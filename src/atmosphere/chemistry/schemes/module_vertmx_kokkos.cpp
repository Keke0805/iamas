// Kokkos translation of module_vertmx_wrf plus its drydep gather/writeback.
// Different (column,species) systems are independent. A system's vertical
// Thomas forward/backward sweeps must remain ordered.
#include <Kokkos_Core.hpp>
#include "flcl-cxx.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>

extern "C" void iamas_vertmx_stage(int stage, int begin);
namespace {
using Exec = Kokkos::DefaultHostExecutionSpace;
using Work2 = Kokkos::View<double**, Kokkos::LayoutLeft, Exec>;
using Work3 = Kokkos::View<double***, Kokkos::LayoutLeft, Exec>;
constexpr int task_chunk = 64;
}

extern "C" void iamas_vertmx_kokkos_batch(
    double dt, double eps, int gas_only,
    flcl_ndarray_t *chem_p, flcl_ndarray_t *alt_p,
    flcl_ndarray_t *zgrid_p, flcl_ndarray_t *kzh_p,
    flcl_ndarray_t *ddvel_p, flcl_ndarray_t *aerosol_p) {
    if (!Kokkos::is_initialized()) Kokkos::abort("vertmx requires initialized Kokkos");
    auto chem = flcl::view_from_ndarray<double***>(*chem_p);
    auto alt = flcl::view_from_ndarray<double**>(*alt_p);
    auto zgrid = flcl::view_from_ndarray<double**>(*zgrid_p);
    auto kzh = flcl::view_from_ndarray<double**>(*kzh_p);
    auto ddvel = flcl::view_from_ndarray<double**>(*ddvel_p);
    auto aerosol = flcl::view_from_ndarray<int*>(*aerosol_p);
    if (chem.extent(0) > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        chem.extent(1) > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        chem.extent(2) > static_cast<size_t>(std::numeric_limits<int>::max()))
        Kokkos::abort("vertmx dimensions exceed supported integer range");
    const int ns = static_cast<int>(chem.extent(0));
    const int nk = static_cast<int>(chem.extent(1));
    const int nc = static_cast<int>(chem.extent(2));
    if (nk < 2) Kokkos::abort("vertmx requires at least two levels");
    if (ns == 0 || nc == 0) return;
    const std::int64_t tasks = static_cast<std::int64_t>(ns) * nc;
    const int slots = static_cast<int>(std::min<std::int64_t>(task_chunk, tasks));
    Exec exec;

    iamas_vertmx_stage(0, 1);
    // Local ownership avoids shared mutable caches when called by host threads.
    // LayoutLeft gives each system a contiguous vertical workspace. The total
    // per-call solve workspace is bounded by 64 systems, irrespective of nc.
    Work2 rho("vertmx_rho", nk, nc), a("vertmx_a", nk, nc);
    Work2 b("vertmx_b", nk, nc), turb("vertmx_turb", nk, nc);
    Work3 work("vertmx_system_workspace", nk, slots, 6);
    exec.fence();
    iamas_vertmx_stage(0, 0);

    iamas_vertmx_stage(1, 1);
    Kokkos::parallel_for("vertmx_prepare_columns",
        Kokkos::RangePolicy<Exec>(exec, 0, nc), KOKKOS_LAMBDA(int c) {
        const double base = zgrid(0,c);
        for (int k = 0; k < nk; ++k) {
            rho(k,c) = 1.0 / alt(k,c);
            turb(k,c) = k == 0 ? 0.0 : (kzh(k,c) < 1.e-6 ? 1.e-6 : kzh(k,c));
            // Preserve the original relative-height subtraction order.
            const double z0 = zgrid(k,c) - base;
            const double z1 = zgrid(k+1,c) - base;
            b(k,c) = 1.0 / (rho(k,c) * (z1-z0));
        }
        a(0,c) = 0.0; // Unused lower coefficient, initialized for diagnostics.
        for (int k = 1; k < nk; ++k) {
            const double h0 = (zgrid(k,c)+zgrid(k-1,c))/2.0 - base;
            const double h1 = (zgrid(k+1,c)+zgrid(k,c))/2.0 - base;
            const double rho_at_w = 0.5 * (rho(k,c)+rho(k-1,c));
            a(k,c) = rho_at_w / (h1-h0);
        }
    });
    exec.fence();
    iamas_vertmx_stage(1, 0);

    iamas_vertmx_stage(2, 1);
    for (std::int64_t first = 0; first < tasks; first += task_chunk) {
        const int length = static_cast<int>(std::min<std::int64_t>(task_chunk, tasks-first));
        Kokkos::parallel_for("vertmx_solve_column_species",
            Kokkos::RangePolicy<Exec>(exec, 0, length), KOKKOS_LAMBDA(int slot) {
            const std::int64_t task = first + slot;
            const int s = static_cast<int>(task % ns);
            const int c = static_cast<int>(task / ns);
            // work planes: 0=phi, 1=lhs1, 2=lhs2, 3=lhs3, 4=rhs, 5=q.
            for (int k = 0; k < nk; ++k)
                work(k,slot,0) = chem(s,k,c) < eps ? eps : chem(s,k,c);
            if (!gas_only || aerosol(s) == 0) {
                const double vd = ddvel(s,c);
                const double alfa = 0.25, beta = 0.75;
                const double a2 = a(1,c)*turb(1,c);
                work(0,slot,4) = (1.0/(dt*b(0,c))-alfa*(vd*rho(0,c)+a2))*work(0,slot,0)
                              + alfa*(a2*work(1,slot,0));
                work(0,slot,1) = 0.0;
                work(0,slot,2) = 1.0/(dt*b(0,c))+beta*(vd*rho(0,c)+a2);
                work(0,slot,3) = -beta*a2;
                for (int k = 1; k < nk-1; ++k) {
                    const double ak1 = a(k,c)*turb(k,c);
                    const double ak2 = a(k+1,c)*turb(k+1,c);
                    work(k,slot,4) = (1.0/(dt*b(k,c))-alfa*(ak1+ak2))*work(k,slot,0)
                        + alfa*(ak1*work(k-1,slot,0)+ak2*work(k+1,slot,0));
                    work(k,slot,1) = -beta*ak1;
                    work(k,slot,2) = 1.0/(dt*b(k,c))+beta*(ak1+ak2);
                    work(k,slot,3) = -beta*ak2;
                }
                const int top = nk-1;
                const double ak1 = a(top,c)*turb(top,c);
                work(top,slot,4) = (1.0/(dt*b(top,c))-alfa*ak1)*work(top,slot,0)
                                + alfa*(ak1*work(top-1,slot,0));
                work(top,slot,1) = -beta*ak1;
                work(top,slot,2) = 1.0/(dt*b(top,c))+beta*ak1;
                work(top,slot,3) = 0.0;
                // Thomas solver: preserve the two ordered vertical sweeps.
                work(0,slot,5) = -work(0,slot,3)/work(0,slot,2);
                work(0,slot,4) = work(0,slot,4)/work(0,slot,2);
                for (int k = 1; k < nk; ++k) {
                    const double p = 1.0/(work(k,slot,2)+work(k,slot,1)*work(k-1,slot,5));
                    work(k,slot,5) = -work(k,slot,3)*p;
                    work(k,slot,4) = (work(k,slot,4)-work(k,slot,1)*work(k-1,slot,4))*p;
                }
                for (int k = nk-2; k >= 0; --k)
                    work(k,slot,4) = work(k,slot,4)+work(k,slot,5)*work(k+1,slot,4);
                for (int k = 0; k < nk; ++k) work(k,slot,0) = work(k,slot,4);
            }
            // Match drydep exactly: top concentration remains unchanged,
            // including the gas-only aerosol path and its lower-level floor.
            for (int k = 0; k < nk-1; ++k) {
                const double value = work(k,slot,0);
                chem(s,k,c) = value < eps ? eps : value;
            }
        });
        exec.fence(); // Workspace slots may only be reused after completion.
    }
    iamas_vertmx_stage(2, 0);
}
