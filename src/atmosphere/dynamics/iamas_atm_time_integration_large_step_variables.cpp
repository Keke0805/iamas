/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop large_step_variables subroutine in kokkos

    * Author: Xin Qu
    Last update: Nov 2025

    Contact:
    Chun Zhao, PhD
    Professor
    University of Science and Technology of China
    Hefei, Anhui, 230026
    Email: ChunZhao@ustc.edu.cn

    Please report any bugs or problems to Chun Zhao

    Additional information:
    aemol.ustc.edu.cn: resources
**/
#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>
#include <flcl-cxx.hpp>
#include <cstring>
#include <iostream>
typedef double real;
using real1d = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real2d = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real3d = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real4d = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real5d = Kokkos::Experimental::OffsetView<real*****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int2d = Kokkos::Experimental::OffsetView<int**, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int5d = Kokkos::Experimental::OffsetView<int*****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int1d = Kokkos::Experimental::OffsetView<int*, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int3d = Kokkos::Experimental::OffsetView<int***, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;

extern "C" void atm_recover_large_step_variables_kokkos(int &nVertLevels,int &num_scalars, int &maxEdges, int &nCells, int &nEdges, int &nCellsSolve, int &nEdgesSolve, real &dt, int &ns, int &rk_step, 
    flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* rw_save_p, flcl_ndarray_t* w_p, flcl_ndarray_t* rw_f, flcl_ndarray_t* rw_p_p, flcl_ndarray_t* rtheta_p_p, flcl_ndarray_t* rtheta_pp_p, flcl_ndarray_t* rtheta_p_save_p, flcl_ndarray_t* rt_diabatic_tend_p, flcl_ndarray_t* rho_p_p, 
    flcl_ndarray_t* rho_p_save_p, flcl_ndarray_t* rho_pp_p, flcl_ndarray_t* rho_zz_p, flcl_ndarray_t* rho_base_p, flcl_ndarray_t* ruAvg_p, flcl_ndarray_t* ru_save_p, flcl_ndarray_t* ru_p_p, flcl_ndarray_t* u_p, flcl_ndarray_t* ru_f, flcl_ndarray_t* exner_p, flcl_ndarray_t* exner_base_p, 
    flcl_ndarray_t* rtheta_base_p, flcl_ndarray_t* pressure_p_p, flcl_ndarray_t* zz_p, flcl_ndarray_t* theta_m_p, flcl_ndarray_t* pressure_b_p, flcl_ndarray_t* scalars_p, flcl_ndarray_t* fzm_p, flcl_ndarray_t* fzp_p, 
    flcl_ndarray_t* phi_p_p, flcl_ndarray_t* phi_base_p, 
    flcl_ndarray_t* zb_p, flcl_ndarray_t* zb3_p, flcl_ndarray_t* zb_cell_p, flcl_ndarray_t* zb3_cell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* nEdgesOnCell_p, 
    real &cf1, real &cf2, real &cf3, 
    flcl_ndarray_t* bdyMaskCell_p,   
    int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd, 
    int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd){
    auto wwAvg_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p); auto wwAvg = real2d(wwAvg_raw,{1,1});
    auto rw_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_save_p); auto rw_save = real2d(rw_save_raw,{1,1});
    auto w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*w_p); auto w = real2d(w_raw,{1,1});
    auto rw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_f); auto rw = real2d(rw_raw,{1,1});
    auto rw_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_p_p); auto rw_p = real2d(rw_p_raw,{1,1});
    auto rtheta_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_p_p); auto rtheta_p = real2d(rtheta_p_raw,{1,1});
    auto rtheta_pp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_pp_p); auto rtheta_pp = real2d(rtheta_pp_raw,{1,1});
    auto rtheta_p_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_p_save_p); auto rtheta_p_save = real2d(rtheta_p_save_raw,{1,1});
    auto rt_diabatic_tend_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rt_diabatic_tend_p); auto rt_diabatic_tend = real2d(rt_diabatic_tend_raw,{1,1});
    auto rho_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_p_p); auto rho_p = real2d(rho_p_raw,{1,1});
    auto rho_p_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_p_save_p); auto rho_p_save = real2d(rho_p_save_raw,{1,1});
    auto rho_pp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_pp_p); auto rho_pp = real2d(rho_pp_raw,{1,1});
    auto rho_zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_p); auto rho_zz = real2d(rho_zz_raw,{1,1});
    auto rho_base_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_base_p); auto rho_base = real2d(rho_base_raw,{1,1});
    auto ruAvg_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ruAvg_p); auto ruAvg = real2d(ruAvg_raw,{1,1});
    auto ru_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_save_p); auto ru_save = real2d(ru_save_raw,{1,1});
    auto ru_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_p_p); auto ru_p = real2d(ru_p_raw,{1,1});
    auto u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*u_p); auto u = real2d(u_raw,{1,1});
    auto ru_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_f); auto ru = real2d(ru_raw,{1,1});
    auto exner_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*exner_p); auto exner = real2d(exner_raw,{1,1});
    auto exner_base_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*exner_base_p); auto exner_base = real2d(exner_base_raw,{1,1});
    auto rtheta_base_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_base_p); auto rtheta_base = real2d(rtheta_base_raw,{1,1});
    auto pressure_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pressure_p_p); auto pressure_p = real2d(pressure_p_raw,{1,1});
    auto zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zz_p); auto zz = real2d(zz_raw,{1,1});
    auto theta_m_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*theta_m_p); auto theta_m = real2d(theta_m_raw,{1,1});
    auto pressure_b_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pressure_b_p); auto pressure_b = real2d(pressure_b_raw,{1,1});
    auto scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_p); auto scalars = real3d(scalars_raw,{1,1,1});
    auto fzm_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzm_p); auto fzm = real1d(fzm_raw,{1});
    auto fzp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzp_p); auto fzp = real1d(fzp_raw,{1});
    auto phi_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*phi_p_p); auto phi_p = real2d(phi_p_raw,{1,1});
    auto phi_base_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*phi_base_p); auto phi_base = real2d(phi_base_raw,{1,1});
    auto zb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*zb_p); auto zb = real3d(zb_raw,{1,1,1});
    auto zb3_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*zb3_p); auto zb3 = real3d(zb3_raw,{1,1,1});
    auto zb_cell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*zb_cell_p); auto zb_cell = real3d(zb_cell_raw,{1,1,1});
    auto zb3_cell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*zb3_cell_p); auto zb3_cell = real3d(zb3_cell_raw,{1,1,1});
    auto edgesOnCell_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p); auto edgesOnCell_sign = real2d(edgesOnCell_sign_raw,{1,1});
    auto cellsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p); auto cellsOnEdge = int2d(cellsOnEdge_raw,{1,1});
    auto edgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p); auto edgesOnCell = int2d(edgesOnCell_raw,{1,1});
    auto nEdgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p); auto nEdgesOnCell = int1d(nEdgesOnCell_raw,{1});
    auto bdyMaskCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskCell_p); auto bdyMaskCell = int1d(bdyMaskCell_raw,{1});
    const int nRelaxZone = 5;
    real rgas = 287.0;
    real cp = 7. * rgas / 2.;
    real cv = cp - rgas;
    real rcv = rgas / (cp - rgas);
    real p0 = 1.0e+05; // this should come from somewhere else...

    // Avoid FP errors caused by a potential division by zero below by 
    // initializing the "garbage cell" of rho_zz to a non-zero value
    for(int k = 1; k <= nVertLevels; k++){
        rho_zz(k, nCells + 1) = 1.0;
    }

    // compute new density everywhere so we can compute u from ru.
    // we will also need it to compute theta_m below

    real invNs = 1.0 / static_cast<real>(ns);

    Kokkos::parallel_for(
        "Cell_kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
        for(int k = 1; k <= nVertLevels; k++){
            rho_p(k, iCell) = rho_p_save(k, iCell) + rho_pp(k, iCell);
            rho_zz(k, iCell) = rho_p(k, iCell) + rho_base(k, iCell);
        }
        w(1, iCell) = 0.0;
        for(int k = 2; k <= nVertLevels; k++){
            wwAvg(k,iCell) = rw_save(k,iCell) + (wwAvg(k,iCell) * invNs);
            rw(k,iCell) = rw_save(k,iCell) + rw_p(k,iCell);

            // pick up part of diagnosed w from omega - divide by density later
            w(k,iCell) = rw(k,iCell)/(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell));
        }
        w(nVertLevels+1,iCell) = 0.0;
        if (rk_step == 3){
            for(int k = 1; k <= nVertLevels; k++){
                rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell) 
                                 - dt * rho_zz(k,iCell) * rt_diabatic_tend(k,iCell);
                theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell);
                exner(k,iCell) = Kokkos::pow((zz(k,iCell)*(rgas/p0)*(rtheta_p(k,iCell)+rtheta_base(k,iCell))),rcv);
                // phi=ln(p/p0)
                phi_p(k,iCell) = cp/cv * log( rgas/p0 * zz(k,iCell) * (rtheta_p(k,iCell)+rtheta_base(k,iCell)) ) - phi_base(k,iCell);
                // pressure_p is perturbation pressure
                // pressure_p(k,iCell) = p0 * exp(phi_p(k,iCell) + phi_base(k,iCell)) - pressure_b(k,iCell)
                pressure_p(k,iCell) = zz(k,iCell) * rgas * (exner(k,iCell)*rtheta_p(k,iCell)+rtheta_base(k,iCell) 
                                                          * (exner(k,iCell)-exner_base(k,iCell)));
            }
        }else{
            for(int k = 1; k <= nVertLevels; k++){
                rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell);
                theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell);
            }
        }
    // }
    });
    Kokkos::fence();

    /** 
        recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).  
        we solved for these in the acoustic-step loop.  
        we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve
    */
   Kokkos::parallel_for(
    "Edges_kokkos", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
    // for(int iEdge = edgeStart; iEdge <= edgeEnd; iEdge++){
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);

        for(int k = 1; k <= nVertLevels; k++){
            ruAvg(k,iEdge) = ru_save(k,iEdge) + (ruAvg(k,iEdge) * invNs);
            ru(k,iEdge) = ru_save(k,iEdge) + ru_p(k,iEdge);
            u(k,iEdge) = 2.*ru(k,iEdge)/(rho_zz(k,cell1)+rho_zz(k,cell2));
        }
    // }
    });
    Kokkos::fence();

    Kokkos::parallel_for(
        "Cell_kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
        /** 
          finish recovering w from (rho*omega)_p.  as when we formed (rho*omega)_p from u and w, we need
          to use the same flux-divergence operator as is used for the horizontal theta transport
          (See Klemp et al 2003).
        */
       if (bdyMaskCell(iCell) <= nRelaxZone){
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iEdge=edgesOnCell(i,iCell);

                real flux = (cf1*ru(1,iEdge) + cf2*ru(2,iEdge) + cf3*ru(3,iEdge));
                w(1,iCell) = w(1,iCell) + edgesOnCell_sign(i,iCell) * 
                                    (zb_cell(1,i,iCell) + Kokkos::copysign(1.0,flux)*zb3_cell(1,i,iCell))*flux;
                for(int k = 2; k <= nVertLevels; k++){
                    real flux = (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge));
                    w(k,iCell) = w(k,iCell) + edgesOnCell_sign(i,iCell) * 
                                        (zb_cell(k,i,iCell)+Kokkos::copysign(1.0,flux)*zb3_cell(k,i,iCell))*flux;
                }
            }
            w(1,iCell) = w(1,iCell)/(cf1*rho_zz(1,iCell)+cf2*rho_zz(2,iCell)+cf3*rho_zz(3,iCell));
            for(int k = 2; k <= nVertLevels; k++){
                w(k,iCell) = w(k,iCell)/(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell));
            }
        }
    // }
    });
    Kokkos::fence();
}
