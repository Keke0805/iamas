/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop advance_scalars_mono_precomm subroutine in kokkos

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
#include <flcl-cxx.hpp>
#include <iostream>
#include <cstdlib>
#include <cstring>
extern "C" void iamas_precomm_stage(int, int);
struct PrecommStage {
 int id;
 explicit PrecommStage(int i):id(i){Kokkos::fence();iamas_precomm_stage(id,1);}
 ~PrecommStage(){Kokkos::fence();iamas_precomm_stage(id,0);}
};
typedef double real;


KOKKOS_INLINE_FUNCTION real flux4(const real &q_im2, const real &q_im1, const real &q_i, const real &q_ip1, const real &ua){
	return ua * (7. * (q_i + q_im1) - (q_ip1 + q_im2) ) / 12.0;
}

KOKKOS_INLINE_FUNCTION real flux3(const real &q_im2, const real &q_im1, const real &q_i, const real &q_ip1, const real &ua, const real &coef3){
	return flux4(q_im2, q_im1, q_i, q_ip1, ua) + coef3 * Kokkos::abs(ua) * ((q_ip1 - q_im2) - 3. * (q_i - q_im1)) / 12.0;
}

KOKKOS_INLINE_FUNCTION real vert_upwind_flux(const real &dt, const real &w, const real &q_below, const real &q_here){
    return dt * (Kokkos::max(0.0, w) * q_below + Kokkos::min(0.0, w) * q_here);
}

static void precomm_baseline(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
	auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_old_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto cellsOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
	auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
	auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
	auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
	auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
	auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
	auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
	auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
	auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
	auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
	auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
	const int SCALE_IN = 0, SCALE_OUT = 1;
	const double eps = 1.e-20;
	const int nRelaxZone = 5;
	bool local_advance_density;
	local_advance_density = advance_density;
    // Each stage is launched once over all scalars. This keeps the algorithmic
    // ordering between stages while avoiding num_scalars_dummy repeated kernel launches.
    { PrecommStage stage(0);
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });
    }

    { PrecommStage stage(1);
    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });
    }

    { PrecommStage stage(2);
    Kokkos::parallel_for("precomm_vert_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            wdtn_arr_scalars(1, iCell, iScalar) =
                wwAvg(1, iCell) *
                (fnm(1) * scalar_new_arr_scalars(1, iCell, iScalar) +
                 fnp(1) * scalar_new_arr_scalars(0, iCell, iScalar));

            for (int k = 2; k < nVertLevels_dummy - 1; k++) {
                wdtn_arr_scalars(k, iCell, iScalar) = flux3(
                    scalar_new_arr_scalars(k-2, iCell, iScalar),
                    scalar_new_arr_scalars(k-1, iCell, iScalar),
                    scalar_new_arr_scalars(k,   iCell, iScalar),
                    scalar_new_arr_scalars(k+1, iCell, iScalar),
                    wwAvg(k, iCell), coef_3rd_order);
            }

            int k = nVertLevels_dummy - 1;
            wdtn_arr_scalars(k, iCell, iScalar) =
                wwAvg(k, iCell) *
                (fnm(k) * scalar_new_arr_scalars(k, iCell, iScalar) +
                 fnp(k) * scalar_new_arr_scalars(k-1, iCell, iScalar));
    });
    }

    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;

            if (cell1 > nCellsSolve - 1 && cell2 > nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                }
            } else {
                const int nAdv = nAdvCellsForEdge(iEdge);
                if (nAdv == 10) {
                    int ica[10];
                    real coef_pos[10];
                    real coef_neg[10];

                    for (int i = 0; i < 10; i++) {
                        ica[i] = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        coef_pos[i] = c + c3;
                        coef_neg[i] = c - c3;
                    }

                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        const real uh = uhAvg(k, iEdge);
                        const real* coef = (uh > 0.0) ? coef_pos : coef_neg;
                        flux_arr_scalars(k, iEdge, iScalar) = uh * (
                            coef[0]*scalar_new_arr_scalars(k, ica[0], iScalar) +
                            coef[1]*scalar_new_arr_scalars(k, ica[1], iScalar) +
                            coef[2]*scalar_new_arr_scalars(k, ica[2], iScalar) +
                            coef[3]*scalar_new_arr_scalars(k, ica[3], iScalar) +
                            coef[4]*scalar_new_arr_scalars(k, ica[4], iScalar) +
                            coef[5]*scalar_new_arr_scalars(k, ica[5], iScalar) +
                            coef[6]*scalar_new_arr_scalars(k, ica[6], iScalar) +
                            coef[7]*scalar_new_arr_scalars(k, ica[7], iScalar) +
                            coef[8]*scalar_new_arr_scalars(k, ica[8], iScalar) +
                            coef[9]*scalar_new_arr_scalars(k, ica[9], iScalar));
                    }
                } else {
                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                    }
                    for (int i = 0; i < nAdv; i++) {
                        int iCell = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        for (int k = 0; k < nVertLevels_dummy; k++) {
                            real uh = uhAvg(k, iEdge);
                            real coef = (uh > 0.0) ? (c + c3) : (c - c3);
                            flux_arr_scalars(k, iEdge, iScalar) +=
                                uh * coef * scalar_new_arr_scalars(k, iCell, iScalar);
                        }
                    }
                }
            }
    });
    }

    { PrecommStage stage(4);
    Kokkos::parallel_for("precomm_vert_upwind_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            for (int k = 0; k < nVertLevels_dummy; k++) {
                scalar_new_arr_scalars(k, iCell, iScalar) =
                    scalar_old_arr_scalars(k, iCell, iScalar) * rho_zz_old(k, iCell);
            }

            for (int k = 0; k < nVertLevels_dummy; k++) {
                real flux_in = 0.0;
                real flux_out = 0.0;
                if (k > 0) {
                    flux_in = vert_upwind_flux(dt, wwAvg(k, iCell),
                        scalar_old_arr_scalars(k-1, iCell, iScalar),
                        scalar_old_arr_scalars(k, iCell, iScalar));
                    wdtn_arr_scalars(k, iCell, iScalar) =
                        dt * wdtn_arr_scalars(k, iCell, iScalar) - flux_in;
                }
                if (k < nVertLevels_dummy - 1) {
                    flux_out = vert_upwind_flux(dt, wwAvg(k+1, iCell),
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }
                scalar_new_arr_scalars(k, iCell, iScalar) +=
                    (flux_in - flux_out) * rdnw(k);
            }

            for (int k = 0; k < nVertLevels_dummy; k++) {
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::min(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::max(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::max(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::min(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
            }
    });
    }

    { PrecommStage stage(5);
    Kokkos::parallel_for("precomm_horiz_lbc_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    flux_arr_scalars(k, iEdge, iScalar) = flux_upwind;
                }
            }
    });
    }

    { PrecommStage stage(6);
    Kokkos::parallel_for("precomm_horiz_div_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            for (int i = 0; i < nEdgesOnCell(iCell); i++) {
                int iEdge = edgesOnCell(i, iCell) - 1;
                real sign = edgesOnCell_sign(i, iCell);
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                const bool lbc_edge =
                    (config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                    (bdyMaskEdge(iEdge) == nRelaxZone - 1);

                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    real flux_tmp_local =
                        lbc_edge ? 0.0 :
                        dt * flux_arr_scalars(k, iEdge, iScalar) - flux_upwind;

                    scalar_new_arr_scalars(k, iCell, iScalar) -=
                        sign * flux_upwind * invAreaCell(iCell);
                    scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) -=
                        Kokkos::max(0.0, sign * flux_tmp_local) * invAreaCell(iCell);
                    scale_arr_scalars(k, SCALE_IN, iScalar, iCell) -=
                        Kokkos::min(0.0, sign * flux_tmp_local) * invAreaCell(iCell);
                }
            }
    });
    }

    { PrecommStage stage(7);
    Kokkos::parallel_for("precomm_limiter_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            for (int k = 0; k < nVertLevels_dummy; k++) {
                real local_smax;
                real local_smin;

                if (k == 0) {
                    local_smax = Kokkos::max(
                        scalar_old_arr_scalars(0, iCell, iScalar),
                        scalar_old_arr_scalars(1, iCell, iScalar));
                    local_smin = Kokkos::min(
                        scalar_old_arr_scalars(0, iCell, iScalar),
                        scalar_old_arr_scalars(1, iCell, iScalar));
                } else if (k == nVertLevels_dummy - 1) {
                    local_smax = Kokkos::max(
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k-1, iCell, iScalar));
                    local_smin = Kokkos::min(
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k-1, iCell, iScalar));
                } else {
                    local_smax = Kokkos::max(
                        Kokkos::max(
                            scalar_old_arr_scalars(k-1, iCell, iScalar),
                            scalar_old_arr_scalars(k, iCell, iScalar)),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                    local_smin = Kokkos::min(
                        Kokkos::min(
                            scalar_old_arr_scalars(k-1, iCell, iScalar),
                            scalar_old_arr_scalars(k, iCell, iScalar)),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }

                for (int i = 0; i < nEdgesOnCell(iCell); i++) {
                    int nbr = cellsOnCell(i, iCell) - 1;
                    real qnbr = scalar_old_arr_scalars(k, nbr, iScalar);
                    local_smax = Kokkos::max(local_smax, qnbr);
                    local_smin = Kokkos::min(local_smin, qnbr);
                }

                real rho = local_advance_density ? rho_zz_int(k, iCell) : rho_zz_new(k, iCell);
                real sf_in =
                    (local_smax * rho - scalar_new_arr_scalars(k, iCell, iScalar)) /
                    (scale_arr_scalars(k, SCALE_IN, iScalar, iCell) + eps);
                real sf_out =
                    (local_smin * rho - scalar_new_arr_scalars(k, iCell, iScalar)) /
                    (scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) - eps);
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    Kokkos::min(1.0, Kokkos::max(0.0, sf_in));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    Kokkos::min(1.0, Kokkos::max(0.0, sf_out));
            }
    });
    }
    Kokkos::fence();
}

static void precomm_fused(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
	auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_old_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto cellsOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
	auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
	auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
	auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
	auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
	auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
	auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
	auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
	auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
	auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
	auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
	const int SCALE_IN = 0, SCALE_OUT = 1;
	const double eps = 1.e-20;
	const int nRelaxZone = 5;
	bool local_advance_density;
	local_advance_density = advance_density;
    // Each stage is launched once over all scalars. This keeps the algorithmic
    // ordering between stages while avoiding num_scalars_dummy repeated kernel launches.
    { PrecommStage stage(0);
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });
    }

    { PrecommStage stage(1);
    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });
    }

    { PrecommStage stage(2);
    Kokkos::parallel_for("precomm_vert_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            wdtn_arr_scalars(1, iCell, iScalar) =
                wwAvg(1, iCell) *
                (fnm(1) * scalar_new_arr_scalars(1, iCell, iScalar) +
                 fnp(1) * scalar_new_arr_scalars(0, iCell, iScalar));

            for (int k = 2; k < nVertLevels_dummy - 1; k++) {
                wdtn_arr_scalars(k, iCell, iScalar) = flux3(
                    scalar_new_arr_scalars(k-2, iCell, iScalar),
                    scalar_new_arr_scalars(k-1, iCell, iScalar),
                    scalar_new_arr_scalars(k,   iCell, iScalar),
                    scalar_new_arr_scalars(k+1, iCell, iScalar),
                    wwAvg(k, iCell), coef_3rd_order);
            }

            int k = nVertLevels_dummy - 1;
            wdtn_arr_scalars(k, iCell, iScalar) =
                wwAvg(k, iCell) *
                (fnm(k) * scalar_new_arr_scalars(k, iCell, iScalar) +
                 fnp(k) * scalar_new_arr_scalars(k-1, iCell, iScalar));
    });
    }

    if (false) {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_species4",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,edgeStart-1}, {(num_scalars_dummy+3)/4,edgeEnd}),
        KOKKOS_LAMBDA(int species_block, int iEdge) {
            const int c1 = cellsOnEdge(0,iEdge)-1;
            const int c2 = cellsOnEdge(1,iEdge)-1;
            const bool active = c1 < nCellsSolve || c2 < nCellsSolve;
            const int nAdv = active ? nAdvCellsForEdge(iEdge) : 0;
            const int sb = species_block*4;
            const int se = Kokkos::min(sb+4,num_scalars_dummy);
            // The validated mesh has a maximum of 10 stencil cells. Retain a
            // generic fallback for larger stencils rather than truncating them.
            if (nAdv <= 10) {
                int cell[10]; real positive[10],negative[10];
                for(int i=0;i<nAdv;++i){
                    cell[i]=advCellsForEdge(i,iEdge)-1;
                    const real c=adv_coefs(i,iEdge), c3=adv_coefs_3rd(i,iEdge);
                    positive[i]=c+c3; negative[i]=c-c3;
                }
                for(int sp=sb;sp<se;++sp){
                    
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);
                        const real* cf=uh>0 ? positive:negative;
                        real value=0.0;
                        if(nAdv==10){
                            // Preserve the original specialized summation order.
                            value=uh*(cf[0]*scalar_new_arr_scalars(k,cell[0],sp)+
                                      cf[1]*scalar_new_arr_scalars(k,cell[1],sp)+
                                      cf[2]*scalar_new_arr_scalars(k,cell[2],sp)+
                                      cf[3]*scalar_new_arr_scalars(k,cell[3],sp)+
                                      cf[4]*scalar_new_arr_scalars(k,cell[4],sp)+
                                      cf[5]*scalar_new_arr_scalars(k,cell[5],sp)+
                                      cf[6]*scalar_new_arr_scalars(k,cell[6],sp)+
                                      cf[7]*scalar_new_arr_scalars(k,cell[7],sp)+
                                      cf[8]*scalar_new_arr_scalars(k,cell[8],sp)+
                                      cf[9]*scalar_new_arr_scalars(k,cell[9],sp));
                        }else{
                            for(int i=0;i<nAdv;++i)value+=uh*cf[i]*scalar_new_arr_scalars(k,cell[i],sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }else{
                for(int sp=sb;sp<se;++sp){
                    
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);real value=0.0;
                        for(int i=0;i<nAdv;++i){
                            const real c=adv_coefs(i,iEdge),c3=adv_coefs_3rd(i,iEdge);
                            const real cf=uh>0 ? c+c3:c-c3;
                            value+=uh*cf*scalar_new_arr_scalars(k,advCellsForEdge(i,iEdge)-1,sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }
    });
    }
    } else {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;

            if (cell1 > nCellsSolve - 1 && cell2 > nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                }
            } else {
                const int nAdv = nAdvCellsForEdge(iEdge);
                if (nAdv == 10) {
                    int ica[10];
                    real coef_pos[10];
                    real coef_neg[10];

                    for (int i = 0; i < 10; i++) {
                        ica[i] = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        coef_pos[i] = c + c3;
                        coef_neg[i] = c - c3;
                    }

                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        const real uh = uhAvg(k, iEdge);
                        const real* coef = (uh > 0.0) ? coef_pos : coef_neg;
                        flux_arr_scalars(k, iEdge, iScalar) = uh * (
                            coef[0]*scalar_new_arr_scalars(k, ica[0], iScalar) +
                            coef[1]*scalar_new_arr_scalars(k, ica[1], iScalar) +
                            coef[2]*scalar_new_arr_scalars(k, ica[2], iScalar) +
                            coef[3]*scalar_new_arr_scalars(k, ica[3], iScalar) +
                            coef[4]*scalar_new_arr_scalars(k, ica[4], iScalar) +
                            coef[5]*scalar_new_arr_scalars(k, ica[5], iScalar) +
                            coef[6]*scalar_new_arr_scalars(k, ica[6], iScalar) +
                            coef[7]*scalar_new_arr_scalars(k, ica[7], iScalar) +
                            coef[8]*scalar_new_arr_scalars(k, ica[8], iScalar) +
                            coef[9]*scalar_new_arr_scalars(k, ica[9], iScalar));
                    }
                } else {
                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                    }
                    for (int i = 0; i < nAdv; i++) {
                        int iCell = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        for (int k = 0; k < nVertLevels_dummy; k++) {
                            real uh = uhAvg(k, iEdge);
                            real coef = (uh > 0.0) ? (c + c3) : (c - c3);
                            flux_arr_scalars(k, iEdge, iScalar) +=
                                uh * coef * scalar_new_arr_scalars(k, iCell, iScalar);
                        }
                    }
                }
            }
    });
    }

    }
    { PrecommStage stage(4);
    Kokkos::parallel_for("precomm_vert_upwind_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            for (int k = 0; k < nVertLevels_dummy; k++) {
                real q = scalar_old_arr_scalars(k, iCell, iScalar) * rho_zz_old(k, iCell);
                real flux_in = 0.0;
                real flux_out = 0.0;
                if (k > 0) {
                    flux_in = vert_upwind_flux(dt, wwAvg(k, iCell),
                        scalar_old_arr_scalars(k-1, iCell, iScalar),
                        scalar_old_arr_scalars(k, iCell, iScalar));
                    wdtn_arr_scalars(k, iCell, iScalar) =
                        dt * wdtn_arr_scalars(k, iCell, iScalar) - flux_in;
                }
                if (k < nVertLevels_dummy - 1) {
                    flux_out = vert_upwind_flux(dt, wwAvg(k+1, iCell),
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }
                q += (flux_in - flux_out) * rdnw(k);
                scalar_new_arr_scalars(k,iCell,iScalar) = q;
            }

            for (int k = 0; k < nVertLevels_dummy; k++) {
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::min(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::max(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::max(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::min(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
            }
    });
    }

    { PrecommStage stage(5);
    Kokkos::parallel_for("precomm_horiz_lbc_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    flux_arr_scalars(k, iEdge, iScalar) = flux_upwind;
                }
            }
    });
    }

    { PrecommStage stage(8);
    Kokkos::parallel_for("precomm_div_limiter_tiled",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,cellSolveStart-1}, {num_scalars_dummy,cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar,int iCell){
            const real inv_area=invAreaCell(iCell);
            const int neighbors=nEdgesOnCell(iCell);
            // A bounded, column-private tile retains k as the unit-stride loop.
            // No shared scratch, allocation, global padding or neighbor reorder.
            for(int kb=0;kb<nVertLevels_dummy;kb+=64){
                const int count=Kokkos::min(64,nVertLevels_dummy-kb);
                real q[64],incoming[64],outgoing[64],smax[64],smin[64];
                
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    q[t]=scalar_new_arr_scalars(k,iCell,iScalar);
                    incoming[t]=scale_arr_scalars(k,SCALE_IN,iScalar,iCell);
                    outgoing[t]=scale_arr_scalars(k,SCALE_OUT,iScalar,iCell);
                    const int below=k>0 ? k-1:k;
                    const int above=k<nVertLevels_dummy-1 ? k+1:k;
                    // Match endpoint and interior max/min operation order.
                    if(k==0){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                    }else if(k==nVertLevels_dummy-1){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                    }else{
                        smax[t]=Kokkos::max(Kokkos::max(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                        smin[t]=Kokkos::min(Kokkos::min(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                    }
                }
                for(int i=0;i<neighbors;++i){
                    const int edge=edgesOnCell(i,iCell)-1;
                    const int nbr=cellsOnCell(i,iCell)-1;
                    const real sign=edgesOnCell_sign(i,iCell);
                    const int c1=cellsOnEdge(0,edge)-1;
                    const int c2=cellsOnEdge(1,edge)-1;
                    const bool lbc=(config_apply_lbcs && bdyMaskEdge(edge)==nRelaxZone) || bdyMaskEdge(edge)==nRelaxZone-1;
                    const real edge_dt=dvEdge(edge)*dt;
                    
                    for(int t=0;t<count;++t){
                        const int k=kb+t;
                        const real uh=uhAvg(k,edge);
                        const real upwind=edge_dt*(Kokkos::max(0.0,uh)*scalar_old_arr_scalars(k,c1,iScalar)+Kokkos::min(0.0,uh)*scalar_old_arr_scalars(k,c2,iScalar));
                        const real correction=lbc ? 0.0:dt*flux_arr_scalars(k,edge,iScalar)-upwind;
                        q[t]-=sign*upwind*inv_area;
                        outgoing[t]-=Kokkos::max(0.0,sign*correction)*inv_area;
                        incoming[t]-=Kokkos::min(0.0,sign*correction)*inv_area;
                        const real neighbor=scalar_old_arr_scalars(k,nbr,iScalar);
                        smax[t]=Kokkos::max(smax[t],neighbor);
                        smin[t]=Kokkos::min(smin[t],neighbor);
                    }
                }
                
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    const real rho=local_advance_density ? rho_zz_int(k,iCell):rho_zz_new(k,iCell);
                    const real si=(smax[t]*rho-q[t])/(incoming[t]+eps);
                    const real so=(smin[t]*rho-q[t])/(outgoing[t]-eps);
                    scalar_new_arr_scalars(k,iCell,iScalar)=q[t];
                    scale_arr_scalars(k,SCALE_IN,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,si));
                    scale_arr_scalars(k,SCALE_OUT,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,so));
                }
            }
    });
    }
    Kokkos::fence();
}

static void precomm_reuse(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
	auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_old_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto cellsOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
	auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
	auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
	auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
	auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
	auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
	auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
	auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
	auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
	auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
	auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
	const int SCALE_IN = 0, SCALE_OUT = 1;
	const double eps = 1.e-20;
	const int nRelaxZone = 5;
	bool local_advance_density;
	local_advance_density = advance_density;
    // Each stage is launched once over all scalars. This keeps the algorithmic
    // ordering between stages while avoiding num_scalars_dummy repeated kernel launches.
    { PrecommStage stage(0);
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });
    }

    { PrecommStage stage(1);
    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });
    }

    { PrecommStage stage(2);
    Kokkos::parallel_for("precomm_vert_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            wdtn_arr_scalars(1, iCell, iScalar) =
                wwAvg(1, iCell) *
                (fnm(1) * scalar_new_arr_scalars(1, iCell, iScalar) +
                 fnp(1) * scalar_new_arr_scalars(0, iCell, iScalar));

            for (int k = 2; k < nVertLevels_dummy - 1; k++) {
                wdtn_arr_scalars(k, iCell, iScalar) = flux3(
                    scalar_new_arr_scalars(k-2, iCell, iScalar),
                    scalar_new_arr_scalars(k-1, iCell, iScalar),
                    scalar_new_arr_scalars(k,   iCell, iScalar),
                    scalar_new_arr_scalars(k+1, iCell, iScalar),
                    wwAvg(k, iCell), coef_3rd_order);
            }

            int k = nVertLevels_dummy - 1;
            wdtn_arr_scalars(k, iCell, iScalar) =
                wwAvg(k, iCell) *
                (fnm(k) * scalar_new_arr_scalars(k, iCell, iScalar) +
                 fnp(k) * scalar_new_arr_scalars(k-1, iCell, iScalar));
    });
    }

    if (true) {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_species4",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,edgeStart-1}, {(num_scalars_dummy+3)/4,edgeEnd}),
        KOKKOS_LAMBDA(int species_block, int iEdge) {
            const int c1 = cellsOnEdge(0,iEdge)-1;
            const int c2 = cellsOnEdge(1,iEdge)-1;
            const bool active = c1 < nCellsSolve || c2 < nCellsSolve;
            const int nAdv = active ? nAdvCellsForEdge(iEdge) : 0;
            const int sb = species_block*4;
            const int se = Kokkos::min(sb+4,num_scalars_dummy);
            // The validated mesh has a maximum of 10 stencil cells. Retain a
            // generic fallback for larger stencils rather than truncating them.
            if (nAdv <= 10) {
                int cell[10]; real positive[10],negative[10];
                for(int i=0;i<nAdv;++i){
                    cell[i]=advCellsForEdge(i,iEdge)-1;
                    const real c=adv_coefs(i,iEdge), c3=adv_coefs_3rd(i,iEdge);
                    positive[i]=c+c3; negative[i]=c-c3;
                }
                for(int sp=sb;sp<se;++sp){
                    
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);
                        const real* cf=uh>0 ? positive:negative;
                        real value=0.0;
                        if(nAdv==10){
                            // Preserve the original specialized summation order.
                            value=uh*(cf[0]*scalar_new_arr_scalars(k,cell[0],sp)+
                                      cf[1]*scalar_new_arr_scalars(k,cell[1],sp)+
                                      cf[2]*scalar_new_arr_scalars(k,cell[2],sp)+
                                      cf[3]*scalar_new_arr_scalars(k,cell[3],sp)+
                                      cf[4]*scalar_new_arr_scalars(k,cell[4],sp)+
                                      cf[5]*scalar_new_arr_scalars(k,cell[5],sp)+
                                      cf[6]*scalar_new_arr_scalars(k,cell[6],sp)+
                                      cf[7]*scalar_new_arr_scalars(k,cell[7],sp)+
                                      cf[8]*scalar_new_arr_scalars(k,cell[8],sp)+
                                      cf[9]*scalar_new_arr_scalars(k,cell[9],sp));
                        }else{
                            for(int i=0;i<nAdv;++i)value+=uh*cf[i]*scalar_new_arr_scalars(k,cell[i],sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }else{
                for(int sp=sb;sp<se;++sp){
                    
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);real value=0.0;
                        for(int i=0;i<nAdv;++i){
                            const real c=adv_coefs(i,iEdge),c3=adv_coefs_3rd(i,iEdge);
                            const real cf=uh>0 ? c+c3:c-c3;
                            value+=uh*cf*scalar_new_arr_scalars(k,advCellsForEdge(i,iEdge)-1,sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }
    });
    }
    } else {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;

            if (cell1 > nCellsSolve - 1 && cell2 > nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                }
            } else {
                const int nAdv = nAdvCellsForEdge(iEdge);
                if (nAdv == 10) {
                    int ica[10];
                    real coef_pos[10];
                    real coef_neg[10];

                    for (int i = 0; i < 10; i++) {
                        ica[i] = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        coef_pos[i] = c + c3;
                        coef_neg[i] = c - c3;
                    }

                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        const real uh = uhAvg(k, iEdge);
                        const real* coef = (uh > 0.0) ? coef_pos : coef_neg;
                        flux_arr_scalars(k, iEdge, iScalar) = uh * (
                            coef[0]*scalar_new_arr_scalars(k, ica[0], iScalar) +
                            coef[1]*scalar_new_arr_scalars(k, ica[1], iScalar) +
                            coef[2]*scalar_new_arr_scalars(k, ica[2], iScalar) +
                            coef[3]*scalar_new_arr_scalars(k, ica[3], iScalar) +
                            coef[4]*scalar_new_arr_scalars(k, ica[4], iScalar) +
                            coef[5]*scalar_new_arr_scalars(k, ica[5], iScalar) +
                            coef[6]*scalar_new_arr_scalars(k, ica[6], iScalar) +
                            coef[7]*scalar_new_arr_scalars(k, ica[7], iScalar) +
                            coef[8]*scalar_new_arr_scalars(k, ica[8], iScalar) +
                            coef[9]*scalar_new_arr_scalars(k, ica[9], iScalar));
                    }
                } else {
                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                    }
                    for (int i = 0; i < nAdv; i++) {
                        int iCell = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        for (int k = 0; k < nVertLevels_dummy; k++) {
                            real uh = uhAvg(k, iEdge);
                            real coef = (uh > 0.0) ? (c + c3) : (c - c3);
                            flux_arr_scalars(k, iEdge, iScalar) +=
                                uh * coef * scalar_new_arr_scalars(k, iCell, iScalar);
                        }
                    }
                }
            }
    });
    }

    }
    { PrecommStage stage(4);
    Kokkos::parallel_for("precomm_vert_upwind_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            for (int k = 0; k < nVertLevels_dummy; k++) {
                real q = scalar_old_arr_scalars(k, iCell, iScalar) * rho_zz_old(k, iCell);
                real flux_in = 0.0;
                real flux_out = 0.0;
                if (k > 0) {
                    flux_in = vert_upwind_flux(dt, wwAvg(k, iCell),
                        scalar_old_arr_scalars(k-1, iCell, iScalar),
                        scalar_old_arr_scalars(k, iCell, iScalar));
                    wdtn_arr_scalars(k, iCell, iScalar) =
                        dt * wdtn_arr_scalars(k, iCell, iScalar) - flux_in;
                }
                if (k < nVertLevels_dummy - 1) {
                    flux_out = vert_upwind_flux(dt, wwAvg(k+1, iCell),
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }
                q += (flux_in - flux_out) * rdnw(k);
                scalar_new_arr_scalars(k,iCell,iScalar) = q;
            }

            for (int k = 0; k < nVertLevels_dummy; k++) {
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::min(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::max(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::max(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::min(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
            }
    });
    }

    { PrecommStage stage(5);
    Kokkos::parallel_for("precomm_horiz_lbc_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    flux_arr_scalars(k, iEdge, iScalar) = flux_upwind;
                }
            }
    });
    }

    { PrecommStage stage(8);
    Kokkos::parallel_for("precomm_div_limiter_tiled",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,cellSolveStart-1}, {num_scalars_dummy,cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar,int iCell){
            const real inv_area=invAreaCell(iCell);
            const int neighbors=nEdgesOnCell(iCell);
            // A bounded, column-private tile retains k as the unit-stride loop.
            // No shared scratch, allocation, global padding or neighbor reorder.
            for(int kb=0;kb<nVertLevels_dummy;kb+=64){
                const int count=Kokkos::min(64,nVertLevels_dummy-kb);
                real q[64],incoming[64],outgoing[64],smax[64],smin[64];
                
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    q[t]=scalar_new_arr_scalars(k,iCell,iScalar);
                    incoming[t]=scale_arr_scalars(k,SCALE_IN,iScalar,iCell);
                    outgoing[t]=scale_arr_scalars(k,SCALE_OUT,iScalar,iCell);
                    const int below=k>0 ? k-1:k;
                    const int above=k<nVertLevels_dummy-1 ? k+1:k;
                    // Match endpoint and interior max/min operation order.
                    if(k==0){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                    }else if(k==nVertLevels_dummy-1){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                    }else{
                        smax[t]=Kokkos::max(Kokkos::max(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                        smin[t]=Kokkos::min(Kokkos::min(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                    }
                }
                for(int i=0;i<neighbors;++i){
                    const int edge=edgesOnCell(i,iCell)-1;
                    const int nbr=cellsOnCell(i,iCell)-1;
                    const real sign=edgesOnCell_sign(i,iCell);
                    const int c1=cellsOnEdge(0,edge)-1;
                    const int c2=cellsOnEdge(1,edge)-1;
                    const bool lbc=(config_apply_lbcs && bdyMaskEdge(edge)==nRelaxZone) || bdyMaskEdge(edge)==nRelaxZone-1;
                    const real edge_dt=dvEdge(edge)*dt;
                    
                    for(int t=0;t<count;++t){
                        const int k=kb+t;
                        const real uh=uhAvg(k,edge);
                        const real upwind=edge_dt*(Kokkos::max(0.0,uh)*scalar_old_arr_scalars(k,c1,iScalar)+Kokkos::min(0.0,uh)*scalar_old_arr_scalars(k,c2,iScalar));
                        const real correction=lbc ? 0.0:dt*flux_arr_scalars(k,edge,iScalar)-upwind;
                        q[t]-=sign*upwind*inv_area;
                        outgoing[t]-=Kokkos::max(0.0,sign*correction)*inv_area;
                        incoming[t]-=Kokkos::min(0.0,sign*correction)*inv_area;
                        const real neighbor=scalar_old_arr_scalars(k,nbr,iScalar);
                        smax[t]=Kokkos::max(smax[t],neighbor);
                        smin[t]=Kokkos::min(smin[t],neighbor);
                    }
                }
                
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    const real rho=local_advance_density ? rho_zz_int(k,iCell):rho_zz_new(k,iCell);
                    const real si=(smax[t]*rho-q[t])/(incoming[t]+eps);
                    const real so=(smin[t]*rho-q[t])/(outgoing[t]-eps);
                    scalar_new_arr_scalars(k,iCell,iScalar)=q[t];
                    scale_arr_scalars(k,SCALE_IN,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,si));
                    scale_arr_scalars(k,SCALE_OUT,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,so));
                }
            }
    });
    }
    Kokkos::fence();
}

static void precomm_simd(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
	auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_old_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto cellsOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
	auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
	auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
	auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
	auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
	auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
	auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
	auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
	auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
	auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
	auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
	const int SCALE_IN = 0, SCALE_OUT = 1;
	const double eps = 1.e-20;
	const int nRelaxZone = 5;
	bool local_advance_density;
	local_advance_density = advance_density;
    // Each stage is launched once over all scalars. This keeps the algorithmic
    // ordering between stages while avoiding num_scalars_dummy repeated kernel launches.
    { PrecommStage stage(0);
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });
    }

    { PrecommStage stage(1);
    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });
    }

    { PrecommStage stage(2);
    Kokkos::parallel_for("precomm_vert_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            wdtn_arr_scalars(1, iCell, iScalar) =
                wwAvg(1, iCell) *
                (fnm(1) * scalar_new_arr_scalars(1, iCell, iScalar) +
                 fnp(1) * scalar_new_arr_scalars(0, iCell, iScalar));

            #pragma omp simd
            for (int k = 2; k < nVertLevels_dummy - 1; k++) {
                wdtn_arr_scalars(k, iCell, iScalar) = flux3(
                    scalar_new_arr_scalars(k-2, iCell, iScalar),
                    scalar_new_arr_scalars(k-1, iCell, iScalar),
                    scalar_new_arr_scalars(k,   iCell, iScalar),
                    scalar_new_arr_scalars(k+1, iCell, iScalar),
                    wwAvg(k, iCell), coef_3rd_order);
            }

            int k = nVertLevels_dummy - 1;
            wdtn_arr_scalars(k, iCell, iScalar) =
                wwAvg(k, iCell) *
                (fnm(k) * scalar_new_arr_scalars(k, iCell, iScalar) +
                 fnp(k) * scalar_new_arr_scalars(k-1, iCell, iScalar));
    });
    }

    if (true) {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_species4",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,edgeStart-1}, {(num_scalars_dummy+3)/4,edgeEnd}),
        KOKKOS_LAMBDA(int species_block, int iEdge) {
            const int c1 = cellsOnEdge(0,iEdge)-1;
            const int c2 = cellsOnEdge(1,iEdge)-1;
            const bool active = c1 < nCellsSolve || c2 < nCellsSolve;
            const int nAdv = active ? nAdvCellsForEdge(iEdge) : 0;
            const int sb = species_block*4;
            const int se = Kokkos::min(sb+4,num_scalars_dummy);
            // The validated mesh has a maximum of 10 stencil cells. Retain a
            // generic fallback for larger stencils rather than truncating them.
            if (nAdv <= 10) {
                int cell[10]; real positive[10],negative[10];
                for(int i=0;i<nAdv;++i){
                    cell[i]=advCellsForEdge(i,iEdge)-1;
                    const real c=adv_coefs(i,iEdge), c3=adv_coefs_3rd(i,iEdge);
                    positive[i]=c+c3; negative[i]=c-c3;
                }
                for(int sp=sb;sp<se;++sp){
                    #pragma omp simd
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);
                        const real* cf=uh>0 ? positive:negative;
                        real value=0.0;
                        if(nAdv==10){
                            // Preserve the original specialized summation order.
                            value=uh*(cf[0]*scalar_new_arr_scalars(k,cell[0],sp)+
                                      cf[1]*scalar_new_arr_scalars(k,cell[1],sp)+
                                      cf[2]*scalar_new_arr_scalars(k,cell[2],sp)+
                                      cf[3]*scalar_new_arr_scalars(k,cell[3],sp)+
                                      cf[4]*scalar_new_arr_scalars(k,cell[4],sp)+
                                      cf[5]*scalar_new_arr_scalars(k,cell[5],sp)+
                                      cf[6]*scalar_new_arr_scalars(k,cell[6],sp)+
                                      cf[7]*scalar_new_arr_scalars(k,cell[7],sp)+
                                      cf[8]*scalar_new_arr_scalars(k,cell[8],sp)+
                                      cf[9]*scalar_new_arr_scalars(k,cell[9],sp));
                        }else{
                            for(int i=0;i<nAdv;++i)value+=uh*cf[i]*scalar_new_arr_scalars(k,cell[i],sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }else{
                for(int sp=sb;sp<se;++sp){
                    #pragma omp simd
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);real value=0.0;
                        for(int i=0;i<nAdv;++i){
                            const real c=adv_coefs(i,iEdge),c3=adv_coefs_3rd(i,iEdge);
                            const real cf=uh>0 ? c+c3:c-c3;
                            value+=uh*cf*scalar_new_arr_scalars(k,advCellsForEdge(i,iEdge)-1,sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }
    });
    }
    } else {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;

            if (cell1 > nCellsSolve - 1 && cell2 > nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                }
            } else {
                const int nAdv = nAdvCellsForEdge(iEdge);
                if (nAdv == 10) {
                    int ica[10];
                    real coef_pos[10];
                    real coef_neg[10];

                    for (int i = 0; i < 10; i++) {
                        ica[i] = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        coef_pos[i] = c + c3;
                        coef_neg[i] = c - c3;
                    }

                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        const real uh = uhAvg(k, iEdge);
                        const real* coef = (uh > 0.0) ? coef_pos : coef_neg;
                        flux_arr_scalars(k, iEdge, iScalar) = uh * (
                            coef[0]*scalar_new_arr_scalars(k, ica[0], iScalar) +
                            coef[1]*scalar_new_arr_scalars(k, ica[1], iScalar) +
                            coef[2]*scalar_new_arr_scalars(k, ica[2], iScalar) +
                            coef[3]*scalar_new_arr_scalars(k, ica[3], iScalar) +
                            coef[4]*scalar_new_arr_scalars(k, ica[4], iScalar) +
                            coef[5]*scalar_new_arr_scalars(k, ica[5], iScalar) +
                            coef[6]*scalar_new_arr_scalars(k, ica[6], iScalar) +
                            coef[7]*scalar_new_arr_scalars(k, ica[7], iScalar) +
                            coef[8]*scalar_new_arr_scalars(k, ica[8], iScalar) +
                            coef[9]*scalar_new_arr_scalars(k, ica[9], iScalar));
                    }
                } else {
                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                    }
                    for (int i = 0; i < nAdv; i++) {
                        int iCell = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        for (int k = 0; k < nVertLevels_dummy; k++) {
                            real uh = uhAvg(k, iEdge);
                            real coef = (uh > 0.0) ? (c + c3) : (c - c3);
                            flux_arr_scalars(k, iEdge, iScalar) +=
                                uh * coef * scalar_new_arr_scalars(k, iCell, iScalar);
                        }
                    }
                }
            }
    });
    }

    }
    { PrecommStage stage(4);
    Kokkos::parallel_for("precomm_vert_upwind_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            #pragma omp simd
            for (int k = 0; k < nVertLevels_dummy; k++) {
                real q = scalar_old_arr_scalars(k, iCell, iScalar) * rho_zz_old(k, iCell);
                real flux_in = 0.0;
                real flux_out = 0.0;
                if (k > 0) {
                    flux_in = vert_upwind_flux(dt, wwAvg(k, iCell),
                        scalar_old_arr_scalars(k-1, iCell, iScalar),
                        scalar_old_arr_scalars(k, iCell, iScalar));
                    wdtn_arr_scalars(k, iCell, iScalar) =
                        dt * wdtn_arr_scalars(k, iCell, iScalar) - flux_in;
                }
                if (k < nVertLevels_dummy - 1) {
                    flux_out = vert_upwind_flux(dt, wwAvg(k+1, iCell),
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }
                q += (flux_in - flux_out) * rdnw(k);
                scalar_new_arr_scalars(k,iCell,iScalar) = q;
            }

            #pragma omp simd
            for (int k = 0; k < nVertLevels_dummy; k++) {
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::min(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::max(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::max(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::min(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
            }
    });
    }

    { PrecommStage stage(5);
    Kokkos::parallel_for("precomm_horiz_lbc_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    flux_arr_scalars(k, iEdge, iScalar) = flux_upwind;
                }
            }
    });
    }

    { PrecommStage stage(8);
    Kokkos::parallel_for("precomm_div_limiter_tiled",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,cellSolveStart-1}, {num_scalars_dummy,cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar,int iCell){
            const real inv_area=invAreaCell(iCell);
            const int neighbors=nEdgesOnCell(iCell);
            // A bounded, column-private tile retains k as the unit-stride loop.
            // No shared scratch, allocation, global padding or neighbor reorder.
            for(int kb=0;kb<nVertLevels_dummy;kb+=64){
                const int count=Kokkos::min(64,nVertLevels_dummy-kb);
                real q[64],incoming[64],outgoing[64],smax[64],smin[64];
                #pragma omp simd
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    q[t]=scalar_new_arr_scalars(k,iCell,iScalar);
                    incoming[t]=scale_arr_scalars(k,SCALE_IN,iScalar,iCell);
                    outgoing[t]=scale_arr_scalars(k,SCALE_OUT,iScalar,iCell);
                    const int below=k>0 ? k-1:k;
                    const int above=k<nVertLevels_dummy-1 ? k+1:k;
                    // Match endpoint and interior max/min operation order.
                    if(k==0){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                    }else if(k==nVertLevels_dummy-1){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                    }else{
                        smax[t]=Kokkos::max(Kokkos::max(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                        smin[t]=Kokkos::min(Kokkos::min(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                    }
                }
                for(int i=0;i<neighbors;++i){
                    const int edge=edgesOnCell(i,iCell)-1;
                    const int nbr=cellsOnCell(i,iCell)-1;
                    const real sign=edgesOnCell_sign(i,iCell);
                    const int c1=cellsOnEdge(0,edge)-1;
                    const int c2=cellsOnEdge(1,edge)-1;
                    const bool lbc=(config_apply_lbcs && bdyMaskEdge(edge)==nRelaxZone) || bdyMaskEdge(edge)==nRelaxZone-1;
                    const real edge_dt=dvEdge(edge)*dt;
                    #pragma omp simd
                    for(int t=0;t<count;++t){
                        const int k=kb+t;
                        const real uh=uhAvg(k,edge);
                        const real upwind=edge_dt*(Kokkos::max(0.0,uh)*scalar_old_arr_scalars(k,c1,iScalar)+Kokkos::min(0.0,uh)*scalar_old_arr_scalars(k,c2,iScalar));
                        const real correction=lbc ? 0.0:dt*flux_arr_scalars(k,edge,iScalar)-upwind;
                        q[t]-=sign*upwind*inv_area;
                        outgoing[t]-=Kokkos::max(0.0,sign*correction)*inv_area;
                        incoming[t]-=Kokkos::min(0.0,sign*correction)*inv_area;
                        const real neighbor=scalar_old_arr_scalars(k,nbr,iScalar);
                        smax[t]=Kokkos::max(smax[t],neighbor);
                        smin[t]=Kokkos::min(smin[t],neighbor);
                    }
                }
                #pragma omp simd
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    const real rho=local_advance_density ? rho_zz_int(k,iCell):rho_zz_new(k,iCell);
                    const real si=(smax[t]*rho-q[t])/(incoming[t]+eps);
                    const real so=(smin[t]*rho-q[t])/(outgoing[t]-eps);
                    scalar_new_arr_scalars(k,iCell,iScalar)=q[t];
                    scale_arr_scalars(k,SCALE_IN,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,si));
                    scale_arr_scalars(k,SCALE_OUT,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,so));
                }
            }
    });
    }
    Kokkos::fence();
}

static void precomm_fused_simd(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
	auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_old_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto cellsOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
	auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
	auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
	auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
	auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
	auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
	auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
	auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
	auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
	auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
	auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
	const int SCALE_IN = 0, SCALE_OUT = 1;
	const double eps = 1.e-20;
	const int nRelaxZone = 5;
	bool local_advance_density;
	local_advance_density = advance_density;
    // Each stage is launched once over all scalars. This keeps the algorithmic
    // ordering between stages while avoiding num_scalars_dummy repeated kernel launches.
    { PrecommStage stage(0);
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });
    }

    { PrecommStage stage(1);
    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });
    }

    { PrecommStage stage(2);
    Kokkos::parallel_for("precomm_vert_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            wdtn_arr_scalars(1, iCell, iScalar) =
                wwAvg(1, iCell) *
                (fnm(1) * scalar_new_arr_scalars(1, iCell, iScalar) +
                 fnp(1) * scalar_new_arr_scalars(0, iCell, iScalar));

            #pragma omp simd
            for (int k = 2; k < nVertLevels_dummy - 1; k++) {
                wdtn_arr_scalars(k, iCell, iScalar) = flux3(
                    scalar_new_arr_scalars(k-2, iCell, iScalar),
                    scalar_new_arr_scalars(k-1, iCell, iScalar),
                    scalar_new_arr_scalars(k,   iCell, iScalar),
                    scalar_new_arr_scalars(k+1, iCell, iScalar),
                    wwAvg(k, iCell), coef_3rd_order);
            }

            int k = nVertLevels_dummy - 1;
            wdtn_arr_scalars(k, iCell, iScalar) =
                wwAvg(k, iCell) *
                (fnm(k) * scalar_new_arr_scalars(k, iCell, iScalar) +
                 fnp(k) * scalar_new_arr_scalars(k-1, iCell, iScalar));
    });
    }

    if (false) {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_species4",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,edgeStart-1}, {(num_scalars_dummy+3)/4,edgeEnd}),
        KOKKOS_LAMBDA(int species_block, int iEdge) {
            const int c1 = cellsOnEdge(0,iEdge)-1;
            const int c2 = cellsOnEdge(1,iEdge)-1;
            const bool active = c1 < nCellsSolve || c2 < nCellsSolve;
            const int nAdv = active ? nAdvCellsForEdge(iEdge) : 0;
            const int sb = species_block*4;
            const int se = Kokkos::min(sb+4,num_scalars_dummy);
            // The validated mesh has a maximum of 10 stencil cells. Retain a
            // generic fallback for larger stencils rather than truncating them.
            if (nAdv <= 10) {
                int cell[10]; real positive[10],negative[10];
                for(int i=0;i<nAdv;++i){
                    cell[i]=advCellsForEdge(i,iEdge)-1;
                    const real c=adv_coefs(i,iEdge), c3=adv_coefs_3rd(i,iEdge);
                    positive[i]=c+c3; negative[i]=c-c3;
                }
                for(int sp=sb;sp<se;++sp){
                    #pragma omp simd
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);
                        const real* cf=uh>0 ? positive:negative;
                        real value=0.0;
                        if(nAdv==10){
                            // Preserve the original specialized summation order.
                            value=uh*(cf[0]*scalar_new_arr_scalars(k,cell[0],sp)+
                                      cf[1]*scalar_new_arr_scalars(k,cell[1],sp)+
                                      cf[2]*scalar_new_arr_scalars(k,cell[2],sp)+
                                      cf[3]*scalar_new_arr_scalars(k,cell[3],sp)+
                                      cf[4]*scalar_new_arr_scalars(k,cell[4],sp)+
                                      cf[5]*scalar_new_arr_scalars(k,cell[5],sp)+
                                      cf[6]*scalar_new_arr_scalars(k,cell[6],sp)+
                                      cf[7]*scalar_new_arr_scalars(k,cell[7],sp)+
                                      cf[8]*scalar_new_arr_scalars(k,cell[8],sp)+
                                      cf[9]*scalar_new_arr_scalars(k,cell[9],sp));
                        }else{
                            for(int i=0;i<nAdv;++i)value+=uh*cf[i]*scalar_new_arr_scalars(k,cell[i],sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }else{
                for(int sp=sb;sp<se;++sp){
                    #pragma omp simd
                    for(int k=0;k<nVertLevels_dummy;++k){
                        const real uh=uhAvg(k,iEdge);real value=0.0;
                        for(int i=0;i<nAdv;++i){
                            const real c=adv_coefs(i,iEdge),c3=adv_coefs_3rd(i,iEdge);
                            const real cf=uh>0 ? c+c3:c-c3;
                            value+=uh*cf*scalar_new_arr_scalars(k,advCellsForEdge(i,iEdge)-1,sp);
                        }
                        flux_arr_scalars(k,iEdge,sp)=value;
                    }
                }
            }
    });
    }
    } else {
    { PrecommStage stage(3);
    Kokkos::parallel_for("precomm_horiz_flux_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;

            if (cell1 > nCellsSolve - 1 && cell2 > nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                }
            } else {
                const int nAdv = nAdvCellsForEdge(iEdge);
                if (nAdv == 10) {
                    int ica[10];
                    real coef_pos[10];
                    real coef_neg[10];

                    for (int i = 0; i < 10; i++) {
                        ica[i] = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        coef_pos[i] = c + c3;
                        coef_neg[i] = c - c3;
                    }

                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        const real uh = uhAvg(k, iEdge);
                        const real* coef = (uh > 0.0) ? coef_pos : coef_neg;
                        flux_arr_scalars(k, iEdge, iScalar) = uh * (
                            coef[0]*scalar_new_arr_scalars(k, ica[0], iScalar) +
                            coef[1]*scalar_new_arr_scalars(k, ica[1], iScalar) +
                            coef[2]*scalar_new_arr_scalars(k, ica[2], iScalar) +
                            coef[3]*scalar_new_arr_scalars(k, ica[3], iScalar) +
                            coef[4]*scalar_new_arr_scalars(k, ica[4], iScalar) +
                            coef[5]*scalar_new_arr_scalars(k, ica[5], iScalar) +
                            coef[6]*scalar_new_arr_scalars(k, ica[6], iScalar) +
                            coef[7]*scalar_new_arr_scalars(k, ica[7], iScalar) +
                            coef[8]*scalar_new_arr_scalars(k, ica[8], iScalar) +
                            coef[9]*scalar_new_arr_scalars(k, ica[9], iScalar));
                    }
                } else {
                    for (int k = 0; k < nVertLevels_dummy; k++) {
                        flux_arr_scalars(k, iEdge, iScalar) = 0.0;
                    }
                    for (int i = 0; i < nAdv; i++) {
                        int iCell = advCellsForEdge(i, iEdge) - 1;
                        real c = adv_coefs(i, iEdge);
                        real c3 = adv_coefs_3rd(i, iEdge);
                        for (int k = 0; k < nVertLevels_dummy; k++) {
                            real uh = uhAvg(k, iEdge);
                            real coef = (uh > 0.0) ? (c + c3) : (c - c3);
                            flux_arr_scalars(k, iEdge, iScalar) +=
                                uh * coef * scalar_new_arr_scalars(k, iCell, iScalar);
                        }
                    }
                }
            }
    });
    }

    }
    { PrecommStage stage(4);
    Kokkos::parallel_for("precomm_vert_upwind_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars_dummy, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            wdtn_arr_scalars(0, iCell, iScalar) = 0.0;
            wdtn_arr_scalars(nVertLevels_dummy, iCell, iScalar) = 0.0;

            #pragma omp simd
            for (int k = 0; k < nVertLevels_dummy; k++) {
                real q = scalar_old_arr_scalars(k, iCell, iScalar) * rho_zz_old(k, iCell);
                real flux_in = 0.0;
                real flux_out = 0.0;
                if (k > 0) {
                    flux_in = vert_upwind_flux(dt, wwAvg(k, iCell),
                        scalar_old_arr_scalars(k-1, iCell, iScalar),
                        scalar_old_arr_scalars(k, iCell, iScalar));
                    wdtn_arr_scalars(k, iCell, iScalar) =
                        dt * wdtn_arr_scalars(k, iCell, iScalar) - flux_in;
                }
                if (k < nVertLevels_dummy - 1) {
                    flux_out = vert_upwind_flux(dt, wwAvg(k+1, iCell),
                        scalar_old_arr_scalars(k, iCell, iScalar),
                        scalar_old_arr_scalars(k+1, iCell, iScalar));
                }
                q += (flux_in - flux_out) * rdnw(k);
                scalar_new_arr_scalars(k,iCell,iScalar) = q;
            }

            #pragma omp simd
            for (int k = 0; k < nVertLevels_dummy; k++) {
                scale_arr_scalars(k, SCALE_IN, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::min(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::max(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
                scale_arr_scalars(k, SCALE_OUT, iScalar, iCell) =
                    -rdnw(k) *
                    (Kokkos::max(0.0, wdtn_arr_scalars(k+1, iCell, iScalar)) -
                     Kokkos::min(0.0, wdtn_arr_scalars(k, iCell, iScalar)));
            }
    });
    }

    { PrecommStage stage(5);
    Kokkos::parallel_for("precomm_horiz_lbc_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars_dummy, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                int cell1 = cellsOnEdge(0, iEdge) - 1;
                int cell2 = cellsOnEdge(1, iEdge) - 1;
                for (int k = 0; k < nVertLevels_dummy; k++) {
                    real flux_upwind =
                        dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar) +
                         Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    flux_arr_scalars(k, iEdge, iScalar) = flux_upwind;
                }
            }
    });
    }

    { PrecommStage stage(8);
    Kokkos::parallel_for("precomm_div_limiter_tiled",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,cellSolveStart-1}, {num_scalars_dummy,cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar,int iCell){
            const real inv_area=invAreaCell(iCell);
            const int neighbors=nEdgesOnCell(iCell);
            // A bounded, column-private tile retains k as the unit-stride loop.
            // No shared scratch, allocation, global padding or neighbor reorder.
            for(int kb=0;kb<nVertLevels_dummy;kb+=64){
                const int count=Kokkos::min(64,nVertLevels_dummy-kb);
                real q[64],incoming[64],outgoing[64],smax[64],smin[64];
                #pragma omp simd
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    q[t]=scalar_new_arr_scalars(k,iCell,iScalar);
                    incoming[t]=scale_arr_scalars(k,SCALE_IN,iScalar,iCell);
                    outgoing[t]=scale_arr_scalars(k,SCALE_OUT,iScalar,iCell);
                    const int below=k>0 ? k-1:k;
                    const int above=k<nVertLevels_dummy-1 ? k+1:k;
                    // Match endpoint and interior max/min operation order.
                    if(k==0){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(0,iCell,iScalar),scalar_old_arr_scalars(1,iCell,iScalar));
                    }else if(k==nVertLevels_dummy-1){
                        smax[t]=Kokkos::max(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                        smin[t]=Kokkos::min(scalar_old_arr_scalars(k,iCell,iScalar),scalar_old_arr_scalars(k-1,iCell,iScalar));
                    }else{
                        smax[t]=Kokkos::max(Kokkos::max(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                        smin[t]=Kokkos::min(Kokkos::min(scalar_old_arr_scalars(below,iCell,iScalar),scalar_old_arr_scalars(k,iCell,iScalar)),scalar_old_arr_scalars(above,iCell,iScalar));
                    }
                }
                for(int i=0;i<neighbors;++i){
                    const int edge=edgesOnCell(i,iCell)-1;
                    const int nbr=cellsOnCell(i,iCell)-1;
                    const real sign=edgesOnCell_sign(i,iCell);
                    const int c1=cellsOnEdge(0,edge)-1;
                    const int c2=cellsOnEdge(1,edge)-1;
                    const bool lbc=(config_apply_lbcs && bdyMaskEdge(edge)==nRelaxZone) || bdyMaskEdge(edge)==nRelaxZone-1;
                    const real edge_dt=dvEdge(edge)*dt;
                    #pragma omp simd
                    for(int t=0;t<count;++t){
                        const int k=kb+t;
                        const real uh=uhAvg(k,edge);
                        const real upwind=edge_dt*(Kokkos::max(0.0,uh)*scalar_old_arr_scalars(k,c1,iScalar)+Kokkos::min(0.0,uh)*scalar_old_arr_scalars(k,c2,iScalar));
                        const real correction=lbc ? 0.0:dt*flux_arr_scalars(k,edge,iScalar)-upwind;
                        q[t]-=sign*upwind*inv_area;
                        outgoing[t]-=Kokkos::max(0.0,sign*correction)*inv_area;
                        incoming[t]-=Kokkos::min(0.0,sign*correction)*inv_area;
                        const real neighbor=scalar_old_arr_scalars(k,nbr,iScalar);
                        smax[t]=Kokkos::max(smax[t],neighbor);
                        smin[t]=Kokkos::min(smin[t],neighbor);
                    }
                }
                #pragma omp simd
                for(int t=0;t<count;++t){
                    const int k=kb+t;
                    const real rho=local_advance_density ? rho_zz_int(k,iCell):rho_zz_new(k,iCell);
                    const real si=(smax[t]*rho-q[t])/(incoming[t]+eps);
                    const real so=(smin[t]*rho-q[t])/(outgoing[t]-eps);
                    scalar_new_arr_scalars(k,iCell,iScalar)=q[t];
                    scale_arr_scalars(k,SCALE_IN,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,si));
                    scale_arr_scalars(k,SCALE_OUT,iScalar,iCell)=Kokkos::min(1.0,Kokkos::max(0.0,so));
                }
            }
    });
    }
    Kokkos::fence();
}

extern "C" void advance_scalars_mono_precomm_kokkos(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels_dummy, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars_dummy, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p,
		flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p,
		flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p, flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, flcl_ndarray_t* rho_zz_int_p){
    static const int mode=[](){
        const char* text=std::getenv("IAMAS_PRECOMM_MODE");
        if(!text || !*text)return 0;
        char* end=nullptr;const long v=std::strtol(text,&end,10);
        if(*end || v<0 || v>4){std::cerr<<"Invalid IAMAS_PRECOMM_MODE; expected 0..4\n";std::abort();}
        return static_cast<int>(v);
    }();
    if(mode==0)precomm_baseline(maxEdges,nCells,nEdges,nVertLevels_dummy,dt,cellStart,cellEnd,vertexStart,vertexEnd,edgeStart,edgeEnd,cellSolveStart,cellSolveEnd,vertexSolveStart,vertexSolveEnd,edgeSolveStart,edgeSolveEnd,coef_3rd_order,nCellsSolve,num_scalars_dummy,uhAvg_p,wwAvg_p,scalar_tend_p,rho_zz_old_p,rho_zz_new_p,scalars_old_p,scalars_new_p,invAreaCell_p,dvEdge_p,cellsOnEdge_p,cellsOnCell_p,edgesOnCell_p,edgesOnCell_sign_p,nEdgesOnCell_p,fnm_p,fnp_p,rdnw_p,nAdvCellsForEdge_p,advCellsForEdge_p,adv_coefs_p,adv_coefs_3rd_p,scale_arr_p,scalar_old_p,scalar_new_p,s_max_p,s_min_p,wdtn_arr_p,flux_arr_p,flux_upwind_tmp_p,flux_tmp_p,scale_array_scalars_p,scalar_new_arr_scalars_p,scalar_old_arr_scalars_p,flux_array_scalars_p,wdtn_array_scalars_p,bdyMaskCell_p,bdyMaskEdge_p,advance_density,config_apply_lbcs,rho_zz_int_p);
    else if(mode==1)precomm_fused(maxEdges,nCells,nEdges,nVertLevels_dummy,dt,cellStart,cellEnd,vertexStart,vertexEnd,edgeStart,edgeEnd,cellSolveStart,cellSolveEnd,vertexSolveStart,vertexSolveEnd,edgeSolveStart,edgeSolveEnd,coef_3rd_order,nCellsSolve,num_scalars_dummy,uhAvg_p,wwAvg_p,scalar_tend_p,rho_zz_old_p,rho_zz_new_p,scalars_old_p,scalars_new_p,invAreaCell_p,dvEdge_p,cellsOnEdge_p,cellsOnCell_p,edgesOnCell_p,edgesOnCell_sign_p,nEdgesOnCell_p,fnm_p,fnp_p,rdnw_p,nAdvCellsForEdge_p,advCellsForEdge_p,adv_coefs_p,adv_coefs_3rd_p,scale_arr_p,scalar_old_p,scalar_new_p,s_max_p,s_min_p,wdtn_arr_p,flux_arr_p,flux_upwind_tmp_p,flux_tmp_p,scale_array_scalars_p,scalar_new_arr_scalars_p,scalar_old_arr_scalars_p,flux_array_scalars_p,wdtn_array_scalars_p,bdyMaskCell_p,bdyMaskEdge_p,advance_density,config_apply_lbcs,rho_zz_int_p);
    else if(mode==2)precomm_reuse(maxEdges,nCells,nEdges,nVertLevels_dummy,dt,cellStart,cellEnd,vertexStart,vertexEnd,edgeStart,edgeEnd,cellSolveStart,cellSolveEnd,vertexSolveStart,vertexSolveEnd,edgeSolveStart,edgeSolveEnd,coef_3rd_order,nCellsSolve,num_scalars_dummy,uhAvg_p,wwAvg_p,scalar_tend_p,rho_zz_old_p,rho_zz_new_p,scalars_old_p,scalars_new_p,invAreaCell_p,dvEdge_p,cellsOnEdge_p,cellsOnCell_p,edgesOnCell_p,edgesOnCell_sign_p,nEdgesOnCell_p,fnm_p,fnp_p,rdnw_p,nAdvCellsForEdge_p,advCellsForEdge_p,adv_coefs_p,adv_coefs_3rd_p,scale_arr_p,scalar_old_p,scalar_new_p,s_max_p,s_min_p,wdtn_arr_p,flux_arr_p,flux_upwind_tmp_p,flux_tmp_p,scale_array_scalars_p,scalar_new_arr_scalars_p,scalar_old_arr_scalars_p,flux_array_scalars_p,wdtn_array_scalars_p,bdyMaskCell_p,bdyMaskEdge_p,advance_density,config_apply_lbcs,rho_zz_int_p);
    else if(mode==3)precomm_simd(maxEdges,nCells,nEdges,nVertLevels_dummy,dt,cellStart,cellEnd,vertexStart,vertexEnd,edgeStart,edgeEnd,cellSolveStart,cellSolveEnd,vertexSolveStart,vertexSolveEnd,edgeSolveStart,edgeSolveEnd,coef_3rd_order,nCellsSolve,num_scalars_dummy,uhAvg_p,wwAvg_p,scalar_tend_p,rho_zz_old_p,rho_zz_new_p,scalars_old_p,scalars_new_p,invAreaCell_p,dvEdge_p,cellsOnEdge_p,cellsOnCell_p,edgesOnCell_p,edgesOnCell_sign_p,nEdgesOnCell_p,fnm_p,fnp_p,rdnw_p,nAdvCellsForEdge_p,advCellsForEdge_p,adv_coefs_p,adv_coefs_3rd_p,scale_arr_p,scalar_old_p,scalar_new_p,s_max_p,s_min_p,wdtn_arr_p,flux_arr_p,flux_upwind_tmp_p,flux_tmp_p,scale_array_scalars_p,scalar_new_arr_scalars_p,scalar_old_arr_scalars_p,flux_array_scalars_p,wdtn_array_scalars_p,bdyMaskCell_p,bdyMaskEdge_p,advance_density,config_apply_lbcs,rho_zz_int_p);
    else if(mode==4)precomm_fused_simd(maxEdges,nCells,nEdges,nVertLevels_dummy,dt,cellStart,cellEnd,vertexStart,vertexEnd,edgeStart,edgeEnd,cellSolveStart,cellSolveEnd,vertexSolveStart,vertexSolveEnd,edgeSolveStart,edgeSolveEnd,coef_3rd_order,nCellsSolve,num_scalars_dummy,uhAvg_p,wwAvg_p,scalar_tend_p,rho_zz_old_p,rho_zz_new_p,scalars_old_p,scalars_new_p,invAreaCell_p,dvEdge_p,cellsOnEdge_p,cellsOnCell_p,edgesOnCell_p,edgesOnCell_sign_p,nEdgesOnCell_p,fnm_p,fnp_p,rdnw_p,nAdvCellsForEdge_p,advCellsForEdge_p,adv_coefs_p,adv_coefs_3rd_p,scale_arr_p,scalar_old_p,scalar_new_p,s_max_p,s_min_p,wdtn_arr_p,flux_arr_p,flux_upwind_tmp_p,flux_tmp_p,scale_array_scalars_p,scalar_new_arr_scalars_p,scalar_old_arr_scalars_p,flux_array_scalars_p,wdtn_array_scalars_p,bdyMaskCell_p,bdyMaskEdge_p,advance_density,config_apply_lbcs,rho_zz_int_p);
}
