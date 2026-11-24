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
    Kokkos::parallel_for("precomm_init_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, cellStart - 1, 0}, {num_scalars_dummy, cellEnd, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int iCell, int k) {
            scalar_old_arr_scalars(k, iCell, iScalar) = scalars_old(iScalar, k, iCell);
            scalar_new_arr_scalars(k, iCell, iScalar) = scalars_new(iScalar, k, iCell);
    });

    Kokkos::parallel_for("precomm_init_ghost_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_scalars_dummy, nVertLevels_dummy}),
        KOKKOS_LAMBDA(int iScalar, int k) {
            scalar_old_arr_scalars(k, nCells, iScalar) = 0.0;
            scalar_new_arr_scalars(k, nCells, iScalar) = 0.0;
    });

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
    Kokkos::fence();
}
