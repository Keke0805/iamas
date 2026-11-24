/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop advance_scalars_mono_afcomm subroutine in kokkos

    * Author: Xin Qu
    Last update: June 2026

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


extern "C" void advance_scalars_mono_afcomm_kokkos(int &maxEdges, int &nCells, int &nEdges, int &nVertLevels, real &dt,
		int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
		int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
		real &coef_3rd_order, int &nCellsSolve, int &num_scalars, flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_old_p, 
		flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* scalars_old_p, flcl_ndarray_t* scalars_new_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnCell_p,
		flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* nAdvCellsForEdge_p,
		flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* scale_arr_p,flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p, 							
        flcl_ndarray_t* s_max_p, flcl_ndarray_t* s_min_p, flcl_ndarray_t* wdtn_arr_p, flcl_ndarray_t* flux_arr_p, flcl_ndarray_t* flux_upwind_tmp_p, flcl_ndarray_t* flux_tmp_p, 
        flcl_ndarray_t* scale_array_scalars_p, flcl_ndarray_t* scalar_new_arr_scalars_p, flcl_ndarray_t* scalar_old_arr_scalars_p, flcl_ndarray_t* flux_array_scalars_p,flcl_ndarray_t* wdtn_array_scalars_p,
		flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
	    bool &advance_density, bool &config_apply_lbcs, bool &config_adv_mono_comm, flcl_ndarray_t* rho_zz_int_p){
	auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
	auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
	auto scalars_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_new_p);
	auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
	auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
	auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
	auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
	auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
	auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
	auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
    auto scale_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*scale_array_scalars_p);
    auto scalar_new_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_arr_scalars_p);
    auto scalar_old_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_arr_scalars_p);
    auto flux_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*flux_array_scalars_p);
    auto wdtn_arr_scalars = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_array_scalars_p);
	auto bdyMaskCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskCell_p);
	auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);

    const int SCALE_IN = 0, SCALE_OUT = 1;
    const int nSpecZone = 2;
    const int nRelaxZone = 5;
    bool local_advance_density = advance_density;
    // Launch each algorithmic stage once over all scalars instead of launching
    // every stage separately for every scalar.
    Kokkos::parallel_for("afcomm_flux_upwind_scale_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, edgeStart - 1}, {num_scalars, edgeEnd}),
        KOKKOS_LAMBDA(int iScalar, int iEdge) {
            int cell1 = cellsOnEdge(0, iEdge) - 1;
            int cell2 = cellsOnEdge(1, iEdge) - 1;
            if (cell1 <= nCellsSolve - 1 || cell2 <= nCellsSolve - 1) {
                for (int k = 0; k < nVertLevels; k++) {
                    real flux_upwind = dvEdge(iEdge) * dt *
                        (Kokkos::max(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell1, iScalar)
                       + Kokkos::min(0.0, uhAvg(k, iEdge)) *
                             scalar_old_arr_scalars(k, cell2, iScalar));
                    real flux = dt * flux_arr_scalars(k, iEdge, iScalar) - flux_upwind;
                    if ((config_apply_lbcs && (bdyMaskEdge(iEdge) == nRelaxZone)) ||
                        (bdyMaskEdge(iEdge) == nRelaxZone - 1)) {
                        flux = 0.0;
                    } else if (config_adv_mono_comm ||
                               (cell1 <= nCellsSolve - 1 && cell2 <= nCellsSolve - 1)) {
                        flux =
                            Kokkos::max(0.0, flux) *
                                Kokkos::min(scale_arr_scalars(k, SCALE_OUT, iScalar, cell1),
                                            scale_arr_scalars(k, SCALE_IN,  iScalar, cell2))
                          + Kokkos::min(0.0, flux) *
                                Kokkos::min(scale_arr_scalars(k, SCALE_IN,  iScalar, cell1),
                                            scale_arr_scalars(k, SCALE_OUT, iScalar, cell2));
                    }
                    flux_arr_scalars(k, iEdge, iScalar) = flux;
                }
            }
    });

    Kokkos::parallel_for("afcomm_scale_wdtn_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            for (int k = 1; k < nVertLevels; k++) {
                real flux = wdtn_arr_scalars(k, iCell, iScalar);
                flux =
                    Kokkos::max(0.0, flux) *
                        Kokkos::min(scale_arr_scalars(k-1, SCALE_OUT, iScalar, iCell),
                                    scale_arr_scalars(k,   SCALE_IN,  iScalar, iCell))
                  + Kokkos::min(0.0, flux) *
                        Kokkos::min(scale_arr_scalars(k,   SCALE_OUT, iScalar, iCell),
                                    scale_arr_scalars(k-1, SCALE_IN,  iScalar, iCell));
                wdtn_arr_scalars(k, iCell, iScalar) = flux;
            }
    });

    Kokkos::parallel_for("afcomm_update_scalar_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellSolveStart - 1}, {num_scalars, cellSolveEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            for (int i = 0; i < nEdgesOnCell(iCell); i++) {
                int iEdge = edgesOnCell(i, iCell) - 1;
                real sign = edgesOnCell_sign(i, iCell);
                for (int k = 0; k < nVertLevels; k++) {
                    scalar_new_arr_scalars(k, iCell, iScalar) -=
                        sign * flux_arr_scalars(k, iEdge, iScalar) * invAreaCell(iCell);
                }
            }
            for (int k = 0; k < nVertLevels; k++) {
                real rho = local_advance_density ? rho_zz_int(k, iCell) : rho_zz_new(k, iCell);
                scalar_new_arr_scalars(k, iCell, iScalar) =
                    (scalar_new_arr_scalars(k, iCell, iScalar)
                    - rdnw(k) *
                        (wdtn_arr_scalars(k+1, iCell, iScalar) -
                         wdtn_arr_scalars(k, iCell, iScalar)))
                    / rho;
            }
    });

    Kokkos::parallel_for("afcomm_writeback_all",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, cellStart - 1}, {num_scalars, cellEnd}),
        KOKKOS_LAMBDA(int iScalar, int iCell) {
            if (bdyMaskCell(iCell) <= nSpecZone) {
                for (int k = 0; k < nVertLevels; k++) {
                    scalars_new(iScalar, k, iCell) =
                        Kokkos::max(0.0, scalar_new_arr_scalars(k, iCell, iScalar));
                }
            }
    });
    Kokkos::fence();
}
