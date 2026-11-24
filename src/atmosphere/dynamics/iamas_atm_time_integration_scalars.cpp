/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop advance_scalars subroutine in kokkos

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

KOKKOS_INLINE_FUNCTION real flux4(const real q_im2,const real q_im1, const real q_i, const real q_ip1, const real ua){
	return ua * (7. * (q_i + q_im1) - (q_ip1 + q_im2) ) / 12.0;
}

KOKKOS_INLINE_FUNCTION real flux3(const real q_im2,const real q_im1, const real q_i, const real q_ip1, const real ua, const real coef3){
	return flux4(q_im2, q_im1, q_i, q_ip1, ua) + coef3 * Kokkos::abs(ua) * ((q_ip1 - q_im2) - 3. * (q_i - q_im1)) / 12.0;
}

extern "C" void advance_scalars_kokkos(int &maxEdges,int &num_scalars_dummy, int &nCells, int &nVertLevels_dummy, real &dt,
						  int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd,
						  int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd,
						  real &coef_3rd_order, flcl_ndarray_t* scalar_old_p, flcl_ndarray_t* scalar_new_p, flcl_ndarray_t* rho_zz_old_p, flcl_ndarray_t* rho_zz_new_p, flcl_ndarray_t* kdiff_p,
						  flcl_ndarray_t* uhAvg_p, flcl_ndarray_t* wwAvg_p, flcl_ndarray_t* deriv_two_p, flcl_ndarray_t* dvEdge_p, 
						  flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* invAreaCell_p,
						  flcl_ndarray_t* scalar_tend_save_p, flcl_ndarray_t* fnm_p, flcl_ndarray_t* fnp_p, flcl_ndarray_t* rdnw_p, flcl_ndarray_t* meshScalingDel2_p, flcl_ndarray_t* meshScalingDel4_p,
						  flcl_ndarray_t* bdyMaskCell_p, flcl_ndarray_t* bdyMaskEdge_p,
						  flcl_ndarray_t* nAdvCellsForEdge_p, flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, flcl_ndarray_t* rho_edge_p, flcl_ndarray_t* qv_init_p,flcl_ndarray_t* zgrid_p,
						  int &nCellsSolve, int &nEdges, flcl_ndarray_t* horiz_flux_arr_p, int &rk_step, int &config_time_integration_order,
						  bool &advance_density, bool &config_apply_lbcs,flcl_ndarray_t* wdtn_scalar_p, flcl_ndarray_t* scalar_tend_column_p, flcl_ndarray_t* scalar_tend_p, flcl_ndarray_t* rho_zz_int_p){
    auto dvEdge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p);
    auto qv_init = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*qv_init_p);
    auto fnm = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnm_p);
    auto fnp = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fnp_p);
    auto rdnw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdnw_p);
    auto meshScalingDel2 = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*meshScalingDel2_p);
    auto meshScalingDel4 = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*meshScalingDel4_p);
    auto invAreaCell = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p);
    auto rho_zz_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_old_p);
    auto rho_zz_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_new_p);
    auto uhAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*uhAvg_p);
    auto wwAvg = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p);
    auto kdiff = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kdiff_p);
    auto rho_edge = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_edge_p);
    auto zgrid = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zgrid_p);
    auto adv_coefs = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p);
    auto adv_coefs_3rd = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p);
    auto edgesOnCell_sign = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p);
    auto rho_zz_int = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_int_p);
    auto scalar_old = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_old_p);
    auto scalar_new = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_new_p);
    auto scalar_tend_save = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_tend_save_p);
    auto deriv_two = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*deriv_two_p);
    auto horiz_flux_arr = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*horiz_flux_arr_p);
    auto scalar_tend = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_tend_p);
    auto cellsOnEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p);
    auto edgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p);
    auto advCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p);
    auto nAdvCellsForEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p);
    auto nEdgesOnCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p);
    auto bdyMaskCell = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskCell_p);
    auto bdyMaskEdge = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*bdyMaskEdge_p);
	auto wdtn_scalar = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*wdtn_scalar_p);
	auto scalar_tend_column = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalar_tend_column_p);
    bool local_advance_density;
    real weight_time_old, weight_time_new;
    local_advance_density = advance_density;
    const int nRelaxZone = 5;
	/*
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !
      !  horizontal flux divergence, accumulate in scalar_tend


      !  horiz_flux_arr stores the value of the scalar at the edge.
      !  a better name perhaps would be scalarEdge

      !  weights for the time interpolation of the input density
      !
	**/
    if(!advance_density) {
        weight_time_new = 1.;
    }else{
        if((rk_step == 1) && config_time_integration_order == 3) weight_time_new = 1./3;
        if((rk_step == 1) && config_time_integration_order == 2) weight_time_new = 1./2;
        if(rk_step == 2) weight_time_new = 1./2;
        if(rk_step == 3) weight_time_new = 1.;
    }
    weight_time_old = 1. - weight_time_new;
    for(int iEdge = edgeStart - 1; iEdge < edgeEnd; iEdge++){
        for(int k = 0; k < nVertLevels_dummy; k++){
            for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
                horiz_flux_arr(iScalar,k,iEdge) = 0.0;
            }
        }
    }
    Kokkos::parallel_for(
        "Edges_kokkos", Kokkos::RangePolicy<>(edgeStart - 1, edgeEnd), KOKKOS_LAMBDA(int iEdge) {
        if((!config_apply_lbcs) || (bdyMaskEdge(iEdge) < nRelaxZone - 1)){
			for (int j = 0; j < nAdvCellsForEdge(iEdge); j++){
			    int iAdvCell = advCellsForEdge(j,iEdge) - 1;
			    for(int k = 0; k < nVertLevels_dummy; k++){
					real scalar_weight = adv_coefs(j,iEdge) + copysign(1.0, uhAvg(k,iEdge)) * adv_coefs_3rd(j,iEdge);
					for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
						horiz_flux_arr(iScalar,k,iEdge) += scalar_weight * scalar_new(iScalar,k,iAdvCell);
					}
			    }
			}
        }else if(config_apply_lbcs && (bdyMaskEdge(iEdge) >= nRelaxZone-1) && (bdyMaskEdge(iEdge) <= nRelaxZone)){
            int cell1 = cellsOnEdge(0,iEdge) - 1;
            int cell2 = cellsOnEdge(1,iEdge) - 1;
            for(int k = 0; k < nVertLevels_dummy; k++){
                real u_direction = copysign(0.5,uhAvg(k,iEdge));
                real u_positive = dvEdge(iEdge)*abs(u_direction + 0.5);
                real u_negative = dvEdge(iEdge)*abs(u_direction - 0.5);
                for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
                    horiz_flux_arr(iScalar,k,iEdge) = u_positive*scalar_new(iScalar,k,cell1) + u_negative*scalar_new(iScalar,k,cell2);
                }
            }   
        }
    });
	Kokkos::fence();
	// scalar update, for each column sum fluxes over horizontal edges, add physics tendency, and add vertical flux divergence in update.
	Kokkos::parallel_for(
        "Cells_kokkos", Kokkos::RangePolicy<>(cellSolveStart - 1, cellSolveEnd), KOKKOS_LAMBDA(int iCell) {
		if(bdyMaskCell(iCell) <= nRelaxZone){ // specified zone for regional_iAMAS is not updated in this routine
			 for(int k = 0; k < nVertLevels_dummy; k++){
				for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
					scalar_tend_column(iScalar,k,iCell) = 0.0;
				}
			 }
			for(int i = 0; i < nEdgesOnCell(iCell); i++){
				int iEdge = edgesOnCell(i,iCell) - 1;
				//! here we add the horizontal flux divergence into the scalar tendency.
                		//! note that the scalar tendency is modified.
				for(int k = 0; k < nVertLevels_dummy; k++){
					for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
						scalar_tend_column(iScalar,k,iCell) = scalar_tend_column(iScalar,k,iCell) - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*horiz_flux_arr(iScalar,k,iEdge);
					}
				}
			}
			
			for(int k = 0; k < nVertLevels_dummy; k++){
				for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
					scalar_tend_column(iScalar,k,iCell) = scalar_tend_column(iScalar,k,iCell) * invAreaCell(iCell) + scalar_tend_save(iScalar,k,iCell);
				}
			}
			
			// vertical flux divergence and update of the scalars
			for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
				wdtn_scalar(iScalar,0,iCell) = 0.0;
				wdtn_scalar(iScalar,nVertLevels_dummy,iCell) = 0.0;
			}
			int k = 1;
			for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
				wdtn_scalar(iScalar,k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell));
			}
			
			for(int k = 2; k < nVertLevels_dummy-1; k++){
				for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
					wdtn_scalar(iScalar,k,iCell) = flux3( scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  
						scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell),  
						wwAvg(k,iCell), coef_3rd_order );
				}
			}
			
			
			k = nVertLevels_dummy - 1;
			for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
				wdtn_scalar(iScalar,k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell));
			}
			for(int k = 0; k < nVertLevels_dummy; k++){
				real rho_zz_new_inv = 1.0 / (weight_time_old * rho_zz_old(k,iCell) + weight_time_new * rho_zz_new(k,iCell));
				for(int iScalar = 0; iScalar < num_scalars_dummy; iScalar++){
					scalar_new(iScalar,k,iCell) = (scalar_old(iScalar,k,iCell) * rho_zz_old(k,iCell) + dt * (scalar_tend_column(iScalar,k,iCell) - rdnw(k) * (wdtn_scalar(iScalar,k+1,iCell) - wdtn_scalar(iScalar,k,iCell)))) * rho_zz_new_inv;
				}
			}
		}
	});
	Kokkos::fence();
}




