/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop dyn_tend subroutine in kokkos

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

KOKKOS_INLINE_FUNCTION real flux4(const real &q_im2,const real &q_im1, const real &q_i, const real &q_ip1, const real &ua){
	return ua * (7. * (q_i + q_im1) - (q_ip1 + q_im2) ) / 12.0;
}

KOKKOS_INLINE_FUNCTION real flux3(const real &q_im2,const real &q_im1, const real &q_i, const real &q_ip1, const real &ua, const real &coef3){
	return flux4(q_im2, q_im1, q_i, q_ip1, ua) + coef3 * Kokkos::abs(ua) * ((q_ip1 - q_im2) - 3. * (q_i - q_im1)) / 12.0;
}

extern "C" void atm_compute_dyn_tend_kokkos(int &nCells, int &nEdges, int &nVertices, int &nVertLevels, 
    int &nCellsSolve, int &nEdgesSolve, int &vertexDegree, int &maxEdges, int &maxEdges2, int &num_scalars, int &moist_start, int &moist_end, 
    flcl_ndarray_t* fEdge_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* dcEdge_p, flcl_ndarray_t* invDcEdge_p, flcl_ndarray_t* invDvEdge_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* invAreaTriangle_p, flcl_ndarray_t* meshScalingDel2_p, flcl_ndarray_t* meshScalingDel4_p, 
    flcl_ndarray_t* weightsOnEdge_p, flcl_ndarray_t* zgrid_p, flcl_ndarray_t* rho_edge_p, flcl_ndarray_t* rho_zz_p, flcl_ndarray_t* ru_p, flcl_ndarray_t* u_p, flcl_ndarray_t* v_p, flcl_ndarray_t* tend_u_p, 
    flcl_ndarray_t* divergence_p, flcl_ndarray_t* vorticity_p, flcl_ndarray_t* ke_p, flcl_ndarray_t* pv_edge_p, flcl_ndarray_t* theta_m_p, flcl_ndarray_t* rw_p, flcl_ndarray_t* tend_rho_p, 
    flcl_ndarray_t* rt_diabatic_tend_p, flcl_ndarray_t* tend_theta_p, flcl_ndarray_t* tend_w_p, flcl_ndarray_t* w_p, flcl_ndarray_t* cqw_p, flcl_ndarray_t *rb_p, flcl_ndarray_t* rr_p, flcl_ndarray_t* pp_p, flcl_ndarray_t* pressure_b_p, flcl_ndarray_t* phi_p_p, flcl_ndarray_t* phi_base_p, flcl_ndarray_t* zz_p, flcl_ndarray_t* zxu_p, flcl_ndarray_t* cqu_p, 
    flcl_ndarray_t* h_divergence_p, flcl_ndarray_t* kdiff_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* edgesOnVertex_sign_p, flcl_ndarray_t* rw_save_p, flcl_ndarray_t* ru_save_p, 
    flcl_ndarray_t* theta_m_save_p, flcl_ndarray_t* exner_p, flcl_ndarray_t* rr_save_p, flcl_ndarray_t* scalars_p, flcl_ndarray_t* tend_u_euler_p, flcl_ndarray_t* tend_w_euler_p, flcl_ndarray_t* tend_theta_euler_p, flcl_ndarray_t* deriv_two_p, 
    flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* verticesOnEdge_p, flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnEdge_p, flcl_ndarray_t* cellsOnCell_p, flcl_ndarray_t* edgesOnVertex_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* nEdgesOnEdge_p, 
    flcl_ndarray_t* latCell_p, flcl_ndarray_t* latEdge_p, flcl_ndarray_t* angleEdge_p, flcl_ndarray_t* u_init_p, flcl_ndarray_t* v_init_p, flcl_ndarray_t* advCellsForEdge_p, flcl_ndarray_t* nAdvCellsForEdge_p, flcl_ndarray_t* adv_coefs_p, flcl_ndarray_t* adv_coefs_3rd_p, 
    flcl_ndarray_t* rdzu_p, flcl_ndarray_t* rdzw_p, flcl_ndarray_t* fzm_p, flcl_ndarray_t* fzp_p, flcl_ndarray_t* qv_init_p, flcl_ndarray_t* t_init_p, real& cf1, real& cf2, real& cf3, real& r_earth, flcl_ndarray_t* ur_cell_p, flcl_ndarray_t* vr_cell_p, flcl_ndarray_t* defc_a_p, flcl_ndarray_t* defc_b_p, 
    flcl_ndarray_t* tend_w_pgf_p, flcl_ndarray_t* tend_w_buoy_p, real& coef_3rd_order, real &c_s, bool &config_mix_full, char *config_horiz_mixing, real& config_del4u_div_factor, 
    real& config_h_mom_eddy_visc2, real& config_v_mom_eddy_visc2, real& config_h_theta_eddy_visc2, real& config_v_theta_eddy_visc2, 
    real& config_h_theta_eddy_visc4, real& config_h_mom_eddy_visc4, real& config_visc4_2dsmag, real& config_len_disp, int &rk_step, real& dt, 
    real& config_iamas_cam_coef, bool &config_include_kinetic_visc,  
    bool &config_rayleigh_damp_u, real &config_rayleigh_damp_u_timescale_days, int &config_number_rayleigh_damp_u_levels, 
    flcl_ndarray_t* tend_rtheta_adv_p, flcl_ndarray_t* rthdynten_p, flcl_ndarray_t* qtot_p, flcl_ndarray_t* tend_rtheta_physics_p, flcl_ndarray_t* tend_rho_physics_p, flcl_ndarray_t* tend_ru_physics_p,
    flcl_ndarray_t* delsq_theta_p, flcl_ndarray_t* delsq_w_p, flcl_ndarray_t* delsq_divergence_p, flcl_ndarray_t* delsq_u_p, flcl_ndarray_t* delsq_vorticity_p, flcl_ndarray_t* dpdz_p,
    int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd, 
    int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd){
    // std::cout << strcmp(config_horiz_mixing, "2d_smagorinsky") << std::endl;
    auto fEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fEdge_p); auto fEdge = real1d(fEdge_raw,{1});
    auto dvEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p); auto dvEdge = real1d(dvEdge_raw,{1});
    auto dcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dcEdge_p); auto dcEdge = real1d(dcEdge_raw,{1});
    auto invDcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invDcEdge_p); auto invDcEdge = real1d(invDcEdge_raw,{1});
    auto invDvEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invDvEdge_p); auto invDvEdge = real1d(invDvEdge_raw,{1});
    auto invAreaCell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p); auto invAreaCell = real1d(invAreaCell_raw,{1});
    auto invAreaTriangle_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaTriangle_p); auto invAreaTriangle = real1d(invAreaTriangle_raw,{1});
    auto meshScalingDel2_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*meshScalingDel2_p); auto meshScalingDel2 = real1d(meshScalingDel2_raw,{1});
    auto meshScalingDel4_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*meshScalingDel4_p); auto meshScalingDel4 = real1d(meshScalingDel4_raw,{1});
    auto weightsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*weightsOnEdge_p); auto weightsOnEdge = real2d(weightsOnEdge_raw,{1,1});
    auto zgrid_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zgrid_p); auto zgrid = real2d(zgrid_raw,{1,1});
    auto rho_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_edge_p); auto rho_edge = real2d(rho_edge_raw,{1,1});
    auto rho_zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_p); auto rho_zz = real2d(rho_zz_raw,{1,1});
    auto ru_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_p); auto ru = real2d(ru_raw,{1,1});
    auto u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*u_p); auto u = real2d(u_raw,{1,1});
    auto v_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*v_p); auto v = real2d(v_raw,{1,1});
    auto tend_u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_u_p); auto tend_u = real2d(tend_u_raw,{1,1});
    auto divergence_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*divergence_p); auto divergence = real2d(divergence_raw,{1,1});
    auto vorticity_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*vorticity_p); auto vorticity = real2d(vorticity_raw,{1,1});
    auto ke_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ke_p); auto ke = real2d(ke_raw,{1,1});
    auto pv_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pv_edge_p); auto pv_edge = real2d(pv_edge_raw,{1,1});
    auto theta_m_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*theta_m_p); auto theta_m = real2d(theta_m_raw,{1,1});
    auto rw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_p); auto rw = real2d(rw_raw,{1,1});
    auto tend_rho_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rho_p); auto tend_rho = real2d(tend_rho_raw,{1,1});
    auto rt_diabatic_tend_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rt_diabatic_tend_p); auto rt_diabatic_tend = real2d(rt_diabatic_tend_raw,{1,1});
    auto tend_theta_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_theta_p); auto tend_theta = real2d(tend_theta_raw,{1,1});
    auto tend_w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_w_p); auto tend_w = real2d(tend_w_raw,{1,1});
    auto w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*w_p); auto w = real2d(w_raw,{1,1});
    auto cqw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cqw_p); auto cqw = real2d(cqw_raw,{1,1});
    auto rb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rb_p); auto rb = real2d(rb_raw,{1,1});
    auto rr_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rr_p); auto rr = real2d(rr_raw,{1,1});
    auto pp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pp_p); auto pp = real2d(pp_raw,{1,1});
    auto pressure_b_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pressure_b_p); auto pressure_b = real2d(pressure_b_raw,{1,1});
    auto phi_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*phi_p_p); auto phi_p = real2d(phi_p_raw,{1,1});
    auto phi_base_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*phi_base_p); auto phi_base = real2d(phi_base_raw,{1,1});
    auto zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zz_p); auto zz = real2d(zz_raw,{1,1});
    auto zxu_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zxu_p); auto zxu = real2d(zxu_raw,{1,1});
    auto cqu_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cqu_p); auto cqu = real2d(cqu_raw,{1,1});
    auto h_divergence_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*h_divergence_p); auto h_divergence = real2d(h_divergence_raw,{1,1});
    auto kdiff_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kdiff_p); auto kdiff = real2d(kdiff_raw,{1,1});
    auto edgesOnCell_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p); auto edgesOnCell_sign = real2d(edgesOnCell_sign_raw,{1,1});
    auto edgesOnVertex_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnVertex_sign_p); auto edgesOnVertex_sign = real2d(edgesOnVertex_sign_raw,{1,1});
    auto rw_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_save_p); auto rw_save = real2d(rw_save_raw,{1,1});
    auto ru_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_save_p); auto ru_save = real2d(ru_save_raw,{1,1});
    auto theta_m_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*theta_m_save_p); auto theta_m_save = real2d(theta_m_save_raw,{1,1});
    auto exner_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*exner_p); auto exner = real2d(exner_raw,{1,1});
    auto rr_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rr_save_p); auto rr_save = real2d(rr_save_raw,{1,1});
    auto scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_p); auto scalars = real3d(scalars_raw,{1,1,1});
    auto tend_u_euler_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_u_euler_p); auto tend_u_euler = real2d(tend_u_euler_raw,{1,1});
    auto tend_w_euler_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_w_euler_p); auto tend_w_euler = real2d(tend_w_euler_raw,{1,1});
    auto tend_theta_euler_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_theta_euler_p); auto tend_theta_euler = real2d(tend_theta_euler_raw,{1,1});
    auto deriv_two_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*deriv_two_p); auto deriv_two = real3d(deriv_two_raw,{1,1,1});
    auto cellsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p); auto cellsOnEdge = int2d(cellsOnEdge_raw,{1,1});
    auto verticesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*verticesOnEdge_p); auto verticesOnEdge = int2d(verticesOnEdge_raw,{1,1});
    auto edgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p); auto edgesOnCell = int2d(edgesOnCell_raw,{1,1});
    auto edgesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnEdge_p); auto edgesOnEdge = int2d(edgesOnEdge_raw,{1,1});
    auto cellsOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnCell_p); auto cellsOnCell = int2d(cellsOnCell_raw,{1,1});
    auto edgesOnVertex_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnVertex_p); auto edgesOnVertex = int2d(edgesOnVertex_raw,{1,1});
    auto nEdgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p); auto nEdgesOnCell = int1d(nEdgesOnCell_raw,{1});
    auto nEdgesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnEdge_p); auto nEdgesOnEdge = int1d(nEdgesOnEdge_raw,{1});
    auto latCell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*latCell_p); auto latCell = real1d(latCell_raw,{1});
    auto latEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*latEdge_p); auto latEdge = real1d(latEdge_raw,{1});
    auto angleEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*angleEdge_p); auto angleEdge = real1d(angleEdge_raw,{1});
    auto u_init_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*u_init_p); auto u_init = real1d(u_init_raw,{1});
    auto v_init_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*v_init_p); auto v_init = real1d(v_init_raw,{1});
    auto advCellsForEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*advCellsForEdge_p); auto advCellsForEdge = int2d(advCellsForEdge_raw,{1,1});
    auto nAdvCellsForEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nAdvCellsForEdge_p); auto nAdvCellsForEdge = int1d(nAdvCellsForEdge_raw,{1});
    auto adv_coefs_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_p); auto adv_coefs = real2d(adv_coefs_raw,{1,1});
    auto adv_coefs_3rd_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*adv_coefs_3rd_p); auto adv_coefs_3rd = real2d(adv_coefs_3rd_raw,{1,1});
    auto rdzu_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdzu_p); auto rdzu = real1d(rdzu_raw,{1});
    auto rdzw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdzw_p); auto rdzw = real1d(rdzw_raw,{1});
    auto fzm_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzm_p); auto fzm = real1d(fzm_raw,{1});
    auto fzp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzp_p); auto fzp = real1d(fzp_raw,{1});
    auto qv_init_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*qv_init_p); auto qv_init = real1d(qv_init_raw,{1});
    auto t_init_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*t_init_p); auto t_init = real2d(t_init_raw,{1,1});
    auto ur_cell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ur_cell_p); auto ur_cell = real2d(ur_cell_raw,{1,1});
    auto vr_cell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*vr_cell_p); auto vr_cell = real2d(vr_cell_raw,{1,1});
    auto defc_a_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*defc_a_p); auto defc_a = real2d(defc_a_raw,{1,1});
    auto defc_b_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*defc_b_p); auto defc_b = real2d(defc_b_raw,{1,1});
    auto tend_w_pgf_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_w_pgf_p); auto tend_w_pgf = real2d(tend_w_pgf_raw,{1,1});
    auto tend_w_buoy_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_w_buoy_p); auto tend_w_buoy = real2d(tend_w_buoy_raw,{1,1});
    auto tend_rtheta_adv_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rtheta_adv_p); auto tend_rtheta_adv = real2d(tend_rtheta_adv_raw,{1,1});
    auto rthdynten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rthdynten_p); auto rthdynten = real2d(rthdynten_raw,{1,1});
    auto qtot_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*qtot_p); auto qtot = real2d(qtot_raw,{1,1});
    auto tend_rtheta_physics_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rtheta_physics_p); auto tend_rtheta_physics = real2d(tend_rtheta_physics_raw,{1,1});
    auto tend_rho_physics_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rho_physics_p); auto tend_rho_physics = real2d(tend_rho_physics_raw,{1,1});
    auto tend_ru_physics_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_ru_physics_p); auto tend_ru_physics = real2d(tend_ru_physics_raw,{1,1});
    auto delsq_theta_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*delsq_theta_p); auto delsq_theta = real2d(delsq_theta_raw,{1,1});
    auto delsq_w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*delsq_w_p); auto delsq_w = real2d(delsq_w_raw,{1,1});
    auto delsq_divergence_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*delsq_divergence_p); auto delsq_divergence = real2d(delsq_divergence_raw,{1,1});
    auto delsq_u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*delsq_u_p); auto delsq_u = real2d(delsq_u_raw,{1,1});
    auto delsq_vorticity_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*delsq_vorticity_p); auto delsq_vorticity = real2d(delsq_vorticity_raw,{1,1});
    auto dpdz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*dpdz_p); auto dpdz = real2d(dpdz_raw,{1,1});
    real rayleigh_damp_coef_s[nVertLevels];
    real kin_visc_s[nVertLevels * (nCells + 1)];
    auto rayleigh_damp_coef = real1d((real*)rayleigh_damp_coef_s, {1},{nVertLevels + 1});
    auto kin_visc = real2d((real*)kin_visc_s, {1,1},{nVertLevels+1,nCells+2});
    const real prandtl = 1.0;
    const real omega = 7.29212e-5;
    const real gravity = 9.80616;
    const real seconds_per_day = 86400.0;
    real prandtl_inv = 1.0 / prandtl;
    real invDt = 1.0 / dt;
    real inv_r_earth = 1.0 / r_earth;
    real v_mom_eddy_visc2   = config_v_mom_eddy_visc2;
    real v_theta_eddy_visc2 = config_v_theta_eddy_visc2;
    const real rgas = 287.0;
    const real cp = 7. * rgas / 2.0;
    const real cv = cp - rgas;
    const real rv = 461.6;
    real h_mom_eddy_visc4,h_theta_eddy_visc4;
    if(rk_step == 1){
        // jwfeng: Kinetic viscosity, T = theta_m*exner*/((1+rv/rd)*qv)
        // Kokkos::parallel_for(
        //     "Cell_Kokkos",Kokkos::RangePolicy<>(1, nCells + 2), KOKKOS_LAMBDA(int iCell) {
        for(int iCell = 1; iCell <= nCells + 1; iCell++){
            for(int k = 1; k <= nVertLevels; k++){
                kin_visc(k,iCell) = 0.0;
            }
        // });
        }
        Kokkos::fence();
        if(config_include_kinetic_visc){
            Kokkos::parallel_for(
               "Cell_Kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
                for(int k = 1; k <= nVertLevels; k++){
                    kin_visc(k,iCell) = (4.03*0.21 + 3.42*0.78 + 0.01*3.9) * 1.e-7 * 
                     Kokkos::pow(( exner(k,iCell)*theta_m(k,iCell)/(1.0+rv/rgas)/cqu(k,iCell) ),0.69);
                }
            // }
            });
            Kokkos::fence();
        }
        
        // Smagorinsky eddy viscosity, based on horizontal deformation (in this case on model coordinate surfaces).
        // The integration coefficients were precomputed and stored in defc_a and defc_b
        if(!strcmp(config_horiz_mixing,"2d_smagorinsky")){
            Kokkos::parallel_for(
               "Cell_Kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
                real d_diag_s[nVertLevels + 1];
                real d_off_diag_s[nVertLevels + 1];
                auto d_diag = real1d((real*)d_diag_s, {1},{nVertLevels+2});
                auto d_off_diag = real1d((real*)d_off_diag_s, {1},{nVertLevels+2});
            // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
                for(int k = 1; k <= nVertLevels; k++){
                    d_diag(k) = 0.0;
                    d_off_diag(k) = 0.0;
                }
                for(int iEdge = 1; iEdge <= nEdgesOnCell(iCell); iEdge++){
                    for(int k = 1; k <= nVertLevels; k++){
                        d_diag(k)     = d_diag(k)  + defc_a(iEdge,iCell)*u(k,edgesOnCell(iEdge,iCell))  
                                                   - defc_b(iEdge,iCell)*v(k,edgesOnCell(iEdge,iCell));
                        d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,edgesOnCell(iEdge,iCell))  
                                                   + defc_a(iEdge,iCell)*v(k,edgesOnCell(iEdge,iCell));
                    }
                }
                for(int k = 1; k <= nVertLevels; k++){
                    // here is the Smagorinsky formulation, 
                    // followed by imposition of an upper bound on the eddy viscosity
                    kdiff(k,iCell) = Kokkos::min(Kokkos::pow((c_s * config_len_disp),2) * Kokkos::sqrt(Kokkos::pow(d_diag(k),2) + Kokkos::pow(d_off_diag(k),2)),(0.01*Kokkos::pow(config_len_disp,2)) * invDt);
                }
            // }
            });
            Kokkos::fence();
            h_mom_eddy_visc4   = config_visc4_2dsmag * Kokkos::pow(config_len_disp,3);
            h_theta_eddy_visc4 = h_mom_eddy_visc4;
        }else if(!strcmp(config_horiz_mixing,"2d_fixed")){
            Kokkos::parallel_for(
               "Cell_Kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
                for(int k = 1; k <= nVertLevels; k++){
                    kdiff(k,iCell) = config_h_theta_eddy_visc2;
                }
            // }
            });
            Kokkos::fence();
            h_mom_eddy_visc4 = config_h_mom_eddy_visc4;
            h_theta_eddy_visc4 = config_h_theta_eddy_visc4;
        }
        
        if(config_iamas_cam_coef > 0.0){
            Kokkos::parallel_for(
               "Cell_Kokkos", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
                // 2nd-order filter for top absorbing layer as in CAM-SE :  WCS 10 May 2017
                // From iAMAS-CAM V4.0 code, with addition to config-specified coefficient (V4.0_coef = 0.2; SE_coef = 1.0)
                kdiff(nVertLevels-2,iCell) = Kokkos::max(kdiff(nVertLevels-2,iCell),    2.0833*config_len_disp*config_iamas_cam_coef);
                kdiff(nVertLevels-1,iCell) = Kokkos::max(kdiff(nVertLevels-1,iCell),2.0*2.0833*config_len_disp*config_iamas_cam_coef);
                kdiff(nVertLevels  ,iCell) = Kokkos::max(kdiff(nVertLevels  ,iCell),4.0*2.0833*config_len_disp*config_iamas_cam_coef);
            // }
            }); 
            Kokkos::fence();
        }
        
    }
    /**
     * tendency for density.
     * accumulate total water here for later use in w tendency calculation.
     * accumulate horizontal mass-flux
    */
   Kokkos::parallel_for(
    "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
        for(int k = 1; k<= nVertLevels; k++){
            h_divergence(k,iCell) = 0.0;
        }
        for(int i = 1; i <= nEdgesOnCell(iCell); i++){
            int iEdge = edgesOnCell(i,iCell);
            real edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge);
            for(int k = 1; k <= nVertLevels; k++){
                h_divergence(k,iCell) = h_divergence(k,iCell) + edge_sign * ru(k,iEdge);
            }
        }
    // }
    });
    Kokkos::fence();

    // compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho
    Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
        real r = invAreaCell(iCell);
        for(int k = 1; k <= nVertLevels; k++){
            h_divergence(k,iCell) = h_divergence(k,iCell) * r;
        }
    // }
    });
    Kokkos::fence();

    /**
     * dp / dz and tend_rho
     * 
     * only needed on first rk_step with pert variables defined a pert from time t
    */
    if(rk_step == 1){
        real rgas_cprcv = rgas * cp / cv;
        Kokkos::parallel_for(
            "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
            for(int k = 1; k <= nVertLevels; k++){
                tend_rho(k,iCell) = -h_divergence(k,iCell)-rdzw(k)*(rw(k+1,iCell)-rw(k,iCell)) + tend_rho_physics(k,iCell);
                dpdz(k,iCell) = -gravity*(rb(k,iCell)*(qtot(k,iCell)) + rr_save(k,iCell)*(1.+qtot(k,iCell)));
            }
        // }
        });
        Kokkos::fence();
    }
    

    // Compute u (normal) velocity tendency for each edge (cell face)
    Kokkos::parallel_for(
        "Egde_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
    // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
        real wduz_s[nVertLevels + 1];
        auto wduz = real1d((real*)wduz_s, {1},{nVertLevels+2});
        real q_s[nVertLevels];
        auto q = real1d((real*)q_s, {1},{nVertLevels+1});
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        // horizontal pressure gradient 
        if(rk_step == 1){
            for(int k = 1; k <= nVertLevels; k++){
                tend_u_euler(k,iEdge) =  -cqu(k,iEdge)*( (pp(k,cell2)-pp(k,cell1))*invDcEdge(iEdge)/(.5*(zz(k,cell2)+zz(k,cell1))) 
                                            -0.5*zxu(k,iEdge)*(dpdz(k,cell1)+dpdz(k,cell2)) );
            }
        }

        // vertical transport of u
        wduz(1) = 0.;
        int k = 2;
        wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge));
        for(int k = 3; k <= nVertLevels - 1; k++){
            wduz(k) = flux3( u(k-2,iEdge),u(k-1,iEdge),u(k,iEdge),u(k+1,iEdge),0.5*(rw(k,cell1)+rw(k,cell2)), 1.0);
        }
        k = nVertLevels;
        wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge));

        wduz(nVertLevels+1) = 0.;
        for(int k = 1; k <= nVertLevels; k++){
            tend_u(k,iEdge) = - rdzw(k)*(wduz(k+1)-wduz(k)); //  first use of tend_u
        }
        // Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009
        for(int k = 1; k <= nVertLevels; k++){
            q(k) = 0.0;
        }
        for(int j = 1; j <= nEdgesOnEdge(iEdge); j++){
            int eoe = edgesOnEdge(j,iEdge);
            for(int k = 1; k <= nVertLevels; k++){
                real workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe));
                // the original definition of pv_edge had a factor of 1/density.  We have removed that factor
                // given that it was not integral to any conservation property of the system
                q(k) = q(k) + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv;
            }
        }
        for(int k = 1; k <= nVertLevels; k++){
            // horizontal ke gradient and vorticity terms in the vector invariant formulation
            // of the horizontal momentum equation
            tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge)* (q(k) - (ke(k,cell2) - ke(k,cell1)) * invDcEdge(iEdge))                    
                                             - u(k,iEdge)*0.5*(h_divergence(k,cell1)+h_divergence(k,cell2));
            #ifdef CURVATURE
            // curvature terms for the sphere
            tend_u(k,iEdge) = tend_u(k,iEdge) 
                             - 2.*omega*Kokkos::cos(angleEdge(iEdge))*Kokkos::cos(latEdge(iEdge))  
                               *rho_edge(k,iEdge)*.25*(w(k,cell1)+w(k+1,cell1)+w(k,cell2)+w(k+1,cell2))          
                             - u(k,iEdge)*.25*(w(k+1,cell1)+w(k,cell1)+w(k,cell2)+w(k+1,cell2))                  
                               *rho_edge(k,iEdge) * inv_r_earth;
            #endif
        }
    // }
    });
    Kokkos::fence();
    /**
     * horizontal mixing for u
     * mixing terms are integrated using forward-Euler, so this tendency is only computed in the
     * first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.
     */
    if(rk_step == 1){
        // del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
        // First, storage to hold the result from the first del^2 computation.
        Kokkos::parallel_for(
            "Edge_Kokkos",Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
        // for(int iEdge = edgeStart; iEdge <= edgeEnd; iEdge++){
            for(int k = 1; k <= nVertLevels; k++){
                delsq_u(k,iEdge) = 0.0;
            }
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            int vertex1 = verticesOnEdge(1,iEdge);
            int vertex2 = verticesOnEdge(2,iEdge);
            real r_dc = invDcEdge(iEdge);
            real r_dv = Kokkos::min(invDvEdge(iEdge), 4*invDcEdge(iEdge));
            for(int k = 1; k <= nVertLevels; k++){
                // Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                // only valid for h_mom_eddy_visc4 == constant
                real u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) * r_dc  
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) * r_dv;

                delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion;

                real kdiffu = 0.5*(kdiff(k,cell1)+kdiff(k,cell2));

                // jwfeng: Kinetic viscosity
                real u_turb_flux_kin_visc = 0.5*(kin_visc(k,cell1) + kin_visc(k,cell2)) * meshScalingDel2(iEdge) * u_diffusion;
                // include 2nd-orer diffusion here 
                tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) 
                                        + rho_edge(k,iEdge)* kdiffu * u_diffusion * meshScalingDel2(iEdge)
                                        + u_turb_flux_kin_visc;
            }
        // }
        });
        Kokkos::fence();
        if(h_mom_eddy_visc4 > 0.0){ // 4th order mixing is active
            Kokkos::parallel_for(
               "Vertex_Kokkos", Kokkos::RangePolicy<>(vertexStart, vertexEnd + 1), KOKKOS_LAMBDA(int iVertex) {
            // for(int iVertex = vertexStart; iVertex <= vertexEnd; iVertex++){
                for(int k = 1; k <= nVertLevels; k++){
                    delsq_vorticity(k,iVertex) = 0.0;
                }
                for(int i = 1; i <= vertexDegree; i++){
                    int iEdge = edgesOnVertex(i,iVertex);
                    real edge_sign = invAreaTriangle(iVertex) * dcEdge(iEdge) * edgesOnVertex_sign(i,iVertex);
                    for(int k = 1; k <= nVertLevels; k++){
                        delsq_vorticity(k,iVertex) = delsq_vorticity(k,iVertex) + edge_sign * delsq_u(k,iEdge);
                    }
                }
            // }
            });
            Kokkos::fence();
            Kokkos::parallel_for(
                "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
                for(int k = 1; k <= nVertLevels; k++){
                    delsq_divergence(k,iCell) = 0.0;
                }
                real r = invAreaCell(iCell);
                for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                    int iEdge = edgesOnCell(i,iCell);
                    real edge_sign = r * dvEdge(iEdge) * edgesOnCell_sign(i,iCell);
                    for(int k = 1; k <= nVertLevels; k++){
                        delsq_divergence(k,iCell) = delsq_divergence(k,iCell) + edge_sign * delsq_u(k,iEdge);
                    }
                }
            // }
            });
            Kokkos::fence();
            Kokkos::parallel_for(
                "Edge_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
            // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                int vertex1 = verticesOnEdge(1,iEdge);
                int vertex2 = verticesOnEdge(2,iEdge);

                real u_mix_scale = meshScalingDel4(iEdge)*h_mom_eddy_visc4;
                real r_dc = u_mix_scale * config_del4u_div_factor * invDcEdge(iEdge);
                real r_dv = u_mix_scale * Kokkos::min(invDvEdge(iEdge), 4*invDcEdge(iEdge));
                for(int k = 1; k <= nVertLevels; k++){
                    /**
                     * Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                     * only valid for h_mom_eddy_visc4 == constant
                     * Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor 
                     * relative to the rotational part.  The stability constraint on the divergence component is much less
                     * stringent than the rotational part, and this flexibility may be useful.
                     */
                    real u_diffusion =  rho_edge(k,iEdge) *  ( ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) * r_dc  
                                                       -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) * r_dv );
                    tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - u_diffusion;
                }
            // }
            });
            Kokkos::fence();

        } // 4th order mixing is active 
        
        // vertical mixing for u - 2nd order filter in physical (z) space
        if(v_mom_eddy_visc2 > 0.0){
            if (config_mix_full){ // mix full state
                Kokkos::parallel_for(
                    "Edge_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
                // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
                    int cell1 = cellsOnEdge(1,iEdge);
                    int cell2 = cellsOnEdge(2,iEdge);
                    for(int k = 2; k <= nVertLevels - 1; k++){
                        real z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2));
                        real z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2));
                        real z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2));
                        real z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2));

                        real zm = 0.5*(z1+z2);
                        real z0 = 0.5*(z2+z3);
                        real zp = 0.5*(z3+z4);

                        tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  
                                            (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      
                                        -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm));
                    }
                // }
                });
                Kokkos::fence();
                
            }else{ // idealized cases where we mix on the perturbation from the initial 1-D state
                Kokkos::parallel_for(
                    "Edge_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
                // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
                    real u_mix_s[nVertLevels];
                    auto u_mix = real1d((real*)u_mix_s, {1},{nVertLevels+1});
                    int cell1 = cellsOnEdge(1,iEdge);
                    int cell2 = cellsOnEdge(2,iEdge);
                    for(int k = 1; k <= nVertLevels; k++){
                        u_mix(k) = u(k,iEdge) - u_init(k) * Kokkos::cos( angleEdge(iEdge) ) 
                                           - v_init(k) * Kokkos::sin( angleEdge(iEdge) );
                    }
                    for(int k = 2; k <= nVertLevels - 1; k++){
                        real z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2));
                        real z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2));
                        real z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2));
                        real z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2));

                        real zm = 0.5*(z1+z2);
                        real z0 = 0.5*(z2+z3);
                        real zp = 0.5*(z3+z4);

                        tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  
                                            (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      
                                        -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm));
                    }
                // }
                });
                Kokkos::fence();
            } // mix perturbation state
        } // vertical mixing of horizontal momentum
    } // (rk_step 1 test for computing mixing terms)
    
    // add in mixing and physics tendency for u
    // Rayleigh damping on u
    if(config_rayleigh_damp_u){
        real rayleigh_coef_inverse = 1.0 / ( double(config_number_rayleigh_damp_u_levels) 
                                         * (config_rayleigh_damp_u_timescale_days*seconds_per_day) );
        Kokkos::parallel_for(
            "k_Kokkos",Kokkos::RangePolicy<>(nVertLevels-config_number_rayleigh_damp_u_levels+1, nVertLevels + 1), KOKKOS_LAMBDA(int k) {
        // for(int k = nVertLevels-config_number_rayleigh_damp_u_levels+1; k <= nVertLevels; k++){
            rayleigh_damp_coef(k) = real(k - (nVertLevels-config_number_rayleigh_damp_u_levels))*rayleigh_coef_inverse;
        // }
        });
        Kokkos::fence();
        Kokkos::parallel_for(
            "Edge_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
        // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
            for(int k = nVertLevels-config_number_rayleigh_damp_u_levels+1; k <= nVertLevels; k++){
                tend_u(k,iEdge) = tend_u(k,iEdge) - rho_edge(k,iEdge)*u(k,iEdge)*rayleigh_damp_coef(k);
            }
        // }
        });
        Kokkos::fence();
    }
    
    Kokkos::parallel_for(
        "Edge_Kokkos",Kokkos::RangePolicy<>(edgeSolveStart, edgeSolveEnd + 1), KOKKOS_LAMBDA(int iEdge) {
    // for(int iEdge = edgeSolveStart; iEdge <= edgeSolveEnd; iEdge++){
        for(int k = 1; k <= nVertLevels; k++){
            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge) + tend_ru_physics(k,iEdge);
        }
    // }
    });
    Kokkos::fence();

    // rhs for w
    // horizontal advection for w

    Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){ // Technically updating fewer cells than before...
        real ru_edge_w_s[nVertLevels];
        auto ru_edge_w = real1d((real*)ru_edge_w_s, {1},{nVertLevels+1});
        real flux_arr_s[nVertLevels + 1];
        auto flux_arr = real1d((real*)flux_arr_s, {1},{nVertLevels + 2});
        for(int k = 1; k <= nVertLevels + 1; k++){
            tend_w(k,iCell) = 0.0;
        }
        for(int i = 1; i <= nEdgesOnCell(iCell); i++){
            int iEdge = edgesOnCell(i,iCell);
            real edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * 0.5;
            for(int k = 2; k <= nVertLevels; k++){
                ru_edge_w(k) = fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge);
            }
            for(int k = 1; k <= nVertLevels; k++){
                flux_arr(k) = 0.0;
            }
            // flux_arr stores the value of w at the cell edge used in the horizontal transport
            for(int j = 1; j <= nAdvCellsForEdge(iEdge); j++){
                int iAdvCell = advCellsForEdge(j,iEdge);
                for(int k = 2; k <= nVertLevels; k++){
                    real scalar_weight = adv_coefs(j,iEdge) + Kokkos::copysign(1.0,ru_edge_w(k)) * adv_coefs_3rd(j,iEdge);
                    flux_arr(k) = flux_arr(k) + scalar_weight * w(k,iAdvCell);
                }
            }

            for(int k = 2; k <= nVertLevels; k++){
                tend_w(k,iCell) = tend_w(k,iCell) - edgesOnCell_sign(i,iCell) * ru_edge_w(k)*flux_arr(k);
            }
        }
    // }
    });
    Kokkos::fence();

    #ifdef CURVATURE
        Kokkos::parallel_for(
            "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
            for(int k = 2; k <= nVertLevels; k++){
                tend_w(k,iCell) = tend_w(k,iCell) + (rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))*          
                                      ( Kokkos::pow(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell),2.0)             
                                       +Kokkos::pow(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell),2.0) )/r_earth   
                                + 2.*omega*Kokkos::cos(latCell(iCell))                                       
                                       *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))                 
                                       *(rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k));
            }
        // }
        });
        Kokkos::fence();
    #endif

    // horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),

    // but here we can also code in hyperdiffusion if we wish (2nd order at present)
    if(rk_step == 1){
        /**
         * del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         * First, storage to hold the result from the first del^2 computation.
         * we copied code from the theta mixing, hence the theta* names.
         */
        Kokkos::parallel_for(
            "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
            for(int k = 1; k <= nVertLevels; k++){
                delsq_w(k,iCell) = 0.0;
            }
            for(int k = 1; k <= nVertLevels + 1; k++){
                tend_w_euler(k,iCell) = 0.0;
            }
            real r_areaCell = invAreaCell(iCell);
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iEdge = edgesOnCell(i,iCell);
                real edge_sign = 0.5 * r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge);

                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                for(int k = 2; k <= nVertLevels; k++){
                    real w_turb_flux_kin_visc = 2.0*edge_sign*(w(k,cell2)-w(k,cell1));
                    real w_turb_flux =  edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1));
                    delsq_w(k,iCell) = delsq_w(k,iCell) + w_turb_flux;
                    w_turb_flux = w_turb_flux * meshScalingDel2(iEdge) * 0.25 * 
                                    (kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2));
                    w_turb_flux_kin_visc = w_turb_flux_kin_visc * meshScalingDel2(iEdge) * 0.25 * 
                                    (kin_visc(k,cell1)+kin_visc(k,cell2)+kin_visc(k-1,cell1)+kin_visc(k-1,cell2));
                    tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + w_turb_flux + w_turb_flux_kin_visc;
                }
            }
        // }
        });
        Kokkos::fence();
        if (h_mom_eddy_visc4 > 0.0){ // 4th order mixing is active
            Kokkos::parallel_for(
                "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){ // Technically updating fewer cells than before...
                real r_areaCell = h_mom_eddy_visc4 * invAreaCell(iCell);
                for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                    int iEdge = edgesOnCell(i,iCell);
                    int cell1 = cellsOnEdge(1,iEdge);
                    int cell2 = cellsOnEdge(2,iEdge);

                    real edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell) * invDcEdge(iEdge);
                    for(int k = 2; k <= nVertLevels; k++){
                        tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - edge_sign * (delsq_w(k,cell2) - delsq_w(k,cell1));
                    }
                }
            // }
            });
            Kokkos::fence();

        } // 4th order mixing is active 
    } // horizontal mixing for w computed in first rk_step
    
    /**
     * Note for OpenMP parallelization: We could avoid allocating the delsq_w scratch
     * array, and just use the delsq_theta array as was previously done; however,
     * particularly when oversubscribing cores with threads, there is the risk that
     * some threads may reach code further below that re-uses the delsq_theta array, 
     * in which case we would need a barrier somewhere between here and that code 
     * below to ensure correct behavior. 
     */

     // vertical advection, pressure gradient and buoyancy for w
     Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
     // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
        real wdwz_s[nVertLevels + 1];
        auto wdwz = real1d((real*)wdwz_s, {1},{nVertLevels+2});
        wdwz(1) = 0.0;

        int k = 2;
        wdwz(k) =  0.25*(rw(k,iCell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell));
        for(int k = 3; k <= nVertLevels - 1; k++){
            wdwz(k) = flux3( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)), 1.0);
        }
        k = nVertLevels;
        wdwz(k) =  0.25*(rw(k,iCell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell));

        wdwz(nVertLevels+1) = 0.0;
        // next we are also dividing through by the cell area after the horizontal flux divergence
        for(int k = 2; k <= nVertLevels; k++){
            tend_w(k,iCell) = tend_w(k,iCell) * invAreaCell(iCell) -rdzu(k)*(wdwz(k+1)-wdwz(k));
        }
        if(rk_step == 1){
            for(int k = 2; k <= nVertLevels; k++){
                tend_w_pgf(k,iCell) = -cqw(k,iCell)*rdzu(k)*(pp(k,iCell)-pp(k-1,iCell));
                // dpdz is the buoyancy term here.
                tend_w_buoy(k,iCell) = cqw(k,iCell)*(fzm(k)*dpdz(k,iCell) + fzp(k)*dpdz(k-1,iCell));

                tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + tend_w_pgf(k,iCell) + tend_w_buoy(k,iCell);
            }
        }
     // }
     });
     Kokkos::fence();

     if(rk_step == 1){
        if(v_mom_eddy_visc2 > 0.0){
            Kokkos::parallel_for(
                "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
                for(int k = 2; k <= nVertLevels; k++){
                    tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  
                        (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              
                       -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k);
                }
            // }
            });
            Kokkos::fence();
        }
     } // mixing term computed first rk_step
     
     // add in mixing terms for w
     Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
     // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
        for(int k = 2; k <= nVertLevels; k++){
            tend_w(k,iCell) = tend_w(k,iCell) + tend_w_euler(k,iCell);
        }
     // }
     });
     Kokkos::fence();

     // rhs for theta
     // horizontal advection for theta
     Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
     // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){ // Technically updating fewer cells than before...
        real flux_arr_s[nVertLevels + 1];
        auto flux_arr = real1d((real*)flux_arr_s, {1},{nVertLevels + 2});
        for(int k = 1; k <= nVertLevels; k++){
            tend_theta(k,iCell) = 0.0;
        }
        for(int i = 1; i <= nEdgesOnCell(iCell); i++){
            int iEdge = edgesOnCell(i,iCell);
            for(int k = 1; k <= nVertLevels; k++){
                flux_arr(k) = 0.0;
            }
            for(int j = 1; j <= nAdvCellsForEdge(iEdge); j++){
                int iAdvCell = advCellsForEdge(j,iEdge);
                for(int k = 1; k <= nVertLevels; k++){
                    real scalar_weight = adv_coefs(j,iEdge) + Kokkos::copysign(1.0,ru(k,iEdge))*adv_coefs_3rd(j,iEdge);
                    flux_arr(k) = flux_arr(k) + scalar_weight* theta_m(k,iAdvCell);
                }
            }
            for(int k = 1; k <= nVertLevels; k++){
                tend_theta(k,iCell) = tend_theta(k,iCell) - edgesOnCell_sign(i,iCell)*ru(k,iEdge)*flux_arr(k);
            }
        }
     // }
     });
     Kokkos::fence();

     // addition to pick up perturbation flux for rtheta_pp equation
     if(rk_step > 1){
        Kokkos::parallel_for(
            "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iEdge = edgesOnCell(i,iCell);
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                for(int k = 1; k <= nVertLevels; k++){
                    real flux = edgesOnCell_sign(i,iCell)*dvEdge(iEdge)*(ru_save(k,iEdge)-ru(k,iEdge))*0.5*(theta_m_save(k,cell2)+theta_m_save(k,cell1));
                    tend_theta(k,iCell) = tend_theta(k,iCell)-flux;  // division by areaCell picked up down below
                }
            }
        // }
        });
        Kokkos::fence();
     }
     // horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
     // but here we can also code in hyperdiffusion if we wish (2nd order at present)
     if(rk_step == 1){
        Kokkos::parallel_for(
            "Cell_Kokkos",Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; iCell++){
            for(int k = 1; k <= nVertLevels; k++){
                delsq_theta(k,iCell) = 0.0;
            }
            for(int k = 1; k <= nVertLevels; k++){
                tend_theta_euler(k,iCell) = 0.0;
            }
            real r_areaCell = invAreaCell(iCell);
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iEdge = edgesOnCell(i,iCell);
                real edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge);
                real pr_scale = prandtl_inv * meshScalingDel2(iEdge);
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                for(int k = 1; k <= nVertLevels; k++){
                    // we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below
                    real theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge);
                    delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux;
                    theta_turb_flux = theta_turb_flux*0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * pr_scale;
                    tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux;
                }
            }
        // }
        });
        Kokkos::fence();
        if(h_theta_eddy_visc4 > 0.0){ // 4th order mixing is active
            Kokkos::parallel_for(
                "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
            // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){ // Technically updating fewer cells than before...
                real r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell);
                for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                    int iEdge = edgesOnCell(i,iCell);
                    real edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge);

                    int cell1 = cellsOnEdge(1,iEdge);
                    int cell2 = cellsOnEdge(2,iEdge);
                    for(int k = 1; k <= nVertLevels; k++){
                        tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1));
                    }
                }
            // }
            });
            Kokkos::fence();

        } // 4th order mixing is active 
     } // theta mixing calculated first rk_step
     
     /**
      * vertical advection plus diabatic term
      * Note: we are also dividing through by the cell area after the horizontal flux divergence
      */

    Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
        real wdtz_s[nVertLevels + 1];
        auto wdtz = real1d((real*)wdtz_s, {1},{nVertLevels+2});
        wdtz(1) = 0.0;

        int k = 2;
        wdtz(k) =  rw(k,iCell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell));
        wdtz(k) =  wdtz(k) + (rw_save(k,iCell)-rw(k,iCell))*(fzm(k)*theta_m_save(k,iCell)+fzp(k)*theta_m_save(k-1,iCell));
        for(int k = 3; k <= nVertLevels - 1; k++){
            wdtz(k) = flux3( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell), coef_3rd_order );
            wdtz(k) =  wdtz(k) + (rw_save(k,iCell)-rw(k,iCell))*(fzm(k)*theta_m_save(k,iCell)+fzp(k)*theta_m_save(k-1,iCell));  // rtheta_pp redefinition
        }
        k = nVertLevels;
        wdtz(k) =  rw_save(k,iCell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell));  // rtheta_pp redefinition

        wdtz(nVertLevels+1) = 0.0;
        for(int k = 1; k <= nVertLevels; k++){
            tend_theta(k,iCell) = tend_theta(k,iCell)*invAreaCell(iCell) -rdzw(k)*(wdtz(k+1)-wdtz(k));
            tend_rtheta_adv(k,iCell) = tend_theta(k,iCell);   
            rthdynten(k,iCell) = tend_theta(k,iCell)/rho_zz(k,iCell);  
            tend_theta(k,iCell) = tend_theta(k,iCell) + rho_zz(k,iCell)*rt_diabatic_tend(k,iCell);
        }
    // }
    });
    Kokkos::fence();

    // vertical mixing for theta - 2nd order
    if(rk_step == 1){
        if ( v_theta_eddy_visc2 > 0.0 ){ // vertical mixing for theta_m
            if (config_mix_full){
                Kokkos::parallel_for(
                    "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
                // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
                    for(int k = 2; k <= nVertLevels - 1; k++){
                        real z1 = zgrid(k-1,iCell);
                        real z2 = zgrid(k  ,iCell);
                        real z3 = zgrid(k+1,iCell);
                        real z4 = zgrid(k+2,iCell);

                        real zm = 0.5*(z1+z2);
                        real z0 = 0.5*(z2+z3);
                        real zp = 0.5*(z3+z4);

                        tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(
                                                (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 
                                                -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm));
                    }
                // }
                });
                Kokkos::fence();
            }else{ // idealized cases where we mix on the perturbation from the initial 1-D state
                Kokkos::parallel_for(
                    "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
                // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
                    for(int k = 2; k <= nVertLevels - 1; k++){
                        real z1 = zgrid(k-1,iCell);
                        real z2 = zgrid(k  ,iCell);
                        real z3 = zgrid(k+1,iCell);
                        real z4 = zgrid(k+2,iCell);

                        real zm = 0.5*(z1+z2);
                        real z0 = 0.5*(z2+z3);
                        real zp = 0.5*(z3+z4);

                        tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(
                                                ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      
                                                -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm));
                    }
                // }
                });
                Kokkos::fence();
            }
        }
    } // compute vertical theta mixing on first rk_step
    
    Kokkos::parallel_for(
        "Cell_Kokkos",Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){
        for(int k = 1; k <= nVertLevels; k++){
            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell) + tend_rtheta_physics(k,iCell);
        }
    // }
    });
    Kokkos::fence();
}
