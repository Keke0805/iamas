! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_todynamics
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_dmpar

 use iamas_atmphys_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_get_tend


!Interface between the physics parameterizations and the non-hydrostatic dynamical core.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
!
! subroutines in iamas_atmphys_todynamics:
! ---------------------------------------
! physics_get_tend: add and mass-weigh tendencies before being added to dynamics tendencies.
! tend_toEdges   : interpolate wind-tendencies from centers to edges of grid-cells.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added calculation of the advective tendency of the potential temperature due to horizontal
!   and vertical advection, and horizontal mixing (diffusion).
!   Laura D. Fowler (birch.mmm.ucar.edu) / 2013-11-19.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * modified the sourcecode to accomodate the packages "cu_kain_fritsch_in" and "cu_tiedtke_in".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-24.
! * added the calculation of rthdynten which is the tendency of potential temperature due to horizontal and
!   vertical advections needed in the Grell-Freitas scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * added the option bl_mynn for the calculation of the tendency for the cloud ice number concentration.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * in subroutine physics_get_tend_work, added the option cu_ntiedtke in the calculation of rucuten_Edge.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-28.

interface physics_get_tend_interface1
   subroutine physics_get_tend_phase1_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve,&
         tend_u,rublten_Edge,mass_edge,tend_th,tend_scalars,rqvblten,rqcblten,rqiblten,&
         index_qv,index_qc,index_qi,mass,rthblten) bind(C, name="physics_get_tend_phase1_kokkos")
         use, intrinsic :: iso_c_binding
         use :: flcl_ndarray_mod
         integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve
         integer(c_int), intent(in) :: index_qv, index_qc, index_qi
         type(nd_array_t) :: mass_edge
         type(nd_array_t) :: rublten_Edge
         type(nd_array_t) :: tend_u
         type(nd_array_t) :: tend_th
         type(nd_array_t) :: tend_scalars
         type(nd_array_t) :: mass
         type(nd_array_t) :: rthblten
         type(nd_array_t) :: rqvblten
         type(nd_array_t) :: rqcblten
         type(nd_array_t) :: rqiblten
   end subroutine physics_get_tend_phase1_kokkos
   subroutine physics_get_tend_phase2_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
      tend_th,tend_scalars,&
      index_qv,index_qc,index_qi,mass,rthcuten,rqvcuten,rqccuten,rqicuten) bind(C, name="physics_get_tend_phase2_kokkos")
      use, intrinsic :: iso_c_binding
      use :: flcl_ndarray_mod
      integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
      integer(c_int), intent(in) :: index_qv, index_qc, index_qi
      type(nd_array_t) :: tend_th
      type(nd_array_t) :: tend_scalars
      type(nd_array_t) :: mass
      type(nd_array_t) :: rthcuten
      type(nd_array_t) :: rqvcuten
      type(nd_array_t) :: rqccuten
      type(nd_array_t) :: rqicuten
   end subroutine physics_get_tend_phase2_kokkos
   subroutine physics_get_tend_phase3_kokkos(nVertLevels, nCells, nEdges, num_scalars, nCellsSolve,&
      tend_scalars,tend_th,&
      index_qv,mass,rthfddaten,rqvfddaten) bind(C, name="physics_get_tend_phase3_kokkos")
      use, intrinsic :: iso_c_binding
      use :: flcl_ndarray_mod
      integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
      integer(c_int), intent(in) :: index_qv
      type(nd_array_t) :: tend_th
      type(nd_array_t) :: tend_scalars
      type(nd_array_t) :: mass
      type(nd_array_t) :: rthfddaten
      type(nd_array_t) :: rqvfddaten
   end subroutine physics_get_tend_phase3_kokkos
end interface physics_get_tend_interface1
interface physics_get_tend_interface2
   subroutine physics_get_tend_phase4_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
      tend_th,tend_scalars,&
      index_qv,theta_m,tend_theta,scalars) bind(C, name="physics_get_tend_phase4_kokkos")
      use, intrinsic :: iso_c_binding
      use :: flcl_ndarray_mod
      integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
      integer(c_int), intent(in) :: index_qv
      type(nd_array_t) :: tend_th
      type(nd_array_t) :: tend_scalars
      type(nd_array_t) :: theta_m
      type(nd_array_t) :: tend_theta
      type(nd_array_t) :: scalars
   end subroutine physics_get_tend_phase4_kokkos
   subroutine tend_toEdges_kokkos(nCells,nVertLevels,nCellsSolve,nEdges,Ux_tend,Uy_tend,U_tend,edgeNormalVectors,east,north,cellsOnEdge) bind(C, name="tend_toEdges_kokkos")
      use, intrinsic :: iso_c_binding
      use :: flcl_ndarray_mod
      integer(c_int), intent(in) :: nCells, nCellsSolve, nEdges
      type(nd_array_t) :: Ux_tend
      type(nd_array_t) :: Uy_tend
      type(nd_array_t) :: U_tend
      type(nd_array_t) :: edgeNormalVectors
      type(nd_array_t) :: east
      type(nd_array_t) :: north
      type(nd_array_t) :: cellsOnEdge
   end subroutine tend_toEdges_kokkos
end interface physics_get_tend_interface2
 contains
 subroutine physics_get_tend_phase1_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve,&
   tend_u,rublten_Edge,mass_edge,tend_th,tend_scalars,rqvblten,rqcblten,rqiblten,&
   index_qv,index_qc,index_qi,mass,rthblten)
   use, intrinsic :: iso_c_binding
   use :: flcl_ndarray_mod
   integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve
   integer(c_int), intent(in) :: index_qv, index_qc, index_qi
   real(c_double), dimension(nVertLevels,nEdges+1), intent(in) :: mass_edge
   real(c_double), dimension(nVertLevels,nEdges+1), intent(inout) :: rublten_Edge
   real(c_double), dimension(nVertLevels,nEdges+1), intent(inout) :: tend_u
   real(c_double), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
   real(c_double), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: mass
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rthblten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqvblten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqcblten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqiblten
   call physics_get_tend_phase1_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve,&
   to_nd_array(tend_u),to_nd_array(rublten_Edge),to_nd_array(mass_edge),to_nd_array(tend_th),to_nd_array(tend_scalars),to_nd_array(rqvblten),to_nd_array(rqcblten),to_nd_array(rqiblten),&
   index_qv,index_qc,index_qi,to_nd_array(mass),to_nd_array(rthblten))
end subroutine physics_get_tend_phase1_kokkos_f
subroutine physics_get_tend_phase2_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
   tend_th,tend_scalars,&
   index_qv,index_qc,index_qi,mass,rthcuten,rqvcuten,rqccuten,rqicuten)
   use, intrinsic :: iso_c_binding
   use :: flcl_ndarray_mod
   integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
   integer(c_int), intent(in) :: index_qv, index_qc, index_qi
   real(c_double), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
   real(c_double), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: mass
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rthcuten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqvcuten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqccuten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqicuten
   call physics_get_tend_phase2_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
   to_nd_array(tend_th),to_nd_array(tend_scalars),&
   index_qv,index_qc,index_qi,to_nd_array(mass),to_nd_array(rthcuten),to_nd_array(rqvcuten),to_nd_array(rqccuten),to_nd_array(rqicuten))
end subroutine physics_get_tend_phase2_kokkos_f
subroutine physics_get_tend_phase3_kokkos_f(nVertLevels, nCells, nEdges, num_scalars, nCellsSolve,&
   tend_scalars,tend_th,&
   index_qv,mass,rthfddaten,rqvfddaten)
   use, intrinsic :: iso_c_binding
   use :: flcl_ndarray_mod
   integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
   integer(c_int), intent(in) :: index_qv
   real(c_double), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
   real(c_double), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: mass
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rthfddaten
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: rqvfddaten
   call physics_get_tend_phase3_kokkos(nVertLevels, nCells, nEdges,  num_scalars, nCellsSolve,&
   to_nd_array(tend_scalars),to_nd_array(tend_th),&
   index_qv,to_nd_array(mass),to_nd_array(rthfddaten),to_nd_array(rqvfddaten))
end subroutine physics_get_tend_phase3_kokkos_f
subroutine physics_get_tend_phase4_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
   tend_th,tend_scalars,&
   index_qv,theta_m,tend_theta,scalars)
   use, intrinsic :: iso_c_binding
   use :: flcl_ndarray_mod
   integer(c_int), intent(in) :: nCells, nVertLevels, nEdges, num_scalars, nCellsSolve
   integer(c_int), intent(in) :: index_qv
   real(c_double), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
   real(c_double), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
   real(c_double), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta
   real(c_double), dimension(num_scalars,nVertLevels,nCells+1), intent(in) :: scalars
   call physics_get_tend_phase4_kokkos(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
   to_nd_array(tend_th),to_nd_array(tend_scalars),&
   index_qv,to_nd_array(theta_m),to_nd_array(tend_theta),to_nd_array(scalars))
end subroutine physics_get_tend_phase4_kokkos_f
subroutine tend_toEdges_kokkos_f(nCells,nVertLevels,nCellsSolve,nEdges,Ux_tend,Uy_tend,U_tend,edgeNormalVectors,east,north,cellsOnEdge)
   use, intrinsic :: iso_c_binding
   use :: flcl_ndarray_mod
   integer(c_int), intent(in) :: nCells, nCellsSolve, nEdges,nVertLevels
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: Ux_tend
   real(c_double), dimension(nVertLevels,nCells+1), intent(in) :: Uy_tend
   real(c_double), dimension(nVertLevels,nEdges+1), intent(inout) :: U_tend
   real(c_double), dimension(3,nEdges+1), intent(in) :: edgeNormalVectors
   real(c_double), dimension(3,nCells+1), intent(in) :: east
   real(c_double), dimension(3,nCells+1), intent(in) :: north
   integer(c_int), dimension(2,nEdges+1), intent(in) :: cellsOnEdge
   call tend_toEdges_kokkos(nCells,nVertLevels,nCellsSolve,nEdges,to_nd_array(Ux_tend),to_nd_array(Uy_tend),to_nd_array(U_tend),to_nd_array(edgeNormalVectors),to_nd_array(east),to_nd_array(north),to_nd_array(cellsOnEdge))
end subroutine tend_toEdges_kokkos_f
!=================================================================================================================
 subroutine physics_get_tend( block, mesh, state, diag, tend, tend_physics, configs, rk_step, dynamics_substep, &
                              tend_ru_physics, tend_rtheta_physics, tend_rho_physics )
!=================================================================================================================
   
 use iamas_atm_dimensions

!input variables:
 type(block_type),intent(in),target:: block
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: state
 type(iamas_pool_type),intent(in):: configs
 integer, intent(in):: rk_step
 integer, intent(in):: dynamics_substep

!inout variables:
 type(iamas_pool_type),intent(inout):: diag
 type(iamas_pool_type),intent(inout):: tend
 type(iamas_pool_type),intent(inout):: tend_physics

 real(kind=RKIND),intent(inout),dimension(:,:) :: tend_ru_physics, tend_rtheta_physics, tend_rho_physics

!local variables:
 character(len=StrKIND),pointer :: config_pbl_scheme, &
                                   config_convection_scheme, &
                                   config_radt_lw_scheme, &
                                   config_radt_sw_scheme, &
                                   config_fdda_scheme

 integer:: i,iCell,k,n
 integer,pointer:: index_qv, index_qc, index_qr, index_qi, index_qs, index_qg
 integer,pointer:: index_ni
 integer,pointer:: nCells,nCellsSolve,nEdges,nEdgesSolve

 real(kind=RKIND),dimension(:,:),pointer:: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: mass_edge     ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: theta_m       ! time level 1
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rniblten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten, &
                                           rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rthfddaten,rqvfddaten,rufddaten,rvfddaten                                 
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 real(kind=RKIND),dimension(:,:),pointer:: rthdynten
 
 real(kind=RKIND),dimension(:,:),pointer:: tend_rtheta_adv
 real(kind=RKIND),dimension(:,:),pointer:: tend_u_phys !nick
 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_theta_euler,tend_diabatic,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars
 real(kind=RKIND):: coeff

 real(kind=RKIND):: tem
 real(kind=RKIND),dimension(:,:),pointer:: rublten_Edge,rucuten_Edge,rufddaten_Edge

 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th


!=================================================================================================================
 call iamas_pool_get_dimension(mesh, 'nCells', nCells)
 call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
 call iamas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

 call iamas_pool_get_config(configs, 'config_pbl_scheme', config_pbl_scheme)
 call iamas_pool_get_config(configs, 'config_convection_scheme', config_convection_scheme)
 call iamas_pool_get_config(configs, 'config_radt_lw_scheme', config_radt_lw_scheme)
 call iamas_pool_get_config(configs, 'config_radt_sw_scheme', config_radt_sw_scheme)
 call iamas_pool_get_config(configs, 'config_fdda_scheme'   , config_fdda_scheme)

 call iamas_pool_get_array(state, 'theta_m', theta_m, 1)
 call iamas_pool_get_array(state, 'scalars', scalars, 1)
 call iamas_pool_get_array(state, 'rho_zz', mass, 2)
 call iamas_pool_get_array(diag , 'rho_edge', mass_edge)

 call iamas_pool_get_array(diag , 'tend_rtheta_adv', tend_rtheta_adv)

 call iamas_pool_get_array(diag , 'tend_u_phys', tend_u_phys) !nick

 call iamas_pool_get_dimension(state, 'index_qv', index_qv)
 call iamas_pool_get_dimension(state, 'index_qc', index_qc)
 call iamas_pool_get_dimension(state, 'index_qr', index_qr)
 call iamas_pool_get_dimension(state, 'index_qi', index_qi)
 call iamas_pool_get_dimension(state, 'index_qs', index_qs)
 call iamas_pool_get_dimension(state, 'index_qg', index_qg)
 call iamas_pool_get_dimension(state, 'index_ni', index_ni)

 call iamas_pool_get_array(tend_physics, 'rublten', rublten)
 call iamas_pool_get_array(tend_physics, 'rvblten', rvblten)
 call iamas_pool_get_array(tend_physics, 'rublten_Edge', rublten_Edge)
 call iamas_pool_get_array(tend_physics, 'rthblten', rthblten)
 call iamas_pool_get_array(tend_physics, 'rqvblten', rqvblten)
 call iamas_pool_get_array(tend_physics, 'rqcblten', rqcblten)
 call iamas_pool_get_array(tend_physics, 'rqiblten', rqiblten)
 call iamas_pool_get_array(tend_physics, 'rniblten', rniblten)

 call iamas_pool_get_array(tend_physics, 'rucuten', rucuten)
 call iamas_pool_get_array(tend_physics, 'rvcuten', rvcuten)
 call iamas_pool_get_array(tend_physics, 'rucuten_Edge', rucuten_Edge)
 call iamas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
 call iamas_pool_get_array(tend_physics, 'rqvcuten', rqvcuten)
 call iamas_pool_get_array(tend_physics, 'rqccuten', rqccuten)
 call iamas_pool_get_array(tend_physics, 'rqrcuten', rqrcuten)
 call iamas_pool_get_array(tend_physics, 'rqicuten', rqicuten)
 call iamas_pool_get_array(tend_physics, 'rqscuten', rqscuten)
 call iamas_pool_get_array(tend_physics, 'rthdynten', rthdynten)
 call iamas_pool_get_array(tend_physics, 'rufddaten', rufddaten)
 call iamas_pool_get_array(tend_physics, 'rvfddaten', rvfddaten)
 call iamas_pool_get_array(tend_physics, 'rthfddaten', rthfddaten)
 call iamas_pool_get_array(tend_physics, 'rqvfddaten', rqvfddaten)

 call iamas_pool_get_array(tend,'rt_diabatic_tend',tend_diabatic)

 call iamas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
 call iamas_pool_get_array(tend_physics, 'rthratensw', rthratensw)

 call iamas_pool_get_array(tend,'u'           , tend_u         )
 call iamas_pool_get_array(tend,'theta_m'     , tend_theta     )
 call iamas_pool_get_array(tend,'theta_euler' ,tend_theta_euler)
 call iamas_pool_get_array(tend,'scalars_tend',tend_scalars    )

!initialize the tendency for the potential temperature and all scalars due to PBL, convection,
!and longwave and shortwave radiation:
! allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th = 0._RKIND

 tend_scalars(:,:,:) = 0._RKIND

 tend_ru_physics(:,:) = 0._RKIND
 tend_rtheta_physics(:,:) = 0._RKIND
 tend_rho_physics(:,:) = 0._RKIND       ! NB: rho tendency is not currently supplied by physics, but this
                                        !     field may be later filled with IAU or other tendencies

 !
 ! In case some variables are not allocated due to their associated packages,
 ! we need to make their pointers associated here to avoid triggering run-time
 ! checks when calling physics_get_tend_work
 !
 if (.not. associated(rucuten)) allocate(rucuten(0,0))
 if (.not. associated(rvcuten)) allocate(rvcuten(0,0))
 if (.not. associated(rqrcuten)) allocate(rqrcuten(0,0))
 if (.not. associated(rqscuten)) allocate(rqscuten(0,0))
 if (.not. associated(rniblten)) allocate(rniblten(0,0))
 if (.not. associated(rthdynten)) allocate(rthdynten(0,0))
 if (.not. associated(rublten)) allocate(rublten(0,0))
 if (.not. associated(rvblten)) allocate(rvblten(0,0))
 if (.not. associated(rthblten)) allocate(rthblten(0,0))
 if (.not. associated(rqvblten)) allocate(rqvblten(0,0))
 if (.not. associated(rqcblten)) allocate(rqcblten(0,0))
 if (.not. associated(rqiblten)) allocate(rqiblten(0,0))
 if (.not. associated(rthcuten)) allocate(rthcuten(0,0))
 if (.not. associated(rqvcuten)) allocate(rqvcuten(0,0))
 if (.not. associated(rqccuten)) allocate(rqccuten(0,0))
 if (.not. associated(rqicuten)) allocate(rqicuten(0,0))

    !czhao +++++++++++++++++
    ! Zero-out the theta tendency in the top six layers of the model:
    do k = nVertLevels-5, nVertLevels
       rthblten(k,:) = 0.0
    end do
    !czhao ----------------- 
 !! FDDA nudging variables
 if(.not. associated(rufddaten)) allocate(rufddaten(0,0))
 if(.not. associated(rvfddaten)) allocate(rvfddaten(0,0))
 if(.not. associated(rthfddaten)) allocate(rthfddaten(0,0))
 if(.not. associated(rqvfddaten)) allocate(rqvfddaten(0,0))


 call physics_get_tend_work(block, mesh, nCells, nEdges, nCellsSolve, nEdgesSolve, &
                           rk_step, dynamics_substep, &
                           config_pbl_scheme, config_convection_scheme, config_radt_lw_scheme, &
                           config_radt_sw_scheme, config_fdda_scheme, &
                           index_qv, index_qc, index_qr, index_qi, index_qs, index_ni, &
                           rublten, rvblten, mass_edge, rublten_Edge, &
                           tend_ru_physics, &
                           rucuten, rvcuten, rucuten_Edge, &
                           tend_th, tend_scalars, mass, rthblten, rqvblten, rqcblten, rqiblten, rniblten, &
                           rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten, &
                           rthratenlw, rthratensw, rthdynten, &
                           rufddaten, rvfddaten, rthfddaten, rqvfddaten, rufddaten_Edge, &
                           tend_u_phys, tend_rtheta_adv, tend_diabatic, &
                           theta_m, scalars, &
                           tend_rtheta_physics, &
                           tend_theta_euler &
                           )

 !
 ! Clean up any pointers that were allocated with zero size before the call to
 ! physics_get_tend_work
 !
 if (size(rucuten) == 0) deallocate(rucuten)
 if (size(rvcuten) == 0) deallocate(rvcuten)
 if (size(rqrcuten) == 0) deallocate(rqrcuten)
 if (size(rqscuten) == 0) deallocate(rqscuten)
 if (size(rniblten) == 0) deallocate(rniblten)
 if (size(rthdynten) == 0) deallocate(rthdynten)
 if (size(rublten) == 0) deallocate(rublten)
 if (size(rvblten) == 0) deallocate(rvblten)
 if (size(rthblten) == 0) deallocate(rthblten)
 if (size(rqvblten) == 0) deallocate(rqvblten)
 if (size(rqcblten) == 0) deallocate(rqcblten)
 if (size(rqiblten) == 0) deallocate(rqiblten)
 if (size(rthcuten) == 0) deallocate(rthcuten)
 if (size(rqvcuten) == 0) deallocate(rqvcuten)
 if (size(rqccuten) == 0) deallocate(rqccuten)
 if (size(rqicuten) == 0) deallocate(rqicuten)
 if (size(rufddaten) == 0) deallocate(rufddaten)
 if (size(rvfddaten) == 0) deallocate(rvfddaten)
 if (size(rthfddaten)== 0) deallocate(rthfddaten)
 if (size(rqvfddaten)== 0) deallocate(rqvfddaten)

! deallocate(theta)
 deallocate(tend_th)

! if(rk_step .eq. 3) then
!    call iamas_log_write('')
!    call iamas_log_write('--- enter subroutine physics_get_tend:')
!    call iamas_log_write('max rthblten   = $r',realArgs=(/maxval(rthblten(:,1:nCellsSolve))/))
!    call iamas_log_write('min rthblten   = $r',realArgs=(/minval(rthblten(:,1:nCellsSolve))/))
!    call iamas_log_write('max rthcuten   = $r',realArgs=(/maxval(rthcuten(:,1:nCellsSolve))/))
!    call iamas_log_write('min rthcuten   = $r',realArgs=(/minval(rthcuten(:,1:nCellsSolve))/))
!    call iamas_log_write('max rthratenlw = $r',realArgs=(/maxval(rthratenlw(:,1:nCellsSolve))/))
!    call iamas_log_write('min rthratenlw = $r',realArgs=(/minval(rthratenlw(:,1:nCellsSolve))/))
!    call iamas_log_write('max rthratensw = $r',realArgs=(/maxval(rthratensw(:,1:nCellsSolve))/))
!    call iamas_log_write('min rthratensw = $r',realArgs=(/minval(rthratensw(:,1:nCellsSolve))/))
!    call iamas_log_write('--- end subroutine physics_get_tend')
!    call iamas_log_write('')
! endif

 end subroutine physics_get_tend

 !==================================================================================================
 subroutine physics_get_tend_work(block, mesh, nCells, nEdges, nCellsSolve, nEdgesSolve, &
                                 rk_step, dynamics_substep, &
                                 config_pbl_scheme, config_convection_scheme, config_radt_lw_scheme, &
                                 config_radt_sw_scheme, config_fdda_scheme, &
                                 index_qv, index_qc, index_qr, index_qi, index_qs, index_ni, &
                                 rublten, rvblten, mass_edge, rublten_Edge, tend_u, &
                                 rucuten, rvcuten, rucuten_Edge, &
                                 tend_th, tend_scalars, mass, rthblten, rqvblten, rqcblten, rqiblten, rniblten, &
                                 rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten, &
                                 rthratenlw, rthratensw, rthdynten, &
                                 rufddaten, rvfddaten, rthfddaten, rqvfddaten, rufddaten_Edge, &
                                 tend_u_phys, tend_rtheta_adv, tend_diabatic, &
                                 theta_m, scalars, tend_theta, tend_theta_euler &
                                )
!==================================================================================================

    use iamas_atm_dimensions
 
    implicit none

    type(block_type), intent(in) :: block
    type(iamas_pool_type), intent(in) :: mesh
    integer, intent(in) :: nCells, nEdges, nCellsSolve, nEdgesSolve
    integer, intent(in) :: rk_step, dynamics_substep
    character(len=StrKIND), intent(in) :: config_pbl_scheme
    character(len=StrKIND), intent(in) :: config_convection_scheme
    character(len=StrKIND), intent(in) :: config_radt_lw_scheme
    character(len=StrKIND), intent(in) :: config_radt_sw_scheme
    character(len=StrKIND), intent(in) :: config_fdda_scheme !! fdda_scheme
    integer, intent(in) :: index_qv, index_qc, index_qr, index_qi, index_qs, index_ni
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rublten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rvblten
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: mass_edge
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rublten_Edge
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: tend_u
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rucuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rvcuten
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rucuten_Edge
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
    real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: mass
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqvblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqcblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqiblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rniblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqvcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqccuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqrcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqicuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqscuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthratenlw
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthratensw
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: rthdynten
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: tend_u_phys
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: tend_rtheta_adv
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: tend_diabatic
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
    real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(in) :: scalars
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: tend_theta_euler
    !! added for FDDA
    real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: rufddaten
    real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: rvfddaten
    real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: rthfddaten
    real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: rqvfddaten
   
    real(kind=RKIND),dimension(:,:),pointer:: rufddaten_Edge !! commented out on 2024/09/16
    !!real(kind=RKIND),intent(inout),dimension(nVertLevels,nEdges):: rufddaten_Edge


    integer :: i, k
    real (kind=RKIND) :: coeff

    !add coupled tendencies due to PBL processes:
    if (config_pbl_scheme .ne. 'off') then
       if (rk_step == 1 .and. dynamics_substep == 1) then
          call tend_toEdges(block,mesh,rublten,rvblten,rublten_Edge)

          !MGD for PV budget? should a similar line be in the cumulus section below?
          tend_u_phys(1:nVertLevels,1:nEdges) = rublten_Edge(1:nVertLevels,1:nEdges)
       end if

!   !czhao +++++++++++++++++
!   ! Zero-out the theta tendency in the top six layers of the model:
!   do k = nVertLevels-5, nVertLevels
!      rthblten(k,:) = 0.0
!      rthblten_Edge(k,:) = 0.0
!   end do
!   !czhao ----------------- 

       call physics_get_tend_phase1_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve, nEdgesSolve,&
              tend_u,rublten_Edge,mass_edge,tend_th,tend_scalars,rqvblten,rqcblten,rqiblten,&
           index_qv,index_qc,index_qi,mass,rthblten)
      !  do i = 1, nEdgesSolve
      !  do k  = 1, nVertLevels
      !     tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i)
      !  enddo
      !  enddo

      !  do i = 1, nCellsSolve
      !  do k = 1, nVertLevels
      !     tend_th(k,i) = tend_th(k,i) + rthblten(k,i)*mass(k,i)
      !     tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvblten(k,i)*mass(k,i)
      !     tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqcblten(k,i)*mass(k,i)
      !     tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqiblten(k,i)*mass(k,i)
      !  enddo
      !  enddo

       pbl_select: select case (trim(config_pbl_scheme))

          case("bl_mynn")

             do i = 1, nCellsSolve
             do k = 1, nVertLevels
                tend_scalars(index_ni,k,i) = tend_scalars(index_ni,k,i) + rniblten(k,i)*mass(k,i)
             enddo
             enddo
   
          case default         

        end select pbl_select
    endif

    !add coupled tendencies due to convection:
    if (config_convection_scheme .ne. 'off') then
      call physics_get_tend_phase2_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
           tend_th,tend_scalars,&
           index_qv,index_qc,index_qi,mass,rthcuten,rqvcuten,rqccuten,rqicuten)
      !  do i = 1, nCellsSolve
      !  do k = 1, nVertLevels
      !     tend_th(k,i) = tend_th(k,i) + rthcuten(k,i)*mass(k,i)
      !     tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvcuten(k,i)*mass(k,i)
      !     tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqccuten(k,i)*mass(k,i)
      !     tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqicuten(k,i)*mass(k,i)
      !  enddo
      !  enddo

        convection_select: select case(config_convection_scheme)

           case('cu_kain_fritsch')
              do i = 1, nCellsSolve
              do k = 1, nVertLevels
                 tend_scalars(index_qr,k,i) = tend_scalars(index_qr,k,i) + rqrcuten(k,i)*mass(k,i)
                 tend_scalars(index_qs,k,i) = tend_scalars(index_qs,k,i) + rqscuten(k,i)*mass(k,i)
              enddo
              enddo
    
           case('cu_tiedtke','cu_ntiedtke')
              if (rk_step == 1 .and. dynamics_substep == 1) then
                 call tend_toEdges(block,mesh,rucuten,rvcuten,rucuten_Edge)
                 
                 tend_u_phys(1:nVertLevels,1:nEdges) = tend_u_phys(1:nVertLevels,1:nEdges) &
                                                     + rucuten_Edge(1:nVertLevels,1:nEdges)
              end if
              do i = 1, nEdgesSolve
              do k  = 1, nVertLevels
                 tend_u(k,i)=tend_u(k,i)+rucuten_Edge(k,i)*mass_edge(k,i)
              enddo
              enddo

            case default
        end select convection_select
    endif

  !add coupled tendencies due to four-dimensional data assimilation (fdda):
 if(config_fdda_scheme .ne. 'off') then

   ! call iamas_log_write('#### Begin running fdda in todynamics', messageType=iAMAS_LOG_WARN)

   allocate(rufddaten_Edge(nVertLevels,nEdges+1))
   rufddaten_Edge(:,:) = 0.
   call tend_toEdges(block,mesh,rufddaten,rvfddaten,rufddaten_Edge)
   do i = 1, nEdgesSolve
   do k  = 1, nVertLevels
      tend_u(k,i)=tend_u(k,i)+rufddaten_Edge(k,i)*mass_edge(k,i)
   enddo
   enddo
   deallocate(rufddaten_Edge)
   call physics_get_tend_phase3_kokkos_f(nVertLevels,nCells, nEdges, num_scalars, nCellsSolve,&
           tend_scalars,tend_th,&
           index_qv,mass,rthfddaten,rqvfddaten)
  !  do i = 1, nCellsSolve
  !  do k = 1, nVertLevels
  !     tend_th(k,i) = tend_th(k,i) + rthfddaten(k,i)*mass(k,i)
  !     ! tend_qv(k,i) = tend_qv(k,i) + rqvfddaten(k,i)*mass(k,i)
  !     tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvfddaten(k,i)*mass(k,i)
  !  enddo
  !  enddo
endif


    !add coupled tendencies due to longwave radiation:
    if (config_radt_lw_scheme .ne. 'off') then
       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratenlw(k,i)*mass(k,i)
       enddo
       enddo
    endif
    
    !add coupled tendencies due to shortwave radiation:
    if (config_radt_sw_scheme .ne. 'off') then
       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratensw(k,i)*mass(k,i)
       enddo
       enddo
    endif

    !if non-hydrostatic core, convert the tendency for the potential temperature to a
    !tendency for the modified potential temperature:
    call physics_get_tend_phase4_kokkos_f(nCells, nVertLevels, nEdges, num_scalars, nCellsSolve,&
           tend_th,tend_scalars,&
           index_qv,theta_m,tend_theta,scalars)
   !  do i = 1, nCellsSolve
   !  do k = 1, nVertLevels
   !     coeff = (1. + R_v/R_d * scalars(index_qv,k,i))
   !     tend_th(k,i) = coeff * tend_th(k,i) + R_v/R_d * theta_m(k,i) * tend_scalars(index_qv,k,i) / coeff
   !     tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i)
   !  enddo
   !  enddo

 end subroutine physics_get_tend_work

!=================================================================================================================
 subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend)
!=================================================================================================================

 use iamas_atm_dimensions

!input arguments:
 type(block_type),intent(in),target:: block
 type(iamas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:),target:: Ux_tend,Uy_tend 

!output arguments:
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j
 integer:: cell1, cell2
 integer,pointer:: nCells,nCellsSolve,nEdges
 integer,dimension(:,:),pointer:: cellsOnEdge

 real(kind=RKIND),dimension(:,:),pointer:: Ux_tend_halo,Uy_tend_halo
 real(kind=RKIND), dimension(:,:), pointer :: east, north, edgeNormalVectors

 
!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_dimension(mesh, 'nCells', nCells)
 call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)

 call iamas_pool_get_array(mesh, 'east', east)
 call iamas_pool_get_array(mesh, 'north', north)
 call iamas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)

 call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
 Ux_tend_halo => Ux_tend
 Uy_tend_halo => Uy_tend

 tempField => tempFieldTarget 
 tempField % block => block
 tempField % dimSizes(1) = nVertLevels
 tempField % dimSizes(2) = nCellsSolve
 tempField % sendList => block % parinfo % cellsToSend
 tempField % recvList => block % parinfo % cellsToRecv
 tempField % copyList => block % parinfo % cellsToCopy
 tempField % prev => null()
 tempField % next => null()
 tempField % isActive = .true.

 tempField % array => Ux_tend_halo
 call iamas_dmpar_exch_halo_field(tempField)
 
 tempField % array => Uy_tend_halo
 call iamas_dmpar_exch_halo_field(tempField)
 call tend_toEdges_kokkos_f(nCells,nVertLevels,nCellsSolve,nEdges,Ux_tend,Uy_tend,U_tend,edgeNormalVectors,east,north,cellsOnEdge)
!  do iEdge = 1, nEdges
!     cell1 = cellsOnEdge(1,iEdge)
!     cell2 = cellsOnEdge(2,iEdge)
                           
!     U_tend(:,iEdge) =  Ux_tend_halo(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
!                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
!                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
!                      + Uy_tend_halo(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   &
!                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell1)   &
!                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell1))  &
!                      + Ux_tend_halo(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
!                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
!                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
!                      + Uy_tend_halo(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   &
!                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell2)   &
!                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell2))
!  end do
 
 end subroutine tend_toEdges

!=================================================================================================================
 end module iamas_atmphys_todynamics
!=================================================================================================================
