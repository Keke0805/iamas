! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_driver_oml
 use iamas_kind_types
 use iamas_pool_routines

 use iamas_atmphys_constants
 use iamas_atmphys_landuse
 use iamas_atmphys_lsm_noahinit
 use iamas_atmphys_vars
 use iamas_constants

!wrf physics
 use module_sf_oml
 
 implicit none
 private
 public :: driver_oml1d

 integer,private:: i,j,k,n

!iAMAS driver for 1d ocean mixed layer
!Bill Skamarock (send comments to skamaroc@ucar.edu).
!2014-06-11.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added update of the sea-surface temperature below call to subroutine oml1d.
!   Laura D. Fowler (laura@ucar.edu) / 2014-07-24.


 contains


!=================================================================================================================
 subroutine driver_oml1d(configs,mesh,diag,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: diag

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: sfc_input


  logical,pointer:: config_oml1d
  real(kind=RKIND),pointer:: oml_gamma
  real(kind=RKIND),pointer:: oml_relaxation_time

  real(kind=RKIND), dimension(:), pointer:: t_oml
  real(kind=RKIND), dimension(:), pointer:: t_oml_initial
  real(kind=RKIND), dimension(:), pointer:: t_oml_200m_initial
  real(kind=RKIND), dimension(:), pointer:: h_oml
  real(kind=RKIND), dimension(:), pointer:: h_oml_initial
  real(kind=RKIND), dimension(:), pointer:: hu_oml
  real(kind=RKIND), dimension(:), pointer:: hv_oml

  real(kind=RKIND), dimension(:), pointer:: sst
  real(kind=RKIND), dimension(:), pointer:: skintemp
  real(kind=RKIND), dimension(:), pointer:: hfx
  real(kind=RKIND), dimension(:), pointer:: lh
  real(kind=RKIND), dimension(:), pointer:: gsw
  real(kind=RKIND), dimension(:), pointer:: glw
  real(kind=RKIND), dimension(:), pointer:: sfc_emiss
  real(kind=RKIND), dimension(:), pointer:: ust
  real(kind=RKIND), dimension(:), pointer:: xland

  real(kind=RKIND), dimension(:,:), pointer:: uReconstructZonal
  real(kind=RKIND), dimension(:,:), pointer:: uReconstructmeridional
  real(kind=RKIND), dimension(:)  , pointer:: latCell

  integer, pointer:: nCells, nCellsSolve

  integer :: iCell
  !qkyang changed to support Shin-Hong PBL>>>>>>>>>>>>>>>>
  ! real(kind=RKIND):: f_coriolis
  real(kind=RKIND),dimension(:), pointer:: f_coriolis
  !qkyang changed to support Shin-Hong PBL<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('--- enter subroutine driver_oml1d:')

!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>>>>>
call iamas_pool_get_array(diag_physics,'f_coriolis',f_coriolis)
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<<<<<

! namelist parameters for ocean mixed layer model
 call iamas_pool_get_config(configs,'config_oml1d'               ,config_oml1d )
 call iamas_pool_get_config(configs,'config_oml_gamma'           ,oml_gamma )
 call iamas_pool_get_config(configs,'config_oml_relaxation_time' ,oml_relaxation_time )

! state and initial state for ocean mixed layer model
 call iamas_pool_get_array(diag_physics,'t_oml'              ,t_oml )
 call iamas_pool_get_array(diag_physics,'t_oml_initial'      ,t_oml_initial )
 call iamas_pool_get_array(diag_physics,'t_oml_200m_initial' ,t_oml_200m_initial )
 call iamas_pool_get_array(diag_physics,'h_oml'              ,h_oml )
 call iamas_pool_get_array(diag_physics,'h_oml_initial'      ,h_oml_initial )
 call iamas_pool_get_array(diag_physics,'hu_oml'             ,hu_oml )
 call iamas_pool_get_array(diag_physics,'hv_oml'             ,hv_oml )

! state and diagnostics from iamas
 call iamas_pool_get_array(diag_physics,'hfx'       ,hfx )
 call iamas_pool_get_array(diag_physics,'lh'        ,lh )
 call iamas_pool_get_array(diag_physics,'gsw'       ,gsw )
 call iamas_pool_get_array(diag_physics,'glw'       ,glw )
 call iamas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call iamas_pool_get_array(diag_physics,'ust'       ,ust )

 call iamas_pool_get_array(sfc_input,'sst'       ,sst )
 call iamas_pool_get_array(sfc_input,'skintemp'  ,skintemp )
 call iamas_pool_get_array(sfc_input,'xland'     ,xland )

 call iamas_pool_get_array(diag,'uReconstructZonal'      ,uReconstructZonal )
 call iamas_pool_get_array(diag,'uReconstructMeridional' ,uReconstructMeridional )

 call iamas_pool_get_array(mesh,'latCell'      ,latCell )

 call iamas_pool_get_dimension(mesh,'nCells'       ,nCells )
 call iamas_pool_get_dimension(mesh,'nCellsSolve'  ,nCellsSolve )

 if (config_oml1d) then
 
   do iCell = 1, nCellsSolve

! if ocean point, call the 1d ocean mixed layer model 
     if( xland(iCell) .gt. 1.5) then
     !qkyang changed to support Shin-Hong PBL>>>>>>>>>>>>>
      !  f_coriolis = 2.*omega*cos(latCell(iCell))
       f_coriolis(iCell) = 2.*omega*cos(latCell(iCell))
       call oml1d( t_oml(iCell), t_oml_initial(iCell), h_oml(iCell), h_oml_initial(iCell),  &
                   hu_oml(iCell), hv_oml(iCell), skintemp(iCell), hfx(iCell),               &
                   lh(iCell), gsw(iCell), glw(iCell), t_oml_200m_initial(iCell),            &
                   uReconstructZonal(1,iCell), uReconstructMeridional(1,iCell),             &
                   ust(iCell), f_coriolis(iCell), sfc_emiss(iCell), stbolt, gravity, dt_pbl,       &
                   oml_gamma, oml_relaxation_time                                          )
      !qkyang changed to support Shin-Hong PBL<<<<<<<<<<<<
       sst(iCell) = skintemp(iCell)

     end if

   end do

 end if

!call iamas_log_write('--- end subroutine driver_oml1d.')

 end subroutine driver_oml1d

!=================================================================================================================
 end module iamas_atmphys_driver_oml
!=================================================================================================================
