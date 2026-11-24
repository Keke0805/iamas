!************************************************************************
!
! * Author: Jiawang Feng and Chun Zhao from AEMOL USTC
! Last update: Dec 2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao
!
! References: 
!   Copied from WRF-Chem-v3.5:chem/module_mosaic_settling.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

! jwfeng: add diag of gravitational settling velocity 2022/09/02

MODULE module_mosaic_settling
use iamas_atmchem_configs
use iamas_atmchem_vars
use iamas_atmphys_vars
use iamas_atmchem_utilities
use iamas_log
!---------------------------------------------------------------------------------
! Purpose:
! calculate gravitational settling for MOSAIC aerosols 
!
! Author: Chun Zhao 
! v1 : 2016/05/30
!---------------------------------------------------------------------------------

interface
  subroutine f_gocart_settling_scheme(dt,t_phy,moist, &
    chem,rho_phy,dz8w,p8w,p_phy,dx,g,v_grav_3d, &
    its,ite,kts,kte,p_qv,ncomp_aer,massptr_aer,dens_aer,hygro_aer,dcen_sect, &
    num_moist,num_chem,nbin_a) bind(C,name="c_gocart_settling_scheme")
    use, intrinsic :: iso_c_binding
    use :: flcl_mod
    implicit none
    integer(c_int) :: its,ite,kts,kte,p_qv
    integer(c_int) :: num_moist,num_chem,nbin_a
    type(nd_array_t) :: moist
    type(nd_array_t) :: chem,v_grav_3d
    type(nd_array_t) :: t_phy,p_phy,dz8w,rho_phy
    type(nd_array_t) :: p8w
    type(nd_array_t) :: dx
    type(nd_array_t) :: ncomp_aer,massptr_aer,dens_aer,hygro_aer,dcen_sect
    real(c_double) :: dt,g
  end subroutine f_gocart_settling_scheme
end interface

CONTAINS

SUBROUTINE mosaic_settling_driver(dt,config_flags,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
         v_grav_3d, & !jwfeng add
         its,ite, kts,kte                   )

  IMPLICIT NONE

   TYPE(iAMAS_chem_config),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) ::                      &
                                  its,ite, kts,kte
    REAL, DIMENSION( 1:num_moist,kts:kte,its:ite ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( 1:num_chem,kts:kte,its:ite ),                 &
         INTENT(INOUT ) ::                                   chem
    REAL, DIMENSION( 1:num_chem,kts:kte,its:ite ),                 &
         INTENT(INOUT ) ::                                   v_grav_3d
   REAL,  DIMENSION( kts:kte+1,its:ite ),                        &
         INTENT(IN    ) :: p8w
   REAL,  DIMENSION( kts:kte,its:ite ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,rho_phy
  REAL, INTENT(IN   ) :: dt,g
   REAL,  DIMENSION( its:ite ), INTENT(IN) :: dx

!  INTEGER, PARAMETER :: settling_option=1

  settling_scheme: SELECT CASE(config_flags%settling_option)
   
  CASE (0)
      CALL chemistry_debug_message('no gravitational settling scheme called for mosaic aerosol')

  CASE (1)
      CALL chemistry_debug_message('gocart gravitational settling scheme called for mosaic aerosol')
      ! call gocart_settling_scheme(dt,t_phy,moist,  &
      !   chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
      !   v_grav_3d,&
      !   its,ite, kts,kte)
      call gocart_settling_scheme_flcl(dt,t_phy,moist,  &
          chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
          v_grav_3d,&
          its,ite, kts,kte)

  CASE (2)
      ! call cam_settling()

  CASE DEFAULT
  END SELECT settling_scheme


END SUBROUTINE mosaic_settling_driver

SUBROUTINE gocart_settling_scheme_flcl(dt,t_phy,moist,  &
  chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
  v_grav_3d, &
  its,ite, kts,kte                   )

use, intrinsic :: iso_c_binding
use :: flcl_mod
use module_data_binaer_asect, only: ncomp_aer,massptr_aer,hygro_aer,dens_aer,dcen_sect
use module_data_mosaic_therm, only: nbin_a
IMPLICIT NONE

INTEGER(c_int), INTENT(IN) :: its,ite, kts,kte
REAL(c_double), DIMENSION( 1:num_moist,kts:kte,its:ite ), INTENT(IN) :: moist
REAL(c_double), DIMENSION( 1:num_chem,kts:kte,its:ite ), INTENT(INOUT) :: chem
REAL(c_double), DIMENSION( 1:num_chem,kts:kte,its:ite ), INTENT(INOUT) :: v_grav_3d
REAL(c_double), DIMENSION( kts:kte,its:ite ), INTENT(IN) :: t_phy,p_phy,dz8w,rho_phy
REAL(c_double), DIMENSION( kts:kte+1,its:ite ), INTENT(IN) :: p8w
REAL(c_double), INTENT(IN) :: dt,g
REAL(c_double), DIMENSION( its:ite ), INTENT(IN) :: dx

call f_gocart_settling_scheme(dt,to_nd_array(t_phy),to_nd_array(moist), &
to_nd_array(chem),to_nd_array(rho_phy),to_nd_array(dz8w),to_nd_array(p8w), &
to_nd_array(p_phy),to_nd_array(dx),g,to_nd_array(v_grav_3d), &
its,ite,kts,kte,p_qv,to_nd_array(ncomp_aer),to_nd_array(massptr_aer), &
to_nd_array(dens_aer),to_nd_array(hygro_aer),to_nd_array(dcen_sect), &
num_moist,num_chem,nbin_a)

END SUBROUTINE gocart_settling_scheme_flcl


END MODULE module_mosaic_settling
