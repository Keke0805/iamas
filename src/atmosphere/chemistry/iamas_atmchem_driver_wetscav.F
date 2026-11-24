!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute wet scavenging tendencies for chem(:,:,:,:). 
! * Author: Jiawang Feng, Zihan Xia from AEMOL USTC
! Last update: jwfeng 03/2020
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_driver_wetscav
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_configs
 use iamas_atmchem_utilities
 use iamas_atmphys_vars, only: kts, kte
 use iamas_log

!chemistry schemes:
 use module_mosaic_wetscav

 implicit none
 private
 public :: driver_wetscav


REAL, PARAMETER :: mwso4 = 96.00   ! Molecular mass of SO4-- (g/mol)
REAL, PARAMETER :: mwno3 = 62.0    ! Molecular mass of NO3- (g/mol)

CONTAINS

subroutine driver_wetscav(state, diag_physics, time_lev, configs, tend_chemistry, its, ite)

   IMPLICIT NONE
   type(iamas_pool_type),intent(in):: diag_physics
   type(iamas_pool_type),intent(in):: configs
   type(iamas_pool_type),intent(inout):: state
   type(iamas_pool_type),intent(inout):: tend_chemistry
   integer,intent(in):: its,ite,time_lev
   real(kind=RKIND),dimension(:,:),pointer :: qlsink, ilsink
   real(kind=RKIND),dimension(:,:),pointer :: precr, preci, precs, precg
   real(kind=RKIND),dimension(:,:,:),pointer :: tend_wetls
   logical,pointer :: config_do_chem_budget

   !----------------------------------------------------------------------
   !
   !  wet removal by grid-resolved precipitation
   !  scavenging of cloud-phase aerosols and gases by collection, freezing, ...
   !  scavenging of interstitial-phase aerosols by impaction
   !  scavenging of gas-phase gases by mass transfer and reaction
   !
   !  This driver calls subroutines for wet scavenging.
   !       Only mosaic for now
   !
   !----------------------------------------------------------------------
   ! local vars
   real, dimension( kts:kte, its:ite ) :: r2evap, i2evap
   integer :: numgas_aqfrac
   call iamas_pool_get_config(configs,'config_do_chem_budget' ,config_do_chem_budget )
   !r2evap,i2evap are not important,not used further.
   r2evap(kts:kte,its:ite) = 0.0
   i2evap(kts:kte,its:ite) = 0.0
   numgas_aqfrac = num_gas

   ! copy iAMAS arrays to local arrays:
   call iamas_pool_get_array(diag_physics,'qlsink' ,qlsink )
   call iamas_pool_get_array(diag_physics,'ilsink' ,ilsink )
   call iamas_pool_get_array(diag_physics,'precr'  , precr )
   call iamas_pool_get_array(diag_physics,'preci'  , preci )
   call iamas_pool_get_array(diag_physics,'precs'  , precs )
   call iamas_pool_get_array(diag_physics,'precg'  , precg )

   if (config_do_chem_budget) then
      chem_old(1:num_chem,kts:kte,its:ite) = chem(1:num_chem,kts:kte,its:ite)
   end if

   call chemistry_debug_message('wetscav from iAMAS success')

   call wetscav_cbmz_mosaic (                                                          &
        id = 1,ktau = 1,dtstep = dt_dyn,ktauc = 1,configs = configs,         &
        dtstepc = dt_dyn, alt = alt(kts:kte,its:ite), t_phy = t_phy(kts:kte,its:ite),  &
        p8w = pres2_hyd(kts:kte+1,its:ite), t8w = t2(kts:kte+1,its:ite),               &
        p_phy = pres_hyd(kts:kte,its:ite),chem = chem(1:num_chem,kts:kte,its:ite),     &
        rho_phy = rho_phy(kts:kte,its:ite), cldfra = cldfrac(kts:kte,its:ite),         &
        qlsink = qlsink(kts:kte,its:ite), ilsink = ilsink(kts:kte,its:ite),            &
        precr = precr(kts:kte,its:ite),preci = preci(kts:kte,its:ite),                 &
        precs = precs(kts:kte,its:ite),precg = precg(kts:kte,its:ite),                 &
        r2evap = r2evap(kts:kte,its:ite),i2evap = i2evap(kts:kte,its:ite),             & 
        gas_aqfrac = gas_aqfrac(1:numgas_aqfrac,kts:kte,its:ite),                     &
        numgas_aqfrac = numgas_aqfrac,                                                &
        qv = qv(kts:kte,its:ite), qc = qc(kts:kte,its:ite),                            &
        qi = qi(kts:kte,its:ite), qs = qs(kts:kte,its:ite),                            &
        its = its, ite = ite, kts = kts, kte = kte)
      ! Calculate wet deposition of the species by accumulating their (negative)
      ! column change. It is assumed that all precipitation reaches the ground.

   call chemistry_debug_message( 'wetscav scheme success')

   if (config_do_chem_budget) then
      call iamas_pool_get_array(tend_chemistry,'tend_wetls',tend_wetls)
      tend_wetls(1:num_chem,kts:kte,its:ite) = (chem(1:num_chem,kts:kte,its:ite) - chem_old(1:num_chem,kts:kte,its:ite)) / dt_dyn
   end if

   call chemistry_debug_message( 'wetscav to iAMAS success')
end subroutine driver_wetscav

END MODULE iamas_atmchem_driver_wetscav
