! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!

!=================================================================================================================
 module iamas_atmphys_init_microphysics
 use iamas_dmpar
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_atmphys_constants, only: degrad

 use iamas_atmphys_utilities

 use iamas_log, only : iamas_log_write

!use module_mp_thompson, only: is_aerosol_aware,naCCN0,naCCN1,naIN0,naIN1,ntb_arc,ntb_arw,ntb_art,ntb_arr, &
!                              ntb_ark,tnccn_act
 use module_mp_thompson, only: is_aerosol_aware

 implicit none
 private
 public:: init_thompson_clouddroplets_foriAMAS

!iAMAS main initialization of the Thompson parameterization of cloud microphysics with nucleation of cloud
!droplets based on distributions of CCNs and INs (aerosol-aware parameterization).
!Laura D. Fowler (send comments to laura@ucar.edu).
!2016-03-28.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added "use iamas_dmpar" at the top of the module.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-04.


 contains


!=================================================================================================================
 subroutine init_thompson_clouddroplets_foriAMAS(mesh,configs,sfc_input,diag_physics)
!=================================================================================================================
!input variables:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: sfc_input

!inout variables:
 type(iamas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask
 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell

 real(kind=RKIND),dimension(:),pointer:: nt_c,mu_c

 logical,pointer:: config_prescribed_nc, config_aer_cloud_droplet_on
 character(len=StrKIND),pointer:: config_chem_opt

 integer:: iCell
 real(kind=RKIND):: xlat_loc, xlon_loc 

!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('--- enter subroutine init_thompson_clouddroplets_foriAMAS:')

 call iamas_pool_get_config(configs,'config_prescribed_nc',config_prescribed_nc)
!jwfeng
#ifdef DO_CHEMISTRY
 call iamas_pool_get_config(configs,'config_chem_opt',config_chem_opt)
 call iamas_pool_get_config(configs,'config_aer_cloud_droplet_on',config_aer_cloud_droplet_on)
#endif

 call iamas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call iamas_pool_get_array(mesh,'latCell',latCell)
 call iamas_pool_get_array(mesh,'lonCell',lonCell)

 call iamas_pool_get_array(sfc_input,'landmask',landmask)
 call iamas_pool_get_array(diag_physics,'nt_c',nt_c)
 call iamas_pool_get_array(diag_physics,'mu_c',mu_c)

!... initialize the prescribed number of cloud droplets, and mu_c (parameter in the exponential of the generalized
!gamma distribution) as a function of the land-cean mask. as set in the thompson cloud microphysics scheme, nt_c
!is set to 100 per cc (100.E6 m^-3) for maritime cases and 300 per cc (300.E6 m^-3) for continental cases.
!czhao add the control of regional nt_c

!jwfeng
 is_aerosol_aware = .false.

#ifdef DO_CHEMISTRY
 if(trim(config_chem_opt) .ne. "no_chemistry" .and. config_aer_cloud_droplet_on) then
    is_aerosol_aware = .true.
 else
    is_aerosol_aware = .false.
 endif
#endif

 do iCell = 1, nCellsSolve

   xlat_loc = latCell(iCell) / degrad
   xlon_loc = lonCell(iCell) / degrad

   if(landmask(iCell) .eq. 1) then
      nt_c(iCell) = 200.e6
   elseif(landmask(iCell) .eq. 0) then
      nt_c(iCell) = 50.e6
   endif

   !overwrite the NC number if prescribed for some regions
   if (config_prescribed_nc) then 

    if (xlat_loc.ge.10.and.xlat_loc.le.45.and.   &
        xlon_loc.ge.100.and.xlon_loc.le.150) then  !East Asia
     if(landmask(iCell) .eq. 1) then
        nt_c(iCell) = 1000.e6 
     elseif(landmask(iCell) .eq. 0) then
        nt_c(iCell) = 300.e6
     endif

!!   call iamas_log_write('')
!!   call iamas_log_write('--- czhao prescribe thompson clouddroplets for East Asia')
!!   call iamas_log_write('--- iCell = $i',intArgs=(/iCell/))
!!   call iamas_log_write('--- landmask = $i',intArgs=(/landmask(iCell)/))
!!   call iamas_log_write('--- NT_C = $r',realArgs=(/nt_c(iCell)/))
!!   call iamas_log_write('')

    elseif (xlat_loc.ge.5.and.xlat_loc.le.30.and.   &
            xlon_loc.ge.68.and.xlon_loc.le.90) then  !South Asia
     if(landmask(iCell) .eq. 1) then
        nt_c(iCell) = 1000.e6 
     elseif(landmask(iCell) .eq. 0) then
        nt_c(iCell) = 300.e6
     endif

!!   call iamas_log_write('')
!!   call iamas_log_write('--- czhao prescribe thompson clouddroplets for South Asia')
!!   call iamas_log_write('--- NT_C = $r',realArgs=(/nt_c(iCell)/))
!!   call iamas_log_write('')

    endif 

   endif 

   ! jwfeng: just for test, delete later
!    nt_c(iCell) = 250.e6
    mu_c(iCell) = MIN(15., (1000.e6/nt_c(iCell) + 2.))
 enddo

!call iamas_log_write('--- end subroutine init_thompson_clouddroplets_foriAMAS.')

 end subroutine init_thompson_clouddroplets_foriAMAS 

!=================================================================================================================
 end module iamas_atmphys_init_microphysics
!=================================================================================================================

 
 
