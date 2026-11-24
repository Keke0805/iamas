!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute optical arguments for chem(:,:,:,:). 
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 04/2020
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

module iamas_atmchem_driver_optical
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_configs
 use iamas_atmchem_utilities
 use iamas_atmphys_vars

 use module_optical_averaging

 implicit none
 private
 public:: driver_optical

CONTAINS

subroutine driver_optical(curr_secs, diag, diag_chemistry, configs, its, ite)
    use  module_data_mosaic_therm, only: nbin_a
!------------------------------------------------------------------------
 IMPLICIT NONE
 real(kind=RKIND),intent(in   ):: curr_secs
 type(iamas_pool_type),intent(in   ):: diag
 type(iamas_pool_type),intent(in   ):: configs
 type(iamas_pool_type),intent(inout):: diag_chemistry
 integer,intent(in):: its,ite

!         
! local variables
!
 logical :: haveaer = .true.
 integer :: option_method, option_mie
 integer :: i,h,k
 character(len=StrKIND),pointer :: aer_op_opt, config_chem_opt, aerchem_sia_opt
 character(len=StrKIND):: msg
 real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw, absaersw, extaersw, gaersw, waersw, bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_dust, absaersw_dust, gaersw_dust, waersw_dust
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_seas, absaersw_seas, gaersw_seas, waersw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw, extaerlw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw_dust, tauaerlw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
 real(kind=RKIND),dimension(:,:  ),pointer :: aodsw, aodlw
 logical :: use_thermo_water  !xyruan-Controls whether the optical module uses aerosol water calculated by the thermodynamic scheme

 call iamas_pool_get_array(diag_chemistry,'tauaersw'     , tauaersw      )
 call iamas_pool_get_array(diag_chemistry,'absaersw'     , absaersw      )
 call iamas_pool_get_array(diag_chemistry,'gaersw'       , gaersw        )
 call iamas_pool_get_array(diag_chemistry,'waersw'       , waersw        )
 call iamas_pool_get_array(diag_chemistry,'aodsw'        , aodsw         )
 call iamas_pool_get_array(diag_chemistry,'aodlw'        , aodlw         )
 call iamas_pool_get_array(diag_chemistry,'tauaersw_dust', tauaersw_dust )
 call iamas_pool_get_array(diag_chemistry,'absaersw_dust', absaersw_dust )
 call iamas_pool_get_array(diag_chemistry,'gaersw_dust'  , gaersw_dust   )
 call iamas_pool_get_array(diag_chemistry,'waersw_dust'  , waersw_dust   )
 call iamas_pool_get_array(diag_chemistry,'tauaersw_seas', tauaersw_seas )
 call iamas_pool_get_array(diag_chemistry,'absaersw_seas', absaersw_seas )
 call iamas_pool_get_array(diag_chemistry,'gaersw_seas'  , gaersw_seas   )
 call iamas_pool_get_array(diag_chemistry,'waersw_seas'  , waersw_seas   )
 call iamas_pool_get_array(diag_chemistry,'extaersw'     , extaersw      )
 call iamas_pool_get_array(diag_chemistry,'bscoefsw'     , bscoefsw      )
 call iamas_pool_get_array(diag_chemistry,'tauaerlw'     , tauaerlw      )
 call iamas_pool_get_array(diag_chemistry,'tauaerlw_dust', tauaerlw_dust )
 call iamas_pool_get_array(diag_chemistry,'tauaerlw_seas', tauaerlw_seas )
 call iamas_pool_get_array(diag_chemistry,'extaerlw'     , extaerlw      )
 call iamas_pool_get_array(diag_chemistry,'l2aer'        , l2aer         )
 call iamas_pool_get_array(diag_chemistry,'l3aer'        , l3aer         )
 call iamas_pool_get_array(diag_chemistry,'l4aer'        , l4aer         )
 call iamas_pool_get_array(diag_chemistry,'l5aer'        , l5aer         )
 call iamas_pool_get_array(diag_chemistry,'l6aer'        , l6aer         )
 call iamas_pool_get_array(diag_chemistry,'l7aer'        , l7aer         )
 call iamas_pool_get_array(diag          ,'relhum'       ,relhum         )


 call iamas_pool_get_config(configs,'config_chem_opt',config_chem_opt)
 call iamas_pool_get_config(configs,'config_aer_op_opt',aer_op_opt)
 call iamas_pool_get_config(configs,'config_aerchem_sia_opt',aerchem_sia_opt)

 call chemistry_debug_message('optical averaging')

 call chemistry_debug_message("config_aer_op_opt = "//trim(aer_op_opt))

! xyruan 
 use_thermo_water = .false.

 select case (trim(aerchem_sia_opt))
 case ('MOSAIC_mech', 'AIMACI_mosaic')
    ! use water_a* calculated by the thermodynamic model
    use_thermo_water = .true. 
 case ('simple_SIA', 'no_sia')
    ! use diagnosed water from RH-kappa formula
    use_thermo_water = .false. 
 case default

 end select


 aer_op_opt_select: SELECT CASE(trim(aer_op_opt))
 CASE ('volume_approx')
     option_method=1
     option_mie=1
 CASE ('maxwell_approx')
     option_method=2
     option_mie=1
 CASE ('volume_exact')
     option_method=1
     option_mie=2
 CASE ('maxwell_exact')
     option_method=2
     option_mie=2
 CASE ('shell_exact')
     option_method=3
     option_mie=2
 CASE DEFAULT
     call chemistry_message('WARNING: Invalid aer_op_opt. Defaulting to volume_approx.')
     option_method=1
     option_mie=1
 END SELECT aer_op_opt_select

 if(haveaer) then
    ! call optical_from_iAMAS(diag, its, ite)
    tauaersw(1:nswbands,kts:kte,its:ite) = 0.0 
    absaersw(1:nswbands,kts:kte,its:ite) = 0.0 
    gaersw(1:nswbands,kts:kte,its:ite) = 0.0 
    waersw(1:nswbands,kts:kte,its:ite) = 0.0 
    tauaersw_dust(1:nswbands,kts:kte,its:ite) = 0.0 
    absaersw_dust(1:nswbands,kts:kte,its:ite) = 0.0 
    gaersw_dust(1:nswbands,kts:kte,its:ite) = 0.0 
    waersw_dust(1:nswbands,kts:kte,its:ite) = 0.0 
    tauaersw_seas(1:nswbands,kts:kte,its:ite) = 0.0 
    absaersw_seas(1:nswbands,kts:kte,its:ite) = 0.0 
    gaersw_seas(1:nswbands,kts:kte,its:ite) = 0.0 
    waersw_seas(1:nswbands,kts:kte,its:ite) = 0.0 
    extaersw(1:nswbands,kts:kte,its:ite) = 0.0
    bscoefsw(1:nswbands,kts:kte,its:ite) = 0.0
    extaerlw(1:nlwbands,kts:kte,its:ite) = 0.0 
    tauaerlw(1:nlwbands,kts:kte,its:ite) = 0.0
    tauaerlw_dust(1:nlwbands,kts:kte,its:ite) = 0.0
    tauaerlw_seas(1:nlwbands,kts:kte,its:ite) = 0.0
    l2aer(1:nswbands,kts:kte,its:ite) = 0.0
    l3aer(1:nswbands,kts:kte,its:ite) = 0.0
    l4aer(1:nswbands,kts:kte,its:ite) = 0.0
    l5aer(1:nswbands,kts:kte,its:ite) = 0.0
    l6aer(1:nswbands,kts:kte,its:ite) = 0.0
    l7aer(1:nswbands,kts:kte,its:ite) = 0.0
    call chemistry_debug_message('optical from iAMAS success')
    call chemistry_debug_message('optical driver: call optical averaging')

    !BSINGH(PNNL)- The followingvariables had an undefined behavior at the boundaries.
    call optical_averaging(                                                                                                                         &
         id = 1, curr_secs = curr_secs, dtstep = dt_dyn,                                                                                      &
         config_flags = config_flags, nbin_o = nbin_a,                                                                                              &
         haveaer = haveaer, option_method = option_method, option_mie = option_mie, use_thermo_water = use_thermo_water,                             &
         chem = chem(1:num_chem,kts:kte,its:ite), dz8w = dz(kts:kte,its:ite), alt = alt(kts:kte,its:ite), relhum = relhum(kts:kte,its:ite),&
      !  h2oai,h2oaj,h2oai_cz,h2oaj_cz,h2oac_cz,                      &
         tauaersw = tauaersw(1:nswbands,kts:kte,its:ite), absaersw = absaersw(1:nswbands,kts:kte,its:ite),                                          &
         extaersw = extaersw(1:nswbands,kts:kte,its:ite),   gaersw = gaersw(1:nswbands,kts:kte,its:ite),                                            &
         waersw = waersw(1:nswbands,kts:kte,its:ite),   bscoefsw = bscoefsw(1:nswbands,kts:kte,its:ite),                                            &
         l2aer = l2aer(1:nswbands,kts:kte,its:ite), l3aer = l3aer(1:nswbands,kts:kte,its:ite), l4aer = l4aer(1:nswbands,kts:kte,its:ite),           &
         l5aer = l5aer(1:nswbands,kts:kte,its:ite), l6aer = l6aer(1:nswbands,kts:kte,its:ite), l7aer = l7aer(1:nswbands,kts:kte,its:ite),           &
         tauaerlw = tauaerlw(1:nlwbands,kts:kte,its:ite), extaerlw = extaerlw(1:nlwbands,kts:kte,its:ite),                                          &
         tauaersw_dust=tauaersw_dust(1:nswbands,kts:kte,its:ite), absaersw_dust=absaersw_dust(1:nswbands,kts:kte,its:ite),                          &
         gaersw_dust  =gaersw_dust(1:nswbands,kts:kte,its:ite),  waersw_dust   =waersw_dust(1:nswbands,kts:kte,its:ite),                            &
         tauaersw_seas=tauaersw_seas(1:nswbands,kts:kte,its:ite), absaersw_seas=absaersw_seas(1:nswbands,kts:kte,its:ite),                          &
         gaersw_seas  =gaersw_seas(1:nswbands,kts:kte,its:ite),  waersw_seas   =waersw_seas(1:nswbands,kts:kte,its:ite),                            &
         tauaerlw_dust=tauaerlw_dust(1:nlwbands,kts:kte,its:ite), tauaerlw_seas=tauaerlw_seas(1:nlwbands,kts:kte,its:ite),                          &
        !  tauaersw_bc0,   absaersw_bc0,   gaersw_bc0,   waersw_bc0,    &
        !  tauaersw_sulf0, absaersw_sulf0, gaersw_sulf0, waersw_sulf0,  &
        !  tauaersw_oc0,   absaersw_oc0,   gaersw_oc0,   waersw_oc0,    &
        !  tauaersw_other0,absaersw_other0,gaersw_other0,waersw_other0, &
        !  tauaerlw_bc0,tauaerlw_dust0,tauaerlw_sulf0,                  &
        !  tauaerlw_oc0,tauaerlw_other0,                                &
        !  totoa_a01,totoa_a02,totoa_a03,totoa_a04,                     &
        !  totoa_a05,totoa_a06,totoa_a07,totoa_a08,                     &
         its=its, ite=ite, kts=kts, kte=kte                                                                                                         &
         )
    call chemistry_debug_message('optical averaging success')

    ! call optical_to_iAMAS(diag_chemistry, its, ite)
    aodsw(1:nswbands,its:ite) = 0.0
    do h = 1,nswbands
    do k = kts,kte
    do i = its,ite
        aodsw(h,i)                = aodsw(h,i) + tauaersw(h,k,i)
    end do
    end do
    end do
    write(msg,*) 'max/min of tauaersw      =',maxval(tauaersw(1:nswbands,kts:kte,its:ite)),minval(tauaersw(1:nswbands,kts:kte,its:ite))
    call chemistry_debug_message(msg)
    write(msg,*) 'max/min of tauaersw_dust =',maxval(tauaersw_dust(1:nswbands,kts:kte,its:ite)),minval(tauaersw_dust(1:nswbands,kts:kte,its:ite))
    call chemistry_debug_message(msg)
    write(msg,*) 'max/min of tauaersw_seas =',maxval(tauaersw_seas(1:nswbands,kts:kte,its:ite)),minval(tauaersw_seas(1:nswbands,kts:kte,its:ite))
    call chemistry_debug_message(msg)
    aodlw(1:nlwbands,its:ite) = 0.0
    do h = 1,nlwbands
    do k = kts,kte
    do i = its,ite
        aodlw(h,i)                = aodlw(h,i) + tauaerlw(h,k,i)
    end do
    end do
    end do
    call chemistry_debug_message('optical to iAMAS success')
 else
    call chemistry_message("No aerosols, Won't do optical_averaging")
 end if

end subroutine driver_optical

END MODULE iamas_atmchem_driver_optical
