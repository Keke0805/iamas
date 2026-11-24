!************************************************************************
!
! Module to check if chemistry namelists are set up correctly.
! * Author: Jiawang Feng, Zihan Xia and Chun Zhao / AEMOL USTC
! Last update: Sep 2024
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
!   iamas_atmphys_control.F by Laura D. Fowler
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_control
 use iamas_dmpar
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines

 use iamas_atmchem_utilities

 implicit none
 private
 public:: chemistry_namelist_check

 contains


!=================================================================================================================
 subroutine chemistry_namelist_check(configs)
!=================================================================================================================

!===========================================================================
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + Check if chemistry namelists are set correctly. 
!                                         ---jwfeng   12/2019
!NOTE:
!    If you add a new chemistry namelist config in Registry.xml,
!    you'd better edit its "possible_values" in the same time.
!===========================================================================


!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_chem_opt,          &
                                  config_dust_opt,          &
                                  config_convection_scheme, &
                                  config_ctrans_opt,        &
                                  config_phot_opt,          &
                                  config_gaschem_opt,       &
                                  config_aerchem_sia_opt,   &
                                  config_aerchem_soa_opt,   &
                                  config_stra_gaschem_opt,  &
                                  config_cloudchem_opt
        
 integer,pointer ::               config_n2o5_hetchem !xyruan add                              

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine chemistry_namelist_check:'

 call iamas_pool_get_config(configs,'config_chem_opt'    ,config_chem_opt    )
 call iamas_pool_get_config(configs,'config_dust_opt'    ,config_dust_opt    )
 call iamas_pool_get_config(configs,'config_ctrans_opt'  ,config_ctrans_opt  )
 call iamas_pool_get_config(configs,'config_convection_scheme'  ,config_convection_scheme  )
 call iamas_pool_get_config(configs,'config_phot_opt', config_phot_opt  ) !zhxia add this for photolysis
 call iamas_pool_get_config(configs,'config_gaschem_opt', config_gaschem_opt  ) !zhxia add this for photochemistry
 call iamas_pool_get_config(configs,'config_aerchem_sia_opt', config_aerchem_sia_opt ) !xyruan add this for inorganic aerosol chemistry and interactions
 call iamas_pool_get_config(configs,'config_aerchem_soa_opt', config_aerchem_soa_opt ) !xyruan add this for organic aerosol chemistry and interactions
 call iamas_pool_get_config(configs,'config_stra_gaschem_opt', config_stra_gaschem_opt  ) !zhxia add this for stratospheric photochemistry
 call iamas_pool_get_config(configs,'config_cloudchem_opt', config_cloudchem_opt  ) !xyruan add this for cloud aqueous chemistry
 call iamas_pool_get_config(configs,'config_n2o5_hetchem', config_n2o5_hetchem  ) !xyruan add this for N2O5 heterogenous chemistry

!chem options:
 if(.not. (config_chem_opt .eq. 'no_chemistry'         .or. &
           config_chem_opt .eq. 'notIncluded'          .or. &
           config_chem_opt .eq. 'bin10aer_nature_aq'   .or. &
           config_chem_opt .eq. 'bin09aer_nature_aq'   .or. &
           config_chem_opt .eq. 'bin08aer_nature_aq'   .or. &
           config_chem_opt .eq. 'bin05aer_nature_aq'   .or. &
           config_chem_opt .eq. 'bin04aer_nature_aq'   .or. &
           config_chem_opt .eq. 'bin10aer_nature'      .or. &
           config_chem_opt .eq. 'bin09aer_nature'      .or. &
           config_chem_opt .eq. 'bin08aer_nature'      .or. &
           config_chem_opt .eq. 'bin05aer_nature'      .or. &
           config_chem_opt .eq. 'bin04aer_nature'      .or. &
           config_chem_opt .eq. 'bin10dust_aq'         .or. &
           config_chem_opt .eq. 'bin09dust_aq'         .or. &
           config_chem_opt .eq. 'bin08dust_aq'         .or. &
           config_chem_opt .eq. 'bin05dust_aq'         .or. &
           config_chem_opt .eq. 'bin04dust_aq'         .or. &
           config_chem_opt .eq. 'bin10dust'            .or. &
           config_chem_opt .eq. 'bin09dust'            .or. &
           config_chem_opt .eq. 'bin08dust'            .or. &
           config_chem_opt .eq. 'bin05dust'            .or. &
           config_chem_opt .eq. 'bin04dust'            .or. &
           config_chem_opt .eq. 'bin10seas_aq'         .or. &
           config_chem_opt .eq. 'bin09seas_aq'         .or. &
           config_chem_opt .eq. 'bin08seas_aq'         .or. &
           config_chem_opt .eq. 'bin05seas_aq'         .or. &
           config_chem_opt .eq. 'bin04seas_aq'         .or. &
           config_chem_opt .eq. 'bin10seas'            .or. &
           config_chem_opt .eq. 'bin09seas'            .or. &
           config_chem_opt .eq. 'bin08seas'            .or. &
           config_chem_opt .eq. 'bin05seas'            .or. &
           config_chem_opt .eq. 'bin04seas'            .or. &
           config_chem_opt .eq. 'bin04aer_gas_aq'      .or. &
           config_chem_opt .eq. 'bin04aer_gas'     )) then
    write(iamas_err_message,'(A,A20)') 'illegal value for config_chem_opt:', &
          trim(config_chem_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif

!dust scheme:
 if(.not. (config_dust_opt .eq. 'gocart'     .or. &
           config_dust_opt .eq. 'mosaic'     .or. &
           config_dust_opt .eq. 'Kok'        .or. &
           config_dust_opt .eq. 'shao_ver1_schem1' .or. &
           config_dust_opt .eq. 'shao_ver1_schem2' .or. &
           config_dust_opt .eq. 'shao_ver1_schem3' .or. &
           config_dust_opt .eq. 'shao_ver2_schem1' .or. &
           config_dust_opt .eq. 'shao_ver2_schem2' .or. &
           config_dust_opt .eq. 'shao_ver2_schem3' ))    then
    !ver1 means SHAO in wrfchem_v3.7
    !ver2 means SHAO in wrfchem_v4.0
    !schem1 means SHAO_2001
    !schem2 means SHAO_2004
    !schem3 means SHAO_2011  chenjin changed in 202312
    write(iamas_err_message,'(A,A10)') 'illegal value for config_dust_opt:', &
          trim(config_dust_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif

!ctrans scheme:
 if(.not. (config_ctrans_opt .eq. 'no_ctrans' .or. &
           config_ctrans_opt .eq. 'kfeta'     .or. &
           config_ctrans_opt .eq. 'general')) then

    write(iamas_err_message,'(A,A10)') 'illegal value for config_ctrans_opt:', &
          trim(config_ctrans_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif
!Check if ctrans_opt matches convection_scheme
 if(trim(config_ctrans_opt) .eq. 'kfeta' .and. trim(config_convection_scheme) .ne. 'cu_kain_fritsch') then
    iamas_err_message = 'config_ctrans_opt: kfeta only works when config_convection_scheme is cu_kain_fritsch'
    call chemistry_error_fatal(iamas_err_message)
 endif

!zhxia add: photolysis scheme:
 if(.not. (config_phot_opt .eq. 'no_phot'     .or. &
           config_phot_opt .eq. 'photfastj' ))    then
    write(iamas_err_message,'(A,A10)') 'illegal value for config_phot_opt:', &
          trim(config_phot_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif

!zhxia add: photochemistry scheme:
 if(.not. (config_gaschem_opt .eq. 'no_gaschem'     .or. &
           config_gaschem_opt .eq. 'AIPC_saprc99'   .or. &
           config_gaschem_opt .eq. 'KPP_saprc99'))    then
    write(iamas_err_message,'(A,A10)') 'illegal value for config_gaschem_opt:', &
          trim(config_gaschem_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif

!zhxia add: aerosol chemistry and interactions scheme:
!SIA
 if(.not. (config_aerchem_sia_opt .eq. 'no_sia'     .or. &
           config_aerchem_sia_opt .eq. 'AIMACI_mosaic'  .or. &
           config_aerchem_sia_opt .eq. 'simple_SIA'     .or. &
           config_aerchem_sia_opt .eq. 'MOSAIC_mech' ))    then
    write(iamas_err_message,'(A,A10)') 'illegal value for config_aerchem_sia_opt:', &
          trim(config_aerchem_sia_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif
!SOA
 if(.not. (config_aerchem_soa_opt .eq. 'no_soa'     .or. &
           config_aerchem_soa_opt .eq. 'simple_SOA' ))    then
    write(iamas_err_message,'(A,A10)') 'illegal value for config_aerchem_soa_opt:', &
          trim(config_aerchem_soa_opt)
    call chemistry_error_fatal(iamas_err_message)
 endif

!zhxia add: stratospheric photochemistry scheme:
 if(.not. (config_stra_gaschem_opt .eq. 'no_stra_gaschem'     .or. &
           config_stra_gaschem_opt .eq. 'linoz' ))    then
    write(iamas_err_message,'(A,A10)') 'illegal value for config_stra_gaschem_opt:', &
          trim(config_stra_gaschem_opt)
    call chemistry_error_fatal(iamas_err_message)

 endif

!cloud aqueous chemistry scheme:
 if(.not. (config_cloudchem_opt .eq. 'no_cloudchem'     .or. &
           config_cloudchem_opt .eq. 'CMU_cloudchem' ))    then
   write(iamas_err_message,'(A,A20)') 'illegal value for config_cloudchem_opt:', &
         trim(config_cloudchem_opt)
   call chemistry_error_fatal(iamas_err_message)
 endif

 if(trim(config_cloudchem_opt) .eq. 'CMU_cloudchem' .and. &
    trim(config_chem_opt) .ne. 'bin04aer_gas_aq') then
   iamas_err_message = 'ERROR: CMU_cloudchem currently requires config_chem_opt = bin04aer_gas_aq'
   call chemistry_error_fatal(iamas_err_message)
 endif

!xyruan add: N2O5 heterogenous chemistry
 if(.not. (config_n2o5_hetchem .eq. 0 .or. &
           config_n2o5_hetchem .eq. 1  )) then

    write(iamas_err_message,'(A,I0)') 'illegal value for config_n2o5_hetchem:', &
          config_n2o5_hetchem
    call chemistry_error_fatal(iamas_err_message)
 endif
!Check if config_n2o5_hetchem matches aerosol schemes
 if(config_n2o5_hetchem .eq. 1 .and. trim(config_aerchem_sia_opt) .ne. 'MOSAIC_mech') then
    iamas_err_message = 'ERROR: N2O5 heterogenous chemistry (without Cl- pathway) must be run with MOSAIC aerosol'
    call chemistry_error_fatal(iamas_err_message)
 endif
 end subroutine chemistry_namelist_check

!=================================================================================================================
 end module iamas_atmchem_control
!=================================================================================================================

