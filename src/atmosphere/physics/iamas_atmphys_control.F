! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_control
 use iamas_dmpar
 use iamas_kind_types
 use iamas_pool_routines

 use iamas_atmphys_utilities
 use iamas_atmphys_vars, only: l_mp_tables

 implicit none
 private
 public:: physics_namelist_check, &
          physics_registry_init,  &
          physics_tables_init,    &
          physics_compatibility_check

 logical,public:: moist_physics


!iAMAS control and initialization routines.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines called in iamas_atmphys_control:
! -------------------------------------------
! * physics_namelist_check: checks that physics namelist parameters are defined correctly.
! * physics_registry_init : initializes thicknesses of soil layers for NOAH land-surface scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the namelist option config_eddy_scheme and associated sourcecode.
! * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
! * removed controls to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * removed subroutine physics_idealized_init, also available in iamas_init_atm_cases.F in core_init_atmosphere.
!   Laura D. Fowler (laura@ucar.edu) / 2014-08-11.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "wsm6" to "mp_wsm6" and "kessler" to "mp_kessler".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-09.
! * renamed "kain_fritsch" to "cu_kain_fritsch" and "tiedtke" to "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * renamed "ysu" to "bl_ysu", "ysu_gwdo" to "bl_gwdo_ysu", and "monin_obukhov" to "sf_monin_obukhov".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25. 
! * added the option mp_thompson.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the option cu_grell_freitas.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-31.
! * added the options sf_mynn and bl_mynn and for the MYNN parameterization from WRF version 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * added the option cu_ntiedtke for the "new" Tiedtke parameterization of convection from WRF version 3.8.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * added the physics suite "convection_scale_aware" (see below for the physics options used in the suite).
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-28.
! * added the subroutine physics_tables_init which checks if the files containing the lokk-up tables for the
!   Thompson cloud microphysics are available or not.
!   Laura D. Fowler (laura@ucar.edu) / 2016-11-01.
! * modified checking the config_gwdo_scheme option to allow bl_ysu_gwdo to be run when the MYNN pbl and surface
!   layer scheme options are chosen.
!   Laura D. Fowler (laura@ucar.edu) / 2016-12-22.
! * modified logic in subroutine physics_tables_init so that the Thompson microphysics tables are read in each
!   MPI task.
!   Laura D. Fowler (laura@ucar.edu) / 2016-12-30.
! * In the subroutine physics_tables_init which checks if the files containing the lokk-up tables for the
!   P3 cloud microphysics are available or not.
!   J.Gu (jungu@ustc.edu.cn) / 2024-06-02


 contains


!=================================================================================================================
 subroutine physics_namelist_check(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_physics_suite,     &
                                  config_microp_scheme,     &
                                  config_mp_p3_scheme,      &
                                  config_convection_scheme, &
                                  config_lsm_scheme,        &
                                  config_pbl_scheme,        &
                                  config_gwdo_scheme,       &
                                  config_radt_cld_scheme,   &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme,    &
                                  config_sfclayer_scheme,   &
                                  config_fdda_scheme,       &
                                  config_les_scheme,        &
                                  config_horiz_mixing
 real(kind=RKIND),pointer:: config_h_mom_eddy_visc2,   &
                            config_h_mom_eddy_visc4,   &
                            config_v_mom_eddy_visc2,   &
                            config_h_theta_eddy_visc2, &
                            config_h_theta_eddy_visc4, &
                            config_v_theta_eddy_visc2

!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('--- enter subroutine physics_namelist_check:')

 call iamas_pool_get_config(configs,'config_physics_suite'    ,config_physics_suite    )
 call iamas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call iamas_pool_get_config(configs,'config_mp_p3_scheme'     ,config_mp_p3_scheme    )
 call iamas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call iamas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call iamas_pool_get_config(configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call iamas_pool_get_config(configs,'config_gwdo_scheme'      ,config_gwdo_scheme      )
 call iamas_pool_get_config(configs,'config_radt_cld_scheme'  ,config_radt_cld_scheme  )
 call iamas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call iamas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )
 call iamas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call iamas_pool_get_config(configs,'config_fdda_scheme'      ,config_fdda_scheme      )
 call iamas_pool_get_config(configs,'config_les_scheme'       ,config_les_scheme       )
 call iamas_pool_get_config(configs,'config_horiz_mixing'     ,config_horiz_mixing     )
 call iamas_pool_get_config(configs,'config_h_mom_eddy_visc2' ,config_h_mom_eddy_visc2 )
 call iamas_pool_get_config(configs,'config_h_mom_eddy_visc4' ,config_h_mom_eddy_visc4 )
 call iamas_pool_get_config(configs,'config_v_mom_eddy_visc2' ,config_v_mom_eddy_visc2 )
 call iamas_pool_get_config(configs,'config_h_theta_eddy_visc2',config_h_theta_eddy_visc2)
 call iamas_pool_get_config(configs,'config_h_theta_eddy_visc4',config_h_theta_eddy_visc4)
 call iamas_pool_get_config(configs,'config_v_theta_eddy_visc2',config_v_theta_eddy_visc2)

 call iamas_log_write('')
 call iamas_log_write('----- Setting up physics suite '''//trim(config_physics_suite)//''' -----')

 !
 !setup schemes according to the selected physics suite:
 !
 if (trim(config_physics_suite) == 'mesoscale_reference') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'mp_wsm6'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'cu_ntiedtke'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'bl_ysu'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'bl_ysu_gwdo'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'rrtmg_lw'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'rrtmg_sw'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'cld_fraction'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'sf_monin_obukhov'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'noah'

 else if (trim(config_physics_suite) == 'convection_permitting') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'mp_thompson'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'cu_grell_freitas'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'bl_mynn'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'bl_ysu_gwdo'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'rrtmg_lw'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'rrtmg_sw'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'cld_fraction'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'sf_mynn'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'noah'

 else if (trim(config_physics_suite) == 'none') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'off'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'off'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'off'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'off'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'off'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'off'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'off'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'off'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'off'

 else

    write(iamas_err_message,'(A)') 'Unrecognized choice of physics suite: config_physics_suite = '''// &
                                  trim(config_physics_suite)//''''
    call physics_error_fatal(iamas_err_message)

 end if

!fdda scheme:
 if(.not. (config_fdda_scheme .eq. 'off'        .or. &
           config_fdda_scheme .eq. 'analysis'   .or. &
           config_fdda_scheme .eq. 'scaled'     )) then

    write(iamas_err_message,'(A,A20)') 'illegal value for config_fdda_scheme:', &
          trim(config_fdda_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!cloud microphysics scheme:
 if(.not. (config_microp_scheme .eq. 'off'         .or. &
           config_microp_scheme .eq. 'mp_kessler'  .or. &
           config_microp_scheme .eq. 'mp_thompson' .or. &
           config_microp_scheme .eq. 'mp_morrison' .or. &
           config_microp_scheme .eq. 'mp_p3'       .or. &
           config_microp_scheme .eq. 'mp_wsm6')) then

    write(iamas_err_message,'(A,A10)') 'illegal value for config_microp_scheme:', &
          trim(config_microp_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!P3 cloud microphysics scheme:
 if(config_microp_scheme .eq. 'mp_p3' .and. &
         .not. (config_mp_p3_scheme .eq. 'P3_1CATEGORY' .or. &
         config_mp_p3_scheme .eq. 'P3_1CATEGORY_NC'  .or. &
         config_mp_p3_scheme .eq. 'P3_2CATEGORY'     .or. &
         config_mp_p3_scheme .eq. 'P3_1CAT_3MOM')) then

     write(iamas_err_message,'(A,A10)') 'illegal value for config_mp_p3_scheme:', &
         trim(config_mp_p3_scheme)
     call physics_error_fatal(iamas_err_message)

 endif

!convection scheme:
 if(.not. (config_convection_scheme .eq. 'off'                 .or. &
           config_convection_scheme .eq. 'cu_grell_freitas'    .or. &
           config_convection_scheme .eq. 'cu_kain_fritsch'     .or. &
           config_convection_scheme .eq. 'cu_kain_fritsch_sca' .or. &
           config_convection_scheme .eq. 'cu_tiedtke'          .or. &
           config_convection_scheme .eq. 'cu_ntiedtke')) then

    write(iamas_err_message,'(A,A10)') 'illegal value for config_convection_scheme: ', &
          trim(config_convection_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!pbl scheme:
 if(.not. (config_pbl_scheme .eq. 'off'     .or. &
           config_pbl_scheme .eq. 'bl_mynn' .or. &
           config_pbl_scheme .eq. 'bl_ysu'  .or. &
           config_pbl_scheme .eq. 'bl_shinhong')) then!qkyang added 'bl_shinhong' to support Shin-Hong PBL.

    write(iamas_err_message,'(A,A15)') 'illegal value for pbl_scheme: ', &
          trim(config_pbl_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

 if (trim(config_les_scheme) == 'les_deardorff' .and. trim(config_pbl_scheme) /= 'off') then
    call physics_error_fatal('config_les_scheme=les_deardorff requires config_pbl_scheme=off')
 endif

 if (trim(config_les_scheme) == 'les_deardorff') then
    if (trim(config_horiz_mixing) /= '2d_fixed') then
       call physics_error_fatal('LES requires config_horiz_mixing=2d_fixed to disable native Smagorinsky mixing')
    endif

    if (config_h_mom_eddy_visc2 /= 0.0_RKIND .or. &
        config_h_mom_eddy_visc4 /= 0.0_RKIND .or. &
        config_v_mom_eddy_visc2 /= 0.0_RKIND .or. &
        config_h_theta_eddy_visc2 /= 0.0_RKIND .or. &
        config_h_theta_eddy_visc4 /= 0.0_RKIND .or. &
        config_v_theta_eddy_visc2 /= 0.0_RKIND) then
       call physics_error_fatal('LES requires all native momentum/theta eddy-viscosity coefficients to be zero')
    endif
 endif

!gravity wave drag over orography scheme:
 if(.not. (config_gwdo_scheme .eq. 'off' .or. &
           config_gwdo_scheme .eq. 'bl_ysu_gwdo')) then

    write(iamas_err_message,'(A,A10)') 'illegal value for gwdo_scheme: ', &
          trim(config_gwdo_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!lw radiation scheme:
 if(.not. (config_radt_lw_scheme .eq. 'off'    .or. &
           config_radt_lw_scheme .eq. 'cam_lw' .or. &
           config_radt_lw_scheme .eq. 'rrtmg_lw')) then
 
    write(iamas_err_message,'(A,A10)') 'illegal value for longwave radiation scheme: ', &
          trim(config_radt_lw_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!sw radiation scheme:
 if(.not. (config_radt_sw_scheme .eq. 'off'    .or. &
           config_radt_sw_scheme .eq. 'cam_sw' .or. &
           config_radt_sw_scheme .eq. 'rrtmg_sw')) then
 
    write(iamas_err_message,'(A,A10)') 'illegal value for shortwave radiation _scheme: ', &
          trim(config_radt_sw_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!cloud fraction for radiation schemes:
 if(.not. (config_radt_cld_scheme .eq. 'off'           .or. &
           config_radt_cld_scheme .eq. 'cld_incidence' .or. &
           config_radt_cld_scheme .eq. 'cld_fraction'  .or. &
           config_radt_cld_scheme .eq. 'cld_fraction_tompkins'  .or. &
           config_radt_cld_scheme .eq. 'cld_fraction_thompson')) then

    write(iamas_err_message,'(A,A10)') 'illegal value for calculation of cloud fraction: ', &
          trim(config_radt_cld_scheme)
    call physics_error_fatal(iamas_err_message)

 endif
 if((config_radt_lw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off') .or. &
    (config_radt_sw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off')) then

    call iamas_log_write('')
    write(iamas_err_message,'(A,A10)') &
       '    config_radt_cld_scheme is not set for radiation calculation'
    call physics_message(iamas_err_message)
    write(iamas_err_message,'(A,A10)') &
       '    switch calculation of cloud fraction to config_radt_cld_scheme = cld_incidence'
    call physics_message(iamas_err_message)
    config_radt_cld_scheme = "cld_incidence"

 endif

!surface-layer scheme:
 if(.not. (config_sfclayer_scheme .eq. 'off'     .or. &
           config_sfclayer_scheme .eq. 'sf_mynn' .or. &
           config_sfclayer_scheme .eq. 'sf_monin_obukhov')) then
 
    write(iamas_err_message,'(A,A10)') 'illegal value for surface layer scheme: ', &
          trim(config_sfclayer_scheme)
    call physics_error_fatal(iamas_err_message)
 else
    if(config_pbl_scheme == 'bl_mynn') then
       config_sfclayer_scheme = 'sf_mynn'
    elseif(config_pbl_scheme == 'bl_ysu') then
       config_sfclayer_scheme = 'sf_monin_obukhov'
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
    elseif(config_pbl_scheme == 'bl_shinhong') then
       config_sfclayer_scheme = 'sf_monin_obukhov'
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<
    endif
 endif

!land-surface scheme: note that config_sfclayer_scheme must be defined for the land-surface
!scheme to be called:
 if(config_lsm_scheme .ne. 'off' .and. config_sfclayer_scheme .eq. 'off') then
 
    call physics_error_fatal('land surface scheme: ' // &
                             'set config_sfclayer_scheme different than off')
    
 elseif(.not. (config_lsm_scheme .eq. 'off ' .or. &
               config_lsm_scheme .eq. 'noah')) then
 
    write(iamas_err_message,'(A,A10)') 'illegal value for land surface scheme: ', &
          trim(config_lsm_scheme)
    call physics_error_fatal(iamas_err_message)

 endif

!checks if any physics process is called. if not, return:
 moist_physics = .true.
 
 if(config_microp_scheme     .eq. 'off' .and. &
    config_convection_scheme .eq. 'off' .and. &
    config_lsm_scheme        .eq. 'off' .and. &
    config_pbl_scheme        .eq. 'off' .and. &
    config_radt_lw_scheme    .eq. 'off' .and. &
    config_radt_sw_scheme    .eq. 'off' .and. &
    config_sfclayer_scheme   .eq. 'off') moist_physics = .false.

 call iamas_log_write('')
 call iamas_log_write('    config_microp_scheme       = '//trim(config_microp_scheme))
 if (config_microp_scheme .eq. 'mp_p3') call iamas_log_write('    config_mp_p3_scheme       = '//trim(config_mp_p3_scheme))
 call iamas_log_write('    config_convection_scheme   = '//trim(config_convection_scheme))
 call iamas_log_write('    config_pbl_scheme          = '//trim(config_pbl_scheme))
 call iamas_log_write('    config_gwdo_scheme         = '//trim(config_gwdo_scheme))
 call iamas_log_write('    config_radt_cld_scheme     = '//trim(config_radt_cld_scheme))
 call iamas_log_write('    config_radt_lw_scheme      = '//trim(config_radt_lw_scheme))
 call iamas_log_write('    config_radt_sw_scheme      = '//trim(config_radt_sw_scheme))
 call iamas_log_write('    config_sfclayer_scheme     = '//trim(config_sfclayer_scheme))
 call iamas_log_write('    config_lsm_scheme          = '//trim(config_lsm_scheme))
 call iamas_log_write('    config_fdda_scheme         = '//trim(config_fdda_scheme))
 call iamas_log_write('')

 end subroutine physics_namelist_check

!=================================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs

!inout arguments:
 type(iamas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: config_lsm_scheme
 integer,pointer:: nCells
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:,:),pointer:: dzs

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call iamas_pool_get_config(configs,'config_lsm_scheme',config_lsm_scheme)

 call iamas_pool_get_dimension(mesh,'nCells',nCells)

 call iamas_pool_get_array(sfc_input,'landmask',landmask)
 call iamas_pool_get_array(sfc_input,'dzs'     , dzs    )

!initialization of input variables, if needed:

 if(.not. config_do_restart) then

    lsm_select: select case(trim(config_lsm_scheme))

       case("noah")
       !initialize the thickness of the soil layers for the Noah scheme:
          do iCell = 1, nCells
             dzs(1,iCell) = 0.10_RKIND
             dzs(2,iCell) = 0.30_RKIND
             dzs(3,iCell) = 0.60_RKIND
             dzs(4,iCell) = 1.00_RKIND
          enddo

       case default
    
    end select lsm_select
    
 endif

!call iamas_log_write('--- enter subroutine physics_namelist_check.')
!call iamas_log_write('')

 end subroutine physics_registry_init

!=================================================================================================================
 subroutine physics_tables_init(dminfo,configs)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(iamas_pool_type),intent(in):: configs

!local variables:
 character(len=StrKIND),pointer:: config_microp_scheme
 logical:: l_qr_acr_qg,l_qr_acr_qs,l_qi_aut_qs,l_freezeH2O
 logical:: l_1cat_2mom,l_1cat_3mom,l_2cat_2mom !for P3

!-----------------------------------------------------------------------------------------------------------------

 l_mp_tables = .true.

 if(dminfo % my_proc_id == IO_NODE) then

    call iamas_pool_get_config(configs,'config_microp_scheme',config_microp_scheme)
    mp_select: select case(trim(config_microp_scheme))

        case("mp_thompson")
            l_qr_acr_qg = .false.
            l_qr_acr_qs = .false.
            l_qi_aut_qs = .false.
            l_freezeH2O = .false.

            inquire(file='MP_THOMPSON_QRacrQG_DATA.DBL'  ,exist=l_qr_acr_qg)
            inquire(file='MP_THOMPSON_QRacrQS_DATA.DBL'  ,exist=l_qr_acr_qs)
            inquire(file='MP_THOMPSON_QIautQS_DATA.DBL'  ,exist=l_qi_aut_qs)
            inquire(file='MP_THOMPSON_freezeH2O_DATA.DBL',exist=l_freezeH2O)

        !   call iamas_log_write('')
        !   call iamas_log_write('--- enter subroutine physics_tables_init:')
        !   call iamas_log_write('l_qr_acr_qg = $l',logicArgs=(/l_qr_acr_qg/))
        !   call iamas_log_write('l_qr_acr_qs = $l',logicArgs=(/l_qr_acr_qs/))
        !   call iamas_log_write('l_qi_aut_qs = $l',logicArgs=(/l_qi_aut_qs/))
        !   call iamas_log_write('l_freezeH2O = $l',logicArgs=(/l_freezeH2O/))

            if(.not. (l_qr_acr_qg .and. l_qr_acr_qs .and. l_qi_aut_qs .and. l_freezeH2O)) then
            write(iamas_err_message,'(A)') &
                '--- tables to run the Thompson cloud microphysics scheme do not exist: run build_tables first.'
            call physics_error_fatal(iamas_err_message)
            endif
        !   call iamas_log_write('l_mp_tables = $l',logicArgs=(/l_mp_tables/))
        case("mp_p3")
            l_1cat_2mom = .false.
            l_1cat_3mom = .false.
            l_2cat_2mom = .false.

            inquire(file='p3_lookupTable_1.dat-v5.4_2momI'  ,exist=l_1cat_2mom)
            inquire(file='p3_lookupTable_1.dat-v5.4_3momI'  ,exist=l_1cat_3mom)
            inquire(file='p3_lookupTable_2.dat-v5.3'  ,exist=l_2cat_2mom)

        !   call iamas_log_write('')
        !   call iamas_log_write('--- enter subroutine physics_tables_init:')
        !   call iamas_log_write('l_qr_acr_qg = $l',logicArgs=(/l_qr_acr_qg/))
        !   call iamas_log_write('l_qr_acr_qs = $l',logicArgs=(/l_qr_acr_qs/))
        !   call iamas_log_write('l_qi_aut_qs = $l',logicArgs=(/l_qi_aut_qs/))
        !   call iamas_log_write('l_freezeH2O = $l',logicArgs=(/l_freezeH2O/))

            if(.not. (l_1cat_2mom .and. l_1cat_3mom .and. l_2cat_2mom )) then
            write(iamas_err_message,'(A)') &
                '--- tables to run the P3 cloud microphysics scheme do not exist'
            call physics_error_fatal(iamas_err_message)
            endif
        !   call iamas_log_write('l_mp_tables = $l',logicArgs=(/l_mp_tables/))
        case default
    end select mp_select

 endif

 end subroutine physics_tables_init

!=================================================================================================================
!  routine physics_compatibility_check()
!
!> \brief Checks physics input fields and options for compatibility
!> \author Miles Curry and Michael Duda
!> \date   25 October 2018
!> \details
!>  This routine checks the input fields and run-time options provided
!>  by the user for compatibility. For example, two run-time options may
!>  be mutually exclusive, or an option may require that a certain input
!>  field is provided. The checks performed by this routine are only for
!>  physics related fields and options.
!>
!>  A value of 0 is returned if there are no incompatibilities among
!>  the provided input fields and run-time options, and a non-zero value
!>  otherwise.
!>
 subroutine physics_compatibility_check(dminfo, blockList, streamManager, ierr)
!=================================================================================================================

 implicit none

 type (dm_info), pointer :: dminfo
 type (block_type), pointer :: blockList
 type (iAMAS_streamManager_type), pointer :: streamManager
 integer, intent(out) :: ierr

 real (kind=RKIND) :: maxvar2d_local, maxvar2d_global
 real (kind=RKIND), dimension(:), pointer :: var2d
 integer, pointer :: nCellsSolve
 character (len=StrKIND), pointer :: gwdo_scheme
 type (block_type), pointer :: block
 type (iamas_pool_type), pointer :: meshPool
 type (iamas_pool_type), pointer :: sfc_inputPool

 ierr = 0

 call iamas_pool_get_config(blocklist % configs, 'config_gwdo_scheme', gwdo_scheme)

 if (trim(gwdo_scheme) /= 'off') then
     maxvar2d_local = -huge(maxvar2d_local)
     block => blockList
     do while (associated(block))
         call iamas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call iamas_pool_get_subpool(block % structs, 'sfc_input', sfc_inputPool)
         call iamas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call iamas_pool_get_array(sfc_inputPool, 'var2d', var2d)

         maxvar2d_local = max(maxvar2d_local, maxval(var2d(1:nCellsSolve)))

         block => block % next
     end do

     call iamas_dmpar_max_real(dminfo, maxvar2d_local, maxvar2d_global)

     if (maxvar2d_global <= 0.0_RKIND) then
         call iamas_log_write('*******************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('The GWDO scheme requires valid var2d, con, oa{1,2,3,4}, and ol{1,2,3,4} fields,', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('but these fields appear to be zero everywhere in the model input.', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Either set config_gwdo_scheme = ''off'' in the &physics namelist, or generate', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('the GWDO static fields with the init_atmosphere core.', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('*******************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         ierr = ierr + 1
     end if

 end if

 end subroutine physics_compatibility_check

!=================================================================================================================
 end module iamas_atmphys_control
!=================================================================================================================
