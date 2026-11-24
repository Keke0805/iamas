! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_packages
 use iamas_kind_types
 use iamas_derived_types,only : iamas_pool_type,iamas_io_context_type,iAMAS_LOG_ERR
 use iamas_pool_routines,only : iamas_pool_get_config,iamas_pool_get_package
 use iamas_log,only : iamas_log_write

 implicit none
 private
 public:: atmphys_setup_packages

!iamas_atmphys_packages contains the definitions of all physics packages.
!Laura D. Fowler (laura@ucar.edu) / 2016-03-10.


 contains


!=================================================================================================================
 function atmphys_setup_packages(configs,packages,iocontext) result(ierr)
!=================================================================================================================
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>>
 logical,pointer:: bl_shinhong_in
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<

!inout arguments:
 type (iamas_pool_type), intent(inout) :: configs
 type (iamas_pool_type), intent(inout) :: packages
 type (iamas_io_context_type), intent(inout) :: iocontext

!local variables:
 character(len=StrKIND),pointer:: config_microp_scheme
 character(len=StrKIND),pointer:: config_mp_p3_scheme
 character(len=StrKIND),pointer:: config_convection_scheme
 character(len=StrKIND),pointer:: config_pbl_scheme
 character(len=StrKIND),pointer:: config_les_scheme
 character(len=StrKIND),pointer:: config_les_surface_mode
 logical,pointer:: mp_kessler_in,mp_thompson_in,mp_wsm6_in,mp_morrison_in,mp_p3_in
 logical,pointer:: mp_p3_1category_in, mp_p3_1category_nc_in, mp_p3_2category_in, mp_p3_1cat_3mom_in
 logical,pointer:: cu_grell_freitas_in,cu_kain_fritsch_in,cu_tiedtke_in
 logical,pointer:: bl_mynn_in,bl_ysu_in
 logical,pointer:: les_sgs_in

 integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('--- enter subroutine atmphys_setup_packages:')

 ierr = 0

 call iamas_log_write('----- Setting up package variables -----')
 call iamas_log_write('')

!--- initialization of all packages for parameterizations of cloud microphysics:

 call iamas_pool_get_config(configs,'config_microp_scheme',config_microp_scheme)

 nullify(mp_kessler_in)
 call iamas_pool_get_package(packages,'mp_kessler_inActive',mp_kessler_in)

 nullify(mp_thompson_in)
 call iamas_pool_get_package(packages,'mp_thompson_inActive',mp_thompson_in)

 nullify(mp_wsm6_in)
 call iamas_pool_get_package(packages,'mp_wsm6_inActive',mp_wsm6_in)

 nullify(mp_morrison_in)
 call iamas_pool_get_package(packages,'mp_morrison_inActive',mp_morrison_in)

 nullify(mp_p3_in)
 call iamas_pool_get_package(packages,'mp_p3_inActive',mp_p3_in)

 if(.not.associated(mp_kessler_in)  .or. &
    .not.associated(mp_thompson_in) .or. &
    .not.associated(mp_wsm6_in)     .or. &
    .not.associated(mp_morrison_in) .or. &
    .not.associated(mp_p3_in)) then
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('* Error while setting up packages for cloud microphysics options in atmosphere core.',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    ierr = 1
    return
 endif

 mp_kessler_in           = .false.
 mp_thompson_in          = .false.
 mp_wsm6_in              = .false.
 mp_morrison_in          = .false.
 mp_p3_in                = .false.

 if(config_microp_scheme == 'mp_kessler') then
    mp_kessler_in = .true.
 elseif(config_microp_scheme == 'mp_thompson') then
    mp_thompson_in = .true.
 elseif(config_microp_scheme == 'mp_wsm6') then
    mp_wsm6_in = .true.
 elseif(config_microp_scheme == 'mp_morrison') then
    mp_morrison_in = .true.
 elseif(config_microp_scheme == 'mp_p3') then
    mp_p3_in = .true.
 endif

 call iamas_log_write('    mp_kessler_in           = $l', logicArgs=(/mp_kessler_in/))
 call iamas_log_write('    mp_thompson_in          = $l', logicArgs=(/mp_thompson_in/))
 call iamas_log_write('    mp_wsm6_in              = $l', logicArgs=(/mp_wsm6_in/))
 call iamas_log_write('    mp_morrison_in          = $l', logicArgs=(/mp_morrison_in/))
 call iamas_log_write('    mp_p3_in                = $l', logicArgs=(/mp_p3_in/))

 if (config_microp_scheme == 'mp_p3') then
    !--- initialization of all packages for parameterizations of p3:
    call iamas_pool_get_config(configs,'config_mp_p3_scheme',config_mp_p3_scheme)

    nullify(mp_p3_1category_in)
    call iamas_pool_get_package(packages,'mp_p3_1category_inActive',mp_p3_1category_in)

    nullify(mp_p3_1category_nc_in)
    call iamas_pool_get_package(packages,'mp_p3_1category_nc_inActive',mp_p3_1category_nc_in)

    nullify(mp_p3_2category_in)
    call iamas_pool_get_package(packages,'mp_p3_2category_inActive',mp_p3_2category_in)

    nullify(mp_p3_1cat_3mom_in)
    call iamas_pool_get_package(packages,'mp_p3_1cat_3mom_inActive',mp_p3_1cat_3mom_in)

    if(.not.associated(mp_p3_1category_in)  .or. &
        .not.associated(mp_p3_1category_nc_in) .or. &
        .not.associated(mp_p3_2category_in)     .or. &
        .not.associated(mp_p3_1cat_3mom_in)) then
        call iamas_log_write('=======================================================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for P3 cloud microphysics options in atmosphere core.',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=======================================================================================',messageType=iAMAS_LOG_ERR)
        ierr = 1
        return
    endif

    mp_p3_1category_in      = .false.
    mp_p3_1category_nc_in   = .false.
    mp_p3_2category_in      = .false.
    mp_p3_1cat_3mom_in      = .false.

    if(config_mp_p3_scheme == 'P3_1CATEGORY') then
        mp_p3_1category_in = .true.
    elseif(config_mp_p3_scheme == 'P3_1CATEGORY_NC') then
        mp_p3_1category_nc_in = .true.
    elseif(config_mp_p3_scheme == 'P3_2CATEGORY') then
        mp_p3_2category_in = .true.
    elseif(config_mp_p3_scheme == 'P3_1CAT_3MOM') then
        mp_p3_1cat_3mom_in = .true.
    endif

    call iamas_log_write('    mp_p3_1category_in           = $l', logicArgs=(/mp_p3_1category_in/))
    call iamas_log_write('    mp_p3_1category_nc_in        = $l', logicArgs=(/mp_p3_1category_nc_in/))
    call iamas_log_write('    mp_p3_2category_in           = $l', logicArgs=(/mp_p3_2category_in/))
    call iamas_log_write('    mp_p3_1cat_3mom_in           = $l', logicArgs=(/mp_p3_1cat_3mom_in/))
  end if

!--- initialization of all packages for parameterizations of convection:

 call iamas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)

 nullify(cu_grell_freitas_in)
 call iamas_pool_get_package(packages,'cu_grell_freitas_inActive',cu_grell_freitas_in)

 nullify(cu_kain_fritsch_in)
 call iamas_pool_get_package(packages,'cu_kain_fritsch_inActive',cu_kain_fritsch_in)

 nullify(cu_tiedtke_in)
 call iamas_pool_get_package(packages,'cu_tiedtke_inActive',cu_tiedtke_in)

 if(.not.associated(cu_grell_freitas_in) .or. &
    .not.associated(cu_kain_fritsch_in)  .or. &
    .not.associated(cu_tiedtke_in)     ) then
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('* Error while setting up packages for convection options in atmosphere core.',        messageType=iAMAS_LOG_ERR)
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    ierr = 1
    return
 endif

 cu_grell_freitas_in = .false.
 cu_kain_fritsch_in  = .false.
 cu_tiedtke_in       = .false.

 if(config_convection_scheme=='cu_grell_freitas') then
    cu_grell_freitas_in = .true.
 elseif(config_convection_scheme == 'cu_kain_fritsch' .or. &
        config_convection_scheme == 'cu_kain_fritsch_sca') then
    cu_kain_fritsch_in = .true.
 elseif(config_convection_scheme == 'cu_tiedtke' .or. &
        config_convection_scheme == 'cu_ntiedtke') then
    cu_tiedtke_in = .true.
 endif

 call iamas_log_write('    cu_grell_freitas_in     = $l', logicArgs=(/cu_grell_freitas_in/))
 call iamas_log_write('    cu_kain_fritsch_in      = $l', logicArgs=(/cu_kain_fritsch_in/))
 call iamas_log_write('    cu_tiedtke_in           = $l', logicArgs=(/cu_tiedtke_in/))

!--- initialization of all packages for parameterizations of surface layer and planetary boundary layer:

 call iamas_pool_get_config(configs,'config_pbl_scheme',config_pbl_scheme)
 call iamas_pool_get_config(configs,'config_les_scheme',config_les_scheme)
 call iamas_pool_get_config(configs,'config_les_surface_mode',config_les_surface_mode)

 nullify(bl_mynn_in)
 call iamas_pool_get_package(packages,'bl_mynn_inActive',bl_mynn_in)

 nullify(bl_ysu_in)
 call iamas_pool_get_package(packages,'bl_ysu_inActive',bl_ysu_in)
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
 nullify(bl_shinhong_in)
 call iamas_pool_get_package(packages,'bl_shinhong_inActive',bl_shinhong_in)
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<

 if(.not.associated(bl_mynn_in) .or. &
    .not.associated(bl_ysu_in) .or. &
    .not.associated(bl_shinhong_in)) then!qkyang added  bl_shinhong_in to support Shin-Hong PBL
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('* Error while setting up packages for planetary layer  options in atmosphere core.',  messageType=iAMAS_LOG_ERR)
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    ierr = 1
    return
 endif

 bl_mynn_in = .false.
 bl_ysu_in  = .false.
 !qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
 bl_shinhong_in  = .false.
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<
 if(config_pbl_scheme=='bl_mynn') then
    bl_mynn_in = .true.
 elseif(config_pbl_scheme == 'bl_ysu') then
    bl_ysu_in = .true.
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
 elseif(config_pbl_scheme == 'bl_shinhong') then
    bl_shinhong_in = .true.
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<
 endif

 ! The native surface path uses fields carried by the existing YSU allocation
 ! package even when the YSU PBL driver is off.  Preserve the established LES
 ! physics_flux allocation.  The package-init config pool does not carry the
 ! surface-layer selector, so allocate the same fields for LES=off; downstream
 ! surface and PBL execution remain controlled by their original config guards.
 if (trim(config_pbl_scheme) == 'off' .and. &
     ((trim(config_les_scheme) == 'les_deardorff' .and. &
       trim(config_les_surface_mode) == 'physics_flux') .or. &
      trim(config_les_scheme) == 'off')) then
    bl_ysu_in = .true.
 endif

 call iamas_log_write('    bl_mynn_in              = $l', logicArgs=(/bl_mynn_in/))
 call iamas_log_write('    bl_ysu_in               = $l', logicArgs=(/bl_ysu_in/))
 call iamas_log_write('    bl_shinhong_in          = $l', logicArgs=(/bl_shinhong_in/))!qkyang added to support Shin-Hong PBL

!--- initialization of LES SGS package:

 call iamas_pool_get_config(configs,'config_les_scheme',config_les_scheme)

 nullify(les_sgs_in)
 call iamas_pool_get_package(packages,'les_sgs_inActive',les_sgs_in)

 if(.not.associated(les_sgs_in)) then
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('* Error while setting up packages for LES SGS options in atmosphere core.',messageType=iAMAS_LOG_ERR)
    call iamas_log_write('====================================================================================',messageType=iAMAS_LOG_ERR)
    ierr = 1
    return
 endif

 les_sgs_in = .false.
 if(trim(config_les_scheme) == 'les_deardorff') then
    les_sgs_in = .true.
 endif

 call iamas_log_write('    les_sgs_in              = $l', logicArgs=(/les_sgs_in/))
 call iamas_log_write('')

 end function atmphys_setup_packages

!=================================================================================================================
 end module iamas_atmphys_packages
!=================================================================================================================


