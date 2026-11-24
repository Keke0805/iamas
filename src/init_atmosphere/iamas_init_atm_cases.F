! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module init_atm_cases

   use iamas_kind_types
   use iamas_derived_types
   use iamas_pool_routines
   use iamas_constants
   use iamas_dmpar
   use iamas_sort
   use atm_advection
   use iamas_RBF_interpolation
   use iamas_vector_reconstruction
   use iamas_timer
   use iamas_init_atm_static
   use iamas_init_atm_surface
   use iamas_init_atm_fdda ! jwfeng add this for fdda nudging. 1217/2025
   use iamas_init_atm_emission ! zhxia add this for emission. 0414/2025
   use iamas_init_atm_wrfchem ! zhxia add this for reading chemical field from wrfchem. 06/2025
   use iamas_atmphys_constants, only: svpt0,svp1,svp2,svp3
   use iamas_atmphys_functions
   use iamas_atmphys_initialize_real
   use iamas_log, only : iamas_log_write
   !qkyang added to support building turbulence model>>>>>>>>>>>>
   use AddUpStreamBuildParam_mod, only : AddUpStreamBuildParam_fromPools
   !qkyang added to support building turbulence model<<<<<<<<<<<
   !qkyang added to support slope and shading effects>>>>>>>>>>>>>>>>>
   use AddSlopeShadowParam_mod, only : AddSlopeShadowParam_fromPools
   !qkyang added to support slope and shading effects<<<<<<<<<<<<<<<<

   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in iamas_timekeeping
   use iamas_timekeeping !, only: iAMAS_Time_type, iAMAS_TimeInterval_type, iAMAS_Clock_type, &
                        !        iamas_set_time, iamas_set_timeInterval, iamas_get_time, operator(+), add_t_ti



   contains
!MODIFIED BY jwfeng
   subroutine init_chemistry(mesh, state)

      implicit none

      type (iamas_pool_type), intent(in) :: mesh
      type (iamas_pool_type), intent(inout) :: state

      integer, pointer :: nCells
      integer, pointer :: nVertLevels

      integer :: k, iCell

      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer :: lonCell

      real (kind=RKIND), dimension(:,:,:), pointer :: chem
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call iamas_pool_get_dimension(mesh,'nCells',nCells)
      call iamas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
      call iamas_pool_get_array(state,'chem', chem)
      call iamas_pool_get_array(mesh,'latCell',latCell)
      call iamas_pool_get_array(mesh,'lonCell',lonCell)

      chem(:,:,:)= 1.0e-16
! jwfeng from AEMOL USTC:
! These codes can init specific chem in certain area:
!LAT:35-42N,LON:75-90E,nVert=0
!     do iCell = 1, nCells
!         if (lonCell(iCell)*57.2957795 > 75.0 .and. &
!             lonCell(iCell)*57.2957795 < 90.0 .and. &
!             latCell(iCell)*57.2957795 > 35.0 .and. &
!             latCell(iCell)*57.2957795 < 42.0) then
!                 chem(index_dust_a01,1,iCell) = 0.001
!                 chem(index_dust_a02,1,iCell) = 0.001
!         end if
!     end do

   end subroutine init_chemistry
!END MODIFIED BY jwfeng

   subroutine init_atm_setup_case(domain, stream_manager)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the hydrostatic test case
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use iamas_stream_manager
      use iamas_init_atm_gwd, only : compute_gwd_fields

      implicit none

      type (domain_type), intent(inout) :: domain
      type (iAMAS_streamManager_type), intent(inout) :: stream_manager


      integer :: i
      integer :: ierr
      type (block_type), pointer :: block_ptr

      type (iamas_pool_type), pointer :: mesh
      type (iamas_pool_type), pointer :: fg
      type (iamas_pool_type), pointer :: emis_state
      type (iamas_pool_type), pointer :: state
      type (iamas_pool_type), pointer :: diag
      type (iamas_pool_type), pointer :: diag_physics
      type (iamas_pool_type), pointer :: lbc_state

      integer, pointer :: config_init_case
      logical, pointer :: config_static_interp
      logical, pointer :: config_native_gwd_static
      logical, pointer :: config_met_interp
      logical, pointer :: config_blend_bdy_terrain
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_met_prefix
      logical, pointer :: config_input_chem
      logical, pointer :: config_use_fix_gasinit !determine whether we use the photochemistry simulation results from wrfchem as initial values in iAMAS.
      logical, pointer :: config_use_fix_ozone_init !determine whether we use the photochemistry simulation results from wrfchem as initial values in iAMAS.

      character(len=StrKIND), pointer :: mminlu
      character(len=StrKIND), pointer :: xtime

      type (iAMAS_Time_type)  :: curr_time, stop_time
      type (iAMAS_TimeInterval_type)  :: clock_interval, lbc_stream_interval
      character(len=StrKIND) :: timeString

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: nVertLevels

      ! The next four variables are needed in the argument list for blend_bdy_terrain
      ! with the dryrun argument set to true; accordingly, we never actually need to
      ! set these pointers to fields
      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer :: lonCell
      real (kind=RKIND), dimension(:), pointer :: ter
      integer, dimension(:), pointer :: bdyMaskCell


      call iamas_pool_get_config(domain % blocklist % configs, 'config_init_case', config_init_case)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_input_chem', config_input_chem)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_use_fix_gasinit', config_use_fix_gasinit)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_use_fix_ozone_init', config_use_fix_ozone_init)


      if ((config_init_case == 1) .or. (config_init_case == 2) .or. (config_init_case == 3)) then

         call iamas_log_write(' Jablonowski and Williamson baroclinic wave test case ')
         if (config_init_case == 1) call iamas_log_write(' no initial perturbation ')
         if (config_init_case == 2) call iamas_log_write(' initial perturbation included ')
         if (config_init_case == 3) call iamas_log_write(' normal-mode perturbation included ')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call iamas_pool_get_config(block_ptr % configs, 'config_static_interp', config_static_interp)
            call iamas_pool_get_config(block_ptr % configs, 'config_met_interp', config_met_interp)

            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            call iamas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call iamas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            call iamas_log_write(' calling test case setup ')
            call init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, block_ptr % configs, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            call iamas_log_write(' returned from test case setup ')
            block_ptr => block_ptr % next
         end do

      else if ((config_init_case == 4) .or. (config_init_case == 5)) then

         call iamas_log_write(' squall line - super cell test case ')
         if (config_init_case == 4) call iamas_log_write(' squall line test case')
         if (config_init_case == 5) call iamas_log_write(' supercell test case')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call iamas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call iamas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            call iamas_log_write(' calling test case setup ')
            call init_atm_case_squall_line(domain % dminfo, mesh, nCells, nVertLevels, state, diag, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            call iamas_log_write(' returned from test case setup ')
            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 6 ) then

         call iamas_log_write(' mountain wave test case ')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call iamas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call iamas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            call iamas_log_write(' calling test case setup ')
            call init_atm_case_mtn_wave(domain % dminfo, mesh, nCells, nVertLevels, state, diag, block_ptr % configs, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            call iamas_log_write(' returned from test case setup ')
            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 7 ) then


         call iamas_log_write(' real-data GFS test case ')
         block_ptr => domain % blocklist

         do while (associated(block_ptr))

            call iamas_pool_get_config(block_ptr % configs, 'config_static_interp', config_static_interp)
            call iamas_pool_get_config(block_ptr % configs, 'config_native_gwd_static', config_native_gwd_static)
            call iamas_pool_get_config(block_ptr % configs, 'config_met_interp', config_met_interp)
            call iamas_pool_get_config(block_ptr % configs, 'config_blend_bdy_terrain', config_blend_bdy_terrain)

            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'fg', fg)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)

            call iamas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call iamas_pool_get_dimension(block_ptr % dimensions, 'nEdges', nEdges)
            call iamas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            !
            ! Before proceeding with any other processing that takes non-trivial time (e.g., static field interpolation),
            ! check that the intermediate file with terrain information exists if config_blend_bdy_terrain = true.
            !
            ! NB: When calling blend_bdy_terrain(...) with the 'dryrun' argument set, the nCells, latCell, lonCell,
            !     bdyMaskCell, and ter arguments are not used -- only the config_met_prefix and config_start_time
            !     arguments are used.
            !
            if (config_blend_bdy_terrain) then
                call iamas_pool_get_config(block_ptr % configs, 'config_start_time', config_start_time)
                call iamas_pool_get_config(block_ptr % configs, 'config_met_prefix', config_met_prefix)

                call blend_bdy_terrain(config_met_prefix, config_start_time, &
                                       nCells, latCell, lonCell, bdyMaskCell, ter, .true., ierr)
                if (ierr /= 0) then
                    call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
                    call iamas_log_write('Blending of terrain along domain boundaries would fail, and', messageType=iAMAS_LOG_ERR)
                    call iamas_log_write('config_blend_bdy_terrain = true in the namelist.init_atmosphere file.', &
                                        messageType=iAMAS_LOG_ERR)
                    call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
                end if
            end if

            if (config_static_interp) then

               !
               !  Without a convex mesh partition file, interpolating static fields in parallel
               !     will give incorrect results. Since it is very unlikely that typical users
               !     will have convex partitions, it's safer to just stop if multiple MPI tasks are
               !     detected when performing the static_interp step.
               !  
               if (domain % dminfo % nprocs > 1) then !czhao +++ 
!                call iamas_dmpar_global_abort('****************************************************************', deferredAbort=.true.)
!                call iamas_dmpar_global_abort('Error: Interpolation of static fields does not work in parallel.', deferredAbort=.true.)
!                call iamas_dmpar_global_abort('Please run the static_interp step using only a single MPI task.',  deferredAbort=.true.)
!                call iamas_dmpar_global_abort('****************************************************************')
               end if  !czhao -----
               call init_atm_static(mesh, block_ptr % dimensions, block_ptr % configs)
            end if

            if (config_native_gwd_static) then
               call iamas_log_write(' ')
               call iamas_log_write('Computing GWDO static fields on the native iAMAS mesh')
               call iamas_log_write(' ')
               ierr = compute_gwd_fields(domain)
               if (ierr /= 0) then
                  call iamas_log_write('****************************************************************', messageType=iAMAS_LOG_ERR)
                  call iamas_log_write('Error while trying to compute sub-grid-scale orography',           messageType=iAMAS_LOG_ERR)
                  call iamas_log_write('statistics for use with the GWDO scheme.',                         messageType=iAMAS_LOG_ERR)
                  call iamas_log_write('****************************************************************', messageType=iAMAS_LOG_CRIT)
               end if
            end if

            !
            ! If at this point the mminlu variable is blank, we assume that the static interp step was
            !   not run, and that we are working with a static file created before there was a choice
            !   of land use datasets; in this case, the dataset was almost necessarily USGS
            !
            call iamas_pool_get_array(mesh, 'mminlu', mminlu)
            if (len_trim(mminlu) == 0) then
                  call iamas_log_write('****************************************************************')
                  call iamas_log_write('No information on land use dataset is available.')
                  call iamas_log_write('Assume that we are using ''USGS''.')
                  call iamas_log_write('****************************************************************')
                  write(mminlu,'(a)') 'USGS'
            end if

            call init_atm_case_gfs(block_ptr, mesh, nCells, nEdges, nVertLevels, fg, state, &
                 diag, diag_physics, config_init_case, block_ptr % dimensions, block_ptr % configs)

            ! jwfeng from AEMOL USTC:
            if (.not. config_input_chem) then
               call iamas_log_write('****************************************************************')
               call iamas_log_write('Initialize chem to default value (1.0e-16 ug/kg).'               )
               call iamas_log_write('****************************************************************')
               call init_chemistry(mesh,state)

               !Add this to use the chemistry simulation results from WRF-Chem as initial values in iamas. zhxia@lacar. 06/2025
               if (config_use_fix_gasinit) then
                   call init_atm_case_wrfchem(block_ptr, mesh, nCells, nEdges, nVertLevels, fg, state, &
                                              diag, diag_physics, block_ptr % dimensions, block_ptr % configs)
               end if

               !Add this to use the ozone data from TOST as initial values in iamas. zhxia@lacar. 10/2025
               if (config_use_fix_ozone_init) then
                   call init_atm_case_ozone(block_ptr, mesh, nCells, nEdges, nVertLevels, state, &
                                              diag, diag_physics, block_ptr % dimensions, block_ptr % configs)
               end if
            end if

            if (config_met_interp) then
               call physics_initialize_real(mesh, fg, domain % dminfo, block_ptr % dimensions, block_ptr % configs)
            end if

            !qkyang added to support building turbulence model>>>>>>>>>>>>
            call AddUpStreamBuildParam_fromPools(mesh, block_ptr % dimensions, block_ptr % configs)
            !qkyang added to support building turbulence model<<<<<<<<<<<

            !qkyang added to support slope and shading effects>>>>>>>>>>>>>>>>>
            call AddSlopeShadowParam_fromPools(mesh, block_ptr % dimensions, block_ptr % configs)
            !qkyang added to support slope and shading effects<<<<<<<<<<<<<<<<

            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 8 ) then

         call iamas_log_write('real-data surface (SST) update test case ')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'fg', fg)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)

            ! Defined in iamas_init_atm_surface.F
            call init_atm_case_sfc(domain, domain % dminfo, stream_manager, mesh, fg, state, block_ptr % dimensions, block_ptr % configs)
            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 9 ) then

         call iamas_log_write('Lateral boundary conditions case')

         !
         ! Check that the first-guess interval (which is the same as the clock timestep)
         ! matches the output interval of the 'lbc' stream
         !
         clock_interval = iamas_get_clock_timestep(domain % clock, ierr=ierr)
         lbc_stream_interval = iAMAS_stream_mgr_get_stream_interval(stream_manager, 'lbc', iAMAS_STREAM_OUTPUT, ierr)
         if (clock_interval /= lbc_stream_interval) then
            call iamas_log_write('****************************************************************',   messageType=iAMAS_LOG_ERR)
            call iamas_log_write('The intermediate file interval specified by ''config_fg_interval''', messageType=iAMAS_LOG_ERR)
            call iamas_log_write('does not match the output_interval for the ''lbc'' stream.',         messageType=iAMAS_LOG_ERR)
            call iamas_log_write('Please correct the namelist.init_atmosphere and/or',                 messageType=iAMAS_LOG_ERR)
            call iamas_log_write('streams.init_atmosphere files.',                                     messageType=iAMAS_LOG_ERR)
            call iamas_log_write('****************************************************************',   messageType=iAMAS_LOG_CRIT)
         end if

         curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW)
         stop_time = iamas_get_clock_time(domain % clock, iAMAS_STOP_TIME)

         do while (curr_time <= stop_time)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
               call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
               call iamas_pool_get_subpool(block_ptr % structs, 'fg', fg)
               call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
               call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call iamas_pool_get_subpool(block_ptr % structs, 'lbc_state', lbc_state)

               call iamas_pool_get_array(state, 'xtime', xtime)

               call iamas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
               call iamas_pool_get_dimension(block_ptr % dimensions, 'nEdges', nEdges)
               call iamas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

               call iamas_get_time(curr_time, dateTimeString=timeString)
               xtime = timeString   ! Set field valid time, xtime, to the current time in the time loop

               call init_atm_case_lbc(timeString, block_ptr, mesh, nCells, nEdges, nVertLevels, fg, state, &
                                      diag, lbc_state, block_ptr % dimensions, block_ptr % configs)

               block_ptr => block_ptr % next
            end do

            call iamas_stream_mgr_write(stream_manager, streamID='lbc', ierr=ierr)
            call iamas_stream_mgr_reset_alarms(stream_manager, streamID='lbc', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

            call iamas_advance_clock(domain % clock)
            curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW)

         end do

         !
         ! Ensure that no output alarms are still ringing for the 'lbc' stream after
         ! we exit the time loop above; the main run routine may write out all other
         ! output streams with ringing alarms.
         !
         call iamas_stream_mgr_reset_alarms(stream_manager, streamID='lbc', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

      !zhxia add this for emission. 0414/2025
      else if (config_init_case == 10 ) then

         call iamas_log_write('real-data emission update test case ')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'emis_state', emis_state)
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)

            ! Defined in iamas_init_atm_emission.F
            call init_atm_case_emis(domain, domain % dminfo, stream_manager, mesh, emis_state, state, block_ptr % dimensions, block_ptr % configs)
            block_ptr => block_ptr % next
         end do

      !jwfeng add this for fdda nudging. 1217/2025
      else if (config_init_case == 11 ) then

         call iamas_log_write('real-data FDDA update case ')
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            ! Defined in iamas_init_atm_fdda.F
            call init_atm_case_fdda(block_ptr, domain, domain % dminfo, stream_manager)
            block_ptr => block_ptr % next

         end do

      else

         call iamas_log_write(' ****************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(' Only test cases 1 through 11 are currently supported.', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(' ****************************************************', messageType=iAMAS_LOG_CRIT)

      end if


      !initialization of surface input variables technically not needed to run our current set of
      !idealized test cases:
      if (config_init_case < 7)  then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call iamas_pool_get_subpool(block_ptr % structs, 'fg', fg)

            call physics_idealized_init(mesh, fg)

            block_ptr => block_ptr % next
         end do
      end if

   end subroutine init_atm_setup_case

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, configs, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(in) :: configs
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

!      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:), pointer :: u, ru, w, rw, v
      real (kind=RKIND), dimension(:,:), pointer :: rho, theta
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

!.. initialization of moisture:
      integer, pointer :: index_qv
      real (kind=RKIND), parameter :: rh_max = 0.40 ! Maximum relative humidity
!      real (kind=RKIND), parameter :: rh_max = 0.70 ! Maximum relative humidity
      real (kind=RKIND), dimension(nVertLevels, nCells) :: qsat, relhum
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
!.. end initialization of moisture.

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
      integer, pointer :: nz1, nCellsSolve, nEdges, maxEdges, nVertices

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, verticesOnEdge, cellsOnCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND) :: es, qvs, xnutr, znut, ptemp
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(nVertLevels ) :: eta, etav, teta, ppi, tt, temperature_1d

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf
      real (kind=RKIND), pointer :: cf1, cf2, cf3

      !  storage for (lat,z) arrays for zonal velocity calculation

      logical, parameter :: rebalance = .true.
      integer, parameter :: nlat=721
      real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nVertLevels + 1, nlat) :: zgrid_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: u_2d, pp_2d, rho_2d, qv_2d, etavs_2d, zz_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: p_2d, pb_2d, ppb_2d, rr_2d, rb_2d, tb_2d, rtb_2d
      real (kind=RKIND), dimension(nVertLevels, nlat-1) :: zx_2d
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat, hx_1d
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

!      logical, parameter :: moisture = .true.
      logical, parameter :: moisture = .false.

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:), pointer :: latCell, latVertex, lonVertex, latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: fEdge, fVertex

      real (kind=RKIND), pointer :: sphere_radius
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: config_init_case


      call iamas_pool_get_config(configs, 'config_init_case', config_init_case)
      call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call iamas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)


      !
      ! Scale all distances and areas from a unit sphere to one with radius sphere_radius
      !
      call iamas_pool_get_array(mesh, 'xCell', xCell)
      call iamas_pool_get_array(mesh, 'yCell', yCell)
      call iamas_pool_get_array(mesh, 'zCell', zCell)
      call iamas_pool_get_array(mesh, 'xEdge', xEdge)
      call iamas_pool_get_array(mesh, 'yEdge', yEdge)
      call iamas_pool_get_array(mesh, 'zEdge', zEdge)
      call iamas_pool_get_array(mesh, 'xVertex', xVertex)
      call iamas_pool_get_array(mesh, 'yVertex', yVertex)
      call iamas_pool_get_array(mesh, 'zVertex', zVertex)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call iamas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call iamas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call iamas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      xCell(:) = xCell(:) * sphere_radius
      yCell(:) = yCell(:) * sphere_radius
      zCell(:) = zCell(:) * sphere_radius
      xVertex(:) = xVertex(:) * sphere_radius
      yVertex(:) = yVertex(:) * sphere_radius
      zVertex(:) = zVertex(:) * sphere_radius
      xEdge(:) = xEdge(:) * sphere_radius
      yEdge(:) = yEdge(:) * sphere_radius
      zEdge(:) = zEdge(:) * sphere_radius
      dvEdge(:) = dvEdge(:) * sphere_radius
      dcEdge(:) = dcEdge(:) * sphere_radius
      areaCell(:) = areaCell(:) * sphere_radius**2.0
      areaTriangle(:) = areaTriangle(:) * sphere_radius**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * sphere_radius**2.0

      call iamas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call iamas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call iamas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)

      call iamas_pool_get_dimension(mesh, 'nVertLevels', nz1)
      call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz = nz1 + 1

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_array(mesh, 'dzu', dzu)
      call iamas_pool_get_array(mesh, 'rdzu', rdzu)
      call iamas_pool_get_array(mesh, 'fzm', fzm)
      call iamas_pool_get_array(mesh, 'fzp', fzp)
      call iamas_pool_get_array(mesh, 'zxu', zxu)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(mesh, 'hx', hx)
      call iamas_pool_get_array(mesh, 'dss', dss)

      call iamas_pool_get_array(diag, 'exner_base', pb)
      call iamas_pool_get_array(diag, 'rho_base', rb)
      call iamas_pool_get_array(diag, 'rho_p', rr)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag, 'theta_base', tb)
      call iamas_pool_get_array(diag, 'theta', theta)
      call iamas_pool_get_array(diag, 'rtheta_base', rtb)
      call iamas_pool_get_array(diag, 'rtheta_p', rt)
      call iamas_pool_get_array(diag, 'exner', p)
      call iamas_pool_get_array(diag, 'pressure_base', ppb)
      call iamas_pool_get_array(diag, 'pressure_p', pp)
      call iamas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call iamas_pool_get_array(diag, 'ru', ru)
      call iamas_pool_get_array(diag, 'rw', rw)
      call iamas_pool_get_array(diag, 'v', v)

      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(state, 'theta_m', t)
      call iamas_pool_get_array(state, 'scalars', scalars)
      call iamas_pool_get_array(state, 'u', u)
      call iamas_pool_get_array(state, 'w', w)

      call iamas_pool_get_array(mesh, 'latCell', latCell)
      call iamas_pool_get_array(mesh, 'latVertex', latVertex)
      call iamas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call iamas_pool_get_array(mesh, 'latEdge', latEdge)
      call iamas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call iamas_pool_get_array(mesh, 'fEdge', fEdge)
      call iamas_pool_get_array(mesh, 'fVertex', fVertex)

      call iamas_pool_get_array(mesh, 'cf1', cf1)
      call iamas_pool_get_array(mesh, 'cf2', cf2)
      call iamas_pool_get_array(mesh, 'cf3', cf3)

!.. initialization of moisture:
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0
!.. end initialization of moisture.

      surface_pressure(:) = 0.0

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      xnutr = 0.0
      zd = 12000.0
      znut = eta_t

      etavs = (1.0 - 0.252) * pii/2.
      r_earth = sphere_radius
      omega_e = omega
      p0 = 1.0e+05

      call iamas_log_write(' point 1 in test case setup ')

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,nCells
        do k=1,nz
          phi = latCell(iCell)
          hx(k,iCell) = u0/gravity*cos(etavs)**1.5                                   &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        end do
      end do

      !     Metrics for hybrid coordinate and vertical stretching

      str = 1.5
      zt = 45000.
      dz = zt/float(nz1)

      call iamas_log_write(' hx computation complete ')

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            sh(k) = (real(k-1)*dz/zt)**str

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) and define sh(k) = zc(k)/zt
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz

            zw(k) = float(k-1)*dz
!            zw(k) = sh(k)*zt
!
!           ah(k) governs the transition between terrain-following
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate

            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
!            ah(k) = 0.
            call iamas_log_write(' k, sh, zw, ah = $i $r $r $r', intArgs=(/k/), realArgs=(/sh(k),zw(k),ah(k)/))
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2)
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3)
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      call iamas_log_write(' cf1, cf2, cf3 = $r $r $r', realArgs=(/cf1,cf2,cf3/))

      do iCell=1,nCells
        do k=1,nz
          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
                         + ah(k) * sh(k)* zt
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

      !do k=1,nz1
      !  call iamas_log_write(' k, zgrid(k,1),hx(k,1) = $i $r $r', intArgs=(/k/), realArgs=(/zgrid(k,1),hx(k,1)/))
      !end do

      !do k=1,nz1
      !  call iamas_log_write(' k, zxu(k,1) = $i $r', intArgs=(/k/), realArgs=(/zxu(k,1)/))
      !end do

      call iamas_log_write(' grid metrics setup complete ')

!**************  section for 2d (z,lat) calc for zonal velocity

      dlat = 0.5*pii/float(nlat-1)
      do i = 1,nlat

        lat_2d(i) = float(i-1)*dlat
        phi = lat_2d(i)
        hx_1d    = u0/gravity*cos(etavs)**1.5                           &
                   *((-2.*sin(phi)**6                                   &
                         *(cos(phi)**2+1./3.)+10./63.)                  &
                         *(u0)*cos(etavs)**1.5                          &
                    +(1.6*cos(phi)**3                                   &
                         *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)

        do k=1,nz
          zgrid_2d(k,i) = (1.-ah(k))*(sh(k)*(zt-hx_1d)+hx_1d)  &
                         + ah(k) * sh(k)* zt
        end do
        do k=1,nz1
          zz_2d(k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
          ppb_2d(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b))
          pb_2d(k,i) = (ppb_2d(k,i)/p0)**(rgas/cp)
          rb_2d(k,i) = ppb_2d(k,i)/(rgas*t0b*zz_2d(k,i))
          tb_2d(k,i) = t0b/pb_2d(k,i)
          rtb_2d(k,i) = rb_2d(k,i)*tb_2d(k,i)
          p_2d(k,i) = pb_2d(k,i)
          pp_2d(k,i) = 0.0
          rr_2d(k,i) = 0.0
        end do


        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb_2d(k,i)+pp_2d(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k) >= znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do

          phi = lat_2d(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))


            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
            ptemp   = ppb_2d(k,i) + pp_2d(k,i)

            !get moisture
            if (moisture) then
              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
          end do

          do itrp = 1,25
            do k=1,nz1
              rr_2d(k,i)  = (pp_2d(k,i)/(rgas*zz_2d(k,i)) - rb_2d(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                            &
                          *(1.25*(rr_2d(1,i)+rb_2d(1,i))*(1.+qv_2d(1,i))   &
                            -.25*(rr_2d(2,i)+rb_2d(2,i))*(1.+qv_2d(2,i)))

            ppi(1) = ppi(1)-ppb_2d(1,i)
            do k=1,nz1-1

              ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                       &
                            ( (rr_2d(k  ,i)+(rr_2d(k  ,i)+rb_2d(k  ,i))*qv_2d(k  ,i))*fzp(k+1)   &
                            + (rr_2d(k+1,i)+(rr_2d(k+1,i)+rb_2d(k+1,i))*qv_2d(k+1,i))*fzm(k+1) )
            end do

            do k=1,nz1
              pp_2d(k,i) = .2*ppi(k)+.8*pp_2d(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          rho_2d(k,i) = rr_2d(k,i)+rb_2d(k,i)
          etavs_2d(k,i) = ((ppb_2d(k,i)+pp_2d(k,i))/p0 - 0.252)*pii/2.
          u_2d(k,i) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(k,i))**1.5)
        end do

      end do  ! end loop over latitudes for 2D zonal wind field calc

      !SHP-balance:: in case of rebalacing for geostrophic wind component
      if (rebalance) then

        do i=1,nlat-1
          do k=1,nz1
            zx_2d(k,i) = (zgrid_2d(k,i+1)-zgrid_2d(k,i))/(dlat*r_earth)
          end do
        end do

        call init_atm_recompute_geostrophic_wind(u_2d, rho_2d, pp_2d, qv_2d, lat_2d, zz_2d, zx_2d,     &
                                        cf1, cf2, cf3, fzm, fzp, rdzw, nz1, nlat, dlat, sphere_radius)

      end if

!******************************************************************

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, nCells
        do k=1,nz1
          ztemp    = .5*(zgrid(k+1,i)+zgrid(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b))
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do

!       if(i == 1) then
!         do k=1,nz1
!           call iamas_log_write(' k, ppb, pb, rb, tb (k,1) = $i $r $r $r $r', intArgs=(/k/), realArgs=(/ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)/))
!         end do
!       end if

!     iterations to converge temperature as a function of pressure
!
        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb(k,i)+pp(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k).ge.znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do
          phi = latCell(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*scalars(index_qv,k,i))

            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture
            if (moisture) then

                !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )

               if(ptemp < 50000.) then
                  relhum(k,i) = 0.0
               elseif(ptemp > p0) then
                  relhum(k,i) = 1.0
               else
                  relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
               end if
               relhum(k,i) = min(rh_max,relhum(k,i))

               !.. calculation of water vapor mixing ratio:
               if (temperature_1d(k) > 273.15) then
                   es  = 1000.*0.6112*exp(17.67*(temperature_1d(k)-273.15)/(temperature_1d(k)-29.65))
               else
                   es  = 1000.*0.6112*exp(21.8745584*(temperature_1d(k)-273.15)/(temperature_1d(k)-7.66))
               end if
               qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
               if(relhum(k,i) .eq. 0.0) qsat(k,i) = 0.0
               scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*scalars(index_qv,k,i))

          end do

          do itrp = 1,25
            do k=1,nz1
              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1

!              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)   &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))

               ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                                  &
                             ( (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))*fzp(k+1)   &
                             + (rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))*fzm(k+1) )

            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          p (k,i) = ((ppb(k,i)+pp(k,i))/p0)**(rgas/cp)
          t (k,i) = tt(k)/p(k,i)
          rt (k,i) = t(k,i)*rr(k,i)+rb(k,i)*(t(k,i)-tb(k,i))
          rho_zz (k,i) = rb(k,i) + rr(k,i)
        end do

        !calculation of surface pressure:
        surface_pressure(i) = 0.5*dzw(1)*gravity                                    &
                        * (1.25*(rr(1,i) + rb(1,i)) * (1. + scalars(index_qv,1,i))  &
                        -  0.25*(rr(2,i) + rb(2,i)) * (1. + scalars(index_qv,2,i)))
        surface_pressure(i) = surface_pressure(i) + pp(1,i) + ppb(1,i)

      end do  ! end loop over cells

      !call iamas_log_write('')
      !call iamas_log_write('--- initialization of water vapor:')
      !do iCell = 1, nCells
      !   if(iCell == 1 .or. iCell == nCells) then
      !      do k = nz1, 1, -1
      !         call iamas_log_write('$i $i $r $r $r $r', intArgs=(/iCell,k/), realArgs=(/t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)/))
      !      end do
       !     call iamas_log_write('')
      !   end if
      !end do

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,nEdges

         vtx1 = verticesOnEdge(1,iEdge)
         vtx2 = verticesOnEdge(2,iEdge)
         lat1 = latVertex(vtx1)
         lat2 = latVertex(vtx2)
         iCell1 = cellsOnEdge(1,iEdge)
         iCell2 = cellsOnEdge(2,iEdge)
         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)

         if (config_init_case == 2) then
            r_pert = sphere_distance( latEdge(iEdge), lonEdge(iEdge), &
                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1) * sphere_radius / dvEdge(iEdge)

         else if (config_init_case == 3) then
            lon_Edge = lonEdge(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
         else
            u_pert = 0.0
         end if

         if (rebalance) then

           call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,dvEdge(iEdge),sphere_radius,u0,nz1,nlat)
           do k=1,nVertLevels
             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
             u(k,iEdge) = fluxk + u_pert
           end do

         else

           do k=1,nVertLevels
             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
             fluxk = u0*flux*(cos(etavs)**1.5)
             u(k,iEdge) = fluxk + u_pert
           end do

         end if

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nz1
            ru(k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,iEdge)
         end do

      !
      ! Generate rotated Coriolis field
      !

         fEdge(iEdge) = 2.0 * omega_e * &
                                       ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha_grid) + &
                                         sin(latEdge(iEdge)) * cos(alpha_grid) &
                                       )
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 2.0 * omega_e * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha_grid) + &
                                          sin(latVertex(iVtx)) * cos(alpha_grid) &
                                         )
      end do

      !
      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
      !

      !
      !     pre-calculation z-metric terms in omega eqn.
      !
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         
         ! Avoid a potential divide by zero below if areaCell(nCells+1) is used in the denominator
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
            do k = 1, nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)

!  WCS fix 20120711

                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do

                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  else
                     z_edge3 = 0.
                  end if

               end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

            end do
         end if

      end do

      ! for including terrain
      rw = 0.0
      w = 0.0
      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then
               rw(k,cell2) = rw(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do


      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
           end do
         end do
      end do

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

        call iamas_log_write(' i, psurf, lat = $i $r $r', intArgs=(/i/), realArgs=(/psurf,latCell(i)*180./3.1415828/))
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_jw


   subroutine init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1_in,lat2_in,dvEdge,a,u0,nz1,nlat)

      implicit none

      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: u_2d,etavs_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(out) :: flux_zonal
      real (kind=RKIND), intent(in) :: lat1_in, lat2_in, dvEdge, a, u0

      integer :: k,i
      real (kind=RKIND) :: lat1, lat2, w1, w2
      real (kind=RKIND) :: dlat,da,db

      lat1 = abs(lat1_in)
      lat2 = abs(lat2_in)
      if(lat2 <= lat1) then
        lat1 = abs(lat2_in)
        lat2 = abs(lat1_in)
      end if

      do k=1,nz1
        flux_zonal(k) = 0.
      end do

      do i=1,nlat-1
        if( (lat1 <= lat_2d(i+1)) .and. (lat2 >= lat_2d(i)) ) then

        dlat = lat_2d(i+1)-lat_2d(i)
        da = (max(lat1,lat_2d(i))-lat_2d(i))/dlat
        db = (min(lat2,lat_2d(i+1))-lat_2d(i))/dlat
        w1 = (db-da) -0.5*(db-da)**2
        w2 = 0.5*(db-da)**2

        do k=1,nz1
          flux_zonal(k) = flux_zonal(k) + w1*u_2d(k,i) + w2*u_2d(k,i+1)
        end do

        end if

      end do

   !  renormalize for setting cell-face fluxes

      do k=1,nz1
        flux_zonal(k) = sign(1.0_RKIND,lat2_in-lat1_in)*flux_zonal(k)*dlat*a/dvEdge/u0
      end do

   end subroutine init_atm_calc_flux_zonal


   !SHP-balance
   subroutine init_atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                         cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat,rad)

      implicit none

      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(inout) :: u_2d
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: rho_2d, pp_2d, qv_2d, zz_2d
      real (kind=RKIND), dimension(nz1,nlat-1), intent(in) :: zx_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(in) :: fzm, fzp, rdzw
      real (kind=RKIND), intent(in) :: cf1, cf2, cf3, dlat, rad

      !local variable
      real (kind=RKIND), dimension(nz1,nlat-1) :: pgrad, ru, u
      real (kind=RKIND), dimension(nlat-1) :: f
      real (kind=RKIND), dimension(nz1+1)  :: dpzx

   !   real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND) :: rdx, qtot, r_earth, phi
      integer :: k,i, itr

      r_earth  = rad
      omega_e = omega
      rdx = 1./(dlat*r_earth)

      do i=1,nlat-1
        do k=1,nz1
          pgrad(k,i) = rdx*(pp_2d(k,i+1)/zz_2d(k,i+1)-pp_2d(k,i)/zz_2d(k,i))
        end do

        dpzx(:) = 0.

        k=1
        dpzx(k) = .5*zx_2d(k,i)*(cf1*(pp_2d(k  ,i+1)+pp_2d(k  ,i))        &
                                +cf2*(pp_2d(k+1,i+1)+pp_2d(k+1,i))        &
                                +cf3*(pp_2d(k+2,i+1)+pp_2d(k+2,i)))
        do k=2,nz1
           dpzx(k) = .5*zx_2d(k,i)*(fzm(k)*(pp_2d(k  ,i+1)+pp_2d(k  ,i))   &
                                   +fzp(k)*(pp_2d(k-1,i+1)+pp_2d(k-1,i)))
        end do

        do k=1,nz1
           pgrad(k,i) = pgrad(k,i) - rdzw(k)*(dpzx(k+1)-dpzx(k))
        end do
      end do


      !initial value of v and rv -> that is from analytic sln.
      do i=1,nlat-1
         do k=1,nz1
            u(k,i) = .5*(u_2d(k,i)+u_2d(k,i+1))
            ru(k,i) = u(k,i)*(rho_2d(k,i)+rho_2d(k,i+1))*.5
         end do
      end do
   
      call iamas_log_write('MAX U wind before REBALANCING ----> $r', realArgs=(/maxval(abs(u))/))
   
      !re-calculate geostrophic wind using iteration 
      do itr=1,50
      do i=1,nlat-1
         phi = (lat_2d(i)+lat_2d(i+1))/2.
         f(i) = 2.*omega_e*sin(phi)
         do k=1,nz1
            if (f(i).eq.0.) then
              ru(k,i) = 0.
            else
              qtot = .5*(qv_2d(k,i)+qv_2d(k,i+1))
              ru(k,i) = - ( 1./(1.+qtot)*pgrad(k,i) + tan(phi)/r_earth*u(k,i)*ru(k,i) )/f(i)
            end if
              u(k,i) = ru(k,i)*2./(rho_2d(k,i)+rho_2d(k,i+1))
         end do
      end do
      end do
   
      call iamas_log_write('MAX U wind after REBALANCING ----> $r', realArgs=(/maxval(abs(u))/))
   
      !update 2d ru
      do i=2,nlat-1
        do k=1,nz1
          u_2d(k,i) = (ru(k,i-1)+ru(k,i))*.5
        end do
      end do

      i=1
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i+1)-u_2d(k,i+2))*.5
      end do
      i=nlat
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i-1)-u_2d(k,i-2))*.5
      end do

   end subroutine init_atm_recompute_geostrophic_wind


   subroutine init_atm_case_squall_line(dminfo, mesh, nCells, nVertLevels, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup squall line and supercell test case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (iamas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, cell1, cell2
      integer, pointer :: nEdges, nVertices, maxEdges, nCellsSolve
      integer, pointer :: index_qv

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znu, znw, znwc, znwv
      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znuc, znuv

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND), dimension(nVertLevels, nCells) :: relhum, thi, tbi, cqwb

      real (kind=RKIND) ::  r, xnutr
      real (kind=RKIND) ::  ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(nVertLevels ) :: qvb
      real (kind=RKIND), dimension(nVertLevels ) :: t_init_1d

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2
      real (kind=RKIND), pointer :: cf1, cf2, cf3
      real (kind=RKIND) :: ztr, thetar, ttr, thetas, um, us, zts, pitop, pibtop, ptopb, ptop, rcp, rcv, p0
      real (kind=RKIND) :: radx, radz, zcent, xmid, delt, xloc, rad, yloc, ymid, a_scale
      real (kind=RKIND) :: pres, temp, es, qvs

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius

      real (kind=RKIND), dimension(:,:), pointer :: t_init, w, rw, v, rho, theta
      real (kind=RKIND), dimension(:), pointer :: u_init, qv_init, angleEdge, fEdge, fVertex

      call iamas_pool_get_array(mesh, 'xCell', xCell)
      call iamas_pool_get_array(mesh, 'yCell', yCell)
      call iamas_pool_get_array(mesh, 'zCell', zCell)
      call iamas_pool_get_array(mesh, 'xEdge', xEdge)
      call iamas_pool_get_array(mesh, 'yEdge', yEdge)
      call iamas_pool_get_array(mesh, 'zEdge', zEdge)
      call iamas_pool_get_array(mesh, 'xVertex', xVertex)
      call iamas_pool_get_array(mesh, 'yVertex', yVertex)
      call iamas_pool_get_array(mesh, 'zVertex', zVertex)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call iamas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

      call iamas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call iamas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      !
      ! Scale all distances
      !

      a_scale = 1.0

      xCell(:) = xCell(:) * a_scale
      yCell(:) = yCell(:) * a_scale
      zCell(:) = zCell(:) * a_scale
      xVertex(:) = xVertex(:) * a_scale
      yVertex(:) = yVertex(:) * a_scale
      zVertex(:) = zVertex(:) * a_scale
      xEdge(:) = xEdge(:) * a_scale
      yEdge(:) = yEdge(:) * a_scale
      zEdge(:) = zEdge(:) * a_scale
      dvEdge(:) = dvEdge(:) * a_scale
      dcEdge(:) = dcEdge(:) * a_scale
      areaCell(:) = areaCell(:) * a_scale**2.0
      areaTriangle(:) = areaTriangle(:) * a_scale**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * a_scale**2.0

      call iamas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call iamas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call iamas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz1 = nVertLevels
      nz = nz1 + 1

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_array(mesh, 'dzu', dzu)
      call iamas_pool_get_array(mesh, 'rdzu', rdzu)
      call iamas_pool_get_array(mesh, 'fzm', fzm)
      call iamas_pool_get_array(mesh, 'fzp', fzp)
      call iamas_pool_get_array(mesh, 'zxu', zxu)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(mesh, 'hx', hx)
      call iamas_pool_get_array(mesh, 'dss', dss)
      call iamas_pool_get_array(mesh, 't_init', t_init)
      call iamas_pool_get_array(mesh, 'u_init', u_init)
      call iamas_pool_get_array(mesh, 'qv_init', qv_init)
      call iamas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call iamas_pool_get_array(mesh, 'fEdge', fEdge)
      call iamas_pool_get_array(mesh, 'fVertex', fVertex)

      call iamas_pool_get_array(mesh, 'cf1', cf1)
      call iamas_pool_get_array(mesh, 'cf2', cf2)
      call iamas_pool_get_array(mesh, 'cf3', cf3)

      call iamas_pool_get_array(diag, 'pressure_base', ppb)
      call iamas_pool_get_array(diag, 'exner_base', pb)
      call iamas_pool_get_array(diag, 'rho_base', rb)
      call iamas_pool_get_array(diag, 'theta_base', tb)
      call iamas_pool_get_array(diag, 'rtheta_base', rtb)
      call iamas_pool_get_array(diag, 'exner', p)
      call iamas_pool_get_array(diag, 'cqw', cqw)

      call iamas_pool_get_array(diag, 'pressure_p', pp)
      call iamas_pool_get_array(diag, 'rho_p', rr)
      call iamas_pool_get_array(diag, 'rtheta_p', rt)
      call iamas_pool_get_array(diag, 'ru', ru)
      call iamas_pool_get_array(diag, 'rw', rw)
      call iamas_pool_get_array(diag, 'v', v)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag, 'theta', theta)

      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(state, 'theta_m', t)
      call iamas_pool_get_array(state, 'u', u)
      call iamas_pool_get_array(state, 'w', w)
      call iamas_pool_get_array(state, 'scalars', scalars)

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      xnutr = 0.
      zd = 12000.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

     call iamas_log_write(' point 1 in test case setup ')

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,nCells
        do k=1,nz
          hx(k,iCell) = 0.  ! squall line or supercell on flat plane
        end do
      end do

      !     metrics for hybrid coordinate and vertical stretching

      str = 1.0
      zt = 20000.
      dz = zt/float(nz1)

      call iamas_log_write(' hx computation complete ')

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k)
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate

!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!            call iamas_log_write(' k, zc, zw, ah = $i $r $r $r', intArgs=(/k/), realArgs=(/zc(k),zw(k),ah(k)/))
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2)
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3)
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      do iCell=1,nCells
        do k=1,nz
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

!
! convective initialization
!
         ztr    = 12000.
         thetar = 343.
         ttr    = 213.
         thetas = 300.5

!         call iamas_log_write(' rgas, cp, gravity = $r $r $r', realArgs=(/rgas,cp, gravity/))

      if ( test_case == 4) then ! squall line parameters
         um = 12.
         us = 10.
         zts = 2500.
      else if (test_case == 5) then !supercell parameters
         um = 30.
         us = 15.
         zts = 5000.
      end if

         do i=1,nCells
            do k=1,nz1
               ztemp = .5*(zgrid(k,i)+zgrid(k+1,i))
               if(ztemp .gt. ztr) then
                  t (k,i) = thetar*exp(9.8*(ztemp-ztr)/(1003.*ttr))
                  relhum(k,i) = 0.25
               else
                  t (k,i) = 300.+43.*(ztemp/ztr)**1.25
                  relhum(k,i) = (1.-0.75*(ztemp/ztr)**1.25)
                  if(t(k,i).lt.thetas) t(k,i) = thetas
               end if
               tb(k,i) = t(k,i)
               thi(k,i) = t(k,i)
               tbi(k,i) = t(k,i)
               cqw(k,i) = 1.
               cqwb(k,i) = 1.
            end do
         end do

!         relhum(:,:) = 0.

!  set the velocity field - we are on a plane here.

         do i=1, nEdges
            cell1 = cellsOnEdge(1,i)
            cell2 = cellsOnEdge(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               if(ztemp.lt.zts)  then
                  u(k,i) = um*ztemp/zts
               else
                  u(k,i) = um
               end if
               if(i == 1 ) u_init(k) = u(k,i) - us
               u(k,i) = cos(angleEdge(i)) * (u(k,i) - us)
            end do
            end if
         end do

         call iamas_dmpar_bcast_reals(dminfo, nz1, u_init)

!
!    for reference sounding
!
     do itr=1,30

      pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
      pibtop = 1.-.5*dzw(1)*gravity*(1.+qvb(1))/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t(k,1)+t(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))
         pibtop = pibtop-dzu(k)*gravity/(cp*cqwb(k,1)*.5*(tb(k,1)+tb(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))

         !call iamas_log_write('$i $r $r $r $r', intArgs=(/k/), realArgs=(/pitop,tb(k,1),dzu(k),tb(k,1)/))
      end do
      pitop = pitop-.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
      pibtop = pibtop-.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,1)*zz(nz1,1))

      call iamas_dmpar_bcast_real(dminfo, pitop)
      call iamas_dmpar_bcast_real(dminfo, pibtop)

      ptopb = p0*pibtop**(1./rcp)
      call iamas_log_write('ptopb = $r', realArgs=(/0.01_RKIND*ptopb/))

      do i=1, nCells
         pb(nz1,i) = pibtop+.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,i))/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*cqwb(k+1,i)*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*cqw(k+1,i)*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            ppb(k,i) = p0*(zz(k,i)*rgas*rtb(k,i)/p0)**(cp/cv)
         end do
      end do

     !
     ! update water vapor mixing ratio from humidity profile
     !
      do i= 1,nCells
         do k=1,nz1
            temp     = p(k,i)*thi(k,i)
            pres     = p0*p(k,i)**(1./rcp)
            qvs      = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
            scalars(index_qv,k,i) = min(0.014_RKIND,relhum(k,i)*qvs)
         end do
      end do

      do k=1,nz1
!*********************************************************************
!           QVB = QV INCLUDES MOISTURE IN REFERENCE STATE
!            qvb(k) = scalars(index_qv,k,1)

!           QVB = 0 PRODUCES DRY REFERENCE STATE
            qvb(k) = 0.
!*********************************************************************
      end do

      do i= 1,nCells
         do k=1,nz1
            t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
            tb(k,i) = tbi(k,i)*(1.+1.61*qvb(k))
         end do
         do k=2,nz1
            cqw (k,i) = 1./(1.+.5*(scalars(index_qv,k,i)+scalars(index_qv,k-1,i)))
            cqwb(k,i) = 1./(1.+.5*(qvb(k)+qvb(k-1)))
         end do
      end do

      end do !end of iteration loop

      call iamas_log_write(' base state sounding ')
      call iamas_log_write(' k,     pb,     rb,     tb,     rtb,     t,     rr,      p,    qvb')
      do k=1,nVertLevels
         call iamas_log_write('$i $r $r $r $r $r $r $r $r', intArgs=(/k/), realArgs=(/pb(k,1),rb(k,1),tb(k,1),rtb(k,1),t(k,1),rr(k,1),p(k,1),qvb(k)/))
      end do

!
!     potential temperature perturbation
!
!      delt = -10.
!      delt = -0.01
      delt = 3.
      radx  = 10000.
      radz  = 1500.
      zcent = 1500.

      if (test_case == 4) then          ! squall line prameters
         call iamas_dmpar_max_real(dminfo, maxval(xCell(:)), xmid)
         xmid = xmid * 0.5
         ymid = 0.0          ! Not used for squall line
      else if (test_case == 5) then     ! supercell parameters
         call iamas_dmpar_max_real(dminfo, maxval(xCell(:)), xmid)
         call iamas_dmpar_max_real(dminfo, maxval(yCell(:)), ymid)
         xmid = xmid * 0.5
         ymid = ymid * 0.5
      end if

      do i=1, nCells
        xloc = xCell(i) - xmid
        if (test_case == 4) then
           yloc = 0.                            !squall line setting
        else if (test_case == 5) then
           yloc = yCell(i) - ymid !supercell setting
        end if

        do k = 1,nz1
          ztemp     = .5*(zgrid(k+1,i)+zgrid(k,i))
          rad =sqrt((xloc/radx)**2+(yloc/radx)**2+((ztemp-zcent)/radz)**2)
          if(rad.lt.1)  then
            thi(k,i) = thi(k,i) + delt*cos(.5*pii*rad)**2
          end if
           t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
        end do
      end do

      do itr=1,30

        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t (k,1)+t (k-1,1)) &
                                                  *.5*(zz(k,1)+zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)
        call iamas_log_write('ptop  = $r $r', realArgs=(/0.01_RKIND*ptop, 0.01_RKIND*ptopb/))

        call iamas_dmpar_bcast_real(dminfo, ptop)

        do i = 1, nCells

          pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
          do k=nz1-1,1,-1
!             pp(k,i) = pp(k+1,i)+.5*dzu(k+1)*gravity*                   &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)  &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))
               pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*(    &
                            fzm(k+1)*(rb(k+1,i)*(scalars(index_qv,k+1,i)-qvb(k+1))    &
                                     +rr(k+1,i)*(1.+scalars(index_qv,k+1,i)))         &
                           +fzp(k+1)*(rb(k  ,i)*(scalars(index_qv,k  ,i)-qvb(k))      &
                                     +rr(k  ,i)*(1.+scalars(index_qv,k  ,i))))
          end do
          if (itr==1.and.i==1) then
             do k=1,nz1
                call iamas_log_write('pp-check $r', realArgs=(/pp(k,i)/))
             end do
          end if
          do k=1,nz1
             rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                     -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)
             p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
             rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
          end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      do k=1,nz1
         qv_init(k) = scalars(index_qv,k,1)
      end do

      t_init_1d(:) = t(:,1)
      call iamas_dmpar_bcast_reals(dminfo, nz1, t_init_1d)
      call iamas_dmpar_bcast_reals(dminfo, nz1, qv_init)

      do i=1,nCells
         do k=1,nz1
            t_init(k,i) = t_init_1d(k)
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do
      end do

      do i=1,nEdges
        cell1 = cellsOnEdge(1,i)
        cell2 = cellsOnEdge(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)
          end do
        end if
      end do


      !
      !  we are assuming w and rw are zero for this initialization
      !  i.e., no terrain
      !
       rw = 0.0
       w = 0.0

       zb = 0.0
       zb3 = 0.0

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,nEdges
         fEdge(iEdge) = 0.0
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 0.0
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, nVertLevels
                 v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
              end do
            end if
         end do
      end do

     ! call iamas_log_write(' k,u_init, t_init, qv_init ')
     ! do k=1,nVertLevels
     !   call iamas_log_write('$i $r $r $r', intArgs=(/k/), realArgs=(/u_init(k),t_init(k,1),qv_init(k)/))
     ! end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_squall_line


!----------------------------------------------------------------------------------------------------------


   subroutine init_atm_case_mtn_wave(dminfo, mesh, nCells, nVertLevels, state, diag, configs, init_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (iamas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(inout) :: configs
      integer, intent(in) :: init_case

      real (kind=RKIND), parameter :: t0=288., hm=250.

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, deriv_two, zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, cellsOnCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2, nz1
      integer, pointer :: nEdges, maxEdges, nCellsSolve, nVertices
      integer, pointer :: index_qv

      real (kind=RKIND) :: ptop, pitop, ptopb, p0, flux, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND) :: ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(nVertLevels, nCells) :: relhum
      real (kind=RKIND) :: es, qvs, xnutr, ptemp
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2
      real (kind=RKIND) :: um, vm,rcp, rcv
      real (kind=RKIND) :: xmid, temp, pres, a_scale

      real (kind=RKIND) :: xi, xa, xc, xla, zinv, xn2, xn2m, xn2l, sm, dzh, dzht, dzmin, z_edge, z_edge3

      integer, dimension(nCells, 2) :: next_cell
      real (kind=RKIND),  dimension(nCells) :: hxzt
      logical, parameter :: terrain_smooth = .false.

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_theta_adv_order

      real (kind=RKIND), pointer :: cf1, cf2, cf3

      real (kind=RKIND), dimension(:,:), pointer :: t_init, w, rw, v, rho, theta
      real (kind=RKIND), dimension(:), pointer :: u_init, v_init, angleEdge, fEdge, fVertex


      call iamas_pool_get_array(mesh, 'xCell', xCell)
      call iamas_pool_get_array(mesh, 'yCell', yCell)
      call iamas_pool_get_array(mesh, 'zCell', zCell)
      call iamas_pool_get_array(mesh, 'xEdge', xEdge)
      call iamas_pool_get_array(mesh, 'yEdge', yEdge)
      call iamas_pool_get_array(mesh, 'zEdge', zEdge)
      call iamas_pool_get_array(mesh, 'xVertex', xVertex)
      call iamas_pool_get_array(mesh, 'yVertex', yVertex)
      call iamas_pool_get_array(mesh, 'zVertex', zVertex)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call iamas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

      call iamas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call iamas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call iamas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)

      call iamas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call iamas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call iamas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call iamas_pool_get_array(mesh, 't_init', t_init)
      call iamas_pool_get_array(mesh, 'u_init', u_init)
      call iamas_pool_get_array(mesh, 'v_init', v_init)
      call iamas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call iamas_pool_get_array(mesh, 'fEdge', fEdge)
      call iamas_pool_get_array(mesh, 'fVertex', fVertex)

      !
      ! Scale all distances
      !
      a_scale = 1.0

      xCell(:) = xCell(:) * a_scale
      yCell(:) = yCell(:) * a_scale
      zCell(:) = zCell(:) * a_scale
      xVertex(:) = xVertex(:) * a_scale
      yVertex(:) = yVertex(:) * a_scale
      zVertex(:) = zVertex(:) * a_scale
      xEdge(:) = xEdge(:) * a_scale
      yEdge(:) = yEdge(:) * a_scale
      zEdge(:) = zEdge(:) * a_scale
      dvEdge(:) = dvEdge(:) * a_scale
      dcEdge(:) = dcEdge(:) * a_scale
      areaCell(:) = areaCell(:) * a_scale**2.0
      areaTriangle(:) = areaTriangle(:) * a_scale**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * a_scale**2.0


      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call iamas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz1 = nVertLevels
      nz = nz1 + 1

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_array(mesh, 'dzu', dzu)
      call iamas_pool_get_array(mesh, 'rdzu', rdzu)
      call iamas_pool_get_array(mesh, 'fzm', fzm)
      call iamas_pool_get_array(mesh, 'fzp', fzp)
      call iamas_pool_get_array(mesh, 'zxu', zxu)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(mesh, 'hx', hx)
      call iamas_pool_get_array(mesh, 'dss', dss)

      call iamas_pool_get_array(mesh, 'cf1', cf1)
      call iamas_pool_get_array(mesh, 'cf2', cf2)
      call iamas_pool_get_array(mesh, 'cf3', cf3)

      call iamas_pool_get_array(diag, 'pressure_base', ppb)
      call iamas_pool_get_array(diag, 'exner_base', pb)
      call iamas_pool_get_array(diag, 'rho_base', rb)
      call iamas_pool_get_array(diag, 'theta_base', tb)
      call iamas_pool_get_array(diag, 'rtheta_base', rtb)
      call iamas_pool_get_array(diag, 'exner', p)
      call iamas_pool_get_array(diag, 'cqw', cqw)
      call iamas_pool_get_array(diag, 'pressure_p', pp)
      call iamas_pool_get_array(diag, 'rho_p', rr)
      call iamas_pool_get_array(diag, 'rtheta_p', rt)
      call iamas_pool_get_array(diag, 'ru', ru)
      call iamas_pool_get_array(diag, 'rw', rw)
      call iamas_pool_get_array(diag, 'v', v)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag, 'theta', theta)

      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(state, 'theta_m', t)
      call iamas_pool_get_array(state, 'u', u)
      call iamas_pool_get_array(state, 'w', w)
      call iamas_pool_get_array(state, 'scalars', scalars)

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      xnutr = 0.1
      zd = 10500.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

      ! for hx computation
      xa = 5000. !SHP - should be changed based on grid distance
      xla = 4000.
      xc = maxval (xCell(:))/2.

      !     metrics for hybrid coordinate and vertical stretching
      str = 1.0
      zt = 21000.
      dz = zt/float(nz1)
!      call iamas_log_write(' dz = $r', realArgs=(/dz/))

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k)
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate

!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!            call iamas_log_write(' k, zc, zw, ah = $i $r $r $r', intArgs=(/k/), realArgs=(/zc(k),zw(k),ah(k)/))
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      d1  = .5*dzw(1)
      d2  = dzw(1)+.5*dzw(2)
      d3  = dzw(1)+dzw(2)+.5*dzw(3)
      !cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      cof1 = (2.*dzu(2)+dzu(3))/(dzu(2)+dzu(3))*dzw(1)/dzu(2)
      cof2 =     dzu(2)        /(dzu(2)+dzu(3))*dzw(1)/dzu(3)
      cf1  = fzp(2) + cof1
      cf2  = fzm(2) - cof1 - cof2
      cf3  = cof2

! setting for terrain
      do iCell=1,nCells
         xi = xCell(iCell)
         !====1. for pure cosine mountain
         ! if(abs(xi-xc).ge.2.*xa)  then
         !    hx(1,iCell) = 0.
         ! else
         !    hx(1,iCell) = hm*cos(.5*pii*(xi-xc)/(2.*xa))**2.
         ! end if

         !====2. for cosine mountain
         !if(abs(xi-xc).lt.xa)  THEN
         !     hx(1,iCell) = hm*cos(pii*(xi-xc)/xla)**2. *cos(.5*pii*(xi-xc)/xa )**2.
         ! else
         !    hx(1,iCell) = 0.
         ! end if

         !====3. for shock mountain
         hx(1,iCell) = hm*exp(-((xi-xc)/xa)**2)*cos(pii*(xi-xc)/xla)**2.

         hx(nz,iCell) = zt

!***** SHP -> get the temporary point information for the neighbor cell ->> should be changed!!!!!
         do i=1,nCells
            !option 1
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)-sqrt(3.)*dcEdge(1)) next_cell(iCell,1) = i
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)+sqrt(3.)*dcEdge(1)) next_cell(iCell,2) = i
            !option 2
            next_cell(iCell,1) = iCell - 8 ! note ny=4
            next_cell(iCell,2) = iCell + 8 ! note ny=4

            if (xCell(iCell).le. 3.*dcEdge(1)) then
                next_cell(iCell,1) = 1
            else if (xCell(iCell).ge. maxval(xCell(:))-3.*dcEdge(1)) then
                next_cell(iCell,2) = 1
            end if

         end do
      end do
      
      call iamas_log_write(' hx computation complete ')

      if (terrain_smooth) then
         call iamas_log_write('***************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Please contact the iAMAS-A developers for up-to-date terrain-smoothing code.', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Otherwise, set terrain_smooth=.false. in the mountain wave test case',        messageType=iAMAS_LOG_ERR)
         call iamas_log_write('   initialization routine and re-compile.',                                   messageType=iAMAS_LOG_ERR)
         call iamas_log_write('***************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      do iCell=1,nCells
        do k=1,nz
            if (terrain_smooth) then
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)
            else
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(1,iCell)/zt)+hx(1,iCell)) &
                           + (1.-ah(k)) * zc(k)
            end if
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

      call iamas_log_write(' grid metrics setup complete ')

!
! mountain wave initialization
!
         !SHP-original
         !zinv = 1000.
         !SHP-schar case
         zinv = 3000.

         xn2  = 0.0001
         xn2m = 0.0000
         xn2l = 0.0001

         vm = 10.
         um = 0.

         do k=1,nz1
           v_init(k) = vm
           u_init(k) = um
         end do

         do i=1,nCells
            do k=1,nz1
               ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
               tb(k,i) =  t0*(1. + xn2m/gravity*ztemp)
               if(ztemp .le. zinv) then
                  t (k,i) = t0*(1.+xn2l/gravity*ztemp)
               else
                  t (k,i) = t0*(1.+xn2l/gravity*zinv+xn2/gravity*(ztemp-zinv))
               end if
                  relhum(k,i) = 0.
            end do
         end do

!  set the velocity field - we are on a plane here.

         do i=1, nEdges
            cell1 = cellsOnEdge(1,i)
            cell2 = cellsOnEdge(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               u(k,i) =   vm*sin(angleEdge(i)) + um*cos(angleEdge(i))
            end do
            end if
         end do

!
!     reference sounding based on dry atmosphere
!
      pitop = 1.-.5*dzw(1)*gravity/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*(fzm(k)*tb(k,1)+fzp(k)*tb(k-1,1))   &
                                         *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
      end do
      pitop = pitop-.5*dzw(nz1)*gravity/(cp*tb(nz1,1)*zz(nz1,1))
      ptopb = p0*pitop**(1./rcp)

      do i=1, nCells
         pb(nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            cqw(k,i) = 1.
         end do
      end do

       call iamas_log_write(' ***** base state sounding ***** ')
       call iamas_log_write('k       pb        p         rb         rtb         rr          tb          t')
       do k=1,nVertLevels
          call iamas_log_write('$i $r $r $r $r $r $r $r', intArgs=(/k/), realArgs=(/pb(k,1),p(k,1),rb(k,1),rtb(k,1),rr(k,1),tb(k,1),t(k,1)/))
       end do

       scalars(index_qv,:,:) = 0.

!-------------------------------------------------------------------
!     ITERATIONS TO CONVERGE MOIST SOUNDING
      do itr=1,30
        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))

        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*(fzm(k)*t (k,1)+fzp(k)*t (k-1,1)) &
                                                   *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)

        do i = 1, nCells

           pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
           do k=nz1-1,1,-1
              pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*                   &
                            (fzm(k)*(rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))  &
                            +fzp(k)*(rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i)))
           end do
           do k=1,nz1
              rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                      -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)
              p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
              rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
           end do
!
!     update water vapor mixing ratio from humitidty profile
!
           do k=1,nz1
              temp   = p(k,i)*t(k,i)
              pres   = p0*p(k,i)**(1./rcp)
              qvs    = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
              scalars(index_qv,k,i) = min(0.014_RKIND,relhum(k,i)*qvs)
           end do

           do k=1,nz1
              t (k,i) = t(k,i)*(1.+1.61*scalars(index_qv,k,i))
           end do
           do k=2,nz1
              cqw(k,i) = 1./(1.+.5*( scalars(index_qv,k-1,i)  &
                                    +scalars(index_qv,k  ,i)))
           end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      call iamas_log_write(' *** sounding for the simulation ***')
      call iamas_log_write('    z       theta       pres         qv       rho_m        u        rr')
      do k=1,nz1
         call iamas_log_write('$r $r $r $r $r $r $r', realArgs=(/.5*(zgrid(k,1)+zgrid(k+1,1))/1000.,   &
                       t(k,1)/(1.+1.61*scalars(index_qv,k,1)),        &
                       .01*p0*p(k,1)**(1./rcp),                       &
                       1000.*scalars(index_qv,k,1),                   &
                       (rb(k,1)+rr(k,1))*(1.+scalars(index_qv,k,1)),  &
                       u_init(k), rr(k,1)/))
      end do

      do i=1,ncells
         do k=1,nz1
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do

        do k=1,nz1
            t_init(k,i) = t(k,i)
        end do
      end do

      do i=1,nEdges
        cell1 = cellsOnEdge(1,i)
        cell2 = cellsOnEdge(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)
          end do
        end if
      end do

!
!     pre-calculation z-metric terms in omega eqn.
!
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do

                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  else
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

            end do

         end if
       end do

!     for including terrain
      w(:,:) = 0.0
      rw(:,:) = 0.0

!
!     calculation of omega, rw = zx * ru + zz * rw
!

      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then
               rw(k,cell2) = rw(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do


      do iEdge=1,nEdges
         fEdge(iEdge) = 0.
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, nVertLevels
                 v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
              end do
            end if
         end do
      end do

!      do k=1,nVertLevels
!        call iamas_log_write(' k,u_init, t_init = $i $r $r $r', intArgs=(/k/), realArgs=(/u_init(k),t_init(k,1)/))
!      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_mtn_wave


   subroutine init_atm_case_gfs(block, mesh, nCells, nEdges, nVertLevels, fg, state, diag, diag_physics, init_case, dims, configs)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Real-data test case using GFS data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use iamas_dmpar
      use init_atm_read_met
      use init_atm_llxy
      use init_atm_hinterp
      use iamas_hash
      use iamas_atmphys_constants, only : svpt0, svp1, svp2, svp3

      implicit none

      type (block_type), intent(inout), target :: block
      type (iamas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: fg
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(inout):: diag_physics
      integer, intent(in) :: init_case
      type (iamas_pool_type), intent(inout):: dims
      type (iamas_pool_type), intent(inout):: configs

      type (parallel_info), pointer :: parinfo
      type (dm_info), pointer :: dminfo

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

!      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      type (met_data) :: field
      type (proj_info) :: proj

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: vert_level, latPoints, lonPoints, ter
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, pressure_b, exner_b, rho_zz, rho_b, rho_p, phi_b, phi_total
      real (kind=RKIND), dimension(:,:), pointer :: theta_b, rtheta_b, exner, pressure_p, dss, theta_m, rtheta_p
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:), pointer :: destField1d
      real (kind=RKIND), dimension(:,:), pointer :: destField2d
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      real (kind=RKIND) :: target_z
      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
      integer, pointer :: nCellsSolve, nz1
      integer :: nInterpPoints, ndims

      integer :: nfglevels_actual
      integer, pointer :: index_qv

      integer, dimension(5) :: interp_list
      real (kind=RKIND) :: maskval
      real (kind=RKIND) :: msgval
      real (kind=RKIND) :: fillval
      integer :: masked

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, edgesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: v
      real (kind=RKIND), dimension(:,:), pointer :: sorted_arr
      integer, dimension(:), pointer :: bdyMaskCell

      type (field1DReal), pointer :: tempField
      type (field1DReal), pointer :: ter_field
      type (field1DReal), target :: tempFieldTarget

      real(kind=RKIND), dimension(:), pointer :: hs, hs1, sm0
      real(kind=RKIND) :: hm, hm_global, zh, dzmin, dzmina, dzmina_global, dzminf, dzminf_global, sm
      real(kind=RKIND) :: dcsum
      integer :: nsmterrain, kz, sfc_k
      logical :: hybrid, smooth,cam_lev30,cam_lev32

      integer :: it
      real (kind=RKIND) :: p_check

      ! For interpolating terrain and land use
      integer :: nx, ny
      integer :: istatus

      real (kind=RKIND), allocatable, dimension(:,:) :: rslab, maskslab
      integer, dimension(:), pointer :: mask_array
      integer, dimension(nEdges), target :: edge_mask
      character (len=StrKIND) :: fname
      logical :: is_sfc_field

      real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r
      real (kind=RKIND) :: lat, lon, x, y

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND), dimension(nVertLevels, nCells) :: rel_hum, temperature, qv
      real (kind=RKIND) :: ptmp, es, rs, rgas_moist, qvs, xnutr, znut, ptemp, rcv
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(nVertLevels + 1) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(:), pointer :: specified_zw
      real (kind=RKIND), dimension(nVertLevels + 1) :: sh, zw, ah
      real (kind=RKIND), dimension(nVertLevels) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(nVertLevels) :: eta, etav, teta, ppi, tt

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      integer, parameter :: nlat=361
      real (kind=RKIND), dimension(nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
      real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nlat, nVertLevels) :: u_2d, etavs_2d
      real (kind=RKIND), dimension(nVertLevels + 1) :: fsum
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND) :: alt, als, zetal, zl

      !  calculation of the water vapor mixing ratio:
      real (kind=RKIND) :: sh_max,sh_min,global_sh_max,global_sh_min

      character (len=StrKIND), pointer :: config_met_prefix
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_source_file
      character (len=StrKIND), pointer :: config_source_opt

      ! added by qkyang to support anthropogenic heating emissions>>>>>>>>>>>>>>
      character (len=StrKIND), pointer :: config_AH_emis_file
      logical, pointer :: is_ah_emis
      ! added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<<
      logical, pointer :: config_met_interp
      logical, pointer :: config_vertical_grid
      integer, pointer :: config_nsmterrain
      integer, pointer :: config_nsm
      real (kind=RKIND), pointer :: config_dzmin
      real (kind=RKIND), pointer :: config_ztop
      logical, pointer :: config_tc_vertical_grid
      character (len=StrKIND), pointer :: config_specified_zeta_levels
      logical, pointer :: config_use_spechumd
      integer, pointer :: config_nfglevels
      integer, pointer :: config_nfgsoillevels
      logical, pointer :: config_smooth_surfaces
      logical, pointer :: config_cam_lev30 
      logical, pointer :: config_cam_lev32 
      integer, pointer :: config_theta_adv_order
      real (kind=RKIND), pointer :: config_coef_3rd_order
      logical, pointer :: config_blend_bdy_terrain
      logical, pointer :: config_log_pressure ! jwfeng: If we are using log pressure (phi = ln(p/p0))

      character (len=StrKIND), pointer :: config_extrap_airtemp
      integer :: extrap_airtemp

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), pointer :: cf1, cf2, cf3
      integer, dimension(:), pointer :: landmask

      real (kind=RKIND), dimension(:,:), pointer :: dzs_fg
      real (kind=RKIND), dimension(:,:), pointer :: zs_fg

      real (kind=RKIND), dimension(:), pointer :: sst
      real (kind=RKIND), dimension(:), pointer :: seaice
      real (kind=RKIND), dimension(:), pointer :: xice
      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: w
      real (kind=RKIND), dimension(:,:), pointer :: theta
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:), pointer :: relhum
      real (kind=RKIND), dimension(:,:), pointer :: spechum
      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:), pointer :: precipw
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructY
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructZ
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructMeridional

      real (kind=RKIND), dimension(:), pointer :: psfc
      real (kind=RKIND), dimension(:), pointer :: skintemp
      real (kind=RKIND), dimension(:), pointer :: snow
      real (kind=RKIND), dimension(:), pointer :: snowc
      real (kind=RKIND), dimension(:,:), pointer :: u_fg
      real (kind=RKIND), dimension(:,:), pointer :: v_fg
      real (kind=RKIND), dimension(:,:), pointer :: z_fg
      real (kind=RKIND), dimension(:,:), pointer :: t_fg
      real (kind=RKIND), dimension(:,:), pointer :: rh_fg
      real (kind=RKIND), dimension(:,:), pointer :: sh_fg
      real (kind=RKIND), dimension(:,:), pointer :: gfs_z
      real (kind=RKIND), dimension(:,:), pointer :: p_fg
      real (kind=RKIND), dimension(:,:), pointer :: st_fg
      real (kind=RKIND), dimension(:,:), pointer :: sm_fg
      real (kind=RKIND), dimension(:,:), pointer :: erod
      real (kind=RKIND), dimension(:), pointer :: AH_emis! added by qkyang to support anthropogenic heating emissions
      real (kind=RKIND), dimension(:), pointer :: dms_0
      real (kind=RKIND), dimension(:), pointer :: soilz
      type (hashtable) :: level_hash
      logical :: too_many_fg_levs
      logical, allocatable, dimension(:) :: found_u, found_v, found_t, found_moisture, found_z
      logical :: found_soilhgt
      integer :: level_value
      integer :: missing_required_fields

      ! For outputting surface fields u10, v10, q2, rh2, and t2m from first-guess data
      real (kind=RKIND), dimension(:), pointer :: u10
      real (kind=RKIND), dimension(:), pointer :: v10
      real (kind=RKIND), dimension(:), pointer :: q2
      real (kind=RKIND), dimension(:), pointer :: rh2
      real (kind=RKIND), dimension(:), pointer :: t2m

      character (len=StrKIND) :: errstring

      call iamas_pool_get_config(configs, 'config_source_file' , config_source_file)
      call iamas_pool_get_config(configs, 'config_source_opt', config_source_opt)
      ! added by qkyang to support anthropogenic heating emissions>>>>>>>>>>>>>>
      call iamas_pool_get_config(configs, 'config_AH_emis_file' , config_AH_emis_file)
      call iamas_pool_get_config(configs, 'is_ah_emis' , is_ah_emis)
      ! added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<<
      call iamas_pool_get_config(configs, 'config_met_prefix', config_met_prefix)
      call iamas_pool_get_config(configs, 'config_start_time', config_start_time)
      call iamas_pool_get_config(configs, 'config_met_interp', config_met_interp)
      call iamas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call iamas_pool_get_config(configs, 'config_nsmterrain', config_nsmterrain)
      call iamas_pool_get_config(configs, 'config_nsm', config_nsm)
      call iamas_pool_get_config(configs, 'config_dzmin', config_dzmin)
      call iamas_pool_get_config(configs, 'config_ztop', config_ztop)
      call iamas_pool_get_config(configs, 'config_tc_vertical_grid', config_tc_vertical_grid)
      call iamas_pool_get_config(configs, 'config_specified_zeta_levels', config_specified_zeta_levels)
      call iamas_pool_get_config(configs, 'config_use_spechumd', config_use_spechumd)
      call iamas_pool_get_config(configs, 'config_nfglevels', config_nfglevels)
      call iamas_pool_get_config(configs, 'config_nfgsoillevels', config_nfgsoillevels)
      call iamas_pool_get_config(configs, 'config_smooth_surfaces', config_smooth_surfaces)
      call iamas_pool_get_config(configs, 'config_cam_lev30', config_cam_lev30)
      call iamas_pool_get_config(configs, 'config_cam_lev32', config_cam_lev32)
      call iamas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)
      call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call iamas_pool_get_config(configs, 'config_blend_bdy_terrain', config_blend_bdy_terrain)
      call iamas_pool_get_config(configs, 'config_log_pressure', config_log_pressure)
      
      call iamas_pool_get_config(configs, 'config_extrap_airtemp', config_extrap_airtemp)
      if (trim(config_extrap_airtemp) == 'constant') then
         extrap_airtemp = 0
      else if (trim(config_extrap_airtemp) == 'linear') then
         extrap_airtemp = 1
      else if (trim(config_extrap_airtemp) == 'lapse-rate') then
         extrap_airtemp = 2
      else
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('* Invalid value for namelist variable config_extrap_airtemp *', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
      call iamas_log_write("Using option '" // trim(config_extrap_airtemp) // "' for vertical extrapolation of temperature")

      parinfo => block % parinfo
      dminfo => block % domain % dminfo

      call iamas_pool_get_field(mesh, 'ter', ter_field)

      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call iamas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call iamas_pool_get_array(mesh, 'angleEdge', angleEdge)

      call iamas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_array(mesh, 'dzu', dzu)
      call iamas_pool_get_array(mesh, 'rdzu', rdzu)
      call iamas_pool_get_array(mesh, 'fzm', fzm)
      call iamas_pool_get_array(mesh, 'fzp', fzp)
      call iamas_pool_get_array(mesh, 'zxu', zxu)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(mesh, 'hx', hx)
      call iamas_pool_get_array(mesh, 'ter', ter)
      call iamas_pool_get_array(mesh, 'dss', dss)

      call iamas_pool_get_array(diag, 'exner_base', exner_b)
      call iamas_pool_get_array(diag, 'rho_base', rho_b)
      call iamas_pool_get_array(diag, 'theta_base', theta_b)
      call iamas_pool_get_array(diag, 'rtheta_base', rtheta_b)
      call iamas_pool_get_array(diag, 'exner', exner)
      call iamas_pool_get_array(diag, 'pressure_base', pressure_b)
      call iamas_pool_get_array(diag, 'pressure_p', pressure_p)
      call iamas_pool_get_array(diag, 'phi_base', phi_b)
      call iamas_pool_get_array(diag, 'phi_total', phi_total)
      call iamas_pool_get_array(diag, 'pressure', pressure)
      call iamas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call iamas_pool_get_array(diag, 'relhum', relhum)
      call iamas_pool_get_array(diag, 'spechum', spechum)
      call iamas_pool_get_array(diag, 'ru', ru)
      call iamas_pool_get_array(diag, 'rw', rw)

      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(diag, 'rho_p', rho_p)
      call iamas_pool_get_array(state, 'theta_m', theta_m)
      call iamas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call iamas_pool_get_array(state, 'scalars', scalars)
      call iamas_pool_get_array(state, 'u', u)
      call iamas_pool_get_array(state, 'w', w)
      call iamas_pool_get_array(diag, 'theta', theta)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag_physics, 'precipw', precipw)
      call iamas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call iamas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call iamas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      call iamas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call iamas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

      call iamas_pool_get_array(mesh, 'latCell', latCell)
      call iamas_pool_get_array(mesh, 'lonCell', lonCell)
      call iamas_pool_get_array(mesh, 'latEdge', latEdge)
      call iamas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call iamas_pool_get_array(mesh, 'cf1', cf1)
      call iamas_pool_get_array(mesh, 'cf2', cf2)
      call iamas_pool_get_array(mesh, 'cf3', cf3)
      call iamas_pool_get_array(mesh, 'landmask', landmask)

      call iamas_pool_get_array(fg, 'dzs_fg', dzs_fg)
      call iamas_pool_get_array(fg, 'zs_fg', zs_fg)
      call iamas_pool_get_array(fg, 'sst', sst)
      call iamas_pool_get_array(fg, 'xice', xice)
      call iamas_pool_get_array(fg, 'seaice', seaice)
      call iamas_pool_get_array(fg, 'st_fg', st_fg)
      call iamas_pool_get_array(fg, 'sm_fg', sm_fg)
      call iamas_pool_get_array(fg, 'erod', erod)
      call iamas_pool_get_array(fg, 'dms_0', dms_0)
      call iamas_pool_get_array(fg, 'AH_emis', AH_emis)! added by qkyang to support anthropogenic heating emissions
      call iamas_pool_get_array(fg, 'psfc', psfc)
      call iamas_pool_get_array(fg, 'skintemp', skintemp)
      call iamas_pool_get_array(fg, 'snow', snow)
      call iamas_pool_get_array(fg, 'snowc', snowc)
      call iamas_pool_get_array(fg, 'u', u_fg)
      call iamas_pool_get_array(fg, 'v', v_fg)
      call iamas_pool_get_array(fg, 'z', z_fg)
      call iamas_pool_get_array(fg, 't', t_fg)
      call iamas_pool_get_array(fg, 'rh', rh_fg)
      call iamas_pool_get_array(fg, 'sh', sh_fg)
      call iamas_pool_get_array(fg, 'gfs_z', gfs_z)
      call iamas_pool_get_array(fg, 'p', p_fg)

      call iamas_pool_get_dimension(dims, 'nVertLevels', nz1)
      call iamas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
      nz = nz1 + 1

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      rcv = rgas/(cp-rgas)
      omega_e = omega
      p0 = 1.e+05

      scalars(:,:,:) = 0.

      
      !
      ! If requested, blend the terrain along the domain boundaries with terrain from
      ! an intermediate file. For global domains, this routine will have no effect even
      ! if called, since terrain is only blended for cells with bdyMaskCell > 0.
      !
      if (config_blend_bdy_terrain) then
         call iamas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
         call iamas_pool_get_array(mesh, 'ter', ter)

         call blend_bdy_terrain(config_met_prefix, config_start_time, nCells, latCell, lonCell, bdyMaskCell, ter, .false., istatus)
         if (istatus /= 0) then
             call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
             call iamas_log_write('* Blending of terrain along domain boundaries failed!       *', messageType=iAMAS_LOG_ERR)
             call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
         end if
      end if

      if (config_vertical_grid) then

      !
      ! Vertical grid setup
      !
      allocate(hs (nCells+1))
      allocate(hs1(nCells+1))

!     Fourth order smoother for terrain

      nsmterrain = config_nsmterrain

      do i=1,nsmterrain

         do iCell=1,nCells
            hs(iCell) = 0.
            if(ter(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)

                  ! For smoothing at cells along the boundary of the mesh, set the terrain value
                  ! for non-existent neighbors, which map to the "garbage cell", to the same as
                  ! the terrain in the cell being smoothed
                  if (cellsOnCell(j,iCell) == nCells+1) then
                     ter(nCells+1) = ter(iCell)
                  end if

                  hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                        / dcEdge(edgesOnCell(j,iCell))    &
                                        * (ter(cellsOnCell(j,iCell))-ter(iCell))
               end do
            end if
            hs(iCell) = ter(iCell) + 0.216*hs(iCell)
         end do

         do iCell=1,nCells
            ter(iCell) = 0.
            if(hs(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)

                  ! For smoothing at cells along the boundary of the mesh, set the terrain value
                  ! for non-existent neighbors, which map to the "garbage cell", to the same as
                  ! the terrain in the cell being smoothed
                  if (cellsOnCell(j,iCell) == nCells+1) then
                     hs(nCells+1) = hs(iCell)
                  end if

                  ter(iCell) = ter(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                          / dcEdge(edgesOnCell(j,iCell))    &
                                          * (hs(cellsOnCell(j,iCell))-hs(iCell))
               end do
            end if
            ter(iCell) = hs(iCell) - 0.216*ter(iCell)
         end do

         ! note that ther variable ter used throughout this section is a pointer to grid % ter % array, here we are passing ter's parent field
         call iamas_dmpar_exch_halo_field(ter_field)

      end do

      do iCell=1,nCells
         hx(:,iCell) = ter(iCell)
      end do

      hm = maxval(ter(1:nCellsSolve))
      call iamas_dmpar_max_real(dminfo, hm, hm_global)
      hm = hm_global
      call iamas_log_write('max ter = $r', realArgs=(/hm/))

!     Metrics for hybrid coordinate and vertical stretching

      !
      ! If a the name of a file with vertical coordinate values has been specified,
      ! use those values to setup the vertical grid
      !
      if (len_trim(config_specified_zeta_levels) > 0) then

         call iamas_log_write('Setting up vertical grid using levels from '''//trim(config_specified_zeta_levels)//'''')

         if (read_text_array(dminfo, trim(config_specified_zeta_levels), specified_zw) /= 0) then
            call iamas_log_write('Failed to read vertical levels from '''//trim(config_specified_zeta_levels)//'''', &
                                 messageType=iAMAS_LOG_CRIT)
         end if

         if (size(specified_zw) /= nz) then
            call iamas_log_write('In the namelist.init_atmosphere file, config_nvertlevels = $i, but ', intArgs=(/nz1/), &
                                 messageType=iAMAS_LOG_ERR)
            call iamas_log_write('but '''//trim(config_specified_zeta_levels)//''' has $i values.', intArgs=(/size(specified_zw)/), &
                                 messageType=iAMAS_LOG_ERR)
            call iamas_log_write(''''//trim(config_specified_zeta_levels)//''' must contain nVertLevels+1 ($i) values.', intArgs=(/nz/), &
                                 messageType=iAMAS_LOG_CRIT)
         end if

         zw(:) = specified_zw(:)
         zt = zw(nz)

         deallocate(specified_zw)

      !
      ! Otherwise, see if the user has requested to set up the vertical grid as in the iAMAS TC configuration
      !
      else if (config_tc_vertical_grid) then
!     Below can be modified to specify the vertical coordinate distribution -czhao
!++++++++++++++czhao
      str = 1.5
!      str = 1.
         call iamas_log_write('Setting up vertical levels as in 2014 TC experiments')

         zt = config_ztop
         dz = zt/float(nz1)

!... Laura D. Fowler: change the values als,alt,and zetal valid for nVertLevels=41 to values
!    needed for nVertLevels equal or greater than 55.
         if (nVertLevels >= 55) then
            als   =  .075  ! 55 levels.
            alt   = 1.70   ! 55 levels.
            zetal =  .75   ! 55 levels.
         else
            als   =  .075  ! 41 levels.
            alt   = 1.23   ! 41 levels.
            zetal =  .31   ! 41 levels.
         endif
         call iamas_log_write('')
         call iamas_log_write('--- config_tc_vertical_grid = $l', logicArgs=(/config_tc_vertical_grid/))
         call iamas_log_write('--- als   = $r', realArgs=(/als/))
         call iamas_log_write('--- alt   = $r', realArgs=(/alt/))
         call iamas_log_write('--- zetal = $r', realArgs=(/zetal/))
         if (nVertLevels /= 55 .and. nVertLevels /= 41) then
            call iamas_log_write('********************************************************************************')
            call iamas_log_write('* Note: als, alt, and zetal have not been tested for nVertLevels = $i', intArgs=(/nVertLevels/))
            call iamas_log_write('*       The values of these parameters have only been tested with 41 and 55 layers')
            call iamas_log_write('********************************************************************************')
         end if
         call iamas_log_write('')
!... end Laura D. Fowler / 2016-04-12.

         do k=1,nz
            zl = 1.-alt*(1.-zetal)
            if ((real(k-1)/real(nz1)).LT.zetal) then
                 zw(k) =              (als*real(K-1)/real(nz1)      &
                     +(3.*(1.-alt)+2.*(alt-als)*zetal)              &
                                     *((K-1)*dz/(zt*zetal))**2      &
                     -(2.*(1.-alt) +  (alt-als)*zetal)              &
                                     *(real(k-1)*dz/(zt*zetal))**3) * zt
            else
                 zw(K) = (zl+alt*(real(K-1)/real(nz1)-zetal))*zt
            end if
            if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
         end do

      !
      ! Otherwise, use the vertical level configuration from iAMAS v2.0
      !
      else

         call iamas_log_write('Setting up vertical levels as in iAMAS 2.0 and earlier')

         str = 1.5
         zt = config_ztop
         dz = zt/float(nz1)

         do k=1,nz
            zw(k) = (real(k-1)/real(nz1))**str*zt
            if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
         end do
      end if !config_tc_vertical_grid

      cam_lev30 = config_cam_lev30
      if (cam_lev30) then 
!We use specified 30 vertical grids for CAM5.3 ++++++czhao 
      !CAM-L30
      zw(1) = 0.
      zw(2) = 109.69314
      zw(3) = 241.32851
      zw(4) = 394.90024
      zw(5) = 570.41372
      zw(6) = 767.86412
      zw(7) = 987.25538
      zw(8) = 1228.58545
      zw(9) = 1665.17067
      zw(10) = 2297.01418
      zw(11) = 3124.11476
      zw(12) = 4146.47283
      zw(13) = 5334.83465
      zw(14) = 6523.19676
      zw(15) = 7711.55805
      zw(16) = 8899.92028
      zw(17) = 10088.28224
      zw(18) = 11276.64419
      zw(19) = 12465.00619
      zw(20) = 13653.36722
      zw(21) = 14841.72987
      zw(22) = 16030.09128
      zw(23) = 17218.45255
      zw(24) = 18406.81464
      zw(25) = 20180.21649
      zw(26) = 22538.65857
      zw(27) = 25482.14044
      zw(28) = 29156.92293
      zw(29) = 33563.00532
      zw(30) = 38700.38866
      zw(31) = 44569.07186
      zt = zw(31)
!We use specified 30 vertical grids for CAM5.3 --------czhao 
      endif

      cam_lev32 = config_cam_lev32
      if (cam_lev32) then 
!We use specified 32 vertical grids for CAM5.4 ++++++czhao 
      !CAM-L32
      zw(1) = 0.
      zw(2) = 109.69314  
      zw(3) = 241.32851  
      zw(4) = 394.90024  
      zw(5) = 570.41372  
      zw(6) = 767.86412  
      zw(7) = 987.25538  
      zw(8) = 1228.58545 
      zw(9) = 1665.17067 
      zw(10)= 2297.01418 
      zw(11)= 3124.11476 
      zw(12)= 4146.47283 
      zw(13)= 5334.83465 
      zw(14)= 6523.19676 
      zw(15)= 7711.55805 
      zw(16)= 8899.92028 
      zw(17)= 10088.28224
      zw(18)= 11276.64419
      zw(19)= 12465.00618
      zw(20)= 13653.36722
      zw(21)= 14841.72987
      zw(22)= 16030.09128
      zw(23)= 17218.45255
      zw(24)= 18497.99874
      zw(25)= 19789.30527
      zw(26)= 21086.79176
      zw(27)= 22373.60067
      zw(28)= 23673.77514
      zw(29)= 25708.42285
      zw(30)= 29156.92293
      zw(31)= 33563.00532
      zw(32)= 38700.38866
      zw(33)= 44569.07186
      zt = zw(33)
!We use specified 32 vertical grids for CAM5.4 --------czhao 
      endif 

!--------------------------czhao

!     ah(k) governs the transition between terrain-following
!        and pure height coordinates
!           ah(k) = 1           is a smoothed terrain-following coordinate
!           ah(k) = 1.-zw(k)/zt is the basic terrain-following coordinate
!           ah(k) = 0           is a height coordinate

      hybrid = .true.
!      hybrid = .false.

      kz = nz
      if (hybrid) then

         zh = 30000.0
!         zh = 0.5*zt

         do k=1,nz
            if (zw(k) < zh) then
               ah(k) = cos(.5*pii*zw(k)/zh)**6

!!!               ah(k) = ah(k)*(1.-zw(k)/zt)

            else
               ah(k) = 0.
               kz = min(kz,k)
            end if
         end do

      else

         do k=1,nz
            ah(k) = 1.-zw(k)/zt
         end do

      end if

      do k=1,nz
         call iamas_log_write('$i $r $r', intArgs=(/k/), realArgs=(/zw(k), ah(k)/))
      end do

      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2)
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3)
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))


!     Smoothing algorithm for coordinate surfaces 

      smooth = config_smooth_surfaces

      if (smooth) then

         dzmin = config_dzmin

         allocate(sm0(nCells+1))

         do iCell=1,nCells
            dcsum = 0.0
            do j=1,nEdgesOnCell(iCell)
               dcsum = dcsum + dcEdge(edgesOnCell(j,iCell))
            end do
            dcsum = dcsum / real(nEdgesOnCell(iCell))
            sm0(iCell) = max(0.01_RKIND, 0.125 * min(1.0_RKIND, 3000.0_RKIND/dcsum))
         end do

         do k=2,kz-1
            hx(k,:) = hx(k-1,:)
            dzminf = zw(k)-zw(k-1)

            do i=1,config_nsm + k
               do iCell=1,nCells

                  sm = sm0(iCell) * min((3.0*zw(k)/zt)**2.0, 1.0_RKIND)

                  hs1(iCell) = 0.
                  do j = 1,nEdgesOnCell(iCell)

                     ! For smoothing at cells along the boundary of the mesh, set the hx value
                     ! for non-existent neighbors, which map to the "garbage cell", to the same as
                     ! the hx in the cell being smoothed
                     if (cellsOnCell(j,iCell) == nCells+1) then
                        hx(k,nCells+1) = hx(k,iCell)
                     end if

                     hs1(iCell) = hs1(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                           / dcEdge(edgesOnCell(j,iCell))    &
                                           *  (hx(k,cellsOnCell(j,iCell))-hx(k,iCell))
                  end do
                  hs(iCell) = hx(k,iCell) + sm*hs1(iCell)

               end do

               tempField => tempFieldTarget
               tempField % block => block
               tempField % dimSizes(1) = nCells
               tempField % sendList => parinfo % cellsToSend
               tempField % recvList => parinfo % cellsToRecv
               tempField % copyList => parinfo % cellsToCopy
               tempField % array => hs
               tempField % isActive = .true.
               tempField % prev => null()
               tempField % next => null()

               call iamas_dmpar_exch_halo_field(tempField)

               do iCell=1,nCells
                  dzmina = (zw(k) + ah(k)*hs(iCell)) - (zw(k-1) + ah(k-1)*hx(k-1,iCell))
                  if (dzmina > dzmin*(zw(k)-zw(k-1))) then
                     hx(k,iCell) = hs(iCell)
                     if (dzmina < dzminf) dzminf = dzmina
                  end if
               end do

            end do
            call iamas_dmpar_min_real(dminfo, dzminf, dzminf_global)
            call iamas_log_write('$i $i $r $r', intArgs=(/k,i/), realArgs=(/sm,dzminf_global/(zw(k)-zw(k-1))/))
         end do

         deallocate(sm0)

         do k=kz,nz
               hx(k,:) = 0.
         end do

      else

         do k=2,nz1
            dzmina = minval(zw(k)+ah(k)*hx(k,:)-zw(k-1)-ah(k-1)*hx(k-1,:))
            call iamas_log_write('$i $r', intArgs=(/k/), realArgs=(/dzmina/(zw(k)-zw(k-1))/))
         end do

      end if

      deallocate(hs )
      deallocate(hs1)

!     Height of coordinate levels (calculation of zgrid)

      do iCell=1,nCells
         do k=1,nz
            zgrid(k,iCell) = zw(k) + ah(k)*hx(k,iCell)
         end do
         do k=1,nz1
            zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
         end do
      end do

      do i=1, nEdges
         iCell1 = cellsOnEdge(1,i)
         iCell2 = cellsOnEdge(2,i)
         do k=1,nz1
            zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
         end do
      end do
      do i=1, nCells
         do k=1,nz1
           ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
           dss(k,i) = 0.
           ztemp = zgrid(k,i)
           if (ztemp.gt.zd+.1)  then
               dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
           end if
         end do
      end do


      ! For z-metric term in omega equation
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! Avoid referencing the garbage cell for exterior edges
         if (cell1 == nCells+1) then
            cell1 = cell2
         end if
         if (cell2 == nCells+1) then
            cell2 = cell1
         end if

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                        d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                        d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do

                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  else
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

            end do

         end if
      end do

      call iamas_log_write(' grid metrics setup complete ')

      end if    ! config_vertical_grid


      if (config_met_interp) then

      !ldf (2011-11-19): added initialization of the sea-surface temperature, seaice fraction, and
      !seaice flag:
       sst = 0.0
       xice = 0.0
       seaice = 0.0
      !ldf end.

      !
      ! First, try to locate the LANDSEA field for use as an interpolation mask
      !
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error opening initial meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                   messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)
         if (trim(field % field) == 'LANDSEA') then

            allocate(maskslab(-2:field % nx+3, field % ny))
            maskslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            maskslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            maskslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            maskslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            maskslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            maskslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            maskslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
            call iamas_log_write('minval, maxval of LANDSEA = $r $r', realArgs=(/minval(maskslab), maxval(maskslab)/))

         end if

         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()

      if (.not. allocated(maskslab)) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('LANDSEA field not found in meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                   messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      edge_mask(:) = 1


      !
      ! Horizontally interpolate meteorological data
      !
      allocate(vert_level(config_nfglevels))
      vert_level(:) = -1.0
      allocate(found_u(config_nfglevels), found_v(config_nfglevels), found_t(config_nfglevels), &
               found_moisture(config_nfglevels), found_z(config_nfglevels))
      found_u(:) = .false.
      found_v(:) = .false.
      found_t(:) = .false.
      found_moisture(:) = .false.
      found_z(:) = .false.
      found_soilhgt = .false.

      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error opening initial meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                   messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      call iamas_hash_init(level_hash)
      too_many_fg_levs = .false.

      call read_next_met_field(field, istatus)

      do while (istatus == 0)

         interp_list(1) = FOUR_POINT
         interp_list(2) = SEARCH
         interp_list(3) = 0

         maskval = -1.0
         masked = -1
         fillval = 0.0
         msgval = -1.e30

         mask_array => landmask

         if (trim(field % field) == 'UU' .or. &
             trim(field % field) == 'VV' .or. &
             trim(field % field) == 'TT' .or. &
             trim(field % field) == 'RH' .or. &
             trim(field % field) == 'SPECHUMD' .or. &
             trim(field % field) == 'GHT' .or. &
             trim(field % field) == 'PMSL' .or. &
             trim(field % field) == 'PSFC' .or. &
             trim(field % field) == 'SOILHGT' .or. &
             trim(field % field) == 'SM000010' .or. &
             trim(field % field) == 'SM010040' .or. &
             trim(field % field) == 'SM040100' .or. &
             trim(field % field) == 'SM100200' .or. &
             trim(field % field) == 'SM010200' .or. &
             trim(field % field) == 'SM000007' .or. &
             trim(field % field) == 'SM007028' .or. &
             trim(field % field) == 'SM028100' .or. &
             trim(field % field) == 'SM100255' .or. &
             trim(field % field) == 'SM100289' .or. &
             trim(field % field) == 'ST000010' .or. &
             trim(field % field) == 'ST010040' .or. &
             trim(field % field) == 'ST040100' .or. &
             trim(field % field) == 'ST100200' .or. &
             trim(field % field) == 'ST010200' .or. &
             trim(field % field) == 'ST000007' .or. &
             trim(field % field) == 'ST007028' .or. &
             trim(field % field) == 'ST028100' .or. &
             trim(field % field) == 'ST100255' .or. &
             trim(field % field) == 'ST100289' .or. &
             trim(field % field) == 'PRES' .or. &
             trim(field % field) == 'PRESSURE' .or. &
             trim(field % field) == 'SNOW' .or. &
             trim(field % field) == 'SEAICE' .or. &
             trim(field % field) == 'SKINTEMP') then

            if (trim(field % field) == 'SM000010' .or. &
                trim(field % field) == 'SM010040' .or. &
                trim(field % field) == 'SM040100' .or. &
                trim(field % field) == 'SM100200' .or. &
                trim(field % field) == 'SM010200' .or. &
                trim(field % field) == 'SM000007' .or. &
                trim(field % field) == 'SM007028' .or. &
                trim(field % field) == 'SM028100' .or. &
                trim(field % field) == 'SM100255' .or. &
                trim(field % field) == 'SM100289' .or. &
                trim(field % field) == 'ST000010' .or. &
                trim(field % field) == 'ST010040' .or. &
                trim(field % field) == 'ST040100' .or. &
                trim(field % field) == 'ST100200' .or. &
                trim(field % field) == 'ST010200' .or. &
                trim(field % field) == 'ST000007' .or. &
                trim(field % field) == 'ST007028' .or. &
                trim(field % field) == 'ST028100' .or. &
                trim(field % field) == 'ST100255' .or. &
                trim(field % field) == 'ST100289' .or. &
                trim(field % field) == 'SNOW' .or. &
                trim(field % field) == 'SEAICE' .or. &
                trim(field % field) == 'SKINTEMP') then
               k = 1
            else if (trim(field % field) /= 'PMSL' .and. &
                     trim(field % field) /= 'PSFC' .and. &
                     trim(field % field) /= 'SOILHGT') then

               ! Since the hash table can only store integers, transfer the bit pattern from
               ! the real-valued xlvl into an integer; that the result is not an integer version
               ! of the level is not important, since we only want to test uniqueness of levels
               level_value = transfer(field % xlvl, level_value)
               if (.not. iamas_hash_search(level_hash, level_value)) then
                  call iamas_hash_insert(level_hash, level_value)
                  if (iamas_hash_size(level_hash) > config_nfglevels) then
                     too_many_fg_levs = .true.
                  end if
               end if

               !
               ! In case we have more than config_nfglevels levels, just keep cycling through
               ! the remaining fields in the intermediate file for the purpose of counting how
               ! many unique levels are found using the code above
               !
               if (too_many_fg_levs) then
                  call read_next_met_field(field, istatus)
                  cycle
               end if

               do k=1,config_nfglevels
                  if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
               end do
               if (vert_level(k) == -1.0) vert_level(k) = field % xlvl
            else
               k = 1
            end if

            if (trim(field % field) == 'UU') then
               found_u(k) = .true.
            else if (trim(field % field) == 'VV') then
               found_v(k) = .true.
            else if (trim(field % field) == 'TT') then
               found_t(k) = .true.
            else if (trim(field % field) == 'RH' .or. trim(field % field) == 'SPECHUMD') then
               found_moisture(k) = .true.
            else if (trim(field % field) == 'GHT') then
               found_z(k) = .true.
            else if (trim(field % field) == 'SOILHGT') then
               found_soilhgt = .true.
            end if

            !
            ! Set up projection
            !
            call map_init(proj)

            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat,RKIND), &
                            loninc = real(field % deltalon,RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
!                            nxmax = nint(360.0 / field % deltalon), &
            end if


            !
            ! Horizontally interpolate the field at level k
            !
            if (trim(field % field) == 'UU') then
               call iamas_log_write('Interpolating U at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))

               ! For U10, interpolate to cell centers
               if (vert_level(k) == 200100.0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call iamas_pool_get_array(fg, 'u10', destField1d)
                  ndims = 1

               ! otherwise to edges
               else
                  mask_array => edge_mask

                  nInterpPoints = nEdges
                  latPoints => latEdge
                  lonPoints => lonEdge
                  call iamas_pool_get_array(fg, 'u', destField2d)
                  ndims = 2
               end if

            else if (trim(field % field) == 'VV') then
               call iamas_log_write('Interpolating V at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))

               ! For V10, interpolate to cell centers
               if (vert_level(k) == 200100.0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call iamas_pool_get_array(fg, 'v10', destField1d)
                  ndims = 1

               ! otherwise to edges
               else
                  mask_array => edge_mask

                  nInterpPoints = nEdges
                  latPoints => latEdge
                  lonPoints => lonEdge
                  call iamas_pool_get_array(fg, 'v', destField2d)
                  ndims = 2
               end if

            else if (trim(field % field) == 'TT') then
               call iamas_log_write('Interpolating TT at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 't', destField2d)
               ndims = 2
            else if (trim(field % field) == 'RH') then
               call iamas_log_write('Interpolating RH at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'rh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'SPECHUMD') then
               call iamas_log_write('Interpolating SPECHUMD at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'GHT') then
               call iamas_log_write('Interpolating GHT at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'z', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRES') then
               call iamas_log_write('Interpolating PRES at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRESSURE') then
               call iamas_log_write('Interpolating PRESSURE at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PMSL') then
               call iamas_log_write('Interpolating PMSL')
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'pmsl', destField1d)
               ndims = 1
            else if (trim(field % field) == 'PSFC') then
               call iamas_log_write('Interpolating PSFC')
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'psfc', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SOILHGT') then
               call iamas_log_write('Interpolating SOILHGT')
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'soilz', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SM000010') then
               call iamas_log_write('Interpolating SM000010')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 10.
               zs_fg(k,:) = 10.
            else if (trim(field % field) == 'SM010200') then
               call iamas_log_write('Interpolating SM100200')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 200.-10.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'SM010040') then
               call iamas_log_write('Interpolating SM010040')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 40.-10.
               zs_fg(k,:) = 40.
            else if (trim(field % field) == 'SM040100') then
               call iamas_log_write('Interpolating SM040100')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-40.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'SM100200') then
               call iamas_log_write('Interpolating SM100200')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 200.-100.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'SM000007') then
               call iamas_log_write('Interpolating SM000007')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 7.
               zs_fg(k,:) = 7.
            else if (trim(field % field) == 'SM007028') then
               call iamas_log_write('Interpolating SM007028')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 28.-7.
               zs_fg(k,:) = 28.
            else if (trim(field % field) == 'SM028100') then
               call iamas_log_write('Interpolating SM028100')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-28.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'SM100255') then
               call iamas_log_write('Interpolating SM100255')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 255.-100.
               zs_fg(k,:) = 255.
            else if (trim(field % field) == 'SM100289') then
               call iamas_log_write('Interpolating SM100289')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 289.-100.
               zs_fg(k,:) = 289.
            else if (trim(field % field) == 'ST000010') then
               call iamas_log_write('Interpolating ST000010')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 10.
               zs_fg(k,:) = 10.
            else if (trim(field % field) == 'ST010200') then
               call iamas_log_write('Interpolating ST010200')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 200.-10.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'ST010040') then
               call iamas_log_write('Interpolating ST010040')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 40.-10.
               zs_fg(k,:) = 40.
            else if (trim(field % field) == 'ST040100') then
               call iamas_log_write('Interpolating ST040100')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-40.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'ST100200') then
               call iamas_log_write('Interpolating ST100200')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 200.-100.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'ST000007') then
               call iamas_log_write('Interpolating ST000007')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 7.
               zs_fg(k,:) = 7.
            else if (trim(field % field) == 'ST007028') then
               call iamas_log_write('Interpolating ST007028')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 28.-7.
               zs_fg(k,:) = 28.
            else if (trim(field % field) == 'ST028100') then
               call iamas_log_write('Interpolating ST028100')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-28.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'ST100255') then
               call iamas_log_write('Interpolating ST100255')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 255.-100.
               zs_fg(k,:) = 255.
            else if (trim(field % field) == 'ST100289') then
               call iamas_log_write('Interpolating ST100289')

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'st_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 289.-100.
               zs_fg(k,:) = 289.
            else if (trim(field % field) == 'SNOW') then
               call iamas_log_write('Interpolating SNOW')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0

               masked = 0
               fillval = 0.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'snow', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SEAICE') then
               call iamas_log_write('Interpolating SEAICE')

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 1.0
               masked = 1
               fillval = 0.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'xice', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SKINTEMP') then
               call iamas_log_write('Interpolating SKINTEMP')
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'skintemp', destField1d)
               ndims = 1
            end if

            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               if (mask_array(i) /= masked) then
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field%ny)
                  end if
                  if (ndims == 1) then
                     if (maskval /= -1.0) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  else if (ndims == 2) then
                     if (maskval /= -1.0) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  end if
               else
                  if (ndims == 1) then
                     destField1d(i) = fillval
                  else if (ndims == 2) then
                     destField2d(k,i) = fillval
                  end if
               end if
            end do

            !
            ! In addition to interpolating wind fields to cell edges, we should
            ! also interpolate to cell centers at the surface in order to
            ! produce U10 and V10 fields
            !
            is_sfc_field = .false.
            mask_array => landmask
            if (index(field % field, 'UU') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'u10', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'VV') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'v10', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'TT') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 't2m', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'RH') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'rh2', destField1d)
               ndims = 1
               is_sfc_field = .true.
            end if

            if (is_sfc_field) then
               do i=1,nInterpPoints
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, &
                                                   msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
               end do
            end if

            deallocate(rslab)

         end if

         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()
      level_value = iamas_hash_size(level_hash)
      call iamas_hash_destroy(level_hash)

      missing_required_fields = 0
      do k=1,config_nfglevels
         if (vert_level(k) == -1.0 .or. vert_level(k) == 200100.0) cycle
         if (.not. found_u(k) .or. .not. found_v(k) .or. .not. found_t(k) .or. &
             .not. found_moisture(k) .or. .not. found_z(k)) then
            missing_required_fields = missing_required_fields + 1
            call iamas_log_write('Incomplete required field set at meteorological level $r.', &
                                realArgs=(/vert_level(k)/), messageType=iAMAS_LOG_ERR)
         end if
      end do
      ! SOILHGT is consumed only by optional regional boundary-terrain blending.
      if (config_blend_bdy_terrain .and. .not. found_soilhgt) then
         missing_required_fields = missing_required_fields + 1
         call iamas_log_write('Missing required SOILHGT field in meteorological input.', messageType=iAMAS_LOG_ERR)
      end if
      if (missing_required_fields > 0) then
         call iamas_log_write('Refusing to interpolate incomplete meteorological input.', messageType=iAMAS_LOG_CRIT)
      end if

      if (too_many_fg_levs) then
         write(errstring,'(a,i4)') '       Please increase config_nfglevels to at least ', level_value
         call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error: The meteorological data file has more than config_nfglevels.', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(trim(errstring),                                                       messageType=iAMAS_LOG_ERR)
         call iamas_log_write('       in the namelist and re-run.',                                  messageType=iAMAS_LOG_ERR)
         call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_CRIT)
      end if


      !
      ! Check how many distinct levels we actually found in the meteorological data
      !
      do k=1,config_nfglevels
         if (vert_level(k) == -1.0) exit
      end do
      nfglevels_actual = k-1
      call iamas_log_write('*************************************************')
      call iamas_log_write('Found $i levels in the first-guess data', intArgs=(/nfglevels_actual/))
      call iamas_log_write('*************************************************')


      !
      ! Extract surface fields from first-guess
      !
!      call iamas_pool_get_array(fg, 'u10', u10)
!      call iamas_pool_get_array(fg, 'v10', v10)
      call iamas_pool_get_array(fg, 'q2', q2)
      call iamas_pool_get_array(fg, 'rh2', rh2)
      call iamas_pool_get_array(fg, 't2m', t2m)
      call iamas_pool_get_array(fg, 'u', u_fg)
      call iamas_pool_get_array(fg, 'v', v_fg)
      call iamas_pool_get_array(fg, 't', t_fg)
      call iamas_pool_get_array(fg, 'rh', rh_fg)
!      u10(:) = 0.0
!      v10(:) = 0.0
      q2(:) = 0.0
      rh2(:) = 0.0
      t2m(:) = 0.0

      do k=1,config_nfglevels
         if (vert_level(k) == 200100.0) then
!            u10(:) = u_fg(k,:)
!            v10(:) = v_fg(k,:)
            t2m(:) = t_fg(k,:)
            rh2(:) = rh_fg(k,:)
            do iCell = 1, nCells
               es = svp1 * 10.0_RKIND * exp(svp2 * (t2m(iCell)-svpt0) / (t2m(iCell)-svp3))
               es = min(es, 0.99_RKIND * 0.01_RKIND * psfc(iCell))
               rs = 0.622_RKIND * es * 100.0_RKIND / (psfc(iCell) - es * 100.0_RKIND)
               q2(iCell) = 0.01_RKIND * rs * rh2(iCell)
             end do
         end if
      end do



      !
      ! For isobaric data, fill in the 3-d pressure field; otherwise, ensure
      ! that the surface pressure and height fields are filled in
      !
      if (minval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0 .and. &
          maxval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0) then
         call iamas_log_write('Setting pressure field for isobaric data')
         do k=1,config_nfglevels
            if (vert_level(k) /= 200100.0) then
               p_fg(k,:) = vert_level(k)
            else
               p_fg(k,:) = psfc(:)
            end if
         end do
      else
         call iamas_pool_get_array(fg, 'z', z_fg)
         call iamas_pool_get_array(fg, 'soilz', soilz)
         call iamas_log_write('Assuming model-level input data')
         do k=1,config_nfglevels
            if (vert_level(k) == 200100.0) then
               p_fg(k,:) = psfc(:)
               z_fg(k,:) = soilz(:)
            end if
         end do
      end if

      ! Set SST based on SKINTEMP field if it wasn't found in input data
      if (minval(sst) == 0.0 .and. maxval(sst) == 0.0) then
         call iamas_log_write('Setting SST from SKINTEMP')
         !where (landmask == 0) sst = skintemp
         sst = skintemp
      end if

      ! Set SNOWC (snow-cover flag) based on SNOW
      snowc(:) = 0.0
      where (snow > 0.0) snowc = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MGD CHECK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iCell=1,nCells
   if (landmask(iCell) == 1) then

      do k = 1, config_nfgsoillevels
         if (st_fg(k,iCell) <= 0.0) call iamas_log_write('Bad st_fg $i $i', intArgs=(/k, iCell/))
      end do

      do k = 1, config_nfgsoillevels
         if (sm_fg(k,iCell) <= 0.0) then
            call iamas_log_write('Bad sm_fg $r $i $i', intArgs=(/k, iCell/), realArgs=(/sm_fg(k,iCell)/))
            sm_fg(k,iCell) = 0.001
         end if
      end do
      !LDF end.

   end if
end do
call iamas_log_write('Done with soil consistency check')

! jwfeng from AEMOL USTC: 05/2019
! Get EROD and DMS_0 from a separate file
! This separate file is generated by a Fortran file.
      erod(:,:) = 0._RKIND
      dms_0(:) = 0._RKIND
!     call read_met_init('metdata/EROD_DMS0', .true., config_start_time(1:13), istatus)
      call read_met_init(trim(config_source_file), .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('****************'//trim(config_source_file)//' is not found! ****************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(         '********************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            !
            ! Set up projection
            !
            call map_init(proj)
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat, RKIND), &
                            loninc = real(field % deltalon, RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat, RKIND), &
                            lon1 = real(field % startlon, RKIND))

            else if (field % iproj == PROJ_PS) then
               call map_set(PROJ_PS, proj, &
                            dx = real(field % dx,RKIND), &
                            truelat1 = real(field % truelat1,RKIND), &
                            stdlon = real(field % xlonc,RKIND), &
                            knowni = real(field % nx / 2.0,RKIND), &
                            knownj = real(field % ny / 2.0,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else
               call iamas_log_write('Unsupported source-file projection, iproj = $i', &
                                   intArgs=(/ field % iproj /), messageType=iAMAS_LOG_ERR)
               call iamas_dmpar_global_abort('Unsupported source-file projection in config_source_file')
 !                          nxmax = nint(360.0 / field % deltalon), &
            end if

         if (trim(config_source_opt) == 'erod') then
   
            if (trim(field % field) == 'EROD01') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'erod', destField2d)
               masked = 0 !landmask=0 reprsent ocean. 1 present land.
               ndims = 2
               k     = 1
            else if (trim(field % field) == 'EROD02') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'erod', destField2d)
               masked = 0
               ndims = 2
               k     = 2
            else if (trim(field % field) == 'EROD03') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'erod', destField2d)
               masked = 0
               ndims = 2
               k     = 3
            else if (trim(field % field) == 'DMS_0') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'dms_0', destField1d)
               masked = 1
               ndims = 1
            else
               deallocate(field % slab)
               call read_next_met_field(field, istatus)
               cycle
            end if

         else if (trim(config_source_opt) == 'ssm') then

            if (trim(field % field) == 'SSM') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'erod', destField2d)
               masked = 0
               ndims = 2
               k = 1

            else if (trim(field % field) == 'DMS_0') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'dms_0', destField1d)
               masked = 1
               ndims = 1

            else
               deallocate(field % slab)
               call read_next_met_field(field, istatus)
               cycle
            end if
         else
            call iamas_log_write('Unsupported config_source_opt = '//trim(config_source_opt), &
                                messageType=iAMAS_LOG_CRIT)
         end if

            interp_list(1) = FOUR_POINT
            interp_list(2) = W_AVERAGE4
            interp_list(3) = 0

!           maskval = 1.0
!           masked = 0
!           fillval = 0.0
!           mask_array => landmask
            fillval = 0.0
            msgval = 0.0
            mask_array => landmask


            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               
               if (mask_array(i) /= masked) then
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field%ny)
                  end if
                  if (ndims == 1) then
                     if (maskval /= -1.0) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  else if (ndims == 2) then
                     if (maskval /= -1.0) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  end if
               else

                  if (ndims == 1) then
                     destField1d(i) = fillval
                  else if (ndims == 2) then
                     destField2d(k,i) = fillval
                  end if
               end if
            end do
            deallocate(rslab)
            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()

      if (trim(config_source_opt) == 'ssm') then
         where (erod(1,:) < 0._RKIND) erod(1,:) = 0._RKIND
         where (erod(1,:) > 1._RKIND) erod(1,:) = 1._RKIND
         erod(2,:) = 0._RKIND
         erod(3,:) = 0._RKIND
      end if

! added by qkyang to support anthropogenic heating emissions>>>>>>>>>>>>>>
! Get AH_emis 
if (is_ah_emis) then
      AH_emis(:) = 0. !Default value
      call read_met_init(trim(config_AH_emis_file), .true., config_start_time(1:13), istatus)!config_start_time support emissions time

      if (istatus /= 0) then
         call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('****************'//trim(config_AH_emis_file)//' is not found! ****************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(         '********************************************************', messageType=iAMAS_LOG_CRIT)
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            !
            ! Set up projection
            !
            call map_init(proj)
            if (field % iproj == PROJ_PS) then
               call map_set(PROJ_PS, proj, &
                            dx = real(field % dx,RKIND), &
                            truelat1 = real(field % truelat1,RKIND), &
                            stdlon = real(field % xlonc,RKIND), &
                            knowni = real(field % nx / 2.0,RKIND), &
                            knownj = real(field % ny / 2.0,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
 !                          nxmax = nint(360.0 / field % deltalon), &
            end if
            if (trim(field % field) == 'AH_emis') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'AH_emis', destField1d)
               masked = 0 
               ndims = 1
            end if

            interp_list(1) = FOUR_POINT
            ! interp_list(2) = W_AVERAGE4
            interp_list(2) = W_AVERAGE16
            interp_list(3) = 0

            fillval = 0.0 !missing value
            msgval = 0.0
            mask_array => landmask !landmask=0 is water

            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               if (mask_array(i) /= masked) then
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field%ny)
                  end if
                  if (ndims == 1) then
                     if (maskval /= -1.0) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  else if (ndims == 2) then
                     if (maskval /= -1.0) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  end if
               else
                  if (ndims == 1) then
                     destField1d(i) = fillval
                  else if (ndims == 2) then
                     destField2d(k,i) = fillval
                  end if
               end if
            end do
            deallocate(rslab)
            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()
end if
! added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<<

      !
      ! Get SEAICE from a separate file
      !
      call read_met_init('SEAICE_FRACTIONAL', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('SEAICE_FRACTIONAL file not found...')
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (trim(field % field) == 'SEAICE') then

               call iamas_log_write('PROCESSING SEAICE')

               !
               ! Set up projection
               !
               call map_init(proj)

               if (field % iproj == PROJ_PS) then
                  call map_set(PROJ_PS, proj, &
                               dx = real(field % dx,RKIND), &
                               truelat1 = real(field % truelat1,RKIND), &
                               stdlon = real(field % xlonc,RKIND), &
                               knowni = real(field % nx / 2.0,RKIND), &
                               knownj = real(field % ny / 2.0,RKIND), &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               end if

               if (trim(field % field) == 'SEAICE') then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call iamas_pool_get_array(fg, 'xice', destField1d)
                  ndims = 1
               end if

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0

               masked = 1
               fillval = 0.0
               msgval = 1.01
               mask_array => landmask


               allocate(rslab(field % nx, field % ny))
               rslab(:,:) = field % slab(:,:)
               do i=1,nInterpPoints
                  if (mask_array(i) /= masked) then
                     lat = latPoints(i)*DEG_PER_RAD
                     lon = lonPoints(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     else if (x >= real(field%nx)+0.5) then
                        lon = lon - 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (y < 0.5) then
                        y = 1.0
                     else if (y >= real(field%ny)+0.5) then
                        y = real(field%ny)
                     end if
                     if (ndims == 1) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField1d(i) == msgval) destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField2d(k,i) == msgval) destField2d(k,i) = fillval
                     end if
                  else
                     if (ndims == 1) then
                        destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = fillval
                     end if
                  end if
               end do
               deallocate(rslab)

            end if

            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()


      !
      ! Get OMLD climatology from a separate file
      !
      call read_met_init('OMLD', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('OMLD file not found...')
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (index(field % field, 'OMLD') /= 0) then

               call iamas_log_write('PROCESSING OMLD')

               !
               ! Set up projection
               !
               call map_init(proj)

               if (field % iproj == PROJ_LATLON) then
                  call map_set(PROJ_LATLON, proj, &
                               latinc = real(field % deltalat,RKIND), &
                               loninc = real(field % deltalon,RKIND), &
                               knowni = 1.0_RKIND, &
                               knownj = 1.0_RKIND, &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               else
                  call iamas_log_write('We were expecting OMLD field to be on a lat-lon projection...', messageType=iAMAS_LOG_CRIT)
               end if

               if (index(field % field, 'OMLD') /= 0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call iamas_pool_get_array(state, 'h_oml_initial', destField1d)
                  ndims = 1
               end if

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               masked = 1
               fillval = 0.0
               msgval = 0.0
               mask_array => landmask


               allocate(rslab(field % nx, field % ny))
               rslab(:,:) = field % slab(:,:)
               do i=1,nInterpPoints
                  if (mask_array(i) /= masked) then
                     lat = latPoints(i)*DEG_PER_RAD
                     lon = lonPoints(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     else if (x >= real(field%nx)+0.5) then
                        lon = lon - 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (y < 0.5) then
                        y = 1.0
                     else if (y >= real(field%ny)+0.5) then
                        y = real(field%ny)
                     end if
                     if (ndims == 1) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField1d(i) == msgval) destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField2d(k,i) == msgval) destField2d(k,i) = fillval
                     end if
                  else
                     if (ndims == 1) then
                        destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = fillval
                     end if
                  end if
               end do
               deallocate(rslab)

            end if

            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()


      if (allocated(maskslab)) deallocate(maskslab)

      ! Freeze really cold ocean
      where (sst < 271.0 .and. landmask == 0) xice = 1.0

      ! Limit XICE to values between 0 and 1. Although the input meteorological field is between 0.
      ! and 1., interpolation to the iAMAS grid can yield values of XiCE less than 0. and greater
      ! than 1.:
      where (xice < 0._RKIND) xice = 0._RKIND
      where (xice > 1._RKIND) xice = 1._RKIND

      ! Set SEAICE (0/1 flag) based on XICE (fractional ice coverage)
      seaice(:) = 0.0
      where (xice >= 0.5) seaice = 1.0


      !
      ! Compute normal wind component and store in fg % u
      !
      do iEdge=1,nEdges
         do k=1,nfglevels_actual
            u_fg(k,iEdge) = cos(angleEdge(iEdge)) * u_fg(k,iEdge) &
                          + sin(angleEdge(iEdge)) * v_fg(k,iEdge)
         end do
      end do

      !
      ! Vertically interpolate meteorological data
      !
      allocate(sorted_arr(2,nfglevels_actual))

      do iCell=1,nCells

         ! T
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 9999999.0 ! jwfeng: avoid err when ztop is larger than 99999.0
            sorted_arr(2,k) = t_fg(k,iCell)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            ! jwfeng: Use temperature instead of theta_m
            temperature(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                         sorted_arr(:,1:nfglevels_actual-1), order=1, &
                                         extrap=extrap_airtemp, ierr=istatus)
            if (istatus /= 0) then
               write(errstring,'(a,i4,a,i10)') 'Error in interpolation of t(k,iCell) for k=', k, ', iCell=', iCell
               call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_ERR)
               call iamas_log_write(trim(errstring),                                                     messageType=iAMAS_LOG_ERR)
               call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_CRIT)
            end if
         end do


         ! RH
         sorted_arr(:,:) = -999.0
         relhum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
            sorted_arr(2,k) = rh_fg(k,iCell)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            relhum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                       sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
         end do


         ! SPECHUM: if first-guess values are negative, set those values to zero before
         ! vertical interpolation.
         sorted_arr(:,:) = -999.0
         spechum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
            sorted_arr(2,k) = max(0._RKIND,sh_fg(k,iCell))
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            spechum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                        sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
            if (target_z < z_fg(1,iCell) .and. k < nVertLevels) spechum(k,iCell) = spechum(k+1,iCell)
         end do


         ! GHT
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
            sorted_arr(2,k) = z_fg(k,iCell)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            gfs_z(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                      sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
         end do


         ! PRESSURE
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) then
               sorted_arr(1,k) = 99999999.0
               sfc_k = k
            end if
            sorted_arr(2,k) = log(p_fg(k,iCell))
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            pressure(k,iCell) = exp(vertical_interp(target_z, nfglevels_actual-1, &
                                    sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1))
         end do

      end do


      do iEdge=1,nEdges

         ! U
         sorted_arr(:,:) = -999.0
         do k=1,nfglevels_actual
            sorted_arr(1,k) = 0.5 * (z_fg(k,cellsOnEdge(1,iEdge)) + z_fg(k,cellsOnEdge(2,iEdge)))
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 0.5 * (fg % soilz % array(cellsOnEdge(1,iEdge)) + fg % soilz % array(cellsOnEdge(2,iEdge)))
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
            sorted_arr(2,k) = u_fg(k,iEdge)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k=1,nVertLevels
            target_z = 0.25 * (zgrid(k,cellsOnEdge(1,iEdge)) + zgrid(k+1,cellsOnEdge(1,iEdge)) + zgrid(k,cellsOnEdge(2,iEdge)) + zgrid(k+1,cellsOnEdge(2,iEdge)))
!           u(k,iEdge) = vertical_interp(target_z, nfglevels_actual, sorted_arr, order=1, extrap=0)
            u(k,iEdge) = vertical_interp(target_z, nfglevels_actual-1, sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
         end do

      end do


      !
      ! Reconstruct zonal and meridional winds for diagnostic puposes:
      !
      call iamas_rbf_interp_initialize(mesh)
      call iamas_init_reconstruct(mesh)
      call iamas_reconstruct(mesh, u,                 &
                            uReconstructX,           &
                            uReconstructY,           &
                            uReconstructZ,           &
                            uReconstructZonal,       &
                            uReconstructMeridional   &
                           )


      !
      ! Adjust surface pressure for difference in topography
      !
      do sfc_k=1,nfglevels_actual
         if (vert_level(sfc_k) == 200100.) exit
      end do
      do iCell=1,nCells

         ! We need to extrapolate
            sorted_arr(:,:) = -999.0
            do k=1,nfglevels_actual
               sorted_arr(1,k) = z_fg(k,iCell)
               if (vert_level(k) == 200100.0) then
!NOSFC                  sorted_arr(1,k) = fg % soilz % array(iCell)
                  sorted_arr(1,k) = 99999999.0
               end if
               sorted_arr(2,k) = log(p_fg(k,iCell))
            end do
            call iamas_quicksort(nfglevels_actual, sorted_arr)
            target_z = zgrid(1,iCell)
            psfc(iCell) = exp(vertical_interp(target_z, nfglevels_actual, sorted_arr, order=1, extrap=1))

      end do

      deallocate(sorted_arr)


      ! Diagnose the water vapor mixing ratios:
      global_sh_min = 0._RKIND
      global_sh_max = 0._RKIND
      if(config_use_spechumd) then
         sh_min = minval(spechum(:,1:nCellsSolve))
         sh_max = maxval(spechum(:,1:nCellsSolve))
         call iamas_dmpar_min_real(dminfo,sh_min,global_sh_min)
         call iamas_dmpar_max_real(dminfo,sh_max,global_sh_max)
      endif
      call iamas_log_write('')
      call iamas_log_write('--- global_sh_min = $r', realArgs=(/global_sh_min/))
      call iamas_log_write('--- global_sh_max = $r', realArgs=(/global_sh_max/))
      call iamas_log_write('')

      call iamas_log_write('--- config_use_spechumd = $l', logicArgs=(/config_use_spechumd/))
      if(.not. config_use_spechumd .or. (global_sh_min==0._RKIND .and. global_sh_max==0._RKIND)) then
         !--- calculate the saturation mixing ratio and interpolated first-guess relative humidity:
         if (config_use_spechumd) then
            call iamas_log_write('config_use_spechumd=T, but specific humidity was not found in '//trim(config_met_prefix)//':'//config_start_time(1:13), messageType=iAMAS_LOG_WARN)
         end if
         call iamas_log_write(' *** initializing water vapor mixing ratio using first-guess relative humidity')
         call iamas_log_write('')

         do k = 1, nVertLevels
            do iCell = 1, nCells
!              es = svp1*10.*exp(svp2*(t(k,iCell)-svpt0)/(t(k,iCell)-svp3))
!              es = min(es,0.99*0.01*pressure(k,iCell))
!              rs = 0.622*es*100. / (pressure(k,iCell)-es*100.)

               !
               ! Note: the RH field provided by ungrib should always be with respect to liquid water,
               !       hence, we can always call rslf; see the routine fix_gfs_rh in WPS/ungrib/src/rrpr.F .
               !
               rs = rslf(pressure(k,iCell),temperature(k,iCell)) ! Use temperature instead of theta_m, jwfeng
               scalars(index_qv,k,iCell) = 0.01_RKIND*rs*relhum(k,iCell)
               qv(k,iCell) = scalars(index_qv,k,iCell)
             enddo
         enddo
      else
         !--- use the interpolated first-guess specific humidity:
         call iamas_log_write(' *** initializing water vapor mixing ratio using first-guess specific humidity')
         call iamas_log_write('')
         do k = 1, nVertLevels
            do iCell = 1, nCells
               scalars(index_qv,k,iCell) = spechum(k,iCell)/(1._RKIND-spechum(k,iCell))
               qv(k,iCell) = scalars(index_qv,k,iCell)
            enddo
         enddo
      endif

      !
      ! After RH has been used to compute qv (unless config_use_spechumd = T and a valid spechum field
      ! is available), modify the RH field to be with respect to ice for temperatures below freezing.
      ! NB: Here we pass in 1:nCells explicitly, since computations involving the "garbage cell" could
      ! trigger FPEs.
      !
      call convert_relhum_wrt_ice(temperature(:,1:nCells), relhum(:,1:nCells)) ! Use temperature instead of theta_m, jwfeng

      !
      ! Diagnose fields needed in initial conditions file (u, w, rho, theta)
      ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
      !
      do iCell=1,nCells

         ! Q2
         es = 6.112 * exp((17.27*(t2m(iCell) - 273.16))/(t2m(iCell) - 35.86))
         rs = 0.622 * es * 100. / (psfc(iCell) - es * 100.)
         q2(iCell) = 0.01 * rs * rh2(iCell)


         do k=1,nVertLevels
            ! PI
            exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)

            ! THETA - can compute this using PI instead
            ! t(k,iCell) = t(k,iCell) / exner(k,iCell)
            ! jwfeng: theta_m of RHS is temperature, will be converted into theta_m here
            theta_m(k,iCell) = temperature(k,iCell) / exner(k,iCell) * (1.0 + 1.61*qv(k,iCell))

            ! RHO_ZZ, couple with vertical metric here, jwfeng
            rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
         end do
      end do


      !
      ! Calculation of the initial precipitable water:
      !
      do iCell = 1,nCells
         precipw(iCell) = 0.0
         do k = 1,nVertLevels
            precipw(iCell) = precipw(iCell) + rho_zz(k,iCell)*zz(k,iCell)*qv(k,iCell)*(zgrid(k+1,iCell)-zgrid(k,iCell))
         end do
      end do

      !
      ! Reference state based on a dry isothermal atmosphere
      !
      do iCell=1,nCells
         do k=1,nz1
            ztemp    = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
            pressure_b(k,iCell) = p0*exp(-gravity*ztemp/(rgas*t0b))     ! pressure_base
            exner_b(k,iCell) = (pressure_b(k,iCell)/p0)**(rgas/cp)      ! exner_base
            rho_b(k,iCell) = pressure_b(k,iCell)/(rgas*t0b)/zz(k,iCell) ! rho_base
            theta_b (k,iCell) = t0b/exner_b(k,iCell)                    ! theta_base
            rtheta_b(k,iCell) = rho_b(k,iCell)*theta_b(k,iCell)         ! rtheta_base
            exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)    ! exner
            pressure_p(k,iCell) = pressure(k,iCell)-pressure_b(k,iCell) ! pressure_p
            rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)           ! rho_p

            ! jwfeng: Use phi = ln(p/p0) as log pressure.
            phi_b(k,iCell) = -gravity*ztemp/(rgas*t0b)
            phi_total(k,iCell) = log(pressure(k,iCell)/p0)
         end do
      end do

      if (config_log_pressure) then
         ! jwfeng: Use phi instead of pressure_p to balance the initial thermodynamic fields.
         ! d(phi)/dz = -g*rho_zz*(1+qv)/pressure = -g*(1+qv)/(rgas*T_m)/zz, T_m = (1+1.61*qv)*T = theta_m * exner
         ! To be simplified, use phi without perturbation. Klemp et al., (2021) Adapting the iAMAS ...
         do iCell=1,nCells
         do k=2,nVertLevels
            it = 0
            p_check = 2.0 * 0.000001
            do while ( (it < 100) .and. (p_check > 0.000001) )
               p_check = phi_total(k,iCell)

               phi_total(k,iCell) = phi_total(k-1,iCell) - gravity*dzu(k)*( &
                  fzm(k)*( (1.0+qv(k  ,iCell))/(rgas*theta_m(k  ,iCell)*exner(k  ,iCell))/zz(k  ,iCell) ) + &
                  fzp(k)*( (1.0+qv(k-1,iCell))/(rgas*theta_m(k-1,iCell)*exner(k-1,iCell))/zz(k-1,iCell) )   &
               )
               pressure(k,iCell) = p0 * exp(phi_total(k,iCell))
               exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)
               ! jwfeng: pressure is adjusted above, so we need to re-compute theta here 
               theta_m(k,iCell) = temperature(k,iCell) / exner(k,iCell) * (1.0 + 1.61*qv(k,iCell))
               rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
               rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)

               p_check = abs(p_check - phi_total(k,iCell))

               pressure_p(k,iCell) = pressure(k,iCell) - pressure_b(k,iCell)

               it = it + 1
            end do
         end do
         end do
      else
         do iCell=1,nCells
         do k=2,nVertLevels
            it = 0
            p_check = 2.0 * 0.0001
            do while ( (it < 30) .and. (p_check > 0.0001) )

               p_check = pressure_p(k,iCell)
               !  WCS 20130821 - iAMAS hydrostatic relation
               pressure_p(k,iCell) = pressure_p(k-1,iCell) - &
                                        (fzm(k)*rho_p(k,iCell) + fzp(k)*rho_p(k-1,iCell))*gravity*dzu(k) &
                                       -(fzm(k)*rho_zz(k  ,iCell)*qv(k  ,iCell) &
                                       + fzp(k)*rho_zz(k-1,iCell)*qv(k-1,iCell))*gravity*dzu(k)
               pressure(k,iCell) = pressure_p(k,iCell) + pressure_b(k,iCell)
               exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)
               !  WCS 20130821 - couple with vertical metric
               rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
               rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)

               p_check = abs(p_check - pressure_p(k,iCell))

               it = it + 1
            end do
         end do
         end do
      end if ! jwfeng: config_log_pressure

      do iEdge=1,nEdges
         do k=1,nVertLevels
            ru(k,iEdge) = u(k,iEdge) * 0.5*(rho_zz(k,cellsOnEdge(1,iEdge)) + rho_zz(k,cellsOnEdge(2,iEdge)))
         end do
      end do


      rw(:,:) = 0.0

      do iCell=1,nCellsSolve

         do i=1,nEdgesOnCell(iCell)
            iEdge=edgesOnCell(i,iCell)

            do k = 2, nVertLevels
               flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
               if (iCell == cellsOnEdge(1,iEdge)) then
                  rw(k,iCell) = rw(k,iCell) - (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb(k,1,iEdge)*flux
               else
                  rw(k,iCell) = rw(k,iCell) + (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb(k,2,iEdge)*flux
               end if

               if (config_theta_adv_order ==3) then 
                  if (iCell == cellsOnEdge(1,iEdge)) then
                     rw(k,iCell) = rw(k,iCell)    &
                                  + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                    (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb3(k,1,iEdge)*flux
                  else
                     rw(k,iCell) = rw(k,iCell)    &
                                  - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                    (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb3(k,2,iEdge)*flux
                  end if
               end if

            end do

         end do

      end do


      ! Compute w from rho_zz and rw
      do iCell=1,nCellsSolve
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do

      deallocate(vert_level)
      deallocate(found_u, found_v, found_t, found_moisture, found_z)


      ! Calculate surface pressure (This is an ad-hoc calculation. The actual surface pressure is actually re-calculated at
      !the top of the subroutine iAMAS_to_physics in ../core_atmos_physics/iamas_atmphys_interface_nhyd.F
      do iCell=1,nCells
         surface_pressure(iCell) = 0.5*gravity/rdzw(1)                                              &
                                 * (1.25* rho_zz(1,iCell) * (1. + qv( 1, iCell))  &
                                 -  0.25* rho_zz(2,iCell) * (1. + qv( 2, iCell)))
         surface_pressure(iCell) = surface_pressure(iCell) + pressure_p(1,iCell) + pressure_b(1,iCell)
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = theta_m(k,iCell) / (1.0 + 1.61 * qv(k,iCell))
         end do
      end do


      end if    ! config_met_interp

   end subroutine init_atm_case_gfs


   !-----------------------------------------------------------------------
   !  routine init_atm_case_lbc
   !
   !> \brief Computes lbc_{rho,theta,u,w,qx} fields for lateral boundary conditions
   !> \author Michael Duda
   !> \date   22 April 2019
   !> \details
   !>  This routine is similar to the init_atm_case_gfs routine in that it reads
   !>  atmospheric fields from "intermediate" files and horizontally and vertically
   !>  interpolates them to an iAMAS mesh. However, rather than producing model
   !>  initial conditions, this routine is responsible for producing only those
   !>  fields that are needed as model lateral boundary conditions.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_case_lbc(timestamp, block, mesh, nCells, nEdges, nVertLevels, fg, state, diag, lbc_state, dims, configs)

      use iamas_dmpar, only : iamas_dmpar_min_real, iamas_dmpar_max_real
      use init_atm_read_met, only : met_data, read_met_init, read_met_close, read_next_met_field
      use init_atm_llxy, only : proj_info, map_init, map_set, latlon_to_ij, PROJ_LATLON, PROJ_GAUSS, DEG_PER_RAD
      use init_atm_hinterp, only : interp_sequence, FOUR_POINT, SIXTEEN_POINT, W_AVERAGE4, SEARCH
      use iamas_hash, only : hashtable, iamas_hash_init, iamas_hash_destroy, iamas_hash_search, iamas_hash_size, iamas_hash_insert

      implicit none

      character(len=*), intent(in) :: timestamp
      type (block_type), intent(inout), target :: block
      type (iamas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: fg
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(inout) :: lbc_state
      type (iamas_pool_type), intent(inout):: dims
      type (iamas_pool_type), intent(inout):: configs

      type (dm_info), pointer :: dminfo

      real (kind=RKIND), parameter :: t0b = 250.0

      type (met_data) :: field
      type (proj_info) :: proj

      real (kind=RKIND), dimension(:), pointer :: dzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: vert_level, latPoints, lonPoints
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zz
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, t, rt
      real (kind=RKIND), dimension(:), pointer :: destField1d
      real (kind=RKIND), dimension(:,:), pointer :: destField2d
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND) :: target_z
      integer :: iCell, iEdge, i, k, nVertLevelsP1
      integer, pointer :: nCellsSolve
      integer :: nInterpPoints, ndims

      integer :: nfglevels_actual
      integer, pointer :: index_qv

      integer, dimension(5) :: interp_list
      real (kind=RKIND) :: msgval

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: sorted_arr

      integer :: sfc_k

      integer :: it
      real (kind=RKIND) :: p_check

      integer :: istatus

      real (kind=RKIND), allocatable, dimension(:,:) :: rslab

      real (kind=RKIND) :: flux
      real (kind=RKIND) :: lat, lon, x, y

      real (kind=RKIND) :: p0

      real (kind=RKIND) :: etavs, ztemp

      real (kind=RKIND) :: rs, rcv

      real (kind=RKIND), dimension(nVertLevels + 1) :: sh

      !  calculation of the water vapor mixing ratio:
      real (kind=RKIND) :: sh_max,sh_min,global_sh_max,global_sh_min

      character (len=StrKIND), pointer :: config_met_prefix
      logical, pointer :: config_use_spechumd
      logical, pointer :: config_blend_bdy_terrain
      integer, pointer :: config_nfglevels
      integer, pointer :: config_theta_adv_order
      real (kind=RKIND), pointer :: config_coef_3rd_order

      character (len=StrKIND), pointer :: config_extrap_airtemp
      integer :: extrap_airtemp

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: angleEdge

      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: w
      real (kind=RKIND), dimension(:,:), pointer :: theta
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:), pointer :: relhum
      real (kind=RKIND), dimension(:,:), pointer :: spechum
      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: rw

      real (kind=RKIND), dimension(:,:), pointer :: u_fg
      real (kind=RKIND), dimension(:,:), pointer :: v_fg
      real (kind=RKIND), dimension(:,:), pointer :: z_fg
      real (kind=RKIND), dimension(:,:), pointer :: t_fg
      real (kind=RKIND), dimension(:,:), pointer :: rh_fg
      real (kind=RKIND), dimension(:,:), pointer :: sh_fg
      real (kind=RKIND), dimension(:,:), pointer :: p_fg
      real (kind=RKIND), dimension(:), pointer :: soilz

      type (hashtable) :: level_hash
      logical :: too_many_fg_levs
      logical, allocatable, dimension(:) :: found_u, found_v, found_t, found_moisture, found_z
      logical :: found_soilhgt
      integer :: level_value
      integer :: missing_required_fields

      character (len=StrKIND) :: errstring

      real (kind=RKIND) :: max_zgrid_local, max_zgrid_global
      integer, pointer :: index_tke_sgs
      real (kind=RKIND), pointer :: config_les_tke_min
      character (len=StrKIND), pointer :: config_les_scheme
      character (len=StrKIND), pointer :: config_les_lbc_tke_mode


      call iamas_log_write('Interpolating LBCs at time '//trim(timestamp))

      call iamas_pool_get_config(configs, 'config_met_prefix', config_met_prefix)
      call iamas_pool_get_config(configs, 'config_use_spechumd', config_use_spechumd)
      call iamas_pool_get_config(configs, 'config_blend_bdy_terrain', config_blend_bdy_terrain)
      call iamas_pool_get_config(configs, 'config_nfglevels', config_nfglevels)
      call iamas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)
      call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      call iamas_pool_get_config(configs, 'config_les_lbc_tke_mode', config_les_lbc_tke_mode)
      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)

      call iamas_pool_get_config(configs, 'config_extrap_airtemp', config_extrap_airtemp)
      if (trim(config_extrap_airtemp) == 'constant') then
         extrap_airtemp = 0
      else if (trim(config_extrap_airtemp) == 'linear') then
         extrap_airtemp = 1
      else if (trim(config_extrap_airtemp) == 'lapse-rate') then
         extrap_airtemp = 2
      else
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('* Invalid value for namelist variable config_extrap_airtemp *', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
      call iamas_log_write("Using option '" // trim(config_extrap_airtemp) // "' for vertical extrapolation of temperature")

      dminfo => block % domain % dminfo

      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'angleEdge', angleEdge)

      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'dzu', dzu)
      call iamas_pool_get_array(mesh, 'fzm', fzm)
      call iamas_pool_get_array(mesh, 'fzp', fzp)
      call iamas_pool_get_array(mesh, 'zz', zz)

      call iamas_pool_get_array(diag, 'exner_base', pb)
      call iamas_pool_get_array(diag, 'rho_base', rb)
      call iamas_pool_get_array(diag, 'theta_base', tb)
      call iamas_pool_get_array(diag, 'rtheta_base', rtb)
      call iamas_pool_get_array(diag, 'exner', p)
      call iamas_pool_get_array(diag, 'pressure_base', ppb)
      call iamas_pool_get_array(diag, 'pressure_p', pp)
      call iamas_pool_get_array(diag, 'pressure', pressure)
      call iamas_pool_get_array(diag, 'relhum', relhum)
      call iamas_pool_get_array(diag, 'spechum', spechum)
      call iamas_pool_get_array(diag, 'ru', ru)
      call iamas_pool_get_array(diag, 'rw', rw)

      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(diag, 'rho_p', rr)
      call iamas_pool_get_array(state, 'theta_m', t)
      call iamas_pool_get_array(diag, 'rtheta_p', rt)
      call iamas_pool_get_array(lbc_state, 'lbc_scalars', scalars)
      call iamas_pool_get_array(lbc_state, 'lbc_u', u)
      call iamas_pool_get_array(lbc_state, 'lbc_w', w)
      call iamas_pool_get_array(lbc_state, 'lbc_theta', theta)
      call iamas_pool_get_array(lbc_state, 'lbc_rho', rho)

      call iamas_pool_get_array(mesh, 'latCell', latCell)
      call iamas_pool_get_array(mesh, 'lonCell', lonCell)
      call iamas_pool_get_array(mesh, 'latEdge', latEdge)
      call iamas_pool_get_array(mesh, 'lonEdge', lonEdge)

      call iamas_pool_get_array(fg, 'u', u_fg)
      call iamas_pool_get_array(fg, 'v', v_fg)
      call iamas_pool_get_array(fg, 'z', z_fg)
      call iamas_pool_get_array(fg, 't', t_fg)
      call iamas_pool_get_array(fg, 'rh', rh_fg)
      call iamas_pool_get_array(fg, 'sh', sh_fg)
      call iamas_pool_get_array(fg, 'p', p_fg)

      call iamas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
      nVertLevelsP1 = nVertLevels + 1

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)
      call iamas_pool_get_dimension(lbc_state, 'index_lbc_tke_sgs', index_tke_sgs)

      etavs = (1.0_RKIND - 0.252_RKIND) * pii / 2.0_RKIND
      rcv = rgas / (cp - rgas)
      p0 = 1.0e+05_RKIND

      scalars(:,:,:) = 0.0_RKIND
      if (associated(config_les_scheme) .and. associated(config_les_lbc_tke_mode) .and. &
          associated(config_les_tke_min) .and. associated(index_tke_sgs)) then
         if (trim(config_les_scheme) == 'les_deardorff' .and. index_tke_sgs > 0) then
            select case (trim(config_les_lbc_tke_mode))
            case ('diagnose', 'minimum')
               scalars(index_tke_sgs,:,:) = config_les_tke_min
            case ('input')
               ! The generic meteorological intermediate format has no SGS-TKE field.
               ! Keep the safe lower bound unless a future input reader fills this field.
               scalars(index_tke_sgs,:,:) = config_les_tke_min
            case default
               call iamas_log_write('Invalid config_les_lbc_tke_mode: '// &
                                   trim(config_les_lbc_tke_mode), messageType=iAMAS_LOG_CRIT)
            end select
         end if
      end if

      !
      ! Check that we have what looks like a valid zgrid field. If the max value for zgrid is zero,
      ! the input file likely does not contain vertical grid information.
      !
      max_zgrid_local = maxval(zgrid(:,1:nCellsSolve))
      call iamas_dmpar_max_real(dminfo, max_zgrid_local, max_zgrid_global)
      if (max_zgrid_global == 0.0_RKIND) then
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('The maximum value of the zgrid field is 0. Please ensure that the ''input'' stream ', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('contains valid vertical grid information.', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_CRIT)
      end if


      !
      ! Horizontally interpolate meteorological data
      !
      allocate(vert_level(config_nfglevels))
      vert_level(:) = -1.0
      allocate(found_u(config_nfglevels), found_v(config_nfglevels), found_t(config_nfglevels), &
               found_moisture(config_nfglevels), found_z(config_nfglevels))
      found_u(:) = .false.
      found_v(:) = .false.
      found_t(:) = .false.
      found_moisture(:) = .false.
      found_z(:) = .false.
      found_soilhgt = .false.

      ! TODO: We should check that timestamp is actually of length >= 13
      call read_met_init(trim(config_met_prefix), .false., timestamp(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error opening initial meteorological data file '//trim(config_met_prefix)//':'//timestamp(1:13), &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_CRIT)
      end if

      call iamas_hash_init(level_hash)
      too_many_fg_levs = .false.

      call read_next_met_field(field, istatus)

      do while (istatus == 0)

         interp_list(1) = FOUR_POINT
         interp_list(2) = SEARCH
         interp_list(3) = 0

         msgval = -1.e30

         if (trim(field % field) == 'UU' .or. &
             trim(field % field) == 'VV' .or. &
             trim(field % field) == 'TT' .or. &
             trim(field % field) == 'RH' .or. &
             trim(field % field) == 'SPECHUMD' .or. &
             trim(field % field) == 'GHT' .or. &
             trim(field % field) == 'SOILHGT' .or. &
             trim(field % field) == 'PRES' .or. &
             trim(field % field) == 'PRESSURE') then

            if (trim(field % field) /= 'SOILHGT') then

               ! Since the hash table can only store integers, transfer the bit pattern from
               ! the real-valued xlvl into an integer; that the result is not an integer version
               ! of the level is not important, since we only want to test uniqueness of levels
               level_value = transfer(field % xlvl, level_value)
               if (.not. iamas_hash_search(level_hash, level_value)) then
                  call iamas_hash_insert(level_hash, level_value)
                  if (iamas_hash_size(level_hash) > config_nfglevels) then
                     too_many_fg_levs = .true.
                  end if
               end if

               !
               ! In case we have more than config_nfglevels levels, just keep cycling through
               ! the remaining fields in the intermediate file for the purpose of counting how
               ! many unique levels are found using the code above
               !
               if (too_many_fg_levs) then
                  call read_next_met_field(field, istatus)
                  cycle
               end if

               do k=1,config_nfglevels
                  if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
               end do
               if (vert_level(k) == -1.0) vert_level(k) = field % xlvl
            else
               k = 1
            end if

            if (trim(field % field) == 'UU') then
               found_u(k) = .true.
            else if (trim(field % field) == 'VV') then
               found_v(k) = .true.
            else if (trim(field % field) == 'TT') then
               found_t(k) = .true.
            else if (trim(field % field) == 'RH' .or. trim(field % field) == 'SPECHUMD') then
               found_moisture(k) = .true.
            else if (trim(field % field) == 'GHT') then
               found_z(k) = .true.
            else if (trim(field % field) == 'SOILHGT') then
               found_soilhgt = .true.
            end if

            !
            ! Set up projection
            !
            call map_init(proj)

            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat,RKIND), &
                            loninc = real(field % deltalon,RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            end if


            !
            ! Horizontally interpolate the field at level k
            !
            if (trim(field % field) == 'UU') then
               call iamas_log_write('Interpolating U at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nEdges
               latPoints => latEdge
               lonPoints => lonEdge
               call iamas_pool_get_array(fg, 'u', destField2d)
               ndims = 2
            else if (trim(field % field) == 'VV') then
               call iamas_log_write('Interpolating V at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nEdges
               latPoints => latEdge
               lonPoints => lonEdge
               call iamas_pool_get_array(fg, 'v', destField2d)
               ndims = 2
            else if (trim(field % field) == 'TT') then
               call iamas_log_write('Interpolating TT at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 't', destField2d)
               ndims = 2
            else if (trim(field % field) == 'RH') then
               call iamas_log_write('Interpolating RH at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'rh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'SPECHUMD') then
               call iamas_log_write('Interpolating SPECHUMD at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'sh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'GHT') then
               call iamas_log_write('Interpolating GHT at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'z', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRES') then
               call iamas_log_write('Interpolating PRES at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRESSURE') then
               call iamas_log_write('Interpolating PRESSURE at $i $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'SOILHGT') then
               call iamas_log_write('Interpolating SOILHGT')
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'soilz', destField1d)
               ndims = 1
            end if

            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               lat = latPoints(i)*DEG_PER_RAD
               lon = lonPoints(i)*DEG_PER_RAD
               call latlon_to_ij(proj, lat, lon, x, y)
               if (x < 0.5) then
                  lon = lon + 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               else if (x >= real(field%nx)+0.5) then
                  lon = lon - 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               end if
               if (y < 0.5) then
                  y = 1.0
               else if (y >= real(field%ny)+0.5) then
                  y = real(field%ny)
               end if
               if (ndims == 1) then
                  destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, 1, 1, msgval, interp_list, 1)
               else if (ndims == 2) then
                  destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, 1, 1, msgval, interp_list, 1)
               end if
            end do

            deallocate(rslab)

         end if

         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()
      level_value = iamas_hash_size(level_hash)
      call iamas_hash_destroy(level_hash)

      missing_required_fields = 0
      do k=1,config_nfglevels
         if (vert_level(k) == -1.0 .or. vert_level(k) == 200100.0) cycle
         if (.not. found_u(k) .or. .not. found_v(k) .or. .not. found_t(k) .or. &
             .not. found_moisture(k) .or. .not. found_z(k)) then
            missing_required_fields = missing_required_fields + 1
            call iamas_log_write('Incomplete required field set at meteorological level $r.', &
                                realArgs=(/vert_level(k)/), messageType=iAMAS_LOG_ERR)
         end if
      end do
      if (config_blend_bdy_terrain .and. .not. found_soilhgt) then
         missing_required_fields = missing_required_fields + 1
         call iamas_log_write('Missing required SOILHGT field in meteorological input.', messageType=iAMAS_LOG_ERR)
      end if
      if (missing_required_fields > 0) then
         call iamas_log_write('Refusing to interpolate incomplete meteorological input.', messageType=iAMAS_LOG_CRIT)
      end if

      if (too_many_fg_levs) then
         write(errstring,'(a,i4)') '       Please increase config_nfglevels to at least ', level_value
         call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error: The meteorological data file has more than config_nfglevels.', messageType=iAMAS_LOG_ERR)
         call iamas_log_write(trim(errstring),                                                       messageType=iAMAS_LOG_ERR)
         call iamas_log_write('       in the namelist and re-run.',                                  messageType=iAMAS_LOG_ERR)
         call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_CRIT)
      end if


      !
      ! Check how many distinct levels we actually found in the meteorological data
      !
      do k=1,config_nfglevels
         if (vert_level(k) == -1.0) exit
      end do
      nfglevels_actual = k-1
      call iamas_log_write('*************************************************')
      call iamas_log_write('Found $i levels in the first-guess data', intArgs=(/nfglevels_actual/))
      call iamas_log_write('*************************************************')


      !
      ! For isobaric data, fill in the 3-d pressure field; otherwise, ensure
      ! that the surface pressure and height fields are filled in
      !
      if (minval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0 .and. &
          maxval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0) then
         call iamas_log_write('Setting pressure field for isobaric data')
         do k=1,config_nfglevels
            if (vert_level(k) /= 200100.0) then
               p_fg(k,:) = vert_level(k)
            end if
         end do
      else
         call iamas_pool_get_array(fg, 'z', z_fg)
         call iamas_pool_get_array(fg, 'soilz', soilz)
         call iamas_log_write('Assuming model-level input data')
         do k=1,config_nfglevels
            if (vert_level(k) == 200100.0) then
               z_fg(k,:) = soilz(:)
            end if
         end do
      end if


      !
      ! Compute normal wind component and store in fg % u
      !
      do iEdge=1,nEdges
         do k=1,nfglevels_actual
            u_fg(k,iEdge) = cos(angleEdge(iEdge)) * u_fg(k,iEdge) &
                          + sin(angleEdge(iEdge)) * v_fg(k,iEdge)
         end do
      end do

      !
      ! Vertically interpolate meteorological data
      !
      allocate(sorted_arr(2,nfglevels_actual))

      do iCell=1,nCells

         ! T
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = t_fg(k,iCell)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            t(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                         sorted_arr(:,1:nfglevels_actual-1), order=1, &
                                         extrap=extrap_airtemp, ierr=istatus)
            if (istatus /= 0) then
               write(errstring,'(a,i4,a,i10)') 'Error in interpolation of t(k,iCell) for k=', k, ', iCell=', iCell
               call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_ERR)
               call iamas_log_write(trim(errstring),                                                     messageType=iAMAS_LOG_ERR)
               call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_CRIT)
            end if
         end do


         ! RH
         sorted_arr(:,:) = -999.0
         relhum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = rh_fg(k,iCell)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            relhum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                       sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
         end do


         ! SPECHUM: if first-guess values are negative, set those values to zero before
         ! vertical interpolation.
         sorted_arr(:,:) = -999.0
         spechum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = max(0._RKIND,sh_fg(k,iCell))
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            spechum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                        sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
            if (target_z < z_fg(1,iCell) .and. k < nVertLevels) spechum(k,iCell) = spechum(k+1,iCell)
         end do


         ! PRESSURE
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) then
               sorted_arr(1,k) = 99999.0
               sfc_k = k
               p_fg(k,iCell) = 1.0   ! Any value that has valid log is fine...
            end if
            sorted_arr(2,k) = log(p_fg(k,iCell))
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            pressure(k,iCell) = exp(vertical_interp(target_z, nfglevels_actual-1, &
                                    sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1))
         end do

      end do


      do iEdge=1,nEdges

         ! U
         sorted_arr(:,:) = -999.0
         do k=1,nfglevels_actual
            sorted_arr(1,k) = 0.5 * (z_fg(k,cellsOnEdge(1,iEdge)) + z_fg(k,cellsOnEdge(2,iEdge)))
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = u_fg(k,iEdge)
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         do k=1,nVertLevels
            target_z = 0.25 * (zgrid(k,cellsOnEdge(1,iEdge)) + zgrid(k+1,cellsOnEdge(1,iEdge)) &
                             + zgrid(k,cellsOnEdge(2,iEdge)) + zgrid(k+1,cellsOnEdge(2,iEdge)))
            u(k,iEdge) = vertical_interp(target_z, nfglevels_actual-1, sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
         end do

      end do

      deallocate(sorted_arr)


      ! Diagnose the water vapor mixing ratios:
      global_sh_min = 0._RKIND
      global_sh_max = 0._RKIND
      if(config_use_spechumd) then
         sh_min = minval(spechum(:,1:nCellsSolve))
         sh_max = maxval(spechum(:,1:nCellsSolve))
         call iamas_dmpar_min_real(dminfo,sh_min,global_sh_min)
         call iamas_dmpar_max_real(dminfo,sh_max,global_sh_max)
      endif
      call iamas_log_write('')
      call iamas_log_write('--- global_sh_min = $r', realArgs=(/global_sh_min/))
      call iamas_log_write('--- global_sh_max = $r', realArgs=(/global_sh_max/))
      call iamas_log_write('')

      call iamas_log_write('--- config_use_spechumd = $l', logicArgs=(/config_use_spechumd/))
      if(.not. config_use_spechumd .or. (global_sh_min==0._RKIND .and. global_sh_max==0._RKIND)) then
         !--- calculate the saturation mixing ratio and interpolated first-guess relative humidity:
         if (config_use_spechumd) then
            call iamas_log_write('config_use_spechumd=T, but specific humidity was not found in ' &
                                //trim(config_met_prefix)//':'//timestamp(1:13), messageType=iAMAS_LOG_WARN)
         end if
         call iamas_log_write(' *** initializing water vapor mixing ratio using first-guess relative humidity')
         call iamas_log_write('')

         do k = 1, nVertLevels
            do iCell = 1, nCells
               !
               ! Note: the RH field provided by ungrib should always be with respect to liquid water,
               !       hence, we can always call rslf; see the routine fix_gfs_rh in WPS/ungrib/src/rrpr.F .
               !
               rs = rslf(pressure(k,iCell),t(k,iCell))
               scalars(index_qv,k,iCell) = 0.01_RKIND*rs*relhum(k,iCell)
             enddo
         enddo
      else
         !--- use the interpolated first-guess specific humidity:
         call iamas_log_write(' *** initializing water vapor mixing ratio using first-guess specific humidity')
         call iamas_log_write('')
         do k = 1, nVertLevels
            do iCell = 1, nCells
               scalars(index_qv,k,iCell) = spechum(k,iCell)/(1._RKIND-spechum(k,iCell))
            enddo
         enddo
      endif

      !
      ! Diagnose fields needed in initial conditions file (u, w, rho, theta)
      ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
      !
      do iCell=1,nCells

         do k=1,nVertLevels
            ! PI
            p(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)

            ! THETA - can compute this using PI instead
            t(k,iCell) = t(k,iCell) * (p0 / pressure(k,iCell)) ** (rgas / cp)

            ! RHO_ZZ
            rho_zz(k,iCell) = pressure(k,iCell) / rgas / (p(k,iCell) * t(k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / (1.0 + scalars(index_qv,k,iCell))
         end do
      end do


      !
      ! Reference state based on a dry isothermal atmosphere
      !
      do iCell=1,nCells
         do k=1,nVertLevels
            ztemp    = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
            ppb(k,iCell) = p0*exp(-gravity*ztemp/(rgas*t0b))      ! pressure_base
            pb (k,iCell) = (ppb(k,iCell)/p0)**(rgas/cp)           ! exner_base
            rb (k,iCell) = ppb(k,iCell)/(rgas*t0b)                ! rho_base
            tb (k,iCell) = t0b/pb(k,iCell)                        ! theta_base
            rtb(k,iCell) = rb(k,iCell)*tb(k,iCell)                ! rtheta_base
            p  (k,iCell) = pb(k,iCell)                            ! exner
            pp (k,iCell) = 0.                                     ! pressure_p
            rr (k,iCell) = 0.                                     ! rho_p
         end do
      end do

      do iCell=1,nCells
         do k=1,nVertLevels

            ! couple with vertical metric
            rb(k,iCell) = rb(k,iCell) / zz(k,iCell)
            rho_zz(k,iCell) = rho_zz(k,iCell) / zz(k,iCell)

            pp(k,iCell) = pressure(k,iCell) - ppb(k,iCell)
            rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         end do
      end do

      do iCell=1,nCells
         k = 1

         ! couple with vertical metric, note: rr is coupled here
         rho_zz(k,iCell) = ((pressure(k,iCell) / p0)**(cv / cp)) * (p0 / rgas) &
                            / (t(k,iCell)*(1.0 + 1.61*scalars(index_qv,k,iCell))) / zz(k,iCell)
         rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         do k=2,nVertLevels
            it = 0
            p_check = 2.0 * 0.0001
            do while ( (it < 30) .and. (p_check > 0.0001) )

               p_check = pp(k,iCell)

               ! iAMAS hydrostatic relation
               pp(k,iCell) = pp(k-1,iCell) - (fzm(k)*rr(k,iCell) + fzp(k)*rr(k-1,iCell))*gravity*dzu(k) &
                                           - (fzm(k)*rho_zz(k,iCell)*scalars(index_qv,k,iCell) &
                                                  + fzp(k)*rho_zz(k-1,iCell)*scalars(index_qv,k-1,iCell))*gravity*dzu(k)
               pressure(k,iCell) = pp(k,iCell) + ppb(k,iCell)
               p(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)

               ! couple with vertical metric
               rho_zz(k,iCell) = pressure(k,iCell) / rgas &
                     / (p(k,iCell)*t(k,iCell)*(1.0 + 1.61*scalars(index_qv,k,iCell)))/zz(k,iCell)
               rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

               p_check = abs(p_check - pp(k,iCell))

               it = it + 1
            end do
         end do
      end do

      ! Compute theta_m and rho-tilde
      do iCell=1,nCells
         do k=1,nVertLevels
            t(k,iCell) = t(k,iCell) * (1.0 + 1.61*scalars(index_qv,k,iCell))
            rr(k,iCell) = rr(k,iCell)*zz(k,iCell)
         end do
      end do

      do iEdge=1,nEdges
         do k=1,nVertLevels
            ru(k,iEdge) = u(k,iEdge) * 0.5*(rho_zz(k,cellsOnEdge(1,iEdge)) + rho_zz(k,cellsOnEdge(2,iEdge)))
         end do
      end do


      rw(:,:) = 0.0

      do iCell=1,nCellsSolve

         do i=1,nEdgesOnCell(iCell)
            iEdge=edgesOnCell(i,iCell)

            do k = 2, nVertLevels
               flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
               if (iCell == cellsOnEdge(1,iEdge)) then
                  rw(k,iCell) = rw(k,iCell) - (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb(k,1,iEdge)*flux
               else
                  rw(k,iCell) = rw(k,iCell) + (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb(k,2,iEdge)*flux
               end if

               if (config_theta_adv_order ==3) then
                  if (iCell == cellsOnEdge(1,iEdge)) then
                     rw(k,iCell) = rw(k,iCell)    &
                                  + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                    (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb3(k,1,iEdge)*flux
                  else
                     rw(k,iCell) = rw(k,iCell)    &
                                  - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                    (fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))*zb3(k,2,iEdge)*flux
                  end if
               end if

            end do

         end do

      end do


      ! Compute w from rho_zz and rw
      do iCell=1,nCellsSolve
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do

      deallocate(vert_level)
      deallocate(found_u, found_v, found_t, found_moisture, found_z)


      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_lbc


   integer function nearest_edge(target_lat, target_lon, &
                                 start_edge, &
                                 nCells, nEdges, maxEdges, nEdgesOnCell, edgesOnCell, cellsOnEdge, latCell, lonCell, latEdge, lonEdge)

      implicit none

      real (kind=RKIND), intent(in) :: target_lat, target_lon
      integer, intent(in) :: start_edge
      integer, intent(in) :: nCells, nEdges, maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
      real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
      real (kind=RKIND), dimension(nEdges), intent(in) :: latEdge, lonEdge

      integer :: i, cell1, cell2, iCell
      integer :: iEdge
      integer :: current_edge
      real (kind=RKIND) :: cell1_dist, cell2_dist
      real (kind=RKIND) :: current_distance, d
      real (kind=RKIND) :: nearest_distance

      nearest_edge = start_edge
      current_edge = -1

      do while (nearest_edge /= current_edge)
         current_edge = nearest_edge
         current_distance = sphere_distance(latEdge(current_edge), lonEdge(current_edge), target_lat, target_lon, 1.0_RKIND)
         nearest_edge = current_edge
         nearest_distance = current_distance
         cell1 = cellsOnEdge(1,current_edge)
         cell2 = cellsOnEdge(2,current_edge)
         cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), target_lat, target_lon, 1.0_RKIND)
         cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), target_lat, target_lon, 1.0_RKIND)
         if (cell1_dist < cell2_dist) then
            iCell = cell1
         else
            iCell = cell2
         end if
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            if (iEdge <= nEdges) then
               d = sphere_distance(latEdge(iEdge), lonEdge(iEdge), target_lat, target_lon, 1.0_RKIND)
               if (d < nearest_distance) then
                  nearest_edge = iEdge
                  nearest_distance = d
               end if
            end if
         end do
      end do

   end function nearest_edge


   real (kind=RKIND) function vertical_interp(target_z, nz, zf, order, extrap, surface_val, sealev_val, ierr)

      implicit none

      real (kind=RKIND), intent(in) :: target_z
      integer, intent(in) :: nz
      real (kind=RKIND), dimension(2,nz), intent(in) :: zf      ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
      integer, intent(in), optional :: order
      integer, intent(in), optional :: extrap                   ! can take values 0 = constant, 1 = linear (default), 2 = lapse-rate
      real (kind=RKIND), intent(in), optional :: surface_val
      real (kind=RKIND), intent(in), optional :: sealev_val
      integer, intent(out), optional :: ierr

      integer :: k, lm, lp
      real (kind=RKIND) :: wm, wp
      real (kind=RKIND) :: slope

      integer :: interp_order, extrap_type
      real (kind=RKIND) :: surface, sealevel

      if (present(ierr)) ierr = 0

      if (present(order)) then
         interp_order = order
      else
         interp_order = 2
      end if

      if (present(extrap)) then
         extrap_type = extrap
      else
         extrap_type = 1
      end if

      if (present(surface_val)) then
         surface = surface_val
      else
         surface = 200100.0
      end if

      if (present(sealev_val)) then
         sealevel = sealev_val
      else
         sealevel = 201300.0
      end if

      !
      ! Extrapolation required
      !
      if (target_z < zf(1,1)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,1)
         else if (extrap_type == 1) then
            slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
            vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
         else if (extrap_type == 2) then
            vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
         end if
         return
      end if
      if (target_z >= zf(1,nz)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,nz)
         else if (extrap_type == 1) then
            slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
            vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
         else if (extrap_type == 2) then
             call iamas_log_write('extrap_type == 2 not implemented for target_z >= zf(1,nz)', messageType=iAMAS_LOG_ERR)
             if (present(ierr)) ierr = 1
             return
         end if
         return
      end if


      !
      ! No extrapolation required
      !
      do k=1,nz-1
         if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
            lm = k
            lp = k+1
            wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
            wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
            exit
         end if
      end do

      vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

      return

   end function vertical_interp


!----------------------------------------------------------------------------------------------------------

   real (kind=RKIND) function env_qv( z, temperature, pressure, rh_max )

      implicit none
      real (kind=RKIND) :: z, temperature, pressure, ztr, es, qvs, p0, rh_max

      p0 = 100000.

!      ztr = 5000.
!
!      if(z .gt. ztr) then
!         env_qv = 0.
!      else
!         if(z.lt.2000.) then
!            env_qv = .5
!         else
!            env_qv = .5*(1.-(z-2000.)/(ztr-2000.))
!         end if
!      end if

       if (pressure .lt. 50000. ) then
           env_qv = 0.0
       else
           env_qv = (1.-((p0-pressure)/50000.)**1.25)
       end if

       env_qv = min(rh_max,env_qv)

! env_qv is the relative humidity, turn it into mixing ratio
       if (temperature .gt. 273.15) then
           es  = 1000.*0.6112*exp(17.67*(temperature-273.15)/(temperature-29.65))
       else
           es  = 1000.*0.6112*exp(21.8745584*(temperature-273.16)/(temperature-7.66))
       end if
       qvs = (287.04/461.6)*es/(pressure-es)

       ! qvs =  380.*exp(17.27*(temperature-273.)/(temperature-36.))/pressure

        env_qv = env_qv*qvs

   end function env_qv


   subroutine physics_idealized_init(mesh, fg)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      type (iamas_pool_type), intent(inout) :: fg

      !local variables:
      integer :: iCell, iMonth, iSoil,iLanduse,iSoilType
      integer, pointer :: nCells, nSoilLevels, nMonths
      integer, dimension(:), pointer :: landmask, lu_index, soilcat_top
      real (kind=RKIND), dimension(:), pointer :: ter, xice, shdmin, shdmax, vegfra, sfc_albbck, xland, seaice
      real (kind=RKIND), dimension(:), pointer :: snow, snowc, snoalb, snowh, skintemp, sst, tmn
      real (kind=RKIND), dimension(:,:), pointer :: tslb, smcrel, sh2o, smois, dzs, albedo12m, greenfrac,soilctop,landusef

      !---------------------------------------------------------------------------------------------

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nSoilLevels', nSoilLevels)
      call iamas_pool_get_dimension(mesh, 'nMonths', nMonths)

      call iamas_pool_get_array(mesh, 'ter', ter)
      call iamas_pool_get_array(mesh, 'landmask', landmask)
      call iamas_pool_get_array(mesh, 'lu_index', lu_index)
      call iamas_pool_get_array(mesh, 'soilcat_top', soilcat_top)
      call iamas_pool_get_array(mesh, 'soilctop', soilctop)
      call iamas_pool_get_array(mesh, 'landusef', landusef)
      call iamas_pool_get_array(mesh, 'shdmin', shdmin)
      call iamas_pool_get_array(mesh, 'shdmax', shdmax)
      call iamas_pool_get_array(mesh, 'snoalb', snoalb)
      call iamas_pool_get_array(mesh, 'albedo12m', albedo12m)
      call iamas_pool_get_array(mesh, 'greenfrac', greenfrac)

      call iamas_pool_get_array(fg, 'xice', xice)
      call iamas_pool_get_array(fg, 'vegfra', vegfra)
      call iamas_pool_get_array(fg, 'sfc_albbck', sfc_albbck)
      call iamas_pool_get_array(fg, 'xland', xland)
      call iamas_pool_get_array(fg, 'seaice', seaice)
      call iamas_pool_get_array(fg, 'snow', snow)
      call iamas_pool_get_array(fg, 'snowc', snowc)
      call iamas_pool_get_array(fg, 'snowh', snowh)
      call iamas_pool_get_array(fg, 'skintemp', skintemp)
      call iamas_pool_get_array(fg, 'sst', sst)
      call iamas_pool_get_array(fg, 'tmn', tmn)
      call iamas_pool_get_array(fg, 'tslb', tslb)
      call iamas_pool_get_array(fg, 'smcrel', smcrel)
      call iamas_pool_get_array(fg, 'sh2o', sh2o)
      call iamas_pool_get_array(fg, 'smois', smois)
      call iamas_pool_get_array(fg, 'dzs', dzs)

      !initialization of surface input variables that are not needed if we run the current set of
      !idealized test cases:


      do iCell = 1, nCells

         !terrain,soil type, and vegetation:
         ter(iCell) = 0.0
         xice(iCell) = 0.0
         landmask(iCell) = 0
         lu_index(iCell) = 0
         soilcat_top(iCell) = 0
         shdmin(iCell) = 0.0
         shdmax(iCell) = 0.0
         vegfra(iCell) = 0.0
         sfc_albbck(iCell) = 0.0
         xland(iCell) = 0.0
         seaice(iCell) = 0.0

         do iLanduse = 1,24
            landusef(iLanduse,iCell) = 0.0
         end do

         do iSoilType = 1,16
            soilctop(iSoilType,iCell) = 0.0
         end do
         !snow coverage:
         snow(iCell) = 0.0
         snowc(iCell) = 0.0
         snoalb(iCell) = 0.08
         snowh(iCell) = 0.0

         !surface and sea-surface temperatures:
         skintemp(iCell) = 288.0
         sst(iCell) = 288.0

         !soil layers:
         tmn(iCell) = 288.0
         do iSoil = 1, nSoilLevels
            tslb(iSoil,iCell)   = 288.0
            smcrel(iSoil,iCell) =   0.0
            sh2o(iSoil,iCell) =   0.0
            smois(iSoil,iCell) =   0.0
            dzs(iSoil,iCell) =   0.0
         end do

         !monthly climatological surface albedo and greeness fraction:
         do iMonth = 1, nMonths
            albedo12m(iMonth,iCell) = 0.08
            greenfrac(iMonth,iCell) = 0.0
         end do

      end do

   end subroutine physics_idealized_init


   subroutine decouple_variables(mesh, nCells, nVertLevels, state, diag)

      implicit none

      type (iamas_pool_type), intent(in) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (iamas_pool_type), intent(inout) :: state
      type (iamas_pool_type), intent(inout) :: diag

      integer :: iCell, iEdge, k

      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: rdzw
      real (kind=RKIND), dimension(:,:), pointer :: zz, pp, ppb, rho, rho_zz, theta, theta_m
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer, pointer :: index_qv

      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(diag, 'pressure_p', pp)
      call iamas_pool_get_array(diag, 'pressure_base', ppb)
      call iamas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag, 'theta', theta)
      call iamas_pool_get_array(state, 'rho_zz', rho_zz)
      call iamas_pool_get_array(state, 'theta_m', theta_m)
      call iamas_pool_get_array(state, 'scalars', scalars)

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      ! Compute surface pressure
      do iCell=1,nCells
         surface_pressure(iCell) = 0.5*gravity/rdzw(1)                                        &
                                   * (1.25* rho_zz(1,iCell) * (1. + scalars(index_qv, 1, iCell))  &
                                      -  0.25* rho_zz(2,iCell) * (1. + scalars(index_qv, 2, iCell)))
         surface_pressure(iCell) = surface_pressure(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do


      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = theta_m(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine decouple_variables


   !-----------------------------------------------------------------------
   !  routine blend_bdy_terrain
   !
   !> \brief Combines first-guess terrain with static terrain along regional domain boundaries
   !> \author Michael Duda
   !> \date   25 April 2019
   !> \details
   !>  This routine combines terrain from the first-guess dataset provided in an intermediate
   !>  file with the terrain field produced by the init_atmosphere core's "static interpolation"
   !>  stage in the boundary cells of a regional mesh. Specifically, where the value of
   !>  the bdyMaskCell field is nBdyLayers or nBdyLayers-1, the terrain field, ter, is interpolated
   !>  directly from the first-guess terrain; where the value of bdy MaskCell is between nBdyLayers-2
   !>  and 1, the terrain field is a combination of the first-guess terrain and the high-resolution
   !>  "static" terrain field.
   !>
   !>  When dryrun=true, the nCells, latCell, lonCell, bdyMaskCell, and ter arguments are not used
   !>   -- only the config_met_prefix and config_start_time arguments are used. The dryrun argument
   !>  allows calling code to determine if blend_bdy_terrain will succeed without actually blending
   !>  the boundary terrain.
   !>
   !>  For global meshes, where bdyMaskCell == 0 everywhere, this routine will have no impact
   !>  on the model terrain field.
   !
   !-----------------------------------------------------------------------
   subroutine blend_bdy_terrain(config_met_prefix, config_start_time, nCells, latCell, lonCell, bdyMaskCell, ter, dryrun, ierr)

      use init_atm_read_met, only : read_met_init, read_met_close, read_next_met_field, met_data
      use init_atm_llxy, only : map_init, map_set, proj_info, latlon_to_ij, PROJ_LATLON, PROJ_GAUSS, DEG_PER_RAD
      use init_atm_hinterp, only : interp_sequence, FOUR_POINT

      implicit none

      character(len=*), intent(in) :: config_met_prefix
      character(len=*), intent(in) :: config_start_time
      integer, intent(in) :: nCells
      real (kind=RKIND), dimension(:), intent(in) :: latCell  ! These four variables (latCell, lonCell, bdyMaskCell, and ter)
      real (kind=RKIND), dimension(:), intent(in) :: lonCell  ! may actually have more than nCells elements, for example,
      integer, dimension(:), intent(in) :: bdyMaskCell        ! if the arrays include a "garbage cell".
      real (kind=RKIND), dimension(:), intent(inout) :: ter   ! 
      logical, intent(in) :: dryrun
      integer, intent(out) :: ierr

      integer, parameter :: nBdyLayers = 7   ! The number of relaxation layers plus the number of specified layers
      integer, parameter :: nSpecLayers = 2  ! The number of specified layers

      integer :: i
      integer :: istatus
      integer, dimension(2) :: interp_list
      real (kind=RKIND) :: weight
      real (kind=RKIND) :: lat, lon
      real (kind=RKIND) :: x, y
      real (kind=RKIND), allocatable, dimension(:,:) :: rslab
      type (met_data) :: field
      type (proj_info) :: proj


      ierr = 0

      if (.not. dryrun) then
         call iamas_log_write('Blending first-guess terrain field along domain boundary')
      end if

      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error opening file with terrain field, '//trim(config_met_prefix)//':'//config_start_time(1:13), &
                             messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', &
                             messageType=iAMAS_LOG_ERR)
         ierr = 1
         return
      end if

      !
      ! Loop over fields in the intermediate file looking for the SOILHGT field
      !
      call read_next_met_field(field, istatus)
      do while (istatus == 0)
         if (trim(field % field) == 'SOILHGT') then

            if (.not. dryrun) then
               interp_list(1) = FOUR_POINT
               interp_list(2) = 0

               !
               ! Set up map projection - currently, only the regular lat-lon projection is handled
               !
               call map_init(proj)
          
               if (field % iproj == PROJ_LATLON) then
                  call map_set(PROJ_LATLON, proj, &
                               latinc = real(field % deltalat,RKIND), &
                               loninc = real(field % deltalon,RKIND), &
                               knowni = 1.0_RKIND, &
                               knownj = 1.0_RKIND, &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               else if (field % iproj == PROJ_GAUSS) then
                  call map_set(PROJ_GAUSS, proj, &
                               nlat = nint(field % deltalat), &
                               loninc = 360.0_RKIND / real(field % nx,RKIND), &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               end if

               !
               ! Copy the first-guess terrain field into an array that includes some periodic points
               !
               allocate(rslab(-2:field % nx+3, field % ny))
               rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
               rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
               rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
               rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
               rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
               rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
               rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

               !
               ! For each cell in the iAMAS mesh, perform terrain blending if the cell is a boundary cell
               !
               do i=1,nCells
                  if (bdyMaskCell(i) > 0) then
                     lat = latCell(i)*DEG_PER_RAD
                     lon = lonCell(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     else if (x >= real(field%nx)+0.5) then
                        lon = lon - 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (y < 0.5) then
                        y = 1.0
                     else if (y >= real(field%ny)+0.5) then
                        y = real(field%ny)
                     end if

                     !
                     ! Is this a specified cell?
                     !
                     if (bdyMaskCell(i) > (nBdyLayers - nSpecLayers)) then
                        ter(i) = interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, &
                                                 1, 1, -1.0E30_RKIND, interp_list, 1)

                     !
                     ! Or a relaxation cell?
                     !
                     else
                        weight = real(bdyMaskCell(i),kind=RKIND) / real(nBdyLayers - nSpecLayers,kind=RKIND)
                        ter(i) = weight * interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, &
                                                          1, 1, -1.0E30_RKIND, interp_list, 1) &
                               + (1.0_RKIND - weight) * ter(i)

                     end if
                  end if
               end do

               deallocate(rslab)
            end if

            !
            ! At this point, we have found and processed the first-guess terrain field, so we can return
            !
            deallocate(field % slab)
            call read_met_close()
            return

         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()

      !
      ! If we have reached this point, no first-guess terrain field was found...
      !
      ierr = 1
      call iamas_log_write('********************************************************************************', &
                          messageType=iAMAS_LOG_ERR)
      call iamas_log_write('SOILHGT field not found in intermediate file ' &
                                     //trim(config_met_prefix)//':'//config_start_time(1:13)                , &
                          messageType=iAMAS_LOG_ERR)
      call iamas_log_write('********************************************************************************', &
                          messageType=iAMAS_LOG_ERR)

   end subroutine blend_bdy_terrain


   !-----------------------------------------------------------------------
   !  routine convert_relhum_wrt_ice
   !
   !> \brief Converts an RH field given w.r.t. liquid water to RH w.r.t. ice below freezing
   !> \author Wei Wang
   !> \date   11 May 2019
   !> \details
   !>  This routine takes as input a temperature field (in K) and an RH field (in percent),
   !>  which is assumed to be with respect to liquid water everywhere. Upon return, the
   !>  relative humidity for temperatures below 253.15 K has been modified to be with
   !>  respect to ice. For temperatures in the range (253.15, 273.15] the relative humidity
   !>  uses a blend of the saturation mixing ratios for liquid water and ice.
   !>
   !>  This routine uses the formula from the WPS ungrib program to re-compute RH with
   !>  respect to ice in an attempt to re-construct RH that is more like the ungrib input.
   !
   !-----------------------------------------------------------------------
   subroutine convert_relhum_wrt_ice(t, relhum)

      implicit none

      real (kind=RKIND), dimension(:,:), intent(in) :: t
      real (kind=RKIND), dimension(:,:), intent(inout) :: relhum

      integer :: iCell, k
      integer :: nVertLevels
      integer :: nCells
      real (kind=RKIND) :: eis, ews, r1


      nVertLevels = size(t, dim=1)
      nCells = size(t, dim=2)

      call iamas_log_write('')
      call iamas_log_write('Recomputing RH w.r.t. ice below freezing')
      call iamas_log_write('')

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if ( t(k,iCell) <= 273.15_RKIND ) then

               ! use formula in ungrib to reconstruct RH 
               eis = 0.01_RKIND * exp (9.550426_RKIND - (5723.265_RKIND / t(k,iCell)) + (3.53068_RKIND * log(t(k,iCell))) &
                   - (0.00728332_RKIND * t(k,iCell)))
               ews = 6.112_RKIND * exp(17.67_RKIND * (t(k,iCell)-273.15_RKIND) / ((t(k,iCell)-273.15_RKIND)+243.5_RKIND))

               ! A linear approximation to the GFS blending region ( -20 C > T < 0 C )
               if ( t(k,iCell) > 253.15_RKIND ) then
                  r1 = ((273.15_RKIND - t(k,iCell)) / 20.0_RKIND)
                  r1 = (r1 * eis) + ((1.0_RKIND-r1)*ews)
               else
                  r1 = eis
               end if
               r1 = max(r1, 1.0e-12_RKIND)
               ews = max(ews, 0.0_RKIND)
               relhum(k,iCell) = ews / r1 * relhum(k,iCell)
               relhum(k,iCell) = min(relhum(k,iCell), 100.0_RKIND)
               relhum(k,iCell) = max(relhum(k,iCell), 0.0_RKIND)
            end if
         end do
      end do

   end subroutine convert_relhum_wrt_ice


   !-----------------------------------------------------------------------
   !  routine read_text_array
   !
   !> \brief Reads a real-valued array from a text file and broadcasts to all tasks
   !> \author Michael Duda
   !> \date   11 May 2019
   !> \details
   !>  This routine reads a list of real values from the specified text file on
   !>  the IO_NODE task and broadcasts the values to all other MPI tasks. Upon
   !>  successful return, the array pointer xarray will have been allocated with
   !>  a size equal to the number of lines in the text file, the array will be filled
   !>  with the values from the file, and a value of 0 will be returned.
   !>
   !>  This routine will print an error message and return a non-zero value if
   !>  any of the following conditions occur:
   !>
   !>  1) The text file does not exist
   !>  2) The text file cannot be opened for reading
   !>  3) The text does not contain readable real values
   !>
   !>  If the return value of this function is non-zero, the xarray pointer
   !>  will be unassociated.
   !
   !-----------------------------------------------------------------------
   function read_text_array(dminfo, filename, xarray) result(ierr)

      use iamas_io_units, only : iamas_new_unit, iamas_release_unit
      use iamas_log, only : iamas_log_write
      use iamas_derived_types, only : iAMAS_LOG_ERR

      implicit none

      type (dm_info), intent(in) :: dminfo
      character(len=*), intent(in) :: filename
      real (kind=RKIND), dimension(:), pointer :: xarray

      integer :: ierr

      integer :: i
      integer :: nlines
      integer :: iunit
      integer :: iexists
      logical :: exists
      real (kind=RKIND) :: rtemp


      ierr = 1
      nullify(xarray)

      !
      ! Check whether the file exists
      !
      if (dminfo % my_proc_id == IO_NODE) then
         inquire(file=filename, exist=exists)
         if (exists) then
            iexists = 1
         else
            iexists = 0
         end if
      end if
      call iamas_dmpar_bcast_int(dminfo, iexists)

      if (iexists == 0) then
         call iamas_log_write('Text file '''//filename//''' does not exist.', messageType=iAMAS_LOG_ERR)
         return
      end if

      !
      ! Count the number of lines in the file
      !
      if (dminfo % my_proc_id == IO_NODE) then
         call iamas_new_unit(iunit)
         open(unit=iunit, file=filename, form='formatted', status='old', iostat=ierr)
         if (ierr /= 0) then
            nlines = -1
         else
            nlines = 0
            read(unit=iunit, fmt=*, iostat=ierr) rtemp
            do while (ierr == 0)
               nlines = nlines + 1
               read(unit=iunit, fmt=*, iostat=ierr) rtemp
            end do
         end if
         close(unit=iunit)
         call iamas_release_unit(iunit)
      end if
      call iamas_dmpar_bcast_int(dminfo, nlines)

      if (nlines <= 0) then
         if (nlines < 0) then
            call iamas_log_write('Text file '''//filename//''' could not be opened for reading.', messageType=iAMAS_LOG_ERR)
         else
            call iamas_log_write('Text file '''//filename//''' contains no readable real values.', messageType=iAMAS_LOG_ERR)
         end if
         ierr = 1
         return
      end if

      !
      ! Allocate output array, read, and broadcast
      !
      allocate(xarray(nlines))
      if (dminfo % my_proc_id == IO_NODE) then
         call iamas_new_unit(iunit)
         open(unit=iunit, file=filename, form='formatted', status='old', iostat=ierr)
         do i=1,nlines
            read(unit=iunit, fmt=*, iostat=ierr) xarray(i)
         end do
         close(unit=iunit)
         call iamas_release_unit(iunit)
      end if

      call iamas_dmpar_bcast_reals(dminfo, nlines, xarray)
      ierr = 0

   end function read_text_array


end module init_atm_cases
