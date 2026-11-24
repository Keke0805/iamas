! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module atm_core

   use iamas_derived_types
   use iamas_pool_routines
   use iamas_dmpar
   use iamas_log, only : iamas_log_write, iamas_log_info

   type (iAMAS_Clock_type), pointer :: clock


   contains


   function atm_core_init(domain, startTimeStamp) result(ierr)

      use iamas_timekeeping
      use iamas_kind_types
      use iamas_stream_manager
      use iamas_atm_dimensions, only : iamas_atm_set_dims
      use iamas_atm_diagnostics_manager, only : iamas_atm_diag_setup
      use iamas_atm_threading, only : iamas_atm_threading_init
      use iamas_timer
      use flcl_util_kokkos_mod
      use iamas_les_init, only : iamas_les_initialize

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
      integer :: ierr

      real (kind=RKIND), pointer :: dt
      type (block_type), pointer :: block

      integer :: i
      logical, pointer :: config_do_restart

      type (iamas_pool_type), pointer :: state
      type (iamas_pool_type), pointer :: mesh
      type (iamas_pool_type), pointer :: diag
      type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: initial_time1, initial_time2
      type (iAMAS_Time_Type) :: startTime

      integer, pointer :: nVertLevels, maxEdges, maxEdges2, num_scalars, num_chem

      real (kind=R8KIND) :: timer_start_time, timer_stop_time !jwfeng
      call kokkos_initialize()
      ierr = 0


      !
      ! Setup threading
      !
      call iamas_atm_threading_init(domain % blocklist, ierr)
      if ( ierr /= 0 ) then
         call iamas_log_write('Threading setup failed for core '//trim(domain % core % coreName), messageType=iAMAS_LOG_CRIT)
      end if


      !
      ! Set up inner dimensions used by arrays in optimized dynamics routines
      !
      call iamas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call iamas_pool_get_dimension(state, 'nVertLevels', nVertLevels)
      call iamas_pool_get_dimension(state, 'maxEdges', maxEdges)
      call iamas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
      call iamas_pool_get_dimension(state, 'num_scalars', num_scalars)
      call iamas_pool_get_dimension(state, 'num_chem', num_chem)
      call iamas_atm_set_dims(nVertLevels, maxEdges, maxEdges2, num_scalars, num_chem)

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock
      iamas_log_info => domain % logInfo


      call iamas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)


      !
      ! If this is a restart run, read the restart stream, else read the input
      ! stream.
      ! Regardless of which stream we read for initial conditions, reset the
      ! input alarms for both input and restart before reading any remaining
      ! input streams.
      !
      if (config_do_restart) then
         call iAMAS_stream_mgr_read(domain % streamManager, streamID='restart', ierr=ierr)
      else
         call iamas_dmpar_get_time(timer_start_time) !jwfeng
         call iamas_timer_start('Stream input in init') !jwfeng
         call iAMAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=ierr)
         call iamas_timer_stop('Stream input in init') !jwfeng
         call iamas_dmpar_get_time(timer_stop_time) !jwfeng
         call iamas_log_write('--- Timing for stream input in init: $r s', &
              realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng
      end if
      if (ierr /= iAMAS_STREAM_MGR_NOERR) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error reading initial conditions',                                                 messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', direction=iAMAS_STREAM_INPUT, ierr=ierr)
      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=iAMAS_STREAM_INPUT, ierr=ierr)

      !
      ! Read all other inputs
      ! For now we don't do this here to match results with previous code; to match requires 
      ! that we read in SST and seaice fields after the call to atm_iamas_init_block()
      !
!      call iAMAS_stream_mgr_read(domain % streamManager, ierr=ierr)
!      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, direction=iAMAS_STREAM_INPUT, ierr=ierr) 

      if (.not. config_do_restart) then
         block => domain % blocklist
         do while (associated(block))
            call iamas_pool_get_subpool(block % structs, 'state', state)
            call iamas_pool_initialize_time_levels(state)
            block => block % next
         end do
      end if

      call iamas_les_initialize(domain)

      !
      ! Read a new data stream for Incremental Analysis Update (IAU), if config_IAU_option /= 'off'    : HA (June-15-2016)
      ! FIXME: should I check xtime for the IAU fields? Maybe not.
      ! Note: Because the 'iau' stream has the 'iau' package attached to it, the iAMAS_stream_mgr_read( )
      !       call here will actually try to read the stream only if IAU is being used in the run.
      !
      call iAMAS_stream_mgr_read(domain % streamManager, streamID='iau', whence=iAMAS_STREAM_NEAREST, ierr=ierr)
      if (ierr /= iAMAS_STREAM_MGR_NOERR) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error reading IAU files',                                                          messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='iau', ierr=ierr)

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = iamas_get_clock_time(clock, iAMAS_START_TIME, ierr)
      call iamas_get_time(startTime, dateTimeString=startTimeStamp) 


      call iamas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call iamas_pool_get_field(state, 'u', u_field, 1)

      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_dmpar_exch_halo_field(u_field)
      call iamas_dmpar_get_time(timer_stop_time) !jwfeng
      call iamas_log_write('--- Timing for exchanging halo u_field: $r s', &
           realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng


      !
      ! Perform basic compatibility checks among the fields that were read and the run-time options that were selected
      !
      call iamas_atm_run_compatibility(domain % dminfo, domain % blocklist, domain % streamManager, ierr)
      if (ierr /= 0) then
          call iamas_log_write('Please correct issues with the model input fields and/or namelist.', messageType=iAMAS_LOG_ERR)
          return
      end if


      call iamas_timer_start('Blocks init') !jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      block => domain % blocklist
      do while (associated(block))
         call iamas_pool_get_subpool(block % structs, 'mesh', mesh)
         call iamas_pool_get_subpool(block % structs, 'state', state)

         call atm_iamas_init_block(domain % dminfo, domain % streamManager, block, mesh, dt)

         call iamas_pool_get_array(state, 'xtime', xtime, 1)
         xtime = startTimeStamp

         ! Initialize initial_time in second time level. We need to do this because initial state
         ! is read into time level 1, and if we write output from the set of state arrays that
         ! represent the original time level 2, the initial_time field will be invalid.
         call iamas_pool_get_array(state, 'initial_time', initial_time1, 1)
         call iamas_pool_get_array(state, 'initial_time', initial_time2, 2)
         initial_time2 = initial_time1

         block => block % next
      end do
      call iamas_dmpar_get_time(timer_stop_time) !jwfeng
      call iamas_timer_stop('Blocks init') !jwfeng
      call iamas_log_write('--- Timing for init blocks: $r s', &
           realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng

      call iamas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call iamas_pool_get_field(diag, 'pv_edge', pv_edge_field)
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_dmpar_exch_halo_field(pv_edge_field)
      call iamas_dmpar_get_time(timer_stop_time) !jwfeng
      call iamas_log_write('--- Timing for exchanging halo pv_edge: $r s', &
           realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng

      call iamas_pool_get_field(diag, 'ru', ru_field)
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_dmpar_exch_halo_field(ru_field)
      call iamas_dmpar_get_time(timer_stop_time) !jwfeng
      call iamas_log_write('--- Timing for exchanging halo ru_field: $r s', &
           realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng

      call iamas_pool_get_field(diag, 'rw', rw_field)
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_dmpar_exch_halo_field(rw_field)
      call iamas_dmpar_get_time(timer_stop_time) !jwfeng
      call iamas_log_write('--- Timing for exchanging halo rw_field: $r s', &
           realArgs=(/real(timer_stop_time - timer_start_time, kind=RKIND)/)) !jwfeng

      call iamas_atm_diag_setup(domain % streamManager, domain % blocklist % configs, &
                               domain % blocklist % structs, domain % clock, domain % dminfo)


   end function atm_core_init


   subroutine atm_simulation_clock_init(core_clock, configs, ierr)

      use iamas_timekeeping

      implicit none

      type (iAMAS_Clock_type), intent(inout) :: core_clock
      type (iamas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (iAMAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (iAMAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: local_err
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_run_duration
      character (len=StrKIND), pointer :: config_stop_time
      character (len=StrKIND) :: startTimeStamp


      ierr = 0

      call iamas_pool_get_config(configs, 'config_dt', config_dt)
      call iamas_pool_get_config(configs, 'config_start_time', config_start_time)
      call iamas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call iamas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call iamas_pool_get_config(configs, 'config_stop_time', config_stop_time)

      if(trim(config_start_time) == 'file') then
         open(22,file=trim(config_restart_timestamp_name),form='formatted',status='old')
         read(22,*) startTimeStamp
         close(22)
      else
        startTimeStamp = config_start_time
      end if
      call iamas_set_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=local_err)
      call iamas_set_timeInterval(timeStep, dt=config_dt, ierr=local_err)

      if (trim(config_run_duration) /= "none") then
         call iamas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
         call iamas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)

         if (trim(config_stop_time) /= "none") then
            call iamas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
            if(startTime + runduration /= stopTime) then
               call iamas_log_write('config_run_duration and config_stop_time are inconsistent: using config_run_duration.', messageType=iAMAS_LOG_WARN)
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call iamas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
         call iamas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=local_err)
      else
          call iamas_log_write('Neither config_run_duration nor config_stop_time were specified.', messageType=iAMAS_LOG_ERR)
          ierr = 1
      end if

      !TODO: set phyics alarms here...
      !....
      !....

   end subroutine atm_simulation_clock_init


   subroutine atm_iamas_init_block(dminfo, stream_manager, block, mesh, dt)
   
      use atm_time_integration
      use iamas_rbf_interpolation
      use iamas_vector_reconstruction
      use iamas_stream_manager
      use iamas_atm_boundaries, only : iamas_atm_setup_bdy_masks
#ifdef DO_PHYSICS
!     use iamas_atmphys_aquaplanet
      use iamas_atmphys_control
      use iamas_atmphys_init
      use iamas_atmphys_manager
#endif

#ifdef DO_CHEMISTRY
      use iamas_atmchem_vars
      use iamas_atmchem_gas_init
      use iamas_atmchem_binaer_init
      use module_dep_simple
#endif

      implicit none
   
      type (dm_info), intent(in) :: dminfo
      type (iAMAS_streamManager_type), intent(inout) :: stream_manager
      type (block_type), intent(inout) :: block
      type (iamas_pool_type), intent(inout) :: mesh     !MGD does this need to be a pointer?
      real (kind=RKIND), intent(in) :: dt

      type (iamas_pool_type), pointer :: state
      type (iamas_pool_type), pointer :: diag
      type (iamas_pool_type), pointer :: tend
      type (iamas_pool_type), pointer :: sfc_input
      type (iamas_pool_type), pointer :: diag_physics
#ifdef DO_CHEMISTRY
      type (iamas_pool_type), pointer :: diag_chemistry, tend_chemistry !zhxia add this for photolysis
      integer :: itimestep
#endif
      type (iamas_pool_type), pointer :: atm_input

      integer :: iCell,iEdge,iVertex
      
      real (kind=RKIND), dimension(:,:), pointer :: u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:), pointer :: meshScalingRegionalCell, meshScalingRegionalEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invDvEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge, invDcEdge
      real (kind=RKIND), dimension(:), pointer :: areaTriangle, invAreaTriangle
      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      integer :: thread
      character(len=StrKIND), pointer :: mminlu

      integer, pointer :: nThreads
      integer, dimension(:), pointer :: cellThreadStart, cellThreadEnd
      integer, dimension(:), pointer :: cellSolveThreadStart, cellSolveThreadEnd
      integer, dimension(:), pointer :: edgeThreadStart, edgeThreadEnd
      integer, dimension(:), pointer :: edgeSolveThreadStart, edgeSolveThreadEnd
      integer, dimension(:), pointer :: vertexThreadStart, vertexThreadEnd
      integer, dimension(:), pointer :: vertexSolveThreadStart, vertexSolveThreadEnd


      logical, pointer :: config_do_restart, config_do_DAcycling
      character(len=StrKIND), pointer :: config_chem_opt
      logical, pointer :: config_use_fix_gasinit !determine whether we use the photochemistry simulation results from wrfchem as initial values in iAMAS. zhxia@lacar
      logical, pointer :: config_input_chem !determine whether we use the specific values of chem array in iAMAS. zhxia@lacar

      call atm_compute_signs(mesh)
   
      call iamas_pool_get_subpool(block % structs, 'diag', diag)
      call iamas_pool_get_subpool(block % structs, 'state', state)

      call iamas_pool_get_subpool(block % structs, 'state', state)

      call iamas_pool_get_config(block % configs, 'config_chem_opt', config_chem_opt)
      call iamas_pool_get_config(block % configs, 'config_use_fix_gasinit', config_use_fix_gasinit)
      call iamas_pool_get_config(block % configs, 'config_input_chem', config_input_chem)
      call iamas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)
      call iamas_pool_get_config(block % configs, 'config_do_DAcycling', config_do_DAcycling)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!! Compute inverses to avoid divides later
      
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)      
      call iamas_pool_get_array(mesh, 'invDvEdge', invDvEdge)
      call iamas_pool_get_array(mesh, 'invDcEdge', invDcEdge)      
      
      call iamas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call iamas_pool_get_array(mesh, 'invAreaTriangle', invAreaTriangle)
      
      call iamas_pool_get_dimension(mesh, 'nCells', nCells)      
      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
 
      do iCell=1,nCells 
         invAreaCell(iCell) = 1.0_RKIND / areaCell(iCell)
      end do
        
      do iEdge=1,nEdges 
         invDvEdge(iEdge) = 1.0_RKIND / dvEdge(iEdge)
      end do
      
      do iEdge=1,nEdges 
         invDcEdge(iEdge) = 1.0_RKIND / dcEdge(iEdge)
      end do
      
      do iVertex=1,nVertices
         invAreaTriangle(iVertex) = 1.0_RKIND / areaTriangle(iVertex)
      end do
      
      !!!!! End compute inverses
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call atm_adv_coef_compression(mesh)

      call atm_couple_coef_3rd_order(mesh, block % configs)
      
      call iamas_pool_get_dimension(state, 'nVertices', nVertices)
      call iamas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

      allocate(ke_vertex(nVertLevels,nVertices+1))  ! ke_vertex is a module variable defined in iamas_atm_time_integration.F
      ke_vertex(:,nVertices+1) = 0.0_RKIND
      allocate(ke_edge(nVertLevels,nEdges+1))       ! ke_edge is a module variable defined in iamas_atm_time_integration.F
      ke_edge(:,nEdges+1) = 0.0_RKIND

      call iamas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

      call iamas_pool_get_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
      call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

      call iamas_pool_get_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
      call iamas_pool_get_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

      call iamas_pool_get_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
      call iamas_pool_get_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
      call iamas_pool_get_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

!$OMP PARALLEL DO
      do thread=1,nThreads
         if (.not. config_do_restart .or. (config_do_restart .and. config_do_DAcycling)) then
            call atm_init_coupled_diagnostics(state, 1, diag, mesh, block % configs, &
                                              cellThreadStart(thread), cellThreadEnd(thread), &
                                              vertexThreadStart(thread), vertexThreadEnd(thread), &
                                              edgeThreadStart(thread), edgeThreadEnd(thread), &
                                              cellSolveThreadStart(thread), cellSolveThreadEnd(thread), &
                                              vertexSolveThreadStart(thread), vertexSolveThreadEnd(thread), &
                                              edgeSolveThreadStart(thread), edgeSolveThreadEnd(thread))
         end if

         call atm_compute_solve_diagnostics(dt, state, 1, diag, mesh, block % configs, &
                                            cellThreadStart(thread), cellThreadEnd(thread), &
                                            vertexThreadStart(thread), vertexThreadEnd(thread), &
                                            edgeThreadStart(thread), edgeThreadEnd(thread))
      end do
!$OMP END PARALLEL DO

      deallocate(ke_vertex)
      deallocate(ke_edge)

      call iamas_rbf_interp_initialize(mesh)
      call iamas_init_reconstruct(mesh)

      call iamas_pool_get_array(state, 'u', u, 1)
      call iamas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call iamas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call iamas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      call iamas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call iamas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call iamas_reconstruct(mesh, u,                   &
                            uReconstructX,             &
                            uReconstructY,             &
                            uReconstructZ,             &
                            uReconstructZonal,         &
                            uReconstructMeridional     &
                           )
   
#ifdef DO_PHYSICS
      !proceed with initialization of physics parameterization if moist_physics is set to true:
      call iamas_pool_get_subpool(block % structs, 'sfc_input', sfc_input)

      ! Before calling physics_init, ensure that mminlu contains the name of the land use dataset
      call iamas_pool_get_array(sfc_input, 'mminlu', mminlu)
      if (len_trim(mminlu) == 0) then
            call iamas_log_write('****************************************************************')
            call iamas_log_write('No information on land use dataset is available.')
            call iamas_log_write('Assume that we are using ''USGS''.')
            call iamas_log_write('****************************************************************')
            write(mminlu,'(a)') 'USGS'
      end if


      if (moist_physics) then
         !initialization of some input variables in registry:
         call iamas_pool_get_subpool(block % structs, 'tend', tend)
         call iamas_pool_get_subpool(block % structs, 'diag_physics', diag_physics)
         call iamas_pool_get_subpool(block % structs, 'atm_input', atm_input)
         call physics_tables_init(dminfo, block % configs)
         call physics_registry_init(mesh, block % configs, sfc_input)
         call physics_run_init(block % configs, mesh, state, clock, stream_manager)

         !initialization of all physics:
         call physics_init(dminfo, clock, block % configs, mesh, diag, tend, state, 1, diag_physics, &
                           atm_input, sfc_input)
      endif
#endif
!jwfeng: 03/21/2020
#ifdef DO_CHEMISTRY
      !initialization of all chemistry:
      if (trim(config_chem_opt) .ne. "no_chemistry") then
          call iamas_pool_get_subpool(block % structs, 'diag_chemistry', diag_chemistry)
          call iamas_pool_get_subpool(block % structs, 'tend_chemistry', tend_chemistry) !zhxia
          call atmchem_vars_init(block%configs, mesh, state, diag_chemistry, tend_chemistry, clock) !zhxia
          call iamas_pool_get_array(sfc_input, 'mminlu', mminlu)
          itimestep = 1
          call init_data_binaer_asect( state, tend_chemistry, trim(config_chem_opt), is_aerosol )
          if (config_input_chem .or. config_use_fix_gasinit) then
             call iamas_log_write("gas init using specific values of chem array") !zhxia
          else 
             call iamas_log_write("gas init using default profile")
             call gas_init( block%configs, mesh ) !zhxia
          endif
          call dep_init( state,id = itimestep , mminlu_loc = trim(mminlu))
      end if
#else
      call iamas_log_write("**************************" )
      call iamas_log_write("*** Don't do chemistry ***" )
      call iamas_log_write("**************************" )
#endif

      call atm_compute_mesh_scaling(mesh, block % configs)

      call atm_compute_damping_coefs(mesh, block % configs)

      call atm_compute_vert_implicit_coefs(mesh, block % configs)
      !
      ! Set up mask fields used in limited-area simulations
      !
      call iamas_atm_setup_bdy_masks(mesh, block % configs)

      call iamas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call iamas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call iamas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)

      call iamas_log_write('min/max of meshScalingDel2 = $r $r', &
                          realArgs=(/minval(meshScalingDel2(1:nEdgesSolve)), maxval(meshScalingDel2(1:nEdgesSolve))/))
      call iamas_log_write('min/max of meshScalingDel4 = $r $r', &
                          realArgs=(/minval(meshScalingDel4(1:nEdgesSolve)), maxval(meshScalingDel4(1:nEdgesSolve))/))

   end subroutine atm_iamas_init_block
   
   
   function atm_core_run(domain) result(ierr)
   
      use iamas_timekeeping
      use iamas_kind_types
      use iamas_stream_manager
      use iamas_derived_types, only : iAMAS_STREAM_LATEST_BEFORE, iAMAS_STREAM_INPUT, iAMAS_STREAM_INPUT_OUTPUT
      use iamas_timer
      use iamas_atm_boundaries, only : iamas_atm_update_bdy_tend
      use iamas_atm_diagnostics_manager, only : iamas_atm_diag_update, iamas_atm_diag_compute, iamas_atm_diag_reset
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: ierr
   
      real (kind=RKIND), pointer :: dt
      logical, pointer :: config_do_restart
      logical, pointer :: config_apply_lbcs
      logical, pointer :: config_emis_update
      type (block_type), pointer :: block_ptr

      type (iAMAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_fdda_scheme
      integer :: itimestep

      integer :: stream_dir
      character(len=StrKIND) :: input_stream, read_time

      type (iamas_pool_type), pointer :: state, diag, mesh, diag_physics, tend, tend_physics

      ! For high-frequency diagnostics output
      character (len=StrKIND) :: tempfilename

      ! For timing information
      real (kind=R8KIND) :: integ_start_time, integ_stop_time
      real (kind=R8KIND) :: diag_start_time, diag_stop_time
      real (kind=R8KIND) :: input_start_time, input_stop_time
      real (kind=R8KIND) :: output_start_time, output_stop_time
      
      ierr = 0

      clock => domain % clock
      iamas_log_info => domain % logInfo
      
      ! Eventually, dt should be domain specific
      call iamas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Avoid writing a restart file at the initial time
      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

      ! Also, for restart runs, avoid writing the initial history or diagnostics fields to avoid overwriting those from the preceding run
      if (config_do_restart) then
         call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='output', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)
         call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='diagnostics', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)
      end if

      call iamas_dmpar_get_time(diag_start_time)

      if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call iamas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
            call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call atm_compute_output_diagnostics(state, 1, diag, mesh)

            block_ptr => block_ptr % next
         end do
      end if
      call iamas_atm_diag_reset()
      call iamas_atm_diag_update()
      call iamas_atm_diag_compute()

      call iamas_dmpar_get_time(diag_stop_time)

      call iamas_dmpar_get_time(output_start_time)
      call iamas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call iamas_dmpar_get_time(output_stop_time)
      if (ierr /= iAMAS_STREAM_MGR_NOERR .and. &
          ierr /= iAMAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
          ierr /= iAMAS_STREAM_MGR_ERR_CLOBBER_REC) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('Error writing one or more output streams',                                         messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
      if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
         call iamas_log_write('--- Timing for diagnostic computation: $r s', realArgs=(/real(diag_stop_time - diag_start_time, kind=RKIND)/))
         call iamas_log_write('--- Timing for stream output: $r s', realArgs=(/real(output_stop_time - output_start_time, kind=RKIND)/))
      end if

      call iamas_atm_diag_reset()

      call iamas_stream_mgr_reset_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

      block_ptr => domain % blocklist
      call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
      call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
      call iamas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)

      call iamas_pool_get_config(domain % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_emis_update', config_emis_update)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_fdda_scheme', config_fdda_scheme)

      !
      ! Read initial boundary state
      !
      if (config_apply_lbcs .and. &
          iAMAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='lbc_in', direction=iAMAS_STREAM_INPUT, ierr=ierr)) then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call iamas_atm_update_bdy_tend(clock, domain % streamManager, block_ptr, .true., ierr)
            if (ierr /= 0) then
               currTime = iamas_get_clock_time(clock, iAMAS_NOW, ierr)
               call iamas_get_time(curr_time=currTime, dateTimeString=timeStamp)
               call iamas_log_write('Failed to process LBC data on or before '//trim(timeStamp), messageType=iAMAS_LOG_ERR)
               return
            end if

            block_ptr => block_ptr % next
         end do
      end if

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 1
      do while (.not. iamas_is_clock_stop_time(clock))

         currTime = iamas_get_clock_time(clock, iAMAS_NOW, ierr)
         call iamas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         call iamas_log_write('')
         call iamas_log_write('Begin timestep '//trim(timeStamp))

         !
         ! Read future boundary state and compute boundary tendencies
         !
         if (config_apply_lbcs .and. &
             iAMAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='lbc_in', direction=iAMAS_STREAM_INPUT, ierr=ierr)) then
            block_ptr => domain % blocklist
            do while (associated(block_ptr))
               call iamas_atm_update_bdy_tend(clock, domain % streamManager, block_ptr, .false., ierr)
               if (ierr /= 0) then
                  call iamas_log_write('Failed to process LBC data at next time after '//trim(timeStamp), messageType=iAMAS_LOG_ERR)
                  return
               end if

               block_ptr => block_ptr % next
            end do
         end if

         ! Regardless of whether boundary tendencies were updated, above, we do not want to read the 'lbc_in' stream
         ! as a general input stream, below.
         call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='lbc_in', direction=iAMAS_STREAM_INPUT, ierr=ierr)


         ! Don't read emission_update file if we set config_emis_update=false. @zhxia 0420/2025
         if (.not. config_emis_update) then
             call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='emission', direction=iAMAS_STREAM_INPUT, ierr=ierr)
         endif 

         ! Don't read fdda file if we set config_fdda_scheme="off". @jwfeng 1218/2025
         if (trim(config_fdda_scheme) .eq. "off") then
             call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='fdda', direction=iAMAS_STREAM_INPUT, ierr=ierr)
         endif 

         !
         ! Read external field updates
         !
         call iAMAS_stream_mgr_begin_iteration(domain % streamManager, ierr=ierr)
         do while (iAMAS_stream_mgr_get_next_stream(domain % streamManager, streamID=input_stream, directionProperty=stream_dir))
            if (stream_dir == iAMAS_STREAM_INPUT .or. stream_dir == iAMAS_STREAM_INPUT_OUTPUT) then
               if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, streamID=input_stream, &
                                                  direction=iAMAS_STREAM_INPUT, ierr=ierr)) then
                  call iamas_dmpar_get_time(input_start_time)
                  call iAMAS_stream_mgr_read(domain % streamManager, streamID=input_stream, whence=iAMAS_STREAM_LATEST_BEFORE, &
                                            actualWhen=read_time, ierr=ierr)
                  call iamas_dmpar_get_time(input_stop_time)
                  if (ierr /= iAMAS_STREAM_MGR_NOERR) then
                     call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
                     call iamas_log_write('Error reading input stream '//trim(input_stream),                                  messageType=iAMAS_LOG_ERR)
                     call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
                  end if

                  call iamas_log_write('----------------------------------------------------------------------')
                  call iamas_log_write('  Read '''//trim(input_stream)//''' input stream valid at '//trim(read_time))
                  call iamas_log_write('--- Timing for stream input: $r s', realArgs=(/real(input_stop_time - input_start_time, kind=RKIND)/))
                  call iamas_log_write('----------------------------------------------------------------------')

                  call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID=input_stream, direction=iAMAS_STREAM_INPUT, ierr=ierr)
               end if
            end if
         end do

         call iamas_timer_start("time integration")
         call iamas_dmpar_get_time(integ_start_time)
         call atm_do_timestep(domain, dt, itimestep)
         call iamas_dmpar_get_time(integ_stop_time)
         call iamas_timer_stop("time integration")   
         call iamas_log_write('--- Timing for integration step: $r s', realArgs=(/real(integ_stop_time - integ_start_time, kind=RKIND)/))

         ! Move time level 2 fields back into time level 1 for next time step
         call iamas_pool_get_subpool(domain % blocklist % structs, 'state', state)
         call iamas_pool_shift_time_levels(state)
         
         ! Advance clock before writing output
         itimestep = itimestep + 1
         call iamas_advance_clock(clock)
         currTime = iamas_get_clock_time(clock, iAMAS_NOW, ierr)


         call iamas_dmpar_get_time(diag_start_time)

         !
         ! Write any output streams that have alarms ringing, after computing diagnostics fields
         !
         call iamas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
           block_ptr => domain % blocklist
           do while (associated(block_ptr))
              call iamas_pool_get_subpool(block_ptr % structs, 'state', state)
              call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
              call iamas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
              call iamas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
              call iamas_pool_get_subpool(block_ptr % structs, 'tend', tend)
              call iamas_pool_get_subpool(block_ptr % structs, 'tend_physics', tend_physics)
              call atm_compute_output_diagnostics(state, 1, diag, mesh)

              block_ptr => block_ptr % next
           end do
         end if

         call iamas_atm_diag_update()
         call iamas_atm_diag_compute()
         call iamas_dmpar_get_time(diag_stop_time)

         call iamas_timer_start("write output files")
         call iamas_dmpar_get_time(output_start_time)
         call iamas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call iamas_dmpar_get_time(output_stop_time)
		 call iamas_timer_stop("write output files")

         if (ierr /= iAMAS_STREAM_MGR_NOERR .and. &
             ierr /= iAMAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
             ierr /= iAMAS_STREAM_MGR_ERR_CLOBBER_REC) then
            call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
            call iamas_log_write('Error writing one or more output streams',                                         messageType=iAMAS_LOG_ERR)
            call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
         end if
         if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
            call iamas_log_write('--- Timing for diagnostic computation: $r s', realArgs=(/real(diag_stop_time - diag_start_time, kind=RKIND)/))
            call iamas_log_write('--- Timing for stream output: $r s', realArgs=(/real(output_stop_time - output_start_time, kind=RKIND)/))
         end if

         ! reset any diagnostics here

         if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='diagnostics', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
            block_ptr => domain % blocklist
            do while (associated(block_ptr))

               call iamas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call iamas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
               call atm_reset_diagnostics(diag, diag_physics)

               block_ptr => block_ptr % next
            end do
         end if


         ! Only after we've successfully written the restart file should we we
         !    write the restart_timestamp file
         if (iAMAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then
            if (domain % dminfo % my_proc_id == 0) then
               open(22,file=trim(config_restart_timestamp_name),form='formatted',status='replace')
               write(22,*) trim(timeStamp)
               close(22)
            end if
         end if

         call iamas_atm_diag_reset()

         call iamas_stream_mgr_reset_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

      end do
   
   end function atm_core_run

   
   subroutine atm_compute_output_diagnostics(state, time_lev, diag, mesh)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to history files
   !
   ! Input: state - contains model prognostic fields
   !        mesh  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use iamas_constants
   
      implicit none
   
      type (iamas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(in) :: mesh
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND), parameter :: p_ref = 100000.0_RKIND
      real (kind=RKIND), dimension(:,:), pointer :: theta, temperature, rho, theta_m, rho_zz, zz
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base, pressure_p, pressure
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call iamas_pool_get_dimension(state, 'index_qv', index_qv)

      call iamas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call iamas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
      call iamas_pool_get_array(state, 'scalars', scalars, time_lev)

      call iamas_pool_get_array(diag, 'theta', theta)
      call iamas_pool_get_array(diag, 'temperature', temperature)
      call iamas_pool_get_array(diag, 'rho', rho)
      call iamas_pool_get_array(diag, 'pressure_p', pressure_p)
      call iamas_pool_get_array(diag, 'pressure_base', pressure_base)
      call iamas_pool_get_array(diag, 'pressure', pressure)

      call iamas_pool_get_array(mesh, 'zz', zz)

      do iCell=1,nCells
         do k=1,nVertLevels
            theta(k,iCell) = theta_m(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            pressure(k,iCell) = pressure_base(k,iCell) + pressure_p(k,iCell)
            temperature(k,iCell) = theta(k,iCell) * (pressure(k,iCell) / p_ref)**(rgas / cp)
         end do
      end do

   end subroutine atm_compute_output_diagnostics
   
   
   subroutine atm_reset_diagnostics(diag, diag_physics)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! reset some diagnostics after output
   !
   ! Input: diag          - contains dynamics diagnostic fields
   !        daig_physics  - contains physics diagnostic fields
   !
   ! Output: whatever diagnostics need resetting after output
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none
   
      type (iamas_pool_type), intent(inout) :: diag
      type (iamas_pool_type), intent(inout) :: diag_physics
   
      real (kind=RKIND), dimension(:), pointer :: refl10cm_1km_max

      call iamas_pool_get_array(diag_physics, 'refl10cm_1km_max', refl10cm_1km_max)
      if(associated(refl10cm_1km_max)) then
         refl10cm_1km_max(:) = 0.
      endif
   
   end subroutine atm_reset_diagnostics


   subroutine atm_do_timestep(domain, dt, itimestep)
   
      use iamas_timekeeping
      use iamas_kind_types
      use atm_time_integration

#ifdef DO_CHEMISTRY
      use iamas_atmchem_driver
#endif

#ifdef DO_PHYSICS
      use iamas_atmphys_control
      use iamas_atmphys_driver
      use iamas_atmphys_manager
      use iamas_atmphys_update
#endif
  
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: itimestep
      
      type (iAMAS_Time_Type) :: startTime, currTime
      type (iAMAS_TimeInterval_Type) :: xtimeTime
      integer :: s, s_n, s_d, dd
      real (kind=RKIND) :: xtime_s, xtime_d, curr_secs !zhxia add curr_secs for total simulation seconds from start time.
      integer :: ierr
      character(len=StrKIND), pointer :: config_chem_opt

      clock => domain % clock
      iamas_log_info => domain % logInfo

      startTime = iamas_get_clock_time(clock, iAMAS_START_TIME, ierr)
      currTime = iamas_get_clock_time(clock, iAMAS_NOW, ierr)
         
      xtimeTime = currTime - startTime
      call iamas_get_timeInterval(interval=xtimeTime, DD=dd, S=s, S_n=s_n, S_d=s_d, dt=curr_secs, ierr=ierr)         
      xtime_s = (s + s_n / s_d) !total simulation seconds - total simulation days * 60 * 60 *24
      xtime_d = dd !total simulation days, eg. xtimeTime=21d 6h 3min 5s, so xtime_d=21, xtime_s=6*3600+3*60+5


#ifdef DO_PHYSICS
      !proceed with physics if moist_physics is set to true:
      if(moist_physics) then
         call physics_timetracker(domain,dt,clock,itimestep,xtime_s)
         call physics_driver(domain,itimestep,xtime_s)

      endif
#endif

!AEMOL USTC 02/2020
!jwfeng: You can add CHEMISTRY=-DDO_CHEMISTRY to activate chemistry_driver 
#ifdef DO_CHEMISTRY
      call iamas_pool_get_config(domain % blocklist % configs, 'config_chem_opt', config_chem_opt)
      if (trim(config_chem_opt) .ne. "no_chemistry") then
         call chemistry_driver(domain,itimestep,xtime_d,xtime_s,curr_secs)
      end if
#endif

      call atm_timestep(domain, dt, currTime, itimestep)

!jwfeng diag
#ifdef DO_CHEMISTRY
      if (trim(config_chem_opt) .ne. "no_chemistry") then
         call driver_chem_diag(domain)
      end if
#endif

   end subroutine atm_do_timestep
   
   
   function atm_core_finalize(domain) result(ierr)
   
      use iamas_decomp
      use iamas_timekeeping
      use iamas_atm_diagnostics_manager, only : iamas_atm_diag_cleanup
      use iamas_atm_threading, only : iamas_atm_threading_finalize
   
#ifdef DO_PHYSICS
      use iamas_atmphys_finalize
#endif

      use flcl_util_kokkos_mod

      implicit none

#ifdef DO_PHYSICS
   interface
      subroutine rrtmg_lw_workspace_finalize() bind(C,name="rrtmg_lw_workspace_finalize")
         use iso_c_binding
      end subroutine rrtmg_lw_workspace_finalize
      subroutine rrtmg_sw_workspace_finalize() bind(C,name="rrtmg_sw_workspace_finalize")
         use iso_c_binding
      end subroutine rrtmg_sw_workspace_finalize
   end interface
#endif

#ifdef DO_CHEMISTRY
   interface
      subroutine mixactivate_workspace_finalize() bind(C,name="mixactivate_workspace_finalize")
         use iso_c_binding
      end subroutine mixactivate_workspace_finalize
      subroutine mosaic_settling_workspace_finalize() bind(C,name="mosaic_settling_workspace_finalize")
         use iso_c_binding
      end subroutine mosaic_settling_workspace_finalize
      subroutine optical_averaging_workspace_finalize() bind(C,name="optical_averaging_workspace_finalize")
         use iso_c_binding
      end subroutine optical_averaging_workspace_finalize
      subroutine ctrans_general_workspace_finalize() bind(C,name="ctrans_general_workspace_finalize")
         use iso_c_binding
      end subroutine ctrans_general_workspace_finalize
   end interface
#endif

      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr 
      integer :: ierr

      ierr = 0

      clock => domain % clock
      iamas_log_info => domain % logInfo


      call iamas_atm_diag_cleanup()

      call iamas_destroy_clock(clock, ierr)
      call iamas_decomp_destroy_decomp_list(domain % decompositions)

#ifdef DO_PHYSICS
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call atmphys_finalize(block_ptr%configs)

            block_ptr => block_ptr%next
         end do
#endif
#ifdef DO_PHYSICS
      call rrtmg_lw_workspace_finalize()
      call rrtmg_sw_workspace_finalize()
#endif
#ifdef DO_CHEMISTRY
      call mixactivate_workspace_finalize()
      call mosaic_settling_workspace_finalize()
      call optical_averaging_workspace_finalize()
      call ctrans_general_workspace_finalize()
#endif
      !
      ! Finalize threading
      !
      call iamas_atm_threading_finalize(domain % blocklist)

      call iamas_log_write('')
      call iamas_log_write('********************************************************')
      call iamas_log_write('   Finished running the atmosphere core')
      call iamas_log_write('********************************************************')
      call kokkos_finalize()
   
   end function atm_core_finalize


   subroutine atm_compute_mesh_scaling(mesh, configs)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      type (iamas_pool_type), intent(in) :: configs

      integer :: iCell,iEdge, cell1, cell2
      integer, pointer :: nEdges, nCells
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:), pointer :: meshScalingRegionalCell, meshScalingRegionalEdge
      logical, pointer :: config_h_ScaleWithMesh

      call iamas_pool_get_array(mesh, 'meshDensity', meshDensity)
      call iamas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call iamas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)
      call iamas_pool_get_array(mesh, 'meshScalingRegionalCell', meshScalingRegionalCell)
      call iamas_pool_get_array(mesh, 'meshScalingRegionalEdge', meshScalingRegionalEdge)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'nCells', nCells)

      call iamas_pool_get_config(configs, 'config_h_ScaleWithMesh', config_h_ScaleWithMesh)

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.25
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.75
         end do
      end if

      !
      ! Compute the scaling factors to be used in relaxation zone of regional configuration
      !
      meshScalingRegionalCell(:) = 1.0
      meshScalingRegionalEdge(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            meshScalingRegionalEdge(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.25
         end do

         do iCell=1,nCells
            meshScalingRegionalCell(iCell) = 1.0 / (meshDensity(iCell))**0.25
         end do
      end if

   end subroutine atm_compute_mesh_scaling


   subroutine atm_compute_signs(mesh)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh

      integer :: i, j, iCell, iVtx
      integer, pointer :: nCells, nVertices, nEdges, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnVertex, verticesOnEdge, cellsOnEdge, edgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, kiteForCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnVertex_sign, edgesOnCell_sign
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell


      call iamas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call iamas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call iamas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call iamas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)
      call iamas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call iamas_pool_get_array(mesh, 'kiteForCell', kiteForCell)
      call iamas_pool_get_array(mesh, 'zb', zb)
      call iamas_pool_get_array(mesh, 'zb3', zb3)
      call iamas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call iamas_pool_get_array(mesh, 'zb3_cell', zb3_cell)

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call iamas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)


      do iVtx=1,nVertices
         do i=1,vertexDegree
            if (edgesOnVertex(i,iVtx) <= nEdges) then
               if (iVtx == verticesOnEdge(2,edgesOnVertex(i,iVtx))) then
                  edgesOnVertex_sign(i,iVtx) = 1.0
               else
                  edgesOnVertex_sign(i,iVtx) = -1.0
               end if
            else
               edgesOnVertex_sign(i,iVtx) = 0.0
            end if
         end do
      end do

      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            if (edgesOnCell(i,iCell) <= nEdges) then
               if (iCell == cellsOnEdge(1,edgesOnCell(i,iCell))) then
                  edgesOnCell_sign(i,iCell) = 1.0
                  zb_cell(:,i,iCell) = zb(:,1,edgesOnCell(i,iCell))
                  zb3_cell(:,i,iCell) = zb3(:,1,edgesOnCell(i,iCell))
               else
                  edgesOnCell_sign(i,iCell) = -1.0
                  zb_cell(:,i,iCell) = zb(:,2,edgesOnCell(i,iCell))
                  zb3_cell(:,i,iCell) = zb3(:,2,edgesOnCell(i,iCell))
               end if
            else
               edgesOnCell_sign(i,iCell) = 0.0
            end if
         end do
      end do

      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iVtx = verticesOnCell(i,iCell)
            if (iVtx <= nVertices) then
               do j=1,vertexDegree
                  if (iCell == cellsOnVertex(j,iVtx)) then
                     kiteForCell(i,iCell) = j
                     exit
                  end if
               end do
               if (j > vertexDegree) then
                  call iamas_log_write('Unexpected error while identifying kiteForCell', messageType=iAMAS_LOG_ERR)
               end if
            else
               kiteForCell(i,iCell) = 1
            end if
         end do
      end do

   end subroutine atm_compute_signs


   subroutine atm_compute_damping_coefs(mesh, configs)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      type (iamas_pool_type), intent(in) :: configs

      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      real (kind=RKIND), pointer :: config_xnutr, config_zd
      real (kind=RKIND) :: z, zt, m1, pii
      real (kind=RKIND), dimension(:,:), pointer :: dss, zgrid
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      real (kind=RKIND) :: dx_scale_power

      m1 = -1.0
      pii = acos(m1)

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      call iamas_pool_get_array(mesh, 'meshDensity', meshDensity)
      call iamas_pool_get_array(mesh, 'dss', dss)
      call iamas_pool_get_array(mesh, 'zgrid', zgrid)

      call iamas_pool_get_config(configs, 'config_zd', config_zd)
      call iamas_pool_get_config(configs, 'config_xnutr', config_xnutr)

      dx_scale_power = 1.0
      dss(:,:) = 0.0
      do iCell=1,nCells
         zt = zgrid(nVertLevels+1,iCell)
         do k=1,nVertLevels
            z = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
            if (z > config_zd) then
               dss(k,iCell) = config_xnutr*sin(0.5*pii*(z-config_zd)/(zt-config_zd))**2.0
               dss(k,iCell) = dss(k,iCell) / meshDensity(iCell)**(0.25*dx_scale_power)
            end if
         end do
      end do

   end subroutine atm_compute_damping_coefs

   ! jwfeng calculating adjusted epssm
   subroutine atm_compute_vert_implicit_coefs(mesh, configs)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      type (iamas_pool_type), intent(in) :: configs

      ! local variables
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, epssm_adj
      real (kind=RKIND), pointer :: epssm, zb, zt ! for adjustment of epssm
      integer, pointer :: nCells, nVertLevels

      real :: z, pii
      integer :: iCell, k

      pii = acos(-1.0)

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'epssm_adj', epssm_adj)

      call iamas_pool_get_config(configs, 'config_epssm', epssm)
      call iamas_pool_get_config(configs, 'config_epssm_zb', zb)
      call iamas_pool_get_config(configs, 'config_epssm_zt', zt)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            z = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
            ! jwfeng: epssm_adj >= epssm, <= 0.5
            if (z > zt) then
               epssm_adj(k,iCell) = 0.5
            else if (z < zb) then
               epssm_adj(k,iCell) = epssm
            else
               epssm_adj(k,iCell) = epssm + (0.5-epssm)*sin(0.5*pii*(z-zb)/(zt-zb))**2
            end if
         end do
      end do

   end subroutine atm_compute_vert_implicit_coefs


   subroutine atm_adv_coef_compression( mesh )

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh


      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge

      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell
      integer, pointer :: nCells, nEdges
      integer :: cell_list(20)
      logical :: addcell


      call iamas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call iamas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call iamas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)
      call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nEdges', nEdges)

      do iEdge = 1, nEdges
         nAdvCellsForEdge(iEdge) = 0
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         !
         ! do only if this edge flux is needed to update owned cells
         !
         if (cell1 <= nCells .or. cell2 <= nCells) then
 
            cell_list(1) = cell1
            cell_list(2) = cell2
            n = 2 
  
          !  add cells surrounding cell 1.  n is number of cells currently in list
            do i = 1, nEdgesOnCell(cell1)
               if (cellsOnCell(i,cell1) /= cell2) then
                  n = n + 1
                  cell_list(n) = cellsOnCell(i,cell1)
               end if
            end do
  
          !  add cells surrounding cell 2 (brute force approach)
            do iCell = 1, nEdgesOnCell(cell2)
               addcell = .true.
               do i=1,n
                  if (cell_list(i) == cellsOnCell(iCell,cell2)) addcell = .false.
               end do
               if (addcell) then
                  n = n+1
                  cell_list(n) = cellsOnCell(iCell,cell2)
               end if
            end do
  
  
            nAdvCellsForEdge(iEdge) = n
            do iCell = 1, nAdvCellsForEdge(iEdge)
               advCellsForEdge(iCell,iEdge) = cell_list(iCell)
            end do
  
          ! we have the ordered list, now construct coefficients
  
            adv_coefs(:,iEdge) = 0.
            adv_coefs_3rd(:,iEdge) = 0.
          
          ! pull together third and fourth order contributions to the flux
          ! first from cell1
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,1,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(1,1,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell1)
               j_in = 0
               do j=1, n
                 if( cell_list(j) == cellsOnCell(iCell,cell1) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
            end do
  
          ! pull together third and fourth order contributions to the flux
          ! now from cell2
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,2,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(1,2,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell2)
               j_in = 0
               do j=1, n
                  if( cell_list(j) == cellsOnCell(iCell,cell2) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,2,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(iCell+1,2,iEdge)
            end do
  
            do j = 1,n
               adv_coefs    (j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs    (j,iEdge) / 12.
               adv_coefs_3rd(j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end do
  
          ! 2nd order centered contribution - place this in the main flux weights
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
          !  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply
  
            do j=1,n
               adv_coefs    (j,iEdge) = dvEdge(iEdge) * adv_coefs    (j,iEdge)
               adv_coefs_3rd(j,iEdge) = dvEdge(iEdge) * adv_coefs_3rd(j,iEdge)
            end do
 
         end if  ! only do for edges of owned-cells
         
      end do ! end loop over edges

   end subroutine atm_adv_coef_compression


   subroutine atm_couple_coef_3rd_order(mesh, configs)

      implicit none

      type (iamas_pool_type), intent(inout) :: mesh
      type (iamas_pool_type), intent(in) :: configs

      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs_3rd
      real (kind=RKIND), dimension(:,:,:), pointer :: zb3_cell
      real (kind=RKIND), pointer :: config_coef_3rd_order

      call iamas_pool_get_array(mesh, 'zb3_cell', zb3_cell)
      call iamas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)

      adv_coefs_3rd(:,:) = config_coef_3rd_order * adv_coefs_3rd(:,:)
      zb3_cell(:,:,:) = config_coef_3rd_order * zb3_cell(:,:,:)

   end subroutine atm_couple_coef_3rd_order


   !-----------------------------------------------------------------------
   !  routine iamas_atm_run_compatibility
   !
   !> \brief Checks input fields and options for compatibility
   !> \author Michael Duda
   !> \date   22 October 2018
   !> \details
   !>  This routine checks the input fields and run-time options provided
   !>  by the user for compatibility. For example, two run-time options may
   !>  be mutually exclusive, or an option may require that a certain input
   !>  field is provided.
   !>
   !>  A value of 0 is returned if there are no incompatibilities among
   !>  the provided input fields and run-time options, and a non-zero value
   !>  otherwise.
   !
   !-----------------------------------------------------------------------
   subroutine iamas_atm_run_compatibility(dminfo, blockList, streamManager, ierr)

      use iamas_atmphys_control, only : physics_compatibility_check
      use iamas_atm_boundaries, only : iamas_atm_bdy_checks

      implicit none

      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: blockList
      type (iAMAS_streamManager_type), pointer :: streamManager
      integer, intent(out) :: ierr

      integer :: local_ierr

      ierr = 0

      !
      ! Physics specific checks found in physics/iamas_atmphys_control.F
      !
      call physics_compatibility_check(dminfo, blockList, streamManager, local_ierr)
      ierr = ierr + local_ierr

      !
      ! Checks for limited-area simulations
      !
      call iamas_atm_bdy_checks(dminfo, blockList, streamManager, local_ierr)
      ierr = ierr + local_ierr

   end subroutine iamas_atm_run_compatibility


end module atm_core
