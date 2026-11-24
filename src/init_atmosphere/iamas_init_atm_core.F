! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module init_atm_core

   use iamas_dmpar
   use iamas_abort, only : iamas_dmpar_global_abort

   contains


   function init_atm_core_init(domain, startTimeStamp) result(ierr)
   
      use iamas_derived_types
      use iamas_stream_manager
      use iamas_io_streams, only : iAMAS_STREAM_NEAREST
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block
      type (iamas_pool_type), pointer :: state, mesh 
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: initial_time
      character (len=StrKIND), pointer :: config_start_time
      integer, pointer :: config_init_case
      logical, pointer :: config_input_chem
      real (kind=RKIND), pointer :: sphere_radius
      integer :: ierr


      ierr = 0

      block => domain % blocklist
      do while (associated(block))
         call iamas_pool_get_subpool(block % structs, 'state', state)
         call iamas_pool_get_subpool(block % structs, 'mesh', mesh)
         call iamas_pool_get_array(state, 'xtime', xtime)
         call iamas_pool_get_array(state, 'initial_time', initial_time)
         call iamas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
         call iamas_pool_get_config(block % configs, 'config_start_time', config_start_time)

         startTimeStamp = config_start_time
         xtime = config_start_time
         initial_time = config_start_time
         domain % sphere_radius = a      ! Appears in output files
         sphere_radius = a               ! Used in setting up test cases

         block => block % next
      end do 

      call iAMAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius, streamID='output', ierr=ierr)
      call iAMAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius, streamID='surface', ierr=ierr)

      !
      ! We don't actually expect the time in the (most likely 'static') file to
      !    match the time in the namelist, so just read whatever time we find in
      !    the input file.
      !
      call iAMAS_stream_mgr_read(domain % streamManager, whence=iAMAS_STREAM_NEAREST, ierr=ierr)
      call iAMAS_stream_mgr_reset_alarms(domain % streamManager, direction=iAMAS_STREAM_INPUT, ierr=ierr)

      !jwfeng from AEMOL USTC: 06/01/2020
      call iamas_pool_get_config(domain % blocklist % configs, 'config_init_case' , config_init_case)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_input_chem', config_input_chem)
      if (config_input_chem .and. config_init_case .eq. 7) then
        call iamas_log_write('****************************************************************')
        call iamas_log_write('Initialize chem from input file.'                                )
        call iamas_log_write('****************************************************************')
        call iAMAS_stream_mgr_read(domain % streamManager, streamID='chem_input', ierr=ierr)
        if (ierr /= iAMAS_STREAM_MGR_NOERR) then
           call iamas_dmpar_global_abort('********************************************************************************', &
                                        deferredAbort=.true.)
           call iamas_dmpar_global_abort('Error reading chemistry initial conditions', &
                                        deferredAbort=.true.)
           call iamas_dmpar_global_abort('********************************************************************************')
        end if
        call iAMAS_stream_mgr_reset_alarms(domain % streamManager, streamID='chem_input', direction=iAMAS_STREAM_INPUT, ierr=ierr)
      end if

   end function init_atm_core_init
   
   
   function init_atm_core_run(domain) result(ierr)
   
      use iamas_derived_types
      use iamas_stream_manager
      use iamas_timer
      use init_atm_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: ierr


      ierr = 0
      
      call init_atm_setup_case(domain, domain % streamManager)
   
   !
   ! Note: The following initialization calls have been moved to the iamas_init_atm_case_*() subroutines,
   !       since values computed by these routines are needed to produce initial fields
   !
   !   call atm_initialize_advection_rk(mesh)
   !   call atm_initialize_deformation_weights(mesh)
  
      call iamas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call iamas_stream_mgr_reset_alarms(domain % streamManager, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)
   
   end function init_atm_core_run
   
   
   function init_atm_core_finalize(domain) result(ierr)
   
      use iamas_derived_types
      use iamas_decomp
      use iamas_stream_manager
      use iamas_log, only : iamas_log_write
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer :: ierr


      ierr = 0

      call iamas_decomp_destroy_decomp_list(domain % decompositions)

      call iamas_log_write('')
      call iamas_log_write('********************************************************')
      call iamas_log_write('   Finished running the init_atmosphere core')
      call iamas_log_write('********************************************************')
   
   end function init_atm_core_finalize
   
end module init_atm_core
