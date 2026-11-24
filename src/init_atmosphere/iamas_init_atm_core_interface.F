! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module init_atm_core_interface

   use iamas_attlist
   use iamas_derived_types
   use iamas_pool_routines
   use iamas_dmpar
   use iamas_io_units
   use iamas_log, only : iamas_log_write

   contains


   !***********************************************************************
   !
   !  routine init_atm_setup_core
   !
   !> \brief   Init atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_core(core)

      use iamas_derived_types, only : core_type
      use init_atm_core, only : init_atm_core_init, init_atm_core_run, init_atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => init_atm_core_init
      core % core_run => init_atm_core_run
      core % core_finalize => init_atm_core_finalize
      core % define_packages => init_atm_define_packages
      core % setup_packages => init_atm_setup_packages
      core % setup_decompositions => init_atm_setup_decompositions
      core % setup_clock => init_atm_setup_clock
      core % setup_log => init_atm_setup_log
      core % get_mesh_stream => init_atm_get_mesh_stream
      core % setup_immutable_streams => init_atm_setup_immutable_streams
      core % setup_derived_dimensions => init_atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => init_atm_setup_decomposed_dimensions
      core % setup_block => init_atm_setup_block
      core % setup_namelist => init_atm_setup_namelists

      core % Conventions = 'iAMAS'
      core % source = 'iAMAS'

#include "core_variables.inc"

   end subroutine init_atm_setup_core


   !***********************************************************************
   !
   !  routine init_atm_setup_domain
   !
   !> \brief   Init atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_domain(domain)

      use iamas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "domain_variables.inc"

   end subroutine init_atm_setup_domain


   !***********************************************************************
   !
   !  function init_atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  init_atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_packages(configs, packages, iocontext) result(ierr)

      use iamas_derived_types, only : iamas_pool_type, iamas_io_context_type
      use iamas_pool_routines, only : iamas_pool_get_config, iamas_pool_get_package

      implicit none

      type (iamas_pool_type), intent(inout) :: configs
      type (iamas_pool_type), intent(inout) :: packages
      type (iamas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      logical, pointer :: initial_conds, sfc_update, lbcs, lbc_chemistry
      logical, pointer :: emis_update, fdda_update !zhxia add this for emission. 0414/2025
      logical, pointer :: gwd_stage_in, gwd_stage_out, vertical_stage_in, vertical_stage_out, met_stage_in, met_stage_out
      logical, pointer :: config_native_gwd_static, config_static_interp, config_vertical_grid, config_met_interp, config_input_chem
      integer, pointer :: config_init_case


      ierr = 0

      call iamas_pool_get_config(configs, 'config_init_case', config_init_case)
      call iamas_pool_get_config(configs, 'config_static_interp', config_static_interp)
      call iamas_pool_get_config(configs, 'config_native_gwd_static', config_native_gwd_static)
      call iamas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call iamas_pool_get_config(configs, 'config_met_interp', config_met_interp)
      call iamas_pool_get_config(configs, 'config_input_chem', config_input_chem)

      nullify(initial_conds)
      call iamas_pool_get_package(packages, 'initial_condsActive', initial_conds)

      nullify(sfc_update)
      call iamas_pool_get_package(packages, 'sfc_updateActive', sfc_update)

      nullify(lbcs)
      call iamas_pool_get_package(packages, 'lbcsActive', lbcs)

      nullify(lbc_chemistry)
      call iamas_pool_get_package(packages, 'lbc_chemistryActive', lbc_chemistry)

      nullify(emis_update)
      call iamas_pool_get_package(packages, 'emis_updateActive', emis_update)

      nullify(fdda_update)
      call iamas_pool_get_package(packages, 'fdda_updateActive', fdda_update)

      nullify(gwd_stage_in)
      call iamas_pool_get_package(packages, 'gwd_stage_inActive', gwd_stage_in)

      nullify(gwd_stage_out)
      call iamas_pool_get_package(packages, 'gwd_stage_outActive', gwd_stage_out)

      nullify(vertical_stage_in)
      call iamas_pool_get_package(packages, 'vertical_stage_inActive', vertical_stage_in)

      nullify(vertical_stage_out)
      call iamas_pool_get_package(packages, 'vertical_stage_outActive', vertical_stage_out)

      nullify(met_stage_in)
      call iamas_pool_get_package(packages, 'met_stage_inActive', met_stage_in)

      nullify(met_stage_out)
      call iamas_pool_get_package(packages, 'met_stage_outActive', met_stage_out)

      if (.not. associated(initial_conds) .or. &
          .not. associated(sfc_update) .or. &
          .not. associated(lbcs) .or. &
          .not. associated(lbc_chemistry) .or. &
          .not. associated(emis_update) .or. &
          .not. associated(fdda_update) .or. &
          .not. associated(gwd_stage_in) .or. &
          .not. associated(gwd_stage_out) .or. &
          .not. associated(vertical_stage_in) .or. &
          .not. associated(vertical_stage_out) .or. &
          .not. associated(met_stage_in) .or. &
          .not. associated(met_stage_out)) then
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         call iamas_log_write('* Error while setting up packages for init_atmosphere core.',                      messageType=iAMAS_LOG_ERR)
         call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (config_init_case == 8) then
         initial_conds = .false.
         sfc_update = .true.
      else
         initial_conds = .true.
         sfc_update = .false.
      end if

      if (config_init_case == 9) then
         lbcs = .true.
      else
         lbcs = .false.
      end if

      if (config_init_case == 9 .and. config_input_chem) then
         lbc_chemistry = .true.
      else
         lbc_chemistry = .false.
      end if

      ! zhxia add this for emission. 0414/2025
      if (config_init_case == 10) then
         emis_update = .true.
      else
         emis_update = .false.
      end if

      ! jwfeng add this for fdda. 1218/2025
      if (config_init_case == 11) then
         fdda_update = .true.
      else
         fdda_update = .false.
      end if



      if (config_init_case == 7) then

         !
         ! The logic here is a little convoluted
         ! For input, we want to read in fields from all earlier stages, except if those earlier stages are being run now
         ! For output, we want to output the fields that were computed in a stage and all of those from earlier stages
         !
         gwd_stage_in = config_native_gwd_static .and. &
                        (.not. config_static_interp)
         gwd_stage_out = config_native_gwd_static
         vertical_stage_in = config_vertical_grid .and. &
                             (.not. config_native_gwd_static) .and. &
                             (.not. config_static_interp)
         vertical_stage_out = config_vertical_grid
         met_stage_in = config_met_interp .and. &
                        (.not. config_native_gwd_static) .and. &
                        (.not. config_static_interp) .and. &
                        (.not. config_vertical_grid)
         met_stage_out = config_met_interp

      else if (config_init_case == 8) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         vertical_stage_in = .true.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.

      !
      ! When interpolating LBC fields, we need all inputs that would be needed for the interpolation
      ! of ICs, so met_stage_in = .true.
      !
      else if (config_init_case == 9) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .true.
         met_stage_out = .true.

         initial_conds = .false.   ! Also, turn off the initial_conds package to avoid writing the IC "output" stream

      ! zhxia add this for emission. 0414/2025
      else if (config_init_case == 10) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.

         initial_conds = .false.

      else if (config_init_case == 11) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         vertical_stage_in = .true.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.
         initial_conds = .false.

      else
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .true.
      end if

   end function init_atm_setup_packages


   !***********************************************************************
   !
   !  function init_atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_clock(core_clock, configs) result(ierr)

      use iamas_timekeeping, only : iamas_set_time, iamas_set_timeInterval, iamas_create_clock
      use iamas_derived_types, only : iAMAS_Clock_type, iAMAS_Time_type, iAMAS_Timeinterval_type, iamas_pool_type
      use iamas_pool_routines, only : iamas_pool_get_config
      use iamas_kind_types, only : StrKIND

      implicit none

      type (iAMAS_Clock_type), intent(inout) :: core_clock
      type (iamas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_start_time, config_stop_time
      integer, pointer :: config_fg_interval

      type (iAMAS_Time_type)         :: start_time, stop_time
      type (iAMAS_TimeInterval_type) :: dt


      ierr = 0

      call iamas_pool_get_config(configs, 'config_start_time', config_start_time)
      call iamas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call iamas_pool_get_config(configs, 'config_fg_interval', config_fg_interval)

      call iamas_set_time(start_time, dateTimeString=trim(config_start_time))
      call iamas_set_time(stop_time, dateTimeString=trim(config_stop_time))
      call iamas_set_timeInterval(dt, S=config_fg_interval)
      call iamas_create_clock(core_clock, start_time, dt, stopTime=stop_time)

   end function init_atm_setup_clock


   !***********************************************************************
   !
   !  function init_atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_log(logInfo, domain) result(iErr)!{{{

      use iamas_derived_types, only : iamas_log_type, domain_type
      use iamas_log, only : iamas_log_init, iamas_log_open

      implicit none

      type (iamas_log_type), intent(inout), pointer :: logInfo  !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain  !< domain object to provide info for setting up log manager
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call iamas_log_init(logInfo, domain, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call iamas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call iamas_log_write('')
#ifdef SINGLE_PRECISION
      call iamas_log_write('Using default single-precision reals')
#else
      call iamas_log_write('Using default double-precision reals')
#endif
      call iamas_log_write('')


   end function init_atm_setup_log!}}}


   !***********************************************************************
   !
   !  function init_atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_get_mesh_stream(configs, stream) result(ierr)

      use iamas_kind_types, only : StrKIND
      use iamas_derived_types, only : iamas_pool_type
      use iamas_pool_routines, only : iamas_pool_get_config

      implicit none

      type (iamas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      ierr = 0

      write(stream,'(a)') 'input'

   end function init_atm_get_mesh_stream


   !***********************************************************************
   !
   !  function init_atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_decompositions(decompList) result(ierr)

      use iamas_derived_types, only : iamas_decomp_list, iamas_decomp_function, iAMAS_DECOMP_NOERR
      use iamas_decomp, only : iamas_decomp_create_decomp_list, iamas_decomp_register_method, &
                              iamas_uniform_decomp 

      implicit none

      type (iamas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (iamas_decomp_function), pointer :: decompFunc

      ierr = 0

      call iamas_decomp_create_decomp_list(decompList)

      decompFunc => iamas_uniform_decomp

      call iamas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == iAMAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function init_atm_setup_decompositions


   !***********************************************************************
   !
   !  function init_atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  init atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_block(block) result(ierr)

      use iamas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call init_atm_generate_structs(block, block % structs, block % dimensions, block % packages)

   end function init_atm_setup_block


#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module init_atm_core_interface
