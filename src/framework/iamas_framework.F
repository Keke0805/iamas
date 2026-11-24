! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!-----------------------------------------------------------------------
!  iamas_framework
!
!> \brief iAMAS Framework routines
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module contains all routines related to the general iAMAS framework interface.
!
!-----------------------------------------------------------------------
module iamas_framework

   use iamas_dmpar
   use iamas_derived_types
   use iamas_domain_routines
   use iamas_pool_routines
   use iamas_timer
   use iamas_timekeeping
   use iamas_io
   use iamas_io_units
   use iamas_block_decomp


   contains


!-----------------------------------------------------------------------
!  routine iamas_framework_init_phase1
!
!> \brief iAMAS framework initialization phase 1 routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine initializes the first parts of iAMAS framework. It initializes
!>  MPI, the log unit numbers.
!
!-----------------------------------------------------------------------
   subroutine iamas_framework_init_phase1(dminfo, mpi_comm)!{{{

      implicit none

      type (dm_info), pointer :: dminfo
      integer, intent(in), optional :: mpi_comm

      allocate(dminfo)
      call iamas_dmpar_init(dminfo, mpi_comm)

   end subroutine iamas_framework_init_phase1!}}}

!-----------------------------------------------------------------------
!  routine iamas_framework_init_phase2
!
!> \brief iAMAS framework initialization phase 2 routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine finalizes the initialization of the iAMAS framework. It calls initializes 
!>  the time keeper, and the IO infrastructure.
!
!-----------------------------------------------------------------------
   subroutine iamas_framework_init_phase2(domain, io_system, calendar)!{{{
      use mpi

      implicit none

      type (domain_type), pointer :: domain

      type (iosystem_desc_t), optional, pointer :: io_system
      character(len=*), intent(in), optional :: calendar

      character(len=StrKIND), pointer :: config_calendar_type
      integer, pointer :: config_pio_num_iotasks, config_pio_stride
      integer, pointer :: config_pio_local_num_iotasks, config_pio_local_stride
      logical, pointer :: config_io_use_multiple_files, config_io_do_output_indices
      logical, pointer :: config_is_init_atm
      integer :: pio_num_iotasks
      integer :: pio_stride
      
      integer :: local_pio_num_iotasks
      integer :: local_pio_stride

      integer :: mpi_io_rank, mpi_io_size
      integer :: mpi_io_comm, mpi_ierr
     
      integer :: mpi_local_rank, mpi_local_size
      integer :: mpi_local_comm
      integer :: group_size
      integer :: color

      call iamas_timer_init(domain)

#ifdef iAMAS_DEBUG
      call iamas_pool_set_error_level(iAMAS_POOL_WARN)
#endif

      call iamas_pool_get_config(domain % configs, 'config_calendar_type', config_calendar_type)
      call iamas_pool_get_config(domain % configs, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call iamas_pool_get_config(domain % configs, 'config_pio_stride', config_pio_stride)

      call iamas_pool_get_config(domain % configs, 'config_pio_local_num_iotasks', config_pio_local_num_iotasks)
      call iamas_pool_get_config(domain % configs, 'config_pio_local_stride', config_pio_local_stride)
      
      call iamas_pool_get_config(domain % configs, 'config_io_use_multiple_files', config_io_use_multiple_files)
      call iamas_pool_get_config(domain % configs, 'config_io_do_output_indices', config_io_do_output_indices)
      
      call iamas_pool_get_config(domain % configs, 'config_is_init_atm', config_is_init_atm)

      
      domain % ioContext % io_use_mulfiles => config_io_use_multiple_files
      domain % ioContext % do_output_indices => config_io_do_output_indices
      
      domain % ioContext % is_init_atm => config_is_init_atm
      call iamas_log_write( 'config_is_init_atm = $l', logicArgs=(/config_is_init_atm/) )

      if(domain % ioContext % is_init_atm) then
         call iamas_log_write('Init atmosphere ...')
      else
         call iamas_log_write('Atmosphere ...')
         if(domain % ioContext % io_use_mulfiles)then
            call iamas_log_write('I/O using multiple files...')
         else
            call iamas_log_write('I/O using single file...')
         end if
      end if


      if (present(calendar)) then
         call iamas_timekeeping_init(calendar)
      else
         call iamas_timekeeping_init(config_calendar_type)
      end if

      pio_num_iotasks = config_pio_num_iotasks
      pio_stride = config_pio_stride

      group_size = pio_stride
      color = domain%dminfo%my_proc_id/group_size
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, mod(domain%dminfo%my_proc_id, group_size), mpi_local_comm, mpi_ierr)
      call MPI_Comm_dup(mpi_local_comm, domain%dminfo % local_comm, mpi_ierr)
      
      call MPI_Comm_rank(mpi_local_comm, mpi_local_rank, mpi_ierr)
      call MPI_Comm_size(mpi_local_comm, mpi_local_size, mpi_ierr)

      domain%dminfo % local_nprocs = mpi_local_size
      domain%dminfo % local_proc_id = mpi_local_rank

      color = mod(domain%dminfo%my_proc_id, group_size)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, domain%dminfo%my_proc_id / group_size, mpi_io_comm, mpi_ierr)
      call MPI_Comm_dup(mpi_io_comm, domain%dminfo % io_comm, mpi_ierr)
      
      call MPI_Comm_rank(mpi_io_comm, mpi_io_rank, mpi_ierr)
      call MPI_Comm_size(mpi_io_comm, mpi_io_size, mpi_ierr)

      domain%dminfo % io_nprocs = mpi_io_size
      domain%dminfo % io_proc_id = mpi_io_rank

      if (pio_num_iotasks == 0) then
         pio_num_iotasks = domain % dminfo % nprocs
      end if


      domain % ioContext % dminfo => domain % dminfo
   
      if(domain % ioContext % io_use_mulfiles)then
          local_pio_num_iotasks = config_pio_local_num_iotasks
          local_pio_stride = config_pio_local_stride
          
          if (local_pio_num_iotasks == 0) then
             local_pio_num_iotasks = 1
          end if
          call iAMAS_io_init_mul(domain % ioContext, pio_num_iotasks, pio_stride, io_system, local_pio_num_iotasks, local_pio_stride)
      else
          call iAMAS_io_init(domain % ioContext, pio_num_iotasks, pio_stride, io_system)
      end if

   end subroutine iamas_framework_init_phase2!}}}


!-----------------------------------------------------------------------
!  routine iamas_framework_finalize
!
!> \brief iAMAS framework finalization routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine finalizes the iAMAS framework. It calls routines related to finalizing different parts of iAMAS, that are housed within the framework.
!
!-----------------------------------------------------------------------  
   subroutine iamas_framework_finalize(dminfo, domain, io_system)!{{{
  
      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain
      type (iosystem_desc_t), optional, pointer :: io_system

      call iAMAS_io_finalize(domain % ioContext, .false.)

      call iamas_deallocate_domain(domain)

      call iamas_dmpar_finalize(dminfo)

      call iamas_finish_block_proc_list(dminfo)

      call iamas_timekeeping_finalize()

   end subroutine iamas_framework_finalize!}}}

end module iamas_framework
