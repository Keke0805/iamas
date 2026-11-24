! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!***********************************************************************
!
!  iamas_timer
!
!> \brief   iAMAS Timer module
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This module provides developers with internal timer routines. 
!> Additionally it provides standard interfaces to additional timer libarries,
!> such as tau or gptl.
!> Timers are stored as trees of lists of timers.
!
!-----------------------------------------------------------------------
      module iamas_timer

        use iamas_kind_types
        use iamas_derived_types
        use iamas_dmpar
        use iamas_threading
        use iamas_log

#ifdef iAMAS_PERF_MOD_TIMERS
        use perf_mod
#endif

#ifdef iAMAS_GPTL_TIMERS
        use gptl
#endif

        implicit none
        save

#ifdef _PAPI
        include 'f90papi.h'
#endif

        type (iamas_timer_root), pointer :: timer_root => null()

        public :: iamas_timer_set_context, &
                  iamas_timer_start, &
                  iamas_timer_stop, &
                  iamas_timer_write, &
                  iamas_timer_init, &
                  iamas_timer_finalize

        contains

!***********************************************************************
!
!  routine iamas_timer_set_context
!
!> \brief   iAMAS Timer set context routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine sets the timer_root for the timer infrastructure to point to a specific domain's timer_root.
!> This allows multiple cores to have iAMAS timers, and keep their timing information separate.
!
!-----------------------------------------------------------------------
        subroutine iamas_timer_set_context(domain)!{{{

           type (domain_type), intent(in) :: domain

           timer_root => domain % timer_root

        end subroutine iamas_timer_set_context!}}}

!***********************************************************************
!
!  routine iamas_timer_start
!
!> \brief   iAMAS Timer start routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine starts a timer named 'timer_name' as a child of whatever the
!> most recently started timer that is still running is.
!
!-----------------------------------------------------------------------
        subroutine iamas_timer_start(timer_name, clear_timer_in)!{{{

          character (len=*), intent (in) :: timer_name !< Input: name of timer, stored as name of timer
          logical, optional, intent(in) :: clear_timer_in !< Input: flag to clear timer

          logical :: setup_timer, timer_found, clear_timer, check_flag

          type (iamas_timer_node), pointer :: current_timer, prev_timer

          character (len=len_trim(timer_name)) :: trimmed_name

          integer :: usecs, nlen, iErr, threadNum, numThreads

          numThreads = iamas_threading_get_max_threads()
          threadNum = iamas_threading_get_thread_num()

          trimmed_name = trim(timer_name)
          nlen = len(trimmed_name)

          clear_timer = .false.

          if ( present(clear_timer_in) ) then
             clear_timer = clear_timer_in
          end if

#ifdef iAMAS_TAU_TIMERS
          call tau_start(trimmed_name)
#endif

#ifdef iAMAS_PERF_MOD_TIMERS
          call t_startf(trimmed_name)
#endif

#ifdef iAMAS_GPTL_TIMERS
          iErr = gptlstart(trimmed_name)
#endif

#ifdef iAMAS_NATIVE_TIMERS
             if ( threadNum == 0 ) then
                setup_timer = .false.
   
                ! If the root_timer in timer_root has not been set, create it and use it for the current timer
                if ( .not. associated(timer_root % root_timer) ) then
                   allocate(timer_root % root_timer)
                   current_timer => timer_root % root_timer
                   setup_timer = .true.
                ! If root_timer is set, current_timer should be as well. So look for this timer under current_timer, or 
                ! create it as a child.
                else
                   ! If current_timer is not set, set it to root_timer and search...
                   if ( .not. associated(timer_root % current_timer) ) then
                      timer_root % current_timer => timer_root % root_timer
                   end if
   
                   ! If current_timer doesn't have any children, make a new child that will be this timer
                   if ( .not. associated(timer_root % current_timer % child) ) then
                      allocate(timer_root % current_timer % child)
                      current_timer => timer_root % current_timer % child
                      current_timer % parent => timer_root % current_timer
                      setup_timer = .true.
                   else ! Search through current_timer's children for this timer, or create it as a new child
                      timer_found = .false.
                      current_timer => timer_root % current_timer % child
                      do while (associated(current_timer) .and. .not. timer_found)
                         if ( current_timer % nlen == nlen) then
                            if ( current_timer % timer_name(1:current_timer % nlen) == trimmed_name(1:nlen) ) then
                               timer_found = .true.
                            end if
                         end if
   
                         if ( .not. timer_found ) then
                            prev_timer => current_timer
                            current_timer => current_timer % next
                         end if
                      end do
   
                      if (.not. timer_found) then
                         allocate(prev_timer % next)
                         current_timer => prev_timer % next
                         current_timer % parent => timer_root % current_timer
                         setup_timer = .true.
                      end if
                   end if
                end if
   
                ! Setup timer if needed
                if ( setup_timer ) then
                   current_timer % timer_name = trimmed_name
                   current_timer % nlen = nlen
                   current_timer % printed = .false.
                   allocate(current_timer % running(numThreads))
                   allocate(current_timer % start_time(numThreads))
                   allocate(current_timer % end_time(numThreads))
                   allocate(current_timer % total_time(numThreads))
                   allocate(current_timer % max_time(numThreads))
                   allocate(current_timer % min_time(numThreads))
                   allocate(current_timer % avg_time(numThreads))
                   current_timer % running(:) = .false.
                   current_timer % start_time(:) = 0.0_R8KIND
                   current_timer % end_time(:) = 0.0_R8KIND
                   current_timer % total_time(:) = 0.0_R8KIND
                   current_timer % max_time(:) = -huge(RKIND)
                   current_timer % min_time(:) = huge(RKIND)
                   current_timer % avg_time(:) = 0.0_RKIND
                   current_timer % calls = 0
                end if

                ! Set current timer to be current timer (so timers started after will be children of this timer instance)
                timer_root % current_timer => current_timer
             end if
             call iamas_threading_barrier()

             current_timer => timer_root % current_timer

             if ( clear_timer) then
                if ( threadNum == 0 ) then
                   current_timer % calls = 0
                end if
                current_timer % running(threadNum + 1) = .false.
                current_timer % start_time(threadNum + 1) = 0.0_R8KIND
                current_timer % end_time(threadNum + 1) = 0.0_R8KIND
                current_timer % total_time(threadNum + 1) = 0.0_R8KIND
                current_timer % max_time(threadNum + 1) = -huge(RKIND)
                current_timer % min_time(threadNum + 1) = huge(RKIND)
                current_timer % avg_time(threadNum + 1) = 0.0_RKIND
             end if

             ! Set start time
             ! TODO after making sure the timer structures work as expected...
             if ( threadNum == 0 ) then
                current_timer % calls = current_timer % calls + 1
             end if
             current_timer % running(threadNum + 1) = .true.
#ifdef _PAPI
             call PAPIF_get_real_usec(usecs, check_flag)
             current_timer % start_time(threadNum + 1) = usecs/1.0e6
#else
             call iamas_dmpar_get_time(current_timer % start_time(threadNum + 1))
#endif
#endif

        end subroutine iamas_timer_start!}}}

!***********************************************************************
!
!  routine iamas_timer_stop
!
!> \brief   iAMAS Timer stop routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine stops a timer named 'timer_name'. It should be the most
!> recently started timer, so an error is thrown if it is not.
!
!-----------------------------------------------------------------------
        subroutine iamas_timer_stop(timer_name)!{{{

          character (len=*), intent(in) :: timer_name !< Input: name of timer to stop

          character (len=len_trim(timer_name)) :: trimmed_name !< Trimmed timer name

          real (kind=R8KIND) :: temp_time

          integer :: usecs, nlen, iErr, threadNum

          threadNum = iamas_threading_get_thread_num()

          trimmed_name = trim(timer_name)
          nlen = len(trimmed_name)

#ifdef iAMAS_TAU_TIMERS
          call tau_stop(trimmed_name)
#endif

#ifdef iAMAS_PERF_MOD_TIMERS
          call t_stopf(trimmed_name)
#endif

#ifdef iAMAS_GPTL_TIMERS
          iErr = gptlstop(trimmed_name)
#endif

#ifdef iAMAS_NATIVE_TIMERS
          ! Timer to stop should be timer_root % current_timer, since you should only be allowed to stop the most 
          ! recently started timer
          if ( .not. associated(timer_root % current_timer) ) then
             call iamas_log_write('Trying to stop a timer when no timer has been started.', iAMAS_LOG_CRIT)
          end if

          if ( .not. timer_root % current_timer % nlen == nlen ) then
             if ( .not. timer_root % current_timer % timer_name(1:timer_root % current_timer % nlen) == trimmed_name(1:nlen) ) then
                call iamas_log_write('Trying to stop timer ' // trim(trimmed_name) // ' when ' // &
                           trim(timer_root % current_timer % timer_name) // ' is the most recently started timer.', iAMAS_LOG_CRIT)
             end if
          end if

          ! Set stop time
          ! TODO: Compute timer statistics...
          timer_root % current_timer % running(threadNum + 1) = .false.
#ifdef _PAPI
          call PAPIF_get_real_usec(usecs, check_flag)
          current_timer % end_time(threadNum + 1) = usecs/1.0e6
#else
          call iamas_dmpar_get_time(timer_root % current_timer % end_time(threadNum + 1))
#endif
          temp_time = timer_root % current_timer % end_time(threadNum + 1) &
                    - timer_root % current_timer % start_time(threadNum + 1)
          timer_root % current_timer % total_time(threadNum + 1) = timer_root % current_timer % total_time(threadNum + 1) &
                                                                 + temp_time
          timer_root % current_timer % max_time(threadNum + 1) = max(timer_root % current_timer % max_time(threadNum + 1), &
                                                                 real(temp_time, kind=RKIND))
          timer_root % current_timer % min_time(threadNum + 1) = min(timer_root % current_timer % min_time(threadNum + 1), &
                                                                 real(temp_time, kind=RKIND))

          ! Set current_timer to the parent of current_timer
          call iamas_threading_barrier()
          if ( threadNum == 0 ) then
             timer_root % current_timer => timer_root % current_timer % parent
          end if
          call iamas_threading_barrier()
#endif

        end subroutine iamas_timer_stop!}}}

!***********************************************************************
!
!  routine iamas_timer_write_header
!
!> \brief   iAMAS Timer write header routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine writes header information about what will be printed for each timer_root.
!
!-----------------------------------------------------------------------
       subroutine iamas_timer_write_header()!{{{

#ifdef iAMAS_NATIVE_TIMERS

          character(len=StrKIND) :: msg

          call iamas_log_write('')
          call iamas_log_write('')
          call iamas_log_write(' Timer information: ')
          call iamas_log_write('    Globals are computed across all threads and processors')
          call iamas_log_write('')
          call iamas_log_write(' Columns:')
          call iamas_log_write('    calls(call ranks): total start/stop calls summed over MPI ranks')
          call iamas_log_write('    local(s): accumulated time on rank 0')
          call iamas_log_write('    avg/min/max(s): accumulated time averaged/minimum/maximum over MPI ranks')
          call iamas_log_write('    total%: maximum accumulated section time divided by the root maximum')
          call iamas_log_write('')
          call iamas_log_write('')

          write(msg,'(1x, a100, 1x, a16, 1x, a14, 1x, a14, 1x, a14, 1x, a14, 1x, a10)') &
                    'section', 'calls(call ranks)', 'local(s)', 'avg(s)', 'min(s)', 'max(s)', 'total%'
          call iamas_log_write(msg)
#endif

        end subroutine iamas_timer_write_header!}}}

!***********************************************************************
!
!  routine iamas_timer_write
!
!> \brief   iAMAS Timer write routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine writes all timer output to stdout for the currently set
!> timer_root. It prints timer information in depth first format.
!
!-----------------------------------------------------------------------
       subroutine iamas_timer_write()!{{{
          type (iamas_timer_node), pointer :: current_timer, parent_timer
          character (len=StrKIND) :: indentation, section
          logical :: next_timer_found
          integer :: iErr, levels, i, numThreads, global_calls

          real (kind=RKIND) :: percent_inc, percent_exc, efficiency
          real (kind=RKIND) :: global_time, global_max_total, global_ave_total
          real (kind=RKIND) :: local_total, global_sum_total, global_min_total, global_max_total_accum, root_total
          character(len=StrKIND) :: msg

#ifdef iAMAS_GPTL_TIMERS
#ifdef iAMAS_DEBUG
          iErr = gptlpr(timer_root % dminfo % my_proc_id)
#else
          if ( timer_root % dminfo % my_proc_id == 0 ) then
             iErr = gptlpr(timer_root % dminfo % my_proc_id)
          end if
#endif
#endif

#ifdef iAMAS_NATIVE_TIMERS
          numThreads = iamas_threading_get_max_threads()
          current_timer => timer_root % root_timer
          levels = 1

          ! Initialize indentation with all spaces.
          indentation = ''
          do i = 1, StrKIND
             indentation(i:i) = ' '
          end do

          do while ( associated(current_timer) )
             if ( current_timer % running(1) ) then
                call iamas_log_write('Timer ' // trim(current_timer % timer_name) // ' is still running.', iAMAS_LOG_ERR)
             end if

             if ( .not. current_timer % printed ) then
                ! Aggregate accumulated timer time across OpenMP threads and MPI ranks.
                ! The output intentionally matches the requested hotspot table:
                ! local is rank-0 accumulated time; avg/min/max are across-rank accumulated times.
                local_total = real(sum(current_timer % total_time(:)), kind=RKIND)
                call iamas_dmpar_sum_real(timer_root % dminfo, local_total, global_sum_total)
                call iamas_dmpar_min_real(timer_root % dminfo, local_total, global_min_total)
                call iamas_dmpar_max_real(timer_root % dminfo, local_total, global_max_total_accum)
                call iamas_dmpar_sum_int(timer_root % dminfo, current_timer % calls, global_calls)
                global_ave_total = global_sum_total / real(timer_root % dminfo % nprocs, kind=RKIND)
                if (levels == 1) root_total = global_max_total_accum
                percent_inc = global_max_total_accum / max(1.0e-30_RKIND, root_total)

                ! Keep the original timer fields consistent for the traversal and parent percentages.
                current_timer % total_time(1) = global_max_total_accum
                current_timer % avg_time(1) = global_ave_total / max(1, current_timer % calls)
                if ( associated(current_timer % parent) ) then
                   percent_exc = current_timer % total_time(1) / max(1.0e-30_RKIND, &
                                      real(current_timer % parent % total_time(1), kind=RKIND))
                else
                   percent_exc = 0.0_RKIND
                end if

                section = trim(current_timer % timer_name)
                parent_timer => current_timer % parent
                do while ( associated(parent_timer) )
                   section = trim(parent_timer % timer_name) // ':' // trim(section)
                   parent_timer => parent_timer % parent
                end do
                write(msg,'(1x, a100, 1x, i16, 1x, 5f14.4, 1x, f9.2)') trim(section), global_calls, &
                          local_total, global_ave_total, global_min_total, global_max_total_accum, percent_inc * 100.0_RKIND
                call iamas_log_write(msg)
                current_timer % printed = .true.
             end if


             ! If current_timer has children, move down a level (and increment levels)
             ! If current_timer has no children, print siblings (without incrementing levels)
             ! If current_timer has no children, move back to the parent timer, and decrement levels.
             next_timer_found = .false.
             if ( associated(current_timer % child) ) then
                if ( .not. current_timer % child % printed ) then
                    current_timer => current_timer % child
                    levels = levels + 1
                    next_timer_found = .true.
                end if
             end if

             if ( .not. next_timer_found ) then
                if ( associated(current_timer % next) ) then
                   current_timer => current_timer % next
                else
                   if ( associated(current_timer % parent) ) then
                      current_timer => current_timer % parent
                      levels = levels - 1
                   else
                      nullify(current_timer)
                   end if
                end if
             end if
          end do
#endif

        end subroutine iamas_timer_write!}}}

!***********************************************************************
!
!  routine iamas_timer_init
!
!> \brief   iAMAS Timer init routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine initializes the timers within a given domain.
!> Additionally, it sets the context to the timers within the domain, which can
!> be changed later using the iamas_timer_set_context routine.
!
!-----------------------------------------------------------------------
        subroutine iamas_timer_init(domain)!{{{
          type (domain_type), intent(inout) :: domain !< Input/Output: Domain structure

          integer :: iErr

          allocate(domain % timer_root)
          timer_root => domain % timer_root
          timer_root % dminfo => domain % dminfo

#ifdef iAMAS_GPTL_TIMERS
          iErr = gptlsetoption(gptloverhead, 0)
          iErr = gptlsetoption(gptlpercent, 0)
          iErr = gptlsetoption(gptlsync_mpi, 1)

          iErr = gptlinitialize()
#endif

        end subroutine iamas_timer_init!}}}

!***********************************************************************
!
!  routine iamas_timer_finalize
!
!> \brief   iAMAS Timer finalize routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine destroys all timers for the given domain.
!> Timers are destroyed in depth first format
!
!-----------------------------------------------------------------------
        subroutine iamas_timer_finalize(domain)!{{{
          type (domain_type), intent(inout) :: domain !< Input/Output: Domain structure

#ifdef iAMAS_NATIVE_TIMERS
          ! Destroy all children of domain's root_timer
          call iamas_timer_destroy_children(domain % timer_root % root_timer)

          ! Deallocate all arrays of the root_timer, and the timer itself
          deallocate(domain % timer_root % root_timer % start_time)
          deallocate(domain % timer_root % root_timer % end_time)
          deallocate(domain % timer_root % root_timer % total_time)
          deallocate(domain % timer_root % root_timer % max_time)
          deallocate(domain % timer_root % root_timer % min_time)
          deallocate(domain % timer_root % root_timer % avg_time)
          deallocate(domain % timer_root % root_timer % running)
          deallocate(domain % timer_root % root_timer)
#endif

          deallocate(domain % timer_root)

        end subroutine iamas_timer_finalize!}}}

!***********************************************************************
!
!  recursive routine iamas_timer_destroy_chldren
!
!> \brief   iAMAS Timer Destroy Children Routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine destroys all children timers of parent_timer.
!> It does not destroy parent_timer itself.
!
!-----------------------------------------------------------------------
        recursive subroutine iamas_timer_destroy_children(parent_timer)!{{{
           type (iamas_timer_node), pointer :: parent_timer

           type (iamas_timer_node), pointer :: current_timer

           if (.not. associated(parent_timer % child)) then
              return
           end if

           ! First, destroy all children of any sibling in this generation
           current_timer => parent_timer % child
           do while ( associated(current_timer) )
              if ( associated(current_timer % child) ) then
                 call iamas_timer_destroy_children(current_timer)
              end if
              current_timer => current_timer % next
           end do

           ! This generation should not have any children now, so destroy all siblings
           do while ( associated(parent_timer % child) )
              current_timer => parent_timer % child

              if ( associated(current_timer) ) then
                 ! If there is more than one sibling, delete the first one, and set child to be it's sibling
                 if ( associated(current_timer % next) ) then
                    parent_timer % child => current_timer % next
                 else
                    nullify(parent_timer % child)
                 end if

                 ! Deallocate all arrays, and the timer itself
                 deallocate(current_timer % start_time)
                 deallocate(current_timer % end_time)
                 deallocate(current_timer % total_time)
                 deallocate(current_timer % max_time)
                 deallocate(current_timer % min_time)
                 deallocate(current_timer % avg_time)
                 deallocate(current_timer % running)
                 deallocate(current_timer)
              else
                 nullify(parent_timer % child)
              end if
           end do

        end subroutine iamas_timer_destroy_children!}}}

      end module iamas_timer

! vim: foldmethod=marker et ts=2
