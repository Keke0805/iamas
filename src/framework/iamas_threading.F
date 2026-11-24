! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!-----------------------------------------------------------------------
!  iamas_threading
!
!> \brief iAMAS Threading Support
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This module will provide interfaces to support functions / routines for OpenMP threading.
!
!-----------------------------------------------------------------------

module iamas_threading

   use iamas_kind_types

   implicit none
   private

   public :: iamas_threading_get_num_threads, iamas_threading_set_num_threads, iamas_threading_in_parallel
   public :: iamas_threading_get_thread_num, iamas_threading_barrier, iamas_threading_get_max_threads 
   public :: iamas_threading_get_thread_limit

   contains

!-----------------------------------------------------------------------
!  function iamas_threading_get_num_threads
!
!> \brief iAMAS Threading number of threads function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns the number of threads currently available.
!
!-----------------------------------------------------------------------
   function iamas_threading_get_num_threads() result(numThreads)!{{{
      integer :: numThreads
      integer :: omp_get_num_threads

      numThreads = 1

#ifdef iAMAS_OPENMP
      numThreads = omp_get_num_threads()
#endif

   end function iamas_threading_get_num_threads!}}}


!-----------------------------------------------------------------------
!  routine iamas_threading_set_num_threads
!
!> \brief iAMAS Threading set number of threads routine
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This routine sets the number of threads for the next parallel region.
!
!-----------------------------------------------------------------------
   subroutine iamas_threading_set_num_threads(numThreads)!{{{
      integer, intent(in) :: numThreads

#ifdef iAMAS_OPENMP
      call omp_set_num_threads(numThreads)
#endif

   end subroutine iamas_threading_set_num_threads!}}}

!-----------------------------------------------------------------------
!  function iamas_threading_in_parallel
!
!> \brief iAMAS Threading in parallel function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns a logical where true means it was called within a
!>  parallel region, and false means it was not.
!
!-----------------------------------------------------------------------
   function iamas_threading_in_parallel() result(parallelRegion)!{{{
      logical :: parallelRegion
      logical :: omp_in_parallel

      parallelRegion = .false.

#ifdef iAMAS_OPENMP
      parallelRegion = omp_in_parallel()
#endif

   end function iamas_threading_in_parallel!}}}

!-----------------------------------------------------------------------
!  function iamas_threading_get_thread_num
!
!> \brief iAMAS Threading get thread number function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns current thread's number
!
!-----------------------------------------------------------------------
   function iamas_threading_get_thread_num() result(threadNum)!{{{
      integer :: threadNum
      integer :: omp_get_thread_num

      threadNum = 0

#ifdef iAMAS_OPENMP
      threadNum = omp_get_thread_num()
#endif

   end function iamas_threading_get_thread_num!}}}

!-----------------------------------------------------------------------
!  routine iamas_threading_barrier
!
!> \brief iAMAS Threading barrier routine
!> \author Doug Jacobsen
!> \date   10/15/2015
!> \details
!>  This routine implements an OpenMP barrier to synchronize all threads.
!
!-----------------------------------------------------------------------
   subroutine iamas_threading_barrier()!{{{

#ifdef iAMAS_OPENMP
      !$omp barrier
#endif

   end subroutine iamas_threading_barrier!}}}

!-----------------------------------------------------------------------
!  function iamas_threading_get_max_threads
!
!> \brief iAMAS Threading maximum number of threads function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns maximum number of threads a single MPI process can use.
!
!-----------------------------------------------------------------------
   function iamas_threading_get_max_threads() result(maxThreads)!{{{
      integer :: maxThreads
      integer :: omp_get_max_threads

      maxThreads = 1

#ifdef iAMAS_OPENMP
      maxThreads = omp_get_max_threads()
#endif

   end function iamas_threading_get_max_threads!}}}

!-----------------------------------------------------------------------
!  function iamas_threading_get_thread_limit
!
!> \brief iAMAS Threading thread limit function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns limit on the total number of threads.
!
!-----------------------------------------------------------------------
   function iamas_threading_get_thread_limit() result(threadLimit)!{{{
      integer :: threadLimit
      integer :: omp_get_thread_limit

      threadLimit = 1

#ifdef iAMAS_OPENMP
      threadLimit = omp_get_thread_limit()
#endif

   end function iamas_threading_get_thread_limit!}}}

end module iamas_threading
