! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module iamas_atm_threading


    contains


    !-----------------------------------------------------------------------
    !  routine iamas_atm_threading_init
    !
    !> \brief  Pre-computes thread loop bounds for cell, edge, and vertex elements
    !> \author Michael Duda
    !> \date   6 July 2015
    !> \details
    !>  This routine is responsible for computing thread loop bounds for cell,
    !>  edge, and vertex elements in each block of the input blocklist argument.
    !>  Starting and ending loop bounds are computed for these three element
    !>  types for all elements (e.g., nCells) as well as owned elements (e.g.,
    !>  nCellsSolve).
    !>
    !>  When iAMAS is compiled without OpenMP support, this routine computes loop
    !>  bounds as though there is just a single thread; otherwise, it is assumed
    !>  that all threads (given by OMP_get_num_threads()) will be used to
    !>  decompose each of the element ranges.
    !>
    !>  At present, a return value of 0 is always returned in the optional
    !>  output argument, ierr.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_threading_init(blocklist, ierr) 

        use iamas_derived_types, only : block_type
        use iamas_pool_routines, only : iamas_pool_get_dimension, iamas_pool_add_dimension
#ifdef iAMAS_OPENMP
        use omp_lib
#endif

        implicit none

        type (block_type), pointer :: blocklist
        integer, intent(out), optional :: ierr

        type (block_type), pointer :: block
        integer :: threadid
        integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve, nVertices, nVerticesSolve

        integer :: nThreads
        integer, dimension(:), pointer :: cellThreadStart, cellThreadEnd
        integer, dimension(:), pointer :: cellSolveThreadStart, cellSolveThreadEnd
        integer, dimension(:), pointer :: edgeThreadStart, edgeThreadEnd
        integer, dimension(:), pointer :: edgeSolveThreadStart, edgeSolveThreadEnd
        integer, dimension(:), pointer :: vertexThreadStart, vertexThreadEnd
        integer, dimension(:), pointer :: vertexSolveThreadStart, vertexSolveThreadEnd

        
        if (present(ierr)) ierr = 0

        block => blocklist
        do while (associated(block))
#ifdef iAMAS_OPENMP
!$OMP PARALLEL
!$OMP MASTER
            nThreads = OMP_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
#else 
            nThreads = 1
#endif 

            allocate(cellThreadStart(nThreads))
            allocate(cellThreadEnd(nThreads))
            allocate(cellSolveThreadStart(nThreads))
            allocate(cellSolveThreadEnd(nThreads))
            allocate(edgeThreadStart(nThreads))
            allocate(edgeThreadEnd(nThreads))
            allocate(edgeSolveThreadStart(nThreads))
            allocate(edgeSolveThreadEnd(nThreads))
            allocate(vertexThreadStart(nThreads))
            allocate(vertexThreadEnd(nThreads))
            allocate(vertexSolveThreadStart(nThreads))
            allocate(vertexSolveThreadEnd(nThreads))

            call iamas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call iamas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call iamas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call iamas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
            call iamas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
            call iamas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

#ifdef iAMAS_OPENMP
!$OMP PARALLEL PRIVATE(threadid)
            threadid = OMP_get_thread_num()

            cellThreadStart(threadid+1) = (threadid * nCells / nThreads) + 1
            cellThreadEnd(threadid+1)   = ((threadid+1) * nCells / nThreads)
            cellSolveThreadStart(threadid+1) = (threadid * nCellsSolve / nThreads) + 1
            cellSolveThreadEnd(threadid+1)   = ((threadid+1) * nCellsSolve / nThreads)
            edgeThreadStart(threadid+1) = (threadid * nEdges / nThreads) + 1
            edgeThreadEnd(threadid+1)   = ((threadid+1) * nEdges / nThreads)
            edgeSolveThreadStart(threadid+1) = (threadid * nEdgesSolve / nThreads) + 1
            edgeSolveThreadEnd(threadid+1)   = ((threadid+1) * nEdgesSolve / nThreads)
            vertexThreadStart(threadid+1) = (threadid * nVertices / nThreads) + 1
            vertexThreadEnd(threadid+1)   = ((threadid+1) * nVertices / nThreads)
            vertexSolveThreadStart(threadid+1) = (threadid * nVerticesSolve / nThreads) + 1
            vertexSolveThreadEnd(threadid+1)   = ((threadid+1) * nVerticesSolve / nThreads)
!$OMP END PARALLEL
#else 
            cellThreadStart(1) = 1
            cellThreadEnd(1)   = nCells
            cellSolveThreadStart(1) = 1
            cellSolveThreadEnd(1)   = nCellsSolve
            edgeThreadStart(1) = 1
            edgeThreadEnd(1)   = nEdges
            edgeSolveThreadStart(1) = 1
            edgeSolveThreadEnd(1)   = nEdgesSolve
            vertexThreadStart(1) = 1
            vertexThreadEnd(1)   = nVertices
            vertexSolveThreadStart(1) = 1
            vertexSolveThreadEnd(1)   = nVerticesSolve
#endif 

            call iamas_pool_add_dimension(block % dimensions, 'nThreads', nThreads)

            call iamas_pool_add_dimension(block % dimensions, 'cellThreadStart', cellThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'cellThreadEnd', cellThreadEnd)
            call iamas_pool_add_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

            call iamas_pool_add_dimension(block % dimensions, 'edgeThreadStart', edgeThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'edgeThreadEnd', edgeThreadEnd)
            call iamas_pool_add_dimension(block % dimensions, 'edgeSolveThreadStart', edgeSolveThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'edgeSolveThreadEnd', edgeSolveThreadEnd)

            call iamas_pool_add_dimension(block % dimensions, 'vertexThreadStart', vertexThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'vertexThreadEnd', vertexThreadEnd)
            call iamas_pool_add_dimension(block % dimensions, 'vertexSolveThreadStart', vertexSolveThreadStart)
            call iamas_pool_add_dimension(block % dimensions, 'vertexSolveThreadEnd', vertexSolveThreadEnd)

            !
            ! Because pools make an internal copy of dimensions, we can now
            ! delete our copies of the thread bounds arrays
            !
            deallocate(cellThreadStart)
            deallocate(cellThreadEnd)
            deallocate(cellSolveThreadStart)
            deallocate(cellSolveThreadEnd)
            deallocate(edgeThreadStart)
            deallocate(edgeThreadEnd)
            deallocate(edgeSolveThreadStart)
            deallocate(edgeSolveThreadEnd)
            deallocate(vertexThreadStart)
            deallocate(vertexThreadEnd)
            deallocate(vertexSolveThreadStart)
            deallocate(vertexSolveThreadEnd)

            block => block % next
        end do

    end subroutine iamas_atm_threading_init


    !-----------------------------------------------------------------------
    !  routine iamas_atm_threading_finalize
    !
    !> \brief  Deallocates memory associated with threading in iAMAS
    !> \author Michael Duda
    !> \date   6 July 2015
    !> \details
    !>  This routine deallocates any memory that was allocated in the call to
    !>  iamas_atm_threading_init() for each block in the input block list.
    !>
    !>  At present, a return value of 0 is always returned in the optional
    !>  output argument, ierr.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_threading_finalize(blocklist, ierr) 

        use iamas_derived_types, only : block_type
        use iamas_pool_routines, only : iamas_pool_remove_dimension

        implicit none

        type (block_type), pointer :: blocklist
        integer, intent(out), optional :: ierr

        type (block_type), pointer :: block


        if (present(ierr)) ierr = 0

        block => blocklist
        do while (associated(block))

            call iamas_pool_remove_dimension(block % dimensions, 'nThreads')

            call iamas_pool_remove_dimension(block % dimensions, 'cellThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'cellThreadEnd')
            call iamas_pool_remove_dimension(block % dimensions, 'cellSolveThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'cellSolveThreadEnd')

            call iamas_pool_remove_dimension(block % dimensions, 'edgeThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'edgeThreadEnd')
            call iamas_pool_remove_dimension(block % dimensions, 'edgeSolveThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'edgeSolveThreadEnd')

            call iamas_pool_remove_dimension(block % dimensions, 'vertexThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'vertexThreadEnd')
            call iamas_pool_remove_dimension(block % dimensions, 'vertexSolveThreadStart')
            call iamas_pool_remove_dimension(block % dimensions, 'vertexSolveThreadEnd')

            block => block % next
        end do

    end subroutine iamas_atm_threading_finalize
 
end module iamas_atm_threading
