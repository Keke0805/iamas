! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module column_diagnostics

    use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
    use iamas_kind_types, only : RKIND
    use iamas_dmpar
    use iamas_pool_routines
    use iamas_constants
    use iamas_log, only : iamas_log_write

    type (iAMAS_pool_type), pointer :: mesh
    type (iAMAS_pool_type), pointer :: state
    type (iAMAS_pool_type), pointer :: diag

    type (iAMAS_clock_type), pointer :: clock  

    public :: column_diagnostics_setup, &
              column_diagnostics_compute

    private

    logical :: need_qc_column, need_qi_column, need_qs_column, need_qg_column


    contains


    !-----------------------------------------------------------------------
    !  routine column_dignostics_setup
    !
    !> \brief Set up the column diagnostics module
    !> \author Ziyu Zhang
    !> \date  20 June 2025
    !> \details
    !>  Compute Column diagnostic before model output is written.
    !>  To simplify column post-handling
    !
    !-----------------------------------------------------------------------
    subroutine column_diagnostics_setup(all_pools, simulation_clock)

        use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
        use iamas_pool_routines, only : iamas_pool_get_subpool

        implicit none

        type (iAMAS_pool_type), pointer :: all_pools
        type (iAMAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call iamas_pool_get_subpool(all_pools, 'mesh', mesh)
        call iamas_pool_get_subpool(all_pools, 'state', state)
        call iamas_pool_get_subpool(all_pools, 'diag', diag)
   
    end subroutine column_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine column_diagnostics_compute
    !
    !> \brief Compute column diagnostic before model output is written
    !> \author Ziyu Zhang
    !> \date  20 June 2025
    !> \details
    !>  Compute column diagnostic before model output is written
    !
    !-----------------------------------------------------------------------
    subroutine column_diagnostics_compute()

        use iamas_atm_diagnostics_utils, only : iAMAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_qc_column = iAMAS_field_will_be_written('qc_column')
        need_any_diags = need_any_diags .or. need_qc_column
        need_qi_column = iAMAS_field_will_be_written('qi_column')
        need_any_diags = need_any_diags .or. need_qi_column
        need_qs_column = iAMAS_field_will_be_written('qs_column')
        need_any_diags = need_any_diags .or. need_qs_column
        need_qg_column = iAMAS_field_will_be_written('qg_column')
        need_any_diags = need_any_diags .or. need_qg_column

        if (need_any_diags) then
            call sum_diagnostics(mesh, state, 1, diag)
        end if
   
    end subroutine column_diagnostics_compute


!==================================================================================================
    subroutine sum_diagnostics(mesh, state, time_lev, diag)
!==================================================================================================

       !input arguments:
        type (iamas_pool_type), intent(in)  :: mesh
        type (iamas_pool_type), intent(in) :: state
        integer, intent(in) :: time_lev              ! which time level to use from state
       
       !inout arguments:
        type (iamas_pool_type), intent(inout) :: diag

       !local variables:
        integer :: iCell, iVert
        integer, pointer :: nCells, nVertLevels
        integer, pointer :: index_qc, index_qi, index_qs, index_qg
        
        real (kind=RKIND), dimension(:,:), pointer :: height
        real (kind=RKIND), dimension(:,:), pointer :: qc, qi, qs, qg
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
       
        real (kind=RKIND), dimension(:), pointer :: qc_column
        real (kind=RKIND), dimension(:), pointer :: qi_column
        real (kind=RKIND), dimension(:), pointer :: qs_column
        real (kind=RKIND), dimension(:), pointer :: qg_column

       !--------- local fields---------------

        real (kind=RKIND), dimension(:,:), allocatable :: height_d

       !-------------------------------------

        call iamas_pool_get_dimension(mesh, 'nCells', nCells)
        call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)


        call iamas_pool_get_dimension(state, 'index_qc', index_qc)
        call iamas_pool_get_dimension(state, 'index_qi', index_qi)
        call iamas_pool_get_dimension(state, 'index_qs', index_qs)
        call iamas_pool_get_dimension(state, 'index_qg', index_qg)

        call iamas_pool_get_array(mesh, 'zgrid', height)

        call iamas_pool_get_array(state, 'scalars', scalars, time_lev)
        qc => scalars(index_qc,:,:)
        qi => scalars(index_qi,:,:)
        qs => scalars(index_qs,:,:)
        qg => scalars(index_qg,:,:)

        call iamas_pool_get_array(diag, 'qc_column', qc_column)
        call iamas_pool_get_array(diag, 'qi_column', qi_column)
        call iamas_pool_get_array(diag, 'qs_column', qs_column)
        call iamas_pool_get_array(diag, 'qg_column', qg_column)

        if(.not. allocated(height_d)) allocate(height_d(nVertLevels,nCells))

        height_d = height(2:nVertLevels+1,:) - height(1:nVertLevels,:)

        if (need_qc_column) then
            qc_column = 0.0
            do iCell = 1, nCells
            do iVert = 1, nVertLevels
                qc_column(iCell) = qc_column(iCell) + (height_d(iVert,iCell)*qc(iVert,iCell))*1000000.
            enddo
            enddo
        end if

        if (need_qi_column) then
            qi_column = 0.0
            do iCell = 1, nCells
            do iVert = 1, nVertLevels
                qi_column(iCell) = qi_column(iCell) + (height_d(iVert,iCell)*qi(iVert,iCell))*1000000.
            enddo
            enddo
        end if

        if (need_qs_column) then
            qs_column = 0.0
            do iCell = 1, nCells
            do iVert = 1, nVertLevels
                qs_column(iCell) = qs_column(iCell) + (height_d(iVert,iCell)*qs(iVert,iCell))*1000000.
            enddo
            enddo
        end if

        if (need_qg_column) then
            qg_column = 0.0
            do iCell = 1, nCells
            do iVert = 1, nVertLevels
                qg_column(iCell) = qg_column(iCell) + (height_d(iVert,iCell)*qg(iVert,iCell))*1000000.
            enddo
            enddo
        end if

        if (allocated(height_d)) deallocate(height_d)

    end subroutine sum_diagnostics

end module column_diagnostics
