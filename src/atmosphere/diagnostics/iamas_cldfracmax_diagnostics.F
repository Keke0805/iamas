! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module cldfracmax_diagnostics

    use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
    use iamas_kind_types, only : RKIND
    use iamas_dmpar
    use iamas_pool_routines
    use iamas_constants
    use iamas_log, only : iamas_log_write

    type (iAMAS_pool_type), pointer :: mesh
    type (iAMAS_pool_type), pointer :: diag_physics
    type (iAMAS_pool_type), pointer :: diag

    type (iAMAS_clock_type), pointer :: clock  

    public :: cldfracmax_diagnostics_setup, &
              cldfracmax_diagnostics_compute

    private

    logical :: need_cldfracmax


    contains


    !-----------------------------------------------------------------------
    !  routine cldfracmax_dignostics_setup
    !
    !> \brief Set up the cldfracmax diagnostics module
    !> \author Ziyu Zhang
    !> \date  20 June 2025
    !> \details
    !>  Compute maximun cldfrac diagnostic before model output is written.
    !>  To simplify cldfrac post-handling
    !
    !-----------------------------------------------------------------------
    subroutine cldfracmax_diagnostics_setup(all_pools, simulation_clock)

        use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
        use iamas_pool_routines, only : iamas_pool_get_subpool

        implicit none

        type (iAMAS_pool_type), pointer :: all_pools
        type (iAMAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call iamas_pool_get_subpool(all_pools, 'mesh', mesh)
        call iamas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call iamas_pool_get_subpool(all_pools, 'diag', diag)

    end subroutine cldfracmax_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine cldfracmax_diagnostics_compute
    !
    !> \brief Compute maximun cldfrac diagnostic before model output is written
    !> \author Ziyu Zhang
    !> \date   20 June 2025
    !> \details
    !>  Compute maximun cldfrac diagnostic before model output is written
    !
    !-----------------------------------------------------------------------
    subroutine cldfracmax_diagnostics_compute()

        use iamas_atm_diagnostics_utils, only : iAMAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_cldfracmax = iAMAS_field_will_be_written('cldfrac_max')
        need_any_diags = need_any_diags .or. need_cldfracmax

        if (need_any_diags) then
            call max_diagnostics(mesh, diag_physics, diag)
        end if
   
    end subroutine cldfracmax_diagnostics_compute


!==================================================================================================
    subroutine max_diagnostics(mesh, diag_physics, diag)
!==================================================================================================

       !input arguments:
        type (iamas_pool_type), intent(in)  :: mesh
        type (iamas_pool_type), intent(in) :: diag_physics

       !inout arguments:
        type (iamas_pool_type), intent(inout) :: diag

       !local variables:
        integer :: iCell
        integer, pointer :: nCells
        
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac
       
        real (kind=RKIND), dimension(:), pointer :: cldfrac_max

       !-------------------------------------

        call iamas_pool_get_dimension(mesh, 'nCells', nCells)

        call iamas_pool_get_array(diag_physics, 'cldfrac', cldfrac)

        call iamas_pool_get_array(diag, 'cldfrac_max', cldfrac_max)

        if (need_cldfracmax) then
            do iCell = 1, nCells
                cldfrac_max(iCell) = maxval(cldfrac(:,iCell))
            enddo
        end if

    end subroutine max_diagnostics

end module cldfracmax_diagnostics
