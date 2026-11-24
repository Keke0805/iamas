! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module iamas_atm_diagnostics_manager

    private

    public :: iamas_atm_diag_setup, &
              iamas_atm_diag_update, &
              iamas_atm_diag_compute, &
              iamas_atm_diag_reset, &
              iamas_atm_diag_cleanup

    contains


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_setup
    !
    !> \brief Initialize the diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Initialize the diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_setup(stream_mgr, configs, structs, clock, dminfo)

        use iamas_atm_diagnostics_utils, only : iamas_atm_diag_utils_init
        use iamas_derived_types, only : iAMAS_streamManager_type, iAMAS_pool_type, iAMAS_clock_type, dm_info
        use diagnostic_template, only : diagnostic_template_setup
        use isobaric_diagnostics, only : isobaric_diagnostics_setup
        use convective_diagnostics, only : convective_diagnostics_setup
        use pv_diagnostics, only : pv_diagnostics_setup
        use soundings, only : soundings_setup
        use aerosols_diagnostics, only : aerosols_diagnostics_setup
        use column_diagnostics, only : column_diagnostics_setup
        use isohypsic_diagnostics, only : isohypsic_diagnostics_setup
        use cldfracmax_diagnostics, only : cldfracmax_diagnostics_setup

        implicit none

        type (iAMAS_streamManager_type), target, intent(inout) :: stream_mgr
        type (iAMAS_pool_type), pointer :: configs
        type (iAMAS_pool_type), pointer :: structs
        type (iAMAS_clock_type), pointer :: clock
        type (dm_info), intent(in) :: dminfo


        !
        ! Prepare the diagnostics utilities module for later use by diagnostics
        !
        call iamas_atm_diag_utils_init(stream_mgr)

        call diagnostic_template_setup(configs, structs, clock)
        call isobaric_diagnostics_setup(structs, clock)
        call convective_diagnostics_setup(structs, clock)
        call pv_diagnostics_setup(configs, structs, clock)
        call soundings_setup(configs, structs, clock, dminfo)
        call aerosols_diagnostics_setup(structs, clock)
        call column_diagnostics_setup(structs, clock)
        call isohypsic_diagnostics_setup(structs, clock)
        call cldfracmax_diagnostics_setup(structs, clock)

    end subroutine iamas_atm_diag_setup


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_update
    !
    !> \brief Handle diagnostics accumulation at the end of each timestep
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  iAMAS_atm_diag_update.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_update()

        use diagnostic_template, only : diagnostic_template_update
        use convective_diagnostics, only : convective_diagnostics_update

        implicit none


        call diagnostic_template_update()
        call convective_diagnostics_update()

    end subroutine iamas_atm_diag_update


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_compute
    !
    !> \brief Compute diagnostics before they are written to output streams
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  iAMAS_atm_diag_compute.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_compute()

        use diagnostic_template, only : diagnostic_template_compute
        use isobaric_diagnostics, only : isobaric_diagnostics_compute
        use convective_diagnostics, only : convective_diagnostics_compute
        use pv_diagnostics, only : pv_diagnostics_compute
        use soundings, only : soundings_compute
        use aerosols_diagnostics, only : aerosols_diagnostics_compute
        use column_diagnostics, only : column_diagnostics_compute
        use isohypsic_diagnostics, only : isohypsic_diagnostics_compute 
        use cldfracmax_diagnostics, only : cldfracmax_diagnostics_compute

        implicit none


        call diagnostic_template_compute()
        call isobaric_diagnostics_compute()
        call convective_diagnostics_compute()
        call pv_diagnostics_compute()
        call soundings_compute()
        call aerosols_diagnostics_compute()
        call column_diagnostics_compute()
        call isohypsic_diagnostics_compute()
        call cldfracmax_diagnostics_compute()

    end subroutine iamas_atm_diag_compute


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_reset
    !
    !> \brief Resets a diagnostic after it has been computed and written
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  iAMAS_atm_diag_reset.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_reset()

        use diagnostic_template, only : diagnostic_template_reset
        use convective_diagnostics, only : convective_diagnostics_reset

        implicit none


        call diagnostic_template_reset()
        call convective_diagnostics_reset()

    end subroutine iamas_atm_diag_reset


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_cleanup
    !
    !> \brief Finalizes diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_cleanup()

        use iamas_atm_diagnostics_utils, only : iamas_atm_diag_utils_finalize
        use diagnostic_template, only : diagnostic_template_cleanup
        use soundings, only : soundings_cleanup

        implicit none


        call diagnostic_template_cleanup()
        call soundings_cleanup()

        !
        ! Take care of any needed cleanup in the diagnostics utility module
        !
        call iamas_atm_diag_utils_finalize()

    end subroutine iamas_atm_diag_cleanup

end module iamas_atm_diagnostics_manager
