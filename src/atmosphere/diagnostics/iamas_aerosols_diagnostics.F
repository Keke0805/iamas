! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module aerosols_diagnostics

    use iamas_dmpar
    use iamas_kind_types
    use iamas_derived_types
    use iamas_pool_routines
    use iamas_constants
    use iamas_log, only : iamas_log_write

    type (iAMAS_pool_type), pointer :: mesh
    type (iAMAS_pool_type), pointer :: state
    type (iAMAS_pool_type), pointer :: diag

    type (iAMAS_clock_type), pointer :: clock

    public :: aerosols_diagnostics_setup, &
              aerosols_diagnostics_compute

    private

    logical :: need_dust_column, need_seas_column
    logical :: need_any_column

    contains


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine aerosols_diagnostics_setup(all_pools, simulation_clock)

        use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
        use iamas_pool_routines, only : iamas_pool_get_subpool

        implicit none

        type (iAMAS_pool_type), pointer :: all_pools
        type (iAMAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call iamas_pool_get_subpool(all_pools, 'mesh', mesh)
        call iamas_pool_get_subpool(all_pools, 'state', state)
        call iamas_pool_get_subpool(all_pools, 'diag', diag)
   
    end subroutine aerosols_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in iamas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine aerosols_diagnostics_compute()

        use iamas_atm_diagnostics_utils, only : iAMAS_field_will_be_written
        use module_data_binaer_asect, only : lptr_dust_aer, lptr_cl_aer, lptr_na_aer, nsize_aer, ai_phase

        implicit none

        real (kind=RKIND), dimension(:), pointer :: dust_column, seas_column
        real (kind=RKIND), dimension(:,:), pointer :: zgrid, rho
        real (kind=RKIND), dimension(:,:,:), pointer :: chem
        real (kind=RKIND), dimension(:,:), allocatable :: rhoZ, field
        integer, pointer :: nCells, nVertLevels
        integer :: time_lev              ! which time level to use from state
        integer :: l_dust, l_cl, l_na, iphase, itype, n


        need_any_column = .False.

        need_dust_column = iAMAS_field_will_be_written('dust_column')
        need_seas_column = iAMAS_field_will_be_written('seas_column')

        need_any_column = need_dust_column .or. need_seas_column
        
        itype = 1
        iphase = ai_phase

        if (need_any_column) then
            call iamas_pool_get_dimension(mesh, 'nCells', nCells)
            call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
            call iamas_pool_get_array(mesh, 'zgrid', zgrid)
            call iamas_pool_get_array(diag, 'rho', rho)
            time_lev = 1
            call iamas_pool_get_array(state, 'chem', chem, time_lev)
            if (.not. allocated(rhoZ)) allocate(rhoZ(nVertLevels, nCells))
            if (.not. allocated(field)) allocate(field(nVertLevels, nCells))
            rhoZ(1:nVertLevels,1:nCells) = rho(1:nVertLevels,1:nCells) * &
                                           (zgrid(2:nVertLevels+1, 1:nCells) - zgrid(1:nVertLevels, 1:nCells))
        end if

        if (need_dust_column) then
            field(1:nVertLevels,1:nCells) = 0.0_RKind
            call iamas_pool_get_array(diag, 'dust_column', dust_column)
            do n = 1, nsize_aer(itype)
               l_dust = lptr_dust_aer(n,itype,iphase)
               field(1:nVertLevels,1:nCells) = field(1:nVertLevels,1:nCells) + chem(l_dust,1:nVertLevels,1:nCells)
            enddo
            dust_column(1:nCells) = sum(rhoZ(1:nVertLevels, 1:nCells) * field(1:nVertLevels, 1:nCells), dim=1)
        end if

        if (need_seas_column) then
            field(1:nVertLevels,1:nCells) = 0.0_RKind
            call iamas_pool_get_array(diag, 'seas_column', seas_column)
            do n = 1, nsize_aer(itype)
               l_cl = lptr_cl_aer(n,itype,iphase)
               l_na = lptr_na_aer(n,itype,iphase)
               if (l_cl > 1) field(1:nVertLevels,1:nCells) = field(1:nVertLevels,1:nCells) &
                   + chem(l_cl,1:nVertLevels,1:nCells)
               if (l_na > 1) field(1:nVertLevels,1:nCells) = field(1:nVertLevels,1:nCells) &
                   + chem(l_na,1:nVertLevels,1:nCells)
            enddo
            seas_column(1:nCells) = sum(rhoZ(1:nVertLevels, 1:nCells) * field(1:nVertLevels, 1:nCells), dim=1)
        end if

    end subroutine aerosols_diagnostics_compute
 
end module aerosols_diagnostics
