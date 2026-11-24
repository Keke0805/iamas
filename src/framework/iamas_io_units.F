! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!***********************************************************************
!
!  iamas_io_units
!
!> \brief   iAMAS Logging module
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This module contains routines for handling log and error files.
!
!-----------------------------------------------------------------------

module iamas_io_units

   use iamas_kind_types

   integer, parameter, private :: maxUnits = 99
   logical, dimension(0:maxUnits), private, save :: unitsInUse

   contains

!***********************************************************************
!
!  routine iamas_new_unit
!
!> \brief   iAMAS New unit routine
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine determines a new unit that is not yet in use, and returns
!> the unit number
!
!-----------------------------------------------------------------------
    subroutine iamas_new_unit(newUnit)!{{{
        integer, intent(inout) :: newUnit

        integer :: i

        logical :: opened

        do i = 1, maxUnits
            if (.not. unitsInUse(i)) then
                inquire(i, opened=opened)
                if (opened) then
                    unitsInUse(i) = .true.
                else
                    newUnit = i
                    unitsInUse(newUnit) = .true.
                    return
                endif
            end if
        end do

    end subroutine iamas_new_unit!}}}

!***********************************************************************
!
!  routine iamas_release_unit
!
!> \brief   iAMAS Release unit routine
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine releases a unit that is in use.
!
!-----------------------------------------------------------------------
    subroutine iamas_release_unit(releasedUnit)!{{{
        integer, intent(in) :: releasedUnit

        unitsInUse(releasedUnit) = .false.

    end subroutine iamas_release_unit!}}}

end module iamas_io_units
