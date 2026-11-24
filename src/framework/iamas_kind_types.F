! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!***********************************************************************
!
!  iamas_kind_types
!
!> \brief   iAMAS Kind definition module
!> \author  Michael Duda, Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This module defines the kind types for basic fortran data types within iAMAS.
!
!-----------------------------------------------------------------------

module iamas_kind_types

   integer, parameter :: R4KIND = selected_real_kind(6)
   integer, parameter :: R8KIND = selected_real_kind(12)
#ifdef SINGLE_PRECISION
   integer, parameter :: RKIND  = selected_real_kind(6)
#else
   integer, parameter :: RKIND  = selected_real_kind(12)
#endif

   integer, parameter :: I8KIND = selected_int_kind(18)

   integer, parameter :: StrKIND = 512
   integer, parameter :: ShortStrKIND = 64

   contains

!***********************************************************************
!
!  routine dummy_kinds
!
!> \brief   iAMAS Dummy kind routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This is a dummy routine that doesn't do anything.
!
!-----------------------------------------------------------------------
   subroutine dummy_kinds()

   end subroutine dummy_kinds

end module iamas_kind_types
