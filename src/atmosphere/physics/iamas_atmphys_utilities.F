! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_utilities

 use iamas_kind_types

 implicit none
 private
 public:: physics_error_fatal, &
          physics_message

 character(len=StrKIND),public:: iamas_err_message

!Simple utility subroutines to broadcast messages or abort physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.


 contains


!=================================================================================================================
 subroutine physics_message(str)
!=================================================================================================================

 use iamas_log, only : iamas_log_write

!input arguments:
 character(len=*),intent(in):: str

!-----------------------------------------------------------------------------------------------------------------

 call iamas_log_write(trim(str))

 end subroutine physics_message

!=================================================================================================================
 subroutine physics_error_fatal(str)
!=================================================================================================================

 use iamas_derived_types, only : iAMAS_LOG_ERR, iAMAS_LOG_CRIT
 use iamas_log, only : iamas_log_write

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

 call iamas_log_write(' ',                                                                          messageType=iAMAS_LOG_ERR)
 call iamas_log_write('------------------------------ FATAL CALLED ------------------------------', messageType=iAMAS_LOG_ERR)
 call iamas_log_write(trim(str),                                                                    messageType=iAMAS_LOG_ERR)
 call iamas_log_write('iAMAS core_physics abort', messageType=iAMAS_LOG_CRIT)
 
 end subroutine physics_error_fatal

!=================================================================================================================
 end module iamas_atmphys_utilities
!=================================================================================================================
