!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Some chemistry utilities.
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 11/2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!     iAMAS-v5.2:src/core_atmosphere/physics/iamas_atmphys_utilities.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_utilities

 use iamas_kind_types
! use iamas_abort, only : iamas_dmpar_global_abort
 use iamas_log,only : iamas_log_write
 use iamas_derived_types, only : iAMAS_LOG_ERR, iAMAS_LOG_CRIT



 implicit none
 private
 public:: chemistry_error_fatal, &
          chemistry_error_message, &
          chemistry_message, &
          chemistry_debug_message

 character(len=StrKIND),public:: iamas_err_message

!jwfeng from AEMOL.
!Simple utility subroutines to broadcast messages or abort chemistry parameterizations.



 contains


!=================================================================================================================
 subroutine chemistry_message(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!-----------------------------------------------------------------------------------------------------------------

 call iamas_log_write(trim(str))

 end subroutine chemistry_message

!=================================================================================================================
 subroutine chemistry_debug_message(str)
!=================================================================================================================
 use iamas_atmchem_vars, only: detailed_debug_on

!input arguments:
 character(len=*),intent(in):: str

!-----------------------------------------------------------------------------------------------------------------

 if (detailed_debug_on) call iamas_log_write(trim(str))

 end subroutine chemistry_debug_message


!=================================================================================================================
 subroutine chemistry_error_message(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

 call iamas_log_write('   ',                                                                        messageType=iAMAS_LOG_ERR)
 call iamas_log_write('--------------------------iAMAS core_chemistry ERROR-----------------------', messageType=iAMAS_LOG_ERR)
 call iamas_log_write(trim(str),                                                                    messageType=iAMAS_LOG_ERR)

 end subroutine chemistry_error_message


!=================================================================================================================
 subroutine chemistry_error_fatal(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

!call iamas_dmpar_global_abort(' ',                                                                          deferredAbort=.true.)
!call iamas_dmpar_global_abort('------------------------------ FATAL CALLED ------------------------------', deferredAbort=.true.)
!call iamas_dmpar_global_abort(trim(str),                                                                    deferredAbort=.true.)
!call iamas_dmpar_global_abort('ERROR: iAMAS core_chemistry abort')

 call iamas_log_write('   ',                                                                        messageType=iAMAS_LOG_ERR)
 call iamas_log_write('------------------------------ FATAL CALLED ------------------------------', messageType=iAMAS_LOG_ERR)
 call iamas_log_write(trim(str),                                                                    messageType=iAMAS_LOG_ERR)
 call iamas_log_write('iAMAS core_chemistry abort', messageType=iAMAS_LOG_CRIT)
 
 end subroutine chemistry_error_fatal

!=================================================================================================================
 end module iamas_atmchem_utilities
!=================================================================================================================
