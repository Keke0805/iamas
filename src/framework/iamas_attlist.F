! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!***********************************************************************
!
!  iamas_attlist
!
!> \brief   iAMAS Attribute list module
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This module provides type definitions and subroutines for working with attribute lists.
!
!-----------------------------------------------------------------------

module iamas_attlist

   use iamas_kind_types
   use iamas_derived_types

   interface iamas_add_att
      module procedure iamas_add_att_int0d
      module procedure iamas_add_att_int1d
      module procedure iamas_add_att_real0d
      module procedure iamas_add_att_real1d
      module procedure iamas_add_att_text
   end interface iamas_add_att

   interface iamas_get_att
      module procedure iamas_get_att_int0d
      module procedure iamas_get_att_int1d
      module procedure iamas_get_att_real0d
      module procedure iamas_get_att_real1d
      module procedure iamas_get_att_text
   end interface iamas_get_att


contains

!***********************************************************************
!
!  routine iamas_add_att_int0d
!
!> \brief   iAMAS Add 0D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 0D integer attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_add_att_int0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = iAMAS_ATT_INT
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueInt = attValue

   end subroutine iamas_add_att_int0d!}}}

!***********************************************************************
!
!  routine iamas_add_att_int1d
!
!> \brief   iAMAS Add 1D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 1D integer attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_add_att_int1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, dimension(:), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = iAMAS_ATT_INTA
      allocate(cursor % attValueIntA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueIntA(:) = attValue(:)

   end subroutine iamas_add_att_int1d!}}}

!***********************************************************************
!
!  routine iamas_add_att_real0d
!
!> \brief   iAMAS Add 0D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 0D real attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_add_att_real0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = iAMAS_ATT_REAL
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueReal = attValue

   end subroutine iamas_add_att_real0d!}}}

!***********************************************************************
!
!  routine iamas_add_att_real1d
!
!> \brief   iAMAS Add 1D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a 1D real attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_add_att_real1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), dimension(:), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = iAMAS_ATT_REALA
      allocate(cursor % attValueRealA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueRealA(:) = attValue(:)

   end subroutine iamas_add_att_real1d!}}}

!***********************************************************************
!
!  routine iamas_add_att_text
!
!> \brief   iAMAS Add text attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine adds a text attribute the attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_add_att_text(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      character (len=*), intent(in) :: attValue !< Input: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = iAMAS_ATT_TEXT
      write(cursor % attName,'(a)') trim(attName)
      write(cursor % attValueText,'(a)') trim(attValue)

   end subroutine iamas_add_att_text!}}}

!***********************************************************************
!
!  routine iamas_get_att_int0d
!
!> \brief   iAMAS get 0D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 0D integer attribute.
!
!-----------------------------------------------------------------------
   subroutine iamas_get_att_int0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= iAMAS_ATT_INT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueInt
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_get_att_int0d!}}}

!***********************************************************************
!
!  routine iamas_get_att_int1d
!
!> \brief   iAMAS get 1D integer attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 1D integer attribute.
!
!-----------------------------------------------------------------------
   subroutine iamas_get_att_int1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, dimension(:), pointer :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= iAMAS_ATT_INTA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueIntA)))
               attValue(:) = cursor % attValueIntA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_get_att_int1d!}}}

!***********************************************************************
!
!  routine iamas_get_att_real0d
!
!> \brief   iAMAS get 0D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 0D real attribute.
!
!-----------------------------------------------------------------------
   subroutine iamas_get_att_real0d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= iAMAS_ATT_REAL) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueReal
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_get_att_real0d!}}}

!***********************************************************************
!
!  routine iamas_get_att_real1d
!
!> \brief   iAMAS get 1D real attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a 1D real attribute.
!
!-----------------------------------------------------------------------
   subroutine iamas_get_att_real1d(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      real (kind=RKIND), dimension(:), pointer :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= iAMAS_ATT_REALA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueRealA)))
               attValue(:) = cursor % attValueRealA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_get_att_real1d!}}}

!***********************************************************************
!
!  routine iamas_get_att_text
!
!> \brief   iAMAS get text attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine returns the attribute value of a text attribute.
!
!-----------------------------------------------------------------------
   subroutine iamas_get_att_text(attList, attName, attValue, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      character (len=*), intent(out) :: attValue !< Output: Attribute value
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= iAMAS_ATT_TEXT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               write(attValue,'(a)') trim(cursor % attValueText)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_get_att_text!}}}

!***********************************************************************
!
!  routine iamas_duplicate_attlist
!
!> \brief   iAMAS duplicate attribute list routine
!> \author  Doug Jacobsen
!> \date    04/26/2016
!> \details 
!> This routine creates a copy of an attribute list, and returns it as destAttList.
!> This routine assumes that destAttList is empty, however it will remove every
!> attribute from destAttList before it attempts to copy into it.
!
!-----------------------------------------------------------------------
   subroutine iamas_duplicate_attlist(srcAttList, destAttList, ierr)!{{{

      implicit none

      type ( att_list_type ), pointer :: srcAttList !< Input: Source attribute list
      type ( att_list_type ), pointer :: destAttList !< Output: Destination attribute list
      integer, intent(out), optional :: ierr !< Optional Output: Error code

      type ( att_list_type ), pointer :: srcCursor

      if ( present(ierr) ) ierr = 0

      if ( associated(destAttList) ) then
         call iamas_deallocate_attlist(destAttList)
      end if

      allocate(destAttList)

      destAttList % attType = -1
      destAttList % attName = ''
      nullify(destAttList % next)
      nullify(destAttList % attValueIntA)
      nullify(destAttList % attValueRealA)

      srcCursor => srcAttList
      do while ( associated(srcCursor) )
         if ( srcCursor % attType == iAMAS_ATT_INT ) then
            call iamas_add_att(destAttList, srcCursor % attName, srcCursor % attValueInt)
         else if ( srcCursor % attType == iAMAS_ATT_INTA ) then
            call iamas_add_att(destAttList, srcCursor % attName, srcCursor % attValueIntA)
         else if ( srcCursor % attType == iAMAS_ATT_REAL ) then
            call iamas_add_att(destAttList, srcCursor % attName, srcCursor % attValueReal)
         else if ( srcCursor % attType == iAMAS_ATT_REALA ) then
            call iamas_add_att(destAttList, srcCursor % attName, srcCursor % attValueRealA)
         else if ( srcCursor % attType == iAMAS_ATT_TEXT ) then
            call iamas_add_att(destAttList, srcCursor % attName, srcCursor % attValueText)
         end if

         srcCursor => srcCursor % next
      end do

   end subroutine iamas_duplicate_attlist!}}}

!***********************************************************************
!
!  routine iamas_remove_att
!
!> \brief   iAMAS remove attribute routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine removes an attribute from an attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_remove_att(attList, attName, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      character (len=*), intent(in) :: attName !< Input: Attribute name
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor, cursor_prev

      if (present(ierr)) ierr = 0

      cursor => attList

      ! Item is at the head of the list
      if (trim(attName) == trim(cursor % attName)) then
         attList => cursor % next
         if (cursor % attType == iAMAS_ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == iAMAS_ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         return
      end if

      cursor_prev => cursor
      cursor => cursor % next
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            cursor_prev % next => cursor % next

            if (cursor % attType == iAMAS_ATT_REALA) then
               deallocate(cursor % attValueRealA)
            else if (cursor % attType == iAMAS_ATT_INTA) then
               deallocate(cursor % attValueIntA)
            end if
            deallocate(cursor)
            
            return
         end if 

         cursor_prev => cursor
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine iamas_remove_att!}}}

!***********************************************************************
!
!  routine iamas_deallocate_attlist
!
!> \brief   iAMAS attribute list deallocation routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This routine deallocates an attribute list.
!
!-----------------------------------------------------------------------
   subroutine iamas_deallocate_attlist(attList, ierr)!{{{

      implicit none

      type (att_list_type), pointer :: attList !< Input/Output: Attribute list
      integer, intent(out), optional :: ierr !< Output: Error flag

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         attList => attList % next
         if (cursor % attType == iAMAS_ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == iAMAS_ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         cursor => attList
      end do

   end subroutine iamas_deallocate_attlist!}}}

end module iamas_attlist
