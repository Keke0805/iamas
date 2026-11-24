module iamas_c_interfacing


    contains


    !-----------------------------------------------------------------------
    !  routine iamas_sanitize_string
    !
    !> \brief Converts C null characters in a Fortran string to spaces
    !> \author Michael Duda
    !> \date   19 February 2019
    !> \details
    !>  Converts all C null characters in a Fortran string to spaces.
    !>  This may be useful for strings that were provided by C code through other
    !>  Fortran code external to iAMAS.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_sanitize_string(str)
    
        use iso_c_binding, only : c_null_char
    
        implicit none
    
        character(len=*), intent(inout) :: str
    
        integer :: i
    
        do i=1,len(str)
            if (str(i:i) == c_null_char) then
                str(i:i) = ' '
            end if
        end do
    
    end subroutine iamas_sanitize_string


    !-----------------------------------------------------------------------
    !  routine iamas_c_to_f_string
    !
    !> \brief Converts a C null-terminated string to a Fortran string
    !> \author Michael Duda
    !> \date   9 July 2014
    !> \details
    !>  Converts a C null-terminated string to a Fortran string.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_c_to_f_string(cstring, fstring)
    
        
        use iso_c_binding, only : c_char, c_null_char
    
        implicit none
    
        character(kind=c_char), dimension(*), intent(in) :: cstring
        character(len=*), intent(out) :: fstring
    
        integer :: i, j
    
        j = len(fstring)
        do i=1,j
            if (cstring(i) == c_null_char) exit
        end do
        if (i > j) then
            i = j
        else
            i = i - 1
        end if
        fstring(1:i) = transfer(cstring(1:i), fstring)
        fstring = fstring(1:i)
    
    end subroutine iamas_c_to_f_string
    
    
    !-----------------------------------------------------------------------
    !  routine iamas_f_to_c_string
    !
    !> \brief Converts a Fortran string to a C null-terminated string
    !> \author Michael Duda
    !> \date   9 July 2014
    !> \details
    !>  Converts a Fortran string to a C null-terminated string.
    !>  The output argument is an assumed-size array that must be large enough
    !>  to contain the Fortran string, plus a c_null_char character.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_f_to_c_string(fstring, cstring)
    
        use iso_c_binding, only : c_char, c_null_char
        
        implicit none
        
        character(len=*), intent(in) :: fstring
        character(kind=c_char), dimension(*), intent(out) :: cstring
            
        integer :: i
                
        do i=1,len_trim(fstring)
            cstring(i) = fstring(i:i)
        end do
        cstring(i) = c_null_char 
 
    end subroutine iamas_f_to_c_string

end module iamas_c_interfacing
