module iamas_stream_list

#define COMMA ,
#define LIST_DEBUG_WRITE(M) ! call iamas_log_write(M)
#define LIST_WARN_WRITE(M) call iamas_log_write( M , messageType=iAMAS_LOG_WARN)
#define LIST_ERROR_WRITE(M) call iamas_log_write( M , messageType=iAMAS_LOG_ERR)

    use iamas_kind_types, only : StrKIND
    use iamas_log
    use iamas_derived_types


    contains


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_create
    !
    !> \brief Initialize a new iAMAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Instantiates and initializes a stream_list type, to store all active streams.
    !
    !-----------------------------------------------------------------------
    subroutine iAMAS_stream_list_create(list, ierr) !{{{

        implicit none

        type (iAMAS_stream_list_type), pointer :: list
        integer, intent(out), optional :: ierr


        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_create()')

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR

        allocate(list)

        list % nItems = 0
        nullify(list % head)

    end subroutine iAMAS_stream_list_create !}}}


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_destroy
    !
    !> \brief Free all memory associated with an iAMAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Destroys a stream list type, freeing all memory that was created as
    !>  part of the list.
    !
    !-----------------------------------------------------------------------
    subroutine iAMAS_stream_list_destroy(list, ierr) !{{{

        implicit none

        type (iAMAS_stream_list_type), pointer:: list
        integer, intent(out), optional :: ierr

        type (iAMAS_stream_list_type), pointer :: node


        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_destroy()')

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR

        if (.not. associated(list)) return

        if (list % nItems == 0) then
            deallocate(list)
            return
        end if

        node => list % head
        do while (associated(node))
            list % head => list % head % next
            deallocate(node)
            node => list % head
        end do

        deallocate(list)

    end subroutine iAMAS_stream_list_destroy !}}}


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_insert
    !
    !> \brief Add a stream to a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Adds a stream to the list of streams, first makes sure the stream doesn't exist in the stream list.
    !
    !-----------------------------------------------------------------------
    subroutine iAMAS_stream_list_insert(list, stream, ierr) !{{{

        implicit none

        type (iAMAS_stream_list_type), intent(inout) :: list
        type (iAMAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        type (iAMAS_stream_list_type), pointer :: node


        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_insert()')

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR

        nullify(stream % next)

        if (.not. associated(list % head)) then
            list % head => stream
        else
            node => list % head
            do while (associated(node % next))
                if (node % name == stream % name) then
                    if (present(ierr)) ierr = iAMAS_STREAM_LIST_DUPLICATE
                    LIST_ERROR_WRITE('Found duplicate item '//trim(stream % name)//' in list.')
                    return
                end if
                node => node % next
            end do
            node % next => stream
        end if

        list % nItems = list % nItems + 1

    end subroutine iAMAS_stream_list_insert !}}}


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_remove
    !
    !> \brief Remove a stream from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Removes a stream from the list of streams. Exits if the stream doesn't exist.
    !
    !-----------------------------------------------------------------------
    subroutine iAMAS_stream_list_remove(list, streamName, stream, ierr) !{{{

        implicit none

        type (iAMAS_stream_list_type), intent(inout) :: list
        character (len=*), intent(in) :: streamName
        type (iAMAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        type (iAMAS_stream_list_type), pointer :: node, nodePrev


        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_remove()')

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR

        ! Return if no streams exist in stream list
        if ( list % nItems == 0 ) then
            if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOT_FOUND
            LIST_ERROR_WRITE('Item '//trim(streamName)//' not found in list.')
            nullify(stream)
            return
        end if

        ! Check the head of the stream list
        node => list % head
        if (associated(node)) then
            if (node % name == streamName) then
                list % head => node % next
                stream => node
                list % nItems = list % nItems - 1
                return
            end if
        end if

        ! Loop through all streams until we find the one with StreamName
        nodePrev => node
        node => node % next
        do while (associated(node))
            if (node % name == streamName) then
                nodePrev % next => node % next
                stream => node
                list % nItems = list % nItems - 1
                return
            end if

            nodePrev => node
            node => node % next
        end do

        ! If the routine hasn't returned yet, the stream was not found. Return an error.
        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOT_FOUND
        LIST_ERROR_WRITE('Item '//trim(streamName)//' not found in list.')
        nullify(stream)

    end subroutine iAMAS_stream_list_remove !}}}


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_query
    !
    !> \brief Get streams from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  This function allows one to search through a stream list for a set of
    !>  matching streams. In the first call to this function, the stream
    !>  argument should be a null pointer. When called a single time, this
    !>  function sets the stream argument to the first stream in the list that
    !>  matches the streamPattern argument. Upon subsequent calls, searching
    !>  begins at the first stream *after* the stream argument, to allow one to
    !>  find all streams that match streamPattern.  Each time a valid stream is
    !>  found, this function returns .true. and if no more streams are found, it
    !>  will return .false.
    !
    !-----------------------------------------------------------------------
    logical function iAMAS_stream_list_query(list, streamPattern, stream, ierr) result(found) !{{{

        use iso_c_binding, only: c_char, c_int
        use iamas_c_interfacing, only : iamas_f_to_c_string, iamas_c_to_f_string

        implicit none

        type (iAMAS_stream_list_type), intent(in) :: list
        character (len=*), intent(in) :: streamPattern
        type (iAMAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        character(kind=c_char), dimension(StrKIND+1) :: c_pattern
        character(kind=c_char), dimension(StrKIND+1) :: c_test_string
        integer(kind=c_int) :: c_match

        type (iAMAS_stream_list_type), pointer :: node

        interface
           subroutine check_regex_match(pattern, string, imatch) bind(c)
              use iso_c_binding, only: c_char, c_int
              character(kind=c_char), dimension(*), intent(in) :: pattern, string
              integer(kind=c_int), intent(out) :: imatch
           end subroutine check_regex_match
        end interface

        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_query()')

        call iamas_f_to_c_string(streamPattern, c_pattern)

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR
        found = .false.

        ! Return if no streams exist in stream list
        if ( list % nItems == 0 ) then
            LIST_DEBUG_WRITE(' -- No items matching '//trim(streamPattern)//' found in list.')
            nullify(stream)
            return
        end if

        if ( associated(stream) ) then
           if ( associated(stream % next) ) then
              node => stream % next
           else
              nullify(node)
           end if
        else
           node => list % head
        end if

        do while (associated(node))
            call iamas_f_to_c_string(node % name, c_test_string)

            call check_regex_match(c_pattern, c_test_string, c_match)

            if ( c_match == -1 ) then
                call iamas_log_write('Regular expression matching failed.', iAMAS_LOG_ERR)
                call iamas_log_write('       Pattern was: ' // trim(streamPattern), iAMAS_LOG_ERR)
                call iamas_log_write('       Test string was: ' // trim(node % name), iAMAS_LOG_ERR)
            end if

            if ( c_match == 1 ) then
                found = .true.
                stream => node
                return
            end if
            node => node % next
        end do

        LIST_DEBUG_WRITE(' -- No items matching '//trim(streamPattern)//' found in list.')
        nullify(stream)

    end function iAMAS_stream_list_query !}}}


    !-----------------------------------------------------------------------
    !  routine printlist
    !
    !> \brief Prints the contents of a list.
    !> \author Michael Duda
    !> \date   25 August 2014
    !> \details
    !>  Traverses a list, printing the 'name' component from each node.
    !
    !-----------------------------------------------------------------------
    subroutine printlist(list)!{{{

        implicit none

        type (iAMAS_stream_list_type), intent(in) :: list

        type (iAMAS_stream_list_type), pointer :: node
        integer :: i

        i = 1
        call iamas_log_write('----------------------------')
        call iamas_log_write('List contains:')
        node => list % head
        do while (associated(node))
            call iamas_log_write('    $i) '//trim(node % name), intArgs=(/i/))
            i = i + 1
            node => node % next
        end do
        call iamas_log_write('----------------------------')

    end subroutine printlist!}}}


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_list_length
    !
    !> \brief Returns the length of a stream list
    !> \author Michael Duda
    !> \date   25 August 2014
    !> \details
    !>  Returns the number of items stored in a stream list.
    !
    !-----------------------------------------------------------------------
    integer function iAMAS_stream_list_length(list, ierr) result(nItems) !{{{

        implicit none

        type (iAMAS_stream_list_type), intent(in) :: list
        integer, intent(out), optional :: ierr


        LIST_DEBUG_WRITE(' -- Called iAMAS_stream_list_length()')

        if (present(ierr)) ierr = iAMAS_STREAM_LIST_NOERR

        nItems = list % nItems

    end function iAMAS_stream_list_length !}}}


end module iamas_stream_list
