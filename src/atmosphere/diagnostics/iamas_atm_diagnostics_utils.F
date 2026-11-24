! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module iamas_atm_diagnostics_utils

    use iamas_derived_types, only : iAMAS_streamManager_type 

    private

    public :: iamas_atm_diag_utils_init, &
              iamas_atm_diag_utils_finalize, &
              iamas_field_will_be_written, &
              iamas_stream_inclusion_count


    type (iAMAS_streamManager_type), pointer :: streamManager

    contains


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_utils_init
    !
    !> \brief Prepares diagnostics utilities for use by diagnostics modules
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This routine prepares the diagnostics utilities for use by diagnostics
    !>  modules.
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_utils_init(stream_mgr)

        use iamas_derived_types, only : iAMAS_streamManager_type

        implicit none

        type (iAMAS_streamManager_type), target, intent(inout) :: stream_mgr

        streamManager => stream_mgr

    end subroutine iamas_atm_diag_utils_init


    !-----------------------------------------------------------------------
    !  routine iAMAS_atm_diag_utils_finalize
    !
    !> \brief Performs cleanup after diagnostics utilities will no longer be used
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This routine performs cleanup after diagnostics utilities will no 
    !>  longer be used
    !
    !-----------------------------------------------------------------------
    subroutine iamas_atm_diag_utils_finalize()

        implicit none

        ! Nothing to do here at present...

    end subroutine iamas_atm_diag_utils_finalize


    !-----------------------------------------------------------------------
    !  routine iAMAS_field_will_be_written
    !
    !> \brief Decide if a field will be written in next call to iamas_stream_mgr_write
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This function queries the stream manager to see whether there are any
    !>  streams that contain the field 'fieldName' and whose output alarms
    !>  are also ringing. If so, the function returns .true..
    !>  The assumption is that, between the call to this function and the next
    !>  call to write all streams with iamas_stream_mgr_write(), the stream
    !>  (or streams) containing the named field will not have their alarms
    !>  externally reset.
    !
    !-----------------------------------------------------------------------
    logical function iamas_field_will_be_written(fieldName)

        use iamas_kind_types, only : StrKIND
        use iamas_derived_types, only : iAMAS_STREAM_OUTPUT, iAMAS_STREAM_INPUT_OUTPUT
        use iamas_stream_manager, only : iamas_stream_mgr_begin_iteration, iamas_stream_mgr_get_next_stream, &
                                        iAMAS_stream_mgr_ringing_alarms, iamas_stream_mgr_get_next_field

        implicit none

        character(len=*), intent(in) :: fieldName

        character (len=StrKIND) :: streamNameItr
        character (len=StrKIND) :: fieldNameItr
        integer :: streamDirection
        logical :: streamActive
        logical :: fieldActive
        integer :: ierr

        iamas_field_will_be_written = .false.

        call iamas_stream_mgr_begin_iteration(streamManager)
        do while (iamas_stream_mgr_get_next_stream(streamManager, streamID = streamNameItr, &
                                                  directionProperty = streamDirection, activeProperty = streamActive))

            if (streamActive .and. ( streamDirection == iAMAS_STREAM_OUTPUT .or. streamDirection == iAMAS_STREAM_INPUT_OUTPUT )) then

                if (iAMAS_stream_mgr_ringing_alarms(streamManager, streamID=streamNameItr, &
                                                   direction=iAMAS_STREAM_OUTPUT, ierr=ierr)) then

                    call iamas_stream_mgr_begin_iteration(streamManager, streamID=streamNameItr)
                    do while (iamas_stream_mgr_get_next_field(streamManager, streamNameItr, fieldNameItr, isActive=fieldActive))

                        if (fieldActive .and. (fieldNameItr == fieldName)) then
                            iamas_field_will_be_written = .true.
                            return
                        end if

                    end do
                end if

            end if

        end do

    end function iamas_field_will_be_written


    !-----------------------------------------------------------------------
    !  routine iAMAS_stream_inclusion_count
    !
    !> \brief Returns the number of streams containing the specified field
    !> \author Michael Duda
    !> \date   18 October 2016
    !> \details
    !>  This function queries the stream manager to determine how many streams
    !>  contain the specified field. The optional argument 'direction' can be
    !>  used to limit the count to only input streams, only output streams, or
    !>  both input and output streams. By default, the function only counts
    !>  streams that have a non-'none' input/output interval, but even streams
    !>  with input_interval or output_interval equal to 'none' can be considered
    !>  by setting the optional argument 'includeInactive' to .true..
    !
    !-----------------------------------------------------------------------
    integer function iamas_stream_inclusion_count(fieldName, direction, includeInactive)

        use iamas_kind_types, only : StrKIND
        use iamas_derived_types, only : iAMAS_STREAM_INPUT, iAMAS_STREAM_OUTPUT, iAMAS_STREAM_INPUT_OUTPUT, &
                                       iAMAS_STREAM_PROPERTY_RECORD_INTV
        use iamas_stream_manager, only : iamas_stream_mgr_begin_iteration, iamas_stream_mgr_get_next_stream, &
                                        iAMAS_stream_mgr_ringing_alarms, iamas_stream_mgr_get_next_field, &
                                        iamas_stream_mgr_get_property

        implicit none

        character(len=*), intent(in) :: fieldName
        integer, intent(in), optional :: direction
        logical, intent(in), optional :: includeInactive

        character (len=StrKIND) :: streamNameItr
        character (len=StrKIND) :: fieldNameItr
        character (len=StrKIND) :: recordIntervalIn
        character (len=StrKIND) :: recordIntervalOut
        integer :: streamDirection
        logical :: streamActive
        logical :: fieldActive
        integer :: ierr

        integer :: local_direction
        logical :: local_includeInactive


        if (present(direction)) then
            local_direction = direction
        else
            local_direction = iAMAS_STREAM_INPUT_OUTPUT
        end if

        if (present(includeInactive)) then
            local_includeInactive = includeInactive
        else
            local_includeInactive = .false.
        end if


        iamas_stream_inclusion_count = 0

        call iamas_stream_mgr_begin_iteration(streamManager)
        STREAM_LOOP: do while (iamas_stream_mgr_get_next_stream(streamManager, streamID = streamNameItr, &
                                                               directionProperty = streamDirection, activeProperty = streamActive))

            call iAMAS_stream_mgr_get_property(streamManager, trim(streamNameItr), iAMAS_STREAM_PROPERTY_RECORD_INTV, &
                                              recordIntervalIn, direction=iAMAS_STREAM_INPUT, ierr=ierr)

            call iAMAS_stream_mgr_get_property(streamManager, trim(streamNameItr), iAMAS_STREAM_PROPERTY_RECORD_INTV, &
                                              recordIntervalOut, direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

            ! Determine whether this stream is "active" for the purposes of consideration here
            if (.not. local_includeInactive) then 
                if (streamActive) then
                    streamActive = ((local_direction == iAMAS_STREAM_INPUT .and. trim(recordIntervalIn) /= 'none') &
                                    .or. (local_direction == iAMAS_STREAM_OUTPUT .and. trim(recordIntervalOut) /= 'none') &
                                    .or. (local_direction == iAMAS_STREAM_INPUT_OUTPUT .and. ((trim(recordIntervalIn) /= 'none') &
                                                                                             .or. (trim(recordIntervalOut) /= 'none'))))
                end if
            else
                streamActive = .true.
            end if

            if (streamActive .and. ((local_direction == iAMAS_STREAM_INPUT .and. streamDirection == iAMAS_STREAM_INPUT) &
                                    .or. (local_direction == iAMAS_STREAM_OUTPUT .and. streamDirection == iAMAS_STREAM_OUTPUT) &
                                    .or. (local_direction == iAMAS_STREAM_INPUT_OUTPUT) &
                                    .or. (streamDirection == iAMAS_STREAM_INPUT_OUTPUT))) then

                call iamas_stream_mgr_begin_iteration(streamManager, streamID=streamNameItr)
                do while (iamas_stream_mgr_get_next_field(streamManager, streamNameItr, fieldNameItr, isActive=fieldActive))

                    if (fieldActive .and. (fieldNameItr == fieldName)) then
                        iamas_stream_inclusion_count = iamas_stream_inclusion_count + 1
                        cycle STREAM_LOOP
                    end if

                end do

            end if

        end do STREAM_LOOP

    end function iamas_stream_inclusion_count


end module iamas_atm_diagnostics_utils
