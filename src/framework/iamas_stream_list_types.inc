    integer, parameter :: iAMAS_STREAM_LIST_NOERR     = 0, &
                          iAMAS_STREAM_LIST_DUPLICATE = 1, &
                          iAMAS_STREAM_LIST_NOT_FOUND = 2
    
    ! jwfeng
    integer, public, parameter :: maxVarArrays = 32    , &
                                  maxVarsInArray = 256

    type iAMAS_arrayField_list_type
        ! jwfeng from AEMOL: 06/09/2020
        ! Add these vars to customly determine which vars in array fields should be outputed
        character(len=StrKIND) :: arrayFields(maxVarArrays)
        character(len=StrKIND) :: varsInArray(maxVarsInArray, maxVarArrays)
    end type iAMAS_arrayField_list_type

    type iAMAS_stream_list_type

        ! Used by list head
        integer :: nItems = 0
        type (iAMAS_stream_list_type), pointer :: head => null()

        ! Used by streams
        integer :: direction
        logical :: valid = .false.
        logical :: immutable = .false.
        logical :: active_stream = .true.
        logical :: blockWrite = .false.
        character(len=StrKIND) :: filename
        character(len=StrKIND) :: filename_template
        character(len=StrKIND) :: filename_interval
        type (iAMAS_Stream_type), pointer :: stream => null()
        integer :: timeLevel = 0
        integer :: nRecords
        integer :: precision = iAMAS_IO_NATIVE_PRECISION
        integer :: clobber_mode
        integer :: io_type
        type (iAMAS_TimeInterval_type), pointer :: recordInterval => null()
        type (iAMAS_stream_list_type), pointer :: alarmList_in => null()
        type (iAMAS_stream_list_type), pointer :: alarmList_out => null()
        type (iamas_pool_type), pointer :: att_pool => null()
        type (iamas_pool_type), pointer :: field_pool => null()
        type (iamas_pool_type), pointer :: field_pkg_pool => null()
        type (iamas_pool_type), pointer :: pkg_pool => null()
        type (iAMAS_Time_type), pointer :: referenceTime => null()

        ! Used by alarms
        type (iAMAS_stream_list_type), pointer :: streamList => null()

        ! Used by streams and alarms
        character(len=StrKIND) :: name
        type (iAMAS_stream_list_type), pointer :: xref => null()
        type (iAMAS_stream_list_type), pointer :: next => null()

        type (iAMAS_arrayField_list_type) :: arrayFieldList

    end type iAMAS_stream_list_type

