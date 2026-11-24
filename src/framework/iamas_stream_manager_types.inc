    integer, public, parameter :: iAMAS_STREAM_ERR_FATAL  = 1000, &
                                  iAMAS_STREAM_ERR_WARN   = 1001, &
                                  iAMAS_STREAM_ERR_SILENT = 1002

    integer, public, parameter :: iAMAS_STREAM_MGR_NOERR            =  0, &
                                  iAMAS_STREAM_MGR_ERR_CLOBBER_FILE = -1, &
                                  iAMAS_STREAM_MGR_ERR_CLOBBER_REC  = -2, &
                                  iAMAS_STREAM_MGR_ERROR            = -3

    integer, public, parameter :: iAMAS_STREAM_INPUT        = 1, &
                                  iAMAS_STREAM_OUTPUT       = 2, &
                                  iAMAS_STREAM_INPUT_OUTPUT = 3, &
                                  iAMAS_STREAM_NONE         = 4

    integer, public, parameter :: iAMAS_STREAM_PROPERTY_ACTIVE        =  5, &
                                  iAMAS_STREAM_PROPERTY_IMMUTABLE     =  6, &
                                  iAMAS_STREAM_PROPERTY_FILENAME      =  7, &
                                  iAMAS_STREAM_PROPERTY_REF_TIME      =  8, &
                                  iAMAS_STREAM_PROPERTY_RECORD_INTV   =  9, &
                                  iAMAS_STREAM_PROPERTY_PRECISION     = 10, &
                                  iAMAS_STREAM_PROPERTY_FILENAME_INTV = 11, &
                                  iAMAS_STREAM_PROPERTY_CLOBBER       = 12, &
                                  iAMAS_STREAM_PROPERTY_IOTYPE        = 13

    integer, public, parameter :: iAMAS_STREAM_CLOBBER_NEVER     = 100, &
                                  iAMAS_STREAM_CLOBBER_APPEND    = 101, &
                                  iAMAS_STREAM_CLOBBER_TRUNCATE  = 102, &
                                  iAMAS_STREAM_CLOBBER_OVERWRITE = 103



    type iAMAS_streamManager_type

        integer :: numStreams = 0
        integer :: errorLevel

        type (iAMAS_IO_context_type), pointer :: ioContext
        type (iAMAS_Clock_type), pointer :: streamClock
        type (iAMAS_Pool_type), pointer :: allFields
        type (iAMAS_Pool_type), pointer :: allPackages
        type (iAMAS_Pool_type), pointer :: allStructs
        type (iAMAS_Pool_type), pointer :: defaultAtts

        type (iAMAS_stream_list_type), pointer :: currentStream => NULL()
        type (iAMAS_stream_list_type), pointer :: streams
        type (iAMAS_stream_list_type), pointer :: alarms_in
        type (iAMAS_stream_list_type), pointer :: alarms_out

    end type iAMAS_streamManager_type


