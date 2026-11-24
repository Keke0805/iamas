   integer, parameter :: iAMAS_MAX_ALARMS = 40

   integer, parameter :: iAMAS_NOW = 0, &
                         iAMAS_START_TIME = 1, &
                         iAMAS_STOP_TIME = 2

   integer, parameter :: iAMAS_FORWARD = 1, &
                         iAMAS_BACKWARD = -1

   integer, parameter :: iAMAS_GREGORIAN = 0, &
                         iAMAS_GREGORIAN_NOLEAP = 1, &
                         iAMAS_360DAY = 2

   type iAMAS_Time_type
      type (ESMF_Time) :: t
   end type

   type iAMAS_TimeInterval_type
      type (ESMF_TimeInterval) :: ti
   end type

   type iAMAS_Alarm_type
      character (len=ShortStrKIND) :: alarmID
      logical :: isRecurring
      logical :: isSet
      type (iAMAS_Time_type) :: ringTime
      type (iAMAS_Time_type) :: prevRingTime
      type (iAMAS_TimeInterval_type) :: ringTimeInterval
      type (iAMAS_Alarm_type), pointer :: next => null()
   end type
   
   type iAMAS_Clock_type
      integer :: direction
      integer :: nAlarms
      type (ESMF_Clock) :: c
      type (iAMAS_Alarm_type), pointer :: alarmListHead => null()
   end type

