   integer, parameter :: iAMAS_POOL_TABLE_SIZE = 128

   integer, parameter :: iAMAS_POOL_SILENT = 1001, &
                         iAMAS_POOL_WARN   = 1002, &
                         iAMAS_POOL_FATAL  = 1003

   integer, parameter :: iAMAS_POOL_FIELD     = 1004, &
                         iAMAS_POOL_CONFIG    = 1005, &
                         iAMAS_POOL_DIMENSION = 1006, &
                         iAMAS_POOL_SUBPOOL   = 1007, &
                         iAMAS_POOL_PACKAGE   = 1008

   integer, parameter :: iAMAS_POOL_REAL      = 1009, &
                         iAMAS_POOL_INTEGER   = 1010, &
                         iAMAS_POOL_LOGICAL   = 1011, &
                         iAMAS_POOL_CHARACTER = 1012

   type iamas_pool_data_type
      integer :: contentsType
      integer :: contentsDims
      integer :: contentsTimeLevs

      ! For storing fields
      type (field0DReal), pointer :: r0 => null()
      type (field1DReal), pointer :: r1 => null()
      type (field2DReal), pointer :: r2 => null()
      type (field3DReal), pointer :: r3 => null()
      type (field4DReal), pointer :: r4 => null()
      type (field5DReal), pointer :: r5 => null()
      type (field0DReal), dimension(:), pointer :: r0a => null()
      type (field1DReal), dimension(:), pointer :: r1a => null()
      type (field2DReal), dimension(:), pointer :: r2a => null()
      type (field3DReal), dimension(:), pointer :: r3a => null()
      type (field4DReal), dimension(:), pointer :: r4a => null()
      type (field5DReal), dimension(:), pointer :: r5a => null()
      type (field0DInteger), pointer :: i0 => null()
      type (field1DInteger), pointer :: i1 => null()
      type (field2DInteger), pointer :: i2 => null()
      type (field3DInteger), pointer :: i3 => null()
      type (field0DInteger), dimension(:), pointer :: i0a => null()
      type (field1DInteger), dimension(:), pointer :: i1a => null()
      type (field2DInteger), dimension(:), pointer :: i2a => null()
      type (field3DInteger), dimension(:), pointer :: i3a => null()
      type (field0DChar), pointer :: c0 => null()
      type (field1DChar), pointer :: c1 => null()
      type (field0DChar), dimension(:), pointer :: c0a => null()
      type (field1DChar), dimension(:), pointer :: c1a => null()
      type (field0DLogical), pointer :: l0 => null()
      type (field0DLogical), dimension(:), pointer :: l0a => null()
      type (iamas_pool_type), pointer :: p => null()
 
      ! For storing config options, dimensions, and packages
      integer, pointer :: simple_int => null()
      integer, dimension(:), pointer :: simple_int_arr => null()
      real(kind=RKIND), pointer :: simple_real => null()
      logical, pointer :: simple_logical => null()
      character(len=StrKIND), pointer :: simple_char => null()
   end type iamas_pool_data_type

   type iamas_pool_member_type
      character (len=StrKIND) :: key
      integer :: keyLen
      integer :: contentsType
      type (iamas_pool_data_type), pointer :: data => null()
      type (iamas_pool_member_type), pointer :: next => null()
      type (iamas_pool_member_type), pointer :: iteration_next => null()
      type (iamas_pool_member_type), pointer :: iteration_prev => null()
   end type iamas_pool_member_type

   type iamas_pool_head_type
      type (iamas_pool_member_type), pointer :: head => null()
   end type iamas_pool_head_type

   type iamas_pool_type
      integer :: size
      type (iamas_pool_head_type), dimension(:), pointer :: table => null()
      type (iamas_pool_member_type), pointer :: iterator => null()
      type (iamas_pool_member_type), pointer :: iteration_head => null()
      type (iamas_pool_member_type), pointer :: iteration_tail => null()
   end type iamas_pool_type

   type iamas_pool_iterator_type
      character (len=StrKIND) :: memberName
      integer :: memberType
      integer :: dataType
      integer :: nDims
      integer :: nTimeLevels
   end type iamas_pool_iterator_type

   type iamas_pool_field_info_type
      integer :: fieldType
      integer :: nDims
      integer :: nTimeLevels
      integer :: nHaloLayers
      logical :: isActive
   end type iamas_pool_field_info_type
