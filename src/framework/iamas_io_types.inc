#ifdef USE_PIO2
   integer, parameter :: iAMAS_IO_OFFSET_KIND = PIO_OFFSET_KIND
#else
   integer, parameter :: iAMAS_IO_OFFSET_KIND = PIO_OFFSET
#endif

   ! File access modes
   integer, parameter :: iAMAS_IO_READ  = 1, &
                         iAMAS_IO_WRITE = 2

   ! I/O formats
   integer, parameter :: iAMAS_IO_NETCDF   = 3, &
                         iAMAS_IO_PNETCDF  = 4, &
                         iAMAS_IO_NETCDF4  = 5, &
                         iAMAS_IO_PNETCDF5 = 6

   ! Field and attribute types
   integer, parameter :: iAMAS_IO_REAL     =  7,  &
                         iAMAS_IO_DOUBLE   =  8,  &
                         iAMAS_IO_INT      =  9,  &
                         iAMAS_IO_LOGICAL  = 10,  &
                         iAMAS_IO_CHAR     = 11

   ! Field precision
   integer, parameter :: iAMAS_IO_SINGLE_PRECISION = 12, &
                         iAMAS_IO_DOUBLE_PRECISION = 13, &
#ifdef SINGLE_PRECISION
                         iAMAS_IO_NATIVE_PRECISION = iAMAS_IO_SINGLE_PRECISION
#else
                         iAMAS_IO_NATIVE_PRECISION = iAMAS_IO_DOUBLE_PRECISION
#endif

   ! Unlimited / record dimension
   integer, parameter :: iAMAS_IO_UNLIMITED_DIM = -123456

   ! Error codes
   integer, parameter :: iAMAS_IO_NOERR              =  0, &
                         iAMAS_IO_ERR_INVALID_MODE   = -1, &
                         iAMAS_IO_ERR_INVALID_FORMAT = -2, &
                         iAMAS_IO_ERR_LONG_FILENAME  = -3, &
                         iAMAS_IO_ERR_UNINIT_HANDLE  = -4, &
                         iAMAS_IO_ERR_PIO            = -5, &
                         iAMAS_IO_ERR_DATA_MODE      = -6, &
                         iAMAS_IO_ERR_NOWRITE        = -7, &
                         iAMAS_IO_ERR_REDEF_DIM      = -8, &
                         iAMAS_IO_ERR_REDEF_VAR      = -9, &
                         iAMAS_IO_ERR_UNDEFINED_DIM  = -10, &
                         iAMAS_IO_ERR_UNDEFINED_VAR  = -11, &
                         iAMAS_IO_ERR_REDEF_ATT      = -12, &
                         iAMAS_IO_ERR_WRONG_ATT_TYPE = -13, &
                         iAMAS_IO_ERR_NO_DECOMP      = -14, &
                         iAMAS_IO_ERR_TWO_UNLIMITED_DIMS = -15, &
                         iAMAS_IO_ERR_WRONG_MODE         = -16, &
                         iAMAS_IO_ERR_NO_UNLIMITED_DIM   = -17, &
                         iAMAS_IO_ERR_UNIMPLEMENTED      = -18, &
                         iAMAS_IO_ERR_WOULD_CLOBBER      = -19, &
                         iAMAS_IO_ERR_NOEXIST_READ       = -20, &
                         iAMAS_IO_ERR_MISSING_DIM        = -21

   type iAMAS_IO_Handle_type
      logical :: initialized = .false.
      logical :: preexisting_file = .false.
      logical :: data_mode = .false.
      type (file_desc_t) :: pio_file
      character (len=StrKIND) :: filename
      integer :: iomode
      integer :: ioformat
      integer :: pio_unlimited_dimid
      integer :: preexisting_records = 0
      integer (kind=iAMAS_IO_OFFSET_KIND) :: frame_number = 1
      type (dimlist_type), pointer :: dimlist_head => null()
      type (dimlist_type), pointer :: dimlist_tail => null()
      type (fieldlist_type), pointer :: fieldlist_head => null()
      type (fieldlist_type), pointer :: fieldlist_tail => null()
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (iamas_io_context_type), pointer :: ioContext => null()
   end type iAMAS_IO_Handle_type

   type decomphandle_type
      integer :: field_type
      integer, dimension(:), pointer :: dims
      integer, dimension(:), pointer :: indices
      type (io_desc_t) :: pio_iodesc
   end type decomphandle_type

   type atthandle_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer :: precision
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
   end type atthandle_type

   type dimhandle_type
      character (len=StrKIND) :: dimname
      logical :: is_unlimited_dim = .false.
      integer :: dimsize
      integer :: dimid
   end type dimhandle_type

   type fieldhandle_type
      character (len=StrKIND) :: fieldname
      integer :: fieldid
      type (Var_desc_t) :: field_desc
      integer :: field_type
      logical :: has_unlimited_dim = .false.
      integer :: ndims
      integer :: precision
      type (dimhandle_type), pointer, dimension(:) :: dims
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (decomphandle_type), pointer :: decomp => null()

      type (decomphandle_type), pointer :: decomp_mul => null()
   end type fieldhandle_type

   type decomplist_type
      type (decomphandle_type), pointer :: decomphandle
      type (decomplist_type), pointer :: next => null()
   end type decomplist_type

   type attlist_type
      type (atthandle_type), pointer :: atthandle
      type (attlist_type), pointer :: next => null()
   end type attlist_type

   type dimlist_type
      type (dimhandle_type), pointer :: dimhandle
      type (dimlist_type), pointer :: next => null()
   end type dimlist_type

   type fieldlist_type
      type (fieldhandle_type), pointer :: fieldhandle
      type (fieldlist_type), pointer :: next => null()
   end type fieldlist_type

   type iamas_io_context_type
      type (decomplist_type), pointer :: decomp_list => null()
      type (iosystem_desc_t), pointer :: pio_iosystem => null()
      integer :: master_pio_iotype = -999
      type (dm_info), pointer :: dminfo => null()
      !added for multiple files
      type (decomplist_type), pointer :: decomp_list_mul => null()
      type (iosystem_desc_t), pointer :: pio_iosystem_mul => null()
      logical, pointer :: io_use_mulfiles => null()
      logical, pointer :: do_output_indices => null()
      logical, pointer :: is_init_atm => null()
   end type iamas_io_context_type
