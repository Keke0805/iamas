   integer, parameter :: iAMAS_DMPAR_NOERR = 0
   integer, parameter :: iAMAS_DMPAR_MISSING_GROUP = 1
   integer, parameter :: iAMAS_DMPAR_EXISTING_GROUP = 2
   integer, parameter :: iAMAS_DMPAR_MISSING_FIELD = 3
   integer, parameter :: iAMAS_DMPAR_FIELD_TIMELEVEL_ERR = 4
   integer, parameter :: iAMAS_DMPAR_FIELD_HALO_ERR = 5
   integer, parameter :: iAMAS_DMPAR_BUFFER_EXISTS = 6

   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: initialized_mpi

     ! For multi-file PIO
     integer :: local_nprocs, local_proc_id, local_comm
     integer :: io_nprocs, io_proc_id, io_comm
     ! Add variables specific to block decomposition. {{{
     ! These are used in iamas_block_decomp.F
     integer :: total_blocks
     logical :: explicitDecomp
     integer, dimension(:), allocatable :: block_proc_list
     integer, dimension(:), allocatable :: block_local_id_list
     !}}}
   end type dm_info


   type iamas_exchange_list
     integer :: endPointID
     integer :: nlist
     integer, dimension(:), pointer :: srcList
     integer, dimension(:), pointer :: destList
     type (iamas_exchange_list), pointer :: next
   end type iamas_exchange_list


   type iamas_exchange_list_pointer
     type (iamas_exchange_list), pointer :: exchList
   end type iamas_exchange_list_pointer


   type iamas_multihalo_exchange_list
     type (iamas_exchange_list_pointer), dimension(:), pointer :: halos
     ! Pointers to the mulithalo exchange lists for this variable on the prev and next blocks on this processor
     type (iamas_multihalo_exchange_list), pointer :: prev, next
   end type iamas_multihalo_exchange_list


   type iamas_communication_list
     integer :: procID
     integer :: nlist
     integer :: bufferOffset
     real (kind=RKIND), dimension(:), pointer :: rbuffer => null()
     integer, dimension(:), pointer :: ibuffer => null()
     integer :: reqID
     type (iamas_communication_list), pointer :: next => null()
     integer :: commListSize
     logical :: received
     logical :: unpacked
   end type iamas_communication_list

   type iamas_exchange_field_list
      integer :: nLen
      character (len=StrKIND) :: fieldName
      integer :: fieldType, nDims
      logical, dimension(:), pointer :: haloLayers => null()
      logical, dimension(:), pointer :: timeLevels => null()
      type (iamas_exchange_field_list), pointer :: next => null()
   end type iamas_exchange_field_list

   type iamas_exchange_group
      integer :: nLen
      character (len=StrKIND) :: groupName
      type (iamas_pool_type), pointer :: fieldPool => null()
      type (iamas_exchange_field_list), pointer :: fieldList => null()

      ! Each group should have it's own communication lists, since they have the buffers for the communications
      type (iamas_communication_list), pointer :: sendList => null()
      type (iamas_communication_list), pointer :: recvList => null()
      type (iamas_exchange_group), pointer :: next => null()
      logical :: allReceived
   end type iamas_exchange_group
