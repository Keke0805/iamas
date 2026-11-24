   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
      type (iamas_pool_type), pointer :: configs, packages

      type (iAMAS_Clock_type), pointer :: clock
      type (iAMAS_log_type), pointer :: logInfo => null()
      type (iAMAS_streamManager_type), pointer :: streamManager
      type (iamas_decomp_list), pointer :: decompositions => null()
      type (iamas_io_context_type), pointer :: ioContext => null()

      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo

      ! Store exchange group information here
      type (iamas_exchange_group), pointer :: exchangeGroups => null()

      ! Domain specific constants
      logical :: on_a_sphere = .true.
      logical :: is_periodic = .false.
      real (kind=RKIND) :: sphere_radius = 1.0_RKIND
      real (kind=RKIND) :: x_period = 0.0_RKIND
      real (kind=RKIND) :: y_period = 0.0_RKIND
      character (len=StrKIND) :: namelist_filename = '' !< Constant: Name of namelist file
      character (len=StrKIND) :: streams_filename = '' !< Constant: Name of stream configuration file
      character (len=StrKIND) :: mesh_spec = '' !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id = '' !< parent_id attribute, read in from input file.

      ! Pointer to timer root
      type (iamas_timer_root), pointer :: timer_root => null()

      ! Back pointer to core
      type (core_type), pointer :: core => null()

      ! Domain_type is a linked list
      type (domain_type), pointer :: next => null()
   end type domain_type
