    ! main particle type (used by particle_framework)
    type iamas_particle_type
       ! IO block number
       !integer :: ioBlock = 0
       ! pool containing particle data needed for communication
       ! between computational cells
       type (iamas_pool_type), pointer :: haloDataPool => NULL()
       ! pool containing particle data not needing to be passed
       ! between processors
       type (iamas_pool_type), pointer :: nonhaloDataPool => NULL()
    end type iamas_particle_type

    ! just a list of particles (main unit that code operates on)
    type iamas_particle_list_type
       type (iamas_particle_type), pointer :: particle => NULL()
       ! next item in the list
       type (iamas_particle_list_type), pointer :: next => NULL()
       ! previous item in the list
       type (iamas_particle_list_type), pointer :: prev => NULL()
    end type iamas_particle_list_type

    type iamas_list_of_particle_list_type
       type (iamas_particle_list_type), pointer :: list => NULL()
    end type iamas_list_of_particle_list_type
