        type iamas_timer_root
           type (iamas_timer_node), pointer :: root_timer => null()
           type (iamas_timer_node), pointer :: current_timer => null()
           type (dm_info), pointer :: dminfo => null()
        end type iamas_timer_root

        type iamas_timer_node
          character (len=StrKIND) :: timer_name
          logical :: printed
          integer :: nlen
          logical, pointer, dimension(:) :: running
          integer :: calls
          real (kind=R8KIND), pointer, dimension(:) :: start_time, end_time, total_time
          real (kind=RKIND), pointer, dimension(:) :: max_time, min_time, avg_time
          type (iamas_timer_node), pointer :: next => null()
          type (iamas_timer_node), pointer :: child => null()
          type (iamas_timer_node), pointer :: parent => null()
        end type iamas_timer_node
