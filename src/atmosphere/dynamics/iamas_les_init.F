module iamas_les_init

   use iamas_derived_types, only : domain_type, block_type
   use iamas_kind_types, only : RKIND, StrKIND
   use iamas_pool_routines
   use iamas_log, only : iamas_log_write

   implicit none
   private

   public :: iamas_les_initialize

contains

   subroutine iamas_les_initialize(domain)
      type(domain_type), intent(inout) :: domain

      type(block_type), pointer :: block
      type(iamas_pool_type), pointer :: state, mesh
      real(kind=RKIND), pointer :: affine_condition(:)
      character(len=StrKIND), pointer :: config_les_scheme
      real(kind=RKIND), pointer :: config_les_tke_min
      logical, pointer :: config_do_restart
      integer, pointer :: index_tke_sgs
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars
      integer :: time_lev

      nullify(config_les_scheme, config_les_tke_min, config_do_restart)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_les_scheme', config_les_scheme)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_les_tke_min', config_les_tke_min)
      call iamas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return
      if (.not. associated(config_les_tke_min)) return
      if (.not. associated(config_do_restart)) return

      ! Geometry caches are reconstructed on every cold/restart initialization.
      ! Prognostic state is untouched; each block owns its mesh cache.
      block => domain % blocklist
      do while (associated(block))
         nullify(mesh,affine_condition)
         call iamas_pool_get_subpool(block % structs,'mesh',mesh)
         call iamas_pool_get_array(mesh,'les_affine_condition',affine_condition)
         if (.not.associated(affine_condition)) error stop 'C affine cache unavailable at initialization'
         affine_condition=0.0_RKIND
         block => block % next
      end do

      if (config_do_restart) then
         ! The restart stream has already restored the prognostic state into
         ! the main time level.  Leave every stored value untouched here;
         ! atm_rk_integration_setup refreshes its work level from that state.
         ! Non-finite or out-of-range values belong in the read-only restart
         ! QC gate, not in a silent mutation during model initialization.
         call iamas_log_write('Preserved LES SGS TKE exactly as read from restart')
         return
      end if

      block => domain % blocklist
      do while (associated(block))
         nullify(state, index_tke_sgs)
         call iamas_pool_get_subpool(block % structs, 'state', state)
         call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)

         if (associated(index_tke_sgs)) then
            if (index_tke_sgs > 0) then
               do time_lev = 1, 2
                  nullify(scalars)
                  call iamas_pool_get_array(state, 'scalars', scalars, time_lev)
                  if (associated(scalars)) then
                     scalars(index_tke_sgs,:,:) = config_les_tke_min
                  end if
               end do
            end if
         end if

         block => block % next
      end do

      call iamas_log_write('Initialized cold-start LES SGS TKE to $r', &
           realArgs=(/config_les_tke_min/))
   end subroutine iamas_les_initialize

end module iamas_les_init
