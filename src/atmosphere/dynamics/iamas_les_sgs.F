! C3 shared-stress kernel; independently verified source embedded in this candidate.
module iamas_les_cross_kernel
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: cross_prepare, cross_gather
contains
  subroutine cross_prepare(nc,ne,nz,cells,gids,z,dc,dv,h,coef,volume,offset,elist,slist,ierr) bind(C)
    integer(c_int), value :: nc,ne,nz
    integer(c_int), intent(in) :: cells(2,ne),gids(ne)
    real(c_double), intent(in) :: z(nz+1,nc),dc(ne),dv(ne)
    real(c_double), intent(out) :: h(nz-1,ne),coef(3,nz-1,2,ne),volume(nz-1,ne)
    integer(c_int), intent(out) :: offset(nc+1),elist(2*ne),slist(2*ne),ierr
    integer :: e,c,c1,c2,s,j,p,q,hold_e,hold_s
    integer, allocatable :: counts(:),cursor(:)
    real(c_double) :: low,high,a,b,slope,sgn
    ierr=1
    if (nc<2 .or. ne<1 .or. nz<3) return
    ierr=2
    if (.not.all(ieee_is_finite(z)) .or. .not.all(ieee_is_finite(dc)) .or. &
        .not.all(ieee_is_finite(dv))) return
    if (any(dc<=0) .or. any(dv<=0) .or. any(z(2:nz+1,:)<=z(1:nz,:))) return
    ierr=3
    if (any(cells<1) .or. any(cells>nc) .or. any(gids<=0)) return
    if (any(cells(1,:)==cells(2,:))) return
    allocate(counts(nc),cursor(nc));counts=0
    do e=1,ne
      c1=cells(1,e);c2=cells(2,e)
      counts(c1)=counts(c1)+1;counts(c2)=counts(c2)+1
      do j=1,nz-1
        low=.5_c_double*(.5_c_double*(z(j,c1)+z(j+1,c1))+.5_c_double*(z(j,c2)+z(j+1,c2)))
        high=.5_c_double*(.5_c_double*(z(j+1,c1)+z(j+2,c1))+.5_c_double*(z(j+1,c2)+z(j+2,c2)))
        h(j,e)=high-low
        volume(j,e)=dc(e)*dv(e)*h(j,e)
        slope=(z(j+1,c2)-z(j+1,c1))/dc(e)
        do s=1,2
          c=cells(s,e);a=z(j+1,c)-z(j,c);b=z(j+2,c)-z(j+1,c)
          coef(1,j,s,e)=-.5_c_double*slope*(-b/(a*(a+b)))
          coef(2,j,s,e)=-.5_c_double*slope*((b-a)/(a*b))
          coef(3,j,s,e)=-.5_c_double*slope*(a/(b*(a+b)))
          sgn=-1._c_double
          if (s==2) sgn=1._c_double
          coef(2,j,s,e)=coef(2,j,s,e)+sgn/dc(e)
        end do
      end do
    end do
    ierr=4
    if (any(h<=0) .or. .not.all(ieee_is_finite(coef)) .or. .not.all(ieee_is_finite(volume))) return
    offset(1)=1
    do c=1,nc
      offset(c+1)=offset(c)+counts(c)
    end do
    cursor=offset(1:nc)
    do e=1,ne
      do s=1,2
        c=cells(s,e);p=cursor(c);elist(p)=e;slist(p)=s;cursor(c)=p+1
      end do
    end do
    ! Sort each small cell star by global ID, independent of local numbering.
    do c=1,nc
      do p=offset(c)+1,offset(c+1)-1
        hold_e=elist(p);hold_s=slist(p);q=p-1
        do while(q>=offset(c))
          if (gids(elist(q))<=gids(hold_e)) exit
          elist(q+1)=elist(q);slist(q+1)=slist(q);q=q-1
        end do
        elist(q+1)=hold_e;slist(q+1)=hold_s
      end do
      do p=offset(c)+1,offset(c+1)-1
        if (gids(elist(p))==gids(elist(p-1))) return
      end do
    end do
    ierr=0
  end subroutine

  subroutine cross_gather(nc,ne,nz,noc,noe,cells,h,coef,volume,offset,elist,slist, &
      rho,u,w,kv,owned_c,owned_e,shear,fu,fw,prod,loss,ierr) bind(C)
    integer(c_int), value :: nc,ne,nz,noc,noe
    integer(c_int), intent(in) :: cells(2,ne),offset(nc+1),elist(2*ne),slist(2*ne)
    integer(c_int), intent(in) :: owned_c(noc),owned_e(noe)
    real(c_double), intent(in) :: h(nz-1,ne),coef(3,nz-1,2,ne),volume(nz-1,ne)
    real(c_double), intent(in) :: rho(nz,nc),u(nz,ne),w(nz+1,nc),kv(nz-1,ne)
    real(c_double), intent(out) :: shear(nz-1,ne),fu(nz,ne),fw(nz+1,nc),prod(nz,nc),loss
    integer(c_int), intent(out) :: ierr
    real(c_double), allocatable :: impulse(:,:),energy(:,:)
    real(c_double) :: rho_low,rho_high,weight,nan
    integer :: e,c,c1,c2,j,s,t,p,owner,level,center
    logical, allocatable :: seen_c(:),seen_e(:)
    ierr=1;loss=0
    if (nc<2 .or. ne<1 .or. nz<3 .or. noc<0 .or. noe<0) return
    nan=ieee_value(0._c_double,ieee_quiet_nan);fu=nan;fw=nan;prod=nan;shear=nan
    ierr=2
    if (.not.all(ieee_is_finite(rho)) .or. .not.all(ieee_is_finite(u)) .or. &
        .not.all(ieee_is_finite(w)) .or. .not.all(ieee_is_finite(kv))) return
    if (any(rho<=0) .or. any(kv<0)) return
    ierr=3
    if (any(owned_c<1) .or. any(owned_c>nc) .or. any(owned_e<1) .or. any(owned_e>ne)) return
    allocate(seen_c(nc),seen_e(ne));seen_c=.false.;seen_e=.false.
    do owner=1,noc
      c=owned_c(owner)
      if (seen_c(c)) return
      seen_c(c)=.true.
    end do
    do owner=1,noe
      e=owned_e(owner)
      if (seen_e(e)) return
      seen_e(e)=.true.
    end do
    allocate(impulse(nz-1,ne),energy(nz-1,ne))
    do e=1,ne
      c1=cells(1,e);c2=cells(2,e)
      do j=1,nz-1
        shear(j,e)=(u(j+1,e)-u(j,e))/h(j,e)
        do s=1,2
          c=cells(s,e)
          do t=1,3
            shear(j,e)=shear(j,e)+coef(t,j,s,e)*w(j+t-1,c)
          end do
        end do
        rho_low=.5_c_double*(rho(j,c1)+rho(j,c2))
        rho_high=.5_c_double*(rho(j+1,c1)+rho(j+1,c2))
        weight=volume(j,e)*(.5_c_double*(rho_low+rho_high))*kv(j,e)
        impulse(j,e)=weight*shear(j,e);energy(j,e)=impulse(j,e)*shear(j,e)
      end do
    end do
    ierr=4
    if (.not.all(ieee_is_finite(energy))) return
    do owner=1,noe
      e=owned_e(owner)
      fu(1,e)=impulse(1,e)/h(1,e)
      do j=2,nz-1
        fu(j,e)=impulse(j,e)/h(j,e)-impulse(j-1,e)/h(j-1,e)
      end do
      fu(nz,e)=-impulse(nz-1,e)/h(nz-1,e)
    end do
    do owner=1,noc
      c=owned_c(owner);fw(:,c)=0;prod(:,c)=0
      do p=offset(c),offset(c+1)-1
        e=elist(p);s=slist(p)
        prod(1:nz-1,c)=prod(1:nz-1,c)+.25_c_double*energy(:,e)
        prod(2:nz,c)=prod(2:nz,c)+.25_c_double*energy(:,e)
        do level=1,nz+1
          do center=max(2,level-1),min(nz,level+1)
            fw(level,c)=fw(level,c)-impulse(center-1,e)*coef(level-center+2,center-1,s,e)
          end do
        end do
      end do
    end do
    loss=sum(energy);ierr=0
  end subroutine
end module

module iamas_les_sgs

   use iamas_kind_types, only : RKIND, StrKIND
   use iamas_constants, only : gravity, rvord, cp
   use iamas_atm_boundaries, only : nRelaxZone
   use iamas_pool_routines
   use iamas_les_cross_kernel, only: cross_prepare, cross_gather

   implicit none
   private

   public :: iamas_les_sgs_stage1_interface
   public :: iamas_les_compute_coefficients
   public :: iamas_les_add_sgs_tendencies
   public :: iamas_les_enforce_tke_floor
   public :: iamas_les_apply_lbc_tke_minimum
   public :: iamas_les_prepare_lbc_tke_driving

   ! Geometry is immutable; caches are keyed by the actual block zgrid target
   ! and owned ranges. The current production build is serial within each rank.
   type les_cross_cache_type
      real(kind=RKIND), pointer :: zkey(:,:) => null()
      integer :: cs,ce,es,ee,nc,ne,nz
      integer, allocatable :: cm(:),em(:),cells(:,:),gids(:),oc(:),oe(:),offset(:),elist(:),slist(:)
      real(kind=RKIND), allocatable :: h(:,:),coef(:,:,:,:),volume(:,:),dc(:),dv(:),ru(:),rv(:)
      real(kind=RKIND), allocatable :: rho(:,:),u(:,:),w(:,:),kv(:,:),shear(:,:),fu(:,:),fw(:,:),prod(:,:)
      type(les_cross_cache_type), pointer :: next => null()
   end type
   type(les_cross_cache_type), pointer, save :: les_cross_head => null()

contains

   subroutine iamas_les_sgs_stage1_interface(dt)
      real(kind=RKIND), intent(in) :: dt

      ! Stage-1 link interface. The prognostic-TKE closure and SGS
      ! tendencies are connected in the following implementation stages.
      if (dt < 0.0_RKIND) error stop 'Invalid LES timestep'
   end subroutine iamas_les_sgs_stage1_interface

   pure real(kind=RKIND) function iamas_les_safe_finite(value) result(safe_value)
      real(kind=RKIND), intent(in) :: value

      if (value /= value .or. abs(value) > 1.0e30_RKIND) then
         safe_value = 0.0_RKIND
      else
         safe_value = value
      end if
   end function iamas_les_safe_finite

   pure real(kind=RKIND) function iamas_les_safe_nonnegative(value) result(safe_value)
      real(kind=RKIND), intent(in) :: value

      if (value /= value .or. value < 0.0_RKIND .or. abs(value) > 1.0e30_RKIND) then
         safe_value = 0.0_RKIND
      else
         safe_value = value
      end if
   end function iamas_les_safe_nonnegative

   subroutine iamas_les_enforce_tke_floor(state, configs, nVertLevels, cellStart, cellEnd)
      type(iamas_pool_type), intent(in) :: state
      type(iamas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels, cellStart, cellEnd

      character(len=StrKIND), pointer :: config_les_scheme
      real(kind=RKIND), pointer :: config_les_tke_min
      integer, pointer :: index_tke_sgs
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_old
      integer :: iCell, k

      if (nVertLevels < 1) return

      nullify(config_les_scheme, config_les_tke_min, index_tke_sgs, scalars, scalars_old)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return

      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)
      call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)
      call iamas_pool_get_array(state, 'scalars', scalars, 2)
      call iamas_pool_get_array(state, 'scalars', scalars_old, 1)
      if (.not. associated(config_les_tke_min)) return
      if (.not. associated(index_tke_sgs)) return
      if (index_tke_sgs <= 0) return

      if (associated(scalars)) then
         do iCell = cellStart, cellEnd
            do k = 1, nVertLevels
               scalars(index_tke_sgs,k,iCell) = max(iamas_les_safe_nonnegative(scalars(index_tke_sgs,k,iCell)), &
                                                    config_les_tke_min)
            end do
         end do
      end if
      if (associated(scalars_old)) then
         do iCell = cellStart, cellEnd
            do k = 1, nVertLevels
               scalars_old(index_tke_sgs,k,iCell) = max(iamas_les_safe_nonnegative(scalars_old(index_tke_sgs,k,iCell)), &
                                                        config_les_tke_min)
            end do
         end do
      end if
   end subroutine iamas_les_enforce_tke_floor

   subroutine iamas_les_apply_lbc_tke_minimum(state, mesh, configs, nVertLevels, cellStart, cellEnd)
      type(iamas_pool_type), intent(in) :: state
      type(iamas_pool_type), intent(in) :: mesh
      type(iamas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels, cellStart, cellEnd

      character(len=StrKIND), pointer :: config_les_scheme, config_les_lbc_tke_mode
      real(kind=RKIND), pointer :: config_les_tke_min
      integer, pointer :: index_tke_sgs
      integer, dimension(:), pointer :: bdyMaskCell
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_old
      integer :: iCell

      if (nVertLevels < 1) return

      nullify(config_les_scheme, config_les_lbc_tke_mode, config_les_tke_min)
      nullify(index_tke_sgs, bdyMaskCell, scalars, scalars_old)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      call iamas_pool_get_config(configs, 'config_les_lbc_tke_mode', config_les_lbc_tke_mode)
      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)
      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return
      if (.not. associated(config_les_lbc_tke_mode)) return
      if (trim(config_les_lbc_tke_mode) /= 'minimum') return

      call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)
      call iamas_pool_get_array(state, 'scalars', scalars, 2)
      call iamas_pool_get_array(state, 'scalars', scalars_old, 1)
      call iamas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
      if (.not. associated(config_les_tke_min)) return
      if (.not. associated(index_tke_sgs)) return
      if (.not. associated(bdyMaskCell)) return
      if (index_tke_sgs <= 0) return

      do iCell = cellStart, cellEnd
         if (bdyMaskCell(iCell) > 0) then
            if (associated(scalars)) then
               scalars(index_tke_sgs,1:nVertLevels,iCell) = config_les_tke_min
            end if
            if (associated(scalars_old)) then
               scalars_old(index_tke_sgs,1:nVertLevels,iCell) = config_les_tke_min
            end if
         end if
      end do
   end subroutine iamas_les_apply_lbc_tke_minimum

   subroutine iamas_les_prepare_lbc_tke_driving(state, configs, scalars_driving, nVertLevels, nCells)
      type(iamas_pool_type), intent(in) :: state
      type(iamas_pool_type), intent(in) :: configs
      real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars_driving
      integer, intent(in) :: nVertLevels, nCells

      character(len=StrKIND), pointer :: config_les_scheme, config_les_lbc_tke_mode
      real(kind=RKIND), pointer :: config_les_tke_min
      integer, pointer :: index_tke_sgs

      if (nVertLevels < 1 .or. nCells < 1) return

      nullify(config_les_scheme, config_les_lbc_tke_mode, config_les_tke_min, index_tke_sgs)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      call iamas_pool_get_config(configs, 'config_les_lbc_tke_mode', config_les_lbc_tke_mode)
      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)
      call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)

      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return
      if (.not. associated(config_les_lbc_tke_mode)) return
      if (trim(config_les_lbc_tke_mode) /= 'minimum') return
      if (.not. associated(config_les_tke_min)) return
      if (.not. associated(index_tke_sgs)) return
      if (index_tke_sgs <= 0 .or. index_tke_sgs > size(scalars_driving, 1)) return

      scalars_driving(index_tke_sgs,1:nVertLevels,1:nCells+1) = config_les_tke_min
   end subroutine iamas_les_prepare_lbc_tke_driving

   subroutine iamas_les_compute_coefficients(state, diag, mesh, configs, nVertLevels, dt, &
                                              cellStart, cellEnd, rk_step)
      type(iamas_pool_type), intent(in) :: state
      type(iamas_pool_type), intent(in) :: diag
      type(iamas_pool_type), intent(in) :: mesh
      type(iamas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels, cellStart, cellEnd, rk_step
      real(kind=RKIND), intent(in) :: dt

      character(len=StrKIND), pointer :: config_les_scheme
      character(len=StrKIND), pointer :: config_les_coeff_mode
      real(kind=RKIND), pointer :: config_les_ck
      real(kind=RKIND), pointer :: config_les_prandtl
      real(kind=RKIND), pointer :: config_les_tke_min
      real(kind=RKIND), pointer :: config_les_mix_upper_bound
      integer, pointer :: index_qv, index_tke_sgs
      real(kind=RKIND), dimension(:), pointer :: areaCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid, theta_m
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars
      real(kind=RKIND), dimension(:,:), pointer :: les_mix_length
      real(kind=RKIND), dimension(:,:), pointer :: les_mix_length_h, les_mix_length_v
      real(kind=RKIND), dimension(:,:), pointer :: les_km_h, les_km_v
      real(kind=RKIND), dimension(:,:), pointer :: les_kh_h, les_kh_v
      real(kind=RKIND), dimension(:,:), pointer :: les_bv_freq2
      integer :: iCell, k
      real(kind=RKIND) :: dz, delta, theta_here, theta_low, theta_high
      real(kind=RKIND) :: z_low, z_high, dtheta_dz, n2
      real(kind=RKIND) :: tke, mix_length, mix_length_h, mix_length_v
      real(kind=RKIND) :: velocity_scale, pr_inv
      real(kind=RKIND) :: cap_h, cap_v

      if (rk_step /= 1) return
      if (nVertLevels < 2) return
      if (dt <= 0.0_RKIND) return

      nullify(config_les_scheme)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return

      call iamas_pool_get_config(configs, 'config_les_ck', config_les_ck)
      call iamas_pool_get_config(configs, 'config_les_prandtl', config_les_prandtl)
      call iamas_pool_get_config(configs, 'config_les_coeff_mode', config_les_coeff_mode)
      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)
      call iamas_pool_get_config(configs, 'config_les_mix_upper_bound', config_les_mix_upper_bound)
      if (.not. associated(config_les_ck) .or. .not. associated(config_les_prandtl)) return
      if (.not. associated(config_les_coeff_mode)) return
      if (.not. associated(config_les_tke_min) .or. .not. associated(config_les_mix_upper_bound)) return

      call iamas_pool_get_dimension(state, 'index_qv', index_qv)
      call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)
      if (.not. associated(index_qv) .or. .not. associated(index_tke_sgs)) return
      if (index_qv <= 0 .or. index_tke_sgs <= 0) return

      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(state, 'theta_m', theta_m, 2)
      call iamas_pool_get_array(state, 'scalars', scalars, 2)
      call iamas_pool_get_array(diag, 'les_mix_length', les_mix_length)
      call iamas_pool_get_array(diag, 'les_mix_length_h', les_mix_length_h)
      call iamas_pool_get_array(diag, 'les_mix_length_v', les_mix_length_v)
      call iamas_pool_get_array(diag, 'les_km_h', les_km_h)
      call iamas_pool_get_array(diag, 'les_km_v', les_km_v)
      call iamas_pool_get_array(diag, 'les_kh_h', les_kh_h)
      call iamas_pool_get_array(diag, 'les_kh_v', les_kh_v)
      call iamas_pool_get_array(diag, 'les_bv_freq2', les_bv_freq2)

      if (.not. associated(areaCell) .or. .not. associated(zgrid)) return
      if (.not. associated(theta_m) .or. .not. associated(scalars)) return
      if (.not. associated(les_mix_length) .or. .not. associated(les_km_h)) return
      if (.not. associated(les_mix_length_h) .or. .not. associated(les_mix_length_v)) return
      if (.not. associated(les_km_v) .or. .not. associated(les_kh_h)) return
      if (.not. associated(les_kh_v) .or. .not. associated(les_bv_freq2)) return

      do iCell = cellStart, cellEnd
         do k = 1, nVertLevels
            dz = max(abs(zgrid(k+1,iCell) - zgrid(k,iCell)), 1.0_RKIND)
            delta = max((areaCell(iCell) * dz)**(1.0_RKIND / 3.0_RKIND), 1.0_RKIND)

            theta_here = theta_m(k,iCell) / &
                         (1.0_RKIND + rvord * scalars(index_qv,k,iCell))

            if (k == 1) then
               theta_low = theta_here
               theta_high = theta_m(k+1,iCell) / &
                            (1.0_RKIND + rvord * scalars(index_qv,k+1,iCell))
               z_low = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
               z_high = 0.5_RKIND * (zgrid(k+1,iCell) + zgrid(k+2,iCell))
            else if (k == nVertLevels) then
               theta_low = theta_m(k-1,iCell) / &
                           (1.0_RKIND + rvord * scalars(index_qv,k-1,iCell))
               theta_high = theta_here
               z_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               z_high = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
            else
               theta_low = theta_m(k-1,iCell) / &
                           (1.0_RKIND + rvord * scalars(index_qv,k-1,iCell))
               theta_high = theta_m(k+1,iCell) / &
                            (1.0_RKIND + rvord * scalars(index_qv,k+1,iCell))
               z_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               z_high = 0.5_RKIND * (zgrid(k+1,iCell) + zgrid(k+2,iCell))
            end if

            dtheta_dz = (theta_high - theta_low) / max(abs(z_high - z_low), 1.0_RKIND)
            n2 = gravity * dtheta_dz / max(theta_here, 100.0_RKIND)
            tke = max(iamas_les_safe_nonnegative(scalars(index_tke_sgs,k,iCell)), config_les_tke_min)
            velocity_scale = sqrt(tke)

            mix_length = delta
            if (n2 > 1.0e-6_RKIND) then
               mix_length = min(delta, 0.76_RKIND * velocity_scale / sqrt(n2))
               mix_length = max(mix_length, 0.001_RKIND * delta)
            end if

            cap_h = config_les_mix_upper_bound * areaCell(iCell) / dt
            cap_v = config_les_mix_upper_bound * dz * dz / dt
            les_mix_length(k,iCell) = mix_length
            les_bv_freq2(k,iCell) = iamas_les_safe_finite(n2)

            select case (trim(config_les_coeff_mode))
            case ('legacy_volume')
               mix_length_h = mix_length
               mix_length_v = mix_length
               pr_inv = (1.0_RKIND + 2.0_RKIND * mix_length / delta) / &
                        max(config_les_prandtl, 1.0e-6_RKIND)
               les_km_h(k,iCell) = iamas_les_safe_nonnegative( &
                    min(config_les_ck * velocity_scale * mix_length, cap_h))
               les_km_v(k,iCell) = iamas_les_safe_nonnegative( &
                    min(config_les_ck * velocity_scale * mix_length, cap_v))
               les_kh_h(k,iCell) = iamas_les_safe_nonnegative( &
                    min(les_km_h(k,iCell) * pr_inv, cap_h))
               les_kh_v(k,iCell) = iamas_les_safe_nonnegative( &
                    min(les_km_v(k,iCell) * pr_inv, cap_v))

            case ('wrf_anisotropic')
               mix_length_h = max(sqrt(areaCell(iCell)), 1.0_RKIND)
               mix_length_v = dz
               if (n2 > 0.0_RKIND) then
                  mix_length_v = min(dz, 0.76_RKIND * velocity_scale / &
                                     sqrt(max(n2, 1.0e-30_RKIND)))
               end if
               mix_length_v = max(mix_length_v, 0.0_RKIND)

               les_km_h(k,iCell) = iamas_les_safe_nonnegative(min(max( &
                    config_les_ck * velocity_scale * mix_length_h, &
                    1.0e-6_RKIND * mix_length_h * mix_length_h), cap_h))
               les_km_v(k,iCell) = iamas_les_safe_nonnegative(min(max( &
                    config_les_ck * velocity_scale * mix_length_v, &
                    1.0e-6_RKIND * dz * dz), cap_v))
               les_kh_h(k,iCell) = iamas_les_safe_nonnegative( &
                    les_km_h(k,iCell) / max(config_les_prandtl, 1.0e-6_RKIND))
               les_kh_v(k,iCell) = iamas_les_safe_nonnegative( &
                    les_km_v(k,iCell) * (1.0_RKIND + 2.0_RKIND * mix_length_v / dz))

            case default
               error stop 'Unsupported config_les_coeff_mode'
            end select

            les_mix_length_h(k,iCell) = mix_length_h
            les_mix_length_v(k,iCell) = mix_length_v
         end do
      end do
   end subroutine iamas_les_compute_coefficients

   subroutine iamas_les_add_sgs_tendencies(tend, state, diag, diag_physics, mesh, configs, &
                                            nVertLevels, dt, rk_step, dynamics_substep, &
                                            cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)
      type(iamas_pool_type), intent(inout) :: tend
      type(iamas_pool_type), intent(in) :: state
      type(iamas_pool_type), intent(in) :: diag
      type(iamas_pool_type), intent(in) :: diag_physics
      type(iamas_pool_type), intent(in) :: mesh
      type(iamas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels, rk_step, dynamics_substep
      real(kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: edgeSolveStart, edgeSolveEnd

      character(len=StrKIND), pointer :: config_les_scheme
      character(len=StrKIND), pointer :: config_les_surface_mode
      character(len=StrKIND), pointer :: config_les_tke_diffusion_mode
      real(kind=RKIND), pointer :: config_les_ck
      real(kind=RKIND), pointer :: config_les_ce0
      real(kind=RKIND), pointer :: config_les_scalar_schmidt
      real(kind=RKIND), pointer :: config_les_tke_min
      real(kind=RKIND), pointer :: config_les_heat_flux
      real(kind=RKIND), pointer :: config_les_drag_coefficient
      logical, pointer :: config_les_output_budget
      integer, pointer :: num_scalars
      integer, pointer :: index_qv, index_tke_sgs
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      real(kind=RKIND), dimension(:), pointer :: invDcEdge, invDvEdge, dvEdge, rdzw
      real(kind=RKIND), dimension(:), pointer :: invAreaCell, meshScalingDel2
      real(kind=RKIND), dimension(:), pointer :: areaCell
      real(kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign, zgrid, zz
      real(kind=RKIND), dimension(:,:), pointer :: affine_d1u, affine_d1v, affine_d2u, affine_d2v
      real(kind=RKIND), dimension(:,:), pointer :: rho_zz, u, w, theta_m
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars
      real(kind=RKIND), dimension(:,:), pointer :: rho_edge, divergence, vorticity, v
      real(kind=RKIND), dimension(:,:), pointer :: ur_cell, vr_cell
      real(kind=RKIND), dimension(:,:), pointer :: les_km_h, les_km_v
      real(kind=RKIND), dimension(:,:), pointer :: les_kh_h, les_kh_v
      real(kind=RKIND), dimension(:,:), pointer :: les_mix_length, les_bv_freq2
      real(kind=RKIND), dimension(:,:), pointer :: les_shear_prod, les_buoy_prod
      real(kind=RKIND), dimension(:,:), pointer :: les_shear_prod_h
      real(kind=RKIND), dimension(:,:), pointer :: les_shear_prod_v_internal
      real(kind=RKIND), dimension(:,:), pointer :: les_surface_prod_current
      real(kind=RKIND), dimension(:,:), pointer :: les_surface_prod_stress_shear
      real(kind=RKIND), dimension(:,:), pointer :: les_surface_resolved_ke_sink_edge
      real(kind=RKIND), dimension(:,:), pointer :: les_dissipation, les_tke_diffusion
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_diffusion_h, les_tke_diffusion_v
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_h_flux, les_tke_v_flux
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_tend_total
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_source_pre_limiter
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_source_post_limiter
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_limiter_residual
      real(kind=RKIND), dimension(:,:), pointer :: les_tke_budget_closure_residual
      real(kind=RKIND), dimension(:,:), pointer :: les_tend_un, les_tend_w
      real(kind=RKIND), dimension(:,:), pointer :: les_tend_theta_m, les_tend_theta, les_tend_qv
      real(kind=RKIND), dimension(:,:), pointer :: les_unw_flux
      real(kind=RKIND), dimension(:,:), pointer :: les_uw_flux, les_vw_flux
      real(kind=RKIND), dimension(:,:), pointer :: les_wtheta_m_flux, les_wtheta_flux, les_wqv_flux
      real(kind=RKIND), dimension(:), pointer :: hfx, qfx, ust, cpm
      real(kind=RKIND), dimension(:,:), pointer :: tend_u, tend_w, tend_theta
      real(kind=RKIND), dimension(:,:), pointer :: tend_u_euler, tend_w_euler, tend_theta_euler
      real(kind=RKIND), dimension(:,:,:), pointer :: tend_scalars
      real(kind=RKIND), allocatable :: vertical_flux(:), cross_prod(:,:)
      integer :: i, iCell, iEdge, iScalar, k
      integer :: cell1, cell2, vertex1, vertex2
      real(kind=RKIND) :: edge_sign, rho_int, coeff, coeff_up, coeff_down, grad_up, grad_down
      real(kind=RKIND) :: zc_low, zc_here, zc_high, dz_up, dz_down
      real(kind=RKIND) :: diffusion, scalar_diffusion, schmidt_inv
      real(kind=RKIND) :: d_diag, d_off_diag, du_dz, dv_dz, dw_dz
      real(kind=RKIND) :: km_h_safe, km_v_safe, kh_v_safe, mix_safe, bv_safe
      real(kind=RKIND) :: shear_prod, buoy_prod, dissipation, source
      real(kind=RKIND) :: shear_prod_h_value, shear_prod_v_internal_value
      real(kind=RKIND) :: surface_prod_stress_shear
      real(kind=RKIND) :: source_raw, budget_sum
      real(kind=RKIND) :: dz, delta, tke, ce1, ce2, ce
      real(kind=RKIND) :: speed, surface_prod, theta_dry
      real(kind=RKIND) :: edge_ust, rho_surface, cpm_surface
      real(kind=RKIND) :: rho_phys_low, rho_phys_high
      real(kind=RKIND) :: heat_flux_kinematic, moisture_flux
      real(kind=RKIND) :: rho_tend, velocity_scale_cell
      logical :: use_wrf_tke_diffusion, apply_scalar_tendencies
      logical, pointer :: config_apply_lbcs
      integer, pointer :: affine_bdy_mask(:)

      if (rk_step /= 1) return
      if (nVertLevels < 2) return
      if (dt <= 0.0_RKIND) return

      ! Scalar transport is outside the dry-dynamics substep loop and consumes
      ! scalars_tend once with the full outer timestep.  Dynamic tendencies are
      ! integrated inside that loop with dt/dynamics_split, so only persistent
      ! scalar-tendency writes are restricted to its first substep.
      apply_scalar_tendencies = (dynamics_substep == 1)

      nullify(config_les_scheme)
      call iamas_pool_get_config(configs, 'config_les_scheme', config_les_scheme)
      if (.not. associated(config_les_scheme)) return
      if (trim(config_les_scheme) /= 'les_deardorff') return

      call iamas_pool_get_config(configs, 'config_les_ck', config_les_ck)
      call iamas_pool_get_config(configs, 'config_les_ce0', config_les_ce0)
      call iamas_pool_get_config(configs, 'config_les_scalar_schmidt', config_les_scalar_schmidt)
      call iamas_pool_get_config(configs, 'config_les_tke_diffusion_mode', config_les_tke_diffusion_mode)
      call iamas_pool_get_config(configs, 'config_les_tke_min', config_les_tke_min)
      call iamas_pool_get_config(configs, 'config_les_surface_mode', config_les_surface_mode)
      call iamas_pool_get_config(configs, 'config_les_heat_flux', config_les_heat_flux)
      call iamas_pool_get_config(configs, 'config_les_drag_coefficient', config_les_drag_coefficient)
      call iamas_pool_get_config(configs, 'config_les_output_budget', config_les_output_budget)
      call iamas_pool_get_dimension(state, 'num_scalars', num_scalars)
      call iamas_pool_get_dimension(state, 'index_qv', index_qv)
      call iamas_pool_get_dimension(state, 'index_tke_sgs', index_tke_sgs)

      nullify(zz, rdzw)
      call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call iamas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call iamas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call iamas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call iamas_pool_get_array(mesh, 'invDcEdge', invDcEdge)
      call iamas_pool_get_array(mesh, 'invDvEdge', invDvEdge)
      call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call iamas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call iamas_pool_get_array(mesh, 'areaCell', areaCell)
      call iamas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(mesh, 'zz', zz)
      call iamas_pool_get_array(mesh, 'rdzw', rdzw)
      call iamas_pool_get_config(configs,'config_apply_lbcs',config_apply_lbcs)
      call iamas_pool_get_array(mesh,'bdyMaskCell',affine_bdy_mask)
      if (.not.associated(config_apply_lbcs) .or. .not.associated(affine_bdy_mask)) &
         call iamas_les_affine_abort('C2 LES boundary contract absent')
      call iamas_les_prepare_affine(mesh,cellSolveStart,cellSolveEnd,config_apply_lbcs)
      call iamas_pool_get_array(mesh,'les_affine_d1u',affine_d1u)
      call iamas_pool_get_array(mesh,'les_affine_d1v',affine_d1v)
      call iamas_pool_get_array(mesh,'les_affine_d2u',affine_d2u)
      call iamas_pool_get_array(mesh,'les_affine_d2v',affine_d2v)

      if (.not. associated(zz) .or. .not. associated(rdzw)) then
         error stop 'E6 LES requires mesh fields zz and rdzw'
      end if
      if (size(zz,1) < nVertLevels .or. size(zz,2) < cellSolveEnd .or. &
          size(rdzw) < nVertLevels) then
         error stop 'E6 LES mesh fields zz or rdzw have invalid bounds'
      end if

      call iamas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call iamas_pool_get_array(state, 'u', u, 2)
      call iamas_pool_get_array(state, 'w', w, 2)
      call iamas_pool_get_array(state, 'theta_m', theta_m, 2)
      call iamas_pool_get_array(state, 'scalars', scalars, 2)

      call iamas_pool_get_array(diag, 'rho_edge', rho_edge)
      call iamas_pool_get_array(diag, 'divergence', divergence)
      call iamas_pool_get_array(diag, 'vorticity', vorticity)
      call iamas_pool_get_array(diag, 'v', v)
      call iamas_pool_get_array(diag, 'uReconstructZonal', ur_cell)
      call iamas_pool_get_array(diag, 'uReconstructMeridional', vr_cell)
      call iamas_pool_get_array(diag, 'les_mix_length', les_mix_length)
      call iamas_pool_get_array(diag, 'les_bv_freq2', les_bv_freq2)
      call iamas_pool_get_array(diag, 'les_km_h', les_km_h)
      call iamas_pool_get_array(diag, 'les_km_v', les_km_v)
      call iamas_pool_get_array(diag, 'les_kh_h', les_kh_h)
      call iamas_pool_get_array(diag, 'les_kh_v', les_kh_v)
      call iamas_pool_get_array(diag, 'les_shear_prod', les_shear_prod)
      call iamas_pool_get_array(diag, 'les_shear_prod_h', les_shear_prod_h)
      call iamas_pool_get_array(diag, 'les_shear_prod_v_internal', les_shear_prod_v_internal)
      call iamas_pool_get_array(diag, 'les_surface_prod_current', les_surface_prod_current)
      call iamas_pool_get_array(diag, 'les_surface_prod_stress_shear', les_surface_prod_stress_shear)
      call iamas_pool_get_array(diag, 'les_surface_resolved_ke_sink_edge', &
                                les_surface_resolved_ke_sink_edge)
      call iamas_pool_get_array(diag, 'les_buoy_prod', les_buoy_prod)
      call iamas_pool_get_array(diag, 'les_dissipation', les_dissipation)
      call iamas_pool_get_array(diag, 'les_tke_diffusion', les_tke_diffusion)
      call iamas_pool_get_array(diag, 'les_tke_diffusion_h', les_tke_diffusion_h)
      call iamas_pool_get_array(diag, 'les_tke_diffusion_v', les_tke_diffusion_v)
      call iamas_pool_get_array(diag, 'les_tke_h_flux', les_tke_h_flux)
      call iamas_pool_get_array(diag, 'les_tke_v_flux', les_tke_v_flux)
      call iamas_pool_get_array(diag, 'les_tke_tend_total', les_tke_tend_total)
      call iamas_pool_get_array(diag, 'les_tke_source_pre_limiter', les_tke_source_pre_limiter)
      call iamas_pool_get_array(diag, 'les_tke_source_post_limiter', les_tke_source_post_limiter)
      call iamas_pool_get_array(diag, 'les_tke_limiter_residual', les_tke_limiter_residual)
      call iamas_pool_get_array(diag, 'les_tke_budget_closure_residual', les_tke_budget_closure_residual)
      call iamas_pool_get_array(diag, 'les_tend_un', les_tend_un)
      call iamas_pool_get_array(diag, 'les_tend_w', les_tend_w)
      call iamas_pool_get_array(diag, 'les_tend_theta_m', les_tend_theta_m)
      call iamas_pool_get_array(diag, 'les_tend_theta', les_tend_theta)
      call iamas_pool_get_array(diag, 'les_tend_qv', les_tend_qv)
      call iamas_pool_get_array(diag, 'les_unw_flux', les_unw_flux)
      call iamas_pool_get_array(diag, 'les_uw_flux', les_uw_flux)
      call iamas_pool_get_array(diag, 'les_vw_flux', les_vw_flux)
      call iamas_pool_get_array(diag, 'les_wtheta_m_flux', les_wtheta_m_flux)
      call iamas_pool_get_array(diag, 'les_wtheta_flux', les_wtheta_flux)
      call iamas_pool_get_array(diag, 'les_wqv_flux', les_wqv_flux)

      nullify(hfx, qfx, ust, cpm)
      if (trim(config_les_surface_mode) == 'physics_flux') then
         call iamas_pool_get_array(diag_physics, 'hfx', hfx)
         call iamas_pool_get_array(diag_physics, 'qfx', qfx)
         call iamas_pool_get_array(diag_physics, 'ust', ust)
         call iamas_pool_get_array(diag_physics, 'cpm', cpm)
         if (.not. associated(hfx) .or. .not. associated(qfx) .or. &
             .not. associated(ust)) return
      end if

      call iamas_pool_get_array(tend, 'u', tend_u)
      call iamas_pool_get_array(tend, 'w', tend_w)
      call iamas_pool_get_array(tend, 'theta_m', tend_theta)
      call iamas_pool_get_array(tend, 'u_euler', tend_u_euler)
      call iamas_pool_get_array(tend, 'w_euler', tend_w_euler)
      call iamas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call iamas_pool_get_array(tend, 'scalars_tend', tend_scalars)

      if (.not. associated(num_scalars) .or. .not. associated(index_tke_sgs)) return
      if (.not. associated(index_qv) .or. index_qv <= 0) return
      if (.not. associated(tend_scalars) .or. index_tke_sgs <= 0) return
      if (.not. associated(config_les_tke_diffusion_mode)) return
      if (.not. associated(les_tke_diffusion) .or. .not. associated(les_tke_diffusion_h)) return
      if (.not. associated(les_tke_diffusion_v) .or. .not. associated(les_tke_h_flux)) return
      if (.not. associated(les_tke_v_flux)) return
      if (.not. associated(les_tend_un) .or. .not. associated(les_tend_w)) return
      if (.not. associated(les_tend_theta_m) .or. .not. associated(les_tend_theta)) return
      if (.not. associated(les_tend_qv)) return
      if (.not. associated(les_unw_flux) .or. .not. associated(les_uw_flux)) return
      if (.not. associated(les_vw_flux) .or. .not. associated(les_wtheta_m_flux)) return
      if (.not. associated(les_wtheta_flux)) return
      if (.not. associated(les_wqv_flux)) return
      if (.not. associated(les_shear_prod_h) .or. &
          .not. associated(les_shear_prod_v_internal)) return
      if (.not. associated(les_surface_prod_current) .or. &
          .not. associated(les_surface_prod_stress_shear)) return
      if (.not. associated(les_surface_resolved_ke_sink_edge)) return
      schmidt_inv = 1.0_RKIND / max(config_les_scalar_schmidt, 1.0e-6_RKIND)
      select case (trim(config_les_tke_diffusion_mode))
      case ('legacy_km_schmidt')
         use_wrf_tke_diffusion = .false.
      case ('wrf_kh_x2')
         use_wrf_tke_diffusion = .true.
      case default
         error stop 'Unsupported config_les_tke_diffusion_mode'
      end select

      do iCell = cellSolveStart, cellSolveEnd
         les_shear_prod(:,iCell) = 0.0_RKIND
         les_shear_prod_h(:,iCell) = 0.0_RKIND
         les_shear_prod_v_internal(:,iCell) = 0.0_RKIND
         les_surface_prod_current(:,iCell) = 0.0_RKIND
         les_surface_prod_stress_shear(:,iCell) = 0.0_RKIND
         les_buoy_prod(:,iCell) = 0.0_RKIND
         les_dissipation(:,iCell) = 0.0_RKIND
         les_tke_diffusion(:,iCell) = 0.0_RKIND
         les_tke_diffusion_h(:,iCell) = 0.0_RKIND
         les_tke_diffusion_v(:,iCell) = 0.0_RKIND
         les_tke_v_flux(:,iCell) = 0.0_RKIND
         if (associated(les_tke_tend_total)) les_tke_tend_total(:,iCell) = 0.0_RKIND
         if (associated(les_tke_source_pre_limiter)) les_tke_source_pre_limiter(:,iCell) = 0.0_RKIND
         if (associated(les_tke_source_post_limiter)) les_tke_source_post_limiter(:,iCell) = 0.0_RKIND
         if (associated(les_tke_limiter_residual)) les_tke_limiter_residual(:,iCell) = 0.0_RKIND
         if (associated(les_tke_budget_closure_residual)) then
            les_tke_budget_closure_residual(:,iCell) = 0.0_RKIND
         end if
         les_tend_w(:,iCell) = 0.0_RKIND
         les_tend_theta_m(:,iCell) = 0.0_RKIND
         les_tend_theta(:,iCell) = 0.0_RKIND
         les_tend_qv(:,iCell) = 0.0_RKIND
         les_uw_flux(:,iCell) = 0.0_RKIND
         les_vw_flux(:,iCell) = 0.0_RKIND
         les_wtheta_m_flux(:,iCell) = 0.0_RKIND
         les_wtheta_flux(:,iCell) = 0.0_RKIND
         les_wqv_flux(:,iCell) = 0.0_RKIND
      end do
      do iEdge = edgeSolveStart, edgeSolveEnd
         les_tend_un(:,iEdge) = 0.0_RKIND
         les_unw_flux(:,iEdge) = 0.0_RKIND
         les_tke_h_flux(:,iEdge) = 0.0_RKIND
         les_surface_resolved_ke_sink_edge(:,iEdge) = 0.0_RKIND
      end do

      ! C3: full edge-normal xz/yz stress, physical-density flux and its
      ! coordinate-density divergence, with the same TKE energy transfer.
      allocate(cross_prod(nVertLevels,cellSolveEnd-cellSolveStart+1))
      call iamas_les_cross_apply(mesh,configs,nVertLevels,cellSolveStart,cellSolveEnd, &
           edgeSolveStart,edgeSolveEnd,rho_zz,u,w,zz,rdzw,les_km_v, &
           tend_u,tend_u_euler,tend_w,tend_w_euler,les_tend_un,les_tend_w, &
           les_unw_flux,les_uw_flux,les_vw_flux,cross_prod)

      ! Horizontal momentum diffusion on edges. This follows the native iAMAS
      ! vector-Laplacian discretization, but uses the LES eddy viscosity.
      do iEdge = edgeSolveStart, edgeSolveEnd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         if (cell1 <= 0 .or. cell2 <= 0 .or. vertex1 <= 0 .or. vertex2 <= 0) cycle
         do k = 1, nVertLevels
            diffusion = 2.0_RKIND * (divergence(k,cell2) - divergence(k,cell1)) * invDcEdge(iEdge) - &
                        (vorticity(k,vertex2) - vorticity(k,vertex1)) * &
                        min(invDvEdge(iEdge), 4.0_RKIND * invDcEdge(iEdge))
            coeff = iamas_les_safe_nonnegative(0.5_RKIND * (les_km_h(k,cell1) + les_km_h(k,cell2)))
            diffusion = rho_edge(k,iEdge) * coeff * meshScalingDel2(iEdge) * diffusion
            tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + diffusion
            tend_u(k,iEdge) = tend_u(k,iEdge) + diffusion
            les_tend_un(k,iEdge) = les_tend_un(k,iEdge) + &
                                    diffusion / max(rho_edge(k,iEdge), 1.0e-12_RKIND)
         end do
      end do

      ! Cell-centered horizontal diffusion for vertical momentum, heat, and
      ! all active atmospheric scalars.
      do iCell = cellSolveStart, cellSolveEnd
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= 0 .or. cell2 <= 0) cycle
            edge_sign = invAreaCell(iCell) * edgesOnCell_sign(i,iCell) * &
                        dvEdge(iEdge) * invDcEdge(iEdge) * meshScalingDel2(iEdge)

            do k = 1, nVertLevels
               coeff = iamas_les_safe_nonnegative(0.5_RKIND * (les_kh_h(k,cell1) + les_kh_h(k,cell2)))
               diffusion = edge_sign * rho_edge(k,iEdge) * coeff * &
                           (theta_m(k,cell2) - theta_m(k,cell1))
               tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + diffusion
               tend_theta(k,iCell) = tend_theta(k,iCell) + diffusion
               les_tend_theta_m(k,iCell) = les_tend_theta_m(k,iCell) + &
                                            diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)

               coeff = iamas_les_safe_nonnegative(0.5_RKIND * (les_km_h(k,cell1) + les_km_h(k,cell2))) * schmidt_inv
               do iScalar = 1, num_scalars
                  if (use_wrf_tke_diffusion .and. iScalar == index_tke_sgs) cycle
                  coeff_up = coeff
                  if (iScalar == index_qv) then
                     coeff_up = iamas_les_safe_nonnegative(0.5_RKIND * &
                                (les_kh_h(k,cell1) + les_kh_h(k,cell2)))
                  end if
                  scalar_diffusion = edge_sign * rho_edge(k,iEdge) * coeff_up * &
                                     (scalars(iScalar,k,cell2) - scalars(iScalar,k,cell1))
                  if (apply_scalar_tendencies) then
                     tend_scalars(iScalar,k,iCell) = tend_scalars(iScalar,k,iCell) + scalar_diffusion
                  end if
                  if (iScalar == index_qv) then
                     les_tend_qv(k,iCell) = les_tend_qv(k,iCell) + &
                                             scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
                  end if
                  if (iScalar == index_tke_sgs) then
                     diffusion = scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
                     les_tke_diffusion_h(k,iCell) = les_tke_diffusion_h(k,iCell) + diffusion
                     les_tke_diffusion(k,iCell) = les_tke_diffusion(k,iCell) + diffusion
                     les_tke_h_flux(k,iEdge) = -coeff * meshScalingDel2(iEdge) * invDcEdge(iEdge) * &
                                               (scalars(index_tke_sgs,k,cell2) - &
                                                scalars(index_tke_sgs,k,cell1))
                  end if
               end do

               if (use_wrf_tke_diffusion) then
                  coeff = 2.0_RKIND * iamas_les_safe_nonnegative(0.5_RKIND * &
                          (les_kh_h(k,cell1) + les_kh_h(k,cell2)))
                  scalar_diffusion = edge_sign * rho_edge(k,iEdge) * coeff * &
                                     (scalars(index_tke_sgs,k,cell2) - &
                                      scalars(index_tke_sgs,k,cell1))
                  if (apply_scalar_tendencies) then
                     tend_scalars(index_tke_sgs,k,iCell) = &
                          tend_scalars(index_tke_sgs,k,iCell) + scalar_diffusion
                  end if
                  diffusion = scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
                  les_tke_diffusion_h(k,iCell) = les_tke_diffusion_h(k,iCell) + diffusion
                  les_tke_diffusion(k,iCell) = les_tke_diffusion(k,iCell) + diffusion
                  les_tke_h_flux(k,iEdge) = -coeff * meshScalingDel2(iEdge) * invDcEdge(iEdge) * &
                                            (scalars(index_tke_sgs,k,cell2) - &
                                             scalars(index_tke_sgs,k,cell1))
               end if
            end do

            ! Vertical-w horizontal mixing is part of the common C3 cross stress.
         end do
      end do

      allocate(vertical_flux(nVertLevels+1))

      ! C3 vertical u diffusion is included in the common cross stress above.

      ! E6 candidate: vertical heat/scalar fluxes use physical density
      ! rho=zz*rho_zz, while their divergence remains the coordinate-density
      ! tendency rdzw*(G_up-G_down).
      ! Vertical diffusion of heat and scalars at cell centers.
      do iCell = cellSolveStart, cellSolveEnd
         vertical_flux = 0.0_RKIND
         do k = 2, nVertLevels
            zc_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
            zc_here = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
            coeff = iamas_les_safe_nonnegative(0.5_RKIND * (les_kh_v(k-1,iCell) + les_kh_v(k,iCell)))
            rho_phys_low = iamas_les_safe_nonnegative(zz(k-1,iCell) * rho_zz(k-1,iCell))
            rho_phys_high = iamas_les_safe_nonnegative(zz(k,iCell) * rho_zz(k,iCell))
            rho_int = 0.5_RKIND * (rho_phys_low + rho_phys_high)
            vertical_flux(k) = rho_int * coeff * (theta_m(k,iCell) - theta_m(k-1,iCell)) / &
                               max(abs(zc_here-zc_low), 1.0_RKIND)
            les_wtheta_m_flux(k,iCell) = -vertical_flux(k) / max(rho_int, 1.0e-12_RKIND)
         end do
         do k = 1, nVertLevels
            diffusion = rdzw(k) * (vertical_flux(k+1) - vertical_flux(k))
            tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + diffusion
            tend_theta(k,iCell) = tend_theta(k,iCell) + diffusion
            les_tend_theta_m(k,iCell) = les_tend_theta_m(k,iCell) + &
                                         diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
         end do

         do iScalar = 1, num_scalars
            if (use_wrf_tke_diffusion .and. iScalar == index_tke_sgs) cycle
            vertical_flux = 0.0_RKIND
            do k = 2, nVertLevels
               zc_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               zc_here = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
               if (iScalar == index_qv) then
                  coeff = iamas_les_safe_nonnegative(0.5_RKIND * &
                          (les_kh_v(k-1,iCell) + les_kh_v(k,iCell)))
               else
                  coeff = iamas_les_safe_nonnegative(0.5_RKIND * &
                          (les_km_v(k-1,iCell) + les_km_v(k,iCell))) * schmidt_inv
               end if
               rho_phys_low = iamas_les_safe_nonnegative(zz(k-1,iCell) * rho_zz(k-1,iCell))
               rho_phys_high = iamas_les_safe_nonnegative(zz(k,iCell) * rho_zz(k,iCell))
               rho_int = 0.5_RKIND * (rho_phys_low + rho_phys_high)
               vertical_flux(k) = rho_int * coeff * &
                                  (scalars(iScalar,k,iCell) - scalars(iScalar,k-1,iCell)) / &
                                  max(abs(zc_here-zc_low), 1.0_RKIND)
               if (iScalar == index_qv) then
                  les_wqv_flux(k,iCell) = -vertical_flux(k) / max(rho_int, 1.0e-12_RKIND)
               end if
               if (iScalar == index_tke_sgs) then
                  les_tke_v_flux(k,iCell) = -vertical_flux(k) / max(rho_int, 1.0e-12_RKIND)
               end if
            end do
            do k = 1, nVertLevels
               scalar_diffusion = rdzw(k) * (vertical_flux(k+1) - vertical_flux(k))
               if (apply_scalar_tendencies) then
                  tend_scalars(iScalar,k,iCell) = tend_scalars(iScalar,k,iCell) + scalar_diffusion
               end if
               if (iScalar == index_qv) then
                  les_tend_qv(k,iCell) = les_tend_qv(k,iCell) + &
                                          scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
               end if
               if (iScalar == index_tke_sgs) then
                  diffusion = scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
                  les_tke_diffusion_v(k,iCell) = les_tke_diffusion_v(k,iCell) + diffusion
                  les_tke_diffusion(k,iCell) = les_tke_diffusion(k,iCell) + diffusion
               end if
            end do
         end do

         if (use_wrf_tke_diffusion) then
            vertical_flux = 0.0_RKIND
            do k = 2, nVertLevels
               zc_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               zc_here = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
               coeff = 2.0_RKIND * iamas_les_safe_nonnegative(0.5_RKIND * &
                       (les_kh_v(k-1,iCell) + les_kh_v(k,iCell)))
               rho_phys_low = iamas_les_safe_nonnegative(zz(k-1,iCell) * rho_zz(k-1,iCell))
               rho_phys_high = iamas_les_safe_nonnegative(zz(k,iCell) * rho_zz(k,iCell))
               rho_int = 0.5_RKIND * (rho_phys_low + rho_phys_high)
               vertical_flux(k) = rho_int * coeff * &
                                  (scalars(index_tke_sgs,k,iCell) - &
                                   scalars(index_tke_sgs,k-1,iCell)) / &
                                  max(abs(zc_here-zc_low), 1.0_RKIND)
               les_tke_v_flux(k,iCell) = -vertical_flux(k) / max(rho_int, 1.0e-12_RKIND)
            end do
            do k = 1, nVertLevels
               scalar_diffusion = rdzw(k) * (vertical_flux(k+1) - vertical_flux(k))
               if (apply_scalar_tendencies) then
                  tend_scalars(index_tke_sgs,k,iCell) = &
                       tend_scalars(index_tke_sgs,k,iCell) + scalar_diffusion
               end if
               diffusion = scalar_diffusion / max(rho_zz(k,iCell), 1.0e-12_RKIND)
               les_tke_diffusion_v(k,iCell) = les_tke_diffusion_v(k,iCell) + diffusion
               les_tke_diffusion(k,iCell) = les_tke_diffusion(k,iCell) + diffusion
            end do
         end if
      end do

      ! C3 normal stress tau_zz=2*rho*Kv*dw/dz. As in MPAS, physical
      ! density enters the flux and reference coordinate thickness its divergence.
      do iCell = cellSolveStart, cellSolveEnd
         if (.not.iamas_les_affine_cell_active(config_apply_lbcs,affine_bdy_mask(iCell))) cycle
         do k = 1, nVertLevels
            dz = zgrid(k+1,iCell)-zgrid(k,iCell)
            if (dz<=0.0_RKIND) call iamas_les_affine_abort('C3 invalid physical layer thickness')
            vertical_flux(k) = 2.0_RKIND * rho_zz(k,iCell)*zz(k,iCell) * les_km_v(k,iCell) * &
                               (w(k+1,iCell)-w(k,iCell))/dz
         end do
         do k = 2, nVertLevels
            diffusion = (vertical_flux(k)-vertical_flux(k-1)) * &
                        2.0_RKIND/(1.0_RKIND/rdzw(k-1)+1.0_RKIND/rdzw(k))
            rho_int = 0.5_RKIND*(rho_zz(k,iCell)+rho_zz(k-1,iCell))
            tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + diffusion
            tend_w(k,iCell) = tend_w(k,iCell) + diffusion
            les_tend_w(k,iCell) = les_tend_w(k,iCell) + diffusion/max(rho_int,1.0e-12_RKIND)
         end do
      end do

      ! Lower-boundary fluxes. Prescribed mode is retained for idealized LES.
      ! In physics_flux mode, Noah and the surface-layer scheme diagnose HFX,
      ! QFX, and UST while the PBL driver remains off; SGS applies the fluxes
      ! to the lowest atmospheric layer.
      if (trim(config_les_surface_mode) == 'prescribed') then
         do iEdge = edgeSolveStart, edgeSolveEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= 0 .or. cell2 <= 0) cycle
            dz = max(0.5_RKIND * (abs(zgrid(2,cell1)-zgrid(1,cell1)) + &
                                  abs(zgrid(2,cell2)-zgrid(1,cell2))), 1.0_RKIND)
            speed = 0.5_RKIND * (sqrt(ur_cell(1,cell1)**2 + vr_cell(1,cell1)**2) + &
                                  sqrt(ur_cell(1,cell2)**2 + vr_cell(1,cell2)**2))
            diffusion = -rho_edge(1,iEdge) * max(config_les_drag_coefficient, 0.0_RKIND) * &
                        speed * u(1,iEdge) / dz
            tend_u_euler(1,iEdge) = tend_u_euler(1,iEdge) + diffusion
            tend_u(1,iEdge) = tend_u(1,iEdge) + diffusion
            les_tend_un(1,iEdge) = les_tend_un(1,iEdge) + &
                                    diffusion / max(rho_edge(1,iEdge), 1.0e-12_RKIND)
            les_unw_flux(1,iEdge) = -max(config_les_drag_coefficient, 0.0_RKIND) * speed * u(1,iEdge)
            les_surface_resolved_ke_sink_edge(1,iEdge) = iamas_les_safe_finite( &
                 -u(1,iEdge) * diffusion / max(rho_edge(1,iEdge), 1.0e-12_RKIND))
         end do

         do iCell = cellSolveStart, cellSolveEnd
            dz = max(abs(zgrid(2,iCell)-zgrid(1,iCell)), 1.0_RKIND)
            diffusion = rho_zz(1,iCell) * config_les_heat_flux / dz
            tend_theta_euler(1,iCell) = tend_theta_euler(1,iCell) + diffusion
            tend_theta(1,iCell) = tend_theta(1,iCell) + diffusion
            les_tend_theta_m(1,iCell) = les_tend_theta_m(1,iCell) + &
                                         diffusion / max(rho_zz(1,iCell), 1.0e-12_RKIND)
            les_wtheta_m_flux(1,iCell) = config_les_heat_flux
         end do
      else if (trim(config_les_surface_mode) == 'physics_flux') then
         do iEdge = edgeSolveStart, edgeSolveEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= 0 .or. cell2 <= 0) cycle
            dz = max(0.5_RKIND * (abs(zgrid(2,cell1)-zgrid(1,cell1)) + &
                                  abs(zgrid(2,cell2)-zgrid(1,cell2))), 1.0_RKIND)
            speed = 0.5_RKIND * (sqrt(ur_cell(1,cell1)**2 + vr_cell(1,cell1)**2) + &
                                  sqrt(ur_cell(1,cell2)**2 + vr_cell(1,cell2)**2))
            edge_ust = 0.5_RKIND * (iamas_les_safe_nonnegative(ust(cell1)) + &
                                    iamas_les_safe_nonnegative(ust(cell2)))
            diffusion = -rho_edge(1,iEdge) * edge_ust**2 * u(1,iEdge) / &
                        (max(speed, 1.0e-6_RKIND) * dz)
            tend_u_euler(1,iEdge) = tend_u_euler(1,iEdge) + diffusion
            tend_u(1,iEdge) = tend_u(1,iEdge) + diffusion
            les_tend_un(1,iEdge) = les_tend_un(1,iEdge) + &
                                    diffusion / max(rho_edge(1,iEdge), 1.0e-12_RKIND)
            les_unw_flux(1,iEdge) = -edge_ust**2 * u(1,iEdge) / max(speed, 1.0e-6_RKIND)
            les_surface_resolved_ke_sink_edge(1,iEdge) = iamas_les_safe_finite( &
                 -u(1,iEdge) * diffusion / max(rho_edge(1,iEdge), 1.0e-12_RKIND))
         end do

         ! E6 candidate: HFX/CPM and QFX are physical surface flux targets.
         ! Convert them once to a lowest-layer coordinate-density source with
         ! rdzw(1); saved kinematic fluxes divide by physical density.
         do iCell = cellSolveStart, cellSolveEnd
            rho_surface = max(iamas_les_safe_nonnegative(zz(1,iCell) * &
                              rho_zz(1,iCell)), 1.0e-6_RKIND)
            if (associated(cpm)) then
               cpm_surface = max(iamas_les_safe_nonnegative(cpm(iCell)), 100.0_RKIND)
            else if (associated(index_qv) .and. index_qv > 0) then
               cpm_surface = cp * (1.0_RKIND + 0.8_RKIND * &
                                    max(scalars(index_qv,1,iCell), 0.0_RKIND))
            else
               cpm_surface = cp
            end if
            heat_flux_kinematic = iamas_les_safe_finite(hfx(iCell)) / &
                                  (rho_surface * cpm_surface)
            moisture_flux = iamas_les_safe_finite(qfx(iCell))

            diffusion = rho_surface * heat_flux_kinematic
            if (associated(index_qv) .and. index_qv > 0) then
               theta_dry = theta_m(1,iCell) / &
                           (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND))
               diffusion = diffusion * &
                           (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND)) + &
                           rvord * theta_dry * moisture_flux
               scalar_diffusion = moisture_flux * rdzw(1)
               if (apply_scalar_tendencies) then
                  tend_scalars(index_qv,1,iCell) = tend_scalars(index_qv,1,iCell) + &
                                                   scalar_diffusion
               end if
               les_tend_qv(1,iCell) = les_tend_qv(1,iCell) + &
                                       scalar_diffusion / max(rho_zz(1,iCell), 1.0e-12_RKIND)
            end if
            diffusion = diffusion * rdzw(1)
            tend_theta_euler(1,iCell) = tend_theta_euler(1,iCell) + diffusion
            tend_theta(1,iCell) = tend_theta(1,iCell) + diffusion
            les_tend_theta_m(1,iCell) = les_tend_theta_m(1,iCell) + &
                                         diffusion / max(rho_zz(1,iCell), 1.0e-12_RKIND)
            les_wqv_flux(1,iCell) = moisture_flux / rho_surface
            les_wtheta_m_flux(1,iCell) = heat_flux_kinematic * &
                                         (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND)) + &
                                         rvord * theta_dry * les_wqv_flux(1,iCell)
         end do
      end if

      ! C3 interior fluxes were reconstructed from the common edge stress.
      ! Surface exchange and dry-theta conversions retain their existing paths.
      do iCell = cellSolveStart, cellSolveEnd
         if (trim(config_les_surface_mode) == 'physics_flux') then
            velocity_scale_cell = max(sqrt(ur_cell(1,iCell)**2 + vr_cell(1,iCell)**2), 1.0e-6_RKIND)
            les_uw_flux(1,iCell) = -iamas_les_safe_nonnegative(ust(iCell))**2 * &
                                    ur_cell(1,iCell) / velocity_scale_cell
            les_vw_flux(1,iCell) = -iamas_les_safe_nonnegative(ust(iCell))**2 * &
                                    vr_cell(1,iCell) / velocity_scale_cell
         else if (trim(config_les_surface_mode) == 'prescribed') then
            velocity_scale_cell = sqrt(ur_cell(1,iCell)**2 + vr_cell(1,iCell)**2)
            les_uw_flux(1,iCell) = -max(config_les_drag_coefficient, 0.0_RKIND) * &
                                    velocity_scale_cell * ur_cell(1,iCell)
            les_vw_flux(1,iCell) = -max(config_les_drag_coefficient, 0.0_RKIND) * &
                                    velocity_scale_cell * vr_cell(1,iCell)
         end if

         do k = 1, nVertLevels
            theta_dry = theta_m(k,iCell) / &
                        (1.0_RKIND + rvord * max(scalars(index_qv,k,iCell), 0.0_RKIND))
            les_tend_theta(k,iCell) = (les_tend_theta_m(k,iCell) - &
                                      rvord * theta_dry * les_tend_qv(k,iCell)) / &
                                      (1.0_RKIND + rvord * max(scalars(index_qv,k,iCell), 0.0_RKIND))
         end do

         theta_dry = theta_m(1,iCell) / &
                     (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND))
         les_wtheta_flux(1,iCell) = (les_wtheta_m_flux(1,iCell) - &
                                     rvord * theta_dry * les_wqv_flux(1,iCell)) / &
                                     (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND))
         do k = 2, nVertLevels
            theta_dry = 0.5_RKIND * ( &
               theta_m(k-1,iCell) / (1.0_RKIND + rvord * max(scalars(index_qv,k-1,iCell), 0.0_RKIND)) + &
               theta_m(k,iCell) / (1.0_RKIND + rvord * max(scalars(index_qv,k,iCell), 0.0_RKIND)))
            les_wtheta_flux(k,iCell) = (les_wtheta_m_flux(k,iCell) - &
                                        rvord * theta_dry * les_wqv_flux(k,iCell)) / &
                                        (1.0_RKIND + rvord * 0.5_RKIND * &
                                         (max(scalars(index_qv,k-1,iCell), 0.0_RKIND) + &
                                          max(scalars(index_qv,k,iCell), 0.0_RKIND)))
         end do
      end do

      ! E6MSA2C: native-affine horizontal deformation, with geometry and
      ! tangential reconstruction treated as one measurement operator.
      ! C3 adds the native divergence normal term, paired with 2*grad(div).
      ! Buoyancy and dissipation follow the WRF Deardorff closure.
      ce1 = 0.19_RKIND * config_les_ck / 0.10_RKIND
      ce2 = max(0.0_RKIND, config_les_ce0 - ce1)
      do iCell = cellSolveStart, cellSolveEnd
         ! C2: match the scalar integrator's specified-zone exclusion. These
         ! diagnostic source slots remain inactive, not physical zero evidence.
         if (.not.iamas_les_affine_cell_active(config_apply_lbcs,affine_bdy_mask(iCell))) cycle
         do k = 1, nVertLevels
            d_diag = 0.0_RKIND
            d_off_diag = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               d_diag = d_diag + affine_d1u(i,iCell)*u(k,iEdge) + affine_d1v(i,iCell)*v(k,iEdge)
               d_off_diag = d_off_diag + affine_d2u(i,iCell)*u(k,iEdge) + affine_d2v(i,iCell)*v(k,iEdge)
            end do

            if (k == 1) then
               zc_here = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
               zc_high = 0.5_RKIND * (zgrid(k+1,iCell) + zgrid(k+2,iCell))
               du_dz = (ur_cell(k+1,iCell) - ur_cell(k,iCell)) / max(abs(zc_high-zc_here), 1.0_RKIND)
               dv_dz = (vr_cell(k+1,iCell) - vr_cell(k,iCell)) / max(abs(zc_high-zc_here), 1.0_RKIND)
            else if (k == nVertLevels) then
               zc_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               zc_here = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
               du_dz = (ur_cell(k,iCell) - ur_cell(k-1,iCell)) / max(abs(zc_here-zc_low), 1.0_RKIND)
               dv_dz = (vr_cell(k,iCell) - vr_cell(k-1,iCell)) / max(abs(zc_here-zc_low), 1.0_RKIND)
            else
               zc_low = 0.5_RKIND * (zgrid(k-1,iCell) + zgrid(k,iCell))
               zc_high = 0.5_RKIND * (zgrid(k+1,iCell) + zgrid(k+2,iCell))
               du_dz = (ur_cell(k+1,iCell) - ur_cell(k-1,iCell)) / max(abs(zc_high-zc_low), 1.0_RKIND)
               dv_dz = (vr_cell(k+1,iCell) - vr_cell(k-1,iCell)) / max(abs(zc_high-zc_low), 1.0_RKIND)
            end if

            dz = max(abs(zgrid(k+1,iCell) - zgrid(k,iCell)), 1.0_RKIND)
            d_diag = iamas_les_safe_finite(d_diag)
            d_off_diag = iamas_les_safe_finite(d_off_diag)
            du_dz = iamas_les_safe_finite(du_dz)
            dv_dz = iamas_les_safe_finite(dv_dz)
            dw_dz = iamas_les_safe_finite((w(k+1,iCell) - w(k,iCell)) / dz)
            km_h_safe = iamas_les_safe_nonnegative(les_km_h(k,iCell))
            km_v_safe = iamas_les_safe_nonnegative(les_km_v(k,iCell))
            kh_v_safe = iamas_les_safe_nonnegative(les_kh_v(k,iCell))
            bv_safe = iamas_les_safe_finite(les_bv_freq2(k,iCell))
            shear_prod_h_value = km_h_safe * (d_diag*d_diag + d_off_diag*d_off_diag + divergence(k,iCell)**2)
            shear_prod_v_internal_value = cross_prod(k,iCell-cellSolveStart+1) + &
                 2.0_RKIND*km_v_safe*dw_dz*dw_dz
            shear_prod = shear_prod_h_value + shear_prod_v_internal_value
            les_shear_prod_h(k,iCell) = shear_prod_h_value
            les_shear_prod_v_internal(k,iCell) = shear_prod_v_internal_value
            buoy_prod = -kh_v_safe * bv_safe

            if (k == 1 .and. trim(config_les_surface_mode) == 'prescribed') then
               speed = sqrt(ur_cell(1,iCell)**2 + vr_cell(1,iCell)**2)
               surface_prod = max(config_les_drag_coefficient, 0.0_RKIND) * speed**3 / dz
               surface_prod_stress_shear = -les_uw_flux(1,iCell) * du_dz - &
                                            les_vw_flux(1,iCell) * dv_dz
               les_surface_prod_current(1,iCell) = surface_prod
               les_surface_prod_stress_shear(1,iCell) = surface_prod_stress_shear
               shear_prod = shear_prod + surface_prod
               if (associated(index_qv) .and. index_qv > 0) then
                  theta_dry = theta_m(1,iCell) / &
                              (1.0_RKIND + rvord * scalars(index_qv,1,iCell))
                  buoy_prod = 0.5_RKIND * (buoy_prod + gravity * config_les_heat_flux / &
                                          max(theta_dry, 100.0_RKIND))
               end if
            else if (k == 1 .and. trim(config_les_surface_mode) == 'physics_flux') then
               speed = sqrt(ur_cell(1,iCell)**2 + vr_cell(1,iCell)**2)
               surface_prod = iamas_les_safe_nonnegative(ust(iCell))**2 * speed / dz
               surface_prod_stress_shear = -les_uw_flux(1,iCell) * du_dz - &
                                            les_vw_flux(1,iCell) * dv_dz
               les_surface_prod_current(1,iCell) = surface_prod
               les_surface_prod_stress_shear(1,iCell) = surface_prod_stress_shear
               ! E1 minimum-change candidate: retain the legacy positive-definite
               ! surface work as a diagnostic, but use the stress--shear transfer
               ! in the SGS TKE source. The existing total-source limiter remains
               ! responsible for TKE positivity when this local transfer is negative.
               surface_prod = surface_prod_stress_shear
               shear_prod = shear_prod + surface_prod
               rho_surface = max(iamas_les_safe_nonnegative(zz(1,iCell) * &
                                   rho_zz(1,iCell)), 1.0e-6_RKIND)
               if (associated(cpm)) then
                  cpm_surface = max(iamas_les_safe_nonnegative(cpm(iCell)), 100.0_RKIND)
               else
                  cpm_surface = cp
               end if
               heat_flux_kinematic = iamas_les_safe_finite(hfx(iCell)) / &
                                     (rho_surface * cpm_surface)
               if (associated(index_qv) .and. index_qv > 0) then
                  theta_dry = theta_m(1,iCell) / &
                              (1.0_RKIND + rvord * max(scalars(index_qv,1,iCell), 0.0_RKIND))
                  buoy_prod = 0.5_RKIND * (buoy_prod + gravity * heat_flux_kinematic / &
                                          max(theta_dry, 100.0_RKIND))
               end if
            end if

            delta = max((areaCell(iCell) * dz)**(1.0_RKIND / 3.0_RKIND), 1.0_RKIND)
            if (k == 1 .or. k == nVertLevels) then
               ce = 3.9_RKIND
            else
               mix_safe = max(iamas_les_safe_nonnegative(les_mix_length(k,iCell)), 1.0_RKIND)
               ce = ce1 + ce2 * mix_safe / delta
            end if
            mix_safe = max(iamas_les_safe_nonnegative(les_mix_length(k,iCell)), 1.0_RKIND)
            tke = max(iamas_les_safe_nonnegative(scalars(index_tke_sgs,k,iCell)), config_les_tke_min)
            dissipation = -ce * tke**1.5_RKIND / mix_safe

            if (config_les_output_budget) then
               les_shear_prod(k,iCell) = shear_prod
               les_buoy_prod(k,iCell) = buoy_prod
               les_dissipation(k,iCell) = dissipation
            end if

            source_raw = iamas_les_safe_finite(shear_prod + buoy_prod + dissipation)
            source = max(source_raw, -(tke-config_les_tke_min) / dt)
            if (config_les_output_budget .and. associated(les_tke_tend_total)) then
               les_tke_tend_total(k,iCell) = iamas_les_safe_finite(les_tke_diffusion(k,iCell) + source)
            end if
            if (config_les_output_budget .and. associated(les_tke_source_pre_limiter)) then
               les_tke_source_pre_limiter(k,iCell) = source_raw
            end if
            if (config_les_output_budget .and. associated(les_tke_source_post_limiter)) then
               les_tke_source_post_limiter(k,iCell) = source
            end if
            if (config_les_output_budget .and. associated(les_tke_limiter_residual)) then
               les_tke_limiter_residual(k,iCell) = source - source_raw
            end if
            if (config_les_output_budget .and. associated(les_tke_budget_closure_residual) .and. &
                associated(les_tke_tend_total)) then
               budget_sum = les_tke_diffusion_h(k,iCell) + les_tke_diffusion_v(k,iCell) + &
                            shear_prod + buoy_prod + dissipation + (source - source_raw)
               les_tke_budget_closure_residual(k,iCell) = &
                    iamas_les_safe_finite(les_tke_tend_total(k,iCell) - budget_sum)
            end if
            if (apply_scalar_tendencies) then
               tend_scalars(index_tke_sgs,k,iCell) = tend_scalars(index_tke_sgs,k,iCell) + &
                                                      rho_zz(k,iCell) * source
            end if
         end do
      end do

      deallocate(vertical_flux,cross_prod)
   end subroutine iamas_les_add_sgs_tendencies


   ! The numerical scalar integrator uses this same nRelaxZone boundary.
   logical function iamas_les_affine_cell_active(apply_lbcs,mask)
      logical, intent(in) :: apply_lbcs
      integer, intent(in) :: mask
      iamas_les_affine_cell_active = (.not.apply_lbcs) .or. (mask<=nRelaxZone)
   end function iamas_les_affine_cell_active

   subroutine iamas_les_affine_abort(message)
#ifdef _MPI
      include 'mpif.h'
      integer :: ierr
#endif
      character(len=*), intent(in) :: message
      write(*,*) trim(message)
#ifdef _MPI
      ! All ranks in this standalone atmosphere executable must terminate.
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#endif
      error stop 'C2 LES affine fatal'
   end subroutine iamas_les_affine_abort

   ! E6MSA2C: per-block mesh cache, initialized once for each solved cell.
   ! This routine is called only after the les_deardorff early-return guard.
   subroutine iamas_les_prepare_affine(mesh,cellStart,cellEnd,apply_lbcs)
      type(iamas_pool_type), intent(in) :: mesh
      integer, intent(in) :: cellStart,cellEnd
      logical, intent(in) :: apply_lbcs
      integer, pointer :: bdyMask(:)
      integer, pointer :: nEdges
      integer, pointer :: edgesOnCell(:,:),nEdgesOnCell(:),edgesOnEdge(:,:),nEdgesOnEdge(:)
      real(kind=RKIND), pointer :: weightsOnEdge(:,:),areaCell(:),latCell(:),lonCell(:)
      real(kind=RKIND), pointer :: xCell(:),yCell(:),zCell(:),xEdge(:),yEdge(:),zEdge(:)
      real(kind=RKIND), pointer :: latEdge(:),lonEdge(:),angleEdge(:)
      real(kind=RKIND), pointer :: d1u(:,:),d1v(:,:),d2u(:,:),d2v(:,:),ready(:)
      real(kind=RKIND), allocatable :: edge_geometry(:,:),neighbor_geometry(:,:,:),weights(:,:)
      integer, allocatable :: neighbors(:)
      real(kind=RKIND) :: center(3),condition,max_condition
      integer :: iCell,i,j,edge,neighbor,n,neigh,number_built,number_specified,maxCellEdges,maxEdgeNeighbors
      logical :: ok

      call iamas_pool_get_array(mesh,'les_affine_condition',ready)
      if (.not.associated(ready)) call iamas_les_affine_abort('C LES affine cache absent')
      if (cellStart<1 .or. cellEnd>size(ready)) call iamas_les_affine_abort('C LES affine cell bounds')
      if (any(ready(cellStart:cellEnd)/=ready(cellStart:cellEnd))) call iamas_les_affine_abort('C LES affine cache NaN')
      call iamas_pool_get_array(mesh,'bdyMaskCell',bdyMask)
      if (.not.associated(bdyMask)) call iamas_les_affine_abort('C2 LES boundary mask absent')
      if (cellEnd>size(bdyMask)) call iamas_les_affine_abort('C2 LES boundary mask bounds')
      if (all(ready(cellStart:cellEnd)/=0.0_RKIND)) return
      call iamas_pool_get_dimension(mesh,'nEdges',nEdges)
      call iamas_pool_get_array(mesh,'edgesOnCell',edgesOnCell)
      call iamas_pool_get_array(mesh,'nEdgesOnCell',nEdgesOnCell)
      call iamas_pool_get_array(mesh,'edgesOnEdge',edgesOnEdge)
      call iamas_pool_get_array(mesh,'nEdgesOnEdge',nEdgesOnEdge)
      call iamas_pool_get_array(mesh,'weightsOnEdge',weightsOnEdge)
      call iamas_pool_get_array(mesh,'areaCell',areaCell)
      call iamas_pool_get_array(mesh,'latCell',latCell)
      call iamas_pool_get_array(mesh,'lonCell',lonCell)
      call iamas_pool_get_array(mesh,'xCell',xCell)
      call iamas_pool_get_array(mesh,'yCell',yCell)
      call iamas_pool_get_array(mesh,'zCell',zCell)
      call iamas_pool_get_array(mesh,'latEdge',latEdge)
      call iamas_pool_get_array(mesh,'lonEdge',lonEdge)
      call iamas_pool_get_array(mesh,'angleEdge',angleEdge)
      call iamas_pool_get_array(mesh,'xEdge',xEdge)
      call iamas_pool_get_array(mesh,'yEdge',yEdge)
      call iamas_pool_get_array(mesh,'zEdge',zEdge)
      call iamas_pool_get_array(mesh,'les_affine_d1u',d1u)
      call iamas_pool_get_array(mesh,'les_affine_d1v',d1v)
      call iamas_pool_get_array(mesh,'les_affine_d2u',d2u)
      call iamas_pool_get_array(mesh,'les_affine_d2v',d2v)
      maxCellEdges=size(edgesOnCell,1)
      maxEdgeNeighbors=size(edgesOnEdge,1)
      allocate(edge_geometry(6,maxCellEdges),neighbor_geometry(7,maxEdgeNeighbors,maxCellEdges))
      allocate(neighbors(maxCellEdges),weights(2,2*maxCellEdges))
      number_built=0
      number_specified=0
      max_condition=0.0_RKIND
      do iCell=cellStart,cellEnd
         if (.not.iamas_les_affine_cell_active(apply_lbcs,bdyMask(iCell))) then
            ! -1 explicitly marks a non-prognostic specified-boundary cell.
            ! It is NOT a fitted geometry condition number.
            d1u(:,iCell)=0.0_RKIND; d1v(:,iCell)=0.0_RKIND
            d2u(:,iCell)=0.0_RKIND; d2v(:,iCell)=0.0_RKIND
            ready(iCell)=-1.0_RKIND
            number_specified=number_specified+1
            cycle
         end if
         if (ready(iCell)<0.0_RKIND) call iamas_les_affine_abort('C2 active cell marked specified')
         if (ready(iCell)>0.0_RKIND) cycle
         n=nEdgesOnCell(iCell)
         if (n<3 .or. n>maxCellEdges) call iamas_les_affine_abort('C LES affine invalid cell stencil')
         if (areaCell(iCell)/=areaCell(iCell) .or. areaCell(iCell)<=0.0_RKIND) &
            call iamas_les_affine_abort('C LES affine invalid area')
         center=[xCell(iCell),yCell(iCell),zCell(iCell)]
         edge_geometry=0.0_RKIND
         neighbor_geometry=0.0_RKIND
         neighbors=0
         do i=1,n
            edge=edgesOnCell(i,iCell)
            if (edge<1 .or. edge>nEdges) call iamas_les_affine_abort('C LES affine edge halo missing')
            edge_geometry(:,i)=[xEdge(edge),yEdge(edge),zEdge(edge),latEdge(edge),lonEdge(edge),angleEdge(edge)]
            neigh=nEdgesOnEdge(edge)
            if (neigh<1 .or. neigh>maxEdgeNeighbors) call iamas_les_affine_abort('C LES affine neighbor count')
            neighbors(i)=neigh
            do j=1,neigh
               neighbor=edgesOnEdge(j,edge)
               if (neighbor<1 .or. neighbor>nEdges) then
                  write(*,*) 'C_AFFINE_MISSING_HALO cell edge neighbor=',iCell,edge,neighbor
                  call iamas_les_affine_abort('C LES affine neighbor halo missing')
               end if
               neighbor_geometry(:,j,i)=[xEdge(neighbor),yEdge(neighbor),zEdge(neighbor), &
                  latEdge(neighbor),lonEdge(neighbor),angleEdge(neighbor),weightsOnEdge(j,edge)]
            end do
         end do
         call iamas_les_affine_native(center,latCell(iCell),lonCell(iCell),sqrt(areaCell(iCell)), &
            edge_geometry(:,1:n),neighbors(1:n),neighbor_geometry(:,:,1:n),1,weights(:,1:2*n),condition,ok)
         if (.not.ok) then
            write(*,*) 'C_AFFINE_REJECT cell lat lon condition=',iCell,latCell(iCell),lonCell(iCell),condition
            call iamas_les_affine_abort('C LES affine geometry gate failed')
         end if
         do i=1,n
            d1u(i,iCell)=weights(1,2*i-1)
            d1v(i,iCell)=weights(1,2*i)
            d2u(i,iCell)=weights(2,2*i-1)
            d2v(i,iCell)=weights(2,2*i)
         end do
         ready(iCell)=condition
         max_condition=max(max_condition,condition)
         number_built=number_built+1
      end do
      write(*,*) 'C2_AFFINE_CACHE_PASS active specified max_condition=',number_built,number_specified,max_condition
      deallocate(edge_geometry,neighbor_geometry,neighbors,weights)
   end subroutine iamas_les_prepare_affine

! Geometry-only candidate kernel. RKIND is provided by the containing module.
! Parameters: U0,V0,L*Ux,L*Uy,L*Vx,L*Vy; observations are native edge u,v.
! Tangential rows MUST include weightsOnEdge and its neighboring normals.
! Coefficients reconstruct (Ux-Vy,Uy+Vx), not the full 3-D strain tensor.
! Rejected geometry is a fatal preflight condition, never a zero-strain fallback.
subroutine iamas_les_affine_weights(rows, length, weights, condition, ok)
   real(kind=RKIND), intent(in) :: rows(:,:), length
   real(kind=RKIND), intent(out) :: weights(:,:), condition
   logical, intent(out) :: ok
   real(kind=RKIND) :: gram(6,6), augmented(6,12), swap_row(12)
   real(kind=RKIND) :: inverse(6,6), selector(6,2), dual(6,2)
   real(kind=RKIND) :: largest, pivot_value, factor, gram_norm, inverse_norm
   integer :: i,k,pivot,m

   ok=.false.
   weights=0.0_RKIND
   condition=huge(1.0_RKIND)
   if (size(rows,1)/=6 .or. size(rows,2)<6) return
   m=size(rows,2)
   if (size(weights,1)/=2 .or. size(weights,2)/=m) return
   if (length/=length .or. length<=0.0_RKIND .or. length>1.0e9_RKIND) return
   if (any(rows/=rows) .or. any(abs(rows)>1.0e15_RKIND)) return

   gram=matmul(rows,transpose(rows))
   largest=maxval(abs(gram))
   if (largest<=tiny(1.0_RKIND)) return
   augmented=0.0_RKIND
   augmented(:,1:6)=gram
   do i=1,6
      augmented(i,i+6)=1.0_RKIND
   end do
   do k=1,6
      pivot=k
      do i=k+1,6
         if (abs(augmented(i,k))>abs(augmented(pivot,k))) pivot=i
      end do
      if (abs(augmented(pivot,k))<1.0e-12_RKIND*largest) return
      swap_row=augmented(k,:)
      augmented(k,:)=augmented(pivot,:)
      augmented(pivot,:)=swap_row
      pivot_value=augmented(k,k)
      augmented(k,:)=augmented(k,:)/pivot_value
      do i=1,6
         if (i==k) cycle
         factor=augmented(i,k)
         augmented(i,:)=augmented(i,:)-factor*augmented(k,:)
      end do
   end do
   inverse=augmented(:,7:12)
   gram_norm=0.0_RKIND
   inverse_norm=0.0_RKIND
   do i=1,6
      gram_norm=max(gram_norm,sum(abs(gram(i,:))))
      inverse_norm=max(inverse_norm,sum(abs(inverse(i,:))))
   end do
   condition=gram_norm*inverse_norm
   if (condition/=condition .or. condition>1.0e8_RKIND) return
   selector=0.0_RKIND
   selector(3,1)=1.0_RKIND/length
   selector(6,1)=-1.0_RKIND/length
   selector(4,2)=1.0_RKIND/length
   selector(5,2)=1.0_RKIND/length
   dual=matmul(inverse,selector)
   weights=matmul(transpose(dual),rows)
   if (any(weights/=weights)) return
   ! Independent affine moment contract on both output rows.
   if (maxval(abs(matmul(weights,transpose(rows))-transpose(selector))) &
       >1.0e-10_RKIND/length) return
   ok=.true.
end subroutine iamas_les_affine_weights

! Native geometry adapter for the offline candidate kernel; no model state writes.
! Arguments contain only one cell's existing edge stencil, so no global array copies.
subroutine iamas_les_affine_frame(lat,lon,angle,normal,tangent)
   real(kind=RKIND), intent(in) :: lat,lon,angle
   real(kind=RKIND), intent(out) :: normal(3),tangent(3)
   real(kind=RKIND) :: east(3),north(3)
   east=[-sin(lon),cos(lon),0.0_RKIND]
   north=[-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat)]
   normal=cos(angle)*east+sin(angle)*north
   tangent=-sin(angle)*east+cos(angle)*north
end subroutine iamas_les_affine_frame

subroutine iamas_les_affine_row(position,direction,center,east,north,length,row)
   real(kind=RKIND), intent(in) :: position(3),direction(3),center(3),east(3),north(3),length
   real(kind=RKIND), intent(out) :: row(6)
   real(kind=RKIND) :: x,y,e,n
   x=dot_product(position-center,east)/length
   y=dot_product(position-center,north)/length
   e=dot_product(direction,east)
   n=dot_product(direction,north)
   row=[e,n,e*x,e*y,n*x,n*y]
end subroutine iamas_les_affine_row

subroutine iamas_les_affine_native(center,cell_lat,cell_lon,length,edge_geometry, &
                                  neighbors,neighbor_geometry,mode,weights,condition,ok)
   real(kind=RKIND), intent(in) :: center(3),cell_lat,cell_lon,length
   ! Each geometry column: x,y,z,latitude,longitude,angle; neighbor adds weight.
   real(kind=RKIND), intent(in) :: edge_geometry(:,:),neighbor_geometry(:,:,:)
   integer, intent(in) :: neighbors(:),mode
   real(kind=RKIND), intent(out) :: weights(:,:),condition
   logical, intent(out) :: ok
   real(kind=RKIND) :: east(3),north(3),normal(3),tangent(3),row(6)
   real(kind=RKIND), allocatable :: rows(:,:)
   integer :: i,j,n
   ok=.false.
   weights=0.0_RKIND
   condition=huge(1.0_RKIND)
   n=size(edge_geometry,2)
   if (size(edge_geometry,1)/=6 .or. n<3) return
   if (size(neighbors)/=n .or. size(neighbor_geometry,1)/=7) return
   if (size(neighbor_geometry,3)/=n) return
   if (size(weights,1)/=2 .or. size(weights,2)/=2*n) return
   if (mode/=0 .and. mode/=1) return
   if (length/=length .or. length<=0.0_RKIND) return
   if (any(neighbors<0) .or. any(neighbors>size(neighbor_geometry,2))) return
   call iamas_les_affine_frame(cell_lat,cell_lon,0.0_RKIND,east,north)
   allocate(rows(6,2*n))
   rows=0.0_RKIND
   do i=1,n
      call iamas_les_affine_frame(edge_geometry(4,i),edge_geometry(5,i),edge_geometry(6,i),normal,tangent)
      call iamas_les_affine_row(edge_geometry(1:3,i),normal,center,east,north,length,rows(:,2*i-1))
      if (mode==0) then
         call iamas_les_affine_row(edge_geometry(1:3,i),tangent,center,east,north,length,rows(:,2*i))
      else
         do j=1,neighbors(i)
            call iamas_les_affine_frame(neighbor_geometry(4,j,i),neighbor_geometry(5,j,i), &
                                         neighbor_geometry(6,j,i),normal,tangent)
            call iamas_les_affine_row(neighbor_geometry(1:3,j,i),normal,center,east,north,length,row)
            rows(:,2*i)=rows(:,2*i)+neighbor_geometry(7,j,i)*row
         end do
      end if
   end do
   call iamas_les_affine_weights(rows,length,weights,condition,ok)
   deallocate(rows)
end subroutine iamas_les_affine_native

   subroutine iamas_les_cross_apply(mesh,configs,nz,cs,ce,es,ee,rhozz,u,w,zz,rdzw,kmv, &
       tu,tue,tw,twe,dun,dw,unw,uw,vw,production)
      type(iamas_pool_type), intent(in) :: mesh,configs
      integer, intent(in) :: nz,cs,ce,es,ee
      real(kind=RKIND), intent(in) :: rhozz(:,:),u(:,:),w(:,:),zz(:,:),rdzw(:),kmv(:,:)
      real(kind=RKIND), intent(inout) :: tu(:,:),tue(:,:),tw(:,:),twe(:,:),dun(:,:),dw(:,:)
      real(kind=RKIND), intent(inout) :: unw(:,:),uw(:,:),vw(:,:)
      real(kind=RKIND), intent(out) :: production(nz,ce-cs+1)
      type(les_cross_cache_type), pointer :: a
      integer, pointer :: nc0,ne0,coe(:,:),eoc(:,:),nec(:),gids(:),bmc(:),bme(:)
      real(kind=RKIND), pointer :: z(:,:),dc(:),dv(:),area(:),latc(:),lonc(:),late(:),lone(:),angle(:)
      logical, pointer :: lbcs
      logical, allocatable :: wantc(:),wante(:)
      integer, allocatable :: reverse_c(:),reverse_e(:)
      real(kind=RKIND), allocatable :: local_z(:,:)
      integer :: i,c,e,k,j,p,c1,c2,noc,noe,status
      real(kind=RKIND) :: delta,loss,dzetaw,rhoi,normal(3),tangent(3),east(3),north(3)
      real(kind=RKIND) :: mxx,mxy,myy,det,nx,ny,weight,kinematic
      production=0.0_RKIND
      call iamas_pool_get_array(mesh,'zgrid',z)
      a=>les_cross_head
      do while(associated(a))
         if (associated(a%zkey,z)) then
            if (a%cs==cs.and.a%ce==ce.and.a%es==es.and.a%ee==ee) exit
         end if
         a=>a%next
      end do
      if (.not.associated(a)) then
         call iamas_pool_get_dimension(mesh,'nCells',nc0)
         call iamas_pool_get_dimension(mesh,'nEdges',ne0)
         call iamas_pool_get_array(mesh,'cellsOnEdge',coe)
         call iamas_pool_get_array(mesh,'edgesOnCell',eoc)
         call iamas_pool_get_array(mesh,'nEdgesOnCell',nec)
         call iamas_pool_get_array(mesh,'indexToEdgeID',gids)
         call iamas_pool_get_array(mesh,'dcEdge',dc)
         call iamas_pool_get_array(mesh,'dvEdge',dv)
         call iamas_pool_get_array(mesh,'areaCell',area)
         call iamas_pool_get_array(mesh,'bdyMaskCell',bmc)
         call iamas_pool_get_array(mesh,'bdyMaskEdge',bme)
         call iamas_pool_get_config(configs,'config_apply_lbcs',lbcs)
         if (.not.associated(lbcs).or..not.associated(bmc).or..not.associated(bme)) &
            call iamas_les_affine_abort('C3 cross boundary contract absent')
         allocate(wantc(nc0),wante(ne0),reverse_c(nc0),reverse_e(ne0))
         wantc=.false.;wante=.false.;reverse_c=0;reverse_e=0;noc=0;noe=0
         do c=cs,ce
            if (.not.iamas_les_affine_cell_active(lbcs,bmc(c))) cycle
            noc=noc+1
            do i=1,nec(c)
               e=eoc(i,c)
               if (e<1.or.e>ne0) call iamas_les_affine_abort('C3 cross required edge halo absent')
               wante(e)=.true.
            end do
         end do
         do e=es,ee
            if (.not.iamas_les_affine_cell_active(lbcs,bme(e))) cycle
            wante(e)=.true.;noe=noe+1
         end do
         if (count(wante)==0) return
         do e=1,ne0
            if (.not.wante(e)) cycle
            if (any(coe(:,e)<1).or.any(coe(:,e)>nc0)) &
               call iamas_les_affine_abort('C3 cross required cell halo absent')
            wantc(coe(:,e))=.true.
         end do
         allocate(a);a%zkey=>z;a%cs=cs;a%ce=ce;a%es=es;a%ee=ee;a%nz=nz
         a%nc=count(wantc);a%ne=count(wante)
         allocate(a%cm(a%nc),a%em(a%ne),a%cells(2,a%ne),a%gids(a%ne),a%oc(noc),a%oe(noe))
         j=0
         do c=1,nc0
            if (.not.wantc(c)) cycle
            j=j+1;a%cm(j)=c;reverse_c(c)=j
         end do
         j=0
         do e=1,ne0
            if (.not.wante(e)) cycle
            j=j+1;a%em(j)=e;reverse_e(e)=j;a%cells(:,j)=reverse_c(coe(:,e));a%gids(j)=gids(e)
         end do
         j=0
         do c=cs,ce
            if (.not.iamas_les_affine_cell_active(lbcs,bmc(c))) cycle
            j=j+1;a%oc(j)=reverse_c(c)
         end do
         j=0
         do e=es,ee
            if (.not.iamas_les_affine_cell_active(lbcs,bme(e))) cycle
            j=j+1;a%oe(j)=reverse_e(e)
         end do
         allocate(a%h(nz-1,a%ne),a%coef(3,nz-1,2,a%ne),a%volume(nz-1,a%ne),a%dc(a%ne),a%dv(a%ne))
         allocate(a%offset(a%nc+1),a%elist(2*a%ne),a%slist(2*a%ne),a%ru(2*a%ne),a%rv(2*a%ne))
         allocate(local_z(nz+1,a%nc));local_z=z(:,a%cm);a%dc=dc(a%em);a%dv=dv(a%em)
         call cross_prepare(a%nc,a%ne,nz,a%cells,a%gids,local_z,a%dc,a%dv,a%h,a%coef,a%volume, &
                            a%offset,a%elist,a%slist,status)
         if (status/=0) call iamas_les_affine_abort('C3 cross geometry preparation failed')
         call iamas_pool_get_array(mesh,'latCell',latc)
         call iamas_pool_get_array(mesh,'lonCell',lonc)
         call iamas_pool_get_array(mesh,'latEdge',late)
         call iamas_pool_get_array(mesh,'lonEdge',lone)
         call iamas_pool_get_array(mesh,'angleEdge',angle)
         a%ru=0.0_RKIND;a%rv=0.0_RKIND
         do i=1,size(a%oc)
            c=a%oc(i);j=a%cm(c)
            call iamas_les_affine_frame(latc(j),lonc(j),0.0_RKIND,east,north)
            mxx=0.0_RKIND;mxy=0.0_RKIND;myy=0.0_RKIND
            do p=a%offset(c),a%offset(c+1)-1
               e=a%em(a%elist(p))
               call iamas_les_affine_frame(late(e),lone(e),angle(e),normal,tangent)
               nx=dot_product(normal,east);ny=dot_product(normal,north);weight=dc(e)*dv(e)
               a%ru(p)=weight*nx;a%rv(p)=weight*ny
               mxx=mxx+weight*nx*nx;mxy=mxy+weight*nx*ny;myy=myy+weight*ny*ny
            end do
            det=mxx*myy-mxy*mxy
            if (det<=1.e-12_RKIND*mxx*myy) call iamas_les_affine_abort('C3 flux reconstruction singular')
            do p=a%offset(c),a%offset(c+1)-1
               nx=a%ru(p);ny=a%rv(p)
               a%ru(p)=(myy*nx-mxy*ny)/det;a%rv(p)=(mxx*ny-mxy*nx)/det
            end do
         end do
         allocate(a%rho(nz,a%nc),a%u(nz,a%ne),a%w(nz+1,a%nc),a%kv(nz-1,a%ne))
         allocate(a%shear(nz-1,a%ne),a%fu(nz,a%ne),a%fw(nz+1,a%nc),a%prod(nz,a%nc))
         a%next=>les_cross_head;les_cross_head=>a
         write(*,*) 'C3_CROSS_CACHE_READY columns edges owners=',a%nc,a%ne,noc,noe
      end if
      call iamas_pool_get_array(mesh,'areaCell',area)
      do i=1,a%nc
         c=a%cm(i);a%rho(:,i)=rhozz(:,c)*zz(:,c);a%w(:,i)=w(:,c)
      end do
      do i=1,a%ne
         e=a%em(i);a%u(:,i)=u(:,e);c1=a%cm(a%cells(1,i));c2=a%cm(a%cells(2,i))
         do k=1,nz-1
            a%kv(k,i)=.25_RKIND*(kmv(k,c1)+kmv(k+1,c1)+kmv(k,c2)+kmv(k+1,c2))
         end do
      end do
      call cross_gather(a%nc,a%ne,nz,size(a%oc),size(a%oe),a%cells,a%h,a%coef,a%volume, &
          a%offset,a%elist,a%slist,a%rho,a%u,a%w,a%kv,a%oc,a%oe,a%shear,a%fu,a%fw,a%prod,loss,status)
      if (status/=0) call iamas_les_affine_abort('C3 cross dynamic input/force invalid')
      do i=1,size(a%oe)
         j=a%oe(i);e=a%em(j);c1=a%cm(a%cells(1,j));c2=a%cm(a%cells(2,j))
         do k=1,nz
            delta=a%fu(k,j)*rdzw(k)/(a%dc(j)*a%dv(j))
            tu(k,e)=tu(k,e)+delta;tue(k,e)=tue(k,e)+delta
            rhoi=.5_RKIND*(rhozz(k,c1)+rhozz(k,c2))
            dun(k,e)=dun(k,e)+delta/max(rhoi,1.e-12_RKIND)
         end do
         unw(2:nz,e)=-a%kv(:,j)*a%shear(:,j)
      end do
      do i=1,size(a%oc)
         j=a%oc(i);c=a%cm(j)
         do k=1,nz
            production(k,c-cs+1)=a%prod(k,j)*rdzw(k)/(area(c)*max(rhozz(k,c),1.e-12_RKIND))
         end do
         do k=2,nz
            dzetaw=.5_RKIND*(1.0_RKIND/rdzw(k-1)+1.0_RKIND/rdzw(k))
            delta=a%fw(k,j)/(area(c)*dzetaw)
            tw(k,c)=tw(k,c)+delta;twe(k,c)=twe(k,c)+delta
            rhoi=.5_RKIND*(rhozz(k-1,c)+rhozz(k,c))
            dw(k,c)=dw(k,c)+delta/max(rhoi,1.e-12_RKIND)
            uw(k,c)=0.0_RKIND;vw(k,c)=0.0_RKIND
            do p=a%offset(j),a%offset(j+1)-1
               e=a%elist(p);kinematic=-a%kv(k-1,e)*a%shear(k-1,e)
               uw(k,c)=uw(k,c)+a%ru(p)*kinematic;vw(k,c)=vw(k,c)+a%rv(p)*kinematic
            end do
         end do
      end do
   end subroutine iamas_les_cross_apply

end module iamas_les_sgs
