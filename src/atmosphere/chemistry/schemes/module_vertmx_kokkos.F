! Batched Kokkos entry for the existing dry-deposition vertmx loop.
! The original module_vertmx_wrf remains the numerical reference.
module module_vertmx_kokkos
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use flcl_ndarray_mod, only: nd_array_t, to_nd_array
#ifndef VERTMX_STANDALONE
  use iamas_timer, only: iamas_timer_start, iamas_timer_stop
#endif
  implicit none
  private
  public :: vertmx_kokkos_batch, vertmx_kokkos_enabled, vertmx_kokkos_mode
  interface
    subroutine c_vertmx_batch(dt, eps, gas_only, chem, alt, zgrid, kzh, ddvel, aerosol) &
        bind(C, name='iamas_vertmx_kokkos_batch')
      import :: c_double, c_int, nd_array_t
      real(c_double), value :: dt, eps
      integer(c_int), value :: gas_only
      type(nd_array_t), intent(in) :: chem, alt, zgrid, kzh, ddvel, aerosol
    end subroutine
  end interface
contains
  integer function vertmx_kokkos_mode()
    integer, save :: mode = -1
    character(len=16) :: text
    integer :: status, length
    if (mode == -1) then
      call get_environment_variable('IAMAS_VERTMX_MODE', text, length=length, status=status)
      mode = 1
      if (status == 0 .and. length > 0) then
        select case (trim(text))
        case ('0')
          mode = 0
        case ('1')
          mode = 1
        case ('2')
          mode = 2
        case default
          error stop 'IAMAS_VERTMX_MODE must be 0 (Fortran), 1 (Kokkos), or 2 (validate)'
        end select
      else if (status /= 1 .and. status /= 0) then
        error stop 'Cannot read IAMAS_VERTMX_MODE'
      end if
    end if
    vertmx_kokkos_mode = mode
  end function

  logical function vertmx_kokkos_enabled()
    vertmx_kokkos_enabled = vertmx_kokkos_mode() /= 0
  end function

  subroutine vertmx_kokkos_batch(dt, eps, gas_only, chem, alt, zgrid, kzh, ddvel, is_aerosol)
    real(c_double), intent(in) :: dt, eps
    logical, intent(in) :: gas_only
    real(c_double), intent(inout) :: chem(:,:,:)
    real(c_double), intent(in) :: alt(:,:), zgrid(:,:), kzh(:,:), ddvel(:,:)
    logical, intent(in) :: is_aerosol(:)
    integer(c_int) :: aerosol(size(is_aerosol)), gas
    if (size(chem,2) < 2) error stop 'vertmx requires at least two vertical levels'
    if (size(chem,1) /= size(is_aerosol)) error stop 'vertmx species mask mismatch'
    if (size(alt,1) /= size(chem,2) .or. size(alt,2) /= size(chem,3)) error stop 'vertmx alt shape mismatch'
    if (size(zgrid,1) /= size(chem,2)+1 .or. size(zgrid,2) /= size(chem,3)) error stop 'vertmx zgrid shape mismatch'
    if (size(kzh,1) /= size(chem,2) .or. size(kzh,2) /= size(chem,3)) error stop 'vertmx kzh shape mismatch'
    if (size(ddvel,1) /= size(chem,1) .or. size(ddvel,2) /= size(chem,3)) error stop 'vertmx ddvel shape mismatch'
    aerosol = merge(1_c_int, 0_c_int, is_aerosol)
    gas = merge(1_c_int, 0_c_int, gas_only)
    call c_vertmx_batch(dt, eps, gas, to_nd_array(chem), to_nd_array(alt), &
         to_nd_array(zgrid), to_nd_array(kzh), to_nd_array(ddvel), to_nd_array(aerosol))
  end subroutine

  subroutine c_vertmx_stage(stage, begin) bind(C, name='iamas_vertmx_stage')
    integer(c_int), value :: stage, begin
#ifndef VERTMX_STANDALONE
    character(len=32) :: label
    select case(stage)
    case (0)
      label = 'vertmx_workspace'
    case (1)
      label = 'vertmx_geometry'
    case (2)
      label = 'vertmx_solve_writeback'
    case default
      error stop 'Invalid vertmx timing stage'
    end select
    if (begin /= 0) then
      call iamas_timer_start(trim(label))
    else
      call iamas_timer_stop(trim(label))
    end if
#endif
  end subroutine
end module module_vertmx_kokkos
