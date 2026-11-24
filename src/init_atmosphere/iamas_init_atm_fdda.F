! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!==================================================================================================
module iamas_init_atm_fdda

   use iamas_kind_types
   use iamas_derived_types
   use iamas_pool_routines
   use iamas_constants
   use iamas_dmpar
   use iamas_sort
   use atm_advection
   use iamas_RBF_interpolation
   use iamas_vector_reconstruction
   use iamas_timer
   use iamas_atmphys_constants, only: svpt0,svp1,svp2,svp3
   use iamas_atmphys_functions
   use iamas_atmphys_initialize_real
   use iamas_log, only : iamas_log_write

   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in iamas_timekeeping
   use iamas_timekeeping !, only: iAMAS_Time_type, iAMAS_TimeInterval_type, iAMAS_Clock_type, &
                        !        iamas_set_time, iamas_set_timeInterval, iamas_get_time, operator(+), add_t_ti



implicit none
private
public :: init_atm_case_fdda


contains


!==================================================================================================
subroutine init_atm_case_fdda(block, domain, dminfo, stream_manager)
!==================================================================================================

   use iamas_stream_manager
   
   implicit none
   
   !input arguments:
   type (block_type), intent(inout), target :: block
   type (domain_type), intent(inout) :: domain
   type (dm_info), intent(in)        :: dminfo
   type (iAMAS_streamManager_type), intent(inout) :: stream_manager
   
   !local variables:
   type (iAMAS_Time_type)  :: curr_time, stop_time
   type (iamas_pool_type), pointer :: mesh, fg, state, diag, fdda
   integer, pointer :: nCells, nEdges, nVertLevels
   character(len=StrKIND) :: timeString
   character(len=StrKIND), pointer :: xtime
   integer :: ierr
   
   call iamas_pool_get_subpool(block % structs, 'mesh', mesh)
   call iamas_pool_get_subpool(block % structs, 'fg', fg)
   call iamas_pool_get_subpool(block % structs, 'state', state)
   call iamas_pool_get_subpool(block % structs, 'diag', diag)
   call iamas_pool_get_subpool(block % structs, 'fdda', fdda)

   call iamas_pool_get_dimension(block % dimensions, 'nCells', nCells)
   call iamas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
   call iamas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
   
   !==================================================================================================
   call iamas_pool_get_array(state, 'xtime', xtime)
   
   !loop over all times:
   curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW) 
   stop_time = iamas_get_clock_time(domain % clock, iAMAS_STOP_TIME) 
   
   do while (curr_time <= stop_time)
      call iamas_get_time(curr_time, dateTimeString=timeString)
      xtime = timeString
      call iamas_log_write('')
      call iamas_log_write('*************************************************')
      call iamas_log_write('--- Timing for FDDA case at: '//trim(timeString))

      !data to the iAMAS grid:
      call interp_fdda_to_iAMAS(timeString, block, mesh, nCells, nEdges, nVertLevels, fg, state, &
                               diag, fdda, block % dimensions, block % configs)
   
      !write the interpolated SST/SKINTEMP field as a new time slice in the iAMAS output file:
      call iamas_stream_mgr_write(stream_manager, streamID='fdda', ierr=ierr)
      call iamas_stream_mgr_reset_alarms(stream_manager, streamID='fdda', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)
   
      call iamas_advance_clock(domain % clock)
      curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW) 
   end do
   !
   ! Ensure that no output alarms are still ringing for the 'fdda' stream after
   ! we exit the time loop above; the main run routine may write out all other
   ! output streams with ringing alarms.
   !
   call iamas_stream_mgr_reset_alarms(stream_manager, streamID='fdda', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

end subroutine init_atm_case_fdda

subroutine interp_fdda_to_iAMAS(timeString, block, mesh, nCells, nEdges, nVertLevels, fg, state, diag, fdda, dims, configs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real-data test case using GFS data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use iamas_dmpar
   use init_atm_read_met
   use init_atm_llxy
   use init_atm_hinterp
   use iamas_hash
   use iamas_atmphys_constants, only : svpt0, svp1, svp2, svp3

   implicit none

   character(len=StrKIND), intent(in) :: timeString
   type (block_type), intent(inout), target :: block
   type (iamas_pool_type), intent(inout) :: mesh
   integer, intent(in) :: nCells
   integer, intent(in) :: nEdges
   integer, intent(in) :: nVertLevels
   type (iamas_pool_type), intent(inout) :: fg
   type (iamas_pool_type), intent(inout) :: state
   type (iamas_pool_type), intent(inout) :: diag
   type (iamas_pool_type), intent(inout) :: fdda
   type (iamas_pool_type), intent(inout):: dims
   type (iamas_pool_type), intent(inout):: configs

   type (parallel_info), pointer :: parinfo
   type (dm_info), pointer :: dminfo

   real (kind=RKIND), parameter :: u0 = 35.0
   real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

   real (kind=RKIND) :: omega_e

   real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
   real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
   real (kind=RKIND), parameter :: theta_c = pii/4.0
   real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
   real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
   real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

   type (met_data) :: field
   type (proj_info) :: proj

   real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
   real (kind=RKIND), dimension(:), pointer :: vert_level, latPoints, lonPoints, ter
   real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
   real (kind=RKIND), dimension(:,:), pointer :: pressure, pressure_b, exner_b, rho_zz, rho_b, rho_p, phi, phi_b, phi_total
   real (kind=RKIND), dimension(:,:), pointer :: theta_b, rtheta_b, exner, pressure_p, dss, theta_m, rtheta_p
   real (kind=RKIND), dimension(:), pointer :: surface_pressure
   real (kind=RKIND), dimension(:), pointer :: destField1d
   real (kind=RKIND), dimension(:,:), pointer :: destField2d
   real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
   real (kind=RKIND), dimension(:,:,:), pointer :: scalars
   real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

   real (kind=RKIND) :: target_z
   integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
   integer, pointer :: nCellsSolve, nz1
   integer :: nInterpPoints, ndims

   integer :: nfglevels_actual
   integer, pointer :: index_qv

   integer, dimension(5) :: interp_list
   real (kind=RKIND) :: maskval
   real (kind=RKIND) :: msgval
   real (kind=RKIND) :: fillval
   integer :: masked

   !This is temporary variable here. It just need when calculate tangential velocity v.
   integer :: eoe, j
   integer, dimension(:), pointer :: nEdgesOnCell
   integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, edgesOnCell, cellsOnCell
   real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
   real (kind=RKIND), dimension(:,:), pointer :: v
   real (kind=RKIND), dimension(:,:), pointer :: sorted_arr
   integer, dimension(:), pointer :: bdyMaskCell

   type (field1DReal), pointer :: tempField
   type (field1DReal), pointer :: ter_field
   type (field1DReal), target :: tempFieldTarget

   real(kind=RKIND), dimension(:), pointer :: hs, hs1, sm0
   real(kind=RKIND) :: hm, hm_global, zh, dzmin, dzmina, dzmina_global, dzminf, dzminf_global, sm
   real(kind=RKIND) :: dcsum
   integer :: nsmterrain, kz, sfc_k
   logical :: hybrid, smooth,cam_lev30,cam_lev32

   integer :: it
   real (kind=RKIND) :: p_check

   ! For interpolating terrain and land use
   integer :: nx, ny
   integer :: istatus

   real (kind=RKIND), allocatable, dimension(:,:) :: rslab, maskslab
   integer, dimension(:), pointer :: mask_array
   integer, dimension(nEdges), target :: edge_mask
   character (len=StrKIND) :: fname
   logical :: is_sfc_field

   real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r
   real (kind=RKIND) :: lat, lon, x, y

   real (kind=RKIND) :: ptop, p0
   real (kind=RKIND) :: lon_Edge

   real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

   real (kind=RKIND), dimension(nVertLevels, nCells) :: rel_hum, temperature, qv
   real (kind=RKIND) :: ptmp, es, rs, rgas_moist, qvs, xnutr, znut, ptemp, rcv
   integer :: iter

   real (kind=RKIND), dimension(nVertLevels + 1) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
   real (kind=RKIND), dimension(nVertLevels + 1) :: znuc, znuv, bn, divh, dpn

   real (kind=RKIND), dimension(:), pointer :: specified_zw
   real (kind=RKIND), dimension(nVertLevels + 1) :: sh, zw, ah
   real (kind=RKIND), dimension(nVertLevels) :: zu, dzw, rdzwp, rdzwm
   real (kind=RKIND), dimension(nVertLevels) :: eta, etav, teta, ppi, tt

   real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf

   !  storage for (lat,z) arrays for zonal velocity calculation

   integer, parameter :: nlat=361
   real (kind=RKIND), dimension(nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
   real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
   real (kind=RKIND), dimension(nlat, nVertLevels) :: u_2d, etavs_2d
   real (kind=RKIND), dimension(nVertLevels + 1) :: fsum
   real (kind=RKIND), dimension(nlat) :: lat_2d
   real (kind=RKIND) :: dlat
   real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2
   real (kind=RKIND) :: alt, als, zetal, zl

   !  calculation of the water vapor mixing ratio:
   real (kind=RKIND) :: sh_max,sh_min,global_sh_max,global_sh_min

   character (len=StrKIND), pointer :: config_fdda_prefix
   character (len=StrKIND), pointer :: config_start_time
   ! added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<<
   logical, pointer :: config_met_interp
   logical, pointer :: config_log_pressure
   logical, pointer :: config_vertical_grid
   integer, pointer :: config_nsmterrain
   integer, pointer :: config_nsm
   real (kind=RKIND), pointer :: config_dzmin
   real (kind=RKIND), pointer :: config_ztop
   logical, pointer :: config_tc_vertical_grid
   character (len=StrKIND), pointer :: config_specified_zeta_levels
   logical, pointer :: config_use_spechumd
   integer, pointer :: config_nfglevels
   integer, pointer :: config_nfgsoillevels
   logical, pointer :: config_smooth_surfaces
   integer, pointer :: config_theta_adv_order
   real (kind=RKIND), pointer :: config_coef_3rd_order
   logical, pointer :: config_blend_bdy_terrain

   character (len=StrKIND), pointer :: config_extrap_airtemp
   integer :: extrap_airtemp

   real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
   real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
   real (kind=RKIND), dimension(:), pointer :: angleEdge
   real (kind=RKIND), pointer :: cf1, cf2, cf3
   integer, dimension(:), pointer :: landmask

   real (kind=RKIND), dimension(:,:), pointer :: dzs_fg
   real (kind=RKIND), dimension(:,:), pointer :: zs_fg

   real (kind=RKIND), dimension(:), pointer :: sst
   real (kind=RKIND), dimension(:), pointer :: seaice
   real (kind=RKIND), dimension(:), pointer :: xice
   real (kind=RKIND), dimension(:,:), pointer :: u
   real (kind=RKIND), dimension(:,:), pointer :: w
   real (kind=RKIND), dimension(:,:), pointer :: theta
   real (kind=RKIND), dimension(:,:), pointer :: rho
   real (kind=RKIND), dimension(:,:), pointer :: relhum
   real (kind=RKIND), dimension(:,:), pointer :: spechum
   real (kind=RKIND), dimension(:,:), pointer :: ru
   real (kind=RKIND), dimension(:,:), pointer :: rw
   real (kind=RKIND), dimension(:), pointer :: precipw
   real (kind=RKIND), dimension(:,:), pointer :: uReconstructX
   real (kind=RKIND), dimension(:,:), pointer :: uReconstructY
   real (kind=RKIND), dimension(:,:), pointer :: uReconstructZ
   real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal
   real (kind=RKIND), dimension(:,:), pointer :: uReconstructMeridional
   real (kind=RKIND), dimension(:,:), pointer :: qv_fdda_new
   real (kind=RKIND), dimension(:,:), pointer :: th_fdda_new
   real (kind=RKIND), dimension(:,:), pointer :: u_fdda_new
   real (kind=RKIND), dimension(:,:), pointer :: v_fdda_new

   real (kind=RKIND), dimension(:), pointer :: psfc
   real (kind=RKIND), dimension(:), pointer :: skintemp
   real (kind=RKIND), dimension(:), pointer :: snow
   real (kind=RKIND), dimension(:), pointer :: snowc
   real (kind=RKIND), dimension(:,:), pointer :: u_fg
   real (kind=RKIND), dimension(:,:), pointer :: v_fg
   real (kind=RKIND), dimension(:,:), pointer :: z_fg
   real (kind=RKIND), dimension(:,:), pointer :: t_fg
   real (kind=RKIND), dimension(:,:), pointer :: rh_fg
   real (kind=RKIND), dimension(:,:), pointer :: sh_fg
   real (kind=RKIND), dimension(:,:), pointer :: gfs_z
   real (kind=RKIND), dimension(:,:), pointer :: p_fg
   real (kind=RKIND), dimension(:,:), pointer :: erod
   real (kind=RKIND), dimension(:), pointer :: AH_emis! added by qkyang to support anthropogenic heating emissions
   real (kind=RKIND), dimension(:), pointer :: dms_0
   real (kind=RKIND), dimension(:), pointer :: soilz
   type (hashtable) :: level_hash
   logical :: too_many_fg_levs
   integer :: level_value

   character (len=StrKIND) :: errstring

   call iamas_pool_get_config(configs, 'config_fdda_prefix', config_fdda_prefix)
   call iamas_pool_get_config(configs, 'config_met_interp', config_met_interp)
   call iamas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
   call iamas_pool_get_config(configs, 'config_nsmterrain', config_nsmterrain)
   call iamas_pool_get_config(configs, 'config_nsm', config_nsm)
   call iamas_pool_get_config(configs, 'config_dzmin', config_dzmin)
   call iamas_pool_get_config(configs, 'config_ztop', config_ztop)
   call iamas_pool_get_config(configs, 'config_tc_vertical_grid', config_tc_vertical_grid)
   call iamas_pool_get_config(configs, 'config_specified_zeta_levels', config_specified_zeta_levels)
   call iamas_pool_get_config(configs, 'config_use_spechumd', config_use_spechumd)
   call iamas_pool_get_config(configs, 'config_nfglevels', config_nfglevels)
   call iamas_pool_get_config(configs, 'config_nfgsoillevels', config_nfgsoillevels)
   call iamas_pool_get_config(configs, 'config_smooth_surfaces', config_smooth_surfaces)
   call iamas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)
   call iamas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
   call iamas_pool_get_config(configs, 'config_blend_bdy_terrain', config_blend_bdy_terrain)
   call iamas_pool_get_config(configs, 'config_log_pressure', config_log_pressure)
   
   call iamas_pool_get_config(configs, 'config_extrap_airtemp', config_extrap_airtemp)
   if (trim(config_extrap_airtemp) == 'constant') then
      extrap_airtemp = 0
   else if (trim(config_extrap_airtemp) == 'linear') then
      extrap_airtemp = 1
   else if (trim(config_extrap_airtemp) == 'lapse-rate') then
      extrap_airtemp = 2
   else
       call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
       call iamas_log_write('* Invalid value for namelist variable config_extrap_airtemp *', messageType=iAMAS_LOG_ERR)
       call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
   end if

   parinfo => block % parinfo
   dminfo => block % domain % dminfo

   call iamas_pool_get_field(mesh, 'ter', ter_field)

   call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
   call iamas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
   call iamas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
   call iamas_pool_get_array(mesh, 'dvEdge', dvEdge)
   call iamas_pool_get_array(mesh, 'dcEdge', dcEdge)
   call iamas_pool_get_array(mesh, 'areaCell', areaCell)
   call iamas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
   call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
   call iamas_pool_get_array(mesh, 'angleEdge', angleEdge)

   call iamas_pool_get_array(mesh, 'deriv_two', deriv_two)
   call iamas_pool_get_array(mesh, 'zb', zb)
   call iamas_pool_get_array(mesh, 'zb3', zb3)

   call iamas_pool_get_array(mesh, 'zgrid', zgrid)
   call iamas_pool_get_array(mesh, 'rdzw', rdzw)
   call iamas_pool_get_array(mesh, 'dzu', dzu)
   call iamas_pool_get_array(mesh, 'rdzu', rdzu)
   call iamas_pool_get_array(mesh, 'fzm', fzm)
   call iamas_pool_get_array(mesh, 'fzp', fzp)
   call iamas_pool_get_array(mesh, 'zxu', zxu)
   call iamas_pool_get_array(mesh, 'zz', zz)
   call iamas_pool_get_array(mesh, 'hx', hx)
   call iamas_pool_get_array(mesh, 'ter', ter)
   call iamas_pool_get_array(mesh, 'dss', dss)

   call iamas_pool_get_array(diag, 'exner_base', exner_b)
   call iamas_pool_get_array(diag, 'rho_base', rho_b)
   call iamas_pool_get_array(diag, 'theta_base', theta_b)
   call iamas_pool_get_array(diag, 'rtheta_base', rtheta_b)
   call iamas_pool_get_array(diag, 'exner', exner)
   call iamas_pool_get_array(diag, 'pressure_base', pressure_b)
   call iamas_pool_get_array(diag, 'pressure_p', pressure_p)
   call iamas_pool_get_array(diag, 'pressure', pressure)
   call iamas_pool_get_array(diag, 'phi_base', phi_b)
   call iamas_pool_get_array(diag, 'phi_total', phi_total)
   call iamas_pool_get_array(diag, 'surface_pressure', surface_pressure)
   call iamas_pool_get_array(diag, 'relhum', relhum)
   call iamas_pool_get_array(diag, 'spechum', spechum)
   call iamas_pool_get_array(diag, 'ru', ru)
   call iamas_pool_get_array(diag, 'rw', rw)

   call iamas_pool_get_array(state, 'rho_zz', rho_zz)
   call iamas_pool_get_array(diag, 'rho_p', rho_p)
   call iamas_pool_get_array(state, 'theta_m', theta_m)
   call iamas_pool_get_array(diag, 'rtheta_p', rtheta_p)
   call iamas_pool_get_array(state, 'scalars', scalars)
   call iamas_pool_get_array(state, 'u', u)
   call iamas_pool_get_array(state, 'w', w)
   call iamas_pool_get_array(diag, 'theta', theta)
   call iamas_pool_get_array(diag, 'rho', rho)
   call iamas_pool_get_array(diag, 'uReconstructX', uReconstructX)
   call iamas_pool_get_array(diag, 'uReconstructY', uReconstructY)
   call iamas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
   call iamas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
   call iamas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

   !jwfeng: For fdda nudging 1218/2025
   call iamas_pool_get_array(fdda, 'qv_fdda_new', qv_fdda_new)
   call iamas_pool_get_array(fdda, 'th_fdda_new', th_fdda_new)
   call iamas_pool_get_array(fdda, 'u_fdda_new', u_fdda_new)
   call iamas_pool_get_array(fdda, 'v_fdda_new', v_fdda_new)

   call iamas_pool_get_array(mesh, 'latCell', latCell)
   call iamas_pool_get_array(mesh, 'lonCell', lonCell)
   call iamas_pool_get_array(mesh, 'latEdge', latEdge)
   call iamas_pool_get_array(mesh, 'lonEdge', lonEdge)
   call iamas_pool_get_array(mesh, 'cf1', cf1)
   call iamas_pool_get_array(mesh, 'cf2', cf2)
   call iamas_pool_get_array(mesh, 'cf3', cf3)
   call iamas_pool_get_array(mesh, 'landmask', landmask)

   call iamas_pool_get_array(fg, 'dzs_fg', dzs_fg)
   call iamas_pool_get_array(fg, 'zs_fg', zs_fg)
   call iamas_pool_get_array(fg, 'sst', sst)
   call iamas_pool_get_array(fg, 'xice', xice)
   call iamas_pool_get_array(fg, 'seaice', seaice)
   call iamas_pool_get_array(fg, 'psfc', psfc)
   call iamas_pool_get_array(fg, 'skintemp', skintemp)
   call iamas_pool_get_array(fg, 'snow', snow)
   call iamas_pool_get_array(fg, 'snowc', snowc)
   call iamas_pool_get_array(fg, 'u', u_fg)
   call iamas_pool_get_array(fg, 'v', v_fg)
   call iamas_pool_get_array(fg, 'z', z_fg)
   call iamas_pool_get_array(fg, 't', t_fg)
   call iamas_pool_get_array(fg, 'rh', rh_fg)
   call iamas_pool_get_array(fg, 'sh', sh_fg)
   call iamas_pool_get_array(fg, 'gfs_z', gfs_z)
   call iamas_pool_get_array(fg, 'p', p_fg)

   call iamas_pool_get_dimension(dims, 'nVertLevels', nz1)
   call iamas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
   nz = nz1 + 1

   call iamas_pool_get_dimension(state, 'index_qv', index_qv)

   xnutr = 0.
   zd = 12000.
   znut = eta_t

   etavs = (1.-0.252)*pii/2.
   rcv = rgas/(cp-rgas)
   omega_e = omega
   p0 = 1.e+05

   scalars(:,:,:) = 0.

   !
   ! If requested, blend the terrain along the domain boundaries with terrain from
   ! an intermediate file. For global domains, this routine will have no effect even
   ! if called, since terrain is only blended for cells with bdyMaskCell > 0.
   !
   if (config_blend_bdy_terrain) then
      call iamas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
      call iamas_pool_get_array(mesh, 'ter', ter)

      call blend_bdy_terrain(config_fdda_prefix, config_start_time, nCells, latCell, lonCell, bdyMaskCell, ter, .false., istatus)
      if (istatus /= 0) then
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('* Blending of terrain along domain boundaries failed!       *', messageType=iAMAS_LOG_ERR)
          call iamas_log_write('*************************************************************', messageType=iAMAS_LOG_CRIT)
      end if
   end if

   if (config_vertical_grid) then

   !
   ! Vertical grid setup
   !
   allocate(hs (nCells+1))
   allocate(hs1(nCells+1))

   ! Fourth order smoother for terrain

   nsmterrain = config_nsmterrain

   do i=1,nsmterrain

      do iCell=1,nCells
         hs(iCell) = 0.
         if(ter(iCell) .ne. 0.) then
            do j = 1,nEdgesOnCell(iCell)

               ! For smoothing at cells along the boundary of the mesh, set the terrain value
               ! for non-existent neighbors, which map to the "garbage cell", to the same as
               ! the terrain in the cell being smoothed
               if (cellsOnCell(j,iCell) == nCells+1) then
                  ter(nCells+1) = ter(iCell)
               end if

               hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                     / dcEdge(edgesOnCell(j,iCell))    &
                                     * (ter(cellsOnCell(j,iCell))-ter(iCell))
            end do
         end if
         hs(iCell) = ter(iCell) + 0.216*hs(iCell)
      end do

      do iCell=1,nCells
         ter(iCell) = 0.
         if(hs(iCell) .ne. 0.) then
            do j = 1,nEdgesOnCell(iCell)

               ! For smoothing at cells along the boundary of the mesh, set the terrain value
               ! for non-existent neighbors, which map to the "garbage cell", to the same as
               ! the terrain in the cell being smoothed
               if (cellsOnCell(j,iCell) == nCells+1) then
                  hs(nCells+1) = hs(iCell)
               end if

               ter(iCell) = ter(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                       / dcEdge(edgesOnCell(j,iCell))    &
                                       * (hs(cellsOnCell(j,iCell))-hs(iCell))
            end do
         end if
         ter(iCell) = hs(iCell) - 0.216*ter(iCell)
      end do

      ! note that ther variable ter used throughout this section is a pointer to grid % ter % array, here we are passing ter's parent field
      call iamas_dmpar_exch_halo_field(ter_field)

   end do

   do iCell=1,nCells
      hx(:,iCell) = ter(iCell)
   end do

   hm = maxval(ter(1:nCellsSolve))
   call iamas_dmpar_max_real(dminfo, hm, hm_global)
   hm = hm_global

   ! Metrics for hybrid coordinate and vertical stretching

   !
   ! If a the name of a file with vertical coordinate values has been specified,
   ! use those values to setup the vertical grid
   !
   if (len_trim(config_specified_zeta_levels) > 0) then

      if (read_text_array(dminfo, trim(config_specified_zeta_levels), specified_zw) /= 0) then
         call iamas_log_write('Failed to read vertical levels from '''//trim(config_specified_zeta_levels)//'''', &
                              messageType=iAMAS_LOG_CRIT)
      end if

      if (size(specified_zw) /= nz) then
         call iamas_log_write('In the namelist.init_atmosphere file, config_nvertlevels = $i, but ', intArgs=(/nz1/), &
                              messageType=iAMAS_LOG_ERR)
         call iamas_log_write('but '''//trim(config_specified_zeta_levels)//''' has $i values.', intArgs=(/size(specified_zw)/), &
                              messageType=iAMAS_LOG_ERR)
         call iamas_log_write(''''//trim(config_specified_zeta_levels)//''' must contain nVertLevels+1 ($i) values.', intArgs=(/nz/), &
                              messageType=iAMAS_LOG_CRIT)
      end if

      zw(:) = specified_zw(:)
      zt = zw(nz)

      deallocate(specified_zw)

   !
   ! Otherwise, see if the user has requested to set up the vertical grid as in the iAMAS TC configuration
   !
   else if (config_tc_vertical_grid) then
      str = 1.5
      zt = config_ztop
      dz = zt/float(nz1)

   !... Laura D. Fowler: change the values als,alt,and zetal valid for nVertLevels=41 to values
   !    needed for nVertLevels equal or greater than 55.
      if (nVertLevels >= 55) then
         als   =  .075  ! 55 levels.
         alt   = 1.70   ! 55 levels.
         zetal =  .75   ! 55 levels.
      else
         als   =  .075  ! 41 levels.
         alt   = 1.23   ! 41 levels.
         zetal =  .31   ! 41 levels.
      endif
      if (nVertLevels /= 55 .and. nVertLevels /= 41) then
         call iamas_log_write('********************************************************************************')
         call iamas_log_write('* Note: als, alt, and zetal have not been tested for nVertLevels = $i', intArgs=(/nVertLevels/))
         call iamas_log_write('*       The values of these parameters have only been tested with 41 and 55 layers')
         call iamas_log_write('********************************************************************************')
      end if
   !... end Laura D. Fowler / 2016-04-12.

      do k=1,nz
         zl = 1.-alt*(1.-zetal)
         if ((real(k-1)/real(nz1)).LT.zetal) then
              zw(k) =              (als*real(K-1)/real(nz1)      &
                  +(3.*(1.-alt)+2.*(alt-als)*zetal)              &
                                  *((K-1)*dz/(zt*zetal))**2      &
                  -(2.*(1.-alt) +  (alt-als)*zetal)              &
                                  *(real(k-1)*dz/(zt*zetal))**3) * zt
         else
              zw(K) = (zl+alt*(real(K-1)/real(nz1)-zetal))*zt
         end if
         if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
      end do

   !
   ! Otherwise, use the vertical level configuration from iAMAS v2.0
   !
   else

      call iamas_log_write('Setting up vertical levels as in iAMAS 2.0 and earlier')

      str = 1.5
      zt = config_ztop
      dz = zt/float(nz1)

      do k=1,nz
         zw(k) = (real(k-1)/real(nz1))**str*zt
         if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
      end do
   end if !config_tc_vertical_grid

   hybrid = .true.

   kz = nz
   if (hybrid) then
      zh = 30000.0

      do k=1,nz
         if (zw(k) < zh) then
            ah(k) = cos(.5*pii*zw(k)/zh)**6
         else
            ah(k) = 0.
            kz = min(kz,k)
         end if
      end do

   else

      do k=1,nz
         ah(k) = 1.-zw(k)/zt
      end do

   end if

   do k=1,nz1
      dzw (k) = zw(k+1)-zw(k)
      rdzw(k) = 1./dzw(k)
      zu(k  ) = .5*(zw(k)+zw(k+1))
   end do
   do k=2,nz1
      dzu (k)  = .5*(dzw(k)+dzw(k-1))
      rdzu(k)  =  1./dzu(k)
      fzp (k)  = .5* dzw(k  )/dzu(k)
      fzm (k)  = .5* dzw(k-1)/dzu(k)
      rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
      rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
   end do

   !**********  how are we storing cf1, cf2 and cf3?

   COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2)
   COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3)
   CF1  = FZP(2) + COF1
   CF2  = FZM(2) - COF1 - COF2
   CF3  = COF2

  
   !     Smoothing algorithm for coordinate surfaces 

   smooth = config_smooth_surfaces

   if (smooth) then

      dzmin = config_dzmin

      allocate(sm0(nCells+1))

      do iCell=1,nCells
         dcsum = 0.0
         do j=1,nEdgesOnCell(iCell)
            dcsum = dcsum + dcEdge(edgesOnCell(j,iCell))
         end do
         dcsum = dcsum / real(nEdgesOnCell(iCell))
         sm0(iCell) = max(0.01_RKIND, 0.125 * min(1.0_RKIND, 3000.0_RKIND/dcsum))
      end do

      do k=2,kz-1
         hx(k,:) = hx(k-1,:)
         dzminf = zw(k)-zw(k-1)

         do i=1,config_nsm + k
            do iCell=1,nCells

               sm = sm0(iCell) * min((3.0*zw(k)/zt)**2.0, 1.0_RKIND)

               hs1(iCell) = 0.
               do j = 1,nEdgesOnCell(iCell)

                  ! For smoothing at cells along the boundary of the mesh, set the hx value
                  ! for non-existent neighbors, which map to the "garbage cell", to the same as
                  ! the hx in the cell being smoothed
                  if (cellsOnCell(j,iCell) == nCells+1) then
                     hx(k,nCells+1) = hx(k,iCell)
                  end if

                  hs1(iCell) = hs1(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                        / dcEdge(edgesOnCell(j,iCell))    &
                                        *  (hx(k,cellsOnCell(j,iCell))-hx(k,iCell))
               end do
               hs(iCell) = hx(k,iCell) + sm*hs1(iCell)

            end do

            tempField => tempFieldTarget
            tempField % block => block
            tempField % dimSizes(1) = nCells
            tempField % sendList => parinfo % cellsToSend
            tempField % recvList => parinfo % cellsToRecv
            tempField % copyList => parinfo % cellsToCopy
            tempField % array => hs
            tempField % isActive = .true.
            tempField % prev => null()
            tempField % next => null()

            call iamas_dmpar_exch_halo_field(tempField)

            do iCell=1,nCells
               dzmina = (zw(k) + ah(k)*hs(iCell)) - (zw(k-1) + ah(k-1)*hx(k-1,iCell))
               if (dzmina > dzmin*(zw(k)-zw(k-1))) then
                  hx(k,iCell) = hs(iCell)
                  if (dzmina < dzminf) dzminf = dzmina
               end if
            end do

         end do
         call iamas_dmpar_min_real(dminfo, dzminf, dzminf_global)
      end do

      deallocate(sm0)

      do k=kz,nz
            hx(k,:) = 0.
      end do

   else

      do k=2,nz1
         dzmina = minval(zw(k)+ah(k)*hx(k,:)-zw(k-1)-ah(k-1)*hx(k-1,:))
      end do

   end if

   deallocate(hs )
   deallocate(hs1)

   !     Height of coordinate levels (calculation of zgrid)

   do iCell=1,nCells
      do k=1,nz
         zgrid(k,iCell) = zw(k) + ah(k)*hx(k,iCell)
      end do
      do k=1,nz1
         zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
      end do
   end do

   do i=1, nEdges
      iCell1 = cellsOnEdge(1,i)
      iCell2 = cellsOnEdge(2,i)
      do k=1,nz1
         zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
      end do
   end do
   do i=1, nCells
      do k=1,nz1
        ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
        dss(k,i) = 0.
        ztemp = zgrid(k,i)
        if (ztemp.gt.zd+.1)  then
            dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
        end if
      end do
   end do


   ! For z-metric term in omega equation
   do iEdge = 1,nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)

      ! Avoid referencing the garbage cell for exterior edges
      if (cell1 == nCells+1) then
         cell1 = cell2
      end if
      if (cell2 == nCells+1) then
         cell2 = cell1
      end if

      if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

         do k = 1, nVertLevels

            if (config_theta_adv_order == 2) then

               z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

            else !theta_adv_order == 3 or 4

               d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
               do i=1, nEdgesOnCell(cell1)
                  if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
               end do
               do i=1, nEdgesOnCell(cell2)
                  if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
               end do

               z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                             - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               if (config_theta_adv_order == 3) then
                  z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
               else
                  z_edge3 = 0.
               end if

            end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

         end do

      end if
   end do

   end if    ! config_vertical_grid


   if (config_met_interp) then

   !ldf (2011-11-19): added initialization of the sea-surface temperature, seaice fraction, and
   !seaice flag:
    sst = 0.0
    xice = 0.0
    seaice = 0.0
   !ldf end.

   !
   ! First, try to locate the LANDSEA field for use as an interpolation mask
   !
   call read_met_init(trim(config_fdda_prefix), .false., timeString(1:13), istatus)

   if (istatus /= 0) then
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
      call iamas_log_write('Error opening fdda file ' &
                                   //trim(config_fdda_prefix)//':'//timeString(1:13),                   messageType=iAMAS_LOG_ERR)
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
   end if

   call read_next_met_field(field, istatus)

   do while (istatus == 0)
      if (trim(field % field) == 'LANDSEA') then

         allocate(maskslab(-2:field % nx+3, field % ny))
         maskslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
         maskslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
         maskslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
         maskslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
         maskslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
         maskslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
         maskslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

      end if

      deallocate(field % slab)
      call read_next_met_field(field, istatus)
   end do

   call read_met_close()

   if (.not. allocated(maskslab)) then
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
      call iamas_log_write('LANDSEA field not found in fdda file ' &
                                   //trim(config_fdda_prefix)//':'//timeString(1:13),                   messageType=iAMAS_LOG_ERR)
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
   end if

   edge_mask(:) = 1


   !
   ! Horizontally interpolate meteorological data
   !
   allocate(vert_level(config_nfglevels))
   vert_level(:) = -1.0

   call read_met_init(trim(config_fdda_prefix), .false., timeString(1:13), istatus)

   if (istatus /= 0) then
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
      call iamas_log_write('Error opening initial fdda file ' &
                                   //trim(config_fdda_prefix)//':'//timeString(1:13),                   messageType=iAMAS_LOG_ERR)
      call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
   else
      call iamas_log_write('    Opening fdda file '//trim(config_fdda_prefix)//':'//timeString(1:13))
   end if

   call iamas_hash_init(level_hash)
   too_many_fg_levs = .false.

   call read_next_met_field(field, istatus)

   do while (istatus == 0)

      interp_list(1) = FOUR_POINT
      interp_list(2) = SEARCH
      interp_list(3) = 0

      maskval = -1.0
      masked = -1
      fillval = 0.0
      msgval = -1.e30

      mask_array => landmask

      if (trim(field % field) == 'UU' .or. &
          trim(field % field) == 'VV' .or. &
          trim(field % field) == 'TT' .or. &
          trim(field % field) == 'RH' .or. &
          trim(field % field) == 'SPECHUMD' .or. &
          trim(field % field) == 'GHT' .or. &
          trim(field % field) == 'PMSL' .or. &
          trim(field % field) == 'PSFC' .or. &
          trim(field % field) == 'SOILHGT' .or. &
          trim(field % field) == 'SM000010' .or. &
          trim(field % field) == 'SM010040' .or. &
          trim(field % field) == 'SM040100' .or. &
          trim(field % field) == 'SM100200' .or. &
          trim(field % field) == 'SM010200' .or. &
          trim(field % field) == 'SM000007' .or. &
          trim(field % field) == 'SM007028' .or. &
          trim(field % field) == 'SM028100' .or. &
          trim(field % field) == 'SM100255' .or. &
          trim(field % field) == 'SM100289' .or. &
          trim(field % field) == 'ST000010' .or. &
          trim(field % field) == 'ST010040' .or. &
          trim(field % field) == 'ST040100' .or. &
          trim(field % field) == 'ST100200' .or. &
          trim(field % field) == 'ST010200' .or. &
          trim(field % field) == 'ST000007' .or. &
          trim(field % field) == 'ST007028' .or. &
          trim(field % field) == 'ST028100' .or. &
          trim(field % field) == 'ST100255' .or. &
          trim(field % field) == 'ST100289' .or. &
          trim(field % field) == 'PRES' .or. &
          trim(field % field) == 'PRESSURE' .or. &
          trim(field % field) == 'SNOW' .or. &
          trim(field % field) == 'SEAICE' .or. &
          trim(field % field) == 'SKINTEMP') then

         if(trim(field % field) == 'SM000010' .or. &
            trim(field % field) == 'SM010040' .or. &
            trim(field % field) == 'SM040100' .or. &
            trim(field % field) == 'SM100200' .or. &
            trim(field % field) == 'SM010200' .or. &
            trim(field % field) == 'SM000007' .or. &
            trim(field % field) == 'SM007028' .or. &
            trim(field % field) == 'SM028100' .or. &
            trim(field % field) == 'SM100255' .or. &
            trim(field % field) == 'SM100289' .or. &
            trim(field % field) == 'ST000010' .or. &
            trim(field % field) == 'ST010040' .or. &
            trim(field % field) == 'ST040100' .or. &
            trim(field % field) == 'ST100200' .or. &
            trim(field % field) == 'ST010200' .or. &
            trim(field % field) == 'ST000007' .or. &
            trim(field % field) == 'ST007028' .or. &
            trim(field % field) == 'ST028100' .or. &
            trim(field % field) == 'ST100255' .or. &
            trim(field % field) == 'ST100289' .or. &
            trim(field % field) == 'SNOW' .or. &
            trim(field % field) == 'SEAICE' .or. &
            trim(field % field) == 'SKINTEMP') then
            k = 1
         else if (trim(field % field) /= 'PMSL' .and. &
                  trim(field % field) /= 'PSFC' .and. &
                  trim(field % field) /= 'SOILHGT') then

            ! Since the hash table can only store integers, transfer the bit pattern from
            ! the real-valued xlvl into an integer; that the result is not an integer version
            ! of the level is not important, since we only want to test uniqueness of levels
            level_value = transfer(field % xlvl, level_value)
            if (.not. iamas_hash_search(level_hash, level_value)) then
               call iamas_hash_insert(level_hash, level_value)
               if (iamas_hash_size(level_hash) > config_nfglevels) then
                  too_many_fg_levs = .true.
               end if
            end if

            ! In case we have more than config_nfglevels levels, just keep cycling through
            ! the remaining fields in the intermediate file for the purpose of counting how
            ! many unique levels are found using the code above
            if (too_many_fg_levs) then
               call read_next_met_field(field, istatus)
               cycle
            end if

            do k=1,config_nfglevels
               if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
            end do
            if (vert_level(k) == -1.0) vert_level(k) = field % xlvl
         else
            k = 1
         end if
         ! Set up projection
         call map_init(proj)

         if (field % iproj == PROJ_LATLON) then
            call map_set(PROJ_LATLON, proj, &
                         latinc = real(field % deltalat,RKIND), &
                         loninc = real(field % deltalon,RKIND), &
                         knowni = 1.0_RKIND, &
                         knownj = 1.0_RKIND, &
                         lat1 = real(field % startlat,RKIND), &
                         lon1 = real(field % startlon,RKIND))
         else if (field % iproj == PROJ_GAUSS) then
            call map_set(PROJ_GAUSS, proj, &
                         nlat = nint(field % deltalat), &
                         loninc = 360.0_RKIND / real(field % nx,RKIND), &
                         lat1 = real(field % startlat,RKIND), &
                         lon1 = real(field % startlon,RKIND))
         end if

         ! Horizontally interpolate the field at level k
         if (trim(field % field) == 'UU') then
            ! For U10, interpolate to cell centers
            if (vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'u10', destField1d)
               ndims = 1
            ! otherwise to edges
            else
               mask_array => edge_mask
               nInterpPoints = nEdges
               latPoints => latEdge
               lonPoints => lonEdge
               call iamas_pool_get_array(fg, 'u', destField2d)
               ndims = 2
            end if
         else if (trim(field % field) == 'VV') then
            ! For V10, interpolate to cell centers
            if (vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call iamas_pool_get_array(fg, 'v10', destField1d)
               ndims = 1
            ! otherwise to edges
            else
               mask_array => edge_mask

               nInterpPoints = nEdges
               latPoints => latEdge
               lonPoints => lonEdge
               call iamas_pool_get_array(fg, 'v', destField2d)
               ndims = 2
            end if

         else if (trim(field % field) == 'TT') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 't', destField2d)
            ndims = 2
         else if (trim(field % field) == 'RH') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'rh', destField2d)
            ndims = 2
         else if (trim(field % field) == 'SPECHUMD') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'sh', destField2d)
            ndims = 2
         else if (trim(field % field) == 'GHT') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'z', destField2d)
            ndims = 2
         else if (trim(field % field) == 'PRES') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'p', destField2d)
            ndims = 2
         else if (trim(field % field) == 'PRESSURE') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'p', destField2d)
            ndims = 2
         else if (trim(field % field) == 'PMSL') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'pmsl', destField1d)
            ndims = 1
         else if (trim(field % field) == 'PSFC') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'psfc', destField1d)
            ndims = 1
         else if (trim(field % field) == 'SOILHGT') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'soilz', destField1d)
            ndims = 1
         else if (trim(field % field) == 'SKINTEMP') then
            nInterpPoints = nCells
            latPoints => latCell
            lonPoints => lonCell
            call iamas_pool_get_array(fg, 'skintemp', destField1d)
            ndims = 1
         end if

         allocate(rslab(-2:field % nx+3, field % ny))
         rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
         rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
         rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
         rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
         rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
         rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
         rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

         do i=1,nInterpPoints
            if (mask_array(i) /= masked) then
               lat = latPoints(i)*DEG_PER_RAD
               lon = lonPoints(i)*DEG_PER_RAD
               call latlon_to_ij(proj, lat, lon, x, y)
               if (x < 0.5) then
                  lon = lon + 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               else if (x >= real(field%nx)+0.5) then
                  lon = lon - 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               end if
               if (y < 0.5) then
                  y = 1.0
               else if (y >= real(field%ny)+0.5) then
                  y = real(field%ny)
               end if
               if (ndims == 1) then
                  if (maskval /= -1.0) then
                     destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                  else
                     destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                  end if
               else if (ndims == 2) then
                  if (maskval /= -1.0) then
                     destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                  else
                     destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                  end if
               end if
            else
               if (ndims == 1) then
                  destField1d(i) = fillval
               else if (ndims == 2) then
                  destField2d(k,i) = fillval
               end if
            end if
         end do

         deallocate(rslab)

      end if

      deallocate(field % slab)
      call read_next_met_field(field, istatus)
   end do

   call read_met_close()
   level_value = iamas_hash_size(level_hash)
   call iamas_hash_destroy(level_hash)

   if (too_many_fg_levs) then
      write(errstring,'(a,i4)') '       Please increase config_nfglevels to at least ', level_value
      call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_ERR)
      call iamas_log_write('Error: The meteorological data file has more than config_nfglevels.', messageType=iAMAS_LOG_ERR)
      call iamas_log_write(trim(errstring),                                                       messageType=iAMAS_LOG_ERR)
      call iamas_log_write('       in the namelist and re-run.',                                  messageType=iAMAS_LOG_ERR)
      call iamas_log_write('*******************************************************************', messageType=iAMAS_LOG_CRIT)
   end if

   !
   ! Check how many distinct levels we actually found in the meteorological data
   !
   do k=1,config_nfglevels
      if (vert_level(k) == -1.0) exit
   end do
   nfglevels_actual = k-1
   call iamas_log_write('    Found $i levels in the first-guess data', intArgs=(/nfglevels_actual/))

   !
   ! Extract surface fields from first-guess
   !
   call iamas_pool_get_array(fg, 'u', u_fg)
   call iamas_pool_get_array(fg, 'v', v_fg)
   call iamas_pool_get_array(fg, 't', t_fg)
   call iamas_pool_get_array(fg, 'rh', rh_fg)
   !
   ! For isobaric data, fill in the 3-d pressure field; otherwise, ensure
   ! that the surface pressure and height fields are filled in
   !
   if (minval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0 .and. &
       maxval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0) then
      call iamas_log_write('    Setting pressure field for isobaric data')
      do k=1,config_nfglevels
         if (vert_level(k) /= 200100.0) then
            p_fg(k,:) = vert_level(k)
         else
            p_fg(k,:) = psfc(:)
         end if
      end do
   else
      call iamas_pool_get_array(fg, 'z', z_fg)
      call iamas_pool_get_array(fg, 'soilz', soilz)
      call iamas_log_write('    Assuming model-level input data')
      do k=1,config_nfglevels
         if (vert_level(k) == 200100.0) then
            p_fg(k,:) = psfc(:)
            z_fg(k,:) = soilz(:)
         end if
      end do
   end if
   ! Compute normal wind component and store in fg % u
   do iEdge=1,nEdges
      do k=1,nfglevels_actual
         u_fg(k,iEdge) = cos(angleEdge(iEdge)) * u_fg(k,iEdge) &
                       + sin(angleEdge(iEdge)) * v_fg(k,iEdge)
      end do
   end do

   ! Vertically interpolate meteorological data
   allocate(sorted_arr(2,nfglevels_actual))

   do iCell=1,nCells
      ! T
      sorted_arr(:,:) = -999.0
      do k = 1, nfglevels_actual
         sorted_arr(1,k) = z_fg(k,iCell)
         if (vert_level(k) == 200100.0) sorted_arr(1,k) = 9999999.0 ! jwfeng: avoid err when ztop is larger than 99999.0
         sorted_arr(2,k) = t_fg(k,iCell)
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k = 1, nVertLevels
         target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
         ! jwfeng: Use temperature instead of theta_m
         temperature(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                      sorted_arr(:,1:nfglevels_actual-1), order=1, &
                                      extrap=extrap_airtemp, ierr=istatus)
         if (istatus /= 0) then
            write(errstring,'(a,i4,a,i10)') 'Error in interpolation of t(k,iCell) for k=', k, ', iCell=', iCell
            call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_ERR)
            call iamas_log_write(trim(errstring),                                                     messageType=iAMAS_LOG_ERR)
            call iamas_log_write('*****************************************************************', messageType=iAMAS_LOG_CRIT)
         end if
      end do

      ! RH
      sorted_arr(:,:) = -999.0
      relhum(:,iCell) = 0._RKIND
      do k = 1, nfglevels_actual
         sorted_arr(1,k) = z_fg(k,iCell)
         if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
         sorted_arr(2,k) = rh_fg(k,iCell)
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k = nVertLevels, 1, -1
         target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
         relhum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                    sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
         if (target_z < z_fg(1,iCell) .and. k < nVertLevels) relhum(k,iCell) = relhum(k+1,iCell)
      end do

      ! SPECHUM: if first-guess values are negative, set those values to zero before
      ! vertical interpolation.
      sorted_arr(:,:) = -999.0
      spechum(:,iCell) = 0._RKIND
      do k = 1, nfglevels_actual
         sorted_arr(1,k) = z_fg(k,iCell)
         if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
         sorted_arr(2,k) = max(0._RKIND,sh_fg(k,iCell))
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k = nVertLevels, 1, -1
         target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
         spechum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                     sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
         if (target_z < z_fg(1,iCell) .and. k < nVertLevels) spechum(k,iCell) = spechum(k+1,iCell)
      end do

      ! GHT
      sorted_arr(:,:) = -999.0
      do k = 1, nfglevels_actual
         sorted_arr(1,k) = z_fg(k,iCell)
         if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
         sorted_arr(2,k) = z_fg(k,iCell)
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k = 1, nVertLevels
         target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
         gfs_z(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                   sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
      end do

      ! PRESSURE
      sorted_arr(:,:) = -999.0
      do k = 1, nfglevels_actual
         sorted_arr(1,k) = z_fg(k,iCell)
         if (vert_level(k) == 200100.0) then
            sorted_arr(1,k) = 99999999.0
            sfc_k = k
         end if
         sorted_arr(2,k) = log(p_fg(k,iCell))
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k = 1, nVertLevels
         target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
         pressure(k,iCell) = exp(vertical_interp(target_z, nfglevels_actual-1, &
                                 sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1))
      end do
   end do


   do iEdge=1,nEdges
      ! U
      sorted_arr(:,:) = -999.0
      do k=1,nfglevels_actual
         sorted_arr(1,k) = 0.5 * (z_fg(k,cellsOnEdge(1,iEdge)) + z_fg(k,cellsOnEdge(2,iEdge)))
         if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999999.0
         sorted_arr(2,k) = u_fg(k,iEdge)
      end do
      call iamas_quicksort(nfglevels_actual, sorted_arr)
      do k=1,nVertLevels
         target_z = 0.25 * (zgrid(k,cellsOnEdge(1,iEdge)) + zgrid(k+1,cellsOnEdge(1,iEdge)) + zgrid(k,cellsOnEdge(2,iEdge)) + zgrid(k+1,cellsOnEdge(2,iEdge)))
         u(k,iEdge) = vertical_interp(target_z, nfglevels_actual-1, sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
      end do
   end do

   ! Reconstruct zonal and meridional winds for diagnostic puposes:
   call iamas_rbf_interp_initialize(mesh)
   call iamas_init_reconstruct(mesh)
   call iamas_reconstruct(mesh, u,                 &
                         uReconstructX,           &
                         uReconstructY,           &
                         uReconstructZ,           &
                         uReconstructZonal,       &
                         uReconstructMeridional   &
                        )
   ! Adjust surface pressure for difference in topography
   do sfc_k=1,nfglevels_actual
      if (vert_level(sfc_k) == 200100.) exit
   end do
   do iCell=1,nCells
      ! We need to extrapolate
         sorted_arr(:,:) = -999.0
         do k=1,nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) then
               sorted_arr(1,k) = 99999999.0
            end if
            sorted_arr(2,k) = log(p_fg(k,iCell))
         end do
         call iamas_quicksort(nfglevels_actual, sorted_arr)
         target_z = zgrid(1,iCell)
         psfc(iCell) = exp(vertical_interp(target_z, nfglevels_actual, sorted_arr, order=1, extrap=1))
   end do

   deallocate(sorted_arr)


   ! Diagnose the water vapor mixing ratios:
   global_sh_min = 0._RKIND
   global_sh_max = 0._RKIND
   if(config_use_spechumd) then
      sh_min = minval(spechum(:,1:nCellsSolve))
      sh_max = maxval(spechum(:,1:nCellsSolve))
      call iamas_dmpar_min_real(dminfo,sh_min,global_sh_min)
      call iamas_dmpar_max_real(dminfo,sh_max,global_sh_max)
   endif
   call iamas_log_write('    config_use_spechumd = $l', logicArgs=(/config_use_spechumd/))
   if(.not. config_use_spechumd .or. (global_sh_min==0._RKIND .and. global_sh_max==0._RKIND)) then
      !--- calculate the saturation mixing ratio and interpolated first-guess relative humidity:
      if (config_use_spechumd) then
         call iamas_log_write('config_use_spechumd=T, but specific humidity was not found in '//trim(config_fdda_prefix)//':'//timeString(1:13), messageType=iAMAS_LOG_WARN)
      end if
      call iamas_log_write('    initializing water vapor mixing ratio using first-guess relative humidity')

      do k = 1, nVertLevels
         do iCell = 1, nCells
            !
            ! Note: the RH field provided by ungrib should always be with respect to liquid water,
            !       hence, we can always call rslf; see the routine fix_gfs_rh in WPS/ungrib/src/rrpr.F .
            !
            rs = rslf(pressure(k,iCell),temperature(k,iCell)) ! Use temperature instead of theta_m, jwfeng
            scalars(index_qv,k,iCell) = 0.01_RKIND*rs*relhum(k,iCell)
            qv(k,iCell) = scalars(index_qv,k,iCell)
          enddo
      enddo
   else
      !--- use the interpolated first-guess specific humidity:
      call iamas_log_write(' *** initializing water vapor mixing ratio using first-guess specific humidity')
      do k = 1, nVertLevels
         do iCell = 1, nCells
            scalars(index_qv,k,iCell) = spechum(k,iCell)/(1._RKIND-spechum(k,iCell))
            qv(k,iCell) = scalars(index_qv,k,iCell)
         enddo
      enddo
   endif

   !
   ! After RH has been used to compute qv (unless config_use_spechumd = T and a valid spechum field
   ! is available), modify the RH field to be with respect to ice for temperatures below freezing.
   ! NB: Here we pass in 1:nCells explicitly, since computations involving the "garbage cell" could
   ! trigger FPEs.
   !
   call convert_relhum_wrt_ice(temperature(:,1:nCells), relhum(:,1:nCells)) ! Use temperature instead of theta_m, jwfeng

   !
   ! Diagnose fields needed in initial conditions file (u, w, rho, theta)
   ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
   !
   do iCell=1,nCells
      do k=1,nVertLevels
         ! PI
         exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)

         ! THETA - can compute this using PI instead
         ! t(k,iCell) = t(k,iCell) / exner(k,iCell)
         ! jwfeng: theta_m of RHS is temperature, will be converted into theta_m here
         theta_m(k,iCell) = temperature(k,iCell) / exner(k,iCell) * (1.0 + 1.61*qv(k,iCell))

         ! RHO_ZZ, couple with vertical metric here, jwfeng
         rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
      end do
   end do

   !
   ! Reference state based on a dry isothermal atmosphere
   !
   do iCell=1,nCells
      do k=1,nz1
         ztemp    = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
         pressure_b(k,iCell) = p0*exp(-gravity*ztemp/(rgas*t0b))     ! pressure_base
         exner_b(k,iCell) = (pressure_b(k,iCell)/p0)**(rgas/cp)      ! exner_base
         rho_b(k,iCell) = pressure_b(k,iCell)/(rgas*t0b)/zz(k,iCell) ! rho_base
         theta_b (k,iCell) = t0b/exner_b(k,iCell)                    ! theta_base
         rtheta_b(k,iCell) = rho_b(k,iCell)*theta_b(k,iCell)         ! rtheta_base
         exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)    ! exner
         pressure_p(k,iCell) = pressure(k,iCell)-pressure_b(k,iCell) ! pressure_p
         rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)           ! rho_p

         ! jwfeng: Use phi = ln(p/p0) as log pressure.
         phi_b(k,iCell) = -gravity*ztemp/(rgas*t0b)
         phi_total(k,iCell) = log(pressure(k,iCell)/p0)
      end do
   end do

   if (config_log_pressure) then
      ! jwfeng: Use phi instead of pressure_p to balance the initial thermodynamic fields.
      ! d(phi)/dz = -g*rho_zz*(1+qv)/pressure = -g*(1+qv)/(rgas*T_m)/zz, T_m = (1+1.61*qv)*T = theta_m * exner
      ! To be simplified, use phi without perturbation. Klemp et al., (2021) Adapting the iAMAS ...
      do iCell=1,nCells
      do k=2,nVertLevels
         it = 0
         p_check = 2.0 * 0.000001
         do while ( (it < 100) .and. (p_check > 0.000001) )
            p_check = phi_total(k,iCell)

            phi_total(k,iCell) = phi_total(k-1,iCell) - gravity*dzu(k)*( &
               fzm(k)*( (1.0+qv(k  ,iCell))/(rgas*theta_m(k  ,iCell)*exner(k  ,iCell))/zz(k  ,iCell) ) + &
               fzp(k)*( (1.0+qv(k-1,iCell))/(rgas*theta_m(k-1,iCell)*exner(k-1,iCell))/zz(k-1,iCell) )   &
            )
            pressure(k,iCell) = p0 * exp(phi_total(k,iCell))
            exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)
            ! jwfeng: pressure is adjusted above, so we need to re-compute theta here 
            theta_m(k,iCell) = temperature(k,iCell) / exner(k,iCell) * (1.0 + 1.61*qv(k,iCell))
            rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
            rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)

            p_check = abs(p_check - phi_total(k,iCell))

            pressure_p(k,iCell) = pressure(k,iCell) - pressure_b(k,iCell)

            it = it + 1
         end do
      end do
      end do
   else
      do iCell=1,nCells
      do k=2,nVertLevels
         it = 0
         p_check = 2.0 * 0.0001
         do while ( (it < 30) .and. (p_check > 0.0001) )

            p_check = pressure_p(k,iCell)
            !  WCS 20130821 - iAMAS hydrostatic relation
            pressure_p(k,iCell) = pressure_p(k-1,iCell) - &
                                     (fzm(k)*rho_p(k,iCell) + fzp(k)*rho_p(k-1,iCell))*gravity*dzu(k) &
                                    -(fzm(k)*rho_zz(k  ,iCell)*qv(k  ,iCell) &
                                    + fzp(k)*rho_zz(k-1,iCell)*qv(k-1,iCell))*gravity*dzu(k)
            pressure(k,iCell) = pressure_p(k,iCell) + pressure_b(k,iCell)
            exner(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)
            !  WCS 20130821 - couple with vertical metric
            rho_zz(k,iCell) = pressure(k,iCell) / rgas / (exner(k,iCell)*theta_m(k,iCell)) / zz(k,iCell)
            rho_p(k,iCell) = rho_zz(k,iCell) - rho_b(k,iCell)

            p_check = abs(p_check - pressure_p(k,iCell))

            it = it + 1
         end do
      end do
      end do
   end if ! jwfeng: config_log_pressure

   deallocate(vert_level)

   ! Compute fdda vars
   do iCell=1,nCells
      do k=1,nVertLevels
         th_fdda_new(k,iCell) = theta_m(k,iCell) / (1.0 + 1.61 * qv(k,iCell))
         qv_fdda_new(k,iCell) = qv(k,iCell)
         u_fdda_new(k,iCell) = uReconstructZonal(k,iCell)
         v_fdda_new(k,iCell) = uReconstructMeridional(k,iCell)
      end do
   end do

   end if    ! config_met_interp

end subroutine interp_fdda_to_iAMAS

real (kind=RKIND) function vertical_interp(target_z, nz, zf, order, extrap, surface_val, sealev_val, ierr)

   implicit none

   real (kind=RKIND), intent(in) :: target_z
   integer, intent(in) :: nz
   real (kind=RKIND), dimension(2,nz), intent(in) :: zf      ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
   integer, intent(in), optional :: order
   integer, intent(in), optional :: extrap                   ! can take values 0 = constant, 1 = linear (default), 2 = lapse-rate
   real (kind=RKIND), intent(in), optional :: surface_val
   real (kind=RKIND), intent(in), optional :: sealev_val
   integer, intent(out), optional :: ierr

   integer :: k, lm, lp
   real (kind=RKIND) :: wm, wp
   real (kind=RKIND) :: slope

   integer :: interp_order, extrap_type
   real (kind=RKIND) :: surface, sealevel

   if (present(ierr)) ierr = 0

   if (present(order)) then
      interp_order = order
   else
      interp_order = 2
   end if

   if (present(extrap)) then
      extrap_type = extrap
   else
      extrap_type = 1
   end if

   if (present(surface_val)) then
      surface = surface_val
   else
      surface = 200100.0
   end if

   if (present(sealev_val)) then
      sealevel = sealev_val
   else
      sealevel = 201300.0
   end if

   !
   ! Extrapolation required
   !
   if (target_z < zf(1,1)) then
      if (extrap_type == 0) then
         vertical_interp = zf(2,1)
      else if (extrap_type == 1) then
         slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
         vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
      else if (extrap_type == 2) then
         vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
      end if
      return
   end if
   if (target_z >= zf(1,nz)) then
      if (extrap_type == 0) then
         vertical_interp = zf(2,nz)
      else if (extrap_type == 1) then
         slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
         vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
      else if (extrap_type == 2) then
          call iamas_log_write('extrap_type == 2 not implemented for target_z >= zf(1,nz)', messageType=iAMAS_LOG_ERR)
          if (present(ierr)) ierr = 1
          return
      end if
      return
   end if


   !
   ! No extrapolation required
   !
   do k=1,nz-1
      if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
         lm = k
         lp = k+1
         wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
         wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
         exit
      end if
   end do

   vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

   return

end function vertical_interp
function read_text_array(dminfo, filename, xarray) result(ierr)

   use iamas_io_units, only : iamas_new_unit, iamas_release_unit
   use iamas_log, only : iamas_log_write
   use iamas_derived_types, only : iAMAS_LOG_ERR

   implicit none

   type (dm_info), intent(in) :: dminfo
   character(len=*), intent(in) :: filename
   real (kind=RKIND), dimension(:), pointer :: xarray

   integer :: ierr

   integer :: i
   integer :: nlines
   integer :: iunit
   integer :: iexists
   logical :: exists
   real (kind=RKIND) :: rtemp


   ierr = 1
   nullify(xarray)

   !
   ! Check whether the file exists
   !
   if (dminfo % my_proc_id == IO_NODE) then
      inquire(file=filename, exist=exists)
      if (exists) then
         iexists = 1
      else
         iexists = 0
      end if
   end if
   call iamas_dmpar_bcast_int(dminfo, iexists)

   if (iexists == 0) then
      call iamas_log_write('Text file '''//filename//''' does not exist.', messageType=iAMAS_LOG_ERR)
      return
   end if

   !
   ! Count the number of lines in the file
   !
   if (dminfo % my_proc_id == IO_NODE) then
      call iamas_new_unit(iunit)
      open(unit=iunit, file=filename, form='formatted', status='old', iostat=ierr)
      if (ierr /= 0) then
         nlines = -1
      else
         nlines = 0
         read(unit=iunit, fmt=*, iostat=ierr) rtemp
         do while (ierr == 0)
            nlines = nlines + 1
            read(unit=iunit, fmt=*, iostat=ierr) rtemp
         end do
      end if
      close(unit=iunit)
      call iamas_release_unit(iunit)
   end if
   call iamas_dmpar_bcast_int(dminfo, nlines)

   if (nlines <= 0) then
      if (nlines < 0) then
         call iamas_log_write('Text file '''//filename//''' could not be opened for reading.', messageType=iAMAS_LOG_ERR)
      else
         call iamas_log_write('Text file '''//filename//''' contains no readable real values.', messageType=iAMAS_LOG_ERR)
      end if
      ierr = 1
      return
   end if

   !
   ! Allocate output array, read, and broadcast
   !
   allocate(xarray(nlines))
   if (dminfo % my_proc_id == IO_NODE) then
      call iamas_new_unit(iunit)
      open(unit=iunit, file=filename, form='formatted', status='old', iostat=ierr)
      do i=1,nlines
         read(unit=iunit, fmt=*, iostat=ierr) xarray(i)
      end do
      close(unit=iunit)
      call iamas_release_unit(iunit)
   end if

   call iamas_dmpar_bcast_reals(dminfo, nlines, xarray)
   ierr = 0

end function read_text_array
subroutine blend_bdy_terrain(config_met_prefix, config_start_time, nCells, latCell, lonCell, bdyMaskCell, ter, dryrun, ierr)

   use init_atm_read_met, only : read_met_init, read_met_close, read_next_met_field, met_data
   use init_atm_llxy, only : map_init, map_set, proj_info, latlon_to_ij, PROJ_LATLON, PROJ_GAUSS, DEG_PER_RAD
   use init_atm_hinterp, only : interp_sequence, FOUR_POINT

   implicit none

   character(len=*), intent(in) :: config_met_prefix
   character(len=*), intent(in) :: config_start_time
   integer, intent(in) :: nCells
   real (kind=RKIND), dimension(:), intent(in) :: latCell  ! These four variables (latCell, lonCell, bdyMaskCell, and ter)
   real (kind=RKIND), dimension(:), intent(in) :: lonCell  ! may actually have more than nCells elements, for example,
   integer, dimension(:), intent(in) :: bdyMaskCell        ! if the arrays include a "garbage cell".
   real (kind=RKIND), dimension(:), intent(inout) :: ter   ! 
   logical, intent(in) :: dryrun
   integer, intent(out) :: ierr

   integer, parameter :: nBdyLayers = 7   ! The number of relaxation layers plus the number of specified layers
   integer, parameter :: nSpecLayers = 2  ! The number of specified layers

   integer :: i
   integer :: istatus
   integer, dimension(2) :: interp_list
   real (kind=RKIND) :: weight
   real (kind=RKIND) :: lat, lon
   real (kind=RKIND) :: x, y
   real (kind=RKIND), allocatable, dimension(:,:) :: rslab
   type (met_data) :: field
   type (proj_info) :: proj


   ierr = 0

   if (.not. dryrun) then
      call iamas_log_write('Blending first-guess terrain field along domain boundary')
   end if

   call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

   if (istatus /= 0) then
      call iamas_log_write('********************************************************************************', &
                          messageType=iAMAS_LOG_ERR)
      call iamas_log_write('Error opening file with terrain field, '//trim(config_met_prefix)//':'//config_start_time(1:13), &
                          messageType=iAMAS_LOG_ERR)
      call iamas_log_write('********************************************************************************', &
                          messageType=iAMAS_LOG_ERR)
      ierr = 1
      return
   end if

   !
   ! Loop over fields in the intermediate file looking for the SOILHGT field
   !
   call read_next_met_field(field, istatus)
   do while (istatus == 0)
      if (trim(field % field) == 'SOILHGT') then

         if (.not. dryrun) then
            interp_list(1) = FOUR_POINT
            interp_list(2) = 0

            !
            ! Set up map projection - currently, only the regular lat-lon projection is handled
            !
            call map_init(proj)
       
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat,RKIND), &
                            loninc = real(field % deltalon,RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            end if

            !
            ! Copy the first-guess terrain field into an array that includes some periodic points
            !
            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            !
            ! For each cell in the iAMAS mesh, perform terrain blending if the cell is a boundary cell
            !
            do i=1,nCells
               if (bdyMaskCell(i) > 0) then
                  lat = latCell(i)*DEG_PER_RAD
                  lon = lonCell(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field%ny)
                  end if

                  !
                  ! Is this a specified cell?
                  !
                  if (bdyMaskCell(i) > (nBdyLayers - nSpecLayers)) then
                     ter(i) = interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, &
                                              1, 1, -1.0E30_RKIND, interp_list, 1)

                  !
                  ! Or a relaxation cell?
                  !
                  else
                     weight = real(bdyMaskCell(i),kind=RKIND) / real(nBdyLayers - nSpecLayers,kind=RKIND)
                     ter(i) = weight * interp_sequence(x, y, 1, rslab, -2, field%nx + 3, 1, field%ny, &
                                                       1, 1, -1.0E30_RKIND, interp_list, 1) &
                            + (1.0_RKIND - weight) * ter(i)

                  end if
               end if
            end do

            deallocate(rslab)
         end if

         !
         ! At this point, we have found and processed the first-guess terrain field, so we can return
         !
         deallocate(field % slab)
         call read_met_close()
         return

      end if

      deallocate(field % slab)
      call read_next_met_field(field, istatus)
   end do

   call read_met_close()

   !
   ! If we have reached this point, no first-guess terrain field was found...
   !
   ierr = 1
   call iamas_log_write('********************************************************************************', &
                       messageType=iAMAS_LOG_ERR)
   call iamas_log_write('SOILHGT field not found in intermediate file ' &
                                  //trim(config_met_prefix)//':'//config_start_time(1:13)                , &
                       messageType=iAMAS_LOG_ERR)
   call iamas_log_write('********************************************************************************', &
                       messageType=iAMAS_LOG_ERR)

end subroutine blend_bdy_terrain

subroutine convert_relhum_wrt_ice(t, relhum)

   implicit none

   real (kind=RKIND), dimension(:,:), intent(in) :: t
   real (kind=RKIND), dimension(:,:), intent(inout) :: relhum

   integer :: iCell, k
   integer :: nVertLevels
   integer :: nCells
   real (kind=RKIND) :: eis, ews, r1


   nVertLevels = size(t, dim=1)
   nCells = size(t, dim=2)

   call iamas_log_write('    Recomputing RH w.r.t. ice below freezing')

   do iCell = 1, nCells
      do k = 1, nVertLevels
         if ( t(k,iCell) <= 273.15_RKIND ) then

            ! use formula in ungrib to reconstruct RH 
            eis = 0.01_RKIND * exp (9.550426_RKIND - (5723.265_RKIND / t(k,iCell)) + (3.53068_RKIND * log(t(k,iCell))) &
                - (0.00728332_RKIND * t(k,iCell)))
            ews = 6.112_RKIND * exp(17.67_RKIND * (t(k,iCell)-273.15_RKIND) / ((t(k,iCell)-273.15_RKIND)+243.5_RKIND))

            ! A linear approximation to the GFS blending region ( -20 C > T < 0 C )
            if ( t(k,iCell) > 253.15_RKIND ) then
               r1 = ((273.15_RKIND - t(k,iCell)) / 20.0_RKIND)
               r1 = (r1 * eis) + ((1.0_RKIND-r1)*ews)
            else
               r1 = eis
            end if
            r1 = max(r1, 1.0e-12_RKIND)
            ews = max(ews, 0.0_RKIND)
            relhum(k,iCell) = ews / r1 * relhum(k,iCell)
            relhum(k,iCell) = min(relhum(k,iCell), 100.0_RKIND)
            relhum(k,iCell) = max(relhum(k,iCell), 0.0_RKIND)
         end if
      end do
   end do

end subroutine convert_relhum_wrt_ice

!==================================================================================================
end module iamas_init_atm_fdda
!==================================================================================================

