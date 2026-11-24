!==================================================================================================
! Module: AddSlopeShadowParam
! Purpose: Precompute terrain slope and directional horizon fields for slope/shadow radiation
! Author: qkyang
! Date: 2026-04-09
!==================================================================================================
module AddSlopeShadowParam_mod

use iamas_derived_types, only : RKIND, iamas_pool_type
use iamas_pool_routines
use iamas_log, only : iamas_log_write
implicit none

public :: AddSlopeShadowParam
public :: AddSlopeShadowParam_fromPools

real(kind=RKIND), parameter :: PI = 3.141592653589793_RKIND

contains

!==================================================================================================
subroutine AddSlopeShadowParam_fromPools(mesh, dims, configs)
!==================================================================================================
! Wrapper for iAMAS pool-style call site:
! keep call signature compact in driver code and resolve concrete arrays here.
implicit none

type (iamas_pool_type), intent(inout) :: mesh
type (iamas_pool_type), intent(in) :: dims
type (iamas_pool_type), intent(in) :: configs

integer, pointer :: config_slope_rad, config_topo_shading
real(kind=RKIND), pointer :: config_shadlen
integer, pointer :: nCells, nCellsSolve, maxEdges, nShadowDirBins
integer, dimension(:), pointer :: nEdgesOnCell
integer, dimension(:,:), pointer :: cellsOnCell
real(kind=RKIND), dimension(:), pointer :: xCell, yCell, ter
real(kind=RKIND), dimension(:), pointer :: ter_slope, ter_hx, ter_hy
real(kind=RKIND), dimension(:,:), pointer :: ter_horizon_angle

call iamas_pool_get_config(configs, 'config_slope_rad', config_slope_rad)
call iamas_pool_get_config(configs, 'config_topo_shading', config_topo_shading)
if (.not. (config_slope_rad == 1 .or. config_topo_shading == 1)) return

call iamas_pool_get_config(configs, 'config_shadlen', config_shadlen)
call iamas_pool_get_dimension(dims, 'nCells', nCells)
call iamas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
call iamas_pool_get_dimension(dims, 'maxEdges', maxEdges)
call iamas_pool_get_dimension(dims, 'nShadowDirBins', nShadowDirBins)

call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
call iamas_pool_get_array(mesh, 'xCell', xCell)
call iamas_pool_get_array(mesh, 'yCell', yCell)
call iamas_pool_get_array(mesh, 'ter', ter)
call iamas_pool_get_array(mesh, 'ter_slope', ter_slope)
call iamas_pool_get_array(mesh, 'ter_hx', ter_hx)
call iamas_pool_get_array(mesh, 'ter_hy', ter_hy)
call iamas_pool_get_array(mesh, 'ter_horizon_angle', ter_horizon_angle)

call iamas_log_write('Start call AddSlopeShadowParam>>>>>>>>>')

! Call the subroutine
call AddSlopeShadowParam(nCells, nCellsSolve, maxEdges, nEdgesOnCell, cellsOnCell, &
                         xCell, yCell, ter, nShadowDirBins, config_shadlen, &
                         ter_slope, ter_hx, ter_hy, ter_horizon_angle)
call iamas_log_write('AddSlopeShadowParam done<<<<<<<<<')

end subroutine AddSlopeShadowParam_fromPools

!==================================================================================================
subroutine AddSlopeShadowParam(nCells, nCellsSolve, maxEdges, nEdgesOnCell, cellsOnCell, &
                               xCell, yCell, ter, nShadowDirBins, shadlen, &
                               ter_slope, ter_hx, ter_hy, ter_horizon_angle)
!==================================================================================================
!
! This routine precomputes static terrain quantities:
! 1) ter_hx, ter_hy and ter_slope from local least-squares terrain plane fit
! 2) min_shadow_dist as local graph step scale
! 3) ter_horizon_angle(dir,cell) by directional neighbor-ring ray marching up to shadlen
!
! Notes on parallel behavior:
! - The directional search is performed as iterative ring expansion on the mesh graph.
! - In decomposed runs, the effective search range depends on available halo depth in nCells.
! - This routine is written to support iterative extension in future communication hooks.
!
implicit none

integer, intent(in) :: nCells, nCellsSolve, maxEdges, nShadowDirBins
integer, dimension(nCells), intent(in) :: nEdgesOnCell
integer, dimension(maxEdges, nCells), intent(in) :: cellsOnCell
real(kind=RKIND), dimension(nCells), intent(in) :: xCell, yCell, ter
real(kind=RKIND), intent(in) :: shadlen

real(kind=RKIND), dimension(nCells), intent(out) :: ter_slope, ter_hx, ter_hy
real(kind=RKIND), dimension(nShadowDirBins, nCells), intent(out) :: ter_horizon_angle

real(kind=RKIND), dimension(nCells) :: min_shadow_dist
integer :: iCell, m, nb, nnbr, iBin, step, max_steps, ibest
real(kind=RKIND) :: sxx, sxy, syy, sxz, syz, det
real(kind=RKIND) :: dx, dy, dz, grad_mag
real(kind=RKIND) :: step_dist, dist_walk, dist_from_src, min_dist_local
real(kind=RKIND) :: bin_angle, dirx, diry, proj, best_proj
real(kind=RKIND) :: topoelev, horizon_max

ter_slope(:) = 0.0_RKIND
ter_hx(:) = 0.0_RKIND
ter_hy(:) = 0.0_RKIND
min_shadow_dist(:) = 0.0_RKIND ! minimum neighboring cell-center distance used as local shadow step scale
!min_shadow_dist approximately equal to the grid spacing
ter_horizon_angle(:,:) = 0.0_RKIND

if (nShadowDirBins <= 0) return

bin_angle = 2.0_RKIND * PI / real(nShadowDirBins, RKIND)

do iCell = 1, nCellsSolve

   nnbr = nEdgesOnCell(iCell)
   if (nnbr <= 0) cycle

   ! Build local least-squares plane:
   ! ter(x,y) = c0 + ter_hx*dx + ter_hy*dy
   sxx = 0.0_RKIND
   sxy = 0.0_RKIND
   syy = 0.0_RKIND
   sxz = 0.0_RKIND
   syz = 0.0_RKIND
   min_dist_local = huge(1.0_RKIND)

   do m = 1, nnbr
      nb = cellsOnCell(m, iCell)
      if (nb <= 0 .or. nb > nCells) cycle

      dx = xCell(nb) - xCell(iCell)
      dy = yCell(nb) - yCell(iCell)
      dz = ter(nb) - ter(iCell)
      step_dist = sqrt(dx*dx + dy*dy)

      if (step_dist > 1.0e-6_RKIND) min_dist_local = min(min_dist_local, step_dist)

      sxx = sxx + dx*dx
      sxy = sxy + dx*dy
      syy = syy + dy*dy
      sxz = sxz + dx*dz
      syz = syz + dy*dz
   end do

   if (min_dist_local < huge(1.0_RKIND)) then
      min_shadow_dist(iCell) = min_dist_local
   else
      min_shadow_dist(iCell) = 0.0_RKIND
   end if

   det = sxx*syy - sxy*sxy
   if (abs(det) > 1.0e-18_RKIND) then
      ter_hx(iCell) = (sxz*syy - syz*sxy) / det
      ter_hy(iCell) = (syz*sxx - sxz*sxy) / det
      grad_mag = sqrt(ter_hx(iCell)*ter_hx(iCell) + ter_hy(iCell)*ter_hy(iCell))
      ter_slope(iCell) = atan(grad_mag)
   end if

   if (shadlen <= 0.0_RKIND) cycle
   if (min_shadow_dist(iCell) <= 1.0e-6_RKIND) cycle

   max_steps = max(2, min(500, int(shadlen / min_shadow_dist(iCell)) + 2))

   do iBin = 1, nShadowDirBins
      ! Direction bin center angle in radians.
      dirx = sin((real(iBin, RKIND) - 0.5_RKIND) * bin_angle)
      diry = cos((real(iBin, RKIND) - 0.5_RKIND) * bin_angle)

      horizon_max = -0.5_RKIND * PI
      dist_walk = 0.0_RKIND
      ibest = iCell

      ! Directional ring expansion along graph-connected neighbors.
      do step = 1, max_steps
         best_proj = 0.05_RKIND
         nb = -1
         step_dist = 0.0_RKIND

         do m = 1, nEdgesOnCell(ibest)
            if (cellsOnCell(m, ibest) <= 0 .or. cellsOnCell(m, ibest) > nCells) cycle
            dx = xCell(cellsOnCell(m, ibest)) - xCell(ibest)
            dy = yCell(cellsOnCell(m, ibest)) - yCell(ibest)
            dist_from_src = sqrt(dx*dx + dy*dy)
            if (dist_from_src <= 1.0e-6_RKIND) cycle
            proj = (dx*dirx + dy*diry) / dist_from_src
            if (proj > best_proj) then
               best_proj = proj
               nb = cellsOnCell(m, ibest)
               step_dist = dist_from_src
            end if
         end do

         if (nb <= 0) exit

         dist_walk = dist_walk + step_dist
         if (dist_walk > shadlen) exit

         dx = xCell(nb) - xCell(iCell)
         dy = yCell(nb) - yCell(iCell)
         dist_from_src = sqrt(dx*dx + dy*dy)
         if (dist_from_src > 1.0e-6_RKIND) then
            topoelev = atan((ter(nb) - ter(iCell)) / dist_from_src)
            if (topoelev > horizon_max) horizon_max = topoelev
         end if

         ibest = nb
      end do

      ! Negative horizon angle is set to zero for shading test stability.
      ter_horizon_angle(iBin, iCell) = max(0.0_RKIND, horizon_max)
   end do
end do

end subroutine AddSlopeShadowParam

end module AddSlopeShadowParam_mod
