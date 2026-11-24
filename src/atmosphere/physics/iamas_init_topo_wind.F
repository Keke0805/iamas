!qkyang added to support topo_wind=1,2
!Date: September 25, 2025
!Reference code: https://github.com/wrf-model/WRF/blob/master/dyn_em/start_em.F
!Topo Wind Intercomparsion: https://dtcenter.ucar.edu/eval/meso_mod/topo_wind/index.php
!=================================================================================================================
module iamas_init_topo_wind

   use iamas_kind_types
   use iamas_pool_routines
   use iamas_log

   implicit none
   private
   public :: init_topo_wind_fields

   contains

!=================================================================================================================
   subroutine init_topo_wind_fields(configs, mesh, sfc_input, diag_physics) 
!=================================================================================================================

   implicit none

   type(iamas_pool_type), intent(in) :: configs
   type(iamas_pool_type), intent(in) :: mesh
   type(iamas_pool_type), intent(inout) :: sfc_input
   type(iamas_pool_type), intent(inout) :: diag_physics

   ! Local pointers
   integer, pointer :: nCells, nCellsSolve, iswater
   integer, pointer :: config_topo_wind
   integer, dimension(:), pointer :: nEdgesOnCell
   integer, dimension(:,:), pointer :: cellsOnCell
   real(kind=RKIND), dimension(:), pointer :: ter
   real(kind=RKIND), dimension(:), pointer :: var2d
   real(kind=RKIND), dimension(:), pointer :: ctopo, ctopo2
   real(kind=RKIND), dimension(:), pointer :: lap_hgt

   ! Local variables
   integer :: iCell, iNeighbor, neighborIdx
   real(kind=RKIND) :: alpha, vfac
   real(kind=RKIND) :: ter_sum, ter_center
   integer :: valid_neighbors
   integer, dimension(:), pointer :: ivgtyp

   ! Get configuration
   call iamas_pool_get_config(configs, 'config_topo_wind', config_topo_wind)

   ! Get mesh information
   call iamas_pool_get_dimension(mesh, 'nCells', nCells)
   call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
   call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
   call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
   call iamas_pool_get_array(sfc_input, 'var2d', var2d)!"variance of orography" "m^{2}"
   call iamas_pool_get_array(sfc_input, 'ter', ter)

   ! Get surface input arrays
   call iamas_pool_get_array(sfc_input, 'ivgtyp', ivgtyp)
   call iamas_pool_get_array(sfc_input, 'iswater', iswater) 

   ! Get or allocate ctopo and ctopo2 arrays
   call iamas_pool_get_array(diag_physics, 'ctopo', ctopo)
   call iamas_pool_get_array(diag_physics, 'ctopo2', ctopo2)

   call iamas_log_write('Initializing topo_wind fields...')

   ! Initialize arrays
   ctopo(:) = 1.0_RKIND
   ctopo2(:) = 1.0_RKIND

   ! Calculate ctopo and ctopo2 based on topo_wind option
   if (config_topo_wind == 1) then
      call iamas_log_write('Computing ctopo and ctopo2 for topo_wind = 1')

      ! Allocate lap_hgt for laplacian calculation
      allocate(lap_hgt(nCells))
      lap_hgt(:) = 0.0_RKIND
   
      ! Calculate laplacian of height (lap_hgt) for each cell
      do iCell = 1, nCellsSolve
         ter_sum = 0.0_RKIND
         ter_center = ter(iCell)
         valid_neighbors = 0
   
         ! Sum terrain heights of neighboring cells
         do iNeighbor = 1, nEdgesOnCell(iCell)
            neighborIdx = cellsOnCell(iNeighbor, iCell)
            if (neighborIdx > 0 .and. neighborIdx <= nCells) then
               ter_sum = ter_sum + ter(neighborIdx)
               valid_neighbors = valid_neighbors + 1
            end if
         end do
   
         ! Calculate laplacian: average of neighbors minus center
         if (valid_neighbors > 0) then
            lap_hgt(iCell) = ter_sum / real(valid_neighbors, RKIND) - ter_center
         else
            lap_hgt(iCell) = 0.0_RKIND
         end if
      end do

      do iCell = 1, nCellsSolve
         ! Only apply over land (xland < 1.5)
        !  if (xland(iCell) < 1.5_RKIND) then
         if (ivgtyp(iCell) /= iswater) then
            ctopo(iCell) = sqrt(var2d(iCell))
            
            if (ctopo(iCell) <= 2.718_RKIND) then
               ctopo(iCell) = 1.0_RKIND
            else
               ctopo(iCell) = log(ctopo(iCell))
            end if

            ! Apply laplacian height correction
            if (lap_hgt(iCell) > -10.0_RKIND) then
               ! Keep ctopo as is
               !  continue
               cycle
            else if (lap_hgt(iCell) >= -20.0_RKIND) then
               alpha = (lap_hgt(iCell) + 20.0_RKIND) / 10.0_RKIND
               ctopo(iCell) = alpha * ctopo(iCell) + (1.0_RKIND - alpha)
            else if (lap_hgt(iCell) >= -30.0_RKIND) then
               ctopo(iCell) = (lap_hgt(iCell) + 30.0_RKIND) / 10.0_RKIND
               ctopo2(iCell) = (lap_hgt(iCell) + 30.0_RKIND) / 10.0_RKIND
            else
               ctopo(iCell) = 0.0_RKIND
               ctopo2(iCell) = 0.0_RKIND
            end if
         end if
      end do

      ! Clean up
      deallocate(lap_hgt)
      
   else if (config_topo_wind == 2) then
      call iamas_log_write('Computing ctopo and ctopo2 for topo_wind = 2')
      
      do iCell = 1, nCellsSolve
         if (ivgtyp(iCell) /= iswater) then
            ! Calculate "STANDARD DEVIATION OF SUBGRID-SCALE OROGRAPHY"   "m"  based on sqrt(var2d(iCell))
            vfac = min(1.575_RKIND, (sqrt(var2d(iCell)) * 0.4_RKIND / 200.0_RKIND + 1.175_RKIND))
            vfac = vfac * vfac
         else
            ! Over water, leave it alone
            vfac = 1.0_RKIND
         end if
         
         ctopo(iCell) = ctopo(iCell) * vfac
      end do

   end if



   call iamas_log_write('Topo_wind fields initialization completed.')

   end subroutine init_topo_wind_fields

!=================================================================================================================
end module iamas_init_topo_wind
