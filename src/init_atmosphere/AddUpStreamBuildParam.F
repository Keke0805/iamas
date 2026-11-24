!==================================================================================================
! Module: AddUpStreamBuildParam
! Purpose: Calculate upstream building parameters for each grid cell
! Author: qkyang
! Date: September 22, 2025
!==================================================================================================
module AddUpStreamBuildParam_mod

use iamas_derived_types, only : RKIND, StrKIND, iamas_pool_type
use iamas_pool_routines
use iamas_log, only : iamas_log_write
implicit none

! Make the subroutine public
public :: AddUpStreamBuildParam
public :: AddUpStreamBuildParam_fromPools

! Building parameter constants (corresponding to Python USBP dictionary)
integer, parameter :: N_LCZ_TYPES = 11
integer, parameter :: IVGTYP_MIN = 51
integer, parameter :: IVGTYP_MAX = 61
integer, parameter :: N_SECTORS = 8
! If it is not an LCZ, the parameters of the 5th LCZ type are used by default.
integer, parameter :: NonLCZ_TYPE = 5

! Building parameters for LCZ types 1-11 (indices 51-61)
real(kind=RKIND), parameter, dimension(N_LCZ_TYPES) :: &
    ZR_PARAMS    = (/37.5, 17.5, 6.5, 37.5, 17.5, 6.5, 3.0, 6.5, 6.5, 10.0, 10.0/), &
    ZR_STD_PARAMS= (/4.0,  3.0,  1.0, 1.0,  1.0,  1.0, 1.0, 1.0, 1.0, 1.0,  1.0/), &
    RF_PARAMS    = (/22.2, 22.0, 9.6, 42.86,26.25,13.0,25.0,28.9,43.33,23.8, 5.0/), &
    RD_PARAMS    = (/20.0, 14.0, 5.2, 50.0, 35.0, 13.0,3.33,32.5,43.3, 28.6, 100.0/)

real(kind=RKIND), parameter :: PI = 3.141592653589793_RKIND
real(kind=RKIND), parameter :: DEG_TO_RAD = PI / 180.0_RKIND
real(kind=RKIND), parameter :: RAD_TO_DEG = 180.0_RKIND / PI

contains

!==================================================================================================
subroutine AddUpStreamBuildParam_fromPools(mesh, dims, configs)
!==================================================================================================
! Wrapper for iAMAS pool-style call site:
! keep call signature compact in driver code and resolve concrete arrays here.
implicit none

type (iamas_pool_type), intent(inout) :: mesh
type (iamas_pool_type), intent(in) :: dims
type (iamas_pool_type), intent(in) :: configs

logical, pointer :: config_build_turb_n, config_build_turb_d
logical, pointer :: config_use_lcz
character(len=StrKIND), pointer :: config_landuse_data

integer, pointer :: nCells, maxEdges, nUpstreamSect
integer, dimension(:), pointer :: lu_index, nEdgesOnCell
integer, dimension(:,:), pointer :: cellsOnCell
real(kind=RKIND), dimension(:), pointer :: latCell, lonCell, UpstreamAngles
real(kind=RKIND), dimension(:,:), pointer :: Upstream_zr, Upstream_zr_std, Upstream_rf, Upstream_rd
real(kind=RKIND), dimension(:,:), pointer :: Upstream_zd, Upstream_RSL, Upstream_Z0C

call iamas_pool_get_config(configs, 'config_build_turb_n', config_build_turb_n)
call iamas_pool_get_config(configs, 'config_build_turb_d', config_build_turb_d)
if (.not. (config_build_turb_n .or. config_build_turb_d)) return

call iamas_pool_get_config(configs, 'config_use_lcz', config_use_lcz)
call iamas_pool_get_config(configs, 'config_landuse_data', config_landuse_data)
call iamas_pool_get_dimension(dims, 'nCells', nCells)
call iamas_pool_get_dimension(dims, 'maxEdges', maxEdges)
call iamas_pool_get_dimension(dims, 'nUpstreamSect', nUpstreamSect)

call iamas_pool_get_array(mesh, 'lu_index', lu_index)
call iamas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
call iamas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
call iamas_pool_get_array(mesh, 'latCell', latCell)
call iamas_pool_get_array(mesh, 'lonCell', lonCell)
call iamas_pool_get_array(mesh, 'UpstreamAngles', UpstreamAngles)
call iamas_pool_get_array(mesh, 'Upstream_zr', Upstream_zr)
call iamas_pool_get_array(mesh, 'Upstream_zr_std', Upstream_zr_std)
call iamas_pool_get_array(mesh, 'Upstream_rf', Upstream_rf)
call iamas_pool_get_array(mesh, 'Upstream_rd', Upstream_rd)
call iamas_pool_get_array(mesh, 'Upstream_zd', Upstream_zd)
call iamas_pool_get_array(mesh, 'Upstream_RSL', Upstream_RSL)
call iamas_pool_get_array(mesh, 'Upstream_Z0C', Upstream_Z0C)

call iamas_log_write('Start call AddUpStreamBuildParam>>>>>>>>>')

! Call the subroutine
call AddUpStreamBuildParam(lu_index, config_use_lcz, config_landuse_data, &
                           nCells, maxEdges, nEdgesOnCell, cellsOnCell, &
                           latCell, lonCell, nUpstreamSect, &
                           UpstreamAngles, Upstream_zr, Upstream_zr_std, &
                           Upstream_rf, Upstream_rd, &
                           Upstream_zd, Upstream_RSL, Upstream_Z0C)

call iamas_log_write('AddUpStreamBuildParam done<<<<<<<<<')

end subroutine AddUpStreamBuildParam_fromPools

!==================================================================================================
subroutine AddUpStreamBuildParam(lu_index, config_use_lcz, config_landuse_data, &
                                nCells, maxEdges, nEdgesOnCell, cellsOnCell, &
                                latCell, lonCell, nUpstreamSect, &
                                UpstreamAngles, Upstream_zr, Upstream_zr_std, &
                                Upstream_rf, Upstream_rd, &
                                Upstream_zd, Upstream_RSL, Upstream_Z0C)
!==================================================================================================

implicit none

! Input arguments
integer, dimension(nCells), intent(in) :: lu_index
logical, intent(in) :: config_use_lcz
character(len=*), intent(in) :: config_landuse_data
integer, intent(in) :: nCells, maxEdges
integer, dimension(nCells), intent(in) :: nEdgesOnCell
integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell  
real(kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
integer, intent(in) :: nUpstreamSect

! Output arguments
real(kind=RKIND), dimension(nUpstreamSect), intent(out) :: UpstreamAngles
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_zr
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_zr_std
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_rf
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_rd
! New output arguments for pre-calculated parameters
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_zd
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_RSL
real(kind=RKIND), dimension(nUpstreamSect, nCells), intent(out) :: Upstream_Z0C

! Local variables
integer :: iCell, iSector, iNeighbor, neighborIdx
integer :: sector_counts(nUpstreamSect)
real(kind=RKIND) :: sector_angle, bearing
real(kind=RKIND) :: current_zr, current_zr_std, current_rf, current_rd
real(kind=RKIND) :: neighbor_zr, neighbor_zr_std, neighbor_rf, neighbor_rd
real(kind=RKIND) :: sector_sums_zr(nUpstreamSect), sector_sums_zr_std(nUpstreamSect)
real(kind=RKIND) :: sector_sums_rf(nUpstreamSect), sector_sums_rd(nUpstreamSect)
logical :: sector_has_data(nUpstreamSect)

! Parameters for calculating zd, RSL and Z0C (following module_sf_urban.F)
real(kind=RKIND), parameter :: VonK = 0.4_RKIND
real(kind=RKIND), parameter :: Cd = 1.2_RKIND
real(kind=RKIND), parameter :: alpha_macd = 4.43_RKIND
real(kind=RKIND), parameter :: beta_macd = 1.0_RKIND

! Initialize sector angles (only once, not per cell)
sector_angle = 360.0_RKIND / real(nUpstreamSect, RKIND)
do iSector = 1, nUpstreamSect
    UpstreamAngles(iSector) = real(iSector-1, RKIND) * sector_angle + sector_angle * 0.5_RKIND
end do

! Process each cell
do iCell = 1, nCells
    
    ! Get current cell building parameters
    call get_building_params(lu_index(iCell), config_use_lcz, config_landuse_data, &
                           current_zr, current_zr_std, current_rf, current_rd)
    
    ! Initialize sector accumulators
    sector_counts(:) = 0
    sector_sums_zr(:) = 0.0_RKIND
    sector_sums_zr_std(:) = 0.0_RKIND
    sector_sums_rf(:) = 0.0_RKIND
    sector_sums_rd(:) = 0.0_RKIND
    sector_has_data(:) = .false.
    
    ! Process neighbors
    do iNeighbor = 1, nEdgesOnCell(iCell)
        neighborIdx = cellsOnCell(iNeighbor, iCell)
        
        if (neighborIdx > 0 .and. neighborIdx <= nCells) then
            ! Calculate bearing from current cell to neighbor
            bearing = calculate_bearing(latCell(iCell), lonCell(iCell), &
                                      latCell(neighborIdx), lonCell(neighborIdx))
            
            ! Determine sector
            iSector = get_upstream_sector(bearing, nUpstreamSect)
            
            ! Get neighbor building parameters
            call get_building_params(lu_index(neighborIdx), config_use_lcz, config_landuse_data, &
                                   neighbor_zr, neighbor_zr_std, neighbor_rf, neighbor_rd)
            
            ! Accumulate in sector
            sector_counts(iSector) = sector_counts(iSector) + 1
            sector_sums_zr(iSector) = sector_sums_zr(iSector) + neighbor_zr
            sector_sums_zr_std(iSector) = sector_sums_zr_std(iSector) + neighbor_zr_std
            sector_sums_rf(iSector) = sector_sums_rf(iSector) + neighbor_rf
            sector_sums_rd(iSector) = sector_sums_rd(iSector) + neighbor_rd
            sector_has_data(iSector) = .true.
        end if
    end do
    
    ! Calculate sector averages and interpolate missing sectors
    call process_sectors(nUpstreamSect, sector_counts, sector_has_data, &
                        sector_sums_zr, sector_sums_zr_std, sector_sums_rf, sector_sums_rd, &
                        current_zr, current_zr_std, current_rf, current_rd, &
                        Upstream_zr(:,iCell), Upstream_zr_std(:,iCell), &
                        Upstream_rf(:,iCell), Upstream_rd(:,iCell))
    
    ! Calculate derived parameters for each sector
    do iSector = 1, nUpstreamSect
        call calculate_derived_params(Upstream_zr(iSector,iCell), &
                                    Upstream_rf(iSector,iCell), &
                                    Upstream_rd(iSector,iCell), &
                                    Upstream_zd(iSector,iCell), &
                                    Upstream_RSL(iSector,iCell), &
                                    Upstream_Z0C(iSector,iCell))
    end do
    
end do

end subroutine AddUpStreamBuildParam

!==================================================================================================
subroutine get_building_params(ivgtyp_value, config_use_lcz, config_landuse_data, &
                              zr, zr_std, rf, rd)
!==================================================================================================
implicit none

integer, intent(in) :: ivgtyp_value
logical, intent(in) :: config_use_lcz
character(len=*), intent(in) :: config_landuse_data
real(kind=RKIND), intent(out) :: zr, zr_std, rf, rd

integer :: idx

if (config_use_lcz) then
    ! Use LCZ classification (11 types, indices 51-61)
    if (ivgtyp_value >= IVGTYP_MIN .and. ivgtyp_value <= IVGTYP_MAX) then
        idx = ivgtyp_value - IVGTYP_MIN + 1
        zr = ZR_PARAMS(idx)
        zr_std = ZR_STD_PARAMS(idx)
        rf = RF_PARAMS(idx)
        rd = RD_PARAMS(idx)
    else
        zr = 0.0_RKIND
        zr_std = 0.0_RKIND
        rf = 0.0_RKIND
        rd = 0.0_RKIND
    end if
else
    ! Use single urban category based on landuse data
    select case(trim(config_landuse_data))
    case('USGS')
        ! Urban index = 1, use LCZ type 5 parameters (index NonLCZ_TYPE in array, corresponds to ivgtyp=55)
        if (ivgtyp_value == 1) then
            zr = ZR_PARAMS(NonLCZ_TYPE)      ! 17.5
            zr_std = ZR_STD_PARAMS(NonLCZ_TYPE)  ! 1.0
            rf = RF_PARAMS(NonLCZ_TYPE)      ! 26.25
            rd = RD_PARAMS(NonLCZ_TYPE)      ! 35.0
        else
            zr = 0.0_RKIND
            zr_std = 0.0_RKIND
            rf = 0.0_RKIND
            rd = 0.0_RKIND
        end if
    case('MODIFIED_IGBP_MODIS_NOAH', 'MODIFIED_IGBP_MODIS_NOAH_15s')
        ! Urban index = 13, use LCZ type 5 parameters (index NonLCZ_TYPE in array, corresponds to ivgtyp=55)
        if (ivgtyp_value == 13) then
            zr = ZR_PARAMS(NonLCZ_TYPE)      ! 17.5
            zr_std = ZR_STD_PARAMS(NonLCZ_TYPE)  ! 1.0
            rf = RF_PARAMS(NonLCZ_TYPE)      ! 26.25
            rd = RD_PARAMS(NonLCZ_TYPE)      ! 35.0
        else
            zr = 0.0_RKIND
            zr_std = 0.0_RKIND
            rf = 0.0_RKIND
            rd = 0.0_RKIND
        end if
    case default
        ! Default case: no urban parameters
        zr = 0.0_RKIND
        zr_std = 0.0_RKIND
        rf = 0.0_RKIND
        rd = 0.0_RKIND
    end select
end if

end subroutine get_building_params

!==================================================================================================
function calculate_bearing(lat1, lon1, lat2, lon2) result(bearing_deg)
!==================================================================================================
implicit none

real(kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2
real(kind=RKIND) :: bearing_deg

real(kind=RKIND) :: dlon, y, x, bearing_rad

! Convert to radians for calculation
dlon = (lon2 - lon1) * DEG_TO_RAD
y = sin(dlon) * cos(lat2 * DEG_TO_RAD)
x = cos(lat1 * DEG_TO_RAD) * sin(lat2 * DEG_TO_RAD) - &
    sin(lat1 * DEG_TO_RAD) * cos(lat2 * DEG_TO_RAD) * cos(dlon)
bearing_rad = atan2(y, x)
bearing_deg = bearing_rad * RAD_TO_DEG

! Convert to 0-360 range
bearing_deg = mod(bearing_deg + 360.0_RKIND, 360.0_RKIND)

end function calculate_bearing

!==================================================================================================
function get_upstream_sector(bearing, n_sectors) result(sector)
!==================================================================================================
implicit none

real(kind=RKIND), intent(in) :: bearing
integer, intent(in) :: n_sectors
integer :: sector

real(kind=RKIND) :: sector_angle

sector_angle = 360.0_RKIND / real(n_sectors, RKIND)
sector = int(bearing / sector_angle) + 1

! Ensure sector is within valid range
if (sector > n_sectors) sector = n_sectors
if (sector < 1) sector = 1

end function get_upstream_sector

!==================================================================================================
subroutine process_sectors(nUpstreamSect, sector_counts, sector_has_data, &
                          sector_sums_zr, sector_sums_zr_std, sector_sums_rf, sector_sums_rd, &
                          current_zr, current_zr_std, current_rf, current_rd, &
                          result_zr, result_zr_std, result_rf, result_rd)
!==================================================================================================
implicit none

integer, intent(in) :: nUpstreamSect
integer, dimension(nUpstreamSect), intent(in) :: sector_counts
logical, dimension(nUpstreamSect), intent(in) :: sector_has_data
real(kind=RKIND), dimension(nUpstreamSect), intent(in) :: sector_sums_zr, sector_sums_zr_std
real(kind=RKIND), dimension(nUpstreamSect), intent(in) :: sector_sums_rf, sector_sums_rd
real(kind=RKIND), intent(in) :: current_zr, current_zr_std, current_rf, current_rd
real(kind=RKIND), dimension(nUpstreamSect), intent(out) :: result_zr, result_zr_std
real(kind=RKIND), dimension(nUpstreamSect), intent(out) :: result_rf, result_rd

integer :: iSector, left_idx, right_idx, j
real(kind=RKIND) :: upstream_zr, upstream_zr_std, upstream_rf, upstream_rd
real(kind=RKIND), parameter :: MISSING_VALUE = -999.0_RKIND

! First pass: calculate averages for sectors with data
do iSector = 1, nUpstreamSect
    if (sector_has_data(iSector) .and. sector_counts(iSector) > 0) then
        upstream_zr = sector_sums_zr(iSector) / real(sector_counts(iSector), RKIND)
        upstream_zr_std = sector_sums_zr_std(iSector) / real(sector_counts(iSector), RKIND)
        upstream_rf = sector_sums_rf(iSector) / real(sector_counts(iSector), RKIND)
        upstream_rd = sector_sums_rd(iSector) / real(sector_counts(iSector), RKIND)
        
        ! Average with current cell parameters
        result_zr(iSector) = (current_zr + upstream_zr) * 0.5_RKIND
        result_zr_std(iSector) = (current_zr_std + upstream_zr_std) * 0.5_RKIND
        result_rf(iSector) = (current_rf + upstream_rf) * 0.5_RKIND
        result_rd(iSector) = (current_rd + upstream_rd) * 0.5_RKIND
    else
        ! Mark as needing interpolation
        result_zr(iSector) = MISSING_VALUE
        result_zr_std(iSector) = MISSING_VALUE
        result_rf(iSector) = MISSING_VALUE
        result_rd(iSector) = MISSING_VALUE
    end if
end do

! Second pass: interpolate missing sectors
do iSector = 1, nUpstreamSect
    if (result_zr(iSector) < -900.0_RKIND) then  ! Needs interpolation
        
        ! Find nearest valid sectors
        left_idx = 0
        right_idx = 0
        
        ! Search left
        do j = 1, nUpstreamSect-1
            left_idx = iSector - j
            if (left_idx < 1) left_idx = left_idx + nUpstreamSect
            if (result_zr(left_idx) > -900.0_RKIND) exit
            left_idx = 0
        end do
        
        ! Search right
        do j = 1, nUpstreamSect-1
            right_idx = iSector + j
            if (right_idx > nUpstreamSect) right_idx = right_idx - nUpstreamSect
            if (result_zr(right_idx) > -900.0_RKIND) exit
            right_idx = 0
        end do
        
        ! Interpolate or use current cell values
        if (left_idx > 0 .and. right_idx > 0) then
            ! Average of left and right
            result_zr(iSector) = (result_zr(left_idx) + result_zr(right_idx)) * 0.5_RKIND
            result_zr_std(iSector) = (result_zr_std(left_idx) + result_zr_std(right_idx)) * 0.5_RKIND
            result_rf(iSector) = (result_rf(left_idx) + result_rf(right_idx)) * 0.5_RKIND
            result_rd(iSector) = (result_rd(left_idx) + result_rd(right_idx)) * 0.5_RKIND
        else if (left_idx > 0) then
            ! Use left
            result_zr(iSector) = result_zr(left_idx)
            result_zr_std(iSector) = result_zr_std(left_idx)
            result_rf(iSector) = result_rf(left_idx)
            result_rd(iSector) = result_rd(left_idx)
        else if (right_idx > 0) then
            ! Use right
            result_zr(iSector) = result_zr(right_idx)
            result_zr_std(iSector) = result_zr_std(right_idx)
            result_rf(iSector) = result_rf(right_idx)
            result_rd(iSector) = result_rd(right_idx)
        else
            ! Use current cell values
            result_zr(iSector) = current_zr
            result_zr_std(iSector) = current_zr_std
            result_rf(iSector) = current_rf
            result_rd(iSector) = current_rd
        end if
    end if
end do

end subroutine process_sectors

!==================================================================================================
subroutine calculate_derived_params(zr, rf, rd, zd, H_RSL, Z0C)
!==================================================================================================
! Calculate zero-plane displacement (zd), roughness sublayer height (H_RSL), 
! and roughness length (Z0C) based on building parameters
! Reference: module_sf_urban.F Z0C_TBL calculation
!==================================================================================================
implicit none

real(kind=RKIND), intent(in) :: zr, rf, rd
real(kind=RKIND), intent(out) :: zd, H_RSL, Z0C

! Local variables (only keep necessary ones)
real(kind=RKIND) :: HGT, R_TBL

! Parameters for Z0C calculation (from module_sf_urban.F)
real(kind=RKIND), parameter :: VonK = 0.4_RKIND
real(kind=RKIND), parameter :: Cd = 1.2_RKIND
real(kind=RKIND), parameter :: alpha_macd = 4.43_RKIND
real(kind=RKIND), parameter :: beta_macd = 1.0_RKIND

! Check for valid building parameters
if (zr <= 0.0_RKIND .or. rf <= 0.0_RKIND .or. rd <= 0.0_RKIND) then
    zd = 0.0_RKIND
    H_RSL = 0.0_RKIND
    Z0C = 0.1_RKIND  ! Default roughness length for non-urban areas
    return
end if

! Calculate normalized parameters following module_sf_urban.F (directly use rd and rf)
HGT = zr / (rd + rf)     ! Frontal area index (normalized height for 2-d canyon)
R_TBL = rf / (rd + rf)   ! Plan areal fraction (building coverage ratio for 2-d canyon)

! Calculate zero-plane displacement following Macdonald's (1998) formulation
zd = zr * (1.0_RKIND + (alpha_macd ** (-R_TBL)) * (R_TBL - 1.0_RKIND))

! Ensure zd is reasonable (0 <= zd < zr)
if (zd < 0.0_RKIND) then
    zd = 0.0_RKIND
else if (zd >= zr) then
    zd = 0.8_RKIND * zr  ! Limit to 80% of building height
end if

! Calculate roughness sublayer height (3 times building height)
H_RSL = 3.0_RKIND * zr
if (H_RSL < 10.0_RKIND) then
    H_RSL = 10.0_RKIND  ! Minimum value
end if

! Calculate roughness length following module_sf_urban.F
if (zr > zd) then
    Z0C = zr * (1.0_RKIND - zd/zr) * &
          exp(-(0.5_RKIND * beta_macd * Cd / (VonK**2) * &
               (1.0_RKIND - zd/zr) * HGT)**(-0.5_RKIND))
    
    ! Ensure Z0C is reasonable
    if (Z0C <= 0.0_RKIND .or. Z0C > zr) then
        Z0C = 0.1_RKIND * zr  ! Default to 10% of building height
    end if
else
    Z0C = 0.1_RKIND  ! Default for flat areas
end if

end subroutine calculate_derived_params

end module AddUpStreamBuildParam_mod
