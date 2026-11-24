!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to generate emisson_update file.
! * Author: Zihan Xia from AEMOL USTC, xiazh@mail.ustc.edu.cn
! Last update: zhxia 04/2026
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! References:
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
 module iamas_init_atm_emission
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_timekeeping
 use iamas_timer
 use iamas_log, only : iamas_log_write
 
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met
 use iamas_kd_tree

 implicit none
 private
 public :: init_atm_case_emis

 ! vars related to kdtree
 type (iamas_kd_type), pointer :: kdTree => null()
 type (iamas_kd_type), dimension(:), pointer :: kdPoints
 type (iamas_kd_type), pointer :: res
 logical :: kdTree_initialized = .false.

 ! map_table used in area_average_interp method. zhxia@2025
  type map_table
    integer, dimension(:,:), allocatable :: map_index_i  ! (nCells, max_points)
    integer, dimension(:,:), allocatable :: map_index_j  ! (nCells, max_points)
    integer, dimension(:), allocatable :: point_counts ! (nCells)
    integer :: max_points = 0
  end type map_table
  
  type(map_table), save :: anthro_map, bioms_map
  logical, save :: anthro_map_initialized = .false.
  logical, save :: bioms_map_initialized = .false.

  ! interp_weight_table: precomputed bilinear weights for point_counts<=1 cells
  ! avoids repeated rslab alloc/copy for every species in 'variable' interp
  ! path. zhxia@202604
  type interp_weight_table
    integer,          dimension(:),   allocatable :: cell_indices ! which iCells need interp
    real(kind=RKIND), dimension(:,:), allocatable :: weights      ! (4, n_cells) bilinear weights
    integer,          dimension(:,:), allocatable :: src_i        ! (4, n_cells) source i index
    integer,          dimension(:,:), allocatable :: src_j        ! (4, n_cells) source j index
    logical,          dimension(:),   allocatable :: use_fillval  ! .true. if masked
    integer :: n_cells = 0
  end type interp_weight_table

  type(interp_weight_table), save :: anthro_wt, bioms_wt
  logical, save :: anthro_wt_initialized  = .false.
  logical, save :: bioms_wt_initialized   = .false.

  ! cached proj_info to avoid rebuilding every species call (optimization point 3). zhxia@202604
  type(proj_info), save :: anthro_proj, bioms_proj
  logical, save :: anthro_proj_initialized = .false.
  logical, save :: bioms_proj_initialized  = .false.

 contains

!==================================================================================================
 subroutine init_atm_case_emis(domain, dminfo, stream_manager, mesh, emis_state, state, dims, configs)
!==================================================================================================

 use iamas_stream_manager

 implicit none

!input arguments:
 type (domain_type), intent(inout) :: domain
 type (dm_info), intent(in)        :: dminfo
 type (iAMAS_streamManager_type), intent(inout) :: stream_manager
 type (iamas_pool_type), intent(inout) :: mesh
 type (iamas_pool_type), intent(inout) :: emis_state
 type (iamas_pool_type), intent(inout) :: state
 type (iamas_pool_type), intent(in)    :: dims
 type (iamas_pool_type), intent(in)    :: configs

!local variables:
 type (iAMAS_Time_type)  :: curr_time, stop_time
 character(len=StrKIND) :: timeString

 character(len=StrKIND), pointer :: config_anthro_emis_prefix
 character(len=StrKIND), pointer :: config_bioms_emis_prefix
 character(len=StrKIND), pointer :: config_emis_interp_type
 logical, pointer :: config_use_anthro_emis, config_use_bioms_emis
 character(len=StrKIND), pointer :: xtime
 integer :: ierr
 

!==================================================================================================

 call iamas_pool_get_config(configs, 'config_anthro_emis_prefix', config_anthro_emis_prefix)
 call iamas_pool_get_config(configs, 'config_bioms_emis_prefix', config_bioms_emis_prefix)
 call iamas_pool_get_config(configs, 'config_use_anthro_emis', config_use_anthro_emis)
 call iamas_pool_get_config(configs, 'config_use_bioms_emis', config_use_bioms_emis)
 call iamas_pool_get_config(configs, 'config_emis_interp_type', config_emis_interp_type)
 call iamas_pool_get_array(state, 'xtime', xtime)

!loop over all times:
 curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW) 
 stop_time = iamas_get_clock_time(domain % clock, iAMAS_STOP_TIME) 

 do while (curr_time <= stop_time)
    call iamas_get_time(curr_time, dateTimeString=timeString)
    xtime = timeString

!   call iamas_log_write('Processing '//trim(config_anthro_emis_prefix)//':'//timeString(1:13))

    !read the anthropogenic and biomass data from the emission file, and interpolate the data to the iAMAS grid:
    call interp_emis_to_iAMAS(timeString(1:13), mesh, emis_state, dims, dminfo, config_anthro_emis_prefix, config_bioms_emis_prefix, config_use_anthro_emis, config_use_bioms_emis, config_emis_interp_type)

    !write the interpolated anthropogenic/biomass field as a new time slice in the iAMAS output file:
    call iamas_stream_mgr_write(stream_manager, streamID='emission', ierr=ierr)
    call iamas_stream_mgr_reset_alarms(stream_manager, streamID='emission', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

    call iamas_advance_clock(domain % clock)
    curr_time = iamas_get_clock_time(domain % clock, iAMAS_NOW) 

 end do

!
! Deallocate and free the KD Tree
!
 if (kdTree_initialized) then
     call iamas_kd_free(kdTree)
     deallocate(kdPoints)
 end if

 if (allocated(anthro_map%map_index_i)) deallocate(anthro_map%map_index_i)
 if (allocated(anthro_map%map_index_j)) deallocate(anthro_map%map_index_j)
 if (allocated(anthro_map%point_counts)) deallocate(anthro_map%point_counts)
 if (allocated(bioms_map%map_index_i)) deallocate(bioms_map%map_index_i)
 if (allocated(bioms_map%map_index_j)) deallocate(bioms_map%map_index_j)
 if (allocated(bioms_map%point_counts)) deallocate(bioms_map%point_counts)

 ! free interp_weight_table caches. zhxia@202604
 if (allocated(anthro_wt%cell_indices)) deallocate(anthro_wt%cell_indices)
 if (allocated(anthro_wt%weights))      deallocate(anthro_wt%weights)
 if (allocated(anthro_wt%src_i))        deallocate(anthro_wt%src_i)
 if (allocated(anthro_wt%src_j))        deallocate(anthro_wt%src_j)
 if (allocated(anthro_wt%use_fillval))  deallocate(anthro_wt%use_fillval)
 if (allocated(bioms_wt%cell_indices))  deallocate(bioms_wt%cell_indices)
 if (allocated(bioms_wt%weights))       deallocate(bioms_wt%weights)
 if (allocated(bioms_wt%src_i))         deallocate(bioms_wt%src_i)
 if (allocated(bioms_wt%src_j))         deallocate(bioms_wt%src_j)
 if (allocated(bioms_wt%use_fillval))   deallocate(bioms_wt%use_fillval)
 anthro_wt_initialized  = .false.
 bioms_wt_initialized   = .false.
 anthro_proj_initialized = .false.
 bioms_proj_initialized  = .false.

 !
 ! Ensure that no output alarms are still ringing for the 'emission' stream after
 ! we exit the time loop above; the main run routine may write out all other
 ! output streams with ringing alarms.
 !
 call iamas_stream_mgr_reset_alarms(stream_manager, streamID='emission', direction=iAMAS_STREAM_OUTPUT, ierr=ierr)

 end subroutine init_atm_case_emis

!==================================================================================================
 subroutine interp_emis_to_iAMAS(timeString, mesh, emis_state, dims, dminfo, config_anthro_emis_prefix, config_bioms_emis_prefix, config_use_anthro_emis, config_use_bioms_emis, config_emis_interp_type)
!==================================================================================================

 use iamas_dmpar 

 implicit none

!input arguments:
 character(len=*), intent(in) :: timeString
 type (iamas_pool_type), intent(in) :: mesh
 type (iamas_pool_type), intent(in) :: dims
 type (dm_info), intent(in)   :: dminfo
 character(len=*), intent(in) :: config_anthro_emis_prefix, config_bioms_emis_prefix, config_emis_interp_type
 logical, intent(in) :: config_use_anthro_emis, config_use_bioms_emis

!inout arguments:
 type (iamas_pool_type), intent(inout) :: emis_state


!local variables:
 type(met_data) :: field !real*4 emission data.

 integer :: istatus
 integer :: masked
 integer, dimension(:), pointer :: mask_array

 real(kind=RKIND) :: fillval, maskval, msgval

 real(kind=RKIND), dimension(:), pointer :: destField1d => null()
 real(kind=RKIND), dimension(:), pointer :: e_nh3     ,e_co       ,e_no      ,e_no2     ,e_so2     , &
                                            e_bc      ,e_oc       ,e_pm25    ,e_pm10    ,e_so4     , &
                                            e_no3     ,e_c2h6     ,e_c3h8    ,e_c2h2    ,e_alk3    , &
                                            e_alk4    ,e_alk5     ,e_ethene  ,e_c3h6    ,e_ole1    , &
                                            e_ole2    ,e_aro1     ,e_aro2    ,e_hcho    ,e_ccho    , &
                                            e_rcho    ,e_acet     ,e_mek     ,e_isoprene,e_terp    , &
                                            e_sesq    ,e_phen     ,e_cres    ,e_meoh    ,e_gly     , &
                                            e_mgly    ,e_bacl     ,e_isoprod ,e_methacro,e_mvk     , &
                                            e_prod2   ,e_ch4      ,e_bald    ,e_hcooh   ,e_cco_oh  , &
                                            e_rco_oh  ,e_co2
 real(kind=RKIND), dimension(:), pointer :: ebu_bc    ,ebu_oc     ,ebu_pm25  ,ebu_pm10  ,ebu_ch4   , &
                                            ebu_co2   ,ebu_co     ,ebu_no    ,ebu_no2   ,ebu_so2   , &
                                            ebu_acet  ,ebu_alk3   ,ebu_alk4  ,ebu_alk5  ,ebu_aro1  , &
                                            ebu_aro2  ,ebu_bald   ,ebu_ccho  ,ebu_c2o2h ,ebu_ethen , &
                                            ebu_hcho  ,ebu_hcooh  ,ebu_hono  ,ebu_isopr ,ebu_mek   , &
                                            ebu_meoh  ,ebu_metha  ,ebu_mgly  ,ebu_mvk   ,ebu_ole1  , &
                                            ebu_ole2  ,ebu_phen   ,ebu_prod2 ,ebu_rcho
 integer, pointer :: nCells

!==================================================================================================
 call iamas_pool_get_dimension(dims, 'nCells', nCells)

 ! initialize mask and fill values early (needed by precompute_interp_weights). zhxia@202604
 call iamas_pool_get_array(mesh, 'landmask', mask_array)
 msgval  = 0.0_RKIND
 masked  = -1
 maskval = -1.0_RKIND
 fillval =  0.0_RKIND

 ! initialize kdtree and map_table
 if (trim(config_emis_interp_type) .ne. 'coarse2fine') then 
    if (.not. kdTree_initialized) then
       call iamas_log_write("Create kdtree for map_index")
       call init_kd_tree(mesh,nCells)
    endif
    if (config_use_anthro_emis) then
       if (.not. anthro_map_initialized) then
          call read_met_init(trim(config_anthro_emis_prefix),.false.,timeString(1:13),istatus)
          call read_next_met_field(field,istatus)
          if (istatus == 0) then
             if (field % iproj .ne. PROJ_LATLON) then
                 call iamas_log_write("Error: config_emis_interp_type="//trim(config_emis_interp_type)//" but field % iproj != PROJ_LATLON", messageType=iAMAS_LOG_ERR)
             endif
             call iamas_log_write("Create map_table for anthro")
             call precompute_map_table(mesh, field, nCells, anthro_map)
             anthro_map_initialized = .true.
             ! precompute proj and interp weight table for 'variable' mode. zhxia@202604
             if (trim(config_emis_interp_type) == 'variable' .and. .not. anthro_wt_initialized) then
                call iamas_log_write("Precompute interp_weight_table for anthro")
                call precompute_interp_weights(mesh, field, nCells, anthro_map, mask_array, masked, &
                                               fillval, anthro_proj, anthro_wt)
                anthro_wt_initialized   = .true.
                anthro_proj_initialized = .true.
             end if
             deallocate(field%slab)
          endif
          call read_met_close()
       endif
    endif
    if (config_use_bioms_emis) then
       if (.not. bioms_map_initialized) then
          call read_met_init(trim(config_bioms_emis_prefix),.false.,timeString(1:13),istatus)
          call read_next_met_field(field,istatus)
          if (istatus == 0) then
             if (field % iproj .ne. PROJ_LATLON) then
                 call iamas_log_write("Error: config_emis_interp_type="//trim(config_emis_interp_type)//" but field % iproj != PROJ_LATLON", messageType=iAMAS_LOG_ERR)
             endif
             call iamas_log_write("Create map_table for bioms")
             call precompute_map_table(mesh, field, nCells, bioms_map)
             bioms_map_initialized = .true.
             ! precompute proj and interp weight table for 'variable' mode. zhxia@202604
             if (trim(config_emis_interp_type) == 'variable' .and. .not. bioms_wt_initialized) then
                call iamas_log_write("Precompute interp_weight_table for bioms")
                call precompute_interp_weights(mesh, field, nCells, bioms_map, mask_array, masked, &
                                               fillval, bioms_proj, bioms_wt)
                bioms_wt_initialized   = .true.
                bioms_proj_initialized = .true.
             end if
             deallocate(field%slab)
          endif
          call read_met_close()
       endif
    endif
 endif 

 call iamas_pool_get_array(emis_state, 'e_co2'     , e_co2)
 call iamas_pool_get_array(emis_state, 'e_nh3'     , e_nh3)
 call iamas_pool_get_array(emis_state, 'e_co'      , e_co)
 call iamas_pool_get_array(emis_state, 'e_no'      , e_no)
 call iamas_pool_get_array(emis_state, 'e_no2'     , e_no2)
 call iamas_pool_get_array(emis_state, 'e_so2'     , e_so2)
 call iamas_pool_get_array(emis_state, 'e_bc'      , e_bc)
 call iamas_pool_get_array(emis_state, 'e_oc'      , e_oc)
 call iamas_pool_get_array(emis_state, 'e_pm25'    , e_pm25)
 call iamas_pool_get_array(emis_state, 'e_pm10'    , e_pm10)
 call iamas_pool_get_array(emis_state, 'e_so4'     , e_so4)
 call iamas_pool_get_array(emis_state, 'e_no3'     , e_no3)
 call iamas_pool_get_array(emis_state, 'e_c2h6'    , e_c2h6)
 call iamas_pool_get_array(emis_state, 'e_c3h8'    , e_c3h8)
 call iamas_pool_get_array(emis_state, 'e_c2h2'    , e_c2h2)
 call iamas_pool_get_array(emis_state, 'e_alk3'    , e_alk3)
 call iamas_pool_get_array(emis_state, 'e_alk4'    , e_alk4)
 call iamas_pool_get_array(emis_state, 'e_alk5'    , e_alk5)
 call iamas_pool_get_array(emis_state, 'e_ethene'  , e_ethene)
 call iamas_pool_get_array(emis_state, 'e_c3h6'    , e_c3h6)
 call iamas_pool_get_array(emis_state, 'e_ole1'    , e_ole1)
 call iamas_pool_get_array(emis_state, 'e_ole2'    , e_ole2)
 call iamas_pool_get_array(emis_state, 'e_aro1'    , e_aro1)
 call iamas_pool_get_array(emis_state, 'e_aro2'    , e_aro2)
 call iamas_pool_get_array(emis_state, 'e_hcho'    , e_hcho)
 call iamas_pool_get_array(emis_state, 'e_ccho'    , e_ccho)
 call iamas_pool_get_array(emis_state, 'e_rcho'    , e_rcho)
 call iamas_pool_get_array(emis_state, 'e_acet'    , e_acet)
 call iamas_pool_get_array(emis_state, 'e_mek'     , e_mek)
 call iamas_pool_get_array(emis_state, 'e_isoprene', e_isoprene)
 call iamas_pool_get_array(emis_state, 'e_terp'    , e_terp)
 call iamas_pool_get_array(emis_state, 'e_sesq'    , e_sesq)
 call iamas_pool_get_array(emis_state, 'e_phen'    , e_phen)
 call iamas_pool_get_array(emis_state, 'e_cres'    , e_cres)
 call iamas_pool_get_array(emis_state, 'e_meoh'    , e_meoh)
 call iamas_pool_get_array(emis_state, 'e_gly'     , e_gly)
 call iamas_pool_get_array(emis_state, 'e_mgly'    , e_mgly)
 call iamas_pool_get_array(emis_state, 'e_bacl'    , e_bacl)
 call iamas_pool_get_array(emis_state, 'e_isoprod' , e_isoprod)
 call iamas_pool_get_array(emis_state, 'e_methacro', e_methacro)
 call iamas_pool_get_array(emis_state, 'e_mvk'     , e_mvk)
 call iamas_pool_get_array(emis_state, 'e_prod2'   , e_prod2)
 call iamas_pool_get_array(emis_state, 'e_ch4'     , e_ch4)
 call iamas_pool_get_array(emis_state, 'e_bald'    , e_bald)
 call iamas_pool_get_array(emis_state, 'e_hcooh'   , e_hcooh)
 call iamas_pool_get_array(emis_state, 'e_cco_oh'  , e_cco_oh)
 call iamas_pool_get_array(emis_state, 'e_rco_oh'  , e_rco_oh)

 call iamas_pool_get_array(emis_state, 'ebu_oc'    , ebu_oc)
 call iamas_pool_get_array(emis_state, 'ebu_bc'    , ebu_bc)
 call iamas_pool_get_array(emis_state, 'ebu_pm25'  , ebu_pm25)
 call iamas_pool_get_array(emis_state, 'ebu_pm10'  , ebu_pm10)
 call iamas_pool_get_array(emis_state, 'ebu_ch4'   , ebu_ch4)
 call iamas_pool_get_array(emis_state, 'ebu_co2'   , ebu_co2)
 call iamas_pool_get_array(emis_state, 'ebu_co'    , ebu_co)
 call iamas_pool_get_array(emis_state, 'ebu_no'    , ebu_no)
 call iamas_pool_get_array(emis_state, 'ebu_no2'   , ebu_no2)
 call iamas_pool_get_array(emis_state, 'ebu_so2'   , ebu_so2)
 call iamas_pool_get_array(emis_state, 'ebu_acet'  , ebu_acet)
 call iamas_pool_get_array(emis_state, 'ebu_alk3'  , ebu_alk3)
 call iamas_pool_get_array(emis_state, 'ebu_alk4'  , ebu_alk4)
 call iamas_pool_get_array(emis_state, 'ebu_alk5'  , ebu_alk5)
 call iamas_pool_get_array(emis_state, 'ebu_aro1'  , ebu_aro1)
 call iamas_pool_get_array(emis_state, 'ebu_aro2'  , ebu_aro2)
 call iamas_pool_get_array(emis_state, 'ebu_bald'  , ebu_bald)
 call iamas_pool_get_array(emis_state, 'ebu_ccho'  , ebu_ccho)
 call iamas_pool_get_array(emis_state, 'ebu_c2o2h' , ebu_c2o2h) !ebu_c2o2h=ebu_cco_oh
 call iamas_pool_get_array(emis_state, 'ebu_ethen' , ebu_ethen) !ebu_ethen=ebu_ethene
 call iamas_pool_get_array(emis_state, 'ebu_hcho'  , ebu_hcho)
 call iamas_pool_get_array(emis_state, 'ebu_hcooh' , ebu_hcooh)
 call iamas_pool_get_array(emis_state, 'ebu_hono'  , ebu_hono)
 call iamas_pool_get_array(emis_state, 'ebu_isopr' , ebu_isopr) !ebu_isopr=ebu_isoprene
 call iamas_pool_get_array(emis_state, 'ebu_mek'   , ebu_mek)
 call iamas_pool_get_array(emis_state, 'ebu_meoh'  , ebu_meoh)
 call iamas_pool_get_array(emis_state, 'ebu_metha' , ebu_metha) !ebu_metha=ebu_methacro
 call iamas_pool_get_array(emis_state, 'ebu_mgly'  , ebu_mgly)
 call iamas_pool_get_array(emis_state, 'ebu_mvk'   , ebu_mvk)
 call iamas_pool_get_array(emis_state, 'ebu_ole1'  , ebu_ole1)
 call iamas_pool_get_array(emis_state, 'ebu_ole2'  , ebu_ole2)
 call iamas_pool_get_array(emis_state, 'ebu_phen'  , ebu_phen)
 call iamas_pool_get_array(emis_state, 'ebu_prod2' , ebu_prod2)
 call iamas_pool_get_array(emis_state, 'ebu_rcho'  , ebu_rcho)

 !interpolation to the iAMAS grid:
 call iamas_log_write('Interpolating with interp type:'//trim(config_emis_interp_type))

 call iamas_log_write('Processing emission file for date:'//timeString(1:13))
 if (config_use_anthro_emis) then
     !read anthropogenic emission data. open intermediate file:
     call read_met_init(trim(config_anthro_emis_prefix),.false.,timeString(1:13),istatus)
     call iamas_log_write('open anthropogenic emission data')
     if(istatus /= 0) then
        call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
        call iamas_log_write('Error opening anthropogenic emission file '//trim(config_anthro_emis_prefix)//':'//timeString(1:13), messageType=iAMAS_LOG_ERR)
        call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
     end if
  
     !scan through all fields in the file, looking for the anthropogenic emission field:
     call read_next_met_field(field,istatus)
     do while (istatus == 0)
      !anthropogenic emission data:
        if(trim(field % field) == 'e_co2') then
            e_co2(1:nCells) = 0.0_RKIND
            destField1d => e_co2
 
        else if(trim(field % field) == 'e_nh3') then
            e_nh3(1:nCells) = 0.0_RKIND
            destField1d => e_nh3

        else if(trim(field % field) == 'e_co') then
            e_co(1:nCells) = 0.0_RKIND
            destField1d => e_co

        else if(trim(field % field) == 'e_no') then
            e_no(1:nCells) = 0.0_RKIND
            destField1d => e_no

        else if(trim(field % field) == 'e_no2') then
            e_no2(1:nCells) = 0.0_RKIND
            destField1d => e_no2

        else if(trim(field % field) == 'e_so2') then
            e_so2(1:nCells) = 0.0_RKIND
            destField1d => e_so2

        else if(trim(field % field) == 'e_bc') then
            e_bc(1:nCells) = 0.0_RKIND
            destField1d => e_bc

        else if(trim(field % field) == 'e_oc') then
            e_oc(1:nCells) = 0.0_RKIND
            destField1d => e_oc

        else if(trim(field % field) == 'e_pm25') then
            e_pm25(1:nCells) = 0.0_RKIND
            destField1d => e_pm25

        else if(trim(field % field) == 'e_pm10') then
            e_pm10(1:nCells) = 0.0_RKIND
            destField1d => e_pm10

        else if(trim(field % field) == 'e_so4') then
            e_so4(1:nCells) = 0.0_RKIND
            destField1d => e_so4

        else if(trim(field % field) == 'e_no3') then
            e_no3(1:nCells) = 0.0_RKIND
            destField1d => e_no3

        else if(trim(field % field) == 'e_c2h6') then
            e_c2h6(1:nCells) = 0.0_RKIND
            destField1d => e_c2h6

        else if(trim(field % field) == 'e_c3h8') then
            e_c3h8(1:nCells) = 0.0_RKIND
            destField1d => e_c3h8

        else if (trim(field % field) == 'e_c2h2') then
            e_c2h2(1:nCells) = 0.0_RKIND
            destField1d => e_c2h2

        else if (trim(field % field) == 'e_alk3') then
            e_alk3(1:nCells) = 0.0_RKIND
            destField1d => e_alk3

        else if (trim(field % field) == 'e_alk4') then
            e_alk4(1:nCells) = 0.0_RKIND
            destField1d => e_alk4

        else if (trim(field % field) == 'e_alk5') then
            e_alk5(1:nCells) = 0.0_RKIND
            destField1d => e_alk5

        else if (trim(field % field) == 'e_ethene') then
            e_ethene(1:nCells) = 0.0_RKIND
            destField1d => e_ethene

        else if (trim(field % field) == 'e_c3h6') then
            e_c3h6(1:nCells) = 0.0_RKIND
            destField1d => e_c3h6

        else if (trim(field % field) == 'e_ole1') then
            e_ole1(1:nCells) = 0.0_RKIND
            destField1d => e_ole1

        else if (trim(field % field) == 'e_ole2') then
            e_ole2(1:nCells) = 0.0_RKIND
            destField1d => e_ole2

        else if (trim(field % field) == 'e_aro1') then
            e_aro1(1:nCells) = 0.0_RKIND
            destField1d => e_aro1

        else if (trim(field % field) == 'e_aro2') then
            e_aro2(1:nCells) = 0.0_RKIND
            destField1d => e_aro2

        else if (trim(field % field) == 'e_hcho') then
            e_hcho(1:nCells) = 0.0_RKIND
            destField1d => e_hcho

        else if (trim(field % field) == 'e_ccho') then
            e_ccho(1:nCells) = 0.0_RKIND
            destField1d => e_ccho

        else if (trim(field % field) == 'e_rcho') then
            e_rcho(1:nCells) = 0.0_RKIND
            destField1d => e_rcho

        else if (trim(field % field) == 'e_acet') then
            e_acet(1:nCells) = 0.0_RKIND
            destField1d => e_acet

        else if (trim(field % field) == 'e_mek') then
            e_mek(1:nCells) = 0.0_RKIND
            destField1d => e_mek

        else if (trim(field % field) == 'e_isopren') then !e_isopren=e_isoprene
            e_isoprene(1:nCells) = 0.0_RKIND
            destField1d => e_isoprene

        else if (trim(field % field) == 'e_terp') then
            e_terp(1:nCells) = 0.0_RKIND
            destField1d => e_terp

        else if (trim(field % field) == 'e_sesq') then
            e_sesq(1:nCells) = 0.0_RKIND
            destField1d => e_sesq

        else if (trim(field % field) == 'e_phen') then
            e_phen(1:nCells) = 0.0_RKIND
            destField1d => e_phen

        else if (trim(field % field) == 'e_cres') then
            e_cres(1:nCells) = 0.0_RKIND
            destField1d => e_cres

        else if (trim(field % field) == 'e_meoh') then
            e_meoh(1:nCells) = 0.0_RKIND
            destField1d => e_meoh

        else if (trim(field % field) == 'e_gly') then
            e_gly(1:nCells) = 0.0_RKIND
            destField1d => e_gly

        else if (trim(field % field) == 'e_mgly') then
            e_mgly(1:nCells) = 0.0_RKIND
            destField1d => e_mgly

        else if (trim(field % field) == 'e_bacl') then
            e_bacl(1:nCells) = 0.0_RKIND
            destField1d => e_bacl

        else if (trim(field % field) == 'e_isoprod') then
            e_isoprod(1:nCells) = 0.0_RKIND
            destField1d => e_isoprod

        else if (trim(field % field) == 'e_methacr') then !e_methacr=e_methacro
            e_methacro(1:nCells) = 0.0_RKIND
            destField1d => e_methacro

        else if (trim(field % field) == 'e_mvk') then
            e_mvk(1:nCells) = 0.0_RKIND
            destField1d => e_mvk

        else if (trim(field % field) == 'e_prod2') then
            e_prod2(1:nCells) = 0.0_RKIND
            destField1d => e_prod2

        else if (trim(field % field) == 'e_ch4') then
            e_ch4(1:nCells) = 0.0_RKIND
            destField1d => e_ch4

        else if (trim(field % field) == 'e_bald') then
            e_bald(1:nCells) = 0.0_RKIND
            destField1d => e_bald

        else if (trim(field % field) == 'e_hcooh') then
            e_hcooh(1:nCells) = 0.0_RKIND
            destField1d => e_hcooh

        else if (trim(field % field) == 'e_cco_oh') then
            e_cco_oh(1:nCells) = 0.0_RKIND
            destField1d => e_cco_oh

        else if (trim(field % field) == 'e_rco_oh') then
            e_rco_oh(1:nCells) = 0.0_RKIND
            destField1d => e_rco_oh

        else
            ! unrecognized field, skip without interpolating
            deallocate(field%slab)
            call read_next_met_field(field, istatus)
            cycle
 
        end if
    
        !call iamas_log_write(trim(field % field)//' processing')
        if (anthro_wt_initialized .and. anthro_proj_initialized) then
           call interp_to_iAMAS(mesh,nCells,field,destField1d,config_emis_interp_type,msgval,masked,maskval,fillval, &
                                mask_array,anthro_map,cached_proj=anthro_proj,cached_wt=anthro_wt)
        else
           call interp_to_iAMAS(mesh,nCells,field,destField1d,config_emis_interp_type,msgval,masked,maskval,fillval, &
                                mask_array,anthro_map)
        end if

        deallocate(field%slab)
        call read_next_met_field(field,istatus)
     end do

     !close intermediate file 
     call read_met_close()
 endif

 if (config_use_bioms_emis) then
     !read biomass burning emission data. open intermediate file:
     call read_met_init(trim(config_bioms_emis_prefix),.false.,timeString(1:13),istatus)
     call iamas_log_write('open biomass burning emission data')
     if(istatus /= 0) then
        call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
        call iamas_log_write('Error opening biomass burning emission file '//trim(config_bioms_emis_prefix)//':'//timeString(1:13), messageType=iAMAS_LOG_ERR)
        call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
     end if

     !scan through all fields in the file, looking for the biomass burning emission field:
     call read_next_met_field(field,istatus)
     do while (istatus == 0)
      !biomass burning emission data:
        if(trim(field % field) == 'ebu_bc') then
            ebu_bc(1:nCells) = 0.0_RKIND
            destField1d => ebu_bc

        else if(trim(field % field) == 'ebu_oc') then
            ebu_oc(1:nCells) = 0.0_RKIND
            destField1d => ebu_oc

        else if(trim(field % field) == 'ebu_pm25') then
            ebu_pm25(1:nCells) = 0.0_RKIND
            destField1d => ebu_pm25

        else if(trim(field % field) == 'ebu_pm10') then
            ebu_pm10(1:nCells) = 0.0_RKIND
            destField1d => ebu_pm10

        else if(trim(field % field) == 'ebu_ch4') then
            ebu_ch4(1:nCells) = 0.0_RKIND
            destField1d => ebu_ch4

        else if(trim(field % field) == 'ebu_co2') then
            ebu_co2(1:nCells) = 0.0_RKIND
            destField1d => ebu_co2

        else if(trim(field % field) == 'ebu_co') then
            ebu_co(1:nCells) = 0.0_RKIND
            destField1d => ebu_co

        else if(trim(field % field) == 'ebu_no') then
            ebu_no(1:nCells) = 0.0_RKIND
            destField1d => ebu_no

        else if(trim(field % field) == 'ebu_no2') then
            ebu_no2(1:nCells) = 0.0_RKIND
            destField1d => ebu_no2

        else if(trim(field % field) == 'ebu_so2') then
            ebu_so2(1:nCells) = 0.0_RKIND
            destField1d => ebu_so2

        else if(trim(field % field) == 'ebu_acet') then
            ebu_acet(1:nCells) = 0.0_RKIND
            destField1d => ebu_acet

        else if(trim(field % field) == 'ebu_alk3') then
            ebu_alk3(1:nCells) = 0.0_RKIND
            destField1d => ebu_alk3

        else if(trim(field % field) == 'ebu_alk4') then
            ebu_alk4(1:nCells) = 0.0_RKIND
            destField1d => ebu_alk4

        else if(trim(field % field) == 'ebu_alk5') then
            ebu_alk5(1:nCells) = 0.0_RKIND
            destField1d => ebu_alk5

        else if(trim(field % field) == 'ebu_aro1') then
            ebu_aro1(1:nCells) = 0.0_RKIND
            destField1d => ebu_aro1

        else if(trim(field % field) == 'ebu_aro2') then
            ebu_aro2(1:nCells) = 0.0_RKIND
            destField1d => ebu_aro2

        else if(trim(field % field) == 'ebu_bald') then
            ebu_bald(1:nCells) = 0.0_RKIND
            destField1d => ebu_bald

        else if(trim(field % field) == 'ebu_ccho') then
            ebu_ccho(1:nCells) = 0.0_RKIND
            destField1d => ebu_ccho

        else if(trim(field % field) == 'ebu_c2o2h') then !ebu_c2o2h=ebu_cco_oh
            ebu_c2o2h(1:nCells) = 0.0_RKIND
            destField1d => ebu_c2o2h

        else if(trim(field % field) == 'ebu_ethen') then !ebu_ethen=ebu_ethene
            ebu_ethen(1:nCells) = 0.0_RKIND
            destField1d => ebu_ethen

        else if(trim(field % field) == 'ebu_hcho') then
            ebu_hcho(1:nCells) = 0.0_RKIND
            destField1d => ebu_hcho

        else if(trim(field % field) == 'ebu_hcooh') then
            ebu_hcooh(1:nCells) = 0.0_RKIND
            destField1d => ebu_hcooh

        else if(trim(field % field) == 'ebu_hono') then
            ebu_hono(1:nCells) = 0.0_RKIND
            destField1d => ebu_hono

        else if(trim(field % field) == 'ebu_isopr') then !ebu_isopr=ebu_isoprene
            ebu_isopr(1:nCells) = 0.0_RKIND
            destField1d => ebu_isopr

        else if(trim(field % field) == 'ebu_mek') then
            ebu_mek(1:nCells) = 0.0_RKIND
            destField1d => ebu_mek

        else if(trim(field % field) == 'ebu_meoh') then
            ebu_meoh(1:nCells) = 0.0_RKIND
            destField1d => ebu_meoh

        else if(trim(field % field) == 'ebu_metha') then !ebu_metha=ebu_methacro
            ebu_metha(1:nCells) = 0.0_RKIND
            destField1d => ebu_metha

        else if(trim(field % field) == 'ebu_mgly') then
            ebu_mgly(1:nCells) = 0.0_RKIND
            destField1d => ebu_mgly

        else if(trim(field % field) == 'ebu_mvk') then
            ebu_mvk(1:nCells) = 0.0_RKIND
            destField1d => ebu_mvk

        else if(trim(field % field) == 'ebu_ole1') then
            ebu_ole1(1:nCells) = 0.0_RKIND
            destField1d => ebu_ole1

        else if(trim(field % field) == 'ebu_ole2') then
            ebu_ole2(1:nCells) = 0.0_RKIND
            destField1d => ebu_ole2

        else if(trim(field % field) == 'ebu_phen') then
            ebu_phen(1:nCells) = 0.0_RKIND
            destField1d => ebu_phen

        else if(trim(field % field) == 'ebu_prod2') then
            ebu_prod2(1:nCells) = 0.0_RKIND
            destField1d => ebu_prod2

        else if(trim(field % field) == 'ebu_rcho') then
            ebu_rcho(1:nCells) = 0.0_RKIND
            destField1d => ebu_rcho

        else
            ! unrecognized field, skip without interpolating. zhxia@202606
            deallocate(field%slab)
            call read_next_met_field(field, istatus)
            cycle

        end if

        !call iamas_log_write(trim(field % field)//' processing')
        if (bioms_wt_initialized .and. bioms_proj_initialized) then
           call interp_to_iAMAS(mesh,nCells,field,destField1d,config_emis_interp_type,msgval,masked,maskval,fillval, &
                                mask_array,bioms_map,cached_proj=bioms_proj,cached_wt=bioms_wt)
        else
           call interp_to_iAMAS(mesh,nCells,field,destField1d,config_emis_interp_type,msgval,masked,maskval,fillval, &
                                mask_array,bioms_map)
        end if

        deallocate(field%slab)
        call read_next_met_field(field,istatus)
     end do

     !close intermediate file
     call read_met_close()
 endif

 end subroutine interp_emis_to_iAMAS

!==================================================================================================
subroutine init_kd_tree(mesh,nCells)
!==================================================================================================
  type (iamas_pool_type), intent(in) :: mesh
  integer, intent(in) :: nCells 

  real(kind=RKIND), dimension(:), pointer :: latCell, lonCell, xCell, yCell, zCell
  integer :: iCell

  if (kdTree_initialized) return

  call iamas_pool_get_array(mesh, 'latCell', latCell)
  call iamas_pool_get_array(mesh, 'lonCell', lonCell)
  call iamas_pool_get_array(mesh, 'xCell', xCell)
  call iamas_pool_get_array(mesh, 'yCell', yCell)
  call iamas_pool_get_array(mesh, 'zCell', zCell)


  ! Initialize the KD-Tree
  allocate(kdPoints(nCells))
  do iCell = 1, nCells
    allocate(kdPoints(iCell) % point(3))
    kdPoints(iCell) % point = (/xCell(iCell), yCell(iCell), zCell(iCell)/)
    kdPoints(iCell) % id = iCell
  end do
  kdTree => iamas_kd_construct(kdPoints, 3)
  if (.not. associated(kdTree)) then
     call iamas_log_write('Error creating the KD-Tree for emission interplation', messageType=iAMAS_LOG_CRIT)
  endif
  kdTree_initialized = .true.
end subroutine init_kd_tree

!==================================================================================================
subroutine precompute_map_table(mesh, field, nCells, map_tbl)
!==================================================================================================
  type (iamas_pool_type), intent(in) :: mesh
  type (met_data), intent(in) :: field
  integer, intent(in) :: nCells
  type(map_table), intent(inout) :: map_tbl
  
  real(kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
  real(kind=RKIND) :: lat_deg, lon_deg, lat_rad, lon_rad, x, y, z, dist
  real(kind=RKIND) :: lat_max, lat_min, lon_max, lon_min
  integer :: i, j, iCell
  real(kind=RKIND), parameter :: deg2rad = 3.141592653589793_RKIND / 180.0_RKIND
  real(kind=RKIND), parameter :: buffer = 0.017 !1 deg
  
  ! allocate array
  if (allocated(map_tbl%point_counts)) deallocate(map_tbl%point_counts)
  allocate(map_tbl%point_counts(nCells))
  map_tbl%point_counts = 0
 
  ! calculate lat/lon range in this process
  call iamas_pool_get_array(mesh, 'latEdge', latEdge)
  call iamas_pool_get_array(mesh, 'lonEdge', lonEdge)
  lat_max = maxval(latEdge)
  lat_min = minval(latEdge)
  lon_max = maxval(lonEdge)
  lon_min = minval(lonEdge)

  ! Step1: count the number of lat-lon grids contained by each iamas grid
  ! loop all lat-lon grids
  do j = 1, field % ny
    ! calculate lat-lon of the current grid
    lat_deg = field % startlat + (j-1) * field % deltalat
    ! transfer to Cartesian coordinate
    lat_rad = lat_deg * deg2rad
    ! skip grid which in not included by this process
    if (lat_rad < lat_min - buffer .or. lat_rad > lat_max  + buffer) cycle

    do i = 1, field % nx
      ! calculate lat-lon of the current grid   
      lon_deg = field % startlon + (i-1) * field % deltalon      
      ! transfer to Cartesian coordinate
      lon_rad = lon_deg * deg2rad
      ! skip grid which in not included by this process
      if (lon_rad < lon_min - buffer .or. lon_rad > lon_max  + buffer) cycle

      x = cos(lat_rad) * cos(lon_rad)
      y = cos(lat_rad) * sin(lon_rad)
      z = sin(lat_rad)
     
      ! search iamas grid in kdtree
      call iamas_kd_search(kdTree, [x, y, z], res, dist)
    
      ! accumulate and count
      if (associated(res)) then
        iCell = res%id
        if (iCell > 0 .and. iCell <= nCells) then
          map_tbl%point_counts(iCell) = map_tbl%point_counts(iCell) + 1
        end if
      end if
    end do
  end do
  
  ! allocate map_table using thr max point_counts
  map_tbl%max_points = maxval(map_tbl%point_counts)
  if (allocated(map_tbl%map_index_i)) deallocate(map_tbl%map_index_i)
  if (allocated(map_tbl%map_index_j)) deallocate(map_tbl%map_index_j)
  allocate(map_tbl%map_index_i(nCells, map_tbl%max_points))
  allocate(map_tbl%map_index_j(nCells, map_tbl%max_points))
  map_tbl%map_index_i = 0
  map_tbl%map_index_j = 0
  map_tbl%point_counts = 0  ! reset counter
  
  ! Step2: write map_table
  ! loop all lat-lon grids
  do j = 1, field % ny
    ! calculate lat-lon of the current grid
    lat_deg = field % startlat + (j-1) * field % deltalat
    ! transfer to Cartesian coordinate
    lat_rad = lat_deg * deg2rad
    ! skip grid which in not included by this process
    if (lat_rad < lat_min - buffer .or. lat_rad > lat_max  + buffer) cycle

    do i = 1, field % nx
      ! calculate lat-lon of the current grid   
      lon_deg = field % startlon + (i-1) * field % deltalon      
      ! transfer to Cartesian coordinate
      lon_rad = lon_deg * deg2rad
      ! skip grid which in not included by this process
      if (lon_rad < lon_min - buffer .or. lon_rad > lon_max  + buffer) cycle

      x = cos(lat_rad) * cos(lon_rad)
      y = cos(lat_rad) * sin(lon_rad)
      z = sin(lat_rad)
      
      call iamas_kd_search(kdTree, [x, y, z], res, dist)
      if (associated(res)) then
        iCell = res%id
        if (iCell > 0 .and. iCell <= nCells) then
          map_tbl%point_counts(iCell) = map_tbl%point_counts(iCell) + 1
          map_tbl%map_index_i(iCell, map_tbl%point_counts(iCell)) = i
          map_tbl%map_index_j(iCell, map_tbl%point_counts(iCell)) = j
        end if
      end if
    end do
  end do
end subroutine precompute_map_table

!==================================================================================================
 subroutine interp_to_iAMAS(mesh,nCells,field,destField1d,config_emis_interp_type,msgval,masked,maskval,fillval, &
                            mask_array, current_map, cached_proj, cached_wt)
!==================================================================================================

!input arguments:
 type (iamas_pool_type), intent(in) :: mesh
 integer, intent(in) :: nCells
 type (met_data), intent(in)  :: field !real*4 meteorological data.
 character(len=*), intent(in) :: config_emis_interp_type

 integer, intent(in) :: masked
 integer, dimension(:), intent(in), pointer :: mask_array

 real(kind=RKIND), intent(in) :: fillval, maskval, msgval

 type(map_table), intent(in) :: current_map

 ! optional cached proj and weight table (optimization: avoid rebuilding every species). zhxia@202604
 type(proj_info),           intent(in), optional :: cached_proj
 type(interp_weight_table), intent(in), optional :: cached_wt

!inout arguments:
 real(kind=RKIND), intent(inout), dimension(:), pointer :: destField1d

!local variables:
 type(proj_info) :: proj
 integer, dimension(5) :: interp_list
 integer :: i, k, nInterpPoints, iCell
 real(kind=RKIND) :: lat,lon,x,y
 real(kind=RKIND), dimension(:,:), allocatable :: rslab

 real(kind=RKIND), dimension(:), pointer :: latPoints, lonPoints
 real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
 

!--------------------------------------------------------------------------------------------------


 ! zhxia@202604: reuse cached proj if provided, otherwise build it once
 if (present(cached_proj)) then
    proj = cached_proj
 else
    call map_init(proj)   
    if (field % iproj == PROJ_LATLON) then
       call map_set(PROJ_LATLON, proj, &
                    latinc = real(field % deltalat,RKIND), &
                    loninc = real(field % deltalon,RKIND), &
                    knowni = 1.0_RKIND, &
                    knownj = 1.0_RKIND, &
                    lat1 = real(field % startlat,RKIND), &
                    lon1 = real(field % startlon,RKIND))
       call iamas_log_write('--- The projection is PROJ_LATLON.')
    else if (field % iproj == PROJ_GAUSS) then
       call map_set(PROJ_GAUSS, proj, &
                    nlat = nint(field % deltalat), &
                    loninc = 360.0_RKIND / real(field % nx,RKIND), &
                    lat1 = real(field % startlat,RKIND), &
                    lon1 = real(field % startlon,RKIND))
       call iamas_log_write('--- The projection is PROJ_GAUSS.')
    else if (field % iproj == PROJ_PS) then
       call map_set(PROJ_PS, proj, &
                    dx = real(field % dx,RKIND), &
                    truelat1 = real(field % truelat1,RKIND), &
                    stdlon = real(field % xlonc,RKIND), &
                    knowni = real(field % nx / 2.0,RKIND), &
                    knownj = real(field % ny / 2.0,RKIND), &
                    lat1 = real(field % startlat,RKIND), &
                    lon1 = real(field % startlon,RKIND))
       call iamas_log_write('--- The projection is PROJ_PS.')
    end if
 end if
 
 if (field % iproj == PROJ_LATLON) then

      if (trim(config_emis_interp_type)=='fine2coarse') then
           !applicable to interpolating the original fine-resolution data to the coarse-resolution grid of IAMAS, using the regional average method
           call area_average_interp(field, nCells, destField1d, mask_array, masked, fillval, current_map)          
      else if (trim(config_emis_interp_type)=='variable') then
           !applicable to interpolating the original fine-resolution data to the variational grid of IAMAS - both fine and coarse

           ! zhxia@202604: area_average_interp skips point_counts<=1 cells (see subroutine)
           call area_average_interp(field, nCells, destField1d, mask_array, masked, fillval, current_map)

           if (present(cached_wt) .and. cached_wt%n_cells > 0) then
              ! zhxia@202604: use precomputed bilinear weights, no rslab alloc/copy ---
              do k = 1, cached_wt%n_cells
                 iCell = cached_wt%cell_indices(k)
                 if (cached_wt%use_fillval(k)) then
                    destField1d(iCell) = fillval
                 else
                    destField1d(iCell) = &
                         cached_wt%weights(1,k) * real(field%slab(cached_wt%src_i(1,k), cached_wt%src_j(1,k)), RKIND) + &
                         cached_wt%weights(2,k) * real(field%slab(cached_wt%src_i(2,k), cached_wt%src_j(2,k)), RKIND) + &
                         cached_wt%weights(3,k) * real(field%slab(cached_wt%src_i(3,k), cached_wt%src_j(3,k)), RKIND) + &
                         cached_wt%weights(4,k) * real(field%slab(cached_wt%src_i(4,k), cached_wt%src_j(4,k)), RKIND)
                 end if
              end do
           else
              ! --- fallback: original rslab path (used when weight table not available) zhxia@202604---
              interp_list(1) = FOUR_POINT
              interp_list(2) = AVERAGE4
              interp_list(3) = SEARCH
              interp_list(4) = 0

              call iamas_pool_get_array(mesh, 'latCell', latCell)
              call iamas_pool_get_array(mesh, 'lonCell', lonCell)

              nInterpPoints = nCells
              latPoints => latCell
              lonPoints => lonCell

              allocate(rslab(-2:field % nx+3, field % ny))
              rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
              rslab( 0, 1:field % ny) = field % slab(field % nx  , 1:field % ny)
              rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
              rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
              rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
              rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
              rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

              do i = 1,nInterpPoints
                 if (current_map%point_counts(i) .le. 1) then
                    if(mask_array(i) /= masked) then
                       lat = latPoints(i) * DEG_PER_RAD
                       lon = lonPoints(i) * DEG_PER_RAD
                       call latlon_to_ij(proj, lat, lon, x, y)
                       if(y <= 0.5) then
                          y = 1.0
                       else if(y >= real(field%ny)+0.5) then
                          y = real(field % ny)
                       end if
                       if(x < 0.5) then
                          lon = lon + 360.0
                          call latlon_to_ij(proj, lat, lon, x, y)
                       else if (x >= real(field%nx)+0.5) then
                          lon = lon - 360.0
                          call latlon_to_ij(proj, lat, lon, x, y)
                       end if
                       destField1d(i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                                        msgval,interp_list,1,maskval=maskval)
                    else
                       destField1d(i) = fillval
                    end if
                 endif
              end do
              deallocate(rslab)
           end if

      else if (trim(config_emis_interp_type)=='coarse2fine') then
           !applicable to the original  coarse-resolution data to the fine-resolution grid of IAMAS

           interp_list(1) = FOUR_POINT
           interp_list(2) = AVERAGE4
           interp_list(3) = SEARCH
           interp_list(4) = 0
          
           call iamas_pool_get_array(mesh, 'latCell', latCell)
           call iamas_pool_get_array(mesh, 'lonCell', lonCell)
          
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
          
           allocate(rslab(-2:field % nx+3, field % ny))
           rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
           rslab( 0, 1:field % ny) = field % slab(field % nx  , 1:field % ny)
           rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
           rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
           rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
           rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
           rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
          
           do i = 1,nInterpPoints
              if(mask_array(i) /= masked) then
                 lat = latPoints(i) * DEG_PER_RAD
                 lon = lonPoints(i) * DEG_PER_RAD
                 call latlon_to_ij(proj, lat, lon, x, y)
                 if(y <= 0.5) then
                    y = 1.0
                 else if(y >= real(field%ny)+0.5) then
                    y = real(field % ny)
                 end if
                 if(x < 0.5) then
                    lon = lon + 360.0
                    call latlon_to_ij(proj, lat, lon, x, y)
                 else if (x >= real(field%nx)+0.5) then
                    lon = lon - 360.0
                    call latlon_to_ij(proj, lat, lon, x, y)
                 end if
          
                 destField1d(i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                                  msgval,interp_list,1,maskval=maskval)
              else
                 destField1d(i) = fillval
              end if
           end do
           deallocate(rslab)

      endif

 else !other proj must use default interp method

      interp_list(1) = FOUR_POINT
      interp_list(2) = AVERAGE4
      interp_list(3) = SEARCH
      interp_list(4) = 0
     
      call iamas_pool_get_array(mesh, 'latCell', latCell)
      call iamas_pool_get_array(mesh, 'lonCell', lonCell)
     
      nInterpPoints = nCells
      latPoints => latCell
      lonPoints => lonCell
     
      allocate(rslab(-2:field % nx+3, field % ny))
      rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
      rslab( 0, 1:field % ny) = field % slab(field % nx  , 1:field % ny)
      rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
      rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
      rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
      rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
      rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
     
      do i = 1,nInterpPoints
         if(mask_array(i) /= masked) then
            lat = latPoints(i) * DEG_PER_RAD
            lon = lonPoints(i) * DEG_PER_RAD
            call latlon_to_ij(proj, lat, lon, x, y)
            if(y <= 0.5) then
               y = 1.0
            else if(y >= real(field%ny)+0.5) then
               y = real(field % ny)
            end if
            if(x < 0.5) then
               lon = lon + 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            else if (x >= real(field%nx)+0.5) then
               lon = lon - 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            end if
     
            destField1d(i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                             msgval,interp_list,1,maskval=maskval)
         else
            destField1d(i) = fillval
         end if
      end do
      deallocate(rslab)

 endif

 end subroutine interp_to_iAMAS

!==================================================================================================
subroutine area_average_interp(field, nCells, destField1d, mask_array, masked, fillval, current_map)
! applicable to interpolating the original fine-resolution data to the coarse-resolution grid of IAMAS, using the regional average method. zhxia@lacar. 2025/06/26
!==================================================================================================
  type (met_data), intent(in) :: field
  integer, intent(in) :: nCells
  real(kind=RKIND), dimension(:), intent(inout) :: destField1d
  integer, dimension(:), intent(in) :: mask_array
  integer, intent(in) :: masked
  real(kind=RKIND), intent(in) :: fillval
  type(map_table), intent(in) :: current_map
  
  integer :: i, j, k, iCell
  
  ! calculate area average
  ! zhxia@202604: skip point_counts<=1 cells — those will be filled by the
  ! interp_sequence / precomputed-weight path afterwards (avoids doing wasted work)
  do iCell = 1, nCells
    if (mask_array(iCell) == masked) then
      destField1d(iCell) = fillval
    else if (current_map%point_counts(iCell) > 1) then
      do k = 1, current_map%point_counts(iCell)
         i = current_map%map_index_i(iCell, k)
         j = current_map%map_index_j(iCell, k)
         destField1d(iCell) = destField1d(iCell) + real(field%slab(i,j), RKIND)
      enddo
      destField1d(iCell) = destField1d(iCell) / real(current_map%point_counts(iCell), RKIND)
    else if (current_map%point_counts(iCell) == 0) then
      destField1d(iCell) = fillval
    end if
    ! point_counts == 1: leave as-is (will be overwritten by weight-table / interp_sequence path)
  end do
end subroutine area_average_interp

!==================================================================================================
subroutine precompute_interp_weights(mesh, field, nCells, current_map, mask_array, masked, &
                                     fillval, out_proj, out_wt)
! Precomputes bilinear interpolation weights for all cells where point_counts<=1.
! Called once per emission source (anthro/bioms) during initialization.
! Eliminates the repeated rslab alloc/copy and latlon_to_ij calls in the 'variable' path.
! zhxia@202604
!==================================================================================================
  use init_atm_llxy

  type(iamas_pool_type),      intent(in)  :: mesh
  type(met_data),            intent(in)  :: field
  integer,                   intent(in)  :: nCells
  type(map_table),           intent(in)  :: current_map
  integer, dimension(:),     intent(in)  :: mask_array
  integer,                   intent(in)  :: masked
  real(kind=RKIND),          intent(in)  :: fillval
  type(proj_info),           intent(out) :: out_proj
  type(interp_weight_table), intent(out) :: out_wt

  real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
  real(kind=RKIND) :: lat, lon, x, y
  real(kind=RKIND) :: wx, wy                   ! fractional offsets within source cell
  integer :: i0, j0                            ! lower-left source cell indices (1-based)
  integer :: iCell, n, cnt
  real(kind=RKIND), parameter :: DEG_PER_RAD_L = 180.0_RKIND / 3.141592653589793_RKIND

  ! --- build proj (saved into out_proj so caller can cache it) ---
  call map_init(out_proj)
  if (field % iproj == PROJ_LATLON) then
     call map_set(PROJ_LATLON, out_proj, &
                  latinc = real(field % deltalat, RKIND), &
                  loninc = real(field % deltalon, RKIND), &
                  knowni = 1.0_RKIND, &
                  knownj = 1.0_RKIND, &
                  lat1   = real(field % startlat, RKIND), &
                  lon1   = real(field % startlon, RKIND))
  else if (field % iproj == PROJ_GAUSS) then
     call map_set(PROJ_GAUSS, out_proj, &
                  nlat   = nint(field % deltalat), &
                  loninc = 360.0_RKIND / real(field % nx, RKIND), &
                  lat1   = real(field % startlat, RKIND), &
                  lon1   = real(field % startlon, RKIND))
  else if (field % iproj == PROJ_PS) then
     call map_set(PROJ_PS, out_proj, &
                  dx       = real(field % dx,       RKIND), &
                  truelat1 = real(field % truelat1, RKIND), &
                  stdlon   = real(field % xlonc,    RKIND), &
                  knowni   = real(field % nx / 2.0, RKIND), &
                  knownj   = real(field % ny / 2.0, RKIND), &
                  lat1     = real(field % startlat, RKIND), &
                  lon1     = real(field % startlon, RKIND))
  end if

  call iamas_pool_get_array(mesh, 'latCell', latCell)
  call iamas_pool_get_array(mesh, 'lonCell', lonCell)

  ! --- first pass: count how many cells need interpolation ---
  cnt = 0
  do iCell = 1, nCells
     if (current_map%point_counts(iCell) <= 1) cnt = cnt + 1
  end do

  out_wt%n_cells = cnt
  if (cnt == 0) return

  allocate(out_wt%cell_indices(cnt))
  allocate(out_wt%weights(4, cnt))
  allocate(out_wt%src_i(4, cnt))
  allocate(out_wt%src_j(4, cnt))
  allocate(out_wt%use_fillval(cnt))

  ! --- second pass: compute and store weights for each cell ---
  n = 0
  do iCell = 1, nCells
     if (current_map%point_counts(iCell) > 1) cycle
     n = n + 1
     out_wt%cell_indices(n) = iCell

     if (mask_array(iCell) == masked) then
        out_wt%use_fillval(n) = .true.
        out_wt%weights(:, n)  = 0.0_RKIND
        out_wt%src_i(:, n)    = 1
        out_wt%src_j(:, n)    = 1
        cycle
     end if

     out_wt%use_fillval(n) = .false.

     lat = latCell(iCell) * DEG_PER_RAD_L
     lon = lonCell(iCell) * DEG_PER_RAD_L
     call latlon_to_ij(out_proj, lat, lon, x, y)

     ! clamp y to valid range
     if (y <= 0.5_RKIND) then
        y = 1.0_RKIND
     else if (y >= real(field%ny, RKIND) + 0.5_RKIND) then
        y = real(field%ny, RKIND)
     end if
     ! wrap x if needed
     if (x < 0.5_RKIND) then
        lon = lon + 360.0_RKIND
        call latlon_to_ij(out_proj, lat, lon, x, y)
     else if (x >= real(field%nx, RKIND) + 0.5_RKIND) then
        lon = lon - 360.0_RKIND
        call latlon_to_ij(out_proj, lat, lon, x, y)
     end if

     ! bilinear stencil: lower-left corner
     !i0 = max(1, min(int(x), field%nx - 1))
     i0 = floor(x)
     if (i0 < 1) i0 = field%nx
     if (i0 > field%nx) i0 = 1
     j0 = max(1, min(int(y), field%ny - 1))

     wx = x - real(i0, RKIND)   ! fractional offset in i (0..1)
     wy = y - real(j0, RKIND)   ! fractional offset in j (0..1)
     ! clamp to [0,1]
     wx = max(0.0_RKIND, min(1.0_RKIND, wx))
     wy = max(0.0_RKIND, min(1.0_RKIND, wy))

     ! four corner indices (with periodic wrap in i)
     out_wt%src_i(1, n) = i0;                       out_wt%src_j(1, n) = j0
     out_wt%src_i(2, n) = mod(i0, field%nx) + 1;   out_wt%src_j(2, n) = j0
     out_wt%src_i(3, n) = i0;                       out_wt%src_j(3, n) = min(j0+1, field%ny)
     out_wt%src_i(4, n) = mod(i0, field%nx) + 1;   out_wt%src_j(4, n) = min(j0+1, field%ny)

     ! bilinear weights
     out_wt%weights(1, n) = (1.0_RKIND - wx) * (1.0_RKIND - wy)
     out_wt%weights(2, n) = wx               * (1.0_RKIND - wy)
     out_wt%weights(3, n) = (1.0_RKIND - wx) * wy
     out_wt%weights(4, n) = wx               * wy

  end do

end subroutine precompute_interp_weights

!==================================================================================================
 end module iamas_init_atm_emission
!==================================================================================================
