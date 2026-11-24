!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to generate gas_init for iamas.
! * Author: Zihan Xia from LACAR USTC, xiazh@mail.ustc.edu.cn
! Last update: zhxia 06/2025
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
 module iamas_init_atm_wrfchem
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_sort
 use iamas_timekeeping
 use iamas_timer
 use iamas_log, only : iamas_log_write
 
 use init_atm_vinterp
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met

 implicit none
 private
 public :: init_atm_case_wrfchem
 public :: init_atm_case_ozone


 contains


!==================================================================================================
 subroutine init_atm_case_wrfchem(block, mesh, nCells, nEdges, nVertLevels, fg, state, diag, diag_physics, dims, configs)
!==================================================================================================

 use iamas_stream_manager

 implicit none

!input arguments:
 type (block_type), intent(inout), target :: block
 type (iamas_pool_type), intent(inout) :: mesh
 integer, intent(in) :: nCells
 integer, intent(in) :: nEdges
 integer, intent(in) :: nVertLevels
 type (iamas_pool_type), intent(inout) :: fg
 type (iamas_pool_type), intent(inout) :: state
 type (iamas_pool_type), intent(inout) :: diag
 type (iamas_pool_type), intent(inout):: diag_physics
 type (iamas_pool_type), intent(inout):: dims
 type (iamas_pool_type), intent(inout):: configs

 type (parallel_info), pointer :: parinfo
 type (dm_info), pointer :: dminfo

!local variables:
 character(len=StrKIND), pointer :: config_fix_gasinit_file
 integer, pointer :: config_nfggaslevels
 real (kind=RKIND), dimension(:,:,:), pointer :: chem
 integer,pointer ::                                                  &
     p_o3      , p_h2o2   , p_no     , p_no2     , p_no3   , &
     p_n2o5    , p_hono   , p_hno3   , p_hno4    , p_so2   , &
     p_h2so4   , p_co     , p_hcho   , p_ccho    , p_rcho  , &
     p_acet    , p_mek    , p_hcooh  , p_meoh    , p_etoh  , &
     p_cco_oh  , p_rco_oh , p_gly    , p_mgly    , p_bacl  , &
     p_cres    , p_bald   , p_isoprod, p_methacro, p_mvk   , &
     p_prod2   , p_dcb1   , p_dcb2   , p_dcb3    , p_ethene, &
     p_isoprene, p_c2h6   , p_c3h8   , p_c2h2    , p_c3h6  , &
     p_alk3    , p_alk4   , p_alk5   , p_aro1    , p_aro2  , &
     p_ole1    , p_ole2   , p_terp   , p_rno3    , p_nphe  , &
     p_phen    , p_pan    , p_pan2   , p_pbzn    , p_ma_pan, &
     p_co2     , p_cco_ooh, p_rco_o2 , p_rco_ooh , p_xn    , &
     p_xc      , p_ho     , p_ho2    , p_c_o2    , p_cooh  , &
     p_rooh    , p_ro2_r  , p_r2o2   , p_ro2_n   , p_cco_o2, &
     p_bzco_o2 , p_ma_rco3, p_sesq   , p_ch4
 real(kind=RKIND), dimension(:,:), pointer ::                     & 
     o3_fg      , h2o2_fg   , no_fg     , no2_fg     , no3_fg   , &
     n2o5_fg    , hono_fg   , hno3_fg   , hno4_fg    , so2_fg   , &
     h2so4_fg   , co_fg     , hcho_fg   , ccho_fg    , rcho_fg  , &
     acet_fg    , mek_fg    , hcooh_fg  , meoh_fg    , etoh_fg  , &
     cco_oh_fg  , rco_oh_fg , gly_fg    , mgly_fg    , bacl_fg  , &
     cres_fg    , bald_fg   , isoprod_fg, methacro_fg, mvk_fg   , &
     prod2_fg   , dcb1_fg   , dcb2_fg   , dcb3_fg    , ethene_fg, &
     isoprene_fg, c2h6_fg   , c3h8_fg   , c2h2_fg    , c3h6_fg  , &
     alk3_fg    , alk4_fg   , alk5_fg   , aro1_fg    , aro2_fg  , &
     ole1_fg    , ole2_fg   , terp_fg   , rno3_fg    , nphe_fg  , &
     phen_fg    , pan_fg    , pan2_fg   , pbzn_fg    , ma_pan_fg, &
     co2_fg     , cco_ooh_fg, rco_o2_fg , rco_ooh_fg , xn_fg    , &
     xc_fg      , ho_fg     , ho2_fg    , c_o2_fg    , cooh_fg  , &
     rooh_fg    , ro2_r_fg  , r2o2_fg   , ro2_n_fg   , cco_o2_fg, &
     bzco_o2_fg , ma_rco3_fg, sesq_fg   , ch4_fg

 real (kind=RKIND), dimension(:), pointer :: vert_level !model layer in wrfchem data
 real (kind=RKIND), dimension(:,:), pointer :: height
 real (kind=RKIND), dimension(:,:), pointer :: zgrid
 real (kind=RKIND) :: target_z
 type (iAMAS_Time_type)  :: curr_time
 character(len=StrKIND) :: timeString
 
 real (kind=RKIND), dimension(:,:), pointer :: destField2d
 real (kind=RKIND), dimension(:,:), pointer :: sorted_arr
 type (met_data) :: field
 integer, dimension(5) :: interp_list
 integer :: masked
 real(kind=RKIND) :: fillval, maskval, msgval
 integer, dimension(:), pointer :: mask_array 
 integer :: istatus
 integer ::iCell, k, l

 call iamas_pool_get_config(configs, 'config_fix_gasinit_file', config_fix_gasinit_file)
 call iamas_pool_get_config(configs, 'config_nfggaslevels', config_nfggaslevels)
 curr_time = iamas_get_clock_time(block % domain % clock, iAMAS_NOW)
 call iamas_get_time(curr_time, dateTimeString=timeString)
 call iamas_pool_get_array(mesh, 'zgrid', zgrid)

 call iamas_log_write('Adjust gas initial values to wrfchem simulation results.')

 call iamas_pool_get_array(mesh, 'landmask' , mask_array)
 call iamas_pool_get_array(state, 'chem'    , chem)
 call iamas_pool_get_dimension(state, 'index_o3'       , p_o3)
 call iamas_pool_get_dimension(state, 'index_h2o2 '    , p_h2o2)
 call iamas_pool_get_dimension(state, 'index_no'       , p_no)
 call iamas_pool_get_dimension(state, 'index_no2'      , p_no2)
 call iamas_pool_get_dimension(state, 'index_no3'      , p_no3)
 call iamas_pool_get_dimension(state, 'index_n2o5'     , p_n2o5)
 call iamas_pool_get_dimension(state, 'index_hono'     , p_hono)
 call iamas_pool_get_dimension(state, 'index_hno3'     , p_hno3)
 call iamas_pool_get_dimension(state, 'index_hno4'     , p_hno4)
 call iamas_pool_get_dimension(state, 'index_so2'      , p_so2)
 call iamas_pool_get_dimension(state, 'index_h2so4'    , p_h2so4)
 call iamas_pool_get_dimension(state, 'index_co'       , p_co)
 call iamas_pool_get_dimension(state, 'index_hcho'     , p_hcho)
 call iamas_pool_get_dimension(state, 'index_ccho'     , p_ccho)
 call iamas_pool_get_dimension(state, 'index_rcho'     , p_rcho)
 call iamas_pool_get_dimension(state, 'index_acet'     , p_acet)
 call iamas_pool_get_dimension(state, 'index_mek'      , p_mek)
 call iamas_pool_get_dimension(state, 'index_hcooh'    , p_hcooh)
 call iamas_pool_get_dimension(state, 'index_meoh'     , p_meoh)
 call iamas_pool_get_dimension(state, 'index_etoh'     , p_etoh)
 call iamas_pool_get_dimension(state, 'index_cco_oh'   , p_cco_oh)
 call iamas_pool_get_dimension(state, 'index_rco_oh'   , p_rco_oh)
 call iamas_pool_get_dimension(state, 'index_gly'      , p_gly)
 call iamas_pool_get_dimension(state, 'index_mgly'     , p_mgly)
 call iamas_pool_get_dimension(state, 'index_bacl'     , p_bacl)
 call iamas_pool_get_dimension(state, 'index_cres'     , p_cres)
 call iamas_pool_get_dimension(state, 'index_bald'     , p_bald)
 call iamas_pool_get_dimension(state, 'index_isoprod'  , p_isoprod)
 call iamas_pool_get_dimension(state, 'index_methacro' , p_methacro)
 call iamas_pool_get_dimension(state, 'index_mvk'      , p_mvk)
 call iamas_pool_get_dimension(state, 'index_prod2'    , p_prod2)
 call iamas_pool_get_dimension(state, 'index_dcb1'     , p_dcb1)
 call iamas_pool_get_dimension(state, 'index_dcb2'     , p_dcb2)
 call iamas_pool_get_dimension(state, 'index_dcb3'     , p_dcb3)
 call iamas_pool_get_dimension(state, 'index_ethene'   , p_ethene)
 call iamas_pool_get_dimension(state, 'index_isoprene' , p_isoprene)
 call iamas_pool_get_dimension(state, 'index_c2h6'     , p_c2h6)
 call iamas_pool_get_dimension(state, 'index_c3h8'     , p_c3h8)
 call iamas_pool_get_dimension(state, 'index_c2h2'     , p_c2h2 )
 call iamas_pool_get_dimension(state, 'index_c3h6'     , p_c3h6)
 call iamas_pool_get_dimension(state, 'index_alk3'     , p_alk3)
 call iamas_pool_get_dimension(state, 'index_alk4'     , p_alk4)
 call iamas_pool_get_dimension(state, 'index_alk5'     , p_alk5)
 call iamas_pool_get_dimension(state, 'index_aro1'     , p_aro1)
 call iamas_pool_get_dimension(state, 'index_aro2'     , p_aro2)
 call iamas_pool_get_dimension(state, 'index_ole1'     , p_ole1)
 call iamas_pool_get_dimension(state, 'index_ole2'     , p_ole2)
 call iamas_pool_get_dimension(state, 'index_terp'     , p_terp)
 call iamas_pool_get_dimension(state, 'index_rno3'     , p_rno3)
 call iamas_pool_get_dimension(state, 'index_nphe'     , p_nphe)
 call iamas_pool_get_dimension(state, 'index_phen'     , p_phen)
 call iamas_pool_get_dimension(state, 'index_pan'      , p_pan) 
 call iamas_pool_get_dimension(state, 'index_pan2'     , p_pan2 )
 call iamas_pool_get_dimension(state, 'index_pbzn'     , p_pbzn)
 call iamas_pool_get_dimension(state, 'index_ma_pan'   , p_ma_pan)
 call iamas_pool_get_dimension(state, 'index_co2'      , p_co2)
 call iamas_pool_get_dimension(state, 'index_cco_ooh'  , p_cco_ooh)
 call iamas_pool_get_dimension(state, 'index_rco_o2'   , p_rco_o2)
 call iamas_pool_get_dimension(state, 'index_rco_ooh'  , p_rco_ooh)
 call iamas_pool_get_dimension(state, 'index_xn'       , p_xn)
 call iamas_pool_get_dimension(state, 'index_xc'       , p_xc)
 call iamas_pool_get_dimension(state, 'index_ho'       , p_ho)
 call iamas_pool_get_dimension(state, 'index_ho2'      , p_ho2)
 call iamas_pool_get_dimension(state, 'index_c_o2'     , p_c_o2)
 call iamas_pool_get_dimension(state, 'index_cooh'     , p_cooh)
 call iamas_pool_get_dimension(state, 'index_rooh'     , p_rooh)
 call iamas_pool_get_dimension(state, 'index_ro2_r'    , p_ro2_r)
 call iamas_pool_get_dimension(state, 'index_r2o2'     , p_r2o2)
 call iamas_pool_get_dimension(state, 'index_ro2_n'    , p_ro2_n)
 call iamas_pool_get_dimension(state, 'index_cco_o2'   , p_cco_o2)
 call iamas_pool_get_dimension(state, 'index_bzco_o2'  , p_bzco_o2)
 call iamas_pool_get_dimension(state, 'index_ma_rco3'  , p_ma_rco3)
 call iamas_pool_get_dimension(state, 'index_sesq'     , p_sesq)
 call iamas_pool_get_dimension(state, 'index_ch4'      , p_ch4)

 allocate(vert_level(config_nfggaslevels))
 allocate(height(config_nfggaslevels,nCells))
 vert_level(:) = -1.0
 height(:,:) = 0.0

 !
 ! Horizontally interpolate chemical data
 !

 call iamas_log_write('Processing fix_gasinit file:'//trim(config_fix_gasinit_file))
 !read wrfchem simulation data. open intermediate file:
 call read_met_init(trim(config_fix_gasinit_file),.true.,timeString(1:13),istatus)
 call iamas_log_write('open gas simulation data')
 if(istatus /= 0) then
    call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
    call iamas_log_write('Error opening fix gasinit file '//trim(config_fix_gasinit_file)                  , messageType=iAMAS_LOG_ERR)
    call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
 end if
 
 !scan through all fields in the file, looking for the gas species field:
 call read_next_met_field(field,istatus)

 do while (istatus == 0)

  interp_list(1) = FOUR_POINT
  interp_list(2) = SEARCH
  interp_list(3) = 0

  msgval  = 0.0_RKIND !missing value
  masked  = -1             !no mask
  maskval = -1.0_RKIND
  fillval =  0.0_RKIND

  ! gas simulation data:
    do k=1,config_nfggaslevels
       if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
    end do
    if (vert_level(k) == -1.0) vert_level(k) = field % xlvl

    if(trim(field % field) == 'o3') then
        call iamas_log_write('Interpolating o3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'o3_fg'       , destField2d)
 
    else if(trim(field % field) == 'h2o2') then
        call iamas_log_write('Interpolating h2o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'h2o2_fg'    , destField2d)

    else if(trim(field % field) == 'no') then
        call iamas_log_write('Interpolating no at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'no_fg'       , destField2d)

    else if(trim(field % field) == 'no2') then
        call iamas_log_write('Interpolating no2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'no2_fg'    , destField2d)

    else if(trim(field % field) == 'no3') then
        call iamas_log_write('Interpolating no3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'no3_fg'       , destField2d)

    else if(trim(field % field) == 'n2o5') then
        call iamas_log_write('Interpolating n2o5 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'n2o5_fg'    , destField2d)

    else if(trim(field % field) == 'hono') then
        call iamas_log_write('Interpolating hono at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'hono_fg'       , destField2d)

    else if(trim(field % field) == 'hno3') then
        call iamas_log_write('Interpolating hno3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'hno3_fg'    , destField2d)

    else if(trim(field % field) == 'hno4') then
        call iamas_log_write('Interpolating hno4 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'hno4_fg'       , destField2d)

    else if(trim(field % field) == 'so2') then
        call iamas_log_write('Interpolating so2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'so2_fg'    , destField2d)

    else if(trim(field % field) == 'h2so4') then
        call iamas_log_write('Interpolating h2so4 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'h2so4_fg'       , destField2d)

    else if(trim(field % field) == 'co') then
        call iamas_log_write('Interpolating co at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'co_fg'    , destField2d)

    else if(trim(field % field) == 'hcho') then
        call iamas_log_write('Interpolating hcho at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'hcho_fg'       , destField2d)

    else if(trim(field % field) == 'ccho') then
        call iamas_log_write('Interpolating ccho at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ccho_fg'    , destField2d)

    else if(trim(field % field) == 'rcho') then
        call iamas_log_write('Interpolating rcho at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rcho_fg'       , destField2d)

    else if(trim(field % field) == 'acet') then
        call iamas_log_write('Interpolating acet at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'acet_fg'    , destField2d)

    else if(trim(field % field) == 'mek') then
        call iamas_log_write('Interpolating mek at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'mek_fg'       , destField2d)

    else if(trim(field % field) == 'hcooh') then
        call iamas_log_write('Interpolating hcooh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'hcooh_fg'    , destField2d)

    else if(trim(field % field) == 'meoh') then
        call iamas_log_write('Interpolating meoh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'meoh_fg'       , destField2d)

    else if(trim(field % field) == 'etoh') then
        call iamas_log_write('Interpolating etoh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'etoh_fg'    , destField2d)

    else if(trim(field % field) == 'cco_oh') then
        call iamas_log_write('Interpolating cco_oh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'cco_oh_fg'       , destField2d)

    else if(trim(field % field) == 'rco_oh') then
        call iamas_log_write('Interpolating rco_oh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rco_oh_fg'    , destField2d)

    else if(trim(field % field) == 'gly') then
        call iamas_log_write('Interpolating gly at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'gly_fg'       , destField2d)

    else if(trim(field % field) == 'mgly') then
        call iamas_log_write('Interpolating mgly at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'mgly_fg'    , destField2d)

    else if(trim(field % field) == 'bacl') then
        call iamas_log_write('Interpolating bacl at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'bacl_fg'       , destField2d)

    else if(trim(field % field) == 'cres') then
        call iamas_log_write('Interpolating cres at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'cres_fg'    , destField2d)

    else if(trim(field % field) == 'bald') then
        call iamas_log_write('Interpolating bald at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'bald_fg'       , destField2d)

    else if(trim(field % field) == 'isoprod') then
        call iamas_log_write('Interpolating isoprod at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'isoprod_fg'    , destField2d)

    else if(trim(field % field) == 'methacro') then
        call iamas_log_write('Interpolating methacro at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'methacro_fg'       , destField2d)

    else if(trim(field % field) == 'mvk') then
        call iamas_log_write('Interpolating mvk at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'mvk_fg'    , destField2d)

    else if(trim(field % field) == 'prod2') then
        call iamas_log_write('Interpolating prod2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'prod2_fg'       , destField2d)

    else if(trim(field % field) == 'dcb1') then
        call iamas_log_write('Interpolating dcb1 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'dcb1_fg'    , destField2d)

    else if(trim(field % field) == 'dcb2') then
        call iamas_log_write('Interpolating dcb2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'dcb2_fg'       , destField2d)

    else if(trim(field % field) == 'dcb3') then
        call iamas_log_write('Interpolating dcb3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'dcb3_fg'       , destField2d)

    else if(trim(field % field) == 'ethene') then
        call iamas_log_write('Interpolating ethene at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ethene_fg'    , destField2d)

    else if(trim(field % field) == 'isoprene') then
        call iamas_log_write('Interpolating isoprene at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'isoprene_fg'       , destField2d)

    else if(trim(field % field) == 'c2h6') then
        call iamas_log_write('Interpolating c2h6 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'c2h6_fg'    , destField2d)

    else if(trim(field % field) == 'c3h8') then
        call iamas_log_write('Interpolating c3h8 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'c3h8_fg'       , destField2d)

    else if(trim(field % field) == 'c2h2') then
        call iamas_log_write('Interpolating c2h2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'c2h2_fg'    , destField2d)

    else if(trim(field % field) == 'c3h6') then
        call iamas_log_write('Interpolating c3h6 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'c3h6_fg'       , destField2d)

    else if(trim(field % field) == 'alk3') then
        call iamas_log_write('Interpolating alk3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'alk3_fg'    , destField2d)

    else if(trim(field % field) == 'alk4') then
        call iamas_log_write('Interpolating alk4 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'alk4_fg'       , destField2d)

    else if(trim(field % field) == 'alk5') then
        call iamas_log_write('Interpolating alk5 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'alk5_fg'    , destField2d)

    else if(trim(field % field) == 'aro1') then
        call iamas_log_write('Interpolating aro1 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'aro1_fg'       , destField2d)

    else if(trim(field % field) == 'aro2') then
        call iamas_log_write('Interpolating aro2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'aro2_fg'    , destField2d)

    else if(trim(field % field) == 'ole1') then
        call iamas_log_write('Interpolating ole1 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ole1_fg'       , destField2d)

    else if(trim(field % field) == 'ole2') then
        call iamas_log_write('Interpolating ole2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ole2_fg'    , destField2d)

    else if(trim(field % field) == 'terp') then
        call iamas_log_write('Interpolating terp at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'terp_fg'       , destField2d)

    else if(trim(field % field) == 'rno3') then
        call iamas_log_write('Interpolating rno3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rno3_fg'    , destField2d)

    else if(trim(field % field) == 'nphe') then
        call iamas_log_write('Interpolating nphe at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'nphe_fg'       , destField2d)

    else if(trim(field % field) == 'phen') then
        call iamas_log_write('Interpolating phen at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'phen_fg'    , destField2d)

    else if(trim(field % field) == 'pan') then
        call iamas_log_write('Interpolating pan at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'pan_fg'       , destField2d)

    else if(trim(field % field) == 'pan2') then
        call iamas_log_write('Interpolating pan2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'pan2_fg'    , destField2d)

    else if(trim(field % field) == 'pbzn') then
        call iamas_log_write('Interpolating pbzn at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'pbzn_fg'       , destField2d)

    else if(trim(field % field) == 'ma_pan') then
        call iamas_log_write('Interpolating ma_pan at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ma_pan_fg'    , destField2d)

    else if(trim(field % field) == 'co2') then
        call iamas_log_write('Interpolating co2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'co2_fg'       , destField2d)

    else if(trim(field % field) == 'cco_ooh') then
        call iamas_log_write('Interpolating cco_ooh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'cco_ooh_fg'    , destField2d)

    else if(trim(field % field) == 'rco_o2') then
        call iamas_log_write('Interpolating rco_o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rco_o2_fg'       , destField2d)

    else if(trim(field % field) == 'rco_ooh') then
        call iamas_log_write('Interpolating rco_ooh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rco_ooh_fg'    , destField2d)

    else if(trim(field % field) == 'xn') then
        call iamas_log_write('Interpolating xn at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'xn_fg'       , destField2d)

    else if(trim(field % field) == 'xc') then
        call iamas_log_write('Interpolating xc at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'xc_fg'    , destField2d)

    else if(trim(field % field) == 'ho') then
        call iamas_log_write('Interpolating ho at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ho_fg'       , destField2d)

    else if(trim(field % field) == 'ho2') then
        call iamas_log_write('Interpolating ho2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ho2_fg'       , destField2d)

    else if(trim(field % field) == 'c_o2') then
        call iamas_log_write('Interpolating c_o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'c_o2_fg'    , destField2d)

    else if(trim(field % field) == 'cooh') then
        call iamas_log_write('Interpolating cooh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'cooh_fg'       , destField2d)

    else if(trim(field % field) == 'rooh') then
        call iamas_log_write('Interpolating rooh at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'rooh_fg'    , destField2d)

    else if(trim(field % field) == 'ro2_r') then
        call iamas_log_write('Interpolating ro2_r at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ro2_r_fg'       , destField2d)

    else if(trim(field % field) == 'r2o2') then
        call iamas_log_write('Interpolating r2o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'r2o2_fg'    , destField2d)

    else if(trim(field % field) == 'ro2_n') then
        call iamas_log_write('Interpolating ro2_n at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ro2_n_fg'       , destField2d)

    else if(trim(field % field) == 'cco_o2') then
        call iamas_log_write('Interpolating cco_o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'cco_o2_fg'    , destField2d)

    else if(trim(field % field) == 'bzco_o2') then
        call iamas_log_write('Interpolating bzco_o2 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'bzco_o2_fg'       , destField2d)

    else if(trim(field % field) == 'ma_rco3') then
        call iamas_log_write('Interpolating ma_rco3 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ma_rco3_fg'       , destField2d)

    else if(trim(field % field) == 'sesq') then
        call iamas_log_write('Interpolating sesq at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'sesq_fg'    , destField2d)

    else if (trim(field % field) == 'ch4') then
        call iamas_log_write('Interpolating ch4 at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call iamas_pool_get_array(fg, 'ch4_fg'    , destField2d)

    else if (trim(field % field) == 'height') then
        call iamas_log_write('Interpolating height at $i wrfchem model layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        destField2d => height

    end if

    call interp_to_iAMAS(mesh,nCells,k,field,destField2d,interp_list,msgval,masked,maskval,fillval,mask_array)

    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do

 !close intermediate file 
 call read_met_close()

 !
 ! Vertically interpolate chemical  data
 !

 call iamas_log_write('Vertically interpolating all gas initial field.')
 
 call iamas_pool_get_array(fg, 'o3_fg'       , o3_fg)
 call iamas_pool_get_array(fg, 'h2o2_fg'    , h2o2_fg)
 call iamas_pool_get_array(fg, 'no_fg'       , no_fg)
 call iamas_pool_get_array(fg, 'no2_fg'      , no2_fg)
 call iamas_pool_get_array(fg, 'no3_fg'      , no3_fg)
 call iamas_pool_get_array(fg, 'n2o5_fg'     , n2o5_fg)
 call iamas_pool_get_array(fg, 'hono_fg'     , hono_fg)
 call iamas_pool_get_array(fg, 'hno3_fg'     , hno3_fg)
 call iamas_pool_get_array(fg, 'hno4_fg'     , hno4_fg)
 call iamas_pool_get_array(fg, 'so2_fg'      , so2_fg)
 call iamas_pool_get_array(fg, 'h2so4_fg'    , h2so4_fg)
 call iamas_pool_get_array(fg, 'co_fg'       , co_fg)
 call iamas_pool_get_array(fg, 'hcho_fg'     , hcho_fg)
 call iamas_pool_get_array(fg, 'ccho_fg'     , ccho_fg)
 call iamas_pool_get_array(fg, 'rcho_fg'     , rcho_fg)
 call iamas_pool_get_array(fg, 'acet_fg'     , acet_fg)
 call iamas_pool_get_array(fg, 'mek_fg'      , mek_fg)
 call iamas_pool_get_array(fg, 'hcooh_fg'    , hcooh_fg)
 call iamas_pool_get_array(fg, 'meoh_fg'     , meoh_fg)
 call iamas_pool_get_array(fg, 'etoh_fg'     , etoh_fg)
 call iamas_pool_get_array(fg, 'cco_oh_fg'   , cco_oh_fg)
 call iamas_pool_get_array(fg, 'rco_oh_fg'   , rco_oh_fg)
 call iamas_pool_get_array(fg, 'gly_fg'      , gly_fg)
 call iamas_pool_get_array(fg, 'mgly_fg'     , mgly_fg)
 call iamas_pool_get_array(fg, 'bacl_fg'     , bacl_fg)
 call iamas_pool_get_array(fg, 'cres_fg'     , cres_fg)
 call iamas_pool_get_array(fg, 'bald_fg'     , bald_fg)
 call iamas_pool_get_array(fg, 'isoprod_fg'  , isoprod_fg)
 call iamas_pool_get_array(fg, 'methacro_fg' , methacro_fg)
 call iamas_pool_get_array(fg, 'mvk_fg'      , mvk_fg)
 call iamas_pool_get_array(fg, 'prod2_fg'    , prod2_fg)
 call iamas_pool_get_array(fg, 'dcb1_fg'     , dcb1_fg)
 call iamas_pool_get_array(fg, 'dcb2_fg'     , dcb2_fg)
 call iamas_pool_get_array(fg, 'dcb3_fg'     , dcb3_fg)
 call iamas_pool_get_array(fg, 'ethene_fg'   , ethene_fg)
 call iamas_pool_get_array(fg, 'isoprene_fg' , isoprene_fg)
 call iamas_pool_get_array(fg, 'c2h6_fg'     , c2h6_fg)
 call iamas_pool_get_array(fg, 'c3h8_fg'     , c3h8_fg)
 call iamas_pool_get_array(fg, 'c2h2_fg'     , c2h2_fg)
 call iamas_pool_get_array(fg, 'c3h6_fg'     , c3h6_fg)
 call iamas_pool_get_array(fg, 'alk3_fg'     , alk3_fg)
 call iamas_pool_get_array(fg, 'alk4_fg'     , alk4_fg)
 call iamas_pool_get_array(fg, 'alk5_fg'     , alk5_fg)
 call iamas_pool_get_array(fg, 'aro1_fg'     , aro1_fg)
 call iamas_pool_get_array(fg, 'aro2_fg'     , aro2_fg)
 call iamas_pool_get_array(fg, 'ole1_fg'     , ole1_fg)
 call iamas_pool_get_array(fg, 'ole2_fg'     , ole2_fg)
 call iamas_pool_get_array(fg, 'terp_fg'     , terp_fg)
 call iamas_pool_get_array(fg, 'rno3_fg'     , rno3_fg)
 call iamas_pool_get_array(fg, 'nphe_fg'     , nphe_fg)
 call iamas_pool_get_array(fg, 'phen_fg'     , phen_fg)
 call iamas_pool_get_array(fg, 'pan_fg'      , pan_fg) 
 call iamas_pool_get_array(fg, 'pan2_fg'     , pan2_fg)
 call iamas_pool_get_array(fg, 'pbzn_fg'     , pbzn_fg)
 call iamas_pool_get_array(fg, 'ma_pan_fg'   , ma_pan_fg)
 call iamas_pool_get_array(fg, 'co2_fg'      , co2_fg)
 call iamas_pool_get_array(fg, 'cco_ooh_fg'  , cco_ooh_fg)
 call iamas_pool_get_array(fg, 'rco_o2_fg'   , rco_o2_fg)
 call iamas_pool_get_array(fg, 'rco_ooh_fg'  , rco_ooh_fg)
 call iamas_pool_get_array(fg, 'xn_fg'       , xn_fg)
 call iamas_pool_get_array(fg, 'xc_fg'       , xc_fg)
 call iamas_pool_get_array(fg, 'ho_fg'       , ho_fg)
 call iamas_pool_get_array(fg, 'ho2_fg'      , ho2_fg)
 call iamas_pool_get_array(fg, 'c_o2_fg'     , c_o2_fg)
 call iamas_pool_get_array(fg, 'cooh_fg'     , cooh_fg)
 call iamas_pool_get_array(fg, 'rooh_fg'     , rooh_fg)
 call iamas_pool_get_array(fg, 'ro2_r_fg'    , ro2_r_fg)
 call iamas_pool_get_array(fg, 'r2o2_fg'     , r2o2_fg)
 call iamas_pool_get_array(fg, 'ro2_n_fg'    , ro2_n_fg)
 call iamas_pool_get_array(fg, 'cco_o2_fg'   , cco_o2_fg)
 call iamas_pool_get_array(fg, 'bzco_o2_fg'  , bzco_o2_fg)
 call iamas_pool_get_array(fg, 'ma_rco3_fg'  , ma_rco3_fg)
 call iamas_pool_get_array(fg, 'sesq_fg'     , sesq_fg)
 call iamas_pool_get_array(fg, 'ch4_fg'      , ch4_fg)

 allocate(sorted_arr(2,config_nfggaslevels))  

 do iCell=1,nCells

    do k = 1, nVertLevels
       target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))

       ! o3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = o3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_o3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! h2o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = h2o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_h2o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! no
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = no_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_no,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! no2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = no2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_no2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! no3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = no3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_no3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! n2o5
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = n2o5_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_n2o5,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! hono
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = hono_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_hono,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! hno3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = hno3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_hno3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! hno4
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = hno4_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_hno4,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! so2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = so2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_so2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! h2so4
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = h2so4_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_h2so4,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! co
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = co_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_co,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! hcho
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = hcho_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_hcho,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ccho
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ccho_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ccho,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rcho
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rcho_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rcho,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! acet
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = acet_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_acet,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! mek
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = mek_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_mek,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! hcooh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = hcooh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_hcooh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! meoh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = meoh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_meoh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! etoh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = etoh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_etoh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! cco_oh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = cco_oh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_cco_oh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rco_oh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rco_oh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rco_oh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! gly
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = gly_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_gly,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! mgly
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = mgly_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_mgly,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! bacl
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = bacl_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_bacl,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! cres
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = cres_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_cres,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! bald
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = bald_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_bald,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! isoprod
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = isoprod_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_isoprod,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! methacro
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = methacro_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_methacro,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! mvk
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = mvk_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_mvk,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! prod2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = prod2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_prod2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! dcb1
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = dcb1_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_dcb1,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! dcb2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = dcb2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_dcb2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! dcb3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = dcb3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_dcb3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ethene
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ethene_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ethene,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! isoprene
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = isoprene_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_isoprene,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! c2h6
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = c2h6_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_c2h6,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! c3h8
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = c3h8_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_c3h8,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! c2h2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = c2h2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_c2h2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! c3h6
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = c3h6_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_c3h6,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! alk3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = alk3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_alk3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! alk4
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = alk4_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_alk4,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! alk5
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = alk5_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_alk5,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! aro1
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = aro1_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_aro1,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! aro2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = aro2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_aro2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ole1
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ole1_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ole1,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ole2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ole2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ole2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! terp
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = terp_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_terp,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rno3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rno3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rno3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! nphe
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = nphe_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_nphe,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! phen
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = phen_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_phen,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! pan
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = pan_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_pan,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! pan2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = pan2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_pan2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! pbzn
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = pbzn_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_pbzn,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ma_pan
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ma_pan_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ma_pan,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! co2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = co2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_co2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! cco_ooh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = cco_ooh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_cco_ooh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rco_o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rco_o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rco_o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rco_ooh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rco_ooh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rco_ooh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! xn
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = xn_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_xn,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! xc
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = xc_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_xc,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ho
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ho_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ho,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ho2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ho2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ho2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! c_o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = c_o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_c_o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! cooh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = cooh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_cooh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! rooh
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = rooh_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_rooh,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ro2_r
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ro2_r_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ro2_r,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! r2o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = r2o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_r2o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ro2_n
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ro2_n_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ro2_n,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! cco_o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = cco_o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_cco_o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! bzco_o2
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = bzco_o2_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_bzco_o2,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ma_rco3
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ma_rco3_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ma_rco3,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! sesq
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = sesq_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_sesq,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, & 
                                     extrap=0) !extrap=0 represents no extrapolation.

       ! ch4
       sorted_arr(:,:) = -999.0
       do l = 1, config_nfggaslevels
          sorted_arr(1,l) = height(l,iCell)
          sorted_arr(2,l) = ch4_fg(l,iCell)
       enddo
       call iamas_quicksort(config_nfggaslevels, sorted_arr)
       chem(p_ch4,k,iCell) = vertical_interp(target_z, config_nfggaslevels-1, &
                                     sorted_arr(:,1:config_nfggaslevels-1),order=1, &         
                                     extrap=0) !extrap=0 represents no extrapolation.

    enddo !nVert

 enddo !nCell
 deallocate(sorted_arr)

 end subroutine init_atm_case_wrfchem

!==================================================================================================
 subroutine init_atm_case_ozone(block, mesh, nCells, nEdges, nVertLevels, state, diag, diag_physics, dims, configs)
!==================================================================================================

 use iamas_stream_manager

 implicit none

!input arguments:
 type (block_type), intent(inout), target :: block
 type (iamas_pool_type), intent(inout) :: mesh
 integer, intent(in) :: nCells
 integer, intent(in) :: nEdges
 integer, intent(in) :: nVertLevels
 type (iamas_pool_type), intent(inout) :: state
 type (iamas_pool_type), intent(inout) :: diag
 type (iamas_pool_type), intent(inout):: diag_physics
 type (iamas_pool_type), intent(inout):: dims
 type (iamas_pool_type), intent(inout):: configs

 type (parallel_info), pointer :: parinfo
 type (dm_info), pointer :: dminfo

!local variables:
 character(len=StrKIND), pointer :: config_fix_ozone_init_file
 logical, pointer :: config_use_fix_ozone_init !determine whether we use the ozone dataset from TOST as initial values in iAMAS.
 real (kind=RKIND), dimension(:,:,:), pointer :: chem
 integer,pointer :: p_o3

 integer :: nmonth=12,nlev=26,nlat=36,nlon=72 !data dims in ozone dataset from TOST
 real (kind=RKIND), dimension(:), allocatable :: vert_level !vertical layer in ozone dataset from TOST
 real (kind=RKIND), dimension(:,:), allocatable :: height
 real (kind=RKIND), dimension(:,:), pointer :: zgrid
 real (kind=RKIND) :: target_z
 type (iAMAS_Time_type)  :: curr_time
 character(len=StrKIND) :: timeString
 
 real (kind=RKIND), dimension(:,:), allocatable :: destField2d
 real (kind=RKIND), dimension(:,:), allocatable :: sorted_arr
 type (met_data) :: field
 integer, dimension(5) :: interp_list
 integer :: masked
 real(kind=RKIND) :: fillval, maskval, msgval
 integer, dimension(:), pointer :: mask_array 
 integer :: istatus
 integer ::iCell, k, l

 call iamas_pool_get_config(configs, 'config_fix_ozone_init_file', config_fix_ozone_init_file)
 curr_time = iamas_get_clock_time(block % domain % clock, iAMAS_NOW)
 call iamas_get_time(curr_time, dateTimeString=timeString)
 call iamas_pool_get_array(mesh, 'zgrid', zgrid)

 call iamas_log_write('Adjust ozone initial values to ozone dataset from TOST.')

 call iamas_pool_get_array(mesh, 'landmask' , mask_array)
 call iamas_pool_get_array(state, 'chem'    , chem)
 call iamas_pool_get_dimension(state, 'index_o3'       , p_o3)
 
 allocate(vert_level(nlev))
 allocate(height(nlev,nCells))
 allocate(destField2d(nlev,nCells))
 vert_level(:) = -1.0
 height(:,:) = 0.0
 destField2d(:,:) = 0.0

 !
 ! Horizontally interpolate chemical data
 !

 call iamas_log_write('Processing fix_ozone init file:'//trim(config_fix_ozone_init_file))
 !read wrfchem simulation data. open intermediate file:
 call read_met_init(trim(config_fix_ozone_init_file),.true.,timeString(1:13),istatus)
 call iamas_log_write('open ozone clim data')
 if(istatus /= 0) then
    call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_ERR)
    call iamas_log_write('Error opening fix ozone init file '//trim(config_fix_ozone_init_file)                  , messageType=iAMAS_LOG_ERR)
    call iamas_log_write('********************************************************************************', messageType=iAMAS_LOG_CRIT)
 end if
 
 !scan through all fields in the file, looking fo ozone field:
 call read_next_met_field(field,istatus)

 do while (istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = SEARCH
    interp_list(3) = 0

    msgval  = 0.0_RKIND !missing value
    masked  = -1             !no mask
    maskval = -1.0_RKIND
    fillval =  0.0_RKIND

    ! gas simulation data:
    do k=1,nlev
       if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
    end do
    if (vert_level(k) == -1.0) vert_level(k) = field % xlvl

    if(trim(field % field) == 'o3') then
        call iamas_log_write('Interpolating o3 at $i vertical layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call interp_to_iAMAS(mesh,nCells,k,field,destField2d,interp_list,msgval,masked,maskval,fillval,mask_array)
    else if (trim(field % field) == 'height') then
        call iamas_log_write('Interpolating height at $i vertical layer $r', intArgs=(/k/), realArgs=(/vert_level(k)/))
        call interp_to_iAMAS(mesh,nCells,k,field,height,interp_list,msgval,masked,maskval,fillval,mask_array)
    end if

    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do

 !close intermediate file
 call read_met_close()

 !
 ! Vertically interpolate chemical  data
 !

 call iamas_log_write('Vertically interpolating ozone initial field.')

 allocate(sorted_arr(2,nlev))
 
 do iCell=1,nCells

     do k = 1, nVertLevels
        target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))

        ! o3
        sorted_arr(:,:) = -999.0
        do l = 1, nlev
           sorted_arr(1,l) = height(l,iCell)
           sorted_arr(2,l) = destField2d(l,iCell)
        enddo
        call iamas_quicksort(nlev, sorted_arr)
        chem(p_o3,k,iCell) = vertical_interp(target_z, nlev, &
                                      sorted_arr(:,1:nlev),order=1, &
                                      extrap=0) !extrap=0 represents no extrapolation.
     enddo !nVert
 enddo !nCell
 deallocate(sorted_arr)


 deallocate(vert_level)
 deallocate(height)
 deallocate(destField2d)

  end subroutine init_atm_case_ozone

!==================================================================================================
 subroutine interp_to_iAMAS(mesh,nCells,k,field,destField2d,interp_list,msgval,masked,maskval,fillval, &
                           mask_array)
!==================================================================================================

!input arguments:
 type (iamas_pool_type), intent(in) :: mesh
 integer, intent(in) :: nCells, k
 type (met_data), intent(in)  :: field !real*4 meteorological data.

 integer, intent(in) :: masked
 integer, dimension(5), intent(in) :: interp_list
 integer, dimension(:), intent(in), pointer :: mask_array

 real(kind=RKIND), intent(in) :: fillval, maskval, msgval

!inout arguments:
 real(kind=RKIND), intent(inout), dimension(:,:) :: destField2d

!local variables:
 type(proj_info) :: proj
 integer :: i, nInterpPoints
 real(kind=RKIND) :: lat,lon,x,y
 real(kind=RKIND), dimension(:,:), allocatable :: rslab

 real(kind=RKIND), dimension(:), pointer :: latPoints, lonPoints
 real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
 
!--------------------------------------------------------------------------------------------------

 call iamas_pool_get_array(mesh, 'latCell', latCell)
 call iamas_pool_get_array(mesh, 'lonCell', lonCell)

 call map_init(proj)   
 if(field % iproj == PROJ_LATLON) then
    call map_set(PROJ_LATLON, proj, &
                 latinc = real(field % deltalat,RKIND), &
                 loninc = real(field % deltalon,RKIND), &
                 knowni = 1.0_RKIND, &
                 knownj = 1.0_RKIND, &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   call iamas_log_write('--- The projection is PROJ_LATLON.')
 else if(field % iproj == PROJ_GAUSS) then
    call map_set(PROJ_GAUSS, proj, &
                 nlat = nint(field % deltalat), &
                 loninc = 360.0_RKIND / real(field % nx,RKIND), &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   call iamas_log_write('--- The projection is PROJ_GAUSS.')
 else if(field % iproj == PROJ_PS) then
    call map_set(PROJ_PS, proj, &
                 dx = real(field % dx,RKIND), &
                 truelat1 = real(field % truelat1,RKIND), &
                 stdlon = real(field % xlonc,RKIND), &
                 knowni = real(field % nx / 2.0,RKIND), &
                 knownj = real(field % ny / 2.0,RKIND), &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   call iamas_log_write('--- The projection is PROJ_PS.')
 end if

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

       destField2d(k,i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                        msgval,interp_list,1,maskval=maskval)
    else
       destField2d(k,i) = fillval
    end if
 end do
 deallocate(rslab)

 end subroutine interp_to_iAMAS

!==================================================================================================
 end module iamas_init_atm_wrfchem
!==================================================================================================

