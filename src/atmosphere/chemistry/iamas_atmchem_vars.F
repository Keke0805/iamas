!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to declare or initialize variables needed by chemistry processes,and add vertical gas initial field.
! * Author: Jiawang Feng, Chen Jin, Zihan Xia and Chun Zhao
! Last update: Sep 2024
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao

! References:
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

 module iamas_atmchem_vars
 use iamas_kind_types
 use iamas_timekeeping
 use iamas_pool_routines
 use module_data_binaer_asect
 use iamas_atmchem_configs
 use iamas_atmphys_vars
 use iamas_constants , only: pii
 use iamas_log,only : iamas_log_write
 use iso_c_binding,only: c_int
 use flcl_mod
 implicit none
 public :: atmchem_vars_init
 save

!----------------------------------------------------------------------------!
!---------------------------Wrf vars-----------------------------------------!
!----------------------------------------------------------------------------!

 type(iAMAS_chem_config),public :: config_flags

! Whether to print detailed debug information.
 logical,public :: detailed_debug_on
! Number of different fields.
 integer,public :: num_species
 integer,public :: num_scalars
 integer,public :: num_chem
 integer,public :: num_scalar
 integer,public :: num_vert_mix
 integer,public :: num_gas
 integer,public :: num_phr !zhxia
!  integer,public :: nswbands, nlwbands

 integer,public :: ntot_amode_cam_mam = 3
 integer,public :: kemit
!  integer,public :: nbin_a

 real(kind=RKIND),public :: emis_dust_factor,dust_emass_jfrac,dust_emass_cfrac

 integer,public :: curr_year            !Current year.
 integer,public :: curr_month           !Current month.
 integer,public :: julday               !Initial Julian day.
 logical,public :: LeapYear             !Is leap year or not.
 logical,dimension(:),allocatable,public :: is_aerosol

 real(kind=RKIND),public :: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND),public :: gmt         !Greenwich mean time hour of model start (hr)
! Species indents used in wrf-chem.

! Used by interface module

 integer,pointer:: index_dust_a01,index_qv,index_qc,index_qi,index_qs
 integer,dimension(:),pointer:: ivgtyp,isltyp
 real(kind=RKIND),                 pointer :: len_disp
 real(kind=RKIND),dimension(:)    ,pointer :: meshDensity,u10,v10,ust,xland,dms_0
 real(kind=RKIND),dimension(:)    ,pointer :: rdzw,surface_pressure
 real(kind=RKIND),dimension(:)    ,pointer :: raincv,skintemp,gsw,vegfra,rmol
 real(kind=RKIND),dimension(:)    ,pointer :: snowh,snowc,hfx,hpbl,znt
 real(kind=RKIND),dimension(:,:)  ,pointer :: pressure,smois,zgrid,erod
 real(kind=RKIND),dimension(:,:)  ,pointer :: u_phy,v_phy
 real(kind=RKIND),dimension(:,:)  ,pointer :: rho_zz,zz,w,exner,pressure_p,pressure_b
 real(kind=RKIND),dimension(:,:)  ,pointer :: theta_m,kzh
 real(kind=RKIND),dimension(:,:)  ,pointer :: cldfrac,cldfrac_old
 real(kind=RKIND),dimension(:,:,:),pointer :: scalars,chem
 real(kind=RKIND),dimension(:,:)  ,pointer :: soilctop
 real(kind=RKIND),dimension(:,:)  ,pointer :: landusef
 real(kind=RKIND),dimension(:)    ,pointer :: lai
 real(kind=RKIND),dimension(:,:)  ,pointer :: sh2o

 real(kind=RKIND),dimension(:)    ,allocatable :: psfc,psfc_hyd,psfc_hydd,xlat,xlong
 real(kind=RKIND),dimension(:,:)  ,allocatable :: alt,pres_hydd,th,tke_pbl,viscosity,f_ice_local,f_rain_local
 real(kind=RKIND),dimension(:,:)  ,allocatable :: qv,qc,qi,qs,pres,znu_hyd,qndrop3d,qnice3d,nimey3d
 real(kind=RKIND),dimension(:,:)  ,allocatable :: t2,pres2,pres2_hydd ! nVertLevels + 1
 real(kind=RKIND),dimension(:,:,:),allocatable :: chem_old
 real(kind=RKIND),dimension(:,:,:),allocatable :: gas_aqfrac

 real(kind=RKIND),dimension(:)    ,allocatable     :: dx
 real(kind=RKIND),dimension(:,:)  ,allocatable     :: t_phy,rho_phy,dz,pres2_hyd,pres_hyd
 real(kind=RKIND),dimension(:,:,:),allocatable     :: moist

! type(dualview_r64_1d_t) :: dv_dx
! type(dualview_r64_2d_t) :: dv_t_phy,dv_chem,dv_rho_phy,dv_dz,dv_pres2_hyd,dv_pres_hyd
! type(dualview_r64_3d_t) :: dv_moist

! How to get p_p25???
 integer,parameter :: p_p25              = 1
 integer,public,parameter :: param_first_scalar = 1
 real   ,parameter :: pi                 = pii
 integer,public,pointer :: &
        p_qndrop, p_qniceaer, p_qnicemix                                 , &
        p_qv, p_qc,  p_qr, p_qi, p_qs                                    , &
        p_dust_a01,  p_dust_a02,  p_dust_a03,  p_dust_a04,  p_dust_a05   , &
        p_dust_a06,  p_dust_a07,  p_dust_a08,  p_dust_a09,  p_dust_a10   , &
        p_water_a01, p_water_a02, p_water_a03, p_water_a04, p_water_a05  , &
        p_water_a06, p_water_a07, p_water_a08, p_water_a09, p_water_a10  , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05    , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10    , &
        p_dust_cw01, p_dust_cw02, p_dust_cw03, p_dust_cw04, p_dust_cw05  , &
        p_dust_cw06, p_dust_cw07, p_dust_cw08, p_dust_cw09, p_dust_cw10  , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05   , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10   , &
        p_oc_a01,    p_oc_a02,    p_oc_a03,    p_oc_a04                  , &
        p_bc_a01,    p_bc_a02,    p_bc_a03,    p_bc_a04                  , &
        p_so4_a01,   p_so4_a02,   p_so4_a03,   p_so4_a04                 , &
        p_nh4_a01,   p_nh4_a02,   p_nh4_a03,   p_nh4_a04                 , &
        p_msa_a01,   p_msa_a02,   p_msa_a03,   p_msa_a04                 , &
        p_no3_a01,   p_no3_a02,   p_no3_a03,   p_no3_a04                 , &
        p_oin_a01,   p_oin_a02,   p_oin_a03,   p_oin_a04                 , &
        p_soa_a01,   p_soa_a02,   p_soa_a03,   p_soa_a04                 , &
        p_na_a01,    p_na_a02,    p_na_a03,    p_na_a04                  , &
        p_na_a05,    p_na_a06,    p_na_a07,    p_na_a08,    p_na_a09,    p_na_a10, &
        p_cl_a01,    p_cl_a02,    p_cl_a03,    p_cl_a04                  , &
        p_cl_a05,    p_cl_a06,    p_cl_a07,    p_cl_a08,    p_cl_a09,    p_cl_a10, &
        p_ca_a01,    p_ca_a02,    p_ca_a03,    p_ca_a04                  , &
        p_co3_a01,   p_co3_a02,   p_co3_a03,   p_co3_a04                 , &
        p_hysw_a01,  p_hysw_a02,  p_hysw_a03,  p_hysw_a04                , &
        p_oc_cw01,    p_oc_cw02,    p_oc_cw03,    p_oc_cw04              , &
        p_bc_cw01,    p_bc_cw02,    p_bc_cw03,    p_bc_cw04              , &
        p_so4_cw01,   p_so4_cw02,   p_so4_cw03,   p_so4_cw04             , &
        p_nh4_cw01,   p_nh4_cw02,   p_nh4_cw03,   p_nh4_cw04             , &
        p_msa_cw01,   p_msa_cw02,   p_msa_cw03,   p_msa_cw04             , &
        p_no3_cw01,   p_no3_cw02,   p_no3_cw03,   p_no3_cw04             , &
        p_oin_cw01,   p_oin_cw02,   p_oin_cw03,   p_oin_cw04             , &
        p_soa_cw01,   p_soa_cw02,   p_soa_cw03,   p_soa_cw04             , &
        p_na_cw01,    p_na_cw02,    p_na_cw03,    p_na_cw04              , &
        p_na_cw05,    p_na_cw06,    p_na_cw07,    p_na_cw08,    p_na_cw09,    p_na_cw10, &
        p_cl_cw01,    p_cl_cw02,    p_cl_cw03,    p_cl_cw04              , &
        p_cl_cw05,    p_cl_cw06,    p_cl_cw07,    p_cl_cw08,    p_cl_cw09,    p_cl_cw10, &
        p_ca_cw01,    p_ca_cw02,    p_ca_cw03,    p_ca_cw04              , &
        p_co3_cw01,   p_co3_cw02,   p_co3_cw03,   p_co3_cw04


! zhxia/chenjin
 integer,public,pointer :: &
  p_co2,       p_no,        p_no2,         p_so2                   , &
  p_nh3,       p_co,        p_o3,          p_h2o2                  , &
  p_no3,       p_n2o5,      p_hono,        p_hno3                  , &
  p_hno4,      p_h2so4,     p_hcho,        p_ccho                  , &
  p_rcho,      p_acet,      p_mek,         p_hcooh                 , &
  p_meoh,      p_etoh,      p_cco_oh,      p_rco_oh                , &
  p_gly,       p_mgly,      p_bacl,        p_cres                  , &
  p_bald,      p_isoprod,   p_methacro,    p_mvk                   , &
  p_prod2,     p_dcb1,      p_dcb2,        p_dcb3                  , &
  p_ethene,    p_isoprene,  p_c2h6,        p_c3h8                  , &
  p_c2h2,      p_c3h6,      p_alk3,        p_alk4                  , &
  p_alk5,      p_aro1,      p_aro2,        p_ole1                  , &
  p_ole2,      p_terp,      p_rno3,        p_nphe                  , &
  p_phen,      p_pan,       p_pan2,        p_pbzn                  , &
  p_ma_pan,    p_cco_ooh,   p_rco_o2,      p_rco_ooh               , &
  p_xn,        p_xc,        p_ho,          p_ho2                   , &
  p_c_o2,      p_cooh,      p_rooh,        p_ro2_r                 , &
  p_r2o2,      p_ro2_n,     p_cco_o2,      p_bzco_o2               , &
  p_ma_rco3,   p_sesq,      p_ch4,         p_hcl                   , &
  p_msa                                                            , &
  p_alk1,      p_alk2

! zhxia
 integer,public,pointer :: &
  p_ph_o2,       p_ph_o31d,     p_ph_o33p,      p_ph_no2           , &
  p_ph_no3o2,    p_ph_no3o,     p_ph_hno2,      p_ph_hno3          , &
  p_ph_hno4,     p_ph_h2o2,     p_ph_ch2or,     p_ph_ch2om         , &
  p_ph_ch3cho,   p_ph_ch3coch3, p_ph_ch3coc2h5, p_ph_hcocho        , &
  p_ph_ch3cocho, p_ph_hcochest, p_ph_ch3o2h,    p_ph_ch3coo2h      , &
  p_ph_ch3ono2,  p_ph_hcochob,  p_ph_n2o5

 integer,public,pointer :: &
  p_sulf,     p_cro    , & ! p_qr, p_hcl
  p_ch3o2, p_ethp,     p_ch3oh,   p_c2h5oh,    p_open,    p_op3    , &
  p_to2,   p_c2o3,     p_ro2,     p_ano2,      p_nap,     p_xo2    , &
  p_xpar,  p_isoprd,   p_isopp,   p_isopn,     p_isopo2            , & !p_dms 
  p_dmso,     p_dmso2,   p_ch3so2h,   p_ch3sch2oo,p_ch3so2         , & ! p_msa
  p_ch3so2oo,              p_ch3so2ch2oo,          p_xyl           , &
  p_mtf,   p_soag,     p_cvasoa1, p_cvasoa2,   p_cvasoa3, p_cvasoa4, &
  p_cvbsoa1,      p_op2,   p_eto2,  p_mbo,     p_ch3co2h           , &
  p_ald2,  p_iso,    p_macr,             p_aco3                    , & !p_alk1
  p_api1,  p_api2,     p_lim1,    p_lim2,      p_api,     p_lim    , &
  p_mpan,  p_mvkooh,   p_isooh,   p_xooh,      p_o,       p_o1d_cb4, &
  p_n2o,   p_h2,       p_ald,     p_c2h4                           , &
  p_eo,    p_pro2,     p_po2,     p_cvbsoa3,   p_aceto2,  p_bigene , &
  p_bigalk,p_eneo2,    p_alko2,   p_isopr,     p_iso2,    p_mvko2  , &
  p_c10h16,p_cvbsoa2,  p_terpo2,  p_tol,       p_cvbsoa4, p_xoh    , &
  p_meko2, p_voca,     p_vocbb,   p_smpa,      p_smpbb,   p_hso3   , &
  p_so3,               p_nc4h10,  p_c5h8,      p_benzene, p_toluene, &
  p_oxyl,  p_npropol,      p_tbut2ene,  p_ipropol, p_apinene       , &
  p_bpinene,p_ch3cl,   p_ch2cl2,  p_chcl3,     p_ch3ccl3, p_cdicleth, &
  p_tdicleth,p_tricleth,p_tce,    p_noa,       p_aroh14,  p_raroh14, &
  p_arnoh14,p_aroh17,  p_raroh17, p_arnoh17,   p_anhy,    p_tm123b , &
  p_tm124b,p_tm135b,   p_oethtol, p_methtol,   p_pethtol, p_dime35eb, &
  p_ch3s,  p_ch3so,    p_ch3so3,  p_hoch2ch2o2,p_ic3h7o2, p_rn10o2 , &
  p_rn13o2,p_rn16o2,   p_rn19o2,  p_rn9o2,     p_rn12o2,  p_rn15o2 , &
  p_rn18o2,p_nrn6o2,   p_nrn9o2,  p_nrn12o2,   p_rn11o2,  p_rn14o2 , &
  p_rn8o2, p_rn17o2,   p_rn13ao2, p_rn16ao2,   p_rn15ao2, p_rn18ao2, &
  p_ru14o2,p_ru12o2,   p_ru10o2,  p_nru14o2,   p_nru12o2, p_ra13o2 , &
  p_ra16o2,p_ra19ao2,  p_ra19co2, p_rtn28o2,   p_rtn26o2, p_nrtn28o2, &
  p_rtn25o2,p_rtn24o2, p_rtn23o2, p_rtn14o2,   p_rtn10o2, p_rtx28o2, &
  p_rtx24o2,p_rtx22o2, p_nrtx28o2,p_ch3o2no2,  p_ra22ao2, p_ra22bo2, &
  p_ra25o2,       p_ch3ooh,        p_c3h6ooh, p_etooh  , &
  p_prooh, p_aceth,    p_onit,    p_onitr   , &
  p_clno2,     p_msia,    p_ch3oo  , &
  p_c2h5o2,p_ch3cooh,  p_acetp,   p_acetol,    p_glyald,  p_hydrald, &
  p_alkooh,p_mekooh,   p_tolooh,  p_terpooh,   p_ket,     p_c2h5co3h, &
  p_hoch2co3h,p_ru12pan,p_rtn26pan,p_phan,     p_ppn,     p_c2h5ooh, &
  p_hoc2h4ooh,p_rn13ooh, p_rn10ooh,   p_rn16ooh, p_rn19ooh, &
  p_rn8ooh,p_rn11ooh,  p_rn14ooh, p_rn17ooh,   p_rn9ooh,  p_rn12ooh, &
  p_rn15ooh,p_rn18ooh, p_nrn6ooh, p_nrn12ooh,  p_ru14ooh, p_ru12ooh
!jwfeng from AERMOL USTC 10/2021
!Not included species, but will be included in future.
!Set to -1 for now, remember to change them into pointer vars.
 integer,public :: &
        p_ora1=-1,   p_op1=-1,    p_so4aj=-1,  p_nh4aj=-1  , p_dms=-1    , &
        p_no3aj=-1,  p_bc1=-1,    p_bc2=-1,    p_oc1=-1,     p_oc2=-1    , &
        p_seas_1=-1, p_seas_2=-1, p_seas_3=-1, p_seas_4=-1,  p_paa=-1


 contains
!=================================================================================================================
subroutine atmchem_vars_init(configs,mesh,state,diag_chemistry,tend_chemistry,clock)
!=================================================================================================================
  type(iamas_pool_type),intent(in) :: mesh
  type(iamas_pool_type),intent(in) :: configs
  type(iamas_pool_type),intent(in) :: state
  type(iamas_pool_type),intent(inout) :: diag_chemistry
  type(iamas_pool_type),intent(inout) :: tend_chemistry !zhxia

  integer,pointer :: nCellsSolve,nVertLevels,nSoils,nChem,nScalars
  integer,pointer :: species_start,species_end
  integer,pointer :: moist_start,moist_end
  integer,pointer :: gas_start,gas_end
  integer,pointer :: phr_start,phr_end !zhxia
  integer         :: nSpecies,nMoist,nPhr !zhxia
  real(kind=RKIND),pointer :: config_dt
  type(iAMAS_Time_Type) :: currTime
  type(iAMAS_Time_Type) :: startTime
  type(iAMAS_Clock_type),intent(in) :: clock

  character(len=StrKIND) :: timeStamp
  integer :: ierr
  real(kind=RKIND) :: utc_s,utc_h
  integer :: h,m,s,s_n,s_d,DoY,yr,month
  integer :: n,nhi,itype,isize
  real    :: dum

!local config variables
  integer,pointer :: gas_drydep_opt,aer_drydep_opt,kdepvel
!   integer,pointer :: nswbands_pointer, nlwbands_pointer
  character(len=StrKIND),pointer :: chem_opt,phot_opt,gaschem_opt,stra_gaschem_opt
  logical, pointer :: config_print_debug_info

!jwfeng: 2020/10/09 add for aerosol optical variables
!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw, absaersw, extaersw, gaersw, waersw, bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw, extaerlw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_dust, absaersw_dust, gaersw_dust, waersw_dust
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_seas, absaersw_seas, gaersw_seas, waersw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw_dust, tauaerlw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer

!zhxia: add for photolysis variables
!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: ph

!zhxia: add for linoz variables
!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: linoz_tlstt

!zhxia: add for radicals variables
!local pointers:
 real(kind=RKIND),dimension(:,:),pointer :: o1d
 real(kind=RKIND),dimension(:,:),pointer :: o3p
 real(kind=RKIND),dimension(:,:),pointer :: tbu_o
 real(kind=RKIND),dimension(:,:),pointer :: hocoo
 real(kind=RKIND),dimension(:,:),pointer :: bzno2_o
 real(kind=RKIND),dimension(:,:),pointer :: bz_o

 integer,pointer :: nCells
 real (kind=RKIND), dimension(:), pointer :: latCell
 real (kind=RKIND), dimension(:), pointer :: lonCell
 integer:: time_lev

  time_lev = 1 

  call iamas_pool_get_dimension(mesh,'nCells',nCells)
  call iamas_pool_get_dimension(mesh,'nVertLevels', nVertLevels)
  call iamas_pool_get_array(state,'chem', chem, time_lev)
  call iamas_pool_get_array(mesh,'latCell',latCell)
  call iamas_pool_get_array(mesh,'lonCell',lonCell)

  call iamas_pool_get_config(configs, 'config_print_debug_info', config_print_debug_info)
  call iamas_pool_get_config(configs  ,     'config_dt'     ,      config_dt      )
  call iamas_pool_get_dimension(state ,     'species_start' ,      species_start  )
  call iamas_pool_get_dimension(state ,     'species_end'   ,      species_end    )
  call iamas_pool_get_dimension(state ,     'moist_start'   ,      moist_start    )
  call iamas_pool_get_dimension(state ,     'moist_end'     ,      moist_end      )
  call iamas_pool_get_dimension(state ,     'gas_start'     ,      gas_start      )
  call iamas_pool_get_dimension(state ,     'gas_end'       ,      gas_end        )
  call iamas_pool_get_dimension(tend_chemistry ,     'phr_start'     ,     phr_start      ) !zhxia
  call iamas_pool_get_dimension(tend_chemistry ,     'phr_end'       ,     phr_end        ) !zhxia
  call iamas_pool_get_dimension(state ,     'num_chem'      ,      nChem          )
  call iamas_pool_get_dimension(state ,     'num_scalars'   ,      nScalars       )
  call iamas_pool_get_dimension(mesh  ,     'nSoilLevels'   ,      nSoils         )
  call iamas_pool_get_dimension(mesh  ,     'nCellsSolve'   ,      nCellsSolve    )

  call iamas_pool_get_config(configs,'config_chem_opt'        ,chem_opt           )
  call iamas_pool_get_config(configs,'config_gas_drydep_opt'  ,gas_drydep_opt     )
  call iamas_pool_get_config(configs,'config_aer_drydep_opt'  ,aer_drydep_opt     )
  call iamas_pool_get_config(configs,'config_kdepvel'         ,kdepvel            )
  call iamas_pool_get_config(configs,'config_phot_opt'        ,phot_opt) !zhxia
  call iamas_pool_get_config(configs,'config_stra_gaschem_opt',stra_gaschem_opt) !zhxia
  call iamas_pool_get_config(configs,'config_gaschem_opt',gaschem_opt) !zhxia
  
!jwfeng: 2020/10/09 add for aerosol optical variables
  call iamas_pool_get_array(diag_chemistry,'tauaersw'  , tauaersw )
  call iamas_pool_get_array(diag_chemistry,'absaersw'  , absaersw )
  call iamas_pool_get_array(diag_chemistry,'gaersw'    , gaersw   )
  call iamas_pool_get_array(diag_chemistry,'waersw'    , waersw   )
  call iamas_pool_get_array(diag_chemistry,'tauaersw_dust', tauaersw_dust )
  call iamas_pool_get_array(diag_chemistry,'absaersw_dust', absaersw_dust )
  call iamas_pool_get_array(diag_chemistry,'gaersw_dust'  , gaersw_dust   )
  call iamas_pool_get_array(diag_chemistry,'waersw_dust'  , waersw_dust   )
  call iamas_pool_get_array(diag_chemistry,'tauaersw_seas', tauaersw_seas )
  call iamas_pool_get_array(diag_chemistry,'absaersw_seas', absaersw_seas )
  call iamas_pool_get_array(diag_chemistry,'gaersw_seas'  , gaersw_seas   )
  call iamas_pool_get_array(diag_chemistry,'waersw_seas'  , waersw_seas   )
  call iamas_pool_get_array(diag_chemistry,'extaersw'  , extaersw )
  call iamas_pool_get_array(diag_chemistry,'bscoefsw'  , bscoefsw )
  call iamas_pool_get_array(diag_chemistry,'tauaerlw'  , tauaerlw )
  call iamas_pool_get_array(diag_chemistry,'tauaerlw_dust'  , tauaerlw_dust )
  call iamas_pool_get_array(diag_chemistry,'tauaerlw_seas'  , tauaerlw_seas )
  call iamas_pool_get_array(diag_chemistry,'extaerlw'  , extaerlw )
  call iamas_pool_get_array(diag_chemistry,'l2aer'     , l2aer    )
  call iamas_pool_get_array(diag_chemistry,'l3aer'     , l3aer    )
  call iamas_pool_get_array(diag_chemistry,'l4aer'     , l4aer    )
  call iamas_pool_get_array(diag_chemistry,'l5aer'     , l5aer    )
  call iamas_pool_get_array(diag_chemistry,'l6aer'     , l6aer    )
  call iamas_pool_get_array(diag_chemistry,'l7aer'     , l7aer    )
!zhxia: add for photolysis variables
  call iamas_pool_get_array(tend_chemistry,'ph'        , ph          )
!zhxia: add for linoz variables
  call iamas_pool_get_array(state         ,'linoz_tlstt'     , linoz_tlstt )
!zhxia: add for radicals variables
  call iamas_pool_get_array(state         ,'o1d'     , o1d ,time_lev)
  call iamas_pool_get_array(state         ,'o3p'     , o3p ,time_lev)
  call iamas_pool_get_array(state         ,'tbu_o'     , tbu_o ,time_lev)
  call iamas_pool_get_array(state         ,'hocoo'     , hocoo ,time_lev)
  call iamas_pool_get_array(state         ,'bzno2_o'     , bzno2_o ,time_lev)
  call iamas_pool_get_array(state         ,'bz_o'     , bz_o ,time_lev)

  config_flags%chem_opt         = chem_opt
  config_flags%chem_in_opt      = 1
  config_flags%gas_drydep_opt   = gas_drydep_opt
  config_flags%aer_drydep_opt   = aer_drydep_opt
!  config_flags%num_land_cat     = 24
  config_flags%kdepvel          = kdepvel
  config_flags%settling_option  = kdepvel

!!! All of these will be passed from Registry.xml
!!! and are belonged to phys in WRF-chem

! if(radt_sw_scheme == ...)
! jwfeng: Fixed a bug
  config_flags%ra_sw_physics    = RRTMG_SWSCHEME
! end if
  config_flags%bl_pbl_physics   = 1
  config_flags%wsubscheme       = 0
  config_flags%ice_para         = 1
  config_flags%mixcloud         = 0

  nSpecies = species_end - species_start + 1
  nMoist   =   moist_end -   moist_start + 1
  num_gas  =   gas_end   -     gas_start + 1
  nPhr = phr_end   -  phr_start + 1 !zhxia
  if (gas_start .eq. -1) num_gas = 0
!   nswbands = nswbands_pointer
!   nlwbands = nlwbands_pointer

! Get start gmt and julday:
  startTime   = iamas_get_clock_time(clock,iAMAS_START_TIME, ierr)
  call iamas_get_time(curr_time=startTime,YYYY=yr,MM=month,H=h,M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
  gmt         = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
  julday      = DoY

! Update the current julian day and current year:
  currTime    = iamas_get_clock_time(clock,iAMAS_NOW,ierr)
  call iamas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,MM=month,H=h,M=m,&
                     S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
  utc_h       = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
  utc_s       = real(s + s_n/s_d)
  curr_year   = yr
  curr_month  = month
  curr_julday = real(DoY-1) + utc_h / 24.0
  LeapYear    = isLeapYear(yr)

  num_species      = nSpecies
  num_chem         = nChem
  num_scalars      = nScalars
  num_scalar       = num_scalars
  num_moist        = nMoist ! num_moist is declared in iamas_atmphys_vars.F
  num_vert_mix     = num_gas
  num_phr          = nPhr !zhxia

  detailed_debug_on = config_print_debug_info

  if(.not.allocated(is_aerosol)) allocate(is_aerosol(1:num_chem))
  is_aerosol(:)    = .false.

  dust_emass_jfrac = 0.15
  dust_emass_cfrac = 0.85
!jwfeng: 2020/10/09 add for aerosol optical variables
  tauaersw(:,:,:) = 0.0
  absaersw(:,:,:) = 0.0
  extaersw(:,:,:) = 0.0
    gaersw(:,:,:) = 0.0
    waersw(:,:,:) = 0.0
  bscoefsw(:,:,:) = 0.0
  extaerlw(:,:,:) = 0.0
  tauaerlw(:,:,:) = 0.0
     l2aer(:,:,:) = 0.0
     l3aer(:,:,:) = 0.0
     l4aer(:,:,:) = 0.0
     l5aer(:,:,:) = 0.0
     l6aer(:,:,:) = 0.0
     l7aer(:,:,:) = 0.0
!zhxia: add for photolysis variables
  if (phot_opt .ne. 'no_phot') ph(:,:,:) = 0.0

!zhxia: add for linoz variables
  if (stra_gaschem_opt .ne. 'no_stra_gaschem') linoz_tlstt(:,:,:) = 0.0

!zhxia: add for radicals variables
  if (gaschem_opt .ne. 'no_gaschem') o1d(:,:) = 0.0
  if (gaschem_opt .ne. 'no_gaschem') o3p(:,:) = 0.0
  if (gaschem_opt .ne. 'no_gaschem') tbu_o(:,:) = 0.0
  if (gaschem_opt .ne. 'no_gaschem') hocoo(:,:) = 0.0
  if (gaschem_opt .ne. 'no_gaschem') bzno2_o(:,:) = 0.0
  if (gaschem_opt .ne. 'no_gaschem') bz_o(:,:) = 0.0

!init chem array.   updated on 2025/05/22. @zhxia
!Delete now. Because atmosphere will firstly read init_final.nc which has chem array.
  chem(:,:,:)= max(chem(:,:,:), 1.0e-16)

  call iamas_log_write( '' )
  call iamas_log_write( '--- Basic settings of chemistry species' )
  call iamas_log_write( '       chem_opt = '//trim(config_flags%chem_opt) )
  call iamas_log_write( '         nMoist = $i',intArgs = (/num_moist  /) )
  call iamas_log_write( '       nScalars = $i',intArgs = (/num_scalars/) )
  call iamas_log_write( '           nGas = $i',intArgs = (/num_gas    /) )
  call iamas_log_write( '       nSpecies = $i',intArgs = (/num_species/) )
  call iamas_log_write( '          nChem = $i',intArgs = (/num_chem   /) )
  call iamas_log_write( '      gas_start = $i',intArgs = (/gas_start/) )
  call iamas_log_write( '        gas_end = $i',intArgs = (/gas_end  /) )
  call iamas_log_write( '  species_start = $i',intArgs = (/species_start/) )
  call iamas_log_write( '    species_end = $i',intArgs = (/species_end   /) )
  call iamas_log_write( '         nCells = $i',intArgs = (/nCells   /) )
  call iamas_log_write( '    nCellsSolve = $i',intArgs = (/nCellsSolve   /) )
  call iamas_log_write( '   Size of chem = $i,$i,$i',intArgs = (/size(chem,dim=1), size(chem,dim=2), size(chem,dim=3)/) )
  call iamas_log_write( '' )

end subroutine atmchem_vars_init

 end module iamas_atmchem_vars
!=================================================================================================================
