MODULE module_sf_noahdrv
!-------------------------------
  USE module_sf_noahlsm,  only: SFLX, XLF, XLV, CP, R_D, RHOWATER, NATURAL, SHDTBL, LUTYPE, SLTYPE, STBOLT, &
      &                         KARMAN, LUCATS, NROTBL, RSTBL, RGLTBL, HSTBL, SNUPTBL, MAXALB, LAIMINTBL,   &
      &                         LAIMAXTBL, Z0MINTBL, Z0MAXTBL, ALBEDOMINTBL, ALBEDOMAXTBL, EMISSMINTBL,     &
      &                         EMISSMAXTBL, TOPT_DATA, CMCMAX_DATA, CFACTR_DATA, RSMAX_DATA, BARE, NLUS,   &
      &                         SLCATS, BB, DRYSMC, F11, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, WLTSMC, QTZ, &
      &                         NSLTYPE, SLPCATS, SLOPE_DATA, SBETA_DATA, FXEXP_DATA, CSOIL_DATA,           &
      &                         SALP_DATA, REFDK_DATA, REFKDT_DATA, FRZK_DATA, ZBOT_DATA, CZIL_DATA,        &
      &                         SMLOW_DATA, SMHIGH_DATA, LVCOEF_DATA, NSLOPE, &
      &                         FRH2O,ZTOPVTBL,ZBOTVTBL, &
      &                         LOW_DENSITY_RESIDENTIAL, HIGH_DENSITY_RESIDENTIAL, HIGH_INTENSITY_INDUSTRIAL, &
      !qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      & LCZ_1,LCZ_2,LCZ_3,LCZ_4,LCZ_5,LCZ_6,LCZ_7,LCZ_8,LCZ_9,LCZ_10,LCZ_11!,err_message !err_message qkyang test-
      !qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  USE module_sf_urban,    only: urban, oasis, IRI_SCHEME
  USE module_sf_noahlsm_glacial_only, only: sflx_glacial
  USE module_sf_bep,      only: bep
  USE module_sf_bep_bem,  only: bep_bem
#if defined(iamas)
use iamas_atmphys_date_time, only: cal_mon_day
use iamas_atmphys_utilities, only: physics_error_fatal!, physics_message ! physics_message qkyang test-
#define FATAL_ERROR(M) call physics_error_fatal( M )
#else
 use module_ra_gfdleta,  only: cal_mon_day
 use module_wrf_error
#define FATAL_ERROR(M) call wrf_error_fatal( M )
#endif
#if ( WRF_CHEM == 1 )
  USE module_data_gocart_dust
#endif

!-------------------------------

!
CONTAINS
!
!----------------------------------------------------------------
! Urban related variable are added to arguments - urban
!----------------------------------------------------------------
   SUBROUTINE lsm(DZ8W,QV3D,P8W3D,T3D,TSK,                      &
                  HFX,QFX,LH,GRDFLX, QGH,GSW,SWDOWN,GLW,SMSTAV,SMSTOT, &
                  SFCRUNOFF, UDRUNOFF,IVGTYP,ISLTYP,ISURBAN,ISICE,iswater,VEGFRA,    &
                  !qkyang added "iswater" to support LCZ landuse dataset  
                  ALBEDO,ALBBCK,ZNT,Z0,TMN,XLAND,XICE,EMISS,EMBCK,   &
                  SNOWC,QSFC,RAINBL,MMINLU,                     &
                  num_soil_layers,DT,DZS,ITIMESTEP,             &
                  SMOIS,TSLB,SNOW,CANWAT,                       &
                  CHS,CHS2,CQS2,CPM,ROVCP,SR,chklowq,lai,qz0,   & !H
                  myj,frpcpn,                                   &
                  SH2O,SNOWH,                                   & !H
                  U_PHY,V_PHY,                                  & !I
                  SNOALB,SHDMIN,SHDMAX,                         & !I
                  SNOTIME,                                      & !?
                  ACSNOM,ACSNOW,                                & !O
                  SNOPCX,                                       & !O
                  POTEVP,                                       & !O
                  SMCREL,                                       & !O
                  XICE_THRESHOLD,                               &
                  RDLAI2D,USEMONALB,                            &
                  RIB,                                          & !?
                  NOAHRES,opt_thcnd,                            &
! Noah UA changes
                  ua_phys,flx4_2d,fvb_2d,fbur_2d,fgsn_2d,       &
                  ids,ide, jds,jde, kds,kde,                    &
                  ims,ime, jms,jme, kms,kme,                    &
                  its,ite, jts,jte, kts,kte,                    &
                  config_urban_scheme,                          & !qkyang replace sf_urban_physics with config_urban_scheme
                  !qkyang added to support anthropogenic heating emissions+++++++++
                  AH_emis,                                      &
                  is_ah_emis,                                   &
                  AH_URB2D,                                     &
                  !qkyang added to support anthropogenic heating emissions----------
                  CMR_SFCDIF,CHR_SFCDIF,CMC_SFCDIF,CHC_SFCDIF,  &
                  CMGR_SFCDIF,CHGR_SFCDIF,                      &
!Optional Urban
                  TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,QC_URB2D, & !H urban
                  UC_URB2D,                                     & !H urban
                  XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D,  & !H urban
                  TRL_URB3D,TBL_URB3D,TGL_URB3D,                & !H urban
                  SH_URB2D,LH_URB2D,G_URB2D,RN_URB2D,TS_URB2D,  & !H urban
                  PSIM_URB2D,PSIH_URB2D,U10_URB2D,V10_URB2D,    & !O urban
                  GZ1OZ0_URB2D,  AKMS_URB2D,                    & !O urban
                  TH2_URB2D,Q2_URB2D, UST_URB2D,                & !O urban
                  DECLIN_URB,COSZ_URB2D,OMG_URB2D,              & !I urban
                  XLAT_URB2D,                                   & !I urban
                  num_roof_layers, num_wall_layers,             & !I urban
                  num_road_layers, DZR, DZB, DZG,               & !I urban
                  CMCR_URB2D,TGR_URB2D,TGRL_URB3D,SMR_URB3D,    & !H urban
                  DRELR_URB2D,DRELB_URB2D,DRELG_URB2D,          & !H urban
                  FLXHUMR_URB2D,FLXHUMB_URB2D,FLXHUMG_URB2D,    & !H urban
                  julian, julyr,                                & !H urban
                  FRC_URB2D,UTYPE_URB2D,                        & !O
                  num_urban_layers,                             & !I multi-layer urban
                  num_urban_hi,                                 & !I multi-layer urban
                  tsk_rural_bep,                                & !H multi-layer urban
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,      & !H multi-layer urban
                  tlev_urb3d,qlev_urb3d,                        & !H multi-layer urban
                  tw1lev_urb3d,tw2lev_urb3d,                    & !H multi-layer urban
                  tglev_urb3d,tflev_urb3d,                      & !H multi-layer urban
                  sf_ac_urb3d,lf_ac_urb3d,cm_ac_urb3d,          & !H multi-layer urban
                  sfvent_urb3d,lfvent_urb3d,                    & !H multi-layer urban
                  sfwin1_urb3d,sfwin2_urb3d,                    & !H multi-layer urban
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,    & !H multi-layer urban
                  lp_urb2d,hi_urb2d,lb_urb2d,hgt_urb2d,         & !H multi-layer urban
                  mh_urb2d,stdh_urb2d,lf_urb2d,                 & !SLUCM
                  th_phy,rho,p_phy,ust,                         & !I multi-layer urban
                  gmt,julday,xlong,xlat,                        & !I multi-layer urban
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,              & !O multi-layer urban
                  a_e_bep,b_u_bep,b_v_bep,                      & !O multi-layer urban
                  b_t_bep,b_q_bep,b_e_bep,dlg_bep,              & !O multi-layer urban
                  dl_u_bep,sf_bep,vl_bep,sfcheadrt,INFXSRT, soldrain      &   !O multi-layer urban         
                 ,SDA_HFX, SDA_QFX, HFX_BOTH, QFX_BOTH, QNORM, fasdas     &   !fasdas
                 ,RC2,XLAI2                                               &
                  )

!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!----------------------------------------------------------------
! --- atmospheric (WRF generic) variables
!-- DT          time step (seconds)
!-- DZ8W        thickness of layers (m)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- FLHC        exchange coefficient for heat (m/s)
!-- FLQC        exchange coefficient for moisture (m/s)
!-- PSFC        surface pressure (Pa)
!-- XLAND       land mask (1 for land, 2 for water)
!-- QGH         saturated mixing ratio at 2 meter
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- History variables
!-- CANWAT      canopy moisture content (mm)
!-- TSK         surface temperature (K)
!-- TSLB        soil temp (k)
!-- SMOIS       total soil moisture content (volumetric fraction)
!-- SH2O        unfrozen soil moisture content (volumetric fraction)
!                note: frozen soil moisture (i.e., soil ice) = SMOIS - SH2O
!-- SNOWH       actual snow depth (m)
!-- SNOW        liquid water-equivalent snow depth (m)
!-- ALBEDO      time-varying surface albedo including snow effect (unitless fraction)
!-- ALBBCK      background surface albedo (unitless fraction)
!-- CHS          surface exchange coefficient for heat and moisture (m s-1);
!-- CHS2        2m surface exchange coefficient for heat  (m s-1);
!-- CQS2        2m surface exchange coefficient for moisture (m s-1);
! --- soil variables
!-- num_soil_layers   the number of soil layers
!-- ZS          depths of centers of soil layers   (m)
!-- DZS         thicknesses of soil layers (m)
!-- SLDPTH      thickness of each soil layer (m, same as DZS)
!-- TMN         soil temperature at lower boundary (K)
!-- SMCWLT      wilting point (volumetric)
!-- SMCDRY      dry soil moisture threshold where direct evap from
!               top soil layer ends (volumetric)
!-- SMCREF      soil moisture threshold below which transpiration begins to
!                   stress (volumetric)
!-- SMCMAX      porosity, i.e. saturated value of soil moisture (volumetric)
!-- NROOT       number of root layers, a function of veg type, determined
!               in subroutine redprm.
!-- SMSTAV      Soil moisture availability for evapotranspiration (
!                   fraction between SMCWLT and SMCMXA)
!-- SMSTOT      Total soil moisture content frozen+unfrozen) in the soil column (mm)
! --- snow variables
!-- SNOWC       fraction snow coverage (0-1.0)
! --- vegetation variables
!-- SNOALB      upper bound on maximum albedo over deep snow
!-- SHDMIN      minimum areal fractional coverage of annual green vegetation
!-- SHDMAX      maximum areal fractional coverage of annual green vegetation
!-- XLAI        leaf area index (dimensionless)
!-- XLAI2       leaf area index (same as XLAI) passed to output (dimensionless)
!-- Z0BRD       Background fixed roughness length (M)
!-- Z0          Background vroughness length (M) as function
!-- ZNT         Time varying roughness length (M) as function
!-- ALBD(IVGTPK,ISN) background albedo reading from a table
! --- LSM output
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          upward moisture flux at the surface (W m-2)
!-- grdflx(I) ground heat flux (W m-2)
!-- FDOWN       radiation forcing at the surface (W m-2) = SOLDN*(1-alb)+LWDN
!----------------------------------------------------------------------------
!-- EC          canopy water evaporation ((W m-2)
!-- EDIR        direct soil evaporation (W m-2)
!-- ET          plant transpiration from a particular root layer (W m-2)
!-- ETT         total plant transpiration (W m-2)
!-- ESNOW       sublimation from (or deposition to if <0) snowpack (W m-2)
!-- DRIP        through-fall of precip and/or dew in excess of canopy
!                 water-holding capacity (m)
!-- DEW         dewfall (or frostfall for t<273.15) (M)
!-- SMAV        Soil Moisture Availability for each layer, as a fraction
!                 between SMCWLT and SMCMAX (dimensionless fraction)
! ----------------------------------------------------------------------
!-- BETA        ratio of actual/potential evap (dimensionless)
!-- ETP         potential evaporation (W m-2)
! ----------------------------------------------------------------------
!-- FLX1        precip-snow sfc (W m-2)
!-- FLX2        freezing rain latent heat flux (W m-2)
!-- FLX3        phase-change heat flux from snowmelt (W m-2)
! ----------------------------------------------------------------------
!-- ACSNOM      snow melt (mm) (water equivalent)
!-- ACSNOW      accumulated snow fall (mm) (water equivalent)
!-- SNOPCX      snow phase change heat flux (W/m^2)
!-- POTEVP      accumulated potential evaporation (m)
!-- RIB         Documentation needed!!!
! ----------------------------------------------------------------------
!-- RUNOFF1     surface runoff (m s-1), not infiltrating the surface
!-- RUNOFF2     subsurface runoff (m s-1), drainage out bottom of last
!                  soil layer (baseflow)
!  important note: here RUNOFF2 is actually the sum of RUNOFF2 and RUNOFF3
!-- RUNOFF3     numerical trunctation in excess of porosity (smcmax)
!                  for a given soil layer at the end of a time step (m s-1).
!SFCRUNOFF     Surface Runoff (mm)
!UDRUNOFF      Total Underground Runoff (mm), which is the sum of RUNOFF2 and RUNOFF3
! ----------------------------------------------------------------------
!-- RC          canopy resistance (s m-1)
!-- RC2         canopy resistance (same as RC) passed to output
!-- PC          plant coefficient (unitless fraction, 0-1) where PC*ETP = actual transp
!-- RSMIN       minimum canopy resistance (s m-1)
!-- RCS         incoming solar rc factor (dimensionless)
!-- RCT         air temperature rc factor (dimensionless)
!-- RCQ         atmos vapor pressure deficit rc factor (dimensionless)
!-- RCSOIL      soil moisture rc factor (dimensionless)

!-- EMISS       surface emissivity (between 0 and 1)
!-- EMBCK       Background surface emissivity (between 0 and 1)

!-- ROVCP       R/CP
!               (R_d/R_v) (dimensionless)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!
!-- SR          fraction of frozen precip (0.0 to 1.0)
!----------------------------------------------------------------

! IN only
   !Changed by qkyang to support urban model+++++++++++++++++++++
   !INTEGER,  INTENT(IN)      :: sf_urban_physics               !original code
   INTEGER, INTENT(IN)      :: config_urban_scheme !replace sf_urban_physics with config_urban_scheme
   INTEGER                     :: sf_urban_physics !added by qkyang
   INTEGER,  INTENT(IN   )   ::  julyr!added by qkyang
   REAL,  INTENT(IN   )   ::  julian!added by qkyang
   !Changed by qkyang to support urban model--------------------- 
   !Changed by qkyang to speed iAMAS model+++++++++++++++++++++
   REAL,    DIMENSION(its:ite), INTENT(IN)::                    & 
                                                          XICE, &
                                                         XLAND, &
                                                           TMN, &
                                                        VEGFRA, &
                                                        SHDMIN, &
                                                        SHDMAX, &
                                                        SNOALB, &
                                                           GLW, &
                                                           GSW, &
                                                        SWDOWN, &
                                                        RAINBL, &
                                                         EMBCK, &
                                                            SR, &
                                                           CPM, &
                                                           QZ0
   INTEGER, DIMENSION(its:ite), INTENT(IN)::                    &
                                                        IVGTYP, &
                                                        ISLTYP

   REAL,    DIMENSION( its:ite ), INTENT(INOUT)::      &
                                                           TSK, &
                                                           QGH
                                                           
   REAL,    DIMENSION(its:ite), INTENT(INOUT)::                 &
                                                           HFX, &
                                                           QFX, &
                                                            LH, &
                                                        GRDFLX, &
                                                          QSFC, &
                                                          CQS2, &
                                                           CHS, &
                                                          CHS2, &
                                                          SNOW, &
                                                         SNOWC, &
                                                         SNOWH, & !new
                                                        CANWAT, &
                                                        SMSTAV, &
                                                        SMSTOT, &
                                                     SFCRUNOFF, &
                                                      UDRUNOFF, &
                                                        ACSNOM, &
                                                        ACSNOW, &
                                                       SNOTIME, &
                                                        SNOPCX, &
                                                         EMISS, &
                                                           RIB, &
                                                        POTEVP, &
                                                        ALBEDO, &
                                                           ZNT, &
                                                        ALBBCK, &
                                                            Z0, &
                                                           LAI
REAL,    DIMENSION(its:ite), INTENT(OUT)::                      &
                                                       NOAHRES, &
                                                       CHKLOWQ

REAL,DIMENSION(1:num_soil_layers,its:ite),INTENT(INOUT) ::      &
                                                         SMOIS, & ! total soil moisture
                                                          SH2O, & ! new soil liquid
                                                          TSLB    ! TSLB     STEMP

REAL,DIMENSION(1:num_soil_layers,its:ite),INTENT(OUT)::         &
                                                         SMCREL

REAL,OPTIONAL, DIMENSION(its:ite), INTENT(IN)::                 &
                                                          XLAT, &
                                                         XLONG, &
                                                    XLAT_URB2D

!===Yang,2014/10/08,hydrological variable for single layer UCM===                   
REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) ::          &
                                                   DRELR_URB2D, &
                                                   DRELB_URB2D, &
                                                   DRELG_URB2D, &
                                                 FLXHUMR_URB2D, &
                                                 FLXHUMB_URB2D, &
                                                 FLXHUMG_URB2D, &
                                                    CMCR_URB2D, &
                                                     TGR_URB2D
REAL,OPTIONAL,DIMENSION(1:num_roof_layers,its:ite),INTENT(INOUT)::   &
                                                         TGRL_URB3D, &
                                                          SMR_URB3D
! state variable surface_driver <--> lsm <--> urban
 REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) ::        &
                                                     TR_URB2D, &
                                                     TB_URB2D, &
                                                     TG_URB2D, &
                                                     TC_URB2D, &
                                                     QC_URB2D, &
                                                     UC_URB2D, &
                                                   XXXR_URB2D, &
                                                   XXXB_URB2D, &
                                                   XXXG_URB2D, &
                                                   XXXC_URB2D, &
                                                     SH_URB2D, &
                                                     LH_URB2D, &
                                                      G_URB2D, &
                                                     RN_URB2D, &
                                                     TS_URB2D
REAL, OPTIONAL, DIMENSION( 1:num_roof_layers, its:ite ), INTENT(INOUT) :: TRL_URB3D
REAL, OPTIONAL, DIMENSION( 1:num_wall_layers, its:ite ), INTENT(INOUT) :: TBL_URB3D
REAL, OPTIONAL, DIMENSION( 1:num_road_layers, its:ite ), INTENT(INOUT) :: TGL_URB3D
! output variable lsm --> surface_driver
 REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) ::        &
                                                 PSIM_URB2D, &
                                                 PSIH_URB2D, &
                                               GZ1OZ0_URB2D, &
                                                  U10_URB2D, &
                                                  V10_URB2D, &
                                                  TH2_URB2D, &
                                                   Q2_URB2D, &
                                                 AKMS_URB2D
 REAL, DIMENSION( its:ite ), INTENT(OUT) :: UST_URB2D
 REAL, DIMENSION( its:ite ), INTENT(IN) :: FRC_URB2D
 INTEGER, DIMENSION( its:ite ), INTENT(IN) :: UTYPE_URB2D

 ! Variables for multi-layer UCM (Martilli et al. 2002)
   REAL, OPTIONAL, DIMENSION(                     its:ite ), INTENT(INOUT) :: tsk_rural_bep
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: trb_urb4d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tw1_urb4d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tw2_urb4d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tgb_urb4d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tlev_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: qlev_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tw1lev_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tw2lev_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tglev_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: tflev_urb3d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: lf_ac_urb3d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: sf_ac_urb3d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: cm_ac_urb3d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: sfvent_urb3d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: lfvent_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfwin1_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfwin2_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfw1_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfw2_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfr_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_layers, its:ite ), INTENT(INOUT) :: sfg_urb3d
   REAL, OPTIONAL, DIMENSION( 1:num_urban_hi, its:ite ), INTENT(IN) :: hi_urb2d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(IN) :: lp_urb2d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(IN) :: lb_urb2d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(IN) :: hgt_urb2d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(IN) :: mh_urb2d
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(IN) :: stdh_urb2d
   REAL, OPTIONAL, DIMENSION( 4, its:ite ), INTENT(IN) :: lf_urb2d
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::a_u_bep   !Implicit momemtum component X-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::a_v_bep   !Implicit momemtum component Y-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::a_t_bep   !Implicit component pot. temperature
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::a_q_bep   !Implicit momemtum component X-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::a_e_bep   !Implicit component TKE
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::b_u_bep   !Explicit momentum component X-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::b_v_bep   !Explicit momentum component Y-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::b_t_bep   !Explicit component pot. temperature
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::b_q_bep   !Implicit momemtum component Y-direction
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::b_e_bep   !Explicit component TKE
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::vl_bep    !Fraction air volume in grid cell
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::dlg_bep   !Height above ground
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::sf_bep    !Fraction air at the face of grid cell
   REAL, OPTIONAL, DIMENSION( kms:kme, its:ite ), INTENT(INOUT) ::dl_u_bep  !Length scale

! Local variables for multi-layer UCM (Martilli et al. 2002)
   REAL,    DIMENSION( its:ite ) :: HFX_RURAL,LH_RURAL,GRDFLX_RURAL ! ,RN_RURAL
   REAL,    DIMENSION( its:ite ) :: QFX_RURAL ! ,QSFC_RURAL,UMOM_RURAL,VMOM_RURAL
   REAL,    DIMENSION( its:ite ) :: ALB_RURAL,EMISS_RURAL,TSK_RURAL ! ,UST_RURAL
!   REAL,    DIMENSION( its:ite ) :: QSFC_URB
   REAL,    DIMENSION( its:ite ) :: HFX_URB,UMOM_URB,VMOM_URB
   REAL,    DIMENSION( its:ite ) :: QFX_URB
!    REAL, OPTIONAL, DIMENSION( its:ite ) :: ALBEDO_URB,EMISS_URB,UMOM,VMOM,UST 
   REAL, DIMENSION(its:ite) ::EMISS_URB
   REAL, DIMENSION(its:ite) ::RL_UP_URB
   REAL, DIMENSION(its:ite) ::RS_ABS_URB
   REAL, DIMENSION(its:ite) ::GRDFLX_URB
   REAL, OPTIONAL, DIMENSION(its:ite), INTENT(IN) :: COSZ_URB2D
   REAL, OPTIONAL, DIMENSION(its:ite), INTENT(IN) :: OMG_URB2D
   REAL, DIMENSION(its:ite), INTENT(OUT) :: FLX4_2D,FVB_2D,FBUR_2D,FGSN_2D
   !Changed by qkyang to speed iAMAS model--------------------- 

!qkyang added to support anthropogenic heating emissions+++++++
   REAL,    DIMENSION(its:ite), INTENT(IN):: AH_emis
   logical,pointer:: is_ah_emis
   REAL,    DIMENSION(its:ite), INTENT(OUT):: AH_URB2D
!qkyang added to support anthropogenic heating emissions-------

   INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte

!    INTEGER,  INTENT(IN   )   ::  sf_urban_physics               !urban !qkyang comment
   INTEGER,  INTENT(IN   )   ::  isurban
   INTEGER,  INTENT(IN   )   ::  isice,iswater !qkyang added "iswater" to support LCZ landuse dataset
   !INTEGER,  INTENT(IN   )   ::  julian, julyr                  !urban !qkyang comment

!Changed by qkyang to speed iAMAS model+++++++++++++++++++++ 
!added by Wei Yu  for routing
    REAL,    DIMENSION( its:ite )                     , &
             INTENT(INOUT)  :: sfcheadrt,INFXSRT,soldrain 
!Changed by qkyang to speed iAMAS model--------------------- 
    real :: etpnd1
!end added

   CHARACTER(LEN=*), INTENT(IN   )    ::                 MMINLU

   REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )            , &
            INTENT(IN   )    ::                           QV3D, &
                                                         p8w3D, &
                                                          DZ8W, &
                                                          T3D                

   INTEGER, INTENT(IN)       ::     num_soil_layers,ITIMESTEP

   REAL,     INTENT(IN   )   ::     DT,ROVCP

   REAL,     DIMENSION(1:num_soil_layers), INTENT(IN)::DZS

! IN and OUT
      

   INTEGER, INTENT(IN)       ::                         OPT_THCND

! Noah UA changes
   LOGICAL,                                INTENT(IN)  :: UA_PHYS
   REAL                                                :: FLX4,FVB,FBUR,FGSN
   !Changed by qkyang to speed iAMAS model+++++++++++++++++++++ 
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) ::  RC2, XLAI2
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMR_SFCDIF
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHR_SFCDIF
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMGR_SFCDIF
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHGR_SFCDIF
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMC_SFCDIF
   REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHC_SFCDIF
   !Changed by qkyang to speed iAMAS model--------------------- 
! Local variables (moved here from driver to make routine thread safe, 20031007 jm)

      REAL, DIMENSION(1:num_soil_layers) ::  ET

      REAL, DIMENSION(1:num_soil_layers) ::  SMAV

      REAL  ::  BETA, ETP, SSOIL,EC, EDIR, ESNOW, ETT,        &
                FLX1,FLX2,FLX3, DRIP,DEW,FDOWN,RC,PC,RSMIN,XLAI,  &
!                RCS,RCT,RCQ,RCSOIL
                RCS,RCT,RCQ,RCSOIL,FFROZP

    LOGICAL,    INTENT(IN   )    ::     myj,frpcpn

! DECLARATIONS - LOGICAL
! ----------------------------------------------------------------------
      LOGICAL, PARAMETER :: LOCAL=.false.
      LOGICAL :: FRZGRA, SNOWNG

      LOGICAL :: IPRINT

! ----------------------------------------------------------------------
! DECLARATIONS - INTEGER
! ----------------------------------------------------------------------
      INTEGER :: I,J, ICE,NSOIL,SLOPETYP,SOILTYP,VEGTYP
      INTEGER :: NROOT
      INTEGER :: KZ ,K
      INTEGER :: NS
! ----------------------------------------------------------------------
! DECLARATIONS - REAL
! ----------------------------------------------------------------------

      REAL  :: SHMIN,SHMAX,DQSDT2,LWDN,PRCP,PRCPRAIN,                    &
               Q2SAT,Q2SATI,SFCPRS,SFCSPD,SFCTMP,SHDFAC,SNOALB1,         &
               SOLDN,TBOT,ZLVL, Q2K,ALBBRD, ALBEDOK, ETA, ETA_KINEMATIC, &
               EMBRD,                                                    &
               Z0K,RUNOFF1,RUNOFF2,RUNOFF3,SHEAT,SOLNET,E2SAT,SFCTSNO,   &
! mek, WRF testing, expanded diagnostics
               SOLUP,LWUP,RNET,RES,Q1SFC,TAIRV,SATFLG
! MEK MAY 2007
      REAL ::  FDTLIW
! MEK JUL2007 for pot. evap.
      REAL :: RIBB
      REAL ::  FDTW

      REAL  :: EMISSI

      REAL  :: SNCOVR,SNEQV,SNOWHK,CMC, CHK,TH2

      REAL  :: SMCDRY,SMCMAX,SMCREF,SMCWLT,SNOMLT,SOILM,SOILW,Q1,T1
      REAL  :: SNOTIME1    ! LSTSNW1 INITIAL NUMBER OF TIMESTEPS SINCE LAST SNOWFALL

      REAL  :: DUMMY,Z0BRD
!
      REAL  :: COSZ, SOLARDIRECT
!
      REAL, DIMENSION(1:num_soil_layers)::  SLDPTH, STC,SMC,SWC
!
      REAL, DIMENSION(1:num_soil_layers) ::     ZSOIL, RTDIS
      REAL, PARAMETER  :: TRESH=.95E0, A2=17.67,A3=273.15,A4=29.65,   &
                          T0=273.16E0, ELWV=2.50E6,  A23M4=A2*(A3-A4)
! MEK MAY 2007
      REAL, PARAMETER  :: ROW=1.E3,ELIW=XLF,ROWLIW=ROW*ELIW

! ----------------------------------------------------------------------
! DECLARATIONS START - urban
! ----------------------------------------------------------------------

! input variables surface_driver --> lsm
     INTEGER, INTENT(IN) :: num_roof_layers
     INTEGER, INTENT(IN) :: num_wall_layers
     INTEGER, INTENT(IN) :: num_road_layers
     REAL, OPTIONAL, DIMENSION(1:num_roof_layers), INTENT(IN) :: DZR
     REAL, OPTIONAL, DIMENSION(1:num_wall_layers), INTENT(IN) :: DZB
     REAL, OPTIONAL, DIMENSION(1:num_road_layers), INTENT(IN) :: DZG
     REAL, OPTIONAL, INTENT(IN) :: DECLIN_URB
     !-- u_phy         u-velocity interpolated to theta points (m/s)
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: U_PHY
     !-- v_phy         v-velocity interpolated to theta points (m/s)
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: V_PHY
     !potential temperature  [K]
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: TH_PHY

     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: P_PHY
     REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: RHO
!Changed by qkyang to speed iAMAS model+++++++++++++++++++++ 
     REAL, OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: UST
!Changed by qkyang to speed iAMAS model--------------------- 

     LOGICAL, intent(in) :: rdlai2d
     LOGICAL, intent(in) :: USEMONALB

! input variables lsm --> urban
     INTEGER :: UTYPE_URB ! urban type [urban=1, suburban=2, rural=3]
     REAL :: TA_URB       ! potential temp at 1st atmospheric level [K]
     REAL :: QA_URB       ! mixing ratio at 1st atmospheric level  [kg/kg]
     REAL :: UA_URB       ! wind speed at 1st atmospheric level    [m/s]
     REAL :: U1_URB       ! u at 1st atmospheric level             [m/s]
     REAL :: V1_URB       ! v at 1st atmospheric level             [m/s]
     REAL :: SSG_URB      ! downward total short wave radiation    [W/m/m]
     REAL :: LLG_URB      ! downward long wave radiation           [W/m/m]
     REAL :: RAIN_URB     ! precipitation                          [mm/h]
     REAL :: RHOO_URB     ! air density                            [kg/m^3]
     REAL :: ZA_URB       ! first atmospheric level                [m]
     REAL :: DELT_URB     ! time step                              [s]
     REAL :: SSGD_URB     ! downward direct short wave radiation   [W/m/m]
     REAL :: SSGQ_URB     ! downward diffuse short wave radiation  [W/m/m]
     REAL :: XLAT_URB     ! latitude                               [deg]
     REAL :: COSZ_URB     ! cosz
     REAL :: OMG_URB      ! hour angle
     REAL :: ZNT_URB      ! roughness length                       [m]
     REAL :: TR_URB
     REAL :: TB_URB
     REAL :: TG_URB
     REAL :: TC_URB
     REAL :: QC_URB
     REAL :: UC_URB
     REAL :: XXXR_URB
     REAL :: XXXB_URB
     REAL :: XXXG_URB
     REAL :: XXXC_URB
     REAL, DIMENSION(1:num_roof_layers) :: TRL_URB  ! roof layer temp [K]
     REAL, DIMENSION(1:num_wall_layers) :: TBL_URB  ! wall layer temp [K]
     REAL, DIMENSION(1:num_road_layers) :: TGL_URB  ! road layer temp [K]
     LOGICAL  :: LSOLAR_URB

!===Yang,2014/10/08,hydrological variable for single layer UCM===
     INTEGER :: jmonth, jday, tloc
     INTEGER :: IRIOPTION, USOIL, DSOIL
     REAL :: AOASIS, OMG
     REAL :: DRELR_URB
     REAL :: DRELB_URB
     REAL :: DRELG_URB
     REAL :: FLXHUMR_URB
     REAL :: FLXHUMB_URB
     REAL :: FLXHUMG_URB
     REAL :: CMCR_URB
     REAL :: TGR_URB
     REAL, DIMENSION(1:num_roof_layers) :: SMR_URB  ! green roof layer moisture
     REAL, DIMENSION(1:num_roof_layers) :: TGRL_URB ! green roof layer temp [K]


! output variables urban --> lsm
     REAL :: TS_URB     ! surface radiative temperature    [K]
     REAL :: QS_URB     ! surface humidity                 [-]
     REAL :: SH_URB     ! sensible heat flux               [W/m/m]
     REAL :: LH_URB     ! latent heat flux                 [W/m/m]
     REAL :: LH_KINEMATIC_URB ! latent heat flux, kinetic  [kg/m/m/s]
     REAL :: SW_URB     ! upward short wave radiation flux [W/m/m]
     REAL :: ALB_URB    ! time-varying albedo            [fraction]
     REAL :: LW_URB     ! upward long wave radiation flux  [W/m/m]
     REAL :: G_URB      ! heat flux into the ground        [W/m/m]
     REAL :: RN_URB     ! net radiation                    [W/m/m]
     REAL :: PSIM_URB   ! shear f for momentum             [-]
     REAL :: PSIH_URB   ! shear f for heat                 [-]
     REAL :: GZ1OZ0_URB   ! shear f for heat                 [-]
     REAL :: U10_URB    ! wind u component at 10 m         [m/s]
     REAL :: V10_URB    ! wind v component at 10 m         [m/s]
     REAL :: TH2_URB    ! potential temperature at 2 m     [K]
     REAL :: Q2_URB     ! humidity at 2 m                  [-]
     REAL :: CHS_URB
     REAL :: CHS2_URB
     REAL :: UST_URB
! NUDAPT Parameters urban --> lam
     REAL :: mh_urb
     REAL :: stdh_urb
     REAL :: lp_urb
     REAL :: hgt_urb
     REAL, DIMENSION(4) :: lf_urb
! Variables for multi-layer UCM (Martilli et al. 2002)
   REAL, OPTIONAL, INTENT(IN  )   ::                                   GMT 
   INTEGER, OPTIONAL, INTENT(IN  ) ::                               JULDAY
   INTEGER, INTENT(IN  ) ::                               NUM_URBAN_LAYERS
   INTEGER, INTENT(IN  ) ::                               NUM_URBAN_HI
   REAL :: SIGMA_SB,RL_UP_RURAL,RL_UP_TOT,RS_ABS_TOT,UMOM,VMOM
   REAL :: CMR_URB, CHR_URB, CMC_URB, CHC_URB, CMGR_URB, CHGR_URB
   REAL :: frc_urb,lb_urb
   REAL :: check 
! ----------------------------------------------------------------------
! DECLARATIONS END - urban
! ----------------------------------------------------------------------

  REAL, PARAMETER  :: CAPA=R_D/CP
  REAL :: APELM,APES,SFCTH2,PSFC
  real, intent(in) :: xice_threshold
  character(len=80) :: message_text
!
!  FASDAS
!
   REAL,    DIMENSION( ims:ime, jms:jme ),  OPTIONAL,                     &
            INTENT(INOUT)  ::   SDA_HFX, SDA_QFX, HFX_BOTH, QFX_BOTH, QNORM
   INTEGER, INTENT(IN   )  ::  fasdas
!  local vars
   REAL                    ::  XSDA_HFX, XSDA_QFX, XQNORM
   REAL                    ::  HFX_PHY, QFX_PHY
   REAL                    ::  DZQ
   REAL                    ::  HCPCT_FASDAS

   !Added by qkyang to support urban model+++++++++++++++++++++
   sf_urban_physics =config_urban_scheme
   CALL cal_mon_day(julday,julyr,jmonth,jday) !qkyang replace julian with julday
   !Added by qkyang to support urban model--------------------- 

   HFX_PHY = 0.0   ! initialize
   QFX_PHY = 0.0
   XQNORM  = 0.0
   XSDA_HFX = 0.0
   XSDA_QFX = 0.0
!
!  END FASDAS
!
  FLX4  = 0.0 !BSINGH - Initialized to 0.0
  FVB   = 0.0 !BSINGH - Initialized to 0.0
  FBUR  = 0.0 !BSINGH - Initialized to 0.0
  FGSN  = 0.0 !BSINGH - Initialized to 0.0
  SOILW = 0.0 !BSINGH - Initialized to 0.0

  sigma_sb=5.67e-08

! MEK MAY 2007
      FDTLIW=DT/ROWLIW
! MEK JUL2007
      FDTW=DT/(XLV*RHOWATER)
! debug printout
         IPRINT=.false.

!      SLOPETYP=2
      SLOPETYP=1
!      SHDMIN=0.00

      NSOIL=num_soil_layers

     DO NS=1,NSOIL
     SLDPTH(NS)=DZS(NS)
     ENDDO

     JLOOP : DO J=jts,jte

      IF(ITIMESTEP.EQ.1)THEN
        DO 50 I=its,ite
!*** initialize soil conditions for IHOP 31 May case
!         IF((XLAND(I)-1.5) < 0.)THEN
!            if (I==108.and.j==85) then
!                  DO NS=1,NSOIL
!                      smois(NS,I)=0.10
!                      SH2O(NS,I)=0.10
!                  enddo
!             endif
!         ENDIF

!*** SET ZERO-VALUE FOR SOME OUTPUT DIAGNOSTIC ARRAYS
!           IF((XLAND(I)-1.5).GE.0.)THEN
!qkyang  added ".or. VEGTYP.eq.iswater" to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF(((XLAND(I)-1.5).GE.0.) .or. (ivgtyp(I) == iswater))THEN                                  ! begining of land/sea if block
!qkyang  added ".or. VEGTYP.eq.iswater" to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! check sea-ice point
#if 0
            IF( XICE(I).GE. XICE_THRESHOLD .and. IPRINT ) PRINT*, ' sea-ice at water point, I=',I,'J=',J
#endif
!***   Open Water Case
            smstav(I)=1.0
            smstot(I)=1.0
            DO NS=1,NSOIL
              smois(NS,I)=1.0
              tslb(NS,I)=273.16                                          !STEMP
              smcrel(NS,I)=1.0
            ENDDO
          ELSE
            IF ( XICE(I) .GE. XICE_THRESHOLD ) THEN
!***        SEA-ICE CASE
              smstav(I)=1.0
              smstot(I)=1.0
              DO NS=1,NSOIL
                smois(NS,I)=1.0
                smcrel(NS,I)=1.0
              ENDDO
            ENDIF
          ENDIF
!
   50   CONTINUE
      ENDIF                                                               ! end of initialization over ocean

!-----------------------------------------------------------------------
      ILOOP : DO I=its,ite
! surface pressure
        PSFC=P8w3D(i,1,j)
! pressure in middle of lowest layer
        SFCPRS=(P8W3D(I,KTS+1,j)+P8W3D(i,KTS,j))*0.5
! convert from mixing ratio to specific humidity
         Q2K=QV3D(i,1,j)/(1.0+QV3D(i,1,j))
!
!         Q2SAT=qgh(I,j)
         Q2SAT=qgh(I)/(1.0+qgh(I))        ! Q2SAT is sp humidity
! add check on myj=.true.
!        IF((Q2K.GE.Q2SAT*TRESH).AND.Q2K.LT.qz0(I))THEN
        IF((myj).AND.(Q2K.GE.Q2SAT*TRESH).AND.Q2K.LT.qz0(I))THEN
          SATFLG=0.
          chklowq(I)=0.
        ELSE
          SATFLG=1.0
          chklowq(I)=1.
        ENDIF

        SFCTMP=T3D(i,1,j)
        ZLVL=0.5*DZ8W(i,1,j)

!        TH2=SFCTMP+(0.0097545*ZLVL)
! calculate SFCTH2 via Exner function vs lapse-rate (above)
         APES=(1.E5/PSFC)**CAPA
         APELM=(1.E5/SFCPRS)**CAPA
         SFCTH2=SFCTMP*APELM
         TH2=SFCTH2/APES
!
         EMISSI = emiss(I)
         LWDN=glw(I)*EMISSI
! SOLDN is total incoming solar
        SOLDN=swdown(I)
! GSW is net downward solar
!        SOLNET=gsw(I)
! use mid-day albedo to determine net downward solar (no solar zenith angle correction)
        SOLNET=SOLDN*(1.-albedo(I))
        PRCP=rainbl(I)/DT
        VEGTYP=ivgtyp(I)
        SOILTYP=ISLTYP(I)
        SHDFAC=vegfra(I)/100.
        T1=tsk(I)
        CHK=chs(I)
        SHMIN=SHDMIN(I)/100. !NEW
        SHMAX=shdmax(I)/100. !NEW
! convert snow water equivalent from mm to meter
        SNEQV=snow(I)*0.001
! snow depth in meters
        SNOWHK=snowh(I)
        SNCOVR=snowc(I)

! if "SR" present, set frac of frozen precip ("FFROZP") = snow-ratio ("SR", range:0-1)
! SR from e.g. Ferrier microphysics
! otherwise define from 1st atmos level temperature
       IF(FRPCPN) THEN
          FFROZP=sr(I)
        ELSE
          IF (SFCTMP <=  273.15) THEN
            FFROZP = 1.0
	  ELSE
	    FFROZP = 0.0
	  ENDIF
        ENDIF
!qkyang  added ".or. VEGTYP.eq.iswater" to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF(((XLAND(I)-1.5).GE.0.) .or. (ivgtyp(I) == iswater))THEN                                  ! begining of land/sea if block
!qkyang  added ".or. VEGTYP.eq.iswater" to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Open water points
          TSK_RURAL(I)=tsk(I)
          HFX_RURAL(i)=hfx(I)
          QFX_RURAL(i)=qfx(I)
          LH_RURAL(i)=lh(I)
          emiss_rural(i)=emiss(I)
          GRDFLX_RURAL(i)=grdflx(I)

        ELSE
! Land or sea-ice case

          IF (XICE(I) >= XICE_THRESHOLD) THEN
             ! Sea-ice point
             ICE = 1
          ELSE IF ( VEGTYP == ISICE ) THEN
             ! Land-ice point
             ICE = -1
          ELSE
             ! Neither sea ice or land ice.
             ICE=0
          ENDIF
          DQSDT2=Q2SAT*A23M4/(SFCTMP-A4)**2

          IF(snow(I).GT.0.0)THEN
! snow on surface (use ice saturation properties)
            SFCTSNO=SFCTMP
            E2SAT=611.2*EXP(6174.*(1./273.15 - 1./SFCTSNO))
            Q2SATI=0.622*E2SAT/(SFCPRS-E2SAT)
            Q2SATI=Q2SATI/(1.0+Q2SATI)    ! spec. hum.
            IF (T1 .GT. 273.14) THEN
! warm ground temps, weight the saturation between ice and water according to SNOWC
              Q2SAT=Q2SAT*(1.-snowc(I)) + Q2SATI*snowc(I)
              DQSDT2=DQSDT2*(1.-snowc(I)) + Q2SATI*6174./(SFCTSNO**2)*snowc(I)
            ELSE
! cold ground temps, use ice saturation only
              Q2SAT=Q2SATI
              DQSDT2=Q2SATI*6174./(SFCTSNO**2)
            ENDIF
! for snow cover fraction at 0 C, ground temp will not change, so DQSDT2 effectively zero
! V3.8 add condition for SWDOWN to restrict condition to positive forcing (JD)
            IF(T1 .GT. 273. .AND. snowc(I) .GT. 0. .AND. swdown(I) .GT. 10.)DQSDT2=DQSDT2*(1.-snowc(I))
          ENDIF

          ! Land-ice or land points use the usual deep-soil temperature.
          TBOT=tmn(I)

          IF(ISURBAN.EQ.1) THEN
! assumes these only need to be set for USGS land data
          IF(VEGTYP.EQ.25) SHDFAC=0.0000
          IF(VEGTYP.EQ.26) SHDFAC=0.0000
          IF(VEGTYP.EQ.27) SHDFAC=0.0000
          ENDIF
          IF(SOILTYP.EQ.14.AND.XICE(I).EQ.0.)THEN
#if 0
         IF(IPRINT)PRINT*,' SOIL TYPE FOUND TO BE WATER AT A LAND-POINT'
         IF(IPRINT)PRINT*,i,j,'RESET SOIL in surfce.F'
#endif
            SOILTYP=7
          ENDIF
          SNOALB1 = snoalb(I)
! converts canwat in mm to CMC in meters
          CMC=canwat(I)/1000.

!-------------------------------------------
!*** convert snow depth from mm to meter
!
!          IF(RDMAXALB) THEN
!           SNOALB=ALBMAX(I,J)*0.01
!         ELSE
!           SNOALB=MAXALB(IVGTPK)*0.01
!         ENDIF

!        SNOALB1=0.80
!        SHMIN=0.00
        ALBBRD=albbck(I)
        Z0BRD=z0(I)
        EMBRD=embck(I)
        SNOTIME1 = snotime(I)
        RIBB=rib(I)
!FEI: temporaray arrays above need to be changed later by using SI

          DO NS=1,NSOIL
            SMC(NS)=smois(NS,I)
            STC(NS)=tslb(NS,I)                                          !STEMP
            SWC(NS)=SH2O(NS,I)
          ENDDO
!
          if ( (SNEQV.ne.0..AND.SNOWHK.eq.0.).or.(SNOWHK.le.SNEQV) )THEN
            SNOWHK= 5.*SNEQV
          endif
!

!Fei: urban. for urban surface, if calling UCM, redefine the natural surface in cities as
! the "NATURAL" category in the VEGPARM.TBL

	   IF(SF_URBAN_PHYSICS == 1.OR. SF_URBAN_PHYSICS==2.OR.SF_URBAN_PHYSICS==3 ) THEN
         !qkyang replaced to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!            IF( ivgtyp(I) == ISURBAN .or. ivgtyp(I) == LOW_DENSITY_RESIDENTIAL .or. &
!             ivgtyp(I) == HIGH_DENSITY_RESIDENTIAL .or. ivgtyp(I) == HIGH_INTENSITY_INDUSTRIAL) THEN
         IF(IVGTYP(I) == ISURBAN .or. &
            IVGTYP(I) == low_density_residential .or. &  
            IVGTYP(I) == HIGH_DENSITY_RESIDENTIAL .or. &
            IVGTYP(I) == HIGH_INTENSITY_INDUSTRIAL .or. &
            IVGTYP(I) == LCZ_1 .or. IVGTYP(I) == LCZ_2 .or. &
            IVGTYP(I) == LCZ_3 .or. IVGTYP(I) == LCZ_4 .or. &
            IVGTYP(I) == LCZ_5 .or. IVGTYP(I) == LCZ_6 .or. &
            IVGTYP(I) == LCZ_7 .or. IVGTYP(I) == LCZ_8 .or. &
            IVGTYP(I) == LCZ_9 .or. IVGTYP(I) == LCZ_10 .or. &
            IVGTYP(I) == LCZ_11 )  then !urban
         !qkyang replaced to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		     VEGTYP = NATURAL
           SHDFAC = SHDTBL(NATURAL)
           ALBEDOK =0.2         !  0.2
           ALBBRD  =0.2         !0.2
           EMISSI = 0.98                                 !for VEGTYP=5
		     IF ( frc_urb2d(i) < 0.99 ) THEN
              if(sf_urban_physics.eq.1)then
                T1= ( tsk(I) -frc_urb2d(i) * TS_URB2D(i) )/ (1-frc_urb2d(i))
              elseif((sf_urban_physics.eq.2).OR.(sf_urban_physics.eq.3))then
                T1=TSK_RURAL_BEP(I)
              endif
	        ELSE
		        T1 = tsk(I)
           ENDIF
         ENDIF

     ELSE
           !qkyang replaced to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!           IF( ivgtyp(I) == ISURBAN .or. ivgtyp(I) == LOW_DENSITY_RESIDENTIAL .or. &
!             ivgtyp(I) == HIGH_DENSITY_RESIDENTIAL .or. ivgtyp(I) == HIGH_INTENSITY_INDUSTRIAL) THEN
         IF(IVGTYP(I) == ISURBAN .or. &
            IVGTYP(I) == low_density_residential .or. &  
            IVGTYP(I) == HIGH_DENSITY_RESIDENTIAL .or. &
            IVGTYP(I) == HIGH_INTENSITY_INDUSTRIAL .or. &
            IVGTYP(I) == LCZ_1 .or. IVGTYP(I) == LCZ_2 .or. &
            IVGTYP(I) == LCZ_3 .or. IVGTYP(I) == LCZ_4 .or. &
            IVGTYP(I) == LCZ_5 .or. IVGTYP(I) == LCZ_6 .or. &
            IVGTYP(I) == LCZ_7 .or. IVGTYP(I) == LCZ_8 .or. &
            IVGTYP(I) == LCZ_9 .or. IVGTYP(I) == LCZ_10 .or. &
            IVGTYP(I) == LCZ_11 )  then !urban
         !qkyang replaced to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            VEGTYP = ISURBAN
     	 ENDIF
     ENDIF

!===Yang, 2014/10/08, hydrological processes for urban vegetation in single layer UCM===
           AOASIS = 1.0
           USOIL = 1
           DSOIL = 2
           IRIOPTION=IRI_SCHEME
	  IF(SF_URBAN_PHYSICS == 1) THEN
           OMG= OMG_URB2D(i) !OMG is "SOLAR HOUR ANGLE" ["radians"]!qkyang coment
           tloc=mod(int(OMG/3.14159*180./15.+12.+0.5 ),24)
           if (tloc.lt.0) tloc=tloc+24
           if (tloc==0) tloc=24
           !CALL cal_mon_day(julday,julyr,jmonth,jday) !qkyang replace julian with julday !qkyang coment
           !qkyang replaced to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!                 IF( ivgtyp(I) == ISURBAN .or. ivgtyp(I) == LOW_DENSITY_RESIDENTIAL .or. &
!               ivgtyp(I) == HIGH_DENSITY_RESIDENTIAL .or. ivgtyp(I) == HIGH_INTENSITY_INDUSTRIAL) THEN
           IF(IVGTYP(I) == ISURBAN .or. &
              IVGTYP(I) == low_density_residential .or. &  
              IVGTYP(I) == HIGH_DENSITY_RESIDENTIAL .or. &
              IVGTYP(I) == HIGH_INTENSITY_INDUSTRIAL .or. &
              IVGTYP(I) == LCZ_1 .or. IVGTYP(I) == LCZ_2 .or. &
              IVGTYP(I) == LCZ_3 .or. IVGTYP(I) == LCZ_4 .or. &
              IVGTYP(I) == LCZ_5 .or. IVGTYP(I) == LCZ_6 .or. &
              IVGTYP(I) == LCZ_7 .or. IVGTYP(I) == LCZ_8 .or. &
              IVGTYP(I) == LCZ_9 .or. IVGTYP(I) == LCZ_10 .or. &
              IVGTYP(I) == LCZ_11 )  then !urban
           !qkyang replaced to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            AOASIS = oasis  ! urban oasis effect
                   IF (IRIOPTION ==1) THEN
                   IF (tloc==21 .or. tloc==22) THEN  !irrigation on vegetaion in urban area, MAY-SEP, 9-10pm
                   IF (jmonth==5 .or. jmonth==6 .or. jmonth==7 .or. jmonth==8 .or. jmonth==9) THEN
!                      IF (SMC(USOIL) .LT. SMCREF) SMC(USOIL)= SMCREF
!                      IF (SMC(DSOIL) .LT. SMCREF) SMC(DSOIL)= SMCREF
                      IF (SMC(USOIL) .LT. SMCREF) SMC(USOIL)= REFSMC(ISLTYP(I))                      
                      IF (SMC(DSOIL) .LT. SMCREF) SMC(DSOIL)= REFSMC(ISLTYP(I))                    
                   ENDIF
                   ENDIF
                   ENDIF
                ENDIF
           ENDIF
	  IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN
          IF(AOASIS > 1.0) THEN
          FATAL_ERROR('Urban oasis option is for SF_URBAN_PHYSICS == 1 only')
          ENDIF
          IF(IRIOPTION == 1) THEN
          FATAL_ERROR('Urban irrigation option is for SF_URBAN_PHYSICS == 1 only')
          ENDIF
          ENDIF

#if 0
          IF(IPRINT) THEN
!
       print*, 'BEFORE SFLX, in Noahlsm_driver'
       print*, 'ICE', ICE, 'DT',DT, 'ZLVL',ZLVL, 'NSOIL', NSOIL,   &
       'SLDPTH', SLDPTH, 'LOCAL',LOCAL, 'LUTYPE',&
        LUTYPE, 'SLTYPE',SLTYPE, 'LWDN',LWDN, 'SOLDN',SOLDN,      &
        'SFCPRS',SFCPRS, 'PRCP',PRCP,'SFCTMP',SFCTMP,'Q2K',Q2K,   &
         'TH2',TH2,'Q2SAT',Q2SAT,'DQSDT2',DQSDT2,'VEGTYP', VEGTYP,&
         'SOILTYP',SOILTYP, 'SLOPETYP',SLOPETYP, 'SHDFAC',SHDFAC,&
         'SHMIN',SHMIN, 'ALBBRD',ALBBRD,'SNOALB1',SNOALB1,'TBOT',&
          TBOT, 'Z0BRD',Z0BRD, 'Z0K',Z0K, 'CMC',CMC, 'T1',T1,'STC',&
          STC, 'SMC',SMC, 'SWC',SWC,'SNOWHK',SNOWHK,'SNEQV',SNEQV,&
          'ALBEDOK',ALBEDOK,'CHK',CHK,'ETA',ETA,'SHEAT',SHEAT,      &
          'ETA_KINEMATIC',ETA_KINEMATIC, 'FDOWN',FDOWN,'EC',EC,   &
          'EDIR',EDIR,'ET',ET,'ETT',ETT,'ESNOW',ESNOW,'DRIP',DRIP,&
          'DEW',DEW,'BETA',BETA,'ETP',ETP,'SSOIL',SSOIL,'FLX1',FLX1,&
          'FLX2',FLX2,'FLX3',FLX3,'SNOMLT',SNOMLT,'SNCOVR',SNCOVR,&
          'RUNOFF1',RUNOFF1,'RUNOFF2',RUNOFF2,'RUNOFF3',RUNOFF3, &
          'RC',RC, 'PC',PC,'RSMIN',RSMIN,'XLAI',XLAI,'RCS',RCS,  &
          'RCT',RCT,'RCQ',RCQ,'RCSOIL',RCSOIL,'SOILW',SOILW,     &
          'SOILM',SOILM,'Q1',Q1,'SMCWLT',SMCWLT,'SMCDRY',SMCDRY,&
          'SMCREF',SMCREF,'SMCMAX',SMCMAX,'NROOT',NROOT
           endif
#endif

          IF (rdlai2d) THEN
             IF (SHDFAC > 0.0 .AND. lai(I) <= 0.0) lai(I) = 0.01
             xlai = lai(I)
          endif

    IF ( ICE == 1 ) THEN

       ! Sea-ice case

       DO NS = 1, NSOIL
          SH2O(NS,I) = 1.0
       ENDDO
       lai(I) = 0.01

       CYCLE ILOOP

    ELSEIF (ICE == 0) THEN

       ! Non-glacial land
!
!  FASDAS
!
       IF( fasdas == 1 ) THEN

         DZQ = DZ8W(I,1,J)
         XSDA_HFX= SDA_hfx(I,J)*RHO(I,1,J)*cpm(I)*DZQ  ! W/m^2
         ! TWG2015 Bugfix remove factor of 1000.0 for correct units
         XSDA_QFX= SDA_QFX(I,J)*RHO(I,1,J)*DZQ          ! Kg/m2/s of water
         XQNORM = QNORM(I,J)

       ENDIF
!
!  END FASDAS
!
       CALL SFLX (I,J,FFROZP, ISURBAN, DT,ZLVL,NSOIL,SLDPTH,       &    !C
                 LOCAL,                                            &    !L
                 LUTYPE, SLTYPE,                                   &    !CL
                 LWDN,SOLDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2K,DUMMY,   &    !F
                 DUMMY,DUMMY, DUMMY,                               &    !F PRCPRAIN not used
                 TH2,Q2SAT,DQSDT2,                                 &    !I
                 VEGTYP,SOILTYP,SLOPETYP,SHDFAC,SHMIN,SHMAX,       &    !I
                 ALBBRD, SNOALB1,TBOT, Z0BRD, Z0K, EMISSI, EMBRD,  &    !S
                 CMC,T1,STC,SMC,SWC,SNOWHK,SNEQV,ALBEDOK,CHK,dummy,&    !H
                 ETA,SHEAT, ETA_KINEMATIC,FDOWN,                   &    !O
                 EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,                    &    !O
                 BETA,ETP,SSOIL,                                   &    !O
                 FLX1,FLX2,FLX3,                                   &    !O
	          	 FLX4,FVB,FBUR,FGSN,UA_PHYS,                       &    !UA 
                 SNOMLT,SNCOVR,                                    &    !O
                 RUNOFF1,RUNOFF2,RUNOFF3,                          &    !O
                 RC,PC,RSMIN,XLAI,RCS,RCT,RCQ,RCSOIL,              &    !O
                 SOILW,SOILM,Q1,SMAV,                              &    !D
                 RDLAI2D,USEMONALB,                                &
                 SNOTIME1,                                         &
                 RIBB,                                             &
                 SMCWLT,SMCDRY,SMCREF,SMCMAX,NROOT,                &
                 sfcheadrt(i),                                   &    !I
                 infxsrt(i),ETPND1,OPT_THCND,AOASIS              &    !O
                ,XSDA_QFX, HFX_PHY, QFX_PHY, XQNORM, fasdas, HCPCT_FASDAS    &   ! fasdas
                 )

#ifdef WRF_HYDRO
                 soldrain(i) = RUNOFF2*DT*1000.0
#endif
    ELSEIF (ICE == -1) THEN

       !
       ! Set values that the LSM is expected to update,
       ! but don't get updated for glacial points.
       !
       SOILM = 0.0 !BSINGH(PNNL)- SOILM is undefined for this case, it is used for diagnostics so setting it to zero
       XLAI = 0.01 ! KWM Should this be Zero over land ice?  Does this value matter?
       RUNOFF2 = 0.0
       RUNOFF3 = 0.0
       DO NS = 1, NSOIL
          SWC(NS) = 1.0
          SMC(NS) = 1.0
          SMAV(NS) = 1.0
       ENDDO
!
!  FASDAS
!
       IF( fasdas == 1 ) THEN

         DZQ = DZ8W(I,1,J)
         XSDA_HFX= SDA_HFX(I,J)*RHO(I,1,J)*cpm(I)*DZQ  ! W/m^2
         XSDA_QFX= 0.0          ! Kg/m2/s of water
         XQNORM = 0.0

       ENDIF
!
!  END FASDAS
!

       CALL SFLX_GLACIAL(I,J,ISICE,FFROZP,DT,ZLVL,NSOIL,SLDPTH,   &    !C
            &    LWDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2K,              &    !F
            &    TH2,Q2SAT,DQSDT2,                                &    !I
            &    ALBBRD, SNOALB1,TBOT, Z0BRD, Z0K, EMISSI, EMBRD, &    !S
            &    T1,STC(1:NSOIL),SNOWHK,SNEQV,ALBEDOK,CHK,        &    !H
            &    ETA,SHEAT,ETA_KINEMATIC,FDOWN,                   &    !O
            &    ESNOW,DEW,                                       &    !O
            &    ETP,SSOIL,                                       &    !O
            &    FLX1,FLX2,FLX3,                                  &    !O
            &    SNOMLT,SNCOVR,                                   &    !O
            &    RUNOFF1,                                         &    !O
            &    Q1,                                              &    !D
            &    SNOTIME1,                                        &
            &    RIBB)

    ENDIF

       lai(I) = xlai
       if (present(rc2) .and. present(xlai2)) then
       RC2(i) = RC                                              ! for output
       XLAI2(i) = XLAI
       endif

#if 0
          IF(IPRINT) THEN

       print*, 'AFTER SFLX, in Noahlsm_driver'
       print*, 'ICE', ICE, 'DT',DT, 'ZLVL',ZLVL, 'NSOIL', NSOIL,   &
       'SLDPTH', SLDPTH, 'LOCAL',LOCAL, 'LUTYPE',&
        LUTYPE, 'SLTYPE',SLTYPE, 'LWDN',LWDN, 'SOLDN',SOLDN,      &
        'SFCPRS',SFCPRS, 'PRCP',PRCP,'SFCTMP',SFCTMP,'Q2K',Q2K,   &
         'TH2',TH2,'Q2SAT',Q2SAT,'DQSDT2',DQSDT2,'VEGTYP', VEGTYP,&
          'SOILTYP',SOILTYP, 'SLOPETYP',SLOPETYP, 'SHDFAC',SHDFAC,&
         'SHDMIN',SHMIN, 'ALBBRD',ALBBRD,'SNOALB',SNOALB1,'TBOT',&
          TBOT, 'Z0BRD',Z0BRD, 'Z0K',Z0K, 'CMC',CMC, 'T1',T1,'STC',&
          STC, 'SMC',SMC, 'SWc',SWC,'SNOWHK',SNOWHK,'SNEQV',SNEQV,&
          'ALBEDOK',ALBEDOK,'CHK',CHK,'ETA',ETA,'SHEAT',SHEAT,      &
          'ETA_KINEMATIC',ETA_KINEMATIC, 'FDOWN',FDOWN,'EC',EC,   &
          'EDIR',EDIR,'ET',ET,'ETT',ETT,'ESNOW',ESNOW,'DRIP',DRIP,&
          'DEW',DEW,'BETA',BETA,'ETP',ETP,'SSOIL',SSOIL,'FLX1',FLX1,&
          'FLX2',FLX2,'FLX3',FLX3,'SNOMLT',SNOMLT,'SNCOVR',SNCOVR,&
          'RUNOFF1',RUNOFF1,'RUNOFF2',RUNOFF2,'RUNOFF3',RUNOFF3, &
          'RC',RC, 'PC',PC,'RSMIN',RSMIN,'XLAI',XLAI,'RCS',RCS,  &
          'RCT',RCT,'RCQ',RCQ,'RCSOIL',RCSOIL,'SOILW',SOILW,     &
          'SOILM',SOILM,'Q1',Q1,'SMCWLT',SMCWLT,'SMCDRY',SMCDRY,&
          'SMCREF',SMCREF,'SMCMAX',SMCMAX,'NROOT',NROOT
           endif
#endif

!***  UPDATE STATE VARIABLES
          canwat(I)=CMC*1000.
          snow(I)=SNEQV*1000.
!          snowh(I)=SNOWHK*1000.
          snowh(I)=SNOWHK                   ! SNOWHK in meters
          albedo(I)=ALBEDOK
          ALB_RURAL(i)=ALBEDOK
          albbck(I)=ALBBRD
          z0(I)=Z0BRD
          emiss(I) = EMISSI
          emiss_rural(i) = EMISSI
! Noah: activate time-varying roughness length (V3.3 Feb 2011)
          znt(I)=Z0K
!
! FASDAS
!
! Update Skin Temperature

       IF( fasdas == 1 ) THEN
           XSDA_QFX= XSDA_QFX*ELWV*XQNORM

        !TWG2015 Bugfix to multiply Heat Capacity by Soil Depth for correct
        !units 

           T1 = T1 + (XSDA_HFX-XSDA_QFX)*DT/(HCPCT_FASDAS*DZS(1))

       END IF
!
! END FASDAS
!
          tsk(I)=T1
          TSK_RURAL(I)=T1
          if (present(tsk_rural_bep)) then
             IF(SF_URBAN_PHYSICS == 2 .or. SF_URBAN_PHYSICS == 3) THEN
               TSK_RURAL_BEP(I)=T1
             END IF
          endif
          hfx(I)=SHEAT
          HFX_RURAL(i)=SHEAT
! MEk Jul07 add potential evap accum
        potevp(I)=potevp(I)+ETP*FDTW
          qfx(I)=ETA_KINEMATIC
          QFX_RURAL(i)=ETA_KINEMATIC

#ifdef WRF_HYDRO
!added by Wei Yu
!         qfx(I) = qfx(I) + ETPND1
!         ETA = ETA + ETPND1/2.501E6*dt
!end added by Wei Yu
#endif


          lh(I)=ETA
          LH_RURAL(i)=ETA
          grdflx(I)=SSOIL
          GRDFLX_RURAL(i)=SSOIL
          snowc(I)=SNCOVR
          chs2(I)=cqs2(I)
          snotime(I) = SNOTIME1
!      prevent diagnostic ground q (q1) from being greater than qsat(tsk)
!      as happens over snow cover where the cqs2 value also becomes irrelevant
!      by setting cqs2=chs in this situation the 2m q should become just qv(k=1)
! ww: comment out this change to avoid Q2 drop due to change of radiative flux
!         IF (Q1 .GT. qsfc(I)) THEN
!           cqs2(I) = chs(I)
!         ENDIF
!          qsfc(I)=Q1
! Convertqsfc(I)back to mixing ratio
           qsfc(I)= Q1/(1.0-Q1)
!
           ! QSFC_RURAL(I,J)= Q1/(1.0-Q1)
! Calculate momentum flux from rural surface for use with multi-layer UCM (Martilli et al. 2002)

          DO 80 NS=1,NSOIL
           smois(NS,I)=SMC(NS)
           tslb(NS,I)=STC(NS)                                        !  STEMP
           SH2O(NS,I)=SWC(NS)
   80     CONTINUE
!       ENDIF

        FLX4_2D(i)  = FLX4
	fvb_2d(i)   = FVB
	fbur_2d(i)  = FBUR
	fgsn_2d(i)  = FGSN

     !
     ! Residual of surface energy balance equation terms
     !

     IF ( UA_PHYS ) THEN
         noahres(I) = ( solnet + lwdn ) - sheat + ssoil - eta &
              - ( emissi * STBOLT * (t1**4) ) - flx1 - flx2 - flx3 - flx4

     ELSE
         noahres(I) = ( solnet + lwdn ) - sheat + ssoil - eta &
              - ( emissi * STBOLT * (t1**4) ) - flx1 - flx2 - flx3
     ENDIF

        IF (SF_URBAN_PHYSICS == 1 ) THEN                                              ! Beginning of UCM CALL if block
!--------------------------------------
! URBAN CANOPY MODEL START - urban
!--------------------------------------
! Input variables lsm --> urban

!qkyang replaced to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!                 IF( ivgtyp(I) == ISURBAN .or. ivgtyp(I) == LOW_DENSITY_RESIDENTIAL .or. &
!                   ivgtyp(I) == HIGH_DENSITY_RESIDENTIAL .or. ivgtyp(I) == HIGH_INTENSITY_INDUSTRIAL) THEN
               IF(IVGTYP(I) == ISURBAN .or. &
                  IVGTYP(I) == low_density_residential .or. &  
                  IVGTYP(I) == HIGH_DENSITY_RESIDENTIAL .or. &
                  IVGTYP(I) == HIGH_INTENSITY_INDUSTRIAL .or. &
                  IVGTYP(I) == LCZ_1 .or. IVGTYP(I) == LCZ_2 .or. &
                  IVGTYP(I) == LCZ_3 .or. IVGTYP(I) == LCZ_4 .or. &
                  IVGTYP(I) == LCZ_5 .or. IVGTYP(I) == LCZ_6 .or. &
                  IVGTYP(I) == LCZ_7 .or. IVGTYP(I) == LCZ_8 .or. &
                  IVGTYP(I) == LCZ_9 .or. IVGTYP(I) == LCZ_10 .or. &
                  IVGTYP(I) == LCZ_11 )  then !urban
               !qkyang replaced to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! Call urban

!
            UTYPE_URB = utype_urb2d(i) !urban type (low, high or industrial)
            TA_URB    = SFCTMP           ! [K]
            QA_URB    = Q2K              ! [kg/kg]
            UA_URB    = SQRT(U_PHY(I,1,J)**2.+V_PHY(I,1,J)**2.)
            U1_URB    = U_PHY(I,1,J)
            V1_URB    = V_PHY(I,1,J)
            IF(UA_URB < 1.) UA_URB=1.    ! [m/s]
            SSG_URB   = SOLDN            ! [W/m/m]
            SSGD_URB  = 0.8*SOLDN        ! [W/m/m]
            SSGQ_URB  = SSG_URB-SSGD_URB ! [W/m/m]
            LLG_URB   = glw(I)         ! [W/m/m]
            RAIN_URB  = rainbl(I)      ! [mm]
            RHOO_URB  = SFCPRS / (287.04 * SFCTMP * (1.0+ 0.61 * Q2K)) ![kg/m/m/m]
            ZA_URB    = ZLVL             ! [m]
            DELT_URB  = DT               ! [sec]
            XLAT_URB  = XLAT_URB2D(I)  ! [deg]
            COSZ_URB  = COSZ_URB2D(i)  !
            OMG_URB   = OMG_URB2D(i)   !
            ZNT_URB   = znt(I)
            LSOLAR_URB = .FALSE.
            TR_URB = TR_URB2D(i)
            TB_URB = TB_URB2D(i)
            TG_URB = TG_URB2D(i)
            TC_URB = TC_URB2D(i)
            QC_URB = QC_URB2D(i)
            UC_URB = UC_URB2D(i)
            TGR_URB   = TGR_URB2D(i)
            CMCR_URB  = CMCR_URB2D(i)
            FLXHUMR_URB = FLXHUMR_URB2D(i)
            FLXHUMB_URB = FLXHUMB_URB2D(i)
            FLXHUMG_URB = FLXHUMG_URB2D(i)
            DRELR_URB = DRELR_URB2D(i)
            DRELB_URB = DRELB_URB2D(i)
            DRELG_URB = DRELG_URB2D(i)

            DO K = 1,num_roof_layers
              TRL_URB(K) = TRL_URB3D(K,I)
              SMR_URB(K) = SMR_URB3D(K,I)
              TGRL_URB(K)= TGRL_URB3D(K,I)
            END DO
            DO K = 1,num_wall_layers
              TBL_URB(K) = TBL_URB3D(K,I)
            END DO
            DO K = 1,num_road_layers
              TGL_URB(K) = TGL_URB3D(K,I)
            END DO

            XXXR_URB = XXXR_URB2D(i)
            XXXB_URB = XXXB_URB2D(i)
            XXXG_URB = XXXG_URB2D(i)
            XXXC_URB = XXXC_URB2D(i)
!
!
!      Limits to avoid dividing by small number
            if (chs(I) < 1.0E-02) then
               chs(I)  = 1.0E-02
            endif
            if (chs2(I) < 1.0E-02) then
               chs2(I)  = 1.0E-02
            endif
            if (cqs2(I) < 1.0E-02) then
               cqs2(I)  = 1.0E-02
            endif
!
            CHS_URB  = chs(I)
            CHS2_URB = chs2(I)
            IF (PRESENT(CMR_SFCDIF)) THEN
               CMR_URB = CMR_SFCDIF(i)
               CHR_URB = CHR_SFCDIF(i)
               CMGR_URB = CMGR_SFCDIF(i)
               CHGR_URB = CHGR_SFCDIF(i)
               CMC_URB = CMC_SFCDIF(i)
               CHC_URB = CHC_SFCDIF(i)
            ENDIF
! NUDAPT for SLUCM
            mh_urb = MH_URB2D(i)
            stdh_urb = STDH_URB2D(i)
            lp_urb = LP_URB2D(i)
            hgt_urb = HGT_URB2D(i)
            lf_urb = 0.0
            DO K = 1,4
              lf_urb(K)=lf_urb2d(K,I)
            ENDDO
            frc_urb = frc_urb2d(i)
            lb_urb = LB_URB2D(i)
            check = 0
            if (I.eq.73.and.J.eq.125)THEN
               check = 1
            end if
!
! Call urban
            !CALL cal_mon_day(julday,julyr,jmonth,jday)!qkyang replace julian with julday
            CALL urban(LSOLAR_URB,                                      & ! I
                       num_roof_layers,num_wall_layers,num_road_layers, & ! C
                       DZR,DZB,DZG,                                     & ! C
                       UTYPE_URB,TA_URB,QA_URB,UA_URB,U1_URB,V1_URB,SSG_URB, & ! I
                       SSGD_URB,SSGQ_URB,LLG_URB,RAIN_URB,RHOO_URB,     & ! I
                       ZA_URB,DECLIN_URB,COSZ_URB,OMG_URB,              & ! I
                       XLAT_URB,DELT_URB,ZNT_URB,                       & ! I
                       CHS_URB, CHS2_URB,                               & ! I
                       TR_URB, TB_URB, TG_URB, TC_URB, QC_URB,UC_URB,   & ! H
                       TRL_URB,TBL_URB,TGL_URB,                         & ! H
                       XXXR_URB, XXXB_URB, XXXG_URB, XXXC_URB,          & ! H
                       TS_URB,QS_URB,SH_URB,LH_URB,LH_KINEMATIC_URB,    & ! O
                       SW_URB,ALB_URB,LW_URB,G_URB,RN_URB,PSIM_URB,PSIH_URB, & ! O
                       GZ1OZ0_URB,                                      & !O
                       CMR_URB, CHR_URB, CMC_URB, CHC_URB,              &
                       U10_URB, V10_URB, TH2_URB, Q2_URB,               & ! O
                       UST_URB,mh_urb, stdh_urb, lf_urb, lp_urb,        & ! 0
                       hgt_urb,frc_urb,lb_urb, check,CMCR_URB,TGR_URB,  & ! H
                       TGRL_URB,SMR_URB,CMGR_URB,CHGR_URB,jmonth,       & ! H
                       DRELR_URB,DRELB_URB,                             & ! H
                       DRELG_URB,FLXHUMR_URB,FLXHUMB_URB,FLXHUMG_URB,   &
                       !qkyang added to support anthropogenic heating emissions+++++++++
                       AH_emis(i),is_ah_emis,AH_URB2D(i))
                       !qkyang added to support anthropogenic heating emissions--------

#if 0
          IF(IPRINT) THEN

       print*, 'AFTER CALL URBAN'
       print*,'num_roof_layers',num_roof_layers, 'num_wall_layers',  &
        num_wall_layers,                                             &
       'DZR',DZR,'DZB',DZB,'DZG',DZG,'UTYPE_URB',UTYPE_URB,'TA_URB', &
        TA_URB,                                                      &
        'QA_URB',QA_URB,'UA_URB',UA_URB,'U1_URB',U1_URB,'V1_URB',    &
         V1_URB,                                                     &
         'SSG_URB',SSG_URB,'SSGD_URB',SSGD_URB,'SSGQ_URB',SSGQ_URB,  &
        'LLG_URB',LLG_URB,'RAIN_URB',RAIN_URB,'RHOO_URB',RHOO_URB,   &
        'ZA_URB',ZA_URB, 'DECLIN_URB',DECLIN_URB,'COSZ_URB',COSZ_URB,&
        'OMG_URB',OMG_URB,'XLAT_URB',XLAT_URB,'DELT_URB',DELT_URB,   &
         'ZNT_URB',ZNT_URB,'TR_URB',TR_URB, 'TB_URB',TB_URB,'TG_URB',&
         TG_URB,'TC_URB',TC_URB,'QC_URB',QC_URB,'TRL_URB',TRL_URB,   &
          'TBL_URB',TBL_URB,'TGL_URB',TGL_URB,'XXXR_URB',XXXR_URB,   &
         'XXXB_URB',XXXB_URB,'XXXG_URB',XXXG_URB,'XXXC_URB',XXXC_URB,&
         'TS_URB',TS_URB,'QS_URB',QS_URB,'SH_URB',SH_URB,'LH_URB',   &
         LH_URB, 'LH_KINEMATIC_URB',LH_KINEMATIC_URB,'SW_URB',SW_URB,&
         'ALB_URB',ALB_URB,'LW_URB',LW_URB,'G_URB',G_URB,'RN_URB',   &
          RN_URB, 'PSIM_URB',PSIM_URB,'PSIH_URB',PSIH_URB,          &
         'U10_URB',U10_URB,'V10_URB',V10_URB,'TH2_URB',TH2_URB,      &
          'Q2_URB',Q2_URB,'CHS_URB',CHS_URB,'CHS2_URB',CHS2_URB
           endif
#endif

            TS_URB2D(i) = TS_URB

            albedo(I) = frc_urb2d(i)*ALB_URB+(1-frc_urb2d(i))*ALBEDOK   ![-]
            hfx(I) = frc_urb2d(i)*SH_URB+(1-frc_urb2d(i))*SHEAT         ![W/m/m]
            qfx(I) = frc_urb2d(i)*LH_KINEMATIC_URB &
                     + (1-frc_urb2d(i))*ETA_KINEMATIC                ![kg/m/m/s]
            lh(I) = frc_urb2d(i)*LH_URB+(1-frc_urb2d(i))*ETA            ![W/m/m]
            grdflx(I) = frc_urb2d(i)*G_URB+(1-frc_urb2d(i))*SSOIL       ![W/m/m]
            tsk(I) = frc_urb2d(i)*TS_URB+(1-frc_urb2d(i))*T1            ![K]
            Q1 = frc_urb2d(i)*QS_URB+(1-frc_urb2d(i))*Q1            ![-]
! Convert QSFC back to mixing ratio
            qsfc(I)= Q1/(1.0-Q1)
            ust(i)= frc_urb2d(i)*UST_URB+(1-frc_urb2d(i))*ust(i)      ![m/s]

#if 0
    IF(IPRINT)THEN

    print*, ' FRC_URB2D', FRC_URB2D,                        &
    'ALB_URB',ALB_URB, 'ALBEDOK',ALBEDOK, &
    'albedo(I)',  albedo(I),                  &
    'SH_URB',SH_URB,'SHEAT',SHEAT, 'hfx(I)',hfx(I),  &
    'LH_KINEMATIC_URB',LH_KINEMATIC_URB,'ETA_KINEMATIC',  &
     ETA_KINEMATIC, 'qfx(I)',qfx(I),                  &
    'LH_URB',LH_URB, 'ETA',ETA, 'lh(I)',lh(I),        &
    'G_URB',G_URB,'SSOIL',SSOIL,'grdflx(I)', grdflx(I),&
    'TS_URB',TS_URB,'T1',T1,'tsk(I)',tsk(I),          &
    'QS_URB',QS_URB,'Q1',Q1,'qsfc(I)',qsfc(I)
     endif
#endif



! Renew Urban State Varialbes

            TR_URB2D(i) = TR_URB
            TB_URB2D(i) = TB_URB
            TG_URB2D(i) = TG_URB
            TC_URB2D(i) = TC_URB
            QC_URB2D(i) = QC_URB
            UC_URB2D(i) = UC_URB

            TGR_URB2D(i) =TGR_URB 
            CMCR_URB2D(i)=CMCR_URB
            FLXHUMR_URB2D(i)=FLXHUMR_URB 
            FLXHUMB_URB2D(i)=FLXHUMB_URB
            FLXHUMG_URB2D(i)=FLXHUMG_URB
            DRELR_URB2D(i) = DRELR_URB 
            DRELB_URB2D(i) = DRELB_URB 
            DRELG_URB2D(i) = DRELG_URB

            DO K = 1,num_roof_layers
              TRL_URB3D(K,I) = TRL_URB(K)
              SMR_URB3D(K,I) = SMR_URB(K)
              TGRL_URB3D(K,I)= TGRL_URB(K)
            END DO
            DO K = 1,num_wall_layers
              TBL_URB3D(K,I) = TBL_URB(K)
            END DO
            DO K = 1,num_road_layers
              TGL_URB3D(K,I) = TGL_URB(K)
            END DO
            XXXR_URB2D(i) = XXXR_URB
            XXXB_URB2D(i) = XXXB_URB
            XXXG_URB2D(i) = XXXG_URB
            XXXC_URB2D(i) = XXXC_URB

            SH_URB2D(i)    = SH_URB
            LH_URB2D(i)    = LH_URB
            G_URB2D(i)     = G_URB
            RN_URB2D(i)    = RN_URB
            PSIM_URB2D(i)  = PSIM_URB
            PSIH_URB2D(i)  = PSIH_URB
            GZ1OZ0_URB2D(i)= GZ1OZ0_URB
            U10_URB2D(i)   = U10_URB
            V10_URB2D(i)   = V10_URB
            TH2_URB2D(i)   = TH2_URB
            Q2_URB2D(i)    = Q2_URB
            ust_urb2d(i)   = UST_URB
            AKMS_URB2D(i)  = KARMAN * ust_urb2d(i)/(GZ1OZ0_URB2D(i)-PSIM_URB2D(i))
            IF (PRESENT(CMR_SFCDIF)) THEN
               CMR_SFCDIF(i) = CMR_URB
               CHR_SFCDIF(i) = CHR_URB
               CMGR_SFCDIF(i) = CMGR_URB
               CHGR_SFCDIF(i) = CHGR_URB
               CMC_SFCDIF(i) = CMC_URB
               CHC_SFCDIF(i) = CHC_URB
            ENDIF
          END IF

         ENDIF                                   ! end of UCM CALL if block

!--------------------------------------
! Urban Part End - urban
!--------------------------------------

!***  DIAGNOSTICS
          smstav(I)=SOILW
          smstot(I)=SOILM*1000.
          DO NS=1,NSOIL
          smcrel(NS,I)=SMAV(NS)
          ENDDO
!         Convert the water unit into mm
          sfcrunoff(I)=sfcrunoff(I)+RUNOFF1*DT*1000.0
          udrunoff(I)=udrunoff(I)+RUNOFF2*DT*1000.0
! snow defined when fraction of frozen precip (FFROZP) > 0.5,
          IF(FFROZP.GT.0.5)THEN
            acsnow(I)=acsnow(I)+PRCP*DT
          ENDIF
          IF(snow(I).GT.0.)THEN
            acsnom(I)=acsnom(I)+SNOMLT*1000.
! accumulated snow-melt energy
            snopcx(I)=snopcx(I)-SNOMLT/FDTLIW
          ENDIF

        ENDIF                                                           ! endif of land-sea test

      ENDDO ILOOP                                                       ! of I loop
   ENDDO JLOOP                                                          ! of J loop

      IF (SF_URBAN_PHYSICS == 2) THEN

         do j=jts,jte
         do i=its,ite
            emiss_urb(i)=0.
            rl_up_urb(i)=0.
            RS_ABS_URB(i)=0.
            GRDFLX_URB(i)=0.
            b_q_bep(kms:kme,i)=0.
         end do
         end do
       CALL BEP(frc_urb2d,utype_urb2d,itimestep,dz8w,dt,u_phy,v_phy,   &
                th_phy,rho,p_phy,swdown,glw,                           &
                gmt,julday,xlong,xlat,declin_urb,cosz_urb2d,omg_urb2d, &
                num_urban_layers,num_urban_hi,                         &
                trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,               &
                sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,             &
                lp_urb2d,hi_urb2d,lb_urb2d,hgt_urb2d,                  &
                a_u_bep,a_v_bep,a_t_bep,                               &
                a_e_bep,b_u_bep,b_v_bep,                               &
                b_t_bep,b_e_bep,b_q_bep,dlg_bep,                       &
                dl_u_bep,sf_bep,vl_bep,                                &
                rl_up_urb,rs_abs_urb,emiss_urb,grdflx_urb,             &
                ids,ide, jds,jde, kds,kde,                             &
                ims,ime, jms,jme, kms,kme,                             &
                its,ite, jts,jte, kts,kte )

       ENDIF

       
       IF (SF_URBAN_PHYSICS == 3) THEN


         do j=jts,jte
         do i=its,ite
            emiss_urb(i)=0.
            rl_up_urb(i)=0.
            RS_ABS_URB(i)=0.
            GRDFLX_URB(i)=0.
            b_q_bep(kms:kme,i)=0.
         end do
         end do
        

       CALL BEP_BEM(frc_urb2d,utype_urb2d,itimestep,dz8w,dt,u_phy,v_phy, &
                th_phy,rho,p_phy,swdown,glw,                           &
                gmt,julday,xlong,xlat,declin_urb,cosz_urb2d,omg_urb2d, &
                num_urban_layers,num_urban_hi,                         &
                trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,               &
                tlev_urb3d,qlev_urb3d,tw1lev_urb3d,tw2lev_urb3d,       &
                tglev_urb3d,tflev_urb3d,sf_ac_urb3d,lf_ac_urb3d,       &
                cm_ac_urb3d,sfvent_urb3d,lfvent_urb3d,                 &
                sfwin1_urb3d,sfwin2_urb3d,                             &
                sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,             &
                lp_urb2d,hi_urb2d,lb_urb2d,hgt_urb2d,                  &
                a_u_bep,a_v_bep,a_t_bep,                               &
                a_e_bep,b_u_bep,b_v_bep,                               &
                b_t_bep,b_e_bep,b_q_bep,dlg_bep,                       &
                dl_u_bep,sf_bep,vl_bep,                                &
                rl_up_urb,rs_abs_urb,emiss_urb,grdflx_urb,qv3d,        &
                ids,ide, jds,jde, kds,kde,                             &
                ims,ime, jms,jme, kms,kme,                             &
                its,ite, jts,jte, kts,kte )

       ENDIF

    if((sf_urban_physics.eq.2).OR.(sf_urban_physics.eq.3))then         !Bep begin
         do j=jts,jte
         do i=its,ite
            UMOM_URB(i)=0.
            VMOM_URB(i)=0.
            HFX_URB(i)=0.
            QFX_URB(i)=0.
         do k=kts,kte
            A_U_BEP(k,i)=A_U_BEP(k,i)*frc_urb2d(i)
            A_V_BEP(k,i)=A_V_BEP(k,i)*frc_urb2d(i)
            A_T_BEP(k,i)=A_T_BEP(k,i)*frc_urb2d(i)
            A_Q_BEP(k,i)=0.        
            A_E_BEP(k,i)=0.
            B_U_BEP(k,i)=B_U_BEP(k,i)*frc_urb2d(i)
            B_V_BEP(k,i)=B_V_BEP(k,i)*frc_urb2d(i)
            B_T_BEP(k,i)=B_T_BEP(k,i)*frc_urb2d(i)
            B_Q_BEP(k,i)=B_Q_BEP(k,i)*frc_urb2d(i)
            B_E_BEP(k,i)=B_E_BEP(k,i)*frc_urb2d(i)
            HFX_URB(i)=HFX_URB(i)+B_T_BEP(K,i)*RHO(I,K,J)*CP*                &
                          DZ8W(I,K,J)*VL_BEP(K,i)
            QFX_URB(i)=QFX_URB(i)+B_Q_BEP(K,i)*               &
                          DZ8W(I,K,J)*VL_BEP(K,i)
            UMOM_URB(i)=UMOM_URB(i)+ (A_U_BEP(k,i)*U_PHY(I,K,J)+             &
                          B_U_BEP(K,i))*DZ8W(I,K,J)*VL_BEP(K,i)
            VMOM_URB(i)=VMOM_URB(i)+ (A_V_BEP(K,i)*V_PHY(I,K,J)+             &
                          B_V_BEP(K,i))*DZ8W(I,K,J)*VL_BEP(K,i)
            VL_BEP(k,i)=(1.-frc_urb2d(i))+VL_BEP(k,i)*frc_urb2d(i)
            SF_BEP(k,i)=(1.-frc_urb2d(i))+SF_BEP(k,i)*frc_urb2d(i)
         end do
            A_U_BEP(1,i)=(1.-frc_urb2d(i))*(-ust(i)*ust(i))/dz8w(i,1,j)/   &
                           ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+A_U_BEP(1,i)
            A_V_BEP(1,i)=(1.-frc_urb2d(i))*(-ust(i)*ust(i))/dz8w(i,1,j)/   &
                           ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+A_V_BEP(1,i)
            B_T_BEP(1,i)=(1.-frc_urb2d(i))*HFX_RURAL(i)/dz8w(i,1,j)/rho(i,1,j)/CP+ &
                            B_T_BEP(1,i)
            B_Q_BEP(1,i)=(1.-frc_urb2d(i))*QFX_RURAL(i)/dz8w(i,1,j)/rho(i,1,j)+B_Q_BEP(1,i)
            umom=(1.-frc_urb2d(i))*ust(i)*ust(i)*u_phy(i,1,j)/               &
                         ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+UMOM_URB(i)
            vmom=(1.-frc_urb2d(i))*ust(i)*ust(i)*v_phy(i,1,j)/               &
                         ((u_phy(i,1,j)**2+v_phy(i,1,j)**2.)**.5)+VMOM_URB(i)
            SF_BEP(1,i)=1.
           
! compute upward longwave radiation from the rural part and total
!           rl_up_rural=-emiss_rural(i)*sigma_sb*(TSK_RURAL(I)**4.)-(1.-emiss_rural(i))*glw(I)
!           rl_up_tot=(1.-frc_urb2d(i))*rl_up_rural+frc_urb2d(i)*rl_up_urb(i)   
!           emiss(I)=(1.-frc_urb2d(i))*emiss_rural(i)+frc_urb2d(i)*emiss_urb(i)
! using the emissivity and the total longwave upward radiation estimate the averaged skin temperature
           IF (frc_urb2d(i).GT.0.) THEN
              rl_up_rural=-emiss_rural(i)*sigma_sb*(TSK_RURAL(I)**4.)-(1.-emissi)*glw(I)
              rl_up_tot=(1.-frc_urb2d(i))*rl_up_rural+frc_urb2d(i)*rl_up_urb(i)   
              emiss(I)=(1.-frc_urb2d(i))*emiss_rural(i)+frc_urb2d(i)*emiss_urb(i)
              TS_URB2D(i)=(max(0.,(-rl_up_urb(i)-(1.-emiss_urb(i))*glw(I))/emiss_urb(i)/sigma_sb))**0.25
              tsk(I)=(max(0., (-1.*rl_up_tot-(1.-emiss(I))*glw(I) )/emiss(I)/sigma_sb))**.25
           rs_abs_tot=(1.-frc_urb2d(i))*swdown(I)*(1.-albedo(I))+frc_urb2d(i)*RS_ABS_URB(i)  
          if(swdown(I).gt.0.)then
           albedo(I)=1.-rs_abs_tot/swdown(I)
          else
           albedo(I)=ALB_RURAL(i)
          endif
! rename *_urb to sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d
         grdflx(I)= (1.-frc_urb2d(i))*GRDFLX_RURAL(i)+frc_urb2d(i)*GRDFLX_URB(i) 
         qfx(I)=(1.-frc_urb2d(i))*QFX_RURAL(i)+QFX_URB(i)
!        lh(I)=(1.-frc_urb2d(i))*QFX_RURAL(i)*xlv
         lh(I)=qfx(I)*xlv
         hfx(I) = HFX_URB(i)+(1-frc_urb2d(i))*HFX_RURAL(i)         ![W/m/m]
            SH_URB2D(i)    = HFX_URB(i)/frc_urb2d(i)
            LH_URB2D(i)    = QFX_URB(i)*xlv
            G_URB2D(i)     = GRDFLX_URB(i)
            RN_URB2D(i)    = RS_ABS_URB(i)+emiss_urb(i)*glw(I)-rl_up_urb(i)
            ust(i)=(umom**2.+vmom**2.)**.25
!          if(tsk(I).gt.350)write(*,*)'tsk too big!',i,j,tsk(I)
!          if(tsk(I).lt.260)write(*,*)'tsk too small!',i,j,tsk(I),rl_up_tot,rl_up_urb(i),rl_up_rural
!                    print*,'ivgtyp(I),i,j,sigma_sb',ivgtyp(I),i,j,sigma_sb
!                    print*,'hfx,lh,qfx,grdflx,ts_urb2d',hfx(I),lh(I),qfx(I),grdflx(I),TS_URB2D(i)
!                    print*,'tsk,albedo,emiss',tsk(I),albedo(I),emiss(I)
!                if(i.eq.56.and.j.eq.29)then
!                    print*,'ivgtyp(I), qfx, hfx',ivgtyp(I),HFX_RURAL(i),QFX_RURAL(i)
!                    print*,'emiss_rural,emiss_urb',emiss_rural(i),emiss_urb(i)
!                    print*,'rl_up_rural,rl_up_urb(i)',rl_up_rural,rl_up_urb(i)
!                    print*,'tsk_rural,TS_URB2D(i),tsk',TSK_RURAL(I),TS_URB2D(i),tsk(I)
!                    print*,'reconstruction fei',((emiss(I)*tsk(I)**4.-frc_urb2d(i)*emiss_urb(i)*TS_URB2D(i)**4.)/(emiss_rural(i)*(1.-frc_urb2d(i))))**.25
!                    print*,'ivgtyp(I),hfx,hfx_urb,hfx_rural',hfx(I),HFX_URB(i),HFX_RURAL(i)
!                    print*,'lh,lh_rural',lh(I),LH_RURAL(i)
!                    print*,'qfx',qfx(I)
!                    print*,'ts_urb2d',TS_URB2D(i)
!                    print*,'ust',ust(i)
!                    print*,'swdown,glw',swdown(I),glw(I)
!                endif
            else
              SH_URB2D(i)    = 0.
              LH_URB2D(i)    = 0.
              G_URB2D(i)     = 0.
              RN_URB2D(i)    = 0.
            endif
!          IF( ivgtyp(I) == 1 .or. ivgtyp(I) == LOW_DENSITY_RESIDENTIAL .or. &
!                  ivgtyp(I) == HIGH_DENSITY_RESIDENTIAL .or. ivgtyp(I) == HIGH_INTENSITY_INDUSTRIAL) THEN
!                    print*,'ivgtyp(I), qfx, hfx',ivgtyp(I),HFX_RURAL(i),QFX_RURAL(i)
!                    print*,'ivgtyp,hfx,hfx_urb,hfx_rural',hfx(I),HFX_URB(i),HFX_RURAL(i)
!                    print*,'lh,lh_rural',lh(I),LH_RURAL(i)
!                    print*,'qfx',qfx(I)
!                    print*,'ts_urb2d',TS_URB2D(i)
!                    print*,'ust',ust(i)
!          endif
        enddo
        enddo

       endif                                                           !Bep end

!------------------------------------------------------
   END SUBROUTINE lsm
!------------------------------------------------------

!For iAMAS, the below section of the module is moved to module_physics_lsm_noahinit.F in
!./../core_physics to accomodate differences in the mpi calls between WRF and iAMAS.I thought
!that it would be cleaner to do this instead of adding a lot of #ifdef statements throughout
!the initialization subroutine.

!qkyang delete LSMINIT, lsm_mosaic, lsm_mosaic_init, SOIL_VEG_GEN_PARM that were not used

END MODULE module_sf_noahdrv
