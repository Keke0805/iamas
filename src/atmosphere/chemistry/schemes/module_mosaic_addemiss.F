!************************************************************************
!
! * Author: Zihan Xia, Jiawang Feng, Chen Jin and Chun Zhao from LACAR USTC
! Last update: 04/2025
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao
!
! References: 
!   Copied from WRF-Chem-v3.5:chem/module_mosaic_addemiss.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
module module_mosaic_addemiss

   USE iamas_atmchem_vars
   USE iamas_atmphys_vars
   USE iamas_kind_types
   USE iamas_pool_routines
   USE module_data_binaer_asect
   use iamas_atmchem_gocart_dust_init, only: has_done_gocart_init
   use iamas_atmchem_gocart_dust_init, only: gocart_dust_init_foriAMAS

   USE iamas_atmchem_utilities
 
contains

subroutine mosaic_addemiss_nature(configs,slai,ust,smois,ivgtyp,vegfra,isltyp, &
            dtstep, u10, v10, alt, dz8w, xland, sstemp,                        &
            ktau,num_soil_layers,u_phy,v_phy,rho_phy,p8w,snowc,erod,           &
            xlat,xlong,dx,g,emis_dust_factor,emis_seas_factor,                 &
            tend_emis, e_dust, e_seas,                                         &
            dust_emiss_active, seasalt_emiss_active,                           &
            do_seas_sst_adjust,                                                &
            soilctop,landusef,sh2o,xlai,                                       &
            its,ite, kts,kte                                                   )

       implicit none

       type(iamas_pool_type),intent(in):: configs

       character(len=strkind), intent(in   ) :: dust_emiss_active
       integer, intent(in   )                :: seasalt_emiss_active           
       logical, intent(in   )                :: do_seas_sst_adjust
       integer, intent(in   )                :: ktau, num_soil_layers
       integer, intent(in   )                :: its,ite, kts,kte

       real, intent(in   ) :: dtstep,emis_dust_factor,emis_seas_factor,g

    ! 10-m wind speed components (m/s)
       real,  dimension( its:ite),                                             &
              intent(in   ) ::   u10,v10,xland,slai,ust,dx,xlat,xlong,xlai

       real,  dimension( 1:4, its:ite),                                        &
              intent(inout) ::   sh2o

       integer,  dimension( its:ite),                                          &
              intent(in   ) ::   ivgtyp, isltyp

       real,  dimension( kts:kte, its:ite),                                    &
              intent(in   ) ::   alt, u_phy, v_phy, rho_phy, dz8w

       real,  dimension( kts:kte+1, its:ite),                                  &
              intent(in   ) ::   p8w

       real, dimension( 1:num_soil_layers, its:ite),                           &
              intent(inout) ::   smois
       real,  dimension( 1:3, its:ite),                                        &
          intent(in   ) ::    erod
       real,  dimension( its:ite),                                             &
          intent(in   ) ::    sstemp,snowc,vegfra
       real, dimension( 1:num_chem, kts:kte, its:ite),                         &
              intent(inout) ::   tend_emis
       real, dimension( its:ite),                                              &
              intent(inout) ::   e_dust, e_seas
       real, dimension( 1:16, its:ite),                                        &
              intent(in) ::   soilctop
       real, dimension( 1:24, its:ite),                                        &
              intent(in) ::   landusef
       character(len=StrKIND), pointer :: config_source_opt
       logical :: source_is_ssm
! do seasalt emissions
if( lptr_cl_aer(1,1,1) >= param_first_scalar .or. lptr_na_aer(1,1,1) >= param_first_scalar ) then
    if( seasalt_emiss_active == 3 ) then
          call mosaic_seasalt_emiss2(                               &
          dtstep, u10, v10, u_phy,v_phy, alt, rho_phy, dz8w,        &
          xland, dx,g,tend_emis,sstemp,e_seas,emis_seas_factor,     &
          do_seas_sst_adjust,                                       &
          its,ite, kts,kte, seasalt_emiss_active                    )
    else if( seasalt_emiss_active > 0) then
        call mosaic_seasalt_emiss(                                  &
        dtstep, u10, v10, alt, dz8w, xland,                         &
        tend_emis, e_seas,emis_seas_factor,                         &
        its,ite, kts,kte, seasalt_emiss_active                      )
    else if( seasalt_emiss_active > 4) then
        write(iamas_err_message,'(A,I3)') 'illegal value for config_seas_opt, must chosen from 1~4:', seasalt_emiss_active
        call chemistry_error_fatal(iamas_err_message)
    endif
else
    call chemistry_message("Don't do seasalt emission.")
endif

if( lptr_dust_aer(1,1,1) >= param_first_scalar ) then
    call iamas_pool_get_config(configs, 'config_source_opt', config_source_opt)
    source_is_ssm = (trim(config_source_opt) == 'ssm')

    if (dust_emiss_active == 'mosaic') then
        call mosaic_dust_emiss (                                            &
                 slai = slai,      ust = ust,        smois = smois        , &
                 ivgtyp = ivgtyp,  isltyp = isltyp                        , &
                 dtstep = dtstep,   u10 = u10,      v10 = v10             , &
                 alt = alt,        dz8w = dz8w,     xland = xland         , &
                 tend_emis = tend_emis, e_dust = e_dust                   , &
                 its = its, ite = ite, kts = kts, kte = kte                 &
                )
    elseif (dust_emiss_active == 'gocart') then

        if (.not. has_done_gocart_init) then
            call chemistry_debug_message("Reading from file: SOILPARM.TBL")
            call gocart_dust_init_foriAMAS(configs)
            has_done_gocart_init = .true.
        endif
        call mosaic_dust_gocartemis (&
                 ktau = 0 ,dt = dtstep ,num_soil_layers = num_soils       , &
                 alt  = alt,        u_phy = u_phy, v_phy = v_phy          , &
                 rho_phy = rho_phy, dz8w = dz8w, smois = smois            , &
                 u10 = u10, v10 = v10, p8w = p8w, erod = erod             , &
                 source_is_ssm = source_is_ssm                            , &
                 ivgtyp = ivgtyp, isltyp = isltyp, xland = xland          , &
                 xlat = xlat, xlong = xlong, dx = dx                      , &
                 g = g, emis_dust_factor = emis_dust_factor               , &
                 tend_emis = tend_emis, e_dust = e_dust                   , &
                 its = its, ite = ite, kts = kts, kte = kte                 &
                 )
    elseif (dust_emiss_active == 'Kok') then
        call mosaic_dust_Kokemis (ktau = 0,  dt = dtstep                          , &
            num_soil_layers = num_soils, alt = alt, u_phy = u_phy                 , &
            v_phy = v_phy, rho_phy = rho_phy, dz8w = dz8w, smois = smois          , &
            ust = ust  ,u10 = u10,v10 = v10,p8w = p8w                             , &
            ivgtyp = ivgtyp,isltyp = isltyp,xland = xland                         , &
            xlat = xlat,xlong = xlong,dx = dx,g = g                               , &
            e_dust = e_dust,emis_dust_factor = emis_dust_factor                   , &
            num_land = 24,landf = landusef,num_soil = 16,soilf = soilctop         , &
            tend_emis = tend_emis,sh2o = sh2o ,snowc = snowc,xlai = xlai          , &
            its = its, ite = ite, kts = kts, kte = kte                              &
            )
    elseif (dust_emiss_active == 'shao_ver1_schem1') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 1, shao_v = 1, num_soil = 16                                    , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    elseif (dust_emiss_active == 'shao_ver1_schem2') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 2,shao_v = 1,num_soil = 16                                      , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    elseif (dust_emiss_active == 'shao_ver1_schem3') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 3,shao_v = 1,num_soil = 16                                      , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    elseif (dust_emiss_active == 'shao_ver2_schem1') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 1,shao_v = 2,num_soil = 16                                      , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    elseif (dust_emiss_active == 'shao_ver2_schem2') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 2,shao_v = 2,num_soil = 16                                      , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    elseif (dust_emiss_active == 'shao_ver2_schem3') then
        call mosaic_dust_shaoemis(ktau = 0,dt = dtstep,num_soil_layers = num_soils , &
            rho_phy = rho_phy,dz8w = dz8w,smois = smois,ust = ust                  , &
            isltyp = isltyp,vegfra = vegfra,g = g                                  , &
            imod = 3,shao_v = 2,num_soil = 16                                      , &
            soil_top_cat = soilctop,erod = erod,alt = alt                          , &
            source_is_ssm = source_is_ssm                                          , &
            tend_emis = tend_emis,e_dust = e_dust,xlat = xlat,xlong = xlong        , &
            its = its, ite = ite, kts = kts, kte = kte                               &
            )
    end if 
else
    call chemistry_message("Don't do dust emission.")
endif
end subroutine mosaic_addemiss_nature 


! zhxia add this subroutine for anthropogenic emission update. 04/2025
subroutine mosaic_addemiss_anthropogenic(configs, tend_emis, dz8w, alt              , &
            e_nh3    ,e_co    ,e_no      ,e_no2       ,e_so2                        , &
            e_bc     ,e_oc    ,e_pm25    ,e_pm10      ,e_so4                        , &
            e_no3    ,e_c2h6  ,e_c3h8    ,e_c2h2      ,e_alk3                       , &
            e_alk4   ,e_alk5  ,e_ethene  ,e_c3h6      ,e_ole1                       , &
            e_ole2   ,e_aro1  ,e_aro2    ,e_hcho      ,e_ccho                       , &
            e_rcho   ,e_acet  ,e_mek     ,e_isoprene  ,e_terp                       , &
            e_sesq   ,e_phen  ,e_cres    ,e_meoh      ,e_gly                        , &
            e_mgly   ,e_bacl  ,e_isoprod ,e_methacro  ,e_mvk                        , &
            e_prod2  ,e_ch4   ,e_bald    ,e_hcooh     ,e_cco_oh                     , &
            e_rco_oh ,e_co2   ,its,ite   ,kts,kte                                     )

    implicit none
    type(iamas_pool_type),intent(in):: configs

    real, dimension( 1:num_chem, kts:kte, its:ite)                                 , &
              intent(inout) ::   tend_emis

    real,  dimension( kts:kte, its:ite)                                            , &
           intent(in   ) ::   alt, dz8w

    real,  dimension(its:ite)                                                      , &
              intent(in   ) :: e_nh3    ,e_co    ,e_no      ,e_no2       ,e_so2    , &
                               e_bc     ,e_oc    ,e_pm25    ,e_pm10      ,e_so4    , &
                               e_no3    ,e_c2h6  ,e_c3h8    ,e_c2h2      ,e_alk3   , &
                               e_alk4   ,e_alk5  ,e_ethene  ,e_c3h6      ,e_ole1   , &
                               e_ole2   ,e_aro1  ,e_aro2    ,e_hcho      ,e_ccho   , &
                               e_rcho   ,e_acet  ,e_mek     ,e_isoprene  ,e_terp   , &
                               e_sesq   ,e_phen  ,e_cres    ,e_meoh      ,e_gly    , &
                               e_mgly   ,e_bacl  ,e_isoprod ,e_methacro  ,e_mvk    , &
                               e_prod2  ,e_ch4   ,e_bald    ,e_hcooh     ,e_cco_oh , &
                               e_rco_oh ,e_co2

    INTEGER,      INTENT(IN   ) :: its,ite, kts,kte

    !local vars
    integer i,k,n,l

    !----------------------------------------------------------------------------------
    ! Add vertical weights to adjust emission height. zhxia@2026 
    !----------------------------------------------------------------------------------
    ! NOx vertical distribution weights
    ! NOx source composition in North China (estimated from MEIC 2019):：
    ! Road traffic ~30% + Residential heating ~20% → all at k=1 (surface)
    ! Industry ~35% → distributed across k=1~4
    ! Power plants ~15% → concentrated in k=3~6
    ! Combined weighted proportions for each layer:
    !   k=1(~24m):  50%   k=2(~77m):  15%   k=3(~142m): 20%
    !   k=4(~221m):  8%   k=5(~314m):  5%   k=6(~421m):  2%
    real, save :: w_nox(6)
    data w_nox / 0.50, 0.15, 0.20, 0.08, 0.05, 0.02 /
 
    ! SO2 vertical distribution weights
    ! SO2 has a higher proportion from power plants and industry, and a lower proportion from surface sources compared to NOx:
    !   k=1(~24m):  40%   k=2(~77m):  15%   k=3(~142m): 22%
    !   k=4(~221m): 13%   k=5(~314m):  7%   k=6(~421m):  3%
    real, save :: w_so2(6)
    data w_so2 / 0.40, 0.15, 0.22, 0.13, 0.07, 0.03 /
    !----------------------------------------------------------------------------------

    ! now the anthropogenic and biomass burning emissiosn use same distribution, it can be improve later --czhao
    ! not 100%, because most mass falling down to the size less than the lowest boundary
    ! fraction of sorgam i/aitken mode emissions that go to each of the mosaic 8 "standard" sections
    ! Default
    real, save :: fr8b_aem_mosaic_i(8) = (/0.500,0.500,0.000,0.000,0.000,0.000,0.000,0.000/)
    ! fraction of sorgam j/accum mode emissions that go to each of the mosaic 8 "standard" sections\
    ! Default
    real, save :: fr8b_aem_mosaic_j(8) = (/0.024,0.147,0.350,0.332,0.125,0.019,0.001,0.000/)

    ! fraction of mosaic fine (< 2.5 um) emissions that go to each of the mosaic 8 "standard" sections
    !real, save :: fr8b_aem_mosaic_f(8) = (/ 0.060, 0.045, 0.245, 0.400, 0.100, 0.150, 0., 0./) !10-may-2005
    real :: fr8b_aem_mosaic_f(8)
    
    ! fraction of mosaic coarse (> 2.5 um) emissions that go to each of the mosaic 8 "standard" sections
    real, save :: fr8b_aem_mosaic_c(8) = (/ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.300, 0.700 /) ! 10-may-2005

    real :: fr_aem_mosaic_f(8)
    real :: fr_aem_mosaic_c(8)
    real :: ant_emass_ifrac,ant_emass_jfrac
    real :: aem_oc,aem_bc,aem_no3,aem_so4,aem_oin,aem_num
    integer :: itype,iphase

    !emis(gas): mol km-2 hour-1, emis(aerosol): ug m-2 s-1.
    !alt: m3 kg-1, M_air = 0.0289634 kg mol-1, dz8w: m.
    !tend(gas): ppm s-1, tend(aerosol): ug kg-1 s-1.
    !chem(gas): ppm, chem(aerosol): ug kg-1.

    ! gas (if we use the anthropogenic emission and biomass burning emission at the same time, we need to add them together. Please note that we need to set initial values of tend_emis to 0 at each timestep. @zhxia)
    tend_emis(p_co2,1,its:ite)      = tend_emis(p_co2,1,its:ite)      + e_co2(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_nh3,1,its:ite)      = tend_emis(p_nh3,1,its:ite)      + e_nh3(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_co,1,its:ite)       = tend_emis(p_co,1,its:ite)       + e_co(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_c2h6,1,its:ite)     = tend_emis(p_c2h6,1,its:ite)     + e_c2h6(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_c3h8,1,its:ite)     = tend_emis(p_c3h8,1,its:ite)     + e_c3h8(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_c2h2,1,its:ite)     = tend_emis(p_c2h2,1,its:ite)     + e_c2h2(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_alk3,1,its:ite)     = tend_emis(p_alk3,1,its:ite)     + e_alk3(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_alk4,1,its:ite)     = tend_emis(p_alk4,1,its:ite)     + e_alk4(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_alk5,1,its:ite)     = tend_emis(p_alk5,1,its:ite)     + e_alk5(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_ethene,1,its:ite)   = tend_emis(p_ethene,1,its:ite)   + e_ethene(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_c3h6,1,its:ite)     = tend_emis(p_c3h6,1,its:ite)     + e_c3h6(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_ole1,1,its:ite)     = tend_emis(p_ole1,1,its:ite)     + e_ole1(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_ole2,1,its:ite)     = tend_emis(p_ole2,1,its:ite)     + e_ole2(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_aro1,1,its:ite)     = tend_emis(p_aro1,1,its:ite)     + e_aro1(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_aro2,1,its:ite)     = tend_emis(p_aro2,1,its:ite)     + e_aro2(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_hcho,1,its:ite)     = tend_emis(p_hcho,1,its:ite)     + e_hcho(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_ccho,1,its:ite)     = tend_emis(p_ccho,1,its:ite)     + e_ccho(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_rcho,1,its:ite)     = tend_emis(p_rcho,1,its:ite)     + e_rcho(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_acet,1,its:ite)     = tend_emis(p_acet,1,its:ite)     + e_acet(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_mek,1,its:ite)      = tend_emis(p_mek,1,its:ite)      + e_mek(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_isoprene,1,its:ite) = tend_emis(p_isoprene,1,its:ite) + e_isoprene(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_terp,1,its:ite)     = tend_emis(p_terp,1,its:ite)     + e_terp(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_sesq,1,its:ite)     = tend_emis(p_sesq,1,its:ite)     + e_sesq(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_phen,1,its:ite)     = tend_emis(p_phen,1,its:ite)     + e_phen(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_cres,1,its:ite)     = tend_emis(p_cres,1,its:ite)     + e_cres(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_meoh,1,its:ite)     = tend_emis(p_meoh,1,its:ite)     + e_meoh(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_gly,1,its:ite)      = tend_emis(p_gly,1,its:ite)      + e_gly(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_mgly,1,its:ite)     = tend_emis(p_mgly,1,its:ite)     + e_mgly(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_bacl,1,its:ite)     = tend_emis(p_bacl,1,its:ite)     + e_bacl(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_isoprod,1,its:ite)  = tend_emis(p_isoprod,1,its:ite)  + e_isoprod(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_methacro,1,its:ite) = tend_emis(p_methacro,1,its:ite) + e_methacro(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_mvk,1,its:ite)      = tend_emis(p_mvk,1,its:ite)      + e_mvk(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_prod2,1,its:ite)    = tend_emis(p_prod2,1,its:ite)    + e_prod2(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_ch4,1,its:ite)      = tend_emis(p_ch4,1,its:ite)      + e_ch4(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_bald,1,its:ite)     = tend_emis(p_bald,1,its:ite)     + e_bald(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_hcooh,1,its:ite)    = tend_emis(p_hcooh,1,its:ite)    + e_hcooh(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_cco_oh,1,its:ite)   = tend_emis(p_cco_oh,1,its:ite)   + e_cco_oh(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600
    tend_emis(p_rco_oh,1,its:ite)   = tend_emis(p_rco_oh,1,its:ite)   + e_rco_oh(its:ite) * alt(1,its:ite) * 0.0289634 / dz8w(1,its:ite) / 3600

    ! ── Vertical distribution of NO / NO2 / SO2 (6-layer weighting) ────────────────
    ! Original code assigns all to k=1, now changed to 6-layer weighting with total conservation (sum of weights = 1.0). zhxia@2026 
    do k = 1, 6
        tend_emis(p_no, k,its:ite)  = tend_emis(p_no, k,its:ite)  + &
            w_nox(k) * e_no(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
        tend_emis(p_no2,k,its:ite)  = tend_emis(p_no2,k,its:ite)  + &
            w_nox(k) * e_no2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
        tend_emis(p_so2,k,its:ite)  = tend_emis(p_so2,k,its:ite)  + &
            w_so2(k) * e_so2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
    end do

    ! aerosol (if we use the anthropogenic emission and biomass burning emission at the same time, we need to add them together. Please note that we need to set initial values of tend_emis to 0 at each timestep. @zhxia)
    fr_aem_mosaic_f(:) = 0.0
    fr_aem_mosaic_c(:) = 0.0
    ant_emass_ifrac=0.20
    ant_emass_jfrac=0.80
    fr8b_aem_mosaic_f = fr8b_aem_mosaic_i*ant_emass_ifrac+fr8b_aem_mosaic_j*ant_emass_jfrac    
    ! = (/0.100,0.100,0.000,0.000,0.000,0.000,0.000,0.000/) + (/0.0192,0.1176,0.280,0.2656,0.100,0.0152,0.0008,0.000/) 
    ! = (/0.1192,0.2176,0.280,0.2656,0.100,0.0152,0.0008,0.000/)

    itype = 1
    iphase = ai_phase

    if (nsize_aer(itype) .eq. 8) then
       fr_aem_mosaic_f(:) = fr8b_aem_mosaic_f(:)
       fr_aem_mosaic_c(:) = fr8b_aem_mosaic_c(:)
    else if (nsize_aer(itype) .eq. 4) then
       do n = 1, nsize_aer(itype)
          fr_aem_mosaic_f(n) = fr8b_aem_mosaic_f(2*n-1) + fr8b_aem_mosaic_f(2*n)
          fr_aem_mosaic_c(n) = fr8b_aem_mosaic_c(2*n-1) + fr8b_aem_mosaic_c(2*n)
       end do
    end if 

    do n = 1, nsize_aer(itype)
       do i = its, ite
          aem_so4 = fr_aem_mosaic_f(n)*ant_emass_jfrac*e_so4(i)*alt(1,i)/dz8w(1,i)
          aem_no3 = fr_aem_mosaic_f(n)*ant_emass_jfrac*e_no3(i)*alt(1,i)/dz8w(1,i)
          aem_oc  = fr_aem_mosaic_f(n)*ant_emass_jfrac*e_oc(i)*alt(1,i)/dz8w(1,i)
          aem_bc  = fr_aem_mosaic_f(n)*ant_emass_jfrac*e_bc(i)*alt(1,i)/dz8w(1,i)
          aem_oin = fr_aem_mosaic_f(n)*ant_emass_jfrac*e_pm25(i)*alt(1,i)/dz8w(1,i) + fr_aem_mosaic_c(n)*e_pm10(i)*alt(1,i)/dz8w(1,i)
 
          l = lptr_so4_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_so4
 
          l = lptr_no3_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_no3

          l = lptr_oc_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_oc

          l = lptr_bc_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_bc

          l = lptr_oin_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_oin

          ! compute number emissions
          ! first sum the mass-emissions/density
          aem_num = (aem_so4/dens_so4_aer) + (aem_no3/dens_no3_aer) + (aem_oc/dens_oc_aer ) + (aem_bc /dens_bc_aer ) + (aem_oin/dens_oin_aer)
          ! then multiply by 1.0e-6 to convert ug to g and divide by particle volume at center of section (cm3)
          aem_num = aem_num*1.0e-6/volumcen_sect(n,itype)
          l = numptr_aer(n,itype,iphase)
          tend_emis(l,1,i)  = tend_emis(l,1,i)  + aem_num

       enddo
    enddo

end subroutine mosaic_addemiss_anthropogenic

! zhxia add this subroutine for biomass burning emission update. 04/2025
subroutine mosaic_addemiss_bioms(configs, tend_emis, dz8w, alt                          , &
                                 ebu_bc    ,ebu_oc     ,ebu_pm25  ,ebu_pm10  ,ebu_ch4   , &
                                 ebu_co2   ,ebu_co     ,ebu_no    ,ebu_no2   ,ebu_so2   , &
                                 ebu_acet  ,ebu_alk3   ,ebu_alk4  ,ebu_alk5  ,ebu_aro1  , &
                                 ebu_aro2  ,ebu_bald   ,ebu_ccho  ,ebu_c2o2h ,ebu_ethen , &
                                 ebu_hcho  ,ebu_hcooh  ,ebu_hono  ,ebu_isopr ,ebu_mek   , &
                                 ebu_meoh  ,ebu_metha  ,ebu_mgly  ,ebu_mvk   ,ebu_ole1  , &
                                 ebu_ole2  ,ebu_phen   ,ebu_prod2 ,ebu_rcho             , &
                                 its,ite  ,kts,kte                                   )

    implicit none
    type(iamas_pool_type),intent(in):: configs

    real, dimension( 1:num_chem, kts:kte, its:ite)                                    , &
              intent(inout) ::   tend_emis

    real,  dimension( kts:kte, its:ite)                                               , &
           intent(in   ) ::   alt,dz8w

    real,  dimension(its:ite)                                                         , &
              intent(in   ) :: ebu_bc    ,ebu_oc     ,ebu_pm25  ,ebu_pm10  ,ebu_ch4   , &
                               ebu_co2   ,ebu_co     ,ebu_no    ,ebu_no2   ,ebu_so2   , &
                               ebu_acet  ,ebu_alk3   ,ebu_alk4  ,ebu_alk5  ,ebu_aro1  , &
                               ebu_aro2  ,ebu_bald   ,ebu_ccho  ,ebu_c2o2h ,ebu_ethen , &
                               ebu_hcho  ,ebu_hcooh  ,ebu_hono  ,ebu_isopr ,ebu_mek   , &
                               ebu_meoh  ,ebu_metha  ,ebu_mgly  ,ebu_mvk   ,ebu_ole1  , &
                               ebu_ole2  ,ebu_phen   ,ebu_prod2 ,ebu_rcho
 
    INTEGER,      INTENT(IN   ) :: its,ite, kts,kte

    !local vars
    integer i,k,n,l

    ! now the anthropogenic and biomass burning emissiosn use same distribution, it can be improve later --czhao
    ! not 100%, because most mass falling down to the size less than the lowest boundary
    ! fraction of sorgam i/aitken mode emissions that go to each of the mosaic 8 "standard" sections
    ! Default
    real, save :: fr8b_aem_mosaic_i(8) = (/0.500,0.500,0.000,0.000,0.000,0.000,0.000,0.000/)
    ! fraction of sorgam j/accum mode emissions that go to each of the mosaic 8 "standard" sections\
    ! Default
    real, save :: fr8b_aem_mosaic_j(8) = (/0.024,0.147,0.350,0.332,0.125,0.019,0.001,0.000/)

    ! fraction of mosaic fine (< 2.5 um) emissions that go to each of the mosaic 8 "standard" sections
    !real, save :: fr8b_aem_mosaic_f(8) = (/ 0.060, 0.045, 0.245, 0.400, 0.100, 0.150, 0., 0./) !10-may-2005
    real :: fr8b_aem_mosaic_f(8)
    
    ! fraction of mosaic coarse (> 2.5 um) emissions that go to each of the mosaic 8 "standard" sections
    real, save :: fr8b_aem_mosaic_c(8) = (/ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.300, 0.700 /) ! 10-may-2005

    real :: fr_aem_mosaic_f(8)
    real :: fr_aem_mosaic_c(8)
    real :: ant_emass_ifrac,ant_emass_jfrac
    real :: aem_oc,aem_bc,aem_oin,aem_num
    integer :: itype,iphase

    real, parameter :: w_bio(10) = (/ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 /)

    !emis(gas): mol km-2 hour-1, emis(aerosol): ug m-2 s-1.
    !alt: m3 kg-1, M_air = 0.0289634 kg mol-1, dz8w: m.
    !tend(gas): ppm s-1, tend(aerosol): ug kg-1 s-1.
    !chem(gas): ppm, chem(aerosol): ug kg-1.

    ! gas (if we use the anthropogenic emission and biomass burning emission at the same time, we need to add them together. Please note that we need to set initial values of tend_emis to 0 at each timestep. @zhxia)
    do k = 1, 10
       tend_emis(p_ch4,k,its:ite)      = tend_emis(p_ch4,k,its:ite) + w_bio(k) * ebu_ch4(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_co2,k,its:ite)      = tend_emis(p_co2,k,its:ite) + w_bio(k) * ebu_co2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_co,k,its:ite)       = tend_emis(p_co,k,its:ite)  + w_bio(k) * ebu_co(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_no,k,its:ite)       = tend_emis(p_no,k,its:ite)  + w_bio(k) * ebu_no(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_no2,k,its:ite)      = tend_emis(p_no2,k,its:ite) + w_bio(k) * ebu_no2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_so2,k,its:ite)      = tend_emis(p_so2,k,its:ite) + w_bio(k) * ebu_so2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_acet,k,its:ite)     = tend_emis(p_acet,k,its:ite) + w_bio(k) * ebu_acet(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_alk3,k,its:ite)     = tend_emis(p_alk3,k,its:ite) + w_bio(k) * ebu_alk3(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_alk4,k,its:ite)     = tend_emis(p_alk4,k,its:ite) + w_bio(k) * ebu_alk4(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600 
       tend_emis(p_alk5,k,its:ite)     = tend_emis(p_alk5,k,its:ite) + w_bio(k) * ebu_alk5(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_aro1,k,its:ite)     = tend_emis(p_aro1,k,its:ite) + w_bio(k) * ebu_aro1(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_aro2,k,its:ite)     = tend_emis(p_aro2,k,its:ite) + w_bio(k) * ebu_aro2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_bald,k,its:ite)     = tend_emis(p_bald,k,its:ite) + w_bio(k) * ebu_bald(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_ccho,k,its:ite)     = tend_emis(p_ccho,k,its:ite) + w_bio(k) * ebu_ccho(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_cco_oh,k,its:ite)   = tend_emis(p_cco_oh,k,its:ite)  + w_bio(k) * ebu_c2o2h(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_ethene,k,its:ite)   = tend_emis(p_ethene,k,its:ite)  + w_bio(k) * ebu_ethen(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_hcho,k,its:ite)     = tend_emis(p_hcho,k,its:ite) + w_bio(k) * ebu_hcho(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_hcooh,k,its:ite)    = tend_emis(p_hcooh,k,its:ite) + w_bio(k) * ebu_hcooh(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_hono,k,its:ite)     = tend_emis(p_hono,k,its:ite) + w_bio(k) * ebu_hono(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_isoprene,k,its:ite) = tend_emis(p_isoprene,k,its:ite) + w_bio(k) * ebu_isopr(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_mek,k,its:ite)      = tend_emis(p_mek,k,its:ite)  + w_bio(k) * ebu_mek(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_meoh,k,its:ite)     = tend_emis(p_meoh,k,its:ite)  + w_bio(k) * ebu_meoh(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_methacro,k,its:ite) = tend_emis(p_methacro,k,its:ite) + w_bio(k) * ebu_metha(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_mgly,k,its:ite)     = tend_emis(p_mgly,k,its:ite) + w_bio(k) * ebu_mgly(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_mvk,k,its:ite)      = tend_emis(p_mvk,k,its:ite) + w_bio(k) * ebu_mvk(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_ole1,k,its:ite)     = tend_emis(p_ole1,k,its:ite) + w_bio(k) * ebu_ole1(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_ole2,k,its:ite)     = tend_emis(p_ole2,k,its:ite)  + w_bio(k) * ebu_ole2(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_phen,k,its:ite)     = tend_emis(p_phen,k,its:ite)  + w_bio(k) * ebu_phen(its:ite)  * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_prod2,k,its:ite)    = tend_emis(p_prod2,k,its:ite) + w_bio(k) * ebu_prod2(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
       tend_emis(p_rcho,k,its:ite)     = tend_emis(p_rcho,k,its:ite) + w_bio(k) * ebu_rcho(its:ite) * alt(k,its:ite) * 0.0289634 / dz8w(k,its:ite) / 3600
    enddo

    ! aerosol (if we use the anthropogenic emission and biomass burning emission at the same time, we need to add them together. Please note that we need to set initial values of tend_emis to 0 at each timestep. @zhxia)
    fr_aem_mosaic_f(:) = 0.0
    fr_aem_mosaic_c(:) = 0.0
    ant_emass_ifrac=0.20
    ant_emass_jfrac=0.80
    fr8b_aem_mosaic_f = fr8b_aem_mosaic_i*ant_emass_ifrac+fr8b_aem_mosaic_j*ant_emass_jfrac    
    ! = (/0.100,0.100,0.000,0.000,0.000,0.000,0.000,0.000/) + (/0.0192,0.1176,0.280,0.2656,0.100,0.0152,0.0008,0.000/) 
    ! = (/0.1192,0.2176,0.280,0.2656,0.100,0.0152,0.0008,0.000/)

    itype = 1
    iphase = ai_phase

    if (nsize_aer(itype) .eq. 8) then
       fr_aem_mosaic_f(:) = fr8b_aem_mosaic_f(:)
       fr_aem_mosaic_c(:) = fr8b_aem_mosaic_c(:)
    else if (nsize_aer(itype) .eq. 4) then
       do n = 1, nsize_aer(itype)
          fr_aem_mosaic_f(n) = fr8b_aem_mosaic_f(2*n-1) + fr8b_aem_mosaic_f(2*n)
          fr_aem_mosaic_c(n) = fr8b_aem_mosaic_c(2*n-1) + fr8b_aem_mosaic_c(2*n)
       end do
    end if 

    do n = 1, nsize_aer(itype)
       do i = its, ite
       !do 1810 k =   1, min(configs%kemit,kte)
          do k = 1,10
             aem_oc  = w_bio(k) * fr_aem_mosaic_f(n)*ant_emass_jfrac*ebu_oc(i)*alt(k,i)/dz8w(k,i)
             aem_bc  = w_bio(k) * fr_aem_mosaic_f(n)*ant_emass_jfrac*ebu_bc(i)*alt(k,i)/dz8w(k,i)
             aem_oin = w_bio(k) * fr_aem_mosaic_f(n)*ant_emass_jfrac*ebu_pm25(i)*alt(k,i)/dz8w(k,i) + w_bio(k) * fr_aem_mosaic_c(n)*ebu_pm10(i)*alt(k,i)/dz8w(k,i)

             l = lptr_oc_aer(n,itype,iphase)
             tend_emis(l,k,i)  = tend_emis(l,k,i)  + aem_oc

             l = lptr_bc_aer(n,itype,iphase)
             tend_emis(l,k,i)  = tend_emis(l,k,i)  + aem_bc

             l = lptr_oin_aer(n,itype,iphase)
             tend_emis(l,k,i)  = tend_emis(l,k,i)  + aem_oin

             ! compute number emissions
             ! first sum the mass-emissions/density
             aem_num = (aem_oc/dens_oc_aer ) + (aem_bc /dens_bc_aer ) + (aem_oin/dens_oin_aer)
             ! then multiply by 1.0e-6 to convert ug to g and divide by particle
             ! volume at center of section (cm3)
             aem_num = aem_num*1.0e-6/volumcen_sect(n,itype)
             l = numptr_aer(n,itype,iphase)
             tend_emis(l,k,i)  = tend_emis(l,k,i)  + aem_num
          enddo
       enddo
    enddo

end subroutine mosaic_addemiss_bioms

!----------------------------------------------------------------------

subroutine mosaic_dust_emiss(        slai,ust, smois, ivgtyp, isltyp  , &
            dtstep, u10, v10, alt, dz8w, xland                        , &
            tend_emis, e_dust                                         , &
            its,ite, kts,kte                                            )

!
! adds dust emissions for mosaic aerosol species (i.e. emission tendencies
! over time dtstep are applied to the aerosol mixing ratios)
!
! This is a simple dust scheme based on Shaw et al. (2008) to appear in
! Atmospheric Environment, recoded by Jerome Fast
!
! NOTE: 
! 1) This version only works with the 8-bin version of MOSAIC.
! 2) Dust added to MOSAIC's other inorganic specie, OIN.  If Ca and CO3 are 
!    activated in the Registry, a small fraction also added to Ca and CO3.
! 3) The main departure from Shaw et al., is now alphamask is computed since
!    the land-use categories in that paper and in WRF differ.  WRF currently 
!    does not have that many land-use categories and adhoc assumptions had to
!    be made. This version was tested for Mexico in the dry season.  The main
!    land-use categories in WRF that are likely dust sources are grass, shrub,
!    and savannna (that WRF has in the desert regions of NW Mexico).  Having
!    dust emitted from these types for other locations and other times of the
!    year is not likely to be valid.
! 4) An upper bound on ustar was placed because the surface parameterizations
!    in WRF can produce unrealistically high values that lead to very high
!    dust emission rates.
! 5) Other departures' from Shaw et al. noted below, but are probably not as
!    important as 2) and 3).
! 6) Now the dust added into dust species

!  USE module_configure, only:  grid_config_rec_type
   IMPLICIT NONE


   INTEGER,      INTENT(IN   ) :: its,ite, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( its:ite),                                             &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( its:ite),                                          &
          INTENT(IN   ) ::   ivgtyp, isltyp

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( kts:kte, its:ite),                                    &
          INTENT(IN   ) ::   alt, dz8w
   REAL, DIMENSION( 1:num_soils, its:ite),                                 &
          INTENT(INOUT) ::   smois
   REAL, DIMENSION( 1:num_chem, kts:kte, its:ite),                         &
          INTENT(INOUT) ::   tend_emis
   REAL, DIMENSION( its:ite),                                              &
          INTENT(INOUT) ::   e_dust
         
! local variables
        integer i, k, l, l_oin, l_ca, l_co3, n, ii, l_dust
        integer iphase, itype, izob
        integer p1st

        real dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real factaa, factbb, fracoin, fracca, fracco3, fractot
        real ustart, ustar1, ustart0
        real alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust, mass1part
        real :: c_const

        p1st = param_first_scalar
!
! from Nickovic et al., JGR, 2001 and Shaw et al. 2007
! beta: fraction of clay, small silt, large silt, and sand correcsponding to Zobler class (7)
! beta (1,*) for 0.5-1 um
! beta (2,*) for 1-10 um
! beta (3,*) for 10-25 um
! beta (4,*) for 25-50 um
!
        beta(1,1)=0.12
        beta(2,1)=0.04
        beta(3,1)=0.04
        beta(4,1)=0.80
        beta(1,2)=0.34
        beta(2,2)=0.28
        beta(3,2)=0.28
        beta(4,2)=0.10
        beta(1,3)=0.45
        beta(2,3)=0.15
        beta(3,3)=0.15
        beta(4,3)=0.25
        beta(1,4)=0.12
        beta(2,4)=0.09
        beta(3,4)=0.09
        beta(4,4)=0.70
        beta(1,5)=0.40
        beta(2,5)=0.05
        beta(3,5)=0.05
        beta(4,5)=0.50
        beta(1,6)=0.34
        beta(2,6)=0.18
        beta(3,6)=0.18
        beta(4,6)=0.30
        beta(1,7)=0.22
        beta(2,7)=0.09
        beta(3,7)=0.09
        beta(4,7)=0.60
        gamma(1)=0.08
        gamma(2)=1.00
        gamma(3)=1.00
        gamma(4)=0.12

        itype=1
  if (nsize_aer(itype) .eq. 8) then
!jwfeng:04/2020
!Modified to 10 bins, up to 40um
          sz(1)=0
          sz(2)=1.78751e-06
          sz(3)=0.000273786
          sz(4)=0.00847978
          sz(5)=0.056055
          sz(6)=0.0951896
          sz(7)=0.17
          sz(8)=0.67
  else if (nsize_aer(itype) .eq. 4) then
          !BSINGH(12/11/2013): Based on the suggestions by Manish Shrivastva, sz variable is modified below.
          !Original values are commented out
          !sz(1)=0.0
          !sz(2)=0.01
          !sz(3)=0.07
          !sz(4)=0.92
          sz(1)=1.78751e-06
          sz(2)=0.00875357
          sz(3)=0.1512446
          sz(4)=0.84
          sz(5)=0.0
          sz(6)=0.0
          sz(7)=0.0
          sz(8)=0.0
        endif

!   for now just do itype=1
        itype = 1
        iphase = ai_phase

!   loop over i,j and apply dust emissions
        k = kts
        do 1820 i = its, ite

    if( xland(i) > 1.5 ) cycle

! compute wind speed anyway, even though ustar is used below

        dumlandfrac = 1.
        dumspd10=(u10(i)*u10(i) + v10(i)*v10(i))**(0.5)
        if(dumspd10 >= 5.0) then
           dumspd10 = dumlandfrac* &
         ( dumspd10*dumspd10*(dumspd10-5.0))
         else
            dumspd10=0.
         endif

! part1 - compute vegetation mask
!
! * f8, f50, f51, f52 refer to vegetation classes from the Olsen categories
!   for desert, sand desert, grass semi-desert, and shrub semi-desert
! * in WRF, vegetation types 7, 8 and 10 are grassland, shrubland, and savanna
!   that are dominate types in Mexico and probably have some erodable surface
!   during the dry season
! * currently modified these values so that only a small fraction of cell
!   area is erodable
! * these values are highly tuneable!

         alphamask=0.001
         if (ivgtyp(i) .eq. 7) then
           f8=0.005
           f50=0.00
           f51=0.10
           f51=0.066
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i) .eq. 8) then
           f8=0.010
           f50=0.00
           f51=0.00
           f52=0.15
           f52=0.10
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i) .eq. 10) then
           f8=0.00
           f50=0.00
           f51=0.01
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif

! part2 - zobler
! 
! * in Shaw's paper, dust is computed for 4 size ranges:
!   0.5-1 um 
!    1-10 um  
!   10-25 um  
!   25-50 um
! * Shaw's paper also accounts for sub-grid variability in soil
!   texture, but here we just assume the same soil texture for each
!   grid cell
! * since MOSAIC is currently has a maximum size range up to 10 um,
!   neglect upper 2 size ranges and lowest size range (assume small)
! * map WRF soil classes arbitrarily to Zolber soil textural classes
! * skip dust computations for WRF soil classes greater than 13, i.e. 
!   do not compute dust over water, bedrock, and other surfaces
! * should be skipping for water surface at this point anyway
!
         izob=0
         if(isltyp(i).eq.1) izob=1
         if(isltyp(i).eq.2) izob=1
         if(isltyp(i).eq.3) izob=4
         if(isltyp(i).eq.4) izob=2
         if(isltyp(i).eq.5) izob=2
         if(isltyp(i).eq.6) izob=2
         if(isltyp(i).eq.7) izob=7
         if(isltyp(i).eq.8) izob=2
         if(isltyp(i).eq.9) izob=6
         if(isltyp(i).eq.10) izob=5
         if(isltyp(i).eq.11) izob=2
         if(isltyp(i).eq.12) izob=3
         if(isltyp(i).ge.13) izob=0
         if(izob.eq.0) goto 1840
!
! part3 - dustprod
!
         do ii=1,4
           delta(ii)=0.0
         enddo
         sumdelta=0.0
         do ii=1,4
           delta(ii)=beta(ii,izob)*gamma(ii)
           if(ii.lt.4) then
             sumdelta=sumdelta+delta(ii)
           endif
         enddo
         do ii=1,4
           delta(ii)=delta(ii)/sumdelta
         enddo

! part4 - wetness
!
! * assume dry for now, have passed in soil moisture to this routine
!   but needs to be included here
! * wetfactor less than 1 would reduce dustflux
! * convert model soil moisture (m3/m3) to gravimetric soil moisture
!   (mass of water / mass of soil in %) assuming a constant density 
!   for soil
         pclay=beta(1,izob)*100.
         wp=0.0014*pclay*pclay+0.17*pclay
         smois_grav=(smois(1,i)/2.6)*100.
         if(smois_grav.gt.wp) then
           wetfactor=sqrt(1.0+1.21*(smois_grav-wp)**0.68)
         else
           wetfactor=1.0
         endif
!        wetfactor=1.0

! part5 - dustflux
! lower bound on ustar = 20 cm/s as in Shaw et al, but set upper
! bound to 100 cm/s
         c_const=1.e-14  ! default

         ustar1=ust(i)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=c_const*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0
! units kg m-2 s-1
         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo
! convert to ug m-2 s-1  
         ftot=ftot*1.0e+09
         e_dust(i) = ftot

!   apportion other inorganics only
         factaa = alt(k,i)/dz8w(k,i)
         factbb = factaa * ftot
         fracoin = 0.97
         fracca = 0.03*0.4
         fracco3 = 0.03*0.6
         fractot = fracoin + fracca + fracco3

         do n = 1, nsize_aer(itype)
            l_dust = lptr_dust_aer(n,itype,iphase)
            l_ca = lptr_ca_aer(n,itype,iphase)
            l_co3 = lptr_co3_aer(n,itype,iphase)
            if (l_dust >= p1st) then
                tend_emis(l_dust,k,i) = factbb * sz(n) * fracoin
            end if
            if (l_ca >= p1st) then
                tend_emis(l_ca,k,i) = factbb * sz(n) * fracca
            end if
            if (l_co3 >= p1st) then
                tend_emis(l_co3,k,i) = factbb * sz(n) * fracco3
            end if
! mass1part is mass of a single particle in ug, density of dust ~2.5 g cm-3
               densdust=2.5
               mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
               l = numptr_aer(n,itype,iphase)
               if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + factbb * sz(n) * fractot / mass1part
         end do


1840    continue

1820    continue
1830    continue

        return

   END subroutine mosaic_dust_emiss

!====================================================================================
!add another dust emission scheme following GOCART mechanism  --czhao  09/17/2009
!====================================================================================
!chenjin add a new way to calculate dust emission
   subroutine mosaic_dust_Kokemis (ktau,dt,num_soil_layers,alt,u_phy,&
         v_phy,rho_phy,dz8w,smois,ust,u10,v10,p8w,                   &
         ivgtyp,isltyp,xland,xlat,xlong,dx,g,                        &
         e_dust,emis_dust_factor,                                    &
         num_land,landf,num_soil,soilf,tend_emis,sh2o,snowc,xlai,    &
         its,ite, kts,kte                                   )

    USE module_data_gocart_dust
    USE module_data_binaer_asect

  IMPLICIT NONE

  INTEGER,      INTENT(IN   ) :: ktau,num_soil_layers,num_land,num_soil,     &
                                 its,ite, kts,kte
  INTEGER,DIMENSION( its:ite),                                               &
         INTENT(IN   ) ::   ivgtyp,isltyp

  REAL, DIMENSION( 1:num_chem, kts:kte, its:ite),                            &
         INTENT(INOUT) ::   tend_emis

         REAL, DIMENSION(1:num_land, its:ite),                               &
         INTENT(IN) ::  landf

!  LOGICAL :: prescribed_baresoil = true

  REAL, DIMENSION(1:num_soil, its:ite),                                      &
         INTENT(IN) ::  soilf

  REAL, DIMENSION(its:ite),                                                  &
         INTENT(IN) ::  snowc,xlai


  REAL,  DIMENSION( its:ite), INTENT(INOUT ) ::  e_dust

  REAL, DIMENSION( 1:num_soil_layers, its:ite),                              &
      INTENT(INOUT) ::               smois,sh2o

  REAL,  DIMENSION( ite:ite),                                                &
          INTENT(IN   ) ::  u10, v10, ust, xland, xlat, xlong    

  REAL,  DIMENSION( kts:kte, its:ite),                                       &
          INTENT(IN   ) ::  alt, u_phy, v_phy, rho_phy, dz8w

  REAL,  DIMENSION( kts:kte+1, its:ite),                                     &
          INTENT(IN   ) ::  p8w

   REAL,  DIMENSION( its:ite),                                               &
          INTENT(IN   ) ::  dx
  REAL, INTENT(IN   ) :: dt,g

  REAL, INTENT(IN ) ::     emis_dust_factor

  !local variables
 ! REAL, PARAMETER  :: Tfactor=5.0*1.e-4 !factor for resolution sensitivity, highly tunable  !no longer needed -czhao 
 ! REAL, PARAMETER  :: Sfactor=1.0       !source erodibility factor, place holder !no longer needed -czhao 
  INTEGER, PARAMETER  :: BARESOIL=19    !bare soil with sparse veg. from USGS 24-class data
  INTEGER, PARAMETER  :: WATERBODY=16    !water body from USGS 24-class data
! LOGICAL,PARAMETER :: prescribe_soil=.false.
  real baresoilf
  real f_v
  real,parameter  :: xlai_thresh=0.3    !threshhold for LAI
 ! REAL alpha !sandblasting mass efficiency (m^-1) !no longer needed -czhao
  REAL f_soil    !grid cell fraction of exposed bare soil 
  real mclay,mclay1 !mass fraction of clay particles 
  real f_clay !fraction of clay soil
 ! real index1  !no longer needed -czhao
 ! real qflux  !total horizontally saltating mass flux (kg m^-1 s^-1) !!variable no longer needed, -czhao
  real term1,term2,term3
  real ustar   ! wind friction speed (m/s) accounting for the Owen effect
  real ustar_t ! threshold wind friction speed for saltation (m/s)
  real ustar_t_std ! [m/s] The soil threshold friction speed at standard air density (1.2250 kg/m3) 
  real, parameter :: ustar_t_std_min = 0.16 ! [m/s] minimum standardized soil threshold friction speed 

  real wind10m,wind10m_t

 ! REAL, PARAMETER  :: CSfactor=2.61    !saltation constant !variable no longer needed
  real :: Cd     ! [dimless] The dust emission coefficient, which depends on the soil's standardized threshold friction speed
  real, parameter :: Cd0 = 4.4e-5      ! [dimless] proportionality constant in calculation of dust emission coefficient 
  real, parameter :: Ca = 2.7          ! [dimless] proportionality constant in scaling of dust emission exponent
  real, parameter :: Ce = 2.0          ! [dimless] proportionality constant scaling exponential dependence of dust emission coefficient on standardized soil threshold friction speed
 !now c_tune is controled by namelist.input
 !real, parameter :: C_tune = 0.25  ! [dimless] global tuning constant for vertical dust flux; set to produce ~same global dust flux in control sim (I_2000) as old parameterization 
  real :: c_tune
  real, parameter :: rho_phy_std = 1.2250  ! [kg/m3] density of air at standard pressure (101325) and temperature (293 K) 

  REAL :: sand_f(num_soil) ! sand fraction of 16-class soil type from Noah land surface model 
  REAL :: silt_f(num_soil) ! silt fraction of 16-class soil type from Noah land surface model 
  REAL :: clay_f(num_soil) ! clay fraction of 16-class soil type from Noah land surface model 
  integer isoil,i,k
  real gwet,gwet_t,bulk_dens !bulk density of soil 
  real,parameter :: liquid_dens = 1.e3  !kg/m3 
  real,parameter :: dosp = 75*10e-6  ! m, optimal saltation particles 
  real fwet !factor dependent on soil moisture
  real,parameter :: fzfactor=1.0 !factor dependent on surface roughness, set as a place holder for now
  real re_t_f,re_t

  real :: rhosoil  !kg/m3

  real totalemis,sz(8),mass1part,densdust
  real factaa, factbb, fracoin, fracca, fracco3, fractot
  integer p1st,l_dust,l_ca,l_co3,l,n,l_tr1r1,l_tr1r2,l_tr1r3,l_tr1r4
  integer iphase, itype
  integer,parameter :: dust_size_option = 2

 !**********************************************************
 !1. Initializing parameters 
 !**********************************************************
   !now use namelist to control global tuning constant for vertical dust flux
   c_tune = emis_dust_factor 
   !  for now just do itype=1
   itype = 1
   iphase = ai_phase
   p1st = param_first_scalar
   sz(:) = 0.0
   rhosoil = dens_dust_aer*1.e3 ! g/cm3 to kg/m3

   !-----------------------------------------
   !1.1 Checking land & soil information 
   !-----------------------------------------
    if (num_land.ne.24) then
      write(0,*) 'incorrect number of land category ',num_land,baresoil
      print*, 'incorrect number of land category ',num_land,baresoil
      stop
    endif

    if (num_soil.ne.16) then
      write(0,*) 'incorrect number of soil category ',num_soil
      print*, 'incorrect number of soil category ',num_soil
      stop
    endif

   !-----------------------------------------
   !1.2 Initialize soil categary percentage 
   !It follows CLM in CESM and can also be found (similarly) in clm_varcon of module_sf_clm.F) 
   !-----------------------------------------
    sand_f = (/92,82,58,17,10,43,58,10,32,52, 6,22,0,0,0,0/)
    silt_f = (/ 5,12,32,70,85,39,15,56,34, 6,47,20,0,0,0,0/)
    clay_f = (/ 3, 6,10,13, 5,18,27,34,34,42,47,58,0,0,0,0/)

   !-----------------------------------------
   !1.3 Initialize dust emission fraction of modes 
   !-----------------------------------------
    !accfrac=emass_jfrac1
    !corfrac=emass_cfrac1
    !accfrac=dust_emass_jfrac
    !corfrac=dust_emass_cfrac

 !**********************************************************
 !2. Dust emission calculation 
 !**********************************************************

  k=kts   !surface only
  do i=its,ite

      !zero the emission array
      e_dust(i)=0.0
      totalemis = 0.0

      !baresoilf fraction
!      if (prescribed_baresoil) then
      baresoilf = landf(baresoil,i)
!      else
!      f_v=xlai(i)/xlai_thresh
!      f_v=max(0.0,f_v)
!      f_v=min(1.0,f_v)
!      baresoilf = (1.-snowc(i))*(1.-landf(i,waterbody,j))*(1.-f_v)
!      endif

   !dust emission only over baresoil land 
    if(xland(i).lt.1.5.and.baresoilf.gt.0)then
   !if(landf(i,waterbody,j).lt.0.5.and.baresoilf.gt.0)then


   !-----------------------------------------------------
   !2.1 grid cell fraction of exposed bare soil (f_soil) 
   !-----------------------------------------------------
      if (smois(1,i).gt.0) then
       f_soil = baresoilf*sh2o(1,i)/smois(1,i)
      else
       f_soil = baresoilf
      endif

   !-----------------------------------------------------
   !2.2 sandblasting mass efficiency (alpha) 
   !-----------------------------------------------------
      !2.2.1 clay percentage 
      f_clay=0.0
      do isoil=1,num_soil
       f_clay=f_clay+clay_f(isoil)*soilf(isoil,i)
      enddo

      !2.2.2 mass fraction of clay particles 
      if (f_clay.gt.20.and.f_clay.le.100) then
       mclay = 20*0.01
      else
       if (f_clay.lt.0) f_clay=0.0
       if (f_clay.gt.20) f_clay=20.0
       mclay = f_clay*0.01
      endif
      !2.2.3 alpha  !no longer needed -czhao 
 !     index1= 13.4*mclay-6.0
 !     alpha = 100.0*exp(index1*alog(10.0))

   !-----------------------------------------------------
   !2.3 threshold wind friction speed for saltation (ustar_t)  
   !-----------------------------------------------------
      !2.3.1 soil moisture content (gwet)
      !chenjin : Add porosity through lsm_noah_init.
      porosity(isltyp(i))=maxsmc(isltyp(i))

      bulk_dens=(1.0-porosity(isltyp(i)))*2.7e3
      gwet=smois(1,i)*liquid_dens/bulk_dens

      !2.3.2 soil moisture dependent factor (fwet) 
      mclay1 = f_clay*0.01
      mclay1 = min(max(mclay1,0.0),1.0)
      if (mclay1==0) then
       gwet_t = 0.0
      else
       gwet_t=(0.17*mclay1+0.14*(mclay1**2))/mclay1
      endif

      if (gwet.le.gwet_t) then
       fwet = 1.0
      else
       fwet = (1+1.21*((100*(gwet-gwet_t))**0.68))**0.5
      endif

      !2.3.3 threshold friction Reynolds factor (re_t_f) 
      re_t = 0.38+1331*((100*dosp)**1.56)
      if (re_t.gt.10) then
      re_t_f = 0.0144*(1-0.0858*exp(-0.0617*(re_t-10)))**2
      else
       if (re_t.ge.0.03) then
        !re_t_f = 0.1291*0.1291/(-1+1.928*re_t)   !from CLM document
        re_t_f = 0.1291*0.1291/(-1+1.928*(re_t**0.0922)) !from CLM code
       else
       goto 999
       endif
      endif

      !2.3.4 ustar_t 
      term1=re_t_f*rhosoil*g*dosp*(1+(6.e-7)/(rhosoil*g*(dosp**2.5)))
      ustar_t=fzfactor*(term1**0.5)*(rho_phy(k,i)**(-0.5))*fwet
      ustar_t_std = ustar_t*sqrt(rho_phy(k,i)/rho_phy_std) ! standardized soil threshold friction speed -czhao      
   !-----------------------------------------------------
   !2.4 horizontally saltating mass flux (kg/m/s) (qflux) 
   !-----------------------------------------------------
      !2.4.1 friction velocity accounting for the Owen effect 
      wind10m=sqrt(u10(i)*u10(i)+v10(i)*v10(i))
      ! we don't trust the u10,v10 values, if model layers are very thin near surface
      if(dz8w(kts,i).le.10.) wind10m=sqrt(u_phy(kts,i)*u_phy(kts,i)+v_phy(kts,i)*v_phy(kts,i))
      wind10m_t=ustar_t*wind10m/ust(i)
      if (wind10m.lt.wind10m_t) then
       ustar = ust(i)
      else
       ustar = ust(i)+0.003*(wind10m-wind10m_t)**2
      endif

 !not needed any longer  -czhao
 !     !2.4.2 friction velocity accounting for the Owen effect 
 !     if (ustar.gt.ustar_t) then
 !      term1 = csfactor*rho_phy(i,kts,j)*(ustar**3)/g
 !      term2 = 1-ustar_t/ustar
 !      term3 = (1+ustar_t/ustar)**2
 !      qflux = term1*term2*term3
 !     else
 !      qflux = 0.0
 !     endif

   !-----------------------------------------------------
   !2.5 vertical mass fluxes of dust (kg m-2 s-1) (totalemis) 
   !-----------------------------------------------------
      if (ustar > ustar_t) then
 !     totalemis = tfactor*sfactor*f_soil*alpha*qflux  !czhao
      Cd = Cd0*exp(-Ce*(ustar_t_std-ustar_t_std_min)/ustar_t_std_min) ! the dust emission coefficient -czhao
      totalemis = Cd*mclay*rho_phy(k,i)*((ustar**2.0-ustar_t**2.0)/ustar_t_std)*(ustar/ustar_t)**(Ca*(ustar_t_std-ustar_t_std_min)/ustar_t_std_min)  ! the vertical dust flux -czhao
      totalemis = totalemis*f_soil*C_tune  ! account for bare soil fraction, frozen soil fraction, and apply global tuning parameter -czhao
      else
      totalemis = 0.0
      endif


      !2.5.1 unit convert from kg/m2/s to ug/m2/s
      totalemis = totalemis*1.e9
      e_dust(i)=totalemis

      !2.5.2 convert emission to mass concentration with size distribution
      !ug/m2/s to ug/kg-dry air (alt in m3/kg)
! chenjin: Using tend_emis instead of directly changing chem
      factaa = alt(k,i)/dz8w(k,i)
      factbb = factaa * totalemis 
      !size distribution
      SELECT CASE (dust_size_option)  !option for dust size distribution 
      CASE(1)  !Zhao et al., 2010  ACP 
      !option 1
      !!the size is following the best match2 of modal distribution
      !!we have ~3.6% within 1 um radius and ~7.4% within 1.25 um radius
      !!versus ~6.0% within 1 um radius following GOCART original and ~9.3% within 1.25 um radius
      !!totally ~82% are within the MOSAIC size range 0.03-10 um diameter
      !!now the dust distribution follows the parameter defined in module_data_sorgam.F 
        sz(1)=1.0e-8
        sz(2)=1.0e-6
        sz(3)=3.0e-4
        sz(4)=3.5e-3
        sz(5)=0.018
        sz(6)=0.070
        sz(7)=0.2595
        sz(8)=0.42
    
      CASE(2) !Kok et al., 2011 PNAS 
      !option 2
      !!based on Kok et al., 2011, we can get another size distribution with lambda=12 and Cv=12.62  -czhao
        sz(1)=7.7e-7
        sz(2)=1.5e-5
        sz(3)=1.9e-4
        sz(4)=1.8e-3
        sz(5)=0.0119
        sz(6)=0.0574
        sz(7)=0.2018
        sz(8)=0.4564
        !sz(9)=0.2702   !10um-20um in diameter
      
      CASE(3) !Kok et al., 2011 PNAS
      !option 3
      !!based on Kok et al., 2011, we can get another size distribution with lambda=8 and Cv=6.83  -czhao
        sz(1)=1.43e-6
        sz(2)=2.71e-5
        sz(3)=3.57e-4
        sz(4)=3.34e-3
        sz(5)=0.0220
        sz(6)=0.1048
        sz(7)=0.3425
        sz(8)=0.4782
        !sz(9)=0.0484   !10um-20um in diameter

      CASE DEFAULT
       write(0,*) " dust_size_option should be 1,2,3"
       print*, " dust_size_option should be 1,2,3"
       stop
      END SELECT !dust_size_option

       do 1850 n = 1, nsize_aer(itype)
!jchenjin: Only dust for now  !xyruan: Add ca and co3
         l_dust = lptr_dust_aer(n,itype,iphase)
         l_ca = lptr_ca_aer(n,itype,iphase)
         l_co3 = lptr_co3_aer(n,itype,iphase)
!         l_tr1r1 = lptr_tr1r1_aer(n,itype,iphase)
!         l_tr1r2 = lptr_tr1r2_aer(n,itype,iphase)
!         l_tr1r3 = lptr_tr1r3_aer(n,itype,iphase)
!         l_tr1r4 = lptr_tr1r4_aer(n,itype,iphase)
!
         if  (l_ca >= p1st) then
           fracoin = 0.99
           fracca = 0.01*0.4
           fracco3 = 0.01*0.6
           fractot = fracoin
         else
           fracoin = 1.0
           fracca = 0.0
           fracco3 = 0.0
           fractot = fracoin
         endif

         if (l_dust >= p1st) then
                if (nsize_aer(itype).eq.4 .or. nsize_aer(itype).eq.5) then
                    tend_emis(l_dust,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracoin
                if (l_ca > p1st) tend_emis(l_ca,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracca
                if (l_co3 > p1st) tend_emis(l_co3,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracco3
                !if (n <= 3) densdust=2.5
                !if (n > 3 ) densdust=2.65
                densdust = rhosoil/1.e3 !kg/m3 to g/cm3
                ! mass1part is mass of a single particle in ug
                mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                l = numptr_aer(n,itype,iphase)
                if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
         else if (nsize_aer(itype).eq.8 .or. nsize_aer(itype).eq.9 .or. nsize_aer(itype).eq.10) then
                tend_emis(l_dust,k,i) =  factbb * sz(n) * fracoin
                if (l_ca > p1st) tend_emis(l_ca,k,i) = factbb * sz(n) * fracca
                if (l_co3 > p1st) tend_emis(l_co3,k,i) = factbb * sz(n) * fracco3
                !if (n <= 5) densdust=2.5
                !if (n > 5 ) densdust=2.65
                densdust = rhosoil/1.e3 !kg/m3 to g/cm3
                ! mass1part is mass of a single particle in ug
                mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                l = numptr_aer(n,itype,iphase)
                if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + factbb * sz(n) * fractot / mass1part
         endif !nsize_aer
        endif !l_dust
!no source tracking for now chenjin
!         !for source tracking
!          if (l_tr1r1 >= p1st) then
!             if (nsize_aer(itype).eq.4) then
!              if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                  (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!              chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!              factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!              else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                       (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
!              chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!              factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!              else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                       (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!              chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!              factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!              else !region 4 for all others 
!              chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!              factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!              endif
!              if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!              factbb*(sz(2*n-1)+sz(2*n))*fracca
!              if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!              factbb*(sz(2*n-1)+sz(2*n))*fracco3
!              !if (n <= 3) densdust=2.5
!              !if (n > 3 ) densdust=2.65
!              densdust = rhosoil/1.e3 !kg/m3 to g/cm3
!              ! mass1part is mass of a single particle in ug    
!              mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!              l = numptr_aer(n,itype,iphase)
!              if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                         factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
!             else if (nsize_aer(itype).eq.8) then
!              if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                  (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!              chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!              factbb * sz(n) * fracoin
!              else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                       (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
!              chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!              factbb * sz(n) * fracoin
!              else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                       (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!              chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!              factbb * sz(n) * fracoin
!              else    ! region 4 for all others
!              chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!              factbb * sz(n) * fracoin
!              endif
!              if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!              factbb*sz(n)*fracca
!              if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!              factbb*sz(n)*fracco3
!              !if (n <= 5) densdust=2.5
!              !if (n > 5 ) densdust=2.65
!              densdust = rhosoil/1.e3 !kg/m3 to g/cm3
!              ! mass1part is mass of a single particle in ug
!              mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!              l = numptr_aer(n,itype,iphase)
!              if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                        factbb * sz(n) * fractot / mass1part
!             endif !nsize_aer
!            endif !l_tr1r1

1850    continue

 999  continue

     endif ! xland
  enddo ! i

 end subroutine mosaic_dust_Kokemis


  subroutine mosaic_dust_gocartemis (ktau,dt,num_soil_layers,alt,u_phy,&
         v_phy,rho_phy,dz8w,smois,u10,v10,p8w,erod,source_is_ssm,      &
         ivgtyp,isltyp,xland,xlat,xlong,dx,g,emis_dust_factor,         &
         tend_emis,e_dust,                                             &
         its,ite, kts,kte                                     )

  USE module_data_gocart_dust

  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,                   &
                                  its,ite, kts,kte
   INTEGER,DIMENSION( its:ite)                  ,                          &
          INTENT(IN   ) ::  ivgtyp, isltyp

   REAL, DIMENSION( 1:num_chem, kts:kte, its:ite) ,                        &
          INTENT(INOUT) ::   tend_emis
   REAL,  DIMENSION( its:ite)    ,                                         &
          INTENT(INOUT   ) ::    e_dust
   !czhao

   REAL, INTENT(IN) :: emis_dust_factor 

   REAL, DIMENSION(1:num_soil_layers, its:ite) ,                           &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION(1:3, its:ite)                ,                         &
          INTENT(IN   ) ::    erod
   LOGICAL, OPTIONAL, INTENT(IN) :: source_is_ssm
   REAL,  DIMENSION( its:ite)                   ,                          &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland,                &
                                                     dx,                   &
                                                     xlat,                 &
                                                     xlong                 
   REAL,  DIMENSION(kts:kte, its:ite),                                     &
          INTENT(IN   ) :: alt,u_phy,v_phy,rho_phy,dz8w

   REAL,  DIMENSION(kts:kte+1, its:ite),                                   &
          INTENT(IN   ) :: p8w

   REAL, INTENT(IN   ) :: dt,g

!
! local variables
!
  integer :: nmx,i,k,ndt,imx,lmx
  integer ilwi, start_month
  real*8, DIMENSION (3) :: erodin
  real*8, DIMENSION (5) :: bems
  real*8, DIMENSION( its:ite) :: dxy
  real*8  w10m,gwet,airden,airmas,srcfac
  real*8  cdustemis,jdustemis,cdustcon,jdustcon
  real*8  cdustdens,jdustdens,mass1part,jdustdiam,cdustdiam,dp_meanvol_tmp
  real factaa, factbb, fracoin, fracca, fracco3, fractot
  real*8  conver,converi
  real dttt
  real soilfacj,rhosoilj,rhosoilc
  integer p1st,l_dust,l_ca,l_co3,l,n,l_tr1r1,l_tr1r2,l_tr1r3,l_tr1r4
  real :: densdust
  real sz(10),ftot,frac_10um
  real totalemis,accfrac,corfrac,rscale
  integer iphase, itype

  rscale=1.01  !  to account for the emitted dust with radius larger than 10 um
  accfrac=0.15   !dust_emass_jfrac assign 7% to accumulation mode
  corfrac=0.85   !dust_emass_cfrac assign 93% to coarse mode
  accfrac=accfrac*rscale
  corfrac=corfrac*rscale
  p1st = param_first_scalar

!czhao give an example of MOSAIC 8 size bins here (unit: um) from Jerome 2005 paper
! 0.039-0.078; 0.078-0.156; 0.156-0.312; 0.312-0.625; 0.625-1.25; 1.25-2.5; 2.5-5.0; 5.0-10.0
!in the model the size bound for 8 bins is calculated in the code using dlo_sec and dhi_sec
!size3
!       sz(1)=0.0
!       sz(2)=0.0
!       sz(3)=0.0003
!       sz(4)=0.0048
!       sz(5)=0.0130
!       sz(6)=0.0250
!       sz(7)=0.1400
!       sz(8)=0.8169
!the size is following the best match2 of modal distribution
!we have ~3.6% within 1 um radius and ~7.4% within 1.25 um radius
!versus ~6.0% within 1 um radius following GOCART original and ~9.3% within 1.25 um radius
! totally ~82% are within the MOSAIC size range 0.03-10 um diameter
! now the dust distribution follows the parameter defined in module_data_sorgam.F 
        sz(1)=1.0e-8
        sz(2)=1.0e-6
        sz(3)=3.0e-4
        sz(4)=3.5e-3
        sz(5)=0.018
        sz(6)=0.070
        sz(7)=0.2595
        sz(8) =0.42
        sz(9) =0.20
        sz(10)=0.0287
!!based on Kok et al., 2011, we can get another size distribution with lambda=12 and Cv=12.62  -czhao
!       sz(1)=7.7e-7
!       sz(2)=1.5e-5
!       sz(3)=1.9e-4
!       sz(4)=1.8e-3
!       sz(5)=0.0119
!       sz(6)=0.0574
!       sz(7)=0.2018
!       sz(8)=0.4564
!       !sz(9)=0.2702   !10um-20um in diameter


!!based on Kok et al., 2011, we can get another size distribution with lambda=8 and Cv=6.83  -czhao
!       sz(1)=1.43e-6
!       sz(2)=2.71e-5
!       sz(3)=3.57e-4
!       sz(4)=3.34e-3
!       sz(5)=0.0220
!       sz(6)=0.1048
!       sz(7)=0.3425
!       sz(8)=0.4782
!       !sz(9)=0.0484   !10um-20um in diameter


!        frac_10um=sum(sz(1:8))
!        sz(:)=sz(:)/sum(sz(1:8))  ! relative ratio in diameter range of 0-10 um

!  for now just do itype=1
        itype = 1
        iphase = ai_phase

  conver=1.e-9
  converi=1.e9
!
! number of dust bins
  nmx=5
  k=kts
  do i=its,ite


!
! don't do dust over water!!!
     if(xland(i).lt.1.5)then

     ilwi=1
     start_month = 3   ! it doesn't matter, ch_dust is not a month dependent now, a constant
     w10m=sqrt(u10(i)*u10(i)+v10(i)*v10(i))
!    airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g   ! kg 
     airmas=1.0  ! czhao initialize here but not used any more

! we don't trust the u10,v10 values, if model layers are very thin near surface
     if(dz8w(kts,i).lt.12.)w10m=sqrt(u_phy(kts,i)*u_phy(kts,i)+v_phy(kts,i)*v_phy(kts,i))
    !erodin(1)=erod(i,j,1)/dx/dx   ! czhao erod shouldn't be scaled to the area, because it's a fraction
    !erodin(2)=erod(i,j,2)/dx/dx
    !erodin(3)=erod(i,j,3)/dx/dx
    if (present(source_is_ssm) .and. source_is_ssm) then  ! ybxue added dust source function according to config_source_opt
       erodin(1) = max(0.d0, min(1.d0, dble(erod(1,i))))
       erodin(2) = max(0.d0, min(1.d0, dble(erod(1,i))))
       erodin(3) = max(0.d0, min(1.d0, dble(erod(1,i))))
    else
       erodin(1) = erod(1,i)
       erodin(2) = erod(2,i)
       erodin(3) = erod(3,i)
    end if
!  volumetric soil moisture over porosity

!jwfeng : Add porosity through lsm_noah_init.
    porosity(isltyp(i))=maxsmc(isltyp(i))

     gwet=smois(1,i)/porosity(isltyp(i))
     ndt=ifix(dt)
     airden=rho_phy(kts,i)
     dxy(i)=dx(i)*dx(i)

    call mosaic_source_du( nmx, dt, &
                     erodin, ilwi, dxy(i), w10m, gwet, airden, airmas, &
                     bems,start_month,g)

!bems: kg/timestep/cell
    !sum up the dust emission from 0.1-10 um in radius 
    ! unit change from kg/timestep/cell to ug/m2/s
    totalemis=(sum(bems(1:5))/dt)*converi/dxy(i)
    !tune total dust emission, previously it's done for ch_dust
    totalemis=totalemis*emis_dust_factor
    !totalemis=totalemis*rscale !to account for the particles larger than 10 um 
                              ! based on assumed size distribution
    jdustemis = totalemis*accfrac   ! accumulation mode
    cdustemis = totalemis*corfrac   ! coarse mode 

    e_dust(i)=cdustemis+jdustemis
    ftot=jdustemis+cdustemis  ! total dust emission in ug/m2/s 
!    ftot=ftot*frac_10um  ! only in diameter range of 0-10um

!   apportion other inorganics only
! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!        factaa = (dt/dz8w(i,k,j))*alt(i,k,j) 
         factaa = alt(k,i)/dz8w(k,i)
         factbb = factaa * ftot

!         fracoin = 0.99
!         fracca = 0.01*0.4
!         fracco3 = 0.01*0.6
!         fractot = fracoin

         do 1810 n = 1, nsize_aer(itype)

            l_dust = lptr_dust_aer(n,itype,iphase)
!jwfeng: Only dust for now !xyruan: add ca and co3
            l_ca = lptr_ca_aer(n,itype,iphase)
            l_co3 = lptr_co3_aer(n,itype,iphase)
!l_tr1r1 = lptr_tr1r1_aer(n,itype,iphase)
!l_tr1r2 = lptr_tr1r2_aer(n,itype,iphase)
!l_tr1r3 = lptr_tr1r3_aer(n,itype,iphase)
!l_tr1r4 = lptr_tr1r4_aer(n,itype,iphase)

           if  (l_ca >= p1st) then
            fracoin = 0.99
            fracca = 0.01*0.4
            fracco3 = 0.01*0.6
            fractot = fracoin
           else 
             fracoin = 1.0
             fracca = 0.0
             fracco3 = 0.0
             fractot = fracoin
           endif

            if (l_dust >= p1st) then
! jwfeng:04/2020 Add bin 5,9,10
                if (nsize_aer(itype).eq.4 .or. nsize_aer(itype).eq.5) then
                   tend_emis(l_dust,k,i) = factbb * (sz(2*n-1)+sz(2*n)) * fracoin
                   if (l_ca > p1st) tend_emis(l_ca,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracca
                   if (l_co3 > p1st) tend_emis(l_co3,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracco3
                   if (n <= 3) densdust=2.5
                   if (n > 3 ) densdust=2.65
                   ! mass1part is mass of a single particle in ug
                   mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                   l = numptr_aer(n,itype,iphase)
                   if (l >= p1st) &
                       tend_emis(l,k,i) = tend_emis(l,k,i) + factbb * (sz(2*n-1)+sz(2*n)) * fractot/mass1part

                else if (nsize_aer(itype).eq.8 .or. nsize_aer(itype).eq.9 .or. nsize_aer(itype).eq.10) then
                    tend_emis(l_dust,k,i) = factbb * sz(n) * fracoin
                    if (l_ca > p1st) tend_emis(l_ca,k,i) =  factbb * sz(n) * fracca
                    if (l_co3 > p1st) tend_emis(l_co3,k,i) =  factbb * sz(n) * fracco3
                    if (n <= 5) densdust=2.5
                    if (n > 5 ) densdust=2.65
                 !  mass1part is mass of a single particle in ug
                    mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st)  &
                        tend_emis(l,k,i) = tend_emis(l,k,i) + factbb * sz(n) * fractot / mass1part
                endif !numptr_aer
            endif !l_dust

!JWFeng : No source tracking for now !
            !for source tracking
!           if (l_tr1r1 >= p1st) then
!            if (nsize_aer(itype).eq.4) then
!             if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                 (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                      (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then 
!             chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                      (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!             chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else !region 4 for all others 
!             chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             endif
!             if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!             factbb*(sz(2*n-1)+sz(2*n))*fracca
!             if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!             factbb*(sz(2*n-1)+sz(2*n))*fracco3
!             if (n <= 3) densdust=2.5
!             if (n > 3 ) densdust=2.65
!             ! mass1part is mass of a single particle in ug    
!             mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!             l = numptr_aer(n,itype,iphase)
!             if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                        factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
!            else if (nsize_aer(itype).eq.8) then
!             if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                 (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!             factbb * sz(n) * fracoin
!             else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                      (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!             factbb * sz(n) * fracoin
!             else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                      (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!             chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!             factbb * sz(n) * fracoin
!             else    ! region 4 for all others
!             chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!             factbb * sz(n) * fracoin
!             endif
!             if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!             factbb*sz(n)*fracca
!             if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!             factbb*sz(n)*fracco3
!             if (n <= 5) densdust=2.5
!             if (n > 5 ) densdust=2.65
!             ! mass1part is mass of a single particle in ug
!             mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!             l = numptr_aer(n,itype,iphase)
!             if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                       factbb * sz(n) * fractot / mass1part
!            endif !nsize_aer
!           endif !l_tr1r1

1810    continue

     endif ! xland
  enddo ! i

end subroutine mosaic_dust_gocartemis

  SUBROUTINE mosaic_source_du( nmx, dt1, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

 USE module_data_gocart_dust

  INTEGER, INTENT(IN)    :: nmx
  REAL*8,    INTENT(IN)  :: erod(ndcls)
  INTEGER, INTENT(IN)    :: ilwi,month

  REAL*8,    INTENT(IN)    :: w10m, gwet
  REAL*8,    INTENT(IN)    :: dxy
  REAL*8,    INTENT(IN)    :: airden, airmas
  REAL*8,    INTENT(OUT)   :: bems(nmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, n, m, k

  ! default is 1 ug s2 m-5 == 1.e-9 kg s2 m-5
  !ch_dust(:,:)=0.8D-9   ! ch_dust is defined here instead of in the chemics_ini.F if with SORGAM  -czhao
  !ch_dust(:,:)=1.0D-9  ! default 
!  ch_dust(:,:)=0.50D-9  ! ch_dust is scaled for Sahara 
   ch_dust(:,:)=1.0D-9   ! now totalemis is tuned by emis_dust_factor 


  ! executable statemenst
  DO n = 1, nmx
     ! jwfeng:
     ! Use Marticorena and Bergametti (1995). to calculate threshold velocity instead of the original Bagnold(1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              ! Case of surface dry enough to erode
             IF (gwet < 0.5) THEN  !  Pete's modified value
!              IF (gwet < 0.2) THEN
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)*dxy  ! (m2)
              IF (ilwi == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m**2 &
                      * (w10m - u_ts)*dt1  ! (kg)
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              ! Update dust mixing ratio at first model level.
              !tc(n) = tc(n) + dsrc / airmas    !kg/kg-dryair -czhao
              bems(n) = dsrc     ! kg/timestep/cell

  ENDDO

END SUBROUTINE mosaic_source_du


subroutine mosaic_dust_shaoemis(ktau,dt,  num_soil_layers,        & 
    rho_phy,dz8w,smois,ust,                                       &
    isltyp,vegfra,g,                                              &
    imod,shao_v,num_soil,                                         &
    soil_top_cat, erod,alt,source_is_ssm,                         &
    tend_emis,e_dust,xlat,xlong,                                  &
    its,ite, kts,kte                                         )

    USE module_data_gocart_dust
    USE module_data_binaer_asect
   !phys/module_sf_noahmplsm/NOAHMP_TABLES
   !  USE NOAHMP_TABLES, ONLY: DRYSMC_nmp => SMCDRY_TABLE
   !  USE module_data_sorgam
   !  USE module_data_gocart_dust
   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) :: ktau, imod,shao_v,num_soil, num_soil_layers,   &
                                its,ite, kts,kte
   INTEGER, DIMENSION( its:ite) ,                                                &
        INTENT(IN   ) ::                                 isltyp
   REAL, DIMENSION(1:num_soil,its:ite) ,                                         &
        INTENT(IN   ) ::                           soil_top_cat
   REAL,  DIMENSION(1:3, its:ite) ,                                              &
   INTENT(IN   ) ::                           erod
   LOGICAL, OPTIONAL, INTENT(IN) :: source_is_ssm
   REAL, DIMENSION( 1:num_chem, kts:kte, its:ite) ,                              &
   INTENT(INOUT) ::   tend_emis
   ! REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,                  &
   !       INTENT(INOUT ) ::                              emis_dust
   REAL, DIMENSION(1:num_soil_layers, its:ite) ,                                 &
         INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( its:ite)  ,                                                 &
        INTENT(IN   ) ::                                ust,vegfra
   REAL,  DIMENSION( kts:kte, its:ite),                                          &
        INTENT(IN   ) ::                           rho_phy, alt,dz8w
   REAL, DIMENSION( its:ite), INTENT(INOUT ) :: e_dust
   REAL,  DIMENSION( its:ite)                   ,                                &
        INTENT(IN  )::                              xlat, xlong
   REAL, INTENT(IN   ) :: dt,g
   !----------------------------------------------------------------------------
   ! local variables
   !---------------------------------------------------------------------------
   ! ust_t_min is calculated value from qf03
   REAL,  DIMENSION( its:ite) ::          ust_t_min,                             &
                                          rough_cor,                             &
                                          smois_cor
   integer, parameter :: imax=100     ! No. of particle size intervalsfor psd
   integer, parameter :: jmax=4       ! No. of log-normal distributionsfor constructing psd
   real(8), dimension(0:imax) :: d0
   real(8), dimension(imax)   :: dd
   real(8), dimension(imax)   :: psdm, dpsdm, ppsdm
   real(8), dimension(imax)   :: psdf, dpsdf, ppsdf
   real(8), parameter :: dcut=20.d0            ! dust cutoff particle size
   real(8), parameter :: rhop = 2650.d0        ! particle density [kg/m3]
   
   integer             :: i,j,k,p,idst
   integer, parameter  :: nmx = 5  ! No. of dust bins
   real                :: ust_grid, airden, dz_lowest, smc
   real, dimension(nmx):: tc, bems
   real(8)             :: gwet, cf, srcfac
   real(8)             :: conver,converi
   real(8)             :: ust_min, rough_cor_in, smois_cor_in
   real, dimension(16) :: soilc
   real                :: tot_soilc
   integer             :: domsoilc
   integer             :: cc
   character*1 :: tmp
   real(8)             :: ustart0_out, lambda, sigma
   real(8), dimension(imax) :: ustart
   real, dimension(12) :: thr
   data thr/0.001, 0.003, 0.037, 0.061, 0.072, 0.049,0.084,0.110,0.095,0.126,0.141, 0.156/
   
   real factaa, factbb, fracoin, fracca, fracco3, fractot
   integer l,n,p1st,l_dust,l_ca,l_co3,l_tr1r1,l_tr1r2,l_tr1r3,l_tr1r4
   integer iphase, itype
   real totalemis,ftot
   real sz(10),frac_10um
   real :: densdust, mass1part
   integer,parameter :: dust_size_option = 2
   integer           :: stype        ! No. of soil texture classes
   real(8), dimension(imax,12)   :: psd_m, dpsd_m, ppsd_m
   real(8), dimension(imax,12)   :: psd_f, dpsd_f, ppsd_f
   logical :: is_dust_source
   !-----------------------------------------------------------------------------------------
   if(shao_v .eq.1)then
      stype= 4
   elseif(shao_v .eq.2)then
      stype= 12          ! shao_v==2 represent select WRF-Chemv4.0     
   endif
   !size distribution
   SELECT CASE (dust_size_option)  !option for dust size distribution 
    CASE(1)  !Zhao et al., 2010  ACP 
    !option 1
    !!the size is following the best match2 of modal distribution
    !!we have ~3.6% within 1 um radius and ~7.4% within 1.25 um radius
    !!versus ~6.0% within 1 um radius following GOCART original and ~9.3%
    !within 1.25 um radius
    !!totally ~82% are within the MOSAIC size range 0.03-10 um diameter
    !!now the dust distribution follows the parameter defined in
    !module_data_sorgam.F 
      sz(1)=1.0e-8
      sz(2)=1.0e-6
      sz(3)=3.0e-4
      sz(4)=3.5e-3
      sz(5)=0.018
      sz(6)=0.070
      sz(7)=0.2595
      sz(8)=0.42
      sz(9)=0.20
      sz(10)=max(0.0, 1.0-sum(sz(1:9)))
    CASE(2) !Kok et al., 2011 PNAS 
    !option 2
    !!based on Kok et al., 2011, we can get another size distribution with
    !lambda=12 and Cv=12.62  -czhao
      sz(1)=7.7e-7
      sz(2)=1.5e-5
      sz(3)=1.9e-4
      sz(4)=1.8e-3
      sz(5)=0.0119
      sz(6)=0.0574
      sz(7)=0.2018
      sz(8)=0.4564
      sz(9)=0.2702   !10um-20um in diameter
      sz(10)=max(0.0, 1.0-sum(sz(1:9)))
   
    CASE(3) !Kok et al., 2011 PNAS
    !option 3
    !!based on Kok et al., 2011, we can get another size distribution with
    !lambda=8 and Cv=6.83  -czhao
      sz(1)=1.43e-6
      sz(2)=2.71e-5
      sz(3)=3.57e-4
      sz(4)=3.34e-3
      sz(5)=0.0220
      sz(6)=0.1048
      sz(7)=0.3425
      sz(8)=0.4782
      sz(9)=0.0484   !10um-20um in diameter
      sz(10)=max(0.0, 1.0-sum(sz(1:9)))
    CASE DEFAULT
     write(0,*) " dust_size_option should be 1,2,3"
     print*, " dust_size_option should be 1,2,3"
     stop
    END SELECT !dust_size_option
   !initialization
   conver=1.e-9
   converi=1.e9
   p1st = param_first_scalar
   itype = 1
   iphase = ai_phase
   k=kts   !in the bottom layer
   ! calculate soil-psd once for all 4/12 types
   ! initialize
   psd_m(:,:) = 0.
   psd_f(:,:) = 0.
   dpsd_m(:,:) = 0.
   dpsd_f(:,:) = 0.
   ppsd_m(:,:) = 0.
   ppsd_f(:,:) = 0.
   if(shao_v .eq. 1)then
   do p = 1, stype
    if (p.eq.1) then            ! sand
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csandm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csandf,jmax,shao_v)
    elseif (p.eq.2) then        ! sandy clay loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csloam,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csloaf,jmax,shao_v)
    elseif (p.eq.3) then        ! loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,cloamm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,cloamf,jmax,shao_v)
    else                        ! clay
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,cclaym,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,cclayf,jmax,shao_v)
    endif
       psd_m(:,p) = psdm(:)
       psd_f(:,p) = psdf(:)
       dpsd_m(:,p) = dpsdm(:)
       dpsd_f(:,p) = dpsdf(:)
       ppsd_m(:,p) = ppsdm(:)
       ppsd_f(:,p) = ppsdf(:)
   enddo
   endif
   if(shao_v .eq.2)then
   do p = 1, stype
   !     write(6,*) p
    if (p.eq.1) then            ! sand
       call psd_create(d0, dd, psdm, dpsdm, ppsdm,imax,vcsandm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf,imax,vcsandf,jmax,shao_v)
    elseif (p.eq.2) then        ! loamy sand
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,closam,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,closaf,jmax,shao_v)
    elseif (p.eq.3) then        ! sandy loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csalom,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csalof,jmax,shao_v)
    elseif (p.eq.4) then        ! silt loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csilom,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csilof,jmax,shao_v)
    elseif (p.eq.5) then        ! silt
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csiltm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csiltf,jmax,shao_v)
    elseif (p.eq.6) then        ! loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm,imax,vcloamm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf,imax,vcloamf,jmax,shao_v)
    elseif (p.eq.7) then        ! sandy clay loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csclom,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csclof,jmax,shao_v)
    elseif (p.eq.8) then        ! silty clay loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csiclm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csiclf,jmax,shao_v)
    elseif (p.eq.9) then        ! clay loam
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,ccloam,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,ccloaf,jmax,shao_v)
    elseif (p.eq.10) then        ! sandy clay
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csaclm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csaclf,jmax,shao_v)
    elseif (p.eq.11) then        ! silty clay
       call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax,csilcm,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax,csilcf,jmax,shao_v)
    elseif (p.eq.12) then        ! clay
       call psd_create(d0, dd, psdm, dpsdm, ppsdm,imax,vcclaym,jmax,shao_v)
       call psd_create(d0, dd, psdf, dpsdf, ppsdf,imax,vcclayf,jmax,shao_v)
    endif
       psd_m(:,p) = psdm(:)
       psd_f(:,p) = psdf(:)
       dpsd_m(:,p) = dpsdm(:)
       dpsd_f(:,p) = dpsdf(:)
       ppsd_m(:,p) = ppsdm(:)
       ppsd_f(:,p) = ppsdf(:)
   enddo
   endif
   !----------------------------------------------------------------------------------
   ! Before calculating dust emission, some parameters should be set.
   j = 0
   1  j = j+1
   if ( dd(j) .le. dcut ) then
    idst = j
    goto 1
   endif
   do i=its,ite
   !do dust over dust source area only
   !if (sum(erod(:,i)).gt.0.) then    !use fraction of erodible surface areaas dust source indicator     
   is_dust_source = .false.           ! ybxue added dust source function according to config_source_opt
   if (present(source_is_ssm) .and. source_is_ssm) then
      srcfac = max(0.d0, min(1.d0, dble(erod(1,i))))
      if (srcfac .gt. 0.d0) is_dust_source = .true.
   else
      if (sum(erod(:,i)) .gt. 0.d0) is_dust_source = .true.
   end if
   if (is_dust_source) then
     !tc(1)=chem(i,kts,j,p_dust_1)*conver   ![kg/kg-dryair]
     !tc(2)=chem(i,kts,j,p_dust_2)*conver
     !tc(3)=chem(i,kts,j,p_dust_3)*conver
     !tc(4)=chem(i,kts,j,p_dust_4)*conver
     !tc(5)=chem(i,kts,j,p_dust_5)*conver
      ust_grid=ust(i)   ! u* for one grid
   ! delta z for the lowest layer is needed for unit conversion
      dz_lowest = dz8w(1,i)
   ! for soil moisture using volumetric soil moisture (smois)
      if (smois(1,i).gt.1.) then
         smois(1,i) = 1.
         !CALL wrf_message('UoC CTDE WARNING: vol. soil moisture > 1, reset')
      endif
      gwet=smois(1,i)
      airden=rho_phy(kts,i)                       ![kg/m3]
      cf=vegfra(i)                                ! in [%]
      cf = cf/100.d0                                ! fraction 
   !initialization
       tot_soilc=0.
       do cc = 1, 12
          soilc(cc) = 0.
       enddo
   
       do cc = 1, 12
          soilc(cc) = soil_top_cat(cc,i)
          tot_soilc = tot_soilc + soilc(cc)
       enddo
   !
       domsoilc = isltyp(i)
   !set DRYSMC = thr chenjin
       !if (  config_flags%sf_surface_physics .eq. 3   ) then
       !    DRYSMC = DRYSMC_ruc
      !! elseif (  config_flags%sf_surface_physics .eq. 4   ) then
      !!     DRYSMC = DRYSMC_nmp
       !elseif ( config_flags%sf_surface_physics .eq. 1 .or. &
       !&        config_flags%sf_surface_physics .eq. 5 .or. &
       !&        config_flags%sf_surface_physics .eq. 7 .or. &
       !&        config_flags%sf_surface_physics .eq. 8 .or. &
       !&        config_flags%sf_surface_physics .eq. 0) then
           DRYSMC(1:12) = thr
       !    CALL wrf_message('UoC dust: DRYSMC reset for dust emission')
       !endif
       call qf03_driver( nmx, idst, g, rhop, airden, dt,        &
                ust_grid, gwet, cf, ust_min, imod,shao_v, dz_lowest, &
                soilc, tot_soilc, domsoilc,                   &
                tc, bems, rough_cor_in, smois_cor_in, DRYSMC(1:12),    &
                d0, dd, psd_m, dpsd_m, ppsd_m, psd_f, dpsd_f,ppsd_f,imax,stype)
       !frac_10um=sum(sz(1:8))
       !sz(:)=sz(:)/sum(sz(1:8))
       totalemis=sum(bems(1:5))*converi      !kg/m2/s >> ug/m2/s
       !tot= totalemis*frac_10um
       e_dust(i)=totalemis
       ! chenjin: Using tend_emis instead of directly changing chem
       factaa = alt(kts,i) / dz_lowest !ug/m2/s>>ug/kg-dryair
       factbb = factaa * totalemis
       !factbb = sum(tc(1:5))*converi
   do n = 1, nsize_aer(itype)
   
       l_dust = lptr_dust_aer(n,itype,iphase)
       l_ca = lptr_ca_aer(n,itype,iphase)
       l_co3 = lptr_co3_aer(n,itype,iphase)
   !    l_tr1r1 = lptr_tr1r1_aer(n,itype,iphase)
   !    l_tr1r2 = lptr_tr1r2_aer(n,itype,iphase)
   !    l_tr1r3 = lptr_tr1r3_aer(n,itype,iphase)
   !    l_tr1r4 = lptr_tr1r4_aer(n,itype,iphase)
   
      if  (l_ca >= p1st) then
        fracoin = 0.99
        fracca = 0.01*0.4
        fracco3 = 0.01*0.6
        fractot = fracoin
      else
         fracoin = 1.0
         fracca = 0.0
         fracco3 = 0.0
         fractot = fracoin
       endif
   
       if (l_dust >= p1st) then
        if (nsize_aer(itype).eq.4) then
           tend_emis(l_dust,k,i) = factbb * (sz(2*n-1)+sz(2*n)) * fracoin
           if (l_ca > p1st) tend_emis(l_ca,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracca
           if (l_co3 > p1st) tend_emis(l_co3,k,i) =  factbb * (sz(2*n-1)+sz(2*n)) * fracco3
        if (n <= 3) densdust=2.5
        if (n > 3 ) densdust=2.65
        !densdust = rhosoil/1.e3 !kg/m3 to g/cm3
        ! mass1part is mass of a single particle in ug
        mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
        l = numptr_aer(n,itype,iphase)
        if (l >= p1st) tend_emis(l,k,i) = factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
        else if (nsize_aer(itype).eq.8 .or. nsize_aer(itype).eq.9 .or. nsize_aer(itype).eq.10) then
           tend_emis(l_dust,k,i) = factbb * sz(n) * fracoin
           if (l_ca > p1st) tend_emis(l_ca,k,i) =  factbb * sz(n) * fracca
           if (l_co3 > p1st) tend_emis(l_co3,k,i) =  factbb * sz(n) * fracco3
        if (n <= 5) densdust=2.5
        if (n > 5 ) densdust=2.65
        !densdust = rhosoil/1.e3 !kg/m3 to g/cm3
        ! mass1part is mass of a single particle in ug
        mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
        l = numptr_aer(n,itype,iphase)
       if (l >= p1st) tend_emis(l,k,i) = factbb * sz(n) * fractot / mass1part
        endif !nsize_aer
       endif !l_dust
   enddo ! n
       !for source tracking
   !no source tracking for now chenjin
   !     if (l_tr1r1 >= p1st) then
   !        if (nsize_aer(itype).eq.4) then
   !         if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1for US 
   !             (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
   !         chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
   !         factbb * (sz(2*n-1)+sz(2*n)) * fracoin
   !         else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2for East Asia
   !                  (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
   !         chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
   !         factbb * (sz(2*n-1)+sz(2*n)) * fracoin
   !         else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. &!region3for North Africa 
   !                  (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
   !         chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
   !         factbb * (sz(2*n-1)+sz(2*n)) * fracoin
   !         else !region 4 for all others 
   !         chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
   !         factbb * (sz(2*n-1)+sz(2*n)) * fracoin
   !         endif
   !         if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
   !         factbb*(sz(2*n-1)+sz(2*n))*fracca
   !         if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
   !         factbb*(sz(2*n-1)+sz(2*n))*fracco3
   !         if (n <= 3) densdust=2.5
   !         if (n > 3 ) densdust=2.65
   !         !densdust = rhosoil/1.e3 !kg/m3 to g/cm3
   !         ! mass1part is mass of a single particle in ug    
   !         mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
   !         l = numptr_aer(n,itype,iphase)
   !         if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
   !                    factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
   !        else if (nsize_aer(itype).eq.8) then
   !         if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1for US 
   !             (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
   !         chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
   !         factbb * sz(n) * fracoin
   !         else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2for East Asia
   !                  (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
   !         chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
   !         factbb * sz(n) * fracoin
   !         else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3for North Africa 
   !                  (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
   !         chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
   !         factbb * sz(n) * fracoin
   !         else    ! region 4 for all others
   !         chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
   !         factbb * sz(n) * fracoin
   !         endif
   !         if (l_ca >= p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
   !         factbb*sz(n)*fracca
   !         if (l_co3 >= p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
   !         factbb*sz(n)*fracco3
   !         if (n <= 5) densdust=2.5
   !         if (n > 5 ) densdust=2.65
   !         !densdust = rhosoil/1.e3 !kg/m3 to g/cm3
   !         ! mass1part is mass of a single particle in ug
   !         mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
   !         l = numptr_aer(n,itype,iphase)
   !         if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
   !                   factbb * sz(n) * fractot / mass1part
   !        endif !nsize_aer
   !       endif !l_tr1r1
   
   
     else   ! no dust sourc
          e_dust(i)=0.
          ust_min = -9999.d0
          rough_cor_in = 1.d0
          smois_cor_in = 1.d0
   !          endif
   endif !dsource/erod
      ust_t_min(i) = ust_min
      rough_cor(i) = rough_cor_in
      smois_cor(i) = smois_cor_in
   
   enddo ! i
   
   end subroutine mosaic_dust_shaoemis
!*********************************************************************************
!----------------------------------------------------------------------------
! Yaping Shao, 13 June 2000
!
! - Generate particle size distribution density function 
!   (both minimally-dispersed and fully-dispersed as the 
!   sum of four log-normal distributions.
!
! d(0,imax):    output, particle size at 0, 1, 2, ..., imax points      [um]
! dm(imax):     output, particle size at middle of 0-1, 1-2, etc      [um]
! psd(imax):    output, particle size distribution density at dm        [um^-1]
! dpsd(imax):   output, Delta P for sections 0-1, 1-2, etc.             [ ]
! ppsd(imax):   output, P for sections 0-1, 1-2, etc.           [ ]
! imax:         input, length dm, psd, dpsdm, ppsd, etc.
! cmtrix:       jmaxx coefficient matrix 
!               e.g.
!               w1   = cmtrix(1, 1): weight for first log-normal distribution
!               dln1 = cmtrix(2, 1): mean log-particle size of first log-normal
!               distribution    
!               sig1 = cmtrix(3, 1): sigma of log-particle size for first
!               log-normal distribution
!               careful with the dimension of dln and sig
!----------------------------------------------------------------------------
subroutine psd_create(d, dm, psd, dpsd, ppsd, imax, cmtrix, jmax,ver)
    integer :: i, j, imax, jmax,ver
    real(8), dimension(3, jmax) :: cmtrix
    real(8) :: d(0:imax), dm(imax)
    real(8) :: psd(imax), dpsd(imax), ppsd(imax)   ! for p(d), Delta P(d) andP(d)
    real(8) :: p, pp, w, dln, sig
    real(8) :: cn
    real(8), parameter :: eps=1.d-7
    real(8), parameter :: dref=1000.d0
    real(8) :: fu, fd, phi
    cn = 1.d0/dsqrt(2.d0*3.14159d0)
!
! initialise psd, dpsd, ppsd
!
    psd = 0.d0
    dpsd = 0.d0
    ppsd = 0.d0
! Estimate d using phi scale. phi varies between from 9 to -1
! with increment 0.1. Reference particle size d0 = 1000 um
    fu = 10.d0
    fd = -1.d0
    do i = 0, imax
      phi = fu - i*(fu-fd)/imax
      d(i) = dref/2.d0**phi
    enddo
    do i = 1, imax
    if(ver .eq.1)then
    dm(i) = dexp( (dlog(d(i))+dlog(d(i-1)) )/2.d0 )
    elseif(ver .eq.2)then
    dm(i) = (d(i)-d(i-1))/dlog(d(i)/d(i-1)) ! d = delta(d)/delta(logd),because del (log d) = (del d) / d   ! mklose [17122014]
    endif
      pp = 0.d0
      do j = 1, jmax
        w = cmtrix(1, j)
        dln = cmtrix(2, j)
        sig = cmtrix(3, j)
        if ( (w.gt.eps) .and. (sig.ne.0.) ) then
          p = w*cn/sig*dexp( -(dlog(dm(i))-dln)**2/(2*sig**2) )
        else
          p=0.d0
        endif
        pp = pp + p
      enddo
!
      dpsd(i) = pp*( dlog(d(i)) - dlog(d(i-1)) )        ! Delta P over i
      if (i.eq.1) then
        ppsd(i) = 0.d0 + dpsd(i)                        ! P(d), with P(0) = 0
      else
        ppsd(i) = ppsd(i-1) + dpsd(i)
      endif
      psd(i) = pp/dm(i)           ! p(d), particle sizedistribution density
    enddo
!
! Normalize in case ppsd(imax) is not 1
!
    dpsd = dpsd/ppsd(imax)
    psd  = psd/ppsd(imax)
    ppsd = ppsd/ppsd(imax)

   end subroutine psd_create
!************************************************************************************************
! Y. Shao, 29 Jan 2004 
!         
! JY Kang, 01 Dec 2008
! Modify the code for WRF_chem
!
! M. Klose, 2010-2015 Modifications
!-----------------------------------------------------------------------------------
! Calculate sediment flux for multi-particle size soils as a weighted average of
! Q(d) 
!           dust emission F(d) for covered and moisture soil
!
! Options for dust calculation:
! 1 Shao (2001)
! 2 Shao (2004)
! 3 Shao (2011): simplification of 2; added on 26 Sep 2009
!
!--------------------------------------------------------------------------------------
!
!      input: 
!      n:          number of particle size ranges.
!      dm:         median diameter of each particle size.                   [m]
!      m_fract:    Weight fraction of each particle range. Sum m_fract = 1
!      ustar:      Friciton velocity.
!      [m/s]
!      cf:         fraction area covered by roughness elements/vegetation cover
!      w:          surface soil moisture content [m^3/m^3]
!      wr:         air-dry soil moisture [m^3/m^3] from SOILPARM.TBL
!      c:          Owen's coefficient
!      output:  
!      ustart:     Mean threshold velocity of each particle range.
!      q:          Sand flux from each size range.
!      ffq:        Weighted sand flux from each size range.
!      qtotal:     Total sand flux (weighted average).
!      f:          Weighted dust flux ejected by single sand size range for
!      given d_d. 
!      fff:        Weighted dust flux ejected by all sand size range from each
!      dust size range.
!      ftotal:     Total dust flux (weighted average).
!--------------------------------------------------------------------------------------
subroutine qf03_driver ( nmx, idst, g, rhop, rho, dt, &
                  ustar, w, cf, ust_min, imod,shao_v, dz_lowest, &
                  soilc, tot_soilc, domsoilc,    &
                  tc, bems, rough_cor_in, smois_cor_in, wr, &
                  d0, dd, psd_m, dpsd_m, ppsd_m, psd_f, dpsd_f,ppsd_f,imax,stype)
     INTEGER, INTENT(IN   )                   :: nmx, imod,shao_v, idst,stype
     REAL   , INTENT(IN   )                   :: g, dt, rho, dz_lowest, ustar
     REAL(8), INTENT(IN   )                   :: w
     REAL(8), INTENT(INOUT)                   :: cf
     REAL(8), INTENT(  OUT)                   :: ust_min,rough_cor_in,smois_cor_in
     REAL   , INTENT(INOUT), DIMENSION(nmx)   :: tc
     REAL   , INTENT(  OUT), DIMENSION(nmx)   :: bems
     REAL   , INTENT(IN   ), DIMENSION(12)    :: wr

!kang [2009/01/07] soil class type
     REAL, INTENT(IN   )                      :: tot_soilc
     REAL, DIMENSION(16), INTENT(IN)          :: soilc
     INTEGER,  INTENT(IN  )                   :: domsoilc

     integer :: imax
     real(8), dimension(0:imax), intent(in) :: d0
     real(8), dimension(imax), intent(in)   :: dd
     real(8), dimension(imax,12), intent(in)   :: psd_m, dpsd_m, ppsd_m
     real(8), dimension(imax,12), intent(in)   :: psd_f, dpsd_f, ppsd_f
     real(8), intent(in) :: rhop
!local variables
!      particle-size distributions
     integer, parameter :: nbins=5 
     real(8), dimension(imax)      :: psdm, dpsdm, ppsdm    ! minimally dispersed
     real(8), dimension(imax)      :: psdf, dpsdf, ppsdf    ! fully dispersed
     real(8), dimension(imax)      :: psds, dpsds, ppsds    ! sediment       
     character(4) :: s_type
     real(8) :: smois_correc
     integer :: i, j, n, kk, ij, ffile
     real(8) :: total, qtotal, ftotal
     real(8) :: ftotalb, ftotalp
!
     real(8), dimension(imax)   :: beta_d, beta_s  ! beta1 and beta2 used byShao et al. (1996) dust emission
     real(8), dimension(imax)   :: ustart, q, ffq, fff
     real(8), dimension(imax,imax) :: f
     real(8), parameter :: c_lambda=0.35
     real(8) :: h, lambda
     real(8) :: ghl, fc
     real(8) :: phl
     real(8) :: cys, u0, al0, rhos, smass, omega, rys
     real(8) :: ddm, a1, a2, a3
     real(8) :: zeta, sigma_m                      !u*sqrt(rhos/p),bombardment coefficient
     real(8) :: ustart0_out, qwhite_out, f_mb_out, f_hlys_out,pmass_out,vhlys_out
     real(8) :: sigma
     real    :: cell_ftotal
     integer :: isl, cc
     character :: msg, soil
     real, dimension(nbins)    :: dbin, fbin, cell_fbin,dbin_v4,dbin_v3
     data dbin_v4/2.,3.6,6.,12.,20./ !size cut diameter (um) consistent withGOCART model  v4.0
     data dbin_v3/2.5,5.,10.,20.,0./    !size cut diameter (um)   v3.7.1
     integer, dimension(nbins) :: ibin
!initialize
     cell_ftotal = 0.
     do n = 1, nbins
     cell_fbin(n) = 0.
     enddo
DO cc = 1, 12  ! soil category
    if (soilc(cc).eq.0.) then
            go to 103
    endif
  if(shao_v .eq.1)then
     dbin= dbin_v3
     if (cc.eq.1.or.cc.eq.2) then
      s_type = 'sand'
         psdm(:)=psd_m(:,1)
         psdf(:)=psd_f(:,1)
         dpsdm(:)=dpsd_m(:,1)
         dpsdf(:)=dpsd_f(:,1)
         ppsdm(:)=ppsd_m(:,1)
         ppsdf(:)=ppsd_f(:,1)
     elseif (cc.eq.3.or.cc.eq.4..or.cc.eq.6.or.cc.eq.8.or.cc.eq.9) then
      s_type = 'loam'
         psdm(:)=psd_m(:,3)
         psdf(:)=psd_f(:,3)
         dpsdm(:)=dpsd_m(:,3)
         dpsdf(:)=dpsd_f(:,3)
         ppsdm(:)=ppsd_m(:,3)
         ppsdf(:)=ppsd_f(:,3)
     elseif (cc.eq.7) then
      s_type = 'sloa'
         psdm(:)=psd_m(:,2)
         psdf(:)=psd_f(:,2)
         dpsdm(:)=dpsd_m(:,2)
         dpsdf(:)=dpsd_f(:,2)
         ppsdm(:)=ppsd_m(:,2)
         ppsdf(:)=ppsd_f(:,2)
     elseif (cc.eq.5.or.cc.eq.10.or.cc.eq.11.or.cc.eq.12) then
      s_type = 'clay'
         psdm(:)=psd_m(:,4)
         psdf(:)=psd_f(:,4)
         dpsdm(:)=dpsd_m(:,4)
         dpsdf(:)=dpsd_f(:,4)
         ppsdm(:)=ppsd_m(:,4)
         ppsdf(:)=ppsd_f(:,4)
     else
      go to 103
     endif
elseif(shao_v .eq.2)then
     dbin=dbin_v4
         psdm(:)=psd_m(:,cc)
         psdf(:)=psd_f(:,cc)
         dpsdm(:)=dpsd_m(:,cc)
         dpsdf(:)=dpsd_f(:,cc)
         ppsdm(:)=ppsd_m(:,cc)
         ppsdf(:)=ppsd_f(:,cc)
! frontal area index
       if (cf .gt. 0.95) then
          write(6,*) 'cover fraction too large, reset'
          cf = 0.95
       endif! as larger values lead to large lambda, which is problematic in Raupach scheme
endif

! default settings: 
     rhos = 1000.d0      ! bulk density of soil [kg/m3] 
     phl = 30000.        ! plastic pressure [N/m2]
     cys = 0.00001       ! cys : parameter

     sigma = rhop/rho    ! particle-air density ratio
!-------------------
! frontal area index
!-------------------
       lambda = - c_lambda*dlog( 1.d0 - cf )
       call r_c(lambda, rough_cor_in,shao_v)
! Matching WRF_soil class with Shao's class for moisture correction of Fecan
! (cc:WRF_soil class, isl:Shao's class)
! cc
! 1:sand, 2:loamy sand, 3:sandy loam, 4:silt loam, 5:silt, 6:loam, 7:sandy clay
! loam,
! 8:silty clay loam, 9:clay loam, 10:sandy clay, 11:silty clay, 12:clay
! isl
! 1:sand, 2:loamy sand, 3:sandy loam, 4:loam, 5:silt loam, 6: silt, 7:sandy clay
! loam, 
! 8:clay loam, 9:silty clay loam, 10:sandy clay, 11:silty clay, 12:clay
      if (cc.eq.1.or.cc.eq.2.or.cc.eq.3.or.cc.eq.7.or.cc.eq.10.or.  &
      &   cc.eq.11.or.cc.eq.12) then
       isl = cc
      elseif (cc.eq.4) then
       isl = 5
      elseif (cc.eq.5) then
       isl = 6
      elseif (cc.eq.6) then
       isl = 4
      elseif (cc.eq.8) then
       isl = 9
      elseif (cc.eq.9) then
       isl = 8
      endif

       call h_c(w, wr(cc), isl, smois_correc)
!----------------------------------------------
! for each particle size group, estimate ustart
!----------------------------------------------
!
       ust_min = 999.0
       do i = 1, imax
          call ustart0(dd(i), sigma, g, rho, ustart0_out)
          ustart(i) = ustart0_out
          ustart(i) = rough_cor_in*smois_correc*ustart(i)
          ust_min = dmin1(ust_min, ustart(i))
          call qwhite(ustart(i), ustar, rho, g, qwhite_out,shao_v)
          q(i) = qwhite_out
          q(i) = (1.d0-cf)*q(i)
       enddo
       if (cc.eq.domsoilc) then
        smois_cor_in = smois_correc
       endif
! 
! 
       IF ( ustar .le. ust_min ) THEN    ! no erosion goto 102
          q = 0.d0
          ffq = 0.d0
          qtotal = 0.d0
          fff = 0.d0
          ftotal = 0.d0
          fbin = 0.d0
          goto 102
       ELSE
          ghl = dexp( -(ustar - ust_min)**3.d0 )
          dpsds = ghl*dpsdm + (1.-ghl)*dpsdf
          psds  = ghl*psdm + (1.-ghl)*psdf
          ppsds = ghl*ppsdm + (1.-ghl)*ppsdf
!
          ffq = q*dpsds
          qtotal = sum(ffq)

!--------------
! dust emission
!--------------
! 
! size bin
          do n=1,nbins
           ibin(n)=0
           do i=imax,1,-1
            if(d0(i).ge.dbin(n)) ibin(n)=i
           enddo
           if(ibin(n).eq.0) stop 'wrong dust classes'
          enddo
!--------------------------------
! Shao (2001) dust emission model
!--------------------------------
      IF (imod .eq. 1) THEN
         do i = idst+1, imax
            ddm = dd(i)*1.d-6
            call pmass(rhop, ddm, pmass_out)                       ! mass of saltating particles
            smass = pmass_out
            u0 = 10*ustar
            al0 = 13.d0*3.14159d0/180.d0
            call vhlys(phl, 2, smass, al0, u0, ddm, vhlys_out)
             omega = vhlys_out     ![m3]
!               
            do j = 1, idst
               rys = psdm(j)/psdf(j)
               a1 = cys*( (1.-ghl) + ghl*rys )
               a2 = ffq(i)*g/ustar**2/smass

               if ( dpsdf(j) .lt. dpsdm(j) ) then
                  a3 = dpsdf(j)*rhos*omega
               else
                  a3 = dpsdf(j)*rhos*omega + (dpsdf(j)-dpsdm(j))*smass
               endif
               f(i,j) = a1*a2*a3    ![kg/m2/s]
            enddo
         enddo
         ftotal = 0.0
         do j = 1, idst
            fff(j) = 0.
            do i = idst+1, imax
               fff(j) = fff(j) + f(i,j)
            enddo
            fff(j) = (1-cf)*fff(j)
            ftotal = ftotal + fff(j)
         enddo
! 
         do n=1,nbins
          j0=1
          if(n.gt.1) j0=ibin(n-1)+1
          fbin(n)=0
          do j=j0,ibin(n)
           fbin(n)=fbin(n)+fff(j)
          enddo
         enddo
!--------------------------------
! Shao (2004) dust emission model
!--------------------------------
      ELSEIF (imod .eq. 2) THEN
        zeta = ustar*dsqrt( rhos/phl )
        sigma_m = 12.d0*zeta*zeta*(1.d0+14.d0*zeta)
! 
        do i = idst+1, imax
           do j = 1, idst
              rys = psdm(j)/psdf(j)
              a1 = cys*dpsdf(j)*( (1.-ghl) + ghl*rys )
              a2 = (1.+sigma_m)
              a3 = ffq(i)*g/ustar**2
              f(i,j) = a1*a2*a3
           enddo
        enddo
!           
        ftotal = 0.0
        do j = 1, idst
           fff(j) = 0.
           do i = idst+1, imax
              fff(j) = fff(j) + f(i,j)
           enddo
           fff(j) = (1-cf)*fff(j)
           ftotal = ftotal + fff(j)
        enddo
        do n=1,nbins
         j0=1
         if(n.gt.1) j0=ibin(n-1)+1
         fbin(n)=0
         do j=j0,ibin(n)
          fbin(n)=fbin(n)+fff(j)
         enddo
        enddo
!--------------------------------------------------------------------------
! Shao (2011) minimal version, ghl = 1, Q independent of sand particle size
! 
! See Eq. (34) in 
! Shao, Y., M. Ishizuka, M. Mikami, J. Leys (2011): Parameterization of size-
! resolved dust emission and validation with measurements, JGR, 116, D08203, 
! doi: 10.1029/2010JD014527
!--------------------------------------------------------------------------
      ELSEIF (imod .eq. 3) THEN
        zeta = ustar*dsqrt( rhos/phl )
        sigma_m = 12.d0*zeta*zeta*(1.d0+14.d0*zeta)
! 
        ftotal = 0.0
! 
        do j = 1, idst
           a1 = cys*dpsdm(j)
           a2 = (1.+sigma_m)
           a3 = qtotal*g/ustar**2
           fff(j) = a1*a2*a3
           fff(j) = (1-cf)*fff(j)
           ftotal = ftotal + fff(j)
        enddo
!           
        do n=1,nbins
         j0=1
         if(n.gt.1) j0=ibin(n-1)+1
         fbin(n)=0
         do j=j0,ibin(n)
          fbin(n)=fbin(n)+fff(j)
         enddo
        enddo
! 
      ENDIF  ! dust scheme
    ENDIF    ! ust < ust_t

102   continue
    do n = 1, nbins
    cell_fbin(n) = cell_fbin(n) + (soilc(cc)/tot_soilc)*fbin(n)
    enddo
    cell_ftotal = cell_ftotal + ftotal*soilc(cc)/tot_soilc

103   continue

    ENDDO  ! cc, soil category


      do n = 1, nbins
! fbin : [kg/m2/s],  dz_lowest : [m],  rho : [kg/m3],  dt : [s] -> tc :
! [kg/kg-dryair]
       tc(n) = tc(n) + cell_fbin(n)/dz_lowest/rho*dt     ![kg/kg-dryair]
       bems(n) = cell_fbin(n)     ![kg/m2/s]
      enddo
      if(shao_v .eq.1)then
       tc(5) = tc(5) + cell_ftotal/dz_lowest/rho*dt    ![kg/kg-dryair]
       bems(5) = cell_ftotal         ![kg/m2/s]
      endif
    END subroutine qf03_driver

!*****************************************************************************
    subroutine ustart0(dum, sigma, g, rho, ustart0_out)
!
! Y. Shao, 13 June 2000
!
! Calculate ustar0(d) using Shao and Lu (2000) for uncovered
! dry surface
!
! dum:    particle diameter                     [um]
! ustar0: threshold friction velocity           [m/s]
!
    real, intent(in)     :: g, rho
    real(8), intent(in)  :: dum, sigma
    real(8), intent(out) :: ustart0_out
    real(8) :: dm
    real(8), parameter :: gamma = 1.65d-4      ! a constant
    real(8), parameter :: f = 0.0123

    dm = dum*1d-6

    ustart0_out = f*(sigma*g*dm + gamma/(rho*dm) )
    ustart0_out = dsqrt( ustart0_out )
!      end function
    end subroutine ustart0
!*****************************************************************************
    subroutine qwhite(ust, ustar, rho, g, qwhite_out,shao_v)
!
! Yaping Shao 17-07-99!
!
! White (1979) Sand Flux Equation
! Q = c*rho*u_*^3 over g (1 - u_*t over u_*)(1 + u_*t^2/u_*^2)
! qwhite: Streamwise Sand Flux;       [kg m-1 s-1]
! c     : 2.6
! ust   : threhold friction velocity  [m/s]
! ustar : friction velocity           [m/s]
!
    real(8) :: c
    integer,intent(in) :: shao_v
    real, intent(in) :: ustar, rho, g
    real(8), intent(in) :: ust
    real(8), intent(out) :: qwhite_out
    real(8) :: a, b
!      if(shao_v.eq.1)then
!      c = 0.5d0
!      elseif(shao_v.eq.2)then
    c = 2.3d0
!      endif
    a = rho/g
!      IF (ustar.lt.ust) THEN 
!        qwhite = 0.d0
!      ELSE
!        b = ust/ustar 
!        qwhite = c*a*ustar**3.*(1.-b)*(1.+b*b)
!      ENDIF
    IF (ustar.lt.ust) THEN
      qwhite_out = 0.d0
    ELSE
      b = ust/ustar
!      if(shao_v.eq.1)then
!        qwhite_out = c*a*ustar**3.*(1.-b)*(1.+b*b)
!     elseif(shao_v.eq.2)then
      qwhite_out = c*a*ustar**3.*(1.-b)*(1.+b)**2.
!      endif
    ENDIF
    END subroutine qwhite
!*****************************************************************************
    subroutine vhlys(p, k, xm, alpha, u, d, vhlys_out)
!
! Volume removal according to Lu and Shao (1999), Equation (8)
! alpha: impact angle             [^o]
! u    : impact velocity          [m/s]
! p    : plastic pressure         [N/m^2]
! xm   : particle mass            [kg]
! d    : particle diameter        [m]
!
!
    REAL(8),intent(in) :: alpha, xm, u, d, p
    REAL(8) :: beta
    REAL(8), PARAMETER :: pi=3.1415927d0
    REAL(8) :: t1, t2, t3
    INTEGER,intent(in) :: k
    real(8),intent(out) :: vhlys_out

    beta = dsqrt( p*k*d/xm )
    t1 = u*u/(beta*beta)*( dsin(2.d0*alpha) - 4.d0*dsin(alpha)*dsin(alpha) )
    t2 = u*dsin(alpha)/beta
    t3 = 7.5d0*pi*t2**3.d0/d
    vhlys_out = d*( t1 + t3 )

    END subroutine vhlys
!*****************************************************************************
! A routine for correction of ust for soil moisture content
!
! w : volumetric soil moisture
! isl: soil texture type, ranging from 1 to 12
!
! Author: Yaping Shao, 5/05/2001
! Reference: Fecan et al. (1999), Ann. Geophysicae,17,149-157
!
! Data based on Shao and Jung, 2000, unpublished manuscript
! Data invented for sand, loamy sand, sandy loam, loam, clay loam, and clay
! Soil classes:
! 1 sand, 2 loamy sand, 3 sandy loam, 4 loam, 5 silty loam, 6 silt, 7 sandy clay
! loam, 8 clay loam, 
! 9 silty clay loam, 10 sandy clay, 11 silty clay, 12 clay
!----------------------------------------------------------------------
    subroutine h_c (w, wr, isl, h)

    real(8) :: a(12), b(12)!, thr(12) 
    real(8), intent(in)  :: w
    real(8), intent(out) :: h
    real, intent(in)     :: wr
    integer, intent(in)  :: isl
    character*100 :: msg  ! error message string

!       data thr/0.001, 0.003, 0.037, 0.049, 0.061, 0.072, 0.084, 0.095, 0.110,
!       0.126, 0.141, 0.156/
    data a /21.19, 33.03, 44.87, 17.79, 20.81, 23.83, 26.84, 29.86, 27.51,25.17, 22.82, 20.47/
    data b / 0.68,  0.71,  0.85,  0.61,  0.66,  0.71,  0.75,  0.80,  0.75,0.70,  0.64,  0.59/

    if ( w.lt.0. ) then
      write(msg, *) 'soil moisture correction (h_c): w = ', w, ' <  0'
      call wrf_error_fatal(msg)
!         stop
    endif

    if ( w.le.wr ) then
       h = 1.0
    else
       h = sqrt( 1 + a(isl)*( w-wr )**b(isl) )
    endif

    END subroutine h_c
!*****************************************************************************
    subroutine pmass(rhop, d, pmass_out)
!
!       Particle Mass   
!       rhop: particle density          [kg m^-3]
!       d   : particle size             [m]
!       
    REAL(8), PARAMETER :: pi=3.1415927d0
    REAL(8),intent(in) :: rhop, d
    real(8),intent(out) :: pmass_out

    pmass_out = (pi*rhop*d**3.d0)/6.d0

    END subroutine pmass

!*****************************************************************************
    subroutine r_c (x,r,shao_v)
!
!   Y. Shao 17-07-92
!   CORRECTION FUNCTION FOR UST(D) BASED ON Raupach et al. (1992)
!   x = frontal area index
!
!   R_C = (1 - sig m x)^{1/2} (1 + m beta x)^{1/2}
!         Note I deife R_C = u_{*tR}/u_{*tS}
!         While Raupach et al. defined 
!                      R_C = u_{*tS}/u_{*tR} and their R function is
!   R_C = (1 - sig m x)^{-1/2} (1 + m beta x)^{-1/2}
!   sig   : basal to frontal area; about 1
!   m     : parameter less than 1; about 0.5
!   beta  : a ratio of drag coef.; about 90.
!
    real(8) :: xc,beta
    integer,intent(in) :: shao_v
    real(8), intent(in) :: x
    real(8), intent(out) :: r
    real(8), parameter :: sig=1., m=0.5
    if(shao_v.eq.1)then
    beta=90.
    elseif(shao_v.eq.2)then
    beta=200.
    endif
!
    xc = 1./(sig*m)
    IF (x.ge.xc) THEN
      r   = 999.           ! Full covered surface
    ELSE
      r   = dsqrt(1.-sig*m*x)*dsqrt(1.+m*beta*x)
    ENDIF
!
    END subroutine r_c
!*************************************************************************************
!----------------------------------------------------------------------
subroutine mosaic_seasalt_emiss(                                       &
           dtstep, u10, v10, alt, dz8w, xland,                         &
           tend_emis, e_seas, emis_seas_factor,                        &
           its,ite, kts,kte, seasalt_emiss_active                      )
!
! adds seasalt emissions for mosaic aerosol species
! (i.e., seasalt emissions tendencies over time dtstep are applied 
! to the aerosol mixing ratios)
!

   USE iamas_atmchem_vars, only:  num_chem, param_first_scalar

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: its,ite, kts,kte

   REAL, INTENT(IN   ) ::    dtstep, emis_seas_factor
   INTEGER, INTENT(IN) ::    seasalt_emiss_active
! 10-m wind speed components (m/s)
   REAL,  DIMENSION( its:ite),                                         &
          INTENT(IN   ) ::   u10, v10, xland
! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( 1:num_chem, kts:kte, its:ite),                     &
         INTENT(INOUT ) ::   tend_emis
   REAL, DIMENSION( its:ite), INTENT(INOUT  ) :: e_seas
! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( kts:kte, its:ite)                       ,         &
          INTENT(IN   ) ::   alt,dz8w


! local variables
    integer i, k, l, n, l_oc, l_cl, l_na
    integer iphase, itype
    integer p1st

    real dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
    real factaa, factbb, fracna, fraccl, fracorg

    real :: ssemfact_numb( maxd_asize, maxd_atype )
    real :: ssemfact_mass( maxd_asize, maxd_atype )


    ! Factors for Fuentes et al (ACP, 2011, doi:10.5194/acp-11-2585-2011) 
    !                                              seasalt emission scheme.
    ! These are for calculating the seawater_OC content dependent factors.
    real, save :: alpha_f1(4) =   &
                 (/ 12.328, 38.077, 102.31, 281.65 /)
    real, save :: alpha_f2(4) =   &
                 (/ 2.2958, 8.0935, 25.251, 46.80 /)
    real, save :: alpha_f3(4) =   &
                 (/ 0.00452, 0.00705, 0.00080, 0.000761 /)
    real, save :: beta_f1(4) =   &
                 (/ 0.0311, -0.031633, 0.013154, -0.0017762 /)
    real, save :: beta_f2(4) =   &
                 (/ -13.916, 35.73, -9.7651, 1.1665 /)
    real, save :: beta_f3(4) =   &
                 (/ 4747.8, 12920.0, 7313.4, 6610.0 /)
    real :: nti(4), dp0gi(4)
    ! seawater OC<0.2um concentration (in uM) - first is for low activity, the 2nd for high activity
    !    High activity conc is from the average for RHaMBLe cruise in high activity regions
    real, save :: oc02um(2) = (/ 0.0, 280.0 /)        
    ! Organic fraction for seasalt emissions (by size bin).
    ! These are estimated from Figure 5 of Fuentes et al (ACP, 2011), assuming a value of 
    !        2 ug/l for Chl-a for the high activity
    real, save :: org_frac_low_activity(8) =   &
       (/ 0.05,  0.05,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 /)
    real, save :: org_frac_high_activity(8) =   &
       (/ 0.2,  0.2,  0.1,  0.01,  0.01,  0.01,  0.01,  0.01 /)
    

    p1st = PARAM_FIRST_SCALAR

!   for now just do itype=1
    itype = 1
    iphase = ai_phase

    ! if using Feuntes et al (2011) then calculate the seawater OC dependent factors
    if(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
        do i=1,4
            nti(i)   = beta_f1(i) * oc02um(seasalt_emiss_active-2)**2.0 + beta_f2(i) &
                            * oc02um(seasalt_emiss_active-2) + beta_f3(i)
            dp0gi(i) = alpha_f1(i) + alpha_f2(i) &
                            * exp(-alpha_f3(i)*oc02um(seasalt_emiss_active-2))
        end do
    end if


!   compute emissions factors for each size bin
    do n = 1, nsize_aer(itype)
!     changed the lower bound of the emission size limit to match lowest section edge, Qing.Yang@pnl.gov
        ! DL - 30/3/2012 - select emission scheme (1=original, 3&4=Feuntes et al, 2011)
        if(seasalt_emiss_active == 1)then
!             changed the lower bound of the emission size limit to match lowest section edge, Qing.Yang@pnl.gov
!            dumdlo = max( dlo_sect(n,itype), 0.1e-4 )
!            dumdhi = max( dhi_sect(n,itype), 0.1e-4 )
            dumdlo = max( dlo_sect(n,itype), 0.02e-4 )
            dumdhi = max( dhi_sect(n,itype), 0.02e-4 )
            call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
            ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype) )
        elseif(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
            call seasalt_emit_feuntes_1bin( dlo_sect(n,itype), dhi_sect(n,itype),  &
                ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype), nti, dp0gi, oc02um(seasalt_emiss_active-2) )
        endif

!   convert mass emissions factor from (g/m2/s) to (ug/m2/s)
        ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*1.0e6
        !jwfeng
        ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*emis_seas_factor
        ssemfact_numb(n,itype) = ssemfact_numb(n,itype)*emis_seas_factor
    end do


!   loop over i,j and apply seasalt emissions
    k = kts
    do i = its, ite
    
        ! Skip this point if over land. xland=1 for land and 2 for water.
        ! Also, there is no way to differentiate fresh from salt water.
        ! Currently, this assumes all water is salty.
        if( xland(i) < 1.5 ) cycle
    
        !wig: As far as I can tell, only real.exe knows the fractional breakdown
        !     of land use. So, in wrf.exe, dumoceanfrac will always be 1.
        dumoceanfrac = 1. !fraction of grid i,j that is salt water
        dumspd10 = dumoceanfrac* &
             ( (u10(i)*u10(i) + v10(i)*v10(i))**(0.5*3.41) )

!store the emissions  (ug/m2/s)
        e_seas(i)=0.0
        do n = 1, nsize_aer(itype)
            e_seas(i) = e_seas(i) + &
                          dumspd10*ssemfact_mass(n,itype)
        enddo
        ! factaa is (s*m2/kg-air)
        ! factaa*ssemfact_mass(n) is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
        ! factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!       factaa = (dt/dz8w(i,k,j))*alt(i,k,j) 
        factaa = alt(k,i)/dz8w(k,i)
        factbb = factaa * dumspd10

        if(seasalt_emiss_active == 1)then
            ! apportion seasalt mass emissions assumming that seasalt is pure nacl
            fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
            fraccl = 1.0 - fracna
            do n = 1, nsize_aer(itype)
                l_cl = lptr_cl_aer(n,itype,iphase)
                l_na = lptr_na_aer(n,itype,iphase)
                if ( (l_cl >= p1st) .and. (l_na >= p1st) ) then
                    tend_emis(l_cl,k,i) = factbb * ssemfact_mass(n,itype) * fraccl
                    tend_emis(l_na,k,i) = factbb * ssemfact_mass(n,itype) * fracna
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + &
                                                        factbb * ssemfact_numb(n,itype)
                end if
            end do ! n = 1, nsize_aer(itype)

        elseif (seasalt_emiss_active.eq.3 .or. seasalt_emiss_active.eq.4) then
            do n = 1, nsize_aer(itype)
                ! apportion seasalt mass emissions assumming that seasalt is a 
                ! simple mix of pure nacl and a single primary organic 
                if(seasalt_emiss_active.eq.3)then
                    fracorg = org_frac_low_activity(n)
                elseif(seasalt_emiss_active.eq.4)then
                    fracorg = org_frac_high_activity(n)
                endif
                fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
                fraccl = 1.0 - fracna
                fracna = (1.0-fracorg)*fracna
                fraccl = (1.0-fracorg)*fraccl
                
                ! only apply emissions if bin has both na and cl species
                l_oc   = lptr_oc_aer(n,itype,iphase)
                l_cl = lptr_cl_aer(n,itype,iphase)
                l_na = lptr_na_aer(n,itype,iphase)
                if ((l_cl >= p1st) .and. (l_na >= p1st) ) then
                    tend_emis(l_cl,k,i) = factbb * ssemfact_mass(n,itype) * fraccl
                    tend_emis(l_na,k,i) = factbb * ssemfact_mass(n,itype) * fracna
                   ! tend_emis(l_oc,k,i) = factbb * ssemfact_mass(n,itype) * fracorg

                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + &
                                                        factbb * ssemfact_numb(n,itype)
                end if
            end do ! n = 1, nsize_aer(itype)
        endif

    end do ! i = its, ite

    return

END subroutine mosaic_seasalt_emiss

!c----------------------------------------------------------------------
!c   following is from gong06b.f in
!c    /net/cirrus/files1/home/rce/oldfiles1/box/seasaltg
!c----------------------------------------------------------------------
subroutine seasalt_emitfactors_1bin( ireduce_smallr_emit,    &
           dpdrylo_cm, dpdryhi_cm,      &
           emitfact_numb, emitfact_surf, emitfact_mass )
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c    dpdrylo_cm  = lower dry diameter (cm)
!c    dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c    number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c    dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c    dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c   uses bubble emissions formula (eqn 5a) from 
!c    Gong et al. [JGR, 1997, p 3805-3818]
!c
!c   *** for rdry < rdry_star, this formula overpredicts emissions.
!c    A strictly ad hoc correction is applied to the formula,
!c    based on sea-salt size measurements of
!c    O'Dowd et al. [Atmos Environ, 1997, p 73-80]
!c
!c   *** the correction is only applied when ireduce_smallr_emit > 0
!c
    implicit none

!c   subr arguments
    integer ireduce_smallr_emit
    real dpdrylo_cm, dpdryhi_cm,                &
                emitfact_numb, emitfact_surf, emitfact_mass

!c   local variables
    integer isub_bin, nsub_bin

    real alnrdrylo
    real drydens, drydens_43pi_em12, x_4pi_em8
    real dum, dumadjust, dumb, dumexpb
    real dumsum_na, dumsum_ma, dumsum_sa
    real drwet, dlnrdry
    real df0drwet, df0dlnrdry, df0dlnrdry_star
    real relhum
    real rdry, rdrylo, rdryhi, rdryaa, rdrybb
    real rdrylowermost, rdryuppermost, rdry_star
    real rwet, rwetaa, rwetbb
    real rdry_cm, rwet_cm
    real sigmag_star
    real xmdry, xsdry

    real pi
    parameter (pi = 3.1415936536)

!c   c1-c4 are constants for seasalt hygroscopic growth parameterization
!c   in Eqn 3 and Table 2 of Gong et al. [1997]
    real c1, c2, c3, c4, onethird
    parameter (c1 = 0.7674)
    parameter (c2 = 3.079)
    parameter (c3 = 2.573e-11)
    parameter (c4 = -1.424)
    parameter (onethird = 1.0/3.0)


!c   dry particle density (g/cm3)
    drydens = 2.165
!c   factor for radius (micrometers) to mass (g)
    drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
!c   factor for radius (micrometers) to surface (cm2)
    x_4pi_em8 = 4.0*pi*1.0e-8
!c   bubble emissions formula assume 80% RH
    relhum = 0.80

!c   rdry_star = dry radius (micrometers) below which the
!c   dF0/dr emission formula is adjusted downwards
    rdry_star = 0.1
    if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
!c   sigmag_star = geometric standard deviation used for
!c   rdry < rdry_star
    sigmag_star = 1.9

!c   initialize sums
    dumsum_na = 0.0
    dumsum_sa = 0.0
    dumsum_ma = 0.0

!c   rdrylowermost, rdryuppermost = lower and upper 
!c   dry radii (micrometers) for overall integration
        rdrylowermost = dpdrylo_cm*0.5e4
        rdryuppermost = dpdryhi_cm*0.5e4

!c
!c   "section 1"
!c   integrate over rdry > rdry_star, where the dF0/dr emissions 
!c   formula is applicable
!c   (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20,
!c   and the entire integration is done here)
!c
    if (rdryuppermost .le. rdry_star) goto 2000

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

    nsub_bin = 1000

    alnrdrylo = log( rdrylo )
    dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

!c   compute rdry, rwet (micrometers) at lowest size
    rdrybb = exp( alnrdrylo )
    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

    do 1900 isub_bin = 1, nsub_bin

!c   rdry, rwet at sub_bin lower boundary are those
!c   at upper boundary of previous sub_bin
    rdryaa = rdrybb
    rwetaa = rwetbb

!c   compute rdry, rwet (micrometers) at sub_bin upper boundary
    dum = alnrdrylo + isub_bin*dlnrdry
    rdrybb = exp( dum )

    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

!c   geometric mean rdry, rwet (micrometers) for sub_bin
    rdry = sqrt(rdryaa * rdrybb)
    rwet = sqrt(rwetaa * rwetbb)
    drwet = rwetbb - rwetaa

!c   xmdry is dry mass in g
    xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
    xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumb is "B" in Gong's Eqn 5a
!c   df0drwet is "dF0/dr" in Gong's Eqn 5a
    dumb = ( 0.380 - log10(rwet) ) / 0.650
    dumexpb = exp( -dumb*dumb)
    df0drwet = 1.373 * (rwet**(-3.0)) *             &
              (1.0 + 0.057*(rwet**1.05)) *             &
              (10.0**(1.19*dumexpb))

    dumsum_na = dumsum_na + drwet*df0drwet
    dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
    dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

1900    continue


!c
!c   "section 2"
!c   integrate over rdry < rdry_star, where the dF0/dr emissions 
!c   formula is just an extrapolation and predicts too many emissions
!c
!c   1.  compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry) 
!c    at rdry_star
!c   2.  for rdry < rdry_star, assume dF0/dln(rdry) is lognormal,
!c    with the same lognormal parameters observed in 
!c    O'Dowd et al. [1997]
!c
2000    if (rdrylowermost .ge. rdry_star) goto 3000

!c   compute dF0/dln(rdry) at rdry_star
    rdryaa = 0.99*rdry_star
    rdry_cm = rdryaa*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetaa = rwet_cm*1.0e4

    rdrybb = 1.01*rdry_star
    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

    rwet = 0.5*(rwetaa + rwetbb)
    dumb = ( 0.380 - log10(rwet) ) / 0.650
    dumexpb = exp( -dumb*dumb)
    df0drwet = 1.373 * (rwet**(-3.0)) *             &
              (1.0 + 0.057*(rwet**1.05)) *             &
              (10.0**(1.19*dumexpb))

    drwet = rwetbb - rwetaa
    dlnrdry = log( rdrybb/rdryaa )
    df0dlnrdry_star = df0drwet * (drwet/dlnrdry)


!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

    nsub_bin = 1000

    alnrdrylo = log( rdrylo )
    dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

    do 2900 isub_bin = 1, nsub_bin

!c   geometric mean rdry (micrometers) for sub_bin
    dum = alnrdrylo + (isub_bin-0.5)*dlnrdry
    rdry = exp( dum )

!c   xmdry is dry mass in g
    xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
    xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumadjust is adjustment factor to reduce dF0/dr
    dum = log( rdry/rdry_star ) / log( sigmag_star )
    dumadjust = exp( -0.5*dum*dum )

    df0dlnrdry = df0dlnrdry_star * dumadjust

    dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
    dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
    dumsum_sa = dumsum_sa + dlnrdry*df0dlnrdry*xsdry

2900    continue


!c
!c  all done
!c
3000    emitfact_numb = dumsum_na
    emitfact_mass = dumsum_ma
    emitfact_surf = dumsum_sa

    return
end subroutine seasalt_emitfactors_1bin

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!c----------------------------------------------------------------------
!c   Following is an adaption of the subroutine above.
!c   It has been modified to include flux from Fuentes et al, ACP, 2010
!c   for dealing with smaller particle emissions.
!c----------------------------------------------------------------------
subroutine seasalt_emit_feuntes_1bin( &
           dpdrylo_cm, dpdryhi_cm,      &
           emitfact_numb, emitfact_surf, emitfact_mass, nti, dp0gi, oc02um  )
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c    dpdrylo_cm  = lower dry diameter (cm)
!c    dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c    number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c    dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c    dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c  Uses bubble emissions formula (eqn 5a) from 
!c    Gong et al. [JGR, 1997, p 3805-3818] for particles 
!c  with dry diameters of 0.45 um or more
!c
!c  For smaller particles we use the parameterisation of 
!c  Fuentes et al, ACP, 2010.
!c


    implicit none

    ! subr arguments
    !integer ireduce_smallr_emit
    real dpdrylo_cm, dpdryhi_cm,                &
                emitfact_numb, emitfact_surf, emitfact_mass
    real, intent(in) :: nti(4), dp0gi(4), oc02um

    ! local variables
    integer isub_bin, nsub_bin, jd, nsub_lower_bin, nsub_upper_bin

    real alnrdrylo
    real drydens, drydens_43pi_em12, x_4pi_em8, drydens_16pi_em12, x_pi_em8
    real dum, dumadjust, dumb, dumexpb
    real dumsum_na, dumsum_ma, dumsum_sa
    real drwet, dlnrdry, ddwet, ddry, dwet
    real df0drwet, df0dlnrdry, df0dlnrdry_star
    real relhum
    real rdry, rdrylo, rdryhi, rdryaa, rdrybb
    real rdrylowermost, rdryuppermost, rdry_star
    real rwet, rwetaa, rwetbb
    real rdry_cm, rwet_cm
    real sigmag_star
    real xmdry, xsdry
    
    real ddrylo, ddryhi, alogddrylo
    real ddrybb, ddry_cm, dwet_cm, dwetbb
    real ddryaa, dwetaa, dlogddry, df0dlogddry

    real pi
    parameter (pi = 3.1415936536)

    ! c1-c4 are constants for seasalt hygroscopic growth parameterization
    ! in Eqn 3 and Table 2 of Gong et al. [1997]
    real c1, c2, c3, c4, onethird
    parameter (c1 = 0.7674)
    parameter (c2 = 3.079)
    parameter (c3 = 2.573e-11)
    parameter (c4 = -1.424)
    parameter (onethird = 1.0/3.0)


    !  constants for seasalt/organic aerosol distribution parameterisation
    !    from Fuentes et al, ACP, 2010
    real, save :: width_ssinorg_f(4) = &
                 (/ 1.55, 1.7, 1.5, 1.7 /)
    real, save :: width_ssorg_f(4) = &
                 (/ 1.55, 1.9, 1.5, 1.7 /)
    real :: width_ss_f(4), frac
    ! scale_factor is a combination of:
    !   1) (Q)  sweep air flow (58.3 cm3/s)  
    !   2) (Ab) total surface area covered by bubbles (0.0146 m2)
    !   4) scaling factor of 3.84e-6 from the whitecap coverage formulation of 
    !          Monahan and O'Muircheartaigh (1980) - as the Gong et al formulation
    !          includes a whitecap coverage factor too
    ! this is for converting from dNt to dFp for Fuentes et al (ACP, 2010)
    real, parameter :: scale_factor = 58.3/(0.0146)*3.84e-6 !0.01533


    ! select which distribution width to use for Fuentes et al,
    if(oc02um.gt.0e0)then
        width_ss_f = width_ssorg_f
    else
        width_ss_f = width_ssinorg_f
    end if


    ! dry particle density (g/cm3)
    drydens = 2.165
    ! factor for radius (micrometers) to mass (g)
    drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
    ! factor for diameter (micrometers) to mass (g)
    drydens_16pi_em12 = drydens*(1.0/6.0)*pi*1.0e-12
    ! factor for radius (micrometers) to surface (cm2)
    x_4pi_em8 = 4.0*pi*1.0e-8
    ! factor for diameter (micrometers) to surface (cm2)
    x_pi_em8 = pi*1.0e-8
    ! bubble emissions formula assume 80% RH
    relhum = 0.80

    ! rdry_star = dry radius (micrometers) below which the
    ! dF0/dr emission formula from Fuentes et al, ACP, 2010
    ! will be used.
    rdry_star = 0.45 / 2.0
    !if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
    ! sigmag_star = geometric standard deviation used for
    ! rdry < rdry_star
    sigmag_star = 1.9

    ! initialize sums
    dumsum_na = 0.0
    dumsum_sa = 0.0
    dumsum_ma = 0.0

    ! rdrylowermost, rdryuppermost = lower and upper 
    ! dry radii (micrometers) for overall integration
    rdrylowermost = dpdrylo_cm*0.5e4
    rdryuppermost = dpdryhi_cm*0.5e4


    ! "section 1"
    ! integrate over rdry > rdry_star, where the dF0/dr emissions 
    ! formula from the original subroutine is applicable
    !  (so using Gong et al for all emission profile)
    if (rdrylowermost .ge. rdry_star) then

        ! rdrylo, rdryhi = lower and upper dry radii (micrometers) 
        ! for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = rdryuppermost

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        ! compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! rdry, rwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            rdryaa = rdrybb
            rwetaa = rwetbb

            ! compute rdry, rwet (micrometers) at sub_bin upper boundary
            dum = alnrdrylo + isub_bin*dlnrdry
            rdrybb = exp( dum )

            rdry_cm = rdrybb*1.0e-4
            rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                      ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
            rwetbb = rwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            rdry = sqrt(rdryaa * rdrybb)
            rwet = sqrt(rwetaa * rwetbb)
            drwet = rwetbb - rwetaa

            ! xmdry is dry mass in g
            xmdry = drydens_43pi_em12 * (rdry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_4pi_em8 * (rdry**2.0)

            ! dumb is "B" in Gong's Eqn 5a
            ! df0drwet is "dF0/dr" in Gong's Eqn 5a
            dumb = ( 0.380 - log10(rwet) ) / 0.650
            dumexpb = exp( -dumb*dumb)
            df0drwet = 1.373 * (rwet**(-3.0)) *             &
                      (1.0 + 0.057*(rwet**1.05)) *             &
                      (10.0**(1.19*dumexpb))

            dumsum_na = dumsum_na + drwet*df0drwet
            dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
            dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

        end do


    ! "section 2"
    ! integrate over rdry < rdry_star for old parameterisation, and use
    ! Fuentes et al below that
    !
    ! 1.  for rdry < rdry_star, use Fuentes et al parameterisation
    ! 2.  for rdry > rdry_star, use Gong et al parameterisation, as above
    elseif (rdryuppermost .gt. rdry_star) then

        ! determine the fraction of size bin below rdry_star
        frac = (log(rdry_star)-log(rdrylowermost)) / (log(rdryuppermost)-log(rdrylowermost))
        nsub_lower_bin = floor(frac*1000.0)     ! calc number of size bins to use for section below rdry_star
        nsub_upper_bin = 1000-nsub_lower_bin ! use remaining size bins for section above rdry_star


    ! 1.................
        ! ddrylo, ddryhi = lower and upper dry diameter (micrometers) 
        ! for this part of the integration
        ddrylo = rdrylowermost*2.0
        ddryhi = rdry_star*2.0

        nsub_bin = nsub_lower_bin

        alogddrylo = log10( ddrylo )
        dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

        ! compute ddry, dwet (micrometers) at lowest size
        ddrybb = 10.0**( alogddrylo )
        ddry_cm = ddrybb*1.0e-4
        dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
                  ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
        dwetbb = dwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! ddry, dwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            ddryaa = ddrybb
            !dwetaa = dwetbb

            ! compute ddry, dwet (micrometers) at sub_bin upper boundary
            dum = alogddrylo + isub_bin*dlogddry
            ddrybb = 10.0**( dum )

            ddry_cm = ddrybb*1.0e-4
            !dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
              !        ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
            !dwetbb = dwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            ddry = sqrt(ddryaa * ddrybb)
            !dwet = sqrt(dwetaa * dwetbb)

            ! xmdry is dry mass in g
            xmdry = drydens_16pi_em12 * (ddry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_pi_em8 * (ddry**2.0)

            ! Equation 2 from Fuentes et al (ACP, 2011). Wet diameter is scaled by a factor
            ! of 1e3 to convert from um to nm for this calculation
            df0dlogddry = 0.0
            do jd = 1, 4
                df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
                            * exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
            end do

            dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
            dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
            dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

        end do

    
    ! 2................
        ! rdrylo, rdryhi = lower and upper dry radii (micrometers) 
        ! for this part of the integration
        rdrylo = rdry_star
        rdryhi = rdryuppermost

        nsub_bin = nsub_upper_bin

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        ! compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! rdry, rwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            rdryaa = rdrybb
            rwetaa = rwetbb

            ! compute rdry, rwet (micrometers) at sub_bin upper boundary
            dum = alnrdrylo + isub_bin*dlnrdry
            rdrybb = exp( dum )

            rdry_cm = rdrybb*1.0e-4
            rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                      ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
            rwetbb = rwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            rdry = sqrt(rdryaa * rdrybb)
            rwet = sqrt(rwetaa * rwetbb)
            drwet = rwetbb - rwetaa

            ! xmdry is dry mass in g
            xmdry = drydens_43pi_em12 * (rdry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_4pi_em8 * (rdry**2.0)

            ! dumb is "B" in Gong's Eqn 5a
            ! df0drwet is "dF0/dr" in Gong's Eqn 5a
            dumb = ( 0.380 - log10(rwet) ) / 0.650
            dumexpb = exp( -dumb*dumb)
            df0drwet = 1.373 * (rwet**(-3.0)) *             &
                      (1.0 + 0.057*(rwet**1.05)) *             &
                      (10.0**(1.19*dumexpb))

            dumsum_na = dumsum_na + drwet*df0drwet
            dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
            dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

        end do



    ! "section 3"
    ! where rdry < rdry_star for the whole bin we use Fuentes et al 
    ! for all emissions
    !
    else
    
        ! ddrylo, ddryhi = lower and upper dry diameter (micrometers) 
        ! for this part of the integration
        ddrylo = rdrylowermost*2.0
        ddryhi = rdryuppermost*2.0

        nsub_bin = 1000

        alogddrylo = log10( ddrylo )
        dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

        ! compute ddry, dwet (micrometers) at lowest size
        ddrybb = 10.0**( alogddrylo )
        ddry_cm = ddrybb*1.0e-4
        dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
                  ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
        dwetbb = dwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! ddry, dwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            ddryaa = ddrybb
            !dwetaa = dwetbb

            ! compute ddry, dwet (micrometers) at sub_bin upper boundary
            dum = alogddrylo + isub_bin*dlogddry
            ddrybb = 10.0**( dum )

            ddry_cm = ddrybb*1.0e-4
            !dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
              !        ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
            !dwetbb = dwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            ddry = sqrt(ddryaa * ddrybb)
            !dwet = sqrt(dwetaa * dwetbb)

            ! xmdry is dry mass in g
            xmdry = drydens_16pi_em12 * (ddry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_pi_em8 * (ddry**2.0)

            ! Equation 2 from Fuentes et al (ACP, 2011). Wet diameter is scaled by a factor
            ! of 1e3 to convert from um to nm for this calculation
            df0dlogddry = 0.0
            do jd = 1, 4
                df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
                            * exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
            end do

            dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
            dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
            dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

        end do
    


    end if
    !c
    !c  all done
    !c    
    emitfact_numb = dumsum_na
    emitfact_mass = dumsum_ma
    emitfact_surf = dumsum_sa


    return
end subroutine seasalt_emit_feuntes_1bin

subroutine mosaic_seasalt_emiss2(                            &
           dtstep, u10, v10,u_phy,v_phy,alt,rho_phy,dz8w,    &
           xland,dx,g,tend_emis,sstemp,e_seas,               &
           emis_seas_factor,do_seas_sst_adjust,              &
           its,ite, kts,kte,seasalt_emiss_active    )

! Routine to calculate seasalt emissions for MOSAIC over the time dtstep... 
! Chun.Zhao@pnnl.gov; 10-Aug-2011
! New formula including correction of particle radius < 0.2 um and SST correction
! Reference:
! Gong 2003 and Jeagle et al. 2011
!------------------------------------------------------------------------

    use iamas_atmchem_vars, only:  num_chem, param_first_scalar


    implicit none

    integer, intent(in   ) :: its,ite, kts,kte
    real, intent(in   ) :: dtstep, g, emis_seas_factor
    logical, intent(in   ) :: do_seas_sst_adjust
    real, dimension(its:ite), intent(inout) :: e_seas
    real, dimension(its:ite), intent(in   ) :: sstemp
    real, dimension(its:ite), intent(in   ) :: dx
    real, dimension(its:ite), intent(in   ) :: u10       , &
                                               v10       , &
                                               xland
    real, dimension(kts:kte,its:ite), intent(in   ) :: alt, u_phy, v_phy, rho_phy, dz8w
    real, dimension(1:num_chem, kts:kte, its:ite), intent(inout) :: tend_emis

! local variables
    integer :: i,k,n,l
    real dumdlo,dumdhi
    !jwfeng: There is no need to use 1000 sub-bins, use 5 instead
    integer, parameter :: nsub_bin = 5 ! sub-bin for integration
    real bems_mass(nsub_bin),bems_num(nsub_bin)  ! sea-salt emission for bins 
    real total_mass,total_num
    real w10m,gwet,airden,airmas
    real*8, DIMENSION( its:ite) :: dxy
    real fracna,fraccl,factaa,factbb,dumspd10,fracorg,dumoceanfrac
    real conver,converi
    real ssemfact_numb( maxd_asize, maxd_atype )
    real ssemfact_mass( maxd_asize, maxd_atype )
    real sst_c
    integer l_cl,l_na,iphase,itype,p1st,l_oc
    character(len=strkind) :: debugString
    INTEGER, INTENT(IN) ::    seasalt_emiss_active
    real, save :: org_frac_low_activity(8) =   &
       (/ 0.05,  0.05,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 /)
    real, save :: org_frac_high_activity(8) =   &
       (/ 0.2,  0.2,  0.1,  0.01,  0.01,  0.01,  0.01,  0.01 /)

    p1st    = PARAM_FIRST_SCALAR
    itype   = 1
    iphase  = ai_phase
    conver  = 1.d-9
    converi = 1.d9

    k=kts
    do i=its,ite
        e_seas(i)=0.0
        !sea-salt only over water
        if ( xland(i).gt.1.5 ) then
        !-----------------------------------------------------
        !1. Pre-setting inputs 
        !-----------------------------------------------------
            !1.1 10-m wind speed
            w10m = sqrt(u10(i)*u10(i)+v10(i)*v10(i))
            ! we don't trust the u10,v10 values, for thin model surface layer
            if (dz8w(kts,i).le.10.) then
                w10m = sqrt(u_phy(kts,i)*u_phy(kts,i)+v_phy(kts,i)*v_phy(kts,i))
            end if
            !1.2 Grid area 
            dxy = dx(i)*dx(i)
            !1.3 Sea-surface temperature in degree C
            sst_c=sstemp(i)-273.15
        !-----------------------------------------------------
        !2. Emissions 
        !-----------------------------------------------------
            !2.1 Emission for bins 
            do n = 1, nsize_aer(itype)
                dumdlo = max(dlo_sect(n,itype), 0.02e-4 )
                dumdhi = max(dhi_sect(n,itype), 0.02e-4 )
                dumdlo = dumdlo * 1.e4 ! cm to um
                dumdhi = dumdhi * 1.e4 ! cm to um 
                call seasalt_source(dtstep,w10m,sst_c,dumdlo,dumdhi,nsub_bin,bems_mass,bems_num,do_seas_sst_adjust)
                !bems_mass: particles kg m-2 s-1
                !bems_num: particles # m-2 s-1
                ssemfact_mass(n,itype) = sum(bems_mass(1:nsub_bin))*converi !kg to ug
                ssemfact_numb(n,itype) = sum(bems_num(1:nsub_bin))
                !jwfeng: add seasalt emis factor
                ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*emis_seas_factor
                ssemfact_numb(n,itype) = ssemfact_numb(n,itype)*emis_seas_factor
            enddo

            !2.2 Sum the total emission
            total_mass = 0.0
            total_num  = 0.0
            do n = 1, nsize_aer(itype)
                total_mass = total_mass + ssemfact_mass(n,itype)
                total_num  = total_num  + ssemfact_numb(n,itype)
            enddo
            !ug m-2 s-1
            e_seas(i) = total_mass

            !2.3 Apportion seasalt mass emissions assumming that seasalt is pure NaCl
!            fracseas = 1.0

         !-----------------------------------------------------
         !3. Update chemistry tendency (Accumulation and Coarse modes) 
         !-----------------------------------------------------
            !3.1 Factor converting emission to concentration
            !factaa is (s*m2/kg-air)
            !factaa*ssemfact_mass is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
            !factaa*ssemfact_numb is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
            dumoceanfrac = 1. !fraction of grid i,j that is salt water
            dumspd10 = dumoceanfrac* &
                 ( (u10(i)*u10(i) + v10(i)*v10(i))**(0.5*3.41) )
! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!         factaa = (dtstep/dz8w(k,i))*alt(k,i) 
            factaa = alt(k,i)/dz8w(k,i)
            factbb = factaa * dumspd10

        if(seasalt_emiss_active == 1)then
            ! apportion seasalt mass emissions assumming that seasalt is pure nacl
            fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
            fraccl = 1.0 - fracna
            do n = 1, nsize_aer(itype)
                l_cl = lptr_cl_aer(n,itype,iphase)
                l_na = lptr_na_aer(n,itype,iphase)
                if ( (l_cl >= p1st) .and. (l_na >= p1st) ) then
                    tend_emis(l_cl,k,i) = factbb * ssemfact_mass(n,itype) * fraccl
                    tend_emis(l_na,k,i) = factbb * ssemfact_mass(n,itype) * fracna
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) + &
                                                        factbb * ssemfact_numb(n,itype)
                end if
            end do ! n = 1, nsize_aer(itype) 
            !3.2 Update 
        elseif (seasalt_emiss_active.eq.3 .or. seasalt_emiss_active.eq.4) then
            do n = 1, nsize_aer(itype)

                if(seasalt_emiss_active.eq.3)then
                    fracorg = org_frac_low_activity(n)
                elseif(seasalt_emiss_active.eq.4)then
                    fracorg = org_frac_high_activity(n)
                endif

                fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
                fraccl = 1.0 - fracna
                fracna = (1.0-fracorg)*fracna
                fraccl = (1.0-fracorg)*fraccl

                l_cl = lptr_cl_aer(n,itype,iphase)
                l_na = lptr_na_aer(n,itype,iphase)
                l_oc   = lptr_oc_aer(n,itype,iphase)
                if ( l_cl >= p1st .and. l_na >= p1st ) then
                    tend_emis(l_cl,k,i) = factaa * ssemfact_mass(n,itype) * fraccl
                    tend_emis(l_na,k,i) = factaa * ssemfact_mass(n,itype) * fracna
                   ! tend_emis(l_oc,k,i) = factaa * ssemfact_mass(n,itype) * fracorg
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(l,k,i) = tend_emis(l,k,i) +   &
                                                        factaa * ssemfact_numb(n,itype)
                endif
            enddo
        endif
        endif !xland 

    enddo

end subroutine mosaic_seasalt_emiss2

subroutine seasalt_source(dt,w10m,sst_c,diamlo,diamhi,nsub_bin,bems_mass,bems_num,do_seas_sst_adjust)
! ****************************************************************************
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) The wet radius r is estimated at 80% RH
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr_w/dr_d) * dFn/dr_w
! *    or 
! *    dFM = 4/3*pi*rho_d*r_d^3 * dFn
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
  IMPLICIT NONE

  !input & output
  REAL,    INTENT(IN)    :: dt
  REAL,    INTENT(IN)    :: w10m
  REAL,    INTENT(IN)    :: sst_c
  REAL,    INTENT(IN)    :: diamlo,diamhi
  INTEGER, INTENT(IN)    :: nsub_bin
  LOGICAL, INTENT(IN)    :: do_seas_sst_adjust
  REAL,    INTENT(INOUT) :: bems_mass(nsub_bin),bems_num(nsub_bin)

  !local variables
  !values used in Gong's paper (2003)
  REAL, PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/)
  REAL, PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
  !Change suggested by MC
  !REAL, PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/)
  !REAL, PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
  REAL, PARAMETER :: b_old(2)=(/0.380, 0.650/)
  REAL, PARAMETER :: b_new(2)=(/0.433, 0.433/)
  REAL :: c0(5), b0(2)
 !REAL*8, PARAMETER :: dr=5.0D-2 ! um, radius interval for integration 
  REAL, PARAMETER :: theta=30.0
  REAL, PARAMETER :: den_seas=2.2e3 !sea-salt density, kg/m3
  REAL :: r_to_mass,xmdry
 !Swelling coefficient frh (d rwet / d rd)
 !REAL*8,    PARAMETER :: frh = 1.65
 !REAL*8,    PARAMETER :: frh = 2.d0
  LOGICAL, PARAMETER :: new=.TRUE.  ! new from Gong (2003) or old (1997)
  REAL SSTFACTOR
  REAL :: rho_d, rdrylo, rdryhi
  real :: rd1, rd2, rw1,rw2, rd1_cm,rd2_cm,rw1_cm,rw2_cm
  real :: rdry,rwet,drdry,drwet
  real,parameter :: c1 = 0.7674
  real,parameter :: c2 = 3.079
  real,parameter :: c3 = 2.573e-11
  real,parameter :: c4 = -1.424
  real,parameter :: onethird = 1.0/3.0
  real,parameter :: pi = 3.1415936536
  real,parameter :: relhum = 0.8
  real :: a, b, dfndr,dfn,dfm
  integer isub_bin

 !-----------------------------------------------------
 !1. Initializing
 !-----------------------------------------------------
     bems_mass = 0.0
     bems_num = 0.0
     rho_d = den_seas
     !factor for radius (micrometers) to mass (kg)
     r_to_mass = rho_d*(4.0/3.0)*pi*1.0e-18
     rdrylo = diamlo/2.0
     rdryhi = diamhi/2.0
     drdry = (rdryhi - rdrylo)/nsub_bin

 !-----------------------------------------------------
 !2. Emission source 
 !-----------------------------------------------------
     !2.1 compute rdry,rwet (micrometers) at lowest bound 
     rd2 = rdrylo
     rd2_cm = rdrylo*1.e-4
     rw2_cm = ( rd2_cm**3 +  &
               (c1*(rd2_cm**c2))/((c3*(rd2_cm**c4))-log10(relhum)) &
              )**onethird
     rw2 = rw2_cm*1.0e4

     do 2000 isub_bin = 1, nsub_bin

       !2.2 rdry&rwet at the sub_bin lower and upper bound 
       rd1 = rd2
       rw1 = rw2
       rd2 = rdrylo + isub_bin*drdry
       rd2_cm = rd2*1.e-4
       rw2_cm = ( rd2_cm**3 +  &
                 (c1*(rd2_cm**c2))/((c3*(rd2_cm**c4))-log10(relhum)) &
                )**onethird
       rw2 = rw2_cm*1.0e4

       !2.3 geometric mean rdry&rwet (micrometers) for sub_bin
       rdry = 0.5*(rd1+rd2)
       rwet = 0.5*(rw1+rw2)
       drwet= rw2 - rw1

       !2.4 convert factor (xmdry),from wet number to dry mass in kg
       !xmdry = r_to_mass * (drwet/drdry) * (rdry**3.0)
       xmdry = r_to_mass * (rdry**3.0)

       !2.5 Coefficients for source formula 
       IF (new) THEN
          a = -1.0*(4.7*(1.0 + theta*rwet)**(-0.017*rwet**(-1.44)))
          c0 = c_new
          b0 = b_new
       ELSE
          a = -3.0
          c0 = c_old
          b0 = b_old
       END IF
       b = (b0(1) - LOG10(rwet))/b0(2)

       !2.6 Source function 
       !particles # m-2 s-1 um-1
       dfndr = (c0(1)*(rwet**a))*(w10m**c0(2))*(1.0+c0(3)*rwet**c0(4))* &
              10**(c0(5)*EXP(-(b**2)))
       if (do_seas_sst_adjust) then ! SST adjustment 
        if(sst_c.ge.0.and.sst_c.le.30) then
        sstfactor = (0.3+0.1*sst_c-0.0076*sst_c**2+0.00021*sst_c**3)
        dfndr = dfndr*sstfactor
        endif
       endif

       !particles # m-2 s-1
       dfn = dfndr*drwet
       !particles kg m-2 s-1
       dfm = xmdry*dfn
       bems_mass(isub_bin) = dfm
       bems_num(isub_bin) = dfn

 2000  CONTINUE

       return

END SUBROUTINE seasalt_source


end module module_mosaic_addemiss
