!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute natural, anthropogenic and biomass burning emission. 
! * Author: Zihan Xia, Jiawang Feng, Chen Jin and Chun Zhao from AEMOL USTC
! * Email: xiazh@mail.ustc.edu.cn
! Last update: zhxia 04/2025
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

 module iamas_atmchem_driver_emission
 use iamas_timer, only: iamas_timer_start, iamas_timer_stop
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_constants
 use iamas_atmphys_vars
 use iamas_atmchem_configs
 use iamas_atmchem_interface
 use iamas_atmchem_utilities
!chemistry schemes:
 use module_mosaic_addemiss
 use module_data_binaer_asect, only : lptr_dust_aer, lptr_cl_aer, lptr_na_aer

 implicit none
 private
 public:: driver_emission

!iAMAS driver for chemistry emission.

 contains

!=================================================================================================================
 subroutine driver_emission(itimestep,configs,mesh,state,emis_state,sfc_input,diag,diag_physics,diag_chemistry,tend_chemistry,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: state
 type(iamas_pool_type),intent(in):: emis_state
 type(iamas_pool_type),intent(in):: sfc_input
 type(iamas_pool_type),intent(in):: diag
 type(iamas_pool_type),intent(in):: diag_physics

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep
!inout arguments:
 type(iamas_pool_type),intent(inout):: tend_chemistry
 type(iamas_pool_type),intent(inout):: diag_chemistry
!local variables:
 character(len=StrKIND),pointer:: config_chem_opt, config_dust_opt
 integer,pointer:: config_seas_opt
 logical,pointer:: config_emis_update, config_use_anthro_emis, config_use_bioms_emis
 real,pointer:: config_dust_factor, config_seas_factor
 logical,pointer:: do_seas_sst_adjust
 integer:: i,k,j
 character(len=StrKIND):: debugString
 real(kind=RKIND),dimension(:),allocatable :: slai
 real(kind=RKIND),dimension(:,:),allocatable :: landusef_temp !zhxia add this for avoiding temp array warning. 04/2025
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_emis
 real(kind=RKIND),dimension(:),pointer    :: sst,e_dust,e_seas
 real(kind=RKIND),dimension(:,:),pointer  :: smois,erod
 real(kind=RKIND),dimension(its:ite),target :: not_used_pointer

! Note, the shape and order of species array must be same between tend_emis and chem (belong to different var struct) because they use the same index of species and have matrix addition calculation. @zhxia
 real(kind=RKIND), dimension(:), pointer :: e_nh3    ,e_co    ,e_no      ,e_no2       ,e_so2        , &
                                            e_bc     ,e_oc    ,e_pm25    ,e_pm10      ,e_so4        , &
                                            e_no3    ,e_c2h6  ,e_c3h8    ,e_c2h2      ,e_alk3       , &
                                            e_alk4   ,e_alk5  ,e_ethene  ,e_c3h6      ,e_ole1       , &
                                            e_ole2   ,e_aro1  ,e_aro2    ,e_hcho      ,e_ccho       , &
                                            e_rcho   ,e_acet  ,e_mek     ,e_isoprene  ,e_terp       , &
                                            e_sesq   ,e_phen  ,e_cres    ,e_meoh      ,e_gly        , &
                                            e_mgly   ,e_bacl  ,e_isoprod ,e_methacro  ,e_mvk        , &
                                            e_prod2  ,e_ch4   ,e_bald    ,e_hcooh     ,e_cco_oh     , &
                                            e_rco_oh ,e_co2
 real(kind=rkind), dimension(:), pointer :: ebu_bc    ,ebu_oc     ,ebu_pm25  ,ebu_pm10  ,ebu_ch4   , &
                                            ebu_co2   ,ebu_co     ,ebu_no    ,ebu_no2   ,ebu_so2   , &
                                            ebu_acet  ,ebu_alk3   ,ebu_alk4  ,ebu_alk5  ,ebu_aro1  , &
                                            ebu_aro2  ,ebu_bald   ,ebu_ccho  ,ebu_c2o2h ,ebu_ethen , &
                                            ebu_hcho  ,ebu_hcooh  ,ebu_hono  ,ebu_isopr ,ebu_mek   , &
                                            ebu_meoh  ,ebu_metha  ,ebu_mgly  ,ebu_mvk   ,ebu_ole1  , &
                                            ebu_ole2  ,ebu_phen   ,ebu_prod2 ,ebu_rcho

 call iamas_pool_get_config(configs      ,'config_chem_opt'        ,config_chem_opt)
 call iamas_pool_get_config(configs      ,'config_dust_opt'        ,config_dust_opt) 
 call iamas_pool_get_config(configs      ,'config_seas_opt'        ,config_seas_opt)
 call iamas_pool_get_config(configs      ,'config_emis_update'     ,config_emis_update)
 call iamas_pool_get_config(configs      ,'config_use_anthro_emis' ,config_use_anthro_emis)
 call iamas_pool_get_config(configs      ,'config_use_bioms_emis'  ,config_use_bioms_emis)
 call iamas_pool_get_config(configs      ,'config_dust_factor'     ,config_dust_factor)
 call iamas_pool_get_config(configs      ,'config_seas_factor'     ,config_seas_factor)
 call iamas_pool_get_config(configs      ,'do_seas_sst_adjust'     ,do_seas_sst_adjust)

 if(.not.allocated(       slai) )    allocate(         slai(its:ite) )
 if(.not.allocated(       landusef_temp) )    allocate(         landusef_temp(1:24,its:ite) )
 landusef_temp = landusef(1:24, its:ite)
 call iamas_pool_get_array(sfc_input     ,'smois'      ,smois)      !defined in iamas_atmchem_vars.F, don't get in iamas_atmchem_interface.F
 call iamas_pool_get_array(sfc_input,'sst'             ,sst  )
 call iamas_pool_get_array(sfc_input,'erod'            ,erod )
 call iamas_pool_get_array(tend_chemistry,'tend_emis'  ,tend_emis)
 not_used_pointer(:) = 0.0
 if ( lptr_dust_aer(1,1,1)>=param_first_scalar ) then
    call iamas_pool_get_array(diag_chemistry,'e_dust'  ,e_dust)
 else
    e_dust => not_used_pointer
 endif
 if ( lptr_cl_aer(1,1,1)>=param_first_scalar .and. lptr_na_aer(1,1,1)>=param_first_scalar ) then
    call iamas_pool_get_array(diag_chemistry,'e_seas'  ,e_seas)
 else
    e_seas => not_used_pointer
 endif

 if (config_emis_update) then
    if (trim(config_chem_opt) == 'bin04aer_gas_aq' .or.   &   
        trim(config_chem_opt) == 'bin04aer_gas') then

        !anthropogenic.    zhxia@202504
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
        !biomass burning.  zhxia@2026
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

    else
        call chemistry_error_fatal('Error: Emission update fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can set config_emis_update = True.')
    endif
 endif


 ! init set values =0
 tend_emis(:,:,:) = 0.

 ! Dust and seasalt emission.
 call iamas_timer_start('emis_natural')
 call mosaic_addemiss_nature(configs = configs                                                                                        , &
        slai = slai(its:ite)                  ,ust = ust(its:ite)                     ,smois = smois(1:num_soils, its:ite)            , &
        ivgtyp = ivgtyp(its:ite)              ,vegfra = vegfra(its:ite)               ,isltyp = isltyp(its:ite)                       , &
        dtstep = dt_dyn                       ,u10 = u10(its:ite)                     ,v10 = v10(its:ite)                             , &
        alt = alt(kts:kte, its:ite)           ,dz8w = dz(kts:kte, its:ite)            ,xland = xland(its:ite)                         , & 
        sstemp = sst(its:ite)                 ,ktau = 1                               ,num_soil_layers = num_soils                    , &
        u_phy = u_phy(kts:kte, its:ite)       ,v_phy = v_phy(kts:kte, its:ite)        ,rho_phy = rho_phy(kts:kte, its:ite)            , &
        p8w = pres2_hyd(kts:kte+1, its:ite)   ,snowc = snowc(its:ite)                 ,erod = erod(1:3, its:ite)                      , &
        xlat = xlat(its:ite)                  ,xlong = xlong(its:ite)                 ,dx = dx(its:ite)                               , &
        g = gravity                           ,emis_dust_factor = config_dust_factor  ,emis_seas_factor = config_seas_factor          , &
        tend_emis = tend_emis(1:num_chem, kts:kte, its:ite),e_dust = e_dust(its:ite)  ,e_seas = e_seas(its:ite)                       , &
        dust_emiss_active = config_dust_opt   ,seasalt_emiss_active = config_seas_opt ,do_seas_sst_adjust = do_seas_sst_adjust        , &
        soilctop = soilctop(1:16, its:ite)    ,landusef = landusef_temp(1:24,its:ite) ,sh2o = sh2o(1:4, its:ite),xlai = lai(its:ite)  , &
        its = its, ite = ite, kts = kts, kte = kte)
 call iamas_timer_stop('emis_natural')

 call chemistry_debug_message('Nature aerosol emission success')

 !Anthropogenic and biomass burning emission. @zhxia
 if (config_emis_update) then
     if (config_use_anthro_emis) then
 call iamas_timer_start('emis_anthropogenic')
         call mosaic_addemiss_anthropogenic(configs = configs                                                                                    , &
                  tend_emis = tend_emis(1:num_chem, kts:kte, its:ite), dz8w = dz(kts:kte, its:ite), alt = alt(kts:kte, its:ite)                , &
                  e_nh3 = e_nh3(its:ite)          , e_co = e_co(its:ite)          , e_no = e_no(its:ite)            , e_no2 = e_no2(its:ite)     , &
                  e_so2 = e_so2(its:ite)          , e_bc = e_bc(its:ite)          , e_oc = e_oc(its:ite)            , e_pm25 = e_pm25(its:ite)   , &
                  e_pm10 = e_pm10(its:ite)        , e_so4 = e_so4(its:ite)        , e_no3 = e_no3(its:ite)          , e_c2h6 = e_c2h6(its:ite)   , &
                  e_c3h8 = e_c3h8(its:ite)        , e_c2h2 = e_c2h2(its:ite)      , e_alk3 = e_alk3(its:ite)        , e_alk4 = e_alk4(its:ite)   , &
                  e_alk5 = e_alk5(its:ite)        , e_ethene = e_ethene(its:ite)  , e_c3h6 = e_c3h6(its:ite)        , e_ole1 = e_ole1(its:ite)   , &
                  e_ole2 = e_ole2(its:ite)        , e_aro1 = e_aro1(its:ite)      , e_aro2 = e_aro2(its:ite)        , e_hcho = e_hcho(its:ite)   , &
                  e_ccho = e_ccho(its:ite)        , e_rcho = e_rcho(its:ite)      , e_acet = e_acet(its:ite)        , e_mek = e_mek(its:ite)     , &
                  e_isoprene = e_isoprene(its:ite), e_terp = e_terp(its:ite)      , e_sesq = e_sesq(its:ite)        , e_phen = e_phen(its:ite)   , &
                  e_cres = e_cres(its:ite)        , e_meoh = e_meoh(its:ite)      , e_gly = e_gly(its:ite)          , e_mgly = e_mgly(its:ite)   , &
                  e_bacl = e_bacl(its:ite)        , e_isoprod = e_isoprod(its:ite), e_methacro = e_methacro(its:ite), e_mvk = e_mvk(its:ite)     , &
                  e_prod2 = e_prod2(its:ite)      , e_ch4 = e_ch4(its:ite)        , e_bald = e_bald(its:ite)        , e_hcooh = e_hcooh(its:ite) , &
                  e_cco_oh = e_cco_oh(its:ite)    , e_rco_oh = e_rco_oh(its:ite)  , e_co2 = e_co2(its:ite)                                       , &
                  its = its, ite = ite, kts = kts, kte = kte)
 call iamas_timer_stop('emis_anthropogenic')

         call chemistry_debug_message('Anthropogenic emission success')
     endif
     if (config_use_bioms_emis)  then
 call iamas_timer_start('emis_biomass')
         call mosaic_addemiss_bioms(configs = configs                                                                                , &
                  tend_emis = tend_emis(1:num_chem, kts:kte, its:ite), dz8w = dz(kts:kte, its:ite), alt = alt(kts:kte, its:ite)      , &
                  ebu_bc   = ebu_bc(its:ite)        , ebu_oc   = ebu_oc(its:ite)       , ebu_pm25 = ebu_pm25(its:ite)                , &
                  ebu_pm10 = ebu_pm10(its:ite)      , ebu_ch4  = ebu_ch4(its:ite)      , ebu_co2  = ebu_co2(its:ite)                 , &
                  ebu_co   = ebu_co(its:ite)        , ebu_no   = ebu_no(its:ite)       , ebu_no2  = ebu_no2(its:ite)                 , &
                  ebu_so2  = ebu_so2(its:ite)       , ebu_acet = ebu_acet(its:ite)     , ebu_alk3 = ebu_alk3(its:ite)                , &
                  ebu_alk4 = ebu_alk4(its:ite)      , ebu_alk5 = ebu_alk5(its:ite)     , ebu_aro1 = ebu_aro1(its:ite)                , &
                  ebu_aro2 = ebu_aro2(its:ite)      , ebu_bald = ebu_bald(its:ite)     , ebu_ccho = ebu_ccho(its:ite)                , &
                  ebu_c2o2h= ebu_c2o2h(its:ite)     , ebu_ethen= ebu_ethen(its:ite)    , ebu_hcho = ebu_hcho(its:ite)                , &
                  ebu_hcooh= ebu_hcooh(its:ite)     , ebu_hono = ebu_hono(its:ite)     , ebu_isopr= ebu_isopr(its:ite)               , &
                  ebu_mek  = ebu_mek(its:ite)       , ebu_meoh = ebu_meoh(its:ite)     , ebu_metha= ebu_metha(its:ite)               , &
                  ebu_mgly = ebu_mgly(its:ite)      , ebu_mvk  = ebu_mvk(its:ite)      , ebu_ole1 = ebu_ole1(its:ite)                , &
                  ebu_ole2 = ebu_ole2(its:ite)      , ebu_phen = ebu_phen(its:ite)     , ebu_prod2= ebu_prod2(its:ite)               , &
                  ebu_rcho = ebu_rcho(its:ite)                                                                                       , &
                  its = its, ite = ite, kts = kts, kte = kte)
 call iamas_timer_stop('emis_biomass')
         call chemistry_debug_message('Biomass burning emission success')
     endif
 endif

 end subroutine driver_emission

!=================================================================================================================
 end module iamas_atmchem_driver_emission
!=================================================================================================================
