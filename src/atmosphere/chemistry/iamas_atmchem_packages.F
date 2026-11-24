!************************************************************************
!
! Module to Set Up Definitions of All Chemistry Packages Through Configs.
! * Author: Jiawang Feng, Zihan Xia and Chun Zhao / AEMOL USTC
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
!
! References: 
!   iamas_atmphys_packages.F by Laura D. Fowler
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_packages
 use iamas_kind_types
 use iamas_derived_types,only : iamas_pool_type,iamas_io_context_type,iAMAS_LOG_ERR
 use iamas_pool_routines,only : iamas_pool_get_config,iamas_pool_get_package
!! use iamas_io_units,only : stderrUnit
 use iamas_log,only : iamas_log_write


 implicit none
 private
 public:: atmchem_setup_packages

 contains


!=================================================================================================================
 function atmchem_setup_packages(configs,packages,iocontext) result(ierr)
!=================================================================================================================

!===========================================================================
!IN THIS FUNCTION, IT INCLUDES FOLLOWING PARTS:
!  + Definitions of all chemistry packages. 
!                                         ---jwfeng   12/2019
!NOTE:
!    If you want to add a package,you need to add a "package" in
!    Registry.xml and add a new config_chem_opt relating to your
!    new package.Don't forget to update package under chem vars.
!===========================================================================

!inout arguments:
type (iamas_pool_type), intent(inout) :: configs
type (iamas_pool_type), intent(inout) :: packages
type (iamas_io_context_type), intent(inout) :: iocontext

!local variables:
character(len=StrKIND),pointer:: config_chem_opt
character(len=StrKIND),pointer:: config_phot_opt !zhxia add for photolysis
character(len=StrKIND),pointer:: config_gaschem_opt !zhxia add for photochemistry
character(len=StrKIND),pointer:: config_aerchem_sia_opt !xyruan add for inorganic aerosol chemistry and interactions
character(len=StrKIND),pointer:: config_aerchem_soa_opt !xyruan add for organic aerosol chemistry and interactions
character(len=StrKIND),pointer:: config_stra_gaschem_opt !zhxia add for stratospheric photochemistry
character(len=StrKIND),pointer:: config_cloudchem_opt !xyruan add for cloud aqueous chemistry
logical,pointer:: no_chemistry, notIncluded
logical,pointer:: bin10aer_nature_aq, bin09aer_nature_aq, bin08aer_nature_aq, bin05aer_nature_aq, bin04aer_nature_aq
logical,pointer:: bin10aer_nature   , bin09aer_nature   , bin08aer_nature   , bin05aer_nature   , bin04aer_nature   
logical,pointer:: bin10dust_aq, bin09dust_aq, bin08dust_aq, bin05dust_aq, bin04dust_aq
logical,pointer:: bin10dust   , bin09dust   , bin08dust   , bin05dust   , bin04dust   
logical,pointer:: bin10seas_aq, bin09seas_aq, bin08seas_aq, bin05seas_aq, bin04seas_aq
logical,pointer:: bin10seas   , bin09seas   , bin08seas   , bin05seas   , bin04seas   
logical,pointer:: bin04aer_gas_aq, bin04aer_gas
logical,pointer:: no_phot, photfastj !zhxia
logical,pointer:: no_gaschem, AIPC_saprc99, KPP_saprc99 !zhxia
logical,pointer:: no_stra_gaschem, linoz !zhxia
logical,pointer:: no_cloudchem, CMU_cloudchem !xyruan
logical,pointer:: no_sia, AIMACI_mosaic, simple_SIA, MOSAIC_mech !zhxia,xyruan
logical,pointer:: no_soa, simple_SOA !xyruan
integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
!Debug Information,Not Necessary.
call iamas_log_write( ' ')
call iamas_log_write( '--- enter subroutine atmchem_setup_packages:')

ierr = 0

call iamas_log_write( '----- Setting up package variables -----')
call iamas_log_write( ' ')

!--- zhxia:initialization of all packages for photolysis shemes:

call iamas_pool_get_config(configs,'config_phot_opt',config_phot_opt)
nullify(no_phot)
nullify(photfastj)
call iamas_pool_get_package(packages, 'no_photActive' ,  no_phot )
call iamas_pool_get_package(packages, 'photfastjActive' ,  photfastj )

if(.not.associated(no_phot) .or. &
   .not.associated(photfastj)) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for photolysis schemes options(config_phot_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_phot = .false.
photfastj = .false.

PHOT_OPT:select case (config_phot_opt)
    case ('no_phot')
        no_phot = .true.
    case ('photfastj')
        photfastj = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for photolysis schemes options(config_phot_opt). *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select PHOT_OPT

call iamas_log_write( ' no_phot  = $l', logicArgs=(/no_phot/) )
call iamas_log_write( ' photfastj  = $l', logicArgs=(/photfastj/) )

!--- zhxia:initialization of all packages for photochemistry shemes:

call iamas_pool_get_config(configs,'config_gaschem_opt',config_gaschem_opt)
nullify(no_gaschem)
nullify(AIPC_saprc99)
nullify(KPP_saprc99)
call iamas_pool_get_package(packages, 'no_gaschemActive' ,  no_gaschem )
call iamas_pool_get_package(packages, 'AIPC_saprc99Active' ,  AIPC_saprc99 )
call iamas_pool_get_package(packages, 'KPP_saprc99Active' ,  KPP_saprc99 )

if(.not.associated(no_gaschem) .or. &
   .not.associated(AIPC_saprc99) .or. &
   .not.associated(KPP_saprc99)) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for photochemistry schemes options(config_gaschem_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_gaschem = .false.
AIPC_saprc99 = .false.
KPP_saprc99 = .false.

GASCHEM_OPT:select case (config_gaschem_opt)
    case ('no_gaschem')
        no_gaschem = .true.
    case ('AIPC_saprc99')
        AIPC_saprc99 = .true.
    case ('KPP_saprc99')
        KPP_saprc99 = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for photochemistry schemes options(config_gaschem_opt). *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select GASCHEM_OPT

call iamas_log_write( ' no_gaschem  = $l', logicArgs=(/no_gaschem/) )
call iamas_log_write( ' AIPC_saprc99  = $l', logicArgs=(/AIPC_saprc99/) )
call iamas_log_write( ' KPP_saprc99  = $l', logicArgs=(/KPP_saprc99/) )

!--- zhxia:initialization of all packages for aerosol chemistry and interactions shemes:
! SIA 
call iamas_pool_get_config(configs,'config_aerchem_sia_opt',config_aerchem_sia_opt)
nullify(no_sia)
nullify(AIMACI_mosaic)
nullify(simple_SIA)
nullify(MOSAIC_mech)
call iamas_pool_get_package(packages, 'no_siaActive' ,  no_sia )
call iamas_pool_get_package(packages, 'AIMACI_mosaicActive' ,  AIMACI_mosaic )
call iamas_pool_get_package(packages, 'simple_SIAActive' ,  simple_SIA )
call iamas_pool_get_package(packages, 'MOSAIC_mechActive' ,  MOSAIC_mech )

if(.not.associated(no_sia) .or. &
   .not.associated(AIMACI_mosaic) .or. &
   .not.associated(simple_SIA) .or. &
   .not.associated(MOSAIC_mech )) then

   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for secondary inorganic aerosol chemistry and interactions schemes options(config_aerchem_sia_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_sia = .false.
AIMACI_mosaic = .false.
simple_SIA = .false.
MOSAIC_mech = .false.

AERCHEM_SIA_OPT:select case (config_aerchem_sia_opt)
    case ('no_sia')
        no_sia = .true.
    case ('AIMACI_mosaic')
        AIMACI_mosaic = .true.
    case ('simple_SIA')
        simple_SIA = .true.
    case ('MOSAIC_mech')
        MOSAIC_mech = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for inorganic aerosol chemistry and interactions schemes options(config_aerchem_sia_opt). *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select AERCHEM_SIA_OPT

call iamas_log_write( ' no_sia  = $l', logicArgs=(/no_sia/) )
call iamas_log_write( ' AIMACI_mosaic  = $l', logicArgs=(/AIMACI_mosaic/) )
call iamas_log_write( ' simple_SIA  = $l', logicArgs=(/simple_SIA/) )
call iamas_log_write( ' MOSAIC_mech  = $l', logicArgs=(/MOSAIC_mech/) )

! SOA 
call iamas_pool_get_config(configs,'config_aerchem_soa_opt',config_aerchem_soa_opt)
nullify(no_soa)
nullify(simple_SOA)
call iamas_pool_get_package(packages, 'no_soaActive' ,  no_soa )
call iamas_pool_get_package(packages, 'simple_SOAActive' ,  simple_SOA )

if(.not.associated(no_soa) .or. &
   .not.associated(simple_SOA)) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for organic aerosol chemistry and interactions schemes options(config_aerchem_soa_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_soa = .false.
simple_SOA = .false.

AERCHEM_SOA_OPT:select case (config_aerchem_soa_opt)
    case ('no_soa')
        no_soa = .true.
    case ('simple_SOA')
        simple_SOA = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for organic aerosol chemistry and interactions schemes options(config_aerchem_soa_opt). *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select AERCHEM_SOA_OPT

call iamas_log_write( ' no_soa  = $l', logicArgs=(/no_soa/) )
call iamas_log_write( ' simple_SOA  = $l', logicArgs=(/simple_SOA/) )


!--- zhxia:initialization of all packages for stratospheric photochemistry shemes:

call iamas_pool_get_config(configs,'config_stra_gaschem_opt',config_stra_gaschem_opt)
nullify(no_stra_gaschem)
nullify(linoz)
call iamas_pool_get_package(packages, 'no_stra_gaschemActive' ,  no_stra_gaschem )
call iamas_pool_get_package(packages, 'linozActive' ,  linoz )

if(.not.associated(no_stra_gaschem) .or. &
   .not.associated(linoz)) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for photochemistry schemes options(config_stra_gaschem_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_stra_gaschem = .false.
linoz = .false.

STRA_GASCHEM_OPT:select case (config_stra_gaschem_opt)
    case ('no_stra_gaschem')
        no_stra_gaschem = .true.
    case ('linoz')
        linoz = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for stratospheric photochemistry schemes options(config_stra_gaschem_opt). *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select STRA_GASCHEM_OPT

call iamas_log_write( ' no_stra_gaschem  = $l', logicArgs=(/no_stra_gaschem/) )
call iamas_log_write( ' linoz  = $l', logicArgs=(/linoz/) )

!---xyruan: initialization of all packages for cloud aqueous chemistry schemes:

call iamas_pool_get_config(configs,'config_cloudchem_opt',config_cloudchem_opt)
nullify(no_cloudchem)
nullify(CMU_cloudchem)
call iamas_pool_get_package(packages, 'no_cloudchemActive' ,  no_cloudchem )
call iamas_pool_get_package(packages, 'CMU_cloudchemActive' ,  CMU_cloudchem )

if(.not.associated(no_cloudchem) .or. &
   .not.associated(CMU_cloudchem)) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for cloudchem options. *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_cloudchem = .false.
CMU_cloudchem = .false.

CLOUDCHEM_OPT:select case (config_cloudchem_opt)
    case ('no_cloudchem')
        no_cloudchem = .true.
    case ('CMU_cloudchem')
        CMU_cloudchem = .true.
    case default
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('* Error while setting up packages for cloudchem options. *',messageType=iAMAS_LOG_ERR)
        call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select CLOUDCHEM_OPT

call iamas_log_write( ' no_cloudchem  = $l', logicArgs=(/no_cloudchem/) )
call iamas_log_write( ' CMU_cloudchem  = $l', logicArgs=(/CMU_cloudchem/) )

!--- initialization of all packages for chem_opt:

call iamas_pool_get_config(configs,'config_chem_opt',config_chem_opt)

nullify(no_chemistry, notIncluded)
nullify(bin10aer_nature_aq , bin09aer_nature_aq , bin08aer_nature_aq , bin05aer_nature_aq , bin04aer_nature_aq )
nullify(bin10aer_nature    , bin09aer_nature    , bin08aer_nature    , bin05aer_nature    , bin04aer_nature    )
nullify(bin10dust_aq, bin09dust_aq, bin08dust_aq, bin05dust_aq, bin04dust_aq)
nullify(bin10dust   , bin09dust   , bin08dust   , bin05dust   , bin04dust   )
nullify(bin10seas_aq, bin09seas_aq, bin08seas_aq, bin05seas_aq, bin04seas_aq)
nullify(bin10seas   , bin09seas   , bin08seas   , bin05seas   , bin04seas   )
nullify(bin04aer_gas_aq, bin04aer_gas )
call iamas_pool_get_package(packages, 'no_chemistryActive'       ,  no_chemistry       )
call iamas_pool_get_package(packages, 'notIncludedActive'        ,  notIncluded        )
call iamas_pool_get_package(packages, 'bin10aer_nature_aqActive' ,  bin10aer_nature_aq )
call iamas_pool_get_package(packages, 'bin09aer_nature_aqActive' ,  bin09aer_nature_aq )
call iamas_pool_get_package(packages, 'bin08aer_nature_aqActive' ,  bin08aer_nature_aq )
call iamas_pool_get_package(packages, 'bin05aer_nature_aqActive' ,  bin05aer_nature_aq )
call iamas_pool_get_package(packages, 'bin04aer_nature_aqActive' ,  bin04aer_nature_aq )
call iamas_pool_get_package(packages, 'bin10aer_natureActive'    ,  bin10aer_nature    )
call iamas_pool_get_package(packages, 'bin09aer_natureActive'    ,  bin09aer_nature    )
call iamas_pool_get_package(packages, 'bin08aer_natureActive'    ,  bin08aer_nature    )
call iamas_pool_get_package(packages, 'bin05aer_natureActive'    ,  bin05aer_nature    )
call iamas_pool_get_package(packages, 'bin04aer_natureActive'    ,  bin04aer_nature    )
call iamas_pool_get_package(packages, 'bin10dust_aqActive'       ,  bin10dust_aq       )
call iamas_pool_get_package(packages, 'bin09dust_aqActive'       ,  bin09dust_aq       )
call iamas_pool_get_package(packages, 'bin08dust_aqActive'       ,  bin08dust_aq       )
call iamas_pool_get_package(packages, 'bin05dust_aqActive'       ,  bin05dust_aq       )
call iamas_pool_get_package(packages, 'bin04dust_aqActive'       ,  bin04dust_aq       )
call iamas_pool_get_package(packages, 'bin10dustActive'          ,  bin10dust          )
call iamas_pool_get_package(packages, 'bin09dustActive'          ,  bin09dust          )
call iamas_pool_get_package(packages, 'bin08dustActive'          ,  bin08dust          )
call iamas_pool_get_package(packages, 'bin05dustActive'          ,  bin05dust          )
call iamas_pool_get_package(packages, 'bin04dustActive'          ,  bin04dust          )
call iamas_pool_get_package(packages, 'bin10seas_aqActive'       ,  bin10seas_aq       )
call iamas_pool_get_package(packages, 'bin09seas_aqActive'       ,  bin09seas_aq       )
call iamas_pool_get_package(packages, 'bin08seas_aqActive'       ,  bin08seas_aq       )
call iamas_pool_get_package(packages, 'bin05seas_aqActive'       ,  bin05seas_aq       )
call iamas_pool_get_package(packages, 'bin04seas_aqActive'       ,  bin04seas_aq       )
call iamas_pool_get_package(packages, 'bin10seasActive'          ,  bin10seas          )
call iamas_pool_get_package(packages, 'bin09seasActive'          ,  bin09seas          )
call iamas_pool_get_package(packages, 'bin08seasActive'          ,  bin08seas          )
call iamas_pool_get_package(packages, 'bin05seasActive'          ,  bin05seas          )
call iamas_pool_get_package(packages, 'bin04seasActive'          ,  bin04seas          )
call iamas_pool_get_package(packages, 'bin04aer_gas_aqActive'    ,  bin04aer_gas_aq    )
call iamas_pool_get_package(packages, 'bin04aer_gasActive'       ,  bin04aer_gas       )

if(.not.associated(no_chemistry       ) .or. &
   .not.associated(notIncluded        ) .or. &
   .not.associated(bin10aer_nature_aq ) .or. &
   .not.associated(bin09aer_nature_aq ) .or. &
   .not.associated(bin08aer_nature_aq ) .or. &
   .not.associated(bin05aer_nature_aq ) .or. &
   .not.associated(bin04aer_nature_aq ) .or. &
   .not.associated(bin10aer_nature    ) .or. &
   .not.associated(bin09aer_nature    ) .or. &
   .not.associated(bin08aer_nature    ) .or. &
   .not.associated(bin05aer_nature    ) .or. &
   .not.associated(bin04aer_nature    ) .or. &
   .not.associated(bin10dust_aq       ) .or. &
   .not.associated(bin09dust_aq       ) .or. &
   .not.associated(bin08dust_aq       ) .or. &
   .not.associated(bin05dust_aq       ) .or. &
   .not.associated(bin04dust_aq       ) .or. &
   .not.associated(bin10dust          ) .or. &
   .not.associated(bin09dust          ) .or. &
   .not.associated(bin08dust          ) .or. &
   .not.associated(bin05dust          ) .or. &
   .not.associated(bin04dust          ) .or. &
   .not.associated(bin10seas_aq       ) .or. &
   .not.associated(bin09seas_aq       ) .or. &
   .not.associated(bin08seas_aq       ) .or. &
   .not.associated(bin05seas_aq       ) .or. &
   .not.associated(bin04seas_aq       ) .or. &
   .not.associated(bin10seas          ) .or. &
   .not.associated(bin09seas          ) .or. &
   .not.associated(bin08seas          ) .or. &
   .not.associated(bin05seas          ) .or. &
   .not.associated(bin04seas          ) .or. &
   .not.associated(bin04aer_gas_aq    ) .or. &
   .not.associated(bin04aer_gas       )) then
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('* Error while setting up packages for chem options1(config_chem_opt). *',messageType=iAMAS_LOG_ERR)
   call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
   ierr = 1
   return
endif

no_chemistry        = .false.
notIncluded         = .false.
bin10aer_nature_aq  = .false.
bin09aer_nature_aq  = .false.
bin08aer_nature_aq  = .false.
bin05aer_nature_aq  = .false.
bin04aer_nature_aq  = .false.
bin10aer_nature     = .false.
bin09aer_nature     = .false.
bin08aer_nature     = .false.
bin05aer_nature     = .false.
bin04aer_nature     = .false.
bin10dust_aq        = .false.
bin09dust_aq        = .false.
bin08dust_aq        = .false.
bin05dust_aq        = .false.
bin04dust_aq        = .false.
bin10dust           = .false.
bin09dust           = .false.
bin08dust           = .false.
bin05dust           = .false.
bin04dust           = .false.
bin10seas_aq        = .false.
bin09seas_aq        = .false.
bin08seas_aq        = .false.
bin05seas_aq        = .false.
bin04seas_aq        = .false.
bin10seas           = .false.
bin09seas           = .false.
bin08seas           = .false.
bin05seas           = .false.
bin04seas           = .false.
bin04aer_gas_aq     = .false.
bin04aer_gas        = .false.

CHEM_OPT:select case (config_chem_opt)
    case ('no_chemistry')
        no_chemistry = .true.
    case ('notIncluded')
        notIncluded = .true.
    case ('bin10aer_nature_aq')
        bin10aer_nature_aq  = .true.
    case ('bin09aer_nature_aq')
        bin09aer_nature_aq  = .true.
    case ('bin08aer_nature_aq')
        bin08aer_nature_aq  = .true.
    case ('bin05aer_nature_aq')
        bin05aer_nature_aq  = .true.
    case ('bin04aer_nature_aq')
        bin04aer_nature_aq  = .true.
    case ('bin10aer_nature')
        bin10aer_nature     = .true.
    case ('bin09aer_nature')
        bin09aer_nature     = .true.
    case ('bin08aer_nature')
        bin08aer_nature     = .true.
    case ('bin05aer_nature')
        bin05aer_nature     = .true.
    case ('bin04aer_nature')
        bin04aer_nature     = .true.
    case ('bin10dust_aq')
        bin10dust_aq  = .true.
    case ('bin09dust_aq')
        bin09dust_aq  = .true.
    case ('bin08dust_aq')
        bin08dust_aq  = .true.
    case ('bin05dust_aq')
        bin05dust_aq  = .true.
    case ('bin04dust_aq')
        bin04dust_aq  = .true.
    case ('bin10dust')
        bin10dust     = .true.
    case ('bin09dust')
        bin09dust     = .true.
    case ('bin08dust')
        bin08dust     = .true.
    case ('bin05dust')
        bin05dust     = .true.
    case ('bin04dust')
        bin04dust     = .true.
    case ('bin10seas_aq')
        bin10seas_aq  = .true.
    case ('bin09seas_aq')
        bin09seas_aq  = .true.
    case ('bin08seas_aq')
        bin08seas_aq  = .true.
    case ('bin05seas_aq')
        bin05seas_aq  = .true.
    case ('bin04seas_aq')
        bin04seas_aq  = .true.
    case ('bin10seas')
        bin10seas     = .true.
    case ('bin09seas')
        bin09seas     = .true.
    case ('bin08seas')
        bin08seas     = .true.
    case ('bin05seas')
        bin05seas     = .true.
    case ('bin04seas')
        bin04seas     = .true.
    case ('bin04aer_gas_aq')
        bin04aer_gas_aq = .true.
    case ('bin04aer_gas')
        bin04aer_gas = .true.
    case default
     call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
     call iamas_log_write('* Error while setting up packages for chem options2(config_chem_opt). *',messageType=iAMAS_LOG_ERR)
     call iamas_log_write('=================================================',messageType=iAMAS_LOG_ERR)
end select CHEM_OPT

call iamas_log_write( '    no_chemistry    = $l', logicArgs=(/no_chemistry/) )
call iamas_log_write( '    notIncluded     = $l', logicArgs=(/notIncluded/) )
call iamas_log_write( '    bin10aer_nature_aq     = $l', logicArgs=(/bin10aer_nature_aq /) )
call iamas_log_write( '    bin09aer_nature_aq     = $l', logicArgs=(/bin09aer_nature_aq /) )
call iamas_log_write( '    bin08aer_nature_aq     = $l', logicArgs=(/bin08aer_nature_aq /) )
call iamas_log_write( '    bin05aer_nature_aq     = $l', logicArgs=(/bin05aer_nature_aq /) )
call iamas_log_write( '    bin04aer_nature_aq     = $l', logicArgs=(/bin04aer_nature_aq /) )   
call iamas_log_write( '    bin10aer_nature        = $l', logicArgs=(/bin10aer_nature    /) )  
call iamas_log_write( '    bin09aer_nature        = $l', logicArgs=(/bin09aer_nature    /) )  
call iamas_log_write( '    bin08aer_nature        = $l', logicArgs=(/bin08aer_nature    /) )  
call iamas_log_write( '    bin05aer_nature        = $l', logicArgs=(/bin05aer_nature    /) )  
call iamas_log_write( '    bin04aer_nature        = $l', logicArgs=(/bin04aer_nature    /) )  
call iamas_log_write( '    bin10dust_aq    = $l', logicArgs=(/bin10dust_aq/) )  
call iamas_log_write( '    bin09dust_aq    = $l', logicArgs=(/bin09dust_aq/) )  
call iamas_log_write( '    bin08dust_aq    = $l', logicArgs=(/bin08dust_aq/) )  
call iamas_log_write( '    bin05dust_aq    = $l', logicArgs=(/bin05dust_aq/) )  
call iamas_log_write( '    bin04dust_aq    = $l', logicArgs=(/bin04dust_aq/) )  
call iamas_log_write( '    bin10dust       = $l', logicArgs=(/bin10dust   /) )  
call iamas_log_write( '    bin09dust       = $l', logicArgs=(/bin09dust   /) )  
call iamas_log_write( '    bin08dust       = $l', logicArgs=(/bin08dust   /) )  
call iamas_log_write( '    bin05dust       = $l', logicArgs=(/bin05dust   /) )  
call iamas_log_write( '    bin04dust       = $l', logicArgs=(/bin04dust   /) )  
call iamas_log_write( '    bin10seas_aq    = $l', logicArgs=(/bin10seas_aq/) )  
call iamas_log_write( '    bin09seas_aq    = $l', logicArgs=(/bin09seas_aq/) )  
call iamas_log_write( '    bin08seas_aq    = $l', logicArgs=(/bin08seas_aq/) )  
call iamas_log_write( '    bin05seas_aq    = $l', logicArgs=(/bin05seas_aq/) )  
call iamas_log_write( '    bin04seas_aq    = $l', logicArgs=(/bin04seas_aq/) )  
call iamas_log_write( '    bin10seas       = $l', logicArgs=(/bin10seas   /) )  
call iamas_log_write( '    bin09seas       = $l', logicArgs=(/bin09seas   /) )  
call iamas_log_write( '    bin08seas       = $l', logicArgs=(/bin08seas   /) )  
call iamas_log_write( '    bin05seas       = $l', logicArgs=(/bin05seas   /) )  
call iamas_log_write( '    bin04seas       = $l', logicArgs=(/bin04seas   /) )  
call iamas_log_write( '    bin04aer_gas_aq    = $l', logicArgs=(/bin04aer_gas_aq  /) )  
call iamas_log_write( '    bin04aer_gas    = $l', logicArgs=(/bin04aer_gas  /) ) 

end function atmchem_setup_packages

!=================================================================================================================
 end module iamas_atmchem_packages
!=================================================================================================================
