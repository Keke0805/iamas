! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_driver_lsm

!Added by qkyang to support urban model+++++++++++++++++++++
 use iamas_atmphys_manager, only: gmt,curr_julday,julday,year
!Added by qkyang to support urban model---------------------



 use iamas_kind_types
 use iamas_pool_routines

 use iamas_atmphys_constants
 use iamas_atmphys_landuse, only: isurban
 use iamas_atmphys_lsm_noahinit
 use iamas_atmphys_vars

!wrf physics
 use module_sf_noahdrv
 use module_sf_noah_seaice_drv
 use module_sf_sfcdiags
 
 implicit none
 private
 public:: init_lsm,       &
 !qkyang comment to speed iAMAS model++++++++++++
!           allocate_lsm,   &
!           deallocate_lsm, &
!qkyang comment to speed iAMAS model-------------
          driver_lsm

!qkyang comment to support urban model++++++++++++
!urban physics: iAMAS does not plan to run the urban physics option. 
!integer,parameter:: sf_urban_physics = 0 !activate urban canopy model (=0: no urban canopy) 
!qkyang replaced sf_urban_physics with config_urban_scheme
!qkyang comment to support urban model-------------
 
 logical,parameter:: rdmaxalb = .false. !use snow albedo from geogrid;false use table values
 logical,parameter:: myj      = .false. !true if using Mellor-Yamada PBL scheme.
 logical,parameter:: frpcpn   = .false.
 logical,parameter:: rdlai2d  = .false.


!iAMAS driver for parameterization of land surface processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_driver_lsm:
! ---------------------------------------
! allocate_lsm  : allocate local arrays for land surface parameterization.
! deallocate_lsm: deallocate local arrays for land surface parameterization.
! driver_lsm    : main driver (called from subroutine physics_driver).
! lsm_from_iAMAS : initialize local arrays.
! lsm_to_iAMAS   : copy local arrays to iAMAS arrays.
!
! WRF physics called from driver_lsm:
! ------------------------ ----------
!    * module_sf_noahdrv : NOAH 4-layers land surface scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine lsm.
!   Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the definition of isurban to landuse_init_foriAMAS in iamas_atmphys_landuse.F.
!   isurban is now defined as a function of the input landuse data file.
!   Dominikus Heinzeller (IMK) / 2014-07-24.
! * removed the global variable qz0 and initialized the local variable qz0_p to 0. qz0 is only
!   used in the MYJ PBL parameterization which is not available in iAMAS.
!   Laura D. Fowler (laura@ucar.edu) / 2015-03-05.
! * in subroutine lsm_from_iAMAS, modified the initialization of local variable dzs_p.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-11.
! * in subroutine lsm_to_iAMAS, removed updating isltyp, ivgtyp, shdmin, and shdmax as they
!   are constant in the Noah lsm.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-11.
! * in subroutine lsm_from_iAMAS, modified the calculation of rain_bl, now that the convection
!   and cloud microphysics parameterizations are contained in "packages." The original
!   calculation failed when configuration_convection_scheme or config_microphysics_scheme
!   was set of 'off'.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-13.
! * added call to the subroutine sfcdiags to update the 2-meter temperature, potential temperature, and
!   water vapor mixing ratio after call to lsm.
!   Laura D. Fowler (laura@ucar.edu) / 2016-05-11.
! * added the calculation of surface variables over seaice cells when config_frac_seaice is set to true.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-03.
! * now use isice and iswater initialized in the init file instead of initialized in iamas_atmphys_landuse.F.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-13.
! * updated the call to subroutine lsm as we updated the Noah LSM scheme to WRF version 3.9.0.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-30.
! * since we removed the local variable lsm_scheme from iamas_atmphys_vars.F, now defines lsm_scheme as a
!   pointer to config_lsm_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * added call to seaice_noah to include the parameterization of seaice for the updated Noah land surface
!   scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-19.

!
! DOCUMENTATION:
! ./physics_wrf/module_sf_noahdrv.F: main driver for the "NOAH" land-surface parameterization.
! In the argument list,I added "OPTIONAL" to the declaration of a few arrays to avoid compiling
! with the "urban physics" option. These arrays are:
! .. num_roof_layers; num_wall_layers; num_road_layers;num_urban_layers.
! .. ust_urb2d;frc_urb2d;utype_urb2d.
! Laura D. Fowler (01-18-2011).


 contains
 
!=================================================================================================================
 subroutine init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(iamas_pool_type):: mesh
 type(iamas_pool_type),intent(in):: configs

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: sfc_input
!local pointers:
 character(len=StrKIND),pointer:: lsm_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_lsm_scheme',lsm_scheme)

 lsm_select: select case (trim(lsm_scheme))

    case ("noah")
       call noah_init_foriAMAS(dminfo,mesh,configs,diag_physics,sfc_input)
    
    case default
 
 end select lsm_select

 end subroutine init_lsm

!=================================================================================================================
 subroutine driver_lsm(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!added by qkyang to support urban model++++++++++++
integer,pointer:: config_urban_scheme
!qkyang changed to support BEP BEM model>>>>>>>>>>>>>>>>>      
integer,pointer:: num_urban_layers,num_urban_hi
real(kind=RKIND),dimension(:),pointer:: TSK_RURAL_BEP
!qkyang changed to support BEP BEM model<<<<<<<<<<<<<<<<<<
real(kind=RKIND),dimension(:),pointer:: DZR,DZB,DZG
real(kind=RKIND),pointer:: declin
integer,intent(in):: its,ite!,JTS,JTE,IMS,IME,JMS,JME
!added by qkyang to support urban model-------------

 !Changed by qkyang to speed iAMAS model+++++++++++++++++++++
!input arguments:
!Added by qkyang to support urban model+++++++++++++++++++++
 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D
 real(kind=RKIND),dimension(:),pointer:: TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,QC_URB2D,UC_URB2D
 real(kind=RKIND),dimension(:,:),pointer:: TRL_URB3D,TBL_URB3D,TGL_URB3D
 real(kind=RKIND),dimension(:),pointer:: SH_URB2D,LH_URB2D,G_URB2D,RN_URB2D,TS_URB2D
 real(kind=RKIND),dimension(:,:),pointer:: TRB_URB4D,TW1_URB4D,TW2_URB4D,TGB_URB4D
 real(kind=RKIND),dimension(:,:),pointer:: TLEV_URB3D,QLEV_URB3D,TW1LEV_URB3D,TW2LEV_URB3D,TGLEV_URB3D,TFLEV_URB3D
 real(kind=RKIND),dimension(:),pointer:: SF_AC_URB3D,LF_AC_URB3D,CM_AC_URB3D,SFVENT_URB3D,LFVENT_URB3D
 real(kind=RKIND),dimension(:,:),pointer:: SFWIN1_URB3D,SFWIN2_URB3D,SFW1_URB3D,SFW2_URB3D,SFR_URB3D,SFG_URB3D
 real(kind=RKIND),dimension(:),pointer:: LP_URB2D,LB_URB2D,HGT_URB2D,MH_URB2D,STDH_URB2D
 real(kind=RKIND),dimension(:,:),pointer:: HI_URB2D,LF_URB2D
 real(kind=RKIND),dimension(:),pointer:: CMCR_URB2D,TGR_URB2D
 real(kind=RKIND),dimension(:,:),pointer:: TGRL_URB3D,SMR_URB3D
 real(kind=RKIND),dimension(:),pointer:: DRELR_URB2D,DRELB_URB2D,DRELG_URB2D,FLXHUMR_URB2D,FLXHUMB_URB2D,FLXHUMG_URB2D
 real(kind=RKIND),dimension(:,:),pointer:: a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep
 real(kind=RKIND),dimension(:,:),pointer:: b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep
 real(kind=RKIND),dimension(:,:),pointer:: dlg_bep,dl_u_bep,sf_bep,vl_bep
 real(kind=RKIND),dimension(:),pointer:: FRC_URB2D,UST_URB2D,ust
 real(kind=RKIND),dimension(:),pointer:: u10,v10,psih,psim,gz1oz0,AKHS,AKMS
 real(kind=RKIND),dimension(:),pointer:: CMR_SFCDIF,CHR_SFCDIF,CMGR_SFCDIF,CHGR_SFCDIF,CMC_SFCDIF,CHC_SFCDIF
 integer,dimension(:),pointer:: UTYPE_URB2D
 real(kind=RKIND),dimension(:),pointer:: coszr,hrang
 integer:: k
!Added by qkyang to support urban model---------------------

!Changed by qkyang to speed iAMAS model+++++++++++++++++++++
 REAL,    DIMENSION(its:ite)::                                & 
                                                      swdown, &
                                                      rainbl, &
                                                         qz0, &
                                                     flxsnow, &
                                                     fvbsnow, &
                                                    fbursnow, &
                                                    fgsnsnow, &
                                                     infxsrt, &
                                                   sfcheadrt, &
                                                    soldrain, &
                                                    icedepth, &
                                                      snowsi, &
                                                  PSIM_URB2D, &
                                                  PSIH_URB2D, &
                                                GZ1OZ0_URB2D, &
                                                   U10_URB2D, &
                                                   V10_URB2D, &
                                                   TH2_URB2D, &
                                                    Q2_URB2D, &
                                                  AKMS_URB2D, &
                                                     tsk_sea, &
                                                     tsk_ice, &
                                                        xlat, &
                                                       xlong
!Changed by qkyang to speed iAMAS model--------------------- 

!qkyang added to support anthropogenic heating emissions++++++
 real(kind=RKIND),dimension(:),pointer  :: AH_emis
 logical,pointer:: is_ah_emis
 real(kind=RKIND),dimension(:),pointer  :: AH_URB2D
!qkyang added to support anthropogenic heating emissions------
!local pointers:

 character(len=StrKIND),pointer:: config_microp_scheme,    &
                                  config_convection_scheme

 integer,dimension(:),pointer:: isltyp,ivgtyp

 real(kind=RKIND),dimension(:),pointer  :: acsnom,acsnow,canwat,chs,chs2,chklowq,cpm,cqs2,glw,   &
                                           grdflx,gsw,hfx,lai,lh,noahres,potevp,qfx,qgh,qsfc,    &
                                           br,sfc_albedo,sfc_albedo_seaice,sfc_emibck,sfc_emiss, &
                                           sfcrunoff,smstav,smstot,snotime,snopcx,sr,udrunoff,   &
                                           z0,znt
 real(kind=RKIND),dimension(:),pointer  :: shdmin,shdmax,snoalb,sfc_albbck,snow,snowc,snowh,tmn, &
                                           skintemp,vegfra,xice,xland
 real(kind=RKIND),dimension(:),pointer  :: t2m,th2m,q2
 real(kind=RKIND),dimension(:),pointer  :: raincv,rainncv
 real(kind=RKIND),dimension(:,:),pointer:: sh2o,smcrel,smois,tslb,dzs

!local variables and arrays:
 logical:: do_fill
 integer:: ip,iEdg
 integer:: i,j,n

!-----------------------------------------------------------------------------------------------------------------
 !Changed by qkyang to speed iAMAS model---------------------

!input arguments:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs

 integer,intent(in):: itimestep

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_sfc_albedo,config_frac_seaice
 character(len=StrKIND),pointer:: lsm_scheme
 character(len=StrKIND),pointer:: mminlu
 integer,pointer:: isice,iswater !qkyang added iswater to support LCZ landuse dataset

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_sfc_albedo' ,config_sfc_albedo )
 call iamas_pool_get_config(configs,'config_frac_seaice',config_frac_seaice)
 call iamas_pool_get_config(configs,'config_lsm_scheme',lsm_scheme)

!added by qkyang to support urban model++++++++++++
 call iamas_pool_get_config(configs,'config_urban_scheme',config_urban_scheme)
 !qkyang added to support BEP BEM model>>>>>>>>>>>>>>>>>
 call iamas_pool_get_config(configs,'num_urban_layers',num_urban_layers)
 call iamas_pool_get_config(configs,'num_urban_hi',num_urban_hi)
 call iamas_pool_get_array(diag_physics,'TSK_RURAL_BEP' , TSK_RURAL_BEP )
 !qkyang added to support BEP BEM model<<<<<<<<<<<<<<<<<<
 call iamas_pool_get_array(diag_physics,'DZR' , DZR )
 call iamas_pool_get_array(diag_physics,'DZB' , DZB )
 call iamas_pool_get_array(diag_physics,'DZG' , DZG )
 call iamas_pool_get_array(diag_physics,'declin' , declin )
!added by qkyang to support urban model-------------
 call iamas_pool_get_array(sfc_input,'mminlu',mminlu)
 call iamas_pool_get_array(sfc_input,'isice' ,isice )
 call iamas_pool_get_array(sfc_input,'iswater',iswater )!qkyang added to support LCZ landuse dataset

 !Changed by qkyang to speed iAMAS model+++++++++++++++++++++
 !added by qkyang to support urban model+++++++++++++++++++++
!  call iamas_pool_get_config(configs,'config_urban_scheme',config_urban_scheme)
! IF ((config_urban_scheme.eq.1).OR.(config_urban_scheme.EQ.2).OR.(config_urban_scheme.EQ.3)) THEN
 call iamas_pool_get_array(mesh,'latCell',latCell)
 call iamas_pool_get_array(mesh,'lonCell',lonCell)
 call iamas_pool_get_array(diag_physics,'XXXR_URB2D' , XXXR_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXB_URB2D' , XXXB_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXG_URB2D' , XXXG_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXC_URB2D' , XXXC_URB2D )
 call iamas_pool_get_array(diag_physics,'TR_URB2D' , TR_URB2D )
 call iamas_pool_get_array(diag_physics,'TB_URB2D' , TB_URB2D )
 call iamas_pool_get_array(diag_physics,'TG_URB2D' , TG_URB2D )
 call iamas_pool_get_array(diag_physics,'TC_URB2D' , TC_URB2D )
 call iamas_pool_get_array(diag_physics,'QC_URB2D' , QC_URB2D )
 call iamas_pool_get_array(diag_physics,'UC_URB2D' , UC_URB2D )
 call iamas_pool_get_array(diag_physics,'TRL_URB3D' , TRL_URB3D )
 call iamas_pool_get_array(diag_physics,'TBL_URB3D' , TBL_URB3D )
 call iamas_pool_get_array(diag_physics,'TGL_URB3D' , TGL_URB3D )
 call iamas_pool_get_array(diag_physics,'SH_URB2D' , SH_URB2D )
 call iamas_pool_get_array(diag_physics,'LH_URB2D' , LH_URB2D )
 call iamas_pool_get_array(diag_physics,'G_URB2D' , G_URB2D )
 call iamas_pool_get_array(diag_physics,'RN_URB2D' , RN_URB2D )
 call iamas_pool_get_array(diag_physics,'TS_URB2D' , TS_URB2D )
 call iamas_pool_get_array(diag_physics,'TRB_URB4D' , TRB_URB4D )
 call iamas_pool_get_array(diag_physics,'TW1_URB4D' , TW1_URB4D )
 call iamas_pool_get_array(diag_physics,'TW2_URB4D' , TW2_URB4D )
 call iamas_pool_get_array(diag_physics,'TGB_URB4D' , TGB_URB4D )
 call iamas_pool_get_array(diag_physics,'TLEV_URB3D' , TLEV_URB3D )
 call iamas_pool_get_array(diag_physics,'QLEV_URB3D' , QLEV_URB3D )
 call iamas_pool_get_array(diag_physics,'TW1LEV_URB3D' , TW1LEV_URB3D )
 call iamas_pool_get_array(diag_physics,'TW2LEV_URB3D' , TW2LEV_URB3D )
 call iamas_pool_get_array(diag_physics,'TGLEV_URB3D' , TGLEV_URB3D )
 call iamas_pool_get_array(diag_physics,'TFLEV_URB3D' , TFLEV_URB3D )
 call iamas_pool_get_array(diag_physics,'SF_AC_URB3D' , SF_AC_URB3D )
 call iamas_pool_get_array(diag_physics,'LF_AC_URB3D' , LF_AC_URB3D )
 call iamas_pool_get_array(diag_physics,'CM_AC_URB3D' , CM_AC_URB3D )
 call iamas_pool_get_array(diag_physics,'SFVENT_URB3D' , SFVENT_URB3D )
 call iamas_pool_get_array(diag_physics,'LFVENT_URB3D' , LFVENT_URB3D )
 call iamas_pool_get_array(diag_physics,'SFWIN1_URB3D' , SFWIN1_URB3D )
 call iamas_pool_get_array(diag_physics,'SFWIN2_URB3D' , SFWIN2_URB3D )
 call iamas_pool_get_array(diag_physics,'SFW1_URB3D' , SFW1_URB3D )
 call iamas_pool_get_array(diag_physics,'SFW2_URB3D' , SFW2_URB3D )
 call iamas_pool_get_array(diag_physics,'SFR_URB3D' , SFR_URB3D )
 call iamas_pool_get_array(diag_physics,'SFG_URB3D' , SFG_URB3D )
 call iamas_pool_get_array(diag_physics,'LP_URB2D' , LP_URB2D )
 call iamas_pool_get_array(diag_physics,'HI_URB2D' , HI_URB2D )
 call iamas_pool_get_array(diag_physics,'LB_URB2D' , LB_URB2D )
 call iamas_pool_get_array(diag_physics,'HGT_URB2D' , HGT_URB2D )
 call iamas_pool_get_array(diag_physics,'MH_URB2D' , MH_URB2D )
 call iamas_pool_get_array(diag_physics,'STDH_URB2D' , STDH_URB2D )
 call iamas_pool_get_array(diag_physics,'LF_URB2D' , LF_URB2D )
 call iamas_pool_get_array(diag_physics,'CMCR_URB2D' , CMCR_URB2D )
 call iamas_pool_get_array(diag_physics,'TGR_URB2D' , TGR_URB2D )
 call iamas_pool_get_array(diag_physics,'TGRL_URB3D' , TGRL_URB3D )
 call iamas_pool_get_array(diag_physics,'SMR_URB3D' , SMR_URB3D )
 call iamas_pool_get_array(diag_physics,'DRELR_URB2D' , DRELR_URB2D )
 call iamas_pool_get_array(diag_physics,'DRELB_URB2D' , DRELB_URB2D )
 call iamas_pool_get_array(diag_physics,'DRELG_URB2D' , DRELG_URB2D )
 call iamas_pool_get_array(diag_physics,'FLXHUMR_URB2D' , FLXHUMR_URB2D )
 call iamas_pool_get_array(diag_physics,'FLXHUMB_URB2D' , FLXHUMB_URB2D )
 call iamas_pool_get_array(diag_physics,'FLXHUMG_URB2D' , FLXHUMG_URB2D )
 call iamas_pool_get_array(diag_physics,'A_U_BEP' , A_U_BEP )
 call iamas_pool_get_array(diag_physics,'A_V_BEP' , A_V_BEP )
 call iamas_pool_get_array(diag_physics,'A_T_BEP' , A_T_BEP )
 call iamas_pool_get_array(diag_physics,'A_Q_BEP' , A_Q_BEP )
 call iamas_pool_get_array(diag_physics,'A_E_BEP' , A_E_BEP )
 call iamas_pool_get_array(diag_physics,'B_U_BEP' , B_U_BEP )
 call iamas_pool_get_array(diag_physics,'B_V_BEP' , B_V_BEP )
 call iamas_pool_get_array(diag_physics,'B_T_BEP' , B_T_BEP )
 call iamas_pool_get_array(diag_physics,'B_Q_BEP' , B_Q_BEP )
 call iamas_pool_get_array(diag_physics,'B_E_BEP' , B_E_BEP )
 call iamas_pool_get_array(diag_physics,'DLG_BEP' , DLG_BEP )
 call iamas_pool_get_array(diag_physics,'DL_U_BEP' , DL_U_BEP )
 call iamas_pool_get_array(diag_physics,'SF_BEP' , SF_BEP )
 call iamas_pool_get_array(diag_physics,'VL_BEP' , VL_BEP )
 call iamas_pool_get_array(diag_physics,'FRC_URB2D' , FRC_URB2D )
 call iamas_pool_get_array(diag_physics,'UST_URB2D' , UST_URB2D )
 call iamas_pool_get_array(diag_physics,'ust' , ust )
 call iamas_pool_get_array(diag_physics,'u10'   ,u10   )
 call iamas_pool_get_array(diag_physics,'v10'   ,v10   )
 call iamas_pool_get_array(diag_physics,'psih'  ,psih    )
 call iamas_pool_get_array(diag_physics,'psim'  ,psim    )
 call iamas_pool_get_array(diag_physics,'gz1oz0'  ,gz1oz0  )
 call iamas_pool_get_array(diag_physics,'AKHS'    ,AKHS    )
 call iamas_pool_get_array(diag_physics,'AKMS'    ,AKMS    )
 call iamas_pool_get_array(diag_physics,'UTYPE_URB2D' , UTYPE_URB2D )
 call iamas_pool_get_array(diag_physics,'CMR_SFCDIF',CMR_SFCDIF)
 call iamas_pool_get_array(diag_physics,'CHR_SFCDIF',CHR_SFCDIF)
 call iamas_pool_get_array(diag_physics,'CMGR_SFCDIF',CMGR_SFCDIF)
 call iamas_pool_get_array(diag_physics,'CHGR_SFCDIF',CHGR_SFCDIF)
 call iamas_pool_get_array(diag_physics,'CMC_SFCDIF',CMC_SFCDIF)
 call iamas_pool_get_array(diag_physics,'CHC_SFCDIF',CHC_SFCDIF)
 call iamas_pool_get_array(diag_physics,'hrang'     ,hrang     )
 call iamas_pool_get_array(diag_physics,'coszr'     ,coszr     )
! ENDIF
 !added by qkyang to support urban model---------------------
 
!qkyang added to support anthropogenic heating emissions++++++
 call iamas_pool_get_array(sfc_input,'AH_emis'   ,AH_emis   )
 call iamas_pool_get_config(configs,'is_ah_emis' ,is_ah_emis)
 call iamas_pool_get_array(diag_physics,'AH_URB2D' ,AH_URB2D)
!qkyang added to support anthropogenic heating emissions------

 call iamas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call iamas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )

 call iamas_pool_get_array(diag_physics,'acsnom'           ,acsnom           )
 call iamas_pool_get_array(diag_physics,'acsnow'           ,acsnow           )
 call iamas_pool_get_array(diag_physics,'canwat'           ,canwat           )
 call iamas_pool_get_array(diag_physics,'chs'              ,chs              )
 call iamas_pool_get_array(diag_physics,'chs2'             ,chs2             )
 call iamas_pool_get_array(diag_physics,'chklowq'          ,chklowq          )
 call iamas_pool_get_array(diag_physics,'cpm'              ,cpm              )
 call iamas_pool_get_array(diag_physics,'cqs2'             ,cqs2             )
 call iamas_pool_get_array(diag_physics,'glw'              ,glw              )
 call iamas_pool_get_array(diag_physics,'grdflx'           ,grdflx           )
 call iamas_pool_get_array(diag_physics,'gsw'              ,gsw              )
 call iamas_pool_get_array(diag_physics,'hfx'              ,hfx              )
 call iamas_pool_get_array(diag_physics,'lai'              ,lai              )
 call iamas_pool_get_array(diag_physics,'lh'               ,lh               )
 call iamas_pool_get_array(diag_physics,'noahres'          ,noahres          )
 call iamas_pool_get_array(diag_physics,'potevp'           ,potevp           )
 call iamas_pool_get_array(diag_physics,'qfx'              ,qfx              )
 call iamas_pool_get_array(diag_physics,'qgh'              ,qgh              )
 call iamas_pool_get_array(diag_physics,'qsfc'             ,qsfc             )
 call iamas_pool_get_array(diag_physics,'br'               ,br               )
 call iamas_pool_get_array(diag_physics,'sfc_albedo'       ,sfc_albedo       )
 call iamas_pool_get_array(diag_physics,'sfc_albedo_seaice',sfc_albedo_seaice)
 call iamas_pool_get_array(diag_physics,'sfc_emibck'       ,sfc_emibck       )
 call iamas_pool_get_array(diag_physics,'sfc_emiss'        ,sfc_emiss        )
 call iamas_pool_get_array(diag_physics,'sfcrunoff'        ,sfcrunoff        )
 call iamas_pool_get_array(diag_physics,'smstav'           ,smstav           )
 call iamas_pool_get_array(diag_physics,'smstot'           ,smstot           )
 call iamas_pool_get_array(diag_physics,'snotime'          ,snotime          )
 call iamas_pool_get_array(diag_physics,'snopcx'           ,snopcx           )
 call iamas_pool_get_array(diag_physics,'udrunoff'         ,udrunoff         )
 call iamas_pool_get_array(diag_physics,'z0'               ,z0               )
 call iamas_pool_get_array(diag_physics,'znt'              ,znt              )
 call iamas_pool_get_array(diag_physics,'t2m'              ,t2m              )
 call iamas_pool_get_array(diag_physics,'th2m'             ,th2m             )
 call iamas_pool_get_array(diag_physics,'q2'               ,q2               )

 call iamas_pool_get_array(sfc_input,'isltyp'    ,isltyp    )
 call iamas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp    )
 call iamas_pool_get_array(sfc_input,'shdmin'    ,shdmin    )
 call iamas_pool_get_array(sfc_input,'shdmax'    ,shdmax    )
 call iamas_pool_get_array(sfc_input,'snoalb'    ,snoalb    )
 call iamas_pool_get_array(sfc_input,'sfc_albbck',sfc_albbck)
 call iamas_pool_get_array(sfc_input,'snow'      ,snow      )
 call iamas_pool_get_array(sfc_input,'snowc'     ,snowc     )
 call iamas_pool_get_array(sfc_input,'snowh'     ,snowh     )
 call iamas_pool_get_array(sfc_input,'tmn'       ,tmn       )
 call iamas_pool_get_array(sfc_input,'skintemp'  ,skintemp  )
 call iamas_pool_get_array(sfc_input,'vegfra'    ,vegfra    )
 call iamas_pool_get_array(sfc_input,'xice'      ,xice      )
 call iamas_pool_get_array(sfc_input,'xland'     ,xland     )
 call iamas_pool_get_array(sfc_input,'dzs'       ,dzs       )
 call iamas_pool_get_array(sfc_input,'sh2o'      ,sh2o      )
 call iamas_pool_get_array(sfc_input,'smcrel'    ,smcrel    )
 call iamas_pool_get_array(sfc_input,'smois'     ,smois     )
 call iamas_pool_get_array(sfc_input,'tslb'      ,tslb      )

!In Registry.xml, dzs is a function of nCells. In the Noah lsm scheme, dzs is independent
!of cell locations:

 do i = its,ite
    xlat(i)   = latCell(i) / degrad
    xlong= lonCell(i) / degrad
    swdown(i) = gsw(i) / (1._RKIND - sfc_albedo(i))
    rainbl(i) = 0._RKIND
    qz0(i)    = 0._RKIND
    !initialization of arrays to run the UA Noah LSM snow cover parameterization:
    flxsnow(i)    = 0._RKIND
    fvbsnow(i)    = 0._RKIND
    fbursnow(i)   = 0._RKIND
    fgsnsnow(i)   = 0._RKIND
!Changed by qkyang to support urban model+++++++++++++++
!    ust_urb2d_p(i,j)    = ust_urb2d(i)!added by qkyang
     IF (config_urban_scheme.eq.0) THEN !Added by qkyang
      ust_urb2d(i)    = 0._RKIND!qkyang replace ust_urb_p with ust_urb2d_p
      FRC_URB2D(i)    = 0._RKIND!qkyang replace frc_urb_p with frc_urb2d_p
      UTYPE_URB2D(i)  = low_density_residential!qkyang replace UTYPE_URB_p with UTYPE_URB2d_p
     ENDIF !Added by qkyang
!Changed by qkyang to support urban model---------------
    !initialization of arrays to run the Noah LSM hydrological parameterization (not currently
    !implemented in iAMAS):
    infxsrt(i)    = 0._RKIND
    sfcheadrt(i)  = 0._RKIND
    soldrain(i)   = 0._RKIND
 enddo

! do j = jts,jte
! do i = its,ite
!     tsk_p(i,j)        = skintemp(i)
!     qgh_p(i,j)        = qgh(i)
! enddo
! enddo


if(config_frac_seaice) then
    do i = its,ite
       !modify the surface albedo and surface emissivity, and surface temperatures over sea-ice points:
       if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
          sfc_albedo(i) = (sfc_albedo(i) - 0.08_RKIND*(1._RKIND-xice(i))) / xice(i)
          sfc_emiss(i)  = (sfc_emiss(i) - 0.98_RKIND*(1._RKIND-xice(i))) / xice(i)
       endif
    enddo
    !calculate sea-surface and sea-ice temperatures over sea-ice grid cells:
    call correct_tsk_over_seaice(ims,ime,jms,jme,its,ite,jts,jte, &
      xice_threshold,xice(its:ite),skintemp(its:ite),tsk_sea,tsk_ice)
    do i = its,ite
!       do j = jts,jte
!        skintemp(i) = tsk_ice(i,j)!qkyang consider jte=1, thus tsk_ice(i,1)
!     tsk_p(i,j) = tsk_ice(i,j)!qkyang consider jte=1, thus tsk_ice(i,1)
    skintemp(i) = tsk_ice(i)!qkyang consider jte=1, thus tsk_ice(i,1)
    !initialize the surface albedo, the surface albedo over snow-covered seaice, and the
    !seaice thickness.
       icedepth(i) = seaice_thickness_default
       snowsi(i)   = seaice_snowdepth_min
!       enddo
    enddo
endif


 call iamas_pool_get_array(diag_physics,'sr',sr)
 if(config_microp_scheme .ne. 'off') then
    call iamas_pool_get_array(diag_physics,'rainncv',rainncv)
    do i = its,ite
       rainbl(i) = rainbl(i) + rainncv(i)
    enddo
 else
    do i = its,ite
       sr(i)   = 0._RKIND
    enddo
 endif

 if(config_convection_scheme .ne. 'off') then
    call iamas_pool_get_array(diag_physics,'raincv',raincv)
    do i = its,ite
       rainbl(i) = rainbl(i) + raincv(i)
    enddo
 endif
!Changed by qkyang to speed iAMAS model---------------------

!call to land-surface scheme:

 lsm_select: select case (trim(lsm_scheme))
    case("noah")
       call lsm( &
                dz8w        = dz_p        , p8w3d       = pres2_hyd_p  , t3d       = t_p          , &
                qv3d        = qv_p        , xland       = xland(its:ite),xice      = xice(its:ite), &
                ivgtyp      = ivgtyp(its:ite),isltyp    = isltyp(its:ite),tmn      = tmn(its:ite) , &
                vegfra      = vegfra(its:ite),shdmin    = shdmin(its:ite),shdmax   = shdmax(its:ite), &
                snoalb      = snoalb(its:ite),glw       = glw(its:ite)  , gsw      = gsw(its:ite), &
                swdown      = swdown(its:ite),rainbl    = rainbl(its:ite),embck    = sfc_emibck(its:ite), &
                sr          = sr(its:ite) , qgh         = qgh(its:ite), cpm       = cpm(its:ite) , &
                qz0         = qz0(its:ite), tsk         = skintemp(its:ite),hfx    = hfx(its:ite) , &
                qfx         = qfx(its:ite), lh          = lh(its:ite)  , grdflx    = grdflx(its:ite), &
                qsfc        = qsfc(its:ite),cqs2        = cqs2(its:ite), chs       = chs(its:ite) , &
                chs2        = chs2(its:ite),snow        = snow(its:ite), snowc     = snowc(its:ite), &
                snowh       = snowh(its:ite),canwat     = canwat(its:ite),smstav   = smstav(its:ite), &
                smstot      = smstot(its:ite),sfcrunoff = sfcrunoff(its:ite),udrunoff=udrunoff(its:ite), &
                acsnom      = acsnom(its:ite),acsnow    = acsnow(its:ite), snotime = snotime(its:ite), &
                snopcx      = snopcx(its:ite),emiss     = sfc_emiss(its:ite),rib   = br(its:ite) , &
                potevp      = potevp(its:ite),albedo    = sfc_albedo(its:ite),albbck= sfc_albbck(its:ite), &
                z0          = z0(its:ite) , znt         = znt(its:ite) , lai       = lai(its:ite) , &
                noahres     = noahres(its:ite), chklowq = chklowq(its:ite),sh2o    = sh2o(1:num_soils,its:ite), &
                smois       = smois(1:num_soils,its:ite),tslb= tslb(1:num_soils,its:ite),smcrel=smcrel(1:num_soils,its:ite), &
                dzs         = dzs(1:num_soils,its), isurban= isurban      , isice  = isice        , &
                iswater     = iswater, & !qkyang added "iswater" to support LCZ landuse dataset
                rovcp       = rcp         , dt          = dt_pbl       , myj       = myj          , &
                itimestep   = itimestep   , &
                !changed by qkyang to support urban model+++++++++++++++++++
! #IF ((config_urban_scheme.eq.1).OR.(config_urban_scheme.EQ.2).OR.(config_urban_scheme.EQ.3))
                U_PHY       = u_p         , V_PHY       = v_p          , th_phy   = th_p          , &
                !qkyang added to support BEP BEM model>>>>>>>>>>>>>>>>>
                rho         = rho_p       , P_PHY      =pres_p         , &
                TSK_RURAL_BEP =TSK_RURAL_BEP(its:ite), &
                !qkyang added to support BEP BEM model<<<<<<<<<<<<<<<<
                gmt         = gmt         , & 
                xlat        = xlat(its:ite), XLAT_URB2D=xlat(its:ite),xlong=xlong(its:ite), &
                DZR         = DZR         , DZB         = DZB          , DZG       = DZG          , &
                TRL_URB3D   = TRL_URB3D(1:num_soils,its:ite),TBL_URB3D=TBL_URB3D(1:num_soils,its:ite), TGL_URB3D = TGL_URB3D(1:num_soils,its:ite), &
                TRB_URB4D   = TRB_URB4D(1:num_urban_layers,its:ite),TW1_URB4D=TW1_URB4D(1:num_urban_layers,its:ite), TW2_URB4D = TW2_URB4D(1:num_urban_layers,its:ite), &
                TGB_URB4D   = TGB_URB4D(1:num_urban_layers,its:ite),TLEV_URB3D=TLEV_URB3D(1:num_urban_layers,its:ite),QLEV_URB3D=QLEV_URB3D(1:num_urban_layers,its:ite), &
                TW1LEV_URB3D=TW1LEV_URB3D(1:num_urban_layers,its:ite),TW2LEV_URB3D=TW2LEV_URB3D(1:num_urban_layers,its:ite),TGLEV_URB3D=TGLEV_URB3D(1:num_urban_layers,its:ite), &
                TFLEV_URB3D =TFLEV_URB3D(1:num_urban_layers,its:ite),SFWIN1_URB3D=SFWIN1_URB3D(1:num_urban_layers,its:ite),SFWIN2_URB3D=SFWIN2_URB3D(1:num_urban_layers,its:ite),&
                SFW1_URB3D  = SFW1_URB3D(1:num_urban_layers,its:ite),SFW2_URB3D=SFW2_URB3D(1:num_urban_layers,its:ite), SFR_URB3D = SFR_URB3D(1:num_urban_layers,its:ite), &
                SFG_URB3D   = SFG_URB3D(1:num_urban_layers,its:ite), TGRL_URB3D=TGRL_URB3D(1:num_soils,its:ite), SMR_URB3D = SMR_URB3D(1:num_soils,its:ite), &
                A_U_BEP     = A_U_BEP(kms:kme,its:ite), A_V_BEP=A_V_BEP(kms:kme,its:ite),A_T_BEP= A_T_BEP(kms:kme,its:ite), &
                A_Q_BEP     = A_Q_BEP(kms:kme,its:ite), A_E_BEP=A_E_BEP(kms:kme,its:ite),B_U_BEP= B_U_BEP(kms:kme,its:ite), &
                B_V_BEP     = B_V_BEP(kms:kme,its:ite), B_T_BEP=B_T_BEP(kms:kme,its:ite),B_Q_BEP= B_Q_BEP(kms:kme,its:ite), &
                B_E_BEP     = B_E_BEP(kms:kme,its:ite), DLG_BEP=DLG_BEP(kms:kme,its:ite),DL_U_BEP=DL_U_BEP(kms:kme,its:ite), &
                SF_BEP      = SF_BEP(kms:kme,its:ite), VL_BEP= VL_BEP(kms:kme,its:ite), XXXR_URB2D= XXXR_URB2D(its:ite), &
                XXXB_URB2D  = XXXB_URB2D(its:ite), XXXG_URB2D=XXXG_URB2D(its:ite), XXXC_URB2D= XXXC_URB2D(its:ite), &
                TR_URB2D    = TR_URB2D(its:ite), TB_URB2D=TB_URB2D(its:ite), TG_URB2D= TG_URB2D(its:ite), &
                TC_URB2D    = TC_URB2D(its:ite), QC_URB2D=QC_URB2D(its:ite), UC_URB2D= UC_URB2D(its:ite), &
                SH_URB2D    = SH_URB2D(its:ite), &
                LH_URB2D    = LH_URB2D(its:ite),G_URB2D = G_URB2D(its:ite),RN_URB2D = RN_URB2D(its:ite), &
                TS_URB2D    = TS_URB2D(its:ite),SF_AC_URB3D=SF_AC_URB3D(its:ite),LF_AC_URB3D=LF_AC_URB3D(its:ite), &
                CM_AC_URB3D = CM_AC_URB3D(its:ite),SFVENT_URB3D =SFVENT_URB3D(its:ite),LFVENT_URB3D=LFVENT_URB3D(its:ite),&
                LP_URB2D    = LP_URB2D(its:ite),HI_URB2D=HI_URB2D(1:num_urban_hi,its:ite), LB_URB2D=LB_URB2D(its:ite), &
                HGT_URB2D   = HGT_URB2D(its:ite),&
                MH_URB2D    = MH_URB2D(its:ite),STDH_URB2D=STDH_URB2D(its:ite),LF_URB2D=LF_URB2D(1:num_soils,its:ite), &
                CMCR_URB2D  = CMCR_URB2D(its:ite),TGR_URB2D=TGR_URB2D(its:ite),DRELR_URB2D=DRELR_URB2D(its:ite), &
                DRELB_URB2D = DRELB_URB2D(its:ite),DRELG_URB2D =DRELG_URB2D(its:ite),FLXHUMR_URB2D=FLXHUMR_URB2D(its:ite),&
                FLXHUMB_URB2D=FLXHUMB_URB2D(its:ite),FLXHUMG_URB2D=FLXHUMG_URB2D(its:ite),&
                PSIM_URB2D  = PSIM_URB2D(its:ite),PSIH_URB2D= PSIH_URB2D(its:ite),GZ1OZ0_URB2D=GZ1OZ0_URB2D(its:ite), &
                U10_URB2D   = U10_URB2D(its:ite), V10_URB2D= V10_URB2D(its:ite), TH2_URB2D=TH2_URB2D(its:ite), &
                Q2_URB2D    = Q2_URB2D(its:ite), AKMS_URB2D= AKMS_URB2D(its:ite), &
                CHR_SFCDIF  = CHR_SFCDIF(its:ite),CMR_SFCDIF=CMR_SFCDIF(its:ite),CMGR_SFCDIF=CMGR_SFCDIF(its:ite), &
                CHGR_SFCDIF = CHGR_SFCDIF(its:ite),CMC_SFCDIF=CMC_SFCDIF(its:ite),CHC_SFCDIF=CHC_SFCDIF(its:ite), &
                ust         = ust(its:ite), DECLIN_URB  = declin       , COSZ_URB2D = coszr(its:ite), &
                OMG_URB2D   = hrang(its:ite), &
                config_urban_scheme = config_urban_scheme ,  &!qkyang replaced sf_urban_physics with config_urban_scheme
                !changed by qkyang to support urban model--------------
                !qkyang added to support anthropogenic heating emissions++++++
                AH_emis     = AH_emis(its:ite), &
                is_ah_emis  = is_ah_emis, &
                AH_URB2D    = AH_URB2D, &
                !qkyang added to support anthropogenic heating emissions------
                frpcpn      = frpcpn      , rdlai2d     = rdlai2d      , &
                opt_thcnd   = opt_thcnd   , ua_phys     = ua_phys      , flx4_2d   = flxsnow(its:ite), &
                fvb_2d      = fvbsnow(its:ite), fbur_2d = fbursnow(its:ite), fgsn_2d= fgsnsnow(its:ite), &
                utype_urb2d = utype_urb2d(its:ite), frc_urb2d = frc_urb2d(its:ite), ust_urb2d = ust_urb2d(its:ite), &
                !qkyang replaced utype_urb_p,frc_urb_p,ust_urb_p with utype_urb2d_p, frc_urb2d_p,ust_urb2d_p
                infxsrt     = infxsrt(its:ite), sfcheadrt= sfcheadrt(its:ite), soldrain  = soldrain(its:ite), &
                fasdas      = fasdas      , & !qkyang comment julian      = 0            , julyr     = 0            , &
                julday      = julday      , julian      = curr_julday  , julyr     = year         , &
                !qkyang added julday      = julday      , julian      = curr_julday  , julyr     = year
                xice_threshold   = xice_threshold    ,                                              &
                usemonalb        = config_sfc_albedo ,                                              &
                mminlu           = mminlu ,                                                         &
                num_soil_layers  = num_soils ,                                                      &
                num_roof_layers  = num_soils ,                                                      &
                num_wall_layers  = num_soils ,                                                      &
                num_road_layers  = num_soils ,                                                      &
                !qkyang changed to support BEP BEM model>>>>>>>>>>>>>>>>>                                                     &
                num_urban_layers = num_urban_layers ,                                               &
                num_urban_hi     = num_urban_hi ,                                                   &
                !qkyang changed to support BEP BEM model>>>>>>>>>>>>>>>>>                                                    &
!                 sf_urban_physics = sf_urban_physics ,                                               &qkyang comment
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,             &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,             &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte               &
               )

       if(config_frac_seaice) then
          call seaice_noah( &
                dz8w      = dz_p        , p8w3d   = pres2_hyd_p  , t3d      = t_p      , &
                qv3d      = qv_p        , xice    = xice(its:ite), snoalb2d = snoalb(its:ite), &
                glw       = glw(its:ite), swdown  = swdown(its:ite),rainbl  = rainbl(its:ite), &
                sr        = sr(its:ite) , qgh     = qgh(its:ite), tsk      = skintemp(its:ite), &
                hfx       = hfx(its:ite), qfx     = qfx(its:ite) , lh       = lh(its:ite), &
                grdflx    = grdflx(its:ite), potevp  = potevp(its:ite), qsfc= qsfc(its:ite), &
                emiss     = sfc_emiss(its:ite), albedo=sfc_albedo(its:ite),rib=br(its:ite), &
                cqs2      = cqs2(its:ite),chs     = chs(its:ite) , chs2     = chs2(its:ite), &
                z02d      = z0(its:ite) , znt     = znt(its:ite) , tslb     = tslb(1:num_soils,its:ite), &
                snow      = snow(its:ite), snowc  = snowc(its:ite), snowh2d = snowh(its:ite), &
                snopcx    = snopcx(its:ite),acsnow= acsnow(its:ite), acsnom = acsnom(its:ite), &
                sfcrunoff = sfcrunoff(its:ite), albsi   = sfc_albedo_seaice(its:ite), snowsi   = snowsi(its:ite), &
                icedepth  = icedepth(its:ite), noahres = noahres(its:ite), dt       = dt_pbl   , &
                frpcpn    = frpcpn      ,                                                &
                seaice_albedo_opt        = seaice_albedo_opt        ,                    &
                seaice_albedo_default    = seaice_albedo_default    ,                    &
                seaice_thickness_opt     = seaice_thickness_opt     ,                    &
                seaice_thickness_default = seaice_thickness_default ,                    &
                seaice_snowdepth_opt     = seaice_snowdepth_opt     ,                    &
                seaice_snowdepth_max     = seaice_snowdepth_max     ,                    &
                seaice_snowdepth_min     = seaice_snowdepth_min     ,                    &
                xice_threshold           = xice_threshold           ,                    &
                num_soil_layers          = num_soils                ,                    &
                sf_urban_physics         = config_urban_scheme      ,                    &
                !qkyang added to support BEP BEM model>>>>>>>>>>>>>>>>>
                B_T_BEP=B_T_BEP(kms:kme,its:ite),B_Q_BEP= B_Q_BEP(kms:kme,its:ite),rho= rho_p, &
                !qkyang added to support BEP BEM model<<<<<<<<<<<<<<<<
                !qkyang replace sf_urban_physics with config_urban_scheme
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,  &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,  &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                          )
       endif
       call sfcdiags( &
                hfx   = hfx(its:ite), qfx= qfx(its:ite), tsk  =skintemp(its:ite), qsfc = qsfc(its:ite), chs = chs(its:ite), &
                chs2  = chs2(its:ite),cqs2=cqs2(its:ite), t2  = t2m(its:ite), th2  = th2m(its:ite), q2  = q2(its:ite), &
                psfc  = psfc_p , t3d     = t_p     , qv3d = qv_p  , cp   = cp     , R_d = R_d   , &
                rovcp = rcp    , ua_phys = ua_phys                                              , &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,           &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,           &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                    )
    case default

 end select lsm_select

!Changed by qkyang to speed iAMAS model+++++++++++++++++++++
if(config_frac_seaice) then
    do i = its,ite
      do j = jts,jte!
       if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
          chs(i)        = xice(i)*chs(i)       +(1._RKIND-xice(i))*chs_sea(i,j)
          chs2(i)       = xice(i)*chs2(i)      +(1._RKIND-xice(i))*chs2_sea(i,j)
          cqs2(i)       = xice(i)*cqs2(i)      +(1._RKIND-xice(i))*cqs2_sea(i,j)
          cpm(i)        = xice(i)*cpm(i)       +(1._RKIND-xice(i))*cpm_sea(i,j)
          hfx(i)        = xice(i)*hfx(i)       +(1._RKIND-xice(i))*hfx_sea(i,j)
          lh(i)         = xice(i)*lh(i)        +(1._RKIND-xice(i))*lh_sea(i,j)
          qfx(i)        = xice(i)*qfx(i)       +(1._RKIND-xice(i))*qfx_sea(i,j)
          qgh(i)        = xice(i)*qgh(i)       +(1._RKIND-xice(i))*qgh_sea(i,j)
          qsfc(i)       = xice(i)*qsfc(i)      +(1._RKIND-xice(i))*qsfc_sea(i,j)
        skintemp(i)     = xice(i)*skintemp(i)  +(1._RKIND-xice(i))*tsk_sea(i)
!           skintemp(i)   = xice(i)*tsk_p(i,j)   +(1._RKIND-xice(i))*tsk_sea(i,j)
          sfc_albedo(i) = xice(i)*sfc_albedo(i)+(1._RKIND-xice(i))*0.08_RKIND
          sfc_emiss(i)  = xice(i)*sfc_emiss(i) +(1._RKIND-xice(i))*0.98_RKIND
       endif
       enddo
    enddo
 endif

!added by qkyang to support urban model+++++++++++++++++++++
! output variable lsm --> surface_driver !Based on WRF
 !urban

     IF(config_urban_scheme.eq.1) THEN
         do i = its,ite!DO i=i_start(ij),i_end(ij)              !urban
          !qkyang replaced to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

             U10(I)  = U10_URB2D(I)      !U10(I,J)  = U10_URB2D(I,J)      !urban
             V10(I)  = V10_URB2D(I)      !V10(I,J)  = V10_URB2D(I,J)      !urban
             PSIM(I) = PSIM_URB2D(I)     !PSIM(I,J) = PSIM_URB2D(I,J)     !urban
             PSIH(I) = PSIH_URB2D(I)     !PSIH(I,J) = PSIH_URB2D(I,J)     !urban
             GZ1OZ0(I) = GZ1OZ0_URB2D(I) !GZ1OZ0(I,J) = GZ1OZ0_URB2D(I,J) !urban
!m             AKHS(I) = AKHS_URB2D_p(I,J)   !!AKHS(I,J) = AKHS_URB2D(I,J)    !urban
!-- ...(AKHS and AKMS generated by urban model, probably not be used since no MYJ, qnse yet):
             !akhs          sfc exchange coefficient of heat/moisture
             AKHS(I) = CHS(I)                !AKHS(I,J) = CHS(I,J)             !urban
             !akms          sfc exchange coefficient of momentum
             AKMS(I) = AKMS_URB2D(I)     !AKMS(I,J) = AKMS_URB2D(I,J)     !urban
           END IF                                               !urban
         ENDDO                                                  !urban
     ENDIF
! urban BEP
     IF((config_urban_scheme.eq.2).OR.(config_urban_scheme.eq.3)) THEN
        do i = its,ite!DO i=i_start(ij),i_end(ij)              !urban
          do j = jts,jte!DO j=j_start(ij),j_end(ij)              !urban
          !qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF(IVGTYP(I) == ISURBAN .or. &
             IVGTYP(I) == low_density_residential .or. &  
             IVGTYP(I) == HIGH_DENSITY_RESIDENTIAL .or. &
             IVGTYP(I) == HIGH_INTENSITY_INDUSTRIAL .or. &
             IVGTYP(I) == LCZ_1 .or. IVGTYP(I) == LCZ_2 .or. &
             IVGTYP(I) == LCZ_3 .or. IVGTYP(I) == LCZ_4 .or. &
             IVGTYP(I) == LCZ_5 .or. IVGTYP(I) == LCZ_6 .or. &
             IVGTYP(I) == LCZ_7 .or. IVGTYP(I) == LCZ_8 .or. &
             IVGTYP(I) == LCZ_9 .or. IVGTYP(I) == LCZ_10 .or. &
             IVGTYP(I) == LCZ_11 ) then !urban
             !qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            T2M(I)   = th_p(i,1,j)/((1.E5/psfc_p(I,J))**RCP)!TH_PHY(i,1,j)/((1.E5/PSFC(I,J))**RCP)   !urban
            TH2M(I)  = th_p(i,1,j)!TH_PHY(i,1,j)                 !urban
            Q2(I)   = qv_p(i,1,j)!qv_curr(i,1,j)                !urban
            U10(I)  = u_p(I,1,J)!U_phy(I,1,J)                   !urban
            V10(I)  = v_p(I,1,J)!V_phy(I,1,J)                   !urban
           END IF                                               !urban
         ENDDO                                                  !urban
        ENDDO                                                   !urban
     ENDIF
!added by qkyang to support urban model---------------------


!Changed by qkyang to speed iAMAS model---------------------

!copy local arrays to iAMAS grid:
! call lsm_to_iAMAS(configs,mesh,diag_physics,sfc_input,its,ite)!qkyang comment 

!Changed by qkyang to speed iAMAS model+++++++++++++++++++++
 if(config_frac_seaice) then
    if(allocated(chs_sea)   ) deallocate(chs_sea   )
    if(allocated(chs2_sea)  ) deallocate(chs2_sea  )
    if(allocated(cqs2_sea)  ) deallocate(cqs2_sea  )
    if(allocated(cpm_sea)   ) deallocate(cpm_sea   )
    if(allocated(hfx_sea)   ) deallocate(hfx_sea   )
    if(allocated(qfx_sea)   ) deallocate(qfx_sea   )
    if(allocated(qgh_sea)   ) deallocate(qgh_sea   )
    if(allocated(qsfc_sea)  ) deallocate(qsfc_sea  )
    if(allocated(lh_sea)    ) deallocate(lh_sea    )
 endif
!Changed by qkyang to speed iAMAS model---------------------
 end subroutine driver_lsm

!Changed by qkyang to support urban model+++++++++++++++++++++
!=================================================================================================================
 subroutine correct_tsk_over_seaice(ims,ime,jms,jme,its,ite,jts,jte,xice_thresh,xice,tsk,tsk_sea,tsk_ice)
!=================================================================================================================

!input arguments:
 integer,intent(in):: ims,ime,its,ite,jms,jme,jts,jte
 real(kind=RKIND),intent(in):: xice_thresh
 real(kind=RKIND),intent(in),dimension(its:ite):: tsk,xice

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite):: tsk_sea,tsk_ice

!local variables:
 integer:: i,j

!-----------------------------------------------------------------------------------------------------------------

!initialize the local sea-surface temperature and local sea-ice temperature to the local surface
!temperature:
 do i = its,ite
!    do j = jts,jte
    tsk_sea(i) = tsk(i)
    tsk_ice(i) = tsk(i)

    if(xice(i).ge.xice_thresh .and. xice(i).le.1._RKIND) then
       !over sea-ice grid cells, limit sea-surface temperatures to temperatures warmer than 271.4:
       tsk_sea(i) = max(tsk_sea(i),271.4_RKIND)

       !over sea-ice grid cells, avoids unphysically too cold sea-ice temperatures for grid cells
       !with small sea-ice fractions:
       if(xice(i).lt.0.2_RKIND .and. tsk_ice(i).lt.253.15_RKIND) tsk_ice(i) = 253.15_RKIND
       if(xice(i).lt.0.1_RKIND .and. tsk_ice(i).lt.263.15_RKIND) tsk_ice(i) = 263.15_RKIND
    endif
!    enddo
 enddo

 end subroutine correct_tsk_over_seaice
!Changed by qkyang to support urban model---------------------
!=================================================================================================================
 end module iamas_atmphys_driver_lsm
!=================================================================================================================
