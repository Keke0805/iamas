! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
#define DM_BCAST_CHAR(A)     call iamas_dmpar_bcast_char(dminfo,A)
#define DM_BCAST_INTEGER(A)  call iamas_dmpar_bcast_int(dminfo,A)
#define DM_BCAST_INTEGERS(A) call iamas_dmpar_bcast_ints(dminfo,size(A),A)
#define DM_BCAST_REAL(A)     call iamas_dmpar_bcast_real(dminfo,A)
#define DM_BCAST_REALS(A)    call iamas_dmpar_bcast_reals(dminfo,size(A),A)

!=================================================================================================================
 module iamas_atmphys_lsm_noahinit
 use iamas_dmpar
 use iamas_kind_types
 use iamas_pool_routines

 use iamas_atmphys_constants
 use iamas_atmphys_utilities
!wrf physics
 use module_sf_noahlsm

 implicit none
 private
 public:: noah_init_foriAMAS


!Initialization of the NOAH 4-layer land surface scheme.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_lsm_noahinit:
! -----------------------------------------
! noah_init_foriAMAS: call lsminit from subroutine init_lsm (module iamas_atmphyse_driver_lsm.F).
! lsminit          : main initialization subroutine for the NOAH 4-layer land surface scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * replaced the variable grav with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * added "use iamas_kind_types" at the top of the module.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * in subroutine soil_veg_gen_parm, modified reading the updated file VEGPARM.TBL so that we can update the NOAH
!   land surface scheme.added the categories low_density_residential,high_density_residential,and high_intensity_
!   industrial.added the variables ztopvtbl and zbotvtbl.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-25.

 contains


!=================================================================================================================
 subroutine noah_init_foriAMAS(dminfo,mesh,configs,diag_physics,sfc_input)
!=================================================================================================================
!input arguments:
 type(dm_info):: dminfo
 type(iamas_pool_type):: mesh
 type(iamas_pool_type),intent(in):: configs
!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: sfc_input

!-----------------------------------------------------------------------------------------------------------------

!read formatted files needed for land-surface model:
 call lsminit(dminfo,mesh,configs,diag_physics,sfc_input)

 end subroutine noah_init_foriAMAS

!=================================================================================================================
 subroutine lsminit(dminfo,mesh,configs,diag_physics,sfc_input)
!=================================================================================================================
 !added by qkyang to support urban model++++++++++++
  USE module_sf_urban
  use iamas_atmphys_landuse, only: isurban
!   USE module_sf_bep
!   USE module_sf_bep_bem
!qkyang changed to support BEP BEM model>>>>>>>>>>>>>>>>>     
 character(len=StrKIND),pointer:: pbl_scheme
 logical,pointer:: flag_bep 
 integer,pointer:: num_urban_layers,num_urban_hi,nVertLevelsP1
!qkyang changed to support BEP BEM model<<<<<<<<<<<<<<<<<<
 integer,pointer:: config_urban_scheme,nVertLevels
 real(kind=RKIND),dimension(:),pointer:: DZR,DZB,DZG
 real(kind=RKIND),dimension(:),pointer:: skintemp,tmn
 real(kind=RKIND),dimension(:),pointer:: XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D
 real(kind=RKIND),dimension(:),pointer:: TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,QC_URB2D
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
 real(kind=RKIND),dimension(:),pointer:: FRC_URB2D
 integer,dimension(:),pointer:: UTYPE_URB2D
 !added by qkyang to support urban model------------

!qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
logical,pointer:: config_use_lcz
!qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: sfc_input

!local pointers::
 logical,pointer:: input_sfc_albedo,restart
 character(len=StrKIND),pointer:: mminlu,mminsl
 integer,pointer:: nCells,nSoilLevels
 integer,dimension(:),pointer:: ivgtyp,isltyp
 real(kind=RKIND),dimension(:),pointer:: snoalb,snow,snowh
 real(kind=RKIND),dimension(:,:),pointer:: tslb,smois,sh2o

!local variables:
 logical,parameter:: fndsnowh = .true.
 integer:: iCell
 integer:: errflag,ns
 real(kind=RKIND):: bx,fk,smcmax,psisat,free
 real(kind=RKIND),parameter:: blim = 5.5
 real(kind=RKIND),parameter:: hlice = 3.335e5
 real(kind=RKIND),parameter:: t0 = 273.15

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_array(sfc_input,'mminlu'              ,mminlu          )
 call iamas_pool_get_config(configs,'input_soil_data'      ,mminsl          )
 call iamas_pool_get_config(configs,'config_sfc_snowalbedo',input_sfc_albedo)
 call iamas_pool_get_config(configs,'config_do_restart'    ,restart         )
 call iamas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call iamas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call iamas_pool_get_array(sfc_input,'isltyp', isltyp)
 call iamas_pool_get_array(sfc_input,'ivgtyp', ivgtyp)
 call iamas_pool_get_array(sfc_input,'sh2o'  , sh2o  )
 call iamas_pool_get_array(sfc_input,'smois' , smois )
 call iamas_pool_get_array(sfc_input,'tslb'  , tslb  )
 call iamas_pool_get_array(sfc_input,'snoalb', snoalb)
 call iamas_pool_get_array(sfc_input,'snow'  , snow  )
 call iamas_pool_get_array(sfc_input,'snowh' , snowh )

!qkyang added config_use_lcz to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
call iamas_pool_get_config(configs, 'config_use_lcz', config_use_lcz)
write(err_message,*) "iamas_atmphys_lsm_noahinit.F: config_use_lcz= ", config_use_lcz
call physics_message(err_message)
!qkyang added config_use_lcz to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!added by qkyang to support urban model++++++++++++
 call iamas_pool_get_config(configs,'config_urban_scheme',config_urban_scheme)
  !qkyang added to support BEP BEM model>>>>>>>>>>>>>>>>>
 call iamas_pool_get_config(configs,'num_urban_layers',num_urban_layers)
 call iamas_pool_get_config(configs,'num_urban_hi',num_urban_hi)
 call iamas_pool_get_dimension(mesh,'nVertLevelsP1',nVertLevelsP1)
 !qkyang added to support BEP BEM model<<<<<<<<<<<<<<<<<<
 call iamas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
  !qkyang added to support BEP and BEM models>>>>>>>>>>>>
 call iamas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)
 call iamas_pool_get_config(configs,'flag_bep' , flag_bep )
 !qkyang added to support BEP and BEM models<<<<<<<<<<<<
 call iamas_pool_get_array(diag_physics,'DZR' , DZR )
 call iamas_pool_get_array(diag_physics,'DZB' , DZB )
 call iamas_pool_get_array(diag_physics,'DZG' , DZG )
 call iamas_pool_get_array(diag_physics,'XXXR_URB2D' , XXXR_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXB_URB2D' , XXXB_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXG_URB2D' , XXXG_URB2D )
 call iamas_pool_get_array(diag_physics,'XXXC_URB2D' , XXXC_URB2D )

 call iamas_pool_get_array(diag_physics,'TR_URB2D' , TR_URB2D )
 call iamas_pool_get_array(diag_physics,'TB_URB2D' , TB_URB2D )
 call iamas_pool_get_array(diag_physics,'TG_URB2D' , TG_URB2D )
 call iamas_pool_get_array(diag_physics,'TC_URB2D' , TC_URB2D )
 call iamas_pool_get_array(diag_physics,'QC_URB2D' , QC_URB2D )

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
 call iamas_pool_get_array(diag_physics,'UTYPE_URB2D' , UTYPE_URB2D )

 call iamas_pool_get_array(sfc_input,'skintemp' , skintemp )
 call iamas_pool_get_array(sfc_input,'tmn' , tmn )
 !added by qkyang to support urban model------------

 !reads the NOAH LSM tables:
 call iamas_log_write('')
 call physics_message('--- initialize NOAH LSM tables' )
 call soil_veg_gen_parm(dminfo,mminlu,mminsl)
 
 !qkyang added config_use_lcz to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 write(err_message,*) "iamas_atmphys_lsm_noahinit.F: config_urban_scheme= ", &
                               config_urban_scheme
 call physics_message(err_message)
 !qkyang added to support BEP and BEM models>>>>>>>>>>>>>>>
 if ((config_urban_scheme.EQ.2).OR.(config_urban_scheme.EQ.3)) THEN
      flag_bep=.true.
      if (trim(pbl_scheme).NE.'bl_ysu') then
          call physics_message('iamas_atmphys_lsm_noahinit.F: config_pbl_scheme=' // trim ( pbl_scheme ))
          call physics_message('Please use config_pbl_scheme=ysu with BEP/BEM urban scheme')
          call physics_error_fatal("Please use config_pbl_scheme=ysu with BEP/BEM urban scheme")
      endif
 else
      flag_bep=.false.
 endif
 !qkyang added to support BEP and BEM models<<<<<<<<<<<<<<<<<<<<<<
 IF ((config_urban_scheme.eq.1).OR.(config_urban_scheme.EQ.2).OR.(config_urban_scheme.EQ.3)) THEN
   call iamas_log_write('')

   call physics_message('--- start read URBPARM.TBL tables' )
   !qkyang added config_use_lcz to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   call urban_param_init(DZR,DZB,DZG,nSoilLevels,nCells,config_urban_scheme,config_use_lcz)
   !qkyang added config_use_lcz to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   call physics_message('--- end read URBPARM.TBL tables' )
   call physics_message('--- start call urban_var_init' )
   CALL urban_var_init(&!ISURBAN,TSK,TSLB,TMN,IVGTYP,                & !urban&
              ISURBAN,skintemp,TSLB,TMN,IVGTYP,                & !urban
              !ims,ime,jms,jme,kms,kme,nSoilLevels,         & !urban
              nCells,nVertLevels,nVertLevelsP1,nSoilLevels,    & !qkyang added nVertLevelsP1 to support BEP BEM model
              LOW_DENSITY_RESIDENTIAL,                        &
              HIGH_DENSITY_RESIDENTIAL,                       &
              HIGH_INTENSITY_INDUSTRIAL,                      &
              !qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              LCZ_1,                                          &
              LCZ_2,                                          &
              LCZ_3,                                          &
              LCZ_4,                                          &
              LCZ_5,                                          &
              LCZ_6,                                          &
              LCZ_7,                                          &
              LCZ_8,                                          &
              LCZ_9,                                          &
              LCZ_10,                                         &
              LCZ_11,                                         &
              config_use_lcz,                                 &
              !qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              !num_roof_layers,num_wall_layers,num_road_layers, & !urban
              !restart,sf_urban_physics,                        & !urban
              restart,config_urban_scheme,                     & !urban
              XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D,     & !urban
              TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,QC_URB2D,    & !urban
              TRL_URB3D,TBL_URB3D,TGL_URB3D,                   & !urban
              SH_URB2D,LH_URB2D,G_URB2D,RN_URB2D, TS_URB2D,    & !urban
              num_urban_layers,                                & !urban
              num_urban_hi,                                    & !urban
              TRB_URB4D,TW1_URB4D,TW2_URB4D,TGB_URB4D,         & !urban

              TLEV_URB3D,QLEV_URB3D,                           & !urban
              TW1LEV_URB3D,TW2LEV_URB3D,                       & !urban
              TGLEV_URB3D,TFLEV_URB3D,                         & !urban

              SF_AC_URB3D,LF_AC_URB3D,CM_AC_URB3D,             & !urban
              SFVENT_URB3D,LFVENT_URB3D,                       & !urban
              SFWIN1_URB3D,SFWIN2_URB3D,                       & !urban

              SFW1_URB3D,SFW2_URB3D,SFR_URB3D,SFG_URB3D,       & !urban
              LP_URB2D,HI_URB2D,LB_URB2D,                      & !urban
              HGT_URB2D,MH_URB2D,STDH_URB2D,LF_URB2D,          & !urban

              CMCR_URB2D,TGR_URB2D,TGRL_URB3D,SMR_URB3D,       & ! inout
              DRELR_URB2D,DRELB_URB2D,DRELG_URB2D,             & ! inout
              FLXHUMR_URB2D, FLXHUMB_URB2D, FLXHUMG_URB2D,     & ! inout

              A_U_BEP,A_V_BEP,A_T_BEP,A_Q_BEP,                 & !multi-layer urban
              A_E_BEP,B_U_BEP,B_V_BEP,                         & !multi-layer urban
              B_T_BEP,B_Q_BEP,B_E_BEP,DLG_BEP,                 & !multi-layer urban
              DL_U_BEP,SF_BEP,VL_BEP,                          & !multi-layer urban
              FRC_URB2D, UTYPE_URB2D)                            !urban
   call physics_message('--- end call urban_var_init' )
 ENDIF
 !added by qkyang to support urban model------------

 call physics_message('--- end initialize NOAH LSM tables' )
 call iamas_log_write('')

 if(.not.restart) then

    errflag = 0
    do iCell = 1, nCells
       if(isltyp(iCell) .lt. 1) then
          errflag = 1
          write(err_message,*) "module_sf_noahlsm.F: lsminit: out of range ISLTYP ", &
                               iCell,isltyp(iCell)
          call physics_message(err_message)
       endif
       if(.not. input_sfc_albedo)  then
           IF(ivgtyp(iCell) == ISURBAN .or. &
              ivgtyp(iCell) == low_density_residential .or. &  
              ivgtyp(iCell) == HIGH_DENSITY_RESIDENTIAL .or. &
              ivgtyp(iCell) == HIGH_INTENSITY_INDUSTRIAL .or. &
              ivgtyp(iCell) == LCZ_1 .or. ivgtyp(iCell) == LCZ_2 .or. &
              ivgtyp(iCell) == LCZ_3 .or. ivgtyp(iCell) == LCZ_4 .or. &
              ivgtyp(iCell) == LCZ_5 .or. ivgtyp(iCell) == LCZ_6 .or. &
              ivgtyp(iCell) == LCZ_7 .or. ivgtyp(iCell) == LCZ_8 .or. &
              ivgtyp(iCell) == LCZ_9 .or. ivgtyp(iCell) == LCZ_10 .or. &
              ivgtyp(iCell) == LCZ_11 )THEN !urban
             snoalb(iCell) = maxalb(ISURBAN)*0.01
           else
            snoalb(iCell) = maxalb(ivgtyp(iCell))*0.01 !default
           endif
        endif
    enddo
    if(errflag .eq. 1) &
       call physics_error_fatal("module_sf_noahlsm.F: lsminit: out of range value "// &
                            "of ISLTYP. Is this field in the input?" )

!initializes soil liquid water content SH2O:
    do iCell = 1, nCells
       bx = bb(isltyp(iCell))
       smcmax = maxsmc(isltyp(iCell))
       psisat = satpsi(isltyp(iCell))
       if((bx > 0.0).and.(smcmax > 0.0).and.(psisat > 0.0)) then
          do ns = 1, nSoilLevels
! ----------------------------------------------------------------------
!SH2O  <= SMOIS for T < 273.149K (-0.001C)
             if(tslb(ns,iCell) < 273.149) then
! ----------------------------------------------------------------------
! first guess following explicit solution for Flerchinger Eqn from Koren
! et al, JGR, 1999, Eqn 17 (KCOUNT=0 in FUNCTION FRH2O).
! ISLTPK is soil type
                bx = bb(isltyp(iCell))
                smcmax = maxsmc(isltyp(iCell))
                psisat = satpsi(isltyp(iCell))
                if(bx >  blim) bx = blim
                fk = (((hlice/(gravity*(-psisat))) * &
                     ((tslb(ns,iCell)-t0)/tslb(ns,iCell)) )**(-1/bx) )*smcmax
                if (fk < 0.02) fk = 0.02
                sh2o(ns,iCell) = min(fk,smois(ns,iCell))
! ----------------------------------------------------------------------
! now use iterative solution for liquid soil water content using
! FUNCTION FRH2O with the initial guess for SH2O from above explicit
! first guess.
                call frh2o(free,tslb(ns,iCell),smois(ns,iCell),sh2o(ns,iCell), &
                           smcmax,bx,psisat)
                sh2o(ns,iCell) = free
             else         ! of if (tslb(i,ns,j)
! ----------------------------------------------------------------------
! SH2O = SMOIS ( for T => 273.149K (-0.001C)
                sh2o(ns,iCell)=smois(ns,iCell)
! ----------------------------------------------------------------------
             endif        ! of if (tslb(i,ns,j)
          enddo
       else                  ! of if ((bx > 0.0)
          do ns = 1, nSoilLevels
             sh2o(ns,iCell)=smois(ns,iCell)
          enddo
       endif ! of if ((bx > 0.0)
    enddo    ! do iCell

!initialize physical snow height SNOWH:
    if(.not.fndsnowh)then
!if no snowh do the following:
       call physics_message( 'SNOW HEIGHT NOT FOUND - VALUE DEFINED IN LSMINIT' )
       do iCell = 1, nCells
          snowh(iCell)=snow(iCell)*0.005 ! snow in mm and snowh in m
       enddo
    endif
    110 continue

 endif

 end subroutine lsminit

!=================================================================================================================
 subroutine soil_veg_gen_parm(dminfo,mminlu,mminsl)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 character(len=*),intent(inout):: mminlu, mminsl

!local variables:
 character*128:: mess , message

 integer,parameter:: open_ok = 0
 integer:: lumatch,iindex,lc,num_slope
 integer:: istat

!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!ALBBCK: SFC albedo (in percentage)
!Z0    : Roughness length (m)
!SHDFAC: Green vegetation fraction (in percentage)

!Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!ALBEDO,and Z0 are specified in LAND-USE TABLE. SHDFAC is the monthly green vegetation data.
!CMXTBL: MAX CNPY Capacity (m)
!NROTBL: Rooting depth (layer)
!RSMIN : Mimimum stomatal resistance (s m-1)
!RSMAX : Max. stomatal resistance (s m-1)
!RGL   : Parameters used in radiation stress function
!HS    : Parameter used in vapor pressure deficit function
!TOPT  : Optimum transpiration air temperature. (K)
!CMCMAX: Maximum canopy water capacity
!CFACTR: Parameter used in the canopy inteception calculation
!SNUP  : Threshold snow depth (in water equivalent m) that implies 100% snow cover.
!LAI   : Leaf area index (dimensionless)
!MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETATION PROPERTIES FROM VEGPARM.TBL

!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('       enter subroutine soil_veg_gen_parm:')

!read in the vegetation properties from vegparm.tbl:

 if(dminfo % my_proc_id == IO_NODE) then
    open(16,file='VEGPARM.TBL',form='FORMATTED',status='OLD',iostat=istat)
    if(istat /= open_ok) &
       call physics_error_fatal('subroutine soil_veg_gen_arm: ' // &
                                'failure opening VEGPARM.TBL')

    lumatch=0
    find_lutype : do while (lumatch == 0)
       read(16,*,end=2002)
       read(16,*,end=2002) lutype
       read(16,*) lucats,iindex

       if(lutype.eq.trim(mminlu))then
          write(mess,*) 'iamas_atmphys_lsm_noahinit.F: landuse type = ' // trim ( lutype ) // ' found', &
                        lucats, ' categories'
          call physics_message(mess)
          lumatch=1
       else
          call physics_message('iamas_atmphys_lsm_noahinit.F: skipping over lutype = ' // trim ( lutype ))
          do lc = 1, lucats+20+22!qkyang +22 to support LCZ
             read(16,*)
          enddo
       endif
    enddo find_lutype

!prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008:
    if(size(shdtbl)       < lucats .or. &
       size(nrotbl)       < lucats .or. &
       size(rstbl)        < lucats .or. &
       size(rgltbl)       < lucats .or. &
       size(hstbl)        < lucats .or. &
       size(snuptbl)      < lucats .or. &
       size(maxalb)       < lucats .or. &
       size(laimintbl)    < lucats .or. &
       size(laimaxtbl)    < lucats .or. &
       size(z0mintbl)     < lucats .or. &
       size(z0maxtbl)     < lucats .or. &
       size(albedomintbl) < lucats .or. &
       size(albedomaxtbl) < lucats .or. &
       size(ztopvtbl)     < lucats .or. &
       size(zbotvtbl)     < lucats .or. &
       size(emissmintbl ) < lucats .or. &
       size(emissmaxtbl ) < lucats) then
       call physics_error_fatal('table sizes too small for value of lucats in module_sf_noahdrv.f')
    endif

    if(lutype.eq.mminlu)then
       do lc = 1, lucats
          read(16,*) iindex,shdtbl(lc),nrotbl(lc),rstbl(lc),rgltbl(lc),hstbl(lc),snuptbl(lc),  &
                     maxalb(lc),laimintbl(lc),laimaxtbl(lc),emissmintbl(lc),emissmaxtbl(lc),   &
                     albedomintbl(lc),albedomaxtbl(lc),z0mintbl(lc),z0maxtbl(lc),ztopvtbl(lc), &
                     zbotvtbl(lc)
       enddo
       read (16,*)
       read (16,*)topt_data
       read (16,*)
       read (16,*)cmcmax_data
       read (16,*)
       read (16,*)cfactr_data
       read (16,*)
       read (16,*)rsmax_data
       read (16,*)
       read (16,*)bare
       read (16,*)
       read (16,*)natural
       read (16,*)
       read (16,*)
       read (16,*)
       read (16,*)low_density_residential
       read (16,*)
       read (16,*)high_density_residential
       read (16,*)
       read (16,*)high_intensity_industrial
       !qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       read (16,*)
       read (16,*)LCZ_1
       read (16,*)
       read (16,*)LCZ_2
       read (16,*)
       read (16,*)LCZ_3
       read (16,*)
       read (16,*)LCZ_4
       read (16,*)
       read (16,*)LCZ_5
       read (16,*)
       read (16,*)LCZ_6
       read (16,*)
       read (16,*)LCZ_7
       read (16,*)
       read (16,*)LCZ_8
       read (16,*)
       read (16,*)LCZ_9
       read (16,*)
       read (16,*)LCZ_10
       read (16,*)
       read (16,*)LCZ_11
       !qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    endif

    2002 continue
    close (16)
    if(lumatch == 0) &
       call physics_error_fatal ('land use dataset '''//mminlu//''' not found in VEGPARM.TBL.')
      
 endif ! end dminfo

!distribute data to all processors: 
 DM_BCAST_CHAR(lutype)
 DM_BCAST_INTEGER(lucats)
 DM_BCAST_INTEGER(iindex)
 DM_BCAST_INTEGER(lumatch)
 DM_BCAST_REALS(shdtbl)
 DM_BCAST_INTEGERS(nrotbl)
 DM_BCAST_REALS(rstbl)
 DM_BCAST_REALS(rgltbl)
 DM_BCAST_REALS(hstbl)
 DM_BCAST_REALS(snuptbl)
 DM_BCAST_REALS(laimintbl)
 DM_BCAST_REALS(laimaxtbl)
 DM_BCAST_REALS(z0mintbl)
 DM_BCAST_REALS(z0maxtbl)
 DM_BCAST_REALS(emissmintbl)
 DM_BCAST_REALS(emissmaxtbl)
 DM_BCAST_REALS(albedomintbl)
 DM_BCAST_REALS(albedomaxtbl)
 DM_BCAST_REALS(ztopvtbl)
 DM_BCAST_REALS(zbotvtbl)
 DM_BCAST_REALS(maxalb)
 DM_BCAST_REAL(topt_data)
 DM_BCAST_REAL(cmcmax_data)
 DM_BCAST_REAL(cfactr_data)
 DM_BCAST_REAL(rsmax_data)
 DM_BCAST_INTEGER(bare)
 DM_BCAST_INTEGER(natural)
 DM_BCAST_INTEGER(low_density_residential)
 DM_BCAST_INTEGER(high_density_residential)
 DM_BCAST_INTEGER(high_intensity_industrial)
 !qkyang added to support LCZ landuse dataset>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 DM_BCAST_INTEGER(LCZ_1)
 DM_BCAST_INTEGER(LCZ_2)
 DM_BCAST_INTEGER(LCZ_3)
 DM_BCAST_INTEGER(LCZ_4)
 DM_BCAST_INTEGER(LCZ_5)
 DM_BCAST_INTEGER(LCZ_6)
 DM_BCAST_INTEGER(LCZ_7)
 DM_BCAST_INTEGER(LCZ_8)
 DM_BCAST_INTEGER(LCZ_9)
 DM_BCAST_INTEGER(LCZ_10)
 DM_BCAST_INTEGER(LCZ_11)
!qkyang added to support LCZ landuse dataset<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!call iamas_log_write(' LUTYPE  = '//trim(lutype))
!call iamas_log_write(' LUCATS  = $i',intArgs=(/lucats/))
!call iamas_log_write(' IINDEX  = $i',intArgs=(/iindex/))
!call iamas_log_write(' LUMATCH = $i',intArgs=(/lumatch/))

!call iamas_log_write(' TOPT_DATA   = $r',realArgs=(/topt_data/))
!call iamas_log_write(' CMCMAX_DATA = $r',realArgs=(/cmcmax_data/))
!call iamas_log_write(' CFACTR_DATA = $r',realArgs=(/cfactr_data/))
!call iamas_log_write(' RSMAX_DATA  = $r',realArgs=(/rsmax_data/))
!call iamas_log_write(' BARE        = $i',intArgs=(/bare/))
!call iamas_log_write(' NATURAL     = $i',intArgs=(/natural/))
call iamas_log_write(' LOW_DENSITY_RESIDENTIAL = $i', intArgs=(/low_density_residential/))
call iamas_log_write(' HIGH_DENSITY_RESIDENTIAL= $i', intArgs=(/high_density_residential/))
call iamas_log_write(' HIGH_INTENSITY_INDUSTRIAL = $i', intArgs=(/high_intensity_industrial/))
call iamas_log_write(' LCZ_1  = $i', intArgs=(/LCZ_1/))
call iamas_log_write(' LCZ_2  = $i', intArgs=(/LCZ_2/))
call iamas_log_write(' LCZ_3  = $i', intArgs=(/LCZ_3/))
call iamas_log_write(' LCZ_4  = $i', intArgs=(/LCZ_4/))
call iamas_log_write(' LCZ_5  = $i', intArgs=(/LCZ_5/))
call iamas_log_write(' LCZ_6  = $i', intArgs=(/LCZ_6/))
call iamas_log_write(' LCZ_7  = $i', intArgs=(/LCZ_7/))
call iamas_log_write(' LCZ_8  = $i', intArgs=(/LCZ_8/))
call iamas_log_write(' LCZ_9  = $i', intArgs=(/LCZ_9/))
call iamas_log_write(' LCZ_10 = $i', intArgs=(/LCZ_10/))
call iamas_log_write(' LCZ_11 = $i', intArgs=(/LCZ_11/))
call iamas_log_write('')
! do lc = 1, lucats
!   call iamas_log_write('$i $r $r $r $r $r $r $r $r $r $r $r $r $r $r $r $r $r', intArgs=(/lc/),     &
!                realArgs=(/shdtbl(lc),float(nrotbl(lc)),rstbl(lc),rgltbl(lc),hstbl(lc),snuptbl(lc), &
!                maxalb(lc),laimintbl(lc),laimaxtbl(lc),emissmintbl(lc),emissmaxtbl(lc),             &
!                albedomintbl(lc),albedomaxtbl(lc),z0mintbl(lc),z0maxtbl(lc),ztopvtbl(lc),           &
!                zbotvtbl(lc)/))
! enddo
 call iamas_log_write('    end read VEGPARM.TBL')

!read in soil properties from soilparm.tbl:

 if(dminfo % my_proc_id == IO_NODE) then
    open(16,file='SOILPARM.TBL',form='FORMATTED',status='OLD',iostat=istat)
    if(istat /= open_ok) &
       call physics_error_fatal('module_sf_noahlsm.F: soil_veg_gen_parm: ' // &
                                'failure opening SOILPARM.TBL' )

    write(mess,*) '   input soil texture classification = ', trim (mminsl)
    call physics_message(mess)

    lumatch=0
    read(16,*)
    read(16,2000,end=2003) sltype
    2000 format(a4)
    read(16,*)slcats,iindex
    if(sltype.eq.mminsl)then
       write(mess,*) '   soil texture classification = ', trim ( sltype ) , ' found', &
                  slcats,' categories'
       call physics_message ( mess )
       lumatch=1
    endif

!prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008:
    if(size(bb) < slcats .or. &
       size(drysmc) < slcats .or. &
       size(f11   ) < slcats .or. &
       size(maxsmc) < slcats .or. &
       size(refsmc) < slcats .or. &
       size(satpsi) < slcats .or. &
       size(satdk ) < slcats .or. &
       size(satdw ) < slcats .or. &
       size(wltsmc) < slcats .or. &
       size(qtz   ) < slcats) then
!      call wrf_error_fatal('table sizes too small for value of slcats in module_sf_noahdrv.f')
    endif
    if(sltype.eq.mminsl) then
       do lc = 1, slcats
          read(16,*) iindex,bb(lc),drysmc(lc),f11(lc),maxsmc(lc),refsmc(lc),satpsi(lc), &
                      satdk(lc),satdw(lc),wltsmc(lc),qtz(lc)
       enddo
    endif

    2003 continue
    close(16)
    if(lumatch.eq.0)then
       call physics_message( 'soil texture in input file does not ' )
       call physics_message( 'match soilparm table'                 )
       call physics_error_fatal( 'inconsistent or missing soilparm file' )
    endif

 endif

!distribute data to all processors: 
 DM_BCAST_INTEGER(lumatch)
 DM_BCAST_CHAR(sltype)
 DM_BCAST_CHAR(mminsl)
 DM_BCAST_INTEGER(slcats)
 DM_BCAST_INTEGER(iindex)
 DM_BCAST_REALS(bb)
 DM_BCAST_REALS(drysmc)
 DM_BCAST_REALS(f11)
 DM_BCAST_REALS(maxsmc)
 DM_BCAST_REALS(refsmc)
 DM_BCAST_REALS(satpsi)
 DM_BCAST_REALS(satdk)
 DM_BCAST_REALS(satdw)
 DM_BCAST_REALS(wltsmc)
 DM_BCAST_REALS(qtz)

!call iamas_log_write(' LUMATCH=$i',intArgs=(/lumatch/))
!call iamas_log_write(' SLTYPE ='//trim(sltype))
!call iamas_log_write(' MMINSL ='//mminsl)
!call iamas_log_write(' SLCATS =$i',intArgs=(/slcats/))
!call iamas_log_write(' IINDEX =$i',intArgs=(/iindex/))
!call iamas_log_write('')
!do lc = 1, slcats
!   call iamas_log_write('$i $r $r $r $r $r $r $r $r $r $r', intArgs=(/lc/), &
!                realArgs=(/bb(lc),drysmc(lc),f11(lc),maxsmc(lc),refsmc(lc),satpsi(lc), &
!                satdk(lc),satdw(lc),wltsmc(lc),qtz(lc)/))
!enddo
 call iamas_log_write('    end read SOILPARM.TBL')

!read in general parameters from genparm.tbl:

 if(dminfo % my_proc_id == IO_NODE) then
    open(16,file='GENPARM.TBL',form='FORMATTED',status='OLD',iostat=istat)
    if(istat /= open_ok) &
       call physics_error_fatal('module_sf_noahlsm.F: soil_veg_gen_parm: ' // &
                                'failure opening GENPARM.TBL' )
    read(16,*)
    read(16,*)
    read(16,*) num_slope

    slpcats=num_slope
!prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008:
    if(size(slope_data) < num_slope) &
       call physics_error_fatal('NUM_SLOPE too large for slope_data array' // &
                                'in module_sf_noahdrv')

    do lc = 1, slpcats
       read(16,*)slope_data(lc)
    enddo
    read(16,*)
    read(16,*)sbeta_data
    read(16,*)
    read(16,*)fxexp_data
    read(16,*)
    read(16,*)csoil_data
    read(16,*)
    read(16,*)salp_data
    read(16,*)
    read(16,*)refdk_data
    read(16,*)
    read(16,*)refkdt_data
    read(16,*)
    read(16,*)frzk_data
    read(16,*)
    read(16,*)zbot_data
    read(16,*)
    read(16,*)czil_data
    read(16,*)
    read(16,*)smlow_data
    read(16,*)
    read(16,*)smhigh_data
    read(16,*)
    read(16,*)lvcoef_data
    close(16)
 endif

 DM_BCAST_INTEGER(num_slope)
 DM_BCAST_INTEGER(slpcats)
 DM_BCAST_REALS(slope_data)
 DM_BCAST_REAL(sbeta_data)
 DM_BCAST_REAL(fxexp_data)
 DM_BCAST_REAL(csoil_data)
 DM_BCAST_REAL(salp_data)
 DM_BCAST_REAL(refdk_data)
 DM_BCAST_REAL(refkdt_data)
 DM_BCAST_REAL(frzk_data)
 DM_BCAST_REAL(zbot_data)
 DM_BCAST_REAL(czil_data)
 DM_BCAST_REAL(smlow_data)
 DM_BCAST_REAL(smhigh_data)
 DM_BCAST_REAL(lvcoef_data)

!call iamas_log_write('')
!write(mess,*) 'input general parameters'
!call physics_message(mess)
!call iamas_log_write('NUM_SLOPE=$i',intArgs=(/num_slope/))
!do lc = 1, slpcats
!   call iamas_log_write('$i $r', intArgs=(/lc/),realArgs=(/slope_data(lc)/))
!enddo
!call iamas_log_write('       end subroutine soil_veg_gen_parm:')
 call iamas_log_write('    end read GENPARM.TBL')

 end subroutine soil_veg_gen_parm

!=================================================================================================================
 end module iamas_atmphys_lsm_noahinit
!=================================================================================================================
