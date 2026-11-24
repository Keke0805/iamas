!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute dry deposition tendencies for chem(:,:,:,:). 
! * Author: Jiawang Feng,ChenJin and Chun Zhao
! Last update: chenjin 11/2023
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
! Terms of Use:
!  1) Users are requested to consult the primary author prior to
!     modifying this module or incorporating it or its submodules in
!     another code. This is meant to ensure that the any linkages and/or
!     assumptions will not adversely affect the operation of this module.
!  2) The source code in this module is intended for research and
!     educational purposes. Users are requested to contact the primary
!     author regarding the use of the code for any commercial application.
!  3) Users preparing publications resulting from the usage of this code
!     are requested to cite one or more of the references below
!     (depending on the application) for proper acknowledgement.
!
! References: 
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

!===========================================================================
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + Dry deposition calculation. 
!  + Tendency calculation for budget analysis 
!                                         ---jwfeng   12/2019
!=============================================================================

MODULE iamas_atmchem_driver_drydep
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_configs
 use iamas_atmchem_utilities
 use iamas_atmphys_vars
!chemistry schemes:
 use module_dep_simple
 use module_gocart_drydep
 use module_aer_drydep
 use module_vertmx_wrf
 use module_vertmx_kokkos, only: vertmx_kokkos_batch, vertmx_kokkos_mode
 use module_mosaic_settling
 use module_mixactivate_wrappers
 use iamas_log
 use iamas_timer
!use module_mosaic_drydep

 implicit none
 private
 public:: driver_drydep

!iAMAS driver for chemistry drydep.

! subroutines in iamas_atmchem_driver_drydep:
! ---------------------------------------
! driver_drydep    : main driver (called from subroutine chemistry_driver).


 contains

!=================================================================================================================
 subroutine driver_drydep(itimestep,configs,mesh,state,sfc_input,diag,diag_physics,tend,tend_chemistry,its,ite)
!=================================================================================================================
!input arguments:
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep
!inout arguments:
 type(iamas_pool_type),intent(inout):: diag
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend
 type(iamas_pool_type),intent(inout):: tend_chemistry
 type(iamas_pool_type),intent(inout):: state
 
!local variables:
 integer :: i,k,j,h,nv,ll,l
 integer :: vertmx_mode
 real(kind=RKIND), allocatable :: vertmx_candidate(:,:,:)
 real(kind=RKIND) :: vertmx_error, vertmx_scale
 character(len=256) :: vertmx_message
 integer :: itf,ktf,jtf
 integer :: k_a, k_c, kmax, m_mam
 integer :: idrydep_onoff,aer_mech_id
 integer :: iphase, itype
 integer,pointer :: gas_drydep_option
 real,parameter  :: m2cm = 100.
 real,dimension(1:num_chem,1:config_flags%kdepvel,its:ite) :: dep_vel 
 real,dimension(kts:kte) :: dryrho_1d
 real,dimension(its:ite) :: frac_removed
 real,dimension(kts:kte,its:ite) :: dryrho_phy
 real(kind=RKIND),dimension(:,:), pointer :: ddvel
 real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
 real(kind=RKIND),dimension(:), pointer :: lonCell
 real(kind=RKIND),dimension(:)  ,pointer ::  latCell
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn1
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn2
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn3
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn4
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn5
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn6
 real(kind=RKIND),dimension(:,:)  ,pointer :: nsource
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nihf3d
 real(kind=RKIND),dimension(:,:)  ,pointer :: niimm3d
 real(kind=RKIND),dimension(:,:)  ,pointer :: nidep3d
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst1
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst2
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst3
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst4
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon1
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon2
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon3
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon4
 real(kind=RKIND),dimension(:,:)  ,pointer ::  wsub3d
 real(kind=RKIND),dimension(:,:)  ,pointer :: tend_qndrop_activation
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_mixdry
 real(kind=RKIND),dimension(:,:,:),pointer :: v_grav_3d

 !zhxia add this for aer_drydep. zhxia@2026
 character(len=StrKIND), pointer :: mminlu

 !zhxia add this for determining which scheme is used for vertical mixing calculation.
 character(len=StrKIND),pointer :: chem_opt

!turbulent transport
 real :: pblst(kts:kte),ekmfull(kts:kte+1),zzfull(kts:kte+1),zz(kts:kte)
 logical,pointer :: config_do_chem_budget
 real, dimension(its:ite) ::                      &
                              aer_res,            &
                              aer_res_def,        &
                              aer_res_zcen
 real, dimension(kts:kte,its:ite ) ::             &
                                    h2oai,        &
                                    h2oaj
 character(len=StrKIND):: debugString
!-----------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------

!jwfeng: Whether to calculate budget (tendencies)
 call iamas_timer_start('dd_get_fields')
 call iamas_pool_get_config(configs,'config_gas_drydep_opt' ,gas_drydep_option )
 call iamas_pool_get_config(configs,'config_do_chem_budget' ,config_do_chem_budget )
 call iamas_pool_get_config(configs,'config_chem_opt'       ,chem_opt      )
 call iamas_pool_get_array(mesh          ,'latCell'         ,latCell       )
 call iamas_pool_get_array(diag          ,'relhum'          ,relhum        )
 call iamas_pool_get_array(mesh          ,'lonCell'         ,lonCell       )
 call iamas_pool_get_array(diag          ,'nihf3d'         ,nihf3d        )
 call iamas_pool_get_array(diag          ,'niimm3d'        ,niimm3d       )
 call iamas_pool_get_array(diag          ,'nidep3d'        ,nidep3d       )
 call iamas_pool_get_array(diag          ,'wsub3d'         ,wsub3d        )
 call iamas_pool_get_array(diag          ,'ccn1'           ,ccn1          )
 call iamas_pool_get_array(diag          ,'ccn2'           ,ccn2          )
 call iamas_pool_get_array(diag          ,'ccn3'           ,ccn3          )
 call iamas_pool_get_array(diag          ,'ccn4'           ,ccn4          )
 call iamas_pool_get_array(diag          ,'ccn5'           ,ccn5          )
 call iamas_pool_get_array(diag          ,'ccn6'           ,ccn6          )
 call iamas_pool_get_array(diag          ,'nsource'        ,nsource       )
 call iamas_pool_get_array(diag          ,'rndst1'         ,rndst1        )
 call iamas_pool_get_array(diag          ,'rndst2'         ,rndst2        )
 call iamas_pool_get_array(diag          ,'rndst3'         ,rndst3        )
 call iamas_pool_get_array(diag          ,'rndst4'         ,rndst4        )
 call iamas_pool_get_array(diag          ,'nacon1'         ,nacon1        )
 call iamas_pool_get_array(diag          ,'nacon2'         ,nacon2        )
 call iamas_pool_get_array(diag          ,'nacon3'         ,nacon3        )
 call iamas_pool_get_array(diag          ,'nacon4'         ,nacon4        )
 call iamas_pool_get_array(tend_chemistry,'ddvel'          ,ddvel         )
 call iamas_pool_get_array(tend_chemistry,'v_grav_3d'      ,v_grav_3d     )
 call iamas_pool_get_array(tend          ,'tend_qndrop_activation'      ,tend_qndrop_activation     )
 call iamas_pool_get_array(sfc_input     , 'mminlu'        ,mminlu) !ensure that mminlu contains the name of the land use dataset. zhxia@2026
 call iamas_timer_stop('dd_get_fields')
 call iamas_timer_start('dd_init_arrays')
 if (config_do_chem_budget) then
      chem_old(1:num_chem,kts:kte,its:ite) = chem(1:num_chem,kts:kte,its:ite)
 end if
  !initialization:
  v_grav_3d(1:num_chem,kts:kte,its:ite)  = 0._RKIND
  ddvel(1:num_chem,its:ite)      = 0._RKIND
!-----------------------------------------------------------------------------------------------------------------
 aer_res_def     = 0
 aer_res_zcen    = 0
 aer_res         = 0

 jtf=min0(jte,jde-1)
 ktf=min0(kte,kde-1)
 itf=min0(ite,ide-1)
 call iamas_timer_stop('dd_init_arrays')
!-----------------------------------------------------------------------------------------------------------------
 
!
! compute dry deposition velocities = ddvel
!
! 28-jun-2005 rce - initialize ddvel=0; call aerosol drydep routine
!           only when gas_drydep_opt == WESELY
!       the wesely_driver routine computes aer_res, and currently
!	    you cannot compute aerosol drydep without it !!
! 08-jul-2005 rce - pass idrydep_onoff to mixactivate routines
!
 idrydep_onoff = 0
 call chemistry_debug_message('drydep_driver success ')

 gas_drydep_select: SELECT CASE(gas_drydep_option)
     CASE (1)
!1 means use wesely module
!
! gas_drydep_opt == WESELY means 
!     wesely for gases 
!     other (appropriate) routine for aerosols
!
!      CALL wrf_debug(15,'DOING DRY DEP VELOCITIES WITH WESELY METHOD')
         
 call iamas_timer_start('dd_gas_velocity')
         call wesely_driver(id=itimestep,ktau=itimestep,dtstep=dt_dyn,                                                                          &
              config_flags=config_flags,gmt=gmt,                                                                                                &
              julday=int(curr_julday),t_phy=t_phy(kts:kte,its:ite),moist=moist(1:num_moist,kts:kte,its:ite),                                    &
              p8w=pres2_hyd(kts:kte+1,its:ite),t8w=t2(kts:kte+1,its:ite),raincv=raincv(its:ite),                                                &
              p_phy=pres_hyd(kts:kte,its:ite),chem=chem(1:num_chem,kts:kte,its:ite),rho_phy=rho_phy(kts:kte,its:ite),dz8w=dz(kts:kte,its:ite),  &
              ddvel=ddvel(1:num_chem,its:ite),aer_res_def=aer_res_def(its:ite),                                                                 &
              aer_res_zcen=aer_res_zcen(its:ite),ivgtyp=ivgtyp(its:ite),tsk=skintemp(its:ite),                                                  &
              gsw=gsw(its:ite),vegfra=vegfra(its:ite),pbl=hpbl(its:ite),rmol=rmol(its:ite),ust=ust(its:ite),                                    &
              znt=znt(its:ite),xlat=xlat(its:ite),xlong=xlong(its:ite),z=zgrid(kts:kte+1,its:ite),z_at_w=zgrid(kts:kte+1,its:ite),              &
              snowh=snowh(its:ite),numgas=num_gas,                                                                                              &
              its=its,ite=ite,kts=kts,kte=kte,                                                                                                  &
              latCell=latCell(its:ite),lonCell=lonCell(its:ite)                  )
 call iamas_timer_stop('dd_gas_velocity')
              call chemistry_debug_message('wesely_driver success')
       
     CASE DEFAULT 
                                                     
 END SELECT gas_drydep_select 

!       IF(config_flags%chem_opt /= CHEM_TRACER    .and.                   &
!          config_flags%chem_opt /= CHEM_TRACE2    .and.                   &
!          config_flags%chem_opt /= CO2_TRACER     .and.                   & 
!          config_flags%chem_opt /= GHG_TRACER     .and.                   &
!          config_flags%chem_opt /= CHEM_VASH      .and.                   &
!          config_flags%chem_opt /= CHEM_VOLC_4BIN .and.                   &
!          config_flags%chem_opt /= DUST           .and.                   &
!          config_flags%chem_opt /= GOCART_SIMPLE  .and.                   &
!          config_flags%chem_opt /= GOCARTRACM_KPP )THEN     
!         call wesely_driver(id,ktau,dtstep,                              &
!              config_flags,current_month,                                &
!              gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
!              p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
!              ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
!              snowh,numgas,                                              &
!              ids,ide, jds,jde, kds,kde,                                 &
!              ims,ime, jms,jme, kms,kme,                                 &
!              its,ite, jts,jte, kts,kte                                  )
!-----------------------------------------------------------------
! this does aerosol species (dust,seas, bc,oc) for gocart only
!-----------------------------------------------------------------
!          IF(config_flags%chem_opt == MOZCART_KPP ) then
!            call gocart_drydep_driver( dtstep,                            &
!                  config_flags, numgas,                                   &
!                  t_phy, moist, p8w, t8w, rmol,aer_res_def,               &
!                  p_phy, chem, rho_phy, dz8w, ddvel, xland, hfx,          &
!                  ivgtyp, tsk, vegfra, pbl, ust, znt, xlat, xlong,        &
!                  ids,ide, jds,jde, kds,kde,                              &
!                  ims,ime, jms,jme, kms,kme,                              &
!                  its,ite, jts,jte, kts,kte )
!          ENDIF
!       ELSEIF ( config_flags%chem_opt == GOCART_SIMPLE ) then
!         call wesely_driver(id,ktau,dtstep,                              &
!              config_flags,current_month,                                &
!              gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
!              p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    & 
!              ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
!              snowh, numgas,                                             &
!              ids,ide, jds,jde, kds,kde,                                 &
!              ims,ime, jms,jme, kms,kme,                                 &
!              its,ite, jts,jte, kts,kte                                  )
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
!         call gocart_drydep_driver(dtstep,                                &
!               config_flags,numgas,                                       &
!               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
!               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
!               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
!               ids,ide, jds,jde, kds,kde,                                 &
!               ims,ime, jms,jme, kms,kme,                                 &
!               its,ite, jts,jte, kts,kte                                  )
!	   ELSEIF ( config_flags%chem_opt == DUST ) then
!
! this does aerosol species (dust) for gocart only
!
!         call gocart_drydep_driver(dtstep,                                &
!              config_flags,numgas,                                        &
!              t_phy,moist,p8w,t8w,rmol,aer_res_def,                       &
!              p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                    &
!              ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                   &
!              ids,ide, jds,jde, kds,kde,                                  &
!              ims,ime, jms,jme, kms,kme,                                  &
!              its,ite, jts,jte, kts,kte                                   )
!
!
!       ELSEIF ( config_flags%chem_opt == GOCARTRACM_KPP) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
!         call wesely_driver(id,ktau,dtstep,                              &
!              config_flags,current_month,                                &
!              gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
!              p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
!              ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
!              snowh, numgas,                                             &
!              ids,ide, jds,jde, kds,kde,                                 &
!              ims,ime, jms,jme, kms,kme,                                 &
!              its,ite, jts,jte, kts,kte                                  )
!         call gocart_drydep_driver(dtstep,                                &
!               config_flags,numgas,                                       &
!               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
!               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
!               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
!               ids,ide, jds,jde, kds,kde,                                 &
!               ims,ime, jms,jme, kms,kme,                                 &
!               its,ite, jts,jte, kts,kte                                  )

!       ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
!          ddvel(its:ite,jte:jte,1:num_chem) = 0.
!       END IF
!       ddvel temporarily close,maybe this is wrong ;chenjin
!      if (config_flags%aer_aerodynres_opt == 2) then
!         ! use aerodynamic resistance from center of layer kts to surface
!         aer_res(its:ite,jts:jte) = aer_res_zcen(its:ite,jts:jte)
!      else
!         ! this is the default -- use aerodynamic resistance from 
!         ! "default reference height" (currently 2 m) to surface
 call iamas_timer_start('dd_aer_resistance')
       aer_res(its:ite) = aer_res_def(its:ite)
       write(debugString, *) "--- Min/Max of aer_res = ", maxval(aer_res(its:ite)),"/",minval(aer_res(its:ite))
       call chemistry_debug_message(debugString)
!      end if

       idrydep_onoff = 1
       aer_mech_id = 3
 call iamas_timer_stop('dd_aer_resistance')

! !  config_flags%aer_drydep_opt <= 0   -- aerosol depositon velocities are set to zero
! !  config_flags%aer_drydep_opt >= 100 -- aerosol depvels are calculated in subr aer_depvel_driver
! !                                        (see module_aer_drydep.F for details)
! !  config_flags%aer_drydep_opt == 1   -- SORGAM chem packages use subr vdvg_2 of module_aerosols_sorgam.F
! !                                        MOSAIC chem packages subr aerosol_depvel_2 of module_mosaic_drydep.F
! !  config_flags%aer_drydep_opt == 11  -- SORGAM chem packages use subr vdvg of module_aerosols_sorgam.F

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!next is the calculation of the aerosol ddvel
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if ((config_flags%aer_drydep_opt <= 0) .or. (aer_mech_id <= 0)) then
          call chemistry_debug_message('AEROSOL DRY DEP VELOCITIES  = 0.0')
       else if (config_flags%aer_drydep_opt <= 99) then
                                                    
! !    END SELECT adrydep_select
       else ! (config_flags%aer_drydep_opt > 99)
          call chemistry_debug_message('DOING DRY DEP VELOCITIES THRU AER_DRYDEP_DRIVER')
 call iamas_timer_start('dd_aer_velocity')
          call aer_drydep_driver(                                                                                                                           &
                  id=itimestep, ktau=itimestep, dtstep=dt_dyn,                                                                                              &
                  config_flags=config_flags, aer_mech_id=aer_mech_id,                                                                                       & 
                  gmt=gmt, julday=int(curr_julday), mminlu = trim(mminlu),                                                                                  & 
                  t_phy=t_phy(kts:kte,its:ite),rho_phy=rho_phy(kts:kte,its:ite), p_phy=pres_hyd(kts:kte,its:ite),                                           &
                  alt=alt(kts:kte,its:ite),p8w=pres2_hyd(kts:kte+1,its:ite),t8w=t2(kts:kte+1,its:ite),dz8w=dz(kts:kte,its:ite),z=zgrid(kts:kte+1,its:ite),  &
                  z_at_w=zgrid(kts:kte+1,its:ite), ust=ust(its:ite),aer_res=aer_res(its:ite), ivgtyp=ivgtyp(its:ite),                                       &
                  vegfra=vegfra(its:ite), pbl=hpbl(its:ite),rmol=rmol(its:ite), znt=znt(its:ite),                                                           &
                  moist=moist(1:num_moist,kts:kte,its:ite),chem=chem(1:num_chem,kts:kte,its:ite), ddvel=ddvel(1:num_chem,its:ite),                          &
                  h2oai=h2oai(kts:kte,its:ite),h2oaj=h2oaj(kts:kte,its:ite), numgas=num_gas,                                                                &
                  its=its,ite=ite,kts=kts,kte=kte)
 call iamas_timer_stop('dd_aer_velocity')

      end if
 call iamas_timer_start('dd_velocity_prepare')
      if (config_flags%aer_drydep_opt > 0) then
         if ((aer_mech_id > 0) .and. (aer_mech_id <= 4)) then
             ! limit aerosol ddvels to <= 0.5 m/s
             ! drydep routines occasionally produce unrealistically-large particle
             !    diameter leading to unrealistically-large sedimentation velocity
          ddvel(num_gas+1:num_chem,its:ite) =min(0.50,ddvel(num_gas+1:num_chem,its:ite) )
         end if
      end if

      !write(*,*) ddvel(num_gas+1,its:ite)
! !      ddvel(its:ite,jts:jte,p_o3) = 0.1*ddvel(its:ite,jts:jte,p_o3)
! !      ddvel(its:ite,jts:jte,1:num_chem) = 0
! !       !For MAM aersols dry deposition
! !       if ((aer_mech_id == 4)) then
! !          do m_mam = 1, num_chem
! !             !Look for species with _a or _c and adjust these species for dry deposition
! !             k_a   = index(trim(adjustl(chem_dname_table(1,m_mam))),'_a') !To include *_a* species
! !             k_c   = index(trim(adjustl(chem_dname_table(1,m_mam))),'_c') !To include *_c* species
! !             kmax = max(k_a, k_c)
! !             if(kmax > 0 ) then !kmax>0 means a species has _a or _c
! !                frac_removed(its:ite,jts:jte) = max( 0.0, min( 1.0, ddvel(its:ite,jts:jte,m_mam)*dtstep/dz8w(its:ite,1,jts:jte) ) )
! !                chem(its:ite,1,jts:jte,m_mam) = chem(its:ite,1,jts:jte,m_mam)*(1.0 - frac_removed(its:ite,jts:jte)) !update lowest layer
! !             endif
! !          enddo
! !       endif
                             

! ! Add dep_vel here for all aerosol and gas species in the dvel array 
! !ll = max( 1, min( config_flags%ndepvel, num_vert_mix ) )
      ll = num_chem    ! not confirmed now
      dep_vel(:,:,:) = 0.
      do i=its,ite
      do k=1,config_flags%kdepvel
      do l=1,ll
         dep_vel(l,k,i) = ddvel(l,i)
      enddo
      enddo
      enddo
 call iamas_timer_stop('dd_velocity_prepare')

! !   This will be called later from subgrd_transport_driver.F !!!!!!!!
! !
! !
! !     dep_vel_o3=0.
! !      if (num_vert_mix == 1) then
     call iamas_timer_start('vertmx')
     vertmx_mode = vertmx_kokkos_mode()
     if (vertmx_mode == 1) then
        call vertmx_kokkos_batch(dt_dyn, epsilc, index(trim(chem_opt), "aq") /= 0, &
             chem(1:num_chem,kts:kte,its:ite), alt(kts:kte,its:ite), &
             zgrid(kts:kte+1,its:ite), kzh(kts:kte,its:ite), &
             ddvel(1:num_chem,its:ite), is_aerosol(1:num_chem))
     else
        if (vertmx_mode == 2) then
           allocate(vertmx_candidate(num_chem,kte-kts+1,ite-its+1))
           vertmx_candidate = chem(1:num_chem,kts:kte,its:ite)
           call vertmx_kokkos_batch(dt_dyn, epsilc, index(trim(chem_opt), "aq") /= 0, &
                vertmx_candidate, alt(kts:kte,its:ite), zgrid(kts:kte+1,its:ite), &
                kzh(kts:kte,its:ite), ddvel(1:num_chem,its:ite), is_aerosol(1:num_chem))
        endif
     do i=its,ite
          pblst=0.
!
!====================================================================================
!Start with vertical turbulent mixing    -czhao
!====================================================================================
!
      do k=kts,kte+1
         zzfull(k)=zgrid(k,i)-zgrid(kts,i)
      enddo
      do k=kts,kte
         ekmfull(k)=max(1.e-6,kzh(k,i))
      enddo
      ekmfull(kts)=0.
      ekmfull(kte+1)=0.
      do k=kts,kte
         zz(k)=(zgrid(k+1,i)+zgrid(k,i))/2.0 -zgrid(kts,i)
      enddo
! !     dep_vel_o3(i,j)=ddvel(i,j,p_o3)

! !Jiawang Feng:nv start at 2 in wrf-chem.
      do nv=1,num_chem
!        if(is_CAMMGMP_used .and. .not.vertMixAero(nv))cycle !Balwinder.Singh@pnnl.gov: Do mix constituents which are already mixed by CAMMGMP microphysics
         do k=kts,kte
            pblst(k)=max(epsilc,chem(nv,k,i))
            dryrho_1d(k) = 1./alt(k,i)
         enddo

         ! zhxia add this for determining which scheme is used for vertical scheme calculation.
         if (index(trim(chem_opt), "aq") /= 0) then
            !Vertical turbulent mixing with dry deposition for gases only. zhxia@2026
            if(.not.is_aerosol(nv))then ! mix gases not aerosol
                call vertmx(dt_dyn,pblst,ekmfull,dryrho_1d, &
                            zzfull,zz,ddvel(nv,i),kts,kte)
            end if

         else
            !Vertical turbulent mixing with dry deposition for gases and aerosols. zhxia@2026
            call vertmx(dt_dyn,pblst,ekmfull,dryrho_1d, &
                        zzfull,zz,ddvel(nv,i),kts,kte)
         endif

         do k=kts,kte-1
         chem(nv,k,i)=max(epsilc,pblst(k))
         enddo
      enddo
     enddo !i loop
        if (vertmx_mode == 2) then
           call iamas_timer_start('vertmx_validate')
           if (.not. all(ieee_is_finite(vertmx_candidate)) .or. &
               .not. all(ieee_is_finite(chem(1:num_chem,kts:kte,its:ite)))) then
              error stop 'vertmx validation: nonfinite concentration'
           endif
           vertmx_error = maxval(abs(vertmx_candidate-chem(1:num_chem,kts:kte,its:ite)))
           vertmx_scale = maxval(abs(vertmx_candidate-chem(1:num_chem,kts:kte,its:ite)) / &
                max(abs(chem(1:num_chem,kts:kte,its:ite)),epsilc))
           write(vertmx_message,'(a,i0,a,i0,a,es12.4,a,es12.4)') &
                'VERTMX_VALIDATE step=',itimestep,' values=',size(vertmx_candidate), &
                ' max_abs=',vertmx_error,' max_rel=',vertmx_scale
           call iamas_log_write(trim(vertmx_message), masterOnly=.false.)
           if (vertmx_scale > 2.e-13_RKIND) error stop 'vertmx validation: numerical mismatch'
           chem(1:num_chem,kts:kte,its:ite) = vertmx_candidate
           deallocate(vertmx_candidate)
           call iamas_timer_stop('vertmx_validate')
        endif
     endif
     call iamas_timer_stop('vertmx')
! !      endif ! num_vert_mix = 0
! !
! !  vertical mixing and activation of aerosol
! !
 call iamas_timer_start('dd_density_prepare')
     where( alt(kts:kte,its:ite) /= 0. )  !get dry density to conserve mass in mixactivate, wig, 24-apr-2006
        dryrho_phy(kts:kte,its:ite) = 1./alt(kts:kte,its:ite)
     elsewhere
        dryrho_phy(kts:kte,its:ite) = 0.
     end where
 call iamas_timer_stop('dd_density_prepare')
     
 
     ! zhxia add this for determining which scheme is used for vertical scheme calculation. 
     if (index(trim(chem_opt), "aq") /= 0) then
        !----------------------------------------
        !Vertical turbulent mixing with dry deposition and aerosol activation for aerosols only. zhxia@2026
        !----------------------------------------
 call iamas_timer_start('dd_mixactivate')
        call mosaic_mixactivate (                                                                                                                                     &
             id=itimestep, ktau=itimestep, dtstep=dt_dyn, config_flags=config_flags,                                                                                  &
             idrydep_onoff=idrydep_onoff, p_phy=pres_hyd(kts:kte,its:ite),rho_phy=dryrho_phy(kts:kte,its:ite),                                                        &
             t_phy=t_phy(kts:kte,its:ite), w=w(kts:kte+1,its:ite),cldfra=cldfrac(kts:kte,its:ite), cldfra_old=cldfrac_old(kts:kte,its:ite),                           &
             ddvel=ddvel(1:num_chem,its:ite), z=zgrid(kts:kte+1,its:ite),dz8w=dz(kts:kte,its:ite), p_at_w=pres2_hyd(kts:kte+1,its:ite),t_at_w=t2(kts:kte+1,its:ite),  &
             exch_h=kzh(kts:kte,its:ite),qv=qv(kts:kte,its:ite),qc=qc(kts:kte,its:ite),                                                                               &
             qi=qi(kts:kte,its:ite),qs=qs(kts:kte,its:ite),f_ice_phy=f_ice_local(kts:kte,its:ite),f_rain_phy=f_rain_local(kts:kte,its:ite),qndrop3d=qndrop3d(kts:kte,its:ite),    &
             qnice3d=qnice3d(kts:kte,its:ite),tend_qndrop_activation=tend_qndrop_activation(kts:kte,its:ite),f_qv=f_qv,f_qc=f_qc,                                     &
             f_qi=f_qi,f_qs=f_qs, chem=chem(1:num_chem,kts:kte,its:ite), ccn1=ccn1(kts:kte,its:ite), ccn2=ccn2(kts:kte,its:ite), ccn3=ccn3(kts:kte,its:ite),          &
             ccn4=ccn4(kts:kte,its:ite), ccn5=ccn5(kts:kte,its:ite), ccn6=ccn6(kts:kte,its:ite), nsource=nsource(kts:kte,its:ite),                                    &
             nihf3d=nihf3d(kts:kte,its:ite),niimm3d=niimm3d(kts:kte,its:ite),nidep3d=nidep3d(kts:kte,its:ite),nimey3d=nimey3d(kts:kte,its:ite),                       &
             rndst1=rndst1(kts:kte,its:ite),rndst2=rndst2(kts:kte,its:ite),rndst3=rndst3(kts:kte,its:ite),rndst4=rndst4(kts:kte,its:ite),                             &
             nacon1=nacon1(kts:kte,its:ite),nacon2=nacon2(kts:kte,its:ite),nacon3=nacon3(kts:kte,its:ite),nacon4=nacon4(kts:kte,its:ite),                             &
             tke_pbl=tke_pbl(kts:kte,its:ite),relhum3d=relhum(kts:kte,its:ite),wsub3d=wsub3d(kts:kte,its:ite),                                                        &                              
             its=its,ite=ite,  kts=kts,kte=kte)
 call iamas_timer_stop('dd_mixactivate')
   
     endif  

 call iamas_timer_start('dd_positive_clip')
     chem(1:num_chem,kts:kte,its:ite) =max(epsilc,chem(1:num_chem,kts:kte,its:ite))   
 call iamas_timer_stop('dd_positive_clip')
     call chemistry_debug_message('call mosaic settling routine')
 call iamas_timer_start('dd_settling')
     call mosaic_settling_driver(dt=dt_dyn,config_flags=config_flags,t_phy=t_phy(kts:kte,its:ite),moist=moist(1:num_moist,kts:kte,its:ite),&
          chem=chem(1:num_chem,kts:kte,its:ite),rho_phy=rho_phy(kts:kte,its:ite),dz8w=dz(kts:kte,its:ite)                                  &
          ,p8w=pres2_hyd(kts:kte+1,its:ite),p_phy=pres_hyd(kts:kte,its:ite),dx=dx(its:ite),g=g,                                            &
          v_grav_3d=v_grav_3d(1:num_chem,kts:kte,its:ite),its=its,ite=ite, kts=kts,kte=kte) ! jwfeng add v_grav_3d
 call iamas_timer_stop('dd_settling')


     call chemistry_debug_message('drydep driver success')

 call iamas_timer_start('dd_budget_update')
     if (config_do_chem_budget) then
       call iamas_pool_get_array(tend_chemistry,'tend_mixdry',tend_mixdry)
       tend_mixdry(1:num_chem,kts:kte,its:ite) = (chem(1:num_chem,kts:kte,its:ite) - chem_old(1:num_chem,kts:kte,its:ite)) / dt_dyn
     end if
 call iamas_timer_stop('dd_budget_update')
 
 call iamas_timer_start('dd_scalar_writeback')
     scalars(p_qndrop,kts:kte,its:ite) = qndrop3d(kts:kte,its:ite)
     scalars(p_qniceaer,kts:kte,its:ite) = qnice3d(kts:kte,its:ite)
     scalars(p_qnicemix,kts:kte,its:ite) = nimey3d(kts:kte,its:ite)
 call iamas_timer_stop('dd_scalar_writeback')
     call chemistry_debug_message('drydep to iAMAS success')

 end subroutine driver_drydep

!=================================================================================================================
END MODULE iamas_atmchem_driver_drydep
!=================================================================================================================
