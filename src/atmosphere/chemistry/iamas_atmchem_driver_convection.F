!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute convection for chem.
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 03/2021
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

module iamas_atmchem_driver_convection

   use iamas_kind_types
   use iamas_pool_routines
   use iamas_atmchem_vars
   use iamas_atmchem_utilities
   use iamas_atmphys_vars
   use iamas_atmchem_constants, only: xlv, g, r_v
   use iamas_constants, only : cp
   !chemistry schemes:
   use module_ctrans_kfeta, only: kfeta_cup_transport
   use module_ctrans_general, only: general_cup_transport
   use iamas_log
   use iamas_derived_types

   implicit none

   private

   public :: driver_convection

contains

   subroutine driver_convection(state, diag_physics, time_lev, configs, diag_chemistry, tend_chemistry, its, ite)

      type(iamas_pool_type), intent(in) :: configs
      type(iamas_pool_type), intent(inout) :: state
      type(iamas_pool_type), intent(inout) :: tend_chemistry
      type(iamas_pool_type), intent(inout) :: diag_physics
      type(iamas_pool_type), intent(inout) :: diag_chemistry
      integer, intent(in) :: its, ite, time_lev
      
      real(kind = RKIND), dimension(:, :), pointer :: umf, uer, udr
      real(kind = RKIND), dimension(:, :), pointer :: dmf, der, ddr, wu
      real(kind = RKIND), dimension(:), pointer :: cutop, cubot, nca
      integer, dimension(:), pointer :: k22_shalc, kbcon_shalc, ktop_shalc
      real(kind = RKIND), dimension(:), pointer :: xmb_shalc
      integer, dimension(:), pointer :: ldb, let, ksl, lfs, klc
      integer, dimension(:), pointer :: ishall, check_kfeta

      real(kind = RKIND), dimension(:, :, :), pointer :: tend_ctrans
      real(kind = RKIND), dimension(:, :, :), pointer :: tend_wetconv
      real(kind = RKIND), dimension(:, :), pointer :: aq_aer_ratio
      integer, dimension(:), pointer :: ktop_deepc, kbcon_deepc

      !local vars
      logical, pointer :: scav_on
      character(len = StrKIND), pointer :: convection_scheme, chem_opt, ctrans_opt
      logical,pointer :: config_do_chem_budget

      call iamas_pool_get_config(configs, 'config_chem_opt', chem_opt)
      call iamas_pool_get_config(configs, 'config_cu_scav_on', scav_on)
      call iamas_pool_get_config(configs, 'config_convection_scheme', convection_scheme)
      call iamas_pool_get_config(configs, 'config_ctrans_opt', ctrans_opt)
      call iamas_pool_get_config(configs,'config_do_chem_budget' ,config_do_chem_budget )

      call iamas_pool_get_array(diag_physics, 'umf'        , umf         )
      call iamas_pool_get_array(diag_physics, 'uer'        , uer         )
      call iamas_pool_get_array(diag_physics, 'udr'        , udr         )
      call iamas_pool_get_array(diag_physics, 'dmf'        , dmf         )
      call iamas_pool_get_array(diag_physics, 'der'        , der         )
      call iamas_pool_get_array(diag_physics, 'ddr'        , ddr         )
      call iamas_pool_get_array(diag_physics, 'wu'         , wu          )
      call iamas_pool_get_array(diag_physics, 'ldb'        , ldb         )
      call iamas_pool_get_array(diag_physics, 'let'        , let         )
      call iamas_pool_get_array(diag_physics, 'ksl'        , ksl         )
      call iamas_pool_get_array(diag_physics, 'lfs'        , lfs         )
      call iamas_pool_get_array(diag_physics, 'klc'        , klc         )
      call iamas_pool_get_array(diag_physics, 'cutop'      , cutop       )
      call iamas_pool_get_array(diag_physics, 'cubot'      , cubot       )
      call iamas_pool_get_array(diag_physics, 'nca'        , nca         )
      call iamas_pool_get_array(diag_physics, 'ishall'     , ishall      )
      call iamas_pool_get_array(diag_physics, 'check_kfeta', check_kfeta )

      call iamas_pool_get_array(tend_chemistry, 'tend_ctrans', tend_ctrans)
      call iamas_pool_get_array(tend_chemistry, 'tend_wetconv', tend_wetconv)
      call iamas_pool_get_array(diag_chemistry, 'aq_aer_ratio', aq_aer_ratio)
      call iamas_pool_get_array(diag_physics, 'kbcon_deepc', kbcon_deepc)
      call iamas_pool_get_array(diag_physics, 'ktop_deepc', ktop_deepc)

      call iamas_pool_get_array(diag_physics, 'k22_shalc',   k22_shalc  )
      call iamas_pool_get_array(diag_physics, 'kbcon_shalc', kbcon_shalc)
      call iamas_pool_get_array(diag_physics, 'ktop_shalc',  ktop_shalc )
      call iamas_pool_get_array(diag_physics, 'xmb_shalc',   xmb_shalc  )

      kbcon_deepc(its:ite) = 0
      ktop_deepc(its:ite) = 0
      aq_aer_ratio(kts:kte,its:ite) = 0
      tend_wetconv(1:num_chem,kts:kte,its:ite) = 0.0

      convection_select: select case(trim(ctrans_opt))
      case ("kfeta")
         call chemistry_debug_message('convection_DRIVER CALLING kfeta')
         call kfeta_cup_transport(                                     &
              dt=dt_dyn, dx=dx(its:ite), rho_phy=rho_phy(kts:kte,its:ite), alt=alt(kts:kte,its:ite)     , &
              chem=chem(1:num_chem,kts:kte,its:ite), u=u_phy(kts:kte,its:ite), v=v_phy(kts:kte,its:ite), t_phy=t_phy(kts:kte,its:ite)       , &
              dz8w=dz(kts:kte,its:ite), p_phy=pres_hyd(kts:kte,its:ite), z=zgrid(kts:kte+1,its:ite)                     , &
              umf3d=umf(kts:kte,its:ite), uer3d=uer(kts:kte,its:ite), udr3d=udr(kts:kte,its:ite), dmf3d=dmf(kts:kte,its:ite)     , &
              der3d=der(kts:kte,its:ite), ddr3d=ddr(kts:kte,its:ite),  wu3d= wu(kts:kte,its:ite), ldb2d=ldb(its:ite)     , &
              let2d=let(its:ite), ksl2d=ksl(its:ite), lfs2d=lfs(its:ite), klc2d=klc(its:ite)     , &
              cutop=int(cutop(its:ite)), cubot=int(cubot(its:ite)),nca2d=nca(its:ite)     , &
              ishall2d=ishall(its:ite), check2d=check_kfeta(its:ite)               , &
              num_chem=num_chem, chem_opt=chem_opt, scav_opt=scav_on , &
              its=its, ite=ite, kts=kts, kte=kte )
      case ("general")
         !jwfeng: *_shallow vars only present in cu_grell_freitas scheme
         !        if you want to use general ctrans module, you must set
         !        ishallow to 0 to avoid using *_shallow vars.
         if (trim(convection_scheme) .ne. "cu_grell_freitas") then
            call general_cup_transport(                                     &
                 dt=dt_dyn, rho_phy=rho_phy(kts:kte,its:ite), raincv=raincv(its:ite)            , &
                 chem=chem(1:num_chem,kts:kte,its:ite), wetdep_3d = tend_wetconv(1:num_chem,kts:kte,its:ite)          , &
                 u=u_phy(kts:kte,its:ite), v=v_phy(kts:kte,its:ite), t_phy=t_phy(kts:kte,its:ite)                      , &
                 qv=qv(kts:kte,its:ite), dz8w=dz(kts:kte,its:ite), p_phy=pres_hyd(kts:kte,its:ite), xlv=xlv            , &
                 cp=cp, g=g, r_v=r_v, z=zgrid(kts:kte+1,its:ite)                               , &
                 tracert_3d = tend_ctrans(1:num_chem,kts:kte,its:ite), aq_aer_ratio = aq_aer_ratio(kts:kte,its:ite), &
                 ktop_2d = ktop_deepc(its:ite), kbcon_2d = kbcon_deepc(its:ite)         , &
                 ishallow=0, num_gas=num_gas, num_chem=num_chem           , &
                 chem_opt=chem_opt, cu_scav_on=scav_on, conv_tr_aqchem=0  , &
                 its=its, ite=ite, kts=kts, kte=kte)
         else
            ! jwfeng: GF can calculate shallow convection through setting ishallow=1.
            !         But we don't do shallow for now to compared to KF
            call general_cup_transport(                                     &
                 dt=dt_dyn, rho_phy=rho_phy(kts:kte,its:ite), raincv=raincv(its:ite)            , &
                 chem=chem(1:num_chem,kts:kte,its:ite), wetdep_3d = tend_wetconv(1:num_chem,kts:kte,its:ite)          , &
                 u=u_phy(kts:kte,its:ite), v=v_phy(kts:kte,its:ite), t_phy=t_phy(kts:kte,its:ite)                      , &
                 qv=qv(kts:kte,its:ite), dz8w=dz(kts:kte,its:ite), p_phy=pres_hyd(kts:kte,its:ite), xlv=xlv            , &
                 cp=cp, g=g, r_v=r_v, z=zgrid(kts:kte+1,its:ite)                              , &
                 tracert_3d = tend_ctrans(1:num_chem,kts:kte,its:ite), aq_aer_ratio = aq_aer_ratio(kts:kte,its:ite), &
                 ktop_2d = ktop_deepc(its:ite), kbcon_2d = kbcon_deepc(its:ite)         , &
                 k22_shallow=k22_shalc(its:ite), kbcon_shallow=kbcon_shalc(its:ite)     , &
                 ktop_shallow=ktop_shalc(its:ite), xmb_shallow=xmb_shalc(its:ite)       , &
                 ishallow=0, num_gas=num_gas, num_chem=num_chem           , &
                 chem_opt=chem_opt, cu_scav_on=scav_on, conv_tr_aqchem=0  , &
                 its=its, ite=ite, kts=kts, kte=kte)
         end if
      case default
         call chemistry_error_fatal('Wrong ctrans_opt'//trim(ctrans_opt)//']' )
      end select convection_select

      call chemistry_debug_message('convection scheme success')
      if (config_do_chem_budget) then
         call iamas_pool_get_array(tend_chemistry, 'tend_ctrans', tend_ctrans)
         to_iamas: select case(trim(ctrans_opt))
         case ("kfeta")
            tend_ctrans(1:num_chem,kts:kte,its:ite) = (chem(1:num_chem,kts:kte,its:ite) - chem_old(1:num_chem,kts:kte,its:ite)) / dt_dyn
         end select to_iamas
      end if

   end subroutine driver_convection

end module iamas_atmchem_driver_convection
