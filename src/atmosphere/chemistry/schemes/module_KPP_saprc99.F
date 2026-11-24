!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute KPP_saprc99 subroutine for chem(:,:,:,:). 
! * Author: Zihan Xia
! * Email: xiazh@mail.ustc.edu.cn
! Last update: zhxia 10/2025
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
 module module_KPP_saprc99 
    use iamas_kind_types
    use iamas_atmchem_vars
    use iamas_log, only : iamas_log_write
    use saprc99_Model

 contains
 subroutine KPP_saprc99_driver(itimestepc,mesh,state,                    &
                               chem,dtstepc,                             &
                               rho_phy,t_phy,p_phy,qv,                   &
                               ph, is_trop,                              &
                               its,ite,kts,kte)

 IMPLICIT NONE

 integer, intent(in) :: itimestepc
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(inout):: state
 integer, intent(in) :: its,ite,kts,kte
 real, dimension(1:num_chem, kts:kte, its:ite), intent(inout) :: chem
 real(kind=dp), intent(in) :: dtstepc
 real, dimension(kts:kte, its:ite), intent(in) :: rho_phy, t_phy, p_phy, qv
 real, dimension(1:num_phr, kts:kte, its:ite), intent(in) :: ph
 logical, dimension(kts:kte, its:ite), intent(in) :: is_trop


!local vars
 real(kind=dp) :: t, dval(nspec)
 integer, dimension(20) :: icntrl
 integer, dimension(20) :: istate
 real(kind=dp), dimension(20) :: rcntrl
 real(kind=dp), dimension(20) :: rstate
 real(kind=dp) :: conv, oconv
 integer :: i,k,n,time_lev
 real(kind=dp),dimension(:,:),pointer :: o1d,o3p,tbu_o,hocoo,bzno2_o,bz_o

 ! a few constants 
 real(kind=rkind), parameter  ::  navgdro = 6.022e23   ! molecules/mol
 real(kind=rkind), parameter  ::   &   !molecular weights
      mwh = 1.0079,  mwo = 15.9994!,  mwair = 28.97
 real(kind=rkind), parameter :: mwh2o = 2*mwh + mwo
 !dens2con air
 real(kind=rkind), parameter ::          dens2con_a = 1.e-3     & 
                              * (1./28.97)          &
                              * navgdro              
 !dens2con water
 real(kind=rkind), parameter ::          dens2con_w = 1.e-3     &
                               * (1./mwh2o)         &
                               * navgdro             

 ! initialization
 if (.not. allocated(jv)) then
    allocate( jv(num_phr) )
 end if
 jv = 0.0_dp 

 icntrl  = 0
 rcntrl  = 0.d0

 do i=1,nvar
    rtol(i) = 1.0d-3
    atol(i) = 1.0d0
 end do

 tstart = 0.0d0
 t = tstart
 tend   = real(dtstepc, kind=dp)
 time = t !current time, used to update sun. should be delete after we switch to j-values @zhxia

 icntrl(1) = 1
 icntrl(3) = 2
 !icntrl(15) = 5 !call update_sun and update_rconst from w/in the int. it can be adjust to change simulation speed. @zhxia
 icntrl(15) = 1 !Call Update_RCONST from within the integrator
 rcntrl(3) = 0.01_dp * real(dtstepc, kind=dp)

 if (itimestepc==1) then
     c = 0.0_dp
 endif
 
 time_lev = 1
 call iamas_pool_get_array(state         ,'o1d'         , o1d    , time_lev )
 call iamas_pool_get_array(state         ,'o3p'         , o3p    , time_lev )
 call iamas_pool_get_array(state         ,'tbu_o'       , tbu_o  , time_lev )
 call iamas_pool_get_array(state         ,'hocoo'       , hocoo  , time_lev )
 call iamas_pool_get_array(state         ,'bzno2_o'     , bzno2_o, time_lev )
 call iamas_pool_get_array(state         ,'bz_o'        , bz_o   , time_lev )

 do i=its, ite
    do k=kts, kte
       if (is_trop(k,i)) then
           conv=1.e-6_dp*dens2con_a*rho_phy(k,i) !ppm->molecules/cm³
           oconv = 1.e0_dp/conv
           temp = real(t_phy(k,i), kind=dp)
           cfactor = conv !impact fall function

           ! photolysis rate
           jv(p_ph_o2)           = real(max( ph(p_ph_o2          ,k,i),0.)/60., kind=dp)
           jv(p_ph_o31d)         = real(max( ph(p_ph_o31d        ,k,i),0.)/60., kind=dp)
           jv(p_ph_o33p)         = real(max( ph(p_ph_o33p        ,k,i),0.)/60., kind=dp)
           jv(p_ph_no2)          = real(max( ph(p_ph_no2         ,k,i),0.)/60., kind=dp)
           jv(p_ph_no3o2)        = real(max( ph(p_ph_no3o2       ,k,i),0.)/60., kind=dp)
           jv(p_ph_no3o)         = real(max( ph(p_ph_no3o        ,k,i),0.)/60., kind=dp)
           jv(p_ph_hno2)         = real(max( ph(p_ph_hno2        ,k,i),0.)/60., kind=dp)
           jv(p_ph_hno3)         = real(max( ph(p_ph_hno3        ,k,i),0.)/60., kind=dp)
           jv(p_ph_hno4)         = real(max( ph(p_ph_hno4        ,k,i),0.)/60., kind=dp)
           jv(p_ph_h2o2)         = real(max( ph(p_ph_h2o2        ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch2or)        = real(max( ph(p_ph_ch2or       ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch2om)        = real(max( ph(p_ph_ch2om       ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3cho)       = real(max( ph(p_ph_ch3cho      ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3coch3)     = real(max( ph(p_ph_ch3coch3    ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3coc2h5)    = real(max( ph(p_ph_ch3coc2h5   ,k,i),0.)/60., kind=dp)
           jv(p_ph_hcocho)       = real(max( ph(p_ph_hcocho      ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3cocho)     = real(max( ph(p_ph_ch3cocho    ,k,i),0.)/60., kind=dp)
           jv(p_ph_hcochest)     = real(max( ph(p_ph_hcochest    ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3o2h)       = real(max( ph(p_ph_ch3o2h      ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3coo2h)     = real(max( ph(p_ph_ch3coo2h    ,k,i),0.)/60., kind=dp)
           jv(p_ph_ch3ono2)      = real(max( ph(p_ph_ch3ono2     ,k,i),0.)/60., kind=dp)
           jv(p_ph_hcochob)      = real(max( ph(p_ph_hcochob     ,k,i),0.)/60., kind=dp)
           jv(p_ph_n2o5)         = real(max( ph(p_ph_n2o5        ,k,i),0.)/60., kind=dp)
           ! var
           c(ind_h2so4)    = real(max(chem(p_h2so4   ,k,i),0.), kind=dp) *conv
           c(ind_hcooh)    = real(max(chem(p_hcooh   ,k,i),0.), kind=dp) *conv
           c(ind_cco_oh)   = real(max(chem(p_cco_oh  ,k,i),0.), kind=dp) *conv
           c(ind_rco_oh)   = real(max(chem(p_rco_oh  ,k,i),0.), kind=dp) *conv
           c(ind_cco_ooh)  = real(max(chem(p_cco_ooh ,k,i),0.), kind=dp) *conv
           c(ind_rco_ooh)  = real(max(chem(p_rco_ooh ,k,i),0.), kind=dp) *conv
           c(ind_xn)       = real(max(chem(p_xn      ,k,i),0.), kind=dp) *conv
           c(ind_xc)       = real(max(chem(p_xc      ,k,i),0.), kind=dp) *conv
           c(ind_so2)      = real(max(chem(p_so2     ,k,i),0.), kind=dp) *conv
           c(ind_o1d)      = real(max(o1d(k,i),0.), kind=dp) *conv!real(max(chem(p_o1d     ,k,i),0.), kind=dp) *conv
           c(ind_c2h6)     = real(max(chem(p_c2h6    ,k,i),0.), kind=dp) *conv
           c(ind_bacl)     = real(max(chem(p_bacl    ,k,i),0.), kind=dp) *conv
           c(ind_pan)      = real(max(chem(p_pan     ,k,i),0.), kind=dp) *conv
           c(ind_pan2)     = real(max(chem(p_pan2    ,k,i),0.), kind=dp) *conv
           c(ind_pbzn)     = real(max(chem(p_pbzn    ,k,i),0.), kind=dp) *conv
           c(ind_ma_pan)   = real(max(chem(p_ma_pan  ,k,i),0.), kind=dp) *conv
           c(ind_h2o2)     = real(max(chem(p_h2o2    ,k,i),0.), kind=dp) *conv
           c(ind_n2o5)     = real(max(chem(p_n2o5    ,k,i),0.), kind=dp) *conv
           c(ind_hono)     = real(max(chem(p_hono    ,k,i),0.), kind=dp) *conv
           c(ind_c2h2)     = real(max(chem(p_c2h2    ,k,i),0.), kind=dp) *conv
           c(ind_alk3)     = real(max(chem(p_alk3    ,k,i),0.), kind=dp) *conv
           c(ind_tbu_o)    = real(max(tbu_o(k,i),0.), kind=dp) *conv!real(max(chem(p_tbu_o   ,k,i),0.), kind=dp) *conv
           c(ind_alk5)     = real(max(chem(p_alk5    ,k,i),0.), kind=dp) *conv
           c(ind_aro2)     = real(max(chem(p_aro2    ,k,i),0.), kind=dp) *conv
           c(ind_hno4)     = real(max(chem(p_hno4    ,k,i),0.), kind=dp) *conv
           c(ind_cooh)     = real(max(chem(p_cooh    ,k,i),0.), kind=dp) *conv
           c(ind_hocoo)    = real(max(hocoo(k,i),0.), kind=dp) *conv!real(max(chem(p_hocoo   ,k,i),0.), kind=dp) *conv
           c(ind_bzno2_o)  = real(max(bzno2_o(k,i),0.), kind=dp) *conv!real(max(chem(p_bzno2_o ,k,i),0.), kind=dp) *conv
           c(ind_meoh)     = real(max(chem(p_meoh    ,k,i),0.), kind=dp) *conv
           c(ind_alk4)     = real(max(chem(p_alk4    ,k,i),0.), kind=dp) *conv
           c(ind_aro1)     = real(max(chem(p_aro1    ,k,i),0.), kind=dp) *conv
           c(ind_dcb2)     = real(max(chem(p_dcb2    ,k,i),0.), kind=dp) *conv
           c(ind_dcb3)     = real(max(chem(p_dcb3    ,k,i),0.), kind=dp) *conv
           c(ind_cres)     = real(max(chem(p_cres    ,k,i),0.), kind=dp) *conv
           c(ind_dcb1)     = real(max(chem(p_dcb1    ,k,i),0.), kind=dp) *conv
           c(ind_nphe)     = real(max(chem(p_nphe    ,k,i),0.), kind=dp) *conv
           c(ind_rooh)     = real(max(chem(p_rooh    ,k,i),0.), kind=dp) *conv
           c(ind_bald)     = real(max(chem(p_bald    ,k,i),0.), kind=dp) *conv
           c(ind_phen)     = real(max(chem(p_phen    ,k,i),0.), kind=dp) *conv
           c(ind_co)       = real(max(chem(p_co      ,k,i),0.), kind=dp) *conv
           c(ind_mgly)     = real(max(chem(p_mgly    ,k,i),0.), kind=dp) *conv
           c(ind_acet)     = real(max(chem(p_acet    ,k,i),0.), kind=dp) *conv
           c(ind_hno3)     = real(max(chem(p_hno3    ,k,i),0.), kind=dp) *conv
           c(ind_ethene)   = real(max(chem(p_ethene  ,k,i),0.), kind=dp) *conv
           c(ind_gly)      = real(max(chem(p_gly     ,k,i),0.), kind=dp) *conv
           c(ind_bz_o)     = real(max(bz_o(k,i),0.), kind=dp) *conv!real(max(chem(p_bz_o    ,k,i),0.), kind=dp) *conv
           c(ind_isoprene) = real(max(chem(p_isoprene,k,i),0.), kind=dp) *conv
           c(ind_r2o2)     = real(max(chem(p_r2o2    ,k,i),0.), kind=dp) *conv
           c(ind_terp)     = real(max(chem(p_terp    ,k,i),0.), kind=dp) *conv
           c(ind_methacro) = real(max(chem(p_methacro,k,i),0.), kind=dp) *conv
           c(ind_ole1)     = real(max(chem(p_ole1    ,k,i),0.), kind=dp) *conv
           c(ind_isoprod)  = real(max(chem(p_isoprod ,k,i),0.), kind=dp) *conv
           c(ind_ole2)     = real(max(chem(p_ole2    ,k,i),0.), kind=dp) *conv
           c(ind_mvk)      = real(max(chem(p_mvk     ,k,i),0.), kind=dp) *conv
           c(ind_ccho)     = real(max(chem(p_ccho    ,k,i),0.), kind=dp) *conv
           c(ind_hcho)     = real(max(chem(p_hcho    ,k,i),0.), kind=dp) *conv
           c(ind_rno3)     = real(max(chem(p_rno3    ,k,i),0.), kind=dp) *conv
           c(ind_o3p)      = real(max(o3p(k,i),0.), kind=dp) *conv!real(max(chem(p_o3p     ,k,i),0.), kind=dp) *conv
           c(ind_rcho)     = real(max(chem(p_rcho    ,k,i),0.), kind=dp) *conv
           c(ind_mek)      = real(max(chem(p_mek     ,k,i),0.), kind=dp) *conv
           c(ind_prod2)    = real(max(chem(p_prod2   ,k,i),0.), kind=dp) *conv
           c(ind_o3)       = real(max(chem(p_o3      ,k,i),0.), kind=dp) *conv
           c(ind_ho2)      = real(max(chem(p_ho2     ,k,i),0.), kind=dp) *conv
           c(ind_ro2_n)    = real(max(chem(p_ro2_n   ,k,i),0.), kind=dp) *conv
           c(ind_ma_rco3)  = real(max(chem(p_ma_rco3 ,k,i),0.), kind=dp) *conv
           c(ind_c_o2)     = real(max(chem(p_c_o2    ,k,i),0.), kind=dp) *conv
           c(ind_bzco_o2)  = real(max(chem(p_bzco_o2 ,k,i),0.), kind=dp) *conv
           c(ind_ro2_r)    = real(max(chem(p_ro2_r   ,k,i),0.), kind=dp) *conv
           c(ind_no)       = real(max(chem(p_no      ,k,i),0.), kind=dp) *conv
           c(ind_no2)      = real(max(chem(p_no2     ,k,i),0.), kind=dp) *conv
           c(ind_no3)      = real(max(chem(p_no3     ,k,i),0.), kind=dp) *conv
           c(ind_cco_o2)   = real(max(chem(p_cco_o2  ,k,i),0.), kind=dp) *conv
           c(ind_rco_o2)   = real(max(chem(p_rco_o2  ,k,i),0.), kind=dp) *conv
           c(ind_oh)       = real(max(chem(p_ho      ,k,i),0.), kind=dp) *conv !oh in kpp, but ho in chem
           c(ind_c3h6)     = real(max(chem(p_c3h6    ,k,i),0.), kind=dp) *conv
           c(ind_sesq)     = real(max(chem(p_sesq    ,k,i),0.), kind=dp) *conv
           c(ind_co2)      = real(max(chem(p_co2     ,k,i),0.), kind=dp) *conv
           ! fix
           c(ind_ch4) = real(max(chem(p_ch4,k,i),0.), kind=dp) *conv
           c(ind_air) = real(dens2con_a * rho_phy(k,i), kind=dp) 
           c(ind_o2) = 0.20946e0 * c(ind_air) ! o2(oxygen) makes up approximately 21% of the air. 
           c(ind_h2o) = real(dens2con_w * qv(k,i) * rho_phy(k,i), kind=dp)
           c(ind_h2) = 5.40e-7 * c(ind_air) !h2(hydrogen )makes up approximately 0.000055% of the air.

           !call update_sun() ! if we want to use j-values, we need to adjust model.eqn @zhxia
           call update_rconst()

           !CALL GetMass( C, DVAL ) !can be used to check mass balance. @zhxia

           call integrate( tin       = tstart,   &
                           tout      = tend,     &
                           icntrl_u  = icntrl,   &
                           rcntrl_u  = rcntrl,   &
                           istatus_u = istate,   &
                           rstatus_u = rstate   )

           !CALL GetMass( C, DVAL ) !can be used to check mass balance. @zhxia

           chem(p_h2so4   ,k,i) = real(max(c(ind_h2so4),0.), kind   =sp) *oconv
           chem(p_hcooh   ,k,i) = real(max(c(ind_hcooh),0.), kind   =sp) *oconv
           chem(p_cco_oh  ,k,i) = real(max(c(ind_cco_oh),0.), kind  =sp) *oconv
           chem(p_rco_oh  ,k,i) = real(max(c(ind_rco_oh),0.), kind  =sp) *oconv
           chem(p_cco_ooh ,k,i) = real(max(c(ind_cco_ooh),0.), kind =sp) *oconv
           chem(p_rco_ooh ,k,i) = real(max(c(ind_rco_ooh),0.), kind =sp) *oconv
           chem(p_xn      ,k,i) = real(max(c(ind_xn),0.), kind      =sp) *oconv
           chem(p_xc      ,k,i) = real(max(c(ind_xc),0.), kind      =sp) *oconv
           chem(p_so2     ,k,i) = real(max(c(ind_so2),0.), kind     =sp) *oconv
            o1d(           k,i) = real(max(c(ind_o1d),0.), kind     =sp) *oconv
           chem(p_c2h6    ,k,i) = real(max(c(ind_c2h6),0.), kind    =sp) *oconv
           chem(p_bacl    ,k,i) = real(max(c(ind_bacl),0.), kind    =sp) *oconv
           chem(p_pan     ,k,i) = real(max(c(ind_pan),0.), kind     =sp) *oconv
           chem(p_pan2    ,k,i) = real(max(c(ind_pan2),0.), kind    =sp) *oconv
           chem(p_pbzn    ,k,i) = real(max(c(ind_pbzn),0.), kind    =sp) *oconv
           chem(p_ma_pan  ,k,i) = real(max(c(ind_ma_pan),0.), kind  =sp) *oconv
           chem(p_h2o2    ,k,i) = real(max(c(ind_h2o2),0.), kind    =sp) *oconv
           chem(p_n2o5    ,k,i) = real(max(c(ind_n2o5),0.), kind    =sp) *oconv
           chem(p_hono    ,k,i) = real(max(c(ind_hono),0.), kind    =sp) *oconv
           chem(p_c2h2    ,k,i) = real(max(c(ind_c2h2),0.), kind    =sp) *oconv
           chem(p_alk3    ,k,i) = real(max(c(ind_alk3),0.), kind    =sp) *oconv
           tbu_o(          k,i) = real(max(c(ind_tbu_o),0.), kind   =sp) *oconv
           chem(p_alk5    ,k,i) = real(max(c(ind_alk5),0.), kind    =sp) *oconv
           chem(p_aro2    ,k,i) = real(max(c(ind_aro2),0.), kind    =sp) *oconv
           chem(p_hno4    ,k,i) = real(max(c(ind_hno4),0.), kind    =sp) *oconv
           chem(p_cooh    ,k,i) = real(max(c(ind_cooh),0.), kind    =sp) *oconv
           hocoo(          k,i) = real(max(c(ind_hocoo),0.), kind   =sp) *oconv
           bzno2_o(        k,i) = real(max(c(ind_bzno2_o),0.), kind =sp) *oconv
           chem(p_meoh    ,k,i) = real(max(c(ind_meoh),0.), kind    =sp) *oconv
           chem(p_alk4    ,k,i) = real(max(c(ind_alk4),0.), kind    =sp) *oconv
           chem(p_aro1    ,k,i) = real(max(c(ind_aro1),0.), kind    =sp) *oconv
           chem(p_dcb2    ,k,i) = real(max(c(ind_dcb2),0.), kind    =sp) *oconv
           chem(p_dcb3    ,k,i) = real(max(c(ind_dcb3),0.), kind    =sp) *oconv
           chem(p_cres    ,k,i) = real(max(c(ind_cres),0.), kind    =sp) *oconv
           chem(p_dcb1    ,k,i) = real(max(c(ind_dcb1),0.), kind    =sp) *oconv
           chem(p_nphe    ,k,i) = real(max(c(ind_nphe),0.), kind    =sp) *oconv
           chem(p_rooh    ,k,i) = real(max(c(ind_rooh),0.), kind    =sp) *oconv
           chem(p_bald    ,k,i) = real(max(c(ind_bald),0.), kind    =sp) *oconv
           chem(p_phen    ,k,i) = real(max(c(ind_phen),0.), kind    =sp) *oconv
           chem(p_co      ,k,i) = real(max(c(ind_co),0.), kind      =sp) *oconv
           chem(p_mgly    ,k,i) = real(max(c(ind_mgly),0.), kind    =sp) *oconv
           chem(p_acet    ,k,i) = real(max(c(ind_acet),0.), kind    =sp) *oconv
           chem(p_hno3    ,k,i) = real(max(c(ind_hno3),0.), kind    =sp) *oconv
           chem(p_ethene  ,k,i) = real(max(c(ind_ethene),0.), kind  =sp) *oconv
           chem(p_gly     ,k,i) = real(max(c(ind_gly),0.), kind     =sp) *oconv
           bz_o(           k,i) = real(max(c(ind_bz_o),0.), kind    =sp) *oconv
           chem(p_isoprene,k,i) = real(max(c(ind_isoprene),0.), kind=sp) *oconv
           chem(p_r2o2    ,k,i) = real(max(c(ind_r2o2),0.), kind    =sp) *oconv
           chem(p_terp    ,k,i) = real(max(c(ind_terp),0.), kind    =sp) *oconv
           chem(p_methacro,k,i) = real(max(c(ind_methacro),0.), kind=sp) *oconv
           chem(p_ole1    ,k,i) = real(max(c(ind_ole1),0.), kind    =sp) *oconv
           chem(p_isoprod ,k,i) = real(max(c(ind_isoprod),0.), kind =sp) *oconv
           chem(p_ole2    ,k,i) = real(max(c(ind_ole2),0.), kind    =sp) *oconv
           chem(p_mvk     ,k,i) = real(max(c(ind_mvk),0.), kind     =sp) *oconv
           chem(p_ccho    ,k,i) = real(max(c(ind_ccho),0.), kind    =sp) *oconv
           chem(p_hcho    ,k,i) = real(max(c(ind_hcho),0.), kind    =sp) *oconv
           chem(p_rno3    ,k,i) = real(max(c(ind_rno3),0.), kind    =sp) *oconv
            o3p(           k,i) = real(max(c(ind_o3p),0.), kind     =sp) *oconv
           chem(p_rcho    ,k,i) = real(max(c(ind_rcho),0.), kind    =sp) *oconv
           chem(p_mek     ,k,i) = real(max(c(ind_mek),0.), kind     =sp) *oconv
           chem(p_prod2   ,k,i) = real(max(c(ind_prod2),0.), kind   =sp) *oconv
           chem(p_o3      ,k,i) = real(max(c(ind_o3),0.), kind      =sp) *oconv
           chem(p_ho2     ,k,i) = real(max(c(ind_ho2),0.), kind     =sp) *oconv
           chem(p_ro2_n   ,k,i) = real(max(c(ind_ro2_n),0.), kind   =sp) *oconv
           chem(p_ma_rco3 ,k,i) = real(max(c(ind_ma_rco3),0.), kind =sp) *oconv
           chem(p_c_o2    ,k,i) = real(max(c(ind_c_o2),0.), kind    =sp) *oconv
           chem(p_bzco_o2 ,k,i) = real(max(c(ind_bzco_o2),0.), kind =sp) *oconv
           chem(p_ro2_r   ,k,i) = real(max(c(ind_ro2_r),0.), kind   =sp) *oconv
           chem(p_no      ,k,i) = real(max(c(ind_no),0.), kind      =sp) *oconv
           chem(p_no2     ,k,i) = real(max(c(ind_no2),0.), kind     =sp) *oconv
           chem(p_no3     ,k,i) = real(max(c(ind_no3),0.), kind     =sp) *oconv
           chem(p_cco_o2  ,k,i) = real(max(c(ind_cco_o2),0.), kind  =sp) *oconv
           chem(p_rco_o2  ,k,i) = real(max(c(ind_rco_o2),0.), kind  =sp) *oconv
           chem(p_ho      ,k,i) = real(max(c(ind_oh),0.), kind      =sp) *oconv !oh in kpp, but ho in chem
           chem(p_c3h6    ,k,i) = real(max(c(ind_c3h6),0.)  , kind  =sp) *oconv
           chem(p_sesq    ,k,i) = real(max(c(ind_sesq),0.)  , kind  =sp) *oconv
           chem(p_co2     ,k,i) = real(max(c(ind_co2),0.)   , kind  =sp) *oconv
       endif                           
    enddo
 enddo

 t = rstate(1)

 time = t

end subroutine KPP_saprc99_driver
end module module_KPP_saprc99
