!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Interface between iAMAS and chemistry schemes for all common varaiables in
! different chemistry processes.
! * Author: Jiawang Feng, Zihan Xia from AEMOL USTC
! Last update: jwfeng 11/2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!     iAMAS-v5.2:src/core_atmosphere/physics/iamas_atmphys_interface.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module utils
  use, intrinsic :: iso_c_binding
  use flcl_mod
  use flcl_dualview_mod

  implicit none
  private

  public kokkos_allocate_dualview_wb

  interface kokkos_allocate_dualview_wb
  module procedure kokkos_allocate_dv_r64_1d_wb
  module procedure kokkos_allocate_dv_r64_2d_wb
  module procedure kokkos_allocate_dv_r64_3d_wb
  module procedure kokkos_allocate_dv_r64_4d_wb
  end interface kokkos_allocate_dualview_wb

  contains

  subroutine kokkos_allocate_dv_r64_1d_wb(A, v_A, n_A, b1)
    real(flcl_dualview_r64_f_t), pointer, dimension(:), intent(inout) :: A
    type(dualview_r64_1d_t), intent(out) :: v_A
    character(len=*), intent(in) :: n_A
    integer, intent(in) :: b1(2)
    real(flcl_dualview_r64_f_t), pointer, dimension(:) :: tmp_A
    integer :: lb, ub

    lb = b1(1); ub = b1(2)
    if (ub < lb) error stop 'Invalid 1d bounds'
    call kokkos_allocate_dualview(tmp_A, v_A, n_A, int(ub-lb+1, c_size_t))
    A(lb:ub) => tmp_A(:)
  end subroutine kokkos_allocate_dv_r64_1d_wb

  subroutine kokkos_allocate_dv_r64_2d_wb(A, v_A, n_A, b1, b2)
    real(flcl_dualview_r64_f_t), pointer, dimension(:,:), intent(inout) :: A
    type(dualview_r64_2d_t), intent(out) :: v_A
    character(len=*), intent(in) :: n_A
    integer, intent(in) :: b1(2),b2(2)
    real(flcl_dualview_r64_f_t), pointer, contiguous, dimension(:,:) :: tmp_A
    integer :: lb1,ub1,lb2,ub2

    lb1=b1(1); ub1=b1(2)
    lb2=b2(1); ub2=b2(2)
    if (ub1 < lb1 .or. ub2 < lb2) error stop 'Invalid 2d bounds'
    call kokkos_allocate_dualview(tmp_A, v_A, n_A, int(ub1-lb1+1,c_size_t), int(ub2-lb2+1,c_size_t))
    A(lb1:ub1,lb2:ub2) => tmp_A(:,:)
  end subroutine kokkos_allocate_dv_r64_2d_wb

  subroutine kokkos_allocate_dv_r64_3d_wb(A, v_A, n_A, b1, b2, b3)
    real(flcl_dualview_r64_f_t), pointer, dimension(:,:,:), intent(inout) :: A
    type(dualview_r64_3d_t), intent(out) :: v_A
    character(len=*), intent(in) :: n_A
    integer, intent(in) :: b1(2),b2(2),b3(2)
    real(flcl_dualview_r64_f_t), pointer, contiguous, dimension(:,:,:) :: tmp_A
    integer :: lb1,ub1,lb2,ub2,lb3,ub3

    lb1=b1(1); ub1=b1(2)
    lb2=b2(1); ub2=b2(2)
    lb3=b3(1); ub3=b3(2)
    if (ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3) error stop 'Invalid 3d bounds'
    call kokkos_allocate_dualview(tmp_A, v_A, n_A, int(ub1-lb1+1,c_size_t), int(ub2-lb2+1,c_size_t), int(ub3-lb3+1,c_size_t))
    A(lb1:ub1,lb2:ub2,lb3:ub3) => tmp_A(:,:,:)
  end subroutine kokkos_allocate_dv_r64_3d_wb

  subroutine kokkos_allocate_dv_r64_4d_wb(A, v_A, n_A, b1, b2, b3, b4)
    real(flcl_dualview_r64_f_t), pointer, dimension(:,:,:,:), intent(inout) :: A
    type(dualview_r64_4d_t), intent(out) :: v_A
    character(len=*), intent(in) :: n_A
    integer, intent(in) :: b1(2),b2(2),b3(2),b4(2)
    real(flcl_dualview_r64_f_t), pointer, contiguous, dimension(:,:,:,:) :: tmp_A
    integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4

    lb1=b1(1); ub1=b1(2)
    lb2=b2(1); ub2=b2(2)
    lb3=b3(1); ub3=b3(2)
    lb4=b4(1); ub4=b4(2)
    if (ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 .or. ub4 < lb4) error stop 'Invalid 4d bounds'
    call kokkos_allocate_dualview(tmp_A, v_A, n_A, int(ub1-lb1+1,c_size_t), &
    & int(ub2-lb2+1,c_size_t), int(ub3-lb3+1,c_size_t), int(ub4-lb4+1,c_size_t))
    A(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) => tmp_A(:,:,:,:)
  end subroutine kokkos_allocate_dv_r64_4d_wb
end module
MODULE iamas_atmchem_interface
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_constants
 use iamas_atmphys_vars
 use iamas_atmchem_utilities

 implicit none
 private
 public:: allocate_forall_chemistry,   &
          deallocate_forall_chemistry, &
          chemistry_from_iAMAS

! iAMAS interface for all chemistry.

! subroutines:
! ---------------------------------------
! allocate_forall_chemistry  : allocate local arrays for all chemistry.
! deallocate_forall_chemistry: deallocate local arrays for all chemistry.
! chemistry_from_iAMAS        : initialize local arrays.
! chemistry_to_iAMAS          : copy local arrays to iAMAS arrays.
! --------------------------------------

 contains

!=================================================================================================================
 subroutine allocate_forall_chemistry
!=================================================================================================================
    use utils
    if(.not.allocated(           dx) )    allocate(           dx(its:ite) )
    !if(.not.associated(          dx) )    call kokkos_allocate_dualview_wb(dx, dv_dx, "dx", [its,ite])
    if(.not.allocated(         xlat) )    allocate(         xlat(its:ite) )
    if(.not.allocated(        xlong) )    allocate(        xlong(its:ite) )
    if(.not.allocated(         psfc) )    allocate(         psfc(its:ite) )
    if(.not.allocated(     psfc_hyd) )    allocate(     psfc_hyd(its:ite) )
    if(.not.allocated(    psfc_hydd) )    allocate(    psfc_hydd(its:ite) )
    if(.not.allocated(           qv) )    allocate(           qv(kts:kte, its:ite) )
    if(.not.allocated(           qc) )    allocate(           qc(kts:kte, its:ite) )
    if(.not.allocated(           qi) )    allocate(           qi(kts:kte, its:ite) )
    if(.not.allocated(           qs) )    allocate(           qs(kts:kte, its:ite) )
    if(.not.allocated(           th) )    allocate(           th(kts:kte, its:ite) )
    if(.not.allocated(        t_phy) )    allocate(        t_phy(kts:kte, its:ite) )
    !if(.not.associated(        t_phy))    call kokkos_allocate_dualview_wb(t_phy, dv_t_phy, "t_phy", [kts,kte], [its,ite])
    if(.not.allocated(         pres) )    allocate(         pres(kts:kte, its:ite) )
    if(.not.allocated(      rho_phy) )    allocate(      rho_phy(kts:kte, its:ite) )
    !if(.not.associated(     rho_phy) )    call kokkos_allocate_dualview_wb(rho_phy, dv_rho_phy, "rho_phy", [kts,kte], [its,ite])
    if(.not.allocated(          alt) )    allocate(          alt(kts:kte, its:ite) )
    if(.not.allocated(           dz) )    allocate(           dz(kts:kte, its:ite) )
    !if(.not.associated(          dz) )    call kokkos_allocate_dualview_wb(dz, dv_dz, "dz", [kts,kte], [its,ite])
    if(.not.allocated(     pres_hyd) )    allocate(     pres_hyd(kts:kte, its:ite) )
    !if(.not.associated(    pres_hyd) )    call kokkos_allocate_dualview_wb(pres_hyd, dv_pres_hyd, "pres_hyd", [kts,kte], [its,ite])
    if(.not.allocated(    pres_hydd) )    allocate(    pres_hydd(kts:kte, its:ite) )
    if(.not.allocated(      znu_hyd) )    allocate(      znu_hyd(kts:kte, its:ite) )
    if(.not.allocated(           t2) )    allocate(           t2(kts:kte+1, its:ite) )
    if(.not.allocated(        pres2) )    allocate(        pres2(kts:kte+1, its:ite) )
    if(.not.allocated(    pres2_hyd) )    allocate(    pres2_hyd(kts:kte+1, its:ite) )
    !if(.not.associated(   pres2_hyd) )    call kokkos_allocate_dualview_wb(pres2_hyd, dv_pres2_hyd, "pres2_hyd", [kts,kte+1], [its,ite])
    if(.not.allocated(   pres2_hydd) )    allocate(   pres2_hydd(kts:kte+1, its:ite) )
    if(.not.allocated(        moist) )    allocate(        moist(1:num_moist, kts:kte, its:ite) )
    !if(.not.associated(       moist) )    call kokkos_allocate_dualview_wb(moist, dv_moist, "moist", [1,num_moist], [kts,kte], [its,ite])
    if(.not.allocated(        tke_pbl) )  allocate(        tke_pbl(kts:kte, its:ite) )
    if(.not.allocated(        f_ice_local) )  allocate(        f_ice_local(kts:kte, its:ite) )
    if(.not.allocated(        f_rain_local) )  allocate(        f_rain_local(kts:kte, its:ite) )
    if(.not.allocated(     chem_old) )    allocate(     chem_old(1:num_chem, kts:kte, its:ite) )
    if(.not.allocated(   gas_aqfrac) )    allocate(   gas_aqfrac(1:num_gas, kts:kte, its:ite) )

    gas_aqfrac(1:num_gas,kts:kte,its:ite) = 0.0_RKIND


 end subroutine allocate_forall_chemistry

!=================================================================================================================
 subroutine deallocate_forall_chemistry
!=================================================================================================================
    if(allocated(           dx) )   deallocate(           dx)
    !if(associated(          dx) )   call kokkos_deallocate_dualview(dx, dv_dx)
    if(allocated(         xlat) )   deallocate(         xlat)
    if(allocated(        xlong) )   deallocate(        xlong)
    if(allocated(     psfc_hyd) )   deallocate(     psfc_hyd)
    if(allocated(    psfc_hydd) )   deallocate(    psfc_hydd)
    if(allocated(      qndrop3d) )   deallocate(    qndrop3d)
    if(allocated(      qnice3d) )   deallocate(      qnice3d)
    if(allocated(      nimey3d) )   deallocate(      nimey3d)
    if(allocated(           qv) )   deallocate(           qv)
    if(allocated(           qc) )   deallocate(           qc)
    if(allocated(           qi) )   deallocate(           qi)
    if(allocated(           qs) )   deallocate(           qs)
    if(allocated(           th) )   deallocate(           th)
    if(allocated(        t_phy) )   deallocate(        t_phy)
    !if(associated(       t_phy) )   call kokkos_deallocate_dualview(t_phy, dv_t_phy)
    if(allocated(      tke_pbl) )   deallocate(      tke_pbl)
    if(allocated(  f_ice_local) )   deallocate(  f_ice_local)
    if(allocated( f_rain_local) )   deallocate( f_rain_local)
    if(allocated(      rho_phy) )   deallocate(      rho_phy)
    !if(associated(     rho_phy) )   call kokkos_deallocate_dualview(rho_phy, dv_rho_phy)
    if(allocated(          alt) )   deallocate(          alt)
    if(allocated(           dz) )   deallocate(           dz)
    !if(associated(          dz) )   call kokkos_deallocate_dualview(dz, dv_dz)
    if(allocated(         pres) )   deallocate(         pres)
    if(allocated(     pres_hyd) )   deallocate(     pres_hyd)
    !if(associated(    pres_hyd) )   call kokkos_deallocate_dualview(pres_hyd, dv_pres_hyd)
    if(allocated(    pres_hydd) )   deallocate(    pres_hydd)
    if(allocated(      znu_hyd) )   deallocate(      znu_hyd)
    if(allocated(           t2) )   deallocate(           t2)
    if(allocated(        pres2) )   deallocate(        pres2)
    if(allocated(    pres2_hyd) )   deallocate(    pres2_hyd)
    !if(associated(   pres2_hyd) )   call kokkos_deallocate_dualview(pres2_hyd, dv_pres2_hyd)
    if(allocated(   pres2_hydd) )   deallocate(   pres2_hydd)
    if(allocated(        moist) )   deallocate(        moist)
    !if(associated(       moist) )   call kokkos_deallocate_dualview(moist, dv_moist)
    if(allocated(     chem_old) )   deallocate(     chem_old)
    if(allocated(   gas_aqfrac) )   deallocate(   gas_aqfrac)

 end subroutine deallocate_forall_chemistry

!=================================================================================================================
 subroutine chemistry_from_iAMAS(configs,mesh,state,time_lev,sfc_input,diag,diag_physics,tend_chemistry,its,ite)
!=================================================================================================================

!input variables:
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: state
 type(iamas_pool_type),intent(in):: sfc_input
 type(iamas_pool_type),intent(in):: diag_physics
 type(iamas_pool_type),intent(in):: diag
 type(iamas_pool_type),intent(in):: tend_chemistry
 integer,intent(in):: its,ite,time_lev
 real(kind=RKIND),dimension(:)    ,pointer :: latCell,lonCell

!local variables:
 integer :: i,k,j,h,n
 real(kind=RKIND) :: z0,z1,z2,w1,w2,fzm,fzp
 real(kind=RKIND) :: rho_a,rho1,rho2,tem1,tem2
 character(len=StrKIND) :: CHEM_MSG


!-----------------------------------------------------------------------------------------------------------------
!Jiawang Feng:
!Gocart use a constantly dxy as area.Can we use array meshDensity to reach the same results?
!jwfeng: 2022/08/11, Use len_disp/meshDensity**0.25 instead of areaCell**0.5 to cal dx
 call iamas_pool_get_config(configs,'config_len_disp',len_disp)
 call iamas_pool_get_array(mesh,'meshDensity',meshDensity)
 call iamas_pool_get_array(mesh,'latCell' ,latCell    )
 call iamas_pool_get_array(mesh,'lonCell' ,lonCell   )
 call iamas_pool_get_array(mesh,'zgrid'   ,zgrid   )
 call iamas_pool_get_array(mesh,'zz'      ,zz      )
 call iamas_pool_get_array(mesh,'rdzw'    ,rdzw    )

 call iamas_pool_get_array(diag_physics,'u10'         ,u10           )
 call iamas_pool_get_array(diag_physics,'v10'         ,v10           )
 call iamas_pool_get_array(diag_physics,'ust'         ,ust           )
 call iamas_pool_get_array(diag_physics,'gsw'         ,gsw           )
 call iamas_pool_get_array(diag_physics,'hfx'         ,hfx           )
 call iamas_pool_get_array(diag_physics,'znt'         ,znt           )
 call iamas_pool_get_array(diag_physics,'rmol'        ,rmol          )
 call iamas_pool_get_array(diag_physics,'hpbl'        ,hpbl          )
 call iamas_pool_get_array(diag_physics,'raincv'      ,raincv        )
 call iamas_pool_get_array(diag_physics,'kzh'         ,kzh           )
 call iamas_pool_get_array(diag_physics,'cldfrac'     ,cldfrac       )
 call iamas_pool_get_array(diag_physics,'cldfrac_old' ,cldfrac_old   )
 call iamas_pool_get_array(diag_physics,'lai'         ,lai           )

 call iamas_pool_get_array(diag,'surface_pressure'      ,surface_pressure)
 call iamas_pool_get_array(diag,'exner'                 ,exner           )
 call iamas_pool_get_array(diag,'uReconstructZonal'     ,u_phy           )
 call iamas_pool_get_array(diag,'uReconstructMeridional',v_phy           )
 call iamas_pool_get_array(diag,'pressure_base'         ,pressure_b      )
 call iamas_pool_get_array(diag,'pressure_p'            ,pressure_p      )

 call iamas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call iamas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call iamas_pool_get_array(state,'w      ',w      ,time_lev)
 call iamas_pool_get_array(state,'scalars',scalars,time_lev)
 call iamas_pool_get_array(state,'chem'   ,chem   ,time_lev)

 call iamas_pool_get_array(sfc_input,'xland'    ,xland    )
 call iamas_pool_get_array(sfc_input,'ivgtyp'   ,ivgtyp   )
 call iamas_pool_get_array(sfc_input,'isltyp'   ,isltyp   )
 call iamas_pool_get_array(sfc_input,'skintemp' ,skintemp )
 call iamas_pool_get_array(sfc_input,'vegfra'   ,vegfra   )
 call iamas_pool_get_array(sfc_input,'snowh'    ,snowh    )
 call iamas_pool_get_array(sfc_input,'snowc'    ,snowc    )
 call iamas_pool_get_array(sfc_input,'soilctop' ,soilctop )
 call iamas_pool_get_array(sfc_input,'landusef' ,landusef )
 call iamas_pool_get_array(sfc_input,'sh2o'     ,sh2o     )

 call iamas_pool_get_dimension(state,'index_qv' ,index_qv )
 call iamas_pool_get_dimension(state,'index_qc' ,index_qc )
 call iamas_pool_get_dimension(state,'index_qi' ,index_qi )
 call iamas_pool_get_dimension(state,'index_qs' ,index_qs )
 qv = scalars(index_qv,kts:kte,its:ite)
 qc = scalars(index_qc,kts:kte,its:ite)
 qi = scalars(index_qi,kts:kte,its:ite)
 qs = scalars(index_qs,kts:kte,its:ite)
 qndrop3d = scalars(p_qndrop,kts:kte,its:ite)
 qnice3d = scalars(p_qniceaer,kts:kte,its:ite)
 nimey3d = scalars(p_qnicemix,kts:kte,its:ite)
 moist = scalars(1:num_moist,kts:kte,its:ite)

 if(associated(meshDensity      )  .and. associated(latCell         )  .and. &
    associated(lonCell          )  .and. associated(zgrid           )  .and. &
    associated(zz               )  .and. associated(rdzw            )  .and. &
    associated(u10              )  .and. associated(v10             )  .and. &
    associated(ust              )  .and. associated(gsw             )  .and. &
    associated(hfx              )  .and. associated(znt             )  .and. &   
    associated(rmol             )  .and. associated(hpbl            )  .and. &  
    associated(raincv           )  .and. associated(surface_pressure)  .and. & 
    associated(exner            )  .and. associated(u_phy           )  .and. &
    associated(v_phy            )  .and. associated(pressure_b      )  .and. &
    associated(pressure_p       )  .and. associated(rho_zz          )  .and. &
    associated(theta_m          )  .and. associated(w               )  .and. &
    associated(scalars          )  .and. associated(chem            )  .and. &
    associated(xland            )  .and. associated(ivgtyp          )  .and. &
    associated(isltyp           )  .and. associated(skintemp        )  .and. &
    associated(vegfra           )  .and. associated(snowh           )  .and. &
    associated(cldfrac          )  .and. associated(cldfrac_old     )  .and. &
    associated(kzh              )  .and. associated(snowc           )  .and. &
    associated(soilctop         )  .and. associated(landusef        )  .and. &
    associated(lai              )  .and. associated(sh2o            )        & 
    ) then
 else
    call chemistry_error_fatal('Something wrong in chemistry interface:some arguments are not correctly associated!!!!!')
 endif

 !calculation of the surface pressure using hydrostatic assumption down to the surface::
 do i = its,ite
                   tem1 = zgrid(2,i)-zgrid(1,i)
                   tem2 = zgrid(3,i)-zgrid(2,i)
                   rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv(1,i))
                   rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv(2,i))
    surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                             *(rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo

!   scalars(1:num_scalars,kts:kte,its:ite) = max(0._RKIND, scalars(1:num_scalars,kts:kte,its:ite))
!   moist(1:num_moist,kts:kte,its:ite) = max(0._RKIND, moist(1:num_moist,kts:kte,its:ite))
!   chem(1:num_chem,kts:kte,its:ite) = max(epsilc,chem(1:num_chem,kts:kte,its:ite))
!   qv(kts:kte,its:ite) = max(0.,qv(kts:kte,its:ite))

 rho_phy(kts:kte,its:ite) = zz(kts:kte,its:ite) * rho_zz(kts:kte,its:ite) * (1._RKIND + qv(kts:kte,its:ite))
 th(kts:kte,its:ite) = theta_m(kts:kte,its:ite) / (1._RKIND + R_v/R_d * qv(kts:kte,its:ite))
 t_phy(kts:kte,its:ite) = th(kts:kte,its:ite) * exner(kts:kte,its:ite)
 alt(kts:kte,its:ite) = 1._RKIND / rho_phy(kts:kte,its:ite)
 pres(kts:kte,its:ite) = pressure_p(kts:kte,its:ite) + pressure_b(kts:kte,its:ite)
 dz(kts:kte,its:ite) = zgrid(kts+1:kte+1,its:ite) - zgrid(kts:kte,its:ite) 
 f_ice_local(kts:kte,its:ite)  = 0.0_RKIND
 f_rain_local(kts:kte,its:ite) = 0.0_RKIND
 tke_pbl(kts:kte,its:ite) = 0.0_RKIND

 xlat(its:ite) = latCell(its:ite) * 57.2957795
 xlong(its:ite) = lonCell(its:ite) * 57.2957795
 dx(its:ite) = len_disp / (meshDensity(its:ite)**0.25)

!interpolation of pressure and temperature from theta points to w points:
 do i = its,ite
 do k = kts+1, kte
              tem1 = 1./(zgrid(k+1,i)-zgrid(k-1,i))
               fzm = (zgrid(k,i)-zgrid(k-1,i)) * tem1
               fzp = (zgrid(k+1,i)-zgrid(k,i)) * tem1
       t2(k,i) = fzm*t_phy(k,i) + fzp*t_phy(k-1,i)
    pres2(k,i) = fzm*pres(k,i)  + fzp*pres(k-1,i)
 end do
 end do

!interpolation of pressure and temperature from theta points to the top-of-the-model: follows
!the calculation of the top-of-the-model pressure and temperature in WRF (subroutine phy_prep
!in ./dyn_em/module_big_step_utilities.F):
 k = kte+1
 do i = its,ite
                z0 = zgrid(k,i)
                z1 = 0.5*(zgrid(k,i)+zgrid(k-1,i)) 
                z2 = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
                w1 = (z0-z2)/(z1-z2)
                w2 = 1.-w1
       t2(k,i) = w1*t_phy(k-1,i) + w2*t_phy(k-2,i)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2(k,i) = exp(w1*log(pres(k-1,i))+w2*log(pres(k-2,i)))
 end do

 !ldf (2012-06-22): recalculates the pressure at the surface as an extrapolation of the
 !pressures in the 2 layers above the surface, as was originally done:
 k = kts
 do i = its,ite
                z0 = zgrid(k,i)
                z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i)) 
                z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
                w1 = (z0-z2)/(z1-z2)
                w2 = 1.-w1
       t2(k,i) = w1*t_phy(k,i) + w2*t_phy(k+1,i)
    pres2(k,i) = w1* pres(k,i) + w2* pres(k+1,i)
       psfc(i) = surface_pressure(i)
 enddo

 !calculation of the hydrostatic pressure:
 do i = its,ite
    !pressure at w-points:
    k = kte+1
    pres2_hyd(k,i)  = pres2(k,i)
    pres2_hydd(k,i) = pres2(k,i)
    do k = kte,1,-1
       rho_a = rho_phy(k,i) / (1.+qv(k,i))
       pres2_hyd(k,i)  = pres2_hyd(k+1,i)  + gravity*rho_phy(k,i)*dz(k,i)
       pres2_hydd(k,i) = pres2_hydd(k+1,i) + gravity*rho_a*dz(k,i)
    enddo
    !pressure at theta-points:
    do k = kte,1,-1
       pres_hyd(k,i)  = 0.5*(pres2_hyd(k+1,i)+pres2_hyd(k,i))
       pres_hydd(k,i) = 0.5*(pres2_hydd(k+1,i)+pres2_hydd(k,i))
    enddo
    !surface pressure:
     psfc_hyd(i) = pres2_hyd(1,i)
    psfc_hydd(i) = pres2_hydd(1,i)
    !znu:
    do k = kte,1,-1
       znu_hyd(k,i) = pres_hyd(k,i) / psfc_hyd(i) 
    enddo
 enddo

 end subroutine chemistry_from_iAMAS

END MODULE iamas_atmchem_interface

