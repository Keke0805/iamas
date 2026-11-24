!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to Compute Convective Transport following KFeta scheme
! * Author: Jiawang Feng and Chun Zhao
! Last update: jwfeng 11/2021
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

module module_ctrans_general

   use iamas_atmchem_utilities, only: chemistry_message, chemistry_error_fatal
   use module_cu_gd, only: cup_dd_edt, cup_dd_he, cup_dd_moisture      , &
                           cup_dd_nms, cup_env, cup_env_clev, cup_kbcon, &
                           cup_ktop, cup_maximi, cup_minimi, cup_up_he , &
                           cup_up_moisture, cup_up_nms
   use iamas_atmchem_vars, only: p_co, p_so2, p_hno3, p_hno4, p_dms     , &
                                p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1    , &
                                p_op1, p_paa, p_h2so4, p_so4aj, p_nh4aj , &
                                p_no3aj, p_bc1, p_bc2, p_oc1, p_oc2    , &
                                p_seas_1, p_seas_2, p_seas_3, p_seas_4
   use iamas_log
   use iamas_kind_types
   use iamas_timer
   implicit none

   interface
      subroutine f_cup_ct(ktop,kbcon,k23,kbcon3,ktop3,xmb3,tracer,j,aaeq,t,q,z1,  &
            pre,p,tracert,tracert3,hstary,dtime,psur,us,vs,tcrit,    &
            wetdep_3d,rho_phy,dz8w,xl,rv,cp,g,ipr,jpr,npr,num_chem,chem_opt,         &
            cu_scav_on,conv_tr_aqchem, chem, aq_aer_ratio,           &
            ishallow,num_gas,                                        &
            its,ite, kts,kte,                                        &
            nsize_aer, ncomp_aer, massptr_aer, dens_aer,             &
            mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer, &
            mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer,&
            mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer, &
            p_co, p_so2, p_hno3, p_hno4, p_dms,p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1 , &
            p_op1, p_paa, p_h2so4, p_so4aj, p_nh4aj , p_no3aj) bind(C,name="c_cup_ct")
         use, intrinsic :: iso_c_binding
         use :: flcl_mod
         implicit none
         integer(c_int) :: j,ipr,jpr,npr,num_chem
         real(c_double) :: dtime,tcrit,xl,rv,cp,g
         integer(c_int) :: conv_tr_aqchem,ishallow,num_gas
         integer(c_int) :: its, ite, kts, kte
         character(c_char),intent(in), dimension(80) :: chem_opt
         type(nd_array_t) :: ktop, kbcon, k23, kbcon3, ktop3
         type(nd_array_t) :: xmb3,tracer,aaeq,t,q,z1,pre
         type(nd_array_t) :: p,tracert,tracert3,hstary,psur,us,vs
         type(nd_array_t) :: wetdep_3d,rho_phy,dz8w,aq_aer_ratio
         logical(c_bool), intent(in), value :: cu_scav_on
         type(nd_array_t) :: chem
         type(nd_array_t) :: nsize_aer, ncomp_aer, massptr_aer, dens_aer
         integer(c_int) :: mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer
         integer(c_int) :: mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer
         integer(c_int) :: mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer
         integer(c_int) :: p_co, p_so2, p_hno3, p_hno4, p_dms, p_n2o5
         integer(c_int) :: p_nh3, p_h2o2, p_o3, p_ora1, p_op1, p_paa
         integer(c_int) :: p_h2so4, p_so4aj, p_nh4aj, p_no3aj
      end subroutine f_cup_ct  
   end interface

   private

   real, parameter :: qcldwtr_cutoff = 1.0e-6 ! kg / kg
   real, parameter :: mwdry = 28.966  ! molecular mass of dry air (g/mol)
   real, parameter :: mwso4 = 96.00   ! molecular mass of so4-- (g/mol)
   real, parameter :: mwno3 = 62.0    ! molecular mass of no3- (g/mol)
   real, parameter :: mwnh4 = 18.0985 ! molecular mass of nh4+ (g/mol)

   public :: general_cup_transport

contains
   subroutine c_ctrans_aqchem_bridge(temp, pres_pa, taucld, prcrate, &
      wcavg, wtavg, airm, alfa0, alfa2, alfa3, gas, aerosol, liquid, &
      gaswdep, aerwdep, hpwdep) bind(C, name="c_ctrans_aqchem_bridge")
   use, intrinsic :: iso_c_binding
   use module_ctrans_aqchem, only: aqchem
   implicit none

   real(c_double) :: temp, pres_pa, taucld, prcrate
   real(c_double) :: wcavg, wtavg, airm, alfa0, alfa2, alfa3, hpwdep
   real(c_double), dimension(12) :: gas, gaswdep
   real(c_double), dimension(36) :: aerosol, aerwdep
   real(c_double), dimension(41) :: liquid

   call aqchem(temp, pres_pa, taucld, prcrate, wcavg, wtavg, airm, &
      alfa0, alfa2, alfa3, gas, aerosol, liquid, gaswdep, aerwdep, hpwdep)
   end subroutine c_ctrans_aqchem_bridge

   subroutine cup_ct_flcl(ktop,kbcon,k23,kbcon3,ktop3,xmb3,tracer,j,aaeq,t,q,z1,  &
      pre,p,tracert,tracert3,hstary,dtime,psur,us,vs,tcrit,    &
      wetdep_3d,rho_phy,dz8w,xl,rv,cp,g,ipr,jpr,npr,num_chem,chem_opt,         &
      cu_scav_on,conv_tr_aqchem, chem, aq_aer_ratio,           &
      ishallow,num_gas,                                        &
      its,ite, kts,kte)
   use, intrinsic :: iso_c_binding
   use :: flcl_mod

   use module_data_binaer_asect, only: nsize_aer, ncomp_aer, massptr_aer, dens_aer, &
          mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer, &
          mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer, &
          mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer

   use iamas_atmchem_vars, only: p_co, p_so2, p_hno3, p_hno4, p_dms     , &
                     p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1    , &
                     p_op1, p_paa, p_h2so4, p_so4aj, p_nh4aj , &
                     p_no3aj, p_bc1, p_bc2, p_oc1, p_oc2    , &
                     p_seas_1, p_seas_2, p_seas_3, p_seas_4

   implicit none

   integer(c_int) :: j,ipr,jpr,npr,num_chem
   real(c_double) :: dtime,tcrit,xl,rv,cp,g
   integer(c_int) :: conv_tr_aqchem,ishallow,num_gas
   integer(c_int) :: its, ite, kts, kte
   character(len=80), intent(in) :: chem_opt

   integer(c_int), dimension (its:ite), intent(inout) :: ktop, kbcon, k23, kbcon3, ktop3
   real(c_double), dimension (its:ite), intent (in) :: xmb3,aaeq,z1,psur
   real(c_double), dimension (its:ite), intent (inout) :: pre
   real(c_double), dimension(kts:kte,its:ite), intent(in) :: t,p,hstary,us,vs
   real(c_double), dimension(kts:kte,its:ite), intent(inout) :: q,aq_aer_ratio
   real(c_double), dimension (1:num_chem, kts:kte, its:ite), intent(inout) :: tracert
   real(c_double), dimension (1:num_chem, kts:kte, its:ite), intent(inout) :: tracer
   real(c_double), dimension (1:num_chem, kts:kte, its:ite), intent(inout) :: tracert3
   real(c_double), dimension (1:num_chem, kts:kte, its:ite), intent(inout) :: chem
   real(c_double), dimension(1:num_chem,kts:kte,its:ite), intent(out) :: wetdep_3d
   real(c_double), dimension(kts:kte,its:ite), intent(in) :: dz8w
   real(c_double), dimension(kts:kte,its:ite), intent(in) :: rho_phy
   logical, intent(in) :: cu_scav_on
   logical(c_bool) :: cu_scav_on_c

   ! integer(c_int), dimension (maxd_atype) :: nsize_aer
   ! integer(c_int), dimension (maxd_atype) :: ncomp_aer
   ! integer(c_int), dimension (maxd_acomp, maxd_asize, maxd_atype, maxd_aphase) :: massptr_aer
   ! real(c_double), dimension (maxd_acomp, maxd_atype) :: dens_aer

   ! integer(c_int) :: mastercompindx_seas_aer,mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer,&
   !                   mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer
   
   ! integer(c_int) :: p_co, p_so2, p_hno3, p_hno4, p_dms,       &
   !                   p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1,      &
   !                   p_op1, p_paa, p_h2so4, p_so4aj, p_nh4aj,  &
   !                   p_no3aj, p_bc1, p_bc2, p_oc1, p_oc2    ,  &
   !                   p_seas_1, p_seas_2, p_seas_3, p_seas_4

   cu_scav_on_c = cu_scav_on

   call f_cup_ct(to_nd_array(ktop),to_nd_array(kbcon),to_nd_array(k23),to_nd_array(kbcon3),to_nd_array(ktop3), &
      to_nd_array(xmb3),to_nd_array(tracer),j,to_nd_array(aaeq),to_nd_array(t),to_nd_array(q),to_nd_array(z1),  &
      to_nd_array(pre),to_nd_array(p),to_nd_array(tracert),to_nd_array(tracert3),to_nd_array(hstary),        &
      dtime,to_nd_array(psur),to_nd_array(us),to_nd_array(vs),tcrit,    &
      to_nd_array(wetdep_3d),to_nd_array(rho_phy),to_nd_array(dz8w),xl,rv,cp,g,ipr,jpr,npr,num_chem,chem_opt,         &
      cu_scav_on_c,conv_tr_aqchem, to_nd_array(chem), to_nd_array(aq_aer_ratio),         &
      ishallow,num_gas,its,ite, kts,kte , &         
      to_nd_array(nsize_aer),to_nd_array(ncomp_aer), to_nd_array(massptr_aer), to_nd_array(dens_aer),        &
      mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer, &
      mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer,&
      mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer, &
      p_co, p_so2, p_hno3, p_hno4, p_dms,p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1 , &
      p_op1, p_paa, p_h2so4, p_so4aj, p_nh4aj , p_no3aj)

   end subroutine cup_ct_flcl
   subroutine general_cup_transport(                  &
              dt, rho_phy, raincv, chem, wetdep_3d  , &
              u, v, t_phy, qv, dz8w, p_phy          , &
              xlv, cp, g, r_v, z                    , &
              tracert_3d, aq_aer_ratio              , &
              ktop_2d, kbcon_2d                     , &
              k22_shallow, kbcon_shallow            , &
              ktop_shallow, xmb_shallow, ishallow   , &
              num_gas, num_chem, chem_opt           , &
              cu_scav_on, conv_tr_aqchem            , &
              its, ite, kts, kte          )

      integer, intent(in) :: num_gas, ishallow, num_chem
      integer, intent(in) :: its, ite, kts, kte
      integer, intent(in) :: conv_tr_aqchem
      logical, intent(in) :: cu_scav_on
      character(len=*), intent(in) :: chem_opt

      real, intent(in) :: dt
      real, intent(in) :: xlv, r_v
      real, intent(in) :: cp, g
      real, dimension(kts:kte,its:ite), intent(in) :: qv
      real, dimension(kts:kte,its:ite), intent(in) :: u
      real, dimension(kts:kte,its:ite), intent(in) :: v
      real, dimension(kts:kte,its:ite), intent(in) :: t_phy
      real, dimension(kts:kte,its:ite), intent(in) :: z
      real, dimension(kts:kte,its:ite), intent(in) :: p_phy
      real, dimension(kts:kte,its:ite), intent(in) :: dz8w
      real, dimension(kts:kte,its:ite), intent(in) :: rho_phy
      real, dimension(its:ite), intent(in) :: raincv
      real, dimension(1:num_chem,kts:kte,its:ite), intent(inout) :: chem
      real, dimension(1:num_chem,kts:kte,its:ite), intent(out) :: wetdep_3d
      real, dimension(its:ite), intent(in), optional :: xmb_shallow

      integer, dimension(its:ite), intent(in), optional :: k22_shallow
      integer, dimension(its:ite), intent(in), optional :: kbcon_shallow
      integer, dimension(its:ite), intent(in), optional :: ktop_shallow
      integer, dimension(its:ite), intent(out), optional :: kbcon_2d
      integer, dimension(its:ite), intent(out), optional :: ktop_2d
      real, dimension(1:num_chem,kts:kte,its:ite), intent(out) :: tracert_3d
      real, dimension(kts:kte,its:ite), intent (inout) :: aq_aer_ratio ! Output to diagnosticss
      ! local vars
      real, dimension(kts:kte,its:ite) :: outt, outq, outqc
      real, dimension(its:ite) :: pret, ter11
      ! basic environmental input includes moisture convergence (mconv)
      ! omega (omeg), windspeed (us, vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      real, dimension(kts:kte,its:ite) :: T, TN, q, qo, PO, P, US, VS, hstary
      real, dimension(1:num_chem,kts:kte,its:ite) :: tracer, tracert, tracert3
      real, dimension(its:ite) :: Z1, PSUR, AAEQ, xmb3
      integer, dimension(its:ite) :: ktop, kbcon, k23, kbcon3, ktop3

      integer :: nv, i, j, k, icldck, ipr, jpr, npr
      integer :: itf, jtf, ktf, iopt
      real :: tcrit, dp, dq, epsilc
      character(len=300) :: message
      real :: tracert_ratio

      character(len=StrKIND):: debugString

      ! jwfeng: *_shallow vars should present if ishallow equals 1
      if (ishallow .eq. 1) then
         if (.not. (present(xmb_shallow)   .and. &
                    present(k22_shallow)   .and. &
                    present(kbcon_shallow) .and. &
                    present(ktop_shallow)) ) then
            call chemistry_error_fatal("*_shallow vars should present if ishallow equals 1")
         end if
      end if
      
      call iamas_timer_start("assignment")
      epsilc = 1.e-30
      wetdep_3d(:,:,:) = 0.0
      ktop_2d(:) = 0.0
      kbcon_2d(:) = 0.0
      tracert_3d(:,:,:) = 0.0
      ipr = 0
      jpr = 0
      npr = 1
      if (p_co .gt. 0) npr = p_co
      tcrit = 258.
      iopt = 0
      call iamas_timer_stop("assignment")

         call iamas_timer_start("calculate1")
         do i = its, ite
            ktop(i) = 0
            kbcon(i)= 0
            psur(i) = p_phy(kts, i) * .01
            ter11(i) = z(kts, i)
            aaeq(i) = 0.
            ! rainrate is input for this transport/wet-deposition routine
            pret(i) = raincv(i) / dt
            if (pret(i) .le. 0.) aaeq(i) = 20.
         enddo
         if (ishallow .eq. 1) then
            do i = its, ite
               xmb3(i) = xmb_shallow(i)
               ktop3(i) = ktop_shallow(i)
               k23(i) = k22_shallow(i)
               kbcon3(i) = kbcon_shallow(i)
            enddo
         else
            do i = its, ite
               xmb3(i) = 0.
               ktop3(i) = 0
               k23(i) = 0
               kbcon3(i) = 0
            enddo
         endif
         do i = its, ite
            do k = kts, kte
               po(k, i) = p_phy(k, i) * .01
               p(k, i) = po(k, i)
               us(k, i) = u(k, i)
               vs(k, i) = v(k, i)
               t(k, i) = t_phy(k, i)
               q(k, i) = qv(k, i)
               if (q(k, i) .lt. 1.e-08) q(k, i) = 1.e-08
            enddo
         enddo
         do i = its,ite
            do K=kts,kte
               do nv=1,num_chem
                  tracer(nv,k,i)=max(epsilc,chem(nv,k,i))
                  tracert(nv,k,i)=0.
                  tracert3(nv,k,i)=0.
               enddo
            enddo
         enddo
         do i=its,ite
            do k=kts,kte
               ! cu_co_ten(i,k,j)=0.
               ! hstary(i,k)=hstar4(nv)*exp(dhr(nv)*(1./t(i,k)-1./298.))
               if(i.eq.ipr)then
                  write(message, "(I3,x,4F9.3)") k, pret(i), tracer(npr,k, i), p(k, i), z(k, i)
                  call chemistry_message(trim(message))
               endif
            enddo
         enddo
         call iamas_timer_stop("calculate1")
         
         call iamas_timer_start("cup_ct")
         ! call cup_ct(ktop, kbcon, k23, kbcon3, ktop3, xmb3, &
         !      tracer, j, aaeq, t, q, ter11, pret, p, tracert, tracert3, &
         !      hstary, dt, psur, us, vs, tcrit, wetdep_1d, &
         !      xlv, r_v, cp, g, ipr, jpr, npr, num_chem, chem_opt, &
         !      cu_scav_on, conv_tr_aqchem, chem, aq_aer_ratio, &
         !      ishallow, num_gas, &
         !      its, ite, kts, kte)

         call cup_ct_flcl(ktop, kbcon, k23, kbcon3, ktop3, xmb3, &
              tracer, j, aaeq, t, q, ter11, pret, p, tracert, tracert3, &
              hstary, dt, psur, us, vs, tcrit, wetdep_3d, rho_phy, dz8w, &
              xlv, r_v, cp, g, ipr, jpr, npr, num_chem, chem_opt, &
              cu_scav_on, conv_tr_aqchem, chem, aq_aer_ratio, &
              ishallow, num_gas, &
              its, ite, kts, kte)
         call iamas_timer_stop("cup_ct")

         call iamas_timer_start("calculate2")
         do i = its, ite
            if (pret(i) .le. 0.) then
               do k = kts, kte
                  do nv = 1, num_chem
                     tracert(nv, k, i) = 0.
                  enddo
               enddo
            endif
         enddo

         if (ishallow .eq. 1) then
            call neg_check_ct('shallow', pret, ktop3, epsilc, dt, &
                 tracer, tracert3, iopt, num_chem, &
                 its, ite, kts, kte, ipr, jpr, npr, j)
            do i = its, ite
               do k = kts, kte
                  do nv = 1, num_chem
                     chem(nv, k, i) = max(epsilc, chem(nv, k, i) &
                                       + tracert3(nv, k, i) * dt)
                     tracert_3d(nv, k, i) = tracert_3d(nv, k, i) + tracert3(nv, k, i)
                  enddo
               enddo
            enddo
         endif
         call iamas_timer_stop("calculate2")

         call iamas_timer_start("calculate3")
         ! now deep convection
         call neg_check_ct('deep', pret, ktop, epsilc, dt, &
              tracer, tracert, iopt, num_chem, &
              its, ite, kts, kte, ipr, jpr, npr, j)
         do i = its, ite
            do k = kts, kte
                  if (pret(i) .gt. 0.) then
                     do nv = 1, num_chem
                        chem(nv, k, i) = max(epsilc, chem(nv, k, i) &
                                          + tracert(nv, k, i) * dt)
                        ! if (nv .eq. npr) then
                        !    cu_co_ten(i, k, j) = tracert(i, k, npr) * dt
                        ! endif
                        tracert_3d(nv, k, i) = tracert_3d(nv, k, i) + tracert(nv, k, i)
                     enddo
                  else
                     do nv = 1, num_chem
                        tracert(nv, k, i) = 0.
                        ! if (nv .eq. npr)cu_co_ten(i, k, j) = 0.
                     enddo
                  endif
            enddo
         enddo
         ! jwfeng: Check if tracert is conserved vertically when cu_scav_on==.false.
         if (.not. cu_scav_on) then
            do i = its, ite
               do nv = 1, num_chem
                  tracert_ratio = dt*sum(tracert_3d(nv,:,i)*dz8w(:,i)*rho_phy(:,i))/ &
                                     sum(      chem(nv,:,i)*dz8w(:,i)*rho_phy(:,i))*100.0
                  if (abs(tracert_ratio) .ge. 1.0 .and. sum(chem(nv,:,i)) .ge. 1.0e-10) then
                     write(message, "(a,I5,x,a,I3,x,a,F9.4,a)") "iCell=", i, "iChem=", nv, "diff_ratio=", tracert_ratio, "%"
                     call chemistry_message("Tendencies of ctrans are not conserved at:")
                     call chemistry_message(trim(message))
                  end if
               enddo
            enddo
         endif
         ktop_2d(its:ite) = ktop(its:ite)
         kbcon_2d(its:ite)= kbcon(its:ite)
         call iamas_timer_stop("calculate3")

   end subroutine general_cup_transport

   subroutine neg_check_ct(name, pret, ktop, epsilc, dt, q, outq, iopt, &
              num_chem, its, ite, kts, kte, ipr, jpr, npr, j)

      integer, intent(in) :: iopt, num_chem, its, ite, kts, kte
      integer, intent(in) :: ipr, jpr, npr, j

      real, dimension(1:num_chem, kts:kte, its:ite), intent(inout) :: q, outq
      real, dimension(its:ite), intent(in) :: pret
      integer, dimension(its:ite), intent(in) :: ktop
      real, intent(in) :: dt, epsilc
      real :: tracermin, tracermax, thresh, qmem, qmemf, qmem2, qtest, qmem1
      character*(*) name
      character(len=300) message
      integer nv, i, k
      ! check whether routine produces negative q's. This can happen, since
      ! tendencies are calculated based on forced q's. This should have no
      ! influence on conservation properties, it scales linear through all
      ! tendencies. Use iopt=0 to test for each tracer seperately, iopt=1
      ! for a more severe limitation...
      thresh = epsilc
      !thresh=1.e-30
      if (iopt .eq. 0) then
         do nv = 1, num_chem
            do i = its, ite
               if (pret(i) .le. 0 .and. name .eq. 'deep') cycle
               tracermin = q(nv, kts, i)
               tracermax = q(nv, kts, i)
               do k = kts + 1, kte - 1
                  tracermin = min(tracermin, q(nv, k, i))
                  tracermax = max(tracermax, q(nv, k, i))
               enddo
               tracermin = max(tracermin, thresh)
               qmemf = 1.
               ! first check for minimum restriction
               do k = kts, ktop(i)
               ! tracer tendency
                  qmem = outq(nv, k, i)
                  ! only necessary if there is a tendency
                  if (qmem .lt. 0.) then
                     qtest = q(nv, k, i) + outq(nv, k, i) * dt
                     if (qtest .lt. tracermin) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(nv, k, i)
                        qmem2 = (tracermin - q(nv, k, i)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        if (qmemf .gt. 1.) then
                           write(message, "(a,x,2F9.3)") 'something wrong in negct_1', qmem2, qmem1
                           call chemistry_message(trim(message))
                        endif
                        if (i .eq. ipr .and. nv .eq. npr) then
                           write(message, "(I3,x,4F9.3)") k, qtest, qmem2, qmem1, qmemf
                           call chemistry_message(trim(message))
                        endif
                        qmemf = max(qmemf, 0.)
                     endif
                  endif
               enddo
               do k = kts, ktop(i)
                  outq(nv, k, i) = outq(nv, k, i) * qmemf
               enddo
               ! now check max
               qmemf = 1.
               do k = kts, ktop(i)
                  ! tracer tendency
                  qmem = outq(nv, k, i)
                  ! only necessary if there is a tendency
                  if (qmem .gt. 0.) then
                     qtest = q(nv, k, i) + outq(nv, k, i) * dt
                     if (qtest .gt. tracermax) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(nv, k, i)
                        qmem2 = (tracermax - q(nv, k, i)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        if (qmemf .gt. 1.) then
                           write(message, "(a,x,2F9.3)") 'something wrong in negct_2', qmem2, qmem1
                           call chemistry_error_fatal(trim(message))
                        endif
                        if (i .eq. ipr .and. nv .eq. npr) then
                           write(message, "(I3,x,4F9.3)") k, qtest, qmem2, qmem1, qmemf
                           call chemistry_message(trim(message))
                        endif
                        qmemf = max(qmemf, 0.)
                     endif
                  endif
               enddo
               do k = kts, ktop(i)
                  outq(nv, k, i) = outq(nv, k, i) * qmemf
               enddo
            enddo
         enddo
      ! ELSE
      elseif (iopt .eq. 1) then
         do i = its, ite
            qmemf = 1.
            do k = kts, ktop(i)
               do nv = 1, num_chem
                  ! tracer tendency
                  qmem = outq(nv, k, i)
                  ! only necessary if tendency is larger than zero
                  if (qmem .lt. 0.) then
                     qtest = q(nv, k, i) + outq(nv, k, i) * dt
                     if (qtest .lt. thresh) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(nv, k, i)
                        qmem2 = (thresh - q(nv, k, i)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        qmemf = max(0., qmemf)
                     endif
                  endif
               enddo
            enddo
            do nv = 1, num_chem
               do k = kts, ktop(i)
                  outq(nv, k, i) = outq(nv, k, i) * qmemf
               enddo
            enddo
         enddo
      endif
   end subroutine neg_check_ct
!-------------------------------------------------------
END MODULE module_ctrans_general
