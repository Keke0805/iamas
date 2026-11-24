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

module module_cu_gd

   implicit none

contains

   subroutine cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,maxens2,edtc,               &
              its,ite, kts,kte                     )

      integer, intent(in) :: its, ite, kts, kte
      integer, intent(in) :: maxens2
       ! ierr error value, maybe modified in this routine
      real, dimension (kts:kte,its:ite), intent(in) :: us, vs, z, p
      real, dimension (1:maxens2,its:ite), intent(out) :: edtc
      real, dimension (its:ite), intent(out) :: edt
      real, dimension (its:ite), intent(in) :: pwav, pwev
      real, intent(in) :: edtmax, edtmin
      integer, dimension (its:ite), intent(in) :: ktop, kbcon
      integer, dimension (its:ite), intent(inout) :: ierr

      !  local variables in this routine
      integer i, k, kk
      real    einc, pef, pefb, prezk, zkbc
      real, dimension (its:ite) :: vshear, sdp, vws
      !--- determine downdraft strength in terms of windshear
      ! */ calculate an average wind shear over the depth of the cloud
      do i = its, ite
         edt(i) = 0.
         vws(i) = 0.
         sdp(i) = 0.
         vshear(i) = 0.
      enddo
      do i = its, ite
         do kk = kts, kte - 1
            if (ierr(i) .ne. 0) cycle
            if (kk .le. min0(kte,ktop(i)) .and. kk .ge. kbcon(i)) then
               vws(i) = vws(i) + &
                  (abs((us(kk+1,i)-us(kk,i))/(z(kk+1,i)-z(kk,i))) &
                  +   abs((vs(kk+1,i)-vs(kk,i))/(z(kk+1,i)-z(kk,i)))) * &
                  (p(kk,i) - p(kk+1,i))
               sdp(i) = sdp(i) + p(kk,i) - p(kk+1,i)
            endif
            if (kk .eq. kte-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
         enddo
      enddo
      do i = its, ite
         if (ierr(i) .eq. 0) then
            pef=(1.591-.639*vshear(i)+.0953*(vshear(i)**2) &
               -.00496*(vshear(i)**3))
            if (pef .gt. edtmax) pef = edtmax
            if (pef .lt. edtmin) pef = edtmin
!
!--- cloud base precip efficiency
!
            zkbc=z(kbcon(i),i)*3.281e-3
            prezk = .02
            if (zkbc .gt. 3.) then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
                  *(- 1.2569798e-2+zkbc*(4.2772e-4-zkbc*5.44e-6))))
            endif
            if (zkbc .gt. 25) then
               prezk = 2.4
            endif
            pefb = 1. / (1. + prezk)
            if (pefb .gt. edtmax) pefb = edtmax
            if (pefb .lt. edtmin) pefb = edtmin
            edt(i) = 1. - .5 * (pefb + pef)
            !--- edt here is 1-precipeff!
            !           einc=(1.-edt(i))/float(maxens2)
            !           einc=edt(i)/float(maxens2+1)
            !--- 20 percent
            einc = .2 * edt(i)
            do k = 1, maxens2
               edtc(k,i) = edt(i) + float(k - 2) * einc
            enddo
         endif
      enddo
      do i = its, ite
         if (ierr(i) .eq. 0) then
            do k = 1, maxens2
               edtc(k,i) = -edtc(k,i) * pwav(i) / pwev(i)
               if (edtc(k,i) .gt. edtmax)edtc(k,i) = edtmax
               if (edtc(k,i) .lt. edtmin)edtc(k,i) = edtmin
            enddo
         endif
      enddo

   end subroutine cup_dd_edt

   subroutine cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr, &
         jmin,ierr,he,dby,he_cup,                 &
         its,ite, kts,kte                )

      integer, intent (in) :: its, ite, kts, kte
      ! hcd = downdraft moist static energy
      ! he = moist static energy on model levels
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! dby = buoancy term
      ! cdd = detrainment function
      ! z_cup = heights of model cloud levels
      ! entr = entrainment rate
      ! zd = downdraft normalized mass flux
      !
      real, dimension (kts:kte,its:ite), intent (in) :: he, he_cup, hes_cup, z_cup, cdd, zd
      ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: jmin
      ! input and output
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (kts:kte,its:ite), intent (out) :: hcd, dby
      !  local variables in this routine
      integer :: i, k, ki
      real :: dz

      do i = its, ite
         do k =  kts + 1, kte
            dby(k,i) = 0.
            if (ierr(i) .eq. 0) then
               hcd(k,i) = hes_cup(k,i)
            endif
         enddo
      enddo

      do i = its, ite
         if (ierr(i) .eq. 0) then
            k = jmin(i)
            hcd(k,i) = hes_cup(k,i)
            dby(k,i) = hcd(jmin(i),i) - hes_cup(k,i)
!
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(ki + 1,i) - z_cup(ki,i)
               hcd(ki,i)=(hcd(ki+1,i)*(1.-.5*cdd(ki,i)*dz) &
                  +entr*dz*he(ki,i) &
                  )/(1.+entr*dz-.5*cdd(ki,i)*dz)
               dby(ki,i) = hcd(ki,i) - hes_cup(ki,i)
            enddo
         endif
      enddo
   end subroutine cup_dd_he

   subroutine cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup, &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,         &
              gamma_cup,pwev,bu,qrcd,                     &
              q,he,t_cup,iloop,xl,                        &
              its,ite, kts,kte                   )

      integer, intent (in) :: its, ite, kts, kte
      ! cdd = detrainment function
      ! q = environmental q on model levels
      ! q_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! hes_cup = saturation h on model cloud levels
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      ! qcd = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd = evaporate at that level
      ! pwev = total normalized integrated evaoprate (i2)
      ! entr = entrainment rate
      real, dimension(kts:kte,its:ite), intent(in) :: zd, t_cup, hes_cup
      real, dimension(kts:kte,its:ite), intent(in) :: hcd, qes_cup, q_cup
      real, dimension(kts:kte,its:ite), intent(in) :: z_cup, cdd, gamma_cup
      real, dimension(kts:kte,its:ite), intent(in) :: q, he
      real, intent (in) :: entr, xl
      integer, intent (in) :: iloop
      integer, dimension (its:ite), intent (in) :: jmin
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (kts:kte,its:ite), intent (out) :: qcd, qrcd, pwd
      real, dimension (its:ite), intent (out) :: pwev, bu
      !  local variables in this routine
      integer :: i, k, ki
      real :: dh, dz, dqeva

      do i = its, ite
         bu(i) = 0.
         pwev(i) = 0.
      enddo
      do i = its, ite
         do k = kts, kte
            qcd(k,i) = 0.
            qrcd(k,i) = 0.
            pwd(k,i) = 0.
         enddo
      enddo

      do i = its, ite
         if (ierr(i) .eq. 0) then
            k = jmin(i)
            dz = z_cup(k + 1,i) - z_cup(k,i)
            qcd(k,i) = q_cup(k,i)
            ! qcd(k,i) = .5 * (qes_cup(k,i) + q_cup(k,i))
            qrcd(k,i) = qes_cup(k,i)
            pwd(jmin(i),i) = min(0., qcd(k,i) - qrcd(k,i))
            pwev(i) = pwev(i) + pwd(jmin(i),i)
            qcd(k,i) = qes_cup(k,i)
            dh = hcd(k,i) - hes_cup(k, i)
            bu(i) = dz * dh
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(ki + 1,i) - z_cup(ki,i)
               qcd(ki,i)=(qcd(ki+1,i)*(1.-.5*cdd(ki,i)*dz) &
                  +entr*dz*q(ki,i) &
                  )/(1.+entr*dz-.5*cdd(ki,i)*dz)
               !--- to be negatively buoyant, hcd should be smaller than hes!
               dh = hcd(ki,i) - hes_cup(ki,i)
               bu(i) = bu(i) + dz * dh
               qrcd(ki,i)=qes_cup(ki,i)+(1./xl)*(gamma_cup(ki,i) &
                  /(1.+gamma_cup(ki,i)))*dh
               dqeva = qcd(ki,i) - qrcd(ki,i)
               if (dqeva .gt. 0.) dqeva = 0.
               pwd(ki,i) = zd(ki,i) * dqeva
               qcd(ki,i) = qrcd(ki,i)
               pwev(i) = pwev(i) + pwd(ki,i)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),hes_cup(i,ki),dh,dqeva
!        endif
            enddo
            if (pwev(i) .eq. 0 .and. iloop .eq. 1) then
               ierr(i) = 7
            endif
            if (bu(i) .ge. 0 .and. iloop .eq. 1) then
               ierr(i) = 7
            endif
         endif
      enddo

   end subroutine cup_dd_moisture

   subroutine cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr, &
              itest,kdet,z1,                          &
              its,ite, kts,kte               )

      integer, intent (in) :: its, ite, kts, kte
      ! z_cup = height of cloud model level
      ! z1 = terrain elevation
      ! entr = downdraft entrainment rate
      ! jmin = downdraft originating level
      ! kdet = level above ground where downdraft start detraining
      ! itest = flag to whether to calculate cdd

      real, dimension (kts:kte,its:ite), intent (in) :: z_cup
      real, dimension (its:ite), intent (in) :: z1
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: jmin, kdet
      integer, intent (in) :: itest
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      ! zd is the normalized downdraft mass flux
      ! cdd is the downdraft detrainmen function

      real, dimension (kts:kte,its:ite), intent (out) :: zd
      real, dimension (kts:kte,its:ite), intent (inout) :: cdd
      !  local variables in this routine
      integer :: i, k, ki
      real :: a, perc, dz
      !--- perc is the percentage of mass left when hitting the ground
      perc = .03
      do i = its, ite
         do k = kts, kte
            zd(k, i) = 0.
            if (itest .eq. 0)cdd(k, i) = 0.
         enddo
      enddo
      a = 1. - perc

      do i = its, ite
         if (ierr(i) .eq. 0) then
            zd(jmin(i),i) = 1.
            !--- integrate downward, specify detrainment(cdd)!
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(ki + 1,i) - z_cup(ki,i)
               if (ki .le. kdet(i) .and. itest .eq. 0) then
                  cdd(ki,i)=entr+(1.- (a*(z_cup(ki,i)-z1(i)) &
                     +perc*(z_cup(kdet(i),i)-z1(i)) ) &
                     /(a*(z_cup(ki+1,i)-z1(i)) &
                     +perc*(z_cup(kdet(i),i)-z1(i))))/dz
               endif
               zd(ki,i)=zd(ki+1,i)*(1.+(entr-cdd(ki,i))*dz)
            enddo
         endif
      enddo

   end subroutine cup_dd_nms

   subroutine cup_env(z,qes,he,hes,t,q,p,z1, &
         psur,ierr,tcrit,itest,xl,cp,        &
         its,ite, kts,kte           )

      integer, intent (in) :: its, ite, kts, kte
      !
      ! ierr error value, maybe modified in this routine
      ! q = environmental mixing ratio
      ! qes = environmental saturation mixing ratio
      ! t = environmental temp
      ! tv = environmental virtual temp
      ! p = environmental pressure
      ! z = environmental heights
      ! he = environmental moist static energy
      ! hes = environmental saturation moist static energy
      ! psur = surface pressure
      ! z1 = terrain elevation
      !
      real, dimension (kts:kte, its:ite), intent (in) :: p, t
      real, dimension (kts:kte, its:ite), intent (out) :: he, hes, qes
      real, dimension (kts:kte, its:ite), intent (inout) :: z, q
      real, dimension (its:ite), intent (in) :: psur, z1
      real, intent (in) :: xl, cp
      integer, dimension (its:ite), intent (inout) :: ierr
      integer, intent (in) :: itest

      !  local variables in this routine
      integer :: i, k, iph
      real, dimension (1:2) :: ae,be,ht
      real, dimension (kts:kte, its:ite) :: tv
      real :: tcrit,e,tvbar

      ht(1) = xl / cp
      ht(2) = 2.834e6 / cp
      be(1) = .622 * ht(1) / .286
      ae(1) = be(1) / 273. + alog(610.71)
      be(2) = .622 * ht(2) / .286
      ae(2) = be(2) / 273. + alog(610.71)
      do  i = its, ite
         do k = kts, kte
            if (ierr(i) .eq. 0) then
               !csgb - iph is for phase, dependent on tcrit (water or ice)
               iph = 1
               if (t(k, i) .le. tcrit) iph = 2
               e = exp(ae(iph) - be(iph) / t(k, i))
               qes(k, i) = .622 * e / (100. * p(k, i) - e)
               if (qes(k, i) .le. 1.e-08)qes(k, i) = 1.e-08
               if (q(k, i) .gt. qes(k, i))q(k, i) = qes(k, i)
               tv(k, i) = t(k, i) + .608 * q(k, i) * t(k, i)
            endif
         enddo
      enddo
      !--- z's are calculated with changed h's and q's and t's
      !--- if itest=2
      if (itest .ne. 2) then
         do i = its, ite
            if (ierr(i) .eq. 0) then
               z(1, i) = max(0., z1(i)) - (alog(p(1, i)) - &
                  alog(psur(i))) * 287. * tv(1, i) / 9.81
            endif
         enddo

         ! --- calculate heights
         do i = its , ite
            do k = kts+1, kte
               if (ierr(i) .eq. 0) then
                  tvbar = .5 * tv(k, i) + .5 * tv(k-1, i)
                  z(k, i) = z(k-1, i) - (alog(p(k, i)) - &
                     alog(p(k-1, i))) * 287. * tvbar / 9.81
               endif
            enddo
         enddo
      else
         do i = its, ite
            do k = kts, kte
               if (ierr(i) .eq. 0) then
                  z(k,i)=(he(k,i)-1004.*t(k,i)-2.5e6*q(k,i))/9.81
                  z(k,i) = max(1.e-3, z(k,i))
               endif
            enddo
         enddo
      endif
      !--- calculate moist static energy - he
      !    saturated moist static energy - hes
      do i = its, ite
         do k = kts, kte
            if (ierr(i) .eq. 0) then
               if(itest .eq. 0) he(k,i)=9.81*z(k,i)+1004.*t(k,i)+2.5e06*q(k,i)
               hes(k,i)=9.81*z(k,i)+1004.*t(k,i)+2.5e06*qes(k,i)
               if (he(k,i) .ge. hes(k,i))he(k,i) = hes(k,i)
            endif
         enddo
      enddo

   end subroutine cup_env

   subroutine cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              its,ite, kts,kte                        )

      integer, intent (in) :: its, ite, kts, kte
      !
      ! ierr error value, maybe modified in this routine
      ! q = environmental mixing ratio
      ! q_cup = environmental mixing ratio on cloud levels
      ! qes = environmental saturation mixing ratio
      ! qes_cup = environmental saturation mixing ratio on cloud levels
      ! t = environmental temp
      ! t_cup = environmental temp on cloud levels
      ! p = environmental pressure
      ! p_cup = environmental pressure on cloud levels
      ! z = environmental heights
      ! z_cup = environmental heights on cloud levels
      ! he = environmental moist static energy
      ! he_cup = environmental moist static energy on cloud levels
      ! hes = environmental saturation moist static energy
      ! hes_cup = environmental saturation moist static energy on cloud levels
      ! gamma_cup = gamma on cloud levels
      ! psur = surface pressure
      ! z1 = terrain elevation
      !
      real, dimension(kts:kte, its:ite), intent(in) :: qes, q, he
      real, dimension(kts:kte, its:ite), intent(in) :: hes, z, p, t
      real, dimension(kts:kte, its:ite), intent(out) :: qes_cup, q_cup
      real, dimension(kts:kte, its:ite), intent(out) :: he_cup, hes_cup
      real, dimension(kts:kte, its:ite), intent(out) :: z_cup, p_cup
      real, dimension(kts:kte, its:ite), intent(out) :: gamma_cup, t_cup
      real, dimension(its:ite), intent(in) :: psur, z1
      real, intent (in) :: xl, rv, cp
      integer, dimension (its:ite), intent (inout) :: ierr
      !  local variables in this routine
      integer :: i, k

      do i = its , ite
         do k = kts + 1, kte
            if (ierr(i) .eq. 0) then
               qes_cup(k,i)=.5*(qes(k-1,i)+qes(k,i))
               q_cup(k,i)=.5*(q(k-1,i)+q(k,i))
               hes_cup(k,i)=.5*(hes(k-1,i)+hes(k,i))
               he_cup(k,i)=.5*(he(k-1,i)+he(k,i))
               if (he_cup(k, i) .gt. hes_cup(k, i))he_cup(k, i) = hes_cup(k, i)
               z_cup(k,i)=.5*(z(k-1,i)+z(k,i))
               p_cup(k,i)=.5*(p(k-1,i)+p(k,i))
               t_cup(k,i)=.5*(t(k-1,i)+t(k,i))
               gamma_cup(k,i)=(xl/cp)*(xl/(rv*t_cup(k,i) &
                  *t_cup(k,i)))*qes_cup(k,i)
            endif
         enddo
      enddo
      do i=its,ite
         if(ierr(i).eq.0)then
            qes_cup(1,i)=qes(1,i)
            q_cup(1,i)=q(1,i)
            hes_cup(1,i)=hes(1,i)
            he_cup(1,i)=he(1,i)
            z_cup(1,i)=.5*(z(1,i)+z1(i))
            p_cup(1,i)=.5*(p(1,i)+psur(i))
            t_cup(1,i)=t(1,i)
            gamma_cup(1,i)=xl/cp*(xl/(rv*t_cup(1,i) &
               *t_cup(1,i)))*qes_cup(1,i)
         endif
      enddo

   end subroutine cup_env_clev

   subroutine cup_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
                        ierr, kbmax, p_cup, cap_max, &
                        its, ite, kts, kte)

      integer, intent (in) :: its, ite, kts, kte
      real, dimension (kts:kte,its:ite), intent (in) :: he_cup, hes_cup, p_cup
      real, dimension (its:ite), intent (in) :: cap_max, cap_inc
      integer, dimension (its:ite), intent (in) :: kbmax
      integer, dimension (its:ite), intent (inout) :: kbcon, k22, ierr
      integer, intent (in) :: iloop

      !local variables
      integer :: i
      real :: pbcdif, plus
      ! --- determine the level of convective cloud base  - kbcon
      do i = its, ite
         kbcon(i)=1
         if(ierr(i).ne.0)cycle
         kbcon(i)=k22(i)
         go to 32
31 continue
         kbcon(i) = kbcon(i) + 1
         if (kbcon(i) .gt. kbmax(i) + 2) then
            if (iloop .lt. 4)ierr(i) = 3
!        if (iloop .lt. 4)ierr(i) = 997
            cycle
         endif
32 continue
         if (he_cup(k22(i),i) .lt. hes_cup(kbcon(i),i)) go to 31

         ! cloud base pressure and max moist static energy pressure
         ! i.e., the depth (in mb) of the layer of negative buoyancy
         if (kbcon(i) - k22(i) .eq. 1) cycle
         pbcdif = - p_cup(kbcon(i),i) + p_cup(k22(i),i)
         plus = max(25., cap_max(i) - float(iloop - 1) * cap_inc(i))
         if (iloop .eq. 4) plus = cap_max(i)
         if (pbcdif .gt. plus) then
            k22(i) = k22(i) + 1
            kbcon(i) = k22(i)
            go to 32
         endif
      enddo

   end subroutine cup_kbcon

   subroutine cup_ktop(ilo,dby,kbcon,ktop,ierr, &
         its,ite, kts,kte         )

      integer, intent (in) :: its, ite, kts, kte
      ! dby = buoancy term
      ! ktop = cloud top (output)
      ! ilo = flag
      ! ierr error value, maybe modified in this routine
      real, dimension (kts:kte,its:ite), intent (inout) :: dby
      integer, dimension (its:ite), intent (in) :: kbcon
      integer, intent (in) :: ilo
      integer, dimension (its:ite), intent (out) :: ktop
      integer, dimension (its:ite), intent (inout) :: ierr
      !  local variables
      integer :: i, k

      do 42 i = its, ite
         ktop(i) = 1
         if (ierr(i) .eq. 0) then
            do 40 k = kbcon(i) + 1, kte - 1
               if (dby(k,i) .le. 0.) then
                  ktop(i) = k - 1
                  go to 41
               endif
            40 continue

               if (ilo .eq. 1) ierr(i) = 5
               !if (ilo .eq. 2)ierr(i) = 998
               go to 42
41 continue
            do k = ktop(i) + 1, kte
               dby(k,i) = 0.
            enddo
         endif
      42 continue

   end subroutine cup_ktop

   subroutine cup_maximi(array,ks,ke,maxx,ierr, &
              its,ite, kts,kte         )

      integer, intent (in) :: its, ite, kts, kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks, kend  check - range
      real, dimension (kts:kte,its:ite), intent (in) :: array
      integer, dimension (its:ite), intent (in) :: ierr, ke
      integer, intent (in) :: ks
      integer, dimension (its:ite), intent (out) :: maxx
      real, dimension (its:ite) :: x
      real :: xar
      integer :: i, k

      do i = its, ite
         maxx(i) = ks
         if (ierr(i) .eq. 0) then
            x(i) = array(ks, i)

            do k = ks, ke(i)
               xar = array(k, i)
               if (xar .ge. x(i)) then
                  x(i) = xar
                  maxx(i) = k
               endif
            enddo
         endif
      enddo

   end subroutine cup_maximi

   subroutine cup_minimi(array,ks,kend,kt,ierr, &
              its,ite, kts,kte         )

      integer, intent (in) :: its, ite, kts, kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks, kend  check - range
      real, dimension (kts:kte,its:ite), intent (in) :: array
      integer, dimension (its:ite), intent (in) :: ierr, ks, kend
      integer, dimension (its:ite), intent (out) :: kt
      real, dimension (its:ite) :: x
      integer :: i, k, kstop

      do i = its, ite
         kt(i) = ks(i)
         if (ierr(i) .eq. 0) then
            x(i) = array(ks(i),i)
            kstop = max(ks(i) + 1, kend(i))

            do k = ks(i) + 1, kstop
               if (array(k, i) .lt. x(i)) then
                  x(i) = array(k, i)
                  kt(i) = k
               endif
            enddo
         endif
      enddo

   end subroutine cup_minimi

   subroutine cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc, &
              kbcon,ierr,dby,he,hes_cup,                 &
              its,ite, kts,kte                  )

      integer, intent (in) :: its, ite, kts, kte
      ! hc = cloud moist static energy
      ! hkb = moist static energy at originating level
      ! he = moist static energy on model levels
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! dby = buoancy term
      ! cd = detrainment function
      ! z_cup = heights of model cloud levels
      ! entr = entrainment rate
      real, dimension (kts:kte,its:ite), intent (in) :: he, he_cup
      real, dimension (kts:kte,its:ite), intent (in) :: hes_cup, z_cup, cd
  ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: kbcon, k22
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (kts:kte,its:ite), intent (out) :: hc, dby
      real, dimension (its:ite), intent (out) :: hkb

      ! local variables
      integer :: i, k
      real :: dz

      !--- moist static energy inside cloud
      do i = its, ite
         if (ierr(i) .eq. 0.) then
            hkb(i) = he_cup(k22(i),i)
            do k = 1, k22(i)
               hc(k, i) = he_cup(k, i)
               dby(k, i) = 0.
            enddo
            do k = k22(i), kbcon(i) - 1
               hc(k, i) = hkb(i)
               dby(k, i) = 0.
            enddo
            k = kbcon(i)
            hc(k, i) = hkb(i)
            dby(kbcon(i),i) = hkb(i) - hes_cup(k, i)
         endif
      enddo
      do i = its , ite
         do k = kts + 1, kte
            if (k .gt. kbcon(i) .and. ierr(i) .eq. 0.) then
               dz = z_cup(k, i) - z_cup(k - 1,i)
               hc(k, i)=(hc(k-1,i)*(1.-.5*cd(k, i)*dz)+entr* &
                  dz*he(k-1,i))/(1.+entr*dz-.5*cd(k ,i)*dz)
               dby(k, i) = hc(k, i) - hes_cup(k, i)
            endif
         enddo
      enddo

   end subroutine cup_up_he

   subroutine cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,    &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,         &
              q,gamma_cup,zu,qes_cup,k22,qe_cup,xl,         &
              its,ite, kts,kte                     )

      integer, intent (in) :: its, ite, kts, kte
      ! cd = detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd = detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      real, dimension (kts:kte,its:ite), intent (in) :: q, zu
      real, dimension (kts:kte,its:ite), intent (in) ::gamma_cup, qe_cup
      real, dimension (kts:kte,its:ite), intent (in) ::dby, qes_cup
      real, dimension (kts:kte,its:ite), intent (in) ::z_cup, cd
      ! entr = entrainment rate
      real, intent (in) :: mentr_rate, xl
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22
      integer, dimension (its:ite), intent (inout) :: ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (i1)
      ! c0 = conversion rate (cloud to rain)
      real, dimension (kts:kte,its:ite), intent (out) :: qc, qrc
      real, dimension (kts:kte,its:ite), intent (out) :: pw, clw_all
      real, dimension (its:ite), intent (out) :: pwav

      ! local variables
      integer :: iall, i, k
      real :: dh, qrch, c0, dz, radius

      iall = 0
      c0 = .002
      !--- no precip for small clouds
      if (mentr_rate .gt. 0.) then
         radius = .2 / mentr_rate
         if (radius .lt. 900.) c0 = 0.
!         if (radius .lt. 900.) iall = 0
      endif
      do i = its, ite
         pwav(i) = 0.
      enddo
      do i = its, ite
         do k = kts, kte
            pw(k,i) = 0.
            if (ierr(i) .eq. 0) qc(k,i) = qes_cup(k,i)
            clw_all(k,i) = 0.
            qrc(k,i) = 0.
         enddo
      enddo
      do i = its, ite
         if (ierr(i) .eq. 0.) then
            do k = k22(i), kbcon(i) - 1
               qc(k,i) = qe_cup(k22(i),i)
            enddo
         endif
      enddo

      do i = its , ite
         do k = kts + 1, kte
            if (ierr(i) .ne. 0) cycle
            if (k .lt. kbcon(i)) cycle
            if (k .gt. ktop(i)) cycle
            dz = z_cup(k,i) - z_cup(k - 1,i)
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
            qc(k,i)=(qc(k-1,i)*(1.-.5*cd(k,i)*dz)+mentr_rate* &
               dz*q(k-1,i))/(1.+mentr_rate*dz-.5*cd(k,i)*dz)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
            qrch=qes_cup(k,i)+(1./xl)*(gamma_cup(k,i) &
               /(1.+gamma_cup(k,i)))*dby(k,i)
!
!------- liquid water content in cloud after rainout
!
            clw_all(k,i) = qc(k,i) - qrch
            qrc(k,i)=(qc(k,i)-qrch)/(1.+c0*dz*zu(k,i))
            if (qrc(k,i) .lt. 0.) then
               qrc(k,i) = 0.
            endif
!
!-------   3.condensation
!
            pw(k,i) = c0 * dz * qrc(k,i) * zu(k,i)
            if (iall .eq. 1) then
               qrc(k,i) = 0.
               pw(k,i) = (qc(k,i) - qrch) * zu(k,i)
               if (pw(k,i) .lt. 0.)pw(k,i) = 0.
            endif
!
!----- set next level
!
            qc(k,i) = qrc(k,i) + qrch
!
!--- integrated normalized ondensate
!
            pwav(i) = pwav(i) + pw(k,i)
         enddo
      enddo

   end subroutine cup_up_moisture

   subroutine cup_up_nms(zu, z_cup, entr, cd, kbcon, ktop, ierr, k22, &
              its, ite, kts, kte)

      IMPLICIT NONE
      ! on input
      integer, intent (in) :: its, ite, kts, kte
      ! cd = detrainment function
      real, dimension (kts:kte,its:ite), intent (in) :: z_cup, cd
      ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22
      ! input and output
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      ! zu is the normalized mass flux
      real, dimension (kts:kte,its:ite), intent (out) :: zu
      ! local variables in this routine

      integer :: i, k
      real :: dz
      ! initialize for this go around
      do i = its, ite
          do k = kts, kte
              zu(k, i) = 0.
          enddo
      enddo
      ! do normalized mass budget
      do i = its, ite
          if (ierr(i) .eq. 0) then
              do k = k22(i), kbcon(i)
                  zu(k, i) = 1.
              end do
              do k = kbcon(i) + 1, ktop(i)
                  dz = z_cup(k, i) - z_cup(k-1, i)
                  zu(k, i) = zu(k - 1, i) * (1. + (entr - cd(k ,i)) * dz)
              end do
          end if
      end do

   end subroutine cup_up_nms

!-------------------------------------------------------
end module module_cu_gd
