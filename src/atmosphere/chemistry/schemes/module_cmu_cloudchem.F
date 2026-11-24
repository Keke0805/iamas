
!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!**********************************************************************************  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!      xyruan modify this code for iamas. 06/2026          !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_cmu_cloudchem
 use iamas_kind_types, only: RKIND
 use iamas_atmchem_vars, only: num_chem, param_first_scalar, p_qc,        &
      p_nh3, p_hno3, p_hcl, p_sulf, p_h2so4, p_hcho, p_ora1,             &
      p_so2, p_h2o2, p_o3, p_ho, p_ho2, p_no3, p_no, p_no2,              &
      p_hono, p_pan, p_op1
 use module_data_binaer_other, only: l2maxd

 implicit none
 private

 integer, parameter :: l_so4_aqyy  = 1
 integer, parameter :: l_no3_aqyy  = 2
 integer, parameter :: l_cl_aqyy   = 3
 integer, parameter :: l_nh4_aqyy  = 4
 integer, parameter :: l_na_aqyy   = 5
 integer, parameter :: l_oin_aqyy  = 6
 integer, parameter :: l_bc_aqyy   = 7
 integer, parameter :: l_oc_aqyy   = 8
 integer, parameter :: l_dust_aqyy = 9
 integer, parameter :: l_tr1r1_aqyy = 10
 integer, parameter :: l_tr1r2_aqyy = 11
 integer, parameter :: l_tr1r3_aqyy = 12
 integer, parameter :: l_tr1r4_aqyy = 13
 integer, parameter :: nyyy = 13

 public :: cmu_cloudchem_driver

contains

!=======================================================================
 subroutine cmu_cloudchem_driver( id, ktau, ktauc, dtstepc,          &
                                     num_moist,                         &
                                     t_phy, rho_phy, p_phy,             &
                                     moist, chem, gas_aqfrac,           &
                                     numgas_aqfrac, cldfra, ph_no2,     &
                                     its, ite, kts, kte )
!=======================================================================
 ! NOTE: num_moist is passed in as an argument rather than use-associated
 ! from iamas_atmphys_vars. That module declares num_moist as a
 ! BIND(C, name="num_moist") variable; pulling it into this module makes the
 ! C external symbol "num_moist" travel with this .mod and collides with the
 ! same symbol re-exported through other modules in any unit that USEs both.
 use iamas_atmchem_utilities, only: chemistry_debug_message, chemistry_error_fatal
 use module_data_binaer_asect, only: cw_phase, nphase_aer
 use module_data_binaer_other, only: k_pegbegin
 use module_mosaic_driver, only: mapaer_tofrom_host

 implicit none

 integer, intent(in) :: id, ktau, ktauc
 integer, intent(in) :: num_moist, numgas_aqfrac
 integer, intent(in) :: its, ite, kts, kte
 real(kind=RKIND), intent(in) :: dtstepc
 real(kind=RKIND), dimension(kts:kte, its:ite), intent(in) :: t_phy, rho_phy, p_phy
 real(kind=RKIND), dimension(kts:kte, its:ite), intent(in) :: cldfra, ph_no2
 real(kind=RKIND), dimension(1:num_moist, kts:kte, its:ite), intent(in) :: moist
 real(kind=RKIND), dimension(1:num_chem,  kts:kte, its:ite), intent(inout) :: chem
 real(kind=RKIND), dimension(1:numgas_aqfrac, kts:kte, its:ite), intent(inout) :: gas_aqfrac
!   gas_aqfrac - fraction (0-1) of gas that is dissolved in cloud water
 
 integer :: i, k
 integer :: icase, kpeg, k_pegshift, mpeg
 integer :: iphotol_onoff, iradical_onoff
 real(kind=RKIND) :: gas_aqfrac_box(numgas_aqfrac)
 real :: ph_aq_box
 real(kind=RKIND), parameter :: qcldwtr_cutoff = 1.0e-6_RKIND
 character(len=256) :: msg

 if ((cw_phase <= 0) .or. (cw_phase > nphase_aer)) then
    call chemistry_error_fatal( &
         'Error: CMU_cloudchem requires active cloud-water aerosol phase.')
 endif

 if ((p_qc < 1) .or. (p_qc > num_moist)) then
    call chemistry_error_fatal( &
         'Error: CMU_cloudchem cannot find valid cloud-water moisture index.')
 endif

 iphotol_onoff = 1
!   following line turns aqueous radical chem off unconditionally
 iradical_onoff = 0

 icase = 0
 k_pegshift = k_pegbegin - kts
 mpeg = 1
 gas_aqfrac_box(:) = 0.0_RKIND
 gas_aqfrac(1:numgas_aqfrac,kts:kte,its:ite) = 0.0_RKIND
 ph_aq_box = 0.0

 do i = its, ite
    do k = kts, kte
       if (moist(p_qc,k,i) > qcldwtr_cutoff) then
          icase = icase + 1
          kpeg = k + k_pegshift

          call mapaer_tofrom_host( imap = 0,                         &
               its = its, ite = ite, kts = kts, kte = kte,            &
               it = i, ktmaps = k, ktmape = k,                        &
               moist = moist, chem = chem,                            &
               t_phy = t_phy, p_phy = p_phy, rho_phy = rho_phy        &
               )


          call cmu_cloudchem_1box(                                  &
               id, ktau, ktauc, dtstepc,                               &
               iphotol_onoff, iradical_onoff, ph_no2(k,i),             &
               ph_aq_box, gas_aqfrac_box, numgas_aqfrac,               &
               i, 1, k, kpeg, mpeg, icase )


          call mapaer_tofrom_host( imap = 1,                         &
               its = its, ite = ite, kts = kts, kte = kte,            &
               it = i, ktmaps = k, ktmape = k,                        &
               moist = moist, chem = chem,                            &
               t_phy = t_phy, p_phy = p_phy, rho_phy = rho_phy        &
               )

          gas_aqfrac(1:numgas_aqfrac,k,i) = gas_aqfrac_box(1:numgas_aqfrac)

       endif
    enddo
 enddo

 write(msg, '(A,I0)') 'leaving cmu_cloudchem_driver - ktau =', ktau
 call chemistry_debug_message(trim(msg))

 end subroutine cmu_cloudchem_driver

!=======================================================================
 subroutine cmu_cloudchem_1box( id, ktau, ktauc, dtstepc,            &
                                   iphotol_onoff, iradical_onoff,       &
                                   photol_no2_box,                      &
                                   ph_aq_box, gas_aqfrac_box,           &
                                   numgas_aqfrac,                       &
                                   it, jt, kt, kpeg, mpeg, icase )
!=======================================================================
 use module_data_binaer_asect, only: msectional, maxd_asize, maxd_atype,&
      cw_phase, lptr_so4_aer, lptr_no3_aer, lptr_cl_aer,               &
      lptr_nh4_aer, lptr_na_aer, lptr_oin_aer, lptr_bc_aer,            &
      lptr_oc_aer, lptr_dust_aer, lptr_tr1r1_aer, lptr_tr1r2_aer,      &
      lptr_tr1r3_aer, lptr_tr1r4_aer
 use module_data_binaer_other, only: ltot2, rsub
 use module_data_cmu_bulkaqchem, only: meqn1max

 implicit none

 integer, intent(in) :: id, ktau, ktauc
 integer, intent(in) :: iphotol_onoff, iradical_onoff
 integer, intent(in) :: numgas_aqfrac, it, jt, kt, kpeg, mpeg, icase
 real(kind=RKIND), intent(in) :: dtstepc, photol_no2_box
 real, intent(inout) :: ph_aq_box
 real(kind=RKIND), intent(inout), dimension(1:numgas_aqfrac) :: gas_aqfrac_box

 integer :: icase_in, idecomp_hmsa_hso5, iphase, iradical_in, istat_aqop
 integer :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)
 real :: co2_mixrat_in, ph_cmuaq_cur, photol_no2_in, xprescribe_ph
 real :: yaq_beg(meqn1max), yaq_end(meqn1max)
 real :: rbox(l2maxd), rbox_sv1(l2maxd)
 real :: rbulk_cwaer(nyyy,2)
 real, dimension(maxd_asize,maxd_atype) :: fr_partit_cw

 ph_cmuaq_cur = 0.0
 yaq_beg(:) = 0.0
 yaq_end(:) = 0.0

 iphase = cw_phase
 lptr_yyy_cwaer(:,:,l_so4_aqyy)  = lptr_so4_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_no3_aqyy)  = lptr_no3_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_cl_aqyy )  = lptr_cl_aer( :,:,iphase)
 lptr_yyy_cwaer(:,:,l_nh4_aqyy)  = lptr_nh4_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_na_aqyy )  = lptr_na_aer( :,:,iphase)
 lptr_yyy_cwaer(:,:,l_oin_aqyy)  = lptr_oin_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_bc_aqyy )  = lptr_bc_aer( :,:,iphase)
 lptr_yyy_cwaer(:,:,l_oc_aqyy )  = lptr_oc_aer( :,:,iphase)
 lptr_yyy_cwaer(:,:,l_dust_aqyy) = lptr_dust_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_tr1r1_aqyy) = lptr_tr1r1_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_tr1r2_aqyy) = lptr_tr1r2_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_tr1r3_aqyy) = lptr_tr1r3_aer(:,:,iphase)
 lptr_yyy_cwaer(:,:,l_tr1r4_aqyy) = lptr_tr1r4_aer(:,:,iphase)

!
!   xfer from rsub to rbox
!
 rbox(:) = 0.0
 rbox_sv1(:) = 0.0
 if ((ltot2 > 0) .and. (ltot2 <= l2maxd)) then
    rbox(1:ltot2) = max(0.0, rsub(1:ltot2,kpeg,mpeg))
    rbox_sv1(1:ltot2) = rbox(1:ltot2)
 endif

!
!
!   do bulk cloud-water chemistry
!
!
 icase_in = icase
 iradical_in = 1
 idecomp_hmsa_hso5 = 1
 co2_mixrat_in = 350.0
 photol_no2_in = photol_no2_box

!   when xprescribe_ph >= 0, ph is forced to its value
 xprescribe_ph = -1.0e31

!   turn off aqueous phase photolytic and radical chemistry
!   if either of the iphotol_onoff and iradical_onoff flags are 0
 if ((iphotol_onoff <= 0) .or. (iradical_onoff <= 0)) then
    photol_no2_in = 0.0
    iradical_in = 0
 endif

 gas_aqfrac_box(:) = 0.0_RKIND

 call interface_to_aqoperator1(                                      &
      istat_aqop, dtstepc, rbox, gas_aqfrac_box,                     &
      rbulk_cwaer, lptr_yyy_cwaer,                                   &
      co2_mixrat_in, photol_no2_in, xprescribe_ph,                   &
      iradical_in, idecomp_hmsa_hso5,                                &
      yaq_beg, yaq_end, ph_cmuaq_cur,                                &
      numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase_in )

 ph_aq_box = ph_cmuaq_cur

!
!   calculate fraction of cloud-water associated with each activated aerosol bin
!
 call partition_cldwtr(rbox, fr_partit_cw, it, jt, kt, kpeg, mpeg, icase_in)

!
!   distribute changes in bulk cloud-water composition among size bins
!
 call distribute_bulk_changes(                                      &
      rbox, rbox_sv1, fr_partit_cw, rbulk_cwaer, lptr_yyy_cwaer,     &
      it, jt, kt, kpeg, mpeg, icase_in )

!
!   xfer back to rsub
!
 if ((ltot2 > 0) .and. (ltot2 <= l2maxd)) then
    rsub(1:ltot2,kpeg,mpeg) = max(0.0, rbox(1:ltot2))
 endif

!
!   do move-sections
!
 if (msectional < 1000000000) then
    call cloudchem_apply_move_sections(rbox, rbox_sv1, it, jt, kt, kpeg, mpeg, icase_in)
 endif

 end subroutine cmu_cloudchem_1box

!=======================================================================
 subroutine interface_to_aqoperator1(                                  &
      istat_aqop, dtstepc, rbox, gas_aqfrac_box,                       &
      rbulk_cwaer, lptr_yyy_cwaer,                                     &
      co2_mixrat_in, photol_no2_in, xprescribe_ph,                     &
      iradical_in, idecomp_hmsa_hso5,                                  &
      yaq_beg, yaq_end, ph_cmuaq_cur,                                  &
      numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase )
!=======================================================================
 use module_data_cmu_bulkaqchem, only: meqn1max, naers, ngas,           &
      na4, naa, nac, nae, nahmsa, nahso5, nan, nao, nar, nas, nad,     &
      nt1r1, nt1r2, nt1r3, nt1r4,                                      &
      ng4, nga, ngc, ngch3o2h, ngh2o2, nghcho, nghcooh, nghno2,        &
      ngho2, ngn, ngno, ngno2, ngno3, ngo3, ngoh, ngpan, ngso2
 use module_cmu_bulkaqchem, only: aqoperator1
 use module_data_binaer_asect, only: maxd_asize, maxd_atype, cw_phase, &
      nsize_aer, ntype_aer, mw_cl_aer, mw_na_aer, mw_nh4_aer,          &
      mw_no3_aer, mw_so4_aer
 use module_data_binaer_other, only: aboxtest_units_convert, cairclm,  &
      ktemp, ptotclm, rcldwtr_sub

 implicit none

 integer, intent(in) :: iradical_in, idecomp_hmsa_hso5
 integer, intent(in) :: numgas_aqfrac, id, it, jt, kt, kpeg, mpeg, ktau, icase
 integer, intent(inout) :: istat_aqop
 integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)
 real(kind=RKIND), intent(in) :: dtstepc
 real, intent(in) :: co2_mixrat_in, photol_no2_in, xprescribe_ph
 real, intent(inout) :: ph_cmuaq_cur
 real, intent(inout), dimension(1:l2maxd) :: rbox
 real, intent(inout), dimension(nyyy,2) :: rbulk_cwaer
 real(kind=RKIND), intent(inout), dimension(1:numgas_aqfrac) :: gas_aqfrac_box
 real, intent(inout), dimension(meqn1max) :: yaq_beg, yaq_end

 integer :: i, isize, itype, l, lyyy
 integer :: iaq, istat_fatal, istat_warn
 real, parameter :: eps = 0.622
 real :: cair_moleperm3
 real :: dum, dumb
 real :: factgas, factlwc, factpatm, factphoto
 real :: factaerbc, factaercl, factaerna, factaernh4
 real :: factaerno3, factaeroc, factaeroin, factaerdust, factaerso4
 real :: factaertr1r1, factaertr1r2, factaertr1r3, factaertr1r4
 real :: lwc, p_atm, photo_in, rh, temp
 real :: tstep_beg_sec, tstep_end_sec
 real :: totsulf_beg, totsulf_end
 real(kind=RKIND) :: gas_aqfrac_cmu(ngas)
 real :: gas(ngas), aerosol(naers)
 double precision :: tstep_beg_sec_dp, tstep_end_sec_dp
 double precision :: temp_dp, p_atm_dp, lwc_dp, rh_dp
 double precision :: co2_mixrat_in_dp, photo_in_dp, ph_cmuaq_cur_dp
 double precision :: xprescribe_ph_dp
 double precision :: gas_dp(ngas), gas_aqfrac_cmu_dp(ngas)
 double precision :: aerosol_dp(naers), yaq_beg_dp(meqn1max), yaq_end_dp(meqn1max)

 gas_aqfrac_cmu(:) = 0.0_RKIND

!
!   units conversion factors 
!   'cmuaq-bulk' value = pegasus value X factor
!
!   [pres in atmospheres] = [pres in dynes/cm2] * factpatm
 factpatm = 1.0/1.01325e6
!   [cldwtr in g-h2o/m3-air] = [cldwtr in mole-h2o/mole-air] * factlwc
 factlwc = 28.966*eps*1.0e6*cairclm(kpeg)
!   [aq photolysis rate scaling factor in --] = [jno2 in 1/min] * factphoto
 factphoto = 1.6

!   [gas in ppm] = [gas in mole/mole-air] * factgas
 factgas = 1.0e6

!   [aerosol in ug/m3-air] = [aerosol in mole/mole-air] * factaer
 dum = cairclm(kpeg)*1.0e12
 factaerso4 = dum*mw_so4_aer
 factaerno3 = dum*mw_no3_aer
 factaercl  = dum*mw_cl_aer
 factaernh4 = dum*mw_nh4_aer
 factaerna  = dum*mw_na_aer
 factaeroin = dum
 factaerdust = dum
 factaertr1r1 = dum
 factaertr1r2 = dum
 factaertr1r3 = dum
 factaertr1r4 = dum
 factaeroc = dum
 factaerbc = dum

!   If aboxtest_units_convert=10, turn off units conversions both here
!   and in module_mosaic.  This is for testing, to allow exact agreements.
 if (aboxtest_units_convert == 10) then
    factpatm = 1.0
    factlwc = 1.0
    factphoto = 1.0
    factgas = 1.0
    factaerso4 = 1.0
    factaerno3 = 1.0
    factaercl = 1.0
    factaernh4 = 1.0
    factaerna = 1.0
    factaeroin = 1.0
    factaerdust = 1.0
    factaertr1r1 = 1.0
    factaertr1r2 = 1.0
    factaertr1r3 = 1.0
    factaertr1r4 = 1.0
    factaeroc = 1.0
    factaerbc = 1.0
 endif

!
!   map from rbox to gas,aerosol
!
 temp = rbox(ktemp)

 lwc = rcldwtr_sub(kpeg,mpeg) * factlwc
 p_atm = ptotclm(kpeg) * factpatm

!   set p_atm so that cmu code's cair will match cairclm
 p_atm = cairclm(kpeg)*1.0e3*0.082058*temp
 photo_in = photol_no2_in * factphoto
 rh = 1.0
 iaq = 1
 tstep_beg_sec = 0.0
 tstep_end_sec = real(dtstepc)

!   map gases and convert to ppm
 gas(:) = 0.0
 if (p_nh3 >= param_first_scalar .and. p_nh3 <= l2maxd) &
    gas(nga) = rbox(p_nh3) * factgas
 if (p_hno3 >= param_first_scalar .and. p_hno3 <= l2maxd) &
    gas(ngn) = rbox(p_hno3) * factgas
 if (p_hcl >= param_first_scalar .and. p_hcl <= l2maxd) &
    gas(ngc) = rbox(p_hcl) * factgas
 if (p_sulf >= param_first_scalar .and. p_sulf <= l2maxd) &
    gas(ng4) = rbox(p_sulf) * factgas
 if (p_h2so4 >= param_first_scalar .and. p_h2so4 <= l2maxd) &
    gas(ng4) = rbox(p_h2so4) * factgas
 if (p_hcho >= param_first_scalar .and. p_hcho <= l2maxd) &
    gas(nghcho) = rbox(p_hcho) * factgas
 if (p_ora1 >= param_first_scalar .and. p_ora1 <= l2maxd) &
    gas(nghcooh) = rbox(p_ora1) * factgas
 if (p_so2 >= param_first_scalar .and. p_so2 <= l2maxd) &
    gas(ngso2) = rbox(p_so2) * factgas
 if (p_h2o2 >= param_first_scalar .and. p_h2o2 <= l2maxd) &
    gas(ngh2o2) = rbox(p_h2o2) * factgas
 if (p_o3 >= param_first_scalar .and. p_o3 <= l2maxd) &
    gas(ngo3) = rbox(p_o3) * factgas
 if (p_ho >= param_first_scalar .and. p_ho <= l2maxd) &
    gas(ngoh) = rbox(p_ho) * factgas
 if (p_ho2 >= param_first_scalar .and. p_ho2 <= l2maxd) &
    gas(ngho2) = rbox(p_ho2) * factgas
 if (p_no3 >= param_first_scalar .and. p_no3 <= l2maxd) &
    gas(ngno3) = rbox(p_no3) * factgas
 if (p_no >= param_first_scalar .and. p_no <= l2maxd) &
    gas(ngno) = rbox(p_no) * factgas
 if (p_no2 >= param_first_scalar .and. p_no2 <= l2maxd) &
    gas(ngno2) = rbox(p_no2) * factgas
 if (p_hono >= param_first_scalar .and. p_hono <= l2maxd) &
    gas(nghno2) = rbox(p_hono) * factgas
 if (p_pan >= param_first_scalar .and. p_pan <= l2maxd) &
    gas(ngpan) = rbox(p_pan) * factgas
 if (p_op1 >= param_first_scalar .and. p_op1 <= l2maxd) &
    gas(ngch3o2h) = rbox(p_op1) * factgas

!   compute bulk activated-aerosol mixing ratios
 aerosol(:) = 0.0
 rbulk_cwaer(:,:) = 0.0

 do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)
       do lyyy = 1, nyyy
          l = lptr_yyy_cwaer(isize,itype,lyyy)
          if (l >= param_first_scalar .and. l <= l2maxd) then
             rbulk_cwaer(lyyy,1) = rbulk_cwaer(lyyy,1) + rbox(l)
          endif
       enddo
    enddo
 enddo

!   map them to 'aerosol' array and convert to ug/m3
 aerosol(na4)   = rbulk_cwaer(l_so4_aqyy,1)  * factaerso4
 aerosol(nan)   = rbulk_cwaer(l_no3_aqyy,1)  * factaerno3
 aerosol(nac)   = rbulk_cwaer(l_cl_aqyy,1)   * factaercl
 aerosol(naa)   = rbulk_cwaer(l_nh4_aqyy,1)  * factaernh4
 aerosol(nas)   = rbulk_cwaer(l_na_aqyy,1)   * factaerna
 aerosol(nar)   = rbulk_cwaer(l_oin_aqyy,1)  * factaeroin
 aerosol(nad)   = rbulk_cwaer(l_dust_aqyy,1) * factaerdust
 aerosol(nt1r1) = rbulk_cwaer(l_tr1r1_aqyy,1) * factaertr1r1
 aerosol(nt1r2) = rbulk_cwaer(l_tr1r2_aqyy,1) * factaertr1r2
 aerosol(nt1r3) = rbulk_cwaer(l_tr1r3_aqyy,1) * factaertr1r3
 aerosol(nt1r4) = rbulk_cwaer(l_tr1r4_aqyy,1) * factaertr1r4
 aerosol(nae)   = rbulk_cwaer(l_bc_aqyy,1)   * factaerbc
 aerosol(nao)   = rbulk_cwaer(l_oc_aqyy,1)   * factaeroc

!   convert arguments to double prec
 tstep_beg_sec_dp = dble(tstep_beg_sec)
 tstep_end_sec_dp = dble(tstep_end_sec)
 temp_dp = dble(temp)
 p_atm_dp = dble(p_atm)
 lwc_dp = dble(lwc)
 rh_dp = dble(rh)
 co2_mixrat_in_dp = dble(co2_mixrat_in)
 photo_in_dp = dble(photo_in)
 xprescribe_ph_dp = dble(xprescribe_ph)
 ph_cmuaq_cur_dp = dble(ph_cmuaq_cur)

 gas_dp(:) = dble(gas(:))
 aerosol_dp(:) = dble(aerosol(:))
 gas_aqfrac_cmu_dp(:) = dble(gas_aqfrac_cmu(:))
 yaq_beg_dp(:) = dble(yaq_beg(:))
 yaq_end_dp(:) = dble(yaq_end(:))

!   total sulfur species conc as sulfate (ug/m3)
 if (temp_dp > 0.0d0) then
    cair_moleperm3 = 1.0e3*p_atm_dp/(0.082058d0*temp_dp)
 else
    cair_moleperm3 = 0.0
 endif
 totsulf_beg = ( aerosol_dp(na4)/96.0d0                                 &
              + aerosol_dp(nahso5)/113.0d0 + aerosol_dp(nahmsa)/111.0d0 &
              + (gas_dp(ngso2) + gas_dp(ng4))*cair_moleperm3 )*96.0d0

 call aqoperator1(                                                       &
      istat_fatal, istat_warn,                                           &
      tstep_beg_sec_dp, tstep_end_sec_dp,                                &
      gas_dp, aerosol_dp, gas_aqfrac_cmu_dp,                             &
      temp_dp, p_atm_dp, lwc_dp, rh_dp,                                  &
      co2_mixrat_in_dp, photo_in_dp, xprescribe_ph_dp,                   &
      iradical_in, idecomp_hmsa_hso5, iaq,                               &
      yaq_beg_dp, yaq_end_dp, ph_cmuaq_cur_dp )

 totsulf_end = ( aerosol_dp(na4)/96.0d0                                  &
              + aerosol_dp(nahso5)/113.0d0 + aerosol_dp(nahmsa)/111.0d0 &
              + (gas_dp(ngso2) + gas_dp(ng4))*cair_moleperm3 )*96.0d0

!   convert arguments back to single prec
 temp = real(temp_dp)
 p_atm = real(p_atm_dp)
 lwc = real(lwc_dp)
 rh = real(rh_dp)
 ph_cmuaq_cur = real(ph_cmuaq_cur_dp)

 gas(:) = real(gas_dp(:))
 aerosol(:) = real(aerosol_dp(:))
 gas_aqfrac_cmu(:) = real(gas_aqfrac_cmu_dp(:), RKIND)
 yaq_beg(:) = real(yaq_beg_dp(:))
 yaq_end(:) = real(yaq_end_dp(:))

!
!   warning message when status flags are non-zero
!
 istat_aqop = 0
 if (istat_fatal /= 0) then
    write(6,*) '*** cmu_cloudchem_driver, interface_to_aqoperator1'
    write(6,'(a,3i5,2i10)') '    id,it,kt, istat_fatal, warn =', &
         id, it, kt, istat_fatal, istat_warn
    istat_aqop = -10
 endif

!
!   warning message when sulfur mass balance error exceeds the greater
!   of (1.0e-3 ug/m3) OR (1.0e-3 X total sulfur mixing ratio)
!

 dum = totsulf_end - totsulf_beg
 dumb = max(totsulf_beg, totsulf_end)
 if (abs(dum) > max(1.0e-3, 1.0e-3*dumb)) then
    write(6,*) '*** cmu_cloudchem_driver, sulfur balance warning'
    write(6,'(a,3i5,1p,3e12.4)') &
         '    id,it,kt, total_sulfur_beg, _end, _error =', &
         id, it, kt, totsulf_beg, totsulf_end, dum
 endif


!
!   map from gas,aerosol to rbox
!
 if (p_nh3 >= param_first_scalar .and. p_nh3 <= l2maxd) &
    rbox(p_nh3) = gas(nga) / factgas
 if (p_hno3 >= param_first_scalar .and. p_hno3 <= l2maxd) &
    rbox(p_hno3) = gas(ngn) / factgas
 if (p_hcl >= param_first_scalar .and. p_hcl <= l2maxd) &
    rbox(p_hcl) = gas(ngc) / factgas
 if (p_sulf >= param_first_scalar .and. p_sulf <= l2maxd) &
    rbox(p_sulf) = gas(ng4) / factgas
 if (p_h2so4 >= param_first_scalar .and. p_h2so4 <= l2maxd) &
    rbox(p_h2so4) = gas(ng4) / factgas
 if (p_hcho >= param_first_scalar .and. p_hcho <= l2maxd) &
    rbox(p_hcho) = gas(nghcho) / factgas
 if (p_ora1 >= param_first_scalar .and. p_ora1 <= l2maxd) &
    rbox(p_ora1) = gas(nghcooh) / factgas
 if (p_so2 >= param_first_scalar .and. p_so2 <= l2maxd) &
    rbox(p_so2) = gas(ngso2) / factgas
 if (p_h2o2 >= param_first_scalar .and. p_h2o2 <= l2maxd) &
    rbox(p_h2o2) = gas(ngh2o2) / factgas
 if (p_o3 >= param_first_scalar .and. p_o3 <= l2maxd) &
    rbox(p_o3) = gas(ngo3) / factgas
 if (p_ho >= param_first_scalar .and. p_ho <= l2maxd) &
    rbox(p_ho) = gas(ngoh) / factgas
 if (p_ho2 >= param_first_scalar .and. p_ho2 <= l2maxd) &
    rbox(p_ho2) = gas(ngho2) / factgas
 if (p_no3 >= param_first_scalar .and. p_no3 <= l2maxd) &
    rbox(p_no3) = gas(ngno3) / factgas
 if (p_no >= param_first_scalar .and. p_no <= l2maxd) &
    rbox(p_no) = gas(ngno) / factgas
 if (p_no2 >= param_first_scalar .and. p_no2 <= l2maxd) &
    rbox(p_no2) = gas(ngno2) / factgas
 if (p_hono >= param_first_scalar .and. p_hono <= l2maxd) &
    rbox(p_hono) = gas(nghno2) / factgas
 if (p_pan >= param_first_scalar .and. p_pan <= l2maxd) &
    rbox(p_pan) = gas(ngpan) / factgas
 if (p_op1 >= param_first_scalar .and. p_op1 <= l2maxd) &
    rbox(p_op1) = gas(ngch3o2h) / factgas

 gas_aqfrac_box(:) = 0.0_RKIND

 if (p_nh3 >= param_first_scalar .and. p_nh3 <= numgas_aqfrac) &
    gas_aqfrac_box(p_nh3) = gas_aqfrac_cmu(nga)
 if (p_hno3 >= param_first_scalar .and. p_hno3 <= numgas_aqfrac) &
    gas_aqfrac_box(p_hno3) = gas_aqfrac_cmu(ngn)
 if (p_hcl >= param_first_scalar .and. p_hcl <= numgas_aqfrac) &
    gas_aqfrac_box(p_hcl) = gas_aqfrac_cmu(ngc)
 if (p_sulf >= param_first_scalar .and. p_sulf <= numgas_aqfrac) &
    gas_aqfrac_box(p_sulf) = gas_aqfrac_cmu(ng4)
 if (p_h2so4 >= param_first_scalar .and. p_h2so4 <= numgas_aqfrac) &
    gas_aqfrac_box(p_h2so4) = gas_aqfrac_cmu(ng4)
 if (p_hcho >= param_first_scalar .and. p_hcho <= numgas_aqfrac) &
    gas_aqfrac_box(p_hcho) = gas_aqfrac_cmu(nghcho)
 if (p_ora1 >= param_first_scalar .and. p_ora1 <= numgas_aqfrac) &
    gas_aqfrac_box(p_ora1) = gas_aqfrac_cmu(nghcooh)
 if (p_so2 >= param_first_scalar .and. p_so2 <= numgas_aqfrac) &
    gas_aqfrac_box(p_so2) = gas_aqfrac_cmu(ngso2)
 if (p_h2o2 >= param_first_scalar .and. p_h2o2 <= numgas_aqfrac) &
    gas_aqfrac_box(p_h2o2) = gas_aqfrac_cmu(ngh2o2)
 if (p_o3 >= param_first_scalar .and. p_o3 <= numgas_aqfrac) &
    gas_aqfrac_box(p_o3) = gas_aqfrac_cmu(ngo3)
 if (p_ho >= param_first_scalar .and. p_ho <= numgas_aqfrac) &
    gas_aqfrac_box(p_ho) = gas_aqfrac_cmu(ngoh)
 if (p_ho2 >= param_first_scalar .and. p_ho2 <= numgas_aqfrac) &
    gas_aqfrac_box(p_ho2) = gas_aqfrac_cmu(ngho2)
 if (p_no3 >= param_first_scalar .and. p_no3 <= numgas_aqfrac) &
    gas_aqfrac_box(p_no3) = gas_aqfrac_cmu(ngno3)
 if (p_no >= param_first_scalar .and. p_no <= numgas_aqfrac) &
    gas_aqfrac_box(p_no) = gas_aqfrac_cmu(ngno)
 if (p_no2 >= param_first_scalar .and. p_no2 <= numgas_aqfrac) &
    gas_aqfrac_box(p_no2) = gas_aqfrac_cmu(ngno2)
 if (p_hono >= param_first_scalar .and. p_hono <= numgas_aqfrac) &
    gas_aqfrac_box(p_hono) = gas_aqfrac_cmu(nghno2)
 if (p_pan >= param_first_scalar .and. p_pan <= numgas_aqfrac) &
    gas_aqfrac_box(p_pan) = gas_aqfrac_cmu(ngpan)
 if (p_op1 >= param_first_scalar .and. p_op1 <= numgas_aqfrac) &
    gas_aqfrac_box(p_op1) = gas_aqfrac_cmu(ngch3o2h)

 rbulk_cwaer(l_so4_aqyy,2)  = aerosol(na4)   / factaerso4
 rbulk_cwaer(l_no3_aqyy,2)  = aerosol(nan)   / factaerno3
 rbulk_cwaer(l_cl_aqyy,2)   = aerosol(nac)   / factaercl
 rbulk_cwaer(l_nh4_aqyy,2)  = aerosol(naa)   / factaernh4
 rbulk_cwaer(l_na_aqyy,2)   = aerosol(nas)   / factaerna
 rbulk_cwaer(l_oin_aqyy,2)  = aerosol(nar)   / factaeroin
 rbulk_cwaer(l_dust_aqyy,2) = aerosol(nad)   / factaerdust
 rbulk_cwaer(l_tr1r1_aqyy,2) = aerosol(nt1r1) / factaertr1r1
 rbulk_cwaer(l_tr1r2_aqyy,2) = aerosol(nt1r2) / factaertr1r2
 rbulk_cwaer(l_tr1r3_aqyy,2) = aerosol(nt1r3) / factaertr1r3
 rbulk_cwaer(l_tr1r4_aqyy,2) = aerosol(nt1r4) / factaertr1r4
 rbulk_cwaer(l_bc_aqyy,2)   = aerosol(nae)   / factaerbc
 rbulk_cwaer(l_oc_aqyy,2)   = aerosol(nao)   / factaeroc

 end subroutine interface_to_aqoperator1

!=======================================================================
 subroutine partition_cldwtr(rbox, fr_partit_cw, it, jt, kt, kpeg, mpeg, icase)
!=======================================================================
 use module_data_binaer_asect, only: maxd_asize, maxd_atype, cw_phase,  &
      nsize_aer, ntype_aer, ncomp_aer, massptr_aer, numptr_aer,         &
      dens_aer, mw_aer, volumlo_sect, volumhi_sect

 implicit none

 integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase
 real, intent(inout), dimension(1:l2maxd) :: rbox
 real, intent(inout), dimension(maxd_asize,maxd_atype) :: fr_partit_cw

 integer :: iphase, isize, itype
 integer :: jdone_mass, jdone_numb, jpos, jpos_mass, jpos_numb
 integer :: l, ll
 real, parameter :: partit_wght_mass = 0.5
 real :: dum, dummass, dumnumb, dumvolu
 real :: tmass, tnumb, umass, unumb, wmass, wnumb
 real, dimension(maxd_asize,maxd_atype) :: fmass, fnumb, xmass, xnumb

 iphase = cw_phase
 fmass(:,:) = 0.0
 fnumb(:,:) = 0.0
 xmass(:,:) = 0.0
 xnumb(:,:) = 0.0
 tmass = 0.0
 tnumb = 0.0
 umass = 0.0
 unumb = 0.0

!   compute
!     xmass, xnumb = mass, number mixing ratio for a bin
!     tmass, tnumb = sum over all bins of xmass, xnumb
!     umass, unumb = max over all bins of xmass, xnumb
!   set xmass, xnumb = 0.0 if bin mass, numb < 1.0e-37
!   constrain xnumb so that mean particle volume is
!     within bin boundaries
 do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)
       dummass = 0.0
       dumvolu = 0.0
       do ll = 1, ncomp_aer(itype)
          l = massptr_aer(ll,isize,itype,iphase)
          if (l >= param_first_scalar .and. l <= l2maxd) then
             dum = max(0.0, rbox(l))*mw_aer(ll,itype)
             dummass = dummass + dum
             dumvolu = dumvolu + dum/dens_aer(ll,itype)
          endif
       enddo

       l = numptr_aer(isize,itype,iphase)
       dumnumb = 0.0
       if (l >= param_first_scalar .and. l <= l2maxd) then
          dumnumb = max(0.0, rbox(l))
          if (dumnumb > dumvolu/volumlo_sect(isize,itype)) then
             dumnumb = dumvolu/volumlo_sect(isize,itype)
             rbox(l) = dumnumb
          else if (dumnumb < dumvolu/volumhi_sect(isize,itype)) then
             dumnumb = dumvolu/volumhi_sect(isize,itype)
             rbox(l) = dumnumb
          endif
       endif

       if (dummass < 1.0e-37) dummass = 0.0
       xmass(isize,itype) = dummass
       if (dumnumb < 1.0e-37) dumnumb = 0.0
       xnumb(isize,itype) = dumnumb

       tmass = tmass + xmass(isize,itype)
       tnumb = tnumb + xnumb(isize,itype)
       umass = max(umass, xmass(isize,itype))
       unumb = max(unumb, xnumb(isize,itype))
    enddo
 enddo

!   compute
!     fmass, fnumb = fraction of total mass, number that is in a bin
!   if tmass<1e-35 and umass>0, set fmass=1 for bin with largest xmass
!   if tmass<1e-35 and umass=0, set fmass=0 for all
 jdone_mass = 0
 jdone_numb = 0
 jpos_mass = 0
 jpos_numb = 0
 do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)
       fmass(isize,itype) = 0.0
       if (tmass >= 1.0e-35) then
          fmass(isize,itype) = xmass(isize,itype)/tmass
       else if (umass > 0.0) then
          if ((jdone_mass == 0) .and. (xmass(isize,itype) == umass)) then
             jdone_mass = 1
             fmass(isize,itype) = 1.0
          endif
       endif
       if (fmass(isize,itype) > 0.0) jpos_mass = jpos_mass + 1

       fnumb(isize,itype) = 0.0
       if (tnumb >= 1.0e-35) then
          fnumb(isize,itype) = xnumb(isize,itype)/tnumb
       else if (unumb > 0.0) then
          if ((jdone_numb == 0) .and. (xnumb(isize,itype) == unumb)) then
             jdone_numb = 1
             fnumb(isize,itype) = 1.0
          endif
       endif
       if (fnumb(isize,itype) > 0.0) jpos_numb = jpos_numb + 1
    enddo
 enddo

!   if only 1 bin has fmass or fnumb > 0, set value to 1.0 exactly
 if ((jpos_mass == 1) .or. (jpos_numb == 1)) then
    do itype = 1, ntype_aer
       do isize = 1, nsize_aer(itype)
          if ((jpos_mass == 1) .and. (fmass(isize,itype) > 0.0)) then
             fmass(isize,itype) = 1.0
          endif
          if ((jpos_numb == 1) .and. (fnumb(isize,itype) > 0.0)) then
             fnumb(isize,itype) = 1.0
          endif
       enddo
    enddo
 endif

!
!   compute fr_partit_cw as weighted average of fmass & fnumb, except
!   if tmass<1e-35 and umass=0, use only fnumb
!   if tnumb<1e-35 and unumb=0, use only fmass
!   if tmass,tnumb<1e-35 and umass,unumb=0, 
!       set fr_partit_cw=1 for center bin of itype=1
!
 fr_partit_cw(:,:) = 0.0
 if ((jpos_mass == 0) .and. (jpos_numb == 0)) then
    itype = 1
    isize = (nsize_aer(itype)+1)/2
    fr_partit_cw(isize,itype) = 1.0
 else if (jpos_mass == 0) then
    fr_partit_cw(:,:) = fnumb(:,:)
 else if (jpos_numb == 0) then
    fr_partit_cw(:,:) = fmass(:,:)
 else
    wmass = max(0.0, min(1.0, partit_wght_mass))
    wnumb = 1.0 - wmass
    fr_partit_cw(:,:) = wmass*fmass(:,:) + wnumb*fnumb(:,:)

    jpos = 0
    do itype = 1, ntype_aer
       do isize = 1, nsize_aer(itype)
          if (fr_partit_cw(isize,itype) > 0.0) jpos = jpos + 1
       enddo
    enddo

!   if only 1 bin has fr_partit_cw > 0, set value to 1.0 exactly
    if (jpos == 1) then
       do itype = 1, ntype_aer
          do isize = 1, nsize_aer(itype)
             if (fr_partit_cw(isize,itype) > 0.0) fr_partit_cw(isize,itype) = 1.0
          enddo
       enddo
    endif
 endif

 end subroutine partition_cldwtr

!=======================================================================
 subroutine distribute_bulk_changes(                                    &
      rbox, rbox_sv1, fr_partit_cw, rbulk_cwaer, lptr_yyy_cwaer,        &
      it, jt, kt, kpeg, mpeg, icase )
!=======================================================================
 use module_data_binaer_asect, only: maxd_asize, maxd_atype, nsize_aer, ntype_aer

 implicit none

 integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase
 integer, intent(in) :: lptr_yyy_cwaer(maxd_asize,maxd_atype,nyyy)
 real, intent(inout), dimension(1:l2maxd) :: rbox, rbox_sv1
 real, intent(in), dimension(maxd_asize,maxd_atype) :: fr_partit_cw
 real, intent(in), dimension(nyyy,2) :: rbulk_cwaer

 integer :: icount, idone, isize, itype, jpos, jpos_sv, l, lyyy, ncount
 real :: del_r_current, del_r_remain, fr, frsum_cur
 real :: duma, dumb
 real :: del_rbulk_cwaer(nyyy)
 real :: fr_cur(maxd_asize,maxd_atype)

 do lyyy = 1, nyyy
    del_rbulk_cwaer(lyyy) = rbulk_cwaer(lyyy,2) - rbulk_cwaer(lyyy,1)
 enddo

 jpos = 0
 do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)
       if (fr_partit_cw(isize,itype) > 0.0) jpos = jpos + 1
    enddo
 enddo
 jpos_sv = jpos

!
!   distribution is trivial when only 1 bin has fr_partit_cw > 0
!
 if (jpos_sv == 1) then
    do lyyy = 1, nyyy
       do itype = 1, ntype_aer
          do isize = 1, nsize_aer(itype)
             fr = fr_partit_cw(isize,itype)
             if (fr == 1.0) then
                l = lptr_yyy_cwaer(isize,itype,lyyy)
                if (l >= param_first_scalar .and. l <= l2maxd) rbox(l) = rbulk_cwaer(lyyy,2)
             endif
          enddo
       enddo
    enddo
    return
 endif

 do lyyy = 1, nyyy

!
!   distribution is simple when del_rbulk_cwaer(lyyy) >= 0
!
    if (del_rbulk_cwaer(lyyy) == 0.0) then
       cycle
    else if (del_rbulk_cwaer(lyyy) > 0.0) then
       do itype = 1, ntype_aer
          do isize = 1, nsize_aer(itype)
             fr = fr_partit_cw(isize,itype)
             if (fr > 0.0) then
                l = lptr_yyy_cwaer(isize,itype,lyyy)
                if (l >= param_first_scalar .and. l <= l2maxd) then
                   rbox(l) = rbox(l) + fr*del_rbulk_cwaer(lyyy)
                endif
             endif
          enddo
       enddo
       cycle
    endif

!
!   distribution is complicated when del_rbulk_cwaer(lyyy) < 0, 
!   because you cannot produce any negative mixrats
!
    del_r_remain = del_rbulk_cwaer(lyyy)
    fr_cur(:,:) = fr_partit_cw(:,:)
    ncount = max(1, jpos_sv*2)
    icount = 0

!   iteration loop
    do while (icount <= ncount)
       icount = icount + 1
       del_r_current = del_r_remain
       jpos = 0
       frsum_cur = 0.0

       do itype = 1, ntype_aer
          do isize = 1, nsize_aer(itype)
             fr = fr_cur(isize,itype)
             if (fr > 0.0) then
                l = lptr_yyy_cwaer(isize,itype,lyyy)
                if (l >= param_first_scalar .and. l <= l2maxd) then
                   duma = fr*del_r_current
                   dumb = rbox(l) + duma
                   if (dumb > 0.0) then
                      jpos = jpos + 1
                   else if (dumb == 0.0) then
                      fr_cur(isize,itype) = 0.0
                   else
                      duma = -rbox(l)
                      dumb = 0.0
                      fr_cur(isize,itype) = 0.0
                   endif
                   del_r_remain = del_r_remain - duma
                   rbox(l) = dumb
                   frsum_cur = frsum_cur + fr_cur(isize,itype)
                else
                   fr_cur(isize,itype) = 0.0
                endif
             endif
          enddo ! isize = 1, nsize_aer
       enddo ! itype = 1, ntype_aer

!   done if jpos = jpos_sv, because bins reached zero mixrat
       if (jpos == jpos_sv) then
          idone = 1
!   del_r_remain starts as negative, so done if non-negative
       else if (del_r_remain >= 0.0) then
          idone = 2
!   del_r_remain starts as negative, so done if non-negative
       else if (abs(del_r_remain) <= 1.0e-7*abs(del_rbulk_cwaer(lyyy))) then
          idone = 3
!   done if all bins have fr_cur = 0
       else if (frsum_cur <= 0.0) then
          idone = 4
!   same thing basically
       else if (jpos <= 0) then
          idone = 5
       else
          idone = 0
       endif

!   check for done, and (conditionally) print message
       if (idone > 0) exit

!   rescale fr_cur for next iteration
       fr_cur(:,:) = fr_cur(:,:)/frsum_cur
    enddo
 enddo

 end subroutine distribute_bulk_changes

!=======================================================================
 subroutine cloudchem_apply_move_sections(rbox, rbox_sv1, it, jt, kt, kpeg, mpeg, icase)
!=======================================================================
 use module_data_binaer_asect, only: msectional, maxd_asize, maxd_atype,&
      cw_phase, nsize_aer, ntype_aer, ncomp_aer, massptr_aer,           &
      mw_aer, dens_aer, drymass_aftgrow, drymass_pregrow,              &
      drydens_aftgrow, drydens_pregrow
 use module_mosaic_movesect, only: move_sections

 implicit none

 integer, intent(in) :: it, jt, kt, kpeg, mpeg, icase
 real, intent(inout), dimension(1:l2maxd) :: rbox, rbox_sv1

 integer :: idum_msect, iphase, isize, itype, l, ll
 real :: densdefault, dmaft, dmpre, dvaft, dvpre, duma
 real :: smallmassbb

 iphase = cw_phase

!
!   compute drymass before and after growth
!   set drydens = -1.0, so it will be calculated
!
 densdefault = 2.0
 smallmassbb = 1.0e-30

 do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)
       dmaft = 0.0
       dmpre = 0.0
       dvaft = 0.0
       dvpre = 0.0

       do ll = 1, ncomp_aer(itype)
          l = massptr_aer(ll,isize,itype,iphase)
          if (l >= param_first_scalar .and. l <= l2maxd) then
             duma = mw_aer(ll,itype)
             dmaft = dmaft + duma*rbox(l)
             dmpre = dmpre + duma*rbox_sv1(l)
             duma = duma/dens_aer(ll,itype)
             dvaft = dvaft + duma*rbox(l)
             dvpre = dvpre + duma*rbox_sv1(l)
          endif
       enddo

       drymass_aftgrow(isize,itype) = dmaft
       drymass_pregrow(isize,itype) = dmpre

       if (min(dmaft,dvaft) <= smallmassbb) then
          drydens_aftgrow(isize,itype) = densdefault
       else
          drydens_aftgrow(isize,itype) = dmaft/dvaft
       endif
       if (min(dmpre,dvpre) <= smallmassbb) then
          drydens_pregrow(isize,itype) = densdefault
       else
          drydens_pregrow(isize,itype) = dmpre/dvpre
       endif
    enddo
 enddo

!   apply move sections routine
!   (and conditionally turn on movesect diagnostics)
 idum_msect = msectional
 call move_sections(2, it, kpeg, mpeg)
 msectional = idum_msect

end subroutine cloudchem_apply_move_sections

end module module_cmu_cloudchem
