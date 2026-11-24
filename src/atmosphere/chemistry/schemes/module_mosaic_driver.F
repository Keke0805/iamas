!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module for MOSAIC (Model for Simulating Aerosol Interactions & Chemistry)
! * Author: Xueyin Ruan from LACAR USTC, ftrxyyg@mail.ustc.edu.cn
! Last update: xyruan 01/2026
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

 module module_mosaic_driver

 use iamas_atmchem_vars
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_atmchem_utilities
 use iamas_atmchem_configs
 use iamas_log, only : iamas_log_write

!
!   *** NOTE - when the cw species are NOT in the registry, then
!   then the p_xxx_cwnn variables are not in module_state_description,
!   and the following cpp directive should be commented out
!
#define cw_species_are_in_registry


 contains
!***********************************************************************
 subroutine MOSAIC_mech_driver(                         &
        id, curr_secs, ktau, dtstep, ktauc, dtstepc, configs, &
        t_phy, rho_phy, p_phy,                                &
        moist, chem, vbs_nbin,                                &
        its, ite, kts, kte                            )

 use module_data_binaer_asect
 use module_data_binaer_other
 use module_mosaic_therm, only:  aerchemistry, print_mosaic_stats
 use module_data_mosaic_therm, only:                     &
         iprint_mosaic_fe1, iprint_mosaic_perform_stats, &
         iprint_mosaic_diag1, iprint_mosaic_input_ok,    &
         n2o5_flag
 use module_mosaic_newnuc, only:  mosaic_newnuc_1clm
 use module_mosaic_coag, only:  mosaic_coag_1clm

 implicit none

!-----------------------------------------------------------------------
! DESCRIPTION
!
! mosaic_aerchem_driver is the interface between iAMAS and the
!   mosaic aerosol-chemistry routine cat computes condensation/evaporation
!   of trace gases to/from aerosol particles (AP).  It currently treats
!   water vapor and the 4 inorganic trace gases (nh3, h2so4, hno3, and hcl).
!   The aerosol-chemistry routine can work with either a sectional
!   (multiple size bins) or modal (multiple modes) representation.  
!
!   In both cases, condensation/evaporation to/from each bins/mode is 
!   first computed.  For sectional representation, AP mass and number 
!   are then transferred between size bins as a result of AP 
!   positive/negative growth.  Either a moving-center or two-moment
!   algorithm can be used to compute this transfer.
!
! mosaic_aerchem_driver is organized as follows
!   loop over i
!	call mapaer_tofrom_host to map 1 column of gas and aerosol mixing 
!	    ratios from the chem array to the rsub array (and convert units)
!	call aerchemistry to do the aerosol chemistry calculations
!	    for timestep = dtstepc
!	call mapaer_tofrom_host to map 1 column of gas and aerosol mixing 
!	    ratios from the rsub array back to the chem array
!
!-----------------------------------------------------------------------

 !subr argumentos
 type(iamas_pool_type),intent(in):: configs

 integer, intent(in) :: id, ktau, ktauc,      &
                        its, ite, kts, kte,   &
                        vbs_nbin(1)
                    
 !   id - domain index
 !   ktau - time step number
 !   ktauc - gas and aerosol chemistry time step number

 real(kind=8), intent(in) :: curr_secs
 real, intent(in) :: dtstep, dtstepc
 !   dtstep - main model time step (s)
 !   dtstepc - time step for gas and aerosol chemistry(s)

 real, dimension( kts:kte, its:ite), intent(in)::   &
                               t_phy, rho_phy, p_phy
 !   t_phy - temperature (K)
 !   rho_phy - air density (kg/m^3)
 !   p_phy - air pressure (Pa)

 real, dimension( 1:num_moist, kts:kte, its:ite), intent(in) :: moist
 !   moist - mixing ratios of moisture species (water vapor, 
 !	cloud water, ...) (kg/kg for mass species, #/kg for number species)
 
 real, dimension( 1:num_chem, kts:kte, its:ite ), intent(inout ) :: chem
 !   chem - mixing ratios of trace gase (ppm) and aerosol species
 !	(ug/kg for mass species, #/kg for number species)


!-----------------------------------------------------------------------
!   local variables
 integer :: i, idum, istat, it, j, jt, k, l, n
 integer :: k_pegshift, kclm_calcbgn, kclm_calcend
 integer :: ktmaps, ktmape
 integer :: levdbg_err, levdbg_info
 integer :: mode_force_dump
 integer :: idiagbb_dum
 integer, parameter :: debug_level=0
 integer, parameter :: aercoag_onoff = 1
 integer, parameter :: aernewnuc_onoff = 1

 real :: dtchem, dtcoag, dtnuc
 real :: dum
 real :: rsub0(l2maxd,kmaxd,nsubareamaxd)

 integer,pointer :: config_n2o5_hetchem !xyruan add

 character*100 msg

! xyruan- copy over N2O5 control flag
    call iamas_pool_get_config(configs,'config_n2o5_hetchem',config_n2o5_hetchem)
    n2o5_flag = config_n2o5_hetchem

! The default values for these informational printout settings are set
! in module_data_mosaic_therm.F.
    if (debug_level .lt. 15) then
       iprint_mosaic_fe1 = 1
       iprint_mosaic_perform_stats = 0
       iprint_mosaic_diag1 = 0
       iprint_mosaic_input_ok = 0
    end if

!   ktmaps,ktmape = first/last wrf kt for which aer chem is done
    ktmaps = kts
    ktmape = kte

!   rce 2005-mar-09 - added kclm_calcbgn/end 
!   kclm_calcbgn,kclm_calcend = first/last pegasus array k
!   for which aer chem is done
    k_pegshift = k_pegbegin - kts
    kclm_calcbgn = kts + k_pegshift
    kclm_calcend = kte + k_pegshift

!   set some variables to their wrf-chem "standard" values
    mode_force_dump = 0
    levdbg_err = 0
    levdbg_info = 15

!   eventually iymdcur & ihmscur should be set to the correct date/time 
!   using wrf timing routines
    iymdcur = 1 + int( curr_secs/86400._8 )
    ihmscur = nint( mod( curr_secs, 86400._8 ) )

    t = curr_secs
    ncorecnt = ktau - 1

!   set "pegasus" grid size variables
    itot = ite
    nsubareas = 1


    call print_mosaic_stats( 0 )

    do 2910 it = its, ite

    dtchem = dtstepc


    call mapaer_tofrom_host( 0,              &
        its,ite, kts,kte,                    &
        it, ktmaps,ktmape,                   &
        moist, chem,                         &
        t_phy, p_phy, rho_phy                 )

!rce 22-jul-2006 - save initial mixrats
    rsub0(:,:,:) = rsub(:,:,:)
    
    !idiagaa_dum = 0
    idiagbb_dum = 110

    call aerchemistry( it, kclm_calcbgn, kclm_calcend, dtchem, vbs_nbin )

!  note units for aerosol is now ug/m3

    call chemistry_debug_message("mosaic_aerchem_driver: back from aerchemistry")


     if (aernewnuc_onoff > 0) then
         dtnuc = dtchem
         call mosaic_newnuc_1clm( istat,      &
         it, kclm_calcbgn, kclm_calcend,      &
         idiagbb_dum, dtchem, dtnuc, rsub0,   &
         id, ktau, ktauc, its, ite, kts, kte )
     end if

     if (aercoag_onoff > 0) then
         dtcoag = dtchem
         call mosaic_coag_1clm( istat,        &
         it, kclm_calcbgn, kclm_calcend,      &
         idiagbb_dum, dtchem, dtcoag,         &
         id, ktau, ktauc, its, ite, kts, kte )
     end if

    call mapaer_tofrom_host( 1,              &
        its,ite, kts,kte,                    &
        it, ktmaps,ktmape,                   &
        moist, chem,                         &
        t_phy, p_phy, rho_phy                 )


2910    continue

    call print_mosaic_stats( 1 )

    write(msg, '(A,I0)') 'leaving mosaic_aerchem_driver - ktau =', ktau
    call chemistry_debug_message( msg )

    return
    end subroutine MOSAIC_mech_driver

    subroutine mapaer_tofrom_host( imap,        &
        its,ite, kts,kte,                       &
        it, ktmaps,ktmape,                      &
        moist, chem,                            &
        t_phy, p_phy, rho_phy                    )

    use module_data_binaer_asect
    use module_data_binaer_other
    use module_mosaic_csuesat, only:  esat_gchm

    implicit none

!   subr arguments

!   imap determines mapping direction (chem-->rsub if <=0, rsub-->chem if >0)
    integer, intent(in) :: imap
!   wrf array dimensions
    integer, intent(in) :: its, ite, kts, kte
!   do mapping for iAMAS k,i grid points = [ktmaps:ktmape,it]
    integer, intent(in) :: it, ktmaps, ktmape
!   

    real, dimension( 1:num_moist, kts:kte, its:ite), intent(in) :: moist
    real, dimension( 1:num_chem, kts:kte, its:ite ), intent(inout) :: chem
    real, dimension( kts:kte, its:ite), intent(in)::   &
                                               t_phy, rho_phy, p_phy

!   local variables
    integer ido_l, idum, iphase, itype,   &
        k, k1, k2, kt, kt1, kt2, k_pegshift, l, n
    integer p1st
    real dum, dumesat, dumrsat, dumrelhum, onemeps
    real factdens, factpres, factmoist, factgas,   &
        factaerso4, factaerno3, factaercl, factaermsa,   &
        factaerco3, factaernh4, factaerna, factaerca,   &
        factaeroin, factaeroc, factaerbc,factaerdust,   &
        factaertr1r1,factaertr1r2,factaertr1r3,factaertr1r4, &
       factaerpcg1_b_c,factaerpcg2_b_c,factaerpcg3_b_c,factaerpcg4_b_c,factaerpcg5_b_c,factaerpcg6_b_c, &
       factaerpcg7_b_c,factaerpcg8_b_c,factaerpcg9_b_c,factaeropcg1_b_c,factaeropcg2_b_c,factaeropcg3_b_c, &
       factaeropcg4_b_c,factaeropcg5_b_c,factaeropcg6_b_c,factaeropcg7_b_c,factaeropcg8_b_c,factaerpcg1_b_o,&
       factaerpcg2_b_o,factaerpcg3_b_o,factaerpcg4_b_o,factaerpcg5_b_o,factaerpcg6_b_o,factaerpcg7_b_o,&
       factaerpcg8_b_o,factaerpcg9_b_o,factaeropcg1_b_o,factaeropcg2_b_o,factaeropcg3_b_o,factaeropcg4_b_o, &
       factaeropcg5_b_o,factaeropcg6_b_o,factaeropcg7_b_o,factaeropcg8_b_o,factaerpcg1_f_c,factaerpcg2_f_c, &
       factaerpcg3_f_c,factaerpcg4_f_c,factaerpcg5_f_c,factaerpcg6_f_c,factaerpcg7_f_c,factaerpcg8_f_c, &
       factaerpcg9_f_c, factaeropcg1_f_c,factaeropcg2_f_c,factaeropcg3_f_c,factaeropcg4_f_c, &
       factaeropcg5_f_c,factaeropcg6_f_c,factaeropcg7_f_c,factaeropcg8_f_c,factaerpcg1_f_o, &
       factaerpcg2_f_o,factaerpcg3_f_o,factaerpcg4_f_o,factaerpcg5_f_o,factaerpcg6_f_o,factaerpcg7_f_o, &
       factaerpcg8_f_o,factaerpcg9_f_o,factaeropcg1_f_o,factaeropcg2_f_o,factaeropcg3_f_o,factaeropcg4_f_o,&
       factaeropcg5_f_o,factaeropcg6_f_o,factaeropcg7_f_o,factaeropcg8_f_o,&
       factaersmpa,factaersmpbb, &
       factaerant1_c,factaerant2_c,factaerant3_c,factaerant4_c, &
       factaerant1_o,factaerant2_o,factaerant3_o,factaerant4_o, &
       factaerbiog1_c,factaerbiog2_c,factaerbiog3_c,factaerbiog4_c, &
       factaerbiog1_o,factaerbiog2_o,factaerbiog3_o,factaerbiog4_o, &
        factaerhysw, factaerwater, factaernum

    real, parameter :: eps=0.622

    character*80 msg


!
!   units conversion factors 
!   wrf-chem value = pegasus value X factor
!
    factdens = 28.966e3      ! moleair/cm3 --> kgair/m3
    factpres = 0.1           ! dyne/cm2 --> pa
    factmoist = eps          ! moleh2o/moleair --> kgh2o/kgair
    factgas = 1.0e6          ! mole/moleair --> ppm

!wig 9-Nov-2004: Change to converting from concentration to converting
!                from mixing ratio.
!   factaernum = 40.9        ! #/moleair --> #/m3 at STP
!! at 1 atm & 298 k,  1 m3 = 40.9 moleair,  1 liter = 0.0409 moleair
    factaernum = 1000./28.966 ! 1 kg air = (1000/28.966) moleair

    dum = factaernum*1.0e6   ! g/moleair --> ug/m3 at STP
    factaerso4   = dum*mw_so4_aer
    factaerno3   = dum*mw_no3_aer
    factaercl    = dum*mw_cl_aer
    factaermsa   = dum*mw_msa_aer
    factaerco3   = dum*mw_co3_aer
    factaernh4   = dum*mw_nh4_aer
    factaerna    = dum*mw_na_aer
    factaerca    = dum*mw_ca_aer
    factaeroin   = dum*mw_oin_aer
    factaerdust  = dum*mw_dust_aer
    factaertr1r1  = dum*mw_tr1r1_aer
    factaertr1r2  = dum*mw_tr1r2_aer
    factaertr1r3  = dum*mw_tr1r3_aer
    factaertr1r4  = dum*mw_tr1r4_aer
    factaeroc    = dum*mw_oc_aer
    factaerbc    = dum*mw_bc_aer
    factaerhysw  = dum*mw_water_aer
    factaerwater = dum*mw_water_aer


!   rce 2005-mar-09 - set ktot in mapaer_tofrom_host;
!   use k_pegshift for calc of ktot and k (=k_peg)
!   k_pegshift = k index shift between wrf-chem and pegasus arrays
    k_pegshift = k_pegbegin - kts

!   set ktot = highest k index for pegasus arrays
!   since kts=1 and k_pegbegin=1, get k_pegshift=0 and ktot=kte-1
    ktot = (kte-1) + k_pegshift
!   *** check that ktot and kte <= kmaxd ***
    if ((kte > kmaxd) .or. (ktot > kmaxd) .or. (ktot <= 0)) then
        write( msg, '(a,4i5)' )   &
        '*** subr mapaer_tofrom_host -- ' //   &
        'ktot, kmaxd, kts, kte', ktot, kmaxd, kts, kte
        call chemistry_message( msg )
        msg = '*** subr mosaic_aerchem_driver -- ' //   &
        'kte>kmaxd OR ktot>kmaxd OR ktot<=0'
        call chemistry_error_fatal( msg )
    end if

!   rce 2005-apr-28 - changed mapping loops to improve memory access
!   now do rsub(l,k1:k2,m) <--> chem(l,kt1:kt2,it) for each species
    kt1 = ktmaps
    kt2 = ktmape
    k1 = kt1 + k_pegshift
    k2 = kt2 + k_pegshift

    if (imap .gt. 0) goto 2000
!   imap==0 -- map species and state variables from host arrays 
!              to rsub, cairclm, ptotclm

!   first zero everything (except relhumclm)
    rsub(:,:,:) = 0.0
    cairclm(:) = 0.0
    ptotclm(:) = 0.0
    afracsubarea(:,:) = 0.0
    relhumclm(:) = aboxtest_min_relhum
    rcldwtr_sub(:,:) = 0.0

    adrydens_sub( :,:,:,:) = 0.0
    aqmassdry_sub(:,:,:,:) = 0.0
    aqvoldry_sub( :,:,:,:) = 0.0

    p1st = param_first_scalar

            if (p_h2so4 .ge. p1st) then
                rsub(kh2so4,k1:k2,1) = chem(p_h2so4,kt1:kt2,it)/factgas
            elseif (p_sulf .ge. p1st)   then
                rsub(kh2so4,k1:k2,1) = chem(p_sulf,kt1:kt2,it)/factgas
            endif

        if (p_hno3 .ge. p1st)   &
        rsub(khno3,k1:k2,1)  = chem(p_hno3,kt1:kt2,it)/factgas
        if (p_hcl .ge. p1st)   &
        rsub(khcl,k1:k2,1)   = chem(p_hcl,kt1:kt2,it)/factgas
        if (p_nh3 .ge. p1st)   &
        rsub(knh3,k1:k2,1)   = chem(p_nh3,kt1:kt2,it)/factgas
        if (p_msa .ge. p1st)   &
        rsub(kmsa,k1:k2,1)   = chem(p_msa,kt1:kt2,it)/factgas
        if (p_n2o5 .ge. p1st)   &
        rsub(kn2o5,k1:k2,1)   = chem(p_n2o5,kt1:kt2,it)/factgas
        !xyruan: mosaic_port: No clno2 now!
        !if (p_clno2 .ge. p1st)   &
        !rsub(kclno2,k1:k2,1)   = chem(p_clno2,kt1:kt2,it)/factgas

!   rce 2005-apr-12 - added following species for cldchem, here and below:
!   ko3, kso2, kh2o2, khcho, khcooh, koh, kho2, 
!   kno3, kno, kno2, khono, kpan, kch3o2, kch3oh, kch3ooh
        if (p_o3 .ge. p1st)   &
        rsub(ko3,k1:k2,1)   = chem(p_o3,kt1:kt2,it)/factgas
        if (p_so2 .ge. p1st)   &
        rsub(kso2,k1:k2,1)   = chem(p_so2,kt1:kt2,it)/factgas
        if (p_h2o2 .ge. p1st)   &
        rsub(kh2o2,k1:k2,1)   = chem(p_h2o2,kt1:kt2,it)/factgas
        if (p_hcho .ge. p1st)   &
        rsub(khcho,k1:k2,1)   = chem(p_hcho,kt1:kt2,it)/factgas
        if (p_ho .ge. p1st)   &
        rsub(koh,k1:k2,1)   = chem(p_ho,kt1:kt2,it)/factgas
        if (p_ho2 .ge. p1st)   &
        rsub(kho2,k1:k2,1)   = chem(p_ho2,kt1:kt2,it)/factgas
        if (p_no3 .ge. p1st)   &
        rsub(kno3,k1:k2,1)   = chem(p_no3,kt1:kt2,it)/factgas
        if (p_no .ge. p1st)   &
        rsub(kno,k1:k2,1)   = chem(p_no,kt1:kt2,it)/factgas
        if (p_no2 .ge. p1st)   &
        rsub(kno2,k1:k2,1)   = chem(p_no2,kt1:kt2,it)/factgas
        if (p_hono .ge. p1st)   &
        rsub(khono,k1:k2,1)   = chem(p_hono,kt1:kt2,it)/factgas
        if (p_pan .ge. p1st)   &
        rsub(kpan,k1:k2,1)   = chem(p_pan,kt1:kt2,it)/factgas

        do iphase=1,nphase_aer
        do itype=1,ntype_aer
        do n = 1, nsize_aer(itype)
        rsub(lptr_so4_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_so4_aer(n,itype,iphase),kt1:kt2,it)/factaerso4
        rsub(numptr_aer(n,itype,iphase),k1:k2,1) =   &
            chem(numptr_aer(n,itype,iphase),kt1:kt2,it)/factaernum

        if (lptr_no3_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_no3_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_no3_aer(n,itype,iphase),kt1:kt2,it)/factaerno3
        if (lptr_cl_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_cl_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_cl_aer(n,itype,iphase),kt1:kt2,it)/factaercl
        if (lptr_msa_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_msa_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_msa_aer(n,itype,iphase),kt1:kt2,it)/factaermsa
        if (lptr_co3_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_co3_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_co3_aer(n,itype,iphase),kt1:kt2,it)/factaerco3
        if (lptr_nh4_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_nh4_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_nh4_aer(n,itype,iphase),kt1:kt2,it)/factaernh4
        if (lptr_na_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_na_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_na_aer(n,itype,iphase),kt1:kt2,it)/factaerna
        if (lptr_ca_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_ca_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_ca_aer(n,itype,iphase),kt1:kt2,it)/factaerca

        if (lptr_oin_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_oin_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_oin_aer(n,itype,iphase),kt1:kt2,it)/factaeroin
        if (lptr_dust_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_dust_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_dust_aer(n,itype,iphase),kt1:kt2,it)/factaerdust
        !if (lptr_tr1r1_aer(n,itype,iphase) .ge. p1st)   &
        !    rsub(lptr_tr1r1_aer(n,itype,iphase),k1:k2,1) =   &
        !    chem(lptr_tr1r1_aer(n,itype,iphase),kt1:kt2,it)/factaertr1r1
        !if (lptr_tr1r2_aer(n,itype,iphase) .ge. p1st)   &
        !    rsub(lptr_tr1r2_aer(n,itype,iphase),k1:k2,1) =   &
        !    chem(lptr_tr1r2_aer(n,itype,iphase),kt1:kt2,it)/factaertr1r2
        !if (lptr_tr1r3_aer(n,itype,iphase) .ge. p1st)   &
        !    rsub(lptr_tr1r3_aer(n,itype,iphase),k1:k2,1) =   &
        !    chem(lptr_tr1r3_aer(n,itype,iphase),kt1:kt2,it)/factaertr1r3
        !if (lptr_tr1r4_aer(n,itype,iphase) .ge. p1st)   &
        !    rsub(lptr_tr1r4_aer(n,itype,iphase),k1:k2,1) =   &
        !   chem(lptr_tr1r4_aer(n,itype,iphase),kt1:kt2,it)/factaertr1r4
        if (lptr_oc_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_oc_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_oc_aer(n,itype,iphase),kt1:kt2,it)/factaeroc
        if (lptr_bc_aer(n,itype,iphase) .ge. p1st)   &
            rsub(lptr_bc_aer(n,itype,iphase),k1:k2,1) =   &
            chem(lptr_bc_aer(n,itype,iphase),kt1:kt2,it)/factaerbc
        if (hyswptr_aer(n,itype) .ge. p1st)   &
            rsub(hyswptr_aer(n,itype),k1:k2,1) =   &
            chem(hyswptr_aer(n,itype),kt1:kt2,it)/factaerhysw
        if (waterptr_aer(n,itype) .ge. p1st)   &
            rsub(waterptr_aer(n,itype),k1:k2,1) =   &
            chem(waterptr_aer(n,itype),kt1:kt2,it)/factaerwater

        end do ! size
        end do ! type
        end do ! phase

!   map state variables
    afracsubarea(k1:k2,1) = 1.0
    rsub(ktemp,k1:k2,1) = t_phy(kt1:kt2,it)
    rsub(kh2o,k1:k2,1) = moist(p_qv,kt1:kt2,it)/factmoist
    cairclm(k1:k2) = rho_phy(kt1:kt2,it)/factdens
    ptotclm(k1:k2) = p_phy(kt1:kt2,it)/factpres
    if (p_qc .ge. p1st)   &
        rcldwtr_sub(k1:k2,1) = moist(p_qc,kt1:kt2,it)/factmoist

!   compute relhumclm from water vapor, temperature, and pressure
!   *** force relhum to between aboxtest_min/max_relhum
        do kt = ktmaps, ktmape
        k = kt + k_pegshift
        onemeps = 1.0 - 0.622
        dumesat = esat_gchm( rsub(ktemp,k,1) )
        dumrsat = dumesat / (ptotclm(k) - onemeps*dumesat)
        dumrelhum = rsub(kh2o,k,1) / max( dumrsat, 1.e-20 )
        dumrelhum = max( 0.0, min( 0.99, dumrelhum ) )

       relhumclm(k) = dumrelhum

        relhumclm(k) = max( relhumclm(k), aboxtest_min_relhum )
        relhumclm(k) = min( relhumclm(k), aboxtest_max_relhum )
        end do

!   *** force temperature to be > aboxtest_min_temp
    do kt = ktmaps, ktmape
        k = kt + k_pegshift
        rsub(ktemp,k,1) =   &
        max( rsub(ktemp,k,1), aboxtest_min_temp )
    end do


    return

!   imap==1 -- map species from rsub back to host arrays 

2000    continue

    p1st = param_first_scalar

            if (p_h2so4 .ge. p1st)   then
                chem(p_h2so4,kt1:kt2,it) = rsub(kh2so4,k1:k2,1)*factgas
           elseif (p_sulf .ge. p1st)   then   
                chem(p_sulf,kt1:kt2,it) = rsub(kh2so4,k1:k2,1)*factgas 
            endif       
        if (p_hno3 .ge. p1st)   &
        chem(p_hno3,kt1:kt2,it)  = rsub(khno3,k1:k2,1)*factgas
        if (p_hcl .ge. p1st)   &
        chem(p_hcl,kt1:kt2,it)   = rsub(khcl,k1:k2,1)*factgas 
        if (p_nh3 .ge. p1st)   &
        chem(p_nh3,kt1:kt2,it)  = rsub(knh3,k1:k2,1)*factgas
        if (p_msa .ge. p1st)   &
        chem(p_msa,kt1:kt2,it)  = rsub(kmsa,k1:k2,1)*factgas
        if (p_n2o5 .ge. p1st)   &
        chem(p_n2o5,kt1:kt2,it)  = rsub(kn2o5,k1:k2,1)*factgas
        !xyruan: mosaic_port: No clno2 now!
        !if (p_clno2 .ge. p1st)   &
        !chem(p_clno2,kt1:kt2,it)  = rsub(kclno2,k1:k2,1)*factgas

        if (p_o3 .ge. p1st)   &
        chem(p_o3,kt1:kt2,it)  = rsub(ko3,k1:k2,1)*factgas
        if (p_so2 .ge. p1st)   &
        chem(p_so2,kt1:kt2,it)  = rsub(kso2,k1:k2,1)*factgas
        if (p_h2o2 .ge. p1st)   &
        chem(p_h2o2,kt1:kt2,it)  = rsub(kh2o2,k1:k2,1)*factgas
        if (p_hcho .ge. p1st)   &
        chem(p_hcho,kt1:kt2,it)  = rsub(khcho,k1:k2,1)*factgas
        if (p_ho .ge. p1st)   &
        chem(p_ho,kt1:kt2,it)  = rsub(koh,k1:k2,1)*factgas
        if (p_ho2 .ge. p1st)   &
        chem(p_ho2,kt1:kt2,it)  = rsub(kho2,k1:k2,1)*factgas
        if (p_no3 .ge. p1st)   &
        chem(p_no3,kt1:kt2,it)  = rsub(kno3,k1:k2,1)*factgas
        if (p_no .ge. p1st)   &
        chem(p_no,kt1:kt2,it)  = rsub(kno,k1:k2,1)*factgas
        if (p_no2 .ge. p1st)   &
        chem(p_no2,kt1:kt2,it)  = rsub(kno2,k1:k2,1)*factgas
        if (p_hono .ge. p1st)   &
        chem(p_hono,kt1:kt2,it)  = rsub(khono,k1:k2,1)*factgas
        if (p_pan .ge. p1st)   &
        chem(p_pan,kt1:kt2,it)  = rsub(kpan,k1:k2,1)*factgas

        do iphase=1,nphase_aer
        do itype=1,ntype_aer
        do n = 1, nsize_aer(itype)
        chem(lptr_so4_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_so4_aer(n,itype,iphase),k1:k2,1)*factaerso4
        chem(numptr_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(numptr_aer(n,itype,iphase),k1:k2,1)*factaernum

        if (lptr_no3_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_no3_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_no3_aer(n,itype,iphase),k1:k2,1)*factaerno3
        if (lptr_cl_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_cl_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_cl_aer(n,itype,iphase),k1:k2,1)*factaercl
        if (lptr_msa_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_msa_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_msa_aer(n,itype,iphase),k1:k2,1)*factaermsa
        if (lptr_co3_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_co3_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_co3_aer(n,itype,iphase),k1:k2,1)*factaerco3
        if (lptr_nh4_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_nh4_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_nh4_aer(n,itype,iphase),k1:k2,1)*factaernh4
        if (lptr_na_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_na_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_na_aer(n,itype,iphase),k1:k2,1)*factaerna
        if (lptr_ca_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_ca_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_ca_aer(n,itype,iphase),k1:k2,1)*factaerca
        if (lptr_oin_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_oin_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_oin_aer(n,itype,iphase),k1:k2,1)*factaeroin
        if (lptr_dust_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_dust_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_dust_aer(n,itype,iphase),k1:k2,1)*factaerdust
        ! if (lptr_tr1r1_aer(n,itype,iphase) .ge. p1st)   &
        !     chem(lptr_tr1r1_aer(n,itype,iphase),kt1:kt2,it) =   &
        !     rsub(lptr_tr1r1_aer(n,itype,iphase),k1:k2,1)*factaertr1r1
        ! if (lptr_tr1r2_aer(n,itype,iphase) .ge. p1st)   &
        !     chem(lptr_tr1r2_aer(n,itype,iphase),kt1:kt2,it) =   &
        !     rsub(lptr_tr1r2_aer(n,itype,iphase),k1:k2,1)*factaertr1r2
        ! if (lptr_tr1r3_aer(n,itype,iphase) .ge. p1st)   &
        !     chem(lptr_tr1r3_aer(n,itype,iphase),kt1:kt2,it) =   &
        !     rsub(lptr_tr1r3_aer(n,itype,iphase),k1:k2,1)*factaertr1r3
        ! if (lptr_tr1r4_aer(n,itype,iphase) .ge. p1st)   &
        !     chem(lptr_tr1r4_aer(n,itype,iphase),kt1:kt2,it) =   &
        !     rsub(lptr_tr1r4_aer(n,itype,iphase),k1:k2,1)*factaertr1r4
        if (lptr_oc_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_oc_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_oc_aer(n,itype,iphase),k1:k2,1)*factaeroc
        if (lptr_bc_aer(n,itype,iphase) .ge. p1st)   &
            chem(lptr_bc_aer(n,itype,iphase),kt1:kt2,it) =   &
            rsub(lptr_bc_aer(n,itype,iphase),k1:k2,1)*factaerbc
        if (hyswptr_aer(n,itype) .ge. p1st)   &
            chem(hyswptr_aer(n,itype),kt1:kt2,it) =   &
            rsub(hyswptr_aer(n,itype),k1:k2,1)*factaerhysw
        if (waterptr_aer(n,itype) .ge. p1st)   &
            chem(waterptr_aer(n,itype),kt1:kt2,it) =   &
            rsub(waterptr_aer(n,itype),k1:k2,1)*factaerwater

        end do ! size
        end do ! type
        end do ! phase

    return

    end subroutine mapaer_tofrom_host

!-----------------------------------------------------------------------
    end module module_mosaic_driver
