!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to declare or initialize variables needed by chemistry processes. 
! * Author: Jiawang Feng, Zihan Xia, Xueyin Ruan and Chun Zhao
! Last update: Jan 2026
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao

! References: 
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************


module iamas_atmchem_binaer_init
use iamas_kind_types
use iamas_pool_routines

implicit none
public :: init_data_binaer_asect
save

contains
!-----------------------------------------------------------------------
! *** note - eventually is_aerosol will be a subr argument
subroutine init_data_binaer_asect( state, tend_chemistry, chem_opt, is_aerosol )
!   subroutine init_data_binaer_asect( id,n2o5_hetchem, is_aerosol )
!    subroutine init_data_mosaic_asect( )

    use module_data_binaer_asect
    use module_data_binaer_other, only:  lunerr, lunout,   &
        aboxtest_testmode, aboxtest_units_convert,   &
        aboxtest_rh_method, aboxtest_map_method,   &
        aboxtest_gases_fixed, aboxtest_min_temp,   &
        aboxtest_min_relhum, aboxtest_max_relhum
    use module_data_mosaic_therm, only:  nbin_a, nbin_a_maxd, n2o5_flag
    use module_mosaic_csuesat, only:  init_csuesat
    use module_mosaic_movesect, only:  move_sections, test_move_sections
    ! use module_peg_util, only:  peg_error_fatal
    use iamas_atmchem_utilities, only: chemistry_error_fatal, chemistry_message

    use iamas_atmchem_vars, only:  param_first_scalar, num_chem , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05   , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10   , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05  , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10

    implicit none

! *** note - eventually is_aerosol will be a subr argument
    type(iamas_pool_type),intent(in) :: state
    type(iamas_pool_type),intent(in) :: tend_chemistry
    character(len=*),intent(IN) :: chem_opt
    logical, intent(out) :: is_aerosol(num_chem)

!   local variables
    integer idum, itype, l, ldum, n, nhi, nsize_aer_dum
    real dum
    real, parameter :: pi = 3.14159265
    character*200 msg


!
!   set indexs to iAMAS chem-array species
!
    call init_data_binaer_index(state, tend_chemistry)

!
!   set some "pegasus" control variables
!
    msectional = 20
    maerocoag = -2
    maerchem = 1
    maeroptical = 1
    maerchem_boxtest_output = -1

!
!   set ntype_aer = 1
!
    ntype_aer = 1

!
!   set number of aerosol bins using the wrf-chem sulfate pointers
!
    nsize_aer(:) = 0
    itype=1
! jwfeng from AEMOL USTC
! Assume that all schemes have num_a0*.p_num_a0* will be -1 if num_a0* is not included in.
    if (p_num_a01 .ge. param_first_scalar) nsize_aer(itype) = 1
    if (p_num_a02 .ge. param_first_scalar) nsize_aer(itype) = 2
    if (p_num_a03 .ge. param_first_scalar) nsize_aer(itype) = 3
    if (p_num_a04 .ge. param_first_scalar) nsize_aer(itype) = 4
    if (p_num_a05 .ge. param_first_scalar) nsize_aer(itype) = 5
    if (p_num_a06 .ge. param_first_scalar) nsize_aer(itype) = 6
    if (p_num_a07 .ge. param_first_scalar) nsize_aer(itype) = 7
    if (p_num_a08 .ge. param_first_scalar) nsize_aer(itype) = 8
    if (p_num_a09 .ge. param_first_scalar) nsize_aer(itype) = 9
    if (p_num_a10 .ge. param_first_scalar) nsize_aer(itype) = 10

    if (nsize_aer(itype) .le. 0) then
       call chemistry_error_fatal(&
             'init_data_binaer_asect - nsize_aer = 0')
    else if (nsize_aer(itype) .gt. maxd_asize) then
       call chemistry_error_fatal(&
             'init_data_binaer_asect - nsize_aer > maxd_asize' )
    end if

!
!   set nbin_a to total number of aerosol bins (for all types)
!
    nbin_a = 0
    do itype = 1, ntype_aer
        nbin_a = nbin_a + nsize_aer(itype)
    end do
    write(msg,'(a,(1x,I4))') '--- nbin_a = ' , nbin_a
    call chemistry_message(msg)
    if (nbin_a .gt. nbin_a_maxd) then
        call chemistry_error_fatal(&
            'init_data_binaer_asect - nbin_a > nbin_a_maxd' )
    end if

!
!   set nphase_aer (number of active aerosol species phases),
!   the xx_phase, and maerosolincw
!
    nphase_aer = 0
    maerosolincw = 0
    if (nsize_aer(1) .gt. 0) then
        nphase_aer = 1
        ai_phase = 1

        ! jwfeng from AEMOL USTC:
        ! Assume that all schemes have num_a0*.p_num_a0* will be -1 if num_a0* is not included in.
        if (p_num_cw01 .ge. param_first_scalar) then
            nphase_aer = 2
            cw_phase = 2
            maerosolincw = 1
        end if
    end if


! #if defined ( aboxtest_box_testing_active )
! ! *** start of "box testing" code section ***
! !     these code lines should be inactive when running wrf-chem
! !
! !   set some variables to "box test" values
!      call aboxtest_get_extra_args( 10,   &
!          msectional, maerosolincw, maerocoag,   &
!          maerchem, maeroptical, maerchem_boxtest_output,   &
!          lunerr, lunout, idum, dum, dum )
!      call aboxtest_get_extra_args( 11,   &
!          aboxtest_testmode, aboxtest_units_convert,   &
!          aboxtest_rh_method, aboxtest_map_method,   &
!          aboxtest_gases_fixed, nsize_aer_dum,   &
!          idum, idum, idum, dum, dum )
 
!      itype = 1
!      if (nsize_aer_dum > 0) nsize_aer(itype) = nsize_aer_dum
 
!      aboxtest_min_temp = 0.0
!      aboxtest_min_relhum = 0.0
!      aboxtest_max_relhum = 1.0
! ! ***  end  of "box testing" code section ***
! #endif


!
!   set master aerosol chemical types
!
    !ntot_mastercomp_aer = 96  !czhao add dust or tr1r1-4
    ntot_mastercomp_aer = 13 !zhxia :When we add new aerosol species, we need to consider whether to update this and massptr. 10/2024

        l = 1 
        mastercompindx_oc_aer    = l    
        name_mastercomp_aer( l ) = 'oc'
        dens_mastercomp_aer( l ) =  dens_oc_aer
        mw_mastercomp_aer(   l ) =    mw_oc_aer
        hygro_mastercomp_aer(l ) = hygro_oc_aer

        l = 2 
        mastercompindx_bc_aer    = l    
        name_mastercomp_aer( l ) = 'bc'
        dens_mastercomp_aer( l ) =  dens_bc_aer
        mw_mastercomp_aer(   l ) =    mw_bc_aer
        hygro_mastercomp_aer(l ) = hygro_bc_aer

        l = 3 
        mastercompindx_so4_aer    = l    
        name_mastercomp_aer( l ) = 'so4'
        dens_mastercomp_aer( l ) =  dens_so4_aer
        mw_mastercomp_aer(   l ) =    mw_so4_aer
        hygro_mastercomp_aer(l ) = hygro_so4_aer

        l = 4 
        mastercompindx_no3_aer    = l    
        name_mastercomp_aer( l ) = 'no3'
        dens_mastercomp_aer( l ) =  dens_no3_aer
        mw_mastercomp_aer(   l ) =    mw_no3_aer
        hygro_mastercomp_aer(l ) = hygro_no3_aer

        l = 5 
        mastercompindx_oin_aer    = l    
        name_mastercomp_aer( l ) = 'oin'
        dens_mastercomp_aer( l ) =  dens_oin_aer
        mw_mastercomp_aer(   l ) =    mw_oin_aer
        hygro_mastercomp_aer(l ) = hygro_oin_aer

        l = 6 
        mastercompindx_dust_aer = l    
        name_mastercomp_aer( l ) = 'dust'
        dens_mastercomp_aer( l ) =  dens_dust_aer
        mw_mastercomp_aer(   l ) =    mw_dust_aer
        hygro_mastercomp_aer(l ) = hygro_dust_aer

        l = 7
        mastercompindx_nh4_aer = l
        name_mastercomp_aer( l ) = 'nh4'
        dens_mastercomp_aer( l ) =  dens_nh4_aer
        mw_mastercomp_aer(   l ) =    mw_nh4_aer
        hygro_mastercomp_aer(l ) = hygro_nh4_aer

        l = 8
        mastercompindx_msa_aer = l
        name_mastercomp_aer( l ) = 'msa'
        dens_mastercomp_aer( l ) =  dens_msa_aer
        mw_mastercomp_aer(   l ) =    mw_msa_aer
        hygro_mastercomp_aer(l ) = hygro_msa_aer

        l = 9
        mastercompindx_soa_aer = l
        name_mastercomp_aer( l ) = 'soa'
        dens_mastercomp_aer( l ) =  dens_soa_aer
        mw_mastercomp_aer(   l ) =    mw_soa_aer
        hygro_mastercomp_aer(l ) = hygro_soa_aer

        l = 10
        mastercompindx_na_aer = l
        name_mastercomp_aer( l ) = 'sodium'
        dens_mastercomp_aer( l ) =  dens_na_aer
        mw_mastercomp_aer(   l ) =    mw_na_aer
        hygro_mastercomp_aer(l ) = hygro_na_aer

        l = 11
        mastercompindx_cl_aer = l
        name_mastercomp_aer( l ) = 'chloride'
        dens_mastercomp_aer( l ) =  dens_cl_aer
        mw_mastercomp_aer(   l ) =    mw_cl_aer
        hygro_mastercomp_aer(l ) = hygro_cl_aer

        l = 12
        mastercompindx_ca_aer = l
        name_mastercomp_aer( l ) = 'calcium'
        dens_mastercomp_aer( l ) =  dens_ca_aer
        mw_mastercomp_aer(   l ) =    mw_ca_aer
        hygro_mastercomp_aer(l ) = hygro_ca_aer

        l = 13
        mastercompindx_co3_aer = l
        name_mastercomp_aer( l ) = 'carbonate'
        dens_mastercomp_aer( l ) =  dens_co3_aer
        mw_mastercomp_aer(   l ) =    mw_co3_aer
        hygro_mastercomp_aer(l ) = hygro_co3_aer



!
!   set section size arrays
! 
    do itype = 1, ntype_aer 
        nhi = nsize_aer(itype) 
        dlo_sect(1,itype) = 3.90625e-6 
        if (nhi.eq.4 .or. nhi.eq.8) then
            dhi_sect(nhi,itype) = 10.0e-4 
        else if (nhi.eq.5 .or. nhi.eq.10) then
            dhi_sect(nhi,itype) = 40.0e-4 
        else if (nhi.eq.9) then
            dhi_sect(nhi,itype) = 20.0e-4 
        else
            write( msg,*) 'Wrong nsize_aer:', nhi
            call chemistry_error_fatal(msg)
        end if
        dum = log( dhi_sect(nhi,itype)/dlo_sect(1,itype) ) / nhi !replaced alog by log by Manish Shrivastava on 11/28/2011. alog denoted natural log in fortran 77. log(x) is natural log in fortran 90
        do n = 2, nhi
                dlo_sect(n,itype) = dlo_sect(1,itype) * exp( (n-1)*dum )
                dhi_sect(n-1,itype) = dlo_sect(n,itype)
        end do
        do n = 1, nhi
            dcen_sect(n,itype) = sqrt( dlo_sect(n,itype)*dhi_sect(n,itype) )
            volumlo_sect(n,itype) = (pi/6.) * (dlo_sect(n,itype)**3)
            volumhi_sect(n,itype) = (pi/6.) * (dhi_sect(n,itype)**3)
            volumcen_sect(n,itype) = (pi/6.) * (dcen_sect(n,itype)**3)
            sigmag_aer(n,itype) = (dhi_sect(n,itype)/dlo_sect(n,itype))**0.289
        end do
        write(msg,'(a,10(1x,f8.3))') '--- dlo_sect  (um):' , dlo_sect(:,itype) *1.0e4
        call chemistry_message(msg)
        write(msg,'(a,10(1x,f8.3))') '--- dhi_sect  (um):' , dhi_sect(:,itype) *1.0e4
        call chemistry_message(msg)
        write(msg,'(a,10(1x,f8.3))') '--- dcen_sect (um):' , dcen_sect(:,itype)*1.0e4
        call chemistry_message(msg)
    end do

!
!   set pointers to iAMAS chem-array species
!
    call init_data_binaer_ptr( is_aerosol )

!
!   csuesat initialization
!
    call init_csuesat

!
!   move_sect initialization (and testing)
!
!    subr move_sections( iflag, iclm, k, m )
    call move_sections(    -1,    1,    1, 1 )

    call test_move_sections( 1,   1,    1, 1 )
    

end subroutine init_data_binaer_asect

subroutine init_data_binaer_index(state, tend_chemistry)
    use iamas_kind_types
    use iamas_pool_routines
    use iamas_atmchem_vars, only: &
        p_qv,        p_qc,    p_qr,    p_qi,        p_qs       , &
        p_qndrop,    p_qniceaer,  p_qnicemix                   , &
        p_dust_a01,  p_dust_a02,  p_dust_a03,  p_dust_a04,  p_dust_a05   , &
        p_dust_a06,  p_dust_a07,  p_dust_a08,  p_dust_a09,  p_dust_a10   , &
        p_water_a01, p_water_a02, p_water_a03, p_water_a04, p_water_a05  , &
        p_water_a06, p_water_a07, p_water_a08, p_water_a09, p_water_a10  , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05    , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10    , &
        p_dust_cw01, p_dust_cw02, p_dust_cw03, p_dust_cw04, p_dust_cw05  , &
        p_dust_cw06, p_dust_cw07, p_dust_cw08, p_dust_cw09, p_dust_cw10  , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05   , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10   , &
        p_oc_a01,    p_oc_a02,    p_oc_a03,    p_oc_a04                  , &
        p_bc_a01,    p_bc_a02,    p_bc_a03,    p_bc_a04                  , &
        p_so4_a01,   p_so4_a02,   p_so4_a03,   p_so4_a04                 , &
        p_nh4_a01,   p_nh4_a02,   p_nh4_a03,   p_nh4_a04                 , &
        p_msa_a01,   p_msa_a02,   p_msa_a03,   p_msa_a04                 , &
        p_no3_a01,   p_no3_a02,   p_no3_a03,   p_no3_a04                 , &
        p_oin_a01,   p_oin_a02,   p_oin_a03,   p_oin_a04                 , &
        p_soa_a01,   p_soa_a02,   p_soa_a03,   p_soa_a04                 , &
        p_na_a01,    p_na_a02,    p_na_a03,    p_na_a04                  , &
        p_na_a05,    p_na_a06,    p_na_a07,    p_na_a08,    p_na_a09,    p_na_a10, &
        p_cl_a01,    p_cl_a02,    p_cl_a03,    p_cl_a04                  , &
        p_cl_a05,    p_cl_a06,    p_cl_a07,    p_cl_a08,    p_cl_a09,    p_cl_a10, &
        p_ca_a01,    p_ca_a02,    p_ca_a03,    p_ca_a04                  , &
        p_co3_a01,   p_co3_a02,   p_co3_a03,   p_co3_a04                 , &
        p_hysw_a01,  p_hysw_a02,  p_hysw_a03,  p_hysw_a04                , &
        p_oc_cw01,    p_oc_cw02,    p_oc_cw03,    p_oc_cw04              , &
        p_bc_cw01,    p_bc_cw02,    p_bc_cw03,    p_bc_cw04              , &
        p_so4_cw01,   p_so4_cw02,   p_so4_cw03,   p_so4_cw04             , &
        p_nh4_cw01,   p_nh4_cw02,   p_nh4_cw03,   p_nh4_cw04             , &
        p_msa_cw01,   p_msa_cw02,   p_msa_cw03,   p_msa_cw04             , &
        p_no3_cw01,   p_no3_cw02,   p_no3_cw03,   p_no3_cw04             , &
        p_oin_cw01,   p_oin_cw02,   p_oin_cw03,   p_oin_cw04             , &
        p_soa_cw01,   p_soa_cw02,   p_soa_cw03,   p_soa_cw04             , &
        p_na_cw01,    p_na_cw02,    p_na_cw03,    p_na_cw04              , &
        p_na_cw05,    p_na_cw06,    p_na_cw07,    p_na_cw08,    p_na_cw09,    p_na_cw10, &
        p_cl_cw01,    p_cl_cw02,    p_cl_cw03,    p_cl_cw04              , &
        p_cl_cw05,    p_cl_cw06,    p_cl_cw07,    p_cl_cw08,    p_cl_cw09,    p_cl_cw10, &
        p_ca_cw01,    p_ca_cw02,    p_ca_cw03,    p_ca_cw04              , &
        p_co3_cw01,   p_co3_cw02,   p_co3_cw03,   p_co3_cw04             , &
        p_co2,       p_no,        p_no2,         p_so2                   , &
        p_nh3,       p_co,        p_o3,          p_h2o2                  , &
        p_no3,       p_n2o5,      p_hono,        p_hno3                  , &
        p_hno4,      p_h2so4,     p_hcho,        p_ccho                  , &
        p_rcho,      p_acet,      p_mek,         p_hcooh                 , &
        p_meoh,      p_etoh,      p_cco_oh,      p_rco_oh                , &
        p_gly,       p_mgly,      p_bacl,        p_cres                  , &
        p_bald,      p_isoprod,   p_methacro,    p_mvk                   , &
        p_prod2,     p_dcb1,      p_dcb2,        p_dcb3                  , &
        p_ethene,    p_isoprene,  p_c2h6,        p_c3h8                  , &
        p_c2h2,      p_c3h6,      p_alk3,        p_alk4                  , &
        p_alk5,      p_aro1,      p_aro2,        p_ole1                  , &
        p_ole2,      p_terp,      p_rno3,        p_nphe                  , &
        p_phen,      p_pan,       p_pan2,        p_pbzn                  , &
        p_ma_pan,    p_cco_ooh,   p_rco_o2,      p_rco_ooh               , &
        p_xn,        p_xc,        p_ho,          p_ho2                   , &
        p_c_o2,      p_cooh,      p_rooh,        p_ro2_r                 , &
        p_r2o2,      p_ro2_n,     p_cco_o2,      p_bzco_o2               , &
        p_ma_rco3,   p_sesq,      p_ch4,         p_hcl                   , &
        p_msa,       p_sulf                                              , &
        p_alk1,      p_alk2                                              , &
        p_ph_o2,       p_ph_o31d,     p_ph_o33p,      p_ph_no2           , &
        p_ph_no3o2,    p_ph_no3o,     p_ph_hno2,      p_ph_hno3          , &
        p_ph_hno4,     p_ph_h2o2,     p_ph_ch2or,     p_ph_ch2om         , &
        p_ph_ch3cho,   p_ph_ch3coch3, p_ph_ch3coc2h5, p_ph_hcocho        , &
        p_ph_ch3cocho, p_ph_hcochest, p_ph_ch3o2h,    p_ph_ch3coo2h      , &
        p_ph_ch3ono2,  p_ph_hcochob,  p_ph_n2o5



    implicit none
    type(iamas_pool_type),intent(in) :: state
    type(iamas_pool_type),intent(in) :: tend_chemistry

!jwfeng: 12/31/2019 
!call iamas subroutine to pass index_* to p_* 
!physics
    call iamas_pool_get_dimension(state, 'index_qv'      , p_qv           )
    call iamas_pool_get_dimension(state, 'index_qc'      , p_qc           )
    call iamas_pool_get_dimension(state, 'index_qr'      , p_qr           )
    call iamas_pool_get_dimension(state, 'index_qi'      , p_qi           )
    call iamas_pool_get_dimension(state, 'index_qs'      , p_qs           )
    call iamas_pool_get_dimension(state, 'index_qndrop'  , p_qndrop       )
    call iamas_pool_get_dimension(state, 'index_qniceaer', p_qniceaer     )
    call iamas_pool_get_dimension(state, 'index_qnicemix', p_qnicemix     )

!chemistry
    call iamas_pool_get_dimension(state, 'index_dust_a01', p_dust_a01     )
    call iamas_pool_get_dimension(state, 'index_dust_a02', p_dust_a02     )
    call iamas_pool_get_dimension(state, 'index_dust_a03', p_dust_a03     )
    call iamas_pool_get_dimension(state, 'index_dust_a04', p_dust_a04     )
    call iamas_pool_get_dimension(state, 'index_dust_a05', p_dust_a05     )
    call iamas_pool_get_dimension(state, 'index_dust_a06', p_dust_a06     )
    call iamas_pool_get_dimension(state, 'index_dust_a07', p_dust_a07     )
    call iamas_pool_get_dimension(state, 'index_dust_a08', p_dust_a08     )
    call iamas_pool_get_dimension(state, 'index_dust_a09', p_dust_a09     )
    call iamas_pool_get_dimension(state, 'index_dust_a10', p_dust_a10     )

    call iamas_pool_get_dimension(state, 'index_water_a01', p_water_a01   )
    call iamas_pool_get_dimension(state, 'index_water_a02', p_water_a02   )
    call iamas_pool_get_dimension(state, 'index_water_a03', p_water_a03   )
    call iamas_pool_get_dimension(state, 'index_water_a04', p_water_a04   )
    call iamas_pool_get_dimension(state, 'index_water_a05', p_water_a05   )
    call iamas_pool_get_dimension(state, 'index_water_a06', p_water_a06   )
    call iamas_pool_get_dimension(state, 'index_water_a07', p_water_a07   )
    call iamas_pool_get_dimension(state, 'index_water_a08', p_water_a08   )
    call iamas_pool_get_dimension(state, 'index_water_a09', p_water_a09   )
    call iamas_pool_get_dimension(state, 'index_water_a10', p_water_a10   )

    call iamas_pool_get_dimension(state, 'index_num_a01', p_num_a01     )
    call iamas_pool_get_dimension(state, 'index_num_a02', p_num_a02     )
    call iamas_pool_get_dimension(state, 'index_num_a03', p_num_a03     )
    call iamas_pool_get_dimension(state, 'index_num_a04', p_num_a04     )
    call iamas_pool_get_dimension(state, 'index_num_a05', p_num_a05     )
    call iamas_pool_get_dimension(state, 'index_num_a06', p_num_a06     )
    call iamas_pool_get_dimension(state, 'index_num_a07', p_num_a07     )
    call iamas_pool_get_dimension(state, 'index_num_a08', p_num_a08     )
    call iamas_pool_get_dimension(state, 'index_num_a09', p_num_a09     )
    call iamas_pool_get_dimension(state, 'index_num_a10', p_num_a10     )
 
    call iamas_pool_get_dimension(state, 'index_oc_a01',   p_oc_a01       )
    call iamas_pool_get_dimension(state, 'index_oc_a02',   p_oc_a02       )
    call iamas_pool_get_dimension(state, 'index_oc_a03',   p_oc_a03       )
    call iamas_pool_get_dimension(state, 'index_oc_a04',   p_oc_a04       )

    call iamas_pool_get_dimension(state, 'index_bc_a01',   p_bc_a01       )
    call iamas_pool_get_dimension(state, 'index_bc_a02',   p_bc_a02       )
    call iamas_pool_get_dimension(state, 'index_bc_a03',   p_bc_a03       )
    call iamas_pool_get_dimension(state, 'index_bc_a04',   p_bc_a04       )

    call iamas_pool_get_dimension(state, 'index_so4_a01',  p_so4_a01      )
    call iamas_pool_get_dimension(state, 'index_so4_a02',  p_so4_a02      )
    call iamas_pool_get_dimension(state, 'index_so4_a03',  p_so4_a03      )
    call iamas_pool_get_dimension(state, 'index_so4_a04',  p_so4_a04      )

    call iamas_pool_get_dimension(state, 'index_nh4_a01',  p_nh4_a01      )
    call iamas_pool_get_dimension(state, 'index_nh4_a02',  p_nh4_a02      )
    call iamas_pool_get_dimension(state, 'index_nh4_a03',  p_nh4_a03      )
    call iamas_pool_get_dimension(state, 'index_nh4_a04',  p_nh4_a04      )

    call iamas_pool_get_dimension(state, 'index_msa_a01',  p_msa_a01      )
    call iamas_pool_get_dimension(state, 'index_msa_a02',  p_msa_a02      )
    call iamas_pool_get_dimension(state, 'index_msa_a03',  p_msa_a03      )
    call iamas_pool_get_dimension(state, 'index_msa_a04',  p_msa_a04      )

    call iamas_pool_get_dimension(state, 'index_no3_a01',  p_no3_a01      )
    call iamas_pool_get_dimension(state, 'index_no3_a02',  p_no3_a02      )
    call iamas_pool_get_dimension(state, 'index_no3_a03',  p_no3_a03      )
    call iamas_pool_get_dimension(state, 'index_no3_a04',  p_no3_a04      )

    call iamas_pool_get_dimension(state, 'index_oin_a01',  p_oin_a01      )
    call iamas_pool_get_dimension(state, 'index_oin_a02',  p_oin_a02      )
    call iamas_pool_get_dimension(state, 'index_oin_a03',  p_oin_a03      )
    call iamas_pool_get_dimension(state, 'index_oin_a04',  p_oin_a04      )

    call iamas_pool_get_dimension(state, 'index_soa_a01',  p_soa_a01      )
    call iamas_pool_get_dimension(state, 'index_soa_a02',  p_soa_a02      )
    call iamas_pool_get_dimension(state, 'index_soa_a03',  p_soa_a03      )
    call iamas_pool_get_dimension(state, 'index_soa_a04',  p_soa_a04      )

    call iamas_pool_get_dimension(state, 'index_na_a01',   p_na_a01       )
    call iamas_pool_get_dimension(state, 'index_na_a02',   p_na_a02       )
    call iamas_pool_get_dimension(state, 'index_na_a03',   p_na_a03       )
    call iamas_pool_get_dimension(state, 'index_na_a04',   p_na_a04       )
    call iamas_pool_get_dimension(state, 'index_na_a05',   p_na_a05       )
    call iamas_pool_get_dimension(state, 'index_na_a06',   p_na_a06       )
    call iamas_pool_get_dimension(state, 'index_na_a07',   p_na_a07       )
    call iamas_pool_get_dimension(state, 'index_na_a08',   p_na_a08       )
    call iamas_pool_get_dimension(state, 'index_na_a09',   p_na_a09       )
    call iamas_pool_get_dimension(state, 'index_na_a10',   p_na_a10       )

    call iamas_pool_get_dimension(state, 'index_cl_a01',   p_cl_a01       )
    call iamas_pool_get_dimension(state, 'index_cl_a02',   p_cl_a02       )
    call iamas_pool_get_dimension(state, 'index_cl_a03',   p_cl_a03       )
    call iamas_pool_get_dimension(state, 'index_cl_a04',   p_cl_a04       )
    call iamas_pool_get_dimension(state, 'index_cl_a05',   p_cl_a05       )
    call iamas_pool_get_dimension(state, 'index_cl_a06',   p_cl_a06       )
    call iamas_pool_get_dimension(state, 'index_cl_a07',   p_cl_a07       )
    call iamas_pool_get_dimension(state, 'index_cl_a08',   p_cl_a08       )
    call iamas_pool_get_dimension(state, 'index_cl_a09',   p_cl_a09       )
    call iamas_pool_get_dimension(state, 'index_cl_a10',   p_cl_a10       )

    call iamas_pool_get_dimension(state, 'index_ca_a01',   p_ca_a01       )
    call iamas_pool_get_dimension(state, 'index_ca_a02',   p_ca_a02       )
    call iamas_pool_get_dimension(state, 'index_ca_a03',   p_ca_a03       )
    call iamas_pool_get_dimension(state, 'index_ca_a04',   p_ca_a04       )

    call iamas_pool_get_dimension(state, 'index_co3_a01',  p_co3_a01      )
    call iamas_pool_get_dimension(state, 'index_co3_a02',  p_co3_a02      )
    call iamas_pool_get_dimension(state, 'index_co3_a03',  p_co3_a03      )
    call iamas_pool_get_dimension(state, 'index_co3_a04',  p_co3_a04      )

    call iamas_pool_get_dimension(state, 'index_hysw_a01', p_hysw_a01     )
    call iamas_pool_get_dimension(state, 'index_hysw_a02', p_hysw_a02     )
    call iamas_pool_get_dimension(state, 'index_hysw_a03', p_hysw_a03     )
    call iamas_pool_get_dimension(state, 'index_hysw_a04', p_hysw_a04     )

    call iamas_pool_get_dimension(state, 'index_dust_cw01', p_dust_cw01     )
    call iamas_pool_get_dimension(state, 'index_dust_cw02', p_dust_cw02     )
    call iamas_pool_get_dimension(state, 'index_dust_cw03', p_dust_cw03     )
    call iamas_pool_get_dimension(state, 'index_dust_cw04', p_dust_cw04     )
    call iamas_pool_get_dimension(state, 'index_dust_cw05', p_dust_cw05     )
    call iamas_pool_get_dimension(state, 'index_dust_cw06', p_dust_cw06     )
    call iamas_pool_get_dimension(state, 'index_dust_cw07', p_dust_cw07     )
    call iamas_pool_get_dimension(state, 'index_dust_cw08', p_dust_cw08     )
    call iamas_pool_get_dimension(state, 'index_dust_cw09', p_dust_cw09     )
    call iamas_pool_get_dimension(state, 'index_dust_cw10', p_dust_cw10     )

    call iamas_pool_get_dimension(state, 'index_num_cw01', p_num_cw01     )
    call iamas_pool_get_dimension(state, 'index_num_cw02', p_num_cw02     )
    call iamas_pool_get_dimension(state, 'index_num_cw03', p_num_cw03     )
    call iamas_pool_get_dimension(state, 'index_num_cw04', p_num_cw04     )
    call iamas_pool_get_dimension(state, 'index_num_cw05', p_num_cw05     )
    call iamas_pool_get_dimension(state, 'index_num_cw06', p_num_cw06     )
    call iamas_pool_get_dimension(state, 'index_num_cw07', p_num_cw07     )
    call iamas_pool_get_dimension(state, 'index_num_cw08', p_num_cw08     )
    call iamas_pool_get_dimension(state, 'index_num_cw09', p_num_cw09     )
    call iamas_pool_get_dimension(state, 'index_num_cw10', p_num_cw10     )
 
    call iamas_pool_get_dimension(state, 'index_oc_cw01',   p_oc_cw01       )
    call iamas_pool_get_dimension(state, 'index_oc_cw02',   p_oc_cw02       )
    call iamas_pool_get_dimension(state, 'index_oc_cw03',   p_oc_cw03       )
    call iamas_pool_get_dimension(state, 'index_oc_cw04',   p_oc_cw04       )

    call iamas_pool_get_dimension(state, 'index_bc_cw01',   p_bc_cw01       )
    call iamas_pool_get_dimension(state, 'index_bc_cw02',   p_bc_cw02       )
    call iamas_pool_get_dimension(state, 'index_bc_cw03',   p_bc_cw03       )
    call iamas_pool_get_dimension(state, 'index_bc_cw04',   p_bc_cw04       )

    call iamas_pool_get_dimension(state, 'index_so4_cw01',  p_so4_cw01      )
    call iamas_pool_get_dimension(state, 'index_so4_cw02',  p_so4_cw02      )
    call iamas_pool_get_dimension(state, 'index_so4_cw03',  p_so4_cw03      )
    call iamas_pool_get_dimension(state, 'index_so4_cw04',  p_so4_cw04      )

    call iamas_pool_get_dimension(state, 'index_nh4_cw01',  p_nh4_cw01      )
    call iamas_pool_get_dimension(state, 'index_nh4_cw02',  p_nh4_cw02      )
    call iamas_pool_get_dimension(state, 'index_nh4_cw03',  p_nh4_cw03      )
    call iamas_pool_get_dimension(state, 'index_nh4_cw04',  p_nh4_cw04      )

    call iamas_pool_get_dimension(state, 'index_msa_cw01',  p_msa_cw01      )
    call iamas_pool_get_dimension(state, 'index_msa_cw02',  p_msa_cw02      )
    call iamas_pool_get_dimension(state, 'index_msa_cw03',  p_msa_cw03      )
    call iamas_pool_get_dimension(state, 'index_msa_cw04',  p_msa_cw04      )

    call iamas_pool_get_dimension(state, 'index_no3_cw01',  p_no3_cw01      )
    call iamas_pool_get_dimension(state, 'index_no3_cw02',  p_no3_cw02      )
    call iamas_pool_get_dimension(state, 'index_no3_cw03',  p_no3_cw03      )
    call iamas_pool_get_dimension(state, 'index_no3_cw04',  p_no3_cw04      )

    call iamas_pool_get_dimension(state, 'index_oin_cw01',  p_oin_cw01      )
    call iamas_pool_get_dimension(state, 'index_oin_cw02',  p_oin_cw02      )
    call iamas_pool_get_dimension(state, 'index_oin_cw03',  p_oin_cw03      )
    call iamas_pool_get_dimension(state, 'index_oin_cw04',  p_oin_cw04      )

    call iamas_pool_get_dimension(state, 'index_soa_cw01',  p_soa_cw01      )
    call iamas_pool_get_dimension(state, 'index_soa_cw02',  p_soa_cw02      )
    call iamas_pool_get_dimension(state, 'index_soa_cw03',  p_soa_cw03      )
    call iamas_pool_get_dimension(state, 'index_soa_cw04',  p_soa_cw04      )

    call iamas_pool_get_dimension(state, 'index_na_cw01',   p_na_cw01       )
    call iamas_pool_get_dimension(state, 'index_na_cw02',   p_na_cw02       )
    call iamas_pool_get_dimension(state, 'index_na_cw03',   p_na_cw03       )
    call iamas_pool_get_dimension(state, 'index_na_cw04',   p_na_cw04       )
    call iamas_pool_get_dimension(state, 'index_na_cw05',   p_na_cw05       )
    call iamas_pool_get_dimension(state, 'index_na_cw06',   p_na_cw06       )
    call iamas_pool_get_dimension(state, 'index_na_cw07',   p_na_cw07       )
    call iamas_pool_get_dimension(state, 'index_na_cw08',   p_na_cw08       )
    call iamas_pool_get_dimension(state, 'index_na_cw09',   p_na_cw09       )
    call iamas_pool_get_dimension(state, 'index_na_cw10',   p_na_cw10       )

    call iamas_pool_get_dimension(state, 'index_cl_cw01',   p_cl_cw01       )
    call iamas_pool_get_dimension(state, 'index_cl_cw02',   p_cl_cw02       )
    call iamas_pool_get_dimension(state, 'index_cl_cw03',   p_cl_cw03       )
    call iamas_pool_get_dimension(state, 'index_cl_cw04',   p_cl_cw04       )
    call iamas_pool_get_dimension(state, 'index_cl_cw05',   p_cl_cw05       )
    call iamas_pool_get_dimension(state, 'index_cl_cw06',   p_cl_cw06       )
    call iamas_pool_get_dimension(state, 'index_cl_cw07',   p_cl_cw07       )
    call iamas_pool_get_dimension(state, 'index_cl_cw08',   p_cl_cw08       )
    call iamas_pool_get_dimension(state, 'index_cl_cw09',   p_cl_cw09       )
    call iamas_pool_get_dimension(state, 'index_cl_cw10',   p_cl_cw10       )

    call iamas_pool_get_dimension(state, 'index_ca_cw01',   p_ca_cw01       )
    call iamas_pool_get_dimension(state, 'index_ca_cw02',   p_ca_cw02       )
    call iamas_pool_get_dimension(state, 'index_ca_cw03',   p_ca_cw03       )
    call iamas_pool_get_dimension(state, 'index_ca_cw04',   p_ca_cw04       )

    call iamas_pool_get_dimension(state, 'index_co3_cw01',  p_co3_cw01      )
    call iamas_pool_get_dimension(state, 'index_co3_cw02',  p_co3_cw02      )
    call iamas_pool_get_dimension(state, 'index_co3_cw03',  p_co3_cw03      )
    call iamas_pool_get_dimension(state, 'index_co3_cw04',  p_co3_cw04      )

! zhxia/chenjin add these for photochemistry
    call iamas_pool_get_dimension(state, 'index_co2', p_co2     )
    call iamas_pool_get_dimension(state, 'index_so2', p_so2     )
    call iamas_pool_get_dimension(state, 'index_nh3', p_nh3    )
    call iamas_pool_get_dimension(state, 'index_co', p_co    )
    call iamas_pool_get_dimension(state, 'index_no2', p_no2     )
    call iamas_pool_get_dimension(state, 'index_no', p_no    )
    call iamas_pool_get_dimension(state, 'index_o3', p_o3    )
    call iamas_pool_get_dimension(state, 'index_h2o2', p_h2o2    )
    call iamas_pool_get_dimension(state, 'index_no3', p_no3    )
    call iamas_pool_get_dimension(state, 'index_n2o5', p_n2o5    )
    call iamas_pool_get_dimension(state, 'index_hono', p_hono    )
    call iamas_pool_get_dimension(state, 'index_hno3', p_hno3    )
    call iamas_pool_get_dimension(state, 'index_hno4', p_hno4    )
    call iamas_pool_get_dimension(state, 'index_h2so4', p_h2so4    )
    call iamas_pool_get_dimension(state, 'index_h2so4', p_sulf    )
    call iamas_pool_get_dimension(state, 'index_hcho', p_hcho    )
    call iamas_pool_get_dimension(state, 'index_ccho', p_ccho    )
    call iamas_pool_get_dimension(state, 'index_rcho', p_rcho    )
    call iamas_pool_get_dimension(state, 'index_acet', p_acet    )
    call iamas_pool_get_dimension(state, 'index_mek', p_mek    )
    call iamas_pool_get_dimension(state, 'index_hcooh', p_hcooh    )
    call iamas_pool_get_dimension(state, 'index_meoh', p_meoh    )
    call iamas_pool_get_dimension(state, 'index_etoh', p_etoh    )
    call iamas_pool_get_dimension(state, 'index_cco_oh', p_cco_oh    )
    call iamas_pool_get_dimension(state, 'index_rco_oh', p_rco_oh    )
    call iamas_pool_get_dimension(state, 'index_gly', p_gly    )
    call iamas_pool_get_dimension(state, 'index_mgly', p_mgly    )
    call iamas_pool_get_dimension(state, 'index_bacl', p_bacl    )
    call iamas_pool_get_dimension(state, 'index_cres', p_cres    )
    call iamas_pool_get_dimension(state, 'index_bald', p_bald    )
    call iamas_pool_get_dimension(state, 'index_isoprod', p_isoprod    )
    call iamas_pool_get_dimension(state, 'index_methacro', p_methacro    )
    call iamas_pool_get_dimension(state, 'index_mvk', p_mvk    )
    call iamas_pool_get_dimension(state, 'index_prod2', p_prod2    )
    call iamas_pool_get_dimension(state, 'index_dcb1', p_dcb1    )
    call iamas_pool_get_dimension(state, 'index_dcb2', p_dcb2    )
    call iamas_pool_get_dimension(state, 'index_dcb3', p_dcb3    )
    call iamas_pool_get_dimension(state, 'index_ethene', p_ethene    )
    call iamas_pool_get_dimension(state, 'index_isoprene', p_isoprene    )
    call iamas_pool_get_dimension(state, 'index_c2h6', p_c2h6    )
    call iamas_pool_get_dimension(state, 'index_c3h8', p_c3h8    )
    call iamas_pool_get_dimension(state, 'index_c2h2', p_c2h2    )
    call iamas_pool_get_dimension(state, 'index_c3h6', p_c3h6    )
    call iamas_pool_get_dimension(state, 'index_alk3', p_alk3    )
    call iamas_pool_get_dimension(state, 'index_alk4', p_alk4    )
    call iamas_pool_get_dimension(state, 'index_alk5', p_alk5    )
    call iamas_pool_get_dimension(state, 'index_aro1', p_aro1    )
    call iamas_pool_get_dimension(state, 'index_aro2', p_aro2    )
    call iamas_pool_get_dimension(state, 'index_ole1', p_ole1    )
    call iamas_pool_get_dimension(state, 'index_ole2', p_ole2    )
    call iamas_pool_get_dimension(state, 'index_terp', p_terp    )
    call iamas_pool_get_dimension(state, 'index_rno3', p_rno3    )
    call iamas_pool_get_dimension(state, 'index_nphe', p_nphe    )
    call iamas_pool_get_dimension(state, 'index_phen', p_phen    )
    call iamas_pool_get_dimension(state, 'index_pan', p_pan    )
    call iamas_pool_get_dimension(state, 'index_pan2', p_pan2    )
    call iamas_pool_get_dimension(state, 'index_pbzn', p_pbzn    )
    call iamas_pool_get_dimension(state, 'index_ma_pan', p_ma_pan    )
    call iamas_pool_get_dimension(state, 'index_cco_ooh', p_cco_ooh    )
    call iamas_pool_get_dimension(state, 'index_rco_o2', p_rco_o2    )
    call iamas_pool_get_dimension(state, 'index_rco_ooh', p_rco_ooh    )
    call iamas_pool_get_dimension(state, 'index_xn', p_xn    )
    call iamas_pool_get_dimension(state, 'index_xc', p_xc    )
    call iamas_pool_get_dimension(state, 'index_ho', p_ho    )
    call iamas_pool_get_dimension(state, 'index_ho2', p_ho2    )
    call iamas_pool_get_dimension(state, 'index_c_o2', p_c_o2    )
    call iamas_pool_get_dimension(state, 'index_cooh', p_cooh    )
    call iamas_pool_get_dimension(state, 'index_rooh', p_rooh    )
    call iamas_pool_get_dimension(state, 'index_ro2_r', p_ro2_r    )
    call iamas_pool_get_dimension(state, 'index_r2o2', p_r2o2    )
    call iamas_pool_get_dimension(state, 'index_ro2_n', p_ro2_n    )
    call iamas_pool_get_dimension(state, 'index_cco_o2', p_cco_o2    )
    call iamas_pool_get_dimension(state, 'index_bzco_o2', p_bzco_o2    )
    call iamas_pool_get_dimension(state, 'index_ma_rco3', p_ma_rco3    )
    call iamas_pool_get_dimension(state, 'index_sesq', p_sesq    )
    call iamas_pool_get_dimension(state, 'index_ch4', p_ch4    )
    call iamas_pool_get_dimension(state, 'index_hcl', p_hcl    )
    call iamas_pool_get_dimension(state, 'index_msa', p_msa    )

    call iamas_pool_get_dimension(state, 'index_alk1', p_alk1    )
    call iamas_pool_get_dimension(state, 'index_alk2', p_alk2    )

! zhxia add these for photolysis
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_o2', p_ph_o2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_o31d', p_ph_o31d)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_o33p', p_ph_o33p)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_no2', p_ph_no2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_no3o2', p_ph_no3o2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_no3o', p_ph_no3o)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hno2', p_ph_hno2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hno3', p_ph_hno3)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hno4', p_ph_hno4)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_h2o2', p_ph_h2o2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch2or', p_ph_ch2or)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch2om', p_ph_ch2om)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3cho', p_ph_ch3cho)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3coch3', p_ph_ch3coch3)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3coc2h5',p_ph_ch3coc2h5)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hcocho', p_ph_hcocho)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3cocho', p_ph_ch3cocho)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hcochest', p_ph_hcochest)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3o2h', p_ph_ch3o2h)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3coo2h', p_ph_ch3coo2h)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_ch3ono2', p_ph_ch3ono2)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_hcochob', p_ph_hcochob)
    call iamas_pool_get_dimension(tend_chemistry,'index_ph_n2o5', p_ph_n2o5)

end subroutine init_data_binaer_index


!-----------------------------------------------------------------------
    subroutine init_data_binaer_ptr( is_aerosol )

    ! use module_configure
    use iamas_atmchem_vars, only:  param_first_scalar,num_chem  , &
        p_qv,        p_qc,       p_qr,    p_qi,        p_qs    , &
        p_qndrop,    p_qniceaer,  p_qnicemix                   , &
        p_dust_a01,  p_dust_a02,  p_dust_a03,  p_dust_a04,  p_dust_a05   , &
        p_dust_a06,  p_dust_a07,  p_dust_a08,  p_dust_a09,  p_dust_a10   , &
        p_water_a01, p_water_a02, p_water_a03, p_water_a04, p_water_a05  , &
        p_water_a06, p_water_a07, p_water_a08, p_water_a09, p_water_a10  , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05    , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10    , &
        p_dust_cw01, p_dust_cw02, p_dust_cw03, p_dust_cw04, p_dust_cw05  , &
        p_dust_cw06, p_dust_cw07, p_dust_cw08, p_dust_cw09, p_dust_cw10  , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05   , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10   , &
        p_oc_a01,    p_oc_a02,    p_oc_a03,    p_oc_a04                  , &
        p_bc_a01,    p_bc_a02,    p_bc_a03,    p_bc_a04                  , &
        p_so4_a01,   p_so4_a02,   p_so4_a03,   p_so4_a04                 , &
        p_nh4_a01,   p_nh4_a02,   p_nh4_a03,   p_nh4_a04                 , &
        p_msa_a01,   p_msa_a02,   p_msa_a03,   p_msa_a04                 , &
        p_no3_a01,   p_no3_a02,   p_no3_a03,   p_no3_a04                 , &
        p_oin_a01,   p_oin_a02,   p_oin_a03,   p_oin_a04                 , &
        p_soa_a01,   p_soa_a02,   p_soa_a03,   p_soa_a04                 , &
        p_na_a01,    p_na_a02,    p_na_a03,    p_na_a04                  , &
        p_na_a05,    p_na_a06,    p_na_a07,    p_na_a08,    p_na_a09,    p_na_a10, &
        p_cl_a01,    p_cl_a02,    p_cl_a03,    p_cl_a04                  , &
        p_cl_a05,    p_cl_a06,    p_cl_a07,    p_cl_a08,    p_cl_a09,    p_cl_a10, &
        p_ca_a01,    p_ca_a02,    p_ca_a03,    p_ca_a04                  , &
        p_co3_a01,   p_co3_a02,   p_co3_a03,   p_co3_a04                 , &
        p_hysw_a01,  p_hysw_a02,  p_hysw_a03,  p_hysw_a04                , &
        p_oc_cw01,    p_oc_cw02,    p_oc_cw03,    p_oc_cw04              , &
        p_bc_cw01,    p_bc_cw02,    p_bc_cw03,    p_bc_cw04              , &
        p_so4_cw01,   p_so4_cw02,   p_so4_cw03,   p_so4_cw04             , &
        p_nh4_cw01,   p_nh4_cw02,   p_nh4_cw03,   p_nh4_cw04             , &
        p_msa_cw01,   p_msa_cw02,   p_msa_cw03,   p_msa_cw04             , &
        p_no3_cw01,   p_no3_cw02,   p_no3_cw03,   p_no3_cw04             , &
        p_oin_cw01,   p_oin_cw02,   p_oin_cw03,   p_oin_cw04             , &
        p_soa_cw01,   p_soa_cw02,   p_soa_cw03,   p_soa_cw04             , &
        p_na_cw01,    p_na_cw02,    p_na_cw03,    p_na_cw04              , &
        p_na_cw05,    p_na_cw06,    p_na_cw07,    p_na_cw08,    p_na_cw09,    p_na_cw10, &
        p_cl_cw01,    p_cl_cw02,    p_cl_cw03,    p_cl_cw04              , &
        p_cl_cw05,    p_cl_cw06,    p_cl_cw07,    p_cl_cw08,    p_cl_cw09,    p_cl_cw10, &
        p_ca_cw01,    p_ca_cw02,    p_ca_cw03,    p_ca_cw04              , &
        p_co3_cw01,   p_co3_cw02,   p_co3_cw03,   p_co3_cw04             , &
        p_co2,       p_no,        p_no2,         p_so2                   , &
        p_nh3,       p_co,        p_o3,          p_h2o2                  , &
        p_no3,       p_n2o5,      p_hono,        p_hno3                  , &
        p_hno4,      p_h2so4,     p_hcho,        p_ccho                  , &
        p_rcho,      p_acet,      p_mek,         p_hcooh                 , &
        p_meoh,      p_etoh,      p_cco_oh,      p_rco_oh                , &
        p_gly,       p_mgly,      p_bacl,        p_cres                  , &
        p_bald,      p_isoprod,   p_methacro,    p_mvk                   , &
        p_prod2,     p_dcb1,      p_dcb2,        p_dcb3                  , &
        p_ethene,    p_isoprene,  p_c2h6,        p_c3h8                  , &
        p_c2h2,      p_c3h6,      p_alk3,        p_alk4                  , &
        p_alk5,      p_aro1,      p_aro2,        p_ole1                  , &
        p_ole2,      p_terp,      p_rno3,        p_nphe                  , &
        p_phen,      p_pan,       p_pan2,        p_pbzn                  , &
        p_ma_pan,    p_cco_ooh,   p_rco_o2,      p_rco_ooh               , &
        p_xn,        p_xc,        p_ho,          p_ho2                   , &
        p_c_o2,      p_cooh,      p_rooh,        p_ro2_r                 , &
        p_r2o2,      p_ro2_n,     p_cco_o2,      p_bzco_o2               , &
        p_ma_rco3,   p_sesq,      p_ch4,         p_hcl                   , &
        p_msa,       p_sulf                                              , &
        p_alk1,      p_alk2



    use module_data_binaer_asect
    use module_data_binaer_other, only:   &
        kh2so4, khno3, khcl, knh3,kmsa, ko3, kh2o, ktemp,   &
        kso2, kh2o2, khcho, koh, kho2,   &
        kno3, kno, kno2, khono, kpan,   &
                lmaxd, l2maxd, ltot, ltot2, lunout, lunerr, &
                name,kpcg1_b_c,kpcg2_b_c,kpcg3_b_c,kpcg4_b_c, &
                kpcg5_b_c,kpcg6_b_c,kpcg7_b_c,kpcg8_b_c, &
                kpcg9_b_c,kpcg1_b_o,kpcg2_b_o,kpcg3_b_o, &
                kpcg4_b_o,kpcg5_b_o,kpcg6_b_o,kpcg7_b_o, &
                kpcg8_b_o,kpcg9_b_o,kopcg1_b_c,kopcg2_b_c,&
                kopcg3_b_c, kopcg4_b_c,kopcg5_b_c,kopcg6_b_c,&
                kopcg7_b_c,kopcg8_b_c,kopcg1_b_o,kopcg2_b_o,&
                kopcg3_b_o,kopcg4_b_o,kopcg5_b_o,kopcg6_b_o,&
                kopcg7_b_o,kopcg8_b_o,&
                kpcg1_f_c,kpcg2_f_c,kpcg3_f_c,kpcg4_f_c, &
                kpcg5_f_c,kpcg6_f_c,kpcg7_f_c,kpcg8_f_c, &
                kpcg9_f_c,kpcg1_f_o,kpcg2_f_o,kpcg3_f_o, &
                kpcg4_f_o,kpcg5_f_o,kpcg6_f_o,kpcg7_f_o, &
                kpcg8_f_o,kpcg9_f_o,kopcg1_f_c,kopcg2_f_c,&
                kopcg3_f_c, kopcg4_f_c,kopcg5_f_c,kopcg6_f_c,&
                kopcg7_f_c,kopcg8_f_c,kopcg1_f_o,kopcg2_f_o,&
                kopcg3_f_o,kopcg4_f_o,kopcg5_f_o,kopcg6_f_o,&
                kopcg7_f_o,kopcg8_f_o,                      &
                ksmpa,ksmpbb,                      &
                kant1_c,kant2_c,kant3_c,kant4_c,kant1_o,kant2_o, &
                kant3_o,kant4_o,                            &
                kbiog1_c,kbiog2_c,kbiog3_c,kbiog4_c,kbiog1_o,kbiog2_o, &
                kbiog3_o,kbiog4_o, &
                kn2o5, kclno2                            

    ! use module_peg_util, only:  peg_error_fatal, peg_message
    use iamas_atmchem_utilities, only: chemistry_error_fatal, chemistry_message
    use module_mosaic_wetscav, only: initwet
        ! use module_scalar_tables, only:  chem_dname_table

    implicit none

!   subr arguments
        logical, intent(out) :: is_aerosol(num_chem)
!   local variables
    integer l, ll, n, p1st
    integer iaddto_ncomp, iaddto_ncomp_plustracer
    integer l_mastercomp, lptr_dum
    integer mcindx_dum
    integer isize, itype, iphase
    integer nphasetxt, nsizetxt, nspectxt, ntypetxt
    integer ncomp_dum(maxd_asize,maxd_aphase)
    integer ncomp_plustracer_dum(maxd_asize,maxd_aphase)

    integer y_oc, y_bc, y_so4, y_no3, y_oin, y_dust, y_nh4, y_msa, y_soa, &
            y_cl, y_na, y_ca, y_co3, &
            y_hysw, y_water, y_num
    integer y_cw_dust, y_cw_num

    character*200 msg
    character*20 phasetxt, sizetxt, spectxt, typetxt

    p1st = param_first_scalar
!
!   set up pointers to aerosol species in the wrf-chem "chem" array
!   note:  lptr=1 points to the first chem species which is "unused"
!
    itype=1
    lptr_bc_aer(:,itype, :)       = -1
    lptr_oc_aer(:,itype,:)        = -1
    lptr_so4_aer(:,itype, :)      = -1
    lptr_msa_aer(:,itype, :)      = -1
    lptr_no3_aer(:,itype, :)      = -1
    lptr_cl_aer(:,itype, :)       = -1
    lptr_co3_aer(:,itype, :)      = -1
    lptr_nh4_aer(:,itype, :)      = -1
    lptr_na_aer(:,itype, :)       = -1
    lptr_ca_aer(:,itype, :)       = -1
    lptr_oin_aer(:,itype, :)      = -1
    lptr_dust_aer(:,itype, :)     = -1
    lptr_soa_aer(:,itype, :)     = -1
    lptr_tr1r1_aer(:,itype, :)    = -1
    lptr_tr1r2_aer(:,itype, :)    = -1
    lptr_tr1r3_aer(:,itype, :)    = -1
    lptr_tr1r4_aer(:,itype, :)    = -1
    lptr_pcg1_b_c_aer(:,itype, :) = -1
    lptr_pcg2_b_c_aer(:,itype, :) = -1
    lptr_pcg3_b_c_aer(:,itype, :) = -1
    lptr_pcg4_b_c_aer(:,itype, :) = -1
    lptr_pcg5_b_c_aer(:,itype, :) = -1
    lptr_pcg6_b_c_aer(:,itype, :) = -1
    lptr_pcg7_b_c_aer(:,itype, :) = -1
    lptr_pcg8_b_c_aer(:,itype, :) = -1
    lptr_pcg9_b_c_aer(:,itype, :) = -1
    lptr_pcg1_b_o_aer(:,itype, :) = -1
    lptr_pcg2_b_o_aer(:,itype, :) = -1
    lptr_pcg3_b_o_aer(:,itype, :) = -1
    lptr_pcg4_b_o_aer(:,itype, :) = -1
    lptr_pcg5_b_o_aer(:,itype, :) = -1
    lptr_pcg6_b_o_aer(:,itype, :) = -1
    lptr_pcg7_b_o_aer(:,itype, :) = -1
    lptr_pcg8_b_o_aer(:,itype, :) = -1
    lptr_pcg9_b_o_aer(:,itype, :) = -1
    lptr_opcg1_b_c_aer(:,itype, :)= -1
    lptr_opcg2_b_c_aer(:,itype, :)= -1
    lptr_opcg3_b_c_aer(:,itype, :)= -1
    lptr_opcg4_b_c_aer(:,itype, :)= -1
    lptr_opcg5_b_c_aer(:,itype, :)= -1
    lptr_opcg6_b_c_aer(:,itype, :)= -1
    lptr_opcg7_b_c_aer(:,itype, :)= -1
    lptr_opcg8_b_c_aer(:,itype, :)= -1
    lptr_opcg1_b_o_aer(:,itype, :)= -1
    lptr_opcg2_b_o_aer(:,itype, :)= -1
    lptr_opcg3_b_o_aer(:,itype, :)= -1
    lptr_opcg4_b_o_aer(:,itype, :)= -1
    lptr_opcg5_b_o_aer(:,itype, :)= -1
    lptr_opcg6_b_o_aer(:,itype, :)= -1
    lptr_opcg7_b_o_aer(:,itype, :)= -1
    lptr_opcg8_b_o_aer(:,itype, :)= -1
    lptr_pcg1_f_c_aer(:,itype, :) = -1
    lptr_pcg2_f_c_aer(:,itype, :) = -1
    lptr_pcg3_f_c_aer(:,itype, :) = -1
    lptr_pcg4_f_c_aer(:,itype, :) = -1
    lptr_pcg5_f_c_aer(:,itype, :) = -1
    lptr_pcg6_f_c_aer(:,itype, :) = -1
    lptr_pcg7_f_c_aer(:,itype, :) = -1
    lptr_pcg8_f_c_aer(:,itype, :) = -1
    lptr_pcg9_f_c_aer(:,itype, :) = -1
    lptr_pcg1_f_o_aer(:,itype, :) = -1
    lptr_pcg2_f_o_aer(:,itype, :) = -1
    lptr_pcg3_f_o_aer(:,itype, :) = -1
    lptr_pcg4_f_o_aer(:,itype, :) = -1
    lptr_pcg5_f_o_aer(:,itype, :) = -1
    lptr_pcg6_f_o_aer(:,itype, :) = -1
    lptr_pcg7_f_o_aer(:,itype, :) = -1
    lptr_pcg8_f_o_aer(:,itype, :) = -1
    lptr_pcg9_f_o_aer(:,itype, :) = -1
    lptr_opcg1_f_c_aer(:,itype, :)= -1
    lptr_opcg2_f_c_aer(:,itype, :)= -1
    lptr_opcg3_f_c_aer(:,itype, :)= -1
    lptr_opcg4_f_c_aer(:,itype, :)= -1
    lptr_opcg5_f_c_aer(:,itype, :)= -1
    lptr_opcg6_f_c_aer(:,itype, :)= -1
    lptr_opcg7_f_c_aer(:,itype, :)= -1
    lptr_opcg8_f_c_aer(:,itype, :)= -1
    lptr_opcg1_f_o_aer(:,itype, :)= -1
    lptr_opcg2_f_o_aer(:,itype, :)= -1
    lptr_opcg3_f_o_aer(:,itype, :)= -1
    lptr_opcg4_f_o_aer(:,itype, :)= -1
    lptr_opcg5_f_o_aer(:,itype, :)= -1
    lptr_opcg6_f_o_aer(:,itype, :)= -1
    lptr_opcg7_f_o_aer(:,itype, :)= -1
    lptr_opcg8_f_o_aer(:,itype, :)= -1
    lptr_smpa_aer(:,itype, :)     = -1
    lptr_smpbb_aer(:,itype, :)    = -1
    lptr_ant1_c_aer(:,itype, :)   = -1
    lptr_ant2_c_aer(:,itype, :)   = -1
    lptr_ant3_c_aer(:,itype, :)   = -1
    lptr_ant4_c_aer(:,itype, :)   = -1
    lptr_ant1_o_aer(:,itype, :)   = -1
    lptr_ant2_o_aer(:,itype, :)   = -1
    lptr_ant3_o_aer(:,itype, :)   = -1
    lptr_ant4_o_aer(:,itype, :)   = -1
    lptr_biog1_c_aer(:,itype, :)  = -1
    lptr_biog2_c_aer(:,itype, :)  = -1
    lptr_biog3_c_aer(:,itype, :)  = -1
    lptr_biog4_c_aer(:,itype, :)  = -1
    lptr_biog1_o_aer(:,itype, :)  = -1
    lptr_biog2_o_aer(:,itype, :)  = -1
    lptr_biog3_o_aer(:,itype, :)  = -1
    lptr_biog4_o_aer(:,itype, :)  = -1


    waterptr_aer(:,itype)         = -1
    hyswptr_aer(:,itype)          = -1
    numptr_aer(:,itype,:)         = -1
 

    if (nsize_aer(itype) .ge. 1) then
        if (p_oc_a01    .ge. p1st) lptr_oc_aer(01,itype,ai_phase)   = p_oc_a01
        if (p_bc_a01    .ge. p1st) lptr_bc_aer(01,itype,ai_phase)   = p_bc_a01
        if (p_so4_a01   .ge. p1st) lptr_so4_aer(01,itype,ai_phase)  = p_so4_a01
        if (p_nh4_a01   .ge. p1st) lptr_nh4_aer(01,itype,ai_phase)  = p_nh4_a01
        if (p_msa_a01   .ge. p1st) lptr_msa_aer(01,itype,ai_phase)  = p_msa_a01
        if (p_no3_a01   .ge. p1st) lptr_no3_aer(01,itype,ai_phase)  = p_no3_a01
        if (p_oin_a01   .ge. p1st) lptr_oin_aer(01,itype,ai_phase)  = p_oin_a01
        if (p_cl_a01    .ge. p1st) lptr_cl_aer(01,itype,ai_phase)   = p_cl_a01
        if (p_na_a01    .ge. p1st) lptr_na_aer(01,itype,ai_phase)   = p_na_a01
        if (p_ca_a01    .ge. p1st) lptr_ca_aer(01,itype,ai_phase)   = p_ca_a01
        if (p_co3_a01   .ge. p1st) lptr_co3_aer(01,itype,ai_phase)  = p_co3_a01
        if (p_dust_a01  .ge. p1st) lptr_dust_aer(01,itype,ai_phase) = p_dust_a01
        if (p_soa_a01   .ge. p1st) lptr_soa_aer(01,itype,ai_phase)  = p_soa_a01
        if (p_hysw_a01  .ge. p1st) hyswptr_aer(01,itype)            = p_hysw_a01
        if (p_water_a01 .ge. p1st) waterptr_aer(01,itype)           = p_water_a01
        if (p_num_a01   .ge. p1st) numptr_aer(01,itype,ai_phase)    = p_num_a01
    end if

    if (nsize_aer(itype) .ge. 2) then
        if (p_oc_a02    .ge. p1st) lptr_oc_aer(02,itype,ai_phase)   = p_oc_a02
        if (p_bc_a02    .ge. p1st) lptr_bc_aer(02,itype,ai_phase)   = p_bc_a02
        if (p_so4_a02   .ge. p1st) lptr_so4_aer(02,itype,ai_phase)  = p_so4_a02
        if (p_nh4_a02   .ge. p1st) lptr_nh4_aer(02,itype,ai_phase)  = p_nh4_a02
        if (p_msa_a02   .ge. p1st) lptr_msa_aer(02,itype,ai_phase)  = p_msa_a02
        if (p_no3_a02   .ge. p1st) lptr_no3_aer(02,itype,ai_phase)  = p_no3_a02
        if (p_oin_a02   .ge. p1st) lptr_oin_aer(02,itype,ai_phase)  = p_oin_a02
        if (p_cl_a02    .ge. p1st) lptr_cl_aer(02,itype,ai_phase)   = p_cl_a02
        if (p_na_a02    .ge. p1st) lptr_na_aer(02,itype,ai_phase)   = p_na_a02
        if (p_ca_a02    .ge. p1st) lptr_ca_aer(02,itype,ai_phase)   = p_ca_a02
        if (p_co3_a02   .ge. p1st) lptr_co3_aer(02,itype,ai_phase)  = p_co3_a02
        if (p_dust_a02  .ge. p1st) lptr_dust_aer(02,itype,ai_phase) = p_dust_a02
        if (p_soa_a02   .ge. p1st) lptr_soa_aer(02,itype,ai_phase)  = p_soa_a02
        if (p_hysw_a02  .ge. p1st) hyswptr_aer(02,itype)            = p_hysw_a02
        if (p_water_a02 .ge. p1st) waterptr_aer(02,itype)           = p_water_a02
        if (p_num_a02   .ge. p1st) numptr_aer(02,itype,ai_phase)    = p_num_a02
    end if

    if (nsize_aer(itype) .ge. 3) then
        if (p_oc_a03    .ge. p1st) lptr_oc_aer(03,itype,ai_phase)   = p_oc_a03
        if (p_bc_a03    .ge. p1st) lptr_bc_aer(03,itype,ai_phase)   = p_bc_a03
        if (p_so4_a03   .ge. p1st) lptr_so4_aer(03,itype,ai_phase)  = p_so4_a03
        if (p_nh4_a03   .ge. p1st) lptr_nh4_aer(03,itype,ai_phase)  = p_nh4_a03
        if (p_msa_a03   .ge. p1st) lptr_msa_aer(03,itype,ai_phase)  = p_msa_a03
        if (p_no3_a03   .ge. p1st) lptr_no3_aer(03,itype,ai_phase)  = p_no3_a03
        if (p_oin_a03   .ge. p1st) lptr_oin_aer(03,itype,ai_phase)  = p_oin_a03
        if (p_cl_a03    .ge. p1st) lptr_cl_aer(03,itype,ai_phase)   = p_cl_a03
        if (p_na_a03    .ge. p1st) lptr_na_aer(03,itype,ai_phase)   = p_na_a03
        if (p_ca_a03    .ge. p1st) lptr_ca_aer(03,itype,ai_phase)   = p_ca_a03
        if (p_co3_a03   .ge. p1st) lptr_co3_aer(03,itype,ai_phase)  = p_co3_a03
        if (p_dust_a03  .ge. p1st) lptr_dust_aer(03,itype,ai_phase) = p_dust_a03
        if (p_soa_a03   .ge. p1st) lptr_soa_aer(03,itype,ai_phase)  = p_soa_a03
        if (p_hysw_a03  .ge. p1st) hyswptr_aer(03,itype)            = p_hysw_a03
        if (p_water_a03 .ge. p1st) waterptr_aer(03,itype)           = p_water_a03
        if (p_num_a03   .ge. p1st) numptr_aer(03,itype,ai_phase)    = p_num_a03
    end if

    if (nsize_aer(itype) .ge. 4) then
        if (p_oc_a04    .ge. p1st) lptr_oc_aer(04,itype,ai_phase)   = p_oc_a04
        if (p_bc_a04    .ge. p1st) lptr_bc_aer(04,itype,ai_phase)   = p_bc_a04
        if (p_so4_a04   .ge. p1st) lptr_so4_aer(04,itype,ai_phase)  = p_so4_a04
        if (p_nh4_a04   .ge. p1st) lptr_nh4_aer(04,itype,ai_phase)  = p_nh4_a04
        if (p_msa_a04   .ge. p1st) lptr_msa_aer(04,itype,ai_phase)  = p_msa_a04
        if (p_no3_a04   .ge. p1st) lptr_no3_aer(04,itype,ai_phase)  = p_no3_a04
        if (p_oin_a04   .ge. p1st) lptr_oin_aer(04,itype,ai_phase)  = p_oin_a04
        if (p_cl_a04    .ge. p1st) lptr_cl_aer(04,itype,ai_phase)   = p_cl_a04
        if (p_na_a04    .ge. p1st) lptr_na_aer(04,itype,ai_phase)   = p_na_a04
        if (p_ca_a04    .ge. p1st) lptr_ca_aer(04,itype,ai_phase)   = p_ca_a04
        if (p_co3_a04   .ge. p1st) lptr_co3_aer(04,itype,ai_phase)  = p_co3_a04
        if (p_dust_a04  .ge. p1st) lptr_dust_aer(04,itype,ai_phase) = p_dust_a04
        if (p_soa_a04   .ge. p1st) lptr_soa_aer(04,itype,ai_phase)  = p_soa_a04
        if (p_hysw_a04  .ge. p1st) hyswptr_aer(04,itype)            = p_hysw_a04
        if (p_water_a04 .ge. p1st) waterptr_aer(04,itype)           = p_water_a04
        if (p_num_a04   .ge. p1st) numptr_aer(04,itype,ai_phase)    = p_num_a04
    end if

!jwfeng:
!I don't know why wrfchem has no "if" begore hysw,water and num pointer in bin 5,6,7,8
!So I just copy that.        
    if (nsize_aer(itype) .ge. 5) then
        if (p_cl_a05  .ge. p1st) lptr_cl_aer(05,itype,ai_phase) = p_cl_a05
        if (p_na_a05  .ge. p1st) lptr_na_aer(05,itype,ai_phase) = p_na_a05
        if (p_dust_a05  .ge. p1st) lptr_dust_aer(05,itype,ai_phase) = p_dust_a05
                                   waterptr_aer(05,itype)           = p_water_a05
                                   numptr_aer(05,itype,ai_phase)    = p_num_a05
    end if

    if (nsize_aer(itype) .ge. 6) then
        if (p_cl_a06  .ge. p1st) lptr_cl_aer(06,itype,ai_phase) = p_cl_a06
        if (p_na_a06  .ge. p1st) lptr_na_aer(06,itype,ai_phase) = p_na_a06
        if (p_dust_a06  .ge. p1st) lptr_dust_aer(06,itype,ai_phase) = p_dust_a06
                                   waterptr_aer(06,itype)           = p_water_a06
                                   numptr_aer(06,itype,ai_phase)    = p_num_a06
    end if

    if (nsize_aer(itype) .ge. 7) then
        if (p_cl_a07  .ge. p1st) lptr_cl_aer(07,itype,ai_phase) = p_cl_a07
        if (p_na_a07  .ge. p1st) lptr_na_aer(07,itype,ai_phase) = p_na_a07
        if (p_dust_a07  .ge. p1st) lptr_dust_aer(07,itype,ai_phase) = p_dust_a07
                                   waterptr_aer(07,itype)           = p_water_a07
                                   numptr_aer(07,itype,ai_phase)    = p_num_a07
    end if

    if (nsize_aer(itype) .ge. 8) then
        if (p_cl_a08  .ge. p1st) lptr_cl_aer(08,itype,ai_phase) = p_cl_a08
        if (p_na_a08  .ge. p1st) lptr_na_aer(08,itype,ai_phase) = p_na_a08
        if (p_dust_a08  .ge. p1st) lptr_dust_aer(08,itype,ai_phase) = p_dust_a08
                                   waterptr_aer(08,itype)           = p_water_a08
                                   numptr_aer(08,itype,ai_phase)    = p_num_a08
    end if

    if (nsize_aer(itype) .ge. 9) then
        if (p_cl_a09  .ge. p1st) lptr_cl_aer(09,itype,ai_phase) = p_cl_a09
        if (p_na_a09  .ge. p1st) lptr_na_aer(09,itype,ai_phase) = p_na_a09
        if (p_dust_a09  .ge. p1st) lptr_dust_aer(09,itype,ai_phase) = p_dust_a09
                                   waterptr_aer(09,itype)           = p_water_a09
                                   numptr_aer(09,itype,ai_phase)    = p_num_a09
    end if

    if (nsize_aer(itype) .ge. 10) then
        if (p_cl_a10  .ge. p1st) lptr_cl_aer(10,itype,ai_phase) = p_cl_a10
        if (p_na_a10  .ge. p1st) lptr_na_aer(10,itype,ai_phase) = p_na_a10
        if (p_dust_a10  .ge. p1st) lptr_dust_aer(10,itype,ai_phase) = p_dust_a10
                                   waterptr_aer(10,itype)           = p_water_a10
                                   numptr_aer(10,itype,ai_phase)    = p_num_a10
    end if

if (cw_phase .gt. 0) then
!   this code is "active" only when cw species are in the registry
    if (nsize_aer(itype) .ge. 1) then
        if(p_oc_cw01  .ge.p1st) lptr_oc_aer(01,itype,cw_phase)   = p_oc_cw01
        if(p_bc_cw01  .ge.p1st) lptr_bc_aer(01,itype,cw_phase)   = p_bc_cw01
        if(p_so4_cw01  .ge.p1st) lptr_so4_aer(01,itype,cw_phase)   = p_so4_cw01
        if(p_nh4_cw01  .ge.p1st) lptr_nh4_aer(01,itype,cw_phase)   = p_nh4_cw01
        if(p_msa_cw01  .ge.p1st) lptr_msa_aer(01,itype,cw_phase)   = p_msa_cw01
        if(p_no3_cw01  .ge.p1st) lptr_no3_aer(01,itype,cw_phase)   = p_no3_cw01
        if(p_oin_cw01  .ge.p1st) lptr_oin_aer(01,itype,cw_phase)   = p_oin_cw01
        if(p_cl_cw01   .ge.p1st) lptr_cl_aer(01,itype,cw_phase)    = p_cl_cw01
        if(p_na_cw01   .ge.p1st) lptr_na_aer(01,itype,cw_phase)    = p_na_cw01
        if(p_ca_cw01   .ge.p1st) lptr_ca_aer(01,itype,cw_phase)    = p_ca_cw01
        if(p_co3_cw01  .ge.p1st) lptr_co3_aer(01,itype,cw_phase)   = p_co3_cw01
        if(p_dust_cw01.ge.p1st) lptr_dust_aer(01,itype,cw_phase) = p_dust_cw01
        if(p_soa_cw01  .ge.p1st) lptr_soa_aer(01,itype,cw_phase)   = p_soa_cw01
        numptr_aer(01,itype,cw_phase)        = p_num_cw01
    end if

    if (nsize_aer(itype) .ge. 2) then
        if(p_oc_cw02  .ge.p1st) lptr_oc_aer(02,itype,cw_phase)   = p_oc_cw02
        if(p_bc_cw02  .ge.p1st) lptr_bc_aer(02,itype,cw_phase)   = p_bc_cw02
        if(p_so4_cw02  .ge.p1st) lptr_so4_aer(02,itype,cw_phase)   = p_so4_cw02
        if(p_nh4_cw02  .ge.p1st) lptr_nh4_aer(02,itype,cw_phase)   = p_nh4_cw02
        if(p_msa_cw02  .ge.p1st) lptr_msa_aer(02,itype,cw_phase)   = p_msa_cw02
        if(p_no3_cw02  .ge.p1st) lptr_no3_aer(02,itype,cw_phase)   = p_no3_cw02
        if(p_oin_cw02  .ge.p1st) lptr_oin_aer(02,itype,cw_phase)   = p_oin_cw02
        if(p_cl_cw02   .ge.p1st) lptr_cl_aer(02,itype,cw_phase)    = p_cl_cw02
        if(p_na_cw02   .ge.p1st) lptr_na_aer(02,itype,cw_phase)    = p_na_cw02
        if(p_ca_cw02   .ge.p1st) lptr_ca_aer(02,itype,cw_phase)    = p_ca_cw02
        if(p_co3_cw02  .ge.p1st) lptr_co3_aer(02,itype,cw_phase)   = p_co3_cw02
        if(p_dust_cw02.ge.p1st) lptr_dust_aer(02,itype,cw_phase) = p_dust_cw02
        if(p_soa_cw02  .ge.p1st) lptr_soa_aer(02,itype,cw_phase)   = p_soa_cw02
        numptr_aer(02,itype,cw_phase)        = p_num_cw02
    end if

    if (nsize_aer(itype) .ge. 3) then
        if(p_oc_cw03  .ge.p1st) lptr_oc_aer(03,itype,cw_phase)   = p_oc_cw03
        if(p_bc_cw03  .ge.p1st) lptr_bc_aer(03,itype,cw_phase)   = p_bc_cw03
        if(p_so4_cw03  .ge.p1st) lptr_so4_aer(03,itype,cw_phase)   = p_so4_cw03
        if(p_nh4_cw03  .ge.p1st) lptr_nh4_aer(03,itype,cw_phase)   = p_nh4_cw03
        if(p_msa_cw03  .ge.p1st) lptr_msa_aer(03,itype,cw_phase)   = p_msa_cw03
        if(p_no3_cw03  .ge.p1st) lptr_no3_aer(03,itype,cw_phase)   = p_no3_cw03
        if(p_oin_cw03  .ge.p1st) lptr_oin_aer(03,itype,cw_phase)   = p_oin_cw03
        if(p_cl_cw03   .ge.p1st) lptr_cl_aer(03,itype,cw_phase)    = p_cl_cw03
        if(p_na_cw03   .ge.p1st) lptr_na_aer(03,itype,cw_phase)    = p_na_cw03
        if(p_ca_cw03   .ge.p1st) lptr_ca_aer(03,itype,cw_phase)    = p_ca_cw03
        if(p_co3_cw03  .ge.p1st) lptr_co3_aer(03,itype,cw_phase)   = p_co3_cw03
        if(p_dust_cw03.ge.p1st) lptr_dust_aer(03,itype,cw_phase) = p_dust_cw03
        if(p_soa_cw03  .ge.p1st) lptr_soa_aer(03,itype,cw_phase)   = p_soa_cw03
        numptr_aer(03,itype,cw_phase)        = p_num_cw03
    end if

    if (nsize_aer(itype) .ge. 4) then
        if(p_oc_cw04.ge.p1st)   lptr_oc_aer(04,itype,cw_phase)   = p_oc_cw04
        if(p_bc_cw04  .ge.p1st) lptr_bc_aer(04,itype,cw_phase)   = p_bc_cw04
        if(p_so4_cw04  .ge.p1st) lptr_so4_aer(04,itype,cw_phase)   = p_so4_cw04
        if(p_nh4_cw04  .ge.p1st) lptr_nh4_aer(04,itype,cw_phase)   = p_nh4_cw04
        if(p_msa_cw04  .ge.p1st) lptr_msa_aer(04,itype,cw_phase)   = p_msa_cw04
        if(p_no3_cw04  .ge.p1st) lptr_no3_aer(04,itype,cw_phase)   = p_no3_cw04
        if(p_oin_cw04  .ge.p1st) lptr_oin_aer(04,itype,cw_phase)   = p_oin_cw04
        if(p_cl_cw04   .ge.p1st) lptr_cl_aer(04,itype,cw_phase)    = p_cl_cw04
        if(p_na_cw04   .ge.p1st) lptr_na_aer(04,itype,cw_phase)    = p_na_cw04
        if(p_ca_cw04   .ge.p1st) lptr_ca_aer(04,itype,cw_phase)    = p_ca_cw04
        if(p_co3_cw04  .ge.p1st) lptr_co3_aer(04,itype,cw_phase)   = p_co3_cw04
        if(p_dust_cw04.ge.p1st) lptr_dust_aer(04,itype,cw_phase) = p_dust_cw04
        if(p_soa_cw04  .ge.p1st) lptr_soa_aer(04,itype,cw_phase)   = p_soa_cw04
        numptr_aer(04,itype,cw_phase)        = p_num_cw04
    end if

    if (nsize_aer(itype) .ge. 5) then
        if(p_cl_cw05.ge.p1st) lptr_cl_aer(05,itype,cw_phase) = p_cl_cw05
        if(p_na_cw05.ge.p1st) lptr_na_aer(05,itype,cw_phase) = p_na_cw05
        if(p_dust_cw05.ge.p1st) lptr_dust_aer(05,itype,cw_phase) = p_dust_cw05
        numptr_aer(05,itype,cw_phase)        = p_num_cw05
    end if

    if (nsize_aer(itype) .ge. 6) then
        if(p_cl_cw06.ge.p1st) lptr_cl_aer(06,itype,cw_phase) = p_cl_cw06
        if(p_na_cw06.ge.p1st) lptr_na_aer(06,itype,cw_phase) = p_na_cw06
        if(p_dust_cw06.ge.p1st) lptr_dust_aer(06,itype,cw_phase) = p_dust_cw06
        numptr_aer(06,itype,cw_phase)        = p_num_cw06
    end if

    if (nsize_aer(itype) .ge. 7) then
        if(p_cl_cw07.ge.p1st) lptr_cl_aer(07,itype,cw_phase) = p_cl_cw07
        if(p_na_cw07.ge.p1st) lptr_na_aer(07,itype,cw_phase) = p_na_cw07
        if(p_dust_cw07.ge.p1st) lptr_dust_aer(07,itype,cw_phase) = p_dust_cw07
        numptr_aer(07,itype,cw_phase)        = p_num_cw07
    end if

    if (nsize_aer(itype) .ge. 8) then
        if(p_cl_cw08.ge.p1st) lptr_cl_aer(08,itype,cw_phase) = p_cl_cw08
        if(p_na_cw08.ge.p1st) lptr_na_aer(08,itype,cw_phase) = p_na_cw08
        if(p_dust_cw08.ge.p1st) lptr_dust_aer(08,itype,cw_phase) = p_dust_cw08
        numptr_aer(08,itype,cw_phase)        = p_num_cw08
    end if

    if (nsize_aer(itype) .ge. 9) then
        if(p_cl_cw09.ge.p1st) lptr_cl_aer(09,itype,cw_phase) = p_cl_cw09
        if(p_na_cw09.ge.p1st) lptr_na_aer(09,itype,cw_phase) = p_na_cw09
        if(p_dust_cw09.ge.p1st) lptr_dust_aer(09,itype,cw_phase) = p_dust_cw09
        numptr_aer(09,itype,cw_phase)        = p_num_cw09
    end if

    if (nsize_aer(itype) .ge. 10) then
        if(p_cl_cw10.ge.p1st) lptr_cl_aer(10,itype,cw_phase) = p_cl_cw10
        if(p_na_cw10.ge.p1st) lptr_na_aer(10,itype,cw_phase) = p_na_cw10
        if(p_dust_cw10.ge.p1st) lptr_dust_aer(10,itype,cw_phase) = p_dust_cw10
        numptr_aer(10,itype,cw_phase)        = p_num_cw10
    end if

end if

!
!   define the massptr_aer and mastercompptr_aer pointers
!   and the name() species names
!

!   first initialize
    do l = 1, l2maxd
        write( name(l), '(a,i4.4,15x)' ) 'r', l
    end do
    massptr_aer(:,:,:,:) = -999888777
    mastercompptr_aer(:,:) = -999888777

do itype = 1, ntype_aer
    if (itype .eq. 1) then
        typetxt = ' '
        ntypetxt = 1
        if (ntype_aer .gt. 1) then
        typetxt = '_t1'
        ntypetxt = 3
        end if
    else if (itype .le. 9) then
        write(typetxt,'(a,i1)') '_t', itype
        ntypetxt = 3
    else if (itype .le. 99) then
        write(typetxt,'(a,i2)') '_t', itype
        ntypetxt = 4
    else
        typetxt = '_t?'
        ntypetxt = 3
    end if

    ncomp_dum(:,:) = 0
    ncomp_plustracer_dum(:,:) = 0

    do isize = 1, nsize_aer(itype)
        n =isize

        if (isize .le. 9) then
            write(sizetxt,'(i1)') isize
            nsizetxt = 1
        else if (isize .le. 99) then
            write(sizetxt,'(i2)') isize
            nsizetxt = 2
        else if (isize .le. 999) then
            write(sizetxt,'(i3)') isize
            nsizetxt = 3
        else
            sizetxt = 's?'
            nsizetxt = 2
        end if


        do iphase = 1, nphase_aer

            if (iphase .eq. ai_phase) then
                phasetxt = 'a'
                nphasetxt = 1
            else if (iphase .eq. cw_phase) then
                phasetxt = 'cw'
                nphasetxt = 2
            else 
                phasetxt = 'p?'
                nphasetxt = 2
            end if

            do l_mastercomp = -2, ntot_mastercomp_aer

                iaddto_ncomp = 1
                iaddto_ncomp_plustracer = 1

                if (l_mastercomp .eq. -2) then
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = numptr_aer(n,itype,iphase)
                    mcindx_dum = -2
                    spectxt = 'num_'
                    nspectxt = 4
                
                else if (l_mastercomp .eq. -1) then
                    if (iphase .ne. ai_phase) cycle
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = waterptr_aer(n,itype)
                    mcindx_dum = -1
                    spectxt = 'water_'
                    nspectxt = 6
                else if (l_mastercomp .eq. 0) then
                    if (iphase .ne. ai_phase) cycle
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = hyswptr_aer(n,itype)
                    mcindx_dum = 0
                    spectxt = 'hysw_'
                    nspectxt = 5

                else if (l_mastercomp .eq. mastercompindx_oc_aer) then
                    lptr_dum = lptr_oc_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_oc_aer
                    spectxt = 'oc_'
                    nspectxt = 3

                else if (l_mastercomp .eq. mastercompindx_bc_aer) then
                    lptr_dum = lptr_bc_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_bc_aer
                    spectxt = 'bc_'
                    nspectxt = 3
                
                else if (l_mastercomp .eq. mastercompindx_so4_aer) then
                    lptr_dum = lptr_so4_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_so4_aer
                    spectxt = 'so4_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_nh4_aer) then
                    lptr_dum = lptr_nh4_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_nh4_aer
                    spectxt = 'nh4_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_msa_aer) then
                    lptr_dum = lptr_msa_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_msa_aer
                    spectxt = 'msa_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_no3_aer) then
                    lptr_dum = lptr_no3_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_no3_aer
                    spectxt = 'no3_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_oin_aer) then
                    lptr_dum = lptr_oin_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_oin_aer
                    spectxt = 'oin_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_na_aer) then
                    lptr_dum = lptr_na_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_na_aer
                    spectxt = 'na_'
                    nspectxt = 3

                else if (l_mastercomp .eq. mastercompindx_cl_aer) then
                    lptr_dum = lptr_cl_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_cl_aer
                    spectxt = 'cl_'
                    nspectxt = 3

                else if (l_mastercomp .eq. mastercompindx_ca_aer) then
                    lptr_dum = lptr_ca_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_ca_aer
                    spectxt = 'ca_'
                    nspectxt = 3

                else if (l_mastercomp .eq. mastercompindx_co3_aer) then
                    lptr_dum = lptr_co3_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_co3_aer
                    spectxt = 'co3_'
                    nspectxt = 4

                else if (l_mastercomp .eq. mastercompindx_dust_aer) then
                    lptr_dum = lptr_dust_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_dust_aer
                    spectxt = 'dust_'
                    nspectxt = 5
                else if (l_mastercomp .eq. mastercompindx_soa_aer) then
                    lptr_dum = lptr_soa_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_soa_aer
                    spectxt = 'soa_'
                    nspectxt = 4
                else
                    ! goto 2500
                    cycle
                end if


                if (lptr_dum .gt. lmaxd) then
                    ! rce 2005-mar-14 - added check for lptr_dum > lmaxd
                    write( msg, '(a,3(1x,i4))' ) 'itype, isize, iphase =',   &
                    itype, isize, iphase
                    call chemistry_message( msg )
                    write( msg, '(a,3(1x,i4))' ) 'l_mastercomp, lptr_dum, lmaxd =',   &
                    l_mastercomp, lptr_dum, lmaxd
                    call chemistry_message( msg )
                    msg = '*** subr init_data_binaer_ptr error - lptr_dum > lmaxd'
                    call chemistry_error_fatal( msg )
                
                else if (lptr_dum .ge. p1st) then
                
                    ncomp_dum(isize,iphase) = ncomp_dum(isize,iphase) + iaddto_ncomp
                    ncomp_plustracer_dum(isize,iphase) =   &
                    ncomp_plustracer_dum(isize,iphase) + iaddto_ncomp_plustracer
                
                    name(lptr_dum) =   &
                    spectxt(1:nspectxt) // phasetxt(1:nphasetxt) //   &
                    sizetxt(1:nsizetxt) //  typetxt(1:ntypetxt)
                
                    if (l_mastercomp .eq. -2) then
                        !(numptr_aer is already set)
                        mprognum_aer(n,itype,iphase) = 1
                    
                    else if (l_mastercomp .eq. -1) then
                        !(waterptr_aer is already set)
                        continue
                    
                    else if (l_mastercomp .eq. 0) then
                        !(hyswptr_aer is already set)
                        continue
                    
                    else if (l_mastercomp .gt. 0) then
                        ll = ncomp_plustracer_dum(isize,iphase)
                        ! Note massptr_aer holds the species type l1, size bin n and is used at many places in movesect and binaer drydep code
                        massptr_aer(ll,n,itype,iphase) = lptr_dum
                        mastercompptr_aer(ll,itype) = mcindx_dum
                    
                        name_aer(ll,itype) = name_mastercomp_aer(mcindx_dum)
                        dens_aer(ll,itype) = dens_mastercomp_aer(mcindx_dum)
                        mw_aer(ll,itype) = mw_mastercomp_aer(mcindx_dum)
                        hygro_aer(ll,itype) = hygro_mastercomp_aer(mcindx_dum)
                    
                    end if
                end if ! if for lptr_dum
            end do    ! l_mastercomp = -1, ntot_mastercomp_aer 2500
        end do    ! iphase = 1, nphase_aer 2600
    end do    ! isize = 1, nsize_aer(itype) 2700


!   now set ncomp_aer and ncomp_plustracer_aer, 
!   *** and check that the values computed for each size and phase all match
    ncomp_aer(itype) = ncomp_dum(1,ai_phase)
    ncomp_plustracer_aer(itype) = ncomp_plustracer_dum(1,ai_phase)
    do iphase = 1, nphase_aer
        do isize = 1, nsize_aer(itype)
            if (ncomp_aer(itype) .ne. ncomp_dum(isize,iphase)) then
                msg =  '*** subr init_data_binaer_ptr - ' //   &
                'ncomp_aer .ne. ncomp_dum'
            call chemistry_message( msg )
            write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
            call chemistry_message( msg )
            write(msg,9350) 'ncomp_aer, ncomp_dum =',   &
                ncomp_aer(itype), ncomp_dum(isize,iphase)
            call chemistry_error_fatal( msg )
            end if
            if (ncomp_plustracer_aer(itype) .ne.   &
                ncomp_plustracer_dum(isize,iphase)) then
                msg = '*** subr init_data_binaer_ptr - ' //   &
                'ncomp_plustracer_aer .ne. ncomp_plustracer_dum'
            call chemistry_message( msg )
            write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
            call chemistry_message( msg )
            write(msg,9350)   &
                'ncomp_plustracer_aer, ncomp_plustracer_dum =',   &
                ncomp_plustracer_aer(itype),   &
                ncomp_plustracer_dum(isize,iphase)
            call chemistry_error_fatal( msg )
            end if
        end do
    end do


end do    ! itype = 1, ntype_aer 2800


9320    format( a, i1, i1, a, 8x )

!
!   output wrfch pointers
!

!jwfeng: 03/24/2020 
!modify to a nicer output format
9350    format( a, 32(1x,i4) )
    ! msg = ' '
    ! call chemistry_message( msg )
    ! msg = 'output from subr init_data_binaer_ptr'
    ! call chemistry_message( msg )
    write(msg,9350) '    nphase_aer = ', nphase_aer
    call chemistry_message( msg )

    do iphase=1,nphase_aer
    write(msg,'(a,i2)') '    Information of phase', iphase
    call chemistry_message( msg )
    msg = '    -------------------------------------------------------------------'
    call chemistry_message( msg )
    write(msg,9350) '            ntype_aer = ', ntype_aer
    call chemistry_message( msg )
    write(msg,9350) '            ncomp_aer = ', ncomp_aer
    call chemistry_message( msg )
      
    do itype=1,ntype_aer

    write(msg,9350) '                itype = ', itype
    call chemistry_message( msg )

    write(msg,9350) '            nsize_aer = ', nsize_aer(itype)
    call chemistry_message( msg )

    write(msg,9350) '          lptr_oc_aer = ',   &
        (lptr_oc_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '          lptr_bc_aer = ',   &
        (lptr_bc_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_so4_aer = ',   &
        (lptr_so4_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_no3_aer = ',   &
        (lptr_no3_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_oin_aer = ',   &
        (lptr_oin_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_cl_aer = ',   &
        (lptr_cl_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_na_aer = ',   &
        (lptr_na_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_ca_aer = ',   &
        (lptr_ca_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         lptr_co3_aer = ',   &
        (lptr_co3_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '        lptr_dust_aer = ',   &
        (lptr_dust_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '        lptr_nh4_aer = ',   &
        (lptr_nh4_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '        lptr_msa_aer = ',   &
        (lptr_msa_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '        lptr_soa_aer = ',   &
        (lptr_soa_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '          hyswptr_aer = ',   &
        (hyswptr_aer(n,itype), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         waterptr_aer = ',   &
        (waterptr_aer(n,itype), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '           numptr_aer = ',   &
        (numptr_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

9352    format( 5a )

    do ll = 1, ncomp_plustracer_aer(itype)
    ! write(msg,'(a,i3,a)') 'for comp No.', ll,':'
    write(msg,9350) '          massptr_aer = ',   &
        (massptr_aer(ll,n,itype,iphase), n=1,nsize_aer(itype)),ll
        call chemistry_message( msg )

            ! write(0,9352) 'chem_dname_table =  ', &
        !   (chem_dname_table(id,(massptr_aer(ll,n,itype,iphase)))(1:14),n=1,4)


    end do
    end do ! type
    msg = '    -------------------------------------------------------------------'
    call chemistry_message( msg )
    end do ! phase

!
!   check aerosol species pointers for "validity"
!
    do iphase=1,nphase_aer
    do itype=1,ntype_aer

    y_oc = 0
    y_bc = 0
    y_so4 = 0
    y_no3 = 0
    y_oin = 0
    y_dust = 0
    y_nh4 = 0
    y_msa = 0
    y_soa = 0
    y_cl = 0
    y_na = 0
    y_ca = 0
    y_co3 = 0
    y_hysw = 0
    y_water = 0
    y_num = 0

    do n = 1, nsize_aer(itype)
        if (lptr_oc_aer(n,itype,iphase)   .ge. p1st) y_oc   = y_oc + 1
        if (lptr_bc_aer(n,itype,iphase)   .ge. p1st) y_bc   = y_bc + 1
        if (lptr_so4_aer(n,itype,iphase)  .ge. p1st) y_so4   = y_so4 + 1
        if (lptr_no3_aer(n,itype,iphase)  .ge. p1st) y_no3   = y_no3 + 1
        if (lptr_oin_aer(n,itype,iphase)  .ge. p1st) y_oin   = y_oin + 1
        if (lptr_dust_aer(n,itype,iphase) .ge. p1st) y_dust = y_dust + 1
        if (lptr_nh4_aer(n,itype,iphase)  .ge. p1st) y_nh4 = y_nh4 + 1
        if (lptr_msa_aer(n,itype,iphase)  .ge. p1st) y_msa = y_msa + 1
        if (lptr_soa_aer(n,itype,iphase)  .ge. p1st) y_soa = y_soa + 1
        if (lptr_cl_aer(n,itype,iphase)   .ge. p1st) y_cl  = y_cl + 1
        if (lptr_na_aer(n,itype,iphase)   .ge. p1st) y_na  = y_na + 1
        if (lptr_ca_aer(n,itype,iphase)   .ge. p1st) y_ca  = y_ca + 1
        if (lptr_co3_aer(n,itype,iphase)  .ge. p1st) y_co3 = y_co3 + 1
        if (hyswptr_aer(n,itype)          .ge. p1st) y_hysw = y_hysw + 1
        if (waterptr_aer(n,itype)         .ge. p1st) y_water = y_water + 1
        if (numptr_aer(n,itype,iphase)    .ge. p1st) y_num = y_num + 1
    end do

!   these must be defined for all aerosol bins
!    if (y_so4 .ne. nsize_aer(itype)) then
    ! if (y_dust .ne. nsize_aer(itype)) then
    !     msg = '*** subr init_data_binaer_ptr - ptr error for dust'
    !     call chemistry_message( msg )
    !     write(msg,9350) 'phase, type=', iphase,itype
    !     call chemistry_error_fatal( msg )
    ! else if (y_seas .ne. nsize_aer(itype)) then
    !     msg = '*** subr init_data_binaer_ptr - ptr error for seas'
    !     call chemistry_message( msg )
    !     write(msg,9350) 'phase, type=', iphase,itype
    !     call chemistry_error_fatal( msg )
    ! ! else if (y_oc .ne. nsize_aer(itype)) then
    ! !     msg = '*** subr init_data_binaer_ptr - ptr error for oc'
    ! !     call chemistry_message( msg )
    ! !     write(msg,9350) 'phase, type=', iphase,itype
    !     call chemistry_error_fatal( msg )
    ! else if (y_hysw .ne. nsize_aer(itype)) then
    !     msg = '*** subr init_data_binaer_ptr - ptr error for hysw'
    !     call chemistry_message( msg )
    !     write(msg,9350) 'phase, type=', iphase,itype
    !     call chemistry_error_fatal( msg )
    ! else if (y_water .ne. nsize_aer(itype)) then
    !     msg = '*** subr init_data_binaer_ptr - ptr error for water'
    !     call chemistry_message( msg )
    !     write(msg,9350) 'phase, type=', iphase,itype
    !     call chemistry_error_fatal( msg )
    if (y_num .ne. nsize_aer(itype)) then
        msg = '*** subr init_data_binaer_ptr - ptr error for num'
        call chemistry_message( msg )
        write(msg,9350) 'phase, type=', iphase,itype
        call chemistry_error_fatal( msg )
    end if


!   these must be defined for all aerosol bins
!       or else undefined for all aerosol bins
!    if ((y_no3 .ne. 0) .and. (y_no3 .ne. nsize_aer(itype))) then
!        msg = '*** subr init_data_binaer_ptr - ptr error for no3'
!        call chemistry_message( msg )
!        write(msg,9350) 'phase, type=', iphase,itype
!        call chemistry_error_fatal( msg )
!   else if ((y_dust .ne. 0) .and. (y_dust .ne. nsize_aer(itype))) then
!   if ((y_dust .ne. 0) .and. (y_dust .ne. nsize_aer(itype))) then
!       msg = '*** subr init_data_binaer_ptr - ptr error for dust'
!       call chemistry_message( msg )
!       write(msg,9350) 'phase, type=', iphase,itype
!       call chemistry_error_fatal( msg )
!   end if

    enddo ! type
    enddo ! phase

!   set pointers for gases
!  zhxia :When we add new gas species, we need to consider whether to update this. 10/2024

    if (p_h2so4 .ge. p1st) then
        kh2so4 = p_h2so4
    endif

    if (p_hno3 .ge. p1st) then
         khno3 = p_hno3
    endif

    if (p_hcl .ge. p1st) then
          khcl = p_hcl
    endif

    if (p_nh3 .ge. p1st) then
          knh3 = p_nh3
    endif

    if (p_msa .ge. p1st) then
          kmsa = p_msa
    endif

    if (p_n2o5 .ge. p1st) then
         kn2o5 = p_n2o5
    endif

    if (p_o3 .ge. p1st) then
           ko3 = p_o3
    endif

!   rce 2005-apr-12 - added following species for cldchem, here and below:
!   kso2, kh2o2, khcho, khcooh, koh, kho2,
!   kno3, kno, kno2, khono, kpan, kch3o2, kch3oh, kch3ooh
    if (p_so2    .ge. p1st) kso2    = p_so2
    if (p_h2o2   .ge. p1st) kh2o2   = p_h2o2
    if (p_hcho   .ge. p1st) khcho   = p_hcho
    if (p_ho     .ge. p1st) koh     = p_ho
    if (p_ho2    .ge. p1st) kho2    = p_ho2
    if (p_no3    .ge. p1st) kno3    = p_no3
    if (p_no     .ge. p1st) kno     = p_no
    if (p_no2    .ge. p1st) kno2    = p_no2
    if (p_hono   .ge. p1st) khono   = p_hono
    if (p_pan    .ge. p1st) kpan    = p_pan
!
!   calc ltot, ltot2, kh2o, ktemp
!
    is_aerosol(:) = .false.
    ltot = 0
    ltot = max( ltot, kh2so4 )
    ltot = max( ltot, khno3 )
    ltot = max( ltot, khcl )
    ltot = max( ltot, knh3 )
    ltot = max( ltot, kmsa )
    ltot = max( ltot, kn2o5 )
    ltot = max( ltot, ko3 )
    ltot = max( ltot, kso2    )
    ltot = max( ltot, kh2o2   )
    ltot = max( ltot, khcho   )
    ltot = max( ltot, koh     )
    ltot = max( ltot, kho2    )
    ltot = max( ltot, kno3    )
    ltot = max( ltot, kno     )
    ltot = max( ltot, kno2    )
    ltot = max( ltot, khono   )
    ltot = max( ltot, kpan    )

    do iphase=1,nphase_aer
        do itype=1,ntype_aer
        do n = 1, nsize_aer(itype)
            do ll = 1, ncomp_plustracer_aer(itype)
               ltot = max( ltot, massptr_aer(ll,n,itype,iphase) )
               is_aerosol(massptr_aer(ll,n,itype,iphase))=.true.
            end do
            ltot = max( ltot, hyswptr_aer(n,itype) )
            ltot = max( ltot, waterptr_aer(n,itype) )
            ltot = max( ltot, numptr_aer(n,itype,iphase) )
            l = hyswptr_aer(n,itype)
            if (l .ge. p1st) is_aerosol(l)=.true.
            l = waterptr_aer(n,itype)
            if (l .ge. p1st) is_aerosol(l)=.true.
            l = numptr_aer(n,itype,iphase)
            if (l .ge. p1st) is_aerosol(l)=.true.
        end do
        end do
    end do


    kh2o = ltot + 1
    ktemp = ltot + 2
    ltot2 = ktemp

    write( msg, '(a,4(1x,i4))' ) 'ltot, ltot2, lmaxd, l2maxd =',   &
        ltot, ltot2, lmaxd, l2maxd
    call chemistry_message( msg )
    if ((ltot .gt. lmaxd) .or. (ltot2 .gt. l2maxd)) then
        msg = '*** subr init_data_binaer_ptr - ltot/ltot2 too big'
        call chemistry_error_fatal( msg )
    end if

    if (p_h2so4   .ge. p1st) name(kh2so4 ) = 'h2so4'
    if (p_hno3    .ge. p1st) name(khno3  ) = 'hno3'
    if (p_hcl    .ge. p1st) name(khcl   ) = 'hcl'
    if (p_nh3    .ge. p1st) name(knh3   ) = 'nh3'
    if (p_msa    .ge. p1st) name(kmsa   ) = 'msa'
    if (p_n2o5   .ge. p1st) name(kn2o5  ) = 'n2o5'
    if (p_o3     .ge. p1st) name(ko3    ) = 'o3'
    if (p_so2    .ge. p1st) name(kso2   ) = 'so2'
    if (p_h2o2   .ge. p1st) name(kh2o2  ) = 'h2o2'
    if (p_hcho   .ge. p1st) name(khcho  ) = 'hcho'
    if (p_ho     .ge. p1st) name(koh    ) = 'oh'
    if (p_ho2    .ge. p1st) name(kho2   ) = 'ho2'
    if (p_no3    .ge. p1st) name(kno3   ) = 'no3'
    if (p_no     .ge. p1st) name(kno    ) = 'no'
    if (p_no2    .ge. p1st) name(kno2   ) = 'no2'
    if (p_hono   .ge. p1st) name(khono  ) = 'hono'
    if (p_pan    .ge. p1st) name(kpan   ) = 'pan'
    name(ktemp)  = 'temp'
    name(kh2o)   = 'h2o'

    call initwet(   &
    ntype_aer, nsize_aer, ncomp_aer,   &
    massptr_aer, dens_aer, numptr_aer,           &
    maxd_acomp, maxd_asize,maxd_atype, maxd_aphase,   &
    dcen_sect, sigmag_aer, &
    waterptr_aer, dens_water_aer, &
    scavimptblvol, scavimptblnum, nimptblgrow_mind,   &
    nimptblgrow_maxd, dlndg_nimptblgrow)

    return
end subroutine init_data_binaer_ptr

end module iamas_atmchem_binaer_init
