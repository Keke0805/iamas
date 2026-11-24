!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to set gases concentrations for chem(:,:,:,:).
! * Author: Zihan Xia from LACAR USTC, xiazh@mail.ustc.edu.cn
! Last update: zhxia 05/2025
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! References:
! module_input_chem_data.F in wrfchem.
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

MODULE iamas_atmchem_gas_init
   USE iamas_atmchem_vars
   USE iamas_atmphys_vars
   USE iamas_kind_types
   USE iamas_pool_routines
   USE iamas_atmchem_utilities
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial atmospheric chemistry profile data
    INTEGER :: k_loop          ! Used for loop index
    INTEGER :: lo              ! number of chemicals in inital profile
    INTEGER :: logg            ! number of final chemical species (nch-1)
    INTEGER :: kx              ! number of vertical levels in temp profile
    INTEGER :: kxm1

    PARAMETER( kx=16, kxm1=kx-1, logg=350, lo=34) ! DL (6/2/2013) Changed value of logg from 200 to 350 for additional gas species
    INTEGER, DIMENSION(logg)                     :: iref
    INTEGER, DIMENSION(logg)                     :: gas_index 
    REAL, DIMENSION(logg)                        :: fracref
    REAL, DIMENSION(kx)                          :: dens
    REAL, DIMENSION(kx+1)                        :: zfa
    REAL, DIMENSION(kx+1)                        :: zfa_bdy
    REAL, DIMENSION(lo  ,kx)                     :: xl
    REAL                                         :: so4vaptoaer
    DATA so4vaptoaer/.999/
    CHARACTER (LEN=20), DIMENSION(logg) :: ggnam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The idealized profile is converted from the NALROM chemistry model.
!
!  The variable "iref" is the reference index
!    "fracref" is the reference fraction correpsonding to iref
!
!  For example: WRF-Chem species number 1 for RADM is SO2.  iref(1) is 12, and XL(12,K) is the profile for SO2.
!
!  Note: NALROM calculates lumped OX (XL(1) =O3+NO2+HNO3+...) and a lumped NOX (XL(2)=NO+NO2+NO3+2N2O5+HO2NO2+HONO).  But XL(33) is strictly O3, and XL(34)=NOx=(NO+NO2 only).
!
!  Short-lived species are initialized to steady-state equilibrium since they are short-lived.  The short-lived species within a lumped category (Ox , NOx, or NO3+N2O5 in our case) would be renormalized to the lumped class after the steady-state equilibrium concentrations are determined.
!
!  The following is a list of long-lived species
!
!      NAMEL( 1)='OX        '
!      NAMEL( 2)='NOX       '
!      NAMEL( 3)='HNO3      '
!      NAMEL( 4)='H2O2      '
!      NAMEL( 5)='CH3OOH    '
!      NAMEL( 6)='CO        '
!      NAMEL( 7)='ISOPRENE  '
!      NAMEL( 8)='CH2O      '
!      NAMEL( 9)='CH3CHO    '
!      NAMEL(10)='PAN       '
!      NAMEL(11)='OTHER ALKA'
!      NAMEL(12)='SO2       '
!      NAMEL(13)='BUTANE    '
!      NAMEL(14)='ETHENE    '
!      NAMEL(15)='PROPENE+   '
!      NAMEL(16)='PPN       '
!      NAMEL(17)='MEK       '
!      NAMEL(18)='RCHO      '
!      NAMEL(19)='SO4=      '
!      NAMEL(20)='MVK       '
!      NAMEL(21)='MACR      '
!      NAMEL(22)='HAC       '
!      NAMEL(23)='MGLY      '
!      NAMEL(24)='HPAN      '
!      NAMEL(25)='MPAN      '
!      NAMEL(26)='PROPANE   '
!      NAMEL(27)='ACETYLENE '
!      NAMEL(28)='OH        '
!      NAMEL(29)='HO2       '
!      NAMEL(30)='NO3+N2O5  '
!      NAMEL(31)='HO2NO2    '
!      NAMEL(32)='SUM RO2   '
!      NAMEL(33)='OZONE     '
!      NAMEL(34)='NOX       '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DATA dens/ 2.738E+18, 5.220E+18, 7.427E+18, 9.202E+18, &
             1.109E+19, 1.313E+19, 1.525E+19, 1.736E+19, &
             1.926E+19, 2.074E+19, 2.188E+19, 2.279E+19, &
             2.342E+19, 2.384E+19, 2.414E+19, 2.434E+19  /

  DATA ZFA_BDY/    0.,   85.,  212.,  385.,  603.,  960., 1430., 2010., &
            2850., 4010., 5340., 6900., 8510.,10200.,12100.,16000., &
           21000./

! Profile heights in meters
  DATA ZFA/    0.,   85.,  212.,  385.,  603.,  960., 1430., 2010., &
            2850., 4010., 5340., 6900., 8510.,10200.,12100.,16000., &
           21000./

!wig: To match the xl profile to the correct species, match WRF's p_<species> flag with iref(p_<species>-1) to get the value of the first index in xl, e.g. p_o3=6, iref(6-1)=1, so xl(1,:) is the ozone profile. See gasprofile_init_pnnl for an explination of what height each index represents.

!O3-czhao
      DATA (xl(1,k_loop),k_loop=1,kx) &
      / 1.68E-07, 1.68E-07, 5.79E-08, 5.24E-08, 5.26E-08, &
       5.16E-08, 4.83E-08, 4.50E-08, 4.16E-08, 3.80E-08, 3.56E-08, &
       3.35E-08, 3.15E-08, 3.08E-08, 3.06E-08, 3.00E-08/
!NO2 & NO -czhao
      DATA (xl(2,k_loop),k_loop=1,kx) &
      / 4.06E-10, 4.06E-10, 2.16E-10, 1.37E-10, 9.47E-11, &
       6.95E-11, 5.31E-11, 4.19E-11, 3.46E-11, 3.01E-11, 2.71E-11, &
       2.50E-11, 2.35E-11, 2.26E-11, 2.20E-11, 2.16E-11/
!HNO3-czhao
      DATA (xl(3,k_loop),k_loop=1,kx) &
      / 9.84E-10, 9.84E-10, 5.66E-10, 4.24E-10, 3.26E-10, &
       2.06E-10, 1.12E-10, 7.33E-11, 7.03E-11, 7.52E-11, 7.96E-11, &
       7.56E-11, 7.27E-11, 7.07E-11, 7.00E-11, 7.00E-11/
!H2O2-czhao
      DATA (xl(4,k_loop),k_loop=1,kx) &
      / 8.15E-10, 8.15E-10, 8.15E-10, 8.15E-10, 8.15E-10, &
       8.65E-10, 1.07E-09, 1.35E-09, 1.47E-09, 1.47E-09, 1.47E-09, &
       1.47E-09, 1.45E-09, 1.43E-09, 1.40E-09, 1.38E-09/
!op1,op2-czhao
      DATA (xl(5,k_loop),k_loop=1,kx) &
      / 4.16E-10, 4.16E-10, 4.16E-10, 4.16E-10, 4.16E-10, &
       4.46E-10, 5.57E-10, 1.11E-09, 1.63E-09, 1.63E-09, 1.63E-09, &
       1.63E-09, 1.61E-09, 1.59E-09, 1.57E-09, 1.54E-09/

!  CO is 70 ppbv at top, 80 throughout troposphere
!ORA1,ORA2,NH3,ol2-czhao
      DATA (xl(6,k_loop),k_loop=1,kx)  / 7.00E-08, kxm1*8.00E-08/
!iso-czhao, isoprene
      DATA (xl(7,k_loop),k_loop=1,kx) &
      / 8.33E-29, 8.33E-29, 8.33E-29, 8.33E-29, 8.33E-29, &
       1.33E-28, 3.54E-28, 1.85E-28, 1.29E-29, 1.03E-30, 1.72E-31, &
       7.56E-32, 1.22E-31, 2.14E-31, 2.76E-31, 2.88E-31/
!hcho-czhao
      DATA (xl(8,k_loop),k_loop=1,kx) &
      / 9.17E-11, 9.17E-11, 9.17E-11, 9.17E-11, 9.17E-11, &
       1.03E-10, 1.55E-10, 2.68E-10, 4.47E-10, 4.59E-10, 4.72E-10, &
       4.91E-10, 5.05E-10, 5.13E-10, 5.14E-10, 5.11E-10/
!ALD-czhao
      DATA (xl(9,k_loop),k_loop=1,kx) &
      / 7.10E-12, 7.10E-12, 7.10E-12, 7.10E-12, 7.10E-12, &
       7.36E-12, 1.02E-11, 2.03E-11, 2.98E-11, 3.01E-11, 3.05E-11, &
       3.08E-11, 3.08E-11, 3.06E-11, 3.03E-11, 2.99E-11/
!hc3-czhao
      DATA (xl(10,k_loop),k_loop=1,kx) &
      / 4.00E-11, 4.00E-11, 4.00E-11, 3.27E-11, 2.51E-11, &
       2.61E-11, 2.20E-11, 1.69E-11, 1.60E-11, 1.47E-11, 1.37E-11, &
       1.30E-11, 1.24E-11, 1.20E-11, 1.18E-11, 1.17E-11/
!eth-czhao
      DATA (xl(11,k_loop),k_loop=1,kx) &
      / 1.15E-16, 1.15E-16, 2.46E-15, 2.30E-14, 1.38E-13, &
       6.25E-13, 2.31E-12, 7.32E-12, 1.87E-11, 3.68E-11, 6.10E-11, &
       9.05E-11, 1.22E-10, 1.50E-10, 1.70E-10, 1.85E-10/
!SO2-czhao
      DATA (xl(12,k_loop),k_loop=1,kx) &
      / 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, &
       1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, &
       1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10, 1.00E-10/
!hc8-czhao
      DATA (xl(13,k_loop),k_loop=1,kx) &
      / 1.26E-11, 1.26E-11, 2.02E-11, 2.50E-11, 3.02E-11, &
       4.28E-11, 6.62E-11, 1.08E-10, 1.54E-10, 2.15E-10, 2.67E-10, &
       3.24E-10, 3.67E-10, 3.97E-10, 4.16E-10, 4.31E-10/
!olt-czhao
      DATA (xl(14,k_loop),k_loop=1,kx) &
      / 1.15E-16, 1.15E-16, 2.46E-15, 2.30E-14, 1.38E-13, &
       6.25E-13, 2.31E-12, 7.32E-12, 1.87E-11, 3.68E-11, 6.10E-11, &
       9.05E-11, 1.22E-10, 1.50E-10, 1.70E-10, 1.85E-10/
!oli,tol-czhao
      DATA (xl(15,k_loop),k_loop=1,kx) &
      / 1.00E-20, 1.00E-20, 6.18E-20, 4.18E-18, 1.23E-16, &
       2.13E-15, 2.50E-14, 2.21E-13, 1.30E-12, 4.66E-12, 1.21E-11, &
       2.54E-11, 4.47E-11, 6.63E-11, 8.37E-11, 9.76E-11/
      DATA (xl(16,k_loop),k_loop=1,kx) &
      / 1.23E-11, 1.23E-11, 1.23E-11, 1.23E-11, 1.23E-11, &
       1.20E-11, 9.43E-12, 3.97E-12, 1.19E-12, 1.11E-12, 9.93E-13, &
       8.66E-13, 7.78E-13, 7.26E-13, 7.04E-13, 6.88E-13/
!ket-czhao
      DATA (xl(17,k_loop),k_loop=1,kx) &
      / 1.43E-12, 1.43E-12, 1.43E-12, 1.43E-12, 1.43E-12, &
       1.50E-12, 2.64E-12, 8.90E-12, 1.29E-11, 1.30E-11, 1.32E-11, &
       1.32E-11, 1.31E-11, 1.30E-11, 1.29E-11, 1.27E-11/
!RCHO-czhao
      DATA (xl(18,k_loop),k_loop=1,kx) &
       / 3.61E-13, 3.61E-13, 3.61E-13, 3.61E-13, 3.61E-13, &
       3.58E-13, 5.22E-13, 1.75E-12, 2.59E-12, 2.62E-12, 2.64E-12, &
       2.66E-12, 2.65E-12, 2.62E-12, 2.60E-12, 2.57E-12/
!SULF-czhao
      DATA (xl(19,k_loop),k_loop=1,kx) &
       / 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, &
       5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, &
       5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11, 5.00E-11/

!xyl,aco3,hono-czhao
      DATA (xl(20,k_loop),k_loop=1,kx)/kx*1.E-20/
      DATA (xl(21,k_loop),k_loop=1,kx)/kx*1.E-20/
      DATA (xl(22,k_loop),k_loop=1,kx)/kx*1.E-20/
!gly,mgly,dcb,
      DATA (xl(23,k_loop),k_loop=1,kx)/kx*1.E-20/
      DATA (xl(24,k_loop),k_loop=1,kx)/kx*1.E-20/
      DATA (xl(25,k_loop),k_loop=1,kx)/kx*1.E-20/

! Propane - Gregory PEM-West A 25 ppt median marine boundary layer
!hc5-czhao
      DATA (xl(26,k_loop),k_loop=1,kx) &
      /5.00E-13, 1.24E-12, 2.21E-12, 3.27E-12, 4.71E-12, &
       6.64E-12, 9.06E-12, 1.19E-11, 1.47E-11, 1.72E-11, &
       1.93E-11, 2.11E-11, 2.24E-11, 2.34E-11, 2.42E-11, 2.48E-11/
! Acetylene - Gregory PEM-West A 53 ppt median marine boundary layer
      DATA (xl(27,k_loop),k_loop=1,kx) &
      /1.00E-12, 2.48E-12, 4.42E-12, 6.53E-12, 9.42E-12, &
       1.33E-11, 1.81E-11, 2.37E-11, 2.95E-11, 3.44E-11, &
       3.85E-11, 4.22E-11, 4.49E-11, 4.69E-11, 4.84E-11, 4.95E-11/
! OH
      DATA (xl(28,k_loop),k_loop=1,kx) &
       / 9.80E+06, 9.80E+06, 4.89E+06, 2.42E+06, 1.37E+06, &
       9.18E+05, 7.29E+05, 6.26E+05, 5.01E+05, 4.33E+05, 4.05E+05, &
       3.27E+05, 2.54E+05, 2.03E+05, 1.74E+05, 1.52E+05/
! HO2
      DATA (xl(29,k_loop),k_loop=1,kx) &
       / 5.74E+07, 5.74E+07, 7.42E+07, 8.38E+07, 8.87E+07, &
       9.76E+07, 1.15E+08, 1.34E+08, 1.46E+08, 1.44E+08, 1.40E+08, &
       1.36E+08, 1.31E+08, 1.28E+08, 1.26E+08, 1.26E+08/
! NO3+N2O5
      DATA (xl(30,k_loop),k_loop=1,kx) &
       / 5.52E+05, 5.52E+05, 3.04E+05, 2.68E+05, 2.32E+05, &
       1.66E+05, 1.57E+05, 1.72E+05, 1.98E+05, 2.22E+05, 2.43E+05, &
       2.75E+05, 3.00E+05, 3.18E+05, 3.32E+05, 3.39E+05/
! HO2NO2-czhao
      DATA (xl(31,k_loop),k_loop=1,kx) &
       / 7.25E+07, 7.25E+07, 6.36E+07, 5.55E+07, 4.94E+07, &
       3.66E+07, 2.01E+07, 9.57E+06, 4.75E+06, 2.37E+06, 1.62E+06, &
       9.86E+05, 7.05E+05, 5.63E+05, 4.86E+05, 4.41E+05/
! Sum of RO2 &
      DATA (xl(32,k_loop),k_loop=1,kx) &
       / 9.14E+06, 9.14E+06, 1.46E+07, 2.14E+07, 2.76E+07, &
       3.62E+07, 5.47E+07, 1.19E+08, 2.05E+08, 2.25E+08, 2.39E+08, &
       2.58E+08, 2.82E+08, 2.99E+08, 3.08E+08, 3.15E+08/
! O3 <--This is not the O3 used for RADM2 or CBMZ (wig)
!in molecule number concentration
      DATA (xl(33,k_loop),k_loop=1,kx) &
       / 8.36E+11, 8.36E+11, 4.26E+11, 4.96E+11, 6.05E+11, &
       6.93E+11, 7.40E+11, 7.74E+11, 7.82E+11, 7.75E+11, 7.69E+11, &
       7.59E+11, 7.54E+11, 7.50E+11, 7.47E+11, 7.45E+11/
! NOx (NO+NO2) in molecule number concentration
      DATA (xl(34,k_loop),k_loop=1,kx) &
       / 1.94E+09, 1.94E+09, 1.53E+09, 1.24E+09, 1.04E+09, &
       8.96E+08, 7.94E+08, 7.11E+08, 6.44E+08, 6.00E+08, 5.70E+08, &
       5.49E+08, 5.35E+08, 5.28E+08, 5.24E+08, 5.23E+08/

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gas_init(configs, mesh)
 type(iamas_pool_type),intent(in) :: configs
 type(iamas_pool_type),intent(in):: mesh 

 ! local variables
 character(len=StrKIND),pointer :: chem_opt,gaschem_opt
 real(kind=RKIND), dimension(:,:),allocatable :: zgrid_t 
 Integer :: k,i

 call iamas_pool_get_config(configs,'config_chem_opt' ,chem_opt) 
 call iamas_pool_get_config(configs,'config_gaschem_opt' ,gaschem_opt)
 call iamas_pool_get_array(mesh,'zgrid'   ,zgrid   )
 
 if (trim(chem_opt) == 'bin04aer_gas_aq' .OR. &
     trim(chem_opt) == 'bin04aer_gas') then

     allocate(zgrid_t(kts:kte, its:ite))

     do k = kts, kte
        do i = its, ite
           zgrid_t(k,i) = 0.5*(zgrid(k,i)+zgrid(k+1,i))
        enddo
     enddo

     call setup_gasprofile_maps(gaschem_opt,num_gas)
     call make_chem_profile (nx1 = its, nx2 = ite, nz1 = kts, nz2 = kte, &
                             nch = num_chem, numgas = num_gas,           &
                             zgrid=zgrid_t, chem = chem(1:num_chem,kts:kte,its:ite)   )

     do i = its, ite
        do k = kts, kte
           chem(p_co2,k,i) = 380
           chem(p_ch4,k,i) = 1.7
        enddo
     enddo

     deallocate(zgrid_t)
 endif

END SUBROUTINE gas_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE setup_gasprofile_maps(gaschem_opt, numgas)
  character(len=StrKIND), intent(in) :: gaschem_opt 
  integer, intent(in) :: numgas

   !gaschem_opt_select: SELECT CASE(trim(gaschem_opt))
   !  CASE ("AIPC_saprc99")
          call setup_gasprofile_map_saprcnov
   !  END SELECT gaschem_opt_select

END SUBROUTINE setup_gasprofile_maps 

SUBROUTINE setup_gasprofile_map_saprcnov
  iref(:)    = 7 !default value
  iref(1:49) = (/1,4,2,2,30,30,23,3,31,12, &
                  19,6,8,9,18,17,6,6,6,23, &
                  23,23,21,20,14,7,26,27,15,13, &
                  11,11,23,23,6,6,23,10,16,25, &
                  32,32,32,28,29,5,5,32,32/)

  fracref(:)    = 1. !default value
  fracref(1:49) = (/1.00E+00,1.00E+00,2.50E-01,7.50E-01,9.00E-01,1.00E-01, &
                  1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00, &
                  1.00E+00,1.00E+00,1.00E+00,1.00E+00,6.25E-04,5.00E-04, &
                  2.50E-04,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00, &
                  1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00, &
                  5.00E-01,5.00E-01,1.00E+00,1.00E+00,8.00E-03,8.00E-03, &
                  1.00E+00,1.00E+00,1.00E+00,1.00E+00,3.00E-01,2.00E-01, &
                  2.00E-01,1.00E+00,1.00E+00,1.00E+00,5.00E-01,1.00E+00, &
                  3.00E-01/)
  ggnam(:) = 'None' !default value
  ggnam(1:49) = (/'O3  ','H2O2','NO  ','NO2 ','NO3 ','N2O5', &
                  'HONO','HNO3','HNO4','SO2 ','H2SO','CO  ', & !H2SO really is H2SO4
                  'HCHO','CCHO','RCHO','MEK ','HCOH','CCOH', & !HCOH=HCOOH, CCOH=CCO_OH
                  'RCOH','GLY ','MGLY','CRES','MEAC','MVK ', & !RCOH=RCO_OH, MEAC=METHACRO
                  'ETHE','ISPR','C3H8','C2H2','C3H6','ALK3', & !ETHE=ETHENE, ISPR=ISOPRENE
                  'ALK4','ALK5','ARO1','ARO2','OLE1','OLE2', &
                  'RNO3','PAN ','PAN2','MAPN','COOH','RCO2', & !MAPN=MA_PAN, COOH=CCO_OOH, RCO2=RCO_O2
                  'ROOH','OH  ','HO2 ','COOH','ROOH','RO2R', & !ROOH=RCO_OOH, RO2R=RO2_R
                  'RO2N'/)                                     !RO2R=RO2_R,RO2N=RO2_N

  ! There is a bug. The most of vocs species is mismatch. @zhxia
  gas_index(1:49) = (/p_o3,p_h2o2,p_no,p_no2,p_no3,p_n2o5, &
                  p_hono,p_hno3,p_hno4,p_so2,p_h2so4,p_co, &
                  p_hcho,p_ccho,p_rcho,p_mek,p_hcooh,p_cco_oh, & 
                  p_rco_oh,p_gly,p_mgly,p_cres,p_methacro,p_mvk, & 
                  p_ethene,p_isoprene,p_c3h8,p_c2h2,p_c3h6,p_alk3, &
                  p_alk4,p_alk5,p_aro1,p_aro2,p_ole1,p_ole2, &
                  p_rno3,p_pan,p_pan2,p_ma_pan,p_cco_ooh,p_rco_o2, &
                  p_rco_ooh,p_ho,p_ho2,p_cooh,p_rooh,p_ro2_r, &
                  p_ro2_n/)  
END SUBROUTINE setup_gasprofile_map_saprcnov


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE make_chem_profile ( nx1, nx2, nz1, nz2, nch, numgas, &
                               zgrid, chem )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx1, nz1
  INTEGER, INTENT(IN) :: nx2, nz2
  INTEGER, INTENT(IN) :: nch, numgas
  REAL, DIMENSION(nz1:nz2,nx1:nx2), INTENT(IN) :: zgrid
  REAL, DIMENSION(nch,nz1:nz2,nx1:nx2),INTENT(INOUT) :: chem

  CHARACTER (LEN=80) :: message
  INTEGER :: i, k, l, is, ix

  REAL, DIMENSION(lo,nz1:kx,nx1:nx2):: chprof
  REAL, DIMENSION(nz1:kx,nx1:nx2)     :: zprof
  REAL, DIMENSION(lo,nz1:nz2,nx1:nx2 ) :: stor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!   do k=1,kde-1
!     si_zsig(:,k,:) = 0.5 * ( si_zsigf(:,k,:) + si_zsigf(:,k+1,:) )
!   enddo
!   si_zsig(:,kde,:) = 0.5 * ( 3. * si_zsigf(:,kde,:) - si_zsigf(:,kde-1,:) )

      ! Vertically flip the chemistry data as it is given top down and
      ! heights are bottom up. Fill temp 3D chemical and profile array,
      ! keep chem slot 1 open as vinterp_chem assumes there is no data.
      DO k=  1,kx
      DO i=nx1,nx2
         chprof(1:lo,k,i) = xl(1:lo,kx-k+1)
         zprof(k,i) = 0.5*(zfa(k)+zfa(k+1))
      ENDDO
      ENDDO
!
! return xl to previous value for next time...
!   34 chemicals (lo), 16 vertical levels (kx)
!     DO i=lo-6,lo
!       xl(i,1:kx)=xl(i,1:kx)*dens(1:kx)
!     ENDDO

! Change number concentrations to mixing ratios for short-lived NALROM species
     do k=1,kx
        chprof(lo-6:lo,k,:) = chprof(lo-6:lo,k,:)/dens(k)
     end do

     ! Interpolate temp 3D chemical and profile array to iAMAS grid
     call vinterp_chem(nx1, nx2, nz1, kx, nz2, nch, zprof, zgrid, &
                         chprof, stor,.false.)

     ! Here is where the chemistry profile is constructed
     chem(1:numgas,nz1:nz2,nx1:nx2)  = 1. * SPREAD(stor(7,nz1:nz2,nx1:nx2), dim=1, ncopies=numgas) *1.E6 !default values
     DO l=1,49 !saprcnov=49
        is=iref(l)
        ix=gas_index(l) 
        chem(ix,nz1:nz2,nx1:nx2)=fracref(l)*stor(is,nz1:nz2,nx1:nx2)*1.E6
     ENDDO
     RETURN
  END SUBROUTINE make_chem_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_chem(nx1, nx2, nz1, nz_in, nz_out, nch, z_in, z_out, &
                 data_in, data_out, extrapolate)

    ! Interpolates columns of chemistry data from one set of height surfaces to
    ! another
 
    INTEGER, INTENT(IN)                :: nx1, nx2
    INTEGER, INTENT(IN)                :: nz1
    INTEGER, INTENT(IN)                :: nz_in
    INTEGER, INTENT(IN)                :: nz_out
    INTEGER, INTENT(IN)                :: nch
    REAL, INTENT(IN)                   :: z_in (nz1:nz_in,nx1:nx2 )
    REAL, INTENT(IN)                   :: z_out(nz1:nz_out,nx1:nx2)
    REAL, INTENT(IN)                   :: data_in (lo,nz1:nz_in,nx1:nx2)
    REAL, INTENT(INOUT)                :: data_out(lo,nz1:nz_out,nx1:nx2)
    LOGICAL, INTENT(IN)                :: extrapolate
    INTEGER                            :: i,l
    INTEGER                            :: k,kk
    REAL                               :: desired_z
    REAL                               :: dvaldz
    REAL                               :: wgt0
  
!   Loop over the number of chemical species
    chem_loop: DO l = 1, lo
!czhao     data_out(nx1:nx2,nz1:nz_out,ny1:ny2,l) = -99999.9
    data_out(l,:,:) = -99999.9

        DO i = nx1, nx2

          output_loop: DO k = nz1, nz_out
            desired_z = z_out(k,i)
            IF (desired_z .LT. z_in(1,i)) THEN

              IF ((desired_z - z_in(1,i)).LT. 0.0001) THEN
                 data_out(l,k,i) = data_in(l,1,i)
              ELSE
                IF (extrapolate) THEN
                  ! Extrapolate downward because desired height level is below
                  ! the lowest level in our input data.  Extrapolate using simple
                  ! 1st derivative of value with respect to height for the bottom 2
                  ! input layers.
  
                  ! Add a check to make sure we are not using the gradient of 
                  ! a very thin layer
                  IF ( (z_in(1,i) - z_in(2,i)) .GT. 0.001) THEN
                    dvaldz = (data_in(l,1,i) - data_in(l,2,i)) / &
                              (z_in(1,i)  - z_in(2,i) )
                  ELSE
                    dvaldz = (data_in(l,1,i) - data_in(l,3,i)) / &
                              (z_in(1,i)  - z_in(3,i) )
                  ENDIF
                  data_out(l,k,i) = MAX( data_in(l,1,i) + &
                                dvaldz * (desired_z-z_in(1,i)), 0.)
                ELSE
                  data_out(l,k,i) = data_in(l,1,i)
                ENDIF
              ENDIF
            ELSE IF (desired_z .GT. z_in(nz_in,i)) THEN
              IF ( (z_in(nz_in,i) - desired_z) .LT. 0.0001) THEN
                 data_out(l,k,i) = data_in(l,nz_in,i)
              ELSE
                IF (extrapolate) THEN
                  ! Extrapolate upward
                  IF ( (z_in(nz_in-1,i)-z_in(nz_in,i)) .GT. 0.0005) THEN
                    dvaldz = (data_in(l,nz_in,i) - data_in(l,nz_in-1,i)) / &
                               (z_in(nz_in,i)  - z_in(nz_in-1,i))
                  ELSE
                    dvaldz = (data_in(l,nz_in,i) - data_in(l,nz_in-2,i)) / &
                               (z_in(nz_in,i)  - z_in(nz_in-2,i)) 
                  ENDIF
                  data_out(l,k,i) =  MAX( data_in(l,nz_in,i) + &
                           dvaldz * (desired_z-z_in(nz_in,i)), 0.)
                ELSE
                  data_out(l,k,i) = data_in(l,nz_in,i)
                ENDIF
              ENDIF
            ELSE
              ! We can trap between two levels and linearly interpolate
  
              input_loop:  DO kk = 1, nz_in-1
                IF (desired_z .EQ. z_in(kk,i) )THEN
                  data_out(l,k,i) = data_in(l,kk,i)
                  EXIT input_loop
                ELSE IF (desired_z .EQ. z_in(kk+1,i) )THEN
                  data_out(l,k,i) = data_in(l,kk+1,i)
                  EXIT input_loop
                ELSE IF ( (desired_z .GT. z_in(kk,i)) .AND. &
                          (desired_z .LT. z_in(kk+1,i)) ) THEN
                  wgt0 = (desired_z - z_in(kk+1,i)) /  (z_in(kk,i)-z_in(kk+1,i))
                  data_out(l,k,i) = MAX( wgt0*data_in(l,kk,i) + (1.-wgt0)*data_in(l,kk+1,i), 0.)
                  EXIT input_loop
                ENDIF        
              ENDDO input_loop

            ENDIF
          ENDDO output_loop
        ENDDO 
    ENDDO chem_loop

    RETURN
  END SUBROUTINE vinterp_chem

END MODULE iamas_atmchem_gas_init
