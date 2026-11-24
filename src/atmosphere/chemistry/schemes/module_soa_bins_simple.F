!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module for simple secondary organic aerosol (SOA) formation for 
! size-resolved aerosols
! * Author: Xueyin Ruan from AEMOL USTC, ftrxyyg@mail.ustc.edu.cn
! Last update: xyruan 1/2026
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

module module_soa_bins_simple
  use iamas_atmchem_vars
  use module_data_mosaic_therm, only: nbin_a

  implicit none
  private
  public :: simple_SOA_driver

  !=====================================================================
  ! SIMPLE SOA MODULE (GROUPED-VOC, SIZE-RESOLVED)
  !
  ! PURPOSE
  !   Provide a compact, robust SOA (Secondary Organic Aerosol) source
  !   for a sectional aerosol model (NB=4 or 8 bins), using:
  !     (1) a kinetic OH-oxidation proxy: dVOC ~ kOH * [OH] * [VOC] * dt
  !     (2) a mass yield Y to convert oxidized VOC to condensable organics (SOG)
  !     (3) a simple absorptive partition factor: OA / (OA + Cstar)
  !     (4) bin distribution proportional to existing OA (POA + SOA), with a
  !         lognormal fallback when OA is near zero.
  !
  !
  ! UNITS
  !   - rho_air: kg/m3
  !   - OH and grouped VOCs: ppmv
  !   - POA and SOA in each size bin: ug/kg-air
  !   - Bin edges dlo/dhi: um (micrometers)
  !
  ! DEFAULT PARAMETERS (HIGH UNCERTAINTY)
  !   The parameters below are intended as reasonable starting points:
  !
  !   Kinetic OH rate constants (effective group values):
  !     kOH_aro  ~ 2e-11 cm3/molec/s
  !     kOH_terp ~ 1e-10
  !     kOH_isop ~ 1e-10
  !     kOH_oth  ~ 5e-12
  !
  !   Mass yields Y (dimensionless):
  !     Y_aro  = 0.15
  !     Y_terp = 0.25
  !     Y_isop = 0.08
  !     Y_oth  = 0.05
  !
  !   Effective condensable product molecular weight:
  !     MW_SOG = 200 g/mol
  !
  !   Absorptive partition strength:
  !     Cstar_ugkg = 25 ug/kg-air 
  !     (Typical useful range 10–50; larger values reduce partitioning.)
  !
  ! These defaults can be adjusted later if observational constraints become available.
  !
  !
  !   dSOA_total = dSOG * OA_total / (OA_total + Cstar_ugkg)
  !   (optional)
  !   where OA_total = sum_i (POA_i + SOA_i) [ug/kg]
  !
  !   Bin distribution: w_i ∝ (POA_i + SOA_i), fallback to a lognormal
  !   profile centered in the accumulation mode if OA_total is nearly zero.
  !
  !=====================================================================

  ! --- constants ---
  real(kind=8), parameter :: MW_AIR   = 28.97d-3        ! [kg/mol]
  real(kind=8), parameter :: AVOG     = 6.02214076d23   ! [molec/mol]
  real(kind=8), parameter :: PPMV2VMR = 1.0d-6          ! ppmv -> mol/mol

  ! --- default kinetic parameters (effective group values) ---
  real(kind=8), parameter :: KOH_ARO_DEF  = 2.0d-11      ! [cm3/molec/s]
  real(kind=8), parameter :: KOH_TERP_DEF = 1.0d-10
  real(kind=8), parameter :: KOH_ISOP_DEF = 1.0d-10
  real(kind=8), parameter :: KOH_OTH_DEF  = 5.0d-12

  ! --- default yields (dimensionless) ---
  real(kind=8), parameter :: Y_ARO_DEF  = 0.15d0
  real(kind=8), parameter :: Y_TERP_DEF = 0.25d0
  real(kind=8), parameter :: Y_ISOP_DEF = 0.08d0
  real(kind=8), parameter :: Y_OTH_DEF  = 0.05d0

  ! --- default effective product MW ---
  real(kind=8), parameter :: MW_SOG_DEF = 200.0d0        ! [g/mol]

  ! --- default partition strength ---
  real(kind=8), parameter :: CSTAR_DEF  = 25.0d0         ! [ug/kg-air]

  ! --- default control ---
  ! apply OA/(OA+Cstar) if true; else all SOG becomes SOA
  logical :: use_partition_DEF = .true.

contains

!***********************************************************************
  subroutine simple_SOA_driver (dtstepc, t_phy, rho_phy, chem,  &
                                its,ite, kts,kte              )

  use module_data_binaer_asect, only: dlo_sect, dhi_sect, &
                                      ai_phase,           &
                                      lptr_soa_aer

  implicit none
  real, intent(in) :: dtstepc
  real, dimension( kts:kte, its:ite), intent(in) :: t_phy, rho_phy
  real, dimension( 1:num_chem, kts:kte, its:ite ), intent(inout) :: chem
  integer, intent(in) :: its,ite, kts,kte

  ! local vars
  integer:: iphase,itype
  integer:: i,k
  real(kind=8) :: T, rho_air
  real(kind=8), dimension(1:nbin_a) :: dlo, dhi
  real(kind=8) :: OH_ppmv, VOC_ARO_ppmv, VOC_TERP_ppmv, VOC_ISOP_ppmv, VOC_OTH_ppmv
  real(kind=8), dimension(1:nbin_a) :: POA_ugkg, SOA_ugkg

  itype = 1
  iphase = ai_phase

  do i = its, ite
  do k = kts, kte

   T = t_phy(k,i)              ! K
   rho_air = rho_phy(k,i)      ! kg/m3 
   dlo(1:nbin_a) = dlo_sect(1:nbin_a,itype)*10000 !cm-->um   
   dhi(1:nbin_a) = dhi_sect(1:nbin_a,itype)*10000 !cm-->um

   OH_ppmv = chem(p_ho,k,i)

  ! --- build VOC groups (ppmv) ---
   ! VOC_ARO  = ARO1 + ARO2 + PHEN + NPHE + CRES + BALD + BACL + PBZN
   VOC_ARO_ppmv = chem(p_aro1,k,i) + chem(p_aro2,k,i) + &
                  chem(p_phen,k,i) + chem(p_nphe,k,i) + &
                  chem(p_cres,k,i) + chem(p_bald,k,i) + &
                  chem(p_bacl,k,i) + chem(p_pbzn,k,i) 
   
   ! VOC_TERP = TERP + SESQ
   VOC_TERP_ppmv = chem(p_terp,k,i) + chem(p_sesq,k,i)

   ! VOC_ISOP = ISOPRENE   ! or + ISOPROD + MVK + METHACRO + PROD2
   VOC_ISOP_ppmv =  chem(p_isoprene,k,i)

   ! VOC_OTH  = ALK3 + ALK4 + ALK5 + OLE1 + OLE2 + ETHENE + C2H6 + C2H2 + C3H6
   VOC_OTH_ppmv = chem(p_alk3,k,i) + chem(p_alk4,k,i) + &
                  chem(p_alk5,k,i) + chem(p_ole1,k,i) + &
                  chem(p_ole2,k,i) + chem(p_ethene,k,i) + &
                  chem(p_c2h6,k,i) + chem(p_c2h2,k,i) + &
                  chem(p_c3h6,k,i)

   ! POA
   POA_ugkg(1:nbin_a) = chem(lptr_oc_aer(1:nbin_a,itype,iphase),k,i)

   ! SOA
   SOA_ugkg(1:nbin_a) = chem(lptr_soa_aer(1:nbin_a,itype,iphase),k,i)

   ! --- call SOA; uses internal default parameters ---
   call soa_step_bins_groups (dtstepc = dtstepc, T = T, rho_air = rho_air, &
                              dlo = dlo, dhi = dhi, &
                              OH_ppmv = OH_ppmv, &
                              VOC_ARO_ppmv = VOC_ARO_ppmv, VOC_TERP_ppmv = VOC_TERP_ppmv, &
                              VOC_ISOP_ppmv = VOC_ISOP_ppmv, VOC_OTH_ppmv = VOC_OTH_ppmv, &
                              POA_ugkg = POA_ugkg, SOA_ugkg = SOA_ugkg  )


   chem(lptr_soa_aer(1:nbin_a,itype,iphase),k,i) = SOA_ugkg(1:nbin_a)

  end do ! k
  end do ! i


  end subroutine simple_SOA_driver

  !---------------------------------------------------------------------
  ! SUBROUTINE: soa_step_bins_groups
  ! PURPOSE
  !   Compute SOA production and update size-bin SOA (ug/kg-air) 
  !   using grouped VOCs and OH.
  !
  ! INPUTS
  !   dtstepc, T      : timestep [s], temperature [K] (T currently unused)
  !   rho_air        : air density [kg/m3]
  !   dlo(:), dhi(:) : bin edges [um], length NB=4 or 8
  !   OH_ppmv        : OH [ppmv]
  !   VOC_*_ppmv     : grouped VOC precursors [ppmv]
  !   POA_ugkg(:)    : primary OA per bin [ug/kg-air]
  !
  ! INPUT/OUTPUT
  !   SOA_ugkg(:)    : secondary OA per bin [ug/kg-air]
  !
  ! OPTIONAL PARAMETERS (all have defaults if not provided)
  !   kOH_*        : [cm3/molec/s]
  !   Y_*          : [-]
  !   MW_SOG       : [g/mol]
  !   Cstar_ugkg   : partition strength [ug/kg-air]
  !   use_partition: apply OA/(OA+Cstar) if true; else all dSOG becomes SOA
  !---------------------------------------------------------------------
  subroutine soa_step_bins_groups(dtstepc, T, rho_air, dlo, dhi, &
                                  OH_ppmv, &
                                  VOC_ARO_ppmv, VOC_TERP_ppmv, &
                                  VOC_ISOP_ppmv, VOC_OTH_ppmv, &
                                  POA_ugkg, SOA_ugkg, &
                                  kOH_aro, kOH_terp, kOH_isop, kOH_oth, &
                                  Y_aro,  Y_terp,  Y_isop,  Y_oth, &
                                  MW_SOG_gmol, Cstar_ugkg,         &
                                  use_partition )
    implicit none
    real(kind=8), intent(in)    :: dtstepc, T, rho_air
    real(kind=8), intent(in)    :: dlo(:), dhi(:)
    real(kind=8), intent(in)    :: OH_ppmv
    real(kind=8), intent(in)    :: VOC_ARO_ppmv, VOC_TERP_ppmv, &
                                   VOC_ISOP_ppmv, VOC_OTH_ppmv
    real(kind=8), intent(in)    :: POA_ugkg(:)
    real(kind=8), intent(inout) :: SOA_ugkg(:)

    ! Optional overrides
    real(kind=8), intent(in), optional :: kOH_aro, kOH_terp, kOH_isop, kOH_oth
    real(kind=8), intent(in), optional :: Y_aro,  Y_terp,  Y_isop,  Y_oth
    real(kind=8), intent(in), optional :: MW_SOG_gmol, Cstar_ugkg
    logical,      intent(in), optional :: use_partition

    integer :: ibin
    real(kind=8) :: Nair
    real(kind=8) :: OH_c, VOCa_c, VOCt_c, VOCi_c, VOCo_c
    real(kind=8) :: dVOCa_molec, dVOCt_molec, dVOCi_molec, dVOCo_molec
    real(kind=8) :: dVOCa_molkg, dVOCt_molkg, dVOCi_molkg, dVOCo_molkg
    real(kind=8) :: dSOG_ugkg, dSOA_total
    real(kind=8) :: OA_sum, frac
    real(kind=8) :: w(size(dlo))

    ! Effective parameters (defaults unless overridden)
    real(kind=8) :: kA, kT, kI, kO
    real(kind=8) :: yA, yT, yI, yO
    real(kind=8) :: MWs, cS
    logical      :: use_p

    kA  = KOH_ARO_DEF
    kT  = KOH_TERP_DEF
    kI  = KOH_ISOP_DEF
    kO  = KOH_OTH_DEF
    yA  = Y_ARO_DEF
    yT  = Y_TERP_DEF
    yI  = Y_ISOP_DEF
    yO  = Y_OTH_DEF
    MWs = MW_SOG_DEF
    cS  = CSTAR_DEF
    use_p = use_partition_DEF

    if (present(kOH_aro))   kA = kOH_aro
    if (present(kOH_terp))  kT = kOH_terp
    if (present(kOH_isop))  kI = kOH_isop
    if (present(kOH_oth))   kO = kOH_oth
    if (present(Y_aro))     yA = Y_aro
    if (present(Y_terp))    yT = Y_terp
    if (present(Y_isop))    yI = Y_isop
    if (present(Y_oth))     yO = Y_oth
    if (present(MW_SOG_gmol)) MWs = MW_SOG_gmol
    if (present(Cstar_ugkg))  cS  = Cstar_ugkg
    if (present(use_partition))  use_p = use_partition

    !---------------------------------------------------------------
    ! (1) Convert ppmv -> molec/cm3 using air density
    !---------------------------------------------------------------
    Nair  = air_number_density_cm3(rho_air)

    OH_c  = ppmv_to_molec_cm3(OH_ppmv,       Nair)
    VOCa_c= ppmv_to_molec_cm3(VOC_ARO_ppmv,  Nair)
    VOCt_c= ppmv_to_molec_cm3(VOC_TERP_ppmv, Nair)
    VOCi_c= ppmv_to_molec_cm3(VOC_ISOP_ppmv, Nair)
    VOCo_c= ppmv_to_molec_cm3(VOC_OTH_ppmv,  Nair)

    !---------------------------------------------------------------
    ! (2) Kinetic consumption of VOC by OH (molec/cm3)
    !     dVOC = kOH * [OH] * [VOC] * dtstepc
    !---------------------------------------------------------------
    dVOCa_molec = max0(kA) * OH_c * VOCa_c * max0(dtstepc)
    dVOCt_molec = max0(kT) * OH_c * VOCt_c * max0(dtstepc)
    dVOCi_molec = max0(kI) * OH_c * VOCi_c * max0(dtstepc)
    dVOCo_molec = max0(kO) * OH_c * VOCo_c * max0(dtstepc)

    !---------------------------------------------------------------
    ! (3) Convert dVOC (molec/cm3) -> mol/kg-air
    !---------------------------------------------------------------
    dVOCa_molkg = dVOC_molec_to_molkg(dVOCa_molec, Nair)
    dVOCt_molkg = dVOC_molec_to_molkg(dVOCt_molec, Nair)
    dVOCi_molkg = dVOC_molec_to_molkg(dVOCi_molec, Nair)
    dVOCo_molkg = dVOC_molec_to_molkg(dVOCo_molec, Nair)

    !---------------------------------------------------------------
    ! (4) Convert oxidized VOC to condensable organics (SOG) mass [ug/kg-air]
    !     dSOG = sum( dVOC_molkg * Y * MW_SOG[g/mol] ) * 1e6 [ug/g]
    !---------------------------------------------------------------
    dSOG_ugkg = ( dVOCa_molkg * max0(yA) + &
                 dVOCt_molkg * max0(yT) + &
                 dVOCi_molkg * max0(yI) + &
                 dVOCo_molkg * max0(yO) ) * max0(MWs) * 1.0d6
    dSOG_ugkg = max0(dSOG_ugkg)

    !---------------------------------------------------------------
    ! (5) Optional absorptive partitioning: OA/(OA + Cstar)
    !---------------------------------------------------------------
    if (use_p) then
      OA_sum = 0d0
      do ibin = 1, nbin_a
        OA_sum = OA_sum + max0(POA_ugkg(ibin) + SOA_ugkg(ibin))
      end do
      if (OA_sum > 0d0) then
        frac = OA_sum / (OA_sum + max0(cS))
      else
        frac = 0d0
      end if
      dSOA_total = dSOG_ugkg * frac
    else
      dSOA_total = dSOG_ugkg
    end if

    !---------------------------------------------------------------
    ! (6) Bin distribution weights
    !     - Primary: proportional to existing OA_i = POA_i + SOA_i
    !     - Fallback: lognormal (accumulation-centered) if OA is tiny
    !---------------------------------------------------------------
    OA_sum = 0d0
    do ibin = 1, nbin_a
      w(ibin) = max0(POA_ugkg(ibin) + SOA_ugkg(ibin))
      OA_sum = OA_sum + w(ibin)
    end do

    if (OA_sum <= 1d-12) then
      call lognorm_weights(dlo, dhi, w, d0=0.35d0, sigln=0.8d0)
    else
      call normalize_weights(w)
    end if

    !---------------------------------------------------------------
    ! (7) Update SOA per bin
    !---------------------------------------------------------------
    if (dSOA_total > 0d0) then
      do ibin = 1, nbin_a
        SOA_ugkg(ibin) = max0(SOA_ugkg(ibin) + w(ibin) * dSOA_total)
      end do
    end if

  end subroutine soa_step_bins_groups

  !---------------------------------------------------------------------
  ! FUNCTION: max0
  ! PURPOSE
  !   Clamp negatives to 0 for numerical safety.
  !---------------------------------------------------------------------
  real(kind=8) function max0(x)
    implicit none
    real(kind=8), intent(in) :: x
    if (x < 0d0) then
      max0 = 0d0
    else
      max0 = x
    end if
  end function max0

  !---------------------------------------------------------------------
  ! SUBROUTINE: normalize_weights
  ! PURPOSE
  !   Normalize a non-negative weight vector so that sum(w)=1.
  !   If the sum is zero, assign all weight to the first bin.
  !---------------------------------------------------------------------
  subroutine normalize_weights(w)
    implicit none
    real(kind=8), intent(inout) :: w(:)
    integer :: ibin
    real(kind=8) :: s
    s = 0d0
    do ibin = 1, nbin_a
      w(ibin) = max0(w(ibin))
      s = s + w(ibin)
    end do
    if (s > 0d0) then
      do ibin = 1, nbin_a
        w(ibin) = w(ibin) / s
      end do
    else
      w = 0d0
      w(1) = 1d0
    end if
  end subroutine normalize_weights

  !---------------------------------------------------------------------
  ! SUBROUTINE: lognorm_weights
  ! PURPOSE
  !   Fallback lognormal weights for distributing new SOA when OA is
  !   near zero (OA-based weights become ill-defined).
  !
  ! INPUT
  !   dlo, dhi : bin edges [um]
  !   d0       : center diameter [um]
  !   sigln    : width in ln-space
  !
  ! OUTPUT
  !   w(:)     : normalized weights
  !---------------------------------------------------------------------
  subroutine lognorm_weights(dlo, dhi, w, d0, sigln)
    implicit none
    real(kind=8), intent(in)  :: dlo(:), dhi(:)
    real(kind=8), intent(out) :: w(:)
    real(kind=8), intent(in)  :: d0, sigln
    integer :: ibin
    real(kind=8) :: dmid, x

    do ibin = 1, nbin_a
      dmid = sqrt(dlo(ibin) * dhi(ibin))
      x    = log(dmid / d0)
      w(ibin) = exp( -0.5d0 * (x/sigln)**2 )
    end do
    call normalize_weights(w)
  end subroutine lognorm_weights

  !---------------------------------------------------------------------
  ! FUNCTION: air_number_density_cm3
  ! PURPOSE
  !   Compute air number density Nair [molec/cm3] from air density rho_air [kg/m3].
  !
  !   mol/m3  = rho_air / MW_AIR
  !   molec/m3 = mol/m3 * Avog
  !   molec/cm3 = molec/m3 / 1e6
  !---------------------------------------------------------------------
  real(kind=8) function air_number_density_cm3(rho_air)
    implicit none
    real(kind=8), intent(in) :: rho_air
    air_number_density_cm3 = max0(rho_air) / MW_AIR * AVOG / 1.0d6
  end function air_number_density_cm3

  !---------------------------------------------------------------------
  ! FUNCTION: ppmv_to_molec_cm3
  ! PURPOSE
  !   Convert ppmv to concentration [molec/cm3] using Nair [molec/cm3].
  !---------------------------------------------------------------------
  real(kind=8) function ppmv_to_molec_cm3(ppmv, Nair)
    implicit none
    real(kind=8), intent(in) :: ppmv, Nair
    ppmv_to_molec_cm3 = max0(ppmv) * PPMV2VMR * max0(Nair)
  end function ppmv_to_molec_cm3

  !---------------------------------------------------------------------
  ! FUNCTION: dVOC_molec_to_molkg
  ! PURPOSE
  !   Convert a change in concentration (molec/cm3) to mol/kg-air.
  !
  ! LOGIC
  !   dVOC_molec_cm3 / Nair  -> delta mixing ratio (mol/mol)
  !   (mol/mol) / MW_AIR     -> mol/kg-air
  !---------------------------------------------------------------------
  real(kind=8) function dVOC_molec_to_molkg(dVOC_molec_cm3, Nair)
    implicit none
    real(kind=8), intent(in) :: dVOC_molec_cm3, Nair
    if (Nair > 0d0) then
      dVOC_molec_to_molkg = max0(dVOC_molec_cm3) / Nair / MW_AIR
    else
      dVOC_molec_to_molkg = 0d0
    end if
  end function dVOC_molec_to_molkg


end module module_soa_bins_simple

