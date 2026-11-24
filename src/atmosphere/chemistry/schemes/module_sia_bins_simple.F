!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module for simple secondary inorganic aerosol (SIA) formation for 
! size-resolved aerosols
! * Author: Xueyin Ruan from AEMOL USTC, ftrxyyg@mail.ustc.edu.cn
! Last update: xyruan 12/2025
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

module module_sia_bins_simple 
  use iamas_atmchem_vars, only: num_chem, p_h2so4, p_hno3, p_nh3
  use module_data_mosaic_therm, only: nbin_a

  implicit none
  private
  public :: simple_SIA_driver

  !=====================================================================
  ! PURPOSE
  !   A minimal secondary inorganic aerosol (SIA) module for size-resolved
  !   aerosol mass bins. It updates sulfate (SO4), nitrate (NO3), and
  !   ammonium (NH4) in each aerosol size bin using:
  !
  !     (1) Condensation of gas-phase H2SO4 (already computed in gas chemistry)
  !     (2) Ammonia neutralization of sulfate (2 NH3 per 1 SO4 in molar units)
  !     (3) Formation of ammonium nitrate: NH3 + HNO3 -> NH4NO3 (1:1)
  !     (4) Optional volatilization of NH4NO3 back to NH3 + HNO3 controlled by T/RH
  !
  ! DESIGN GOALS
  !   - Very short and robust (does not require pressure or air density)
  !   - Consistent with your units:
  !       Gas species: ppmv (molar mixing ratio in ppm)
  !       Aerosol species in bins: ug/kg-air (mass mixing ratio)
  !   - Automatically adapts to either 4-bin or 8-bin configurations by
  !     receiving bin edges dlo(:), dhi(:) with size(dlo)=NB.
  !
  ! IMPORTANT LIMITATIONS / ASSUMPTIONS
  !   - All H2SO4_g is condensed each call (fast irreversible condensation).
  !     If you need a finite condensation timescale, we can add tau_cond later.
  !   - Nitrate is treated only as NH4NO3. We do NOT include coarse-mode
  !     nitrate (e.g., HNO3 uptake on dust/sea-salt) in this simple scheme.
  !   - Size distribution of newly formed mass is parameterized, not
  !   prognostic:
  !     Positive increments are distributed via a lognormal weight
  !     centered in the accumulation mode.
  !
  !=====================================================================

  !--- constants for unit conversion (per kg of air) ---
  real(kind=8), parameter :: MW_AIR = 28.97d-3      ! [kg/mol] mean molar mass of dry air
  real(kind=8), parameter :: PPMV_TO_MOLFRAC = 1d-6 ! 1 ppmv = 1e-6 mol/mol

  !--- molar masses (g/mol) for aerosol ions (as carried in your aerosol tracers) ---
  real(kind=8), parameter :: MW_SO4 = 96.06d0
  real(kind=8), parameter :: MW_NO3 = 62.00d0
  real(kind=8), parameter :: MW_NH4 = 18.04d0

contains
!***********************************************************************
  subroutine simple_SIA_driver (dtstepc, t_phy, relhum, chem, &
                                do_volatilize,               &
                                its,ite, kts,kte              )

  use module_data_binaer_asect, only: dlo_sect, dhi_sect, &
                                      ai_phase,           &
                                      lptr_so4_aer,       &
                                      lptr_no3_aer,       &
                                      lptr_nh4_aer
  
  implicit none
  real, intent(in) :: dtstepc
  real, dimension( kts:kte, its:ite), intent(in) :: t_phy,relhum
  real, dimension( 1:num_chem, kts:kte, its:ite ), intent(inout) :: chem
  logical, intent(in) :: do_volatilize
  integer, intent(in) :: its,ite, kts,kte


  ! local vars
  integer:: iphase,itype
  integer:: i,k
  real(kind=8) :: T, RH
  real(kind=8), dimension(1:nbin_a) :: dlo, dhi
  real(kind=8) :: H2SO4_g_ppmv, HNO3_g_ppmv, NH3_g_ppmv
  real(kind=8), dimension(1:nbin_a) :: SO4_ugkg, NO3_ugkg, NH4_ugkg


  itype = 1
  iphase = ai_phase
  
  do i = its, ite
  do k = kts, kte

   T = t_phy(k,i)          ! K
   RH = relhum(k,i)/100.   ! (%)-->(0—1)
   dlo(1:nbin_a) = dlo_sect(1:nbin_a,itype)*10000 !cm-->um   
   dhi(1:nbin_a) = dhi_sect(1:nbin_a,itype)*10000 !cm-->um

   H2SO4_g_ppmv = chem(p_h2so4,k,i)
   HNO3_g_ppmv = chem(p_hno3,k,i)
   NH3_g_ppmv = chem(p_nh3,k,i)

   SO4_ugkg(1:nbin_a) = chem(lptr_so4_aer(1:nbin_a,itype,iphase),k,i)
   NO3_ugkg(1:nbin_a) = chem(lptr_no3_aer(1:nbin_a,itype,iphase),k,i)
   NH4_ugkg(1:nbin_a) = chem(lptr_nh4_aer(1:nbin_a,itype,iphase),k,i)

  call sia_step_bins(dtstepc, T, RH, dlo, dhi, &
                     H2SO4_g_ppmv, HNO3_g_ppmv, NH3_g_ppmv, &
                     SO4_ugkg, NO3_ugkg, NH4_ugkg, &
                     do_volatilize )
   
    chem(p_h2so4,k,i) = H2SO4_g_ppmv
    chem(p_hno3,k,i) = HNO3_g_ppmv
    chem(p_nh3,k,i) = NH3_g_ppmv

    chem(lptr_so4_aer(1:nbin_a,itype,iphase),k,i) = SO4_ugkg(1:nbin_a)
    chem(lptr_no3_aer(1:nbin_a,itype,iphase),k,i) = NO3_ugkg(1:nbin_a)
    chem(lptr_nh4_aer(1:nbin_a,itype,iphase),k,i) = NH4_ugkg(1:nbin_a)


  end do ! k
  end do ! i

  end subroutine simple_SIA_driver

  !---------------------------------------------------------------------
  ! SUBROUTINE: sia_step_bins
  ! PURPOSE
  !   Advance SIA (SO4/NO3/NH4) over one timestep for a size-resolved
  !   aerosol.
  !
  ! INPUT/OUTPUT (UNITS)
  !   Gas species (in/out): ppmv
  !     H2SO4_g_ppmv  : gas-phase sulfuric acid (already computed
  !     elsewhere)
  !     HNO3_g_ppmv   : gas-phase nitric acid
  !     NH3_g_ppmv    : gas-phase ammonia
  !
  !   Aerosol species (in/out): ug/kg-air per bin
  !     SO4_ugkg(:), NO3_ugkg(:), NH4_ugkg(:)
  !
  !   Bin edges:
  !     dlo(:), dhi(:) in um. size(dlo)=NB defines NB automatically (4
  !     or 8).
  !
  ! METEO
  !   T: K, RH: 0–1
  !
  ! LOGIC OVERVIEW
  !   (A) Convert gas ppmv -> mol/kg-air.
  !   (B) Convert aerosol ug/kg -> mol/kg (per bin and bulk).
  !   (C) Condense all H2SO4_g to sulfate (bulk increment).
  !   (D) Neutralize sulfate using NH3; then form NH4NO3 with remaining
  !   NH3+HNO3.
  !   (E) Optionally volatilize NH4NO3 based on fp(T,RH).
  !   (F) Apply resulting bulk increments back to bins:
  !       - Positive increments: distributed by lognormal weights w(i)
  !       - Negative increments (from volatilization): applied as
  !       proportional scaling
  !         to preserve existing bin shape.
  !
  ! WHY BULK THEN BIN?
  !   It keeps stoichiometry stable and code short, while still updating
  !   bins.
  !
  ! dt
  !   Present for interface consistency; current implementation does not
  !   explicitly use dt because H2SO4 is assumed to fully condense each
  !   call.
  !---------------------------------------------------------------------
  subroutine sia_step_bins(dt, T, RH, dlo, dhi, &
                           H2SO4_g_ppmv, HNO3_g_ppmv, NH3_g_ppmv, &
                           SO4_ugkg, NO3_ugkg, NH4_ugkg, &
                           do_volatilize)
    implicit none
    real(kind=8), intent(in)    :: dt, T, RH
    real(kind=8), intent(in)    :: dlo(:), dhi(:)
    real(kind=8), intent(inout) :: H2SO4_g_ppmv, HNO3_g_ppmv, NH3_g_ppmv
    real(kind=8), intent(inout) :: SO4_ugkg(:), NO3_ugkg(:), NH4_ugkg(:)
    logical,      intent(in)    :: do_volatilize

    integer :: ibin
    real(kind=8) :: w(size(dlo))

    real(kind=8), allocatable :: so4_bin_mol(:), no3_bin_mol(:), nh4_bin_mol(:)
    real(kind=8) :: so4_old, so4_new
    real(kind=8) :: no3_old, no3_new
    real(kind=8) :: nh4_old, nh4_new

    real(kind=8) :: h2so4_molkg, hno3_molkg, nh3_molkg
    real(kind=8) :: nh3_bulk, hno3_bulk, nh4_bulk, no3_bulk

    real(kind=8) :: add_so4, need_nh4_for_s, use_nh3_for_s
    real(kind=8) :: form_no3, fp, no3_pre, evap_no3
    real(kind=8), parameter :: tiny = 1d-20

    allocate(so4_bin_mol(nbin_a), no3_bin_mol(nbin_a), nh4_bin_mol(nbin_a))

    ! Weight profile used for distributing newly formed secondary mass
    call compute_lognorm_weights(dlo, dhi, w, d0=0.35d0, sigln=0.8d0)

    ! Gas conversion: ppmv -> mol/kg-air
    h2so4_molkg = ppmv_to_molkg(H2SO4_g_ppmv)
    hno3_molkg  = ppmv_to_molkg(HNO3_g_ppmv)
    nh3_molkg   = ppmv_to_molkg(NH3_g_ppmv)

    ! Aerosol conversion: ug/kg-air -> mol/kg-air
    so4_old = 0d0; no3_old = 0d0; nh4_old = 0d0
    do ibin = 1, nbin_a
      so4_bin_mol(ibin) = ugkg_to_molkg(SO4_ugkg(ibin), MW_SO4)
      no3_bin_mol(ibin) = ugkg_to_molkg(NO3_ugkg(ibin), MW_NO3)
      nh4_bin_mol(ibin) = ugkg_to_molkg(NH4_ugkg(ibin), MW_NH4)
      so4_old = so4_old + so4_bin_mol(ibin)
      no3_old = no3_old + no3_bin_mol(ibin)
      nh4_old = nh4_old + nh4_bin_mol(ibin)
    end do

    !-------------------------------------------------------------------
    ! (1) Condense all H2SO4_g -> sulfate
    !     1 mol H2SO4 is treated as producing 1 mol SO4 in aerosol
    !     tracers.
    !-------------------------------------------------------------------
    add_so4      = max0(h2so4_molkg)
    so4_new      = so4_old + add_so4
    h2so4_molkg  = 0d0   ! fully condensed

    !-------------------------------------------------------------------
    ! (2) Bulk stoichiometry for NH3/HNO3 with aerosol
    ! sulfate/nitrate/ammonium
    !-------------------------------------------------------------------
    nh3_bulk = nh3_molkg
    hno3_bulk= hno3_molkg
    nh4_bulk = nh4_old
    no3_bulk = no3_old

    ! (2a) Neutralize sulfate first: 2 NH4 needed per 1 SO4 (molar)
    need_nh4_for_s = max0(2d0 * so4_new - nh4_bulk)
    use_nh3_for_s  = min(nh3_bulk, need_nh4_for_s)
    nh4_bulk = nh4_bulk + use_nh3_for_s
    nh3_bulk = nh3_bulk - use_nh3_for_s

    ! (2b) Form ammonium nitrate: NH3 + HNO3 -> NH4NO3
    form_no3  = min(nh3_bulk, hno3_bulk)
    no3_pre   = no3_bulk + form_no3
    nh4_bulk  = nh4_bulk + form_no3
    nh3_bulk  = nh3_bulk - form_no3
    hno3_bulk = hno3_bulk - form_no3

    ! (2c) Optional volatilization of NH4NO3
    if (do_volatilize) then
      fp        = fp_nh4no3(T, RH, 298d0, 0.5d0, 0.15d0, 4d0)
      no3_new   = fp * no3_pre
      evap_no3  = no3_pre - no3_new
      nh4_new   = nh4_bulk - evap_no3

      nh3_bulk  = nh3_bulk  + evap_no3
      hno3_bulk = hno3_bulk + evap_no3

      no3_bulk  = max0(no3_new)
      nh4_bulk  = max0(nh4_new)
    else
      no3_bulk  = no3_pre
      nh4_bulk  = nh4_bulk
    end if

    ! Final bulk after chemistry
    so4_new = max0(so4_new)
    no3_new = max0(no3_bulk)
    nh4_new = max0(nh4_bulk)

    !-------------------------------------------------------------------
    ! (3) Convert bulk changes into bin updates
    !     - Sulfate increment is always non-negative here.
    !     - Nitrate/ammonium can be negative if volatilization
    !     dominates.
    !-------------------------------------------------------------------
    add_so4 = max0(so4_new - so4_old)
    form_no3 = no3_new - no3_old
    use_nh3_for_s = nh4_new - nh4_old   ! net ammonium change (includes
                                        ! both sulfate and nitrate pairing)

    ! (3a) Sulfate: distribute positive increment by w
    if (add_so4 > 0d0) then
      do ibin = 1, nbin_a
        so4_bin_mol(ibin) = so4_bin_mol(ibin) + w(ibin) * add_so4
      end do
    end if

    ! (3b) Nitrate: positive -> distribute by w; negative -> scale down
    ! proportionally
    if (form_no3 > 0d0) then
      do ibin = 1, nbin_a
        no3_bin_mol(ibin) = no3_bin_mol(ibin) + w(ibin) * form_no3
      end do
    else if (form_no3 < 0d0 .and. no3_old > tiny) then
      do ibin = 1, nbin_a
        no3_bin_mol(ibin) = no3_bin_mol(ibin) * (no3_new / no3_old)
      end do
    end if

    ! (3c) Ammonium: positive -> distribute by w; negative -> scale down
    ! proportionally
    if (use_nh3_for_s > 0d0) then
      do ibin = 1, nbin_a
        nh4_bin_mol(ibin) = nh4_bin_mol(ibin) + w(ibin) * use_nh3_for_s
      end do
    else if (use_nh3_for_s < 0d0 .and. nh4_old > tiny) then
      do ibin = 1, nbin_a
        nh4_bin_mol(ibin) = nh4_bin_mol(ibin) * (nh4_new / nh4_old)
      end do
    end if

    !-------------------------------------------------------------------
    ! (4) Convert updated per-bin mol/kg back to ug/kg and write back
    ! gases
    !-------------------------------------------------------------------
    do ibin = 1, nbin_a
      SO4_ugkg(ibin) = molkg_to_ugkg(so4_bin_mol(ibin), MW_SO4)
      NO3_ugkg(ibin) = molkg_to_ugkg(no3_bin_mol(ibin), MW_NO3)
      NH4_ugkg(ibin) = molkg_to_ugkg(nh4_bin_mol(ibin), MW_NH4)
      SO4_ugkg(ibin) = max0(SO4_ugkg(ibin))
      NO3_ugkg(ibin) = max0(NO3_ugkg(ibin))
      NH4_ugkg(ibin) = max0(NH4_ugkg(ibin))
    end do

    H2SO4_g_ppmv = molkg_to_ppmv(h2so4_molkg)   ! will be ~0 by design
    HNO3_g_ppmv  = molkg_to_ppmv(hno3_bulk)
    NH3_g_ppmv   = molkg_to_ppmv(nh3_bulk)

    deallocate(so4_bin_mol, no3_bin_mol, nh4_bin_mol)

  end subroutine sia_step_bins

  !---------------------------------------------------------------------
  ! FUNCTION: max0
  ! PURPOSE
  !   Numerical safety helper: clamp any negative value to 0.
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
  ! FUNCTION: clamp01
  ! PURPOSE
  !   Clamp a value into [0, 1]. Used for partition fractions.
  !---------------------------------------------------------------------
  real(kind=8) function clamp01(x)
    implicit none
    real(kind=8), intent(in) :: x
    if (x < 0d0) then
      clamp01 = 0d0
    else if (x > 1d0) then
      clamp01 = 1d0
    else
      clamp01 = x
    end if
  end function clamp01

  !---------------------------------------------------------------------
  ! FUNCTION: ppmv_to_molkg
  ! PURPOSE
  !   Convert gas mixing ratio in ppmv (mol/mol * 1e6) to mol/kg-air.
  !
  ! WHY THIS WORKS WITHOUT PRESSURE/DENSITY
  !   ppmv is essentially mol/mol. For a dry-air basis:
  !     mol species per kg air = (mol/mol) / (kg/mol air)
  !
  ! INPUT
  !   ppmv : [ppmv]
  !
  ! OUTPUT
  !   mol/kg-air
  !---------------------------------------------------------------------
  real(kind=8) function ppmv_to_molkg(ppmv)
    implicit none
    real(kind=8), intent(in) :: ppmv
    real(kind=8) :: vmr
    vmr = max0(ppmv) * PPMV_TO_MOLFRAC   ! mol/mol
    ppmv_to_molkg = vmr / MW_AIR         ! mol/kg-air
  end function ppmv_to_molkg

  !---------------------------------------------------------------------
  ! FUNCTION: molkg_to_ppmv
  ! PURPOSE
  !   Convert mol/kg-air back to ppmv.
  !---------------------------------------------------------------------
  real(kind=8) function molkg_to_ppmv(molkg)
    implicit none
    real(kind=8), intent(in) :: molkg
    real(kind=8) :: vmr
    vmr = max0(molkg) * MW_AIR           ! mol/mol
    molkg_to_ppmv = vmr / PPMV_TO_MOLFRAC
  end function molkg_to_ppmv

  !---------------------------------------------------------------------
  ! FUNCTION: ugkg_to_molkg
  ! PURPOSE
  !   Convert aerosol mass mixing ratio (ug/kg-air) to mol/kg-air using
  !   MW.
  !
  ! INPUT
  !   ugkg    : [ug/kg-air]
  !   MW_gmol : [g/mol]
  !
  ! OUTPUT
  !   mol/kg-air
  !---------------------------------------------------------------------
  real(kind=8) function ugkg_to_molkg(ugkg, MW_gmol)
    implicit none
    real(kind=8), intent(in) :: ugkg, MW_gmol
    ugkg_to_molkg = max0(ugkg) * 1d-6 / MW_gmol   ! ug/kg -> g/kg ->mol/kg
  end function ugkg_to_molkg

  !---------------------------------------------------------------------
  ! FUNCTION: molkg_to_ugkg
  ! PURPOSE
  !   Convert mol/kg-air to aerosol mass mixing ratio (ug/kg-air).
  !---------------------------------------------------------------------
  real(kind=8) function molkg_to_ugkg(molkg, MW_gmol)
    implicit none
    real(kind=8), intent(in) :: molkg, MW_gmol
    molkg_to_ugkg = max0(molkg) * MW_gmol * 1d6   ! mol/kg -> g/kg ->ug/kg
  end function molkg_to_ugkg

  !---------------------------------------------------------------------
  ! SUBROUTINE: compute_lognorm_weights
  ! PURPOSE
  !   Construct a normalized weight vector w(i) over size bins for
  !   distributing newly formed secondary mass (positive increments). 
  !   We use a lognormal-like shape centered at diameter d0 (um) with 
  !   width sigln (ln-space std dev).
  !
  ! DESIGN CHOICE
  !   For SIA (sulfate/nitrate/ammonium), new mass typically grows the
  !   accumulation mode (~0.1–1 um). This weight is a compact
  !   parameterization.
  !
  ! INPUT
  !   dlo(:), dhi(:): lower/upper bin diameters in um
  !   d0            : geometric center (um), e.g. 0.35
  !   sigln         : width in ln-space, e.g. 0.8
  !
  ! OUTPUT
  !   w(:)          : normalized weights summing to 1
  !---------------------------------------------------------------------
  subroutine compute_lognorm_weights(dlo, dhi, w, d0, sigln)
    implicit none
    real(kind=8), intent(in)  :: dlo(:), dhi(:)
    real(kind=8), intent(out) :: w(size(dlo))
    real(kind=8), intent(in)  :: d0, sigln
    integer :: i, nb
    real(kind=8) :: dmid, x, sumw

    nb   = size(dlo)
    sumw = 0d0
    do i = 1, nb
      dmid = sqrt(dlo(i) * dhi(i))    ! geometric mean diameter
      x    = log(dmid / d0)
      w(i) = exp( -0.5d0 * (x/sigln)**2 )
      sumw = sumw + w(i)
    end do

    if (sumw > 0d0) then
      do i = 1, nb
        w(i) = w(i) / sumw
      end do
    else
      w = 0d0
      w(1) = 1d0
    end if
  end subroutine compute_lognorm_weights

  !---------------------------------------------------------------------
  ! FUNCTION: fp_nh4no3
  ! PURPOSE
  !   A smooth (logistic) partition factor for NH4NO3 aerosol stability.
  !   fp ~ 1  => nitrate stays in aerosol
  !   fp ~ 0  => nitrate evaporates to NH3 + HNO3
  !
  ! NOTE
  !   This is a crude surrogate for full thermodynamics (e.g.,
  !   ISORROPIA).
  !   Parameters (T0, RH0, betaT, betaRH) are tunable.
  !---------------------------------------------------------------------
  real(kind=8) function fp_nh4no3(T, RH, T0, RH0, betaT, betaRH)
    implicit none
    real(kind=8), intent(in) :: T, RH, T0, RH0, betaT, betaRH
    real(kind=8) :: x
    x = betaT*(T - T0) - betaRH*(RH - RH0)
    fp_nh4no3 = 1d0 / (1d0 + exp(x))
    fp_nh4no3 = clamp01(fp_nh4no3)
  end function fp_nh4no3

end module module_sia_bins_simple 

