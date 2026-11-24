!=============================================================================
! iamas_atmphys_tompkins_cloud.F
!
! Prognostic subgrid cloud fraction scheme following:
!   Tompkins, A. M. (2002). A Prognostic Parameterization for the
!   Subgrid-Scale Variability of Water Vapor and Clouds in Large-Scale
!   Models and Its Use to Diagnose Cloud Cover.
!   J. Atmos. Sci., 59(12), 1917-1942.
!
! Interfaces with:
!   Cumulus scheme
!     -> rqccuten, rqicuten : convective QC/QI tendencies [kg/kg/s]
!        used as proxy for d(M_cu * rc_cu)/dz in Eq.(30)
!   PBL scheme
!     -> tke_pbl: TKE [m2/s2]
!     -> kzq    : PBL moisture diffusivity [m2/s]
!     -> kpbl   : diagnosed PBL top level index
!
! Call sequence (from iamas_atmphys_driver.F):
!   1. PBL scheme      -> tke_pbl, kzq updated
!   2. Cumulus scheme  -> rqccuten, rqicuten updated
!   3. Microphysics    -> snapshot rc_before = qc+qi, call micro, compute delta_rc
!   4. tompkins_cloud_driver() -> update q_tom, bdw_tom, cldfrac
!   5. Radiation       -> receives cldfrac
!
! Registry variables required (add to Registry.xml):
!   q_tom    (nVertLevels, nCells, Time)  real  [1]        shape param q
!   bdw_tom  (nVertLevels, nCells, Time)  real  [kg/kg]    PDF width b-a
!   cldfrac(nVertLevels, nCells, Time) real  [1]        cloud fraction
!
! Author: Jiawang Feng from LACAR
! Date: 2026/05/12
!=============================================================================

module module_tompkins_cloud

   use iamas_kind_types, only : RKIND
   use iamas_log,        only : iamas_log_write

   implicit none
   private

   public :: tompkins_cloud_init
   public :: tompkins_cloud_driver

   !---------------------------------------------------------------------------
   ! Physical and tunable constants
   !---------------------------------------------------------------------------

   ! Fixed beta distribution shape parameter p (Tompkins 2002, Sec. 5e)
   real(RKIND), parameter :: p_fix    = 2.0_RKIND

   ! Symmetric equilibrium value for q (= p_fix, gives p=q=2 bell shape)
   real(RKIND), parameter :: q_min    = 2.0_RKIND

   ! Convective skewness rate constant K (dimensionless, tuned in paper)
   real(RKIND), parameter :: K_conv   = 10.0_RKIND

   ! 3D turbulent dissipation constant kappa (Garratt 1992, Table 8.2)
   real(RKIND), parameter :: kappa_tv = 7.44_RKIND

   ! Smagorinsky constant Cs for horizontal subgrid mixing (Eq. 24)
   real(RKIND), parameter :: Cs_smag  = 0.23_RKIND

   ! von Karman constant
   real(RKIND), parameter :: vonk     = 0.40_RKIND

   ! Blackadar (1962) asymptotic mixing length in PBL [m] (Eq. 19)
   real(RKIND), parameter :: lam_inf  = 150.0_RKIND

   ! Maximum horizontal mixing timescale = 20 days [s] (Appendix)
   real(RKIND), parameter :: tau_h_max = 1.728e6_RKIND

   ! Minimum distribution width = bdw_frac * qv (Appendix safeguard)
   real(RKIND), parameter :: bdw_frac  = 0.10_RKIND

   ! Rd / Rv ratio for saturation mixing ratio calculation
   real(RKIND), parameter :: eps_mv    = 0.6220_RKIND

   ! Small positive number to avoid division by zero
   real(RKIND), parameter :: eps_tiny  = 1.0e-16_RKIND

contains

!=============================================================================
! tompkins_cloud_init
!
! Initialize q_tom and bdw_tom at model start.
! Called once from iamas_atmphys_driver_init() after scalars are allocated.
!
! Input:
!   nCells, nVertLevels  -- grid dimensions
!   qv(nVertLevels,nCells) -- water vapor mixing ratio [kg/kg]
!
! Output:
!   q_tom(nVertLevels,nCells)  -- initialized to symmetric value q_min=2
!   bdw_tom(nVertLevels,nCells)-- initialized to 10% of local qv
!=============================================================================
   subroutine tompkins_cloud_init(nCells, nVertLevels, qv, q_tom, bdw_tom)

      integer,     intent(in)  :: nCells, nVertLevels
      real(RKIND), intent(in)  :: qv(nVertLevels, nCells)
      real(RKIND), intent(out) :: q_tom(nVertLevels, nCells)
      real(RKIND), intent(out) :: bdw_tom(nVertLevels, nCells)

      integer :: k, iCell

      do iCell = 1, nCells
         do k = 1, nVertLevels
            q_tom(k, iCell)   = q_min
            bdw_tom(k, iCell) = bdw_frac * max(qv(k, iCell), 1.0e-8_RKIND)
         end do
      end do

      call iamas_log_write('tompkins_cloud_init: q_tom and bdw_tom initialized.')

   end subroutine tompkins_cloud_init


!=============================================================================
! tompkins_cloud_driver
!
! Main driver, called every time step after PBL, convection, and microphysics.
! Updates prognostic q_tom, bdw_tom and diagnoses cldfrac.
!
! Input:
!   dt                           -- model time step [s]
!   nCells, nVertLevels          -- grid dimensions
!   qv(nVertLevels,nCells)       -- water vapor mix ratio      [kg/kg]
!   qc(nVertLevels,nCells)       -- cloud liquid mix ratio     [kg/kg]
!   qi(nVertLevels,nCells)       -- cloud ice mix ratio        [kg/kg]
!   t_phys(nVertLevels,nCells)   -- temperature                [K]
!   p_phys(nVertLevels,nCells)   -- pressure                   [Pa]
!   rho_phys(nVertLevels,nCells) -- air density                [kg/m3]
!   dz_phys(nVertLevels,nCells)  -- layer geometric thickness  [m]
!   z_phys(nVertLevels,nCells)   -- layer mid-point height AGL [m]
!   k_pbltop(nCells)             -- PBL top level index (from kpbl)
!
!   --- From PBL scheme (after PBL call) ---
!   tke_pbl(nVertLevels,nCells) -- TKE from PBL scheme         [m2/s2]
!   kzq(nVertLevels,nCells)     -- PBL moisture diffusivity    [m2/s]
!
!   --- From Cumulus scheme (after convection call) ---
!   rqccuten(nVertLevels,nCells) -- convective QC tendency     [kg/kg/s]
!   rqicuten(nVertLevels,nCells) -- convective QI tendency     [kg/kg/s]
!
!   --- From microphysics ---
!   delta_rc(nVertLevels,nCells) -- change in (qc+qi) due to microphysics [kg/kg]
!                                   = (qc+qi)_after - (qc+qi)_before
!                                   pass array of zeros if not using micro feedback
!
! Input/Output:
!   q_tom(nVertLevels,nCells)    -- beta PDF shape parameter q (prognostic)
!   bdw_tom(nVertLevels,nCells)  -- beta PDF width b-a         (quasi-prognostic)
!
! Output:
!   cldfrac(nVertLevels,nCells) -- diagnosed cloud fraction [0,1]
!=============================================================================
   subroutine tompkins_cloud_driver(                                       &
         dt, nCells, nVertLevels,                                          &
         qv, qc, qi,                                                       &
         T_phys, p_phys, rho_phys, dz_phys, z_phys,                        &
         k_pbltop,                                                         &
         tke_pbl, kzq,                                                     &
         rqccuten, rqicuten,                                               &
         delta_rc,                                                         &
         kdiff, len_disp,                                                  &
         q_tom, bdw_tom,                                                   &
         cldfrac)

      integer,     intent(in)  :: nCells, nVertLevels
      real(RKIND), intent(in)  :: dt, len_disp

      real(RKIND), intent(in), dimension(nVertLevels, nCells) ::           &
         qv, qc, qi,                                                       &
         T_phys, p_phys, rho_phys, dz_phys, z_phys,                        &
         tke_pbl, kzq,                                                     &
         rqccuten, rqicuten,                                               &
         delta_rc,                                                         &
         kdiff

      integer,     intent(in),    dimension(nCells)           :: k_pbltop

      real(RKIND), intent(inout), dimension(nVertLevels, nCells) ::        &
         q_tom, bdw_tom

      real(RKIND), intent(out),   dimension(nVertLevels, nCells) ::        &
         cldfrac

      !--- Local variables ---------------------------------------------------
      integer     :: iCell, k
      real(RKIND) :: rt, rs, rv_bar, rc_bar
      real(RKIND) :: a_val, b_val, cfrac
      real(RKIND) :: tke_local, mix_len
      real(RKIND) :: tau_v, tau_h
      real(RKIND) :: dqrt_dz
      real(RKIND) :: prod_bdw
      real(RKIND) :: eta
      real(RKIND) :: dq_conv, dq_turb, dq_micro
      real(RKIND) :: dbw_conv, dbw_diss, dbw_turb
      real(RKIND) :: Drc_Dz
      real(RKIND) :: bdw_min_local
      real(RKIND) :: q_new

      do iCell = 1, nCells
         do k = 1, nVertLevels

            !-----------------------------------------------------------------
            ! [1] Basic thermodynamic quantities
            !-----------------------------------------------------------------
            rt     = qv(k,iCell) + qc(k,iCell) + qi(k,iCell)
            rs     = sat_mix_ratio(T_phys(k,iCell), p_phys(k,iCell))
            rv_bar = qv(k,iCell)
            rc_bar = qc(k,iCell) + qi(k,iCell)

            ! Physical lower bound for distribution width (Appendix)
            bdw_min_local = bdw_frac * max(qv(k,iCell), 1.0e-8_RKIND)

            !-----------------------------------------------------------------
            ! [2] Diagnose cloud fraction and PDF bounds a, b
            !     using current prognostic q_tom and bdw_tom
            !-----------------------------------------------------------------
            call diagnose_pdf(rt, rs, rv_bar, rc_bar,                      &
                              q_tom(k,iCell), bdw_tom(k,iCell),            &
                              a_val, b_val, cfrac)
            cldfrac(k,iCell) = cfrac

            !-----------------------------------------------------------------
            ! [3] Compute turbulence timescales
            !
            ! tau_v (Eq. 22): 3D turbulent eddy timescale
            !   Applied only within the PBL (k <= k_pbltop).
            !   Above PBL: 3D turbulence negligible -> tau_v = very large.
            !
            ! tau_h (Eq. 25): horizontal Smagorinsky mixing timescale
            !   Applied at all levels.
            !-----------------------------------------------------------------
            tke_local = max(tke_pbl(k,iCell), 1.0e-6_RKIND)

            if (k <= k_pbltop(iCell)) then
               ! Blackadar (1962) mixing length (Eq. 19)
               mix_len = vonk * z_phys(k,iCell) /                          &
                         (1.0_RKIND + vonk * z_phys(k,iCell) / lam_inf)
               ! 3D turbulent dissipation timescale (Eq. 22)
               tau_v   = kappa_tv * mix_len / sqrt(tke_local)
               tau_v   = max(tau_v, dt)   ! floor at one time step
            else
               tau_v = 1.0e30_RKIND   ! effectively no 3D dissipation above PBL
            end if

            ! Horizontal Smagorinsky mixing timescale (Eq. 25)
            !   tau_h^{-1} = Cs^2 * sqrt( (du/dx)^2 + (dv/dy)^2 )
            !   Use tau_h^{-1} = kdiff / config_len_disp^2 instead.
            ! tau_h = 1.0_RKIND / max(                                        &
            !            Cs_smag**2 * sqrt(dudx(k,iCell)**2 + dvdy(k,iCell)**2), &
            !            1.0_RKIND / tau_h_max)
            tau_h = 1.0_RKIND / max(                 &
                       kdiff(k,iCell) / (len_disp**2), &
                       1.0_RKIND / tau_h_max)
            tau_h = min(tau_h, tau_h_max)

            !-----------------------------------------------------------------
            ! [4] Compute q_tom tendencies from three processes
            !
            !  dq/dt = dq_conv + dq_turb + dq_micro
            !-----------------------------------------------------------------

            ! [4a] Convective skewness source (Eq. 30)
            !      dq/dt = (K / rs) * d(M_cu * rc_cu) / dz
            !
            !      Grell-Freitas does not directly output M_cu.
            !      Approximation:  d(M_cu*rc_cu)/dz ~ rho * (rqccuten + rqicuten)
            !      Only positive values (detrainment = skewness source).
            Drc_Dz  = rho_phys(k,iCell) *                                  &
                      max(0.0_RKIND, rqccuten(k,iCell) + rqicuten(k,iCell))
            dq_conv = K_conv / max(rs, eps_tiny) * Drc_Dz * dt

            ! [4b] Turbulent + horizontal mixing relaxation toward q0 (Eq. 26)
            !      dq/dt = (q0 - q) * (1/tau_v + 1/tau_h)
            dq_turb = (q_min - q_tom(k,iCell)) *                           &
                      (1.0_RKIND/tau_v + 1.0_RKIND/tau_h) * dt

            ! [4c] Microphysics feedback (Eq. 34, linearized)
            !      Precipitation removes condensate -> q relaxes toward q0.
            !      Only applied when condensate decreases (delta_rc < 0).
            dq_micro = 0.0_RKIND
            if (rc_bar > eps_tiny .and. delta_rc(k,iCell) < 0.0_RKIND) then
               dq_micro = (b_val + (delta_rc(k,iCell) / rc_bar) *          &
                          (b_val - rs)) * p_fix /                           &
                          max(rt - a_val, eps_tiny)                        &
                          - (p_fix + q_tom(k,iCell))
               ! Limit to avoid overshoot within a single time step
               dq_micro = dq_micro * min(dt / max(tau_v, 1.0_RKIND), 1.0_RKIND)
            end if

            ! Update q_tom; enforce floor at q_min = 2 (symmetric distribution)
            q_new = q_tom(k,iCell) + dq_conv + dq_turb + dq_micro
            q_tom(k,iCell) = max(q_min, q_new)

            !-----------------------------------------------------------------
            ! [5] Update bdw_tom (distribution width b-a)
            !
            ! Partially cloudy:  bdw is DIAGNOSED from iteration (Sec. 4).
            !   The prognostic value is slaved to the diagnosed value.
            !
            ! Clear sky / overcast:  prognostic equation is active (Eq. 36).
            !   Sources: convective widening (Eq. 31), turbulent production (Eq. 20).
            !   Sinks:   turbulent + horizontal mixing dissipation.
            !-----------------------------------------------------------------
            if (cfrac > 0.0_RKIND .and. cfrac < 1.0_RKIND) then

               !--- Partial cloud: slave bdw to iteration-diagnosed value ---
               bdw_tom(k,iCell) = b_val - a_val

            else

               !--- Clear sky or overcast: integrate prognostic equation ---

               ! eta = (p+q)^2*(p+q+1) / (p*q)  [shape factor for Eq. 20]
               eta = (p_fix + q_tom(k,iCell))**2 *                         &
                     (p_fix + q_tom(k,iCell) + 1.0_RKIND) /                &
                     (p_fix * q_tom(k,iCell))

               ! Turbulent production of bdw (Eq. 20 term 1):
               !   -(eta/bdw) * w'qt' * d(qt)/dz
               !   w'qt' approximated as -Kq * d(qt)/dz  (PBL diffusivity)
               !   d(qt)/dz via upwind finite difference
               if (k < nVertLevels) then
                  dqrt_dz = ((qv(k+1,iCell) + qc(k+1,iCell) + qi(k+1,iCell)) - &
                             (qv(k,  iCell) + qc(k,  iCell) + qi(k,  iCell))) / &
                             max(dz_phys(k,iCell), 1.0_RKIND)
               else
                  dqrt_dz = 0.0_RKIND
               end if
               ! w'qt' = -Kq * d(qt)/dz; production = -(eta/bdw) * w'qt' * d(qt)/dz
               !       = (eta/bdw) * Kq * (d(qt)/dz)^2    [always non-negative]
               prod_bdw = (eta / max(bdw_tom(k,iCell), eps_tiny)) *        &
                           kzq(k,iCell) * dqrt_dz**2

               ! Turbulent transport term (Eq. 20 term 2, -d(w'bdw')/dz):
               !   Handled implicitly by iAMAS scalar transport.
               !   Set to 0 here to avoid double-counting.
               dbw_turb = 0.0_RKIND

               ! Convective widening from updated q (Eq. 31):
               !   bdw -> (rt-a)*(p+q)/p  as q increases
               if (rt > a_val + eps_tiny) then
                  dbw_conv = ((rt - a_val) * (p_fix + q_tom(k,iCell)) /    &
                              p_fix - bdw_tom(k,iCell)) *                  &
                              min(dt / max(tau_v, 1.0_RKIND), 1.0_RKIND)
               else
                  dbw_conv = 0.0_RKIND
               end if

               ! Turbulent + horizontal dissipation of bdw (Eq. 20 last term):
               !   -bdw * (1/tau_v + 1/tau_h)
               dbw_diss = -bdw_tom(k,iCell) *                              &
                           (1.0_RKIND/tau_v + 1.0_RKIND/tau_h) * dt

               ! Forward Euler step for bdw
               bdw_tom(k,iCell) = bdw_tom(k,iCell)                        &
                                 + (prod_bdw + dbw_turb + dbw_conv + dbw_diss) * dt

            end if  ! partial / clear / overcast

            ! Enforce physical lower bound (Appendix)
            bdw_tom(k,iCell) = max(bdw_tom(k,iCell), bdw_min_local)

         end do  ! k
      end do     ! iCell

   end subroutine tompkins_cloud_driver


!=============================================================================
! diagnose_pdf  (private)
!
! Given current prognostic state (q, bdw) and water amounts, determine:
!   (1) PDF bounds a, b
!   (2) Cloud fraction cfrac via analytic beta integral
!
! Three regimes handled:
!   rc_bar = 0              -> clear sky,  cfrac = 0,  bounds from bdw
!   rv_bar >= rs            -> overcast,   cfrac = 1,  bounds from bdw
!   0 < rc_bar < rv_bar     -> partial cloud, bisect Eq.(16) for a
!=============================================================================
   subroutine diagnose_pdf(rt, rs, rv_bar, rc_bar, q, bdw, a_out, b_out, cfrac)

      real(RKIND), intent(in)  :: rt, rs, rv_bar, rc_bar, q, bdw
      real(RKIND), intent(out) :: a_out, b_out, cfrac

      real(RKIND) :: xs

      if (rc_bar <= eps_tiny) then
         !--- Clear sky ---
         cfrac = 0.0_RKIND
         ! Eq.(12): rt = (b-a)*p/(p+q) + a  =>  a = rt - bdw*p/(p+q)
         a_out = rt - bdw * p_fix / (p_fix + q)
         a_out = max(0.0_RKIND, a_out)
         b_out = a_out + bdw

      else if (rv_bar >= rs - eps_tiny) then
         !--- Fully overcast ---
         cfrac = 1.0_RKIND
         a_out = rt - bdw * p_fix / (p_fix + q)
         b_out = a_out + bdw

      else
         !--- Partial cloud: solve Eq.(16) for lower bound a ---
         call bisect_for_a(rt, rs, rv_bar, q, a_out, b_out)

         ! Cloud fraction from analytic I_x(p,q) with p=2 (Eq. 2 + Eq. 15):
         !   C = 1 - I_{xs}(2, q)
         !     = (1 - xs)^q * (1 + q*xs)
         ! where xs = p*(rs-a) / [(rt-a)*(p+q)]  is the normalised saturation point
         xs = p_fix * (rs - a_out) / ((rt - a_out) * (p_fix + q))
         xs = max(0.0_RKIND, min(1.0_RKIND, xs))

         cfrac = (1.0_RKIND - xs)**q * (1.0_RKIND + q * xs)
         cfrac = max(0.0_RKIND, min(1.0_RKIND, cfrac))

      end if

   end subroutine diagnose_pdf


!=============================================================================
! bisect_for_a  (private)
!
! Bisection solver for the PDF lower bound a satisfying Eq.(16):
!
!   F(a) = (rt-a)*I_{x0}(3,q) + (a-rs)*I_{x0}(2,q) + rs - rv_bar = 0
!
! where  x0 = 2*(rs-a) / [(rt-a)*(2+q)]   (normalised sat. point with p=2)
!
! F is monotonically decreasing in a -> bisection is globally convergent.
! Search interval:  a in [0, min(rt, rs))
! Tolerance: 1e-13 kg/kg (sub-nanogram level, well within float64 precision)
!=============================================================================
   subroutine bisect_for_a(rt, rs, rv_bar, q, a_out, b_out)

      real(RKIND), intent(in)  :: rt, rs, rv_bar, q
      real(RKIND), intent(out) :: a_out, b_out

      integer,     parameter :: max_iter = 80
      real(RKIND), parameter :: tol_a    = 1.0e-13_RKIND

      real(RKIND) :: a_lo, a_hi, a_mid
      real(RKIND) :: f_lo, f_hi, f_mid
      integer     :: iter

      ! Search bounds: a in [0, min(rt,rs))
      a_lo =  0.0_RKIND
      a_hi = min(rt, rs) * (1.0_RKIND - 1.0e-9_RKIND)

      f_lo = residual_eq16(a_lo, rt, rs, rv_bar, q)
      f_hi = residual_eq16(a_hi, rt, rs, rv_bar, q)

      ! Fallback: if no sign change, snap to the nearer bound
      if (f_lo * f_hi > 0.0_RKIND) then
         if (abs(f_lo) <= abs(f_hi)) then
            a_out = a_lo
         else
            a_out = a_hi
         end if
         b_out = a_out + (rt - a_out) * (p_fix + q) / p_fix
         return
      end if

      do iter = 1, max_iter
         a_mid = 0.5_RKIND * (a_lo + a_hi)
         if ((a_hi - a_lo) < tol_a) exit
         f_mid = residual_eq16(a_mid, rt, rs, rv_bar, q)
         if (f_lo * f_mid <= 0.0_RKIND) then
            a_hi = a_mid
            f_hi = f_mid
         else
            a_lo = a_mid
            f_lo = f_mid
         end if
      end do

      a_out = 0.5_RKIND * (a_lo + a_hi)

      ! Recover b from Eq.(12): b = a + (rt-a)*(p+q)/p
      b_out = a_out + (rt - a_out) * (p_fix + q) / p_fix

   end subroutine bisect_for_a


!=============================================================================
! residual_eq16  (private, pure)
!
! Evaluates the residual of Tompkins (2002) Eq.(16):
!
!   F(a) = (rt-a) * I_{x0}(p+1, q) + (a-rs) * I_{x0}(p, q) + rs - rv_bar
!
! With p=2 fixed, the incomplete beta ratios are evaluated analytically:
!
!   I_x(2, q) = 1 - (1-x)^q * (1 + q*x)
!   I_x(3, q) = 1 - (1-x)^q * [1 + q*x + q*(q+1)/2 * x^2]
!
! Derivation of I_x(3,q):  use the recurrence relation
!   I_x(p+1, q) = I_x(p, q) - x^p*(1-x)^q / [p*B(p,q)]
! with B(2,q) = 1 / [q*(q+1)] to get the closed form above.
!=============================================================================
   pure function residual_eq16(a, rt, rs, rv_bar, q) result(f)

      real(RKIND), intent(in) :: a, rt, rs, rv_bar, q
      real(RKIND) :: f, x0, Ix_2q, Ix_3q, omx

      ! Normalised saturation point: x0 = p*(rs-a) / [(rt-a)*(p+q)]
      x0 = p_fix * (rs - a) / ((rt - a) * (p_fix + q))
      x0 = max(0.0_RKIND, min(1.0_RKIND, x0))

      omx = 1.0_RKIND - x0   ! (1 - x0), computed once

      ! Analytic I_x(2, q) with p=2
      Ix_2q = 1.0_RKIND - omx**q * (1.0_RKIND + q * x0)

      ! Analytic I_x(3, q) with p=2
      Ix_3q = 1.0_RKIND - omx**q *                                         &
              (1.0_RKIND + q * x0 + 0.5_RKIND * q * (q + 1.0_RKIND) * x0**2)

      ! Eq.(16) residual
      f = (rt - a) * Ix_3q + (a - rs) * Ix_2q + rs - rv_bar

   end function residual_eq16


!=============================================================================
! sat_mix_ratio  (private, pure)
!
! Saturation water vapor mixing ratio via Tetens approximation.
!   T >= -23 C : over liquid water (Bolton 1980)
!   T <  -23 C : over ice         (Murray 1967)
!
! Arguments:
!   T_K  [K]  -- temperature
!   p_Pa [Pa] -- pressure
!=============================================================================
   pure function sat_mix_ratio(T_K, p_Pa) result(rs)

      real(RKIND), intent(in) :: T_K, p_Pa
      real(RKIND) :: rs, es, Tc

      Tc = T_K - 273.15_RKIND

      if (Tc >= -23.0_RKIND) then
         ! Bolton (1980), valid for liquid
         es = 611.2_RKIND * exp(17.67_RKIND * Tc / (Tc + 243.5_RKIND))
      else
         ! Murray (1967), valid for ice
         es = 611.2_RKIND * exp(21.874_RKIND * Tc / (Tc + 265.5_RKIND))
      end if

      ! rs = eps * es / (p - es)
      rs = eps_mv * es / max(p_Pa - es, 1.0_RKIND)

   end function sat_mix_ratio

end module module_tompkins_cloud
