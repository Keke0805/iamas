!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute linoz subroutine for stratospheric ozone chemistry.
! * Author: Zihan Xia
! * Email: xiazh@mail.ustc.edu.cn
! Last update: zhxia 10/2025
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module module_linoz 
  use iamas_kind_types
  use iamas_derived_types
  use iamas_pool_routines
  use iamas_atmchem_vars
  use iamas_atmchem_configs
  use iamas_atmchem_utilities
  use iamas_atmphys_vars
  use iamas_log, only : iamas_log_write
 
  implicit none
  private
  public  :: linoz_driver
  public  :: linoz_read

  private :: linoz_chem3
  private :: linoz_stratl
  private :: linoz_strt2m
  private :: linoz_somlfq

  real(kind=RKIND), dimension(:,:,:,:), allocatable, private :: linoz_tparm

! !REMARKS:
!  LINOZ Climatology:
!  ============================================================================
!  The LINOZ stratospheric chemistry tables (eg., Linoz_March2007.dat) for ozone consist of:
!
!  7 tables, each a function of:
!    12 months,
!    18 latitudes (-85 to 85 in 10 deg. increments)
!    25 altitudes ( z*=10-58 km in 2 km increments)
!
!  The 7 data fields are:
!    1- ozone (Logan climatology), v/v
!    2- Temperature climatology, K
!    3- Column ozone climatology, Logan ozone integrated above box, DU
!    4- ozone (P-L) for climatological ozone, v/v/s
!    5- d(P-L) / dO3, 1/s
!    6- d(P-L) / dT, v/v/s/K
!    7- d(P-L) / d(column O3), v/v/s/DU
!------------------------------------------------------------------------------
CONTAINS

  subroutine linoz_driver( dtstepc, diag, state, configs, mesh  , &
                           its, ite, kts, kte                     &
                         )
  
    IMPLICIT NONE
    ! input/ouput vars
    real(kind=RKIND), intent(in) :: dtstepc 
    type(iamas_pool_type),intent(in):: diag
    type(iamas_pool_type),intent(inout):: state
    type(iamas_pool_type),intent(in):: mesh
    type(iamas_pool_type),intent(in):: configs
    integer, intent(in) :: its,ite,kts,kte

    !local variables:
    integer, save :: lastmonth = -99
    integer :: month
    real,pointer :: config_chemdt

    ! current month
    month = curr_month !zhxia compute this in chem_driver

    ! if new month, get new parameters
    if ( month /= lastmonth ) then
       call linoz_stratl(state = state, mesh = mesh      , &
                         its=its, ite=ite, kts=kts, kte=kte)                    

       lastmonth =  month
    endif

    ! linoz needs time step in seconds
    call iamas_pool_get_config(configs, 'config_chemdt', config_chemdt)

    ! call the linoz chemistry
    call linoz_chem3( dtstepc = dtstepc, diag = diag, state = state, mesh = mesh  , &
                      its=its, ite=ite, kts=kts, kte=kte)

  end subroutine linoz_driver

!---------------------------------------------------------------------

  subroutine linoz_chem3( dtstepc, diag, state, mesh, its, ite, kts, kte)
  !description: applies linearized chemistry based on tables from pratmo model using climatological t, o3, time of year

    IMPLICIT NONE
    ! input/ouput vars
    real(kind=RKIND),intent(in)    :: dtstepc      ! chemical time step in current time[seconds]
    type(iamas_pool_type),intent(inout):: state
    type(iamas_pool_type),intent(in):: mesh
    type(iamas_pool_type),intent(in):: diag
    integer, intent(in) :: its,ite,kts,kte !must kts=1,kte=nVertLevels. zhxia.

    !local variables:
    ! scalars
    integer,pointer :: nCells, nVertLevels
    integer, dimension(:), pointer :: is_trop_layer
    real(kind=RKIND),dimension(:,:,:),pointer :: linoz_tlstt   ! K,i,linoz_nfields
    real(kind=RKIND), dimension(:), pointer :: areaCell
    integer  :: iCell, l, lbot, lpos 
    real(kind=RKIND) :: climo3,  climpml,    dco3, dero3,     dertmp
    real(kind=RKIND) :: derco3,  dmass,      dtmp, sso3
    real(kind=RKIND) :: air_mass_in_cell, o3_mass_in_cell, o3_molecules, area_cm2

    ! we need to define local arrays to hold corresponding values
    real(kind=RKIND), allocatable :: o3_conc(:,:),dcolo3(:,:),colo3(:,:),out_data(:,:)

    ! constant
    real(kind=RKIND), parameter :: AIR_MW = 28.97_RKIND !g/mol
    real(kind=RKIND), parameter :: O3_MW = 48.0_RKIND !g/mol
    real(kind=RKIND), parameter :: DU_CONST = 2.687e16_RKIND !molecules/cm²
    real(kind=RKIND), parameter :: avo = 6.022e23_RKIND !molecules/mol
 
    !debug zhxia
    real(kind=RKIND), dimension(:),pointer :: latCell, lonCell
    real(kind=RKIND) :: lat, lon
    real, parameter :: pi = 3.14159265358979323846

    !=================================================================
    ! linoz_chem3 begins here!
    !=================================================================

    ! dimensions
    call iamas_pool_get_dimension(mesh, 'nCells', nCells)
    call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

    ! point to  array
    call iamas_pool_get_array(diag, 'is_trop_layer', is_trop_layer)
    call iamas_pool_get_array(mesh, 'areaCell', areaCell)
    call iamas_pool_get_array(state, 'linoz_tlstt', linoz_tlstt)
  
    ! echo info
    call chemistry_debug_message('linoz_chem3: doing linoz')

    !=================================================================
    ! perform stratospheric chemistry
    !=================================================================
    allocate(o3_conc(kts:kte,its:ite))
    allocate(dcolo3(kts:kte,its:ite))
    allocate(colo3(kts:kte,its:ite))
    allocate(out_data(kts:kte,its:ite))
    o3_conc(kts:kte,its:ite) = chem(p_o3,kts:kte,its:ite) / 1000000.0 !ppmv->v/v

    ! start at top layer and continue to lowest layer for strat. chem
    out_data = 0e+0

    do iCell = its,ite

       lbot = is_trop_layer(iCell) + 1
       lpos = 1
       if (pres2_hyd(kte,iCell)*0.01 .ge. 0.3e+0) then !pa->hPa
          lpos = kte
       else
          do while ((lpos < nVertLevels) .and. (pres2_hyd(lpos+1,iCell)*0.01 .ge. 0.3e+0)) !pa->hPa
             lpos = lpos +1
          enddo
          lpos = lpos-1
       endif

       ! loop over levels
       do l = kte, lbot, -1

          ! skip if tracer is negative
          if ( o3_conc(l,iCell) .lt. 0.e+0 ) cycle
  
          !air_mass_in_cell, unit:kg
          air_mass_in_cell = rho_phy(l,iCell)*areaCell(iCell)*dz(l,iCell)
          !o3_mass_in_cell, unit:kg
          o3_mass_in_cell = o3_conc(l,iCell)*air_mass_in_cell*(O3_MW / AIR_MW)
          !o3_molecules, unit:molecules
          o3_molecules = (o3_mass_in_cell * 1000._RKIND / O3_MW) * avo
          !area_cm2, unit:cm2
          area_cm2 = areaCell(iCell) * 1.0e4_RKIND
          
          ! calculate ozone column above box (and save)
          ! dcolo3 = ozone column (in du) in given layer
          ! colo3 =  ozone column above layer + half of column in layer

          ! bdf spc is in v/v, make conversion to du
          dcolo3(l,iCell) = (o3_molecules / area_cm2) / DU_CONST
          if (l.eq.nVertLevels) then !top model layer
             colo3(l,iCell) = dcolo3(l,iCell)*0.5
          else
             colo3(l,iCell) = colo3(l+1,iCell) + (dcolo3(l,iCell)+dcolo3(l+1,iCell))*0.5
          endif
          out_data(l,iCell) = colo3(l,iCell)

          ! ++++++ climatological p-l:   ++++++
          climpml=linoz_tlstt(l,iCell,4)      ! climatological p-l = (p-l)^o

          ! ++++++ local ozone feedback: ++++++
          dero3=linoz_tlstt(l,iCell,5)               ! derivative w.r.t. o3. 
                                             !  dero3=-1/(time constant)
          if (dero3.eq.0) cycle              ! skip linoz if lifetime
                                             !  is infinite.
          climo3=linoz_tlstt(l,iCell,1)              ! climatological o3 = f^o
          derco3=linoz_tlstt(l,iCell,7)              ! derivative w.r.t. column o3
          dco3=(colo3(l,iCell)-linoz_tlstt(l,iCell,3)) ! deviation from o3 climatology.
          ! ++++++ temperature feedback: ++++++
          dertmp=linoz_tlstt(l,iCell,6)              ! derivative w.r.t. temperature
          dtmp=(t_phy(l,iCell) - &       ! deviation in temperature
                   linoz_tlstt(l,iCell,2))           !  from climatology.

          ! ++++++ calculate steady-state ozone: ++++++
          sso3=climo3 - (climpml+dtmp*dertmp+dco3*derco3)/dero3

          ! ++++++ change in ozone mass due to chemistry: ++++++
          !sso3 = f^*
          dmass=(sso3-o3_conc(l,iCell))*(1.0-exp(dero3*dtstepc))

          ! ++++++ update ozone mass ++++++
          ! linox valid only up to 58 km, so do not use above 0.3 hpa
          ! dbj: impose exponential fall off of mixing ratio
          ! between 0.3 and 0.01 hpa (with fall off of a scale height)
          if (pres2_hyd(l,iCell)*0.01 .le. 0.3e+0) then !pa->hPa
             o3_conc(l,iCell) = &
                  (pres_hyd(l,iCell)/pres_hyd(lpos-1,iCell))  &
                  * o3_conc(lpos-1,iCell)
          else
             o3_conc(l,iCell) = o3_conc(l,iCell)+dmass
          endif

          chem(p_o3,l,iCell) =  o3_conc(l,iCell) * 1000000.0 !v/v->ppmv

       enddo       ! loop over l

    enddo          ! loop over iCell

    !write our calculated column o3 maximum
    !write(6,*) 'max of columns= ',maxval(out_data)

    deallocate(o3_conc)
    deallocate(dcolo3)
    deallocate(colo3)
    deallocate(out_data)

    call chemistry_debug_message('done linoz_chem3')

  end subroutine linoz_chem3

!---------------------------------------------------------------------

  subroutine linoz_stratl( state, mesh, its, ite, kts, kte)
  !description: performs a monthly fixup of chemistry parameters for the linoz stratospheric ozone chemistry.

    IMPLICIT NONE

!  the linoz climatology array is linoz_tparm(25,18,12,n),
!  which has the following dimensions
!    * 25 layers from 58 km to 10 km by 2 km intervals
!    * 18 latitudes (85s, 75s, ...85n)
!    * 12 months
!    *  n fields (currently n=7)
!

    ! input/ouput vars
    type(iamas_pool_type),intent(inout):: state
    type(iamas_pool_type),intent(in):: mesh
    integer, intent(in) :: its,ite,kts,kte !must kts=1,kte=nVertLevels. zhxia.

    !local variables:
    real(kind=RKIND),dimension(:,:,:),pointer :: linoz_tlstt
    real(kind=RKIND), dimension(:), pointer :: latCell
    integer, pointer :: nfields_linoz,nlat_linoz,nlevels_linoz,nmonths_linoz
    integer,pointer :: nCells, nVertLevels
    integer :: iCell, ii, ixxx, k, l, lr, n, month
    real(kind=RKIND) :: current_lat_rad, current_lat_deg
    real(kind=RKIND) :: pi = 3.141592653589793

    ! arrays
    real(kind=RKIND),dimension(:),allocatable  :: p0l
    real(kind=RKIND),dimension(:),allocatable  :: strtx
    real(kind=RKIND),dimension(:),allocatable  :: strt0l, strt1l, strt2l

    ! month names
    character(len=3)  :: cmonth(12) = (/'jan', 'feb', 'mar', 'apr', &
                                        'may', 'jun', 'jul', 'aug', &
                                        'sep', 'oct', 'nov', 'dec'/)

    !=================================================================
    ! linoz_stratl begins here!
    !=================================================================

    ! initialization
    call iamas_pool_get_dimension(mesh, 'nCells', nCells)
    call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels) 
    call iamas_pool_get_dimension(state, 'linoz_Nlevels', nlevels_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nlat', nlat_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nmonths',nmonths_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nfields',nfields_linoz)
    call iamas_pool_get_array(state, 'linoz_tlstt', linoz_tlstt)
    call iamas_pool_get_array(mesh,'latCell',latCell)
    IF (.NOT. ASSOCIATED(linoz_tlstt)) THEN
        call chemistry_error_fatal('Error: Failed to get pointer to "linoz_tlstt" in linoz_stratl!')
    ENDIF

    month = curr_month
    allocate(p0l(kte+1))
    allocate(strtx(nlevels_linoz))
    allocate(strt0l(kte+1))
    allocate(strt1l(kte+1))
    allocate(strt2l(kte+1))

    ! echo info to stdout
    call chemistry_debug_message('Interpolating linoz fields for '//cmonth(month))

    !=================================================================
    ! lookup data in the linoz climatology
    !=================================================================

    ! loop over the # of fields in the linoz climatology
    do n = 1, nfields_linoz

       ! loop over iAMAS latitudes
       do iCell = its, ite

          ! ii is the index of the nearest linoz data column corresponding to iAMAS latitude index iCell
          current_lat_rad = latCell(iCell)
          current_lat_deg = current_lat_rad * (180.0 / pi)
          ixxx = int( 0.1e+0 * current_lat_deg + 10.e+0 )
          ii   = min( 18, max( 1, ixxx ) )

          ! loop over the # of levels in the linoz climatology
          do k = 1, nlevels_linoz

             ! strtx(k) is the column of data from the linoz climatology for the given month, latitude (jj), level (k), and species (n).

             ! nearest neighbor, no interpolation
             strtx(k) = linoz_tparm(k,ii,month,n)
          enddo

          ! p0l are the pressure at the level edges of each iAMAS grid box
          do l = kts, kte+1
             p0l(l) = pres2_hyd(kte+2-l,iCell)*0.01 !inverse index and transform unit. Pa -> hPa
          enddo

          ! interpolate and calculate moments of column distribution
          call linoz_strt2m( nlevels_linoz = nlevels_linoz,        &
                             nVertLevels = kte, strtx = strtx, p0l = p0l,    &
                             strt0l = strt0l, strt1l = strt1l, strt2l = strt2l)

          ! store loss freq/yields & moments in tlstt/swt/sww for exact ctm layers lm down order reversed from c.mclinden version
          do lr = kts,kte
             linoz_tlstt(lr, iCell, n) = strt0l(kte+1-lr)
          enddo

       enddo   ! loop over iCell

    enddo      ! loop over n

    deallocate(p0l)
    deallocate(strtx)
    deallocate(strt0l)
    deallocate(strt1l)
    deallocate(strt2l)

  end subroutine linoz_stratl

!--------------------------------------------------------------------

  subroutine linoz_strt2m( nlevels_linoz,         &
                           nVertLevels,     strtx,      p0l,    &
                           strt0l,    strt1l,     strt2l  )
  !description: interpolates quantities from the linoz vertical grid to the iAMAS vertical grid. it also computes the 1st \& 2nd moments of the distribution.

    IMPLICIT NONE
    ! input/ouput vars
    integer, intent(in):: nlevels_linoz
    integer,          intent(in)  :: nVertLevels        ! number of grid levels
    real(kind=RKIND), intent(in)  :: p0l(nVertLevels+1) ! pressure edges !nVertLevels=kte

    ! quantity on the linoz vertical grid (i.e. fields #1-7 of the linoz clim.
    real(kind=RKIND), intent(in)  :: strtx(nlevels_linoz)
!
! !output parameters:
!
    ! 0th moment of distribution, on geos-chem grid edges
    real(kind=RKIND), intent(out) :: strt0l(nVertLevels+1)

    ! 1st moment of distribution, on geos-chem grid edges
    real(kind=RKIND), intent(out) :: strt1l(nVertLevels+1)

    ! 2nd moment of distribution, on geos-chem grid edges
    real(kind=RKIND), intent(out) :: strt2l(nVertLevels+1)

!  call somlfq(p1,p2,f0,f1,f2,ps,f,nl)
!  - p1,p2 are the pressure edges for the ctm layer onto which the
!    coefficients will be mapped. [p1>p2 i believe {pjc}]
!  - f0,f1,f2 are the ctm layer vertical moments determined in somlfq
!  - ps are the pressure layer edges of the original [ie linox] grid
!  - f is the column of coefficients (on the original grid); note
!    f is flipped relative to strtx and since the coefficients begin
!    at z*=10, f(1)=f(2)=...=f(5)=0
!  - nl is 30; size of f()
!                                                                             .
!   the box model calculations were performed at z*=10km, 12km, ... and
!   so these would represent the centres with the corresponding edges at
!   9,11km ; 11,13km; ...
!   ps() represents the edges (although ps(1) is set to 1000mb).
!   the first few values are:
!     ps(1)=1000
!     ps(2)=874.947105    (note ps(2) is not quite 1000 exp(-1/16) as the
!     ps(3)=656.117767     the average pressure is used - not the pressure
!     ps(4)=492.018914     at the average z*)
!     ps(5)=368.96213
!     ps(6)=276.68257
!     ps(7)=207.48266
!     ...
!     ps(30)=0.276682568
!     ps(31)=0.0
!                                                                             .
!     f(1) spans ps(1)-ps(2)
!     f(2) spans ps(2)-ps(3)
!     ...
!     f(30) spans ps(30)-ps(31)

    !local variables:
    ! scalars
    integer  :: nl,  ncbox, l,      k,    nx
    real(kind=RKIND) :: p1,  p2,    f0,     f1
    real(kind=RKIND) :: f2,  xpsd,  xpslm1, xpsl

    ! arrays
    real(kind=RKIND) :: ps(nlevels_linoz+6)
    real(kind=RKIND) :: f (nlevels_linoz+5) 

    !=================================================================
    ! initialization
    !=================================================================

    ! now make nx a local variable, since we already pass
    nx         = nlevels_linoz

    ! now make nl a local variable instead of a parameter
    nl         = nlevels_linoz + 5

    !=================================================================
    ! set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
    !
    ! assume that stratospheric chemical parameters always start at
    ! 52 km (n=27).  scan downward from 52 km to 14 km (nx=20) by
    ! 2 km.
    !
    ! 58 km (n=30) scan downward from 58 km to 10 km (nx=25) by 2 km
    ! intervals, constant >58km
    !
    !  n.b. f(@30km) assumed to be constant from 29-31 km (by mass)
    !=================================================================

    xpsd       = 10.e+0 **(-0.125e+0)
    xpslm1     = 1000.e+0
    ps(1)      = 1000.e+0
    do l = 2,nl
       xpsl    = xpslm1 *xpsd
       ps(l)   = 0.5e+0 *(xpslm1 +xpsl)
       xpslm1  = xpsl
    enddo
    ps(nl+1)   = 0.e+0
    do l = 1,nl-nx
       f(l)     = 0.e+0
    enddo

    ! k=1 is at the top of atmosphere
    do k = 1,nx
       f(nl+1-k)= strtx(k) !strtx has increasing preasure.
    enddo

    do k = 1,nVertLevels
       p1       = p0l(k+1)
       p2       = p0l(k)
       call linoz_somlfq(p1,p2,f0,f1,f2,ps,f,nl)
       strt0l(k)= f0
       strt1l(k)= f1
       strt2l(k)= f2
    enddo

  end subroutine linoz_strt2m

!--------------------------------------------------------------------

  subroutine linoz_somlfq(p1,p2,f0,f1,f2,ps,f,nl)

  !description: calculates loss freq moments from a set of loss frequencies at std z*, given a ctm model interval pressu pressure range: p1 > p2 (decreasing up)

    IMPLICIT NONE
    ! input/ouput vars
    integer,   intent(in)  :: nl
    real(kind=RKIND),  intent(in)  :: f(nl)
    real(kind=RKIND),  intent(in)  :: ps(nl+1)
    real(kind=RKIND),  intent(in)  :: p1
    real(kind=RKIND),  intent(in)  :: p2
!
! !output parameters:
!
    real(kind=RKIND),  intent(out) :: f0
    real(kind=RKIND),  intent(out) :: f1
    real(kind=RKIND),  intent(out) :: f2
!
! remarks:
! the pressure levels between z* values are:
!      ps(i) > ps(i+1) bounds z*(i)
!                                                                             .
! nl:  z* levels, ==> ps(nl+1) = 0  (extrapolate chemical loss to top)
!      z1 = 16.d0*log10(1000.d0/p1)
!      z2 = 16.d0*log10(1000.d0/p2)
!                                                                             .
! the moments for a square-wave or 'bar': f(x)=f0  b<=x<=c, =0.0 else
!      s0 =   f0 (x)                      [from x=b to x=c]
!      s1 = 3 f0 (x^2 - x)                [from x=b to x=c]
!      s2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
!

    !local variables:
    integer  i
    real(kind=RKIND)   xb,xc,pc,pb,third,sgnf0

    f0     = 0.e+0
    f1     = 0.e+0
    f2     = 0.e+0
    do i = 1,nl
       pc   = min(p1,ps(i))
       pb   = max(p2,ps(i+1))
       if (pc .gt. pb)  then

          ! have condition:  p1>=pc > pb>=p2, 0<=xb < xc<=1
          xc = (pc-p2)/(p1-p2)
          xb = (pb-p2)/(p1-p2)

          ! assume that the loss freq, f, is constant over interval [xb,xc],
          ! f0: (c-b),
          ! f1: 6((c2-c)-(b2-b)),
          ! f2: 5((2c3-3c2+c)-(2b3-3b2+b))
          ! calculate its contribution to the moments in the interval [0,1]
          f0 = f0 +f(i) *(xc -xb)
          f1 = f1 +f(i) *3.e+0*((xc *xc -xc) - (xb *xb -xb))
          f2 = f2 +f(i) *5.e+0* &
               ((xc+xc-1.e+0)*(xc*xc -xc) - &
                (xb+xb-1.e+0)*(xb*xb -xb))
       endif
    enddo

    ! restrain moments: force monotonicity & positive at min end pt

    ! cam: tables can be + or -
    if (f0.ne.0.0) then
       sgnf0=f0 / abs(f0)
    else
       sgnf0=1.0
    endif
    f0=abs(f0)

    !f0 = max(f0, 0.d0)
    third = 1.e+0/3.e+0
    if (f2 .gt. 0.e+0)  then

       ! do not allow reversal of curvature: f2 > 0
       f2   = min(f2, abs(f1)*third, 5.e-1*f0)
       if (f1 .lt. 0.e+0)  then
          f1 = max(-(f0+f2), f1)
       else
          f1 = min(+(f0+f2), f1)
       endif
    else

       ! f2 < 0 = curved down at ends, allow if f1 < f0
       f1  = min(f0,max(-f0,f1))
       f2  = max(f2,(abs(f1)-f0),(-abs(f1)*third))
    endif

    ! cam: apply sign
    f0=sgnf0 * f0
    f1=sgnf0 * f1
    f2=sgnf0 * f2

  end subroutine linoz_somlfq

!---------------------------------------------------------------------

  subroutine linoz_read(state , &
                        linoz_data_file )
    
    !description: read linoz data from data file.

    use iamas_io_units

    IMPLICIT NONE

    ! input/ouput vars
    type(iamas_pool_type), intent(inout) :: state
    character(len=StrKIND), intent(in):: linoz_data_file !path of linoz data file

    !local variables:

    ! scalars
    logical            :: fileexists
    integer            :: iu_file, k,     j
    integer            :: m,       ntbls, ios
    real(kind=RKIND)   :: tmax,    tmin

    ! pointer
    integer,pointer    :: nfields_linoz,nlat_linoz,nlevels_linoz,nmonths_linoz

    ! strings
    character(len=80)  :: heading, titl1
    character(len=255) :: filename
    character(len=255) :: filemsg
    character(len=255) :: thisloc
    character(len=StrKIND):: debugString

    !=================================================================
    ! linoz_read begins here!
    !=================================================================

    ! assume success
    thisloc        = ' -> at linoz_read '

    ! get fields
    call iamas_pool_get_dimension(state, 'linoz_Nlevels', nlevels_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nlat', nlat_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nmonths',nmonths_linoz)
    call iamas_pool_get_dimension(state, 'linoz_Nfields',nfields_linoz)

    allocate(linoz_tparm(nlevels_linoz, nlat_linoz, nmonths_linoz, nfields_linoz))
    
    ! filename
    filename       = trim(linoz_data_file)

    ! test if the file exists
    inquire( file=trim( filename ), exist=fileexists )

    ! test if the file exists and define an output string
    if ( .not. fileexists ) then
       filemsg = 'linoz (linoz_read): required file not found'
       call chemistry_error_fatal('Error:'//trim( filemsg )//trim( filename ))
    endif

    !=================================================================
    ! read climatological data from the linoz lookup tables
    !=================================================================

    ! find a free file lun
    call iamas_new_unit(iu_file)

    ! new std z*=2km levels from model:  z*=10,12,...(25*2)+8 km
    open( iu_file, file=trim( filename ), status='old', &
          form='formatted',      iostat=ios )

    ! return if there was an error opening the file
    if ( ios /= 0 ) then
       call chemistry_error_fatal('Error: can not open'//trim(filename)//trim(thisloc))
    endif

    ! read header
    read ( iu_file, '(a)' ) heading
    call chemistry_debug_message(trim( heading ))

    ! loop over # of fields in the linoz climatology
    do ntbls = 1, nfields_linoz

       ! zero min & max values
       tmin = +1.e+30
       tmax = -1.e+30

       ! skip header line
       read( iu_file, '(a)' ) titl1

       ! loop over # of months in the linoz climatology
       do m = 1, nmonths_linoz

          ! loop over # of latitudes in the linoz climatology
          do j = 1, nlat_linoz

             ! read data into input_opt%linoz_tparm
             read( iu_file, '(20x,6e11.4/(8e11.4))', iostat=ios ) &
                   ( linoz_tparm(k,j,m,ntbls), &
                     k=nlevels_linoz,1,-1                 )

             ! stop on error
             if ( ios > 0 ) then
                call chemistry_error_fatal('error reading ' // trim( filename)//trim(thisloc))
             endif

             ! loop over # of levels in the linoz climatology
             ! and compute the overall min & max
             do k = 1, nlevels_linoz
                tmax = max( tmax, linoz_tparm(k,j,m,ntbls) )
                tmin = min( tmin, linoz_tparm(k,j,m,ntbls) )
             enddo
          enddo
       enddo

       ! write overall min & max
       write (debugString, '(a, a80, 1p, 2e10.3)') '  linoz data:  ', titl1,tmin,tmax
       call chemistry_debug_message(trim( debugString ))

    enddo

    ! echo info
    call chemistry_debug_message('finished reading linoz data')

    ! close the files
    close( iu_file )
    call iamas_release_unit(iu_file)

  end subroutine linoz_read

end module module_linoz
