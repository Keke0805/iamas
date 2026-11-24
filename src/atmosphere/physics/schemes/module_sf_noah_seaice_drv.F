module module_sf_noah_seaice_drv
#if defined(iamas)
use iamas_atmphys_utilities, only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal( M )
#define WRITE_MESSAGE(M) call physics_message( M )
#else
use module_wrf_error
#define FATAL_ERROR(M) call wrf_error_fatal( M )
#define WRITE_MESSAGE(M) call wrf_message( M )
#endif
  use module_sf_noah_seaice
  implicit none
contains
  subroutine seaice_noah( SEAICE_ALBEDO_OPT, SEAICE_ALBEDO_DEFAULT, SEAICE_THICKNESS_OPT, &
       &                  SEAICE_THICKNESS_DEFAULT, SEAICE_SNOWDEPTH_OPT,               &
       &                  SEAICE_SNOWDEPTH_MAX, SEAICE_SNOWDEPTH_MIN,                   &
       &                  T3D, QV3D, P8W3D, DZ8W, NUM_SOIL_LAYERS, DT, FRPCPN, SR,      &
       &                  GLW, SWDOWN, RAINBL, SNOALB2D, QGH, XICE, XICE_THRESHOLD,     &
       &                  ALBSI, ICEDEPTH, SNOWSI,                                      &
       &                  TSLB, EMISS, ALBEDO, Z02D, TSK, SNOW, SNOWC, SNOWH2D, &
       &                  CHS, CHS2, CQS2,                                              &
       &                  RIB, ZNT, LH, HFX, QFX, POTEVP, GRDFLX, QSFC, ACSNOW,         &
       &                  ACSNOM, SNOPCX, SFCRUNOFF, NOAHRES,                           &
       &                  SF_URBAN_PHYSICS, B_T_BEP, B_Q_BEP, RHO,                      &
       &                  IDS, IDE, JDS, JDE, KDS, KDE,                                 &
       &                  IMS, IME, JMS, JME, KMS, KME,                                 &
       &                  ITS, ITE, JTS, JTE, KTS, KTE  )

#if defined(wrfmodel)
#if (NMM_CORE != 1)
    USE module_state_description, ONLY : NOAHUCMSCHEME
    USE module_state_description, ONLY : BEPSCHEME
    USE module_state_description, ONLY : BEP_BEMSCHEME
#endif
#endif

    implicit none

    INTEGER, INTENT(IN)       ::               SEAICE_ALBEDO_OPT
    REAL   , INTENT(IN)       ::               SEAICE_ALBEDO_DEFAULT
    INTEGER, INTENT(IN)       ::               SEAICE_THICKNESS_OPT
    REAL,    INTENT(IN)       ::               SEAICE_THICKNESS_DEFAULT
    INTEGER, INTENT(IN)       ::               SEAICE_SNOWDEPTH_OPT
    REAL,    INTENT(IN)       ::               SEAICE_SNOWDEPTH_MAX
    REAL,    INTENT(IN)       ::               SEAICE_SNOWDEPTH_MIN

    INTEGER, INTENT(IN)       ::                            IDS, &
         &                                                  IDE, &
         &                                                  JDS, &
         &                                                  JDE, &
         &                                                  KDS, &
         &                                                  KDE

    INTEGER, INTENT(IN)       ::                            IMS, &
         &                                                  IME, &
         &                                                  JMS, &
         &                                                  JME, &
         &                                                  KMS, &
         &                                                  KME

    INTEGER, INTENT(IN)       ::                            ITS, &
         &                                                  ITE, &
         &                                                  JTS, &
         &                                                  JTE, &
         &                                                  KTS, &
         &                                                  KTE
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )            , &
         &   INTENT (IN)      ::                            T3D, &
         &                                                 QV3D, &
         &                                                P8W3D, &
         &                                                 DZ8W
         
!Changed by qkyang to speed iAMAS model+++++++++++++++++++++
    REAL,    DIMENSION( its:ite )                              , &
         &   INTENT (IN)      ::                             SR, &
         &                                                  GLW, &
         &                                               SWDOWN, &
         &                                               RAINBL, &
         &                                             SNOALB2D, &
         &                                                 XICE, &
         &                                                  RIB
!Changed by qkyang to speed iAMAS model---------------------

    LOGICAL, INTENT (IN)      :: FRPCPN
    REAL   , INTENT (IN)      :: DT
    INTEGER, INTENT (IN)      :: NUM_SOIL_LAYERS
    REAL   , INTENT (IN)      :: XICE_THRESHOLD

!Changed by qkyang to speed iAMAS model+++++++++++++++++++++
    REAL,     DIMENSION( 1:num_soil_layers, its:ite )          , &
         INTENT(INOUT)   ::                                TSLB
   REAL,    DIMENSION( its:ite ), INTENT(INOUT)::       &
                                                            TSK, &
         &                                                  QGH
    REAL,    DIMENSION( its:ite )                              , &
         &   INTENT (INOUT)   ::                          EMISS, &
         &                                               ALBEDO, &
         &                                                ALBSI, &
         &                                                 Z02D, &
         &                                                 SNOW, &
         &                                                SNOWC, &
         &                                              SNOWH2D, &
         &                                                  CHS, &
         &                                                 CQS2

    REAL,    DIMENSION( its:ite )                              , &
         &   INTENT (OUT)     ::                            HFX, &
         &                                                   LH, &
         &                                                  QFX, &
         &                                                  ZNT, &
         &                                               POTEVP, &
         &                                               GRDFLX, &
         &                                                 QSFC, &
         &                                               ACSNOW, &
         &                                               ACSNOM, &
         &                                               SNOPCX, &
         &                                            SFCRUNOFF, &
         &                                              NOAHRES, &
         &                                                 CHS2

    REAL,    DIMENSION( its:ite )                              , &
         &   INTENT(INOUT)    ::                         SNOWSI

    REAL,    DIMENSION( its:ite )                              , &
         &   INTENT (INOUT)   ::                        ICEDEPTH

    INTEGER, INTENT (IN)      ::               SF_URBAN_PHYSICS
    REAL,    OPTIONAL, DIMENSION( kms:kme, its:ite )  , &
         &   INTENT (INOUT)   ::                        B_Q_BEP, &
         &                                              B_T_BEP
    !Changed by qkyang to speed iAMAS model---------------------
    REAL,    OPTIONAL, DIMENSION(  ims:ime, kms:kme, jms:jme  )  , &
         &   INTENT (IN)      ::                            RHO
    INTEGER :: I
    INTEGER :: J
    REAL    :: FFROZP
    REAL    :: ZLVL
    INTEGER :: NSOIL
    REAL    :: LWDN
    REAL    :: SOLNET
    REAL    :: SFCPRS
    REAL    :: PRCP
    REAL    :: SFCTMP
    REAL    :: Q2
    REAL    :: TH2
    REAL    :: Q2SAT
    REAL    :: DQSDT2
    REAL    :: SNOALB
    REAL    :: TBOT
    REAL    :: SITHICK

    REAL    :: ALBEDOK
    REAL    :: ALBBRD
    REAL    :: Z0BRD
    REAL    :: EMISSI
    REAL    :: T1
    REAL, DIMENSION(1:NUM_SOIL_LAYERS)::  STC
    REAL    :: SNOWH
    REAL    :: SNEQV
    REAL    :: CH
    REAL    :: SNCOVR
    REAL    :: RIBB

    REAL    :: Z0
    REAL    :: ETA
    REAL    :: SHEAT
    REAL    :: ETA_KINEMATIC
    REAL    :: FDOWN
    REAL    :: ESNOW
    REAL    :: DEW
    REAL    :: ETP
    REAL    :: SSOIL
    REAL    :: FLX1
    REAL    :: FLX2
    REAL    :: FLX3
    REAL    :: SNOMLT
    REAL    :: RUNOFF1
    REAL    :: Q1

    REAL    :: APES
    REAL    :: APELM
    REAL    :: PSFC
    REAL    :: SFCTSNO
    REAL    :: E2SAT
    REAL    :: Q2SATI
    INTEGER :: NS
    REAL    :: FDTW
    REAL    :: FDTLIW
    REAL    :: ALBEDOSI
    REAL    :: SNOWONSI
    REAL, PARAMETER  :: CAPA   = R_D / CP
    REAL, PARAMETER  :: A2     = 17.67
    REAL, PARAMETER  :: A3     = 273.15
    REAL, PARAMETER  :: A4     = 29.65
    REAL, PARAMETER  :: A23M4  = A2 * ( A3 - A4 )
    REAL, PARAMETER  :: ROW    = 1.E3
    REAL, PARAMETER  :: ELIW   = XLF
    REAL, PARAMETER  :: ROWLIW = ROW * ELIW

    CHARACTER(len=80) :: message

    FDTLIW = DT / ROWLIW
    FDTW   = DT / ( XLV * RHOWATER )

    NSOIL  = NUM_SOIL_LAYERS

    SEAICE_JLOOP : do J = JTS, JTE
       SEAICE_ILOOP : do I = ITS, ITE

          ! Skip the points that are not sea-ice points.
          IF ( xice(i) < XICE_THRESHOLD ) THEN
              IF ( SEAICE_THICKNESS_OPT == 1 ) THEN
                  icedepth(i) = 0.0
              ENDIF
              IF ( SEAICE_SNOWDEPTH_OPT == 1 ) THEN
                  snowsi(i) = 0.0
              ENDIF
              CYCLE SEAICE_ILOOP
          ENDIF

          SELECT CASE ( SEAICE_THICKNESS_OPT )
          CASE DEFAULT
              WRITE(message,'("Namelist value for SEAICE_THICKNESS_OPT not recognized: ",I6)') SEAICE_THICKNESS_OPT
              FATAL_ERROR(message)
          CASE (0)
              ! Use uniform sea-ice thickness.
              SITHICK = SEAICE_THICKNESS_DEFAULT
          CASE (1)
              ! Use the sea-ice as read in from the input files.
              ! Limit the to between 0.10 and 10.0 m.
              IF ( icedepth(i) < -1.E6 ) THEN
                  WRITE_MESSAGE("Field ICEDEPTH not found in input files.")
                  WRITE_MESSAGE(".... Namelist SEAICE_THICKNESS_OPT=1 requires ICEDEPTH field.")
                  WRITE_MESSAGE(".... Try namelist option SEAICE_THICKNESS_OPT=0.")
                  FATAL_ERROR("SEAICE_THICKNESS_OPT")
              ENDIF
              SITHICK = MIN ( MAX ( 0.10 , icedepth(i) ) , 10.0 )
              icedepth(i) = SITHICK
          END SELECT

          SFCTMP = T3D(I,1,J)
          T1     = tsk(i)
          IF ( SEAICE_ALBEDO_OPT == 2 ) THEN
              IF ( albsi(i) < -1.E6 ) THEN
                  FATAL_ERROR("Field ALBSI not found in input.  Field ALBSI is required if SEAICE_ALBEDO_OPT=2")
              ENDIF
              SNOALB = albsi(i)
              albedo(i) = albsi(i)
              ALBEDOK = albsi(i)
              ALBBRD = albsi(i)
              ALBEDOSI = albsi(i)
          ELSE
              SNOALB = snoalb2d(i)
          ENDIF
          ZLVL   = 0.5 * DZ8W(I,1,J)
          EMISSI = emiss(i)               ! But EMISSI might change in SFLX_SEAICE
          LWDN   = glw(i) * EMISSI        ! But EMISSI might change in SFLX_SEAICE

          ! convert snow water equivalent from mm to meter
          SNEQV = snow(i) * 0.001

          ! snow depth in meters
          SNOWH = snowh2d(i)
          SNCOVR = snowc(i)

          ! Use mid-day albedo to determine net downward solar (no solar zenith angle correction)
          SOLNET = swdown(i) * (1.-albedo(i))   ! But ALBEDO might change after SFLX_SEAICE

          ! Pressure in middle of lowest layer.  Why don't we use the true surface pressure?
          ! Are there places where we would need to use the true surface pressure?
          SFCPRS = ( P8W3D(I,KTS+1,j) + P8W3D(I,KTS,J) ) * 0.5

          ! surface pressure
          PSFC   = P8W3D(I,1,J)

          ! Convert lowest model level humidity from mixing ratio to specific humidity
          Q2     = QV3D(I,1,J) / ( 1.0 + QV3D(I,1,J) )

          ! Calculate TH2 via Exner function
          APES   = ( 1.E5 / PSFC )   ** CAPA
          APELM  = ( 1.E5 / SFCPRS ) ** CAPA
          TH2    = ( SFCTMP * APELM ) / APES

          ! Q2SAT is specific humidity
          Q2SAT  = qgh(i) / ( 1.0 + qgh(i) )
          DQSDT2 = Q2SAT * A23M4 / ( SFCTMP - A4 ) ** 2

          SELECT CASE ( SEAICE_SNOWDEPTH_OPT )
          CASE DEFAULT
              
              WRITE(message,'("Namelist value for SEAICE_SNOWDEPTH_OPT not recognized: ",I6)') SEAICE_SNOWDEPTH_OPT
              FATAL_ERROR(message)

          CASE ( 0 )

              ! Minimum and maximum bounds on snow depth are enforced in SFLX_SEAICE

          CASE ( 1 ) 

              ! Snow depth on sea ice comes from a 2D array, SNOWSI, bounded by user-specified
              ! minimum and maximum values.  No matter what anybody else says about snow 
              ! accumulation and melt, we want the snow depth on sea ice to be specified
              ! as SNOWSI (bounded by SEAICE_SNOWDEPTH_MIN and SEAICE_SNOWDEPTH_MAX).
              SNOWONSI = MAX ( SEAICE_SNOWDEPTH_MIN , MIN ( snowsi(i) , SEAICE_SNOWDEPTH_MAX ) )
              SNEQV = SNOWONSI * 0.3
              snowh2d(i) = SNOWONSI

          END SELECT

          IF ( snow(i) .GT. 0.0 ) THEN
             ! If snow on surface, use ice saturation properties
             SFCTSNO = SFCTMP ! Lowest model Air temperature
             E2SAT = 611.2 * EXP ( 6174. * ( 1./273.15 - 1./SFCTSNO ) )
             Q2SATI = 0.622 * E2SAT / ( SFCPRS - E2SAT )
             Q2SATI = Q2SATI / ( 1.0 + Q2SATI )    ! Convert to specific humidity
             ! T1 is skin temperature
             IF (T1 .GT. 273.14) THEN
                ! Warm ground temps, weight the saturation between ice and water according to SNOWC
                Q2SAT = Q2SAT * (1.-snowc(i)) + Q2SATI * snowc(i)
                DQSDT2 = DQSDT2 * (1.-snowc(i)) + Q2SATI * 6174. / (SFCTSNO**2) * snowc(i)
             ELSE
                ! Cold ground temps, use ice saturation only
                Q2SAT = Q2SATI
                DQSDT2 = Q2SATI * 6174. / (SFCTSNO**2)
             ENDIF
             IF ( ( T1 .GT. 273. ) .AND. ( snowc(i) .GT. 0.0 ) ) THEN   
                ! If (SNOW > 0) can we have (SNOWC <= 0) ?  Perhaps not, so the check on 
                ! SNOWC here might be superfluous.
                DQSDT2 = DQSDT2 * ( 1. - snowc(i) )
             ENDIF
          ENDIF

          PRCP = rainbl(i) / DT

          ! If "SR" is present, set frac of frozen precip ("FFROZP") = snow-ratio ("SR", range:0-1)
          ! SR from e.g. Ferrier microphysics
          ! otherwise define from 1st atmos level temperature

          IF (FRPCPN) THEN
             FFROZP = sr(i)
          ELSE
             IF (SFCTMP <=  273.15) THEN
                FFROZP = 1.0
             ELSE
                FFROZP = 0.0
             ENDIF
          ENDIF

          ! Sea-ice point has deep-level temperature of about -1.8 C
          TBOT = 271.36
          ! TBOT=273.15  ! appropriate value for lake ice.

          ! INTENT(IN) for SFLX_SEAICE, values unchanged by SFLX_SEAICE
          !       I           --
          !       J           --
          !       FFROZP      --
          !       DT          --
          !       ZLVL        --
          !       NSOIL       --
          !       LWDN        --
          !       SOLNET      --
          !       SFCPRS      --
          !       PRCP        --
          !       SFCTMP      --
          !       Q2          --
          !       TH2         --
          !       Q2SAT       --
          !       DQSDT2      --
          !       SNOALB      --
          !       TBOT        --

          Z0BRD  = z02d(i)

          DO NS = 1, NSOIL
             STC(NS) = TSLB(NS,I)
          ENDDO

          CH = chs(i)
          RIBB = rib(i)

          ! INTENT(INOUT) for SFLX_SEAICE, values updated by SFLX_SEAICE
          !       Z0BRD       --
          !       EMISSI      --
          !       T1          --
          !       STC         --
          !       SNOWH       --
          !       SNEQV       --
          !       SNCOVR      --
          !       CH          -- but the result isn't used for anything.
          !                      Might as well be intent in to SFLX_SEAICE and changed locally in 
          !                      that routine?
          !       RIBB        -- but the result isn't used for anything.  
          !                      Might as well be intent in to SFLX_SEAICE and changed locally in 
          !                      that routine?

          ! INTENT(OUT) for SFLX_SEAICE.  Input value should not matter.
          Z0               = -1.E36
          ETA              = -1.E36
          SHEAT            = -1.E36
          ETA_KINEMATIC    = -1.E36
          FDOWN            = -1.E36  ! Returned value unused.  Might as well be local to SFLX_SEAICE ?
          ESNOW            = -1.E36  ! Returned value unused.  Might as well be local to SFLX_SEAICE ?
          DEW              = -1.E36  ! Returned value unused.  Might as well be local to SFLX_SEAICE ?
          ETP              = -1.E36
          SSOIL            = -1.E36
          FLX1             = -1.E36
          FLX2             = -1.E36
          FLX3             = -1.E36
          SNOMLT           = -1.E36
          RUNOFF1          = -1.E36
          Q1               = -1.E36

          call sflx_seaice(I, J, SEAICE_ALBEDO_OPT, SEAICE_ALBEDO_DEFAULT,  &    !C
               &           SEAICE_SNOWDEPTH_OPT, SEAICE_SNOWDEPTH_MAX,      &    !C
               &           SEAICE_SNOWDEPTH_MIN,                            &    !C
               &           FFROZP, DT, ZLVL, NSOIL,                         &    !C
               &           SITHICK,                                         &
               &           LWDN, SOLNET, SFCPRS, PRCP, SFCTMP, Q2,          &    !F
               &           TH2, Q2SAT, DQSDT2,                              &    !I
               &           SNOALB, TBOT, Z0BRD, Z0, EMISSI,                 &    !S
               &           T1, STC, SNOWH, SNEQV, ALBEDOK, CH,              &    !H
               &           ALBEDOSI, SNOWONSI,                              &
               &           ETA, SHEAT, ETA_KINEMATIC, FDOWN,                &    !O
               &           ESNOW, DEW, ETP, SSOIL, FLX1, FLX2, FLX3,        &    !O
               &           SNOMLT, SNCOVR,                                  &    !O
               &           RUNOFF1, Q1, RIBB)

          ! Update our 2d arrays with results from SFLX_SEAICE
          albedo(i)  = ALBEDOK
          emiss(i)   = EMISSI
          tsk(i)     = T1
          z02d(i)    = Z0BRD
          snowh2d(i) = SNOWH
          snowc(i)   = SNCOVR

          ! Convert snow water equivalent from (m) back to (mm)
          snow(i)    = SNEQV * 1000.

          ! Update our ice temperature array with results from SFLX_SEAICE
          DO NS = 1,NSOIL
             TSLB(NS,I) = STC(NS)
          ENDDO

          ! Intent (OUT) from SFLX_SEAICE
          znt(i)    = Z0
          lh(i)     = ETA
          hfx(i)    = SHEAT
          qfx(i)    = ETA_KINEMATIC
          potevp(i) = potevp(i) + ETP*FDTW
          grdflx(i) = SSOIL

          ! Exchange Coefficients
          chs2(i) = cqs2(i)
          IF (Q1 .GT. qsfc(i)) THEN
             cqs2(i) = chs(i)
          ENDIF

          ! Convert QSFC term back to Mixing Ratio.
          qsfc(i)   = Q1 / ( 1.0 - Q1 )

          IF ( SEAICE_SNOWDEPTH_OPT == 1 ) THEN
              snowsi(i) = SNOWONSI
          ENDIF

          ! Accumulated snow precipitation.
          IF ( FFROZP .GT. 0.5 ) THEN
             acsnow(i) = acsnow(i) + PRCP * DT
          ENDIF

          ! Accumulated snow melt.
          acsnom(i) = acsnom(i) + SNOMLT * 1000.

          ! Accumulated snow-melt energy.
          snopcx(i) = snopcx(i) - SNOMLT/FDTLIW

          ! Surface runoff
          sfcrunoff(i) = sfcrunoff(i) + RUNOFF1 * DT * 1000.0

          !
          ! Residual of surface energy balance terms
          !
          noahres(i) = ( SOLNET + LWDN ) &
               &         - SHEAT + SSOIL - ETA &
               &         - ( EMISSI * STBOLT * (T1**4) ) &
               &         - FLX1 - FLX2 - FLX3
               !qkyang changed to support BEP and BEP_BEM model>>>>>>>>>>>>>>>>
! #if defined(wrfmodel)
! #if (NMM_CORE != 1)
        !   IF ( ( SF_URBAN_PHYSICS == NOAHUCMSCHEME ) .OR. &
        !        (SF_URBAN_PHYSICS == BEPSCHEME )      .OR. &
        !        ( SF_URBAN_PHYSICS == BEP_BEMSCHEME ) ) THEN
        !      if ( PRESENT (B_T_BEP) ) then
        !         B_T_BEP(I,1,J)=hfx(i)/dz8w(i,1,j)/rho(i,1,j)/CP
        !      endif
        !      if ( PRESENT (B_Q_BEP) ) then
        !         B_Q_BEP(I,1,J)=qfx(i)/dz8w(i,1,j)/rho(i,1,j)
        !      endif
        !   ENDIF
! #endif
! #endif
          IF ( (SF_URBAN_PHYSICS == 2 ) .OR. ( SF_URBAN_PHYSICS == 3 ) ) THEN
             if ( PRESENT (B_T_BEP) ) then
                B_T_BEP(1,i)=hfx(i)/dz8w(i,1,j)/rho(i,1,j)/CP
             endif
             if ( PRESENT (B_Q_BEP) ) then
                B_Q_BEP(1,i)=qfx(i)/dz8w(i,1,j)/rho(i,1,j)
             endif
          ENDIF
              !qkyang changed to support BEP and BEP_BEM model<<<<<<<<<<<<<<<<<

       enddo SEAICE_ILOOP
    enddo SEAICE_JLOOP

  end subroutine seaice_noah

end module module_sf_noah_seaice_drv
