!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_sfcdiags

CONTAINS

   SUBROUTINE SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2,TH2,Q2,       &
                     PSFC,CP,R_d,ROVCP,CHS,T3D,QV3D,UA_PHYS,       &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
      !Changed by qkyang to speed iAMAS model+++++++++++++++++++++
      REAL,     DIMENSION( its:ite )                             , &
                INTENT(IN)                  ::                HFX, &
                                                              QFX, &
                                                             QSFC

      REAL,     DIMENSION( its:ite )                   , &
                INTENT(IN)                  ::                TSK

      REAL,     DIMENSION( its:ite )                            , &
                INTENT(INOUT)               ::                Q2, &
                                                             TH2, &
                                                              T2
      REAL,     DIMENSION( its:ite )                             , &
                INTENT(IN)                  ::               CHS2, &
                                                             CQS2

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN)                  ::               PSFC

      REAL,     INTENT(IN   )               ::       CP,R_d,ROVCP

! UA changes
      LOGICAL, INTENT(IN) :: UA_PHYS   ! UA: flag for UA option
      REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )            , &
            INTENT(IN   )    ::                           QV3D,T3D
      REAL,     DIMENSION( its:ite )                    , &
                INTENT(IN)                  ::               CHS

! LOCAL VARS
      INTEGER ::  I,J
      REAL    ::  RHO

      DO J=jts,jte
        DO I=its,ite
          RHO = PSFC(I,J)/(R_d * TSK(I))
          if(CQS2(I).lt.1.E-5) then
             Q2(I)=QSFC(I)
          else
              IF ( UA_PHYS ) THEN
                  Q2(I) = QSFC(I) - CHS(I)/CQS2(I)*(QSFC(I) - QV3D(i,1,j))
              ELSE
                  Q2(I) = QSFC(I) - QFX(I)/(RHO*CQS2(I))
              ENDIF
          endif
          if(CHS2(I).lt.1.E-5) then
             T2(I) = TSK(I) 
          else
              IF ( UA_PHYS ) THEN
                  T2(I) = TSK(I) - CHS(I)/CHS2(I)*(TSK(I) - T3D(i,1,j))
              ELSE
                  T2(I) = TSK(I) - HFX(I)/(RHO*CP*CHS2(I))
              ENDIF
          endif
          TH2(I) = T2(I)*(1.E5/PSFC(I,J))**ROVCP
        ENDDO
      ENDDO

  END SUBROUTINE SFCDIAGS
      !Changed by qkyang to speed iAMAS model----------------------
END MODULE module_sf_sfcdiags
