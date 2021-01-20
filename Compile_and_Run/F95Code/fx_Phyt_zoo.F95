C Eutrophication model Inner Oslofjord - File  PHYT_ZOO.FOR
C DIAGNOSTIC OUTPUT CONTROLLED BY PPTEST: number of layers debugged
C for times > TTRIG

      Module fx_Phyt_Zoo
      use ModelParam_RunControl
      use ModelParam_Plankton
      use ModelParam_Topography
      use ModelParam_Decomposition
      use ModelVar_Topography
      use ModelVar_RunControl
      use ModelVar_HydroBioChem
      use fx_RunControl
      use fx_Stoichiometry

      implicit none

$undefine DEBUG_PRPROD
$undefine DEBUG_ZOO
$undefine TEST_STEP
         ! traces time step limitation

$if defined DEBUG_PRPROD
$define Debug_CheckPoint
$define DEBUG_ACTIVE
$endif

$if defined DEBUG_ZOO
$define DEBUG_ACTIVE
$endif

$if defined TEST_STEP
$define DEBUG_ACTIVE
$endif

C ==========================================================

      contains 


      
      SUBROUTINE PRPROD


!      INCLUDE 'EUTRO.INC'




C !!!!!! NB Should rather have death_rate inv. prop. to growth?
C           would mean that old organisms died faster.
C           and perhaps limit population at great depths.



C =============== Local variables:
      real*8 TEMP_L, NO3_L, NH4_L, PO4_L, SiO2_L, OXYG_L
      real*8 NO3DV_L, NH4DV_L, PO4DV_L, SiO2DV_L
      real*8 OXYGDV_L, DOCDV_L, CZOODV_L
      real*8 BACT_L, CZOO_L, SFYT_L
      real*8 CFYT_ACTIVE, CFYT_LG, NFYT_LG, PFYT_LG, CHL_LG
      real*8 FYT_LOSS, FYT_TEMP_FACTOR
      LOGICAL GRAZING
      real*8 ZOO_FILTER_COMPENSATION
      real*8 ZOO_TEMP_FACTOR, ZOO_O2_FACTOR
      real*8 ZOO_FACTOR, ZOO_FACTOR_ABOVE
      integer I, I_B, LSURF, LMAX, L, IS, FG, XG
      real*8 X, Y, G_TEMP, G_TN, R_CRIT, R_CRIT_ADAPT
      real*8 NUTR_MIN
      real*8 N_SURPLUS,     P_SURPLUS,     S_SURPLUS
      real*8 N_SURPLUS_OPT(2), P_SURPLUS_OPT(2), S_SURPLUS_OPT
      real*8 G_D_OPT
      real*8 CFYT_SUM,G_D,CFF
      real*8 nutr_fac, LIGHT_LIMITATION, Uptake_Light_Lim
      real*8 G_RATE, F_RESP, NET_PROD
      real*8 CHL_OPT
      real*8 NH4_UPTAKE, NO3_INHIBIT, NO3_UPTAKE, SiO2_UPTAKE
      real*8 N_UPTAKE_MAX, N_UPTAKE, N_UPTAKE_F
      real*8 PO4_UPTAKE

      real*8 NEXT_AREA, VL
      real*8 FYT_SED_RATE, FYT_DEATH_RATE
      real*8 FROM_ABOVE
      real*8 CONTINUES_TO_NEXT


      real*8 RESPIRED (4)
      real*8 C_RESPIRED, N_RESPIRED, P_RESPIRED
      real*8 S_RESPIRED
      EQUIVALENCE (RESPIRED (1), C_RESPIRED),
     &            (RESPIRED (2), N_RESPIRED),
     &            (RESPIRED (3), P_RESPIRED),
     &            (RESPIRED (4), S_RESPIRED)

      real*8 REMINERALIZED (4)
      real*8 C_REMINERALIZED, N_REMINERALIZED, P_REMINERALIZED
      real*8 Si_REMINERALIZED
      EQUIVALENCE (REMINERALIZED (1), C_REMINERALIZED),
     &            (REMINERALIZED (2), N_REMINERALIZED),
     &            (REMINERALIZED (3), P_REMINERALIZED),
     &            (REMINERALIZED (4), Si_REMINERALIZED)

      real*8 CSED_SUM, NSED_SUM, PSED_SUM, SSED_SUM
      real*8 xSED_SUM(4)
      EQUIVALENCE (xSED_SUM (1), CSED_SUM),
     &            (xSED_SUM (2), NSED_SUM),
     &            (xSED_SUM (3), PSED_SUM),
     &            (xSED_SUM (4), SSED_SUM)

      real*8 C_DET, N_DET, P_DET, S_DET, R_DET
      real*8 x_DET(5)
      EQUIVALENCE (x_DET(1), C_DET),
     &            (x_DET(2), N_DET),
     &            (x_DET(3), P_DET),
     &            (x_DET(4), S_DET),
     &            (x_DET(5), R_DET)


      real*8 C_TO_BOTTOM, N_TO_BOTTOM
      real*8 P_TO_BOTTOM, S_TO_BOTTOM, R_TO_BOTTOM
      real*8 C_SINKING_SED, N_SINKING_SED
      real*8 P_SINKING_SED, S_SINKING_SED, R_SINKING_SED
      real*8 C_ADD_TO_DET, N_ADD_TO_DET
      real*8 P_ADD_TO_DET, S_ADD_TO_DET, R_ADD_TO_DET
      real*8 x_TO_BOTTOM(5), x_SINKING_SED(5), x_ADD_TO_DET(5)
      EQUIVALENCE
     &   (x_TO_BOTTOM  (1), C_TO_BOTTOM  ),
     &   (x_SINKING_SED(1), C_SINKING_SED),
     &   (x_ADD_TO_DET (1), C_ADD_TO_DET ),
     &   (x_TO_BOTTOM  (2), N_TO_BOTTOM  ),
     &   (x_SINKING_SED(2), N_SINKING_SED),
     &   (x_ADD_TO_DET (2), N_ADD_TO_DET ),
     &   (x_TO_BOTTOM  (3), P_TO_BOTTOM  ),
     &   (x_SINKING_SED(3), P_SINKING_SED),
     &   (x_ADD_TO_DET (3), P_ADD_TO_DET ),
     &   (x_TO_BOTTOM  (4), S_TO_BOTTOM  ),
     &   (x_SINKING_SED(4), S_SINKING_SED),
     &   (x_ADD_TO_DET (4), S_ADD_TO_DET ),
     &   (x_TO_BOTTOM  (5), R_TO_BOTTOM  ),
     &   (x_SINKING_SED(5), R_SINKING_SED),
     &   (x_ADD_TO_DET (5), R_ADD_TO_DET )


C ------------------------ PARAMETERS --------------------------

C -------- Stochiometric ratios:  OX_C, OX_NITR, DENITR_C
!      INCLUDE 'STOICHIOM.inc'

      real*8 RAD_ADAPT_COEFF / 40.0/ ! W/m2, coefficient in light adaption
C                  - reduced from 75W/m2 total radiation,
C                    is applied to visible fraction.
      real*8 R_CRIT_MIN  / 5 /   ! W/m2, limit for light adaption,
C                                applies to visible fraction.

      real*8 CFYT_SEED/1.e-12/ ! Helps growth recover from
                               ! accidental negative values, which
                               ! has happened


C Uptake limited by minimum turnover time, set to 1-2 hours :
C (limits problems with large phytoplankton stocks
C  and small nutrient pools)

      real*8 MINIMUM_NH4_TURNOVER_TIME,
     &     MINIMUM_NO3_TURNOVER_TIME,
     &     MINIMUM_PO4_TURNOVER_TIME,
     &     MIN_SIO2_TURNOVER_TIME
      PARAMETER ( MINIMUM_NH4_TURNOVER_TIME = 2./24. ,
     &            MINIMUM_NO3_TURNOVER_TIME = 4./24. ,
     &            MINIMUM_PO4_TURNOVER_TIME = 2/24.  ,
     &            MIN_SiO2_TURNOVER_TIME = 4/24.    )


      real*8 BOTTOM_AREA, SETTLE_AREA, SINKING_VELOCITY

C     ---  Temperature response (Eppley)
C      real*8 TEMP_RESP
C      PARAMETER ( TEMP_RESP = 0.0633 )
C             (erstattet av variable FTRESP(2) og ZTRESP)

C     ----------- PRELIMINARY ZOOPLANKTON MODEL:
      real*8 CFOOD ,NFOOD, PFOOD, COEFF, CFOOD_EFF
      real*8 GRAZING_C, GRAZING_N, GRAZING_P, GRAZING_S, INGEST_C
      real*8 R_MAX, R_ZOO, ZOO_DEATH, ZOO_RESP, ZOO_SINK
      real*8 EXCR_PROD, Zoo_Conditions, Zoo_Conditions_above


C -------- Time step limit:
      real*8 MAX_BIO_TSTEP
      INTEGER step_lim_layer, STEP_LIM_EVENT
      COMMON /BIO_TSTEP/ MAX_BIO_TSTEP, step_lim_layer, step_lim_event
C          coding here depends on DGRADE being performed before PHYT_ZOO
C          since result from DGRADE is used as upper limit in here

      real*8 new_step


C \\\\\\\\\\\\\\\\\\\\\\\\\\\\

$if defined DEBUG_ZOO
      LOGICAL TEST_ZOO
$endif


$if defined DEBUG_ACTIVE
      character*1 ANSWER
      logical ppt, PPT_ON
      PPT_ON = (PPTEST.NE.0) .AND. T.GE. TTRIG
      IF (PPT_ON) THEN
           WRITE( debug_unit, '('' ####### PHYT_ZOO.FOR'')')
      END IF
$endif
C ////////////////////////////


C ==============================================================
C  1. GENERAL PREPARATION:
C ==============================================================

C Start value for biological limit on time step, given by DGRADE

      MXTBIO = MAX_BIO_TSTEP

      DO I_B = 1,NBI
         NFIX(I_B) = 0.0  ! Nitrogen fixation (mg/day)
      enddo

      IF ( .NOT. BIOACT ) RETURN


C Nutrient surplus at optimal ratios:
      DO I=1,2
         N_SURPLUS_OPT(I) = (NCOPT(I)-NCMIN(I))/NCMIN(I)
         P_SURPLUS_OPT(I) = (PCOPT(I)-PCMIN(I))/PCMIN(I)
      END DO
      S_SURPLUS_OPT = (SCOPT-SCMIN)/SCMIN

C Nutrient poor food compensated by increased filtering in zooplankton,
C ensure that response exponent is kept within range 0..1:
      ZOO_FILTER_COMPENSATION = MAX(0.0, Min(1.0, ZFCOMP ) )
C ==================================================================
C  2. LOOP THROUGH BASINS
      DO I_B = 1,NBI
C ==================================================================

$if defined DEBUG_PRPROD
       IF (PPT_ON) THEN
          WRITE( debug_unit, '('' ======= Basin nr.'',I3)') I_B
       END IF
       ppt = ppt_ON !Can be turned off until next basin below
$endif

C ....... Index limits:
       LSURF = INDXI(I_B)+1
       LMAX  = INDXI(I_B+1)

C ....... Initiate sedimentation terms accumulating through layers:
C   ..... Sedimentation flux to bottom:
       CSED_SUM = 0.0
       NSED_SUM = 0.0
       PSED_SUM = 0.0
       SSED_SUM = 0.0

C   ..... Sedimentation flux from layer above (mg/day):
       C_SINKING_SED = 0.0
       N_SINKING_SED = 0.0
       P_SINKING_SED = 0.0
       S_SINKING_SED = 0.0
       ZOO_SINK = 0.0


C   ============================================================
C    b. loop through layers
       DO L = LSURF, LMAX
C   ============================================================

!            IF (l.EQ.1) CALL RC_PRINT (L, 'PHYT_ZOO _1')

$if defined DEBUG_PRPROD

         if (PPT .AND. PPTEST.GT.0 .AND. L.gT.LSURF+PPTEST-1 ) THEN
             PPT=.FALSE.
         ENDIF

         if (PPT_on .AND. PPTEST.LT.0 .AND. L.EQ.LSURF-PPTEST-1 ) THEN
             PPT=.TRUE.
         ENDIF

         IF (PPT) THEN
            WRITE(DEBUG_UNIT,'(''  --------- layer '',I5)') L
         ENDIF

$endif


$if defined DEBUG_ACTIVE
         if (PPT) THEN
            WRITE(*,'('' Continue. diagnostics? (Y[=ENTER]/N'')')
            READ (*,'(A1)' ) ANSWER
            IF ( ANSWER.EQ.'N' .OR. ANSWER .EQ.'n' ) PPT = .false.
         END IF
$endif


C      ......... Volume and areas of current layer:
         VL = VLAYER(L)
C      ......... Area at top of the layer below:
         IF( L .LT. LMAX ) THEN
             NEXT_AREA = AREA(L+1)
         ELSE
            NEXT_AREA = 0.0
         END IF
C      ......... Bottom area in contact with current layer:
         BOTTOM_AREA = (AREA(L)-NEXT_AREA)

C   ---------------------------------------------------------
C   Store array values in single variables for
C   more compact reference later:

         TEMP_L = TEMP(L)
          NO3_L = MAX(0.0D0, NO3(L))
          NH4_L = MAX(0.0D0, NH4(L))
          PO4_L = MAX(0.0D0, PO4(L))
         SiO2_L = MAX(0.0D0, SiO2(L))
         OXYG_L = OXYG(L)
         BACT_L = MAX(0.0D0, BACT(L))
         CZOO_L = MAX(0.0D0, CZOO(L))

$if defined DEBUG_ZOO
         IF (CZOO(L) .LT. 0.0) THEN
            TEST_ZOO = .TRUE.
         ELSE
            TEST_ZOO = .FALSE.
         ENDIF
$endif

$if defined DEBUG_ACTIVE
         IF (PPT) TEST_ZOO=.TRUE.
$endif

         SFYT_L = MAX(0.0D0, SFYT(L))

C   ---------------------------------------------------------
C   Accumulators for local changes, will be updated and
C   restored into arrays when processing of this layer is completed:
         NO3DV_L  = NO3DV(L)
         NH4DV_L  = NH4DV(L)
         PO4DV_L  = PO4DV(L)
         SiO2DV_L = SiO2DV(L)
         OXYGDV_L = OXYGDV(L)
         DOCDV_L  = DOCDV(L)

C ===================================================
C Phytoplankton preparatory loop
C      Sum concentrations for uptake and grazing,
C      and control debug printing in case of negative CFYT:
         CFYT_SUM = 0.0
         CFOOD = BACT_L *GRZBAC
         NFOOD = CFOOD  *NCBACT
         PFOOD = CFOOD  *PCBACT
         DO FG = 1, FYTGRP
C ===================================================

            CFYT_LG = MAX( 0.0D0, CFYT(L,FG))
            NFYT_LG = MAX( 0.0D0, NFYT(L,FG))
            PFYT_LG = MAX( 0.0D0, PFYT(L,FG))
C Sum of fytoplankton, used below in uptake calculations:
            CFYT_SUM = CFYT_SUM + max( cfyt_seed, CFYT_LG )

C Sum of available food concentration
C for zooplankton, weighted by grazing efficiency.
            XG = MIN(2,FG)
            CFOOD = CFOOD + CFYT_LG*GRZFYT (XG)
            NFOOD = NFOOD + NFYT_LG*GRZFYT (XG)
            PFOOD = PFOOD + PFYT_LG*GRZFYT (XG)


C    ============================================
         ENDDO ! (preparatory phytoplankton loop)
C    ============================================





$if defined DEBUG_PRPROD
         if (PPT) THEN
            WRITE(DEBUG_UNIT,*)
     &         ' ------- input values:'
            WRITE( DEBUG_UNIT,'('' in water:'', 3A15)' )
     &         'CONC:','DERIV:','VL*DERIV:'
            WRITE(DEBUG_UNIT,'(1X,A8,'':'',3G15.7)')
     &         'NO3'  , NO3_L  , NO3DV_L  , VL*NO3DV_L ,
     &         'NH4'  , NH4_L  , NH4DV_L  , VL*NH4DV_L ,
     &         'PO4'  , PO4_L  , PO4DV_L  , VL*PO4DV_L ,
     &         'SiO2' , SiO2_L , SiO2DV_L , VL*SiO2DV_L,
     &         'OXYG' , OXYG_L , OXYGDV_L , VL*OXYGDV_L,
     &         'BACT'  , BACT_L  , BACTDV(L)  , VL*BACTDV(L) ,
     &         'CZOO'  , CZOO_L  , CZOODV(L)  , VL*CZOODV(L) ,
     &         'CDET'  , CDET(L) , CDETDV(L)  , VL*CDETDV(L) ,
     &         'NDET'  , NDET(L) , NDETDV(L)  , VL*NDETDV(L) ,
     &         'PDET'  , PDET(L) , PDETDV(L)  , VL*PDETDV(L) ,
     &         'SDET'  , SDET(L) , SDETDV(L)  , VL*SDETDV(L) ,
     &         'SFYT'  , SFYT(L) , SFYTDV(L)  , VL*SFYTDV(L)
            WRITE( DEBUG_UNIT,'('' in sediment:'', 2A15)' )
     &         '..DV1','BOTTOM_AREA*..DV1'
            WRITE( DEBUG_UNIT,'(1X,A12,'':'',2G15.7)')
     &         'CSED'   , CSEDDV(L), BOTTOM_AREA*CSEDDV(L),
     &         'NSED'   , NSEDDV(L), BOTTOM_AREA*NSEDDV(L),
     &         'PSED'   , PSEDDV(L), BOTTOM_AREA*PSEDDV(L),
     &         'RSED'   , RSEDDV(L), BOTTOM_AREA*RSEDDV(L),
     &         'SSED'   , SSEDDV(L), BOTTOM_AREA*SSEDDV(L)
!             IF(DEBUG_UNIT.eq.6) PAUSE
         END IF
$endif

$if defined DEBUG_ZOO
         if (TEST_ZOO) THEN
            if ( CZOO(L).lt. 0.0) 
     &         WRITE(DEBUG_UNIT,*)' NEGATIVE ZOOPLANKTON:'
            WRITE( DEBUG_UNIT,'(10X, 3A15)' )
     &         'CONC:','DERIV:','VL*DERIV:'
            WRITE(DEBUG_UNIT,'(1X,A8,'':'',3G15.7)')
     &         'CZOO(L)'  , CZOO(L)  , CZOODV(L)  , VL*CZOODV(L)
         ENDIF
$endif

C --------------------------------------------------------------
C Microzooplankton; active vertical movements now included.
C --------------------------------------------------------------

C ---- Activity factor given by temperature:
         ZOO_TEMP_FACTOR = EXP (ZTRESP*(TEMP_L-20.0))

C --------- include oxygen level in activity factor,
C           low oxygen levels reduces both growth and respiration,
C           and induces movement into layers with better conditions.
         X = ( OXYG_L - MAX(0.0,ZOXMIN) )
         if ( X .gt. 0.0 ) THEN ! above lower limit:
            ZOO_O2_FACTOR = X/( X + MAX(0.0,ZOXOPT - ZOXMIN) )
            ZOO_FACTOR = ZOO_TEMP_FACTOR * ZOO_O2_FACTOR

C      Death assumed to be effect of grazing from higher trophic levels
C      not included in the model as explicit component. It is assumed
C      that the grazing stops completely at some low concentrations,
C      and is proportional to concentration**2 above this, up to a
C      level where saturation occurs and it is then proportional to
C      concentration. It is implicitly assumed that biomass of higher
C      grazing levels is approximately proportional to total zooplankton
C      biomass, and that grazing pr. grazing biomass is proportional to
C      concentration of prey up to some saturation level.


C      Scaled zooplankton concentration,
C      controls 'Auto-control' by grazing:
            X =  MAX( 0.0D0, CZOO_L-ZCCRIT(1) )
            X =  X / (X + MAX ( 1.e-6, ZCCRIT(2) ) )
                 ! ---> 0 for CZOO_L --> ZCCRIT(1)
                 ! ---> 1 for CZOO_L --> infinity

C      Combine effect of self-control and oxygen-related death:
            ZOO_DEATH =  ZOODR(1)*X*ZOO_FACTOR
     &                 + ZOODR(2)*(1.0-ZOO_O2_FACTOR)*ZOO_TEMP_FACTOR
                             ! 1/day
         ELSE  ! max. oxygen-related reduction, but no grazing
            ZOO_DEATH = ZOODR(2) * ZOO_TEMP_FACTOR
            ZOO_FACTOR = 0.0
         ENDIF
         ZOO_DEATH = ZOO_DEATH*CZOO_L


C ---- respiration from living zooplankton:
C            (grazing dependent. resp. included in growth efficiency)
         ZOO_resp  = CZOO_L * ZRESP* ZOO_FACTOR

C      Accumulator for zoplankton changes due to biological processes,
C      will be combined with transports below.
         CZOODV_L = - ZOO_DEATH - ZOO_RESP

C -------  Sedimentation of dead zooplankton (mg/m3/day),
C          initiates cumulative sedimentation from within layer:

         X = MAX(0.0, Min (1.0, ZDCYCL))
              ! = Fraction immediately recycled

              ! The rest enters detritus with specified
              ! initial degradation rate:
         C_ADD_TO_DET = (1.0-X)*ZOO_DEATH
            R_ADD_TO_DET = C_ADD_TO_DET*DGRATE(1)

!         IF (l.EQ.1) write(555,"(2(A/),4g15.7)")
!     &         ' PRPROD zooplankton ',
!     &         ' R_ADD_TO_DET, C_ADD_TO_DET, DGRATE',
!     &           R_ADD_TO_DET, C_ADD_TO_DET, DGRATE

         N_ADD_TO_DET = C_ADD_TO_DET*NCZOO
         P_ADD_TO_DET = C_ADD_TO_DET*PCZOO
         S_ADD_TO_DET = 0.0
C             ! = sediment flux (will add phytoplankton terms later)

C ------- effect on water concentrations of nutrients:
         X = ZOO_RESP + X*ZOO_DEATH
         OXYGDV_L = OXYGDV_L - X*OX_C
         NH4DV_L  =  NH4DV_L + X*NCZOO
         PO4DV_L  =  PO4DV_L + X*PCZOO
C                 should be differentiated between components?
C                 (could compensate for nutrients lacking in food)
C                 should be dependent on activity?
C                 (can also be included in net assimilation rate)


         X = - 0.5*CZOO_L   ! Always <=0
         IF ( CZOODV_L .lt. X ) then ! avoid generating CZOO<0
              ! Should not happen if X=0,
              ! because then CZOODV_L = 0 too

            NEW_STEP = X/CZOODV_L
            IF ( NEW_STEP.GT.0 .AND. NEW_STEP.LT.mxtbio) THEN
                 !FIRST COND. IS SAFETY MEASURE, NOT REALLY NEEDED

                mxtbio=NEW_STEP

$if defined TEST_STEP
                step_lim_layer = L
                STEP_LIM_EVENT = 1
$endif
          ENDIF
       ENDIF

C -----------------------------------------------------------
C Feeding conditions:

C        ......  Effective food concentration according to nutrient content:
         CFOOD_EFF =  Min (CFOOD, NFOOD/NCZOO, PFOOD/PCZOO)

C     ...... Max. normal ration (before filter_compensation), depends on
C            temperature (Laundry & Lorenzen) & oxygen (in ZOO_FACTOR)
C            and food concentration as carbon:
C              only particulate food triggers activity,
C              but also bacteria are preyed upon:
         if ( CFOOD - BACT_L*GRZBAC .gt. ZCFMIN ) THEN
            R_max = ZOO_FACTOR * ZFMX20 *
     &              MAX(0.0D0, CFOOD /(ZCFSAT + CFOOD) )
C             Grazing limited by biomass availability as carbon
C             (could increase somewhat at high CFYT concentrations
C              if N:C or P:C ratios in phytoplankton where low?)
         else
            R_max = 0.0
         endif


$if defined DEBUG_ZOO
         if (TEST_ZOO)
     &      WRITE( DEBUG_UNIT, '(1X, A12, '':'', G15.7:)' )
     &      'CZOO', CZOO_L,
     &      'ZOO_DEATH', ZOO_DEATH,
     &      'ZOO_resp', ZOO_resp,
     &      'CFOOD', CFOOD,'CFOOD_EFF',CFOOD_EFF,
     &      'R_MAX', R_MAX,
     &      'R_zoo', R_zoo
$endif



C        ....... Grazing control:
         GRAZING =      CZOO_L    .gt. 0.0
     &            .and. CFOOD_EFF .gt. 0.0

C ---- Actual grazing:
         IF ( GRAZING) THEN

C        ...... Available relative ration:
            R_zoo = R_max * (CFOOD / CFOOD_Eff)**Zoo_filter_compensation
C               nutrient-poor food compensated by increased filtering,
C               totally if Zoo_filter_compensation=1,
C               not at all if Zoo_filter_compensation=0.
C               Exponent set above within limits 0..1 from
C               externally specified exponent.

C        ..... Phytoplankton+bacterial concencentration (carbon)
C              removed by grazing pr. time unit:
            GRAZING_C = CZOO_L * R_ZOO                   ! mg/m3/day
C        ...... Effect on bacteria population, with implicit
C               fixed content of N and P:
            X = GRAZING_C * (BACT_L*GRZBAC)/CFOOD        ! mg/m3/day
            BACTDV(L) = BACTDV(L) - X
C        ...... Initiate N, P an Si available in grazed material
C               with contribution from bacteria
            GRAZING_N = X * NCBACT
            GRAZING_P = X * PCBACT
            GRAZING_S = 0.0
         END IF

C    Grazing on each phytoplankton group will be handled in the main
C    phytoplankton loop below, with effect on N, P and Si rations,
C    and with corresponding effects on phytoplankton concentrations.

C    The ingestion/excretion of grazed material will be handled
C    by a final zooplankton section after the phytoplankton loop



C    ================================================================
C    Main loop through phytoplankton groups:
         DO FG = 1,FYTGRP
C    ================================================================


C    ----------------------------------------------------------------
C    Preparations:
            XG = MIN(2,FG)  !index for growth characteristics
C       Store into single variables for effectivity/compact code:
                    ! Fytoplankton measure for passive and loss terms
            CFYT_LG = MAX( 0.0D0, CFYT(L,FG))
                    ! Fytoplankton measure for active terms:
            CFYT_ACTIVE = MAX( CFYT_SEED, CFYT_LG)
                    ! controls uptake & growth, is always >0, this
                    ! helps recovery from accidental small value <0

            NFYT_LG = MAX(0.0D0, NFYT(L,FG))
            PFYT_LG = MAX(0.0D0, PFYT(L,FG))
            CHL_LG  = MAX(0.0D0, CHL(L,FG))

$if defined DEBUG_PRPROD
            IF(PPT) THEN
               WRITE(DEBUG_UNIT,'('' Phytoplankton group '',I2)') FG
               WRITE(DEBUG_UNIT,'(1X,A8,'':'',3G15.7:))')
     &         'CFYT' , CFYT_LG , CFYTDV(L,FG) , VL*CFYTDV(L,FG),
     &         'NFYT' , NFYT_LG , NFYTDV(L,FG) , VL*NFYTDV(L,FG),
     &         'PFYT' , PFYT_LG , PFYTDV(L,FG) , VL*PFYTDV(L,FG),
     &         'CHL'  , CHL _LG , CHLDV (L,FG) , VL*CHLDV (L,FG)
            END IF
$endif


C    -----------------------------------------------------------------
C    Primary production as function of light, temperature
C    and internal nutrient content.

C    ---------- Temperature limited obtainable growth rate:
            FYT_temp_FACTOR = exp(FTRESP(XG)*(TEMP_L-20.0))
            G_TEMP = GMX20(XG)* FYT_TEMP_FACTOR
            GTEMP(L,FG) = G_TEMP


C   ---------- Nutrient limitation:

C   .......... Relative N, P surplus:
C               (via real*8 NUTR_MIN to avoid overflow for small CFYT:
            NUTR_MIN  = CFYT_ACTIVE*NCMIN(XG)
            N_SURPLUS = MIN(N_SURPLUS_OPT(XG), NFYT_LG/NUTR_MIN -1.0)
            NUTR_MIN  = CFYT_ACTIVE*PCMIN(XG)
            P_SURPLUS = MIN(P_SURPLUS_OPT(XG), PFYT_LG/NUTR_MIN -1.0)

$if defined DEBUG_PRPROD
      IF(PPT) THEN
        WRITE( DEBUG_UNIT, * )' ------- Nutrient limitation:'
        WRITE( DEBUG_UNIT,'(3(1X,A9,'':'',G15.7:))')
     &      'NCMIN'  , NCMIN(XG) ,
     &      'N_SURPLUS', N_SURPLUS,
     &      'PCMIN'  , PCMIN(XG) ,
     &      'P_SURPLUS', P_SURPLUS
      END IF
$endif

C   ........... Silicate surplus for diatoms:
            IF (FG.EQ.1) THEN  ! Diatoms with silicate limitation
               NUTR_MIN  = CFYT_ACTIVE*SCMIN
              S_SURPLUS = MIN( S_SURPLUS_OPT, SFYT_L/NUTR_MIN -1.0 )

$if defined DEBUG_PRPROD
              IF(PPT) THEN
                 WRITE( DEBUG_UNIT, '(3( 1X, A9,'':'',G15.7:))')
     &               'SCMIN'  , SCMIN ,
     &               'S_SURPLUS', S_SURPLUS
              END IF
$endif

            END IF

C   ...... Nutrient limitation relative to situation with optimal rates:
            IF ( N_SURPLUS .le. 0.0 .OR. P_SURPLUS .le. 0.0 ) THEN
               nutr_fac = 0.0
            ELSE
C          ...... Growth "resistance" at optimal nutrient ratios:
               G_D_OPT =   1.0
     &                   + 1.0/N_SURPLUS_OPT(XG)
     &                   + 1.0/P_SURPLUS_OPT(XG)
C          ...... Growth "resistance" at current nutrient ratios:
               G_D     = 1.0 + 1.0/N_SURPLUS     + 1.0/P_SURPLUS
               IF (FG.eq.1) THEN   ! Diatoms - silicate dependent
                  IF (S_SURPLUS .le. 0.0 ) THEN
                     nutr_fac = 0.0
                  ELSE
                     nutr_fac =  ( G_D_OPT +1.0/S_SURPLUS_OPT)
     &                        /( G_D     +1.0/S_SURPLUS    )
                  END IF
               ELSE                ! Other groups - only N and P
                 nutr_fac = G_D_OPT/G_D
               END IF
            END IF
C    NOTE! values N_SURPLUS, P_SURPLUS and S_SURPLUS
C          are used to calculate uptake below.


C   ----------- Combined temperature and nutrient limited growth rate:
            G_TN = G_TEMP * nutr_fac
            GTN(L,FG) = G_TN


C  ------------ Light limitation:
C        ........ critical radiation R_Crit:
            X =  ABS(G_TN *CFYT_ACTIVE) ! Max. production
            IF ( CHL_LG .GT. X/1000.0 ) THEN
               R_CRIT = MAX( DBLE(R_Crit_Min/1000), X/CHL_LG )  ! W/m2
C                  lower limit only to avoid problems with RAD=0
C                  in case CFYT=0 in expression for light limitation.
C                  R_Crit should not get much below R_crit_Min
            ELSE
               R_CRIT = 1000
C                  in fact a lower limit to effective quantum yield,
C                  helps the simulation get away from a situation with
C                  chlorophyll measure initiated to zero.
            END IF

C   ------------ Specific growth rate is a function of
C                combined temp/nutrient & light limitation:
            LIGHT_LIMITATION = RAD(L)/(R_CRIT**3+RAD(L)**3)**(1./3.)
            G_RATE =  G_TN * LIGHT_LIMITATION
            GRATE(L,FG) = G_RATE


C    .......... Dark respiration (containing cell death?)
            F_RESP = RESP20*FYT_TEMP_FACTOR
            FRESP(L,FG) = F_RESP
C               (Laws & Bannister 1980, Aruga 1965a,b)
C               Sakshaug et al. 1987 (Mar.Ecol. Prog: Ser., 41,275-282
C            approx 40% reduction in 10 days at low light: rate -0.04/d
C               Light dependent respiration assumed to be included in
C               growth equation, which thus gives net production:


C    .......... Net assimilated carbon production (growth):
            NET_PROD = G_RATE*CFYT_ACTIVE

C    .......... Carbon production excreted as dissolved organic carbon:
C               = EXCRF(1) at nutrient sufficient gross growth,
C               + EXCRF(2) of nutrient limitation of gross growth.
            EXCR_PROD = G_TEMP*LIGHT_LIMITATION
     &                  *( EXCRF(1)
     &                    +EXCRF(2)*(1.0D0-MIN(1.0D0,nutr_fac)))
     &                  * CFYT_ACTIVE
            DOCDV_L = DOCDV_L + EXCR_PROD

C    .......... Oxygen consumption connected to carbon production:
            OXYGDV_L = OXYGDV_L + (NET_PROD+EXCR_PROD)*OX_C
C                effect of nitrate reduction is included below,
C                in connection with nitrate uptake.

C    ---------------------------------------------------------------
C    Light adaption:
C    ......... Optimum chlorophyll/carbon measure:
C                    ( long term adaption )
C                  fixed by amount of light and nutrient limitation
            R_CRIT_ADAPT = max( R_CRIT_MIN,
     &                          RAD_ADAPT_COEFF**0.6 *RAD(L)**0.4 )
            CHL_OPT = G_TN /R_CRIT_ADAPT * CFYT_ACTIVE
C                NOTE:
C                 CHL, CHL_OPT = amount of "chlorophyll" / water volume:

            CHLREL(L,FG) = CHL_LG / MAX(1.D-20,CHL_LG*1.0D-10,CHL_OPT)
C              Optimum used in adaptive term when updating CHLDV below


$if defined DEBUG_PRPROD
      IF(PPT) THEN
        WRITE( DEBUG_UNIT , '(2(1X,A12,'':'',G17.7:))' )
     &      'G_D_OPT',       G_D_OPT,
     &      'nutr_fac',      nutr_fac,
     &      'CHL ',          CHL_LG,
     &      'TEMP',          TEMP_L,
     &      'RAD ',          RAD(L) ,
     &      'G_TEMP',        G_TEMP,
     &      'G_TN',          G_TN,
     &      'R_CRIT',        R_CRIT,
     &      'G_RATE',        G_RATE,
     &      'R_CRIT_ADAPT',  R_CRIT_ADAPT,
     &      'CHL_OPT',       CHL_OPT,
     &      'F_RESP',        F_RESP,
     &      'OXYGDV',        OXYGDV_L
        END IF
$endif

$if defined Debug_CheckPoint
        IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 1'
$endif

C -----------------------------------------------------------------
C Nutrient uptake
C        Currently: Equal maximum specific uptake capacity
C                   for all phytoplankton groups:

C    .......... Light limitation factor of uptake:
            Uptake_Light_Lim = 0.5+0.5*LIGHT_LIMITATION

C    .......... Fraction of phytoplankton carbon in group FG,
C               used to limit total uptake:
            CFF = CFYT_ACTIVE/CFYT_SUM
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 2'
$endif
C    .......... Nitrogen uptake as one way process: ..............
C         .... Uptake capacity as function of concentration,
C              relative rates gN/gC/day:
C              .... Ammonium:
            IF (NH4_L.GT.0.0) THEN
                NH4_UPTAKE = Uptake_Light_Lim
     &                      * VMNH4(XG)*NH4_L/(KSNH4(XG)+NH4_L)
            ELSE
                NH4_UPTAKE = 0.0
            END IF

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 3'
$endif

C               .... Nitrate:
            NO3_INHIBIT = 1.0 / ( 1.0 + (NH4_L/KSNH4(XG))**NH4EXP(XG) )
            IF (NO3_L.GT.0) THEN
                NO3_UPTAKE = NO3_INHIBIT * Uptake_Light_Lim
     &              * VMNO3(XG)*NO3_L/(KSNO3(XG)+NO3_L)
            ELSE
                NO3_UPTAKE = 0.0
            END IF
C         .... Total uptake capacity:
            N_UPTAKE_MAX = (NH4_UPTAKE+NO3_UPTAKE)

C    ........ If capacity is larger than needed for sustaining growth
C             at optimal ratio, reduce uptake down to what is needed
C             as linear function of nutrient surplus:
            N_UPTAKE =   N_UPTAKE_MAX
     &               + MIN( 0.0D0,( G_RATE*NCOPT(XG) - N_UPTAKE_MAX)
     &                          * N_SURPLUS/N_SURPLUS_OPT(XG) )
        
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 4'
$endif

C     ....... Convert to absolute values and distribute on NO3 and NH4,
C             applying limit on turnover:
            IF (N_UPTAKE_MAX .GT. 0.0) THEN
              N_UPTAKE_F = CFYT_ACTIVE*  N_UPTAKE/N_UPTAKE_MAX
              NH4_UPTAKE = MIN ( NH4_UPTAKE * N_UPTAKE_F,
     &                           NH4_L/MINIMUM_NH4_TURNOVER_TIME*CFF )
              NO3_UPTAKE = MIN ( NO3_UPTAKE * N_UPTAKE_F,
     &                           NO3_L/MINIMUM_NO3_TURNOVER_TIME*CFF )
              N_UPTAKE = NH4_UPTAKE + NO3_UPTAKE
            END IF

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 5'
$endif


C    ........ Phosphorus uptake as one way process: ............
C             uptake capacity scaled by relative deficiency,
C             partly by needs for sustaining growth with optimal ratio:

            IF(PO4_L.GT.0.0) THEN
C       .... Uptake capacity as function of concentrations,
C            relative rates gP/gC/day:
               PO4_UPTAKE =   Uptake_Light_Lim
     &                       *VMPO4(XG)*PO4_L/(KSPO4(XG)+PO4_L)
        
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 6'
$endif
C       .... If capacity is larger than needed for sustaining growth
C            at optimal ratio, reduce uptake down to what is needed
C            as linear function of nutrient surplus:
               PO4_UPTAKE =   PO4_UPTAKE
     &             + MIN( 0.0D0, (G_RATE*PCOPT(XG)*PLUXURY-PO4_UPTAKE)
     &                                * P_SURPLUS/P_SURPLUS_OPT(XG) )

          ! Introduced: factor PLUXURY which will allow uptake of up to
          ! PLUXURY times the growth optimal ratio, without affecting
          ! growth limitation


$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 7'
$endif

C       .... Convert to absolute rates, apply limit on turnover:
                PO4_UPTAKE = MIN( PO4_UPTAKE * CFYT_ACTIVE,
     &                            PO4_L/MINIMUM_PO4_TURNOVER_TIME*CFF )
C                                       as mgP/m3/day
            ELSE
                PO4_UPTAKE = 0.0
            END IF

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 8'
$endif

C    .......... Limit timestep to avoid negative nutrients < 0:
            IF (CFF.GT.1.0) THEN
               WRITE (*,*) 'IN phyt_zoo: CFF=', cFF,'>1 !!!!'
            endif

            IF (NO3_UPTAKE .GT. NO3_L ) THEN
               NEW_STEP = CFF*NO3_L/NO3_UPTAKE
               IF ( NEW_STEP.LT.mxtbio) THEN
                   mxtbio=NEW_STEP
$if defined TEST_STEP
                   step_lim_layer = L
                   STEP_LIM_EVENT = 4
$endif
               ENDIF
            ENDIF

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 9'
$endif

            IF (NH4_UPTAKE .GT. NH4_L ) THEN
               NEW_STEP = CFF*NH4_L/NH4_UPTAKE
               IF ( NEW_STEP.LT.mxtbio) THEN
                   mxtbio=NEW_STEP
$if defined TEST_STEP
                   step_lim_layer = L
                   STEP_LIM_EVENT = 5
$endif
               ENDIF
            ENDIF

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 10'
$endif

            IF (PO4_UPTAKE .GT. PO4_L ) THEN
               NEW_STEP = CFF*PO4_L/PO4_UPTAKE
               IF ( NEW_STEP.LT.mxtbio) THEN
                   mxtbio=NEW_STEP
$if defined TEST_STEP
                   step_lim_layer = L
                   STEP_LIM_EVENT = 6
$endif
               ENDIF
            ENDIF
        
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 11'
$endif

$if defined DEBUG_PRPROD
      if ( PPT) THEN
         WRITE( DEBUG_UNIT, '(2(2X,A12,'':'',G15.7:))' )
     &     'NH4'         ,NH4_L      ,
     &     'NO3'         ,NO3_L      ,
     &     'NO3_INHIBIT' ,NO3_INHIBIT ,
     &     'N_UPTAKE_MAX',N_UPTAKE_MAX,
     &     'NH4_UPTAKE'  ,NH4_UPTAKE  ,
     &     'NO3_UPTAKE'  ,NO3_UPTAKE  ,
     &     'N_UPTAKE'    ,N_UPTAKE    ,
     &     'PO4'         ,PO4_L      ,
     &     'PO4_UPTAKE'  ,PO4_UPTAKE ,
     &     'MXTBIO'       ,MXTBIO
!            IF(DEBUG_UNIT.eq.6) PAUSE
       END IF
$endif

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 12'
$endif

C    ... Net change in P and N pool due to uptake:
            NO3DV_L = NO3DV_L - NO3_UPTAKE
            NH4DV_L = NH4DV_L - NH4_UPTAKE
            PO4DV_L = PO4DV_L - PO4_UPTAKE
C     (Effect on phytoplankton at end of loop)


C    ... NOTE!!!! Reduction of NO3 upon uptake affects oxygen balance:
            OXYGDV_L = OXYGDV_L + NO3_UPTAKE*OX_NITR


C --------------------------------------------------------------
C Fytoplankton lysis (death) and respiration, works on CFYT_LG,
C and modifies nutrient uptake working on fytoplankton:


          !     Maximum death rate for given concentration:
            if (CFYT_LG.gt.0.0) then
               FYT_DEATH_RATE = FDEATH(XG)*fyt_temp_factor*
     &                        CFYT_LG/(ABS(FDCSAT(XG))+CFYT_LG)
            else
               FYT_DEATH_RATE = 0.0
            Endif

          !     Minimum fraction of N and P lost in death process:
            Y = MAX(0.0,Min(1.0,FDNUTR(XG))) !(Also used further below)

          !     value of Y modifies death rate:
            FYT_DEATH_RATE = FYT_DEATH_RATE*(1.0 - (1.-Y)*NUTR_FAC)
          !     Always full rate if Y=1 (N& P lost in proportion)
          !        or NUTR_FAC = 0 (growth stopped due to nutr. lim)
          !     Factor 1-NUTRFAC if Y=0 (only nutrient competetion)
          !     Partly dependent on NUTR_FAC for intermediate Y



C Initiate accumulator for net phytoplankton loss rate
C (F_RESP added later):
            FYT_LOSS = FYT_DEATH_RATE
       ! Phytoplankton loss terms affecting C,N,P (and Si) equally
       ! will be accumulated as specific rates in FYT_LOSS.
       ! The combined effect on phytoplankton measures is handled
       ! in the last section of the phytoplankton loop.

          !    Oxidation rate:
            X = FYT_DEATH_RATE + F_resp
          !    Effect on oxygen:
            OXYGDV_L = OXYGDV_L - X*CFYT_LG*OX_C


          !     but the effective value is proportional to nutrient
          !     sufficiency, because the nutrient competiton aspect
          !     is dependent on that:         (uses value of Y above)

          !     Nitrogen:
            X = Y + (1.-Y)*N_SURPLUS/N_SURPLUS_OPT(XG)
            X = X*FYT_DEATH_RATE
                ! X = Specific loss as nutrient:
            NH4DV_L  = NH4DV_L + (X+F_RESP)*NFYT_LG
                ! The rest accounted for as uptake in phytoplankton balance:
            N_UPTAKE   = N_UPTAKE   + (FYT_DEATH_RATE-X)*NFYT_LG
                ! (Will counteract effect of FYT_LOSS on NFYT)

          !     Phosphorus:
            X = Y + (1.-Y)*P_SURPLUS/P_SURPLUS_OPT(XG)
            X = X*FYT_DEATH_RATE
                ! X = Specific loss as nutrient:
            PO4DV_L  = PO4DV_L + (X+F_RESP)*PFYT_LG
                ! The rest accounted for as uptake in phytoplankton balance:
            PO4_UPTAKE = PO4_UPTAKE + (FYT_DEATH_RATE-X)*PFYT_LG
                ! (Will counteract effect of FYT_LOSS on NFYT)

C (Silicate for diatoms: FYT_DEATH_RATE --> sedimenting matter,
C  see below, while respiration does not affect silicate.)


$if defined DEBUG_PRPROD
            IF (PPT) WRITE( DEBUG_UNIT,'(3( 1X, A6, '':'', G15.7:))')
     &        'NO3DV' ,  NO3DV_L ,
     &        'NH4DV' ,  NH4DV_L ,
     &        'PO4DV' ,  PO4DV_L ,
     &        'OXYGDV' , OXYGDV_L
$endif



C ------------------------------------------------------------
            IF (FG.eq.1) then !Silicate uptake and loss for diatoms:
C ------------------------------------------------------------

              IF( SiO2_L .GT. 0.0) THEN
C             .... Uptake capacity as function of concentrations,
C                  relative rates gP/gC/day:
                 SiO2_UPTAKE = Uptake_Light_Lim
     &            *VMSiO2*SiO2_L/(KSSiO2+SiO2_L)

C             .... If capacity is larger than needed for sustaining
C                  growth at optimal ratio, reduce uptake down to what
C                  is needed as linear function of nutrient surplus:
                 SiO2_UPTAKE =   SiO2_UPTAKE
     &                        + MIN(0.0D0, (G_RATE*SCOPT-SiO2_UPTAKE)
     &                                     * S_SURPLUS/S_SURPLUS_OPT )

C             .... Convert to absolute rates, apply limit on turnover:
                 SiO2_UPTAKE = MIN( SiO2_UPTAKE * CFYT_ACTIVE,
     &                              SiO2_L/MIN_SIO2_TURNOVER_TIME )
C                   as mgSi/m3/day
              ELSE
                 SiO2_UPTAKE = 0.0
              END IF
              IF (SiO2_UPTAKE .GT. SiO2_L ) then
                  NEW_STEP = CFF*SiO2_L/SiO2_UPTAKE
                  IF ( NEW_STEP.LT.mxtbio) THEN
                       mxtbio=NEW_STEP
$if defined TEST_STEP
                       step_lim_layer = L
                       STEP_LIM_EVENT = 7
$endif
                  ENDIF
              endif

C    ............ Update state derivatives:
              SiO2DV_L = SiO2DV_L - SiO2_UPTAKE
C                no change in oxygen balance assumed for silicate

$if defined DEBUG_PRPROD
      if ( PPT) THEN
          WRITE ( DEBUG_UNIT,'('' SiO2='',G15.7,'' UPTAKE='',G15.7)' )
     &              SiO2_L, SiO2_UPTAKE
          WRITE( DEBUG_UNIT,'(3( 1X, A6, '':'', G15.7:))')
     &       'SFYTDV',  SFYTDV(L),
     &       'SiO2DV' , SiO2DV_L
          WRITE(DEBUG_UNIT,
     &       '('' L='',I5,'' MXTBIO='',G15.7,'' CFYT_ACTIVE ='',G15.7)')
     &                         L,  MXTBIO, CFYT_ACTIVE
      END IF
$endif

C --------------------------------------------------------------------
            ENDIF
C --------------------------------------------------------------------


$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 13'
$endif

C    ---------- If no phytoplankton, skip metabolic processes.
            IF (CFYT(L,FG) .LE. 0.0) THEN ! zero irrelevant info
               SEDVF(L) = 0.0
               CYCLE ! to next phytoplankton group
            END IF


C --------------------------------------------------------------------
            if (FG.eq.2) Then
                   ! N_fixation and vertical movement in search of
                   ! optimal balance between light and nutrients
                   ! for fytoplankton group 2:
C --------------------------------------------------------------------

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.1'
$endif

C           ......... N fixation:
               X = MIN ( CFYT_active,
     &                   PFYT_LG/PCOPT(XG))*NCOPT(XG) - NFYT_LG
C                = N deficit relative to P or C, as concentration
               X = MAX(0.0D0, X*NFIXRR/365)
                    ! NFIXRR = fraction of deficit compensated pr. year
C                = compensating fixation pr. day
C                  added to N_uptake, balanced as import to system
               N_UPTAKE = N_UPTAKE + X           ! pr. m3 water volume
               NFIX (I_B) = NFIX (I_B) + X*VL ! total import by fixation
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.2'
$endif


C           ..... Net vertical movement (m/day) controlled by difference
C                 between nutrient and light limitation:
               X = LIGHT_LIMITATION - nutr_fac
               iF ( X.GT.0 ) then ! Nutrient limitation dominates
                  X =  F2SINK*(1.0 - nutr_fac)  ! sinking
               ELSE               ! Light limitation dominates
                  X = -F2RIZE*(1.0 - LIGHT_LIMITATION) ! rizing
               ENDIF
C Changed compared to explanation below:
               ! *nutr_fac is skipped
               ! Not sum of both any more.
C                 Sinks if nutrients are strongly limiting, without
C                 regard to light conditions (first part of expression)
C                 Rizes if there is light limitation and the nutrient
C                 limitation is not too strong.

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.3'
$endif


C light_limitation:   0                  0.5                   1
C nutr_fac:+ --------------------------------------------------------
C     0    |        F2SINK              F2SINK               F2SINK
C     0.5  | 0.5(F2SINK-F2RIZE)   0.5(F2SINK-0.5*F2RIZE)   0.5*F2SINK
C     1    |       -F2RIZE           -0.5*F2RIZE               0

               X = X ! * DEPENDED on temperature. SKIPPED
$if defined Debug_CheckPoint
      IF (PPT) THEN
         WRITE(DEBUG_UNIT,'(A,4(I4,A))')
     &     ' PHYT_ZOO Checkpoint FG2.5, X=',X,
     &     ' L=', L, ' LSURF=',LSURF,' LMAX=',LMAX,' ZMID='
         WRITE(DEBUG_UNIT,'(7G11.4)') ZMID
      ENDIF
$endif

C                 Convert to fraction/day removed to next layer (I):
               if ( X.lt.0 .and. L.gt. LSURF ) THEN
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.6a, X=',X
$endif

                  X = ABS(X)/(ZMID(L)-ZMID(L-1))
                  I = L-1
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.6b, X=',X
$endif

               else if (L.lt.LMAX) then
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.7a, X=',X
$endif

                  X = ABS(X)/(ZMID(L+1)-ZMID(L))
                  I = L+1
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.7b, X=',X
$endif

               else
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.8a, X=',X
$endif

                  X = 0.0
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.8b, X=',X
$endif

               endif
        
$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint FG2.9, I,X=',I,X
$endif


C                 The effect is a decrease in this layer:
               FYT_LOSS = FYT_LOSS + X
C                 and a gain in layer I:
               X = X*VL/VLAYER(I)
               CFYTDV(I,FG) = CFYTDV(I,FG) + X*CFYT_LG
               NFYTDV(I,FG) = NFYTDV(I,FG) + X*NFYT_LG
               PFYTDV(I,FG) = PFYTDV(I,FG) + X*PFYT_LG
               CHLDV (I,FG) = CHLDV (I,FG) + X* CHL_LG

C --------------------------------------------------------------------
            END IF
C --------------------------------------------------------------------

$if defined Debug_CheckPoint
      IF (PPT) WRITE(DEBUG_UNIT,*)'PHYT_ZOO Checkpoint 14'
$endif


C --------------------------------------------------------------
C Diatom sedimentation:
C Increased sedimentation in dark weather not accounted for (RIEBESEL)

            IF (FG.eq.1) THEN

C  ..... rate of transformation into sinking diatoms increases as a
C        function of nutrient limitation for nutr_fac < 1./DSNINV,
C        more suddenly the higher DSNEXP is.  Cutoff at low
C        population density:

               X = nutr_fac*DSNINV
               if ( X .lt. 1.0) then
                  X = 1./ max (1.0D-10, X**ABS(DSNEXP) ) - 1.0
               Else
                  X = 0.0
               Endif
                 ! = 0.0 for nutr_fac >= 1/DSNINV
                 ! > 0.0 for smaller nutr_fac

               X = 1.0+X*MAX( 0.0D0,(CFYT_ACTIVE-DSCLIM)/CFYT_ACTIVE)
                 ! = 1.0 always for CFYT_ACTIVE <= DSCLIM


               FYT_SED_RATE = (DSRATE(1) + DSRATE(2)*(X-1))/X
                 ! transition from one value to the other

               SEDVF (L) = FYT_SED_RATE

C                   returned to calling program as information

C  ......... Include in fytoplankton loss:
               FYT_LOSS = FYT_LOSS + FYT_SED_RATE

C     ....... Add to amount of C,N,P settling from this layer (mg/m3/d)

               C_ADD_TO_DET = C_ADD_TO_DET + FYT_SED_RATE*CFYT_LG
               ! Carbon also accumulated as degradation potential:
                  R_ADD_TO_DET = R_ADD_TO_DET +
     &                           FYT_SED_RATE*CFYT_LG*DGRATE(1)

!         IF (l.EQ.1) WRITE (555,"(2(A/),4g15.7)")
!     &         ' PRPROD DIATOM SINKING ',
!     &         ' R_ADD_TO_DET, C_ADD_TO_DET, DGRATE',
!     &           R_ADD_TO_DET, C_ADD_TO_DET, DGRATE

c
               N_ADD_TO_DET = N_ADD_TO_DET + FYT_SED_RATE*NFYT_LG
               P_ADD_TO_DET = P_ADD_TO_DET + FYT_SED_RATE*PFYT_LG

               S_ADD_TO_DET = S_ADD_TO_DET +
     &                        (FYT_SED_RATE+FYT_DEATH_RATE)*SFYT_L
C                 Silicate part assumed to sink upon lysis of cells
C                 (other components remineralized, see above)


$if defined DEBUG_PRPROD
               IF (PPT) THEN
                  WRITE( DEBUG_UNIT,'(2(1X,A14,'':'',G14.7:))')
     &               'SEDVF(L)'    , SEDVF(L)  ,
     &               'C_ADD_TO_DET'     , C_ADD_TO_DET   ,
     &               'N_ADD_TO_DET'     , N_ADD_TO_DET   ,
     &               'P_ADD_TO_DET'     , P_ADD_TO_DET   ,
     &               'S_ADD_TO_DET'     , S_ADD_TO_DET
               END IF
$endif


            END IF  ! ---------- end of sedimentation code

$if defined DEBUG_PRPROD
            IF (PPT) THEN
               WRITE( DEBUG_UNIT,'(3(1X,A10,'':'',G14.7:))')
     &               'FYT_LOSS'   ,  FYT_LOSS
!               IF(DEBUG_UNIT.eq.6) PAUSE
            END IF
$endif

C -----------------------------------------------------------



C Effect of zooplankton grazing on phytoplankton:
            IF (GRAZING) THEN
C      ....... relative grazing rate added to loss term:
               COEFF= GRAZING_C * GRZFYT(XG) / CFOOD
               FYT_LOSS = FYT_LOSS + COEFF
C      ....... accumulate N and P ration:
               GRAZING_N = GRAZING_N + NFYT_LG * COEFF
               GRAZING_P = GRAZING_P + PFYT_LG * COEFF
               IF (FG.EQ.1) THEN
                  GRAZING_S = GRAZING_S + SFYT_L * COEFF
               END IF
            END IF

C --------------------------------------------------------------


C Combined effect of loss terms on phytoplankton concentrations:
            X = FYT_LOSS + F_RESP
            CFYTDV(L,FG) = CFYTDV(L,FG) + NET_PROD   -X*CFYT_LG
            NFYTDV(L,FG) = NFYTDV(L,FG) + N_UPTAKE   -X*NFYT_LG
            PFYTDV(L,FG) = PFYTDV(L,FG) + PO4_UPTAKE -X*PFYT_LG

C        ....... Adaptive change of chlorophyll contents:
            CHLDV(L,FG)  = CHLDV(L,FG)  + G_RATE*( CHL_OPT - CHL_LG )
     &                                  - X*CHL_LG

C ------------------------------------------------------------
C Corresponding changes in silicate for first group:
            IF (FG.EQ.1) THEN
               SFYTDV(L) = SFYTDV(L) + SiO2_Uptake - FYT_LOSS*SFYT_L
            END IF

C Ensure against numeric instabilities due to too fast loss:
            ! (used X value from above)
            IF (X .GT. 0.5) THEN
               NEW_STEP = 0.5 / X
               IF ( NEW_STEP.LT.mxtbio) THEN
                    mxtbio=NEW_STEP
$if defined TEST_STEP
                    step_lim_layer = L
                    STEP_LIM_EVENT = 8
$endif
               ENDIF

            ENDIF

$if defined DEBUG_PRPROD
      IF (PPT) THEN
         WRITE(DEBUG_UNIT,'('' Phytoplankton group '',I2)') FG
         WRITE( DEBUG_UNIT,'('' finally: '', 2A15)' )
     &      'DERIV:','VL*DERIV:'
         WRITE(DEBUG_UNIT,'(1X,A8,'':'',2G15.7:))')
     &         'CFYT' , CFYTDV(L,FG) , VL*CFYTDV(L,FG),
     &         'NFYT' , NFYTDV(L,FG) , VL*NFYTDV(L,FG),
     &         'PFYT' , PFYTDV(L,FG) , VL*PFYTDV(L,FG),
     &         'CHL'  , CHLDV (L,FG) , VL*CHLDV (L,FG)
      END IF
$endif


C ==================================================================
         ENDDO   ! End of loop through phytoplankton groups
C ==================================================================


C ------------------------------------------------------------------
C Microzooplankton:
C    Ingestion and feces pellets production
         IF (GRAZING) THEN

C    ..... Ingested amount of C,N,P limited by ZOOEFF(i),I=1,2,3

            INGEST_C =  MIN( GRAZING_C *ZOOEFF(1)
     &                       ,GRAZING_P / PCZOO *ZOOEFF(2)
     &                       ,GRAZING_N / NCZOO *ZOOEFF(3) )
C             (efficiency could be decreasing for high food conc.)
            CZOODV_L = CZOODV_L + INGEST_C
C               Proportional increase of N and P zooplankton pool
C               in N:C and P:C ratios NCZOO and PCZOO is implicit,
C               and is applied in mass-balance calculations.

C    ..... Noningested part of grazed material:
            C_respired  = GRAZING_C - INGEST_C
            N_respired  = GRAZING_N - INGEST_C * NCZOO
            P_respired  = GRAZING_P - INGEST_C * PCZOO
            S_respired  = GRAZING_S


            DO IS = 1,4
C            ......... partly sedimentation as fecal pellets:
               X_ADD_TO_DET(IS) =  X_ADD_TO_DET(IS)
     &                   + RESPIRED(IS)*(1.0-ZGCYCL)
C            ......... partly directly remineralizated:
               RESPIRED(IS) = RESPIRED(IS)*ZGCYCL
            ENDDO

        ! carbon also accumulated as degradation potential:
            R_ADD_TO_DET = R_ADD_TO_DET +
     &                     C_RESPIRED*(1.0-ZGCYCL)*DGRATE(1)

!         IF (l.EQ.1) WRITE (555,"(2(A/),4g15.7)")
!     &         ' PRPROD FECAL PELLETS ',
!     &         ' R_ADD_TO_DET, C_ADD_TO_DET, DGRATE',
!     &           R_ADD_TO_DET, C_ADD_TO_DET, DGRATE



         ELSE ! No grazing, zero respiration terms
            do IS = 1,4
               respired(IS)=0.0
            enddo
         ENDIF  ! End of zooplankton grazing

$if defined DEBUG_ZOO
           if (TEST_ZOO) THEN
              write( DEBUG_UNIT,*) '------ zooplankton grazing -------'
              WRITE( DEBUG_UNIT, '( 2( 1X, A10, '':'', G15.7:) )' )
     &           'R_max'    ,R_max    ,
     &           'CFOOD' ,CFOOD ,
     &           'R_ZOO'    ,R_ZOO    ,
     &           'GRAZING_C',GRAZING_C,
     &           'GRAZING_P',GRAZING_P,
     &           'GRAZING_N',GRAZING_N,
     &           'GRAZING_S',GRAZING_S,
     &           'INGEST_C' ,INGEST_C ,
     &           'INGEST_P' ,INGEST_C*NCZOO ,
     &           'INGEST_N' ,INGEST_C*PCZOO
            END IF
$endif
       !     Effect on oxygen and nutrient concentrations:
         OXYGDV_L =  OXYGDV_L - C_Respired*OX_C
         NH4DV_L  =  NH4DV_L  + N_Respired
         PO4DV_L  =  PO4DV_L  + P_Respired
         SiO2DV_L = SiO2DV_L  + S_Respired





C -------------------------------------------------------------------
C     Zooplankton migration:
C     Vertical movement to optimize apparent growth conditions, without
C     considering quality of food (unintelligent optimization).
C     Controlled by relative vertical gradient in growth conditions up
C     to value 1/ZMIGRH, where max. migration velocity ZMIGRV (m/day)
C     is achieved:


         ! Find measure of relative fullfillment of growth potential:
         X = ZFMX20 *ZOOEFF(1)
            ! = Maximum relative growth rate at 20 deg C and enough food

         if (CZOO_L.gt.0.0) then
             Zoo_Conditions = (CZOODV_L/CZOO_L)/X
         Else
             Zoo_Conditions = 0
         Endif  ! = Realized part, including effect
                !   of respiration and death


      ! ........... migration downwards toward better conditions:
         if ( L.gt.LSURF ) then

          ! Set X = dimensionless factor on migration velocity:
            X = ZMIGRH * (Zoo_Conditions-Zoo_Conditions_above)/
     &                       (ZMID(L)-ZMID(L-1))
            X = max( -1.0D0, min( 1.0D0, X) )
          !      if large difference over distance <=ZMIGRH
          !      give full migration, weaker gradient
          !      lower migration velocity.


            IF (X .NE. 0.0) THEN

             ! Change to transition rate as water volume/day :
             ! m3/day    m/day  *    m2
               X     = ZMIGRV*X*Area(L)

             ! 1.Limit time step according to relative rate:
             ! 2 Change to transit rate as mgC/day (+/-) = (down/up)
             !   (Activity factor influences migration)

               if ( X .gt. 0.0) then  ! From layer above
                  NEW_STEP = 0.5*VLAYER(L-1)/X
                  IF ( NEW_STEP.LT.mxtbio) THEN
                       mxtbio=NEW_STEP
$if defined TEST_STEP
                       step_lim_layer = L
                       STEP_LIM_EVENT = 2
$endif
                  ENDIF
                  X = X*MAX(0.0D0,CZOO(L-1)-1.e-10)*ZOO_FACTOR_ABOVE
                                ! Turn off for small concentrations

               elseif (X.lt.0.0) then         ! To layer above
                  NEW_STEP = -0.5*VL/x
                  IF ( NEW_STEP.LT.mxtbio) THEN
                       mxtbio=NEW_STEP
$if defined TEST_STEP
                       step_lim_layer = L
                       STEP_LIM_EVENT = 3
$endif
                  ENDIF
                  X = X*max(0.0D0, CZOO_L-1.e-10)*ZOO_FACTOR
               endif


             ! Include effect of migration in concentrations:
               CZOODV(L-1) = CZOODV(L-1) - X/VLAYER(L-1)
               CZOODV_L    = CZOODV_L  + X/VL

            endif

$if defined DEBUG_ZOO
         if (TEST_ZOO)
     &      WRITE( DEBUG_UNIT, '(1X, A12, '':'', G15.7:)' )
     &      'ZOO_Conditions', ZOO_Conditions,
     &      'ZOO_Conditions_above', ZOO_Conditions_above,
     &      'Migration (mg/day):', X,
     &      'CZOODV(L-1):', CZOODV(L-1),'VLAYER(L-1):',VLAYER(L-1),
     &      'CZOODV_L:',CZOODV_L,' VL=',VL
$endif


         endif





       !     Effect of local processes in time derivative
       !     for zooplankton concentration:
         CZOODV(L) = CZOODV(L) + CZOODV_L


       !     Stored for calculation of
       !     migration between this and next layer:
         ZOO_FACTOR_ABOVE = ZOO_FACTOR
         Zoo_Conditions_above = Zoo_Conditions


C ==================================================================
C Distribution of sedimenting matter:

C ------------------------------------------------------------------
C Old code: distributed throughout column instantaneously
C           but with remineralization on the way:

         C_DET = CDET(L)
         N_DET = NDET(L)
         P_DET = PDET(L)
         S_DET = SDET(L)
         R_DET = RDET(L)

C Fluxes x_SINKING_SED from layers above and xD_Sed from current layer
C will partly sediment on bottom within layer, the rest continuing
C as x_SINKING_SED to lower layers (unit mg/d):

C  Effective sedimentation area (reduced due to resuspension)
         SETTLE_AREA = (1.0 - RESUSP)*BOTTOM_AREA

C Store particulate flux downwards as mg/m2/day:
         X = AREA(L)
         CDFLUX(L) = C_SINKING_SED/X
         NDFLUX(L) = N_SINKING_SED/X
         PDFLUX(L) = PDFLUX(l) + P_SINKING_SED/X
                    ! PDFLUX initiated in DGRADE.FOR
         SDFLUX(L) = S_SINKING_SED/X

         SINKING_VELOCITY = SEDVEL(1) + SEDVEL(2)*ZMID(L)

! ............. Update derivatives for dead organic matter:
         DO IS = 1,5  !  using equivalenced arrays (C,N,P,S,R):

         ! ---------- Vertical fluxes
           FROM_ABOVE =   x_SINKING_SED(IS)  ! mg/day

           X = SINKING_VELOCITY  * x_DET(IS) ! mg/m2/day
           CONTINUES_TO_NEXT = X*NEXT_AREA  ! mg/day
           x_TO_BOTTOM(IS) = ( FROM_ABOVE/AREA(L) +X ) * 0.5*SETTLE_AREA
                              ! always = 0 if BOTTOM_AREA=0

         ! Accumulate total sedimentation of C,N,P & S :
           if (IS.le.4) xSED_SUM(IS) = xSED_SUM(IS) + x_TO_BOTTOM(IS)


         ! ---------- Net effect on water concentration of detritus:
           x_ADD_TO_DET(IS) =   x_ADD_TO_DET(IS)
     &       + (FROM_ABOVE - x_TO_BOTTOM(IS) - CONTINUES_TO_NEXT)/VL

         ! ---------- send on to next layer (mg/day):
           x_SINKING_SED(IS) = CONTINUES_TO_NEXT
         ENDDO


! ............... Detritus component updated with additions
!                 from living matter and net effect of sinking:
            CDETDV(L) = CDETDV(L) + C_ADD_TO_DET ! mg/m3/day
            NDETDV(L) = NDETDV(L) + N_ADD_TO_DET
            PDETDV(L) = PDETDV(L) + P_ADD_TO_DET
            SDETDV(L) = SDETDV(L) + S_ADD_TO_DET
            RDETDV(L) = RDETDV(L) + R_ADD_TO_DET

! ............... Accumulated on bottom: update with contribution
!                 from sedimentation of pelagic component:
            if (BOTTOM_AREA.gt.0) then
               X = 1.0/BOTTOM_AREA
               CSEDDV(L) = CSEDDV(L) + C_TO_BOTTOM*X
               NSEDDV(L) = NSEDDV(L) + N_TO_BOTTOM*X
               PSEDDV(L) = PSEDDV(L) + P_TO_BOTTOM*X
               SSEDDV(L) = SSEDDV(L) + S_TO_BOTTOM*X
               RSEDDV(L) = RSEDDV(L) + R_TO_BOTTOM*X
            endif

!            IF (l.EQ.1) CALL RC_PRINT (L,'PHYT_ZOO_2')




$if defined DEBUG_PRPROD
      IF (PPT) THEN
        WRITE( DEBUG_UNIT,'(2(1X,A16,'':'',G14.7:))')
     &     'NEXT_AREA'   , NEXT_AREA ,
     &     'BOTTOM_AREA'    , BOTTOM_AREA  ,
     &     'SETTLE_AREA'    , SETTLE_AREA  ,
     &     'VL'          , VL        ,
     &     'C_ADD_TO_DET'     , C_ADD_TO_DET   ,
     &     'N_ADD_TO_DET'     , N_ADD_TO_DET   ,
     &     'P_ADD_TO_DET'     , P_ADD_TO_DET   ,
     &     'S_ADD_TO_DET'     , S_ADD_TO_DET   ,
     &     'C_TO_BOTTOM'        , C_TO_BOTTOM      ,
     &     'N_TO_BOTTOM'        , N_TO_BOTTOM      ,
     &     'P_TO_BOTTOM'        , P_TO_BOTTOM      ,
     &     'S_TO_BOTTOM'        , S_TO_BOTTOM      ,
     &     'C_SINKING_SED'     , C_SINKING_SED   ,
     &     'N_SINKING_SED'     , N_SINKING_SED   ,
     &     'P_SINKING_SED'     , P_SINKING_SED   ,
     &     'S_SINKING_SED'     , S_SINKING_SED
!        IF(DEBUG_UNIT.eq.6) PAUSE

        write( DEBUG_UNIT,*)
     &          '>>>>>>>layer processed, output values of derivatives:'
        WRITE( DEBUG_UNIT,'('' in water:'', 2A15)' )
     &      'DERIV:','VL*DERIV:'
        WRITE(DEBUG_UNIT,'(1X,A8,'':'',2G15.7:))')
     &     'NO3D',NO3DV_L  , VL*NO3DV_L ,
     &     'NH4D',NH4DV_L  , VL*NH4DV_L ,
     &     'PO4D',PO4DV_L  , VL*PO4DV_L ,
     &     'SiO2',SiO2DV_L , VL*SiO2DV_L,
     &     'OXYG',OXYGDV_L , VL*OXYGDV_L,
     &     'CZOO',CZOODV(L) , VL*CZOODV(L),
     &     'SFYT',SFYTDV(L) , VL*SFYTDV(L)
        WRITE( DEBUG_UNIT,'('' in sediment:'', 2A15)' )
     &          '..DV1','BOTTOM_AREA*..DV1'
        WRITE(DEBUG_UNIT,'(1X,A12,'':'',2G15.7)')
     &     'CSED'   , CSEDDV(L), BOTTOM_AREA*CSEDDV(L) ,
     &     'NSED'   , NSEDDV(L), BOTTOM_AREA*NSEDDV(L) ,
     &     'PSED'   , PSEDDV(L), BOTTOM_AREA*PSEDDV(L)
!        IF(DEBUG_UNIT.eq.6) PAUSE 'Press ENTER to continue'
      END IF
$endif

$if defined DEBUG_PRPROD
      IF (TEST_ZOO) THEN
        WRITE( DEBUG_UNIT,'('' in water:'', 2A15)' )
     &      'DERIV:','VL*DERIV:'
        WRITE(DEBUG_UNIT,'(1X,A8,'':'',2G15.7:))')
     &     'CZOO',CZOODV(L) , VL*CZOODV(L)
      END IF
$endif

C      ........ Restore updated derivatives:
            NO3DV(L)  = NO3DV_L
            NH4DV(L)  = NH4DV_L
            PO4DV(L)  = PO4DV_L
            SiO2DV(L) = SiO2DV_L
            OXYGDV(L) = OXYGDV_L
            DOCDV(L)  = DOCDV_L

C ================================
         END DO ! NEXT LAYER
C ================================


$if defined DEBUG_PRPROD
      if (PPT) THEN
        WRITE(DEBUG_UNIT,*)
     &      '>>>>>>>>>>>>> basin processed, updated imports: '
        WRITE(DEBUG_UNIT,'(3(1X,A8,'':'',G15.7:))')
     &     'CSEDXP'   , CSEDXP(I_B) ,
     &     'NSEDXP'   , NSEDXP(I_B) ,
     &     'PSEDXP'   , PSEDXP(I_B) ,
     &     'SSEDXP'   , SSEDXP(I_B)
      END IF
$endif

C ===============================================================
      END DO  ! NEXT BASIN
C ===============================================================

$if defined TEST_STEP
      IF ( PPTEST.NE.0 .AND. T.GE. TTRIG)
     &     CALL MXTBIO_MESSAGE
     &         ( MXTBIO, DEBUG_UNIT, step_lim_layer, STEP_LIM_EVENT)
$endif

      end subroutine

$if defined TEST_STEP
      SUBROUTINE MXTBIO_MESSAGE(
     &   MXTBIO, DEBUG_UNIT, step_lim_layer, STEP_LIM_EVENT)
      

      real*8 MXTBIO
      INTEGER DEBUG_UNIT, step_lim_layer, STEP_LIM_EVENT

      CHARACTER*16 TEXT(-4:8)  /
     &  'PO4_USED(BACT)','NH4_USED(BACT)',
     &  'NO3_USED(BACT)','DOC_USED(BACT)','DGRADE (BACT)',
     &  'CZOO_DV','CZOO_MIGR_DOWN','CZOO_MIGR_UP',
     &  'NO3_UPTAKE','NH4_UPTAKE','PO4_UPTAKE','SIO2_UPTAKE',
     &  'FYT_LOSS' /

      WRITE( DEBUG_UNIT,'('' MXTBIO = '',G9.2,'' IN LAYER '',I3,
     &                    '' BY '',I2,'':'',A)' ) 
     &     MXTBIO, STEP_LIM_LAYER,
     &     STEP_LIM_EVENT, TEXT(STEP_LIM_EVENT)

      end subroutine
$endif


      
		end Module
