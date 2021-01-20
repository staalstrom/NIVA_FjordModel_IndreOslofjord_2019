C Eutrophication model Inner Oslofjord - File  DGRADE.FOR
      Module fx_Decomposition

      use ModelParam_Plankton
      use ModelParam_Decomposition
      use Modelparam_InitState
      use ModelVar_Topography
      use ModelVar_HydroBioChem
      use ModelVar_RunControl
      use fx_Stoichiometry
      use fx_RunControl
      implicit none

      contains
          
$DEFINE test_illegal_value


$DEFINE IMPORT_BY_COMMON
         ! Change from subroutine with all arguments in subroutine call
         ! to subroutine with most arguments imported by common


C DIAGNOSTIC OUTPUT CONTROLLED BY TEST_MODE TURNED ON (1) OR OFF (0)
C           =1: Debug_print to DEBUG_UNIT

$unDEFINE TEST_STEP

      
$unDEFINE DEBUG_ACTIVE

         ! Extra debug print:
$IF DEFINED DEBUG_ACTIVE
$DEFINE DEBUG_PHOSPHORUS
$DEFINE DEBUG_ASED
$DEFINE DEBUG_LEVEL 1 
$ENDIF
   		! at least some debugging printout




C #####################################################################
$if defined IMPORT_BY_COMMON

      SUBROUTINE DGRADE( TEST_MODE, MLI )
      
      
      LOGICAL TEST_MODE
      integer MLI
!
!      ! Rest of arguments in:
!!      INCLUDE 'eutro.inc'
!
!
!C #####################################################################
!$ELSE
!C #####################################################################
!
!C ====================================================================
!      SUBROUTINE DGRADE( BIOACT, TEST_MODE, ND, DEPTH,
!     &            MLI, NBI, INDXI, VLAYER, AREA,
!     &            TEMP, OXYG, NO3, NH4, PO4,
!     &            DOC, BACT, NCBACT, PCBACT,
!     &            GMX20B, BTRESP, TTURNB,
!     &            DGRATE, DGCMAX, DGNFAC, DGPFAC, DGSFAC, DGWSF,
!     &            DOXBRR,  DOXKB, DOXKM, DOXKS, DOXLIM,
!     &            DNITRR, DNITKS, KOXN, DNITF, DNITXP,
!     &            SULFRR, SULFOX,
!     &            RAMMOX, KAMMOX,
!     &            BURIAL,
!     &            PPAMAX, PPRMAX, PPOXEX, PPOXMX,
!     &            PADRET, PADMAX, PADRLS, PADASD,
!     &            ASEDLR, ASEDOX, ASEDMX,
!     &            CDET, NDET, PDET, SDET, RDET,
!     &            CSED, NSED, PSED, SSED, RSED,
!     &            XSED, XBUR, ASED, PADS,
!     &            PO4DV, NH4DV, NO3DV, OXYGDV, SiO2DV,
!     &            DOCDV, BACTDV,
!     &            CDETDV, NDETDV, PDETDV, SDETDV, RDETDV,
!     &    CSEDDV, NSEDDV, PSEDDV, SSEDDV, RSEDDV,
!     &    XSEDDV, XBURDV, ASEDDV, PADSDV,
!     &    DENITR, CSEDXP, NSEDXP, PSEDXP, SSEDXP, ASEDXP, DGRV )  !=
!      
!
!
!
!C Bacterial utilisation of DOC,
!C and biodegradation of sedimented material,
!C distributed on two degradable fractions.
!
!C ============================================================
!C                     Input arguments
!C ============================================================
!
!      LOGICAL BIOACT, TEST_MODE
!
!C   Dimensions:
!      INTEGER ND
!      INTEGER MLI        !  Sum of number of layers, all basins
!      INTEGER NBI        !  Number of basins
!      INTEGER INDXI(NBI+1)  !  Layer index for basins
!
!C   Topography:
!      real*8 DEPTH(ND)
!      real*8 VLAYER (MLI)
!      real*8 AREA  (MLI)
!
!C   External factors controlling DOC and sediment degradation:
!      real*8 TEMP  (MLI)   !  Temperature in middle of each layer
!      real*8 OXYG  (MLI)   !  Oxygen concentration
!      real*8 NO3   (MLI)   !  Nitrate concentration
!      real*8 NH4(MLI), PO4(MLI)
!      real*8 DOC(MLI)
!
!      real*8 BACT(MLI)
!C   Stochiometric requirements of bacteria:
!      real*8 NCBACT, PCBACT
!
!C   Growth limitations of pelagic bacteria:
!      real*8 GMX20B ! Max specific growth rate at 20 oC
!      real*8 BTRESP ! Temperature dependence exp(BTRESP*(T-20))
!
!      real*8 TTURNB ! Minimum turnover time of
!                  ! DOC, N or P pools by bacteria
!
!
!C   Internal factors controlling degradation:
!
! ! Detritus in water: (mg/m3)
!      real*8 CDET (MLI)  ! Carbon
!      real*8 NDET (MLI)  ! Nitrogen
!      real*8 PDET (MLI)  ! Phosphorus
!      real*8 SDET (MLI)  ! Silicate
!      real*8 RDET (MLI)  ! C*degradation rate
!
! ! On bottom sediments: (mg/m2)
!      real*8 CSED (MLI)  ! Carbon
!      real*8 NSED (MLI)  ! Nitrogen
!      real*8 PSED (MLI)  ! Phosphorus
!      real*8 SSED (MLI)  ! Silicate
!      real*8 RSED (MLI)  ! Carbon*degradation rate
!
!      real*8 ASED (MLI)    ! Oxygen debt in sediments (as sulfide,
!                         ! most of it buffered as iron sulphide.
!
!      real*8 PADS (MLI)    ! Phosphorus adsorbed to sediment mg/m2
!
!
!! Cumulative effect of loss rates for sediment components.
!! Values are used to correct sediment matter accumulated
!! since t=0 into assumed total matter including
!! residual of matter present at t>0.
!      real*8 XSED (MLI)   ! Degradation of carbon, fraction 1.
!      real*8 XBUR (NBI)   ! Effect of burial, common to all components
!        ! Separate terms can be modified & combined:
!        ! Multiplying a rate by a factor corresponds to using
!        ! the factor as an exponent on the X... value.
!        ! Adding two rates corresponds to multiplying X... values.
!
!        ! Resulting X is residual fraction of original amounts that
!        ! should be left at the current time. Under assumption of
!        ! constant yearly pattern of added new material, the material
!        ! accumulated since start of simulation can be multiplied
!        ! by 1/(1-X) to estimate approximate total amount:
!
!        !   Accumulated mass without regard to inital mass:
!        !      M_acc(t) = INTEGRAL[ {Supply(t)-k(t)*M_acc(t)}*dt ]
!        !   Total mass present:
!        !      M_tot(t) = M_acc(t) + M_init*exp[-INTEGRAL(k(t)*dt)]
!        !      M_tot(t) = M_acc(t) + M_init*X(t)
!        !   Assuming that M_tot is approximately equal to M_init:
!        !             (quasi-stationary state)
!        !      M_tot*(1-X(t)) = M_acc(t)
!        !   equivalent to:
!        !      M_tot = M_acc(t)/(1-X(t))
!
!        ! X-->0 means CORR_FAC --> 1.0 as simulation proceeds
!        ! to high T values.
!
!        ! Can reduce (or set to zero ) initial X... values to handle
!        ! cases where simulation starts with realistic
!        ! values >0 of sed. matter
!
!
!
!  ! aerobic degradation:
!      real*8 DOXBRR(2)
!               ! (1): fraction of max. aerobic degradation 
!               !      continuing without macrofauna
!               ! (2): fraction of fanua-related degradation
!               !      which is bacterial 
!      real*8 DGRATE (2) ! Max. specific. degradation rates at T=20 C (1/d)
!                      ! for detritus (1) and for sediment matter (2).
!      real*8 DGCMAX (2) ! Max. absolute gradation rates at T=20 C:
!                      ! 1: as mgC/m3/day for detritus
!                      ! 2: as mgC/m2/day for bottom sediments
!  ! Relative factors for degradation of components relative to C:
!      real*8 DGNFAC     ! Nitrogen
!      real*8 DGPFAC     ! Phosphorus
!      real*8 DGSFAC     ! Silicate
!      real*8 DGWSF      ! Correction factor for conc. of oxygen and NO3
!                      ! When controlling bottom degradation
!                      ! Should be <1.
!      real*8 DOXKB      ! Half saturation O2-conc. for bacterial contrib.
!      real*8 DOXKM      ! Half saturation for macrofanual contribution
!      real*8 DOXKS      ! Coefficient of addition to O2 half.saturation
!                      !   for macrofanual contribution pr. increase of
!                      !   sulphide in sediments

!      real*8 DOXLIM     ! Lower oxygen limit for bottom fauna
!
!      real*8 DNITRR     ! Relative rate of degr. by denitrification
!      real*8 DNITKS     ! NO3 half-saturation concentration for
!                      !   denitrification from external NO3 (µg N/l)
!                      !   without oxygen present in water.
!      real*8 KOXN       ! Oxygen inhibition of denitrification of external
!                      ! nitrate: max. increase of nitrate half-saturation 
!                      ! conc. due to oxic degradation zone as buffer
!      real*8 DNITF      ! Degree of bacterial aerobic degradation that
!                      ! gives 50% nitrification of ammonium
!                      ! released throughout the oxic zone.
!      real*8 DNITXP     ! Exponent in expression for how nitrification
!                      ! and denitrification of released ammonium
!                      ! depends on the extent of aerobic degradation.
!
!      real*8 SULFRR     ! maximum rate of sulfate reduction
!      real*8 SULFOX(2)  ! 1. Oxygen concentration where sulfate reduction
!                      !    becomes possible.
!                      ! 2. Half saturation constant for response
!                      !    to oxygen below limit
!
!      real*8 RAMMOX        ! Max.rate of ammonium oxidation
!      real*8 KAMMOX        ! Oxygen half saturation for ammonium oxidation
!
!      real*8 BURIAL (NBI)  ! Sediment burial rates (1/year)
!
!C Inorganic phosphorus precipitation (with Fe mainly)
!      real*8 PPAMAX        ! abs. max rate (mg/m2/day) at oxygen >=PPOXMX
!      real*8 PPRMAX        ! max. relative rate (m/day)
!      real*8 PPOXEX        ! exponent of dependence on oxygen levels.
!      real*8 PPOXMX        ! max. oxygen level for oxygen dependence
!
!      real*8 PADRET        ! Fraction of remineralized P retained
!C                          in oxic sediments.
!      real*8 PADMAX        ! Maximum limit of retained P (mg/m2)
!      real*8 PADRLS        ! Rel. release rate for excess retained P
!      real*8 PADASD        ! Sulphide content giving full release of
!                         ! P buffered in sediments.
!
!      real*8 ASEDLR(2)     ! Leakage rate of sulphide into sediments (1/year)
!      real*8 ASEDOX        ! Ratio between oxygen content in sediments
!                         ! and oxygen levels in water that gives no change
!                         ! ( = distribution depth of sulphide buffer)
!      real*8 ASEDMX        ! Maximum sulphide content of sediment
!      real*8 ASOXTL        ! (day/m)
!                         ! Conversion from oxygen debt leakage /l/m2/day)
!                         ! to oxygen debt conentration in sediment (l/m3)
!
!C ============================================================
!C                  Output arguments
!C ============================================================
!
!C   Updated derivatives in water (mg/m3, except Oxygen):
!      real*8 PO4DV  (MLI)    !  Nitrate concentration
!      real*8 NH4DV  (MLI)    !  Ammonium (+urea) concentration
!      real*8 NO3DV  (MLI)    !  Nitrate concentration
!      real*8 SiO2DV (MLI)    !  Silicate concentration
!      real*8 DOCDV  (MLI)    !  Dissolved organic carbon
!      real*8 BACTDV (MLI)    !  Pelagic bacteria
!      real*8 OXYGDV (MLI)    !  Oxygen concentration (ml/l = liter/m3)
!      real*8 CDETDV (MLI)  ! Carbon    detritus
!      real*8 NDETDV (MLI)  ! Nitrogen  detritus
!      real*8 PDETDV (MLI)  ! Phosphorus detritus
!      real*8 SDETDV (MLI)  ! Silicate  detritus
!      real*8 RDETDV (MLI)  ! C*degradation rate
!
!C   Initiated derivatives of sediments contents (mg/m2/day)
!      real*8 CSEDDV (MLI)  ! Carbon
!      real*8 NSEDDV (MLI)  ! Nitrogen
!      real*8 PSEDDV (MLI)  ! Phosphorus
!      real*8 SSEDDV (MLI)  ! Silicate
!      real*8 RSEDDV (MLI)  ! C*degradation rate
!      real*8 ASEDDV (MLI)  ! Oxygen debt in sediments (liter O2/m2)
!      real*8 PADSDV (MLI)  ! Phosphorus adsorbed in sediment
!      real*8 XSEDDV (MLI)  ! Residuals of initial amounts.
!      real*8 XBURDV (NBI)  !  "       
!
!C   Net export of C, N, P and Si into sediments from model
!C   (including import of C,N and P from sediments assumed present
!C    at start of simulation.)
!      real*8 CSEDXP (NBI)
!      real*8 NSEDXP (NBI)
!      real*8 PSEDXP (NBI)
!      real*8 SSEDXP (NBI)
!      real*8 ASEDXP (NBI)
!      real*8 DENITR (NBI)    !  Nitrogen removed by denitrification (mg/d)
!
!
!C #####################################################################
!$endif
           ! Subroutine - for both interfaces
C #####################################################################


C ===============================================================
C                       LOCAL VARIABLES
C ===============================================================

       real*8 Days_pr_Year /365./

C -------- Stochiometric ratios:  OX_C, OX_NITR, DENITR_C
!      INCLUDE 'STOICHIOM.INC'

      integer IB, LSURF, LMAX, L, N_FRAC
      real*8 OXYG_L, ASED_L, BACT_L, F_OX_B_REL, F_OX_B_ABS
      real*8 RR, MAX_RATE_C, DEPTH_FACTOR
      real*8 F_OXIC, F_DENITR, F_SULFRED, F_SUM, RB
      real*8 C_REMIN, N_REMIN, P_REMIN, S_REMIN, R_Remin
      real*8 C_remin_tot, C_degr_spes, XB, XC
      real*8 N_Nitr, NH3_release, NH3_removed, NO3_removed, N_removed

      real*8 OXYG_EQUIV, NO3_EFF, OX_ASED_INHIB, R_C_OX_MAX

      real*8 VL, NEXT_AREA, SED_AREA, AREA_ON_VOL, DZ, volume_factor
      real*8 C_loss, N_loss, P_loss, S_loss
      real*8 X, Y, P_Precip, P_net_retained
      real*8 BACT_GROWTH, TEMP_FACTOR, DOC_USED
      real*8 NH4_L, NO3_L, N_Used
      real*8 NH4_USED, NH4_ADDED, one_FDNH3
      real*8 NO3_USED, NO3_ADDED
      real*8 PO4_USED, PO4_ADDED
      real*8 OXYG_USED
      real*8 DOC_EFFICIENCY /0.5/


      real*8 OXYG_DEBT_IN, OXYG_DEBT_OUT
      real*8 R_C_ANOX_MAX, DENITR_C_Corr
      real*8 N_MAX, N_NO3, Q_D_OX, N_OX_NITR, N_ANOX_NITR
      real*8 THETA_S, K_D, H_D, U_S, W_S
      real*8 A_1, A_2, B, BETA_1, BETA_2, E
      real*8 AB_BetaE, A_B, Div, A_Beta, BetaA_BetaE

      real*8 DGRV(20,2)  ! Test values exported to main module 


$if defined DEBUG_ACTIVE
!      INCLUDE 'DEBUG.INC'
      LOGICAL DEBUG_LAYER
$endif

$IF DEBUG_LEVEL == 2
      CHARACTER Error_String*74
      INTEGER ERR_EQ
$endif


C ===============================================================
C                    Time step limit:
C ===============================================================
      real*8 MAX_BIO_TSTEP, MAX_BIO_RATE
      integer step_lim_layer, step_lim_event

      COMMON /BIO_TSTEP/ MAX_BIO_TSTEP, step_lim_layer, step_lim_event

      real*8 new_rate

C          coding here depends on DGRADE being performed before PHYT_ZOO
C          since result here will be used as upper limit in PHYT_ZOO.





C ===============================================================
      ! Statement function which translates X (1-->0)
      ! into correction factor for accumulated sedimented matter:
      real*8 CORR_FAC
      CORR_FAC ( X ) = 1.0/(1.0-X*(1.0-0.2*X))
                     ! Term -0.2*X keeps factor <5 in the start
C ===============================================================




C ===============================================================
C                     Executable part:


$IF DEBUG_LEVEL == 1
      if (TEST_MODE) WRITE(DEBUG_UNIT,*)
     &     ' >>>>>>>>>>>> subroutine DGRADE'
$endif

      MAX_BIO_RATE = 0.001 ! will accumulate max value
                           ! for time step control
$if defined TEST_STEP
      step_lim_event = 0
$endif

      ! FDNH3: Fraction of NH3 in denitrified organic matter which is
      ! converted to molecular nitrogen in combination with nitrate.
      ! complement value, part released as NH4+ :
      one_FDNH3 = max( 0.0, min( 1.0,1.0-FDNH3 ) )


C ===============================================================
C     Zero unused derivatives, since all are integrated in ACSL:
C ===============================================================

      if ( BIOACT ) THEN
         LMAX = INDXI(NBI+1)
      ELSE         ! Biological processes are turned off;
         LMAX = 0  ! just initiate minimally and return.
      ENDIF

      DO L = LMAX+1,MLI
         CSEDDV( L ) = 0.0
         NSEDDV( L ) = 0.0
         PSEDDV( L ) = 0.0
         SSEDDV( L ) = 0.0
         ASEDDV( L ) = 0.0
         PADSDV( L ) = 0.0
      ENDDO


C ================================================================
      DO IB = 1,NBI  ! Process each basin
C ================================================================


$IF DEBUG_LEVEL == 1
         if (TEST_MODE) WRITE(DEBUG_UNIT,*) ' ========= basin:', IB
$ENDIF


         DENITR (IB) = 0.0  ! Initiate integrated transports out of
         CSEDXP (IB) = 0.0  ! system by denitrification or export.
         NSEDXP (IB) = 0.0  ! Export values will below accumulate
         PSEDXP (IB) = 0.0  ! permanent loss into sediment, and will
         SSEDXP (IB) = 0.0  ! eventually include other exports as well.
         ASEDXP (IB) = 0.0

         if (.NOT.BIOACT) CYCLE  ! to next basin


         RB = BURIAL(IB)/Days_pr_year   ! =: 1/day
             !  relative rate of permanent disappearance (burial)
             ! (dimension 1/time = D/(v**2)),
             !  with    D=sediment diffusion
             !          v=sedimentation rate

         ! Reduction factor of initial sedimented matter
         ! due to burial (prepared for differentiating between layers):
         XBURDV(IB) = -RB*XBUR(IB)
         XB = XBUR(IB)


         P_Precip = 0  !  Phosphorus precipitation from layer above
                       !  will be updated for each layer.


         LSURF = INDXI(IB)+1
         LMAX  = INDXI(IB+1)


$IF DEBUG_LEVEL == 2
         if (TEST_MODE) WRITE(DEBUG_UNIT,*) '      RB = ', RB
$ENDIF



C ===============================================================
         DO L = LSURF, LMAX  ! scan layers
C ===============================================================


         ! Volume, area and depth thickness:
            VL = VLAYER(L)
            !  Bottom area of current layer, and area per volume:
            IF( L .LT. LMAX ) THEN
               NEXT_AREA = AREA(L+1)
            ELSE
               NEXT_AREA = 0.0
            ENDIF
            SED_AREA = AREA(L)-NEXT_AREA
            AREA_ON_VOL = SED_AREA/ VL
            DZ = Depth( L-LSurf+2 ) - Depth( L-LSurf+1 )


$IF DEBUG_LEVEL == 1

      DEBUG_LAYER = TEST_MODE

      if (DEBUG_LAYER) THEN
        WRITE( DEBUG_UNIT,'(3(1X,A8,'':'',G15.7:))')
     &     'VL'          , VL,  'DZ',   DZ,
     &     'NEXT_AREA'   , NEXT_AREA,
     &     'SED_AREA'    , SED_AREA,
     &     'AREA_ON_VOL' , AREA_ON_VOL
      ENDIF
$ENDIF


            OXYG_L = MAX(0.0D0,OXYG(L))      ! 02>=0
            NH4_L  = MAX(0.0D0, NH4(L)) ! into single variables
            NO3_L  = MAX(0.0D0, NO3(L)) ! for later reference.
            ASED_L = MAX(0.0D0, ASED(L))
            BACT_L = MAX(0.0D0, BACT(L))

            TEMP_FACTOR = exp(BTRESP*(TEMP(L)-20.0))
                 ! modifies all specified biological rates
                 ! and realisation rate of oxygen demand

C ....... Specified oxygen comsumption:            
            OXYG_USED = ODmRat*TEMP_FACTOR *MAX(0.0D0,ODM(L))
            ODMDev(L) = ODMDev(L) - OXYG_USED

$IF DEBUG_LEVEL == 1
            if (DEBUG_LAYER) THEN
                WRITE( DEBUG_UNIT,'('' -------- layer:'',I5)') L
                WRITE( DEBUG_UNIT,'('' --------------- input :'')' )
                WRITE(DEBUG_UNIT,'(1X,A8,1X,3A15)')
     &            'VARIABLE','DERIVATIVE','DERIV*VOLUME','VALUE'
                WRITE(DEBUG_UNIT,'(1X,A8,'':'',3G15.7:))')
     &            'DOC ', DOCDV(L)  , VL*DOCDV(L),   DOC(L),
     &            'BACT', bactDV(L) , VL*bactDV(L),  bact(L),
     &            'OXYG', OXYGDV(L) , VL*OXYGDV(L),  OXYG(L),
     &            'NO3' , NO3DV(L)  , VL*NO3DV(L) ,  NO3(L),
     &            'NH4' , NH4DV(L)  , VL*NH4DV(L),   NH4(L),
     &            'PO4' , PO4DV(L)  , VL*PO4DV(L),   PO4(L)
                WRITE(DEBUG_UNIT,'(1X,A8,'':'',2G15.7:))')
     &            'SiO2', SiO2DV(L) , VL*SiO2DV(L)
                WRITE(DEBUG_UNIT,'(2(1x,A,G11.5))')
     &        'TEMP', TEMP(L), 'temp_factor', TEMP_FACTOR
            ENDIF
$ENDIF


C ---------------------------------------------------------------
C      Pelagic bacteria using DOC

            !      Was assumed to operate only in oxic conditions,
            !      prevents oxygen getting below 0.0 except due to
            !      sulphate reduction in sediments.
            !  if (OXYG_L .gt. 0.0) then


            ! --------- Bacteria growth at 20 oC:
               ! ------- limited by C, N or P turnover
               !         ( with 50% efficiency on nitrate uptake):
               BACT_GROWTH =  MIN( DOC(L)*DOC_Efficiency,
     &                             (NH4_L + NO3_L/2.0)/NCBACT,
     &                             DBLE(PO4(L)/PCBACT)
     &                           ) / TTURNB
               ! ------- or by specific growth rate:
               Bact_Growth = MIN (Bact_Growth,GMX20B*BACT_L)
               Bact_Growth = MAX( 0.0D0, Bact_Growth )

            ! ---------- at ambient temperature and oxygen:
               Bact_Growth = Bact_Growth*Temp_Factor

            !else
            !   Bact_Growth = 0.0
            !endif



         ! Transit to detritus:
            X = BACDET * BACT_L

         ! --------- total change:
            BACTDV(L) =  BACTDV(L) + Bact_Growth - X

         ! ---------- effects on organic carbon pools and on
         !            nitrogen, phosphorus and oxygen:
            DOC_USED =  Bact_Growth/DOC_Efficiency

            DOCDV(L)  =  DOCDV(L) - DOC_USED

            CDETDV(L) = CDETDV(L) + X
            RDETDV(L) = RDETDV(L) + X*DGRATE(1)
            NDETDV(L) = NDETDV(L) + X*NCBACT
            PDETDV(L) = PDETDV(L) + X*PCBACT

            if (Bact_Growth .gt. 0) then
                NEW_RATE = DOC_USED/DOC(L)
                IF ( NEW_RATE.GT.MAX_BIO_RATE ) THEN
                   Max_bio_RATE=NEW_RATE
$if defined TEST_STEP
                   step_lim_layer = L
                   STEP_LIM_EVENT = -1
$ENDIF
                ENDIF
            ENDIF

            N_USED = Bact_Growth * NCBACT
            NH4_USED = min( N_used, NH4_L*Temp_Factor/TTURNB )
            NO3_USED = N_USED - NH4_USED ! (primarily assimilate ammonium)
                 ! Note !  No pelagic denitrification included yet.

            PO4_USED = Bact_Growth*PCBACT

            OXYG_USED = OXYG_USED + (DOC_USED - Bact_Growth)*OX_C
     &                  - NO3_USED*OX_NITR
                 ! Carbon burned off consumes oxygen equivalents,
                 ! and oxygen content influenced by nitrate reduction
                 ! in bacterial uptake: (Only during oxic conditions)


$IF DEBUG_LEVEL == 2
      if (DEBUG_LAYER) THEN
         WRITE( DEBUG_UNIT,'(3(1X,A11,'':'',G13.7))')
     &      'BACT_GROWTH', BACT_GROWTH,
     &      'DOC_USED', DOC_USED,'OXYG_USED', OXYG_USED,
     &      '     NO3_USED', NO3_USED,
     &      'NH4_USED', NH4_USED, 'PO4_USED', PO4_USED
      ENDIF
$ENDIF



! -----------------------------------------------------------------
!  Free NH4 is transformed to NO3 under oxic conditions,
!  (not related to sediments):
! -----------------------------------------------------------------            

            if (OXYG_L .gt. 0.0) then
               X = RammOX*TEMP_FACTOR * NH4_L*OXYG_L/(KammoX+OXYG_L)
            Else
               X = 0.0
            Endif
            NH4_USED  = NH4_USED + X
            NO3_ADDED = X
            OXYG_USED =  OXYG_USED + OX_NITR*X


C ----------------------------------------------------------------
C       Leak-out of oxygen debt from sediments:
C ----------------------------------------------------------------

        ! Leaks out of sediment with a rate which depends on
        ! excess sulphide vs. saturation of sediment:

            X = ASED(L) + ASEDOX*OXYG(L) 
              ! (liter O2 /m2 = m*(ml/litre))
				  ! ==> ASEDOX has unit meter:
              !     and can be considered to represent
              !     distribution depth in sediment
              !     of active sulphide buffer

            if ( X .le. ASEDMX ) THEN
               OXYG_DEBT_OUT = MIN(0.0d0,X*ASEDLR(1)/365.)
            ELSE
               OXYG_DEBT_OUT = ( ASEDMX*ASEDLR(1)
     &                           + (X-ASEDMX)*ASEDLR(2) )/365.0
            ENDIF
                      ! liter O2/m2/day
                ! ASEDLR = specific leakage rate (1/year), of
                !          amount up to and above limit,
                !          use precipitation and diffusion
                !          consideration to set value.
                ! ASEDOX = Ratio between oxygen debt in sediments
                !    and water oxygen concentration affecting
                !    transport equally. Also equal to ratio between
                !    negative ASED (oxygen in sediments) and oxygen
                !    content in water at equilibrium,
                !  = height of water column containing same amount of
                !    oxygen as in fully oxygenated sediments, thus
                !    can be set approximtely = the depth
                !    of sediments taken into consideration.

$if defined DEBUG_ACTIVE
$if defined DEBUG_ASED
      if (DEBUG_LAYER) THEN
        WRITE ( DEBUG_UNIT, '(''Before sulphide in/out of sediment:'')')
        WRITE ( DEBUG_UNIT, '(1X,2(A16,'':'',G15.7:))')
     &     'ASED(L)', ASED(L)
     &    ,'OXYG_USED',OXYG_USED
     &    ,'OXYG_DEBT_OUT',OXYG_DEBT_OUT
     &    ,'ASEDOX', ASEDOX
     &    ,'ASEDMX', ASEDMX
     &    ,'ASEDLR(1)', ASEDLR(1),' (2)',ASEDLR(1) 
     &    ,'Excess sulphide', X
      ENDIF
$ENDIF
$ENDIF
            OXYG_USED = OXYG_USED + OXYG_DEBT_OUT * AREA_ON_VOL
                !liter/m3/day


       ! ----------- Initiate phosphorus and ammonium increase:
            NH4_ADDED = 0.0
            PO4_ADDED = 0.0



C ------------------------------------------------------------------
C           DEGRADATION OF ORGANIC MATTER:
C ------------------------------------------------------------------


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DO N_FRAC = 1,2  ! Detritus or sediment
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               SELECT CASE (N_FRAC)
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
               CASE (1)
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%

$IF DEBUG_LEVEL == 1
                 if (DEBUG_LAYER) THEN
                     WRITE( DEBUG_UNIT,*) ' *** Detritus degradation:'
                 ENDIF
$ENDIF
            ! Organic matter as mg/m3:
                 C_remin = max(0.0D0,CDET(L))
                 N_remin = max(0.0D0,NDET(L))
                 P_remin = max(0.0D0,PDET(L))
                 S_remin = max(0.0D0,SDET(L))
                 R_remin = max(0.0D0,RDET(L))

            ! Additional oxygen consumption: 
                 C_degr_spes = 0.0

            ! Effective limiting concentrations:
                 OXYG_EQUIV = OXYG(L)  ! Net oxygen debt included
                 
                 NO3_EFF = NO3_L
            
            ! Sulfide in sediments has reduced inhibiting effect,
            ! (weight factor DGWSF and area ratio):
                 OX_ASED_INHIB = 1.0 + DOXKS*ASED_L
     &                                 *DGWSF*SED_AREA/AREA(L)
            ! Depth dependence for sinking detritus:
                 DEPTH_FACTOR = ZMID(L)/(DGDETZ+ZMID(L))

            ! NO CONVERSION TO GET RATE PR. VOLUME: 
                 VOLUME_FACTOR = 1.0

            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               CASE (2)
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

$IF DEBUG_LEVEL == 1
                 if (DEBUG_LAYER) THEN
                     WRITE( DEBUG_UNIT,*) ' *** Sediment degradation:'
                 ENDIF
$ENDIF
            ! Organic matter as mg/m2:
                 C_remin = MAX(0.0D0,CSED(L))
                 N_remin = MAX(0.0D0,NSED(L))
                 P_remin = MAX(0.0D0,PSED(L))
                 S_remin = MAX(0.0D0,SSED(L))
                 R_remin = MAX(0.0D0,RSED(L))

            ! Additional oxygen consumption, 
            ! specified carbon degradation in sediments : 
                 X = DEPTH(L-LSURF+2)-CDRDEPTH(IB)
                 IF (X.gt.0.0) then
                    C_degr_spes = MAX(0.0,CDRSED(IB))*min(1.0D0,X/DZ) 
                 ELSE 
                    C_degr_spes = 0.0
                 ENDIF

            ! Effective limiting concentrations:
                 OXYG_EQUIV = OXYG(L)*DGWSF - OXYG_DEBT_OUT*ASOXTL 
                 ! Note: in current model formulation: 
                 ! H2S conc. >0 affects processes,  
                 ! so negative OXYG_EQUIV has to be noted 

                 NO3_EFF =  NO3_L* DGWSF
                    ! (Oxygen and NO3 concentrations changed by factor
                    !  DGWSF in saturation expressions below, to take
                    !  account of transport barrier water <--> bottom
                    !  DGWSF should be < 1.0
            
            ! Full inhibitory effect of sulphide in sediments:
                 OX_ASED_INHIB = 1.0 + DOXKS*ASED_L
                    ! Anoxic bottom only affects degradation in sediment
                 
            ! No explicit depth dependence for sediment:
                 DEPTH_FACTOR = 1.0
            
            ! CONVERSION FROM RATE PR. AREA TO RATE PR. VOLUME: 
                 VOLUME_FACTOR = AREA_ON_VOL


            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               END SELECT
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      ! ------ Convert to max. degradation pr. volume or area
      !        under in situ temperatures and with excess oxygen:
      !        assumed under no kind of external kinetic/transport
      !        limitations: intrinsic degradation potential.

      ! DGRATE can be adjusted at this point to take account of
      ! reduction in degradability of organic matter:

      ! Effective specific C rate:
               if ( C_remin.gt.0.0 ) then
                   RR = R_REMIN/C_REMIN ! By def of RDET or RSED
               else
                   RR = 0.0
               endif
               MAX_RATE_C = RR*TEMP_FACTOR*DEPTH_FACTOR

               C_remin = C_remin*MAX_RATE_C
               N_remin = N_remin*MAX_RATE_C*DGNFAC
               P_remin = P_remin*MAX_RATE_C*DGPFAC
               S_remin = S_remin*MAX_RATE_C*DGSFAC

      ! Reduction rate for R... (DET or SED) is increased to give a
      ! higher specific rate than for C..., the effect is a steady
      ! decrease of degradation rates of residual C.

               X = max ( 1.D-10, MAX_rate_C )  ! < DGRATE
               R_remin = MAX_RATE_C*R_remin*
     &                  (1.0D0 + ACCLRC * X**ACCLXP )

!               if (L.eq.1) then
!                  WRITE(*,"(1x,4A15)") 
!                       'C_degr_spes', 'C_remin', 'R_remin','Max-rate_C'
!                  WRITE(*,*) C_degr_spes, C_REMIN, R_REMIN, max_rate_c
!                  PAUSE
!               endif

       
       ! Add specified background degradation of carbon:       
               C_remin_tot = C_remin + C_degr_spes

		   !--------------------------------------------------------------------
		   ! Reduction factors for aerobic degradation,
		   ! composed of bacterial and fauna contribution:
                  ! ( Bacterial oxygen requirement could depend on
                  !   amount of organic matter to be degraded )

                  ! Bacterial oxic degradation:
            if (OXYG_EQUIV.gt.0.0) then
               if ( DOXKB .gt. 0.0 ) then
                  F_OX_B_REL = OXYG_EQUIV/(DOXKB+OXYG_EQUIV)
               else
                  F_OX_B_REL = 1.0
               endif
            else
               F_OX_B_REL = 0.0
            endif

                  ! Fauna-related degradation as function of
                  ! oxygen content above lower limit DOXLIM,
                  ! with half-saturation concentrartion DOXKM

            Y = MAX( 0.0D0, OXYG_L - DOXLIM )
            if ( DOXKM .gt. 0.0 ) then
               X = Y/(DOXKM+Y)
            elseif (Y.gt.0.0) then
               X = 1.0
            else
               X = 0.0
            endif

                  ! Further inhibition due to sulphide in sediment:

            X = X/OX_ASED_INHIB

            Y = (1-X)*DOXBRR  ! Fraction of max. degradation
                              ! continuing without fanua 

            F_OX_B_ABS = Y*F_OX_B_REL
            F_OXIC = F_OX_B_ABS + X
                  ! total factor on rate, always 0 <= F_OXIC < =1


      ! --------- The amount of aerobic degradation is limited by oxygen:
            R_C_OX_MAX = F_OXIC * C_remin_tot ! mg/m3/day or mg/m2/day

      !           but also by absolute limit pr. volume or area:
            R_C_OX_MAX =   R_C_OX_MAX * DGCMAX(N_FRAC)/
     &                   ( R_C_OX_MAX + DGCMAX(N_FRAC) )


		   ! ---------------------------------------------------------------
		   !                  Anoxic degradation:
		   ! ---------------------------------------------------------------

	   ! Measure of available organic matter for anoxic degradation:
            R_C_ANOX_MAX = Y*(1.0-F_OX_B_REL)*C_remin_tot
                ! ( Anoxic degradation should not replace
                !   the assumed fauna-related part )
			
			! KUTTET UT:
			!            if (F_OXIC .gt. 0) then
			!                R_C_ANOX_MAX = R_C_ANOX_MAX*F_OX_B_ABS/F_OXIC
			!            endif
			               
               !  Fauna contribution increases aerobic degradation,
               !  and also reduces possible anoxic degradation.
               !  (Implement ordinary saturation curves for for F_Oxic ?


	   ! Nitrified fraction of ammonium released in oxic zone:
$if defined test_illegal_value          
            N_OX_NITR = (F_OX_B_REL/(F_OX_B_REL+DNITF))
            if ( N_OX_NITR .lt.0.0 ) then
               write (6,*)'N_OX_NITR =', N_OX_NITR,' <0' 
!               Pause
            endif
            N_OX_NITR = N_OX_NITR**DNITXP
$ELSE
            N_OX_NITR = (F_OX_B_REL/(F_OX_B_REL+DNITF))**DNITXP
$endif
               ! Factor describing how large part of ammonium generated
               ! in oxic zone that can be nitrified before being
               ! released to the water. (0.5 because generation of
               ! ammonium is distributed throughout the oxic zone)


			   ! Nitrified fraction of ammonium released in anoxic zone:

$if defined test_illegal_value           
            N_ANOX_NITR = (F_OX_B_REL/(F_OX_B_REL+0.5*DNITF))
            if ( N_ANOX_NITR .lt.0.0 ) then
               write (6,*)'N_ANOX_NITR =', N_ANOX_NITR,' <0' 
!               Pause
            endif
            N_ANOX_NITR = N_ANOX_NITR**DNITXP
$else
            N_ANOX_NITR = (F_OX_B_REL/(F_OX_B_REL+0.5*DNITF))**DNITXP
$endIf               
               ! Factor as above, but for ammonium generated
               ! in the anoxic zone. Thin oxic zones have
               ! double effect compared to above,
               ! because all ammonium passes the whole zone.


		   ! ------------------------------------------------------------------
            if (R_C_ANOX_MAX .le. 0.0) then
                 ! No potential for anoxic degradation
		   ! ------------------------------------------------------------------

               F_DENITR  = 0.0
               F_SULFRED = 0.0


		   ! ------------------------------------------------------------------
            else ! Distribute on denitrification and sulfate reduction:
		   ! ------------------------------------------------------------------


               IF ( OXYG_EQUIV .gt. SULFOX(1) ) THEN
                   THETA_S = 0.0
               ELSEIF ( SULFOX(2) .gt. 0.0 ) THEN 
                   ! increases to max value = 1. for oxyg--> - inf.
                   X = (SULFOX(1)-OXYG_EQUIV)**SULFXP
                   THETA_S = SULFRR* X/(SULFOX(2)**SULFXP + X )
               ELSE
                   THETA_S = 1.0
               ENDIF  ! Sulfate reduction inhibited by oxygen.


               X = R_C_ANOX_MAX / C_remin_tot
               U_S = THETA_S*X
               K_D = DNITRR * X


$IF DEBUG_LEVEL == 2
               if (DEBUG_LAYER) Write (DEBUG_UNIT,'(4(1X,A,1x,G12.6))')
     &            'THETA_S', THETA_S, 'U_S', U_S, 'X', X,  'K_D', K_D
$ENDIF


      ! -------------------------------------------------
               IF ( K_D .le. 0 ) then
                      ! No denitrification, any
                      ! anoxic degradation
                      ! happens by sulfate reduction:
      ! -------------------------------------------------

                   F_DENITR = 0.0
                   F_SULFRED = U_S
                      ! FOR DEBUG PURPOSES (DGRV), SEE BELOW:
                   N_max = 0.0
                   N_NO3 = 0.0
                   Q_D_OX  = 0.0


      ! -------------------------------------------------
               ELSE  ! Potential for denitrification:
                      ! including possible feedback
                      ! from nitrification:
      ! -------------------------------------------------


      ! Maximum nitrate consumption if all is denitrified:
              ! Stoichiometric N:C ratio corrected for NH3 -> N2:     
                   DENITR_C_CORR = DENITR_C + 
     &                     DENITR_NH3*(1-one_FDNH3)*N_REMIN/C_REMIN_TOT
                              ! one_FDNH3 = 1.0-fraction of NH3 reacting
                              ! with nitrate to molecular nitrogen
                              ! as part of denitrification
              ! Abs. value     
                   N_max = DNITRR * R_C_ANOX_MAX * DENITR_C_CORR
     

      ! External nitrate available for denitrification:
                   
!                   N_NO3 = N_MAX * NO3_EFF/ (DNITKS+NO3_EFF )
               
                    N_NO3 = N_MAX * NO3_EFF
     &                   /(DNITKS+F_OX_B_REL*KOXN+NO3_EFF)
                                  ! Oxic zone as barrier
                                  ! Depends on demand


      ! Nitrified ammonium from (bacterial) aerobic degradation
      ! available for denitrification:

                   Q_D_OX  = DNOXFR * N_OX_NITR * F_OXIC*N_Remin
               
               ! All released nitrogen included here.
               
               ! DNOXFR puts a limit to how large fraction can be
               ! denitrifed. (50% reasonable if particles have
               ! open surfaces, could presumably be higher behind
               ! some kind of membrane



$IF DEBUG_LEVEL == 2
      if (DEBUG_LAYER) THEN
            write (DEBUG_UNIT,'(3(1X,A12,G10.4))')
     &          'C_remin_tot', C_remin_tot, 'N_Remin', N_Remin,
     &          'R_C_ANOX_MAX', R_C_ANOX_MAX,
     &          'F_OX_B_REL', F_OX_B_REL,
     &          'F_OX_B_ABS', F_OX_B_ABS, 'F_Oxic', F_Oxic,
     &          'N_max',N_max, 'N_NO3',N_NO3,
     &          'N_OX_NITR', N_OX_NITR, 'Q_D_OX', Q_D_OX,
     &          'N_ANOX_NITR', N_ANOX_NITR
      endif
$ENDIF


                   H_D = ( DNOXFR + (1.0-DNOXFR)*N_OX_NITR)
     &                     * N_ANOX_NITR*N_remin
                     
                     ! Factor that multiplied with (F_DENITR+F_SULFRED)
                     ! will give amount of nitrified ammonium from
                     ! anoxic zone that is available for
                     ! denitrification (Up to 100 % will be
                     ! denitrified if N_OX_NITR = 1 )
                     ! At low values of N_OX_NITR: 
                       !   same fraction as Q_D_OX
                     ! At high values of N_OX_NITR: 100%



                   W_S = THETA_S / DNITRR
                   A_1 = N_MAX + Q_D_OX
                   A_2 = N_NO3 + Q_D_OX

$IF DEBUG_LEVEL == 2
                if (DEBUG_LAYER) THEN
                     write (DEBUG_UNIT,'(3(1X,A5,G10.4))')
     &                   'N_NO3', N_NO3, 'Q_D_OX', Q_D_OX, 'H_D', H_D,
     &                   'W_S', W_S, 'A_1', A_1, 'A_2', A_2
                endif
$ENDIF

                   IF ( H_D .LE. 0.0 ) THEN   ! Nitrate only from
                      F_DENITR = K_D*A_2/A_1  ! NO3 and aerobic degradation,
                      F_SULFRED = U_S - W_S*F_DENITR
                                              ! will always be >=0.

                   ELSE  ! Amonium from anoxic degr. may be nitrified
                         ! and denitrified: feedback and interaction.

                         ! First, attempt solution F_DENITR>=0, F_SULFRED>=0:
                      BETA_2 = H_D*K_D ! (also used below)
                      BETA_1 = (one_FDNH3 - W_S)*BETA_2
                      B = H_D*U_S+BETA_1
                      E = N_MAX - N_NO3
                      A_B = A_1 - B
                      AB_BetaE = A_1*B - BETA_1*E
                      X = A_B*A_B + 4*AB_BetaE 

$IF DEBUG_LEVEL == 2
                if (DEBUG_LAYER) THEN
                     write (DEBUG_UNIT,'(5(1X,A,''='',G10.4))')
     &                     'BETA_1', BETA_1, 'BETA_2', BETA_2,
     &                     'B', B, 'E', E, 'X=Rootsquare', X
                endif
$ENDIF
      
                      if ( X .ge. 0.0 ) then
                         
                         X = 0.5*( SQRT(X) - A_B )
                         
                         Div = X + A_B
                         if (Div .ne. 0.0D0) then !Iterative correction:
                             X = AB_BetaE / Div
                         Endif
                         
                         if ( X .ge. 0.0 ) THEN
                              ! X = H_D*(one_FDNH3*F_DENITR+F_SULFRED)
                            F_DENITR = K_D*(A_2+X)/(A_1+X)
                            X = X/H_D - one_FDNH3*F_DENITR
                         ENDIF !(Was only skipped if X.lt. 0.0)
                      endif

                      
                      if (X .ge. 0.0) then ! (F_DENITR was set above)
                         F_SULFRED = X
                      
                      ELSE ! Enough nitrate to denitrify all
                           ! available carbon: no sulfate reduction.
                           ! Solve 2. order equation for F_DENITR
                         A_Beta = A_1 - Beta_1
                         BetaA_BetaE = Beta_1*(A_1 - E)
                      
                         X = A_Beta*A_Beta + 4*BetaA_BetaE 
                         
                         X = 0.5* (SQRT(X)- A_Beta )
                         Div = X + A_Beta
                         if (Div .ne. 0.0D0) then !Iterative correction:
                             X = BetaA_BetaE / Div
                         Endif
                         
                         if (one_FDNH3.gt.0.0) then
                            F_DENITR = X/H_D/one_FDNH3
                         else
                            F_DENITR = K_D*A_2/A_1
                         Endif
                         F_SULFRED = 0.0
                      ENDIF
                   ENDIF
      

       ! ============= Controls the two equations: ================
$IF DEBUG_LEVEL == 2

                   Error_string=' '
                   do err_eq=1, 2
                      SELECT CASE (ERR_EQ)
                      CASE (1)
                         X = H_D*(one_FDNH3*F_DENITR+F_SULFRED)
                         if ( N_MAX + Q_D_OX + X .gt. 0.0 ) THEN
                            X = K_D*(N_NO3+Q_D_OX+X)/(N_MAX+Q_D_OX+X)
                         endif ! (else X=0 already)
                         A_1 = F_DENITR
                      CASE (2)
                         X = MAX ( 0.0D0, U_S -W_S*F_DENITR )
                         A_1 = F_SULFRED
                      END SELECT

                      IF ( ABS(X-A_1) .GT. 1.0E-11 ) THEN
                         write( Error_string,
     &                    '('' DGRADE: Feil i lign.'',I1,2(A,G19.11))')
     &                     Err_Eq,' v.s=',X,' h.s.=',A_1
                         WRITE ( Debug_unit, *) Error_string
                      Endif
                   enddo
$ENDIF


      ! -------------------------------------------------
               ENDIF
      ! -------------------------------------------------

	   ! ---------------------------------------------------------------
            ENDIF
	   ! ---------------------------------------------------------------


C #####################################################################
C #####################################################################
       ! End of code to be tested by DGRADE_test

C Spesiell test-kode: NB! Lagres før omregning
            if ( L .eq. LDGRV) then
               DGRV(1,N_FRAC) = F_OX_B_REL  ! relative values:
               DGRV(2,N_FRAC) = F_OX_B_ABS  !  "
               DGRV(3,N_FRAC) = F_OXIC      !  "
               DGRV(4,N_FRAC) = F_DENITR    !  "
               DGRV(5,N_FRAC) = F_SULFRED   !  "
               DGRV(6,N_FRAC) = MAX_RATE_C  ! Specific rate
                 ! Converted to volume basis always
               DGRV(7,N_FRAC) = C_remin_tot *VOLUME_FACTOR
               DGRV(8,N_FRAC) = R_Remin     *VOLUME_FACTOR
               DGRV(9,N_FRAC) = R_C_anox_max*VOLUME_FACTOR
               DGRV(10,N_FRAC)= N_max       *VOLUME_FACTOR
               DGRV(11,N_FRAC)= N_No3       *VOLUME_FACTOR
               DGRV(12,N_FRAC)= Q_D_OX      *VOLUME_FACTOR
               DGRV(13,N_FRAC)= 
     &                (F_DENITR*one_FDNH3 + F_SULFRED)*N_remin
     &                 *N_ANOX_NITR*VOLUME_FACTOR
                 !  = nitrate from anoxic degradation
                 !    available for denitrification
               DGRV(14,N_FRAC)= N_OX_NITR   ! Fraction
               DGRV(15,N_FRAC)= N_ANOX_NITR ! Fraction
               DGRV(16,N_FRAC) =  OX_ASED_INHIB
               DGRV(17,N_FRAC) =  OXYG_EQUIV
            endif



C -------------- Adjust distribution factors to sum to 1, and set
C                F_SUM as a reduction factor. This is done to be
C                able to use the factors as fractions of total
C                degradation below when including effect on
C                derivatives of water concentrations.
            F_SUM = F_OXIC + F_DENITR + F_SULFRED
                ! = Reduction factor for actual total degradation



            if (F_SUM .gt. 0.0) then       !Distribution factors:
               F_OXIC    = F_OXIC/F_SUM       ! adjusted to sum to 1.0
               F_DENITR  = F_DENITR/F_SUM
               F_SULFRED = F_SULFRED/F_SUM
            endif

$IF DEBUG_LEVEL == 1
      if (DEBUG_LAYER) THEN
        WRITE(DEBUG_UNIT,'(2(1X,A12,G13.7:))')
     &     'F_sum   '     , F_sum,
     &     'F_OXIC   '    , F_OXIC,
     &     'F_DENITR '    , F_DENITR,
     &     'F_SULFRED'    , F_SULFRED
      ENDIF
$ENDIF

       ! ----------- Actual sum of remineralization after
       !             external limiting factors are applied:

            C_Remin_tot = C_Remin_tot*F_SUM
            C_Remin = C_Remin*F_SUM   
                    ! both C terms are used below and must be rescaled
            N_Remin = N_Remin*F_SUM
            P_Remin = P_Remin*F_SUM
            S_Remin = S_Remin*F_SUM
            R_remin = R_remin*F_SUM


$IF DEBUG_LEVEL == 1
      if (DEBUG_LAYER) THEN
        WRITE( DEBUG_UNIT,
     &  '('' Actual values''/2( 1X, A, '':'', G15.7:))')
     &     'C_remin' , C_remin,
     &     'C_remin_tot' , C_remin_tot,
     &     'N_remin' , N_remin,
     &     'P_remin' , P_remin,
     &     'S_remin' , S_remin,
     &     'R_remin' , R_remin
      endif
$ENDIF



            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               SELECT CASE (N_FRAC)
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               CASE (1) ! Detritus in water
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            ! Update time derivatives of accumulated detritus:
            ! (was initiated by MTRANS, called previously by EUTRO.SRC)
                  CDETDV(L) = CDETDV(L) - C_remin
                  NDETDV(L) = NDETDV(L) - N_remin
                  PDETDV(L) = PDETDV(L) - P_remin
                  SDETDV(L) = SDETDV(L) - S_remin
                  RDETDV(L) = RDETDV(L) - R_remin
            ! Store "degradability" for monitoring/control purposes
                  RRDET(L)  = RR

$IF DEBUG_LEVEL == 1
      if (DEBUG_LAYER) THEN
        WRITE( DEBUG_UNIT, '( 2(1X, A9, '':'', G15.7:))')
     &     'CDETDV(L)' , CDETDV(L),
     &     'NDETDV(L)' , NDETDV(L),
     &     'PDETDV(L)' , PDETDV(L),
     &     'SDETDV(L)' , SDETDV(L),
     &     'RETDV(L)'  , RDETDV(L)
      endif
$ENDIF

            OXYG_USED = OXYG_USED + C_Remin_tot*OX_C*F_Sulfred
               !     Any sulphate reduction affects oxygen
               !     directly.


            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               CASE (2) ! Bottom sediment
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               ! Total loss (exclusive of specified organic carbon
               ! decomposition) including permanent burial in sediment:

                  C_loss = C_remin + RB*CSED( L )
                  N_loss = N_remin + RB*NSED( L )
                  P_loss = P_remin + RB*PSED( L )
                  S_loss = S_remin + RB*SSED( L )
                          ! [mg/m2/d = (1/d)*mg/m2]

                  R_remin= R_remin + RB*RSED( L )


          ! Initiate time derivatives of accumulated sediment matter:
          ! (new sedimentation is added by PRPROD, called later).
                  CSEDDV(L) = - C_loss
                  NSEDDV(L) = - N_loss
                  PSEDDV(L) = - P_loss
                  SSEDDV(L) = - S_loss
                  RSEDDV(L) = - R_remin ! R_Loss term not needed


            ! Store "degradability" for monitoring/control purposes
                  RRSED(L)  = RR


$IF DEBUG_LEVEL == 2
      if (DEBUG_LAYER) THEN
        WRITE(DEBUG_UNIT,*) ' Before residual correction'
        WRITE(DEBUG_UNIT,'(2(1X,A8,'':'',G15.8:))')
     &     'C_loss' , C_loss,  'C_remin' , C_remin,
     &     'N_loss' , N_loss,  'N_remin' , N_remin,
     &     'P_loss' , P_loss,  'P_remin' , P_remin,
     &     'S_loss' , S_loss,  'S_remin' , S_remin
      ENDIF
$ENDIF

      ! Include contribution from assumed residual of sediment
      ! present at start of simulation (per area).
      ! The contribution from assumed residual of sedimented matter
      ! is introduced as import to the system (summed over area)

               XC = XSED(L) ! Reduction factor for initial carbon
                            ! in sediment due to degradation

$if defined TEST_ILLEGAL_VALUE            
               IF (XC.LT.0.0) THEN
                  WRITE(6,*)' xc=',xc,' <0'
               ENDIF
$ENDIF
               XSEDDV(L) = -F_SUM*MAX_RATE_C*XSED(L)
                            ! Derivative of reduction factor:



               C_Remin = CORR_FAC( XB*XC) * C_remin
               C_Remin_tot =  C_remin + C_degr_spes*F_SUM
               CSEDXP(IB) = CSEDXP(IB) 
     &                     + (C_loss - C_remin_tot)*SED_AREA

               N_remin = CORR_FAC( XB*(XC**DGNFAC) ) * N_remin
               NSEDXP(IB) = NSEDXP(IB) + (N_loss - N_remin)*SED_AREA

               P_remin = CORR_FAC( XB*(XC**DGPFAC) ) * P_remin
               PSEDXP(IB) = PSEDXP(IB) + (P_loss - P_remin)*SED_AREA

               S_remin = CORR_FAC( XB*(XC**DGSFAC) ) * S_remin
               SSEDXP(IB) = SSEDXP(IB) + (S_loss - S_remin)*SED_AREA

$IF DEBUG_LEVEL == 2
      if (DEBUG_LAYER) THEN
        WRITE(DEBUG_UNIT,*) ' After residual correction'
        WRITE(DEBUG_UNIT,'(2(1X,A8,'':'',G15.8:))')
     &     'CSEDXP(IB)' , CSEDXP(IB),  'C_remin' , C_remin,
     &     'NSEDXP(IB)' , NSEDXP(IB),  'N_remin' , N_remin,
     &     'PSEDXP(IB)' , PSEDXP(IB),  'P_remin' , P_remin,
     &     'SSEDXP(IB)' , SSEDXP(IB),  'S_remin' , S_remin
      ENDIF
$ENDIF


        ! Increase in oxygen debt in sediment
        ! (only for degradation of bottom sediment):
            OXYG_DEBT_IN = C_Remin_tot*OX_C*F_Sulfred


	    ! Phosphorus degraded at bottom is stored in sediments
	    ! PADS(..) or released as function of sulphide in sediments
            if ( PADASD .gt. ASED_L ) then
               X = MAX ( 0.0D0, 1.0D0-ASED_L/PADASD )
            ELSEif ( ASED_L .gt. 0 ) then
               X = 0.0
            Else
               X = 1.0
            Endif

            X = PADMAX*X - MAX(0.0D0, PADS(L) )
                  ! first term = adjusted storage capacity as mg/m2

            if ( X .gt. 0) then ! Unused capacity:
                  ! Remineralized P retained in sediment (oxic cond.)
               P_net_retained = P_remin*PADRET* X/PADMAX
                                ! as fraction of mineralized P.
            else  ! X<=0: too much stored.
                  ! Excess in storage is released to water:
               P_net_retained = PADRLS* X ! mg/m2
            endif

	    ! P-buffer is buried with same rate as sulphide
            P_loss      = PADS(L) * RB * PSBURF  ! Factor relative to
                                                 ! general rate
C        Net change of buffer in sediment:
            PADSDV(L) = P_net_retained - P_loss

C        Net release to water:
            P_remin = P_remin - P_net_retained
            !    (no P is actively drawn from water to sediment,
            !     except by precipitation, see below )

C        Permanently exported:
            PSEDXP(IB) = PSEDXP(IB) + P_loss * Sed_AREA

$if defined DEBUG_ACTIVE
$if defined DEBUG_PHOSPHORUS
      if ( DEBUG_LAYER ) THEN
        WRITE(DEBUG_UNIT,'(3(1X,A,'':'',G13.7))')
     &     'OXYG_DEBT_IN', OXYG_DEBT_IN,
     &     'PADMAX ',PADMAX,
     &     'PADS(L)',PADS(L),
     &     'X=unused cap.',X,
     &     'P_remin'   ,  P_remin,
     &     'PO4_ADDED' , PO4_ADDED,
     &     'PO4_USED ' , PO4_USED,
     &     'P_precip'  , P_precip,
     &     'P_loss'     , P_loss,
     &     'PADSDV(L)',PADSDV(L),
     &     'PSEDXP(IB)' ,PSEDXP(IB)
      ENDIF
$ENDIF
$ENDIF

        ! Convert to pr. volume basis for computing effect on
        ! water concentrations (P already handled above):

               C_remin_tot = C_remin_tot*Area_on_Vol
               N_remin = N_remin*Area_on_Vol
               S_remin = S_remin*Area_on_Vol
               P_remin = P_remin*Area_on_Vol



        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           END SELECT
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C ----------------------------------------------------------------
C         Effect of degradation on water concentrations
C ----------------------------------------------------------------

	    ! ........... Nitrogen, N_remin is total amount degraded:

	    ! Removed by denitrification:
            NH3_removed = (1.0-one_FDNH3)*F_DENITR*N_remin   ! mg/m3/day
            NO3_removed = F_Denitr*Denitr_C_Corr*C_remin_tot ! mg/m3/day

	    ! Nitrified:
            N_Nitr = (   N_ANOX_NITR*(one_FDNH3*F_DENITR + F_SULFRED) 
     &                +  N_OX_NITR*F_OXIC  )*N_Remin
                         
                         ! assumes local nitrification to affect
                         ! all remineralized nitrogen, necessary to
                         ! reproduce observed low NH4 concentrations
                         ! in deep oxic waters
            
            NO3_ADDED = NO3_ADDED + N_Nitr
               ! (mg/m3/day)
            

            NH3_release = N_remin - N_Nitr - NH3_removed
            NH4_ADDED = NH4_ADDED + NH3_release
                        ! (mg/m3/day)


            N_removed = (NH3_removed+NO3_removed)*VL        ! mg/day  
            DENITR(IB) = DENITR(IB) + N_removed     ! sum total in basin
            
            NO3_USED =  NO3_USED + NO3_removed


            if ( L.eq.LDGRV ) then
               DGRV(18,N_FRAC) = N_remin ! Remineralized   N mg/m3/day
               DGRV(19,N_FRAC) = N_Nitr  ! Nitrified part     "
               DGRV(20,N_FRAC) = NH3_removed + NO3_removed      
                                            ! Denitrified part   "
            endif


C  .... Oxygen used by aerobic degradation and ammonium nitrification:
            OXYG_USED =   OXYG_USED + C_Remin_tot*OX_C*F_OXIC 
     &                  + N_NITR * OX_NITR
           !          (ml/liter/day)= (liter/m3/day)


C  .... Silicate:
            SiO2DV(L) = SiO2DV(L)  + S_remin

C  .... Phosphorus:
            PO4_added = PO4_added + P_remin


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ENDDO  ! Next fraction
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



C -------------------------------------------------------------------
C      Precipitation of phosphorus from oxygenated waters:
C -------------------------------------------------------------------

         ! Contribution from layer above (initiated above for top layer)
            PO4_ADDED = PO4_ADDED + P_Precip/VL

         ! Note P flux (for comparison with sediment traps):
            PDFLUX(L) = P_Precip/AREA(L)  ! MG/M2/DAY

         ! New equilibrium in precipitation from this layer:
            if (OXYG_L .le. 0.0 ) then
               X = 0.0
            else
               if ( OXYG_L .ge. PPOXMX ) then
                  X = 1.0
               else
                  X = (OXYG_L/PPOXMX)**PPOXEX
               endif
               X = min( DBLE(PPAMAX), PO4(L)*PPRMAX*X )
                    ! mg/m2/day, mg/m3*m/day * dim.less
               X = max ( 0.0D0, X )
            endif

            PADSDV(L) = PADSDV(L) + X  ! Settling on bottom (mg/m2/day)

            P_Precip = X*Next_Area ! mg/day sinking to next layer
                    ! stored to use for next layer

            PO4_USED  = PO4_USED + X*Area(L) /VL
                    ! Sinking out from this layer in total.


C --------------------------------------------------------
C     Balance oxygen debt in sediment
C --------------------------------------------------------

            X =  ASED(L) * RB 
            ASEDXP(IB) = ASEDXP(IB) + X*SED_AREA
            ASEDDV(L)  = OXYG_DEBT_IN - OXYG_DEBT_OUT - X
                ! literO2/m2/day = (mgC/m2/day)*(literO2/mgC)



$if defined DEBUG_ACTIVE
$if defined DEBUG_ASED
      if (DEBUG_LAYER) THEN
        WRITE ( DEBUG_UNIT, '(1X,2(A14,'':'',G15.7:))')
     &     'ASEDDV(L)', ASEDDV(L),
     &     'ASEDXP(IB)', ASEDXP(IB),
     &     'X=burial ', X,
     &     'OXYG_DEBT_IN',OXYG_DEBT_IN,
     &     'OXYG_DEBT_OUT',OXYG_DEBT_out,
     &     'OXYG_USED',OXYG_USED
      ENDIF
$ENDIF                                   
$ENDIF                                   


C ----------------------------------------------------------------
C   Update total derivatives of oxygen and free nutrients:
C ----------------------------------------------------------------

            OXYGDV(L) = OXYGDV(L) - OXYG_USED
            NH4DV(L)  = NH4DV(L)  + NH4_ADDED - NH4_USED
            NO3DV(L)  = NO3DV(L)  + NO3_ADDED - NO3_USED
            PO4DV(L)  = PO4DV(L)  + PO4_ADDED - PO4_USED

$IF DEBUG_LEVEL == 1
      if (DEBUG_LAYER) THEN
        WRITE ( DEBUG_UNIT, '(1X,2(A14,'':'',G15.7:))')
     &     'N_nitr', N_nitr,
     &     'NH3_release', NH3_release,
     &     'NH3_removed', NH3_removed,
     &     'NO3_removed', NO3_removed,
     &     'N_removed', N_removed,
     &     'NH4_ADDED', NH4_ADDED,
     &     'NH4_USED',  NH4_USED,
     &     'NO3_ADDED', NO3_ADDED,
     &     'NO3_USED',  NO3_USED,
     &     'DENITR(IB)',DENITR(IB)
        WRITE(DEBUG_UNIT,*)
     &      '     effect of biomass degradation:'
        WRITE( DEBUG_UNIT,'(1x, A12, 4A15)' )'after degr.',
     &          'CONS.','DERIV.', 'SED_AREA*DERIV'
        WRITE( DEBUG_UNIT, '(1X, A12, '':'', 3G15.7))' )
     &      'CSED', CSED(L), CSEDDV(L), SED_AREA*CSEDDV(L),
     &      'NSED', NSED(L), NSEDDV(L), SED_AREA*NSEDDV(L),
     &      'PSED', PSED(L), PSEDDV(L), SED_AREA*PSEDDV(L),
     &      'SSED', SSED(L), SSEDDV(L), SED_AREA*SSEDDV(L),
     &      'RSED', RSED(L), RSEDDV(L), SED_AREA*RSEDDV(L),
     &      'PADS', PADS(L), PADSDV(L), SED_AREA*PADSDV(L),
     &      'ASED', ASED(L), ASEDDV(L), SED_AREA*ASEDDV(L),
     &      'XSED', XSED(L), XSEDDV(L)
      ENDIF
$ENDIF


C ---------------------------------------------------------------
C   Update max.rate for limiting time-step:
C ---------------------------------------------------------------

            NEW_RATE = NO3_USED/(0.1+NO3_L)
            IF ( NEW_RATE.GT.MAX_BIO_RATE ) THEN
                 Max_bio_RATE=NEW_RATE
$if defined TEST_STEP
                 step_lim_layer = L
                 STEP_LIM_EVENT = -2
$ENDIF
            ENDIF
            
            NEW_RATE = NH4_USED/(0.1+NH4_L)
            IF ( NEW_RATE.GT.MAX_BIO_RATE ) THEN
                 Max_bio_RATE=NEW_RATE
$if defined TEST_STEP
                 step_lim_layer = L
                 STEP_LIM_EVENT = -3
$ENDIF
            ENDIF
            
            NEW_RATE = PO4_USED/MAX(0.1D0,PO4(L))
            IF ( NEW_RATE.GT.MAX_BIO_RATE ) THEN
                 Max_bio_RATE=NEW_RATE
$if defined TEST_STEP
                 step_lim_layer = L
                 STEP_LIM_EVENT = -4
$ENDIF
            ENDIF


$IF DEBUG_LEVEL == 1
      if (DEBUG_LAYER) THEN
        WRITE( DEBUG_UNIT,'('' ----------------- output :'')' )
        WRITE(DEBUG_UNIT,'(1X,A8,1X,2A15)')
     &     'VARIABLE','DERIVATIVE','DERIV*VOLUME'
        WRITE(DEBUG_UNIT,'(1X,A8,'':'',2G15.7:))')
     &     'DOC ', DOCDV(L)  , VL*DOCDV(L),
     &     'BACT', bactDV(L) , VL*bactDV(L),
     &     'OXYG', OXYGDV(L) , VL*OXYGDV(L),
     &     'NO3' , NO3DV(L)  , VL*NO3DV(L),
     &     'NH4' , NH4DV(L)  , VL*NH4DV(L),
     &     'PO4' , PO4DV(L)  , VL*PO4DV(L),
     &     'SiO2', SiO2DV(L) , VL*SiO2DV(L)
        WRITE(DEBUG_UNIT,'(''   DENITR='',G15.7)') DENITR(IB)
      ENDIF
$ENDIF
C ========================================================
         END DO  ! next layer
C ========================================================


C ========================================================
      END DO  ! next basin
C ========================================================


      MAX_BIO_TSTEP = 0.5/ MAX_BIO_RATE
!       Kan endres til: (kan evt. innføre en egen variabel)
!       MXTBIO = 0.5/ MAX_BIO_RATE  ! Initial value used by PRPROD.FOR


      end Subroutine

      end Module