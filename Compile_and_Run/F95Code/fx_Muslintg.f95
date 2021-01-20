C======================================================================
C Eutrophication model for inner Oslo fjord
C Norwegian Institute for Water Research  (NIVA)
C Birger Bjerkeng
C File: MUSLINTG.FOR
C   contains integrating part of submodel for mussels (Mytilus Edilis)
C   The initiating part is found in file MUSLINIT.FOR
C =====================================================================

      Module fx_MusselIntegrate
      use ModelParam_Mussels
      use ModelParam_Plankton
      use ModelParam_RunControl
      use ModelParam_InitState
      use ModelParam_Plankton
      use ModelParam_Decomposition
      use ModelVar_RunControl
      use ModelVar_Topography
      use ModelVar_Mussels
      use ModelVar_HydroBioChem
      use fx_RunControl
      use fx_Stoichiometry
      implicit none

$undefine DEBUG_DERIV
$if defined DEBUG_DERIV
$define DEBUG_DERIV_GROUPS
$endif

$undefine DEBUG_INTEGRATE
$if defined DEBUG_INTEGRATE
$define DEBUG_INTEGRATE_GROUPS
$endif

$define MUSLDV_VERSION 2
         ! =1: Subroutine arguments transferred in subroutine calls
         ! =2: Arguments defined in included Common block in MUSL.INC

      contains

C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C Subroutine call to calculate derivatives
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

$if MUSLDV_VERSION == 1
      SUBROUTINE MUSLDV ( Tinteg, BIOACT, DEBUG_ON, MLI, NBI, INDXI,
     &      VLAYER, AREA, BOTTOM, TEMP,
     &      MFYTG, FYTGRP, CFYT, NFYT, PFYT, SFYT, CHL,
     &      BACT, CZOO, NCZOO, PCZOO, NCBACT, PCBACT,
     &      PO4DV,  NH4DV, OXYGDV, SiO2DV,
     &      BACTDV, CZOODV,
     &      NCMUSL, PCMUSL, TMSPWN, TMSETL,
     &      MSVC,   MSERMX, MSREXP, MSINDW, MSWR, MSQW, MSBW,
     &      MSEASS, MSCREQ, MRSP15, MTRESP,
     &      UFRIC3, MFWFAC, MFILTM,
     &      GRMFYT, GRMBCT, GRMZOO, MCFMIN, MUSLDR, DGRATE,
     &      MXDETR, MSAGES, MSLAYR,
     &      M_NR, M_WT, M_WA, M_LGT, MUSLWM, MUSLBT, MCOVER,
     &    CFYTDV, NFYTDV, PFYTDV, SFYTDV, CHLDV,
     &    CDETDV, NDETDV, PDETDV, SDETDV, RDETDV,
     &    CSEDDV, NSEDDV, PSEDDV, SSEDDV, RSEDDV,
     &    CSEDXP, NSEDXP, PSEDXP, SSEDXP )
      

      real*8 TInteg
      LOGICAL BIOACT, DEBUG_ON
      INTEGER NBI, INDXI(NBI+1)
      INTEGER MLI
      real*8 BOTTOM (MLI)   ! Total bottom area within depth layer
      real*8 NCMUSL, PCMUSL ! Stochiometric ratio of mussels
      real*8 TMSPWN, TMSETL ! Time constants for spawning and settling

  ! Specific data on blue mussel individuals:

  ! Max. filtering capacity for individual of 1g dW:
      real*8 MSVC    ! l/h

  ! Upper limit to fraction of net growth used for reproduction.
      real*8 MSERMX

  ! Exponent of weight dependence for reproductive effort:
      real*8 MSREXP

  ! Weight relative to 1g dry weight soft tissue:
  !  (1): limits between low and high range for weight dependence
  !       of filtering and respiration
  !  (2): maximum effective weight in filtering relation
      real*8 MSINDW(2) ! dimensionless ratio
  
  !  (3): Target value of active weight in regard of
  !       switch to reproduction instead of tissue growth:
      real*8 MSWR

  ! Exponents of weight relation:
      real*8 MSQW (2)  ! on filtering
      real*8 MSBW (2)  ! on respiration

  ! Maximum ingestion efficiency for carbon, nitrogen and phosporus:
      real*8 MSEASS(3)

  ! Food concentration where unrestricted effective filtering ...
  ! equals physiological needs for ingested material
  ! (defines max. growth rate, could be set directly):
      real*8 MSCREQ    ! mgC/m3

  ! Routine respiration for individual of weight 1g Dw at 15 deg.C
      real*8 MRSP15    ! l O2/hour

! Temperature coeff. for exponent variation of respiration:
      real*8 MTRESP

      real*8 GRMFYT(MFYTG), GRMBCT, GRMZOO ! Grazing coefficients

      real*8 MCFMIN
      real*8 MUSLDR(4) ! 'Mussel mortality as rates pr. year 
                     ! (1): proportional to excess vs. area capacity
                     !      (=rate at 100% excess, only age >0)
                     !   (in settling class: included in TMSETL)
                     ! (2): intrinsic rate at settling
                     ! (3): intrinsic rate at 1cm length
                     ! (2) & (3) is combined into a term
                     !  decreasing exponentially with length
                     ! (4): rate increasing with high age (inversely
                     !      proportional to years left to MSAGMX)
                     ! Total rate is sum of terms 1 + (2&3) + 4
      real*8 MSAGMX    ! Age where rate (4) applies


      real*8 DGRATE(2) ! Degradation rates of fresh detritus or
                     ! sediment matter (1): from plankton (2): mussels

      real*8 MXDETR    ! Fraction of dead organic material that is
                     ! entered into the detritus phase,
                     ! the rest go directly into the sediment


      integer MSAGES      ! Number of age groups in external arrays
                          ! ( = Max age +1)
      integer MSLAYR      ! Max number of layers with Mussels

      real*8 M_NR   ( 0:MSAGES-1, MSLayr, NBI ) ! Number of mussels
      real*8 M_WT   ( 0:MSAGES-1, MSLayr, NBI ) ! ä total  weight mgC
      real*8 M_WA   ( 0:MSAGES-1, MSLayr, NBI ) ! ä active weight mgC
      real*8 M_LGT  ( 0:MSAGES-1, MSLayr, NBI ) ! Mean length (cm)
            ! (Addressed 1:MSAGES outside)

      real*8 MUSLWM  ( MSLayr, NBI )
           ! Max. individual weight of age 0. Only used as information
           ! to help interpret results and see if growth is reasonable,
           ! since mean weight is driven down by settling of larvae
           ! throughout growth season.
      real*8 MUSLBT  ( MSLayr, NBI ) ! Bottom area available
      real*8 MCOVER  (NBI)           ! Max. fraction used
      real*8 TEMP (MLI), VLAYER (MLI), AREA(MLI)
      real*8 UFRIC3(NBI) ! Friction velocity, measure of wind circulation,
                       ! low values may limit effective filtering
      real*8 MFWFAC      ! wind circulation effectivity factor,
C
      real*8 MFILTM(NBI) ! Volume exchange between water masses and
                       ! water close to mussels relative to filtered volume
      real*8 BACT (MLI)
      real*8 CZOO (MLI)
      Integer MFYTG, FYTGRP ! Array dimensions
      real*8 CFYT (MLI,MFYTG), NFYT(MLI, MFYTG), PFYT(MLI, MFYTG)
      real*8 SFYT(MLI)
      real*8 CHL  (MLI, MFYTG)
      real*8 NCZOO,  PCZOO
      real*8 NCBACT, PCBACT

C Water concentration derivatives, updated with effect of mussels:
      real*8 BACTDV (MLI)
      real*8 CFYTDV (MLI,MFYTG), NFYTDV(MLI,MFYTG), PFYTDV(MLI,MFYTG)
      real*8 SFYTDV (MLI)
      real*8 CDETDV (MLI), NDETDV(MLI), PDETDV(MLI), SDETDV(MLI)
      real*8 RDETDV (MLI)
      real*8 CSEDDV (MLI), NSEDDV(MLI), PSEDDV(MLI), SSEDDV(MLI)
      real*8 RSEDDV (MLI)
      real*8 CSEDXP (NBI), NSEDXP(NBI), PSEDXP(NBI), SSEDXP(NBI)
      real*8 CHLDV  (MLI, MFYTG )
      real*8 OXYGDV (MLI), NH4DV (MLI), PO4DV (MLI), SiO2DV(MLI)
      real*8 CZOODV (MLI)

$else

      SUBROUTINE MUSLDV  ( DEBUG_ON, MFYTG, MSAGES, MSLAYR, MBI,
     &                     M_NR, M_WT, M_WA, M_LGT )

      LOGICAL DEBUG_ON
      INTEGER MFYTG, MSAGES, MSLAYR, MBI
      real*8 M_NR   ( 0:MSAGES-1, MSLayr, MBI ) ! Number of mussels
      real*8 M_WT   ( 0:MSAGES-1, MSLayr, MBI ) ! Sum total weight mgC
      real*8 M_WA   ( 0:MSAGES-1, MSLayr, MBI ) ! Sum active weight mgC
      real*8 M_LGT  ( 0:MSAGES-1, MSLayr, MBI ) ! Mean length (cm)
        ! (NOTE: Addressed as 1:MSAGES outside of subroutine)


!      include 'eutro.inc' ! connects to model

$endif

C =============================================================
C ---- Internal model description (states, derivatives)
!      include 'MUSL.INC'

C ---- Internal Work variables:
!      INCLUDE 'MUSL_VAR.INC'

      Integer IB

C =============================================================
C Stochiometric ratios OX_C (+OX_NITR and DENITR_C, not used):
!      include 'STOICHIOM.inc'


$if DEBUG_DERIV
      Logical DEBUG_NOW, DEBUG_LAYER
      real*8 T_LAST
      save T_LAST
      real*8 Sum_W_T_Deriv, Sum_W_A_Deriv
$endif

$if defined DEBUG_DERIV_GROUPS
      logical Test_Group(0:Max_Age_Group),Test_Some_Group
      common /group_debug_control/
     &   Test_Group, Test_Some_Group

$endif

      real*8 Assim_respired

      real*8 sec_per_day
      parameter (sec_per_day=24.*3600.)


C  ...... Statement functions for relation between size measures:
      real*8 INDIVID_AREA, Length
      real*8 LENGTH_OF_W, Wa
      Individ_area(Length) = Area_per_sqlength*Length*Length ! m2
      ! W_of_length(Length ) = W_1cm*(Length**3.0) 
             ! Weight pr. individual, not used
      Length_of_W (Wa)     = (Wa/W_1cm)**(1./3.)  ! Length

$if defined DEBUG_DERIV
      DEBUG_NOW = DEBUG_ON .and. TEST_THIS_TIME (T, T_last)

      IF ( DEBUG_NOW ) THEN
           WRITE( debug_unit,
     &       '('' ####### MUSLDV at Tinteg = '',G14.7,
     &         '' Initiated='',L3,'' BIOACT='',L3 )' ) 
     &              Tinteg, Initiated, BIOACT
      END IF
$endif

      if ( .NOT. Initiated) RETURN

      INACTIVE = .NOT. BIOACT
          ! defined in musl.cmn, controls integrating call

      if (INACTIVE) THEN
         RETURN
      endif

      Time_in_year = MOD( TInteg, DAYS_PER_YEAR )
      Settling = Time_in_year .gt. Spawning_start + 30 
     &     .and. Time_in_year .lt. Spawning_start + 90
      

$IF DEFINED DEBUG_DERIV
      IF ( DEBUG_NOW ) THEN
         write(*,*)'Time_in_year, Settling:',Time_in_year, Settling

      endif
$endif

      Settling_Length = Length_of_W (Wa_Settling)
      Settling_footprint = Individ_area(Settling_Length)

      Layers = Min( MSLayr , MSNLay ) 
                   ! Dim.,   set during init,  
                   ! To avoid trouible with RESTOR command

$if defined DEBUG_DERIV
      IF ( DEBUG_NOW ) THEN
           WRITE( debug_unit,
     &      '(1x,A,L2,2(2(1x,A,G10.3:)/))' )
     &        'Settling='           , SETTLING,
     &        'Time_in_year='       , Time_in_year,
     &        'Spawning_start='     , Spawning_start,
     &        'Settling_Length='    , Settling_Length,
     &        'Settling_Footprint=' , Settling_Footprint
      END IF
$endif


      Assim_respired = min(1.0,max(0.0,MRASSF))
           ! fraction of gross tissue growth which is respired 


C ================================================================
C Calculate new derivatives for mussel stock,
C and update derivatives of affected concentrations
C of food, nutrients and oxygen
C=================================================================


C ################################################################
      Do IB = 1, Basins
C ################################################################

         L_BASE = INDXI(IB)

$if defined DEBUG_DERIV
       IF (DEBUG_NOW) THEN
          WRITE( debug_unit,
     &           '('' --------- Basin IB = '',I4, ''  L_Base='',I4 )')
     &                IB, L_Base
       END IF
$endif

         CSED_SUM = 0.0
         NSED_SUM = 0.0
         PSED_SUM = 0.0
         SSED_SUM = 0.0

C ################################################################
         Do L = 1, Layers
C ################################################################


            L_X = L_BASE + L  ! = global water layer index
            Tempr = TEMP   (L_X)  ! = temperature
            VL  = VLAYER (L_X)  ! = water volume

$if defined DEBUG_DERIV
       DEBUG_LAYER =  DEBUG_NOW .And. Test_this_layer(IB,L)
       IF (DEBUG_LAYER) THEN
          WRITE( DEBUG_UNIT,
     &          '('' ......... Layer '',2(2x,A,I3),2(2x,A,G12.5))')
     &           'L=',L,'L_X=',L_X,' Tempr',Tempr, 'VL', VL
       END IF
$endif

C     ............ Add up available food concentration:
C                  (All sources assumed equally available)
            BACT_L = MAX(0.0D0,BACT(L_X))*GRMBCT
            CZOO_L = MAX(0.0D0,CZOO(L_X))*GRMZOO
            C_FOOD = BACT_L        + CZOO_L
            N_FOOD = BACT_L*NCBACT + CZOO_L*NCZOO
            P_FOOD = BACT_L*PCBACT + CZOO_L*PCZOO
            DO FG = 1, MIN(MFYTG,FYTGRP)
               C_FOOD = C_FOOD + MAX(0.0D0,CFYT (L_X,FG))*GRMFYT(FG)
               N_FOOD = N_FOOD + MAX(0.0D0,NFYT (L_X,FG))*GRMFYT(FG)
               P_FOOD = P_FOOD + MAX(0.0D0,PFYT (L_X,FG))*GRMFYT(FG)
            ENDDO

C     ........... Growth efficient part of food concentration,
C                 including limiting assimilation efficiency:
            Assim_Conc = min( C_FOOD*MSEASS(1),
     &                        N_FOOD/NCMUSL*MSEASS(2),
     &                        P_FOOD/PCMUSL*MSEASS(3)  )


$if defined DEBUG_DERIV
            IF (DEBUG_LAYER) THEN
               WRITE( debug_unit, '(1x,4(A,'':'',G12.5))' )
     &           'C_FOOD', C_FOOD, 'N_FOOD', N_FOOD,
     &           'P_FOOD', P_FOOD, 'Assim_Conc', Assim_CONC
            END IF
$endif



C     ..................... Filtering rate .......................
C     ....... Max. filtering rate (m3/d) for individuals with
C             dry weight of soft tissue = 1g as function of Tempr:
C                   (valid from Food conc. from 0.03 to 10 mgDW/l?)
                v_c = 0.024  *  MSVC*exp(-0.3/max(0.01D0,Tempr+1.0) )
C               m3/d  m3/l*h/d   filtering in l/h
C                   (   no reduction during spawning or at high temp.
C                     - could be introduced )


$if defined DEBUG_DERIV_GROUPS
      Test_Some_Group =  Test_Groups
     &                  ( DEBUG_LAYER, Age_Groups, Test_Group, 0 )

$endif


C .... Find unrestricted volume filtering summed over population,
C      and sum occupied area:
            Sum_V = 0.0
            Sum_Area = 0.0
            Do Age = Age_Groups, 0, -1

               MS = M_NR (Age, L, IB) !Number of individuals

$if defined DEBUG_DERIV_GROUPS
                      ! trigger debug print for groups with mass <0
                      ! in addition to groups specified for debugging
               Test_Group(Age) =
     &             Test_Group(Age)
     &           .or. (MS .lt. 0.0)
     &           .or. (M_WA(Age, L, IB) .lt. 0.0)
     &           .or. (M_WA(Age, L, IB) .lt. 0.0)
               Test_Some_Group = Test_Some_Group .or. Test_Group(Age)
$endif

               if (MS .gt. 0.0) THEN
                  Wa = max( Wa_Settling, (M_WA(Age, L, IB))/MS )
                       ! = Active soft DW/individual as mg carbon
                       !   always at least settling weight

                  Wa_rel = Wa/W_1g ! Converts weight Wa as mg carbon 
                                   ! to total dry weight in gram for
                                   ! comparison with model parameters
                                   ! given as g dry weight.
                                   ! (W_1g = mgC per gram dry weight)

C           ........ Unrestricted filtering pr. individual,
C                    log_linear relation with breakpoint at MSINDW(1),
C                    upper limit at weight MSINDW(2):
                  if (Wa_rel.lt.MSINDW(1)) then
                     V = v_c * (Wa_rel/MSINDW(1))**MSQW(1)
                  else
                     V = v_c *( min(DBLE(MSINDW(2)),Wa_rel)
     &                          / MSINDW(1))**MSQW(2)
                  endif
                  X_Age  (Age) = V
                  Sum_V = Sum_V + V*MS !Total age-class max. filtering


                  X = MS* Individ_area( DBLE(M_LGT( Age, L, IB)) )
                  Sum_area = Sum_area + X
                      ! Accumulated sum will influence mortality



$if defined DEBUG_DERIV_GROUPS

       IF (Test_Group(Age) ) THEN
          WRITE( debug_unit,
     &          '((2x,A,I4),3(2x,A,G13.7:))')
     &           'Age =',Age, 'Wa_rel=',Wa_rel,'V=',V, 'äArea',X
       END IF
$endif

               Endif
            Enddo

C .........  General "overpopulation" death rate:            
            X = MUSLBT(L,IB)*MCOVER(IB)
            if ( X .gt. 0 ) THEN
               DEATH_RATE_AREA = MUSLDR(1) *
     &             Max (0.0D0, Sum_area/X-1.0 )
            Else
               DEATH_RATE_AREA = MUSLDR(1) ! (not important)
            endif


C     ....... Exchange of water volume to filter at given wind circulation
            V =  MFILTM(IB)*Sum_V
              ! Lower limit if no wind.
              ! (MFILTM = fraction of filtered volume pr. day)
            if (UFRIC3(IB).gt.0.0) then
               X =   MFWFAC*UFRIC3(IB)**(1./3.)/SQRT(AREA(L_X))
     &             * sec_per_day  ! velocity over travel distance
               V = MAX( V, VL*X)
            endif

$if defined DEBUG_DERIV
            IF (DEBUG_LAYER) THEN
               WRITE (DEBUG_UNIT,'( 3(1X,A,G10.3)/ 4(1X,A,G10.3) )')
     &               'FILTER-STAGE 1: MFILTM:',    MFILTM(IB),
     &               'MFWFAC:', MFWFAC, 'UFRIC3:', UFRIC3(IB),
     &               'AREA:', AREA(L_X), 'V:',V, 'Sum_V:',Sum_V,
     &               'äarea used: ', Sum_Area
            END IF
$endif

            if (Sum_V.gt.0) then
               X = V/(V+Sum_V)
            else
               X = 1.0
            endif
                ! = reduction of concentration at full filtration:
                !      V = flushing of water past mussels,
                !  Sum_V = filtering through mussels.
            
            if( L.eq.1) then 
               MVRED(IB) = X ! Export as info to user interface
                ! Can be used to adjust mfiltm to reasonable values
            endif


C     ....... Filtering is reduced if necessary to keep edible particle
C             concentration above cutoff value (assumed to be net result of
C             filtering being turned on/off around this value, with random
C             variation over population):
C                   This will reduce respiration (see below)
            if ( (C_FOOD-BACT_L)*X .lt. MCFMIN ) THEN
               if (C_FOOD-BACT_L .gt. MCFMIN .and. SUM_V .gt. 0.0 ) then
                   X = MIN(1.0D0, MCFMIN/(C_FOOD-BACT_L) ) ! ( >X above)
                   V_red = (V/Sum_V) *(1.0-X)/X ! < 1 because:
C                            [ V(1-X)/X<Sum_V   when   X>V/(V+Sum_V)  ]
               else
                   X = 0.0
                   V_red = 0.0
               endif
            else
               V_red = 1.0
            ENDIF

$if defined DEBUG_DERIV
            IF (DEBUG_LAYER) THEN
               WRITE (DEBUG_UNIT,'( 4(1X,A,G10.3))')
     &               '*STAGE 2: MCFMIN:', MCFMIN,
     &               'BACT_L:', BACT_L,  'X:', X, 'V_red:', V_red
            END IF
$endif


C            NOTE: ONLY PARTICULATE FOOD TRIGGERS ACTIVITY,
C                  ALTHOUGH EVEN BACTERIA ARE FILTERED TO SOME EXTENT

            Food_availability = X
C               :  Adjusted reduction in concentration close to mussels
C           (V_red: fraction of time used filtering)

C     ....... Basic (starving) respiration rate for 
C             mussel individs as mgC/day:
            R_c     = (MRSP15/ Ox_C)   * 24. * exp(MTRESP*(Tempr-15.))
C           mgC/day = [literO2/h /(literO2/mgC)*h/day


$if defined DEBUG_DERIV
            IF (DEBUG_LAYER) THEN
               WRITE (DEBUG_UNIT,'(4(1X,a,'':'',G10.3:))')
     &               'v_c', v_c, 'R_c(starving)', R_c,
     &               'Food_avail', Food_availability
                 Sum_W_A_DERIV= 0.0
                 Sum_W_T_DERIV= 0.0
            END IF
$endif


C     ............ activity of each age group:
            Sum_Spawning  = 0.0
            Sum_V = 0.0
            Sum_R = 0.0
            Sum_A = 0.0
            Sum_D = 0.0


         ! Loop down to age 0 to have values set for
         ! for special handling of age classe 0 below loop

         !-------------------------------------------
            Do Age = Age_Groups, 0, -1
         !-------------------------------------------

               MS = M_NR (Age, L, IB)

$if defined DEBUG_DERIV_GROUPS
                IF (Test_Group(Age)) THEN
                   WRITE (DEBUG_UNIT, '('' ... Age='',I2,
     &                                  3(1X,a,''='',G10.3))' )
     &               Age, 'MS', MS,
     &               'M_WA(...)', M_WA(Age, L, IB),
     &               'M_WT(...)', M_WT(Age, L, IB)
                END IF
$endif


               if (MS .le. 0.0) THEN

                  dWa_dt = 0.0
                  DWT_DT = 0.0
                  DR    = 0.0

               else

               ! ----------- filtering and growth --------------------

               ! ..... Mean soft tissue weight/individual as mg carbon.

                  WT = M_WT( Age, L, IB)/MS
                       ! = Total weight (in mass balance control)
                       !   Controls spawning 
                       
                  WR = Max( 0.0D0, Min(M_WA(Age,L,IB)/MS, WT ) )
                       ! = Effective respiration weight

                  WA = Max ( Wa_settling , WR )
                       ! = Growth active weight, 
                       !   controls filtering, growth and size



               ! The settling process (below) increase numbers,
               ! without increasing M_WT of Age class 0.

               ! Biomass increase M_WT is taken care of by the ordinary
               ! growth process below, to avoid separate code for
               ! biomass transfer connected to settling. This retards
               ! the growth a bit, but since the initial weight
               ! is insignificant, this is not important.

               ! However, active biomass M_WA for age class 0 is
               ! increased with WA_Settling pr. individual, if not,
               ! the mean active weight would be kept down, reducing
               ! growth of age classe 0 to substantially lower levels.
               ! As an extra safeguard, settling weight is a lower
               ! limit to the active weight.

               ! The addition to M_WA will ensure that age class 0
               ! will grow normally, at about the same absolute rate
               ! as if settling weight had been added into M_WT as well.

               ! Thus Wa will be > Wt for Age class 0. When production
               ! of spawning material starts Wt will get > Wa,
               ! denoting storage of spawning material.


                  Wa_rel = Wa/W_1g ! Converts weight Wa as mg carbon 
                                   ! to total dry weight in gram for
                                   ! comparison with model parameters
                                   ! given as g dry weight.
                                   ! (W_1g = mgC per gram dry weight)

C           .... Averaged activity rates pr. individual:
C              Volume effectively cleared (m3/d) at mean concentration
                  V = Food_Availability*V_red*X_Age(Age)
C

C              Food mass filtered (mgC/day):
                  F = V*C_FOOD

C              Asymptotic max. absorbtion rate (mgC/day) pr. individual:
C              proportional to unrestricted filtering:
                  A_Asymp = MSCREQ*WA_rel**MSCWXP
     &                     * MSEASS(1) * X_Age(Age)
C                           mgC/m3 * m3/day
                          ! X_Age stores unrestricted V pr. individual

$if defined DEBUG_DERIV_GROUPS

                  IF (Test_Group(Age) ) THEN
                      WRITE (DEBUG_UNIT, '(4(1x,A,''='',G10.3))' )
     &               'WA_rel', WA_rel, 'V' , V  , 'F'  , F,
     &               'A_Asymp'  , A_Asymp
                  END IF
$endif

C Assimilation:
                  A = min ( A_Asymp, V*Assim_CONC ) ! <= F*MSEASS
C                                   ( see def. of Assim_Conc)

C Respiration (mgC/d) controlled by accumulated weight W:
C Breakpoint around (W/W_1g)=(MSINDW(1)), with different
C slope of respiration/weight relation above and below.

                  X = (WR/W_1g)/MSINDW(1)
               
               ! Effective weight <=WT, WA
               ! This is to avoid unreasonably fast respiration
               ! and negative biomass for age class 0, but also to
               ! avoid increased respiration due to stored spawning
               ! material in higher age classes.
                  
                  if ( X .lt. 1.0 ) then
                     R = R_c * X**MSBW(1)
                  else
                     R = R_c * X**MSBW(2)
                  endif
                  R = max(0.0D0,(A - R)*Assim_respired) + R    


C              .... Net producton (growth) (mgC/d)
                  P = A - R

C              .... Sum effect over age groups, for calculating
C                   accumulated effect in water and sediment:
                  SUM_V = Sum_V + V*MS   ! Total volume filtered
                  SUM_R = Sum_R + R*MS   ! Total respiration as carbon
                  SUM_A = Sum_A + A*MS   ! Assimilation  (mgC/d)


C   ... Net growth (of Wt) pr. individual (mg/d):
                  
                  DWT_DT = P
                  dWA_dt = P

                  if ( Age.gt.0 ) then
                    if ( DWT_DT .gt. 0.0 ) then
                         ! Part of net positive growth is used
                         ! for producing spawning material:
                       X = DWT_DT *
     &                     Min ( 1.0D0, MSERMX* (Wa_rel/MSWR)**MSREXP)
                       dWa_dt = DWt_dt - X
                    endif

                     
                    if ( Tempr .gt. 7.0 ) then
                       if (Time_7_deg .gt. 999.0 ) then
                          ! note start time
                           Time_7_deg = Time_in_year
                       endif
                       
                       ! Activated after 15 days of T>7.0 degrees:
                       
                       if ( Time_in_Year .gt. Time_7_deg+15) then
                          if ( Spawning_start .gt. 999.0) then 
                               Spawning_start = Time_in_year
                          endif

                          ! Released from mussel biomass (mg/day):
                          Spawning = (Tempr-7)/(Tempr-5)
                          SPAWNING = max(0.0D0,(Wt-Wa)/TMSPWN*Spawning)  
                          DWT_DT = DWT_DT - Spawning
                          SUM_SPAWNING = Sum_Spawning + Spawning*MS
                       endif
                    else
                       Spawning = 0
                    endif
                  endif



C              .... Length is updated, increasing monotonously
C                   as function of max. active weight of age class
C                   that has occurred until now. 
C                   If length is increasing, death rate will have
C                   a lower limit given by individual size growth
C                   in each class:

                  LENGTH = Length_of_W(Wa) ! cm
                  if ( LENGTH .gt. M_LGT( Age, L, IB ) ) then
                      DR_LOW_LIM = 0.6777*dWa_dt/Wa ! (/day)
                  else
                      LENGTH = M_LGT( Age, L, IB )
                      DR_LOW_LIM = 0.0
                  endif    
                  if ( Length.lt.0.0 .or. 
     &                 abs(Length) .gt. 1.e10 ) Length = 0.0
                         ! Ensure against insufficient initiation
                  M_LGT( Age, L, IB ) = Length

                  ! Length is used to limit numbers versus area,
                  ! mean shell size can never be reduced,
                  ! except for age group 0 during settling
                  ! due to added smaller shells.
                  ! Effect of settling on mean length in age group 0
                  ! is taken care of in MUSSEL_INTEGRATE



      ! ---------- Mortality ------------------------
      ! Mortality rate is a sum of three rates,
      ! all given as 1/year:

      ! 1. given by area restrictions, computed above as
      !    Death_rate_Area, applies equally to all older age classes.
                  if (Age.gt.0) then
                     DR = Death_rate_Area
                  else
                     DR = 0
                  endif

      ! 2: MUSLDR(2) for yngel, MUSLDR(3) for eldre skjell:
                  if ( Age .eq. 0) then
                      DR = DR + MUSLDR(2)
                  else
                      DR = DR + MUSLDR(3)
                  Endif

      ! 3. A rate <increasing at high age, to give linear
      !    decrease in numbers for high ages:

                  X = MUSLDR(4)/Max( 1.0,MSAGMX-Age )
                  DR = (DR + X)/365. !  (/year)  --->  (/day)
                  DR = MAX ( DR_LOW_LIM, DR )
                  

$if defined DEBUG_DERIV_GROUPS

                IF (Test_Group(Age) ) THEN
                    WRITE (DEBUG_UNIT, '(5(1x,A,'':'',G10.3))' )
     &                'WA', WA, 'WT', WT, 'V'  , V  ,
     &                'F'  , F, 'A'  , A  , 'R' , R,  'P' , P,
     &               'SUM_V', SUM_V, 'SUM_R', SUM_R, 'SUM_A', SUM_A,
     &               'DR', DR, 'LENGTH', LENGTH,
     &               'DWT_DT', DWT_DT, 'dWa_dt', dWa_dt
                    WRITE (DEBUG_UNIT, '(2(1x,A,'':'',G10.3))' )
     &               'SUM_AREA', SUM_AREA,
     &               'SPAWNING', SPAWNING
!                   if (debug_unit.eq.6) Pause
                END IF
$endif
               ENDIF

C -------------------------------------------------
                  
               X = DR*MAX ( 0.0D0, M_WT(Age, L, IB) )
                 ! = biomass loss by death rates, ensures against 
                 !   "negative" death, if M_WT should become <0
                 !    by unforseen event

C            Add effect of mortality as total biomass died pr. day,
C            converted to sedimentary matter below, and
C            remineralized elsewhere in model.
               SUM_D = SUM_D + X

C            Store derivative information for use in MUSSEL_INTEGRATE:
               DEATH_RATE (Age,L,IB ) = DR
               
               W_T_DERIV (Age,L,IB) = dWT_DT*MS - X

               W_A_DERIV (Age,L,IB) = dWa_dt*MS - M_WA(Age,L,IB)*DR
                     ! not (MAX(0.0D0, ..... ) here,
                     ! M_WA is not included in mass balance, 
                     ! so there is no need to avoid drawing


$if defined DEBUG_DERIV_GROUPS
               IF (Test_Group(Age) ) THEN
                  WRITE (DEBUG_UNIT,'('' ... Result for age='',I5)') Age
                  WRITE (DEBUG_UNIT, '(4(1X,A,'':'',G10.3))' )
     &               'MS', MS,
     &               'W_A_DERIV', W_A_DERIV(Age, L, IB),
     &               'W_T_DERIV', W_T_DERIV(Age, L, IB),
     &               'DEATH_RATE', DEATH_RATE(Age,L,IB)
!                  IF(DEBUG_UNIT.eq.6) PAUSE
                END IF
$endif

$if defined DEBUG_DERIV
               IF (DEBUG_LAYER) THEN
                 Sum_W_A_DERIV= Sum_W_A_DERIV +  W_A_DERIV (Age, L, IB)
                 Sum_W_T_DERIV= Sum_W_T_DERIV +  W_T_DERIV (Age, L, IB)
                    ! Does not include settling addition to M_WA
               END IF
$endif




         !-------------------------------------------
            Enddo  ! Next age group.
         !-------------------------------------------


$if defined DEBUG_DERIV
            IF ( DEBUG_LAYER ) THEN
              WRITE ( DEBUG_UNIT,
     &        '('' ... Sum over age groups of:'',2(1x,A,''='',G10.3))')
     &               'W_A_DERIV', Sum_W_A_DERIV,
     &                'W_T_DERIV', Sum_W_T_DERIV
               END IF
$endif



! Accumulate max. individual weight for age 0:

	   !   Unrestricted filtering, and respiration:
            X =  MAX(0.0D0,MUSLWM(L,IB))/W_1g / MSINDW(1)
            if (X.lt. 0.0) then 
               write(6,*) 'Feil i linje 874 i MUSLINTG, X=',X
               stop
            endif

            if ( X .lt. 1.0 ) then
               V = v_c * X**MSQW(1)
               R = R_c * X**MSBW(1)
            else
               V = v_c *( min( DBLE(MSINDW(2)/MSINDW(1)), X)**MSQW(2))
               R = R_c *X**MSBW(2)
            endif
   
	   !   Net growth due to restrictions on filtering and assimilation:
            A = V*min( Assim_Conc*Food_Availability*V_red, 
     &                 MSCREQ*MSEASS(1)*(MUSLWM(L,IB)/W_1g)**MSCWXP )
            R = max(0.0D0,(A-R)*Assim_respired) + R
            W_max_DERIV (L, IB) = A - R


C     .......... Settling into age class 0, starts 30 days after
C                spawning started. Would start on empty area by
C                filling up available area with time constant TMSETL,
C                decreases with increasing occupied area, to 50%
C                at specified max. area.

            if ( SETTLING ) then
                X = MUSLBT(L,IB)*MCOVER(IB)
                if (X.gt.0.0) then
                   Settling_rate (L,IB) 
     &               = max( 0.0D0, X*X / (X + Sum_area)
     &                               / Settling_footprint / TMSETL )
                endif
            else
                Settling_rate (L,IB) = 0.0
            endif

$if defined DEBUG_DERIV
                IF (DEBUG_LAYER) THEN
                    WRITE (DEBUG_UNIT, '(3(1x,a,'':'',G10.3))' )
     &               'W_Max_DERIV:  '  , W_Max_DERIV(L,IB),
     &               'Settling_rate:'  , Settling_rate(L,IB),
     &               'SUM_SPAWNING', SUM_SPAWNING
!                   if (debug_unit.eq.6) Pause
                END IF
$endif


C     ......... Effect on water and sediment biological and chemical
C               water concentrations:

            Sum_area = Bottom(L_X)  ! Now used for total bottom area


$if defined DEBUG_DERIV
         if (DEBUG_LAYER ) THEN
            WRITE(DEBUG_UNIT,*)
            WRITE( DEBUG_UNIT,'('' before :'', 3A15)' )
     &         'CONC:', 'DERIV:','VL*DERIV:'
            DO FG = 1, MIN(MFYTG,FYTGRP)
               WRITE(DEBUG_UNIT,'(1X,A4,''(..,'',I3,'')'',3G15.7)')
     &        'CFYT',FG, CFYT(L_X,FG),CFYTDV(L_X,FG), VL*CFYTDV(L_X,FG),
     &        'PFYT',FG, PFYT(L_X,FG),PFYTDV(L_X,FG), VL*PFYTDV(L_X,FG),
     &        'NFYT',FG, NFYT(L_X,FG),NFYTDV(L_X,FG), VL*NFYTDV(L_X,FG),
     &        'CHL' ,FG, CHL (L_X,FG),CHLDV (L_X,FG), VL*CHLDV (L_X,FG)
            ENDDO
            WRITE(DEBUG_UNIT,'(5X,A8,3G15.7)')
     &        'BACT' , BACT(L_X) , BACTDV(L_X) , VL*BACTDV(L_X),
     &        'CZOO' , CZOO(L_X) , CZOODV(L_X) , VL*CZOODV(L_X),
     &        'SFYT' , SFYT(L_X) , SFYTDV(L_X) , VL*SFYTDV(L_X)
           WRITE(DEBUG_UNIT,'(5X,A8,15X,2G15.7:))')
     &         'NH4'  , NH4DV(L_X)  , VL*NH4DV(L_X) ,
     &         'PO4'  , PO4DV(L_X)  , VL*PO4DV(L_X) ,
     &         'OXYG' , OXYGDV(L_X) , VL*OXYGDV(L_X)
!            if (debug_unit.eq.6) Pause
            WRITE( DEBUG_UNIT,'('' in sediment:'', A10, 3A15)' )
     &         '..DV1','SUM_AREA*..DV1', '..DV2','SUM_AREA*..DV2'
            WRITE( DEBUG_UNIT,'(1X,A12,'':'',2G15.7:)')
     &         'CSED', CSEDDV(L_X), Sum_AREA*CSEDDV(L_X),
     &         'NSED', NSEDDV(L_X), Sum_AREA*NSEDDV(L_X),
     &         'PSED', PSEDDV(L_X), Sum_AREA*PSEDDV(L_X),
     &         'SSED', SSEDDV(L_X), Sum_AREA*SSEDDV(L_X),
     &         'RSED', RSEDDV(L_X), Sum_AREA*RSEDDV(L_X)
!            IF(DEBUG_UNIT.eq.6) PAUSE
         END IF
$endif 

C        ..... Effect on bacterial and phytoplankton biomass:
            V = Sum_V/VL  ! Clearing converted to relative value
            DO FG = 1, MIN(MFYTG,FYTGRP)
               VF = V*GRMFYT(FG) ! (grazing availability/selectivity)
               CFYTDV(L_X,FG) = CFYTDV(L_X,FG)
     &                          - MAX(0.0D0,CFYT(L_X,FG))*VF
               NFYTDV(L_X,FG) = NFYTDV(L_X,FG)
     &                          - MAX(0.0D0,NFYT(L_X,FG))*VF
               PFYTDV(L_X,FG) = PFYTDV(L_X,FG)
     &                          - MAX(0.0D0,PFYT(L_X,FG))*VF
               CHLDV (L_X,FG) = CHLDV (L_X,FG)
     &                          - MAX(0.0D0,CHL (L_X,FG))*VF
            ENDDO

            Filtered_Si = Sum_V*GRMFYT(1)*MAX(0.0D0,SFYT(L_X)) !(mgSi/d)
            SFYTDV(L_X) = SFYTDV(L_X) - Filtered_Si/VL
            BACTDV(L_X) = BACTDV(L_X) - BACT_L*V


C         ..... Excretion of faeces and mussel mortality
C               increases sediment concentrations (as mg/m2).
C               A fraction MXDETR of excreted faeces
C               go into pelagic detritus,
C               the rest enter bottom sedinment directly.

            Excr_C  = Sum_V*C_FOOD - Sum_A         ! mgC/d)
            Excr_N  = Sum_V*N_FOOD - Sum_A*NCMUSL  ! (mgN/d)
            Excr_P  = Sum_V*P_FOOD - Sum_A*PCMUSL  ! (mgP/d)
            Excr_Si = Filtered_Si                 ! (mgSi/d)


            if ( Sum_Area .gt. 0 ) THEN ! (if =0, Excr=0 also)

C           .... C, N and P in two degradable and one residual fraction:

               X = 1.0-MXDETR  ! = Fraction of faeces entering sediment
                               !   the rest goes to detritus

        !  Total amount to sediment (mg/day):

               C_TO_SED = Excr_C*X + Sum_D
               R_TO_SED = Excr_C*X*DGRATE(1) + Sum_D*DGRATE(2)
               N_TO_SED = Excr_N*X + Sum_D*NCMUSL
               P_TO_SED = Excr_P*X + Sum_D*PCMUSL

        !  Accumulated over layers:

               CSED_SUM = CSED_SUM + C_TO_SED
               NSED_SUM = NSED_SUM + N_TO_SED
               PSED_SUM = PSED_SUM + P_TO_SED
               SSED_SUM = SSED_SUM + Excr_Si*X
        !  Impact on bottom conc. of organic matter:
               XX = 1.0/Sum_Area
               CSEDDV(L_X) = CSEDDV(L_X) + C_TO_SED*XX
               RSEDDV(L_X) = RSEDDV(L_X) + R_TO_SED*XX
               NSEDDV(L_X) = NSEDDV(L_X) + N_TO_SED*XX
               PSEDDV(L_X) = PSEDDV(L_X) + P_TO_SED*XX
               SSEDDV(L_X) = SSEDDV(L_X) + Excr_Si*X*XX
        !  Residal feces to pelagic detritus:
               XX = MXDETR/VL
               CDETDV(L_X) = CDETDV(L_X) + Excr_C*XX
               RDETDV(L_X) = RDETDV(L_X) + Excr_C*XX*DGRATE(1)
               NDETDV(L_X) = NDETDV(L_X) + Excr_N*XX
               PDETDV(L_X) = PDETDV(L_X) + Excr_P*XX
               SDETDV(L_X) = SDETDV(L_X) + Excr_Si*XX
            END IF

C        ......... Zooplankton affected by both grazing and spawning:

            Sum_Spawning= Sum_Spawning / VL
C                  If higher nutrient requirements in zooplankton,
C                  effective spawning is reduced acordingly,
C                  and the excess carbon respired,
C                  to conserve N and P:
            Spawning = Sum_Spawning
     &                 * Min ( 1.0, NCMUSL/NCZOO, PCMUSL/PCZOO )
            X = Spawning - CZOO_L*V

$if defined DEBUG_DERIV
                IF (DEBUG_LAYER) THEN
                    WRITE (DEBUG_UNIT, '(3(1x,a,'':'',G10.3))' )
     &                'finally : SUM_SPAWNING',  SUM_SPAWNING,
     &                'W_Max_DERIV:  '  ,   W_Max_DERIV(L,IB),
     &                'SPAWNING',  SPAWNING,
     &                'CZOO_L', CZOO_L, 
     &                'V',V,'X',X,
     &                'X is added to CZOODV(L_X) = ',CZOODV(L_X)
!                   if (debug_unit.eq.6) Pause
                END IF
$endif

            CZOODV(L_X) = CZOODV(L_X) + X

C        ......... Respiration and superfluous carbon in spawning
C                  affects oxygen concentrations:
            Sum_R = Sum_R/VL  ! Effect as mgC/day/m3
            XX =   - ( Sum_R + Sum_Spawning - Spawning)*OX_C
            OXYGDV(L_X) = OXYGDV(L_X) + XX

C         ..... Respiration changes nutrient concentrations,
C               together with any excess nutrient in spawned material
C               over what is needed according to zooplankton ratios:
            XX = (Sum_R+Sum_Spawning)*NCMUSL-Spawning*NCZOO
            NH4DV (L_X) = NH4DV (L_X) + XX
            XX = (Sum_R+Sum_Spawning)*PCMUSL-Spawning*PCZOO
            PO4DV (L_X) = PO4DV (L_X) + XX
            SiO2DV (L_X) = SiO2DV (L_X)

$if defined DEBUG_DERIV
         if (DEBUG_LAYER ) THEN
            WRITE(DEBUG_UNIT,*)
            WRITE( DEBUG_UNIT,'('' after :'', 2(2A15:10X))' )
     &         ('DERIV:','VL*DERIV:',FG=1,2)
            DO FG = 1, MIN(MFYTG,FYTGRP)
               WRITE(DEBUG_UNIT,'(2(2X,A4,''(.,'',I2,'')'',2G14.7))')
     &        'CFYT',FG, CFYTDV(L_X,FG), VL*CFYTDV(L_X,FG),
     &        'PFYT',FG, PFYTDV(L_X,FG), VL*PFYTDV(L_X,FG),
     &        'NFYT',FG, NFYTDV(L_X,FG), VL*NFYTDV(L_X,FG),
     &        'CHL' ,FG, CHLDV (L_X,FG), VL*CHLDV (L_X,FG)
            ENDDO
            WRITE(DEBUG_UNIT,'(2(2x,A8,2G15.7))')
     &        'BACT' , BACTDV(L_X) , VL*BACTDV(L_X),
     &        'CZOO' , CZOODV(L_X) , VL*CZOODV(L_X),
     &        'SFYT' , SFYTDV(L_X) , VL*SFYTDV(L_X),
     &        'OXYG' , OXYGDV(L_X) , VL*OXYGDV(L_X),
     &        'NH4'  , NH4DV (L_X)  , VL*NH4DV(L_X),
     &        'PO4'  , PO4DV (L_X)  , VL*PO4DV(L_X)
            WRITE( DEBUG_UNIT,'('' Detritus:'', 2A15)' )
     &         '..DV1','Volume*..DV1'
            WRITE( DEBUG_UNIT,'(1X,A12,'':'',2G15.7:)')
     &         'CDET', CDETDV(L_X), VL*CDETDV(L_X),
     &         'NDET', NDETDV(L_X), VL*NDETDV(L_X),
     &         'PDET', PDETDV(L_X), VL*PDETDV(L_X),
     &         'SDET', SDETDV(L_X), VL*SDETDV(L_X),
     &         'RDET', RDETDV(L_X), VL*RDETDV(L_X)
            WRITE( DEBUG_UNIT,'('' in sediment:'', 2A15)' )
     &         '..DV1','Sum_AREA*..DV1'
            WRITE( DEBUG_UNIT,'(1X,A12,'':'',2G15.7:)')
     &         'CSED', CSEDDV(L_X), Sum_AREA*CSEDDV(L_X),
     &         'NSED', NSEDDV(L_X), Sum_AREA*NSEDDV(L_X),
     &         'PSED', PSEDDV(L_X), Sum_AREA*PSEDDV(L_X),
     &         'SSED', SSEDDV(L_X), Sum_AREA*SSEDDV(L_X),
     &         'RSED', RSEDDV(L_X), Sum_AREA*RSEDDV(L_X)
!            IF(DEBUG_UNIT.eq.6) PAUSE
         END IF
$endif



C ################################################################
         End do    ! Next layer
C ################################################################



C ################################################################
      End do       ! Next basin
C ################################################################


      DERIVATIVE =.true.

      END SUBROUTINE

C ================================================================
C Advance state values to new time T:

      SUBROUTINE MUSSEL_INTEGRATE ( MSTEST, T, T_Prev, MBI,
     &     MSAGES, MSLAYR,
     &     M_NR, M_WT, M_WA, M_LGT, MUSLWM,
     &     CTMUSL, CAMUSL, C0MUSL, CMUSDV,
     &     MUSLWA, MUSLWT )


      LOGICAL MSTEST
      real*8 T, T_prev
      INTEGER MBI

      integer MSAGES      ! Number of age groups in external arrays
                          ! ( = Max age +1)
      integer MSLAYR      ! Max number of layers with Mussels

      real*8 M_NR   ( 0:MSAGES-1, MSLayr, MBI )
      real*8 M_WT   ( 0:MSAGES-1, MSLayr, MBI )
      real*8 M_WA   ( 0:MSAGES-1, MSLayr, MBI )
      real*8 M_LGT   ( 0:MSAGES-1, MSLayr, MBI )
            ! (Addressed 1:MSAGES outside)
      real*8 MUSLWM               ( MSLayr, MBI )

      ! weight as mg Carbon:
      real*8 CTMUSL (MBI)  ! Sum total soft tissue weight allclasses
      real*8 C0MUSL (MBI)  ! Sum total weight of age class 0
      real*8 CAMUSL (MBI)  ! Sum active weight over all classes
      real*8 CMUSDV (MBI)  ! Derivative in time of CTMUSL (pr.day)

C For each class:
      real*8 MUSLWT (* ) ! Mean total weight of individuals
      real*8 MUSLWA (* ) ! Mean total weight of individuals


$if DEBUG_INTEGRATE
      Logical DEBUG_UNIT_OPEN
      Logical DEBUG_NOW, DEBUG_LAYER
      real*8 T_LAST_DEBUG
      SAVE T_LAST_DEBUG
      real*8 SUM_W_A_OLD, SUM_W_T_OLD
$endif


C =============================================================
C Internal model description (states, derivatives):

!        include 'MUSL.INC'


C =============================================================
C Stochiometric ratios OX_C (+OX_NITR and DENITR_C, not used):
!      include 'STOICHIOM.INC'
C Internal variables:


      real*8 Sum_Weight, Sum_active_weight, Sum_Zero_class
      real*8 Sum_Weight_Increase


      real*8 TSTEP
      Integer New_Year, K1,K2
      Integer IB, L, Age, NewAge, Age_G_Export
      real*8 MS, Wa, Wt, dWt_dt, MS_Prev, Wa_Prev, Wt_Prev
      real*8 X, Length_Prev, Length, X1, X2
      real*8 T_LAST
      SAVE T_LAST


$if defined DEBUG_INTEGRATE_GROUPS
      logical Test_Group(0:Max_Age_Group),Test_Some_Group
      common /group_debug_control/
     &   Test_Group, Test_Some_Group
$endif


      TSTEP = MAX(0.0D0, T - T_Prev )


$if defined DEBUG_INTEGRATE
      debug_now = mstest . AND. Test_this_time ( T, T_last_debug )


      IF (DEBUG_NOW) THEN
         WRITE (debug_unit, 
     &      '(//'' ####### MUSLINIT.FOR, MUSSEL_INTEGRATE'')')
         WRITE( debug_unit,'('' Initiated='',L8,'' Inactive='',L8 )')
     &              Initiated, Inactive
         WRITE( debug_unit,'(1x,3(A,G15.8))')
     &      ' T_Prev=', T_Prev,' T=', T,' TSTEP=', TSTEP
      ENDIF

      call HELLO (' has entered MUSSEL_integrate')

$endif

C -------------------------------------------------------------



      if (INACTIVE) then
         Do IB = 1, Basins
             Cmusdv(ib) = 0.0
         Enddo
         goto 300
      Endif

      if ( .not. DERIVATIVE ) GO TO 200

C ============================================================
C Update values with change over TSTEP.
C ============================================================

C If a new year has been entered, (New_Year=1) updated values for
C age class 0 to (Age_Groups - 2) are moved to next age class,
C (Age_Group - 1) is added to last class (Age_Group),
C and the new values for age class 0 are set to zero.
C For layers with MUSLBT = 0, no settling or growth
C will occur, because derivatives will be set to zero
C (no initial biomass, see MUSLINIT).

      if (T_last .ne. T_prev ) then
          write(*,*)' T_last=',T_last, '<> T_prev=', T_prev


$if defined DEBUG_INTEGRATE

          IF ( Debug_unit.gt.10 ) THEN
             inquire(DEBUG_UNIT, OPENED = DEBUG_UNIT_OPEN )
             IF ( DEBUG_UNIT_OPEN ) THEN
                  WRITE( debug_unit,*)
     &              ' T_last=',T_last, ' T_Last-T_prev=',T_Last-T_prev

             END IF
         END IF

$endif


          NEW_YEAR = 0 ! To avoid trouble
      Else

         X1 =  T / DAYS_PER_YEAR
         K1 = INT ( X1 )
         X2 =  T_Prev / DAYS_PER_YEAR
         K2 = INT ( X2 )


         NEW_YEAR = MIN(1, MAX( 0,  K1 - K2 ) )


$if defined DEBUG_INTEGRATE
         IF (DEBUG_NOW) THEN
!          write(*,*) 'X(T)     ',X1
!          write(*,*) 'X(T_Prev)',X2
           WRITE( debug_unit,
     &       '(1x,3A6/1x,3I6)' )
     &      'K1','K2', 'New_year',
     &       K1,  K2,   New_year
         END IF
$endif

      endif


      T_last = T


$if defined DEBUG_INTEGRATE
         IF (DEBUG_NOW) THEN
           WRITE( debug_unit,
     &      '(1x,3A14/1x,2G14.7,I10)' )
     &     'T', 'T_Prev', 'New_Year',
     &      T ,  T_Prev ,  New_Year
        END IF
$endif




      Do IB = 1, Basins

$if defined DEBUG_INTEGRATE
       IF (DEBUG_NOW) THEN
          WRITE( debug_unit, '('' .......UPDATE Basin IB = '',I3)') IB
          Sum_W_A_Old = 0.0
          Sum_W_T_Old = 0.0
       END IF
$endif

         Sum_Weight = 0.0
         Sum_Active_Weight = 0.0
         Sum_Weight_Increase = 0.0
         Sum_Zero_class = 0.0

         Do L = 1, Layers


$if defined DEBUG_INTEGRATE

       DEBUG_layer =  DEBUG_NOW .And. Test_this_layer(IB,L)
       IF ( DEBUG_LAYER ) THEN
           WRITE( debug_unit, '('' Layer L = '',I3)') L
       ENDIF
$endif



$if defined DEBUG_INTEGRATE_GROUPS
         Test_Some_Group =  Test_Groups ( DEBUG_LAYER,
     &                  Age_Groups, Test_Group, New_Year )
$endif


            Do Age = Age_Groups, -New_Year, -1
               NewAge = Age + New_Year ! always = 0 last time


               if ( Age.ge.0 ) then


$if defined DEBUG_INTEGRATE_GROUPS
         IF ( Test_Group(Age) ) THEN
             WRITE(DEBUG_UNIT, '(1x, A11, 4A14/ 13x, I3, 4G14.7)')
     &               'Before: Age', 'M_NR', 'M_WT', 'M_WA', 'M_LGT',
     &               Age, M_NR(Age,L,IB), M_WT (Age,L,IB),
     &                    M_WA(Age,L,IB), M_LGT(Age,L,IB)
         END IF
$endif


              ! ....... Updated values due to death and net growth:

                  MS =   M_NR(Age,L,IB)
     &                 * (1.0-DEATH_RATE(Age,L,IB)*TSTEP)

              !  Growth controlling weight:
                  Wa     = M_WA(Age,L,IB) + W_A_DERIV(Age,L,IB)*TSTEP

              !  Mass conservative weight;
                  dWt_dt = W_T_DERIV(Age,L,IB)
                  Wt     = M_WT(Age,L,IB) + dWt_dt*TSTEP
                  Length = M_LGT(Age,L,IB)


$if defined DEBUG_INTEGRATE_GROUPS
         IF ( Test_Group(Age) ) THEN
             WRITE(DEBUG_UNIT, '(1x, 2A14/ 1x, 2G14.7)')
     &               'W_T_Deriv','W_A_Deriv',
     &               W_T_Deriv(Age,L,IB), W_A_Deriv (Age,L,IB)
         END IF
$endif

              ! Sum values for export:

                  Sum_Weight = Sum_Weight + Wt
                  Sum_Active_Weight = Sum_Active_Weight + Wa
                  Sum_Weight_Increase = Sum_Weight_increase + dWt_dt

$if defined DEBUG_INTEGRATE
                  IF(MSTEST) THEN
                     Sum_W_A_Old = Sum_W_A_Old + M_WA(Age,L,IB)
                     Sum_W_T_Old = Sum_W_T_Old + M_WT(Age,L,IB)
                  END IF
$endif



               else  ! Happens only if New_Year>0:

                  MS = 0.0  ! Initiates values to be stored
                  Wa = 0.0  ! as new class 0.
                  Wt = 0.0
                  Length = 0.0
                  MUSLWM(L,IB) = Wa_settling
                       ! Reset maximum obtainable individ weight
                       ! for new age group 0.
               endif



C           .....Store into same or higher age groups, depending on
C                whether a new year was entered:

               if (NewAge .gt. Age_Groups) then

                  MS_Prev = Ms  ! Only if New_Year =1; save values to be
                  Wa_Prev = Wa  ! added to younger group in next pass
                  Wt_Prev = Wt
                  Length_Prev = Length

               else

$if defined DEBUG_INTEGRATE_GROUPS
               IF(Test_Group(NewAge)) THEN
                  IF ( MS .GT. 0.0 ) THEN
                     WRITE(DEBUG_UNIT,'(1x,a6,6a12/1x,I6,6G12.6)')
     &               'NewAge','MS','Length','WT','WA','WA/MS','WT/MS',
     &                NewAge, MS, Length, WT, WA, WA/MS, WT/MS
                  ELSE
                     WRITE(DEBUG_UNIT,'(1x,a6,4a12/1x,I6,4G12.6)')
     &                'NewAge', 'MS', 'Length', 'WT','WA',
     &                 NewAge, MS, Length, WT, WA
                  ENDIF
               END IF
$endif

                  if (New_Year.eq.1 .and. NewAge .eq. Age_Groups) then

                     MS = MS_Prev + Ms ! Add two last age groups
                     Wa = Wa_Prev + Wa
                     Wt = Wt_Prev + Wt
                     X = MS + MS_Prev
                     if (X.gt.0.0) then
                        X = (  MS*Length*Length 
     &                       + MS_Prev*Length_Prev*Length_Prev) / X
                        if (X.gt.0.0) then
                           Length = SQRT(X)
                        endif
                     endif

                  endif

                  M_NR (NewAge,L,IB) = MS
                  M_WT (NewAge,L,IB) = Wt
                  M_WA (NewAge,L,IB) = Wa
                  M_LGT (NewAge,L,IB) = Length


$if defined DEBUG_INTEGRATE_GROUPS
         IF ( Test_Group(NewAge) ) THEN
             WRITE(DEBUG_UNIT, '(1x, A11, 4A14/ 13x, I3, 4G14.7)')
     &               'After: NewAge', 'M_NR', 'M_WT', 'M_WA', 'M_LGT',
     &               NewAge, M_NR(NewAge,L,IB), M_WT (NewAge,L,IB),
     &                    M_WA(NewAge,L,IB), M_LGT(NewAge,L,IB)
         END IF
$endif



C NB! HVA MED SKALL-PRODUKSJON - ER OGS ENDEL AV ORGANISK PRODUKSJON,
C     ca. 30-70% av produksjon av soft tissue if›lge Dare 76.!
C     hhv. 8, 15 og 24 % av total for P, C og N
C     iflg. m†linger i Oslofjorden.


               endif




            End do  ! next age group

! Special handling og age group 0,
C presumes that age index above loops from maximum to age 0:

            Sum_Zero_class = Sum_Zero_class + Wa



      !  Add settling of larvae:

            X = SETTLING_RATE (L, IB)*TSTEP
            M_NR (0,L,IB) = MS + X  ! :increase in number

            if (MS + X.gt.0) then

            ! Adjust mean shell length for settling:
               M_LGT(0,L,IB) = (Length*MS + Settling_Length*X)/
     &                          M_NR (0,L,IB )

            ! Increase nominal weight controlling filtering and growth:
               M_WA(0,L,IB) =  M_WA(0,L,IB) + X*Wa_Settling

!        ...... Adjust max. mussel weight for age class 0
               if (MS.eq.0) then
                  MUSLWM(L,IB) = Wa_Settling
               else
                  MUSLWM(L,IB) = MUSLWM(L,IB) + W_max_deriv(L,IB)*TSTEP
               endif
            endif



$if defined DEBUG_INTEGRATE
$if defined DEBUG_INTEGRATE_GROUPS
            IF( Test_Group(0) ) THEN
$else
            IF( Debug_Layer ) THEN
$endif
               WRITE(DEBUG_UNIT,
     &            '('' Age 0: (Settling):'',
     &               2(3(1x,A,'':'',G10.3:)/))' )
     &                'MS',MS, 'X',X,
     &                'M_LGT',  M_LGT(0,L,IB),
     &                'M_NR',   M_NR (0,L,IB),
     &                'M_WA',   M_WA (0,L,IB),
     &                'MUSLWM', MUSLWM (L,IB)
            END IF
$endif


                  ! Reduction given by developing d(n*Lmean)/dt
            ! only numbers and mean length are updated here,
            ! biomass is not added here, this is instead taken
            ! care of by the growth process.



$if defined DEBUG_INTEGRATE
            IF(DEBUG_LAYER) THEN
               WRITE(DEBUG_UNIT,
     &            '('' M_NR(0,L,IB (settling):'',G14.7,
     &              '' MUSLWM(L,IB):'',G14.7 )' )
     &                     M_NR(0,L,IB), MUSLWM(L,IB)
            END IF
$endif

         End do ! next layer

         CTMUSL(IB) = Sum_Weight
         CAMUSL(IB) = Sum_Active_Weight
         C0MUSL(IB) = Sum_zero_class
         CMUSDV(IB) = Sum_weight_increase

$if defined DEBUG_INTEGRATE
            IF ( DEBUG_LAYER.and. TSTEP.ne. 0 ) THEN
                WRITE(DEBUG_UNIT, '(1x,2A21,A26)')
     &            'Old','New', 'Derivative calculated:'

                WRITE(DEBUG_UNIT, '(2(/1x,A4,'':'',2E21.12,G15.7))')
     &               'M_WT', Sum_W_T_Old, Sum_Weight,
     &                  (Sum_Weight-Sum_W_T_Old)/Tstep,
     &               'M_WA', Sum_W_A_Old, Sum_Active_Weight,
     &                  (Sum_Active_Weight-Sum_W_A_Old)/Tstep
               END IF
$endif

      End do

      if (New_year.gt.0) then 
          Spawning_start = 1000.0
          Time_7_deg = 1000.0
      endif


  200 CONTINUE


C ........... Sum over layers and basins, and calculate mean weights
C             within age_group:

$if defined DEBUG_INTEGRATE_GROUPS
            IF (Test_Some_Group ) THEN
               WRITE(DEBUG_UNIT, * ) ' Individ data by age class:'
               WRITE(DEBUG_UNIT,'(1x,'' UPDATED VALUES:''/1X,A5,5A14)')
     &                'Age', 'MS=MUSSEL', 'M_WA',
     &                'M_WT','Wa/MS','Wt/MS'
            END IF
$endif

      AGE_G_EXPORT = Min (Age_Groups, Age_groups_external_info-1 )
      Do Age = 0, Age_G_Export
         NewAge = Age + 1  ! Index 0:max -> 1:max+1 in exported arrays
         MS = 0.0
         WT = 0.0
         WA = 0.0
         Do IB = 1, Basins
            Do L = 1, Layers
                MS = MS + M_NR(Age,L,IB)
                Wa = WA + M_WA(Age,L,IB)
                Wt = Wt + M_WT(Age,L,IB)
            Enddo
         Enddo
         if ( MS .gt. 0.0 ) then
             MUSLWT (NewAge) = WT/MS
             MUSLWA (NewAge) = WA/MS
         else
             MUSLWT (NewAge) = 0.0
             MUSLWA (NewAge) = 0.0
         endif

$if defined DEBUG_INTEGRATE_GROUPS
         IF ( Test_Group(Age) ) THEN
             WRITE(DEBUG_UNIT, '(2x,I6,5G15.7)')
     &               Age, MS, WT, WA, MUSLWT(NewAge), MUSLWA(NewAge)
         END IF
$endif

      Enddo

 300  CONTINUE
      DERIVATIVE =.FALSE. ! NO LONGER VALID AT PRESENT TIME
      
      call HELLO (' leaves MUSSEL_integrate')

      end Subroutine


! =====================================================================
! Only defined if debug of derivative or integration is specified
! May restrict which points in time or water layers gets debug print:

$if DEBUG_DERIV || DEBUG_INTEGRATE  
      Logical function Test_this_time ( T, T_PREV )

          real*8 T
          real*8 T_Prev
          ! Test_this_time = ( NINT(T/7.0) .gt. NINT(T_PREV/7.0) )
          T_prev = T
          Test_this_time = .true.
       End Function


      Logical function Test_this_layer( IB, L )

          integer IB, L
          ! Test_This_layer= (ib.eq.1.and.L.eq.1)
          Test_This_layer = .true.
       End function

$endif


! =====================================================================
! Only defined if debug of age groups is specified.
! Controls which age-groups gets debug print:

$if defined DEFINED_DERIV_GROUPS || defined DEBUG_INTEGRATE_GROUPS

      Logical function Test_Groups ( MSTEST,
     &        Age_Groups, Test_Group, NSHIFT )
      

      Logical MSTEST
      Integer Age_Groups
      logical Test_Group (0:Age_Groups)
      integer NSHIFT

      integer Ages_Tested
      parameter ( Ages_Tested = 2 )
      integer debug_spec (Ages_Tested) /0,1/

      integer Index, AgeToDebug, AgeInvolved

      do Index = 0, Age_Groups
         Test_Group(Index) = .FALSE.
      Enddo
      Test_Groups = .false.

      if (MSTEST) then
         do Index = 1, Ages_Tested
            AgeToDebug =  debug_spec (Index)
            do AgeInvolved = AgeToDebug, AgeToDebug + NSHIFT
               if ( AgeInvolved .ge.0 .and.
     &              AgeInvolved .le. Age_Groups ) then
                  write(*,*)' Debug group:',AgeToDebug
                  Test_Group ( AgeToDebug ) = .true.
                  Test_Groups = .true.
               endif
            Enddo
         enddo
      Endif


      end function

$endif

      end Module
