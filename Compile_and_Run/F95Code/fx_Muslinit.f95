C =====================================================================
C Eutrophication model for inner Oslo fjord
C Norwegian Institute for Water Research  (NIVA)
C Birger Bjerkeng
C File: MUSLINIT.FOR
C   contains initiating part of submodel for mussels (Mytilus Edilis)
C   The integrating part is found in file MUSLINTG.FOR
C =====================================================================

      Module fx_MusselInitialise
      use ModelParam_Mussels
      use ModelParam_InitState
      use ModelParam_RunControl
      use ModelVar_Mussels
      use ModelVar_HydroBioChem
      use ModelVar_Topography
      use fx_RunControl
      implicit none


C DIAGNOSTIC OUTPUT CONTROLLED BY MSTEST TURNED ON OR OFF:
$undefine DEBUG_PRINT
C           =1: Debug_print to DEBUG_UNIT if requested through MSTEST



      contains


C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C Subroutine call to initiate new simulation:
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE MSINIT (MBI, MSAGES, MUSL_DIM_layers,
     &                   M_NR, M_WT, M_WA, M_LGT )


      INTEGER MBI, MSAGES, MUSL_DIM_layers
      real*8 M_NR   ( 0:MSAGES-1, MUSL_DIM_layers, MBI )
      real*8 M_WT   ( 0:MSAGES-1, MUSL_DIM_layers, MBI )
      real*8 M_WA   ( 0:MSAGES-1, MUSL_DIM_layers, MBI )
      real*8 M_LGT  ( 0:MSAGES-1, MUSL_DIM_layers, MBI )
!      INCLUDE 'EUTRO.INC'


C =============================================================
C Internal model description (dimensions, states, derivatives):
!        include 'MUSL.INC'
C =============================================================

      real*8 Prob, P_Weight, P_Area, Number_Factor
      real*8 Sum_Weight, Sum_zero_class

      real*8 Length
      Integer IB, L, Age
      real*8 MS, Wa
      Integer L_Base, MAXLAYER
      real*8 Sum_area
      real*8 X

C  ...... Statement functions for relation between size measures:
      real*8 INDIVID_AREA, W_OF_LENGTH, Length_of_W

      Individ_area(Length) = Area_per_sqlength* Length*Length ! m2
      W_of_length(Length ) = W_1cm*(Length**3.0) ! Weight pr. individual
      Length_of_W (Wa)     = (Wa/W_1cm)**(1./3.)  ! Length


C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C EXECUTE initiation, setting up arrays and start values:
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

$if defined DEBUG_PRINT
      IF (MSTEST) THEN
C          CALL OPN(Debug_unit)
         WRITE( debug_unit, '('' ####### MSINIT.FOR'')')
       ENDIF
$endif

      Age_groups_external_info = MSAGES ! Used in Mussel_integrate
C      TRENGS DENNE N


C ----------------- Set up storage arrays: -----------------------

C .......  Find number of layers with mussels for each basin:
      L = 1
      Do While ( L .lt. ND )
          Last_Fraction = (MUSLDP-Depth(L))/
     &                (Depth(L+1)-Depth(L))
          if (Last_Fraction .gt. 1.0 ) then
             L = L + 1
          else
             Exit
          endif
      Enddo  ! L = last required layer for mussels

      Last_Fraction = min (1.0D0, Last_Fraction)

C ........ check storage arrays:

      if ( L .gt. MUSL_DIM_layers ) THEN
         WRITE(*,*)
     &     'Musl_Dim_layers=',MUSL_DIM_layers, 
     &     'cannot include mussels to MUSLDP=',MUSLDP
!         PAUSE
      endif

      if (      L .gt. MLayers .or. NBI .gt. MBasins
     &     .or. MSAGES-1 .gt. Max_Age_Group ) then
          WRITE(*, '(1x,A/2(10x,2(A,I6)))' )
     &     'EUTRO.SRC dim. requires increased dim. in MUSL.CMN:',
     &     'MLayers=',MLayers, ', requires ', L,
     &     'MBasins=',MBasins, ', requires ', NBI,
     &     'Max_Age_Group =', Max_Age_Group, ', requires ', MSAGES-1
          Stop
      endif

      Layers = Min ( L, MUSL_DIM_layers, MLayers)
      MSNLAY = Layers  ! Part of external model description
      Basins = Min ( NBI, MBasins)
      Age_Groups = Min( MSAGES-1,Max_Age_Group )
                        ! MSAGES Counted from 1

$if defined DEBUG_PRINT
         IF (MSTEST) THEN
C           CALL OPN(Debug_unit)
            WRITE( debug_unit,*) ' Layers, Basins, Age_Groups:',
     &                             Layers, Basins, Age_Groups
         END IF
$endif


C ------------------ Store bottom areas: ------------------------
      Sum_Area = 0.0
      DO IB = 1, BASINS
         L_BASE = INDXI(IB)
         MaxLayer = INDXI (IB+1) - L_BASE
         X = LAST_FRACTION
C  Available bottom:
         Do L = Layers,1,-1
            if ( L.GT. MaxLayer) then
               MUSLBT (L,IB) = 0.0 ! Prevents all settling
            else
               MUSLBT (L,IB) = X * Bottom(L+L_BASE)
               X = 1.0
               Sum_Area= Sum_Area + MUSLBT(L,IB)*MCOVER(IB)
            endif
         Enddo
      Enddo

C -------------- fill arrays with initial values: -------------------

      Prob     = 1.0
      P_Weight = 0.0
      P_Area   = 0.0

      Do Age = 0, Age_Groups

C     ......... length, area and weight for individual of given age:

         if ( Age.eq.0 ) then  ! (No age class 0 initially)

            Length = 0.0
            Wa = 0.0
            X_Age (Age) = 0.0

         else

            Length = MSL_MAX*(1.0-exp(-MSL_K*Age)) ! cm
            Wa = W_of_length(Length) ! Weight pr. individual
            X_Age (Age) = Prob  ! Survival probability

         endif

      ! Sum probable mean area and weight occupied for each
      ! individual in class 0, assuming equilibrium population,
      ! and given the probability to survive until Age.

         P_weight = P_Weight + Wa*Prob
         P_Area   = P_Area   + Individ_area(Wa)*Prob


C           Store individual weights and lengths:

         Do IB = 1, Basins
            Do L = 1, Layers
               M_WA ( Age, L, IB)  = Wa  ! Pr. individual
               M_LGT ( Age, L, IB) = Length
            Enddo
         Enddo

         X = 1.0 ! Mortality
         Prob = Prob * exp( - X )


$if defined DEBUG_PRINT
         IF ( MSTEST ) THEN
            WRITE ( DEBUG_UNIT, '(1x, A,I5,3(1x,A,'':'',G10.3:) )' )
     &         'Age', Age, 'WA', WA, 'Length', length, 'Prob', Prob
            END IF
$endif


      Enddo


C     ....... Find number in each basin from area considerations, using
C             area and weight summed over age classes from above:

C     ....... Reduce to keep within specified total biomass in basin,
C             and convert weights to sum over individuals:
C             all weights as Carbon:
      if (Sum_Area .gt. 0.0 .and. .not. BIOOFF ) then
        X = Min(1.0D0, CMUSIN/(Sum_Area/P_area*P_weight))
                ! Correction factor used below,
                ! Upper limit to avoid unreasonable densities
                ! even if CMUSIN is unreasonably specified 
      else
        X = 0.0 ! (No mussel population possible, will be zeroed
                !  by MUSLBT=0, X value does not matter)
      endif

$if defined DEBUG_PRINT
       IF (MSTEST) THEN
          WRITE( debug_unit, '('' '',3(A,'':'',G15.8))')  
     &             'Sum_Area', Sum_Area, 
     &             'P_Weight',P_Weight,
     &             'X',X
       END IF
$endif


      Do IB = 1, Basins

$if defined DEBUG_PRINT
       IF (MSTEST) THEN
          WRITE( debug_unit, '('' ....... Basin IB = '',I3)') IB
       END IF
$endif


         Sum_Weight = 0.0
         Sum_zero_class = 0.0

         Do L = 1, Layers
            Number_Factor = X * MUSLBT (L, IB )*MCOVER(IB)/P_Area
                                
$if defined DEBUG_PRINT
            IF(MSTEST) THEN
               WRITE(DEBUG_UNIT,'('' ----- Layer L='',I5, 2(A,G15.7:))')
     &               L, 'MUSLBT=', MUSLBT(L,IB),
     &                  'Number_factor', Number_factor
               WRITE(DEBUG_UNIT,'(1x, A5,3A14)')
     &                  'Age', 'MS=M_NR', 'M_WA=M_WT', 'Wa/MS'
            END IF
$endif


            Do Age = Age_Groups, 0 , -1
               MS = X_age (Age) * Number_Factor
               Wa = MS * M_WA( Age, L, IB)
                    ! Individual ->total weight
               Sum_Weight = Sum_Weight + Wa
               M_WA ( Age, L, IB) = Wa  ! Sum over age class
               M_WT ( Age, L, IB) = Wa  ! No initial spawning material
               M_NR( Age, L, IB) = Ms


$if defined DEBUG_PRINT
               IF(MSTEST)THEN
                  if ( MS.gt.0 ) THEN
                     WRITE(DEBUG_UNIT,'(1x, I5, 3G14.7)')
     &                   Age, MS, Wa, Wa/MS
                  ELSE
                     WRITE(DEBUG_UNIT,'(1x, I5, 2G14.7)')
     &                   Age, MS, Wa
                  ENDIF
               END IF
$endif


            Enddo

            Sum_zero_class = Sum_zero_class + MS*M_WA(0,L,IB)


C Initial value of max. attainable weight first year:
            MUSLWM ( L, IB) = Wa_Settling

$if defined DEBUG_PRINT
            IF ( MSTEST ) THEN
                WRITE(DEBUG_UNIT,'(1x, A10,G14.7)')
     &                       'MUSLWM(L,IB):', MUSLWM ( L,IB )
            END IF
$endif

         Enddo

C     Corrected sum over age classes and layers, returned as carbon:

         CTMUSL(IB) = Sum_Weight
         CAMUSL(IB) = Sum_Weight
         C0MUSL(IB) = Sum_zero_class


$if defined DEBUG_PRINT
         IF ( MSTEST ) THEN
             WRITE(DEBUG_UNIT,'(1x, A10,G14.7)')
     &                       'CTMUSL(IB):', CTMUSL( IB ),
     &                       'CAMUSL(IB):', CAMUSL( IB ),
     &                       'C0MUSL(IB):', C0MUSL( IB )
!            if (debug_unit.eq.6) Pause
         END IF
$endif

      Enddo

C --------- Define area occupied by individual at settling size,
C           for later use in MUSLINTG:
      Settling_Length = Length_of_W (Wa_Settling)
      Settling_footprint = Individ_area(Settling_Length)

C ------------ Zero unused positions:
      Do IB = Basins+1, MBI
          CTMUSL(IB)=0.0
          CAMUSL(IB)=0.0
          C0MUSL(IB)=0.0
      Enddo

C ------------ Initiate control variables:
      INITIATED = .TRUE.
      DERIVATIVE = .FALSE.
      Spawning_start=1000.0
      Time_7_deg  = 1000.0

      end subroutine

      end Module