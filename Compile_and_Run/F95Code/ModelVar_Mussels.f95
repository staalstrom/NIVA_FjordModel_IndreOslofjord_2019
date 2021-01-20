! =============================================================
! Mussel model variables
! Definition of commonblock MUSSEL_MODEL
! Accessed by both MUSLINIT & MUSLINTG
! Note! should be last included after all other definitions
! =============================================================

      Module ModelVar_Mussels

      use ModelDimensions
      implicit none

            ! Max. model size, dimensioning parameters:
      INTEGER Max_Age_Group

      Integer MBasins, MLayers
      PARAMETER ( Max_Age_Group = 9, MBasins = dimMBI, Mlayers = 11 )

            ! Actual size of model in current run:
      Integer Basins, Layers, Age_Groups

! Internal model description (states, derivatives):
!         For each age class: individual count and weight:



      real*8 Last_fraction

            ! Rates, derivatives:
      real*8 DEATH_RATE     ( 0:Max_Age_Group, MLayers, MBasins )
      real*8 W_T_DERIV      ( 0:Max_Age_Group, MLayers, MBasins )
      real*8 W_A_DERIV      ( 0:Max_Age_Group, MLayers, MBasins )

            !         For class 0:

      real*8 SETTLING_RATE  ( MLayers, MBasins )

            !         For class 0:

      real*8 W_MAX_DERIV    ( MLayers, MBasins )

            !         Work array over age groups:

      real*8 X_Age(0:Max_Age_Group)

            ! Model state control variables, with defined initial values:

      LOGICAL :: INITIATED = .false.
      LOGICAL :: DERIVATIVE = .FALSE.
      LOGICAL INACTIVE
      real*8 Spawning_start, Time_7_deg

            ! Fixed model parameters:

      real*8 :: DAYS_PER_YEAR = 365.0

            ! Growth characteristics:

      real*8 :: MSL_MAX = 6.0  ! Reported values from 6 to 10
      real*8 :: MSL_K = 0.2    ! Reported values from 0.2 to 0.8

            ! Settling characteristics:

      real*8 Settling_footprint, Settling_length

            ! Dimension of age group arrays in MUSSEL_INTEGRATE in MSLINTG.FOR:

      Integer Age_groups_external_info

				! Data on Mussel indivual basis:
				!  ...... Bottom area occupied related to square length:

      real*8, PARAMETER :: AREA_PER_SQLENGTH = 0.447*0.619/2/10000. ! (m2/cm2)

            !  ...... Reference weights as mg Carbon ( 45% of dry weight):

      real*8, parameter :: W_1g         = 450          ! mg Carbon/1g dry weight
      real*8, parameter ::  W_1cm       = 0.00592*W_1g ! mg Carbon at 1 cm length
      real*8, parameter ::  Wa_Settling = 1.26e-6*W_1g ! mg Carbon at settling
     
      
      real*8 TIME_IN_YEAR
      LOGICAL Settling
      real*8 Length
      Integer L, Age, IS
      real*8 MS, WR, Wa, Wt, Wm, Wa_rel
      Integer L_BASE, L_X
      real*8 Tempr
      real*8 BACT_L, CZOO_L, C_FOOD, N_FOOD, P_FOOD
      Integer FG
      real*8 VL, Assim_Conc
      real*8 Sum_Spawning, Sum_area, Sum_V, Sum_R, Sum_A, Sum_D
      real*8 Death_Rate_Area, DR, DR_LOW_LIM
      real*8 v_c, r_c
      real*8 V, F, A, R, P, X, XX , VF, A_Asymp
      real*8 Spawning, dWt_dt, dWa_dt
      real*8 Excr_C, Excr_N, Excr_P, Filtered_Si, Excr_Si
      real*8 C_to_Sed, N_to_Sed, P_to_Sed, R_TO_SED
      real*8 CSED_SUM, NSED_SUM, PSED_SUM, SSED_SUM
      real*8 Food_availability, V_red


      end Module
