! =============== stochiometric ratios: ===================
! called once when model run starts or continues 
! to initiate these stoichiometric ratios:
 
      Module fx_Stoichiometry
      use ModelParam_Decomposition, only: OXCFAC
                  ! Correction factor for oxygen vs. carbon
                  ! relative to "standard" composition CH2O

      implicit none

      real*4 OX_C, OX_NITR, DENITR_C, DENITR_NH3

! used in PHYT_ZOO.FOR, DGRADE.FOR, EutroSUB.FOR and Muslintg
      

      contains

      SUBROUTINE STOICH


      real*4 OX_C_v, OX_NITR_v, DENITR_C_v, DENITR_NH3_v

! aerobic degradation of carbohydrate:
       !  CH2O + O2  -> CO2 + H2O
       !     literO2/mgC = molO2/molC *literO2/molO2 /(mgC/molC)
      PARAMETER ( OX_C_v =     1.0    *32.0/1.429    /12000.   )

! Nitrification of ammonia:
       !  NH3 + 2O2 -> H2O + HNO3
       !         literO2/mgN = molO2/molN *literO2/molO2 /(mgN/molN)
      PARAMETER ( OX_NITR_v  = 2.0        *32.0/1.429    /14000.   )

! Degradation of CH2O by denitrification:
       !  5 CH2O + 4 HNO3 -> 7 H2O + 2 N2 + 5 CO2
       !               gN/gC = molN/molC  *gN/molN    /(gC/molC)
      PARAMETER ( DENITR_C_v = 4.0/5.0    *14.0       /12.        )

! Denitrification combined with oxidation of organic NH3 to N2:
       !  5 NH3 + 3 HNO3 -> 4 N2 + 9 H2O
       !               gN/gN:
      PARAMETER ( DENITR_NH3_v = 3.0/5.0 )

      OX_C       = OX_C_v*OXCFAC
      DENITR_C   = DENITR_C_v*OXCFAC
      OX_NITR    = OX_NITR_v
      DENITR_NH3 = DENITR_NH3_v

      END SUBROUTINE
      
      end Module

