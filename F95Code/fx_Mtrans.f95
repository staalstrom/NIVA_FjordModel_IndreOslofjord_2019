! ===================================================================
! NIVA Fjord model
! File: MTRANS.FoR
! Birger Bjerkeng, NIVA.
! fetched from former ACSL code for main program of model
! ===================================================================

      Module fx_MTRANS
      use ModelDimensions
      use ModelVar_HydroBioCHem
      use ModelVar_Topography
      use ModelVar_RunControl 
      use ModelParam_RunControl
      use ModelParam_Physics, only : HTRMIX
      use fx_Transp_1     

      implicit none
      
      contains



! Derivative and imports due to transports,
! transfer subroutine MTRANS to TRANSP_1:

      SUBROUTINE MTRANS( DERIV, IMPORT, V_GROUPS, V_CONC, EXTERNAL,
     &                   DEBUG_INDEX, STATE_NAME )


      INTEGER V_GROUPS ! = Number of subgroups for this variable,
C                 applies to water concentrations and derivatives.
C                 MASS_TANSPORT called once for each group,
C                 all groups collects into the same import terms,
C                 (IMPORT zeroed in MASS_TRANSP only for GROUP = 1)

!      INCLUDE 'EUTRODIM.inc' !D

C  IMPORTS ACSL COMMON BLOCK, for use in call to MASS_TRANSPORT:

!      INCLUDE 'EUTRO.INC'  !  ACSL Common block



      real*8 DERIV(dimMLI, V_GROUPS), IMPORT (dimMBI)
 
      real*8 V_CONC(dimMLI,V_GROUPS), EXTERNAL(dimMLE,V_GROUPS)
C           internal and external concentrations

      !! LOGICAL ADJUST_STATUS not needed (controlled ACSL code sorting)
C                whether negative net diffusion is allowed
C                to compensate numerical diffusion

      INTEGER DEBUG_INDEX
C                controls printing of debug information.
C
      CHARACTER*6 STATE_NAME  ! Hollerith string with variable name
C                 used in debug print within MASS_TRANSPORT


C ------------------- local -----------------
      LOGICAL M_DEBUG
      INTEGER I_GROUP


      M_DEBUG = (MDEBUG(DEBUG_INDEX)) ! .or. STATE_NAME(1:1).EQ.'P')
     &               .and. T.ge.TTRIG

      DO I_GROUP = 1 , V_GROUPS
         IF ( TRCALC ) THEN
C     ----------- real mass transport according to TRANSP subroutine:
            CALL  MASS_TRANSPORT( M_DEBUG,
     &                      I_GROUP,
     &                      V_CONC(1,I_GROUP), EXTERNAL(1,I_GROUP),
     &                      DERIV(1,I_GROUP),
     &                      IMPORT, STATE_NAME)
         ENDIF
      ENDDO

      END SUBROUTINE

      END Module
