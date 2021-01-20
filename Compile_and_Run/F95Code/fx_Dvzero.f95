! ===================================================================
! NIVA Fjord model, File: fx_DVZERO.F95,  Birger Bjerkeng, NIVA.
! fetched from former ACSL code for main program of model
! ===================================================================

      Module fx_DVZero
      
      use ModelDimensions

      implicit none

      contains
      
      
C -------- Phase 1 of mass transport calculations:

      SUBROUTINE DVZERO( V_GROUPS, NLI, DERIV, IMPORT  )


      INTEGER V_GROUPS ! = Number of subgroups for this variable,
                       !   applies to water concentrations
                       !   and derivatives.
      INTEGER NLI


      real*8 DERIV(dimMLI, V_GROUPS), IMPORT (dimMBI)


C ------------------- local -----------------
      INTEGER I_GROUP,I_LAYER, I_BASIN

      DO I_GROUP = 1 , V_GROUPS
         DO I_LAYER = 1,NLI
            DERIV(I_LAYER,I_GROUP) = 0.0
         ENDDO
         DO I_BASIN = 1,dimMBI
            IMPORT(I_BASIN) = 0.0
         ENDDO
      ENDDO

      END Subroutine

      End Module
