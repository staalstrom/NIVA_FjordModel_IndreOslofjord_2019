      Module fx_TrCalc

      use fx_RunControl, only: DEBUG_UNIT, T
      use fx_Transp_1
      use fx_Transp_u
      
      implicit none

C Eutrofimodell Indre Oslofjord - Fil  TRCALC.FOR

$define DEBUG_MODE

C          Includes code for debug dump
C          Activated at run-time by argument DEBUG (one of parameter array MDEBUG(i))
C =====================================================================
C Fjord model
C   Transport calculations and mass-->concentration calculations
C               May 1990, Birger Bjerkeng, NIVA.
C
C ================================================================



      Contains

C Subroutine ADJ_CONC adjusts values for volume changes,
C to conserve mass within truncation errors.
C This corrects errors due to concentration changes due to transports
C being integrated as for unchanging volume.

      SUBROUTINE ADJ_CONC ( TRCALC, DEBUG,
     &        NAME, NBI, DVTOT, DVTOTL, VTOTZ, INDXI, NLI,
     &        VLAYER, NLVOPN, VFROPN,
     &        TSTEP, DERIV, CONC, CTRL_C )


C In:  NAME, NBI, DVTOT, DVTOTL, VTOTZ, INDXI, NLI
C In/Modified out: CONC

C ----------------------- External arguments: ------------------------
      LOGICAL TRCALC
C              - Indicates whether physical transports are active
      LOGICAL DEBUG
C              - Controls debug output
      CHARACTER*(*) NAME
C              - Name of substance, used in debug dumps
      INTEGER NBI
C              - actual number of internal basins
      real*8  DVTOT   (NBI)
C              - Deviation of current volumes from start volumes
      real*8  DVTOTL  (NBI)
C              - Deviation of previous volumes from start volumes
      real*8  VTOTZ  (NBI)
C              - Start volumes
      INTEGER INDXI(NBI+1)
C              - layer index limits for internal basins.
      INTEGER NLI
C              - number of layers in internal and external basins
      real*8 VLAYER(NLI)
C              - Volume of each layer
      INTEGER NLVOPN(NBI)
C              - Number of layers with open connection to boundary
      real*8    VFROPN(NBI)
C              - fraction of (normal) total volume in open layers

      real*8 TSTEP

      real*8 DERIV(NLI)
C              - time derivatives for given concentrations
      real*8 CONC(NLI)
C           - Concentrations in internal basins
C                   in:  preliminary stepwise integrated values
C                   out: corrected values.
      real*8 CTRL_C  (NLI)
           ! Controlling substance for applying neg. diffusion
           ! as correction of diffusive effect of advection
           ! see TRANSP_V.FOR


C =================== LOCAL VARIABLES ====================
      INTEGER IB, LSURF, LBOTTOM, L
      real*8 C_CORR, norm_volume, NEW_VOLUME, OLD_VOLUME

$IF !DEFINED DEBUG_MODE
C To avoid warning:
      CHARACTER*1 CNAME
      if (DEBUG) CNAME=NAME
$else
!      INCLUDE 'DEBUG.INC'
      INTEGER I,K,N
      if (DEBUG) THEN
         WRITE(DEBUG_UNIT, '(A,A)' ) ' ===== ADJ_CONC:', NAME
         WRITE(DEBUG_UNIT,'(1x, A5, 3A15)')
     &            'I','DERIV(I)','CONC(I)','VLAYER(I)'
         DO I = 1, INDXI(NBI+1)
            WRITE(DEBUG_UNIT,'(1x, I5, 4G15.7))')
     &           I, DERIV(I), CONC(I), VLAYER(I)
         END DO
      endif
$endif

C    ------ Update concentrations for last time-step,
C           possibly with shorter time-step for vertical diffusion:
      CALL UPDATE_CONC ( TRCALC, DEBUG, TSTEP, NAME,
     &        DERIV, CONC, CTRL_C )

C    This will also take care of changes in volume over past time-step.


C ---------- For each internal basin:
      DO IB = 1,NBI
         LSURF = INDXI(IB)+1
         LBOTTOM = INDXI(IB+1)


C     --------- Adjust concentrations in layers with varying volume
C               to achieve mass conservation between basins:
         NORM_VOLUME = VTOTZ(IB)*VFROPN(IB)
         NEW_VOLUME  = NORM_VOLUME + DVTOT(IB)
         OLD_VOLUME  = NORM_VOLUME + DVTOTL(IB)
         C_CORR = OLD_VOLUME/NEW_VOLUME

$if defined DEBUG_MODE
         if (DEBUG) THEN
            WRITE(DEBUG_UNIT,'('' ==== Basin '',I5)') IB
            WRITE(DEBUG_UNIT,*)' After UPDATE_CONC:'
            WRITE(DEBUG_UNIT,'(1x,A5,2A15)')
     &              'I','CONC(I)'
            DO I = LSURF, LBOTTOM, 4
               N = MIN ( 3, LBOTTOM-I)
              WRITE(DEBUG_UNIT,'(4(I5,1X,G14.7))') (I+K,CONC(I+K),K=0,N)
           END DO
         endif
$endif
         DO L = LSURF, LSURF+NLVOPN(IB)-1
            CONC(L) = CONC(L)*C_CORR
         END DO

$if defined DEBUG_MODE
         if (DEBUG) THEN
            WRITE(DEBUG_UNIT,*)' After volume correction C_Corr ', C_Corr, ":"
            DO I = LSURF, LBOTTOM, 4
               N = MIN ( 3, LBOTTOM-I)
              WRITE(DEBUG_UNIT,'(4(I5,1X,G14.7))') (I+K,CONC(I+K),K=0,N)
           END DO
         endif
$endif


      END DO

      END SUBROUTINE

C ----------------------------------------------------
C Subroutine MIX_CONC homogenizes concentrations
C for a given substance over well-mixed layers,
C indicated by homogeneous density.

      SUBROUTINE MIX_CONC ( DEBUG, NAME,
     &               NBI, INDXI, XMIX, LMIX, NLI, VLAYER, DENS, CONC )

C In:  NBI, XMIX, INDXI, NLI, VLAYER
C In/Modified out: CONC

C --------------- External arguments: -----------------
      LOGICAL DEBUG
C            - Controls debug output
      CHARACTER*(*) NAME
C            - Name of substance, used in debug dumps
      INTEGER NBI
C            - maximum and actual number of internal basins
      real*8  XMIX (NBI)
C            - number of lower layer in well mixed volumes,
C              with fractional part indicating
C              degree of mixing for last layer (see SURFMX)
      Integer LMIX (NBI)
C            - index of lowest layer involved in homogenizing
      INTEGER INDXI(NBI+1)
C            - layer index limits for internal basins.
      INTEGER NLI
C            - number of layers in internal and external basins + conn
C                           (unit m3*CONSI)
      real*8 VLAYER(NLI)
C            - volume of each layer
      real*8 DENS(NLI)
C            - density of each layer
      real*8 CONC(NLI)
C            - Concentrations in internal basins
C                   in:  preliminary stepwise integrated values
C                   out: corrected values.



C =================== LOCAL VARIABLES ====================
      INTEGER IB, LSURF, LBOTTOM, L, LT, LB
      real*8 SUM_OF_MASS, V_SUM, X, XR
      real*8 CONC_MIXED
      LOGICAL HOMOGENEOUS

$if defined DEBUG_MODE
      real*8 SUM_BEFORE, SUM_AFTER
      INTEGER I,K,N
      LOGICAL WARNING

      WARNING = .false.

      if (DEBUG)
     &   WRITE(DEBUG_UNIT, '(A,A)' ) ' ===== MIX_CONC:', NAME
$else
C          To avoid compiler warning:
      CHARACTER*1 CNAME
      IF (DEBUG) CNAME=NAME
$endif

C ===============================================
C                Scan basins:
      DO IB = 1,NBI
C ===============================================

$if defined DEBUG_MODE
      if (DEBUG)
     &    WRITE(DEBUG_UNIT,'('' Basin nr. '',I5)') IB
$endif

         X = XMIX(IB)
         IF (X .le. 1.0) THEN

$if defined DEBUG_MODE
      IF (DEBUG)
     &    WRITE( DEBUG_UNIT,*)'  - unchanged, no mixing'
$endif
             CYCLE
         ENDIF

         LSURF = INDXI(IB)+1
         LBOTTOM = INDXI(IB+1)

$if defined DEBUG_MODE
      CALL SUM_MASS(NBI, IB, NLI, INDXI, VLAYER, CONC, SUM_BEFORE)
      if (DEBUG) THEN
         DO I = LSURF, LBOTTOM, 4
             N = MIN ( 3, LBOTTOM-I)
             WRITE(DEBUG_UNIT,'(4(I5,1X,G14.7))') (I+K,CONC(I+K),K=0,N)
         END DO
         WRITE(DEBUG_UNIT,'('' Sum before mixing:'',G16.8)' ) SUM_BEFORE
      ENDIF
$endif


         L = LSURF
         SUM_OF_MASS = 0.0
         V_SUM = 0.0

C     ...... sum over completely mixed layers at surface:
         DO WHILE ( X .GT. 1.0 .and. L .le. LBOTTOM  )
            SUM_OF_MASS = SUM_OF_MASS + CONC(L)*VLAYER(L)
            V_SUM = V_SUM + VLAYER(L)
            X = X - 1.0     ! Remaining number of layers
            L = L + 1       ! Next layer to mix into surface volume
         END DO  ! On EXIT: X = mixing fraction >0 of last layer

C     ...... handle partially mixed layer L
C            if mixing goes below boundary of surface layer:
         IF ( L . GT. LSURF ) THEN
            CONC_MIXED = SUM_OF_MASS / V_SUM ! preliminary value
C        ...... mix in the specified fraction of next layer by adjusting
C           a. ... concentration C1 in well-mixed layers:
C             C1'= C1 + X*V2*(C2-C1)/(V1+V2) = C1*(1-X*R)+C2*X*R
C           b. ... concentration C2 in partially mixed layer:
C             C2'= C2 + X*V1*(C1-C2)/(V1+V2) = C1*X*R1   +C2*(1-X*R1)
C                                            = C1*X*(1-R)+C2*(1-X*(1-R))
C                  with R = V2/(V1+V2), R1 = V1/(V1+V2) = 1-R
            XR = X*VLAYER(L)/(V_SUM + VLAYER(L))
            CONC(L-1)  = CONC_MIXED*(1.0-XR) + CONC(L)*XR
            CONC(L)    = CONC_MIXED*(X-XR)   + CONC(L)*(1.0-X+XR)
C        ...... reset values in completely mixed layers:
            L = L-1
            DO WHILE (L .GT. LSURF)
                 L = L-1
                 CONC(L) = CONC(L+1)
            ENDDO
         ENDIF

         DO WHILE ( X .GT. 1.0 .and. L .le. LBOTTOM  )
            SUM_OF_MASS = SUM_OF_MASS + CONC(L)*VLAYER(L)
            V_SUM = V_SUM + VLAYER(L)
            X = X - 1.0     ! Remaining number of layers
            L = L + 1       ! Next layer to mix into surface volume
         END DO  ! On EXIT: X = mixing fraction >0 of last layer


$if defined DEBUG_MODE
      CALL SUM_MASS( NBI, IB, NLI, INDXI, VLAYER, CONC, SUM_AFTER )
      IF ( ABS( SUM_AFTER - SUM_BEFORE) .gt. 1.0E-6*MAX(ABS(SUM_BEFORE),ABS(SUM_AFTER)) ) THEN
         WARNING = .true.
         WRITE(DEBUG_UNIT,'(/" !!!!! non conservative mixing at T=",G16.8)') T
         WRITE(DEBUG_UNIT,'('' XMIX='',g15.7)' ) XMIX(IB)
         WRITE(DEBUG_UNIT, '(A,A," basin ",I5)')' in MIX_CONC:', NAME,IB
         WRITE(DEBUG_UNIT,*) ' Sum before mixing: ', SUM_BEFORE
         WRITE(DEBUG_UNIT,*) ' Sum after mixing:  ' , SUM_AFTER
      ENDIF
      if (DEBUG .or. WARNING) THEN
         DO I = LSURF, LBOTTOM, 4
             N = MIN ( 3, LBOTTOM-I)
             WRITE(DEBUG_UNIT,'(4(I5,1X,G14.7))') (I+K,CONC(I+K),K=0,N)
         END DO
      ENDIF

$endif


C ===== HOMOGENIZE INSTABILITES ================
C   Scan layers from surface to bottom layer
      LT = LSURF
      DO WHILE ( LT .lt. MIN(LMIX(IB),LBOTTOM ) )
C ==============================================
         LB = LT + 1
         HOMOGENEOUS = DENS(LB) .LE. DENS(LT)
         IF ( HOMOGENEOUS ) THEN
            V_SUM = VLAYER(LT)
            SUM_OF_MASS = V_SUM*CONC(LT)
            DO WHILE ( HOMOGENEOUS )
               SUM_OF_MASS = SUM_OF_MASS + CONC(LB)*VLAYER(LB)
               V_SUM = V_SUM + VLAYER(LB)
               LB = LB+1
               if ( LB .gt. MIN(LMIX(IB),LBOTTOM ) ) EXIT
               HOMOGENEOUS = DENS(LB) .LE. DENS(LT)
            ENDDO
            CONC_MIXED = SUM_OF_MASS / V_SUM
            DO L = LT,LB-1
               CONC(L)  = CONC_MIXED
            ENDDO
         ENDIF
C ===========================================
         LT = LB
      ENDDO   ! Next layer/mixed region
C ===========================================


C ===========================================
      END DO  ! Next basin
C ===========================================
      END SUBROUTINE



$if defined DEBUG_MODE
      SUBROUTINE SUM_MASS(NBI, IB, NLI, INDXI, VLAYER, CONC, SUM)
      
      
      INTEGER NBI, IB
      INTEGER INDXI(NBI+1)
      INTEGER NLI
      real*8 VLAYER(NLI)
      real*8 CONC(NLI)
      real*8 SUM

      INTEGER I

      SUM = 0.0
      DO I = INDXI(IB)+1,INDXI(IB+1)
         SUM = SUM + VLAYER(i)*CONC(i)
      ENDDO
      END SUBROUTINE
$endif



      End Module