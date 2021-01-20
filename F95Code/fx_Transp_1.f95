      Module fx_Transp_1
      
      use ModelDimensions
      use ModelParam_Inputs
      use fx_RunControl
      use fx_Transp_2
      use fx_Transp_u
      use fx_Transp_v
      
      implicit none


C =====================================================================
C Fjord model
C   Transport calculations and mass --> concentration calculations
C               May 1990, Birger Bjerkeng, NIVA.
C               Last revision: August 1990.
C =====================================================================
C   Top Module in transport calculation.
C   Contains two entry points called from ACSL program EUTRO:

C   Subroutine WATER_TRANSPORTS allocates transport arrays (at first call),
C   calculates transports of water within model at given time
C   by calling TRANSP_SETUP in Module TRANSP_2
C   returns necessary model variables to ACSL derivative section,
C   and stores values internally for use when entry MASS_TRANSPORT
C   is called later in the same derivative evaluation.

C   Entry point MASS_TRANSPORT (called via MTRANS subroutine in EUTRO)
C   will call TR_CONC_DERIV in Module TRANSP_2
C   with transport description derived in TRANSP
C   and external parameters giving concentration of some substance,
C   to calculate time derivative of the concentration.
C =====================================================================


C ====================================================================
C          Internal allocated arrays for transport description
C ====================================================================

            ! Water transport arrays later used by TR_CONC_DERIV:

            ! Dimension parameters:

      integer MLV  ! >=INDXI(NBI+1) required
          parameter (MLV = dimMLI)

      integer MLH
            ! >=(4*sum over connections of (NL_INVOLVED+1)
            ! where NL_INVOLvED for each connection
            !     = Sum of layer connections (def. by INDXC)
            !     + sum of layers in each of involved basin
            !       or boundary (defined by INDXI, INDXE).
            ! Requirement exceeded by:
          parameter (MLH = 4*(dimMLC+dimMC*(dimMLI+dimMLE)))

      INTEGER M_TRCF
            ! max. number of transfer coefficients, must exceed
C           ! sum over basins of squared number of layers, always
            ! within square of total number of layers.
            ! (Used by TRNADJ and TR_UPDATE_CONC below)
          parameter (M_TRCF = dimMLI*dimMLI+dimMBI+1)

      integer MW
            ! requirements:
            ! >= total number of layers over basins
            ! >= the number of layers in any boundary
            ! >= the number of layers in any basin connection
            ! >= 2*the number of layers in any basin
          parameter (MW = 2*dimMLI+dimMLC+dimMLE)


      integer ML_C_MAX
            ! > number of layers in any connection:
          parameter (ML_C_MAX=dimMLC+1)


      integer MT_TEMP
            ! Length of intermediate description,
            ! only used for outflows during Setup
            ! TR_TEMP used to store import terms during CONC_DERIV
            ! Required: >=number of basins
            !           >=2*(number of layers in any connection+1)
          parameter (MT_TEMP = dimMBI+2*ML_C_MAX)


C "actual" transport values after taking buffering volumes
C  into consideration: ==> effective mass transport between basins:

      real*8 VTRANS          (2,MLV)
C           - Vertical transports internal basins,
C             Final values:
C                   (1,I) m3/day diffusive flux between layer I and I-1
C                   (2,I) m3/day advective flux up from layer I to I-1
C                         (see TRANV1 and TR_CONC_DERIV)

      real*8 HTRANS          (MLH)
C           - Horizontal transport across connections
C                 Ref. subroutine TRANH and MTRAN1

      INTEGER HTR_L          (MLH,2)
C           - Index of participating layers for flows HTRANS.
C                      Cfr. subroutine TRANH, TRANC, and MTRAN1.
C               Indexes are set local to each basin,
C               with 1 for surface layer.

      INTEGER INDXH          (dimMC+1)     ! MC+1
C           - Buffer segmenting of HTRANS and HTR_L
C             between connections, set by TRANH.


      real*8 TRCF          (M_TRCF)
C           - Transfer coefficients for diffusive mixing,
C             set up by TRANV2, called below.

      INTEGER TRCF_RANGE   (0:1,MLV)
C           - active range of TRCF for the layers

C -------- work arrays for TRANSP_SETUP:

      real*8 VMIX            (2,MLV)
C           - Energy budget for vertical mixing, TRANSP_2



      real*8  WORK_ARRAY       (0:MW)
C            - Scratch array for storing flux values locally in TRANH.
C              and as dummy array in calls to MTRAN1.
C              Length set to 0..N,
C                with N >= total sum of internal layers
C                       >= sum of layer connections,
C                giving more than enough space for it's various uses.

      real*8  TR_TEMP          (MT_TEMP)
      INTEGER TR_SOURCE_L      (MT_TEMP)
      real*8  CONNECT          (ML_C_MAX,2)
C            - Intermediate transport description used in TRANH

      real*8  ZSURF_W          (dimMBI)
C              double precision internal values for surface levels

      real*8  DZ_DT            (dimMBI)
C           internal values for surface level changes from surface
C           level adjustment in ITERATE_SURFACES

      real*8 T_QTRANSP_CALC

      LOGICAL HTR_CALC

      real*8 TIME_STEP_USED

      INTEGER, Private:: I, K, L, IC, IS, NLH, NL_C, NL_AB_MAX, NL_C_MAX, NW, NQ_TR
      
      
      contains


$define VERSION_2
C
C Eutrofimodell Indre Oslofjord - Fil  TRANSP_1.FOR

C NB!  Burde ha beskyttet beskrivelsen bedre mot 
C      utillatelige endringer i variable av bruker,
C      kan gi inkonsistent beskrivelse.


$undefine DEBUG_MODE

$if defined DEBUG_MODE
$define DEBUG_MODE_GT_1
$endif

C                           also controlled by argument TransportTest
C                  IF =2: EXTENDED CONTROL + TIME-STEP DEBUG

!!      SUBROUTINE WATER_TRANSPORTS( HTROFF, TransportTest, ITRZ,
!!     &          DEPTH, MBI, NBI, VDYN, ZBOTMI,
!!     &          DTJETM, QWSURF, MS, BASINQ, DEPTHQ,
!!     &          QJET, QTEMP, QSPP, dDens_dSPP, RNFNDX, QDIAM, NHOLES,
!!     &          XMIX, INDXI, NLI, AREA, SAL, TEMP, SPP, DENSI,
!!     &          VLAYER, VFRAC, NLVOPN, VFROPN, VLCORF,
!!     &          MIXFAC, MIXEXP, N2LIM,
!!     &          SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
!!     &          NBE, INDXE, NLE, DENSE, ZSURFE, DZDTX,
!!     &          NC, INDXC, NLC, BConn1, BConn2, 
!!     &          ZSILL, WIDTH, DPEFF, HTRMIX,
!!     &          MLC, VBUF, VBUFMX, VBUFTR, TCVBUF, MAX_DT, T,
!!     &  UFLOW, VTOTDV, VDYNDV, VBUFDV, BSFLUX, TSTEP, RCQNDX,
!!     &  Ambient_Volume_Flux, Neutral_Depth,
!!     &  BWFREQ, ZMID )


      SUBROUTINE WATER_TRANSPORTS( TransportTest, MBI, MS, 
     &          QJET, DENSE, MLC, MAX_DT, TSTEP)

      
C -------------------- External arguments: ------------------------
C                      transferred to TRANSP_SETUP, see definition there
!!      LOGICAL HTROFF, ITRZ

      LOGICAL TransportTest

!!      real*8 DEPTH(*)

      INTEGER MBI

!!      INTEGER NBI
!!      real*8 VDYN(2,NBI)
!!      real*8 ZBOTMI(NBI)
!!
!!      real*8 DTJETM
!!      real*8 QWSURF(NBI)

      INTEGER MS

!!      INTEGER BASINQ(MS)
!!      real*8 DEPTHQ(MS,2)

      real*8 QJET(MS,2) ! volume flux of discharge (1)
                        ! and recipient water mixed into the discharge (2)

!!      real*8 QTEMP(MS)  ! Temperature flux (deC*m3/s)
!!      real*8 QSPP(MS)   ! flux of particles (G/S)
!!      real*4 dDens_dSPP ! net effect of particle conc. on density
!!      INTEGER RNFNDX(MS,2)
!!      real*4 QDIAM(MS)
!!      INTEGER NHOLES(MS)
!!
!!      real*8 BSFLUX (NBI)
!!C          - buoyancy flux into surface layer from horisontal transports
!!C            in  unit (kg/m2/s), set by TRANH
!!
!!      real*8 XMIX(NBI)
!!
!!      INTEGER INDXI(NBI+1)
!!      INTEGER NLI
!!      real*8    AREA(NLI)
!!      real*8    SAL    (NLI)
!!      real*8    TEMP   (NLI)
!!      real*8    SPP    (NLI)
!!      real*8    DENSI(NLI)
!!      real*8    VLAYER(NLI)
!!      real*8    VFRAC  (NLI)
!!      INTEGER NLVOPN(NBI)
!!      real*8    VFROPN (NBI), VLCORF(NBI)
!!      real*8    MIXFAC
!!      real*8    MIXEXP, N2LIM
!!      real*8    SFMIXC(MBI), SFMIXZ(MBI,2)
!!      real*8    GMIXFR(NBI)
!!      real*8    GMIXDC, GMIXDX
!!
!!      INTEGER NBE
!!      INTEGER INDXE(NBE+1)
!!      INTEGER NLE

      real*8    DENSE(NLE)

!!      real*8    ZSURFE(NBE), DZDTX(NBE)
!!
!!      INTEGER NC
!!      INTEGER INDXC(NC+1)
!!      INTEGER NLC
!!      INTEGER BConn1(NC), BConn2(NC)
!!      real*8    ZSILL (NC)
!!      real*8    WIDTH(NLC), DPEFF(NC)
!!      real*4    HTRMIX(NC)! DEGREE OF MIXING BETWEEN CONTIGUOUS TRANSPORTS
!!                        ! IN SAME DIRECTION
!!
      INTEGER MLC

!!      real*8    VBUFMX(2,MLC), VBUFTR(2,MLC)
!!      real*8    VBUF(2,MLC)
!!      real*8    TCVBUF(2,NC)

      real*8 MAX_DT

!!      real*8 T

!!      real*8 UFLOW (MLC)
!!
!!      real*8 VTOTDV (MBI)
!!C Change of dynamic volume related to surface level:
!!      real*8 VDYNDV (2, MBI)
!!C             (1,I): as net result of final transport calculation
!!C             (2,I): as result of surface level iteration.
!!C         Deviation is a measure of the accuracy and consistency of
!!C         the transport/surface level relation in the two Modules.
!!      real*8 VBUFDV (2, MLC)

      real*8 TSTEP

!!      INTEGER RCQNDX(MS) ! Final receiving layer for outlet QJET
!!C                          taking mixing/entrainment into consideration
!!
!!      real*8 Ambient_Volume_Flux (MS) ! volume mixed into dived outlets
!!      real*8 Neutral_Depth       (MS) ! neutral depth for dived outlets
!!
!!      real*8 BWFREQ (NLI) ! Stability, as Brunt_w„isela frequency
!!      real*8 ZMID   (NLI) ! Mean depth of layers (volume weighted)



C ----------- additional arguments in ENTRIES
C               TRNADJ, MASS_TRANSPORT, TR_UPDATE_CONC
!!      LOGICAL TRCALC, VTRNEG
!!      real*8 CTRLC  (NLI) ! substance used to control adv/diff. compensation
!!                        ! tool to ensure paralell transports of fytoplankton
!!                        ! components
!!
!!



C ===================================================================

C local variables:
      logical OK

      INTEGER LENGTH_TRCF, NT_TEMP, N1,N2,NL_INVOLVED

      
      SAVE   ! (All allowable variables)

$if defined DEBUG_MODE
!      INCLUDE 'DEBUG.INC'
      INTEGER I
      INTEGER NDEBUG /0/
$endif

      LOGICAL DEBUG
      DEBUG = TransportTest


$if defined DEBUG_MODE_GT_1
      WRITE(DEBUG_UNIT,*) 'DEBUG:', DEBUG
      IF (TransportTest) THEN
        if (NDEBUG.LE.0) THEN
          WRITE(*,*) '---------- TRANSP -------------'
  10      WRITE(*,*) ' Test print? (0=no, N>0=N steps, -1=turn off)'
          READ (*,*,ERR=10) NDEBUG
          if (ndebug.lt.0) TransportTest = .FALSE.
          if (ndebug.eq.0) then
             DEBUG = .FALSE.
          ENDIF
        endif
        ndebug = max(0,ndebug-1)
      ELSE
          DEBUG = .FALSE.
          ndebug = 0
      ENDIF
      IF (DEBUG) WRITE(DEBUG_UNIT,*)' >>>>>>> TRANSP at T = ',T
$endif
C        USED BELOW TO CONTROL DIAGNOSTICS IN OTHER ENTRIES


C Execute TRANSP: =================================================


C     ..... Rescale tidal mixing coefficient to apply to N^2=1,
C           and multiply by energy factor:
C                       (used in subroutine TRANV1)
         DO I = 1,NBI
C        ......... tidal mixing:
            MIXCONST(I) =  EMIXRL*MIXCF(I)*(N2SCAL**(MIXEXP/2.0))
         ENDDO



C ...... Calculate required length of internal, allocated arrays:
C        Only needed for first call after initiation, as
C        as safeguard against limits of arrays


C ...........Connected to horizontal transports:
      NL_C_MAX = 0
      NLH = 0
      NW = INDXI(NBI+1)  ! <=dimMLI

C     ..... if horisontal transport calculation is not inactivated:
      DO IC = 1,NC
            N1 = BCONN1(IC)
            N2 = BCONN2(IC)
            NL_C = (INDXC(IC+1)-INDXC(IC))                 !<=dimMLC
            IF (N2.GT.0) THEN
               NL_AB_MAX = MAX ( INDXI(N1+1)-INDXI(N1),
     &                           INDXI(N2+1)-INDXI(N2) )   !<=dimMLI
            ELSE
               N2 = ABS(N2)
               NL_AB_MAX = MAX ( INDXI(N1+1)-INDXI(N1),
     &                           INDXE(N2+1)-INDXE(N2) )   !<=max(dimMLI, dimMLE)
            ENDIF
            NL_INVOLVED = NL_C + NL_AB_MAX          !<=dimMLC+max(dimMLI,dimMLE)
            NL_C_MAX = MAX ( NL_C_MAX, NL_C )       !<=dimMLC
            NW = MAX ( NW, NL_AB_MAX, NL_C )        !<=max(dimMLI,dimMLE, dimMLC)
                           ! (should always have NL_AB_MAX>=NL_C)
            NLH = NLH + NL_INVOLVED +1              !<=dimMLC+dimMC*max(dimMLI,dimMLE)
      ENDDO
      NL_C_MAX = NL_C_MAX+1                         !<=dimMLC+1
                         ! May connect 1. layer below sill

C  max. number of transfer coefficients:
C          sum of sqared number of layers within basins:
C  (Used by TRNADJ and TR_UPDATE_CONC below)
      LENGTH_TRCF = 0
      DO K=1,NBI
         L = INDXI(K+1)-INDXI(K)                    ! <=dimMLI
         LENGTH_TRCF = LENGTH_TRCF + L*L            ! <=dimMLI*dimMLI
                                   ! sum(xi*xi)<=N*N whensum(xi)=N                                                    !
         NW = MAX(NW, 2*L ) !<=max(dimMLI,dimMLC,2*dimMLI)
                            ! dvs. <=max(dimMBI,dimMLI)+max(dimMLI,dimMLE)
      ENDDO

C   Max. 2 transports pr. involved layer for each connection:

      NT_TEMP = MAX( NBI, 2*NL_C_MAX )    !<=max(dimMBI,dimMLC+1)
C         - Intermediate description, only outflows during Setup
C         - TR_TEMP used to store import terms during CONC_DERIV

      NLH = 4*NLH     ! Final HTRANS description, all involved flows,
C                       and possible residual mixing flows.
C                            (cfr. TRANH )
                  !<=4*( dimMLC+dimMC*max(dimMLI,dimMLE) )


C .......... Connected to dived outlets (assumed to rise to surface):
C            Sum of possible depth entrainment transports for dived
C            outlets with outlet density < ambient density.
      NQ_TR = 1  ! To avoid error in allocate
      DO IS = 1,MS
          k = BASINQ(IS)
          if( k.gt.0. and. k.le.NBI) then 
             NQ_TR = NQ_TR + RNFNDX(IS,1) -(INDXI(k)+1) + 1
          endif
      ENDDO

      OK=.true.
      OK = CheckDimension("MT_TEMP" , MT_TEMP , NT_TEMP     )
      if (OK) OK=CheckDimension("ML_C_MAX", ML_C_MAX, NL_C_MAX    )
      if (OK) OK=CheckDimension("MW"      , MW      , NW          )
      if (OK) OK=CheckDimension("MLH"     , MLH     , NLH         )
      if (OK) OK=CheckDimension("MQ_TR"   , MQ_TR   , NQ_TR       )
      if (OK) OK=CheckDimension("LMV"     , MLV     , INDXI(NBI)+1)
      if (OK) OK=CheckDimension("M_TRCF"  , M_TRCF  , LENGTH_TRCF )
      if (.not.OK) then
         WRITE (*,*) '==>Insufficient array dimensions in TRANS_1'
         WRITE (*,*) 'run aborted'
         STOP
      ELSE
             ! Signal that TRANSP_SETUP should calculate QTRANSP after allocation,
             ! at restart ( T decreased from last calculation),
             ! and later about two times a day (triggered if T_QTRANSP_CALC = T)
         if ( T_QTRANSP_CALC.gt. T .or. T_QTRANSP_CALC.le. T - DTJETM )
     &        T_QTRANSP_CALC = T
      ENDIF

      HTR_CALC = .NOT. HTROFF

C ----------- Call subroutine in TRANSP_2 to calculate water transports:
C             Deactivated horisontal transports signalled by NHL = 0.

      TSTEP = MAX_DT ! Further reduced in TRANSP_SETUP if necessary
      CALL HELLO ('In TRANSP_1: Call to TRANSP_Setup:' )

      CALL TRANSP_SETUP( HTR_CALC, DEBUG, ITRZ,
     &         DEPTH, MBI, NBI, VDYN, ZBOTMI,
     &         T_QTRANSP_CALC - T, QWSURF, MS, BASINQ, DEPTHQ,
     &          QJET, QTEMP, QSPP, dDens_dSPP, RNFNDX, QDIAM, NHOLES,
     &          MQ_TR, XMIX,
     &          INDXI, NLI, AREA, SAL, TEMP, SPP, DENSI, VLAYER, VFRAC,
     &          NLVOPN, VFROPN, VLCORF,
     &         MIXFAC, MIXEXP, N2LIM,
     &         SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
     &         NBE, INDXE, NLE, DENSE, ZSURFE, DZDTX,
     &         NC, INDXC, NLC, BConn1, BConn2, 
     &         ZSILL, WIDTH, DPEFF, HTRMIX,
     &         MLC, VBUF, VBUFMX, VBUFTR, TCVBUF, UFLOW,
     &         VTOTDV, VDYNDV, VBUFDV, TSTEP,
     &  MLV, VTRANS, VMIX, MLH, HTRANS, HTR_L, dimMC, INDXH,
     &  NQ_TR, QJET_TRANSP, QJET_TRANSP_INDEX, RCQNDX,
     &  Ambient_Volume_Flux, Neutral_Depth,
     &  MW, WORK_ARRAY, MT_TEMP ,TR_TEMP, TR_SOURCE_L,
     &  ML_C_MAX, CONNECT, dimMBI, ZSURF_W, DZ_DT, BSFLUX,
     &  BWFREQ, ZMID )

!!      CALL TRANSP_SETUP( MLV, VTRANS, VMIX, MLH, HTRANS, HTR_L, dimMC, INDXH,
!!     &  RCQNDX,
!!     &  Ambient_Volume_Flux, Neutral_Depth,
!!     &  MW, WORK_ARRAY, MT_TEMP ,TR_TEMP, TR_SOURCE_L,
!!     &  ML_C_MAX, CONNECT, dimMBI, ZSURF_W, DZ_DT, BSFLUX,
!!     &  BWFREQ, ZMID )

      CALL HELLO ('In TRANSP_1: back from TRANSP_Setup:' )


        NLH = INDXH(NC+1)

$if defined DEBUG_MODE_GT_1
        write(DEBUG_UNIT,*) 
     &    ' after Transp_Setup: NLH=',NLH,' NC=',NC,
     &    ' INDXH=',(INDXH(I),I=1,NC+1)
$endif

$if defined DEBUG_MODE_GT_1
       if (debug) then
          write(*,*) 'From TR_CONC_DERIV with TSTEP=',TSTEP
       endif
$endif

      RETURN  ! ---------------------------------- From call to TRANSP

      end subroutine

C =====================================================================
      subroutine TRNADJ (NBI, INDXI, NLI, VLAYER,
     &                   TSTEP, TRCALC, VTRNEG )
      INTEGER NBI
      INTEGER INDXI(NBI+1)
      INTEGER NLI
      real*8  VLAYER(NLI)
      real*8  TSTEP
      LOGICAL TRCALC, VTRNEG


C Called from main Module, before the subroutine TR_CONC_DERIV
C in TRANSP_U.FOR is used for any substance by calls to MASS_TRANSPORT.

      IF (TRCALC) THEN
C           - transports are active:

         CALL TRANV2( DEBUG, TSTEP, NBI, INDXI, NLI, VLAYER,
     &        MLV, VTRANS, TRCF_RANGE, M_TRCF, TRCF, VTRNEG )

C TRANV2 in TRANSP_V.FOR makes final calculations and adjustments of
C transport description to work for a given time_step in TR_CONC_DERIV.
C Returns VTRNEG = .TRUE. if negative net diffusion occurred.


         TIME_STEP_USED = TSTEP
C Stores actual time-step, will be used in calculating
C mean time derivatives in subsequent calls to mass_transport,
C for Version .ge.2


$if defined DEBUG_MODE
       IF(DEBUG) THEN
          write(DEBUG_UNIT,*)' ==== TRNADJ, endelig transportbeskr.'
          WRITE (DEBUG_UNIT,'(1X,A5,2A15)' )
     &     'K','VTRANS(1,K)','VTRANS(2,K)'
          WRITE(DEBUG_UNIT,'(1X,I5,2G15.7)')
     &     (K,(VTRANS(I,K),I=1,2),K=1,MLV)
          IF ( HTR_CALC .AND. NLH.GT.0 ) THEN
             WRITE(DEBUG_UNIT,'(1X,A5,3A15)' )
     &          'K','HTR_L(K,1)','HTR_L(K,2)','HTRANS(K)'
             WRITE(DEBUG_UNIT,'(1X,I5,2I15,G15.7)' )
     &          (K,(HTR_L(K,I),I=1,2),HTRANS(K), K=1,NLH)
          ENDIF
       ENDIF
$endif

      ENDIF

      end subroutine


C =====================================================================
C Initiate time derivatives for water concentrations:
C If version.ge.2: not including effect of diffusion.
      subroutine MASS_TRANSPORT( MDEBUG, 
     &                      GROUP_NUMBER,
     &                      CONCI, CONCE, CDERIV, IMPORT, NAME )

      LOGICAL MDEBUG
      INTEGER GROUP_NUMBER
      real*8 CONCI  (NLI)
      real*8 CONCE  (NLE)
      real*8 CDERIV (NLI)
      real*8 IMPORT (NBI)
      CHARACTER*(*) NAME
      
      

!!      subroutine MASS_TRANSPORT()
!! MDEBUG, NBI, INDXI, NLI,
!!     &                      VLAYER, NBE, INDXE, NLE,
!!     &                      NC, INDXC, BConn1, BConn2, HTRMIX,
!!     &                      GROUP_NUMBER,
!!     &                      CONCI, CONCE, CDERIV, IMPORT, NAME )


C ----- For all calls:
C Transfers external arguments together with transport terms into
C subroutine TR_CONC_DERIV in Module TRANSP_U to
C initiate time derivative of concentration, see TRANSP_U for details.
C MDEBUG controls debug output for each variable separately,
C contingent upon TransportTEST being activated (see above).

$if defined DEBUG_MODE_GT_1
      IF (DEBUG) 
     &WRITE(DEBUG_UNIT,'(A/6A4,A8,A4,A6,2A4/6I4,I8,I4,I6,2I4)') 
     & 'TR_CONC_DERIV called with ',
     & 'NBI','NLI','NBE','NLE','NC','MW',
     & 'MT_TEMP','MLV','NQ_TR','MLH','MC',
     &  NBI,  NLI,  NBE,  NLE,  NC,  MW, 
     &  MT_TEMP,  MLV,  NQ_TR,  MLH,  dimMC
$endif

      CALL TR_CONC_DERIV (
     &       HTR_CALC, MDEBUG,
     &       NBI, INDXI, NLI,
     &       VLAYER, NBE, INDXE, NLE,
     &       NC, INDXC, BConn1, BConn2, HTRMIX,
     &       GROUP_NUMBER,
     &       CONCI, CONCE, CDERIV, IMPORT,
     &       MW, WORK_ARRAY, MT_TEMP, TR_TEMP, MLV, VTRANS,
     &       MQ_TR, NQ_TR, QJET_TRANSP_INDEX, QJET_TRANSP,
     &       MLH, HTRANS, HTR_L, dimMC, INDXH, NAME )


      end subroutine


$if defined VERSION_2
C =============================================================
C Integrate concentrations over given step:
      SUBROUTINE UPDATE_CONC (TRCALC, MDEBUG, TSTEP, 
     &     NAME, CDERIV, CONCI, CTRLC )

      LOGICAL TRCALC, MDEBUG
      real*8 TSTEP
      CHARACTER*(*) NAME
      real*8 CDERIV(NLI)
C              - time derivatives for given concentrations
      real*8 CONCI(NLI)
C           - Concentrations in internal basins
C                   in:  preliminary stepwise integrated values
C                   out: corrected values.
      real*8 CTRLC  (NLI)
           ! Controlling substance for applying neg. diffusion
           ! as correction of diffusive effect of advection
           ! see TRANSP_V.FOR
      
      L = NW/2
      CALL TR_UPDATE_CONC ( TRCALC, MDEBUG, NAME,
     &     NBI, INDXI, NLI, VLAYER, VTRANS,
     &     TSTEP, M_TRCF, TRCF, TRCF_RANGE,
     &     L, WORK_ARRAY(0), WORK_ARRAY(L),
     &     CDERIV, CONCI, CTRLC )
      END Subroutine
      
      
$endif


      End Module
