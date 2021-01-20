      Module fx_Transp_h

      use fx_RunControl, only: DEBUG_UNIT
      
      implicit none
C ===================================================================
C Eutrophication model   - File:    TRANSP_H.FOR
C                                   Birger Bjerkeng, NIVA.

C Contains subroutines:
C     TRANH   called from TRANSP_SETUP in TRANSP_2.FOR
C     TRANH2  called from TRANH in this Module.
C     MTRAN1  called from TR2_CONC_DERIV in TRANSP_U.FOR
C           No calls to other Modules

C Mulige forbedringer:
C     1. La fordeling av VBUF f›lge med vertikale forflytninger
C        og vertikal blanding i hovedvannmassen.
C        Dette er antagelig en 2. ordens korreksjon,
C        fordi VBUF er tenkt som et forholdsvis hurtig fluktuerende
C        fenomen, men det kan ha en viss betydning ved langvarige inn-
C        str›mninger. Kan ikke bare la transport av VBUF(IB,L) bestemmes
C        av vertikaltransporter i det bassenget buffervolumet ligger i,
C        det vil ikke fungere for ytre bassenger, hvor det ikke er noen
C        egen beskrivelse av vertikaltransporter eller massebalanse.
C        Mulig l›sning: La VBUF i forbindelse med kobling til ytre
C        bassenger ha indeks knyttet til det indre bassenget, og ha
C        vertikal utveksling knyttet til transporter i dette bassenget.
C        VBUF kunne omfordeles ut fra CONNECT-array, slik at en fikk
C        Fordeling i det ytre bassenget ogs† for bruk her.
C        Kanskje b›r VBUF forenkles til et enkelt array,
C        hvor verdi > 0 bety at vann fra basseng A er flyttet
C        reversibelt til B, og verdi <0 motsatt.
C        VBUF-index kan referere til lag i mottagende volum, dvs.
C        basseng B hvis VBUF>0, A hvis VBUF<0, og bruken av VBUF
C        m† da styres CONNECT. Vertikale endringer av VBUF m†tte
C        styres av en kombinasjon av endringer i bassengene
C        for kobling mellom indre bassenger.

C        Naturlig † se p† dette i forbindelse med en overgang til
C        lagmodell med flytende grenseniv†er mellom lagene?


$undefine DEBUG_TRANH
$undefine DEBUG_MTRAN1
$undefine DEBUG_WIND

$if defined DEBUG_TRANH
$define DEBUG_TRANH_GT_1
$endif

$if defined DEBUG_TRANH
$undefine DEBUG_TRANH_GT_2
$endif

C         IF >0: WRITES DEBUGGING OUTPUT


      contains


C ======================================================================
      SUBROUTINE TRANH( DEBUG, DEPTH, NL_C, WIDTH, HTRMIX,
     &                  ZSURFA, NLA, DENSA, BSFLUX_A, XMIX_A,
     &                  ZSURFB, NLB, DENSB, BSFLUX_B, XMIX_B,
     &                  MLH, ILH, LLH, HTRANS, HTR_L, TRANS_SUM,
     &                  VBUF, VBUFDV, VBUFMX, VBUFTR, TCVBUF,
     &                  MW, WK_ARRAY, ZSILL, UFLOW,
     &                  MTR, TR_TEMP, LBUF_SOURCE, MLCONN, CONNECT )

C Calculates horizontal transport in connection between basins A and B.
C Influx of lighter water to surface of basin accumulated in BSFLUX_x

C ----------------------- External arguments: ------------------------

      LOGICAL DEBUG

      real*8 DEPTH(*)
C       (I)    Limiting depth of layers (dimension > max(NL_C,NLA,NLB)).
      INTEGER NL_C
C       (I)    Number of connected layers (WIDTH dimension)
      real*8 WIDTH(NL_C)
C       (I)    Mean (central) width of each layer
      real*4 HTRMIX
C       (I)    Degree of mixing between homogenized inflows.
      real*8 ZSURFA, ZSURFB
C       (I)    Surface depths of basins
      INTEGER NLA, NLB
C       (I)    Number of layers in connected basin (>=NL_C).
      real*8 DENSA(NLA), DENSB(NLB)
C       (I)    Density (sigma-t) of layers in connected basins
      real*8 BSFLUX_A, BSFLUX_B            ! kg/s
C       (I/O)  Buoyance influx to each basin to be updated
C              (Only used for internal basins: not for B if EXTCON=.T.)
      real*8 XMIX_A, XMIX_B
C       (I)    Number of layers in wellmixed surface volume,
C              (>0, integral part is nmber of completely mixed layers
C                   fractional part: partial mixing of layer below.

      INTEGER MLH, ILH, LLH
C       (I,I,O): Dim., lower and upper index of HTRANS and HTR_L
      real*8 HTRANS(MLH)                   ! m3/day
C       (I/O)  Horisontal transport terms ( >0:  A--->B, <0: B---->A )
      INTEGER HTR_L(MLH,2)
C       (I/O)  Index of participating layers in basins A and B
C              for corresponding flows in HTRANS.
      real*8 TRANS_SUM                     ! m3/day
C       (O)    Sum of dynamic transports (surface level control)
      real*8   VBUF    ( 2, NL_C )
C       (I)    Current values of buffer volume levels:
C              VBUF(K,I) = part of volume of water basin K, layer I,
C              which is filled with water from the other
C              basin.  Increases to VBUFMX by inflow,
C              decreases to 0.0 by outflow and mixing into
C              main volume. May become <0 due to abrupt change, but
C              the derivative is then set as if the value were zero

      real*8   VBUFDV  (2, NL_C )
C       (O)    Time derivatives of VBUF
      real*8   VBUFMX  ( 2, NL_C )
C       (I)    Maximum value of buffer volume.  If VBUF >= VBUFMX, all
C              all further inflow goes directly to main volume.
      real*8   VBUFTR  ( 2, NL_C )
C       (I)    Transition volume:
C              Outflows becomes effective as transit into the other
C              basin to the degree that VBUF < VBUFTR in source basin
C              Inflows enters as mass transport into the main volume
C              to the degree that VBUF > VBUFMX - VBUFTR in destination

      real*8   TCVBUF  ( 2 )
C       (I)    Mixing time constant for buffer volume.
      integer MW
      real*8  WK_ARRAY (0:MW)
C            - Scratch array for storing flux pr. width,
C              and distribution weights for transports.
C              MW >= 2*MAX(NLA,NLB) >= 2*NL_C
      real*8    ZSILL  ! max. sill depth of connection.

      real*8    UFLOW (NL_C)
C            - Flow velocities, calculated already
C              by call to subroutine FLOW_PROF in Module FLOW_PROF.for 
      INTEGER MTR
      real*8  TR_TEMP(MTR)
C            - Temporary storage of directed transports out of basins
      INTEGER LBUF_SOURCE(MTR), MLCONN
      real*8  CONNECT(MLCONN,2)
C            - Temporary storage of layer index for transport terms


C --------------------- LOCAL VARIABLES -----------------------------

C Constant for calculating pressure from density difference:

      INTEGER IL, NTBUF, K_S, K_D
      real*8 QT, QB, ZT, ZB, ZB_LIM
      real*8 Q_EFF,Q_MIXING
      INTEGER LD_MIN, LD_MAX, LLH_BASE, ID, IS, L, IT, K, IW
      real*8 D_CENTER, D_MIN, D_MAX

      INTEGER KMD, LD_BOTTOM, LD_C
      real*8 TR, W_SUM, D_SUM, W, WW, DD, X
      real*8 DIRECTION

 

C BUOYANCY DATA:
      real*8 BF, BSF(2)
C Temporary storage of mixed and total number of layers:
      INTEGER LMIX(2), LBOTT(2), UPPER_CONNECTED(2)

      real*8 sec_per_day
      parameter ( sec_per_day = 24.*3600.)

C ------------- debug printout control:
      LOGICAL DUMMY


$if defined DEBUG_TRANH
!!      INCLUDE 'DEBUG.INC'
      IF (DEBUG) THEN
          WRITE(DEBUG_UNIT,*) ' ===================== TRANH '
          WRITE(DEBUG_UNIT,'(''ZSURF: '',G16.8,''|'',G16.8)')
     &            ZSURFA,ZSURFB
      ENDIF
$endif
      if ( MW .LT. NLA + NLB) THEN
          STOP 'NW too small in TRANSP_H.TRANH'
      ENDIF



C Internal storage of bouyancy data for use in loop below:
      BSF(1) = BSFLUX_A
      BSF(2) = BSFLUX_B
      LBOTT(1)  = NLA
      LBOTT(2)  = NLB

      DUMMY = DEBUG !to avoid insignificant compiler message (not used)


      TRANS_SUM = 0.0

C :::::::::: Initiate control variables for homogenizing outflows ::::::
      W_SUM = 0
      D_SUM = 0
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


C ===============================================================
C Scan layers, and use the flow profile to set up
C a preliminary description of transports in
C       TR_TEMP = signed volume flux (+: A---> B)   (-: A <---B )
C       LBUF_SOURCE = source layer number relative to surface

$if defined DEBUG_TRANH_GT_1
         IF (DEBUG) WRITE(DEBUG_UNIT, '(/A)' )
     &          '######## set up TR_TEMP description'
$endif

C ===============================================================

      NTBUF = 0  ! Will count number of transports

      DO IL = 1,NL_C


C     --------- Zero buffer volume derivatives to be set below:
         DO K=1,2
            VBUFDV (K,IL) = 0.0
         END DO

C     --------- Prepare flow setup:
         ZB_LIM = MIN( ZSILL, DEPTH(IL+1))
         ZB = DEPTH(IL)
         QB = UFLOW (IL) *WIDTH(IL)


$if defined DEBUG_TRANH_GT_2
         IF (DEBUG) THEN
            WRITE(DEBUG_UNIT,'( '' ----- layer '',I5,A,G12.6)' )
     &            IL,' HTRMIX:', HTRMIX
         ENDIF
$endif

C :::::::: Initiate outflowing density values in CONNECT :::::::::::::::
         CONNECT(IL,1) = DENSA(IL)
         CONNECT(IL,2) = DENSB(IL)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

         DO WHILE ( ZB .LT. ZB_LIM )
C        ...... Outflows:
            ZT = ZB
            QT = QB
            ZB = ZB_LIM
            if (IL.lt.NL_C) then
               QB = UFLOW(IL+1) *WIDTH(IL+1)           ! m2/s
            else
               QB = 0.0
            endif

C               If flow reverses within layer: split in two flows:
            IF (QB*QT .LT. 0 ) THEN  ! Interpolate to zero flow depth:
                ZB = (ZB*QT - ZT*QB)/(QT-QB)
                QB = 0.0
            ENDIF
C         .... Dynamic transport:
            TR = (ZB-ZT)*(QT+QB)/2.0 * sec_per_day   ! m3/day
C                +: from A to B,   -: from B to A
C                ---- converted from m3/s to m3/day ----

C              Summed for surface levels derivatives:
            TRANS_SUM = TRANS_SUM +TR



C          ..... Set basin index of source basin:
            IF ( TR .GT. 0.0 ) THEN
C              ..... FROM A TO B:
               K_S = 1
            ELSE IF ( TR .LT. 0.0 ) THEN
C              ..... FROM B TO A:
               K_S = 2
            ELSE
               CYCLE
            ENDIF

C ::::::::: Sum for calculation of homogenized densities: :::::::::::::
            if (W_SUM.eq.0) LD_MIN = IL
            LD_MAX = IL
            W_SUM = W_SUM + TR
            D_SUM = D_SUM + TR*CONNECT(IL,K_S)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C ::::::::: If last of series of contiguous transports ::::::::::::::
C           Set up homogenized densities in outflowing layers

            if ( ( QB .eq. 0.0 .or. ZB .eq. DEPTH(NL_C+1) )
     &            .and. W_SUM.ne.0.0 ) then

               X = D_SUM/W_SUM*HTRMIX

$if defined DEBUG_TRANH_GT_2
               IF (DEBUG) THEN
                 WRITE(DEBUG_UNIT, '( 2(1x,A,I6:))' )
     &             '  LD_MIN=',LD_MIN,
     &             '  LD_MAX=',LD_MAX
                 WRITE(DEBUG_UNIT, '( 4(1x,A,G14.8:))' )
     &             '  W_SUM=',W_SUM, 'D_SUM=',D_SUM,' X=',X
               ENDIF
$endif
               DO L = LD_MIN, LD_MAX
                  CONNECT(L,K_S) = CONNECT(L,K_S)*(1.-HTRMIX)+X

$if defined DEBUG_TRANH_GT_2
                  IF (DEBUG) THEN
                     WRITE(DEBUG_UNIT, '( A,I6,A,G14.8))' )
     &             '  L:',L, ' density CONNECT(L,K_S)=', CONNECT(L,K_S)
                  ENDIF
$endif
               ENDDO
               W_SUM = 0.0
               D_SUM = 0.0
            ENDIF
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


$if defined DEBUG_TRANH_GT_2
            IF (DEBUG) THEN
              WRITE(DEBUG_UNIT, '( 4(1x,A,G14.8:))' )
     &          '  ZT=',ZT,'QT=',QT, 'ZB=',ZB,'QB=',QB,
     &          '  TR (m3/d)=',TR,
     &          '  W_SUM=',W_SUM, 'D_SUM=',D_SUM
              WRITE(DEBUG_UNIT, '( 2(1x,A,I6:))' )
     &          '  LD_MIN=',LD_MIN,
     &          '  LD_MAX=',LD_MAX
            ENDIF
$endif

C           ..... Part of outflow may be return of water
C                 recently flowed in from the other basin (VBUF).
C                 contributing to decreasing VBUF values.
C                 Rest of flow is stored as real*8 transport into
C                 buffer volume of the other basin:

            if ( VBUF(K_S,IL) .gt. 0.0) then
               if (VBUF(K_S,IL) .lt. VBUFTR(K_S,IL) ) then
                  Q_EFF = TR * ( 1.0 - VBUF(K_S,IL)/VBUFTR(K_S,IL) )
               else ! buffer volume filled above transition volume
                  Q_EFF = 0.0 ! All flow goes to emptying buffer volume
               endif
            else
               Q_eff = TR
            endif

            VBUFDV(K_S,IL) = VBUFDV(K_S,IL) - ABS(TR-Q_EFF)  ! m3/day

$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) WRITE( DEBUG_UNIT,'(A,I1,A,I3,A,2(2(A,E16.10:)/21x))')
     &        ' Source buffer (', K_S, ',', IL, '):'
     &      , 'VBUF  =' ,VBUF  (K_S,IL)
     &      , ' VBUFTR=',VBUFTR(K_S,IL)
     &      , ' VBUFDV=',VBUFDV(K_S,IL)
     &      , ' Q_EFF=', Q_EFF
$endif

C         ..... Store transport and note source layer for outflow:
            IF ( Q_EFF.NE.0 ) THEN
               NTBUF = NTBUF+1
               LBUF_SOURCE (NTBUF) = IL
               TR_TEMP(NTBUF) = Q_EFF

$if defined DEBUG_TRANH_GT_1
       IF (DEBUG) WRITE( DEBUG_UNIT, '( 1x, A,I3, A,E18.11,A,I5)' )
     &         ' ........ TR_TEMP(',NTBUF,')=',TR_TEMP(NTBUF),
     &         ', local layer ',IL
$endif

            ENDIF

         ENDDO ! Next section of layer
      ENDDO ! Next layer

$if defined DEBUG_TRANH_GT_1
      IF (DEBUG) THEN
        write ( DEBUG_UNIT, '(A,E16.10,A)/)' )
     &     ' Net sum of transports:', trans_sum, ' m3/d, ',
     &     ' ', trans_sum/sec_per_day, ' m3/s'
      endif
$endif


C ==================================================================
C Set depth index for connections between basins according to density:
C     CONNECT(IL,K) at this point contains density profiles on both sides,
C         homogenized for contiguous outflows.
C     Below CONNECT(IL,K) is changed into
C         depth index (with fraction) for layer in basin 3-k
C         with outflowing density >= density in layer IL in basin K.
C ==================================================================
$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) write (debug_unit,*)' ....... TRANH2, basin A'
$endif
      CALL TRANH2( DEBUG, NL_C, CONNECT(1,1), NLB, DENSB )
$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) write (debug_unit,*)' ....... TRANH2, basin B'
$endif
      CALL TRANH2( DEBUG, NL_C, CONNECT(1,2), NLA, DENSA )

$if defined DEBUG_TRANH
      IF (DEBUG) THEN
        DO IL=1,NL_C
         WRITE (DEBUG_UNIT,
     &     '('' L='',I4,'':'',G18.11,''-->'',F8.4,
     &       '' | '',F8.4,''<--'',G18.11)')
     &          IL,DENSA(IL),CONNECT(IL,1),CONNECT(IL,2),DENSB(IL)
        END DO
        DO IL=NL_C+1,MIN(NLA, NLB)
         WRITE (DEBUG_UNIT,
     &     '('' L='',I4,'':'',G18.11,12X, ''|'',12X,G18.11)')
     &          IL,DENSA(IL),DENSB(IL)
        END DO
        IF (NLA.GT.NLB) THEN
          DO IL=NLB+1,NLA
            WRITE (DEBUG_UNIT,'('' L='',I4,'':'',G18.11,4X)')
     &          IL,DENSA(IL)
          END DO
        ELSEIF (NLB.GT.NLA) THEN
          DO IL=NLA+1,NLB
            WRITE (DEBUG_UNIT,
     &       '('' L='',I4,'':'',18X,12X,''|'',12X,G18.11)')
     &             IL,DENSB(IL)
          END DO
        ENDIF

      ENDIF
$endif

C ==================================================================
C Depth of well-mixed layer,
C save integral part used in distribution scheme below:
C ==================================================================
      LMIX(1)  = XMIX_A
      LMIX(2)  = XMIX_B
      LMIX(1)  = MAX(1,LMIX(1))
      LMIX(2)  = MAX(1,LMIX(2))

$if defined DEBUG_TRANH
      IF (DEBUG) THEN
         WRITE(DEBUG_UNIT,*) ' LMIX(1) (2):', LMIX
      ENDIF
$endif

C ============================================================
C Begin to set up final transport description to be used for updating
C volume and mass distribution over next time step.
C In this phase the direct irreversible effect of flow at the moment
C are stored as flows in HTRANS, with layer indexes noted in HTR_L.
C ============================================================

      LLH = ILH-1  ! Base index = last occupied position in HTRANS.
                   ! (ILH was entered as input argument to subroutine)

$if defined DEBUG_TRANH_GT_1
      if (debug)
     &    WRITE(DEBUG_UNIT, '(/A)')
     &         '######## direct transports in HTRANS:'
$endif

C ---------------------------------------------------------
      DO IT = 1, NTBUF

$if defined DEBUG_TRANH
  100    CONTINUE     ! RESTART ON ERROR CONDITION BELOW
$endif
C ---------------------------------------------------------

         TR = TR_TEMP (IT)
         IF ( TR .EQ. 0.0) CYCLE


C     ....... Direction indices:
         IF ( TR .GT. 0.0 ) THEN
             K_S = 1  ! from A (=1) to B (=2)
             K_D = 2
             LD_BOTTOM = NLB
         ELSE
             K_S = 2  ! from B to A
             K_D = 1
             LD_BOTTOM = NLA
         ENDIF

         IS = LBUF_SOURCE(IT)         ! Source layer number
         D_CENTER   = CONNECT(IS,K_S) ! Mode of destination distribution
         LD_C   = AINT ( D_CENTER )   ! Upper main layer


C -------------------------------------------------------------------
C Distribution range for TR_TEMP over destination layers
C is extended so that adjacent outflows going to destination layers
C separated from each other are interleaved into intermediate layers.
C Care is taken not to extend to layers beyond mode of the
C depth distribution of the adjacent flow.
C --------------------------------------------------------------------

         D_MIN = D_CENTER
         LD_MIN = LD_C
         IF ( IT.GT.1 ) THEN
            IF (TR_TEMP(IT-1)*TR .GT. 0.0 ) THEN
                D_MIN  = MIN(  D_MIN, CONNECT(IS-1,K_S) )
                LD_MIN = MIN( LD_MIN, LD_BOTTOM-INT(LD_BOTTOM-D_MIN) )
            endif
         ENDIF ! ( ensures LD_MIN = LD_C or >= D_MIN )

         D_MAX  = CONNECT(IS,K_S)
         LD_MAX = LD_BOTTOM - INT(LD_BOTTOM-D_MAX)
         IF ( IT.LT.NTBUF .AND. TR_TEMP(IT+1)*TR .GT. 0.0 ) THEN
             D_MAX = MAX(  D_MAX, CONNECT(IS+1,K_S) )
            LD_MAX = MAX( LD_MAX, INT( D_MAX ) )
         ENDIF ! ( ensures LD_MAX <=LD_C+1 or <=D_MAX )

C ------------ Suspended: ----------------------------------------
C           inflow below well_mixed layer is not distributed
C           into well_mixed layer:
C        IF ( LD_C.gt.LMIX(K_D) .and. D_MIN .le. LMIX(K_D)) THEN
C            D_MIN = LMIX(K_D)+1
C        ENDIF
C           and inflow into wellmixed layer is not distributed
C           to layers below:
C        IF (LD_C.le.LMIX(K_D) .and. D_MAX .gt. LMIX(K_D) ) THEN
C            D_MAX = LMIX(K_D)
C        ENDIF
C -----------------------------------------------------------------


$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) WRITE( DEBUG_UNIT,
     &       '(/1x,A,I3,A,G16.10,3(1x,A,I3)/T15,3(1x,A,F10.5))')
     &    '--- TR_TEMP(',IT,')=',TR
     &   ,'LD_MIN=', LD_MIN, 'LD_C=', LD_C, 'LD_MAX=', LD_MAX
     &   ,'D_MIN= ', D_MIN, 'D_CENTER=', D_CENTER,' D_MAX=', D_MAX
$endif



C -------------------------------------------------------------
C         Distribute TR on destination layers:
C -------------------------------------------------------------


C ---------- Part of flow to layers <= LD_C: ------------------
         W_SUM = TR*(1.0 - MOD(D_CENTER,1.0D0))


C   Distributed according to a triangle from D_MIN to D_CENTER
C   with height growing linearly (2:1) from 0 at D_MIN.
C   Sum of flow to layers LD_min ... ID is fraction of W_SUM
C   given by the fraction of the triangle between D_MIN and ID+0.5
C       (ID always >= D_MIN).

         WW = 0.0
         DD = D_CENTER-D_MIN         ! Length of whole triangle
         DD = DD*DD                  ! Area of whole triangle
         DO ID = LD_MIN, LD_C-1, 1
            X = (ID+0.5-D_MIN)       ! Length of triangle to ID+0.5
            X = X*X/DD*W_SUM         ! Sum of flows to layers <=ID
            WK_ARRAY(ID) = (X-WW)    ! Last flow by difference
            WW = X
         ENDDO

C   Layer LD_C gets the residue (= W_SUM if LD_MIN=LD_C) :

         WK_ARRAY(LD_C) = W_SUM - WW

$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) WRITE( DEBUG_UNIT, '(3(1x,A,G14.8))')
     &       'layers <= LD_C: W_SUM=', W_SUM,'DD=',DD, 'WW=', WW
$endif


C ----------- Part of flow to layers > LD_C: ------------------
         W_SUM = TR - W_SUM
              ! ( if >0: Should always have LD_MAX > LD_C
              !   according to code above )

         if ( LD_max .gt. LD_C) then

C   Distributed according to a triangle from D_MAX to D_CENTER
C   and height growing linearly (2:1) from 0 at D_MAX.
C   Sum of flow to layers (ID ... LD_max) is fraction of W_SUM.
C    = fraction of triangle between ID-0.5 and D_max
C       (ID always <= D_Max).

            WW = 0.0
            DD = D_MAX- D_CENTER       ! Length of whole triangle
            DD = DD*DD                 ! Area of whole triangle
            ID = LD_MAX
            DO while ( ID .gt. LD_C+1 )
               X = (D_max-ID+0.5)      ! Length of triangle to ID-0.5
               X = X*X/DD*W_SUM        ! Sum of flows to layers >=ID
               WK_ARRAY(ID) = (X-WW)   ! Flow to ID by difference
               WW = X
               ID = ID -1
            ENDDO


$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) WRITE( DEBUG_UNIT, '(3(1x,A,G14.8))')
     &       'layers >LD_C: W_SUM=', W_SUM,'DD=',DD, 'WW=', WW
$endif

             WK_ARRAY(ID) = W_SUM - WW
         ENDIF

$if defined DEBUG_TRANH_GT_2
         W_SUM = 0.0
         do id = ld_min, ld_max
             W_SUM = W_SUM + WK_ARRAY(ID)
         enddo
         IF (DEBUG .or. ( ABS(W_SUM-TR) .gt. 1.0e-10*ABS(TR)) ) THEN
            WRITE( DEBUG_UNIT, '(1x,A6,A18)') 'ID', 'WK_ARRAY(ID)'
            WRITE( DEBUG_UNIT, '(1x,I6,G18.11)')
     &            (ID,WK_ARRAY(ID), ID=LD_MIN,LD_MAX )
            WRITE( DEBUG_UNIT, '(3(1x,A6,G18.11))')
     &          'W_SUM=', W_SUM,' vs. TR ',TR,' diff:', W_SUM - TR
         ENDIF
$endif


C     ...... Calculate flows and store in HTRANS:

C ---------------------------------------------------------
         DO ID = LD_MIN,LD_MAX
C ---------------------------------------------------------

            TR = WK_ARRAY(ID) ! signed transp. into buffer vol.

            IF ( ID .LE. NL_C) THEN
C                Inflow to layer above sill depth:
C           ..... Part of inflow passing buffer volume:
               X = VBUF(K_D,ID) - ( VBUFMX(K_D,ID) - VBUFTR(K_D,ID) )
C                  : part of outer transition volume filled
C                    already with flow in same direction.
               if ( X .GE. 0 ) then
                  if ( X .lt. VBUFTR(K_D,ID) ) then
                      Q_EFF = TR * X/VBUFTR(K_D,ID) !Gradual transition
                  else
                      Q_EFF = TR ! Buffer volume full: flow to basin
                  endif
               else
                  Q_EFF = 0.0 ! Last part of buffer volume empty
               endif
C           ..... Effective mass transport between main volumes:
C                 Inflow reduced by buffer accumulation.
               VBUFDV(K_D,ID) = VBUFDV(K_D,ID) + ABS(TR-Q_EFF)

$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) WRITE ( DEBUG_UNIT,
     &          '(A,I3,2(A,G15.9)/A,I1,A,I3,2(3(A,G15.9:)/))' )
     &        ' WK_ARRAY(', ID, '):', WK_ARRAY(ID), ' TR=', TR
     &      , ' VBUF(',K_D,',',ID
     &      , ')='      ,VBUF(K_D,ID)
     &      , ' VBUFMX=',VBUFMX(K_D,ID)
     &      , ' VBUFTR=',VBUFTR(K_D,ID)
     &      , ' VBUFDV + ',ABS(TR-Q_EFF), '=' ,VBUFDV(K_D,ID)
     &      , ' Q_EFF=', Q_EFF
$endif

               TR =  Q_EFF
            ENDIF

            IF (TR .NE.0 ) THEN
               LLH = LLH+1
               HTR_L(LLH,K_S) = LBUF_SOURCE(IT)
               HTR_L(LLH,K_D) = ID
               HTRANS(LLH) = TR

$if defined DEBUG_TRANH_GT_2
      IF (DEBUG)
     &   WRITE ( DEBUG_UNIT, '( '' >>>: HTRANS('', I3, '')='', G15.7,
     &                2('' HTR_L(_,'',i1, '')='', i3)   )' )
     &               LLH, HTRANS(LLH), (LLH,K,HTR_L(LLH,K),K=1,2)
$endif

            ENDIF


C ---------------------------------------------------------------------
         ENDDO  ! Next layer within destination range for one TR_TEMP
C ---------------------------------------------------------------------


C ------------------------------------------------------------
      ENDDO     ! Next TR_TEMP
C ------------------------------------------------------------


C ================================================================
C Now complete the transport description in HTRANS.
C For each side of the connection, pass through layers, and
C add mixing of water previously flowed in from the other side,
C but not yet effectively mixed into the receiving basin.
C These volumes are noted in VBUF, mixing rates are specified in
C array TCVBUF.  The mixing transfers water from the basin "owning"
C the water in the buffer volume, and reduces the buffer volume
C accordingly. The flow is distributed proportionally over
C direct transports entering the layer, or if no such inflows exist
C set up as a separate transport.
C ================================================================

$if defined DEBUG_TRANH_GT_1
         if (DEBUG) WRITE ( DEBUG_UNIT, '(/A/)')
     &           '######## modify with buffer mixing:',
     &           '         first: Sum transports into each layer'
$endif

C 1. Scan through HTRANS, and sum inflow to each layer on each side
C    above sill depth in WK_ARRAY (basin .
C 2. Scan through buffer volumes, calculate

      DO L = 1, NLA+NLB
         WK_ARRAY(L) = 0.0
      ENDDO

$if defined DEBUG_TRANH_GT_2
         IF (DEBUG .and. ilh.lt.llh )
     &      WRITE ( DEBUG_UNIT, '(1x,2(A5,A18))')
     &          'L', 'HTRANS(L)', 'K', '--> WK_ARRAY(K)'
$endif

      DO L = ILH, LLH
         TR = HTRANS(L)
         if (TR.gt.0) THEN
            K = NLA + HTR_L(L,2)   ! into basin B
         ELSE
            K = HTR_L(L,1)   ! into basin A
         ENDIF
         WK_ARRAY(K) = WK_ARRAY(K) + ABS(TR)

$if defined DEBUG_TRANH_GT_2
         IF (DEBUG) WRITE ( DEBUG_UNIT, '(1x,2(i5,g18.11))')
     &       L, TR, K, WK_ARRAY(K)
$endif

      ENDDO

$if defined DEBUG_TRANH_GT_2
      IF (DEBUG) THEN
         WRITE ( DEBUG_UNIT, '(/5x,A5,A18)' ) 'L', 'WK_ARRAY(L)'
         DO L = 1, NLA+NLB
            write ( debug_unit, '(5x,I5,G18.11)') L, WK_ARRAY(L)
         ENDDO
      ENDIF
$endif

$if defined DEBUG_TRANH_GT_1
      IF (DEBUG) WRITE ( DEBUG_UNIT, '(/A/)')
     &           ' ####### then: Calculate mixing correction:'
$endif

C ------------------------------------------------------------
C Scan through layers, and correct or add transports
C to take care of postponed effect of inflow to buffer volumes:
C ------------------------------------------------------------

      LLH_BASE = LLH


      DO KMD = 1,2                   ! Limit to index of layers
         UPPER_CONNECTED (KMD) = 0   ! in other basin connecting to
      ENDDO                          ! current layer in KMD.

      DO L = 1, NL_C
         DIRECTION = 1
         DO KMD = 1,2
            DIRECTION = - DIRECTION  ! -1 FOR KMD=1, +1 FOR KMD=2
            IF ( VBUF(KMD,L) .LT. 0.0 .or. TCVBUF(KMD).le.0.0 ) then
               Q_MIXING = 0.0
            else
               Q_MIXING = VBUF(KMD,L)/TCVBUF(KMD)
            endif
            IW = (KMD-1)*NLA + L
            W = WK_ARRAY(IW)       ! Transports into layer

            if (W .gt. 0.0) THEN
               ! ------------- Include in existing transports
               WK_ARRAY(IW) = - 1.0 - Q_MIXING/W
                                 ! - Correction factor to be used below,
                                 !    Will amount to Q_Mixing in total.
                                 !    Must always be done.
                                 !    Minus for check below
$if defined DEBUG_TRANH_GT_2
               IF (DEBUG) WRITE( DEBUG_UNIT,
     &              '(''    Corr. factor WK_ARRAY('',I5,'')='',G18.11)')
     &              IW, WK_ARRAY(IW)
$endif

            else  if (Q_mixing.gt.0.0) then
              ! -------------  New transports added

              ! First establish which layers to draw water from:

               LD_MIN = UPPER_CONNECTED(KMD)
               DO WHILE ( LD_MIN .LT. NL_C )
                   IF ( CONNECT (LD_MIN+1,3-KMD) .gt. L-1 ) EXIT
                   LD_MIN = LD_MIN + 1
               ENDDO ! LD_MIN Normally points to first layer -1,
                     ! = NL_C+1 if all layers connected above
               UPPER_CONNECTED(KMD) = LD_MIN

               LD_MAX = LD_MIN
               W_SUM = 0.0 ! will accumulate distribution weights
               DO WHILE ( LD_MAX .LT. NL_C )
                   X = L+1 - CONNECT (LD_MAX+1,3-KMD)
                   IF ( X. le. 0.0 ) EXIT
                   LD_MAX = LD_MAX + 1
                   X = X*(2.0-X)
                   TR_TEMP (LD_MAX) = X ! used as weight below
                   W_SUM = W_SUM + X
               ENDDO ! LD_MAX now normally points to last layer

$if defined DEBUG_TRANH_GT_2
               IF (DEBUG) WRITE( DEBUG_UNIT,
     &              '('' Draw transports from '',2(1X,A,I5),A,G18.11)')
     &                  'LD_MIN+1', LD_MIN+1, 'LD_MAX',LD_MAX,
     &                  ', W_SUM=', W_SUM
$endif

            ! Draw water according to which are most closely connected:
               if (W_SUM.gt.0) then
                  X = DIRECTION*Q_MIXING/W_SUM
                  TR = 0.0   ! Sums transports
                  DO IS = LD_MIN+1, LD_MAX
                     IF ( LLH .ge. MLH ) THEN
                      Write (*,*)'Dimension of HTRANS exceeded in TRANH'
!                      Pause
                     ELSE
                        LLH = LLH +1
                        HTRANS(LLH) = X*TR_TEMP(IS)
                        TR = TR + HTRANS(LLH)
                        HTR_L(LLH,KMD)   = L  ! Destination layer
                        HTR_L(LLH,3-KMD) = IS ! Source layer
                     ENDIF

$if defined DEBUG_TRANH_GT_2
                     IF (DEBUG) WRITE( DEBUG_UNIT,
     &                  '('' Added new HTRANS('', I4, '')='', G15.7,
     &                    2('' HTR_L(_,''i1, '')='', i3)   )' )
     &                    LLH, HTRANS(LLH), (K,HTR_L(LLH,K),K=1,2)
$endif
                  ENDDO
                  Q_MIXING = ABS(TR) ! May be only part of specified
               ELSE
                  Q_MIXING = 0.0
               ENDIF
            ENDIF

            VBUFDV(KMD,L) = VBUFDV(KMD,L) - Q_MIXING
                  ! Transports included in changed rate of VBUF

$if defined DEBUG_TRANH_GT_2
            IF (DEBUG) WRITE(DEBUG_UNIT,
     &           '( 2(A,I3),A,G16.10,A,I5,A,G16.10/3(A,G16.10))' )
     &              ' Buffer (', KMD, ',', L,
     &              ' ): Q_MIXING=', Q_MIXING,
     &              ' WK_ARRAY(',IW,')',W,
     &              ' VBUF=',VBUF(KMD,L),
     &              ' TCVBUF=',TCVBUF(KMD),
     &              ' ->VBUFDV=',VBUFDV(KMD,L)
$endif


         ENDDO
      ENDDO

C -------------------------------------------------------------
C       Scan previous transports accumulated in WK_array,
C       and correct for postponed mixing according to
C       correction factors. All transports will be related to
C       non_zero WK-array, which all have been converted
C       to correction factors above
C -------------------------------------------------------------


$if defined DEBUG_TRANH_GT_2
         IF (DEBUG) WRITE( DEBUG_UNIT,
     &     '('' Corrected for buffer mixing:''/1x,3(A15:A10))')
     &         'Old HTRANS','(number)','corr.factor',
     &         '(number)', '--> New value'
$endif
      DO L = ILH, LLH_BASE
         TR = HTRANS(L)
            ! Set IW to index of WK_ARRAY for destination layer
         if (TR.gt.0) THEN
            IW = HTR_L(L,2)
            if ( iw .gt. NL_C) CYCLE
            IW = NLA + IW
         ELSE
            IW = HTR_L(L,1)
            if ( IW .gt. NL_C) CYCLE
         ENDIF

$if defined DEBUG_TRANH_GT_2
         IF (WK_ARRAY(IW).ge.0.0) then
             WRITE( *,*) 'Error in TRANSP_H postponed mixing:',
     &             'WK_ARRAY(IW)=',WK_ARRAY(IW),' should be <0'
         ENDIF
$endif

         HTRANS(L) = - WK_ARRAY(IW)*TR

$if defined DEBUG_TRANH_GT_2
         IF (DEBUG) WRITE( DEBUG_UNIT, '(1x,3(G15.7:I10))')
     &        TR, L, WK_ARRAY(IW),IW, HTRANS(L)
$endif
      ENDDO


C ===============================================================
C Update buoyancy flux to top or bottom of mixed layer:
C ===============================================================


      DO L = ILH, LLH
C     ....... Direction indexes:
         IF( HTRANS(L) .gt. 0.0 ) then
            K_S = 1
            K_D = 2
         ELSE
            K_S = 2
            K_D = 1
         ENDIF
         IS = HTR_L(L,K_S)
         ID = HTR_L(L,K_D)
C     ....... Skip if destination below well-mixed layer,
C                Terminate loop when below both well_mixed layers
         IF ( ID .GT. LMIX(K_D) ) THEN
            IF ( IS .GT. LMIX(K_S) )  THEN
               EXIT
            ELSE
               CYCLE
            ENDIF
         ELSE
C     ....... Transport into well_mixed layer:
C                Update buoyancy flux to surface,
C                or move to layer below well-mixed layer:
            IF ( K_S .EQ.1 ) THEN
               BF =  ABS(HTRANS(L))*(DENSA(IS)-DENSB(ID))
            ELSE
               BF =  ABS(HTRANS(L))*(DENSB(IS)-DENSA(ID))
            ENDIF      ! BF is given in unit kg/days
            IF ( BF .LT. 0.0 ) THEN
               BSF(K_D) = BSF(K_D) - BF / sec_per_day
                       ! accumulated as kg/s
            ELSE
               HTR_L(L,K_D) = MIN( LBOTT(K_D), MAX( ID, LMIX(K_D)+1 ) )
            ENDIF
         ENDIF
      ENDDO



$if defined DEBUG_TRANH
      IF (DEBUG) THEN
        DO L = ILH, LLH
          WRITE ( DEBUG_UNIT, '('' L='',I3,A,E16.10,2(A,I3))' )  L
     &                , '  HTRANS(L) =', HTRANS (L)
     &                , '  HTR_L(L,1)=', HTR_L (L,1)
     &                , '  HTR_L(L,2)=', HTR_L (L,2)
        ENDDO
      ENDIF
$endif


C ----------------- Restore updated buoyancy flux values:
      BSFLUX_A =  BSF(1)
      BSFLUX_B =  BSF(2)

      END SUBROUTINE



C  ===========================================================
C  Set up depth index for connections between basins,
C  according to density, called from TRANH.

      SUBROUTINE TRANH2 ( DEBUG, NL_C, CONNECT, NLB, DENSB)
      
      
      LOGICAL DEBUG
      INTEGER NL_C  ! Number of connected layers:
      real*8 CONNECT(NL_C)
      INTEGER NLB
      real*8 DENSB(NLB)

C ---------- local variables -------------
      INTEGER LDB, IL
      real*8 DA, DD_AB, DD_B

$if defined DEBUG_TRANH_GT_2
      real*8 X
!!      INCLUDE 'DEBUG.INC'
$endif


C ------------ EXECUTE ------------------
      LDB = 1
      DA = CONNECT(1)

$if defined DEBUG_TRANH_GT_2
      if (debug) WRITE (debug_unit, '(2A6,2A14, A)' )
     &    ' il', ' ldb', ' dd_AB', ' dd_B',
     &    ' CONNECT(IL): density -->index'
$endif


      DO IL = 1,NL_C

$if defined DEBUG_TRANH_GT_2
         X = CONNECT(IL)
$endif

         DA = MAX( DA, CONNECT (IL) )
C     ...... Find depth index in basin B for uppermost layer
C            with density >= density of layer IL in basin A,
C            and a density gradient to the layer below,
C            (i.e. at bottom layer of any homogeneous region in B
C             consisting of more than layer)
         DO WHILE ( LDB.LT.NLB)
             DD_AB = DA - DENSB(LDB)
             if ( DD_AB .le. 0 .and.  DA .lt. DENSB(LDB+1) ) EXIT
             LDB = LDB +1
         ENDDO
         CONNECT(IL) = LDB

         if ( LDB .gt. 1 .and. DD_AB.lt.0 ) THEN
             DD_B = DENSB(LDB) - DENSB(LDB-1)
             if (DD_B .gt. 0.0 ) THEN
                 CONNECT(IL) = FLOAT(LDB)+DD_AB/DD_B

             endif
         endif
C  (Note that CONNECT may point to layer below sill)

$if defined DEBUG_TRANH_GT_2
         IF (DEBUG) WRITE ( DEBUG_UNIT, '(2I6,4G14.7)' )
     &       il, ldb, dd_AB, dd_B, X, CONNECT(IL)
$endif

      ENDDO

      END SUBROUTINE



C ====================================================================
      SUBROUTINE MTRAN1( DEBUG, NLH_DIM, NLH, HTRANS, LAYER_A, LAYER_B,
     &                   NLA, NLB, CONSA, CONSB, HTRMIX,
     &                   NMD_A, NMD_B, MASS_DERIV_A, MASS_DERIV_B,
     &                   IMPORT_A, IMPORT_B )
      

C Update derivative terms with transport across connection from A to B:
C External parameters:

      LOGICAL DEBUG

      INTEGER NLH_DIM, NLH
C           - Number of flows across connection
      real*8  HTRANS (NLH_DIM)
C           - Water transports accross connection:
C                 >0:  A -> B  ,  <0: B->A
      INTEGER LAYER_A(NLH_DIM), LAYER_B(NLH_DIM)
C           - Index for participating layers in the two basins.
      INTEGER NLA, NLB
C          - Number of possible involved layers in basin A and B.
      real*8 CONSA( NLA ), CONSB( NLB )
C           - Concentrations in outflowing layers.
C             Because of residual effects of buffer volume mixing,
C             inflow from layers below sill depth
C             in internal basins is possible and allowed.
      real*4 HTRMIX
            ! Degree of mixing between continguous transports in
            ! the same direction ( intended range 0-1 )

      INTEGER NMD_A, NMD_B
C          - Number of layers in basin A and B with time integration
C            of state variables  (=0 for external basins)
C
      real*8 MASS_DERIV_A (NMD_A), MASS_DERIV_B (*)
C           - Time derivatives in basins (only internal)
      real*8 IMPORT_A, IMPORT_B
C           - Accumulated imports to basins


C ------------------ local variables: ---------------------
      INTEGER L, L1, Source, IA, IB
      real*8 Q, MASS_TRANSP, SUM_Q, SUM_M_TR, C

C ---------------------------------------------------------

$if defined DEBUG_MTRAN1
!!      INCLUDE 'DEBUG.INC'
      character*72 OUTPUT_LINE

      IF (DEBUG) THEN

          write( debug_unit,*) '------------- MTRAN1:'
          WRITE( debug_unit,'(1x,5A6/1x,5I6)')
     &         'NLA', 'NLB', 'NMD_A', 'NMD_B', 'NLH',
     &          NLA, NLB, NMD_A, NMD_B, NLH
          WRITE( debug_unit,'(1x,4A18)')
     &        'CONSA', 'MASS_DERIV_A', 'CONSB', 'MASS_DERIV_B'
          L1 = MAX( NLA, NLB, NMD_A, NMD_B)
          do L = 1, L1
             Output_line = ' '
             if (L .le. NLA) then
                 write ( OUTPUT_LINE(2:18),'(G17.11)') CONSA(L)
             endif
             if (L .le. NMD_A) then
                 write ( OUTPUT_LINE(20:36),'(G17.11)') MASS_DERIV_A(L)
             endif
             if (L .le. NLB) then
                 write ( OUTPUT_LINE(38:54),'(G17.11)') CONSB(L)
             endif
             if (L .le. NMD_B) then
                 write ( OUTPUT_LINE(56:72),'(G17.11)') MASS_DERIV_B(L)
             endif
             write( debug_unit, '(1x,I4,A)') L, Output_line
          enddo

          write( debug_unit,'('' Import _A:'',G17.11, '' _B:'',G17.11)')
     &           IMPORT_A, IMPORT_B

          WRITE( debug_unit,'(4X,3A6,A18)' )
     &           'L', 'IA', 'IB', 'HTRANS(L)'
          WRITE( debug_unit,'(1X,3I6,G18.11)' )
     &           (L, LAYER_A(L), LAYER_B(L), HTRANS(L), L=1,NLH )

      ENDIF
$endif



      L = 1

      DO WHILE ( L .le. NLH)


         L1 = L
         Q = HTRANS(L)


         ! ---- Process any transports A --> B by handling effects in A,
         !      and accumulate sums for handling effects in B:
         do while ( Q .gt. 0.0 )
            IA = LAYER_A(L)
            MASS_TRANSP = Q*CONSA(IA) ! Transports >0 from A to B
            MASS_DERIV_A(IA) = MASS_DERIV_A(IA) - MASS_TRANSP
            IMPORT_A = IMPORT_A - MASS_TRANSP
            if ( NMD_B .gt. 0 ) then
               if (L.eq.L1) then
                  SUM_M_TR = MASS_TRANSP
                  SUM_Q = Q
               else
                  SUM_M_TR = SUM_M_TR + MASS_TRANSP
                  SUM_Q = SUM_Q + Q
               endif
            endif


$if defined DEBUG_MTRAN1
            if (DEBUG) THEN
               write (debug_unit,
     &            '('' Tr. A-->B:'',3(1x,A,I3),3(2(A,G16.10:)/))')
     &             'L1', L1, 'L', L, 'IA', IA,
     &             ' Q', Q, ' MASS_TRANSP', MASS_TRANSP,
     &             ' SUM_M_TR', SUM_M_TR, ' SUM_Q', SUM_Q,
     &             ' MASS_DERIV_A', MASS_DERIV_A(IA),
     &             ' IMPORT_A', IMPORT_A
            ENDIF
$endif


            ! ----- Check whether to homogenize with next transport:
            L = L+1
            if ( L .le. NLH ) THEN
               Q = HTRANS(L)
               if ( NMD_B .le. 0 ) EXIT
               Source = LAYER_A(L)
               if ( Source .gt. IA+1 .or. Source .lt. IA ) EXIT
            else
               EXIT
            endif
         ENDDO   ! on exit from loop, L=L1+1 if single transport


         ! ------- Handle effects on B of processed flows A --> B
         if (L .gt. L1) then

$if defined DEBUG_MTRAN1
         if (DEBUG) THEN
            write ( debug_unit,'(3(1x,A,I5),1x,A,G18.11)')
     &           'On exit: L1=', L1, 'L=', L, 'Source=', Source, 'Q=',Q
         endif
$endif

            IF ( L .gt. L1+1 ) then
            ! ------- transports L1 to L-1 to be partly homogenized.
            !         this will only happen if NMD_B>0, i.e.
            !         if B is internal basin, see above.
               C = (SUM_M_TR/SUM_Q)*HTRMIX

$if defined DEBUG_MTRAN1
         if (DEBUG) THEN
            write ( debug_unit,'(1x,2(A,G15.7)/1X,3A6,3A18)')
     &          ' HTRMIX:', HTRMIX, ' C:',C,
     &            'L1', 'IA', 'IB',
     &            'MASS_TRANSP', 'MASS_DERIV', 'IMPORT_B'
         endif
$endif

               do while ( L1 .lt. L )
                  IB = LAYER_B(L1)
                  IA = LAYER_A(L1)
                  MASS_TRANSP = HTRANS(L1)*(CONSA(IA)*(1.-HTRMIX)+C)
                  MASS_DERIV_B(IB) = MASS_DERIV_B(IB) + MASS_TRANSP
                  IMPORT_B = IMPORT_B + MASS_TRANSP

$if defined DEBUG_MTRAN1
                 if (DEBUG) THEN
                    write ( debug_unit,'(1x,3I6,3G18.11)')
     &                  L1, IA, IB, MASS_TRANSP,
     &                  MASS_DERIV_B(IB), IMPORT_B
                 ENDIF
$endif

                  L1 = L1+1
               ENDDO
            else

            ! ------- single transport L1: Use MASS_TRANSP from above
               IF (NMD_B .GT.0 ) THEN
                  IB = LAYER_B(L1)
                  MASS_DERIV_B(IB) = MASS_DERIV_B(IB) + MASS_TRANSP
                  IMPORT_B = IMPORT_B + MASS_TRANSP
$if defined DEBUG_MTRAN1
                  if (DEBUG) THEN
                     write ( debug_unit,'(1x,2A6,3A18)') 'L1', 'IB',
     &                  'MASS_TRANSP', 'MASS_DERIV', 'IMPORT_B'
                     write ( debug_unit,'(1x,2I6,3G18.11)') L1,  IB,
     &                   MASS_TRANSP, MASS_DERIV_B(IB), IMPORT_B
                  endif
$endif
               ENDIF

               L1 = L
            endif
         endif

         ! ---- Process any transports A <-- B by handling effects in B
         !      if B is internal basin,
         !      and accumulate sums for later handling of effects in A.
         do while ( Q .lt. 0.0 )
            IB = LAYER_B(L)
            MASS_TRANSP = -Q*CONSB(IB)  ! >0 : TRANSPORT A<---B
            IF ( NMD_B .GE. IB ) THEN
               MASS_DERIV_B(IB) = MASS_DERIV_B(IB) - MASS_TRANSP
               IMPORT_B = IMPORT_B - MASS_TRANSP
            ENDIF

            if (L.eq.L1) then
               SUM_M_TR = MASS_TRANSP
               SUM_Q = -Q
            else
               SUM_M_TR = SUM_M_TR + MASS_TRANSP
               SUM_Q = SUM_Q - Q
            endif

$if defined DEBUG_MTRAN1
            if (DEBUG) THEN
               write (debug_unit,
     &            '('' Tr. A<--B:'',3(1x,A,I3),2(2(A,G16.10:)/))')
     &             'L1', L1, 'L', L, 'IB', IB,
     &             ' Q', Q, ' MASS_TRANSP', MASS_TRANSP,
     &             ' SUM_M_TR', SUM_M_TR, ' SUM_Q', SUM_Q
               IF ( NMD_B .GE. IB )
     &            write (debug_unit,
     &            '( 2(A,G16.10))')
     &             ' MASS_DERIV_B', MASS_DERIV_B(IB),
     &             ' IMPORT_B', IMPORT_B
            ENDIF
$endif

            ! ----- Check whether to homogenize with next transport:
            L = L+1
            if (L .le. NLH) then
                Q = HTRANS(L)
                Source = LAYER_B(L)
                if ( Source .gt. IB+1 .or. Source .lt. IB ) EXIT
            else
                EXIT
            endif
         ENDDO   ! on exit from loop, L=L1+1 if single transport

$if defined DEBUG_MTRAN1
         if (DEBUG) THEN
            write ( debug_unit,'(3(1x,A,I5),1x,A,G18.11)')
     &           'On exit: L1=', L1, 'L=', L, 'Source=', Source, 'Q=',Q
         endif
$endif


         ! ------- Handle effects on A of processed flows A <-- B
         if (L .gt. L1) then
            IF ( L .gt. L1+1 ) then
            ! ------- transports L1 to L-1 to be partly homogenized.
            !         this will only happen if NMD_B>0, i.e.
            !         if B is internal basin, see above.
               C = (SUM_M_TR/SUM_Q)*HTRMIX

$if defined DEBUG_MTRAN1
         if (DEBUG) THEN
            write ( debug_unit,'(1x,2(A,G15.7)/1X,3A6,3A18)')
     &          ' HTRMIX:', HTRMIX, ' C:',C,
     &            'L1', 'IA', 'IB',
     &            'MASS_TRANSP', 'MASS_DERIV_A', 'IMPORT_A'
         endif
$endif

               do while ( L1 .lt. L )
                  IB = LAYER_B(L1)
                  IA = LAYER_A(L1)
                  MASS_TRANSP =-HTRANS(L1)*(CONSB(IB)*(1.-HTRMIX)+C)
                  MASS_DERIV_A(IA) = MASS_DERIV_A(IA) + MASS_TRANSP
                  IMPORT_A = IMPORT_A + MASS_TRANSP

$if defined DEBUG_MTRAN1
                 if (DEBUG) THEN
                    write ( debug_unit,'(1x,3I6,3G18.11)')
     &                  L1, IA, IB, MASS_TRANSP,
     &                  MASS_DERIV_A(IA), IMPORT_A
                 ENDIF
$endif

                  L1 = L1+1
               ENDDO
            else
            ! ------- single transport L1: Use MASS_TRANSP from above
               IA = LAYER_A(L1)
               MASS_DERIV_A(IA) = MASS_DERIV_A(IA) + MASS_TRANSP
               IMPORT_A = IMPORT_A + MASS_TRANSP

$if defined DEBUG_MTRAN1
               if (DEBUG) THEN
                  write ( debug_unit,'(1x,2A6,3A18)')
     &            'L1','IA','MASS_TRANSP','MASS_DERIV_A','IMPORT_A'
                  write ( debug_unit,'(1x,2I6,3G18.11)')
     &              L1,  IA, MASS_TRANSP, MASS_DERIV_A(IA), IMPORT_A
               endif
$endif

               L1 = L
            endif
         endif

      END DO

C Check on indexes above serves to prohibit violation of index range,
C as a security measure, and to prohibit updating of inflows to
C external basins, for which NMD_B should be set to 0.

      END SUBROUTINE


      end Module