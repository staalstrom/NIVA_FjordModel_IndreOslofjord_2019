      Module fx_Transp_v

      use fx_RunControl, only: DEBUG_UNIT
      use ModelDimensions, only: dimMBI

      implicit none
      
C >>>>> moved from old main program:
      real*8 MIXCONST(dimMBI)
C <<<<<



C ===================================================================
C Eutrophication model   - File:    TRANSP_V.FOR
C                                   Birger Bjerkeng, NIVA.
$undefine vmix_version_distr
			!    undefined: adds mixing energy to receiving layer only
			!      defined: distributes mixing energy over layers according to depth,
			!               from below inflow layer down to the receiving layer.
			!      NOTE! May also be redistributed according to density.
			!            Does not use volume, since it is effect of localized inflows


$undefine diff_comp
         ! undefined: never negative diffusion
         !   defined: allows negative diffusion to compensate for advection
         !            (extent of compensation defined in O_5, see below)

$undefine kinetic_energy
         ! determines if kinetic energy of inflows is included as mixing energy
         !   (to receiving layer only)

C Contains subroutines:
C     TRANC  preliminary sums up effects of in/outflows within basin
C     TRANV1 sets up arrays describing vertical mixing intensity,
C            and establishes limit on time step.
C                 - both called from TRANSP_SETUP in TRANSP_2.FOR
C     TRANV2 sets up final description of vertical mixing,
C            scaled for the actual timestep taken.
C                 - called from entry point TRNADJ in TRANSP_1.FOR
C No calls are made to other Modules.


!------------- compile options:
!              turns on/off DEBUG-DUMP TO UNIT DEBUG_UNIT
$undefine DEBUG_TRANC

$undefine DEBUG_TRANV1
$undefine DEBUG_TRANV1_SFMIX
$undefine DEBUG_TRANV1_DETAIL
$undefine DEBUG_TRANV1_TSTEP_LIM

$undefine DEBUG_TRANV2
$undefine DEBUG_TRANV2_VTRNEG

$undefine DEBUG_TRCF 


!------------- For consistency of code inclusion:
$if defined DEBUG_TRANV1
$define SOME_DEBUG_OF_TRANV1
$endif

$if defined DEBUG_TRANV1_SFMIX
$define SOME_DEBUG_OF_TRANV1
$endif

$if defined DEBUG_TRANV1_DETAIL
$define SOME_DEBUG_OF_TRANV1
$endif

$if defined DEBUG_TRANV1_TSTEP_LIM
$define SOME_DEBUG_OF_TRANV1
$endif

      Contains

      
C ==================================================================

      SUBROUTINE TRANC( DEBUG, DEPTH, UFLOW, NL_C, MLCONN, CONNECT,
     &                  NLH_DIM, NLH, HTRANS, HTRL_A, HTRL_B,
     &                  NLA, DENS_A, VTRANS_A, VMIX_A,
     &                  NLB, DENS_B, VTRANS_B, VMIX_B )
      

C Accumulates transports to/from other basins for each layer
C in internal basins, stored temporarily in VTRANS_A(B) (K,L):
C    K=1: Net influx from other basins (sum of signed values)
C    K=2: 2*throughflow (sum of absolute values)
C
C The accumulated values will later be used by TRANV1 to set up
C final transport description, and to evaluate permissible time step
C for integration, to ensure stable solution.

C This subroutine adds to current values VTRANS, initiated in HTRANS,
C and may thus be called more than once for each basin.
C
C The subroutine may be called twice for each connection;
C once for each internal basin connected.

C ------------------- formal arguments: -------------------
      LOGICAL DEBUG

      real*8 DEPTH(*) ! Layer depth limits
      real*8 UFLOW(*) ! velocity in flows - gives kinetic energy
C           local index of UFLOW: 1=surface layer in connection

      INTEGER NL_C, NLH_DIM, NLH
C         (I) Number of connected layers and num. of horizontal transports
      INTEGER MLCONN
      real*8  CONNECT(MLCONN,2)
C         (I) Density defined connections, see TRANSP_H.TRANH2
      real*8  HTRANS(NLH_DIM)            ! m3/day
C         (I) Flows between basins: >0 from A to B, <0: from B to A.
      INTEGER HTRL_A(NLH_DIM), HTRL_B(NLH_DIM)
C         (I) Local index (surf.=1) of connected layers in basin A and B
      INTEGER NLA, NLB
C         (I) Number of layers in basin A and B
      real*8 DENS_A(NLA), DENS_B(*)
C         (I) Density, DENS_B is an array and can be used
C             for occurring values of IB even if NLB is zero
      real*8 VTRANS_A(2,NLA), VTRANS_B(2,*)
C        (I/O) (1,L) : Updated net inflow into layer
C                        before adjusting for changing volume.
      real*8 VMIX_A(2,NLA), VMIX_B(2,*)
C              (1,L) : Input: Depth integral of density
C              (2,L) : I/O: Net released potential energy
C                           connected to advective flows between basins

C ---------- LOCAL VARIABLES:
      INTEGER L,IA,IB, I_MIN, IC, I
      real*8 Grav /9.81 /  ! gravity constant, m/s2
      real*8 Q
      real*8 DD, X

$if defined DEBUG_TRANC
!      INCLUDE 'DEBUG.INC'
      IF (DEBUG) THEN
         write(DEBUG_UNIT,*) '------------------- tranc -------'
         WRITE(DEBUG_UNIT,*)  'NLH, NLA, NLB, NL_C, MLCONN:',
     &                         NLH, NLA, NLB, NL_C, MLCONN
      ENDIF
$endif

C ---------------------------------------------------------
C  VTRANS(1,..) accumulates net transport into layer,
C               used when adjusting to volume conservation in TRANV1.
C  VTRANS(2,..) accumulates outgoing transports,
C               used to limit time-step in TRANV1.
C ---------------------------------------------------------

      DO L = 1, NLH       ! note: transport index, not index of layer 

         IA = HTRL_A(L)
         IB = HTRL_B(L)
c            USED AS INDEX TO ACCUMULATING MIXING ENERGY, WILL
c            BE REDISTRIBUTED IN TRANV1 BELOW

         Q = HTRANS(L)

$if defined DEBUG_TRANC
      IF (DEBUG) THEN
         write( DEBUG_UNIT,'(3(3X,A,I3),3X,A,G18.11)')
     &           'L',L,'IA',IA,'IB',IB,'Q',Q
      ENDIF
$endif

C .......... effects in basin A.
         VTRANS_A(1,IA) = VTRANS_A(1,IA) - Q

         IF ( Q .gt. 0.0D0 ) THEN ! accumulate outflows from A:
            VTRANS_A(2,IA) = VTRANS_A(2,IA) + Q

         ELSE    ! inflow to A from layer IB in B

$IF !DEFINED VMIX_VERSION_DISTR
            I_MIN = IA
$else
            I_MIN = MIN (IA, IB+1)
$endif

$if defined DEBUG_TRANC
            IF (DEBUG) write( DEBUG_UNIT,'(1x,A,3I6)')
     &                        'IA, IB, I_MIN=',IA, IB, I_MIN
            if ( IB.lt.1 ) then
                write(debug_unit,'(1x,A/4I5)')
     &           'mlconn, L, IA, IB',
     &            mlconn, L, IA, IB
                STOP 'feil i tranc, see DEBUG_UNIT'
            endif
$endif

            if ( IB .le. nl_c) then
               IC = CONNECT(IB,2) !layer in A connected to IB in basin B
$if defined DEBUG_TRANC
               if ( IC .Lt. 1 .or. IC .gt. NLA ) then
                   write( DEBUG_UNIT,*)
     &                  'l',L,'ib',ib,'CONNECT(IB,2)=',connect(ib,2)
                   STOP 'feil i tranc, see DEBUG_UNIT'
               endif
$endif

               IF ( IC .gt. NL_C ) THEN ! Inflow to A descending by density:
C               Accumulate potential energy dissipated by inflow from
C               density difference between inflow and recipient
C               integrated over depth interval given by density connection:
                  X = DENS_A(IC)-DENS_A(IC-1)
                  IF (X.GT.0.0) THEN
                     X = (DENS_B(IB)-DENS_A(IC-1))/X
                  ELSE
                     X = 1.0
                  ENDIF
C                  X = depth interpolation measure
C                      between layers IC and IC-1 in basin A.
                  DD = ( X*( DEPTH(IC+1)-DEPTH(IC) )
     &                  +DEPTH(IC)-DEPTH(IB+1)     ) * DENS_B(IB)
     &                -( X*( VMIX_A(1,IC)-VMIX_A(1,IC-1) )
     &                  +VMIX_A(1,IC-1)-VMIX_A(1,IB)     ) !unit kg/m2

                  IF ( DD .gt. 0.0 ) then  ! (adds , since Q<0)
                     DD = DD/(IA-I_MIN+1)
                     do I = I_MIN, IA, 1
                        VMIX_A(2,I) = VMIX_A(2,I) - grav*DD*Q
                     enddo
                  ENDIF    ! Effect = Newton*m/day = kg*m2/s2/day

$if defined DEBUG_TRANC
      IF (DEBUG) THEN
         write( DEBUG_UNIT,'(1x,2(A,I3)/3(1X,A,G17.11))')
     &         'DESCENDING FLOW INTO A: I_MIN=',I_MIN,' IC=',IC,
     &           'X',X,'DD',DD, 'VMIX_A(2,I_MIN)', VMIX_A(2,I_MIN)
      ENDIF
$endif

               ENDIF ! Could subtract negative contributions ?
C                    Would then have to handle VMIX_A(2,..)<0
           ENDIF

$if defined KINETIC_ENERGY
C        ...  Kinetic horizontal energy in inflow is also included:
            VMIX_A(2,IA) =
     &            VMIX_A(2,IA) +ABS(1000.*(UFLOW(IB)**2)*Q/2.0)
                                ! kg/m3  *  m2/s2       *m3/day    
$endif
         ENDIF

         IF (NLB.le.0) CYCLE
            !  otherwise B is also an internal basin, 
            ! so calculate effects in basin B.
         VTRANS_B(1,IB) = VTRANS_B(1,IB) + Q
         
         IF ( Q .lt. 0.0D0 ) THEN ! accumulate outflows from B:
            VTRANS_A(2,IB) = VTRANS_A(2,IB) - Q
         ELSE

$IF !DEFINED VMIX_VERSION_DISTR
            I_MIN = IB
$else
            I_MIN = MIN (IA+1, IB)
$endif

$if defined DEBUG_TRANC
            IF (DEBUG) write( DEBUG_UNIT,'(1x,A,3I6)')
     &                        'IA, IB, I_MIN=',IA, IB, I_MIN
            if ( IA.lt.1 ) then
                write(debug_unit,'(1x,A/4I5)')
     &           'mlconn, L, IA, IB',
     &            mlconn, L, IA, IB
                STOP 'feil i tranc, see DEBUG_UNIT'
            endif
$endif

            if ( IA .le. nl_c) then
             IC = CONNECT(IA,1) !layer in B connected to IA in basin A
$if defined DEBUG_TRANC
             if ( IC .Lt. 1 .or. IC .gt. NLB ) then
                write( DEBUG_UNIT,*)
     &               'l',L,'iA',iA,'CONNECT(IA,1)=',connect(iA,1)
                STOP 'feil i tranc, see DEBUG_UNIT'
            endif
$endif

             IF ( IC .gt. NL_C ) THEN !Inflow to B descending by density:
C               Accumulate potential energy dissipated by inflow from
C               density difference between inflow and recipient
C               integrated over depth interval given by density connection:
               X = DENS_B(IC)-DENS_B(IC-1)
               IF (X.GT.0.0) THEN
                  X = (DENS_A(IA)-DENS_B(IC-1))/X
               ELSE
                  X = 1.0
               ENDIF
C                  X = depth interpolation measure
C                      between layers IC and IC-1 in basin B.
               DD =  ( X*( DEPTH(IC+1)-DEPTH(IC) )
     &                 + DEPTH(IC)-DEPTH(IA+1)      ) * DENS_A(IA)
     &              -( X*( VMIX_B(1,IC)-VMIX_B(1,IC-1) )
     &                 +VMIX_B(1,IC-1)-VMIX_B(1,IA)      ) !unit kg/m2

                IF ( DD .gt. 0.0 ) then
                  DD = DD/(IB-I_MIN+1)
                  do I = I_MIN, IB, 1
                     VMIX_B(2,I) = VMIX_B(2,I) + grav*DD*Q
                  enddo
C                    Effect: (m/s2)*[(kg/m3)*m]*m3/day = Newton*m/day
                ENDIF

$if defined DEBUG_TRANC
      IF (DEBUG) THEN
         write( DEBUG_UNIT,'(1x,2(A,I3)/3(1X,A,G17.11))')
     &         'DESCENDING FLOW INTO B: I_MIN=',I_MIN,' IC=',IC,
     &           'X',X,'DD',DD, 'VMIX_B(2,I_MIN)', VMIX_B(2,I_MIN)
      ENDIF
$endif

              ENDIF
             ENDIF

$if defined KINETIC_ENERGY
C        Kinetic horizontal energy in inflow is also included:
            VMIX_B(2,IB) = VMIX_B(2,IB)
     &                       + ABS(1000.*(UFLOW(IA)**2)*Q/2.0)
$endif
         ENDIF
      END DO


      END Subroutine

C ==================================================================
      SUBROUTINE TRANV1( DEBUG, DEPTH, MB, NB,
     &   MIXFAC, MIXEXP, N2LIM,
     &   SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
     &   ZBOTMI, INDXI, NL, AREA, VLAYER, VLCORF, NLVOPN,
     &   DENS, VTRANS, VMIX,
     &   MW, WK_ARRAY, BWFREQ, ZMID,
     &   TSTEP_ADV )
      

C         Initiate vertical mixing in internal basins:
C         adjusted to integration time by TRANV2 below.
C ==================================================================


C ----------------------- External arguments: ------------------------
      LOGICAL DEBUG

      real*8 DEPTH(*)
C        (I) Depth for layer limits (common for all basins)
      INTEGER MB, NB
C        (I) Number of basins
      real*4 MIXFAC
      real*4 MIXEXP, N2LIM
C        (I) Mixing parameters controlling tidal mixing
      real*4 SFMIXC(NB), SFMIXZ(MB,2)
C         Surface mixing constant (as energy, m2/s3)
C         and depth scale for dampening of surface mixing,
C                   EXP(-(Z-SCALE1)/SCALE2)

      real*4 GMIXFR (NB)
C        (I) Fraction of released gravitation energy from inflows
C            used to work against gravitational field in vertical mixing
      real*4 GMIXDC
C        (I) vertical specific reduction rate of gravitational energy
C            at stability BW_FREQ = 1.0 (per meter)
      real*4 GMIXDX
C        (I) Stability dependence in vertical reduction:
C            exponent for BW_FREQ
      real*8 ZBOTMI(NB)
C        (I) Bottom depths in each basin
      INTEGER INDXI(NB+1)
C        (I) Layer index limits for each basin (NB+1 values)
      INTEGER NL
      real*8 AREA(NL)
C        (I) Horizontal area at top of layers
      real*8 VLAYER(NL)
C        (I) Volume of layers
      real*8 VLCORF(NB)
C        (I) Correction factor for volume open across som connection
      INTEGER NLVOPN(NB)
C        (I) Number of layers open across some connection
      real*8 DENS(NL)
C        (I) Density (sigma-t) in each layer
      real*8 VTRANS(2,NL)  !m3/day
C        (I): (1,L): Vertical flux up from layer L to layer L-1,
C                    through adjusted limit (constant VLAYER),
C                    (Values from TRANC, summed from bottom up,
C                     and adjusted for variation in total volume)
C             (2,L) Sum of outflows from layer L.
C                   (used to limit TSTEP_ADV)
C        (O): Prelim. description of vertical transports between layers
C             (1,L): Unchanged.
C             (2,L): Diffusive volume exchange rate between layers.
      real*8 VMIX(2,NL)
C        (I): (1,L): not used, available for temporary storage.
C             (2,L): Potential (+ kinetic) energy involved in
C                    horisontal inflows into layer L (Nm/day)
C                        ( set up by subroutine TRANC)
C        (O): Work done against gravity field by vertical mixing:
C             (1,L): Total mixing energy doing work against gravity
C             (2,L): From inflows (redistributed downwards by TRANV1)

      INTEGER MW ! >= 2*number of layers in each basin. (TRANSP_1.FOR)
      real*8  WK_ARRAY     (0:MW)
C            - Scratch array for storing values locally.

      real*8 BWFREQ (NL) ! Stability, as Brunt_w„isela frequency

      real*8 ZMID   (NL) ! Mean depth of layers (volume weighted)
C              NB! Note ZMID values used by PHYT_ZOO.

      real*8 TSTEP_ADV
C       (I/O): max. limit for next time-step from advective transports


C ------------------------ LOCAL VARIABLES ------------------------
      INTEGER IB, LVBASE, LAYERS, LL, LG
      real*8 ZBOT, DZ, NSQUARED, BW_FREQ, DFCOEF, EFFECTIVE_AREA
      real*8 ADV_TRANSP, DIFF_TRANSP, OUT_FLUX
      real*8 ZGRAV, ZGRAV_BELOW, AREA_BELOW, DIFF_BELOW
      real*8 SUM1, SUM2, E_ON_F, P, X, W, F, SFX, SFX_ABOVE
      real*8 TSTEP_NEW

$if defined DEBUG_TRANV1_TSTEP_LIM
      LOGICAL ADV_STEP_DECREASED
      INTEGER LIMIT_ADV
      real*8 LIMIT_STORE(4)
$endif

      real*8 grav       /9.81  / ! gravitational accel. m/s2
      real*8 rho_Approx /1000.0/ ! Approx. density in kg/m3

      real*8 sec_per_day
      parameter ( sec_per_day = 24.*3600. )

C     real*8 MAX_FRACTION_REPLACED /0.7/ !Mossesundet
      real*8 MAX_FRACTION_REPLACED /0.3/ !Indre Oslofjord


C ------------------------------------------------------------------

$if defined SOME_DEBUG_OF_TRANV1
!      INCLUDE 'DEBUG.INC'
      IF (DEBUG) WRITE(DEBUG_UNIT,*) '----------- TRANV1 -----------'
$endif

$if defined DEBUG_TRANV1
      IF (DEBUG) WRITE(DEBUG_UNIT,'(2(3X,A,G16.8))')
     &   'GMIXDC=',GMIXDC,'GMIXDX=',GMIXDX
$endif

C ----------------------------------------------
      DO IB = 1,NB  !  All internal basins:
C ----------------------------------------------

C   ...... layer index limits:
         LVBASE = INDXI(IB)

         LAYERS = INDXI(IB+1)-LVBASE

$if defined DEBUG_TRANV1
         if (MW.lt.2*LAYERS) Then
             WRITE(*,'(1x,2(A,I5),A)' )
     &          ' Error in TRANSP_V.FOR/TRANV1: MW=',MW,
     &         '< 2*layers=', 2*LAYERS, ': too small work array'
!             Pause
         endif
         IF (DEBUG)
     &       WRITE( DEBUG_UNIT,'(4(3X,A,I5))') 'IB',IB,'LVBASE',LVBASE,
     &             'LAYERS',LAYERS
$endif


C ===================================================================
C Density stratification will influence both vertical distribution
C of mixing energy and conversion of energy into mixing coefficients.
C This code sequence will set up the necessary description:
C ===================================================================


C -------------------------------------------------------------------
C Calculate stratification parameters and mixing reduction factor p:
C -------------------------------------------------------------------

C ............ Initiate for bottom layer:
C     ........ Effective bottom depth, ensure that it is
C              at least 50% of full layer thickness below top of layer,
C              to avoid unreasonable mixing
         ZBOT = MAX(ZBOTMI(IB),(DEPTH(LAYERS)+DEPTH(LAYERS+1))/2.0)
C     ........ Area at lower limit of layer:
         AREA_BELOW = 0.0

$if defined DEBUG_TRANV1_DETAIL
      IF (DEBUG) THEN
          write(DEBUG_UNIT,*)'Depth:',(Depth(LL),LL=1,Layers+1)
      ENDIF
$endif

         X = 1.0

C ............Process layer:
         DO LL = LAYERS, 1, -1   ! Local index, used for depth
            LG = LL + LVBASE ! global layer index

C        ..... Gravity center of layer LL (globally LG),
C               (assumes linearly increasing area):
            X = ( 1.0 + 1.0/(1.0+AREA_BELOW/AREA(LG)) ) /3.0
            ZGRAV = ZBOT*(1.0-X) + DEPTH(LL)*X
            ZMID(LG)   = ZGRAV
                  ! Could be moved outside this subroutine
                  ! as long as depths are fixed.

$if defined DEBUG_TRANV1_DETAIL
            IF (DEBUG)
     &        write(DEBUG_UNIT,'(1X,A5,I5/3(3(2x,A12,''='',G13.6)/))')
     &           'LL:',         LL,
     &           'AREA_BELOW',  AREA_BELOW,
     &           'AREA(LG)',    AREA(LG),
     &           'X',           X,
     &           'ZBOT',        ZBOT,
     &           'DEPTH(LL)',   DEPTH(LL),
     &           'ZGRAV_BELOW', ZGRAV_BELOW,
     &           'ZGRAV',       ZGRAV
$endif

            IF (LL .LT. LAYERS) THEN
C           ..... Effective depth interval;
C                 distance between gravity centers of layers:
               DZ = ( ZGRAV_BELOW - ZGRAV )
C           ..... Density gradient between layer LG+1 and LG:
               NSQUARED = grav*(DENS(LG+1) - DENS(LG)) /Rho_Approx /DZ
C           .... lower limit (--> upper limit to vertical mixing):
               BW_FREQ = SQRT( MAX( 1.0d-20, 
     &          ABS(DBLE(N2LIM)), NSQUARED ) )
C                    Unit 1/s2 = ((m/s2)/(kg/m3))*(kg/m3)/m

C           ...... Stability information used in next sequence,
C              and exported:
               BWFREQ(LG+1) = BW_FREQ

C           ...... Redistribution of mixing energy:
C                  Reduction factor between layers LG+1 and LG:
               P = ABS( GMIXDC*(BW_FREQ**GMIXDX)*DZ )

               if (P.lt.100.0) then
                   P = Exp(-P)
               Else
                   P = 0.0
               Endif
               WK_ARRAY(LL) = P ! Reduction factor between LL and LL+1

$if defined DEBUG_TRANV1_DETAIL
               IF (DEBUG)
     &             write(DEBUG_UNIT,'(1X,3(2X,A,'':'',G18.11))')
     &                'DZ',          DZ,
     &                'NSQUARED',    NSQUARED,
     &                'BW_FREQ',     BW_FREQ,
     &                'P',           P
$endif
            ENDIF

C       ...... Save depth, gravity and surface area
C              as new info. about lower layer:
            ZGRAV_BELOW = ZGRAV
            AREA_BELOW  = AREA(LG)  ! at top
            ZBOT = DEPTH(LL)        ! of layer LG

         END DO ! move to layer above
         
         BWFREQ(LG) = BW_FREQ ! value or this index not used, 
                              ! just included to get reasonable value



$if defined DEBUG_TRANV1
      IF (DEBUG)
     &   WRITE( DEBUG_UNIT,*)' mixing energy redistribution scheme:'
$endif


C --------------------------------------------------------------
C Set values at endpoints to initiate recursions:
         WK_ARRAY(0) = 0.0       ! P(0) (surface)
         WK_Array(Layers) = 0.0  ! P(N) (bottom)

C --------------------------------------------------------------
C Calculate
C   f(i,minus) = SUM(j=1,..i-1 of DZ(j)*Product(s=j...i-1 of P(s) ) )
C      by recursion starting from f(1,minus) = 0
C      and store in second part of WK_ARRAY:

$if defined DEBUG_TRANV1
      IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,A5,3A18)' )
     &        'LL', ' E=VMIX(2,LG) ', ' f(LL,minus) ', ' DZ(next)'
$endif
         DZ = 0
         DO LL = 1, LAYERS, 1
            X = WK_array(LL-1)*(DZ+X) ! = f(LL,minus) ( =0 for LL=1)
            WK_Array(Layers+LL) = X
            DZ = Depth(LL+1)-Depth(LL) ! For use in next layer

$if defined DEBUG_TRANV1
            IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,I5,3G18.11)' )
     &                LL, VMIX(2,LL+LVBASE), X, DZ
$endif

         ENDDO


C --------------------------------------------------------------
C Calculate:
C   f(i,plus) = SUM(j=i...n of Product(s=i...j-1 of P(s) ) )
C      by recursion starting from f(n,plus) = 1.0
C      and store E(i)/F(i) = E(I)/[ (f(i,minus)+f(i,plus))*DZ ]
C      in second part of WK_ARRAY. (DZ = layer thickness)
C   w(i,plus) = SUM(j=i+1...n of E(i)/F(i)*Product(s=i...j-1 of P(s)))
C      by recursion starting from w(n,plus)=E(n)/F(n)
C      and store in VMIX(1,..) (global layer index)

$if defined DEBUG_TRANV1
         IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,A3,2A14,3A16)' )
     &       'LL',' f(LL,plus) ',' DZ ',
     &            ' F(LL) ',' E/F ',' w(LL,plus) '
$endif

         SUM1 = 0.0
         DO LL = LAYERS, 1, -1
            LG = LL + LVBASE ! global layer index
            SUM1 = VMIX(2,LG) + SUM1
            LG = LL + LVBASE
            DZ = Depth(LL+1)-Depth(LL)
            X = DZ + WK_array(LL)*X ! = f(LL,plus) (=1 for L= Layers)
            F = X + WK_Array(Layers+LL)
            E_ON_F = VMIX(2,LG)/F
C                   ( denominator  F(i) >=1.0 always )
            WK_Array(Layers+LL) = E_ON_F
            W = E_ON_F + WK_array(LL)*W ! = w(i,plus)
            VMIX(1,LG) = W

$if defined DEBUG_TRANV1
            IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,I3,2F14.10,3G16.10)' )
     &                    LL, X, DZ, F, E_ON_F, W
$endif

         ENDDO

C --------------------------------------------------------------
C Calculate
C   w(i,minus) = SUM(j=1,..i-1 of E(i)/F(i)*Product(s=1...i-1 of P(s)))
C      by recursion starting from w(1,minus) = 0
C      and store W(I,*) = DZ*(w(i,plus)+w(i,minus))
C      in VMIX(2,LG).
C
C VMIX(1,LG) becomes free for other uses, and is used to store 
C temporarily the additional mixing due to surface mixing 
C specified in SFMIXC and SFMIXZ. The result is used below.

$if defined DEBUG_TRANV1
      IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,A5,4A18)' )
     &      ' LL ', ' E_ON_F(LG-1) ',' w(LL,minus) ',
     &      ' W(LL,*) ', ' Energy=VMIX(2,LG) '
$endif

         SUM2 = 0.0
         E_ON_F = 0.0
         
	    ! variables F and X used as temporary variables:     
         SFX_ABOVE =0.0
         F = 1.0

         DO LL = 1, LAYERS, 1
            LG = LL + LVBASE
            W = WK_array(LL-1)*( E_ON_F + W )
                          ! initiated by WK-array(0)=0
            DZ=Depth(LL+1)-DEPTH(LL)
            VMIX(2,LG) = (W+VMIX(1,LG))* Dz
                   ! Redistributed
            SUM2 = VMIX(2,LG) + SUM2

$if defined DEBUG_TRANV1
            IF (DEBUG) WRITE( DEBUG_UNIT,'(1x,I5,4G18.10)' )
     &            LL, E_ON_F, W, (W+VMIX(1,LG)), VMIX(2,LG)
$endif
            
            E_ON_F = WK_ARRAY(LAYERS+LL)

            IF (LL .EQ.1 ) cycle

            ! -------- Top layer mixing reduced with depth:
            SFX = - (Depth(LL)-SFMIXZ(IB,1))/SFMIXZ(IB,2) 
                           ! (decreases with depth)
            X = max(-30.0D0, min( 30.0D0, SFX) - SFX_ABOVE )
                           ! >-30, =30 - SFX_ABOVE if SFX>30 
                           !       = max(-30,min(30,SFX)) first time
            SFX_ABOVE = X+SFX_ABOVE ! for next layer 
                           ! = 30 if SFX>=30, else reduced
            F = F*EXP( X )
            VMIX(1,LG-1) = F/(1.0+F)
C                       ( Reduction could be stability dependent,
C                         would be included in VMIX(2,LG) above )

$if defined DEBUG_TRANV1_SFMIX
            IF (DEBUG) WRITE( DEBUG_UNIT,
     &                        '(1x,A3,A5,5a12/1X,2I5,5G12.4)' )
     &         'LL', 'LG', 'SFX', 'X', 'SFX_ABOVE', 'F', 'VMIX(1,LG)',
     &         LL, LG, SFX, X, SFX_ABOVE, F, VMIX(1,LG-1)
$endif
         
         ENDDO

$if defined DEBUG_TRANV1
         IF (DEBUG .or. ABS(SUM1-SUM2) .gt. 1.e-6*SUM1 )
     &      WRITE( DEBUG_UNIT,'(1x,A,i5,2(A,G15.8)/20x,A,G15.8)')
     &         'Energy conservation basin',IB,
     &         '  Sum before:',SUM1,' after: ',SUM2,
     &         'Difference:', SUM1-SUM2
$endif

C ===============================================================
C Set up vertical mixing as water transports in VTRANS:
C                                        down:         up:
C ------- AREA(LG-1), DEPTH(LL-1), VTRANS(1,LG-1), VTRANS(2,LG-1)
C              Layer LL-1: DENS(LG-1)
C ------- AREA(LG),   DEPTH(LL),   VTRANS(1,LG),   VTRANS(2,LG)
C              Layer LL  : DENS(LG)
C ===============================================================

C   ........... Scan layers from bottom up by local index LL

$if defined DEBUG_TRANV1_TSTEP_LIM
      LIMIT_ADV  = 0
$endif

         DIFF_BELOW = 0.0
         ZBOT= DEPTH(LAYERS) ! TOP OF BOTTOM LAYER




         DO LL = LAYERS, 2, -1   ! Local index, used for depth
                    ! ( tar alle lag, OK med endret kode i TRANSP_U.for )   
         
            LG = LL + LVBASE           ! Global index
            if (LL.gt. NLVOPN(IB)) then
               EFFECTIVE_AREA = AREA(LG)
            else   
               EFFECTIVE_AREA = AREA(LG)*VLCORF(IB)
            endif ! Reduced area if volume is smaller 
                  ! (could use VBUF also)

C           ..... Effective depth interval;
C                 distance between gravity centers of layers:
            DZ  = ZMID(LG)-ZMID(LG-1)

C           ..... Density gradient:
            BW_FREQ = BWFREQ(LG)
            NSQUARED = BW_FREQ*BW_FREQ

C  ............ Tidal mixing:
            DFCOEF  = MIXFAC* MIXCONST(IB)  ! = Energy
                 ! MIXCONST is scaled with N2SCAL in calling subroutine


$if defined DEBUG_TRANV1_DETAIL
      IF ( DEBUG )THEN
          WRITE (debug_unit, '(2(1x,A12,'':'',G17.11))')
     &         'MIXFAC'       , MIXFAC
     &        ,'MIXCONST(IB)' , MIXCONST(IB)
     &        ,'---> DFCOEF'       , DFCOEF

      ENDIF
$endif

            DFCOEF = DFCOEF/(BW_FREQ**MIXEXP)  ! = mixing coeff.
     &               + VMIX(1,LG-1)*SFMIXC(IB)/NSQUARED
            ! Unit m2/s  (SFMIXC has unit m2/s3, 
            !             VMIX(1,LG-1) dim.less, from loop above

$if defined DEBUG_TRANV1_DETAIL
      IF ( DEBUG )THEN
          WRITE (debug_unit, '(2(1x,A12,'':'',G17.11))')
     &         'BW_FREQ'     , BW_FREQ
     &        ,'VMIX(1,LG-1)'  , VMIX(1,LG-1)
     &        ,'---> DFCOEF' , DFCOEF

      ENDIF
$endif


C     Constant to transform between diffusive transport (m3/day)
C     and effect against gravity: (kg m/s2) m/day.
            X = NSQUARED *(DZ**2) * Rho_Approx ! (kgm/s2/m2)
            if (X.le.0.0) then
               write(*,'(1x,A,2i3,3(1x,A,'':''G12.6))')
     &            'LL,LG:',LL, LG, 'NSQUARED', NSQUARED,
     &            'Rho_Approx', Rho_approx, 'DZ', DZ
!               Pause
            endif


C    ....... Tidal mixing + surface mixing due to constant energy:
            DIFF_TRANSP = EFFECTIVE_AREA*DFCOEF/DZ * sec_per_day
C             m3/day    =   m2    * m2/s / m * s/day

C        ---> effect used against gravity is stored as info:
                  !    VMIX(1,LG)  = DIFF_TRANSP * X
                  !    (kgm/s2)*m/day = (m3/day) * (kgm/s2/m2)
                  !    endret til:    
                  !  VMIX(1,LG-1) = DFCOEF * NSQUARED * RHO_APPROX /  GRAV
C              kg/m2/s   =  m2/s  *  (1/s2)  *   (kg/m3)  / (m/s2)
                  !    NB forsk›vet indeks, slik at LSURF blir satt


C    ....... Add mixing due to gravitationally released energy:
            DIFF_TRANSP =   DIFF_TRANSP 
     &                    + GMIXFR(IB) * VMIX(2,LG)   / X
            !  m3/day =                  (kgm2/s2/dag) / (m2/s2*kg/m3)


C     Utvidet VMIX:
            VMIX(1,LG-1) = DIFF_TRANSP *NSQUARED * RHO_APPROX / GRAV 
     &                           / EFFECTIVE_AREA * DZ / sec_per_day
            ! kg/m2/s    =  m3/day       * (1/s2)  *  (kg/m3) / (m/s2)
            !                    /   m2     * m  /  (s/day)



C     Effect applied to each layer is used for
C     mixing with layer above.
C                GMIXFR = fraction of total effect used to work against 
C                         gravitational field (assumed constant).

C     ..... Vertical advection used to limit time-step:
            ADV_TRANSP = VTRANS(1,LG)
C                         >0 for flow from L to L-1,
C                         <0 for flow from L-1 to L
C          Add to gross outflow term from horisontal transports:
            IF (ADV_TRANSP.lt.0.0) THEN
                VTRANS(2,LG-1) = VTRANS(2,LG-1) - ADV_TRANSP
                     ! flow down, increases gross throughflow
                     ! for layer above
                OUT_FLUX = VTRANS(2,LG) 
            ELSE
                OUT_FLUX = VTRANS(2,LG) + ADV_TRANSP
            ENDIF    ! Increases gross through-flow from this layer

C      ........... Limit time-steps:
$if defined DEBUG_TRANV1_DETAIL

      IF ( DEBUG )THEN
          write (debug_unit, *)
     &      '----- transports between layer LG=', LG, ' and ', LG-1
          WRITE (debug_unit, '(2(1x,A12,'':'',G17.11))')
     &         'Nsquared'    , NSQUARED
     &        ,'VMIX(1,LG)'  , VMIX(1,LG)
     &        ,'VMIX(2,_)'   , (VMIX(2,LG)+VMIX(2,LG-1))/2
     &        ,'AREA(LG)'    , AREA(LG)
     &        ,'EFFECTIVE_AREA', EFFECTIVE_AREA
     &        ,'DFCOEF'      , DFCOEF
     &        ,'DIFF_TRANSP' , DIFF_TRANSP
     &        ,'ADV_TRANSP'  , ADV_TRANSP
     &        ,'OUT_FLUX'    , OUT_FLUX

      ENDIF
$endif

$if defined DEBUG_TRANV1_TSTEP_LIM
            ADV_STEP_DECREASED = .FALSE.
$endif

C      ....... General step according to gross outflux from layer LG:
C              giving no more than 30% of water renewal
C              to ensure stability and reasonable accuracy.

C            write(*,*)'lg, vlayer, OUT_FLUX',
C     &             lg, vlayer(lg), out_flux

            if (out_flux.gt.0.0) then
               TSTEP_NEW = MAX_FRACTION_REPLACED*VLAYER(LG)/ 
     &                     MAX(0.1D0*VLAYER(LG),OUT_FLUX)
C                       minimum limit in divisor
C                       to avoid floating point overflow
               IF (TSTEP_ADV .GT. TSTEP_NEW) THEN
                   TSTEP_ADV = TSTEP_NEW
$if defined DEBUG_TRANV1_TSTEP_LIM
                   ADV_STEP_DECREASED = .true.
$endif
               ENDIF
            Endif

            DIFF_BELOW = MAX(0.0D0,DIFF_TRANSP) ! (Used for next layer)


$if defined DEBUG_TRANV1_TSTEP_LIM
      IF ( DEBUG )THEN
         IF( ADV_STEP_DECREASED ) THEN
            LIMIT_ADV = LG
            LIMIT_STORE (1) = VLAYER(LG)
            LIMIT_STORE (2) = ADV_TRANSP
            LIMIT_STORE (3) = OUT_FLUX
            LIMIT_STORE (4) = TSTEP_ADV
         ENDIF
       ENDIF
$endif

C       ........ Store diffusive transports in VTRANS(2,LG):
            VTRANS(2,LG) = DIFF_TRANSP
C       ...... Save depth og gravity and surface area
C              as new info. about lower layer:
            ZGRAV_BELOW = ZGRAV
            AREA_BELOW  = AREA(LG)
            ZBOT = DEPTH(LL-1) ! TOP OF LAYER LG-1
         END DO


$if defined DEBUG_TRANV1_TSTEP_LIM
      IF ( DEBUG )THEN

        IF (LIMIT_ADV.GT.0) THEN
          write (debug_unit, *)
     &      '----- adv. time step control from layer LG=', LIMIT_ADV
          WRITE (debug_unit, '(3(1x,A12,'':'',G18.11))')
     &         'VLAYER(LG)', LIMIT_STORE (1)
     &        ,'ADV_TRANSP', LIMIT_STORE (2)
     &        ,'OUT_FLUX'  , LIMIT_STORE (3)
     &        ,'TSTEP_ADV' , LIMIT_STORE (4)
        ENDIF
      ENDIF
$endif


C -----------------------------------------------------
      END DO    ! Next basin
C -----------------------------------------------------

      END Subroutine



C ==================================================================
      SUBROUTINE TRANV2( DEBUG, TSTEP, NB, INDXI, NL, VLAYER,
     &        MLV, VTRANS, TRCF_RANGE, LENGTH_TRCF, TRCF, VTRNEG )
      

C         Adjusts vertical mixing in internal basins for given TSTEP
C ==================================================================
C ----------------------- External arguments: ------------------------
C INPUT:
      LOGICAL DEBUG
C           Controls debug printouts
      real*8 TSTEP
C           time-step.
      INTEGER NB
C           Number of basins
      INTEGER INDXI(NB+1)
C           Layer index limits for each basin (NB+1 values)
      INTEGER NL
C           Number of layers in array from outside TRANSP Modules
      real*8 VLAYER(NL)
C           Volume of layers
      INTEGER MLV  ! Allocated length of VTRANS and TRCF_RANGES

C Input/Output:
      real*8 VTRANS(2,MLV)  ! m3/day
C        I: (1,L): Vertical flux up from layer L to layer L-1,
C                  through adjusted limit (constant fraction of volume),
C                  (Values from TRANC, summed from bottom up,
C                   and adjusted for variation in total volume)
C           (2,L) Diffusive flux between layers L and L-1,
C                   indicator value -1 for well_mixed layer.
C        O: Final description of vertical transports between layers:
C             (1,L): (m3/day) Unchanged
C             (2,L): (m3/day) Adjusted diffusive flux

C           To avoid artificial increase of diffusion by advection,
C           part of the advective transport is subtracted from
C           diffusive transport in opposite direction instead of
C           being added to diffusive transport in advective direction.
C           See how VTRANS is used in TRANSP_APPLY.

C Output:
      INTEGER TRCF_RANGE (0:1, MLV)
      INTEGER LENGTH_TRCF
      real*8  TRCF (LENGTH_TRCF) ! Transfer coefficients

      LOGICAL VTRNEG
C          TRUE if some of VTRANS terms are < 0, FALSE if not.


C -------------------- local variables ----------------
      INTEGER IB, LSURF, LBOTTOM, LNUM, SOURCE, LG
      real*8 ADV_TRANSP, DIFF_TRANSP, ABS_ADV, A

      INTEGER TRCF_BASE_INDEX, SOURCE_INDEX, TARGET_INDEX
      INTEGER LIMIT(0:1), START_DIRECTION, DIR, LSHIFT, L, COUNT

      EQUIVALENCE (LIMIT(0),LSURF),(LIMIT(1),LBOTTOM)

      LOGICAL GO_ON, REVERSED
      real*8 VSUM, C_RED, RED_FAC, X

$if defined DIFF_COMP
      real*8 O_5 / 0.5 / ! =0.5 theoretically, set lower to try to
                       ! countereffect numerical instabilities
$endif
C ----------------------------------------------------------------------
C Final vertical transport description is set up:

$if defined DEBUG_TRANV2 || defined DEBUG_TRANV2_VTRNEG
!      INCLUDE 'DEBUG.INC'
      INTEGER I,K
$endif

$if defined DEBUG_TRANV2
      IF (DEBUG)
     &  WRITE(DEBUG_UNIT,*) ' ================= TRANV2,  TSTEP =', TSTEP
$endif

      VTRNEG = .FALSE.

      TRCF_BASE_INDEX = 0

C ------------------------------------------
      DO IB = 1,NB   ! Loop through basins
C ------------------------------------------

         LSURF   = INDXI(IB)+1
         LBOTTOM = INDXI(IB+1)
         LNUM    = INDXI(IB+1)-INDXI(IB)

$if defined DEBUG_TRANV2
         IF (DEBUG) THEN
      write (DEBUG_UNIT,*)'======= Basin ',IB,' at entry:'
      WRITE(DEBUG_UNIT,'(1X,A5,3A16)' )
     &      'K','VLAYER(K)','VTRANS(1,K)','VTRANS(2,K)'
      WRITE(DEBUG_UNIT,'(1X,I5,3E16.8)')
     &      (K,VLAYER(K),(VTRANS(I,K),I=1,2),K=LSURF,LBOTTOM)
      WRITE(DEBUG_UNIT,'(1X,A5,3A16)' )
     &           'LG', 'A', 'ADV_TRANSP', 'DIFF_TRANSP'
         ENDIF
$endif

C -------------------------------------------------------------------
         DO LG = LSURF+1, LBOTTOM ! Loop through layer interfaces
C -------------------------------------------------------------------

C               Transports at top of surface layer not processed,
C               not used if horisontal transports are active,
C               carries special withdrawal if horisontal transports
C               are turned off.
C
C        ...... Advective transport = accumulated
C               net advective flux from layer LG to layer LG-1
C               across (adjustable) layer boundary
            ADV_TRANSP  = VTRANS(1,LG)       ! <0 down, >0 up
C      ....... Diffusive transports between layer LG and LG-1:
            DIFF_TRANSP = VTRANS(2,LG)

            IF (DIFF_TRANSP .GE. 0.0 ) THEN

C      --------------------------------------------------------------
C      Below well_mixed layer, combine DIFF & ADV transp.
         ! N† for alle lag
C        Q1 = DIFF + |ADV|*0.5*(1+a) in direction of ADV.
C        Q2 = DIFF - |ADV|*0.5*(1-a) in opposite direction.
C      Directed sum of transports Q1-Q2 = |ADV|*(0.5+0.5a+0.5-0.5a)=|ADV|
C        where   a = v*dt/dz, with v= adv. velocity = Qadv/Area,
C        so that a = Qadv*dt/Volume
C          (dt = time step, dz = vertical distance between layers.)
C      This is effected by subtracting |ADV|*0.5*(1-a)
C      from diffusive transport in VTRANS(2,..) (working both ways)
C      and keeping advective transport in VTRANS(1,..) unchanged.

C xxxx   However, a is always set to give Q2>=0, that is a>1-DIFF/|ADV|,
C xxxx   to avoid negative diffusion which can give instabilities.

C        With full correction this largely eliminates diffusive effect of
C        advection, which would otherwise spread out a peak distribution
C             0 0 1 0 0 0 0 0 .....
C        into a binomial distribution after n steps:
C             f(i) = n!/i!(n-i)!*p**i*(1-p)**(n-i), i=0,n
C        with variance np(1-p)   ( p = flux*dt/layer_volume)
C
C      This can be seen by direct numerical experiment,
C      and from the fact that distribution of a mass of 1.0
C      from a single point over time T with diffusion coeff. D
C      has variance 2*D*T. The variance of the binomial distribution
C       = np*(1-p) = Qadv*(n*dt)*(1-p) = Qadv*T*(1-p),
C      thus to counteract this the diffusion must be reduced by
C         dD = Qadv*T*(1-p)/(2*T) = 0.5*Qadv*(1-p)
C      where Qadv now is relative flushing (unit 1/time)

                ABS_ADV = ABS(ADV_TRANSP)
                A = TSTEP*ABS_ADV/
     &               ( (VLAYER(LG-1)+VLAYER(LG)) )
                     ! (factors 0.5 on both TSTEP and äVLAYER cancel)

$if defined DEBUG_TRANV2
               IF (DEBUG) THEN
                  WRITE(DEBUG_UNIT,'(1X,I6,3E16.8)')
     &                 LG,   A,   ADV_TRANSP,   DIFF_TRANSP
               ENDIF
$endif

C ------- Diffusion scaled down to compensate for numerical diffusion:
               A = MAX( 0.0D0, DIFF_TRANSP-0.5*ABS_ADV*(1-A) )
$if defined DIFF_COMP
               if ( A .eq. 0.0) then
                   A = min( 0.0D0, DIFF_TRANSP - O_5*ABS_ADV*(1-A) )
                        ! O_5 : theoretically 0.5, adjusted value above.
               endif
$endif 

               VTRANS(2,LG) = A

               if ( A .lt. 0.0) then  ! Signal negative values:

$if defined DEBUG_TRANV2_VTRNEG
                   IF ( DEBUG. and. .not. VTRNEG ) THEN
                       WRITE(DEBUG_UNIT,*) ' In TRANV2: VTRNEG=.t.'
                   ENDIF
$endif
                   VTRNEG = .TRUE.
               Endif
C      --------------------------------------------------------------
            ELSE
C      Well-mixed layer:
C           No internal diffusion, to minimize numerical noise.
C           Advection retained to conserve substances.
               VTRANS(2,LG) = 0.0 ! (could be done in the first place?)
            ENDIF

C -------------------------------------------------
         END DO  ! Next layer interface
C -------------------------------------------------


$if defined DEBUG_TRANV2
         IF (DEBUG) THEN
            write (DEBUG_UNIT,*)'======= At end of TRANV2, part 1:'
            WRITE(DEBUG_UNIT,'(1X,A5,2A16)' )
     &             'K','VTRANS(1,K)','VTRANS(2,K)'
            WRITE(DEBUG_UNIT,'(1X,I5,2E16.8)')
     &             (K,(VTRANS(I,K),I=1,2),K=LSURF,LBOTTOM)
         ENDIF
$endif

$if defined DEBUG_TRANV2_VTRNEG
         IF (DEBUG .and. VTRNEG ) THEN
            write (DEBUG_UNIT,*)'======= At end of TRANV2, part 1:'
            WRITE(DEBUG_UNIT,'(1X,A5,2A16)' )
     &             'K','VTRANS(1,K)','VTRANS(2,K)'
            WRITE(DEBUG_UNIT,'(1X,I5,2E16.8)')
     &             (K,(VTRANS(I,K),I=1,2),K=LSURF,LBOTTOM)
         ENDIF
$endif


C    --------- Set up transfer coefficients TRCF for diffusion
C              over half of time step, to be used in 2-step updating
C              in TRANV2 in TRANSP_V.FOR.

         if (length_TRCF .lt. TRCF_BASE_INDEX + LNUM*LNUM ) then
             STOP 'TRCF too short, see TRANV2 in TRANSP_V.FOR'
         endif

$if defined DEBUG_TRANV2
      IF (DEBUG) write (DEBUG_UNIT,*)'      ======= TRCF =====:'
$endif

C ensured by equivalence above:  LIMIT(0) = LSURF
C                                LIMIT(1) = LBOTTOM

C ----------------------------------------------------------
         DO SOURCE = LSURF, LBOTTOM  ! Loop through layers
C ----------------------------------------------------------

            SOURCE_INDEX = TRCF_BASE_INDEX + SOURCE-LSURF+1
            TRCF(SOURCE_INDEX) = 1.0
            TRCF_RANGE(0,SOURCE) = 0  ! relative to source layer
            TRCF_RANGE(1,SOURCE) = 0

$if defined DEBUG_TRCF
            IF (DEBUG) WRITE(DEBUG_UNIT,'(1X,A,I5:A,G18.11)' )
     &            '>>>>>> LAYER=', SOURCE, ' vlayer=', VLAYER(SOURCE),
     &            ' SOURCE_INDEX=',SOURCE_INDEX
$endif

C        ...... scan in each applicable direction:
            DO START_DIRECTION = -1,+1,2
C           ......... initiate scan:
               DIR = START_DIRECTION
               LSHIFT = (DIR+1)/2  ! 0/1 for directions -1/+1
               IF ( SOURCE .EQ. LIMIT(LSHIFT) ) CYCLE
               GO_ON = .TRUE.
               REVERSED = .FALSE.

               VSUM = 0.0
               L = 0
               COUNT = 0
               RED_FAC = 1.0

$if defined DEBUG_TRCF
               IF (DEBUG) WRITE(DEBUG_UNIT,* ) 'Scans direction ',DIR,
     &              ' accumulating TRCF:'
$endif

               DO WHILE ( GO_ON )

                  DO WHILE ( L+source.NE.LIMIT(LSHIFT) .and. GO_ON
     &                       .and. COUNT.le.2*LNUM )
                     COUNT = COUNT + 1

C         ........ find degree of transfer to layer L+source+DIR:
                     LG = L + SOURCE

C It has been found by analytical solution and confirmed by
C numerical integration of simple 1-dimensional diffusion
C over equidistant layers, with constant cross-sectional area,
C that a substance contained within a single node (k) at time t=0
C will be distributed at time t=T after several integration steps
C according to the following recursions from C(k):
C
C    For i>k: C(i) = C(i-1)*exp[-{2(i-k)-1}/4K(i,i-1)t]
C    For i<k: C(i) = C(i+1)*exp[-{2(k-i)-1}/4K(i,i+1)t]

C where K(i,j) = exchange coefficient between i and j; j=i+1 or i-1.
C with unit 1/time. Investigation with numerical integration shows
C that the recursions hold in most cases even if K varies with i.

C For i>k:
C By putting K(i,j) = Q(i,j)/(V(i)+V(j))/2 where
C    Q = volume exchange between nodes and
C    V = volume of each node

C replacing (2(i-k)-1) = ( i-k-1 + i-k )  by (VSUM1 + VSUM2)
C with VSUM1 = SUM from k to i-1 of V, with k and i-1 counting 50%
C      VSUM2 = SUM from k to i   of V, with k and i   counting 50%
C            = VSUM1 + 0.5*(V(i-1)+V(i))

C and putting the time t = TSTEP/2.0,

C we can write the formula:
C    i>k: C(i) = C(i-1)*exp [ - (VSUM1+VSUM2)/(2*Q*TSTEP) ]

C similarly for i<k, but with i+1 instead of i-1 in VSUM expressions:
C    i<k: C(i) = C(i+1)*exp [ - (VSUM1+VSUM2)/(2*Q*TSTEP) ]

C The exponent arguments given in brackets can be generalized to a case
C where V and Q may vary over the nodes:
C    [] = [ - VSUM * (|k-i|-0.5)/|k-i| / (Q*TSTEP) ]
C with VSUM = sum of V(j) for J = i,...,k

                     X = VSUM
                     VSUM = VSUM + (VLAYER(LG)+VLAYER(LG+DIR))/2.0
                                  ! starts at approx. 1*V
                     X = X + VSUM ! = VSUM1 + VSUM2
                                  ! starts at approx. 1*V

                     DIFF_TRANSP = Max(0.0D0,VTRANS(2,LG+LSHIFT))*TSTEP
                     A = 2.0*DIFF_TRANSP/X

                     IF ( A .GT. 0.002 ) THEN
                        RED_FAC = EXP(-1/A)
                     ELSE
                        RED_FAC = 0.0
                     ENDIF

      ! For the first transition (i=k+1), 
      ! the first order alternative to
      !        C(k+1) = C(k)*exp [ - 1/A ] , 
      !        with A = 2*Q*TSTEP/ [VSUM1+VSUM2] 
      !               = 2*Q*TSTEP/ [0.5*(V(k)+V(k+1))] 
      ! is     C(k+1) = C(k)*(A')/[1.+A'] 
      !        with A'= A/4 = Q*TSTEP/[V(k)+V(k+1)]
      ! reasonable for exchange of volume Q*TSTEP/2 
      ! between volumes V(k) and V(k+1), with mean [V(k)+V(k+1)]/2
                     
                     if (L.eq.0) then
                        C_RED = MAX( RED_FAC, A/4./(1.0+A/4.) )
                     ELSE
                        C_RED = C_RED*RED_FAC
                     ENDIF
                     


$if defined DEBUG_TRCF
         IF (DEBUG) WRITE(DEBUG_UNIT, '(1x,3A15/1x,3G15.7)')
     &           'VSUM', 'DIFF_TRANSP', 'C_RED',
     &           VSUM, DIFF_TRANSP, C_RED
$endif

C                 ......... store concentrations relative to
C                           remaining conc. in SOURCE in array TRCF:
                    IF ( C_RED*LNUM .GE. 0.001 .or. L.eq.0 ) THEN
                        L = L + DIR
                        TARGET_INDEX = SOURCE_INDEX + L

$if defined DEBUG_TRCF
         IF (DEBUG) WRITE(DEBUG_UNIT,'(1x,A,I5,A,I5)')
     &         ' L',L, ' TARGET_INDEX', TARGET_INDEX
$endif

                        IF ( .NOT. REVERSED ) THEN
                           TRCF_RANGE (LSHIFT,SOURCE) = L
                           TRCF( TARGET_INDEX ) = 0.0
                        ENDIF
                        TRCF(TARGET_INDEX) = MIN( 1.0D0,
     &                                   TRCF(TARGET_INDEX)+C_RED )
C                           (max. conc. = remaining in source layer )
                     ELSE
                        GO_ON=.FALSE.
                     ENDIF

$if defined DEBUG_TRCF
         IF (DEBUG) WRITE(DEBUG_UNIT, '(1x,A5,L4)') 'GO_ON', GO_ON
$endif

                  ENDDO

C Reached surface or bottom: the propagation is reflected,
C and contributions added, in order to achieve good results
C when there is high mixing intensity close to boundary:
                  if ( .not. REVERSED) THEN
                     DIR = -DIR
                     REVERSED = .TRUE.
                  ELSE  ! but only once!
                     GO_ON = .FALSE.
                  ENDIF
               ENDDO ! Continue after scan reversal
            ENDDO ! Next START DIRECTION

$if defined DEBUG_TRCF
          IF (DEBUG) THEN
             WRITE(DEBUG_UNIT,'(1x,A:2I5)')
     &          'TARGET RANGE:',
     &          (TRCF_RANGE(L,SOURCE), L=0,1),'SUM TRCF:'
          ENDIF
$endif

C ----------------------------------------------------------------
C TRCF now contains resulting concentrations in layers in RANGE
C in units relative to remaining concentration in SOURCE layer.
C Now converts TRCF to fraction of original concentration:

            VSUM = 0.0
            DO L = TRCF_RANGE(0,SOURCE), TRCF_RANGE(1,SOURCE)
                TARGET_INDEX = SOURCE_INDEX + L
C                                L is now relative to SOURCE

$if defined DEBUG_TRCF
          IF (DEBUG) THEN
             WRITE(DEBUG_UNIT,'(1x,A,I5,A,I5)')
     &         'L:',L, ' TARGET_INDEX', TARGET_INDEX
          ENDIF
$endif

                TRCF(TARGET_INDEX) = TRCF(TARGET_INDEX)*VLAYER(L+SOURCE)
                VSUM = VSUM + TRCF(TARGET_INDEX)
            ENDDO
C TRCF = amount of substance from SOURCE in each TARGET layer
C VSUM = total amount of substance,
C both in concentration units relative to remaining conc. in SOURCE

            C_RED = VLAYER(SOURCE)/VSUM != reduction in source layer

$if defined DEBUG_TRCF
          IF (DEBUG) THEN
             WRITE(DEBUG_UNIT,'(1x,2A15/1x,2G15.7)')
     &          'VSUM','C_RED',
     &           VSUM, C_RED , ' FINAL VALUES:'
          ENDIF
          VSUM = 0.0
$endif

C Now TRCF is converted to amount of substance found in each layer
C relative to initial total amount in SOURCE layer:

            DO L = TRCF_RANGE(0,SOURCE), TRCF_RANGE(1,SOURCE)
               TARGET_INDEX = SOURCE_INDEX + L
               TRCF(TARGET_INDEX) = TRCF(TARGET_INDEX)*C_RED

$if defined DEBUG_TRCF
          IF (DEBUG) WRITE(DEBUG_UNIT,'(2(1x,A,I3),2(1x,A,G18.11))')
     &         'L:',L , 'TARGET_INDEX:', TARGET_INDEX,
     &         'TRCF:', TRCF(TARGET_INDEX),
     &         'VLAYER:', VLAYER(source+L)
          VSUM = VSUM + TRCF(TARGET_INDEX)
$endif

            ENDDO
C that means: TRCF is now net fraction of
C initial volume in SOURCE layer transferred to TARGET layers.

$if defined DEBUG_TRCF
          IF (DEBUG) WRITE(DEBUG_UNIT,'(2(1x,A,G18.11))')
     &         'SUM of TRCF:',VSUM, 'VLAYER:', VLAYER(SOURCE),
     &         'diff:', VSUM - VLAYER(SOURCE)
$endif

C --------------------------------------------------
            TRCF_BASE_INDEX = TRCF_BASE_INDEX + LNUM
         ENDDO ! Next SOURCE layer
C --------------------------------------------------

C --------------------------
      END DO  ! NEXT BASIN.
C --------------------------

      RETURN

      END Subroutine

      End Module