      Module fx_Transp_u

      use fx_RunControl, only: DEBUG_UNIT
      use fx_Transp_h
      
      implicit none
      
      
C Group_number no longer needed

C ==================================================================
C Eutrophication model   - File:    TRANSP_U.FOR
C                                   Birger Bjerkeng, NIVA.

C  CONTROLS DEBUG DUMP OF CALCULATIONS ON FILE:
$define DEBUG_UPDATE 0
C    0 : NO DUMP
C    1 : DUMP END RESULTS
C    2 : DUMP INTERMEDIATE AND END RESULTS


$undefine DEBUG_DERIV

$if defined DEBUG_DERIV
$define DEBUG_DERIV_GT_1
$endif

$define DEBUG_neg_diff

      Contains
      
C =====================================================================
      SUBROUTINE TR_CONC_DERIV (
     &                 HTR_CALC, DEBUG, NBI, INDXI, NLI,
     &                 VLAYER, NBE, INDXE, NLE,
     &                 NC, INDXC, BConn1, BConn2, HTRMIX,
     &                 GROUP_NUMBER,
     &                 CONCI, CONCE, CDERIV, IMPORT,
     &                 MW, CDERIV_R8, MT_TEMP, IMPORT_R8, MLV, VTRANS,
     &                 MQ_TR, NQ_TR, QJET_TRANSP_INDEX, QJET_TRANSP,
     &                 MLH, HTRANS, HTR_L, MC, INDXH, NAME )
      
      
C ----------- Subroutine arguments as defined for TRANSP_SETUP above:
      LOGICAL HTR_CALC   ! turns off calculation of horisontal transports
      LOGICAL DEBUG
      INTEGER NBI
      INTEGER INDXI(NBI+1)
      INTEGER NLI
      real*8    VLAYER(NLI)
      INTEGER NBE
      INTEGER INDXE(NBE+1)
      INTEGER NLE
      INTEGER NC
      INTEGER INDXC(NC+1)
      INTEGER BConn1(NC), BConn2(NC)
      real*4    HTRMIX(NC)
              ! Degree of mixing between contiguous transports
              ! in the same direction

C --------- additional arguments:
C ------ In:

      INTEGER GROUP_NUMBER ! =1 or 2:
C           if 2, do not initiate IMPORT, only update input values
C           (for phytoplankton groups, import summed for total)
      real*8    CONCE  (NLE)
C              - Concentr. of substance in layers in external basins.
      real*8    CONCI (NLI)
C              - Concentrations in internal basins

C ------ Out:
      real*8 CDERIV(NLI)
C              - Time derivative of concentration in internal layers
C                not taking vertical diffusion into account
C                       (unit cons/day)
      real*8 IMPORT(NBI)
C              - Total IMPORT to each basin from other internal
C                and external basins.   (unit m3*cons/day)
C                (-export by surface withdrawal
C                  if horizontal transports are turned off)
C --------- as defined above:

      INTEGER MW
      real*8  CDERIV_R8   (0:MW)
      INTEGER MT_TEMP
      real*8  IMPORT_R8   (MT_TEMP)
C       ...... Work array used for better precision
C              when adding transports.

      INTEGER MLV, MLH, MC
      real*8  VTRANS        (2,MLV)

      INTEGER MQ_TR, NQ_TR
      INTEGER QJET_TRANSP_INDEX(2,MQ_TR)
      real*8 QJET_TRANSP (MQ_TR)

      real*8 HTRANS        (MLH)
      INTEGER HTR_L        (MLH,2)
      INTEGER INDXH        (MC+1)

      CHARACTER*(*) NAME



C -------------------- local variables -------------------------

      INTEGER IB, L, ILH, NLH, NLH_DIM, ILC, IC
      INTEGER IBA, IBB, ILA, ILB, NLA, NLB

      INTEGER LSURF
      real*8 MTDOWN

            ! Dummy arg. in subroutine - not used:
      real*8 CDERIV_EXT(1)
            ! declared as arrays to get
            ! compatibility of argument type

      real*8 EXPORT

      real*8 TRANSPORT



C Initiates time derivative of conc. in each layer with transport terms.
C Total net material exports are accumulated for each basin.

C Effect of volume changes on conc. derivative handled by
C separate adjustment between integration steps in subroutine CCONC.

C Called from ACSL with NLI = dimensioned number of layers,
C                           = number of layers integrated by ACSL.


C Unused MDERIV values are zeroed ( NLI = ACSL DIMENSION ):

$if defined DEBUG_DERIV
!      INCLUDE 'DEBUG.INC'
      real*8 MASS_INCREASE
      LOGICAL DEBUG_PRINT
      DEBUG_PRINT = DEBUG

      IF (DEBUG_PRINT) THEN
         WRITE( DEBUG_UNIT, '(1X,A,A)') ' TR_CONC_DERIV for ', NAME
         WRITE( DEBUG_UNIT, '('' NBI='',I5,'' NLI='',I5,''  INDXI:'')')  
     &           NBI,NLI
         WRITE( DEBUG_UNIT, '(10X,10I6)' )  INDXI
      ENDIF
$endif


C ----- A. Calculate transports as mass flux (amount/day):

C   ------------ Two-way transport between layers:

$if defined DEBUG_DERIV
            IF (DEBUG_PRINT)
     &      WRITE ( DEBUG_UNIT, '(1X,2A6,3A15)' )
     &              'IB', 'L', 'MTDOWN', 'CDERIV_R8(L-1)','CDERIV_R8(L)'
$endif

      DO IB = 1,NBI

C     ..............
         IMPORT_R8(IB) = IMPORT(IB) ! initiate 8-byte import term
C     ..............

         LSURF = INDXI(IB)+1
         CDERIV_R8 (LSURF) = 0.0
         
         DO L = LSURF+1, INDXI(IB+1)
C        ...... Net downwards mass transport from layer above:
C    Advective and diffusive transports separately in VTRANS(1..)
C    and VTRANS(2..) respectively:

C    ... Vertical advection:

            IF ( VTRANS(1,L).lt.0.0 ) then
               MTDOWN = - VTRANS(1,L) * CONCI(L-1)
            ELSE
               MTDOWN =  - VTRANS(1,L) * CONCI(L)
            ENDIF
C    ... Positive diffusive mixing not included here, in TRCF instead.

C    Update derivatives with combined effect:
            CDERIV_R8(L-1) = - MTDOWN + CDERIV_R8(L-1)
            CDERIV_R8(L)   = + MTDOWN ! initiated

$if defined DEBUG_DERIV
             IF (DEBUG_PRINT)
     &      WRITE ( DEBUG_UNIT, '( 1X, 2I6, 3E18.11)' )
     &              IB, L, MTDOWN, CDERIV_R8(L-1), CDERIV_R8(L)
$endif

         END DO

C     ------- when horisontal transports are turned off: balance
C             net volume influx stored in VTRANS(1,LSURF) with
C             a redrawal of water from the surface layer to take care
C             of mass and volume balance with unchanged total volumes:
         IF (.not. HTR_CALC) THEN
             TRANSPORT =  VTRANS(1,LSURF)*CONCI(LSURF)
             CDERIV_R8(LSURF) = CDERIV_R8(LSURF) - TRANSPORT
             IMPORT_R8(IB) = IMPORT_R8(IB)    - TRANSPORT
         ENDIF
      END DO
C            NOTE! Bottom exchange not included here!
C                  Between model layers within well mixed volume:
C                  mass transport <>0, proportional to volume transport.



C   ------------- Add effect of vertical mixing from dived jets:

$if defined DEBUG_DERIV_GT_1
            IF (DEBUG_PRINT .and. nq_tr .gt. 0 )
     &      WRITE ( DEBUG_UNIT, '(1X,3A4,3A17)' )
     &        'IC', 'ILA',
     &        'ILB','TRANSPORT','CDERIV_R8(ILA)','CDERIV_R8(ILB)'
$endif

      DO IC = 1, NQ_TR
         ILA = QJET_TRANSP_INDEX(1,IC)
         ILB = QJET_TRANSP_INDEX(2,IC)
         TRANSPORT = CONCI(ILA) * QJET_TRANSP (IC)
         if (ABS(TRANSPORT) .gt. 1.0d300 ) THEN
             WRITE(*,*) ' Transport very large for ',
     &                  'IC, ILA, ILB, CONCI(ILA), QJET_TRANSP(IC)'
             WRITE(*,'(18x,3I5,2G18.12)')
     &                   IC, ILA, ILB, CONCI(ILA), QJET_TRANSP(IC)
!             Pause
         endif
         CDERIV_R8 (ILA) = CDERIV_R8 (ILA) - TRANSPORT
         CDERIV_R8 (ILB) = CDERIV_R8 (ILB) + TRANSPORT

$if defined DEBUG_DERIV_GT_1
             IF (DEBUG_PRINT .and. nq_tr .gt. 0 )
     &      WRITE ( DEBUG_UNIT, '( 1X, 3I4, 3E18.11)' )
     &              IC, ILA, ILB, TRANSPORT,
     &              CDERIV_R8(ILA), CDERIV_R8(ILB)
$endif

      ENDDO
C        NOTE! Effect of outlet itself on mass balance is taken care of
C              in subroutine QSCALC by calls directly from main program.



C   ------------- Add horizontal transports through connections
C                 if horizontal transports are not deactivated:
      IF ( HTR_CALC ) THEN
          Cderiv_ext = 0.0
          EXPORT = 0.0
          DO IC = 1,NC
             ILH = INDXH(IC)+1
C            write(*,*)'nc, ic:', nc, ic
             NLH = INDXH(IC+1)-INDXH(IC)
             NLH_DIM = MAX(1,NLH)
             ILC = INDXC(IC)+1
             IBA = BCONN1(IC)
             ILA = INDXI(IBA)+1
             NLA = INDXI(IBA+1) - INDXI(IBA)
             IBB = ABS(BCONN2(IC))
             IF(BCONN2(IC).GT.0) THEN
C           ........... Between internal basins:
                ILB = INDXI(IBB)+1
                NLB = INDXI(IBB+1) - INDXI(IBB)
C            write(*,*)'MTRAN1(A):ilh,ila,ilb,mlh',ilh,ila,ilb,mlh
                CALL MTRAN1( DEBUG,
     &            NLH_DIM, NLH, HTRANS(ILH), HTR_L(ILH,1),HTR_L(ILH,2),
     &            NLA, NLB, CONCI(ILA), CONCI(ILB), HTRMIX(IC),
     &            NLA, NLB, CDERIV_R8(ILA), CDERIV_R8(ILB),
     &            IMPORT_R8(IBA), IMPORT_R8(IBB)  )
             ELSE
C           ........... Between internal and external basin:
                ILB = INDXE(IBB)+1
                NLB = INDXE(IBB+1)-INDXE(IBB)
                CALL MTRAN1( DEBUG,
     &            NLH_DIM, NLH, HTRANS(ILH), HTR_L(ILH,1),HTR_L(ILH,2),
     &            NLA, NLB, CONCI(ILA), CONCE(ILB), HTRMIX(IC),
     &            NLA, 0, CDERIV_R8(ILA), CDERIV_EXT,
     &            IMPORT_R8(IBA), EXPORT )

C                 CDERIV_EXT, EXPORT not used, are dummy argument
C                 for net flow external basin.
             ENDIF
          END DO
      ENDIF


c     if (Cderiv_ext(1) .ne. 0.0 ) then
c        write(*,*)
c    &     'Error in MASS_TRANSPORT/MTRAN1: CDERIV_EXT(1) is changed'
c     endif

C ==================================================================
C CDERIV_R8 now contains mass increase pr. time unit.
C In the following section, this is converted to time derivatives
C of concentration, assuming constant volume within each layer in time.
C and included in total real*8 derivatives returned to calling program.
C    ( corrections for varying volumes will be made in CNCADJ,
C      called from main program )

$if defined DEBUG_DERIV
            IF (DEBUG_PRINT)
     &      WRITE ( DEBUG_UNIT, '(1X,2A6,2A18)' )
     &              'IB', 'L', 'MASS_INCREASE', 'CDERIV_R8(L)'

$endif

      DO IB = 1,NBI


         DO L = INDXI(IB)+1,INDXI(IB+1)

$if defined DEBUG_DERIV
            MASS_INCREASE = CDERIV_R8(L)
$endif

            CDERIV(L)  = CDERIV(L) + CDERIV_R8(L)/VLAYER(L)  

$if defined DEBUG_DERIV
            IF (DEBUG_PRINT)
     &      WRITE ( DEBUG_UNIT, '(1X,2I6,2E18.11)' )
     &              IB, L, MASS_INCREASE, CDERIV(L)
$endif

         END DO

C     ....... Move updated import values back to external variables:
         IMPORT(IB) = IMPORT_R8(IB)

      END DO
C            NOTE! Surface and bottom exchange not included!


      END Subroutine



C ===========================================================
      SUBROUTINE TR_UPDATE_CONC ( TRCALC, XDEBUG, NAME,
     &     NBI, INDXI, NLI, VLAYER, VTRANS,
     &     TSTEP, LENGTH_TRCF, TRCF, TRCF_RANGE,
     &     NCW, C_W1, C_W2, CDERIV, CONC, CTRLC )
      

C Update water concentrations over specified timestep, using
C time derivatives at input value of concentrations in CDERIV.


      LOGICAL TRCALC
      LOGICAL XDEBUG
      CHARACTER*(*) NAME

      INTEGER NBI
      INTEGER INDXI(NBI+1)
      INTEGER NLI
      real*8 VLAYER(NLI)

      real*8 VTRANS(2,NLI)

      real*8 TSTEP

      INTEGER TRCF_RANGE (0:1, NLI )
      INTEGER LENGTH_TRCF
      real*8  TRCF (LENGTH_TRCF) ! Transfer coefficients

      INTEGER NCW
      real*8 C_W1(NCW), C_W2(NCW)

      real*8 CDERIV(NLI), CONC(NLI), CTRLC(NLI)


C ----------------------- local variables ------------------------
      integer IB, lbase, LSURF, LBOTTOM, LNUM, SOURCE, L, LG1, LG2
      integer TRCF_BASE_INDEX, SOURCE_INDEX, TARGET_INDEX

      integer StepPhase

$if DEBUG_UPDATE >1
      character*5 :: StepPhaseText(2) = (/"First","Last "/)
$endif

      real*8 ACCUM, X, DIFF, V

$if DEBUG_UPDATE > 0 || defined DEBUG_NEG_DIFF
      real*8 CHK_SUM1, CHK_SUM2, CHK_SUM_SCALE
      logical mdebug
      mdebug = xdebug
$endif

$if DEBUG_UPDATE>2
      IF ( NAME.NE.'TEMP' .and. NAME.ne.'OXYG' ) THEN
         DO L= 1, INDXI(NBI+1)  ! All basins:
            if ( CONC(L).lt.0.0 ) THEN
                 MDEBUG = .true.
            ENDIF
         ENDDO
      ENDIF

$endif

      TRCF_BASE_INDEX = 0 ! Now stays constant through all layers,
C                           is updated between each basin

$if DEBUG_UPDATE>0
      IF (MDEBUG) WRITE ( DEBUG_UNIT,'(1x,3A,G14.8)')
     &      '=========== TR_CONC_UPDATE for variable ', NAME,
     &          ' TSTEP=',TSTEP
$endif

C --------------------------------------------------------
      DO IB = 1,NBI  ! Loop through basins
C --------------------------------------------------------

$if DEBUG_UPDATE>0
 10      CONTINUE
         IF (MDEBUG) WRITE ( DEBUG_UNIT,'(1x,A,I5,10(1H=))')
     &          '      ======== Basin:',IB
$endif

         LBASE = INDXI(IB)
         LSURF = LBASE+1
         LBOTTOM = INDXI(IB+1)
         LNUM  = LBOTTOM - LBASE

$if DEBUG_UPDATE>0
         IF (MDEBUG ) THEN
             WRITE ( DEBUG_UNIT,'(4(2X,A,'':'',I5))')
     &          'LBASE', LBASE, 'LSURF', LSURF,
     &          'LBOTTOM', LBOTTOM, 'LNUM', LNUM
             WRITE ( DEBUG_UNIT, '(1X,A5,3A18)' )
     &           'Input: ', 'CONC(L)', 'CDERIV(L)'
             WRITE ( DEBUG_UNIT, '(1X,I5,2G18.12)' )
     &           (L, CONC(L), CDERIV(L),L=LSURF,LBOTTOM)
         ENDIF
$endif


C ------------------------------------------------------
         if ( TRCALC ) then ! Transports are active
C ------------------------------------------------------

            if (NCW.lt.LNUM) STOP 'NCW too short in TR_UPDATE_CONC'

$if DEBUG_UPDATE>0
      CHK_SUM1 = 0.0
      CHK_SUM_SCALE = 0.0
      IF (MDEBUG) write(DEBUG_UNIT,'(/1x:60(''='')/1x,2A)')
     &       ' A. Store concentration and mass of each layer',
     &       ' in double precision arrays:'
$endif

$if DEBUG_UPDATE >1
          if (MDEBUG) WRITE ( DEBUG_UNIT, '(1x,2A6,2A18)')
     &               'L','L-LBASE','conc.(C_W1)','mass (C_W2)'
$endif
            DO L = LSURF,LBOTTOM

               X = CONC(L)
               C_W1 (L-LBASE) = X ! conc. driving exchange of matter
               C_W2 (L-LBASE) = X*VLAYER(L) ! mass to be updated

$if DEBUG_UPDATE>0
               CHK_SUM1 = CHK_SUM1 + X*VLAYER(L)
               CHK_SUM_SCALE = CHK_SUM_SCALE + ABS(X*VLAYER(L))
$endif

$if DEBUG_UPDATE >1
               IF (MDEBUG) WRITE ( DEBUG_UNIT, '(1x,2I6,2(1x,G18.12))')
     &               L, L-LBASE, C_W1(L-LBASE), C_W2(L-LBASE)
$endif

            ENDDO

            StepPhase = 1

 100        CONTINUE

$if DEBUG_UPDATE>0
      IF (MDEBUG) write(DEBUG_UNIT,'(/1x,60(''=''),3(/1x,A),L4,A)')
     &    ' B. Propagate initial concentrations through transfer',
     &    '    coefficients in TRCF, giving diffusive action over',
     &    StepPhaseText(StepPhase),' half of interval TSTEP'
$endif

            SOURCE_INDEX = TRCF_BASE_INDEX + 1 ! valid for first layer

$if DEBUG_UPDATE >1
               IF (MDEBUG) THEN
                  WRITE(DEBUG_UNIT,'(1x,2(3A7,A17)/8x,A7,31X,A7)')
     &              'SOURCE','SOURCE','LG1','C_W1(LG1)',
     &                        'L','TARGET','LG2','TRCF(TARGET_INDEX)',
     &                       '_INDEX','_INDEX'
               ENDIF
$endif
            DO SOURCE = LSURF, LBOTTOM
               LG1 = SOURCE-LBASE

$if DEBUG_UPDATE >1
               IF (MDEBUG) THEN
                  WRITE(DEBUG_UNIT,'(1x,3I7,G17.11)')
     &              SOURCE, SOURCE_INDEX, LG1, C_W1(LG1)
               ENDIF
$endif

               DO L = TRCF_RANGE(0,SOURCE), TRCF_RANGE(1,SOURCE)
                  if (L.eq.0) then
                     CYCLE
                  End if
                  TARGET_INDEX = SOURCE_INDEX + L
                  LG2 = L+LG1

$if DEBUG_UPDATE >1
                  IF (MDEBUG) THEN
                      WRITE(DEBUG_UNIT,'(1x,38(''.''),3I7,G17.11)')
     &                   L, TARGET_INDEX, LG2, TRCF(TARGET_INDEX)
                  ENDIF
$endif

C Fraction of volume transferred in TRCF is now converted
C into exchange of matter:
                  X = 0.5* TRCF(TARGET_INDEX)*(C_W1(LG1)-C_W1(LG2))
                  C_W2(LG2)= C_W2(LG2)+X
                  C_W2(LG1)= C_W2(LG1)-X
C One-way transfer only might give problems with 'decreasing entropy',
C that is; increasing concentration in high-concentration layers.
C The two-way transfer means that each pair of layers are mixed twice,
C once as SOURCE->TARGET, and once as TARGET<-SOURCE. To avoid doubling
C the effective diffusion, a factor 0.5 is used above. C_W2 now
C should contain the resulting amount of substance in each layer
C after applying diffusive exchange specified in TRCF.

$if DEBUG_UPDATE >1
                  IF (MDEBUG) 
     &               WRITE (DEBUG_UNIT,'(3X,3(1x,A,1x,G17.11))')
     &               'x', X, 'c_w2(lg1):', C_W2(LG1),
     &               'c_w2(lg2):', C_W2(LG2)
$endif
               ENDDO

               SOURCE_INDEX = SOURCE_INDEX + LNUM + 1
            ENDDO

            IF ( StepPhase.eq.1 ) THEN

$if DEBUG_UPDATE>0
            IF (MDEBUG) write(DEBUG_UNIT,'(/1x:60(''='')/1x,A)')
     &         'C. add net supply + production for whole interval TSTEP'
$endif

$if DEBUG_UPDATE >1
            IF (MDEBUG)  WRITE ( DEBUG_UNIT, '(1x,A5,4A18)')
     &               ' L=', '  X:', 'C_W2:', 'C_W1:'
$endif
C
                DO L= 1, LNUM
                   X = CDERIV(L+LBASE)*TSTEP*VLAYER(L+LBASE)
                   C_W2(L)=C_W2(L) + X
                   C_W1(L)=C_W2(L)/VLAYER(L+LBASE)

$if DEBUG_UPDATE>0
                   CHK_SUM1 = CHK_SUM1 + X
                      ! Add time change to initial amount
$endif

$if DEBUG_UPDATE >1
                IF (MDEBUG)
     &             WRITE ( DEBUG_UNIT, '(1x,I5,3G18.12)/))')
     &               L, X, C_W2(L), C_W1(L)
$endif

                ENDDO
                StepPhase = 2
                GO TO 100  ! FOR LAST HALF STEP
            ENDIF


C ............. Negative diffusion to compensate diffusive effects
C               of advection. Operates on old concentrations, to
C               let advection work on result of sources and sinks before
C               negative diffusion is applied.
C               Care has been taken to avoid increase of local minima
C               by drawing matter out of a layer.

            DO L = LSURF+1, LBOTTOM ! exchange between L-1 and L
              LG2 = L-LBASE
              V = VTRANS(2,L)*TSTEP ! diffusion: volume exchanged
              if ( V .ge. 0.0 ) CYCLE

$if defined DEBUG_NEG_DIFF
              IF (MDEBUG)
     &           WRITE( DEBUG_UNIT, '(//3(A,I5)/3(2(1x,A,G15.9:)/))')
     &             ' Negative diffusion: IB=',IB,' L=', L,' LG2=',LG2,
     &             'V=VTRANS(2,L)*TSTEP=', V,
     &             'CTRLC(L-1)=',CTRLC(L-1), 'CTRLC(L)=',CTRLC(L),
     &             'C_W2(LG2-1)=',C_W2(LG2-1),'C_W2(LG2)=',C_W2(LG2)
$endif

        ! Concentration difference of controlling substance:
              DIFF = ( CTRLC(L-1) - CTRLC(L))
        ! Set X = limit on allowed transfer to avoid creating minimum:
              if ( DIFF .lt. 0.0 ) THEN ! draw matter from L-1 to L
                 if ( L .le. LSURF+1 ) CYCLE
                 X = (CTRLC(L-1) - CTRLC(L-2))*VLAYER(L-1)
                 if ( X .le. 0.0 ) CYCLE ! >0 if monotonous L-2,L-1,L
              elseif( DIFF .gt. 0.0 ) THEN ! draw matter from L to L-1
                 if ( L .ge. LBOTTOM ) CYCLE
                 X = ( CTRLC(L+1) - CTRLC(L) ) *VLAYER(L)
                 if ( X .ge. 0.0 ) CYCLE ! <0 if monotonous L-1,L,L+1
              else
                 CYCLE
              endif  ! Now DIFF and X has opposite signs.
              if ( abs(V*DIFF) .gt. abs(X) ) then
                  V = X/DIFF  ! V<0 and reduced in abs. value
              endif           ! to avoid creating local minimum

              V = V*( CONC(L-1) - CONC(L) ) ! Substance transport down
              C_W2(LG2-1) = C_W2(LG2-1) - V
              C_W2(LG2)   = C_W2(LG2)   + V

$if defined DEBUG_NEG_DIFF
              IF (MDEBUG)
     &           WRITE( DEBUG_UNIT, '(//3(2(1x,A,G15.9:)/))')
     &             'X=', X, 'DIFF=', DIFF, ' transp. down V=',V,
     &             'after: C_W2(LG2-1)=',C_W2(LG2-1),
     &             'C_W2(LG2)=',C_W2(LG2),
     &             'Final Conc.(L-1)=',C_W2(LG2-1)/VLAYER(L-1),
     &             'Final Conc.(L)=',C_W2(LG2)/VLAYER(L)

$endif

            ENDDO


$if DEBUG_UPDATE>2
            IF ( .not. MDEBUG .AND. NAME.NE.'TEMP'
     &                        .and. NAME.ne.'OXYG') THEN
               DO L= 1, lnum  ! All basins:
                  if ( C_W2(L).lt.0.0 ) THEN
                     MDEBUG = .true.
                     GOTO 10  ! Must only repeat from current basin
                  ENDIF
               ENDDO
            ENDIF
$endif


$if DEBUG_UPDATE >0
         !  Calculate check-sum for mass conservation

            CHK_SUM2 = 0.0
            DO L = LSURF,LBOTTOM
               LG2 = L-LBASE
               CHK_SUM2 = CHK_SUM2 + C_W2(LG2)
               CHK_SUM_SCALE = CHK_SUM_SCALE + ABS( C_W2(LG2) )
            ENDDO

            X = ABS(CHK_SUM1-CHK_SUM2)
            if ( X .gt. 0.5E-6*CHK_SUM_SCALE ) then
               WRITE (DEBUG_UNIT,'('' Check_sum error:'')')
               WRITE (DEBUG_UNIT,'(/2(1X,A,G17.11))')
     &             '***** CHK_SUM1:', CHK_SUM1,
     &             'CHK_SUM2:', CHK_SUM2,
     &             'CHK_SUM_SCALE:', CHK_SUM_SCALE,
     &             'Relative difference:', 2.0*X/CHK_SUM_SCALE


$if DEBUG_UPDATE>2
               IF (.not.MDEBUG) THEN
                  WRITE (DEBUG_UNIT,'('' Repeats:'')')
                  mdebug=.true.
                  goto 10
               endif
$endif
            
            endif
$endif


C       Scale down and store result in CONC:

$if DEBUG_UPDATE >1
                IF (MDEBUG) WRITE (DEBUG_UNIT,
     &              '(1x,A/1X,2a5,2a19, A15)')
     &              'Finally:', 'L', 'LG2',
     &              'C_W2(LG2)', 'C_W1(LG2)',
     &              'CONC(L)'

$endif

            DO L = LSURF,LBOTTOM
               LG2 = L-LBASE
               X = C_W2(LG2)/VLAYER(L)
               CONC(L) = X

$if DEBUG_UPDATE >1
                IF (MDEBUG) WRITE (DEBUG_UNIT,
     &            '(1X,2I5,2G19.12,2G15.8)')
     &            L, LG2, C_W2(LG2), C_W1(LG2),
     &                    CONC(L)
$endif
            ENDDO


C ------------------------------------
         ELSE  ! TRANSPORTS INACTIVATED:
C ------------------------------------


$if DEBUG_UPDATE>0
            IF (MDEBUG) WRITE ( DEBUG_UNIT, '(1X,A,G15.8)' )
     &          'No transports, only local change over TSTEP=',TSTEP
$endif

            DO L= LSURF, LBOTTOM
               ACCUM    = CONC(L) + CDERIV(L)*TSTEP
               CONC(L)  = ACCUM
            ENDDO

         ENDIF

$if DEBUG_UPDATE>0
         IF (MDEBUG ) THEN
             WRITE ( DEBUG_UNIT, '(1X,A10,A18,1x,A18)' )
     &           'Output: L ', 'CONC(L)'
             WRITE ( DEBUG_UNIT, '(1X,I10,1x,G18.12)' )
     &           (L, CONC(L), L=LSURF,LBOTTOM)
         ENDIF
$endif

C -----------------------------------------------------
         TRCF_BASE_INDEX = TRCF_BASE_INDEX + LNUM*LNUM
      ENDDO  ! Next basin
C -----------------------------------------------------

      RETURN
      END subroutine
      
      End Module
      
