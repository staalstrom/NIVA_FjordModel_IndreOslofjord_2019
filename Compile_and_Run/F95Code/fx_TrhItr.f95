! ===================================================================
      Module fx_TRHITR

      use fx_RunControl, only: DEBUG_UNIT
      use ModelDimensions
      use fx_EigJacob
      use fx_WindCurr
      
      implicit none

! ===================================================================
! Eutrophication model   - Fixed-format file:    fx_TRHITR.f95

$define ALLOCATE_ARRAYS  
			! use dimensioning parameters instead (ABSOFT)


C --------- to try both:
$undefine test_stabilize

C ------------  Compile variables which control debug printout:

$define TEST_ITR 0

		!  controls monitoring of iteration process in ITERATE_SURFACE
		!  on unit nf ( = free unit number from 901 and up )
		
		!   =0: no debug prints

		!   =1: Monitor iteration process only if it did not converge to
		!       specified accuracy, (computing continues with results achieved)
		
		!  >=2: Monitor iterations also if DEBUG=.TRUE., even if they converge
		!       Reports iteration counts on screen.
		
		!  >=3: As 2, but does not close output file between each step,
		!       uses same debug_unit as TRANSP_x Modules


$undefine TEST_GAUSS
C     activates output from GAUSS_ELIM on unit nf,
C     when iteration is monitored, as controlled above.

$undefine TEST_TRHFNC
C     activates dump of info from TRHFNC_CONN on unit nf
C     when iteration is monitored.


$if defined TEST_ITR || DEFINED TEST_GAUSS || DEFINED TEST_TRHFNC
      integer :: NF = 1
      logical DEBUG
$endif

C --------------------- local variables ------------------------


      INTEGER, parameter, private  :: M =dimMBI

      real*8, private  ::   AR(M), Z(M), DZ(M)
      real*8, private  ::   VD(M), DVD_DZ(M,M), DVD_DTX(M)
      real*8, private  ::   VD_INIT(M), VD_TARGET(M), DVD_INF_INIT(M)
      real*8, private  ::   VD_PREV(M), VDDEV(M)
      real*8, private  ::   QARR(M,M+1), WJ(M,M)    ! Work arrays
      real*8, private  ::   VD_REDUCTION(M)
      real*8, private  ::   EIGVMAT(M,M), EIGVEC(M,M), ACCUR
      real*8, private  ::   COEFF(M)
      real*8, private  ::   EigvTstep
      real*8, private  ::   X, Y

C                 Relative error in DZDT permitted:
      real*8, private, parameter  ::   WABS_CONV =1.D-2
      real*8, private ::   W_MAX, EPS
      real*8, private, parameter  ::   EPSILON8 = 1.0D-14

      INTEGER, private  ::  FNC_COUNT, I, K, IA, FNC_PHASE1
      INTEGER, private, parameter ::  FNC_MAX = 300    ! Typically 10-15 iterations required

C Logical loop control variables:
      LOGICAL, private  ::  CONVERGED, CHANGED, SLOW_CHANGE, PHASE2
C       .FALSE. PHASE 1: Iteration to equilibrium values.
C       .TRUE.  PHASE 2: Iteration to final values to be used for
C                        transport calculations (exported through ZSURF)

      LOGICAL, private  ::  EIGVEC_INITIATED

      real*8, private  ::  DZ_MAX

C  Time constant for adjusting transports:
      real*8, private, parameter  ::  TSTEP_LOW=0.0005
      real*8, private  :: TSTEP_CONV
C           TSTEP_LOW could be imported and controlled externally.

$if defined  test_stabilize
      logical, private  ::  STABILIZE_TARGET
$endif


			! Reinitiated by entry point TRHITR_INIT
			! called at start of simulation, to get identical runs

      real*8, private :: DZDT_SUM =0.0
      real*8, private :: DZDT_WGT =0.0 


         ! Error tolerance scale:
      real*8, private :: DZDT_NORM       ! in order to get identical runs
         ! 48 hour running mean of DZ_DT, used to scale accuracy tolerance in

C ====================================================================
C Internally allocated or local arrays:
      INTEGER, private :: NBI_ALLOC =0,
     &         NBE_ALLOC = 0, MLC_ALLOC = 0, NC_ALLOC = 0

$if defined  ALLOCATE_ARRAYS

!         - dimensions, reset to values>0 on allocation of arrays:
      real*8, private , Allocatable ::  DDZERO  (:)
      real*8, private , Allocatable ::    TRANSP_WIDTH        (:)
      real*8, private , Allocatable ::    HLAYER              (:)
      real*8, private , Allocatable ::    TRANSP_CONSTANT     (:)
      real*8, private , Allocatable ::    VTOTD_LAND_SURFACE  (:)
      real*8, private , Allocatable ::    ZSURF_EXT           (:)
      real*8, private , Allocatable ::    DZDT_EXT            (:)
      real*8, private , Allocatable ::    WIND_CURRENT        (:)
      real*8, private , Allocatable ::    I_MAX_WINDC         (:)

      INTEGER, private , Allocatable ::   ILC_V               (:)
      INTEGER, private , Allocatable ::    NL_V               (:)
      INTEGER, private , Allocatable ::   IBA_V               (:)
      INTEGER, private , Allocatable ::   ILA_V               (:)
      INTEGER, private , Allocatable ::   IBB_V               (:)
      INTEGER, private , Allocatable ::   ILB_V               (:)
      Integer, private ::  Alloc_chk=-1
      SAVE Alloc_chk
$else
      real*8, private ::  DDZERO              (dimMLC)
      real*8, private ::  TRANSP_WIDTH        (dimMLC)
      real*8, private ::  HLAYER              (dimMLC)
      real*8, private ::  TRANSP_CONSTANT     (dimMC)
      real*8, private ::  VTOTD_LAND_SURFACE  (dimMBI)
      real*8, private ::  ZSURF_EXT           (dimMBE)
      real*8, private ::  DZDT_EXT            (dimMBE)
      real*8, private ::  WIND_CURRENT        (dimMBI)
      real*8, private ::  I_MAX_WINDC         (dimMBI)
      INTEGER, private :: ILC_V               (dimMBI)
      INTEGER, private ::  NL_V               (dimMBI)
      INTEGER, private :: IBA_V               (dimMBI)
      INTEGER, private :: ILA_V               (dimMBI)
      INTEGER, private :: IBB_V               (dimMBI)
      INTEGER, private :: ILB_V               (dimMBI)
$endif


      contains


C --------------------------------------------------------------------
C For each time step, surface levels to be used in transport calculation
C will be adjusted relative to time integrated values so that the
C further integration will avoid "overshooting" equilibrium values.
C If external level variations are damped in internal basins, this has
C little effect, but it is important for cases with fast reponse and
C small damping of level fluctuations.  In that case it permits much
C larger timesteps than given by the surface dynamics.

C First the subroutine iterates to equilibrium surface levels and rate
C of change, i.e levels giving constant transport values in time.

C Then a weighted mean of these transport values and the values given by
C the imported "actual" surface levels are calculated.  The weighting is
C done to avoid instability under the imported maximum limit on the
C timestep: If the "actual" change rates would integrate to levels
C beyond equilibrium values, the change rate is adjusted so that the
C integration of transports over time will lead the solution in
C in direction of the (linearly changing) equilibrium.  The adjusted
C rates thus represents anticipated mean values over such a timestep.

C This is handled by the following subroutines,
C called from TRANSP_SETUP in TRANSP_2.FOR:

C  Logical function
C      ITERATE_SURFACE:  Called from TRANSP_SETUP after TRHITR_PREP,
C                        performs iteration, returns .TRUE. if OK

C  Subroutine
C      TRHITR_PREP:      Called from TRANSP_SETUP, prepares for iteration.
C                        Calls TRHITR_CONN for each connection.

C  Subroutine TRHFNC:
C                        Called from ITERATE_SURFACE,
C                        evaluating VD as function of Z.
C                        Calls TRHFNC_CONN for each connection.


C =====================================================================



      LOGICAL FUNCTION ITERATE_SURFACE
     &       (DEBUG_ACTIVE, 
     &        ITRZ, TSTEP, NBI, ZSURF, INDXI, NLI, AREA, DZ_DT,
     &        NLC, UFLOW )
      
C In:
      LOGICAL DEBUG_ACTIVE ! .TRUE. activates debug printout
      LOGICAL ITRZ  ! .TRUE. activates iteration

      real*8    TSTEP ! Upper limit to timestep, controls iteration.
      INTEGER NBI          ! Number of basins
      INTEGER INDXI(NBI+1) ! INDXI(I)+1 is index to area OF BASIN I
      integer NLI          ! Dimension of AREA
      real*8 AREA (NLI)      ! Surface areas (m2)
C Out:
      real*8 ZSURF (NBI)   ! Surface levels (in:guess/out:results)
      real*8 DZ_DT (NBI)   ! Iterated change rates for levels (m/s)
      integer NLC
      real*8  UFLOW (NLC)
C              - flow velocity, exported as user info. to main model

      logical ALREADY_OPENED
      
C -------------------- Diagnostic output control:

$if TEST_ITR == 1
      DEBUG = .FALSE.
$else
      DEBUG = DEBUG_ACTIVE
$endif

$if defined  test_stabilize
      STABILIZE_TARGET = .false.
$endif


C -------------------------------------------------------
C                PREPARE FOR ITERATION:
C -------------------------------------------------------

      EPS = 3.0*EPSILON8 ! Lowest relative sign. change in z
C           Used to avoid infinite loop due to too strict accuray req.

      TSTEP_CONV = MAX ( TSTEP_LOW, TSTEP)*24.*3600.  ! days -> seconds
C           Time horizon for projected mean transports in future step

C ------- Check dimensions of work-array:
      IF (NBI .GT. M ) THEN
         STOP 'NBI > M in subroutine TRHITR2'
      ENDIF

  10  CONTINUE

$if defined TEST_ITR
      call print_unit  ! Find unused unit number in nf (otherwise =0)
      if (DEBUG.and.nf.gt.0 ) then
         write(*,*) 'Iteration attempt recorded on file unit ', nf
         if (nf.gt.0) then
             INQUIRE(UNIT = nf, OPENED = ALREADY_OPENED )
             if (.not.ALREADY_OPENED) then
                open (nf,FILE='ITR_SURF.LST', ACCESS='SEQUENTIAL' )
                write(*,*) 'to special file ITR_SURF.LST'
             endif
         endif
         write(nf,*) ' >>>>>>>>>> STARTS ITERATE_SURFACES <<<<<<<<<<'
         write(nf,*) '           TSTEP =',TSTEP
      endif
$endif

C ------- Transfer surface levels and areas to work arrays:
      DO I =1,NBI
         Z(I)  = ZSURF(I)
         DZ(I) = 0.0
         IA = INDXI(I)+1
         AR(I) = AREA(IA)
         VD_REDUCTION(I) = 1.0
      ENDDO


C ----------- Prepare LOOP control variables:
      FNC_COUNT = 0
      CONVERGED = .not.ITRZ
          ! Considered converged already if iteration is turned off,
          ! will then only make one evaluation and exit the loop

      SLOW_CHANGE = .true.
      EIGVEC_INITIATED = .false.
      PHASE2 = .FALSE.



C --------------------------------------------------------------------
C * * * * * * * * * * * ITERATION LOOP: * * * * * * * * * * * * * * *
C --------------------------------------------------------------------

      DO WHILE ( FNC_COUNT .LT. FNC_MAX )


$if defined TEST_ITR
         if (DEBUG) then
            write ( nf,'(/'' --------- Iteration loop ----------'')')
            write ( nf,*) ' CONVERGED, SLOW_CHANGE, PHASE2, FNC_COUNT:'
            write( nf, '(5X,3L10,I5)')
     &                      CONVERGED, SLOW_CHANGE, PHASE2, FNC_COUNT
         endif
$endif


C --------------------------------------------------------
C * * * GET TRANSPORT DESCRIPTION FOR GIVEN Z VALUES * * *
C --------------------------------------------------------

C   Update minimum Z-scale for VD derivative in TRHFNC:
C   Max. 1 meter scale, otherwise divergence may be accelerated,
C   since a large scale gives smaller derivatives,
C   and larger fluctuations.
         IF (FNC_COUNT .GT.0) THEN
            DZ_MAX = 0.0
            DO I = 1, NBI
               IF( .NOT. phase2 ) THEN
                  DZ_MAX = MAX  ( DZ_MAX, ABS( DZ(I) ) )
               ELSE
                  DZ_MAX = MAX  ( DZ_MAX, DBLE(ABS( ZSURF(I)-Z(I)) ) )
               ENDIF  ! More flat response in 2. phase (mean value)
            ENDDO
         ENDIF
C   (stores initial value to be used in next call to ITERATE_SURFACES)

$if defined TEST_ITR
         if (DEBUG) then
            write ( nf,*)
            write ( nf,*) ' -------- calls TRHFNC'
         endif
$endif


C -----------------------------------------------------------
C Store previous VD values for calculating correction factor
C to off_target measure below:
         DO I=1,M
            VD_PREV(I) = VD(I)
         ENDDO

C then get new transports (VD), Jacobian dVD/dZ, and
C time derivative DVD_DTX due to change rate of external Z's
C for last set of internal Z values.

         CALL TRHFNC ( Z, VD, DVD_DZ, DVD_DTX, NBI,
     &                 NLC, UFLOW )
         FNC_COUNT = FNC_COUNT + 1

! ----------------------------------------------------------
         if (CONVERGED) EXIT   ! After final call
! ----------------------------------------------------------

$if defined TEST_ITR
         if (DEBUG) then
            write ( nf, '(''  Z= '',10(3E17.10:/))' ) (Z(I),I=1,NBI)
            write ( nf, '('' VD= '',10(3E17.10:/))' ) (VD(I),I=1,NBI)
            write ( nf,*) 'DVD_DZ(I,K):'
            do I = 1, NBI
                write(nf, '('' I='',I2,'':'')' ) I
                write(nf, '(5x,3('' K='',I3,'':'',G16.10:))' )
     &                           ( K, DVD_DZ(I,K), K=1,NBI)
            enddo
            write ( nf, '('' DVD_DTX '',10(3E16.10:/))' )
     &                     (DVD_DTX(I),I=1,NBI)
         endif
$endif


C -------------------------------------------------------------------
C  Save/evaluate results:

         if ( FNC_COUNT .eq.1 ) THEN
C        ------ First pass: Store volume changes for input values of ZSURF:
            DO I=1,NBI
                VD_INIT(I) = VD(I)
            ENDDO

         ELSE
C        -------------- later passes:

$if defined TEST_ITR
            if (DEBUG) then
               WRITE (NF,*) ' CHECKS CONVERGENCE:'
               Write(nf,'(1x,A5, 3A17)') 'I',' VD','VD_TARGET','AR'
               Write(nf,'(1x,I5, 3E17.10)')
     &                (I,VD(I),VD_TARGET(I),AR(I),I=1,NBI)
            endif
$endif

C        ....... a. Check convergence:
            W_max = 0.0
            DO I =1,NBI
               W_MAX = MAX ( W_MAX, ABS(VD(I)-VD_TARGET(I)) / AR(I) )
            ENDDO
            IF (W_MAX .LT. WABS_CONV*DZDT_NORM ) THEN
                CONVERGED = .TRUE.
            ENDIF

$if defined TEST_ITR
            if (DEBUG) THEN
                 WRITE (NF,*) ' CHECKED CONVERGENCE:'
                 WRITE (NF,'(3(1X,A,'':'',G16.10)/A,L10)' )
     &              'w_max',     w_max,  'wabs_conv', wabs_conv,
     &              'dzdt_norm', dzdt_norm, 'CONVERGED=',CONVERGED
            ENDIF
$endif

            IF( CONVERGED ) THEN
                IF ( SLOW_CHANGE .OR. PHASE2 ) CYCLE  ! Finished
C               .... otherwise: PHASE 1 converged,
C                    initialize phase 2 iteration:
                DO I = 1,NBI
                   DVD_INF_INIT(I) = (VD(I) - VD_INIT(I))/SQRT(AR(I))
                ENDDO
C      Rescales transport values to get symmetric eigenvalue problem,
C      see below.
                PHASE2 = .TRUE.
                CONVERGED = .FALSE.
                FNC_PHASE1 = FNC_COUNT

$if defined TEST_ITR
                if (DEBUG) THEN
        WRITE (NF,'(/'' ############## Switch to phase 2 iteration'')')
        WRITE (NF,'('' DVD_INF_INIT('',I3,'')/SQRT(AR): '',E17.10)' )
     &                 ( I, DVD_INF_INIT(I),I = 1,NBI)
                ENDIF
$endif
            ENDIF


C         ...... b. Look at difference between intended and achieved change,
C                   use difference to update a reduction factor applied to
C                   VD deviation in next step to hit target better:
C                    corr = (intended - achieved)*specified/achieved
C                         = (intended/achieved-1)*specified
            DO I = 1, NBI
               X = VD(I)        - VD_PREV(I)    ! Achieved change
               Y = VD_TARGET(I) - VD_PREV(I)    ! Intended change
               if (Y*X .gt. 0.0 ) THEN
                   VD_REDUCTION(I) =  VD_REDUCTION(I)*Y/X ! Upper limit?
               else
                   VD_REDUCTION(I) = 1.0
               ENDIF
            ENDDO
$if defined TEST_ITR
            if (DEBUG) THEN
                Write ( nf, '('' vd_reduction = '',10(3E17.10:/))' )
     &                (vd_reduction(I),I=1,NBI)
            ENDIF
$endif

         ENDIF


C ---------------------------------------------------------------------
C   * * * GET TARGET VALUES FOR VD BY DIFFERENT STRATEGIES: * * *
C ---------------------------------------------------------------------

         IF (SLOW_CHANGE) THEN
C        ------- Initially, assume small change in VD throughout step,
C                set target value equal to projected 1.order mean:
            DO I=1,NBI
               X = DVD_DTX(I)
               DO K=1,NBI
                  X = X + DVD_DZ(I,K)*VD(K)/AR(K)
               ENDDO
               X = X*TSTEP_CONV
               VDDEV(I) = X
               VD_TARGET(I) = VD_INIT(I) + X/2.0
               IF ( X .GT. SQRT(WABS_CONV)*VD_TARGET(I) ) THEN
                   SLOW_CHANGE=.FALSE.   ! leaves this mode
               ENDIF
            ENDDO

C                Estimate target accuracy by inconsistency
C                in first order assumption VD= A + B*t :
            DO I=1,NBI
               X = 0.0
               DO K=1,NBI
                  X = X + DVD_DZ(I,K)*VDDEV(K)/AR(K)
               ENDDO
               X = X*(TSTEP_CONV/2.0)
               IF ( ABS(X)/AR(I) .GT. WABS_CONV*DZDT_NORM ) THEN
                   SLOW_CHANGE=.FALSE.   ! leaves this mode
               ENDIF
            ENDDO
         ENDIF
         IF ( SLOW_CHANGE ) GOTO 300  ! to find DZ, below


         IF ( .NOT. PHASE2 ) THEN
C -------------------------------------------------------------------
C                        Phase 1 iteration:
C -------------------------------------------------------------------

C    Find values DZ_DT(I) giving constant volume change rate
C    by balancing the changes due to time derivatives of
C    external surface levels, i.e. solving
C    linear equation system:
C                DVD_DZ(I,K)*DZ_DT(K) = - DVD_DTX(I)
C    From TRHFNC it is known that diagonal elements of DVD_DZ
C    are at least as large in absolute value as other elements
C    in the same row. The equation systems can then be solved by
C    simple Gauss elimination, without pivoting of rows.
C    It is presupposed that all basins are connected, directly
C    or indirectly, and that at least one basin is connected to
C    an external area with specified surface levels.
C    If none of the basins are connected to external sources,
C    there are in general no solution, and the model is stopped
C    with a warning.

            CALL GAUSS_ELIM( .TRUE.,
     &                       M, NBI, DVD_DZ, DVD_DTX, DZ_DT, QARR)
C           1.argument .TRUE. activates triangularization of matrix.
C         Triangularized matrix, vector and triangulation coefficents
C         are stored in QARR, for use below when solving for the same
C         equations with another constant vector instead of DVD_DTX.

C        ----- Convert DZ_DT(I) into target volume change rates VD
            DO I=1,NBI
               VD_TARGET(I) = DZ_DT(I) * AR(I)
            ENDDO

$if defined TEST_ITR
            if (DEBUG) then
               Write( nf, *)
               Write( nf, *)' --------- after GAUSS_ELIM: DZ_DT= '
               Write( nf, '(5x,4E17.10)' ) (DZ_DT(I),I=1,NBI)
            endif
$endif



         ELSE
C ----------------------------------------------------------------------
C                     Phase 2 iteration:
C ----------------------------------------------------------------------

C     Finds analytical, transient part of solution
C     for the change from initial current time values
C     VD_INIT, towards the established equilibrium values.
C     Uses rescaled changes DVD_INF_INIT set up above at transition
C     from PHASE 1 to PHASE 2, and iterated values of DVD_DZ,
C     that is, solve the set of equations given by:
C           A(i,k)*X(k) = lambda(i)*X(i)
C     with:
C           A(i,k) = DVD_DZ(I,K)/SQRT(Ar(k)*Ar(i))
C             X(k) = X(k)[t=0]*exp(lambda(k)*t)
C             X(k)[t=0] = VD(k)[t=0]/sqrt(Ar(k) = dZk/dt*sqrt(Ar(k))
C     which is a rescaled version of the equation:
C           DVD_DZ(I,K)*dZk/dt = dVDi/dt

C     .............................................
C     1. Get Eigenvalues and Eigenvector of DVD_DZ:
C     .............................................
            do I=1,nbi
               do K=1,nbi
                  EigvMat(i,k) = DVD_DZ(i,k)/SQRT(Ar(K)*Ar(i)) !(1/sec)
               enddo
            enddo
            CALL JACOBI( EIGVEC_INITIATED,
     &           M, NBI, 20, 1.D-9, EIGVMAT, EIGVEC, WJ, I, K, ACCUR)
C             Returns Eigenvalues in diagonal of EIGVMAT (rotated),
C                = exponential specific adaption rates
C             and eigenvectors by columns in Eigvec.
$if defined TEST_ITR
            if (DEBUG) then
              WRITE (NF,*) ' ----- has called JACOBI'
              WRITE (NF,'(1x,I5,'' iterations, error status:'',I5)') I,K
              WRITE (NF,'('' Accuracy:'',G12.5)') Accur
              DO I=1,NBI
                 DO K=1,NBI
                    Write(nf,'('' '',I3,'',''I3,'' : '',2(A,E17.10))' )
     &                  I,K, 'EIGVMAT', EIGVMAT(I,K),
     &                       'EIGVEC',  EIGVEC (I,K)
                 ENDDO
              ENDDO
            ENDIF
            If ( Accur .gt. 2.0E-7) then
               write(*,"(2A, G12.3,A/2A)")' Warning: Inaccurate eigenvalue solution ',
     &              '(Accur=',Accur, ',required to be <=2.0E-8)',
     &              'in subroutine ITERATE_SURFACE in Module fx_TrhItr ',
     &              'called from TRANSP_SETUP in Module fx_TRANSP_2.'
!               Pause 'Press Enter to continue'
            Endif
$endif


C     ..............................................................
C     2. Find transient solution by solving linear equation system
C        for COEFF(k):   Transients (t=0) + [(VD(t=ì) - VD(t=0))]= 0
C           Transients at t=0:  SUM k ( Eigvec(i,k)*Coeff(k) )
C           Constant vector:    [ ]= DVD_INF_INIT
C        ::> Sum of transients = VD_INIT - VD_INF at t=0,
C                                -->0 for t--> inf
C     ..............................................................
          CALL GAUSS_ELIM(
     &               .true., M, NBI, EIGVEC, DVD_INF_INIT, COEFF, QARR )
$if defined TEST_ITR
            if (DEBUG) then
               Write(nf,*)' COEFF from GAUSS_ELIM: '
               Write(nf,*) (COEFF(I),I=1,NBI)
               WRITE (NF,'('' given DVD_INF_INIT('',I3,''): '',E17.10)')
     &                 ( I, DVD_INF_INIT(I),I = 1,NBI)
            endif
$endif

C     ..............................................................
C     3. Calculate target values for VD over timestep TSTEP_CONV
C        as means of analytical solution to VD(i) over t:
C               VD(i) + SUM [Eigvec(i,k)*Coeff(k)*exp(Eigval(k)*t)]
C                        k
C        or:
C                                                exp(Eigval(k)*t)-1
C         mean= VD(i) + SUM [Eigvec(i,k)*Coeff(k)* ------------------ ]
C                        k                           Eigval(k)*t
C     ..............................................................


C        ..... Transform COEFF(K) to Coeff*mean of exponential term:
            Do K = 1, NBI
               EigvTstep = EigvMat(K,K)*TSTEP_CONV
               if (EigvTstep.ne.0) then
                  Coeff (K) = Coeff(K)* ( exp(EigvTstep)-1.0 )/EigvTstep
               endif  ! if Eigvstep -->0, expression -->1
            END DO
$if defined TEST_ITR
            if (DEBUG) then
               Write(nf,*)' COEFF transformed to mean value factors'
               Write(nf,*) (COEFF(I),I=1,NBI)
            endif
$endif

C        ..... Add terms to get mean of transient solution over TSTEP
            DO I = 1, NBI
               X = DVD_INF_INIT(I)  ! Rescaled change in transport
               DO K = 1, NBI
                  X = X + Eigvec(I,K)*Coeff(K)
               END DO  ! + mean of [DVD_inf_init]*exp(-lambda*t)
               X =  VD_INIT(I) + X*SQRT(Ar(I))
$if defined  test_stabilize
               if (STABILIZE_TARGET) then
$endif
               if (  FNC_COUNT .gt. FNC_PHASE1+6 ) THEN
                   X = ( X + VD_TARGET(I) )/2.0D0
C                     ..... After iteration 6 in PHASE2:
C                     Updates VD_TARGET with 50% of estimated change,
C                     to stabilize against observed fluctuations in
C                     VD_TARGET, which otherwise can lead to
C                     Z fluctuating between 2 sets of values,
C                     and not converging.
               endif
$if defined  test_stabilize
               endif
$endif
               VD_TARGET(I) = X
            ENDDO

C        ..... First time in PHASE 2:
C              find the best starting point for the new iteration:
            IF ( .not. EIGVEC_INITIATED ) THEN
               Y = 0.0
               X = 0.0
               DO I = 1, NBI
                  IF ( VD_TARGET(I) .ne. 0.0 ) THEN
                      X = MAX( X, ABS(VD_INIT(I)/VD_TARGET(I)) )
                  ENDIF
                  IF ( VD(I) .NE. 0.0 ) THEN
                      Y = MAX (Y, ABS(VD_TARGET(I)/VD(I)) )
                  ENDIF
               ENDDO
               X = (1.0-MIN(1.0D0,X*X))*MIN(1.0D0,Y*Y)
C                       --> 0 if X approx. 1 or if Y<<1
C                       --> Y*Y if X <<1
               Y = min(1.0d0, Y**2)  ! Interpolating factor
               DO I = 1, NBI
                  Z(I) = ( X*Z(I) + (1.0-X)*ZSURF(I) )
               ENDDO

$if defined TEST_ITR
        if (DEBUG) then
           write (nf,* )'---- Interpolated as start values of Z : ',
     &            ' for phase 2 iteration:'
           write (nf, '(20(3E20.11:/4x))' ) ( Z (I),i=1,nbi )
        endif
$endif


C           ..... Set status for repeated passes
C                 of the jacobi transient calculation:
               EIGVEC_INITIATED = .True.
C                 will use EIGVEC as initial rotation matrix
C                 in next pass of loop.
               CYCLE  ! Start with finding VD for this Z value.
            ENDIF

         ENDIF



C ------------------------------------------------------------------
C  * * *   FIND NEW APPROX. Z+DZ GIVING VD CLOSER TO TARGET   * * *
C ------------------------------------------------------------------

  300    CONTINUE


$if defined TEST_ITR
            if (DEBUG) then
               Write(nf,*)' VD_TARGET: '
               Write(nf,'(1x4E17.10)') (VD_TARGET(I), I=1,NBI)
            endif
$endif


C    .............................................................
C    1.  Find surface level adjustments DZ(I,..) necessary to get
C        transport values VD_TARGET, i.e. solve system of linear
C        equations:  DVD_DZ(I,K)*DZ(I) = -(VD(I)-VD_TARGET(I))
C    .............................................................


         DO I=1,NBI
            VDDEV(I) = ( VD(I) - VD_TARGET(I) ) * VD_REDUCTION(I)
         ENDDO

         CALL GAUSS_ELIM( SLOW_CHANGE .OR. PHASE2,
     &                    M, NBI, DVD_DZ, VDDEV, DZ(1), QARR)
C         1. ARG. .TRUE.: 1. time DVD_DZ is used, triangularize.
C                 .FALSE.: ( phase 1 in two phase iteration )
C                    DVD_DZ used in GAUSS_ELIM before,
C                    triangularization of matrix already done.
C                    Coefficients QARR which were set then
C                    are used again on new vector VDDEV.

$if defined TEST_ITR
      if (DEBUG) then
         write (nf, '('' VDDEV:'',20(3E20.11:/4x))' )
     &                       ( VDDEV (I), i=1,nbi )
         write (nf, '('' DZ : '',20(3E20.11:/4x))' )
     &                       ( DZ (I),i=1,nbi )
      endif
$endif

C     .........................................................
C     2. Update Z, store new values as new set:
C        Stops if change is close to insignificant,
C        to avoid 'infinite' loop.
C     .........................................................

         changed = .false.
         DO I = 1, NBI
            Z(I)  = Z(I) + DZ(I)
            if ( abs(DZ(I)) .ge. abs(Z(I))*EPS ) THEN
                 changed =  .true.
            endif
         ENDDO
         converged = .not. changed


$if defined TEST_ITR
         if (DEBUG) then
            write (nf, '('' --> Z : '',20(3E20.11:/4x))' )
     &                       ( Z (I),i=1,nbi )
         endif
$endif


      END DO   ! (while .NOT. CONVERGED)


C ------------------------------------------------------------
C * * * * * * * *  END OF ITERATION LOOP  * * * * * * * * * *
C ------------------------------------------------------------


$if TEST_ITR >= 2
        write(*,*)'           FINAL FNC_COUNT = ',FNC_COUNT
$endif


$if defined  test_stabilize
$if defined TEST_ITR
      if ( DEBUG .AND. .NOT. STABILIZE_TARGET ) then
           stabilize_target = .true.
           write(nf,*) '>>>> Repeat with target stabilization'
           goto 10
      endif
$endif
$endif


C --------------  Check results, perform diagnostics if specified:

      if( .not. converged ) then
         write(*,*)' FNC_COUNT  ',FNC_COUNT,'     FNC_MAX=',FNC_MAX
         write(*,*)'Surface level iteration did not converge'
$if defined TEST_ITR
         if( .not. DEBUG ) then
             DEBUG = .true.
             write(nf,*) 'Surface level iteration did not converge'
             write(nf,*) 'Repeats with diagnostics:'
             go to 10
         endif
$if defined  test_stabilize
         if ( .NOT. STABILIZE_TARGET ) then
             stabilize_target = .true.
             write(nf,*) '  ####### Repeats with target stabilization'
             goto 10
         endif
$endif
$endif

      endif


$if defined TEST_ITR
      if ( DEBUG ) then
        if (ITRZ) then
         write(nf,*) '>>>>>>> Iteration process finished  <<<<<<<'
        else
         write(nf,*) '>>>>>>>  No iteration: ITRZ=.false. <<<<<<<'
        endif
      endif
$endif

$if TEST_ITR >=1 && TEST_ITR <=2
      if ( DEBUG .and. nf.gt.0 ) then
             CLOSE(nf)
             nf=-1
!             PAUSE 'Inspect surface iteration or continue with ENTER'
      endif
$endif


C -------------------------- Export result:
      DO I =1,NBI
         ZSURF(I) = Z(I)
         DZ_DT(I) = VD(I)/AR(I)
      ENDDO
      ITERATE_SURFACE = CONVERGED

      END FUNCTION



C ===============================================================
      SUBROUTINE GAUSS_ELIM ( initial, M, N, A, B, X, QARR )
      
   
      LOGICAL Initial
      INTEGER M,N
      real*8 A(M,M), B(M), X(M), QARR(M,M+1)

C  Solves linear system A*X + B = 0 by simple Gauss elimination:
C  No pivoting, and no check on singularity of matrix is done, since
C  the model is assumed to avoid this, see comment at subroutine call.

      real*8 Q
      INTEGER I,K,L

$if defined TEST_GAUSS

      if (debug) then
         WRITE (nf,*) ' -----------  GAUSS-ELIMINATION: ',
     &                ' N=',N,' INITIAL=', initial
         DO I = 1, N
            DO K = 1, N
               WRITE (nf,'('' A('',I3,'',''I3,''):'', E17.10)') 
     &          I, K, A(I,K)
            ENDDO
            WRITE (nf, '(20x,''B('',I3,'')='', E17.10)' ) I, B(I)
         ENDDO
      endif
$endif


      IF ( Initial ) THEN
C     -------------------- First call for new matrix A:


C     --------- Store in work array:

         DO I = 1,N
             DO K = 1,N
                 QARR (I,K) = A(I,K)
             ENDDO
         ENDDO

C     --------- Triangularize matrix, and store coefficients:

         DO I = 1,N-1
             DO L = I+1,N

C        ....... Eliminate column I from eq. L=I+1,..., N
C                by adding equation I multiplied by suitable coeff. Q:
                Q = - QARR(L,I)/QARR(I,I)
                DO K = I+1,N
                   QARR(L,K) = QARR(L,K) + Q*QARR(I,K)
                END DO ! Transformed A stored in QARR(L,K), L>I, K>I
C                        Final transformed values for eq. L=I+1  (L<=K).

C        ....... Store coefficients in unused lower triangular elements
                QARR(L,I)= Q      ! (L>I)
             ENDDO
         ENDDO

      ENDIF


C    --------------------- Always, store transformed B in last QARR col.
      DO L = 1,N
          Q = B(L)
          DO I = 1, L-1
             Q = Q + QARR(L,I)*QARR(I,N+1)
          ENDDO
          QARR (L,N+1 ) = Q  ! transformed B(I)
      ENDDO

$if defined TEST_GAUSS
      if (debug) then
         WRITE (nf, '(/''   Transformed system in QARR:'')')
         DO I = 1, N
            DO K = 1, I-1
               WRITE (nf,'('' Coeff QARR('',I3,'',''I3,''):'', E17.10)')
     &                        I, K, QARR(I,K)
            ENDDO
            DO K = I,N
               WRITE (nf,'('' A~ IN QARR('',I3,'',''I3,''):'', E17.10)')
     &                        I, K, QARR(I,K)
            ENDDO
            WRITE (nf,'('' B~ in QARR('',I3,'',N+1):'', E17.10)')
     &                        I, QARR(I,N+1)
         ENDDO
      endif
$endif

C  -------- Backsolve for X:
      DO I = N,1,-1
         Q = -QARR(I,N+1)
         DO K = N,I+1,-1
             Q = Q - QARR(I,K)*X(K)
         END DO
         if( QARR(I,I) .ne. 0 ) then
             X(I) = Q/QARR(I,I)
$if defined TEST_GAUSS
                  if (debug) then
               WRITE(nf,'('' I='',I3,'' Q='',E17.10,'' X('',
     &                    I3,'')='',E17.10)')
     &                           I,         Q,         I,  X(I)
                  endif
$endif
         elseif (initial) then
            write(*,*) 'Error in ITERATE_SURFACE:',
     &          'Gauss-elimation =>  diagonal',I,' = 0'
         endif
      ENDDO
      END SUBROUTINE



C ================================================================

      SUBROUTINE TRHITR_PREP (  DEPTH, NBI,  QINPUT, XMIX,
     &          INDXI, NLI, DENSI,
     &          NBE, INDXE, NLE, DENSE, ZSURFE, DZDTX,
     &          NC, INDXC, ZSILL, NLC, BConn1, BConn2, WIDTH, DPEFF )
      

C ---------- External arguments TRHINIT:

C ------ Out:
      real*8    DEPTH(*)
C              - Depth values defining layer division
      INTEGER NBI
C              - dimensional and actual number of internal basins
      INTEGER NLI,NLE,NLC
C              - number of layers in internal and external basins + conn
      real*8    QINPUT(NBI)
C              - net water influx to basin from land and atmosphere
C                  (unit m3/s)
      real*8    XMIX (NBI)
C              - number of wind-mixed layers, with fraction
      INTEGER INDXI(NBI+1)
C              - layer index limits for internal basins.
      real*8    DENSI(NLI)
C              - Density (sigma-t) of internal layers
      INTEGER NBE
C              - actual number of external basins
      INTEGER INDXE(NBE+1)
C              - layer index limits for external basins.
      real*8    DENSE(NLE)
C              - Density (sigma-t) of layers in external basins.
      real*8    ZSURFE(NBI), DZDTX(NBI)
C              - Surface level of each external basin (high level >0)
C                with time derivative (m/day)
      INTEGER NC
C              - actual number of connections.
      INTEGER INDXC(NC+1)
C              - width array index limits for connections.

      real*8  ZSILL (NC)    ! Sill depths

      INTEGER BConn1(NC), BConn2(NC)
C              - Connected basin numbers.
C                +n: internal basin
C                -n: external basin
      real*8    WIDTH(NLC)
C              - Transport widths of connected layers
C                Values from I1=INDXC(IC)+1 through I2=INDXC(IC+1)
C                specifies width for layer 1 to I2-I1+1
C                in connection number IC.

      real*4 DPEFF (NC)



!      SAVE NBI_ALLOC, NBE_ALLOC, MLC_ALLOC, NC_ALLOC
!      SAVE DDZERO, TRANSP_WIDTH, HLAYER,
!     &     ILC_V, NL_V, IBA_V, ILA_V, IBB_V, ILB_V,
!     &     VTOTD_LAND_SURFACE, ZSURF_EXT, DZDT_EXT



C ---------------------- local variables ------------------------


      INTEGER IB, IC, ILC, NL, IBA, IBB, ILA, ILB

      real*8 SEC_PER_DAY
      parameter (SEC_PER_DAY = 24.*3600.)


C ------------------- EXECUTE TRHITR_PREP --------------------------

C -------- If necessary, (re)allocate storage arrays:
$if defined  ALLOCATE_ARRAYS
      IF (      INDXC(NC+1).NE.MLC_ALLOC
     &     .OR. NBI.NE.NBI_ALLOC
     &     .OR. NBE.GT.NBE_ALLOC
     &     .OR. NC.NE.NC_ALLOC      ) THEN
         DEALLOCATE ( DDZERO, TRANSP_WIDTH,
     &        HLAYER, TRANSP_CONSTANT,
     &        VTOTD_LAND_SURFACE,
     &        ZSURF_EXT, DZDT_EXT,
     &        WIND_CURRENT, I_MAX_WINDC,
     &        ILC_V, NL_V, IBA_V,
     &        ILA_V, IBB_V, ILB_V,
     &     STAT = Alloc_chk )
         MLC_ALLOC  = INDXC(NC+1)
         NBI_ALLOC  = NBI
         NBE_ALLOC  = NBE
         NC_ALLOC   = NC
         ALLOCATE ( DDZERO(MLC_ALLOC), TRANSP_WIDTH(MLC_ALLOC),
     &        HLAYER(MLC_ALLOC), TRANSP_CONSTANT (NC_ALLOC),
     &        VTOTD_LAND_SURFACE (NBI_ALLOC),
     &        ZSURF_EXT (NBE_ALLOC), DZDT_EXT(NBE_ALLOC),
     &        WIND_CURRENT(NC_ALLOC), I_MAX_WINDC(NC_ALLOC),
     &        ILC_V(NC_ALLOC),  NL_V(NC_ALLOC), IBA_V(NC_ALLOC),
     &        ILA_V(NC_ALLOC), IBB_V(NC_ALLOC), ILB_V(NC_ALLOC),
     &     STAT = Alloc_chk)

         If ( Alloc_chk .ne. 0 ) then
            WRITE(*,*)
     &       'Error in allocation of work arrays in sub routine TRHITR_PREP: ',
     &            Alloc_chk
            STOP
         ENDIF
      ENDIF
$else

         MLC_ALLOC  = INDXC(NC+1)
         NBI_ALLOC  = NBI
         NBE_ALLOC  = NBE
         NC_ALLOC   = NC
$endif


C ------------- Store indexes (not really necessary for each step)
      DO IC = 1,NC_ALLOC
         ILC_V(IC) = INDXC(IC)+1
          NL_V(IC) = INDXC(IC+1)-INDXC(IC)
         IBA = BCONN1(IC)
         IBA_V(IC) = IBA
         ILA_V(IC) = INDXI(IBA)+1
         IBB = BCONN2(IC)
         IBB_V(IC) = IBB
         IF(IBB.GT.0) THEN
C   ........... Between internal basins:
            ILB_V(IC) = INDXI(IBB)+1
         ELSE
C   ........... Between internal (A) and external (B) basin:
            IBB = ABS(IBB)
            ILB_V(IC) = INDXE(IBB)+1
         ENDIF
      END DO


C ------ store volume change rates due to net land/surface input:
      DO IB = 1,NBI_ALLOC
         VTOTD_LAND_SURFACE(IB) =  QINPUT(IB) ! m3/s
      ENDDO

C ------ store fixed external levels with time derivatives,
C        and update mean value for error tolerance scale:
      DZDT_NORM  = 0.0
      DO IB = 1,NBE_ALLOC
         ZSURF_EXT(IB) = ZSURFE(IB)
         DZDT_EXT(IB)  = DZDTX(IB)/SEC_PER_DAY ! m/s ( from m/day)
         DZDT_NORM  = DZDT_NORM + ABS(DZDT_EXT(IB))
      ENDDO
      DZDT_SUM  = 0.979*DZDT_SUM + DZDT_NORM
      DZDT_WGT  = 0.979*DZDT_WGT + NBE_ALLOC
      DZDT_NORM = DZDT_SUM/DZDT_WGT
C           accessible to ITERATE_SURFACES by COMMON block

C -------------- Initiate:
      DO IC = 1,NC_ALLOC


C  ......... Set constant for calculating effective pressure
C            from density difference:
         TRANSP_CONSTANT(IC) =  2.0*MAX(0.0001,DPEFF(IC))*9.81/1000
                                ! 2*DPEFF*grav/rho(sigma-t)
C            DPEFF = fraction of energy converted into kinetic energy.
C                    With no friction losses: DPEFF =1.


C  ......... Prepare the controls for horisontal flow which are
C            independent of surface levels: wind and depth profile.

         ILC  = ILC_V(IC)
         NL   = NL_V(IC)
         IBA  = IBA_V(IC)
         ILA  = ILA_V(IC)
         IBB  = IBB_V(IC)
         IF(IBB.GT.0) THEN
C   ........... Between internal basins:
            ILB = ILB_V(IC)
            CALL TRHITR_CONN ( ZSILL(IC), NL, DEPTH,
     &             DENSI(ILA), DENSI(ILB), WIDTH(ILC),
     &             DDZERO(ILC), TRANSP_WIDTH(ILC), HLAYER(ILC),
     &             IC, XMIX(IBA), XMIX(IBB),
     &             WIND_CURRENT(IC), I_MAX_WINDC(IC) )
         ELSE
C   ........... Between internal (A) and external (B) basin:
            IBB = ABS(IBB)
            ILB = ILB_V(IC)
            CALL TRHITR_CONN ( ZSILL(IC), NL, DEPTH,
     &             DENSI(ILA), DENSE(ILB), WIDTH(ILC),
     &             DDZERO(ILC), TRANSP_WIDTH(ILC), HLAYER(ILC),
     &             IC, XMIX(IBA), XMIX(IBA),
     &             WIND_CURRENT(IC), I_MAX_WINDC(IC) )
         ENDIF   ! (Use internal basin mixing status for directions)
      END DO
      END SUBROUTINE
      
      

      SUBROUTINE TRHITR_INIT
      DZDT_SUM =0.0   ! Reinitiated here to get identical
      DZDT_WGT =0.0   ! successive runs, called by model program
      END SUBROUTINE


C ===========================================================
      SUBROUTINE TRHFNC ( ZSURFI, VTOTD, DVD_DZ, DVD_DTX, NBI,
     &               NLC, UFLOW )
C Note : Arguments returned in double precision.
C ------------ External arguments to entry-point TRHFNC below:
C  In:
      real*8  ZSURFI (NBI)  ! Internal surface levels
      INTEGER NBI               ! Number of internal basins
      INTEGER NLC               ! number of layers in connections
C Out:
      real*8  VTOTD  (NBI)      ! Rate of volume change, internal basins
      real*8  DVD_DZ (NBI,NBI)  ! Effect of perturbation of ZSURFI.
      real*8  DVD_DTX (NBI)     ! Time derivative of VTOTD due to DZDTX
      real*8  UFLOW(NLC)


C ----- Dummy values for external basins:
      real*8  VTOTD_EXT, DVD_DZ_EXT
      real*8  DVDIBA_DZEXT


C        Array DVD_DZ (IBA,IBB) Accumulates the Jacobian matrix for
C        VTOTD(IBA) as function of ZSURFI(IBB),
C        used in calling subroutine for iteration.

      INTEGER IA, IB, IC, ILC, NL, IBA, IBB, ILA, ILB


C ------------ For given surface levels, calculate transports:

C  1. ---  initiate volume change rates with net land/surface input:
      DO IB = 1,NBI
         VTOTD(IB) = VTOTD_LAND_SURFACE(IB)
         DO IA = 1,NBI
            DVD_DZ (IA,IB) = 0.0
         END DO
         DVD_DTX (IB) = 0.0
      END DO
      VTOTD_EXT = 0.0    ! Accumulated external basin values not used,
      DVD_DZ_EXT = 0.0   ! accumulators zeroed to avoid overflow.
C  ( not really necessary, since the additions becomes insignificant
C    if accumulator is 10**11 times greater than addition )

C  2. ---  calculate transports for each connection:
      DO IC = 1,NC_ALLOC
         ILC  = ILC_V(IC)
         NL   =  NL_V(IC)
         IBA  = IBA_V(IC)
         ILA  = ILA_V(IC)
         IBB  = IBB_V(IC)
         IF(IBB.GT.0) THEN
C   ........... Between internal basins:
            ILB = ILB_V(IC)
            CALL TRHFNC_CONN( NL, DDZERO(ILC), TRANSP_WIDTH(ILC),
     &         HLAYER(ILC), ZSURFI(IBA), ZSURFI(IBB), 
     &         TRANSP_CONSTANT(IC),
     &         WIND_CURRENT(IC), I_MAX_WINDC(IC),
     &         VTOTD(IBA), VTOTD(IBB),
     &         DVD_DZ(IBA,IBA), DVD_DZ(IBA,IBB),
     &         DVD_DZ(IBB,IBA), DVD_DZ(IBB,IBB),
     &         UFLOW(ILC))
         ELSE
C   ........... Between internal (A) and external (B) basin:
            IBB = ABS(IBB)
            ILB = ILB_V(IC)
            DVDIBA_DZEXT = 0.0
            CALL TRHFNC_CONN ( NL, DDZERO(ILC), TRANSP_WIDTH(ILC),
     &         HLAYER(ILC), ZSURFI(IBA), ZSURF_EXT(IBB),
     &         TRANSP_CONSTANT(IC),
     &         WIND_CURRENT(IC), I_MAX_WINDC(IC),
     &         VTOTD(IBA), VTOTD_EXT,
     &         DVD_DZ(IBA,IBA), DVDIBA_DZEXT,
     &         DVD_DZ_EXT,      DVD_DZ_EXT,
     &         UFLOW(ILC))
C                    DVD_DZ for external basins are not used,
C                    and accumulates in dummy variable.

C             Accumulate time derivative of VTOTD due to DZDT_EXT:
                 DVD_DTX(IBA) = DVD_DTX(IBA)
     &                           + DVDIBA_DZEXT*DZDT_EXT(IBB)
C                                          UNIT (m3/s)/s
         ENDIF
      END DO

      END SUBROUTINE


C ======================================================================
      SUBROUTINE TRHITR_CONN( ZSILL, NL_C, DEPTH,
     &                        DENSA, DENSB, WIDTH,
     &                        DDZERO, TRANSP_WIDTH, HLAYER,
     &                        IC, XMIX_A, XMIX_B,
     &                        WIND_CURRENT, I_MAX_WINDC )
      

C Called first in each transport iteration:
C Calculates constant values for transport calculation across connection

C Input:
      real*8 ZSILL     ! sill depth
      INTEGER NL_C
C           Number of connected layers
      real*8    DEPTH(NL_C+1)
C           Limiting depth of layers
      real*8    DENSA(NL_C), DENSB(NL_C)
C           Density (sigma-t) of connected layers.
      real*8    WIDTH(NL_C)
C              - Mean transport widths of connected layers
C Output:
C     (Double precision to get accurate estimates in TRANH_ITR_VTOTD)
      real*8  DDZERO (NL_C)
C           Depth integrated density differences
      real*8  TRANSP_WIDTH (NL_C), HLAYER (NL_C)
C           Width and depth of each layer.
      integer IC
      real*8  XMIX_A, XMIX_B
      real*8  WIND_CURRENT, I_MAX_WINDC


C --------------------------------------
      INTEGER IL
      real*8 DDINTG, DDINCR


C ............ Get surface current and depth extension of
C              wind-induced current:
      CALL GET_WIND_CURRENT( nf, DEBUG, IC, XMIX_A, XMIX_B,
     &                   WIND_CURRENT, I_MAX_WINDC )


      DDINTG = 0.0
      DO IL =  1,NL_C

         HLAYER(IL) = ( MIN(ZSILL,DEPTH(IL+1)) - DEPTH(IL) )
         TRANSP_WIDTH(IL) = WIDTH (IL)

         DDZERO (IL) = DDINTG   ! Store value at top of layer
         DDINCR  = (DENSA(IL)-DENSB(IL)) * HLAYER(IL)
         DDINTG  = DDINTG + DDINCR

      END DO
      END SUBROUTINE

C =======================================

      SUBROUTINE TRHFNC_CONN ( NL_C, DDZERO, TRANSP_WIDTH,
     &       HLAYER, ZSURFA, ZSURFB, TRANSP_CONSTANT,
     &       WIND_CURRENT, I_MAX_WINDC,
     &       VDA, VDB, DVDA_DZA, DVDA_DZB, DVDB_DZA, DVDB_DZB,
     &       UFLOW)
      
C Called from function evaluator during iteration:
C Results and other critical intermediate values in double precision


C Input:
      INTEGER NL_C                ! Number of connected layers
      real*8  DDZERO (NL_C)       ! Depth integrated density diff.:
      real*8  TRANSP_WIDTH (NL_C) ! at top of each layer
      real*8  HLAYER (NL_C)       ! Layer thickness of each layer
      real*8  ZSURFA, ZSURFB      ! Surface levels of basins A and B
      real*8  TRANSP_CONSTANT     ! Pressure efficiency; fraction of
                                  ! energy converted into kinetic energy
      real*8  WIND_CURRENT
      real*8  I_MAX_WINDC

C Output:
      real*8  VDA, VDB            ! Accumulates volume change (m3/s)
      real*8  DVDA_DZA, DVDA_DZB, DVDB_DZA,DVDB_DZB
C           Accumulates derivative of DVx on Zy, x,y = A,B.
      real*8 UFLOW(NL_C)

C --------------------------------------
      INTEGER IL
      real*8 DDSURF, DDINTG
      real*8 U, UU
      real*8 DU_DZ, X
      real*8 TRANSP, DTRANSP, Q_BOTTOM, Q_TOP, DQ_BOTTOM, DQ_TOP, TR

C     -------- set pressure difference at surface:
      DDSURF  = 1000.0 * (ZSURFA-ZSURFB)

$if defined TEST_TRHFNC
         if( debug ) then
           write(nf,'(/''    ------------- TRHFNC_CONN:'')')
           write(nf,'('' DDSURF :'',E17.10)' ) DDSURF
         endif
$endif

C     -------- accumulators for transport from A to B
C              and for the derivative due to changes in ZSURFA:
      TRANSP  = 0.0
      DTRANSP = 0.0

C    --------- No flux at bottom of lower layer.
       Q_BOTTOM = 0.0
      DQ_BOTTOM = 0.0

C    --------- Scans layers from bottom up:
      DO IL = NL_C, 1, -1
C        ........ Value at top of each layer:
C             u(z) = e(e*2*alpha*P/rho)**0.5;
C                     e= +1 for DP>0, e=-1 for DP<0.
C                       is directed velocity, given by Stigebrandt ().
C                     P = DDINTG*g is pressure difference.
         DDINTG = DDZERO(IL) + DDSURF



! --------------------------------------------------------------
!  Directed flow velocity U at top of each layer:

         X = TRANSP_CONSTANT * DDINTG
             ! X = Squared pressure induced velocity with sign
             ! TRANSP_CONSTANT = 2*alpha*grav/rho(sigma-t), see above.
             ! --> X has unit m/s2/(kg/m3) * kg/m2 = m2/s2

c       ........ ADD CONTRIBUTION FROM WIND:
         IF (IL .LT. I_MAX_WINDC) THEN
             U = WIND_CURRENT*( 1.0 - FLOAT(IL) / I_MAX_WINDC )
                               ! Wind induced velocity (m/s)
             X = X + U*ABS(U)  ! Combined square of velocity, with sign
                               ! = U *abs(U), with U=velocity

$if defined TEST_TRHFNC
             IF (DEBUG) then
               WRITE ( nf , '(A, I6, 2(A,G12.5))' )
     &          '  ****** TRANH, IL=',IL, ' X=', X, ' U=',U
             endif
$endif

         Endif

         U = SIGN( DSQRT(ABS(X)), X ) ! signed velocity



! --------------------------------------------------------------
!        Derivative of velocity U
!        with Z = ZURFA-ZSURFB=surface level difference :
!               dX/dZ = dU/dZ*abs(U) + U*d(abs(U))/dZ
!                     =      "       + U*sign(U)*dU/dZ
!                     = 2*abs(U)*dU/dZ
!        ----->
!               dU/dZ = 0.5*(dX/dZ)/abs(U)
!                     = 0.5*1000*TRANSP_CONSTANT/abs(U)
!                     = dpeff*9.81/abs(U):

         UU = MAX( abs(U), 1.d-10 )       ! (m/s)
         DU_DZ = 500.0*TRANSP_CONSTANT/UU ! (s-1)


C        ....... accumulate transport & derivative in each layer:
         Q_TOP   =  U    * TRANSP_WIDTH(IL)

               ! Export value to rest of model:
         UFLOW(IL) = U

         DQ_TOP  = DU_DZ * TRANSP_WIDTH(IL)
         TR = ( Q_TOP +  Q_BOTTOM)*HLAYER(IL)/2.0
         TRANSP  = TRANSP  + TR
         DTRANSP = DTRANSP + (DQ_TOP + DQ_BOTTOM)*HLAYER(IL)/2.0

$if defined TEST_TRHFNC
         if( debug ) then
           write(nf,'('' IL='',i3, 5( 2(1X,A,''='',E17.10)/ ) )' )
     &          IL, 'DDZERO', DDZERO(IL), 'DDINTG',DDINTG,
     &              'U',U, 'Q_TOP', Q_TOP, 'Q_bottom',Q_bottom,
     &              'HLAYER', HLAYER(IL), 'TR',TR,
     &              'DU_DZ',DU_DZ

         endif
$endif
         Q_BOTTOM  = Q_TOP  ! for next layer
         DQ_BOTTOM = DQ_TOP
      END DO


C -------- Update transport values + derivative coeffients (Jacobian)
      VDA = VDA - TRANSP                 ! m3/s
      VDB = VDB + TRANSP

$if defined TEST_TRHFNC
         if( debug ) then
             write(nf,'(2(3x,A,'': '',E17.11))')
     &          'TRANSP', TRANSP, 'DTRANSP', DTRANSP,
     &          '  VDA', vda,    '  VDB',  VDB
         endif
$endif

      DVDA_DZA = DVDA_DZA - DTRANSP      ! m3/m/s
      DVDA_DZB = DVDA_DZB + DTRANSP
      DVDB_DZA = DVDB_DZA + DTRANSP
      DVDB_DZB = DVDB_DZB - DTRANSP
      END SUBROUTINE


C ===================================================================
$if defined TEST_ITR || DEFINED TEST_GAUSS || DEFINED TEST_TRHFNC

      subroutine print_unit
      
$if TEST_ITR >= 3
      nf = DEBUG_UNIT
$else
      integer i
      LOGICAL OPENSTATUS
      if (nf .lt. 0 ) then  ! Search for unused unit number:
         NF = 0
         do i=901,999
            INQUIRE( UNIT = I, OPENED = OPENSTATUS  )
            IF (.NOT. OPENSTATUS) THEN
                NF = I   ! Unused unit number
                EXIT
            ENDIF
         END DO
      ENDIF
$endif
      END SUBROUTINE
$endif

      END Module