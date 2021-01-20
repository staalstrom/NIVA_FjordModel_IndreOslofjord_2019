! ===================================================================
      Module fx_JETMIX

      use fx_RunControl, only: DEBUG_UNIT
      use fx_sigmat

      implicit none

$undefine DEBUG_PRINT

$if defined  DEBUG_PRINT
      logical, private :: debug_j
      INTEGER, private :: DBGUNIT

         ! should only be defined together with DEBUG_PRINT:
$undefine DEBUG_TEST_PRINT
$undefine debug_mode
         ! Debug output from JETMIX to file 987
$undefine debug_trajectory
         ! Debug output of profile to file 989

$endif

			! Debugging is done only if external argument DEBUG = .true.
			! or for repeated calculation if process does not converge.


      Contains

! ===================================================================
! Eutrophication model   - File:    JETMIX.FOR
        ! 13. june 2000:
        ! Corrected error in transport calculations, leading to small
        ! or even negative transports for large number of holes.
        
C Subroutines calculating mixing of discharge from dived outlets,
C using jet integration:
C ===================================================================



C ===================================================================


      SUBROUTINE JETMIX( DEBUG, Layers, DEPTH, SAL, TEMP, DENSITY,
     &      INDEX_BASE, OUTLET_LAYER, OUTFLUX,
     6      QSPP, dDensity_dSPP, OUTFLUX_TEMP,
     &      OUTLET_DEPTH, INTAKE_LAYER, INTAKE_FLUX,
     &      JET_DIAM, NUM_OF_HOLES,
     &   MQ_TR, QJET_TRANSP, QJET_TRANSP_INDEX, NQ_TR, DEST_LAYER,
     &   Ambient_Volume_Flux, Neutral_Depth )
 

C --------------------- In:
      LOGICAL DEBUG   !  .TRUE. triggers debug printout if DEBUG_MODE>0:

      INTEGER Layers        ! Number of density Layers

      real*8 DEPTH(Layers+1)  ! Depth limits for Layers
      real*8 SAL    (Layers)  ! Salinity profile
      real*8 TEMP   (Layers)  ! Temperature profile
      real*8 DENSITY(Layers)  ! Density profile (in sigma-t units)

      INTEGER INDEX_BASE    ! Base of global layer indices

      INTEGER OUTLET_LAYER  ! Global index for outlet
      real*8 OUTFLUX          ! Outlet volume flux in m3/s
      real*8 QSPP             ! Density of particles (g/m3)
      real*4 dDensity_dSPP    ! Fraction of particle concentration
                             ! added to density of water with particles

                            ! Due to contents (particles): (sigma-t)
      real*8 OUTFLUX_TEMP     ! Temperature in deg. C of outflux volume

      real*4 OUTLET_DEPTH     ! Outlet depth (m)

      INTEGER INTAKE_LAYER  ! Global index for intake layer
      real*8 INTAKE_FLUX      ! Volume flux in m3/s for recipient water
                            ! intake

      real*4 JET_DIAM         ! Jet diameter (m)
      INTEGER NUM_OF_HOLES  ! Number of holes (jets) for outlet
                            ! (if <=0: will use 1)
                            
      INTEGER MQ_TR         ! Dimension of transport arrays

C --------------------- Out:
      real*8 QJET_TRANSP (MQ_TR)          ! Transport (m3/s)
      INTEGER QJET_TRANSP_INDEX (2,MQ_TR) ! From/to global layer index
      INTEGER NQ_TR                       ! Number of transports
C                   ( accumulated value <=MQ_TR )
      INTEGER DEST_LAYER                  ! Final receiving global layer
      real*8 Ambient_Volume_Flux  ! volume mixed into dived outlets
      real*8 Neutral_Depth        ! neutral depth for dived outlets


C ---------------------- local variables ------------------------

C ------------- integration states:
      integer states
      parameter (states=5)
      real*8  V(STATES), VN(STATES), D(STATES,4)
      real*8  dVNORM(STATES)
      real*8 S  ! position along trajectory

C --------------- Index values for jet state variables:
      Integer VolumeFlux, MomentumZ, BuoyancyFlux, X, Z
      Parameter ( VolumeFlux=1, MomentumZ=2, BuoyancyFlux=3, X=4, Z=5)

C ----------------------- local variables -----------------------
      Logical First_step
      real*8   Jet_Area
      real*8 Theta
      real*8 Velocity
      real*8 MomentumX    ! Constant throughout integration
      real*8 f, Outlet_Sal(1), Outlet_Temp(1), Outlet_Density(1)
      real*8 Amb_dens
      real*8 dRho_dz

      real*8 Pi

      Parameter ( PI = 3.141 592 653 590 )

      LOGICAL NEUTRAL_POINT
      INTEGER L ,LABOVE ,LBELOW ,HALFINDEX , IQ , LIMIT_EXCEEDED ,K
      Integer Intake_Local_Layer
       
      real*8 MIDDLEDEPTH
      real*8 DEPTHLIMIT(2) ,AMB_DEPTH
      real*8 AMB_DENS_I
      real*8 STEP_NORM , JET_FLUX, PREVIOUSFLUX
      INTEGER Num_Holes_Used

      INTEGER MAX_STEPS / 200 /
      INTEGER STEP_COUNT, SUM_STEPS

      real*8 DEPTHTOLERANCE /0.01/  ! m
      real*8 STEP_MAX /10.0/ ,STEP_MIN  /0.001/
      real*8 ACCNORM / 0.001/
C            global error limit 0.1% sufficient to give accurate results
      real*8  ACCURACY(STATES)

$if defined Debug_Print
      debug_j = debug
$endif

      DO K=1,STATES
         ACCURACY(K) = ACCNORM
      ENDDO

    1 CONTINUE

$if defined DEBUG_MODE
      if (DEBUG_j) then
          dbgunit = 987
          WRITE ( DBGUNIT,* )' ---------- performs JETMIX;',
     &        ' INDEX_BASE:'       ,  INDEX_BASE        ,
     &        ' Layers:'           ,  Layers
          WRITE ( DBGUNIT,* )
     &        ' DEPTH(1)..(Layers)  :',  DEPTH(1),Depth(Layers)
          WRITE ( DBGUNIT,* )
     &        ' DENSITY(1)..(Layers):',  DENSITY(1),DENSITY(Layers)
          WRITE ( DBGUNIT,* )
     &        ' OUTLET_LAYER:'     ,  OUTLET_LAYER      ,
     &        ' OUTFLUX:'          ,  OUTFLUX           ,
     &        ' OUTLET_DEPTH:'     ,  OUTLET_DEPTH
          WRITE ( DBGUNIT,* )
     &        ' INTAKE_LAYER:'     ,  INTAKE_LAYER      ,
     &        ' INTAKE_FLUX:'      ,  INTAKE_FLUX
          WRITE ( DBGUNIT,* )
     &        ' JET_DIAM:'         ,  JET_DIAM          ,
     &        ' NUM_OF_HOLES:'     ,  NUM_OF_HOLES
          WRITE ( DBGUNIT,* )
     &        ' MQ_TR:'            ,  MQ_TR             ,
     &        ' NQ_TR:'            ,  NQ_TR

      endif
$endif

C --------------------------------------------------------------------

      Num_Holes_Used = max(1,NUM_OF_HOLES)

             ! Total jet flux, as sum over holes:

      f = MAX(0.0D0, INTAKE_FLUX)
      JET_FLUX = MAX(0.0D0, OUTFLUX) + f
      IF (JET_FLUX .EQ. 0.0) RETURN
      f = f/JET_FLUX  ! discharge fraction of outlet

      Intake_Local_Layer = Min(Layers,Max(1,INTAKE_LAYER-INDEX_BASE))



          ! Temp. and salinity of discharge + mixed-in recipient water:
      Outlet_Temp(1) = OUTFLUX_TEMP*f + TEMP(Intake_Local_Layer)*(1.0-f)
      Outlet_Sal(1)  = SAL(Intake_Local_Layer)*(1.0-f)
          ! Note: OUTFLUX as freshwater.

          ! Density of mixed outlet before particles:
!      write(999,*) 'calls sigmat with Outlet_Sal, Outlet_TEMP from JETMIX'
      CALL SIGMAT (  Outlet_Sal, Outlet_TEMP, 1, Outlet_Density )
          ! Add particle effect (convert from g/m3 to sigma-t = kg/m3)
      Outlet_Density(1) = Outlet_Density(1) + 1.0e-3*dDensity_dSPP*QSPP

      IQ = NQ_TR

C -----------------------------------

$if defined  DEBUG_TEST_PRINT
      WRITE(dbgunit, '(1x,3A20/1x,3E20.12)')
     &     'Outlet_Temp', 'Outlet_Sal', 'Outlet_Density',
     &      Outlet_Temp, Outlet_Sal, Outlet_Density
      WRITE(dbgunit,'('' ************ Accuracy:''/1x,5E12.4)') Accuracy
$endif

          !  Jet specification, each hole:

      Jet_area  = (PI*JET_DIAM*JET_DIAM/4.0)
      Velocity  = JET_FLUX/Num_Holes_Used/ Jet_area
      MomentumX = Jet_area * Velocity * Velocity ! constant along trajectory
      Theta = 0.0

          !  Horizontal direction initially is assumed,
          !  but code for initial state is prepared for other angles.

C --------------------------------------------------------------------
C 1. Calculate initial state for one jet at end of establishment zone:
C --------------------------------------------------------------------
      V(VolumeFlux) = 2.0*Jet_area*Velocity
                    ! Volume flux doubled by establishment of Gauss profile
      V(MomentumZ ) = 0.0
      V(X) = 6.2*JET_DIAM*Cos(THETA)
      V(Z) = OUTLET_DEPTH-6.2*JET_DIAM*sin(THETA)
           ! BuoyancyFlux is initiated below.

C ----- Initialize layer index (convert to local index; 1=surface layer)
      L = MIN(Layers,Max(1,OUTLET_LAYER-INDEX_BASE))
      STEP_NORM = 1.0
      SUM_STEPS = 0
      S = 0.0

      PREVIOUSFLUX = JET_FLUX/Num_Holes_Used

      First_step = .true.
        

C ---------------------------------------------------------------
C 2. Check/update location of depth interval,
C    and establish local density gradient.
C ---------------------------------------------------------------

C ---------------
  200 continue
C ---------------

         ! ... Search for layer with lower depth below outlet depth:
      do while ( L.lt.Layers .and. Depth(L+1) .lt. V(Z) )
         L = L + 1
      end do
      
         ! ... Search for lowest layer with upper depth above outlet depth:
      do while ( L.gt.1 .and. Depth(L) .ge. V(Z) )
         L = L - 1
      end do

         ! Returns with L = number of lowest layer with top above outlet depth.
         ! L = local index for top of bottom layer if
         ! depth is larger than depth range of density profile.

         ! ... Define depth interval consisting of one half of layer,
         !     and set local linear density profile for use in integration.
      MiddleDepth = (Depth(L+1) + Depth(L))/2.0

$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
         WRITE ( DBGUNIT, * )'       -------- Found depth layer'
         WRITE ( DBGUNIT, * )' L=',L,'  MiddleDepth=',MiddleDepth
      ENDIF
$endif

C ---------------
  210 continue
C ---------------

         ! ... Density gradient for half layer containing current depth V(Z):

         ! Upper half:     ------- Depth (L-1) --------
         !     Labove=L-1
         !                        Densi(L-1)
         !
         !               x -------  Depth (L)  -------- Lower half:
         !     Lbelow=L  x                                 Labove=L 
         !               x       Densi(L)               x 
         !                                              x
         !                 -------  Depth(L+1) -------- x
         !                                                 Lbelow=L+1 
         !                        Densi(L+1)
         !
         !                 -------  Depth(L+2) --------
         
      if ( V(Z) .le. MiddleDepth) then ! upper half:
         Labove = max(1,L-1)
         Lbelow = L
         HalfIndex = 0
         DepthTolerance = 0.01
      else                             ! lower half:
         Labove = L
         Lbelow = min(Layers,L+1)
         HalfIndex = 1
         DepthTolerance = -0.01
      endif

         ! ... Layer depth limits (interval slightly expanded to ensure that
         !     step interpolation at depth limits gets beyond actual limit.)
      DepthLimit(1) = MiddleDepth + DepthTolerance
      DepthLimit(2) = Depth(L+Halfindex) - DepthTolerance

         ! ... Density gradient in current half of layer:
      dRho_dz = ( Density(Labove) - Density(Lbelow)       ) /
     &          ( ( Depth(Labove) - Depth(Lbelow+1) )/2.0 )
         ! (For depths above and below density profile, Labove=Lbelow,
         !  so density gradient is zero)

         ! ... Current ambient density calculated from profile in
         !     current half layer:
      Amb_Dens = Density(L) + dRho_Dz * (V(Z) - MiddleDepth)

$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
        WRITE ( DBGUNIT, * )'       --------- Selected half layer:'
        WRITE ( DBGUNIT, '(1x,3A10)' ) 'Labove','Lbelow','HalfIndex'
        WRITE ( DBGUNIT, '(1x,3I10)' )  Labove,  Lbelow,  HalfIndex
        WRITE ( DBGUNIT, '(1x,2a18)'   ) 'DepthTolerance', 'DepthLimits'
        WRITE ( DBGUNIT, '(1x,3E18.10)' )  DepthTolerance,   DepthLimit
        WRITE ( DBGUNIT, '(1x,a18,a)' ) 'dRho_dz     ',
     &           ' Amb_dens from profile'
        WRITE ( DBGUNIT, '(1x,2E18.10)' ) dRho_dz,  Amb_dens
      ENDIF
$endif


      if (First_Step) Then

             ! ... Initialize jet and ambient density state:
         V(BuoyancyFlux) = V(VolumeFlux)/2.0 *
     &            ( Amb_dens - Outlet_Density(1))
                  ! jet outflux(=half of total flux)
                  ! * density difference between jet and water

         if ( V(BuoyancyFlux).eq.0.0 .and. V(MomentumZ).eq.0.0 )
     &       goto 910  ! horisontal jet with no buoyancy - no mixing
                       ! between layers.  

         Amb_Dens_I  = Amb_Dens ! Stored for use in later steps below
         Amb_Depth = V(Z)

             ! ... Initialize derivative norms for error checking:
         dVNorm(VolumeFlux)   = V(VolumeFlux)/10.
         dVNorm(MomentumZ)    = MomentumZ/10.
         dVNorm(BuoyancyFlux) = V(BuoyancyFlux)/10.
         dVNorm(X) = 1.0
         dVNorm(Z) = 1.0


$if defined DEBUG_MODE
      IF (DEBUG_j) THEN
         WRITE ( DBGUNIT, *)'   --------- first step:'
         WRITE ( DBGUNIT, '(1x,2a20)' ) 'Amb_depth', 'V(BuoyancyFlux):'
         WRITE ( DBGUNIT, '(1x,2E20.4)') Amb_depth ,  V(BuoyancyFlux)
      ENDIF
$endif

      else

             ! ... Adjust gradient to correct for deviations
             !     between ambient density as integrated along trajectory
             !     and the correct ambient density calculated from density
             !     profile in the current half layer
             !     Differences arise because RK3 integrates slightly beyond
             !     depth intervals in DEPTH array.
             !     The adjustment will tend to give correct the error
             !     over one half layer thickness in direction of movement
             !     of jet (sign of D(Z,1) below),
             !     preventing accumulation of errors.
         dRho_dz = dRho_dz
     &            + (Amb_Dens - Amb_Dens_I) /
     &               Sign( DepthLimit(1)-DepthLimit(2), D(Z,1))
     
$if defined  DEBUG_MODE
         IF (DEBUG_j) THEN
      WRITE ( DBGUNIT, *)'   ---------- not first step, adjusts dRHo_dz'
      WRITE ( DBGUNIT, '(1x,3a18)' ) 'Amb_dens_I','Amb_depth', 'dRho_dz'
      WRITE ( DBGUNIT, '(1x,3E18.11)') Amb_dens_I , Amb_Depth, dRho_dz
         ENDIF
$endif


      endif

C ---------------------------------------------------------------
C 3. Integrate until Z is close to top of current layer,
C    or until jet turns (Theta changes sign),
C    taking note of depth of neutral density.
C ---------------------------------------------------------------



$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
         WRITE ( DBGUNIT, '('' V:'',5E14.7)' ) V
         WRITE ( DBGUNIT, * )' ---------- calls RK3 '
      ENDIF
$endif

C  ...... integrate through step by 3.order RK method with error control
      CALL RK3(  DepthLimit, dRho_dz,
     &           STEP_NORM, STEP_MIN, STEP_MAX, MAX_STEPS, S,
     &           States, V, Vn, D, dVNorm, ACCURACY, MomentumX,
     &           Limit_Exceeded, Neutral_point, STEP_COUNT )
      SUM_STEPS = SUM_STEPS + STEP_COUNT
      First_step = .false.

C        Integrated ambient density describes effective density profile
C        during integration of buoyancy flux.  Used to correct density
C        gradient above to avoid accumulating deviations between
C        density given by input table and integrated value.
      Amb_Dens_I  = Amb_dens_I + dRho_Dz * (V(Z) - Amb_Depth)
      Amb_Depth = V(Z)

$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
        WRITE ( DBGUNIT, * )' ---------- returns from RK3 '
        WRITE ( DBGUNIT, '('' V:'',5E14.7)' ) V
        WRITE ( DBGUNIT, '('' D:'',5E14.7))' ) (D(k,1),k=1,5)
        WRITE ( DBGUNIT, '(1x,2a20)' ) 'Limit_Exceeded', 'Neutral_point'
        WRITE ( DBGUNIT, '(1x,I20,L20)') Limit_Exceeded, Neutral_point
        WRITE ( DBGUNIT, '(10X,''STEP_COUNT:'',I8)' ) STEP_COUNT
        WRITE ( DBGUNIT, '(10X,''SUM_STEPS:'',I8)' ) SUM_STEPS
      ENDIF
$endif


C  ....... Check return value:
      IF (NEUTRAL_POINT) THEN
          GOTO 900  ! End of calculation, return code below
      ENDIF


C  ....... Prepare for continued integration:

      IF (LIMIT_EXCEEDED.EQ.1) THEN  ! Passed middepth

$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
         WRITE ( DBGUNIT, * ) ' passed middle depth of layer ',L
      ENDIF
$endif

          GOTO 210      ! continue with next half of interval

      ENDIF


      IF (LIMIT_EXCEEDED.EQ.2) THEN  ! Passed to next layer
             ! Checks if endpoint of trajectory:         
         if (     ( DepthLimit(2).le.Depth(1)
     &              .and. V(BuoyancyFlux).ge.0.0 )
     &       .or. ( DepthLimit(2).ge.Depth(Layers+1)
     &              .and. V(BuoyancyFlux).le.0.0 )
     &      ) then
              GOTO 900   ! Stop calculation
          ENDIF

          IQ = IQ+1     ! store accumulated entrainment as transport:
          IF ( IQ .LE. MQ_TR ) THEN
              QJET_TRANSP(IQ) = ( V(VOLUMEFLUX) - PREVIOUSFLUX) *
     &                            Num_Holes_Used
              QJET_TRANSP_INDEX(1,IQ) = L + INDEX_BASE

$if defined  DEBUG_MODE
      IF (DEBUG_j) THEN
         WRITE (DBGUNIT, * )' Entrainment transport from layer  ',L
         WRITE (DBGUNIT,
     &   '('' QJET_TRANSP('',I3,
     &                 '')= (V(VOLUMEFLUX)-PREVIOUSFLUX)*NHOLE''/
     &     1x, E14.7,'' = ('',E14.7,'' - '',E14.7,'')*'',I5)' )
     &     IQ, QJET_TRANSP(IQ), V(VOLUMEFLUX), PREVIOUSFLUX,
     &                                      Num_Holes_Used
      ENDIF
$endif
              PREVIOUSFLUX = V(VOLUMEFLUX)
          ENDIF

          GOTO 200      ! continue with next layer
      ENDIF

      WRITE(*,*) 'JETMIX calculation stopped prematurely (max steps?)'

$if defined  DEBUG_MODE
      IF (DEBUG_J) THEN
         WRITE (*,*) 'Returns from RK3 with LIMIT_EXCEEDED <> 1,2',
     &                ' and NEUTRAL_POINT = .FALSE.'
         WRITE (*,*) ' THIS SHOULD NEVER OCCUR !!'
!         PAUSE
      ELSE
         DEBUG_J=.true.
         goto 1
      ENDIF
$endif

 900  CONTINUE

$if defined  DEBUG_TEST_PRINT
      IF ( ACCURACY(1) .EQ. ACCNORM ) THEN
           DO K=1,states
              ACCURACY(k) = 10.0*ACCURACY(K)
           ENDDO
           write(DBGUNIT,*)
     &        '** RECALCULATES WITH 10*ORIGINAL TOLERANCE'
           GOTO 1  
      ENDIF
$endif


           ! End of calculation.
           ! and other information:

  910 Ambient_Volume_Flux = V(VolumeFlux)*Num_Holes_Used - JET_FLUX
      Neutral_Depth       = V(Z)


$if defined  DEBUG_TEST_PRINT
      write(DBGUNIT,'('' IN JETMIX: ''/1x,2A15)')
     &         'Amb_Vol_Flux','Neutral_D'
      write(DBGUNIT,'(1x,2G15.6)')
     &         Ambient_Volume_Flux, Neutral_Depth
$endif

           ! NOTE: last part of volume flux, added within destination layer,
           !       will not show in QJET_TRANSP


           ! store index to destination Layers

      DEST_LAYER = L + INDEX_BASE

      IF (INTAKE_FLUX .GT.0.0 ) then

        IF (Intake_Local_Layer + INDEX_BASE .ne. Dest_layer ) THEN
                        ! store recipient intake into jet 
                        ! as internal transport.
          IQ = IQ+1 
          IF ( IQ .LE. MQ_TR ) THEN
            QJET_TRANSP(IQ)         = INTAKE_FLUX
            QJET_TRANSP_INDEX(1,IQ) = Intake_Local_Layer + INDEX_BASE

$if defined  DEBUG_MODE
            IF (DEBUG_j) THEN
               WRITE (DBGUNIT,
     &           '('' QJET_TRANSP('',I3,'') = INTAKE_FLUX ='', E14.7)' )
     &                             IQ, QJET_TRANSP(IQ)
            ENDIF
$endif

          ENDIF
        ENDIF

$if defined  DEBUG_MODE
        IF (DEBUG_j) THEN
           WRITE ( DBGUNIT, * )' + Intake from local layer nr. ',
     &                             Intake_Local_Layer
        ENDIF
$endif

      ENDIF

      DO K = NQ_TR+1, min(MQ_TR, IQ)
           QJET_TRANSP_INDEX(2,K) = Dest_layer
      ENDDO
      NQ_TR = min(MQ_TR, IQ)

      RETURN


      END Subroutine


C ===================================================================
C Subroutine for 3.order Runge Kutta integration with error control:
C Returns with integrated values
C      - at one of Ambient.Depths,
C      - at neutral depth,

      SUBROUTINE RK3(  DepthLimit, dRho_dz,
     &                 STEP_NORM, STEP_MIN, STEP_MAX, MAX_STEPS, S,
     &                 N, V, Vn, D, dVNorm, ACCURACY, MomentumX,
     &                 Limit_Exceeded, Neutral_point, STEP_COUNT )
      

C Relevant recipient data:
      real*8 DepthLimit(2) ! (I)   limits of current depth interval
      real*8 dRho_dz       ! (I)   Vertical density gradient

C Integration step control
      real*8 STEP_Norm          ! (I/O) Estimate of normal stepsize
      real*8 STEP_MIN, STEP_MAX ! (I)   Absolute limits on allowed step
      INTEGER MAX_STEPS       ! (I)   Maximum number of evaluated steps,
C                                     both accepted and rejected
      real*8 S  ! accumulated position along trajectory

C Jet data with work arrays:
      INTEGER N  ! (I)   Number of variables
      real*8 V(N)  ! (I/O) Accepted, error-free jet state values
      real*8 Vn(N) ! (O)   Work array, contains new estimate of V
      real*8 D(N,3)! (O)   Work array, derivatives
      real*8 dVNorm(N)   ! (I)   derivative Norm for checking error
      real*8 ACCURACY(N) ! (I)   Relative accuracy
      real*8 MomentumX   ! (I)   Momentum in X direction, constant.

C Out: result status
      Integer Limit_Exceeded  ! = index of DepthLimit exceeded
      Logical Neutral_point   ! neutral point reached
      INTEGER STEP_COUNT      ! Number of steps required

C --------------- Index values for jet state variables V and D:
      Integer BuoyancyFlux, Z
      Parameter ( BuoyancyFlux=3, Z=5)

C -------------------- local variables ------------------------
      real*8 COEFF(2) / 0.5, 0.75 /
      real*8 STEP, STEP_FACTOR, ERR3, ERRL1, ERRL2, STEP_DONE, SNEW
      INTEGER I,K, POINT_ITERATION
      

$if defined  DEBUG_PRINT
      if (debug_j) WRITE ( DBGUNIT, * )'      ----- performs RK3 '
$endif

      STEP_COUNT = 0
      SNEW = S

C ------- Starts here for each new integration step:
   10 S = SNEW
      CALL DERIVATIVES (.true., S, dRho_dz, MOMENTUMX, V, D(1,1) )

C   Accumulate maximum derivative scale for 3 first states:
C                  VolumeFlux,  MomentumZ,  BuoyancyFlux:
         do k=1,3
            dVNorm(k) = max(abs(V(k))/10. , dVNorm(k) )
         enddo

$if defined  DEBUG_PRINT
      if (debug_j) then
         WRITE ( DBGUNIT, * ) ' >>>>>>>>>> New integration Step::'
         WRITE ( DBGUNIT, '('' dVNorm:'',5E14.7)' ) dVNorm
         WRITE ( DBGUNIT, '(''  V:'',5E14.7)' ) V
         WRITE ( DBGUNIT, '('' D(I,1):'',5E14.7)' )  (D(I,1),I=1,N)
      endif
$endif

C ------- Restarts here when repeating due to error control:
   20 CONTINUE
      STEP = MAX( STEP_MIN, MIN(STEP_MAX, STEP_NORM) )

      POINT_ITERATION = 0

C ------- Restarts here with step set by point iteration:
   30 IF ( STEP_COUNT .GE. MAX_STEPS ) RETURN
      STEP_COUNT = STEP_COUNT+1

$if defined  DEBUG_PRINT
      if (debug_j) WRITE ( DBGUNIT, '(10X,''Try STEP:'',E14.7)' ) STEP
$endif

C ------- intermediate derivative values for 2. and 3.order integration.
      DO K=1,2
         DO I=1,N
             VN(I) = V(I) + COEFF(K)*STEP*D(I,K)
         ENDDO
         CALL DERIVATIVES 
     &       (.false., -1.0D0, dRho_dz, MOMENTUMX, VN, D(1,K+1) )
C                          negative S: suppresses printout of values

$if defined  DEBUG_PRINT
      if (debug_j) then
         WRITE ( DBGUNIT, '(1X,''VN:'',5E14.7)' ) VN
         WRITE ( DBGUNIT,
     &           '('' D(I,'',I1,'')'',5E14.7)') K+1, (D(I,K+1),I=1,N)
      endif
$endif

      ENDDO

C  ----------- estimate new 2.order state variable values and
C              error by comparing with 3.order integration.
C              Then decide action:
C                - reject due to error and repeat with smaller step:
C                - reset step to defined termination point and repeat
C                - accept and continue

$if defined  DEBUG_PRINT
      if (debug_j) then
         WRITE (DBGUNIT,*) ' ------- ERROR CONTROL:'
         WRITE (DBGUNIT,'(1X,A10,4A14)')
     &       'STATE NR', 'ERR3', 'ERRL1', 'ERRl2', 'STEP_FACTOR'
      endif
$endif

      STEP_FACTOR = 12.5*STEP_NORM/STEP
      SNEW = S + STEP

      DO I=1,N
C        ------- 2.order estimate of new value:
         VN(I) = V(I) + STEP*D(I,2)
C        ------- estimate 3.order error over step:
         ERR3 = STEP* ( 2*D(I,1) - 6*D(I,2) + 4*D(I,3) )/9.0
C        ------- Scale for permitted local error from maximum of:
C             - change of value over step or by norm, linear with step:
         ERRL1 = MAX( ABS(VN(I)-V(I)), ABS(DVNORM(I)*STEP))*ACCURACY(I)
C             - double integrated 2. derivative, 2. order in step:
         ERRL2 = ABS( ( VN(I) - (V(I)+D(I,1)*STEP) ) * ACCURACY(I) )

C        ------- adjust step to as large a possible value,
C                but with all errors within limits:
         if (ERR3 .gt. 0.0 ) THEN
            STEP_FACTOR = MIN ( STEP_FACTOR,
     &                          MAX ( ( ERRL1/ERR3 )**0.5,
     &                                ( ERRL2/ERR3 )       ) )
         Endif

$if defined  DEBUG_PRINT
      if (debug_j) WRITE (DBGUNIT,'(1X,I10,4E14.7)')
     &       I, ERR3, ERRL1, ERRl2, STEP_FACTOR
$endif

      ENDDO

$if defined  DEBUG_PRINT
      if (debug_j)
     &   WRITE ( DBGUNIT, '(10X,'' 2.order estimate VN:''/1x,5E14.7)' ) VN
$endif


C         In preparation for next step:
      STEP_NORM = 0.9*STEP_FACTOR * STEP
C              Factor 0.9 to reduce rejections and avoid
C              infinite loop due to minimum fluctuation in step.

$if defined  DEBUG_PRINT
      if (debug_j)
     &   WRITE ( DBGUNIT, '(10X,''STEP_NORM:'',E14.7)' ) STEP_NORM
$endif

      if (STEP_FACTOR .lt. 1.0 .AND. STEP .GT. STEP_MIN ) THEN
          IF (POINT_ITERATION .EQ.0 ) GOTO 20
      ENDIF

$if defined  DEBUG_PRINT
      if (debug_j) WRITE ( DBGUNIT, '(10X,'' step was accepted'')')
$endif

C ------------ HAS PERFORMED STEP WITH ACCEPTABLE ERROR --------------
      STEP_DONE = STEP


         ! --- Check if some of specified termination points are inside STEP:
         !     If found, use lowest non-zero positive step ending in such
         !     point for renewed calculation, with maximum 3 iterations.
         !     Variable STEP is reset in the process to successively
         !     smaller values if limitations are found: 

         ! --- Check against limits of depth interval, slightly expanded:
      Limit_Exceeded = 0
      DO K=1,2
        IF (ITP_STEP( Vn(Z), V(Z), D(Z,1), DepthLimit(K),
     &                 STEP, STEP_DONE ) ) THEN
           Limit_Exceeded = K
        endif
      enddo
         
         ! --- Neutral depth (Buoyancy flux = 0)
      Neutral_point = ITP_STEP( Vn(BuoyancyFlux), V(BuoyancyFlux),
     &                          D(BuoyancyFlux,1), 0.0D0,
     &                          STEP, STEP_DONE )

      if (Neutral_point) Limit_Exceeded=0

$if defined  DEBUG_PRINT
      if (debug_j) then
        if (Limit_Exceeded .gt. 0)
     &     WRITE(DBGUNIT, '(10x,'' Depth limit '',I4,'' exceeded'')')
     &                                    Limit_Exceeded
        if (Neutral_point)
     &   WRITE ( DBGUNIT, '(10x,'' Passed neutral point'')')
      endif
$endif


C   ------ repeat up to 3 times to reach point more accurately:
      if ( STEP .lt. STEP_DONE-STEP_MIN .and.
     &     POINT_ITERATION.le.2 ) THEN
          POINT_ITERATION = POINT_ITERATION + 1

$if defined  DEBUG_PRINT
      if (debug_j)
     &     WRITE ( DBGUNIT, 
     &             '(10x,'' Reiterate with interpolated step'')')
$endif

          goto 30
      endif

C ------------ step accepted, store new values:
      DO I=1,N
          V(I) = Vn(I)
      ENDDO


C ----------- continue integration until stop condition is encountered:
      if (Limit_Exceeded.eq.0 .and. .not. Neutral_point ) goto 10

C Final call to update derivatives and print final values
C in debug output:
$if defined  DEBUG_PRINT
      if (debug_j)
     &      CALL DERIVATIVES (.true., S, dRho_dz, MOMENTUMX, V, D(1,1) )
$endif

      END Subroutine RK3


C ========= Find reduced step for interpolating to specified value:
      LOGICAL FUNCTION ITP_STEP ( VNew, V, B, TARGET, STEP, STEP_DONE )
      
      real*8 VNew, V, B, TARGET, STEP, STEP_DONE

C ------------ local variables ----------
      real*8 A, C, ROOTARGUMENT, ROOT, STEP_LIMIT
      INTEGER ROOTSIGN

      ITP_STEP = .false.
 
C -----  2. order approximation over STEP_DONE
C             V^(dS) = V + B*dS + (Vn-V-B*STEP_DONE)*(dS/STEP_DONE)**2
      A = (VNew-V-B*STEP_DONE)/(STEP_DONE**2)
      C = V-TARGET
      IF ( A .NE. 0.0) THEN
         ROOTARGUMENT = B*B-4*A*C
         IF ( ROOTARGUMENT .GE. 0.0 ) THEN
             ROOT = SQRT(ROOTARGUMENT)
             DO ROOTSIGN = -1,+1,2
                STEP_LIMIT = (-B + ROOTSIGN*ROOT)/(2.0*A)
                IF (STEP_LIMIT.GT.0.0 .AND. STEP_LIMIT.LE.STEP ) THEN
                    ITP_STEP = .TRUE.
                    STEP = STEP_LIMIT
                ENDIF
             ENDDO
         ENDIF
      ENDIF

C ------ 1.order interpolation?
      IF (.NOT. ITP_STEP) THEN
         IF ( VNew.ne.V .and. (VNew-TARGET)*(TARGET-V) .GE. 0.0 ) THEN
             STEP_LIMIT = STEP_DONE*ABS( (TARGET-V)/(VNew-V) )
             IF (STEP_LIMIT.GT.0.0 .AND. STEP_LIMIT.LE.STEP ) THEN
                 ITP_STEP = .TRUE.
                 STEP = STEP_LIMIT
             ENDIF
         ENDIF
      ENDIF

      END Function


C ================================================================
      SUBROUTINE DERIVATIVES( ValidState, S, dRho_dz, MomentumX, V, D )
      

! In:              
      Logical ValidState ! .True. if  valid state
                         ! i.e. not intermediate calculation
      real*8 S           ! Position along trajectory (if <1: intermediate values)
      real*8 dRho_dz     ! Ambient vertical density gradient
      real*8 MomentumX   ! Constant horizontal momentum of jet.
      real*8 V(5)        ! Integrated state variables at position S
! Out:
      real*8 D(5)        ! Derivatives of V in same structure


C --------------- Index values for jet state variables V and D:
      Integer VolumeFlux, MomentumZ, BuoyancyFlux, X, Z
      Parameter ( VolumeFlux=1, MomentumZ=2, BuoyancyFlux=3, X=4, Z=5)


C ---------------------- local variables --------------------------

C -------- Constant values and constant expressions:
      real*8
     .  Gravity, RhoScale,
     .  Alpha,  Lambda,  LambdaSq,  LambdaFormFactor,
     .  PI

      Parameter
     .( Gravity = 9.81 ,  RhoScale = 1000.0,
     .  Alpha = 0.057  ,  Lambda  = 1.18  ,  LambdaSq = Lambda*Lambda,
     .                LambdaFormFactor = (LambdaSq + 1.0) / LambdaSq ,
     .  PI = 3.141 592 653 590
     .)

      real*8 DASIN

C ----------- Jet description variables:
      real*8 Momentum, Velocity
      real*8 BSq, B, RhoDiff
      real*8 Theta, SinTheta, CosTheta
      real*8 dRho_ds
      real*8 Dummy1
      logical Dummy2
C ----------- Calculate variables given by state:
C Jet described by Gauss distribution for velocity and density dev.
C   u(r) = u(0)*exp(-(r/b)**2)
C   RhoDiff(r) = RhoDiff(0)*exp(-(r/(lambda*b))**2)
C which leads to:
C   Momentum     = PI *b**2 * u(0)**2 / 2
C   VolumeFlux   = PI *b**2 * u(0)
C   BuoyancyFlux = PI *b**2 * u(0)* RhoDiff(0) / Lambdaformfactor
C the variable b is related to the nominal width of the jet by
C           Width = sqrt(2)*b

      Momentum = SQRT(V(MomentumZ)*V(MomentumZ) + MomentumX*MomentumX )
     &            + 1.0e-20  ! To handle initial phase if MomentumX = 0
      SinTheta  = - V(MomentumZ) / Momentum     ! >0 for upward jet
      Theta     = DASin (SinTheta)
      CosTheta  = Cos (Theta)
      Velocity  = 2.0 * Momentum / V(VolumeFlux)  ! in center of jet
      BSq   = V(VolumeFlux) / Velocity/PI
      B     = Sqrt(BSq)
      Rhodiff   = - V(BuoyancyFlux) /V(VolumeFlux) *LambdaFormFactor
      dRho_ds   = - dRho_dz *SinTheta

      Dummy1 = S
      Dummy2 = ValidState
      
$if defined DEBUG_MODE && defined DEBUG_TRAJECTORY
      if (debug_j .and. ValidState) then
        if( S.eq.0.0 ) then
           WRITE(989,'(1x,A8,5A14)' )
     &       '   S  ', ' VolumeFlux ', ' MomentumZ ', ' BuoyancyFlux ',
     &       '   X   ', '    Z '
        endif
        if (S.ge.0.0) then
c           WRITE(988,"(1x,6F8.3,3E10.3)" )
c     &        S, V(X), V(Z), Velocity, b, Theta,
c     &        Rhodiff, Momentum, dRho_ds
          WRITE(989,'(1x,F8.3,5E14.7)' ) S, V
        endif

      endif
$endif

C ----------- Calculate derivatives of jet state variables:
      D(VolumeFlux  ) = Alpha * 2.0 * PI * B * Velocity
      D(MomentumZ   ) = PI* Gravity * LambdaSq * BSq * RhoDiff /RhoScale
      D(BuoyancyFlux) = V(VolumeFlux) * dRho_ds
      D(Z)            = - SinTheta
      D(X)            =   CosTheta

      End Subroutine

      End Module