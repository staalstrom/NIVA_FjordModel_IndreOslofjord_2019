C Eutrophication model  - File  RANDFNCT.FOR, 
C creating random variations in time series for model EUTRO.

      Module fx_RandFnct
      use fx_RunControl
      use fx_Rand_Vec
      implicit none


$undefine DEBUG_RANDOM
         ! 1,2 or 3 for various degrees of detail in debugging
      

C It is based on an algorithm according to Marsaglia and Zaman 
C stored on C:\ftntools\random\random.for, modified for 
C parallell update of several sequences. 

C The randomization is started by a value RANDOM_SEED_VALUE
C which can be specified, read from a file and kept constant,       
C or created from the computer clock. 

C For Marsaglia and Zaman algorithm, this is converted into
C two values IJ1 and IJ2, conforming to requirements of the 
C subroutine RANDOM, 
C       ( 0<=IJ1<=31328), 
C         0<=IJ2<=30081).

C The RANDFNCT subroutine can handle a series of arrays with
C independent random sequences from this starting point.
C The initial values of IJ1,IJ2 are used for the first sequence, 
C i.e. first element of first array. Starting value of IJ1 is increased 
C by 1 between arrays and wrapped within limits, while IJ2 is increased
C for subsequent elements within an array. 

C The initiation is done for the different arrays during the first call 
C where time is > initial time for that array. In order for some program   
C to reproduce a run exactly, the call sequence at this first time must
C be reproduced. Otherwise, variable call sequences can be used.

      
      INTEGER*4 RANDOM_SEED_VALUE
      LOGICAL   RANDFNCT_INIT
      INTEGER   IJ1, IJ2
      
      contains

C *************************************************************
C Initiate series and reset random seed number, to be used in 
C subsequent calls to RANDOM_SIGNALS:
      
      SUBROUTINE RANDOM_SEED(NSEED)

      
      INTEGER*4 NSEED  ! Variable: Input and/or output, see below
      
      INTEGER ERROR

      character*20 filename /'RandSeed.Dat'/
      
      INTEGER*2 OLD_SEED / -1 /
      SAVE OLD_SEED

      character Cdate*(8), Ctime*(10), Czone*(5) 
      integer TimeArray(1:8)

!      INCLUDE 'RANDFNCT.INC'


C   -------- GET SEED VALUE FROM CALLING PROGRAM,
C            FROM FILE OR FROM INTERNAL CLOCK:

      IF ( NSEED .lt. 0 ) THEN ! READ OLD VALUE FROM FILE:

         open (2,file=filename,IOSTAT=ERROR)
         if (ERROR.EQ.0) THEN
            READ(2,*,IOSTAT=ERROR) Random_seed_value
            close(2)
            OLD_SEED = RANDOM_SEED_VALUE
         endif

         IF(ERROR.NE.0) THEN
            Random_seed_VALUE = 5555555 ! Default value first time
         ENDIF

      elseif ( NSEED .EQ. 0 ) then ! Use time through F95 intrinsic function Date_and_time

         call date_and_time(Cdate, Ctime, Czone,  TimeArray)
         random_seed_VALUE = TimeArray(6)*TimeArray(7)*TimeArray(8) ! minutes*seconds*1000th-part of second
         RANDOM_SEED_VALUE = TimeArray(8) * (RANDOM_SEED_VALUE+1)   ! further modified

      ELSE  ! Use externally specified value

         Random_seed_VALUE = NSeed

      ENDIF


! ****** Obsolete, modify later:

C -------- EXPORT SEED VALUE TO BE USED ALSO BY ACSL RANDOM GENERATOR:

		  ! MS:   CALL SEED ( MS_FORTRAN_SEED(2) ) ! USES LEAST SIGN. PART
      NSEED = RANDOM_SEED_VALUE

C -------- STORE FOR REPEATED RUNS LATER.

      WRITE(*,*)
      WRITE(*,*) 'Start value for random numbers:',Random_seed_value
      if ( random_seed_value .ne. Old_seed ) then
         open (2,file=filename)
         write(2,*) Random_seed_VALUE
         WRITE(2,*) ' SEED FOR NIVA Fjord model'
         close(2)
         WRITE(*,*) 'saved on file ',filename,' for later use'
         WRITE(*,*)
         OLD_SEED = RANDOM_SEED_VALUE
      endif

C ......... first round of calls to RANDOM_SIGNAL will initiate series:      
      RANDFNCT_INIT = .true.

C        Produce Random Seed used for selecting seeds for data series: 
      IJ1 = MOD( ABS(RANDOM_SEED_VALUE), 31329 )
      IJ2 = MOD( ABS(RANDOM_SEED_VALUE/31329), 30082 )

$if defined DEBUG_RANDOM
      write(Debug_Unit,*) ' Called RANDOM_SEED in RANDFNT.FOR, seed=',
     &     RANDOM_SEED_VALUE, '  IJ1, IJ2=', IJ1, IJ2 
$endif
      
      END SUBROUTINE

C ***************************************************************
C Initiate/Update random data for new series:
      SUBROUTINE RANDOM_SIGNALS (T_NOW, T_LAST, NSeries, RESP_FREQ,
     &   IJ, UC, R_TRANSFORM, RANDOM_VALUE )
      
      
      real*8 T_NOW
      INTEGER NSeries
      real*8 T_LAST(NSeries), RESP_FREQ(NSeries)
      INTEGER IJ (2,NSeries )
      real*8 UC(98,NSeries )
      
      EXTERNAL R_TRANSFORM  ! function transforming perturbations

      real*8 RANDOM_VALUE (NSeries)
      
C -----------------------------------------------------      
      real*8 T_INIT, TSTEP/0.02/
      SAVE T_INIT

C --------------- local variables ---------------------
      INTEGER I
      real*8 RV(3), R

      real*8 ALPHA,BETA
      
!      INCLUDE 'RANDFNCT.INC'


$if DEBUG_RANDOM >=2
      integer II 
$endif

$if DEBUG_RANDOM >=1
      write(Debug_Unit,*) ' Called RANDOM_SIGNALS for array at loc. ',
     &     LOC(RANDOM_VALUE), ' T_NOW=',T_NOW
$endif

C     ...... initiate random functions at start of each run:
    
      IF ( RANDFNCT_INIT ) THEN ! First call, initiate time variables:
         T_INIT = T_NOW
         RANDFNCT_INIT = .FALSE.
$if DEBUG_RANDOM >= 2
         write(Debug_Unit,*) ' -- First call in first round'
$endif
      ENDIF
      
      
      IF ( T_INIT .eq. T_NOW ) THEN

$if DEBUG_RANDOM >= 2
         write(Debug_Unit,*) ' -- Call for initial TIME'
$endif
         DO I= 1, NSeries
           RANDOM_VALUE (I) = 0.0
           T_Last(I) = T_NOW
         ENDDO
         IJ(1,1) = -1 ! signals uninitiated
      

      ELSE IF ( IJ(1,1) .eq. -1 ) THEN
$if DEBUG_RANDOM >= 2
         write(Debug_Unit,*) ' -- First call for T > initial TIME',
     &     ' IJ1, IJ2:', IJ1, IJ2
$endif
         DO I= 1, NSeries
           IJ (1,I) = IJ1
           IJ (2,I) = IJ2
           IJ2 = MOD( IJ2+1, 30082)
         ENDDO
         IJ1 = MOD( IJ1+1, 31329 )
         CALL RMARIN ( IJ, UC, NSeries ) 

$if DEBUG_RANDOM >= 2
         DO I= 1, NSeries
            write(Debug_Unit,'(A,I5,A,2I5,A)' ) 
     &        ' Series ', I, ' IJ=', (IJ(II,I),II=1,2), ' UC:'
            write(Debug_Unit,'(5G15.7)')(UC(II,I),II=1,98) 
         ENDDO
$endif
      
      ENDIF


$if DEBUG_RANDOM >= 2
         write ( Debug_Unit,* ) ' T_INIT, T_NOW:',
     &        T_INIT, T_NOW
$endif


C     ...  Use fixed steps in updating random functions, independently
C          of integration steps, which can vary between simulations.
C          Random functions will only depend on seed in start of
C          simulation, ensuring that random functions can be repeated.
C          ORNSTEIN-UHLENBECK function
      
      
      DO I = 1, NSeries

$if DEBUG_RANDOM >= 2
         write(Debug_Unit,'(1x,A,I5,2(A,G12.6))') 
     &    ' Signal nr. ',I,' Old value=', 
     &      RANDOM_VALUE(I),'  Response frequency:', RESP_FREQ(I)
$endif
         
         ALPHA = EXP(-TSTEP*RESP_FREQ(I))
         BETA  = SQRT(1.0-ALPHA**2)

$if DEBUG_RANDOM >= 3
         write(Debug_Unit,'(2(A,G12.6))') 
     &    ' ALPHA: ',ALPHA,' BETA:',BETA 
$endif
         
         DO WHILE (T_NOW-T_LAST(I).gt.1.0e-6)
            
            T_LAST(I) = T_LAST(I)+ TSTEP

C           .. Create 3 random numbers of series:       
C              [0..1>, average 0.5, variance = (0.5**2)/3      
            call ranmar ( IJ(1,I), UC(1,I), 1, RV, 3, 3 )
               ! in \ftntools\random\rand_Vec.for

C           .. Random variable with peaked distribution,  
C              average 0.0, variance = (2*0.5)**2 = 1.0
            R = 2.0*( RV(1)+RV(2)+RV(3) )-3.0

C           .. Transform perturbation (changing R) in externally provided 
C              function R_TRANSFORM:
            CALL R_TRANSFORM (I, R)

C           .. Update series:            
            RANDOM_VALUE(I) = RANDOM_VALUE(I)*ALPHA + R*BETA

$if DEBUG_RANDOM >= 3
         write(DEBUG_UNIT, '(1x,3(A,G12.6),2G12.6,A,I5,A,g12.6,A)' )  
     &    'at T=', T_Last(I), ' R=',R, '  (RV=', RV, 
     %    ') v(', I, ')=', RANDOM_VALUE(I) 
$endif
         
         ENDDO
      ENDDO

$if DEBUG_RANDOM >=1
         write(Debug_Unit,'(3(1x,A,I5,A,G12.6))') 
     &    (' v(', I,')=', RANDOM_VALUE(I), I=1, NSeries)
$endif
      
      END SUBROUTINE

C NON-TRANSFORMING FUNCTION TO BE USED AS R_TRANSFORM ABOVE:
      SUBROUTINE ONE_TO_ONE(i,r)
      INTEGER I
      real*8 R
      END SUBROUTINE

      end Module