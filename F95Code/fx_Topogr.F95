! ===================================================================
! Eutrophication model   - File:    TOPOGR.95
!                                   Birger Bjerkeng, NIVA.
!
! To be compiled as F95 code, but in fixed F77 format,
! although with 80 character line length (ABSOFT options)


      module fx_Topogr

      use ModelDimensions
      use ModelParam_Topography
      use sub_TextInput

      implicit none

$DEFINE PRINT_DESC
       ! Activates print of topographical description
       ! to file Topo_Dsc.txt

$unDEFINE PRINT_INPUT
        ! Activates echo of input and input processing
        ! results to log file

$unDEFINE PRINT_ADAPT_LOG
        ! Activates diagnostic output of adaption process
        ! to log file

$unDEFINE TRACE_CALLS
        ! activates call trace output to log file
      
      integer, parameter :: DiagFile = 999
      character, private :: ErrorMessage*132
      


C ================== Description of basins: =========================

C Maximum number of internal and external basins:
      INTEGER MBasin, MBasex
         PARAMETER ( MBasin = 30, MBasex=2 )
C Maximum total number of (depth,area) data points i internal basins:
      INTEGER MZB
         PARAMETER ( MZB = 200 )

C Names and status (INTERNAL/EXTERNAL BOUNDARY) of basins:
      CHARACTER*40 BasinName(-MBasex:MBasin)
!      LOGICAL BINTRN(MBas)

C Count and index ranges for internal and external basins:
      INTEGER NBASIN,NBASEX
C        -  Index of topographic data arrays in internal basins:
      INTEGER INDXZB(0:MBasin)
C        -  internal Basin nr. I from INDXZB(I-1)+1 to INDXZB(I).

C Depth and area values entered from file:
      real*8 ZBasin(MZB),ZBasinArea(MZB)
C Calculated values:
      real*8 Vol(MZB),DADZ(MZB)
C        - Vol(k) = calculated volume from bottom.
C        - DADZ(k) = the derivative dArea/dZ between k and k+1.
C               [ DADZ(Kmin) = DADZ(Kmin+1) ]
C   -------------- ZBasin(k),   Area(k),   Vol(k)   -----------
C     layer k       DADZ(k)
C   -------------- ZBasin(k+1), Area(k+1), Vol(k+1) -----------

C ================ Description of connections =====================

C  Max. number of connections:
      INTEGER MCON
C  Max total number of width/depth data pairs:
      INTEGER MZCON
         PARAMETER (MCON=30, MZCON = 300)

C  Actual number of connections entered from file:
      INTEGER NCON
C  Basin numbers connected, two for each connection:
      INTEGER NBCON(2,MCON)
C  Index limits to arrays ZCon, WCon:
      INTEGER NDXCON(0:MCON)
C              (IX = NDXCon(IC-1)+1 ... NDXCon(IC) for connection IC)
C  Depth,width values entered from file:
      real*8 ZCON(MZCON), WCON(0:2,MZCON)
C   WCON(I,K):
C      K=0 transport width,
C      K=1,2: width of buffer zone on each side of transport section.
C
      real*8 BUFFER_LENGTH(2,MCON), TRANSITION_LENGTH(2,MCON)
C           length of buffer zone and mixing zone on each side of
C           transport section.
      real*8 MIXING_TIME(2,MCON)
C           mixing time for each buffer volume
      
		
      contains

               ! log file has name For.nnn,
					! with nnn=Diagnstic_File, opens automatically


!    Topographic description, file TOPOGR.FOR
!    Containing subroutine TOPO,
!    performing topographic setup of model:

!    1. Calls TPRD in module fx_TopoRead:

!           Read topographic description from file.
!           Get number of basins and connections,
!           and bottom and zill depths.

!    2. Set layer division and basin and connection indexes.

!    3. Calls TPBAS:
!           Get volume and area data for layers.

!    4. Calls TPCON:
!           Get basin numbers and transport width of connections.

!    A testprogram exists in TESTTOPO.FOR

! ==================================================================
      SUBROUTINE TOPO( TOPOGRAPHY_FILE,  
     & MD, MBI, MLI_DIM, MLI, MBE, MLE, MC, MLC,
     &                 ND, DEPTH,
     &                 NBI, VTOTZ, ZBOTMI, LSHORE, NLVOPN, VFROPN,
     &                 INDXI, AREA, BOTTOM, VFRAC,
     &                 NBE, ZBOTME, INDXE,
     &                 NC, WVDIR, ZSill, INDXC, BConn1,BConn2, WIDTH,
     &                 VBUFMX, VBUFTR, TCVBUF )

! -------------------------- Arguments -----------------------------

! In:
      CHARACTER*(*) TOPOGRAPHY_FILE ! specified as model parameter
      INTEGER MD, MBI, MLI_DIM, MLI, MBE, MLE, MC, MLC

!   MD      = max. number of depths,
!   MBI     = 1+max. number of internal basins
!   MLI     = max. number of layers summed over internal basins
!   MLI_DIM = array dimension for AREA and VFRAC
!   MBE     = 1+max. number of external basins
!   MLE     = max. number of layers summed over external basins
!   MC      = 1+max. number of connections
!   MLC     = max. number of layers summed over connections


! Out:
      INTEGER ND
!              - number of depths
      real*8    DEPTH(MD)
!              - Depth values defining layer division
      INTEGER NBI
!              - actual number of internal basins
      real*8    VTOTZ(MBI), ZBOTMI(MBI)
!              - Total initial volume and bottom depth of basins
      real*8    LSHORE(MBI)
!              - Length of shoreline (on 10m fractal dimension)
      INTEGER NLVOPN (MBI)
!              - Number of layers open to boundary,
!                directly or through other basins.
      real*8 VFROPN(MBI)
!              - summed fraction of volume of open layers
      INTEGER INDXI(MBI+1)
!              - layer index limits for internal basins.
      real*8    AREA(MLI_DIM)
!              - Horizontal areas at top of layer
      real*8    BOTTOM(MLI_DIM)
!              - Bottom area within layer
      real*8    VFRAC(MLI_DIM)
!              - volume layer as fraction of total volume
      INTEGER NBE
!              - actual number of external basins
      INTEGER INDXE(MBE+1)
!              - layer index limits for external basins.
      real*8    ZBOTME(MBE)
!              - Lower limit for zill depth of external basins
      INTEGER NC
!              - actual number of connections.
      
      real*8 WVDIR(MC)
!              - direction of the connection from 1 to 2:      
!                0: north south, 90: west - east etc. 
!

      real*8    ZSill(MC)
!              - Sill depths of connections
      INTEGER INDXC(MC+1)
!              - width array index limits for connections.
      INTEGER BConn1(MC), BConn2(MC)
!              - Connected basin numbers.
!                +n: internal basin
!                -n: external basin (only in BConn2?)
      real*8    WIDTH(MLC)
!              - Transport widths of connections,
!                for layers down to ZSILL.
!                Values from I1=INDXC(IC)+1 through I2=INDXC(IC+1)
!                specifies mean width for layer 1 to I2-I1+1
!                in connection number IC.
      real*8 VBUFMX (2,MLC ), VBUFTR(2,MLC)
!           Buffer volume for mass transport: Is filled up from
!           inflows, emptied by outflows.
!           Mass transport into receiving basin
!           does not start until the buffer volume is emptied,
!           cfr. module TRAN_SUB.FOR
      real*8 TCVBUF (2, MC )
!           Time constant for mixing between buffer volume
!           and main volume of each layer in each connection

! ---------------- Local definitions: -----------------------------
      integer IB, NLI, NLE, N1, COUNT, ID, IC, I, K
      real*8 DDP, DPF, D, ADJ, X
      real*8 DDPHIGH, DDPLOW
      real*8 ZBOTTOM_MAX, OLD_AREA, vFrac_Average
      INTEGER IBC (2)


!      INCLUDE 'EutroDim.inc'
      
                      !  WORK AREAS:
      integer lgtWorkArray
      parameter (lgtWorkArray =dimMLI + 2*dimMBI)
      real*8  AREA_BUF (lgtWorkArray), VFRAC_BUF  (lgtWorkArray)

!              - Horizontal areas and volume fraction of layers
!      INTEGER ALLOC_CHK

      INTEGER NF /777/    ! Unit number for file Topography_Report.txt
      CHARACTER*64 TopographyReportFile /'Topography_Report.txt'/

      integer NLC,I1,I2,N
      logical OK

$if defined PRINT_DESC
      character*80 BasinDescr
$endif

$if defined TRACE_CALLS
      WRITE(DiagFile,*) ' >>>>>> Called TOPO  '
$ENDIF


! Ensure array dimension >=specified limit to number of layers:
      MLI = MIN(MLI,MLI_DIM)


!      DEALLOCATE ( AREA_BUF, VFRAC_BUF, STAT = ALLOC_chk)
!      I = MLI + 2*MBI
!      K = 2
!      ALLOCATE ( AREA_BUF(I), VFRAC_BUF(I), STAT = ALLOC_chk)
!      If ( Alloc_chk .ne. 0 ) then
!          WRITE(*,*)
!     &      'Feil ved allokering av arbeidsområde i subrutine TOPOGR',
!     &       Alloc_chk
!          STOP
!      ENDIF

! -----------------------------------------------------------------
! 1. Read topographic data for basins and connections:
      OK= TPRD (TOPOGRAPHY_FILE, TOPO_NR, MBI, NBI, ZBOTMI, LSHORE,
     &           MBE, NBE, ZBOTME, 
     &           MC, NC, WVDIR, ZSill, BConn1, BConn2)
!         - output from call:
!                NBI: number of internal basins < MBI
!                ZBOTMI(1..NBI) max. depth in each internal basin
!                LSHORE(1..NBI) effective length of shore line
!                NBE: Number of external boundary basins, < MBE
!                ZBOTME(1..NBE) max. depth in each external basin,
!                      (defined from connections specified in TPRD)
!                NC: number of connections <MC
!                ZSill(1..NC) max sill depths
!                BConn1, BConn2: connected basin numbers
!           NBx<=MBx (x=I,E) and NC<=MC checked by TPRD,

      IF(.NOT.OK) then
         write(*,*) 'Stops due to error in topographic description'
!         PAUSE 'Press ENTER to terminate program'
         STOP
      endif

$if defined PRINT_INPUT
      call Write_Topography_Input_Data(DiagFile)
$endif
      
! ------------------------------------------------------------------
! 2. Allocate approx. MLI number of layers to NBI basins:
      if(MLI.LT.NBI) then
         write(*,*) 'Stops because number of layers < number of basins'
!         PAUSE 'Press ENTER to terminate program'
         STOP
      endif


! ----------------- Iterating process:
!          User estimates of parameters:
      DDP = MAX(0.1, ABS(DDPMIN))
!              : Minimum thickness approx. DDPMIN,
      DPF = MAX(1.0001,ABS(DPFACT))
!              : Maximum thickness (DPFACT-1)*Depth
      COUNT=0
      INDXI(1) = 0

! CONTROLS BINARY SEARCH FOR DDP VALUES:
      DDPLOW  = 0.0
      DDPHIGH = -1.0 ! Undefined

! -------- start iterating loop:
    5 NLI = 0
      ND = 2

$if defined PRINT_ADAPT_LOG
      WRITE(DiagFile,'('' COUNT '',I5,'' DDP, DPF:'',2G15.7)') count, DDP,DPF
$ENDIF


! --- Find number of layers according to DDP,DPF:
!         ( Depth(1)=0.0 )


      DO IB = 1,NBI
! Upper limit for number and max. depth of equidistant layers:
         N = MAX(1.0D0, MIN( (ZBOTMI(IB)+DDP-0.01)/DDP, 1.0/(DPF-1.0)))
         D = N*DDP  ! Always <= max depth
         if ( D .lt. ZBOTMI(IB) ) THEN
!     Number of layers with increasing thickness, defined by
!     by increasing depth with factor DPF for each layer:
!     Add smallest number >=0 so that   D1*(DPF**N) >= ZBOTMI(IB)
!     or:              N >= log(ZBOTMI(IB)/D1) / log(DPF)
             N = N+ MAX(0, NINT( log(ZBOTMI(IB)/D) / log(DPF) ) )
         ENDIF
         ND = MAX( ND, N+1)
         NLI = NLI + N
         INDXI(IB+1) = NLI
      ENDDO

! ...... Calculate depths and count:
      D = 0.0
      DO I=1,MIN(ND,MD)
         DEPTH(I)=D
         D = MAX ( D+DDP, D*DPF)
      ENDDO

$if defined PRINT_ADAPT_LOG
      WRITE(DiagFile,
     &  '('' NLI='',I5,'' MLI='',I5,'' ND='',I5,'' MD='',I5,3x,A)')
     &              NLI, MLI, ND, MD,' Depth:'
      WRITE(DiagFile,'( 1X, 5G15.7)') (DEPTH(I),I=1,MIN(ND,MD))
$ENDIF

! ----------------------------------------------------------------
! 3. Get area and volume data for each layer,
!         Store temporarily in extended work space:
!      (Necessary to get the adjustment to work properly)
      N1 = NLI
      IF ( NLI .LE. MLI + 2*MBI .AND. ND .LE. MD) THEN
        ID = 0
        ZBOTTOM_MAX = 0.0
        DO IB = 1,NBI
           I1 = INDXI(IB)+1
           I2 = INDXI(IB+1)
           N  = I2-I1+1
           vFrac_Average = 1.0/N
           CALL TPBas( IB, N, Depth, VTOTZ(IB),
     &                 Area_BUF(I1), VFrac_BUF(I1) )

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '<<<<<< returned from TPBAS'
$ENDIF

!        ..... Skip small bottom volume to avoid
!              too strict numeric stability restrictions:

           DO WHILE ( I2.gt.I1 .and. 
     &                VFRAC_BUF (I2) .Lt. vFrac_Average*0.5 )

              VFRAC_BUF(I2-1) = VFRAC_BUF (I2-1) +VFRAC_BUF(I2)

$if defined PRINT_ADAPT_LOG
              WRITE(DiagFile,*)
     &          ' Basin ',IB,'Combine last layer (no. ',N,
     &            ', volume fraction ', VFRAC_BUF(I2),
     &              ') into layer above, to a total volume fration',
     &              VFRAC_BUF (I2-1)
$ENDIF

                     ! Update counters and indices:
             DO I = IB,NBI
               INDXI(I+1) = INDXI(I+1)-1
             END DO
             NLI = NLI-1
             N = N-1
             I2 = I2-1
           END DO

$if defined PRINT_ADAPT_LOG
              WRITE(DiagFile,*)' Finished with checking bottom volume ',
     &          ' for Basin ',IB, ';  NLI=',NLI
$ENDIF

!    accumulates maximum number of layers and max depth
           ID = MAX (N+1,ID)
           ZBOTTOM_MAX = MAX ( ZBOTTOM_MAX, ZBOTMI(IB) )
        END DO

!   Drop any superfluous depth value, and
!   adjust lower depth value to include all bottom layers:

        ND = MIN (ID, ND)
        DEPTH (ND) = ZBOTTOM_MAX
      END IF
!   Zero unused depth values (cosmetics only)
      DO ID = ND+1,MD
        DEPTH (ID) = 0.0
      ENDDO

$if defined PRINT_ADAPT_LOG
      if (N1.ne.NLI) then
         WRITE(DIAGFILE,*) 'Merged bottom volumes:'
         WRITE(DiagFile,'('' reduced NLI from '',I5,'' to '',I5)') N1,NLI
      endif
$ENDIF

! --- Count number of connection layers with this division,
!     and set index limits for WIDTH array:
      NLC = 0
      INDXC(1)=0
      DO IC=1,NC

!     ------ limit number of connected layers to lowest number of
!            actual layers in each basin:
         N = MIN (ND, MD)
         IBC(1) = BConn1(IC)
         IBC(2) = BConn2(IC)
         DO K=1,2
            IB = ABS(IBC(K))
            if ( IBC(K) .gt. 0 ) THEN
                N = min (N, INDXI(IB+1)-INDXI(IB) )
            ENDIF
         ENDDO

         DO WHILE ( ZSILL(IC) .le. DEPTH(N) )
             N = N-1
             IF (N.LE.0 ) EXIT
         ENDDO ! N = lowest layer at least partly connected
         NLC = NLC + N
         INDXC(IC+1) = NLC
      ENDDO

! --- Count number of external basin layers with this division:
      NLE = 0
      INDXE(1)=0
      DO 60 IB=1,NBE
         DO 50 ID = 1,MIN(ND,MD)
            IF( ZBOTME(IB).le.DEPTH(ID) ) THEN
                GOTO 55
            ENDIF
   50    CONTINUE
         ID = ND+1
   55    NLE = NLE + ID-1
         INDXE(IB+1) = NLE
   60 CONTINUE

! ------ Adjust if necessary according to available array dimensions:
!        Requires that no dimensions are exceeded, and tries to
!        use as many of the available internal layers as possible.
      COUNT = COUNT +1
      OK =  ND.le.MD .and. NLI.le.MLI .and. NLE.le.MLE .and.NLC.le.MLC
$if defined PRINT_ADAPT_LOG
      WRITE(DiagFile,'(1X,8a7)')'ND','MD','NLI','MLI','NLE','MLE','NLC','MLC'
      WRITE(DiagFile,'(8i7)') ND,  MD,  NLI,  MLI,  NLE,  MLE,  NLC,  MLC
      WRITE(DiagFile,'('' OK='', L4)') OK
      WRITE(DiagFile,'(1X,4A10)')
     &      'DDP','DDPLOW','DDPHIGH','Count'
      WRITE(DiagFile,'( 1X, 3G16.7,I6)')
     &       DDP, DDPLOW, DDPHIGH, Count
$ENDIF

      IF( ( .not.OK .and. count.lt.100 ) .or.
     &    ( OK .and.  NLI.NE.MLI  .AND. COUNT.LT.100
     &         .and. abs(DDPHIGH-DDPLOW)*1.0e6 .gt. DDPHIGH )  ) THEN
         ADJ = MAX(  FLOAT(ND)/FLOAT(MD),
     &               FLOAT(NLI)/FLOAT(MLI),
     &               FLOAT(NLE)/FLOAT(MLE),
     &               FLOAT(NLC)/FLOAT(MLC)  )

$if defined PRINT_ADAPT_LOG
      WRITE(DiagFile,'('' ADJ='', G15.7)') ADJ
$ENDIF

! ------- Store lower and upper limits to DDP and DPF:
         IF ( ADJ .GT. 1.0 ) THEN
!         ... DDP too small: update lower limit
             DDPLOW = MAX( DDPLOW, DDP)
         ELSEIF ( ADJ .LT. 1.0) THEN
!         ... DDP too large: update or initiate high limit
             IF ( DDPHIGH.gt. 0.0 ) THEN
                 DDPHIGH = MIN( DDPHIGH, DDP)
             ELSE
                DDPHIGH = DDP
             ENDIF
         ENDIF

         IF(ADJ.NE.1.0) THEN
             IF ( DDPHIGH .gt. DDPLOW ) THEN
                DDP = DDPHIGH - (DDPHIGH-DDPLOW)/2.0
                     ! if change is not numerically significant:
                     !    selects DDPHIGH.
             ELSE
                DDP = DDP * ADJ
             ENDIF

$if defined PRINT_ADAPT_LOG
      WRITE(DiagFile,'(1X,3A10)')
     &      'DDP','DDPLOW','DDPHIGH'
      WRITE(DiagFile,'( 1X, 3G16.7)')
     &       DDP, DDPLOW, DDPHIGH
$ENDIF

             GOTO 5
         ENDIF
      ENDIF

! ----- Check result:
      IF(.not.OK) THEN
         WRITE(*,*)'Impossible to adapt depth layers - program error?'
         WRITE(*,*)'Check subroutine TOPO'
         WRITE(*,*)'       NLI,MLI,ND,MD,NLC,MLC:',
     *                     NLI,MLI,ND,MD,NLC,MLC
         ND=0
!         Pause 'Press ENTER to terminate program'
         STOP
      ENDIF

      IF ( NLI .lt. MLI-NBI ) THEN
         WRITE(*,
     &      '('' Warning:'',I5,'' layers < expected lower limit '',I5/
     &        ''   (= max. number of layers - number of basins)'')')
     &            NLI, MLI-NBI
         WRITE(*,'(''    possibly because of other limitation:'')')
         I = 0
         IF (ND .EQ. MD) THEN
            Write(*,'('' Number of depths = max '',I5)' ) MD
            I = I+1
         ENDIF
         IF (NLE .EQ. MLE) THEN
            write(*,'('' Total number of external layers = max '',I5)') MLE
            I = I+1
         ENDIF
         if (nlc .eq.mlc) then
            write(*,'('' Total number of layers connected = max '',I5)') MLC
            I = I+1
         ENDIF
         if (i.eq.0) Write(*,*) ' No obvious reason - check program'

!         Pause 'Press ENTER to continue'

      ENDIF
! ----- Return adjusted values of division parameters:
      DDPMIN = DDP
      DPFACT = DPF

! --------------------------------------------------------------------
! 4. Get transport width and buffer volume description of connections,
!    find number of layers open towards boundary,
!    transfer areas and volumes to external arrays,
!    and free allocated work space:

! ........ Transport width and buffer volume description:
      DO I1 = 1,NbI
         NLVOPN(I1) = 0 
      ENDDO ! Assume all layers closed initially
      
      DO IC = 1,NC
         I1 = INDXC(IC)+1
         I2 = INDXC(IC+1)
         N  = I2-I1+1
         CALL TPCON(IC, N, DEPTH, WIDTH(I1),
     &              VBUFMX(1,I1), VBUFTR(1,I1), TCVBUF(1,IC), OK )

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '<<<<<< returned from TPCON'
$ENDIF

      END DO

      DO I = I2+1,MLC
         DO K=1,2
            VBUFMX(K,I) = 0.0
            VBUFTR(K,I) = 0.0
         ENDDO
      ENDDO

! Find number of layers open towards boundary,
! (i.e. with volume varying according to in/outflows)
      N1 = -1
      COUNT = 0
      DO While (Count.gt.N1) ! Continues as long as the count changes:
         N1 = Count
         DO IC = 1, NC
            I = INDXC(IC+1)-INDXC(IC) ! = num of connected layers
            I1= BConn1(IC)
            I2= BConn2(IC)
            IF ( I2 .le. 0 ) THEN ! Open directly towards boundary
               NLVOPN(I1) = MAX ( NLVOPN(I1), I)
            ELSE ! or through another basin
               I2 = ABS(I2)
               NLVOPN(I1) = MAX ( NLVOPN(I1), MIN (NLVOPN(I2),I)  )
               NLVOPN(I2) = MAX ( NLVOPN(I2), MIN (NLVOPN(I1),I)  )
            ENDIF
         ENDDO
         COUNT = 0 ! Counts total number of open layers
         DO IB = 1, NBI
            COUNT = COUNT + NLVOPN(IB)
                  !mulig for test-form†l: setter alle lag som variable
                  !  nlvopn(IB) = indxi(ib+1) - indxi(ib)
         ENDDO
      ENDDO

      DO IB=1,NBI
               ! Note summed volume fraction of open layers :

         X = 0.0
         DO I = INDXI(IB)+1, INDXI(IB)+NLVOPN(IB), 1
            X  = X + VFRAC_BUF(I)
         ENDDO
         VFROPN(IB) = X

               !  Transfer areas and volumes to external arrays,

         OLD_AREA = 0.0
         DO I = INDXI(IB+1), INDXI(IB)+1,-1
            AREA(I)   = AREA_BUF(I)
            BOTTOM(I) = AREA(I) - OLD_AREA
               OLD_AREA  = AREA(I)
            VFRAC(I)  = VFRAC_BUF(I)
         ENDDO
      ENDDO

               ! Zero unused positions in arrays:

      DO I = NLI+1,MLI_DIM
            AREA(I) = 0.0
            BOTTOM(I) = 0.0
            VFRAC(I) = 0.0
      ENDDO


                     !      DEALLOCATE ( AREA_BUF, VFRAC_BUF)

! ----------- Print topography description for reference -----------

$if defined  PRINT_DESC

      Write(*,'(A/5x,A)') 
     &    'Description of topography data as used by model to file',
     &    TopographyReportFile
      
      OPEN(NF, FILE = TopographyReportFile,ACCESS='SEQUENTIAL')
      WRITE(nf,'(1X,a,I6)')'Topography alternative nr. ',TOPO_NR

      WRITE(nf,'(6A7)')
     &     'NBI', 'MBI', 'NBE', 'MBE', 'NC', 'MC'
      WRITE(nf,'(6i7)')
     &      NBI,   MBI,   NBE,   MBE,   NC,   MC

      WRITE(nf,'(1X,8A7)')
     &    'ND','MD','NLI','MLI','NLE','MLE','NLC','MLC'
      WRITE(nf,'(8i7)')
     &     ND,  MD,  NLI,  MLI,  NLE,  MLE,  NLC,  MLC

      Write(NF,*) 'Number of depths:', ND
      WRITE(NF,'(1x,I4,G15.7)' ) (I,DEPTH(I),I=1,ND)


      WRITE(NF,*) 'Number of internal basins:', NBI
      DO IB=1,NBI
        Write(NF,*)
        BasinDescr = getBasinDescr( IB)
            ! function involves write,
            ! cannot be part of Write statement in Watcom Fortran
        Write(NF,*) BasinDescr
        Write(NF,*) '  VTOTZ:',  VTOTZ(IB)
        Write(NF,*) '  ZBOTMI:', ZBOTMI(IB)
        Write(NF,*) '  Layer index limits:'
        Write(NF,*) INDXI(IB)+1,INDXI(IB+1)
        Write(NF,*) '  Basin topography:'
        Write(NF,'(1x,A12,3A16)')
     &        'Depth', 'Total Area', 'Bottom Area','Volume fraction'
        Write(NF,'(1x,2('' Ã----'',A18,''----´ ''))' )
     &        ' at top of layer ','of layer'
        DO I  = INDXI(IB)+1,INDXI(IB+1)
           ID = I -INDXI(IB)
           WRITE( NF, '( 1X, 4G15.7)' )
     &          DEPTH(ID), AREA(I), BOTTOM(I), VFRAC(I)
        END DO
      END DO


      WRITE(NF,*) 'Number of external basins:', NBE
      DO IB=1,NBE
        BasinDescr = getBasinDescr(-IB)
        Write(NF,'(/A)') BasinDescr
        Write(NF,*) '  Layer index limits:'
        Write(NF,*) INDXE(IB)+1,INDXE(IB+1)
        Write(NF,*) '  ZBOTME:', ZBOTME(IB)
      END DO

      Write(NF,'(/A,I6)')    'Number of connections:', NC
      DO IC=1,NC
        Write(NF,'(/A,I6)') 'Connection nr. ',IC
        Write(NF,*) ' Connected basins:'

        BasinDescr = getBasinDescr(BConn1(IC))
        Write(NF,'(A)') BasinDescr

        BasinDescr = getBasinDescr(BConn2(IC))
        Write(NF,'(A)') BasinDescr

        Write(NF,*) '  WVDIR:', WVDIR(IC)
        Write(NF,*) '  ZSill:', ZSill(IC)
        Write(NF,*) '  Buffer mixing (days):', (TCVBUF(K,IC),K=1,2)

        Write(NF,*)'  Layer index limits for connections:'
        Write(NF,*) INDXC(IC)+1, INDXC(IC+1)
        Write(NF,*)'  Transport width of connected layers:'
        Write(NF,'( 1X,2A10,  ''  ['',A40, ''  ]'')' )
     &       'DEPTH','MEAN','   -------- BUFFER VOLUMES ---------'
        Write(NF, '( 1X,2A10, ''  ['',2A20,''  ]'')' )
     &       '(AT TOP)','WIDTH','      BASIN 1','      BASIN 2'
        Write(NF, '( 21X,4A13 )' )
     &       'VBUFMX','VBUFTR','VBUFMX','VBUFTR'
        DO I = INDXC(IC)+1,INDXC(IC+1)
          WRITE(NF,'( 1X,G10.5, G10.2, 4G13.5 )') 
     &       depth(i-INDXC(IC)),width(I),
     &       ( VBUFMX(K,I),VBUFTR(K,I), K=1,2 )
        ENDDO
      END DO

      CLOSE(NF)
$ENDIF

! --------------  Ferdig med utskrift  -------------------

      END Subroutine


! ============= LOWER LEVEL SUBROUTINES CALLED BY TOPO =================



! =====================================================================
      SUBROUTINE TPBAS(BNUM, NZ, ZLTOP, VLSUM, ALTOP, VLFRAC )

! Called by TOPO in section 3.
! Gets area and volume data for a basin, for given depths

!  Input:
      INTEGER BNUM
!          - Basin number
      INTEGER NZ
!          - Number of specified layers from TP
      real*8 ZLTOP(NZ)
!          - Specified depths of top of layers
!  Output:
      real*8 VLSUM
!          - Total volume over all layers in basin
      real*8 ALTOP(NZ), VLFRAC(NZ)
!          - Calculated areas at given depths,
!            and fraction of total volumes in each layer



! ------------------ local variables: -------------------------
      INTEGER INTV,I
      real*8 DZ,VPSUM,VLAYER
! ---------------------------------------------------------------
$if defined TRACE_CALLS
      WRITE(DiagFile,*) '>>>>>> Called TPBAS'
$ENDIF

! ------- Interpolate in topographic data for specified depths:
!     --- Search topographic data from top:
      INTV = INDXZB(BNUM-1)+2
      DO 20 I = 1, NZ
!     --- Position to interval containing Zspec(I):
   10    IF( ZLTOP(I).GT.ZBASIN(INTV) .AND. INTV.LT.INDXZB(BNUM) ) THEN
             INTV = INTV + 1
             GOTO 10
         ENDIF  !  ZLTOP(I)<=ZBASIN(INTV) or INTV = maximum.
!     --- Interpolate area, and extend volume by trapez integration
!         from ZBasin(Intv) to ZLTOP(I)
         DZ = ZBASIN(INTV)-ZLTOP(I)
         ALTOP(I) = MAX(0.0D0, ZBasinArea(INTV) + DADZ(INTV)*DZ )
!         stores accumulated volumes from bottom up:
         VLFRAC(I) = MAX (0.0D0, VOL(INTV)
     &             + (ZBasinArea(INTV)+ALTOP(I))*DZ/2.0)
!              - area and volume = 0 for depth below deepest bottom
   20 CONTINUE

! ------- Store total volume in VLSUM,
!         and recalculate VLFRAC into fractions of total volume VLSUM:
      VLSUM = VLFRAC(1)
      if (VLSUM.gt.0) then
        VPSUM = 0.0
        DO I=NZ,1,-1
            VLAYER = VLFRAC(I)-VPSUM
            VPSUM = VLFRAC(I)
            VLFRAC(I) = VLAYER/VLSUM
        ENDDO
      endif
      
      END Subroutine
! =================================================================

      SUBROUTINE TPCON( CNUM, NZ, ZLTOP, WIDTH,
     &  VBUFMX, VBUFTR, TCVBUF, OK )

! Called from TOPO in section 4.
! Gets basin numbers and transport widths for a given connection.

! -----------------------------------------------------------------
! Input (request for data):
      INTEGER CNUM
!        - Connection nr.,
      INTEGER NZ
!        - number of layers down to sill
      real*8 ZLTOP(NZ)
!        - Depths at top of model layers
! Output:
      real*8 WIDTH(NZ)
!        - mean width of each layer
      real*8 VBUFMX(2,NZ), VBUFTR(2,NZ)
!        - buffer volumes in each layer, max.lim and transition zone
      real*8 TCVBUF(2)
!       -  mixing time for buffer volume

      logical OK  ! set to false if errors are encountered


! ------------------- local variables ----------------------------
      INTEGER INTV,IZ, K
      real*8 ZHIGH, ZLOW, ZLOW_INTV, ZMID, TRANSP_AREA, W_INTP
      real*8 DZ_Intp_Interval, DZ_Layer
! ----------------------------------------------------------------

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '>>>>>> Called TPCON:'
$ENDIF

! ------- Interpolate in topographic data for specified depths:
!     --- Search topographic data from top:
      INTV = NDXCON(CNUM-1)+2
      DO IZ = 1, NZ
!     --- Depth limits of layer, but only down to sill depth:
         ZHIGH = ZLTOP(IZ)
         IF( IZ .LT. NZ ) THEN
            ZLOW = MIN( ZCON(NDXCON(CNUM)), ZLTOP(IZ+1) )
         ELSE
            ZLOW = ZCON(NDXCON(CNUM))  ! Sill depth
         ENDIF

!     --- Move down to first topographical interval containing ZHIGH(I):
         DO WHILE ( ZHIGH .GE. ZCON(INTV) .AND. INTV .LT. NDXCON(CNUM) )
             INTV = INTV + 1
         ENDDO

!     --- Prepare for accumulating cross-section area and buffer volumes
         TRANSP_AREA= 0.0
         DO K=1,2
            VBUFMX(K,IZ) = 0.0
            VBUFTR(K,IZ) = 0.0
         END DO

!     --- Process topographical layers down to ZLOW, with ZHIGH
!         carrying upper depth limit of layer not yet included in
!         accumulated variables:
         DO WHILE ( INTV .le. NDXCON(CNUM) )
!        ---  Define depth interval in topographical layer INTV
!             within model layer IZ:
            ZLOW_INTV = MIN (ZLOW, ZCON(INTV) )
            DZ_Layer = ZLOW_INTV - ZHIGH  ! Within sill depth
            ZMID = (ZLOW_INTV+ZHIGH)/2.0
!        --- Interpolate to mid-depth and accumulate area and volumes:
            DO K = 0,2
               DZ_Intp_Interval = ZCON(INTV)  - ZCON(INTV-1)
               if (DZ_Intp_Interval.le.0.0) then
                 write(DiagFile,*) 
     &                'ERROR: Interpolating depth interval<=0',
     &                ' in connection no.', CNUM
                 OK = .false.

                 W_INTP = 0.0
               else
                 W_INTP = (  WCON(K,INTV-1)*( ZCON(INTV) - ZMID )
     &                   + WCON(K,INTV)  *( ZMID       -  ZCON(INTV-1) )
     &                     )  / DZ_Intp_Interval
                 W_INTP = MAX( W_INTP, 0.0D0 )
               endif

               IF ( K .EQ. 0 )  THEN
                  TRANSP_AREA = TRANSP_AREA + DZ_Layer *W_INTP
               ELSE
                  VBUFMX(K,IZ) = VBUFMX(K,IZ)
     &                    + DZ_Layer* W_INTP * BUFFER_LENGTH(K,CNUM)
                  VBUFTR(K,IZ) = VBUFTR(K,IZ)
     &                    + DZ_Layer* W_INTP * TRANSITION_LENGTH(K,CNUM)
               ENDIF

            ENDDO

            ZHIGH = ZLOW_INTV

            IF (ZHIGH .GE. ZLOW) EXIT
            INTV = INTV + 1
         ENDDO

!      ------  Model layer IZ processed, store mean width,
!              that is area/depth interval, both within sill depth 
         IF ( ZLOW_INTV .GT. ZLTOP(IZ) ) THEN
             WIDTH(IZ) = TRANSP_AREA / (ZLOW_INTV - ZLTOP(IZ) )
         ELSE 
             WIDTH(IZ) = 0.0
         ENDIF
      ENDDO

! ------ Common to all buffer layers;  minimum limit for mixing times:
      DO K=1,2
         TCVBUF(K) = MAX(0.01D0, MIXING_TIME(K,CNUM) )
      END DO

      END Subroutine



! ======================================================================
      Logical function TPRD (TOPOGRAPHY_FILE, TOPO_NR,
     &                 MBI, NBI, ZBOTMI, LSHORE,
     &                 MBE, NBE, ZBOTME,
     &                 MC, NC, WVDIR, ZSILL, BConn1, BConn2)


! In:
      CHARACTER*(*) TOPOGRAPHY_FILE
      INTEGER TOPO_NR
      INTEGER MBI, MBE,  MC
!              = External dimensions:
!                    max+1 number of basins and connections.
! Out:
      INTEGER NBI,  NBE, NC
!                 = Actual number of internal and external basins,
!                   and number of  connections.
      real*8 ZBOTMI(MBI), ZBOTME(MBE)
!                 = Max. depth of internal basins
      real*8 WVDIR(MC)
!                 = direction from 1 to 2
      real*8 ZSILL(MC)
!                 = sill depth of connections.
      Integer BConn1(MC),BConn2(MC)
!             basin numbers of connections (+n internal, -1 external)
      real*8 LSHORE(MBI)
!                 = effective length of shoreline


! Local:
      CHARACTER WORD(2)*40
      LOGICAL OK, FoundLabel
      LOGICAL INTRNL
!      INTEGER IOCHECK

      integer IB, IZB, MCEFF
  !, NP
      integer IZ, I, KMIN, KMAX, K, NS(2),IS
      real*4    BUF_LGT(2), TR_LGT(2), MIX_T(2)
      real*8 DZ
      real*4 LENGTH_OF_SHORE
      real*4 CONNECTION_DIRECTION
      character*40 label(3)
      integer nlabel

      DATA LABEL(1) /'BASINS'/
      DATA LABEL(2) /'BOUNDARY'/
      DATA LABEL(3) /'CONNECTIONS'/

     
      type (TextInputFileDescriptor) TopoFileDescr      


$if defined  TRACE_CALLS
      WRITE(DiagFile,*) ' >>>>>> called TPRD  '
$ENDIF


				! ********* Open input file for topographic data

      if (TOPOGRAPHY_FILE .eq. ' ') then
         write(*,*) ' TOPOGRAPHY_FILE is not specified'
         TPRD = .False.
         return
      endif
      
	   IF (.not.TextInp_OpenFile ( TopoFileDescr,11,
     &         "Topography/basin description",TOPOGRAPHY_FILE,',|*' ))
     &   goto 101

      write(DiagFile,*)'Reads topography/basin description from file ',
     &                 trim(TOPOGRAPHY_FILE)

$if defined PRINT_INPUT
      call TextInp_DiagFile(TopoFileDescr,DiagFile,.true., .true.)
$endif      
      
             ! Initially parse for only blanks between words:

      CALL TextInp_DelimitersUsed(TopoFileDescr, 1) 


				! ********* Combine external and internal dimension limits,
				!            and initiate empty descriptions:

      BasinName(0)=' '

            ! ------  Basins:
      NBASIN = 0
      NBASEX = 0
      IZB    = 0
      INDXZB(0) = 0

            ! ------  Connections:

      MCEFF= MIN(MC,MCon)
      NCON = 0
      IZ   = 0
      NDXCon(0) = 0

				! ************** Enter data from file:

				         ! Start reading basin data first by default:
      NLABEL = 1
      OK = .true.

      DO
				         ! Read next nonblank line, and fetch first word:
			
			if (.not.TextInp_Get_Line(TopoFileDescr)) exit

   		if(TextInp_WORD(TopoFileDescr,WORD(1))) then

                     ! Check if label; in that case, reset label status 
                     ! and continue looking for name

            FoundLabel = .false.
				Do I = 1,3

	            IF(WORD(1).EQ.LABEL(I)) THEN
	               IF(I.EQ.1.AND.NLABEL.gt.1) THEN
	
	                     ! IF OUT OF ORDER:
	                     ! STILL PROCESSES INPUT, BUT RESULT IS NOT CORRECT
	                     ! (ALL INTERNAL BASINS SHOULD BE ENTERED
	                     !  BEFORE ANY EXTERNAL BASIN):
	
	                  Write(DiagFile,*)' ERROR: Input out of order:'
	                  Write(DiagFile,*)'  ',label(I),
     &                                     ' after ',label(nlabel)
	                  OK = .FALSE.
	
	               ENDIF
                  FoundLabel = .true.
	               NLABEL = I
	            endif

            End Do
            if (FoundLabel) cycle
         else
            OK=.false.
            exit
         endif

				         ! not a label: branch to data input
                     ! according to current label status

		   select case (NLABEL)
		   
			   case (1) ; call ReadBasinDescr()
			   case (2) ; call ReadExternalBasinName
			   case (3) ; call ReadConnectionDescr()
	
	      end select
	      
$if defined TRACE_CALLS
      WRITE(DiagFile,*) '<<<<<< returned to TPRD'
$ENDIF
	      
      end do

            ! should only exit at end of file:
            
      OK = OK .and. TextInp_EOF(TopoFileDescr)

      if (OK) THEN

					! Return dimension and max.depths of basins:

        do i=1,nbasin
            ZBOTMI(i) = ZBasin(indxzb(i))
        end do
        NBI = NBASIN

        do i=1,nbasex
            ZBOTME(i) = 0.0
        end do
        NBE = NBASEX

               ! Note max depth and basin number for connections:

        DO I=1,NCON
           BConn1(I) = NBCon(1,I)
           BConn2(I) = NBCon(2,I)
           ZSILL(I) = ZCON(NDXCON(I))  ! Upper limit

					! Get necessary maximum depth of external basins,
					! used for calculating number of external layers NLE
               ! in calling subroutine TOPO:

           IB = NBCon(2,I)
           IF(IB.LT.0) THEN
              IB = ABS(IB)
              ZBOTME(IB) = MAX( ZBOTME(IB), 
     &                           min( ZBOTMI(BConn1(I)), ZSILL(I) ) )

               !  (Idea: Could multiply ZSILL by factor 1.5
               !   to take withdrawal into account)

           ENDIF
        END DO 

        NC = NCON

      endif

               ! Finished reading topographic description
               ! without errors:

      call TextInp_CloseFile(TopofileDescr)
      
      if (OK) then
         TPRD = .true.
         return
      endif

               ! Reading of topography data aborted by error: 
  101 continue

      call TextInp_StandardMessage(TopofileDescr)
      TPRD = .false.
      return

      contains
!*****************************************************************_
!Internal subroutines embedded in TPRD:
!*****************************************************************_


      Subroutine ReadBasinDescr()

$if Defined TRACE_CALLS
      write(DiagFile,*)'>>>>>> Called ReadBasinDescr'
$endif
      
				!    ---- Internal basin:
				!    -------- first read length of shore-line:
      if(TextInp_Get_Line(TopoFileDescr)) then
         IF (TextInp_Value(TopoFileDescr,LENGTH_OF_SHORE)) then
            IF(NBASIN.LT.MBI) THEN
               LSHORE(NBASIN+1) = LENGTH_OF_SHORE
            ENDIF
         ENDIF
      else
         OK = .false.
         return
      ENDIF

!    --------- then read depth/area data:
         Kmin = IZB+1
         CALL TPRD2(TopofileDescr,TOPO_NR, ZBasin, 1, 1, ZBasinArea,
     &               IZB, MZB,'Depth, area of '//WORD(1), OK )

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '<<<<<< returned from TPRD2'
$ENDIF

!             calc. accumulated volume and change of area pr.depth:

         Kmax = IZB

         if (Kmax.le.MZB) then

          Vol(Kmax) = 0.0

          DO k = Kmax-1, Kmin,-1
            DZ = Zbasin(k+1)-Zbasin(k)
            if (DZ.le.0.0) then ! ensure against program crash,
                                ! will end with error status
               WRITE(DiagFile,*)
     &          'Error in depth/area table for basin no.',NBASIN + 1,
     &          ': non-increasing depth sequence ',Zbasin(k),Zbasin(k+1)
               OK = .false.

               DADZ(k+1) = 0.0
            else
               DADZ(k+1) = ( ZBasinArea(k)-ZBasinArea(k+1) ) / DZ
            endif
            Vol(k)  = Vol(k+1) + ( ZBasinArea(k)+ZBasinArea(k+1) ) * DZ /2.0
          END DO

          DADZ(Kmin) = DADZ(Kmin+1)

         Endif

         INTRNL=.TRUE.
         NBASIN = NBASIN + 1

         IF (NBASIN.EQ.MBasin+1) THEN
            write(ErrorMessage,'(A,I6)')
     &           ' ERROR: Too many internal basins, max. ',MBI
            call TextInp_Message(TopoFileDescr, ErrorMessage)
           OK = .false.
         ENDIF

         BasinName(NBasin) = WORD(1)
         INDXZB(NBasin) = IZB
         
         end subroutine


!--------------------------------------------------------------------------
! ---- External basin: just store in list, and give out number:

         subroutine ReadExternalBasinName

$if Defined TRACE_CALLS
      write(DiagFile,*)'>>>>>> Called ReadExternalBasinName'
$endif
         
         INTRNL=.false.
         NBASEX = NBASEX + 1

         IF (NBASEX.EQ.MBasex+1) THEN
           WRITE(ErrorMessage,'(A,I6)')
     &       ' ERROR: Too many external basins, max. ',MBE
           call TextInp_Message(TopoFileDescr, ErrorMessage)
           OK = .false.
         ENDIF

         BasinName(-NBasex) = WORD(1)
         
         end subroutine


! ----------------------------------------------------------------------
! Internal subroutine: Read description of connection for two basins:
      subroutine ReadConnectionDescr

$if Defined TRACE_CALLS
      write(DiagFile,*)'>>>>>> Called ReadConnectionDescr'
$endif
   
!         ---- Buffer volume length + mixing time:

   20 DO IB=1,2
         CALL TextInp_DelimitersUsed(TopoFileDescr,2) ! blanks or | 
         OK=OK .and. TextInp_ValueSelect( TopoFileDescr,TOPO_NR,BUF_LGT(IB))
         OK=OK .and. TextInp_ValueSelect( TopoFileDescr,TOPO_NR,TR_LGT(IB))
         OK=OK .and. TextInp_ValueSelect( TopoFileDescr,TOPO_NR,MIX_T(IB))

!         --- Read name of second basin on new line:

         IF (IB.EQ.1) THEN
            CALL TextInp_DelimitersUsed(TopoFileDescr,1) ! Only blanks
            if( TextInp_Get_Line(TopoFileDescr)) then
               OK = OK .and. TextInp_WORD(TopoFileDescr,WORD(2))
            else
               return
            endif
         endif
      end do

           ! -------- Check and identify names:
      DO IS = 1,2
         NS(IS) = 0
         DO IB = NBasin, -NBasex, -1
            IF ( WORD(IS).eq.BasinName(IB)) THEN
              ! Stores basin number (>0 for internal, <0 for external)
              ! Note: external basin should always be specified as second
                NS(IS) = IB
                EXIT
            ENDIF
         END DO

         IF( NS(IS) .EQ. 0 ) THEN  ! (Basin name 0 is not used)
            WRITE(DiagFile,*)' ERROR: Undefined basin ',WORD(IS)
            OK = .FALSE.
         ENDIF
      END DO

      IF(NS(1).LT.0 .AND. NS(2).LT.0 ) THEN
         call TextInp_Message(TopoFileDescr,
     &         ' ERROR: two external basins connected')
         OK = .FALSE.
      ELSEIF(NS(1).LT.0) THEN
              ! Switch basin numbers,
              ! external basin always specified as second:
         IB = NS(1)
         NS(1) = NS(2)
         NS(2) = IB
      ENDIF

!         --- Read direction of the connection ( 0.0,...,360.0):
      OK = TextInp_Get_Line(TopoFileDescr)
      if (OK) then
         OK = OK .and. TextInp_Value(TopoFileDescr,CONNECTION_DIRECTION)
      else
         return
      endif

!         --- Read depth/width data for connection:'
      CALL TPRD2(TopoFileDescr, TOPO_NR,
     & ZCon, 0, 2, WCon, IZ, MZCon, 'Dyp, Bredde', OK )

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '<<<<<< returned from TPRD2'
$ENDIF

!         --- Update list of connections:
      NCON = NCON + 1
      IF(NCON.LE.MCEFF) THEN
         WVDIR(NCON) = CONNECTION_DIRECTION
         DO K=1,2
            NBCon (K,NCON) = NS (K)
            BUFFER_LENGTH (K,NCON) = BUF_LGT(K)
            TRANSITION_LENGTH (K,NCON) = TR_LGT(K)
            MIXING_TIME (K,NCON)   = MIX_T(K)
         ENDDO
         NDXCon (NCON)  = IZ
      ELSEIF(NCON.EQ.MCEFF+1) THEN
         WRITE(ErrorMessage,'(A,I6)') ': Too many connections,  max. ', MCEFF
         call TextInp_Message(TopoFileDescr, ErrorMessage)
         OK = .false.
      ENDIF

      end subroutine


      End Function TPRD


! ==============================================================
      subroutine TPRD2( TopoFileDescr, Topo_Nr,
     &                  X, IY, NY, Y, I, M, TERM, OK )

      type (TextInputFileDescriptor) TopoFileDescr      

!  CALLED BY TPRD.
!            Reads one or more lines of X, Y(IY..NY) data
!            terminated by line (Xmax, 0.0) or keyword 'END':
!            Each line can contain  X, Y  or X, Y, Z
!            In the last case, Z defines a scaling factor for Y
!            and should be prepended by delimiter *.
!            The factor is used as multiplicator for Y in this
!            and all subsequent lines until redefined (only within data set)
!            The series should be terminated by value Y = 0.0
!            or line beginning with 'end'
!            X-values should be increasing from zero.
!            each Y can be entered with a number of alternative values 
!            selectable through argument Topo_Nr

      INTEGER TOPO_NR

      INTEGER IY,NY,I,M
      real*8 X(M), Y(IY:NY,M)
      CHARACTER TERM*(*)
      LOGICAL OK

! OK HAS INPUT VALUE .TRUE./.FALSE. - ACCUMULATES ERROR INDICATOR
! FROM SEVERAL CALLS

! --------------- LOCAL VARIABLES -----------------------
      CHARACTER STRING*80
      real*4 XLAST, YFCT, V_Input

      INTEGER IStart, K
      Character*1 DelimiterChar
      LOGICAL LINEOK
      LOGICAL Profile_Ended

! -------------------------------------------------------

$if defined TRACE_CALLS
      WRITE(DiagFile,*) '>>>>>> Called TPRD2:'
$ENDIF

      Xlast = -1.e-36
      YFCT = 1.0
      Profile_Ended = .false.
      IStart = I
      
      DO
         
         CALL TextInp_DelimitersUsed(TopoFileDescr,1) ! BLANKS only
         if (TextInp_Get_Line( TopoFileDescr)) then
            IF (TextInp_Word(TopoFileDescr, STRING) ) THEN
               if (trim(STRING) == 'END') Exit
!               if (indexCaseNeutral(STRING,'END') == 1) Exit
            else
               OK = .false.
            endif
         else
            OK = .false.
            exit
         endif

               ! More data: Check available array space

         I = I + 1
         IF(I.EQ.M+1) THEN
            WRITE(DiagFile,*) 'Too many values entered into ',TERM
            OK = .false.
         ENDIF

         IF (I.GT.M) CYCLE    ! keep reading lines,
                              ! but do not process and store info

               ! Arrays large enough:
               ! read first word as single real*8 value X(I):
         
         LINEOK = .TRUE.      ! assumed
         call TextInp_Repeat(TopoFileDescr)
         if (.not.TextInp_Value(TopoFileDescr,V_Input)) THEN
            LINEOK = .false.
         else
            X(I) = V_Input      
            if(X(I).le.Xlast)then
               Write(ErrorMessage,*) ' Depths not sorted: ...', XLAST,X(I),'...'
               call TextInp_Message(TopoFileDescr, ErrorMessage)
               OK = .false. ! Line may be OK; but is out of order
            endif

            Xlast = X(I)

               ! read from sequence of Y values delimited by |,
               ! picking value number Topo_Nr:            

            DelimiterChar =" "
            call TextInp_DelimitersUsed(TopoFileDescr,3) ! (comma, | , * or blank)
            DO K = IY,NY
               IF (TextInp_EOL(TopoFileDescr) .or. DelimiterChar == "*") then
                          ! no more input of Y values,
                          ! fill in with last y value
                  IF (K.GT.IY) THEN
                     Y(K,I)=Y(K-1,I)

                  else
                     call TextInp_Message(TopoFileDescr, 
     &                          "expected at least two values before '*'")
                     LINEOK = .false.
                  endif

               else
                  LINEOK = TextInp_ValueSelect ( TopoFileDescr,Topo_Nr,V_Input)
                  Y(K,I) = V_Input
                  DelimiterChar= TextInp_DelimiterFound(TopoFileDescr)
               ENDIF
            ENDDO
            
$if defined PRINT_INPUT            
            write(DiagFile,*)'has read X and Y; LINEOK=', LINEOK,' I=',I,' M=',M
            if (I.le.M) then
               write(DiagFile,*)'X(I)=',X(I),' Y(I,:)',(Y(I,K),K=IY, NY)
            else
               write(DiagFile,*)'I > M'
            endif            
$endif          

            IF ( LINEOK .and. (.not. TextInp_EOL(TopoFileDescr)) ) then

$if defined PRINT_INPUT            
               write(DiagFile,*)'DelimiterChar  ',DelimiterChar
$ENDIF            

               IF (DelimiterChar.EQ."*") THEN
                     ! delimiter * should be followed by multiplication factor:
                  call TextInp_DelimitersUsed(TopoFileDescr,1)
                  LINEOK = TextInp_Value(TopoFileDescr,YFCT)
               ELSE
                  write(DiagFile,*)
     &               'Expected delimiter * in front of factor'
                  LINEOK = .false.

               ENDIF
            ENDIF
            
         !    Skip all records after first with area (or width) = 0   
            if ( Profile_Ended ) then             
               I=I-1
            else
               if ( Y(IY,I).le.0) Profile_Ended=.true. 
               DO K=IY, NY
                  Y (K,I) = Y(K,I) * YFCT
               ENDDO
            Endif
         ENDIF   

         IF (.NOT.LINEOK) THEN
            call TextInp_StandardMessage(TopoFileDescr)
            OK = .false.
         ENDIF

      ENDDO

      if (I.eq.Istart) then
         call TextInp_Message(TopoFileDescr,  'ERROR: empty profile')
      endif

$if defined PRINT_INPUT            
      write(DiagFile,*)'exits TPRD2 with I=',I,'; start value of I was ', IStart
$endif

      end subroutine


! ====================================================================
! get basin description (IB>0 for internal, <0 for external basins)

      character*80 function getBasinDescr(IB)

      character BasinNumber*6
      integer IB

      write(BasinNumber,'(I6)') abs(IB)
      if (IB.gt.0 .and. IB .LE. NBasin) then
         getBasinDescr = '=== Internal basin nr. '// BasinNumber
     &                    // ' name:' // BasinName(IB)   
      elseif (IB.lt.0 .and. IB .GE. -NBasex) then
         getBasinDescr = '=== External basin nr. '// BasinNumber
     &                    // ' name:' // BasinName(IB)   
      else
         getBasinDescr = '=== Undefined basin, nr.'// BasinNumber
         return
      endif

      end function


! =====================================================================
! Write topography description as read from topography file:

      subroutine Write_Topography_Input_Data(nf)
      integer :: nf  ! Unit number to write to

      integer :: nb, nl, nc

      write(nf,*)' Number of basins: NBASIN=',NBASIN, '  NBASEX=', NBASEX

      do nb = 1,NBASEX
         write(nf,*)'External Basin no. ',nb,
     &              ': BasinName(',-nb,')=',trim(BasinName(-nb))
      enddo

      do nb=1, NBASIN
         write(nf,*)
         write(nf,*) 'Internal Basin; BasinName(',nb,') = ',trim(BasinName(nb))
         write(nf,*) 'Table index INDXZB from ',indxzb(nb-1)+1,' to ', indxzb(nb)
         write(nf,'(A6,4A14)') ' Row ', ' ZBasin ',' ZBasinArea ',' Vol ',' DADZ '
         do nl = indxzb(nb-1)+1, indxzb(nb)
            write(nf,'(I6,4G14.4)') nl,ZBasin(nl),ZBasinArea(nl),Vol(nl),DADZ(nl)
         end do
      enddo

      write(nf,*)' Number of connections NCON=', NCON
      do nc=1, NCON
         write(nf,*)
         write(nf,*)'Connection no. ', nc
         write(nf,*)'Basin numbers (NBCON(:,',nc,'):', (NBCON(nb,nc),nb=1,2)
         write(nf,'(3A20)')  'BUFFER_LENGTH', 'TRANSITION_LENGTH', 'MIXING_TIME'
         write(nf,'(3G20.5)')
     &   ( BUFFER_LENGTH(nb,nc), TRANSITION_LENGTH(nb,nc), MIXING_TIME(nb,nc),
     &     nb=1,2 )
         write(nf,*)'Table index NDXCON from ',NDXCon(nc-1)+1,' to ',  NDXCon(nc)
         write(nf,'(A6,2A14)') 'Row ', ' ZCON ',' WCON '
         do nl = NDXCon(nc-1)+1, NDXCon(nc)
            write(nf,'(I6,4G14.4)') nl,ZCON(nl),(WCON(nb,nl), nb=0,2)
         end do

      end do

      end subroutine

      
		end Module
