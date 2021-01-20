      Module fx_SurfExch
      use fx_OxygSat, only: Oxygen_Saturation
      use ModelParam_Physics, only: OXBUBL, OXSFAC

      implicit none
      

      CHARACTER*132 :: METFILE_NAME = "Met_Data.txt"


C =====================================================================
C Eutrophication model - inner Oslo fjord
C File  SURFEXCH.FOR
C Meteorological forcing of model.
C Includes both input of data and subroutines for influence
C of meteorology conditions on model through surface
C =====================================================================

$define met_first_version
$undefine test_metinp
$undefine test_lgtrad
$undefine test_radABS

$if defined test_radABS
$define test_radABS_detail
$endif


      contains

C =================================================================
C                    Input of metorological data
C =================================================================


C -----------------------------------------------------------------
C                Initiate input of meteorological data
C -----------------------------------------------------------------

C ============ FORCES REINITIATING IN NEXT CALL TO METINP
      SUBROUTINE METINI
      
      save

$if defined met_first_version
      INTEGER METFILE, ERRNUM
      LOGICAL FORCE_INIT
      COMMON /MET_CONTROL/METFILE, ERRNUM, FORCE_INIT
      metfile = 101
      FORCE_INIT = .true.
$else
      INTEGER I
      LOGICAL OCCUPIED
      INTEGER METFILE_LIST_UNIT, METFILE, ERRNUM
      LOGICAL FORCE_INIT
      COMMON /MET_CONTROL/METFILE, ERRNUM, FORCE_INIT

      COMMON /MET_CONTROL/ METFILE, ERRNUM, FORCE_INIT
     &
      METFILE_LIST_UNIT = 0
      do for i=101, 700
         inquire (unit=i, opened = occupied )
         if (.not. occupied )then
           if( METFILE_list_unit.eq.0 ) then
              Metfile_list_unit = i
              open(UNIT=i, files='METFILES.', IOCHECK=ERRNUM )
              if ( ERRNUM .ne. 0 ) then
                   write(*,'('' Error '',I5,2A)') ERRNUM,
     &                ' when opening list of metorological files',
     &                ' in METFILES.'
                   Stop
              endif
           else
              METFILE = i ! will be used for actual datafiles below
              FORCE_INIT = .true.
              Exit
           endif
         endif
      END

      if ( metfile_list_unit.eq.0) then
         STOP 'No free unit numbers for METFILES'
      Endif
$endif

      End Subroutine



C =========== Read metorological data from source file:
      SUBROUTINE METINP( T,
     &           AIRTMP, AIRP, HUMGM3, CLOUDS, PRECIP, WINDN, WINDE )
      
      SAVE 

C IN:
      real*8 T
C OUT:
      real*8 AIRTMP  ! Air temperature
      real*8 AIRP    ! Air pressure, mb
      real*8 HUMGM3  ! Humidity in g/cm3
      real*8 CLOUDS  ! CLOUD COVER 0-8, 9=fog
      real*8 PRECIP  ! Precipitation in m/s
      real*8 WINDN   ! North component of wind speed
      real*8 WINDE   ! East component of wind speed
C                       linear interpolation in 6 hour data.


C -------- Input file:
      CHARACTER*80 FORMAT

C -------- File status shared with METINI:
      INTEGER METFILE, ERRNUM
      LOGICAL FORCE_INIT
      COMMON /MET_CONTROL/METFILE, ERRNUM, FORCE_INIT

c -------- Local file status:
      LOGICAL OPEN_STATUS /.FALSE./
      integer nrec   ! number of records entered
      logical First_record
      integer ERR_COUNT, LAST_ERROR

C -------- data buffer: -----------
      INTEGER NV,MT
      PARAMETER (NV=10, MT=10)

           ! Input line buffer:
      character*80 InputLine(2)
      integer Current, Previous

C Primary buffer receiving input values for one observation:
      real*8 INP_BUF(NV)

C Buffers storing data to be used for sequence of observations:
      real*8 TIME       (0:MT)
      real*8 T_PREV
      real*8 PRECIP_SUM, PRECIP_PREV
      real*8 AIR_TEMP    (MT)
      real*8 HUMID_GM3   (MT)
      real*8 CLOUD_COVER (MT)
      real*8 WIND_NORTH  (MT)
      real*8 WIND_EAST   (MT)
      real*8 AIR_PRESSURE(MT)
C Last 6 buffers collected as one matrix for interpolation:
      real*8 DATABUF(0:MT,6)
      EQUIVALENCE (DATABUF(1,1),AIR_TEMP   (1))
      EQUIVALENCE (DATABUF(1,2),HUMID_GM3  (1))
      EQUIVALENCE (DATABUF(1,3),CLOUD_COVER(1))
      EQUIVALENCE (DATABUF(1,4),WIND_NORTH (1))
      EQUIVALENCE (DATABUF(1,5),WIND_EAST  (1))
      EQUIVALENCE (DATABUF(1,6),AIR_Pressure (1))

      
      INTEGER NT         ! Number of observation currently stored
      INTEGER LAST       ! Index to first observation in series
      INTEGER PREV       !          previous observation
      INTEGER NEXT       !          new observation
                         !      (ring buffer)

      real*8 T_ZERO        ! Time in data file corresponding to 
                         ! model time T = 0
      real*8 T_BASE        ! 365*n, shifts time to reuse data for one year
      INTEGER NYEAR_BASE ! used in setting T_BASE
      real*8 DT
      real*8 PI
      PARAMETER ( PI=3.141593 )
      real*8 RADIANS_PER_DEGREE, ANGLE
      PARAMETER ( RADIANS_PER_DEGREE = PI/180. )
      real*8 TX            ! effective temperature for interpolation

      INTEGER I

      LOGICAL FIXED_FORMAT, SKIPLINE
      INTEGER Num_Header_Records

C ----------- opens file with meterological data if required -----------
      INQUIRE(FILE=METFILE_NAME, OPENED = OPEN_STATUS)

      IF ( force_init .OR. .NOT. OPEN_status ) THEN
         CLOSE(METFILE, IOSTAT = ERRNUM )

         OPEN(METFILE, FILE=METFILE_NAME, IOSTAT = ERRNUM,
     &          status='OLD',ACTION='READ')
         IF ( ERRNUM.NE.0 ) THEN
            WRITE(*,
     &         '('' Error status '',I6,'' on opening file ''/1x,A)')
     &                            ERRNUM, trim(METFILE_NAME)
             STOP
         ELSE
            WRITE(*,*)'Reading meteorological data from file "',trim(METFILE_NAME),'":'
         ENDIF

         NREC = 0

                    ! read until line not starting with # is found,
                    ! this is assumed to contain format string
         FORMAT='#'
         DO WHILE (FORMAT(1:1).eq.'#')
            READ (METFILE, '(A)', IOSTAT = ERRNUM) FORMAT
            NREC = NREC +1      ! Counts number of accepted input lines
            write(*,*) FORMAT   ! FORMAT may be asterisk * to specify free format
         END DO
         FIXED_FORMAT = (FORMAT(1:1).ne.'*')

         IF ( ERRNUM.NE.0 ) THEN
            WRITE (*,'(" Error status ",I6,
     &                 " when reading data format string from file "/1x,A)')
     &                         ERRNUM, trim(METFILE_NAME)
            STOP
         ENDIF

         Num_Header_Records = NREC
             ! notes number of lines to skip after Wrap-around with REWIND

         NT = 0
         LAST = 0
         
         TIME(LAST)=0.0
         ERR_COUNT = 0
         T_base = 0.0
         First_record = .true.


C        ....... INITIALIZE PRECIPITATION ACCOUNT
         PRECIP_SUM = 0.0
         PRECIP_PREV= 0.0
         T_PREV = 0.0
         force_init = .false.
         Current = 1
         InputLine(1)=' '
         InputLine(2)=' '
      ENDIF

               ! Fill up buffer until current time is in range of observations
               ! and at least two records have been read (possibly through 
               ! recycling the file)

      DO WHILE ( T .GT. TIME(LAST) .or. NT .LE. 1 )

         READ (metfile,'(A)',IOSTAT=ERRNUM) InputLine(Current)
         NREC = NREC +1

         IF ( ERRNUM .LT. 0 ) THEN  ! End of file,
            IF( NT.eq.0 ) then
               WRITE(*,*) " No data found in file ", trim(METFILE_NAME)
               STOP
            ELSE
C              ..... reread file, shift time scale by number of years
C                    (always at least one 365 day year:)
               NYEAR_base = MAX( 1, NINT((T-T_base)/365.0))
               
               T_base = T_base + 365.*NYEAR_base
               REWIND (METFILE, IOSTAT = ERRNUM)

                     ! read past header records counted above, including format record:
               NREC = 0
               DO WHILE (NREC.lt.Num_Header_Records)
                  READ (METFILE, '(A)', IOSTAT = ERRNUM) InputLine(Current)
                  NREC = NREC +1
               END DO
                     ! read first data line:
               READ (metfile,'(A)',IOSTAT=ERRNUM) InputLine(Current)
               NREC = NREC +1
               ERR_COUNT = 0
            ENDIF
         ENDIF

         SKIPLINE = (InputLine(Current).eq.' ')
         IF (.not. SKIPLINE) then
            DO I = 1, Len(InputLine(Current))
               IF (InputLine(Current)(I:I).ne.' ') then
                  SKIPLINE = ( (InputLine(Current)(I:I).eq.'!')
     &                    .or. (InputLine(Current)(I:I).eq.'#'))
                  exit  ! after first non-blank
               END IF
            END DO
         endif
         
         if (SKIPLINE) CYCLE

         IF ( ERRNUM  .GT. 0 ) THEN
             WRITE(*,'(2(A,I6),2A)') " Error status ", ERRNUM,
     &                 " on reading line number ", NREC,
     &                 " from file ", trim(METFILE_NAME)

         ELSE
            NEXT = MOD( LAST, MT) + 1   ! Index in ring buffer
            if (FIXED_FORMAT) THEN
               READ (InputLine(Current), FORMAT, IOSTAT=ERRNUM) INP_BUF
            else
               READ (InputLine(Current), *, IOSTAT=ERRNUM) INP_BUF
            endif
            NREC = NREC +1

C Expects following input:
C(* : used)
C    1. Date YYMMDD (or YYYYMMDD - only for data source description)
C    2. Hours       (not used)
C *  3. T-value, from start of first year (days)
C                to be calculated in advance)
C *  4. Precipitation in mm since last input
C *  5. Air pressure
C *  6. Air temperature
C *  7. Relative humidity in %
C *  8. Cloud cover 0-8 (9=fog)
C *  9. Wind direction in 360-degr.
C * 10. Wind speed in m/s
            
            IF ( ERRNUM  .NE. 0 ) THEN
               ERR_COUNT = ERR_COUNT + 1
               Last_error = errnum

            ELSE

               IF (ERR_COUNT.GT.0) THEN
                  
                  WRITE(*,*)
     &                ' error in interpretation of ', ERR_COUNT,
     &                ' lines up and including line number ', NREC
!                  Pause
                  WRITE(*,*) ' in file ', trim(METFILE_NAME),
     &                ' , last error status=', Last_error
!                  Pause

                  WRITE(*,* ) ' Format:', trim(FORMAT)

                  WRITE(*,'(A/A)') ' Content of last line:',
     &                  InputLine(Previous)
                  ERR_COUNT = 0

               endif

C           ------ store values to be used:
C              ...... Time to compare with model time,  
C               ...... first year in file is always equal to first year in model
C               ...... during initial cycle: 
               
               if (First_record) then
                  T_Zero = inp_buf(3)-mod(inp_buf(3),365D0)
                  First_record = .false.
               endif
               
               TIME (NEXT) = INP_BUF(3) -T_Zero + T_Base
                     ! Accumulate precipitation in mm:
               PRECIP_SUM = PRECIP_SUM + INP_BUF(4)
               AIR_TEMP  (NEXT) = INP_BUF(6)
C              ...... transform relative into absolute humidity:
               HUMID_GM3 (NEXT) = INP_BUF(7)/100.
     &                         * 4.8 * (AIR_TEMP(NEXT)/109.7 + 1)**7.594
               CLOUD_COVER(NEXT) = INP_BUF(8)
C              ...... transform wind into components:
               ANGLE = INP_BUF(9)*RADIANS_PER_DEGREE
               WIND_NORTH(NEXT) = INP_BUF(10)*COS(ANGLE)
               WIND_EAST (NEXT) = INP_BUF(10)*SIN(ANGLE)
               AIR_PRESSURE(NEXT) = INP_BUF(5)

C           ------ include new obs. as last in sequence:
               NT = MIN(MT,NT+1)
               PREV = LAST
               LAST = NEXT
            ENDIF
         ENDIF

                ! update pointers to preserve previous line for 
                ! possible error mssages (see above)

         Previous= Current
         Current = 3 - Current  ! swap 1 <---> 2
      ENDDO

$if defined test_metinp
      write(999,'(A,2(/1x,A))') "InputLines (Previous and Current): ",
     &                       InputLine(Previous),InputLine(Current)
$endif

C ----------- at least one record entered, and T =< T_LAST.
C ----------- establish output values at time T:
C ........ Interpolate momentaneous values:
      DT = TIME(LAST) - TIME(PREV)
      if (DT .LE. 0.0 ) THEN
          WRITE (*,'(2(A,G12.3)/A)') 
     &       " Time series out of order; time ",
     &         TIME(LAST), " after ", TIME(PREV),
     &        " in meteorological data"
!          PAUSE
      ENDIF

C       INTERPOLATE:
      TX = MAX(T,TIME(PREV))  ! NO EXTRAPOLATION
      DO I = 1,6
         DATABUF(0,I) =  (   DATABUF(PREV,I)*(TIME(LAST)-TX)
     &                 +   DATABUF(LAST,I)*(TX-TIME(PREV))
     &               )/ DT
      END DO


C ------------------------------------------
C Export values, and prepare for next call:
C ------------------------------------------

C ........ Transfer interpolated values to subroutine arguments:

      AIRTMP = DATABUF(0,1)
      HUMGM3 = DATABUF(0,2)
      CLOUDS = DATABUF(0,3)
      WINDN  = DATABUF(0,4)
      WINDE  = DATABUF(0,5)
      AIRP   = DATABUF(0,6)


          ! Return precipitation of current interval (on per day basis:

$if defined test_metinp
      write(999,'(1x,5A13)')
     & 'PRECIP_SUM', 'PRECIP_PREV', 'T', 'T_PREV', 'TIME(LAST)'
      write(999,'(1x,5G13.6)')
     & PRECIP_SUM, PRECIP_PREV, T, T_PREV, TIME(LAST)
$endif

           ! Correct accumulated precipitation
           ! by the amount that was applied over previous
           ! integration step by the model:
           ! (returned rate*integration step)

      PRECIP_SUM = MAX(0.0D0, PRECIP_SUM - PRECIP_PREV*(T-T_PREV) )
      
           ! Set a precipitation rate that will correspond to
           ! accumulated precipitation being applied 
           ! up to next input time, or at least 0.1 day
           ! (should be changed?)
      PRECIP = PRECIP_SUM / MAX(0.1D0,TIME(LAST)-T )    ! mm/day
           ! Save current precipitation rate:
      PRECIP_PREV = PRECIP

           ! Convert precipitation rate to m/s:
      PRECIP = PRECIP /1000. / 24./3600.         ! m3/s/m2

           ! save current time for next call:
      T_PREV = T

$if defined test_metinp
      write(999,'('' metinp:'',A12,6A7)')
     &    'T', 'AIRTMP', 'HUMGM3', 'CLOUDS', 'PRECIP', 'WINDN', 'WINDE'
          write(999,'(8x,f12.2,2f7.2,f7.1,3f7.1)')
     &              T, AIRTMP, HUMGM3, CLOUDS, PRECIP, WINDN, WINDE
$endif

      End Subroutine



C =================================================================
C Solar radiation, also fixes integration step:
      SUBROUTINE LGTRAD( T, CLOUDS, DTMAX, DAYDIV,
     &      SINHS, SUNRAD, DIFRAD, TSTEP )
      
      
      SAVE

C IN:
      real*8 T
      real*8 CLOUDS ! CLOUD COVER 0-8, 9=fog
      real*8 DTMAX  ! Maximum time-step in unit days set outside subr.
      real*4 DAYDIV ! Minimum number of steps per. daylight period
C OUT:
      real*8 SINHS   ! Sinus(h), where h=height of sun in radians
      real*8 SUNRAD  ! Solar radiation pr. horizontal area
      real*8 DIFRAD  ! Diffusive radiation pr. horizontal area.
C         Mean values for the time-step about to be integrated.
C         Final length of time-step is set to ensure
C         a specified resolution of the daylight period.
      real*8 TSTEP   ! Final time-step to be used for integration


C ------------------ local variables ----------------------
      real*8 SUN_CONST  / 1350.   /     ! W/m2
      real*8 PI, PI2
      PARAMETER ( PI  = 3.141593 )
      PARAMETER ( PI2 = PI*2)

      real*8 T_SUN , T_SHIFT
      real*8 GLOBRAD, INCL, LATID, RAD_MAX, R_CLOUD
      real*8 OMEGA, DAY_LENGTH, HALF_DAY
      real*8 DTMAX_R, DTMAX_LIGHT, T1, T2, DT_R
      real*8 COS_MEAN


C ----------------------------------------------------------------
C Initial/constant values:

      LATID = 60.*PI/180  ! Northern latitude 60 o
      DTMAX_R = DTMAX*PI2 ! Upper limit on time-step in radians

C   ......  Solar time in radians relative to closest 12:00h at time T
      T_SUN = (MOD(T,1.0D0)-0.5)*PI2

C   ......  Solar inclination and length of closest daylight period
C   ...... Inclination of sun in radians:
      INCL =  -(23.45*PI/180)
     &            * COS( PI2*MOD( (T+10.5)/365 , 1.0D0 ) )
C   ...... Length of daylight period:
      OMEGA = MAX( -1.0D0, MIN( 1.0D0, -TAN(INCL)*TAN(LATID) ) )
C              = COS(t), t=solar time in radians at sunrize and sundown
C   ...... Daylight-period related intervals in radians:
      HALF_DAY = ACOS(OMEGA)
      DAY_LENGTH= HALF_DAY*2.0
      DTMAX_LIGHT = DAY_LENGTH / MAX( 1.0, ABS( DAYDIV) )
C           (always divides day period in at least DAYDIV intervals)

C   ...... Shifts time reference to following day
C          if T is just before end of night:
      T_SHIFT = 0.0
      DO WHILE (T_SHIFT .le. 1.0 )
         IF (T_SUN .LE. (T_SHIFT+1)*PI2-HALF_DAY-DTMAX_LIGHT/4.0 ) EXIT
         T_SHIFT = T_SHIFT + 1.0
      END DO
      T_SUN = T_SUN - T_SHIFT*PI2

C ------------- Select time step and set mean daylight sun height,
C               measured as mean of COS(t), t=solar time in radians.
C               This gives give mean value of SINHS used in RAD_MAX,
C               since SINHS is linear function of COS(t).
C                For dark periods, COS=OMEGA is used, since this
C               corresponds to SINHS = 0, that is, no light.
C
      IF (T_SUN .LT. -HALF_DAY - MIN(DTMAX_R,DTMAX_LIGHT/4.0 ) ) THEN
C     ...... Sufficiently long time before sunrize to implement
C            dark period, limited by max. step or until sunrize.
         DT_R = MIN( DTMAX_R, -HALF_DAY-T_SUN )
         COS_MEAN = OMEGA  ! Gives SINHS exactly = 0

      ELSEIF (T_SUN .LT. HALF_DAY ) THEN
C     ...... Mainly daylight period:
C       ..... Daylight from T1 to T2:
         T1 = MAX ( - half_DAY, T_SUN)
         T2 = MIN( half_day+DTMAX_LIGHT/4.0,
     &             T_SUN+DTMAX_R, T1+DTMAX_LIGHT*1.25 )
C                    ( T2>=T1 because T_SUN >= -HALF_DAY-DTMAX_R
C                      according to first check on T_SUN )
$if defined test_lgtrad
      if ( T2 .lt. T1 ) then
         write(*,*) 'subroutine LGTRAD, T2=', T2, '< T1=',T1
!         Pause
         T2=T1
      endif
$endif

C       .... Total time interval
C             (possibly including short dark period before T1)
         DT_R = T2 - T_SUN
C       .... Mean COS-value, using COS=OMEGA over dark period:
         IF (DT_R .GT. 0.0 ) THEN
            COS_MEAN = (OMEGA*MAX(0.0D0,T1-T_SUN)+SIN(T2)-SIN(T1))/DT_R
         ELSE
            COS_MEAN = COS (T_SUN)
         ENDIF

      ELSE

            !     ....... Dark period following sundown,
            !             limit timestep to approx. next sunrise:

         DT_R = MIN( DTMAX_R, PI2 - HALF_DAY - T_SUN )
         COS_MEAN = OMEGA

      ENDIF


C    ...... Time-step returned in unit (days):
      TSTEP = MIN( DT_R/PI2, DBLE(DTMAX))

C    ......... mean effective sun height as sin(h), =0 in dark periods
      SINHS = SIN(INCL)*SIN(LATID) + COS(INCL)*COS(LATID)*COS_MEAN
      SINHS = MAX(0.0D0, SINHS)
C           limit 0.0 imposed to avoid problems when used in functions.
C           Daylight values used for albedo outside this subroutine.

$if defined test_lgtrad
      WRITE(999,'(6(1x,A10),1x,A15)' )
     &   'LGTRAD: T','TSTEP','OMEGA','COS_MEAN','INCL','LATID','SINHS'
      WRITE(999,'(6(1x,6F10.6),1x,G15.7)' )
     &            T, TSTEP, OMEGA, COS_MEAN, INCL, LATID, SINHS
$endif

$if defined test_lgtrad
      if (sinhs .ge. 1.0) then
         write(*,*) 'subroutine LGTRAD, sinhs=', sinhs
!         Pause
         sinhs = 1.0
      endif
$endif

C    ....... theoretical mean solar radiation pr. horizontal area,
C            without absorbtion and scattering in atmosphere:
      RAD_MAX  = SUN_CONST* SINHS
      IF (CLOUDS.EQ.9) THEN
         R_CLOUD = 7.0
      ELSE
         R_CLOUD = CLOUDS
      ENDIF
C    ....... empirical global radiation (direct+diffuse):
      GLOBRAD = (0.67-0.0011*(MIN(8.0D0,R_CLOUD)**3))*RAD_MAX
C            radiation after sundown is neglected
C            fog treated as cloud cover 7, approx. according to data

C    ....Split into:
C       .... direct light
C           .... in clear skies:
      if ( sinhs .ge. 0.005 ) then
         SUNRAD = (0.67*RAD_MAX)*(0.45+0.52*exp(-0.14/sinHS))
      else
         SUNRAD = (0.67*RAD_MAX)*(0.45)
      endif
C           .... proportionally reduced with cloud cover:
      SUNRAD = SUNRAD*(1.0-CLOUDS/8.0)

C       .... residual as diffusive radiation:
      DIFRAD = GLOBRAD - SUNRAD
C NOTE! Radiation terms as W/m2 horizontal surface just above water

$if defined test_lgtrad
c       if (t .gt. 100) write(999,'(A/F8.2,7G8.2)')
c     &  ' T   ,T_SUN, INCL, SINHS, RAD_MAX, GLOBRAD, SUNRAD, DIFRAD',
c     &    T, T_SUN, INCL, SINHS, RAD_MAX,GLOBRAD,SUNRAD,DIFRAD
$endif


      End Subroutine


C =================================================================
C ======== OXYGEN EXCHANGE THROUGH SURFACE
C ======== added to derivative of oxygen contents in upper layer
C ======== of each basin, assumed to have been initialized
C ======== by transport calculation

      SUBROUTINE OXEXCH( WNDSPD, T, S, AREA,
     &                   ND, DEPTH, NLI, INDXI, NBI, VLAYER,
     &                   OXCONC, OXYGDV, OXYGMP, OXSAT, OXEXCF )
      
      
      SAVE


      real*8 WNDSPD          !Effective wind speed in m/s
      INTEGER NLI          !Number of layers
      real*8 T(NLI)          !temperature in deg C
      real*8 S(NLI)          !salinity in o/oo
      real*8 AREA(NLI)       !Horizontal area at top of layers (m2)
      INTEGER ND
      real*8 DEPTH (ND)      ! Depth at layer interfaces (m)
      INTEGER NBI          !Number of basins
      INTEGER INDXI(NBI+1) !layer index limits for internal basins.
      real*8 VLAYER (NLI)    !Volume in each layer
      real*8 OXCONC (NLI)    !Oxygen concentration (ml/l= liter/m3)
      real*8 OXYGDV (NLI)    !Oxygen content time derivative (ml/l/day)
      real*8 OXYGMP (NBI)    !Oxygen input from atmosphere (l/day)
      real*8 OXSAT  (NBI)    !Oxygen saturation concentration
      real*8 OXEXCF          !Exchange coefficient (m/day)


C ------------------------ local variables --------------------------
      INTEGER IB, LBASE, LSURF, LNUM, L

      real*8 OX_SAT, P_FACTOR, DP_HALF
      real*8 OX_FLUX

      ! Bubble formation: Increase in required O2 pressure per meter:
      real*8 DEPTH_FACTOR / 0.48 /

C--------------------------------------------------------------------

C Calculate for each basin, only surface layer:
      DO IB = 1,NBI
C           WRITE(*,*) 'OXEXCH, BASIN ',IB
          LBASE = INDXI(IB)
          LSURF = LBASE+1
          LNUM = INDXI(IB+1)-LBASE

C .............. two_way exchange through surface:

C     Solubility in ml/l at 760 mm hg pressure:
          OXSAT(IB) = OXYGEN_SATURATION ( S(LSURF), T(LSURF) )

C     Oxygen specific exchange, given by wind (m/day):
          OXEXCF = OXSFAC* ( 0.04 + 0.67*MAX(0.0D0,WNDSPD-3)
     &                            + 1.07*MAX(0.0D0,WNDSPD-13) )
     &             * exp( 0.029*(T(LSURF)-20) )
C           temperature dependence according to Stigebrandt (1991)

C ------- Change of oxygen content derivative:
          OX_FLUX = - OXEXCF*(OXCONC(LSURF)-OXSAT(IB))*AREA(LSURF)
C           liter/day =    m/day*(liter/m3)           *m2
          OXYGDV(LSURF) = OXYGDV(LSURF) + OX_FLUX/VLAYER(LSURF)
          OXYGMP(IB)  = OXYGMP(IB) + OX_FLUX
C            WRITE(*,*) 'SURFACE EXCHANGE: OX_FLUX/VLAYER ',
C     &                  OX_FLUX/VLAYER(LSURF)
C            WRITE(*,*) 'OXYGDV(Lsurf) ',OXYGDV(Lsurf)

C .............. one_way release by bubbles:
          L=1
          DP_HALF = (DEPTH_FACTOR/2.0)*DEPTH(2) ! Increase to middepth
          OX_SAT = OXSAT(IB) * ( 1.0+ DP_HALF ) ! in middle of layer

C            WRITE(*,*) 'BUBBLES: , possibly in ',LNUM,' layers'
          DO while ( OXCONC(L+LBASE) .gt. OX_SAT )

C                 it is assumed that oxygen super-saturation does not
C                 occur below sub-saturation
             OX_FLUX = - OXSFAC*OXBUBL*(OXCONC(L+LBASE)-OX_SAT)
             OXYGDV(L+LBASE) = OXYGDV(L+LBASE) + OX_FLUX
             OXYGMP(IB)  = OXYGMP(IB) + OX_FLUX*VLAYER(L+LBASE)
C                WRITE(*,*) 'LAYER ',L,' (from surface)',
C    &               'OX_SAT, OXCONC, OX_FLUX, OXYGDV:',
C    &               OX_SAT, OXCONC(L+LBASE), OX_FLUX, OXYGDV(L+LBASE)
             IF ( L.GE. LNUM) EXIT

             L = L+1 ! Find oxygen saturation for next layer:
             OX_SAT =   OXYGEN_SATURATION ( S(L), T(L) ) ! at 1 atm.
             P_FACTOR = 1.0 + DP_HALF
             DP_HALF = (DEPTH_FACTOR/2.0)*DEPTH(L+1)
             P_FACTOR = p_factor + DP_HALF
                    ! = 1+ Depth_factor*Middle_depth

C                WRITE(*,*) 'LAYER ',L,' (from surface)',
C    &               'P_FACTOR, OX_SAT at 1 atm:',
C    &                P_FACTOR, OX_SAT
!                Pause
             OX_SAT = OX_SAT * P_FACTOR  ! in middle of layer
          ENDDO
!                Pause
      END DO

      End Subroutine


C ==========================================================
C Heat and kinetic energy exchange between surface water and atmosphere
      SUBROUTINE ENEXCH( WNDSPD, AIRTMP, HUMGM3, CLOUDS,
     &                   CDFAC, CEFAC, MBI, NBI, NLI, INDXI, TEMP,
     &           UFRIC3, EVAP, QOUT )
      
      
      SAVE

C  In:
      real*8 WNDSPD           ! effective wind speed in m/s
      real*8 AIRTMP           ! air temperature in deg C
      real*8 HUMGM3           ! humidity in g/cm3
      real*8 CLOUDS           ! Cloud cover in 0-9 (8'ths, 9=fog)
      real*4 CDFAC            ! Factor for C_D, controls wind friction
      real*4 CEFAC            ! Factor for C_e, controls evaporation
C                             and heat loss. If =0, turns off back rad.
      INTEGER MBI, NBI      ! Max. and actual number of basins
      INTEGER NLI           ! Number of layers
      INTEGER INDXI(NBI+1)  ! layer index limits for internal basins.
      real*8 TEMP (NLI)       ! Water temperature
C  Out:
      real*8 UFRIC3 (NBI)     ! Frictional velocity,
      real*8 EVAP   (NBI)     ! Evaporation (m3/s/m2)
      real*8 QOUT   (MBI,5)   ! Heat loss at surface (.,2..5) components

C ------------------------ local variables --------------------------
      INTEGER IB, LSURF

      real*8 ZKELVIN /273.15/            ! KELVIN at 0 Celsius
      real*8 TKA, TKS

      real*8 HUMID_AIR_MB, HUMID_SURF_GM3
      real*8 RHOAIR
      real*8 RHOWAT /1020.0 / ! approx. density of seawater (kg/m3)

C Heat flux terms:
      real*8 Q_LW_DOWN, Q_LW_UP, Q_EVAP, Q_COND
      real*8 Q_CLOUD_MAX / 60.0/          ! W/m2

      real*8 STEFAN_BOLZMAN / 5.67e-8 /  ! W/m2/K4

      real*8 WF_CD, WF_CE, DT_EFF
      real*8 STAB    ! Indicates stability regime (+1,-1)
      INTEGER SX   ! Index to coefficient table (1,2):

C Influence of stability of air column described by formulas below,
C with coefficients:

C ...... for wind friction:
      real*8 A_CD     (2) / 0.313, 0.023 /
      real*8 ALPHA_CD (2) / 0.842, 5.673 /
      real*8 BETA_CD  (2) / 0.968, 2.634 /

C ...... for sensible heat exchange and evaporation:
      real*8 A_CE     (2) / 0.371, 2.855 /
      real*8 ALPHA_CE (2) / 0.807, 1.648 /
      real*8 BETA_CE  (2) / 0.922, 1.722 /

      real*8 C_D, C_e
      real*8 EVAP_G
      real*8 LT

      real*8 HEAT_CAPACITY_AIR /1.0e3/ ! J/kg/K.
C This is heat capacity under constant pressure (Cp):
C According to Chemical Engineers handbook at p=1 atm, T=0 degC:
C         (0.238 cal/g/deg = 0.995 J/g/deg)
C 'Thermodynamics', Sears(Addison Wesley) p. 248 (Table 12-3):
C     Nitrogen: Cp = 3.5*R, Oxygen: Cp = 3.52*R,
C     Universial gas constant R = 8.3149 *1000 J/kgmol/deg
C     --> Dry air, with 80 % N and 20 % oxygen:
C         Cp = 8.3149*1000*(0.8*3.5/(14.008*2)+0.2*3.52/(2*16))
C            = 1.014J/g/deg
C NOTE!   Handbook of physics lists wrong value 4.18J/g/deg !!


C ============ Derived from air specs., common to all basins:

C ............ Air temperature in KELVIN:
      TKA = AIRTMP + ZKELVIN

C ............ Air humidity in g/m3 and millibar:
      HUMID_AIR_MB = HUMGM3 *
     &                  ( 1.26 + MIN(0.0151*AIRTMP, 0.00456*AIRTMP) )
C                                first const if T<0, second if T>0

C ........... Air density (should depend on temp. and humidity?)
      RHOAIR = 1.3 !kg/m3

C ........... Longwave radiation from atmosphere  (W/m2):
      Q_LW_DOWN = (0.61+0.05*HUMID_AIR_MB**0.5)*(TKA**4)*STEFAN_BOLZMAN
     &            + Q_CLOUD_MAX*CLOUDS/8.0


C ......... wind speed factors used below:
      WF_CE = WNDSPD**4.3
      WF_CD = WNDSPD**8


C =================== exchange for each basin:
      DO IB = 1, NBI

         LSURF = INDXI(IB)+1

C ........... Surface water temperatures in degrees KELVIN:
         TKS = TEMP(LSURF) + ZKELVIN

C ........... Stability (effect of humidity gradient not included):
         IF (TKS .gt. TKA) THEN
            STAB = 1    ! for unstable conditions
            SX = 1
         ELSE
            STAB = -1   ! for stable conditions
            SX = 2
         ENDIF
         DT_EFF = MIN(5.0D0, ABS(TKS-TKA))

C ........... transfer of kinetic energy:
C ... coefficient for wind friction, exported to be used outside
C     this subroutine:
         C_D = (0.8 + 0.9*WF_CD/(1.0e8+WF_CD )) ! Neutral coefficient
     &           * (1.0E-3)              ! derived from Smith (1980)
     &           *( 1.0 + A_CD(SX)*(DT_EFF**ALPHA_CD(SX))
     &                    /MAX(1.5D0,WNDSPD)**BETA_CD(SX)   )**STAB
C                                     modified for temperature stability
C                                     according to Bunker (1976).
         C_D = C_D* CDFAC ! Scaled by external factor

C ... Wind mixing energy as frictional velocity:'
        UFRIC3(IB)  = ( SQRT( RHOAIR/RHOWAT*C_D ) * WNDSPD )**3

C ............ transfer of heat energy:
C ... coefficient for evaporation and sensible heat loss:
         C_e = (1.3 + 0.6*WF_CE/(12.2**4.3+WF_CE)) ! Neutral coefficient
     &         * (1.0E-3) * 0.8
     &         *( 1 + A_CE(SX)*(DT_EFF**ALPHA_CE(SX))
     &                /MAX(1.5D0,WNDSPD)**BETA_CE(SX) )**STAB
C                neutral value after Bunker 1976, reduced by 20%,
C                and modified for unstable conditions.

         C_e = C_e* CEFAC ! Scaled by external factor

C ............ Saturated humidity in g/m3 at surface temperature:
         HUMID_SURF_GM3   = 4.8 * (TEMP(LSURF)/109.7 + 1)**7.594

C ... evaporation pr. surface     (g/m2/s):
         EVAP_G  = C_e * (HUMID_SURF_GM3 - HUMGM3)*WNDSPD

C ... specific evaporation energy (J/g):
         LT = 2494-2.2*TKS

C ... Energy loss by evaporation (J/m2/s=W/m2):
         Q_EVAP = EVAP_G * LT

C ........... Conductive (sensible) heat transfer (W/m2):
         Q_COND = RHOAIR * HEAT_CAPACITY_AIR * C_e * WNDSPD* (TKS-TKA)

C ........... Longwave radiation from water (W/m2):
         if (CEFAC.gt.0.0) then
            Q_LW_UP = 0.97 * STEFAN_BOLZMAN * TKS**4
         ELSE
            Q_LW_UP = 0.0  ! Can turn off back rad. for test purposes.
         ENDIF

C ----------------- Store output values ------------------------
C .......... Net heat loss from surface (W/m2):
         QOUT(IB,1) = Q_LW_UP - Q_LW_DOWN + Q_EVAP + Q_COND
C            store components for test purposes:
         QOUT(IB,2) = Q_LW_UP
         QOUT(IB,3) = - Q_LW_DOWN
         QOUT(IB,4) = Q_EVAP
         QOUT(IB,5) = Q_COND

C .......... Evaporation in m3/s/m2:
         EVAP(IB) = EVAP_G/1.0e6   ! Divide by g/m3 for fresh water
      ENDDO

      End Subroutine


C ==============================================================

C Radiative heat transport as light
C within water column.

    
      SUBROUTINE RADABS ( NBI, INDXI, MLI, WNDSPD, ND, DEPTH,
     &            QOUT, SUNRAD, DIFRAD, SINHS, IRFRAC, ICEFAC,
     &            RADFAC, PARTC, ATTNCF,
     &            AREA, Vlayer, TEMP, SAL, CEFAC,
     &            TEMPDV, HEATMP, RAD, QABS )
      
      

      SAVE

C IN:
      INTEGER NBI, MLI, ND, INDXI(NBI+1)
      real*8 WNDSPD, DEPTH(ND)
      real*8 QOUT(NBI), SUNRAD, DIFRAD, SINHS
	  real*4 IRFRAC, ICEFAC
      real*4 RADFAC(3)
      real*8 PARTC (MLI)
      real*4 ATTNCF(2)  ! ATTENUATION COEFFICIENTS FOR LIGHT:
                      ! (1) : BASE VALUE (NO ORGANIC PARTICLES)
                      ! (2) : INCREASE PR. mg/l Carbon

      real*8 Vlayer(MLI), AREA(MLI)
      real*8 TEMP(MLI), SAL(MLI)
      real*4 CEFAC
C UPDATED:
      real*8 TEMPDV (MLI)
      real*8 HEATMP (NBI)  ! Heat import to basin, unit [oC*m3/day]

C OUT:
      real*8 RAD (MLI), QABS(2)

C   All radiation terms should be in W/m2

      INTEGER IB, L, LD, LSURF
      real*8 IRF, S_RAD_in, D_RAD_in, S_RAD, D_RAD, albedo, HEAT_IN
      real*8 RAD_BELOW, RAD_ABOVE , RED_FAC, S_FAC, COSHS, SINHS_R
      real*8 EPSILON4 /1.e-6/
      real*8 AREA_BELOW
      real*8 ABSCFF

            ! BURDE V'RE BRUKERSTYRTE KOEFFISIENTER?



c       real*8 PI
c       PARAMETER ( PI=3.141593 )

C Specific heat of water; approx.  4.2 (joule/g/degC = Ws/cm3/degC)
C                         unit transformed to      (W*day/m3/degC)'
      real*8 HEAT_CAP
      parameter ( HEAT_CAP = 4.2 * 1.0E6 / ( 24 * 3600 ) )
C                                 cm3/m3    s/day

$if defined test_radabs
      write(999,'('' RADABS entered:'')')
$endif


C ------ cosinus to sun height = sinus to zenith angle:
      COSHS = SQRT(1.0-SINHS*SINHS)            ! INCIDENT ANGLE
      SINHS_R = SQRT(1.0-(COSHS/1.33)**2)    ! AFTER REFRACTION
C                        (----- SHOULD TAKE WIND INTO ACCOUNT HERE)

C ------ Path pr. length of vertical descent for direct solar radiation:
      S_FAC = 1.0/(SINHS_R+EPSILON4)  


C ------ Effective reflection at surface:
      if (sinhs.lt. 0.0 ) then
          albedo = 1.0 ! to avoid problems in expr. below, 
                       ! radiation=0 anyway
      else
          albedo = 1/(1 + 0.48*WNDSPD**0.62 + 47.0*SINHS**1.95)
      endif
C            SINHS = COS(zenith angle)
C       cfr. Spreadsheet \oslofj\klima\albedo,
C       based on PREISENDORFER & Mobley 1986 J. of Phys. Ocean.,16, 1293
C            ( Reflection at surface (Parsons et al. 1977, p.74) )


C ------ Total radiation entering the water column, after reflection:
      S_RAD_IN = SUNRAD*(1-albedo)
      D_RAD_IN = 0.95*DIFRAD
             ! Will be corrected for use below


C ------- Part of radiation assumed to be IR-radiation which is absorbed
C         in surface layer, and does not contribute to photosynthesis
      IRF = MAX(0.0, MIN(1.0, IRFRAC) )
      QABS(1) = IRF*(S_RAD_IN+D_RAD_IN)*RADFAC(3)
          ! QABS used below, RADFAC(3) controls heat effect

C ------- visual light penetrates water column:
       S_RAD_IN = (1-IRF)*S_RAD_IN*RADFAC(1)
       D_RAD_IN = (1-IRF)*D_RAD_IN*RADFAC(2)
          ! RADFAC(1) (2) controls light and effect, see below.

C         absorbed visual radiation energy:
       QABS(2) = S_RAD_IN + D_RAD_IN
C         ! Only for info, effect on temperature included below.

C ----- For each basin:
      DO IB = 1,NBI
         LSURF = INDXI(IB)+1

C     ------ Radiation entering the water column, after reflection:
         S_RAD = S_RAD_IN
         D_RAD = D_RAD_IN

C     ------ Update heat balance of surface layer with exchange terms:

C <<<<<<<<<<<<<<<< preliminary treatment of lower limit to temperature
       ! No further heat loss if temperature reaches freezing point.
       ! (lowered by salinity, with coefficient 0.055 according to
       ! table in HORNE 1966, Marine Chemistry (Wiley)
       IF (TEMP(LSURF) .GT. 0.0-0.055*SAL(LSURF) ) THEN
         HEAT_IN = (QABS(1) - QOUT(IB))*AREA(LSURF) /HEAT_CAP
         HEATMP(IB) = HEATMP(IB) + HEAT_IN
         TEMPDV(LSURF) = TEMPDV(LSURF) + HEAT_IN/Vlayer(LSURF)
       ELSE  ! Heat insulation:
         QOUT(IB) = 0.0
         S_RAD = ICEFAC*S_RAD ! Possible ice,
         D_RAD = ICEFAC*D_RAD ! reduce light by factor
       ENDIF

C     ------ Radiation integrated over horizontal area:
         RAD_BELOW = (S_RAD + D_RAD)*AREA(LSURF)

C     ------ Light intensity at top of each layer into RAD(L):
         DO L = LSURF, INDXI(IB+1)
             RAD(L) = S_RAD*S_FAC+D_RAD   ! Photosynthetically eff. rad.
             LD = L -INDXI(IB)
C          ...... sum of organic carbon, particulate and dissolved:
             ABSCFF = ATTNCF(1) + ATTNCF(2)*MAX(0.0D0,PARTC(L))
             RED_FAC = -ABSCFF*(DEPTH(LD+1) - DEPTH(LD))
             if (RED_FAC.GE. -100) THEN
                 RED_FAC = exp( MIN(0.0D0,RED_FAC))
             else
                 red_fac = 0.0
             ENDIF

             S_RAD = S_RAD*RED_FAC**S_FAC  ! at angle
C     Strictly: should transfer from s_rad to D_rad with increasing depth,
C     to take account of scattering.

             D_RAD = D_RAD*RED_FAC         ! assumed straight down
             RAD_ABOVE = RAD_BELOW   ! Rad. pr horiz. area
             IF (L.LT. INDXI(IB+1)) THEN
                AREA_BELOW = AREA(L+1)
             ELSE
                AREA_BELOW = 0.0
             ENDIF
             RAD_BELOW = (S_RAD + D_RAD)*AREA_BELOW

C Absorbed radiation heats layer if not turned off:
             HEAT_IN = (RAD_ABOVE-RAD_BELOW) /HEAT_CAP
             if (CEFAC.gt.0.0) then
                HEATMP(IB) = HEATMP(IB) + HEAT_IN
                TEMPDV(L) = TEMPDV(L) + HEAT_IN/Vlayer(L)
C                  unit degC/day
             endif

         END DO

C  ------- for phytoplankton growth control:
C          RAD(L) is converted to geometric mean intensity in layer:
         DO L = INDXI(IB)+1, INDXI(IB+1)-1   ! From surface down
            RAD(L) = SQRT(RAD(L)*RAD(L+1))   ! to next last layer.
         END DO
                                             ! Last layer:
         L = INDXI(IB+1)
         RAD(L) =SQRT( RAD(L) * (S_RAD*S_FAC+D_RAD) )
      END DO

$if defined test_radabs
      write(999,'('' exiting from RADABS'')')
$endif

      End Subroutine

      end Module
