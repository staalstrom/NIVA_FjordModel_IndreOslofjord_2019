! NIVA Fjord Model, file: fx_RunControl_DebugPrint
!                         fixed F77 format

! Revised by Birger Bjerkeng 2019-10-13:
!   Subroutines out_xxxx for writing to ASCII files changed:
!       Output files are only opened at start of run.
!       Formats now handle any number of layers, using ADVANCE='no'/'yes',
!       each time record is pre-fixed with "T=", and a space is inserted between each value.

      Module fx_RunControl
      
      use ModelParam_RunControl, only: DBGDEV, DEBUG_STEPS, TTRIG, TRACE, NPRINT, REINTG
      use ModelVar_RunControl, only: T, TINTEG, TSTEP, TDERIV, NDERIV
      use ModelVar_Topography, only: NBI, INDXI, BasinName
      use ModelVar_HydroBioChem
!      use ModelVar_Mussels, only: MUSLMT
      
      implicit None

      Integer :: DEBUG_UNIT =999
      LOGICAL isOPEN ! used in subroutines out_XXXXX below
      

      Contains


C ------- initiate step, make ready for using dump files
C         in subroutines:

      SUBROUTINE DEBUGF(MODE)

      INTEGER MODE ! Controls file handling at start and end of simulation

       
      LOGICAL OPEN_UNIT


      if ( DBGDEV.ne.6 .and. DBGDEV.ne.9 .and. DBGDEV .lt. 100 ) THEN
         write(*,*)'New debug unit DBGDEV changed from illegal value ',
     &              DBGDEV,' to 999'
         DBGDEV = 999 ! from illegal value
      ENDIF

      INQUIRE (UNIT=DEBUG_UNIT, OPENED=OPEN_UNIT)

      IF (MODE .le.0 .or. DBGDEV.ne.debug_unit ) THEN
C          Initial or final call, or unit number changed by user:
!          STEP_NUMBER = -1  ! New file for each simulation
          if ( OPEN_UNIT .and. DEBUG_UNIT.ne.6 .and. debug_unit.ne.9 )
     &       CLOSE(DEBUG_UNIT)
          DEBUG_UNIT = DBGDEV

      ELSEIF (OPEN_UNIT) THEN
          IF ( DEBUG_UNIT.NE.6 .AND. DEBUG_UNIT.NE. 9) THEN
             IF(DEBUG_STEPS.LE.0) THEN
                CLOSE(DEBUG_UNIT)
             ENDIF
!             IF( STEP_NUMBER.LE.0) THEN
!                WRITE(*,*) 'Num. of steps to dump on unit ',
!     &                      DEBUG_UNIT,'?'
!                READ(*,*,IOSTAT=errchk) STEP_NUMBER
!                WRITE(*,*) 'errchk, STEP_NUMBER=', errchk, STEP_NUMBER
!                IF ( ERRCHK .NE. 0 ) step_number = 3
!             ENDIF
!             STEP_NUMBER = MAX( 0, STEP_NUMBER -1 )
          ENDIF
      ENDIF
  
          ! count down debugging steps towards zero:
      if (DEBUG_STEPS.gt.0) DEBUG_STEPS = DEBUG_STEPS-1
              ! if zero: suppress any debug printing:
      if (DEBUG_STEPS.le.0) TTRIG=T+1000


      END SUBROUTINE


C ================================================================
C Subroutine called from ACSL model to trace execution.
C Can be controlled by logical variable TRACE at run_time.

      SUBROUTINE HELLO( SECTION_NAME)
      
      Character*(*) SECTION_NAME

!      INCLUDE 'EUTRO.INC'

      INTEGER PRINT_UNIT(3)/3*6/
      INTEGER UNITS, K
!      INCLUDE 'DEBUG.CMN'
      
      
      IF (TRACE) THEN
!!! special ACSL function:  CALL AGETI ( 'PRN', Print_UNIT(2) )
         PRINT_UNIT(3) = DEBUG_UNIT
         UNITS = 3
         if (PRINT_UNIT(3).eq.PRINT_UNIT(2)) UNITS = 2
         if (PRINT_UNIT(2).eq.PRINT_UNIT(1)) THEN
            DO K = 3, UNITS
               PRINT_UNIT(K-1)=PRINT_UNIT(K)
            END DO
         END IF

         DO K=1,UNITS
            WRITE( PRINT_UNIT(K), '(" ### Trace:  ",A,
     &                           "  at",3(1x,A,G13.7))')
     &               SECTION_NAME,
     &              'T:',T    !, 'TINTEG:',TINTEG, 'TDERIV:',TDERIV
         END DO
      ENDIF
      
      END Subroutine
      
      ! ===============================================================
! Controls model progress monitoring report
! to debug unit and standard output unit:

      SUBROUTINE TIMPRT(ProgressQualifier)
      character*(*) :: ProgressQualifier


      integer UNIT_now, UNIT_LAST
      integer I
      logical FILE_opened

      IF ( MOD(NDERIV-1, int(NPRINT,kind(NDERIV)) ) .NE. 0 ) RETURN

      inquire(DEBUG_UNIT, OPENED = file_opened )

      if( file_opened ) THEN
         UNIT_LAST = DEBUG_UNIT
      ELSE
         UNIT_LAST = 0
      ENDIF

      I = INDXI(2)
      UNIT_now = 6

         WRITE ( UNIT_now, '(3A,G14.6/4(4x,A,G14.6),A,I4,A)' )
     &          ' *** TIMPRT(',ProgressQualifier,'): T=',T,
     &          ', TSTEP=', TSTEP, ', OXYG(I)', OXYG(I),
     &          ', NO3(I)',NO3(I), ', NH4(I)', NH4(I),'(I=',I,')'
!         write(UNIT_now,*) 'MAXTTR, MXTBIO: ', MAXTTR, MXTBIO
         write( Unit_Now, *)' NDERIV: ',NDERIV, 'NPRINT: ',NPRINT
! 	      DO I = 1,NBI
!	         ZSURFI (I) = VDYN(1,I)/area(INDXI(I)+1)
!	      end do     
!         WRITE ( UNIT_now, * ) '       ZSURFI: ',(ZSURFI(I),I=1,NBI)
!         WRITE ( UNIT_now, * ) '       ZSURFE: ',(ZSURFE(I),I=1,NBE)
         UNIT_now = UNIT_LAST
         UNIT_LAST = 0

      END SUBROUTiNE

! ===============================================================
! ANS@NIVA.NO Subroutines for writing to ASCII files
! Called from Calc_Derivatives in module m3_RunModel
!

! Revisions by Birger Bjerkeng 2019-10-13:

!  Files are opened if REINTG=.true., which is the case after
!  subroutine StartRun has called InitiateModelState.
!  Files are left open for later calls during a run, when REINTG=.false.
!  The files are automatically closed when the program exits.

!  NOTE: Argument ProgressQualifier is not used. Is this intended?
! ===============================================================   
      
      SUBROUTINE out_MUSSEL(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 37

      call Check_OutFile(UNIT_now,'out_MUS_1.dat')

      WRITE ( UNIT_now, '("T= ",f14.4, 9(" ",F14.3))' )
     &          T, MUSLMT(1:9,1,1)
     
      END SUBROUTINE

! ===============================================================   

      SUBROUTINE out_OXYG(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 20
      
      call Check_OutFile(UNIT_now,'out_OXYG.dat')
      call out_LayerArray(UNIT_now,OXYG)     
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_SAL(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 21
      
      call Check_OutFile(UNIT_now,'out_SAL.dat')
      call out_LayerArray(UNIT_now,SAL)     
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_TEMP(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 22
      
      call Check_OutFile(UNIT_now,'out_TEMP.dat')
      
      call out_LayerArray(UNIT_now,TEMP)     
     
      END SUBROUTINE

! ===============================================================      
      
      SUBROUTINE out_TOTC(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 23
      
      call Check_OutFile(UNIT_now,'out_TOTC.dat')
      
      call out_LayerArray(UNIT_now,TOTC)     
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_DOC(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 24
      
      call Check_OutFile(UNIT_now,'out_DOC.dat')
      call out_LayerArray(UNIT_now,DOC)     

      END SUBROUTINE   
      
! ===============================================================      
      
      SUBROUTINE out_TOTN(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 25
      
      call Check_OutFile(UNIT_now,'out_TOTN.dat')
      call out_LayerArray(UNIT_now,TOTN)     
     
      END SUBROUTINE
    
! ===============================================================      
      
      SUBROUTINE out_NO3(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 26
      
      call Check_OutFile(UNIT_now,'out_NO3.dat')
      
      call out_LayerArray(UNIT_now,NO3)     
     
      END SUBROUTINE

! ===============================================================      
      
      SUBROUTINE out_NH4(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 27
      
      call Check_OutFile(UNIT_now,'out_NH4.dat')
      
      call out_LayerArray(UNIT_now,NH4)     
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_TOTP(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 28
      
      call Check_OutFile(UNIT_now,'out_TOTP.dat')
      
      call out_LayerArray(UNIT_now,TOTP)
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_PO4(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 29
      
      call Check_OutFile(UNIT_now,'out_PO4.dat')
      
      call out_LayerArray(UNIT_now,PO4)     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_CZOO(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 30
      
      call Check_OutFile(UNIT_now,'out_CZOO.dat')
      
      call out_LayerArray(UNIT_now,CZOO)
 
      END SUBROUTiNE
      
! ===============================================================      
      
      SUBROUTINE out_CFYT1(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 31
      
      call Check_OutFile(UNIT_now,'out_CFYT1.dat')
      call out_LayerArray(UNIT_now,CFYT(1,1))
      
      END SUBROUTINE

! ===============================================================      
      
      SUBROUTINE out_CFYT2(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 32
      
      call Check_OutFile(UNIT_now,'out_CFYT2.dat')
      
      call out_LayerArray(UNIT_now,CFYT(1,2))
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_DEEPW(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 33
      
      call Check_OutFile(UNIT_now,'out_DEEPW.dat')
      
      WRITE ( UNIT_now, '("T= ",f11.4, 2(" ",f10.3))' )
     &          T, SAL(39), OXYG(39)
     
      END SUBROUTINE
      
! ===============================================================      
      
      SUBROUTINE out_SURW(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 34
      
      call Check_OutFile(UNIT_now,'out_SURW.dat')
      
      WRITE ( UNIT_now, '("T= ",f11.4, 9(" ",f10.3))' )
     &          T, TOTN(1), NO3(1), NH4(1), TOTP(1), PO4(1), CFYT(1,1), CFYT(1,2), CHL(1,1), CHL(1,2)
     
      END SUBROUTINE      

! ===============================================================      
      
      SUBROUTINE out_CHL1(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 35
      
      call Check_OutFile(UNIT_now,'out_CHL1.dat')
      call out_LayerArray(UNIT_now,CHL(1,1))

      END SUBROUTINE

! ===============================================================      

      
      SUBROUTINE out_Si(ProgressQualifier)
      character*(*) :: ProgressQualifier

      integer UNIT_now
      
      UNIT_now = 36
      
      call Check_OutFile(UNIT_now,'out_Si.dat')
      
      WRITE ( UNIT_now, '("T= ",f11.4, 4(" ",f10.3))' )
     &          T, SiO2(1), SFYT(1), SDFLUX(1), SSED(1)

      END SUBROUTINE

! ===============================================================  
! Check status of output files before writing in subroutines out_Xxxxx above:
      Subroutine Check_OutFile(Unit, FileName)
      integer Unit
      character*(*) FileName

      IF(REINTG) then ! In Start of run: close if open from earlier run
                      ! then open as new file:
         INQUIRE(Unit, OPENED=isOPEN)
         IF(isOPEN) CLOSE(Unit) 
         OPEN (Unit, FILE=FileName, STATUS='UNKNOWN')
         call Hello(" Opened " // FileName)
      ENDIF

      End Subroutine

! ===============================================================  

      Subroutine out_LayerArray (Unit, Arr)

      integer Unit
      real*8 Arr(*)
      
      Character*4 AdvanceMode
      integer I, N

      WRITE ( Unit, '("T= ",f11.4)', ADVANCE='no') T
      AdvanceMode='No'
      N=19
      do I =1,INDXI(NBI+1),20
         if (I.gt.INDXI(NBI+1)-20) then
            AdvanceMode='Yes'
            N=INDXI(NBI+1)-I
         endif
         WRITE ( Unit, '(20(" ",F10.3))', ADVANCE=AdvanceMode ) Arr(I:I+N)
      enddo

      END SUBROUTINE
      
! ===============================================================      


      END Module
