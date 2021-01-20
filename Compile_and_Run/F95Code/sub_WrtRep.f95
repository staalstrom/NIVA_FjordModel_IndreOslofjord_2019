Module sub_WrtRep

      use ModelVar_RunControl
      use ModelVar_HydroBioChem
      use ModelVar_Topography
      
      implicit none

      character*78 Values
      integer IOCHECK
      real*8, pointer :: BALANCE (:)
      integer :: ReportUnit

	contains
	
      Subroutine WRTREP(FileUnit)
      integer :: FileUnit
      

      real*8 D_T, X, CV0, CV1, CV2
      PARAMETER (CV0 = 1.0E-9)

      Integer I_B, I_L, I_D

      Character*4 t_y/ 't/yr' /, tons/'Tons'/

      INTEGER ALLOC_CHECK, ALLOC_SPACE /0/


! ------------ Calculate and print info:

      ReportUnit = FileUnit

      if ( alloc_space.lt. NBI ) then
         if (alloc_space.gt.0) DEALLOCATE(BALANCE, STAT=ALLOC_CHECK)
         ALLOCATE  ( BALANCE(NBI), STAT= ALLOC_CHECK )
         if (ALLOC_CHECK.eq.0) THEN
             ALLOC_SPACE = NBI
         ELSE
             ALLOC_SPACE = 0
         ENDIF
      ENDIF
      write(*,*)' TINTEG:', TINTEG, ' TINTGZ:', TINTGZ
      D_T = (TINTEG - TINTGZ)/365.  ! Integrating time-step in years
      CALL PRT_NAME_AND_V ('End Time   (year)', TINTEG/365.)
      CALL PRT_NAME_AND_V ('StartTime  (year)', TINTGZ/365.)
      CALL PRT_NAME_AND_V ('T interval (year)', D_T)

      if ( D_T .le. 0.0 ) THEN

          CALL PRT_LIN(' Integration time <= 0')
          RETURN

      Endif


      CV1 = CV0/D_T ! From mg to metric ton/year


      !===============================
      CALL PRT_LIN (' C input:')
      CALL REPORT_line ( NBI, CV1, CLOADI,    ' From land    ',T_Y)
      CALL REPORT_line ( NBI, CV1, CSEDXI,    '- Perman. sed.',T_Y)

      !===============================
      CALL PRT_LIN (' Ox.demand input:')
      DO I_B = 1, NBI
         Balance (I_B) =   Odmloadi(I_B)*1.429*1.e+3  ! from l to mg
      Enddo
      CALL REPORT_line ( NBI, CV1, Balance,    ' From land    ',T_Y)

      !===============================
      CALL PRT_LIN (' N budget:')

      CALL REPORT_line ( NBI, CV1, NLOADI,    '+ From land   ',T_Y)
      CALL REPORT_line ( NBI, CV1, NSEDXI,    '- Perman. sed.',T_Y)
      CALL REPORT_line ( NBI, CV1, DNITRI,    '- Denitrified ',T_Y)
      CALL REPORT_line ( NBI, CV1, NFIXI ,    '+ N_fixation  ',T_Y)
      IF ( ALLOC_CHECK .EQ. 0 ) THEN
         DO I_B = 1, NBI
            Balance (I_B) =   NITRMI(I_B,1) - NITRMI(I_B,2)
         Enddo
         CALL REPORT_line ( NBI, CV1, Balance,'- Increase    ',T_Y)

      DO I_B = 1, NBI
            Balance (I_B) =  NLOADI(I_B) - NSEDXI(I_B) - Balance(I_B) &
                          - DNITRI(I_B) + NFIXI(I_B)
         Enddo
         CALL REPORT_line ( NBI, CV1, Balance,'= Exported    ',T_Y)
      ENDIF

      CALL REPORT_line ( NBI, CV0, NITRMI ,'= Within model',Tons )


      !===============================
      CALL PRT_LIN (' P budget:')

      CALL REPORT_line ( NBI, CV1, PLOADI,    '+ From land   ',T_Y)
      CALL REPORT_line ( NBI, CV1, PSEDXI,    '- Perman. sed.',T_Y)

      IF ( ALLOC_CHECK .EQ. 0 ) THEN
         DO I_B = 1, NBI
            Balance (I_B) =   PHOSMI(I_B,1) - PHOSMI(I_B,2)
         Enddo
         CALL REPORT_line ( NBI, CV1, Balance,'- Increase    ',T_Y)
         DO I_B = 1, NBI
            Balance (I_B) =   PLOADI (I_B) - PSEDXI(I_B) - Balance(I_B)
         Enddo
         CALL REPORT_line ( NBI, CV1, Balance,'= Exported    ',T_Y)
      ENDIF
      CALL REPORT_line    ( NBI, CV0, PHOSMI ,'= Within model',Tons)

      IF (ALLOC_CHECK .NE. 0) THEN
         WRITE(*,*) ' ALLOCATION ERROR IN WRTREP.FOR'
      ENDIF


      CALL PRT_LIN (' Sediment flux (C,N,P,S) in tons/year:')

      DO I_B = 1, NBI
         CALL PRT_NAME_AND_N ( ' Basin :', I_B)
         I_D = 2
         DO I_L = INDXI(I_B)+2, INDXI(I_B+1)
            X = AREA(I_L)*CV1  ! from mg/m2 to ton/year
            CALL REPORT_FLUX ( I_L, X, DEPTH(I_D), &
                               CDFLXI, NDFLXI, PDFLXI, SDFLXI )
            I_D = I_D +1
         Enddo

      Enddo

      CALL PRT_LIN (' Mean sediment flux (C,N,P,S) in g/m2/year:')
      CV2 = 1.e-3/D_T  ! from mg/m2 to g/m2/year

      DO I_B = 1, NBI
         CALL PRT_NAME_AND_N ( ' Basin :', I_B)

         I_D = 2
         DO I_L = INDXI(I_B)+2, INDXI(I_B+1)
            CALL REPORT_FLUX ( I_L, CV2, DEPTH(I_D), &
                               CDFLXI, NDFLXI, PDFLXI, SDFLXI )
            I_D = I_D +1
         Enddo

      Enddo


      END Subroutine


! ==================================================================



      SUBROUTINE REPORT_line ( N, X, V, V_Name,V_Unit )

      integer N
      real*8 X, V(N)
      character*(*) V_Name,V_unit

      integer L_Name, L_Unit, V_ptr, p_EndHead
      
      integer LField
      parameter (LField=18) 

      integer K, KK, I, V_pr_line, IOCHECK
      real*8 X_V, SUM

      L_name = Len (V_Name)
      Values = V_name
      L_unit = Len (V_Unit)
      p_EndHead = L_Name + 1+ L_Unit
      Values( L_Name+2 : p_EndHead ) = V_Unit
      V_pr_line = ( 78 - p_EndHead ) / LField
      Sum = 0.0
      DO K = 0, N, V_pr_line
         KK= MIN ( V_pr_line, N-K )
         V_ptr = p_EndHead + 1
         DO I = K+1, K+KK
            X_V = X*V(I)
            if (V_ptr.gt.78-LField) then
               write(*,*)'WRTREP error:V_Ptr=',V_ptr
            else
               WRITE ( Values(V_ptr:V_ptr+LField), "(1x,I3,':',G13.6)",iostat=IOCHECK) I, X_V
            Endif
            V_ptr = V_ptr + LField
            Sum = Sum + X_V
         Enddo
         if (KK.lt.V_pr_Line) then ! must happen last time (K from 0)
            WRITE ( Values(V_ptr:78), "(1x,'SUM:',G13.6)", iostat=IOCHECK) SUM
         endif
         IF (IOCHECK.NE.0) WRITE(*,*)'REPORT_LINE WRITE ERROR:',IOCHECK
         CALL PRT_LIN ( Values)
         Values(1:L_name) = ' '
      ENDDO

      end Subroutine

!=========================================================================
      Subroutine PRT_NAME_AND_V ( V_Name, X )
      Character*(*) V_Name
      real*8 X
      integer L_name
      L_name = Len_trim(V_Name)
      Values = V_name(1:L_Name)
      write ( Values(L_name+1:L_name+12),'(G12.5)', IOSTAT=IOCHECK) X
      IF (IOCHECK.NE.0) WRITE(*,*)'PRT_NAME_AND_V WRITE ERROR:',IOCHECK
      CALL PRT_LIN ( Values )
      end Subroutine

!=========================================================================
      Subroutine PRT_NAME_AND_N ( V_Name, N )
      Character*(*) V_Name
      Integer N
      integer L_name
      L_name = Len_trim(V_Name)
      Values = V_name(1:L_Name)
      write ( Values(L_name+1:L_name+5),'(I5)',IOSTAT=IOCHECK ) N
      IF (IOCHECK.NE.0) WRITE(*,*)'PRT_NAME_AND_N WRITE ERROR:',IOCHECK
      CALL PRT_LIN ( Values )
      END subroutine


! ==================================================================

      SUBROUTINE REPORT_FLUX ( I_L, X, DEPTH, CDFLXI, NDFLXI, PDFLXI, SDFLXI )
      INTEGER I_L
      real*8 X, DEPTH, CDFLXI(I_L), NDFLXI(I_L), PDFLXI(I_L), SDFLXI(I_L)

      Values = ' '
      WRITE( Values, '(" Layer:",I3,", depth:",F6.1,"m:",G13.6,3G12.5)', &
            IOSTAT=IOCHECK ) I_L, DEPTH, &
            CDFLXI(I_L)*X, NDFLXI(I_L)*X, PDFLXI(I_L)*X, SDFLXI(I_L)*X
      IF (IOCHECK.NE.0) WRITE(*,*) 'WRITE ERROR:',IOCHECK,'in subroutine REPORT_FLUX'
      CALL PRT_LIN( Values )
      END subroutine

! ==================================================================

      Subroutine PRT_LIN( String)
      Character String*(*)

      WRITE ( ReportUnit  ,'(1x,A)') String

      End subroutine





!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$IFDEF  (xxxxxx) ! Unfinished code:
!C  Subroutine performed when ACSL command WRTREP is entered
!
!      INCLUDE eutro.inc
!
!      WRITE(*,'(A)' )' Chose operations:',
!     &   ' > Filename: Store state description on file'
!     &   ' <2. Retrieve state description from file'
!      read (
!
!
!      open (STORAGE_FILE, FILE=FILE_NAME )
!C  Dimensions:
!      WRITE (STORAGE_FILE)
!     &   T, NBI, NLI, NC, NLC, NBE, NLE, NC, NLC
!
!C  Internal state variables:
!      WRITE (STORAGE_FILE)
!     &   ( SAL(I),  TEMP(I),  OXYG(I),
!     &     DOC(I),  BACT(I),  CZOO(I),
!     &     PO4(I),  NO3(I), NH4(I),
!     &     SIO2(I), SFYT(I), SSED(I), I=1,NLI ),
!     &   ( ( CFYT(I,K), NFYT(I,K), PFYT(I,K), CHL(I,K), AGEF(I,K)
!     &        CSED(I,K), NSED(I,K), I=1,NLI ), K=1,2 )
!     &   ( ( XSED(I,K), I=1,NLI ), K=1,2 )
!
!
!      END

!$ENDIF

end Module