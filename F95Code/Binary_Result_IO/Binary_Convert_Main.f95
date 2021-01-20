!============================================================================
! ABSOFT FORTRAN 95 main program.
!
! Uses module BinRead_Absoft.F95 to read a binary file
! created by NIVA_Fjord_Model)
! and sends the data to a specified fileformat:
!   In this version: to a Statistica data file
!                    by calls to dll library STA32.dll
!                       (STADEV functions defined for C++).
!   
!============================================================================

! Start with command-line arguments:
!       ProgramName BinaryFilename OutputFileName, RunId_ShortLabel, "Long Runid descr. text" RunidNumber
! Argument 1-3 are mandatory, 4 and 5 optional.         

 
      program Binary_Convert
      
      use Binary_Read
      use Binary_Convert_Sub_Library
      use Binary_Convert_Branching
      
      implicit none


             ! Controls test print statements:
      logical DEBUG /.false./ 
             ! set to true to turn on
             ! also activates some Pause statements


            ! File for printing error messages and debug output:
            ! (If 0: to Standard unit, no file is opened)


      !-------------------------------------------------------------------------
      !declarations:
      !-------------------------------------------------------------------------

          ! F77 Functions in BinRead.for:

!	      logical, external :: BinaryReadOpenFile, BinaryReadNames, BinaryReadValues
!	      integer, external :: BinaryReadErrCode
      





         ! File specifications (program command line arguments):
         ! =====================================================

      integer narg
      integer, external :: iargc 

      integer        BinaryUnitNumber/7/
      character*1028 BinaryFileName, OutputFileName
      integer*4 BinaryReadError
      
      integer ErrStatus
      
      integer iascii


          ! Buffer for Run identification strings
          ! (one space extra for null Char):
      
      character*9  RunId_ShortLabel/' '/
      character*41 RunId_LongLabel/' '/
      character*12 RunId_Value_String
      integer*4 RunId_Specified_Value

      integer IOStatus        




				! Arrays for transferring info from Binary file:
          	! They are allocated using vCount, 
			 	! and vNameArray string length is checked against NameLength
			 	! with both vCount and Namelength declared in Binary_Read
				! and read from the Binary file.

      character*8 , allocatable :: vNameArray (:)
      real*8      , allocatable :: vArray     (:,:)


      integer RecordBatch/64/

!      integer*2 last
      
      logical EOF_Found
      

      real*8 FileTime 



!-------------------------------------------------------------------------
!      executable code 
!-------------------------------------------------------------------------


!========================================================================================
! Open file for writing diagnostics:

      if (DiagUnit .ne.0) then 
         open(DiagUnit,File=DiagFile, IOSTAT=IOStatus)   
	      if (IOStatus .ne.0) then
	         if (DEBUG) Pause 'Error on opening Diagnostics file, press Enter to close'
	         Stop
         endif
      endif


!========================================================================================
! Get program arguments:

      Narg = iargc()
      RunId_Specified_Value =-9999 ! default if not specified as 5th program argument.
                                   ! (for adding to existing Statistica files: 
                                   ! if specified value is already in use, 
                                   ! the largest used value+1 is used instead.

      if (Narg.lt.3) then
         Stop 'Requires 3 command line arguments: BinaryFileName, OutputFileName and RunId_ShortLabel'
      endif   
      call getarg(1, BinaryFilename)
      write(DiagUnit,*)'BinaryFilename:',trim(BinaryFilename)
      call getarg(2, OutputFileName)
      write(DiagUnit,*)'OutputFileName:',trim(OutputFileName)
      call getarg(3, RunId_ShortLabel)
      write(DiagUnit,*)'RunId_ShortLabel:',trim(RunId_ShortLabel)

      if (DEBUG) write(DiagUnit,*) ' Narg:', Narg

      if (Narg.le.3) then
         RunId_LongLabel=''
      else
         call getarg(4,RunId_LongLabel)
         write(DiagUnit,*)'RunId_LongLabel:',trim(RunId_LongLabel)

         if (Narg.gt.4) then
            call getarg(5,RunId_Value_String)
            write(DiagUnit,*) 'RunId_Value_String:', trim(RunId_Value_String)

            read ( RunId_Value_String, *,iostat = IOStatus) RunId_Specified_Value 
            if (IOStatus.ne.0) then
                write(DiagUnit,*)'Error in command line argument Runid_Value'
            else
               write(DiagUnit,*) ' RunId_Specified_Value:', RunId_Specified_Value
            endif
         endif
      endif


      if (DEBUG) then
         write(DiagUnit,'(1x,2A)') &
                 ' BinaryFilename:'       , trim(BinaryFilename), &
                 ' StatisticaFileName:'   , trim(OutputFileName), &
                 ' RunId_ShortLabel:'     , trim(RunId_ShortLabel), &
                 ' RunId_LongLabel:'      , trim(RunId_LongLabel)
         write(DiagUnit,*) ' RunId_Specified_Value:', RunId_Specified_Value
         pause "Press Enter to continue..."
      endif  


! ===========================================================================
! Open Binary file name for read, and read variable count from file into variable vCount:

      if (.not.BinaryReadOpenFile(BinaryUnitNumber, trim(BinaryFileName))) then
         write(DiagUnit,'(1x,2A)') ' Error on opening binary file ', &
                              BinaryFileName
         goto 998
      endif

      call t_Days_FileUnit (BinaryUnitNumber, FileTime)

      write(DiagUnit,*) ' FileTime :', FileTime 

!========================================================================================
! Set up work arrays for data from binary file:

      if (DEBUG) write(DiagUnit,'(1x,3(A,I6))') ' vCount=',vCount

      Allocate(vNameArray(-1:vCount), vArray(-1:vCount, RecordBatch), STAT=ErrStatus)
      if(ErrStatus.eq.0) then
         write(DiagUnit,'(1x,3(A,I6))') ' Allocation OK'
      else
         write(DiagUnit,'(1x,3(A,I6))') ' Allocation Error ', ErrStatus
         goto 999
      endif

 
!========================================================================================
! Saving Run_ID and FileTime in array elements -1 and 0:
      vNameArray(-1) = 'RunIdent'
      vNameArray (0) = 'FileTime'
     

! add variable names from binary file:

      if (.not.BinaryReadNames (vNameArray(1))) then
         write(DiagUnit,*)' reading of variable names from binary file failed'
         goto 998
      endif

      write(DiagUnit,*)' reading of variable names from binary file succeeded'



! Prepare Output file for storing records:

      		
		
      if (.not. OpenOutputFile(OutputFileName, vNameArray, &
					 Runid_Specified_Value, RunId_ShortLabel, RunId_LongLabel)) goto 999

  ! Added by ans@niva.no - write to ascii file 
      OPEN (49, FILE='eutro.csv', STATUS='new',POSITION="append")
      CLOSE(49)
  ! ans@niva.no - finished 



!========================================================================================
! Read records from binary file,
! and store in new records in Statistica file:


      EOF_Found = .false.

      
      do while (.not.EOF_Found)

            ! ========== accumulate records in program buffer:

         RecordNum = 0
         
         CollectRecords: do while (RecordNum .lt. RecordBatch)
            if(.not. BinaryReadValues (vCount, vArray(1:vCount,RecordNum+1), BinaryReadError)) then
	            EOF_Found=(BinaryReadError.lt.0)
               if (EOF_Found) then
                  exit CollectRecords
               else
                  goto 998
               endif
            endif
            RecordNum = RecordNum + 1
            vArray(-1,RecordNum) = Runid_Value
      		vArray (0,RecordNum) = FileTime

            
            
  !          if (DEBUG) write(DiagUnit,*)' Found record number ', RecordNum,' in Binary file'
         end do CollectRecords
  ! Added by ans@niva.no - write to ascii file       
         write(DiagUnit,*) ' RecordNum ', RecordNum
         CALL write_ascii(vArray,RecordNum)
  ! ans@niva.no - finished 

			if (.not. AddRecordsToOutputfile(RecordNum, vArray)) goto 999

      end do
       

         


 998  if (BinaryReadError.gt.0) then
         write(DiagUnit,*)"Error ", BinaryReadError," in operations on Binary file"
      endif       

        ! Close files and terminate program:

 999  continue
      if (.not. Closeoutputfile()) then
			write(*,*)' Could not close outputfile "', trim(OutputFilename),'"'
      endif

      end program
