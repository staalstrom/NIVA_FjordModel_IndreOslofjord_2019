Module Binary_Write_sub

! ======================================================================
! NIVA Fjord Model, file BinWrite.F95 created in Absoft v11

! Writes data from a model run to a binary file for later transfer to other file formats.
! The file is opened unformatted with sequential access,
! and with the option FORM='BINARY', which is an ABSOFT extension,
! giving an unformatted file without record length info embedded.
! It corresponds to the WATCOM open option RECORDTYPE='FIXED',
! which when used for an unformatted file with sequential access
! gives "Unformatted, sequential, binary data with a fixed record type".
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note that on Windows computers, the file will have Little-Endian format,
! which means that for numbers, the least significant bytes (or bits) are stored first.
! The Open option CONVERT= can be used to change this.
! Programs reading the binary files must read with the same order, or use conversion.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ======================================================================

      use BinaryFile_ControlVariables

      implicit none
      
      character*8, allocatable :: vNameArray (:) 
              ! 8 is max. name length in old Statistica format currently used

      real*8     , allocatable :: vArray     (:)
              ! All values stored as FORTRAN real*8 (=double floating values in C++).

      integer*4 :: vCount       ! used to count variables in subroutine Next_Var
      integer*4 :: NumVarArrays ! Stores max. count, used to allocate arrays

     
      logical   :: Binary_File_Is_Open  = .false.
      logical   :: Arrays_Are_Allocated = .false.

      integer   :: LastOperation = -1
                  ! notes OpCode between calls to Next_Var (see below)
                  ! used to check that call sequence to this module is valid


      logical, parameter:: DEBUG_BIN_RES=.false. ! Controls debug print


   contains 


!=======================================================================
! Open file for binary output and write number of values per record:

      logical function BinaryWriteOpenFile(UnitNumber, FileName)

      integer       UnitNumber   ! File unit to use in write operations 
      character*(*) FileName     ! File name to write (or overwrite)

      !-------------------------------------------------

      
      BinaryUnit = UnitNumber  ! to be used in calls below

      BinaryFormatVersion = 1

      open(UNIT = BinaryUnit, FILE = FileName, STATUS = 'UNKNOWN', &
        ACCESS  = 'SEQUENTIAL', ACTION = 'WRITE', &
        FORM    = 'BINARY', & ! file written without embedded length information
        IOSTAT  = BinFileErrorStatus)

! File header for file format version 1:

      if (BinFileErrorStatus.eq.0) &
         write (BinaryUnit, IOSTAT = BinFileErrorStatus) BinaryFileHeader
      if (BinFileErrorStatus.eq.0) &
         write (BinaryUnit, IOSTAT = BinFileErrorStatus) BinaryFormatVersion

      Binary_File_Is_Open = (BinFileErrorStatus.eq.0)     
      BinaryWriteOpenFile = Binary_File_Is_Open
   end function



      ! =================================================================
      ! Reset counter. Should be called before sequence of calls to Next_Var

   subroutine Init_VarCount()
      vCount=0
      return
   end subroutine
      
      ! Cound move opCode to this call, with checks on valid sequence of operations,
      ! and save it in internal vaiable for use in Next_Var below.


      ! ----------------------------------------------------------

   subroutine Next_Var(opCode, vName, V)
         ! operate on one variable according to argument opCode.

      integer opCode 
         ! =0: count values 
         ! =1: store names in name array
         ! =2: store value in value array

      Character*(*) vName
      real*8 V

      vCount = vCount+1

      if (opCode.le.0) return

      if (.not.Arrays_Are_Allocated) return

      if (vCount.gt.NumVarArrays) then
         write(*,*)'Exceeds allocated array lengths in subroutine Next_Var' 
         return
      endif
  
      Select Case (Opcode)
         Case (1)
            vNameArray(vCount)=vName
         Case (2)
            if (vNameArray(vCount) .eq. vName) vArray(vCount) = V
      end Select
      if (DEBUG_BIN_RES) then
         write(*,'(1x,A,I6,3A,I6)') &
            ' Next_Var: opcode=',opCode, ' vName=',vName, ' vCount=',vCount
      end if
      return
   end subroutine

! ==================================================================

   subroutine AllocateArrays()

		      ! Allocates arrays according to counting done with calls to 
		      ! Init_VarCount and Next_Var. 
            ! Should only be done once for each opened Binary file.
            ! is required before Next_Var can be called with OpCode>0


      ! -----------------------------------------------------------

		      ! execute AllocateArrays: 
		      ! Only once for each program run.

      if (vCount.le.0) Stop ' no variables initiated for binary result storage' 

      if ( .not. Allocated(vNameArray)) then
         NumVarArrays = vCount 
         Allocate (vNameArray(NumVarArrays), vArray(NumVarArrays),STAT=BinFileErrorStatus)
         if(BinFileErrorStatus.eq.0) then
            Arrays_Are_Allocated = .true.
         else
           write(*,'(1x,3(A,I6))') &
          ' Allocation of arrays failed in Binary_Write_sub;', &
          ' BinFileErrorStatus=', BinFileErrorStatus,' vCount=',vCount
           Arrays_Are_Allocated = .false.
         endif
      endif
      
   end subroutine


!=======================================================================
! Write variable headings, preceded by size info:

   subroutine BinaryWriteNames ()
   
      integer*4 NameLength

      if (.not. Arrays_Are_Allocated) return

      NameLength = len(vNameArray(1))

      write (BinaryUnit, IOSTAT = BinFileErrorStatus) vCount

      if (BinFileErrorStatus.eq.0) &
         write(BinaryUnit, IOSTAT = BinFileErrorStatus) NameLength

      if (BinFileErrorStatus.eq.0) then
         write(BinaryUnit, IOSTAT = BinFileErrorStatus) vNameArray
      endif

      if(BinFileErrorStatus.ne.0) then
         write(*,*) 'Error ',BinaryWriteErrCode(), &
                    ' on writing names to binary file ',BinaryUnit
         STOP
      else
	      if (DEBUG_BIN_RES) then
	         write(*,*) 'names to binary file OK'
	      endif
      endif

   end subroutine

!=======================================================================
! Write array of real*8 values:

      subroutine BinaryWriteValues ()

      if (.not. Arrays_Are_Allocated) return
      write(BinaryUnit, IOSTAT = BinFileErrorStatus) vArray

      if (BinFileErrorStatus.ne.0) then     
         write(*,*) 'Error ',BinFileErrorStatus, &
                   ' on writing values to binary file ',BinaryUnit
         STOP
      else
         if (DEBUG_BIN_RES) then
            write(*,*) 'values to binary file OK'
         end if
      endif

      end subroutine
      
!=======================================================================
! Close binary file:

   logical function BinaryWriteClose()
      Endfile BinaryUnit
      close(BinaryUnit,IOSTAT = BinFileErrorStatus)
      BinaryWriteClose = (BinFileErrorStatus.eq.0)
      Binary_File_Is_Open  = .false.     
   end function
      
!=======================================================================
! Return error code:

   integer function BinaryWriteErrCode()
      BinaryWriteErrCode= BinFileErrorStatus
   end function
   
end module Binary_Write_sub
