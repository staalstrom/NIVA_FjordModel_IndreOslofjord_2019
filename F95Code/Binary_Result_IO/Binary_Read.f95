! ======================================================================
! ABSOFT FORTRAN Code for reading from binary file written by Eutro model
! when transferring results to other file format (e.g. Statistica)

! The file is opened for unformatted sequential access and binary format,
! which in ABSOFT Fortran specifies no record delimiters or length information
! ======================================================================

 Module Binary_Read
      
      use BinaryFile_ControlVariables
 
      implicit none 


      character*128 FileHeaderBuffer
      
      real*4 :: MissingValueCode = -1.0e+39
      
      
!      character*8 , allocatable :: vNameArray (:)
!      real*8      , allocatable :: vArray     (:,:)

      integer*4 :: vCount     ! nuhmber of variables in Binary file
      integer*4 :: NameLength ! Length of variable names in binary file
      
      logical   :: BinaryFile_Open =.false.

  contains

! ----------------------------------------------------------------------

      logical function BinaryReadOpenFile(UnitNumber, FileName)
      integer, intent(in)      :: UnitNumber
      character*(*), intent(in):: FileName
      
      integer*4 HeaderLength

      BinaryReadOpenFile = .false. ! reset on success
      BinaryUnit = UnitNumber
      
      open(UNIT = BinaryUnit, FILE = FileName, STATUS = 'OLD', &
        ACCESS  = 'SEQUENTIAL', ACTION = 'READ', &
        FORM    = 'BINARY', & ! file without embedded length information
        IOSTAT  = BinFileErrorStatus)

      if (BinFileErrorStatus.ne.0) return
     
! Check file version, attempting to read file header:

      HeaderLength = len(BinaryFileHeader)
      if(HeaderLength.gt. len(FileHeaderBuffer)) then
         write(*,*) ' too small character array FileHeaderBuffer in module Binary_Read'
         return
      endif

      BinaryFormatVersion = -1     ! undefined so far

            ! Check if file has the defined File Header text,
            ! if so, read format version number:

      read (BinaryUnit, IOSTAT = BinFileErrorStatus) FileHeaderBuffer(1:HeaderLength)

      if (BinFileErrorStatus.eq.0) then
         if(FileHeaderBuffer(1:HeaderLength).eq.BinaryFileHeader) then
            read (BinaryUnit, IOSTAT = BinFileErrorStatus) BinaryFormatVersion
      		if (BinFileErrorStatus.eq.0) then
            	   if (BinaryFormatVersion .ne. 1) then
               	      write(*,*)' Invalid format version:',BinaryFormatVersion, &
                  	       ' in Binary file header'
               	      return
                   endif
                endif
         endif

				! did not find header, failed on reading version, or version=1:
            ! Otherwise; could be old file format from WATCOM version model,
            ! assigned as version 0. Should then have BASSENG and Lag as two first variables,
				! will be checked when reading name array.
            
         if (BinaryFormatVersion .eq. -1) then ! can be version 0 without header
            REWIND(BinaryUnit,IOSTAT = BinFileErrorStatus)
            BinaryFormatVersion = 0
         endif
      endif

      if (BinFileErrorStatus.eq.0) then
         read(BinaryUnit,IOSTAT = BinFileErrorStatus) vCount
         if (BinFileErrorStatus.ne.0) then
            write(*,*) "Failed on reading vCount, BinFileErrorStatus=", &
                    BinFileErrorStatus
            return
         endif
         write(*,*)'Number of variables in Binary file, vCount=',vCount

         if (BinaryFormatVersion.eq.0) then
            NameLength = 8
         else
            read(BinaryUnit, IOSTAT = BinFileErrorStatus) NameLength
         endif
      endif
   
      if (BinFileErrorStatus.ne.0) then
          write(*,*) "BinaryReadOpenFile failed, error status:", &
                    BinFileErrorStatus
      endif   

      BinaryReadOpenFile =  (BinFileErrorStatus.eq.0)
      BinaryFile_Open = BinaryReadOpenFile 
     
      return
   end function
      
      
! ----------------------------------------------------------------------
      
      logical function BinaryReadNames (vNameArray)
      character*(*) vNameArray (vCount)

		integer*4 i

		BinaryReadNames=.false.

		if (.not. BinaryFile_Open) then
			return
		endif

      write(*,*)' BinaryReadNames, vCount=', vCount

      if (NameLength.gt.len(vNameArray(1))) then
			write(*,*) Namelength, 'character variable names in Binary file;', &
					' vNameArray character length in program code', &
               ' must be increased to at least this size'
			return
		endif

		read(BinaryUnit, IOSTAT = BinFileErrorStatus) &
         (vNameArray(i)(1:NameLength),i=1,vCount)

      if (BinFileErrorStatus.ne.0) then
			write(*,*)' Error on reading variable names from binary file, Error status=',BinFileErrorStatus
			return
		endif      
!      if (BinaryReadNames) write (667,*) ' vNameArray:', 
!     &    (vNameArray(i)(1:8),i=1,vCount)

      if (BinaryFormatVersion.eq.0) then
			if (vNameArray(1).ne.'BASSENG' .or. vNameArray(2).ne.'Lag') then
         	write(*,*) 'Binary file had no text header, but names did not conform to version 0'
				BinaryFormatVersion = -1  ! undefined
				return
			endif
		endif
			
      BinaryReadNames=.true.
      end function
      
! ----------------------------------------------------------------------

      logical function BinaryReadValues ( vCount, vArray, BinaryReadError)
      integer*4 vCount
      real*8    vArray (vCount)
      integer*4 BinaryReadError
      
      integer i
      real*4 vSingle

		if (.not. BinaryFile_Open) then
			BinaryReadValues =.false.
			return
		endif

      if (BinFileErrorStatus.eq.0) then
      	if (BinaryFormatVersion.gt.0) then
	      	read(BinaryUnit, IOSTAT = BinFileErrorStatus) vArray(1:vCount)
	else
		do i=1,vCount
		      	read(BinaryUnit, IOSTAT = BinFileErrorStatus) vSingle
					vArray(i) = vSingle
		enddo
        endif	
    	if (BinFileErrorStatus.gt.0) then
      		write(*,*)'Error in call to BinaryReadValues:', &
     			'BinFileErrorStatus=',BinFileErrorStatus
      	endif
      endif
				
      BinaryReadValues = (BinFileErrorStatus.eq.0)     
		BinaryReadError  = BinFileErrorStatus
!      if (BinaryReadValues) write (667,*) ' vArray:',vArray
      return
      end function

end module            
