! =====================================================================
! Toolbox for free-form input from text files.
! Opening call returns a logical argument OK = .true. if successful.
! Later calls that read input updates an argument OK to .false.
! if an error condition is detected, but otherwise leaves it unchanged.
! This permits the calling program to use a logical variable
! as a cumulative error condition status
! (.true. until the first error, staying = .false. afterwards).


    
module sub_TextInput

      implicit none

$define Debug_TextInput

$if defined Debug_TextInput
$undefine Debug_TextInput_Detail
      integer, parameter :: Debug_file=0 ! standard error unit in Absoft FORTRAN
      logical :: Debug = .false.         ! debug messages on/off
$endif
      
      character*2, private :: CommentChars ="#!"

      integer, parameter, private :: lgtDelim = 20
      type TextInputFileDescriptor
         private
            logical   :: IsDefined = .false. 
            integer   :: Unit
            character :: Role       *64 
            character :: Name       *256
            logical   :: IsOpen
            integer   :: LineLength
            logical   :: HasLine
            logical   :: NewLine
            logical   :: EndOfLine
            logical   :: HasWord
            integer   :: LineNumber
            character :: LineContent*256
            integer   :: searchPos, firstPos, lastPos
            logical   :: RepeatWord
            character :: Delimiter_List*(lgtDelim+1)
            integer   :: ActiveDelimiters
            integer   :: Delimiter_Number
            character :: Quote
            integer   :: ReadError
            integer   :: FormatError
            character :: Message*256
            logical   :: HasMessage =.false.
            integer   :: DiagFileUnit = 0
            logical   :: EchoInput = .false.
            logical   :: LineEchoed
            logical   :: Echo_To_StdOut
      end type


contains


   subroutine TextInp_Debug(Debug_Switch)
      logical Debug_Switch
  
$if defined Debug_TextInput
      Debug = Debug_Switch    ! debug messages on/off
      
$endif

   end subroutine

! =====================================================================
! Open file for input if not already open, and initiate file descriptor.
! If the unit number refers to a file that is already open,
! the open status is checked

   logical function TextInp_OpenFile(fd, Unit, Role, Name, Delimiters)
      type (TextInputFileDescriptor) :: fd
      integer        :: Unit 
      character *(*) :: Role, Name, Delimiters

      logical OK


            ! Used to check status of specified unit:
            
      Character*3 sequential, read, formatted

$if defined Debug_TextInput
      if (Debug) write(Debug_file,*)'  >>> called TextInp_OpenFile, argument Name=',trim(Name),' Unit=', Unit
$endif


               ! transfer specifications to descriptor
               ! and check that vital items are complete:



      fd%Unit = Unit
      fd%Role = Role
      fd%Name = Name
      fd%Delimiter_List = trim(Delimiters)
      fd%ActiveDelimiters = len_trim(Delimiters)
      fd%IsOpen = .false.
      fd%IsDefined = .true.  ! even if opening file fails
      OK = .true.
      fd%Message=""
      
      if (fd%Name /= Name) then
         fd%Message='File name "'//Name//'" is too long.'
         OK = .false.
      endif

      if (len_trim(fd%Delimiter_List)>lgtDelim) then
         fd%Message = fd%Message// &
                  ' Delimiter list "'//trim(Delimiters)//'" is too long.'
         OK = .false.
      endif
      
      if (scan(trim(fd%Delimiter_list)," '"""//CommentChars)>0) then
         fd%Message = fd%Message// &
                  ' Delimiter list "'//trim(Delimiters)// &
                  '" is illegal, contains blanks, quotes or Comment characters'
         OK = .false.
      endif

      if (OK) then

               ! Check file status.
               ! The unit should either be already open for formatted sequential input,
               ! (Will be the case for standad input input 5)
               ! or be available for opening the named file.
               ! In any case, the name given as function argument
               ! will be noted and used


         INQUIRE (UNIT = unit, OPENED=fd%IsOpen, &
                  SEQUENTIAL=sequential, READ=read, FORMATTED=formatted,&
                  IOSTAT=fd%ReadError)      

         if (fd%IsOpen .and. fd%ReadError==0) then
            if(.not. (sequential=='YES' .and. read=='YES' .and. formatted=='YES') ) then
               write(fd%Message,*)'Specified unit for Text input already open', &
                  'but not for sequential formatted input'
               OK = .false.
            endif
            
         else


            OPEN( fd%Unit, File = fd%Name, status='OLD', IOSTAT = fd%ReadError)
   
            if (fd%ReadError.eq.0) then
               write(fd%Message,*) 'Opened ', fd%Role, ' ', fd%Name, &
                        ' for text input'
               fd%ActiveDelimiters = len_trim(fd%Delimiter_List)
               fd%LineNumber   = 0   ! nothing read so far 
            else
               write(fd%Message,*) 'Error code =', fd%ReadError, &
                  ' when opening ', fd%Role, ' ', fd%Name, &
                  ' as Fortran unit ', fd%unit, 'for text input.'
               OK = .false.
            endif

         endif
            
      endif

               ! Final error message common to all errors:   

      if (.not. OK) then
         fd%Message =  fd%Message // ' Further read operations will be ignored'
      endif

      fd%HasMessage = .true.
      fd%IsOpen = OK
      fd%HasLine = .false.
      fd%DiagFileUnit = 0
      fd%EchoInput = .false.
      fd%Echo_To_StdOut = .false.
      
      TextInp_OpenFile = OK

$if defined Debug_TextInput
      if (Debug) write(Debug_file,*)'  <<< return with TextInp_OpenFile=',TextInp_OpenFile
$endif

   end function

   
! =====================================================================
! return file name from Text File descriptor

   Character*256 function TextInp_GetFileName( fd )
      type (TextInputFileDescriptor) :: fd
      TextInp_GetFileName = fd%Name
   end function

   integer function TextInp_GetFileUnit( fd )
      type (TextInputFileDescriptor) :: fd
      TextInp_GetFileUnit = fd%Unit
   end function

   
! =====================================================================
! Check if End of file has been reached for open file,
! signal by returning .true.

   logical function TextInp_EOF( fd )
      type (TextInputFileDescriptor) :: fd
      TextInp_EOF = fd%IsOpen .and. (fd%Readerror < 0)
   end function


! =====================================================================
! Check if End of line has been reached for open file,
! signal by returning .true.

   logical function TextInp_EOL (fd)
   type (TextInputFileDescriptor) :: fd
   TextInp_EOL = fd%IsOpen .and. fd%HasLine &
                .and. fd%EndOfLine .and. (.not. fd%repeatWord)
   end function


! =====================================================================
! Close input file and inactivate file descriptor from further reading:

   subroutine TextInp_CloseFile( fd )
      type (TextInputFileDescriptor) :: fd
      if (fd%IsOpen) then
         close(fd%unit,IOSTAT=fd%Readerror)
         fd%IsOpen=.false.
      endif
   end subroutine

! =====================================================================
! read to next non-blank line, counting also skipped lines:

   logical function TextInp_Get_Line(fd)
      type (TextInputFileDescriptor) :: fd
      
      integer :: i, sP

$if defined Debug_TextInput
      if (Debug) write(Debug_file,*)'  >>> called TextInp_Get_Line, fd%Unit=', fd%Unit
$endif

      fd%HasLine = .false.

      if (fd%HasMessage) then
         fd%Message=""
         fd%HasMessage = .false.
      endif

      if (.not.fd%IsOpen) then
         TextInp_Get_Line=.false.
      else

         do
            if(fd%Echo_to_StdOut .or. fd%Unit==5 .or. fd%Unit==0) then
               write(0,'(A)',ADVANCE='NO')">"
            endif
            if(fd%EchoInput .and. fd%DiagFileUnit.ne.0 .and. fd%DiagFileUnit.ne.6) then
               write(fd%DiagFileUnit,'(A)',ADVANCE='NO')">"
            endif
      
            READ(fd%unit,'(A)',IOSTAT = fd%ReadError )  fd%LineContent
   
            IF ( fd%ReadError > 0) THEN
               write(fd%Message,*) ' Error number ', fd%ReadError, &
                                   ' on reading inputline ', fd%LineNumber+1
            elseif ( fd%ReadError < 0) then
               write(fd%Message,*) ' Passed end of file after ', fd%LineNumber,' lines'
            endif
   
            if (fd%ReadError.ne.0) then
               fd%HasMessage = .true.
               fd%HasLine    = .false.
               exit  ! also for EOF
            else
               fd%HasLine    = .true.
               fd%LineEchoed = .false.
            endif
   
         
               !  ........ skip blank lines and comment lines (beginning with
               !           one of characters in string CommentChars (# OR !)
         
            fd%LineNumber = fd%LineNumber + 1 
            fd%LineLength = len_trim(fd%LineContent)
            if (fd%EchoInput) call TextInp_WriteLine(fd)
   
            if (fd%LineLength>0) then 
               sP = verify(fd%LineContent," ") ! start at first nonblank
                  
               i = index(CommentChars,fd%LineContent(sP:sP))
               if (i == 0) then !line is not a comment line:
                  fd%NewLine    = .true.
                  fd%EndOfLine  = .false.
                  fd%RepeatWord = .false.
                  exit
               endif
            endif 
   
         end do
   
         fd%searchPos = sP
         fd%HasWord = .false.
         TextInp_Get_Line = (fd%ReadError.EQ.0)
      endif

$if defined Debug_TextInput
      if (Debug) then
         write(Debug_file,*)'  <<< returns with TextInp_Get_Line=', TextInp_Get_Line, &
            '  fd%Isopen=',fd%IsOpen, '  fd%Readerror=',fd%Readerror, ' fd%HasLine=', fd%HasLine

         if (fd%Readerror.eq.0) write(Debug_file,*) 'LineContent:', trim(fd%LineContent)
      endif
$endif

   END function



! =====================================================================
! Free-form read of real*8 value, selecting number AltNum 
! of values entered as sequence y1|y2|...|yn ,
! with | as delimiter between alternatives.
! Any other delimiter ends sequence and causes the read function to
! to use the last value if there are fewer than AltNum values.
! (Value sequence requires | to be defined as delimiter)

   logical function TextInp_ValueSelect (fd, AltNum, Value)
      type (TextInputFileDescriptor) :: fd
      integer :: AltNum
      real*4 :: Value

! --------- LOCAL VARIABLES:

      integer :: ValueCount
      LOGICAL :: OK
      character*1 :: Delimiter
      real*4    :: V

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_ValueSelect, AltNum=', AltNum
$endif
        
      ValueCount = 0
      OK = .true.

      Do 
         if (TextInp_Value (fd,V)) then

$if defined Debug_TextInput_Detail
            if (Debug) write(Debug_file,*)'      TextInp_Value succeeded, V=', V
$endif
            ValueCount = ValueCount+1
            Delimiter = TextInp_DelimiterFound(fd)

$if defined Debug_TextInput_Detail
            if (Debug) write(Debug_file,*)'    Delimiter, ValueCount = ',Delimiter, ValueCount
$endif
            if ( Delimiter .ne. "|") then        ! last value in sequence of alternatives
               if (ValueCount<=AltNum) Value = V ! use it if not after specified alternaive
               exit
            endif
            if (ValueCount == AltNum) then
               Value = V
            endif
         else
            call TextInp_StandardMessage(fd)
            OK = .false.
         endif
      end do 

      if (ValueCount>0 .and. OK) then
         TextInp_ValueSelect = OK
      else
         TextInp_ValueSelect = .false.
      endif

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  <<< returns with TextInp_ValueSelect=',TextInp_ValueSelect
      if (Debug) write(Debug_file,*)'      OK, Value=', OK, Value
$endif

   end function

! =====================================================================
! Free-form read of integer value from string, starting at current position:
! Updates file descriptor to prepare for next read:

   Logical function TextInp_Logical (fd, LogicalValue)
   type (TextInputFileDescriptor) :: fd
   logical :: LogicalValue

      integer :: fP, lP

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_Logical, fd%Unit=', fd%Unit
$endif
      

      TextInp_Logical = .false.

      if (TextInp_DelimitWord (fd)) then

               ! has found word, possibly empty (or blank)
               ! and cleared any previous error message
         fP = fd%firstPos
         lP = fd%LastPos
         
         if (fP<=lP) then
            if (indexStr_CaseNeutral("true",fd%LineContent(fP:lP))==1) then
               LogicalValue=.true.
               TextInp_Logical = .true.
            elseif (indexStr_CaseNeutral("false",fd%LineContent(fP:lP))==1) then
               LogicalValue=.false.
               TextInp_Logical = .true.
            else
               write(fd%Message,*) ' Invalid item at position ',fd%firstPos, &
                  " in inputline ", fd%LineNumber,", True/False expected"
               fd%HasMessage = .true.
               
            endif
         else
            write(fd%Message,*) ' Empty item at position ',fd%firstPos, &
                  " in inputline ", fd%LineNumber,", True/False expected"
            fd%HasMessage = .true.
         endif
      endif

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  <<< returns with TextInp_Logical=',TextInp_Logical
      if (Debug) write(Debug_file,*)'      LogicalValue=', LogicalValue
$endif


   end function

! =====================================================================
! Free-form read of integer value from string, starting at current position:
! Updates file descriptor to prepare for next read:

   Logical function TextInp_Integer (fd, IntegerValue)
   type (TextInputFileDescriptor) :: fd
   integer :: IntegerValue

      integer :: fP, lP

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_Integer, fd%Unit=', fd%Unit
$endif
      

      if (TextInp_DelimitWord (fd)) then
               ! has found word, possibly empty (or blank)
               ! and cleared any previous error message
         fP = fd%firstPos
         lP = fd%LastPos
         if (fP<=lP) then
            READ ( fd%LineContent(fP:lP),*,IOSTAT=fd%FormatError)  IntegerValue

            if (fd%FormatError == 0) then
               TextInp_Integer = .true. 
               return

            else
               write(fd%Message,*) " Error number ", fd%FormatError, &
                  " on reading '",fd%LineContent(fP:lP), &
                  " as integer in inputline ", fd%LineNumber
               fd%HasMessage = .true.
            endif
         else
            write(fd%Message,*) ' Empty item at position ',fd%firstPos, &
                  " in inputline ", fd%LineNumber,", integer value expected"
            fd%HasMessage = .true.
      endif

      TextInp_Integer = .false.

      endif

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  <<< returns with TextInp_Integer=',TextInp_Integer
      if (Debug) write(Debug_file,*)'      IntegerValue=', IntegerValue
$endif

   end function

! =====================================================================
! Free-form read of real*8 value from string, starting at current position:
! Updates file descriptor to prepare for next read:

   Logical function TextInp_Value (fd, Value)
   type (TextInputFileDescriptor) :: fd
   real*4 :: Value

      integer :: fP, lP


$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_Value, fd%Unit=', fd%Unit
$endif
 
      if (TextInp_DelimitWord (fd)) then
               ! has found word, possibly empty (or blank)
               ! and cleared any previous error message
         fP = fd%firstPos
         lP = fd%LastPos
         if (fP<=lP) then
            READ ( fd%LineContent(fP:lP),*,IOSTAT=fd%FormatError)  Value

            if (fd%FormatError == 0) then
               TextInp_Value = .true. 
               return

            else
               write(fd%Message,*) " Error number ", fd%FormatError, &
                  " on reading '",fd%LineContent(fP:lP), &
                  " as numerical value in inputline ", fd%LineNumber
               fd%HasMessage = .true.
            endif
         else
            write(fd%Message,*) ' Empty item at position ',fd%firstPos, &
                  " in inputline ", fd%LineNumber,", numerical value expected"
            fd%HasMessage = .true.
      endif

      TextInp_Value = .false.

      endif
$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  <<< returns with TextInp_Value=',TextInp_Value
      if (Debug) write(Debug_file,*)'      Value=', Value
$endif

   end function


! =====================================================================
! Return next item as text string, possibly blank

   logical function TextInp_Word (fd, Word)
      type (TextInputFileDescriptor) :: fd
      character*(*) Word
      
      integer fP, lP

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_Word, fd%Unit=', fd%Unit
$endif
      
      if (TextInp_DelimitWord(fd)) then
         fP = fd%firstPos
         lP = fd%LastPos
         if (fP<=lP) then
            Word = fd%LineContent(fP:lP)
         else
            Word = ""
         endif
         TextInp_Word = .true.
      else
         TextInp_Word = .false.
      endif

$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  <<< returns with TextInp_Word=',TextInp_Word
      if (Debug) write(Debug_file,*)'      Word=', Word
$endif

   end function

! =====================================================================
! prepare for repeated reading of last item,
! will skip scanning for new item:

   subroutine TextInp_Repeat(fd)
      type (TextInputFileDescriptor) :: fd
      if (fd%IsOpen .and. fd%HasWord) then
         fd%RepeatWord=.true.
      endif
   end subroutine 

! =====================================================================
!    Find item in current input line, delimited by blanks
!    and/or one of the defined delimiters, 
!    possibly enclosed by single or double quotes,
!    Notes starting and ending position of the word found,
!    and beginning position for next word.
!    returns .true. if item found, .false. if no more items in the line

   logical function TextInp_DelimitWord (fd)
      type (TextInputFileDescriptor) :: fd

      INTEGER sP, eP, i
      CHARACTER FirstChar*1

! -----------------------------------------------------------------


$if defined Debug_TextInput_Detail
      if (Debug) write(Debug_file,*)'  >>> called TextInp_DelimitWord, fd%Unit=', fd%Unit
$endif

      if ( .not. (fd%IsOpen .and. fd%HasLine) ) then
         TextInp_DelimitWord=.false.
         return
      endif

      if (fd%HasMessage) then
         fd%HasMessage=.false.
         fd%Message = ""
      endif
      
      if (fd%HasWord .and. fd%RepeatWord) then
         fd%RepeatWord = .false.
         TextInp_DelimitWord = .true.
         return
      endif


      sP = fd%searchPos
      eP = fd%LineLength
      if (sP >eP) then
         TextInp_DelimitWord = .false.
         return
      endif

      FirstChar = fd%LineContent(sP:sP)

      if (FirstChar=="'" .or. FirstChar == '"') then

         fd%Quote = FirstChar

               !................................................
               ! first character is a quote; scan for ending quote
               ! (use rest of line if ending quote is not found)
               ! Note: the quote can be included in the quoted string,
               ! by entering two consequtive quote characters;
               ! they will be truncated into one quote character

         fd%firstPos = sP+1
         if (sP>=eP) then
            fd%lastPos = sP ! line ended by empty string, indicated by 
                            ! starting quote as last character.
         else
            Do
               i = scan(fd%LineContent(sP+1:eP),FirstChar)
               if (i <= 0) then      ! *** Searched to end of line without finding quote
                  sP = eP+1          ! set search position after end of line.
                  fd%lastPos = eP    ! last character at end of line
               else                  ! *** Quote found,
                  sP = sP+i+1        ! set search position after ending quote
                  if (sP <= eP) then ! if within lneline; check for repeated quote
                     if (fd%LineContent(sP:sP) == FirstChar) then
                                    ! repeated quote found: keep only the first 
                                    ! as literal quote in the string, 
                        fd%LineContent(sP-1:eP) = fd%LineContent(sP:eP)
                        eP = eP-1
                        cycle       ! and continue search for ending quote:
                     endif
                  endif
                  fd%lastPos = sP-2 ! Otherwise note last position before quote as end of word 
               endif
               exit                 ! last position set; exit search loop
            end do

          endif
               ! sP points to position after ending quote,
               ! or just after end of line.
               ! firstPos,lastPos points to beginnng and end of string.
               ! lastPos will be = firstPos-1 for empty string                            

      else

         fd%Quote = " "
 
               !................................................
               ! unquoted string: search for delimiter, blank
               ! or Comment character:

         fd%firstPos = sP

     
         i = scan( fd%LineContent(sP:eP), &
                   fd%Delimiter_List(1:fd%ActiveDelimiters)//" "//CommentChars) - 1
         if (i < 0) then
            sP = eP+1
         else
            sP = sP+i      ! points to position just after end of item
         endif
         fd%lastPos = sP-1
      endif     
          
               !................................................
               ! End of item is established;
               ! now find beginning of next input item:
      fd%HasWord = .true.

      if (sP <= eP) then
      
                     ! check if comment character:
         if (index(CommentChars,fd%LineContent(sP:sP)) > 0) then
            eP = sP-1
         else   
            i = verify(fd%LineContent(sP:eP), " ")  ! search for non-blank
            if (i<=0) then
               sP = eP+1      ! not found: position set to after end of line
            else
               sP = sP+i-1    ! = Position of nonblank
               fd%Delimiter_Number = index( fd%Delimiter_List(1:fd%ActiveDelimiters), &
                        fd%LineContent(sP:sP ) )
               if (fd%Delimiter_Number > 0) then
                              ! Position to next non-blank
                              ! after nonblank delimiter:
                  i = verify(fd%LineContent(sP+1:eP), " ")
                  if (i<=0) then
                     sP = eP+1
                  else
                     sP = sP+i
                  endif
               else 
                  if (i==1 .and. fd%Quote /= " ") then
                     write(fd%Message,*) &
                        "String in quotes must be followed by blank or delimiter"
                     fd%FormatError = .true.   
                     fd%HasMessage = .true.
                  endif
               endif
            endif
         endif
      else
         fd%EndOfLine = .true.
      endif

      fd%searchPos = sP
      fd%NewLine = .false.

      TextInp_DelimitWord = .true.
      
   end function

! =====================================================================
! Specify how many nonblank delimiters should be recognized:
   subroutine TextInp_DelimitersUsed (fd, DelimiterCount) 
      type (TextInputFileDescriptor) :: fd
      integer :: DelimiterCount

      fd%ActiveDelimiters = max(1, min ( DelimiterCount, &
                                         len_trim(fd%Delimiter_List)))
   end subroutine

! =====================================================================
! Returns the delimiter that was found after last input item:
   character*1 function TextInp_DelimiterFound (fd)
      type (TextInputFileDescriptor) :: fd
      integer nd

      nd = fd%Delimiter_Number
      if (nd>0) then
         TextInp_DelimiterFound = fd%Delimiter_List(nd:nd)
      else
         TextInp_DelimiterFound = " "
      endif

   end function

! =====================================================================
! Prepare for writing to diagnostic file:

   subroutine TextInp_DiagFileUnit (fd,DiagFileUnit, EchoInput, Echo_to_StdOut)

      type (TextInputFileDescriptor) :: fd
      integer :: DiagFileUnit
      logical :: EchoInput
      logical :: Echo_to_StdOut
      
      fd%DiagFileUnit = DiagFileUnit
      fd%EchoInput = EchoInput
      fd%LineEchoed = .false.
      fd%Echo_to_StdOut = Echo_to_StdOut

   end subroutine

! =====================================================================
! write last input line to specified file unit if not done already,
! and to standard error unit if that fails:

   subroutine TextInp_WriteLine( fd )
      type (TextInputFileDescriptor) :: fd
      
      integer :: ErrorCode, nf

      if (fd%IsOpen .and. fd%HasLine ) then
         if (fd%LineEchoed) return

         nf = fd%DiagFileUnit
         do
            write(nf,*, IOSTAT=ErrorCode) "Line ",fd%LineNumber,": ",trim(fd%LineContent)
            if (ErrorCode .ne. 0) then
               write(*,*)"TextInp_WriteLine failed to write input line to unit no. ", nf
            endif
            if ( .not. ((fd%Echo_To_StdOut .or. ErrorCode.ne.0) .and. nf.ne.0) ) exit
            nf=0
         end do
   
         fd%LineEchoed = .true.
      endif

   end subroutine

! =====================================================================
! write Standard message created within module

   subroutine TextInp_StandardMessage(fd) 
      type (TextInputFileDescriptor) :: fd

      if (fd%HasMessage ) then
         call TextInp_Message(fd, fd%Message)  
      endif
   end subroutine

     
! =====================================================================
! write external message to specified DiagFileUnit unit,
! and to standard error unit if that fails,
! or if it is specified in second argument:

   subroutine TextInp_Message( fd, MessageText)
      type (TextInputFileDescriptor) :: fd
      character*(*) :: Messagetext

      integer :: ErrorCode, nf

      call TextInp_WriteLine(fd)

      if (fd%IsDefined ) then

         nf = fd%DiagFileUnit
   
         do
            write(nf,'(A,A,A,A/A)', IOSTAT=ErrorCode) &
                 trim(fd%Role)," ",trim(fd%Name),":", trim(MessageText)
            if (ErrorCode .ne. 0) then
               write(*,*)"TextInp_Message failed to write error message to unit no. ", nf
            endif
            if ( .not. ((fd%Echo_To_StdOut .or. ErrorCode.ne.0) .and. nf.ne.0) ) exit
            nf=0
         end do

         fd%HasMessage = .false.  ! can only be written once

      endif
   end subroutine

! =======================================================================
! Utility functions for character array comparison;
! Check if arrays of single characters are equal,
! ignore differences of case and trailing blanks.
! NOTE: Works only if both arrays start at index 1.

   logical function equalCharArray_CaseNeutral( Array1, Array2)
      character*1 :: Array1(:), Array2(:)
   
      integer i, first1, last1, first2, last2, Diff

      first1 = lbound(Array1,1)
      last1  = ubound(Array1,1)
      first2 = lbound(Array2,1)
      last2  = ubound(Array2,1)
      Diff = first2-first1

      equalCharArray_CaseNeutral= .false.  ! assumed until proven true

            ! Scans from start of both arrays, accounting for index shift.
            ! if one array is longer, it is checked that the excess part
            ! only contains blanks.

      i=first1
      do
         if (i <=last1) then  ! within length of Array1
            if (i+Diff<=last2) then  ! within length of Array2   
               if (.not. equal_CaseNeutral(Array1(i), Array2(i+Diff)) ) return
            else
               if (Array1(i).ne." ") return
            endif
         elseif (i+Diff<=last2) then  ! within length of Array2   
            if (Array2(i+Diff).ne." ") return
         else     ! scanned both arrays completely without any differences
                  ! except possibly in length of trailing blank sequence:
            equalCharArray_CaseNeutral= .true.
            return
         endif
         i = i+1
      end do

      return

   end function   
   
! =======================================================================
! Utility functions for string comparison;
! Check if strings are equal, ignore differences of case:

   logical function equalStr_CaseNeutral( Str1, Str2)
      character*(*) :: Str1, Str2
   
      integer i
      integer LgtStr
      
      LgtStr = len_trim(Str1)
      equalStr_CaseNeutral= .false.  ! assumed until proven true
   
      if (LgtStr.eq.len_trim(Str2)) then
         do i=1,LgtStr
            if (.not. equal_CaseNeutral(Str1(i:i), Str2(i:i)) ) return
         end do
         equalStr_CaseNeutral= .true. 
         return
      endif
      return
   end function   


! =======================================================================
! Utility functions for string comparison;
! find substring in string at any position, return index:

   integer function indexStr_CaseNeutral( Str, SubStr)
   character*(*) :: Str, SubStr

   integer i,k
   integer LgtSub
   
   LgtSub = len_trim(SubStr)

   do i=0,len_trim(Str)-LgtSub
      k=1
      do 
         if (.not. equal_CaseNeutral(Str(i+k:i+k), SubStr(k:k)) ) exit
         k = k + 1
         if (k > LgtSub) then
            indexStr_CaseNeutral = i+1
            return
         endif
      end do
   end do
   indexStr_CaseNeutral = 0
   end function


! =======================================================================
! Checks if two characters are equal, ignoring case differences
! for letters in Norwegian alphabet:

   logical function equal_CaseNeutral(Ch1, Ch2)
   character*(1) :: Ch1, Ch2
   integer i1, i2
   if (Ch1==Ch2) then
      equal_CaseNeutral=.true.
   else
      i1 = ichar(Ch1)
      i2 = ichar(Ch2)
      if (i1+32 == i2) then
         equal_CaseNeutral = isCapitalLetter(Ch1)
      elseif (i1 == i2+32) then
         equal_CaseNeutral = isCapitalLetter(Ch2)
      else
         equal_CaseNeutral = .false.
      endif
   endif
   end function


! =======================================================================
! returns true if character in Norwegian alphabet in standard ASCII character set
! as defined in the ABSOFT Pro Fortran 11.1 User Guide.
! Will not work with characters in DOS Nordic Code Page.

   logical function isCapitalLetter(Ch)
   character*1 Ch
   isCapitalLetter = (Ch >= "A" .and. Ch <= "Z") &
              .or. (Ch == "Æ") .or. (Ch == "Ø") .or. (Ch == "Å")
   end function

  end module
