Module Binary_Convert_Sub_Library

   implicit none

      integer DiagUnit /0/  ! will use default output unit
      character*(*), parameter :: diagFile = "BinConvDiagnostics.txt"
      							! Used if DiagUnit>0

   contains




      ! ======================================================================
      ! Subroutines and functions:

      ! ---------------------------------------------------------------------
      ! Get info about modify time of file according to name:

      subroutine t_Days_FileName(FileName, ModifyTime)
      implicit none
      character*(*) FileName
      real*8 ModifyTime

      integer*4 istat, infoArray(13)
      integer*4, external :: stat
      istat=stat(FileName, InfoArray) ! function in Unix library coming with ABSOFT
      call ConvertTime(istat, InfoArray(10), ModifyTime)
      end subroutine



      ! ---------------------------------------------------------------------
      ! Get info about modify time of file unit:

      subroutine t_Days_FileUnit(FileUnit, ModifyTime)
      implicit none
      integer*4 FileUnit
      real*8 ModifyTime

      integer*4 istat, infoArray(13)
      integer*4, external :: fstat, time
!      character*12 dateString, timeString
      
      istat=fstat(FileUnit, InfoArray) ! Unix library coming with ABSOFT
      call ConvertTime(istat, InfoArray(10), ModifyTime)
      end subroutine


      ! ---------------------------------------------------------------------
      ! Convert time to real local Windows time: 

      subroutine ConvertTime(istat, sTime, ModifyTime)
      implicit none
      integer*4 istat ! =0 if successful
      integer*4 sTime ! seconds since 00:00:00 GMT January 1 1970
                      ! (As in function time in UNIX library)
      real*8 ModifyTime


      integer*4 GMTarray(9), localTarray(9), SecsAdded, DateDiff
                !sec, min, hours, day, month(0-11),year(0=1900)

      integer*4 T
                          
      if (istat.eq.0) then

                      ! stat/fstat in calling routine was successful:
                      ! Time from stat and fstat is system GMT time,
                      ! adjusted for locale. Windows file time is
                      ! in local units, and the GMT must be adjusted 
                      ! for the difference to be correctly
                      ! represented in Statistica:

         call gmtime(stime,GMTarray)
         call ltime (stime,localTarray)
         

                      ! Correction of time:
         SecsAdded = localTarray(1)-GMTarray(1)          &
                   + (localTarray(2)-GMTarray(2))*60     &
                   + (localTarray(3)-GMTarray(3))*3600

                      ! Possible shift of date of +/- 1 day.
                      ! Most significant difference deternmines sign:
         DateDiff = localTarray(6)- GMTarray(6)
         if (DateDiff.eq.0) DateDiff = localTarray(5)- GMTarray(5)
         if (DateDiff.eq.0) DateDiff = localTarray(4)- GMTarray(4)
         SecsAdded = SecsAdded + DateDiff*3600*24

                     
                      ! Modifies time from seconds to days, and
                      ! adds the date value in days for jan 1 00:00:00 1970
                      ! used in Windows system.

         T = sTime+SecsAdded
         ModifyTime = float(T)/(3600*24) + 25569.0
!         write(DiagUnit,*) ' sTime:', sTime
!         write(DiagUnit,*) '     T:', T
!         write(DiagUnit,*) ' ConvertedTime :',ModifyTime 
         
      else

         ModifyTime = 0.

      endif


      ! call TimeTest(stime)
      ! call TimeTest(0)
      ! call TimeTest(time())
      ! call jdate (dateString, timeString, "+0100")
      ! write(DiagUnit,'(3(1x,A))') 'now: ',datestring, timestring      

      end subroutine



      ! ---------------------------------------------------------------------
      ! Subroutine for testing time routines:
      
      subroutine TimeTest(stime)
      implicit none
      integer*4 stime
      
      integer*4 tarray(9)
      character*24, external :: ctime
      character*8 timename(2) /'gmtime','localtime'/
      integer i

      write(DiagUnit,*) ' stime=',stime
      do i=1,2
         Select case (i)
           case (1); call gmtime(stime,tarray)
           case (2); call ltime(stime,tarray)
         end Select
         write(DiagUnit,'(1x,A,9I4)') timename(i), tarray
      end do
      write(DiagUnit,*) 'ctime = ',ctime(stime)
      end subroutine

      ! ---------------------------------------------------------------------
      ! Transform 'C' null-terminated string 
      ! into trimmed uppercase Fortran Character:

      subroutine Normalise(TextString)
      implicit none
      character*(*) TextString
      
      integer nullpos, endpos, firstchar

      nullpos = index(TextString,CHAR(0))          ! position of terminating null
      if (nullpos.gt.0) then
         endpos = len_trim(TextString(1:nullpos-1))
      else
         endpos = len_trim(TextString)
      endif

      firstchar = verify(TextString(1:endpos),' ') ! position of first non-blank
      if (firstchar.gt.0) then                  ! (=0 if no blanks)
         TextString=TextString(firstchar:endpos)
         call UpperCase(TextString)
      else
         TextString=' '
      endif      

      end subroutine


      ! ---------------------------------------------------------------------
      ! Transform text string into a trimmed, uppercase 
      ! 'C-type' null-terminated string, (replace any existing null)
      ! return pointer to string as integer*4:

      integer*4 function CString_p(TextString)
      implicit none
      character*(*), intent(inout) :: TextString

      integer nullpos, endpos, firstchar

            ! account for possible existing terminating null:
      nullpos = index(TextString,CHAR(0))

      if (nullpos.gt.0) then
         endpos = len_trim(TextString(1:nullpos-1))
      else
         endpos = len_trim(TextString)
      endif

      firstchar = verify(TextString(1:endpos),' ') ! position of first non-blank
      if (firstchar.gt.0) then                     ! (=0 if no blanks)
         TextString = TextString(firstchar:endpos)
         call UpperCase(TextString)
      else
         TextString=' '
      endif

         ! Truncate last character if needed to insert a Null character:      
         ! and return address to string:

      endpos = min(endpos - firstchar + 1,len(TextString)-1)
      TextString(endpos+1:endpos+1) = CHAR(0)
      CString_p  = LOC(TextString)
      end function




           ! transform lower case letters a-z to uppercase

      subroutine UpperCase(S)
      implicit none
      character S*(*), C*1
      integer i, IC
      do i=1, len_trim(S)
         C=S(i:i)
         if (C .ge.'a' .and. C .lt.'z') then
            IC=ICHAR(C)
            S(i:i)=CHAR(IC-ICHAR('a')+ICHAR('A'))
         else
            Select case (C)
              case ('æ'); S(i:i)='Æ'
              case ('ø'); S(i:i)='Ø'
              case ('å'); S(i:i)='Å'
            end Select
         endif
      end do   
      end subroutine


     
      
end Module Binary_Convert_Sub_Library