Module ModelParamToolbox

$undefine DEBUG_PRINT

		implicit none

		! ------------- Values for Phase argument in ModelParam_Handler
      !           --- for Parameter setup as part of model initiation:
		integer, parameter :: NotStarted = 0, &
                            Counting   = 1, &
                            Building   = 2, & 
                            Complete   = 3
      !           --- for later model operation:
		integer, parameter :: Writing    = 4

      ! ------------- Descriptive texts for diagnostics:
      character*(10)     :: PhaseDescr(NotStarted:Writing) = &
				             (/"NotStarted", &
                           "Counting  ", &
                           "Building  ", &
                           "Complete  ", &
                           "Writing   " /)
		
		! ------------- Derived types and descriptors for building 
      !               the model parameter description 

		type :: ParameterGroupDescr 
			integer            :: FirstParam
			character, pointer :: StructName(:)
			integer               :: NameLength
		end type
		type (ParameterGroupDescr), allocatable :: pGroupDescr(:)

      type :: ParameterReference
			logical     , pointer :: Logical
			integer     , pointer :: Integer
			real*4      , pointer :: real
			logical     , pointer :: LogicalArray(:)
			integer     , pointer :: IntegerArray(:)
			real*4      , pointer :: RealArray(:)
			integer     , pointer :: IntegerMatrix(:,:)
			real*4      , pointer :: RealMatrix(:,:)
			character   , pointer :: NameChar(:)
			integer               :: NameLength
      end type
      type (ParameterReference), allocatable :: ParamRef(:)

      integer, allocatable :: ParamType(:), ParamDim(:)
      integer, allocatable :: Param_LBound(:,:),Param_UBound(:,:)
      integer, parameter :: &
                  pTypeLogical      = 0, &
                  pTypeInteger      = 1, &
                  pTypeReal         = 2, &
                  pTypeLogicalArray = 3, &
                  pTypeIntegerArray = 4, &
                  pTypeRealArray    = 5, &
                  pTypeIntegerMatrix= 6, &
                  pTypeRealMatrix   = 7, &
                  pTypeCharacter    = 8

      ! ------------- Status variables for functions below:

		integer :: Phase = NotStarted
      logical :: ModelParamError = .false.

		integer GroupCount, ParamCount
		integer dimGroupDescr, dimParamRef
      integer Allocate_Status
      
      integer :: ParamFileUnit
      integer :: ParamOutputError
      character, pointer :: ParamFileNameChar(:)
      integer :: ParamFileNameLength

   contains


!=======================================================================
! Procedures used by ModelParameter group handlers ModelParam_xxxxxxx:

! The procedures will count groups or parameters,
! allocate space for name and register description and reference,
! or write to parameter output file, depending on Phase
! as set by master subroutines ModelParam_Setup and ModelParam_Output



! ----------------------------------------------------------------------
! starts new model parameter group:

		subroutine ParamGroup(TextDescr, StructName)
		character *(*) TextDescr, StructName

$if defined (DEBUG_PRINT)
      write(*,*)'Phase=',Phase,'  Param Group ', StructName, &
                ',  ModelParamError:', ModelParamError
$endif
      
		if (ModelParamError) return
		GroupCount = GroupCount + 1

		if (Phase /= Counting .and. GroupCount > dimGroupDescr) then

			write (*,*) "Error in Parameter Setup, ParamGroup in Phase =",Phase, &
                        ": Group ", StructName, &
                        " exceeds allocated space in pGroupDescr" 
			ModelParamError = .true.
			return

		endif	

		if (Phase == Building) then

			pGroupDescr(GroupCount)%FirstParam = ParamCount+1
			call SaveName(StructName, pGroupDescr(GroupCount)%StructName, pGroupDescr(GroupCount)%NameLength)

		elseif(Phase == Writing) then

			write(ParamFileUnit,*, IOSTAT = ParamOutputError) &
                  "###### Parameter Group ", StructName

         call WriteParamDescrString(TextDescr)

	      if (ParamOutputError /= 0) then
			   write(*,*)" Error when writing group heading ", StructName, &
	                   " to file '", ParamFileNameChar, "'", &
	                   "; Error code=", ParamOutputError
			   ModelParamError = .true.

	      endif
      endif

		end subroutine ParamGroup



! ----------------------------------------------------------------------
! Parameter explanation (only active for writing to file):

		subroutine ParamExpl(DescrText)
		character *(*) DescrText

$if defined (DEBUG_PRINT)
      write(*,*)'Phase=',Phase,'   ParamExpl ', DescrText, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		
      if(Phase == Writing) then
			call WriteParamDescrString(DescrText)
	      if (ParamOutputError /= 0) then
			   write(*,*)" Error when writing Parameter explanation ", &
	                   " to file '", ParamFileNameChar, "'", &
	                   "; Error code=", ParamOutputError
			   ModelParamError = .true.
	      endif
      endif

		end subroutine ParamExpl


! -----------------------------------------------------------------
! Logical parameter:

		subroutine ParamDeclLogical( Param, InitValue, Name, DescrText)
		logical, target :: Param  ! Model parameter
		logical :: InitValue      ! Coded initial value for parameter
		character*(*) :: Name, DescrText

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclInteger ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
      
		if (Phase==Building) then
			ParamType(ParamCount) = pTypeLogical
			ParamDim (ParamCount) = 0
			paramRef (ParamCount)%Logical => Param
         Param = InitValue  

	   elseif(Phase == Writing) then

		   write(ParamFileUnit,'(2A,L6)', IOSTAT = ParamOutputError, advance = 'no') &
               trim(Name)," = ", Param      

			call WriteParamDescr(Name,"", DescrText)

	   endif

	end subroutine ParamDeclLogical



! -----------------------------------------------------------------
! Integer parameter:

		subroutine ParamDeclInteger( Param, InitValue, &
		           Name, UnitText, DescrText)
		integer, target :: Param ! Model parameter
		integer :: InitValue     ! Coded initial value for parameter
		character*(*) :: Name, UnitText, DescrText

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclInteger ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
      
		if (Phase==Building) then

			ParamType(ParamCount) = pTypeInteger
			ParamDim (ParamCount) = 0
			paramRef (ParamCount)%Integer => Param
         Param = InitValue  

	   elseif(Phase == Writing) then

		   write(ParamFileUnit,'(2A,I12)', IOSTAT = ParamOutputError, advance = 'no') &
               trim(Name)," = ", Param      

			call WriteParamDescr(Name, UnitText, DescrText)
	   endif

	end subroutine ParamDeclInteger


   
! -----------------------------------------------------------------
! real*8 parameter:

	subroutine ParamDeclReal( Param, InitValue, &
	           Name, UnitText, DescrText)
		real*4, target :: Param ! Model parameter
		real*4 :: InitValue     ! Coded initial value for parameter
		character*(*) :: Name, UnitText, DescrText

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclReal ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		

      if (Phase==Building) then

			ParamType(ParamCount) = pTypeReal
			ParamDim (ParamCount) = 0
			paramRef (ParamCount)%real => Param
         Param = InitValue

		elseif (Phase == Writing) then

		   write(ParamFileUnit,'(2A,G12.6)', IOSTAT = ParamOutputError, advance = 'no') &
               trim(Name)," = ", Param      
			call WriteParamDescr(Name, UnitText, DescrText)

		endif

		end subroutine ParamDeclReal


! ----------------------------------------------------------------------
! 1-dimensional array of logical values:

		subroutine ParamDeclLogicalArray( paramArray, InitArray, &
		           Name, DescrText )
		logical, target:: paramArray(:) ! Parameter array used by model
		logical:: InitArray(:)          ! Initial values
		character*(*) :: Name, DescrText

		integer i, N, LastInit
      integer First, Last, LastInLine

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclLogicalArray ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
		
		if (Phase == Building) then

			ParamType(ParamCount) = pTypeLogical
			ParamDim (ParamCount) = 1
			Param_LBound(1,ParamCount) = lbound(paramArray,1)
			Param_UBound(1,ParamCount) = ubound(paramArray,1)
			paramRef (ParamCount)%LogicalArray => paramArray

         N = size(paramArray)
         LastInit = 	ubound(InitArray,1)
		   ParamArray = reshape(InitArray, (/N/), &
                              (/(InitArray(LastInit),i=1,N)/))

		elseif (Phase == Writing) then

         First = lbound(paramArray,1)
         Last  = ubound(Paramarray,1)
         LastInLine = min(First + 4,Last)
         		   
		   write(ParamFileUnit,'(2A,I3,A,I3,A,5L7)', IOSTAT = ParamOutputError) &
                  trim(Name),"(",First,":",Last,") = ", &
                  ParamArray(First:LastInLine)

		   if (Last > LastInLine) then
            write(ParamFileUnit,'(10x,5L7)', IOSTAT = ParamOutputError) &
                   ParamArray(LastInLine+1:Last)
         endif

         call WriteParamDescr(Name, "", DescrText)

      endif
      end subroutine ParamDeclLogicalarray

! ----------------------------------------------------------------------
! 1-dimensional array of integer values:

		subroutine ParamDeclIntegerArray( paramArray, InitArray, &
		           Name, UnitText, DescrText )
		integer, target:: paramArray(:) ! Parameter array used by model
		integer:: InitArray(:)          ! Initial values
		character*(*) :: Name, UnitText, DescrText

		integer i, N, LastInit
      integer First, Last, LastInLine

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclIntegerArray ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
		
		if (Phase == Building) then

			ParamType(ParamCount) = pTypeInteger
			ParamDim (ParamCount) = 1
			Param_LBound(1,ParamCount) = lbound(paramArray,1)
			Param_UBound(1,ParamCount) = ubound(paramArray,1)
			paramRef (ParamCount)%IntegerArray => paramArray

         N = size(paramArray)
         LastInit = 	ubound(InitArray,1)
		   ParamArray = reshape(InitArray, (/N/), &
                              (/(InitArray(LastInit),i=1,N)/))

		elseif (Phase == Writing) then

         First = lbound(paramArray,1)
         Last  = ubound(Paramarray,1)
         LastInLine = min(First + 4,Last)
         		   
		   write(ParamFileUnit,'(2A,I3,A,I3,A,5I12)', IOSTAT = ParamOutputError) &
                  trim(Name),"(",First,":",Last,") = ", &
                  ParamArray(First:LastInLine)

		   if (Last > LastInLine) then
            write(ParamFileUnit,'(10x,5I12)', IOSTAT = ParamOutputError) &
                   ParamArray(LastInLine+1:Last)
         endif

         call WriteParamDescr(Name, UnitText, DescrText)

      endif
      end subroutine ParamDeclIntegerarray


! ----------------------------------------------------------------------
! 1-dimensional array of real*8 values:

		subroutine ParamDeclRealArray( paramArray, InitArray, &
		           Name, UnitText, DescrText )
		real*4, target:: paramArray(:) ! Parameter array used by model
		real*4:: InitArray(:)          ! Initial values
		character*(*) :: Name, UnitText, DescrText

		integer i, N, LastInit
      integer First, Last, LastInLine

$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclRealArray ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
		
		if (Phase == Building) then

			ParamType(ParamCount) = pTypeReal
			ParamDim (ParamCount) = 1
			Param_LBound(1,ParamCount) = lbound(paramArray,1)
			Param_UBound(1,ParamCount) = ubound(paramArray,1)
			paramRef (ParamCount)%RealArray => paramArray

         N = size(paramArray)
         LastInit = 	ubound(InitArray,1)
		   ParamArray = reshape(InitArray, (/N/), &
                              (/(InitArray(LastInit),i=1,N)/))

		elseif (Phase == Writing) then

         First = lbound(paramArray,1)
         Last  = ubound(Paramarray,1)
         LastInLine = min(First + 4,Last)
         		   
		   write(ParamFileUnit,'(2A,I3,A,I3,A,5(1X,G13.6))', IOSTAT = ParamOutputError) &
                  trim(Name),"(",First,":",Last,") = ", &
                  ParamArray(First:LastInLine)

		   if (Last > LastInLine) then
            write(ParamFileUnit,'(10x,5(1X,G13.6))', IOSTAT = ParamOutputError) &
                   ParamArray(LastInLine+1:Last)
         endif

         call WriteParamDescr(Name, UnitText, DescrText)

      endif

      end subroutine ParamDeclRealarray

! ----------------------------------------------------------------------
! 2-dimensional array of integer values:

		subroutine ParamDeclIntegerMatrix( paramArray, InitArray, &
		           Name, UnitText, DescrText )
		integer, target:: paramArray(:,:) ! Parameter array used by model
		integer:: InitArray(:,:)          ! Initial values
		character*(*) :: Name, UnitText, DescrText

		integer N(2), LastInitColumn, FirstRow
      integer FirstColumn, LastColumn, LastInLine
      integer i, k, kk
      
$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclIntegermatrix ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
		
		if (Phase == Building) then

			ParamType(ParamCount) = pTypeReal
			ParamDim (ParamCount) = 2
			do k=1,2
            Param_LBound(k,ParamCount) = lbound(paramArray,k)
			   Param_UBound(k,ParamCount) = ubound(paramArray,k)
         enddo
			paramRef (ParamCount)%IntegerMatrix => paramArray

            ! Fill each row of ParamArray from corresponding row in InitArray
            ! Padding each row with value in last column of InitArray,
            ! and repeating last row of Initarray if needed.

         N = shape(paramArray)
         LastInitColumn = 	ubound(InitArray,1)
         FirstRow = lbound(ParamArray,2)
         do k = FirstRow, ubound(ParamArray,2)
            kk = min(ubound(InitArray,2),lbound(InitArray,2) + (k-FirstRow))
		      ParamArray(:,k) = reshape( InitArray(:,kk), (/N(1)/), &
                                       (/(InitArray(LastInitColumn,kk),i=1,N(1))/))
         end do

		elseif (Phase == Writing) then

	      FirstColumn = lbound(paramArray,1)
	      LastColumn  = ubound(Paramarray,1)
         LastInLine  = min(FirstColumn + 4,LastColumn)
         		   
         do k = lbound(ParamArray,2), ubound(ParamArray,2)
			   write(ParamFileUnit,'(2A,3(I3,A),5I12)', IOSTAT = ParamOutputError) &
	               trim(Name), "(", FirstColumn, ":", LastColumn, "," , k,") = ", &
                  ParamArray(FirstColumn:LastInLine,k)
			   if (LastColumn > LastInLine) then
	            write(ParamFileUnit,'(5I12)', IOSTAT = ParamOutputError) &
	               ParamArray(LastInLine+1:LastColumn,k)
	         endif
         end do

	      call WriteParamDescr(Name, UnitText, DescrText)

      endif

      end subroutine ParamDeclIntegerMatrix



! ----------------------------------------------------------------------
! 2-dimensional array of real*8 values:

		subroutine ParamDeclRealMatrix( paramArray, InitArray, &
		           Name, UnitText, DescrText )
		real*4, target:: paramArray(:,:) ! Parameter array used by model
		real*4:: InitArray(:,:)          ! Initial values
		character*(*) :: Name, UnitText, DescrText

		integer N(2), LastInitColumn, FirstRow
      integer FirstColumn, LastColumn, LastInLine
      integer i, k, kk
      
$if defined (DEBUG_PRINT)
         write(*,*)'Phase=',Phase,'   ParamDeclRealmatrix ', Name, &
                   ',  ModelParamError:', ModelParamError
$endif

		if (ModelParamError) return
		call Check_dimParamRef (Name)		
		
		if (Phase == Building) then

			ParamType(ParamCount) = pTypeReal
			ParamDim (ParamCount) = 2
			do k=1,2
            Param_LBound(k,ParamCount) = lbound(paramArray,k)
			   Param_UBound(k,ParamCount) = ubound(paramArray,k)
         enddo
			paramRef (ParamCount)%RealMatrix => paramArray

            ! Fill each row of ParamArray from corresponding row in InitArray
            ! Padding each row with value in last column of InitArray,
            ! and repeating last row of Initarray if needed.

         N = shape(paramArray)
         LastInitColumn = 	ubound(InitArray,1)
         FirstRow = lbound(ParamArray,2)
         do k = FirstRow, ubound(ParamArray,2)
            kk = min(ubound(InitArray,2),lbound(InitArray,2) + (k-FirstRow))
		      ParamArray(:,k) = reshape( InitArray(:,kk), (/N(1)/), &
                                       (/(InitArray(LastInitColumn,kk),i=1,N(1))/))
         end do

		elseif (Phase == Writing) then

	      FirstColumn = lbound(paramArray,1)
	      LastColumn  = ubound(Paramarray,1)
         LastInLine  = min(FirstColumn + 4,LastColumn)
         		   
         do k = lbound(ParamArray,2), ubound(ParamArray,2)
			   write(ParamFileUnit,'(2A,3(I3,A),5(1X,G13.6))', IOSTAT = ParamOutputError) &
	               trim(Name), "(", FirstColumn, ":", LastColumn, "," , k,") = ", &
                  ParamArray(FirstColumn:LastInLine,k)
			   if (LastColumn > LastInLine) then
	            write(ParamFileUnit,'(1X,5(1X,G13.6))', IOSTAT = ParamOutputError) &
	               ParamArray(LastInLine+1:LastColumn,k)
	         endif
         end do
	      call WriteParamDescr(Name, UnitText, DescrText)

      endif

      end subroutine ParamDeclRealMatrix

!===================================================================      


! -------------------------------------------------------------------		
! Check for program errors end ensure that the ParamRef array is not
! addressed outside dimensioned range:

		subroutine Check_dimParamRef(Name)
		character*(*) Name

		ParamCount = ParamCount + 1
      if (Phase /= Counting .and. ParamCount > dimParamRef) then
			write (*,*) "Error in Parameter Setup, Phase =",Phase, &
                     ": Parameter no ", ParamCount, " Name ", Name, &
                     " Exceeds allocated space ",dimParamREF," of ParamRef" 
			ModelParamError = .true.
		elseif (Phase == Building) then
			call SaveName(Name, paramRef(ParamCount)%NameChar,paramRef(ParamCount)%NameLength)
      endif	

      end subroutine
      
! --------------------------------------------------------------
! Transfer name string into a pointer array of single characters,
! allocated as needed:

		subroutine SaveName(Name, NameChar, Length)
		character, pointer :: NameChar(:)
      character*(*) :: Name

		integer :: First, last, length, i
		
      Last = len_trim(Name)       ! last  nonblank (or 0)
		First  = verify(Name, " ")  ! First nonblank (or 0)
      if (Last > 0) then
      	Length = Last-First+1
			Allocate( NameChar(1:Length), &
                   STAT = Allocate_Status)
	      if (Allocate_Status == 0) then
	      	do i = First, Last
					NameChar(i-First+1) = Name(i:i)
				enddo
			else
				write(*,*)"Allocation of NameChar array failed, "
				ModelParamError = .true.
			endif
		else
			write(*,*)'Blank parameter name; not allowed'
			ModelParamError = .true.
      endif

      if (ModelParamError) then
         write(*,*)"Error for Parameter Number ",ParamCount
         if (ParamCount>1) then 
            write(*,*)'Previous parameter was ', ParamRef(ParamCount-1)%NameChar
         endif
      endif   
		end subroutine SaveName




! ==================================================================
! For ParamDecl subroutines when Writing to file:
! write additional parameter description
! and error message if error condition occur.
! (not called if error already has occurred during the write sequence)

      subroutine WriteParamDescr(Name, UnitText, DescrText)
		character*(*) :: Name, UnitText, DescrText

		if (ParamOutputError==0) then
         if (UnitText == " ") then
            write(ParamFileUnit,'(10x,A)', IOSTAT = ParamOutputError) &
                          "# Dimensionless."
         else
            write(ParamFileUnit,'(10x,2A)', IOSTAT = ParamOutputError) &
                          "# Unit: ", UnitText
		   endif
         if (ParamOutputError==0) then
            call WriteParamDescrString(DescrText)
		   endif
		endif
		
      if (ParamOutputError /= 0) then
		   write(*,*)" Error when writing parameter ", Name, &
                   " to file '", ParamFileNameChar, "'", &
                   "; Error code=", ParamOutputError
		   ModelParamError = .true.
      endif


      end subroutine WriteParamDescr


! ------------------------------------------------------------------------
! Write description string to file, with word-wrap:

      subroutine WriteParamDescrString(String)
      character*(*) String
      integer Length, i, k

      Length = len_trim(String)

$if defined (DEBUG_PRINT)
      write(ParamFileUnit, *) "|",String(1:Length),"|"
$endif

      i=0
      do while (i < Length .and. ParamOutputError == 0)

          ! Write either to last blank or hyphen within 70 characters
          ! or first blank or hyphen after 70 characters:

         k = Length-i
      	if (k > 60) then
      	   k = scan(String(i+1:i+60)," -", .true.)            ! search backwards
        	   if (k == 0 ) k = scan(String(i+61:Length)," -")+60 ! search forwards
            if (k == 0 ) k = Length-i			
         endif 

$if defined (DEBUG_PRINT)
      write(ParamFileUnit, *) "i,k=",i,k,"|",String(i+1:i+k),"|"
$endif

         write(ParamFileUnit,'(10x,2A)', IOSTAT = ParamOutputError) &
                  "# ", String(i+1:i+k)

         i = max(i+1,i+k)  ! ensurance against infinite loop

      end do

      write(ParamFileUnit,*)  ! ending with blank line
      
      end subroutine


!----------------------------------------------------------------------

      character*(22) function IntToString(n)
      integer n
      
      integer IOError
      character*22 String
      write(String,"(I12)",IOSTAT = IOError) n
      if (IOError /=0) then
         String = "(IntToString failed!)"
      endif
      IntToString = trim(String)
      end function


!--------------------------------------------------------------------
! Unfinished sketch for trimming output of numerical values

!   subroutine TrimValueString(String, first, last)
!   character*(*) String
!   integer first, last   ! positions to print
!   
!   first = verify (String, " ")
!   last  = verify (String, " ", .true.)
!   delimiter = scan (String(first:last),".E")
!   if (String(delimiter:delimiter)=".") then
!      do while 
!   TrimmedValueString = trim(String)
!   end function

!======================================================================
! for debug: Print saved info about parameter number n:

      subroutine PrintParamInfo(n)
      integer :: pInfUnit = 6
      integer :: n, i
      
	      if (n >= lbound (ParamRef,1) .and. n <= ubound(ParamRef,1)) then

$if defined (DEBUG_PRINT)
	         write(pInfUnit,*)"Parameter ", n, "Type ", ParamType(n), &
	                   " Dimensionality:", ParamDim(n), &
                      " Name and value:"
$endif

            do i=1,ParamRef(n)%NameLength
               write(pInfUnit,'(A)',advance='NO') ParamRef(n)%NameChar(i)
	         enddo

	         select case (ParamType(n))


	         	!----------------------
               case (pTypeLogical)

	         	   select case (ParamDim(n))
	         	      case (0)
			         		if (.not. associated(ParamRef(n)%Logical)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,*)'  = ', ParamRef(n)%Logical
								endif	         	
	         	      case (1)
			         		if (.not. associated(ParamRef(n)%LogicalArray)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,*)'  = ', ParamRef(n)%LogicalArray
								endif
                  end select	         	

         	   !----------------------
         	   case (pTypeInteger)

	         	   select case (ParamDim(n))
	         	      case (0)
			         		if (.not. associated(ParamRef(n)%Integer)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,*)' = ', ParamRef(n)%Integer
								endif	         	
	         	      case (1)
			         		if (.not. associated(ParamRef(n)%IntegerArray)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,'(A)',ADVANCE='NO')'  = '
								 	write(pInfUnit,'(5I6,1x)') ParamRef(n)%IntegerArray
								endif	         	
	         	      case (2)
			         		if (.not. associated(ParamRef(n)%IntegerMatrix)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,'(A)',ADVANCE='NO')' = '
								 	do i=Param_LBound(2,n), Param_LBound(2,n)
                              write(pInfUnit,'("(:,",I6,"):",5I6,1x)') ParamRef(n)%IntegerMatrix(:,i)
                           enddo
								endif	         	
                  end select	         	

         	   !----------------------
					case (pTypeReal)

	         	   select case (ParamDim(n))
	         	      case (0)
			         		if (.not. associated(ParamRef(n)%real)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,'(A,G12.4)')  ' =', ParamRef(n)%real
								endif	         	
	         	      case (1)
			         		if (.not. associated(ParamRef(n)%RealArray)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,*)' = ', ParamRef(n)%RealArray
								endif	         	
	         	      case (2)
			         		if (.not. associated(ParamRef(n)%RealMatrix)) then
									write(pInfUnit,*)'pointer not allocated'
								else
								 	write(pInfUnit,*)' = ', ParamRef(n)%RealMatrix
								endif	         	
                  end select	         	

	         end select

	      endif

      end subroutine

end Module ModelParamToolbox

