Module ModelParamHandling

! Contains declarations, subroutines and functions to:
! 1. Build model parameter description during initiation of model
! 2. Write current complete model parameter set to a text file.
! 3. Read arbitrary subset/order of parameter values from text file

! This Module cooperates with parameter group handlers in Modules 
! ModelParam_xxxxxx.f95 with include files ModelParam_xxxxxx.inc
! The parameter group handlers is called by subroutine
! ModelParam_Sequence, which is defined below.

   use ModelParamToolbox

! Refers Parameter group handler Modules:

	use ModelParam_Topography
	use ModelParam_Boundary
   use ModelParam_InitState
   use ModelParam_Physics
	use ModelParam_Plankton
	use ModelParam_Mussels
	use ModelParam_Decomposition
   use ModelParam_Inputs
   use ModelParam_RunControl
   use sub_TextInput
   
	implicit none


      ! Type of return value from function Parameter_Input
      ! which reads parameter values until interrupted
      ! by top-level command:

   type ModelParamInterrupt
      integer :: Num
      logical :: OK
      logical :: EOF
   end type
 

   contains

! ======================================================================
	subroutine ModelParam_Sequence

! Calls all model parameter group handlers in Modules ModelParam_xxxxxx.
! These handlers in turn calls functions defined in ModelParamToolbox:
!      ParamDeclLogical, ParamDeclInteger, ParamDeclReal.
!      ParamDeclRealArray or ParamDeclRealMatrix
! in order to count, register or write parameter values to a file.

! ModelParam_Sequence is called in ModelParam_Setup and ModelParam_Output 
! defined below, where actions in the ParamDecl.... subroutines 
! are controlled by a variable Phase defined in the ModelParamToolbox

		GroupCount = 0
      ParamCount = 0
      call ParamGroup_Topography
      call ParamGroup_Boundary
      call ParamGroup_InitState
      call ParamGroup_Physics
      call ParamGroup_Plankton
      call ParamGroup_Mussels
      call ParamGroup_Decomposition
      call ParamGroup_Inputs
      call ParamGroup_RunControl

	end subroutine


!=======================================================================
! Setting up parameter descriptions as part of initialisation of model.
! Should only be done once during intitialisation of the program:

	logical function ModelParam_Setup()

      if (Phase <> NotStarted) then 
      	write(*,*) 'Repeated call to ModelParam_Setup; '// &
                    'can only be performed once'
         ModelParam_Setup=.false.
         return           
      endif

				! First counts all parameters as preparation for registering

      Phase = Counting 
		call ModelParam_Sequence
      if (ModelParamError) return 

				! Prepares arrays for parameter register:

		dimGroupDescr = GroupCount+1
      dimParamRef   = ParamCount

		allocate( ParamRef (1:dimParamRef)      , &
                ParamType(1:dimParamRef)      , &
                ParamDim (1:dimParamRef)      , &
                Param_LBound(2,1:dimParamRef) , &
                Param_UBound(2,1:dimParamRef) , &
                pGroupDescr(1:dimGroupDescr)  , &
                stat = Allocate_Status)
      if (Allocate_Status /= 0) then
         write(*,*)' Could not allocate parameter description arrays'
         ModelParamError = .true.
         return
		else
        write(*,*)' Parameter description arrays allocated, dimParamRef=', dimParamRef
      endif

				! Perform registration for later use in parameter updates:		

		Phase = Building

		call ModelParam_Sequence

      if (ModelParamError) return 
		
		      ! Finalise group description array:
      if (GroupCount < dimGroupDescr) then
			pGroupDescr(GroupCount+1)%FirstParam = ParamCount+1
			Phase = Complete
		else
         ModelParamError = .true.
      endif


         	! From now on, no more parameters can be added;
         	! the only possible operations are to
            ! write all parameters to a text file 
            !      (sub ModelParam_Output)
            ! or update selected parameters from a text file
            !      (sub ModelParam_Input) 

      ModelParam_Setup = .not. ModelParamError

   end function

!========================================================================      
! Write description of model parameters to specified file
! The result file can be used as template for parameter input files,
! where values can be set for an arbitrary subset of parameters

   logical function ModelParam_Output(FileUnit, FileName)
      integer FileUnit
		character*(*) Filename
      
      if (Phase == Complete) then
         ParamOutputError = 0
         ParamFileUnit = FileUnit
         open (ParamFileUnit, FILE = FileName, IOSTAT = ParamOutputError)
		   call SaveName(FileName, ParamFileNameChar, ParamFileNameLength)
		   if (ParamOutputError == 0) then
				Phase = Writing

			   call ModelParam_Sequence

			   Phase = Complete
			   if (ParamOutputError == 0) then
               close(ParamFileUnit, IOSTAT = ParamOutputError)
               ParamFileUnit = 0
            endif
			endif

			if (ParamOutputError /=0) then
			   	write(*,*)" Parameter output file could not be closed"
			   ModelParamError = .true.
			endif

			Deallocate(ParamFileNameChar, STAT = Allocate_Status)
	      if (Allocate_Status /= 0) then
	         write(*,*)' Could not deallocate file name character buffer'
	         ModelParamError = .true.			
			endif
			ModelParam_Output = .not. ModelParamError
		else
         ModelParam_Output = .false.
      endif
      
   end function


! ==============================================================
! Read and process free-form input of parameter changes
! until interrupted by one of interrupt commands that are
! listed in character array argument.
! Returns error status or Interrupt number in function value.

   type(ModelParamInterrupt) function ParameterInput(fd, InterruptCommands)
		type (TextInputFileDescriptor) :: fd
      Character*(*) :: InterruptCommands(:)

      integer      :: firstCmd, LastCmd
      Character*40 :: Word
      Character*1  :: CharArray(40)
      integer      :: ArrayIndex(1:2)
      integer i, k, lgtWord
      logical OK, InputError, FoundParameter
      Character*132 :: ErrorMessage

      firstCmd = lbound(InterruptCommands,1)
      lastCmd  = ubound(InterruptCommands,1)
      InputError = .false.


      ReadLines: do while (TextInp_Get_Line(fd))

               ! Read first word in line,
               ! should be either Action Command or Parameter name:

		   if (.not.TextInp_Word (fd,Word)) then
            call TextInp_StandardMessage(fd)
            call TextInp_Message(fd,"Expected command or Parameter name")
            InputError=.true.
            cycle ReadLines
         endif
               ! Check against interrupt commands:

         CheckCommands: do i=firstCmd, LastCmd
            if (equalStr_CaseNeutral(Word,InterruptCommands(i))) then
				   ParameterInput%Num = i
				   ParameterInput%OK  = .true.
				   ParameterInput%EOF = .false.
               return
            endif
         enddo CheckCommands

               ! Check against parameter names:
               ! (stored as single-character arrays in table ParamRef 
               ! defined in module ModelParamToolbox and populated
               ! by the function Param_Setup above:

         lgtWord = len_trim(Word)
         do i=1, lgtWord
            CharArray(i)=Word(i:i)
         enddo

         FoundParameter=.false.

         CheckParameters: do i=1, ParamCount

            if (equalCharArray_CaseNeutral(CharArray(1:lgtWord),ParamRef(i)%NameChar)) then

               FoundParameter=.true.
                     ! Parameter name identified; update parameter value:                              

               if (ParamDim(i) == 0) then
                      ! Simple variable;
                      ! if followed by "=": read "=value" and update parameter,
                  if(TextInp_DelimiterFound(fd) == "=") then
	                  select case (ParamType(i))
	                     case (pTypeLogical)
	                        OK = TextInp_Logical(fd,ParamRef(i)%Logical)
	                     case (pTypeInteger)
	                        OK = TextInp_Integer(fd,ParamRef(i)%Integer)
	                     case (pTypeReal)
	                        OK = TextInp_Value(fd,ParamRef(i)%real)
	                  end select
                     if (.not.OK) then
                        InputError=.true.
                        call TextInp_StandardMessage(fd)   
                     endif
                  endif
                 
               else
                              ! Array or matrix parameter:
                              ! *** May have start indices within parentheses:

                  if(TextInp_DelimiterFound(fd) == "(") then
                     OK = ReadArrayIndices(fd,ParamDim(i),ArrayIndex)

                     if (.not.OK) then
                        call TextInp_StandardMessage(fd)
                        call TextInp_Message( fd, &
                              "Error in starting index for parameter array")
                     else
                              !     Check that indices are within array range
                        CheckIndices: do k=1, Paramdim(i)
                           if (    ArrayIndex(k)<Param_LBound(k,i) &
                               .or.ArrayIndex(k)>Param_UBound(k,i)) then
                              write(ErrorMessage,*) &
                                 "Index no. ",k,"=",ArrayIndex(k), ", is outside limits ", &
                                 "[",Param_LBound(k,1),":",Param_UBound(k,1),"]"
                              call TextInp_Message(fd, ErrorMessage)          
                              OK = .false.
                           endif
                        enddo CheckIndices
                     endif
                     if (OK) then ! Skip blank word implied before end of line or delimiter "=":
                        OK= TextInp_Word (fd, Word)
                        if (OK) then
                           OK = (Word=="")
                        endif
                     endif

                   else ! set indices to first item in array:
                      ArrayIndex = 1  ! (Note: lower bound is always 1 for array arguments)
                      OK = .true.
                   endif  


	                if (.not.OK) then
	                   InputError=.true.
	                   call TextInp_StandardMessage(fd)   
	                else

                        ! After name and optional indices in parenthesis:
                        ! require end of line or assignment operator '=':
                      if (.not. TextInp_EOL(fd)) then
                         if (.not. TextInp_DelimiterFound(fd)=="=") then
                           call TextInp_Message(fd, &
                                   "Expected end of Line or assignment operator '=' after array parameter")
                           OK = .false.
                         endif   

                              ! *** Read values and place in succeeding elements of array:
	                     ReadValues: do while( .not. TextInp_EOL(fd))
	
	                        if(ParamDim(i)==1) then
			                     select case (ParamType(i))
			                        case (pTypeLogical)
			                           OK = TextInp_Logical(fd, ParamRef(i)%LogicalArray(ArrayIndex(1)))
			                        case (pTypeInteger)
			                           OK = TextInp_Integer(fd, ParamRef(i)%IntegerArray(ArrayIndex(1)))
			                        case (pTypeReal)
			                           OK = TextInp_Value(fd, ParamRef(i)%RealArray(ArrayIndex(1)))
			                     end select
	                        else
			                     select case (ParamType(i))
			                        case (pTypeInteger)
			                           OK = TextInp_Integer(fd, &
	                                         ParamRef(i)%IntegerMatrix(ArrayIndex(1),ArrayIndex(2)))
			                        case (pTypeReal)
			                           OK = TextInp_Value(fd, &
	                                         ParamRef(i)%RealMatrix(ArrayIndex(1),ArrayIndex(2)))
			                     end select
	                        endif

	                        if (.not.OK) then
	                           InputError=.true.
	                           call TextInp_StandardMessage(fd)   
	                        endif
	                              !    Move to next element in storage order,
	                              ! skip rest of input if too many elements:
	
	                        ArrayIndex(1)=ArrayIndex(1)+1
	                        if (ArrayIndex(1)>Param_UBound(1,i)) then
	                           if(ParamDim(i)>1) then
	                              ArrayIndex(1) = Param_LBound(1,i)
	                              ArrayIndex(2) = ArrayIndex(2)+1
	                              if (ArrayIndex(2)>Param_UBound(2,i)) then
	                                 exit ReadValues
	                              endif
	                           else
	                              exit ReadValues
	                           endif
	                        endif
	      
	                     enddo ReadValues  

                     endif

                  endif

                  if (.not. TextInp_EOL(fd)) then
                     Inputerror = .true.
                     Call TextInp_Message(fd,"Unused items; too many values in line")
                  endif
               endif
               exit Checkparameters

			   endif

         enddo CheckParameters

         if(.not.FoundParameter) then
            Call TextInp_Message(fd,"Command is not recognized")
         else
            call PrintParamInfo(i)
         endif
            
      end do ReadLines

	  call TextInp_StandardMessage(fd)

  	  ParameterInput%OK = .not. Inputerror
	  ParameterInput%EOF = TextInp_EOF(fd)
	  ParameterInput%Num = i
      
   end function


! ===========================================================================
! Expects input of n integers separated by comma and enclosed in parenthesis.
! Returns .False. for problems with arguments after giving a warning,
! .false. if there is a syntax error, leaving messaging to calling code.

   logical function ReadArrayIndices(fd,n,Indices)
		type (TextInputFileDescriptor) :: fd
		integer n
      integer Indices(:)
      
      integer k
      Character*1 delim

      ReadArrayIndices=.false.

      if ((n<=0).or.(n>ubound(Indices,1)).or.(lbound(Indices,1).ne.1)) then
         call TextInp_Message( fd, &
               "Program error; illegal arguments in call to ReadArrayIndices")
         return
      endif

      if(.not. TextInp_DelimiterFound(fd) == "(") return
      do k = 1, n
         if (.not.TextInp_Integer(fd, Indices(k))) return
         delim = TextInp_DelimiterFound(fd)
         if (k<n) then
            if (delim <> ",") return
         endif
      end do   

      if (delim == ")") ReadArrayIndices = .true.

   end function ReadArrayIndices   


end Module ModelParamHandling

