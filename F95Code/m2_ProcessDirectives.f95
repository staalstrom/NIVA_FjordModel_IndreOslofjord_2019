Module m2_ProcessDirectives

! Reads command input from keyboard or files, consisting of:
! - Switching to command files (hierarcy of up to 10 levels)
! - Model Parameter changes
! - Action commands controlling model runs

! The subroutine ReadCommands is called from the main program,
! and processes directives until an action command is encountered.
! It then returns with the command number as return argument.

! Command file changes are done directly in this module.

! Model parameter changes are done in module ModelParamHandling.



   use ModelParamHandling
   use sub_TextInput
   
	implicit none

	      ! Action command values returned to main program:

	   integer, parameter :: &
               NoCommand       =0, &
               cmd_CommandFile = 1, & 
               cmd_StartRun    = 2, &
               cmd_Continue    = 3, &
               cmd_EndRun      = 4, &
               cmd_Save        = 5, &
               cmd_Restore     = 6, &
               cmd_ExitProg    = 7

         ! Strings recognized as Interrupt commands from input files:
	      ! Note that index starts at 1, so that Action command value
         ! and array index will match in subroutine ReadCommands
         ! (Subroutine will handle array arguments as having 
         !  lower bound 1 independently of the bounds in actual argument)

      character*12, parameter :: &
               CommandList(1:cmd_ExitProg) &
                  = (/"CommandFile ", &
                      "Start       ", &
                      "Continue    ", &
                      "End         ", &
                      "Save        ", &
                      "Restore     ", &
                      "Exit        "/)


	      ! Command file control variables

		integer, parameter :: StdInputUnit = 5   , &
									 FirstCmdFile = 1000, &
                            LastCmdFile  = 1010          	                        


         ! File descriptor array dimensioned for 10 levels of command files
         ! (Array not used as argument, so it is OK with lower bound 1000):
         
      type(TextInputFileDescriptor), private, target :: &
            StandardCommandFile, CommandFile(FirstCmdFile:LastCmdFile)
      type(TextInputFileDescriptor), private, pointer :: cmdfile

      integer, private :: cmdFileNum

		integer, private :: NewFileNumber, FileStatus

		integer, private :: CommandFileNumber = StdInputUnit
      
      integer, private :: LastCommand = NoCommand  ! program not initiated
      logical, private :: NewSession = .true.
      
      character Answer

      
   contains

! =====================================================================
! Controls input from command files.
! Initially prepares reading from Standrad input (unit 5)
! and calls function ParamInput in module ModelParamHandling,
! which reads and processes Model Parameter input,
! and returns at end of file or when one of the Command strings above
! are found.

      subroutine ReadCommands(Command)
      integer :: Command


      type(ModelParamInterrupt) mpInterrupt
      character*132 :: FileName

      if (NewSession) then   ! initiate reading from Standard unit:

         cmdFileNum = FirstCmdFile-1
         cmdFile => StandardCommandFile
         if (.not.TextInp_Openfile(cmdFile,StdInputUnit, &
               "Command file" , "Standard input unit", "()=,")) then
            call TextInp_StandardMessage(cmdFile)
            Command = cmd_ExitProg
            return 
         endif
         NewSession=.false.

      endif



      ReadInput: do

               ! Read parameter input until error, EOF or Command:

        mpInterrupt = ParameterInput(cmdFile,CommandList)

                  ! Input error: terminate program

         if (.not.mpInterrupt%OK) then
            Command = cmd_ExitProg
            return
         endif

                  ! End of file: go back to previous command file level

         if (mpInterrupt%EOF) then
         
            write(*,*) 'End of file "', trim(TextInp_GetFileName(cmdfile)),'", Unit=', TextInp_GetFileUnit(cmdfile)
            
            if (cmdFileNum >= FirstCmdFile) then
               write(*,*) 'Closing input file'
               call TextInp_CloseFile(cmdFile)
               cmdFileNum = cmdFileNum-1
               if (cmdFileNum >= FirstCmdFile) then
                  cmdFile => CommandFile(cmdFileNum)
               else
                  cmdFile => StandardCommandFile
               endif
               write(*,*) 'continues with file ', trim(TextInp_GetFileName(cmdfile)),' Unit=', TextInp_GetFileUnit(cmdfile)
 
               cycle ReadInput
            else ! already reading from standard unit
               write(*,*)'already standard unit, ses Command = cmd_ExitProg'
               Command = cmd_ExitProg
               return
            endif

         endif

                  ! New Command file: 
                  ! get name and open as next command file level 	      

         if (mpInterrupt%Num == cmd_CommandFile) then

	         if (cmdFileNum >= LastCmdFile) then
               call TextInp_Message(cmdFile, &
                        "Exceeded maximum command file hierarchy level")

	         else if(.not.TextInp_Word(cmdFile,FileName)) then
                  call TextInp_Message(cmdFile,&
                         "No file name specified")
	         else
	            cmdFile => CommandFile(cmdFilenum+1)
		         if (.not.TextInp_Openfile(cmdFile,cmdFileNum+1, &
		                   "Command file" , FileName, "()=,")) then
		            call TextInp_StandardMessage(cmdFile)

		         else      ! Success:
                  cmdFileNum = cmdFileNum+1
                  cycle ReadInput
               endif            
            endif        ! Failure:
            call TextInp_Message(cmdFile,"CommandFile ignored")
            
          else 

                  ! Other Action Command: return to calling program

            Command = mpInterrupt%Num
            return

         endif

      end do ReadInput


	end subroutine ReadCommands


end Module
