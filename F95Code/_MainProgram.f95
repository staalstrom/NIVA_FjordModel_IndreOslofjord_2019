! NIVA Fjord Model, version 2010, Birger Bjerkeng
! Main program.
! Calls 

Program NIVA_Fjord_Model

   use m1_InitiatedModel
   use m2_ProcessDirectives
   use m3_RunModel
   use m4_ModelState
   use m5_TerminateModel

   implicit none

   integer ActionCommand
   
   integer ProcessingPhase
   integer, parameter :: NoRunStarted=1, Running=2, BetweenRuns=3
   
!----------------------------------------------------------------------   
   if (InitiatedModel()) then

      ProcessingPhase = NoRunStarted

      do
               ! Read commands until return with action command.
               ! Input of model parameters including run control
               ! is handled locally within the ReadCommands subroutine

          call ReadCommands(ActionCommand)  ! in module m2_ProcessDirectives

               ! Perform the action requested by the return argument:
               ! values of constants cmd_... and corresponding input strings
               ! are defined as parameters in m2_ProcessDirectives 

         select case (ActionCommand)

            case (cmd_StartRun)      ! New model run
               if (ProcessingPhase == Running) then
                  call EndRun
               endif
               call StartRun      ! returns on termination or error
               ProcessingPhase = Running

               call TextInp_Debug(.true.) ! .false from start

               call ContinueRun

            case (cmd_EndRun)
               call EndRun
               ProcessingPhase = BetweenRuns

            case (cmd_Continue)   ! continue model run
               call ContinueRun   ! returns on termination or error

            case (cmd_Save)       ! save state variables to file
               call SaveState

            case (cmd_Restore)    ! restore state variables from file
               call RestoreState

            case (cmd_ExitProg)   ! Finalise run and exit processing loop 
               write(*,*)' calling EndRun'
               call EndRun
               ProcessingPhase = BetweenRuns
               exit

         end select

      end do
               write(*,*)' calling TerminateModel'

      call TerminateModel  ! includes any needed final actions

   end if

!  With input from file given in pipe argument (<file ): Pause fails after end of input file
! (Must use write(*,*) Message, and then Read(*,*, IOSTAT=variable) Answer to handle error

   write(*,*)"Program terminated"
end program NIVA_Fjord_Model

