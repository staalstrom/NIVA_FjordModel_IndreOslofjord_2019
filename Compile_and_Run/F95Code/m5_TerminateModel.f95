Module m5_TerminateModel
   use fx_Eutrosub
   implicit none

   contains
   
   subroutine TerminateModel
   
   write (*,*)"Terminates model session and exits program"

   CALL HELLO ('m5_TerminateModel: CLOSING DEBUG UNIT ')
   CALL DEBUGF(-1)  ! close debug unit

   end Subroutine

end Module