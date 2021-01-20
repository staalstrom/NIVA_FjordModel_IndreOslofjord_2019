Module Code_QA_Procedures

   implicit none

   contains
   
! ====================================================================
! Checks one specific array dimension against size needed by model
! (Consistency of dimensions)
! called by function above, and by subroutines in Module transp_1.for

   logical function CheckDimension(dimName, dimValue, neededvalue )
      character*(*) dimName
      integer dimValue
      integer neededvalue 
      
      if (dimValue .ge. neededValue) then
      	CheckDimension=.true.
      else
      	CheckDimension=.false.
      	write(*,*)'array dimension ',dimName,' = ',dimValue, &
      				' should be >=', neededValue
      endif

   end function

end Module