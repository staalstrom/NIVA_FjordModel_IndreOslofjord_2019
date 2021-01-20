Module m1_InitiatedModel

! Performs initial operations when starting model program:
!    - Checking array dimension variables
!    - Initiating Model parameters
! (does not include initiating state variables; 
!  that is done when a model run is started)

   use ModelDimensions
   use ModelParamHandling, only: ModelParam_Setup, ModelParam_Output
   use ModelVar_RunControl, only: STARTD

   implicit none
   logical InteractiveMode
   character*64 Unit5Name
   
   contains


   logical function InitiatedModel()
   
   logical*4,external :: isatty   ! in Unix support library
   
      WRITE(*,'(///1x,60(1H=))')
      WRITE(*,*)'      NIVA  EUTROPHICATION MODEL, VERSION 2019'
      WRITE(*,'(1x,60(1H=)///)')

      
      if (ModelParam_Setup()) then
         if(ModelParam_Output(11,"Parameter_Output.txt")) then
	         InitiatedModel = Dimensions_OK()
         else
            write(*,*)"Modelparameters set up, but parameter output failed"
            InitiatedModel = .false.
         endif
      endif
      
      InteractiveMode = isatty(5)
           ! returns .true. if a terminal device is connected to the specified unit number.
           ! Trial code: does not work, since it returns true
           ! also when model is started in Batch job with input redirected from file

      write(*,*)'InteractiveMode=',InteractiveMode
      
      inquire(unit=5,name=Unit5Name)
      write(*,*)'Unit 5 has name "',Unit5Name,'"'
      
      STARTD = .false.    ! run not started

              
   end function InitiatedModel


end Module m1_InitiatedModel
