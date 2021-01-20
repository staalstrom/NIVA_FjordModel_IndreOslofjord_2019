Module ModelParam_RunControl

  use ModelParamToolbox

  implicit none 

  ! User Parameters for run control


      logical REPORT
              ! Controls report of mass budget

      logical INITC
              ! Controls initiating of hydrographical profiles
              ! done in subroutine ZHYDR if .True.

      LOGICAL TRACE
              ! controls tracing execution 
              ! by calls to subroutine HELLO

      integer MBPRT(6)
              ! controls debug dump of mass balance calculation
              ! for 6 variables: (1): salt, (2): heat, (3): Oxygen,
              !                  (4): N,    (5): P   , (6): Si
              ! <0: No warning of deviation
              ! =0: warns about deviation, stops after 200 occurrences
              ! >0: also gives detailed report for each layer


      real*4    MBINTV
              ! Approx. time interval 
              ! between mass balance calculations

      LOGICAL MBRSET
              ! Triggers reset of mass balance control

      real*4    ACCUR
              ! controls warning message about mass/heat balance

      LOGICAL VPRT
              ! Activates test print of volume balance

      LOGICAL EXTEST
              ! Activates test print of calculation of
              ! external concentrations ...EX

      LOGICAL DGTEST
              ! Activates test print of degradation processes

      logical MSTEST
            ! turns debug print on/off (if compiled for it)


      real*4    TTRIG
              ! General switch for debugging output:
              ! other switches are active only for T>=TTRIG

      INTEGER ACTION
              ! Switch to perform some action
              !  =1: reduce salinity in basin 1

      LOGICAL BIOOFF
              ! Can be used to turn off biological part
              ! of model, but only after the first step
              
		integer dbgdev
					! Debug device number

      integer DEBUG_STEPS
        ! number of steps remaining to give active debug prints for:
        ! can be reset by user

      Logical LOG_ALL_STEPS
         		! Controls calls to LOGD during simulation


      real*4 STPLIM
			! Max limit to integrating time step


      INTEGER NPRINT
         ! Number of time steps between screen message
         ! about simulation time from derivative section

      real*4 TSTART
         ! TSTART can be reset to start simulation at different seasons

! Communication interval
      real*4 DTTERM

! Simulation interval (from TSTART)
      
      real*4 CIMIN, CIMAX, CITFAC, CIBASE
          ! Controls communication interval.
          !  NEVER BEYOND CIMAX OR T = DTTERM+TSTART,
          !  ABOVE CIMIN IF POSSIBLE WITHIN LIMITATION ABOVE


      LOGICAL TRTEST
          ! .TRUE. triggers debug printouts from TRANSP
          ! and TRNADJ and MTRANS if Modules TRANSP_x.FOR
          ! are compiled with TEST_MODE on

      LOGICAL MDEBUG(7)
         !    index 1..7 transferred into MTRANS
         !    as second last input argument.
         !    Also applies to ADJ_CONC and MIX_CONC in TRCALC.FOR


      LOGICAL DBGRNF
			! controls debug printout of runoff
			
      INTEGER PPTEST
         ! Debug printout for primary production 
         ! and sinking of organic matter
         ! restricted; only down to layer PPTEST (from surface)

      LOGICAL PHYTON

      LOGICAL MXTEST
       !  CONTROLS DEBUG OUTPUT FROM SURFBF AND SURFMX

    ! ########### Storage of simulation results to binary file
    !             calling subroutines in Module Bin_Res.for
    
      logical BinReset
      integer BinFile


      INTEGER NSEED ! should be in separate time series Module?


      LOGICAL REINTG
             ! True: forces reinitialization of integrals and
             !       restarts basis for calculation of means in SPARE



   contains

      subroutine ParamGroup_RunControl

   
      call ParamDeclLogical(REPORT,.false., "REPORT", &
               "Controls report of mass budget in output")

      call ParamDeclLogical(INITC,.TRUE., "INITC", &
               "Controls initiation of hydrographical profiles;" // &
               " is done in subroutine ZHYDR if .True.")

      call ParamDeclLogical(TRACE, .FALSE., "TRACE", & 
               "Turns on/off execution tracing by" // &
               " calls to subroutine HELLO")

      call ParamDeclIntegerArray(MBPRT,(/0/),"MBPRT","days", &
               "Controls debug dump of mass balance" // &
               " calculation for 6 state variables;"//&
               " (1):salt, (2):heat, (3):oxygen, (4):N, (5):P, (6):Si" // &
               " for each variable: if value <0: No warning of deviation."//&
               " =0: warns about deviation, stops after 200 occurrences"//&
               " >0: also gives detailed report for each layer")

      call ParamDeclReal(MBINTV, 0.0,"MBINTV","days", &
               "Approximate time interval" // &
               " between mass balance calculations")

      call ParamDeclLogical( MBRSET, .FALSE.,"MBRSET", &
               "Triggers reset of mass balance control")

      call ParamDeclReal(ACCUR,2.0E-5,"ACCUR","", &
               "Required mass mass/heat balance relative accuracy"// &
               " (exceedance triggers warning if spec. with MBPRT")
               

      call ParamDeclLogical(VPRT,.FALSE.,"VPRT", &
               "Activates test print of volume balance")

      call ParamDeclLogical(EXTEST,.FALSE.,"EXTEST", &
              "Activates test print of calculation of"// &
              " external concentrations ...EX")

      call ParamDeclLogical(DGTEST,.FALSE.,"DGTEST", &
              "Activates test print of degradation processes")
      
      call ParamDeclLogical(MSTEST, .false.,"MSTEST", &
               "Turns on/off debug printing of mussel calculations"// &
               " if mussel code is compiled for debug printing")

      call ParamDeclReal(TTRIG,1.0E30,"TTRIG","", &
              "General switch for debugging output:" //&
              " other switches are active only for T>=TTRIG")

      call ParamDeclInteger(ACTION, 0,"ACTION","", &
              "Switch to perform miscellanuous actions."// &
              "  0: no actions;" // &
              " =1: reduce salinity in basin 1")

      call ParamDeclLogical(BIOOFF, .FALSE., "BIOOFF", &
              "Can be used to turn off biological part" //&
              " of model, but only after the first step" //&
              " has been initiated")

      call ParamDeclInteger(dbgdev, 999,"DBGDEV","", &
				 "Debug device number")

      call ParamDeclInteger(DEBUG_STEPS, 0, "DEBUG_STEPS","", &
             "number of steps remaining to give active debug prints for." // &
            " can be reset by user before starting/resuming a model run")


      call ParamDeclLogical(LOG_ALL_STEPS, .False., "LOG_ALL_STEPS", &
         	"Controls calls to subroutine LOGD during simulation")

      call ParamDeclReal(STPLIM,1.0, "STPLIM","", &
			   "Maximum limit to integrating time step")

      call ParamDeclInteger(NPRINT,10,"NPRINT","",&
            "Number of time steps between each "// &
            "progress monitoring message to the screen")

      call ParamDeclReal(TSTART, 0.0,"TSTART","", &
            "Start time for the simulation (Days) from beginning of year;"// &
            " can be reset to start simulation at different points in time"// &
            " in relation to seasons and beginning of input time series" )

      call ParamDeclReal(DTTERM, 1.0,"DTTERM","", &
            "Simulation interval (from TSTART)")

      call Paramexpl("Controls communication interval:"// &
            " (for model state output time series):")

      call ParamDeclReal(CIMIN, 1.0, "CIMIN","", &
            "Lower limit for communication interval")

      call ParamDeclReal(CIMAX, 7.0, "CIMAX","", &
            "Upper limit for communication interval")

      call ParamDeclReal(CITFAC,0.1, "CITFAC","", &
            "Minimum as fraction of simulation time span so far"//&
            "if between lower and upper limits")

      call ParamDeclReal(CIBASE,0.5,"CIBASE","", &
            "Phase of commuication time, as fraction of day:"// &
            " 0.5 will give results at 12:00 noon") 


      call ParamDeclLogical(TRTEST, .FALSE., "TRTEST", &
            ".TRUE. triggers debug printouts from" // &
            " water transport calculations in subroutines" // &
            " TRANSP, TRNADJ and MTRANS if Modules TRANSP_x.FOR" //&
            " have been compiled with TEST_MODE on")


      call ParamDeclLogicalArray(MDEBUG,(/.FALSE./),"MDEBUG", &
            "Controls debug print of mass transport calculations" // &
            "    1. unity/residence time variable C1," // &
            "    2. Sal, Temp," // &
            "    3. Oxygen   ," // &
            "    4. Nutrients," // &
            "    5. Biological variables, except (6) and (7)" // &
            "    6. Zooplankton" // &
            "    7. Inorganic particles")


      call ParamDeclLogical(DBGRNF, .FALSE.,"DBGRNF", &
			   "Controls debug print of land runoff and discharges")

			
      call ParamDeclInteger( PPTEST, 0, "PPTEST", "", &
            "Debug printout for primary production" // & 
            " and sinking of organic matter" // &
            " restricted; only down to layer PPTEST (from surface)")

      call ParamDeclLogical(MXTEST, .FALSE.,"MXTEST", &
            "CONTROLS DEBUG OUTPUT FROM SURFBF AND SURFMX")

      call ParamExpl("Storage of model snapshots to binary file"//&
                     " calling subroutines in Module Bin_Res.for")
    
      call ParamDeclLogical(BinReset, .true.,"BinReset", &
                 "Set it to .true. to (re)start storage of binary snapshot data")

      call ParamDeclInteger(BinFile, 888,"BinFile","", &
             "Unit number for binary output of results") 

      call ParamDeclInteger(NSEED, 0,"NSEED","", &
             "Specification of random seed:"//&
             ">0: use specified value;"//&
             "=0: create seed from start date and time;"//&
             "<0: read from file (5555555 if not found)")


      call ParamDeclLogical(REINTG,.false.,"REINTG", &
             "Set .true. before resuming a model run"//&
             " to trigger reinitialization of integrals"//&
             " from the current point in time." //&
             " (will then be reset to .false.)")

      
      !  ---------------------- BIOLOGICAL PROCESSES ----------------------


 
   end subroutine

end Module
