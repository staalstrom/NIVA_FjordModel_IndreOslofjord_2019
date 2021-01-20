Module ModelDimensions

   implicit none

! Array dimension parameters:
      Integer dimMD       !Max. number of depths >=layers in deepest basin+1
      Integer dimMBI      !Max. number of internal basins
      Integer dimMBIplus1 !MBI+1 (defined as constant due to limitations in ACSL)
      Integer dimMLI      !Max. total number of layers over all internal basins
      Integer dimMLIplus1 !MLI+1
      Integer dimMBE      !Max. number of external basins (boundaries)
      Integer dimMBEplus1 !MBE+1
      Integer dimMLE      !Max. total number of layers over all external basins
      Integer dimMC       !Max. number of connections between basins
      Integer dimMCplus1  !MC+1
      Integer dimMLC      !Max. total number of layers over all connections
      Integer dimMS       !Max number of land sources for water and runoff
      Integer dimMS_sq    !dimMS^2

      Integer dimMFYTG    !Number of different fytoplankton groups
                           !  1: only one group, diatoms and others combined
                           !  2: diatoms and others can be differentiated
                           !     integer constant LFYT decides how many groups
                           !     are actually used in simulation (others =0)
      Integer dimMUSLAGES !Number of age classes in mussels
                           ! (used as 0:MUSLAGES-1 in MUSLDV subr.)
      Integer dimMUSLAYERS!Max. number of depth layers (fixed) 
                           !for mussel population

! Initial constant value parameters for topographical setup
! (see Module TOPO.....FOR for details):
      REAL dimDDPMIN_value    !Initial value of DDPMIN 
                               ! (minimum depth interval)
      REAL dimDPFACT_value    !Initial value of DDPFACT
                               ! (factor for increasing depth interval) 
      Integer dimNLIMAX_value !Initial value of NLIMAX
                               ! (upper limit to number of layers used,
                               !  within dimensioning value)

      Integer dimResidTime    ! Number of defined volumes for calculation
                               ! of residence time of water
! ==============================================================================

! Select different sets of actual value sets
! suited for different fjords:
!       INCLUDE '.\Include\FJORDDIM.INC'

! ===========================================================
! NIVA fjord model
! From file OSLOFJ\EUT_OSLF.INC.
!   gives values to Parameters defined in EUTRODIM.INC 
!   for Inner Oslofjord with 5 basins
      parameter (dimMD             =    41)
      parameter (dimMBI            =     9)
      parameter (dimMBIplus1       =    10)
      parameter (dimMLI            =   180)
      parameter (dimMLIplus1       =   181)
      parameter (dimNLIMAX_value   =   180)
      parameter (dimMBE            =     1)
      parameter (dimMBEplus1       =     2)
      parameter (dimMLE            =    33)
      parameter (dimMC             =    10)
      parameter (dimMCplus1        =    11)
      parameter (dimMLC            =   220)
      parameter (dimMS             =    15)
      parameter (dimMS_sq          = 15*15)
      parameter (dimMFYTG          =     2)
      parameter (dimMUSLAGES       =    10)
      parameter (dimMUSLAYERS      =    10)
      parameter (dimDDPMIN_value   =   1.0)
      parameter (dimDPFACT_value   =   1.1)
      parameter (dimResidTime      =   15 )

! ===========================================================
   
contains

! ==========================================================================
! Checks requirements on dimensions set in EutroDim.inc file,
! and return True/False:

   logical function Dimensions_OK()

	Dimensions_OK = Checkdimension('dimMBIplus1', dimMBIplus1, dimMBI+1) &
	          .and. Checkdimension('dimMLIplus1', dimMLIplus1, dimMLI+1) &
	          .and. Checkdimension('dimMBEplus1', dimMBEplus1, dimMBE+1) &
	          .and. Checkdimension('dimMCplus1' , dimMCplus1 , dimMC+1 )   

   end function

   
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
