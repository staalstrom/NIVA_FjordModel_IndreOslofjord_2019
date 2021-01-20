Module Binary_Convert_Branching

use Binary_Convert_Statistica

implicit none

		! Transfers general calls to outputfile from main program
		! to the specific functions for the selected output format.
		! This version has only calls to Statistica 5.5 file output,
		! but the code can be developed into case select structures
		! with the switch set for instance by a new command line argument
		! or by the output file extention.

contains

	logical function OpenOutputFile(FileName, vNameArray, &
					 Runid_Specified_Value, RunId_ShortLabel, RunId_LongLabel)
		character*(*), intent(inout) :: FileName      ! Path name for outputfile to be opened or created
		character*(*) :: vNameArray(:) ! Array of variable names
	
	
						! Additional specifications connected to first vNameArray element,
                  ! which is a Run identification given as command line arguments.
						! May be edited inside subroutine called below.

	   integer*4, intent (inout) :: Runid_Specified_Value ! Suggested numerical value
                                             ! (may be changed for Statistica files,
                                             !  see Binary_Convert_Statistica)

	   character*9, intent (inout)  :: RunId_ShortLabel ! short text identification
	           ! 8 bytes + ending null byte for Statistica

		character*41, intent (inout) :: RunId_LongLabel  ! longer description
	           ! 40 bytes + ending null byte for Statistica
											 ! 
		
		OpenOutputFile = OpenStatisticaFile(FileName, vNameArray, &
					 Runid_Specified_Value, RunId_ShortLabel, RunId_LongLabel)
	
	end function
	
	logical function AddRecordsToOutputFile(RecordNum, vArray)		
		integer*4 RecordNum
      real*8 vArray(:,:)    ! Assumes first index to have same range as vNameArray

		AddRecordsToOutputFile = AddRecordsToStatisticaFile(RecordNum, vArray)		

	end function

	logical function CloseOutputFile()
	
		CloseOutputFile = CloseStatisticaFile()

	end function
end module

