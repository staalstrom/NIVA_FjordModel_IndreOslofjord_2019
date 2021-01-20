Module BinaryFile_ControlVariables

! Contains variable and constant declarations, common to writing and reading binary result files


      character(LEN=*), parameter :: BinaryFileHeader =    "NIVA fjord model binary result file "
                    ! Written as start of binary files iwhen storing data from model,
                    ! used to check file version when reading binary data.
                    ! Length should be multiple of 4 bytes to make hexdump inspection easier.

      integer*4 :: BinaryFormatVersion 
                    ! Set to current version when writing data,
                    ! and according to file content when reading data for conversion
                    !                        
      integer*4 :: BinaryUnit, BinFileErrorStatus

end module 
