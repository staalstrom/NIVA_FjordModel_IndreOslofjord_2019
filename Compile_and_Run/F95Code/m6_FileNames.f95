! Revised by Birger Bjerkeng 2019-10-19

      Module m6_FileNames      

! Defines file name variables and fills them with default values.

! If no file FILENAMES.TXT exists, it is created with default file names.
!
! If it exists, the file names are read from that file
! for use as source files for topography, boundary and runoff.

! Currently, only topography is read from file, 
! while Boundary and Runoff data are coded
! as datainitiated tables in the model.


      CHARACTER*64     TOPOGRAPHY_FILE, BOUNDARY_FILE, RUNOFF_FILE

               ! contains:

      NAMELIST /OPTIONS/ & 
          TOPOGRAPHY_FILE, BOUNDARY_FILE, RUNOFF_FILE 

      DATA TOPOGRAPHY_FILE /'Topography.DAT'  /, &
           BOUNDARY_FILE   /'BOUNDARY.DAT'/, &
           RUNOFF_FILE     /'RUNOFF.DAT'  /   

   contains

! =============================================================================
! Read file names from master file Filenames.Txt using a namelist input format,
! Or, if not found: create the master file with default file names:

      SUBROUTINE Initate_FileNames

      
      CHARACTER*30 FILENAME /'Filenames.Txt'/
      LOGICAL exists
      INTEGER fstatus
      

      

! Open file with data file names, or create file with default names  
      INQUIRE ( FILE = FILENAME, EXIST = exists)
      OPEN ( 101, FILE = FILENAME, IOSTAT = fstatus )
      
      write(*,*) ' file ',FILENAME,': exists=',exists, &
                 '  fstatus=',fstatus
      
      if ( .not.exists ) then
         if ( fstatus .eq.0 ) then
           WRITE  (101,OPTIONS,IOSTAT=fstatus)
         endif
      else
         if ( fstatus.eq.0) then
            WRITE(*,'(/A)') ' reads input file names from file '//FILENAME
            READ(101, OPTIONS, IOSTAT=fstatus)
         endif
      endif   
      if (fstatus.ne.0) then   
         write(*,*) 'Error' , fstatus, &
           'in read/write of file '//FILENAME
      endif

      CLOSE(101)

      
      END subroutine

      end Module
