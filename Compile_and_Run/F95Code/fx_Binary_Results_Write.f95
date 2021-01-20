      Module Binary_Write

      use ModelVar_HydroBioChem
      use ModelVar_Topography
      use ModelVar_RunControl
      use Binary_Write_sub
      
      implicit none
 
! ==================================================================
! Based on WATCOM Fortran Code - writes data from model run 
! to binary file for later inclusion in other file formats.

! Reorganized for Absoft version into separate subroutines with global private variables
! instead of shared declarations through use of entry points.

! The file is written as a binary file without record length info.
! Compiler-specific details are handled by the module Binary_Write_sub.for
! which contains these functions which are used in the code below:

!       logical BinaryWriteOpenFile
!       logical BinaryWriteNames
!       logical BinaryWriteValues
!       logical BinaryWriteClose
!       integer BinaryWriteErrCode



      integer :: Obs_Number
      logical :: BinaryOpen = .false.
	
      contains

! ======================================================================


$undefine DEBUG_BIN_RES

! ======================================================================
! Called from sub_Snapshot to save results in binary format:


      subroutine STARTBINARYSTORAGE ( BinaryUnit )

      integer BinaryUnit ! unit number for Binary File


      integer Raw/8/          ! ACSL raw data fil number, 
                              ! should not be changed
      character Filename*1024
      logical NameStatus
      integer lfn, extPos, i



$if defined DEBUG_BIN_RES
      write (*,*) ' BinaryUnit:', BinaryUnit
$endif
          ! ... get file name of RRR-file, find and replace extension:
      inquire( unit=Raw, name=FileName, Named=NameStatus)
      if (.not. NameStatus ) FileName='eutro'
      lfn = len_trim( FileName )
      extpos=lfn+1
      if (lfn.ge.4) then
         do i=0, 3
            if ( FileName(lfn-i:lfn-i) .eq. '.' ) then
               extpos=lfn-i
               exit
            end if
         end do
      end if
      FileName = FileName(1:extpos-1)//'.bdf' ! Binary file
      lfn=len_trim(FileName)

$if defined DEBUG_BIN_RES
      write (*,*) ' FileName:', FileName(1:lfn), ' lfn:',lfn
$endif

      call do_Variables (0, 1, 1, Obs_Number )      ! Count variables

$if defined DEBUG_BIN_RES
         write (*,*) 'STARTBINARYSTORAGE debug ', 1
$endif

      call AllocateArrays                           ! Allocate buffer arrays

$if defined DEBUG_BIN_RES
         write (*,*) 'STARTBINARYSTORAGE debug ', 2
$endif

      call do_Variables (1, 1, 1, Obs_Number)       ! Transfer variable names to buffer

$if defined DEBUG_BIN_RES
         write (*,*) 'STARTBINARYSTORAGE debug ', 3
$endif

          ! Open and write start of binary file:

      if(BinaryWriteOpenFile(BinaryUnit, FileName(1:lfn))) then
         call BinaryWriteNames()      ! Write variable names
         BinaryOpen = .true.
 
         Obs_Number=0

$if defined DEBUG_BIN_RES
         write (*,*) ' Binary file opened'
$endif
         
      else
$if defined DEBUG_BIN_RES
         write (*,*) ' Error code ', BinaryWriteErrCode(),
     &               ' on opening binary file'    
$endif
         BinaryOpen = .false.
      endif
      return
      end subroutine


! ======================================================================


      subroutine BINARYRECORDS

!      include 'eutro.inc'
      integer, parameter :: ReportDepthDim = 11
      real*8 ReportDepth(ReportDepthDim)
     &       /1., 5., 10., 20., 30., 40., 50., 60., 80., 100., 130.0/ 
      integer NextDepth
      integer I_B, I_L

$if defined DEBUG_BIN_RES
      integer I_S
      logical Open555
$endif

      if (.not. BinaryOpen) return

      Obs_Number = Obs_Number + 1
                   ! record number written to binary file
      
! ========================================================================
! Commented out by ans@niva.no 16.12.2018
! ========================================================================      
!      do I_B=1,NBI
!         NextDepth=1
!         do I_L = INDXI(I_B)+1, INDXI(I_B+1)
!            if (DEPTH(I_L+1-INDXI(I_B)).ge.ReportDepth(NextDepth)
!     &          .or. I_L.eq.INDXI(I_B+1) ) then
!                  ! Fixed model layer contains ReportDepth,
!                  ! or is the last layer in the basin:
!
!                       ! Transfer variables to buffer
!               call do_Variables (2, I_B, I_L, Obs_Number )
!                       ! Write variables to file
!               call BinaryWriteValues()
!
!               if (NextDepth.lt. ReportDepthDim) then
!                  NextDepth= NextDepth+1
!               else
!                  exit  ! from do I_L =
!               endif
!            endif   
!         enddo
!      enddo

! ========================================================================
! New loop added by ans@niva.no 16.12.2018
! Write every layer
! ========================================================================
      do I_B=1,NBI
         NextDepth=1
         do I_L = INDXI(I_B)+1, INDXI(I_B+1)
                       ! Transfer variables to buffer
               call do_Variables (2, I_B, I_L, Obs_Number )
                       ! Write variables to file
               call BinaryWriteValues()   
         enddo
      enddo
! ========================================================================
! Edits by ans@niva.no finished
! ========================================================================

            ! preliminary solution to documenting dived jets:

$if defined DEBUG_BIN_RES
      do I_S = 1, NS
         if (Ambient_Volume_Flux(I_S).gt.0.0) then
            inquire (555, opened = Open555)
            if (.not. Open555) then
               write(555,'(1x,A,I12)')
     &            'Documenting dived jets in BinaryRecords, NS=', NS
               write(555,'(1x,A10,A3,2A15)')
     &            'YEARS', 'n', 'Amb_Vol_Flux','Neutral_D'
             endif
               write(555,'(1x,F10.6,I3,2G15.6)')
     &            YEARS, I_S, Ambient_Volume_Flux(I_S), Neutral_Depth(I_S)
         endif
      end do      
$endif

      end subroutine


! =====================================================================
! Operate on defined set of model variables according to opCode.
           
      subroutine do_Variables (opCode, I_B, I_L, Obs_Number)
      
      integer opCode 
         ! =0: count values 
         ! =1: store names in name array
         ! =2: store value in name array
         
      integer I_B, I_L        ! Basin nr. and layer number
      integer Obs_Number      ! Index for recording time


      !------------ local variables:
      real*8 rBasinNr, rLayer, rObs_Num, rSesong
      integer I_G
      character*1 G(2) /'1','2'/
      
!      include 'eutro.inc'

      
      call Init_VarCount
      
      rBasinNr = I_B
      rLayer   = I_L - INDXI(I_B)  ! Store local layer nr. within basin
      rObs_Num = Obs_Number
      rSesong  = MOD(YEARS,1.0d0)

! Adjusting the set of variables need only be done here:
! In current version, the variable names should have max. 8 characters;
! as this is the max. limit for Statistica file library STA32.dll
! Variable VNameArray also declared with *8 character elements 
! (stores variable names when writing and reading binary files)

      Call Next_Var ( opCode,  'BASSENG ', rBasinNr     )
      Call Next_Var ( opCode,  'Lag     ', rLayer       )
      Call Next_Var ( opCode,  'ObsNum  ', rObs_Num     )
      Call Next_Var ( opCode,  'T_AAR   ', YEARS        ) 
      Call Next_Var ( opCode,  'SESONG  ', rSesong      ) 
      Call Next_Var ( opCode,  'DYP_M   ', ZMID   (I_L) )
      Call Next_Var ( opCode,  'TEMP    ', TEMP   (I_L) )
      Call Next_Var ( opCode,  'SAL     ', SAL    (I_L) )
      Call Next_Var ( opCode,  'DENSI   ', DENSI  (I_L) )
      Call Next_Var ( opCode,  'OXYG    ', OXYG   (I_L) )
      Call Next_Var ( opCode,  'TOTN    ', TOTN   (I_L) )
      Call Next_Var ( opCode,  'NO3     ', NO3    (I_L) )
      Call Next_Var ( opCode,  'NH4     ', NH4    (I_L) )
      Call Next_Var ( opCode,  'TOTP    ', TOTP   (I_L) )
      Call Next_Var ( opCode,  'PO4     ', PO4    (I_L) )
      Call Next_Var ( opCode,  'SIO2    ', SIO2   (I_L) )
      Call Next_Var ( opCode,  'CZOO    ', CZOO   (I_L) )
      Call Next_Var ( opCode,  'BACT    ', BACT   (I_L) )

      do I_G = 1,2
         Call Next_Var ( opCode, 'CFYT_'//G(I_G),   CFYT (I_L, I_G) )
         Call Next_Var ( opCode, 'NFYT_'//G(I_G),   NFYT (I_L, I_G) )
         Call Next_Var ( opCode, 'PFYT_'//G(I_G),   PFYT (I_L, I_G) )
         if (I_G.eq.1) then
            Call Next_Var ( opCode, 'SFYT' ,   SFYT (I_L ) )
         endif
      enddo
    
      Call Next_Var ( opCode,  'ODM     ', ODM    (I_L) )
      Call Next_Var ( opCode,  'DOC     ', DOC    (I_L) )
      Call Next_Var ( opCode,  'TOTC    ', TOTC   (I_L) )
      Call Next_Var ( opCode,  'PARTC   ', PARTC  (I_L) )
      Call Next_Var ( opCode,  'PARTN   ', PARTN  (I_L) )
      Call Next_Var ( opCode,  'PARTP   ', PARTP  (I_L) )
      Call Next_Var ( opCode,  'CDFLUX  ', CDFLUX (I_L) )
      Call Next_Var ( opCode,  'NDFLUX  ', NDFLUX (I_L) )
      Call Next_Var ( opCode,  'PDFLUX  ', PDFLUX (I_L) )
      Call Next_Var ( opCode,  'SDFLUX  ', SDFLUX (I_L) )
      Call Next_Var ( opCode,  'CSED    ', CSED   (I_L) )
      Call Next_Var ( opCode,  'NSED    ', NSED   (I_L) )
      Call Next_Var ( opCode,  'PSED    ', PSED   (I_L) )
      Call Next_Var ( opCode,  'SSED    ', SSED   (I_L) )
      Call Next_Var ( opCode,  'ASED    ', ASED   (I_L) )
      Call Next_Var ( opCode,  'PADS    ', PADS   (I_L) )

      Call Next_Var ( opCode,  'SPP     ', SPP     (I_L))
      Call Next_Var ( opCode,  'SPPSED  ', SPPSED  (I_L))
      Call Next_Var ( opCode,  'SPPDV   ', SPPDV   (I_L))
      Call Next_Var ( opCode,  'SPPSEDDV', SPPSEDDV(I_L))
      
! ========================================================================
! Added by ans@niva.no 10.5.2019
! ========================================================================       
      do I_G = 1,2
         Call Next_Var ( opCode, 'CHL_'//G(I_G),   CHL (I_L, I_G) )
      enddo
      Call Next_Var ( opCode,  'C1      ', C1      (I_L)) 
      Call Next_Var ( opCode,  'TOTN    ', TOTN    (I_L))
      Call Next_Var ( opCode,  'TOTP    ', TOTP    (I_L))
!      Call Next_Var ( opCode,  'W_T_DERI', W_T_DERIV (1,I_L - INDXI(I_B),I_B))  
!     Call Next_Var ( opCode,  'Basins  ', Basins)
!     Call Next_Var ( opCode,  'Layers  ', Layers)
!     Call Next_Var ( opCode,  'Age_Grou', Age_Groups)
! ========================================================================
! Edits by ans@niva.no finished
! ========================================================================

      end subroutine
      


      end module Binary_Write


