      Module Binary_Convert_Statistica
      
      use Binary_Convert_Sub_Library
     
		implicit none

          ! STADEV functions as defined in the C++ header file:
          ! ===================================================

! HSTAFILE API_EXP StaOpenFile (LPCSTR szFileName);
! HSTAFILE API_EXP StaCreateFile (short NVars, long NCases, LPCSTR szFileName); 
! HRES API_EXP StaCloseFile (HSTAFILE hSF);
! short API_EXP StaGetNVars (HSTAFILE hSF);
! long API_EXP StaGetNCases (HSTAFILE hSF);
! HRES API_EXP StaAddVars (HSTAFILE hSF, short After, short HowMany);
! HRES API_EXP StaAddCases (HSTAFILE hSF, long After, long HowMany);
! HRES API_EXP StaSetVarName (HSTAFILE hSF, short Var, LPCSTR szName);
! HRES API_EXP StaGetVarName (HSTAFILE hSF, short Var, LPSTR szName);

! HRES API_EXP StaGetNumVarLabels (HSTAFILE hSF, short Var, short FAR * lpNumLabels);
! HRES API_EXP StaGetVarLabelByIndex (HSTAFILE hSF, short Var, short Index,
!                double FAR * lpValue, LPSTR szLabel, LPSTR szLongLabel);
! HRES API_EXP StaGetValueForLabel (HSTAFILE hSF, short Var, LPCSTR szLabel,
!                double FAR * lpValue);
! HRES API_EXP StaGetLongLabelForValue (HSTAFILE hSF, short Var, double Value,
!                LPSTR szLLabel, short BL);
! HRES API_EXP StaAddLabel (HSTAFILE hSF, short Var, double Value,
!                LPCSTR szLabel, LPCSTR szLongLabel);

! HRES StaSetVarFormat (HSTAFILE hSF, short Var, short width,
!                       short dec, short categ, short display);
! HRES API_EXP StaSetData (HSTAFILE hSF, short Var, long Case, double Value);


          ! Equivalent FORTRAN types and specifications

      ! API_EXP : Attribute stdcall
      ! HSTAFILE = long  = Integer*4
      ! HRES     = short = Integer*2 (1 OK, 0 failure)
      ! LPCSTR           = pointer to constant null-terminated string
      ! LPSTR            = pointer to null-terminated string
      ! double           = real*8
                                     ! non-pointer Arguments passed as value


          ! STADEV functions declared for FORTRAN:
          ! ======================================

      integer*4, external, stdcall :: StaOpenFile, StaCreateFile
      integer*2, external, stdcall :: StaCloseFile
      integer*2, external, stdcall :: StaGetNVars
      integer*4, external, stdcall :: StaGetNCases
      integer*2, external, stdcall :: StaAddVars, StaAddCases 
      integer*2, external, stdcall :: StaSetVarName, StaGetVarName
      integer*2, external, stdcall :: StaGetNumVarLabels
      integer*2, external, stdcall :: StaGetVarLabelByIndex
      integer*2, external, stdcall :: StaGetValueForLabel
      integer*2, external, stdcall :: StaGetLongLabelForValue
      integer*2, external, stdcall :: StaAddLabel 

      integer*2, external, stdcall :: StaSetVarFormat
      integer*2, external, stdcall :: StaSetData

      dll_import StaOpenFile        , StaCreateFile, StaCloseFile , &
                 StaGetNVars        , StaGetNCases           , &
                 StaAddVars         , StaAddCases            , &
                 StaSetVarName      , StaGetVarName          , &
                 StaGetNumVarLabels , StaGetVarLabelByIndex  , &
                 StaGetValueForLabel, StaGetLongLabelForValue, &
                 StaAddLabel            , &
                 StaSetVarFormat    , StaSetData

!StaGetVarMD        , not used, removed

         ! STADEV function return values:
         !=======================================

         ! Statistica File handle

      integer*4 hSF  
         ! Return value of HRES STADEV functions on success:
      integer*2 HRES
      integer*2 Result_OK/1/

      logical StatFileExisted
      integer*4 p_to_String, p_to_String2

		logical FileIsOpen/.false./
!		character, pointer :: Output_Error_Text

      integer*4 RecordNum, CasesBefore, CasesUsed, kCase

		
		    ! Statistica file dimensions:      
      integer*2 NVars
      integer*4 NCases

          ! Buffer for Statistica variable names:
      character*9, allocatable :: StaVarName(:)
      integer*2  , allocatable :: StaVarUsed(:)
      

      integer*2 dimStaVar   ! Dimension of allocated arrays
                            ! for Statistica variable info
      
      integer*2  , allocatable :: varNum     (:)
      
      integer*4 RunIdValue


	   integer*4, private :: vFirst, vLast, vCount
      
      
      logical, parameter, private :: DEBUG = .false.
      
      logical, private :: Error
      real*8 RunId_Value, Old_LabelValue, MaxLabelValue



	contains


! =======================================================================================
	   Logical function OpenStatisticaFile(StatisticaFileName, vNameArray, &
					 Runid_Specified_Value, RunId_ShortLabel, RunId_LongLabel)

		character*(*), intent(inout) :: StatisticaFileName

		character*(*) :: vNameArray(:) ! Array of variable names
	
	
						! Additional specifications connected to first vNameArray element,
                  ! which is a Run identification given as command line arguments.
						! May be edited by called subroutine below.
	   integer*4, intent (inout) :: Runid_Specified_Value ! Suggested numerical value
                                             ! (may be changed for Statistica files,
                                             !  see Binary_Convert_Statistica)

	   character*9, intent (inout)  :: RunId_ShortLabel ! short text identification
	           ! 8 bytes + ending null byte for Statistica

		character*41, intent (inout) :: RunId_LongLabel  ! longer description
	           ! 40 bytes + ending null byte for Statistica

      character*9  Old_ShortLabel/' '/
      character*41 Old_LongLabel/' '/

			
	   integer*4 FirstUnref, LastUnref, ToBeAdded
	   integer*4 iVar, kVar

		integer*4 iStatVar, nv, i2
		integer*4 ErrStatus		

      integer*2 :: T_width = 16, T_dec = 6, T_categ = 2, T_display=1
      integer*2 NumLabels
		logical NewLabel

!      character LabelBuffer*41
      integer*2 bufferlength/41/
!      integer*4 p_to_Buffer


		OpenStatisticaFile = .false. ! assumed initally, reset to .true. on success
		

			! Note index range of name array:
	   vFirst = lbound(vNameArray,1)  ! RunId variable; the next variable is file time
	   vLast  = ubound(vNameArray,1)
		vCount = vlast-vFirst+1 

      Error=.false.
      do iVar = vFirst, vLast
        call Normalise(vNameArray(iVar))  ! remove null, align left, set UpperCase
        if (len_trim(vNameArray(iVar)).gt.8) then
          write(DiagUnit,*) ' Binary variable name ', vNameArray(iVar), &
				' has more than 8 characters; too long for Statistica 5.5 file format'
          Error=.true.
		  endif
		enddo
		if (Error) return
   

      ! Open or create Statistica file:

		if (FileIsOpen) then
         write(DiagUnit,*) 'Statistica file was opened before - error in program logic?'
         return
      endif   

      INQUIRE (File = StatisticaFileName, exist= StatFileExisted)   
      p_to_String = CString_p (StatisticaFileName)
      if (p_to_string.eq.0) then  ! can this happen?
         write(DiagUnit,*) 'Too long Statistica file path string'
         return
      endif   

      if (StatFileExisted) then

         hSF = StaOpenFile ( %VAL(p_to_String) ) ! passes string address
                         ! %VAL(character pointer) needed to avoid
                         ! the extra hidden arguments giving character lengths 
                         ! that are added by ABSOFT for character arguments
         if (hSF.eq.0) then
            write(DiagUnit,*)  'Statistica file existed, but could not be opened'
            return
         endif

         NVars  = StaGetNVars  ( %VAL(hSF)) ! number of variables in Statistica file before
         dimStaVar = NVars+vCount         ! max number of variables after adding data

         NCases = StaGetNCases ( %VAL(hSF))
         CasesBefore = NCases
         CasesUsed   = NCases

         if (DEBUG) write(DiagUnit,'(A,2(I6,A))') &
             ' Statistica file existed with ', &
               NVars,' variables and ', NCases, ' Cases'

      else
         NVars     = vCount   ! number of Statistica variables needed in new file   
         dimStaVar = NVars
         
         NCases = 1
         CasesBefore = 0
         CasesUsed = 0

         hSF = StaCreateFile (%VAL(NVars), %VAL(NCases), %VAL(p_to_String) )
         if (hSF.eq.0) then
            write(DiagUnit,*)  'Statistica file do not exist, and could not be created'
            return
         else
            if (DEBUG) write(DiagUnit,*) &
                  'Statistica file created with ', &
                   NVars,' variables and ', NCases, ' Cases'
         endif   
      endif

		FileIsOpen = .true. ! From now on, always end program by closing Statistica file

	   
! Set up work arrays for Statistica file:
      Allocate ( varNum(vFirst:vLast), StaVarName(dimStaVar), StaVarUsed(dimStaVar), STAT=ErrStatus )
      if(ErrStatus.eq.0) then
         write(DiagUnit,'(1x,3(A,I6))') ' Allocation OK'
      else
         write(DiagUnit,'(1x,3(A,I6))') ' Allocation Error ', ErrStatus
         return
      endif

! Prepare Statistica file: 

      FirstUnref=1; LastUnref=NVars
      ToBeAdded = vCount
      varNum = -1  ! all variables has as yet undefined variable numbers



      Error = .false.  ! is set to true if preparations fail


         ! Get existing Statistica variable names (are right adjusted in file): 

      if (StatFileExisted) then
         do iVar=1, NVars
            p_to_String = LOC(STaVarName(iVar))

            if (DEBUG) write(DiagUnit,*) ' IVar:', IVar, ' p_to_String:', p_to_String
            HRES = StaGetVarName(%VAL(hSF), %VAL(iVar), %VAL(p_to_String))
            call StaResult
            if (DEBUG) write(DiagUnit,*) ' StaGetVarName: HRES=', HRES, ' Varname=',STaVarName(iVar)

            if (HRES .ne.Result_OK) then
               write(DiagUnit,'(1x,A,I4)')'Error when trying to get name of Statistica Variable ',iVar
               Error=.true.
            else
               if (DEBUG) write(DiagUnit,*)'Variable ',iVar,':',StaVarName(iVar)
            endif

            call Normalise(StaVarName(iVar))  ! remove null, align left, set UpperCase
            StaVarUsed(iVar) = -1
            do kVar = 1, iVar-1
               if (StaVarName(kVar).eq.StaVarName(iVar)) then
                  write(DiagUnit,*) ' Error: Statistica variable nr. ',iVar,' and ',kVar, &
                             ' have same name: ', StaVarName(iVar)
                  Error=.true.
               endif
            end do
         end do
      endif

     if (.not.StatFileExisted) then
        if (DEBUG) write(DiagUnit,*)' StatFile did not exist, all variables are new'
        iStatVar = 0
     else
        if (DEBUG) write(DiagUnit,*)' StatFile existed, check variables'

             ! Search in existing Statistica variables for the names in vNameArray:
        do iVar = vFirst, vLast

           call Normalise(vNameArray(iVar))  ! remove null, align left, set UpperCase

           do kVar = FirstUnref, LastUnref
                 ! each Statistica variable only referenced once, regardless of names:
              if (StaVarUsed(kVar).lt.0) then
                 if (StaVarName(kVar).eq.vNameArray(iVar)) then ! found

                    if (DEBUG) write(DiagUnit,*)' variable ',vNameArray(iVar), &
                                                ' found as no.', kVar,' in Statistica file'

                    varNum(iVar)=kVar
                    StaVarUsed(kVar)=iVar
                    ToBeAdded = ToBeAdded-1 
                    if (kVar.eq.FirstUnref) FirstUnref = kVar+1
                    if (kVar.eq.LastUnref)  LastUnref  = kVar-1
                    exit
                 endif
              endif
           enddo 
        enddo

        if (DEBUG) write(DiagUnit,*)' To be added: ', ToBeAdded       

        if (ToBeAdded.gt.0) then
             ! add space for new variables:
           HRES=StaAddVars( %VAL(hSF), %VAL(NVars), %VAL(ToBeAdded))
           call StaResult
           if (HRES.ne.result_OK) then
              write(DiagUnit,*) ' Could not add enough new variables'
              return
           endif
        endif
			iStatVar = NVars   
     endif

     if (DEBUG) write(DiagUnit,*)' Ready to initate new variables:',&
         'NVars=', NVars,' ToBeAdded:', ToBeAdded,' vcount=', vcount

         ! Set variable names for new variables:

     if (ToBeAdded.gt.0) then
        do iVar = vFirst, vLast
           if (DEBUG) write(DiagUnit,*)' ivar=', ivar,' varNum(iVar)=',varNum(iVar)
           if ( varNum(iVar).gt.0) Cycle
           iStatVar=iStatVar+1
           StaVarName(iStatVar)=vNameArray(iVar)
           p_to_String = CString_p(StaVarName(iStatVar))
           if (p_to_string.eq.0) then
              write(DiagUnit,*) '"',StaVarName(iStatVar),'"', &
              ' is too long Statistica variable name'
              stop
           endif

           if (DEBUG) write(DiagUnit,*)' adds variable ', vNameArray(iVar),' as no.', iStatVar
   
           if (DEBUG) write(DiagUnit,*)' calls StaSetVarName'
           HRES=StaSetVarName(%VAL(hSF), %VAL(iStatVar), %VAL(p_to_String))
           call StaResult
           if (DEBUG) write(DiagUnit,*)' Result:', HRES

           if (HRES.eq.Result_OK) then
              varNum(iVar) = iStatVar
              StaVarUsed(iStatVar)= iVar
                                ! Set time format for variable varNum(0)
                                ! with file modification time:
              if (iVar.eq.0) then
                 if (DEBUG) write(DiagUnit,*)' calls StaSetVarFormat'
                 HRES=StaSetVarFormat ( %VAL(hSF), %VAL(varNum(0)), %VAL(T_width), &
                                        %VAL(T_dec), %VAL(T_categ), %VAL(T_display))
                 call StaResult
                 if (DEBUG) write(DiagUnit,*)' Result:', HRES
              endif                          
           else
              write(DiagUnit,*) ' Error: could not set variable name ',NVars,' to ', StaVarName(NVars)
              Error=.true.
           endif
        end do 
      endif



!========================================================================================
! Establish value with labels for run identification (Variable VarNum(-1)):


      nv = VarNum(vFirst)

              ! Get range of label values in Statistica file:

      MaxLabelValue=0

      if (DEBUG) WRITE(*,*)'RunIdent is variable ',nv, ' in Statistica file'

      if (DEBUG) write(DiagUnit,*)' calls StaGetNumVarLabels '
      HRES = StaGetNumVarLabels (%VAL(hSF),  %VAL(nv), NumLabels )

      call StaResult
      if (DEBUG) write(DiagUnit,*)' Result HRES:', HRES, "NumLabels", NumLabels

      if (NumLabels.le.0) then

         NewLabel = .true.
      else
	      if (HRES.eq.result_OK .and. NumLabels.gt.0) then
	         Old_ShortLabel=' '; Old_ShortLabel(9:9) =CHAR(0) ! ASCII table position 0 as character
	         Old_LongLabel =' '; Old_LongLabel(41:41)=CHAR(0)
	         p_to_String = LOC(Old_ShortLabel)
	         p_to_String2 = LOC(Old_LongLabel)
	         do i2=1, NumLabels
	
	            if (DEBUG) write(DiagUnit,*)' calls StaGetVarLabelByIndex i2= ', i2
	            HRES = StaGetVarLabelByIndex &
	                      ( %VAL(hSF),  %VAL(nv), %VAL(i2), &
	                        Old_LabelValue, %VAL(p_to_string), %VAL(p_to_string2) ) 
	            call StaResult
	            if (DEBUG) write(DiagUnit,*)' Result HRES:', HRES, 'Old_LabelValue',Old_LabelValue
	
	            if (HRES.ne.result_OK) then
	               write(DiagUnit,*) 'Error when calling StaGetVarLabelByIndex for variable number nv= ',nv, &
                             ' for label with index i2=',i2
	               Error=.true. 
	            else
	               MaxLabelValue = max(MaxLabelValue, Old_LabelValue)
	               if (RunId_Specified_Value.ne.-9999) then
	                  if (Old_LabelValue.eq.RunId_Specified_Value ) then
	                  		! label value already used, prepare for using new value instead:
	                     RunId_Specified_Value = -9999
	                  endif
	               endif
	            endif
	         end do
	      endif
	
	      IF (debug) Pause "Has read existing RUNIDENT labels"
	
	              ! Check if new label is in Statistica file already:
	
	      p_to_String = CString_p(RunId_ShortLabel)
	      if (p_to_string.eq.0) then
	         write(DiagUnit,*) '"',RunId_ShortLabel,'"', &
	          ' is too long for Statistica short labels'
	            stop
	      endif
	      
	      if (DEBUG) write(DiagUnit,*) ' calls StaGetValueForLabel  , nv:', nv,' Label: ',RunId_ShortLabel
	      HRES = StaGetValueForLabel (%VAL(hSF), %VAL(nv), %VAL(p_to_String), RunId_Value)
	      call StaResult
	      if (DEBUG) write(DiagUnit,*) ' Result:', HRES
         NewLabel = ( HRES .ne. result_OK)
	   endif

	   
	   if ( NewLabel .and. .not.Error) then

              ! not found, add new label, with specified value or with value=highest +1

         p_to_String  = CString_p(RunId_ShortLabel)
         p_to_String2 = CString_p(RunId_LongLabel)

         if (p_to_string2.eq.0) then
            write(DiagUnit,*) '"',RunId_LongLabel,'"', &
                       ' is too long for Statistica long labels'
            stop
         endif

         if (RunId_Specified_Value.ne.-9999) then
            RunId_Value = RunId_Specified_Value
         else
            RunId_Value = MaxLabelValue +1
         endif    

         if (DEBUG) then
             write(DiagUnit,*)' calls StaAddLabel  , nv:', nv, &
              ' RunIdValue: ',RunId_Value
             write(DiagUnit,*) "RunId_ShortLabel:", RunId_ShortLabel, &
                            ",  RunId_LongLabel:", RunId_LongLabel
         endif
    
         HRES = StaAddLabel (%VAL(hSF),  %VAL(nv), %VAL(RunId_Value), &
                               %VAL(p_to_String), %VAL(p_to_String2))


         call StaResult
         if (HRES.ne.result_OK) then
            write(DiagUnit,*)' Error when calling StaAddLabel  , nv:', nv,' RunIdValue: ',RunId_Value
            Error=.true.
         endif

      elseif ((HRES.ne.result_OK) .and. .not.Error) then

             ! found, check long label:

        Old_LongLabel=' '
        p_to_String2 = LOC(Old_LongLabel)

        if (DEBUG) write(DiagUnit,*)' StaGetLongLabelForValue ', &
                                    ' , nv:', nv,' RunId_Value: ',RunId_Value
        HRES = StaGetLongLabelForValue &
                (%VAL(hSF), %VAL(nv), %VAL(RunId_Value), &
                 %VAL(p_to_String2), %VAL(bufferLength))
                      ! fills string with terminating null
        call StaResult
        if (DEBUG) write(DiagUnit,*)' Result:', HRES
        if (DEBUG) write(DiagUnit,*)' RunIdvalue =', RunId_Value

        if (HRES.eq.result_OK) then
           call Normalise(RunId_LongLabel)
           call Normalise(Old_LongLabel)
           if (RunId_LongLabel .ne. Old_LongLabel) then
              write(DiagUnit,*)' Error: Long labels inconsistent for RunIdent=', &
                             RunId_ShortLabel
              write(DiagUnit,*)' in command line now      : ', RunId_LongLabel
              write(DiagUnit,*)' in Statistica file before: ', Old_LongLabel
              Error = .true.
           endif
        endif
      endif
      
      if (HRES.ne.result_OK) then
         write(DiagUnit,*)' Error during prep. of RunIdent variable'
         Error = .true.
      else
         OpenStatisticaFile = .true.
      endif

		end function OPenStatisticaFile

!===============================================================================

		
		logical function AddRecordsToStatisticaFile(RecordNum, vArray)		
		integer*4 RecordNum
      real*8 vArray(:,:)    ! Assumes first index to have same range as vNameArray
      
      integer*4 newRecords
		integer*4 nRecord, iVar
      
      AddRecordsToStatisticaFile = .false.
      newRecords = RecordNum-(NCases-CasesUsed)
      
		       ! add records to Statistica File to receive new data,
             ! taking into account the inital unused case in a newly created file:

         write(DiagUnit,*)' Adding ', newRecords,' records to Statistica file'

         HRES = StaAddCases(%VAL(hSF),%VAL(NCases),%VAL(newRecords))
         NCases = NCases+newRecords

         call StaResult
         if (DEBUG) write(DiagUnit,*)' Result:', HRES

         if (HRES .ne. Result_OK) then
            write(DiagUnit,*)' Error when adding ',newRecords,' cases after case ', NCases
            Error=.true.
         else
            do nRecord = 1, RecordNum
               kCase = CasesUsed + nRecord ! Case number in Statistica file
         		if (DEBUG) write(DiagUnit,*)' Save data in record ', kCase,' in Statistica file'

               Error = .false.

               do iVar=vFirst, vLast
                  HRES = StaSetData ( %VAL(hSF), %VAL(VarNum(iVar)), &
                              %VAL(kCase), %VAL(vArray(iVar,nRecord)))
                  call StaResult
                  if (HRES .ne. Result_OK) then
                     write(DiagUnit,*)' Error on storing var. nr. ',VarNum(iVar)
                     Error=.true.
                  endif   
               end do

               if (Error) then
                  write(DiagUnit,*) 'Stopped because of error in writing Statistica file'
                  write(DiagUnit,*) ' in record ',kCase
                  return
               else
      				AddRecordsToStatisticaFile = .true.
					endif   

            end do

            CasesUsed = CasesUsed + RecordNum

         endif   
      end function     

! ------------------------------------------------------------------------------

		logical function CloseStatisticaFile()
			if (FileIsOpen) then       
		   	HRES = StaCloseFile (%VAL(hSF))
				if (HRES .ne. result_OK) then
		         write(*,*)' Error on closing Statistica file'
		   	else
					FileIsOpen=.false.
					if (.not.Error) then
			         write(*,*)' Added ', NCases - CasesBefore, ' records'
			         write(*,*)' Statistica file now has ', NCases, ' records'
			      endif
		      endif
	      endif
			CloseStatisticaFile = .not.FileIsOpen
		end function
!===================================================================

         
      subroutine StaResult
      implicit none
!      if (DEBUG) then  
!         write(DiagUnit,*) ' **** HRES=',HRES
!         if (HRES) then 
!            write(DiagUnit,*) ' if(HRES) as logical is true'
!         else 
!            write(DiagUnit,*) ' if(HRES) as logical is false'
!         endif
!         if (.not.HRES) then 
!            write(DiagUnit,*) ' .not.HRES is true'
!         else 
!            write(DiagUnit,*) ' .not.HRES is false'
!         endif
!         if (HRES.eq.Result_OK) then 
!            write(DiagUnit,*) ' HRES .eq. Result_OK'
!         else 
!            write(DiagUnit,*) ' HRES .eq. Result_OK'
!         endif
!     endif
      
      end subroutine StaResult         

  ! Added by ans@niva.no - write to ascii file      
      SUBROUTINE write_ascii(vArray2,RecordNum2)
              real*8 vArray2(:,:)  
              integer UNIT_now
              integer iascii   
              integer RecordNum2  
              UNIT_now = 49   
              iascii=1
              OPEN (UNIT_now, FILE='eutro.csv', STATUS='old',POSITION="append")
              DO while (iascii .lt. RecordNum2)
 !                     write(UNIT_now, '(1x, F, 52(";", F))') vArray2(1:52,iascii)
                      WRITE(UNIT_now, '(f15.6,51(";",f15.6))' )  vArray2(1:52,iascii) 
                      iascii = iascii +1 
              END DO   
              CLOSE (UNIT_now)
      END SUBROUTINE  
  ! ans@niva.no - finished      

      End Module
