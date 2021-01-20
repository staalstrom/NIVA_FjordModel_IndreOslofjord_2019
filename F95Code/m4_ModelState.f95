Module m4_ModelState

! Code Skeleton - Unfinished

	   use ModelVar_HydroBioChem
	   use ModelVar_Topography
	
	   implicit none


      integer, parameter  ::  ModelState_Save     = 1, &
                              ModelState_Restore  = 2, &
                              ModelState_Printout = 3

            ! Status variables set internally
            ! and returned through the main external call:
            
      integer, private :: CurrentAction, FileUnit
      integer, private :: IO_Status
      logical, private :: Action_Aborted

      integer, parameter :: IO_Error=1, Name_Error=2 
      
   contains

!=====================================================================
   
	   subroutine SaveState
	   
	   end subroutine
	   
	
!=====================================================================

	   subroutine RestoreState
	   
	   end subroutine

 
!=====================================================================

      logical Function ModelState_Action( Action_Spec, FileName)
      
      ! Performs specified action on all Model State variables,
      ! including Topography description:
      
      integer :: Action_Spec
      character*(*) FileName


      !----------------------------------------------------
      
     
      CurrentAction = Action_Spec
      
      ! ------------------ Topography specification:
            
            ! ...... Dimensions:
      
      call do_IntegerDim('ND' ,ND )
      call do_IntegerDim('NBI',NBI)
      call do_IntegerDim('NBE',NBE)
      call do_IntegerDim('NLI',NLI)
      call do_IntegerDim('NLE',NLE)
      call do_IntegerDim('NC' ,NC )
      call do_IntegerDim('NLC',NLC)


            ! ...... Index limit arrays for basins and connections:
      
      call do_IntegerArray('INDXI',INDXI(1:NBI+1))
      call do_IntegerArray('INDXE',INDXE(1:NBE+1))
      call do_IntegerArray('INDXC',INDXC(1:NC+1) )


            ! ...... Depth, area and Volume arrays:

      call do_RealArray('DEPTH' ,DEPTH (1:ND ))
      call do_RealArray('VTOTZ' ,VTOTZ (1:NBI))
      call do_RealArray('ZBOTMI',ZBOTMI(1:NBI))
      call do_RealArray('LSHORE',LSHORE(1:NBI))
      call do_IntegerArray('NLVOPN',NLVOPN(1:NBI))
      call do_RealArray('VFROPN',VFROPN(1:NBI))
      call do_RealArray('AREA'  , AREA   (1:NLI))
      call do_RealArray('BOTTOM', BOTTOM (1:NLI))
      call do_RealArray('VFRAC' , VFRAC  (1:NLI))
 

   ModelState_Action = Action_Aborted

   end function


! ==============================================================

   subroutine do_IntegerDim(Name, Value)
   character*(*) Name
   integer Value

   character*32 Text
   if (Action_Aborted) return
   
   select case (CurrentAction)

      case (ModelState_Printout)
         write(FileUnit,*, IOSTAT=IO_Status) Name, Value

      case (ModelState_Save)
         write(FileUnit, IOSTAT=IO_Status) Name, Value

      case (ModelState_Restore)
         read(FileUnit, IOSTAT=IO_Status) Text
         if(IO_STatus==0) then
            if (Text == Name) then
               read(FileUnit, IOSTAT=IO_Status) Value
            else
               Action_Aborted = Name_Error      
            endif
         endif

   end select

   end subroutine

! ==============================================================

   subroutine do_IntegerArray(Name, Array)
   character*(*) Name
   integer Array(:)

   if (Action_Aborted) return
   
   select case (CurrentAction)

      case (ModelState_Printout)

      case (ModelState_Save)

      case (ModelState_Restore)
   end select

   end subroutine

! ==============================================================

   subroutine do_RealArray(Name, Array)
   character*(*) Name
   real*8 Array(:)

   if (Action_Aborted) return
   
   select case (CurrentAction)

      case (ModelState_Printout)

      case (ModelState_Save)

      case (ModelState_Restore)

   end select

   end subroutine


! ==============================================================

   subroutine do_RealArrayStructure(Name, IndexArray, Array)
   character*(*) Name
   Integer IndexArray(:) ! n+1 index values to array
   real*8 Array(:)         ! n segments of values
                         ! segment no. i from IndexArray(i)+1 to IndexArray(i+1) 

   integer MaxLength, NumberofSegments, i

   if (Action_Aborted) return

   select case (CurrentAction)

      case (ModelState_Printout)
         NumberofSegments = ubound(IndexArray,1)-1
         MaxLength=0
         do i = 1, NumberofSegments
            MaxLength = Max(MaxLength, IndexArray(i+1)- IndexArray(i))
         end do
         call PrintOut_RealArrayStructure(MaxLength, NumberofSegments)
            
      case (ModelState_Save)

      case (ModelState_Restore)

   end select

   return
   
   contains

   ! -------------- internal function for printout:
   
	   subroutine PrintOut_RealArrayStructure (MaxRows, Columns)
      integer MaxRows, Columns

      character*12 ValueString(MaxRows, Columns)
      integer i, k, kbase

            ! Fill text string matrix:

      Valuestring = " "
      do i = lbound(IndexArray,1), ubound(IndexArray,1)-1
         kbase = IndexArray(i)         
         do k=IndexArray(i)+1, IndexArray(i+1)
            write(ValueString(k-kbase,i),'(G12.4)', IOSTAT = IO_Status) Array(k)
            if (IO_Status.ne.0) then
               Action_Aborted = IO_ERROR
               return               
            endif
         end do
      end do

            ! Print text string matrix, with subsets beside each other in columns:

      write(FileUnit,*) Name,":"

      do i=1, MaxRows
         write(FileUnit,*, IOSTAT = IO_Status) (Valuestring(i, k), k=1, Columns)
         if (IO_Status.ne.0) then
            Action_Aborted = IO_ERROR
            return               
         endif
      end do       
	         
	   end subroutine
         
   end subroutine
   
end Module   