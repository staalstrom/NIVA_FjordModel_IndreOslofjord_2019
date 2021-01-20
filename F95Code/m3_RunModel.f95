Module m3_RunModel

$undefine Debug_TimeStepControl

   use sub_InitiateModelState
   use sub_Integration
   use sub_SnapShot 
   implicit none

   contains

!=======================================================   
   subroutine StartRun

$if defined Debug_TimeStepControl
   write(*,*)"Starts model run"
$endif

   CALL STOICH ! reset stoichiometric ratios

   T = TSTART
   call InitiateModelState
   call Hello ("After InitiateModelState")

   call Calc_Derivatives
   call Model_Snapshot


   end subroutine
   
   
!=======================================================   
   subroutine ContinueRun
   
   real*8 Model_First_Time, Model_Time_Step
   real*8 cpu_start, cpu_end

   real*8 T_BeginStep


$if defined Debug_TimeStepControl
   write(*,*)"Continue model run"
$endif

   CALL STOICH ! reset stoichiometric ratios

   TTERM = max(T,dble(TSTART + DTTERM))
   
   Model_First_Time = T
   call cpu_time (cpu_start)

   NDERIV = 0                        ! for counting derivative calculations

   do 

      if (T.ge.TTERM) exit 
      
            ! Set communication interval (time between snapshots)
            ! based on user specifications:
      
	   CINTV = MAX( dble(CIMIN), (T-TSTART)*CITFAC)
	   CINTV = MIN (CINTV, dble(CIMAX))
	
	   IF( CINTV.EQ.CIMAX .and. MOD(CIMAX,1.0).eq.0.0 .and. CIMAX.ge.1.0) then
	      CINTV = CIMAX*(1.0+AINT(T/CIMAX)-T/CIMAX) + MOD(CIBASE,1.0)
		      ! conditional expression shortens CINTV so that next snapshot
		      ! will be at a user-specifie time of day:  T+CINTV = n+p
		      ! with n = integral number (days), p = specified part of day:
		      !      p = mod(CIBASE,1.0) in interval [0;1)
      ENDIF

		CINTV = MAX ( 0.01d0, CINTV)
		CINTV = MIN (CINTV, TTERM-T )

				   ! Calculate next time for snapshot (communication):

      TCOMM = T + CINTV

               ! Integrate model up to next communication time:

      do while (T.lt.TCOMM)

         T_BeginStep = T

               ! set integration time step (previously in subroutine STPADJ):         
					! Adjust maximum time step not to go beyond next communication time
					! or final time, but always as a significant increment of current T
		         ! to avoid problems with eternal loop.
					! (SHOULD BE IMPROVED - DOES NOT ALWAYS WORK)

	      MAXINT = STPLIM ! user defined timestep limit

$if defined Debug_TimeStepControl
   write(*,*)"Calc_Derivatives at T=", T
$endif

         call Calc_Derivatives
         NDERIV = NDERIV + 1

	            ! Limit to reach next communication time or stopping time:
	            
	      MAXINT = MIN( MAXINT, TTERM-T, TCOMM-T, MXTBIO,TSTEP )
	
	                      ! Ensure significant time step: bruke nearest(x,s)?
	      MAXINT = MAX( T*1.e-6, MAXINT )

         CALL TIMPRT("After derivative calculation")

         T = MIN(TTERM, T + MAXINT)   ! integration will update status to this time (in CNCADJ)
         
         call Integrate_Step

         if (T_BeginStep == T) then
            write(*,*)'T=',T, 'MAXINT=',MAXINT,'too small TSTEP to change time'
!            Pause 'Ends run; press Enter to finish'
            return
         endif
      end do

      CALL TIMPRT("Before Model Snapshot")
      call Model_Snapshot

   end do

   CALL TIMPRT("At completion")

   call cpu_time (cpu_end)
   if (NDERIV.gt.0 ) then
      Model_Time_Step = (T - Model_First_Time)/NDERIV
      write(*,*) NDERIV,' integration steps from', Model_First_Time, ' to ', T
      write(*,*)' Model time span:', (T - Model_First_Time), ', time step:', Model_Time_Step, '(days)'
      write(*,*)' Cpu time used (seconds): ', cpu_end - cpu_start
      write(*,*)' Model days per Cpu second:' , (T - Model_First_Time)/(cpu_end - cpu_start)
   endif

   end subroutine

!=========================================================================   
   subroutine EndRun
$if defined Debug_TimeStepControl
   write(*,*)"Ends model run (cannot be continued)"
$endif


   CALL HELLO ('End of model run - first calls CNCADJ')
   CALL CNCADJ ( dimMBI, dimMLI, dimMUSLAGES, dimMUSLAYERS )


   end subroutine

end Module
