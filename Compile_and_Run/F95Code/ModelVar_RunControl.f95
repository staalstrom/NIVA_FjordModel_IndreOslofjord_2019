Module ModelVar_RunControl

      implicit none

!------------------- Model run control variables:

      logical STARTD 

      real*8 T      ! independent integration variable (time in days)

      real*8 SUMTIM ! integrated time (derivative =1)
                 ! probably not needed in new model version    

      real*8 YEARS ! T converted to Years

      real*8 TZERO  ! Initial condition for integrated time TINTEG

      real*8 TDERIV
         ! Time for last call to CNCADJ,
         ! set here to avoid updating of state variables
         ! in initial calls to CNCADJ

      real*8 TINTEG
         ! Time as integrated value, used by CNCADJ to detect
         ! when corrections to integrated terms are required

      real*8 TINTGZ  ! Should be hidden inside subroutines


      real*8 TSTEP
			! TSTEP specified in TRANSP to avoid numeric instabilities due
			! to transports in general, below limit MAXINT
			! which isset by other considerations
                        ! and then adjusted further in LGTRAD.

			
		real*8 MAXSTP  ! maximum step, set by Calc_Derivatives
      real*8 MAXTTR  ! maximum step limited by transports

      real*8 TCOMM

      real*8 TTERM
     
      real*8 MAXINT
        ! max. integrating step
     
      real*8 MXTBIO
	      ! Maximum timestep allowed by biological processes
	      ! set by PRPROD after MTRANS, when integration
	      ! step is already fixed.  Can only be used for
	      ! diagnostic step control. The problem could be
	      ! removed if the derivatives where no longer
	      ! initiated by MTRANS, but instead set to 0, and the effect
	      ! of transports included in the updating in TRANSP_U.FOR
	      ! In that case the timestep-dependent adjustment of
	      ! transports could be postponed to after all derivatives
	      ! where calculated, and done just before integrating over
	      ! the time-step, at the end of the derivative section

      LOGICAL BIOACT
    

      real*8 CINTV


      INTEGER*8 NDERIV
      real*8      TIMFAC  !  Calculate ratio between simulated and machine time
      real*8      DTMEAN  !  Mean length of integration step (simulated days)
                      !  Calculated for last communication interval:
      


end Module ModelVar_RunControl