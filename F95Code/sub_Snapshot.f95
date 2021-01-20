Module sub_Snapshot

! Revised code from Dynamic section of old ACSL model.
! --- Prepare and report snapshot of selected model state variables
!     at user-specified intervals

   use ModelParam_RunControl
   use ModelParam_Physics
   use ModelVar_RunControl
   use ModelVar_HydroBioChem
   use ModelVar_Topography
   use Binary_Write
   use fx_Eutrosub
   use sub_WrtRep
   
   implicit none
      
   integer :: FileUnit = 2000
   
contains
   
   subroutine Model_Snapshot

   
   ! Process results at communication intervals,
   ! recording intermediate state of simulation

   integer N_B, N_L
	real*8 AIRTD, WN2D, WE2D, UFR3D(dimMBI)


   CALL HELLO( 'ModelState_Snapshot'      )




      IF (MBRSET) CALL COMCLC ( 0, dimMBI, dimMLI )

     
   ! If .true. : Will print report (now mass budget)
   !             being reset to .false. after print

      IF (REPORT) then
         CALL WRTREP(FileUnit)
        REPORT = .FALSE.
      endif


     CALL HELLO ('final call to CNCADJ before snapshot')
     CALL CNCADJ(dimMBI,dimMLI, dimMUSLAGES, dimMUSLAYERS )

	      ! Before recording/output:
	      ! Make final calculation of wind mixing depth,
	      ! adjust integrated concentrations to changed volumes
	      ! This is also done by section EULER, before calculation of
	      ! transports, to ensure that correct values of water
	      ! concentrations are used for calculating derivatives
	      ! for integrating step.

             ! CALL VLCALC (DVTOT)
      CALL HELLO ('RETURNED FROM cncadj')


			! Calculate conservation check on heat and mass substances,
			! and handle very small VBUF values.

      CALL HELLO ('calls COMCLC')
      CALL COMCLC ( 0, dimMBI, dimMLI )
      CALL HELLO ('returns from COMCLC')


			! Calculate check on volume balance --------------

      CALL HELLO ('VOLCHK'      )
      CALL VOLCHK ( 2, NBI, NC, INDXC, BCONN1, BCONN2, NLC, VBUF,  &
        VDYN, DVTOT, (VPRT.AND.(T.GE.TTRIG)), vdindx,VBFSUM, DVTOT2, VTDIFF)

		    !  VBFSUM is the net sum of volume displacement according to VBUF
		    !  DVTOT2 is sum of dynamic volume and buffer displacement, and
		    !  should equal DVTOT integrated below,
		    !  VTDIFF is the difference between DVTOT2 and DVTOT.
		    !  Values 1.. NBI applies to each basin,
		    !  while last value NBI+1 is for all internal basins


 !  set surface levels for comparison with external level ZSURFE

!	    write(*,*) ' set surface levels in Model_Snapshot'
	    DO N_B = 1,NBI
!	    	  write(*,*) ' N_B=',N_B
!	    	  write(*,*) ' N_L=',N_L
	         N_L = INDXI(N_B)+1
!	    	  write(*,*) ' N_L=',N_L
	         ZSURFI (N_B) = VDYN(1,N_B)/area(N_L)
!	    	  write(*,*) ' ZSURFI (N_B)=',ZSURFI (N_B)
	    end do


			 !  Save mean values of integrated terms
			 !  and reset for integration of next period.
			 !  Not used in model, just included for monitoring purposes

!	    	  write(*,*) ' calling MeanVector'
	    CALL MeanVector ( dimMBI*5 , QOUTI, SUMTIM, QOUT, QOUTD )
	    CALL MeanVector ( 2  , QABSI, SUMTIM, QABS, QABSD )
	    CALL MeanValue ( QGI  , SUMTIM, SUNRAD+DIFRAD, QGD )
	
      	 !  AIR TEMPERATURE AND WIND TERMS


	    CALL MeanValue  ( AIRTI, SUMTIM, AIRTMP , AIRTD)
	    CALL MeanValue  ( WN2I , SUMTIM,WINDN**2, WN2D)
	    CALL MeanValue  ( WE2I , SUMTIM,WINDE**2, WE2D)
	    CALL MeanVector ( NBI , UFR3I, SUMTIM, UFRIC3 , UFR3D )

	    SUMTIM = 0.0  ! Reset integrated time
	
	    CALL VRATIO ( dimMLI, NLI, FYTGRP, CHL, CFYT, CHLCF )
	    CALL VRATIO ( dimMLI, NLI, FYTGRP, NFYT, CFYT, NCFYT )
	    CALL VRATIO ( dimMLI, NLI, FYTGRP, PFYT, CFYT, PCFYT )
	    CALL VRATIO ( dimMLI, NLI, FYTGRP, CFYTDV, CFYT, CFDVR )
	    CALL VRATIO ( dimMLI, NLI, FYTGRP, GTN, GTEMP, NUTLIM )

	
			    ! ########### Storage of simulation results to binary file
			    !             calling subroutines in module Bin_Res.for

	    	  write(*,*) ' Checking Binfile and BinReset:', BinFile, BinReset
	    if (BinFile.gt.0) then
	      if( BinReset) then
	    	  write(*,*) ' calling StartBinaryStorage'
	         call STARTBINARYSTORAGE( BinFile)
	         BinReset = .false.
	      endif
	    	  write(*,*) ' calling BINARYRECORDS'
	      call BINARYRECORDS
	    endif 


    end subroutine Model_Snapshot
    

    SUBROUTINE VRATIO ( dim1, N1, N2, DIVIDEND, DIVISOR, RATIO )
      
      INTEGER dim1, N1, N2
      REAL*8 DIVIDEND(dim1,*), DIVISOR(dim1,*), RATIO(dim1,*)

      REAL*8 BIG / 9.87654321E+20/
      INTEGER I1, I2
      DO I1 = 1, N1
        Do I2 = 1, N2 
          IF ( ABS( DIVIDEND(I1,I2)) .lt. abs(divisor(i1,I2))*BIG ) THEN
             RATIO(I1,I2) = DIVIDEND(I1,I2)/DIVISOR(I1,I2)
          ELSE  ! Avoids overflow or division by zero:
             RATIO(I1,I2) = SIGN( BIG, DIVIDEND(I1,I2))*SIGN(1.0D0,DIVISOR(I1,I2))
          ENDIF ! (implicitly assumes zero to be a small positive value)
        Enddo 
      ENDDO
      END  Subroutine  

end module sub_Snapshot
