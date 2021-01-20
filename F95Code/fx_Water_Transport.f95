! ===================================================================
! NIVA Fjord model
! File: EUTRONEW.CSL
! Birger Bjerkeng, NIVA.
! Main ACSL code for model
! fetched from former ACSL code for main program of model
! ===================================================================
      Module fx_Water_Transports
	  
      use ModelDimensions
      use ModelParam_RunControl
      use ModelParam_Inputs
      use ModelParam_Physics
      use ModelVar_HydroBioChem
      use ModelVar_RunControl
      use ModelVar_Topography
      use fx_Transp_1
      
      implicit none
	  
	  
	  contains
C ------------------- Water transport calculation ----------------

      SUBROUTINE TRANSP

C  IMPORTS ACSL COMMON BLOCK, for use in call to WATER_TRANSPORTS:
!      INCLUDE 'EUTRODIM.inc' !B
!      INCLUDE 'EUTRO.INC'  !  ACSL Common block


      integer I,K

!!         ! Ensures that all ACSL derivatives are initiated,
!!         ! even those not used by the model:
!!      integer BASINS_INITIATED /0/,
!!     &        CONNECTIONS_INITIATED /0/
!!      SAVE BASINS_INITIATED, CONNECTIONS_INITIATED


      IF (.NOT. TROFF) THEN
C     ----- real calculation of transports:

C>>>>>> moved to Transp_1
!!C     ..... Rescale tidal mixing coefficient to apply to N2=1,
!!C           and multiply by energy factor:
!!         DO I = 1,NBI
!!C        ......... tidal mixing:
!!            MIXCONST(I) =  EMIXRL*MIXCF(I)*(N2SCAL**(MIXEXP/2.0))
!!         ENDDO
C<<<<<<<


!!         CALL WATER_TRANSPORTS ( HTROFF, (TRTEST .and. T.ge.TTRIG),
!!     &          ITRZ, DEPTH, dimMBI, NBI, VDYN, ZBOTMI,
!!     &          DTJETM, QWSURF, dimMS, BASINQ, DEPTHQ,
!!     &          QWATER, QTEMP, QSPP, dDens_dSPP, RNFNDX, QDIAM, NHOLES,
!!     &          XMIX, INDXI, NLI, AREA, SAL, TEMP, SPP, DENSI,
!!     &          VLAYER, VFRAC, NLVOPN, VFROPN, VLCORF,
!!     &          MIXFAC, MIXEXP, N2LIM,
!!     &          SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
!!     &          NBE, INDXE, NLE, DENSEX, ZSURFE, DZDTX,
!!     &          NC, INDXC, NLC, BConn1, BConn2,
!!     &          ZSILL, WIDTH, DPEFF, HTRMIX,
!!     &          dimMLC, VBUF, VBUFMX, VBUFTR, TCVBUF, MAXINT, T,
!!     &      UFLOW, VTOTDV, VDYNDV, VBUFDV, BSFLUX, MAXTTR, RCQNDX,
!!     &      Ambient_Volume_Flux, Neutral_Depth,
!!     &      BWFREQ, ZMID )

         CALL WATER_TRANSPORTS ( (TRTEST .and. T.ge.TTRIG),
     &          dimMBI, dimMS, QWATER, DENSEX, dimMLC, MAXINT, MAXTTR)



C        (maximum timestep = 10.0 days, reduced by transport calc.
C         and returned via MAXTTR)
C         now uses MAXINT instead.

C                       Mixing coefficient varies over time
C                       with input of energy through tides.


!!          BASINS_INITIATED      = max( BASINS_INITIATED, NBI)
!!          CONNECTIONS_INITIATED = max( CONNECTIONS_INITIATED, NC)

      ELSE

C    ------ No transports, set necessary values:
         DO I = 1, NBI
            DO K = INDXI(I)+1,INDXI(I+1)
               ZMID(K) = (DEPTH(K-INDXI(I))+DEPTH(K-INDXI(I)+1))
            ENDDO  ! must be set, is used by PHYT_ZOO
                   ! (normally done in TRANSP_V)
         ENDDO
         MAXTTR = MAXINT

      ENDIF

C     ---------------------------------------------------------------
C     ----- Unused output derivatives are zeroed (all if no transports)

!!      DO I = BASINS_INITIATED+1, dimMBI
!!         VTOTDV(I) = 0.0
!!         DO K=1,2
!!            VDYNDV(K,I) = 0.0
!!         ENDDO
!!         BFX (I) = 0.0
!!      ENDDO
!!      BASINS_INITIATED = max( BASINS_INITIATED, dimMBI)
!!
!!
!!      DO I = CONNECTIONS_INITIATED+1, dimMLC
!!         DO K=1,2
!!            VBUFDV(K,I) = 0.0
!!         ENDDO
!!      ENDDO
!!      CONNECTIONS_INITIATED = max( CONNECTIONS_INITIATED, dimMLC)


C>>>>>>> no longer necessary to set unused members
C        (was needed when ACSL was used): dimMS changed to NS

      DO I = 1, NS
         if (TROFF) RCQNDX(I) = RNFNDX(I,1)
              ! No vertical displacement of jets
         NQDIST(I) = 1 + RCQNDX(I) - MIN (RNFNDX(I,1),RCQNDX(I))
              ! Number of layers sharing land runoff, from RCQNDX and up:
              ! (see TRANSP_2.FOR)
      ENDDO


      TRCALC = .NOT. TROFF  ! Will affect TRNADJ and MTRANS
      END subroutine
	  
      end module
	  



