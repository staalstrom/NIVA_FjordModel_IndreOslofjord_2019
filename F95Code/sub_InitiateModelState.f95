Module sub_InitiateModelState

      use ModelParam_InitState
      use ModelParam_RunControl
      use ModelParam_Topography
      use ModelParam_Plankton
      use ModelParam_Decomposition
      use ModelParam_Boundary
      use ModelVar_HydroBioChem
      use ModelVar_Topography
      use ModelVar_RunControl

      use m6_FileNames
      use fx_Trhitr
      use fx_Eutrosub
      use fx_Mbalan
      use fx_SurfExch
      use fx_MusselInitialise
      use fx_Topogr, only: Topo
      use fx_Boundary

      
      implicit none

   contains

   
      subroutine InitiateModelState
      
      integer N_B, N_L


      call Hello("InitiateModelState")
      
      REINTG = .true.

      TZERO  = TSTART ! Initial condition for integrated time TINTEG
      TDERIV = TSTART ! Time for last call to CNCADJ,
                      ! set here to avoid updating of state variables
                      ! in initial calls to CNCADJ

      TINTEG = TSTART  ! Time as integrated value, used by CNCADJ to detect
                      ! when corrections to integrated terms are required


      MAXINT = 0.001  ! Initiate timestep control for each new run:
      MXTBIO = 0.001
                               ! PRVSTP = 0.001 (not used?)

                      ! ( Initiate time_step to small value)
      TCOMM = TSTART

      !  ----- Internal control values for FORTRAN subroutines

       ! Initiates some control values in TRHITR to ensure that
       ! successive runs without changes give identical results,
       ! see Module TRHITR

      CALL TRHITR_INIT


       ! Initiates file names for input data
       ! found in Module INIT.FOR:

      CALL Initate_FileNames


  ! ################################################################
      CALL DEBUGF(0)  ! Close any debug file which may be open

  !  ---------- INITIATE RANDOM GENERATOR
      CALL RANDOM_SEED ( NSEED )  ! In RANDFNCT.FOR

  !  ------------- INITIALIZE INPUT OF METEOROLOGICAL DATA ------------
      CALL METINI

  !  --------------------------- Topography ---------------------------

  ! ################################################################
  ! ---------- Skips calculation of topography if already done
  !            and specs. are unchanged:

      IF (.not.(STARTD .and. NLIMAX .eq. NLMPRV .and. DPFACT .eq. DPFOLD &
          .and. TOPO_NR.eq. TOPO_OLD )) THEN

  ! ---------- Establish topographic description -----------

          NLI_lim = MIN( dimMLI,NLIMAX )
!!          CALL TOPO( TOPOGRAPHY_FILE, &
!!                 TOPO_NR, dimMD, dimMBI, dimMLI, NLI_lim,  &
!!                 dimMBE, dimMLE, dimMC, dimMLC, &
!!                 DDPMIN, DPFACT, ND, DEPTH, &
!!                 NBI, VTOTZ, ZBOTMI, LSHORE, NLVOPN, VFROPN,  &
!!                 INDXI, AREA, BOTTOM, VFRAC, &
!!                 NBE, ZBOTME, INDXE, &
!!                 NC, WVDIR, ZSill, INDXC, BConn1,BConn2, WIDTH,  &
!!                 VBUFMX, VBUFTR, TCVBUF )

          CALL TOPO( TOPOGRAPHY_FILE, &
                 dimMD, dimMBI, dimMLI, NLI_lim,  &
                 dimMBE, dimMLE, dimMC, dimMLC, &
                 ND, DEPTH, &
                 NBI, VTOTZ, ZBOTMI, LSHORE, NLVOPN, VFROPN,  &
                 INDXI, AREA, BOTTOM, VFRAC, &
                 NBE, ZBOTME, INDXE, &
                 NC, WVDIR, ZSill, INDXC, BConn1,BConn2, WIDTH,  &
                 VBUFMX, VBUFTR, TCVBUF )


          NLI = INDXI(NBI+1)   ! Total number of layers
          NLE = INDXE(NBE+1)
          NLC = INDXC(NC+1 )
          NLMPRV = NLIMAX
          DPFOLD = DPFACT
          TOPO_OLD = TOPO_NR
      ENDIF


 ! ---- ZERO VALUE ARRAY, INITIAL VALUES FOR INTEGRATING TOTAL VOLUMES:
      CALL SETR ( 0.0, VTZERO, dimMBI )

 ! ----- Initialize layer volumes used in mass balance calculations
 !       during initialization:
      CALL VLCALC (VTZERO)

 ! ----- Initiate array for lagged total volume
 !       later used and updated in subroutine CNCADJ
      DVTOTL(1:NBI) = VTZERO(1:NBI)

 ! ----- Initiate buffer volumes and dynamic volumes:
      CALL SETR ( 0.0, VBUFZ,  2*dimMLC )
      CALL SETR ( 0.0, VDYNZ,  2*dimMBI )

 ! ----- Initial condition for integrated volume check variables:
      CALL SETR ( 0.0, VBFSZ,  dimMBIplus1 )
      CALL SETR ( 0.0, DVT2Z,  dimMBIplus1 )
      CALL SETR ( 0.0, VTDZ,   dimMBIplus1 )

!  ========== Conservation check on heat and mass substances =========
      CALL SETR (0.0, SALTB,  dimMBI)
      CALL SETR (0.0, HEATB,  dimMBI)
      CALL SETR (0.0, OXYGB,  dimMBI)
      CALL SETR (0.0, NITRB,  dimMBI)
      CALL SETR (0.0, PHOSB,  dimMBI)
      CALL SETR (0.0, SILIB,  dimMBI)
      CALL SETR (0.0, SPPB ,  dimMBI)

!  ====================== Water concentrations ========================

!  *************** First, zero all values

      CALL SETR (0.0, SAL,  dimMLI)
      CALL SETR (0.0, TEMP, dimMLI)
      CALL SETR (0.0, OXYG, dimMLI)
      CALL SETR (0.0, PO4,  dimMLI)
      CALL SETR (0.0, NO3,  dimMLI)
      CALL SETR (0.0, NH4,  dimMLI)
      CALL SETR (0.0, SiO2, dimMLI)
      CALL SETR (0.0, CFYT, dimMLI*dimMFYTG)
      CALL SETR (0.0, NFYT, dimMLI*dimMFYTG)
      CALL SETR (0.0, PFYT, dimMLI*dimMFYTG)
      CALL SETR (0.0, CHL,  dimMLI*dimMFYTG)
      CALL SETR (0.0, SFYT, dimMLI)

          ! Si only in first group of phytoplankton

  !  --------- Total concentrations of C,N and P in water
  !            are computed by CNCADJ

      CALL SETR (0.0, TOTC, dimMLI)
      CALL SETR (0.0, TOTN, dimMLI)
      CALL SETR (0.0, TOTP, dimMLI)


!  ============= Detritus concentrations =================

      CALL SETR ( 0.0, CDFLXI,  dimMLI )
      CALL SETR ( 0.0, NDFLXI,  dimMLI )
      CALL SETR ( 0.0, PDFLXI,  dimMLI )
      CALL SETR ( 0.0, SDFLXI,  dimMLI )

!  ============= Bottom sediment concentrations =================


        DO N_B = 1,NBI
          DO N_L = INDXI(N_B)+1, INDXI(N_B+1)
            CSED(N_L)=CSEDIN(N_B)
            NSED(N_L)=NSEDIN(N_B)
            PSED(N_L)=PSEDIN(N_B)
            SSED(N_L)=SSEDIN(N_B)
            RSED(N_L)=RSEDIN(N_B)
            ASED(N_L)=ASEDIN(N_B)
             ! set amount of retained P in sediment to maximum
            PADS(N_L)=PADMAX
             ! Initial state of sediment
            XSED(N_L)=XSEDZ(N_B)
            XBUR(N_B)=XSEDZ(N_B)
          End do
        End do

! ======================================================================
! Inorganic particles, added March 2001 BBJ
    CALL SETR ( 0.0, SPPFLXI, dimMLI )

    DO N_B = 1,NBI
      DO N_L = INDXI(N_B)+1, INDXI(N_B+1)
        SPPSED(N_L)=0.0
      end do
    end do
    call SETR(0.0, SPP, dimMLI)



 !************* Then initiate the part that is used with real data

   FYTGRP = MAX( 0, MIN (LFYT, dimMFYTG) )

   if (initc) then
      CALL ZHYDR ( EXTEST , &
          TSTART,  &
          NBI, INDXI, dimMLI, &
 ND, &
 DEPTH, &
 BOUND_INFL,  &
          PO4IN, NO3IN, NH4IN, SIO2IN, BIOOFF,  &
          FYTGRP, dimMFYTG, CFYTIN, NFYTIN, PFYTIN, SFYTIN,  &
          DOCIN, CZOOIN, BACTIN,  &
          NCBACT, PCBACT, NCZOO, PCZOO,  &
          BndFac_N, BndFac_P, BndFac_Ox, &
      SAL, TEMP, PO4, NO3, NH4, SiO2, OXYG,  &
      CFYT, NFYT, PFYT, CHL, SFYT,  &
      ODM, DOC, CZOO, BACT, CDET, NDET, PDET, SDET, RDET )

      CALL SIGMAT( SAL, TEMP, INDXI(NBI+1), DENSI)


   endif   



!  ====================== mussel population ============================

         CALL MSINIT (dimMBI, dimMUSLAGES, dimMUSLAYERS,  &
                      MUSLNR, MUSLMT, MUSLMA, MUSLGT )

 !  Initiate integrated values of import/export

        CALL SETR ( 0.0,   DNITRI,   dimMBI )
        CALL SETR ( 0.0,    NFIXI,   dimMBI )
        CALL SETR ( 0.0,   CLOADI,   dimMBI )
        CALL SETR ( 0.0, ODMLOADI,   dimMBI )
        CALL SETR ( 0.0,   NLOADI,   dimMBI )
        CALL SETR ( 0.0,   PLOADI,   dimMBI )
        CALL SETR ( 0.0,   CSEDXI,   dimMBI )
        CALL SETR ( 0.0,   NSEDXI,   dimMBI )
        CALL SETR ( 0.0,   PSEDXI,   dimMBI )

!  ====================== mass and heat conservation ===================
   !  xxxxMZ: 1. Size scale, accumulated throughout simulation,
   !          2. Base value of amount xxxxMI
        CALL MBINIT (dimMBI)  ! IN INIT.FOR

  !   calculate amounts in each basin, use as inital values
  !   for integration imports in time as means to mass balance control.
      CALL MBALAN( dimMBI, dimMLI, .false.,  &
             SAL,  TEMP,  OXYG,  NO3, NH4, PO4, SIO2,  &
             CFYT, NFYT, PFYT, SFYT,  &
             CDET, NDET, PDET, SDET,  &
             CSED, NSED, PSED, SSED, ASED, PADS,  &
             ODM, DOC , BACT, CZOO, CTMUSL, SPP, SPPSED, &
           SALTMI,HEATMI,OXYGMI,NITRMI,PHOSMI,SILIMI,SPPMI )

       ! requires VLAYER to be set


  !  ---- Initiate unity concentrations for continuity check
  !       on volume or for residence time of water within system
         CALL SETR (C1XTRN, C1EX,   dimMLE )
         CALL SETR (C1XTRN, C1,     dimMLI )


  !  --------- Initiate accumulated net buoyant influx effect
  !            and zero indexes (real*8) for well_mixed volumes,
  !            necessary for first call to TRANSP in section EULER.
        CALL SETR (0.0, BFXZ, dimMBI )
        CALL SETR (0.0, XMIX, dimMBI )
                !  zero values handled by TRANSP


   !  ------- initiate radiation and wind friction integration terms
       CALL SETR  (0.0, QOUTZ, dimMBI*5)
       CALL SETR  (0.0, QABSZ, 2   )
       CALL SETR  (0.0, FR3NTZ, dimMBI)
       CALL SETR  (0.0, UFR3I,  dimMBI)



   !----------------------------------------------------------------------
   ! ---- Initiate state variables for volume measures of basins:
  
      DVTOT = VTZERO
      VDYN  = VDYNZ
      VBUF  = VBUFZ

      	 ! values to compare with VBFSUM, DVTOT2 and VTDIFF

      VBFSI = VBFSZ
      DVT2I = DVT2Z
      VTDI  = VTDZ


			   ! integrated time (remnant from old ACSL code, not needed now?)
      SUMTIM =  0.0   

            ! Initiate accumulated wind mix energy/area 
            ! and Wind frictional velocity in 3rd exponent

      FR3INT = FR3NTZ

            !  ------- Buoyancy energy balance:

            ! initiate Accumulated net buoyant influx effect: (m2/s3*day)
      BFXINT = BFXZ

          !  -------------- INTEGRATE CLIMATE VARIABLES -------------

          ! RADIATION BUDGET TERMS:
	   QOUTI = QOUTZ
	   QABSI = QABSZ
	   QGI   = 0.0

          ! AIR TEMPERATURE AND WIND TERMS
      AIRTI = 0.0
      WN2I  = 0.0
      WE2I  = 0.0

      STARTD = .TRUE.

      end Subroutine



! =================================================================
! Moved from old file INIT.FOR

      SUBROUTINE MBINIT (MBI)
      INTEGER MBI

! Initiate mass balance values

      INTEGER I_B, I_2

!   Zero size scale values :

      DO I_B = 1,MBI

         SALTMZ(I_B) = 0.0
         HEATMZ(I_B) = 0.0
         OXYGMZ(I_B) = 0.0
         NITRMZ(I_B) = 0.0
         PHOSMZ(I_B) = 0.0
         SILIMZ(I_B) = 0.0
         SPPMZ (I_B) = 0.0

         DO I_2 = 1,2

             SALTMI(I_B,I_2) = 0.0
             HEATMI(I_B,I_2) = 0.0
             OXYGMI(I_B,I_2) = 0.0
             NITRMI(I_B,I_2) = 0.0
             PHOSMI(I_B,I_2) = 0.0
             SILIMI(I_B,I_2) = 0.0
             SPPMI (I_B,I_2) = 0.0

         ENDDO

      END DO

      END subroutine

end Module
