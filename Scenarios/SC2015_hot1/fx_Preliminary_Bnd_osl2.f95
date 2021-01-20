C Eutrophication model  - File  BOUNDARY.FOR

! ***** Revision 19 June 2006 Birger Bjerkeng *****
!       Factors for adjusting external N, P and Oxygen

$define BOUNDARY_TEST 1
      ! if defined: writes debug info if EXTEST is true.
      ! >1: includes detailed info about variable values
      ! at start and end of subroutine Hydrex


      Module fx_Boundary
      use ModelDimensions, only: dimMFYTG, dimMBE
      use fx_RunControl
      use fx_OxygSat
      use fx_Sigmat
      use fx_BndTidal
      use fx_Randfnct
      use m6_FileNames

      implicit none


C ======================================================================
C  Initiation of tables in common-block HYDREX_TAB for subroutine HYDREX

C Declaration of common-block for subroutine HYDREX in BOUNDARY.FOR,
C with values entered from BOUNDARY_FILE (see INIT.FOR)

C Statistical description of hydrographic conditions
C on monthly basis:
C -------------------------------------------------------------------

      integer NBT ! Number of external basins in table, >= NBE.
      parameter (NBT = 1)

C ------------------ Depths of pyknocline:
C                    based on data for 1973, 1974, 1987, 1988,
C                    years where wind data has been analyzed)
C ..... Wind adjusted means, rounded values:
      real*8 DM(12,2,NBT)   ! DM(i,j,k) month i, layer j, basin k
C                 j=1: depth (m) of upper limit of pyknocline,
C                  =2: =log10(dd), dd= thickness (m) of pyknocline

C ..... Dependence on wind:
      real*8 DBETA(2)
C          regression coeff. for D(i,..) against northerly wind.
C ..... Standard deviation of random variation:
      real*8 DDEV(2)

C ------------------ Water quality above and below pyknocline,
C                    based on data 1973-1988.
C Temperature:
C   means:
      real*8 TEMPM(12,2,NBT)
C   standard deviations:
      real*8 TEMPD(12,2,NBT)

C Salinity:
C   means:
      real*8 SALM(12,2,NBT)
C   standard deviations:
      real*8 SALD(12,2,NBT)

C Oxygen (ml/l):   (NB! convert unit!)
C   means:
      real*8 OXYGM(12,2,NBT)
C   st. deviations:
      real*8 OXYGD(12,2,NBT)

C Total phosphorus (µgP/l):
C   means:
      real*8 TOTPM(12,2,NBT)
C   st. deviations:
      real*8 TOTPD(12,2,NBT)

C  ( Ortho-phosphate (µgP/l)        NOT USED):
C   means:
C      real*8 PO4NM(12,2,NBT)
C   st. deviations:
C      real*8 PO4ND(12,2,NBT)

C Total nitrogen (µg/l):
C   means:
      real*8 TOTNM(12,2,NBT)
C   st. deviations:
      real*8 TOTND(12,2,NBT)

C  ( Nitrate (µgN/l)    not used)
C   means:
C      real*8 NO3NM(12,2,NBT)
C   st. deviations:
C      real*8 NO3ND(12,2,NBT)

C  ( Ammonium (µg/l)     not used)
C   means:
C      real*8 NH4M(12,2,NBT)
C   st. deviations:
C      real*8 NO3ND(12,2,NBT)

C Total silicon (µg/l):
C   means:
      real*8 TOTSM(12,2,NBT)

C   st. deviations:
      real*8 TOTSD(12,2,NBT)


C -------- All arrays collected into general array:
      INTEGER NPAR
      PARAMETER (NPAR=6)
C     --------- mean values:
      real*8 MV  (12,2,NBT,NPAR)
      EQUIVALENCE
     &    ( MV(1,1,1,1), TEMPM )
     &   ,( MV(1,1,1,2), SALM  )
     &   ,( MV(1,1,1,3), OXYGM )
     &   ,( MV(1,1,1,4), TOTPM )
     &   ,( MV(1,1,1,5), TOTNM )
     &   ,( MV(1,1,1,6), TOTSM )

C     &   ,( MV(1,1,1,7), PO4PM )
C     &   ,( MV(1,1,1,8), NO3NM )
C     &   ,( MV(1,1,1,9), NH4NM )

C     ---------- standard deviations:
      real*8 STDEV (12,2,NBT,NPAR)
      EQUIVALENCE
     &    ( STDEV(1,1,1,1), TEMPD )
     &   ,( STDEV(1,1,1,2),  SALD )
     &   ,( STDEV(1,1,1,3), OXYGD )
     &   ,( STDEV(1,1,1,4), TOTPD )
     &   ,( STDEV(1,1,1,5), TOTND )
     &   ,( STDEV(1,1,1,6), TOTSD )

C     &   ,( STDEV(1,1,1,7), PO4PD )


C  ------ Random function for each parameter:
      real*8 BNDR_T_LAST (2,NBT,0:NPAR) 
      real*8 RESP_FREQ   (2,NBT,0:NPAR)
      real*8 BNDR_RANDOM (2,NBT,0:NPAR)
      
      integer dim_Signal      
      parameter ( dim_Signal= 2*NBT*(1+NPAR) )
      Integer BNDR_IJ (2,  dim_Signal )      
      real*8    BNDR_UC (98, dim_Signal )      

C ------- Two-layer description for current time, found from tables,
C         interpolated in time, with wind dep. and random part added.
      real*8 CD(2,NBT)       ! Depth values
      real*8 CV(2,NBT, NPAR) ! Hydrographic values

      real*8 CSAL(2,NBT), CTEMP(2,NBT), CSILICON (2,NBT)
      EQUIVALENCE (CV(1,1,1), CTEMP(1,1)), (CV(1,1,2), CSAL(1,1))
      EQUIVALENCE (CV(1,1,6), CSILICON(1,1) )

C ------- Final current value for each parameter:
      real*8 TEMP_L
      real*8 SAL_L
      real*8 OXYG_L
      real*8 TOTP_L
      real*8 TOTN_L
      real*8 TOTS_L

      real*8 V_L(NPAR)
      EQUIVALENCE
     &    ( V_L(1),TEMP_L )
     &   ,( V_L(2),SAL_L  )
     &   ,( V_L(3),OXYG_L )
     &   ,( V_L(4),TOTP_L )
     &   ,( V_L(5),TOTN_L )
     &   ,( V_L(6),TOTS_L )


C  Statistical description of hydrographic conditions
C  on monthly basis:
C -------------------------------------------------------------------

C ..... Depths of pyknocline:
C          (based on data for 1973, 1974, 1987, 1988,
C           years where wind data has been analyzed)
C  Wind adjusted means, rounded values:
      DATA DM /
     &  9.0,  9.0, 18.0,  4.5,  2.5,  2.5,   !D (m), upper limit
     & 11.0,  6.0,  7.0, 11.0,  9.0,  8.0,
     & 3*1.3, 0.9, 6*0.7, 1.2, 1.3 /         !log10(thickness(m))


C  Dependence on wind:
      DATA DBETA /-1.4,0.0 /
C          regression coeff. for D(i,..) against northerly wind.

C  Response time of random variation:

      DATA RESP_FREQ/ dim_Signal*0.13 / ! 1/day

C  Standard deviation of random variation:
      DATA DDEV /6.0, 0.4/

C ............ water quality, based on data 1973-1988.
C Temperature:
C   means:
      DATA TEMPM /
     &  2.7,  2.3,  1.8,  5.8, 11.6, 15.7,  !upper
     & 19.2, 17.4, 15.4, 10.0,  5.9,  6.3,
     &  7.5,  6.0,  5.5,  5.0,  7.0,  9.0,  !lower
     & 12.0, 14.0, 13.0, 12.0, 10.5,  9.0 /
C   standard deviations:
      DATA TEMPD /  24*1.5 /

C Salinity:
C   means: endret - bedre tilpasningen med Breiangen
      DATA SALM /
     & 28.6, 27.1, 24.6, 23.8, 21.4, 22.3,
     & 20.0, 23.4, 24.2, 24.5, 27.5, 31.4,
     & 33.0, 33.5, 33.5, 33.5, 32.0, 31.0,
     & 31.0, 31.5, 31.5, 31.5, 32.5, 33.0 /

C   standard deviations:
      DATA SALD /
     &  3.0,  3.9,  2.0,  4.5,  3.8,  3.3,
     &  3.0,  2.8,  2.7,  2.7,  2.0,  2.0,
     &  1.0,  1.5,  1.5,  2.0,  4.5,  3.0,
     &  1.0,  1.5,  1.5,  2.0,  1.0,  1.0 /

C Oxygen (ml/l):   (NB! convert unit!)
C   means:
      DATA OXYGM /
     &  6.8,  7.4,  8.0,  8.0,  7.1,  7.0,
     &  6.0,  5.8,  6.2,  6.2,  6.2,  5.9,
     &  5.4,  5.6,  5.3,  5.3,  5.2,  5.0,
     &  4.9,  4.7,  4.0,  4.6,  4.7,  5.2 /
       !ENDRET LITT I FORHOLD TIL OPPRINNELIG SATT VERDI.

C   st. deviations:
      DATA OXYGD /
     &  0.5,  1.0,  0.6,  0.7,  0.8,  1.2,
     &  0.4,  0.5,  0.4,  0.5,  0.2,  0.2,
     &  0.4,  0.4,  0.3,  0.6,  1.4,  1.0,
     &  0.5,  0.5,  0.5,  0.5,  0.6,  0.5 /

C Total phosphorus (µgP/l):
C   means:
      DATA TOTPM /
     & 29.0, 30.0, 27.0, 25.0, 14.0, 11.0,
     & 10.0,  9.0, 13.0, 16.0, 25.0, 26.0,
     & 29.0, 29.0, 29.0, 29.0, 25.0, 17.0,
     & 24.0, 12.0, 13.0, 18.0, 20.0, 28.0 /
C   st. deviations:
      DATA TOTPD /
     & 12.0, 10.0, 10.0, 12.0,  8.0,  3.0,
     &  2.0,  4.0,  7.0, 15.0,  3.0,  8.0,
     & 4*4.0,                  10.0,  3.0,
     & 5.0, 8.0, 4*6.0  /

C  ( Ortho-phosphate (µgP/l)        NOT USED):
C   means:
C      DATA PO4NM /
C     & 20.0, 21.0, 16.0, 12.0,  5.0,  5.0,
C     &  5.0,  5.0,  5.0,  9.0, 20.0, 24.0,
C     & 29.0, 29.0, 29.0, 29.0, 25.0, 17.0,
C     & 24.0, 12.0, 13.0, 18.0, 20.0, 28.0 /
C   st. deviations:
C      DATA PO4ND /
C     & 12.0, 10.0, 10.0, 12.0,  8.0,  3.0,
C     &  2.0,  4.0,  7.0, 15.0,  3.0,  8.0,
C     & 4*4.0,                  10.0,  3.0,
C     & 5.0, 8.0, 4*6.0  /

C Total nitrogen (µg/l):
C   means:
      DATA TOTNM /
     & 320.0, 330.0, 400.0, 440.0, 280.0, 270.0,
     & 190.0, 210.0, 300.0, 300.0, 300.0, 250.0,
     & 3*240.0, 3*260.0, 6*210.0 /
C   st. deviations:
      DATA TOTND /
     & 6*80.0, 2*50.0, 4*100.0,
     & 3*30.0, 9*70.0      /

C  ( Nitrate (µgN/l)    not used)
C   means:
C      DATA NO3NM /
C     & 130.0, 150.0, 190.0, 140.0,  40.0,  20.0,
C     & 3*10.0,               60.0,  80.0, 110.0,
C     & 6*130.0, 6*80.0    /
C   st. deviations:
C      DATA NH3ND /
C     & 4*80.0, 5*20.0,       30.0,  40.0,  60.0,
C     & 6*60.0, 6*40.0      /
C                ----- Set to max(20., 0.5*mean)

C  ( Ammonium (µg/l)     not used)
C   means:
C      DATA NH4M /
C     & 2*15.0, 70.0, 2*40.0, 90.0, 4*30.0, 2*60.0
C     & 12*30.0            /
C   st. deviations:
C      DATA NO3ND /
C     & 2*15.0, 35.0, 2*20.0, 45.0, 4*15.0, 2*30.0
C     & 12*15.0            /
C                ----- Set to max(15., 0.5*mean)

C Silicate (µg Si/l):
      DATA TOTSM /
     & 4*1400., 1000., 700., 300., 100., 100., 500., 1000., 1200.,
     & 12*300.0 /
C          First row contains assumed values for fresh water,
C          will be rescaled according to salinity):
C              Based on following summarized monthly values for rivers
C              Glomma and Dramselva in Norway:
C                   4*3.0,2.5,4*1.5,3*3.0 mg/l SiO2/l
C          + assumptions about further reduction in summer.

C          Could use following relation for surface silicate
C          in the Oslo_fjord:
C              SiO2 = value for fresh water: approx. 1000µg/l
C              transform into function of salinity:
C              SiO2 = MAX(100.0, 1000.*(1-0.04*SAL_L),17*(SAL_L-25.0))
CV .... as for fresh water

C   standard deviations:
      DATA TOTSD / 12*60.0, 12*30.0 /
C Silicate in mgSiO2/l, same in all sources,

      
      contains

C ==================================================================
C Contains subroutines for:
C   1. Subroutines for initial and external
C      hydrophysical and biochemical variables
C   2. Basic hydrographic density relations
C =================================================================


C ===================================================================
C Initiate hydrographic profile:
      SUBROUTINE ZHYDR ( EXTEST, INIT_TIME,
     &                   NBI, INDXI, MLI, ND, DEPTH, BOUND_INFL,
     &                   PO4IN, NO3IN, NH4IN, SIO2IN, BIOOFF,
     &                   FYTGRP, MFYTG, CFYTIN, NFYTIN, PFYTIN, SFYTIN,
     &                   DOCIN, CZOOIN, BACTIN,
     &                   NCBACT, PCBACT, NCZOO, PCZOO,
     &                   BndFac_N, BndFac_P, BndFac_Ox,
     &             SALZ, TEMPZ, PO4Z, NO3Z, NH4Z, SiO2Z, O2Z,
     &             CFYTZ, NFYTZ, PFYTZ, CHLZ, SFYTZ,
     &             ODMZ, DOCZ, CZOOZ, BACTZ,
     &             CDETZ, NDETZ, PDETZ, SDETZ, RDETZ )
      
           
           ! Only first component of sediment conc. are
           ! initiated here.
C ------------------- subroutine arguments --------------------
C In:
      LOGICAL EXTEST
      real*4 INIT_TIME 
      INTEGER NBI, INDXI (NBI+1), MLI, ND
      real*8 DEPTH (ND)
      real*4 BOUND_INFL (NBI)

   ! -------- RELATIVE START CONCENTRATIONS OF NUTRIENTS
   !          AND BIOLOGICAL COMPONENTS (WILL BE SCALED ACCORDING
   !          TO TOTAL CONCENTRATIONS IN EXTERNAL WATERS)
      real*4 PO4IN, NO3IN, NH4IN, SIO2IN
      logical BIOOFF ! Controls whether biology is turned off
      INTEGER FYTGRP ! = active number of phytoplankton groups
      INTEGER MFYTG
      real*4 CFYTIN(MFYTG), NFYTIN (MFYTG), PFYTIN (MFYTG), SFYTIN
      real*4 DOCIN,  CZOOIN, BACTIN
      real*4 NCBACT, PCBACT, NCZOO, PCZOO
      
            ! factors for adjusting external N, P and Oxygen
      real*4 BndFac_N    ! Total nitrogen, factor on excess of 100 µg/l
      real*4 BndFac_P    ! Phosphorus, factor on value
      real*4 BndFac_Ox   ! Oxygen, factor on deviation from saturation  

C Out: initial values of hydrographic variables:
      real*8 SALZ (MLI), TEMPZ (MLI)
      real*8 PO4Z(MLI), NO3Z(MLI), NH4Z(MLI), SIO2Z(MLI)
      real*8 O2Z (MLI)
      real*8 CFYTZ (MLI,MFYTG), NFYTZ (MLI,MFYTG)
      real*8 PFYTZ (MLI,MFYTG), CHLZ  (MLI,MFYTG)
      real*8 SFYTZ (MLI) ! Silicate, only first phytoplankton group
      real*8 ODMZ(MLI), DOCZ(MLI), CZOOZ(MLI), BACTZ(MLI)
      real*8 CDETZ(MLI), NDETZ(MLI), PDETZ(MLI), SDETZ(MLI), RDETZ(MLI)

C --------------------- local variables ----------------------
      integer mfs
      parameter (mfs=2) ! Should be = number of phytoplankton groups

C Seed values to initiate distribution of nutrients,
C only values for one layer, used for all layers in model.

C For nutrients & C, only ratios between forms are important,
C not physical units.

      INTEGER INDXZ(2) / 0,1 /  ! Layer index limits
      real*8 TEMPS(1)    /0.0  /  ! set directly, not specified as model parameter
      real*8 PO4S(1)
      real*8 PFYTS(1,dimMFYTG)
      
C             approx. 5% of P in each PFYT group
      real*8 NO3S(1), NH4S(1)
      real*8 NFYTS(1,dimMFYTG)
C             approx 5% of N in each NFYT group
      real*8 SiO2S(1), SFYTS(1) ! 20% of Si in plankton

C Phytoplankton at approx. Redfield ratio for C:N:P
      real*8 CFYTS(1,dimMFYTG)

C Chlorophyll level corresponding to critical radiance 5W/m2
C at growth rate 0.25/day
      real*8 CHLS(1,dimMFYTG)

      real*8 DOCS(1), CZOOS(1), BACTS(1)
      
C Initiate without any sedimenting matter:
      real*8 :: CDETS(1)=(/0.0/), NDETS(1)=(/0.0/)

      real*8 :: PDETS(1)=(/0.0/), SDETS(1)=(/0.0/), RDETS(1)=(/0.0/)

      real*8 ZSURFE(dimMBE), DZDTX(dimMBE), EMIX_RELATIVE ! Dummy variables
      INTEGER NBSURF
      PARAMETER (NBSURF=1)


      real*8 WIND_N, WIND_E
      INTEGER TEMPLATE_BASINS, TEMPLATE_LAYERS
      real*8 TEMPLATE_AREA(2) /2*1.0/
      real*4 TIDAL_FACTOR
      LOGICAL FIXED_EXT_TEMP

      PARAMETER (
     &  WIND_N = 0.0, WIND_E = 0.0,
     &  TEMPLATE_BASINS = 1, TEMPLATE_LAYERS = 1,
     &  TIDAL_FACTOR = 0.0,
     &  FIXED_EXT_TEMP  = .TRUE. )

C ------------------------------------------------------------
      INTEGER N_F
      real*4 EXTBIO_INIT(2)

C Initiate values for internal basins, using the same subroutine
C as for external basins during simulation, to get approximate
C realistic conditions. The correspondence between external and 
C internal basins are not set up according to connections, but
C purely on basis of index.
C If number of external basins is less than number of internal basins,
C the last external basin is used for the final internal basins

C T<0 initiates random variables suppresses setting of surface levels:


      call Hello("ZHYDR")
      
$if defined  BOUNDARY_TEST
      if (EXTEST) write(debug_unit,*)'ZHYDR, INDXI',INDXI,' MLI',MLI
$endif


         ! Distribute scalar initial value parameters across arrays:
      PO4S = PO4IN
      NO3S = NO3IN
      SIO2S = SIO2IN
      NH4S = NH4IN
      do N_F = 1, MFYTG
         CFYTS(:,N_F) = CFYTIN(N_F)
         PFYTS(:,N_F) = PFYTIN(N_F)
         NFYTS(:,N_F) = NFYTIN(N_F)
         CHLS(:,N_F) = 0.25/5.0*CFYTIN(N_F)
      end do
      SFYTS = SFYTIN

      DOCS  = DOCIN
      CZOOS = CZOOIN
      BACTS = BACTIN

      IF (BIOOFF) THEN 
         EXTBIO_INIT(1) = 0.0
      ELSE   
         EXTBIO_INIT(1) = 1.0
      ENDIF
      EXTBIO_INIT(2) = EXTBIO_INIT(1)
      
      CALL HYDREX( EXTEST, dble(INIT_TIME), 0.0D0, WIND_N, WIND_E,
     &             NBSURF, NBI, INDXI, MLI, TEMPLATE_BASINS, INDXZ,
     &             TEMPLATE_LAYERS, ND, DEPTH, TEMPLATE_AREA,
     &             TIDAL_FACTOR, FIXED_EXT_TEMP,
     &             TEMPS, PO4s, NO3S, NH4S, SiO2S,
     &             BOUND_INFL, EXTBIO_INIT, MFYTG, FYTGRP,
     &             CFYTS, NFYTS, PFYTS, CHLS, SFYTS,
     &             DOCS, CZOOS, BACTS, NCBACT, PCBACT, NCZOO, PCZOO,
     &             CDETS, NDETS, PDETS, SDETS, RDETS,
     &             BndFac_N, BndFac_P, BndFac_Ox,
     &     ZSURFE, DZDTX, EMIX_RELATIVE,
     &     SALZ, TEMPZ, PO4Z, NO3Z, NH4Z, SiO2Z, O2Z,
     &     CFYTZ, NFYTZ, PFYTZ, CHLZ, SFYTZ,
     &     ODMZ, DOCZ, CZOOZ, BACTZ, 
     &     CDETZ, NDETZ, PDETZ, SDETZ, RDETZ )
         ! write(*,*)' tilbake til ZHYDR fra hydrex: CZOOZ=', CZOOZ
      END Subroutine



C ===================================================================
C Set hydrographic profiles in external basins:
      SUBROUTINE HYDREX( EXTEST, T,  AIR_PRESSURE, WINDN, WINDE,
     &             NBSURF, NBE, INDXE, MLE, NBI, INDXI, MLI,
     &             ND, DEPTH, AREA,
     &             TIDFAC, FIXTMP, TEMP, PO4, NO3, NH4, SiO2,
     &             BOUND_INFL, EXTBIO, MFYTG, FYTGRP,
     &             CFYT, NFYT, PFYT, CHL, SFYT,
     &             DOC, CZOO, BACT, NCBACT, PCBACT, NCZOO, PCZOO,
     &             CDET, NDET, PDET, SDET, RDET,
     &             BndFac_N, BndFac_P, BndFac_Ox,
     &      ZSURFE, DZDTX, EMIX_RELATIVE,
     &      SALEX, TEMPEX, PO4EX, NO3EX, NH4EX, SiO2EX, O2EX,
     &      CFYTEX, NFYTEX, PFYTEX, CHLEX, SFYTEX,
     &      ODMEX, DOCEX, CZOOEX, BACTEX,
     &      CDETEX, NDETEX, PDETEX, SDETEX, RDETEX )
      

C ---------------- subroutine arguments -------------------
C In: .... Current integration time (days)
      logical extest
      real*8 T
      real*8 Air_Pressure

C     .... wind speed, components, affects stratification:
      real*8 WINDN, WINDE  ! (m/s)
C     .... Basin/layer description:
      INTEGER NBE, INDXE(NBE+1), MLE   ! external basins
      INTEGER NBI, INDXI(NBI+1), MLI   ! internal basins
      integer nbsurf  ! number of
      INTEGER ND
      real*8 DEPTH(ND)
      real*8 AREA (MLI) ! Horizontal area of internal basins,
C                       used as weights when calculating
C                       nutrient distribution ratios
      real*4 TIDFAC     ! Degree of full tidal variation in ZSURFE
      LOGICAL FIXTMP  ! True if external temperature should be
C                       found in table, otherwise use:
      real*8 TEMP (MLI) ! temperature of external basins.

C     .... Biochemical depth profiles of internal basins,
C          used to distribute total external nutrients on forms:
      real*8 PO4(MLI), NO3(MLI), NH4(MLI), SiO2(MLI)
      real*4 BOUND_INFL(NBI) ! Weighs the influence of internal
                           ! basins on distr. inorganic/biological
                           ! in external (boundary) areas.

      real*4 EXTBIO(2)  ! Turns on/off biological components externally
                      ! (1): in general (2): additional factor for DOC
                      ! effective values kept within range 0-1.
      INTEGER MFYTG   ! dimensioning number of phytoplankton groups
      INTEGER FYTGRP  ! number of active phytoplankton groups
      real*8 CFYT(MLI,MFYTG), NFYT(MLI,MFYTG)
      real*8 PFYT(MLI,MFYTG), CHL(MLI,MFYTG)
      real*8 SFYT(MLI)   ! Silicate, only in first phytoplankton group
      real*8 DOC(MLI), CZOO(MLI), BACT(MLI)
      real*4 NCZOO, PCZOO, NCBACT, PCBACT
      real*8 CDET(MLI), NDET(MLI), PDET(MLI), SDET(MLI), RDET(MLI)

      
            ! factors for adjusting external N, P and Oxygen
      real*4 BndFac_N    ! Total nitrogen, factor on excess of 100 µg/l
      real*4 BndFac_P    ! Phosphorus, factor on value
      real*4 BndFac_Ox   ! Oxygen, factor on deviation from saturation  
      


C Out:.... Depth profiles of external basins:
C     .. Nutrients:
      real*8 SALEX (MLE), TEMPEX (MLE)
      real*8 PO4EX(MLE), NO3EX(MLE), NH4EX(MLE), SiO2EX(MLE)
C     .. Oksygen:
      real*8 O2EX (MLE)
C     .. Phytoplankton:
      real*8 CFYTEX(MLE,MFYTG), NFYTEX(MLE,MFYTG)
      real*8 PFYTEX(MLE,MFYTG), CHLEX(MLE,MFYTG)
      real*8 SFYTEX(MLE) ! Silicate only in first phytoplankton group
C  --------- other biological components:
      real*8 ODMEX(MLE), DOCEX(MLE), CZOOEX(MLE), BACTEX(MLE)
      real*8 CDETEX(MLE), NDETEX(MLE), PDETEX(MLE), SDETEX(MLE)
      real*8 RDETEX(MLE)

C     .. Change of surface level with time derivative:
      real*8 ZSURFE (NBSURF), DZDTX(NBSURF)
C        Units : SALEX = o/oo, TEMPEX = temp, DENSEX = sigma-t'
C     .. Relative input of mixing energy:
      real*8 EMIX_RELATIVE

C ==================================================================
C -------------------------- variables, tables ---------------------
!C Previously in common-block for subroutine HYDREX in BOUNDARY.FOR,
!C with values entered from BOUNDARY_FILE (see INIT.FOR)
!
!C Statistical description of hydrographic conditions
!C on monthly basis:
!C -------------------------------------------------------------------
!
!      integer NBT ! Number of external basins in table, >= NBE.
!      parameter (NBT = 1)
!
!C ------------------ Depths of pyknocline:
!C                    based on data for 1973, 1974, 1987, 1988,
!C                    years where wind data has been analyzed)
!C ..... Wind adjusted means, rounded values:
!      real*8 DM(12,2,NBT)   ! DM(i,j,k) month i, layer j, basin k
!C                 j=1: depth (m) of upper limit of pyknocline,
!C                  =2: =log10(dd), dd= thickness (m) of pyknocline
!
!C ..... Dependence on wind:
!      real*8 DBETA(2)
!C          regression coeff. for D(i,..) against northerly wind.
!C ..... Standard deviation of random variation:
!      real*8 DDEV(2)
!
!C ------------------ Water quality above and below pyknocline,
!C                    based on data 1973-1988.
!C Temperature:
!C   means:
!      real*8 TEMPM(12,2,NBT)
!C   standard deviations:
!      real*8 TEMPD(12,2,NBT)
!
!C Salinity:
!C   means:
!      real*8 SALM(12,2,NBT)
!C   standard deviations:
!      real*8 SALD(12,2,NBT)
!
!C Oxygen (ml/l):   (NB! convert unit!)
!C   means:
!      real*8 OXYGM(12,2,NBT)
!C   st. deviations:
!      real*8 OXYGD(12,2,NBT)
!
!C Total phosphorus (µgP/l):
!C   means:
!      real*8 TOTPM(12,2,NBT)
!C   st. deviations:
!      real*8 TOTPD(12,2,NBT)
!
!C  ( Ortho-phosphate (µgP/l)        NOT USED):
!C   means:
!C      real*8 PO4NM(12,2,NBT)
!C   st. deviations:
!C      real*8 PO4ND(12,2,NBT)
!
!C Total nitrogen (µg/l):
!C   means:
!      real*8 TOTNM(12,2,NBT)
!C   st. deviations:
!      real*8 TOTND(12,2,NBT)
!
!C  ( Nitrate (µgN/l)    not used)
!C   means:
!C      real*8 NO3NM(12,2,NBT)
!C   st. deviations:
!C      real*8 NO3ND(12,2,NBT)
!
!C  ( Ammonium (µg/l)     not used)
!C   means:
!C      real*8 NH4M(12,2,NBT)
!C   st. deviations:
!C      real*8 NO3ND(12,2,NBT)
!
!C Total silicon (µg/l):
!C   means:
!      real*8 TOTSM(12,2,NBT)
!
!C   st. deviations:
!      real*8 TOTSD(12,2,NBT)
!
!
!C -------- All arrays collected into general array:
!      INTEGER NPAR
!      PARAMETER (NPAR=6)
!C     --------- mean values:
!      real*8 MV  (12,2,NBT,NPAR)
!      EQUIVALENCE
!     &    ( MV(1,1,1,1), TEMPM )
!     &   ,( MV(1,1,1,2), SALM  )
!     &   ,( MV(1,1,1,3), OXYGM )
!     &   ,( MV(1,1,1,4), TOTPM )
!     &   ,( MV(1,1,1,5), TOTNM )
!     &   ,( MV(1,1,1,6), TOTSM )
!
!C     &   ,( MV(1,1,1,7), PO4PM )
!C     &   ,( MV(1,1,1,8), NO3NM )
!C     &   ,( MV(1,1,1,9), NH4NM )
!
!C     ---------- standard deviations:
!      real*8 STDEV (12,2,NBT,NPAR)
!      EQUIVALENCE
!     &    ( STDEV(1,1,1,1), TEMPD )
!     &   ,( STDEV(1,1,1,2),  SALD )
!     &   ,( STDEV(1,1,1,3), OXYGD )
!     &   ,( STDEV(1,1,1,4), TOTPD )
!     &   ,( STDEV(1,1,1,5), TOTND )
!     &   ,( STDEV(1,1,1,6), TOTSD )
!
!C     &   ,( STDEV(1,1,1,7), PO4PD )
!C     &   ,( STDEV(1,1,1,8), NO3ND )
!C     &   ,( STDEV(1,1,1,9), NH4ND )
!
!C  ------ Random function for each parameter:
!      real*8 BNDR_T_LAST (2,NBT,0:NPAR) 
!      real*8 RESP_FREQ   (2,NBT,0:NPAR)
!      real*8 BNDR_RANDOM (2,NBT,0:NPAR)
!      
!      integer dim_Signal      
!      parameter ( dim_Signal= 2*NBT*(1+NPAR) )
!      real*8 BNDR_IJ (2,  dim_Signal )      
!      real*8 BNDR_UC (98, dim_Signal )      
!
!C ------- Two-layer description for current time, found from tables,
!C         interpolated in time, with wind dep. and random part added.
!      real*8 CD(2,NBT)       ! Depth values
!      real*8 CV(2,NBT, NPAR) ! Hydrographic values
!
!      real*8 CSAL(2,NBT), CTEMP(2,NBT), CSILICON (2,NBT)
!      EQUIVALENCE (CV(1,1,1), CTEMP(1,1)), (CV(1,1,2), CSAL(1,1))
!      EQUIVALENCE (CV(1,1,6), CSILICON(1,1) )
!
!C ------- Final current value for each parameter:
!      real*8 TEMP_L
!      real*8 SAL_L
!      real*8 OXYG_L
!      real*8 TOTP_L
!      real*8 TOTN_L
!      real*8 TOTS_L
!
!      real*8 V_L(NPAR)
!      EQUIVALENCE
!     &    ( V_L(1),TEMP_L )
!     &   ,( V_L(2),SAL_L  )
!     &   ,( V_L(3),OXYG_L )
!     &   ,( V_L(4),TOTP_L )
!     &   ,( V_L(5),TOTN_L )
!     &   ,( V_L(6),TOTS_L )

      real*8 CDENS(2,NBT)  ! NBT defined in included file

!      EXTERNAL ONE_TO_ONE

C --------------- local variables ---------------------
      INTEGER IB, IBI, IBT, IG, I, L, LI, ID, IP, LD, LSURF, IW
      real*8 A1, A2, A, WPSUM, WNSUM, WSSUM, W_C_SUMFYT, C_SUMFYT_EXT
      real*8 W

      real*8 R
      INTEGER LOOP
      LOGICAL STABLE
      real*8 Z_EXT, DZ_DT_EXT, Z
      real*8 TEMP_SUM, area_sum

      real*8 OX_SAT_EX
       

C ...... Control random fluctuations (time unit 1 day):
      real*8 T_PREV /-1.E30/
      SAVE T_PREV

      real*8 T_M, T_A, T_Y
      INTEGER IT1, IT2

C ........... Weight variable used in calculating distribution between
C             inorganic and biologically bound nutrients:
      INTEGER MAXFYT_GROUPS
      PARAMETER ( MAXFYT_GROUPS = 2 )
      real*8 SUMAREA, WPO4, WNO3, WNH4, WSIO2
      real*8 WPFYT  (MAXFYT_GROUPS)
      real*8 WNFYT  (MAXFYT_GROUPS)
      real*8 WCFYT  (MAXFYT_GROUPS)
      real*8 WCHL   (MAXFYT_GROUPS)
      real*8 WSFYT  ! Only first group
      real*8 WDOC, WCZOO, WBACT, WCDET, WNDET, WPDET, WSDET, WRDET
      COMMON /HYDREX_WEIGHTS/
     &     SUMAREA, WPO4, WNO3, WNH4, WSIO2,
     &     WPFYT, WNFYT, WCFYT, WCHL,
     &     WSFYT, WDOC, WCZOO, WBACT, WCDET, WNDET, WPDET, WSDET, WRDET

C             can be addressed as single array:
      INTEGER W_LENGTH
      PARAMETER ( W_LENGTH = 5 + 4*MAXFYT_GROUPS + 9 )
      real*8 W_ARRAY(W_LENGTH)
      EQUIVALENCE  (W_ARRAY(1), SUMAREA)


$if BOUNDARY_TEST>1
      CHARACTER W_NAMES(W_LENGTH)*8 /
     &     'SUMAREA', 'WPO4', 'WNO3', 'WNH4', 'WSIO2',
     &     'WPFYT1', 'WPFYT2',
     &     'WNFYT1', 'WNFYT2',
     &     'WCFYT1', 'WCFYT2',
     &     'WCHL1',  'WCHL2' ,
     &     'WSFYT', 'WDOC', 'WCZOO', 'WBACT',
     &     'WCDET', 'WNDET', 'WPDET', 'WSDET', 'RWDET'/

C ...............................................................

      call Hello("HYDREX")

        IF ( EXTEST ) THEN
      write(debug_unit  ,*)'------input arguments received by HYDREX'
      write(debug_unit  ,*)'T:',       T
      write(debug_unit  ,*)'Air_pressure:', Air_pressure
      write(debug_unit  ,*)'WINDN:',   WINDN
      write(debug_unit  ,*)'WINDE:',   WINDE

      write(debug_unit  ,*)'NBSURF:',  NBSURF
      write(debug_unit  ,*)'NBE:',     NBE
      write(debug_unit  ,*)'INDXE:',   INDXE
      write(debug_unit  ,*)'MLE:',     MLE
      write(debug_unit  ,*)'NBI:',     NBI
      write(debug_unit  ,*)'INDXI:',   INDXI
      write(debug_unit  ,*)'MLI:',     MLI
      write(debug_unit  ,*)'ND:',      ND

      write(debug_unit  ,*)'DEPTH:',   DEPTH
      write(debug_unit  ,*)'AREA:',    AREA
      write(debug_unit  ,*)'TIDFAC',   TIDFAC
      write(debug_unit  ,*)'FIXTMP',   FIXTMP
      write(debug_unit  ,*)'TEMP:',    TEMP
      write(debug_unit  ,*)'PO4:',     PO4
      write(debug_unit  ,*)'NO3:',     NO3
      write(debug_unit  ,*)'SiO2:',    SiO2
      write(debug_unit  ,*)'NH4:',     NH4
      write(debug_unit  ,*)'EXTBIO:',  EXTBIO
      write(debug_unit  ,*)'MFYTG:',   MFYTG
      write(debug_unit  ,*)'FYTGRP:',  FYTGRP
      write(debug_unit  ,*)'CFYT:',    CFYT
      write(debug_unit  ,*)'NFYT:',    NFYT
      write(debug_unit  ,*)'PFYT:',    PFYT
      write(debug_unit  ,*)'SFYT:',    SFYT
      write(debug_unit  ,*)'CHL:',     CHL
      write(debug_unit  ,*)'DOC',      DOC
      write(debug_unit  ,*)'CZOO',     CZOO
      write(debug_unit  ,*)'BACT',     BACT
      write(debug_unit  ,*)'NCBACT',   NCBACT
      write(debug_unit  ,*)'PCBACT',   PCBACT
      write(debug_unit  ,*)'NCZOO',    NCZOO
      write(debug_unit  ,*)'PCZOO',    PCZOO
      write(debug_unit  ,*)'CDET',     CDET
      write(debug_unit  ,*)'NDET',     NDET
      write(debug_unit  ,*)'PDET',     PDET
      write(debug_unit  ,*)'SDET',     SDET
      write(debug_unit  ,*)'RDET',     RDET
      write(debug_unit  ,*)'ZSURFE',   ZSURFE
      write(debug_unit  ,*)'DZDTX',    DZDTX
      write(debug_unit  ,*)'EMIX_RELATIVE', EMIX_RELATIVE
          ENDIF
$endif

C ----------------------------------------------------------------
C  Only for calls from EUTRO.CSL for external basins:
      if (T. GE. 0.0) then
C  ---- Tidal surface level in external basins (one common level):
         CALL TIDAL( T, Air_Pressure, Z_EXT, DZ_DT_EXT,EMIX_RELATIVE )
         DO IB = 1, NBSURF
C ------- Tidal fluctuations Z_EXT and DZ_DT_EXT from TIDAL subroutine,
C         applied with factor specified as variable TIDFAC:
            ZSURFE(IB) = Z_EXT*TIDFAC
            DZDTX (IB) = DZ_DT_EXT*TIDFAC
         ENDDO
      endif

C ========= random variations in hydrographic conditions =========

      CALL RANDOM_SIGNALS(T, BNDR_T_LAST, dim_Signal, RESP_FREQ, 
     &   BNDR_IJ, BNDR_UC, ONE_TO_ONE, BNDR_RANDOM ) 
                           ! NO TRANSFORMATION OF PERTURBATIONS

$if defined  BOUNDARY_TEST
      if (EXTEST) THEN
         write(debug_unit,*)'Random signal in array BNDR_RANDOM:'
         write(debug_unit,'(1x,3(3(A,I3),A,G11.4))')
     &    ((('(', ID, ',',  IB, ',',  IP, ')=', BNDR_RANDOM(ID,IB,IP), 
     &            ID=1,2), IB=1,NBT), IP=0,NPAR)
      endif
$ENDIF

C   ...... get index for interpolating in tables DM, MV, STDEV:
      T_m = T*12.0/365.    ! Time in months with fractional part
      T_a = ANINT(T_m)-0.5 ! Previous point midway between whole months
      T_y = MOD(T_a+1.0,12.0D0)    ! values 0.5,  1.5, .. 11.5
C           indicating which values to use in tables for interpolating:
      IT1 = MOD(T_y+11.0,12.0D0)+1 ! values 12,    1,  .. 11
      IT2 = INT( T_y+1.01 )      ! values 1 ,    2,  .. 12
      
$if defined  BOUNDARY_TEST
      if (EXTEST) THEN
         write(debug_unit,'(4(A,":",G12.4),2(A,":",I6))') 
     &      'T',T, 'T_m',T_m, 'T_a',T_a, 'T_y',T_y, 'IT1',IT1, 'IT2',IT2

      endif
$endif

C  Interpolating weights:
      A1 = T_a+1 - T_m !Distance from current time to table values IT2
      A2 = 1-A1        !Distance from table values IT1 to current time


C ---------- calculate interpolated values for parameters,
C            with random component:

      DO IB = 1,NBT
         DO ID = 1,2

C        ...... Get parameter values into table of current values:
            DO IP = 1, NPAR
               CV(ID,IB,IP) = A1*MV(IT1,ID,IB,IP)+A2*MV(IT2,ID,IB,IP)
     &           + BNDR_RANDOM(ID,IB,IP)
     &             *( A1*STDEV(IT1,ID,IB,IP) + A2*STDEV(IT2,ID,IB,IP) )
               CV(ID,IB,IP) = MAX(0.0D0, CV(ID,IB,IP) )
               
$if defined  BOUNDARY_TEST
      if (EXTEST) write( debug_unit,
     &                   '('' CV('',I3,2('','',I3),'')='',G13.5)')
     &               ID,IB,IP,CV(ID,IB,IP)
$ENDIF
            
            ENDDO

C        ...... Depth values, with wind dependence:
            CD(ID, IB) = A1*DM(IT1,ID,IB) + A2*DM(IT2,ID,IB)
     &           + BNDR_RANDOM(ID,IB,0)*DDEV(ID)
     &           + DBETA(ID)*WINDN + 0.0*WINDE
C                              LAST TERM INCLUDED TO AVOID WARNING,
C                              NO EFFECT OF EAST/WEST WIND
         ENDDO

C           Ensure depths > 0:
         CD(1,IB) = MAX ( 0.1D0,  CD(1,IB) )

C           Transform CD(2,IB) from log10 of thickness into
C           depth of lower limit of pyknokline:
         R = max( -10.0D0, min (10.0D0,CD(2,ib) ) )
         CD(2,IB) = CD(1,IB) + 10.0**R
C   Now CD(1,IB) = depth at lower limit of surface layer,
C       CD(2,IB) = lower depth of pycnocline

$if defined  BOUNDARY_TEST
            IF ( EXTEST ) THEN
                write(debug_unit,'('' CD(I,'',I2,''),I=1,2:'',2g16.7)')
     &                      Ib, (cd(Id,ib),ID=1,2)
            ENDIF
$endif

      ENDDO

C    ....... in calls during simulation, after initiation through zhydr:
C        Set top layer temperature = temperature within model
C        if externally fixed temperatures are not specified in FIXTMP,
C        this means that surface temperature is determined by exchange
C        with atmosphere only, and assumed equal in internal and extern
C        basins:
      IF ( T.gt.0.0 .and. T_PREV.ge.0.0 .and. .not.FIXTMP ) then
         TEMP_SUM = 0.0
         AREA_SUM = 0.0
         DO IBI = 1, NBI
             LSURF = INDXI(IBI)+1
             W = max(1.e-20, BOUND_INFL(IBI))
             TEMP_SUM = TEMP_SUM + W*AREA(LSURF)*TEMP(LSURF)
             AREA_SUM = AREA_SUM + W*AREA(LSURF)
         ENDDO
         DO IB = 1, MIN( NBE, NBT)
             CTEMP (1,IB) = TEMP_SUM/AREA_SUM
         ENDDO
      ENDIF
      
      T_PREV = T

C    ........... limit bottom salinity within reasonable values:
C                below 34.8.0 0/00,

C                3. dec 1998: also keep temperature above certain
C                limit dependent on salinity, to avoid inflow of
C                too dense and cold water:
      DO IB = 1, MIN( NBE, NBT)

$if defined  BOUNDARY_TEST
            IF ( EXTEST ) THEN
               write(debug_unit,*)'1.  CSAL:',(CSAL(I,IB),I=1,2)
            ENDIF
$endif

         DO I=1,2
           CSAL(I,IB)  = MIN( CSAL(I,IB), 34.8D0 )
                          !  (SHOULD BE IMPROVED)
      
           CTEMP(I,IB) = MAX( CTEMP(I,IB),
     &                        MIN ( 4.75 + 0.7*(CSAL(I,IB)-33.5),
     &                              5.5  + 0.5*(CSAL(I,IB)-35  )  ) )  
      
      ! For Drammensfjorden:
      !      CSAL(I,IB) = CSAL(I,IB)*0.95  
      ! Reduserer salinitet p† randen med 5% (dvs. 30 -> 28.5)
  
         ENDDO


$if defined  BOUNDARY_TEST
            IF ( EXTEST ) THEN
               write(debug_unit,*)'2.  CSAL:',(CSAL(I,IB),I=1,2)
            ENDIF
$endif

      ENDDO


C    ........... adjust top layer salinity to avoid density inversion:
      DO LOOP = 1,100 ! ensure against infinite loop if error in SIGMAT
! ************* old code:
!         write(999,*) 'calls sigmat with CSAL, CTEMP from HYDREX'
         CALL SIGMAT ( CSAL, CTEMP, 2*NBT,  CDENS )
! ************* new code:
!         CALL SIGMAT ( CSAL, CTEMP, CTEMP, 0.0, 2*NBT,  CDENS )
                 ! 2. CTEMP is Dummy argument for particle concentration,
                 ! - effect on density is set to zero by next argument
                 ! water on the boundary is assumed not to have any
                 ! particle content.
         
         STABLE = .TRUE.
         DO IB = 1,NBT
            IF (CDENS(1,IB) .GT. CDENS(2,IB)) THEN
                CSAL(1,IB) = CSAL(1,IB) -1.0
                STABLE = .FALSE.
            ENDIF
         ENDDO
         IF (STABLE) EXIT
      ENDDO


C ------------------ insert values into external argument arrays:

    ! -------------------------------------------
      DO IB = 1,NBE  ! external basins
    ! -------------------------------------------

         IBT = MIN(NBT, IB)  ! ------ which twolayer table to use.

$if defined  BOUNDARY_TEST
         IF (EXTEST) THEN
           write(debug_unit,'('' ext. basin '',I3,'' use table '',i5)') 
     &         IB, IBT
         ENDIF
$endif


C ........... Adjust silicate values: surface layer value is valid
C             for salinity=0, and is adjusted according to a
C             formula based on surface data from inner Oslofjord
C             with minimum around S=25 o/oo:
         CSILICON(1,IBT) =
     &           MAX( 100.0D0,
     &                CSILICON(1,IBT)*(1-0.04*CSAL(1,IBT) ),
     &                MIN( CSILICON(2,IBT),17.*(CSAL(1,IBT)-25.0) ) )

       ! -----------------------------------------------------
         DO L = INDXE (IB)+1, INDXE (IB+1) ! external layers
       ! -----------------------------------------------------


             LD = L - INDXE(IB)

$if defined  BOUNDARY_TEST
                  IF (EXTEST) THEN
                      write(debug_unit,'('' L, LD:'',2i5)') l,LD
                  ENDIF
$endif

C --------- Interpolation into nonlinear profile:
             Z = (DEPTH(LD)+DEPTH(LD+1))/2.0 ! in middle of layer
             IF ( Z .LT. CD(2,IBT) ) THEN
                A1 = ( (CD(2,IBT)-Z) / (CD(2,IBT)-CD(1,IBT)) )**2
             ELSE
                A1 = 0.0
             ENDIF

             IF ( Z .GT. CD(1,IBT) ) THEN
                A2 = ( Z/CD(1,IBT) - 1 )**1.5
             ELSE
                A2 = 0.0
             ENDIF

$if defined  BOUNDARY_TEST
                  IF ( EXTEST ) THEN
                      write(debug_unit,
     &                   '('' Z + 1.phase A1, A2:'',3G17.7)') Z, A1, A2
                  ENDIF
$endif

             A1 = A1 / ( Z/CD(1,IBT) + A1**4 )**0.25
             A2 = A2 / ( 1 + A2**2)**0.5
$if defined  BOUNDARY_TEST
                  IF (EXTEST) THEN
                      write(debug_unit,
     &                   '(''     2.phase A1, A2:'',2G17.7)') A1, A2
                  ENDIF
$endif

C Interpolation coefficients depend on DEPTH relative to CD(1), CD(2):
C DEPTH: <CD(1)   =CD(1)   =CD(1)+CD(2))/2  =CD(2)   >CD(2)

C A1 and A2 combined into monotonously increasing weight A2:
C New A2:  0<--     0.1                      0.85   -->1
             A2 = (1-A1+A2)/2.0

C A1 set to complementary weight:
             A1 = 1 - A2

C Now A1 and A2 corresponds to (1-f) and f in the description
C given in reports on the model.


$if defined  BOUNDARY_TEST
                  IF (EXTEST) THEN
              write(debug_unit,'(''     final A1, A2:'',2G17.7)') A1, A2
                  ENDIF
$endif


C ------------ interpolate all variables to values for layer L:
             DO IP = 1, NPAR
                   V_L(IP) = A1*CV(1,IBT,IP) + A2*CV(2,IBT,IP)

$if defined  BOUNDARY_TEST
            IF ( EXTEST ) THEN
                write(debug_unit,'('' V_L('',I3,'')='',g16.7)')
     &                      IP, V_L(IP)
            ENDIF
$endif
             
             ENDDO  ! Can also be addressed as SAL_L, TEMP_L etc.


$if BOUNDARY_TEST>1
            IF (EXTEST) THEN
                write(debug_unit,'('' total values (unadjusted):'')')
                write(debug_unit,'(1x, A12,G12.5)') 
     &                   'TEMP_L',TEMP_L,
     &                   'SAL_L ',SAL_L ,
     &                   'OXYG_L',OXYG_L,
     &                   'TOTP_L',TOTP_L,
     &                   'TOTN_L',TOTN_L,
     &                   'TOTS_L',TOTS_L
            ENDIF
$endif


C ------------ transfer to external arrays:

C ......... Hydrophysical properties:
             SALEX (L) =  SAL_L
             TEMPEX(L) = TEMP_L

                                        ! Oxygen deficits adjusted by factor BndFac_Ox:
             OX_SAT_EX = OXYGEN_SATURATION ( SAL_L, TEMP_L )
             if (OXYG_L .lt. OX_SAT_EX) then
                 OXYG_L  = OX_SAT_EX + ( OXYG_L-OX_SAT_EX)*BndFac_Ox
             endif
                
             O2EX  (L) = OXYG_L

             IF (SAL_L .GT. 35.0 ) THEN
                WRITE(*,'('' SALEX('',i2,'')='',g12.5,'' csal='',2g12.5, 
     &                    '' ib,ibt='',2i3,'' A1,A2='',2g12.5,'' z='',
     &                    g12.5,''cd(1,ibt),cd(2,ibt)='',2g12.5)')
     &             L,SALEX(L), CSAL(1,IB),csal(2,iB),IB,IBT,a1,a2,z,
     &             CD(1,IBT), CD(2,IBT)
!                Pause 'Press Enter to continue'
             ENDIF

C ......... Nutrients:
C           distributed as close as possible to internal basins,
C           using area weighted mean concentrations.
C           Fytoplankton carbon and chlorofyll is set
C           relative to nitrogen, by analogy with
C           the basin that has lowest ratios.


!================================================================
! Total Nitrogen includes a seemingly inactive component of about
! 100 µg/l. This is subtracted from the value in order to get the
! biologically active part of nitrogen concentration:


             TOTN_L = MAX( 0.0D0, TOTN_L-100 )

!================================================================


! --------- Find area-weighted sums of concentrations
!           over internal basins ( should done once, and then
!                                  used for all external basins)
!     Zero accumulators:
             DO IW = 1, W_LENGTH
                W_ARRAY(IW) = 0.0
             ENDDO
             W_C_SUMFYT = 0.0

!     Sum over internal basins using weight BOUND_INFL:
             DO IBI = 1, NBI
C                  use data for layer closest in depth to external layer
                LI = MIN( INDXI(IBI) + (L-INDXE(IB)) , INDXI(IBI+1) )
                A = AREA(LI)*MAX(1.0e-20,BOUND_INFL(IBI)) 
                             ! Used to weight basins
                SUMAREA = SUMAREA + A

$if defined  BOUNDARY_TEST
                  IF (EXTEST) THEN
                    write(debug_unit,'(2(1x,A,''='',I6),A,''='',G12.5)') 
     &                      'IBI',IBI, 'LI', LI, 'A', A
                  ENDIF
$endif
C            ...... accumulate in values zeroed as elements in W_ARRAY
                WPO4  = WPO4  + (MAX(0.0D0,PO4(LI))+1.0e-16)*A
                WNO3  = WNO3  + (MAX(0.0D0,NO3(LI))+1.0e-16)*A
                WNH4  = WNH4  + MAX(0.0D0, NH4(LI))*A
                WSIO2 = WSIO2 + (MAX(0.0D0,SIO2(LI))+1.0e-16)*A
                 ! Addition of small value ensures against div. by zero

        ! Add biological components into weight, controlled by EXTBIO,
        ! if =0, no biomass is imported, =1: nutrients are distributed
        ! as inside model between inorganic and biomass forms.
                A = A * MAX ( 0.0, MIN( 1.0, EXTBIO(1) ) )
                if ( A .GT. 0.0) then
                   DO IG = 1, FYTGRP
                      WPFYT(IG) = WPFYT(IG) + MAX(0.0D0,PFYT(LI,IG))*A
                      WNFYT(IG) = WNFYT(IG) + MAX(0.0D0,NFYT(LI,IG))*A
                      WCFYT(IG) = WCFYT(IG) + MAX(0.0D0,CFYT(LI,IG))*A
                      WCHL (IG) = WCHL (IG) + MAX(0.0D0,CHL (LI,IG))*A
                      W_C_SUMFYT= W_C_SUMFYT+ MAX(0.0D0,CFYT(LI,IG))*A
                   
                   ENDDO
                   if ( FYTGRP.GE. 1 ) 
     &                    WSFYT = WSFYT + MAX(0.0D0,SFYT(LI))*A
                   WDOC  = WDOC  + MAX(0.0D0,DOC (LI))*A
                   WCZOO = WCZOO + MAX(0.0D0,CZOO(LI))*A
                   WBACT = WBACT + MAX(0.0D0,BACT(LI))*A
                      ! First, detritus:
                   WCDET = WCDET  + MAX(0.0D0,CDET(LI))*A
                   WNDET = WNDET  + MAX(0.0D0,NDET(LI))*A
                   WPDET = WPDET  + MAX(0.0D0,PDET(LI))*A
                   WSDET = WSDET  + MAX(0.0D0,SDET(LI))*A
                   WRDET = WRDET  + MAX(0.0D0,RDET(LI))*A
                ENDIF
             ENDDO

$if BOUNDARY_TEST>1
      IF (EXTEST) WRITE ( 99,'(A/3(1X,A8,'':'',G15.7))')
     &               ' in subroutine BOUNDARY: ',
     &               (W_NAMES(IW), W_ARRAY(IW), IW = 1,W_LENGTH)
$ENDIF

! --------- Compute total sum of weights for N,P and Si:
             WPSUM = WPO4        + WPDET + WCZOO*PCZOO + WBACT*PCBACT
             WNSUM = WNO3 + WNH4 + WNDET + WCZOO*NCZOO + WBACT*NCBACT
             WSSUM = WSIO2       + WSDET + WSFYT

             DO IG = 1, FYTGRP
                WPSUM = WPSUM + WPFYT(IG)
                WNSUM = WNSUM + WNFYT(IG)
             ENDDO

$if BOUNDARY_TEST>1
             if (EXTEST) WRITE ( 99,'(3(1X,A8,'':'',G15.7))')
     &               'WPSUM', WPSUM, 'WNSUM',WNSUM, 'WSSUM', WSSUM
$ENDIF

             C_SUMFYT_EXT = 0.0

             ! ------ transform total nutrients into measure
             !        relative to weighted internal sum,
             !        N and P also adjusted by 
             !        user-specified scenario factors:
             
             TOTN_L = (TOTN_L/WNSUM)*BndFac_N
             TOTP_L = (TOTP_L/WPSUM)*BndFac_P
             TOTS_L = TOTS_L/WSSUM


             ! ------- use transformed values to calculate external
             !         concentrations giving the same distribution
             !         between different forms:
             DO IG = 1, MFYTG
                NFYTEX(L,IG) = TOTN_L * WNFYT(IG)
                        ! Nonconservative components is set
                        ! in same proportion to corresponding N:

                CFYTEX(L,IG) = TOTN_L * WCFYT(IG)
                C_SUMFYT_EXT = C_SUMFYT_EXT + CFYTEX(L,IG)
                CHLEX (L,IG) = TOTN_L * WCHL(IG)
                PFYTEX(L,IG) = TOTP_L * WPFYT(IG)
             ENDDO

             NO3EX  (L) = TOTN_L * WNO3
             NH4EX  (L) = TOTN_L * WNH4



             CZOOEX (L) = TOTN_L * WCZOO
             BACTEX (L) = TOTN_L * WBACT
             NDETEX (L) = TOTN_L * WNDET
             CDETEX (L) = TOTN_L * WCDET ! (N:C proportional)
             RDETEX (L) = TOTN_L * WRDET
             PDETEX (L) = TOTP_L * WPDET
             PO4EX  (L) = TOTP_L * WPO4

             SIO2EX (L) = TOTS_L * WSIO2
             SFYTEX (L) = TOTS_L * WSFYT
             SDETEX (L) = TOTS_L * WSDET

             
             ODMEX  (L) = 0.0 ! always zero oxygen demand from outside
             
             
             if (W_C_SUMFYT.gt.0.0) then
                DOCEX (L) =  MIN( 1.0, MAX( 0.0,EXTBIO(2) ) )
     &                      *WDOC * C_SUMFYT_EXT / W_C_SUMFYT
             else
                DOCEX (L) =  0.0
             endif
      ! --------------------------------------
        ENDDO ! next layer
      ! --------------------------------------


    ! ----------------------------------------
      ENDDO ! Next basin
    ! ----------------------------------------


$if BOUNDARY_TEST > 1
      IF (EXTEST) THEN
	      write(debug_unit  ,*)'------output arguments exported by HYDREX'
	      write(debug_unit  ,*)'SALEX:',   SALEX
	      write(debug_unit  ,*)'TEMPEX:',  TEMPEX
	      write(debug_unit  ,*)'PO4EX:',   PO4EX
	      write(debug_unit  ,*)'NO3EX:',   NO3EX
	      write(debug_unit  ,*)'NH4EX:',   NH4EX
	      write(debug_unit  ,*)'SiO2EX:',  SiO2EX
	      write(debug_unit  ,*)'O2EX:',    O2EX
	      write(debug_unit  ,*)'CFYTEX:',  CFYTEX
	      write(debug_unit  ,*)'NFYTEX:',  NFYTEX
	      write(debug_unit  ,*)'PFYTEX:',  PFYTEX
	      write(debug_unit  ,*)'SFYTEX:',  SFYTEX
	      write(debug_unit  ,*)'CHLEX:',   CHLEX
	      write(debug_unit  ,*)'ZSURFE:',  ZSURFE
	      write(debug_unit  ,*)'DZDTX:',   DZDTX
	      write(debug_unit  ,*)'CZOOEX:',  CZOOEX
	      write(debug_unit  ,*)'BACTEX:',  BACTEX
	      write(debug_unit  ,*)'NDETEX:',  NDETEX
	      write(debug_unit  ,*)'CDETEX:',  CDETEX
	      write(debug_unit  ,*)'PDETEX:',  PDETEX
	      write(debug_unit  ,*)'SDETEX:',  SDETEX
	      write(debug_unit  ,*)'RDETEX:',  RDETEX
	      write(debug_unit  ,*)'ODMEX: ' ,  ODMEX
	      write(debug_unit  ,*)'DOCEX:',   DOCEX
      ENDIF
$endif

      END Subroutine


      end Module

