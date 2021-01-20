Module ModelVar_HydroBioChem

! Declaration of model state variables,
! and a subroutine for printout, binary save and binary restore


      use ModelDimensions
      use ModelVar_Topography

     implicit none


! ========================= Volumes:

     real*8 VTZERO (dimMBI)
     real*8 DVTOTL (dimMBI)
     real*8 VBUFZ (2,dimMLC), VDYNZ (2, dimMBI )
     real*8 VBFSZ (dimMBIplus1),  DVT2Z(dimMBIplus1)
     real*8 VTDZ (dimMBIplus1)

! ================= Physical state variables:

     real*8 SALTB(dimMBI), HEATB(dimMBI), OXYGB(dimMBI)
     real*8 NITRB(dimMBI), PHOSB(dimMBI)
     real*8 SILIB(dimMBI), SPPB(dimMBI)

     real*8 SAL( dimMLI ), TEMP( dimMLI ), OXYG (dimMLI)
        ! SAL=o/oo, TEMP = degC  OXYG ml/l

! ================= Nutrients:

     real*8 PO4(dimMLI), NO3(dimMLI), NH4(dimMLI), SiO2(dimMLI)

! ================= Phytoplankton:

     real*8 CFYT(dimMLI,dimMFYTG), PFYT(dimMLI,dimMFYTG)
     real*8 NFYT(dimMLI,dimMFYTG), CHL(dimMLI,dimMFYTG)
     real*8 SFYT(dimMLI)

     real*8 ODM (dimMLI)   ! Oxygen demand
     real*8 DOC(dimMLI)    ! Dissolved organic carbon
     real*8 BACT(dimMLI)   ! Bacterial carbon   (N and P implicit)
     real*8 CZOO(dimMLI)   ! Zooplankton carbon (N and P implicit)

     real*8 TOTC (dimMLI), TOTN (dimMLI), TOTP (dimMLI)

! ================= Detritus concentrations:

     real*8 CDET(dimMLI), NDET(dimMLI), PDET(dimMLI), SDET(dimMLI)
     real*8 RDET(dimMLI), RRDET(dimMLI)

! ================= Particulate flux downwards (mg/m2/day):
    real*8 CDFLUX (dimMLI), NDFLUX (dimMLI), PDFLUX (dimMLI)
    real*8 SDFLUX(dimMLI)

! ================= Integrated values (mg/m2)
    real*8 CDFLXI (dimMLI), NDFLXI (dimMLI), PDFLXI (dimMLI)
    real*8 SDFLXI(dimMLI)

! ============= Bottom sediment concentrations =================

     real*8 CSED(dimMLI), NSED(dimMLI), PSED(dimMLI), SSED(dimMLI)
     real*8 RSED(dimMLI), RRSED(dimMLI)

     real*8 PADS(dimMLI)  ! Phosphorus adsorbed to sediment
     real*8 ASED(dimMLI)  ! Oxygen debt in sediments

      real*8 XSED(dimMLI), XBUR(dimMBI)

! ================= radiation and wind friction integration terms
    real*8 QOUTI (dimMBI,5), QABSI(2), QGI
    real*8 QOUTZ (dimMBI,5), QABSZ(2)
    real*8 FR3NTZ (dimMBI), UFR3I(dimMBI)
    real*8 WN2I, WE2I  ! integration of 2. order wind components

! ================= Inorganic particles, added March 2001 BBJ

    real*8 SPP     (dimMLI)
    real*8 SPPFLUX (dimMLI)
    real*8 SPPFLXI (dimMLI)
    real*8 SPPSED  (dimMLI)

    integer FYTGRP        ! number of active fytoplankton groups

! ====================== mussel population

! ---------- Detailed state variables:
      real*8 MUSLNR   ( dimMUSLAGES, dimMUSLAYERS, dimMBI )
      real*8 MUSLMT   ( dimMUSLAGES, dimMUSLAYERS, dimMBI )
      real*8 MUSLMA   ( dimMUSLAGES, dimMUSLAYERS, dimMBI )
      real*8 MUSLGT   ( dimMUSLAGES, dimMUSLAYERS, dimMBI )
      real*8 MUSLWM                 ( dimMUSLAYERS, dimMBI )
      real*8 MUSLBT                 ( dimMUSLAYERS, dimMBI )

! ----------- Summed over layers within basins:

      real*8 CTMUSL(dimMBI)  ! Total weight of mussels mgC
      real*8 CAMUSL(dimMBI)  ! Total active weight     mgC
                                 ! sum over ages
      real*8 C0MUSL(dimMBI)  ! Weight of age class 0
      real*8 CMUSDV(dimMBI)  ! Derivative of CTMUSL,
                           ! used for mass balance of derivatives 
                           ! in EUTROSUB mgC/day

  !  For each age class, sum/mean over all basins and layers:
      real*8 MUSLWT (dimMUSLAGES )  ! Mean total weight of individuals
      real*8 MUSLWA (dimMUSLAGES )  ! Mean total weight of individuals

  !  Actual number of layers with mussels:
      INTEGER MSNLAY


        !  Permanent export to sediment as mg/day and mg integrated
      real*8 CSEDXP (dimMBI), CSEDXI (dimMBI)
      real*8 NSEDXP (dimMBI), NSEDXI (dimMBI)
      real*8 PSEDXP (dimMBI), PSEDXI (dimMBI)
      real*8 ASEDXP (dimMBI)
      real*8 SSEDXP (dimMBI)

      ! Nitrate removed by denitrification
      real*8 DENITR (dimMBI)  ! (mg/day)
      real*8 DNITRI (dimMBI)  ! mg: Integrated sum


      ! Nitrogen fixation , controlled by parameter NFIXRR:
      real*8 NFIX (dimMBI), NFIXI(dimMBI)
            ! Total, as mg/day and integrated mg

                     !  C, ODM, N & P from land, flux and integral:
      real*8 CLOAD   (dimMBI), CLOADI  (dimMBI)
      real*8 ODMLOAD (dimMBI), ODMLOADI(dimMBI)
      real*8 NLOAD   (dimMBI), NLOADI  (dimMBI)
      real*8 PLOAD   (dimMBI), PLOADI  (dimMBI)



!  import rates xxxxMP defined within this section,
!  integrated to values xxxxMI
!  used in subroutine COMCLC to check mass conservation
!  xxxxMI(..,1) contains total integral,
!  since start of simulation, xxxMI(...,2) contains values
!  at last restart of flux integrals, see WRTPRT and CNCADJ

      real*8 SALTMI (dimMBI,2), SALTMP( dimMBI)
      real*8 HEATMI (dimMBI,2), HEATMP (dimMBI)
      real*8 OXYGMI (dimMBI,2), OXYGMP (dimMBI)
      real*8 NITRMI (dimMBI,2), NITRMP (dimMBI)
      real*8 PHOSMI (dimMBI,2), PHOSMP (dimMBI)
      real*8 SILIMI (dimMBI,2), SILIMP (dimMBI)
      real*8 SPPMI  (dimMBI,2), SPPMP  (dimMBI)



!  ====================== mass and heat conservation ===================
      real*8 SALTMZ (dimMBI)
      real*8 HEATMZ (dimMBI)
      real*8 OXYGMZ (dimMBI)
      real*8 NITRMZ (dimMBI)
      real*8 PHOSMZ (dimMBI)
      real*8 SILIMZ (dimMBI)
      real*8 SPPMZ  (dimMBI)

          ! xxxxMZ: 1. Size scale, accumulated throughout simulation,
          !         2. Base value of amount xxxxMI


  !  ---- unity concentrations for continuity check on volume 
  !       or for residence time of water within system
      real*8 C1EX ( dimMLE ), C1 ( dimMLI )


      real*8 XMIX (dimMBI)   ! Number of well mixed surface layers,
                         ! with fractional part, set by CNCADJ,
                         ! used by TRANSP. Recorded values are not
                         ! representative, because the value depend
                         ! on step length: small values for final
                         ! small step. value/timestep is relevant
                         ! to indicate mixing

      integer LMIX(dimMBI) ! index for lowest layer which must be
                         ! homogenized due to non-linear density,
                         ! only used in CNCADJ, available as info here.

 ! ================= fra EULER Section ===================================


               ! ----- Density as a function of salinity and temperature:
      real*8 DENSI ( dimMLI )

      real*8 VLAYER(dimMLI)  ! Volume of basin layers
      real*8 VLCORF(dimMBI)
               ! Correction factor for open part of volume,
               ! computed by VLCALC called by CNCADJ in EUTROSUB,
               ! and used in TRANSP


!  ----- hydrographic/biochemical conditions of external basins ----:
      real*8 SALEX ( dimMLE ), TEMPEX ( dimMLE )
      real*8 DENSEX ( dimMLE ), OXYGEX(dimMLE)
      ! Units : SALEX = o/oo, TEMPEX = temp, DENSEX = sigma-t
      real*8 PO4EX(dimMLE), NO3EX(dimMLE), NH4EX(dimMLE),  &
         SIO2EX(dimMLE)
      !         PO4, NO3, NH4, SiO2 : µg Si/l
      real*8 CFYTEX(dimMLE,dimMFYTG)
      real*8 PFYTEX(dimMLE,dimMFYTG), NFYTEX(dimMLE,dimMFYTG)
      real*8 SFYTEX(dimMLE), CHLEX(dimMLE,dimMFYTG)

      real*8 ZSURFE ( dimMBE ), DZDTX( dimMBE )
     ! Surface level (m) with time derivative (m/day) of ext. basins.

      real*8 EMIXRL  ! Input of energy with semidiurnal tides
                  ! relative to mean value over time, used to control
                  ! mixing energy


  ! ---- Total volume of basins, and calculates water
  !      transports as local variables, later used in MTRANS
      real*8 VDYN ( 2, dimMBI) , DVTOT( dimMBI)
      real*8 VTOTDV ( dimMBI ), VDYNDV ( 2, dimMBI )
      real*8 VBUF (2, dimMLC), VBUFDV(2,dimMLC)
 
   ! ------------- Check on volume derivative balance
        ! Value nr. NBI+1 for external basins
      real*8 VBFSDV (dimMBIplus1), DVT2DV(dimMBIplus1),  &
        VTDDV (dimMBIplus1)

  !  integrated values, should compare with VBFSUM, DVTOT2 and VTDIFF
      real*8 VBFSI (dimMBIplus1), DVT2I(dimMBIplus1),  &
         VTDI (dimMBIplus1)

      LOGICAL TRCALC  ! = status from TRANSP call,
                       !  used by TRNADJ and MTRANS


   ! Buoyancy influx to surface from other basins
      real*8 BSFLUX (dimMBI)

   ! Flow velocities through connections at upper limit of layers:
      real*8 UFLOW (dimMLC)  ! Returned as information from TRANSP,
                           ! not used except as monitoring output

      real*8 BWFreq (dimMLI)
           ! Stability as Brunt_W,isel, frequency (1/S)
      real*8 ZMID  (dimMLI)  ! Mid-depth of layers, volume mean

      LOGICAL VTRNEG  ! .TRUE. signals negative vertical transports
         ! due to advection dominating over diffusion between some
         ! model layers ( see subr. TRANV1 in TRANSP_3.FOR)
         ! May suggest using a finer vertical resolution
      ! Final time step TSTEP has been set by LGTRAD <= MAXTTR
      !  VTRNEG is also used to sort MTRANS after TRNADJ



   ! --------- VOLUME CONTINUITY: -----------
      real*8 C1DV(dimMLI), C1MP( dimMBI )
   ! ---- Salinity SAL:
      real*8 SALDV( dimMLI )
   ! ---- Temperature TEMP
      real*8 TEMPDV( dimMLI )
   ! ---- OXYGEN CONTENT   (Oxygen unit ml/l)
      real*8 OXYGDV( dimMLI )
      real*8 OXSAT(dimMBI), OXEXCF
           ! Oxygen saturation & exchange coeff.


!  ---------------------- nutrients and oxygen ----------------
  ! -------- Input from land of water and substances:
  !          table of values contained in subroutine INFLUX,
  !          which returns inflow as functions of time:

  !  combined with scale factors and transfer coefficient declared here


      real*8    QWATER(dimMS,2)  ! Amount of water (m3/s)
                            ! in freshwater outlet          (k,1)
                            ! and intake of recipient water (k,2)

      real*8    QPO4(dimMS), QNO3(dimMS), QNH4(dimMS)
                  ! Nutrients, as kg/d

      real*8    QODM(dimMS)    ! Oxygen demand as kg/d
      real*8    QDOC(dimMS)   ! Dissolved organic carbon as kg/d
      real*8    QCDET(dimMS), QNDET(dimMS), QPDET(dimMS)
                          ! Particulate matter generated from
                          ! freshwater and sewage input
      real*8    QO2(dimMS)     ! Oxygen as kg/d
      real*8    QSIO2(dimMS)   ! Silicate as kg/d
      real*8    QTEMP(dimMS)   ! Heat as degC*m3/d


      LOGICAL QTR_ACTIVE(dimMS) ! true/false Set by subroutine below
      real*8 QTR_TimeOff (dimMS)  ! Possible switchoff time, set by subroutine


      INTEGER NS      ! Number of outlets (sources), used in EUTROSUB.FOR

      INTEGER RNFNDX(dimMS,2)
               ! Global index of outlet layer (k,1)
               ! and intake layer for QMIXIN  (k,2) ...

      INTEGER RCQNDX(dimMS)
               ! Global index of final receiving layer
               ! set by TRANSP, and used by QCALC
      INTEGER NQDIST(dimMS)
               ! Number of layers receiving runoff
               ! set by TRANSP, and used by QCALC

      real*8 Ambient_Volume_Flux (dimMS) ! volume mixed into dived outlets
      real*8 Neutral_Depth       (dimMS) ! neutral depth for dived outlets
            

  ! >>>>>>>>>>>>>>>>>> Ide til forbedring: <<<<<<<<<<<<<<<<<<<<
  ! Her i modellen spesifiseres et antall utslipp som nå, men
  ! antallet (innenfor dimMS) og bassengnr. spesifiseres også her.
  ! i tillegg settes opp en liste over hvilken kilde utslippet
  ! skal komme fra. Kildene er spesifisert som nå i RUNOFF,
  ! QFW kan brukes til b†de † fordele og ›ke minske kildene,
  ! dvs. summen av QFW for utslipp fra en kilde er total
  ! faktor for kilden. Overføring mellomutslippene
  ! kan kobles ved overf›ringskoeffisienter som nå,
  ! tomme utslipp blir i alle fall ikke beregnet i TRANSP_rutinene


      real*8 QWSURF(dimMBI)
      ! Total surface water influx from atmosphere.

        ! %%%%%%%%%%%%%% particle influx %%%%%%%%%%%%%%%%%%%%

      real*8 QSPP (dimMS)


           ! Nutrient derivatives:

   !  ORTO-PHOSPHATE
      real*8 PO4DV ( dimMLI ), PO4MP ( dimMBI )
   !  NITRATE
      real*8 NO3DV ( dimMLI ), NO3MP ( dimMBI )
   !  AMMONIUM
      real*8 NH4DV ( dimMLI ), NH4MP ( dimMBI )
   !  SILICATE
      real*8 SiO2DV( dimMLI ), SiO2MP( dimMBI )



  !  --- Phytoplankton, described by mass substances C, N and P,
  !      and chlorofyll, treated as property, not as mass substance.
  !    Chlorofyll has unit (C-conc)/day/(W/m2) and really is
  !    a measure of the quantum yield
      real*8 CFYTDV ( dimMLI,dimMFYTG ), CFYTMP ( dimMBI )

      real*8 NFYTDV ( dimMLI,dimMFYTG ), NFYTMP ( dimMBI )
      real*8 PFYTDV ( dimMLI,dimMFYTG ), PFYTMP ( dimMBI )
  !  Silicon: only for group one in phytoplankton
      real*8 SFYTDV ( dimMLI), SFYTMP ( dimMBI )
      real*8 CHLDV ( dimMLI,dimMFYTG ), CHLMP ( dimMBI )

!  ------------- OXYGEN DEMAND -----------
      real*8 ODMDev (dimMLI), ODMImp (dimMBI),  &
         ODMExt(dimMLE)

  !  ------------- DISSOLVED ORGANIC CARBON -----------
      real*8 DOCDV (dimMLI), DOCMP (dimMBI), DOCEX(dimMLE)

  !  ------------- BACTERIAL CARBON -------------------
      real*8 BACTDV (dimMLI), BACTMP (dimMBI), BACTEX(dimMLE)

  !  ------------- ZOOPLANKTON ---------------
      real*8 CZOODV (dimMLI), CZOOMP (dimMBI), CZOOEX(dimMLE)

  !  -------------- DETRITUS -----------------

      real*8 CDETDV ( dimMLI ), CDETMP ( dimMBI ), CDETEX ( dimMLE )
      real*8 NDETDV ( dimMLI ), NDETMP ( dimMBI ), NDETEX ( dimMLE )
      real*8 PDETDV ( dimMLI ), PDETMP ( dimMBI ), PDETEX ( dimMLE )
      real*8 SDETDV ( dimMLI ), SDETMP ( dimMBI ), SDETEX ( dimMLE )
      real*8 RDETDV ( dimMLI ), RDETMP ( dimMBI ), RDETEX ( dimMLE )

  !  --------- SEDIMENTARY MATTER ------------
      real*8 CSEDDV ( dimMLI )  ! Carbon
      real*8 NSEDDV ( dimMLI )  ! Nitrogen
      real*8 PSEDDV ( dimMLI )  ! Phosphorus
      real*8 SSEDDV ( dimMLI )  ! Silicon
      real*8 RSEDDV ( dimMLI )  ! Carbon * degradation rate
      real*8 ASEDDV ( dimMLI )  ! Stored oxygen debt as sulphide
      real*8 PADSDV ( dimMLI )  ! Phosphorus adsorbed
      real*8 XSEDDV (dimMLI), XBURDV(dimMBI)  ! residual X values


  !  --------- Sum of organic particulate matter ------------
      real*8 PartC (dimMLI), PartN (dimMLI), PartP (dimMLI)
 !  Sum of phytoplankton, zooplankton and sedimenting matter,
 !  see PRPROD. PartC is used to calculate light attenuation
 !  by subroutine RADABS in Module SURFEXCH.


   !  ------------------- MONITORING OF PRIMARY PRODUCTION -------------
      real*8 GTEMP(dimMLI,dimMFYTG), GTN(dimMLI,dimMFYTG),  &
          GRATE(dimMLI,dimMFYTG), CHLREL(dimMLI,dimMFYTG),  &
          SEDVF(dimMLI), FRESP(dimMLI,dimMFYTG)



! ------------------ Mussels ---------------------

      real*8 MVRED (dimMBI)
      !  exported from MUSLDV: reduction in effective concentration
      !  before shift_down of filtering (see MUSLINTG.FOR)


  ! Inorganic particles -  Added March 2001 - BBJ
      real*8 SPPDV ( dimMLI ), SPPEX ( dimMLE ), SPPSEDDV( dimMLI)


!  ---------- Other physical variables at time T -----------


!  ---------------------- CLIMATIC VARIABLES --------------------------
      real*8 AIRTMP, AIRP, HUMGM3, CLOUDS, WINDN, WINDE, PRECIP, WNDSPD
      !  deg.Cels  mb    g/m3    0-8     m/s    m/s     m/s     m/s

      real*8 FR3INT(dimMBI)
         ! Accumulated wind mix energy/area as Ufric**3
      real*8 UFRIC3(dimMBI)
         ! Wind frictional velocity in 3rd exponent


  !  ------------- solar and diffuse light radiation -------------
      real*8 SINHS, SUNRAD, DIFRAD


  !  Heat and evaporative water exchange across water/air boundary:
  !  also calculates wind friction coefficient:
  ! ################################################################


     real*8 QOUT (dimMBI,5)   ! Heat loss in W/m2

      ! QOUT(I,K),K=1..5  contains for basin I:
      !   K=1: Net heat loss to atmosphere,
      !        used to get surface temp. derivative
      !   for test purposes:
      !   K=2: Incoming longwave radiation
      !   K=3: Outgoing longwave radiation
      !   K=4: Heat loss by evaporation
      !   K=5: Heat loss by conduction

     real*8 EVAP  ( dimMBI )  ! Evaporation in m/s

  ! ----- Radiation and heat balance as a function of depth:
     real*8 RAD ( dimMLI )
         ! Light radiation in middle of each layer W/m2
  ! ################################################################
  
     real*8 AIRTI  ! integrated air temperature (not used?)

     real*8  QABS(2)  ! Absorbed radiation (W/m2)
              ! 1: IR radiation absorbed in surface layer
              ! 2: Visual radiation absorbed at depth

 !  ------- Buoyancy energy balance:
     real*8 BFX (dimMBI)  ! Relative buoyancy influx at surface (m2/s3)
     real*8 SPP_ToSurf (dimMBI)  ! Scratch variables = influx of particles
     real*8 QTemp_ToSurf(dimMBI) ! and temp*volume to basin surface layer
                               ! from land
                               ! Only used locally in SURFBF
     real*8 BFXINT (dimMBI), BFXZ(dimMBI)


  ! ------------- Check on volume balance --------------
     real*8 VBFSUM (dimMBIplus1), DVTOT2(dimMBIplus1),  &
         VTDIFF (dimMBIplus1)
        ! Value nr. NBI+1 for external basins


 !  Surface levels for comparison with external level ZSURFE
     real*8 ZSURFI ( dimMBI )
          ! = Vdyn/Area: integrated dynamic transport
          !   translated to change of surface levels

   !  Radiation terms in W/2
     real*8 QOUTD (dimMBI,5), QABSD(2), QGD

    !  Phytoplankton monitoring variables
     real*8 CHLCF (dimMLI,dimMFYTG), NCFYT(dimMLI,dimMFYTG),  &
          PCFYT (dimMLI,dimMFYTG), CFDVR(dimMLI,dimMFYTG),  &
          NUTLIM(dimMLI,dimMFYTG)


end Module ModelVar_HydroBioChem