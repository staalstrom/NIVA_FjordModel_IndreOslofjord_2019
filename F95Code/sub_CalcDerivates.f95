Module sub_Integration

    use ModelParam_RunControl
    use ModelParam_InitState
    use ModelParam_Inputs
    use ModelParam_Plankton
    use ModelParam_Physics
    use ModelParam_Boundary
    use ModelVar_RunControl
    use ModelVar_HydroBioChem
    use ModelVar_Topography
    use fx_Decomposition
    use fx_Eutrosub
    use fx_DVZero
    use fx_SurfExch
    use fx_Transp_2
    use fx_Mtrans
    use fx_Phyt_Zoo
    use fx_Boundary
    use fx_Runoff
    use fx_Water_Transports
    
    
    implicit none

contains

   !----------------------------------------------------------------------------
   ! Reorganised code from the discrete Euler section of the old ACSL model code
   !----------------------------------------------------------------------------

!=================================================================================

   Subroutine Integrate_Step

   
      CALL HELLO( 'Integrate_Step (Euler)')
 
      CALL DEBUGF(1)           ! Control debug dump files, and count down debugging steps
   
            ! Updates integrated variables for water and sediment:
            ! mixing due to wind and density instability
            ! and adjust concentrations accordingly.   
            !  (Requires DVTOT, VTOTZ to be set.)
            ! Density is set using adjusted values

      CALL HELLO ('First call to CNCADJ')

      CALL CNCADJ(dimMBI,dimMLI, dimMUSLAGES, dimMUSLAYERS )
            ! Integrate state values with derivatives from previous call to Euler_Step.
            ! The same subroutine is called by the DYNAMIC section,
            ! performed just before each output

            
   end Subroutine


!=================================================================================
   
    Subroutine Calc_Derivatives      

      integer indx_MS

      CALL HELLO( 'Calc_Derivatives')
 

            
! dropped:      if ( LOG_ALL_STEPS ) CALL LOGD(.false.)
! ACSL-specific code for printing variables
 
      YEARS = TINTEG/365.0   ! Current point in time that model is integrated to (days)

         ! Adjust maximum time step not to go beyond next communication time
         ! or final time, but always as a significant increment of current T
            ! to avoid problems with eternal loop.
         !  SHOULD BE IMPROVED - DOES NOT ALWAYS WORK)
         !  from suroutine STPADJ in ACSL/WATCOM version

      MAXINT = STPLIM ! user defined timestep limit

                      ! Limit to reach next communication time,
                      ! or just before stopping time:
      MAXINT = MIN( MAXINT, (TTERM-TTERM*1.e-6)-T, TCOMM+CINTV-T )
                      ! Ensure significant time step:
      MAXINT = MAX( T*1.e-6, MAXINT )

      MAXINT = MIN (MAXINT, MXTBIO)

         ! Zero derivatives in all layers and the import sum variables for each basin:
         ! (can be replaced by simple array assignment to scalar zero in F95 code)

      CALL DVZERO( 1, NLI,C1DV  , C1MP   )
         !  Derivative of C1 due to water transports ( C1EX =C1XTRN)
      CALL DVZERO( 1, NLI,SALDV, SALTMP )
      CALL DVZERO( 1, NLI,SALDV, SALTMP )
      CALL DVZERO( 1, NLI,TEMPDV, HEATMP  )
      CALL DVZERO( 1, NLI, OXYGDV, OXYGMP  )

         ! Check Seasonal and Oxygen content conditions for discharge transfer,
         ! end turn flags QTR_Active on or off and sets QTR_TimeOff.
         ! NOTE1: Must be set before discharges are calculated.
         ! (CHECK_Q_TRANSFER is defined in runoff module)

      CALL CHECK_Q_TRANSFER( dimMS, QTR_OXCOND, BASINQ, NLI, OXYG, NBI, INDXI, &
                               T, QTR_CheckInterval, QTR_MinActiveInterval, &
                        QTR_ACTIVE, QTR_TimeOff )
      NS = dimMS
 
      CALL DVZERO( 1,NLI,PO4DV,PO4MP)   
      CALL DVZERO( 1,NLI,NO3DV,NO3MP)   
      CALL DVZERO( 1,NLI,NH4DV,NH4MP)   
      CALL DVZERO( 1,NLI,SIO2DV,SIO2MP)       
      CALL DVZERO( 2,NLI,CFYTDV,CFYTMP) 
      CALL DVZERO( 2,NLI,NFYTDV,NFYTMP) 
      CALL DVZERO( 2,NLI,PFYTDV,PFYTMP) 
      CALL DVZERO( 1,NLI,SFYTDV,SFYTMP) 
      CALL DVZERO( 2,NLI,CHLDV,CHLMP)   
      CALL DVZERO( 1,NLI,ODMDEV,ODMIMP) 
      CALL DVZERO( 1,NLI,DOCDV,DOCMP)   
      CALL DVZERO( 1,NLI,BACTDV,BACTMP) 
      CALL DVZERO( 1,NLI,CZOODV,CZOOMP) 
      CALL DVZERO( 1,NLI,CDETDV,CDETMP) 
      CALL DVZERO( 1,NLI,NDETDV,NDETMP) 
      CALL DVZERO( 1,NLI,PDETDV,PDETMP) 
      CALL DVZERO( 1,NLI,SDETDV,SDETMP) 
      CALL DVZERO( 1,NLI,RDETDV,RDETMP) 

         ! turns biological processes on or off depending on user specification
         ! But always active in the initial calculation to set variables
         
      BIOACT = (.not. BIOOFF) .or.  (T.le.TSTART)

         !  Degradation of organic matter and DOC in water column and bottom sediments.
         ! Sediment biomass derivatives initiated, nutrient derivatives updated:

      CALL DGRADE( (DGTEST .AND. T.ge.TTRIG), dimMLI)

          ! Get values of Climatic variables: 

      CALL METINP( T, AIRTMP, AIRP, HUMGM3, CLOUDS, PRECIP, WINDN, WINDE)
      WNDSPD = SQRT(WINDE**2+WINDN**2)

         ! Get boundary conditions (ZSURFE and subsequent arguments)
         ! Combining specifications of total concentrations
         ! with distribution patterns based on state of internal basins:

      CALL HYDREX(EXTEST, T, AIRP, WINDN, WINDE, NBE, NBE,  &
                INDXE, dimMLE, NBI, INDXI, dimMLI, ND, DEPTH,  &
                AREA, TIDFAC, FIXTMP, TEMP, PO4, NO3, NH4, SIO2, &
                BOUND_INFL, EXTBIO, dimMFYTG, FYTGRP,  &
                CFYT, NFYT, PFYT, CHL, SFYT, &
                DOC, CZOO, BACT, NCBACT, PCBACT, NCZOO, PCZOO, &
                CDET, NDET, PDET, SDET, RDET, &
                BndFac_N, BndFac_P, BndFac_Ox, &
            ZSURFE, DZDTX, EMIXRL,&
            SALEX, TEMPEX, PO4EX, NO3EX, NH4EX, SiO2EX, OXYGEX, &
            CFYTEX, NFYTEX, PFYTEX, CHLEX, SFYTEX, ODMExt, DOCEX, &
            CZOOEX, BACTEX, CDETEX, NDETEX, PDETEX, SDETEX, RDETEX)


         ! external density, with old code:
      CALL SIGMAT ( SALEX, TEMPEX, INDXE(NBE+1), DENSEX )
         ! new call not yet used:
!       CALL SIGMAT ( DENSEX = SALEX, TEMPEX, TEMPEX, 0.0, INDXE(NBE+1) )
              ! 2. TEMPEX in new call is Dummy argument for particle concentration,
              ! - effect on density is set to zero by next argument
              ! water on the boundary is assumed not to have any
              ! particle content.

           
          ! -------- Inputs from land of water and substances.
          !          Table of values contained in runoff subroutine,
          !          Inputs set by interpolation in table
          !          and applying user coefficients for changing inputs:
          !          Subroutine returns RNFNDX and subsequent arguments:

      CALL RUNOFF ( DBGRNF.AND. T.GE.TTRIG, T, NBI, INDXI, ND, DEPTH, AIRTMP,  &
                     dimMS, BASINQ, AIRTQF, QFW, QFP, QFN, QFODM, QFC, QFS, QFNH4,  &
                     QCDETF, QNDETF, QPDETF, QNCDET, QPCDET, QTRF, QTRNDX, QTR_Times, &
                DEPTHQ,  &
            QMIXIN,  &
            MIXQTM,  &
            QTR_ACTIVE,  &
            RANDFACTOR, &
            RNFNDX, QWATER, QTEMP, QPO4, QNO3, QNH4,  &
                QODM, QDOC, QCDET, QNDET, QPDET, QSiO2, QO2 )
         ! values of QPO4, QNO3, QNH4, QDOC, QO2 etc. required in kg/d


            ! Determine Oxygen saturation (OXSAT) 
         ! and Water/Atmosphere exchange coeff. (OXEXCF)
         ! as function of hydrophysical surface conditions
         ! and relevant model control parameters:

      CALL OXEXCH(WNDSPD, TEMP, SAL, AREA, &
               ND, DEPTH, NLI, INDXI, NBI, &
               VLAYER, OXYG, OXYGDV, OXYGMP, &
            OXSAT, OXEXCF )

         !  Calculate Heat and evaporative water exchange 
         !  across water/air boundary, and also wind friction coefficient:

      CALL SETR  (0.0, QOUT, dimMBI*5)
      CALL SETR  (0.0, UFRIC3,  dimMBI)
      CALL ENEXCH( WNDSPD, AIRTMP, HUMGM3, CLOUDS, &
                     CDFAC, CEFAC,  &
                     dimMBI, NBI, NLI, INDXI, TEMP, &
                     UFRIC3, EVAP, QOUT )
            
         ! Transform pr. area values for precipitation and evaporation
         ! into total water flux values (QWSURF) to be used by transport calculations

      CALL QWCALC( NBI, INDXI, NLI, PRECIP, EVAP, AREA, QWSURF )

         ! ----- Water transport calculated
         !       Details hidden in TRANSP subroutine and used in calls
         !       to entry MASS_TRANSPORT in module TRANSP_1 via MTRANS
         !
      
         ! Sets up water transports and related derivatives (UFLOW, VTOTDV, VDYNDV, VBUFDV)
         ! and descriptors (MAXTTR, NQDIST, RCQNDX, TRCALC, BSFLUX, BWFreq, ZMID)
          !    ( Wind induced flow across connections is set in subroutine
         !      GET_WIND_CURRENT in TRANSP_H.FOR, called via TRANSP)

         ! Subroutine TRANSP calls subroutine Water_Transports where arrays 
         ! and other arguments are transmitted to the code doing the calculations. 
            ! Logical parameter TROFF can be used to switch between real transport calculation
            ! and dummy call (may save time when studying biological processes by themselves)
            ! Mass transport subroutine MTRANS acts accordingly.

      CALL HELLO ('TRANSP'      )
      CALL TRANSP
      CALL HELLO ('Back from TRANSP'      )

      ! Transp also sets maximal time step MAXTTR as limited by water transports.
      MAXSTP = MIN ( MAXINT, MAXTTR )
         ! Obsolete comment - can now change code?
         ! Dependence on old MXTBIO value is hidden from ACSL
         ! The radiation in LGTRAD is computed in connection with
         ! the calculating of timestep, and primary production is
         ! then calculated based on that radiation. This may also
         ! limit the timestep further. The dependence on the old
         ! MXTBIO above is included to reduce the final time step
         ! and the time step used in calculating the radiation.
         ! In reality there is a circular dependence that strictly
         ! should be solved by iteration in each step

         
         ! ------ Calculate check on volume derivative balance

      CALL SETR (0.0, VBFSDV,  dimMBIplus1 )
      CALL SETR (0.0, DVT2DV,  dimMBIplus1 )
      CALL SETR (0.0, VTDDV,   dimMBIplus1 )
      CALL HELLO ('calls VOLCHK')
      CALL VOLCHK (1, NBI, NC, INDXC, BCONN1, BCONN2, NLC, VBUFDV,  &
        VDYNDV, VTOTDV, VPRT .AND. t.GE.TTRIG, vdindx, VBFSDV, DVT2DV, VTDDV )
      CALL HELLO ('returns from VOLCHK' )


                 ! >>>>>>>>>>>>>>>>>> Ide til forbedring: <<<<<<<<<<<<<<<<<<<<
                 ! Her i modellen spesifiseres et antall utslipp som nå, men
                 ! antallet (innenfor dim_ MS) og bassengnr. spesifiseres også her.
                 ! i tillegg settes opp en liste over hvilken kilde utslippet
                 ! skal komme fra. Kildene er spesifisert som nå i RUNOFF,
                 ! QFW kan brukes til både å fordele og øke/minske kildene,
                 ! dvs. summen av QFW for utslipp fra en kilde er total
                 ! faktor for kilden. Overføring mellomutslippene
                 ! kan kobles ved overføringskoeffisienter som nå,
                 ! tomme utslipp blir i allefall ikke beregnet i TRANSP_rutinene

                 
         ! If transports are not completely turned off by model parameter TROFF,
         ! Update derivative of unity concentration with volume input
         ! from land and through water/atmosphere surface, cfr. transport modules

      IF (.not.TROFF ) THEN
            !  Land runoff:
          CALL QCALC ( NQDIST, RCQNDX, C1DV,  &
                       QWATER, C1XTRN*24.*3600., 'C1'  , C1MP )
                               ! /s -> /day
            !  Surface water exchange:
          CALL QSCALC ( C1DV, C1MP, QWSURF, C1XTRN*24.*3600.  )
      END IF


         ! Updates biochemical time derivatives due to land runoff,
         ! at the same time updating net import value.
         ! The 5. argument is factor for converting from kg of substance to
         ! the mass unit per m3 used for water concentrations,

         ! Oxygen
      CALL QCALC ( NQDIST, RCQNDX, OXYGDV, QO2, 1.e3/1.429*RNF, 'O2', OXYGMP )
                                    ! kg--> l

         ! Specified explicit chemical oxygen demand:
      CALL SETR  ( 0.0, ODMLOAD, NBI )
      CALL QCALC ( NQDIST, RCQNDX, ODMDev, QODM, 1.e3/1.429*RNF, 'ODM', ODMLOAD )
                                    ! kg--> l

         ! Phosphorus:
      CALL SETR ( 0.0, PLOAD, NBI )
      CALL QCALC( NQDIST, RCQNDX, PO4DV, QPO4, 1.e6*RNF,'PO4' ,PLOAD)
      CALL QCALC( NQDIST, RCQNDX, PDETDV,QPDET,1.e6*RNF,'PDET',PLOAD)
                                    ! kg --> mg

         ! Nitrogen load from land distributed on NITRMP(=NO3)
         ! and NLOAD (ammonium and particulate)
      CALL SETR ( 0.0, NITRMP, NBI )
      CALL SETR ( 0.0, NLOAD, NBI )
      CALL QCALC( NQDIST, RCQNDX, NO3DV, QNO3, 1.e6*RNF,'NO3' ,NITRMP)
      CALL QCALC( NQDIST, RCQNDX, NH4DV, QNH4, 1.e6*RNF,'NH4' ,NLOAD )
      CALL QCALC( NQDIST, RCQNDX, NDETDV,QNDET,1.e6*RNF,'NDET',NLOAD)
         ! NITRMB will be changed in SUMIMP to include all N components

         ! Silicate:
      CALL QCALC ( NQDIST, RCQNDX, SiO2DV, QSiO2, 1.e6*RNF, 'SIO2' , SiO2MP)
                                                ! kg--> mg
      

         ! Organic carbon to detrital fraction CDET and degradability RDET
      CALL QCALC ( NQDIST, RCQNDX, RDETDV, QCDET, DGRATE(3)*1.e6*RNF, 'RDET', CLOAD )
                       !  NOTE: Dummy update of CLOAD here, only because function QCALC
                       !        requires a parameter. CLOAD is really set below:

      CALL SETR ( 0.0, CLOAD, NBI )
      CALL QCALC ( NQDIST, RCQNDX, CDETDV, QCDET, 1.e6*RNF, 'CDET' , CLOAD)
                                    ! CLOAD as mg/day

         ! Organic carbon to DOC:
      CALL QCALC ( NQDIST, RCQNDX, DOCDV, QDOC, 1.e6*RNF, 'DOC' , CLOAD)

         !###### ! Phosphorus detritus included twice? (See above)
      CALL QCALC ( NQDIST, RCQNDX, PDETDV, QPDET, 1.e6*RNF, 'PDET' , PDETMP )

         !  Heat content:
      CALL QCALC ( NQDIST, RCQNDX, TEMPDV, QTEMP, 1.0*RNF, 'TEMP' , HEATMP  )



        ! %%%%%%%%%%%%%% recent addition - inorganic particle influx %%%%%%%%%%%%%%
!       CALL DVZERO(1, NLI, SPPDV, SPPMP)
!       CALL SETR( 0.0, SPPSEDDV, dimMLI) 

      do indx_MS = 1, dimMS
         if ( (QFSPP(indx_MS).gt.0.0) .and. &
             (T.ge.TIME_SPP(1)) .and. (T.lt. TIME_SPP(2)) ) then
            QSPP(indx_MS) = QFSPP(indx_MS)
         else
            QSPP(indx_MS) = 0.0
         endif
      end do  
         ! Update derivatives with effect of influx:
      CALL QCALC ( NQDIST, RCQNDX, SPPDV, QSPP, 24.*3600.*RNF, 'SPP'  , SPPMP )
            ! (SPPDV as g/m3/day, SPPMP as g/day)
                ! g/s -> g/day, i.e. conc. in g/m3, or mg/l

            
         ! Inorganic particles -  Added March 2001 - BBJ
         ! DRAFT CODE -not yet active
!     CALL SETR( 0.0, SPPEX, dimMLI) 
                ! Update SPPDV and set SPPSEDDV, SPPFLUX:
                ! with result of sinking rates:
!     CALL SPSINK( SPP_SINK_VELOCITY)
                ! Update derivatives and basin import
                ! with result of transports:
!     CALL MTRANS ( SPPDV, SPPMP, 1, SPP, SPPEX, VTRNEG, 7, 'SPP')
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

         !  Mussel submodel; Update derivates:
      CALL MUSLDV( MSTEST .and. T.ge.TTRIG, dimMFYTG,  &
                     dimMUSLAGES, dimMUSLAYERS, NBI,  &
                     MUSLNR, MUSLMT, MUSLMA, MUSLGT )
      CALL HELLO ('returned from MUSLDV')

      
         !  ------------- solar and diffuse light radiation -------------
 
      CALL LGTRAD( T, CLOUDS, MAXSTP, DAYDIV, SINHS, SUNRAD, DIFRAD, TSTEP )

         ! MAXTTR = max timestep from TRANSP, modified into TSTEP
         ! to give sufficient resolution of day-light period.

         ! LGTRAD returns mean light radiation for timestep TSTEP

         
         ! -------------- Integration interval:


         !  TSTEP specified in TRANSP to avoid numeric instabilities due
         !  to transports in general, below limit MAXINT set by other considerations
 

      CALL TRNADJ ( NBI, INDXI, NLI, VLAYER, TSTEP, TRCALC, VTRNEG )

         ! Adjusts transports to time step TSTEP, 
         ! after TSTEP has been set by biological subroutines. 
         ! Necessary because diffusive effect of advection
         ! is compensated as negative diffusion
         ! NOTE: Argument TRCALC served to sort TRNADJ after TRANSP in old ACSL Code
         !       Reconsider - is it used?
         
         ! The logical return argument VTRNEG will be set to .true.
         ! to signal occurrence of negative vertical transports occur,
         ! this can happen if advection dominates over diffusion 
         ! between some model layers ( see subr. TRANV1 in TRANSP_3.FOR)
         ! May suggest using a finer vertical resolution ??
      
         ! Final time step TSTEP has been set by LGTRAD <= MAXTTR
         !  VTRNEG was also used to sort MTRANS after TRNADJ in old ACSL Code


         ! Mass transports due to Water transports:
         ! NOTE: will initiate state derivatives:
         ! (calls to DVZERO not needed? Check!)
         

         ! --- Calculate water transports between layers and across boundaries
         !     Volume transports were hidden from the ACSL translator in the old model,
         !     stored in common blocks accessed by FORTRAN subroutines:
         !           TRANSP  - calculates the water fluxes and returns
         !                time derivative of total volume in each basin.
         !           MASS_TRANSPORT - uses the calculated transports,
         !                     and calculates time derivatives of the
         !                     amount of a substance  in each layer due
         !                     to water transports, via subroutine MTRANS


   
      CALL MTRANS ( C1DV , C1MP  , 1, C1 , C1EX , 1, 'C1  ' )
      CALL MTRANS ( SALDV, SALTMP, 1, SAL, SALEX, 2, 'SAL '  )
      CALL MTRANS ( TEMPDV, HEATMP, 1, TEMP, TEMPEX, 2, 'TEMP')
      CALL MTRANS ( OXYGDV, OXYGMP, 1, OXYG, OXYGEX, 3, 'OXYG')
                                 ! -(Oxygen unit ml/l)
      CALL MTRANS ( PO4DV, PO4MP, 1, PO4, PO4EX, 4, 'PO4')
      CALL MTRANS ( NO3DV, NO3MP,1, NO3, NO3EX, 4, 'NO3')
      CALL MTRANS ( NH4DV, NH4MP ,1, NH4, NH4EX, 4, 'NH4'  )
      CALL MTRANS ( SiO2DV, SiO2MP,1, SiO2, SiO2EX, 4, 'SiO2')

         !  -------- passive transport of pelagic biological entities ---------
         !  import terms for water transport effect only, is added to
         !  terms for Oxygen, NH4 and PO4 and SiO2 for integration
         !  and mass conservation control

         !  --- Phytoplankton, described by mass substances C, N and P,
         !      and chlorofyll, treated as property, not as mass substance.
         !    Chlorofyll has unit (C-conc)/day/(W/m2) and really is
         !    a measure of the quantum yield

      CALL MTRANS ( CFYTDV, CFYTMP,dimMFYTG, CFYT, CFYTEX, 5, 'CFYT'  )
      CALL MTRANS ( NFYTDV, NFYTMP,dimMFYTG, NFYT, NFYTEX, 5, 'NFYT' )
      CALL MTRANS ( PFYTDV, PFYTMP,dimMFYTG, PFYT, PFYTEX, 5, 'PFYT' )
      CALL MTRANS ( SFYTDV, SFYTMP,1, SFYT, SFYTEX, 5, 'PFYT' )
      CALL MTRANS ( CHLDV, CHLMP,dimMFYTG, CHL, CHLEX, 5, 'CHL' )

      CALL MTRANS ( ODMDev, ODMImp,1, ODM, ODMExt, 5, 'ODM' )
      CALL MTRANS ( DOCDV, DOCMP,1, DOC, DOCEX, 5, 'DOC' )

         !  ##############################################################
         !    simplified treatment of external concentrations,
         !    should be incorporated in HYDREX later
         !  ##############################################################

      CALL MTRANS ( BACTDV, BACTMP, 1 , BACT, BACTEX, 5, 'BACT' )
      CALL MTRANS ( CZOODV, CZOOMP, 1, CZOO, CZOOEX, 6, 'CZOO' )
      CALL MTRANS ( CDETDV, CDETMP, 1, CDET, CDETEX, 5, 'CDET' )
      CALL MTRANS ( NDETDV, NDETMP, 1, NDET, NDETEX, 5, 'NDET' )
      CALL MTRANS ( PDETDV, PDETMP, 1, PDET, PDETEX, 5, 'PDET' )
      CALL MTRANS ( SDETDV, SDETMP, 1, SDET, SDETEX, 5, 'SDET' )
      CALL MTRANS ( RDETDV, RDETMP, 1, RDET, RDETEX, 5, 'RDET' )

         ! ----- Calculate Radiation (RAD) and heat balance (QABS)
         !       as a function of depth:

      CALL RADABS ( NBI, INDXI, NLI, WNDSPD, ND, DEPTH,   &
                QOUT,    &
                SUNRAD,    &
                DIFRAD,    &
                SINHS,    &
                IRFRAC,    &
                ICEFAC,  &
                RADFAC,    &
                PARTC,    &
                ATTNCF,  &
                AREA,    &
                VLAYER,    &
                TEMP,    &
                SAL,    &
                CEFAC,    &
                TEMPDV,    &
                HEATMP, &
                RAD,    &
                QABS )

               ! If CEFAC >0:
               !   Also updates temperature derivative with effect of solar
               !   radiation at the surface and distributed in the water
               !   column, and with heat balance terms at the surface.
               !   TEMPDV and HEATMP initiated by MTRANS are included as
               !   input arguments to get correct sorting by ACSL



         !  Primary production, planton grazing and biomass sedimentation:
         !  updates derivatives, previously set:
         !  CSEDDV, NSEDDV, PSEDDV, SSEDDV by degradation, others by MTRANS.
         ! also provides monitoring of primary production
       
      CALL PRPROD

         ! Print progress monitoring message each NPRINT time-step
         
      CALL TIMPRT("At end of Calc_Derivatives")

         ! Print selected output to text files:
         
!      CALL out_DEEPW("Deep water output")
!      CALL out_SURW("Surface water output")      
!      CALL out_OXYG("Oxygen output")
!      CALL out_SAL("Salt output")
!      CALL out_TEMP("Temp output")
!      CALL out_TOTC("C output")
!      CALL out_DOC("C output")
!      CALL out_TOTN("N output")
!      CALL out_NO3("N output")
!      CALL out_NH4("N output")
!      CALL out_TOTP("P output")
!      CALL out_PO4("P output")
!      CALL out_CZOO("Bio output")
!      CALL out_CFYT1("Bio output")
!      CALL out_CFYT2("Bio output")
!      CALL out_CZOO("Bio output")
!      CALL out_CHL1("KlfA output")
!      CALL out_Si("Si output")
!      CALL out_MUSSEL("Mussel output")

         !  Combine imports of different forms into total import rates:

      CALL SUMIMP

         !  NITRMP and PHOSMP are defined, and OXYGMP is updated
         !  with effect of water transport/runoff and permanent sedimentation.
         !  Must be done after PRPROD


         !  ------- Buoyancy energy balance:

         ! BFX is used to control surface mixing:
         ! Set BFX to give surface buoyancy increase
         ! due to land runoff, exchange with atmosphere
         ! and horisontal transports from other basins.
         ! Effect of transport from other basins have been
         ! set in BSFLUX by call to subroutine TRANSP,
         ! which has already been called due to ACSL sorting

         ! ********* old call:

      CALL SURFBF ( MXTEST, dimMBI, NBI, QWSURF, BSFLUX, NS,  &
                 QWATER, BASINQ, RNFNDX, &
                 QOUT, QABS(1), INDXI, NLI, AREA, SAL, TEMP, &
                 BFX)

         ! ********* new call (not yet used):

!       CALL SURFBF ( MXTEST, dimMBI, NBI, QWSURF, BSFLUX, NS,  &
!                     QWATER, QSPP, QTEMP, BASINQ, RNFNDX, &
!                     QOUT, QABS(1), INDXI, NLI, AREA, SAL, TEMP, &
!                     dDens_dSPP, DENSI,
!                     BFX,  SPP_ToSurf, QTemp_ToSurf )


      end subroutine Calc_Derivatives
      
end module
