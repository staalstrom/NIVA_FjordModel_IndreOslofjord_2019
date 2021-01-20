Module ModelParam_Plankton

  use ModelParamToolbox
  use ModelDimensions

  implicit none
 
  ! Parameter declarations for plankton processes


         integer LFYT   ! number of fytoplankton groups used
                        ! effective within dimension limit dimMFYTG
   
   
                ! ******** fytoplankton ***********
   
                        ! Process parameters for fytoplankton groups:
         real*4 GMX20 (dimMFYTG)  ! Growth rates
   
         real*4 FTRESP(dimMFYTG)  ! Temperature response coefficient
                                 ! for phytoplankton growth

   
         real*4 FDEATH(dimMFYTG)  ! Asymptotic maximum for inherent death rate

   
         real*4 FDCSAT(dimMFYTG)  ! Half saturation concentrations
                                 ! for death rate
   
         real*4 FDNUTR(dimMFYTG)  ! Fraction of nutrients lost in death process,
                                 ! residual due to within-biomass heterotrophy
   
         real*4 RESP20    ! Dark respiration rate at 20 degC,
                        ! common to all groups
   
         real*4 EXCRF(2)  ! Additional organic carbon excretion relative to
                        ! light- and temperature-limited max. gross growth,
                        ! increasing with nutrient limitation:
                        !  = EXCRF(1) at nutrient sufficient growth,
                        !  + EXCRF(2)*Nutrient limitation factor
   
                   ! Grazing availability for bacteria and phytoplankton
         real*4 GRZBAC    ! for bacteria
         real*4 GRZFYT(2) ! for fytoplankton
   
   
                   ! Diatom sedimentation characteristics:
                   ! depends on nutrient and light growth limitation:
   
         real*4 DSRATE(2) ! Minimum and maximum sedimentation rate (1/day)
   
         real*4 DSNINV    ! inverse of nutrient limitation at full sinking rate
   
         real*4 DSNEXP    ! Exponent in transition (1-NUTLIM*SDNUTL)**DSNEXP
   
         real*4 DSCLIM    ! Threshold diatom density as carbon,
                        ! no increased sedmentation at lower densities,
                        ! increasing by Monod kinetics
   
         real*4 RESUSP    ! Resuspension for sedimenting material
   
   
                   
         real*4 SEDVEL(2) ! Sinking velocity of detritus
                        ! assumed to increase with depth:
                        !   (1): minimum velocity (m/day)
                        !   (2): increase (m/day/m)
   
                   !  Nutrient uptake characteristics:
            ! ---max uptake-----  ---half saturation constants--
         real*4 VMNH4(2), VMNO3(2), KSNH4(2), KSNO3(2), NH4EXP(2)  ! for N
         real*4 VMPO4(2),           KSPO4(2)                       ! for P
         real*4 VMSiO2  ,           KSSiO2                         ! for Si,
                                         ! only first phytoplankton group
   
         real*4 PLUXURY         ! P luxury uptake (factor on optimal P:C ratio)
   
         real*4 NFIXRR          ! Nitrogen fixation ability
   
         real*4 F2SINK, F2RIZE  ! Flagellate ability to move vertically,
                              ! as max. velocities (m/day):
   
   
   
                ! ********* Zooplankton **********:
   
   
         real*4 ZFCOMP          ! ability to compensate lack of nutrients in food
                              ! by increased filtering and/or selective ingestion :
   
                              ! Grazing/growth characteristics:
         real*4 ZFMX20          ! (/day)  Max. relative ration at T=20oC
   
         real*4 ZTRESP          ! Temperature response coefficient
                              ! for zooplankton activity
   
         real*4 ZOOEFF(3)       ! max. assimilation = growth efficiency
                              ! for carbon, nitrogen and phosphorus
   
         real*4 ZCFMIN          ! (microgC/l) Food conc. where grazing stops
   
         real*4 ZCFSAT          ! (microgC/l) Food half saturation conc.
   
         real*4 ZGCYCL          ! Fraction of uningested material recycled
                              ! to water column; the rest will sediment
                              ! as particles
   
                        ! Death/recycling by higher predators:
         real*4 ZOODR (2)       ! (/day) min and max. rate
                              ! Max. rate due to predators
                              ! - concentration dependent
                              ! controlled by critical zooplankton conc.:
   
         real*4 ZCCRIT(2)       ! (1): lower limit (no predator activity)
                              ! (2): 50% saturation level for ZOODR(2)
   
         real*4 ZDCYCL          ! Fraction of dead zooplankton
                              ! recycled without sedimentation
   
                              ! Oxygen tolerance/preference:
         real*4 ZOXMIN          ! (ml/l) oxygen limit for zooplankton
         real*4 ZOXOPT          ! (ml(l) oxygen half saturation value
                              !        for zooplankton activity
   
   
         real*4 ZRESP           ! (1/day) Relative respiration at T=20oC
   
                        ! Vertical migration:
         real*4 ZMIGRV          ! (m/day) Maximum migration velocity
         real*4 ZMIGRH          ! (m) controlling vertical dimension


                ! ******** stochiometric ratios *********

                   !  Minimal and optimal ratios in fytoplankton
                   !     (unit weight:weight)
                   !  possibly different for diatoms 
                   !     versus other phytoplankton.

         real*4 NCMIN(dimMFYTG), NCOPT(dimMFYTG)  ! nitrogen:carbon
         real*4 PCMIN(dimMFYTG), PCOPT(dimMFYTG)  ! phosphorus:carbon


                   ! Silicon requirements only in first group of phytoplankton,
                   ! can be deactivated as nutrient by setting ratios to zero
         real*4 SCMIN, SCOPT

                 !  Fixed ratios for bacteria and zooplankton:
         real*4 NCZOO , PCZOO
         real*4 NCBACT, PCBACT



  contains

    subroutine ParamGroup_Plankton

     !  Initiate parameters with default values, and register them
     !  for actions (user editing, file input, file storage)


       integer i     ! used to assign arrays

           ! Start parameter Group:

       call paramGroup("Plankton parameters","PLANKTON")


           ! Initiate and declare parameters:

	    call ParamDeclInteger( LFYT, min(2, dimMFYTG), "LFYT", "", &
                 "Number of fytoplankton groups used. "//  &
                 "Effective within dimension limit dimMFYTG")

           ! Phytoplanktn growth rate:

       call ParamDeclRealArray(GMX20,(/2.2, (0.7,i=2,dimMFYTG)/), "GMX20", "(1/day)",       &
                 "Maximum obtainable specific growth rates (1/d) at 20 degC "//  &
                 " with optimal nutrient ratios (< asymptotic rates)." // &
                 " One rate for each phytoplankton group.")

	    call ParamDeclRealArray( FTRESP,(/0.063, (0.063,i=2,dimMFYTG)/), &
                     "FTRESP","1/degC)", &
                 "Temperature response coefficients for phytoplankton growth"//  &
                 " as coefficient in relation exp((Temp-20)**FTRESP)."//&
                 " One rate for each phytoplankton group.")
	    call ParamDeclRealArray(FDEATH, (/0.3, (0.5,i=2,dimMFYTG)/), &
                     "FDEATH","1/day", &
	              "Asymptotic maximum for inherent death rates" // &
                 " of phytoplankton groups at 20 deg C" // &
                 " at high population densities." // &
                 " One rate for each phytoplankton group.")

       call ParamDeclRealArray(FDCSAT, (/1000.0, (1000.0 ,i=2,dimMFYTG)/), &
                     "FDCSAT"," micro-g C/l", &
                 "Half saturation fytoplankton concentrations as carbon" // &
                 " for death rate saturation function. Assumed to" // &
                 " describe stress, increase of heterotrophic components" // &
                 " or switch of metabolism for mixotrophic species." // &
                 " One rate for each phytoplankton group.")

      call ParamDeclRealArray(FDNUTR,(/0.2, (0.2,i=2,dimMFYTG)/), "FDNUTR","", &
                 "Fraction of nutrients lost in death process," // &
                 " the rest recovered, connected to within-biomass heterotrophy." // &
                 " One rate for each phytoplankton group.")
                 
             ! NOTE - should have cross-connections between phytoplankton groups 
             ! to cover mixotrophy

      call ParamDeclReal(RESP20, 0.04, "RESP20", "1/day", &
                 "Dark respiration rate at 20 degC " )

      call ParamDeclRealArray(EXCRF, (/0.2, 0.3 /), "EXCRF"," ", &
                 " Additional organic carbon excretion relative to " //&
                 "light- and temperature-limited max. gross growth, " //&
                 "increasing with nutrient limitation: "//&
                 "= EXCRF(1) at nutrient sufficient growth, "//&
                 "+ EXCRF(2)*Nutrient limitation factor")


          ! Grazing availability of bacteria and phytoplankton

      call ParamDeclReal( GRZBAC, 1.0,"GRZBAC", "", &
                 "Grazing availability of bacteria")

      call ParamDeclRealArray(GRZFYT, (/0.1, (1.0,i=2,dimMFYTG)/), &
                      "GRZFYT","", &
                 "Grazing availability of phytoplankton." // &
                 " One rate for each phytoplankton group.")

               ! Sinking of plankton and detritus:

      call ParamDeclRealArray(DSRATE, (/0.04, 0.6/),"DSRATE","1/day",&
                 "Diatom sedimentation characteristics:" //&
                 " Minimum and maximum sedimentation rate"//&
                 " (increases with nutrient and light growth limitation)")

      call ParamDeclReal(DSNINV, 2.0,"DSNINV","day",&
                 "Inverse of nutrient limitation at full sinking rate")

      call ParamDeclReal(DSNEXP, 2.5,"DSNEXP","",&
                 "Exponent in sinking rate factor: (1-NUTLIM*DSNINV)**DSNEXP")

      call ParamDeclReal(DSCLIM, 10.0, "DSCLIM","micro-g C/l",&
                 "Threshold diatom density for sinking rate," //&
                 " increase by Monod kinetics for higher concentrations")

      call ParamDeclReal(RESUSP, 0.1,"RESUSP","",&
                 "Fraction of sedimenting flux which is resuspended"//&
                 " (reduces area-proportionate sedimentation within each"//&
                 " depth interval, and causes sedimentation to be shifted"//&
                 " towards greater depth")

      call ParamDeclRealArray(SEDVEL, (/4.0, 0.5/), "SEDVEL","(m/day)",&
                 "Sinking velocity of detritus" // &
                  " (1): velocity at surface, (2): increase pr. m depth")


              !  Nutrient uptake characteristics:

                   ! ... Maximum relative uptake:

        call ParamDeclRealArray(VMNH4, (/0.9,  (0.6,i=2,dimMFYTG)/), &
                      "VMNH4","gN/gC/day",&
                 "Maximum relative uptake of ammonium in fytoplankton."// &
                 " One rate for each phytoplankton group.")


        call ParamDeclRealArray(VMNO3, (/0.6,  (0.4,i=2,dimMFYTG)/), &
                     "VMNO3","gN/gC/day",&
                 "Maximum relative uptake of nitrate in fytoplankton." // &
                 " One rate for each phytoplankton group.")

        call ParamDeclRealArray(VMPO4, (/0.20, (0.13,i=2,dimMFYTG)/), &
                      "VMPO4","gP/gC/day",&
                 "Maximum relative uptake of phosphate in fytoplankton." // &
                 " One rate for each phytoplankton group.")

        call ParamDeclReal(VMSiO2, 0.5, "VMSiO2","gSi/gC/day",&
                 "Maximum relative uptake of silicate in diatoms.")


                   ! ... Half saturation concentrations:

        call ParamDeclRealArray( KSNO3, (/7.0, (7.0,i=2,dimMFYTG)/), &
                     "KSNO3","micro-g N/l",&
                 "Half saturation concentration in water for nitrate uptake")

        call ParamDeclRealArray(KSNH4, (/7.0, (7.0,i=2,dimMFYTG)/), &
                     "KSNH4","micro-g N/l", &
                 "Half saturation concentration in water for ammonium uptake")

        call ParamDeclRealArray(KSPO4,(/3.0, (3.0,i=2,dimMFYTG)/), &
                      "KSPO4","micro-g P/l",&
                 "Half saturation concentration in water for phosphate uptake")

        call ParamDeclReal(KSSiO2,90.0, "KSSiO2","micro-g Si/l",&
                 "Half saturation concentration in water for silicate uptake")

        call ParamDeclRealArray(NH4EXP, (/3.0, (3.0,i=2,dimMFYTG)/), &
                      "NH4EXP","",&
                 "Exponential power in NO3 inhibition")

        call ParamDeclReal(PLUXURY, 2.0, "PLUXURY","",&
                 "P luxury uptake (factor on optimal P:C ratio)")

        call ParamDeclReal(NFIXRR, 1.0, "NFIXRR","(1/year)",&
                 "Nitrogen fixation ability (to reduce algal N:P deficit)"// &
                 " relative rate of phytoplankton group 2 (flagellates)")

        call ParamDeclReal(F2SINK, 5.0,"F2SINK","(m/day)",&
                 "Flagellate max. downward velocity")

        call ParamDeclReal(F2RIZE, 10.0, "F2RIZE", "(m/day)",&
                 "Flagellate max. upward migration velocity")


             ! ---------- Zooplankton: -----------------

      call ParamDeclReal( ZFCOMP, 0.5, "ZFCOMP","",&
                 "Zooplankton ability to compensate lack of nutrients in food"//&
                 " by increased filtering and/or selective ingestion"// &
                 " [0...1] = [no compensation...full compensationn")

      call ParamDeclReal(ZFMX20, 1.5, "ZFMX20","(1/day)",&
                 "Max. relative ration for zooplankton at T=20oC")

      call ParamDeclReal(ZTRESP, 0.05, "ZTRESP","",&
                 "Temperature response coefficient for zooplankton activity")

      call ParamDeclRealArray(ZOOEFF, (/0.6,0.8,0.87/), &
                       "ZOOEFF","(0...1)",&
                 "Max. fraction of grazed biomass assimilated (=growth efficiency)"// &
                 " for carbon, nitrogen and phosphorus")

      call ParamDeclReal(ZCFMIN, 10.0, "ZCFMIN","(micro-g C/l)", &
                 "Food conc. where grazing stops " )

      call ParamDeclReal(ZCFSAT, 2000.0, "ZCFSAT","(micro-gC/l)", &
                 " Food half saturation conc.")


      call ParamDeclReal(ZGCYCL, 0.3, "ZGCYCL","[0...1]", &
                 "Fraction of uningested material that is recycled."// &
                 " The rest will sediment as particulate matter")


               ! Death/recycling of zooplankton by predators
               ! implemented only as prcoess rates, not as predator biomass:

      call ParamDeclRealArray(ZOODR, (/0.03, 1.0/), "ZOODR","(1/day)", &
                 "Relative death rates for zooplankton" //&
                 " 1: max. rate due to 'auto-predation' at 20 deg.C" //&
                 " and good oxygen conditions")

      call ParamDeclRealArray(ZCCRIT, (/10.0, 1000.0/), &
                       "ZCCRIT","(micro-g C/l)", &
                 "Critical zooplankton concentrations  "// &
                 " controlling predator-related Zooplankton death (ZOODR(2):" // &
                 " 1: lower limit for predator activity within zooplankton" //&
                 " 2: 50% saturation level" )

      call ParamDeclReal( ZDCYCL , 0.3, "ZDCYCL ", "", &
                       "Fraction of dead zooplankton " // &
                       " recycled without sedimentation" )

            ! Oxygen tolerance/preference:

      call ParamDeclReal( ZOXMIN, 1.0, "ZOXMIN", " (ml/l)", &
                       "Oxygen limit for zooplankton" )

      call ParamDeclReal(ZOXOPT, 2.0, "ZOXOPT",  " (ml/l)", &
                       "Oxygen half saturation value" // &
                       " for zooplankton activity" )

      call ParamDeclReal(ZRESP, 0.05, "ZRESP", "(1/day)", &
                       " Relative respiration at T=20oC" )

             ! Vertical zooplankton migration:
      call ParamDeclReal(ZMIGRV, 10.0, "ZMIGRV", "(m/day)", &
                      "Maximum migration velocity" )

      call ParamDeclReal(ZMIGRH , 5.0, "ZMIGRH ", "(m)", &
                      "Controlling vertical dimension")

             !  stochiometric ratios

					   !  Minimal and optimal ratios in fytoplankton (unit weight:weight)
					   !  of nitrogen:carbon and phosphorus:carbon,
					   !  possibly different for diatoms versus other phytoplankton:
					   
      call ParamDeclRealArray(NCMIN, (/0.06  , (0.06,i=2,dimMFYTG)/) , "NCMIN ", "w:w", &
                      "Minimum Nitrogen:Carbon ratio in phytoplankton (weight:weight) ") 

      call ParamDeclRealArray(NCOPT, (/0.18,  (0.18,i=2,dimMFYTG)/), "NCOPT", "w:w",  &
                      "Optimal Nitrogen:Carbon ratio in phytoplankton (weight:weight) ") 

      call ParamDeclRealArray(PCMIN, (/0.0027, (0.0027,i=2,dimMFYTG)/), "PCMIN", "w:w",  &
                      "Minimum Phosphorus:Carbon ratio in phytoplankton (weight:weight) ") 

      call ParamDeclRealArray(PCOPT, (/0.027, (0.027,i=2,dimMFYTG)/), "PCOPT ", "w:w", &
                      "Optimal Nitrogen:Carbon ratio in phytoplankton (weight:weight) ") 


				    ! Silicon requirements only in first group of phytoplankton,
				    ! can be deactivated as nutrient by setting ratios to zero

      call ParamDeclReal(SCMIN, 0.09, "SCMIN", "(weight:weight)", &
                     "Minimum silisium:Carbon ratio in group 1 " // &
                     " of phytoplankton") 
      
      call ParamDeclReal(SCOPT , 0.16, "SCOPT ", "(weight:weight)", &
                     "Optimum silisium:Carbon ratio in group 1 " // &
                     " of phytoplankton") 


              !  Based on atomic weights N=14, C=12, P=31, Si=28

      call ParamDeclReal(NCZOO, 0.18, "NCZOO", "(weight:weight) ", &
                     "Fixed Nitrogen:Carbon ratio for zooplankton")

      call ParamDeclReal(PCZOO, 0.027,  "PCZOO", "(weight:weight)" , &
                     "Fixed Phosphorus:Carbon ratio for zooplankton")
      
      call ParamDeclReal(NCBACT, 0.18, "NCBACT", "(weight:weight)", &
                     "Fixed Nitrogen:Carbon ratio for bacteria")

      call ParamDeclReal(PCBACT, 0.027, "PCBACT", "(weight:weight)", &
                     "Fixed Phosphorus:Carbon ratio for bacteria")


      end subroutine

end Module  ModelParam_Plankton
