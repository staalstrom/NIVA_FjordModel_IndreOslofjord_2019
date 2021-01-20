Module ModelParam_InitState

  use ModelParamToolbox
  use ModelDimensions

  implicit none

      ! Parameter declarations as elements in structure p_DC,
      ! located in common block Common_Param_Decomposition
      ! declared in ModelParam_Decomposition.inc

  
           ! Initial contents of organic matter in sediments:
      
      real*4 CSEDIN(dimMBI), NSEDIN(dimMBI), PSEDIN(dimMBI)
      real*4 RSEDIN(dimMBI), SSEDIN(dimMBI), ASEDIN(dimMBI)

      real*4 XSEDZ(dimMBI)
		      ! Reduction factors for sedimented matter 
            ! assumed to be presentat start of simulation,
            ! but not included in inital value ?SED
		      !   See subroutine DGRADE for details.
		      ! XSEDZ = 0 will turn this feature off

      real*4 CDRSED  (dimMBI)
            ! Background values of carbon degradation
            ! in sediments (mg/m2/day)

      real*4 CDRDEPTH(dimMBI)
            ! Upper limit to bottom area where CDRSED applies.

            ! ----------- Initial concentrations in water.
      
      real*4 PO4IN, NO3IN, NH4IN, SIO2IN
            ! Initial concentrations of nutrients

      real*4 CFYTIN(2), NFYTIN(2), PFYTIN(2), SFYTIN
            ! Initial concentration and composition of
            ! phytoplankton: (1): diatoms, (2): flagellates

      real*4 DOCIN, BACTIN, CZOOIN
            ! Initial concentration of DOC, bacteria and zooplankton


		      ! ------------ Continuity and residence time;
		      ! More flexible setup in eutro.csl of POPAppl version
		      ! Should be implemented.

      real*4 C1XTRN
				! controls the use of state variable C1:
				! 1.0 : Unity concentration  - acts as continuity check
				! 0.0 : will have C1 to act as residence time in system

      real*4 C1ZERO(dimMBI)
            !  Only active if C1XTRN = 0.0; in that case forces C1 =0
	         !  for basins where C1ZERO is .true.


!  ====================== mussel population ============================

      real*4 MCOVER(dimMBI)
        !  Dimensioning fraction of bottom covered,
        !  excess leads to increased mortality:

      real*4 MUSLDP
        ! max. depth of mussel mats:

      real*4 CMUSIN
        ! Total initial biomass of mussels, as softbody carbon

 
  contains

! -------------------------------------------------------------

      subroutine ParamGroup_InitState

			! Initiate parameters with default values
         ! and register them for actions
         ! (user editing, file input, file storage)

      call paramGroup("Model State Initiation","INITIATION")


      call ParamExpl("Initial contents of organic matter in sediments," // &
                     " one value per inner basin:")
                     
         ! Arrays with dimMBI values:            

      call ParamDeclRealArray(CSEDIN,(/0.0/),"CSEDIN","mgC/m2", &
                "Organic carbon")

      call ParamDeclRealArray(NSEDIN,(/0.0/),"NSEDIN","mgN/m2", &
                "Nitrogen")

      call ParamDeclRealArray(PSEDIN,(/0.0/),"PSEDIN","mgP/m2", &
                "Phosphorus")

      call ParamDeclRealArray(RSEDIN,(/0.0/),"RSEDIN","mgC/m2/day", &
                "Remineralisation capacity:" // &
                " Organic carbon*decomposition rate" // &
                " (state variable R in model description)")

      call ParamDeclRealArray(SSEDIN,(/0.0/),"SSEDIN","mgSi/m2", &
                "Silicon")

      call ParamDeclRealArray(ASEDIN,(/0.0/),"ASEDIN","liter O2/m2", &
                "Oxygen debt in sediments (mainly occurs as sulfides)")

      call ParamDeclRealArray(XSEDZ,(/1.0/),"XSEDZ","", &
                "Controls degree of adjustment for sedimented matter" // &
                " assumed to be present at start of simulation," // &
                " but not included in inital values ?SEDIN" // &
                " (Helps model to get realistic long-term conditions" // &
                " sooner after start of model run.)" // &
                " May be varied continously between"// &
                " 1.0: Full adjustment, 0.0: no adjustment")

      call ParamDeclRealArray(CDRSED,(/0.0/),"CDRSED","mg/m2/day", &
                "Background values for organic carbon decomposition" // &
                " in deep sediments in addition to rates based on " // &
                " accumulated organic matter from sinking organic matter" // &
                " to represent old organic load of sediments")

      call ParamDeclRealArray(CDRDEPTH, (/0.0/), "CDRDEPTH", "m", &
                "Upward depth limit for bottom area" // &
                " where CDRSED applies")
      

      call ParamExpl("Initial concentrations of nutrients;" // &
                " one value applied to all basin layers:")

      call ParamDeclReal(PO4IN, 20.0, "PO4IN","µg P/liter", &
                 "Initial concentration of orthophospate")

      call ParamDeclReal(NO3IN, 200., "NO3IN","µg N/liter", &
                 "Initial concentration of nitrate+nitrite")

      call ParamDeclReal(NH4IN, 100., "NH4IN","µg N/liter", &
                 "Initial concentration of ammonium")

      call ParamDeclReal(SIO2IN, 750.0, "SIO2IN","µg Si/liter", &
                 "Initial concentration of silicate")
      
      
      call ParamExpl("Initial amount and composition of phytoplankton;" // &
                " as concentration of unfiltered water, with" // &
                " separate value for diatoms and flagellates." // &
                " The same concentration is applied to all basin layers:")

      call ParamDeclRealArray(CFYTIN, (/42.0, 42.0/), &
                                    "CFYTIN","µg C/liter", &
                "Carbon in (1): diatoms, (2): flagellates")

      call ParamDeclRealArray(NFYTIN, (/7.2, 7.2/), &
                                    "NFYTIN","µg N/liter", &
                "Nitrogen in (1): diatoms, (2): flagellates")

      call ParamDeclRealArray(PFYTIN, (/1.0, 1.0/), &
                                    "PFYTIN","µg P/liter", &
                "Phosphorus in (1): diatoms, (2): flagellates")

      call ParamDeclReal(SFYTIN, 13.5, "SFYTIN","µg Si/liter", &
                "Silicon  in diatoms")

      call ParamExpl("Initial amount of other components in water," // &
                " The same concentration is applied to all basin layers:")


      call ParamDeclReal(DOCIN, 2000.0,"DOCIN","µgC/liter", &
                " Initial conc. of dissolved organic carbon")

      call ParamDeclReal(BACTIN, 10.0,"BACTIN","µgC/liter", &
                "Initial conc. of bacteria," // &
                " must be >0.0 to activate bacteria compartment")

      call ParamDeclReal(CZOOIN, 1.0,"CZOOIN","µgC/liter", &
                " Initial conc. of dissolved organic carbon" // &
                " must be >0.0 to activate zooplankton compartment")



      call ParamExpl("continuity and residence time:"//&
                " (More flexible setup in eutro.csl of POPAppl version" //&
                " Should be implemented.)")

      call ParamDeclReal(C1XTRN, 1.0, "C1XTRN","", &
			       "Controls the use of state variable C1:" // &
			       " = 0.0 : C1 will be the average time the water" // &
                " in each layer has stayed within specified basins" // &
                " in the model area (ref. C1ZERO)." // &
                " not = 0.0: Constant concentration in all water" // &
                " - used for continuity check)")

      call ParamDeclRealArray(C1ZERO, (/0.0/),"C1ZERO","",&
                "Only active if C1XTRN = 0.0; in that case keeps" // &
                " C1 = 0 in basins for which C1ZERO =0," // &
                " so C1 will be residence time within" // &
                " other parts of the model")

      call ParamDeclRealArray(MCOVER, (/0.1/),"MCOVER","", &
                "Critical fraction of bottom area" // &
                " that can be covered by mussels," // &
                "  exceeding this limit causes increased mortality")

      call ParamDeclReal(MUSLDP, 6.0, "MUSLDP","m", &
                "Maximum depth of mussel settling in benthic zone")

      call ParamDeclReal(CMUSIN, 1200.0E+09, "CMUSIN", "mg carbon", &
                "Total initial biomass of mussels, as softbody carbon")

   end subroutine
   
end Module
