Module ModelParam_Decomposition

  use ModelParamToolbox
  use ModelDimensions

  implicit none

      ! Parameter declarations as elements in structure p_DC,
      ! located in common block Common_Param_Decomposition
      ! declared in ModelParam_Decomposition.inc

	   
      real*4 OXCFAC    ! Correction factor for ratio oxygen:carbon ratio
                     ! in primary production and decomposition,
                     ! used when starting model run only, to modify
                     ! stoichiometric factors, see subroutine STOICH.FOR

      real*4 ODMRAT    ! Realization rate of oxygen demand ODM

      real*4 DGDETZ    ! Half saturation depth for degradation of detritus
                     ! assumed to reach 50% of full degr. at depth DGDETZ

      real*4 DGRATE(3) ! max specific rates pr. day for oxic degradation
                     ! of fresh material at 20 degC

      real*4 ACCLRC    ! Factor aR (>=0) in function for progressive"
                     ! decrease of degradability

      real*4 ACCLXP    ! Exponent betaR (>=0) in function for
                     ! decrease of degradability

      real*4 DGWSF     ! Reduction factor of oxygen and nitrate concentration
                     ! for bottom degradation kinetics relative
                     ! to pelagic detritus

      real*4 DGNFAC, DGPFAC, DGSFAC
                     ! Relative levels of N, P and Si specific
                     ! remineralisation rates relative to C
                     ! (may change composition during degradation)

      real*4 DGCMAX(2) ! Maximum limits to absolute rates of degradation of carbon:
                     !    1: detritus in water mg/m3/day
                     !    2: Sediment matter   mg/m2/day

      real*4 DOXBRR    ! Bacterial fractions of oxic degradation:
					      !  max. bacterial rate independent of fauna
					      !  Note: originally array(1:2), with
					      !  (2) = bacterial fraction of fauna-related 
                     !        degradation, but not used
                     !        (1.0 implicit in code)

      real*4 DOXKB, DOXKM
      					!  Half-saturation concentration for oxic degr. (ml/l)
      					!  for bacterial and fauna contribution:


      real*4 DOXKS
                     ! inverse of oxygen. debt in sediment giving 
                     ! 50% reduced fauna-related degradation! 


      real*4 DOXLIM    ! Lower oxygen limit for bottom fauna


   !  ------------ denitrification -----------------
      
      real*4 DNITRR    ! Maximum fraction of anoxic decomposition
                     ! by denitrification

      real*4 DNITKS    ! Half-saturation NO3N concentration
                     ! for denitrification by external NO3
                     ! when there is no oxic degradation
 
      real*4 KOXN      ! Maximum increase of half-saturation
                     ! NO3N concentration for denitrification
                     ! due to oxic sone as transport buffer

      real*4 DNITF     ! (K_nitr, dimensionless)
                     ! Half saturation constant for response of denitrification  
                     ! to degree of oxic bacterial degradation:
                     ! Low value means that even small oxic activity
                     ! (low oxygen conc.) will nitrify released
                     ! ammonium efficiently

      real*4 DNOXFR    ! Max. degree of denitrification of ammonium
                     ! released by oxic degradation

      real*4 DNITXP    ! Exponent for monod function using DNITF, see DGRADE

      real*4 FDNH3     ! Degree of direct removal of remineralized nitrogen
                     ! as part of denitrification:
      

    !  ----------- Sulphide reduction: --------------

      real*4 SULFRR    ! Relative rates compared to bacterial oxic degr.

      real*4 SULFOX(2) ! (1): Upper limit and (2): half saturation
                     ! for oxygen (debt) concentration in water
                     ! in transition to maximal rate SULFRR

      real*4 sulfxp    ! Exponent for response of sulphate reduction rates"
                     ! to oxygen concentrations


        ! Ammonium nitrification of free ammonium discharged
        ! & ammonium mineralized from sinking matter:

      real*4 RAMMOX    ! Max. NH4 specific rate (1/day)
      real*4 KAMMOX    ! Oxygen half saturation concentration

    !  ------------ Sulphide buffering in sediments ----------

      real*4 ASEDMX    ! Limit for amount of sulphide in sediments

      real*4 ASEDLR(2) ! Sulphide leakage rates (1/year)
                     ! below and above limit

      real*4 ASEDOX    ! Factor with unit m for converting 
                     ! oxygen levels in water (ml/liter = liter/m3)
                     ! to units of oxygen debt in sediment (liter/m2)
                     ! for calculation of the oxygen gradient 
                     ! That drives leakage to water.
                     ! ( Represents effective depth scale for
                     ! active sulphide buffer in sediment)

      real*4 ASOXTL    ! Factor for converting from ASED leakage (literO2/m2/day)
                     ! to contribution to reduced effective oxygen 
                     ! concentration in sediment (literO2/m3)


  ! ################################################################
  
      real*4 BURIAL (dimMBI)  
                     ! General burial or disappearance rate of sediment
                     ! ( specific rate 1/year for active sediment layer,
                     !   = inverse residence time in active layer).

      real*4 GMX20B    ! Max. spec. growth rate of bacteria at 20C

      real*4 BTRESP    ! Coefficient for temperature dependence

      real*4 TTURNB    ! minimum turnover time for DOC, N and P pools
                     ! by bacteria

      real*4 BACDET    ! Transit rate from free bacteria to detritus
                     ! (could have two-way transit at different rates)

      integer LDGRV  ! Layer nr. for saving variables in debug array
                     ! real*4 DGRV(20,2), see subroutine DGRADE



				!************************
				! Inorganic phosphorus:
				!************************

      real*4 PPAMAX    ! Max. precipitation rate of P
                     ! in mg/m2/day for OXYG >= PPOXMX

      real*4 PPRMAX    ! Max. precipitation speed of P in m/day

      real*4 PPOXEX    ! Exponent for dependence of precipitation
                     ! on oxygen concentration below PPOXMX
    
      real*4 PPOXMX    ! oxygen limit (ml/l) for phosphorus precipitation,
                     ! see above

      real*4 PADRET    ! Fraction of remineralized P retained
                     ! under oxic conditions

      real*4 PADMAX    ! Max P retained in sediment (mg/m2)

      real*4 PADRLS    ! Relative rate of release of excess P
                     ! retained in sediments (1/day)

      real*4 PADASD    ! Sulphide content giving max. release rate for
                     ! P buffered in sediments.

		real*4 PSBURF     ! Ratio between burial rate for phosphorus
                     ! and organic matter in sediment

   

  contains

! -------------------------------------------------------------

      subroutine ParamGroup_Decomposition

			! Initiate parameters with default values
         ! and register them for actions
         ! (user editing, file input, file storage)

      call paramGroup("Decomposition of organic matter","DECOMPOSITION")


      call ParamDeclReal(OXCFAC,1.2, "OXCFAC", "", &
                 "Adjusting factor for ratio oxygen:carbon ratio" // &
                 " in primary production and decomposition" // &
                 " relative to default values in model description." // &
                 " Used only when starting model run to modify" // &
                 " stoichiometric factors OX_C, and NITR_C," // &
                 " see model code for further details.")

	   call ParamDeclReal(ODMRAT,0.0, "ODMRAT", "per day", &
                 "Realization rate of oxygen demand ODM")

	   call ParamDeclReal(DGDETZ, 0.1,"DGDETZ", "(m)", &
                 "Half saturation depth for degradation of detritus " // &
                 "from water surface")

      call ParamDeclRealArray(DGRATE, (/0.3, 0.05, 0.01/), &
                                    "DGRATE", "(1/day)", &
                 "max. specific rates pr. day for oxic degradation " // &
                 " of organic material at 20 degC "//  &
                 " (1): pelagic components (plankton)" // &
                 " (2): dead mussels " // &
                 " (3): terrestrial organic carbon input")

      call ParamDeclReal(ACCLRC, 2.0, "ACCLRC", "", &
                 "Factor aR (>, 0) in function for progressive decrease " // &
                 "of degradability of residual organic matter")
       
	   call ParamDeclReal(ACCLXP, 0.0, "ACCLXP", "", &
                 "Exponent betaR (>, 0) in function " // &
                 "for progressive decrease degradability " // &
                 "of residual organic matter")
       
	   call ParamDeclReal(DGWSF, 0.5, "DGWSF", "", &
                 "Reduction factor of oxygen and nitrate concentration " // &
                 "for bottom degradation kinetics " // &
                 "relative to pelagic detritus" )

	   call ParamDeclReal(DGNFAC, 1.0,"DGNFAC", "", &
                 "Ratio between specific remineralisation rates for N and C")
                 
	   call ParamDeclReal(DGPFAC, 1.0,"DGPFAC", "", &
                 "Ratio between specific remineralisation rates for P and C")

	   call ParamDeclReal(DGSFAC, 2.0,"DGSFAC", "", &
                 "Ratio between specific remineralisation rates for Si and C")

      call ParamDeclRealArray(DGCMAX, (/ 1000., 10000.0 /),"DGCMAX","", &
                "Maximum limits on absolute rates of degradation of carbon " // &
                " (1): for detritus in water (mg/m3/day) " // &
                " (2): Sediment matter mg/m2/day")

      call ParamDeclReal(DOXBRR, 0.3, "DOXBRR", "", &
              "max. relative bacterial degradation rate independent of fauna")
			        ! Note: originally array (1:2), (2) not used

      call ParamDeclReal(DOXKB, 0.2,"DOXKB", "ml/l", &
              "Half-saturation concentration for oxic degr. (ml/l); " // &
				  "bacterial contribution")
      
      call ParamDeclReal(DOXKM, 0.5,"DOXKM", "ml/l", &
      			"Half-saturation concentration for oxic degr. (ml/l); " // &
      			"macro-fauna contribution")


      call ParamDeclReal(DOXKS, 0.05,"DOXKS", "(m2/liter O2)", &
               " inverse of oxygen debt in sediment giving" // & 
               " 50% reduced fauna-related degradation")

      
      call ParamDeclReal(DOXLIM, 0.1, "DOXLIM", "ml/l", &
      			"Lower oxygen limit for bottom fauna")


      call ParamExpl("--------- Denitrification")
   
      call ParamDeclReal(DNITRR, 1.0,"DNITRR", "", &
               "Maximum fraction of anoxic decomposition " // &
               "by denitrification")

      call ParamDeclReal(DNITKS , 80.0 ,"DNITKS", "µgN/l", &
                 "Half-saturation NO3N concentration " // &
                 "for denitrification by external NO3" // &
                 "when there is no oxic degradation" )

      call ParamDeclReal(KOXN, 300.0 ,"KOXN", "µgN/l", &
                 "Maximum increase of half-saturation " // &
                 "NO3N concentration for denitrification " // &
                 "due to oxic sone as transport buffer" )

      call ParamDeclReal(DNITF, 0.1, "DNITF", "", &
                 "Half saturation constant for response of denitrification" // &
                 "to degree of oxic bacterial decomposition. " // &
                 "Low value means that even small oxic activity " // &
                 "(low oxygen) will nitrify released ammonium efficiently.")

      call ParamDeclReal(DNOXFR, 1.0, "DNOXFR", "", &
                 "Max. degree of denitrifaction of ammonium " // &
                 "released by oxic degradation")

      
      call ParamDeclReal(DNITXP, 0.333, "DNITXP", "", &
                 "Exponent for monod function using DNITF," // &
                 "(refer model description)")

      call ParamDeclReal(FDNH3 , 1.0, "FDNH3", "", &
                 "Degree of direct removal of remineralized nitrogen " // &
                 "as part of denitrification:")


      call ParamExpl("---------- Sulphide reduction:")

      call ParamDeclReal(SULFRR , 1.0, "SULFRR", "", &
                 "Relative rate compared to bacterial oxic decomposition")

      call ParamDeclRealArray(SULFOX , (/0.05,0.2/), &
                                "SULFOX", "ml O2-equiv./l", &
                 "(1): Upper limit and (2): half saturation" // &
                 "for oxygen equivalent concentration in water" // &
                 "in transition to maximal rate SULFRR")

      call ParamDeclReal(SULFXP, 2.0, "SULFXP", "", &
                 "Exponent for response of sulphate reduction rates" // &
                 "to oxygen concentrations")

      call ParamExpl(" Nitrification of ammonium " // &
                 "discharged in runoff or mineralized " // &
                 "from sinking matter:")

      call ParamDeclReal(RAMMOX, 1.0, "RAMMOX", "per day", &
                 "Maximum NH4 specific rate 1/day")

      call ParamDeclReal(KAMMOX, 0.2, "KAMMOX", "ml/l", &
                 "Oxygen half saturation concentration")

      call ParamExpl("--------- Sulphide buffering in sediments")

      call ParamDeclReal(ASEDMX, 300.0, "ASEDMX", "liter O2 /m2", &
                 "Critical value of amount of sulphide in sediments.")

      call ParamDeclRealArray(ASEDLR, (/0.05,0.5/), "ASEDLR", "per year", &
                 "Sulphide leakage rates (1): within critical level " // &
                 "and (2): for excess above critical level.")

      call ParamDeclReal(ASEDOX, 0.1, "ASEDOX", "m", &
                 "Ratio between oxygen contents (litre/m2) in sediments " // &
                 "and oxygen levels (ml/l) in water. " // &
                 "(Represents a sort of depth scale " // &
                 "for active sulphide buffer in sediment)")

      call ParamDeclReal(ASOXTL, 300.0, "ASOXTL", "day/m", &
                 "Factor for converting ASED leakage (literO2/m2/day)" // &
                 "into contribution to reduced effective oxygen" // & 
                 "concentration in sediment (literO2/m3)" // &
                 "(represents thickness of sediment layer over diff. coeff.)")

      call ParamDeclRealArray(BURIAL, (/0.08/),"BURIAL", "per year", &
                 "General burial or disappearance rate of sediment " // &
                 "  (specific rate 1/year for active sediment layer " // &
                 "    = inverse residence time in active layer).")

                                  ! NOTE: All Burial rates set equal initially 


      call ParamDeclReal(GMX20B, 0.5, "GMX20B", "per day", &
                 "Max. spec. growth rate of bacteria at temp. 20C")

      call ParamDeclReal(BTRESP, 0.06, "BTRESP", "per deg.C)", &
                 "coeff. in temperature dependence: function " // &
                 "exp(BTRESP*(T-20))" ) 
    
      call ParamDeclReal(TTURNB, 0.25, "TTURNB", "days", &
                 "minimum turnover time for DOC, N and P pools by bacteria")

      call ParamDeclReal(BACDET, 0.05, "BACDET", "per day", &
                 "Transit rate from free bacteria to detritus")


		call ParamExpl( "--------- Sedimentation and release " // &
                 "of dissolved inorganic phosphorus:" // &
                 "(crude parameterisation of "// &
                 " processes related to particle sinking")

      call ParamDeclReal(PPAMAX, 0.8, "PPAMAX", "mg/m2/day", &
                 "Max. absolute precipitation rate of P " //&
                 "for OXYG >= limit PPOXMX (below)")

      call ParamDeclReal(PPRMAX, 0.5,  "PPRMAX", "m/day", &
                 "Max. effective precipitation speed of dissolved P " // &
                 "through adsorbtion to sinking particles")

      call ParamDeclReal(PPOXEX, 0.5, "PPOXEX", "", &
                 "Exponent for dependence of P precipitation on" // &
                 "oxygen concentrations below critical value PPOXMX")
                 
      call ParamDeclReal(PPOXMX, 5.0, "PPOXMX", "ml/l", &
                 "Critical value of oxygen conc. for P precipitation;" // &
                 "below this value precipitation are reduced" )

      call ParamDeclReal(PADRET, 0.2, "PADRET", "", &
                 "Fraction of remineralized P retained under oxic conditions")
                 
      call ParamDeclReal(PADMAX, 5000.0, "PADMAX", "mg/m2", &
                 "Maximum amount of phosphorus retained " // &
                 "in active part of sediment")
                 
      call ParamDeclReal(PADRLS, 0.01, "PADRLS", "per day", &
                 "Relative release rate of excess P retained in sediments")

      call ParamDeclReal(PADASD, 10.0, &
                               "PADASD", " mg/m2 as oxygen debt", &
                 "Sulphide content giving maximum release rate for" // &
                 " P buffered in sediments.")

		call ParamDeclReal(PSBURF,1.0, "PSBURF","", &
                 "Ratio between sediment burial rate for phosphorus" // &
                 " bound in sediments and burial rate for organic matter" // &
                 " and sulphide (oxygen debt)")

      
   end subroutine
   
end Module