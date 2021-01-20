Module ModelParam_Boundary

   use ModelParamToolbox
   use ModelDimensions

   implicit none

 
      ! Parameter declarations as elements in structure p_PHYS,
      ! located in common block Common_Param_Physics
      ! declared in ModelParam_Physics.inc

      LOGICAL FIXTMP
           ! Controls how surface temperature at boundary is set

      real*4 BOUND_INFL(dimMBI)
           ! Adjusts relative weight of different basins (multiplied
           ! by areas) for propagating the distribution of nutrients
           ! from internal basins to outer basins


         ! ----------- Factors for adjusting external N, P and Oxygen:

      real*4 BndFac_N
           ! Total nitrogen, factor on excess of 100 µg/l

      real*4 BndFac_P
           ! Phosphorus, factor on value

      real*4 BndFac_Ox
       	  ! Oxygen, factor on deviation from saturation


      real*4 EXTBIO(2)
           ! Controls if biological components are included in 
           ! inflow from boundary areas into model basins
           
      real*4 TIDFAC
         ! degree of normal tidal variation (=0: no variations)


contains

   subroutine ParamGroup_Boundary

          ! Initiate parameters with default values, and register them
          ! for actions (user editing, file input, file storage)


      call paramGroup("Boundary conditions","Boundary")

      call ParamDeclLogical(FIXTMP,.true.,"FIXTMP", &
                "Controls how surface temperature at boundary" // &
                " is set during simulation:" // &
                " .true. : use specified boundary values only." // &
                " .false.: set equal to weighted average" // &
                " of inner basins")

      call ParamDeclRealArray(BOUND_INFL, (/1.0/),"BOUND_INFL","", &
                "Adjusts relative weight of different basins (multiplied" // &
                " by areas) for propagating the distribution of nutrients" // &
                " from internal basins to outer basins")

         ! Factors for adjusting external N, P and Oxygen

      call ParamDeclReal(BndFac_N, 1.0,"BndFac_N","", &
                "Factor on total nitrogen inputs " // &
                "(only affects excess of 100 µg/l)")

      call ParamDeclReal(BndFac_P, 1.0,"BndFac_P","", &
                "Factor on total phosphorus inputs")

      call ParamDeclReal(BndFac_Ox, 0.9,"BndFac_Ox","", &
                "Factor on oxygen levels in inputs," // &
                " (applies to deviation from saturation")


      call ParamDeclRealArray(EXTBIO, (/0.9,0.0/),"EXTBIO","", &
                " Controls to what degree biological components are included" // &
                " in inflow from boundary areas into the model basins." // &
                " Specify in range 0...1 for continuous variation between:" // &
                " EXTBIO(1): = 0: nutrients are imported in inorganic form" // &
                " = 1: nutrients are distributed between inorganic " // &
                " and biomass forms as in the model basins" // &
                " EXTBIO(2): = 0: no DOC, =1: ratio DOC/Phytoplankton" // &
                " as inside the model basins")


      call ParamDeclReal(TIDFAC,1.0,"TIDFAC","", &
              "Adjustment factor for normal tidal variation" // &
              " of boundary surface level (=0: no variations)" )

      end subroutine
     
end Module  
