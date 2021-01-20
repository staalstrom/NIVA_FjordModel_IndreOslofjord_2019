Module ModelParam_Mussels

  use ModelParamToolbox
  use ModelDimensions

  implicit none
 
      ! Parameter declarations for mussel modelling


      real*4 PCMUSL, NCMUSL    ! P:C & N:C fixed weight rations in mussel

                             ! Time constants (days):
      real*4 TMSPWN            ! for spawning (CTMUSL-CAMUSL)
      real*4 TMSETL            ! for settling (available settling area)

                             ! Grazing coefficients for:
      real*4 GRMFYT(dimMFYTG) !   Phytoplankton
      real*4 GRMBCT            !   Bacteria
      real*4 GRMZOO            !   Zooplankton

      real*4 MCFMIN            ! Threshold food concentration
                             !  only particulate food counts
  
      real*4 MUSLDR(4)         ! Mussel mortality as rates pr. year
						           ! (1): proportional to excess vs. area capacity
						           !      (=rate at 100% excess)
						           ! (2): intrinsic rate for age class zero
						           ! (3): intrinsic rate for all ages >0
						           ! (4): rate increasing with high age (inversely
						           !      proportional to years left to MSAGMX)
						           ! Total rate is sum of terms 1 + (2 or or 3) + 4,
						           ! but with lower limit given by 
                             ! size growth in each class

      real*4 MSAGMX            !  Age where rate (4) applies

      real*4 MFWFAC            ! Factor for friction velocity
                             ! as measure of circulation velocity
                             ! when calculating food availability
                             ! for mussels

      real*4 MFILTM (dimMBI)  ! At least this fraction of filtered volume is
                             ! considered new, i.e. from main water body

      real*4 MXdetr            ! Fraction of excreted materials
                             ! entered into detritus fraction

  ! -------------- Specific data on blue mussel individuals:

                             ! Filtering:
         
      real*4 MSVC              ! Max. filtering capacity (liter/hour) 
                             ! for individual with soft tissue dry weight MSINDW(1)

      real*4 MSINDW(2)         !  Critical soft body dry weights in gram:
							        !  (1): limit Wc between lower and upper range of 
							        !       weight function for filtering and respiration.
							        !  (2): maximum weight Wm in filtering relation
           
      real*4 MSWR              ! Weight Wr (gram dry weight)
                             ! with all production going to reproduction

                             ! Exponents of weight relation:
      real*4 MSQW (2)          ! on filtering
      real*4 MSBW (2)          ! on respiration

      real*4 MSERMX            ! Upper limit to fraction of net growth
                             ! used for reproduction.
      
      real*4 MSREXP            ! Exponent of weight dependence
                             ! for reproductive effort

      real*4 MSEASS(3)         ! Maximum ingestion efficiency for carbon:

      real*4 MSCREQ            ! Food concentration where unrestricted
                             ! effective filtering equals physiological
                             ! needs for ingested material for mussel
                             ! with softbody dry weight of 1 gram

      real*4 MSCWXP            ! exponent in weight dependence

      real*4 MRSP15            ! Starving respiration for individual
                             !  with dry soft tissue weight MSINDW(1)
                             ! at temperature 15 deg.C

      real*4 MTRESP            ! Temperature coefficient for exponential
                             ! variation of respiration with temperature:

      real*4 MRASSF            !  Additional respiration as fraction of tissue buildup:


  contains

      subroutine ParamGroup_Mussels

     !  Initiate parameters with default values, and register them
     !  for actions (user editing, file input, file storage)


           ! Start parameter Group:

      call paramGroup("Mussel parameters","MUSSELS")


           ! Initiate and declare parameters:

      call ParamDeclReal( PCMUSL, 0.027, "PCMUSL", "weight:weight", &
                "Fixed P:C ratio in mussels") 

      call ParamDeclReal( NCMUSL, 0.18, "NCMUSL", "weight:weight", &
                "Fixed N:C ratio in mussels")

      call ParamDeclReal( TMSPWN, 15.0, "TMSPWN", "days", &
                "Time constants for spawning") 

      call ParamDeclReal( TMSETL,  2000.0, "TMSETL", "", &
                "Time constants for larvae spawning") 

      call ParamExpl("Coefficients for efficiency of mussel grazing:"// &
                     " Multiplication factors [0...1] for different food types") 

      call ParamDeclRealArray( GRMFYT, (/1.0/), "GRMFYT", "", &
                "Mussel grazing coefficients for phytoplankton") 

      call ParamDeclReal( GRMBCT, 0.1, "GRMBCT", "", &
                "Mussel grazing coefficients for bacteria") 

      call ParamDeclReal( GRMZOO, 0.2, "GRMZOO", "", &
                "Mussel grazing coefficients for zooplankton") 

      call ParamDeclReal( MCFMIN,  20.0, "MCFMIN", "mg C/liter", &
                "Threshold food concentration;"// &
                " only particulate food counts")

      call ParamDeclRealArray( MUSLDR,  (/1.0, 1.0, 0.5, 1.0/), &
                                   "MUSLDR", "fraction per year", &
                "Mussel mortality:"// &
                "( 1): proportional to excess population"// &
                " (= rate at 100% excess over area capacity)"// &
                " (2): intrinsic rate for age class zero"// &
                " (3): intrinsic rate for all ages >0 "// &
                " (4): rate increasing with high age (inversely "// &
                " proportional to years left to MSAGMX)"// &
                " Total rate is sum of terms 1 + (2 or 3) + 4,"// &
                " but also with a lower limit based on" // &
                " shell size growth in each class")


      call ParamDeclReal( MSAGMX, 10.0, "MSAGMX", "years", &
                "Age where rate (4) applies") 

      call ParamDeclReal( MFWFAC,  1.0, "MFWFAC", "", &
                "Factor for friction velocity of wind"// &
                " as measure of circulation velocity in basin"// &
                " when calculating food availability for mussels. "// &
                " Parameterisation of exchange between"// &
                " main water body and water along shoreline")

      call ParamDeclRealArray( MFILTM,  (/2.0/), "MFILTM", "", &
                "At least this fraction of filtered volume"// &
                " is considered new, i.e. from main water body"//&
                " and not recycled from previously filtered water")

      call ParamDeclReal( MXDETR,  1.0, "MXDETR", "", &
                "Fraction of excreted materials"// &
                " entered into detritus fraction")

  ! ------------------ Specific data on blue mussel individuals:

      call ParamExpl("Critical mussel weights for size dependence"// &
                     " of filtering and respiration in mussels:")
 
      call ParamDeclRealArray( MSINDW,  (/0.007, 0.35/), &
                                   "MSINDW", "gram dry weight", &
                "(1): limit Wc between lower and upper range"// &
                " of weight dependence of filtering and respiration;"// &
                " (2): weight Wm when reaching maximum filtering level")

      call ParamDeclReal( MSVC,  0.2, "MSVC", "litres/hour", &
                "Max. filtering capacity Vc for individual"// &
                " of soft tissue dry weight MSINDW(1)")

      call ParamDeclReal( MSWR,  3.0, "MSWR", "gram dry weight", &
                "Individual weight Wr when reaching maximum reproduction")    ! 

      call ParamDeclRealArray( MSQW,  (/1.0, 0.667/), "MSQW", "", &
                "Exponents of weight relation on filtering") 

      call ParamDeclRealArray( MSBW,  (/1.333, 0.667/), "MSBW", "", &
                "Exponents of weight relation on respiration")   ! 

      call ParamDeclReal( MSERMX,  0.9, "MSERMX", "", &
                "Upper limit to fraction of net growth" // &
                " used for reproduction") 

        call ParamDeclReal( MSREXP,  0.33, "MSREXP", "", &
                "Exponent of weight dependence for reproductive effort") 

      call ParamDeclRealArray( MSEASS,  (/0.6, 0.8, 0.8/), &
                "MSEASS", "", &
                "Maximum ingestion efficiency for carbon") 

      call ParamExpl("Food concentration where unrestricted " // &
                     "Effective filtering equals physiological needs" // &
                     " for ingested material")

      call ParamDeclReal( MSCREQ,  300.0, "MSCREQ", "mgC/m3", &
                "For 1 individual with 1 gram dry weight of soft body") 

      call ParamDeclReal( MSCWXP,  0.15, "MSCWXP", "", &
                "Exponent in weight dependence") 

      call ParamDeclReal( MRSP15,  0.0000086, &
                             "MRSP15", "liter O2/h", &
                "Starving respiration for individual"// &
                " of dry soft tissue weight MSINDW(1) at 15 deg.C")       ! 

      call ParamDeclReal( MTRESP,  0.065, "MTRESP", "", &
                "Temperature coefficient for exponential"// &
                " variation of respiration with temperature") 

      call ParamDeclReal( MRASSF,  0.15, "MRASSF", "", &
                "Additional respiration as fraction of tissue buildup") 


      end subroutine
     
end Module     
