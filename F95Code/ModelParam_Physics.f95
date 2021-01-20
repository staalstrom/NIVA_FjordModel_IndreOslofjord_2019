Module ModelParam_Physics

  use ModelParamToolbox
  use ModelDimensions

  implicit none

      ! Parameter declarations for physical processes

      integer vdindx
            ! Controls check on volume derivative balance
            ! vdindx =1: from surface iteration
            !        =2: from final transport calc.

      real*4 DPEFF(dimMC)
			   ! Energy efficiency in horizontal transports
			   ! driven byhorizontal pressure gradients
            ! fraction of potential energy converted to 
			   ! effective kinetic energy

      real*4 HTRMIX(dimMC)
			   ! Degree of mixing between contiguous
            ! horisontal transports in same direction
			   ! 0 = no mixing, 1 = full mixing,
			   ! other values: portion of transport being mixed



   !  Vertical mixing, driven by tidal energy
   !  (Now includes variation with tidal energy input)

      real*4 N2SCAL
            ! Stability for which mixing is specified

      real*4 MIXCF(dimMBI)
            ! Mixing coefficients at stability N2SCAL

      real*4 MIXEXP
            ! Exponent alpha in eq. K=C*N**(-alpha)

      real*4 N2LIM
          ! Lower limit to effective stability in formula for K
          ! i.e. upper limit to mixing. Unit (1/s2)
          ! Based on actual observations in inner Oslofjord,
          ! set to avoid diverging mixing coeff., without influencing
          ! normal situations

      real*4 SFMIXC(dimMBI)
          ! Surface mixing constant (m2/s3), surface value

      real*4 SFMIXZ(dimMBI,2)
          ! Depth scale in exp. dampening:
          !   (SFMIXZ(i,1) : depth with approx. constant mixing 
          !   (SFMIXZ(i,2) : depth constant in exponential dampening

  ! ################################################################


      LOGICAL TROFF
          ! Turns off water transports (branch in TRANSP)

      LOGICAL HTROFF
          ! Turns off only horisontal transport calculations

      Logical ITRZ
          ! Controls surface iteration in subroutine TRANSP

      real*4 MIXFAC 
           !  Multiplicator for specified vertical diffusion

      real*4 GMIXFR(dimMBI)
           ! Fraction of released potential energy in deep inflows
           ! that is used against mixing as work against gravitation
           ! TO BE SET EMPIRICALLY

      real*4 GMIXDC
           ! vertical specific reduction rate of gravitational energy
           ! at stability BW_FREQ = 1.0 (per meter)

      real*4 GMIXDX
           ! Stability dependence in vertical reduction
           ! exponent for BW_FREQ


      real*4 WVFAC (dimMC)
      real*4 WVHMIN(dimMC)
            ! WVFAC: ratio between wind speed and surface current
            ! WVHMIN: minimum depth for velocity distribution
            !         (triangular distribution)


      !  Update oxygen derivative with flux through surface layer
      !  and effect of bubble formation during supersaturation
      real*4 OXBUBL
                ! Specific reduction rate of oxygen super saturation
                ! by bubbles, accounted for as export from system

      real*4 OXSFAC
                ! Increase/decrease surface exchange

!***********************************************************************
! Air/water interface conditions:

      real*4 DAYDIV
       ! Split daylight period in at least DAYDIV integration periods

      real*4 CDFAC
                ! Factor for wind friction
      real*4 CEFAC
                ! Factor for evaporation & heat exchange,
                ! if =0, turns off all heat transfer terms

      real*4 IRFRAC
                ! Fraction of IR-radiation in light
      real*4 ICEFAC
                ! reduction factor for light
                ! in possible ice-cover conditions

      real*4 RADFAC(3)
               ! Variation coefficients for heat/light effect of:
               !  1: visual direct solar radiation
               !  2: visual diffuse radiation
               !  3: infra-red radiation

      real*4 ATTNCF(2)
               ! LIGHT ATTENUATION COEFFICIENTS:
               !  UNITS: (1): 1/m, (2): 1/(m*(µgC/l))


  contains

    subroutine ParamGroup_Physics

          ! Initiate parameters with default values, and register them
          ! for actions (user editing, file input, file storage)

      call paramGroup("Physical processes","PHYSICS")

      call ParamDeclInteger(vdindx, 1, "vdindx", " ", &
                 "Controls check on volume derivative balance " // &
                 " =1: from surface iteration" // &
                 " =2: from final transport calc.")

      call ParamDeclRealArray(DPEFF,(/0.5/), "DPEFF"," ", &
                 " Energy efficiency in horizontal transports" // &
                 " driven by horizontal pressure gradients;" // &
                 " fraction of potential energy converted to" // &
                 " effective kinetic energy")

      call ParamDeclRealArray(HTRMIX, (/0.5/), "HTRMIX"," ",&
                 "Degree of mixing between contiguous" // &
                 " horisontal transports in same direction" // &
                 " 0 = no mixing, 1 = full mixing," // &
                 " other values: portion of transport being mixed")

      call ParamExpl( " Vertical mixing, driven by tidal energy:")

      call ParamDeclReal(N2SCAL, 0.000063, "N2SCAL", "1/s2", &
                 "Stability (BW-frequency squared)" // &
                 " for which mixing is specified" // &
                 " (default 0.000063 = 1/(10**4.2)")

      call ParamDeclRealArray(MIXCF, &
                 (/1e-6,1e-4,1e-5,1e-5,1e-5,1e-6,1e-6,1e-5,1e-5/), &
                 "MIXCF", "m2/s", &
                 "Mixing coefficients at stability N2SCAL")
                 
      call ParamDeclReal(MIXEXP, 1.6,"MIXEXP","", &
                 "Exponent alpha in equation for" // &
                 " vertical mixing coefficient: K=C*N**(-alpha)")
      
      call ParamDeclReal(N2LIM, 1.0e-7, "N2LIM","1/s2", &
                 "Lower limit to effective stability" // &
                 " in formula for K, i.e. upper limit to mixing." // &
                 " set to avoid numerical overflow in mixing coeff." // &
                 " without influencing normal situations" )

      call ParamExpl("Surface mixing; surface values and"// &
                     " depth scales in exponential dampening:")
                     
      call ParamDeclRealArray(SFMIXC, (/1.0e-10/),"SFMIXC","m2/s3", &
                 " Mixing energy at the surface")

      call ParamExpl("Twodimensional array (" // IntToString(dimMBI) // ",2):")
      
      call ParamDeclRealMatrix(SFMIXZ, reshape((/8.0,4.0/), (/1,2/)), &
                   "SFMIXZ", "m", &
                   " (SFMIXZ(i,1): Thickness of layer with approx. constant mixing"// &
                   " (SFMIXZ(i,2): Depth constant in exponential dampening")

                  ! Default values set to keep energy 
                  ! approximately constant down to 8m,
                  ! and then with a reduction by approx.
                  ! a factor 10 for each 10-15m


      call ParamDeclLogical(TROFF, .FALSE., "TROFF", &
                "Turns off water transports (branch in TRANSP)")

      call ParamDeclLogical(HTROFF, .FALSE., "HTROFF", &
                "Turns off only horisontal transport calculations")

      call ParamDeclLogical(ITRZ  , .TRUE., "ITRZ", &
                "Controls surface iteration in subroutine TRANSP")

      call ParamDeclReal(MIXFAC , 1.0, "MIXFAC"," ", &
                "Multiplicator for specified vertical diffusion")

      call ParamDeclRealArray(GMIXFR,(/0.0/),"GMIXFR","", &
                "Fraction of released gravitational potential energy" // &
                " in sinking dense inflows giving vertical mixing" // &
                " as work against gravitation; must be set empirically")
                
      call ParamDeclReal(GMIXDC, 1.25, "GMIXDC","(per meter)", &
                "Vertical specific reduction rate of "// &
                " gravitational energy at stability BW_FREQ = 1.0")

      call ParamDeclReal(GMIXDX, 0.4, "GMIXDX", "", &
                "Stability dependence of vertical reduction," // &
                " exponent for BW_FREQ.")

      call ParamExpl("Coefficients for wind-driven transports " // &
                " across connections between basins:")

      call ParamDeclRealArray(WVFAC, (/0.03/),"WVFAC","", & 
                "Ratio between wind speed and surface" // &
                " wind-driven current across connections")

      call ParamDeclRealArray(WVHMIN, (/3.0/), "WVHMIN", "", &
                "Minimum depth range for winddriven surface current" // &
                " (in triangular distribution)")

      call ParamDeclReal(OXBUBL, 1.0, "OXBUBL","(per day)", &
                "Specific reduction rate of oxygen super-saturation"//  &
                " due to primary production" // &
                " (assumed to be bubbled to atmosphere")

      call ParamDeclReal(OXSFAC, 1.0, "OXSFAC","", &
                "Factor to adjust oxygen surface exchange" // &
                " relative to rates built into model") 

      call ParamExpl("Air/water interface conditions:")

      call ParamDeclReal(DAYDIV, 4.0, "DAYDIV","", & 
                "Split daylight period in at least this many"// &
                " integration time steps")

      call ParamDeclReal(CDFAC, 1.0,"CDFAC","", &
                "Factor for wind friction")
                
      call ParamDeclReal(CEFAC, 1.0," CEFAC", "", &
                "Factor for evaporation & heat exchange,"// &
                " if =0, it turns off all heat transfer terms" )

      call ParamDeclReal(IRFRAC,0.4, "IRFRAC", "",&
                "Fraction of IR-radiation in light energy")
  
      call ParamDeclReal(ICEFAC, 1.0, "ICEFAC","", &
               "Reduction factor for light at supposedly ice-cover conditions." // &
               " (i.e. surface temperature below freezing point)" )

      call ParamDeclRealArray(RADFAC,(/1.0/),"RADFAC","", &
               " Variation coefficients for heat/light effect of:" // &
               " (1): visual direct solar radiation" // &
               " (2): visual diffuse radiation" // &
               " (3): infra-red radiation")

      call ParamDeclRealArray(ATTNCF, (/0.25, 0.00025/), &
                               "ATTNCF","1/m and 1/(m*(µgC/l)", & 
               "Light attenuation coefficients:" // &
               "(1): Constant, default 0.25 [1/m] " // &
               "(2): Organic carbon dependence, unit [1/(m*(µgC/l)]" // &
               "     (default value 0.00025 = 0.003/12.)")

      end subroutine
     
end Module  
