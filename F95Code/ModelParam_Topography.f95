Module ModelParam_Topography

  use ModelParamToolbox
  use ModelDimensions    

  implicit none

      ! Parameter declarations for topography

         ! Used by subroutine TOPO for controlling depth division:
                     
         real*4  DDPMIN
               ! minimum depth of layers (at surface)

         real*4  DPFACT
               ! factor of increasing thickness of deeper layers

         INTEGER NLIMAX
               ! upper limit to number of layers NLI,
               ! can be set at run-time to vary vertical resolution
               ! without recompiling the program

         INTEGER TOPO_NR
               ! May be used to select between 
               ! alternative topographies in topography input file


  contains! *******************************************************************


    subroutine ParamGroup_Topography

           ! Start parameter Group:

      call paramGroup("Topography setup parameters","TOPOGRAPHY")

           ! Initiate and declare parameters:

      call ParamDeclReal(DDPMIN, dimDDPMIN_value, "DDPMIN", "(m)", &
                "Minimum depth of layers (at surface)")

      call ParamDeclReal(DPFACT, dimDPFACT_value, "DPFACT", "", &
                "Thickness increase factor for deeper layers")
      
      call ParamDeclInteger(NLIMAX, dimNLIMAX_value, "NLIMAX", "", &
                "Upper limit to number of layers NLI,"// &
                " can be set at run-time to vary vertical resolution"// &
                " without recompiling the program ")
      
      call ParamDeclInteger(TOPO_NR, 1, "TOPO_NR", "", &
                "Option for alternative topographies")

    end subroutine

end Module ModelParam_Topography
