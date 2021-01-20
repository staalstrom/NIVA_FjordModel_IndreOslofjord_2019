Module ModelVar_Topography

   use ModelDimensions
   implicit none

  !  --------------------------- Topography ---------------------------

      integer ND
      real*8 DEPTH( dimMD )

   ! ----------- Internal basins ------------
      integer NBI
      real*8 VTOTZ  (dimMBI)  ! Initial total volume
      real*8 ZBOTMI (dimMBI)  ! Maximum depth in basins
      real*8 LSHORE (dimMBI)   ! Length of shoreline (for benthos)
      integer NLVOPN( dimMBI )  ! Number of layers in each basin
                                ! not closed off by sills
      real*8 VFROPN( dimMBI )  ! Volume fraction of NLVOPN open layers
      integer INDXI ( dimMBIplus1 )  ! Index limits to layers
      real*8 AREA ( dimMLI ), BOTTOM( dimMLI), VFRAC ( dimMLI )

   ! ----------- External basins ------------
      integer NBE
      integer INDXE ( dimMBEplus1 )
      real*8 ZBOTME ( dimMBE )

! Names of external boundaries (-NBE:-1) and internal basins (1:NBI):
      CHARACTER*40 BasinName(-dimMBE:dimMBI)

  ! --------------- Connections ------------
      integer NC
      real*8 WVDIR (dimMC)
            ! WVDIR: direction of connection (from basin A to B):
            !        0: north south, 90 east west, etc.

      real*8 ZSill (  dimMC  )
      integer INDXC (  dimMCplus1 )
      INTEGER BConn1( dimMC  ), BConn2( dimMC  )
      real*8 TCVBUF(2,dimMC)
      real*8 WIDTH ( dimMLC )  ! mean width of each open layer
      real*8 VBUFMX(2,dimMLC) , VBUFTR(2,dimMLC)


      integer NLI_lim
      integer NLI, NLE, NLC
      
  ! Controls if new topography setup is needed:
 	   INTEGER TOPO_OLD	! remembers last topography used
					         ! in case topography alternative has been
					         ! reselected between runs.
    
      integer NLMPRV   ! previous value of NLIMAX
      real*8 DPFOLD      ! previous value of DPFACT

end Module ModelVar_Topography 