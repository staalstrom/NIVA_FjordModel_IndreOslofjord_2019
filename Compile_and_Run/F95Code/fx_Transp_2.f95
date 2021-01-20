      Module fx_Transp_2

      use fx_Jetmix
      use fx_RunControl
      use fx_Trhitr
      use fx_Transp_h
      use fx_Transp_v

      implicit none
      
      INTEGER MQ_TR
            ! Connected to dived outlets (assumed to rise to surface):
C           ! Sum of possible depth entrainment transports for dived
C           ! outlets with outlet density < ambient density.
            ! <= dimMS*dimMLI
          parameter ( MQ_TR = dimMS*dimMLI )

      real*8 QJET_TRANSP     (MQ_TR)
C           - Vertical water transports due to dived outlets QJET

      INTEGER QJET_TRANSP_INDEX  (2,MQ_TR)
C           - Global index for source layer of each QJET_TRANSP


      
      contains

C ==================================================================   .
C Eutrophication model   - File:    TRANSP_2.FOR
C                                   Birger Bjerkeng, NIVA.

$define volume_version_2
        ! skips constant volume layer.


C  CONTROLS DEBUG DUMP OF CALCULATIONS ON FILE:
C    0 : NO DUMP
C    1 : DUMP END RESULTS
C    2 : DUMP INTERMEDIATE AND END RESULTS

$undefine DEBUG_SETUP
$undefine DEBUG_APPLY
$undefine DEBUG_TSTEP
$undefine DEBUG_vmix

$if defined DEBUG_SETUP
$define DEBUG_SETUP_GT_1
$define DEBUG_TRANH_TRANC_CALL
$endif

$IF DEFINED DEBUG_SETUP || DEFINED DEBUG_APPLY || DEFINED DEBUG_VMIX
$define DEBUG_ANY
$endif

C =====================================================================
C   Transport calculations and mass --> concentration calculations
C =====================================================================


C ====================================================================

C Contains subroutines:

C   TRANSP_SETUP is called from WATER_TRANSPORTS in Module TRANSP_1
C   and calculates transport terms.
C        Calls subroutines TRANH in Module fx_TRANSP_H.F95
C               TRANC and TRANV1 in Module fx_TRANSP_V.F95
C                JETMIX in module fx_Jetmix.f95

C =====================================================================

      SUBROUTINE TRANSP_SETUP ( HTR_CALC, DEBUG, ITRZ,
     &         DEPTH, MBI, NBI, VDYN, ZBOTMI,
     &         DT_QTRANSP_CALC, QWSURF, MS, BASINQ, DEPTHQ,
     &          QJET, QTEMP, QSPP, dDens_dSPP, RNFNDX, QDIAM, NHOLES,
     &          MQ_TR, XMIX,
     &          INDXI, NLI, AREA, SAL, TEMP, SPP, DENSI, VLAYER, VFRAC,
     &          NLVOPN, VFROPN, VLCORF,
     &         MIXFAC, MIXEXP, N2LIM,
     &         SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
     &         NBE, INDXE, NLE, DENSE, ZSURFE, DZDTX,
     &         NC, INDXC, NLC, BConn1, BConn2,
     &         ZSILL, WIDTH, DPEFF, HTRMIX,
     &         MLC, VBUF, VBUFMX, VBUFTR, TCVBUF, UFLOW,
     &         VTOTDV, VDYNDV, VBUFDV, MAX_DT,
     &  MLV, VTRANS, VMIX, MLH, HTRANS, HTR_L, MC, INDXH,
     &  NQ_TR, QJET_TRANSP, QJET_TRANSP_INDEX, RCQNDX,
     &  Ambient_Volume_Flux, Neutral_Depth,
     &  MW, WK_ARRAY, MT_TEMP ,TR_TEMP, TR_SOURCE_L,
     &  ML_C_MAX, CONNECT, MB, ZSURF_W, surf_change, BSFLUX,
     &  BWFREQ, ZMID )

!!      SUBROUTINE TRANSP_SETUP (HTR_CALC, 
!!     &  MLV, VTRANS , VMIX, MLH, HTRANS, HTR_L, MC, INDXH,
!!     &  RCQNDX,
!!     &  Ambient_Volume_Flux, Neutral_Depth,
!!     &  MW, WK_ARRAY, MT_TEMP ,TR_TEMP, TR_SOURCE_L,
!!     &  ML_C_MAX, CONNECT, MB, ZSURF_W, surf_change, BSFLUX,
!!     &  BWFREQ, ZMID )

      
C ----------------------- External arguments: ------------------------


C ------ In:
      LOGICAL HTR_CALC ! .TRUE. turns off calc. of horisontal transp.

      LOGICAL DEBUG !  .TRUE. triggers debug printout if DEBUG_SETUP>0:
      LOGICAL ITRZ  !  .true. :Iterate towards equilibrium surfaces.
C

      real*8 DEPTH(*)
C           Depth values defining layer division
      INTEGER MBI, NBI
C           Dimensional and actual number of internal basins
      real*8 VDYN(2, NBI)
C           Dynamic volume change from zero level (m3)
      real*8 ZBOTMI(NBI)
C              - bottom depth of basins

      real*8 DT_QTRANSP_CALC
C              - = 0 signals that QTRANSPO should be recalculated
      real*8 QWSURF(NBI)
C              - net water influx to surface layer from atmosphere
C                (unit m3/s)

C  ----------- Description of outlets:
      INTEGER MS            ! Number of outlets
      INTEGER BASINQ(MS)    ! Number of receiving basin
      real*4 DEPTHQ(MS,2)     ! Depth of outlet (m)
      real*8 QJET(MS,2)       ! Water flux (m3/s)
      real*8 QTEMP(MS)        ! Temperature flux (deC*m3/s)
      real*8 QSPP(MS)         ! flux of particles (g/s)
      real*4 dDens_dSPP       ! effect of particle conc. on density
      INTEGER RNFNDX(MS,2)  ! Global index of receiving layer
      real*4 QDIAM(MS)        ! Jet diameter (m) for submerged outlets
                            ! mixing depth for surface depth
      INTEGER NHOLES(MS)    ! Number of holes in outlet
      INTEGER MQ_TR         ! Max. number of jet induced transports

      real*8 XMIX(NBI)
C          - number of layers in well-mixed surface volume,
C            with possible fractional part for last layer.

      INTEGER INDXI(NBI+1)
C              - layer index limits for internal basins.
      INTEGER NLI
C           Number of layers in internal basins
      real*8    AREA(NLI)
C              - Horizontal area above layers in internal basins (m2)
      real*8    SAL    (NLI)
C              - Salinity of internal layers
      real*8    TEMP   (NLI)
C              - Temperature of internal layers
      real*8    SPP    (NLI)
C              - Particle content of internal layers
      real*8    DENSI  (NLI)
C              - Density (sigma-t) of internal layers
      real*8    VLAYER (NLI)
C              - Basin volume in each layer.
      real*8    VFRAC  (NLI)
C              - Fraction of total basin volume in each layer
C                at initial surface position
      INTEGER NLVOPN(NBI)
C              - Number of open layers in each basin
      real*8    VFROPN (NBI)
C              - Normal ä fraction of open layers
      real*8    VLCORF(NBI)
C              - correction factor for open volumes (above some sill)
      real*4    MIXFAC
      real*4    MIXEXP, N2LIM
      real*4    SFMIXC(MBI), SFMIXZ(MBI,2)
C              - Mixing parameters in each internal basin
C                used in subroutine VTRAN1.
      real*4    GMIXFR (NBI)
      real*4    GMIXDC
      real*4    GMIXDX

      INTEGER NBE
C              - actual number of external basins
      INTEGER INDXE(NBE+1)
C              - layer index limits for external basins.
      INTEGER NLE
C           Number of layers in external basins
      real*8    DENSE(NLE)
C              - Density (sigma-t) of layers in external basins.
      real*8    ZSURFE(NBI), DZDTX(NBI)
C              - Surface level of each external basin (high level >0),
C                with time derivative (m/day).

      INTEGER NC
C              - actual number of connections.
      INTEGER INDXC(NC+1)
C              - width array index limits for connections.
      INTEGER NLC
C              - Number of layers in connections
      INTEGER BConn1(NC), BConn2(NC)
C              - Connected basin numbers.
C                +n: internal basin
C                -n: external basin
      real*8    ZSILL(NC)
C               - Maximum depth of connections

      real*8    WIDTH(NLC)
C              - Transport widths of connected layers
C                Values from I1=INDXC(IC)+1 through I2=INDXC(IC+1)
C                specifies width for layer 1 to I2-I1+1
C                in connection number IC.

      real*4    DPEFF(NC)  ! Pressure efficiency; fraction of
                         ! energy converted into kinetic energy
      real*4 HTRMIX(NC)
C           Degree of mixing between homogenized inflows.

      INTEGER MLC
C              - Dimensioning number of layers in connections
      real*8    VBUFMX(2,MLC), VBUFTR(2,MLC)
C           Maximum limit and transition zone for buffer volume:
C           Before becoming efficient in mass transfer,
C           outflows must empty volume if filled with water from
C           outflowing basin, and inflows must fill the volume
C           with water from other basin. Cfr. TRANH.
      real*8    VBUF(2,MLC)
C           Current status of buffer volume VBUFMX:
C                 Part of volume filled with water from other basin.
      real*8    TCVBUF(2,NC)
C           Mixing time constant (days):
C                 VBUF/TCVBUF mixes pr. time unit


C ------ Out:
      real*8 UFLOW  (MLC)
C              - flow velocity, exported as user info. to main model
      real*8 VTOTDV (MBI)
C              - Time derivative of mass-related volume (m3/day)
      real*8 VDYNDV (2, MBI)
C              - Time derivative of dynamic volume      (m3/day)
      real*8 VBUFDV (2, MLC)
C      Dim: MBI: Even derivatives for unused basin indexes must be set,
C                since ACSL integrates all elements.
      real*8 MAX_dT
C           Maximum allowed time_step according to transports
C           which can ensure stability (approximately?).


C  ----------- Allocated arrays for transport description
C              from calling subroutine:

      INTEGER MLV, MLH, MC

C    Water transport arrays, allocated by calling subroutine TRANSP
      real*8 VTRANS        (2,MLV)
C           - Describes vertical transports internal basins

      real*8 VMIX          (2,MLV)
C           - Energy budget for vertical mixing.
C             (1,L) Preliminary values (in TRANC, before TRANV1):
C                       Depth integral of sigma_t from surface to
C                       bottom of layer L.
C                   Final value (TRANV1):
C                       Tidal mixing energy in units kg/m2/s,
C                       surface layer value are subtracted from BSFLUX
C                       to give net change of buoyancy after
C                       vertical mixing has counteracted effect
C                       of inflow of water of lighter density.
C                       Any surplus mixing intensity (BSFLUX<0)
C                       is used outside in main model to reduce
C                       buildtup bouyancy which otherwise would
C                       countereffect wind induced mixing

C             (2,I) Potential energy released/consumed in reduced
C                   gravity field by advective flows finding their
C                   target depth ( see TRANC, TRANV1 )

      real*8 HTRANS        (MLH)
C           - Horizontal transport across connections
C                 Ref. subroutine TRANH and MTRAN1

      INTEGER HTR_L        (MLH,2)
C           - Index of participating layers for flows HTRANS.
C                      Cfr. subroutine TRANH, TRANC, and MTRAN1.
C               Indexes are set local to each basin,
C               with 1 for surface layer.

      INTEGER INDXH        (MC+1)
C           - Buffer segmenting of HTRANS and HTR_L
C             between connections, set by TRANH.

      INTEGER NQ_TR
      real*8 QJET_TRANSP   (MQ_TR)
C           - Vertical water transports due to dived outlets QJET
C                 as m3/s from JETMIX, transformed to m3/day before use
      INTEGER QJET_TRANSP_INDEX (2,MQ_TR)
C           - Global index for source and destination layer
C             for these transports
      INTEGER RCQNDX (dimMS)
C           - Final receiving layer for all land runoffs QJET,
C             both surface outlets (rivers) and dived jets.

      real*8 Ambient_Volume_Flux (dimMS) ! volume mixed into dived outlets
      real*8 Neutral_Depth       (dimMS) ! neutral depth for dived outlets

C  ------------------ work arrays -----------
      INTEGER MW, MT_TEMP, ML_C_MAX, MB
C           - dimension of work arrays

      real*8  WK_ARRAY     (0:MW)
C            - Scratch array for storing flux values locally in TRANH.
C              and TRANV1, and as dummy array in calls to MTRAN1.
C              Length set to 0..N, with N = sum of layers + sum of
C              layer connections, giving more than enough space.

      real*8  TR_TEMP      (MT_TEMP)
      INTEGER TR_SOURCE_L  (MT_TEMP)
      real*8  CONNECT      (ML_C_MAX,2)
C            - Intermediate transport description used in TRANH

      real*8  ZSURF_W      (MB)
C              double precision internal values for surface levels

      real*8  surf_change        (MB)
C              double precision internal values for storing
C              surface level change and dynamic volume change

      real*8  BSFLUX       (MB)
C          - buoyancy flux into surface layer from horisontal transports
C            exported in  unit (kg/s/m2).

      real*8 BWFREQ (NLI) ! Stability, as Brunt_w„isela frequency
      real*8 ZMID   (NLI) ! Mean depth of layers (volume weighted)



!      LOGICAL  ITERATE_SURFACE
!      EXTERNAL ITERATE_SURFACE


C ===================================================================

C                        local variables:

      INTEGER IB, IS, L, IC, ILH, NLH, NLH_DIM, ILC, NL_C, I
      INTEGER IBA, IBB, ILA, ILB, NLA, NLB

      real*8   TRDYNAM, VOL_INCREASE, DENS_D, QJET_TEMP
!      real*8   VLAYER_SUM
		real*8 X
      INTEGER LSURF

      real*8 SEC_PER_DAY
      parameter (SEC_PER_DAY = 24.*3600.)

C Dummy & intermediate variables in calls to TRANH, TRANC and MTRAN1:
      real*8 BSFLUX_EXT
      real*8 ZSURF_EXT, VTRANS_EXT(1), VMIX_EXT(1)
      real*8 XMIX_EXT /0.0/

$if defined DEBUG_TSTEP
      real*8 CURRENT_TIME_STEP
$endif

$if defined DEBUG_ANY
!      INCLUDE 'DEBUG.INC'
      INTEGER K

      
      IF(DEBUG) THEN
           WRITE(DEBUG_UNIT,*) '=========== TRANSP_2.TRANSP_SETUP:'
           WRITE(DEBUG_UNIT,*) 'HTR_CALC = ',HTR_CALC
      ENDIF
$endif
C ==================================================================

      CALL HELLO ('entered TRANSP_2' )


$if defined DEBUG_TSTEP
      CURRENT_TIME_STEP = MAX_DT
$endif

C ----------------------------------------------------------------------
C   Horisontal transports may have been deactivated.  This is signalled
C   from calling program (WATER_TRANSPORTS in TRANSP_1) by HTR_CALC =.F.

C   In that case:
C       Subroutines calculating surface levels and horisontal transports
C       are not called, and terms affecting derivative
C       of total basin volumes are skipped.

C       Input of water from land is neglected in volume balance,
C       so chemical influx is regarded as injection of substances only.

C   The purpose of this possibility is to make it possible to run
C   faster simulations of biochemical processes, with vertical exchange.
C   The transports can be turned on and off throughout the simulation.
C ----------------------------------------------------------------------



C -------------------------------------------------------------------
C A.                Initiate transport terms :
C -------------------------------------------------------------------

C 1. Buoyancy flux terms BSFLUX from horizontal transports,
C    to be set by TRANH.

C 2. Vertical transport arrays VTRANS,
C    used in TRANSP as work area for accumulating temporarily
C    net effect of horisontal and other transports affecting
C    volume continuity.
C       Subroutines TRANC sets the temporary description, which
C    is then used by subroutine TRANV1 to convert VTRANS into
C    the final transport description and to limit time-step.

C 3. Find depth integral of sigma_t, and initiate accumulator for
C    energy input from advective flows to vertical mixing.
C    Later used to calculate mixing energy from advective inflow
C    of denser water from other basins to greater depths,
C    see TRANC and TRANV1.

C 4. VTOTDV array used to accumulate volume input from
C    land and through surface-water exchange.  This will be
C    used locally in this subroutine during surface iteration and
C    in calculating volume conservation in transport terms,
C    but is reset to derivative of effective "chemical volume"
C    of basins before return from the subroutine.

c ===================================================================
   10 CONTINUE ! RESTARTING POINT IF DEBUG IS TURNED ON BY SHORT STEP
c ===================================================================

      DO IB = 1,NBI

         LSURF = INDXI(IB)+1

         BSFLUX(IB) = 0.0

         DO L = LSURF,INDXI(IB+1)
             VTRANS(1,L) = 0.0  ! Accumulates net inflow to layer
             VTRANS(2,L) = 0.0  ! Accumulates gross advective outflow
         END DO

C             Include effect of water-atmosphere exchange:
         VTRANS(1,LSURF) = QWSURF(IB)*sec_per_day ! m3/day
C             changed total volume only if horisontal transports active:
         IF (HTR_CALC ) THEN
             VTOTDV(IB) = QWSURF(IB)    ! m3/s
         ELSE
             VTOTDV(IB) = 0.0
         ENDIF
      END DO



C -------------------------------------------------------------------
C B.           Effect of land runoff and forced circulation
C -------------------------------------------------------------------


$if defined DEBUG_SETUP_GT_1
         IF(DEBUG)
     &       WRITE(DEBUG_UNIT,*) ' ===== LAND RUNOFF ======,   MS= ',MS
$endif


C ------------- Scan through outlets, set up derived transports,
C               and store effect on volume influx terms:
      NQ_TR = 0
      DO IS = 1,NS

         if (QJET(IS,1).le.0.0 .and. QJET(IS,2).le.0.0) THEN
            RCQNDX(IS) = RNFNDX(IS,1) ! Dummy, to avoid problems elsewhere
            Ambient_Volume_Flux (IS) = 0.0
            Neutral_Depth (IS) = 0.0
            CYCLE
         ENDIF

         IB = BASINQ(IS)     ! Basin number

         if ( DT_QTRANSP_CALC .eq. 0 ) THEN
              ! Time to recalculate dived outlets:

            ILA   = RNFNDX(IS,1)          ! Outlet layer (global index)
            LSURF = INDXI(IB) + 1         ! Surface layer
            NLA = INDXI(IB+1) - INDXI(IB) !Number of layers in basin

C     -------- set up derived transports for dived outlets:
            IF ( ILA. gt. LSURF) THEN

$if defined DEBUG_SETUP_GT_1
              IF(DEBUG) THEN
                 WRITE( DEBUG_UNIT,*) ' IS=',IS,
     &                  ' DIVED OUTLET TO BASIN', IB
                 WRITE( DEBUG_UNIT, '(1X,A8,5A12,A8)' )
     &             'RNFNDX', 'DEPTHQ',' ', 'QJET',' ','QDIAM','NHOLES'
                 WRITE( DEBUG_UNIT, '(1X, I8,5E12.5,I8)' )
     &             ILA, (DEPTHQ(IS,K), QJET(IS,K), K=1,2),
     &                  QDIAM(IS), NHOLES(IS)
              ENDIF
              K = NQ_TR  ! only used in printout below
$endif

C        .......Calculate mixing/entrainment in JETMIX.FOR,
C               which will accumulate entrainment transports,
C               and return final receiving layer index in ILB.
              ILB = ILA
!
!                      ! temperature of discharge (heat flux over volume)
!                      ! (QTEMP has dimension temp*volume/time)
!              if (QJet(IS,1).gt.0.0) then
!                 QJET_Temp = QTEMP(IS)/QJet(IS,1)/24./3600.
!                             ! QTEMP as degC*(m3/d), QJET as m3/s
!
!                 WRITE( DEBUG_UNIT, '(1X, A8,3A16/1X,I8,3E16.6)' )
!     &             'IS','QTEMP(IS)',' QJET(IS,1)',' QJET_Temp',
!     &             IS, QTEMP(IS), QJET(IS,1), QJET_Temp
!                 
!              else
!                 QJET_Temp = 0.0
!              endif
!
!             CALL HELLO ('calling JETMIX:' )
! ************ new code: 
              CALL JETMIX( DEBUG, NLA, DEPTH,
     &               SAL(LSURF), TEMP(LSURF), DENSI(LSURF),
     &               INDXI(IB), RNFNDX(IS,1), QJET(IS,1),
     &               QSPP(IS), dDens_dSPP, QJET_Temp,
     &               DEPTHQ(IS,1),
     &               RNFNDX(IS,2),
     &               QJET(IS,2),
     &               QDIAM(IS), NHOLES(IS),
     &             MQ_TR, QJET_TRANSP, QJET_TRANSP_INDEX, NQ_TR, ILB,
     &             Ambient_Volume_Flux(IS), Neutral_Depth(IS) )
! ************ old code: 
!              CALL JETMIX( DEBUG, NLA, DEPTH, DENSI(LSURF), INDXI(IB),
!     &                RNFNDX(IS,1), QJET(IS,1), DEPTHQ(IS,1),
!     &                RNFNDX(IS,2), QJET(IS,2),
!     &                QDIAM(IS), NHOLES(IS),
!     &           MQ_TR, QJET_TRANSP, QJET_TRANSP_INDEX, NQ_TR, ILB,
!     &           Ambient_Volume_Flux(IS), Neutral_Depth(IS) )


             CALL HELLO ('returned from JETMIX:' )


                      ! will now have ILB < =ILA

! $$$$$$$$$$$$
$if defined DEBUG_SETUP_GT_1
              IF(DEBUG) THEN
                 write(DEBUG_UNIT,'('' IN TRANSP_2 ''/1x,A3,2A15)')
     &             'IS', 'Amb_Vol_Flux','Neutral_D'
                 write(DEBUG_UNIT,'(1x,I3,2G15.6)')
     &                 IS, Ambient_Volume_Flux(IS), Neutral_Depth(IS)
                 IF ( K .LT. NQ_TR ) THEN
                     WRITE( DEBUG_UNIT, '(1X, A8, A18,2A12)' )
     &                  'I:','QJET_TRANSP (m3/s)',
     &                  ' FROM LAYER',' TO LAYER'
                     WRITE( DEBUG_UNIT, '(1X, I8, E21.15, 2I12)' )
     &                  ( I, QJET_TRANSP(I),
     &                       QJET_TRANSP_INDEX(1,I),
     &                       QJET_TRANSP_INDEX(2,I),
     &                    I=K+1,NQ_TR )
                 ELSE
                     WRITE( DEBUG_UNIT, * )
     &                   'NO TRANSPORTS SET UP BETWEEN LAYERS'
                 ENDIF
                 WRITE(DEBUG_UNIT,*)
     &             ' ILA= ', ILA, ' ILB= ', ILB
              ENDIF
$endif
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

            ELSE

C -------- Surface outlet: find lowest layer to receive part of outlet:
               L = 2
               DO WHILE  ( QDIAM(IS) .gt. DEPTH(L) )
                  L = L + 1
                  if ( L .gt. NLA ) EXIT
               ENDDO
               ILB = L+LSURF-2 ! Global index

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4
$if defined DEBUG_SETUP_GT_1
               IF (DEBUG) THEN
                  WRITE ( DEBUG_UNIT,'(1x,4(a,I3))')
     &             ' Outlet nr.',IS,
     &             ' to surface, distributed between layers ', ILA,
     &             ' and ', ILB, ' in basin ',IB
               ENDIF
$endif
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4


               Ambient_Volume_Flux (IS) = 0.0
               Neutral_Depth (IS) = 0.0
                 ! (Only exported for information, not used otherwise.
                 ! (Default assumption of 0, will be reset
                 ! (by JETMIX for dived outlets



            ENDIF

            RCQNDX(IS) = ILB
C      For submerged outlet: Destination layer (entrainment) <=ILA
C      For surface outlet: Lowest layer receiving outlet     >=ILA

         ELSE

C ------------- within time interval from last evaluation:
C               Use previously calculated layer numbers:
            ILA = RNFNDX(IS,1)
            ILB = RCQNDX(IS)

         ENDIF


C     ----- Include effects of the outlet flux itself:
         ILA = MIN(ILB,ILA) ! now = top layer to receive influx.
         X = QJET(IS,1)*sec_per_day / FLOAT(ILB-ILA+1)
C                       Volume influx to each layer:
         DO L = ILB, ILA, -1
            VTRANS(1,L) = VTRANS(1,L) + X

$if defined DEBUG_SETUP_GT_1
             IF(DEBUG) WRITE( DEBUG_UNIT, '(1x,A,I5,2(1x,A,G21.15))' )
     &          'VTRANS(1,' , L, ') +' , X, ' =:', VTRANS(1,L)
$endif

         ENDDO

C           Changes VTOTDV only if horisontal transports are active:
         IF (HTR_CALC) THEN
             VTOTDV(IB)  = VTOTDV(IB) + QJET(IS,1)  ! m3/s

$if defined DEBUG_SETUP_GT_1
             IF(DEBUG) WRITE( DEBUG_UNIT, '(1x,A,I5,A,'':'',G21.15)' )
     &               'VTOTDV(',IB,')=',VTOTDV(IB)
$endif

         ENDIF

      ENDDO

C     ----- After all outlets are processed,
C           include effect of derived internal transports
C           on net influx terms for each layer:
      DO IS = 1, NQ_TR
         QJET_TRANSP(IS) = QJET_TRANSP(IS)*SEC_PER_DAY ! m3/s -> m3/day
         ILA = QJET_TRANSP_INDEX(1,IS)
         ILB = QJET_TRANSP_INDEX(2,IS)
C           net influx terms:
         VTRANS(1,ILA) = VTRANS(1,ILA) - QJET_TRANSP(IS)
         VTRANS(1,ILB) = VTRANS(1,ILB) + QJET_TRANSP(IS)
C           gross outflux term:
         VTRANS(2,ILA) = VTRANS(2,ILA) + QJET_TRANSP(IS)

$if defined DEBUG_SETUP_GT_1
         IF(DEBUG) WRITE( DEBUG_UNIT, '(1X, A, I5,A,G21.15)' )
     &           'QJET_TRANSP(',IS,') (m3/day)=',QJET_TRANSP(IS)
$endif
      ENDDO


$if defined DEBUG_SETUP_GT_1
          IF(DEBUG) THEN
      WRITE(DEBUG_UNIT,'(A)')  '=========== VTRANS, after JETMIX:',
     &   '(1,K) now contains net influx to layer K',
     &   ' due to local inputs from land, including dived jet mixing'
      WRITE(DEBUG_UNIT,'(1X,A5,2A15)' ) 'K','VTRANS(1,K)','VTRANS(2,K)'
      WRITE(DEBUG_UNIT,'(1X,I5,2G21.15)')(K,(VTRANS(I,K),I=1,2),K=1,MLV)
          ENDIF
$endif


C ----------------------------------------------------------------------
C                       Surface iteration
C ----------------------------------------------------------------------

C   Loads surface levels values calculated from integrated dynamic
C   volumes VDYN(1,..) into work_array ZSURF_W,
C   as starting point for iteration:
C


         DO IB = 1,NBI
            LSURF = INDXI(IB)+1
            ZSURF_W(IB) = (VDYN(1,IB))/AREA(LSURF)

$if defined DEBUG_SETUP
         IF(DEBUG) THEN
           WRITE(DEBUG_UNIT,*)' --- VDYN/AREAS:'
           WRITE(DEBUG_UNIT, '(1x,2(A5,I6),3(A12,G15.7))' )
     &        ' IB:',IB,' LSURF:',LSURF, ' VDYN(1,IB):',VDYN(1,IB),
     &        ' AREA(LSURF):',AREA(LSURF),' ZSURF_W(IB):',ZSURF_W(IB) 
         ENDIF
$endif
         END DO ! Surface levels determined by integrating transports




C ------ If horisontal transports is active: ----------------------
C     Calculate horizontal transports in UFLOW (m/s) and WK_ARRAY (m2/s)
C     If activated by logical ITRZ: Iterates to equilibrium surface
C     levels, which by definition will give approximately constant
C     dynamic volume rate of change in time.
C     Otherwise, just uses defined ZSURF_W values.

      IF ( HTR_CALC ) THEN

$if defined DEBUG_SETUP
         IF(DEBUG) THEN
            NLB = 0  ! ACCUMULATES MAX. NUMBER OF LAYERS
            WRITE( DEBUG_UNIT,'(4(A,I5:))')
     &         ' **** Calls TRHITR_PREP;  NBI=', NBI, ' NLI=', NLI,
     &         ' NBE=', NBE, ' NLE=', NLE, 'NC=', NC, 'NLC=',NLC
            DO IB = 1,NBI
                NLB = MAX( NLB, INDXI (IB+1) - INDXI(IB) )
                WRITE( DEBUG_UNIT, '( 3(1X,A,I5),1X,A,g13.7)' )
     &               'INTERNAL BASIN ', IB, 'LAYERS ',INDXI(IB)+1,
     &               ' TO ', INDXI (IB+1), 'VTOTDV:',VTOTDV(IB)
                WRITE( DEBUG_UNIT, '(3('' DENSI('',I3,'')=''g13.7) )')
     &              ( L, DENSI(L), L = INDXI(IB)+1,INDXI(IB+1) )
            ENDDO
            DO IB = 1,NBE
                NLB = MAX( NLB, INDXE (IB+1) - INDXE(IB) )
                WRITE( DEBUG_UNIT, '(3(1X,A,I5),1X,A,g13.7)' )
     &               'EXTERNAL BASIN ', IB, 'LAYERS ',INDXE(IB)+1,
     &               ' TO ', INDXE (IB+1),'ZSURFE=', ZSURFE(IB)
                WRITE( DEBUG_UNIT, '(3('' DENSI('',I3,'')=''g13.7) )')
     &              ( L, DENSE(L), L = INDXE(IB)+1,INDXE(IB+1) )
            ENDDO
            DO IC = 1,NC
                NLB = MAX( NLB, INDXC (IC+1) - INDXC(IC) )
                WRITE( DEBUG_UNIT, '(3(1X,A,I5),1X,A,g13.7 )')
     &               'CONNECTION ',IC, 'BETWEEN ', BConn1(IC),
     &               'AND ', BConn2(IC), 'DPEFF=', DPEFF(IC)
                WRITE( DEBUG_UNIT, '(3(:'' WIDTH('',I3,'')=''g13.7) )')
     &              ( L, WIDTH(L), L = INDXC(IC)+1,INDXC(IC+1) )
            ENDDO

            WRITE( DEBUG_UNIT, '(3(:'' DEPTH('',I3,'')='',g13.7) )')
     &              ( L, DEPTH(L), L = 1, NLB+1 )
         endif
$endif

         CALL HELLO ('calling TRHITR_PREP' )

C  Prepare iteration step:
         CALL  TRHITR_PREP (  DEPTH, NBI, VTOTDV, XMIX,
     &          INDXI, NLI, DENSI,
     &          NBE, INDXE, NLE, DENSE, ZSURFE, DZDTX,
     &          NC, INDXC, ZSILL, NLC, BConn1, BConn2, WIDTH, DPEFF )
         CALL HELLO ('RETURNED FROM TRHITR_PREP' )

C Perform iteration to equilibrium surface levels for given ZSURFE,
C returns adjusted Z_SURFW in m, corresponding surf_change in unit m/s.

         IF ( .NOT. ITERATE_SURFACE ( DEBUG, ITRZ, MAX_DT,
     &                  NBI, ZSURF_W, INDXI, NLI, AREA, surf_change,
     &                  NLC, UFLOW )
     &      ) THEN

$if defined DEBUG_SETUP
                 IF(DEBUG) THEN
                     WRITE(DEBUG_UNIT,*)
     &                  '     Surface iteration did not converge'
!                     PAUSE 'Press Enter to continue'
                 ENDIF
$endif
         ENDIF
         CALL HELLO ('RETURNED FROM ITERATE_SURFACE' )

$if defined DEBUG_TSTEP
        IF (DEBUG)  CALL TSTEP_DEBUG_PRINT
     &                  ( 'TRANSP_2 (after Iterate_surface)',
     &                    CURRENT_TIME_STEP, MAX_DT )
$endif


$if defined DEBUG_SETUP
         IF(DEBUG) THEN
           WRITE(DEBUG_UNIT,  '(1x,A,L6,A)' )
     &        ' --- After surface iteration, ITRZ=', ITRZ, ' ZSURF_W:'
           WRITE(DEBUG_UNIT,'(1x,4E21.15)')  ZSURF_W
         ENDIF
$endif


      ENDIF
C        Otherwise: ZSURF_W contains integrated values VDYN/Area


C     ------- surface volume change rates in m3/day:
      DO IB = 1,NBI
         LSURF = INDXI(IB)+1

$if defined DEBUG_VMIX
      if (DEBUG)
     &     WRITE(DEBUG_UNIT,
     &           '('' Density profile basin '',I3/1x,A5,3A18)')
     &         IB, ' L ', ' DENSI(L) ',' X ',' DENS_D '
$endif
         DENS_D = ZSURF_W(IB)*DENSI(LSURF)
         DO L = LSURF,INDXI(IB+1)
             I = L-LSURF+1
             X = DEPTH(I+1) - DEPTH(I)
             DENS_D = DENS_D + DENSI(L)* X
             VMIX(1,L) = DENS_D
                ! Depth integrated density at bottom of layer
             VMIX(2,L) = 0.0     ! Accumulator of mixing energy

$if defined DEBUG_VMIX
      if (DEBUG)
     &      WRITE(DEBUG_UNIT,'(1x,I5,3G21.15)')
     &           L, DENSI(L), X, DENS_D
$endif

         ENDDO

         IF ( HTR_CALC ) THEN
C            ------- Store dynamic volume change rate from iteration:
             VDYNDV(1,IB) = surf_change(IB) * SEC_PER_DAY *AREA(LSURF)
C            ------- Initiate volume change rate with local inputs,
C                    to be completed below with horisontal transports:
             VDYNDV(2,IB) = VTOTDV(IB)*sec_per_day
         ELSE
             VDYNDV(1,IB) = 0.0
             VDYNDV(2,IB) = 0.0
         ENDIF
      END DO



C ---------- Calculate transports across each connection: --------------
C            at the same time adding buoyancy fluxes
C            into top and bottom of mixed surface volumes.


C #############################
      IF (.NOT.HTR_CALC ) GO TO 100 ! Horisontal transports deactivated.
C #############################

      INDXH(1)=0 ! Transport buffer HTRANS used dynamically.
      DO IC = 1,NC
         ILC = INDXC(IC)+1
         NL_C = INDXC(IC+1)-INDXC(IC)
         IBA = BCONN1(IC)
         ILA = INDXI(IBA)+1
         NLA = INDXI(IBA+1)-INDXI(IBA)
         IBB = BCONN2(IC)
         IF(IBB.GT.0) THEN
C   ........... Between internal basins:
            ILB = INDXI(IBB)+1
            NLB = INDXI(IBB+1)-INDXI(IBB)

            ILH = INDXH(IC)+1

$if defined DEBUG_TRANH_TRANC_CALL
      if (DEBUG)
     &WRITE( DEBUG_UNIT,
     &       '('' --------- calls TRANH: '',A/1x,10A6/1x,10I6)')
     &   ' 2. basin internal',
     &   'IC','ILC','NL_C','IBA','ILA','NLA','IBB','ILB','NLB','ILH',
     &    IC,  ILC,  NL_C,  IBA,  ILA,  NLA,  IBB,  ILB,  NLB,  ILH
$endif
         CALL HELLO ('Calls TRANH' )
            CALL TRANH ( DEBUG, DEPTH, NL_C, WIDTH(ILC), HTRMIX(IC),
     &         ZSURF_W(IBA), NLA, DENSI(ILA), BSFLUX(IBA), XMIX(IBA),
     &         ZSURF_W(IBB), NLB, DENSI(ILB), BSFLUX(IBB), XMIX(IBB),
     &         MLH, ILH, INDXH(IC+1), HTRANS, HTR_L, TRDYNAM,
     &         VBUF(1,ILC), VBUFDV(1,ILC),
     &         VBUFMX(1,ILC), VBUFTR(1,ILC), TCVBUF(1,IC),
     &         MW, WK_ARRAY, ZSILL(IC), UFLOW(ILC),
     &         MT_TEMP, TR_TEMP, TR_SOURCE_L, ML_C_MAX, CONNECT )
C      ..... Update derivative of dynamic volume (m3/day):
            VDYNDV(2,IBA) = VDYNDV(2,IBA) - TRDYNAM
            VDYNDV(2,IBB) = VDYNDV(2,IBB) + TRDYNAM

C      ..... Update VTRANS(K,I..) with horizontal transport terms:
C                          K=1: net influx to basin layers
C                          K=2: sum of outflows from basin layers
C            Update VMIX(2,I) with net potential energy released in
C            reduced gravity field by advective flows to and from layer
            NLH = INDXH(IC+1)-INDXH(IC)
            NLH_DIM =MAX(1,NLH)

$if defined DEBUG_TRANH_TRANC_CALL
      if (debug)
     &WRITE( DEBUG_UNIT,
     &       '('' --------- calls TRANC: ''/1x,9A6/1x,9I6)')
     &   'ILC','NL_C','ML_C_MAX','NLH','ILH','NLA','ILA','NLB','ILB',
     &    ILC,  NL_C,  ML_C_MAX , NLH , ILH , NLA,  ILA , NLB , ILB 
$endif
         CALL HELLO ('  .......... Calls TRANC' )
            CALL TRANC( DEBUG, DEPTH, UFLOW (ILC),
     &        NL_C, ML_C_MAX, CONNECT,
     &        NLH_DIM, NLH, HTRANS(ILH), HTR_L(ILH,1), HTR_L(ILH,2),
     &        NLA, DENSI(ILA), VTRANS(1,ILA), VMIX(1,ILA),
     &        NLB, DENSI(ILB), VTRANS(1,ILB), VMIX(1,ILB)  )


         ELSE
C   ........... Between internal (A) and external (B) basin:
            IBB = ABS(IBB)
            ILB = INDXE(IBB)+1
            NLB = INDXE(IBB+1)-INDXE(IBB)
            ZSURF_EXT = ZSURFE(IBB)
            ILH = INDXH(IC)+1
$if defined DEBUG_TRANH_TRANC_CALL
      if (debug)
     &WRITE( DEBUG_UNIT,
     &       '('' --------- calls TRANH: '',A/1x,10A6/1x,10I6)')
     &   ' 2. basin external',
     &   'IC','ILC','NL_C','IBA','ILA','NLA','IBB','ILB','NLB','ILH',
     &    IC,  ILC,  NL_C,  IBA,  ILA,  NLA,  IBB,  ILB,  NLB,  ILH
$endif

         CALL HELLO ('  .......... Calls TRANC' )

            CALL TRANH ( DEBUG, DEPTH, NL_C, WIDTH(ILC), HTRMIX(IC),
     &         ZSURF_W(IBA), NLA, DENSI(ILA), BSFLUX(IBA), XMIX(IBA),
     &         ZSURF_EXT, NLB, DENSE(ILB), BSFLUX_EXT, XMIX_EXT,
     &         MLH, ILH, INDXH(IC+1), HTRANS, HTR_L, TRDYNAM,
     &         VBUF(1,ILC), VBUFDV(1,ILC),
     &         VBUFMX(1,ILC), VBUFTR(1,ILC),  TCVBUF(1,IC),
     &         MW, WK_ARRAY, ZSILL(IC), UFLOW(ILC),
     &         MT_TEMP, TR_TEMP, TR_SOURCE_L, ML_C_MAX, CONNECT )
C      ..... Update derivative of dynamic volume in inner basin:
            VDYNDV(2,IBA) = VDYNDV(2,IBA) - TRDYNAM
C      ..... Update VTRANS and VMIX as above, but only for inner basin A
            NLH = INDXH(IC+1)-INDXH(IC)
            NLH_DIM = MAX(1,NLH)
            NLB = 0  ! Prohibits update of VTRANS_EXT [VTRANS(K,ILB)]

$if defined DEBUG_TRANH_TRANC_CALL
      if (debug)
     &  WRITE( DEBUG_UNIT,
     &       '('' --------- calls TRANC: ''/1x,9A6/1x,9I6)')
     &   'ILC','NL_C','ML_C_MAX','NLH','ILH','NLA','ILA','NLB','ILB',
     &    ILC,  NL_C,  ML_C_MAX , NLH , ILH , NLA,  ILA , NLB , ILB 
$endif
         CALL HELLO ('  .......... Calls TRANC' )

            CALL TRANC(  DEBUG, DEPTH, UFLOW (ILC),
     &          NL_C, ML_C_MAX, CONNECT,
     &          NLH_DIM, NLH, HTRANS(ILH), HTR_L(ILH,1), HTR_L(ILH,2),
     &          NLA, DENSI(ILA), VTRANS(1,ILA), VMIX(1,ILA),
     &          NLB, DENSE(ILB), VTRANS_EXT   ,  VMIX_EXT     )
         ENDIF

$if defined DEBUG_SETUP_GT_1
        if (debug) then
         WRITE(DEBUG_UNIT,'('' connection nr. '',I5)' )  IC
         WRITE(DEBUG_UNIT,'('' INDXH(IC)='',I6,'' INDXH(IC+1)='',I6)')
     &       INDXH(IC), INDXH(IC+1)
         IF (NLH.GT.0) THEN
            WRITE(DEBUG_UNIT,'('' connection nr. '',I5/1X,A5,3A15)' )
     &          IC, 'K','HTR_L(K,1)','HTR_L(K,2)','HTRANS(K)'
            WRITE(DEBUG_UNIT,'(1X,I5,2I15,G21.15)' )
     &       (K,(HTR_L(K,I),I=1,2),HTRANS(K), K=INDXH(IC)+1,INDXH(IC+1))
         ENDIF
        endif
$endif
      END DO

$if defined DEBUG_SETUP_GT_1
          IF(DEBUG) THEN
      WRITE(DEBUG_UNIT,*)' ======== VTRANS/VMIX etter TRANH1/TRANC:'
      WRITE(DEBUG_UNIT,'(1X,A5,A12,3A18)' )
     &     'K','VTRANS(1,K)','VTRANS(2,K)','VMIX(1,K)','VMIX(2,K)'
      WRITE(DEBUG_UNIT,'(1X,I5,4E21.15)')
     &         ( K, (VTRANS(I,K),I=1,2),(VMIX(I,K),I=1,2),K=1,MLV)
          ENDIF
$endif


C ############################
  100 CONTINUE
C ############################


C ----------------- Ensure volume continuity: -----------------------

C   VTRANS(1,L) contains net inflow(+) or outflow(-) for
C   the layers in all internal basins due to horizontal transports
C   and land runoff, including mixing of dived jets.
C   Below, VTRANS(1,L) is reset to include advective
C   transports between neighboring layers within basins.

      DO IB = 1,NBI

$if defined DEBUG_SETUP_GT_1
      IF (DEBUG) WRITE (DEBUG_UNIT,'('' IB='',I4)') IB
$endif
         LSURF  = INDXI(IB)+1


C   ........ Accumulate total net inflow to basin from bottom up:
         VOL_INCREASE  = 0.0
!        VLAYER_SUM = 0.0
         DO L = INDXI(IB+1), LSURF,-1
!           VLAYER_SUM = VLAYER_SUM + VLAYER(L)
            VOL_INCREASE = VOL_INCREASE + VTRANS(1,L)
            VTRANS(1,L) = VOL_INCREASE
C             = accumulated inflow = flow up through fixed boundary
$if defined DEBUG_SETUP_GT_1
      IF (DEBUG)
     &   WRITE (DEBUG_UNIT,'('' L='',I4,'' SUM(VTRANS)='',G20.11)')
     &          L, VOL_INCREASE

$endif
         END DO

         IF(.NOT. HTR_CALC) THEN
C        ..... Store a virtual redrawal of water from top layer,
C                 used in TR2_CONC_DERIV, see below:
             VTRANS(1,LSURF) = VOL_INCREASE

$if defined VOLUME_VERSION_2
         ELSE
C        ..... Store total (inflow-outflow) as time derivative of
C             "chemical" volume:
            VTOTDV(IB) = VOL_INCREASE

C        ..... Change to relative accumulated inflow and calculate
C              corrected vertical flows through layer limits,
C              corresponding to layer limits being adjusted
C              so that each layer remains a constant fraction
C              of the total (changing) volume:
            VOL_INCREASE = VOL_INCREASE/VFROPN(IB)
            VLAYER_SUM = 0.0
            DO L = INDXI(IB)+NLVOPN(IB), LSURF, -1
                VLAYER_SUM = VLAYER_SUM + VFRAC(L)
                VTRANS(1,L) = VTRANS(1,L) - VOL_INCREASE*VLAYER_SUM
C                          >0 for flow up from L to   L-1
C                          <0 for flow down to L from L-1
            END DO  ! VTRANS(1,LSURF) should now be =0
c            write(*,*)'IB:',IB, 'VTRANS(1,LSURF)', VTRANS(1,LSURF),
c     &                'vlayer_sum', vlayer_sum,
c     &                'vfropn(ib)', vfropn(ib)
$endif
         ENDIF

C   ........ Buoyancy influx term converted to area-specific value:
         BSFLUX(IB) =  BSFLUX(IB)/ AREA(LSURF)
         ! kg/m2/s  =    (kg/s)  / m2

C
$if defined DEBUG_SETUP_GT_1
            IF (DEBUG) THEN
      write (DEBUG_UNIT,*)'=======BSFLUX(IB)=', BSFLUX(IB)
      write (DEBUG_UNIT,*)'=======VTRANS after continuity corrections'
      WRITE(DEBUG_UNIT,'(1X,A5,2A15)' ) 'K','VTRANS(1,K)','VTRANS(2,K)'
      WRITE(DEBUG_UNIT,'(1X,I5,2E21.15)')
     &            (K,(VTRANS(I,K),I=1,2),K=INDXI(IB)+1,INDXI(IB+1))
            ENDIF
      IF (DEBUG) WRITE (DEBUG_UNIT, '('' VTOTDV='',G21.15)') VTOTDV(IB)
$endif

      END DO


C -------- Update VTRANS with diffusive vertical transports
C          inside all basins, and adjust MAX_DT initiated outside
      CALL TRANV1(DEBUG, DEPTH, MBI, NBI,
     &   MIXFAC, MIXEXP, N2LIM,
     &   SFMIXC, SFMIXZ, GMIXFR, GMIXDC, GMIXDX,
     &   ZBOTMI, INDXI, NLI, AREA, VLAYER, VLCORF, NLVOPN,
     &   DENSI, VTRANS, VMIX,
     &   MW, WK_ARRAY, BWFREQ, ZMID,
     &   MAX_DT )


$if defined DEBUG_TSTEP
      if (debug)
     &   CALL TSTEP_DEBUG_PRINT('TRANSP_2 (after TRANV1)',
     &        CURRENT_TIME_STEP, MAX_DT)
$endif


C   ........ Update buoyancy influx term with work against gravity
C            performed in surface layer and stored by TRANV1:
      DO IB = 1,NBI
         LSURF = INDXI(IB)+1
         BSFLUX(IB) =  BSFLUX(IB) - VMIX(1,LSURF)   ! kg/m2/s

$if defined DEBUG_SETUP
      if (debug)
     &    write (DEBUG_UNIT,'('' Basseng '', I4, '':'',2(1x,A,G12.5))')
     &           IB, 'VMIX(1,LSURF)=',VMIX(1,LSURF),
     &           'BSFLUX(IB)', BSFLUX(IB)
$endif


      END DO



C -------- zero unused derivatives, necessary because
C          ACSL integrates all elements in arrays.

      DO IB = NBI+1, MBI
         VTOTDV (IB)  = 0.0
         VDYNDV(1,IB) = 0.0
         VDYNDV(2,IB) = 0.0
      END DO

      IF (HTR_CALC) THEN
         ILC = NLC+1
      ELSE
         ILC = 1 ! no horis. transp., zero all buffer volume derivatives
      ENDIF
      DO L = ILC, MLC
         VBUFDV(1,L) = 0.0
         VBUFDV(2,L) = 0.0
      END DO

$if defined DEBUG_SETUP
      IF(DEBUG) THEN
         write(DEBUG_UNIT,*) ' ====== Tilslutt i kall til TRANSP ======'
         WRITE(DEBUG_UNIT,'(1X,A5,A12,3A18)' )
     &     'K','VTRANS(1,K)','VTRANS(2,K)','VMIX(1,K)','VMIX(2,K)'
         WRITE(DEBUG_UNIT,'(1X,I5,4E21.15)')
     &         ( K, (VTRANS(I,K),I=1,2),(VMIX(I,K),I=1,2),K=1,MLV)
         IF(HTR_CALC) THEN
            WRITE(DEBUG_UNIT, '('' NC='',I6,'' connections'')') NC
            WRITE(DEBUG_UNIT, '('' NC+1 indices INDXH:'')')
            WRITE(DEBUG_UNIT, '(1x,10I7)') (INDXH(IC),IC=1,NC+1)
            IF (NC.gt.0 .and. INDXH(NC+1) .GT. 0) THEN
               WRITE(DEBUG_UNIT, '(/1X,A5,3A15)' )
     &          'K','HTR_L(K,1)','HTR_L(K,2)','HTRANS(K)'
               WRITE(DEBUG_UNIT,'(1X,I5,2I15,G21.15)' )
     &          (K,(HTR_L(K,I),I=1,2),HTRANS(K), K=1,INDXH(NC+1) )
            ELSE
               WRITE(DEBUG_UNIT,*) 'No horizontal transports HTRANS'
            ENDIF

            DO IC = 1, NC
               WRITE(DEBUG_UNIT,'(//1X,A3,4A18)' )
     &                 'ILC', 'VBUF(1,ILC)'  , 'VBUFDV(1,ILC)'  ,
     &                        'VBUF(2,ILC)',   'VBUFDV(2,ILC)'
               DO ILC = INDXC(IC)+1,INDXC(IC+1)
                  WRITE(DEBUG_UNIT,'(1X,I3,4G21.15)' )
     &                 ILC, (VBUF(I,ILC),VBUFDV(I,ILC), I=1,2 )
               ENDDO
            ENDDO

         ENDIF

         WRITE (DEBUG_UNIT, '(/1X,A5,2A18)') 'IB','VTOTDV'
         WRITE (DEBUG_UNIT, '(1X,I5,G21.15)')
     &         (IB,VTOTDV(IB),IB = 1,NBI)
         WRITE (DEBUG_UNIT, '(1X,A5,2A18)') 'IB','VDYN', 'VDYNDV'
         DO IB = 1, NBI
            WRITE (DEBUG_UNIT, '(1X,I5,4G21.15)')
     &         IB, (VDYN(I,IB), VDYNDV(I,IB),I = 1,2)
         ENDDO
         if (ITRZ) THEN
            DO IB = 1,NBI
               lsurf = INDXI(IB)+1
C     Deviation between two values of VDYNDV:
               WRITE( DEBUG_UNIT,*)
     &          'IB: VDYNDV(1,IB) - VDYNDV(2,IB)= difference (m3/day)'
               WRITE( DEBUG_UNIT,'(1X,I3,3E21.15)' )
     &           IB, VDYNDV(1,IB) , VDYNDV(2,IB),
     &               VDYNDV(1,IB) - VDYNDV(2,IB)
            END DO
         ENDIF
      ENDIF
$endif


      END Subroutine



$if defined DEBUG_TSTEP
      SUBROUTINE TSTEP_DEBUG_PRINT(IDENT, CURRENT_TIME_STEP, MAX_DT)
      
      
      CHARACTER*(*) IDENT
      real*8 CURRENT_TIME_STEP, MAX_DT

!      INCLUDE 'DEBUG.INC'

      if (CURRENT_TIME_STEP.ne.MAX_DT) THEN
         WRITE (DEBUG_UNIT,
     &       '(1x,A,'': MAX_DT from '',G15.7,'' to '',G15.7)')
     &       IDENT, CURRENT_TIME_STEP, MAX_DT
         CURRENT_TIME_STEP = MAX_DT
      ENDIF
      END Subroutine
$endif


      End Module