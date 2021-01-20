C Eutrofimodell Indre Oslofjord -  Fil SURF_MIX.FOR            %DATE

      Module fx_Surf_Mix
      use fx_RunControl
      use fx_Sigmat
      
      implicit none


C ...... Variables used in SURFMX and GET_DENS_MIXED
C        for accumulating and calculating properties of mixed layers:

      real*8, private:: VLAYER_MIXED, SAL_MASS_MIXED, TEMP_MASS_MIXED
      real*8, private:: SAL_MIXED(1), TEMP_MIXED(1), DENS_MIXED(1)

      contains
      


C        BFX-values:
$define DEBUG_SURFMX

$define DEBUG_SURFBF
C        ALSO intermediate results

$if defined  DEBUG_SURFMX
$define DEBUG_SURFMX_GT_1
$endif

$if defined DEBUG_SURFBF
$define DEBUG_SURFBF_GT_1
$endif

C ====================================================================
      SUBROUTINE SURFBF  ( MXTEST, MB, NB, QWSURF, BSFLUX,
     &                     NS, QWATER, BASINQ, RNFNDX,
     &                     QOUT, QRAD, INDXI,
     &                     ML, AREA, SALIN, TEMP,
     &               BFX )

      
C Set BFX = relative buoyancy influx pr. horizontal area,
C with effect of freshwater inflow, inflow from other basins,
C precipitation, evaporation and radiative/conductive heat exchange.

C BFX will later be updated by subroutine TRANSP_SETUP
C in Module TRANSP_2, called through transfer routine TRANSP in
C ACSL model, with effect of influx of lighter water from other basins.

C Integral of BFX is later used in SURFMX,
C called from ACSL discrete section ADJUST,
C to calculate depth of well-mixed surface layer.

C BFX has unit m2/s3:
C       energy pr. density and area kgm/s2/(kg/m3)

C INPUT:
      LOGICAL MXTEST      ! Controls debug printing
      INTEGER MB, NB
      real*8 QWSURF   (MB)  ! Precipitation - evaporation   m3/s
      real*8 BSFLUX   (MB)  ! Buoyancy influx from other basins (kg/m2/s)

      INTEGER NS
      real*8    QWATER (NS)   ! Water flux from land runoff   m3/s
      INTEGER BASINQ (NS)   ! Index of receiving basin
      INTEGER RNFNDX (NS)   ! Global index of receiving layer (1..NLI)

      real*8 QOUT     (MB)       ! net heat loss from water to air (W/m2)
      real*8 QRAD                ! absorbed IR radiation at surface (W/m2)

      INTEGER INDXI (NB+1), ML
      real*8 AREA     (ML)       ! m2
      real*8 SALIN    (ML)       ! o/oo
      real*8 TEMP     (ML)       ! oC

C OUTPUT:
      real*8 BFX (MB)         ! m2/s3

C ======================== LOCAL VARIABLES ==========================

      INTEGER IS, IB, L
      real*8 DDENS_DS, DDENS_DT, B

      real*8 Heat_capacity / 4.0e6 /  ! Joule/m3/K
      real*8 GRAV          / 9.81  /  ! Gravitation constant (m/s2)


$if defined DEBUG_SURFBF
!      INCLUDE 'DEBUG.INC'
      if (mxtest)
     &WRITE(DEBUG_UNIT,*) ' =============== SURF_MIX.SURFBF'
$else
      LOGICAL DUMMY
      DUMMY = MXTEST !TO AVOID WARNING
$endif

C ------------ Zero all values (all are integrated in ACSL)
      DO IB = 1, MB
           BFX (IB) = 0.0
      END DO

C ------------ Accumulate effect of freshwater land runoff:
      DO IS = 1,NS
         IB = BASINQ(IS)
         IF ( RNFNDX(IS) .eq. INDXI(IB)+1 ) THEN
             BFX(IB) = BFX(IB) + QWATER(IS)      ! m3/s
         ENDIF
      ENDDO

C ------------ Effect of surface exchange for each basin:

      DO IB = 1, NB
         L = INDXI(IB)+1

$if defined DEBUG_SURFBF_GT_1
      if (mxtest)
     &   WRITE( DEBUG_UNIT,
     &          '('' >>> Basin'',I5,'': Runoff [m3/s]='',G15.7)')
     &             IB, BFX(IB)
$endif

C     ..... Effect of exchange with atmosphere, water and heat:
C                 Using density derivatives on salinity and temp.:
         CALL DDENS ( DDENS_DS, DDENS_DT, SALIN(L), TEMP(L) )

C           Net freshwater influx, including land runoff,
C           will increase buoyancy (ddens_ds >0):
         B = (BFX(IB)+ QWSURF(IB))/AREA(L) * DDENS_DS    * SALIN(L)
C                    (m3/s)       /m2      * (kg/m3/ppt) * ppt
C                                         [DENS=sigma_t]

C           Heat loss, i.e. cooling, will decrease byoyancy
C           if dDensity/dT <0:
         B = B + ( QOUT(IB)-QRAD ) / HEAT_CAPACITY * DDENS_DT
C                (m3/m2s*J/m3 + J/m2s ) /  (J/m3/degC)  * kg/m3/degC

C        Unit of B is now kg/m2/s: [ mass_deficit/area/time]

$if defined DEBUG_SURFBF_GT_1
          IF (MXTEST) THEN
       WRITE(DEBUG_UNIT,'(1x,A/A)')
     &  'B[kg/m2/s]= (BFX(IB)+QWSURF(IB))/AREA(L) * DDENS_DS *SALIN(L)',
     &  '            - ( QOUT(IB)-QRAD ) / HEAT_CAPACITY * DDENS_DT'
       WRITE(DEBUG_UNIT,'(1X,3(G15.7:A:))')
     &    B,'= (', BFX(IB),' + ', QWSURF(IB), ') /',
     &         AREA(L),'*',DDENS_DS,'*',SALIN(L),
     &        '- (',QOUT(IB),'-',QRAD,' ) /',
     &          HEAT_CAPACITY, '*', DDENS_DT
          ENDIF
$endif

C        ..... Add density effect of influx from other basins,
C              rescale to buoyancy in density relative units,
C              and store as value of BFX:
         BFX(IB) = ( B + BSFLUX (IB) ) * GRAV / 1000.0
C          m2/s3 =       kg/m2/s       * m/s2 / (kg/m3)
C              Unit of BFX is now effect/density/mixing depth.

$if defined DEBUG_SURFBF
       if (MXTEST)
     & WRITE(DEBUG_UNIT,'(1x,A/3(G15.7:A:))')
     &   'BFX(IB)[m2/s3]= ( B + BSFLUX (IB) ) * GRAV / 1000.0  :',
     &   BFX(IB),'= (', B,' + ',BSFLUX(IB),') * ',GRAV,'/ 1000.0'
$endif

      END DO


      END Subroutine




C ====================================================================

      SUBROUTINE SURFMX ( MXTEST, FR3INT, DT, ND, DEPTH, NBI,
     &                    INDXI, BFXINT, MLI, VLAYER, DENS, SAL, TEMP,
     &              XMIX, LMIX )


C Calculate depth of well_mixed layer at current time,
C based on energy balance and density profile given by
C Euler integration of previous time step.

C --------------------------------------------------------------
C Called from subroutine CNCADJ included in ACSL model.
C The resulting mixing depth is used to adjust concentrations
C in subroutine CCONC called later in the same section.
C Salt and temperature are adjusted here.
C --------------------------------------------------------------
C Input:
      LOGICAL MXTEST      ! Controls debug printing
      integer NBI
      real*8    FR3INT (NBI)  ! Accumulated wind mixing energy (m3/s3*day)
      real*8    DT            ! past time interval for accumulated values
      integer ND
      real*8    DEPTH (ND)
      integer INDXI  (NBI+1)
      real*8    BFXINT (NBI)  ! Net buoyancy resistance     (m2/s3*day)
C                integral of BFX from SURFBF over previous timestep.
      integer MLI
      real*8    VLAYER  (MLI)
C Updated:
      real*8    DENS   (MLI)
      real*8    SAL    (MLI)
      real*8    TEMP   (MLI)
C Out:
      real*8    XMIX   (NBI)
      integer LMIX   (NBI)
C          - index of lowest layer involved in homogenizing

C =================== mixing model parameters ====================

C      ---  Parameters for wind/buoyancy mixing:
      real*8 HHH_Ek /0.2 / ! Ekman depth constant, Stigebrandt 1985.
C
      real*8 F_Coriolis /1.26E-4/   ! (s-1)
      real*8 M0 /0.5/     ! Wind mixing effect coefficent
      real*8  SEC_PER_DAY
      PARAMETER (SEC_PER_DAY = 24.*3600.)
      real*8 GRAV /9.81/  !M/S2


C  ----------- LOCAL VARIABLES --------------

$if defined old_code
      INTEGER LB, LT, N_MIXED
$endif      

      INTEGER IB, L, ID, LSURF, LBOTTOM
      real*8 UFRIC, H_EKMAN, WD_ENERGY, HMIX_2
      
      
C .....................................................................

      real*8 ENERGY, WORK, H1, H2, DW, X, AVAILABLE_ENERGY, XR

$if defined DEBUG_SURFMX
      INTEGER I
      if (MXTEST) THEN
        WRITE ( DEBUG_UNIT,'('' ===== SURF_MIX.SURFMX, DT='',G15.7)') DT
        do ib=1,NBI
           CALL DENSITY_PRINT( INDXI(IB)+1, INDXI(IB+1),
     &                          SAL, TEMP, DENS )
        ENDDO
      endif
$else
      LOGICAL DUMMY
      DUMMY = MXTEST !TO AVOID WARNING
$endif



C =====================================================
C      Loop through basins:
      DO IB = 1,NBI
C =====================================================

C ---------------- Initialize counters:
         LSURF   = INDXI(IB)+1
         LBOTTOM = INDXI(IB+1)
         L  = LSURF
         ID = 1

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
        WRITE ( DEBUG_UNIT,*) '-------- SURFMX, BASIN:',IB
        CALL DENSITY_PRINT( LSURF, LBOTTOM, SAL, TEMP, DENS )
      endif
$endif


C --------------------------------------------------------
C     Initialize description of energy balance at surface:
C --------------------------------------------------------

C     ......  FR3INT: Wind mixing energy integrated during previous
C             timestep DT as 3'd power of frictional velocity.
         if (DT .LE. 0.0 ) THEN
            UFRIC = 0.0
            WD_energy = 0.0
         ELSE
            UFRIC = (FR3INT(IB)/DT)**(1./3.)
            WD_energy = M0 *FR3INT(IB)   ! m3/s3*day   - per surface
         ENDIF

C     ......  Ekman limit for wind-mixed depth:
         H_EKMAN   = HHH_Ek * UFRIC / F_Coriolis  ! m (Stigebrandt 1985)



C     ...... energy balance within surface layer:
C        .... Net gravitational energy required or available
C             as a result of homogenizing within surface layer:

         HMIX_2 = DEPTH(ID+1)/2.0  ! mean mixing depth  m
         WORK = BFXINT(IB)*HMIX_2  ! required work      m3/s3*day

C             Required work can be negative
C             due to evaporation or cooling at surface:
         IF ( WORK.lt.0.0) WORK = WORK*0.05
C             It is then assumed that 5% of this energy
C             is available to work against gradients further down.

C             The energy involved in mixing surface buoyancy
C             to greater depth will be included automatically
C             through layer densities, and is therefore
C             not explicitly included in the energy balance below.

C        ..... available mixing energy from wind
         ENERGY = WD_ENERGY


$if defined DEBUG_SURFMX
         if (MXTEST) THEN
            WRITE ( DEBUG_UNIT,*) '------ At surface:'
            WRITE ( DEBUG_UNIT,'(1X,2(A15:G15.7:))' )
     &          'FR3INT(IB)',FR3INT(IB), 'UFRIC', UFRIC,
     &          'WD_ENERGY', WD_ENERGY, 'H_EKMAN', H_EKMAN
            WRITE ( DEBUG_UNIT,*) 'BFXINT(IB), ENERGY, HMIX_2, WORK:'
            WRITE ( DEBUG_UNIT,*) BFXINT(IB),ENERGY, HMIX_2, WORK
         endif
$endif


C ----------------------------------------------------------
C    Initiate accumulative descriptors for well-mixed layer:
         VLAYER_MIXED  = 0.0  ! Signal to GET_DENS_MIXED
C ----------------------------------------------------------


C     ...... Layers further down included if energy is available,
C            from above or as instability.
         X = 1.0

C          Surface layer always assumed mixed here, but if X<1 is set
C          below, this has no effect on concentrations, since there is
C          only within-layer mixing, assumed in the model anyway.

C =================================================================
C     Mix layers until energy has been expended in gravity field:
         DO WHILE ( L .LT. LBOTTOM .AND. X .ge. 1.0 )
C =================================================================

C        ..... Calculate properties of well-mixed volume
C              with current layer to be included:
            CALL GET_DENS_MIXED ( L, SAL, TEMP, VLAYER )

$if defined DEBUG_SURFMX_GT_1
      if (mxtest) then
         write ( DEBUG_UNIT,'(1X,A5,3A16)')
     &          'Mixed to L:','SALINITY:','TEMPERATURE:','DENSITY:'
         write ( DEBUG_UNIT,'(1x,I5,3F16.8)' )
     &           L,SAL_MIXED(1), TEMP_MIXED(1), DENS_MIXED(1)
      endif
$endif

C        ...... Check energy balance for including next layer:
             L =  L + 1   ! Next layer to check
            ID = ID + 1

C ---------------------------------------------------------------------
C Find energy dW required against gravity field when mixing
C one layer (H2,r2,Z2) into well-mixed volume above (H1,r1,Z1).
C    [ H: thickness, r: density, Z: vertical position (>0 downwards) ]
C    dW = (energy after mixing) - (energy before mixing)
C       = [  (H1*r1 + H2*r2) * (-Z1*H1 -Z2*H2) / (H1 + H2)
C                           - H1*r1*(-Z1) - H2*r2*(-Z2)  ] * gravity
C reducing to:
C       = H1*H2*(r2-r1)*(z2-z1) / (H1+H2) * g
C or, since z2-z1 = H2/2:

C    dW = H1*H2**2/(H1+H2)/2*(r2-r1) * g
C                unit m2 *kg/m3 *m /m * m/s2 = kg/s2 = ( kgm2/s2/areal)
C scaled to kinetic unit m3/s3*day by /(s/day) /(kg/m3):
C    dW = H1*H2**2/(H1+H2)/2 *(r2-r1)/1000 * g /(s/day)
C ---------------------------------------------------------------------

C       ..... Applies this to layer L below well-mixed layers:
            H1 = DEPTH(ID)
            H2 = DEPTH(ID+1) - H1
            DW = H1 * H2**2 /(H1+H2) /2.0
     &             *(DENS(L)-DENS_MIXED(1))/1000.0 * GRAV /sec_per_day
C             unit:             m2 *  m/s2 /s *day  = m3/s3*day
C             (Density in sigma_t, 1000 = density in kg/m3)

C        .... If there is an energy gain because DENS(L)<DENS_MIXED,
C             assume 5% efficiency in counteracting density gradients
C             further down:
            IF (DW .LT. 0.0 ) THEN
                 DW = DW*0.05
            ENDIF


C             DW is now change in WORK if mixing completely,
C             while ENERGY is total expendable energy from wind.

C             If insufficient energy, the layers are mixed partly
C             to a fraction X between 0 and 1.
C             It is assumed that a fraction X of both layers
C             are mixed with each other, and then remixed
C             into the orginal layers.
C             This involves an energy X*DW, and results in
C             mixed concentrations:
C                       C1' = C1 + X*V2*(C2-C1)/(V1+V2)
C                       C2' = C2 + X*V1*(C1-C2)/(V1+V2)

C        .... Finds mixing fraction without regard to Ekman limitation:
            AVAILABLE_ENERGY = MAX(0.0D0, (ENERGY - WORK) )

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
         WRITE ( DEBUG_UNIT,*) ' Ekman limitation applied:'
         WRITE ( DEBUG_UNIT,'(3(1x,A,G12.5))' )
     &    ' H1', H1,'H2', H2,'H_Ekman',H_Ekman, 
     &      'DW', DW,'AVAILABLE_ENERGY', AVAILABLE_ENERGY, 'WORK',WORK
      endif
$endif
            
            
            IF(  DW .gt. AVAILABLE_ENERGY ) THEN
               X = MAX(0.0D0, MIN(1.0D0, AVAILABLE_ENERGY / DW ) )

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
         WRITE ( DEBUG_UNIT,*) ' Partial mix:  X=', X
      endif
$endif
            ELSE
               X = 1.0
            ENDIF

C        .... Applies Ekman limitation if no surplus energy
C             from negative buoyancy effects is left:
            IF ( H1+X*H2 .GT. H_EKMAN .AND. WORK .GE.0.0 ) THEN
               X = MAX(0.0D0, (H_Ekman - H1)/H2 )   ! >0

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
               WRITE ( DEBUG_UNIT,*) ' Ekman limitation--> X=', X
      endif
$endif
               Energy = 0.0    ! Wind mixing is now ineffective.
            
            ELSE IF ( WORK + X*DW .lt. 0.0 ) THEN
C              ... Positive energy balance without wind:
               IF ( DW .gt. 0.0 ) THEN  ! (WORK <0)
                   X = min (1.0D0, - WORK/DW)        ! --> X>0
                       ! ability to mix due to energy -WORK
               ELSE
                   X = 1.0
               ENDIF

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
               WRITE ( DEBUG_UNIT,*) ' WORK + X*DW>0, --> X=', X
      endif
$endif
            
            ENDIF


C        ..... Update net buoyant energy (Work done so far):
            WORK = WORK + X*DW
           if (X .gt. 1.0 .or. X.lt.0.0 )
     &              WRITE (*,*) 'Feil i SURFMX, IB,L,X=',IB,L,X

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST)
     &    WRITE ( DEBUG_UNIT,'(''  ENERGY:'',G15.7,''  WORK:'',G15.7)')
     &                        ENERGY, WORK
$endif

C ===========================================
         END DO  ! next layer
C        Exits loop with L = next layer (always >LSURF),
C                    and X = mixing fraction.
C        X>=1 only if L=LBOTTOM.
C ===========================================

C     ..... note preliminary index for lowest layer
C           which is completely homogenized with other layers.
         LMIX(IB) = L + MIN(X,1.0D0)-1
         LMIX(IB) = MAX( LSURF, LMIX(IB))
C                           Fractional part truncated.
C           may be further updated below.

C     .... Store number of well-mixed layers,
         XMIX(IB) = L + X - LSURF  !Fractional part included.
C            - full mixing of n layers, with n=integer part of XMIX,
C             n+1'th layer mixed with fraction part of X
C             ( used by MIX_CONC  in TRCALC)

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST)
     &    WRITE ( DEBUG_UNIT, '(1x, 2(A6,I5,a6,G15.7) )' )
     &        ' L:',L,' X:',X,' LMIX:', LMIX(IB),' XMIX:', XMIX(IB)
$endif

C     ....   Update sal, temp, density with last layer,
C            mixing fraction in X (OK even if X=0):
         IF ( L.GT. LSURF ) THEN
C        ...... mix in the specified fraction of next layer by adjusting
C           a. ... concentration C1 in well-mixed layers:
C             C1'= C1 + X*V2*(C2-C1)/(V1+V2) = C1*(1-X*R)+C2*X*R
C           b. ... concentration C2 in partially mixed layer:
C             C2'= C2 + X*V1*(C1-C2)/(V1+V2) = C1*X*R1   +C2*(1-X*R1)
C                                            = C1*X*(1-R)+C2*(1-X*(1-R))
C                  with R = V2/(V1+V2), R1 = V1/(V1+V2) = 1-R
            XR = X*VLAYER(L)/(VLAYER_MIXED + VLAYER(L))
C           ...... salinity:
            SAL(L-1)  = SAL_MIXED(1)*(1.0-XR) + SAL(L)*XR
            SAL(L)    = SAL_MIXED(1)*(X-XR)   + SAL(L)*(1.0-X+XR)
C           ...... temperature:
            TEMP(L-1) = TEMP_MIXED(1)*(1.0-XR) + TEMP(L)*XR
            TEMP(L)   = TEMP_MIXED(1)*(X-XR)   + TEMP(L)*(1.0-X+XR)
C           ...... density in layer L-1 and L, with shift of index:
            L = L-1

            CALL SIGMAT( SAL(L), TEMP(L), 2, DENS(L))
         ENDIF

C        ...... store final values in remaining mixed layers above:
         DO WHILE ( L .GT. LSURF )
            L = L - 1
            SAL(L)  =  SAL(L+1)
            TEMP(L) = TEMP(L+1)
            DENS(L) = DENS(L+1)
         END DO

C ------------ prepare integrated energy terms for next step in time:
C           Accumulated buoyancy energy debt (BFXINT>0)
C           is balanced against accumulated wind mixing energy,
C           gravitational mixing energy (BFXINT<0) is reset to zero:
         BFXINT(IB) = max ( 0.0D0, BFXINT(IB) - WD_ENERGY/HMIX_2 )
C           Integrated wind mixing energy is reset for each step,
C           so unusable energy is lost:
         FR3INT(IB) = 0.0

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
        WRITE ( DEBUG_UNIT,*) '     ---- Forced Mixing, XMIX=',XMIX(IB)
        CALL DENSITY_PRINT( LSURF, LBOTTOM, SAL, TEMP, DENS )
      endif
$endif


C ======================= FORCED MIXING DONE ===================

! ############################################################
$if defined old_code
! Taken out, no longer needed with new code for vertical mixing,
! which now handles instabilities OK
! seems to lead to problems sometimes, may create
! instabilies in other variables than sal, temp 


C ===================== HOMOGENIZE INSTABILITES ================
C   Scan layers from lower layer in wellmixed surface layer
C   to bottom layer.
      LT = LMIX(IB)
      DO WHILE ( LT.lt.LBOTTOM )
         LB = LT + 1
C ===========================================================


         IF ( DENS(LB) .LE. DENS(LT) ) THEN
C  -------- Instability:
C                Find extension of mixed region:
            DENS_MIXED(1) = DENS(LT)  ! Initial mixed density reference
            VLAYER_MIXED = 0.0     ! Initiating signal to GET_DENS_MIXED
            N_MIXED = 0            ! No layers mixed yet
C        ... Repeated attempts to expand in both directions:
            DO WHILE ( N_MIXED .lt. LB-LT )
               N_MIXED = LB-LT
               DO WHILE ( LB .le. LBOTTOM )
                  IF ( DENS(LB) .GT. DENS_MIXED(1) ) EXIT
                  CALL GET_DENS_MIXED ( LB, SAL, TEMP, VLAYER )
                  LB = LB + 1
               ENDDO  !LB=LBOTTOM+1 or layer LB denser than mixed region
               DO WHILE ( LT .ge. LSURF )
                  IF ( DENS(LT) .LT. DENS_MIXED(1) ) EXIT
                  CALL GET_DENS_MIXED (LT, SAL, TEMP, VLAYER )
                  LT = LT - 1
               ENDDO !LT=LSURF-1, or layer LT has density < mixed region
$if defined DEBUG_SURFMX_GT_1
      if (mxtest .and. LB-LT.gt.n_mixed) then
         write ( DEBUG_UNIT,'('' Mixed from '',2A8,3A15)')
     &          'LT+1','to LB-1','SALINITY:','TEMPERATURE:','DENSITY:'
         write ( DEBUG_UNIT,'(12x, 2I8, 3F15.7)' )
     &           LT+1,LB-1, SAL_MIXED(1), TEMP_MIXED(1), DENS_MIXED(1)
      endif
$endif

            ENDDO
C         .... Mixed region established:
C                  Consists of layers LT through LB-1,
C                  and is stable in relation to layers above and below.
C              Store homogenized values:
            DO L = LT+1,LB-1
               SAL(L)  = SAL_MIXED(1)
               TEMP(L) = TEMP_MIXED(1)
               DENS(L) = DENS_MIXED(1)
            ENDDO
C              Take note of lowest layer involved:
            LMIX(IB) = LB-1
         ENDIF


C ===================================================
C   Continue scanning of layers
C   from next layer/below mixed region:
         LT = LB
         N_MIXED = 0
      ENDDO
C ===================================================

$if defined DEBUG_SURFMX_GT_1
      if (MXTEST) THEN
         WRITE ( DEBUG_UNIT,*) '     ---- Homogenized:, LMIX=',LMIX(IB)
      endif
$endif

! ############################################################
$endif
! ############################################################


C ====================================================
      END DO  ! Basin
C ====================================================


$if defined DEBUG_SURFMX
      if (MXTEST) THEN
         WRITE( DEBUG_UNIT, '('' DT:'',G15.7/1x,A5,3a15)') DT,
     &     'LMIX','XMIX','FR3INT', 'BFXINT'
         WRITE( DEBUG_UNIT, '(I10, 3G15.7)')
     &     ( LMIX(I), XMIX(I), FR3INT(I), BFXINT(I), I=1,NBI )
         do ib=1,NBI
            CALL DENSITY_PRINT( INDXI(IB)+1, INDXI(IB+1),
     &                          SAL, TEMP, DENS )
         ENDDO
      endif
$endif


      contains


C ==================================================================
C Accumulate sal & temp of mixed layers and find resulting density.
      SUBROUTINE GET_DENS_MIXED (L,SAL,TEMP,VLAYER )


      integer L
      real*8 SAL(L), TEMP(L), VLAYER(L)

C ............................................................

C ....... initiate new mixing region:
      if ( VLAYER_MIXED .eq. 0.0 ) THEN
         SAL_MASS_MIXED = 0.0
         TEMP_MASS_MIXED = 0.0
      endif

$if defined DEBUG_SURFMX
      if(MXTEST) write(DEBUG_UNIT,*) 'GET_DENS_MIXED: Includes layer ', L, ' in mixed depth range'
$endif 

C ....... include specified layer:
      SAL_MASS_MIXED = SAL_MASS_MIXED  + VLAYER(L)*SAL(L)
      TEMP_MASS_MIXED = TEMP_MASS_MIXED  + VLAYER(L)*TEMP(L)
      VLAYER_MIXED = VLAYER_MIXED + VLAYER(L)

C ....... calculate properties of resulting mixed region:
      SAL_MIXED(1) = SAL_MASS_MIXED/VLAYER_MIXED
      TEMP_MIXED(1) = TEMP_MASS_MIXED/VLAYER_MIXED

      CALL SIGMAT( SAL_MIXED, TEMP_MIXED, 1, DENS_MIXED)

      END Subroutine

      End Subroutine



$if defined DEBUG_SURFMX
      SUBROUTINE DENSITY_PRINT( LSURF, LBOTTOM, SAL, TEMP, DENS )

      
      INTEGER LSURF, LBOTTOM
      real*8 SAL(*),TEMP(*),DENS(*)

      LOGICAL UNSTABLE
      INTEGER STAB, L, IU
      CHARACTER*14 TEXT(0:1)/'  ','Instability!'/

!      INCLUDE 'DEBUG.INC'

      write ( DEBUG_UNIT, '( 1X, A5, 3A16 )' )
     &          'I:', 'SALINITY:', 'TEMPERATURE:', 'DENSITY:'
      STAB = 0
      DO L = LSURF, LBOTTOM
         UNSTABLE =  (L.LT.LBOTTOM) .and. (DENS(L+1).lt.DENS(L))
         IF (UNSTABLE) THEN
             STAB = 2
         ELSE
             STAB = MAX(0,STAB-1)
         ENDIF
         IF ( STAB.LE.0 ) THEN
            IU = DEBUG_UNIT
         ELSE
            IU = 0
         ENDIF
         Do WHILE (IU .LE. DEBUG_UNIT)
            write ( IU,'(1x,I5,3F16.8,A14)' )
     &          L,SAL(L), TEMP(L), DENS(L), TEXT(MIN(STAB,1))
            IU = IU + DEBUG_UNIT
         ENDDO
      ENDDO
      END Subroutine


$endif

      end  Module