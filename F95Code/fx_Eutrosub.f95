      Module fx_Eutrosub
      
      use fx_RunControl
      
      use ModelParam_RunControl
      use ModelParam_Plankton     , only: NCZOO, PCZOO, NCBACT, PCBACT
      use ModelParam_Decomposition
      use ModelParam_InitState
      use ModelParam_Inputs
      use ModelParam_Mussels
      
      use ModelVar_Topography
      use ModelVar_HydroBioChem
      use ModelVar_RunControl 
      use fx_Mbalan
      use fx_Stoichiometry        ! by subroutine SUMIMP
      use fx_SigmaT               ! by subroutine CNCADJ
      use fx_TrCalc
      use fx_Surf_mix
      use fx_MusselIntegrate
      
      implicit none


      Contains

      
C Eutrofimodell Indre Oslofjord - Fil  EUTROSUB.FOR

C FORTRAN Subroutines performing subtasks of ACSL model.
C All subroutines accessing ACSL common block by include are
C found in this Module or in the main ACSL model source file.

$undefine special_check
$define DEBUG_SUMIMP
$define DEBUG_CNCADJ
$define debug_qcalc
$define debug_VOLCHK
$define debug_negative

$define mass_balance_deriv_check
          ! defined:   mass_balance subroutine checks derivatives
          ! undefined: mass_balance subroutine checks integrated values



C ------------------ monitoring integrated variables: -----------------
C Save mean values of integrated terms INTEGRAL as mean values in MEAN,
C and reset INTEGRAL to zero for integration of next period
      SUBROUTINE MeanValue ( INTEGRAL, DIVISOR, CURRENT_VALUE, MEAN )
      
      real*8 INTEGRAL, DIVISOR, CURRENT_VALUE, MEAN

      IF (DIVISOR .LE. 0.0) THEN
         MEAN = CURRENT_VALUE
      ELSE
         MEAN = INTEGRAL/DIVISOR
         INTEGRAL = 0.0
      ENDIF
      END SUBROUTiNE


      SUBROUTINE MeanVector(N, INTEGRAL, DIVISOR, CURRENT_VALUE, MEAN)
      INTEGER N
      real*8 INTEGRAL(N), DIVISOR, CURRENT_VALUE(N), MEAN(N)

      integer I

      DO I=1,N
         IF (DIVISOR .LE. 0.0) THEN
             MEAN(I) = CURRENT_VALUE(I)
         ELSE
             MEAN(I) = INTEGRAL(I)/DIVISOR
             INTEGRAL(I) = 0.0
         ENDIF
      ENDDO

      END SUBROUTiNE


C =================================================================
      SUBROUTINE COMCLC (MODE, MBI, MLI )

C        Calculations performed at the beginning of simulation,
C        and at the end of each communication interval.
C        (Now: at end of each integration interval)

      INTEGER MODE ! Controls when to perform check:

          ! =0 Only if subroutine checks integrated values
          !    ( mass_balance_deriv_check is undefined)
          ! =1 Also if subroutine checks derivatives
          !    ( mass_balance_deriv_check is defined)

      INTEGER MBI, MLI ! Dimensioning number of basins and layers

!       INCLUDE    'EUTRO.INC'

      INTEGER I_L,I_B

C --------- Check conservation of heat and mass substances:
C           Nitrogen combined as total nitrogen NTOT

$IF !DEFINED mass_balance_deriv_check

      LOGICAL SIZE_SCALE

      ! ------- Check accumulated mass vs. net import.
      !         Accumulate amount in system (water + bottom):

      SIZE_SCALE = .true. ! Will update ....MZ as max. sum of abs.val.


 10   CALL MBALAN ( MBI, MLI, SIZE_SCALE,
     &        SAL  , TEMP , OXYG , NO3  , NH4   , PO4 , SiO2,
     &        CFYT , NFYT , PFYT , SFYT ,
     &        CDET , NDET , PDET , SDET ,
     &        CSED , NSED , PSED , SSED , ASED  , PADS,
     &        ODM  , DOC  , BACT , CZOO , CTMUSL,
     &        SPP  , SPPSED,
     &        SALTB, HEATB, OXYGB, NITRB, PHOSB , SILIB, SPPB )

      ! SIZE_SCALE: Updates size scale .......MZ
      ! .not.SIZE_SCALE:  Balances current amount in each basin against
      !      values found by integrating import from initial contents:
      !      resulting values should be zero within accuracy
      CALL MACCUR ( SIZE_SCALE, MBRSET, 1, 'SALT',
     &              SALTMI, SALTMZ,  SALTB )
      CALL MACCUR ( SIZE_SCALE, MBRSET, 2, 'HEAT',
     &              HEATMI, HEATMZ,  HEATB )
      CALL MACCUR ( SIZE_SCALE, MBRSET, 3, 'OXYGEN',
     &              OXYGMI, OXYGMZ,  OXYGB )
      CALL MACCUR ( SIZE_SCALE, MBRSET, 4, 'NITROGEN',
     &              NITRMI, NITRMZ,  NITRB )
      CALL MACCUR ( SIZE_SCALE, MBRSET, 5, 'PHOSPHORUS',
     &              PHOSMI, PHOSMZ,  PHOSB )
      CALL MACCUR ( SIZE_SCALE, MBRSET, 6, 'SILICATE',
     &              SILIMI, SILIMZ,  SILIB )
!      CALL MACCUR ( SIZE_SCALE, MBRSET, 7, 'Particles',
!     &              SPPMI, SPPMZ,  SPPB )


      if (SIZE_SCALE) THEN
          SIZE_SCALE=.false.
          GOTO 10
      ENDIF

      IF (MODE.eq.0) goto 100

$else
           ! Check derivative budget balance, only if called
           ! at end of CNCADJ, where derivatives are consistent.

      IF (MODE.eq.0) goto 100

! .......... Set size scale:
      CALL MBALAN ( MBI, MLI, .true.,
     &        SALDV, TEMPDV, OXYGDV, NO3DV, NH4DV, PO4DV, SiO2DV,
     &        CFYTDV, NFYTDV , PFYTDV , SFYTDV,
     &        CDETDV, NDETDV , PDETDV , SDETDV,
     &        CSEDDV, NSEDDV,  PSEDDV, SSEDDV, ASEDDV, PADSDV,
     &        ODMDev, DOCDV,   BACTDV , CZOODV, CMUSDV,
     &        SPPDV , SPPSEDDV,
     &        SALTMZ, HEATMZ, OXYGMZ, NITRMZ, PHOSMZ, SILIMZ, SPPMZ )

! ......... Accumulate derivatives within system (water + bottom)
      CALL MBALAN ( MBI, MLI, .false.,
     &        SALDV, TEMPDV, OXYGDV, NO3DV, NH4DV, PO4DV, SiO2DV,
     &        CFYTDV , NFYTDV , PFYTDV , SFYTDV,
     &        CDETDV , NDETDV , PDETDV , SDETDV,
     &        CSEDDV,  NSEDDV,  PSEDDV, SSEDDV, ASEDDV, PADSDV,
     &        ODMDev, DOCDV,   BACTDV , CZOODV, CMUSDV,
     &        SPPDV , SPPSEDDV,
     &        SALTB, HEATB, OXYGB, NITRB, PHOSB, SILIB, SPPB )
C           (subroutine MBALAN in MBALAN.FOR)

C   ....... Balance current amount in each basin against
C           values found by integrating import from initial contents:
C           resulting values should be zero within accuracy
      CALL MACCUR ( .false., MBRSET, 1, 'SALT'       ,
     &              SALTMP, SALTMZ,  SALTB )
      CALL MACCUR ( .false., MBRSET, 2, 'HEAT'       ,
     &              HEATMP, HEATMZ,  HEATB )
      CALL MACCUR ( .false., MBRSET, 3, 'OXYGEN'     ,
     &              OXYGMP, OXYGMZ,  OXYGB )
      CALL MACCUR ( .false., MBRSET, 4, 'NITROGEN'   ,
     &              NITRMP, NITRMZ,  NITRB )
      CALL MACCUR ( .false., MBRSET, 5, 'PHOSPHORUS' ,
     &              PHOSMP, PHOSMZ,  PHOSB )
      CALL MACCUR ( .false., MBRSET, 6, 'SILICATE'    ,
     &             SILIMP, SILIMZ,  SILIB )
!      CALL MACCUR ( .false., MBRSET, 7, 'Particles',
!     &              SPPMP, SPPMZ,  SPPB )

$endif

      MBRSET = .FALSE. ! Must be activated each time from main model

C ----------- Reset very small buffer volumes in VBUF to zero,
C             to avoid confusing ACSL plots (values X where 0<X<1.0e-38
C             displays as 'R', plots as 1.0e30)
 100  DO I_L = 1,INDXC(NC+1)
          DO I_B=1,2
             IF ( ABS(VBUF(I_B,I_L)).LT.1E-36 ) VBUF(I_B,I_L) = 0.0
          ENDDO
      ENDDO

      END SUBROUTiNE


C =================================================================
C Control mass balance, updating size scale:
      SUBROUTINE MACCUR( SIZE_SCALE, RESET, B_M, V_NAME,
     &                   V_intg, V_scale, V_accuracy)
      LOGICAL SIZE_SCALE, RESET
      integer B_M
      CHARACTER V_NAME*(*)

!       INCLUDE    'EUTRO.INC'

      real*8 V_intg  (NBI) !  Initial + integrated net imports
      real*8 V_scale (NBI) !  Accumulated size scale
      real*8 V_accuracy (NBI) ! In: Summed amount of component in model
C                             Out: Accuracy measure

      INTEGER I_V
      real*8 Acc, Scale
      integer PRINT_OUT

      logical WARNING_given /.false./
      character*26 Warning_File /'MSB_Warn.log'/ 
      integer warning_UNit/991/, warning_count /0/, warn_error
      
!      INCLUDE   'DEBUG.INC'


      if (SIZE_SCALE) then
         DO I_V=1,NBI
            V_scale(I_V) = Max( V_scale(I_V) ,V_accuracy(I_V) )
         END DO
      else
         if (T.ge.TTRIG) THEN
            PRINT_OUT = MBPRT(B_M)
         else
            PRINT_OUT = 0
         ENDIF

         IF (PRINT_OUT.gt.0) THEN
            WRITE(DEBUG_UNIT,'('' >>>> MASS CONSERVATION: '', A)')V_NAME
            WRITE ( DEBUG_UNIT, '(1x,A5,4A15)' )
     &              'num.','size scale','integrated',
     &              'summed','accuracy'
         ENDIF

         if (RESET) warning_given = .false.
         DO I_V=1,NBI
            if (RESET) then
               V_intg(I_V) = V_accuracy(I_V) ! Reinitialize mass balance
            ELSE
               Acc = (V_intg(I_V) - V_accuracy(I_V))
               Scale = Max ( abs(V_intg(I_V)),
     &                       abs(V_accuracy(I_V)),
     &                       abs(V_scale(I_V))     )
               if ( Acc .ne. 0.0 ) Acc = Acc / Scale
                    ! May have Scale = 0,
                    ! but then Acc should be zero as well.
               IF (PRINT_OUT.gt.0) THEN
                   WRITE ( DEBUG_UNIT, '(1x,I5,4G15.7)' )
     &                  I_V, V_scale(I_V), V_intg(I_V),
     &                       V_accuracy(I_V), Acc
               ENDIF
               V_accuracy(I_V) = Acc

               IF ( ABS(Acc) .gt. ACCUR .and. MBPRT(B_M).ge.0 ) THEN

                  if (.not.warning_given) then
                     warning_given = .true.
                     open (Warning_Unit, File = Warning_File,
     &                     IOSTAT=warn_error)
                         ! if error: will try FORnnn instead
                  endif

                  write(*,'('' Warning:'',A,
     &                      '' conservation deviates by '',
     &                         G15.7)') V_NAME, ACC

                  write(warning_unit,'(3A,G16.8,A,I5,A,G16.8)')
     &              ' Warning: ',
     &              V_NAME, ' conservation deviates by ', Acc,
     &              ' in basin ', I_V,' at T=', T

                  if (warning_count.ge.200) then
                     write (*,'(A,I5,A)') ' Stops after more than ',
     &                                   warning_count,' warnings'
                     stop
                  endif
                  warning_count = Warning_count+1

               ENDIF
            ENDIF
         END DO
      ENDIF

      END SUBROUTiNE


C ==================================================================
C Calculations performed at the start of each integrated interval,
C and at the end of each communication interval:


C ------ Volume of layers:
C Part 2 of these calculations could be confined to initiating part:

      SUBROUTINE VLCALC (D_VTOT)

!       INCLUDE    'EUTRO.INC'


      real*8 D_VTOT(NBI)

      INTEGER I_B, L_X
      real*8  X

      DO I_B = 1,NBI

      !        Relative correction factor to volume of open layers:
         VLCORF(I_B) = 1.0 + D_VTOT(I_B)/(VTOTZ(I_B)*VFROPN(I_B))

      !        Total volume corrected by this fraction:
         X = VTOTZ(I_B)*VLCORF(I_B)

C         write(*,*) 'I_B:', I_B, 'nlvopn:',nlvopn(I_B),
C     &              'vfropn:',VFROPN(I_B),
C     &              'VTOTZ:',VTOTZ(I_B),'VLFAC:', VLCORF(I:B)

C     ..... 1. Volume of open layers, fixed + variable part
         DO L_X = INDXI(I_B)+1, INDXI(I_B)+ NLVOPN(I_B)
            VLAYER(L_X) = VFRAC(L_X)*X

C            write(*,*)'L_X:',L_X, 'VFRAC:', VFRAC(L_X),
C     &                'VLAYER:', VLAYER(L_X)

         ENDDO
C     ..... 2. Volume of closed layers fixed values:
         DO L_X = INDXI(I_B)+NLVOPN(I_B)+1,INDXI(I_B+1)
            VLAYER(L_X) = VFRAC(L_X)*VTOTZ(I_B)
         ENDDO
         ! Sum over all layers should now be VTOTZ(I_B) + D_VTOT(I_B)
      ENDDO


      END SUBROUTiNE


C -------- Update water and sediment concentrations over one time step:
      SUBROUTINE CNCADJ ( MBI, MLI, MSAGES, MSLAYR )
      integer  MBI, MLI, MSAGES, MSLAYR


!       INCLUDE    'EUTRO.INC'

      real*8 TIME_STEP
      LOGICAL RESTART_INTEGRAL
      COMMON / CNCADJ_COMMON/ TIME_STEP, RESTART_INTEGRAL

      real*8 T_BALANCE
      SAVE T_BALANCE

      INTEGER PRINT_UNIT

      INTEGER I_B, I_L, N_GROUP, L_X, F_G
      real*8 Part_C_SUM, Part_N_SUM, Part_P_SUM

      character nchar_group



      CALL HELLO ( 'CNCADJ'    )


      RESTART_INTEGRAL = .false.

      TINTEG = MAX( TINTEG, T)

      IF ( TINTEG .GT. TDERIV ) THEN


C    Difference between integrated time and time at last call to CNCADJ,
C    other state variables are updated through this time-step before
C    derivatives are calculated again:

         TIME_STEP = TINTEG-TDERIV
         call Integrate_Physical_Variables(TIME_STEP)

$if defined  DEBUG_CNCADJ
         IF( TRTEST .and. T.ge.TTRIG ) THEN
            WRITE(DEBUG_UNIT,'(1X,A,G12.5)')  'CNCADJ at T=' ,T
            WRITE(DEBUG_UNIT,'(1X,A5,3A15)')  'IB:' ,'DVTOTL ' ,'DVTOT '
            WRITE( DEBUG_UNIT, '(1X,I5,2G15.7)')
     &            (I_B, DVTOTL(I_B), DVTOT(I_B), I_B = 1,NBI)
         ENDIF
$endif

C     .... Adjust salt, temperature and particles for time and volume changes:
         CALL ADJ_CONC ( TRCALC, MDEBUG(2), 'SAL', NBI,
     &          DVTOT, DVTOTL, VTOTZ, INDXI, NLI,
     &          VLAYER, NLVOPN, VFROPN,
     &          TIME_STEP, SALDV, SAL, SAL )
         CALL ADJ_CONC ( TRCALC, MDEBUG(2), 'TEMP', NBI,
     &          DVTOT, DVTOTL, VTOTZ, INDXI, NLI,
     &          VLAYER, NLVOPN, VFROPN,
     &          TIME_STEP, TEMPDV, TEMP, TEMP )
     
! Effect of particles on density not yet implemented;
! will require SPP to be updated in parallell with Sal and Temp:     
!         CALL ADJ_CONC ( TRCALC, MDEBUG(2), 'SPP', NBI,
!     &          DVTOT, DVTOTL, VTOTZ, INDXI, NLI,
!     &          VLAYER, NLVOPN, VFROPN,
!     &          TIME_STEP, SPPDV, SPP, SPP )
      ENDIF  ! subroutine ADJ_CONC in TRCALC.FOR

! Note: Mixed surface layer due to instability found in SURFMX called below,
!       where SAL, TEMP is also changed accordingly.
!       For the other variables, homogenizing over well-mixed layer is done
!       in subroutine MIX_CONC, called by CCONC.
! When SPP is included in calculatio of density, it must be updated 
! in the same way as SAL and TEMP; vy ADJ_CONC and SURFMX, instead of by CCONC

C --------------------------------------------------------------------
C ALWAYS: Calculate density profile for adjusted salt, temperature
C         and particle concentration:

! ******** new code:
!      CALL SIGMAT( SAL, TEMP, SPP, dDens_dSPP, INDXI(NBI+1), DENSI)
! ******** old code:
      CALL SIGMAT( SAL, TEMP, INDXI(NBI+1), DENSI)

C --------------------------------------------------------------------

C --------------------------------------------------------------------
C ALWAYS: Update mussel model, and return sum of mussel carbon:
C (Necessary checks on repeated calls within subroutine)

      CALL MUSSEL_INTEGRATE (MSTEST.and.(T.ge.TTRIG), TINTEG, TDERIV,
     &     NBI, MSAGES, MSLAYR,
     &     MUSLNR, MUSLMT, MUSLMA, MUSLGT, MUSLWM,
     &     CTMUSL, CAMUSL, C0MUSL, CMUSDV,
     &     MUSLWA, MUSLWT )
C Dimension of NRMUSL, WAMUSL and WTMUSL was transferred in MUSLINIT
C in EUTRO.SRC
C --------------------------------
      
      CALL Hello ('CNCADJ: Returned from Mussel_Integrate')

C If time has not been increased: Return directly:
      IF ( TINTEG .LE. TDERIV ) GOTO 999


C -------------------------------------------------------------------
C Time increased: Update integrated variables:
C -------------------------------------------------------------------
      CALL Hello ('CNCADJ: starts updating variables')

C  .... Add mixing energy given by FR3INT to total sum over communication interval: 
C       in ACSL/WATCOM version: CALL ADD ( NBI, FR3INT, 1.0, UFR3I )
      DO I_B=1,NBI
         UFR3I (I_B) = UFR3I(I_B) + FR3INT(I_B)
      END DO

      
      CALL Hello ('added mixing energy')
      
      IF (TRACE) THEN
$if defined DEBUG_CNCADJ
         DO PRINT_UNIT = 6, DEBUG_UNIT, MAX(1,DEBUG_UNIT - 6)
$else
         DO PRINT_UNIT = 6, 6 , 1
$endif
           write(PRINT_UNIT,'(1x,4a15/1x,4G15.7/1x,3A15)')
     &          'TDERIV','TINTEG:','DT:','WNDSPD:',
     &           TDERIV,  TINTEG , TINTEG - TDERIV, WNDSPD,
     &          'UFRIC3', 'FR3INT:',  'UFR3I:'
           write(PRINT_UNIT,'(1x,3G15.7)')
     &            (UFRIC3(I_B), FR3INT(I_B),  UFR3I(I_B), I_B = 1,NBI)
         ENDDO
      ENDIF


C  .... Get depth of well-mixed volume by balancing accumulated
C       mixing energy against increase in potenial energy.

C       Other depth intervals may require homogenizing because of
C       nonlinear density function of (T,S).
C       This will show by DENSI being homogeneous over such layers,
C       and indexes of the lowest layer involved in such mixing
C       are noted in LMIX(IB), with IB = basin index.

      CALL Hello ('calls SURFMX')

! new call:
!      CALL SURFMX ( MXTEST, FR3INT, T-TDERIV, ND, DEPTH, NBI, INDXI,
!     &              BFXINT, NLI, VLAYER, DENSI, SAL, TEMP,
!     &              SPP, dDens_dSPP, XMIX, LMIX )
! old call:
      CALL SURFMX ( MXTEST, FR3INT, T-TDERIV, ND, DEPTH, NBI, INDXI,
     &              BFXINT, NLI, VLAYER, DENSI, SAL, TEMP,
     &              XMIX, LMIX )
       ! write(*,*) 'etter SURFMX: XMIX=', XMIX
      
      CALL Hello ('returned from SURFMX')

C  NOTE!  FR3INT and BFXINT is reset within SURFMX,
C         i.e. initialized for integrating next step.


C --------  Adjust state variables for changing volumes
C           and by homogenizing according to results in SURFMX above:

C          NOTE! This list of calls should include
C                all pelagic state variables
C                which are integrated concentrations in water


            ! ---------------------------------------------------
            ! Special state variable C1 for continuity check
            ! or accumulation of residence time:

      if (C1XTRN.eq.0.0) THEN
         DO I_l =1,NLI 
            C1DV(I_l) = C1DV(I_L)+1./365. ! (Array operation)
         ENDDO
      ENDIF
             ! if C1XTRN == 0.0, new water has C1 = 0,
             ! increases with time inside the model layers,
             ! giving mean residence time (or "age") of the water
             ! within the model volumes (regardless of where).
             ! may be modified below by C1ZERO.

      CALL Hello ('calls CCONC')
      CALL CCONC ( 1, 'C1',   TIME_STEP, C1DV,   C1,   C1   )

      
             ! if C1XTRN == 0.0, C1 is reset to 0 for basins 
             ! for which C1ZERO = 0, so residence time is within the 
             ! other basins only.

      if (C1XTRN.eq.0.0) THEN
	      DO I_B = 1, NBI
	         if (C1ZERO(I_B).eq.0.0) THEN
	             DO I_L = INDXI(I_B)+1, INDXI(I_B+1), 1
	                C1(I_L) = 0.0
	             Enddo
	         Endif
	      Enddo
      Endif
            ! ---------------------------------------------------


      
      CALL CCONC ( 3, 'OXYG', TIME_STEP, OXYGDV, OXYG, OXYG )

!      IF (T.GT.600) THEN
!          WRITE(999,"(1X,A,4G12.5)")
!     &       'F›r CCONC: PO4(34..37):',(po4(i_p),I_P=34,37),
!     &       'PO4DV(34..37):',(po4DV(i_p),I_P=34,37)
!      ENDIF

      CALL CCONC ( 4, 'PO4',  TIME_STEP, PO4DV,  PO4,  PO4  )

$if defined  SPECIAL_CHECK
      call check( NLI, PO4, T )
$endif

!      IF (T.GT.600) THEN
!          WRITE(999,"(1X,A,4G12.5)")
!     &       'Etter cconc: PO4(34..37):',(po4(i_p),I_P=34,37),
!     &       'PO4DV(34..37):',(po4DV(i_p),I_P=34,37)
!      ENDIF

      CALL CCONC ( 4, 'NO3',  TIME_STEP, NO3DV,  NO3,  NO3  )
      CALL CCONC ( 4, 'NH4',  TIME_STEP, NH4DV,  NH4,  NH4  )
      CALL CCONC ( 4, 'SiO2', TIME_STEP, SiO2DV, SiO2, SiO2 )


      DO N_GROUP = 1, FYTGRP
         NCHAR_GROUP = CHAR(N_GROUP+ICHAR('0') )
         CALL CCONC ( 5,'CFYT'//NCHAR_GROUP, TIME_STEP,
     &          CFYTDV(1,N_GROUP), CFYT(1,N_GROUP), CFYT(1,N_GROUP) )
         CALL CCONC ( 5,'NFYT'//NCHAR_GROUP, TIME_STEP,
     &          NFYTDV(1,N_GROUP), NFYT(1,N_GROUP), CFYT(1,N_GROUP) )
         CALL CCONC ( 5,'PFYT'//NCHAR_GROUP, TIME_STEP,
     &          PFYTDV(1,N_GROUP), PFYT(1,N_GROUP), CFYT(1,N_GROUP) )
         CALL CCONC ( 5,'CHL' //NCHAR_GROUP, TIME_STEP,
     &          CHLDV(1,N_GROUP) ,  CHL(1,N_GROUP), CFYT(1,N_GROUP) )
      ENDDO
      CALL CCONC ( 5,'SFYT', TIME_STEP, SFYTDV, SFYT, CFYT )
      CALL CCONC ( 5,'ODM' , TIME_STEP, ODMDev, 
     &             ODM, ODM   )
      CALL CCONC ( 5,'DOC ', TIME_STEP, DOCDV,  DOC,  DOC  )
      CALL CCONC ( 5,'BACT', TIME_STEP, BACTDV, BACT, BACT )
      
      
      CALL CCONC ( 6,'CZOO', TIME_STEP, CZOODV, CZOO, CZOO )

       ! ........... Detritus in water:

      CALL CCONC ( 5,'CDET', TIME_STEP, CDETDV, CDET, CDET )
      CALL CCONC ( 5,'NDET', TIME_STEP, NDETDV, NDET, CDET )
      CALL CCONC ( 5,'PDET', TIME_STEP, PDETDV, PDET, CDET )
      CALL CCONC ( 5,'SDET', TIME_STEP, SDETDV, SDET, CDET )
      CALL CCONC ( 5,'RDET', TIME_STEP, RDETDV, RDET, CDET )


         ! %%%%%%%%%%%%%%%%% particles %%%%%%%%%%%%%%%%%%%%%
         ! Added March 2001, BBJ
      CALL CCONC ( 7,'SPP', TIME_STEP, SPPDV, SPP, SPP )
         ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C ADJ_CONC includes update due to time derivative, thus
C avoiding overhead connected with defining water concentrations
C as ACSL (INTEGrated) state variables.
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



C Take copy of total volumes to use next time (in subr. CCONC),
C and note time from last call of CNCADJ to avoid double correction.
!!      CALL VCOPY ( NBI, DVTOT, DVTOTL ):
      DO I_B=1,NBI
          DVTOTL (I_B) = DVTOT(I_B)
      END DO



C =======================================================



      L_X = INDXI(NBI+1)

C --------- UPDATE SEDIMENT CONCENTRATIONS:
C     ......... Carbon:
      CALL INTEGRATE ( L_X, CSED, CSEDDV )
C     ......... Carbon*degradation rate:
      CALL INTEGRATE ( L_X, RSED, RSEDDV )
C     ......... Nitrogen:
      CALL INTEGRATE ( L_X, NSED, NSEDDV )
C     ......... Phosphorus:
      CALL INTEGRATE ( L_X, PSED, PSEDDV )
C     ......... Silicate:
      CALL INTEGRATE ( L_X, SSED, SSEDDV )

C %%%%%%%%%%%%%%%% particles:
      CALL INTEGRATE ( L_X, SPPSED, SPPSEDDV )

C     ......... Buffered oxygen depth in sediments:
      CALL INTEGRATE ( L_X, ASED, ASEDDV)
C     ......... Phosphorus (adsorbed in sediments):
      CALL INTEGRATE ( L_X, PADS, PADSDV )

C     .. Residual fraction of organic carbon at start of simulation
C     (mass balance not needed)

      DO I_L = 1, L_X
         XSED(I_L) = XSED(I_L) + XSEDDV(I_L)*TIME_STEP
      ENDDO

      DO I_B = 1, NBI
         XBUR(I_B) = XBUR(I_B) + XBURDV(I_B)*TIME_STEP
      ENDDO

C ........... integrated mass changes and fluxes:

C --------- Integrate particulate matter fluxes as mg/m2:


      if (REINTG) THEN  ! Store new initial value of integral:

         DO I_B = 1, NBI

            SALTMI(I_B,2) = SALTMI(I_B,1) ! USED IN WRTREP.FOR
            HEATMI(I_B,2) = HEATMI(I_B,1)
            OXYGMI(I_B,2) = OXYGMI(I_B,1)
            NITRMI(I_B,2) = NITRMI(I_B,1)
            PHOSMI(I_B,2) = PHOSMI(I_B,1)
            SILIMI(I_B,2) = SILIMI(I_B,1)
         ENDDO

         TINTGZ = TDERIV  ! New start_point of integration

      ENDIF


      CALL INTEGRATE (NBI, SALTMI, SALTMP )
      CALL INTEGRATE (NBI, HEATMI, HEATMP )
      CALL INTEGRATE (NBI, OXYGMI, OXYGMP )
      CALL INTEGRATE (NBI, NITRMI, NITRMP )
      CALL INTEGRATE (NBI, PHOSMI, PHOSMP )
      CALL INTEGRATE (NBI, SILIMI, SILIMP )


      ! Zero the other integrals:
      restart_integral = REINTG ! Control switch from EUTRO.CSL
      REINTG = .FALSE.

      CALL INTEGRATE (NBI, DNITRI, DENITR )
      CALL INTEGRATE (NBI, NFIXI , NFIX   )
      CALL INTEGRATE (NBI, CLOADI, CLOAD  )
      CALL INTEGRATE (NBI, ODMLOADI, ODMLOAD  )
      CALL INTEGRATE (NBI, NLOADI, NLOAD  )
      CALL INTEGRATE (NBI, PLOADI, PLOAD  )
      CALL INTEGRATE (NBI, CSEDXI, CSEDXP )
      CALL INTEGRATE (NBI, NSEDXI, NSEDXP )
      CALL INTEGRATE (NBI, PSEDXI, PSEDXP )

      CALL INTEGRATE ( L_X, CDFLXI, CDFLUX )
      CALL INTEGRATE ( L_X, NDFLXI, NDFLUX )
      CALL INTEGRATE ( L_X, PDFLXI, PDFLUX )
      CALL INTEGRATE ( L_X, SDFLXI, SDFLUX )

$IF !DEFINED mass_balance_deriv_check
      CALL VLCALC (DVTOT)
$endif

      IF ( T_BALANCE + MBINTV .LE. TINTEG) THEN
          CALL COMCLC (1, MBI, MLI )
          T_BALANCE = TINTEG
      ENDIF

$if defined  mass_balance_deriv_check
      CALL VLCALC (DVTOT)
$endif


C      WRITE(999,"(1X,7(G14.7,','))") Tinteg, VLAYER(1),TEMP(1),
C     &               tempdv(1), HEATMP(1), HEATMI(1,1), HEATB(1)

C ........ always reset TDERIV (important after restart of simulation)

  999 TDERIV = TINTEG
      call hello ('reached 999 in CNCADJ')
C          and accumulate particulate and total C, N and P:

      DO L_X = 1, INDXI(NBI+1)
            ! Initiate the sum of particulate N, C and P with
            ! contributions from zooplankton and detritus:
         Part_C_SUM = CZOO(L_X)        + CDET(L_X)
         Part_N_SUM = CZOO(L_X)*NCZOO  + NDET(L_X)
         Part_P_SUM = CZOO(L_X)*PCZOO  + PDET(L_X)
            ! Add contribution from phytoplankton:
         DO F_G = 1, FYTGRP
             Part_C_SUM = Part_C_Sum + CFYT(L_X,F_G)
             Part_N_SUM = Part_N_Sum + NFYT(L_X,F_G)
             Part_P_SUM = Part_P_Sum + PFYT(L_X,F_G)
         ENDDO
            ! Total concentrations include bacteria
            ! and dissolved matter:
         PartC(L_X) = Part_C_SUM
         TotC (L_X) = Part_C_SUM +BACT(L_X)        +DOC(L_X)
         PartN(L_X) = Part_N_SUM
         TotN (L_X) = Part_N_SUM +BACT(L_X)*NCBACT +NO3(L_X) +NH4(L_X)
         PartP(L_X) = Part_P_SUM
         TotP (L_X) = Part_P_SUM +BACT(L_X)*PCBACT +PO4(L_X)
      ENDDO
      call hello ('exits from CNCADJ')

      END SUBROUTiNE


C ==================================================================
      SUBROUTINE INTEGRATE ( N, V_1, DV )
      integer N
      real*8 V_1 (N) ! Integrated value
      real*8 DV  (N) ! Time derivative

      real*8 TIME_STEP
      LOGICAL RESTART_INTEGRAL
      COMMON / CNCADJ_COMMON/ TIME_STEP, RESTART_INTEGRAL

      integer I

      if ( .not. RESTART_INTEGRAL ) THEN
         DO I = 1, N
            V_1(I) = V_1(I) + DV(I)*TIME_STEP
         ENDDO
      Else
         DO I = 1, N
            V_1(I) = DV(I)*TIME_STEP
         ENDDO
      Endif
      END Subroutine



C ------------ ADJUST CONCENTRATIONS (TRCALC.CALC_CONC):
      SUBROUTINE CCONC ( DEBUG_INDEX, V_NAME,
     &           TIME_STEP, V_DV, V_CONC, CTRL_C )
      INTEGER DEBUG_INDEX
      CHARACTER*(*) V_NAME

!       INCLUDE    'EUTRO.INC'

      real*8 TIME_STEP
      real*8 V_DV(NLI), V_CONC(NLI), CTRL_C(NLI)

C  ----- Adjust concentration for time change and volume adjustments:
      CALL ADJ_CONC ( TRCALC, MDEBUG(DEBUG_INDEX),
     &      V_NAME, NBI, DVTOT, DVTOTL, VTOTZ, INDXI, NLI,
     &      VLAYER, NLVOPN, VFROPN,
     &      TIME_STEP, V_DV, V_CONC, CTRL_C )

$if defined  DEBUG_NEGATIVE
      if (V_NAME.ne.'OXYG')
     &        call check_SIGN( T, NLI, V_NAME, V_DV, V_CONC)
$endif

C  ----- Homogenize over well_mixed layers,
C        indicated by homogeneous density.
          CALL MIX_CONC ( MDEBUG(DEBUG_INDEX), V_NAME,
     &                    NBI, INDXI, XMIX, LMIX, NLI, VLAYER,
     &                    DENSI, V_CONC )


      END SUBROUTiNE


$if defined  DEBUG_NEGATIVE
C =================== check procedure =============================

      subroutine check_SIGN(T, NLI, V_NAME, V_DV, V_CONC)

      real*8 T
      integer NLI

      CHARACTER*(*) V_NAME

      real*8 V_DV(NLI), V_CONC(NLI)

      integer I,NERR
      logical name_printed
      logical file_opened /.false./
      integer neg_unit /6/
      integer count, index(3), N
      integer total_count /0/
      real*8 C_NEG(3), DV_NEG(3)

      name_printed =.false.
      count = 0

      DO I = 1, NLI

         if ( V_CONC(I).lt. - 0.01 ) then
             if (.not. file_opened) then
                  open(unit = 998, file='Negative.val',IOSTAT=nerr)
                  write(*,*)'File "Negative.val" is opened for reporting concentrations < 0'
!                  Pause 'Press Enter to continue'

                  if (nerr.eq.0) then
                       neg_unit=998
                  endif
                  file_opened=.true.
             endif

             count = count + 1
             total_count = total_count + 1
             C_neg(Count) = V_CONC(I)
             DV_Neg(Count) = V_DV(I)
             index(count) = I
         endif
         if (count.eq.3.or.(count.gt.0.and.I.eq.NLI)) then
            if (.not. name_printed) then
                write( neg_unit,
     &                    '('' T='',F10.4, A6)' ) T, V_NAME
                name_printed=.true.
            endif
            write( neg_unit,
     &             '(3(''  /'',I2,2('':''G9.2)))' )
     &           ( INDEX(N), C_neg(N),DV_NEG(N) , N=1,COUNT )
            count = 0
!            if (neg_unit.eq.6) Pause
            if (total_count.gt.10000) then
              write(*,*) 'Reported 10000 negative values to file "negative.val"'
!               Pause 'Press Enter to continue'
            Endif
         endif
      enddo
      end subroutine

$endif




C =================================================================
C Process land influx sources, set water influx to surface layer:

      SUBROUTINE QWCALC( NBI, INDXI, NLI, PRECIP, EVAP, AREA,
     &            QWSURF ) !=
C In:
      INTEGER NBI, INDXI(NBI+1), NLI
      real*8 PRECIP, EVAP(NBI), AREA(NLI)
C Out:
      real*8 QWSURF(NBI)    ! Net water influx from land/atmosphere (m3/s)

C -------------- local variables:
      INTEGER IB, IL

C -------- Calculate total net water influx to surface layer
C          from terms for atmospheric interaction pr. area:
      DO IB = 1,NBI
          IL = INDXI(IB)+1
          QWSURF (IB) = (PRECIP-EVAP(IB))*AREA(IL)
      END DO
      END SUBROUTINE


C ===============================================================
C Add effect of land influx to biochemical conc. time derivatives
       SUBROUTINE QCALC( NUM_OF_LAYERS, RCV_Q_INDEX,
     &          C_DEV, Q_SUBST, W_FACTOR, NAME, C_IMP)


C                                result of jet calculations in TRANSP.
!       INCLUDE    'EUTRO.INC'


      INTEGER NUM_OF_LAYERS (*)
      INTEGER RCV_Q_INDEX (NS) ! Global layer index of receiving layer
      real*8 C_DEV(NLI)  ! Conc. derivative in water Cmass/m3/day
      real*8 C_IMP(NBI)  ! Mass import Cmass/day
      real*8 Q_SUBST(NS) ! Qmass/d
      real*4 W_FACTOR    ! factor going from Qmass to Cmass
C                        (stochiometry and/or effectivity)
      CHARACTER*8 NAME


C  IMPORTS ACSL COMMON BLOCK for NS, BASINQ, VLAYER
      INTEGER I_S, I_B, I_L, L_COUNT
      real*8 C_INPUT
      
      logical WARNING_given /.false./

$if defined  DEBUG_QCALC
!      INCLUDE   'DEBUG.INC'
      LOGICAL DEBUG_ON

!        integer*4 TEST_NAME /'C1  '/
!        DEBUG_ON = (NAME.eq.TEST_NAME.AND. T.GE.TTRIG)
      
      DEBUG_ON = T.GE.TTRIG
      if ( DEBUG_ON ) then
         write(debug_unit,'(A,A8,1x,A,I3)')
     &          ' QCALC FOR VARIABLE ',NAME, ',  NS:',NS
         write(debug_unit,'(1x,2A15)')
     &          ' NUM_OF_LAYERS:', ' RCV_Q_INDEX:'
         write(debug_unit,'(2(1x,I10,5x))')
     &           ( NUM_OF_LAYERS(I_S), RCV_Q_INDEX(I_S),I_S=1,NS )
      endif
$endif

      DO I_S = 1,NS
         I_B = BASINQ(I_S)
         I_L = RCV_Q_INDEX(I_S)

C Input to each receiving layer, distributed below:
         C_INPUT = Q_SUBST(I_S)*W_FACTOR / FLOAT(NUM_OF_LAYERS(I_S) )

$if defined  DEBUG_QCALC
         if (DEBUG_ON) then
             write(debug_unit,*)
     &             'Runoff',I_S,' to basin ',I_B,':'
             write(debug_unit,'(2(1X,a,g15.8))')
     &           'C_IMP(I_B) before:'     , C_IMP(I_B)
     &         , 'C_INPUT to each layer:' , C_INPUT
         endif
$endif

C Total input of substance:
         C_IMP(I_B) = C_IMP(I_B) + Q_SUBST(I_S)*W_FACTOR

$if defined  DEBUG_QCALC
      if (DEBUG_ON) then
          write(debug_unit,'(1X,a,g15.8)')
     &           'C_IMP(I_B) after:',C_IMP(I_B)
      endif
$endif

         warning_Given =.false.
         DO L_COUNT = 1, NUM_OF_LAYERS(I_S)


            if ( I_L .gt. INDXI(I_B) .and. I_L .le. INDXI(I_B+1) ) then

$if defined  DEBUG_QCALC
               if (DEBUG_ON) then
                  write(debug_unit,*)
                  write(debug_unit,'(A,I3,2(1X,a,g15.8))')
     &                 ' Layer ',I_L
     &               , 'before: C_DEV=', C_DEV(I_L)
     &               , 'VLAYER=',vlayer(i_l)
               endif
$endif
               C_DEV(I_L) =  C_DEV(I_L) + C_INPUT/VLAYER(I_L)

$if defined  DEBUG_QCALC
               if (DEBUG_ON) then
                   write(debug_unit,'(1X,a,g15.8)')
     &                  ' after: C_DEV=', C_DEV(I_L)
               endif
$endif
               
            Else
               if (.not. warning_given) then
                  write(*,*) 'EUTROSUB/QCALC: Warning - runoff',I_S,
     &              ' to basin ',I_B,
     &              'into illegal layer number ',I_L
               endif
               EXIT

            Endif
            I_L = I_L-1
         ENDDO
      ENDDO
      END SUBROUTiNE


C ===============================================================
C Add effect of surface influx on concentration time derivatives
C         Only used on unity concentration derivative
C         to account for effect of water exchange through surface
      SUBROUTINE QSCALC( C_DEV, C_IMP, Q_SUBST, W_FACTOR  )

C  IMPORTS ACSL COMMON BLOCK for NLI, NBI, INDXI, VLAYER

!       INCLUDE    'EUTRO.INC'


      real*8 C_DEV(NLI)   ! Cmass/M3
      real*8 C_IMP(NBI)  ! Cmass/day
      real*8 Q_SUBST(NBI) ! SUM (Qmass/time unit) to surface layers
      real*4 W_FACTOR     ! factor going from Q to C
C                          (stochiometry and/or effectivity)

      INTEGER I_B, I_L
      real*8 C_INPUT

      DO I_B = 1,NBI
         I_L = INDXI(I_B)+1
         C_INPUT = Q_SUBST(I_B)*W_FACTOR
         C_DEV(I_L) =  C_DEV(I_L) + C_INPUT/VLAYER(I_L)
         C_IMP(I_B) =  C_IMP(I_B) + C_INPUT
      ENDDO
      END SUBROUTiNE


C ===================================================================
C Consolidate import rate values to use in mass balance:
      SUBROUTINE SUMIMP

!       INCLUDE    'EUTRO.INC'

$if defined  DEBUG_SUMIMP
!      INCLUDE   'DEBUG.INC'
      real*8 OXYG_IMPORT_BEFORE, NITRMP_BEFORE
$endif

      INTEGER I_B


      CALL HELLO( 'SUMIMP'   )

      DO I_B = 1, NBI
C ....... oxygen-carbon balance:
C           Organic carbon as -oxygen consumed in oxic degradation
C           NO3 as + oxygen released in ammonification/uptake
C           Denitrification as oxygen export, because nitrate has lower
C           oxygen( -carbon) equivalent here than in ammonification,
C           the difference is counted as oxygen loss:
$if defined  DEBUG_SUMIMP
         OXYG_IMPORT_BEFORE = OXYGMP(I_B)
         NITRMP_BEFORE = NITRMP(I_B) ! NITRATE FROM LAND RUNOFF
$endif

         OXYGMP(I_B) =  OXYGMP(I_B) + ASEDXP(I_B) 
     &        - ODMImp(I_B) -ODMLOAD(I_B)
     &        - ( CFYTMP(I_B) +DOCMP(I_B) + CDETMP(I_B) + CLOAD(I_B)
     &           +CZOOMP(I_B) + BACTMP(I_B) - CSEDXP(I_B) )*OX_C
     &        + ( NO3MP (I_B) + NITRMP(I_B ) )*OX_NITR
     &        - DENITR(I_B)*(OX_NITR - OX_C/DENITR_C)

C ....... combine import of nitrogen in different forms:
C         (NLOAD CONTAINED NH4&PARTICULATE RUNOFF, NOW CHANGED TO TOTAL
c          LAND RUNOFF OF N):
         NLOAD(I_B) =  NITRMP(I_B) + NLOAD(I_B) ! Total from land
         NITRMP(I_B) = NLOAD(I_B) + NFIX(I_B)
     &        + NO3MP (I_B) + NH4MP(I_B)  - DENITR(I_B)
     &        + NFYTMP(I_B) + NDETMP(I_B) - NSEDXP(I_B)
     &        + CZOOMP(I_B) * NCZOO + BACTMP(I_B) * NCBACT

C ....... combine import of phosphorus in different forms:
         PHOSMP(I_B) = PLOAD(I_B)
     &        + PO4MP(I_B) + PFYTMP(I_B) + PDETMP(I_B)
     &        - PSEDXP(I_B) + CZOOMP(I_B)*PCZOO + BACTMP(I_B)*PCBACT


C       Zooplankton carbon incorporates N and P in Redfield ratios,
C       with N and P removed from free nutrient pool.  Therefore
C       zooplankton carbon enters the N and P mass balances.

C       Bacterial carbon are connected to N and P within the
C       nutrient pool, and does not affect N and P mass balance.
C             cfr. PHYT_ZOO.FOR for details.

C       The mass balances do not consider different release rates
C       for different components.

C ....... combine import of silicate in different forms:
         SILIMP(I_B) =  SIO2MP(I_B) + SFYTMP(I_B)
     &                 + SDETMP(I_B) - SSEDXP(I_B)

$if defined  DEBUG_SUMIMP
         IF ( t.ge.TTRIG ) THEN
            WRITE( DEBUG_UNIT,
     &          '('' >>>> SUMIMP, net import to basin:'',I5,'':'')') I_B
            WRITE( DEBUG_UNIT, '(3(1x,A10,'':''G13.7))')
     &        'SALTMP', SALTMP (I_B) ,
     &        'HEATMP', HEATMP (I_B) ,
     &        'NITRMP_BEFORE', NITRMP_BEFORE,
     &        'NLOAD ', NLOAD (I_B)  ,
     &        'NFIX  ', NFIX  (I_B)  ,
     &        'NO3MP ', NO3MP (I_B)  ,
     &        'NH4MP ', NH4MP (I_B)  ,
     &        'NFYTMP', NFYTMP(I_B)  ,
     &        'NDETMP', NDETMP(I_B)  ,
     &        'NSEDXP', NSEDXP(I_B)  ,
     &        'DENITR', DENITR(I_B)  ,
     &        '*** NITRMP', NITRMP(I_B)  ,
     &        'PLOAD ', PLOAD (I_B)  ,
     &        'PO4MP ', PO4MP (I_B)  ,
     &        'PFYTMP', PFYTMP(I_B)  ,
     &        'PDETMP', PDETMP(I_B)  ,
     &        'PSEDXP', PSEDXP(I_B)  ,
     &        '*** PHOSMP', PHOSMP(I_B)  ,
     &        'CFYTMP', CFYTMP(I_B)  ,
     &        'CDETMP', CDETMP(I_B)  ,
     &        'CSEDXP', CSEDXP(I_B)  ,
     &        'CZOOMP', CZOOMP(I_B)  ,
     &        'ODMImp  ', ODMImp  (I_B)  ,
     &        'DOCMP ', DOCMP (I_B)  ,
     &        'BACTMP', BACTMP(I_B)  ,
     &        'ASEDXP', ASEDXP(I_B)  ,
     &        'OXYGMP', OXYG_IMPORT_BEFORE,
     &        '--->  ', OXYGMP(I_B)  ,
     &        'SIO2MP', SIO2MP(I_B)  ,
     &        'SFYTMP', SFYTMP(I_B)  ,
     &        'SDETMP', SDETMP(I_B)  ,
     &        'SSEDXP', SSEDXP(I_B)  ,
     &        '*** SILIMP', SILIMP(I_B)
         ENDIF
$endif
      ENDDO
      END SUBROUTiNE



C ===============================================================
C Calculate control values for volume conservation:
      SUBROUTINE VOLCHK( CONTEXT_NUM, NBI, NC, INDXC, BCONN1, BCONN2,
     &                   NLC, VBUF, VDYN, DVTOT, VPRT, vdindx,
     &        VBFSUM, DVTOT2, VTDIFF )

C In:
      integer CONTEXT_NUM
           ! = 1: derivatives
           ! = 2: values
      INTEGER NBI
      INTEGER NC, INDXC(NC+1), BCONN1(NC), BCONN2(NC)
      INTEGER NLC
      real*8 VBUF(2, NLC), VDYN(2, NBI), DVTOT( NBI)
      LOGICAL VPRT
      integer vdindx
C Out:
      real*8 VBFSUM(NBI+1)
C     = Net sum of "displaced" volume for each basin:
C     = + ä volume of "internal" water  in buffers outside connections
C       - ä volume of water from other basins inside connections.
C     Last value is sum for external basins.
      real*8 DVTOT2 (NBI+1)
C     = Dynamic volume + displacement, should equal integrated DVTOT
      real*8 VTDIFF (NBI+1)
C     = Difference between dvtot2 and dvtot, should equal zero



      INTEGER I_C, IB, IB_2, I_LC
      real*8 VBUF_EXT, VDYN_EXT, DVTOT_EXT, NET_DISPLACED

$if defined  DEBUG_VOLCHK
!      INCLUDE   'DEBUG.INC'
      CHARACTER CONTEXT(2)*12/'DERIVATIVES','VALUES'/
$endif

      VBUF_EXT = 0.0
      VDYN_EXT = 0.0
      DVTOT_EXT = 0.0

      DO IB = 1,NBI
          VBFSUM(IB) = 0.0
          VDYN_EXT = VDYN_EXT - VDYN(vdindx,IB)
          DVTOT_EXT = DVTOT_EXT - DVTOT(IB)
      ENDDO

$if defined  DEBUG_VOLCHK
      IF (VPRT) THEN
          WRITE(DEBUG_UNIT,*)' >>>>>>>>>>>>>>>>>>>>>>>>>>>'
          WRITE (debug_unit,'('' ----VOLCHK FOR '',A,''   NC='',I3)')
     &           CONTEXT(CONTEXT_NUM), NC
      ENDIF
$endif

      DO I_C = 1, NC
          IB = BCONN1(I_C)
          IB_2 = BCONN2(I_C)
          DO I_LC = INDXC(I_C)+1, INDXC(I_C+1)

$if defined  DEBUG_VOLCHK
             IF (VPRT) THEN
                WRITE (DEBUG_UNIT,'('' I_C, IB, IB_2, I_LC: '',4I6)' )
     &             I_C, IB, IB_2, I_LC
C            ...... BASIN IB ALWAYS INTERNAL:
                WRITE (DEBUG_UNIT,
     &               '('' VBUF(1,I_LC),VBUF(2,I_LC)'',2G15.7)')
     &            VBUF(1,I_LC), VBUF(2,I_LC)
             ENDIF
$endif

             NET_DISPLACED =   - VBUF(1,I_LC) + VBUF(2,I_LC)

$if defined  DEBUG_VOLCHK
             IF (VPRT) WRITE(DEBUG_UNIT,'('' NET_DISPLACED'',G15.7)')
     &            NET_DISPLACED
$endif

             VBFSUM(IB) = VBFSUM(IB) + NET_DISPLACED
             IF ( IB_2 .GT. 0 ) THEN
                VBFSUM(IB_2) = VBFSUM(IB_2) - NET_DISPLACED

$if defined  DEBUG_VOLCHK
               IF (VPRT) WRITE(DEBUG_UNIT,
     &                '('' VBFSUM(IB), VBFSUM(IB_2)'',3G15.7)')
     &             VBFSUM(IB), VBFSUM(IB_2)
$endif

             ELSE
                VBUF_EXT = VBUF_EXT - NET_DISPLACED

$if defined  DEBUG_VOLCHK
                IF (VPRT) WRITE(DEBUG_UNIT,
     &                '('' VBFSUM(IB), VBUFEXT'',3G15.7)')
     &            VBFSUM(IB), VBUF_EXT
$endif

            ENDIF
          ENDDO
      ENDDO
      VBFSUM(NBI+1) = VBUF_EXT

$if defined  DEBUG_VOLCHK
      IF (VPRT) then
          WRITE(DEBUG_UNIT,*) ' VBUF  :',VBUF
          WRITE(DEBUG_UNIT,*) ' VBFSUM:', VBFSUM
          WRITE(DEBUG_UNIT,'(4X,4a18)')
     &        'DVTOT:','VDYN(VDINDX,..):','DVTOT2:','VTDIFF:'
      ENDIF
$endif

      DO IB = 1,NBI
          DVTOT2(IB) = VDYN(vdindx,IB) + VBFSUM(IB)
          VTDIFF(IB) = DVTOT2(IB)-DVTOT(IB)  ! should be zero
$if defined  DEBUG_VOLCHK
         IF (VPRT) then
            WRITE(DEBUG_UNIT,'(1X,4g18.8)')
     &           DVTOT(IB),VDYN(VDINDX,IB),DVTOT2(IB),VTDIFF(IB)
         ENDIF
$endif
      ENDDO


      DVTOT2(NBI+1) = VDYN_EXT  + VBUF_EXT
      VTDIFF(NBI+1) = DVTOT2(NBI+1) - DVTOT_EXT  ! should be zero

$if defined  DEBUG_VOLCHK
         IF (VPRT) WRITE(DEBUG_UNIT,*)' exits VOLCHK'
$endif
      
      END SUBROUTINE


$if defined  SPECIAL_CHECK
      subroutine check(N,PO4,T)
      integer N
      real*8 PO4(N), T

      real*8 po4_mem(39,10),T_mem(10)
      integer Last /0/
      integer Count/0/
      logical negative /.false./
      integer I,K,L

      Last = Mod(LAST,10)+1

      Count=min(10,Count+1)
      negative =.false.

      do I=1,N
         PO4_mem(I,Last)=PO4(I)
         if (PO4(I).lt.0.0) negative=.true.
      enddo

      T_mem(Last)=T

      if (negative) then
         do K = Count,1,-1
            L = Mod(Last-K+10,10)+1
            if(K.eq.Count)
     &      write(*,'('' At T='',G17.11, '': PO4<0 (see unit 999'')')
     &                  T_mem(L)
            write(999,'('' At T='',G17.11, '': PO4='')')  T_mem(L)
            write(999,'(1x,5G12.5)') (PO4_mem(I,L),I=1,N)
         enddo
!         Pause
      Endif

      END SUBROUTiNE
$endif



!=================================================================
      Subroutine Integrate_Physical_Variables(TIME_STEP)
      real*8 TIME_STEP

				!  ---------- Integrate Physical state variables ----------

		DVTOT  = DVTOT + VTOTDV*TIME_STEP
		VDYN   = VDYN  + VDYNDV*TIME_STEP
		VBUF   = VBUF  + VBUFDV*TIME_STEP
      VBFSI  = VBFSI + VBFSDV*TIME_STEP
      DVT2I  = DVT2I + DVT2DV*TIME_STEP
      VTDI   = VTDI  + VTDDV *TIME_STEP

				! integrated time (remnant from old ACSL code, not needed now?)
		SUMTIM = SUMTIM + TIME_STEP   

				! Calculate accumulated wind mix energy/area 
				! and Wind frictional velocity in 3rd exponent
		FR3INT = FR3INT + UFRIC3*TIME_STEP

				!  ------- Buoyancy energy balance:
				! Accumulated net buoyant influx effect: (m2/s3*day)
		BFXINT = BFXINT + BFXZ*TIME_STEP

				!  -------------- INTEGRATE CLIMATE VARIABLES -------------
				! RADIATION BUDGET TERMS:
		QOUTI = QOUTI + QOUT*TIME_STEP
		QABSI = QABSI + QABS*TIME_STEP
		QGI   = QGI   + (SUNRAD+DIFRAD)*TIME_STEP

				! AIR TEMPERATURE AND WIND TERMS
		AIRTI = AIRTI + AIRTMP*TIME_STEP
		WN2I  = WN2I  + (WINDN**2)*TIME_STEP
		WE2I  = WE2I  + (WINDE**2)*TIME_STEP

      end subroutine


      

!==================================================================
! Subroutine for setting all elements in array to a specified value.
! Used for historical reasons (original ACSL Metafortran code)
! Could be replaced by simple Aray assignment      
      
      subroutine SETR (SpecifiedValue, Array, N)
      real*4 SpecifiedValue
      integer N
      real*8 Array(1:N)
      Array(1:N) = SpecifiedValue
      end subroutine

      End Module
