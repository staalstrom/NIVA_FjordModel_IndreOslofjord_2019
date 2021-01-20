      Module fx_MBALAN


####### Må gjennomgås grundig - endring fra entrypoints til subrutiner


      use ModelParam_RunControl
      use ModelParam_Plankton
      use ModelParam_Mussels, only: NCMUSL, PCMUSL
      use ModelVar_RunControl
      use ModelVar_HydroBioChem, only: VLAYER, FYTGRP
      use ModelVar_Topography, only: BOTTOM, INDXI, NBI
      
      use fx_RunControl, only: DEBUG_UNIT
      use fx_Stoichiometry
      
      implicit none
      
      integer, private :: SUM_mode
      real*8, allocatable, private :: STORED_SUM (:)
      integer, private :: IS_ALLOCATED = 0

      contains

      
C ====================================================================
C Calculate balance of mass and heat substance:
      SUBROUTINE MBALAN ( MBI, MLI, SIZE_SCALE,
     &        SAL_C , TEMP_C, OXYG_C, NO3_C,  NH4_C, PO4_C, SiO2_C,
     &        CFYT_C, NFYT_C, PFYT_C, SFYT_C,
     &        CDET_C, NDET_C, PDET_C, SDET_C,
     &        CSED_C, NSED_C, PSED_C, SSED_C, ASED_C, PADS_C,
     &        ODM_C , DOC_C , BACT_C, CZOO_C, CMUSL_C,
     &        SPP_C, SPPSED_C,
     &     SALT_M, HEAT_M, OXYG_M, NITR_M, PHOS_M, SiO2_M, SPP_M ) !=


!      INCLUDE 'EUTRO.INC'  !  ACSL Common block

C Input values: volume and concentrations:
      INTEGER MBI, MLI
      LOGICAL SIZE_SCALE
      real*8  SAL_C(MLI), TEMP_C(MLI), OXYG_C(MLI)
      real*8  NO3_C(MLI), NH4_C(MLI), PO4_C(MLI), SiO2_C(MLI)

      real*8  CFYT_C(MLI,*), NFYT_C(MLI,*)
      real*8  PFYT_C (MLI,*), SFYT_C(MLI)

      real*8  CDET_C (MLI), NDET_C(MLI)
      real*8  PDET_C (MLI), SDET_C(MLI)

      real*8  CSED_C (MLI), NSED_C(MLI)
      real*8  PSED_C (MLI), SSED_C(MLI)
      real*8  ASED_C (MLI), PADS_C(MLI)
      real*8  ODM_C  (MLI)
      real*8  DOC_C  (MLI)
      real*8  BACT_C (MLI)
      real*8  CZOO_C (MLI)
      real*8  CMUSL_C(MBI)

      real*8  SPP_C(MLI), SPPSED_C(MLI)

C Output: amount summed over volumes:
      real*8  SALT_M(MBI), HEAT_M(MBI), OXYG_M(MBI)
      real*8  NITR_M(MBI), PHOS_M(MBI), SiO2_M(MBI), SPP_M(MBI)

C       Conversion factors: OX_C, OX_NITR

C       Identifying names:
      CHARACTER*5
     &  N_NITR /'NITR' /,
     &  N_PHOS /'PHOS' /,
     &  N_SILI /'Si  ' /,
     &  N_OXYG /'OXYG' /,
     &  N_NO3  /'NO3 ' /,
     &  N_NH4  /'NH4 ' /,
     &  N_PO4  /'PO4 ' /,
     &  N_CZOO /'CZOO' /,
     &  N_CMUSL/'CMUSL'/,
     &  N_BACT /'BACT' /,
     &  N_SPP  /'SPP'  /
      CHARACTER NAME*4

      integer N_F

!      write(*,"(A,i12,A,z12)") ' i MBALAN: MBI,',MBI,
!     &     ' ADRESSE: ', LOCFAR(MBI)
      CALL ACCUM_INIT (MBI, SIZE_SCALE)

C .......... Physical properties:
      CALL ACCUM( MBI, MLI,  1, VLAYER, 'SAL',  SAL_C,  'SALT',
     &            0.0, 1.0, SALT_M )
      CALL ACCUM( MBI, MLI,  2, VLAYER, 'TEMP', TEMP_C, 'HEAT',
     &            0.0, 1.0, HEAT_M )
C .......... Inorganic forms:
      CALL ACCUM( MBI, MLI,  3, VLAYER, 'OXYG', OXYG_C, N_OXYG,
     &            0.0, 1.0, OXYG_M )
      CALL ACCUM( MBI, MLI,  43, VLAYER, N_NO3, NO3_C,  N_NITR,
     &            0.0, 1.0, NITR_M )
           CALL ACCUM2( MBI,  3, N_OXYG, 1.0, OX_NITR, OXYG_M )
      CALL ACCUM( MBI, MLI,  4, VLAYER, N_NH4,  NH4_C,  N_NITR,
     &            1.0, 1.0, NITR_M )
      CALL ACCUM( MBI, MLI,  5, VLAYER, N_PO4,  PO4_C,  N_PHOS,
     &            0.0, 1.0, PHOS_M )
      CALL ACCUM( MBI, MLI,  6, VLAYER, 'SiO2', SiO2_C, N_SILI,
     &            0.0, 1.0, SiO2_M )


C ........... OXYGEN DEMAND:
      CALL ACCUM( MBI, MLI,  3, VLAYER, 'ODM', ODM_C, N_OXYG,
     &            1.0, -1.0, OXYG_M )

C ......... phytoplankton, FYTGRP fractions:
      DO N_F = 1,FYTGRP
        NAME = 'FYT'//CHAR(N_F+48)
C           organic carbon as potential oxygen consumption:
        CALL ACCUM( MBI, MLI, 3, VLAYER, 'C'//NAME, CFYT_C(1,N_F),
     &               N_OXYG,  1.0, -OX_C, OXYG_M )
        CALL ACCUM( MBI, MLI, 4, VLAYER, 'N'//NAME, NFYT_C(1,N_F),
     &               N_NITR,  1.0, 1.0,   NITR_M )
        CALL ACCUM( MBI, MLI, 5, VLAYER, 'P'//NAME, PFYT_C(1,N_F),
     &               N_PHOS,  1.0, 1.0,   PHOS_M )
      ENDDO
      CALL ACCUM(MBI, MLI,  6,VLAYER, 'SFYT1', SFYT_C,
     &             N_SILI, 1.0, 1.0, SiO2_M )

C ......... dissolved organic carbon as potential oxygen consumption:
      CALL ACCUM(MBI, MLI,  3, VLAYER, 'DOC' ,DOC_C, N_OXYG,
     &             1.0, -OX_C, OXYG_M )

C ......... Bacteria and zooplankton, assumed fixed composition:
      CALL ACCUM( MBI, MLI,  345, VLAYER,
     &                    N_BACT, BACT_C, N_OXYG, 1.0, -OX_C,  OXYG_M )
          CALL ACCUM2( MBI, 4, N_NITR, 1.0, NCBACT, NITR_M )
          CALL ACCUM2( MBI, 5, N_PHOS, 1.0, PCBACT, PHOS_M )
      CALL ACCUM( MBI, MLI,  345, VLAYER,
     &                    N_CZOO, CZOO_C, N_OXYG, 1.0, -OX_C,  OXYG_M )
          CALL ACCUM2( MBI, 4, N_NITR, 1.0, NCZOO,  NITR_M )
          CALL ACCUM2( MBI, 5, N_PHOS, 1.0, PCZOO,  PHOS_M )

C ......... Mussels: carbon summed over basin, signalled by 3.arg <0,
C                    the weights in 4. argument: are not used.
      CALL ACCUM( MBI, MLI,  -345, VLAYER,
     &                  N_CMUSL, CMUSL_C, N_OXYG, 1.0, -OX_C,  OXYG_M )
          CALL ACCUM2( MBI, 4, N_NITR, 1.0, NCMUSL, NITR_M )
          CALL ACCUM2( MBI, 5, N_PHOS, 1.0, PCMUSL, PHOS_M )

C ......... Detritus, C,N,P AND S + particles:
      CALL ACCUM( MBI, MLI,  3, VLAYER, 'CDET', CDET_C,
     &               N_OXYG, 1.0, -OX_C, OXYG_M )
      CALL ACCUM( MBI, MLI,  4, VLAYER, 'NDET', NDET_C,
     &               N_NITR, 1.0, 1.0,   NITR_M )
      CALL ACCUM( MBI, MLI,  5, VLAYER, 'PDET', PDET_C,
     &               N_PHOS, 1.0, 1.0,   PHOS_M )
      CALL ACCUM( MBI, MLI,  6, VLAYER, 'SDET', SDET_C,
     &               N_SILI, 1.0, 1.0,   SiO2_M )
      CALL ACCUM( MBI, MLI,  7, VLAYER, 'SPP', SPP_C, 
     &               N_SPP,  0.0, 1.0, SPP_M )


C ......... Sedimented material, C,N,P AND S + particles:
      CALL ACCUM( MBI, MLI,  3, BOTTOM, 'CSED', CSED_C,
     &               N_OXYG, 1.0, -OX_C, OXYG_M )
      CALL ACCUM( MBI, MLI,  4, BOTTOM, 'NSED', NSED_C,
     &               N_NITR, 1.0, 1.0,   NITR_M )
      CALL ACCUM( MBI, MLI,  5, BOTTOM, 'PSED', PSED_C,
     &               N_PHOS, 1.0, 1.0,   PHOS_M )
      CALL ACCUM( MBI, MLI,  6, BOTTOM, 'SSED', SSED_C,
     &               N_SILI, 1.0, 1.0,   SiO2_M )
      CALL ACCUM( MBI, MLI,  7, BOTTOM, 'SPP', SPPSED_C,
     &               N_SPP,  1.0, 1.0,   SPP_M )


      CALL ACCUM( MBI, MLI,  3, BOTTOM, 'ASED', ASED_C,
     &             N_OXYG, 1.0, -1.0, OXYG_M )
                        ! (Oxygen debt as positive values)
      CALL ACCUM( MBI, MLI,  5, BOTTOM, 'PADS', PADS_C,
     &             N_PHOS, 1.0, 1.0,   PHOS_M )



      END SUBROUTINE

C ===============================================================
C Accumulate variable over layers, use volumes or area:
      SUBROUTINE ACCUM_INIT ( MBI, SIZE_SCALE )


      INTEGER MBI
      LOGICAL SIZE_SCALE ! =.true. :will sum absolute values


 !     INCLUDE 'EUTRO.INC'  !  ACSL Common block

 !     INCLUDE 'DEBUG.CMN'



!     write(*,*) ' ALLOCATED, MBI ', ALLOCATED, MBI
!     write(*,"(A,z12)") ' ADDRESSE MBI: ', LOCFAR(MBI)
      if ( IS_ALLOCATED .lt. MBI ) then
          if (IS_ALLOCATED.gt.0)  DEALLOCATE ( STORED_SUM )
          ALLOCATE  ( STORED_SUM(MBI) )
          IS_ALLOCATED = MBI
      endif

      if (SIZE_SCALE) then
          SUM_MODE = 1
      else
          SUM_MODE = 2
      endif

      End Subroutine ACCUM_INIT



C ----------------------------------------------------------------
      subroutine ACCUM (MBI, MLI, MBPRT_X, WEIGHT, VARIABLE_NAME,
     &      CONS_C,COMPONENT_NAME, FACTOR_OLD, FACTOR_NEW, SUM_BASIN)

      INTEGER MBPRT_X ! Index to MBPRT - controlling debug printout
      INTEGER MBI,MLI

      real*8 WEIGHT (MLI)           ! Weighting factor (volume or area)
      CHARACTER*(*) VARIABLE_NAME, COMPONENT_NAME
      real*8 CONS_C (MLI)           ! Concentration per volume or area

C         (if M<0: total amount in basin)

      real*4 FACTOR_OLD, FACTOR_NEW  ! Controls init/update.

C  Result exported: initiated or updated amount per basin:

      real*8 SUM_BASIN(MBI)


      integer test_print
      INTEGER M_test, TEST_INDEX
      INTEGER I_B,L_X
      real*8 ACCUM_SUM
      character name_prefix(2)*4  /'Abs_',' '/

      test_print=0
      IF (T.GE.TTRIG) THEN
         M_test = ABS(MBPRT_X)
         DO WHILE ( M_test .gt. 0 )
            TEST_INDEX = MIN(6,MAX(1,MOD( M_test,10)))
                         ! ensure that TEST_INDEX is within range og MBPRT
            test_print = MAX (test_print, MBPRT(TEST_INDEX) )
            M_test = M_test/10
         ENDDO  ! Last value of TEST_INDEX is first digit in M_TEST
      ENDIF

      DO I_B = 1,NBI
         IF ( test_print.gt.0 ) THEN
           WRITE( DEBUG_UNIT,
     &            '(/'' +++++ Sum of '',2A,'', basin '',I5 )')
     &           Name_Prefix(Sum_Mode),VARIABLE_NAME, I_B
         Endif
         if (MBPRT_X.gt.0) then
            IF ( test_print.gt.1 ) THEN
               WRITE(DEBUG_UNIT,'(1x,A5,3A14)')
     &         'Layer','Value','Volume(area)','Cumulative sum'
            ENDIF
            ACCUM_SUM = 0.0
            if ( SUM_MODE.eq.2 ) then
               DO L_X = INDXI(I_B)+1, INDXI(I_B+1)
                  ACCUM_SUM = ACCUM_SUM + CONS_C(L_X) * WEIGHT(L_X)
                  IF ( test_print.gt.1 ) THEN
                      WRITE(DEBUG_UNIT,'(1x,I5,3G14.7)' )
     &                   L_X, CONS_C(L_X), WEIGHT(L_X), ACCUM_SUM
                  ENDIF
               ENDDO
            else  ! Sums absolute values for accuracy scale:
               DO L_X = INDXI(I_B)+1, INDXI(I_B+1)
                  ACCUM_SUM = ACCUM_SUM + ABS(CONS_C(L_X)*WEIGHT(L_X))
                  IF ( test_print.gt.1 ) THEN
                      WRITE(DEBUG_UNIT,'(1x,I5,3G14.7)' )
     &                   L_X, CONS_C(L_X), WEIGHT(L_X), ACCUM_SUM
                  ENDIF
               ENDDO
            ENDIF
         else ! Negative MBPRT_X signals that total amount
              ! in basin is entered, and Weight argument not used.
            ACCUM_SUM = CONS_C(I_B)
         Endif
         IF (test_print.eq.1 .or. (test_print.ge.1.and. MBPRT_X.le.0))
     &         WRITE(DEBUG_UNIT,'(1x,A,3G14.7)' )
     &        'Cumulative sum: ',  ACCUM_SUM
         STORED_SUM (I_B) = ACCUM_SUM
      ENDDO

            ! Update primary component sum
            ! Test printout is controlled by 
            ! TEST_INDEX = last digit of MBPRT_X


      CALL ACCUM2 ( MBI, TEST_INDEX, COMPONENT_NAME,
     &              FACTOR_OLD, FACTOR_NEW,
     &              SUM_BASIN ) !=

            ! Contribution to other component sums 
            ! are handled by separate calls to ACCUM2 from MBALAN
      
      END SUBROUTINE ACCUM

C ===============================================================
C Update component sum according to given factors:
      SUBROUTINE ACCUM2 ( MBI, MBPRT_X, COMPONENT_NAME,
     &               FACTOR_OLD, FACTOR_NEW,
     &              SUM_BASIN ) !=

C ----------------- arguments for additional entry points:
      INTEGER MBPRT_X ! Index to MBPRT - controlling debug printout
      INTEGER MBI

      CHARACTER*(*) COMPONENT_NAME
      real*4 FACTOR_OLD
      real*4 FACTOR_NEW  ! Controls init/update.

C  Result exported: initiated or updated amount per basin:
      real*8 SUM_BASIN(MBI)

      Character*80 PRINT_STRING

      integer TEST_INDEX, test_print
      integer I_B, I, L
      real*8 OLD_SUM

C for additional direct calls from MBALN to ACCUM2
C MBPRT_X controls printout:

      TEST_INDEX = Min(6,Max(1,MBPRT_X))
      
      IF (T.GE. TTRIG) THEN
         test_print = MBPRT(TEST_INDEX)
      ELSE
         test_print = 0
      ENDIF

      DO I_B = 1,NBI
         OLD_SUM = SUM_BASIN(I_B)
         if (Sum_Mode.eq.1) then
             SUM_BASIN(I_B) =   ABS(FACTOR_OLD*OLD_SUM)
     &                    + ABS(FACTOR_NEW*STORED_SUM (I_B) )
         else
             SUM_BASIN(I_B) =   FACTOR_OLD*OLD_SUM
     &                    + FACTOR_NEW*STORED_SUM (I_B)
         endif
         IF ( test_print.gt.0 ) THEN
            IF ( FACTOR_OLD .eq. 0.0 ) THEN
               WRITE(DEBUG_UNIT,'(10x,''initiates '',A,'' basin '',I5)')
     &            COMPONENT_NAME, I_B
               WRITE(PRINT_STRING, '(1x,3(G13.7:A))')
     &            STORED_SUM(I_B), '*', FACTOR_NEW, '=', SUM_BASIN(I_B)
            ELSE
               WRITE(DEBUG_UNIT,'(10x,''added to '',A,'' basin '',I5)')
     &            COMPONENT_NAME, I_B
               WRITE( PRINT_STRING, '(1x,5(G13.7:A))')
     &            OLD_SUM, '*', FACTOR_OLD, '+',
     &            STORED_SUM(I_B), '*', FACTOR_NEW, '=', SUM_BASIN(I_B)
            ENDIF
            L = 0  ! Pack and print string:
            DO I = 1,LEN_TRIM(PRINT_STRING)
               if ( PRINT_STRING(I:I) .ne. ' ') then
                   L = L+1
                   PRINT_STRING(L:L) = PRINT_STRING(I:I)
                endif
            ENDDO
            WRITE( DEBUG_UNIT, '(1x,A)') PRINT_STRING(1:L)
         ENDIF
      ENDDO
      END Subroutine ACCUM2
      
      End Module
