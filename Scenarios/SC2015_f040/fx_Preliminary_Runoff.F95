      Module fx_Runoff

      use m6_FileNames
      use fx_RunControl
      use fx_RandFnct
      use fx_OxygSat

      implicit none

C ==================================================================
C Eutrophication model  - File  RNF_ALL.FOR
C Contains subroutines for handling influx from land sources.
C =================================================================

! versjon fra Milj_Bfj: endret Include-setning for bruk i ver2007



$unDEFINE DEBUG_INFLUX

$IF DEFINED DEBUG_INFLUX
$UNDEFINE DEBUG_IDDEFJ_Q
$ENDIF

$unDEFINE CHECK_Q_TRANSFER


      CONTAINS


C =================================================================
C Influx from land:

      SUBROUTINE RUNOFF ( DEBUG, T, NBI, INDXI, ND, DEPTH, AIRTMP,
     &     MS, BASINQ, AIRTQF, QFW, QFP, QFN, QFODM, QFC, QFS, 
     &   QFNH4, QCDETF, QNDETF, QPDETF, QNCDET, QPCDET,
     &   QTRF, QTRNDX, QTR_Times, DEPTHQ, QMIXIN, MIXQTM, QTR_ACTIVE,
     &   RNDFACTOR, RNFNDX, QWATER, QTEMP, QPO4, QNO3, QNH4,
     &   QODM, QDOC, QCDET, QNDET, QPDET, QSi, QO2 )
      

C -----------------------------------------------------------------
C In:
      LOGICAL DEBUG
      real*8 T
      INTEGER NBI
      INTEGER INDXI ( NBI+1)
      INTEGER ND
      real*8 DEPTH(ND)
      real*8 AIRTMP
      INTEGER MS

      INTEGER BASINQ (MS)! Index of receiving basin

C Scale factor for normal fluxes:
      real*4 AIRTQF    (MS)! Degree of influence of air temperature
      real*4 QFW       (MS)! Water, including silicate and heat,
                         ! which are specified by concentration.
      real*4 QFP       (MS)! Phosphorus
      real*4 QFN       (MS)! Nitrogen
      real*4 QFODM     (MS)! Oxygen demand
      real*4 QFC       (MS)! Organic carbon
      real*4 QFS       (MS)! Silicate
      real*4 QFNH4     (MS)! Limit to fraction of N as ammonium.

      real*4 QCDETF   (MS) ! Upper limit to fraction of dissolved org. C 
                         ! entering detrital particulate fraction 
                         ! in the fjord
      real*4 QNDETF   (MS) ! Upper limit to fraction N as detritus
      real*4 QPDETF   (MS) ! Upper limit to fraction P as detritus
      real*4 QNCDET        ! N:C by weight with C to detritus
      real*4 QPCDET        ! P:C by weight with C to detritus

                         ! Transfer of discharges to other outlets: 
      real*4 QTRF      (MS,MS)! Transfer coefficients (fraction)
      INTEGER QTRNDX (MS,MS)! >0: Outlet index to transfer to
                         ! Only used to NS_TAB in second index.
                         ! Usage see code below
      real*4 QTR_Times (MS,2) ! Time of year to turn transfer on/off
      real*4 DEPTHQ    (MS,2) ! Specified outlet and intake depth
      real*4 QMIXIN    (MS)   ! Recipient water to mix into outlet
      real*4 MIXQTM    (MS,2) ! Time of year to turn intake on/off

      logical QTR_ACTIVE(MS)! on/off switch for transfers and mixing
                            ! of recipient waters into outlets.
                            ! transfers into outlet nr. i, including
                            ! mixing of recipient water, is only done
                            ! if QTR_ACTIVE(i) is .true.

      real*4 RNDFACTOR (MS)   ! Factor for random fluctuations

C Out:
      INTEGER RNFNDX(MS,2) ! Global layer index (1..NLI) for
                           ! (k,1): Receiving layer for final outlet
                           ! (k,2): Intake layer for recipient water
                           !        to mix into discharge water.

      real*8 QWATER   (MS,2) ! Sum of water flux over holes (m3/s)
                           ! (..,1): Land runoff (freshwater)
                           ! (..,2): Intake of recipient water
                           
      real*8 QTEMP    (MS)   ! Sum of heat influx over holes (deg.C*m3/d)

                   ! Flux of biological and chemical components as kg/d
      real*8 QPO4  (MS)  ! Phosphate
      real*8 QNO3  (MS)  ! Nitrate
      real*8 QNH4  (MS)  ! Ammonia
      real*8 QODM  (MS)  ! Oxygen demand
      real*8 QDOC  (MS)  ! Organic carbon to DOC
      real*8 QCDET (MS)  ! Organic carbon to Detritus
      real*8 QNDET (MS)  ! Nitrogen   to Detritus
      real*8 QPDET (MS)  ! Phosphorus to Detritus
      real*8 QSi   (MS)  ! Silicate
      real*8 QO2   (MS)  ! Free oxygen


C --------------------- local variables ------------------------
      INTEGER MONTH, IS
      real*8    XW,  xWday, XPO4,  XNO3,  XNH4,  X
      real*8    xODM, XDOC, XPOC, XCDET, XNDET, XPDET, XSi, XTMP
      INTEGER IX, I_PRIM, TYPE
      INTEGER IB,IL,ID
      INTEGER TransferNumber
      real*8    FractionLeft


C -------- description of pollution sources in table:

      integer ns_Oslofj
      parameter (ns_Oslofj = 15)

      integer ns_Tab
      parameter (ns_Tab =ns_Oslofj)


C ==========================================================
C         Discharges data on monthly basis of P,N and C:


      real*8 WATER_TAB(12,ns_tab) ! m3/s
      real*8 TOTP_TAB (12,NS_TAB) ! kg/mnd  or myg/l
      real*8 TOTN_TAB (12,NS_TAB) ! kg/mnd  or myg/l
      real*8 QFNH4_Tab   (NS_TAB) ! Max. ammonium fraction
      real*8 TOC_TAB  (12,NS_TAB) ! total organic carbon kg/d or mg/l
      real*8 POC_TAB  (12,NS_TAB) ! min. part. organic carbon kg/d or mg/l
                                ! (+ fraction CDETF of the rest) 
      real*8 Si_TAB   (12)        ! Silicate, equal in all sources
      real*8 ODM_Tab     (NS_TAB) ! Oxygen demand mg/l


            ! defined units, set in connection with data declarations 

      integer W_Unit(NS_TAB), 
     &        P_Unit(NS_TAB), N_Unit(NS_TAB), C_Unit(NS_TAB)

C ---------- Temperature of discharge:
C  Inherent temperature before air temp. influence:
      real*8 INHERENT_TMP(12,NS_TAB)! deg.C


C ==========================================================
C Random components:
      
      integer nRandom
      parameter (nRandom =ns_Tab) 



      real*8  T_LAST       (    nRandom)
      real*8 RNF_RESP_FREQ (    nRandom)
      INTEGER IJ_RUNOFF  ( 2, nRandom)
      real*8    UC_RUNOFF  (98, nRandom)

      real*8 RNF_RANDOM(nRandom)
      real*8 xRand(nRandom)

      
      SAVE IJ_RUNOFF, UC_RUNOFF, T_LAST, RNF_RANDOM


C ==========================================================



C ==========================================================


C --------------- Discharge data:

      integer I, K
      real*8 Time_of_year

           ! Unit conversion to kg/d or m3/s: 
           ! Note: Index for mass/time units <=0,
           !             for mass/volume units 1:3
           !             and for water flow units 4:5

      integer Myg_l, Mg_L, kg_mnd, ton_mnd, kg_d, m3_s, litres_s, m3_day   
      parameter ( Myg_l = -1, Mg_L = 0, kg_mnd=1, ton_mnd=2, kg_d=3,
     &            m3_s = 4, litres_s = 5, m3_day = 6 )

           ! Unit conversion for water flow:
      real*8 conv_m3_day 
      parameter (conv_m3_day = 1.0/24.0/3600.0)

           ! Unit conversion array:
      real*8 ConvFactor(-1:6) / 1.0e-6, 1.e-3, 0.0333, 33.3, 1.0,
     &                        1.0, 0.001, conv_m3_day  /



C ######################################################################
C ######################################################################

      integer i_Oslofj, n_Oslofj, ns_Oslofj_12
      parameter (i_Oslofj = 1)
      parameter (n_Oslofj = i_Oslofj-1+ns_Oslofj)
      parameter (ns_Oslofj_12 = 12*(n_Oslofj-i_Oslofj+1))

! Representerer utslipp pr. 2004
! satt opp for miljømål Bunnefjorden


! Sources:

! Surface outlets:

!  1,Årungen+Gloslibekken            Q1 = 0.41*MCA005
!  2,Gjersjøelva+Ljanselva           Q2 = 0.59*MCA005+0.10*MCA006
!  3,Alna+Akerselva+Frognerb.        Q3 = 0.90*MCA006
!  4,Lysakere.+Hoffselva+Mærddalsb.  Q4 =      MCA007
!  5,Sandvikselva                    Q5 = 0.92*MCA008
!  6,Neselva                         Q6 = 0.08*MCA008
!  7,Askere.+Nærsnese.+Årosva.+Sætre Q7 =      MCA009

! Dived outlets:

!  8,BRA hovedutslipp,  50 m
!  9,BRA overløp,       25 m
!  10,VEAS hovedutslipp 42 m
!  11,VEAS overløp      25 m
!  12,NFR               25 m
!  13,Buhrestua r.a.    22 m
!  14,Fagerstrand RA,   25 m
!  15,Til Steilene nord  1 m?

! Tables (Month,Source):
! DISCHARGE OF WATER: litres/s

      DATA (W_Unit(i),i=i_Oslofj, n_Oslofj)/ns_Oslofj*litres_s/
      DATA ((WATER_TAB(k,i),k=1,12),i=i_Oslofj, n_Oslofj)
     &  /950, 1274, 1134, 1034, 1046, 869, 
     &        545, 705, 789, 733, 924, 867
     & , 1831, 2587, 2329, 2131, 1979, 1594
     &       ,  957, 1410, 1783, 1767, 2168, 1925
     & , 2904, 3632, 3851, 6740, 7817, 5989
     &       , 4127, 5072, 8620, 5013, 5470, 5443
     & , 2997, 2486, 3898, 7345, 8904, 4137
     &       , 3211, 4196, 4864, 4528, 4989, 3685
     & , 3069, 3180, 5375, 6133, 6105, 4070
     &       , 2709, 3756, 4346, 4498, 4667, 3056
     & , 238,  247,  417,  476,  474,  316
     &      ,  210,  292,  338,  349,  362,  237
     & , 3396, 3836, 4504, 4356, 4636, 3841
     &       , 2483, 3192, 3975, 3416, 3810, 3244	  
     & , 12*1500
     & , 14, 15, 11, 46, 79, 25
     &      ,17, 25, 28, 53, 28, 9
     & , 3520, 3472, 4208, 2880, 3835, 3284
     &       , 2664, 3270, 4702, 2304, 3210, 3386	 
     & , 0, 5, 1, 0, 0, 7, 0, 165, 750, 0, 1, 36
     & , 12*150
     & , 12*24
     & , 12*14
     & , 12*1.0
     &  /

      DATA (P_Unit(i),i=i_Oslofj, n_Oslofj) /ns_Oslofj*kg_mnd/,
     & ((TOTP_TAB(k,i),k=1,12),i=i_Oslofj, n_Oslofj)
     &  /197, 268, 226, 198, 149,  79
     &      ,  45,  79,  57,  93, 160, 173
     & , 117, 276, 181, 199, 213, 160
     &      , 105, 218, 215, 174, 157, 118
     & , 327, 623, 498, 1027, 836, 569
     &      , 488, 1033, 1104, 632, 666, 426
     & , 115, 105, 200, 309, 388, 265
     &      , 144, 544, 354, 369, 248, 144
     & , 185, 129, 252, 267, 158, 201
     &      , 149, 165, 170, 104, 166,  96
     & ,  10,  12,  18,  12,  28,  21
     &      ,  22,  31,  31,  19,  23,  10
     & , 202, 484, 527, 481, 168, 243
     &      , 145, 158, 557, 223, 478, 143	 
     & , 12*1070
     & , 49, 52, 39, 164, 283, 88
     &      ,60, 91, 101, 190, 99, 9
     & , 3467, 3292, 4298, 2178, 3076, 2492
     &       , 1591, 2603, 3279, 2130, 2110, 1944
     & , 0, 7, 2, 0, 0, 10, 0, 502, 2176, 0, 2, 141
     & , 12*98
     & , 12*15
     & , 12*10
     & , 12*0.11
     &  /


      DATA (N_Unit(i),i=i_Oslofj, n_Oslofj) /ns_Oslofj*ton_mnd/,
     & ((TOTN_TAB(k,i),k=1,12),i=i_Oslofj, n_Oslofj)
     &  /7.64, 10.24, 9.12, 8.31, 8.40, 6.98
     &       , 4.38, 5.66, 6.34, 5.89, 7.42, 6.96
     & , 6.56, 9.56, 8.37, 7.97, 7.70, 5.88
     &       , 3.18, 4.90, 5.62, 5.74, 7.05, 6.63
     & , 6.86, 9.94, 9.70, 14.46, 14.83, 11.39
     &       , 8.32, 11.65, 13.75, 10.01, 11.13, 9.82
     & , 4.30, 4.00, 6.11, 10.71, 12.38, 5.52
     &      , 3.72, 7.21, 6.66, 7.42, 7.52, 5.49
     & , 7.88, 8.08, 16.58, 15.10, 9.49, 8.23
     &       , 5.55, 6.91, 7.88, 9.09, 10.28, 8.40
     & , 0.62, 0.75, 1.41, 1.40, 2.06, 0.78
     &       , 0.51, 0.74, 0.80, 1.09, 1.13, 0.62
     & , 10.92, 11.98, 14.97, 15.07, 12.41, 9.68
     &        , 5.81, 5.94, 13.29, 7.84, 18.79, 9.21
     & , 12*35
     & , 0.4, 0.4, 0.3, 1.3, 2.2, 0.7
     &      ,0.5, 0.7, 0.8, 1.5, 0.8, 0.2
     & , 133, 121, 135, 99, 124, 107
     &      , 65, 89, 124, 118, 153, 132
     & , 0, 0.1, 0, 0, 0, 0.1, 0, 4, 17, 0, 0, 1
     & , 12*3.8
     & , 12*0.6
     & , 12*0.4
     & , 12*0.0024
     &  /

      DATA (C_Unit(i),i=i_Oslofj, n_Oslofj) /ns_Oslofj*ton_mnd/
      DATA ((TOC_TAB(k,i),k=1,12),i=i_Oslofj, n_Oslofj)
     &  /18.97, 25.43, 22.64, 20.63, 20.87, 17.34
     &        , 10.87, 14.07, 15.75, 14.62, 18.44, 17.3
     & , 33.37, 47.64, 41.98, 40.50, 39.50, 30.30
     &        , 18.05, 32.21, 41.64, 36.86, 44.32, 35.38
     & , 34.88, 47.66, 46.43, 100.19, 106.18, 74.94
     &        , 51.97, 78.39, 119.71, 64.61, 77.65, 68.90
     & , 64.79, 50.73, 70.81, 151.65, 198.88, 95.49
     &        , 73.65, 103.53, 129.07, 108.16, 128, 89.65
     & , 38.19, 32.97, 59.91, 93.79, 90.20, 60.12
     &        , 44.93, 67.18, 72.09, 73.44, 62.91, 36.44
     & ,  2.29,  3.46,  6.17,  7.53,  7.87,  4.59
     &        ,  4.42,  4.84,  5.51,  5.43,  6.86, 3.38
     & , 42.04, 54.44, 61.87, 57.30, 64.28, 52.02
     &        , 40.07, 49.03, 84.48, 60.44, 68.88, 47.93
     & , 12*47
     & , 0.5, 0.5, 0.4, 1.6, 2.7, 0.9
     &      ,0.6, 0.9, 1.0, 1.8, 1.0, 0.3
     & , 181, 147, 195, 116, 168, 146
     &      , 96, 132, 189, 126, 169, 149
     & , 0, 0.1, 0, 0, 0, 0.1, 0, 4, 17, 0, 0, 1
     & , 12*4.2
     & , 12*3.6
     & , 12*0.9
     & , 12*0.0089
     &  /


             ! Distribution on nitrate and ammonium
             ! (Maximum fraction found as ammonium):
             !     3-7 % for river runoff,
             !     40 % from treatment plants after denitrification
             !     75 % from untreated sewage
	     !     73 % for VEAS

      DATA (QFNH4_TAB(i),i=i_Oslofj, n_Oslofj)
     &  /   0.04, 0.04, 0.07, 0.04, 0.03, 0.03, 0.05
     &          , 0.4, 0.75, 0.73, 0.75, 0.4, 0.4, 0.4, 0.4 /


             ! particulate carbon:
             ! REMAINS TO BE SPECIFIED 

      DATA ((POC_TAB(k,i),k=1,12),i=i_Oslofj, n_Oslofj)/
     &   ns_Oslofj_12*0.0/


             ! TEMPERATURE: oC
             ! month specific temperature before influence of air temp.:

      DATA ((INHERENT_TMP(k,i),k=1,12),i=i_Oslofj, n_Oslofj) /
     &  84*0.0           ! la elver følger lufttemperatur (sett AIRQTF=1)
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7.
     & , 4*7.,9.,12.,3*14.,12.,10.,7./ !RA iflg. dok. fra OVA
                                       !antas felles for alt avløpsvann  


C Specified total (max) chemical oxygen demand in runoff: 
      DATA (ODM_TAB(i),i=i_Oslofj, n_Oslofj)/ ns_Oslofj*0.0/   ! mg/l
                        ! Discharge from paper mill


C -------------- characteristic of random variations (1/days):
      DATA (RNF_RESP_FREQ(i),
     &   i=i_Oslofj, n_Oslofj) / 
     &   15*0.15/ 


C Silicate in mgSiO2/l, same in all sources,
C Based on data for rivers Glomma and Drammenselva in Norway.
      DATA Si_TAB/ 2*3.0, 3*2.5, 4*2.0, 3*3.0 /



C ====================================================================
C - Check dimensions first time:
      LOGICAL DimensionsChecked /.false./, OK /.true./
      SAVE    DimensionsChecked
      integer N

      N = 0
      if (.not. DimensionsChecked) then
         N = N + Ns_Oslofj
         OK = OK .and. (N.eq.N_Oslofj )
         if (.not. OK) STOP 'Error in dimensions in RNF_ALL.FOR'
         DimensionsChecked = .true.
      endif
            ! residual of code in old version RNF_ALL,
            ! where inputs from more than one fjord could be included.
            ! can be simplified here.



C     call hello ('INFLUX'      )
$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) WRITE (DEBUG_UNIT,*)'======== NFLUX'
$ENDIF

C ---------- Check dimension:
      IF (NS_TAB .GT. MS) THEN
         write(*,*) '(NS_TAB=', NS_TAB, ')  > (MS=', MS,
     &              ' in subroutine RUNOFF'
!         Pause 'Press Enter to stop model'
         Stop
      endif
                  ! NS_TAB <= MS required if redirection for QTRNDX
                  ! to work on all sources.
                  ! Part of runoff can be transferred to
                  ! outlets NS_TAB+1 to MS as well.

C ---------- Month number:
      Month = T/365.0*12.0    ! Whole months relative to T = 0
      Month = MOD(Month,12)+1 ! Index 1..12

           ! Intake of recipient water, alone or mixed into freshwater:
      Time_of_year = mod(T,365.0D0)

$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) WRITE (DEBUG_UNIT,*)' T', T, ' MONTH', MONTH,
     &           ' Time_of_year ',Time_of_year
$ENDIF


C ---------- Initiate values:

C       ....... Water amount and quality variables of outlet:
      DO IS = 1,MS


$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) WRITE (DEBUG_UNIT,*) ' IS:', IS,
     &          ' MIXQTM: ', MIXQTM(IS,1), MIXQTM(IS,2),
     &          ' QWATER(IS,2):', QWATER(IS,2)
$ENDIF

         ! Intake of recipient water, alone or mixed into freshwater:

          IF (   QTR_ACTIVE(IS)
     &           .and. Time_of_year.ge.MIXQTM(IS,1) 
     &           .and. Time_of_year.le.MIXQTM(IS,2) ) then
              QWATER (IS,2) = max(0.0,QMIXin(IS))
          ELSE
              QWATER (IS,2) = 0.0
          ENDIF

         ! Initiate freshwater description to be added to below:
          QWATER (IS,1) = 0.0
          QPO4   (IS) = 0.0
          QNO3   (IS) = 0.0
          QNH4   (IS) = 0.0
          QODM   (IS) = 0.0
          QDOC   (IS) = 0.0
          QCDET  (IS) = 0.0
          QNDET  (IS) = 0.0
          QPDET  (IS) = 0.0
          QSi    (IS) = 0.0
          QTEMP  (IS) = 0.0
      END DO



C --------------------------------------------------------------
C Process sources described in tables above,
C and distribute them on outlets, 
C using correction factors QF.. (1=normal value)


C ======= Random variations:
      


      CALL RANDOM_SIGNALS( T, T_LAST, NRandom, RNF_RESP_FREQ, 
     &   IJ_RUNOFF, UC_RUNOFF, ONE_TO_ONE, RNF_RANDOM )          
       !   in RANDFNCT.FOR, propagates RNF_RANDOM 
       !   as an array of NS_TAB stochastic factors, 
       !   created to have response freq. RNF_Resp_FREQ, 
       !   approximately normally distributed, 
       !   with mean = 0 and st.dev = 1,   
       !   as a first order Ornstein Uhlenbeck series.


!=======================================================================
! Special code for Iddefjord:      

!=======================================================================
! Other fjords:
      
      DO IS = 1,NRandom
          ! xRand is now created from RNF_RANDOM
          ! as stochastic factor log-normally distributed 
          ! and scaled to have standard deviation RNDFACTOR on log scale
          ! (before modification for Tista) and mean value 1:
        if (RNDFACTOR(IS).gt.0.0) then
           xRand(IS) = 
     &           EXP( (RNF_RANDOM(IS)-RNDFACTOR(IS)/2)*RNDFACTOR(IS))
        else
           xRand(IS) = 1.0
        endif
      ENDDO

      
      DO IS = 1,NS_TAB

      ! total amount of water

          xW   = xRAND(IS)*QFW(IS)*WATER_TAB (MONTH,IS)
     &         *ConvFactor(W_Unit(IS))  ! convert to m3/s

          xWDay = xW*24.*3600 ! converted to m3/day
          
      ! ...... Scale total pollutant fluxes:

               ! P, N and C converted to kg/d:

          xPO4 = QFP(IS)* TOTP_TAB (MONTH,IS) * ConvFactor(P_Unit(IS))
          if (P_Unit(IS).le.0) xPO4 = xWDay * xPO4   ! from conc.

          xNO3 = QFN(IS)* TOTN_TAB (MONTH,IS) * ConvFactor(N_Unit(IS))
          if (N_Unit(IS).le.0) xNO3 = xWDay * xNO3   ! from conc.

               !  Org. C total and particulate:
          xDOC = QFC(IS) *TOC_TAB (MONTH,IS) * ConvFactor(C_Unit(IS))
          xPOC = min( xDOC, 
     &           QFC(IS) *POC_TAB (MONTH,IS) * ConvFactor(C_Unit(IS)))
          if  (C_Unit(IS).le.0) then !  from conc.
             xDOC = xWDay*xDOC 
             xPOC = xWDay*xPOC
          endif  


               ! Oxygen demand, SPECIFIED AS mg/l
          xODM = xWDay*QFODM(IS) * ODM_Tab(IS) * 0.001
          !kg/d=  factor   *(m3/d)*(g/m3)      *(kg/g)
          !                 (QFW is included in xW.)

          !       Silicate given as mgSiO2/l = g/m3 in Si_TAB,
          !       converted to kg Si/d:
          !                (m3/d)*(g/m3)
          xSi   = xWDay * QFS(IS) * Si_TAB (MONTH)
     &            *28.09/(28.09+32) * 0.001
          !          gSi/gSiO2      * (kg/g)
          !               (QFW is included in xW)

C     Temperature never below fresh water freezing point (0.0):
          xTMP = MAX( 0.0D0,
     &            INHERENT_TMP(MONTH,IS)
     &            + AIRTQF(IS)*(AIRTMP-INHERENT_TMP(MONTH,IS))
     &              ) * xWDay  ! degC *m3/s * s/day
C              = heat flux in units degC*m3/day

$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) THEN
         WRITE (DEBUG_UNIT, '(1X,''IS:'',I5/1x,A)' ) IS,' total fluxes:'
         WRITE (DEBUG_UNIT, '(3(1X,A6,1X,G15.7))')
     &   'xW   ', xW   ,
     &   'xRand', xRand(IS),
     &   'xPO4 ', xPO4 ,
     &   'xNO3 ', xNO3 ,
     &   'xODM ', xODM ,
     &   'xDOC ', xDOC ,
     &   'xSi  ', xSi  ,
     &   'xTMP ', xTMP
      Endif
$ENDIF


      ! ....... Distribute C,N,P between dissolved and particulate:
      !         minimum as specified in POC_TAB + fraction of
      !         dissolved residual as specified in QCDETF
      !         (xDOC contains total organic carbon at this point)
          xCDET = xPOC+(xDOC-xPOC) * QCDETF(IS)    !  detritus
          xDOC = xDOC - xCDET                      !  to DOC

          xNDET = MIN ( xNO3*QNDETF(IS), xCDET*QNCDET )  !  detritus
          xNO3  = xNO3 - xNDET                           !  dissolved

          xPDET = MIN ( xPO4*QPDETF(IS), xCDET*QPCDET )  !  detritus
          xPO4  = xPO4 - xPDET                           !  dissolved

      ! ....... Fraction NH4 of dissolved N:
          x = max( 0.0D0, min(1.0D0, QFNH4_Tab(IS)*QFNH4(IS) ) )


          xNH4 = x*xNO3
          xNO3 = xNO3 - xNH4


$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) THEN
         WRITE (DEBUG_UNIT, '(1X,A)' )
     &    ' redistributed on chemical forms:'
         WRITE (DEBUG_UNIT, '(3(1X,A6,1X,G15.7))')
     &   'QCDETF', QCDETF(IS),
     &   'QNDETF', QNDETF(IS),
     &   'QPDETF', QPDETF(IS),
     &   'QNCDET', QNCDET,
     &   'QPCDET', QPCDET,
     &   'x    ', x    ,
     &   'xODM  ', xODM ,
     &   'xDOC ', xDOC ,
     &   'xCDET', xCDET,
     &   'xNO3 ', xNO3 ,
     &   'xNH4 ', xNH4 ,
     &   'xNDET', xNDET,
     &   'xPO4 ', xPO4 ,
     &   'xPDET', xPDET,
     &   'xSi  ', xSi  ,
     &   'xTMP ', xTMP
      Endif
$ENDIF

C  ............ update influx values
                   ! Can move parts of runoff to other outlets:
          I_prim = max( 1, MIN(MS, IS))
          FractionLeft=1.0
          TransferNumber=0

          do while (FractionLeft.gt.0.0)

            TransferNumber = TransferNumber+1

            if (TransferNumber.le.MS) then 
                   ! receiving outlet:
               Ix = min(MS, QTRNDX (TransferNumber,IS))
                   ! terminate transfers from IS by first index <=0:
               if (Ix.le.0) TransferNumber = MS+1
            endif

            if (TransferNumber.le.MS) then 

                  ! skip to next transfer if the transfer has no effect
                  ! or transfers into Ix are inactive:
               if (Ix .eq. I_prim .or. (.not. QTR_ACTIVE(Ix))) CYCLE

                  ! transfer fraction x to other outlet IX
                  ! within specified time window within year: 
               if (      Time_of_year.ge.QTR_Times(IS,1) 
     &            .and. Time_of_year.lt.QTR_Times(IS,2) ) then
                  x=QTRF(TransferNumber,IS)
                  if ( x.lt. 0.0 ) then
                       ! abs(value)=water flow limit (m3/s), convert to fraction:
                     x = abs(x)/xW
                  endif
                  x = max(0.0D0,min(FractionLeft,x))
               else
                  x = 0.0
               endif
                  
               FractionLeft   = FractionLeft - x

            else ! end of transfer specs., send rest to outlet IS:

               Ix = I_prim
               x = FractionLeft
               FractionLeft=0.0           

            end if

            if (x.gt.0.0) then 
               ! fraction x transferred to discharge Ix:
               QWATER (Ix,1) = max(0.0D0,x * xW)   + QWATER (Ix,1)
               QPO4   (Ix) = x * xPO4 + QPO4   (Ix)
               QNO3   (Ix) = x * xNO3 + QNO3   (Ix)
               QNH4   (Ix) = x * xNH4 + QNH4   (Ix)
               QODM   (Ix) = x * xODM + QODM   (Ix)
               QDOC   (Ix) = x * xDOC + QDOC   (Ix)
               QCDET  (Ix) = x * xCDET+ QCDET  (Ix)
               QNDET  (Ix) = x * xNDET+ QNDET  (Ix)
               QPDET  (Ix) = x * xPDET+ QPDET  (Ix)
               QSi    (Ix) = x * xSi  + QSi    (Ix)
               QTEMP  (Ix) = x * xTMP + QTEMP  (Ix)
            end if 

          enddo
          


C --------------------------------------------------------------
      ENDDO  ! Of one Source
C --------------------------------------------------------------




C --------------------------------------------------------------
C Final processing of outlets:
      DO IS = 1,MS
C --------------------------------------------------------------


C ---------- Check/update depth index of outlet depth and
C            any intake depth for recipient water mix-in:
C                 (works regardless of input value of RNFNDX )

          IB = MAX( 1, MIN( NBI, BASINQ(IS) ) )
          BASINQ(IS) = IB ! Is reset within allowed range,
                          ! Necessary because it is used to address
                          ! Arrays in other parts of the model.

          DO TYPE = 1, 2  ! Outlet depth, and eventual mix-in intake:
             IL = MAX( INDXI(IB)+1,
     &                 MIN( INDXI(IB+1), RNFNDX(IS, TYPE) ) )
             ID = IL-INDXI(IB)
             DO WHILE ( IL .LT. INDXI(IB+1)
     &                  .AND. DEPTHQ(IS, TYPE) .GT. DEPTH(ID+1) )
                IL = IL +1
                ID = ID +1
             ENDDO
             DO WHILE ( IL .GT. INDXI(IB)+1
     &                 .AND. DEPTHQ(IS, TYPE) .LE. DEPTH(ID) )
                IL = IL -1
                ID = ID -1
             ENDDO
           ! .... Global index of layer containing DEPTHQ(IS, TYPE)
           !       (or bottom layer)
             RNFNDX (IS, TYPE) = IL
          ENDDO




C  ---------  Add oxygen contents (assume approximate saturation):
          xW   = QWATER(IS,1)*24.*3600.      ! = m3/d
          if (xW .eq. 0.0 ) then
               xTMP = 0.0 ! Not important, just to avoid strange values
          else
               xTMP = QTEMP (IS) /xW
          endif  ! = (degC*m3/d)/(m3/d) = degC

          QO2 (IS) = xW * OXYGEN_SATURATION ( 0.0D0, xTMP)*1.4290e-3
C          kg/d    = (m3/d) *(l/m3)                     *(kg/l)
          

$IF DEFINED DEBUG_INFLUX
      IF (DEBUG) THEN
         WRITE (DEBUG_UNIT, '(1X,A,I3)' ) ' result, Q... nr.:',IS
         WRITE (DEBUG_UNIT, '(3(1X,A,1X,G15.7))')
     &    'Temperature',xTMP,
     &    'QWATER (IS,1)',QWATER (IS,1),
     &    '       (IS,2)',QWATER (IS,2),
     &    'QPO4   (IS)',QPO4   (IS),
     &    'QNO3   (IS)',QNO3   (IS),
     &    'QNH4   (IS)',QNH4   (IS),
     &    'QODM   (IS)',QODM   (IS),
     &    'QDOC   (IS)',QDOC   (IS),
     &    'QCDET  (IS)',QCDET  (IS),
     &    'QNDET  (IS)',QNDET  (IS),
     &    'QPDET  (IS)',QPDET  (IS),
     &    'QSi    (IS)',QPDET  (IS),
     &    'QTEMP  (IS)',QTEMP  (IS),
     &    'QO2    (IS)',QO2    (IS)
      Endif
$ENDIF

      END DO

C     call hello ('RETURN from RUNOFF')
      
      END Subroutine


$IF DEFINED Iddefjorden
      SUBROUTINE R_TRANSFORM_TISTA(I,R)
      integer I
      real*8 R

      integer I_Iddefj_Value
      COMMON/Index_Iddefj/I_Iddefj_Value
      
      real*8 X

      if (I.eq.I_Iddefj_Value) then 
         x = R
         R = 0.196*R*((1.0+0.43*ABS(R)**3.4)**0.8)
         write(777,'(''R :'',G12.5,''-->'',G12.5)') X, R 
      Endif
      END subroutine

$ENDIF


! ===========================================================
! Called from EUTRO.CSL
! Activates or deactivates transfer and recipient water mix-in
! For specific outlets dependent on bottom water oxygen conditions

      subroutine CHECK_Q_TRANSFER
     &    ( dimMS, QTR_OXCOND, BASINQ,
     &      NLI, OXYG, NBI, INDXI,
     &      T, QTR_CheckInterval, QTR_MinActiveInterval,
     &      QTR_ACTIVE, QTR_TimeOff )

           ! In:
      integer dimMS
      integer BASINQ (dimMS)       ! Receiving basin for outlets
      real*4  QTR_OXCOND (dimMS,2) ! Oxygen trigger and target values
      integer NLI
      real*8  OXYG(NLI)                ! Water layer oxygen levels
      integer NBI
      integer INDXI(NBI+1)


      real*8    T                     ! Current time
      real*4    QTR_CheckInterval     ! Time interval for checking
                                    ! oxygen conditions
      real*4    QTR_MinActiveInterval ! minimum active period (days)
           ! Out:
      logical QTR_ACTIVE  (dimMS)  ! On/Off switch¬
      real*8    QTR_TimeOff (dimMS)  ! Time for next possible switch-off


             ! Time trigger for checking oxygen conditions:
      real*8 TriggerTime /0.0D0/
      logical FirstCall/.true./
      SAVE TriggerTime, FirstCall


      integer outlet, receiving_basin, bottom_layer
      real*8 deep_oxygen

$IF DEFINED CHECK_Q_TRANSFER
      integer DEBUG_UNIT /9/, i
      save DEBUG_UNIT
$ENDIF


      if (FirstCall) then
         do outlet = 1, dimMS
            QTR_ACTIVE(outlet) = (QTR_OXCOND(outlet,1).ge. 10.0)
            QTR_timeOff(outlet)= 0.0
         enddo
         FirstCall=.false.

      endif

             ! Ensure reasonable trigger time
             ! independently of resets in T (re-runs etc.)
      if ( TriggerTime .gt. T + QTR_CheckInterval ) then
        TriggerTime = T - MOD(T,dble(QTR_CheckInterval))
      endif
      if ( T .lt. TriggerTime) return
      TriggerTime = T - MOD(T,dble(QTR_CheckInterval)) + QTR_CheckInterval

$IF DEFINED CHECK_Q_TRANSFER
            write(DEBUG_UNIT,'(3(1x,A,G16.5))') 
     &           ' CHECK_Q_TRANSFER, T=:', T,
     &           ' TriggerTime:',TriggerTime,
     &           ' QTR_CheckInterval:',QTR_CheckInterval
$ENDIF


      do outlet = 1, dimMS

         if (QTR_OXCOND(outlet,1).lt.10.0) then

            receiving_basin = BASINQ(outlet)
            bottom_layer = INDXI(receiving_basin+1)
            Deep_Oxygen = OXYG(bottom_layer)

$IF DEFINED CHECK_Q_TRANSFER
            write(DEBUG_UNIT,'(3(1x,A,I6),A,G15.6)') ' outlet', outlet,
     &        ' receiving_basin:',receiving_basin,
     &        ' bottom_layer:',bottom_layer,
     &        ' Deep_Oxygen:', Deep_Oxygen
     
$ENDIF

            if ( .not. QTR_ACTIVE (outlet)) then

                      ! turn on dependent of oxygen conditions
               if (Deep_Oxygen .lt. QTR_OXCOND (outlet,1) ) then
                  QTR_ACTIVE (outlet) = .true.
                  QTR_TimeOff(outlet) = T + QTR_MinActiveInterval
                    ! and keep going for a specified mimimum interval
               endif

            else       ! turn off after oxygen target is achieved,
                       ! but only after a specified minimum time
               if ( Deep_Oxygen .ge. QTR_OXCOND (outlet,2)
     &              .and. T .ge. QTR_TimeOff(outlet) ) then
                  QTR_ACTIVE (outlet) = .false.
               endif
            endif

$IF DEFINED CHECK_Q_TRANSFER
            write(DEBUG_UNIT,'(2(1x,A,G15.6))')
     &          ' QTR_ACTIVE:'   , QTR_ACTIVE (outlet),
     &          '   QTR_TimeOff:', QTR_TimeOff(outlet)
$ENDIF

         endif

      enddo

      end subroutine


      end Module
