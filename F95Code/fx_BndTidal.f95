      Module fx_BndTidal

      implicit none    


           ! Ratio between change in water level Z and air pressure p:
           ! for Z in m, p in hPa = 100 Pa = mBar

      real*8 Pressure_response / 0.01 / !( dz:dp) in m/mbar)



           ! Water level data with tidal component coefficients
           ! at main station:

      real*8 Z0                     ! average height cm
      INTEGER COMPONENTS
      PARAMETER ( COMPONENTS = 7 )
      real*8 PERIOD   (COMPONENTS)  ! hours
      real*8 AMPLITUDE(COMPONENTS)  ! cm
      real*8 PHASE    (COMPONENTS)  ! 360-degrees

           ! Correction for secondary station:
      real*8 HeightFactor           ! for tidal variation
      real*8 Z_Correction           ! cm, added after Heightfactor is applied

           ! for Oslofjord - from mb to m (Ola M. Johannesen)
           !!!!!!! should have site-specific constant


! Tidal level specifications for Inner Oslofjord
      integer i
      DATA Z0/ 65.0 /            ! cm
      DATA (PERIOD(i), AMPLITUDE(i), PHASE(i), i=1,COMPONENTS)/
     &       12.42       , 13.4, 152.,  !M2
     &       12.00       , 3.4 , 101.,  !S2
     &       12.66       , 3.3 , 101.,  !N2 
     &       11.97       , 0.9 ,  76.,  !K2
     &       23.93       , 0.5 , 172.,  !K1 
     &       25.82       , 2.1 , 299.,  !O1
     &       4380.       , 11.2, 282.   !Sa (Semiannual?)
     &     /

           ! Corrections zeroed:
      DATA HeightFactor / 1 /
      DATA Z_Correction / 0.0 / !cm

      real*8, private    ::  FMEAN   ! amplitude weighted mean angular frequency
                                   ! to be calculated in first call  
      logical, private :: Initiated = .false.   ! set .true. when done


      contains


C -----------------------------------------------------------

         ! New version in Egersund

      SUBROUTINE TIDAL( T, Air_Pressure, Z, DZ_DT, EMIX_RELATIVE )
      
      
      real*8 T              ! time in days
      real*8 Air_Pressure   ! pressure in mBar (1 mBar = 100 Pa = 100 N/m2)

      real*8 Z,DZ_DT
      real*8 EMIX_RELATIVE
C  ---- Tidal surface level (m) and time derivative (m/day):

      INTEGER I
      real*8 pi_2 / 6.283 185 307 180 /
      real*8 Last_time /0/, Last_pressure /-9999./
      Save Last_time, Last_pressure

      real*8 X, F, A_SUM, A_SIN_SUM, A_COS_SUM

      if (.not.Initiated) then

              ! Calculate FMEAN
              ! = amplitude weighted mean angular frequency in 1/day
              !   using only semidiurnal components

         FMEAN = 0.0
	      A_SUM = 0.0
	      DO I = 1, 4
	         FMEAN = FMEAN + AMPLITUDE(I)* (24./PERIOD(I))
	         A_SUM = A_SUM + AMPLITUDE(I)
	      ENDDO
	      FMEAN = PI_2*FMEAN/A_SUM
	      Initiated = .true.

      ENDIF


              ! tidal surface height and time derivative
              ! as sum over all components:


      A_SIN_SUM = 0.0
      A_COS_SUM = 0.0
      A_SUM = 0.0

      Z = Z0
      DZ_DT = 0.0


!            write(333,*) '*****TIDAL for T=', T
!            write(333,'(1x,A4,4(A12))') 'I', 'F', 'Z', 'DZ_DT', 'SIN'

      DO I = 1, COMPONENTS

                ! Sum tidal components as level in cm 
                ! and time derivative as cm/day

         F = Pi_2*24.0/PERIOD(I)     ! Frequency as 1/days
         Z     = Z     + AMPLITUDE(I)*SIN(F*T+PHASE(I))
         DZ_DT = DZ_DT + AMPLITUDE(I)*COS(F*T+PHASE(I))*F


!            write(333,'(1x,I4,4(G12.5))')
!     &                   I, F, Z, DZ_DT, SIN(F*T+PHASE(I))


                ! For Semi-diurnal components, prepare for 
                ! calculation of effective amplitude at current time
                ! to return parameter for tidal energy (see below)

         if ( I .le. 4 ) then
            A_SIN_SUM = A_SIN_SUM + AMPLITUDE(I)*SIN( (F-FMEAN)*T )
            A_COS_SUM = A_COS_SUM + AMPLITUDE(I)*COS( (F-FMEAN)*T )
            A_SUM = A_SUM + AMPLITUDE(I)**2
         endif

      ENDDO

               ! Correct level and rescale from cm to meters:

      Z = (Z_Correction+Z*HeightFactor)/100.0
      DZ_DT = DZ_DT*HeightFactor/100.0



               ! Effect of Air Pressure:
               ! (ignores retardation in pressure response)

                        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                        ! error in old code, corrected October 2004;
                        ! correction extracted from tidal component loop:
                        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      Z = Z - Pressure_response * (Air_Pressure-1013.0)
                                  ! deviation from normal atm. pressure

      if (T .gt. Last_time .and. T .lt. Last_time + 0.3
     &     .and. Last_Pressure .gt.0.0   ) then

               ! Avoid disaster on restart by taking into 
               ! account change rates in air pressure
               ! in the dz_dt calculation

         X = (Air_pressure-Last_pressure)/(T-Last_time)
         DZ_DT = DZ_DT - Pressure_response *X
      Endif

      Last_time = T
      Last_pressure = Air_Pressure



               ! effective mixing energy
               ! by internal waves created by tidal oscillation:

      EMIX_RELATIVE = (A_SIN_SUM**2 + A_COS_SUM**2 ) / A_SUM

               !   [ current effective amplitude**2] / [mean ampl**2]
               ! should represent tidal energy input (for mixing)
               ! relative to mean value.

               ! =========================================================
               ! Effective amplitude is found by writing each term in Z as
               !     A(i)*sin[F_mean*T+D(I)])
               ! where D(I) = (F-F_mean)*T 
               ! is a term varying slowly with time.

               ! Z(t) approx =  [SUM(A(i)*cos(D(i))] * sin(F_mean*T)
               !               +[SUM(A(i)*sin(D(i))] * cos(F_mean*T)
               ! together with
               !    a*sin(x) + b*cos(x) = A_eff*sin[ x + arccos(a/A_eff) ]
               !    where A_eff = sqrt(a**2+b**2)
               ! =========================================================

      END Subroutine


      end module