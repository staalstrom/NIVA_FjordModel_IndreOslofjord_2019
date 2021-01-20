      Module fx_SigmaT

      use fx_RunControl, only: DEBUG_UNIT
      use IEEE_ARITHMETIC
      
      implicit none

$undefine sigma_test

      logical :: Found_Nan = .false.

      Contains

      
C   Basic hydrographic density relations
C ===========================================================
C SIMPLIFIED DENSITY FORMULA (sigma-t units)  (TO BE ADJUSTED):
      SUBROUTINE SIGMAT (SALIN, TEMP, N, DENS )      

      
      integer N
      real*8 SALIN(N), TEMP(N)
      real*8, intent(Out) :: DENS(N)


C --------------- local variables ------------------
      INTEGER K
      real*8 S,T

      DO K = 1,N
         S = SALIN(K)

         if (IEEE_IS_NAN(S) .and. .not.Found_NaN) then
            write(*,*) 'S is NaN in Sigmat'
            Found_Nan = .true.
         endif

         T = TEMP(K)
         DENS(K) = - 0.117918 + S*0.805700
     &             + T*(  0.062968 - 0.003031*S
     &                   + T*( -0.008592 + 0.000035*S + 0.000060*T ) )
      END DO

$if defined  sigma_test
      write (DEBUG_UNIT, * ) 'SUBRUTINE SIGMAT:'
      write ( DEBUG_UNIT,'(1X,A5,3A16)')
     &          'K:','SALINITY:','TEMPERATURE:','DENSITY:'
      write ( DEBUG_UNIT,'(1x,I5,3F16.8)' )
     &         ( K,SALIN(K), TEMP(K), DENS(K),K=1,N)
$endif

      END SUBROUTINE




C ===========================================================
C Stability factors for surface exchange of water and heat,
C only used for top layers:

      SUBROUTINE DDENS ( DENSDS, DENSDT, SALIN, TEMP )
      
      
      real*8 DENSDS, DENSDT, SALIN, TEMP

      real*8 TF

      TF     =   ( TEMP-4.0+0.21*SALIN )
      DENSDS = (0.76 - 0.007*2*TF*0.21)
      DENSDT = -0.007*2*TF

      END SUBROUTINE

      END Module
