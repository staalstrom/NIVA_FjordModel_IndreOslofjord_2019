      Module fx_WindCurr
      use ModelParam_Physics
      use ModelVar_Topography
      use ModelVar_HydroBioChem
      implicit none   
      
      
      contains


C ===================================================================
C Eutrophication model   - File:    WINDCURR.FOR
C                                   Birger Bjerkeng, NIVA.

C Subroutine calculating and exporting wind-induced current at given depth:


$define DEBUG_WIND


      SUBROUTINE GET_WIND_CURRENT( debug_unit, DEBUG, 
     &            I_CONN, XMIX_A, XMIX_B, 
     &            WIND_CURRENT, I_MAX_WINDC )

      
      integer debug_unit
      LOGICAL DEBUG
      INTEGER I_CONN
      real*8 XMIX_A, XMIX_B, WIND_CURRENT
      real*8 I_MAX_WINDC
      
      
      real*8 pi2_ON_360
      PARAMETER (PI2_ON_360 = 2.0*3.141592654/360.0)
      
      INTEGER I_L
      real*8 A_V, X_D


! WIND INDUCED SURFACE CURRENT: >0 FROM A TO B, <0 FROM B TO A:

      A_V = WVDIR(I_CONN)*PI2_ON_360 ! Direction of connection in radians
      WIND_CURRENT = WVFAC(I_CONN) * (WINDN*COS(A_V)+WINDE*SIN(A_V))


! Find layer containing lower depth limit WVHMIN() of wind-driven current,
! and calculate how large part of the layer is below the depth limit:

      DO I_L = 2,ND-1
         IF (DEPTH(I_L).GE.WVHMIN(I_CONN)) EXIT
      ENDDO   ! ends with I_L=ND if WVHMIN(I_CONN) > DEPTH(ND-1)
      
      if (I_L.eq.2) then
         X_D=0
      else
         X_D = MAX(0.0d0,(DEPTH(I_L)-WVHMIN(I_CONN)))/(DEPTH(I_L)-DEPTH(I_L-1))
      endif
 
      I_MAX_WINDC = FLOAT(I_L) - X_D
      
$if defined DEBUG_WIND

!      IF (DEBUG) then

         WRITE ( DEBUG_UNIT , * )
     &   ' ******* subroutine GET_WIND_CURRENT', 
     &   'I_CONN=',I_CONN, 'WVDIR(I_CONN)=',WVDIR(I_CONN),'WINDN=',WINDN,'WINDE=',WINDE,
     &   'WIND_CURRENT=',WIND_CURRENT,'I_L=', I_L, 'I_MAX_WINDC=', I_MAX_WINDC
!      endif
$endif

      END Subroutine

      End Module
