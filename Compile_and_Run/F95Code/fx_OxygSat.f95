      Module fx_OxygSat
      
     
      implicit none

      contains


C--------------------------------------------------------------------
      real*8 FUNCTION OXYGEN_SATURATION ( S, T )

C Saturation oxygen contents C0 in ml/l at pressure = 760 MM HG:

C This subroutine is also called by HYDGRAPH, the
C Freshwater input subroutine.

      real*8 S, T  ! SALINITY (o/oo) and TEMPERATURE (oC)

C       LN(C)=A1+A2(100/T)+A3*LN(T/100)+A4(T/100)
C             +S(B1+B2(T/100)+B3(T/100)**2)
C     FORMULA FROM 'INTERNATIONAL OCEANOGRAPHIC TABLES' VOL 2, N.I.O.
C     AND UNESCO 1973.
C TRANSFORMATION FROM ML/L AT STANDARD T AND P TO MG/L: FACTOR FMG

C     real*8 FMG /1.429/ ! ml/l --> mg/l   (not used)

      real*8 TA


C CONSTANTS:

      real*8 A1,A2,A3,A4
      real*8 B1,B2,B3

      DATA A1,A2,A3,A4 / -173.4292 , 249.6339 , 143.3483 , -21.8492 /
      DATA B1,B2,B3 / -0.033096 , 0.014259 , -0.0017 /


C ABSOLUTE TEMPERATURE T/100:
          TA = (T+273.15)/100.

C Solubility in ml/l at 760 mm hg pressure:
          OXYGEN_SATURATION =
     &          EXP( A1 + (  A2/TA)+A3*ALOG(TA)+A4*TA
     &                      + S*(B1+TA*(B2+B3*TA)) )

      end function

      end Module