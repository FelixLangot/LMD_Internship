Tabs en K
p(iz) en hPa
qsat en g/kg

if (Tabs(it,iz,ix,iy)- 273.15.gt.0.0) then  
            qsat=rslf(p(iz)*100.0,Tabs(it,iz,ix,iy))*1.e3  
       else  
            qsat=rsif(p(iz)*100.0,Tabs(it,iz,ix,iy))*1.e3  
       endif  

!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
! copied from SAM
!
     REAL FUNCTION RSLF(P,T)
     IMPLICIT NONE
     REAL, INTENT(IN):: P, T
       ! P in Pa and T in K; RSLF in kg/kg
     REAL:: ESL,X
     REAL, PARAMETER:: C0= .611583699E03
     REAL, PARAMETER:: C1= .444606896E02
     REAL, PARAMETER:: C2= .143177157E01
     REAL, PARAMETER:: C3= .264224321E-1
     REAL, PARAMETER:: C4= .299291081E-3
     REAL, PARAMETER:: C5= .203154182E-5
     REAL, PARAMETER:: C6= .702620698E-8
     REAL, PARAMETER:: C7= .379534310E-11
     REAL, PARAMETER:: C8=-.321582393E-13
     X=MAX(-80.,T-273.16)
!      ESL=612.2*EXP(17.67*X/(T-29.65))
     ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
     RSLF=.622*ESL/(P-ESL)
!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!        / T - 9.44523 * ALOG(T) + 0.014025 * T))


     END FUNCTION RSLF
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
! copied from SAM
!
     REAL FUNCTION RSIF(P,T)
     IMPLICIT NONE
     REAL, INTENT(IN):: P, T
       ! P in Pa and T in K; RSIF in kg/kg
     REAL:: ESI,X
     REAL, PARAMETER:: C0= .609868993E03
     REAL, PARAMETER:: C1= .499320233E02
     REAL, PARAMETER:: C2= .184672631E01
     REAL, PARAMETER:: C3= .402737184E-1
     REAL, PARAMETER:: C4= .565392987E-3
     REAL, PARAMETER:: C5= .521693933E-5
     REAL, PARAMETER:: C6= .307839583E-7
     REAL, PARAMETER:: C7= .105785160E-9
     REAL, PARAMETER:: C8= .161444444E-12
     X=MAX(-80.,T-273.16)
     ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
     RSIF=.622*ESI/(P-ESI)
!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)
     END FUNCTION RSIF