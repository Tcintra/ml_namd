      SUBROUTINE FORCE(R,VS,F)

      implicit none
      include 'param.h'
      REAL*8 VS,F(nat)
      REAL*8 MASS(nat),WS,R(nat),R0
      include 'solvent.in'

!    Calculate solvent force

      VS = 0.5*mass(1)*(WS**2)*(R(1)-R0)**2      
      F(1) = -mass(1)*(WS**2)*(R(1)-R0)
      RETURN
      END SUBROUTINE FORCE

