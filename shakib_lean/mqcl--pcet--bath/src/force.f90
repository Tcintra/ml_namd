      SUBROUTINE FORCE(qs,c,Wb,Qb,Fb,Fs)

      implicit none
      include 'param.h'
      integer I
      REAL*8 Fb(mode),Fs(nat),sum,ms,ws,s1,s2
      REAL*8 Mb,Wb(mode),Qb(mode)
      REAL*8 qs(nat),c(mode),wc,coup
      include 'sol_bath.in'

!    Calculate solvent force

      sum=0.d0
      do I=1,mode
         Fb(I) = -Mb*Wb(I)**2*(Qb(I)-(c(I)*qs(1)/(Mb*Wb(I)**2)))
         sum = sum + c(I)*(Qb(I)-(c(I)*qs(1)/(Mb*Wb(I)**2)))
      enddo
      Fs(nat)=sum
   
      RETURN
      END SUBROUTINE FORCE

