      PROGRAM integrateS

! This program integrates the quantum Hamiltonian of the 1st submatrix.

      implicit none
      integer i,j,p,r,np,ne,nmax
      integer mp,me,le,lp,f,kf
      real*8 qp,qe,qp1,qp2,qe1,qe2,fact
      real*8 pi
      PARAMETER (nmax=30) 
      PARAMETER (pi=3.14159d0)
      DOUBLE PRECISION :: phip1(0:nmax),phip2(0:nmax),phie1(0:nmax),phie2(0:nmax)
      DOUBLE PRECISION :: pl(0:nmax),S_ac(0:nmax,0:nmax,0:nmax,0:nmax), &
      S_bc(0:nmax,0:nmax,0:nmax,0:nmax),S_ad(0:nmax,0:nmax,0:nmax,0:nmax), &
      S_bd(0:nmax,0:nmax,0:nmax,0:nmax),S_inta(0:nmax,0:nmax,0:nmax,0:nmax), &
      S_intb(0:nmax,0:nmax,0:nmax,0:nmax),S(0:nmax,0:nmax,0:nmax,0:nmax),dpl(0:nmax)
      real*8 b1,b2,xp1,xp2,xe1,xe2,sum
      real*8 a,b,c,d,hp,he,amp1,amp2
                
! Read number of basis functions
      OPEN(unit=1,file="basis.in") 
      READ(1,*) 
      READ(1,*) np,ne,qp1,qe1,qp2,qe2
         
! Parameters
     
      amp1=1.d0
      amp2=1.d0
      b1=7.5d0
      b2=0.65d0
      kf=4
                 
! Integration parameters 
!      a=1.8d0
!      b=4.2d0
      a=2d0
      b=4d0
      c=-6.d0
      d=12.d0
      mp=100
      me=700
      hp=(b-a)/dble(mp)
      he=(d-c)/dble(me)
      
!--------------------------------------------------------------------
! Calculate matrix elements regarding overlap of the harmonic 
! oscillators
!--------------------------------------------------------------------
      
      do i=0,np
       do j=0,ne
        do p=0,np
         do r=0,ne

          qe=c
          xe1=b2*(qe-qe1)
          xe2=b2*(qe-qe2) 

          qp=a
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_ac(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
          
                   
          qp=b
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_bc(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
          

          sum=0.5d0*(S_ac(i,j,p,r)+S_bc(i,j,p,r))
          do lp=1,mp-1
            qp=a+lp*hp
            xp1=b1*(qp-qp1)
            xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          sum=sum+phip1(i)*phie1(j)*phip2(p)*phie2(r)
          enddo
 

          qe=d
          xe1=b2*(qe-qe1)
          xe2=b2*(qe-qe2) 

          qp=a
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_ad(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
         
          qp=b
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_bd(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
          
          sum=sum+0.5d0*(S_ad(i,j,p,r)+S_bd(i,j,p,r))
          do lp=1,mp-1
            qp=a+lp*hp
            xp1=b1*(qp-qp1)
            xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          sum=sum+phip1(i)*phie1(j)*phip2(p)*phie2(r)
          enddo
          

          do le=1,me-1
          qe=c+le*he
          xe1=b2*(qe-qe1)
          xe2=b2*(qe-qe2)

          qp=a
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_inta(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
          
          qp=b
          xp1=b1*(qp-qp1)
          xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          S_intb(i,j,p,r)=phip1(i)*phie1(j)*phip2(p)*phie2(r)
         
          sum=sum+0.5d0*(S_inta(i,j,p,r)+S_intb(i,j,p,r))
          do lp=1,mp-1
            qp=a+lp*hp
            if (qp==qe) cycle
            xp1=b1*(qp-qp1)
            xp2=b1*(qp-qp2)
          CALL mothpl(kf,i,xp1,pl,dpl)
          phip1(i)=amp1*(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i)*dexp(-(b1**2*(qp-qp1)**2)/2)
          CALL mothpl(kf,j,xe1,pl,dpl)
          phie1(j)=amp2*(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j)*dexp(-(b2**2*(qe-qe1)**2)/2)
          CALL mothpl(kf,p,xp2,pl,dpl)
          phip2(p)=amp1*(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p)*dexp(-(b1**2*(qp-qp2)**2)/2)
          CALL mothpl(kf,r,xe2,pl,dpl)
          phie2(r)=amp2*(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r)*dexp(-(b2**2*(qe-qe2)**2)/2)
          sum=sum+phip1(i)*phie1(j)*phip2(p)*phie2(r)
          enddo
          enddo
          S(i,j,p,r)=hp*he*sum
          write(100,*) i,j,p,r,S(i,j,p,r)
         enddo
        enddo
       enddo
      enddo

      END

      FUNCTION fact(n) 
      IMPLICIT NONE 
      INTEGER n,f
      real*8 fact
      fact = 1
      DO f = 1, n
      fact = fact * f
      ENDDO
      RETURN
      END FUNCTION fact

      SUBROUTINE mothpl(kf,n,x,pl,dpl)

      !       ==========================================================
      !       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
      !                or Ln(x) or Hn(x), and their derivatives
      !       Input :  KF --- Function code
      !                       KF=1 for Chebyshev polynomial Tn(x)
      !                       KF=2 for Chebyshev polynomial Un(x)
      !                       KF=3 for Laguerre polynomial Ln(x)
      !                       KF=4 for Hermite polynomial Hn(x)
      !                n ---  Order of orthogonal polynomials
      !                x ---  Argument of orthogonal polynomials
      !       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
      !                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
      !       =========================================================
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER, INTENT(IN)                      :: kf
      INTEGER, INTENT(IN)                      :: n
      INTEGER nmax
      DOUBLE PRECISION, INTENT(IN)             :: x
      PARAMETER (nmax=30) 
      DOUBLE PRECISION :: pl(0:nmax)
      DOUBLE PRECISION :: dpl(0:nmax)

      a=2.0D0
      b=0.0D0
      c=1.0D0
      y0=1.0D0
      y1=2.0D0*x
      dy0=0.0D0
      dy1=2.0D0
      if (n==0) then
      pl(0)=1.0D0
      dpl(0)=0.0D0
      else if (n==1) then
      pl(1)=2.0D0*x
      dpl(1)=2.0D0
      else if (n>=2) then
      do k=2,n
        c=2.0D0*(k-1.0D0)
        yn=(a*x+b)*y1-c*y0
        dyn=a*y1+(a*x+b)*dy1-c*dy0
        pl(n)=yn
        dpl(n)=dyn
        y0=y1
        y1=yn
        dy0=dy1
        dy1=dyn
      ENDDO
      ENDIF
      RETURN
      END SUBROUTINE mothpl


