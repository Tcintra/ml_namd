      SUBROUTINE adiab_states(ad,qp,qe,vlow,sum)

      implicit none
!     include 'param.h'
      integer i,j,p,r,np,ne,k,ad
      integer f
      real*8 pi,dpi,qp,qe,qp1,qp2,qe1,qe2,fact
      PARAMETER (np=24,ne=24)
      PARAMETER (pi=3.14159d0)
      real*8 phip,phie
      real*8 b1,b2,xp1,xp2,xe1,xe2,sum
      real*8 vlow((np+1)*(ne+1),(np+1)*(ne+1)),pl
      
!parameters

      b1=8d0
      b2=1d0
      qp1=3d0
      qp2=3d0
      qe1=3d0
      qe2=3d0
      dpi=dsqrt(pi)
            
!  Calculating adiabatic states based on linear combination of basis functions
  
      sum=0.d0
      k=0
      DO i=0,np
       DO j=0,ne
       k=k+1
       xp1=b1*(qp-qp1)
       xe1=b2*(qe-qe1)
       phip=(2**i*fact(i)*dpi)**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**j*fact(j)*dpi)**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
       sum=sum+vlow(k,ad)*phip*phie
       ENDDO
      ENDDO  

      DO i=0,np
       DO r=0,ne
       k=k+1
       xp1=b1*(qp-qp1)
       xe2=b2*(qe-qe2)
       phip=(2**i*fact(i)*dpi)**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**r*fact(r)*dpi)**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
       sum=sum+vlow(k,ad)*phip*phie
       ENDDO
      ENDDO 

      DO p=0,np
       DO j=0,ne
        k=k+1
        xp2=b1*(qp-qp2)
        xe1=b2*(qe-qe1)
        phip=(2**p*fact(p)*dpi)**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**j*fact(j)*dpi)**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
        sum=sum+vlow(k,ad)*phip*phie
       ENDDO
      ENDDO

      DO p=0,np
       DO r=0,ne
        k=k+1
        xp2=b1*(qp-qp2)
        xe2=b2*(qe-qe2)
        phip=(2**p*fact(p)*dpi)**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**r*fact(r)*dpi)**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
        sum=sum+vlow(k,ad)*phip*phie
       ENDDO
      ENDDO     
      END

         
