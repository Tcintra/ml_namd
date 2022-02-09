      SUBROUTINE coupling(mm,nn,qp,qe,vlow,adiabm,adiabn)

      implicit none
!      include 'param.h'
      integer i,j,p,r,np,ne,km,kn
      integer n,f,nn,mm
      real*8 pi,qp,qe,qp1,qp2,qe1,qe2,fact
      PARAMETER (np=24,ne=24)
      PARAMETER (pi=3.14159d0)
      real*8 phip,phie
      real*8 b1,b2,xp1,xp2,xe1,xe2,summ,sumn
      real*8 adiab1,adiab2,adiab3,adiab4,adiabm,adiabn
      real*8 vlow((np+1)*(ne+1),(np+1)*(ne+1)),pl
      
!parameters

      b1=8d0
      b2=1d0
      qp1=3d0
      qp2=3d0
      qe1=3d0
      qe2=3d0
            
!  Calculating adiabatic states based on linear combination of basis functions
  
      summ=0
      km=0
      DO i=0,np
       DO j=0,ne
       km=km+1
       xp1=b1*(qp-qp1)
       xe1=b2*(qe-qe1)
       phip=(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
       adiab1=vlow(km,mm)*phip*phie
       summ=summ+adiab1
       ENDDO
      ENDDO  

      DO i=0,np
       DO r=0,ne
       km=km+1
       xp1=b1*(qp-qp1)
       xe2=b2*(qe-qe2)
       phip=(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
       adiab2=vlow(km,mm)*phip*phie
       summ=summ+adiab2
       ENDDO
      ENDDO 

      DO p=0,np
       DO j=0,ne
        km=km+1
        xp2=b1*(qp-qp2)
        xe1=b2*(qe-qe1)
        phip=(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
        adiab3=vlow(km,mm)*phip*phie
        summ=summ+adiab3
       ENDDO
      ENDDO

      DO p=0,np
       DO r=0,ne
        km=km+1
        xp2=b1*(qp-qp2)
        xe2=b2*(qe-qe2)
        phip=(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
        adiab4=vlow(km,mm)*phip*phie
        summ=summ+adiab4
       ENDDO
      ENDDO
      adiabm=summ


      sumn=0
      kn=0
      DO i=0,np
       DO j=0,ne
       kn=kn+1
       xp1=b1*(qp-qp1)
       xe1=b2*(qe-qe1)
       phip=(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
       adiab1=vlow(kn,nn)*phip*phie
       sumn=sumn+adiab1
       ENDDO
      ENDDO  

      DO i=0,np
       DO r=0,ne
       kn=kn+1
       xp1=b1*(qp-qp1)
       xe2=b2*(qe-qe2)
       phip=(2**i*fact(i)*dsqrt(pi))**(-0.5)*b1**0.5*pl(i,xp1)*exp(-(b1**2*(qp-qp1)**2)/2)
       phie=(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
       adiab2=vlow(kn,nn)*phip*phie
       sumn=sumn+adiab2
       ENDDO
      ENDDO 

      DO p=0,np
       DO j=0,ne
        kn=kn+1
        xp2=b1*(qp-qp2)
        xe1=b2*(qe-qe1)
        phip=(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**j*fact(j)*dsqrt(pi))**(-0.5)*b2**0.5*pl(j,xe1)*exp(-(b2**2*(qe-qe1)**2)/2)
        adiab3=vlow(kn,nn)*phip*phie
        sumn=sumn+adiab3
       ENDDO
      ENDDO

      DO p=0,np
       DO r=0,ne
        kn=kn+1
        xp2=b1*(qp-qp2)
        xe2=b2*(qe-qe2)
        phip=(2**p*fact(p)*dsqrt(pi))**(-0.5)*b1**0.5*pl(p,xp2)*exp(-(b1**2*(qp-qp2)**2)/2)
        phie=(2**r*fact(r)*dsqrt(pi))**(-0.5)*b2**0.5*pl(r,xe2)*exp(-(b2**2*(qe-qe2)**2)/2)
        adiab4=vlow(kn,nn)*phip*phie
        sumn=sumn+adiab4
       ENDDO
      ENDDO
      adiabn=sumn
      END

        
