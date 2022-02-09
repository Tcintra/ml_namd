      SUBROUTINE ADIAB_NVE(Kin,Pot,Vpes1,Vpes2,S,adiab,adiabm,adiabn, &
                 Hmat,Smat,vect,vlow,elow1,elow2,elow3,HF,d_tot,R)

!====================================================================
! This subroutine will calculate the adiabatic states,the
! Hellmann-Feynman forces,the nonadiabatic coupling matrix elements 
! (v.d), and d[alpha,beta].
!====================================================================

      implicit none
      include 'param.h'
            
      real*8 Csp,Cse,pi
      real*8 R(nat),Rp0,Re0,qp,qe,qp0,qe0
      integer i,j,k,ad,ni,nj,lp,le,mp,me
      integer info,lwork,itype
      character jobz,uplo
      parameter(jobz='V')
      parameter(uplo='L')
      parameter(lwork=13600)
      parameter(itype=1)
      parameter(pi=3.14159d0)
      real*8 Kin(n,n),Vpes1(n,n),Vpes2(n,n),Vpes(n,n),S(n,n),Pot(n,n)
      real*8 Hmat(n,n),Smat(n,n)
      real*8 work(lwork)
      real*8 a,b,c,d
      real*8 hp,he
      real*8 phip,phie
      real*8 elow1,elow2,elow3,w(n),vect(n,n),vlow(n,n), &
      vect_old(n,n),vdot(n)
      DOUBLE PRECISION HF_1,HF_2,HF(nlv),d_tot(nlv,nlv)
      real*8 sum,adiab,adiabm,adiabn,adiabj,adiabk
     
! Set complex potential parameters
      Csp=2.d-3
      Cse=2.d-3
      Rp0=0.d0
      Re0=0.d0
      qp0=3.d0
      qe0=3.d0
           
! Integration parameters 
!      a=1.8d0
!      b=4.2d0
      a=2.d0
      b=4.d0
      c=-6.d0
      d=12.d0
      mp=100
      me=700
      hp=(b-a)/dble(mp)
      he=(d-c)/dble(me)

! Calculate Vpes matrix elements
! Take into account that negative sign of Vpes is already applied in
! Vpes1 and Vpes2 though building the input matrix.
      do i=1,n
        do j=1,n
          Vpes(i,j)=(Csp*(R(1)-Rp0)*Vpes1(i,j))+(Cse*(R(1)-Re0)*Vpes2(i,j))
          Hmat(i,j)=Kin(i,j)+Pot(i,j)+Vpes(i,j)
!          Hmat(i,j)=Kin(i,j)+Vpes(i,j)
!          Hmat(i,j)=Kin(i,j)+Pot(i,j)
!        Hmat(i,j)=Vpes(i,j)
        Smat(i,j)=S(i,j)
        enddo
      enddo
       
! Solve eigenvalue problem Hc=eSc using DSYGV subroutine
      call dsygv(itype,jobz,uplo,n,Hmat,n,Smat,n,w,work,lwork,info)

      if (info.ne.0) then
        write ( *, * ) 'Warning!'
        write ( *, * ) ' The error return flag INFO = ',info
      end if
! Choose the 3 lowest eigenvalues

      elow1=w(1)
      elow2=w(2)
      elow3=w(3)

!      do j=1,20
!        vdot(j)=0
!         do i=1,n
!            vdot(j)=vdot(j)+Hmat(i,1)*Hmat(i,j)
!         enddo
!         write(222,*) j,vdot(j)
!      enddo
!      stop
!-----------------------------------------------------------------
! Define first nlv eigenvectors
      do j=1,nlv
        do i=1,n
          vect_old(i,j)=vect(i,j)
          vect(i,j)=Hmat(i,j)
        enddo
      enddo
! Check for inversion of all signs in components of eigenvectors
! relative to those of previous time-step
      do j=1,nlv
        vdot(j)=0
      enddo
      do j=1,nlv
         do i=1,n
            vdot(j)=vdot(j)+vect_old(i,j)*vect(i,j)
         enddo
         if(vdot(j).lt.0) then
             do i=1,n
               vect(i,j)=-vect(i,j)
             enddo
         endif
      enddo      
!--------------------------------------------------------------------
! Define the eigenvector corresponding to the first nlv eigenvalues
      do j=1,nlv
        do i=1,n
          vlow(i,j)=vect(i,j)
        enddo  
      enddo
          
!     do k=1,n 
!     sum=0.d0
!     do i=1,n
!       do j=1,n
!         sum=sum+vect(i,k)*vect(j,k)*S(i,j)
!       enddo  
!     enddo
!     write(222,*) sum
!     enddo
!----------------------------------------------------------------
! Calculate the Hellmann-Feynman force on R coordinate.
!   Calculate integral(vlow*(dV/dR)*vlow) using trapezoid rule.   

      do ad=1,nlv
         sum=0.d0
         do i=1,n
           do j=1,n
              sum=sum+vect(i,ad)*vect(j,ad)*(Csp*Vpes1(i,j)+Cse*Vpes2(i,j))
           enddo  
         enddo
         HF(ad)=-sum
      enddo  !End loop over adiabatic states

! Calculate nonadiabatic coupling element, d[alpha,beta]
          
         sum=0.d0
         do i=1,n
           do j=1,n
              sum=sum+vect(i,1)*vect(j,2)*(Csp*Vpes1(i,j)+Cse*Vpes2(i,j))
           enddo  
         enddo
         d_tot(mm,nn)=-sum/(elow1-elow2)

      END

