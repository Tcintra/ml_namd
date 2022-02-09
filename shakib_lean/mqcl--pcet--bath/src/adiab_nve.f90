      SUBROUTINE ADIAB_NVE(p,p_e,p_s,d_c,w,Hmat,vect,HF,d,qs,c,Wb,Qb,Vs,Vb) 

!====================================================================
! This subroutine will calculate the adiabatic states,the
! Hellmann-Feynman forces,the nonadiabatic coupling matrix elements 
! (v.d), and d[alpha,beta].
!====================================================================

      implicit none
      include 'param.h'
           
      real*8 qs(nat),s1,s2,C(mode),Wb(mode),Qb(mode),Vs(nat),Vb(mode)
      real*8 ms,ws,mb,wc,coup 
      integer i,j,l,li,lj,ad,k
      integer info,lwork,itype
      character jobz,uplo
      parameter(jobz='V')
      parameter(uplo='L')
      parameter(lwork=1280)
      real*8 p(n,n),p_e(n,n),p_e_11(n,n),p_e_22(n,n),p_s(n,n),d_c(n,n),p_s_1(n,n)
      real*8 Hmat(n,n)
      real*8 work(lwork)
      real*8 qp,np,lp,hp,a,b,phi,mu1
      real*8 w(n),vect(n,n),vect_old(n,n),vdot(n)
      DOUBLE PRECISION HF(nat,nlv),d(nlv,nlv,nat)
      real*8 sum
      include 'sol_bath.in'
     
! Calculate mixed quantum_classic matrix elements

      DO k=1,nat
      DO i=1,n
        DO j=1,n
          p_s_1(i,j)=p_s(i,j)*qs(k)
          if (j.le.0.5d0*n) then  
            if (i.eq.j) then    
              p_e_11(i,j)=p_e(i,j)+0.5d0*(ms*ws**2.d0)*(qs(k)-s1)**2.d0
            else
              p_e_11(i,j)=p_e(i,j)
            endif
          Hmat(i,j)=p(i,j)+p_e_11(i,j)+d_c(i,j)+p_s_1(i,j)
          else
            if (i.eq.j) then    
              p_e_22(i,j)=-p_e(i,j)+0.5d0*(ms*ws**2.d0)*(qs(k)-s2)**2.d0
            else
              p_e_22(i,j)=-p_e(i,j)
            endif
          Hmat(i,j)=p(i,j)+p_e_22(i,j)+d_c(i,j)+p_s_1(i,j)
          endif
        ENDDO
      ENDDO
      ENDDO
      
! Solve eigenvalue problem Hc=ec using DSYEV subroutine
      call dsyev(jobz,uplo,n,Hmat,n,w,work,lwork,info)

      if (info.ne.0) then
        write ( *, * ) 'Warning!'
        write ( *, * ) ' The error return flag INFO = ',info
      end if
!---------------------------------------------------------------
! Define first nlv eigenvectors.
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
! Calculate the Hellmann-Feynman force on R coordinate.
!   Calculate integral(vlow*(dV/dR)*vlow) using trapezoid rule.   

         do k=1,nat
         do ad=1,nlv

         sum=0.d0

         do i=1,n-0.5*n
            sum=sum+vect(i,ad)**2*(ms*ws**2*(qs(k)-s1))
         enddo

         do i=0.5*n+1,n
            sum=sum+vect(i,ad)**2*(ms*ws**2*(qs(k)-s2))
         enddo

            do i=1,n-0.5*n
              do j=1,n-0.5*n
                sum=sum+vect(i,ad)*vect(j,ad)*p_s(i,j)
              enddo
            enddo

            do i=0.5*n+1,n
              do j=0.5*n+1,n
                sum=sum+vect(i,ad)*vect(j,ad)*p_s(i,j)
              enddo
            enddo

           HF(k,ad)=-sum
         enddo 
         enddo

      !End loop over adiabatic states

! Calculate nonadiabatic coupling element, d[alpha,beta]
         
         do k=1,nat
         do mm=1,nlv
           do nn=1,nlv
             if(mm.ne.nn) then
               sum=0.d0
               do i=1,n-0.5*n
                  sum=sum+vect(i,mm)*vect(i,nn)*(ms*ws**2*(qs(k)-s1))
               enddo

               do i=0.5*n+1,n
                  sum=sum+vect(i,mm)*vect(i,nn)*(ms*ws**2*(qs(k)-s2))
               enddo

               do i=1,n-0.5*n
                 do j=1,n-0.5*n
                   sum=sum+vect(i,mm)*vect(j,nn)*p_s(i,j)
                 enddo
               enddo

               do i=0.5*n+1,n
                 do j=0.5*n+1,n
                   sum=sum+vect(i,mm)*vect(j,nn)*p_s(i,j)
                 enddo
               enddo
               d(mm,nn,k)=-sum/(w(mm)-w(nn))
             endif
           enddo
         enddo
         enddo
      END
