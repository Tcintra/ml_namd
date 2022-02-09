      SUBROUTINE RUNNVT_NHC(neqs,HF,qs,Vs,Fs,Qb,Vb,Fb,Wb,c, &
                 p,p_e,p_s,d_c,w,Hmat,vect,d,IA,JB,beta)
      implicit none
      include 'param.h'

      integer it,i,j,k,ia,jb,ialpha,jbeta,neqs,nw,nc,nsy,res,nres
      parameter (nc=2)
      parameter (nsy=7)
      parameter (nres=1)
      real*8 ms,ws,s1,s2,mb,wc,Wb(mode),c(mode),coup
      REAL*8 qs(nat),Vs(nat),Fs(nat)
      REAL*8 Qb(mode),Vb(mode),Fb(mode)
      real*8 dt,dteq,dt2,dt4,dt8,ttt,beta
      real*8 enkib,potb,hamilb,enkis,EM_R,enki,hamil
      real*8 hameta,hampeta,hameta_1,hameta_rest
      real*8 meta1,meta(nc),eta(nc),peta(nc),weight(nsy)
      REAL*8 p(n,n),p_e(n,n),p_s(n,n),d_c(n,n),d(n,n)
      real*8 HF(nat,nlv),w(n),vect(n,n),Hmat(n,n)
      include 'sol_bath.in'
      include 'NHC.in'

!     Parameters for NH dynamics
      dteq=50.d0
      meta(1)=(nat+mode)*(20*dteq)**2/beta
      do j=2,nc
       meta(j)=(20*dteq)**2/beta
      enddo
      do j=1,nc
       eta(j)=0.d0
       peta(j)=0.d0
      enddo
      dt2=0.5d0*dteq
      dt4=0.25d0*dteq
      dt8=0.125d0*dteq
     
      DO I=1,NAT
         Fs(I)=Fs(I)+0.5d0*(HF(i,ia)+HF(i,ia))
      enddo
 
!  BEGIN MOLECULAR DYNAMICS CYCLE
      
      DO it=1,neqs

!Nose-Hoover chain acting on position and momentum of the solvent and
!bath.

       DO nw=nsy,1,-1

       do res=1,nres

        enkis = 0.d0
        DO k=1,NAT
           ENKIS = ENKIS + ms*Vs(k)**2
        ENDDO

        enkib = 0.d0
        DO I=1,mode
           ENKIb = ENKIb+ Mb*Vb(I)**2
        ENDDO

        enki=enkib+enkis
       
        peta(nc) = peta(nc) + (weight(nw)*DT4/nres)*((peta(nc-1)**2/meta(nc-1))-1/beta)
         
        do j=nc-1,1,-1
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
           
           if (j.eq.1) then
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*(enki-(nat+mode)/beta)
           else
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*((peta(j-1)**2/meta(j-1))-1/beta)
           endif
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
        enddo      
  
        do j=nc,1,-1
           eta(j) = eta(j) - (weight(nw)*DT2/nres)*peta(j)/meta(j)
        enddo

        DO k=NAT,1,-1
           Vs(k) = Vs(k)*exp((-weight(nw)*DT2/nres)*peta(1)/meta(1))
        ENDDO
        
        DO I=mode,1,-1
           Vb(i) = Vb(i)*exp((-weight(nw)*DT2/nres)*peta(1)/meta(1))
        ENDDO
        
       do j=1,nc-1
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
           if (j.eq.1) then
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*(enki-(nat+mode)/beta)
           else
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*((peta(j-1)**2/meta(j-1))-1/beta)
           endif
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
       enddo      

       peta(nc) = peta(nc) + (weight(nw)*DT4/nres)*((peta(nc-1)**2/meta(nc-1))-1/beta)
      
       ENDDO
 
       ENDDO
       
! Velocity verlet algorithm

      DO k=1,nat
         Vs(k) = Vs(k) + DT2*Fs(k)/ms
         qs(k) = qs(k) + DTEQ*Vs(k)
      ENDDO

      DO I=1,mode
         Vb(I) = Vb(I) + DT2*(Fb(I)/Mb)
         Qb(i) = Qb(I) + DTEQ*Vb(I)
      ENDDO

      CALL FORCE(qs,c,Wb,Qb,Fb,Fs)
      CALL ADIAB_NVE(p,p_e,p_s,d_c,w,Hmat,vect,HF,d,qs,c,Wb,Qb,Vs,Vb)

      DO k=1,nat
         Fs(k)=Fs(k)+0.5d0*(HF(k,ia)+HF(k,ia))
         Vs(k) = Vs(k) + DT2*Fs(k)/ms
      ENDDO

       DO I=1,mode
         Vb(I) = Vb(I) + DT2*Fb(I)/Mb
       ENDDO

!Second set of Nose-Hoover chain acting on position and momentum of the solvent and
!bath.

       DO nw=nsy,1,-1

       do res=1,nres

        enkis = 0.d0
        DO k=1,NAT
           ENKIS = ENKIS + ms*Vs(k)**2
        ENDDO

        enkib = 0.d0
        DO I=1,mode
           ENKIb = ENKIb+ Mb*Vb(I)**2
        ENDDO

        enki=enkib+enkis
       
        peta(nc) = peta(nc) + (weight(nw)*DT4/nres)*((peta(nc-1)**2/meta(nc-1))-1/beta)
   
        do j=nc-1,1,-1
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
           if (j.eq.1) then
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*(enki-(nat+mode)/beta)
           else
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*((peta(j-1)**2/meta(j-1))-1/beta)
           endif
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
        enddo      
  
        do j=nc,1,-1
           eta(j) = eta(j) - (weight(nw)*DT2/nres)*peta(j)/meta(j)
        enddo

        DO k=NAT,1,-1
           Vs(k) = Vs(k)*exp((-weight(nw)*DT2/nres)*peta(1)/meta(1))
        ENDDO
        
        DO I=mode,1,-1
           Vb(i) = Vb(i)*exp((-weight(nw)*DT2/nres)*peta(1)/meta(1))
        ENDDO
        
       do j=1,nc-1
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
           if (j.eq.1) then
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*(enki-(nat+mode)/beta)
           else
           peta(j) = peta(j) + (weight(nw)*DT4/nres)*((peta(j-1)**2/meta(j-1))-1/beta)
           endif
           peta(j) = peta(j)*exp((-weight(nw)*dt8/nres)*peta(j+1)/meta(j+1))
       enddo      

       peta(nc) = peta(nc) + (weight(nw)*DT4/nres)*((peta(nc-1)**2/meta(nc-1))-1/beta)
      
       ENDDO
 
       ENDDO

       enkis = 0.d0
       DO k=1,NAT
          ENKIS = ENKIS + ms*Vs(k)**2
       ENDDO

       enkib = 0.d0
       DO I=1,mode
          ENKIb = ENKIb+ Mb*Vb(I)**2
       ENDDO
 
       enki=enkib+enkis
 
       POTb = 0.d0
       DO I=1,mode
         POTb = POTb + 0.5d0*Mb*Wb(I)**2*(Qb(I)-(c(I)*qs(nat)/(Mb*Wb(I)**2)))**2
       ENDDO

      
       hameta_1=dble(nat+mode)*eta(1)/beta

       hameta_rest=0.d0
       do j=2,nc
         hameta_rest=hameta_rest+eta(j)/beta
       enddo

       hameta=hameta_1+hameta_rest

       hampeta=0.d0
       do j=1,nc
         hampeta=hampeta+0.5d0*peta(j)**2/meta(j)
       enddo

       enki=0.5d0*enki
       enkis=0.5d0*enkis
       enkib=0.5d0*enkib
       hamilb=enkib+potb
       TTT=2.d0*ENKI/(KB*(NAT+mode))
       EM_R=0.5d0*(W(ia)+W(ia))
       hamil=hamilb+enkis+EM_R+hameta+hampeta

!      write(222,'(i7,8g15.7)') it,hamil,em_r,enkib,potb,enkis,hameta,hampeta,ttt
!      write(222,'(i7,3g15.7)') it,hamil,qs(1),ttt

!      if(mod(it,iprint).eq.0) then
!        write(50,25) it,hamil,enki,pot,hameta,hampeta,ttt,ptot
!      endif
       
      ENDDO

! 25   format(1x,i7,1x,7g15.7)

      RETURN
      END

