      SUBROUTINE RUNNVE(DT,nmds,IPRINT,TOLNCE,MASS, &
                 BOXL,R,V,F,PI,CHARGE,HF,ELOW1,IT,SMAT,VS, &
                 ALPHA,BETA,elow2,elow3,mcs,vpes1,vpes2,Kin,Pot, &
                 de_bar,nb,flux,unbias,d_tot,dedt,wphi, &
                 dcoup,S,hbar,IALPHA,JBETA,WPHI_FACTOR, &
                 IDUM, P_FACT,NJUMP,S_SHIFT_CUT,LREJECT1, &
                 PI_WEIGHT,IA,JB,GNORM,OBS,POP,RHOW,ITBLOCK, &
                 NJUMP_CUT,NSLICE,NJUMPM,NJP_STAT,STAT_WEIGHT, &
                 PHASE_FACT_RE,PHASE_FACT_IM,STAT_SLICE, &
                 WEIGHT_SLICE,RHONORM,qmat,lreject2,kappa, &
                 nrej1,nrej2,nrej3,d)
      implicit none
      include 'param.h'
      INTEGER IALPHA,JBETA,nrej1,nrej2,nrej3(nmdsbx)
      integer IA,JB,ITBLOCK,NJUMP_CUT,mcs,NJUMP,IB
      integer NJUMPM(NMDSBX,nlv,nlv),NJP_STAT(NMDSBX,NJC_MAX)
      integer nmds,nslice,IPRINT,IDUM,it,i,j,k,count
      REAL*8 R(NAT),V(NAT)
      REAL*8 MASS(NAT),CHARGE(NAT)
      real*8 dt,dt2,boxl,tolnce,elow1
      real*8 hamil,ttt,enki,ptot,pi,vs
      real*8 pc,beta,kti,elow2,elow3
      REAL*8 alpha
      real*8 Vcomp,Kin(n,n),Pot(n,n),S(n,n,nat),Hmat(n,n),Smat(n,n)
      real*8 Vpes1(n,n),Vpes2(n,n),adiab,adiabn,adiabm
      real*8 vlow(n,n),vect(n,n),sum 
      real*8 C1,C2,C3,de
      real*8 de_bar,dEdt,w1
      real*8 Nb(NAT),flux(NAT),D,kappa(NAT),unbias
      real*8 GNORM(NMDSBX,nlv,nlv),STAT_WEIGHT(NMDSBX,nlv,nlv),POP(NMDSBX,nlv,nlv)
      real*8 rhonorm(NMDSBX,nlv,nlv),rhow(nlv,nlv)
      real*8 PHASE_FACT_RE(NMDSBX,nlv,nlv),PHASE_FACT_IM(NMDSBX,nlv,nlv)
      real*8 HBAR,stat_weight_fact,PHASE_RE,PHASE_IM
      real*8 WPHI(nlv,nlv),omega(nlv,nlv),E_R(n),EM_R,obs(nlv,nlv)
      real*8 FWij(nat),F(nat),HF(nlv)
      real*8 WS,R0
      real*8 d_tot(nlv,nlv),dcoup(nlv,nlv),listcon
      real*8 dnorm(nlv,nlv),dcoup_norm(nlv,nlv),vd_norm(nlv,nlv)
      real*8 sumd,sumdp,sumdn
      real*8 qmat(nlv,nlv)
      real*8 P_FACT,S_SHIFT_CUT,WEIGHT_SLICE,PI_WEIGHT
      real*8 check,flag,flag1,flag2,flag3
      complex*16 STAT_SLICE(NMDSBX,nlv,nlv),stat_slice_fact,WPHI_FACTOR
      LOGICAL LREJECT1,LREJECT2,NBCUT
      include 'solvent.in'

      count=0
      DT2=0.5d0*DT
      
      DO I=1,NAT
        FWij(I)=F(I)+0.5d0*(HF(IA)+HF(JB))
      END DO
   
! Initialization of variables for calculation of rate
!     if (mcs.eq.1) then
!       unbias=0.d0
!       do i=1,nmds
!         flux(i)=0.d0
!         kappa(i)=0.d0
!       enddo
!     endif
!     do i=1,nmds
!       nb(i)=0.d0
!     enddo

!======================================================================
! RATE AT T=0
!    OPEN(332,FILE='rate.out',FORM='formatted',STATUS='unknown')
      it=1

!   Calculate v.d
      dcoup(mm,nn)=0.d0
      do k=1,nat
        dcoup(mm,nn)=dcoup(mm,nn)+v(k)*d_tot(mm,nn)
      enddo
      dcoup(mm,nn)=dcoup(mm,nn)
      dcoup(nn,mm)=-dcoup(mm,nn)
           
!     CALL HEAVISIDE(dE,dEdt,it,nmds,mcs,de_bar,Nb,tolnce)

! Initialization of variables required in the calculation of the 
! observable

      IB=1
      P_FACT=1.D0
      PI_WEIGHT=1.D0
      STAT_WEIGHT_FACT=1.D0
      STAT_SLICE_FACT=cmplx(1.D0,0.d0)
      WPHI_FACTOR=CMPLX(1.D0,0.D0)
      FLAG=1.d0
      FLAG1=1.d0
      FLAG2=1.d0
      NBCUT=.FALSE.
      CALL OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,OBS,nmds,Nb, &
            IA,JB,POP,IALPHA,JBETA,RHOW,IB,qmat,rhonorm,nbcut)
         
!     CALL RATE(dEdt,it,Nb,flux,D,kappa,w1,flag)
      GNORM(IB,IA,JB)=GNORM(IB,IA,JB)+1.D0  !!! normalization
      
!   Weight of the quantum transitions
      STAT_WEIGHT(IB,IA,JB)=STAT_WEIGHT(IB,IA,JB)+STAT_WEIGHT_FACT
      PHASE_RE=DBLE(WPHI_FACTOR)
      PHASE_IM=IMAG(WPHI_FACTOR)
      PHASE_FACT_RE(IB,IA,JB)=PHASE_FACT_RE(IB,IA,JB)+PHASE_RE
      PHASE_FACT_IM(IB,IA,JB)=PHASE_FACT_IM(IB,IA,JB)+PHASE_IM

!   Weight of the quantum transitions * q-phase
      STAT_SLICE(IB,IA,JB)=STAT_SLICE(IB,IA,JB)+STAT_SLICE_FACT
      NJUMPM(IB,IA,JB)=0
!==================================================================
      
!  BEGIN MOLECULAR DYNAMICS CYCLE

      DO it=2,nmds
        count=1+count
        weight_slice=1.d0
        flag3=1
!       if((flag1.eq.0).OR.(flag2.eq.0)) go to 129

! Velocity verlet algorithm
      DO I=1,NAT
         V(I) = V(I) + DT2*(FWij(I))/MASS(I)
         R(I) = R(I) + DT*V(I)
      ENDDO

!=============== Termination of trajectory =================
      if ((it.gt.10).and.(abs(R(1)).gt.0.4).AND.(IALPHA.eq.1).AND.(JBETA.eq.1)) then
      go to 100
      end if
!====================================================================== 
       CALL FORCE(R,VS,F)
  
       CALL ADIAB_NVE(Kin,Pot,Vpes1,Vpes2,S,adiab,adiabm,adiabn, &
                 Hmat,Smat,vect,vlow,elow1,elow2,elow3,HF,d_tot,r)

!      CALL POLARIZ(NAT,BOXL,X,Y,Z,CHARGE,RCMX,DE,
!    &              RCMY,RCMZ,RRELX,RRELY,RRELZ,rcut)
       E_R(1)=elow1
       E_R(2)=elow2
       E_R(3)=elow3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Wphi  QUANTUM ADIABATIC PHASE         !!!
!!           ALONG THE BATH TRAJECTORY            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       OMEGA(IALPHA,JBETA)=(E_R(IALPHA)-E_R(JBETA))/HBAR
       WPHI(IALPHA,JBETA)=OMEGA(IALPHA,JBETA)*DT
       WPHI_FACTOR=WPHI_FACTOR &
       *CMPLX(COS(WPHI(IALPHA,JBETA)),SIN(WPHI(IALPHA,JBETA)))

!     Finish integrating the equations of motion

       DO I=1,NAT
         FWij(I)=F(I)+0.5d0*(HF(IALPHA)+HF(JBETA))
         V(I) = V(I) + DT2*(FWij(I))/MASS(I)
       ENDDO

! Calculate the velocity/unnormalized nonadiabatic coupling dot product


       dcoup(mm,nn)=0.d0
       do k=1,nat
         dcoup(mm,nn)=dcoup(mm,nn)+v(k)*d_tot(mm,nn)
       enddo
      
!      dcoup(mm,nn)=dcoup(mm,nn)*flag
!      dcoup(nn,mm)=-dcoup(mm,nn)*flag
       dcoup(mm,nn)=dcoup(mm,nn)
       dcoup(nn,mm)=-dcoup(mm,nn)
!       write(222,*) mm,nn,dcoup(mm,nn)
!       Nonadiabatic dynamics

! Calculate the velocity/normalized nonadiabatic coupling dot product and 
! normalized nonadiabatic coupling vector 
       dnorm(mm,nn)=0.d0
       do k=1,nat
         dnorm(mm,nn)=dnorm(mm,nn)+d_tot(mm,nn)**2
       enddo

       do k=1,nat
         dnorm(mm,nn)=d_tot(mm,nn)/dsqrt(dnorm(mm,nn))
!        write(142,'(2i6,f20.10)') it,k,dnormx(mm,nn,k)
       enddo
       
!      vd_norm(mm,nn)=0.d0
!      do k=1,nat
!        vd_norm(mm,nn)=vd_norm(mm,nn)+vx(k)*dnormx(mm,nn,k)+
!    &                  vy(k)*dnormy(mm,nn,k)+vz(k)*dnormz(mm,nn,k)
!      enddo
      
       dcoup_norm(mm,nn)=0.d0
       do k=1,nat
         dcoup_norm(mm,nn)=dcoup_norm(mm,nn)+v(k)*dnorm(mm,nn)
       enddo
       dcoup_norm(nn,mm)=-dcoup_norm(mm,nn)

! Calculate velocity shifts

      sumd=0.d0
      do k=1,nat
        sumd=sumd+(dnorm(mm,nn)**2)/mass(k)
      enddo

      sumdn=-dcoup_norm(mm,nn)-dsqrt(dcoup_norm(mm,nn)**2+(elow1-elow2)* &
      sumd)
      sumdp=-dcoup_norm(mm,nn)+dsqrt(dcoup_norm(mm,nn)**2+(elow1-elow2)* &
      sumd)
      
!     do mm=1,nlv
!      do nn=1,nlv
        do k=1,nat
          if(dcoup_norm(mm,nn).lt.0) then
            S(mm,nn,k)=dnorm(mm,nn)*sumdn/(mass(k)*sumd)
          endif
          if(dcoup_norm(mm,nn).gt.0) then
            S(mm,nn,k)=dnorm(mm,nn)*sumdp/(mass(k)*sumd)
          endif
        enddo

      sumdn=-dcoup_norm(mm,nn)-dsqrt(dcoup_norm(mm,nn)**2+(elow2-elow1)* &
      sumd*c3)
      sumdp=-dcoup_norm(mm,nn)+dsqrt(dcoup_norm(mm,nn)**2+(elow2-elow1)* &
      sumd*c3)
     
        do k=1,nat
          if(dcoup_norm(mm,nn).lt.0) then
            S(nn,mm,k)=dnorm(mm,nn)*sumdn/(mass(k)*sumd)
          endif
          if(dcoup_norm(mm,nn).gt.0) then
            S(nn,mm,k)=dnorm(mm,nn)*sumdp/(mass(k)*sumd)
          endif
        enddo
!      enddo
!     enddo

      check=dcoup_norm(mm,nn)**2/sumd
!=============== QUANTUM JUMP ========================================
!     IF((DE.GE.0.01).AND.(DE.LE.0.016)) THEN
      IF(MOD(IT,NSLICE).EQ.0) THEN
        CALL MAKEJUMP(IDUM,IALPHA,JBETA,MASS,V,S,P_FACT, &
                      NJUMP,S_SHIFT_CUT,CHECK,E_R,LREJECT1, &
                      LREJECT2,DT,PI_WEIGHT,NSLICE, &
                      WEIGHT_SLICE,DCOUP,tolnce, &
                      boxl,R,nrej1,nrej2,it)
        IF(NJUMP.GT.NJUMP_CUT) THEN
!           NJUMP=NJUMP-1
            FLAG=0.d0
!           nrej3(it)=nrej3(it)+1
        END IF

        IF(LREJECT1) GO TO 888  !!! Continue trajectory adiabatically
      END IF

!======================================================================
 888 continue

 
         enki = 0.d0
         DO I=1,NAT
            ENKI = ENKI+MASS(I)*(V(I)**2)
         ENDDO

      ENKI = 0.5D0*ENKI
      TTT=(2.D0/DBLE(3.D0*NAT))*ENKI
      EM_R=0.5d0*(E_R(ialpha)+E_R(jbeta))
      hamil=enki+Vs+EM_R

      go to 122      
! ================ CALCULATION OF RATE ==========================
!     IF(MOD(IT,ITBLOCK).EQ.0) THEN
        IB=IT/ITBLOCK
!     Impose sink 
        if ((de.lt.0.0065).AND.(IALPHA.EQ.1).AND.(JBETA.EQ.1)) then
          flag1=0
        endif
        if ((de.gt.0.0222).AND.(IALPHA.EQ.1).AND.(JBETA.EQ.1)) then
          flag2=0
        endif
        IF(((IALPHA.EQ.1).AND.(JBETA.EQ.1)).OR.((IALPHA.EQ.2).AND. &
      (JBETA.EQ.2))) THEN
          if ((flag1.eq.1).AND.(flag2.eq.1)) then
            CALL HEAVISIDE(dE,dEdt,ib,nmds,mcs,de_bar,Nb,tolnce)
          endif

          if (flag1.eq.0) then
            Nb(it)=0
            flag3=0
          endif
          if (flag2.eq.0) then
            Nb(it)=1
            flag3=0
          endif
          CALL OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,OBS,nmds,Nb, &
                IA,JB,POP,IALPHA,JBETA,RHOW,IB,qmat,rhonorm,nbcut)
!          CALL RATE(dEdt,it,Nb,flux,D,kappa,w1,flag)
        ENDIF

122    continue

        GNORM(IT,IA,JB)=GNORM(IT,IA,JB)+1.d0   !normalization

! Accumulated weights along a trajectory
        STAT_WEIGHT_FACT=STAT_WEIGHT_FACT*WEIGHT_SLICE 

      if(mod(it,iprint).eq.0) then
!       write(60,25) it,hamil,enki,vs,r(1),elow1,elow2,dcoup(mm,nn),ialpha,jbeta
        write(60,25) it,hamil,enki,vs,r(1),elow1,elow2,dcoup(mm,nn),ialpha,jbeta
      endif

!       if(it.eq.nmds) then
!         do i=1,nmds
!           write(332,*) i,kappa(i)/unbias
!         enddo
!       endif

        NJUMPM(IT,IA,JB)=NJUMPM(IT,IA,JB)+NJUMP
        IF(NJUMP.LE.NJC_MAX) THEN
           NJP_STAT(IT,NJUMP+1)=NJP_STAT(IT,NJUMP+1)+1
        END IF

      ENDDO  ! end of time loop
100  continue    
25   format(1x,i7,1x,7g15.7,2i5)
      RETURN
      END



