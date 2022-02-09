      SUBROUTINE RUNNVE(DT,nmds,IPRINT,ms,w,olA,olB,sp_density, &
                 Qb,Vs,Vb,Fs,Fb,PI,HF,IT,vect,sample, &
                 ALPHA,BETA,mcs,p,p_e,p_s,d_c,wphi, &
                 dcoup,S,qs,c,Wb,IALPHA,JBETA,WPHI_FACTOR, &
                 IDUM,P_FACT,NJUMP,LREJECT1, &
                 PI_WEIGHT,IA,JB,GNORM,OBS,RATE,RHOW,ITBLOCK, &
                 NJUMP_CUT,NSLICE,NJUMPM,NJP_STAT,STAT_WEIGHT, &
                 PHASE_FACT_RE,PHASE_FACT_IM,STAT_SLICE, &
                 WEIGHT_SLICE,RHONORM,lreject2, &
                 nrej1,nrej2,nrej3,d,ratea,rateb,ratec,rated, &
                 ratee,ratef,rateg,rateh)
      implicit none
      include 'param.h'
      INTEGER IALPHA,JBETA,nrej1,nrej2,nrej3(nmdsbx)
      integer IA,JB,ITBLOCK,NJUMP_CUT,mcs,NJUMP,IB
      integer NJUMPM(NMDSBX,nlv,nlv),NJP_STAT(NMDSBX,NJC_MAX)
      integer nmds,nslice,IPRINT,IDUM,it,i,j,k,count
      REAL*8 w(n),norm(nlv,nlv),ws,wc,wb(mode),Mb,s1,s2,coup
      REAL*8 qs(NAT),Vs(NAT),Qb(mode),Vb(mode),c(mode)
      REAL*8 ms,CHARGE(NAT),temp,sp_density
      real*8 dt,dt2,boxl,tolnce,elow1
      real*8 hamil,hamilb,ttt,enkis,enkib,potb,ptot,pi,enki
      real*8 pc,beta,kti
      REAL*8 alpha,sample
      real*8 p(n,n),p_e(n,n),Hmat(n,n)
      real*8 p_s(n,n),d_c(n,n)
      real*8 vect(n,n),sum,olA(0.5*n,0.5*n),olB(0.5*n,0.5*n) 
      real*8 C1,C2,C3,de
      real*8 Nb(NAT),flux(NAT),kappa(NAT),unbias
      real*8 GNORM(NMDSBX,nlv,nlv),STAT_WEIGHT(NMDSBX,nlv,nlv), &
             RATE(NMDSBX,nlv,nlv),OBSERV_IM(nmdsbx,nlv,nlv)
      REAL*8 RATEa(NMDSBX,nlv,nlv),RATEb(NMDSBX,nlv,nlv),RATEc(NMDSBX,nlv,nlv), &
             RATEd(NMDSBX,nlv,nlv)
      REAL*8 RATEe(NMDSBX,nlv,nlv),RATEf(NMDSBX,nlv,nlv),RATEg(NMDSBX,nlv,nlv), &
             RATEh(NMDSBX,nlv,nlv)
      real*8 rhonorm(NMDSBX,nlv,nlv),rhow(nlv,nlv)
      real*8 PHASE_FACT_RE(NMDSBX,nlv,nlv),PHASE_FACT_IM(NMDSBX,nlv,nlv)
      real*8 stat_weight_fact,PHASE_RE,PHASE_IM
      real*8 WPHI(nlv,nlv),omega(nlv,nlv),EM_R,obs(nlv,nlv)
      real*8 FWij(nat),Fs(nat),Fb(nat),HF(nat,nlv)
      real*8 d(nlv,nlv,nat),dcoup(nlv,nlv),S(nlv,nlv,nat),listcon
      real*8 dnorm(nlv,nlv,nat),dcoup_norm(nlv,nlv),vd_norm(nlv,nlv)
      real*8 sumd,sumdp,sumdn
      real*8 qmat(nlv,nlv)
      real*8 P_FACT,S_SHIFT_CUT,WEIGHT_SLICE,PI_WEIGHT
      real*8 check(nlv,nlv),flag,flag1,flag2
      complex*16 STAT_SLICE(NMDSBX,nlv,nlv),stat_slice_fact,WPHI_FACTOR
      LOGICAL LREJECT1,LREJECT2
      include 'sol_bath.in'

      count=0
      DT2=0.5d0*DT
      
      DO k=1,nat
         Fs(k)=Fs(k)+0.5d0*(HF(k,IA)+HF(k,JB))
      ENDDO

!======================================================================
      it=1

! Calculate v.d
      do mm=1,nlv
        do nn=1,nlv
            dcoup(mm,nn)=0.d0
        enddo
      enddo
      do mm=1,nlv
        do nn=1,nlv
         if(mm.ne.nn) then
          do k=1,nat
            dcoup(mm,nn)=dcoup(mm,nn)+Vs(k)*d(mm,nn,k)
          enddo
         endif
        enddo
      enddo     
!=============================================================== 
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
      CALL OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,rate,nmds, &
                 IA,JB,sp_density,olB,IALPHA,JBETA,IB,vect, &
                 ratea,rateb,ratec,rated,ratee,ratef,rateg,rateh)
         
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

! Velocity verlet algorithm

      DO k=1,nat
         Vs(k) = Vs(k) + DT2*Fs(k)/ms
         qs(k) = qs(k) + DT*Vs(k)
      ENDDO
    
      DO I=1,mode
         Vb(I) = Vb(I) + DT2*(Fb(I))/Mb
         Qb(i) = Qb(I) + DT*Vb(I)
      ENDDO
        
        CALL FORCE(qs,c,Wb,Qb,Fb,Fs)

        CALL ADIAB_NVE(p,p_e,p_s,d_c,w,Hmat,vect,HF,d,qs,c,Wb,Qb,Vs,Vb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Wphi  QUANTUM ADIABATIC PHASE         !!!
!!           ALONG THE BATH TRAJECTORY            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       OMEGA(IALPHA,JBETA)=(W(IALPHA)-W(JBETA))/HBAR
       WPHI(IALPHA,JBETA)=OMEGA(IALPHA,JBETA)*DT
       WPHI_FACTOR=WPHI_FACTOR &
       *CMPLX(COS(WPHI(IALPHA,JBETA)),SIN(WPHI(IALPHA,JBETA)))

!     Finish integrating the equations of motion

       DO k=1,nat
         Fs(k)=Fs(k)+0.5d0*(HF(k,IALPHA)+HF(k,JBETA))
         Vs(k) = Vs(k) + DT2*(Fs(k))/ms
       ENDDO

       DO I=1,mode
         Vb(I) = Vb(I) + DT2*(Fb(I))/Mb
       ENDDO


!    go to 121
!   Calculate v.d
      do mm=1,nlv
        do nn=1,nlv
         dcoup(mm,nn)=0.d0
         if(mm.ne.nn) then
          do k=1,nat
            dcoup(mm,nn)=dcoup(mm,nn)+Vs(k)*d(mm,nn,k)
          enddo
         endif
        enddo
      enddo
      
! Calculate the velocity/normalized nonadiabatic coupling dot product
! and 
! normalized nonadiabatic coupling vector 
      do mm=1,nlv
        do nn=1,nlv
         if(mm.ne.nn) then
           norm(mm,nn)=0.d0
           dcoup_norm(mm,nn)=0.d0
           do k=1,nat
             norm(mm,nn)=norm(mm,nn)+d(mm,nn,k)**2
          enddo
          do k=1,nat
            dnorm(mm,nn,k)=d(mm,nn,k)/dsqrt(norm(mm,nn))
          enddo
          sumd=0.d0
          do k=1,nat
            dcoup_norm(mm,nn)=dcoup_norm(mm,nn)+Vs(k)*dnorm(mm,nn,k)
            sumd=sumd+(dnorm(mm,nn,k)**2)/ms
          enddo
! Calculate velocity shifts

      temp=dsqrt(dcoup_norm(mm,nn)**2+(w(mm)-w(nn))*sumd)
      sumdn=(-dcoup_norm(mm,nn)-temp)/sumd
      sumdp=(-dcoup_norm(mm,nn)+temp)/sumd
        do k=1,nat
          if(dcoup_norm(mm,nn).lt.0) then
            temp = sumdn/ms
          endif
          if(dcoup_norm(mm,nn).gt.0) then
            temp = sumdp/ms
          endif
          S(mm,nn,k)=dnorm(mm,nn,k)*temp
        enddo

      check(mm,nn)=dcoup_norm(mm,nn)**2/sumd
      endif
      enddo
      enddo
!=============== QUANTUM JUMP ========================================
      IF(MOD(IT,NSLICE).EQ.0) THEN
        CALL MAKEJUMP(IDUM,IALPHA,JBETA,MS,VS,S,P_FACT, &
                      NJUMP,CHECK,W,LREJECT1, &
                      LREJECT2,DT,PI_WEIGHT,NSLICE, &
                      WEIGHT_SLICE,DCOUP, &
                      nrej1,nrej2,it)
        IF(NJUMP.GT.NJUMP_CUT) THEN
            NJUMP=NJUMP-1
            FLAG=0.d0
        END IF

      END IF
!======================================================================
!121   continue   
      enkis=0.d0
         DO k=1,nat
            enkis=enkis+0.5d0*ms*Vs(k)**2
         ENDDO

         enkib = 0.d0
         POTb = 0.d0
         DO I=1,mode
            ENKIb = ENKIb+0.5d0*Mb*Vb(I)**2
            POTb = POTb + 0.5d0*Mb*Wb(I)**2*(Qb(I)-(c(I)*qs(nat)/(Mb*Wb(I)**2)))**2
         ENDDO
         hamilb=enkib+potb

      ENKI=ENKIb+enkis
      TTT=2.D0*ENKI/DBLE(KB*(NAT+mode))
      EM_R=0.5d0*(W(ialpha)+W(jbeta))
      hamil=hamilb+enkis+EM_R
!     write(777,'(i7,3g15.7)') it,hamil,qs(1),ttt
      IB=IT/ITBLOCK
      CALL OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,rate,nmds, &
                 IA,JB,sp_density,olB,IALPHA,JBETA,IB,vect, &
                  ratea,rateb,ratec,rated,ratee,ratef,rateg,rateh)
      GNORM(IT,IA,JB)=GNORM(IT,IA,JB)+1.d0   !normalization

! Accumulated weights along a trajectory
      STAT_WEIGHT_FACT=STAT_WEIGHT_FACT*WEIGHT_SLICE 

        NJUMPM(IT,IA,JB)=NJUMPM(IT,IA,JB)+NJUMP
        IF(NJUMP.LE.NJC_MAX) THEN
           NJP_STAT(IT,NJUMP+1)=NJP_STAT(IT,NJUMP+1)+1
        END IF
      ENDDO  ! end of time loop

      RETURN
      END
