      SUBROUTINE MAKEJUMP(IDUM,IALPHA,JBETA,MASS,V,S,P_FACT, &
                          NJUMP,S_SHIFT_CUT,CHECK,E_R,LREJECT1, &
                          LREJECT2,DT,PI_WEIGHT,NSLICE, &
                          WEIGHT_SLICE,DCOUP,tolnce, &
                          boxl,R,nrej1,nrej2,it)

      IMPLICIT NONE
      INCLUDE 'param.h'
      INTEGER IDUM,IALPHA,JBETA,IALPHA_N,JBETA_N,nrej1,nrej2
      INTEGER MU,NU,IS,I,IPDECIDE,NJUMP,NSLICE
      INTEGER IT
      REAL*8 MASS(nat),R_EXTRA,RANF
      REAL*8 P_FACT,WEIGHT_SLICE,PI_WEIGHT
      REAL*8 PYES_J,PYES_J12,PYES_J21,YES_J12,YES_J21
      REAL*8 PNO_J,PNO_J12,PNO_J21,DTSLICE
      REAL*8 tolnce,boxl
      REAL*8 S_SHIFT_CUT,DT
      REAL*8 DCOUP(nlv,nlv),vdotd,e_r(n),check
      REAL*8 V(nat),c3,R(nat)
      REAL*8 S(n,n,nat),S_SHIFT(nat),P_SHIFT
      LOGICAL LREJECT1,LREJECT2

      DTSLICE=DBLE(NSLICE)*DT

      LREJECT1=.FALSE.
      LREJECT2=.FALSE.
      CALL RANDOM_NUMBER(ranf)
      R_EXTRA=2.D0*RANF
      IF(R_EXTRA.LT.1.d0) THEN
         IPDECIDE=1
      ELSE
         IPDECIDE=2
      END IF
 
      IF(IPDECIDE.EQ.0) IPDECIDE=1

! IPDECIDE=1 --> CHANGE JBETA !
! IPDECIDE=2 --> CHANGE IALPHA !

      IALPHA_N=IALPHA
      JBETA_N=JBETA
      IF(IPDECIDE.LE.1) THEN    !jbeta change
        IF(JBETA.EQ.1) THEN
           JBETA_N=2
           MU=1
           NU=2
           YES_J12=dcoup(1,2)
           vdotd=YES_J12
!          write(333,*) dcoup(1,2)
           YES_J12=DTSLICE*ABS(YES_J12)
           PYES_J12=YES_J12/(1.D0+YES_J12)
           PNO_J12=1.d0-PYES_J12
           PYES_J=PYES_J12
           PNO_J=1.d0-PYES_J
        ELSEIF(JBETA.EQ.2) THEN
          JBETA_N=1
          MU=2
          NU=1
          YES_J21=dcoup(2,1)
          vdotd=YES_J21
          YES_J21=DTSLICE*ABS(YES_J21)
          PYES_J21=YES_J21/(1.D0+YES_J21)
          PNO_J21=1.d0-PYES_J21
          PYES_J=PYES_J21
          PNO_J=1.d0-PYES_J
        ENDIF
      ELSEIF(IPDECIDE.GT.1) THEN    !ialpha change
        IF(IALPHA.EQ.1) THEN
          IALPHA_N=2
          MU=1
          NU=2
          YES_J12=dcoup(1,2)
          vdotd=YES_J12
          YES_J12=DTSLICE*ABS(YES_J12)
          PYES_J12=YES_J12/(1.d0+YES_J12)
          PNO_J12=1.d0-PYES_J12
          PYES_J=PYES_J12
          PNO_J=1.d0-PYES_J
        ELSEIF(IALPHA.EQ.2) THEN
          IALPHA_N=1
          MU=2
          NU=1
          YES_J21=dcoup(2,1)
          vdotd=YES_J21
          YES_J21=DTSLICE*ABS(YES_J21)
          PYES_J21=YES_J21/(1.d0+YES_J21)
          PNO_J21=1.d0-PYES_J21
          PYES_J=PYES_J21
          PNO_J=1.d0-PYES_J
        ENDIF
      ENDIF

      IF((MU.EQ.1).AND.(NU.EQ.2)) THEN  !!!  COMPUTE MOMENTUM SHIFT  !!! 
        if(check.lt.abs(e_r(1)-e_r(2))) then
          nrej1=nrej1+1
          LREJECT1=.TRUE.
          pyes_j=0.d0
          pno_j=1.d0-pyes_j
        endif
      ENDIF

!  !!!!!!!!!            TESTING JUMP           !!!!!!!!!!
      CALL RANDOM_NUMBER(r_extra)
!     write(444,*) pyes_j,r_extra
      IF((PYES_J.GT.R_EXTRA)) THEN !!! MAKE JUMP !!!
         IF((MU.EQ.1).AND.(NU.EQ.2)) THEN  !!!  COMPUTE MOMENTUM SHIFT  !!! 
           P_SHIFT=dcoup(1,2)
           DO IS=1,NAT
            S_SHIFT(IS)=S(1,2,IS)
!            S_SHIFTY(IS)=Sy(1,2,IS)
!            S_SHIFTZ(IS)=Sz(1,2,IS)
            V(IS)=V(IS)+S_SHIFT(IS)
!            VY(IS)=VY(IS)+S_SHIFTY(IS)
!            VZ(IS)=VZ(IS)+S_SHIFTZ(IS)
           END DO
           PI_WEIGHT=PI_WEIGHT*1.D0/PYES_J
         ELSEIF((MU.EQ.2).AND.(NU.EQ.1)) THEN
           P_SHIFT=dcoup(2,1)
           DO IS=1,NAT
            S_SHIFT(IS)=S(2,1,IS)
!            S_SHIFTY(IS)=Sy(2,1,IS)
!            S_SHIFTZ(IS)=Sz(2,1,IS)
            V(IS)=V(IS)+S_SHIFT(IS)
!            VY(IS)=VY(IS)+S_SHIFTY(IS)
!            VZ(IS)=VZ(IS)+S_SHIFTZ(IS)
           END DO
           PI_WEIGHT=PI_WEIGHT*1.D0/PYES_J
         ENDIF

!        CALL RATTLE(NSCONS,TOLNCE,LISTCON,BOXL,DXT,DYT,DZT,
!    &               MASS,DT,DDD,X,Y,Z,VX,VY,VZ)

!          !!! MOMENTUM SHIFT APPLIED  !!!
         P_FACT= P_FACT * 2.D0 * DTSLICE * P_SHIFT
         WEIGHT_SLICE=  2.D0 * DTSLICE * 1.D0/PYES_J
         IALPHA=IALPHA_N    !!! changing state labels  !!!
         JBETA=JBETA_N      !!!     after jump         !!!
         NJUMP=NJUMP+1
      ELSE
           PI_WEIGHT=PI_WEIGHT*1.D0/PNO_J
           WEIGHT_SLICE=1.D0/PNO_J
      ENDIF

      RETURN
      END 


      SUBROUTINE OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,OBS,nmds,Nb, &
                 IA,JB,POP,IALPHA,JBETA,RHOW,IB,qmat,rhonorm, &
                 NBCUT)  

! This general subroutine has been modified for the calculation of the 
! 'diagonal' part of the rate constant.

      IMPLICIT NONE
      INCLUDE 'param.h'
      REAL*8 FACTRE,FACTIM
      REAL*8 POP(NMDSBX,n,n)
      REAL*8 obs(n,n),qmat(n,n)
      REAL*8 OBSERV_RE,Nb(nat)
      REAL*8 OBSERV_IM
      REAL*8 PI_WEIGHT
      REAL*8 P_FACT
      REAL*8 RHOW(n,n),RHONORM(NMDSBX,n,n)
      REAL*8 rhob
      INTEGER ialpha,jbeta,ia,jb,nmds
      INTEGER i,j,ib
      COMPLEX*16 WPHI_FACTOR
      LOGICAL NBCUT

! Calculate observable in adiabatic basis

!     do i=1,nlv
!       do j=1,nlv
!         OBS(i,j)=qmat(i,j)
          OBS(ialpha,jbeta)=Nb(ib)
!       enddo
!     enddo

      FACTRE=P_FACT*DBLE(WPHI_FACTOR)*PI_WEIGHT
      FACTIM=P_FACT*IMAG(WPHI_FACTOR)*PI_WEIGHT

      OBSERV_RE=FACTRE*OBS(IALPHA,JBETA)
      OBSERV_IM=FACTIM*OBS(IALPHA,JBETA)

      Nb(ib)=observ_re

! Place cut-off on observable

      if(Nb(ib).ge.20) then
          Nb(ib)=20.d0
      endif
      if(Nb(ib).le.-20) then
          Nb(ib)=-20.d0
      endif

! Calculate elements of sum which contribute to Pop=Tr(rhow*observ)

!     IF((IA.EQ.1).AND.(JB.EQ.1)) THEN
!       POP(IB,1,1)=POP(IB,1,1)+OBSERV_RE
!       RHONORM(IB,1,1)=RHONORM(IB,1,1)+RHOW(2,2)
!     ELSEIF((IA.EQ.2).AND.(JB.EQ.2)) THEN
!       POP(IB,2,2)=POP(IB,2,2)+RHOW(2,2)*OBSERV_RE
!       RHONORM(IB,2,2)=RHONORM(IB,2,2)+RHOW(2,2)
!     ELSEIF((IA.EQ.1).AND.(JB.EQ.2)) THEN
!       POP(IB,1,2)=POP(IB,1,2)+RHOW(2,1)*OBSERV_RE
!       RHONORM(IB,1,2)=RHONORM(IB,1,2)+RHOW(2,2)
!     ELSEIF((IA.EQ.2).AND.(JB.EQ.1)) THEN
!       POP(IB,2,1)=POP(IB,2,1)+RHOW(1,2)*OBSERV_RE
!       RHONORM(IB,2,1)=RHONORM(IB,2,1)+RHOW(2,2)
!     ENDIF

      RETURN
      END 




