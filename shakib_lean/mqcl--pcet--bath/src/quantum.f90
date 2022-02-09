      SUBROUTINE MAKEJUMP(IDUM,IALPHA,JBETA,MASS,Vs,S,P_FACT, &
                          NJUMP,CHECK,w,LREJECT1, &
                          LREJECT2,DT,PI_WEIGHT,NSLICE, &
                          WEIGHT_SLICE,DCOUP, &
                          R,nrej1,nrej2,it)

      IMPLICIT NONE
      INCLUDE 'param.h'
      INTEGER IDUM,IALPHA,JBETA,IALPHA_N,JBETA_N,nrej1,nrej2
      INTEGER MU,NU,IS,I,IPDECIDE,NJUMP,NSLICE
      INTEGER IT
      REAL*8 intlength,currentlength
      REAL*8 MASS(nat),R_EXTRA,RANF
      REAL*8 P_FACT,WEIGHT_SLICE,PI_WEIGHT
      REAL*8 PYES_J,YES_J
      REAL*8 PNO_J,DTSLICE
      REAL*8 tolnce,boxl
      REAL*8 S_SHIFT_CUT,DT
      REAL*8 DCOUP(nlv,nlv),vdotd,w(n),check(nlv,nlv)
      REAL*8 Vs(nat),c3,R(nat)
      REAL*8 S(nlv,nlv,nat),P_SHIFT
      LOGICAL LREJECT1,LREJECT2

      DTSLICE=DBLE(NSLICE)*DT
      IPDECIDE=3
!     LREJECT1=.FALSE.
!     LREJECT2=.FALSE.
 
      do mu=1,nlv
        do nu=1,nlv
          if(mu.lt.nu) then
          if(check(mu,nu).lt.abs(w(mu)-w(nu))) then
            nrej1=nrej1+1
!           LREJECT1=.TRUE.
            dcoup(mu,nu)=0.d0
          endif
          endif
        enddo
      enddo

     intlength=1.d0
      do i=1,nlv
        intlength=intlength+dt*abs(dcoup(jbeta,i))+dt*abs(dcoup(ialpha,i))
      enddo

      CALL RANDOM_NUMBER(ranf)
      R_EXTRA=intlength*RANF
      IF(R_EXTRA.LE.1.d0) THEN
         IPDECIDE=0
         IALPHA_N=IALPHA
         JBETA_N=JBETA
      ENDIF
 
! IPDECIDE=1 --> CHANGE JBETA !
! IPDECIDE=2 --> CHANGE IALPHA !
     
        currentlength=1.d0
        i=1
        do while ((IPDECIDE.EQ.3).and.(i.le.nlv))
          currentlength=currentlength+dt*abs(dcoup(jbeta,i))
          if (R_EXTRA.le.currentlength) then
            IPDECIDE=1
            IALPHA_N=IALPHA
            JBETA_N=i
          endif
          i=i+1
        enddo 
       
        i=1
        do while ((IPDECIDE.EQ.3).and.(i.le.nlv))
          currentlength=currentlength+dt*abs(dcoup(ialpha,i))
          if (R_EXTRA.le.currentlength) then
            IPDECIDE=2
            IALPHA_N=i
            JBETA_N=JBETA
          endif
          i=i+1
        enddo 


      IF(IPDECIDE.EQ.1) THEN    !jbeta change
          PYES_J=dt*abs(dcoup(jbeta,jbeta_n))/intlength
          P_SHIFT=dcoup(jbeta,jbeta_n)
          DO IS=1,NAT
            Vs(is)=Vs(is)+S(jbeta,jbeta_n,is)
          END DO
          PI_WEIGHT=PI_WEIGHT*1.D0/PYES_J
!         WEIGHT_SLICE=  2.D0 * DTSLICE * 1.D0/PYES_J
          P_FACT= P_FACT * DTSLICE * P_SHIFT
          JBETA=JBETA_N      !!!     after jump         !!!
          NJUMP=NJUMP+1
      ELSEIF(IPDECIDE.eq.2) THEN    !ialpha change
          PYES_J=dt*abs(dcoup(ialpha,ialpha_n))/intlength
          P_SHIFT=dcoup(ialpha,ialpha_n)
          DO IS=1,NAT
            Vs(is)=Vs(is)+S(ialpha,ialpha_n,is)
          END DO
          PI_WEIGHT=PI_WEIGHT*1.D0/PYES_J
!         WEIGHT_SLICE=  2.D0 * DTSLICE * 1.D0/PYES_J
          P_FACT= P_FACT * DTSLICE * P_SHIFT
          IALPHA=IALPHA_N    !!! changing state labels  !!!
          NJUMP=NJUMP+1
      ELSEIF(IPDECIDE.eq.0) THEN    !no change
          PNO_J=1.d0/intlength
          PI_WEIGHT=PI_WEIGHT*1.D0/PNO_J
!         WEIGHT_SLICE=1.D0/PNO_J
      ENDIF

      RETURN
      END 


      SUBROUTINE OBSERV(WPHI_FACTOR,PI_WEIGHT,P_FACT,rate,nmds, &
                 IA,JB,sp_density,olB,IALPHA,JBETA,IB,vect, &
                 ratea,rateb,ratec,rated,ratee,ratef,rateg,rateh)  
! This subroutine calculates the instantaneous value of an observable.

      IMPLICIT NONE
      INCLUDE 'param.h'
      REAL*8 FACTRE,FACTIM,FACTIMa,FACTIMb,FACTIMc,FACTIMd
      REAL*8 FACTIMe,FACTIMf,FACTIMg,FACTIMh
      REAL*8 rate(nmdsbx,nlv,nlv),vect(n,n),olB(0.5*n,0.5*n)
      REAL*8 ratea(nmdsbx,nlv,nlv),rateb(nmdsbx,nlv,nlv)
      REAL*8 ratec(nmdsbx,nlv,nlv),rated(nmdsbx,nlv,nlv)
      REAL*8 ratee(nmdsbx,nlv,nlv),ratef(nmdsbx,nlv,nlv)
      REAL*8 rateg(nmdsbx,nlv,nlv),rateh(nmdsbx,nlv,nlv)
      REAL*8 OBSERV_RE,sp_density
      REAL*8 OBSERV_IM
      REAL*8 PI_WEIGHT
      REAL*8 P_FACT
      INTEGER ialpha,jbeta,ia,jb,nmds
      INTEGER i,j,ib
      COMPLEX*16 WPHI_FACTOR

! Calculate observable in adiabatic basis
  
      FACTRE=P_FACT*DBLE(WPHI_FACTOR)*PI_WEIGHT
      FACTIM=P_FACT*IMAG(WPHI_FACTOR)*PI_WEIGHT
!     write(444,'(3i7,3g20.10)') ib,ialpha,jbeta,p_fact,imag(wphi_factor),pi_weight
      
! Place cut-off on observable

      if(FACTIM.ge.5) then
          FACTIMA=5.d0
      elseif(FACTIM.le.-5) then
          FACTIMA=-5.d0
      else
      FACTIMa=FACTIM
      endif

      if(FACTIM.ge.10) then
          FACTIMb=10.d0
      elseif(FACTIM.le.-10) then
          FACTIMb=-10.d0
      else
      FACTIMb=FACTIM
      endif

      if(FACTIM.ge.20) then
          FACTIMc=20.d0
      elseif(FACTIM.le.-20) then
          FACTIMc=-20.d0
      else
      FACTIMc=FACTIM
      endif

      if(FACTIM.ge.50) then
          FACTIMd=50.d0
      elseif(FACTIM.le.-50) then
          FACTIMd=-50.d0
      else
      FACTIMd=FACTIM
      endif

      if(FACTIM.ge.80) then
          FACTIMe=80.d0
      elseif(FACTIM.le.-80) then
          FACTIMe=-80.d0
      else
      FACTIMe=FACTIM
      endif

      if(FACTIM.ge.100) then
          FACTIMf=100.d0
      elseif(FACTIM.le.-100) then
          FACTIMf=-100.d0
      else
      FACTIMf=FACTIM
      endif

      if(FACTIM.ge.150) then
          FACTIMg=150.d0
      elseif(FACTIM.le.-150) then
          FACTIMg=-150.d0
      else
      FACTIMg=FACTIM
      endif

      if(FACTIM.ge.200) then
          FACTIMh=200.d0
      elseif(FACTIM.le.-200) then
          FACTIMh=-200.d0
      else
      FACTIMh=FACTIM
      endif

      do i=0.5*n+1,n
       do j=0.5*n+1,n
!        write(444,'(5i7,5g20.10)') ib,ialpha,jbeta,i,j,vect(i,ialpha),vect(j,jbeta),olB(i-0.5*n,j-0.5*n),sp_density,factim
         rate(ib,ia,jb)=rate(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factim
         ratea(ib,ia,jb)=ratea(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factima
         rateb(ib,ia,jb)=rateb(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimb
         ratec(ib,ia,jb)=ratec(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimc
         rated(ib,ia,jb)=rated(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimd
         ratee(ib,ia,jb)=ratee(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factime
         ratef(ib,ia,jb)=ratef(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimf
         rateg(ib,ia,jb)=rateg(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimg
         rateh(ib,ia,jb)=rateh(ib,ia,jb)+vect(i,ialpha)*vect(j,jbeta)*olB(i-0.5*n,j-0.5*n)*sp_density*factimh
       enddo
      enddo
!      write(555,'(3i7,g20.10)') ib,ia,jb,rate(ib,ia,jb)

      RETURN
      END 




