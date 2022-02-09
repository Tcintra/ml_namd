      PROGRAM PCET_SIMPLE
!=======================================================
! SIMPLE PCET MODEL - SLICE ALGORITHM                
!=======================================================

      IMPLICIT NONE

      include 'param.h'
      real*8 ms,ws,s1,s2,coup,mb,wc,Wb(mode),c(mode)
      real*8 Qb(mode),Vb(mode),qs(nat),Vs(nat),Fs(nat),Fb(mode),pcx
      real*8 gaussn,sp_density,avg_sp_density(nlv,nlv),fact(nlv)
      real*8 Hmat(n,n),vect(n,n)
      real*8 HF(nat,nlv),d(nlv,nlv,nat),sum,S(nlv,nlv,nat)
      real*8 pi_weight,sam,sample
      real*8 pi,alpha,alpha1,alpha2
      REAL*8 p_fact,weight_slice
      real*8 wphi(nlv,nlv),rhow(nlv,nlv)
      real*8 dcoup(nlv,nlv),qmat(nlv,nlv),w(n)
      real*8 fcs,obs(nlv,nlv)
      real*8 STAT_WEIGHT(NMDSBX,nlv,nlv),OBSERV_IM(nmdsbx,nlv,nlv)
      REAL*8 GNORM(NMDSBX,nlv,nlv),RHONORM(NMDSBX,nlv,nlv),RATE(NMDSBX,nlv,nlv)
      REAL*8 RATEa(NMDSBX,nlv,nlv),RATEb(NMDSBX,nlv,nlv),RATEc(NMDSBX,nlv,nlv), &
             RATEd(NMDSBX,nlv,nlv)
      REAL*8 RATEe(NMDSBX,nlv,nlv),RATEf(NMDSBX,nlv,nlv),RATEg(NMDSBX,nlv,nlv), &
             RATEh(NMDSBX,nlv,nlv)
      REAL*8 krp,krp_a,krp_b,krp_c,krp_d,krp_e,krp_f,krp_g,krp_h
      REAL*8 PHASE_FACT_RE(NMDSBX,nlv,nlv),PHASE_FACT_IM(NMDSBX,nlv,nlv)
      REAL*8 BETA,DT,NA,NAeq
      REAL*8 p(n,n),p_e(n,n),p_s(n,n),d_c(n,n),olA(0.5*n,0.5*n),olB(0.5*n,0.5*n)
      parameter(pi=3.14159d0)
      INTEGER NATOM,modes,NMDS,NMCS,MCS,NJUMP,NJUMP_CUT,NJUMP_TOT,IDUM,time,seed_dimension
      INTEGER, Allocatable :: seed(:)
      INTEGER IALPHA,JBETA,RRELX,neqs
      INTEGER k,I,J,ni,nj,ITBLOCK,IT,nstlim,nrej1,nrej2,nrej3(nmdsbx)
      INTEGER IA,JB,IB,NJJ,NSLICE
      INTEGER NJUMPM(NMDSBX,nlv,nlv)
      INTEGER NJP_STAT(NMDSBX,NJC_MAX)
      INTEGER INCO,INIT,IPRINT
      COMPLEX*16 stat_slice(NMDSBX,nlv,nlv),wphi_factor
      logical lreject1,lreject2
      COMMON/RANDOM/IDUM 
      include 'sol_bath.in'
      

      READ (5,*) 
      READ (5,*) NATOM,modes,DT
      READ (5,*) 
      READ (5,*) BETA,INCO,INIT
      READ (5,*) 
      READ (5,*) NMDS,NSTLIM,IPRINT,neqs
      READ (5,*) 
      READ (5,*) IDUM
      READ (5,*) 
      READ (5,*) NMCS,ITBLOCK,NJUMP_CUT,NSLICE

! Input files
       OPEN(unit=10,file="matrix.in")
       OPEN(unit=20,file="overlapA.in")
       OPEN(unit=30,file="overlapB.in")

! Read in field file      
 
      DO I=1,n
        DO J=1,n
          READ(10,*) ni,nj,p(i,j),p_e(i,j),p_s(i,j),d_c(i,j)
        ENDDO
      ENDDO
     CLOSE(10)

      DO I=1,0.5*n
        DO J=1,0.5*n
          READ(20,*) ni,nj,olA(i,j)
          READ(30,*) ni,nj,olB(i,j)
        ENDDO
      ENDDO
     CLOSE(20)
     close(30)
! Initialize RNG seed
     CALL RANDOM_SEED(size=seed_dimension)
     ALLOCATE (seed(seed_dimension))
     do i=1,seed_dimension
!        seed(i) = IDUM+3*i-1
         seed(i) = time()+3*i-1
     enddo
     CALL RANDOM_SEED(PUT=seed)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                             !!!
!!!             P_shift = momentum shift factor                 !!!
!!!            GNORM = normalization                            !!! 
!!!                                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO JB=1,nlv
       DO IA=1,nlv
        DO IB=1,NMDSBX
          GNORM(IB,IA,JB)=0.d0
          RHONORM(IB,IA,JB)=0.d0
        END DO
       END DO
      END DO
      
      NJUMP_TOT=0
      nrej1=0
      nrej2=0
      DO IB=1,NMDSBX
        NREJ3(IB)=0
      END DO
                         
      DO JB=1,nlv
       DO IA=1,nlv
        DO IB=1,NMDSBX
          RATE(IB,IA,JB)=0.D0
          RATEa(IB,IA,JB)=0.D0
          RATEb(IB,IA,JB)=0.D0
          RATEc(IB,IA,JB)=0.D0
          RATEd(IB,IA,JB)=0.D0
          RATEe(IB,IA,JB)=0.D0
          RATEf(IB,IA,JB)=0.D0
          RATEg(IB,IA,JB)=0.D0
          RATEh(IB,IA,JB)=0.D0
        END DO
       END DO
      END DO

      DO IB=1,NMDSBX 
       DO JB=1,nlv
        DO IA=1,nlv
          NJUMPM(IB,IA,JB)=0
          STAT_WEIGHT(IB,IA,JB)=0.d0
          STAT_SLICE(IB,IA,JB)=CMPLX(0.d0,0.d0)
          PHASE_FACT_RE(IB,IA,JB)=0.d0
          PHASE_FACT_IM(IB,IA,JB)=0.d0
        END DO
       END DO
       DO NJJ=1,NJC_MAX
         NJP_STAT(IB,NJJ)=0
         IF((IB.EQ.1).AND.(NJJ.EQ.1)) NJP_STAT(IB,NJJ)=NMCS
       END DO
      END DO
      
      DO IA=1,nlv
        DO JB=1,nlv

!       if (IA.ne.JB) then
        fact(ia)=0.d0
        avg_sp_density(ia,jb)=0.d0

        DO MCS=1,NMCS            !!! SAMPLING LOOP STARTS HERE !!!
                          !!! PHASE SPACE POINT RANDOM SAMPLING !!!

        sp_density=0.d0
    
! Initialize position and momentum of solvent
       
        sum=0.d0
        DO I=1,mode
          Wb(I) = -wc*log((I-0.5)/mode)
          Qb(I)=gaussn()/(Wb(I)*DSQRT(BETA*Mb))
          Vb(I)=gaussn()/DSQRT(BETA*Mb)
!         c(I) = Wb(I)*(2.d0*coup*Mb*wc/(mode*pi))
          c(I) = Wb(I)*dsqrt(2.d0*coup*Mb*wc/(mode*pi))
!         sum=sum+(Wb(I)/c(I))*DSQRT(Mb/BETA)
        ENDDO

        DO k=1,nat
!         qs(k)=gaussn()*sum
          qs(k)=gaussn()/(ws*DSQRT(BETA*ms))-0.75d0
          Vs(k)=gaussn()/DSQRT(BETA*ms)
!         write(222,*) qs(k),vs(k)
        ENDDO

! Remove COM motion
!        pcx=0.d0
!        do i=1,mode
!          pcx=pcx+mb*vb(i)
!        enddo
!        pcx=pcx+ms*vs(1)
!        pcx=pcx/(nat+mode)
!        do i=1,mode
!          vb(i)=vb(i)-pcx/mb
!        enddo
!        vs(1)=vs(1)-pcx/ms 

        CALL FORCE(qs,c,Wb,Qb,Fb,Fs)
        CALL ADIAB_NVE_I(p,p_e,p_s,d_c,w,Hmat,vect,HF,d,qs,c,Wb,Qb,Vs,Vb)
          
        CALL RUNNVT_NHC(neqs,HF,qs,Vs,Fs,Qb,Vb,Fb,Wb,c, &
                 p,p_e,p_s,d_c,w,Hmat,vect,d,IA,JB,beta)

!         write(666,*) vs(1)
!         do k=1,mode
!           write(777,*) vb(k)
!         enddo
 
        CALL FORCE(qs,c,Wb,Qb,Fb,Fs)
        CALL ADIAB_NVE(p,p_e,p_s,d_c,w,Hmat,vect,HF,d,qs,c,Wb,Qb,Vs,Vb)

!======================================================================
!Calculate unbiasing factors for sampling
        if (ia.eq.1) then
           fact(ia)=fact(ia)+1.d0+exp(-beta*(w(2)-w(1)))
        elseif (ia.eq.2) then
           fact(ia)=fact(ia)+1.d0+exp(-beta*(w(1)-w(2)))
        endif
!=============================================================== 

        NJUMP=0  !!! === jump counter for the single trajectory === !!!
        IALPHA=IA
        JBETA=JB
        
!======================================================================
!Calculate the spectral density related to reactant projection
!operator at time zero

            do i=1,0.5*n
              do j=1,0.5*n
               sp_density=sp_density+vect(i,ialpha)*vect(j,jbeta)*olA(i,j)
              enddo
            enddo
            avg_sp_density(ia,jb)=avg_sp_density(ia,jb)+sp_density 
!=============================================================== 

      CALL RUNNVE(DT,nmds,IPRINT,ms,w,olA,olB,sp_density, &
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

       NJUMP_TOT=NJUMP_TOT+NJUMP
     
      ENDDO    !!! MONTE CARLO LOOP ENDS HERE !!!
        ENDDO  
      ENDDO  ! End loop on adiabatic states

      open(28,file='rate.out')
      open(38,file='ratea.out')
      open(48,file='rateb.out')
      open(58,file='ratec.out')
      open(68,file='rated.out')
      open(78,file='ratee.out')
      open(88,file='ratef.out')
      open(98,file='rateg.out')
      open(108,file='rateh.out')

!      NAeq=0.d0
!      do ia=1,nlv
!        do jb=1,nlv
!          NAeq=NAeq+avg_sp_density(ia,jb)/fact(ia)
!        enddo
!      enddo
!      write(222,*) NAeq

      DO IB=1,NMDS/ITBLOCK
        krp=0.d0
        krp_a=0.d0
        krp_b=0.d0
        krp_c=0.d0
        krp_d=0.d0
        krp_e=0.d0
        krp_f=0.d0
        krp_g=0.d0
        krp_h=0.d0
        do ia=1,nlv
          do jb=1,nlv
!            write(555,'(3i7,2g20.10)') ib,ia,jb,rate(ib,ia,jb),fact(ia)
             krp=krp+rate(ib,ia,jb)/fact(ia)
!            write(555,'(3i7,g20.10)') ib,ia,jb,rate(ib,ia,jb)/fact(ia)
             krp_a=krp_a+ratea(ib,ia,jb)/fact(ia)
             krp_b=krp_b+rateb(ib,ia,jb)/fact(ia)
             krp_c=krp_c+ratec(ib,ia,jb)/fact(ia)
             krp_d=krp_d+rated(ib,ia,jb)/fact(ia)
             krp_e=krp_e+ratee(ib,ia,jb)/fact(ia)
             krp_f=krp_f+ratef(ib,ia,jb)/fact(ia)
             krp_g=krp_g+rateg(ib,ia,jb)/fact(ia)
             krp_h=krp_h+rateh(ib,ia,jb)/fact(ia)
          enddo
       enddo
       krp=(2.d0*krp)/(beta)
       krp_a=(2.d0*krp_a)/(beta)
       krp_b=(2.d0*krp_b)/(beta)
       krp_c=(2.d0*krp_c)/(beta)
       krp_d=(2.d0*krp_d)/(beta)
       krp_e=(2.d0*krp_e)/(beta)
       krp_f=(2.d0*krp_f)/(beta)
       krp_g=(2.d0*krp_g)/(beta)
       krp_h=(2.d0*krp_h)/(beta)
       write(28,'(i7,5g20.10)') ib,rate(ib,1,1)/fact(1),rate(ib,1,2)/fact(1), &
                                rate(ib,2,1)/fact(2),rate(ib,2,2)/fact(2),krp
       write(38,'(i7,5g20.10)') ib,ratea(ib,1,1)/fact(1),ratea(ib,1,2)/fact(1), &
                                ratea(ib,2,1)/fact(2),ratea(ib,2,2)/fact(2),krp_a
       write(48,'(i7,5g20.10)') ib,rateb(ib,1,1)/fact(1),rateb(ib,1,2)/fact(1), &
                                rateb(ib,2,1)/fact(2),rateb(ib,2,2)/fact(2),krp_b
       write(58,'(i7,5g20.10)') ib,ratec(ib,1,1)/fact(1),ratec(ib,1,2)/fact(1), &
                                ratec(ib,2,1)/fact(2),ratec(ib,2,2)/fact(2),krp_c
       write(68,'(i7,5g20.10)') ib,rated(ib,1,1)/fact(1),rated(ib,1,2)/fact(1), &
                                rated(ib,2,1)/fact(2),rated(ib,2,2)/fact(2),krp_d
       write(78,'(i7,5g20.10)') ib,ratee(ib,1,1)/fact(1),ratee(ib,1,2)/fact(1), &
                                ratee(ib,2,1)/fact(2),ratee(ib,2,2)/fact(2),krp_e
       write(88,'(i7,5g20.10)') ib,ratef(ib,1,1)/fact(1),ratef(ib,1,2)/fact(1), &
                                ratef(ib,2,1)/fact(2),ratef(ib,2,2)/fact(2),krp_f
       write(98,'(i7,5g20.10)') ib,rateg(ib,1,1)/fact(1),rateg(ib,1,2)/fact(1), &
                                rateg(ib,2,1)/fact(2),rateg(ib,2,2)/fact(2),krp_g
       write(108,'(i7,5g20.10)') ib,rateh(ib,1,1)/fact(1),rateh(ib,1,2)/fact(1), &
                                 rateh(ib,2,1)/fact(2),rateh(ib,2,2)/fact(2),krp_h
      ENDDO


      END 
  
