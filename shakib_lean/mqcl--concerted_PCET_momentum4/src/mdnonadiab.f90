      PROGRAM PCET_SIMPLE
!=======================================================
! SIMPLE PCET MODEL - SLICE ALGORITHM                
!=======================================================

      IMPLICIT NONE

      include 'param.h'
      integer ni,nj
      real*8 sigv,v(nat),mass(nat),gaussn
      real*8 ws,R(nat),R0
      real*8 adiab,adiabm,adiabn
      real*8 Hmat(n,n),Smat(n,n),vect(n,n)
      real*8 vs,vlow(n,n),elow1,elow2,elow3
      real*8 HF(nlv),d_tot(nlv,nlv),sum
      real*8 pi_weight
      real*8 hbar,pi,alpha,alpha1,alpha2
      REAL*8 p_fact,S_SHIFT_CUT,weight_slice
      real*8 wphi(nlv,nlv),rhow(nlv,nlv)
      real*8 dcoup(nlv,nlv),qmat(nlv,nlv)
      real*8 fcs,obs(nlv,nlv)
      real*8 nb(nat),flux(nat),kappa,D
      real*8 de_bar,dedt,dedx(nat),unbias
      real*8 STAT_WEIGHT(NMDSBX,nlv,nlv)
      REAL*8 GNORM(NMDSBX,nlv,nlv),RHONORM(NMDSBX,nlv,nlv),POP(NMDSBX,nlv,nlv)
      REAL*8 PHASE_FACT_RE(NMDSBX,nlv,nlv),PHASE_FACT_IM(NMDSBX,nlv,nlv)
      REAL*8 BETA,DT,BOXL
      REAL*8 S(n,n),Pot(n,n),Vpes1(n,n),Vpes2(n,n),Kin(n,n)
      real*8 tolnce,ddd,f(nat),charge
      INTEGER NATOM,NMDS,NMCS,MCS,NJUMP,NJUMP_CUT,NJUMP_TOT,IDUM,seed_dimension
      INTEGER, Allocatable :: seed(:)
      INTEGER IALPHA,JBETA,RRELX
      INTEGER I,J,ITBLOCK,IT,nstlim,nrej1,nrej2,nrej3(nmdsbx)
      INTEGER IA,JB,IB,NJJ,NSLICE
      INTEGER NJUMPM(NMDSBX,nlv,nlv)
      INTEGER NJP_STAT(NMDSBX,NJC_MAX)
      INTEGER INCO,INIT,IPRINT
      COMPLEX*16 stat_slice(NMDSBX,nlv,nlv),wphi_factor
      logical lreject1,lreject2
      COMMON/RANDOM/IDUM 
      include 'solvent.in'
      data pi/3.141592653d0/

      READ (5,*) 
      READ (5,*) NATOM,DT,BOXL
      READ (5,*) 
      READ (5,*) BETA,INCO,INIT
      READ (5,*) 
      READ (5,*) NMDS,NSTLIM,IPRINT
      READ (5,*) 
      READ (5,*) S_SHIFT_CUT
      READ (5,*) 
      READ (5,*) NMCS,ITBLOCK,NJUMP_CUT,NSLICE

! Input files
       OPEN(unit=10,file="matrix.in")

! Output files
      OPEN(60,FILE='info.out',FORM='formatted',STATUS='unknown')
      OPEN(70,FILE='njump.out',FORM='formatted',STATUS='unknown')

! Read in field file      
 
      DO I=1,n
        DO J=1,n
          READ(10,*) ni,nj,S(i,j),Pot(i,j),Vpes1(i,j),Vpes2(i,j),Kin(i,j)
        ENDDO
      ENDDO
     CLOSE(10)

! Initialize RNG seed
     CALL RANDOM_SEED(size=seed_dimension)
     ALLOCATE (seed(seed_dimension))
     do i=1,seed_dimension
       seed(i) = 3*2**i-1
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
          POP(IB,IA,JB)=0.D0
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

      DO MCS=1,NMCS            !!! SAMPLING LOOP STARTS HERE !!!
                          !!! PHASE SPACE POINT RANDOM SAMPLING !!!

      ! Initialize position and momentum of solvent
        R(1)=-0.45d0
        V(1)=0.000182d0
!       DO I=1,nat
!         SIGV=DSQRT(BETA/MASS(I))
!         V(I)=SIGV*(GAUSSN()) 
!       ENDDO 
  
        CALL FORCE(R,VS,F)
                  
        CALL ADIAB_NVE(Kin,Pot,Vpes1,Vpes2,S,adiab,adiabm,adiabn, &
             Hmat,Smat,vect,vlow,elow1,elow2,elow3,HF,d_tot,r)

      JB=1
      IA=1
        NJUMP=0  !!! === jump counter for the single trajectory === !!!
        IALPHA=IA
        JBETA=JB
      
      CALL RUNNVE(DT,nmds,IPRINT,TOLNCE,MASS, &
                 BOXL,R,V,F,PI,CHARGE,HF,ELOW1,IT,SMAT,VS, &
                 ALPHA,BETA,elow2,elow3,mcs,vpes1,vpes2,Kin,Pot, &
                 de_bar,nb,flux,unbias,d_tot,dedt,wphi, &
                 dcoup,S,hbar,IALPHA,JBETA,WPHI_FACTOR, &
                 IDUM,P_FACT,NJUMP,S_SHIFT_CUT,LREJECT1, &
                 PI_WEIGHT,IA,JB,GNORM,OBS,POP,RHOW,ITBLOCK, &
                 NJUMP_CUT,NSLICE,NJUMPM,NJP_STAT,STAT_WEIGHT, &
                 PHASE_FACT_RE,PHASE_FACT_IM,STAT_SLICE, &
                 WEIGHT_SLICE,RHONORM,qmat,lreject2,kappa, &
                 nrej1,nrej2,nrej3,D)

       NJUMP_TOT=NJUMP_TOT+NJUMP
       write(70,*) MCS,NJUMP    
      ENDDO    !!! MONTE CARLO LOOP ENDS HERE !!!
      
      END 
  
