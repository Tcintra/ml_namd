      PROGRAM REORDER

      IMPLICIT NONE
      
      Integer i,j,n,l,nb
      parameter(nb=25)
      REAL*8  Potential,P(nb**2,nb**2)
      
! Informative files
      OPEN(unit=4,file="fort.100")
      OPEN(unit=222,file="../fort.222")
      
! Output file
      DO i=1,nb**2
       j=0
         Do l=1,nb**2
          READ(4,*) n,n,n,n,Potential
          j=j+1
          P(i,j)=Potential
          WRITE(222,*) i,j,P(i,j)
        ENDDO
       ENDDO

 END
