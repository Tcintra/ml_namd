      PROGRAM REORDER

      IMPLICIT NONE
      
      Integer i,j,n,l,nb
      parameter(nb=25)
      REAL*8  Vpes1,V1(nb**2,nb**2)
      
! Informative files
      OPEN(unit=4,file="fort.100")
      OPEN(unit=333,file="../fort.333")
      
! Output file
      DO i=1,nb**2
       j=0
         Do l=1,nb**2
          READ(4,*) n,n,n,n,Vpes1
          j=j+1
          V1(i,j)=Vpes1
          WRITE(333,*) i,j,V1(i,j)
        ENDDO
       ENDDO

 END
