      PROGRAM REORDER

      IMPLICIT NONE
      
      Integer i,j,n,l,nb
      parameter(nb=25)
      REAL*8  Overlap,S(nb**2,nb**2)
      
! Informative files
      OPEN(unit=4,file="fort.100")
      OPEN(unit=111,file="../fort.111")
      
! Output file
      DO i=1,nb**2
       j=0
         Do l=1,nb**2
          READ(4,*) n,n,n,n,Overlap
          j=j+1
          S(i,j)=Overlap
          WRITE(111,*) i,j,S(i,j)
        ENDDO
       ENDDO

 END
