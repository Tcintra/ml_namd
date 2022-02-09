      PROGRAM REORDER

      IMPLICIT NONE
      
      Integer i,j,n,l,nb
      parameter(nb=25)
      REAL*8  Vpes2,V2(nb**2,nb**2)
      
! Informative files
      OPEN(unit=4,file="fort.100")
      OPEN(unit=444,file="../fort.444")
      
! Output file
      DO i=1,nb**2
       j=0
         Do l=1,nb**2
          READ(4,*) n,n,n,n,Vpes2
          j=j+1
          V2(i,j)=Vpes2
          WRITE(444,*) i,j,V2(i,j)
        ENDDO
       ENDDO

 END
