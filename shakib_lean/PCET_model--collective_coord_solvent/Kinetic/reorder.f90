      PROGRAM REORDER

! This code reads the output of each submatrix which
! are based on the indices of proton and electron harmonic
! oscillators and wirtes it in a new file with matrix
! indices from 1 to 625

      IMPLICIT NONE
      
      Integer i,j,n,l,nb
      parameter(nb=25)
      REAL*8  Kinetic,K(nb**2,nb**2)
      
! Informative files
      OPEN(unit=4,file="fort.100")
      OPEN(unit=555,file="../fort.555")
      
! Output file
      DO i=1,nb**2
       j=0
         Do l=1,nb**2
          READ(4,*) n,n,n,n,Kinetic
          j=j+1
          K(i,j)=Kinetic
          WRITE(555,*) i,j,K(i,j)
        ENDDO
       ENDDO

 END
