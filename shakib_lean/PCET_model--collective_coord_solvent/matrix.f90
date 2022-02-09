      PROGRAM MATRIX
      IMPLICIT NONE

      Integer i,j,n,l,m
      parameter(n=625)
      REAL*8  S,V,K,Vpes1,Vpes2
     
! Input files      
      
    OPEN(unit=1,file="fort.111")
    OPEN(unit=2,file="fort.222")
    OPEN(unit=3,file="fort.333")
    OPEN(unit=4,file="fort.444")     
    OPEN(unit=5,file="fort.555")     
    OPEN(unit=1000,file="matrix.in")     

! Output file
        
    DO l=1,n
      DO m=1,n
        READ(1,*) i,j,S
        READ(2,*) i,j,V         
        READ(3,*) i,j,Vpes1
        READ(4,*) i,j,Vpes2
        READ(5,*) i,j,K
        WRITE(1000,'(2i10,5g25.15)') l,m,S,V,Vpes1,Vpes2,K 
      ENDDO
    ENDDO
         
   END 

