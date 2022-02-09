      FUNCTION fact(n) 
      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: n 
      INTEGER :: f
      real*8 fact 
      fact = 1
      DO f = 1, n
      fact = fact * f
      ENDDO
      RETURN
      END FUNCTION fact

      FUNCTION pl(n,x)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER, INTENT(IN)                      :: n
      DOUBLE PRECISION, INTENT(IN)             :: x
                
      a=2.0D0
      b=0.0D0
      c=1.0D0
      y0=1.0D0
      y1=2.0D0*x
      if (n==0) then
      pl=1.0D0
      else if (n==1) then
      pl=2.0D0*x
      else if (n>=2) then
      do k=2,n
        c=2.0D0*(k-1.0D0)
        yn=(a*x+b)*y1-c*y0
        pl=yn
        y0=y1
        y1=yn
        ENDDO
      ENDIF
      RETURN
      END FUNCTION pl

    
