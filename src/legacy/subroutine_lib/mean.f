c   4)   SUBROUTINE: mean(V,n,x)

c         Author: Wasim Raja
c         Date: 13-02-2008
c         
c         This subroutine calculates the mean "x" of a vector "V" 
c         given its length "N". 

          subroutine mean(V,N,x)
          real*4 V(*), x
          integer*4 n
          integer*4 i
          real*4 accum

          accum = 0.0
          do i = 1,N
               accum = accum + V(i)
          end do
          x = accum/N
          end

