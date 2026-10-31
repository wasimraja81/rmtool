c         Author: Wasim Raja
c         Date: 14-11-2007
c         Last modified on: 14-11-2007

c         dotproduct_c(A,B,P,m) is a subroutine that 
c         multiplies two complex vectors A and B of dimensions 
c         m, and calculates the dot product P of A and B.
c         
c
          subroutine dotproduct_c(A,B,P,m)
          complex A(*), B(*)
          real*4 P
          integer*4 m
          integer*4 i
          real*4 accum

          accum = 0.0
          do i = 1,m
               accum = accum + A(i)*B(i)
          end do
          P = accum
          end

