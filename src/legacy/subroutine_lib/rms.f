c         SUBROUTINE: rms(V,n,x)

c         Author: Wasim Raja
c         Date: 13-02-2008
c         
c         This subroutine calculates the rms "x" about mean 
c         of a vector "V" given its length "N". 

          subroutine rms(V,N,x)
          real*4 V(*), x, mu
          integer*4 n
          integer*4 i
          real*4 accum

          call mean(V,N,mu)
          accum = 0.0
          do i = 1,N
               accum = accum + (V(i) - mu)**2
          end do
          x = sqrt(accum/N)
          end
          include 'mean.f'
