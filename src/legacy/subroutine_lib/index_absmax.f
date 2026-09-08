c         Author: Wasim Raja
c         Date: 12-05-2009
c         
c         This subroutine locates the "index" of maxima of 
c         the absolute values of the elements of a vector 
c         "V" given its length "N". It also returns the avg 
c         of the absolute values. 

          subroutine index_absmax(V,N,imax,avg_abs)

          implicit none
          real*4    V(*), avg_abs, maxx, absV 
          integer*4 N
          integer*4 i, imax

          maxx = abs(V(1))
          imax = 1
          avg_abs = 0.0d0
          do i = 1,N
               absV = abs(V(i))
               avg_abs = avg_abs + absV
               if(absV.ge.maxx)then
                       maxx = absV
                       imax = i 
               endif
          end do
          avg_abs = avg_abs/dble(N)

          return
          end 

