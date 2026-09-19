c         Author: Wasim Raja
c         Date: 07-08-2007
c         Last modified on: 07-08-2007

c         linspace(BASE,LIMIT,N,V) is a subroutine that 
c         generates a vector V of N linearly spaced elements 
c         between BASE and LIMIT. N has to be greater 
c         than 1. Base and LIMIT are always included in the 
c         range. If LIMIT < BASE, then the elements are 
c         stored in descending order. If N is not a natural,
c         number, a default value of N = 100 is used!
c
c         THE LINSPACE SUBROUTINE DOES NOT MODIFY ANY INPUT 
c         ARGUMENT!
          subroutine linspace(BASE,LIMIT,N,V)
          real*4 BASE, LIMIT, V(*)
          integer*4 N
          integer*4 i
          real*4 h
c         if N has been specified correctly, then N = N 
c         else N = 100
          if (N.lt.1)then
               write(*,*)'----------------- WARNING -------------------'
               write(*,*)' Wrong vector length, N will be forced to 100'
               write(*,*)'---------------------------------------------'
               N = 100
          end if
          h = (LIMIT - BASE)/(N-1)
          do i = 1,N
               V(i) = BASE + (i-1)*h
          end do
          end
