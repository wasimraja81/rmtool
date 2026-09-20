c        This is a demo Code
c        Date: 13-02-2008
c        Wasim Raja


         implicit none
         integer*4 N
         real*4 A(100), mu
         integer*4 i
         N = 10
         do i = 1,N
              A(i) = i
              write(*,*)A(i)
         end do
         call mean(A,N,mu)
         write(*,*)'mean = ',mu
         end

         include '../mean.f'
