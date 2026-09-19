c        This is a demo Code
c        Date: 13-02-2008
c        Wasim Raja


         implicit none
         integer*4 N
         real*4 A(100), sigma
         integer*4 i
         N = 10
         do i = 1,N
              A(i) = i
              write(*,*)A(i)
         end do
         call rms(A,N,sigma)
         write(*,*)'rms = ',sigma
         end

         include '../rms.f'
         !include '../fort_lib.f'
