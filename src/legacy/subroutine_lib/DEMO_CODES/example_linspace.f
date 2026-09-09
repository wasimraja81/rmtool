c        This is a demo Code
c        Date: 07-08-2007
c        Wasim Raja


         implicit none
         real*4 aa, bb
         integer*4 N
         real*4 A(100)
         integer*4 i
         aa = 1.0
         bb = 97.0
         N = 10
         call linspace(aa,bb,N,A)
         do i = 1,N
              write(*,*)A(i)
         end do
c        write(*,*)(A(i),i=1,N)
         end

         include '../linspace.f'
