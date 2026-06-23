c        This is a demo Code
c        Date: 07-08-2007
c        Wasim Raja


         implicit none
         real*4 A(10), B(10), P
         integer*4 m
         integer*4 i
         m = 10

         call linspace(1.0,10.0,10,A)
         call linspace(0.0,10.0,10,B)
         write(*,*)'A     B'
         do i = 1,m
              write(*,*)A(i), B(i)
         end do
         call dotproduct(A,B,P,m)
         write(*,*)'A.B = ', P
c        write(*,*)(A(i),i=1,N)
         end

         include '../dotproduct.f'
         include '../linspace.f'
