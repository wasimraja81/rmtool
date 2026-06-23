chelp+
! Driver program for subroutine "SHIFT"
!
!         --wasim raja, 29 Dec, 2009
!----------------------------------------
chelp-

      implicit none
      integer*4 N,nshift
      parameter(N=6)
      integer*4 i
      real*4 A(N), O(N)

      !nshift = int(N/2) ! for fftshift
      do i = 1,N
         A(i) = i
      enddo
      write(*,*)"Before shift: "
      write(*,*)(A(i),i = 1,N)
      write(*,*)" "
      write(*,*)'No. of positions to shift: '
      read(*,*)nshift
      call shift(A,N,nshift)
      write(*,*)"After shift: "
      write(*,*)(A(i),i = 1,N)
      end

      include '../shift.f'

