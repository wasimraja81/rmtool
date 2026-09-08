chelp+
! This subroutine shifts the elements of a vector  
! in cyclic order, by a specified number of integer 
! index: "nshift"
! A specific use of this routine could be to shift 
! the zeroth component of Fourier Transformed data 
! to the centre of the transformed vector. Yours 
! truly is accustomed to dealing with appropriately 
! SHIFTED data only, and gets throughly confused 
! when viewing UNSHIFTED data after a Fourier Tran-
! sform.  Hence this routine! 
! This would require the user to specify the input
! argument "nshift" as:
!               nshift = floor(N/2), whose 
! F77 implementation would be :
!               nshift = int(N/2)
!      
! USAGE:     call shift(A,N,nshift)
!            A = Input vector
!            N = Total no. of elements in A
!       nshift = No. of position to shift
!            A = Output shifted vector (A is modified)
!     
!      --wasim raja, 29 Dec, 2009
!      
!      Last updated: 
!      This code has been tested to work for:
!      1) "fft_general_lin" Fourier Transform routine
!------------------------------------------------------
chelp-     

      !! Driver program (for testing)
      !implicit none
      !integer*4 N,nshift
      !parameter(N=6)
      !integer*4 i
      !real*4 A(N), O(N)

      !!nshift = 2
      !nshift = int(N/2) ! for fftshift
      !do i = 1,N
      !   A(i) = i
      !enddo
      ! write(*,*)"Before shift: "
      !write(*,*)(A(i),i = 1,N)
      !call shift(A,N,nshift)
      ! write(*,*)"After shift: "
      !write(*,*)(A(i),i = 1,N)
      !end


      subroutine shift(A,N,nshift)

      implicit none
      integer*4 maxdim
      parameter(maxdim = 4194304) ! 4 MB worth
      real*4 A(*), O(maxdim)
      integer*4 N, nshift
      integer*4 i

      if(nshift.ge.N)then
              nshift = mod(nshift,N)
      endif
      !if(nshift.ge.N.or.nshift.lt.0)then
      if(nshift.lt.0)then
              write(*,*)"INVALID 'nshift' argument passed!"
              write(*,*)'Input array retained as it is!'
              do i = 1,N
                 O(i) = A(i)
              enddo
      else if(nshift.eq.0)then
              do i = 1,N
                 O(i) = A(i)
              enddo
      else
              do i = 1,N
                 if(i.le.nshift)then
                         O(i) = A(N - nshift + i)
                 else
                         O(i) = A(i - nshift)
                 endif
              enddo
      endif
      do i = 1,N
         A(i) = O(i)
      enddo
      return

      end

