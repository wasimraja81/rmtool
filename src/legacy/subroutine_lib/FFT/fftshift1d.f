chelp+
! This subroutine is meant to align Fourier Transformed
! data such that the 0-th component is shifted to the 
! central pixel of the transformed array. 
! Currently this code is used for 1D only.
!
! Required: 
!      1) shift.f
!      --wasim raja, 06 Jan, 2010

      subroutine fftshift1d(InArr,N)


      real*4 InArr(*)
      integer*4 N

      nshift = int(N/2)

      call shift(InArr,N,nshift)

      return
      end

