chelp+
!-----------------------------------------------------------
! This subroutine is meant to align Fourier Transformed
! data such that the 0-th component is shifted to the 
! central pixel of the transformed array. 
! Currently this code is used for 1D and 2D data only.
!
! Required: 
!      1) shift.f
!      --wasim raja, 06 Jan, 2010
!-----------------------------------------------------------
chelp-      
      !-----------------------------------------------------
      ! Last modified: include file dependency taken away. 
      !                Now the dimensions of the arrays are 
      !                also passed to the subroutines. 
      !                              wr, 01 June, 2011
      !
      !-----------------------------------------------------

      subroutine fftshift2d(InArr,N1,N2,dimx,dimy)

      implicit none

      integer*4  dimx, dimy 
      real*4 InArr(dimx,dimy)
      integer*4 N1, N2, nshift1, nshift2 
      integer*4 i, j
      real*4 tmp(dimx+dimy)

      nshift1 = int(N1/2)
      nshift2 = int(N2/2)

      ! Shift along columns for a given row:
      do i = 1,N1
          do j = 1,N2
             tmp(j) = InArr(i,j)
          enddo
          call shift(tmp,N2,nshift2)
          do j = 1,N2
             InArr(i,j) = tmp(j)
          enddo
      enddo

      ! Shift along rows now for a given column:
      do i = 1,N2
          do j = 1,N1
             tmp(j) = InArr(j,i)
          enddo
          call shift(tmp,N1,nshift1)
          do j = 1,N1
             InArr(j,i) = tmp(j)
          enddo
      enddo
      return
      end

