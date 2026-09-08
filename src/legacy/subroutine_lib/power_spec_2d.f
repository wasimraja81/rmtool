chelp+
      !-----------------------------------------------------
      ! This subroutine computes the power spectra given an 
      ! image. It uses the FX scheme for the computation of 
      ! the power spectra and hence is a fast routine. 
      !
      ! The dependencies include: 
      ! 1) fft_general_lin.f (central fft-engine from Desh)
      ! 2) fft2d.f           (tailor made to perform 2D FT)
      ! 3) fftshift2d.f      (shift operation for nice )
      ! 4) fort_lib.f        (collection of my subroutines)
      !
      !                       --wr, 01,June, 2011      
      !-----------------------------------------------------
chelp-

      subroutine power_spec_2d(InData,pspecData,acorData,nx,ny,maxdimx,
     -                         maxdimy)


      implicit none 

      integer*4   nx, ny, maxdimx, maxdimy 
      real*4      InData(maxdimx,maxdimy), pspecData(maxdimx,maxdimy), 
     -            acorData(maxdimx,maxdimy)
      real*4      realData(nx,ny), imagData(nx,ny) 
      integer*4   i, j  
      real*4      atmp 

      ! Compute the mean (to be removed from data): 
      atmp = 0.0 
      do i = 1,nx
         do j = 1,ny
            atmp = atmp + InData(i,j)
         enddo
      enddo
      atmp = atmp/real(nx*ny)

      ! Fill up the real and imaginary part of the InArray: 
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = InData(i,j) - atmp 
            imagData(i,j) = 0.0
         enddo
      enddo
!      !-----------------------------------------------
!      ! TEST: 
!      do i = 1,nx
!         write(*,*)(realData(i,j),j = 1,ny),'      ',
!     -             (imagData(i,j),j = 1,ny)
!      enddo
!      !-----------------------------------------------

      ! The F-part: 
      call fft2d(realData,imagData,nx,ny,nx,ny)

      ! Normalise the Fourier Transformed data : 
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = realData(i,j)/real(nx*ny)
            imagData(i,j) = imagData(i,j)/real(nx*ny)
         enddo
      enddo

!      !-----------------------------------------------
!      ! TEST: 
!      do i = 1,nx
!         write(*,*)(realData(i,j),j = 1,ny),'      ',
!     -             (imagData(i,j),j = 1,ny)
!      enddo
!      !-----------------------------------------------

      ! The X-part: 
      do i = 1,nx
         do j = 1,ny
            pspecData(i,j) = realData(i,j)*realData(i,j) + 
     -                     imagData(i,j)*imagData(i,j) 
         enddo
      enddo

      ! Perform fftshift to make the data visibly appealing: 
      call fftshift2d(pspecData,nx,ny,maxdimx,maxdimy)
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = pspecData(i,j)
            imagData(i,j) = 0.0 
         enddo
      enddo
      call ifft2d(realData,imagData,nx,ny,nx,ny)
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = realData(i,j)*real(nx*ny)
            imagData(i,j) = imagData(i,j)*real(nx*ny)
            acorData(i,j) = sqrt(realData(i,j)**2 + 
     -                           imagData(i,j)**2)
!            acorData(i,j) = realData(i,j) 
         enddo
      enddo
      call fftshift2d(acorData,nx,ny,maxdimx,maxdimy)


      return
      end


