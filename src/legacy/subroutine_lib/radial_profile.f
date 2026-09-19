chelp+
      !-----------------------------------------------------
      ! This subroutine computes the azimuthally averaged 
      ! intensity of an image. 
      !                         --wr, 21 Nov, 2011
      !-----------------------------------------------------
chelp-


      subroutine radial_profile(im_arr,nx,ny,dimx,dimy, 
     -                          rad_prf,rms_prf,npts_arr,nout) 

      implicit none 
      integer*4      dimx, dimy, nx, ny 
      real*4         im_arr(dimx,dimy), rad_prf(*), rms_prf(*) 
      integer*4      npts_arr(*) , npts, ix, iy, i, rnow, rmax, k,j   
      integer*4      cx, cy, nout 
      real*4         atmp, tmp_arr(dimx+dimy,7*(dimx+dimy))
      real*4         tmp_arr2(7*(dimx+dimy))



      ! Initiate the arrays: 
      do i = 1,dimx+dimy
        rad_prf(i) = 0.0 
        npts_arr(i) = 0 
      enddo

      ! Locate the centre of the image: 
      if (mod(nx,2).eq.0)then 
              cx = nx/2 + 1 
      else 
              cx = (nx+1)/2 
      endif 

      if (mod(ny,2).eq.0)then 
              cy = ny/2 + 1 
      else 
              cy = (ny+1)/2 
      endif 

      rmax = 1 + int(sqrt((real(cx)-1.0)**2 + (real(cy)-1.0)**2)) 
      do ix = 1,nx 
         do iy = 1,ny 
            rnow = 1 + int(sqrt((real(ix) - real(cx))**2 + 
     -                          (real(iy) - real(cy))**2)) 
            rad_prf(rnow) = rad_prf(rnow) + im_arr(ix,iy) 
            npts_arr(rnow) = npts_arr(rnow) + 1 
            k = npts_arr(rnow) 
            tmp_arr(rnow,k) = im_arr(ix,iy) 
         enddo 
      enddo 

      do i = 1,rmax 
         npts = npts_arr(i) 
         rad_prf(i) = rad_prf(i)/real(npts) 
         do j = 1,npts 
            tmp_arr2(j) = tmp_arr(i,j)
         enddo 
         call rms(tmp_arr2,npts,atmp)
         rms_prf(i) = atmp 
      enddo 
      nout = rmax 

      end
         


