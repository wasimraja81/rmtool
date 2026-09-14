chelp+
      !-------------------------------------------------------------
      ! This code computes AUTO-CORRELATION of polarized vectors 
      ! (LPOL & POLA) FITS images.  
      ! 
      !                                    -- wr, 14 Nov, 2011
      !-------------------------------------------------------------
chelp-



      implicit none

      !============================================================
      ! dimensions: 
      integer*4   xdim, ydim  
      parameter   (xdim = 4096, ydim = 4096) 

      integer*4   cxpix, cypix, nxpix, nypix 
      integer*4   nchar 
      integer*4   i, ix, iy 
      ! arrays: 
      real*4      Q_arr(xdim,ydim), U_arr(xdim,ydim)
      real*4      acorr_arr(xdim,ydim), tmp_arr(xdim*ydim)

      character   infileQ*220, infileU*220, outfile*220 

      !============================================================
      ! For FITS files: 
      integer*4          naxis, naxes(2), fpixels(2), lpixels(2) 
      integer*4          group, blocksize, rwmode, decimals, status 
      integer*4          iunit1, iunit2

      character          bunit*16, ctype3*16  
      real*4             bscale 
      !============================================================
      real*4      atmp, ptmp, deg2rad, pi
      
      pi = acos(-1.0) 
      deg2rad = pi/180.0 

      bscale = 1.0 
      bunit = '(Jy/Bm)**2'
      ctype3 = 'AutoCorr'

!-------------------------------------------------------------------
      if(iargc().lt.2)then
              write(*,*)'  '
              write(*,*)' Usage: '
              write(*,*)'> autocorr_image <infileQ> <infileU>'
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)'  '
              stop
      else if(iargc().eq.2)then
              call getarg(1,infileQ)
              infileQ = infileQ(1:nchar(infileQ))

              call getarg(2,infileU)
              infileU = infileU(1:nchar(infileU))
      endif

      ! Attempt to find the extension (fits or FITS? )
      if (index(infileQ,'.fits').gt.0)then
              i = index(infileQ,'.fits') 
              outfile(1:) = infileQ(1:i-1)//'.acorr.fits'
      else if (index(infileQ,'.FITS').gt.0)then
              i = index(infileQ,'.FITS') 
              outfile(1:) = infileQ(1:i-1)//'.acorr.fits'
      else
              i = nchar(infileQ) + 1  
              outfile(1:) = infileQ(1:i-1)//'.acorr.fits'
      endif 




      !=======================================================
      ! Initialise STATUS to zero:
      status = 0

      ! Read the infiles: 
      call load_fits_image(infileQ,cxpix,cypix,nxpix,nypix,
     -                           Q_arr,xdim,ydim,status)
      call load_fits_image(infileU,cxpix,cypix,nxpix,nypix,
     -                           U_arr,xdim,ydim,status)

      ! Generate the Q and U arrays from the amplitude and 
      ! position angle images: 
      do ix = 1,nxpix
         do iy = 1,nypix 
            ptmp = Q_arr(ix,iy) 
            atmp = U_arr(ix,iy)*2.0*deg2rad  ! conversion from POLA  

            Q_arr(ix,iy) = ptmp*cos(atmp) 
            U_arr(ix,iy) = ptmp*sin(atmp) 

            !Q_arr(ix,iy) = ptmp
            !U_arr(ix,iy) = 0.0 !ptmp*sin(atmp) 
         enddo      
      enddo


      !=======================================================
      ! Main task of the program begins here...

      !----------------------------------------------------
      ! Compute the 2D Auto correlation :  
      call acorr_2d(Q_arr,U_arr,acorr_arr,
     -                   nxpix,nypix,xdim,ydim)
      !--------------------------------------------------

      ! Now write out FITS files 
      ! Use a reference file for copying the headers: 
      iunit1 = 11  
      rwmode = 0 
      blocksize = 0 
      group = 1 
      decimals = 11 
      call FTOPEN(iunit1,infileQ,rwmode,blocksize,status)

      iunit2 = iunit1 + 1
      call FTINIT(iunit2,outfile,blocksize,status) 
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      ! Copy the entire header from one of the inputs 
      ! to the output file: 
      call FTCPHD(iunit1,iunit2,status) 

      ! Modify the BUNIT in the header: 
      call ftmkls(iunit2,"ctype3",ctype3,"Physical Quantity",status)
      call ftmkls(iunit2,"BUNIT",bunit,"Units of Pixel Data",status)


      naxis = 2 
      naxes(1) = nxpix 
      naxes(2) = nypix 

      do ix = 1,nxpix
         if(mod(ix-1,100).eq.0)then
                 write(*,*)"Doing x-plane: ",ix
         endif
         do iy = 1,nypix
            tmp_arr(iy) = acorr_arr(ix,iy) 
         enddo
         fpixels(1) = ix 
         lpixels(1) = ix 
   
         fpixels(2) = 1 
         lpixels(2) = nypix
         !--------------------------------------------------
         ! Write the FITS FILE now: 
         call ftpsse(iunit2,group,naxis,naxes,fpixels,lpixels,
     -               tmp_arr,status)
      enddo
   
      call FTCLOS(iunit1,status)
      call FTCLOS(iunit2,status)
      !------------------------------------------------




      end


      include '/usr/lib/subroutine_lib/nchar.f'
      include 'load_fits_image.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/my_randn.f'
      include '/usr/lib/subroutine_lib/shift.f'
      include '/usr/lib/subroutine_lib/FFT/fft_general_lin.f'
      include '/usr/lib/subroutine_lib/FFT/fft2d.f'
      include '/usr/lib/subroutine_lib/FFT/ifft2d.f'
      include '/usr/lib/subroutine_lib/FFT/fftshift2d.f'

!========================================================================
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

      subroutine acorr_2d(InData_re,InData_im,acorData,nx,ny,maxdimx,
     -                         maxdimy)


      implicit none 

      integer*4   nx, ny, maxdimx, maxdimy 
      real*4      InData_re(maxdimx,maxdimy),InData_im(maxdimx,maxdimy),
     -            acorData(maxdimx,maxdimy)
      real*4      realData(nx,ny), imagData(nx,ny), 
     -            pspecData(maxdimx,maxdimy) 
      integer*4   i, j  
      real*4      atmp, btmp, pi
      integer*4   nzero, n_positive, mx, my, nwin_x, nwin_y 
      real*4      filter_x(maxdimx), filter_y(maxdimy), filter  
      real*4      dx, dy, tmp_arr(ny), rand_arr(nx,ny) 

!      ! Test : Put random numbers
!      do i = 1,nx
!         call my_randn(ny,0.0,1.0,tmp_arr) 
!         do j = 1,ny
!            rand_arr(i,j) = tmp_arr(j) + 
!     -      exp(-((real(i-257))**2 + (real(j-257))**2)/(2.0*256*256))
!         enddo 
!         ! Test : Put random numbers 
!         realData(i,j) = rand_arr(i,j)  
!         imagData(i,j) = 0.0   
!      enddo
      !------------------------------------------
      ! Compute the mean (to be removed from data): 
      atmp = 0.0 
      btmp = 0.0 
      pi = acos(-1.0) 

      do i = 1,nx
         do j = 1,ny
            atmp = atmp + InData_re(i,j)
            btmp = btmp + InData_im(i,j)
         enddo
      enddo
      atmp = atmp/real(nx*ny)
      btmp = btmp/real(nx*ny)

      ! Fill up the real and imaginary part of the InArray: 
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = InData_re(i,j) - atmp 
            imagData(i,j) = InData_im(i,j) !- btmp 
         enddo
      enddo
      !-----------------------------------------------

      ! The F-part: 
      call fft2d(realData,imagData,nx,ny,nx,ny)

      ! Normalise the Fourier Transformed data : 
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = realData(i,j)/real(nx*ny)
            imagData(i,j) = imagData(i,j)/real(nx*ny)
         enddo
      enddo


      ! Add a filter to suppress high frequencies: 
      ! n = window size in angular domain (in pixels)
      ! Nf = Total number of pixels 
      !
      ! m = N/n : Window size (in pixels) in Fourier Domain 
      !==================================== 
      ! Initialization: 
      do i = 1,nx
         filter_x(i) = 0.0 
         !filter_x(i) = 1.0 
      enddo
      do i = 1,ny
         filter_y(i) = 0.0 
         !filter_y(i) = 1.0 
      enddo
      !==================================== 
      nwin_x = 1; 
      nwin_y = 1; 
      mx = int(nx/nwin_x); 
      if(mod(mx,2).eq.0)then
              n_positive = mx/2 
      else
              n_positive = (mx+1)/2 - 1 
      endif
      if(mod(nx,2).eq.0)then
              nzero = nx/2 + 1 
      else
              nzero = (nx+1)/2 
      endif

      dx = pi/real(mx) 
      do i = 1,n_positive
         atmp = real(i)*dx 
         filter_x(nzero+i) = cos(atmp)*cos(atmp)
         filter_x(nzero-i) = cos(-atmp)*cos(-atmp)

         !filter_x(nzero+i) = -1.0*cos(atmp)*cos(atmp) + 1.0 
         !filter_x(nzero-i) = -1.0*cos(-atmp)*cos(-atmp) + 1.0 
         write(*,*)"filter_x: ",filter_x(nzero+i), nzero+i,atmp*180./pi 
      enddo
      filter_x(nzero) = cos(0.0) 

      my = int(ny/nwin_y); 
      if(mod(my,2).eq.0)then
              n_positive = my/2 
      else
              n_positive = (my+1)/2 - 1 
      endif

      if(mod(ny,2).eq.0)then
              nzero = ny/2 + 1 
      else
              nzero = (ny + 1)/2
      endif
      dy = pi/real(my) 
      do i = 1,n_positive
         atmp = real(i)*dy 
         filter_y(nzero+i) = cos(atmp)*cos(atmp) 
         filter_y(nzero-i) = cos(-atmp)*cos(-atmp) 

         !filter_y(nzero+i) = -1.0*cos(atmp)*cos(atmp) + 1.0 
         !filter_y(nzero-i) = -1.0*cos(-atmp)*cos(-atmp)+ 1.0  
      enddo
      filter_y(nzero) = cos(0.0) 

      do i = 1,nx
        write(87,*)filter_x(i), filter_y(i) 
      enddo

      !-----------------------------------------------
      ! Perform fftshift to make the data visibly appealing: 
      !call fftshift2d(realData,nx,ny,nx,ny)
      !call fftshift2d(imagData,nx,ny,nx,ny)

      ! The X-part: 
      do i = 1,nx
         do j = 1,ny
            realData(i,j) = realData(i,j) 
            imagData(i,j) = imagData(i,j) 

            pspecData(i,j) = realData(i,j)*realData(i,j) + 
     -                     imagData(i,j)*imagData(i,j) 

         enddo
      enddo

      !! Perform fftshift to make the data visibly appealing: 
      call fftshift2d(pspecData,nx,ny,maxdimx,maxdimy)
      do i = 1,nx
         do j = 1,ny
            filter = filter_x(i)**2 + filter_y(j)**2
!            filter = exp(-((real(i - 257))**2 + (real(j - 257))**2)/
!     -                             (2.0*256*256))
            write(88,*)filter 
            realData(i,j) = pspecData(i,j) * filter 
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
         enddo
      enddo
      call fftshift2d(acorData,nx,ny,maxdimx,maxdimy)


      return
      end

