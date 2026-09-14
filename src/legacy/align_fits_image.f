chelp+ 
      !-------------------------------------------------------
      ! This code was developed to align 2 images using 
      ! Fourier Transforms -- The phases in the FFT of the
      ! 2nd image are replaced by the phases in the FFT of 
      ! the 1st. The resultant FFT of the 2nd image is 
      ! then inverse Fourier Transformed. 
      ! 
      ! 
      !                                  --wr, 05 Apr, 2013
      !-------------------------------------------------------
chelp- 


      implicit none 

      integer*4         max_axes, maxdimx, maxdimy, maxunit 
      parameter         (max_axes=99,maxdimx = 4096, maxdimy = 4096,
     -                   maxunit = 99 )

      integer*4         active_units(maxunit)
      integer*4         nchar 
      integer*4         ngood  
      integer*4         ix, iy, i  
      character*220     infile_1, infile_2, cfgfile, outfile 
      character*220     path, outpath 
      character*1       junkchar 
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 
      integer*4         nxshift, nyshift 

      integer*4         naxes(max_axes), naxis, 
     -                  fpixels(max_axes), lpixels(max_axes)
      real*4            xcutoff1, xcutoff2, scale_fac 

      real*4            image1(maxdimx,maxdimy), 
     -                  image2(maxdimx,maxdimy), 
     -                  image3(maxdimx,maxdimy), 
     -                  xarr(maxdimx*maxdimy), 
     -                  yarr(maxdimx*maxdimy), 
     -                  tmp_arr(maxdimx*maxdimy) 
      real*4            rfx(maxdimx*maxdimy), ifx(maxdimx*maxdimy) 
      real*4            atemp, btemp, 
     -                  xpix(maxdimx*maxdimy), ypix(maxdimx*maxdimy) 
      integer*4         iunit, iunit1, iunit2 
      integer*4         status 
      character*120     templine

      integer*4         rwmode, blocksize, group 
      logical           align_by_hand
      character*1       yorn 

      
      !--------------------------------------
      ! Some input parameters: 
      if(iargc().ne.1)then 
              write(*,*)"Usage: "
              write(*,*)"You can either use a config file: "
              write(*,*)"    comb_fits_image <config file> "
              write(*,*)" "
              stop 
      else
              call getarg(1,cfgfile) 
              cfgfile = '../CONFIG/'//cfgfile(1:nchar(cfgfile))
      endif

      do i = 1,maxunit 
         active_units(i) = 0 
      enddo

      !call get_unit(iunit, active_units, maxunit ) 
      call get_lun(iunit) 
      open(iunit,file=cfgfile,status='old',err=101)
      goto 102
101   write(*,*)"Error opening file: ",cfgfile(1:nchar(cfgfile))
      write(*,*)"Quitting now..."
      stop 

102   continue 
      read(iunit,*)junkchar   ! comment line 
      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      path = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_1 = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_2 = templine(1:nchar(templine)) 

      read(iunit,*)xcutoff1, xcutoff2 
      read(iunit,*)scale_fac  
      xcutoff1 = xcutoff1*scale_fac 
      xcutoff2 = xcutoff2*scale_fac 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      outpath = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      outfile = templine(1:nchar(templine))//'.FITS'
      outfile = outpath(1:nchar(outpath))//outfile(1:nchar(outfile)) 

      
      read(iunit,*)yorn 
      read(iunit,*)nxshift, nyshift 
      if(yorn.eq.'y'.or.yorn.eq.'Y')then 
              align_by_hand = .true. 
              write(*,*)"2nd image viz. :",infile_2(1:nchar(infile_2))
              write(*,*)"will be shifted by: ",nxshift, nyshift 
              write(*,*)"pixels x and y directions respectively. "
      else
              align_by_hand = .false. 
      endif 

      !call my_close(iunit,active_units)
      close(iunit) 


      naxis = 2   ! Number of axes  (2 for images) 


      infile_1 = path(1:nchar(path))//infile_1(1:nchar(infile_1))
      infile_2 = path(1:nchar(path))//infile_2(1:nchar(infile_2))

      write(*,*)"infile_1: ",infile_1(1:nchar(infile_1))
      write(*,*)"infile_2: ",infile_2(1:nchar(infile_2))
      write(*,*)" "

      !-------------------------------------- 
      ! We wish to load the entire images into 
      ! memory (Ensure that the images are NOT 
      ! very large): 
      !=======================================
      ! Freeze parameters to match criterion for 
      ! ENTIRE image reading: 
      cxpix = 0 
      cypix = 0 
      
      nxpix = 0 
      nypix = 0 
      nxpix1 = 0 
      nypix1 = 0 
      nxpix2= 0 
      nypix2= 0 
      !=======================================


      call load_fits_image(infile_1, cxpix1,cypix1,nxpix1,nypix1,
     -                  image1, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_1(1:nchar(infile_1))
              write(*,*)"Quitting now..."
              stop 
      endif
   
      call load_fits_image(infile_2, cxpix2,cypix2,nxpix2,nypix2,
     -                  image2, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_2(1:nchar(infile_2))
              write(*,*)"Quitting now..."
              stop 
      endif
   
      
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix1.ne.nxpix2.or.nypix1.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      else
              nxpix = nxpix1 
              nypix = nypix1 
      endif

      ! TEST: 
      write(*,*)"naxes(1): ", nxpix 
      write(*,*)"naxes(2): ", nypix 
      write(*,*)"----------------------------"

   
      ngood = 0 
      do ix = 1,nxpix
        do iy = 1,nypix
          image1(ix,iy) = image1(ix,iy) * scale_fac  
          image2(ix,iy) = image2(ix,iy) * scale_fac  
          ! Also initialize the output image array: 
          image3(ix,iy) = 0.0 
          ! Keep the mean of good points computed (will be handy): 
          if(image1(ix,iy).gt.xcutoff1.and.
     -       image2(ix,iy).gt.xcutoff2)then
             if((image1(ix,iy).eq.image1(ix,iy)).and.
     -          (image2(ix,iy).eq.image2(ix,iy)))then
                  ngood = ngood + 1 
                  xarr(ngood) = image1(ix,iy) 
                  yarr(ngood) = image2(ix,iy) 

                  ! Keep the pixel numbers stored as well 
                  xpix(ngood) = ix
                  ypix(ngood) = iy
             endif
          endif
        enddo 
      enddo 
      write(*,*)"Total Number of pixels: ",nxpix*nypix 
      write(*,*)" Number of good pixels: ",ngood
      !Play time starts: 
      ! Lets align the two images: 
      ! First lets try by hand: 
      if(align_by_hand)then
              do i = 1,ngood 
                 ix = xpix(i)
                 iy = ypix(i) 
                 image3(ix,iy) = yarr(i) 
              enddo 
              do ix = 1,nxpix 
                 do iy = 1,nypix
                    tmp_arr(iy) = image3(ix,iy)
                 enddo
                 call shift(tmp_arr,nypix,nyshift)
                 do iy = 1,nypix
                    image3(ix,iy) = tmp_arr(iy) 
                 enddo
              enddo
              do iy = 1,nypix 
                 do ix = 1,nxpix
                    tmp_arr(ix) = image3(ix,iy)
                 enddo
                 call shift(tmp_arr,nxpix,nxshift)
                 do ix = 1,nxpix
                    image3(ix,iy) = tmp_arr(ix) 
                 enddo
              enddo
              !------------------------------------
      else
              !------------------------------------
              ! Now let's try some fancy method: 
              do i = 1,ngood
                 rfx(i) = xarr(i) 
                 ifx(i) = 0.0 
              enddo
              call fft1d(rfx,ifx,ngood) 
              do i = 1,ngood
                 ! compute the phases of image1: 
                 tmp_arr(i) = atan2(ifx(i),rfx(i))  
              enddo

              ! Now Fourier Transform image2
              do i = 1,ngood
                 rfx(i) = yarr(i) 
                 ifx(i) = 0.0 
              enddo
              call fft1d(rfx,ifx,ngood) 
              do i = 1,ngood
                 ! compute the abs val and replace phases 
                 ! of image2 by those of image1: 
                 atemp = sqrt(rfx(i)**2 + ifx(i)**2)
                 btemp = atan2(ifx(i), rfx(i))
                 !rfx(i) = atemp*cos(tmp_arr(i) - btemp )
                 !ifx(i) = atemp*sin(tmp_arr(i) - btemp )
                 rfx(i) = atemp*cos(tmp_arr(i))
                 ifx(i) = atemp*sin(tmp_arr(i))
              enddo
              ! Now inverse Fourier transform to get back 
              ! a shifted 2nd image: 
              call ifft1d(rfx,ifx,ngood) 
              do i = 1,ngood
                 yarr(i) = rfx(i) 
              enddo
        
              do i = 1,ngood 
                 ix = xpix(i)
                 iy = ypix(i) 
                 image3(ix,iy) = yarr(i) 
              enddo 
      endif

      ! Time to write output fits image: 
      !---------------------------------------------
      ! Some fitsio requirements: 
      blocksize = 0   
      group = 1

      iunit1 = 11 
      rwmode = 0 
      call FTOPEN(iunit1,infile_1,rwmode,blocksize,status) 
      call FTGISZ(iunit1,max_axes,naxes,status)
   
      do i = 1,max_axes
         fpixels(i) = 0 
         lpixels(i) = 0 
      enddo
      
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
      !call ftmkls(iunit2,"BUNIT","frac","Units of Pixel Data",status)

      do ix = 1,nxpix
         if(mod(ix-1,100).eq.0)then
                 write(*,*)"Doing x-plane: ",ix
         endif
         do iy = 1,nypix
            tmp_arr(iy) = image3(ix,iy) 
         enddo
   
         fpixels(1) = ix 
         lpixels(1) = ix 
   
         fpixels(2) = 1 
         lpixels(2) = nypix 
         !--------------------------------------------------
         ! Write the FITS CUBES now: 
         call ftpsse(iunit2,group,naxis,naxes,fpixels,lpixels,
     -               tmp_arr,status)
      enddo

      call FTCLOS(iunit1,status) 
      call FTCLOS(iunit2,status)


      end


      include '/usr/lib/subroutine_lib/FFT/fft1d.f'
      include '/usr/lib/subroutine_lib/FFT/ifft1d.f'
      include '/usr/lib/subroutine_lib/FFT/fft_general_lin.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/shift.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
