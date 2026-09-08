chelp+ 
      !-------------------------------------------------------
      ! This code was developed to compare 2 images using 
      ! scatter plots, to see any correlation between the 
      ! two. 
      ! 
      !                                  --wr, 05 Apr, 2013
      !-------------------------------------------------------
chelp- 
      ! Last modification: wr, 05 Apr, 2013.
      !      -> Provision to compute a map of ratio of the two 
      !         maps after the gains of each map has been 
      !         calibrated out.
      !         This is done by first dividing each map by their 
      !         respective maxima -- this takes out the common 
      !         scaling factor from each of the individual maps. 
      !         Then the ratio of the two maps is computed. 
      !         Any "flaring" and/or spectral index distribution 
      !         (if the 2 maps are at different frequencies) 
      !         may be probed. 
      !                                  --wr, 09 Apr, 2013 
      !
      !---------------------------------------------------------


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
      character*1       junkchar, yorn  
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 

      integer*4         naxes(max_axes), naxis, 
     -                  fpixels(max_axes), lpixels(max_axes)
      real*4            xcutoff1, xcutoff2, scale_fac 

      real*4            image1(maxdimx,maxdimy), 
     -                  image2(maxdimx,maxdimy), 
     -                  image3(maxdimx,maxdimy), 
     -                  xarr(maxdimx*maxdimy), 
     -                  yarr(maxdimx*maxdimy), 
     -                  tmp_arr(maxdimx*maxdimy) 
      real*4            xmaxval, ymaxval
      integer*4         xmaxpix1, ymaxpix1, 
     -                  xmaxpix2, ymaxpix2 
      integer*4         xpix(maxdimx*maxdimy), ypix(maxdimx*maxdimy) 
      integer*4         iunit, iunit1, iunit2 
      integer*4         status 
      character*120         templine, templine1, templine2 
      real*4            xmin, xmax, ymin, ymax 
      character         xlabel*120, ylabel*120, title*120 

      real*4            A, B 

      integer*4         rwmode, blocksize, group 

      
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

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      xlabel = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      ylabel = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      title = templine(1:nchar(templine)) 

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

      status = 0 

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
      !------------------------------------
      ! Locate the position of max in the two images: 
      xmaxval = xarr(1) 
      xmaxpix1 = xpix(1)  
      ymaxpix1 = ypix(1)  

      ymaxval = yarr(1) 
      ymaxpix2 = xpix(1)  
      ymaxpix2 = ypix(1)  
      do i = 1,ngood
         if(xarr(i) .gt. xmaxval)then
                 xmaxval = xarr(i) 
                 xmaxpix1 = xpix(i) 
                 ymaxpix1 = ypix(i) 
         endif
         if(yarr(i) .gt. ymaxval)then
                 ymaxval = yarr(i) 
                 xmaxpix2 = xpix(i) 
                 ymaxpix2 = ypix(i) 
         endif
      enddo
      write(*,*)"    maximum value of image 1:",xmaxval 
      write(*,*)"maxima of image 1 located at: ",xmaxpix1,ymaxpix1 
      write(*,*)"    maximum value of image 2:",ymaxval 
      write(*,*)"maxima of image 2 located at: ",xmaxpix2,ymaxpix2 

      call minima(xarr,ngood,xmin)
      call maxima(xarr,ngood,xmax)
      call minima(yarr,ngood,ymin)
      call maxima(yarr,ngood,ymax)
      !xmax = 0.05 
      write(*,*)"    maximum value of image 1 again:",xmax 
      write(*,*)"    maximum value of image 2 again:",ymax 

      ! Some plotting 
      call pgbeg(0,'1/xs',1,1)
!      call pgbeg(0,
!     - 'degree_of_polarization_of_dominant_RM_components.ps/cps',1,1)
      !call pgbeg(0,'source_noise_assessment_casa_10chan.ps/cps',1,1)
!      call pgbeg(0,'galactic_diffused_emission_towards_3C468-1.ps/cps',
!     -             1,1)
      call pgask(.false.)
      call pgsci(1)
      call pgsch(1.6) 

      call pgenv(xmin,xmax,ymin,ymax,0,1)

      call pglabel(xlabel,ylabel,title)
      !call pgline(ngood,xarr,yarr) 
      call pgsci(8)
      call pgpt(ngood,xarr,yarr,1) 

      call pgsci(1)

      call fit_linear(xarr,yarr,ngood,A,B)
      write(*,*)"+++++++++++++++++++++++++"
      write(*,*)"Y = ",A, " + ", B,"*X"

      do iy = 1,ngood
         tmp_arr(iy) = A + B*xarr(iy) 
      enddo

      call pgsci(1)
      write(*,*)"Do you wish a st. line fit to data (y/n): "
      read(*,*)yorn 
      if (yorn.eq.'y'.or.yorn.eq.'Y')then 
              call pgline(ngood,xarr,tmp_arr) 

              write(templine1,'(f6.3)')A 
              templine1 = "Y = "//templine1(1:nchar(templine1))//" + "
              write(templine2,'(f6.3)')B 
              templine2 = " "//templine2(1:nchar(templine2))//" * X"

              templine = templine1(1:nchar(templine1))//
     -                   templine2(1:nchar(templine2))

              write(*,*)templine(1:nchar(templine))
              call pgsch(1.6)
              call pgmtxt ('T', 0.2, 0.4, 0.0,templine)
      endif

      call pgend 

      call mean(xarr,ngood,A)
      write(*,*)"Average Flux density: ",A, "Jy/Beam "

!      do i = 1,ngood
!         write(77,*)xarr(i), yarr(i) 
!      enddo
!      call pgbeg(0,'2/xs',0,0)
!      call pghist(ngood,xarr,-250.0,250.0,50,0)
!      call pgend 

      !===========================================
      ! Some experiments to write various ata out: 
      ! Time to write output fits image: 
      xmaxval = 1.0 
      ymaxval = 1.0 
      do i = 1,ngood 
         xarr(i) = xarr(i)/xmaxval 
         yarr(i) = yarr(i)/ymaxval 

         ix = xpix(i)
         iy = ypix(i) 
         if(xarr(i).eq.0.0)then 
                 image3(ix,iy) = 0.0  
         else
                 image3(ix,iy) = yarr(i)/xarr(i) ! if ratio image wanted 
                 !image3(ix,iy) = yarr(i)-xarr(i) ! if diff image wanted 
                 !A = real(ix-xmaxpix2)          
                 !B = real(iy-ymaxpix2) 
                 !image3(ix,iy) = sqrt(A**2+B**2) ! if image of the distance 
                                                  ! from center of 2nd
                                                  ! image wanted 
         endif
      enddo 
      !===========================================
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
      call ftmkls(iunit2,"BUNIT","frac","Units of Pixel Data",status)

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


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
