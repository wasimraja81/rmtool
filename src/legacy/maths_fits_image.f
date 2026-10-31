chelp+ 
      !-------------------------------------------------------
      ! This code was developed to perform various operations 
      ! on 2D images. Currently operations pertaining to some 
      ! polarimetric calibration is allowed. However the code 
      ! has the provision of adding any number of mathematical 
      ! operations to be performed on images. 
      !                                  --wr, 02 Jun, 2012
      !-------------------------------------------------------
chelp- 
      ! Last modification: wr, 13 Jun, 2012.
      !
      !                --> Independent section to write the 
      !                    statistical properties of specific 
      !                    regions -- viz. regions with +ve 
      !                    I-pixels, or -ve I-pixels etc. 
      !                    This is achieved by keyword "nofits" 
      !                    that should be passed to the variable 
      !                    "add_req" through th econfig file. 
      !       
      !                --> Threshold in Stokes-I is now being 
      !                    used for computing the I-weighted 
      !                    map-mean of the various Stokes. 
      !                    This is to see if the leakages are 
      !                    proportional to Stokes-I ONLY for 
      !                    strong I-pixels. 
      !
      !                --> Output directory can now be different 
      !                    than the input data directory. Did 
      !                    this to avoid unmindful deletion of 
      !                    input files when the intention is 
      !                    actually to delete output files! [Wild 
      !                    cards used in delete commands are 
      !                    injurious to sleep-deprived people!] 
      !---------------------------------------------------------


      implicit none 

      integer*4         max_axes, maxdimx, maxdimy, maxunit 
      parameter         (max_axes=99,maxdimx = 4096, maxdimy = 4096,
     -                   maxunit = 99 )

      integer*4         active_units(maxunit)
      integer*4         nchar 
      integer*4         ix, iy, i, ichan 
      character*220     infileI, infileQ, infileU, infileV, cfgfile 
      character*220     outfileQ, outfileU, outfileV, summary_file 
      character*220     path, out_path 
      character*1       junkchar 
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 
      integer*4         cxpix3, cypix3, nxpix3, nypix3 
      integer*4         cxpix4, cypix4, nxpix4, nypix4 

      integer*4         naxes(max_axes), naxis, 
     -                  fpixels(max_axes), lpixels(max_axes) 

      real*4            image1(maxdimx,maxdimy), 
     -                  image2(maxdimx,maxdimy), 
     -                  image3(maxdimx,maxdimy), 
     -                  image4(maxdimx,maxdimy), 
     -                  tmp_arr_Q(maxdimy), 
     -                  tmp_arr_U(maxdimy), 
     -                  tmp_arr_V(maxdimy), 
     -                  atmp_I, atmp_Q, atmp_U, atmp_V, 
     -                          dpol_Q, dpol_U, dpol_V  
      integer*4         iunit1, iunit2, iunit3, iunit4 
      integer*4         status, rwmode, blocksize, group 
      integer*4         ineg, ipos 
      integer*4         bchan, echan, pixtype 
      integer*4         iunit  
      character         src_tag*16 
      character         I_tag*16, Q_tag*16, U_tag*16, V_tag*16 
      character         extn_tag*16,
     -                  chan_tag*3, out_tag*64 
      character         templine*220 
      character         add_req*120 
      real*8            fnow 
      real*4            thresh_I 

      
      !--------------------------------------
      ! Some input parameters: 
      if(iargc().ne.1)then 
              write(*,*)"Usage: "
              write(*,*)"    maths_fits_image <config file> "
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
      out_path = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      src_tag = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      I_tag = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      Q_tag = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      U_tag = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      V_tag = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      extn_tag = templine(1:nchar(templine)) 

      read(iunit,*)bchan, echan 
      read(iunit,*)pixtype  

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      add_req = templine(1:nchar(templine)) 

      read(iunit,*)thresh_I 

      !call my_close(iunit,active_units)
      close(iunit) 


      naxis = 2   ! Number of axes  (2 for images) 

      if(pixtype.eq.0)then
              out_tag = '.ADHOC_LCOR.ALLPIX.'
      else if(pixtype.gt.0)then
              out_tag = '.ADHOC_LCOR.POSPIX.'
      else if(pixtype.lt.0)then
              out_tag = '.ADHOC_LCOR.NEGPIX.'
      endif
      ! Write to file the statistical property of 
      ! the instrument: 
      summary_file = out_path(1:nchar(out_path))// 
     -               src_tag(1:nchar(src_tag))//
     -               out_tag(1:nchar(out_tag))// 
     -               'SUMMARY'

      call get_lun(iunit) 
      open(iunit,file=summary_file,status='unknown')
      write(iunit,*)"# Stokes-I weighted map-average Q, U and V "
      write(iunit,*)"# Freq (Hz)   <I*I>    <Q*I>    <U*I>   <V*I>" 

      do ichan = bchan,echan
         write(*,*)"Processing Channel: ",ichan 
         write(*,*)"  "
         !infile1 = 'MY_CASA81.QCL001.1.FITS' 
         !infile2 = 'MY_CASA81.ICL001.1.FITS' 
         if(ichan.lt.10)then
                 write(chan_tag,'(I1)')ichan 
         else if(ichan.ge.10.and.ichan.lt.100)then
                 write(chan_tag,'(I2)')ichan 
         else if(ichan.ge.100.and.ichan.lt.1000)then
                 write(chan_tag,'(I3)')ichan 
         else
                 write(*,*)"Too many channels!"
                 write(*,*)"Quitting now..."
                 stop 
         endif
         infileI = path(1:nchar(path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             I_tag(1:nchar(I_tag))// 
     -             '.'//
     -             extn_tag(1:nchar(extn_tag))
         infileQ = path(1:nchar(path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             Q_tag(1:nchar(Q_tag))// 
     -             '.'//
     -             extn_tag(1:nchar(extn_tag))
         infileU = path(1:nchar(path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             U_tag(1:nchar(U_tag))// 
     -             '.'//
     -             extn_tag(1:nchar(extn_tag))
         infileV = path(1:nchar(path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             V_tag(1:nchar(V_tag))// 
     -             '.'//
     -             extn_tag(1:nchar(extn_tag))

         infileI = infileI(1:nchar(infileI)) 
         infileQ = infileQ(1:nchar(infileQ)) 
         infileU = infileU(1:nchar(infileU)) 
         infileV = infileV(1:nchar(infileV)) 

         write(*,*)"infileI: ",infileI(1:nchar(infileI))
         write(*,*)"infileQ: ",infileQ(1:nchar(infileQ))
         write(*,*)"infileU: ",infileU(1:nchar(infileU))
         write(*,*)"infileV: ",infileV(1:nchar(infileV))
         write(*,*)" "

         outfileQ = out_path(1:nchar(out_path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             Q_tag(1:nchar(Q_tag))// 
     -             out_tag(1:nchar(out_tag))//
     -             extn_tag(1:nchar(extn_tag))
         outfileU = out_path(1:nchar(out_path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             U_tag(1:nchar(U_tag))// 
     -             out_tag(1:nchar(out_tag))//
     -             extn_tag(1:nchar(extn_tag))
         outfileV = out_path(1:nchar(out_path))// 
     -             src_tag(1:nchar(src_tag))//
     -             chan_tag(1:nchar(chan_tag))//
     -             '.'//
     -             V_tag(1:nchar(V_tag))// 
     -             out_tag(1:nchar(out_tag))//
     -             extn_tag(1:nchar(extn_tag))
         write(*,*)"outfileQ: ",outfileQ(1:nchar(outfileQ))
         write(*,*)"outfileU: ",outfileU(1:nchar(outfileU))
         write(*,*)"outfileV: ",outfileV(1:nchar(outfileV))
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
         nxpix3= 0 
         nypix3= 0 
         nxpix4= 0 
         nypix4= 0 
         !=======================================


         call load_fits_image(infileI, cxpix1,cypix1,nxpix1,nypix1,
     -                     image1, maxdimx, maxdimy,status)
         if(status.ne.0)then
                 write(*,*)"File for channel: ",ichan," not found..."
                 write(*,*)"Proceeding w/o channel: ",ichan 
                 status = 0 
                 goto 1919  
         endif
   
         call load_fits_image(infileQ, cxpix2,cypix2,nxpix2,nypix2,
     -                     image2, maxdimx, maxdimy,status)
         if(status.ne.0)then
                 write(*,*)"File for channel: ",ichan," not found..."
                 write(*,*)"Proceeding w/o channel: ",ichan 
                 status = 0 
                 goto 1919  
         endif
   
         call load_fits_image(infileU, cxpix3,cypix3,nxpix3,nypix3,
     -                     image3, maxdimx, maxdimy,status)
         if(status.ne.0)then
                 write(*,*)"File for channel: ",ichan," not found..."
                 write(*,*)"Proceeding w/o channel: ",ichan 
                 status = 0 
                 goto 1919  
         endif
   
         call load_fits_image(infileV, cxpix4,cypix4,nxpix4,nypix4,
     -                     image4, maxdimx, maxdimy,status)
         if(status.ne.0)then
                 write(*,*)"File for channel: ",ichan," not found..."
                 write(*,*)"Proceeding w/o channel: ",ichan 
                 status = 0 
                 goto 1919  
         endif
   
         
         ! Bare minimum check for dimensional mismatch: 
         if (nxpix1.ne.nxpix2.or.nypix1.ne.nypix2)then
                 write(*,*)"Image dimensions do not match!"
                 write(*,*)"Quitting now..."
                 stop
         else if (nxpix1.ne.nxpix3.or.nypix1.ne.nypix3)then
                 write(*,*)"Image dimensions do not match!"
                 write(*,*)"Quitting now..."
                 stop
         else if (nxpix1.ne.nxpix4.or.nypix1.ne.nypix4)then
                 write(*,*)"Image dimensions do not match!"
                 write(*,*)"Quitting now..."
                 stop
         else
                 nxpix = nxpix1 
                 nypix = nypix1 
         endif


         ! TODO: 
         ! Check if the frequency channels of the 2 infiles 
         ! match: 
   
         !---------------------------------------------
         ! Some fitsio requirements: 
         blocksize = 0   
         group = 1 
      
         iunit1 = 11  
         rwmode = 0 
         call FTOPEN(iunit1,infileI,rwmode,blocksize,status) 
         call FTGISZ(iunit1,max_axes,naxes,status)
   
         do i = 1,max_axes
            fpixels(i) = 0 
            lpixels(i) = 0 
         enddo
   
         ! TEST: 
         write(*,*)"naxes(1): ", nxpix 
         write(*,*)"naxes(2): ", nypix 
         write(*,*)"----------------------------"
         
         ! Read the frequency information from the reference FITS file
         call FTGKYD(iunit1,"CRVAL3",fnow,templine,status) 
         if (status.ne. 0)then
                 call printerror(status)
                 write(*,*)"Frequency value for chan: ",ichan,
     -                     " could not be determined from fits file"
                 write(*,*)"Quitting now..."
                 stop
         endif

         write(*,*)"Channel number: ",ichan," Freq. : ",fnow 
         !==============================================================
         ! PROBING statistical polarization properties of 
         ! various regions in the image (viz. +ve only 
         ! I-pixels, -ve only I-pixels etc. )
         if (index(add_req,'nofits').gt.0.or.
     -                 index(add_req,'NOFITS').gt.0)then
            atmp_I = 0.0 
            atmp_Q = 0.0 
            atmp_U = 0.0 
            atmp_V = 0.0 
            if (pixtype.gt.0)then 
              do ix = 1,nxpix 
                do iy = 1,nypix 
                   !if(image1(ix,iy).gt.0.0)then 
                   if(image1(ix,iy).gt.thresh_I)then 
                     atmp_I = atmp_I + image1(ix,iy)*image1(ix,iy) 
                     atmp_Q = atmp_Q + image2(ix,iy)*image1(ix,iy) 
                     atmp_U = atmp_U + image3(ix,iy)*image1(ix,iy) 
                     atmp_V = atmp_V + image4(ix,iy)*image1(ix,iy) 
                   endif
                enddo 
              enddo 
            else if (pixtype.lt.0)then 
              do ix = 1,nxpix 
                do iy = 1,nypix 
                   !if(image1(ix,iy).lt.0.0)then 
                   if(image1(ix,iy).lt.thresh_I)then 
                     atmp_I = atmp_I + image1(ix,iy)*image1(ix,iy) 
                     atmp_Q = atmp_Q + image2(ix,iy)*image1(ix,iy) 
                     atmp_U = atmp_U + image3(ix,iy)*image1(ix,iy) 
                     atmp_V = atmp_V + image4(ix,iy)*image1(ix,iy) 
                   endif
                enddo 
              enddo 
            else  ! pixtype = 0; all pixels 
               do ix = 1,nxpix 
                 do iy = 1,nypix 
                   atmp_I = atmp_I + image1(ix,iy)*image1(ix,iy) 
                   atmp_Q = atmp_Q + image2(ix,iy)*image1(ix,iy) 
                   atmp_U = atmp_U + image3(ix,iy)*image1(ix,iy) 
                   atmp_V = atmp_V + image4(ix,iy)*image1(ix,iy) 
                 enddo 
               enddo 
            endif
   
            ! skip writing FITS images and simply 
            ! write out the ASCII file containing 
            ! statistical polarization properties 
            ! of the images: 
            write(iunit,fmt=111)fnow, atmp_I, atmp_Q, atmp_U, atmp_V 
111         format (f14.3,1x,f11.3,1x,f8.3,1x,f8.3,1x,f11.3)

            call FTCLOS(iunit1,status)
            goto 1919 
         endif 
         !==============================================================
   
         atmp_I = 0.0 
         atmp_Q = 0.0 
         atmp_U = 0.0 
         atmp_V = 0.0 
   
         ineg = 0 
         ipos = 0 
         do ix = 1,nxpix
            do iy = 1,nypix
               ! Use only pixels that are > the specfied thresh in I-image) 
               ! for determination of I-weighted map-means: 
               if(image1(ix,iy).gt.thresh_I)then
                       atmp_I = atmp_I + image1(ix,iy) * image1(ix,iy) ! I(i) * I(i) 
                       atmp_Q = atmp_Q + image2(ix,iy) * image1(ix,iy) ! q(i) * I(i)  
                       atmp_U = atmp_U + image3(ix,iy) * image1(ix,iy) ! u(i) * I(i)  
                       atmp_V = atmp_V + image4(ix,iy) * image1(ix,iy) ! v(i) * I(i)  
                       ipos = ipos + 1 
               else
                       ineg = ineg + 1 
               endif 
               ! --------------------------------------------------------
               ! ignore pixels in the output image that are >/< thresh_I 
               ! depending on your choice: 
               if (pixtype.lt.0)then
                       if(image1(ix,iy).gt.thresh_I)then
                               image1(ix,iy) = 0.0 
                               image2(ix,iy) = 0.0 
                               image3(ix,iy) = 0.0 
                               image4(ix,iy) = 0.0 
                       endif
               else if (pixtype.gt.0)then
                       if(image1(ix,iy).lt.thresh_I)then
                               image1(ix,iy) = 0.0 
                               image2(ix,iy) = 0.0 
                               image3(ix,iy) = 0.0 
                               image4(ix,iy) = 0.0 
                       endif
               endif
               ! ---------------------------------------------------
            enddo 
         enddo 
         dpol_Q = atmp_Q/atmp_I 
         dpol_U = atmp_U/atmp_I 
         dpol_V = atmp_V/atmp_I 
         write(*,*)"Number of pixels < Thresh-I: ",ineg 
         write(*,*)"Number of pixels > Thresh-I: ",ipos 
         write(*,*)"   Total number of pixels: ",nxpix*nypix  
         write(*,*)"instrumental pol frac: ",dpol_Q, dpol_U, dpol_V  
         write(*,*)" "

         ! Now dump the data into the output file: 
         iunit2 = iunit1 + 1   
         call FTINIT(iunit2,outfileQ,blocksize,status) 
         if(status.ne.0)then
                 call printerror(status)
                 stop
         endif
         iunit3 = iunit2 + 1   
         call FTINIT(iunit3,outfileU,blocksize,status) 
         if(status.ne.0)then
                 call printerror(status)
                 stop
         endif
         iunit4 = iunit3 + 1   
         call FTINIT(iunit4,outfileV,blocksize,status) 
         if(status.ne.0)then
                 call printerror(status)
                 stop
         endif
   
         ! Copy the entire header from one of the inputs 
         ! to the output file: 
         call FTCPHD(iunit1,iunit2,status) 
         call FTCPHD(iunit1,iunit3,status) 
         call FTCPHD(iunit1,iunit4,status) 
   
         do ix = 1,nxpix
            if(mod(ix-1,100).eq.0)then
                    write(*,*)"Doing x-plane: ",ix
            endif
            do iy = 1,nypix
               tmp_arr_Q(iy) = image2(ix,iy) - 
     -                   dpol_Q*image1(ix,iy) 
               tmp_arr_U(iy) = image3(ix,iy) - 
     -                   dpol_U*image1(ix,iy) 
               tmp_arr_V(iy) = image4(ix,iy) - 
     -                   dpol_V*image1(ix,iy) 
            enddo
            !write(77,*)(tmp_arr(i),i=1,nypix)
   
            fpixels(1) = ix 
            lpixels(1) = ix 
   
            fpixels(2) = 1 
            lpixels(2) = nypix 
            !--------------------------------------------------
            ! Write the FITS CUBES now: 
            call ftpsse(iunit2,group,naxis,naxes,fpixels,lpixels,
     -                  tmp_arr_Q,status)
            call ftpsse(iunit3,group,naxis,naxes,fpixels,lpixels,
     -                  tmp_arr_U,status)
            call ftpsse(iunit4,group,naxis,naxes,fpixels,lpixels,
     -                  tmp_arr_V,status)
         enddo
   
         call FTCLOS(iunit1,status) 
         call FTCLOS(iunit2,status) 
         call FTCLOS(iunit3,status) 
         call FTCLOS(iunit4,status) 
         write(*,*)"final status: ",status
1919     continue 
         write(*,*)"Done doing channel: ",ichan 
      enddo    ! End of channels 

      close(iunit) 

      end


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
