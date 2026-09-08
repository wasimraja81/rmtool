chelp+ 
      !-------------------------------------------------------
      ! This code was developed to perform various operations 
      ! on 2 2D images . 
      ! 
      !                                  --wr, 09 Jul, 2012
      !-------------------------------------------------------
chelp- 
      ! Last modification: wr, 09 Jul, 2012.
      !
      !---------------------------------------------------------


      implicit none 

      integer*4         max_axes, maxdimx, maxdimy, maxunit 
      parameter         (max_axes=99,maxdimx = 4096, maxdimy = 4096,
     -                   maxunit = 99 )

      integer*4         active_units(maxunit)
      integer*4         nchar 
      integer*4         ngood  
      real*4            amean, arms, thresh 
      real*4            A1, A2 
      integer*4         ix, iy, i  
      character*220     infile_1, infile_2, cfgfile 
      character*220     outfile, outfile2, outfile3, outfile4 
      character*220     path, out_path 
      character*1       junkchar 
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 

      integer*4         naxes(max_axes), naxis, 
     -                  fpixels(max_axes), lpixels(max_axes) 

      real*4            image1(maxdimx,maxdimy), 
     -                  image2(maxdimx,maxdimy), 
     -                  image3(maxdimx,maxdimy), 
     -                  tmp_arr(maxdimy) 
      integer*4         iunit, iunit1, iunit2, iunit3, iunit4, iunit5 
      integer*4         status, rwmode, blocksize, group 
      character*4       optype 
      character         out_name*72, out_class*16 
      character         templine*220 

      
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
      optype = templine(1:nchar(templine)) 

      read(iunit,*)A1, A2 

      !---------------------------------------------------
      ! Default OUTCLASS using OPTYPE (overwritten if 
      ! secified in .cfg file): 
      if(optype.eq.'SUM'.or.optype.eq.'sum')then
              out_class = 'SUM'
              optype = 'SUM'
      else if(optype.eq.'DIFF'.or.optype.eq.'diff')then
              out_class = 'DIFF'
              optype = 'DIFF'
      else if(optype.eq.'MULT'.or.optype.eq.'mult')then
              out_class = 'MULT'
              optype = 'MULT'
      else if(optype.eq.'DIV'.or.optype.eq.'div')then
              out_class = 'DIV'
              optype = 'DIV'
      else if(optype.eq.'LIN'.or.optype.eq.'lin')then
              out_class = 'LIN_COMB'
              optype = 'LIN'
      else
              write(*,*)"Unknown OPTYPE: ",optype(1:nchar(optype))
              write(*,*)"Quitting now..."
              stop 
      endif
      !---------------------------------------------------

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      out_path = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      out_name = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 

      ! Overwrite outclass if specified in .cfg file: 
      if(nchar(templine).gt.0)then 
              out_class = templine(1:nchar(templine)) 
      endif


      !call my_close(iunit,active_units)
      close(iunit) 


      naxis = 2   ! Number of axes  (2 for images) 


      infile_1 = path(1:nchar(path))//infile_1(1:nchar(infile_1))
      infile_2 = path(1:nchar(path))//infile_2(1:nchar(infile_2))

      write(*,*)"infile_1: ",infile_1(1:nchar(infile_1))
      write(*,*)"infile_2: ",infile_2(1:nchar(infile_2))
      write(*,*)" "

      outfile = out_path(1:nchar(out_path))// 
     -          out_name(1:nchar(out_name))//'.'//
     -          out_class(1:nchar(out_class))//
     -          '.FITS'
      outfile2 = out_path(1:nchar(out_path))// 
     -          out_name(1:nchar(out_name))//'.'//
     -          out_class(1:nchar(out_class))//
     -          '.I_LPOL.TXT'
      outfile3 = out_path(1:nchar(out_path))// 
     -          out_name(1:nchar(out_name))//'.'//
     -          out_class(1:nchar(out_class))//
     -          '.I_LPOL_NOISE_REGION.TXT'
      outfile4 = out_path(1:nchar(out_path))// 
     -          out_name(1:nchar(out_name))//'.'//
     -          out_class(1:nchar(out_class))//
     -          '.I_LPOL_ALL_REGION.TXT'

      write(*,*)"outfile: ",outfile(1:nchar(outfile))
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
      
   
      amean = 0.0 
      ngood = 0 
      if (optype.eq.'SUM')then
              do ix = 1,nxpix
                do iy = 1,nypix
                  image3(ix,iy) = image1(ix,iy) + image2(ix,iy) 
                  ! Keep the mean of good points computed (will be handy): 
                  if(image3(ix,iy).eq.image3(ix,iy))then
                          amean = amean + image3(ix,iy) 
                          ngood = ngood + 1 
                  endif
                enddo 
              enddo 
      else if (optype.eq.'DIFF')then
              do ix = 1,nxpix
                do iy = 1,nypix
                  image3(ix,iy) = image1(ix,iy) - image2(ix,iy) 
                  ! Keep the mean of good points computed (will be handy): 
                  if(image3(ix,iy).eq.image3(ix,iy))then
                          amean = amean + image3(ix,iy) 
                          ngood = ngood + 1 
                  endif
                enddo 
              enddo 
      else if (optype.eq.'LIN')then
              do ix = 1,nxpix
                do iy = 1,nypix
                  image3(ix,iy) = A1*image1(ix,iy) + A2*image2(ix,iy) 
                  ! Keep the mean of good points computed (will be handy): 
                  if(image3(ix,iy).eq.image3(ix,iy))then
                          amean = amean + image3(ix,iy) 
                          ngood = ngood + 1 
                  endif
                enddo 
              enddo 
      else if (optype.eq.'MULT')then
              do ix = 1,nxpix
                do iy = 1,nypix
                  image3(ix,iy) = image1(ix,iy) * image2(ix,iy) 
                  ! Keep the mean of good points computed (will be handy): 
                  if(image3(ix,iy).eq.image3(ix,iy))then
                          amean = amean + image3(ix,iy) 
                          ngood = ngood + 1 
                  endif
                enddo 
              enddo 
      else if (optype.eq.'DIV')then
              do ix = 1,nxpix
                 do iy = 1,nypix
                   image3(ix,iy) = image1(ix,iy) / image2(ix,iy) 
                   ! Keep the mean of good points computed (will be handy): 
                   if(image3(ix,iy).eq.image3(ix,iy))then
                           amean = amean + image3(ix,iy) 
                           ngood = ngood + 1 
                   endif
                 enddo 
               enddo 
      endif
      amean = amean/real(ngood) 
      ! Now compute the rms of good pixels: 
      arms = 0.0 
      do ix = 1,nxpix
         do iy = 1,nypix
           if(image3(ix,iy).eq.image3(ix,iy))then
                   arms = arms + (image3(ix,iy) - amean)**2 
           endif
         enddo 
      enddo 
      arms = sqrt(arms/real(ngood)) 
      write(*,*)"mean of output: ",amean 
      write(*,*)" rms of output: ",arms 

      !------------------------------------
      ! TODO: 
      ! Make provision to BLANK outliers: 
      ! 
      !------------------------------------

      ! Now dump the data into the output file: 
      iunit2 = iunit1 + 1   
      iunit3 = iunit2 + 1   
      iunit4 = iunit3 + 1   
      iunit5 = iunit4 + 1   

      open(iunit3,file=outfile2,status='unknown')
      write(iunit3,*)"# Stokes-I   Lpol Intensity (RM=0)  dPOL  RA  Dec"
      write(iunit3,*)"  "

      open(iunit4,file=outfile3,status='unknown')
      write(iunit4,*)"# Stokes-I   Lpol Intensity (RM=0)  dPOL  RA  Dec"
      write(iunit4,*)"# Noise pixels  "

      open(iunit5,file=outfile4,status='unknown')
      write(iunit5,*)"# Stokes-I   Lpol Intensity (RM=0)  dPOL  RA  Dec"
      write(iunit5,*)"# ALL pixels  "

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

      !thresh = arms*20.0 
      thresh = 0.0 
   
      do ix = 1,nxpix
         if(mod(ix-1,100).eq.0)then
                 write(*,*)"Doing x-plane: ",ix
         endif
         do iy = 1,nypix
            !if(abs(image3(ix,iy) - amean).le.thresh)then 
            if(image2(ix,iy) .ge.thresh)then ! Use filter criterion on I-image 
                    tmp_arr(iy) = image3(ix,iy) 
                    write(iunit3,*)image2(ix,iy),
     -                             image1(ix,iy),
     -                             image3(ix,iy), 
     -                             real(ix), real(iy) 
            else
                    tmp_arr(iy) = 0.0 
                    ! Ignore pixels producing 
                    ! NaN in dpol
                    if(image3(ix,iy).ne.image3(ix,iy).or.
     -                 image1(ix,iy).gt.image2(ix,iy))then 
                            image3(ix,iy) = 1.2 
                    endif
                    write(iunit4,*)image2(ix,iy),
     -                             image1(ix,iy),
     -                             image3(ix,iy), 
     -                             real(ix), real(iy) 
            endif
            ! Now write all pixels together in a file: 
            write(iunit5,*)image2(ix,iy),
     -                     image1(ix,iy),
     -                     image3(ix,iy), 
     -                     real(ix), real(iy) 
         enddo
         !write(*,*)(tmp_arr(iy),iy=1,nypix)
   
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
      write(*,*)"final status: ",status

      close(iunit3) 
      close(iunit4) 
      close(iunit5) 

      write(*,*)" "
      write(*,*)" OUTFILES WRITTEN IN: ",out_path(1:nchar(out_path))

      end


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
