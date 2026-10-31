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
     -                  maxunit = 99)

      integer*4         active_units(maxunit)
      integer*4         nchar 
      real*4            thresh 
      integer*4         ix, iy, i  
      character*220     infile_1, infile_2, cfgfile 
      character*220     outfile
      character*220     path, out_path 
      character*1       junkchar 
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 

      integer*4         naxis

      real*4            image1(maxdimx,maxdimy), 
     -                  image2(maxdimx,maxdimy) 
      integer*4         iunit, iunit1, iunit2, iunit3  
      integer*4         status   
      character*4       optype 
      character         out_name*16, out_class*32 
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

      read(iunit,*)thresh 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      optype = templine(1:nchar(templine)) 

      !---------------------------------------------------
      ! Default OUTCLASS using OPTYPE (overwritten if 
      ! secified in .cfg file): 
      if(optype.eq.'SCAT'.or.optype.eq.'scat')then
              out_class = 'SCATTER'
              optype = 'SCATTER'
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
     -          '.TXT'

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

   

      ! Now dump the data into the output file: 
      iunit2 = iunit1 + 1   

      open(iunit3,file=outfile,status='unknown')
      write(iunit3,*)"# RA   Dec   Stokes-I   RMS(I)"
      write(iunit3,*)"  "

   
      do ix = 1,nxpix
         if(mod(ix-1,100).eq.0)then
                 write(*,*)"Doing x-plane: ",ix
         endif
         do iy = 1,nypix
            if(image2(ix,iy) .ge.thresh)then ! Use filter criterion on I-image 
                    write(iunit3,*)ix, iy, image2(ix,iy),
     -                             image1(ix,iy) 
            endif
         enddo
      enddo
   
      close(iunit3) 

      write(*,*)" "
      write(*,*)" OUTFILES WRITTEN IN: ",out_path(1:nchar(out_path))

      end


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
