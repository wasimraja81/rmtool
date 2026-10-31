chelp+ 
      !==========================================
      ! This code is intended to convert a FITS 
      ! image into an ASCII file. 
      !
      !             --wr, 17 Apr, 2013
      !==========================================

chelp- 


      implicit none 

      integer*4       max_axes, maxdimx, maxdimy, maxunit 
      parameter       (max_axes=99,maxdimx = 4096, maxdimy = 4096,
     -                 maxunit = 99)


      integer*4       ix, iy 

      character*220     infile, cfgfile, outfile 
      character*220     path, outpath 
      character*1       junkchar 
      character*220     templine, tag 
      integer*4         cxpix, cypix, nxpix, nypix 
      integer*4         nchar 


      real*4            image(maxdimx,maxdimy) 
      integer*4         iunit  
      integer*4         status 

      !--------------------------------------
      ! Some input parameters: 
      if(iargc().ne.1)then 
              write(*,*)"Usage: "
              write(*,*)"You can either use a config file: "
              write(*,*)"   fits_image2ascii <config file> "
              write(*,*)" "
              stop 
      else
              call getarg(1,cfgfile) 
              cfgfile = '../CONFIG/'//cfgfile(1:nchar(cfgfile))
      endif


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
      infile = path(1:nchar(path))//templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      outpath = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      outfile = templine(1:nchar(templine))
      outfile = outpath(1:nchar(outpath))//outfile(1:nchar(outfile)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      tag = templine(1:nchar(templine))

      

      !call my_close(iunit,active_units)
      close(iunit)

      status = 0 
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
      !=======================================      

      call load_fits_image(infile, cxpix,cypix,nxpix,nypix,
     -                  image, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile(1:nchar(infile))
              write(*,*)"Quitting now..."
              stop 
      endif
     
      write(*,*)"naxes(1): ", nxpix 
      write(*,*)"naxes(2): ", nypix 
      write(*,*)"----------------------------"
      
      iunit = iunit + 1 
      open(iunit,file=outfile,status='new',err=201)
202   goto 203 
201   write(*,*)"Error opening new file: ",outfile(1:nchar(outfile))
      write(*,*)"File may exist already..."
      write(*,*)"Quitting now... proceed with care..."
      stop
203   continue 
      write(iunit,*)"# RA     Dec         "//tag(1:nchar(tag))
      do ix = 1,nxpix
        do iy = 1,nypix
           write(iunit,fmt=301)real(ix), real(iy), image(ix,iy)
        enddo 
      enddo 
      write(*,*)"Total Number of pixels: ",nxpix*nypix 
      !------------------------------------
301   format(F6.1,2x,F6.1,1x,F11.4)      

      end 


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'      
