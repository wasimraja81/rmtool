chelp+
      !----------------------------------------------------
      ! This code is a result of the RA-mismatch in output 
      ! RM-images. That the RA-scaling by secant-delta had 
      ! not been performed in AIPS was not obvious to me. 
      ! However there is a one-to-one correspondence with 
      ! the pixel of the spectral-images and that of the 
      ! RM-images (except in cases where a subimage had 
      ! been used). 
      !
      ! The idea here is to simply modify the header of the 
      ! RM-cubes appropriately to match those of the images
      ! in the spectral domain. 
      !                                --wr, 12 July, 2012
      !----------------------------------------------------
chelp-


      implicit none 


      character*220     fitsfile1, fitsfile2 
      integer*4         status, 
     -                  iunit1, iunit2, 
     -                  rwmode, blocksize, 
     -                  decimals 

      real*8            cval, cdelt 
      real*4            cpix 
      character*72      comment 
      integer*4         nchar 


      if(iargc().ne.2)then
              write(*,*)"Usage: "
              write(*,*)" modify_ra_dec_in_header "
              write(*,*)"        <fitsfile_in> <fitsfile_target>"
              stop
      endif


      call getarg(1,fitsfile1) 
      fitsfile1 = fitsfile1(1:nchar(fitsfile1))
      call getarg(2,fitsfile2) 
      fitsfile2 = fitsfile2(1:nchar(fitsfile2))

      write(*,*)"file 1: ",fitsfile1(1:nchar(fitsfile1))
      write(*,*)"file 2: ",fitsfile2(1:nchar(fitsfile2))

      call get_lun(iunit1) 
      write(*,*)"iunit1: ",iunit1 
      rwmode = 0 
      blocksize = 0 
      decimals = 11 
      status = 0 
      call FTOPEN(iunit1,fitsfile1,rwmode,blocksize,status) 

      rwmode = 1 
      iunit2 = iunit1 + 1 
      write(*,*)"iunit2: ",iunit2 
      call FTOPEN(iunit2,fitsfile2,rwmode,blocksize,status) 

      call ftgkyd(iunit1,"CRVAL1",cval,comment,status)
      call ftgkyd(iunit1,"CDELT1",cdelt,comment,status)
      call ftgkye(iunit1,"CRPIX1",cpix,comment,status)

      write(*,*)"In file: "
      write(*,'(a,f15.11)')"    CRVAL1: ",cval 
      write(*,'(a,f15.11)')"    CRPIX1: ",cpix 
      write(*,'(a,f15.11)')"    CDELT1: ",cdelt 
      write(*,*) " "

      call ftmkyd(iunit2,"CRVAL1",cval,decimals,comment,status)
      call ftmkyd(iunit2,"CDELT1",cdelt,decimals,comment,status)
      call ftmkye(iunit2,"CRPIX1",cpix,decimals,comment,status)

      call ftgkyd(iunit1,"CRVAL2",cval,comment,status)
      call ftgkyd(iunit1,"CDELT2",cdelt,comment,status)
      call ftgkye(iunit1,"CRPIX2",cpix,comment,status)

      write(*,'(a,f15.11)')"    CRVAL2: ",cval 
      write(*,'(a,f15.11)')"    CRPIX2: ",cpix 
      write(*,'(a,f15.11)')"    CDELT2: ",cdelt 

      call ftmkyd(iunit2,"crval2",cval,decimals,comment,status)
      call ftmkyd(iunit2,"cdelt2",cdelt,decimals,comment,status)
      call ftmkye(iunit2,"crpix2",cpix,decimals,comment,status)

      call FTCLOS(iunit1)
      call FTCLOS(iunit2)

      end 
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'


