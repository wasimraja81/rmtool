chelp+      
      !-----------------------------------------------------------
      ! This subroutine loads an entire FITS 2D or a sub-section
      ! of the image provided it is supplied as inputs the centre 
      ! (x,y) pixels as well as the number of pixels along the x 
      ! and the y-axes
      !
      !          -- wasim raja, 16 July, 2009
      !-----------------------------------------------------------
chelp-      
      !
      ! The current version also makes checks regarding the 
      ! location of the centre, and will prompt to the user 
      ! if the coordinate of central pixel read from the FITS
      ! file is not consistent with the total number of pixels
      ! ialong the x and the y-dimensions.
      !
      ! It is also forced here that the subimage can have only 
      ! EVEN number of pixels along both axes. So if a user 
      ! inputs odd number of output pixels, the code automati-
      ! cally adjusts it to an even number, prompting about 
      ! the action to the user. This is just to avoid extra 
      ! lines of coding without any serious effect on the 
      ! output image; after all how does an extra pixel in the
      ! output subimage affect adversely to the user?
      ! However if the "subimage" specified exactly matches 
      ! the input image size, no pixels are missed. 
      !
      ! Finally, care also has been taken to re-adjust the input 
      ! values, in case the user inputs values that exceeds the 
      ! image size. 
      ! 
      ! The data is loaded into an array called subim()  
      !
      !     -- wasim raja, 18 Aug, 2009
      !------------------------------------------------------

      subroutine load_fits_image(infile,cxpix,cypix,nxpix,nypix,
     -                           subim,xdim,ydim,status)
      implicit none

      integer*4 xdim, ydim, maxdim, maxkeys, max_inbuff
      parameter(maxdim=99,maxkeys=500, max_inbuff=16000)
      ! max_inbuff is the buffer size for the arrays meant 
      ! for reading the image/cubes. For images with a large 
      ! number of pixels, one may run out of memory if one 
      ! attempts to read the entire image into a single 
      ! array. We thus read and analyse the image/cube part 
      ! by part, taking care to prevent "out of memory" 
      ! situations by proper choice of max_inbuff parameter.
      integer*4 nchar
      character infile*(*) 
      integer*4 rwmode
      integer*4 status
      integer*4 lun
      integer*4 blocksize
      ! FOR FTGHPR:
      integer*4 bitpix, naxis, naxes(maxdim)
      integer*4 i !, j
      character*16 tmpstr
      character*12 keyname(maxkeys)
      character*18 keyword
      !character*1 yorn

      ! Definitions related to reading IMAGE:
      integer group
      real buff_array(max_inbuff),nullval
      real subim(xdim,ydim)
      logical anyflg

      integer*4 ypix_beg, ypix_end, nypix, cypix
      integer*4 xpix_beg, xpix_end, nxpix, cxpix
      integer*4 iy, ix, kx, ky
      integer*4 nx_totpix, ny_totpix 
      integer*4 nbuffer, firstpix, npixels


      real*4 cxval_im, xinc_im, cyval_im, yinc_im
      integer*4 cxpix_im, cypix_im
      integer*4 nx_1st, nx_2nd, ny_1st, ny_2nd
      integer*4 nxc, nyc
    
      real*4 xval(max_inbuff), yval(max_inbuff) 
      real*4 x1, xn, y1, yn

      character*72 comment
      integer*4 data_precision

      ! IT IS IMPORTANT TO INITIALISE "STATUS" TO 
      ! ZERO, FOR ONE TO MAKE SUBROUTINE CALLS 
      ! SUCCESSFUL... IF PASSED WITH A POSITIVE  
      ! VALUE, THE SUBROUTINE WILL EXIT IMMEDIATELY. 
      ! FITSIO ALTERS THE STATUS's VALUE IF IT 
      ! ENCOUNTERS SOME SERIOUS ERRORS DURING THE 
      ! RUN OF THE PROGRAM, SO THAT ANY SUBSEQUENT 
      ! CALLS TO SUBROUTINES ARE NOT EXECUTED...
      ! THUS "STATUS" CAN BE USED TO MONITOR ERRORS
      ! AND MESSAGES RELATED TO THEM...
      status = 0

!-------------------------------------------------------------------
      ! Get the file containing the Image/Cube:
      !


      rwmode = 0 ! read only mode
      blocksize = 1

      !infile = '../DATA/'//infile(1:nchar(infile))
      infile = infile(1:nchar(infile))
      ! Open the Image/Cube Fits file:
      call get_lun(lun) 
      call FTOPEN(lun,infile,rwmode,blocksize,status)
      if(status.ne.0)then
              write(*,*)" "
              write(*,*)"status = ", status
              write(*,*)"Error opening FITS file: ",
     -                       infile(1:nchar(infile))
              return  
              !stop
!      else
!              write(*,*)"infile:",infile(1:nchar(infile))
      endif

      ! Determine the data-type of the image (BITPIX value):
      ! Possible returned values are: 
      !   8 : unsigned byte, 
      !  16 : signed 2-byte integer, 
      !  32 : signed 4-byte integer, 
      !  64 : signed 8-byte integer, 
      ! -32 : real, and 
      ! -64 : double. 

      call FTGIDT(lun,bitpix,status)

      ! Determine the dimension(NUMBER of AXES) in 
      ! the image:
      call FTGIDM(lun,naxis,status)

      if(naxis.gt.99)then
              write(*,*)"------------- WARNING --------------"
              write(*,*)" "
              write(*,*)"NAXIS is > 99"
              write(*,*)"Currently this code expects that the "
              write(*,*)"image dimensions does not exceed 99"
              write(*,*)" "
              write(*,*)"However if you are sure that your im-"
              write(*,*)"age has NAXIS > 99, then modify the "
              write(*,*)"source code and RECOMPILE it..."
              write(*,*)" "
              write(*,*)"To modify the code: "
              write(*,*)"Open the code in your favourite text-"
              write(*,*)"editor, and search for the section: "
              write(*,*)"        'ERROR_NAXIS99' "
              write(*,*)"The section will describe to you the "
              write(*,*)"necessary modifications that needs to"
              write(*,*)"be made..."
              write(*,*)" "
              write(*,*)"-----------------------------------"
              stop
      endif
      ! Determine the KEYWORD character "values":
      ! In this case, we wish to know the "names" 
      ! of the axes...
      tmpstr(1:) = '              '
      do i = 1,naxis
      !
! ERROR_NAXIS99:
      ! The following if-section takes care of upto 
      ! 99 naxis. I do not think any FITS file will have 
      ! naxis > 99 ever.
      ! If at all such an apparently crazy need arises, 
      ! append appropriate number of "else if" sections 
      ! and corresponding "formats" defined as in label 
      ! number 499 and 599.
         if (i.lt.10)then
                 write(tmpstr(1:),fmt=499)i
                 keyword = 'ctype'//tmpstr(1:nchar(tmpstr))
         else if(i.ge.10.and.i.lt.100)then
                 write(tmpstr(1:),fmt=599)i
                 keyword = 'ctype'//tmpstr(1:nchar(tmpstr))
         endif
         call FTGKEY(lun,keyword, keyname(i),comment,status)
      enddo
499   format (i1)
599   format (i2)

      ! Determine the SIZE along each dimension of the image:
      call FTGISZ(lun,maxdim,naxes,status)

      !---------------------------------------------------
      ! Check if image is NOT a 2D image: 
      do i = 3,naxis
         if(naxes(i) .gt. 1)then
                 write(*,*)"Image is of dimension higher than 2D"
                 write(*,*)"Quitting now..."
                 stop
         endif
      enddo
      !---------------------------------------------------
      npixels = naxes(1)*naxes(2)
      nx_totpix = naxes(1)
      ny_totpix = naxes(2)
      !------------------------------------------------------------
      if(bitpix.eq.8)data_precision = 1
      if(bitpix.eq.16)data_precision = 2
      if(bitpix.eq.32)data_precision = 4
      if(bitpix.eq.64)data_precision = 8
      if(bitpix.eq.-32)data_precision = 4
      if(bitpix.eq.-64)data_precision = 8

      group = 1
      firstpix = 1
      !nullval = -999.0
      nullval = 0.0/nullval ! NaN
      nbuffer = naxes(1)
      !write(*,*)"nbuffer derived from naxes(1) = ",nbuffer


      ! Now generate the axis values:

      ! Please be careful with the units of cxval_im, xinc_im etc. 
      ! Here it is assumed that these are in degrees.
      !
      ! Get the RA's first:
      call FTGKYE(lun,"crval1",cxval_im,comment,status) ! central "value" of 
                                                       ! axis(1) or RA-centre
      call FTGKYJ(lun,"crpix1",cxpix_im,comment,status) ! central "pixel" of 
                                                       ! axis(1) or RA-centre
      call FTGKYE(lun,"cdelt1",xinc_im,comment,status)  ! increment in between 
                                                       ! pixels of axis(1)
      ! Now get the Dec's :
      call FTGKYE(lun,"crval2",cyval_im,comment,status) ! central "value" of 
                                                       ! axis(2) or Dec-centre
      call FTGKYJ(lun,"crpix2",cypix_im,comment,status) ! central "pixel" of 
                                                       ! axis(2) or Dec-centre
      call FTGKYE(lun,"cdelt2",yinc_im,comment,status)  ! increment in between 
                                                       ! pixels of axis(2)
  
      ! Check if the reference pixel is indeed at the centre of the 
      ! image array and also find out the number of points leading 
      ! and lagging the reference pixel:
  
      ! If total number of pixel is EVEN, centre must lie at n_totpix/2 
      ! or n_totpix/2 + 1
      if(mod(nx_totpix,2) .eq. 0)then   
              nxc = nx_totpix/2
              if(cxpix_im .eq. nxc)then
                      nx_1st = nxc - 1
                      nx_2nd = nxc 
              else if(cxpix_im .eq. nxc + 1)then
                      nx_1st = nxc
                      nx_2nd = nxc - 1
              else
                      nx_1st = cxpix_im - 1
                      nx_2nd = nx_totpix - cxpix_im
              endif
      ! If total number of pixel is ODD, centre must lie at (n_totpix + 1)/2
      elseif(mod(nx_totpix,2) .eq. 1)then
              nxc = (nx_totpix+1)/2
              if(cxpix_im .eq. nxc)then
                      nx_1st = nxc - 1
                      nx_2nd = nxc - 1
              else
                      nx_1st = cxpix_im - 1
                      nx_2nd = nx_totpix - cxpix_im
              endif
      endif
  
  
      ! If total number of pixel is EVEN, centre must lie at n_totpix/2 
      ! or n_totpix/2 + 1
      if(mod(ny_totpix,2) .eq. 0)then   
              nyc = ny_totpix/2
              if(cypix_im .eq. nyc)then
                      ny_1st = nyc - 1
                      ny_2nd = nyc 
              else if(cypix_im .eq. nyc + 1)then
                      ny_1st = nyc
                      ny_2nd = nyc - 1
              else
                      ny_1st = cypix_im - 1
                      ny_2nd = ny_totpix - cypix_im
              endif
      ! If total number of pixel is ODD, centre must lie at (n_totpix + 1)/2
      elseif(mod(ny_totpix,2) .eq. 1)then
              nyc = (ny_totpix+1)/2
              if(cypix_im .eq. nyc)then
                      ny_1st = nyc - 1
                      ny_2nd = nyc - 1
              else
                      ny_1st = cypix_im - 1
                      ny_2nd = ny_totpix - cypix_im
              endif
      endif
  
      !Refreshing Standard IX Arithmetic-Series concepts: 
      !val_1 = val_c - n_1st*delta_val
      !val_ntot = val_1 + (ntot - 1)*delta_val
  
      ! generating the x-axis values...
      x1 = cxval_im - nx_1st*xinc_im
      xn = x1 + (nx_totpix - 1)*xinc_im
      call linspace(x1,xn,nx_totpix,xval)
  
      ! generating the y-axis values...
      y1 = cyval_im - ny_1st*yinc_im
      yn = y1 + (ny_totpix - 1)*yinc_im
      call linspace(y1,yn,ny_totpix,yval)

       
!      write(*,*)"========================================"
!      write(*,*)"Before Readjustment: "
!      write(*,*)"   "
!      write(*,*)"     cxpix,nxpix: ",cxpix,nxpix
!      write(*,*)"     cypix,nypix: ",cypix,nypix
!      write(*,*)"   "
!      write(*,*)"========================================"
      ! We will now define the subimage of (nxpix x nypix) 
      ! pixels centered around the pixel: (cxpix,cypix)

      if(nxpix .lt. nx_totpix .and. nxpix .gt. 0)then
         if(mod(nxpix,2) .ne. 0)then
              !   write(*,*)"nxpix expected even, found odd!"
                 nxpix = nxpix - 1
              !   write(*,*)"nxpix further re-adjusted to: ",nxpix
         endif
      endif  
      if(nypix .lt. ny_totpix .and. nypix .gt. 0)then
         if(mod(nypix,2) .ne. 0)then
              !   write(*,*)"nypix expected even, found odd!"
                 nypix = nypix - 1 
              !   write(*,*)"nypix further re-adjusted to: ",nypix
         endif
      endif  
   
      if(nxpix .gt. nx_totpix.or.nxpix.le.0)then
              nxpix = nx_totpix 
              if(mod(nxpix,2) .ne. 0)then
                     nxpix = nxpix - 1 
              endif
      endif 
      if(cxpix .lt. nxpix/2 + 1)then
           cxpix = nxpix/2 + 1 
      endif 
      if(cxpix .gt. nx_totpix + 1 - nxpix/2)then
           cxpix = nx_totpix + 1 - nxpix/2 
      endif

      if(nypix .gt. ny_totpix.or.nypix.le.0)then
           nypix = ny_totpix 
           if(mod(nypix,2) .ne. 0)then
                   nypix = nypix - 1 
           endif
      endif
      if(cypix .lt. nypix/2 + 1)then
           cypix = nypix/2 + 1 
      endif
      if(cypix .gt. ny_totpix + 1 - nypix/2)then
           cypix = ny_totpix + 1 - nypix/2 
      endif

      ypix_beg = cypix-nypix/2 
      ypix_end = cypix + nypix/2 - 1
      if(ypix_end.eq.ny_totpix)then
              ypix_end = ypix_end-1
      endif

      xpix_beg = cxpix-nxpix/2
      xpix_end = cxpix + nxpix/2 - 1

!      write(*,*)"========================================"
!      write(*,*)"After Readjustment: "
!      write(*,*)"   "
!      write(*,*)"     cxpix,nxpix: ",cxpix,nxpix
!      write(*,*)"     cypix,nypix: ",cypix,nypix
!      write(*,*)"   "
!      write(*,*)"xpixbeg, xpixend: ",xpix_beg, xpix_end 
!      write(*,*)"ypixbeg, ypixend: ",ypix_beg, ypix_end 
!      write(*,*)"========================================"



      if(bitpix.eq.8)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVB(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> kx
                                             ! col --> ky
            enddo
            firstpix = firstpix + nbuffer
         enddo
         close(31)
      else if(bitpix.eq.16)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVI(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> ix
                                             ! col --> iy
            enddo
            firstpix = firstpix + nbuffer
         enddo
         close(31)
      else if(bitpix.eq.32)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVJ(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> ix
                                             ! col --> iy
            enddo
            firstpix = firstpix + nbuffer
         enddo
         close(31)
      else if(bitpix.eq.64)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVK(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> ix
                                             ! col --> iy
            enddo
            firstpix = firstpix + nbuffer
         enddo
         close(31)
      else if(bitpix.eq.-32)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVE(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end 
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> ix
                                             ! col --> iy
            enddo
            firstpix = firstpix + nbuffer
         enddo
         close(31)
      else if(bitpix.eq.-64)then
         ky = 0
         firstpix = ypix_beg*nbuffer 
         do iy = ypix_beg,ypix_end
            call FTGPVD(lun,group,firstpix,nbuffer,nullval,
     -                   buff_array,anyflg,status)
            if(status.ne.0)then
                    call printerror(status)
                    stop
            endif
            ky = ky + 1
            kx = 0
            do ix = xpix_beg,xpix_end
               kx = kx + 1
               subim(kx,ky) = buff_array(ix) ! row --> ix
                                             ! col --> iy
            enddo
            firstpix = firstpix + nbuffer
         enddo
      endif

      ! CLOSE THE FITS FILE:
      call FTCLOS(lun,status)
      if (status .gt. 0) then 
              call printerror(status)
              stop
      endif


      end


