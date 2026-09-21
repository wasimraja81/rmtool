chelp+
      !-------------------------------------------------------------
      ! This code extracts a set of planes from an image cube. 
      ! The 3rd axis is considered to be the axis perpendicular 
      ! to the image plane. 
      !                                    -- wr, 19 Aug, 2009
      !-------------------------------------------------------------
chelp-



      implicit none
      include '../INCLUDE/myfits_spec2rm.inc'

      
      real*4    specI(max_dec*maxchan) 
      integer*4 bitpix, naxis, naxes(max_axis), naxes_out(max_axis)
      logical simple, extend
      integer*4 decimals

      real*4 cxval_im, cyval_im, czval_im
      integer*4 cxpix_im, cypix_im, czpix_im
      real*4 xinc_im, yinc_im, zinc_im
      integer*4   ixpix_now, ix, nx_use, ny_use, nz_use,
     -            xpix_beg, xpix_end, ypix_beg, ypix_end, 
     -            zpix_beg, zpix_end, ntot_use 

      real*4    z_extract 

      integer*4 nx_totpix, ny_totpix, nz_totpix 
      integer*4 nbuffer, firstpix

      integer*4 fpixels(max_axis), lpixels(max_axis), incs(max_axis)
      character*8 junkchar
      integer*4 status, nchar
      logical anyflg

      character*72 message 
      real*8 pi 

      integer*4 rwmode
      character*272 infile 
      character*272 outfile 
      character*272 cfgfile
      character*172 path  

      integer*4 nx_1st, nx_2nd, ny_1st, ny_2nd, nz_1st, nz_2nd
      integer*4 nxc, nyc, nzc
      real*4    cdelt, cval 
      integer*4 cpix 
      integer*4  npix_use 
      real*4     zval_use 
    

      integer*4 data_precision
      real*4 nullval
      ! various counters and indices:
      integer*4 i 


      ! Some useless fitsio legacy stuff:
      integer*4 group, blocksize

      ! temporary variables: 
      logical cube 


      pi = acos(-1.0d0)
!-------------------------------------------------------------------
      ! SANITY CHECKS:
      ! Compare the files containing the Q and U Cubes
      ! ans see if they are compatible with each other:

      if(iargc().lt.1)then
              write(*,*)'  '
              write(*,*)' Usage: '
              write(*,*)'> extract_image_from_cube <cfgfile> '
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)' You need a config file containing the '
              write(*,*)' the parameters for this run. '
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)'  '
              stop
      else
              call getarg(1,cfgfile)
              cfgfile = cfgfile(1:nchar(cfgfile))
      endif



      cfgfile = '../CONFIG/'//cfgfile(1:nchar(cfgfile))
      open(11,file=cfgfile,status='old',err=101)
      goto 102

101   write(*,*)"Error opening config file: ",cfgfile(1:nchar(cfgfile))
      write(*,*)"Quitting now..."
      write(*,*)" "
      stop

102   continue


      read(11,*)junkchar     ! comment line
      read(11,'(a)')path
      path = path(1:index(path,';')-1)
      path = path(1:nchar(path))

      read(11,'(a)')infile
      infile = infile(1:index(infile,';')-1)
      infile = infile(1:nchar(infile))
      read(11,'(a)')outfile
      outfile = outfile(1:index(outfile,';')-1)
      outfile = outfile(1:nchar(outfile))
      read(11,*)z_extract  ! The z value to extract 

      close(11)


      ! Do not write the additional files if the 
      ! entire cube is being processed:
      infile(1:) = path(1:nchar(path))//infile(1:nchar(infile))
      write(*,*)"Opening file: ",infile(1:nchar(infile))

      outfile(1:) = outfile(1:nchar(outfile))//'.FITS'

      ! Extract Some basic INFO from the FITS files:
      status = 0 
      call myfits_info(infile,
     -           bitpix,naxis,naxes,
     -           cxval_im,cxpix_im,xinc_im,
     -           cyval_im,cypix_im,yinc_im,
     -           czval_im,czpix_im,zinc_im,
     -           cube,message,status)

      if (status.eq.0)then
              write(*,*)"In-cube opened:",infile(1:nchar(infile))
              write(*,*)"      bitpixQ:",bitpix
              write(*,*)"       naxisQ:",naxis
              write(*,*)" "
              write(*,*)"   ref. x-val:",cxval_im
              write(*,*)"   ref. x-pix:",cxpix_im
              write(*,*)"         xinc:",xinc_im
              write(*,*)" "
              write(*,*)"   ref. y-val:",cyval_im
              write(*,*)"   ref. y-pix:",cypix_im
              write(*,*)"         yinc:",yinc_im
              write(*,*)" "
              write(*,*)"   ref. z-val:",czval_im
              write(*,*)"   ref. z-pix:",czpix_im
              write(*,*)"         zinc:",zinc_im
              write(*,*)" "
              write(*,*)"        cube:",cube
              write(*,*)"      message:",message(1:nchar(message))
              do i = 1,naxis
                 write(*,*)"naxes(",i,") = ",naxes(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the In-cube file as infile"
              write(*,*)"Check if the file exists..."
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif



       nx_totpix = naxes(1)
       ny_totpix = naxes(2)
       nz_totpix = naxes(3)

       if(mod(nx_totpix,2) .eq. 0)then   
               nxc = nx_totpix/2
               if(cxpix_im .eq. nxc)then
                       nx_1st = nxc - 1
                       nx_2nd = nxc 
               else if(cxpix_im .eq. nxc + 1)then
                       nx_1st = nxc
                       nx_2nd = nxc - 1
               else
                       if(cxpix_im.eq.0)then
                               nx_1st = 0
                               nx_2nd = nx_totpix - 1
                       else
                               nx_1st = cxpix_im - 1
                               nx_2nd = nx_totpix - cxpix_im
                       endif
               endif
       ! If total number of pixel is ODD, centre must lie at (n_totpix + 1)/2
       elseif(mod(nx_totpix,2) .eq. 1)then
               nxc = (nx_totpix+1)/2
               if(cxpix_im .eq. nxc)then
                       nx_1st = nxc - 1
                       nx_2nd = nxc - 1
               else
                       if(cxpix_im.eq.0)then
                               nx_1st = 0
                               nx_2nd = nx_totpix - 1
                       else
                               nx_1st = cxpix_im - 1
                               nx_2nd = nx_totpix - cxpix_im
                       endif
               endif
       endif
  
  
       ! For the y-axis
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
                       if(cypix_im.eq.0)then
                               ny_1st = 0
                               ny_2nd = ny_totpix - 1
                       else
                               ny_1st = cypix_im - 1
                               ny_2nd = ny_totpix - cypix_im
                       endif
               endif
       ! If total number of pixel is ODD, centre must lie at (n_totpix + 1)/2
       elseif(mod(ny_totpix,2) .eq. 1)then
               nyc = (ny_totpix+1)/2
               if(cypix_im .eq. nyc)then
                       ny_1st = nyc - 1
                       ny_2nd = nyc - 1
               else
                       if(cypix_im.eq.0)then
                               ny_1st = 0
                               ny_2nd = ny_totpix - 1
                       else
                               ny_1st = cypix_im - 1
                               ny_2nd = ny_totpix - cypix_im
                       endif
               endif
       endif
  
 
       ! If total number of pixel is EVEN, centre must lie at n_totpix/2 
       ! or n_totpix/2 + 1
       if(mod(nz_totpix,2) .eq. 0)then
               nzc = nz_totpix/2
               if(czpix_im .eq. nzc)then
                       nz_1st = nzc - 1
                       nz_2nd = nzc 
               else if(czpix_im .eq. nzc + 1)then
                       nz_1st = nzc
                       nz_2nd = nzc - 1
               else
                       if(czpix_im.eq.0)then
                               nz_1st = 0
                               nz_2nd = nz_totpix - 1
                       else
                               nz_1st = czpix_im - 1
                               nz_2nd = nz_totpix - czpix_im
                       endif
               endif
       ! If total number of pixel is ODD, centre must lie at (n_totpix + 1)/2
       elseif(mod(nz_totpix,2) .eq. 1)then
               nzc = (nz_totpix+1)/2
               if(czpix_im .eq. nzc)then
                       nz_1st = nzc - 1
                       nz_2nd = nzc - 1
               else
                       if(czpix_im.eq.0)then ! taking care of C vs. Fortran programmers
                               nz_1st = 0
                               nz_2nd = nz_totpix - 1
                       else
                               nz_1st = czpix_im - 1
                               nz_2nd = nz_totpix - czpix_im
                       endif
               endif
       endif
  
       write(*,*)" "
       write(*,*)"Sanity checks performed successfully..."
       write(*,*)" "
      ! End of sanity checks...
      !=======================================================


      if(bitpix.eq.8)data_precision = 1
      if(bitpix.eq.16)data_precision = 2
      if(bitpix.eq.32)data_precision = 4
      if(bitpix.eq.64)data_precision = 8
      if(bitpix.eq.-32)data_precision = 4
      if(bitpix.eq.-64)data_precision = 8




      !=======================================================
      group = 1
      firstpix = 1
      nullval = -999.0
      nbuffer = naxes(1)
      rwmode = 1
      blocksize = 1
      ! Open the Image/Cube Fits file:

      ! Initialise STATUS to zero:
      status = 0

      call FTOPEN(21,infile,rwmode,blocksize,status)

      if(status.ne.0)then
              write(*,*)" "
              write(*,*)"In-infile chosen:",infile(1:nchar(infile))
              write(*,*)"status = ", status
              write(*,*)"Error opening In-FITS file..."
              stop
      else
              write(*,*)" "
              write(*,*)"In-file chosen:",infile(1:nchar(infile))
      endif

      !  Create the new RM FITS files. The blocksize parameter is a
      !  historical artifact and the value is ignored by FITSIO.

      call FTINIT(41,outfile,blocksize,status)


      !=======================================================
      ! Main task of the program begins here...

      ! Decide whether the entire cubes need to be read or a
      ! part of them...

      do i = 1,naxis
         fpixels(i) = 1
         lpixels(i) = naxes(i)
         incs(i) = 1
      enddo

      ! Which plane does z_extract fall on: 
      npix_use = czpix_im + int((z_extract - czval_im)/zinc_im)
      zval_use = czval_im + (npix_use - czpix_im)*zinc_im 
      write(*,*)"Requested out plane with z: ",z_extract
      write(*,*)"Found nearest plane with z: ",zval_use, ", 
     -nz: ",npix_use
  
      extend= .false.
      simple = .true.
      ! Copy the entire header from one of the inputs 
      ! to the output file: 
      call FTCPHD(21,41,status) 

      ! Modify the appropriate headers for output: 
      ! Modify the values for the 3rd axis in the header: 
      cdelt = zinc_im  
      cpix = 1 
      cval = zval_use 
      decimals = 13
      ! 1. naxis3: 
      call ftmkyj(41,"naxis3",cpix,"Npix along 3rd axis ",status)
      ! 2. crval3
      call ftmkyd(41,"crval3",cval,decimals,
     -            "value of reference pix",status)
      ! 3. crpix3
      call ftmkyj(41,"crpix3",cpix,"location of reference pix",status)
      ! 4. cdelt3
      call ftmkyd(41,"cdelt3",cdelt,decimals,
     -             "increment along 3rd dim",status)

      !call FTPHPR(41,simple,bitpix,3,naxes_out,0,1,extend,status)

!  Put (append) a new keyword of the appropriate datatype into the CHU
!  subroutine ftpky[e,d,f,g](ounit,keywrd,rval,decim,comm,status)
!                              OR
!  subroutine ftpky[j,k,l,s](ounit,keywrd,keyval,comment,status)
!
!*******************************************************************************
!
!! FTPKYE writes a real*4 value to a header record in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)

      ! Problems related to axis mismatch noticed in the RM 
      ! planes along the direction of RA. 
      ! We will keep the RA/Dec reference pixels as it was 
      ! in the original frequency data cube, and shall only 
      ! change the "cdelt" values depending on the "incs" 
      ! input by user: 

      ! CLOSE THE FITS FILES and open them afresh!
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing In-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS In-cube..."
      endif

      ! Open the I-cube as well for calibration purpose : wr, 18 Apr, 2012
      call FTOPEN(21,infile,rwmode,blocksize,status)


      xpix_beg = fpixels(1)
      xpix_end = lpixels(1)
      nx_use = int((xpix_end - xpix_beg)/incs(1)) + 1

      ypix_beg = fpixels(2)
      ypix_end = lpixels(2)
      ny_use = int((ypix_end - ypix_beg)/incs(2)) + 1

      zpix_beg = fpixels(3)
      zpix_end = lpixels(3)
      nz_use = int((zpix_end - zpix_beg)/incs(3)) + 1

      ntot_use = nx_use*ny_use*nz_use

      ixpix_now = 0 
      do ix = xpix_beg,xpix_end,incs(1) 
         ixpix_now = ixpix_now + 1 
         
         fpixels(1) = ix 
         lpixels(1) = ix 

         fpixels(2) = ypix_beg 
         lpixels(2) = ypix_end 

         ! Get the image plane you wish: 
         fpixels(3) = npix_use 
         lpixels(3) = npix_use 

         call FTGSVE(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                 nullval,specI,anyflg,status)


         fpixels(1) = ixpix_now 
         lpixels(1) = ixpix_now 
         fpixels(2) = 1 
         lpixels(2) = naxes(2)  
         fpixels(3) = 1 
         lpixels(3) = 1 
         naxes_out(1) = naxes(1) 
         naxes_out(2) = naxes(2) 
         naxes_out(3) = 1  

         call ftpsse(41,group,naxis,naxes_out,fpixels,lpixels,
     -               specI,status)
         if(status.gt.0)then
                 call printerror(status) 
         endif
      enddo


      call FTCLOS(41,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing out-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS outfile..."
      endif
      !=======================================================


      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing In-fitsfile!!"
              call printerror(status)
      endif

      ! -----------------------------------------------------------------


      end

      include '/usr/lib/subroutine_lib/nchar.f'
      include 'myfits_info.f'
      !include 'extract_general.f'
      !include '/usr/lib/subroutine_lib/extract_general_v3.f'
      include '/usr/lib/subroutine_lib/extract_general_v4.f'
      include '/usr/lib/subroutine_lib/extract_general_setup.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
