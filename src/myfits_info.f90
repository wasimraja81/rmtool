 ! This subroutine reads the FITS files (images and cubes)
 ! and outputs some basic properties like bitpix, naxis,
 ! naxes, etc.
 ! This subroutine provides a handy means of extracting
 ! these properties for other programs that may require
 ! these information for passing them as inputs when calling
 ! other FITSIO subroutines.
 !
 !  -- wasim raja, 19 August, 2009




!      subroutine myfits_info(infile,
!     -           > bitpix,naxis,naxes,
!     -           cxval_im,cxpix_im,xinc_im
!     -           cyval_im,cypix_im,yinc_im
!     -           czval_im,czpix_im,zinc_im
!     -           freq_axis,cube,message,status)
!
subroutine myfits_info(infile,&
&bitpix,naxis,naxes,&
&cxval_im,cxpix_im,xinc_im,&
&cyval_im,cypix_im,yinc_im,&
&czval_im,czpix_im,zinc_im,&
&freq_axis,cube,message,status)
   implicit none
   integer*4 maxdim
   parameter(maxdim = 100)
   integer*4 bitpix, naxis, naxes(maxdim)

   integer*4 status
   logical cube

   integer*4 rwmode
   character infile*(*), message*(*)
   character*172 comment
   character*72 ctype_i
   character*8 keyname

   real*4 cxval_im, cyval_im, czval_im
   integer*4 cxpix_im, cypix_im, czpix_im
   real*4 xinc_im, yinc_im, zinc_im
   integer*4 freq_axis
   integer*4 i, tmp_status

   ! Some useless fitsio legacy stuff:
   integer*4 group, blocksize


   ! Some default values:
   group = 1
   blocksize = 1
   rwmode = 1  ! read only mode

   status = 0
   ! Open the data files:
   call FTOPEN(11,infile,rwmode,blocksize,status)
   if(status.ne.0)then
      message(1:) = "file-open error..."
      return
   endif

   ! Determine the data-type of the image (BITPIX value):
   ! Possible returned values are:
   !   8 : unsigned byte,
   !  16 : signed 2-byte integer,
   !  32 : signed 4-byte integer,
   !  64 : signed 8-byte integer,
   ! -32 : real, and
   ! -64 : double.

   call FTGIDT(11,bitpix,status)

   ! Determine the dimension(NUMBER of AXES) in
   ! the images:
   call FTGIDM(11,naxis,status)

   if(naxis.gt.99)then
      message(1:) = "NAXIS > 99, are you sure?"
      status = -1999 ! forcing status to non-zero value
      !stop
      goto 9999
   endif

   ! Determine the SIZE along each dimension of the image:
   call FTGISZ(11,maxdim,naxes,status)

   cube = .false.
   freq_axis = 0
   do i = 1,naxis
      write(keyname,'("CTYPE",I0)')i
      ctype_i = ' '
      tmp_status = 0
      call FTGKYS(11,keyname,ctype_i,comment,tmp_status)
      if(tmp_status.eq.0)then
         if(index(ctype_i,'FREQ').gt.0 .or.&
         &index(ctype_i,'freq').gt.0)then
            freq_axis = i
            cube = .true.
            exit
         endif
      endif
   enddo
   if(freq_axis.eq.0)then
      message(1:) = 'No frequency axis (CTYPE*=FREQ) found'
      status = -2002
      goto 9999
   endif

   !=======================================================

   ! Get the RA's first:
   call FTGKYE(11,"crval1",cxval_im,comment,status) ! central "value" of
   ! axis(1) or RA-centre
   call FTGKYJ(11,"crpix1",cxpix_im,comment,status) ! central "pixel" of
   ! axis(1) or RA-centre
   call FTGKYE(11,"cdelt1",xinc_im,comment,status)  ! increment in between
   ! pixels of axis(1)
   ! Get the Dec's :
   call FTGKYE(11,"crval2",cyval_im,comment,status) ! central "value" of
   ! axis(2) or Dec-centre
   call FTGKYJ(11,"crpix2",cypix_im,comment,status) ! central "pixel" of
   ! axis(2) or Dec-centre
   call FTGKYE(11,"cdelt2",yinc_im,comment,status)  ! increment in between
   ! pixels of axis(2)

   czval_im = 0.0
   czpix_im = 0
   zinc_im = 0.0
   if (cube)then
      ! Get the Spectral Channels :
      write(keyname,'("CRVAL",I0)')freq_axis
      call FTGKYE(11,keyname,czval_im,comment,status) ! reference "value" of
      ! axis(3) or Dec-centre
      write(keyname,'("CRPIX",I0)')freq_axis
      call FTGKYJ(11,keyname,czpix_im,comment,status) ! reference "pixel" of
      ! axis(3) i.e,Frequency
      write(keyname,'("CDELT",I0)')freq_axis
      call FTGKYE(11,keyname,zinc_im,comment,status)  ! increment in between
      ! pixels of axis(3)
   endif
   if(status.eq.0)then
      message(1:) = 'basic info-extraction okay!'
   endif

9999 continue
   ! CLOSE THE FITS FILE:
   call FTCLOS(11,status)
   ! -----------------------------------------------------------------

end
