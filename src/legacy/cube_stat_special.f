chelp+
      ! Modified to compute the residual of RM-profiles of "Stokes-I" 
      ! from the Dty RM-beam at RM = 0. The RM profiles of Stokes-I 
      ! would replicate the profile of RM = 0. In addition, it will 
      ! also carry the NOISE information, which we desperately desire 
      ! for assessing the significance of RM-components in the RM-cube. 
      ! 
      ! Some facts: 
      ! 1. I_rms = Q_rms = U_rms (for unpol sources). 
      ! 2. Bandshape modulations make it difficult to assess noise in 
      !    I from its spectra. 
      ! 3. Q and U cannot be used for determining the noise, since 
      !    their spectra may already have modulations from polarized 
      !    sky components. 
      ! 4. But the dc in Stokes-I being very high, its RM-spectra 
      !    at the high RM ends will have contributions mostly from 
      !    the dc, even if small residual modulation in band-shape 
      !    remains. 
      ! 5. Subtracting the RM-beam scaled by the dc value of Stokes-I 
      !    from the RM profile of Stokes-I, will allow one to measure 
      !    the rms accurately at the far (high RMs) of the residual. 
      !
      !    USE THIS CODE ONLY FOR COMPUTING THE RESIDUAL AND ITS 
      !    PROPERTIES FOR "RM-CUBES" FOR STOKES "I-ONLY" 
      !    
      !    the SUBIM parameters are ignored for everything else but 
      !    for the computation of the RMS of the edge RM values. 
      ! 
      !                            --wr, 28 July, 2013 
      !-------------------------------------------------------------
      ! This code computes various statistical parameters of data 
      ! in a FITS image cube. It was designed to analyse RM-cubes 
      ! as well as Spectral Data Cubes. 
      ! 
      ! The reading of data is moderately fast. The entire "2-3" 
      ! plane is read at one go for a given pixel along the 1st 
      ! axis. 
      ! 
      ! Care must be taken to interpret the axes. 
      !           -- The code assumes that: 
      !                        naxis(1) = RA
      !                        naxis(2) = Dec
      !                        naxis(3) = RM/Freq, in the input 
      !                                   cubes.
      !
      !                                    -- wr, 14 Nov, 2011
      !
      ! Provision to remove any slope in the 3rd dimension viz. 
      ! linear trend along the frequency axis incorporated. 
      ! The rms across the 3rd dimension should not be affected 
      ! by such trends, hence it may be desirable to remove 
      ! the slope to estimate the rms correctly. 
      !                                    -- wr, 10 Apr, 2013
      !-------------------------------------------------------------
chelp-



      implicit none
      !include '../INCLUDE/myfits_spec2rm.inc'
      integer*4 max_axis               ! maximum number of axis in FITS file
      parameter(max_axis = 100)
      integer*4 max_ra, max_dec, maxchan
      parameter(max_ra = 1024, max_dec = 1024, maxchan=1024)
      integer*4 max_pix               ! maximum number of pixels allowed
      parameter(max_pix = 134217728)  ! equivalent of 512 MB in real*4


      
      real*4    data_arr(maxchan), specRM(max_dec*maxchan),
     -          dty_beam(maxchan), beam_arr_good(maxchan)
      integer*4 bitpix, naxis, naxes(max_axis), naxes_out(max_axis)
      logical   simple, extend
      integer*4 decimals

      real*4 cxval_im, cyval_im, czval_im
      integer*4 cxpix_im, cypix_im, czpix_im
      real*4 xinc_im, yinc_im, zinc_im


      integer*4 xpix_beg, xpix_end
      integer*4 ypix_beg, ypix_end
      integer*4 zpix_beg, zpix_end
      integer*4 npts_beam 


      integer*4 nx_totpix, ny_totpix, nz_totpix 
      integer*4 nx_use, ny_use, nz_use, ntot_use
      integer*4 nbuffer, firstpix

      integer*4 fpixels(max_axis), lpixels(max_axis), incs(max_axis)
      real*4    data_good(maxchan) 
      character*8 junkchar
      integer*4 status, nchar
      logical anyflg
      logical cube

      character*64 ctype 
      character*72 comment
      real*4 cval,cdelt,cpix 

      integer*4 rwmode
      character*272 infile, message, dty_beam_file 
      character*272 outfile, outfileMAXLOC,outfileMAXVAL, outfileTHRESH
      character*272 outfileTOTFLUX,outfileRMS, outfileMEAN
      character*272 subim_parfile, cfgfile
      character*172 path
      character*1 yorn
      real*4      xarr(maxchan), zarr(maxchan), A, B 
      logical     remove_spectral_slope 

      integer*4 nx_1st, nx_2nd, ny_1st, ny_2nd, nz_1st, nz_2nd
      integer*4 nxc, nyc, nzc
    
      real*4 xval(max_ra), yval(max_dec), zval(maxchan),
     -       zval_good(maxchan) 
      real*4 x1, xn, y1, yn, z1, zn

      integer*4 data_precision
      real*4 nullval
      logical subim
      ! various counters and indices:
      integer*4 i, kk, ix, iy, ixpix_now, iypix_now, imax 
      integer*4 icnt, cnt1, tmp_cnt1, tmp_cnt2, tmp_index 


      real*4 maxloc_arr(max_dec),maxval_arr(max_dec),thresh_arr(max_dec)
      real*4 totflux_arr(max_dec), rms_arr(max_dec), mean_arr(max_dec) 
      real*4 maxxval, mu, std 

      ! RFI related (list of bad-channels based on apriori info)
      real*4 bad_chan(maxchan)
      integer*4 nbad_chan, ngood_chan, flag_arr(maxchan), nuse 
      logical remove_badchan
      character*172  badchan_file

      ! processing related:
      integer*4 nz_out 

      ! Some useless fitsio legacy stuff:
      integer*4 group, blocksize

      ! temporary variables: 
      real*4 atmp  


!-------------------------------------------------------------------

      if(iargc().lt.1)then
              write(*,*)'  '
              write(*,*)' Usage: '
              write(*,*)'> cube_stat <cfgfile> '
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)' You need a config file containing the '
              write(*,*)' the parameters for this run. '
              write(*,*)'  '
!              write(*,*)' You can make some additional requests: '
!              write(*,*)' using this string. Valid requests as '
!              write(*,*)' of now are: '
!              write(*,*)' 1) single_cut: to be used when you intend'
!              write(*,*)'                to write out the Q,U and '
!              write(*,*)'                RM-spectra for only a single'
!              write(*,*)'                "cut" in the sky. By "cut"'
!              write(*,*)'                I mean all pixels for eg.,'
!              write(*,*)'                having constant Dec value.'
!              write(*,*)'   NB: The subim_parfile must be appropriately'
!              write(*,*)'       written for this.'
!              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)'  '
              stop
      else if(iargc().eq.1)then
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
      read(11,*)yorn
      if(yorn.eq.'y'.or.yorn.eq.'Y')then
              remove_badchan = .true.

              read(11,'(a)')badchan_file
              badchan_file = badchan_file(1:index(badchan_file,';')-1)
              badchan_file = badchan_file(1:nchar(badchan_file))
      else
              remove_badchan = .false.
              read(11,*)junkchar
      endif
      read(11,*)yorn
      if(yorn.eq.'y'.or.yorn.eq.'Y')then
              subim = .true.

              read(11,'(a)')subim_parfile
              subim_parfile=subim_parfile(1:index(subim_parfile,';')-1)
              subim_parfile=subim_parfile(1:nchar(subim_parfile))
      else
              subim = .false.
              read(11,*)junkchar
      endif

      read(11,*)yorn
      if(yorn.eq.'y'.or.yorn.eq.'Y')then
              remove_spectral_slope = .true.
      else
              remove_spectral_slope = .false.
      endif


      close(11)



      infile(1:) = path(1:nchar(path))//infile(1:nchar(infile))

      outfileMAXLOC(1:) = outfile(1:nchar(outfile))//'.MAXLOC.FITS'
      outfileTHRESH(1:) = outfile(1:nchar(outfile))//'.THRESH.FITS'
      outfileMAXVAL(1:) = outfile(1:nchar(outfile))//'.MAXVAL.FITS'
      outfileTOTFLUX(1:) = outfile(1:nchar(outfile))//'.TOTFLUX.FITS'
      outfileRMS(1:) = outfile(1:nchar(outfile))//'.RMS.FITS'
      outfileMEAN(1:) = outfile(1:nchar(outfile))//'.MEAN.FITS'

      badchan_file(1:) = '../PAR/'//badchan_file(1:nchar(badchan_file))

      ! Read the Dirty Beam file: 
      dty_beam_file = '../PAR/dty_rm_beam_for_casa_at_rm_zero.txt'
      open(61,file=dty_beam_file,status='old',err=801)
      goto 802
801   continue
      write(*,*)" "
      write(*,*)"Error opening file: ",
     -          dty_beam_file(1:nchar(dty_beam_file))
      write(*,*)"Please ensure the file is in ../PAR/ area!!"
      write(*,*)" "
      stop 
802   continue
      npts_beam = 0 
      do while(.true.)
        npts_beam = npts_beam + 1
        read(61,*,end=803)atmp,dty_beam(npts_beam)   ! Reading the BAD-CHANNEL NUMBERS
      enddo
803   continue
      npts_beam = npts_beam - 1
      write(*,*)"Number of points in beam file: ",npts_beam
      close(61)
804   continue


      ! Read the BAD-CHANNELS required for FLAGGING them:
      if(remove_badchan)then
              open(71,file=badchan_file,status='old',err=701)
              goto 702
701           continue
              write(*,*)" "
              write(*,*)"Error opening file: ",
     -                  badchan_file(1:nchar(badchan_file))
              write(*,*)"Please ensure the file is in ../PAR/ area!!"
              write(*,*)" "
              write(*,*)"We shall not attempt to FLAG bad-channels!"
              remove_badchan = .false.
              goto 704

702           continue
              nbad_chan = 0
              do while(.true.)
                nbad_chan = nbad_chan + 1
                read(71,*,end=703)bad_chan(nbad_chan)  ! Reading the BAD-CHANNEL NUMBERS
                write(*,*)"bad-channels: ",bad_chan(nbad_chan)
              enddo
703           continue
              nbad_chan = nbad_chan - 1
              write(*,*)"Number of Bad Channels: ",nbad_chan
              close(71)
      endif
704   continue

      ! Read the parameters required for the RM-EXTRACTION
      open(31,file='../PAR/extract_rm.par',status='old',err=601)
      goto 602
601   continue
      write(*,*)" "
      write(*,*)"Error opening file: extract_rm.par"
      write(*,*)" Ensure that you have defined "
      write(*,*)" the parameters needed for the "
      write(*,*)" extraction in a file named: "
      write(*,*)" 'extract_rm.par kept in the 'PAR/ "
      write(*,*)" directory and re-run the code."
      write(*,*)" "
      write(*,*)"Quitting now... "
      stop

602   continue

      close(31)


      ! Extract Some basic INFO from the FITS files:
      call myfits_info(infile,
     -           bitpix,naxis,naxes,
     -           cxval_im,cxpix_im,xinc_im,
     -           cyval_im,cypix_im,yinc_im,
     -           czval_im,czpix_im,zinc_im,
     -           cube,message,status)

      if (status.eq.0)then
              write(*,*)" Cube opened:",infile(1:nchar(infile))
              write(*,*)"      bitpix:",bitpix
              write(*,*)"       naxis:",naxis
              write(*,*)" "
              write(*,*)"  ref. x-val:",cxval_im
              write(*,*)"  ref. x-pix:",cxpix_im
              write(*,*)"        xinc:",xinc_im
              write(*,*)" "
              write(*,*)"  ref. y-val:",cyval_im
              write(*,*)"  ref. y-pix:",cypix_im
              write(*,*)"        yinc:",yinc_im
              write(*,*)" "
              write(*,*)"  ref. z-val:",czval_im
              write(*,*)"  ref. z-pix:",czpix_im
              write(*,*)"        zinc:",zinc_im
              write(*,*)" "
              write(*,*)"        cube:",cube
              write(*,*)"     message:",message(1:nchar(message))
              do i = 1,naxis
                 write(*,*)"naxes(",i,") = ",naxes(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the inCube file as infile"
              write(*,*)"Check if the file exists..."
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif






      ! Final sanity checks...
      !
      ! Please be careful with the units of cxval_im, xinc_im etc. 
      !
      ! NOTE: It appeared to me that it's a practice to make the central 
      !       RA and DEC pixels as the reference pixels. However, in case 
      !       of the frequency axis, it is the 1st pixel that is taken
      !       as the reference. It was later found that intermediate 
      !       AIPS/MIRIAD tasks, when used to write out FITS file, do 
      !       not religiously follow such conventions! The reference 
      !       pixel can be any pixel -- I am not aware on what
      !       determines the choice of the reference pixel though.
      !       
      !       In this code, we check whether the reference pixels of all 
      !       the axes are indeed the central pixels or not.
      !       We further check for the C-Fortran index conventions in
      !       cases where the reference pixel, if happens to be the 1st
      !       pixel, whether a 0 is assigned to it, or a 1. 
      !       C-programmers usually refer to the 1st pixel as the 
      !       0-th pixel, whereas Fortran programmers assign index 1 to 
      !       to the 1st pixel. Hence we will assume the reference pixels 
      !       to be 1 in cases where we find the reference values tagged 
      !       to pixel number 0 in the FITS file. 
      !


      write(*,*)"! -----------------------------------------------------
     ------------" 
      write(*,*)"! Final sanity checks..."
      write(*,*)"! "      
      write(*,*)"! NOTE: It appeared to me that it's a practice to make 
     -the central" 
      write(*,*)"! RA and DEC pixels as the reference pixels,  whereas, 
     -in case of "
      write(*,*)"! the frequency axis, it is the 1st pixel that is taken
     - as the "  
      write(*,*)"! reference. It was later found that intermediate AIPS/
     -MIRIAD tasks"
      write(*,*)"! when used to write out FITS files, do not religiously
     - follow such" 
      write(*,*)"! conventions! The reference pixel can be any pixel -- 
     -I am not aware"
      write(*,*)"! on what determines the choice of the reference pixel 
     -though. "
      write(*,*)"! Maybe the programmer's bias..."
      write(*,*)"! "
      write(*,*)"! In this code, we check whether the reference pixels o
     -f all the "
      write(*,*)"! axes are indeed the central pixels, and WARN the user
     - if otherwise."
      write(*,*)"! "
      write(*,*)"! We further check for the C-Fortran index conventions 
     -in"
      write(*,*)"! cases where the reference pixel, if happens to be the
     - 1st pixel,"
      write(*,*)"! whether a 0 is assigned to it, or a 1."
      write(*,*)"! C-programmers usually refer to the 1st pixel as the 0
     --th pixel,"
      write(*,*)"! whereas Fortran programmers assign index 1 to the 1st
     - pixel."
      write(*,*)"! Hence we will assume the reference pixels to be 1 in 
     -cases where"
      write(*,*)"! we find the reference values tagged to pixel number 0
     - in the "
      write(*,*)"! FITS file."
      write(*,*)"! Feel free to correct me in case I have missed somethi
     -ng:"
      write(*,*)"!               wasim@rri.res.in"
      write(*,*)"! -----------------------------------------------------
     ------------" 

        ! Check if the reference pixel is indeed at the centre of the 
        ! image array and also find out the number of points leading 
        ! and lagging the reference pixel:
 


        ! For the x-axis
        ! If total number of pixel is EVEN, centre must lie at n_totpix/2 
        ! or n_totpix/2 + 1



       nx_totpix = naxes(1)
       ny_totpix = naxes(2)
       nz_totpix = naxes(3)
       if(nz_totpix.ne.npts_beam)then
               write(*,*)"!!!!!!!!!!!!!   Error   !!!!!!!!!!!!!"
               write(*,*)"  "
               write(*,*)"       npts in beam file: ",npts_beam 
               write(*,*)"npts in 3rd axis of cube: ",nz_totpix 
               write(*,*)"  "
               write(*,*)"The number of points along 3rd axis of CUBE"
               write(*,*)"does not match with the number of beam pts!"
               write(*,*)"Check compatibility between files: "
               write(*,*)dty_beam_file(1:nchar(dty_beam_file))," & "
               write(*,*)infile(1:nchar(infile))
               stop 
       endif


        if(mod(nx_totpix,2) .eq. 0)then   
                nxc = nx_totpix/2
                if(cxpix_im .eq. nxc)then
                        nx_1st = nxc - 1
                        nx_2nd = nxc 
                else if(cxpix_im .eq. nxc + 1)then
                        nx_1st = nxc
                        nx_2nd = nxc - 1
                else
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre x-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total x-pixels in image: ",nx_totpix
                        write(*,*)"Expected x-centre : ",nxc,"or",nxc+1
                        write(*,*)"Found x-centre at : ",cxpix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"x-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
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
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre x-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total x-pixels in image: ",nx_totpix
                        write(*,*)"Expected x-centre : ",nxc
                        write(*,*)"Found x-centre at : ",cxpix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"x-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
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
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre y-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total y-pixels in image: ",ny_totpix
                        write(*,*)"Expected y-centre : ",nyc,"or",nyc+1
                        write(*,*)"Found y-centre at : ",cypix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"y-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
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
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre y-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total y-pixels in image: ",ny_totpix
                        write(*,*)"Expected y-centre : ",nyc
                        write(*,*)"Found y-centre at : ",cypix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"y-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
                        if(cypix_im.eq.0)then
                                ny_1st = 0
                                ny_2nd = ny_totpix - 1
                        else
                                ny_1st = cypix_im - 1
                                ny_2nd = ny_totpix - cypix_im
                        endif
                endif
        endif
  
 
        ! For the z-axis
        !
        ! I observe that usually the z-reference value is tagged to the
        ! 1st pixel (referred in some FITS file as 0-th or 1st pixel)
        ! However I also keep a provision to check if the central 
        ! z-pixel has been taken as the reference pixel for crval(3)
        ! 

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
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre z-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total z-pixels in image: ",nz_totpix
                        write(*,*)"Expected z-centre : ",nzc,"or",nzc+1
                        write(*,*)"Found z-centre at : ",czpix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"z-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
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
                        write(*,*)"------------------------------"
                        write(*,*)"           WARNING            "
                        write(*,*)"centre z-pixel is offset from "
                        write(*,*)"actual centre of the image... "
                        write(*,*)" "
                        write(*,*)"Total z-pixels in image: ",nz_totpix
                        write(*,*)"Expected z-centre : ",nzc
                        write(*,*)"Found z-centre at : ",czpix_im
                        write(*,*)" "
                        write(*,*)"Proceed with the fact that the"
                        write(*,*)"z-reference pix is offset from"
                        write(*,*)"the centre of image... "
                        !write(*,*)"Quitting now... "
                        !goto 9999
                        !stop
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
              write(*,*)"Infile chosen:",infile(1:nchar(infile))
              write(*,*)"status = ", status
              write(*,*)"Error opening FITS file..."
              stop
      else
              write(*,*)" "
              write(*,*)"Infile chosen:",infile(1:nchar(infile))
      endif


      !  Create the new RM FITS files. The blocksize parameter is a
      !  historical artifact and the value is ignored by FITSIO.
      call ftinit(41,outfileMAXLOC,blocksize,status)
      call ftinit(42,outfileTHRESH,blocksize,status)
      call ftinit(43,outfileMAXVAL,blocksize,status)
      call ftinit(44,outfileTOTFLUX,blocksize,status)
      call ftinit(45,outfileRMS,blocksize,status)
      call ftinit(46,outfileMEAN,blocksize,status)


      !=======================================================
      ! Main task of the program begins here...

      ! Decide whether the entire cubes need to be read or a
      ! part of them...

      ! By default read the entire cube...
      do i = 1,naxis
         fpixels(i) = 1
         lpixels(i) = naxes(i)
         incs(i) = 1
      enddo
      if(.not.subim)then
              junkchar(1:) = 'nopar'
              write(*,*)" "
              write(*,*)"Entire cube will be used..."
              
              !do i = 1,naxis
              !   fpixels(i) = 1
              !   lpixels(i) = naxes(i)
              !   incs(i) = 1
              !enddo
      else
              write(*,*)" "
              write(*,*)"Sub-section of cube will be used"
              write(*,*)"for the analysis... "
              write(*,*)" "
              subim_parfile = '../PAR/'//subim_parfile(1:nchar(subim_par
     -file))
              write(*,*)"subim_parfile used: ",subim_parfile(8:nchar(sub
     -im_parfile))
              open(201,file=subim_parfile,status='old')
              kk = 0
              do while(1.ne.2)
                 kk = kk + 1
                 read(201,'(a)',end=501)junkchar(1:1)
              enddo

501           continue
              close(201)

              write(*,*)"number of lines in par-file: ",kk - 1
              write(*,*)"naxis: ",naxis 
              !if (kk-1 .ne. naxis + 1)then
              if (kk-1 .ne. naxis + 2)then
                      write(*,*)" "
                      write(*,*)"Cannot determine the sub-image:"
                      write(*,*)"Incomplete or Incompatible parfile..."
                      write(*,*)"Modify the file:",subim_parfile(8:nchar
     -(subim_parfile))
                      write(*,*)"and then re-run the program again."
                      write(*,*)"Quitting now..."
                      write(*,*)" "
                      goto 9999
              else
                      open(201,file=subim_parfile,status='old')
                      read(201,'(a)')junkchar ! Read the first comment line
                      do i = 1,naxis
                         read(201,*)fpixels(i),lpixels(i),incs(i)
                         ! ----------------------------------------
                         ! Make provision to use some of the axis 
                         ! fully while others partly...
                         ! if 0 is supplied, all pix will be read 
                         if(fpixels(i).eq.0)then
                                 fpixels(i) = 1  
                         endif
                         if(lpixels(i).eq.0)then
                                 lpixels(i) = naxes(i) 
                         endif
                         if(incs(i).eq.0)then
                                 incs(i) = 1 
                         endif
                         ! ----------------------------------------
                         if(lpixels(i).lt.fpixels(i))then
                                 write(*,*)" "
                                 write(*,*)"Error: In parfile: "
                                 write(*,*)"last-pix > first-pix"
                                 write(*,*)"In Line number:",i+1
                                 write(*,*)"Correct the par-file: ",subi
     -m_parfile(8:nchar(subim_parfile))
                                 write(*,*)"Quitting now..."
                                 write(*,*)" "
                                 goto 9999
                         endif
                         if(lpixels(i).gt.naxes(i))then
                                 write(*,*)" "
                                 write(*,*)"Error: In parfile: "
                                 write(*,*)"last-pix > ",naxes(i)
                                 write(*,*)"Output image exceeds max dim
     -ension"
                                 write(*,*)"In Line number:",i+1
                                 write(*,*)"Correct the par-file: ",subi
     -m_parfile(8:nchar(subim_parfile))
                                 write(*,*)"Quitting now..."
                                 write(*,*)" "
                                 goto 9999
                         endif
                         if(fpixels(i).lt.1)then
                                 write(*,*)" "
                                 write(*,*)"Error: In parfile: "
                                 write(*,*)"first-pix < 1"
                                 write(*,*)"In Line number:",i+1
                                 write(*,*)"Correct the par-file: ",subi
     -m_parfile(8:nchar(subim_parfile))
                                 write(*,*)"Quitting now..."
                                 write(*,*)" "
                                 goto 9999
                         endif
                         if(incs(i).lt.1)then
                                 write(*,*)" "
                                 write(*,*)"Error: In parfile: "
                                 write(*,*)"inc < 1"
                                 write(*,*)"In Line number:",i+1
                                 write(*,*)"Correct the par-file: ",subi
     -m_parfile(8:nchar(subim_parfile))
                                 write(*,*)"Quitting now..."
                                 write(*,*)" "
                                 goto 9999
                         endif
                      enddo
                      close(201)
              endif
      endif
      

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

      if (nz_use .gt. maxchan)then
              if (subim)then
                      write(*,*)" "
                      write(*,*)"--------------- WARNING --------------"
                      write(*,*)"Number of z-pixels in sub-image excee-"
                      write(*,*)"eded maxchan defined in include file"
                      write(*,*)" "
                      write(*,*)"You may need to modify the 'maxchan' "
                      write(*,*)"parameter in the include file and then"
                      write(*,*)"recompile the code..."
                      write(*,*)" "
                      write(*,*)"Closing the FITS file and Quitting..."
                      write(*,*)" "
                      goto 9999
              else
                      write(*,*)" "
                      write(*,*)"-------------- WARNING ---------------"
                      write(*,*)"Number of z-pixels in image exceeded "
                      write(*,*)"maxchan defined in include file!"
                      write(*,*)" "
                      write(*,*)"You may need to modify the 'maxchan' "
                      write(*,*)"parameter in the include file and then"
                      write(*,*)"recompile the code..."
                      write(*,*)" "
                      write(*,*)"Closing the FITS file and Quitting..."
                      write(*,*)" "
                      goto 9999

              endif
      endif
      ! Now generate the axis values:

      !Refreshing Standard IX Arithmetic-Series concepts: 
      !val_1 = val_c - n_1st*delta_val
      !val_ntot = val_1 + (ntot - 1)*delta_val

      ! generating the x-axis values for the entire image...
  
      x1 = cxval_im - nx_1st*xinc_im
      xn = x1 + (nx_totpix - 1)*xinc_im
      call linspace(x1,xn,nx_totpix,xval)
  
      ! generating the y-axis values for the entire image...
      y1 = cyval_im - ny_1st*yinc_im
      yn = y1 + (ny_totpix - 1)*yinc_im
      call linspace(y1,yn,ny_totpix,yval)

      ! generating the z-axis values for the entire image...


      z1 = czval_im - nz_1st*zinc_im
      zn = z1 + (nz_totpix - 1)*zinc_im
      call linspace(z1,zn,nz_totpix,zval)

      !----------------------------------------------------
      ! Populate the bad channel flags
      do i = 1,nz_totpix 
         flag_arr(i) = 1
      enddo
      ! Now mark the bad-channel flags with 0
      do i = 1,nbad_chan
         tmp_index = bad_chan(i)
         flag_arr(tmp_index) = 0
      enddo
      ngood_chan = 0 
      do i = 1,nz_totpix 
         if(flag_arr(i).eq.1)then
                 ngood_chan = ngood_chan + 1 
                 zval_good(ngood_chan) = zval(i) 
         endif
      enddo


      write(*,*)"xpix-beg,xpix-end,inc: ",xpix_beg,xpix_end,incs(1)
      write(*,*)"ypix-beg,ypix-end,inc: ",ypix_beg,ypix_end,incs(2)
      write(*,*)"zpix-beg,zpix-end,inc: ",zpix_beg,zpix_end,incs(3)

      !  Initialize parameters about the output FITS IMAGE
      !  The EXTEND = TRUE parameter indicates that the FITS file
      !  may contain extensions following the primary array.
      !  Other parameters like BITPIX, naxis, naxes etc., are taken 
      !  to be the same as derived from the input images.

      !  Write the required header keywords to the file

!  subroutine ftphpr(ounit,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)

!! FTPHPR writes required primary header keywords.
!
!       ounit   i  fortran output unit number
!       simple  l  does file conform to FITS standard?
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       pcount  i  number of group parameters
!       gcount  i  number of random groups
!       extend  l  may extensions be present in the FITS file?
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)

      extend= .false.
      simple = .true.
      ! Modify the appropriate headers for output: 
      naxes_out(1) = nx_use
      naxes_out(2) = ny_use
      nz_out = 1
      naxes_out(3) = nz_out
      call ftphpr(41,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(42,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(43,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(44,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(45,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(46,simple,bitpix,3,naxes_out,0,1,extend,status)

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

      decimals = 11
      call ftgkys(21,"ctype1",ctype,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(41,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(42,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(43,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(44,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(45,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(46,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftgkye(21,"crval1",cval,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftgkye(21,"cdelt1",cdelt,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftgkye(21,"crpix1",cpix,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      atmp = cval - real(cpix - 1)*cdelt
      cval = atmp + real(xpix_beg - 1)*cdelt

      call ftpkye(41,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkye(41,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crpix1",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      !cdelt = real(nx_totpix - 1)/real(nx_use - 1)*cdelt
      cdelt = real(incs(1))*cdelt
      call ftpkye(41,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftgkys(21,"ctype2",ctype,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(41,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(42,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(43,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(44,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(45,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(46,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftgkye(21,"crval2",cval,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftgkye(21,"cdelt2",cdelt,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftgkye(21,"crpix2",cpix,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      atmp = cval - real(cpix - 1)*cdelt
      cval = atmp + real(ypix_beg - 1)*cdelt

      call ftpkye(41,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkye(41,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crpix2",1.0,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      !cdelt = real(ny_totpix - 1)/real(ny_use - 1)*cdelt
      cdelt = real(incs(2))*cdelt
      call ftpkye(41,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkys(41,"ctype3","MAXLOC","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(42,"ctype3","THRESH","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(43,"ctype3","MAXVAL","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(44,"ctype3","TOTFLUX","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(45,"ctype3","RMS","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(46,"ctype3","MEAN","3rd axis type",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkye(41,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crval3",1.0,decimals,"Reference Pixel value",
     -status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkye(41,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"crpix3",1.0,decimals,"Reference Pixel",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkye(41,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(42,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(43,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(44,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(45,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkye(46,"cdelt3",1.0,decimals,"Pixel size in world coordina
     -te units",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif

      call ftpkys(41,"BUNIT","rad/m**2","Units of Pixel Data",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(42,"BUNIT","DIMLESS","Units of Pixel Data",status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftgkys(21,"BUNIT",ctype,comment,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(43,"BUNIT",ctype(1:nchar(ctype)),"Units of Pixel Data"
     -,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(44,"BUNIT",ctype(1:nchar(ctype)),"Units of Pixel Data"
     -,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(45,"BUNIT",ctype(1:nchar(ctype)),"Units of Pixel Data"
     -,status)
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      call ftpkys(46,"BUNIT",ctype(1:nchar(ctype)),"Units of Pixel Data"
     -,status)

      ! TEST 
      if(status.ne.0)then
              call printerror(status)
              stop
      endif
      ! A few more useful header info: 
      ! EPOCH of the coordinates:
      call ftgkye(21,"EPOCH",cval,comment,status)
      if(status.ne.0)then
              write(*,*)"Keyword 'EPOCH' missing in "
              write(*,*)"input files' FITS HEADER!"
              write(*,*)" "
              write(*,*)"Default EPOCH assumed: 2000.0"
              write(*,*)"Shall we proceed with default epoch (y/n)? "
              read(*,"(a,$)")yorn
              if(yorn.eq.'y' .or. yorn .eq.'Y')then
                      cval = 2000.0
              else
                      write(*,*)"Enter correct EPOCH (decimal Year): "
                      read(*,*)cval
              endif
              ! Force status to 0
              status = 0
      endif
      call ftpkye(41,"EPOCH",cval,decimals,comment,status)
      call ftpkye(42,"EPOCH",cval,decimals,comment,status)
      call ftpkye(43,"EPOCH",cval,decimals,comment,status)
      call ftpkye(44,"EPOCH",cval,decimals,comment,status)
      call ftpkye(45,"EPOCH",cval,decimals,comment,status)
      call ftpkye(46,"EPOCH",cval,decimals,comment,status)
      ! Object/Field name: 
      call ftgkys(21,"OBJECT",ctype,comment,status)
      call ftpkys(41,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(43,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(44,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(45,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(46,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      ! Scaling if any required: 
      call ftpkye(41,"BSCALE",1.0,decimals," ",status)
      call ftpkye(42,"BSCALE",1.0,decimals,comment,status)
      call ftpkye(43,"BSCALE",1.0,decimals,comment,status)
      call ftpkye(44,"BSCALE",1.0,decimals,comment,status)
      call ftpkye(45,"BSCALE",1.0,decimals,comment,status)
      call ftpkye(46,"BSCALE",1.0,decimals,comment,status)

      call ftpkye(41,"BZERO",0.0,decimals," ",status)
      call ftpkye(42,"BZERO",0.0,decimals,comment,status)
      call ftpkye(43,"BZERO",0.0,decimals,comment,status)
      call ftpkye(44,"BZERO",0.0,decimals,comment,status)
      call ftpkye(45,"BZERO",0.0,decimals,comment,status)
      call ftpkye(46,"BZERO",0.0,decimals,comment,status)
      ! Observer name: 
      call ftgkys(21,"OBSERVER",ctype,comment,status)
      call ftpkys(41,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(43,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(44,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(45,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(46,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      ! TELESCOPE name: 
      call ftgkys(21,"TELESCOP",ctype,comment,status)
      call ftpkys(41,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(43,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(44,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(45,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(46,"TELESCOP",ctype(1:nchar(ctype)),comment,status)




      ! CLOSE THE FITS FILES and open them afresh!
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Qcube..."
      endif

      call FTOPEN(21,infile,rwmode,blocksize,status)
      outfile = outfile(1:nchar(outfile))//'.meanrms'
      open(39,file=outfile,status='unknown')

      write(39,*)(1+int((xpix_end-xpix_beg)/incs(1)))*
     -           (1+int((ypix_end-ypix_beg)/incs(2))), 4 
      write(39,*)"# mean along channels,  rms along channels  "

      ! Initialize aray for average spectra: 
      do i = 1,maxchan  
         zarr(i) = 0.0 
      enddo 

      tmp_cnt1 = 0
      tmp_cnt2 = 0
      cnt1 = 0
      ixpix_now = 0
      do ix = xpix_beg,xpix_end,incs(1)
         ixpix_now = ixpix_now + 1
         if(mod(ixpix_now-1,100).eq.0)then
                 write(*,*)"Doing x-plane: ",ix
         endif

         fpixels(1) = ix
         lpixels(1) = ix

         fpixels(2) = ypix_beg
         lpixels(2) = ypix_end

         fpixels(3) = zpix_beg
         lpixels(3) = zpix_end

         !write(*,*)"fpixels: ",(fpixels(i),i = 1,naxis)
         !write(*,*)"lpixels: ",(lpixels(i),i = 1,naxis)
         call FTGSVE(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specRM,anyflg,status)

         iypix_now = 0
         !write(*,*)"nz_use: ", nz_use 
         do iy = ypix_beg,ypix_end,incs(2)
            iypix_now = iypix_now + 1
            do i = 1,nz_use
            !do i = zpix_beg,zpix_end,incs(3)
               data_arr(i) = specRM(iypix_now + (i-1)*ny_use)
            enddo
            call maxima(data_arr,nz_use,maxxval) ! This should get the 
                                                 ! maxima of the
                                                 ! spectrum 
            !----------------------------------------------------
            ngood_chan = 0
            icnt = 0
            do i = zpix_beg,zpix_end,incs(3)
               icnt = icnt + 1 
               if(flag_arr(i).eq.1)then
                       ngood_chan = ngood_chan + 1
                       !data_good(ngood_chan) = data_arr(i)
                       data_good(ngood_chan) = data_arr(icnt)
                       xarr(ngood_chan) = real(ngood_chan) 
                       beam_arr_good(ngood_chan) = dty_beam(icnt) 
               endif
            enddo
            !write(*,*)(data_good(i),i=1,ngood_chan) 
            !----------------------------------------------------
            ! Compute the statistic along the RM axis: 
            !
            ! Remove linear trend if asked:   
            if(remove_spectral_slope)then
                    write(*,*)"Code Not geared to remove slope..."
                    write(*,*)"This was written for an ACUTE necessity!"
                    write(*,*)"Not all features work..."
                    write(*,*)"Quitting now..."
                    stop 
                    !call fit_linear(xarr,data_good,ngood_chan,A,B)
                    !do i = 1,ngood_chan
                    !   xarr(i) = data_good(i) - B*xarr(i) 
                    !enddo
            else
                    do i = 1,ngood_chan
                       xarr(i) = data_good(i)  
                    enddo
            endif
            ! Now compute the residual : 
            do i = 1,ngood_chan
               !write(54,*)xarr(i),maxxval*beam_arr_good(i)
               xarr(i) = data_good(i) - maxxval*beam_arr_good(i)
            enddo

            !call index_max(data_good,ngood_chan,imax,maxxval,mu,std) 
            call index_max(xarr,ngood_chan,imax,maxxval,mu,std) 
            nuse = 100 
            call rms(xarr,nuse,std)

            maxloc_arr(iypix_now) = zval_good(imax)
            !thresh_arr(iypix_now) = maxxval/std
            thresh_arr(iypix_now) = mu/std
            maxval_arr(iypix_now) = maxxval 
            totflux_arr(iypix_now) = mu*ngood_chan 
            rms_arr(iypix_now) = std 
            mean_arr(iypix_now) = mu 

            write(39,*)ix,iy,mu,std
            ! Accumulate the spectra to compute the mean spectra: 
            ! Use only ON-source pixel: 
            if(thresh_arr(iypix_now).gt.20.0)then
                    cnt1 = cnt1 + 1 
                    do i = 1,ngood_chan 
                       zarr(i) = zarr(i) + xarr(i) 
                    enddo 
            endif

         enddo     ! end of iy loop
         !--------------------------------------------------
         ! Write the FITS CUBES now: 
         fpixels(1) = ixpix_now
         lpixels(1) = ixpix_now
         !write(*,*)"nx_use: ",nx_use

         fpixels(2) = 1 
         lpixels(2) = ny_use !ypix_end - ypix_beg + 1

         fpixels(3) = 1
         lpixels(3) = nz_out

         !write(*,*)"fpixels: ",(fpixels(i),i = 1,naxis)
         !write(*,*)"lpixels: ",(lpixels(i),i = 1,naxis)

         call ftpsse(41,group,3,naxes_out,fpixels,lpixels,
     -                  maxloc_arr,status)
         call ftpsse(42,group,3,naxes_out,fpixels,lpixels,
     -                  thresh_arr,status)
         call ftpsse(43,group,3,naxes_out,fpixels,lpixels,
     -                  maxval_arr,status)
         call ftpsse(44,group,3,naxes_out,fpixels,lpixels,
     -                  totflux_arr,status)
         call ftpsse(45,group,3,naxes_out,fpixels,lpixels,
     -                  rms_arr,status)
         call ftpsse(46,group,3,naxes_out,fpixels,lpixels,
     -                  mean_arr,status)
         !write(*,*)"TEST I am here..."
         !--------------------------------------------------
      enddo        ! end of ix loop
9999  continue
      do i = 1,ngood_chan 
         zarr(i) = zarr(i)/real(cnt1)  
         write(49,*)zarr(i) 
      enddo 
      ! CLOSE THE FITS FILES:
      close(39) 
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing RM-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed RM cube..."
      endif


      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      endif

      call FTCLOS(41,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing MAXLOC-file"
              call printerror(status)
      endif

      call FTCLOS(42,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing THRESH-file"
              call printerror(status)
      endif

      call FTCLOS(43,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing MAXVAL-file"
              call printerror(status)
      endif

      call FTCLOS(44,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing TOTFLUX-file"
              call printerror(status)
      endif

      call FTCLOS(45,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing RMS-file"
              call printerror(status)
      endif

      call FTCLOS(46,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing MEAN-file"
              call printerror(status)
      endif
      ! -----------------------------------------------------------------


      end

c         This subroutine locates the "index" of maxima of 
c         the absolute values of the elements of a vector 
c         "V" given its length "N". It also returns the avg 
c         of the absolute values. 
c         Wasim Raja
c         Date: 12-05-2009
c         

          subroutine index_max(V, N, imax, maxx, avg, std)

          implicit none
          real*4    V(*), avg, maxx, Vnow, std 
          integer*4 N
          integer*4 i, imax

          maxx = V(1)
          imax = 1
          avg = 0.0d0
          do i = 1,N
               Vnow = V(i)
               avg = avg + Vnow
               if(Vnow.ge.maxx)then
                       maxx = Vnow
                       imax = i 
               endif
          end do
          avg = avg/dble(N)

          std = 0.0 
          do i = 1,N
             std = std + (V(i) - avg)**2 
          enddo
          std = sqrt(std/real(N))

          return
          end

      include '/usr/lib/subroutine_lib/nchar.f'
      include 'myfits_info.f'
      !include 'extract_general.f'
      include '/usr/lib/subroutine_lib/extract_general_v3.f'
      include '/usr/lib/subroutine_lib/extract_general_setup.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
