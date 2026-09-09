chelp+
      !-------------------------------------------------------------
      ! This code does the tomography of an image or a rectangular 
      ! subset of it given as inputs the FITS Q and the U spectral 
      ! cubes. Two FITS cubes are written out, one each for linear
      ! polarized intensity as a function of RA, Dec & RM, and 
      ! Polarization Position Angle as a function of RA, Dec & RM
      !                                    -- wr, 19 Aug, 2009
      !-------------------------------------------------------------
chelp-


      !-------------------------------------------------------------
      ! Currently, the tomography is done by reading the input image 
      ! pixel-by-pixel, rather than in one go. That is, the input 
      ! FITS file has to be accessed as many times as the number of 
      ! image-pixels required in the output image. This, although 
      ! is memory efficient, may be highly inefficient in time. 
      ! An alternative would be to read an optimum number of data 
      ! points from the data-array in the input file and keep it in 
      ! memory. The hassle however, in such an attempt would lie in 
      ! interpreting the sequence in which the data-array has been 
      ! written in the input file, since the order of axes in the 
      ! input data is not guaranteed to be consistently followed by 
      ! authors of FITS files.
      !
      !  -- wasim raja, 19 Aug, 2009
      !-------------------------------------------------------------
      ! TODO: KEYWORDS for the output FITS files in cases 
      !       when only a subimage is required, has to be 
      !       appropriately inserted... Currently the keywords 
      !       are read from the INPUT files and copied to the 
      !       output files. Any mismatches encountered thus (for
      !       example, if the output image does not contain 
      !       the reference pixel as defined in the input file),
      !       has been taken care of by writing the FULL image 
      !       with pixels outside the range specified by the 
      !       subimage remaining UNDEFINED (NaN). This scheme 
      !       unnecessarily makes the output images as huge as 
      !       the input images.
      !       -- wasim, 09 Sep, 2009

      ! LAST MODIFICATION: 
      !       --> Configuration file replaces command line arguments.
      !                      -- wr, 15 Sep, 2010
      !       --> Keyword modification for sub-image case DONE!!
      !           Now the output file size is proportional to the 
      !           region of the image used unlike previously where 
      !           the output cubes were as huge as the input cubes. 
      !                      -- wr, 06 Jul, 2011
      !       --> File reading made efficient: Instead of reading 
      !           the Q & U spectra pixel-by-pixel, I now read the 
      !           pixels in the "Dec-Freq" plane for a given RA at 
      !           one go. Care must be taken to interpret the axes 
      !           -- it is assumed in this code that: 
      !                        naxis(1) = RA
      !                        naxis(2) = Dec
      !                        naxis(3) = Freq, in the input cubes.
      !           and that 
      !                        naxis(1) = RA
      !                        naxis(2) = Dec
      !                        naxis(3) = RM, in the input cubes.
      ! 
      ! LAST MODIFICATION: 
      !       --> New parameters added to Configuration file for 
      !           bias removal
      !       --> Stokes-I spectral data cube is now required as 
      !           an input (specified in cfg file) for bias removal 
      !           from Q and U spectra. 
      !                      -- wr, 16 Apr, 2012
      ! 
      ! TODO: NEEDS ATTENTION 
      ! LAST MODIFICATION: 
      !       --> RA mismatch in RM-cubes, whereas Dec matches perfect. 
      !           This bug has been rectified by correcting for the 
      !           declination-dependant dRA : secant(dec) factor!!
      ! 
      !                      -- wr, 11 Jul, 2012
      !-------------------------------------------------------------
      !

      implicit none
      include '../INCLUDE/myfits_spec2rm.inc'

      
      real*4    data_arrI(maxchan),data_arrQ(maxchan), 
     -          data_arrU(maxchan), specI(max_dec*maxchan), 
     -          specQ(max_dec*maxchan), specU(max_dec*maxchan), 
     -          p_ex_arr(maxnt*max_dec), phi_ex_arr(maxnt*max_dec)
      real*4    resiQ, resiU, slopeQ, slopeU  
      logical   remove_QU_bias
      integer*4 bitpixQ, naxisQ, naxesQ(max_axis)
      integer*4 bitpixU, naxisU, naxesU(max_axis)
      integer*4 bitpix, naxis, naxes(max_axis), naxes_out(max_axis)
      logical simple, extend
      integer*4 decimals

      real*4 cxval_im, cyval_im, czval_im
      integer*4 cxpix_im, cypix_im, czpix_im
      real*4 xinc_im, yinc_im, zinc_im

      real*4 cxval_imQ, cyval_imQ, czval_imQ
      integer*4 cxpix_imQ, cypix_imQ, czpix_imQ 
      real*4 xinc_imQ, yinc_imQ, zinc_imQ

      integer*4 xpix_beg, xpix_end
      integer*4 ypix_beg, ypix_end
      integer*4 zpix_beg, zpix_end

      real*4 cxval_imU, cyval_imU, czval_imU  
      integer*4 cxpix_imU, cypix_imU, czpix_imU 
      real*4 xinc_imU, yinc_imU, zinc_imU

      integer*4 nx_totpix, ny_totpix, nz_totpix 
      integer*4 nx_out, ny_out, nz_out, ntot_out
      integer*4 nbuffer, firstpix

      integer*4 fpixels(max_axis), lpixels(max_axis), incs(max_axis)
      real*4 L_sq(maxchan),Q_now(maxchan),U_now(maxchan)
      character*8 junkchar
      integer*4 status, nchar
      logical anyflg
      logical cubeQ
      logical cubeU

      character*64 ctype 
      character*72 comment
      real*8 cval,cdelt, sec_delta, pi 
      real*4 cpix, dRM

      integer*4 rwmode
      character*272 infileI, infileQ, infileU, message
      character*272 outfile, outfileRM, outfilePA
      character*272 RMfile, QU_linecutfile
      character*272 subim_parfile, cfgfile
      character*172 path, path_I 
      character*1 yorn

      integer*4 nx_1st, nx_2nd, ny_1st, ny_2nd, nz_1st, nz_2nd
      integer*4 nxc, nyc, nzc
    
      real*4 xval(max_ra), yval(max_dec), zval(maxchan)
      real*4 x1, xn, y1, yn, z1, zn

      integer*4 data_precision
      real*4 nullval
      logical subim
      real*4 conv_fac ! freq-to-lambda conversion factor
      logical MHz
      ! various counters and indices:
      integer*4 i, kk, ix, iy, ixpix_now, iypix_now, irm 
      integer*4 cnt1, cnt2, tmp_cnt1, tmp_cnt2, tmp_index 


      ! Variables/Parameters for RM-extraction:
      real*4 fac
      integer*4 ofac, rem_mean 
      real*4 RM(maxnt),p_ex(maxnt),phi_ex(maxnt)
      real*4 cos_arr(maxnt,maxchan), sin_arr(maxnt,maxchan)

      ! RFI related (list of bad-channels based on apriori info)
      real*4 bad_chan(maxchan)
      integer*4 nbad_chan, ngood_chan, flag_arr(maxchan)
      logical remove_badchan
      character*172  badchan_file

      ! processing related:
      logical  line_cut
      character*72 add_req
      integer*4 nrm_out 

      ! Some useless fitsio legacy stuff:
      integer*4 group, blocksize

      ! temporary variables: 
      real*4 atmp  
      real*8 atmp8


      pi = acos(-1.0d0)
!-------------------------------------------------------------------
      ! SANITY CHECKS:
      ! Compare the files containing the Q and U Cubes
      ! ans see if they are compatible with each other:
      line_cut = .false.

      if(iargc().lt.1)then
              write(*,*)'  '
              write(*,*)' Usage: '
              write(*,*)'> myfits_spec2rm <cfgfile> <addreq>'
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)' You need a config file containing the '
              write(*,*)' the parameters for this run. '
              write(*,*)'  '
              write(*,*)' You can make some additional requests: '
              write(*,*)' using this string. Valid requests as '
              write(*,*)' of now are: '
              write(*,*)' 1) single_cut: to be used when you intend'
              write(*,*)'                to write out the Q,U and '
              write(*,*)'                RM-spectra for only a single'
              write(*,*)'                "cut" in the sky. By "cut"'
              write(*,*)'                I mean all pixels for eg.,'
              write(*,*)'                having constant Dec value.'
              write(*,*)'   NB: The subim_parfile must be appropriately'
              write(*,*)'       written for this.'
              write(*,*)'  '
              write(*,*)'------------------------------------------'
              write(*,*)'  '
              stop
      else if(iargc().eq.1)then
              call getarg(1,cfgfile)
              cfgfile = cfgfile(1:nchar(cfgfile))
              add_req = 'norequests'
      else if(iargc().gt.1)then
              call getarg(1,cfgfile)
              cfgfile = cfgfile(1:nchar(cfgfile))

              call getarg(2,add_req)
              add_req = add_req(1:nchar(add_req))
      endif


      if(index(add_req,'single_cut').gt.0)then
              line_cut = .true.
      else
              line_cut = .false.
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

      read(11,'(a)')infileQ
      infileQ = infileQ(1:index(infileQ,';')-1)
      infileQ = infileQ(1:nchar(infileQ))
      read(11,'(a)')infileU
      infileU = infileU(1:index(infileU,';')-1)
      infileU = infileU(1:nchar(infileU))
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
      read(11,*)rem_mean    ! removing mean Q and mean U -- mean of 
                            ! the local Q/U spectra -- local to the 
                            ! current pixel. 
                            ! This will basically remove the zero-RM 
                            ! component from each pixel. This
                            ! subtraction is made to take effect in the
                            ! RM-extraction code. 
      if (rem_mean.gt.0)then
              write(*,*)"Mean will be removed from each Q and U "
              write(*,*)"spectra in the RM-extraction..."
              write(*,*)" This will require us to read the I-cube"
              write(*,*)" as well. "
              write(*,*)" "
      endif

      ! New provision for "systematic (instrument?)" removal: 
      ! 
      remove_QU_bias = .false.   ! Initialization 
      read(11,*)yorn
      if(yorn.eq.'y'.or.yorn.eq.'Y')then
              remove_QU_bias = .true.
      else
              remove_QU_bias = .false.
      endif
      read(11,*)resiQ, slopeQ
      read(11,*)resiU, slopeU 
      read(11,'(a)')path_I 
      path_I = path_I(1:index(path_I,';')-1)

      read(11,'(a)')infileI
      infileI = infileI(1:index(infileI,';')-1)
      infileI = infileI(1:nchar(infileI))

      if(remove_qu_bias)then
              write(*,*)"Removing the bias from Q and U..."
              write(*,*)"bias in Q specified in cfg file: ",resiQ 
              write(*,*)"bias in U specified in cfg file: ",resiU 
      else
              write(*,*)"No bias removal from Q and U..."
      endif

      close(11)


      ! Do not write the additional files if the 
      ! entire cube is being processed:
      if(.not.subim)then
              line_cut = .false.
      endif

      infileQ(1:) = path(1:nchar(path))//infileQ(1:nchar(infileQ))
      infileU(1:) = path(1:nchar(path))//infileU(1:nchar(infileU))
      infileI(1:) = path_I(1:nchar(path_I))//infileI(1:nchar(infileI))
      write(*,*)"I-fitscube in: ",infileI(1:nchar(infileI))

      outfileRM(1:) = outfile(1:nchar(outfile))//'.RMCUBE'
      outfilePA(1:) = outfile(1:nchar(outfile))//'.PACUBE'
      QU_linecutfile(1:) = outfile(1:nchar(outfile))//'.QU.linecut'
      badchan_file(1:) = '../PAR/'//badchan_file(1:nchar(badchan_file))

      ! Read the BAD-CHANNELS required for FLAGGING them:
      if(remove_badchan)then
              open(71,file=badchan_file,status='old',err=701)
              goto 702
701           continue
              write(*,*)" "
              write(*,*)"Error opening file: ",
     -                  badchan_file(1:nchar(badchan_file))
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
      read(31,*)ofac  ! Reading the oversampling factor
      read(31,*)fac  ! Reading the uncertainty factor

      close(31)


      ! Extract Some basic INFO from the FITS files:
      call myfits_info(infileQ,
     -           bitpixQ,naxisQ,naxesQ,
     -           cxval_imQ,cxpix_imQ,xinc_imQ,
     -           cyval_imQ,cypix_imQ,yinc_imQ,
     -           czval_imQ,czpix_imQ,zinc_imQ,
     -           cubeQ,message,status)

      if (status.eq.0)then
              write(*,*)"Q-cube opened:",infileQ(1:nchar(infileQ))
              write(*,*)"      bitpixQ:",bitpixQ
              write(*,*)"       naxisQ:",naxisQ
              write(*,*)" "
              write(*,*)"   ref. x-val:",cxval_imQ
              write(*,*)"   ref. x-pix:",cxpix_imQ
              write(*,*)"         xinc:",xinc_imQ
              write(*,*)" "
              write(*,*)"   ref. y-val:",cyval_imQ
              write(*,*)"   ref. y-pix:",cypix_imQ
              write(*,*)"         yinc:",yinc_imQ
              write(*,*)" "
              write(*,*)"   ref. z-val:",czval_imQ
              write(*,*)"   ref. z-pix:",czpix_imQ
              write(*,*)"         zinc:",zinc_imQ
              write(*,*)" "
              write(*,*)"        cubeQ:",cubeQ
              write(*,*)"      message:",message(1:nchar(message))
              do i = 1,naxisQ
                 write(*,*)"naxesQ(",i,") = ",naxesQ(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the Q-cube file as infile"
              write(*,*)"Check if the file exists..."
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif

      call myfits_info(infileU,
     -           bitpixU,naxisU,naxesU,
     -           cxval_imU,cxpix_imU,xinc_imU,
     -           cyval_imU,cypix_imU,yinc_imU,
     -           czval_imU,czpix_imU,zinc_imU,
     -           cubeU,message,status)

      if (status.eq.0)then
              write(*,*)"U-cube opened:",infileU(1:nchar(infileU))
              write(*,*)"      bitpixU:",bitpixU
              write(*,*)"       naxisU:",naxisU
              write(*,*)" "
              write(*,*)"   ref. x-val:",cxval_imU
              write(*,*)"   ref. x-pix:",cxpix_imU
              write(*,*)"         xinc:",xinc_imU
              write(*,*)" "
              write(*,*)"   ref. y-val:",cyval_imU
              write(*,*)"   ref. y-pix:",cypix_imU
              write(*,*)"         yinc:",yinc_imU
              write(*,*)" "
              write(*,*)"   ref. z-val:",czval_imU
              write(*,*)"   ref. z-pix:",czpix_imU
              write(*,*)"         zinc:",zinc_imU
              write(*,*)" "
              write(*,*)"        cubeU:",cubeU
              write(*,*)"      message:",message(1:nchar(message))
              do i = 1,naxisU
                 write(*,*)"naxesU(",i,") = ",naxesU(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the U-cube file as infile"
              write(*,*)"Check if the file exists..."
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif

      write(*,*)"Beginning sanity checks..."
      write(*,*)" "
      if (.not.cubeQ)then
              write(*,*)'ERROR: Image Type mis-match!'
              write(*,*)'    The Q-file is not a cube'
              write(*,*)' '
              write(*,*)'Please ensure that you have input'
              write(*,*)'the right cube-files! '
              write(*,*)' '
              write(*,*)'Quitting now... '
              stop
              !goto 9999
      else if (.not.cubeU)then
              write(*,*)'ERROR: Image Type mis-match!'
              write(*,*)'    The U-file is not a cube'
              write(*,*)' '
              write(*,*)'Please ensure that you have input'
              write(*,*)'the right cube-files! '
              write(*,*)' '
              write(*,*)'Quitting now... '
              stop
              !goto 9999
      else if (naxisQ.ne.naxisU)then
              write(*,*)'ERROR: NAXIS mis-match!'
              write(*,*)'    Q and U-cubes have different'
              write(*,*)'    number of axes...'
              write(*,*)' '
              write(*,*)'Please ensure that you have input'
              write(*,*)'the right cube-files! '
              write(*,*)' '
              write(*,*)'Quitting now... '
              stop
              !goto 9999
      else
              do i = 1,naxisQ
                 if (naxesQ(i).ne.naxesU(i))then
                         write(*,*)' '
                         write(*,*)"ERROR: Axes dimension mis-match"
                         write(*,*)"naxis(",i,") differ in Q and U "
                         write(*,*)' '
                         write(*,*)'Please ensure that you have input'
                         write(*,*)'the right cube-files! '
                         write(*,*)' '
                         write(*,*)'Quitting now... '
                         stop
                         !goto 9999
                 endif
              enddo
              write(*,*)' '
      endif

      ! Check to see if there is a pixel to pixel matching...
      if (cxval_imQ.ne.cxval_imU)then 
              write(*,*)"Reference x-val in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (cxpix_imQ.ne.cxpix_imU)then 
              write(*,*)"Reference x-pix in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (xinc_imQ.ne.xinc_imU)then 
              write(*,*)"x-increment in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (cyval_imQ.ne.cyval_imU)then 
              write(*,*)"Reference y-val in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (cypix_imQ.ne.cypix_imU)then 
              write(*,*)"Reference y-pix in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (yinc_imQ.ne.yinc_imU)then 
              write(*,*)"y-increment in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (czval_imQ.ne.czval_imU)then 
              write(*,*)"Reference z-val in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (czpix_imQ.ne.czpix_imU)then 
              write(*,*)"Reference z-pix in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else if (zinc_imQ.ne.zinc_imU)then 
              write(*,*)"z-increment in the 2 images differ" 
              write(*,*)"This program does not know how to" 
              write(*,*)"handle such files."
              write(*,*)"Ensure that you have used the correct"
              write(*,*)"files... "
              write(*,*)" "
              write(*,*)"Quitting now..."
              stop
      else
              cxval_im = cxval_imQ
              cyval_im = cyval_imQ
              czval_im = czval_imQ

              cxpix_im = cxpix_imQ
              cypix_im = cypix_imQ
              czpix_im = czpix_imQ

              xinc_im = xinc_imQ
              yinc_im = yinc_imQ
              zinc_im = zinc_imQ
              write(*,*)'Q and U-Cubes seem compatible.'
              write(*,*)'We will proceed with the tomography now...'
              write(*,*)' '
      endif


      naxis = naxisQ
      do i = 1,naxis
         naxes(i) = naxesQ(i)
      enddo

      if (bitpixQ.eq.bitpixU)then
              bitpix = bitpixQ
      else
              write(*,*)" "
              write(*,*)"Data types in the Q and U-files differ"
              write(*,*)"Forcing data type to real*4 format..."
              bitpix = -32  ! force real*4 when discrepancy exist
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

      call FTOPEN(21,infileQ,rwmode,blocksize,status)

      if(status.ne.0)then
              write(*,*)" "
              write(*,*)"Q-infile chosen:",infileQ(1:nchar(infileQ))
              write(*,*)"status = ", status
              write(*,*)"Error opening Q-FITS file..."
              stop
      else
              write(*,*)" "
              write(*,*)"Q-infile chosen:",infileQ(1:nchar(infileQ))
      endif

      call FTOPEN(22,infileU,rwmode,blocksize,status)
      if(status.ne.0)then
              write(*,*)"U-infile chosen:",infileU(1:nchar(infileU))
              write(*,*)"status = ", status
              write(*,*)"Error opening U-FITS file..."
              stop
      else
              write(*,*)"U-infile chosen:",infileU(1:nchar(infileU))
      endif


      !  Create the new RM FITS files. The blocksize parameter is a
      !  historical artifact and the value is ignored by FITSIO.

      call ftinit(41,outfileRM,blocksize,status)
      call ftinit(42,outfilePA,blocksize,status)


      !=======================================================
      ! Main task of the program begins here...

      ! Decide whether the entire cubes need to be read or a
      ! part of them...

      if(.not.subim)then
              junkchar(1:) = 'nopar'
              write(*,*)" "
              write(*,*)"Entire Q and U-cubes will be used..."
              
              do i = 1,naxis
                 fpixels(i) = 1
                 lpixels(i) = naxes(i)
                 incs(i) = 1
              enddo
      else
              write(*,*)" "
              write(*,*)"Sub-section of Q and U-cubes will be used"
              write(*,*)"for the tomography... "
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
              if (kk .ne. naxis + 2)then
                      write(*,*)" "
                      write(*,*)"Cannot determine the sub-image:"
                      write(*,*)"Incomplete or Incompatible parfile..."
                      write(*,*)"Modify the file:",subim_parfile(8:nchar
     -(subim_parfile))
                      write(*,*)"and then re-run the program again."
                      write(*,*)"Quitting now..."
                      write(*,*)" "
                      call FTCLOS(21,status)
                      call FTCLOS(22,status)
                      call FTCLOS(41,status)
                      call FTCLOS(42,status)
                      stop 
                      !goto 9999
              else
                      open(201,file=subim_parfile,status='old')
                      read(201,'(a)')junkchar ! Read the first comment line
                      do i = 1,naxis
                         read(201,*)fpixels(i),lpixels(i),incs(i)
                         if(lpixels(i).lt.fpixels(i))then
                                 write(*,*)" "
                                 write(*,*)"Error: In parfile: "
                                 write(*,*)"last-pix > first-pix"
                                 write(*,*)"In Line number:",i+1
                                 write(*,*)"Correct the par-file: ",subi
     -m_parfile(8:nchar(subim_parfile))
                                 write(*,*)"Quitting now..."
                                 write(*,*)" "
                                 call FTCLOS(21,status)
                                 call FTCLOS(22,status)
                                 call FTCLOS(41,status)
                                 call FTCLOS(42,status)
                                 stop 
                                 !goto 9999
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
                                 call FTCLOS(21,status)
                                 call FTCLOS(22,status)
                                 call FTCLOS(41,status)
                                 call FTCLOS(42,status)
                                 stop 
                                 !goto 9999
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
                                 call FTCLOS(21,status)
                                 call FTCLOS(22,status)
                                 call FTCLOS(41,status)
                                 call FTCLOS(42,status)
                                 stop 
                                 !goto 9999
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
                                 call FTCLOS(21,status)
                                 call FTCLOS(22,status)
                                 call FTCLOS(41,status)
                                 call FTCLOS(42,status)
                                 stop 
                                 !goto 9999
                         endif
                      enddo
                      close(201)
              endif
      endif
      
      !write(*,*)" "
      !write(*,'(a)')"junkchar:",junkchar(1:nchar(junkchar))
      !do i = 1,naxis
      !   write(*,*)fpixels(i),lpixels(i),incs(i)
      !enddo

      xpix_beg = fpixels(1)
      xpix_end = lpixels(1)
      nx_out = int((xpix_end - xpix_beg)/incs(1)) + 1

      ypix_beg = fpixels(2)
      ypix_end = lpixels(2)
      ny_out = int((ypix_end - ypix_beg)/incs(2)) + 1

      zpix_beg = fpixels(3)
      zpix_end = lpixels(3)
      nz_out = int((zpix_end - zpix_beg)/incs(3)) + 1

      ntot_out = nx_out*ny_out*nz_out

      if (nz_out .gt. maxchan)then
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
                      call FTCLOS(21,status)
                      call FTCLOS(22,status)
                      call FTCLOS(41,status)
                      call FTCLOS(42,status)
                      stop 
                      !goto 9999
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
                      call FTCLOS(21,status)
                      call FTCLOS(22,status)
                      call FTCLOS(41,status)
                      call FTCLOS(42,status)
                      stop 
                      !goto 9999

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

      ! However, determine the units of Frequency Hz/MHz here...

      ! In the absence of the Frequency-unit information in the 
      ! FITS file, it is a bit tricky to assume it. Here, I 
      ! use a rather BAD trick to tackle the problem:

      MHz = .false.  ! Default

      if (czval_im.ge.30.and.czval_im.le.1.0e4)then ! MHz units assumed
              MHz = .true.
              conv_fac = 300.0
              write(*,*)" "
              write(*,*)"reference-frequency: ",czval_im
              write(*,*)"Assuming frequency in MHz"
              write(*,*)" "
              write(*,*)" "
      else if (czval_im.ge.30.0e6.and.czval_im.le.10.0e9)then ! Hz units assumed
              MHz = .false.
              conv_fac = 3.0e8
              write(*,*)" "
              write(*,*)'reference-frequency: ',czval_im
              write(*,*)"Assuming frequency in Hz"
              write(*,*)" "
      else
              write(*,*)" "
              write(*,*)'reference-frequency: ',czval_im
              write(*,*)"Confusing magnitude for reference-frequency..."
              write(*,*)" "
              write(*,*)"Currently we assume that Hz and MHz are the"
              write(*,*)"ONLY units allowed for Frequency."
              write(*,*)" "
              write(*,*)"Also the observation frequency band is assumed"
              write(*,*)"to be between 30MHz and 10 GHz -- well within"
              write(*,*)"the range of present day Synthesis Radio "
              write(*,*)"Telescopes! "
              write(*,*)" "
              write(*,*)"It would have been a happy situation to "
              write(*,*)"have been able to determine the units from"
              write(*,*)"the FITS file itself. Unfortunately that  "
              write(*,*)"did not seem to be the case for the files"
              write(*,*)"being analysed during the development of "
              write(*,*)"this code :-( "
              write(*,*)"Programmers writing FITS file seem not to "
              write(*,*)"worry about the UNITS of the Axes!! "
              write(*,*)" "
              write(*,*)"If you encounter this message, or feel that "
              write(*,*)"the bug can be removed in a more appealing way"
              write(*,*)",do drop a few lines at:   wasim@rri.res.in  "
              write(*,*)" "
              write(*,*)"Closing the opened FITS files... "
              write(*,*)"Good bye for now... "
              call FTCLOS(21,status)
              call FTCLOS(22,status)
              call FTCLOS(41,status)
              call FTCLOS(42,status)
              stop 
              !goto 9999
      endif
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
         !write(*,*)"bad_chan: ",bad_chan(i)
         tmp_index = bad_chan(i)
         flag_arr(tmp_index) = 0
      enddo

      !-------------------------------------
      ! Set up block for RM-synthesis: 
      !
      ! 1) Arrange the good channels in 
      !    ascending order of lambda_sq: 
      ! Count the good channels: 
      ngood_chan = 0
      do i = zpix_end,zpix_beg,-incs(3)
         if(flag_arr(i).eq.1)then
                 ngood_chan = ngood_chan + 1
                 L_sq(ngood_chan) = (conv_fac/zval(i))**2
         endif
      enddo
      nrm_out = ngood_chan*ofac
      call extract_general_setup(L_sq,ngood_chan,fac, nrm_out, RM, 
     -                           cos_arr, sin_arr, maxnt, maxchan)
      dRM = (RM(nrm_out) - RM(1))/real(nrm_out - 1)
      open(77,file='sampled_RM.txt',status='unknown')
      write(77,*)"# ofac: ",ofac
      do i = 1,nrm_out
         write(77,*)RM(i)
      enddo
      close(77)

      open(78,file='sampled_L_sq_good.txt',status='unknown')
      write(78,*)"# L_sq (only good ones) "
      do i = 1,ngood_chan
         write(78,*)L_sq(i) 
      enddo
      close(78)

      open(79,file='sampled_freq.txt',status='unknown')
      write(79,*)"# freq       L_sq       flag"
      do i = zpix_end,zpix_beg,-incs(3)
         atmp = (conv_fac/zval(i))**2
         write(79,*)zval(i),"    ",atmp,"   ",flag_arr(i) 
      enddo
      close(79)


      ! Irrespective of the total number of output pixels, 
      ! we will read the spectra in the cube on a pix-by-pix 
      ! basis. That way, the variable array named "spec" 
      ! need only be defined to have dimension maxchan.

      write(*,*)"xpix-beg,xpix-end,inc: ",xpix_beg,xpix_end,incs(1)
      write(*,*)"ypix-beg,ypix-end,inc: ",ypix_beg,ypix_end,incs(2)
      write(*,*)"zpix-beg,zpix-end,inc: ",zpix_beg,zpix_end,incs(3)

      !  Initialize parameters about the output FITS CUBES
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
      naxes_out(1) = nx_out
      naxes_out(2) = ny_out
      naxes_out(3) = nrm_out
      call ftphpr(41,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(42,simple,bitpix,3,naxes_out,0,1,extend,status)

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
      decimals = 13

      ! =================== CTYPE-2 =======================
      call ftgkys(21,"ctype2",ctype,comment,status)
      call ftpkys(41,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      call ftpkys(42,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)

      call ftgkyd(21,"crval2",cval,comment,status)
      call ftgkyd(21,"cdelt2",cdelt,comment,status)
      call ftgkye(21,"crpix2",cpix,comment,status)

      atmp8 = cval - (dble(cpix) - 1.0d0)*cdelt     ! Value of the 1st pix 
                                                    ! in input image 
      cval = atmp8 + (dble(ypix_beg) - 1.0d0)*cdelt ! Value of the 1st pix 
                                                    ! in output image 

      ! Now compute the values of center pixel to be used
      ! as the reference for the output image: 
      if(mod(ny_out,2).eq.0)then 
              cpix = real(ny_out/2)
      else
              cpix = real((ny_out+1)/2)
      endif
      cval = cval + dble(cpix - 1)*cdelt*dble(incs(2)) 
!      write(*,'(a,f15.11)')" Ref val Dec in output: ",cval 
!      write(*,'(a,f15.11)')" Ref pix Dec in output: ",cpix 
!      write(*,'(a,f15.11)')" Ref del Dec in output: ",cdelt 


      call ftpkyd(41,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkyd(42,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)

!      call ftpkye(41,"crpix2",1.0,decimals,comment(1:nchar(comment)),
!     -status)
!      call ftpkye(42,"crpix2",1.0,decimals,comment(1:nchar(comment)),
!     -status)
      call ftpkye(41,"crpix2",cpix,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crpix2",cpix,decimals,comment(1:nchar(comment)),
     -status)

      !cdelt = real(ny_totpix - 1)/real(ny_out - 1)*cdelt
      cdelt = dble(incs(2))*cdelt
      call ftpkyd(41,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkyd(42,"cdelt2",cdelt,decimals,comment(1:nchar(comment)),
     -status)

      ! =================== CTYPE-1 =======================
      ! You need to scale the dRA by secant-delta: 
      cval = cval*pi/180.0d0     ! assuming cval is in degrees (dble
                                 ! precision)
      sec_delta = 1.0d0/cos(cval) ! Computed at the reference dec of
                                  ! output image 
      if(sec_delta.ne.sec_delta)then
              write(*,*)"WARNING: "
              write(*,*)"You may be close to the Poles..."
              write(*,*)"We will quit now, since I'm too busy to "
              write(*,*)"think on dealing with this issue now!!"
              stop
      endif

      call ftgkys(21,"ctype1",ctype,comment,status)
      call ftpkys(41,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      call ftpkys(42,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)

      call ftgkyd(21,"crval1",cval,comment,status)
!      write(*,'(a,f15.11)')" Ref val RA in input: ",cval 
      call ftgkyd(21,"cdelt1",cdelt,comment,status)
!      write(*,'(a,f15.11)')"Ref delt RA in input: ",cdelt 
      call ftgkye(21,"crpix1",cpix,comment,status)
!      write(*,'(a,f15.11)')" Ref pix RA in input: ",cpix 
       
      ! You need to scale the dRA by secant-delta: 

      atmp8 = cval - (dble(cpix) - 1.0d0)*cdelt*sec_delta     ! Value of the 1st pix 
                                                              ! in input image at the 
                                                              ! reference declination 
      cval = atmp8 + (dble(xpix_beg) - 1.0d0)*cdelt*sec_delta ! Value of the 1st pix 
                                                              ! in output image at the 
                                                              ! reference declination 
!      write(*,'(a,f15.11)')" 1st val RA in output: ", cval 

      ! Now compute the values of center pixel to be used
      ! as the reference for the output image: 
      if(mod(nx_out,2).eq.0)then 
              cpix = real(nx_out/2)
      else
              cpix = real((nx_out+1)/2)
      endif
      cval = cval + dble(cpix - 1)*cdelt*incs(1)*sec_delta 

      call ftpkyd(41,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkyd(42,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
!      call ftpkye(41,"crpix1",1.0,decimals,comment(1:nchar(comment)),
!     -status)
!      call ftpkye(42,"crpix1",1.0,decimals,comment(1:nchar(comment)),
!     -status)
      call ftpkye(41,"crpix1",cpix,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crpix1",cpix,decimals,comment(1:nchar(comment)),
     -status)

      !cdelt = real(nx_totpix - 1)/real(nx_out - 1)*cdelt
      cdelt = dble(incs(1))*cdelt
      call ftpkyd(41,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkyd(42,"cdelt1",cdelt,decimals,comment(1:nchar(comment)),
     -status)
!      write(*,'(a,f15.11)')" Ref val RA in output: ",cval 
!      write(*,'(a,f15.11)')"Ref delt RA in output: ",cdelt 
!      write(*,'(a,f15.11)')" Ref pix RA in output: ",cpix 
!      !stop
      ! =================== CTYPE-3 =======================

      call ftpkys(41,"ctype3","RM-rd/m2","3rd axis type",status)
      call ftpkys(42,"ctype3","RM-rd/m2","3rd axis type",status)

      call ftpkye(41,"crval3",RM(1),decimals,"Reference Pixel value",
     -status)
      call ftpkye(42,"crval3",RM(1),decimals,"Reference Pixel value",
     -status)

      call ftpkye(41,"crpix3",1.0,decimals,"Reference Pixel",status)
      call ftpkye(42,"crpix3",1.0,decimals,"Reference Pixel",status)

      call ftpkye(41,"cdelt3",dRM,decimals,"Pixel size in world coordina
     -te units",status)
      call ftpkye(42,"cdelt3",dRM,decimals,"Pixel size in world coordina
     -te units",status)

      call ftgkys(21,"BUNIT",ctype,comment,status)
      call ftpkys(41,"BUNIT",ctype(1:nchar(ctype)),"Units of Pixel Data"
     -,status)
      call ftpkys(42,"BUNIT","radians","Units of Pixel Data",status)

      ! A few more useful header info: 
      ! EPOCH of the coordinates:
      call ftgkyd(21,"EPOCH",cval,comment,status)
      if(status.ne.0)then
              write(*,*)"Keyword 'EPOCH' missing in "
              write(*,*)"input files' FITS HEADER!"
              write(*,*)" "
              write(*,*)"Default EPOCH assumed: 2000.0"
              write(*,*)"Shall we proceed with default epoch (y/n)? "
              read(*,"(a,$)")yorn
              if(yorn.eq.'y' .or. yorn .eq.'Y')then
                      cval = 2000.0d0
              else
                      write(*,*)"Enter correct EPOCH (decimal Year): "
                      read(*,*)cval
              endif
              ! Force status to 0
              status = 0
      endif
      call ftpkyd(41,"EPOCH",cval,decimals,comment,status)
      call ftpkyd(42,"EPOCH",cval,decimals,comment,status)
      ! Object/Field name: 
      call ftgkys(21,"OBJECT",ctype,comment,status)
      call ftpkys(41,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"OBJECT",ctype(1:nchar(ctype)),comment,status)
      ! Scaling if any required: 
      call ftpkye(41,"BSCALE",1.0,decimals," ",status)
      call ftpkye(42,"BSCALE",1.0,decimals,comment,status)
      call ftpkye(41,"BZERO",0.0,decimals," ",status)
      call ftpkye(42,"BZERO",0.0,decimals,comment,status)
      ! Observer name: 
      call ftgkys(21,"OBSERVER",ctype,comment,status)
      call ftpkys(41,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"OBSERVER",ctype(1:nchar(ctype)),comment,status)
      ! TELESCOPE name: 
      call ftgkys(21,"TELESCOP",ctype,comment,status)
      call ftpkys(41,"TELESCOP",ctype(1:nchar(ctype)),comment,status)
      call ftpkys(42,"TELESCOP",ctype(1:nchar(ctype)),comment,status)




      RMfile = outfile(1:nchar(outfile))//'.RMSPEC'

      if (line_cut)then
              open(16,file=QU_linecutfile,status='unknown',
     -            form='unformatted',access='direct',recl=4*ngood_chan)
              open(17,file=RMfile,status='unknown',form='unformatted',
     -              access='direct',recl=4*nrm_out)
      endif
      write(*,*)" "

      ! dimx and dimy are the sizes along x and y of the 
      ! data-cubes as defined in the program. FITSIO needs 
      ! to know these dimensions so as to be able to perhaps 
      ! fill any unfilled array element location with nullval

!      ! Read the Q and U cubes into 3D arrays: 
!      write(*,*)"Reading FITS Qcube..."
!      call FTG3DE(21,group, nullval, max_ra, max_dec, 
!     -      nx_totpix, ny_totpix,nz_totpix, data_cubeQ, anyflg, status)
!      write(*,*)"Reading FITS Ucube..."
!      call FTG3DE(22,group, nullval, max_ra, max_dec, 
!     -      nx_totpix, ny_totpix,nz_totpix, data_cubeU, anyflg, status)

      ! CLOSE THE FITS FILES and open them afresh!
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Qcube..."
      endif
      call FTCLOS(22,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing U-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Ucube..."
      endif

      call FTOPEN(21,infileQ,rwmode,blocksize,status)
      call FTOPEN(22,infileU,rwmode,blocksize,status)
      ! Open the I-cube as well for calibration purpose : wr, 18 Apr, 2012
      call FTOPEN(40,infileI,rwmode,blocksize,status)

      if(line_cut)then
              open(121,file='rm_spec.txt')
      endif

      tmp_cnt1 = 0
      tmp_cnt2 = 0
      cnt1 = 0
      ixpix_now = 0
      do ix = xpix_beg,xpix_end,incs(1)
         ixpix_now = ixpix_now + 1
         write(*,*)"Doing x-plane: ",ix

         fpixels(1) = ix
         lpixels(1) = ix

         fpixels(2) = ypix_beg
         lpixels(2) = ypix_end

         fpixels(3) = zpix_beg
         lpixels(3) = zpix_end

         !write(*,*)"fpixels: ",(fpixels(i),i = 1,naxis)
         !write(*,*)"lpixels: ",(lpixels(i),i = 1,naxis)
         call FTGSVE(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
         call FTGSVE(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
         call FTGSVE(40,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specI,anyflg,status)
         ! TEST
         !do i = 1,ny_out*nz_out
         !    write(88,*)specQ(i), specU(i)
         !enddo

         iypix_now = 0 
         irm = 0 
         do iy = ypix_beg,ypix_end,incs(2)
            iypix_now = iypix_now + 1
            cnt1 = cnt1 + 1
            !write(87,*)"## RA, Dec",ix,iy
            do i = 1,nz_out
               !data_arrQ(i) = specQ(i + (iy-1)*nz_out)
               !data_arrU(i) = specU(i + (iy-1)*nz_out)
               data_arrQ(i) = specQ(iypix_now + (i-1)*ny_out)
               data_arrU(i) = specU(iypix_now + (i-1)*ny_out)
               data_arrI(i) = specI(iypix_now + (i-1)*ny_out)
               ! TEST
               !write(87,*)data_arrQ(i), data_arrU(i)
            enddo
            !----------------------------------------------------
            ngood_chan = 0
            cnt2 = nz_out + 1
            !do i = zpix_end,zpix_beg,-incs(3)
            !   cnt2 = cnt2 - 1
            !   if(flag_arr(i).eq.1)then
            !           ngood_chan = ngood_chan + 1
            !           Q_now(ngood_chan) = data_arrQ(cnt2)
            !           U_now(ngood_chan) = data_arrU(cnt2)
            !   endif
            !enddo
            if(.not.remove_QU_bias)then
                    do i = zpix_end,zpix_beg,-incs(3)
                       cnt2 = cnt2 - 1
                       if(flag_arr(i).eq.1)then
                          ngood_chan = ngood_chan + 1
                          Q_now(ngood_chan) = data_arrQ(cnt2) 
                          U_now(ngood_chan) = data_arrU(cnt2) 
                       endif
                    enddo
            else
                    do i = zpix_end,zpix_beg,-incs(3)
                       cnt2 = cnt2 - 1
                       if(flag_arr(i).eq.1)then
                          ngood_chan = ngood_chan + 1
                          if(data_arrQ(cnt2).ge.resiQ)then
                                  slopeQ = slopeQ
                          else
                                  slopeQ = -slopeQ
                          endif
                          if(data_arrU(cnt2).ge.resiU)then
                                  slopeU = slopeU
                          else
                                  slopeU = -slopeU
                          endif
                          Q_now(ngood_chan) = data_arrQ(cnt2) - 
     -                               (data_arrI(cnt2)*slopeQ + resiQ) 
                          U_now(ngood_chan) = data_arrU(cnt2) -  
     -                               (data_arrI(cnt2)*slopeU + resiU) 
                       endif
                    enddo
            endif
            !----------------------------------------------------

            ! Perform the tomography now:
!            call extract_general(Q_now,U_now,ngood_chan,
!     -                           nrm_out, p_ex,phi_ex,
!     -                           cos_arr,sin_arr, maxnt, maxchan)
            call extract_general(Q_now,U_now,ngood_chan,
     -                           nrm_out, p_ex,phi_ex,
     -                           cos_arr,sin_arr, maxnt, maxchan, 
     -                           rem_mean)

            ! Fill up p_ex and phi_ex arrays for writing out 
            ! later (we wish to optimally access FITS FILES and 
            ! save time:)
            do i = 1,nrm_out
               irm = irm + 1
               p_ex_arr(iypix_now + (i-1)*ny_out) = p_ex(i)
               phi_ex_arr(iypix_now + (i-1)*ny_out) = phi_ex(i)
               !p_ex_arr(irm) = p_ex(i)
               !phi_ex_arr(irm) = phi_ex(i)
            enddo

            !--------------------------------------------------
            ! Write Q and U in ascending order of frequency and 
            ! not ascending order of wavelength-squared: 
            ! Q(i) and U(i) correspond to i-th Freq channel
            ! p_ex(i) and phi_ex(i) correspond to i-th RM-bin
            if(line_cut)then
                    tmp_cnt1 = tmp_cnt1 + 1
                    write(16,rec=tmp_cnt1)(Q_now(i),i=ngood_chan,1,-1)
                    tmp_cnt1 = tmp_cnt1 + 1
                    write(16,rec=tmp_cnt1)(U_now(i),i=ngood_chan,1,-1)

                    tmp_cnt2 = tmp_cnt2 + 1
                    write(17,rec=tmp_cnt2)(p_ex(i),i=1,nrm_out)
                    tmp_cnt2 = tmp_cnt2 + 1
                    write(17,rec=tmp_cnt2)(phi_ex(i),i=1,nrm_out)
                    ! TEST : 
                    ! Readily output the RM-spectra if a single pixel 
                    ! tomography done: 
                    write(121,*)"## ix, iy: ",ix,iy
                    do i = 1,nrm_out
                      write(121,*)p_ex(i), phi_ex(i)
                    enddo
            endif
            !! Write the FITS RM-CUBE now:
            if(mod(cnt1-1,1000).eq.0)then
                 write(*,*)"doing ",cnt1," out of",nx_out*ny_out 
            endif
         enddo     ! end of iy loop
         !--------------------------------------------------
         ! Write the FITS CUBES now: 
         fpixels(1) = ixpix_now
         lpixels(1) = ixpix_now
         !write(*,*)"nx_out: ",nx_out

         fpixels(2) = 1 
         lpixels(2) = ny_out !ypix_end - ypix_beg + 1

         fpixels(3) = 1
         lpixels(3) = nrm_out

         !write(*,*)"fpixels: ",(fpixels(i),i = 1,naxis)
         !write(*,*)"lpixels: ",(lpixels(i),i = 1,naxis)

         call ftpsse(41,group,3,naxes_out,fpixels,lpixels,
     -                  p_ex_arr,status)
         if(status.gt.0)then
                 call printerror(status) 
         endif
         call ftpsse(42,group,3,naxes_out,fpixels,lpixels,
     -                  phi_ex_arr,status)
         if(status.gt.0)then
                 call printerror(status) 
         endif
         !write(*,*)"TEST I am here..."
         !--------------------------------------------------
      enddo        ! end of ix loop
      if(line_cut)then
              close(121) 
      endif
      ! CLOSE THE FITS FILES:
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Qcube..."
      endif
      call FTCLOS(22,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing U-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Ucube..."
      endif
      call FTCLOS(40,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing I-file"
              call printerror(status)
      else
              write(*,*)"Successfully read and closed FITS Icube..."
      endif
      write(*,*)" ================================"
      write(*,*)"      fac :",fac
      write(*,*)"ngood_chan: ", ngood_chan 
      write(*,*)"   nRM_out: ", nRM_out
      write(*,*)"      cnt1: ", cnt1
      write(*,*)"       RM1: ", RM(1)
      write(*,*)"       RM2: ", RM(nrm_out)
      write(*,*)"       dRM: ", dRM
      !=======================================================


9999  continue
      if(line_cut)then
              close(16)
              close(17)
      endif

!      !write(*,*)"---------------------------"
!      !write(*,*)"Current subroutine: FTCLOS "
!      !write(*,*)"STATUS = ",status
!      !write(*,*)" "
!      if (status .gt. 0)then
!              write(*,*)"Problem closing Q-file"
!              call printerror(status)
!      endif
!
!      ! CLOSE THE FITS FILES:
!      if (status .gt. 0)then
!              write(*,*)"Problem closing Q-file"
!              call printerror(status)
!      endif
!      if(.not.line_cut)then
              call FTCLOS(41,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing RM-file"
                      call printerror(status)
              endif
              call FTCLOS(42,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing PA-file"
                      call printerror(status)
              endif
!      endif

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
