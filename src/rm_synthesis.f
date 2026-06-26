chelp+
      !-------------------------------------------------------------
      ! This code does the tomography of an image or a rectangular 
      ! subset of it given as inputs the FITS Q and the U spectral 
      ! cubes. Two FITS cubes are written out, one each for linear
      ! polarized intensity as a function of RA, Dec RM, and 
      ! Polarization Position Angle as a function of RA, Dec RM
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
      !           the Q U spectra pixel-by-pixel, I now read the 
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

      use rm_synthesis_mod
      implicit none
      
      
      real(sp), allocatable :: data_arrI(:)
      real(sp), allocatable :: data_arrQ(:)
      real(sp), allocatable :: data_arrU(:)
      real(sp), allocatable :: specMask(:)
      real(sp), allocatable :: specI(:)
      real(sp), allocatable :: specQ(:)
      real(sp), allocatable :: specU(:)
      real(sp), allocatable :: p_tile_arr(:)
      real(sp), allocatable :: phi_tile_arr(:)
      integer*1, allocatable :: mask_tile_arr(:)
      integer*2, allocatable :: nvalid_tile_arr(:)
      real(sp)  resiQ, resiU, slopeQ, slopeU
      logical   remove_QU_bias
      integer   bitpixQ, naxisQ, naxesQ(max_axis)
      integer   bitpixU, naxisU, naxesU(max_axis)
      integer   bitpixM, naxisM, naxesM(max_axis)
      integer   bitpix, naxis, naxes(max_axis), naxes_out(max_axis)
      logical   simple, extend
      integer   decimals

      real(sp) cxval_im, cyval_im, czval_im
      integer   cxpix_im, cypix_im, czpix_im
      real(sp) xinc_im, yinc_im, zinc_im

      real(sp) cxval_imQ, cyval_imQ, czval_imQ
      integer   cxpix_imQ, cypix_imQ, czpix_imQ 
      real(sp) xinc_imQ, yinc_imQ, zinc_imQ

      integer   xpix_beg, xpix_end
      integer   ypix_beg, ypix_end
      integer   zpix_beg, zpix_end
      integer   subim_ra_blc, subim_ra_trc, subim_ra_inc
      integer   subim_dec_blc, subim_dec_trc, subim_dec_inc
      integer   subim_chan_blc, subim_chan_trc, subim_chan_inc
      integer   tile_ra, tile_dec

      real(sp) cxval_imU, cyval_imU, czval_imU  
      integer   cxpix_imU, cypix_imU, czpix_imU 
      real(sp) xinc_imU, yinc_imU, zinc_imU
      real(sp) cxval_imM, cyval_imM, czval_imM
      integer   cxpix_imM, cypix_imM, czpix_imM
      real(sp) xinc_imM, yinc_imM, zinc_imM

      integer   nx_totpix, ny_totpix, nz_totpix 
      integer   nx_out, ny_out, nz_out, ntot_out
      integer   nbuffer, firstpix

      integer   fpixels(max_axis), lpixels(max_axis), incs(max_axis)
        integer   fpixels_out(3), lpixels_out(3)
      real(sp), allocatable :: L_sq(:)
      real(sp), allocatable :: Q_now(:)
      real(sp), allocatable :: U_now(:)
      character(len=8) :: junkchar
      integer   status
      logical   anyflg
      logical   cubeQ
      logical   cubeU
      logical   cubeM
      logical   out_amp_open, out_ang_open, out_exists
      logical   out_mask_open, out_nvalid_open
      integer   freq_axis, freq_axisQ, freq_axisU
      integer   freq_axisM

      character(len=64) :: ctype 
      character(len=72) :: comment
      real(dp) cval,cdelt, pi 
      real(sp) cpix, dRM

      integer   rwmode
      character(len=272) :: infileI, infileQ, infileU, message
      character(len=272) :: outfile, outfileAMP, outfileANG
      character(len=272) :: outfileMASK, outfileNVALID
      character(len=272) :: mask_cube_file, mask_input_cube_file,
     -                      mask_trust_mode
      character(len=272) :: RMfile, QU_linecutfile
      character(len=272) :: subim_parfile, cfgfile, cfgfile_in
      character(len=172) :: path, path_I 
      character(len=1) :: yorn

      integer   nx_1st, nx_2nd, ny_1st, ny_2nd, nz_1st, nz_2nd
      integer   nxc, nyc, nzc
    
      real(sp), allocatable :: xval(:)
      real(sp), allocatable :: yval(:)
      real(sp), allocatable :: zval(:)
      real(sp) x1, xn, y1, yn, z1, zn

      integer   data_precision
      real(sp) nullval
      logical   subim
      logical   tile_auto, dry_run
      real(sp) conv_fac ! freq-to-lambda conversion factor
      real(sp) tile_mem_frac
      logical   MHz
      ! various counters and indices:
      integer   i, kk, ix, iy, ixpix_now, iypix_now, irm 
      integer   cnt1, cnt2, tmp_cnt1, tmp_cnt2, tmp_index 
      integer   progress_total, progress_step
      integer   progress_next_pct, progress_next_count
      integer   ix_tile_beg, ix_tile_end, iy_tile_beg, iy_tile_end
      integer   ix_loc, iy_loc, iz
      integer   nx_tile, ny_tile
      integer   ix_out_beg, ix_out_end, iy_out_beg, iy_out_end
      integer   cnt_good, nvalid_pix
      integer   fpixels_nvalid(2), lpixels_nvalid(2)
      integer   naxes_mask(3), naxes_nvalid(2)
      logical   nan_check_on, chan_valid
      logical   use_input_mask, in_mask_open
      real(sp)  mask_val
      integer   in_fields
      integer   mem_unit, ios_mem
      integer(kind=int64) :: mem_avail_kb, mem_kb_tmp
      integer(kind=int64) :: mem_safe_bytes, bytes_per_tile_pixel
      integer(kind=int64) :: tile_pixels_max, tile_bytes_est
      integer(kind=int64) :: image_pixels_total
      character(len=256) :: mem_line


      ! Variables/Parameters for RM-extraction:
      real(sp) fac, beg_rm, end_rm
      integer   ofac, rem_mean, nrm_out, nrm_out_par
      integer   use_auto_rm_range
      integer   output_mode
      integer   ap_angle_mode
      real(sp), allocatable :: RM(:)
      real(sp), allocatable :: p_ex(:)
      real(sp), allocatable :: phi_ex(:)
      real(sp), allocatable :: cos_arr(:,:)
      real(sp), allocatable :: sin_arr(:,:)
      real(sp), allocatable :: wts_now(:)

      ! RFI related (list of bad-channels based on apriori info)
        real(sp), allocatable :: bad_chan(:)
      integer   nbad_chan, ngood_chan
      integer, allocatable :: flag_arr(:)
      logical   remove_badchan
      character(len=172) :: global_badchan_file
      character(len=16) :: masksrc_key, nanchk_key

      ! processing related:
      logical   line_cut
      logical   need_icube
      character(len=72) :: add_req

      ! Some useless fitsio legacy stuff:
      integer   group, blocksize

      ! temporary variables: 
      real(sp) atmp  
      real(dp) atmp8


      pi = acos(-1.0d0)
!-------------------------------------------------------------------
      ! SANITY CHECKS:
      ! Compare the files containing the Q and U Cubes
      ! ans see if they are compatible with each other:
      line_cut = .false.

      if(command_argument_count() < 1)then
              write(*,*)'  '
              write(*,*)' Usage: '
              write(*,*)'> rm_synthesis <cfgfile> <addreq>'
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
      else if(command_argument_count() == 1)then
              call get_command_argument(1, cfgfile)
              cfgfile = cfgfile(1:nchar(cfgfile))
              add_req = 'norequests'
      else if(command_argument_count() > 1)then
              call get_command_argument(1, cfgfile)
              cfgfile = cfgfile(1:nchar(cfgfile))

              call get_command_argument(2, add_req)
              add_req = add_req(1:nchar(add_req))
      endif


      if(index(add_req,'single_cut') > 0)then
              line_cut = .true.
      else
              line_cut = .false.
      endif

      cfgfile_in = cfgfile(1:nchar(cfgfile))
      cfgfile = cfgfile_in(1:nchar(cfgfile_in))
      inquire(file=cfgfile(1:nchar(cfgfile)),exist=anyflg)
      if(.not.anyflg)then
              cfgfile = '../cfg/'//cfgfile_in(1:nchar(cfgfile_in))
              inquire(file=cfgfile(1:nchar(cfgfile)),exist=anyflg)
      endif
      if(.not.anyflg)then
              write(*,*)"Error locating config file: "
              write(*,*)cfgfile_in(1:nchar(cfgfile_in))
              write(*,*)"Tried direct path and ../cfg/ only"
              stop
      endif

      call read_cfg_keyval(cfgfile,
     -          path,infileQ,infileU,outfile,
     -          remove_badchan,global_badchan_file,
     -          subim,subim_parfile,
     -          subim_ra_blc,subim_ra_trc,subim_ra_inc,
     -          subim_dec_blc,subim_dec_trc,subim_dec_inc,
     -          subim_chan_blc,subim_chan_trc,subim_chan_inc,
     -          tile_ra,tile_dec,tile_mem_frac,tile_auto,dry_run,
     -          rem_mean,remove_QU_bias,
     -          resiQ,slopeQ,resiU,slopeU,
     -          path_I,infileI,
     -          ofac,fac,beg_rm,end_rm,nrm_out_par,
     -          use_auto_rm_range,output_mode,
     -          ap_angle_mode,mask_cube_file,
     -          mask_input_cube_file,
     -          mask_trust_mode,status)
      if(status.ne.0)then
              write(*,*)"Error opening/parsing config file: "
              write(*,*)cfgfile(1:nchar(cfgfile))
              write(*,*)"Quitting now..."
              stop
      endif

      if (rem_mean.gt.0)then
              write(*,*)"Mean will be removed from each Q and U "
              write(*,*)"spectra in the RM-extraction..."
              write(*,*)" "
      endif

      need_icube = .false.
      if(remove_qu_bias)need_icube = .true.

      if(remove_qu_bias)then
              write(*,*)"Removing the bias from Q and U..."
              write(*,*)"bias in Q specified in cfg file: ",resiQ
              write(*,*)"bias in U specified in cfg file: ",resiU
      else
              write(*,*)"No bias removal from Q and U..."
      endif


      ! Do not write the additional files if the 
      ! entire cube is being processed:
      if(.not.subim)then
              line_cut = .false.
      endif

      infileQ(1:) = path(1:nchar(path))//infileQ(1:nchar(infileQ))
      infileU(1:) = path(1:nchar(path))//infileU(1:nchar(infileU))
      if(nchar(mask_input_cube_file).gt.0)then
              inquire(file=mask_input_cube_file(
     -               1:nchar(mask_input_cube_file)),exist=anyflg)
              if(.not.anyflg)then
                      mask_input_cube_file(1:) =
     -                      path(1:nchar(path))//
     -                      mask_input_cube_file(
     -                      1:nchar(mask_input_cube_file))
              endif
      endif
      if(need_icube)then
              infileI(1:)=path_I(1:nchar(path_I))//
     -       infileI(1:nchar(infileI))
              write(*,*)"I-fitscube in: ",infileI(1:nchar(infileI))
      endif

      outfileAMP(1:) = outfile(1:nchar(outfile))//'.AMP.RMCUBE.FITS'
      if(output_mode.eq.1)then
              outfileAMP(1:) = outfile(1:nchar(outfile))//
     -              '.REAL.RMCUBE.FITS'
              outfileANG(1:) = outfile(1:nchar(outfile))//
     -              '.IMAG.RMCUBE.FITS'
      else
              if(ap_angle_mode.eq.1)then
                      outfileANG(1:) = outfile(1:nchar(outfile))//
     -                      '.POLA.RMCUBE.FITS'
              else
                      outfileANG(1:) = outfile(1:nchar(outfile))//
     -                      '.PHA.RMCUBE.FITS'
              endif
      endif
      outfileMASK(1:) = outfile(1:nchar(outfile))//'.MASK.CUBE.FITS'
      if(nchar(mask_cube_file).gt.0)then
              outfileMASK(1:) = mask_cube_file(1:nchar(mask_cube_file))
      endif
      outfileNVALID(1:) = outfile(1:nchar(outfile))//'.NVALID.MAP.FITS'
      QU_linecutfile(1:) = outfile(1:nchar(outfile))//'.QU.linecut'
        global_badchan_file(1:) = global_badchan_file(
     -                       1:nchar(global_badchan_file))

      ! Bad channels will be read after cube dimensions are known
                        nbad_chan = 0

      write(*,*)' ========================================'
      write(*,*)' RM Extraction Parameters from Config:'
      if(output_mode.eq.1)then
              write(*,*)' output_mode: ri'
      else
              write(*,*)' output_mode: ap'
      endif
      write(*,*)' ofac: ',ofac
      write(*,*)' fac:  ',fac
      write(*,*)' use_auto_rm_range: ',use_auto_rm_range
      if(output_mode.eq.0)then
              if(ap_angle_mode.eq.1)then
                      write(*,*)' ap_angle_mode: pol'
              else
                      write(*,*)' ap_angle_mode: phase'
              endif
      endif
      if(use_auto_rm_range.eq.0)then
              write(*,*)' beg_rm: ',beg_rm
              write(*,*)' end_rm: ',end_rm
              write(*,*)' nrm: ',nrm_out_par
              write(*,*)' nrm_out (nrm*ofac): ',nrm_out_par*ofac
      else
              write(*,*)' beg/end/nrm are auto-derived from data'
      endif
      write(*,*)' ========================================'

      ! Extract Some basic INFO from the FITS files:
      call myfits_info(infileQ,
     -           bitpixQ,naxisQ,naxesQ,
     -           cxval_imQ,cxpix_imQ,xinc_imQ,
     -           cyval_imQ,cypix_imQ,yinc_imQ,
     -           czval_imQ,czpix_imQ,zinc_imQ,
     -           freq_axisQ,cubeQ,message,status)

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
              write(*,*)"   freq-axisQ:",freq_axisQ
              write(*,*)"      message:",message(1:nchar(message))
              do i = 1,naxisQ
                 write(*,*)"naxesQ(",i,") = ",naxesQ(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the Q-cube file as infile"
              write(*,*)"message:",message(1:nchar(message))
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif

      call myfits_info(infileU,
     -           bitpixU,naxisU,naxesU,
     -           cxval_imU,cxpix_imU,xinc_imU,
     -           cyval_imU,cypix_imU,yinc_imU,
     -           czval_imU,czpix_imU,zinc_imU,
     -           freq_axisU,cubeU,message,status)

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
              write(*,*)"   freq-axisU:",freq_axisU
              write(*,*)"      message:",message(1:nchar(message))
              do i = 1,naxisU
                 write(*,*)"naxesU(",i,") = ",naxesU(i)
              enddo
      else
              write(*,*)"status = ",status
              write(*,*)"something went wrong with the "
              write(*,*)"'myfits_info' subroutine call"
              write(*,*)"with the U-cube file as infile"
              write(*,*)"message:",message(1:nchar(message))
              write(*,*)"Quitting now..."
              stop
              !goto 9999
      endif

      write(*,*)"Beginning sanity checks..."
      write(*,*)" "
      if (.not.cubeQ)then
              write(*,*)'ERROR: Missing spectral axis in Q-file!'
              write(*,*)'    No CTYPE*=FREQ axis detected.'
              write(*,*)' '
              write(*,*)'Please ensure that you have input'
              write(*,*)'a FITS file with a frequency axis. '
              write(*,*)' '
              write(*,*)'Quitting now... '
              stop
              !goto 9999
      else if (.not.cubeU)then
              write(*,*)'ERROR: Missing spectral axis in U-file!'
              write(*,*)'    No CTYPE*=FREQ axis detected.'
              write(*,*)' '
              write(*,*)'Please ensure that you have input'
              write(*,*)'a FITS file with a frequency axis. '
              write(*,*)' '
              write(*,*)'Quitting now... '
              stop
              !goto 9999
      else if (freq_axisQ.ne.freq_axisU)then
              write(*,*)'ERROR: Frequency-axis index mis-match!'
              write(*,*)'    Q FREQ axis = ',freq_axisQ
              write(*,*)'    U FREQ axis = ',freq_axisU
              write(*,*)' '
              write(*,*)'Please reorder Q and U cubes to use'
              write(*,*)'the same frequency-axis placement.'
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
      freq_axis = freq_axisQ

      if(use_input_mask)then
              status = 0
              call myfits_info(mask_input_cube_file,
     -           bitpixM,naxisM,naxesM,
     -           cxval_imM,cxpix_imM,xinc_imM,
     -           cyval_imM,cypix_imM,yinc_imM,
     -           czval_imM,czpix_imM,zinc_imM,
     -           freq_axisM,cubeM,message,status)
              if(status.ne.0)then
                      write(*,*)"status = ",status
                      write(*,*)"Mask cube info read failed"
                      write(*,*)mask_input_cube_file(
     -                         1:nchar(mask_input_cube_file))
                      write(*,*)"message:",message(1:nchar(message))
                      stop
              endif
              if(.not.cubeM)then
                      write(*,*)'ERROR: Missing FREQ axis in mask cube!'
                      write(*,*)'No CTYPE*=FREQ axis detected.'
                      stop
              endif
              if(freq_axisM.ne.freq_axis)then
                      write(*,*)'ERROR: Mask/Q-U FREQ axis mismatch!'
                      write(*,*)'Mask FREQ axis = ',freq_axisM
                      write(*,*)'Q/U  FREQ axis = ',freq_axis
                      stop
              endif
              if(naxisM.ne.3 .and. naxisM.ne.naxisQ)then
                      write(*,*)'ERROR: Unsupported mask NAXIS!'
                      write(*,*)'Mask NAXIS = ',naxisM
                      write(*,*)'Expected NAXIS = 3 or ',naxisQ
                      stop
              endif
              if(naxesM(1).ne.naxesQ(1))then
                      write(*,*)'ERROR: Mask RA size mismatch!'
                      write(*,*)'Mask RA length = ',naxesM(1)
                      write(*,*)'Q/U  RA length = ',naxesQ(1)
                      stop
              endif
              if(naxesM(2).ne.naxesQ(2))then
                      write(*,*)'ERROR: Mask Dec size mismatch!'
                      write(*,*)'Mask Dec length = ',naxesM(2)
                      write(*,*)'Q/U  Dec length = ',naxesQ(2)
                      stop
              endif
              if(naxesM(freq_axisM).ne.naxesQ(freq_axis))then
                      write(*,*)'ERROR: Mask FREQ size mismatch!'
                      write(*,*)'Mask FREQ length = ',naxesM(freq_axisM)
                      write(*,*)'Q/U  FREQ length = ',naxesQ(freq_axis)
                      stop
              endif
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


        write(*,*)"! -----------------------------------------------"
        write(*,*)"! Final sanity checks..."
        write(*,*)"! Reference-pixel conventions are validated below."
        write(*,*)"! -----------------------------------------------"

        ! Check if the reference pixel is indeed at the centre of the 
        ! image array and also find out the number of points leading 
        ! and lagging the reference pixel:
 


        ! For the x-axis
        ! If total number of pixel is EVEN, centre must lie at n_totpix/2 
        ! or n_totpix/2 + 1



       nx_totpix = naxes(1)
       ny_totpix = naxes(2)
       nz_totpix = naxes(freq_axis)

      ! Allocate axis and flag arrays now cube dimensions are known
      allocate(xval(nx_totpix))
      allocate(yval(ny_totpix))
      allocate(zval(nz_totpix))
      allocate(flag_arr(nz_totpix))
      allocate(bad_chan(nz_totpix))

      ! Read the bad channel list now that nz_totpix is known
      if(remove_badchan)then
              open(71,file=global_badchan_file(
     -                 1:nchar(global_badchan_file)),
     -             status='old',iostat=ios_mem)
              if(ios_mem .ne. 0)then
                      write(*,*)"Error opening bad channel file:"
                      write(*,*)global_badchan_file(
     -                         1:nchar(global_badchan_file))
                      write(*,*)"Skipping bad channel flagging."
                      remove_badchan = .false.
              else
                      nbad_chan = 0
                      do while(.true.)
                              nbad_chan = nbad_chan + 1
                              if(nbad_chan .gt. nz_totpix)then
                                      write(*,*)"Too many bad channels"
                                      write(*,*)"Max by cube:"
                                      write(*,*)nz_totpix
                                      close(71)
                                      stop
                              endif
                              read(71,*,end=711)bad_chan(nbad_chan)
                              write(*,*)"bad-chan: ",bad_chan(nbad_chan)
                      enddo
711                   continue
                      nbad_chan = nbad_chan - 1
                      write(*,*)"Number of Bad Channels: ",nbad_chan
                      close(71)
              endif
      endif

      masksrc_key = 'generated'
      if(use_input_mask)masksrc_key = 'input'
      if(use_input_mask .and. remove_badchan .and.
     -   nbad_chan.gt.0)masksrc_key = 'combined'

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
      out_amp_open = .false.
      out_ang_open = .false.
      out_mask_open = .false.
      out_nvalid_open = .false.
      in_mask_open = .false.
      use_input_mask = .false.
      masksrc_key = 'generated'
      if(nchar(mask_input_cube_file).gt.0)then
              use_input_mask = .true.
              masksrc_key = 'input'
      endif
      if(remove_badchan .and. use_input_mask)masksrc_key = 'combined'
      nanchk_key = 'on'
      nan_check_on = .true.
      if(index(mask_trust_mode,'strict').gt.0 .or.
     -   index(mask_trust_mode,'STRICT').gt.0)then
              nanchk_key = 'off'
              nan_check_on = .false.
      endif
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


      !  Create the new RM FITS files unless this is a dry-run.
      !  Also pre-check for output file collisions before calling FITSIO.
      if(.not.dry_run)then
              inquire(file=outfileAMP(1:nchar(outfileAMP)),
     -                exist=out_exists)
              if(out_exists)then
                      write(*,*)" "
                      write(*,*)"ERROR: Output file already exists:"
                      write(*,*)outfileAMP(1:nchar(outfileAMP))
                      write(*,*)"Refusing to overwrite existing file."
                      write(*,*)"Please remove/rename it and run again."
                      stop
              endif
              inquire(file=outfileANG(1:nchar(outfileANG)),
     -                exist=out_exists)
              if(out_exists)then
                      write(*,*)" "
                      write(*,*)"ERROR: Output file already exists:"
                      write(*,*)outfileANG(1:nchar(outfileANG))
                      write(*,*)"Refusing to overwrite existing file."
                      write(*,*)"Please remove/rename it and run again."
                      stop
              endif
              inquire(file=outfileMASK(1:nchar(outfileMASK)),
     -                exist=out_exists)
              if(out_exists)then
                      write(*,*)" "
                      write(*,*)"ERROR: Output file already exists:"
                      write(*,*)outfileMASK(1:nchar(outfileMASK))
                      write(*,*)"Refusing to overwrite existing file."
                      write(*,*)"Please remove/rename it and run again."
                      stop
              endif
              inquire(file=outfileNVALID(1:nchar(outfileNVALID)),
     -                exist=out_exists)
              if(out_exists)then
                      write(*,*)" "
                      write(*,*)"ERROR: Output file already exists:"
                      write(*,*)outfileNVALID(1:nchar(outfileNVALID))
                      write(*,*)"Refusing to overwrite existing file."
                      write(*,*)"Please remove/rename it and run again."
                      stop
              endif

              status = 0
              call ftinit(41,outfileAMP,blocksize,status)
              if(status.ne.0)then
                      write(*,*)"Error creating RM output file:"
                      write(*,*)outfileAMP(1:nchar(outfileAMP))
                      call printerror(status)
                      stop
              endif
              out_amp_open = .true.

              status = 0
              call ftinit(42,outfileANG,blocksize,status)
              if(status.ne.0)then
                      write(*,*)"Error creating PA output file:"
                      write(*,*)outfileANG(1:nchar(outfileANG))
                      call printerror(status)
                      stop
              endif
              out_ang_open = .true.

              status = 0
              call ftinit(43,outfileMASK,blocksize,status)
              if(status.ne.0)then
                      write(*,*)"Error creating MASK output file:"
                      write(*,*)outfileMASK(1:nchar(outfileMASK))
                      call printerror(status)
                      stop
              endif
              out_mask_open = .true.

              status = 0
              call ftinit(44,outfileNVALID,blocksize,status)
              if(status.ne.0)then
                      write(*,*)"Error creating NVALID output file:"
                      write(*,*)outfileNVALID(1:nchar(outfileNVALID))
                      call printerror(status)
                      stop
              endif
              out_nvalid_open = .true.
      endif


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

              ! Keep only RA/Dec/frequency varying in extraction.
              do i = 1,naxis
                 if(i.ne.1 .and. i.ne.2 .and. i.ne.freq_axis)then
                        fpixels(i) = 1
                        lpixels(i) = 1
                        incs(i) = 1
                 endif
              enddo
      else
              write(*,*)" "
              write(*,*)"Sub-section of Q and U-cubes will be used"
              write(*,*)"for the tomography... "
              write(*,*)" "

              ! Initialize all axes; this is required for naxis=4 cubes.
              do i = 1,naxis
                 fpixels(i) = 1
                 lpixels(i) = naxes(i)
                 incs(i) = 1
              enddo

              ! Keep only RA/Dec/frequency varying in extraction.
              do i = 1,naxis
                 if(i.ne.1 .and. i.ne.2 .and. i.ne.freq_axis)then
                         fpixels(i) = 1
                         lpixels(i) = 1
                         incs(i) = 1
                 endif
              enddo
              
              ! Use subimage parameters directly from config
              if (subim_ra_blc .eq. 0) then
                  fpixels(1) = 1
              else
                  fpixels(1) = subim_ra_blc
              endif
              if (subim_ra_trc .eq. 0) then
                  lpixels(1) = naxes(1)
              else
                  lpixels(1) = subim_ra_trc
              endif
              incs(1) = subim_ra_inc
              
              if (subim_dec_blc .eq. 0) then
                  fpixels(2) = 1
              else
                  fpixels(2) = subim_dec_blc
              endif
              if (subim_dec_trc .eq. 0) then
                  lpixels(2) = naxes(2)
              else
                  lpixels(2) = subim_dec_trc
              endif
              incs(2) = subim_dec_inc
              
              if (subim_chan_blc .eq. 0) then
                  fpixels(freq_axis) = 1
              else
                  fpixels(freq_axis) = subim_chan_blc
              endif
              if (subim_chan_trc .eq. 0) then
                  lpixels(freq_axis) = naxes(freq_axis)
              else
                  lpixels(freq_axis) = subim_chan_trc
              endif
              incs(freq_axis) = subim_chan_inc
              
              write(*,*)"Using subimage from config:"
              write(*,*)"RA: ",fpixels(1)," to ",lpixels(1),
     -          " step ",incs(1)
              write(*,*)"Dec: ",fpixels(2)," to ",lpixels(2),
     -          " step ",incs(2)
              write(*,*)"Chan(axis",freq_axis,"): ",
     -          fpixels(freq_axis)," to ",lpixels(freq_axis),
     -          " step ",incs(freq_axis)
              
              ! Validate subimage bounds
              do i = 1,naxis
                 if(lpixels(i).lt.fpixels(i))then
                      write(*,*)" "
                      write(*,*)"Error: In config: "
                      write(*,*)"last-pix < first-pix"
                      write(*,*)"In axis number:",i
                      write(*,*)"Quitting now..."
                      write(*,*)" "
                      call FTCLOS(21,status)
                      call FTCLOS(22,status)
                      call FTCLOS(41,status)
                      call FTCLOS(42,status)
                      stop 
                 endif
                 if(lpixels(i).gt.naxes(i))then
                      write(*,*)" "
                      write(*,*)"Error: In config: "
                      write(*,*)"last-pix > ",naxes(i)
                      write(*,*)"Output image exceeds max dimension"
                      write(*,*)"In axis number:",i
                      write(*,*)"Quitting now..."
                      write(*,*)" "
                      call FTCLOS(21,status)
                      call FTCLOS(22,status)
                      call FTCLOS(41,status)
                      call FTCLOS(42,status)
                      stop
                 endif
              enddo
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

      zpix_beg = fpixels(freq_axis)
      zpix_end = lpixels(freq_axis)
      nz_out = int((zpix_end - zpix_beg)/incs(freq_axis)) + 1

      ntot_out = nx_out*ny_out*nz_out

      ! Allocate per-pixel spectrum work buffers sized to actual nz_out
      allocate(data_arrI(nz_out))
      allocate(data_arrQ(nz_out))
      allocate(data_arrU(nz_out))
      allocate(L_sq(nz_out))
      allocate(Q_now(nz_out))
      allocate(U_now(nz_out))

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
      do i = zpix_end,zpix_beg,-incs(freq_axis)
         if(flag_arr(i).eq.1)then
                 ngood_chan = ngood_chan + 1
                 L_sq(ngood_chan) = (conv_fac/zval(i))**2
         endif
      enddo
      ! Use explicit flag to select RM extraction mode
      if (use_auto_rm_range .eq. 1) then
           nrm_out_par = ngood_chan
           nrm_out = nrm_out_par*ofac
      else
           nrm_out = nrm_out_par*ofac
      endif

      ! Allocate RM arrays sized to actual nrm_out and ngood_chan
      allocate(RM(nrm_out))
      allocate(p_ex(nrm_out))
      allocate(phi_ex(nrm_out))
      allocate(cos_arr(nrm_out, ngood_chan))
      allocate(sin_arr(nrm_out, ngood_chan))

      call extract_general_setup(L_sq,ngood_chan,fac,beg_rm,end_rm,
     -  nrm_out,RM,cos_arr,sin_arr,nrm_out,ngood_chan,
     -  use_auto_rm_range,ofac)
      dRM = (RM(nrm_out) - RM(1))/real(nrm_out - 1)
      open(77,file='sampled_RM.txt',status='unknown')
      write(77,*)"# ofac: ",ofac
      write(77,*)"# ngood_chan: ",ngood_chan
      write(77,*)"# nrm_out: ",nrm_out
      write(77,*)"# fac: ",fac
      write(77,*)"# beg_rm: ",beg_rm
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
      do i = zpix_end,zpix_beg,-incs(freq_axis)
         atmp = (conv_fac/zval(i))**2
         write(79,*)zval(i),"    ",atmp,"   ",flag_arr(i) 
      enddo
      close(79)

      !----------------------------------------------------
      ! Tile planning for memory-efficient cube processing.
      mem_avail_kb = 0_int64
      mem_unit = 91
      open(mem_unit,file='/proc/meminfo',status='old',iostat=ios_mem)
      if(ios_mem.eq.0)then
              do
                      read(mem_unit,'(A)',iostat=ios_mem) mem_line
                      if(ios_mem.ne.0)exit
                      if(index(mem_line,'MemAvailable:').eq.1)then
                              read(mem_line(14:),*,
     -                             iostat=ios_mem) mem_kb_tmp
                              if(ios_mem.eq.0)mem_avail_kb = mem_kb_tmp
                              exit
                      endif
              enddo
              close(mem_unit)
      endif
      if(mem_avail_kb.le.0_int64)then
              mem_avail_kb = 4194304_int64
      endif

      in_fields = 2
      if(need_icube)in_fields = 3

      bytes_per_tile_pixel = int(4,kind=int64)*
     -      (int(in_fields,kind=int64)*int(nz_out,kind=int64) +
     -       int(2*nrm_out,kind=int64))
        mem_safe_bytes = int(tile_mem_frac *
     -      real(mem_avail_kb,kind=dp) * 1024.0_dp,
     -      kind=int64)
      if(mem_safe_bytes.le.bytes_per_tile_pixel)then
              mem_safe_bytes = bytes_per_tile_pixel
      endif
      tile_pixels_max = mem_safe_bytes / bytes_per_tile_pixel
      if(tile_pixels_max.lt.1_int64)tile_pixels_max = 1_int64
      image_pixels_total = nx_out
      image_pixels_total = image_pixels_total * ny_out

      if(tile_auto .or. tile_ra.le.0 .or. tile_dec.le.0)then
              if(tile_pixels_max.ge.image_pixels_total)then
                      tile_ra = nx_out
                      tile_dec = ny_out
              else
                      tile_ra = min(nx_out,
     -                   max(1,int(sqrt(real(tile_pixels_max,
     -                   kind=dp)*real(nx_out,kind=dp)/
     -                   real(ny_out,kind=dp)))))
                      if(tile_ra.lt.1)tile_ra = 1
                      tile_dec = min(ny_out,
     -                   max(1,int(tile_pixels_max/
     -                   int(tile_ra,kind=int64))))
                      if(tile_dec.lt.1)tile_dec = 1
              endif
      else
              tile_ra = max(1,min(tile_ra,nx_out))
              tile_dec = max(1,min(tile_dec,ny_out))
      endif

      tile_bytes_est = int(tile_ra,kind=int64) *
     -                 int(tile_dec,kind=int64) * bytes_per_tile_pixel
      do while(tile_bytes_est.gt.mem_safe_bytes .and.
     -         (tile_ra.gt.1 .or. tile_dec.gt.1))
              if(tile_ra.ge.tile_dec .and. tile_ra.gt.1)then
                      tile_ra = max(1,tile_ra/2)
              else if(tile_dec.gt.1)then
                      tile_dec = max(1,tile_dec/2)
              endif
              tile_bytes_est = int(tile_ra,kind=int64) *
     -                 int(tile_dec,kind=int64) * bytes_per_tile_pixel
      enddo

      write(*,*)" "
      write(*,*)"Tile planner (Phase-2):"
      write(*,*)" MemAvailable(kB): ",mem_avail_kb
      write(*,*)" tile_mem_frac: ",tile_mem_frac
      write(*,*)" tile_ra x tile_dec (output px): ",tile_ra,tile_dec
      write(*,*)" Estimated tile memory (MB): ",
     -           real(tile_bytes_est,kind=dp)/(1024.0_dp*1024.0_dp)

      if(dry_run)then
              open(96,file='tile_autotune.cfg',status='unknown')
              write(96,*)"# Autogenerated tile hints"
              write(96,*)"# Copy these KEY=VALUE lines to your cfg"
              write(96,*)"tile_auto=n"
              write(96,*)"tile_ra=",tile_ra
              write(96,*)"tile_dec=",tile_dec
              write(96,*)"tile_mem_frac=",tile_mem_frac
              write(96,*)"# Suggested subimage chunk for one pass"
              write(96,*)"subim_ra_blc=",xpix_beg
              write(96,*)"subim_ra_trc=",min(xpix_end,
     -             xpix_beg + (tile_ra-1)*incs(1))
              write(96,*)"subim_dec_blc=",ypix_beg
              write(96,*)"subim_dec_trc=",min(ypix_end,
     -             ypix_beg + (tile_dec-1)*incs(2))
              close(96)
              call write_runtime_estimate('runtime_estimate.txt',
     -        image_pixels_total,
     -             nz_totpix,ngood_chan,nbad_chan,nrm_out,output_mode,
     -             tile_ra,tile_dec,nx_out,ny_out,tile_bytes_est,
     -             tile_mem_frac,status)
              write(*,*)"Dry-run mode enabled. Wrote tile_autotune.cfg"
              write(*,*)"Dry-run mode enabled."
              write(*,*)"Wrote runtime_estimate.txt"
              write(*,*)"No tomography executed in dry-run mode."
              goto 9999
      endif

      allocate(specQ(tile_ra*tile_dec*nz_out))
      allocate(specU(tile_ra*tile_dec*nz_out))
      if(use_input_mask)allocate(specMask(tile_ra*tile_dec*nz_out))
      if(need_icube)allocate(specI(tile_ra*tile_dec*nz_out))
      allocate(p_tile_arr(tile_ra*tile_dec*nrm_out))
      allocate(phi_tile_arr(tile_ra*tile_dec*nrm_out))
      allocate(mask_tile_arr(tile_ra*tile_dec*nz_out))
      allocate(nvalid_tile_arr(tile_ra*tile_dec))
      allocate(wts_now(ngood_chan))


      ! Irrespective of the total number of output pixels, 
      ! we will read the spectra in the cube on a pix-by-pix 
      ! basis. That way, the variable array named "spec" 
      ! need only be defined to have dimension maxchan.

      write(*,*)"xpix-beg,xpix-end,inc: ",xpix_beg,xpix_end,incs(1)
      write(*,*)"ypix-beg,ypix-end,inc: ",ypix_beg,ypix_end,incs(2)
      write(*,*)"zpix-beg,zpix-end,inc: ",
     -          zpix_beg,zpix_end,incs(freq_axis)

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
      naxes_mask(1) = nx_out
      naxes_mask(2) = ny_out
      naxes_mask(3) = nz_out
      naxes_nvalid(1) = nx_out
      naxes_nvalid(2) = ny_out
      call ftphpr(41,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(42,simple,bitpix,3,naxes_out,0,1,extend,status)
      call ftphpr(43,simple,8,3,naxes_mask,0,1,extend,status)
      call ftphpr(44,simple,16,2,naxes_nvalid,0,1,extend,status)

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
      ! WCS header for output RM cubes (FITS Paper I/II standard).
      ! Axes 1+2: CRVAL unchanged (passthrough), CRPIX offset by
      !           subimage start, CDELT scaled by stride only.
      !           No sec(delta) - SIN/TAN projection handles geometry.
      ! Axis  3: RM synthesised axis; CTYPE and CUNIT set explicitly.
      ! Frame:   RADESYS/EQUINOX preferred; EPOCH as legacy fallback.
      ! Rotation: PC matrix elements passthrough if present in input.
      decimals = 13
      status = 0

      ! --- Axis 1 (RA): CTYPE passthrough ---
      call ftgkys(21,'ctype1',ctype,comment,status)
      call ftpkys(41,'ctype1',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(42,'ctype1',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(43,'ctype1',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(44,'ctype1',ctype(1:nchar(ctype)),' ',status)
      status = 0

      ! --- Axis 1: CRVAL passthrough, CRPIX offset, CDELT scaled ---
      call ftgkyd(21,'crval1',cval,comment,status)
      call ftpkyd(41,'crval1',cval,decimals,' ',status)
      call ftpkyd(42,'crval1',cval,decimals,' ',status)
      call ftpkyd(43,'crval1',cval,decimals,' ',status)
      call ftpkyd(44,'crval1',cval,decimals,' ',status)
      call ftgkyd(21,'crpix1',atmp8,comment,status)
      atmp8 = (atmp8 - dble(xpix_beg)) / dble(incs(1)) + 1.0d0
      call ftpkyd(41,'crpix1',atmp8,decimals,' ',status)
      call ftpkyd(42,'crpix1',atmp8,decimals,' ',status)
      call ftpkyd(43,'crpix1',atmp8,decimals,' ',status)
      call ftpkyd(44,'crpix1',atmp8,decimals,' ',status)
      call ftgkyd(21,'cdelt1',cdelt,comment,status)
      cdelt = dble(incs(1)) * cdelt
      call ftpkyd(41,'cdelt1',cdelt,decimals,' ',status)
      call ftpkyd(42,'cdelt1',cdelt,decimals,' ',status)
      call ftpkyd(43,'cdelt1',cdelt,decimals,' ',status)
      call ftpkyd(44,'cdelt1',cdelt,decimals,' ',status)
      status = 0

      ! --- Axis 1: CUNIT passthrough if present ---
      call ftgkys(21,'cunit1',ctype,comment,status)
      if(status.eq.0)then
              call ftpkys(41,'cunit1',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(42,'cunit1',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(43,'cunit1',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(44,'cunit1',ctype(1:nchar(ctype)),' ',status)
      endif
      status = 0

      ! --- Axis 2 (Dec): CTYPE passthrough ---
      call ftgkys(21,'ctype2',ctype,comment,status)
      call ftpkys(41,'ctype2',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(42,'ctype2',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(43,'ctype2',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(44,'ctype2',ctype(1:nchar(ctype)),' ',status)
      status = 0

      ! --- Axis 2: CRVAL passthrough, CRPIX offset, CDELT scaled ---
      call ftgkyd(21,'crval2',cval,comment,status)
      call ftpkyd(41,'crval2',cval,decimals,' ',status)
      call ftpkyd(42,'crval2',cval,decimals,' ',status)
      call ftpkyd(43,'crval2',cval,decimals,' ',status)
      call ftpkyd(44,'crval2',cval,decimals,' ',status)
      call ftgkyd(21,'crpix2',atmp8,comment,status)
      atmp8 = (atmp8 - dble(ypix_beg)) / dble(incs(2)) + 1.0d0
      call ftpkyd(41,'crpix2',atmp8,decimals,' ',status)
      call ftpkyd(42,'crpix2',atmp8,decimals,' ',status)
      call ftpkyd(43,'crpix2',atmp8,decimals,' ',status)
      call ftpkyd(44,'crpix2',atmp8,decimals,' ',status)
      call ftgkyd(21,'cdelt2',cdelt,comment,status)
      cdelt = dble(incs(2)) * cdelt
      call ftpkyd(41,'cdelt2',cdelt,decimals,' ',status)
      call ftpkyd(42,'cdelt2',cdelt,decimals,' ',status)
      call ftpkyd(43,'cdelt2',cdelt,decimals,' ',status)
      call ftpkyd(44,'cdelt2',cdelt,decimals,' ',status)
      status = 0

      ! --- Axis 2: CUNIT passthrough if present ---
      call ftgkys(21,'cunit2',ctype,comment,status)
      if(status.eq.0)then
              call ftpkys(41,'cunit2',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(42,'cunit2',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(43,'cunit2',ctype(1:nchar(ctype)),' ',status)
              call ftpkys(44,'cunit2',ctype(1:nchar(ctype)),' ',status)
      endif
      status = 0

      ! --- Axis 3: RM synthesised axis ---
      call ftpkys(41,'ctype3','FDEP','Faraday depth',status)
      call ftpkys(42,'ctype3','FDEP','Faraday depth',status)
      call ftpkys(41,'cunit3','rad/m**2','RM axis units',status)
      call ftpkys(42,'cunit3','rad/m**2','RM axis units',status)
      call ftpkyd(41,'crval3',dble(RM(1)),decimals,
     -            'Reference RM (rad/m^2)',status)
      call ftpkyd(42,'crval3',dble(RM(1)),decimals,
     -            'Reference RM (rad/m^2)',status)
      call ftpkyd(41,'crpix3',1.0d0,decimals,'Reference pixel',status)
      call ftpkyd(42,'crpix3',1.0d0,decimals,'Reference pixel',status)
      call ftpkyd(41,'cdelt3',dble(dRM),decimals,'RM spacing',status)
      call ftpkyd(42,'cdelt3',dble(dRM),decimals,'RM spacing',status)
      call ftpkys(43,'ctype3','FREQ','Frequency axis',status)
      call ftpkys(43,'cunit3','Hz','Frequency axis units',status)
      call ftpkyd(43,'crval3',dble(zval(zpix_beg)),decimals,
     -            'Reference frequency',status)
      call ftpkyd(43,'crpix3',1.0d0,decimals,'Reference pixel',status)
      call ftpkyd(43,'cdelt3',dble(incs(freq_axis))*dble(zinc_im),
     -            decimals,'Frequency spacing',status)
      status = 0

      ! --- PC rotation matrix: passthrough if present in input ---
      call ftgkyd(21,'pc1_1',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'pc1_1',cval,decimals,' ',status)
              call ftpkyd(42,'pc1_1',cval,decimals,' ',status)
              call ftpkyd(43,'pc1_1',cval,decimals,' ',status)
              call ftpkyd(44,'pc1_1',cval,decimals,' ',status)
      endif
      status = 0
      call ftgkyd(21,'pc1_2',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'pc1_2',cval,decimals,' ',status)
              call ftpkyd(42,'pc1_2',cval,decimals,' ',status)
              call ftpkyd(43,'pc1_2',cval,decimals,' ',status)
              call ftpkyd(44,'pc1_2',cval,decimals,' ',status)
      endif
      status = 0
      call ftgkyd(21,'pc2_1',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'pc2_1',cval,decimals,' ',status)
              call ftpkyd(42,'pc2_1',cval,decimals,' ',status)
              call ftpkyd(43,'pc2_1',cval,decimals,' ',status)
              call ftpkyd(44,'pc2_1',cval,decimals,' ',status)
      endif
      status = 0
      call ftgkyd(21,'pc2_2',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'pc2_2',cval,decimals,' ',status)
              call ftpkyd(42,'pc2_2',cval,decimals,' ',status)
              call ftpkyd(43,'pc2_2',cval,decimals,' ',status)
              call ftpkyd(44,'pc2_2',cval,decimals,' ',status)
      endif
      status = 0

      ! --- Coordinate frame: RADESYS/EQUINOX preferred, EPOCH fallback ---
      call ftgkys(21,'radesys',ctype,comment,status)
      if(status.eq.0)then
              call ftpkys(41,'radesys',ctype(1:nchar(ctype)),
     -                    ' ',status)
              call ftpkys(42,'radesys',ctype(1:nchar(ctype)),
     -                    ' ',status)
              call ftpkys(43,'radesys',ctype(1:nchar(ctype)),
     -                    ' ',status)
              call ftpkys(44,'radesys',ctype(1:nchar(ctype)),
     -                    ' ',status)
      endif
      status = 0
      call ftgkyd(21,'equinox',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'equinox',cval,decimals,' ',status)
              call ftpkyd(42,'equinox',cval,decimals,' ',status)
              call ftpkyd(43,'equinox',cval,decimals,' ',status)
              call ftpkyd(44,'equinox',cval,decimals,' ',status)
      else
              status = 0
              call ftgkyd(21,'epoch',cval,comment,status)
              if(status.eq.0)then
                      call ftpkyd(41,'epoch',cval,decimals,' ',status)
                      call ftpkyd(42,'epoch',cval,decimals,' ',status)
                      call ftpkyd(43,'epoch',cval,decimals,' ',status)
                      call ftpkyd(44,'epoch',cval,decimals,' ',status)
              else
                      write(*,*)'WCS: no EQUINOX/EPOCH; default J2000'
                      call ftpkyd(41,'equinox',2000.0d0,decimals,
     -                             'Coord equinox',status)
                      call ftpkyd(42,'equinox',2000.0d0,decimals,
     -                             'Coord equinox',status)
                      call ftpkyd(43,'equinox',2000.0d0,decimals,
     -                             'Coord equinox',status)
                      call ftpkyd(44,'equinox',2000.0d0,decimals,
     -                             'Coord equinox',status)
              endif
      endif

      ! --- LONPOLE/LATPOLE: passthrough if present ---
      status = 0
      call ftgkyd(21,'lonpole',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'lonpole',cval,decimals,' ',status)
              call ftpkyd(42,'lonpole',cval,decimals,' ',status)
              call ftpkyd(43,'lonpole',cval,decimals,' ',status)
              call ftpkyd(44,'lonpole',cval,decimals,' ',status)
      endif
      status = 0
      call ftgkyd(21,'latpole',cval,comment,status)
      if(status.eq.0)then
              call ftpkyd(41,'latpole',cval,decimals,' ',status)
              call ftpkyd(42,'latpole',cval,decimals,' ',status)
              call ftpkyd(43,'latpole',cval,decimals,' ',status)
              call ftpkyd(44,'latpole',cval,decimals,' ',status)
      endif

      ! --- BUNIT: passthrough for cube 1 (amp/re); set for cube 2 ---
      status = 0
      call ftgkys(21,'bunit',ctype,comment,status)
      if(status.ne.0)then
              ctype = 'UNKNOWN'
              status = 0
      endif
      call ftpkys(41,'bunit',ctype(1:nchar(ctype)),
     -            'Pixel data units',status)
      call ftpkys(42,'bunit','rad',
     -            'Pixel data units (angle)',status)
      call ftpkys(43,'bunit','FLAG',
     -            'Mask value: 0 bad, 1 good',status)
      call ftpkys(44,'bunit','COUNT',
     -            'Number of valid channels',status)

      ! --- Metadata: OBJECT, OBSERVER, TELESCOP ---
      status = 0
      call ftgkys(21,'object',ctype,comment,status)
      if(status.ne.0)then
              ctype = 'UNKNOWN'
              status = 0
      endif
      call ftpkys(41,'object',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(42,'object',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(43,'object',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(44,'object',ctype(1:nchar(ctype)),' ',status)
      status = 0
      call ftgkys(21,'observer',ctype,comment,status)
      if(status.ne.0)then
              ctype = 'UNKNOWN'
              status = 0
      endif
      call ftpkys(41,'observer',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(42,'observer',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(43,'observer',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(44,'observer',ctype(1:nchar(ctype)),' ',status)
      status = 0
      call ftgkys(21,'telescop',ctype,comment,status)
      if(status.ne.0)then
              ctype = 'UNKNOWN'
              status = 0
      endif
      call ftpkys(41,'telescop',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(42,'telescop',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(43,'telescop',ctype(1:nchar(ctype)),' ',status)
      call ftpkys(44,'telescop',ctype(1:nchar(ctype)),' ',status)
      status = 0

      call ftpkys(41,'MASKSRC',masksrc_key(1:nchar(masksrc_key)),
     -            'Mask source: generated/input/combined',status)
      call ftpkys(42,'MASKSRC',masksrc_key(1:nchar(masksrc_key)),
     -            'Mask source: generated/input/combined',status)
      call ftpkys(43,'MASKSRC',masksrc_key(1:nchar(masksrc_key)),
     -            'Mask source: generated/input/combined',status)
      call ftpkys(44,'MASKSRC',masksrc_key(1:nchar(masksrc_key)),
     -            'Mask source: generated/input/combined',status)
      call ftpkyj(41,'NBADGLOB',nbad_chan,
     -            'No. of globally bad channels',status)
      call ftpkyj(42,'NBADGLOB',nbad_chan,
     -            'No. of globally bad channels',status)
      call ftpkyj(43,'NBADGLOB',nbad_chan,
     -            'No. of globally bad channels',status)
      call ftpkyj(44,'NBADGLOB',nbad_chan,
     -            'No. of globally bad channels',status)
      call ftpkys(41,'NANCHK',nanchk_key(1:nchar(nanchk_key)),
     -            'NaN validity check on/off',status)
      call ftpkys(42,'NANCHK',nanchk_key(1:nchar(nanchk_key)),
     -            'NaN validity check on/off',status)
      call ftpkys(43,'NANCHK',nanchk_key(1:nchar(nanchk_key)),
     -            'NaN validity check on/off',status)
      call ftpkys(44,'NANCHK',nanchk_key(1:nchar(nanchk_key)),
     -            'NaN validity check on/off',status)
      call ftpkys(41,'MASKTRUS',
     -            mask_trust_mode(1:nchar(mask_trust_mode)),
     -            'Mask trust mode: safe/strict',status)
      call ftpkys(42,'MASKTRUS',
     -            mask_trust_mode(1:nchar(mask_trust_mode)),
     -            'Mask trust mode: safe/strict',status)
      call ftpkys(43,'MASKTRUS',
     -            mask_trust_mode(1:nchar(mask_trust_mode)),
     -            'Mask trust mode: safe/strict',status)
      call ftpkys(44,'MASKTRUS',
     -            mask_trust_mode(1:nchar(mask_trust_mode)),
     -            'Mask trust mode: safe/strict',status)
        status = 0




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
      if(need_icube)then
              call FTOPEN(40,infileI,rwmode,blocksize,status)
      endif
      if(use_input_mask)then
              status = 0
              call FTOPEN(45,mask_input_cube_file,
     -                    rwmode,blocksize,status)
              if(status.ne.0)then
                      write(*,*)"Error opening input mask cube:"
                      write(*,*)mask_input_cube_file(
     -                         1:nchar(mask_input_cube_file))
                      call printerror(status)
                      stop
              endif
              in_mask_open = .true.
      endif

      if(line_cut)then
              open(121,file='rm_spec.txt')
      endif

      tmp_cnt1 = 0
      tmp_cnt2 = 0
      cnt1 = 0
      progress_total = nx_out*ny_out
      progress_step = max(1, progress_total/10)
      progress_next_pct = 10
      progress_next_count = progress_step
      do ix_tile_beg = xpix_beg,xpix_end,tile_ra*incs(1)
        ix_tile_end = min(xpix_end,
     -                     ix_tile_beg + (tile_ra-1)*incs(1))
        nx_tile = int((ix_tile_end - ix_tile_beg)/incs(1)) + 1

        do iy_tile_beg = ypix_beg,ypix_end,tile_dec*incs(2)
            iy_tile_end = min(ypix_end,
     -                        iy_tile_beg + (tile_dec-1)*incs(2))
            ny_tile = int((iy_tile_end - iy_tile_beg)/incs(2)) + 1

            write(*,*)"Doing tile x:[",ix_tile_beg,",",ix_tile_end,
     -                "] y:[",iy_tile_beg,",",iy_tile_end,"]"

            fpixels(1) = ix_tile_beg
            lpixels(1) = ix_tile_end
            fpixels(2) = iy_tile_beg
            lpixels(2) = iy_tile_end
            fpixels(freq_axis) = zpix_beg
            lpixels(freq_axis) = zpix_end

            call FTGSVE(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                   nullval,specQ,anyflg,status)
            call FTGSVE(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                   nullval,specU,anyflg,status)
             if(use_input_mask)then
                      call FTGSVE(45,group,naxis,naxes,fpixels,lpixels,
     -                   incs,nullval,specMask,anyflg,status)
             endif
            if(need_icube)then
                    call FTGSVE(40,group,naxis,naxes,fpixels,lpixels,
     -                   incs,nullval,specI,anyflg,status)
            endif

            do iy_loc = 1,ny_tile
               iy = iy_tile_beg + (iy_loc-1)*incs(2)
               do ix_loc = 1,nx_tile
                  ix = ix_tile_beg + (ix_loc-1)*incs(1)
                  cnt1 = cnt1 + 1

                  do iz = 1,nz_out
                     i = iz
                     tmp_index = ix_loc + (iy_loc-1)*nx_tile +
     -                          (iz-1)*nx_tile*ny_tile
                     data_arrQ(i) = specQ(tmp_index)
                     data_arrU(i) = specU(tmp_index)
                     if(need_icube)then
                             data_arrI(i) = specI(tmp_index)
                     endif
                  enddo

                  ngood_chan = 0
                  cnt_good = 0
                  cnt2 = nz_out + 1
                  if(.not.remove_QU_bias)then
                          do i = zpix_end,zpix_beg,
     -                         -incs(freq_axis)
                             cnt2 = cnt2 - 1
                             chan_valid = (flag_arr(i).eq.1)
                             if(use_input_mask)then
                                     tmp_index = ix_loc +
     -                                  (iy_loc-1)*nx_tile +
     -                                  (cnt2-1)*nx_tile*ny_tile
                                     mask_val = specMask(tmp_index)
                                     if(mask_val.le.0.5)then
                                             chan_valid = .false.
                                     endif
                             endif
                             if(nan_check_on)then
                                     if(data_arrQ(cnt2).ne.
     -                                  data_arrQ(cnt2))then
                                             chan_valid = .false.
                                     endif
                                     if(data_arrU(cnt2).ne.
     -                                  data_arrU(cnt2))then
                                             chan_valid = .false.
                                     endif
                             endif
                             if(chan_valid)then
                                     ngood_chan = ngood_chan + 1
                                     cnt_good = cnt_good + 1
                                     Q_now(ngood_chan) = data_arrQ(cnt2)
                                     U_now(ngood_chan) = data_arrU(cnt2)
                                     wts_now(ngood_chan) = 1.0
                             endif
                             tmp_index = ix_loc + (iy_loc-1)*nx_tile +
     -                                (cnt2-1)*nx_tile*ny_tile
                             if(chan_valid)then
                                     mask_tile_arr(tmp_index) = 1
                             else
                                     mask_tile_arr(tmp_index) = 0
                             endif
                          enddo
                  else
                          do i = zpix_end,zpix_beg,
     -                         -incs(freq_axis)
                             cnt2 = cnt2 - 1
                             chan_valid = (flag_arr(i).eq.1)
                             if(use_input_mask)then
                                     tmp_index = ix_loc +
     -                                  (iy_loc-1)*nx_tile +
     -                                  (cnt2-1)*nx_tile*ny_tile
                                     mask_val = specMask(tmp_index)
                                     if(mask_val.le.0.5)then
                                             chan_valid = .false.
                                     endif
                             endif
                             if(nan_check_on)then
                                     if(data_arrQ(cnt2).ne.
     -                                  data_arrQ(cnt2))then
                                             chan_valid = .false.
                                     endif
                                     if(data_arrU(cnt2).ne.
     -                                  data_arrU(cnt2))then
                                             chan_valid = .false.
                                     endif
                             endif
                             if(chan_valid)then
                                     ngood_chan = ngood_chan + 1
                                     cnt_good = cnt_good + 1
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
                                     Q_now(ngood_chan) =
     -                                  data_arrQ(cnt2) -
     -                                  (data_arrI(cnt2)*slopeQ + resiQ)
                                     U_now(ngood_chan) =
     -                                  data_arrU(cnt2) -
     -                                  (data_arrI(cnt2)*slopeU + resiU)
                                     wts_now(ngood_chan) = 1.0
                             endif
                             tmp_index = ix_loc + (iy_loc-1)*nx_tile +
     -                                (cnt2-1)*nx_tile*ny_tile
                             if(chan_valid)then
                                     mask_tile_arr(tmp_index) = 1
                             else
                                     mask_tile_arr(tmp_index) = 0
                             endif
                          enddo
                  endif

                  nvalid_pix = cnt_good
                  tmp_index = ix_loc + (iy_loc-1)*nx_tile
                  nvalid_tile_arr(tmp_index) = nvalid_pix

                  if(ngood_chan.le.0)then
                          do i = 1,nrm_out
                                  p_ex(i) = 0.0
                                  phi_ex(i) = 0.0
                          enddo
                  else
                          if (output_mode .eq. 1) then
                                  call extract_general_ri_w(
     -                               Q_now,U_now,wts_now,ngood_chan,
     -                               nrm_out,p_ex,phi_ex,
     -                               cos_arr,sin_arr,nrm_out,ngood_chan,
     -                               rem_mean)
                          else
                                  call extract_general_w(
     -                               Q_now,U_now,wts_now,ngood_chan,
     -                               nrm_out,p_ex,phi_ex,
     -                               cos_arr,sin_arr,nrm_out,ngood_chan,
     -                               rem_mean)
                                  if (ap_angle_mode .eq. 1) then
                                          do i = 1,nrm_out
                                                  phi_ex(i) =
     -                                          0.5*phi_ex(i)
                                          enddo
                                  endif
                          endif
                  endif

                  do i = 1,nrm_out
                     tmp_index = ix_loc + (iy_loc-1)*nx_tile +
     -                          (i-1)*nx_tile*ny_tile
                     p_tile_arr(tmp_index) = p_ex(i)
                     phi_tile_arr(tmp_index) = phi_ex(i)
                  enddo

                  if(line_cut)then
                          tmp_cnt1 = tmp_cnt1 + 1
                          write(16,rec=tmp_cnt1)
     -                    (Q_now(i),i=ngood_chan,1,-1)
                          tmp_cnt1 = tmp_cnt1 + 1
                          write(16,rec=tmp_cnt1)
     -                    (U_now(i),i=ngood_chan,1,-1)

                          tmp_cnt2 = tmp_cnt2 + 1
                          write(17,rec=tmp_cnt2)(p_ex(i),i=1,nrm_out)
                          tmp_cnt2 = tmp_cnt2 + 1
                          write(17,rec=tmp_cnt2)(phi_ex(i),i=1,nrm_out)
                          write(121,*)"## ix, iy: ",ix,iy
                          do i = 1,nrm_out
                                  write(121,*)p_ex(i), phi_ex(i)
                          enddo
                  endif

                  if(progress_total.gt.0)then
                        do while(cnt1.ge.progress_next_count .and.
     -                     progress_next_pct.le.100)
                           write(*,*)'Progress: ',progress_next_pct,
     -                        '% (',cnt1,' out of ',progress_total,')'
                           progress_next_pct = progress_next_pct + 10
                           progress_next_count = progress_next_count +
     -                        progress_step
                        enddo
                endif
               enddo
            enddo

            ix_out_beg = int((ix_tile_beg - xpix_beg)/incs(1)) + 1
            ix_out_end = ix_out_beg + nx_tile - 1
            iy_out_beg = int((iy_tile_beg - ypix_beg)/incs(2)) + 1
            iy_out_end = iy_out_beg + ny_tile - 1

            fpixels_out(1) = ix_out_beg
            lpixels_out(1) = ix_out_end
            fpixels_out(2) = iy_out_beg
            lpixels_out(2) = iy_out_end
            fpixels_out(3) = 1
            lpixels_out(3) = nrm_out

            call ftpsse(41,group,3,naxes_out,fpixels_out,lpixels_out,
     -                  p_tile_arr,status)
            if(status.gt.0)then
                    call printerror(status)
            endif
            call ftpsse(42,group,3,naxes_out,fpixels_out,lpixels_out,
     -                  phi_tile_arr,status)
            if(status.gt.0)then
                    call printerror(status)
            endif

            fpixels_out(3) = 1
            lpixels_out(3) = nz_out
            call ftpssb(43,group,3,naxes_mask,fpixels_out,lpixels_out,
     -                  mask_tile_arr,status)
            if(status.gt.0)then
                    call printerror(status)
            endif

            fpixels_nvalid(1) = ix_out_beg
            lpixels_nvalid(1) = ix_out_end
            fpixels_nvalid(2) = iy_out_beg
            lpixels_nvalid(2) = iy_out_end
            call ftpssi(44,group,2,naxes_nvalid,fpixels_nvalid,
     -                  lpixels_nvalid,nvalid_tile_arr,status)
            if(status.gt.0)then
                    call printerror(status)
            endif
        enddo
      enddo
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
      if(need_icube)then
              call FTCLOS(40,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing I-file"
                      call printerror(status)
              else
                      write(*,*)"Successfully read and closed "
                      write(*,*)"FITS Icube..."
              endif
      endif
      if(in_mask_open)then
              call FTCLOS(45,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing input MASK-file"
                      call printerror(status)
              endif
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


      ! Deallocate all dynamically allocated arrays
      if(allocated(xval)) deallocate(xval)
      if(allocated(yval)) deallocate(yval)
      if(allocated(zval)) deallocate(zval)
      if(allocated(flag_arr)) deallocate(flag_arr)
      if(allocated(data_arrI)) deallocate(data_arrI)
      if(allocated(data_arrQ)) deallocate(data_arrQ)
      if(allocated(data_arrU)) deallocate(data_arrU)
      if(allocated(L_sq)) deallocate(L_sq)
      if(allocated(Q_now)) deallocate(Q_now)
      if(allocated(U_now)) deallocate(U_now)
      if(allocated(RM)) deallocate(RM)
      if(allocated(p_ex)) deallocate(p_ex)
      if(allocated(phi_ex)) deallocate(phi_ex)
      if(allocated(cos_arr)) deallocate(cos_arr)
      if(allocated(sin_arr)) deallocate(sin_arr)
      if(allocated(specQ)) deallocate(specQ)
      if(allocated(specU)) deallocate(specU)
      if(allocated(specMask)) deallocate(specMask)
      if(allocated(specI)) deallocate(specI)
      if(allocated(p_tile_arr)) deallocate(p_tile_arr)
      if(allocated(phi_tile_arr)) deallocate(phi_tile_arr)
      if(allocated(mask_tile_arr)) deallocate(mask_tile_arr)
      if(allocated(nvalid_tile_arr)) deallocate(nvalid_tile_arr)
      if(allocated(wts_now)) deallocate(wts_now)

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
      if(out_amp_open)then
              status = 0
              call FTCLOS(41,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing RM-file"
                      call printerror(status)
              endif
      endif
      if(out_ang_open)then
              status = 0
              call FTCLOS(42,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing PA-file"
                      call printerror(status)
              endif
      endif
      if(out_mask_open)then
              status = 0
              call FTCLOS(43,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing MASK-file"
                      call printerror(status)
              endif
      endif
      if(out_nvalid_open)then
              status = 0
              call FTCLOS(44,status)
              if (status .gt. 0)then
                      write(*,*)"Problem closing NVALID-file"
                      call printerror(status)
              endif
      endif
!      endif

      ! -----------------------------------------------------------------


      end

      ! Modern Fortran approach: extraction routines, linspace, and nchar 
      ! are now in rm_synthesis_mod module
      ! Remaining utility subroutines included below:
      include 'myfits_info.f'
      include 'printerror.f'



