
      ! This code does the tomography of an image or a rectangular 
      ! subset of it given as inputs the FITS Q and the U spectral 
      ! cubes.
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


      ! TODO: KEYWORDS for the output FITS files in cases 
      !       when only a subimage is required, has to be 
      !       appropriately inserted... Currently the keywords 
      !       are read from the INPUT files and copied to the 
      !       output files. Any mismatches encountered thus (for
      !       example, if the output image does not contain 
      !       the reference pixel as defined in the input file),
      !       has been taken care of by writing the FULL image 
      !       with pixels outside the range specified by the 
      !       subimage remaining UNDEFINED (NaN).
      !       -- wasim, 09 Sep, 2009

      implicit none
      include '../INCLUDE/myfits_spec2rm.inc'

      integer*4 bitpixQ, naxisQ, naxesQ(max_axis)
      integer*4 bitpixU, naxisU, naxesU(max_axis)
      integer*4 bitpix, naxis, naxes(max_axis)
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
      real*4 specQ(maxchan),specU(maxchan)
      real*4 L_sq(maxchan),Q_now(maxchan),U_now(maxchan)
      character*8 junkchar
      integer*4 status, nchar
      logical anyflg
      logical cubeQ
      logical cubeU

      character*64 ctype 
      character*72 comment
      real*4 cval,dRM

      integer*4 rwmode
      character*72 infileQ, infileU, message
      character*72 outfile, outfileRM, outfilePA
      character*72 RMfile, QU_avefile
      character*72 subim_parfile
      character*72 path
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
      integer*4 i, kk, ix, iy
      integer*4 lnum
      integer*4 cnt1, tmp_cnt

      integer*4 cnt2, tmp_index, null_cnt(maxchan)
      real*4 tmp_avgQ, tmp_avgU

      ! Variables/Parameters for RM-extraction:
      real*4 fac
      integer*4 ofac
      real*4 RM(maxchan),p_ex(maxchan),phi_ex(maxchan)
      real*4 rp_ex(maxchan),rphi_ex(maxchan)
      real*4 ip_ex(maxchan),iphi_ex(maxchan)

      ! RFI related (list of bad-channels based on apriori info)
      real*4 bad_chan(maxchan)
      integer*4 nbad_chan
      logical remove_badchan
      character*72  badchan_file

      ! processing related:
      logical  line_cut
      character*72 add_req
      real*4  RM1 
      integer*4 nrm_out 
      logical fullrange , dummy_fullrange

      ! Some useless fitsio legacy stuff:
      integer*4 group, blocksize


!-------------------------------------------------------------------
      ! SANITY CHECKS:
      ! Compare the files containing the Q and U Cubes
      ! ans see if they are compatible with each other:

      if(iargc().lt.3)then
              write(*,*)'------------------------------------------'
              write(*,*)' You can have a blind run command as...'
              write(*,*)' myfits_spec2rm <infileQ> <infileU>  <outfile> 
     -<badchan-req> <badchan-file> <subim-req> <subim-parfile> '
              write(*,*)'  '
              write(*,*)'infileQ [NO DEFAULT]: '
              write(*,*)' FITS Q-CUBE File (without PATH) '
              write(*,*)' The path to this file is expected '
              write(*,*)' to be defined in the "datadir.par"'
              write(*,*)' file kept in the PAR/ directory '
              write(*,*)'  '
              write(*,*)'infileU [NO DEFAULT]: '
              write(*,*)' FITS U-CUBE File (without PATH) '
              write(*,*)' The path to this file is expected '
              write(*,*)' to be defined in the "datadir.par"'
              write(*,*)' file kept in the PAR/ directory '
              write(*,*)'  '
              write(*,*)'outfile [NO DEFAULT]: '
              write(*,*)' Generic name for outfiles (w/o PATH) '
              write(*,*)' Two cubes -- one for the RM and the '
              write(*,*)' other for the Polarisation Angle '
              write(*,*)' will be written out with appropriate '
              write(*,*)' extensions appended to the generic '
              write(*,*)' name. '
              write(*,*)' These files will be written to the '
              write(*,*)' PAR/ directory by default.'
              write(*,*)'  '
              write(*,*)'Optional arguments:'
              write(*,*)'  '
              write(*,*)'badchan-req [DEFAULT "N"]: '
              write(*,*)' [Yes or No] '
              write(*,*)' To get information on BAD-CHANNELS '
              write(*,*)' based on apriori information '
              write(*,*)'  '
              write(*,*)'badchan-file [NO DEFAULT]: '
              write(*,*)'[Required only when badchan-req = Yes]: '
              write(*,*)' File containing a list of BAD '
              write(*,*)' SPECTRAL CHANNELS.'
              write(*,*)' This file is expected by default '
              write(*,*)' to be in the DATA/ directory. '
              write(*,*)'  '
              write(*,*)'subim-req [DEFAULT "N"]: '
              write(*,*)' [Yes or No] '
              write(*,*)' To determine whether a subsection '
              write(*,*)' of the images is required or the '
              write(*,*)' full images. '
              write(*,*)'  '
              write(*,*)'subim-parfile [NO DEFAULT]: '
              write(*,*)'[Required only when subim-req = Yes]: '
              write(*,*)' File containing the parameters of '
              write(*,*)' the required subimage '
              write(*,*)' This file is expected by default '
              write(*,*)' to be in the PAR/ directory. '
              write(*,*)'  '
              write(*,*)'add-request : '
              write(*,*)' You can make some additional requests: '
              write(*,*)' using this string. Valid requests as '
              write(*,*)' of now are: '
              write(*,*)' 1) single_cut: to be used when you intend'
              write(*,*)'                to write out the Q,U and '
              write(*,*)'                RM-spectra for only a single'
              write(*,*)'                "cut" in the sky. By "cut"'
              write(*,*)'                I mean all pixels for eg.,'
              write(*,*)'                having constant Dec value.'
              write(*,*)'  '
              write(*,*)' You can have a blind run command as...'
              write(*,*)' myfits_spec2rm <infileQ> <infileU>  <outfile> 
     -<badchan-req> <badchan-file> <subim-req> <subim-parfile> '
              write(*,*)'------------------------------------------'
              write(*,*)'  '
              write(*,*)'Input Q-Cube Fits File:'
              write(*,*)' '
              read(*,'(a)')infileQ
              write(*,*)'Input U-Cube Fits File:'
              write(*,*)' '
              read(*,'(a)')infileU
              write(*,*)'Name of output Fits File:'
              write(*,*)' '
              read(*,'(a)')outfile
              write(*,*)' '
              write(*,*)' Do you have the list of BADCHANS ? (y/n)'
              read(*,'(a)')yorn
              write(*,*)' '
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      write(*,*)"Badchan file for RFI removal:"
                      read(*,'(a)')badchan_file
                      write(*,*)' '
                      remove_badchan = .true.
              endif
              write(*,*)' Do you want a subsection of the images? (y/n)'
              write(*,*)' '
              read(*,'(a)')yorn
              write(*,*)' '
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      write(*,*)"Par file for subimage region:"
                      read(*,'(a)')subim_parfile
                      write(*,*)' '
                      subim = .true.
              endif
      else if(iargc().eq.4)then
              subim = .false.
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              call getarg(4,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      write(*,*)"Badchan file for RFI removal: "
                      read(*,'(a)')badchan_file
                      write(*,*)' '
                      remove_badchan = .true.
              else
                      remove_badchan = .false.
              endif
      else if(iargc().eq.5)then
              subim = .false.
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              call getarg(4,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(5,badchan_file)
                      remove_badchan = .true.
              else
                      write(*,*)"Ignoring Arg(5) since Arg(4) indicates 
     -that no badchan removal is needed!"
                      write(*,*)" "
                      remove_badchan = .false.
              endif
      else if(iargc().eq.6)then
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              call getarg(4,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(5,badchan_file)
                      remove_badchan = .true.
              else
                      write(*,*)"Ignoring Arg(5) since Arg(4) indicates 
     -that no badchan removal is needed!"
                      write(*,*)" "
                      remove_badchan = .false.
              endif
              call getarg(6,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      write(*,*)"Par file for subimage region:"
                      read(*,'(a)')subim_parfile
                      write(*,*)' '
                      subim = .true.
              else
                      subim = .false.
              endif
      else if(iargc().eq.7)then
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              call getarg(4,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(5,badchan_file)
                      remove_badchan = .true.
              else
                      write(*,*)"Ignoring Arg(5) since Arg(4) indicates 
     -that no badchan removal is needed!"
                      write(*,*)" "
                      remove_badchan = .false.
              endif
              call getarg(6,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(7,subim_parfile)
                      subim = .true.
              else
                      write(*,*)"Ignoring Arg(7) since Arg(6) indicates 
     -that full image is to be used..."
                      subim = .false.
                      write(*,*)" "
              endif
      else if(iargc().eq.8)then
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              call getarg(4,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(5,badchan_file)
                      remove_badchan = .true.
              else
                      write(*,*)"Ignoring Arg(5) since Arg(4) indicates 
     -that no badchan removal is needed!"
                      write(*,*)" "
                      remove_badchan = .false.
              endif
              call getarg(6,yorn)
              if(yorn.eq.'y'.or.yorn.eq.'Y')then
                      call getarg(7,subim_parfile)
                      subim = .true.
              else
                      write(*,*)"Ignoring Arg(7) since Arg(6) indicates 
     -that full image is to be used..."
                      subim = .false.
                      write(*,*)" "
              endif
              call getarg(8,add_req)
              ! Now decide what all additional requests 
              ! are to be met:
              if(index(add_req,'single_cut').gt.0)then
                      line_cut = .true.
              else
                      line_cut = .false.
              endif
      else if(iargc().eq.3)then
              call getarg(1,infileQ)
              call getarg(2,infileU)
              call getarg(3,outfile)
              subim = .false.
              remove_badchan = .false.
      end if
      ! Do not write the additional files if the 
      ! entire cube is being processed:
      if(.not.subim)then
              line_cut = .false.
      endif

      ! Read the PATH-to-DATA from a file:
      open(11,file='../PAR/datadir.par',status='old',err=401)
      goto 402
401   continue
      write(*,*)" "
      write(*,*)"Error opening file: datadir.par"
      write(*,*)" Ensure that you have defined "
      write(*,*)" the path to your data in a file"
      write(*,*)" named 'datadir.par' kept in the"
      write(*,*)" PAR/ directory and re-run the "
      write(*,*)" code."
      write(*,*)" "
      write(*,*)"Quitting now... "
      stop

402   continue
      read(11,'(a)')junkchar  ! Reading the comment line
      read(11,*)lnum
      if (lnum.lt.1)then
              write(*,*)" "
              write(*,*)"Invalid path number specified in "
              write(*,*)"row #2 of the 'datadir.par' file "
              write(*,*)" "
              close(11)
              write(*,*)"Closing the 'datadir.par'... "
              write(*,*)"Quitting now... "
              write(*,*)" "
              stop
      endif
      if (lnum.gt.1)then
              do i = 1,lnum - 1
                 read(11,'(a)')junkchar
              enddo
      endif
      read(11,'(a)')path
      close(11)


      infileQ(1:) = path(1:nchar(path))//infileQ(1:nchar(infileQ))
      infileU(1:) = path(1:nchar(path))//infileU(1:nchar(infileU))

      outfileRM(1:) = '../DATA/'//outfile(1:nchar(outfile))//'.RMCUBE'
      outfilePA(1:) = '../DATA/'//outfile(1:nchar(outfile))//'.PACUBE'
      QU_avefile(1:) = '../DATA/'//outfile(1:nchar(outfile))//'.AVGQU'
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
                      goto 9999
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
              goto 9999
      endif
      z1 = czval_im - nz_1st*zinc_im
      zn = z1 + (nz_totpix - 1)*zinc_im
      call linspace(z1,zn,nz_totpix,zval)

      do i = 1,nz_out
         L_sq(i) = (conv_fac/zval(nz_out-i+1))**2  
         write(13,*)L_sq(i)
         Q_now(i) = 0.0
         U_now(i) = 0.0
      enddo

      RM1 = 50.0
      nrm_out = 10
      fullrange = .false.
      dummy_fullrange = .true.
      if(fullrange)then
              nrm_out = nz_out
      endif
      ! Make a dummy call to extract general to get the 
      ! RM-values. We need to write the reference pixel 
      ! and its corresponding value in the output FITS 
      ! file for the RM-cube...
!      call extract_general(L_sq,Q_now,U_now,nz_out,fac,ofac,
!     -                           RM1,nrm_out,RM,
!     -                           p_ex,phi_ex,
!     -                           rp_ex,rphi_ex,
!     -                           ip_ex,iphi_ex, 
!     -                           dummy_fullrange)
!
!      dRM = (RM(nz_out) - RM(1))/real(nz_out-1)
      dRM = fac/(L_sq(nz_out) - L_sq(1))
      RM(1) = -0.5*real(nz_out - 1)*dRM
      RM(nz_out) = 0.5*real(nz_out - 1)*dRM
      RM(nrm_out) = RM1 + real(nrm_out - 1)*dRM
      write(*,*)"dRM: ",dRM
!      stop

      write(*,*)" "
      if(fullrange)then
              write(17,*)"  n_RM: ",nz_out
              write(17,*)"max RM: ",RM(nz_out)
              write(17,*)"min RM: ",RM(1)
      else
              write(17,*)"  n_RM: ",nrm_out
              write(17,*)"max RM: ",RM(nrm_out) !RM1 + real(nout - 1)*dRM
              write(17,*)"min RM: ",RM1
      endif
      write(17,*)"  d_RM: ",dRM

      ! Read data from the files...
      ! FTGSV[BIJKED](unit,group,naxis,naxes,fpixels,lpixels,incs,nullval, 
      !               > array,anyflg,status)
      ! Get an arbitrary data subsection from the data array. 
      ! Undefined pixels in the array will be set equal to the 
      ! value of 'nullval', unless nullval=0 in which case no 
      ! testing for undefined pixels will be performed. 

      ! Internal definition for subroutine FTGSVE:
      !
      ! subroutine ftgsvd(iunit,colnum,naxis,naxes,blc,trc,inc,nulval,
      ! >array,anynul,status)
      !
      ! iunit   i  fortran unit number
      ! colnum  i  number of the column to read from
      ! naxis   i  number of dimensions in the FITS array
      ! naxes   i  size of each dimension.
      ! blc     i  'bottom left corner' of the subsection to be read
      ! trc     i  'top right corner' of the subsection to be read
      ! inc     i  increment to be applied in each dimension
      ! nulval  i  value that undefined pixels will be set to
      ! array   i  array of data values that are read from the FITS file
      ! anynul  l  set to .true. if any of the returned values are undefined
      ! status  i  output error status

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
      call ftphpr(41,simple,bitpix,naxis,naxes,0,1,extend,status)
      call ftphpr(42,simple,bitpix,naxis,naxes,0,1,extend,status)

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
      call ftpkys(41,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      call ftpkys(42,"ctype1",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)

      call ftgkye(21,"crval1",cval,comment,status)
      call ftpkye(41,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crval1",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftgkye(21,"crpix1",cval,comment,status)
      call ftpkye(41,"crpix1",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crpix1",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftgkye(21,"cdelt1",cval,comment,status)
      call ftpkye(41,"cdelt1",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"cdelt1",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftgkys(21,"ctype2",ctype,comment,status)
      call ftpkys(41,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)
      call ftpkys(42,"ctype2",ctype(1:nchar(ctype)),comment(1:nchar(comm
     -ent)),status)

      call ftgkye(21,"crval2",cval,comment,status)
      call ftpkye(41,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crval2",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftgkye(21,"crpix2",cval,comment,status)
      call ftpkye(41,"crpix2",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"crpix2",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftgkye(21,"cdelt2",cval,comment,status)
      call ftpkye(41,"cdelt2",cval,decimals,comment(1:nchar(comment)),
     -status)
      call ftpkye(42,"cdelt2",cval,decimals,comment(1:nchar(comment)),
     -status)

      call ftpkys(41,"ctype3","RM--rad/m2","3rd axis type",status)
      call ftpkys(42,"ctype3","RM--rad/m2","3rd axis type",status)

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


      RMfile = outfile(1:nchar(outfile))//'.RMSPEC'

      open(16,file=RMfile,status='unknown',form='unformatted',
     -              access='direct',recl=4*nrm_out)
      open(51,file=QU_avefile,status='unknown')
      write(*,*)"The spectral average of Q and U for each image"
      write(*,*)"pixel will be noted in: ",
     -                QU_avefile(1:nchar(QU_avefile))
      write(*,*)" "
      write(51,*)"### This file logs the spectral average of Q and U "
      write(51,*)"### for each pixel in the image. "
      write(51,*)"### format used: RApix, DECpix, meanQ, meanU "
      write(51,*)"####" 

      tmp_cnt = 0
      cnt1 = 0
      do ix = xpix_beg,xpix_end,incs(1)
         fpixels(1) = ix
         lpixels(1) = ix
         do iy = ypix_beg,ypix_end,incs(2)
            cnt1 = cnt1 + 1
            !
            ! initialise the spec*() arrays
            !do i = 1,nz_out
            !   specQ(i) = 0.0
            !   specU(i) = 0.0
            !enddo
            fpixels(2) = iy
            lpixels(2) = iy

            if(bitpix.eq.8)then
              call FTGSVB(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVB(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            else if (bitpix.eq.16)then
              call FTGSVI(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVI(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            else if (bitpix.eq.32)then
              call FTGSVJ(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVJ(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            else if (bitpix.eq.64)then
              call FTGSVK(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVK(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            else if (bitpix.eq.-32)then
              call FTGSVE(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVE(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            else if (bitpix.eq.-64)then
              call FTGSVD(21,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specQ,anyflg,status)
              call FTGSVD(22,group,naxis,naxes,fpixels,lpixels,incs,
     -                    nullval,specU,anyflg,status)
            endif

            do i = 1,nz_out
               Q_now(i) = specQ(nz_out-i+1) 
               U_now(i) = specU(nz_out-i+1)
            enddo

            !----------------------------------------------------
            ! Replace the bad-channel data with null-values
            if(remove_badchan)then
                    do i = 1,nbad_chan
                       !write(*,*)"bad_chan: ",bad_chan(i)
                       tmp_index = bad_chan(i)
                       Q_now(tmp_index) = nullval
                       U_now(tmp_index) = nullval
                    enddo
            endif
            !----------------------------------------------------
            ! Replace the pixels containing nullval with the mean
            ! of the "good" data pixels
            ! Dubious mistake (not serious for GOOD data) caught!
            ! Proper replacement of BAD data was not being done.
            ! Hopefully corrected!
            ! WR, 22 Dec, 2009
            cnt2 = 0

            do i = 1,nz_out
               if(Q_now(i).eq.nullval.or.U_now(i).eq.nullval)then
                       cnt2 = cnt2 + 1
                       null_cnt(cnt2) = i
                       Q_now(i) = 0.0
                       U_now(i) = 0.0
               else if(Q_now(i).eq.0.0.and.U_now(i).eq.0.0)then
                       cnt2 = cnt2 + 1
                       null_cnt(cnt2) = i
               else
                       Q_now(i) = Q_now(i) 
                       U_now(i) = U_now(i)
               endif
            enddo
            !write(*,*)cnt2,' out of',nz_out,' zpix have nullval'
            !write(*,*)' '
            call mean(Q_now,nz_out,tmp_avgQ)
            call mean(U_now,nz_out,tmp_avgU)
            tmp_avgQ = tmp_avgQ*real(nz_out)/real(nz_out-cnt2)
            tmp_avgU = tmp_avgU*real(nz_out)/real(nz_out-cnt2)
            ! Book-keeping for the mean-Q and mean-U
            ! values:
            write(51,*)ix,iy,tmp_avgQ, tmp_avgU

            do i = 1,cnt2
               tmp_index = null_cnt(i)
               Q_now(tmp_index) = tmp_avgQ
               U_now(tmp_index) = tmp_avgU
            enddo
            !do i = 1,nz_out
            !  write(*,*)Q_now(i),U_now(i)
            !enddo
            !----------------------------------------------------

            ! Perform the tomography now:
            ! subroutine extract_general(t,ryt,iyt,npts,fac,ofac, 
!     -               omega1,nout,
!     -               > omega,
!     -               p_ex,phi_ex,
!     -               rp_ex, rphi_ex, 
!     -               ip_ex, iphi_ex)
            call extract_general(L_sq,Q_now,U_now,nz_out,fac,ofac,
     -                           RM1,nrm_out,RM,
     -                           p_ex,phi_ex,
     -                           rp_ex,rphi_ex,
     -                           ip_ex,iphi_ex, 
     -                           fullrange)
            !write(*,*)"   nrm_out : ",nrm_out
            !write(*,*)"  RM(1),RM1: ",RM(1),RM1
            !write(*,*)"RM(nrm_out): ",RM(nrm_out)
            !write(*,*)"        dRM: ",dRM

            if(line_cut)then
                    write(14,*)"### ix,iy: ",ix,iy
                    write(15,*)"### ix,iy: ",ix,iy
                    do i = 1,nrm_out
                       write(14,*)RM(i),p_ex(i),phi_ex(i)
                       write(15,*)L_sq(i),Q_now(i), U_now(i)
                    enddo
                    tmp_cnt = tmp_cnt + 1
                    write(16,rec=tmp_cnt)(Q_now(i),i=1,nrm_out)
                    tmp_cnt = tmp_cnt + 1
                    write(16,rec=tmp_cnt)(U_now(i),i=1,nrm_out)
                    tmp_cnt = tmp_cnt + 1
                    write(16,rec=tmp_cnt)(p_ex(i),i=1,nrm_out)
                    tmp_cnt = tmp_cnt + 1
                    write(16,rec=tmp_cnt)(phi_ex(i),i=1,nrm_out)
            endif
            !! Write the FITS RM-CUBE now:
            call ftpsse(41,group,naxis,naxes,fpixels,lpixels,p_ex,
     -status)
            call ftpsse(42,group,naxis,naxes,fpixels,lpixels,phi_ex,
     -status)

            if(mod(cnt1-1,1000).eq.0)then
                 write(*,*)"doing ",cnt1," out of",nx_out*ny_out 
                 !write(*,*)"------------------------------ "
                 !write(*,*)"Warning: "
                 !write(*,*)"Test mode on..."
                 !write(*,*)"Comment the lines to suppress "
                 !write(*,*)"unnecessary file-writing! "
                 ! REMOVE the lines: Test-1 start to Test-1 stop
                 !write(*,*)"------------------------------ "
            endif
         enddo     ! end of iy loop
      enddo        ! end of ix loop
      !=======================================================


9999  continue
      close(16)
      close(51)

      ! CLOSE THE FITS FILES:
      call FTCLOS(21,status)
      !write(*,*)"---------------------------"
      !write(*,*)"Current subroutine: FTCLOS "
      !write(*,*)"STATUS = ",status
      !write(*,*)" "
      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      endif

      call FTCLOS(22,status)
      !write(*,*)"---------------------------"
      !write(*,*)"Current subroutine: FTCLOS "
      !write(*,*)"STATUS = ",status
      !write(*,*)" "
      if (status .gt. 0)then
              write(*,*)"Problem closing Q-file"
              call printerror(status)
      endif
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

      ! -----------------------------------------------------------------


      end

      include '/usr/lib/subroutine_lib/nchar.f'
      include 'myfits_info.f'
      !include 'extract_general.f'
      include 'extract_general_v2.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
