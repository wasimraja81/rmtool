chelp+
      !---------------------------------------------------------
      ! This code is intended to generate blank FITS images
      ! using an existing image. The blanked images correspond 
      ! to frequency-planes that were affected by RFI and 
      ! hence could not be imaged. 
      ! Their need arises when I try to combine these image
      ! planes at different frequencies into a spectral cube 
      ! using AIPS/MIRIAD tasks. 
      ! If these planes are left missing, AIPS/MIRIAD tasks 
      ! are unable to define the frequency-axis appropriately. 
      ! They require contiguous frequencies to produce the 
      ! cubes with well-defined axes values. 
      ! Of course the Frequency-pixel values are appropriately 
      ! modified in the header. Only the data are blanked. 
      ! 
      ! inputs: * reference fits file-name from where most header 
      !           parameters are to be copied. 
      !
      !         * output file-name to be written -- this file is 
      !           by default assumed to be the file corresponding
      !           to the 1st frequency channel to be used in the
      !           cubes.
      !
      !         * Frequency Channel number for which the output, 
      !           blank image is to be written.
      !
      ! output: * Blank image corresponding to the frequency 
      !           channel specified.
      !
      !     NB: * I have not bothered to change the min-max values 
      !           (KEYWORDS: DATAMIN and DATAMAX) in the header of
      !           the output blanked image.
      !
      !         * The keyvalues are written in double precision.
      ! 
      !                              -- wasim, 26 Aug, 2010
      !---------------------------------------------------------
chelp-      

      implicit none


      integer*4    status, nchar
      integer*4    rwmode
      character*172 infile, outfile
      character*172 tempstring
      real*8       f1, del_f
      integer*4    this_chan

    
      ! Some useless fitsio legacy stuff:
      integer*4    blocksize, decimals
      character*1  comment
      real*8       keyval_new


!-------------------------------------------------------------------


      if(iargc().ne.3)then
              write(*,*)"Usage: "
              write(*,*)"    blank_image <infile> <outfile> <this_chan>"
              write(*,*)"Quitting now..."
              stop
      else
              call getarg(1,infile)
              infile = infile(1:nchar(infile))
              call getarg(2,outfile)
              outfile = outfile(1:nchar(outfile))
              call getarg(3,tempstring)
              read(tempstring(1:),*)this_chan
      endif

      ! Initialise STATUS to zero:
      status = 0

      call FTOPEN(21,infile,rwmode,blocksize,status)

      if(status.ne.0)then
              write(*,*)" "
              write(*,*)"Error opening input FITS file: ",
     -                  infile(1:nchar(infile))
              write(*,*)"status = ", status
              stop
      else
              write(*,*)" "
              write(*,*)"inp-FITS file chosen: ",infile(1:nchar(infile))
      endif
      ! Read the frequency information from the reference FITS file
      ! (Remember we want the reference image to be the 1st channel)
      call FTGKYD(21,"crval3",f1,tempstring,status)
      call FTGKYD(21,"cdelt3",del_f,tempstring,status)
      keyval_new = f1 + dble(this_chan - 1)*del_f
      write(*,"(a,f16.5)")"   f1: ",f1
      write(*,"(a,f16.5)")"del_f: ",del_f
      write(*,"(a,f16.5)")"f_new: ",keyval_new

      call FTINIT(31,outfile,blocksize,status)
      if(status.ne.0)then
              write(*,*)" "
              write(*,*)"Error opening output FITS file: ",
     -                  outfile(1:nchar(outfile))
              write(*,*)"status = ", status
              stop
      else
              write(*,*)" "
              write(*,*)"out-FITS file chosen: ",
     -                   outfile(1:nchar(outfile))
      endif

      ! BEGIN BUSINESS:
      ! Copy the header from infile to outfile:
      call FTCPHD(21, 31, status)
      if (status .gt. 0)then
              write(*,*)"Problem using FTCPHD..."
              call printerror(status)
              write(*,*)"Quitting now..."
              stop
      else
              write(*,*)"output data written to: ",
     -                        outfile(1:nchar(outfile))
      endif

      ! CLOSE THE FITS FILES:
      call FTCLOS(21,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing FITS-file:",
     -                  infile(1:nchar(infile))
              call printerror(status)
              write(*,*)"Quitting now..."
              stop
      endif

      call FTCLOS(31,status)
      if (status .gt. 0)then
              write(*,*)"Problem closing output-FITS file", 
     -                  outfile(1:nchar(outfile))
              call printerror(status)
              write(*,*)"Quitting now w/o modifying headers!"
              stop
      endif

      ! -----------------------------------------------------------------
      ! If all went well, modify the frequency in the header now:
      rwmode = 1
      blocksize = 0
      decimals = 11
      call FTOPEN(41,outfile,rwmode,blocksize,status)
      if (status .gt. 0)then
              write(*,*)"Problem opening file: ",
     -                outfile(1:nchar(outfile))
              call printerror(status)
              write(*,*)"Quitting now..."
              stop
      else
              comment = '&' ! This is to keep the comment same as the 
                            ! input FITS file
              !call FTMKLS(41,"crval3",keyval_new,comment,status)
             call FTMKYD(41,"crval3",keyval_new,decimals,comment,status)
      endif
      call FTCLOS(41,status)
      write(*,*)"All went well..."

      end

      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
