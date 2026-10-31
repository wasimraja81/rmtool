chelp+
      !------------------------------------------------------
      ! This code may be used to select a specified column 
      ! from an ASCII file containing (nrow x ncol) data pts
      ! and detect outliers in the specified column. 
      ! It outputs the index(or, row numbers) of the bad-data 
      ! detected into an outfile specified. 
      !          -- wr, 05 Sep, 2010
      !------------------------------------------------------
chelp-      


      integer*4      maxchan, maxpfit, maxcol
      parameter      (maxchan=8192, maxpfit=128, maxcol=64)
      integer*4      nright, nleft, nchan, order_h, order_p 
      integer*4      i, j, icnt
      integer*4      n_badcan, bad_chan(maxchan), ncol, colnum
      logical        lr_exclude, silent
      character*128  infile, outfile, templine
      character*128  data_tag
      real*4         ydata(maxchan), fit_array(maxchan),
     -               best_fit_param(maxpfit), thresh_fit
      real*4         xdata(maxchan), tmpdata(maxcol), tmp_rms
      character*1    junkchar


      character*120  xlabel, ylabel, title
      real*4         xmin, xmax, ymin, ymax
      real*4         tmpnum





       !--------------------------------------------------------

       if(iargc().ne.4)then
               write(*,*)"Usage: detect_rfi <infile> <tot col> <use col>
     - <outfile>"
               write(*,*)" "
               write(*,*)"infile [No Default] "
               write(*,*)"  This file contains the data containing "
               write(*,*)"  possible outliers. " 
               write(*,*)" "
               write(*,*)"Total Columns in infile [No Default] "
               write(*,*)"  The total number of columns in the "
               write(*,*)"  infile chosen."
               write(*,*)" "
               write(*,*)"Use Column [No Default] "
               write(*,*)"  The column number containing the data to "
               write(*,*)"  be used for outlier detection. This is "
               write(*,*)"  helpful in case you have multiple columns"
               write(*,*)"  in your infile."
               write(*,*)"  It is however expected that the number of "
               write(*,*)"  rows is equal for each column. "
               write(*,*)" "
               write(*,*)"outfile [No Default] "
               write(*,*)"  This file contains the list of row-numbers "
               write(*,*)"  containing the outliers. If equi-spaced "
               write(*,*)"  data was input through infile, then this"
               write(*,*)"  file contains the indices corresponding "
               write(*,*)"  to the outliers. For eg., if a spectra, "
               write(*,*)"  is input, then the output would contain "
               write(*,*)"  the BAD-channel-numbers. "
               write(*,*)" "
               write(*,*)"Quitting now..."
               stop
       else
               call getarg(1,infile)
               infile = infile(1:nchar(infile))
               call getarg(2,templine) 
               templine = templine(1:nchar(templine))
               read(templine,*)ncol
               call getarg(3,templine) 
               templine = templine(1:nchar(templine))
               read(templine,*)colnum
               write(*,*)"ncol, usecol: ", ncol, colnum
               call getarg(4,outfile)
               outfile = outfile(1:nchar(outfile))
       endif
       open(11,file=infile,status='old',err=101)
       goto 102
101    continue
       write(*,*)" "
       write(*,*)"Error opening file: ",
     -           infile(1:nchar(infile))
       write(*,*)" "
       write(*,*)"Quitting now..."
       stop

102    continue
       icnt = 0
       read(11,*)junkchar  ! read the comment line
       do while(.true.)
          icnt = icnt + 1
          read(11,*,err=201)(tmpdata(i),i=1,ncol)
          ydata(icnt) = tmpdata(colnum)
          xdata(icnt) = real(icnt)
          write(*,*)"icnt, ydata: ",icnt,ydata(icnt)
       enddo

201    continue
       nchan = icnt - 1
       write(*,*)"Length of spectra: ",nchan

       !----------------------------------------------
       ! We wish to fit all points, so fill fit_array
       ! appropriately: 
       do i = 1,nchan
         fit_array(i) = 1.0
       enddo
       !----------------------------------------------
      ! Locate the bad channels now:
      nleft = 1
      nright = nchan
      lr_exclude = .true.
      silent = .true.
      data_tag = ' '

      thresh_fit = 2.0
      order_h = 3
      order_p = 3

      call poly_harm_fit(nchan,xdata, ydata, nleft, nright,
     -                   lr_exclude, order_h, order_p,
     -                   thresh_fit, silent, data_tag,fit_array,
     -                   best_fit_param)

      n_fitparams = 2*order_h + order_p + 1
      tmp_rms = best_fit_param(n_fitparams + 1)
      n_badchan = 0
       open(31,file=outfile,status='unknown',err=301)
       goto 302
301    continue
       write(*,*)" "
       write(*,*)"Error opening file: ",
     -           outfile(1:nchar(outfile))
       write(*,*)" "
       write(*,*)"Quitting now..."
       stop

302    continue
      do j = 1,nchan
         if(abs(ydata(j)-fit_array(j)).gt.thresh_fit*tmp_rms)then
                 n_badchan = n_badchan + 1
                 bad_chan(n_badchan) = j
                 if(j.eq.1)then
                         write(*,*)j,"-st channel is bad " 
                 else
                         write(*,*)j,"-th channel is bad " 
                 endif
                 write(31,*)j 
         endif
      enddo
      write(*,*)"Robust RMS: ",tmp_rms
      write(*,*)"Number of bad channels: ",n_badchan
      ! End of RFI detection
      !-------------------------------------------------------------
      ! Some plots: 
      call minima(xdata,nchan,xmin)
      call maxima(xdata,nchan,xmax)


      call minima(ydata,nchan,ymin)
      call minima(fit_array,nchan,tmpnum)

      if(tmpnum.lt.ymin)then
              ymin = tmpnum
      endif

      call maxima(ydata,nchan,ymax)
      call maxima(fit_array,nchan,tmpnum)

      if(tmpnum.gt.ymax)then
              ymax = tmpnum
      endif

      call pgbeg(0,'/xs',1,1)
      call pgenv(xmin,xmax,ymin,ymax,0,1)

      xlabel = 'Channel Number'
      ylabel = 'Linear Polarised Intensity'
      title = 'Bad Channel Detection'
      call pglab(xlabel,ylabel,title)
      call pgsci(3)
      call pgline(nchan, xdata, ydata)
      call pgsci(2)
      call pgline(nchan, xdata, fit_array)

      call pgend

      close(11)
      end
      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/poly_harm_fit.f'
