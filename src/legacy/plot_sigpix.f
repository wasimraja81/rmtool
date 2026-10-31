chelp+
      !--------------------------------------------------------------
      ! Code to plot data from ascii files arrranged in rows and col. 
      ! 
      !                                         --wr, 12 Apr, 2012 
      !--------------------------------------------------------------
chelp-

      implicit none 

      integer*4      maxread_per_pass 
      parameter      (maxread_per_pass = 2048) 
      character      infile*128, junkchar*1, templine*128, psfile*220, 
     -               outfile*220 
      real*4         tmp_arr(maxread_per_pass), atmp 
      integer*4      nchar 
      logical        psreq 

      integer*4      nr, nr1, nr2, nc, nc1, nc2, nread, 
     -               ntot_read, nc_max   
      integer*4      i, ipass, j, npass, nrem  

      ! PLOT related: 
      real*4         xmin, xmax, xmean, ymin, ymax, ymean, xtmp, ytmp 
      real*4         xarr(maxread_per_pass), yarr(maxread_per_pass) 
      character*62  xlabel, ylabel, title 
      
      ! BINNING RELATED: 
      integer*4     maxbin 
      parameter     (maxbin=1000)
      integer*4     ibinx, ibiny, nbinx, nbiny 
      real*4        dnx, dny, dist(maxbin,maxbin), xdist(maxbin)
      integer*4     reclen, data_precision 
      ! 
      ! COFIG FILE generation 
      real*4      begx, endx, begy, endy 
      integer*4   nset, ncomments_passed  
      character*62   xlabel_left, xlabel_bottom, 
     -               ylabel_left, ylabel_bottom, 
     -               title_pltpar 
      character*32   comments(32) 
      character*220  path_for_cfg 



      if(iargc().lt.6)then
              write(*,*)"Usage: "
              write(*,*)"plot_sigpix < infile > < nrow1 > < nrow2 > < nc
     -ol1 > < ncol2 > < psreq(y/n) >"
              write(*,*)" "
              write(*,*)" Quitting now..."
              stop
      else
              call getarg(1,infile)
              infile = infile(1:nchar(infile))
              call getarg(2,templine)
              read(templine,*)nr1
              call getarg(3,templine)
              read(templine,*)nr2
              call getarg(4,templine)
              read(templine,*)nc1
              call getarg(5,templine)
              read(templine,*)nc2
              call getarg(6,junkchar)

              if(junkchar.eq.'y'.or.junkchar.eq.'Y')then
                      psreq = .true.
              else
                      psreq = .false. 
              endif 
      endif



      nbinx = 50  
      nbiny = 50  

      open (21,file=infile,status='old')

      read(21,*)nr,nc    ! total number of rows and col in file 
      read(21,*)junkchar ! comment line 
      !write(*,*)nr,nc    ! total number of rows and col in file 
      close(21) 
      write(*,*)"Total number of rows in data file: ",nr 
      write(*,*)"Total number of cols in data file: ",nc 
      write(*,*)" "
      if (nr2.lt.nr1)then
              write(*,*)"last row < 1st row found!!"
              write(*,*)"Change the order in the command line."
              write(*,*)"Quitting now..."
              stop
      endif
      if(nr1.le.0)then
              nr1 = 1
      endif
      if(nr2.le.0)then
              nr2 = nr
      endif
      ntot_read = nr2 - nr1 + 1 
      nc_max = max(nc1,nc2)
      if (nr1.gt.nr2)then
              write(*,*)"nr2 must be greater than nr1..."
              write(*,*)"Quitting now..."
              stop
      endif
      if (nr1.gt.nr)then
              write(*,*)"The nr1-th to read exceeds maxrows in data ..."
              write(*,*)"Quitting now..."
              stop
      endif
      if (nr2.gt.nr)then
              write(*,*)"The nr2-th to read exceeds maxrows in data ..."
              write(*,*)"Quitting now..."
              stop
      endif
      if (ntot_read.gt.nr)then
              write(*,*)"nrows to read exceeds maxrows in data ..."
              write(*,*)"Quitting now..."
              stop
      endif
      if (nc_max.gt.nc)then
              write(*,*)"ncol to read exceeds maxcols in data ..."
              write(*,*)"Quitting now..."
              stop
      endif

      ! Read the unwanted rows... 
      open (21,file=infile,status='old') 
      read(21,*)junkchar ! comment line 
      read(21,*)junkchar ! info line  
      do i = 1,nr1-1 
         read(21,*)junkchar 
      enddo 

      ! Scan the file for maxima and minima global to the 
      ! section intended for plotting 
      !
      ! 1st row for initialization 
      read(21,*)(tmp_arr(j),j=1,nc_max) 
      xmax = tmp_arr(nc1) 
      xmin = xmax 
      xmean = xmax

      ymax = tmp_arr(nc2) 
      ymin = ymax 
      ymean = ymax


      do i=2,ntot_read 
         read(21,*,err=101,end=101)(tmp_arr(j),j=1,nc_max) 
         xtmp = tmp_arr(nc1) 
         ytmp = tmp_arr(nc2) 

         xmean = xmean + xtmp 
         ymean = ymean + ytmp 

         if(xtmp.gt.xmax)then 
                 xmax = xtmp 
         endif 
         if(xtmp.lt.xmin)then 
                 xmin = xtmp 
         endif 

         if(ytmp.gt.ymax)then
                 ymax = ytmp 
         endif
         if(ytmp.lt.ymin)then
                 ymin = ytmp
         endif
      enddo
101   continue
      close(21) 

      if (ntot_read .gt. maxread_per_pass)then
              npass = int(ntot_read/maxread_per_pass)
              nrem = ntot_read - npass*maxread_per_pass 
              nread = maxread_per_pass
      else
              npass = 1 
              nrem = 0 
              nread = ntot_read 
      endif
      xmean = xmean/ntot_read
      ymean = ymean/ntot_read

      write(*,*)" " 
      write(*,*)"Number of passes in which data wll be read : ",npass
      write(*,*)"     Number of points to read/plot per pass: ",nread 
      write(*,*)"                               ntot to plot:",ntot_read
      write(*,*)"                         nrem for last plot:",nrem
      write(*,*)" "
      write(*,*)" mean of x variate: ",xmean
      write(*,*)" mean of y variate: ",ymean
      write(*,*)"=================================================="
      write(*,*)" " 
      !stop 

      !write(*,*)"    xrange: ",xmin, xmax
      !write(*,*)"    yrange: ",ymin, ymax
      open (21,file=infile,status='old')
      ! Read the unwanted rows...
      read(21,*)junkchar ! info line 
      read(21,*)junkchar ! comment line 
      do i = 1,nr1-1
         read(21,*)junkchar
      enddo
      write(xlabel,*)nc1
!      xlabel = 'Column Number: '//xlabel(1:nchar(xlabel))
      write(ylabel,*)nc2
!      ylabel = 'Column Number: '//ylabel(1:nchar(ylabel))

      if(nc1.eq.1)xlabel='RA'
      if(nc1.eq.2)xlabel='Dec'
      if(nc1.eq.3)xlabel='Freq'
      if(nc1.eq.4)xlabel='Stokes-I'
      if(nc1.eq.5)xlabel='Stokes-Q'
      if(nc1.eq.6)xlabel='Stokes-U'
      if(nc1.eq.7)xlabel='Stokes-V'
      if(nc1.eq.8)xlabel='LPol'


      if(nc2.eq.1)ylabel='RA'
      if(nc2.eq.2)ylabel='Dec'
      if(nc2.eq.3)ylabel='Freq'
      if(nc2.eq.4)ylabel='Stokes-I'
      if(nc2.eq.5)ylabel='Stokes-Q'
      if(nc2.eq.6)ylabel='Stokes-U'
      if(nc2.eq.7)ylabel='Stokes-V'
      if(nc2.eq.8)ylabel='LPol'


      write(*,*)"I am here..."
      outfile = infile(1:nchar(infile))//'_'//
     -                 xlabel(1:nchar(xlabel))//'_vs_'//
     -                 ylabel(1:nchar(ylabel))//'.bin'
      if(psreq)then
              psfile = infile(1:nchar(infile))//'_'//
     -                 xlabel(1:nchar(xlabel))//'_vs_'//
     -                 ylabel(1:nchar(ylabel))//'.ps'//'/cps'
      else
              psfile = '1/xs'
      endif

      title = 'Scatter Plot: '//ylabel(1:nchar(ylabel))//' vs. '//
     -                           xlabel(1:nchar(xlabel))

      write(templine,*)nc1
      xlabel = xlabel(1:nchar(xlabel))//
     -        ' [Col Num: '//templine(1:nchar(templine))//']'
      write(templine,*)nc2
      ylabel = ylabel(1:nchar(ylabel))//
     -        ' [Col Num: '//templine(1:nchar(templine))//']'

      call pgbeg(0,psfile,1,1)
      call pgask(.false.)
      call pgsci(8) 
      call pgenv(xmin,xmax,ymin,ymax,0,1)
      call pglab(xlabel,ylabel,title)

      do i = 1,nbinx
         do j = 1,nbiny
            dist(i,j) = 0.0 
         enddo
         xdist(i) = 0.0 
      enddo
      if(nbinx.gt.maxbin)nbinx=maxbin
      if(nbiny.gt.maxbin)nbiny=maxbin 

      dnx = (xmax-xmin)/real(nbinx)
      dny = (ymax-ymin)/real(nbiny)

      do ipass = 1,npass
         do i = 1,nread
            read(21,*)(tmp_arr(j),j=1,nc_max)
            xarr(i) = tmp_arr(nc1) 
            yarr(i) = tmp_arr(nc2) 
            !write(*,*)"I am here...",ipass,i 
            !write(*,*)(tmp_arr(j),j=1,nc_max)
            !------------------------------------------
            ! DATA BINNING SECTION 
            ! Bin the data for contour plots instead of 
            ! scatter plots: 
            ! an = a0 + (n - 1)d 
            ibinx = 1 + int((xarr(i) - xmin)/dnx) 
            ibiny = 1 + int((yarr(i) - ymin)/dny) 
            dist(ibinx,ibiny) = dist(ibinx,ibiny) + 1.0 
            xdist(ibinx) = xdist(ibinx) + 1.0 ! For normalizing the
                                              ! y-values in dist(ibinx,ibiny)
            !write(*,*)"ibinx,ibiny",ibinx,ibiny
         enddo
         call pgpt(nread,xarr,yarr,1)
      enddo

      do i = 1,nrem
         read(21,*)(tmp_arr(j),j=1,nc_max)
         !write(*,*)(tmp_arr(j),j=1,nc_max)
         xarr(i) = tmp_arr(nc1) 
         yarr(i) = tmp_arr(nc2) 
         !------------------------------------------
         ! DATA BINNING SECTION 
         ! Bin the data for contour plots instead of 
         ! scatter plots: 
         ! an = a0 + (n - 1)d 
         ibinx = 1 + int((xarr(i) - xmin)/dnx) 
         ibiny = 1 + int((yarr(i) - ymin)/dny) 
         dist(ibinx,ibiny) = dist(ibinx,ibiny) + 1.0 
         xdist(ibinx) = xdist(ibinx) + 1.0 ! For normalizing the
                                           ! y-values in dist(ibinx,ibiny)
      enddo
      do ibinx = 1,nbinx
         if (xdist(ibinx).eq.0)then
                 atmp = 1.0 
         else
                 atmp = xdist(ibinx)
         endif
         do ibiny = 1,nbiny
            
            dist(ibinx,ibiny) = dist(ibinx,ibiny)/atmp ! Weighting  the y-distribution 
                                                       ! by the number of x's encountered 
         enddo
      enddo
      call pgpt(nrem,xarr,yarr,1)

      call pgend 
      close(21) 

      data_precision = 4 
      reclen = nbiny*data_precision  
      open(21,file=outfile,form='unformatted',access='direct',
     -        recl=reclen,status='unknown')
      do i = 1,nbinx
         write(21,rec=i)(dist(i,j),j=1,nbiny)
      enddo

      close(21) 

      ! Now generate the cfg files to enable quick plotting: 
            ! Write a handy cfg file for display routines to use: 
      ! Determine the path where the data is to be sent: 
      call system("pwd >path.tmp")
      open(21,file='path.tmp',status='old',err=301)
      goto 302
301   write(*,*)"Error opening file path.tmp"
      write(*,*)"Quitting now..."
      stop
302   continue

      read(21,'(a)')path_for_cfg
      path_for_cfg = path_for_cfg(1:nchar(path_for_cfg))//'/'
      close(21) 


      outfile = outfile(1:nchar(outfile))
      begy = ymin
      endy = ymax 
      
      begx = xmin
      endx = xmax
      nset = 1 
      call generate_cfg_for_display(path_for_cfg,outfile,nset,nbiny,
     -               data_precision,begx,endx,nbinx,begy,endy)
      ! Now the PLTPAR file generation: 

      xlabel_left = ' '
      ylabel_left = ylabel(1:nchar(ylabel))
      xlabel_bottom = xlabel(1:nchar(xlabel))
      ylabel_bottom = ' '
      title_pltpar = title(1:nchar(title))
      comments(1) = 'Obs Band: 327 MHz '
      comments(2) = ' Proj Id: 17_078 '
      ncomments_passed = 2 

      call generate_pltpar_for_display(
     -               outfile,xlabel_left, ylabel_left,
     -               xlabel_bottom, ylabel_bottom, 
     -               title_pltpar,ncomments_passed,comments)


      stop 
      end

      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/generate_cfg_for_display.f'
      include '/usr/lib/subroutine_lib/generate_pltpar_for_display.f'
