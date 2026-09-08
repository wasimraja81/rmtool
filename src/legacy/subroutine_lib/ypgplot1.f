*************************************************************************
ccc  ###################################################################
cccc      subroutine pgplot_2d(x,y,npts_in,opt_str)
ccc  ===================================================================
ccc  
ccc  A generic subroutine to plot 2-dimensional line/curve/scatter-plots
ccc
ccc     Subroutine Arguments : 
ccc           x() : X-axis data array (real*4)
ccc           y() : Y-axis data array (real*4)
ccc          npts : Number of points in above two array (integer*4)
ccc     opt_str() : Array of Character variables (strings) (character*(*))
ccc
ccc       NOTE : Keep at the least 32 variables in the array opt_str()
ccc              (size of the opt_str() array : at least 32).
ccc              For more information on how to use this array to plot
ccc              the data in different 'kinds' and 'styles', read below
ccc              the description of each of the array elements.....
ccc
ccc     Include "-L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11"
ccc     in the command line (or make file) while compiling the program
ccc     which includes & uses this subroutine (you might have to change
ccc     the above paths appropriately).
ccc
ccc  ===================================================================
ccc
ccc  ^^^^^^^^^^^^^^^^^^^^^^^^^ Author's note ^^^^^^^^^^^^^^^^^^^^^^^^^^^
ccc     This subroutine has been tested successfully for various options 
ccc     (and their combinations) below. However, in case it is found NOT
ccc     to be working as expected (from description below) for some 
ccc     combinations below, it would be appreciated to receive a bug-report
ccc     for the same (and it will certainly help in making this subroutine
ccc     reach faster towards its final version).
ccc                                     - Yogesh Maan
ccc                                       (yogesh@rri.res.in)
ccc                                       October 2009
ccc  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ccc
ccc  ===================================================================
ccc  What the various elements of the character-variable array opt_str()
ccc  can be used for : 
cc
cc      ------------------------------------------------------------------
cc      opt_str(1) =    Plot-device     (character string)
cc      ------------------------------------------------------------------
cc      opt_str(2) =    X-label         (character string)
cc      ------------------------------------------------------------------
cc      opt_str(3) =    Y-label         (character string)
cc      ------------------------------------------------------------------
cc      opt_str(4) =    Plot-title      (character string)
cc      ------------------------------------------------------------------
cc      opt_str(5) =    Plot-type       (character string)
cc              'LINE'  --> line-plot/curve  (DEFAULT)
cc              'SCAT'  --> scatter-plot
cc                      a.) 'NMARKS='$1--> Advanced scatter-plot, where
cc                                        several graph markers (not
cc                                        all of them are same), will
cc                                        be used to draw the scatter-
cc                                        plot. The code-numbers of the
cc                                        symbols to be plotted are expected
cc                                        to be present in y-array itself,
cc                                        followed the data. $1 specifies 
cc                                        the number of code-values.
cc                      a.) 'NCOLS='$2--> Advanced scatter-plot, where
cc                                        several color indices (not
cc                                        all of them are same), will
cc                                        be used to draw the scatter-
cc                                        plot. The code-numbers of the
cc                                        colors to be used are expected
cc                                        to be present in x-array itself,
cc                                        followed the data. $2 specifies 
cc                                        the number of code-values.
cc                      c.) otherwise     default scatter-plot, with default
cc                                        (or that specified in opt_str(8))
cc                                        graph marker.
cc              'HIST'  a.) 'BINNED'  --> Histogram of binned input data
cc                      b.) 'RAW'     --> unbinned input data, to be
cc                                        binned first
cc                      c.) otherwise --> Histogram of unbinned data
cc                              (in (b) & (c), x-array will not be used
cc                               and it's presence would be sort of dummy)
cc      
cc     ****It may also be an appropriate combination of the above****

cc      ------------------------------------------------------------------
cc      opt_str(6) = In case multiple plots are needed in the same
cc                   plot-window (multiple data sets are to be plotted
cc                   in the same window), it can be specified whether
cc                   the pg-window is required to be kept open for next
cc                   data-set....       (character string)
cc              'FIRST' -->     first-call of the routine, 
cc                              This will force to use "pgbegin" and then 
cc                              NOT to use "pgend" (to keep it open for 
cc                              next data-set)
cc              'LAST'  -->     last-call of the routine
cc                              This will assume that the plot-window is
cc                              already open (so will not use pgbegin),
cc                              and will force to use "pgend"
cc              'KEEP'  -->     keep the plot-window open
cc                              This will assume that the plot-window
cc                              is already open, and has to keep it open.
cc                              (so will NOT use either of "pgbegin" and
cc                               "pgend").
cc
cc              ****ONLY ONE of the above OPTIONS can be used at a time****
cc              Absence of all the above assumes the default case, where
cc              only the present data set is to be plotted (both FIRST_call
cc              and LAST_call are .true.)
cc
cc            'XNORM'/'YNORM'   Presence of 'XNORM'/'YNORM' along with 'LAST' 
cc                         or 'KEEP' would force to normalize the X/Y-range 
cc                         of present data to the min-max range plotted in 
cc                         the last call. If the present X/Y-range is to
cc                         be limited to some %age of previous min-max 
cc                         range, then use the syntax 'XNORM='$1, 
cc                         'YNORM='$2, where $1 & $2 specifies the
cc                         respective %ages.
cc
cc              'EQUALIZE'        Force the min/max along the two axis
cc                              to be equal to the minimum/maximum
cc                              along the two directions....square
cc                              viewing surface in plotting units.
cc

cc      ------------------------------------------------------------------
cc      opt_str(7) =    Force min-max   (4 float numbers)
cc                      (flaot xmin, xmax, ymin, ymax)
cc              In case partial control is required, fill the
cc              remaining (of above 4) variables by magic_val(-9999.0)....
cc              ..these (remaining) parameters will be appropriately 
cc              calculated in the subroutine.

cc      ------------------------------------------------------------------
cc      opt_str(8), 
cc      opt_str(9) & 
cc      opt_str(10) = EXTRA 'CONTROLS'
cc              (any of the following inputs can go to either
cc               of the above three strings, in any order)
cc              a.) 'SYMBOL='$1 -->  Integer*4; spcifies the 'symbol'-sign 
cc                                      to be used for the points to be 
cc                                      plotted, in case ofSCATTER-plot 
cc                                      (default symbol=1)
cc              b.) 'CENTER='$2 -->  Real*4;  if ($2 > 0.0), then in the
cc                                      'BINNED'-histogram data the
cc                                      x-array positions will be assumed
cc                                      to be the center-positions of the
cc                                      bins, otherwise these will be
cc                                      considered as the lower-edges 
cc                                      of the varisou bins
cc              c.) 'NBINS='$3  -->  Integer*4;  In case of 'UNBINNED'-
cc                                      histogram data, $3 specifies the
cc                                      number of bins to be used. (permitted
cc                                      range : 10 <= nbins <=200)
cc              d.) 'FILL'      -->  Plots the 'UNBINNED'-histogram data
cc                                      in fill-area-style
cc                                      ***default is EMPTY-BIN-STYLE***
cc              e.) 'LINE'      -->  Plots the 'UNBINNED'-histogram data
cc                                      in simple line-plot syle
cc              f.) 'GRID'      -->  Draw grid-lines at major increments of
cc                                      the coordinates (NOT applicable in
cc                                      case of 'UNBINNED'-histogram-data
cc              g.) 'COLOR='$4  -->  Plots the line/scatter/histogram plot 
cc                                      using the color-index specified by
cc                                      the value "$4". This color-setting
cc                                      will NOT be applied to the axes/box
cc                                      /labels/title or anything written
cc                                      outside the box.
cc              h.) various line-syles
cc                  'SOLID'             normal full-line (default)
cc                  'DASHED'            dashed-line
cc                  'DOTTED'            dotted-line
cc                  'DASH-DOT-DASH'     dash-dot-dash-dot line
cc                  'DASH-DOT-DOT'      dash-dot-dot-dot-dash line
cc               or use 'LSTYLE='$5     where range of $5 is (1-5)
cc
cc              i.) 'LWIDTH='           sets the "line-width" of the
cc                                      lines/curves to be plotted.
cc              j.) 'SWAP'              swap the axis, just before 
cc                                      axis-normalizing and plotting 




cc      ------------------------------------------------------------------
cc      opt_str() = 
cc      opt_str() = 
ccc  ===================================================================
ccc
ccc
ccc
ccc  ^^^^^^^^^^^^^^^^^^^^^^^^^ Author's note ^^^^^^^^^^^^^^^^^^^^^^^^^^^
cc        Latest modifications/successful-tests done in November 2009.
cc                                              - Yogesh Maan
ccc  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ccc  ###################################################################


      subroutine pgplot_2d(x,y,npts_in,opt_str)
 
       implicit none
       integer*4  max_pts
       parameter  (max_pts=4096*1024*4)
       integer*4  i,j,k,l,m,n,i1,npts_in,npts,symbol,nchar_local2d
       integer*4  nbins,pgflag,ci,lstyle,nsyms,ncols,lwidth,
     -            sarray(max_pts),carray(max_pts)
       real*4     xmin,xmax,ymin,ymax,x(*),y(*),
     -            old_xmin,old_xmax,old_ymin,old_ymax,
     -            xmin_local,xmax_local,ymin_local,ymax_local,
     -            x_percen,y_percen
       real*4     magic_val,temp1,temp2,temp3,temp4
       real*4     x_plot(max_pts),y_plot(max_pts),temp_arr(max_pts)
       character(len=*) opt_str(*)
       character  xlabel*(*),ylabel*(*),title*(*),plot_dev*(*)
       character*128  templine1,templine2,templine3
       character*512  cat_line
       character  plot_type*64
       logical    line_plot,scat_plot,hist_plot,binned,
     -            first_call,last_call,center,grid_on,
     -            xnorm,ynorm,tobe_binned,multi_marks,
     -            multi_cols,swap_axis,equalize
       logical    log_scale,linear_scale

        npts = npts_in


ccc     Default settings
        magic_val = -9999.0
        symbol    = 1
        nbins     = 100
        pgflag    = 0   !1           ! plot in current window
        line_plot = .true.
        scat_plot = .false.
        hist_plot = .false.
        binned    = .false.
        tobe_binned=.false.
        linear_scale = .true.
        log_scale = .false.
        first_call = .true.
        last_call = .true.
        center    = .true.
        grid_on   = .false.
        xnorm     = .false.
        ynorm     = .false.
        ci        = 1           ! default color-index
        x_percen  = 100.0
        y_percen  = 100.0
        lstyle    = 1
        multi_marks=.false.     ! different graph-markers
        multi_cols= .false.     ! different color-indices
        nsyms = 1      ! no. of graph-marker codes (symbols) provided
        swap_axis = .false.
        lwidth = 1              ! line-width
        equalize = .false.


cc      ============================================================
cc      Check if the plot-device, labels,title and plot-type
cc      are provided
        !...........................................................
        if(nchar_local2d(opt_str(1)).gt.0)then    ! Plot-device
          templine1(1:) = opt_str(1)
          plot_dev(1:) = templine1(1:nchar_local2d(templine1))
        else
          plot_dev(1:) = '/xs'                  ! default
        end if
        !...........................................................
        if(nchar_local2d(opt_str(2)).gt.0)then    ! Label for X-axis
          templine1(1:) = opt_str(2)
          xlabel(1:) = templine1(1:nchar_local2d(templine1))
        else
          xlabel(1:) = 'X-axis'                 ! default
        end if
        !...........................................................
        if(nchar_local2d(opt_str(3)).gt.0)then    ! Label for Y-axis
          templine1(1:) = opt_str(3)
          ylabel(1:) = templine1(1:nchar_local2d(templine1))
        else
          ylabel(1:) = 'Y-axis'                 ! default
        end if
        !...........................................................
        if(nchar_local2d(opt_str(4)).gt.0)then    ! Plot-title
          templine1(1:) = opt_str(4)
          title(1:) = templine1(1:nchar_local2d(templine1))
        else
          title(1:) = ''                        ! default
        end if
        !...........................................................
        if(nchar_local2d(opt_str(5)).gt.0)then    ! Plot-type
          templine1(1:) = opt_str(5)
          plot_type(1:) = templine1(1:nchar_local2d(templine1))
        else
          plot_type(1:) = 'line'                ! default
        end if
        !...........................................................
        if(nchar_local2d(opt_str(6)).gt.0)then    ! keep plot-win open ?
          templine1(1:) = opt_str(6)
          call upcase_local2d(templine1)
          if(index(templine1,'FIRST').gt.0)then
             first_call = .true.
             last_call  = .false.
          else if(index(templine1,'LAST').gt.0)then
             first_call = .false.
             last_call  = .true.
          else if(index(templine1,'KEEP').gt.0)then
             first_call = .false.
             last_call  = .false.
          end if
          if(index(templine1,'NORM').gt.0)then
            if(index(templine1,'X').gt.0) xnorm = .true.
            if(index(templine1,'Y').gt.0) ynorm = .true.
          else
            xnorm = .false.
            ynorm = .false.
          end if
          i1 = index(templine1,'XNORM=')
          if(i1.gt.0)then
             read(templine1(i1+6:),*,err=710,end=710)x_percen
          else
             x_percen = 100.0
          end if
710       if(x_percen.le.0) x_percen = 100.0
          i1 = index(templine1,'YNORM=')
          if(i1.gt.0)then
             read(templine1(i1+6:),*,err=711,end=711)y_percen
          else
             y_percen = 100.0
          end if
711       if(y_percen.le.0) y_percen = 100.0
          if(index(templine1,'EQUALIZE').gt.0) equalize = .true.
        end if
        old_xmin = magic_val
        old_xmax = magic_val
        old_ymin = magic_val
        old_ymax = magic_val
        if(xnorm.or.ynorm)then
          templine1(1:) = opt_str(16)
          read(templine1(1:),*,err=712,end=712)old_xmin,old_xmax,
     -                                        old_ymin,old_ymax
        end if
712     continue
        !...........................................................
cc          -----------------------------------------------------
cc          find out the min-max in the arrays
            xmin = x(1)
            xmax = xmin
            ymin = y(1)
            ymax = ymin
            do i = 1,npts
               if(x(i).ne.magic_val)then
                 if(xmin.gt.x(i)) xmin = x(i)
                 if(xmax.lt.x(i)) xmax = x(i)
               end if
               if(y(i).ne.magic_val)then
                 if(ymin.gt.y(i)) ymin = y(i)
                 if(ymax.lt.y(i)) ymax = y(i)
               end if
            end do
            if(first_call)then
              !xmin = xmin - (xmax-xmin)/20.0d0
              !xmax = xmax + (xmax-xmin)/20.0d0
              temp1 = (ymax-ymin)/20.0d0
              ymin = ymin - temp1
              ymax = ymax + temp1
            end if
cc          -----------------------------------------------------
        !...........................................................
        temp1 = magic_val
        temp2 = magic_val
        temp3 = magic_val
        temp4 = magic_val
        if(nchar_local2d(opt_str(7)).gt.0)then    ! min-max to be forced
          templine1(1:) = opt_str(7)
          read(templine1(1:),*,err=713,end=713)temp1,temp2,temp3,temp4
        end if
713     continue
        if(first_call)then
          if(temp1.ne.magic_val) xmin = temp1
          if(temp2.ne.magic_val) xmax = temp2
          if(temp3.ne.magic_val) ymin = temp3
          if(temp4.ne.magic_val) ymax = temp4
        end if
        templine1(1:) = ''
        templine2(1:) = ''
        templine3(1:) = ''
        !...........................................................
        if(nchar_local2d(opt_str(8)).gt.0)then    ! extra inputs
          templine1(1:) = opt_str(8)
        end if
        if(nchar_local2d(opt_str(9)).gt.0)then    ! extra inputs
          templine2(1:) = opt_str(9)
        end if
        if(nchar_local2d(opt_str(10)).gt.0)then    ! extra inputs
          templine3(1:) = opt_str(10)
        end if
        cat_line(1:) = templine1(1:nchar_local2d(templine1))
     -               //templine2(1:nchar_local2d(templine2))
     -               //templine3(1:nchar_local2d(templine3))
        call upcase_local2d(cat_line)
        if(nchar_local2d(cat_line).gt.0)then
          i1 = index(cat_line,'SYMBOL=')        ! type of points (scatter-plot)
          if(i1.gt.0)then
            read(cat_line(i1+7:),*,err=714,end=714)symbol
          else
            symbol = 1
          end if
714       i1 = index(cat_line,'CENTER=')        ! position of bins (histogram)
          if(i1.gt.0)then
            read(cat_line(i1+7:),*,err=715,end=715)temp1
            if(temp1.lt.0.0) center = .false.
          else
            center = .true.
          end if
715       i1 = index(cat_line,'NBINS=')         ! no. of bins
          if(i1.gt.0)then
            read(cat_line(i1+6:),*,err=716,end=716)nbins
          else
            nbins = 100
          end if
716       continue
          i1 = index(cat_line,'COLOR=')         ! color-index
          if(i1.gt.0)then
            read(cat_line(i1+6:),*,err=717,end=717)ci
            ci = ci - int(ci/16)*16
            if(ci.eq.0)then
              write(*,*)'WARNING:Plotting-COLOR set to Background color'
            end if
          else
            ci = 1
          end if
717       if(ci.lt.0) ci = 1
          if(index(cat_line,'FILL').gt.0)then   ! fill-area style (histogram)
            pgflag = 2    !3
          else
            pgflag = 0    !1
          end if
          if(index(cat_line,'LINE').gt.0)then   ! line-style (histogram)
            pgflag = 4   !5
          else
            pgflag = 0   !1
          end if
          if(index(cat_line,'GRID').gt.0)then   ! grid-ON
            grid_on = .true.
          else
            grid_on = .false.
          end if
          if(index(cat_line,'SWAP').gt.0)then   ! swap the axis ?
            swap_axis = .true.
          else
            swap_axis = .false.
          end if
          if(index(cat_line,'DASH-DOT-DOT').gt.0) lstyle = 5
          if(index(cat_line,'DOTTED').gt.0) lstyle = 4
          if(index(cat_line,'DASH-DOT-DASH').gt.0) lstyle = 3
          if(index(cat_line,'DASHED').gt.0) lstyle = 2
          i1 = index(cat_line,'LSTYLE=')
          if(i1.gt.0)then
             read(cat_line(i1+7:),*,err=718,end=718)lstyle
          end if
718       if(lstyle.lt.1.or.lstyle.gt.5) lstyle = 1
          i1 = index(cat_line,'LWIDTH=')
          if(i1.gt.0)then
             read(cat_line(i1+7:),*,err=719,end=719)lwidth
          end if
719       if(lwidth.le.0.or.lwidth.gt.200) lwidth=1
        end if



        call upcase_local2d(plot_type)
        if(index(plot_type,'LINE').gt.0)then    ! line-plot
           line_plot = .true.
        else
           line_plot = .false.
        end if
        if(index(plot_type,'SCAT').gt.0)then    ! scatter-plot
           scat_plot = .true.
           i1 = index(plot_type,'NMARKS=')
           if(i1.gt.0)then
              multi_marks = .true.
              nsyms = npts
              read(plot_type(i1+7:),*,err=781,end=781)nsyms
           else
              nsyms = 1
           end if
781        if(nsyms.lt.1) nsyms = 1
           if(nsyms.gt.npts) nsyms = npts
           i1 = index(plot_type,'NCOLS=')
           if(i1.gt.0)then
              multi_cols = .true.
              ncols = npts
              read(plot_type(i1+6:),*,err=782,end=782)ncols
           else
              ncols = 1
           end if
782        if(ncols.lt.1) ncols = 1
           if(ncols.gt.npts) ncols = npts
        else
           scat_plot = .false.
        end if
        if(index(plot_type,'HIST').gt.0)then    ! Histogram
           hist_plot = .true.
           if(index(plot_type,'BIN').gt.0)then  ! binned-data
              binned = .true.
              nbins  = npts
           else if(index(plot_type,'RAW').gt.0)then ! unbinned data, but
              tobe_binned = .true.                  ! to be binned
              binned      = .true.
           else                                 ! unbinned-data
              binned = .false.
           end if
        else
           hist_plot = .false.
        end if
        if(hist_plot.and..not.binned)then ! set the others to .false.
                                ! since it would not be compatible with
                                ! any other kind of plotting
           line_plot = .false.
           scat_plot = .false.
        end if
cc      ============================================================
       
cc      -----------------------------------------------------
ccc     Do the common things
        if(xmax.eq.xmin) xmax = xmin + 1.0
        if(ymax.eq.ymin) ymax = ymin + 1.0
        do i = 1,npts
           x_plot(i) = x(i)
           y_plot(i) = y(i)
        end do
        if(tobe_binned)then     ! data to be binned.....so do it here
           call bin_data_local2d(y,npts,nbins,center,x_plot,y_plot)
           npts = nbins
           xmin = x_plot(1) - (x_plot(2)-x_plot(1))/2.0
           xmax = x_plot(nbins) + (x_plot(2)-x_plot(1))/2.0
           ymin = y_plot(1)
           ymax = ymin
           do i = 1,nbins
              if(y_plot(i).ne.magic_val)then
                if(ymin.gt.y_plot(i)) ymin = y_plot(i)
                if(ymax.lt.y_plot(i)) ymax = y_plot(i)
              end if
           end do
           !temp1 = (xmax-xmin)/20.0d0
           !xmin = xmin - temp1
           !xmax = xmax + temp1
           temp1 = (ymax-ymin)/20.0d0
           ymin = ymin - temp1
           ymax = ymax + temp1
           binned = .true.
        end if
        if(hist_plot.and..not.binned)then
          if(nbins.lt.10) nbins = 10
          if(nbins.gt.200) nbins = 200
        end if


        if(equalize)then  ! equalize the min-max along the two axis
          if(xmin.lt.ymin)then
             ymin = xmin
          else
             xmin = ymin
          end if
          if(xmax.gt.ymax)then
            ymax = xmax
          else
            xmax = ymax
          end if
          if(-xmin.gt.xmax)then
            xmax = -xmin
            ymax = xmax
          else
            xmin = -xmax
            ymin = xmin
          end if
        end if

        if(swap_axis)then
           do i = 1,npts
              temp_arr(i) = x_plot(i)
              x_plot(i)   = y_plot(i)
              y_plot(i)   = temp_arr(i)
           end do
           temp1 = xmin
           xmin  = ymin
           ymin  = temp1
           temp1 = xmax
           xmax  = ymax
           ymax  = temp1
        end if
        if(old_xmin.eq.old_xmax
     -     .or.old_xmin.eq.magic_val
     -     .or.old_xmax.eq.magic_val) xnorm = .false.
        if(old_ymin.eq.old_ymax
     -     .or.old_ymin.eq.magic_val
     -     .or.old_ymax.eq.magic_val) ynorm = .false.
        if(xnorm)then
             old_xmax = old_xmin +
     -                       (old_xmax-old_xmin)*x_percen/100.0
        end if
        if(ynorm)then
             temp1 = (old_ymax - old_ymin)/22.0d0
             old_ymax = old_ymax - temp1  ! correction for 5% extension
             old_ymin = old_ymin + temp1  ! on either side
             old_ymax = old_ymin +
     -                       (old_ymax-old_ymin)*y_percen/100.0
        end if
        if(xnorm)then
           do i = 1,npts
              x_plot(i) = (x_plot(i)-xmin)*
     -                   (old_xmax-old_xmin)/(xmax-xmin) + old_xmin
           end do
        end if
        if(ynorm)then
           do i = 1,npts
              y_plot(i) = (y_plot(i)-ymin)*
     -                   (old_ymax-old_ymin)/(ymax-ymin) + old_ymin
           end do
        end if




        if(scat_plot.and.multi_marks)then
           do i = 1,nsyms
              sarray(i) = 1
              sarray(i) = nint(y(i+npts))
              if(sarray(i).lt.1) sarray(i) = 1
           end do
           do i = nsyms+1,npts
              sarray(i) = sarray(1)
           end do
        end if
        if(scat_plot.and.multi_cols)then
           do i = 1,ncols
              carray(i) = 1
              carray(i) = nint(x(i+npts))
              if(carray(i).lt.1) carray(i) = 1
           end do
           do i = ncols+1,npts
              carray(i) = carray(1)
           end do
        end if

cc      -----------------------------------------------------



c      do the initial settings...............
       if(first_call)then
         call pgbegin(0,plot_dev,1,1)
         if(hist_plot.and..not.binned)then
           ! no need to call "pgenv"-subroutine, will be called 
           ! automatically by "pghist"
         else
           if(grid_on)then
             call pgenv(real(xmin),real(xmax),real(ymin),real(ymax),0,2)
           else
             call pgenv(real(xmin),real(xmax),real(ymin),real(ymax),0,0)
           end if
         end if
       end if
       !if(.not.first_call) write(*,*)'***',xmin,xmax,ymin,ymax

c      now the real plotting.............
       call pgsci(ci)
       call pgsls(lstyle)
       call pgslw(lwidth)
       if(line_plot) call pgline(npts,x_plot,y_plot)
       if(scat_plot)then
         if(multi_cols)then
           call pgsci(1)
           do i = 1,npts
              x_plot(1) = x_plot(i)
              y_plot(1) = y_plot(i)
              call pgsci(carray(i))
              if(multi_marks) symbol = sarray(i)
              call pgpt(1,x_plot,y_plot,symbol)
           end do
         else if(multi_marks)then
           call pgpnts(npts,x_plot,y_plot,sarray,nsyms)
         else
           call pgpt(npts,x_plot,y_plot,symbol)
         end if
       end if
       if(hist_plot)then
         if(binned)then
           call pgbin(nbins,x_plot,y_plot,center)
         else
           write(*,*)'nbins : ',nbins
           call pghist(npts,y_plot,ymin,ymax,nbins,pgflag)
         end if
       end if
       call pgsls(1)
       call pgsci(1)
       call pgslw(1)

c      label, any text to be written in the plot....and finish
       call pglabel(xlabel,ylabel,title)
       if(last_call) call pgend


c      send some info back, for future use.....
       if(first_call)then
         if(hist_plot.and..not.binned)then
           write(templine1(1:),*)'0.0 1.0 0.0 1.0'
         else
           write(templine1(1:),*)xmin,xmax,ymin,ymax
         end if
         opt_str(16) = templine1(1:)
       end if


       return
       end

*************************************************************************



cc###################################################################
cc###################################################################
cc############ LOCAL FUNCTIONS/SUBROUTINES ##########################


****************************************************************
cc      A subroutine to bin (for a given 
cc      number of bins; nbins) the data in dat-array. X-array
cc      returns the central (or lower-edge, depending upon the
cc      value of the logical parameter "center")-value of each
cc      of the bins, and the Y-array returns the number of 
cc      occurances for each of the bins.
cc                                              -yogesh
cc                                               October 2009
****************************************************************

        subroutine bin_data_local2d(dat,npts,nbins,center,xarr,yarr)

        implicit none
        integer*4 max_pts
        parameter (max_pts=4096*1024*4)
        integer*4 i,m,npts,nbins
        real*4    dat(*),xarr(max_pts),yarr(max_pts),
     -            ymin_local,ymax_local,temp1,temp2
        logical   center



           ymin_local = dat(1)
           ymax_local = ymin_local
           do i = 1,npts
              if(ymin_local.gt.dat(i)) ymin_local = dat(i)
              if(ymax_local.lt.dat(i)) ymax_local = dat(i)
           end do
           temp1 = (ymax_local - ymin_local)/real(nbins) ! bin-interval
           do i = 1,nbins
              yarr(i) = 0.0
              if(center)then
                xarr(i) = ymin_local + temp1*real(i-1) + temp1/2.0
              else
                xarr(i) = ymin_local + temp1*real(i-1)
              end if
           end do
           temp2 = 0.0
           do i = 1,npts
              m = int((dat(i)-ymin_local)/temp1) + 1
              if(m.gt.nbins)then        ! this will happen only when
                                        ! dat(i) = ymax_local
                m = nbins
                temp2 = temp2 + 1.0
              end if
              yarr(m) = yarr(m) + 1.0
           end do
c_check           write(*,*)'****sanity check ... temp2 = ',temp2

           return
           end
****************************************************************





       Integer function nchar_local2d(string)
C  Routine to count the number of characters in the
C  input string. Looks for the last occurrence of 
C  non-(null, blank or tab character)
      Implicit none
      integer*4 i,ipos
      character*(*)  string
      character blank,tab,null,c

      blank=' '
        tab=char(9)
        null=char(0)
      ipos = 0
      i      = len(string)
      do while (i.gt.0.and.ipos.eq.0)
         c = string(i:i)
         if (c.ne.blank.and.c.ne.tab.and.c.ne.null) ipos = i
         i = i - 1
      end do
      nchar_local2d = ipos
      return
      end

      subroutine upcase_local2d(string)
c   to convert strings to Upper case characters
      implicit  none
      character*(*) string,temp_string*80
      integer*4     nchar,
     -              ichar,
     -              istring,
     -              i
      character     char,blank,tab,null,c
      integer*4     ipos
*** let us make this subroutine stand-alone
*** (no dependence on any other subroutine)

      temp_string = string
      blank=' '
        tab=char(9)
        null=char(0)
      ipos = 0
      i    = len(temp_string)
      do while (i.gt.0.and.ipos.eq.0)
         c = temp_string(i:i)
         if (c.ne.blank.and.c.ne.tab.and.c.ne.null) ipos = i
         i = i - 1
      end do
      nchar = ipos
***  ************************
      do i=1,nchar
        istring = ichar(string(i:i))
        if(istring.gt.96.and.istring.lt.123) then   ! lower case
          string(i:i) = char(istring-32)
        end if
      end do
      return
      end
