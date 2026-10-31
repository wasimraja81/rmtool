
      !--------------------------------------------------
      ! This is a basic Plotting subroutine intended for 
      ! plotting a Y vs. X

      !   -- wasim raja, 01 Aug, 2009
      !   -- added a cfg file provision for various 
      !      user requirements.
      !--------------------------------------------------


          subroutine myplot1_v2(X,Y,npts,xlabel,ylabel,header,cfg_file)

          real*4 X(*), Y(*)
          integer*4 npts
          real*4 x_pts(npts),y_pts(npts), ymin,ymax,xmin,xmax
          !integer*4 max_dim
          !parameter (max_dim = 262144)
          !real*4 x_pts(max_dim),y_pts(max_dim), ymin,ymax,xmin,xmax
          integer*4 myplot_type
          character*120 cfg_file
          !
          ! myplot_type = 1 --> plot points
          !             = 2 --> plot line
          !character*120  header, xlabel, ylabel
          character  header*(*) , xlabel*(*) , ylabel*(*) 
          logical fixed_scale

          open(29,file=cfg_file,status='old',err=101)
          goto 102
101       write(*,*)' cfg file: ',cfg_file(1:nchar(cfg_file))
          write(*,*)'NOT FOUND!! '
          write(*,*)'Quitting Now...'
          stop

102       continue
          read(29,*)myplot_type
          read(29,*)fixed_scale
          read(29,*)xmin,xmax,ymin,ymax
          close(29)

          !call pgbeg(0,'/xs',1,1)
          ! We shall call pgbeg in the calling program and 
          ! then call this subroutine as many number of times 
          ! as may be specified by pgbeg.
          call pgsch(2.0)


          if(.not.fixed_scale)then
                  xmin = 1.0e18
                  xmax = -1.0e18
                  ymin = 1.0e18
                  ymax = -1.0e18
                  do j = 1,npts
                       x_pts(j) = X(j)
                       y_pts(j) = Y(j)
                       if(xmin.gt.x_pts(j)) xmin = x_pts(j)
                       if(xmax.lt.x_pts(j)) xmax = x_pts(j)
                       if(ymin.gt.y_pts(j)) ymin = y_pts(j)
                       if(ymax.lt.y_pts(j)) ymax = y_pts(j)
                  end do
                  if(ymax.eq.ymin)then
                          ymax = ymax + 1.0
                          ymin = ymin - 1.0
                  endif
                  if(xmax.eq.xmin)then
                          xmax = xmax + 1.0
                          xmin = xmin - 1.0
                  endif
          else 
                  do j = 1,npts
                       x_pts(j) = X(j)
                       y_pts(j) = Y(j)
                  end do
                  if(ymax.eq.ymin)then
                          ymax = ymax + 1.0
                          ymin = ymin - 1.0
                  endif
                  if(xmax.eq.xmin)then
                          xmax = xmax + 1.0
                          xmin = xmin - 1.0
                  endif
          endif
          call pgslw(1)
          !call PGENV (float XMIN, float XMAX, float YMIN, float YMAX, integer JUST, integer AXIS)
          call pgenv(xmin,xmax,ymin,ymax,0,1)
          call pglabel(xlabel,ylabel,header) 
          if (myplot_type.eq.1)then
                  !call PGPT (int N, float XPTS, float YPTS, int SYMBOL)
                  call pgpt(npts,x_pts,y_pts,1)
          else if(myplot_type.eq.2)then
                  call pgline(npts,x_pts,y_pts)
          endif

          return
          end

