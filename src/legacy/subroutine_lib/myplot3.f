
      !--------------------------------------------------
      ! This is a basic Plotting subroutine intended for 
      ! plotting a Y vs. X keeping in mind the need to 
      ! overlay. For this purpose, pgenv is not called 
      ! inside this subroutine. The user is expected to 
      ! call pgenv from inside the calling program as 
      ! and when he/she feels like.

      !   -- wasim raja, 30 Apr, 2010
      !--------------------------------------------------


          subroutine myplot3(X,Y,npts,xlabel,ylabel,header,myplot_type)

          real*4 X(*), Y(*)
          integer*4 npts
          real*4 x_pts(npts),y_pts(npts), ymin,ymax,xmin,xmax
          !integer*4 max_dim
          !parameter (max_dim = 262144)
          !real*4 x_pts(max_dim),y_pts(max_dim), ymin,ymax,xmin,xmax
          integer*4 myplot_type
          !
          ! myplot_type = 1 --> plot points
          !             = 2 --> plot line
          !character*120  header, xlabel, ylabel
          character  header*(*) , xlabel*(*) , ylabel*(*) 


          !call pgbeg(0,'/xs',1,1)
          ! We shall call pgbeg in the calling program and 
          ! then call this subroutine as many number of times 
          ! as may be specified by pgbeg.
          call pgsch(2.0)


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
          call pgslw(1)
          call pglabel(xlabel,ylabel,header) 
          if (myplot_type.eq.1)then
                  !call PGPT (int N, float XPTS, float YPTS, int SYMBOL)
                  call pgpt(npts,x_pts,y_pts,1)
          else if(myplot_type.eq.2)then
                  call pgline(npts,x_pts,y_pts)
          endif

          return
          end

