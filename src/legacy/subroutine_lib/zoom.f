chelp+
      ! This subroutine is used for zooming in to 
      ! a particular section of a plot using the 
      ! cursor position. 
chelp-

      subroutine zoom(xarr,yarr,npts,xlabel,ylabel,title,plot_type)

      real*4         xarr(*), yarr(*) 
      integer*4      plot_type 
      character*120  xlabel, ylabel, title 

      real*4         xmin, xmax, ymin, ymax 
      character*1    tmpchar 
      real*4         x, y 
      logical        satisfied 

      satisfied = .false. 

      call pgbeg(0,'1/xs',1,1)
      call myplot1(xarr,yarr,npts,xlabel,ylabel,title,plot_type)
      do while (.not.satisfied)
        write(*,*)"'a' to zoom, 'b' to original, 'd' to stop... "
        read(*,'(a)')tmpchar
        if(tmpchar.eq.'d'.or.tmpchar.eq.'D')then
                satisfied = .true. 
        else if(tmpchar.eq.'a'.or.tmpchar.eq.'A')then 
                write(*,*)"Click to choose BLC"
                call pgcurs(x, y, tmpchar)
                xmin = x 
                ymin = y 
                write(*,*)"Click to choose TRC"
                call pgcurs(x, y, tmpchar)
                xmax = x 
                ymax = y 

                if(xmax.le.xmin.or.ymax.le.ymin)then
                        write(*,*)"Choose proper BLC/TRC..."
                        goto 123 
                endif
                ! Plot zoomed in section 
                call pgask(.false.)
                call pgenv(xmin,xmax,ymin,ymax,0,1)
                call pglabel(xlabel,ylabel,title) 
                if (plot_type.eq.1)then
                        call pgpt(npts,xarr,yarr,1)
                else if(plot_type.eq.2)then
                        call pgline(npts,xarr,yarr)
                endif
        else if(tmpchar.eq.'b'.or.tmpchar.eq.'B')then
                call pgask(.false.)
                call myplot1(xarr,yarr,npts,xlabel,ylabel,title,
     -                       plot_type)
        endif
123     continue
      enddo
      call pgend 
        
      end
