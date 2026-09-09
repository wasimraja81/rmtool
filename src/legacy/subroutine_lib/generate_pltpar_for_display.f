chelp+
      !----------------------------------------------
      ! This routine generates a default pltpar file 
      ! required by the display routines given vital 
      ! parameters such as trc, blc, title, label etc. 
      ! The user may later manually modify this file 
      ! to suit his/her need. 
      !                         --wr, 17 Dec, 2011
      !----------------------------------------------
chelp-

      subroutine generate_pltpar_for_display(
     -               filename,xlabel_left, ylabel_left,
     -               xlabel_bottom, ylabel_bottom, 
     -               title,ncomments_passed,comments)

      implicit none 

      integer*4     nchar, iline, nlines_comment, ncomments_passed    
      character*220 filename, pltparfile 
      character*32  templine, comments(*)
      character*62  xlabel_left, ylabel_left, xlabel_bottom, 
     -              ylabel_bottom, title 



      pltparfile = filename(1:nchar(filename))//".pltpar"

      open(31,file=pltparfile,status='unknown')
      write(31,fmt=9)'/xs       ! soft_device_name'
      write(31,fmt=9)'/ps       ! hard_device_name'
      write(31,fmt=9)'1.        ! compress_fac ??'
      write(31,fmt=9)'4         ! twindows_used ??'
      write(31,fmt=9)'0.30, 0.87, 0.13, 0.30,'
      write(31,fmt=9)'0.13, 0.30, 0.30, 0.82,'
      write(31,fmt=9)'0.30, 0.87, 0.30, 0.82,'
      write(31,fmt=9)'0.76, 0.99, 0.02, 0.78'
      write(31,fmt=9)'0.0,0.0,0.0,                      xmin_disp'
      write(31,fmt=9)'0.0,0.0,0.0,                      xmax_disp'
      write(31,fmt=9)'0.0,0.0,0.0,                      ymin_disp'
      write(31,fmt=9)'0.0,0.0,0.0,                      ymax_disp'
      write(31,fmt=9)'0.0,0.0,0.1,               dmin,dmax,c_step'
      write(31,*)"  "
      write(31,fmt=9)xlabel_left(1:nchar(xlabel_left))
      write(31,fmt=9)title(1:nchar(title))
      write(31,fmt=9)xlabel_bottom(1:nchar(xlabel_bottom))
      write(31,*)"  "
      write(31,*)"  "
      write(31,fmt=9)ylabel_bottom(1:nchar(ylabel_bottom))
      write(31,fmt=9)ylabel_left(1:nchar(ylabel_left))
      write(31,*)"  "
      write(31,fmt=9)"ABINTS,0.0,0    ! xopt,xtick,nxsub | required f
     -or defining"
      write(31,fmt=9)"BCIMTV,0.0,0    ! yopt,ytick,nysub | the pgbox
     -for window 1"
      write(31,fmt=9)"BCIMTV,0.0,0    ! -ve n(x/y)sub will decide the
     - # of major ticks"
      write(31,fmt=9)"ABINTS,0.0,0    ! & a +ve number will decide th
     -e # of minor ticks"
      write(31,fmt=9)"BCV,0.0,0   "
      write(31,fmt=9)"BCV,0.0,0   "
      write(31,*)"  " 
      write(31,*)"  " 
      write(31,*)"  " 
      write(31,fmt=9)"-1,-1,-1 "
      write(31,fmt=9)"T"
      write(31,fmt=9)"1.,                   txt_disp "
      write(31,fmt=9)"0.5,                  txt_coord "
      write(31,fmt=9)"0.5,                  txt_fjust"
      nlines_comment = ncomments_passed + 1
      if (nlines_comment .gt. 1 .and. nlines_comment .lt.10)then
              write(templine,'(i1)')nlines_comment 
      else if (nlines_comment .ge. 10 .and. nlines_comment .lt.100)then
              write(templine,'(i2)')nlines_comment 
      else
              nlines_comment = 1
              write(templine,'(i1)')nlines_comment
      endif
      templine = templine(1:nchar(templine))//"          number lines"
      write(31,fmt=9)templine(1:nchar(templine))
      do iline = 1,nlines_comment-1
          templine = comments(iline)
          write(31,fmt=9)templine(1:nchar(templine))
      enddo
      templine = '---------------------------'
      write(31,fmt=9)templine(1:nchar(templine))
      write(31,fmt=9)"0,             for multi line plots"
      write(31,fmt=9)"0.33,          yblow factor"
      write(31,fmt=9)"2              line_width_index for PQ figures"
      write(31,fmt=9)"0.8,1.,3.0     Normal_char_size,Label_char_size
     -, label_displacement"
      write(31,fmt=9)"!  end of .par"


      close(31)
!9     format (0x,a)
9     format (a)
      return
      end
