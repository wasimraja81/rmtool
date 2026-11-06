chelp+
      !----------------------------------------------
      ! This routine generates a default config file 
      ! required by the display routines given vital 
      ! parameters such as nrec, nchan, nsets etc. 
      ! The user may later manually modify this file 
      ! to suit his/her need. 
      !                         --wr, 02 Jun, 2011
      !----------------------------------------------
chelp-

      subroutine generate_cfg_for_display(
     -               path,filename,nset,nchan,
     -               data_precision,begx,endx,nptsx,
     -               begy,endy)

      implicit none 

      integer*4     nset, nchan, data_precision 
      integer*4     nchar, nptsx  
      real*4        begx, endx, begy, endy 
      character*220 filename, pltparfile, path, cfgfile  

      character     templine*82


      pltparfile = filename(1:nchar(filename))//".pltpar"
      cfgfile = filename(1:nchar(filename))//".cfg"

      open(31,file=cfgfile,status='unknown')

      write(31,fmt=9)"# Config file for use with my_display routines"
      write(31,fmt=7)path(1:nchar(path))," ; PATH to data "
      write(31,fmt=7)filename(1:nchar(filename))," ; infile: name of the
     - input file (w/o path)"
      write(31,*)nset," ; n_set of data in file (eg., for I,Q,U & V,
     - n_set = 4)"
      write(31,fmt=9)"1,0, 1,0 ; beg_rec, nrec: read beg_rec to beg_rec+
     -nrecrecs for nrec > 0; if nrec < 1, read all records beginning fro
     -m beg_rec"
      write(31,fmt=9)"1 ; plot_set: 0<=plot_set<=n_set |0->Plot All data
     -. "
      write(31,fmt=9)"0,0,0,0 ; collapse (x,y) to: 0->mean, 1->rms, 2->m
     -ax, 3->min, 4->median; Then perform operation on (x,y) using: 2->L
     -OG10 " 
      write(31,fmt=9)"8192 ; nrec_per_pass: total number of records YOU
     - want to read in a single pass"
      write(31,*)nchan, ";reclen: length of a single record (number 
     -of data-points)" 
      write(31,*)data_precision, " ; data_precision"
      write(31,fmt=7)pltparfile(1:nchar(pltparfile))," ; Vital plot para
     -ms stored in this file"
      write(31,fmt=9)"# ========================================"
      write(31,fmt=9)"# Math operations to be performed on data:"
      write(31,fmt=9)"n  ; perform math-operation on any 2 sets: useful 
     -for computing AMP, PHA etc of Polarisation data"
      write(31,fmt=9)"AMP  ; OPTYPE: currently allowed values are: AMP, 
     -PHA"
      write(31,fmt=9)"2 3  ; Data Sets to perform operation on. Currentl
     -y only 2 sets allowed"
      write(31,fmt=9)"# ========================================"
      write(31,fmt=9)"# Plot Related Variables (can be ported into MENU)
     -:"
      write(31,fmt=9)"0, 0 ; remove baseline from (ROW,COLUMN)? : 0-> No
     -, else-> Yes"
      write(31,fmt=9)"3, 3    ; order_harmonic (ROW,COLUMN)"
      write(31,fmt=9)"3, 3    ; order_poly (ROW,COLUMN)" 
      write(31,fmt=9)"n       ; hard_copy (ps_req) ? "
      write(31,fmt=9)"n       ; forceval" 
      write(31,fmt=9)"0.75    ; formin: min val to be forced for image"
      write(31,fmt=9)"1.5     ; formax: max val to be forced for image"
      write(31,fmt=9)"n       ; use default plot dev as in .par?"
      write(31,fmt=9)"3       ; DISPLAY type (2:c; 1,0 or -1:line)"
      write(31,fmt=9)"6       ; n_contur  for contour plots"
      write(31,fmt=9)"0.3     ; y_fraction for each  line plot"
      write(31,fmt=9)"/xs     ; SOFT_PLOT dev choice (e.g. /xd)"
      write(31,fmt=9)"AUTO    ; ps-filename"
      write(31,fmt=9)"default ; PS_TAG"
      write(31,fmt=9)"../PLOTS/   ; PS_DIR"
      write(31,fmt=9)"FRANSPOSE   ; Additional Request"
      write(31,fmt=9)"#=================================================
     -=="
      write(31,fmt=9)"# Some more flexibility for book-keeping" 
      write(31,fmt=9)"2           ; nlines below to be used in the displ
     -ay" 
      write(31,fmt=9)"Mean has been computed ; info line 1 will appear i
     -n the RHS of disp-dev"
      write(31,fmt=9)"for collapsed axes ;info line 2"
      write(31,fmt=9)"#=================================================
     -=="
      write(31,fmt=9)"# Provision for scaling the x and y-axis values of
     - line plots:"
      write(31,fmt=9)"1 1         ; scale axes (x,y)? : 0 -> No, else Ye
     -s"
      write(templine,*)begx," ",endx," ",nptsx,"; begx, endx, nptsx"
      write(31,*)templine(1:nchar(templine))
      write(templine,*)begy," ",endy," ",nchan,"; begy, endy, nptsy"
      write(31,*)templine(1:nchar(templine))
      write(31,9)"y       ; bandflip? "
      write(31,fmt=9)"#=================================================
     -=="

!7     format (0x,a,a)
!8     format (0x,i2,a)
!9     format (0x,a)
7     format (a,a)
8     format (i2,a)
9     format (a)
      close(31)

      return
      end
