chelp+ 
      ! ----------------------------------------------
      ! This is a general routine to read fits binary
      ! tables. It can be easily modified for various 
      ! purposes. 
      ! In its present form, this code first displays 
      ! various keywords in the binary table, which 
      ! may be used by the user to create a configure 
      ! file for reading specific columns. 
      ! 
      !                      --wr, 22 Nov, 2011
      ! ------------------------------------------
chelp- 

      implicit none 

      integer*4      maxpts,maxcols 
      parameter      (maxpts=655360,maxcols=10)
      integer*4      nrows, ncols 
      integer*4      nchar 
      integer*4      status, blocksize, rwmode  
      integer*4      n_hdu, hdu_type 
      character*172  infile, outfile, cfgfile 
      character*72   comment, templine 
!      character*80   card 
      integer*4      nkeys, n_addkeys 
      integer*4      bitpix 
      integer*4      i, j, icol    
      character*32   kw_str, kval_str 
      integer*4      use_col 
      real*4         out_arr(maxpts,maxcols),col_arr(maxpts)
      integer*4      frow, felem, nelements 
      logical        anyflg 
      real* 4        null_val
!      integer*4      cam_1(maxpts), cam_2(maxpts), cam_3(maxpts) 
      character*1    junkchar 
      integer*4      colnum(maxcols), ncol_use 
      integer*4      doall 

      null_val = 0.0 

      status = 0 
      rwmode = 0 
      blocksize = 1 
      bitpix = -32  ! default, modify it using input file 
!      write(*,*)"Input file name (w/o path): "
!      write(*,*)"[Inp file expected in DATA/ area]"
!      read(*,*)infile

      if(iargc() .lt. 1)then
              write(*,*)"Usage: "
              write(*,*)"   read_fits_table_general <infile (w/o path)>"
              write(*,*)"   [infile expected in ../DATA/ area.]"
              write(*,*)" "
              write(*,*)" Additional argument 0/1 may be specified"
              write(*,*)" 0 -> display keywords only (no outfile)"
              write(*,*)"Quitting now..."
              stop
      else if(iargc().eq. 1)then
              call getarg(1,infile)
              infile = infile(1:nchar(infile))
              doall = 0 
      else if(iargc().eq. 2)then
              call getarg(1,infile)
              infile = infile(1:nchar(infile))

              call getarg(2,templine)
              read(templine,*)i
              if(i .le.0)then
                      doall = 0 
              else
                      doall = 1 
              endif
      endif
      write(*,*)"Fits File read: ",infile(1:nchar(infile))
      outfile = infile(1:nchar(infile))//".txt"

      cfgfile = '../CONFIG/'//infile(1:nchar(infile))//'.cfg'
      infile = '../DATA/'//infile(1:nchar(infile))

      call FTOPEN(21,infile,rwmode,blocksize,status) 
      if(status .ne. 0)then
              write(*,*)"Error opening fits file: ",
     -                   infile(1:nchar(infile))
              write(*,*)"Status : ",status
              stop
      endif

      ! ---------------------------------------------------
      ! Find out the number of Header Units in the file: 
      call FTTHDU(21,n_hdu,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTTHDU"
      write(*,*)"TOTAL No. OF HDU: ",n_hdu
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
      ! ---------------------------------------------------

      ! Before anything else, JUMP to the relevant HDU : 
      call FTMAHD(21,n_hdu,hdu_type,status)
      write(*,*)" "

      ! Count the number of keywords in the current header unit: 
      call FTGHSP(21, nkeys,n_addkeys,status)
      write(*,*)"Number of existing keywords: ",nkeys
      write(*,*)"Number of additional keywords: ",n_addkeys

      ! --------------------------------------------------------
      ! Display some information: 
      do i = 1,nkeys
          call FTGKYN(21,i, kw_str,kval_str,comment,status)
          kw_str = kw_str(1:nchar(kw_str))
          kval_str = kval_str(1:nchar(kval_str))
          write(*,*)kw_str(1:nchar(kw_str)),"  ",
     -              kval_str(1:nchar(kval_str))
      enddo
      ! Count the total number of rows and cols: 
      call FTGNRW(21,nrows,status)
      call FTGNCL(21,ncols,status)
      write(*,*)"nrows,ncols: ",nrows,",",ncols 
      if(doall .le. 0)then 
              stop 
      endif 
      ! --------------------------------------------------------


      ! Do the real Business: 
      ! 
      !-------------------------------------------------
      open(33,file=cfgfile,status='old',err=101)
      goto 102 
101   write(*,*)"Error opening cfgfile: ",cfgfile(1:nchar(cfgfile))
      stop
102   continue 
      read(33,*)junkchar ! comment line 
      read(33,*)junkchar ! comment line 
      read(33,*)ncol_use    ! number of columns to read 
      if(ncol_use .gt. maxcols)then
              write(*,*)"Error: ncol_use exceeds maxcol!!"
              stop
      endif
      read(33,*)(colnum(j),j=1,ncol_use)    ! number of columns to read 
      close(33) 

      do icol = 1,ncol_use 
         use_col = colnum(icol) ! Col number now 
         frow = 1 
         felem = 1 
         nelements = nrows 
         call FTGCVE(21,use_col,frow,felem,nelements,null_val, 
     -                  col_arr,anyflg,status) 

         do j = 1,nelements 
            out_arr(j,icol) = col_arr(j) 
         enddo 
      enddo 

      write(*,*)"--------------------------" 
      write(*,*)"nrows out file: ",nelements 
      write(*,*)"ncols out file: ",ncol_use 

      !------------------------------------------------- 
      ! TEST: 
      !
      write(*,*)"------------------------------"
      write(*,*)"sample outputs for columns:",(colnum(j),j = 1,ncol_use)
      write(*,*)" "
      do i = 1,10 
         write(*,*)(out_arr(i,j),j=1,ncol_use) 
      enddo 
      ! Write outfile: 
      open(33,file=outfile, status='unknown') 
      do i = 1,nrows 
         write(33,*)(out_arr(i,j),j=1,ncol_use) 
      enddo 
      write(*,*)"Outfile written to: ",outfile(1:nchar(outfile))
      close(33) 
      call FTCLOS(21,status) 


      end 

      include '/usr/lib/subroutine_lib/nchar.f'
