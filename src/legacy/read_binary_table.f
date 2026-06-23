chelp+ 
      ! ------------------------------------------
      ! This routine is intended at familiarizing 
      ! reading of fits binary tables. 
      !                      --wr, 22 Nov, 2011
      ! ------------------------------------------
chelp- 

      implicit none 

      integer*4      maxpts 
      parameter      (maxpts=655360)
      integer*4      nrows, ncols 
      integer*4      nrow_single 
      integer*4      nchar 
      integer*4      status, blocksize, rwmode  
      integer*4      n_hdu, hdu_type 
      character*172  infile, outfile_1, outfile_2, outfile_3, outfile_4 
      character*72   comment  
!      character*80   card 
      integer*4      nkeys, n_addkeys 
      integer*4      bitpix 
      integer*4      i, iread, nbytes_row  
      character*32   kw_str, kval_str 
      integer*4      use_col 
      integer*2      col_arr(maxpts)
      integer*4      frow, felem, nelements 
      logical        anyflg 
      real* 4        null_val
      integer*4      icam_1, icam_2, icam_3, icam_22 
!      integer*4      cam_1(maxpts), cam_2(maxpts), cam_3(maxpts) 
      integer*4      int_arr(64)

      null_val = 0.0 

      status = 0 
      rwmode = 0 
      blocksize = 1 
      bitpix = -32  ! default, modify it using input file 
      write(*,*)"Input file name (w/o path): "
      write(*,*)"[Inp file expected in DATA/ area]"
      read(*,*)infile

      infile = infile(1:nchar(infile))
      write(*,*)"Fits File read: ",infile(1:nchar(infile))
      outfile_1 = infile(1:nchar(infile))//"_cam1"
      outfile_2 = infile(1:nchar(infile))//"_cam2"
      outfile_3 = infile(1:nchar(infile))//"_cam3"
      outfile_4 = infile(1:nchar(infile))//"_cam22"

      call FTOPEN(21,infile,rwmode,blocksize,status) 
      write(*,*)"In file opening Status : ",status

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
      ! --------------------------------------------------------


      ! Do the real Business: 
      ! Open the output FITS files to be written: 
      call FTINIT(22,outfile_1,blocksize,status) 
      if(status.ne.0)then
              write(*,*)"Error opening unit 22 with Status : ",status
              stop
      endif
      call FTINIT(23,outfile_2,blocksize,status) 
      if(status.ne.0)then
              write(*,*)"Error opening unit 23 with Status : ",status
              stop
      endif
      call FTINIT(24,outfile_3,blocksize,status) 
      if(status.ne.0)then
              write(*,*)"Error opening unit 24 with Status : ",status
              stop
      endif
      call FTINIT(25,outfile_4,blocksize,status) 
      if(status.ne.0)then
              write(*,*)"Error opening unit 25 with Status : ",status
              stop
      endif

      ! Copy all the headers to the output files: 
      write(*,*)"n_hdu: ",n_hdu 
      do i = 1,n_hdu
         call FTMAHD(21,i,hdu_type,status)
         write(*,*)"---------------------------"
         write(*,*)"CURRENT HDU type: ",hdu_type
         call FTCPHD(21, 22,status)
         call FTCPHD(21, 23,status)
         call FTCPHD(21, 24,status)
         call FTCPHD(21, 25,status)
      enddo

      write(*,*)"+++++++++++++++++++++++++++++ "
      write(*,*)"FYI, HDUTYPE: 0 =>  IMAGE_HDU "
      write(*,*)"              1 =>  ASCII_TBL "
      write(*,*)"              2 => BINARY_TBL "
      write(*,*)"+++++++++++++++++++++++++++++ "
      write(*,*)" "

      ! ------------------------------------------------
      ! We wish to filter out the data based on the 
      ! value of a particular column. 
      ! Example: For the Cas-A light curve from RXTE-ASM 
      !          (xa_casa_d1_ch2.lc) we wish to separate 
      !          data for the various cameras. 
      ! Define the relevant column to be used for 
      ! filtering: 
      ! 
      !-------------------------------------------------
      ! Incorporate this section in a parfile later: 
      use_col = 5 ! Col number containing camera information 
      frow = 1 
      felem = 1  
      nelements = nrows 
      call FTGCVI(21,use_col,frow,felem,nelements,null_val,
     -                  col_arr,anyflg,status)

      icam_1 = 0 
      icam_2 = 0  
      icam_3 = 0  
      icam_22 = 0  


      ! Find out the length (in BYTES) of a single row: 
      call FTGKYJ(21,"NAXIS1",nbytes_row,comment,status)
      write(*,*)"Length of a single row (in bytes) in Input file: ",
     -            nbytes_row 
      do i = 1,nelements 
         if(col_arr(i) .eq. 1)then
                 icam_1 = icam_1 + 1 
                 iread = i 
                 ! Read the i-th row and write it to outfile_1: 
                 !!call FTGTBB(unit,frow,startchar,nchars, > array,status)
                 !!call FTPTBB(unit,frow,startchar,nchars,array, > status)
                 call FTGTBB(21,iread,1,nbytes_row, int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error reading row from infile "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_1: ",icam_1
                         write(*,*)"status: ",status
                         stop 
                 endif
                 call FTPTBB(22,icam_1,1,nbytes_row,int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error writing row "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_1: ",icam_1
                         write(*,*)"status: ",status
                         stop 
                 endif

         else if(col_arr(i) .eq. 2)then
                 icam_2 = icam_2 + 1 
                 iread = i 
                 ! Read the i-th row and write it to outfile_1: 
                 call FTGTBB(21,iread,1,nbytes_row, int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error reading row from infile "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_2: ",icam_2
                         write(*,*)"status: ",status
                         stop 
                 endif
                 call FTPTBB(23,icam_2,1,nbytes_row,int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error writing row "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_2: ",icam_2
                         write(*,*)"status: ",status
                         stop 
                 endif

         else if(col_arr(i) .eq. 3)then
                 icam_3 = icam_3 + 1 
                 iread = i 
                 ! Read the i-th row and write it to outfile_1: 
                 call FTGTBB(21,iread,1,nbytes_row, int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error reading row from infile "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_3: ",icam_3
                         write(*,*)"status: ",status
                         stop 
                 endif
                 call FTPTBB(24,icam_3,1,nbytes_row,int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error writing row "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_3: ",icam_3
                         write(*,*)"status: ",status
                         stop 
                 endif

         else
                 icam_22 = icam_22 + 1 
                 iread = i 
                 ! Read the i-th row and write it to outfile_1: 
                 call FTGTBB(21,iread,1,nbytes_row, int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error reading row from infile "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_22: ",icam_22
                         write(*,*)"status: ",status
                         stop 
                 endif
                 call FTPTBB(25,icam_22,1,nbytes_row,int_arr,status)
                 if(status .ne. 0)then
                         write(*,*)"Error writing row "
                         write(*,*)"original Row num: ",i
                         write(*,*)"Row num in outfile_22: ",icam_22
                         write(*,*)"status: ",status
                         stop 
                 endif
         endif
      enddo 

      ! Modify the NAXIS2 keywords in the outfiles: 
      call FTMKYJ(22,"NAXIS2",icam_1,"&",status)
      call FTMKYJ(23,"NAXIS2",icam_2,"&",status)
      call FTMKYJ(24,"NAXIS2",icam_3,"&",status)
      call FTMKYJ(25,"NAXIS2",icam_22,"&",status)


      write(*,*)"       nrows 1: ",icam_1  
      write(*,*)"       nrows 2: ",icam_2  
      write(*,*)"       nrows 3: ",icam_3  
      write(*,*)"      nrows 22: ",icam_22  
      write(*,*)"                --------"
      write(*,*)"nrows inp file: ",nelements 




      ! In case we incorporate the following scheme: 
      ! Scheme 2: Read optimal number of rows in one go, 
      !           and then based on the filter, redistribute 
      !           them into the outfiles. 
      ! Determine the optimal number of rows to read in one go: 
      call FTGRSZ(21, nrow_single,status)
      write(*,*)"Optimal # rows to read at a time: ",nrow_single 

!      call FTGCL(21,use_col,frow,felem,nelements,col_arr,status) 
      ! SBIJKEDCM

      !-------------------------------------------------




      call FTCLOS(21,status) 
      call FTCLOS(22,status) 
      call FTCLOS(23,status) 
      call FTCLOS(24,status) 


      end 

      include '/usr/lib/subroutine_lib/nchar.f'
