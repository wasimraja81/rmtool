C This is an example Fortran program to educate one 
C on manipulating FITS files using the FITSIO tools. 
C 
C This code also intends at removing the inertia of 
C "inaction" that a beginner may suffer from due to 
C his/her non-familiarity of FITSIO usage. 
C
C          -- wasim raja, 16 July, 2009

      implicit none

      integer*4 maxdim,maxkeys, max_inbuff
      parameter(maxdim = 100, maxkeys=500, max_inbuff=8192)
      ! max_inbuff is the buffer size for the arrays meant 
      ! for reading the image/cubes. For images with a large 
      ! number of pixels, one may run out of memory if one 
      ! attempts to read the entire image into a single 
      ! array. We thus read and analyse the image/cube part 
      ! by part, taking care to prevent "out of memory" 
      ! situations by proper choice of max_inbuff parameter.
      integer*4 nchar
      character*172 infile
      integer*4 hdu_type
      integer*4 rwmode,iomode,exist
      integer*4 status
      integer*4 n_hdu           ! Total number of header units
      integer*4 current_hdu_num ! Specific header unit number
      real*4 version
      integer*4 blocksize
      integer*4 tmp_num
      integer*4 nrows, tfields, varidat
      character*16 ttype(5), tform(5), tunit(5)   ! FQ table ==> 5
      !character*16 ttype(19), tform(19), tunit(19) ! SU table ==> 19
      !character*16 ttype(12), tform(12), tunit(12) ! AN table ==> 12
      character*16 extname
      ! FOR FTGHPR:
      integer*4 bitpix, naxis, pcount,gcount, naxes(maxdim)
      logical simple,extend
      integer*4 n_keys, n_space
      integer*4 i, j
      character record*80
      character*16 tmpstr
      character*12 keyname(maxkeys)
      character*18 keyword
      character*1 yorn

      character*16 tmpstr1, comment1
      real*4 keyvalE
      ! Definitions related to reading IMAGE:
      integer group,dim1,nx,ny
      real buff_array(max_inbuff),nullval
      logical anyflg,ltemp
      integer fpixel,row

      integer dim2,nz      ! for CUBES
      integer*4 nbuffer, firstpix, npixels

      character*72 comment
      logical cube
      integer*4 data_precision
      cube = .false. ! by default, we assume that the file 
                     ! contains an image and NOT a cube...
                     ! We will automatically upgrade the 
                     ! status of cube later, after reading 
                     ! the file...

      ! IT IS IMPORTANT TO INITIALISE "STATUS" TO 
      ! ZERO, FOR ONE TO MAKE SUBROUTINE CALLS 
      ! SUCCESSFUL... IF PASSED WITH A POSITIVE  
      ! VALUE, THE SUBROUTINE WILL EXIT IMMEDIATELY. 
      ! FITSIO ALTERS THE STATUS's VALUE IF IT 
      ! ENCOUNTERS SOME SERIOUS ERRORS DURING THE 
      ! RUN OF THE PROGRAM, SO THAT ANY SUBSEQUENT 
      ! CALLS TO SUBROUTINES ARE NOT EXECUTED...
      ! THUS "STATUS" CAN BE USED TO MONITOR ERRORS
      ! AND MESSAGES RELATED TO THEM...
      status = 0

      write(*,*)"------------------------------------------- "
      write(*,*)" This code intends at removing the inertia " 
      write(*,*)" of 'inaction' that a beginner may suffer  " 
      write(*,*)" from due to his/her non-familiarity of " 
      write(*,*)" FITSIO usage."
      write(*,*)" You may however be wanting to jump into the" 
      write(*,*)" business of dissecting an image directly."
      write(*,*)" "
      write(*,*)" Choose as per your need whether you wish" 
      write(*,*)" to play around with various FITSIO calls" 
      write(*,*)" to SUBROUTINES and warm yourself up with " 
      write(*,*)" FITSIO usage, or jump straight to the image" 
      write(*,*)" reading section..."
      write(*,*)" "
      write(*,*)"Need warm-up? (y/n)"
      read(*,*)yorn
      write(*,*)" "

      if (yorn.eq.'n'.or.yorn.eq.'N')then
              write(*,*)"--------------------------- "
              write(*,*)"Looks like you chose to go  "
              write(*,*)"straight into business! We  "
              write(*,*)"shall take you to the IMAGE "
              write(*,*)"section..."
              write(*,*)" "
              write(*,*)"--------------------------- "
              goto 2999
      endif

      write(*,*)'Input FITS FILE: '
      read(*,'(a)')infile
      write(*,*)'chosen file: ',infile

      write(*,*)'read-write mode?'
      read(*,*)rwmode

      blocksize = 1

      write(*,*)"Do you wish to know the FITS  "
      write(*,*)"version you are using? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1111
      else
              goto 1112
      endif
      ! CHECK THE VERSION OF THE FITSIO BEING USED:
1111  continue
      call FTVERS(version)
      write(*,*)"---------------------------"
      write(*,*)"FITSIO VERSION: ",version
      write(*,*)" "
      write(*,*)"---------------------------"
1112  continue
      ! CHECK IF THE FILE TO BE OPENED EXISTS:
      !call FTEXIST(infile,exist,status)
      !write(*,*)"---------------------------"
      !write(*,*)"Current subroutine: FTEXIST"
      !write(*,*)"EXIST CODE: ",exist
      !write(*,*)"Here, STATUS = ",status
      !write(*,*)" "
      !write(*,*)"---------------------------"

      ! OPEN THE FITS FILE
      call FTOPEN(21,infile,rwmode,blocksize,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTOPEN"
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"

      write(*,*)"Do you wish to count the total no.  "
      write(*,*)"of headers in the FITS file? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1113
      else
              goto 1114
      endif
1113  continue
      ! COUNT TOTAL NUMBER OF HDU's IN THE FITS FILE:
      call FTTHDU(21,n_hdu,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTTHDU"
      write(*,*)"TOTAL No. OF HDU: ",n_hdu
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"

1114  continue
      ! GET THE CURRENT HDU NUMBER:
      write(*,*)"Do you wish to know the current "
      write(*,*)"Header number? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1115
      else
              goto 1116
      endif
1115  continue
      call FTGHDN(21,current_hdu_num,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTGHDN"
      write(*,*)"CURRENT HDU Number: ",current_hdu_num
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
1116  continue

      ! MOVE TO A "SPECIFIC" HDU NUMBER:
      write(*,*)"Do you wish to move to a   "
      write(*,*)"specific Header unit? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1117
      else
              goto 1118
      endif
1117  continue
      write(*,*)'Specify the header number you wish to move to: '
      read(*,*)tmp_num
      write(*,*)' '
      call FTMAHD(21,tmp_num,hdu_type,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTMAHD"
      write(*,*)"CURRENT HDU type: ",hdu_type
      write(*,*)"+++++++++++++++++++++++++++++ "
      write(*,*)"FYI, HDUTYPE: 0 =>  IMAGE_HDU "
      write(*,*)"              1 =>  ASCII_TBL "
      write(*,*)"              2 => BINARY_TBL "
      write(*,*)"+++++++++++++++++++++++++++++ "
      write(*,*)" "
      call FTGHDN(21,current_hdu_num,status)
      write(*,*)"CURRENT HDU No. changed to: ",current_hdu_num
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
1118  continue

      ! CHECK THE R/W MODE BEING USED, i.e., 
      ! WHETHER YOU HAVE OPENED THE FITS FILE 
      ! IN READ-ONLY OR READ-WRITE MODE...
      write(*,*)"Do you wish to check the read/write  "
      write(*,*)"mode with which the FITS file has "
      write(*,*)"been opened? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1119
      else
              goto 1120
      endif
1119  continue
      call FTFLMD(21,iomode,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTFLMD"
      write(*,*)"STATUS = ",status
      write(*,*)"R/W MODE USED: ",iomode
      write(*,*)" "
      write(*,*)"---------------------------"
1120  continue
      ! READ THE HEADER INFO: 
      ! (One needs to know apriori the type/length of certain 
      ! strings and variables: e.g., tform, ttype etc.)
      ! 1) BINARY TABLE: 
      ! FTGHBN(unit,maxdim, > nrows,tfields,ttype,tform,tunit,extname,varidat,status)
      write(*,*)"Do you wish to get the header info  "
      write(*,*)"from binary table e.g., tfields, " 
      write(*,*)"ttype etc.? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1121
      else
              goto 1122
      endif
1121  continue
      call FTGHBN(21,maxdim,nrows,tfields,ttype,tform,
     -           tunit,extname,varidat,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTGHBN"
      write(*,*)"  NROWS: ",nrows
      write(*,*)"TFIELDS: ",tfields
      write(*,*)"  TTYPE: ",ttype
      write(*,*)"  TFORM: ",tform
      write(*,*)"  TUNIT: ",tunit
      write(*,*)"EXTNAME: ",extname
      write(*,*)"VARIDAT: ",varidat
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
1122  continue
      ! 2) PRIMARY HDU (ASCII TABLE)
      ! FTGHPR(unit,maxdim, > simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
      !--------------------------------------------------
      write(*,*)"Do you wish to get the primary header "
      write(*,*)"info of ASCII table, e.g, SIMPLE, BITPIX " 
      write(*,*)"NAXIS etc.? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1123
      else
              goto 1124
      endif
1123  continue
      ! You however need to go to the 1ST header unit:
        call FTMAHD(21,1,hdu_type,status)
      !--------------------------------------------------
      call FTGHPR(21,maxdim, simple,bitpix,naxis,naxes,pcount,gcount,
     -            extend,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTGHPR"
      write(*,*)"(reading PRIMARY HDU...)"
      write(*,*)"  SIMPLE: ",simple
      write(*,*)"  BITPIX: ",bitpix
      write(*,*)"   NAXIS: ",naxis
      write(*,*)"   NAXES: ",(naxes(i),i=1,naxis)
      write(*,*)"  EXTEND: ",extend
      write(*,*)"  PCOUNT: ",pcount
      write(*,*)"  GCOUNT: ",gcount
      write(*,*)" STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
1124  continue
      ! FIND THE NUMBER OF KEYWORDS IN THE INPUT TABLE HEADER:
      !--------------------------------------------------
      ! SUBROUTINE: FTGHSP
      !
      ! It was found that FTGHSP GETS the KEYWORDS corresponding 
      ! to the CURRENT HDU only! We will write to a file all the 
      ! KEYWORDS corresponding to all existing HDUs...
      !--------------------------------------------------

      write(*,*)"Do you wish to get the number of keywords "
      write(*,*)"in ALL the HDUs and write them into a text" 
      write(*,*)"file? (y/n) "
      read(*,*)yorn
      if (yorn.eq.'y'.or.yorn.eq.'Y') then
              goto 1125
      else
              goto 1126
      endif
1125  continue
      open(unit=31,file='keywords.txt',status='unknown',err=9999)
      call FTTHDU(21,n_hdu,status)
      do i = 1,n_hdu
         call FTMAHD(21,i,hdu_type,status)
         call FTGHSP(21,n_keys,n_space,status)
         write(*,*)"---------------------------"
         write(*,*)"Current subroutine: FTGHSP"
         write(*,*)"   Current HDUTYPE: ",hdu_type
         write(*,*)"Number of KEYWORDS: ",n_keys
         write(*,*)"  Number of SPACES: ",n_space
         write(*,*)" "
         write(*,*)"---------------------------"

         write(31,*)"---------------------------"
         write(31,*)"       Current HDU: ",i
         write(31,*)"   Current HDUTYPE: ",hdu_type
         write(31,*)"Number of KEYWORDS: ",n_keys
         write(31,*)"  Number of SPACES: ",n_space
         write(31,*)" "
         write(*,*)"The KEYWORDS has been "
         write(*,*)"written into the file "
         write(*,*)"       'keywords.txt' "
         write(*,*)"---------------------------"
         ! GET ALL THE KEYWORDS FROM THE CURRENT HDU:
         write(*,*)"Current subroutine: FTGHREC"
         do j = 1,n_keys
            call FTGREC(21,j,record,status)
            write(31,*)"KEYWORD(",j,"): ",record(1:nchar(record))
         enddo
      enddo
      close(31)  ! Closing the keywords.txt file
1126  continue

      ! We can also get the VALUE of a KEYWORD using 
      ! the subroutine "FTGKEY"
      !call FTGKEY(21,"cdelt1",tmpstr1,comment1,status)
      call FTGKYE(21,"cdelt1",keyvalE,comment1,status)
      write(*,*)" "
      !write(*,*)"Current subroutine: FTGKEY"
      !write(*,*)"cdelt1: ",tmpstr1(1:nchar(tmpstr1))
      write(*,*)"Current subroutine: FTGKYE"
      write(*,*)"cdelt1: ",keyvalE,comment1(1:nchar(comment1))
      write(*,*)"STATUS = ",status
      write(*,*)" "
      ! CLOSE THE FITS FILE:
      call FTCLOS(21,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTCLOS"
      write(*,*)"STATUS = ",status
      write(*,*)" "
      write(*,*)"---------------------------"
      if (status .gt. 0)call printerror(status)

      write(*,*)'maxdim: ',maxdim
1999  goto 9999 ! If the previous section was called, 
                ! we shall ignore the image-reading 
                ! section...
                ! Hence go to END of program label.

2999  continue  ! You will be here if you choose to 
                ! play with IMAGE/CUBE etc.

!-------------------------------------------------------------------
      ! Get the file containing the Image/Cube:
      !
      write(*,*)'---------------------------'
      write(*,*)'Input Image/Cube Fits File:'
      write(*,*)' '
      read(*,'(a)')infile

      write(*,*)'read-write mode?'
      read(*,*)rwmode


      blocksize = 1

      ! Open the Image/Cube Fits file:
      call FTOPEN(21,infile,rwmode,blocksize,status)

      ! Determine the data-type of the image (BITPIX value):
      ! Possible returned values are: 
      !   8 : unsigned byte, 
      !  16 : signed 2-byte integer, 
      !  32 : signed 4-byte integer, 
      !  64 : signed 8-byte integer, 
      ! -32 : real, and 
      ! -64 : double. 

      call FTGIDT(21,bitpix,status)
      !call FTGIET(21,bitpix,status)
      write(*,*)"---------------------------"
      write(*,*)" Data-type of image: ",bitpix
      write(*,*)" "
      write(*,*)"Interpret the values as: "
      write(*,*)"!   8(B) : unsigned byte," 
      write(*,*)"!  16(I) : signed 2-byte integer," 
      write(*,*)"!  32(J) : signed 4-byte integer," 
      write(*,*)"!  64(K) : signed 8-byte integer," 
      write(*,*)"! -32(E) : real, and" 
      write(*,*)"! -64(D) : double." 
      write(*,*)" "
      write(*,*)"---------------------------"

      ! Determine the dimension(NUMBER of AXES) in 
      ! the image:
      call FTGIDM(21,naxis,status)
      write(*,*)" Total number of axes :",naxis
      write(*,*)" "
      write(*,*)"---------------------------"

      if(naxis.gt.99)then
              write(*,*)"NAXIS is > 99"
              write(*,*)"Currently this code expects that the "
              write(*,*)"image dimensions does not exceed 99"
              write(*,*)" "
              write(*,*)"However if you are sure that your im-"
              write(*,*)"age has NAXIS > 99, then modify the "
              write(*,*)"source code and RECOMPILE it..."
              write(*,*)" "
              write(*,*)"To modify the code: "
              write(*,*)"Open the code in your favourite text-"
              write(*,*)"editor, and search for the section: "
              write(*,*)"        'ERROR_NAXIS99' "
              write(*,*)"The section will describe to you the "
              write(*,*)"necessary modifications that needs to"
              write(*,*)"be made..."
              write(*,*)" "
              write(*,*)"---------------------------"
              stop
      endif
      ! Determine the KEYWORD character "values":
      ! In this case, we wish to know the "names" 
      ! of the axes...
      tmpstr(1:) = '              '
      do i = 1,naxis
      !
! ERROR_NAXIS99:
      ! The following if-section takes care of upto 
      ! 99 naxis. I do not think any FITS file will have 
      ! naxis > 99 ever.
      ! If at all such an apparently crazy need arises, 
      ! append appropriate number of "else if" sections 
      ! and corresponding "formats" defined as in label 
      ! number 499 and 599.
         if (i.lt.10)then
                 write(tmpstr(1:),fmt=499)i
                 keyword = 'ctype'//tmpstr(1:nchar(tmpstr))
         else if(i.ge.10.and.i.lt.100)then
                 write(tmpstr(1:),fmt=599)i
                 keyword = 'ctype'//tmpstr(1:nchar(tmpstr))
         endif
         call FTGKEY(21,keyword, keyname(i),comment,status)
         write(*,*)keyword(1:nchar(keyword)),": ",
     -keyname(i)(1:nchar(keyname(i)))
      enddo
499   format (i1)
599   format (i2)

      ! Determine the SIZE along each dimension of the image:
      call FTGISZ(21,maxdim,naxes,status)

      write(*,*)"---------------------------------------"
      write(*,*)"AXIS NUM     AXIS NAME      AXIS SIZE"
      do i = 1,naxis
         write(*,*)"     ",i,")     ",keyname(i),"    ",naxes(i)
      enddo
      write(*,*)" "
      if (naxes(1).gt.1.and.naxes(2).gt.1.and.naxes(3).gt.1)then
              npixels = naxes(1)*naxes(2)*naxes(3)
              write(*,*)"You are dealing with a CUBE..."
              write(*,*)"NPIXELS = ",npixels
              cube = .true.
      else if (naxes(1).gt.1.and.naxes(2).gt.1.and.naxes(3).eq.1)then
              npixels = naxes(1)*naxes(2)
              write(*,*)"You are dealing with a 2D IMAGE..."
              write(*,*)"NPIXELS = ",npixels
              cube = .false.
      endif
      !------------------------------------------------------------
      ! Let's read an image now...
      ! Depending on the BITPIX value of the image as determined 
      ! above, we will need to call one of the following routines:
      ! FTG2DB(unit,group,nullval,dim1,nx,ny,array,anyflg,status)
      ! FTG2DI(unit,group,nullval,dim1,nx,ny,array,anyflg,status)
      ! FTG2DK(unit,group,nullval,dim1,nx,ny,array,anyflg,status)
      ! FTG2DE(unit,group,nullval,dim1,nx,ny,array,anyflg,status)
      ! FTG2DD(unit,group,nullval,dim1,nx,ny,array,anyflg,status)
      !
      ! FTG2DE reads a 2-d image of real values from the primary 
      ! array.
      !
      ! Data conversion and scaling will be performed if necessary
      ! (e.g, if the datatype of the FITS array is not the same
      ! as the array being read).
      
      ! unit    i  Fortran output unit number
      ! group   i  number of the data group, if any
      !            [each group of the primary array is a 
      !             row in the table, where the first 
      !             column contains the group parameters
      !             and the second column contains the 
      !             image itself]
      !             Perhaps it will be helpful to imagine the primary 
      !             array as an array of dimension: (ngroups x 2), 
      !             where each element of the array is also 
      !             an array -- element(j,1) = [group params]
      !                      -- element(j,2) = [image], and
      !                         j = 1:ngroups

      ! nullval r  undefined pixels will be set to this value (unless = 0)
      ! dim1    i  actual first dimension of ARRAY
      ! nx      i  size of the image in the x direction
      ! ny      i  size of the image in the y direction
      ! array   r  the array of values to be read
      ! anyflg  l  set to true if any of the image pixels were undefined
      ! status  i  returned error stataus
      !
      !-------------------------------------------------------------
      ! However for a large image, it may not be a good idea 
      ! to try to read the entire image into a single large array
      ! Hence we read the image only part-by-part, with the aid of 
      ! the max_inbuff parameter defined already. The routine used 
      ! here is one of : FTGPV[BIJKED]
      !-------------------------------------------------------------
      !

      !-------------------------------------------------------------
      ! Sometimes it may be profitable to write out a binary 
      ! image... The following section writes a binary 2D image 
      ! from the FITS image file. Each record consists of the 
      ! data along the NAXIS1 axis.
      if(.not.cube)then
               if(bitpix.eq.8)data_precision = 1
               if(bitpix.eq.16)data_precision = 2
               if(bitpix.eq.32)data_precision = 4
               if(bitpix.eq.64)data_precision = 8
               if(bitpix.eq.-32)data_precision = 4
               if(bitpix.eq.-64)data_precision = 8

              group = 1
              firstpix = 1
              nullval = -999.0
              nbuffer = naxes(1)
              write(*,*)"nbuffer derived from naxes(1) = ",nbuffer

              open(31,file='image.bin',status='unknown',
     -                form='unformatted',access='direct',
     -                recl=nbuffer*data_precision)

              do i = 1,naxes(2)
                 if(bitpix.eq.8)then
                         call FTGPVB(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 else if(bitpix.eq.16)then
                         call FTGPVI(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 else if(bitpix.eq.32)then
                         call FTGPVJ(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 else if(bitpix.eq.64)then
                         call FTGPVK(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 else if(bitpix.eq.-32)then
                         call FTGPVE(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 else if(bitpix.eq.-64)then
                         call FTGPVD(21,group,firstpix,nbuffer,nullval,
     -                                  buff_array,anyflg,status)
                 endif
                 write(31,rec=i)(buff_array(j),j=1,nbuffer)
                 firstpix = firstpix + nbuffer
              enddo
      endif
      close(31)
      !
      ! Finished writing BINARY image file
      ! -----------------------------------------------------------------

      open(31,file='image.dat',err=9999,status='unknown')
      group = 1
      firstpix = 1
      nullval = -999.0
      do while (npixels.gt.0)
         nbuffer = min(max_inbuff,npixels)
         if (bitpix.eq.8)then 
                 call FTGPVB(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         else if (bitpix.eq.16)then
                 call FTGPVI(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         else if (bitpix.eq.32)then
                 call FTGPVJ(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         else if (bitpix.eq.64)then
                 call FTGPVK(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         else if (bitpix.eq.-32)then
                 call FTGPVE(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         else if (bitpix.eq.-64)then
                 call FTGPVD(21,group,firstpix,nbuffer,nullval,
     -                        buff_array,anyflg,status)
         endif
         do i = 1,nbuffer
            write(31,*)buff_array(i) ! NAXIS1 values form the inner loop
         enddo

         npixels = npixels - nbuffer
         firstpix = firstpix + nbuffer
      enddo

      close(31)
      ! Let's try using FTG2DE:

      ! CLOSE THE FITS FILE:
      call FTCLOS(21,status)
      write(*,*)"---------------------------"
      write(*,*)"Current subroutine: FTCLOS "
      write(*,*)"STATUS = ",status
      write(*,*)" "
      if (status .gt. 0)call printerror(status)

      write(*,*)"---------------------------"
      write(*,*)'maxdim: ',maxdim
9999  continue      
      end


C *************************************************************************
      subroutine printerror(status)

C  This subroutine prints out the descriptive text corresponding to the
C  error status value and prints out the contents of the internal
C  error message stack generated by FITSIO whenever an error occurs.

      integer status
      character errtext*30,errmessage*80

C  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return

C  The FTGERR subroutine returns a descriptive 30-character text string that
C  corresponds to the integer error status number.  A complete list of all
C  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

C  FITSIO usually generates an internal stack of error messages whenever
C  an error occurs.  These messages provide much more information on the
C  cause of the problem than can be provided by the single integer error
C  status value.  The FTGMSG subroutine retrieves the oldest message from
C  the stack and shifts any remaining messages on the stack down one
C  position.  FTGMSG is called repeatedly until a blank message is
C  returned, which indicates that the stack is empty.  Each error message
C  may be up to 80 characters in length.  Another subroutine, called
C  FTCMSG, is available to simply clear the whole error message stack in
C  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
C *************************************************************************
      include 'nchar.f'
