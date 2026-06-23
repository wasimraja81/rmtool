***********************************************************************
***** This is a general subroutine which uses pgplot's subroutines to 
***** plot a 2-dimensional array as an image.
***** Variables to be passed :
**    1. X-array              --  real*4
**    2. Y-array              --  real*4
**    3. No. of points        --  integer*4
**    4. xlabel,ylabel,title  --  character*80
**    5. plot_dev             --  character*120
**    6. mode                 --  character*120
**       ("default" or "interactive")
***********************************************************************



        subroutine pgplot_imag(plot_arr,xdim,ydim,i1,i2,j1,j2,xmin,xmax,
     -             ymin,ymax,xlabel,ylabel,title,plot_dev,force_val, 
     -             amin, amax, bw_plot)

        implicit none
        integer*4 xdim, ydim
        integer*4 i,j,k,l,m,n,i1,i2,j1,j2,i0,j0
        real*4    plot_arr(xdim,ydim),xmin,xmax,ymin,ymax,
     -            tr(6),amin,amax,xstep,ystep,rms,mean,
     -            temp_arr(xdim*ydim),
     -            temp_plot_arr(xdim,ydim)
        character xlabel*(*),ylabel*(*),title*(*)
        character*80 tmp_str
        character plot_dev*(*) 
        logical bw_plot, force_val 

        ! bw_plot = .false.
        do i = 1,xdim
           do j = 1,ydim
              temp_plot_arr(i,j) = 0.0
           end do
        end do

        m = 1
        do i = i1,i2
           do j = j1,j2
              temp_arr(m) = plot_arr(i,j)
              m = m + 1
           end do
        end do
        call meanrms_1(temp_arr,rms,mean,m-1)


        if(amin.ge.amax)then
                force_val = .false. 
        endif
        if(.not.force_val)then
                amin = plot_arr(1,1)
                amax = amin
                do i = i1,i2
                   i0 = i - i1 + 1
                   do j = j1,j2
                      j0 = j - j1 + 1
c                      if(plot_arr(i,j)-mean.lt.6.0*rms)then
                       if(plot_arr(i,j).ne.0.0)then
                         amin = min(amin,plot_arr(i,j))
                         amax = max(amax,plot_arr(i,j))
                         temp_plot_arr(i0,j0) = plot_arr(i,j)
                        end if
c                      end if
                   end do
                end do
        endif
c        write(*,*)' Max, Min, mean, rms : ',amax,amin,mean,rms
        xstep = (xmax - xmin)/real(i2 - i1 + 1)
        ystep = (ymax - ymin)/real(j2 - j1 + 1)
        tr(1) = xmin - xstep/2.0
        tr(2) = xstep
        tr(3) = 0.0
        tr(4) = ymin - ystep/2.0
        tr(5) = 0.0
        tr(6) = ystep
!        plot_dev = '/xs'
        call pgbegin(0,plot_dev,1,1)
        !call pgbegin(0,'/xd',1,1)
        call pgpap(0.0,1.0)
        call pgenv(xmin,xmax,ymin,ymax,0,0)
        call pglabel(xlabel,ylabel,title)

        call set_colours_1(bw_plot,amin,amax)

        i0 = i2 - i1 + 1
        j0 = j2 - j1 + 1
        call pgimag(temp_plot_arr,xdim,ydim,1,i0,1,j0,amin,amax,
     -tr)
        CALL PGWEDG('RI', 0.5, 4.0, amin, amax,' ')
        call pgend
        return
        end
*************************************************************************
** include files
c        include 'nchar_lin.f'
*************************************************************************

        SUBROUTINE MEANRMS_1(A,RMS,AMEAN,NP)
C       TO COMPUTE MEAN AND RMS OF 'A' BY EXCLUDING WHAT MAY BE
C       SOME CONTRIBUTION FROM INTERFERENCE
        DIMENSION A(*)
        ITER=0
C
101     ITER=ITER+1
        AMEAN=0.0
        AN=0.0
C
        DO I=1,NP
        IF(ITER.EQ.1)GO TO 1
        DIFF=ABS(A(I)-AMEAN0)
        IF(DIFF.LE.(4.*RMS))GO TO 1
        GO TO 2
1       AMEAN=A(I)+AMEAN
        AN=AN+1.
2       CONTINUE
        END DO
C
        RMS0=RMS
        if(an.gt.0)AMEAN=AMEAN/AN
C
        RMS=0.0
        AN=0.0
        DO I=1,NP
        DIFF=ABS(A(I)-AMEAN)
        IF(ITER.EQ.1)GO TO 11
        IF(DIFF.LE.(4.*RMS0))GO TO 11
        GO TO 12
11      RMS=RMS+DIFF*DIFF
        AN=AN+1.
12      CONTINUE
        END DO
C
        if(an.gt.0.0)RMS=SQRT(RMS/AN)
        AMEAN0=AMEAN

        IF(ITER.EQ.1)GO TO 101
        if(rms.eq.0.0)return
        IF(ABS((RMS0/RMS)-1.).GT.0.05)GO TO 101
        RETURN
        END
*************************************************************
      subroutine upcase_1(string)

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
*******************************************************************

c========================================================================
      subroutine set_colours_1(bw_plot,amin,amax)

      implicit none

      real*4 dmin,dmax,amin,amax
      logical*4 bw_plot,col_reverse
      integer*4 idx,jj,k,ibw,m,i0,ii,kk
      real*4 rgb(2,3),csr,csg,csb ,val1
      logical*4 col_file_found
      character*40 temp_string

            dmin = amin
            dmax = amax

            call pgqcol(idx,jj)
c            write(*,*)'max # col:',jj
c first two colour indices used for background & default (txt) colour
            call pgqinf('HARDCOPY',temp_string,k)
c            write(*,'(a,a)')'HARDCOPY_DEVICE ? ....',
c     -                      temp_string(1:10)
            if(temp_string(1:3).eq.'YES')then
              col_reverse = .true.
            else
              col_reverse = .false.
            end if
c Set up the colour scale
            col_file_found = .false.
            if(.not.bw_plot)then
              if((dmin*dmax).ge.0.0)then    ! we have zero at one end
                if(dmin.ge.0.0)then
                  open(unit=35,
     -               file='colour_+ve_only.dat',type='old',
     -               err=1100)
                  col_file_found = .true.
1100              if(.not.col_file_found)then
                    open(unit=35,
     -               file='colour_+ve_only.dat',type='unknown')
c==============================================================      
c==============================================================                    
                      write(35,*)'5,2,255'
                      write(35,*)'20'
                      write(35,*)'0.0,0.0,0.0'
                      write(35,*)'0.0,0.0,1.0'
                      write(35,*)'50'
                      write(35,*)'0.0,0.0,1.0'
                      write(35,*)'0.0,1.0,1.0'
                      write(35,*)'60'
                      write(35,*)'0.0,1.0,1.0'
                      write(35,*)'0.0,1.0,0.0'
                      write(35,*)'60'
                      write(35,*)'0.0,1.0,0.0'
                      write(35,*)'1.0,1.0,0.0'
                      write(35,*)'64'
                      write(35,*)'1.0,1.0,0.0'
                      write(35,*)'1.0,0.0,0.0'
                      write(35,*)' '
       write(35,*)' The first line contains no.of info-sets',
     -            ' (each 3-line)'
       write(35,*)' to follow,the offset colour index & ',
     -            ' the total no.'
       write(35,*)'of levels defined in the following sets.'
       write(35,*)' '
       write(35,*)'Each 3-line info has'
       write(35,*)' '
       write(35,*)'the no. of steps'
       write(35,*)'starting rgb values'
       write(35,*)'ending rgb values '
       write(35,*)' '
       write(35,*)' --- desh'
c==============================================================
c==============================================================                    
                    close(unit=35)
                    open(unit=35,
     -               file='colour_+ve_only.dat',type='old')
                  end if
                else    ! -ve range
                  open(unit=35,
     -               file='colour_-ve_only.dat',type='old',
     -               err=1101)
                  col_file_found = .true.
1101              if(.not.col_file_found)then
                    open(unit=35,
     -               file='colour_-ve_only.dat',type='unknown')
c==============================================================
c==============================================================                    
                      write(35,*)'8,2,255'
                      write(35,*)'9'
                      write(35,*)'0.5,0.0,0.1'
                      write(35,*)'0.7,0.0,0.0'
                      write(35,*)'25  '
                      write(35,*)'0.7,0.0,0.0'
                      write(35,*)'0.7,0.2,0.0'
                      write(35,*)'40'
                      write(35,*)'0.7,0.2,0.0'
                      write(35,*)'1.0,1.0,0.0'
                      write(35,*)'40'
                      write(35,*)'1.0,1.0,0.0'
                      write(35,*)'0.7,0.7,0.0'
                      write(35,*)'40'
                      write(35,*)'0.7,0.7,0.0'
                      write(35,*)'0.0,0.5,0.0'
                      write(35,*)'40'
                      write(35,*)'0.0,0.5,0.0'
                      write(35,*)'0.0,0.7,0.7'
                      write(35,*)'40'
                      write(35,*)'0.0,0.7,0.7'
                      write(35,*)'0.1,0.1,0.7'
                      write(35,*)'20'
                      write(35,*)'0.1,0.1,0.7'
                      write(35,*)'0.0,0.0,0.0'
                      write(35,*)' '
       write(35,*)' The first line contains no.of info-sets',
     -            ' (each 3-line)'
       write(35,*)' to follow,the offset colour index & ',
     -            ' the total no.'
       write(35,*)'of levels defined in the following sets.'
       write(35,*)' '
       write(35,*)'Each 3-line info has'
       write(35,*)' '
       write(35,*)'the no. of steps'
       write(35,*)'starting rgb values'
       write(35,*)'ending rgb values '
       write(35,*)' '
       write(35,*)' --- desh'
c==============================================================
c==============================================================                    
                    close(unit=35)
                    open(unit=35,
     -               file='colour_-ve_only.dat',type='old')
                  end if
                end if

              else     ! we have zero in the middle
c             make the range look symmetric
                val1 = dmax
                if(val1.lt.abs(dmin))val1 = abs(dmin)
                dmin = -val1
                dmax = val1
                open(unit=35,
     -               file='colour_+ve_to_-ve.dat',type='old',
     -               err=1102)
                col_file_found = .true.
1102            if(.not.col_file_found)then
                  open(unit=35,
     -               file='colour_+ve_to_-ve.dat',type='unknown')
c==============================================================
c==============================================================                    
                  write(35,*)'8,2,255'
                  write(35,*)'30'
                  write(35,*)'0.0,0.5,0.0'
                  write(35,*)'0.0,0.9,0.0'
                  write(35,*)'30'
                  write(35,*)'0.0,0.9,0.0'
                  write(35,*)'0.0,1.0,1.0'
                  write(35,*)'40'
                  write(35,*)'0.0,1.0,1.0'
                  write(35,*)'0.0,0.0,1.0'
                  write(35,*)'27'
                  write(35,*)'0.0,0.0,1.0'
                  write(35,*)'0.0,0.0,0.0'
                  write(35,*)'47'
                  write(35,*)'0.0,0.0,0.0'
                  write(35,*)'1.0,1.0,0.0'
                  write(35,*)'50'
                  write(35,*)'1.0,1.0,0.0'
                  write(35,*)'1.0,0.0,0.0'
                  write(35,*)'15'
                  write(35,*)'1.0,0.0,0.0'
                  write(35,*)'0.7,0.0,0.0'
                  write(35,*)'15'
                  write(35,*)'0.7,0.0,0.0'
                  write(35,*)'0.5,0.0,0.1'
                  write(35,*)' '
       write(35,*)' The first line contains no.of info-sets',
     -            ' (each 3-line)'
       write(35,*)' to follow,the offset colour index & ',
     -            ' the total no.'
       write(35,*)' of levels defined in the following sets.'
       write(35,*)' '
       write(35,*)' Each 3-line info has '
       write(35,*)' '
       write(35,*)' the no. of steps'
       write(35,*)' starting rgb values'
       write(35,*)' ending rgb values'
       write(35,*)' -----------------------------desh '
c==============================================================
c==============================================================                    
                  close(unit=35)
                  open(unit=35,
     -               file='colour_+ve_to_-ve.dat',type='old')
                end if
              end if
            else           ! gray scale
                open(unit=35,
     -               file='gray_scale.dat',type='old',
     -               err=1103)
                col_file_found = .true.
1103            if(.not.col_file_found)then
                  open(unit=35,
     -               file='gray_scale.dat',type='unknown')
c==============================================================                    
                  write(35,*)'2,2,255'
                  write(35,*)'128'
                  write(35,*)'0.0,0.0,0.0'
                  write(35,*)'0.5,0.5,0.5'
                  write(35,*)'127'
                  write(35,*)'0.5,0.5,0.5'
                  write(35,*)'1.0,1.0,1.0'
                  write(35,*)' '
       write(35,*)' The first line contains no.of info-sets',
     -            ' (each 3-line)'
       write(35,*)' to follow,the offset colour index & ',
     -            ' the total no.'
       write(35,*)' of levels defined in the following sets.'
       write(35,*)' '
       write(35,*)' Each 3-line info has '
       write(35,*)' '
       write(35,*)' the no. of steps'
       write(35,*)' starting rgb values'
       write(35,*)' ending rgb values'
       write(35,*)' -----------------------------desh '
c==============================================================                    
                  close(unit=35)
                  open(unit=35,
     -               file='gray_scale.dat',type='old')
              end if
            end if

            read(35,*,err=1111)ibw,ii,kk
            if(ii.lt.2)ii = 2
            k = ii
            do m=1,ibw
              read(35,*,err=1111)i0
              i0 = real(i0)*real(jj-2)/real(kk-ii)

              do idx=1,2     ! read the start and stop settings
                read(35,*,err=1111)
     -               rgb(idx,1),rgb(idx,2),rgb(idx,3)
                val1 = rgb(idx,1) + rgb(idx,2) + rgb(idx,3)
                if(val1.eq.0.0.or.val1.eq.3.0)then   !we have black or white
                  if(col_reverse)then
                    rgb(idx,1) = 1. - rgb(idx,1)
                    rgb(idx,2) = 1. - rgb(idx,2)
                    rgb(idx,3) = 1. - rgb(idx,3)
                  end if
                end if
              end do

              if(i0.eq.0)i0=1
              do idx=1,3
                rgb(2,idx) = (rgb(2,idx)-rgb(1,idx))/
     -                       real(i0)
              end do
              do idx=1,i0
                val1 = real(idx-1)
                csr = rgb(1,1) + val1*rgb(2,1)
                csg = rgb(1,2) + val1*rgb(2,2)
                csb = rgb(1,3) + val1*rgb(2,3)
                if(k.gt.1.and.k.le.jj)then
                  call pgscr(k,csr,csg,csb)
                  k = k + 1
                end if
              end do
            end do
            k = k - 1
            if(k.lt.jj)jj = k

            call pgscir(ii,jj)
            close(unit=35,err=1112)
      return
1111  write(*,*)' Error in reading colour .dat files '
      write(*,*)' in the present directory...'
      write(*,*)' Please copy  ....'
      write(*,*)' colour_+ve_only.dat'
      write(*,*)' colour_-ve_only.dat'
      write(*,*)' colour_+ve_to_-ve.dat'
      write(*,*)' gray_scale.dat'
      write(*,*)'================================'
      return
1112  close(unit=35,err=1113)
      return
1113  write(*,*)' Error during closing the colour-code file.'
      return
      end

***********************************************************************
