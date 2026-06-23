chelp+
C This subroutine writes out a one-dimensional 
C array , "A" into a file, "outfile" having 
C n_row rows and n_col columns. The length of A 
C is expected to be: n_row X n_col
C If it is less than that, the code will crib. 
C If the length of A is greater than n_row X n_col, 
C elements of A beyond index n_row X n_col will be 
C ignored in the file!

C     Author: Wasim Raja
C     Date: 06-03-2008

C        LAST MODIFICATION: 
C        Additionally, now, you ought to specify the outfile 
C        type -- ie., ASCII or BINARY using the variable 
C        "aorb" -- 'a' will force ASCII outfile, while 
C        'b' will force binary. Any other character will 
C        force an ASCII output file.
C        The default (and presently only the default exists)
C        precision is 4 Bytes.
C      
C        --wasim raja
C          Date: 09-12-2009
chelp-
          subroutine vec2array_file(A,n_row,n_col,outfile,aorb)
         integer*4 max_dim
         parameter (max_dim = 16777216)
         real*4 A(*)
         integer*4 i,j,k
         integer*4 n_row, n_col
         character*72 outfile
         character*1 aorb
         integer*4 buffsize

         !k = 0
         !do i = 1,n_row
         !   do j = 1,n_col
         !      k = k + 1
         !      tmp(k) = A(i,j)
         !   enddo
         !enddo

         if(aorb.ne.'b'.and.aorb.ne.'B')then
                 open(911, file = outfile,status='new')
                 k = 0
                 do i = 1,n_row
                    j = k+1
                    write(911,*)(A(k),k = j,j+n_col-1)         
                    !write(*,*)(A(k),k = j,j+n_col-1)         
                    k = k-1
                 enddo
                 close(911)
         else
                 buffsize = 4*n_col
                 open(911, file=outfile,form='unformatted',
     -                     access='direct', recl=buffsize,
     -                     status='new')
                 k = 0
                 do i = 1,n_row
                    j = k+1
                    write(911,rec=i)(A(k),k = j,j+n_col-1)         
                    !write(*,*)(A(k),k = j,j+n_col-1)         
                    k = k-1
                 enddo
                 close(911)
         endif
         end


