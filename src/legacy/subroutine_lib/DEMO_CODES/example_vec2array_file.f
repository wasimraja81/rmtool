c        This is a demo Code
c        Date: 13-02-2008
c        Wasim Raja

         implicit none
         real*4 P(4,4),A(16)
         integer*4 n_row, n_col
         integer*4 i, j, k
         character*72 outfile

         do i = 1,72
            outfile(i:i) = ' '
         enddo

         n_row = 4
         n_col = 3

         k = 0
         do i = 1,n_row
            do j = 1,n_col
               k = k + 1
               P(i,j) = k
            enddo
         enddo
         k = 0
         do i = 1,n_row
            do j = 1,n_col
               k = k + 1
               A(k) = P(i,j)
            enddo
         enddo
         outfile = '../DAT_FILES/out_arrayfile.txt'

         call vec2array_file(A,n_row,n_col,outfile,'a')
         end
         include '../vec2array_file.f'
         !include '../fort_lib.f'
