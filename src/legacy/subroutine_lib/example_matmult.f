chelp+
      !-----------------------------------------
      ! This is a driver routine to demonstrate 
      ! usage of the matmult_2d subroutine.
      !                      --wr, 31 Jan, 2013
      !-----------------------------------------
chelp-


      implicit none 

      integer*4  maxdim_x, maxdim_y 
      parameter  (maxdim_x=4,maxdim_y=4)
      real*4     A(maxdim_x,maxdim_y), B(maxdim_x,maxdim_y), 
     -           P(maxdim_x,maxdim_y) 
      integer*4  i, j, nx, ny 

      nx = 4 
      ny = 4 
      do i = 1,nx
         do j = 1,ny
            !if (i .eq. j)then
            !        A(i,j) = 1.0  
            !        B(i,j) = 2.0 
            !else
            !        A(i,j) = 0.0 
            !        B(i,j) = 0.0 
            !endif
            A(i,j) = real(i + j)
            B(i,j) = real(i - j)
         enddo
      enddo
      call matmult_2d(A,nx,ny,maxdim_x,maxdim_y,
     -                B,nx,ny,maxdim_x,maxdim_y,
     -                P)

      do i = 1,nx
         write(*,*)"|",(A(i,j),j=1,ny),"|",
     -             " |",(B(i,j),j=1,ny),"|"
      enddo
      write(*,*)" "
      do i = 1,nx
         write(*,*)"|",(P(i,j),j=1,ny),"|"
      enddo

      end

      include '/usr/lib/subroutine_lib/fort_lib.f'
