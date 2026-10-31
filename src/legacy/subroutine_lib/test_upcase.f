chelp+
      ! Code to demonstarte the usage of the subroutine 
      ! upcase.f
      !                               wr, 23 Oct 2012
chelp-


      implicit none 
      character*128  fname 
      integer*4       i, ichar, iachar 
      integer*4       nchar 

      character*1  alpha(26)


      write(*,*)"What is your first name: "
      read(*,'(a)')fname 
      do i = 1,nchar(fname)
         write(*,*)fname(i:i),' -> ',ichar(fname(i:i))
      enddo

      call upcase(fname)
      write(*,*)"Your first name is: ",fname 

      end

      include 'upcase.f'
      include 'nchar.f'
