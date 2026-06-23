chelp+ 
      !---------------------------------------
      ! This code displays an error message 
      ! about file-open error. Hence if there 
      ! is trouble opening a file from a 
      ! program, this code may be used. 
      ! 
      !                  --wr, 03 Aug, 2013 
      !---------------------------------------
chelp- 

      subroutine file_open_error(filename)

      implicit none 

      character filename*(*)
      integer*4  nchar 
      write(*,*)"Error opening file: ",filename(1:nchar(filename))

      end
