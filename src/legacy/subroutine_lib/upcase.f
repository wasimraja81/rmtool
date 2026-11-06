      subroutine upcase(string)

c   to convert strings to Upper case characters

      implicit  none

      character*(*) string
      integer*4     nchar,
     -              ichar,
     -              istring,
     -              i
      character     char

      
      do i=1,nchar(string)
        istring = ichar(string(i:i))
        if(istring.gt.96.and.istring.lt.123) then   ! lower case
          string(i:i) = char(istring-32)
        end if
      end do
     
      return
      end
