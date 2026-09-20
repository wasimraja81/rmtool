      Integer function NCHAR(string)
C
C  Routine to count the number of characters in the
C  input string. Looks for the last occurrence of 
C  non-(null, blank or tab character)
C
C
C
      Implicit none
      integer*4 i,ipos
      character*(*)  string
      character blank,tab,null,c

C      data blank,tab,null/' ',9,0/
      blank=' '
        tab=char(9)
        null=char(0)

      ipos = 0
      i      = len(string)
      do while (i.gt.0.and.ipos.eq.0)
         c = string(i:i)
         if (c.ne.blank.and.c.ne.tab.and.c.ne.null) ipos = i
         i = i - 1
      end do

      NCHAR = ipos
      return
      end
