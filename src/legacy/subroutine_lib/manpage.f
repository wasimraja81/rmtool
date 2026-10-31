chelp+
! This code is meant to extract and display the "HELP" section of 
! the programs containing a help section at the beginning of the 
! file. 
!                     --wasim raja, rri, 08 Nov, 2009
chelp-


      subroutine manpage(prg_name)

      implicit none
      character*120 prg_name
      character*112 readline
      logical fread
      integer*4 i, k, ibeg, iend
      integer*4 nchar_sr

      write(*,*)" "
      write(*,*)"Displaying HELP for: ",prg_name(1:nchar_sr(prg_name))
      write(*,*)" "

      open(13,file=prg_name,status='old',err=101)
      goto 102

101   write(*,*)"Error opening file: ",prg_name(1:nchar_sr(prg_name))
      write(*,*)"No 'help' for: ",prg_name(1:nchar_sr(prg_name))
      stop

102   continue
      fread = .true.
      ibeg = 0
      do while(fread)
         read(13,'(a)',end=103)readline
         ibeg = ibeg + 1
         if(index(readline,'chelp+').gt.0)then 
                 fread = .false.
         endif
      enddo
      goto 104

103   write(*,*)"Start of Help section undefined in: ",
     -                                 prg_name(1:nchar_sr(prg_name))
      write(*,*)"No 'help' for: ",prg_name(1:nchar_sr(prg_name))
      stop

104   continue
      fread = .true.
      iend = 0
      do while(fread)
         read(13,'(a)',end=105)readline
         iend = iend + 1
         if(index(readline,'chelp-').gt.0)then 
                 fread = .false.
         endif
      enddo
      goto 106

105   write(*,*)"End of Help section undefined in: ",
     -                                 prg_name(1:nchar_sr(prg_name))
      write(*,*)"No 'help' for: ",prg_name(1:nchar_sr(prg_name))
      stop


      ! Now start printing the help page:
106   continue
      close(13)

      open(13,file=prg_name,status='old',err=201)
      goto 202

201   write(*,*)"Error opening file: ",prg_name(1:nchar_sr(prg_name))
      stop

202   continue

      do i = 1,ibeg
         read(13,'(a)')readline
      enddo

      k = 0
      do i = ibeg,iend-1
         k = k + 1
         read(13,'(a)')readline
         write(*,*)readline(1:nchar_sr(readline))
         if(mod(k,15).eq.0)then
                 write(*,*)" "
501              continue
                 write(*,*)'press Q to quit, C to continue with help...'
                 read(*,*)readline
                 if(readline(1:1).eq.'Q'.or.
     -              readline(1:1).eq.'q')then
                     write(*,*)"Done with help..."
                     stop
                 else if(readline(1:1).ne.'C'.and.
     -              readline(1:1).ne.'c')then
                    goto 501
                 endif
         endif
      enddo
      write(*,*)" "
      close(13)

      end


      Integer function NCHAR_SR(string)
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

      NCHAR_SR = ipos
      return
      end
