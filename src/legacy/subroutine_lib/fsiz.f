************************************************************************
**  This subroutine take the file-name (with full path) as the input and
**  uses the command 'ls -ltr' to find out its size(in Bytes).
************************************************************************

        subroutine fsiz(file_name,siz)

           
        character*(*) file_name
        integer*8     siz, file_size, i, ipos  

        integer*4     itemp,n1,date
        character     permissions*10, author*15, group*15, month*10,
     -                time*10, temp_file*120
        character*228 temp_string
        character     blank, tab, nul, c


        siz = -1
        blank = ' '
        tab   = char(9)
        nul   = char(0)
        temp_file = file_name

*** **********************************************************
*** following block finds out the number of characters in the
*** file-name
        ipos = 0
        i    = len(temp_file)
        do while (i.gt.0.and.ipos.eq.0)
           c = temp_file(i:i)
           if (c.ne.blank.and.c.ne.tab.and.c.ne.nul) ipos = i
           i = i - 1
        end do
        n1 = ipos
*** **********************************************************

c        n1 = nchar(temp_file) 
        temp_string = 'ls -ltr '//temp_file
        write(temp_string(n1+9:),"(a)")' >temp_file_list.txt'
        call system(temp_string) 

        open(unit=11,file='temp_file_list.txt',status='unknown')
        read(11,"(a)",err=9111,end=9111)temp_string
        close(unit=11)
        call system('rm -f temp_file_list.txt')

9111    continue 
        read(temp_string,*,err=9112,end=9112)permissions,itemp,author,
     -                       group,file_size,month,date,time,temp_file
9112    continue

        if(file_size.gt.0) then 
                siz = file_size
        else
                siz = 0
        endif


        return
        end

***     ===========================
ccc        include 'codes/nchar_lin.f'
***     ===========================
