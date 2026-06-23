chelp+
      ! Code to test the authenticity of rob_meanrms.f code.
      ! --wasim raja, 07 August, 2010
chelp-

      implicit none
      real*4       x(4096), xmean, xrms, 
     -             rob_mean, rob_rms, accum
      integer*4    i, npts, nchar
      character*72 infile

      accum = 0.0
      !infile = 'example_rob_meanrms.txt'
      infile = 'fort.31'
      open(11,file=infile,status='old',err=101)
      goto 102
101   write(*,*)"Error opening file: ",infile(1:nchar(infile))
      write(*,*)"Quitting now..."
      stop
102   continue
      i = 0
      do while(.true.)
         i = i + 1
         read(11,*,end=201)x(i)
         accum = accum + x(i)
      enddo
      
201   continue
      close(11)
      npts = i
      accum = accum/real(npts)
      write(*,*)"accum = ", accum

      call mean(x,npts,xmean)
      call rms(x,npts,xrms)

      call meanrms(x,rob_rms,rob_mean,npts)
      write(*,*)"------------------------------- "
      write(*,*)"       mean in data: ",xmean
      write(*,*)"        rms in data: ",xrms
      write(*,*)"        npt in data: ",npts
      write(*,*)" "
      write(*,*)"robust mean in data: ",rob_mean
      write(*,*)" robust rms in data: ",rob_rms
      write(*,*)"------------------------------- "


      end
      include '/usr/lib/subroutine_lib/robust_meanrms.f'
      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
