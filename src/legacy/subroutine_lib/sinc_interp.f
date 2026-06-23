chelp+
      !-------------------------------------------
      ! Routine for sinc interpolating functions   
      ! 
      ! Associated codes needed: 
      !     2) fort_lib.f
      !     3) nchar.f
      ! 
      !                    --wr, 24 Nov, 2010.
      !-------------------------------------------
chelp-

      subroutine sinc_interp(InArr,OutArr,npts,ofac,nout)

      implicit none


      integer*4  npts, ofac, nout 
      real*4     InArr(*), OutArr(*) 
      real*4     accum, arg, val_at, h, samp, dh, pi 
      integer*4  i, j 



      pi = acos(-1.0)
      nout = npts*ofac


      dh = real(npts-1)/real(nout-1)
      do i = 1,nout
         accum = 0.0
         !val_at = 1.0 + real(i-1)*dh
         val_at = real(i-1)*dh
         do j = 1,npts
            samp = real(j) - 1.0
            arg = pi*(val_at - samp)
            if(arg.eq.0.0)then
                    accum = accum + InArr(j)
            else
                    h = sin(arg)/arg
                    accum = accum + InArr(j)*h
            endif
         enddo
         OutArr(i) = accum
      enddo


      return

      end

