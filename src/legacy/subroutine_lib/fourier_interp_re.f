chelp+
      !-------------------------------------------
      ! Routine for interpolating functions using  
      ! zero-padding in the Fourier Domain and
      ! inverse transforming back.
      ! Inherently, Complex functions are assumed; 
      ! in case the function is real, supply the 
      ! imaginary part of the array after filling 
      ! with zeros. 
      ! 
      ! The output over-sampled array although is 
      ! supposed to have npts_input x ofac number 
      ! of elements, we force the number of output 
      ! elements to that power of 2 which just 
      ! exceeds npts_input x ofac. 
      !
      ! Associated codes needed: 
      !     1) FFT_GENERAL package
      !         --> fft1d.f
      !         --> ifft1d.f
      !         --> fft_general_lin.f
      !     2) fort_lib.f
      !     3) nchar.f
      ! 
      !                    --wr, 24 Nov, 2010.
      ! 
      ! This code is modified to interpolate real 
      ! data.
      !-------------------------------------------
chelp-


      subroutine fourier_interp_re(InArrRe,OutArrRe,
     -                          npts,ofac,nout)

      implicit none


      integer*4  npts, ofac, prod, nout 
      real*4     InArrRe(*) 
      real*4     OutArrRe(*), OutArrIm(npts*ofac*2)
      real*4     tmpArrRe(npts), tmpArrIm(npts), norm 
      integer*4  i, j, n1, npad 



      nout = npts*ofac

      !-------------------------------
      ! Make nt a multiple of 2**n
      i = 0
      prod = 1
      do while(prod.lt.nout)
        i = i + 1
        prod = prod*2
      enddo
      nout = prod
      !-------------------------------

      do i = 1,nout
         OutArrRe(i) = 0.0
         OutArrIm(i) = 0.0
      enddo

      do i = 1,npts
         tmpArrRe(i) = InArrRe(i)
         tmpArrIm(i) = 0.0 
      enddo

      call fft1d(tmpArrRe,tmpArrIm,npts)

      if(mod(npts,2).eq.0)then
              n1 = npts/2  
      else
              n1 = (npts-1)/2 
      endif
      npad = nout - npts

      !do i = 1,n1
      do i = 1,n1
         OutArrRe(i) = tmpArrRe(i)
         OutArrIm(i) = tmpArrIm(i)
      enddo
      j = n1
      do i = n1+npad+1,nout
         j = j + 1
         OutArrRe(i) = tmpArrRe(j)
         OutArrIm(i) = tmpArrIm(j)
      enddo

      call ifft1d(OutArrRe,OutArrIm,nout)
      norm = real(nout)/real(npts)
      do i = 1,nout
         OutArrRe(i) = norm*OutArrRe(i)
         !OutArrIm(i) = real(nout/npts)*OutArrIm(i)
      enddo

      return

      end

