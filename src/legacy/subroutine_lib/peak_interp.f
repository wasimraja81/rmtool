chelp+
      !-------------------------------------------------
      ! Routine to find the exact location of the peak 
      ! in a spectral bin. 
      ! Several schemes are allowed:
      ! Interp Type = 1 --> Quadratic Interpolation(in use) 
      ! Interp Type = 2 --> Fourier Interolation (amplitude used)
      !                     [This may be used for interpolating 
      !                      real-arrays by supplying the imaginary 
      !                      array with a zero-filled one]
      ! Interp Type = 3 --> Fourier Interolation (complex interp)
      ! Interp Type = 4 --> Sinc Interpolation 
      ! Interp Type = ? --> Parabolic Interpolation 
      ! ------------------------------------------------
chelp-


      subroutine peak_interp(InArrRe,InArrIm,narr,RM1, RM2, 
     -                       ofac,peak_loc, peak_val,OutPha, 
     -                       interp_type)


      implicit none


      real*4      InArrRe(*), InArrIm(*), peak_loc
      integer*4   narr, nout, ofac, interp_type
      real*4      tmpInArrRe(narr), tmpInArrIm(narr)
      real*4      InArr(narr)
      real*4      peak_val 
      real*4      OutArr(narr*ofac*2), dh, RM1, RM2, dh_in, dh_out 
      real*4      OutArrRe(narr*ofac*2), OutArrIm(narr*ofac*2)
      integer*4   i, imax
      
      real*4      OutPha


      do i = 1,narr
         tmpInArrRe(i) = InArrRe(i)
         tmpInArrIm(i) = InArrIm(i)

         InArr(i) = sqrt(tmpInArrRe(i)**2 + tmpInArrIm(i)**2)
      enddo
      
      ! TODO: Interpolation of phase for all other interpolation 
      !       schemes except for Complex Fourier Interpolation 
      !       to be done. 
      if(interp_type.eq.1)then ! Parabolic Interpolation
              dh = (RM2 - RM1)/real(narr-1)
              call quad_interp(InArr,narr,peak_loc,peak_val)
              peak_loc = RM1 + (peak_loc - 1.0)*dh
              write(*,*)"xxxxxxxxxxxxxxx Warning xxxxxxxxxxxxxxxx"
              write(*,*)"Phase at Interpolated RM not calculated!"
              write(*,*)"Incorporate the scheme NOW!!! "
              write(*,*)" "
      else if(interp_type.eq.2)then ! Fourier interpolation on Amp
              call fourier_interp_re(InArr, OutArr, narr, ofac, nout)
              call maxima_index(OutArr,nout,imax)

              dh_in = (RM2 - RM1)/real(narr-1)
              dh_out = dh_in*real(narr-1)/real(nout-1)

              peak_loc = RM1 + real(imax - 1)*dh_out
              peak_val = OutArr(imax)
              write(*,*)"xxxxxxxxxxxxxxx Warning xxxxxxxxxxxxxxxx"
              write(*,*)"Phase at Interpolated RM not calculated!"
              write(*,*)"Incorporate the scheme NOW!!! "
              write(*,*)" "
      else if(interp_type.eq.3)then ! Fourier interpolation Cmplx
              call fourier_interp(tmpInArrRe,tmpInArrIm,
     -                            OutArrRe,OutArrIm,narr, 
     -                            ofac, nout)
              do i = 1,nout
                 OutArr(i) = sqrt(OutArrRe(i)**2 + OutArrIm(i)**2)
              enddo
              call maxima_index(OutArr,nout,imax)

              dh_in = (RM2 - RM1)/real(narr-1)
              dh_out = dh_in*real(narr-1)/real(nout-1)

              peak_loc = RM1 + real(imax - 1)*dh_out
              peak_val = OutArr(imax)
              OutPha = atan2(OutArrIm(imax),OutArrRe(imax))
      else if(interp_type.eq.4)then ! Sinc Interpolation
              call sinc_interp(InArr,OutArr,narr,ofac,nout)
              call maxima_index(OutArr,nout,imax)
              dh = (RM2 - RM1)/real(nout-1)
              peak_loc = RM1 + real(imax - 1)*dh
              peak_val = OutArr(imax)
              write(*,*)"xxxxxxxxxxxxxxx Warning xxxxxxxxxxxxxxxx"
              write(*,*)"Phase at Interpolated RM not calculated!"
              write(*,*)"Incorporate the scheme NOW!!! "
              write(*,*)" "
      else    ! Parabolic Interpolation as default
              dh = (RM2 - RM1)/real(narr-1)
              call quad_interp(InArr,narr,peak_loc,peak_val)
              peak_loc = RM1 + (peak_loc - 1.0)*dh
              write(*,*)"xxxxxxxxxxxxxxx Warning xxxxxxxxxxxxxxxx"
              write(*,*)"Phase at Interpolated RM not calculated!"
              write(*,*)"Incorporate the scheme NOW!!! "
              write(*,*)" "
      endif

      return
      end

