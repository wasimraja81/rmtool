chelp+ 
! This subroutine is meant for inverse  
! Fourier Transforming 1D data.
! It uses a the "fft_general_lin.f" routine.
! 
! Requires:
!          1) "fft_general_lin.f"
!      
! NOTE: This subroutine aligns the data 
!       such that the 0-th component of 
!       of the spectrum is placed at the
!       first pixel of the Fourier-Trans-
!       formed array. To shift the 0-th
!       component to the centre, use:
!       "fftshift.f" subroutine.
!       This shift operation conforms to
!       the "fftshift" operation in Octave.
!       
!       This routine should be used for 
!       only BACKWARD transformations. In 
!       other words, the sign of the 
!       argument of the exponential is 
!       taken as POSITIVE.
! IMPORTANT: The input arrays are modified 
!            to contain the transformed 
!            output. 
!            So be careful to preserve the 
!            input arrays if you NEED them!
!      
!       --wasim raja, 02 Jan, 2010
!-----------------------------------------------------------
chelp-

       subroutine ifft1d(rfx,ifx,npts)


       implicit none

       real*4    rfx(*)               ! value of function(real) 
       real*4    ifx(*)               ! value of function(imag) 
       integer*4 ierr, npts

       !integer*4 nshift


C --------------------------------------------------------------

       !nshift = int(npts/2)
       !call shift(rfx,npts,nshift)
       !call shift(ifx,npts,nshift)

       CALL FFT(rfx,ifx,1,npts,1,1,ierr)

       !call shift(rfx,npts,nshift)
       !call shift(ifx,npts,nshift)

       return

       end

