chelp+
      !----------------------------------------------
      ! This routine uses Quadratic Interpolation 
      ! for locating the peak of an amplitude spectrum 
      ! within a nominally oversampled spectral bin. 
      ! 
      ! The sampled peak and its 2 nearest neighbour 
      ! on either side of it is used for the 
      ! interpolation. 
      ! If the sampled peak falls on the edge bins, 
      ! then the sampled peak is regarded as the 
      ! true peak. 
      !
      ! The interpolation is performed after the 
      ! amplitudes are compressed to their log10. 
      ! This is to take care of extreme values of 
      ! the ordinates. 
      ! 
      !                -- wr, 30 Oct, 2010
      !----------------------------------------------
chelp-

      subroutine peak_find(indata,npts,peak)

      !----------------------------------------------------------
      !     ya = input array whose peak is to be determined
      !   npts = no. of pts in input array
      ! > peak = location of the bin containing the peak. This is 
      !          a fractional bin number. The calling program 
      !          should coonovert it to the abscissa value if so 
      !          needed. 
      !----------------------------------------------------------

      implicit none

      real*8    indata(*), peak, alpha, beta, gama 
      integer*4 npts, i, peak_now 
      real*8    ya(npts), sample_max 

      ! Convert to logarithmic scale to take care 
      ! of dynamic range related problems.
      do i = 1,npts
         ya(i) = log10(indata(i))
         !ya(i) = indata(i)
      enddo

      ! Locate 3 successive points around the highest 
      ! peak sampled
      ! Initialisation: 
      sample_max = ya(1)
      peak_now = 1

      do i = 1,npts
         if(ya(i).gt.sample_max)then
                 sample_max = ya(i)
                 peak_now = i
         endif
      enddo

      if(peak_now.eq.1.or.peak_now.eq.npts)then
              peak = real(peak_now)
      else
              alpha = ya(peak_now-1)
              beta = ya(peak_now)
              gama = ya(peak_now+1)

              peak = peak_now + 0.5*((alpha-gama)/(alpha+gama-2.0*beta))
      endif

      return
      end
