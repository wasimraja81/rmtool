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

      subroutine quad_interp(indata,npts,peak,peak_val)

      !----------------------------------------------------------
      !     ya = input array whose peak is to be determined
      !   npts = no. of pts in input array
      ! > peak = location of the bin containing the peak. This is 
      !          a fractional bin number. The calling program 
      !          should convert it to the abscissa value if so 
      !          needed. 
      !----------------------------------------------------------
      ! Theory: 
      ! We wish to find the peak "p" for the curve: 
      !            y(x) = a(x - p)^2 + b
      ! interpolated at 3 points -- alpha, beta & gamma 
      ! about the sampled maxima and whose abscissa is 
      ! chosen arbitrarily to be : -1, 0 and 1
      !
      ! Therefore,
      !            alpha = a(-1 - p)^2 + b
      !                  = a.p^2 + 2.a.p + a + b
      !             beta = a.p^2 + b, and
      !            gamma = a(1 - p)^2 + b
      !                  = a.p^2 -2.a.p + a + b
      ! Solving for the peak location "p" and the value 
      ! at the peak, "b", we get: 
      !             
      !                p = 0.5(alpha-gamma)/(alpha - 2.beta + gamma)
      !                b = 

      implicit none

      real*4    indata(*), peak, alpha, beta, gama 
      integer*4 npts, i, peak_now 
      real*4    ya(npts), sample_max, peak_val 

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
      peak_val = sample_max

      do i = 1,npts
         if(ya(i).gt.sample_max)then
                 sample_max = ya(i)
                 peak_now = i
                 peak_val = sample_max
         endif
      enddo

      if(peak_now.eq.1.or.peak_now.eq.npts)then
              peak = real(peak_now)
      else
              alpha = ya(peak_now-1)
              beta = ya(peak_now)
              gama = ya(peak_now+1)

              peak = peak_now + 0.5*((alpha-gama)/(alpha+gama-2.0*beta))
              peak_val = beta - 0.25*(alpha-gama)*peak
      endif
      ! convert peak from log to linear:
      peak_val = 10.0**(peak_val)

      return
      end
