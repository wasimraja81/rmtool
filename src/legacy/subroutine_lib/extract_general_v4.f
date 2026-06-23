C --------------------------------------------------------------
! This subroutine extracts the power in components of frequency 
! present in a modulated signal using a simple matched filtering 
! technique.
! It is assumed that the signal is of the form:
!             y(t) = A*exp[i*omega*t]
! This code thus extracts the angular frequency "omega" rather 
! than the rfequency "nu". 

!  -- wasim, 11 Aug, 2009


C --------------------------------------------------------------



          !subroutine extract_general(t,ryt,iyt,npts,fac,
!     -               nout,nu, p_ex, phi_ex,
!     -               order_taylor)
          
          !
          !      t-> input time series
          !    ryt-> real(y(t))
          !    iyt-> imag(y(t))
          !   npts-> number of data samples
          !   nout-> number of output components (= npts x ofac)
          !    fac-> factor in the uncertainty relation between the 2 domains
          !          e.g., 1.0 if nu <--> t
          !                 pi if RM <--> lambda**2 etc.
          !   ofac-> oversampling factor in determining spectral components
          ! omega1-> beginning frequency values at which spectral power is sought
          !   p_ex-> total power extracted from a complex signal
          ! phi_ex-> phase angle between the real and imaginary part of y(w)
          !  rp_ex-> power extracted from the real part of the signal
          !rphi_ex-> phase of the real part of the signal at a given frequency
          !  ip_ex-> power extracted from the imag part of the signal
          !iphi_ex-> phase of the imag part of the signal at a given frequency

          ! Last Modification: 
          !         --> The extraction is done for a user-defined range
          !             of angular frequencies instead of the full range 
          !             of unaliased angular frequencies.
          !         --> Also the number of components "nout" is to be supplied 
          !             to this subroutine.
          ! NOTE: The need to extract the power at only a subset of the 
          !       unaliased range arises due to the fact that the
          !       full-extraction is a rather time-consuming process.
          !       Moreover we may already know WHERE to look for the 
          !       power in the spectral axis. 
          !
          ! Last Modification:
          !         --> Bandwidth depolarisation correction done.
          !             --wr, 04 Nov, 2010
          !         --> BW-depol Correction is an option now, 
          !             governed by the value of the input variable
          !             "order_taylor". If negative, no correction.
          !             --wr, 23 Nov, 2010
          ! Last Modification: 
          !         --> Bandwidth depolarization correction exluded 
          !             from this routine. Th eright thing to do is 
          !             to incorporate BW-depol correction in the 
          !             cleaning routine. 
          !             --wr, 08 Apr, 2011
          !
          !         --> We make this code computationally efficient 
          !             where one needs to make multiple calls to 
          !             this subroutine with same range of lambda's 
          !             and RM's, by passing arrays of pre-computed 
          !             sine and cosine templates.
          !             For this purpose, a set-up subroutine has 
          !             been written for precomputing the templates. 
          !             The setup subroutine may be called only once,  
          !             while the extrcat_general subroutine may be 
          !             called as many times as is required in a given 
          !             run of the main program. 
          !             --wr, 08 Apr, 2011
          !
          ! Last Modification: 
          !         --> Provision for removing the "mean" from the 
          !             Q and U data incorporated. 
          !             --wr, Feb, 2012 
          !--------------------------------------------------------------------

!          subroutine extract_general(t,ryt,iyt,npts,fac,
!     -               nout,nu, p_ex, phi_ex,
!     -               cos_arr, sine_arr)
          subroutine extract_general(ryt_in,iyt_in,npts,
     -               nout, p_ex, phi_ex,
     -               cos_arr, sin_arr,maxout,maxpts,mean_rem)
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          real*4      ryt_in(*), iyt_in(*) 
          real*4      p_ex(*), phi_ex(*) 
          integer*4   npts, nout, maxout, maxpts
          real*4      c_template(npts), s_template(npts)
          real*4      ryt(npts), iyt(npts) 

          real*4      cos_arr(maxout,maxpts), 
     -                sin_arr(maxout,maxpts)

          real*4      rc_cor, ic_cor, rs_cor, is_cor
          real*4      ryw_tmp, iyw_tmp 
          integer*4   mean_rem 
             

C COUNTERS:
          integer*4   i,kk



          ! Remove mean from Q and U if sought: 
          if (mean_rem.gt.0)then
                  call mean(ryt_in,npts,ryw_tmp)
                  call mean(iyt_in,npts,iyw_tmp)
                  do i = 1,npts
                     ryt(i) = ryt_in(i) - ryw_tmp
                     iyt(i) = iyt_in(i) - iyw_tmp
                  enddo
          else 
                  do i = 1,npts
                     ryt(i) = ryt_in(i) 
                     iyt(i) = iyt_in(i) 
                  enddo
          endif
          !-------------------------------------------------
          do i = 1,nout ! number of nu's
             do kk = 1,npts ! number of t's
                c_template(kk) = cos_arr(i,kk)
                s_template(kk) = sin_arr(i,kk)
             enddo
  
             call dotproduct(ryt,c_template,rc_cor,npts)
             call dotproduct(ryt,s_template,rs_cor,npts)
             call dotproduct(iyt,c_template,ic_cor,npts)
             call dotproduct(iyt,s_template,is_cor,npts)
             rc_cor = rc_cor/dble(npts)
             rs_cor = rs_cor/dble(npts)
             ic_cor = ic_cor/dble(npts)
             is_cor = is_cor/dble(npts)

         ! Combine coherently to construct y(w)
             ryw_tmp = rc_cor - is_cor  ! Real-part of y(w)
             iyw_tmp = rs_cor + ic_cor  ! Imag-part of y(w)
             p_ex(i) = sqrt(ryw_tmp*ryw_tmp + iyw_tmp*iyw_tmp)
             phi_ex(i) = atan2(iyw_tmp,ryw_tmp)
          enddo
C -------------------------------------------------          
          end
