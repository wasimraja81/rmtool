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
          !         --> BW-depol correction is done on individual 
          !             stokes (Q or U) since we invert Q and U 
          !             one at a time in this scheme. This is to 
          !             facilitate cleaning of Q and U separately.
          !         --> Only one-half of the RM-spectra is used 
          !             for computation of depol-correction as well 
          !             as for computation the fourier components.
          !             (Hermitian Symmetry is invoked to produce 
          !             the other half of the RM-spectra.)
          !             
          !             --wr, 15 Feb, 2011
          !--------------------------------------------------------------------

          subroutine extract_general_real(t,ryt,iyt,npts,fac,
     -               nout,nu, p_ex, phi_ex,order_taylor,
     -               stokes)
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          real*4      t(*), ryt(*), iyt(*), nu(*) 
          real*4      p_ex(*), phi_ex(*) 
          character*1 stokes 
          real*4      fac
          real*4      f1, f2, Lsq1, Lsq2, dfreq
          integer*4   npts, nout, nzero, n_positive 
          real*4      ryt_corr(npts) 
          real*4      freq_MHz(npts) 
          real*4      imaccum, reaccum 

          real*4      t_span
          real*4      d_nu, nu_span !,beg_nu, end_nu 
          real*4      omega  
          real*4      nu_positive(nout)

          real*4      weight_amp(nout,npts), 
     -                weight_pha(nout,npts)
          real*4      weight, ph_weight
          integer*4   order_taylor

          real*4      rc_cor, rs_cor 
          real*4      h_tmp,phi_tmp 
             
          real*4      c_template(npts), s_template(npts)

C COUNTERS:
          integer*4   i,j,kk

C CONSTANTS:
          real*4      pi,twopi


          pi = acos(-1.0d0)
          twopi = 2.0d0*pi
          !order_taylor = 2

         ! Generate the temporal frequencies from L_sq data
         j = npts + 1
         do kk = 1,npts
            j = j - 1
            freq_MHz(j) = 300.0d0/sqrt(t(kk))
         enddo
         ! Calculate the edge L_sq: 
         dfreq = (freq_MHz(npts) - freq_MHz(1))/dble(npts-1)
         !BW_MHz = (dble(npts)/dble(npts-1))*(freq_MHz(npts)-freq_MHz(1))
         f1 = freq_MHz(1) - 0.5d0*dfreq
         f2 = freq_MHz(npts) + 0.5d0*dfreq
         Lsq2 = (300.0d0/f1)**2
         Lsq1 = (300.0d0/f2)**2
         !TEST: 
         !write(*,*)"BW_MHz: ",f2 - f1
C Relation between the 2 domains:
          !t_span = t(npts) - t(1) 
!!     _             + (t(2)-t(1))/2.0 + (t(npts)-t(npts-1))/2.0
          t_span = Lsq2 - Lsq1
          d_nu = fac/t_span  ! The factor fac is pi for the cases like 
                             ! RM and Lambda**2 kind of extraction
                             ! For nu vs t kind of extraction, fac = 1.0
C Trial frequencies:
          nu_span = dble(npts-1)*d_nu
          h_tmp = d_nu*real((npts-1)/(nout-1))

          ! Ensure to sample the zero (ie., location of mean):
          if(mod(nout,2).eq.0)then
                  nzero = nout/2 + 1
                  !beg_nu = -(nzero - 1)*h_tmp
                  !end_nu = (nout - nzero)*h_tmp
          else
                  nzero = (nout + 1)/2
                  !beg_nu = -(nzero - 1)*h_tmp
                  !end_nu = (nzero - 1)*h_tmp
          endif
          n_positive = nzero - 1 ! Ensures equal number of 
                                 ! components on either side 
                                 ! of zero. Hence the output 
                                 ! no. of comps on one side 
                                 ! (including zero-comp) is 
                                 ! always even. 
                                 ! _-_-_-|-_-_-_

          do i = 1,n_positive
             nu_positive(i) = dble(i)*h_tmp
             nu(i+nzero) = nu_positive(i)
             nu(i) = -nu_positive(i)
          end do
          nu(nzero) = 0.0

          !-------------------------------------------------
          ! BW-depolarisation correction related:
          ! (Obtain the complex correction factors for 
          ! bandwidth-depolarisation effect)
          if(order_taylor.ge.0)then
            call bw_depol_correct(freq_MHz, npts, nu_positive,
     -                          n_positive, order_taylor, 
     -                          weight_amp, weight_pha)
          !
          !-------------------------------------------------

            do i = 1,n_positive ! number of nu's
               ! Now correct the ryt and iyt arrays for bw-depol:
               j = npts + 1 
               do kk = 1,npts
                  ! The weights are in decreasing order of L_sq. 
                  ! So we need to arrange them in increasing order 
                  ! of L_sq here.
                  j = j - 1
                  weight = weight_amp(i,j)
                  ph_weight = weight_pha(i,j)
  
                  imaccum = atan2(iyt(kk),ryt(kk)) ! phase
                  reaccum = sqrt(ryt(kk)*ryt(kk)+iyt(kk)*iyt(kk))
  
                  ! Do the correction depending on the linear Stokes
                  ! type: 
                  if(stokes.eq.'q'.or.stokes.eq.'Q')then
                          ryt_corr(kk) = reaccum*
     -                              cos(imaccum - ph_weight)/weight
                  else if(stokes.eq.'u'.or.stokes.eq.'U')then
!                          iyt_corr(kk) = reaccum*
!     -                              sin(imaccum - ph_weight)/weight
                          ryt_corr(kk) = reaccum*
     -                              sin(imaccum - ph_weight)/weight
                  else
                          write(*,*)"Error: Wrong Stokes chosen!!"
                          write(*,*)" Stokes should be either Q or U"
                          write(*,*)" Quitting Now..."
                          stop
                  endif
               enddo
               omega = nu_positive(i)  ! assuming nu_trial to 
                                       ! be angular frequency
               do kk = 1,npts ! number of t's
                  phi_tmp = omega*t(kk)
                  c_template(kk) = cos(phi_tmp)
                  s_template(kk) = -sin(phi_tmp)
               enddo
  
               call dotproduct(ryt_corr,c_template,rc_cor,npts)
               call dotproduct(ryt_corr,s_template,rs_cor,npts)
               !call dotproduct(iyt_corr,c_template,ic_cor,npts)
               !call dotproduct(iyt_corr,s_template,is_cor,npts)
               rc_cor = rc_cor/dble(npts)
               rs_cor = rs_cor/dble(npts)
               !ic_cor = ic_cor/dble(npts)
               !is_cor = is_cor/dble(npts)

               p_ex(i+nzero) = sqrt(rc_cor**2 + rs_cor**2)
               p_ex(i) = p_ex(i+nzero)
               phi_ex(i+nzero) = atan2(rs_cor,rc_cor)
               phi_ex(i) = -phi_ex(i+nzero)
            enddo
          else   ! No BW-depol correction sought
            if(stokes.eq.'q'.or.stokes.eq.'Q')then
               do kk = 1,npts ! number of t's
                  ryt_corr(kk) = ryt(kk)
               enddo
            else if(stokes.eq.'u'.or.stokes.eq.'U')then
               do kk = 1,npts ! number of t's
                  ryt_corr(kk) = iyt(kk)
               enddo
            else
               write(*,*)"Error: Wrong Stokes chosen!!"
               write(*,*)" Stokes should be either Q or U"
               write(*,*)" Quitting Now..."
               stop
            endif
            do i = 1,n_positive        ! number of nu's
               omega = nu_positive(i)  ! assuming nu_trial to 
                                       ! be angular frequency
               do kk = 1,npts ! number of t's
                  phi_tmp = omega*t(kk)
                  c_template(kk) = cos(phi_tmp)
                  s_template(kk) = -sin(phi_tmp)
               enddo
  
               call dotproduct(ryt_corr,c_template,rc_cor,npts)
               call dotproduct(ryt_corr,s_template,rs_cor,npts)
               !call dotproduct(iyt,c_template,ic_cor,npts)
               !call dotproduct(iyt,s_template,is_cor,npts)
               rc_cor = rc_cor/dble(npts)
               rs_cor = rs_cor/dble(npts)
               !ic_cor = ic_cor/dble(npts)
               !is_cor = is_cor/dble(npts)

               p_ex(i+nzero) = sqrt(rc_cor**2 + rs_cor**2)
               p_ex(i) = p_ex(i+nzero)
               phi_ex(i+nzero) = atan2(rs_cor,rc_cor)
               phi_ex(i) = -phi_ex(i+nzero)
            enddo
          endif
          ! Now compute the zero-component in the spectrum:
          p_ex(nzero) = 0.0d0
          if(stokes.eq.'q'.or.stokes.eq.'Q')then
             do kk = 1,npts ! number of t's
                p_ex(nzero) = p_ex(nzero) + ryt(kk)
             enddo
          else if(stokes.eq.'u'.or.stokes.eq.'U')then
             do kk = 1,npts ! number of t's
                p_ex(nzero) = p_ex(nzero) + iyt(kk)
             enddo
          endif

          p_ex(nzero) = p_ex(nzero)/dble(npts)
          phi_ex(nzero) = 0.0
          nout = 2*n_positive + 1
C -------------------------------------------------          
          end
