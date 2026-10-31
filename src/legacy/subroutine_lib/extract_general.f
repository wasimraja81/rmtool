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



          !subroutine extract_general(t,ryt,iyt,npts,fac,ofac, 
!     -               > nu,
!     -               p_ex,phi_ex,
!     -               rp_ex, rphi_ex, 
!     -               ip_ex, iphi_ex)
          
          !
          !      t-> input time series
          !    ryt-> real(y(t))
          !    iyt-> imag(y(t))
          !   npts-> number of data samples
          !    fac-> factor in the uncertainty relation between the 2 domains
          !          e.g., 1.0 if nu <--> t
          !                 pi if RM <--> lambda**2 etc.
          !   ofac-> oversampling factor in determining spectral components
          !     nt-> number of trial frequencies = npts*ofac
          !     nu-> frequency values at which spectral power is sought
          !   p_ex-> total power extracted from a complex signal
          ! phi_ex-> phase angle between the real and imaginary part of y(w)
          !  rp_ex-> power extracted from the real part of the signal
          !rphi_ex-> phase of the real part of the signal at a given frequency
          !  ip_ex-> power extracted from the imag part of the signal
          !iphi_ex-> phase of the imag part of the signal at a given frequency


          subroutine extract_general(t,ryt,iyt,npts,fac,ofac,nu,
     -               p_ex,phi_ex,
     -               rp_ex, rphi_ex, 
     -               ip_ex, iphi_ex)
          
          implicit none

          include '../INCLUDE/extract_rm.inc'

          real*4 t(*), ryt(*), iyt(*), nu(max_chan)
          real*4 fac
          integer*4 ofac
          integer*4 npts,nt

          real*4 t_span
          real*4 beg_nu, end_nu, d_nu, nu_span
          real*4 nu_trial(maxnt)
          real*4 omega


          real*4 rc_cor, ic_cor, rs_cor, is_cor
          real*4 p_ex(max_chan),   rp_ex(max_chan),   ip_ex(max_chan)
          real*4 phi_ex(max_chan), rphi_ex(max_chan), iphi_ex(max_chan)
          real*4 pmax
          real*4 p_tmp, ryw_tmp, iyw_tmp
             
          real*4 c_template(max_chan), s_template(max_chan)
          real*4 std_c_template, std_s_template
          real*4 std_ryt, std_iyt

C COUNTERS:
          integer*4 i,j,k,kk

C CONSTANTS:
          real*4 pi, vel_light


          pi = acos(-1.0)
          vel_light = 3.0e8
          nt = ofac*npts

C Relation between the 2 domains:
          t_span = t(npts) - t(1)
          d_nu = fac/t_span  ! The factor fac is pi for the cases like 
                             ! RM and Lambda**2 kind of extraction
                             ! For nu vs t kind of extraction, fac = 1.0
          nu_span = npts*d_nu

C Trial frequencies:
          beg_nu = -0.5*nu_span
          end_nu = 0.5*nu_span
          call linspace(beg_nu,end_nu,nt,nu_trial) ! oversampled nu-space
          call linspace(beg_nu,end_nu,npts,nu)     ! nu to be output

          call rms(ryt,npts,std_ryt)
          call rms(iyt,npts,std_iyt)

          k = 0
          do i = 1,npts ! number of nu's
             pmax = -1.0
             do j = 1,ofac
                k = k + 1
                !omega = 2*pi*nu_trial(k)
                omega = nu_trial(k) ! assuming nu_trial to be angular frequency
                do kk = 1,npts ! number of t's
                   c_template(kk) = cos(omega*t(kk))
                   s_template(kk) = sin(omega*t(kk))
                enddo
                call rms(c_template,npts,std_c_template)
                call rms(s_template,npts,std_s_template)

                call corr(ryt,c_template,rc_cor,npts)
                call corr(ryt,s_template,rs_cor,npts)

                call corr(iyt,c_template,ic_cor,npts)
                call corr(iyt,s_template,is_cor,npts)

                rc_cor = rc_cor*std_ryt*std_c_template
                rs_cor = rs_cor*std_ryt*std_s_template

                ic_cor = ic_cor*std_iyt*std_c_template
                is_cor = is_cor*std_iyt*std_s_template
          
         ! Combine coherently to construct y(w)
                ryw_tmp = rc_cor - is_cor  ! Real-part of y(w)
                iyw_tmp = rs_cor + ic_cor  ! Imag-part of y(w)
                p_tmp = sqrt(ryw_tmp**2 + iyw_tmp**2)

         ! Store the maxima of the power extracted from 
         ! over-sampled frequencies within a resolution bin
                if(p_tmp.ge.pmax)then
                        pmax = p_tmp
                        p_ex(i) = pmax  
                        phi_ex(i) = atan2(iyw_tmp,ryw_tmp) 

         ! phase & amp of the real-part of the signal y(t) at 
         ! angular frequency omega:
                        rphi_ex(i) = atan2(rs_cor,rc_cor)
                        rp_ex(i) = sqrt(rc_cor**2 + rs_cor**2)

         ! phase & amp of the imag-part of the signal y(t) at 
         ! angular frequency omega:
                        iphi_ex(i) = atan2(is_cor,ic_cor)
                        ip_ex(i) = sqrt(ic_cor**2 + is_cor**2)
                endif
             enddo
          enddo
C -------------------------------------------------          
          end
