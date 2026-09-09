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
          !         --> Calculation of redundant trigonometric functions
          !             which is time-consuming, is dome away with. 
          !                                  --ymaan, 16 Sep, 2010
          !--------------------------------------------------------------------


          subroutine extract_general(t,ryt,iyt,npts,fac,ofac,
     -               omega1, nout,nu,
     -               p_ex, phi_ex,
     -               rp_ex, rphi_ex, 
     -               ip_ex, iphi_ex,
     -               fullrange)
          
          implicit none

          include '../INCLUDE/extract_rm.inc'

          real*4 t(*), ryt(*), iyt(*), nu(*)
          real*4 fac
          integer*4 ofac
          integer*4 npts,nt, nout

          real*4 t_span
          real*4 beg_nu, end_nu, d_nu, nu_span
          real*4 nu_trial(maxnt)
          real*4 omega, omega1
          logical fullrange


          real*4 rc_cor, ic_cor, rs_cor, is_cor
          real*4 p_ex(*),   rp_ex(*),   ip_ex(*)
          real*4 phi_ex(*), rphi_ex(*), iphi_ex(*)
          real*4 pmax
          real*4 p_tmp,ryw_tmp,iyw_tmp,h_tmp,phi_tmp,
     _           atemp,phi_tmp1
             
          real*4 c_template(maxchan), s_template(maxchan)

C COUNTERS:
          integer*4 i,j,k,kk

C CONSTANTS:
          real*4 pi,twopi


          pi = acos(-1.0)
          twopi = 2.0*pi

C Relation between the 2 domains:
          t_span = t(npts) - t(1) 
     _             + (t(2)-t(1))/2.0 + (t(npts)-t(npts-1))/2.0
          d_nu = fac/t_span  ! The factor fac is pi for the cases like 
                             ! RM and Lambda**2 kind of extraction
                             ! For nu vs t kind of extraction, fac = 1.0
C Trial frequencies:
          if(fullrange)then
                  nout = npts
                  nu_span = real(nout-1)*d_nu
                  beg_nu = -0.5*nu_span
          else
                  nu_span = real(nout-1)*d_nu
                  beg_nu = omega1
          endif
          !nt = ofac*npts
          nt = ofac*nout

          end_nu = beg_nu + nu_span

          !call linspace(beg_nu,end_nu,nt,nu_trial) ! oversampled nu-space
          !!!call linspace(beg_nu,end_nu,npts,nu)     ! nu to be output
          !call linspace(beg_nu,end_nu,nout,nu)     ! nu to be output
          h_tmp = (end_nu - beg_nu)/real(nt-1)
          do i = 1,nt   
             nu_trial(i) = beg_nu + real(i-1)*h_tmp
          end do
          h_tmp = (end_nu - beg_nu)/real(nout-1)
          do i = 1,nt   
             nu(i) = beg_nu + real(i-1)*h_tmp
          end do



          k = 0
          do i = 1,nout ! number of nu's
             pmax = -1.0
             do j = 1,ofac
                k = k + 1
                !omega = 2*pi*nu_trial(k)
                omega = nu_trial(k) ! assuming nu_trial to be angular frequency
                do kk = 1,npts ! number of t's
                   phi_tmp1 = omega*t(kk)
                   phi_tmp = phi_tmp1 - twopi*int(phi_tmp1/twopi)
                   c_template(kk) = cos(phi_tmp)
                   ! following finds out the sin(phi_tmp) using the
                   ! above value of cos(phi_tmp), ymaan, 16Sep, 2010
                   atemp = c_template(kk)
                   s_template(kk)= sqrt(1.0 - atemp*atemp)
                   if(phi_tmp.lt.0) phi_tmp = phi_tmp + twopi
                   if(phi_tmp.ge.pi)then
                     s_template(kk)= -s_template(kk)
                   end if
                enddo

                call dotproduct(ryt,c_template,rc_cor,npts)
                call dotproduct(ryt,s_template,rs_cor,npts)
                call dotproduct(iyt,c_template,ic_cor,npts)
                call dotproduct(iyt,s_template,is_cor,npts)
                rc_cor = rc_cor/real(npts)
                rs_cor = rs_cor/real(npts)
                ic_cor = ic_cor/real(npts)
                is_cor = is_cor/real(npts)

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
