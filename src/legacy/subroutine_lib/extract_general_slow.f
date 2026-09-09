C --------------------------------------------------------------
! This subroutine pre-computes the sine and the cosine templates 
! relevant for extraction. The idea is to save unnecessary comp-
! utation time. These templates are a function of only the 
! sampled x-arguments (eg., t's or Lsq's) and the angular 
! frequenciies at which the extraction is sought. It happens in 
! certain cases that one may need to call the extract_general 
! routine several times without the parameters of observation 
! changing -- eg., tomography using a FITS image cube or along 
! a pulsar's longitude bins. In these cases, the time that would 
! be saved by pre-computing these templates would be enormous. 
!  -- wasim, 08 Apr, 2011


C --------------------------------------------------------------



          !subroutine extract_general_setup(t,ryt,iyt,npts,fac,
!     -               nout,nu, p_ex, phi_ex,
!     -               order_taylor)
          
          !
          !      t-> input time series
          !   npts-> number of data samples
          !   nout-> number of output components (= npts x ofac)
          !    fac-> factor in the uncertainty relation between the 2 domains
          !          e.g., 1.0 if nu <--> t
          !                 pi if RM <--> lambda**2 etc.
          !--------------------------------------------------------------------

          subroutine extract_general_slow(t,ryt,iyt,npts,fac,nout,nu,
     -                                     p_ex,phi_ex) 
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          real*4      t(*), ryt(*), iyt(*) 
          real*4      nu(*), p_ex(*), phi_ex(*) 
          real*4      fac
          real*4      t1, t2, dt , d_nu 
          integer*4   npts, nout 

          real*4      t_span, nu_span 
          real*4      omega, atemp
          real*4      h_tmp, phi_tmp 

          real*4      c_template(npts), s_template(npts)
          real*4      rc_cor, ic_cor, rs_cor, is_cor
          real*4      ryw_tmp, iyw_tmp           
             

C COUNTERS:
          integer*4   i,kk,nzero,n_positive

C CONSTANTS:
          real*4      pi,twopi


          pi = 3.14159265358979 !acos(-1.0d0)
          twopi = 6.28318530717959 !2.0d0*pi

          ! Calculate the edge t: 
          !dt = (t(npts) - t(1))/dble(npts-1) 
          dt = t(2) - t(1) 
          ! BW_MHz = (dble(npts)/dble(npts-1))*(freq_MHz(npts)-freq_MHz(1))
          t1 = t(1) - 0.5d0*dt 
          t2 = t(npts) + 0.5d0*dt 

          t_span = t2 - t1 
          d_nu = fac/t_span   ! The factor fac is pi for the cases like 
                              ! RM and Lambda**2 kind of extraction
                              ! For nu vs t kind of extraction, fac = 1.0
          nu_span = dble(npts)*d_nu 

          h_tmp = d_nu*(real(npts-1)/real(nout-1))

          ! Ensure to sample the zero (ie., location of mean):
          if(mod(nout,2).eq.0)then
                  nzero = nout/2 + 1
          else
                  nzero = (nout + 1)/2
          endif
          n_positive = nzero - 1 ! Ensures equal number of
                                 ! components on either side
                                 ! of zero. Hence the output
                                 ! no. of comps on one side
                                 ! (including zero-comp) is
                                 ! always even.
                                 ! _-_-_-|-_-_-_

          do i = 1,n_positive
             atemp = dble(i)*h_tmp
             nu(nzero+i) = atemp*0.5
             nu(nzero-i) = -atemp*0.5
          end do
          nu(nzero) = 0.0

cc          write(*,*)'#### nout,nzero,n_positive :',nout,nzero,n_positive
cc          write(*,*)'#### nu : ',(nu(i),i=1,nout)


          !-------------------------------------------------
          do i = 1,nout ! number of nu's
             omega = 2.0*pi*nu(i)  ! assuming nu_trial to be angular frequency
             do kk = 1,npts ! number of t's
                phi_tmp = omega*t(kk)
                c_template(kk) = cos(phi_tmp)
                s_template(kk) = -sin(phi_tmp)
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
