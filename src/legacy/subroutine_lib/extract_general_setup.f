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

          subroutine extract_general_setup(t,npts,fac,nout,nu, 
     -               cos_arr, sin_arr,maxout,maxpts)
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          real*4      t(*), nu(*) 
          real*4      fac
          real*4      f1, f2, Lsq1, Lsq2, dfreq
          integer*4   npts, nout, maxout, maxpts
          real*4      freq_MHz(npts)

          real*4      t_span
          !real*4      beg_nu, end_nu 
          real*4      d_nu, nu_span
          real*4      omega, atemp
          real*4      h_tmp, phi_tmp 

          real*4      cos_arr(maxout,maxpts), 
     -                sin_arr(maxout,maxpts)

             

C COUNTERS:
          integer*4   i,j,kk,nzero,n_positive

C CONSTANTS:
          real*4      pi,twopi


          pi = 3.14159265358979 !acos(-1.0d0)
          twopi = 6.28318530717959 !2.0d0*pi
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
          !nu_span = dble(npts-1)*d_nu
          nu_span = dble(npts)*d_nu
c_previous          beg_nu = -0.5d0*nu_span + 0.5d0*d_nu
c_previous          end_nu = 0.5d0*nu_span - 0.5d0*d_nu
c_previous          h_tmp = (end_nu - beg_nu)/dble(nout-1)
c_previous          do i = 1,nout   
c_previous             nu(i) = beg_nu + dble(i-1)*h_tmp
c_previous          end do

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
             omega = 2.0*nu(i)  ! assuming nu_trial to be angular frequency
             do kk = 1,npts ! number of t's
                phi_tmp = omega*t(kk)
                cos_arr(i,kk) = cos(phi_tmp)
                sin_arr(i,kk) = -sin(phi_tmp)
             enddo
          enddo
C -------------------------------------------------          
          end
