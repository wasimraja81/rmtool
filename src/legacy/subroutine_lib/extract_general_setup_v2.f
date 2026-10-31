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
!      
!      Provision is now made to generate the pre-computed 
!      values (cos_arr & sin_arr) for any arbitrary range 
!      of RMs. 
!      -- wasim, 14 Jan, 2014 


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
          integer*4   npts, nout, maxout, maxpts

          real*4      omega 
          real*4      phi_tmp 

          real*4      cos_arr(maxout,maxpts), 
     -                sin_arr(maxout,maxpts)

             

C COUNTERS:
          integer*4   i,kk 

C CONSTANTS:
          real*4      pi,twopi
          real*4      t_ref


          pi = 3.14159265358979 !acos(-1.0d0)
          twopi = 6.28318530717959 !2.0d0*pi

          call mean(t,npts,t_ref)
          t_ref = 0.0 
!          !order_taylor = 2

          !-------------------------------------------------
          do i = 1,nout ! number of nu's
             omega = 2.0*nu(i)  ! assuming nu_trial to be angular frequency
             do kk = 1,npts ! number of t's
                !phi_tmp = omega*t(kk)
                phi_tmp = omega*(t(kk) - t_ref) 
                cos_arr(i,kk) = cos(phi_tmp)
                sin_arr(i,kk) = -sin(phi_tmp)
             enddo
          enddo
C -------------------------------------------------          
          end
