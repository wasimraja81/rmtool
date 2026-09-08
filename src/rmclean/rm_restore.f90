chelp+
      !-------------------------------------------------
      ! This subroutine restores the clean components
      ! after convolving the cleaned RM-spectra with
      ! a gaussian beam of FWHM matching the resolution
      ! in RM.
      !                             --wr, 07 Sep, 2011
      !-------------------------------------------------
chelp-

      subroutine rm_restore(L_sq, nchan, RM_amp_spec,RM_pha_spec,
     -                      RM_val, nrm)


         implicit none

         integer*4     nrm, nchan, iwidth, nzero, ipix
         real*4        L_sq(*), RM_amp_spec(*), RM_val(*), RM_pha_spec(*)
         real*4        dfreq, f1, f2, Lsq1, Lsq2,
     -                 Lsq_span, FWHM_RM, d_RM,
     -                 fac, pi
         real*4        gaus_template(nrm)
         real*4        sigma, two_sigma_sq, norm
         real*4        yarr(nrm)
         real*4        atmp, ptmp

         ! TEST
         !character*120 xlabel, ylabel, title



         pi = 3.14159265358979 !acos(-1.0d0)
         fac = 3.14159265358979 !acos(-1.0d0)

         ! Determine the width of the Gaussian
         ! with which to convolve the RM-spec:

         f1 = 300.0d0/sqrt(L_sq(nchan))
         f2 = 300.0d0/sqrt(L_sq(1))

         ! Calculate the edge L_sq:
         dfreq = (f2 - f1)/dble(nchan-1)

         f1 = f1 - 0.5d0*dfreq
         f2 = f2 + 0.5d0*dfreq

         Lsq2 = (300.0d0/f1)**2
         Lsq1 = (300.0d0/f2)**2

         ! Relation between the 2 domains:
         Lsq_span = Lsq2 - Lsq1
         FWHM_RM = 0.5*fac/Lsq_span  ! The factor fac is pi for the cases like
         ! RM and Lambda**2 kind of extraction
         ! For nu vs t kind of extraction, fac = 1.0
         d_RM = (RM_val(nrm) - RM_val(1))/dble(nrm - 1)
         sigma = (0.42466*FWHM_RM)  ! Width of the Gaussian
         iwidth = int(0.42466*FWHM_RM/d_RM)+1  ! Width of the Gaussian in pixels
         ! TEST:
         !write(*,*)"FWHM: ",FWHM_RM
         !write(*,*)"d_RM: ",d_RM
         !write(*,*)"iwid: ",iwidth
         !stop

         ! Locate the pixel containing RM = 0.0:
         if(mod(nrm,2).eq.0)then
            nzero = nrm/2 + 1
         else
            nzero = (nrm + 1)/2
         endif

         ! Generate the convolving template (centered at 0):
         two_sigma_sq = 2.0*sigma*sigma
         !norm = 1.0/sqrt(two_sigma_sq*pi) ! for normalised Gaussian
         norm = 1.0                        ! For Gaussian of height 1

         ! TEST:
         !write(*,*)"        nzero: ",nzero
         !write(*,*)"         norm: ",norm
         write(*,*)"         FWHM: ",FWHM_RM
         write(*,*)"         d_RM: ",d_RM
         write(*,*)"       iwidth: ",iwidth
         write(*,*)" two_sigma_sq: ",two_sigma_sq

         do ipix = 1,nrm
            gaus_template(ipix) = norm*exp(-(RM_val(ipix))**2/two_sigma_sq)
            write(88,*)gaus_template(ipix)
         enddo
         ! TEST:
         !call pgbeg(0,'5/xs',0,0)
         !xlabel = 'RM '
         !ylabel = 'Power '
         !title = 'TEST from SUBROUTINE '
         !call myplot1(RM_val,gaus_template,nrm,xlabel,ylabel,title,2)
         !call pgend

         ! Perform the convolution (Using FFT):
         do ipix = 1,nrm
            atmp = RM_amp_spec(ipix)
            ptmp = RM_pha_spec(ipix)
            RM_amp_spec(ipix) = atmp*cos(ptmp)  ! Real Part
            RM_pha_spec(ipix) = atmp*sin(ptmp)  ! Imag Part
         enddo

         call fft1d(RM_amp_spec,RM_pha_spec,nrm)
         call fftshift1d(RM_amp_spec,nrm)
         call fftshift1d(RM_pha_spec,nrm)

         do ipix = 1,nrm
            yarr(ipix) = 0.0
         enddo
         call fft1d(gaus_template,yarr,nrm)
         call fftshift1d(gaus_template,nrm)
         call fftshift1d(yarr,nrm)

         ! Do the X-ing in the F-domain:
         do ipix = 1,nrm
            atmp = RM_amp_spec(ipix)
            ptmp = RM_pha_spec(ipix)

            RM_amp_spec(ipix) = atmp*gaus_template(ipix)
     -                   + ptmp*yarr(ipix)
            RM_pha_spec(ipix) = ptmp*gaus_template(ipix)
     -                   - atmp*yarr(ipix)
         enddo

         call ifft1d(RM_amp_spec,RM_pha_spec,nrm)
         call fftshift1d(RM_amp_spec,nrm)
         call fftshift1d(RM_pha_spec,nrm)

         do ipix = 1,nrm
            atmp = RM_amp_spec(ipix)
            ptmp = RM_pha_spec(ipix)

            RM_amp_spec(ipix) = sqrt(atmp**2 + ptmp**2)
            RM_pha_spec(ipix) = atan2(ptmp,atmp)
         enddo


      end
