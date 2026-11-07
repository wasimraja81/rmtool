chelp+
      ! -------------------------------------------- 
      ! Subroutine to clean the DIRTY RM spectra. 
      ! This, in its present form is the crudest 
      ! clean algorithm without any optimisation. 
      ! 
      ! It requires as input :
      ! 1)              nrm --> integer number
      ! 2)  Dirty RM-ampMAP --> array: nrm points
      ! 3)  Dirty RM-phaMAP --> array: nrm points
      !X4) Re-Dirty RM-beam --> array: (nrm x nrm) 
      !X5) Im-Dirty RM-beam --> array: (nrm x nrm) 
      ! 6)            niter --> integer number
      ! 7)           thresh --> real*4 
      ! 8)             gain --> real*4 
      ! 
      ! and outputs: 
      ! 9)  Clean RM-ampMAP --> array: nrm points
      !10)  Clean RM-phaMAP --> array: nrm points
      !
      !                       --wr, 08 Nov, 2010
      ! 
      !11)  Residual ampMAP --> array: nrm points
      !12)   Residual Q-MAP --> array: nrm points
      !13)   Residual U-MAP --> array: nrm points
      ! 
      !                       --wr, 02 Feb, 2014 
      ! -------------------------------------------- 
chelp-


      subroutine rm_clean(L_sq, nchan, RM_in, 
     -                    DtyMAP_RMamp, DtyMAP_RMpha, 
     -                    nrm, niter, gain, thresh,
     -                    ClnMAP_RMamp, ClnMAP_RMpha, 
     -                    ResiQ, ResiU, ResiAmp,
     -                    ofac,
     -                    interp_type,cos_arr,sin_arr,maxrm,maxchan,
     -                    bw_depol_correct,tn, z0,z1,z2,interactive,
     -                    FWHM_RM)

      implicit none

      integer*4  maxrm, maxchan 
      integer*4  nchan, nrm, niter 
      integer*4  tn 
      real*8     z0(*), z1(*), z2(*)
      real*4     L_sq(*), RM_in(*), FWHM_RM 
      real*4     DtyMAP_RMamp(*), DtyMAP_RMpha(*) 
      real*4     ClnMAP_RMamp(*), ClnMAP_RMpha(*)
!      real*4     DtyBEAM_reRM(maxrm,maxrm), DtyBEAM_imRM(maxrm,maxrm)
      real*4     re_beam(maxrm), im_beam(maxrm) 
      integer*4  ClnComp(nrm)
      real*4     ClnFlux(nrm), ResiAmp(*)
      real*4     ResiQ(*), ResiU(*),
     -           fracQ, fracU, ClnFluxQ(nrm),ClnFluxU(nrm)
      real*4     gain, thresh 
      real*4     beam_now(nrm), avg_absMAP, rms_MAP
      real*4     peak_val, peak_loc, phase_val, frac 
      real*4     RM1, RM2 
      real*4     subarr_re(64), subarr_im(64) 
      integer*4  ofac
      real*4     dh 
      real*4     cos_arr(maxrm,maxchan), sin_arr(maxrm,maxchan)

      integer*4  i, imax, iter, n_subarr, interp_type 
      logical    over
      ! Plotting related:
      real*4     xarr(nrm), yarr(nrm), xmax  
      character*120 xlabel, ylabel, title
      character*1 yorn
      logical*4   interactive, bw_depol_correct 
      real*4      pi, dev 



      pi = 3.1415927
      over = .false.
      iter = 0
      !TEST
      ! Initialise the clean component aaray as well 
      ! as the Residual array: 
      do i = 1,nrm
         ClnComp(i) = 0
         ResiAmp(i) = DtyMAP_RMamp(i)

         ClnFlux(i) = 0.0d0
         ClnFluxQ(i) = 0.0d0
         ClnFluxU(i) = 0.0d0

         ResiQ(i) = ResiAmp(i)*cos(DtyMAP_RMpha(i))
         ResiU(i) = ResiAmp(i)*sin(DtyMAP_RMpha(i))
      enddo

      do while(.not.over)
         iter = iter + 1
         ! Locate the peak of the |amplitude|: 
         call rms(ResiAmp,nrm,rms_MAP)
         call index_absmax(ResiAmp,nrm,imax,avg_absMAP)
         ! Get the "exact" location of the peak using 
         ! interpolation. This is necessary since the 
         ! beam itself is RM-dependent: 

         subarr_re(1) = ResiQ(imax-2)
         subarr_im(1) = ResiU(imax-2)

         subarr_re(2) = ResiQ(imax-1)
         subarr_im(2) = ResiU(imax-1)

         subarr_re(3) = ResiQ(imax)
         subarr_im(3) = ResiU(imax)

         subarr_re(4) = ResiQ(imax+1)
         subarr_im(4) = ResiU(imax+1)

         subarr_re(5) = ResiQ(imax+2)
         subarr_im(5) = ResiU(imax+2)

         n_subarr = 5
         RM1 = RM_in(imax-2)
         RM2 = RM_in(imax+2)

         if(interp_type.eq.0)then ! Peak lies at sampled RM
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 dh = (RM2 - RM1)/real(nrm-1)
                 peak_loc = RM_in(1) + real(imax-1)*dh
                 peak_val = ResiAmp(imax)
                 phase_val = atan2(ResiU(imax),ResiQ(imax))
         else if(interp_type.eq.-3)then ! Perform Fourier Interp using 
                                        ! entire Residual for peak-detection
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 call peak_interp(ResiQ,ResiU,nrm,RM1,RM2,ofac,
     -                            peak_loc, peak_val,phase_val, 
     -                            3)
         else
                 call peak_interp(subarr_re,subarr_im,n_subarr,RM1,RM2,
     -                            ofac,peak_loc, peak_val, phase_val, 
     -                            interp_type)
         endif
         phase_val = 0.5*phase_val  ! The "compute_dirty_beam" code
                                    ! needs PolPA and not the phase: 
                                    ! that's why 0.5*phase_val 

         ! interp_type: 1 => Parabolic interpolation
         ! interp_type: 2 => Fourier interpolation Amplitude
         ! interp_type: 3 => Fourier interpolation Complex
         ! interp_type: 4 => Sinc interpolation
         ! interp_type: ? => default -- Parabolic interpolation

         call compute_dirty_rmbeam(L_sq,nchan,peak_loc,phase_val, 
     -                             RM_in, nrm, re_beam, im_beam, 
     -                             cos_arr,sin_arr,maxrm,maxchan,
     -                             bw_depol_correct,tn,z0,z1,z2)

1999     continue
         do i = 1,nrm
            beam_now(i) = sqrt(re_beam(i)**2 + im_beam(i)**2)
         enddo

         peak_val = ResiAmp(imax) ! Forces stop at first encounter of 
                                  ! negative, if cleaning is done using 
                                  ! residual amplitude obtained from 
                                  ! subtracting gain x peak x amplitude.
                                  ! Irrelevant if resiQ and resiU is 
                                  ! used to derive the resiAmp
         !if(peak_val.gt.thresh*avg_absMAP.and.iter.lt.niter)then
         !if(peak_val.gt.thresh*rms_MAP.and.iter.lt.niter)then
         call maxima(ResiAmp,nrm,xmax)
         !! Suppressing comments to boost speed: 
         !write(*,*)"used peakval, actual_peakval: ",peak_val, xmax 
         !write(*,*)"                  avg_absMAP: ",avg_absMAP 
         !write(*,*)"                         dev: ",dev 
         !write(*,*)"                      thresh: ",thresh 
         !write(*,*)"                     rms_MAP: ",rms_MAP  
         !write(*,*)"              thresh*rms_MAP: ",thresh*rms_MAP  
         !!peak_val = xmax 
         dev = abs(peak_val - avg_absMAP) 
!         if(dev.gt.thresh*rms_MAP.and.iter.lt.niter.and.
!     -                            peak_val.gt.0.0)then
         if(dev.gt.thresh*rms_MAP.and.iter.lt.niter)then
                 frac = gain*peak_val 
                 !fracQ = frac*ResiQ(imax)
                 !fracU = frac*ResiU(imax)
                 fracQ = gain*ResiQ(imax)
                 fracU = gain*ResiU(imax)

                 ! Accumulate the Clean Components at peak location:
                 ClnComp(imax) = ClnComp(imax) + 1

                 ClnFluxQ(imax) = ClnFluxQ(imax) + fracQ
                 ClnFluxU(imax) = ClnFluxU(imax) + fracU
                 ClnFlux(imax) = ClnFlux(imax) + frac 
!                 ClnFlux(imax) = sqrt(ClnFluxQ(imax)**2 +  
!     -                                ClnFluxU(imax)**2) 

                 !write(18,*)frac,fracQ,fracU,ResiQ(imax),ResiU(imax)
                 ! Obtain the residual Map: 
                 do i = 1,nrm 
                    ResiQ(i) = ResiQ(i) - frac*re_beam(i)
                    ResiU(i) = ResiU(i) - frac*im_beam(i)
                    ResiAmp(i) = sqrt(ResiQ(i)*ResiQ(i) + 
     -                                ResiU(i)*ResiU(i)) 
!                    ResiAmp(i) = ResiAmp(i) - frac*beam_now(i) ! cause of negatives
                 enddo 
                 ! Some plotting if asked for: 
                 if(.not.interactive)then
                         goto 1002
                 endif
                 call pgbeg(0,'/xs',1,1)
                 !xlabel = 'RM bin num'
                 xlabel = 'RM (rad/m2)'
                 ylabel = 'Residual Amplitude'
                 title = 'Diagnosis of RMclean'
                 do i = 1,nrm
                    !xarr(i) = real(i)
                    xarr(i) = real(RM_in(i))
                    yarr(i) = real(ResiAmp(i))
                 enddo
                 call myplot1(xarr,yarr,nrm,xlabel,ylabel,title,2)
                 call pgend
                 do i = 1,nrm
                    xarr(i) = real(i)
                    yarr(i) = real(ClnFlux(i))
                 enddo
                 write(*,*)"-----------------------------------------"
                 write(*,*)"Iter No: ",iter, ", Max Loc: ", imax 
                 write(*,*)"   peak_loc: ",peak_loc, " rad/m2"
                 write(*,*)"   PosAngle: ",0.5*phase_val*180.0/pi," deg"
                 write(*,*)"    Max Amp: ",peak_val
                 write(*,*)"   Avg resi: ",avg_absMAP
                 write(*,*)"   RMS resi: ",rms_MAP 
                 write(*,*)"-----------------------------------------"
                 write(*,*)" "
                 write(*,*)"Press 'S' to stop clean..."
                 write(*,*)"Any other key to continue clean..."
                 read(*,'(a)')yorn
                 if(yorn.eq.'S'.or.yorn.eq.'s')then
                         goto 1001
                 endif
                 goto 1003 
1002             continue
                 !! Suppressing comments to boost speed: 
                 !write(*,*)"-----------------------------------------"
                 !write(*,*)"Iter No: ",iter, ", Max Loc: ", imax 
                 !write(*,*)"   peak_loc: ",peak_loc, " rad/m2"
                 !write(*,*)"   PosAngle: ",0.5*phase_val*180.0/pi," deg"
                 !write(*,*)"    Max Amp: ",peak_val
                 !write(*,*)"   Avg resi: ",avg_absMAP
                 !write(*,*)"   RMS resi: ",rms_MAP 
                 !write(*,*)"-----------------------------------------"
                 !write(*,*)" "
1003             continue 
         else
                 ! Stop Clean and Restore the Cleaned Map: 
1001             continue
                 do i = 1,nrm
                    ClnFluxQ(i) = ClnFluxQ(i) + ResiQ(i) 
                    ClnFluxU(i) = ClnFluxU(i) + ResiU(i) 
!                    ClnMAP_RMamp(i) = sqrt(ClnFluxQ(i)*ClnFluxQ(i) + 
!     -                                     ClnFluxU(i)*ClnFluxU(i))
                    ClnMAP_RMamp(i) = ClnFlux(i) + ResiAmp(i) 
                    ClnMAP_RMpha(i) = atan2(ClnFluxU(i),ClnFluxQ(i)) 
                    !ClnMAP_RMpha(i) = ResiAmp(i) ! TEST plotting -- remove later
                 enddo
                 over = .true. 
                 !! Suppressing comments to boost speed: 
                 !write(*,*)"     niter: ",iter
                 !write(*,*)" |max-res|: ",peak_val
                 !write(*,*)" |avg_res|: ",avg_absMAP
                 !write(*,*)" |rms_res|: ",rms_MAP
         endif
      enddo
      ! TEST: 
      ! write out the clean amp and phases: 
      do i = 1,nrm
         write(91,*)RM_in(i),ClnMap_RMamp(i),ClnMap_RMpha(i)*
     -                                     180./3.14159265358979
      enddo
      niter = iter 
      ! Restore the clean comp:
      call rm_restore(FWHM_RM,ClnMAP_RMamp,ClnMAP_RMpha,RM_in,nRM)
      !write(*,*)"In subroutine rm_clean.f: ",nrm 
      return

      end



