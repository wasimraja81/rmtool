chelp+
      ! -------------------------------------------- 
      ! Subroutine to clean the DIRTY RM spectra. 
      ! This, in its present form is the crudest 
      ! clean algorithm without any optimisation. 
      ! 
      ! Cleans Q and U separately. 
      ! 
      ! It requires as input : 
      ! 1)               L_sq --> The Lambda-squared array; real array
      ! 2)              nchan --> Number of spectral channels 
      ! 3)                  W --> Weights for the spectral channels;
      !                           real array
      ! 4)              RM_in --> Sampled RM's at which to compute the
      !                           response function; real array 
      ! 5)       DtyMAP_RMQ/U --> Dirty Map Q/U; real array
      ! 6)             max_rm --> amximum dimension of rm-array; needed 
      !                           to take care of definition in subroutines
      ! 7)                nrm --> integer number of RM-samples 
      ! 8)              niter --> integer number of clean cycles 
      ! 9)               gain --> clean loop gain
      !10)             thresh --> threshold for stopping clean 
      !11)OUT    ClnMAP_RMQ/U --> The cleaned map Q/U; real array
      !12)               ofac --> Oversampling factor for fourier interpolation
      !                           needed to detect precise location of peak in 
      !                           the RM-spectrum; integer
      !13)        interp_type --> The method of interpolation to be used; integer
      ! 
      !                       --wr, 08 Nov, 2010
      ! -------------------------------------------- 
chelp-


      subroutine rm_clean_Q(L_sq, nchan,W, RM_in, 
     -                    DtyMAP_RMQ, maxrm, nrm, 
     -                    niter, gain, thresh,
     -                    ClnMAP_RMQ, ofac, interp_type)

      implicit none

      integer*4  maxrm
      integer*4  nchan, nrm, niter, nbeams 
      real*4     L_sq(*), W(*), RM_in(*) 
      real*4     DtyMAP_RMamp(*), DtyMAP_RMpha(*) 
      real*4     ClnMAP_RMamp(*), ClnMAP_RMpha(*)
      real*4     DtyBEAM_RMQ(maxrm,maxrm)
      real*4     re_beam(maxrm)
      integer*4  ClnComp(nrm)
      real*4     ClnFluxQ(nrm), ResiQ(nrm)
      real*4     frac 
      real*4     gain, thresh 
      real*4     beam_now(nrm), avg_absMAP, rms_MAP
      real*4     peak_val, peak_loc  
      real*4     RM1, RM2, stored_peak(maxrm) 
      real*4     subarr_re(64)
      integer*4  ofac
      real*4     dh 

      integer*4  i, imax, iter, bm_indx, n_subarr, interp_type 
      logical    over
      ! Plotting related:
      real*4     xarr(nrm), yarr(nrm)
      character*120 xlabel, ylabel, title
      character*1 yorn


      write(*,*)"In development..."

      over = .false.
      iter = 0
      nbeams = 0
      ! Initialise the clean component aaray as well 
      ! as the Residual array: 
      do i = 1,nrm
         ClnComp(i) = 0

         ClnFluxQ(i) = 0.0
         ResiQ(i) = DtyMAP_RMQ(i)
      enddo

      ! Start cleaning Q first: 
      do while(.not.over)
         iter = iter + 1
         ! Locate the peak of the |amplitude|: 
         !call meanrms(ResiAmp,rms_MAP,avg_absMAP,nrm)
         call rms(ResiQ,nrm,rms_MAP)
         call index_absmax(ResiQ,nrm,imax,avg_absMAP)
         ! Get the "exact" location of the peak using 
         ! interpolation. This is necessary since the 
         ! beam itself is RM-dependent: 

         subarr_re(1) = ResiQ(imax-2)
         subarr_re(2) = ResiQ(imax-1)
         subarr_re(3) = ResiQ(imax)
         subarr_re(4) = ResiQ(imax+1)
         subarr_re(5) = ResiQ(imax+2)

         n_subarr = 5
         RM1 = RM_in(imax-2)
         RM2 = RM_in(imax+2)

         if(interp_type.eq.0)then ! Peak lies at sampled RM
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 dh = (RM2 - RM1)/real(nrm-1)
                 peak_loc = RM_in(1) + real(imax-1)*dh
                 peak_val = ResiAmp(imax)
                 !peak_val = abs(ResiAmp(imax))
         else  ! Use entire Residual for peak-detection
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 do i = 1,nrm
                    tmp_arr(i) = 0.0
                 enddo
                 narr = nrm
                 call fourier_interp_re(ResiQ,tmp_arr,narr,ofac,nout)
                 call index_absmax(tmp_arr,nout,jmax)

                 dh_in = (RM2 - RM1)/real(narr-1)
                 dh_out = dh_in*real(narr)/real(nout)
                 RM0 = 0.5*(RM1 + RM2) ! mid RM
   
                 ! Locate the 1st and last pts in the resampled 
                 ! abscissa: 
                 band_start = RM0 - 0.5*dh_in*real(narr)
                 band_stop = RM0 + 0.5*dh_in*real(narr)
   
                 samp1 = band_start + 0.5*dh_out
                 samp2 = band_stop - 0.5*dh_out
   
                 ! Fix the abscissa to include the midpoint:
                 if(mod(nout,2).eq.0)then
                     RMshift = RM0-(samp1 + real(nout/2)*dh_out)
                     samp1 = samp1-RMshift+0.5*(dh_in-dh_out)-dh_out
                     samp2 = samp2-RMshift+0.5*(dh_in-dh_out)-dh_out
                 else
                     write(*,*)"Warning: nout is odd!!"
                     write(*,*)"We want nout to be even"
                     write(*,*)"Quitting now..."
                     stop
                 endif
                 peak_loc = samp1 + real(jmax - 1)*dh_out
                 peak_val = tmp_arr(jmax)
         endif


         ! Get the RM-dependent beam:
         ! See if the beam has been computed already: 
         do i = 1,nbeams
            if(peak_loc.eq.stored_peak(i))then
                    bm_indx = i
                    
                    goto 1999
            endif
         enddo
         ! Compute the "phase" for the Qbeam from the 
         ! dirty map:
         call get_phaseQ(peakloc,L_sq,ResiQ,nchan,phase_val)
         call compute_dirty_rmbeamQ(L_sq,W,nchan,peak_loc, 
     -                             RM_in, phase_val, nrm, 
     -                             re_beam, maxrm)

         nbeams = nbeams + 1
         write(*,*)"nbeams: ",nbeams
         stored_peak(nbeams) = peak_loc
         do i = 1,nrm
            DtyBEAM_RMQ(nbeams,i) = re_beam(i)
         enddo
         bm_indx = nbeams

1999     continue
         do i = 1,nrm
            beam_now(i) = abs(DtyBEAM_RMQ(bm_indx,i))
         enddo

         !peak_val = abs(ResiAmp(imax))
         !peak_val = ResiAmp(imax) ! Forces stop at first encounter of 
                                  ! negative, if cleaning is done using 
                                  ! residual amplitude obtained from 
                                  ! subtracting gain x peak x amplitude.
                                  ! Irrelevant if resiQ and resiU is 
                                  ! used to derive the resiAmp
         !if(peak_val.gt.thresh*avg_absMAP.and.iter.lt.niter)then
         if(abs(peak_val).gt.thresh*rms_MAP.and.iter.lt.niter)then
                 frac = gain*peak_val 

                 ! Accumulate the Clean Components at peak location:
                 ClnComp(imax) = ClnComp(imax) + 1
                 ClnFluxQ(imax) = ClnFluxQ(imax) + frac 

                 ! Obtain the residual Map: 
                 do i = 1,nrm 
                    ResiQ(i) = ResiQ(i) - frac*DtyBEAM_RMQ(bm_indx,i)
                 enddo 
                 ! Some plotting 
                 call pgbeg(0,'/xs',1,1)
                 xlabel = 'RM bin num'
                 ylabel = 'Residual Q'
                 title = 'Diagnosis of RMclean'
                 do i = 1,nrm
                    xarr(i) = real(i)
                    yarr(i) = real(ResiQ(i))
                 enddo
                 call myplot1(xarr,yarr,nrm,xlabel,ylabel,title,2)
                 call pgend
                 do i = 1,nrm
                    xarr(i) = real(i)
                    yarr(i) = real(ClnFluxQ(i))
                 enddo
                 write(*,*)"Iter No.  Max Loc      peak_loc    peak_val"
                 write(*,*)"   ",iter,"      ",imax,"     ",
     -                    peak_loc," ", peak_val
                 write(*,*)" "
                 write(*,*)"Press 'S' to stop clean..."
                 write(*,*)"Any other key to continue clean..."
                 read(*,*)yorn
                 if(yorn.eq.'S'.or.yorn.eq.'s')then
                         goto 1001
                 endif
         else
                 ! Stop Clean and Restore the Cleaned Map: 
1001             continue
                 do i = 1,nrm
                    ClnFluxQ(i) = ClnFluxQ(i) + ResiQ(i) 
                 enddo
                 over = .true. 
                 write(*,*)"     niter: ",iter
                 write(*,*)" |max-res|: ",peak_val
                 write(*,*)" |avg_res|: ",avg_absMAP
                 write(*,*)" |rms_res|: ",rms_MAP
         endif
      enddo

      return

      end

      subroutine rm_cleanU(L_sq, nchan,W, RM_in, 
     -                    DtyMAP_RMU, maxrm, nrm, 
     -                    niter, gain, thresh,
     -                    ClnMAP_RMU, ofac, interp_type)

      implicit none

      integer*4  maxrm
      integer*4  nchan, nrm, niter, nbeams 
      real*4     L_sq(*), W(*), RM_in(*) 
      real*4     DtyMAP_RMU(*) 
      real*4     ClnMAP_RMU(*)
      real*4     DtyBEAM_RMU(maxrm,maxrm)
      real*4     re_beam(maxrm)
      integer*4  ClnComp(nrm)
      real*4     ClnFluxU(nrm), ResiU(nrm)
      real*4     frac 
      real*4     gain, thresh 
      real*4     beam_now(nrm), avg_absMAP, rms_MAP
      real*4     peak_val, peak_loc  
      real*4     RM1, RM2, stored_peak(maxrm) 
      real*4     subarr_re(64)
      integer*4  ofac
      real*4     dh 

      integer*4  i, imax, iter, bm_indx, n_subarr, interp_type 
      logical    over
      ! Plotting related:
      real*4     xarr(nrm), yarr(nrm)
      character*120 xlabel, ylabel, title
      character*1 yorn


      write(*,*)"In development..."

      over = .false.
      iter = 0
      nbeams = 0
      ! Initialise the clean component aaray as well 
      ! as the Residual array: 
      do i = 1,nrm
         ClnComp(i) = 0

         ClnFluxU(i) = 0.0
         ResiU(i) = DtyMAP_RMU(i)
      enddo

      ! Start cleaning Q first: 
      do while(.not.over)
         iter = iter + 1
         ! Locate the peak of the |amplitude|: 
         !call meanrms(ResiAmp,rms_MAP,avg_absMAP,nrm)
         call rms(ResiU,nrm,rms_MAP)
         call index_absmax(ResiU,nrm,imax,avg_absMAP)
         ! Get the "exact" location of the peak using 
         ! interpolation. This is necessary since the 
         ! beam itself is RM-dependent: 

         subarr_re(1) = ResiU(imax-2)
         subarr_re(2) = ResiU(imax-1)
         subarr_re(3) = ResiU(imax)
         subarr_re(4) = ResiU(imax+1)
         subarr_re(5) = ResiU(imax+2)

         n_subarr = 5
         RM1 = RM_in(imax-2)
         RM2 = RM_in(imax+2)

         if(interp_type.eq.0)then ! Peak lies at sampled RM
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 dh = (RM2 - RM1)/real(nrm-1)
                 peak_loc = RM_in(1) + real(imax-1)*dh
                 peak_val = ResiU(imax)
                 !peak_val = abs(ResiAmp(imax))
         else  ! Use entire Residual for peak-detection
                 RM1 = RM_in(1)
                 RM2 = RM_in(nrm)
                 do i = 1,nrm
                    tmp_arr(i) = 0.0
                 enddo
                 narr = nrm
                 call fourier_interp_re(ResiU,tmp_arr,narr,ofac,nout)
                 call index_absmax(tmp_arr,nout,jmax)

                 dh_in = (RM2 - RM1)/real(narr-1)
                 dh_out = dh_in*real(narr)/real(nout)
                 RM0 = 0.5*(RM1 + RM2) ! mid RM
   
                 ! Locate the 1st and last pts in the resampled 
                 ! abscissa: 
                 band_start = RM0 - 0.5*dh_in*real(narr)
                 band_stop = RM0 + 0.5*dh_in*real(narr)
   
                 samp1 = band_start + 0.5*dh_out
                 samp2 = band_stop - 0.5*dh_out
   
                 ! Fix the abscissa to include the midpoint:
                 if(mod(nout,2).eq.0)then
                     RMshift = RM0-(samp1 + real(nout/2)*dh_out)
                     samp1 = samp1-RMshift+0.5*(dh_in-dh_out)-dh_out
                     samp2 = samp2-RMshift+0.5*(dh_in-dh_out)-dh_out
                 else
                     write(*,*)"Warning: nout is odd!!"
                     write(*,*)"We want nout to be even"
                     write(*,*)"Quitting now..."
                     stop
                 endif
                 peak_loc = samp1 + real(jmax - 1)*dh_out
                 peak_val = tmp_arr(jmax)
         endif


         ! Get the RM-dependent beam:
         ! See if the beam has been computed already: 
         do i = 1,nbeams
            if(peak_loc.eq.stored_peak(i))then
                    bm_indx = i
                    
                    goto 1999
            endif
         enddo
         ! Compute the "phase" for the Qbeam from the 
         ! dirty map:
         call get_phaseU(peakloc,L_sq,ResiU,nchan,phase_val)
         call compute_dirty_rmbeamU(L_sq,W,nchan,peak_loc, 
     -                             RM_in, phase_val, nrm, 
     -                             re_beam, maxrm)

         nbeams = nbeams + 1
         write(*,*)"nbeams: ",nbeams
         stored_peak(nbeams) = peak_loc
         do i = 1,nrm
            DtyBEAM_RMU(nbeams,i) = re_beam(i)
         enddo
         bm_indx = nbeams

1999     continue
         do i = 1,nrm
            beam_now(i) = abs(DtyBEAM_RMU(bm_indx,i))
         enddo

         !peak_val = abs(ResiAmp(imax))
         !peak_val = ResiAmp(imax) ! Forces stop at first encounter of 
                                  ! negative, if cleaning is done using 
                                  ! residual amplitude obtained from 
                                  ! subtracting gain x peak x amplitude.
                                  ! Irrelevant if resiQ and resiU is 
                                  ! used to derive the resiAmp
         !if(peak_val.gt.thresh*avg_absMAP.and.iter.lt.niter)then
         if(abs(peak_val).gt.thresh*rms_MAP.and.iter.lt.niter)then
                 frac = gain*peak_val 

                 ! Accumulate the Clean Components at peak location:
                 ClnComp(imax) = ClnComp(imax) + 1
                 ClnFluxU(imax) = ClnFluxU(imax) + frac 

                 ! Obtain the residual Map: 
                 do i = 1,nrm 
                    ResiU(i) = ResiU(i) - frac*DtyBEAM_RMU(bm_indx,i)
                 enddo 
                 ! Some plotting 
                 call pgbeg(0,'/xs',1,1)
                 xlabel = 'RM bin num'
                 ylabel = 'Residual Q'
                 title = 'Diagnosis of RMclean'
                 do i = 1,nrm
                    xarr(i) = real(i)
                    yarr(i) = real(ResiU(i))
                 enddo
                 call myplot1(xarr,yarr,nrm,xlabel,ylabel,title,2)
                 call pgend
                 do i = 1,nrm
                    xarr(i) = real(i)
                    yarr(i) = real(ClnFluxU(i))
                 enddo
                 write(*,*)"Iter No.  Max Loc      peak_loc    peak_val"
                 write(*,*)"   ",iter,"      ",imax,"     ",
     -                    peak_loc," ", peak_val
                 write(*,*)" "
                 write(*,*)"Press 'S' to stop clean..."
                 write(*,*)"Any other key to continue clean..."
                 read(*,*)yorn
                 if(yorn.eq.'S'.or.yorn.eq.'s')then
                         goto 1001
                 endif
         else
                 ! Stop Clean and Restore the Cleaned Map: 
1001             continue
                 do i = 1,nrm
                    ClnFluxU(i) = ClnFluxU(i) + ResiU(i) 
                 enddo
                 over = .true. 
                 write(*,*)"     niter: ",iter
                 write(*,*)" |max-res|: ",peak_val
                 write(*,*)" |avg_res|: ",avg_absMAP
                 write(*,*)" |rms_res|: ",rms_MAP
         endif
      enddo

      return

      end








c         Author: Wasim Raja
c         Date: 12-05-2009
c         
c         This subroutine locates the "index" of maxima of 
c         the absolute values of the elements of a vector 
c         "V" given its length "N". It also returns the avg 
c         of the absolute values. 

          subroutine index_absmax(V,N,imax,avg_abs)

          implicit none
          real*4    V(*), avg_abs, maxx, absV 
          integer*4 N
          integer*4 i, imax

          maxx = abs(V(1))
          imax = 1
          avg_abs = 0.0d0
          do i = 1,N
               absV = abs(V(i))
               avg_abs = avg_abs + absV
               if(absV.ge.maxx)then
                       maxx = absV
                       imax = i 
               endif
          end do
          avg_abs = avg_abs/dble(N)

          return
          end 


        SUBROUTINE MEANRMS(A,RMS,AMEAN,NP)
C	TO COMPUTE MEAN AND RMS OF 'A' BY EXCLUDING WHAT MAY BE
C	SOME CONTRIBUTION FROM INTERFERENCE
C
        implicit none
        real*4 A(*)
        real*4 RMS,AMEAN, AN, DIFF, AMEAN0, RMS0  
        integer*4 NP, I, ITER  
C
        ITER=0
C
101     ITER=ITER+1
        AMEAN=0.0
        AN=0.0
C
        DO I=1,NP
         IF(ITER.EQ.1)GO TO 1
         DIFF=ABS(A(I)-AMEAN0)
         IF(DIFF.LE.(4.*RMS))GO TO 1
         GO TO 2
1        AMEAN=A(I)+AMEAN
         AN=AN+1.
2        CONTINUE
        END DO
C
        RMS0=RMS
        if(an.gt.0)AMEAN=AMEAN/AN
C
        RMS=0.0
        AN=0.0
        DO I=1,NP
         DIFF=ABS(A(I)-AMEAN)
         IF(ITER.EQ.1)GO TO 11
         IF(DIFF.LE.(4.*RMS0))GO TO 11
         GO TO 12
11       RMS=RMS+DIFF*DIFF
         AN=AN+1.
12       CONTINUE
        END DO

        end

