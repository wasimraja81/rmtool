chelp+
!  This is a wrapper routine meant to compute moving 
!  average of data. 
!  It uses "mean" to perform a "moving average" 
!  of the entire data. Each call to "moving_average" 
!  will fit a section of length N (N < maxpts). 
!  Two successive sections will have (N-n_stagger) 
!  data points in common. 
!  The common data is fitted using the mean of the 
!  fit values of the two successive regions. 
!  The total number of calls to "mean" will 
!  be such that:
!           ncalls = 1 + floor[ (ntot-N)/nstag ]
!
!  USAGE:      call moving_average(InArr,npts,
!              nsamp_per_fit,n_stagger,OutArr)
!      --wasim raja, 05 Jan, 2010
!  
!  CODES REQUIRED: 
!      1) fort_lib.f
chelp-

        subroutine moving_average(InArr,npts,
     -         nsamp_per_fit,n_stagger,OutArr)
        ! -----------------------------------------------------------
        implicit none

        integer*4 maxpts 
        parameter(maxpts=8192*2)
        integer*4 npts 
        real*4 InArr(*)
        real*4 OutArr(*)

        real*4 yarray(maxpts)
        real*4 mean_now
        !real*4 rms 

        integer*4 i, jj, k

        integer*4 nsamp_per_fit, n_stagger
        integer*4 nsamp_per_fit_now, n_stagger_now
        integer*4 n_rem, n_calls, n_used


        ! My way of preventing modification of 
        ! variables INSIDE the subroutine:
        nsamp_per_fit_now = nsamp_per_fit
        n_stagger_now = n_stagger

        ! take care of non-senses:
        if(n_stagger_now.le.0.or.npts.gt.4096)then
                nsamp_per_fit_now = 4096
                n_stagger_now = int(nsamp_per_fit/2)
                write(*,*)'----------------------------- '
                write(*,*)'We will fit 4096 pts at a time '
                write(*,*)'because either n_stagger is  '
                write(*,*)'invalid or npts exceeds maxpts '
                write(*,*)'allowed by poly_harm_fit!! '
                write(*,*)'This may not yield acceptable '
                write(*,*)'results... '
                write(*,*)'----------------------------- '
        endif

        !------------------------------------------------------------

        n_calls = 1 + int((npts-nsamp_per_fit_now)/n_stagger_now)
        n_used = nsamp_per_fit_now + (n_calls - 1)*n_stagger_now
        n_rem = npts - n_used

        ! at this stage we have the data to be fitted
        !--------------------------------------------
        ! First call to "poly_harm_fit" to fill the 
        ! first nsamp_per_fit number of elements in 
        ! the OutArr
        k = 0 
        do jj = 1,nsamp_per_fit_now
           k = k + 1
           yarray(k) = InArr(jj)
        enddo
         
        !call rob_mean(yarray,nsamp_per_fit_now,mean_now)
        !call meanrms(yarray,rms,mean_now,nsamp_per_fit_now)
        call median(yarray,nsamp_per_fit_now,mean_now)

        do jj = 1,nsamp_per_fit_now
           OutArr(jj) = mean_now
        enddo

        !--------------------------------------------
        if(n_calls.lt.2)then
                goto 100
        endif
        ! Now the subsequent calls:
        do i = 2,n_calls ! refers to fitting call. 
           k = 0 
           do jj = 1+(i-1)*n_stagger_now,nsamp_per_fit_now + 
     -                        (i-1)*n_stagger_now
              k = k + 1
              yarray(k) = InArr(jj)
           enddo
           !call rob_mean(yarray,nsamp_per_fit_now,mean_now) 
           !call meanrms(yarray,rms,mean_now,nsamp_per_fit_now)
           call median(yarray,nsamp_per_fit_now,mean_now)
           ! retain the first n_stagger number 
           ! of points in OutArr as it was, replace 
           ! the next (nsamp_per_fit - nstagger) 
           ! elements by the mean of OutArr and the 
           ! first (nsamp_per_fit - n_stagger) 
           ! elements of fit_array, finally replace 
           ! the next n_stagger elements of OutArr
           ! by the last n_stagger elements of 
           ! fit_array
           do jj = 1,(nsamp_per_fit_now - n_stagger_now)
              OutArr(jj+(i-1)*n_stagger_now) = 
     -               0.5*(OutArr(jj+(i-1)*n_stagger_now)+mean_now)
           enddo
           do jj = 1+nsamp_per_fit_now-n_stagger_now,nsamp_per_fit_now
              OutArr(jj+(i-1)*n_stagger_now)= mean_now
           enddo
           ! ----------------------------------------
        enddo
100     continue
        if(n_rem.gt.0)then
                !call rob_mean(yarray, n_rem, mean_now)
                !call meanrms(yarray,rms,mean_now,n_rem)
                call median(yarray,n_rem,mean_now)
                !--------------------------------------------
                do jj = 1+n_used,n_rem+n_used
                   OutArr(jj)=mean_now
                enddo
        endif

        return
        end

      !include '/usr/lib/subroutine_lib/nchar.f'
      !include '/usr/lib/subroutine_lib/poly_harm_fit.f'
