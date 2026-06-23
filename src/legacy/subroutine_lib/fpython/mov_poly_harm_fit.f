chelp+
!  This is a wrapper routine meant to fit data sets
!  exceeding the "maxpts" parameter of the routine
!  "poly_harm_fit.f". Currently it is assumed that 
!  this parameter maxpts = 4096 in poly_harm_fit.f 
!  routine. 
!  It uses "poly_harm_fit" to perform a "moving fit" 
!  of the entire data. Each call to "poly_harm_fit" 
!  will fit a section of length N (N < maxpts). 
!  Two successive sections will have (N-n_stagger) 
!  data points in common. 
!  The common data is fitted using the mean of the 
!  fit values of the two successive regions. 
!  The total number of calls to "poly_harm_fit" will 
!  be such that:
!           ncalls = 1 + floor[ (ntot-N)/nstag ]
!
!  USAGE:      call mov_poly_harm_fit(InArr,npts,flag_arr,
!              nsamp_per_fit,n_stagger,order_p,order_h,thresh,
!              OutArr)
!      --wasim raja, 05 Jan, 2010
!  
!  CODES REQUIRED: 
!      1) poly_harm_fit.f
chelp-

        subroutine mov_poly_harm_fit(InArr,npts,flag_arr,
     -         nsamp_per_fit,n_stagger,order_p,order_h,thresh,
     -         OutArr)
        ! -----------------------------------------------------------
        implicit none

        integer*4 maxpts, maxpfit
        parameter(maxpts=65536)
        parameter(maxpfit=128)
        integer*4, intent(in) :: npts, order_p, order_h
        integer*4, intent(in) :: nsamp_per_fit, n_stagger
        real*4,    intent(in) :: thresh
        real*4,    intent(in), dimension(npts) :: InArr
        real*4,    intent(out),dimension(npts) :: OutArr
        real*4,    dimension(npts) :: flag_arr
        integer*4 order_p_now, order_h_now

        real*4 xarray(maxpts), yarray(maxpts)
        real*4 fit_array(maxpts)
        real*4 best_fit_param(maxpfit)

        integer*4 i, jj, k
        character*128 data_tag
        logical silent, lr_exclude
        integer*4 nright, nleft
        real*4 thresh_now

        integer*4 nsamp_per_fit_now, n_stagger_now
        integer*4 n_rem, n_calls, n_used
        
        external poly_harm_fit


        ! My way of preventing modification of 
        ! variables INSIDE the subroutine:
        nsamp_per_fit_now = nsamp_per_fit
        n_stagger_now = n_stagger
        order_p_now = order_p
        order_h_now = order_h
        thresh_now = thresh

        ! take care of non-senses:
        !if(n_stagger_now.le.0.or.npts.gt.4096)then
        if(n_stagger_now.le.0.or.nsamp_per_fit_now.gt.maxpts)then
                nsamp_per_fit_now = maxpts
                n_stagger_now = int(nsamp_per_fit/2)
                write(*,*)'----------------------------- '
                write(*,*)'We will fit ',maxpts,' pts at a time '
                write(*,*)'because either n_stagger is  '
                write(*,*)'invalid or npts exceeds maxpts '
                write(*,*)'allowed by poly_harm_fit!! '
                write(*,*)'This may not yield acceptable '
                write(*,*)'results... '
                write(*,*)'----------------------------- '
        endif

        !------------------------------------------------------------
        ! Set up the variables for poly_harm_fit routine:
        nleft = 1
        nright = nsamp_per_fit_now
        lr_exclude = .true. !if true,include l-to-r data,else exclude
        silent = .true.
        data_tag = ' '
        if(thresh_now.lt.3.0)thresh_now = 3.0

        !do j = 1,nsamp_per_fit
        !   xarray(j) = j
        !enddo


        n_calls = 1 + int((npts-nsamp_per_fit_now)/n_stagger_now)
        n_used = nsamp_per_fit_now + (n_calls - 1)*n_stagger_now
        !n_rem = npts - nsamp_per_fit - (n_calls - 1)*n_stagger
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
           xarray(k) = real(jj)
           fit_array(k) = flag_arr(jj)
        enddo
!       write(*,*)'I am in mov-fit... '
!       write(*,*)'npts: ',npts
!       write(*,*)'nleft: ',nleft
!       write(*,*)'nright: ',nright
!       write(*,*)'nfit: ',nsamp_per_fit
!       write(*,*)'nsta: ',n_stagger
!       write(*,*)'or_p: ',order_p
!       write(*,*)'or_h: ',order_h
!       write(*,*)'thre: ',thresh
!       write(*,*)'ncalls: ',n_calls
!       write(*,*)'nused: ',n_used
!       write(*,*)'nrem: ',n_rem
!       write(*,*)'I am in mov-fit... '
         
        call poly_harm_fit(nsamp_per_fit_now,xarray,yarray,
     -nleft,nright,lr_exclude,order_h_now,order_p_now,thresh,
     -silent,data_tag,fit_array,best_fit_param)

        do jj = 1,nsamp_per_fit_now
           OutArr(jj) = fit_array(jj)
        enddo

        !--------------------------------------------
        if(n_calls.lt.2)then
                goto 100
        endif
        ! Now the subsequent calls:
        do i = 2,n_calls ! refers to fitting call. 
           !write(*,*)"Call number: ",i
           k = 0 
           do jj = 1+(i-1)*n_stagger_now,nsamp_per_fit_now + 
     -                        (i-1)*n_stagger_now
              k = k + 1
              yarray(k) = InArr(jj)
              xarray(k) = real(jj)
              fit_array(k) = flag_arr(jj)
           enddo
           call poly_harm_fit(nsamp_per_fit_now,xarray,yarray,
     -nleft,nright,lr_exclude,order_h_now,order_p_now,thresh_now,
     -silent,data_tag,fit_array,best_fit_param)
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
     -               0.5*(OutArr(jj+(i-1)*n_stagger_now)+fit_array(jj))
           enddo
           do jj = 1+nsamp_per_fit_now-n_stagger_now,nsamp_per_fit_now
              OutArr(jj+(i-1)*n_stagger_now)=fit_array(jj)
           enddo
           ! ----------------------------------------
        enddo
100     continue
        if(n_rem.gt.0)then
                nright = n_rem
                ! fit the remaining data points at the end:
                ! Depending on the number of remaining pts
                ! this section of the data may contain a 
                ! large number of bad points -- eg., in the 
                ! case of a triangular-train, if this last
                ! section contains only as many points close
                ! to the turning-point, all such points will
                ! be marked as BAD -- because we do not wish 
                ! to fit the high-frequencies caused by the 
                ! points at the corners. If such a situation
                ! is encountered, we will simply NOT fit:
                k = 0
                do jj = 1+n_used,n_rem+n_used
                   k = k + 1                 
                   yarray(k) = InArr(jj)
                   xarray(k) = real(jj)
                   fit_array(k) = flag_arr(jj)
                enddo
                ! count the number of bad points:
                k = 0
                do jj = 1,n_rem
                   if(fit_array(jj).eq.0.0)then
                           k = k + 1
                   endif
                enddo
                !write(*,*)'ngood pts for fitting last sec:  ',n_rem-k
                if(n_rem-k.lt.min(order_p_now,order_h_now))then
                        ! do no fit
                        k = 0
                        do jj = 1+n_used,n_rem+n_used
                           k = k+1
                           fit_array(k) = InArr(jj)
                           OutArr(jj)=fit_array(k)
                        enddo
                        return
                endif
                call poly_harm_fit(n_rem,xarray,yarray,
     -nleft,nright,lr_exclude,order_h_now,order_p_now,thresh_now,
     -silent,data_tag,fit_array,best_fit_param)
                !--------------------------------------------
                k = 0
                do jj = 1+n_used,n_rem+n_used
                   k = k+1
                   OutArr(jj)=fit_array(k)
                enddo
        endif

        return
        end

      !include '/usr/lib/subroutine_lib/nchar.f'
      !include 'poly_harm_fit.f'
