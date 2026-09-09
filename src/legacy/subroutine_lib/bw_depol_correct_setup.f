chelp+
      !---------------------------------------------------
      ! This subroutine is designed to pre-compute certain 
      ! RM-independent variables to be passed to the code 
      ! "bw_depol_correct". 
      !                           --wr, 12 Sep, 2011
      !---------------------------------------------------
chelp-

      subroutine bw_depol_correct_setup(InArr, nchan, lsq_passed,
     -                                  z0, z1, z2 ) 


      implicit none

      !--------------------------------------
      real*4       InArr(*) 
      integer*4    nchan  
      real*8       L_sq(nchan)
      integer*4    i, j   
      real*8       f1, f2  
      real*8       z0(*), z1(*), z2(*) 
      real*8       delta_f 
      !--------------------------------------
      logical      lsq_passed   ! If true, we assume InArr array 
                                ! is passed in ascending order of 
                                ! Lambda-squared; we will resort 
                                ! it here. 
                                ! Else, we assume that the sampled 
                                ! center-frequencies are passed in 
                                ! the InArr array, and we shall 
                                ! compute from it the L_sq array 
                                ! sorted in the descending order 
                                ! as required by this code 



      !---------------------------------------------------------
      ! Sort issues with the L_sq array: 
      if(Lsq_passed)then
              j = nchan + 1
              do i = 1,nchan
                 j = j - 1
                 L_sq(j) = dble(InArr(i))  ! Resort L_sq values in
                                           ! descending order
              enddo
              ! Now generate the InArr array in 
              ! ascending order of frequencies: 
              do i = 1,nchan
                 InArr(i) = real(300.0/sqrt(L_sq(i)))
              enddo
      else
              do i = 1,nchan
                 L_sq(i) = (300.0d0/InArr(i))**2
              enddo
      endif
      !---------------------------------------------------------
      f1 = dble(InArr(1)) 
      f2 = dble(InArr(nchan)) 
      ! f1 and f2 are the frequencies corresponding to 
      ! the centres of the edge bins. Hence f2 - f1 
      ! will have (nchan - 1) number of bins and so: 
      delta_f = (f2 - f1)/dble(nchan-1)

      ! Store the RM-independent quantities: 
      do i = 1,nchan
         ! L_sq(i) = (300.0/InArr(i))**2
         z0(i) = L_sq(i)
         z1(i) = (300.0d0/(InArr(i) + 0.5d0*delta_f))**2 - z0(i)
         z2(i) = (300.0d0/(InArr(i) - 0.5d0*delta_f))**2 - z0(i)
      enddo

      return

      end

      
