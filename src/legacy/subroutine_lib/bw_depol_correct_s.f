chelp+
      !---------------------------------------------------
      ! This code RETURNS the complex "correction factor" 
      ! needed to correct the effect of inherent smearing 
      ! of the linear polarised Stokes intensities due 
      ! to the FINITE channel-width in the spectral axis. 
      !                           --wr, 12 Oct, 2010
      !---------------------------------------------------
chelp-

      subroutine bw_depol_correct(freq_MHz, nchan,
     -                           RM_in, nrm, tn, weight_amp, 
     -                           weight_pha)


      implicit none

      !--------------------------------------
      real*4         freq_MHz(*), RM_in(*) 
      integer*4      nchan, nrm, tn 
      real*4         weight_amp(nrm,nchan), weight_pha(nrm,nchan) 
      real*4         L_sq(nchan)
      integer*4      i, j, ii  
      real*4         RM, f1, f2  
      real*4         z0(nchan), z1(nchan), z2(nchan), 
     -               zz0, zz1, zz2 
      real*4         beta  
      real*4         delta_f 
      real*4         reaccum, imaccum, ddxn_fx, reval1, reval2, 
     -               imval1, imval2, norm(nchan) 
      integer*4      order_now 
      integer*8      factorial 
      !--------------------------------------
      ! Constants: 
      real*4         vel_light
      real*4         pi

      !---------------------------------------
      real*4         tmp_num 



      vel_light = 3.0e8
      pi = acos(-1.0d0)

      f1 = freq_MHz(1) 
      f2 = freq_MHz(nchan) 
      ! f1 and f2 are the frequencies corresponding to 
      ! the centres of the edge bins. Hence f2 - f1 
      ! will have (nchan - 1) number of bins and so: 
      delta_f = (f2 - f1)/dble(nchan-1)

      ! Apply the corrections to the smoothed Q and U: 
      ! beta is the index of lambda-squared appearing 
      ! due to the conversion from freq-to-lambda-squared.
      ! 
      beta = -1.5 !-3.0/2.0

      ! Store the RM-independent quantities: 
      do i = 1,nchan
         L_sq(i) = (300.0/freq_MHz(i))**2
         z0(i) = L_sq(i)
         z1(i) = (300.0/(freq_MHz(i) + 0.5*delta_f))**2 - z0(i)
         z2(i) = (300.0/(freq_MHz(i) - 0.5*delta_f))**2 - z0(i)

         ! Normalization factor (for RM=0)
         zz0 = z0(i)
         zz1 = z1(i)
         zz2 = z2(i)
         order_now = 0
         tmp_num = 0.0d0
         call val_integral(zz2, tmp_num, order_now, reval2, imval2)
         call val_integral(zz1, tmp_num, order_now, reval1, imval1)
         call ddxn_xpowalpha(zz0,beta,order_now,ddxn_fx)
         reaccum = ddxn_fx*(reval2 - reval1) 
         imaccum = ddxn_fx*(imval2 - imval1) 
         norm(i) = sqrt(reaccum*reaccum + imaccum*imaccum)
      enddo

      ! Obtain the weight factors (complex) for the matched 
      ! filter for each RM and each single channel: 
      do ii = 1,nrm
         RM = RM_in(ii) 
         do i = 1,nchan
            ! Term for norder = 0
            order_now = 0
            zz0 = z0(i)
            zz1 = z1(i)
            zz2 = z2(i)

            call val_integral(zz2, RM, order_now, reval2, imval2)
            call val_integral(zz1, RM, order_now, reval1, imval1)
            call ddxn_xpowalpha(zz0,beta,order_now,ddxn_fx)
            reaccum = ddxn_fx*(reval2 - reval1) 
            imaccum = ddxn_fx*(imval2 - imval1) 
   
            factorial = 1

            do j = 1,tn
               order_now = j
               call val_integral(zz2, RM, order_now, reval2, imval2)
               call val_integral(zz1, RM, order_now, reval1, imval1)
   
               call ddxn_xpowalpha(zz0,beta,order_now,ddxn_fx)
               factorial = factorial*j
   
               reaccum = reaccum + 
     -                          ddxn_fx*(reval2-reval1)/dble(factorial) 
               imaccum = imaccum + 
     -                          ddxn_fx*(imval2-imval1)/dble(factorial) 
            enddo
            weight_amp(ii,i) = sqrt(reaccum*reaccum + imaccum*imaccum)/
     -                                                          norm(i)
            weight_pha(ii,i) = atan2(imaccum,reaccum)

            !! Remove the phase introduced by smearing:
            !imaccum = atan2(U_in(i),Q_in(i)) ! phase
            !reaccum = sqrt(Q_in(i)*Q_in(i)+U_in(i)*U_in(i))
            !Q_corr(i) = reaccum*cos(imaccum - ph_weight)/weight
            !U_corr(i) = reaccum*sin(imaccum - ph_weight)/weight
         enddo
      enddo
      return

      end

      
      !---------------------------------------------------------
      !SUBROUTINE:   val_integral(z, RM, norder)
      !
      !              This is specific to only an integral 
      !              of the type:
      !              -- 
      !              |          
      !              | (z**n) exp(i RM z)dz
      !              |
      !             --
      ! 
      ! The solution is obtained for the real and imaginary parts 
      ! separately for reasons specific to the case of Q and U 
      ! smearing that takes place across the finite frequency 
      ! channel-widths.

      subroutine   val_integral(z, RM, norder, reval, imval)
      implicit none
      real*4    RM, z,  
     -          sin_RMz_by_RM, cos_RMz_by_RM, 
     -          reval, imval, Ic_0, Is_0, 
     -          Ic_now, Is_now, Ic_pre, Is_pre 
      integer*4 i, norder

      if(RM.eq.0.0)then
              Ic_0 = (z**(norder+1))/dble(norder+1)
              Is_0 = 0.0
      else
              Ic_0 = sin(RM*z)/RM
              Is_0 = -cos(RM*z)/RM
      endif

      Ic_now = Ic_0
      Is_now = Is_0
      if(norder.gt.0.and.abs(RM).gt.0)then
              sin_RMz_by_RM = sin(RM*z)/(RM)
              cos_RMz_by_RM = cos(RM*z)/(RM)
              do i = 1,norder
                 Is_pre = Is_now
                 Ic_pre = Ic_now

                 Ic_now = (z**i)*sin_RMz_by_RM - (real(i)/RM)*Is_pre
                 Is_now = -(z**i)*cos_RMz_by_RM + (real(i)/RM)*Ic_pre
              enddo
      endif
      reval = Ic_now 
      imval = Is_now 

      return

      end

      !---------------------------------------------------------
      !SUBROUTINE:   TS-expansion of p(x) at x = x0
      !
      !              This is specific to only
      !              p(x) = x^beta, kind of fn.
      !                        -- wr, 15 Oct, 2010

        subroutine texpand(x,x0,beta,tn,fx)

        !> inputs:
        ! x   -> argument at which the function is to be evaluated
        ! x0  -> TS expansion around this value
        ! beta-> value of the exponent
        ! tn  -> number of terms in TS-expansion 
        ! < output:
        ! fx  -> evaluated function

        implicit none

        real*4 x, x0, beta,fx
        integer*4 i, inow 
        integer*4 tn, prod, fact

        real*4 accum
        real*4 ddx_fx0

        accum = 0.0
        prod = 1
        do i = 1,tn
           inow = i
           call ddxn_xpowalpha(x0,beta,inow,ddx_fx0)
           fact = prod*i

           accum = accum + real((ddx_fx0*(x - x0)**i)/fact)
        enddo
        fx = accum + x0**beta

        end


      !---------------------------------------------------------
      !SUBROUTINE:   Derivative of p(x) at x = x0
      !
      !              This is specific to only
      !              p(x) = x^beta, kind of fn.
      !                        -- wr, 15 Oct, 2010

        subroutine ddxn_xpowalpha(x0,alpha,norder,ddxn)

        implicit none

        integer*4 norder
        real*4 x0, alpha, ddxn
        integer*4 i, k

        k = norder
        if(norder.lt.0)then
                write(*,*)"Invalid norder: ",norder
                write(*,*)"Stopping here..."
                stop
        endif
        if(norder.gt.0)then
                ddxn = 1.0d0
                do i = 1,norder
                   ddxn = ddxn*dble(alpha - k + 1)
                   k = k - 1
                enddo
                ddxn = ddxn*(x0**(alpha - norder))
        else
                ddxn = x0**alpha
        endif

        return

        end
