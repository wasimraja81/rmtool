chelp+
      !---------------------------------------------------
      ! This code RETURNS the complex "correction factor" 
      ! needed to correct the effect of inherent smearing 
      ! of the linear polarised Stokes intensities due 
      ! to the FINITE channel-width in the spectral axis. 
      !                           --wr, 12 Oct, 2010
      !---------------------------------------------------
      ! Last modification: 
      !     -> The required arrays z0, z1 and z2 are now 
      !        passed as arguments. These should be pre-
      !        computed using the subroutine : 
      !              bw_depol_correct_setup
      !        in the calling program itself. 
      !                                --wr, 12 Sep, 2011
      !---------------------------------------------------
chelp-

      subroutine bw_depol_correct(z0, z1, z2, nchan,
     -                           RM_in, tn, weight_amp, 
     -                           weight_pha)


      implicit none

      !--------------------------------------
      real*8       z0(*), z1(*), z2(*)
      real*4       RM_in 
      integer*4    nchan, tn 
      real*8       weight_amp(*), weight_pha(*) 
      integer*4    i, j   
      real*8       RM 
      real*8       zz0, zz1, zz2, 
     -             beta, 
     -             reaccum, imaccum, ddxn_fx, reval1, reval2, 
     -             imval1, imval2, norm(nchan) 
      integer*4    order_now 
      integer*8    factorial 
      !--------------------------------------
      ! Constants: 
      !real*8       vel_light
      real*8       pi

      !---------------------------------------
      real*8       tmp_num 


      pi = 3.14159265358979d0

      ! Apply the corrections to the smoothed Q and U: 
      ! beta is the index of lambda-squared appearing 
      ! due to the conversion from freq-to-lambda-squared.
      ! 
      beta = -1.5d0 ! Taking ONLY bw-depol correction 
                    ! related index; To add the bandshape 
                    ! index, and/or the spectral index 
                    ! add those indices as well appropriately.


      !---------------------------------------------------------
      do i = 1,nchan
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
      RM = 2.0d0*RM_in ! DOUBT: Fourier pairs are:
                       !        Lsq <-> RM or 
                       !        Lsq <-> 2RM ?
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
     -                       ddxn_fx*(reval2-reval1)/dble(factorial) 
            imaccum = imaccum + 
     -                       ddxn_fx*(imval2-imval1)/dble(factorial) 
         enddo
         weight_amp(i) = sqrt(reaccum*reaccum + imaccum*imaccum)/
     -                                                       norm(i)
         weight_pha(i) = atan2(imaccum,reaccum)

         !! Remove the phase introduced by smearing:
         !imaccum = atan2(U_in(i),Q_in(i)) ! phase
         !reaccum = sqrt(Q_in(i)*Q_in(i)+U_in(i)*U_in(i))
         !Q_corr(i) = reaccum*cos(imaccum - ph_weight)/weight
         !U_corr(i) = reaccum*sin(imaccum - ph_weight)/weight
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
      real*8    RM, z,  
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

        real*8 x, x0, beta,fx
        integer*4 i, inow 
        integer*4 tn, prod, fact

        real*8 accum
        real*8 ddx_fx0

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
        real*8 x0, alpha, ddxn
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
