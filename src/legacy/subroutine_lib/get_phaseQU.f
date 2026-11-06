chelp+
      !-------------------------------------------
      ! This code finds the phase of the cosine 
      ! or sine component of a function for which 
      ! the cosine and sine functions can be taken 
      ! as orthogonal basis functions. 
      !                --wr, 07 Dec, 2010
      !-------------------------------------------
chelp-

      subroutine get_phaseQ(omega,T,Y,N,phase_val)

      implicit none
      real*4        omega, T(*), Y(*), phase_val  
      integer*4     N 
      real*4        pi
      parameter     (pi = 3.14159265358979)

      !---------------------------------------
      ! Some math: 
      ! Q = cos(wt + phi)
      !   = cos(wt)cos(phi) - sin(wt)sin(phi)
      ! 
      ! So,
      !     <Q.cos(wt)> = cos(phi) = rc_cor
      !     <Q.sin(-wt)> = sin(phi) = rs_cor
      !
      ! hence:
      !       phi = atan2(rs_cor,rc_cor)
      !---------------------------------------


      do i = 1,N
         phi = omega*T(i)
         c_template(i) = cos(phi)
         s_template(i) = -sin(phi)
      enddo
      call dotproduct(Y,c_template,rc_cor,N)
      call dotproduct(Y,s_template,rs_cor,N)

      phase_val = atan2(rs_cor,rc_cor)

      return
      end


      subroutine get_phaseU(omega,T,Y,N,phase_val)

      implicit none
      real*4        omega, T(*), Y(*), phase_val  
      integer*4     N 
      real*4        pi
      parameter     (pi = 3.14159265358979)


      !---------------------------------------
      ! Some math: 
      ! U = sin(wt + phi)
      !   = sin(wt)cos(phi) + cos(wt)sin(phi)
      ! 
      ! So,
      !     <U.cos(wt)> = sin(phi) = rc_cor
      !     <U.sin(wt)> = cos(phi) = rs_cor
      !
      ! hence:
      !       phi = atan2(rc_cor,rs_cor)
      !---------------------------------------

      do i = 1,N
         phi = omega*T(i)
         c_template(i) = cos(phi)
         s_template(i) = sin(phi)
      enddo
      call dotproduct(Y,c_template,rc_cor,N)
      call dotproduct(Y,s_template,rs_cor,N)

      phase_val = atan2(rc_cor,rs_cor)

      return
      end
