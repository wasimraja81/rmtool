Chelp+
      ! ---------------------------------------------
      ! This subroutine computes the dirty rm-beams 
      ! as a function of RM. 
      ! It is assumed that the signal is of the form:
      !         y(t) = A*exp[i*omega*t]
      !                -- wasim, 11 Aug, 2009
      ! ---------------------------------------------
      ! Last modification:  
      !                 -> Incorporated scheme to 
      !                    account for the bandwidth 
      ! depolarization during cleaning. 
      !                    --wr, 09 Sep, 2011
      !----------------------------------------------



      subroutine compute_dirty_rmbeam(L_sq,nchan,RM_in,phase_in, 
     -                        RM_samp, nrm, re_beam, im_beam,
     -                        cos_arr,sin_arr,maxrm,maxchan,
     -                        bw_correct,tn,z0,z1,z2)
      
      ! NOTE: 
      !      -> L_sq and z0 both contain the sampled lambda-squared 
      !         values. However their sort-order is opposite. 
      implicit none

      !include '../INCLUDE/extract_rm.inc'

      integer*4   nchan, nrm
      integer*4   tn 
      integer*4   maxrm, maxchan 
      real*4      L_sq(*), RM_in, phase_in, RM_samp(*) 
      real*4      re_beam(maxrm), im_beam(maxrm) 
      real*4      ryt(nchan), iyt(nchan)

      real*4      rc_cor, ic_cor, rs_cor, is_cor
      real*4      ryw_tmp, iyw_tmp, phi_tmp 
         
      real*4      c_template(nchan), s_template(nchan)
      real*4      fac
      real*4      cos_arr(maxrm,maxchan),sin_arr(maxrm,maxchan)
      logical     bw_correct  

      real*8      z0(*), z1(*), z2(*) 
      real*8      wt_amp(nchan), wt_pha(nchan) 
      !real*4      lsq_ref 

C COUNTERS:
      integer*4   j,kk

C CONSTANTS:
      !fac = acos(-1.0d0)
      fac = 3.14159265358979

      !call mean(L_sq,nchan,lsq_ref)
      !lsq_ref = 0.0 

      do kk = 1,nchan
         wt_amp(kk) = 1.0d0
         wt_pha(kk) = 0.0d0 
      enddo

      if(bw_correct)then
              call bw_depol_correct(z0, z1, z2, nchan, RM_in, tn, 
     -                              wt_amp, wt_pha)

      endif

      do kk = 1,nchan ! number of channels
         !phi_tmp = 2.0*(RM_in*L_sq(kk) + phase_in) + real(wt_pha(kk)) 
         phi_tmp = 2.0*(RM_in*L_sq(kk) + phase_in) + 
     -                                           real(wt_pha(kk)) 
         ! Generate the inputs whose response 
         ! you wish: 
         ryt(kk) = real(wt_amp(kk))*cos(phi_tmp)
         iyt(kk) = real(wt_amp(kk))*sin(phi_tmp)
      enddo
      do j = 1,nrm
         do kk = 1,nchan ! number of channels
            !phi_tmp = RM_samp(j)*L_sq(kk)
            c_template(kk) = cos_arr(j,kk) !cos(phi_tmp)
            s_template(kk) = sin_arr(j,kk) ! -sin(phi_tmp)
         enddo

         call dotproduct(ryt,c_template,rc_cor,nchan)
         call dotproduct(ryt,s_template,rs_cor,nchan)
         call dotproduct(iyt,c_template,ic_cor,nchan)
         call dotproduct(iyt,s_template,is_cor,nchan)
         rc_cor = rc_cor/dble(nchan)
         rs_cor = rs_cor/dble(nchan)
         ic_cor = ic_cor/dble(nchan)
         is_cor = is_cor/dble(nchan)
         ! Combine coherently to construct y(w)
         ryw_tmp = rc_cor - is_cor  ! Real-part of beam
         iyw_tmp = rs_cor + ic_cor  ! Imag-part of beam

         re_beam(j) = ryw_tmp
         im_beam(j) = iyw_tmp
      enddo
      return
      
      end
