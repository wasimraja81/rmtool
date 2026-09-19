C --------------------------------------------------------------
! This subroutine computes the dirty rm-beams as a function of 
! RM. 
! It is assumed that the signal is of the form:
!             y(t) = A*exp[i*omega*t]

!  -- wasim, 11 Aug, 2009


C --------------------------------------------------------------



          !--------------------------------------------------------------------

          subroutine compute_dirty_rmbeamQ(L_sq,W,nchan,RM_in,phase_in, 
     -                        RM_samp, nrm, Q_beam,maxrm)
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          integer*4   nchan, nrm
          integer*4   maxrm
          real*4      L_sq(*), RM_in, phase_in, W(*), RM_samp(*) 
          real*4      Q_beam(maxrm) 
          real*4      ryt(nchan) 

          real*4      rc_cor 
          real*4      ryw_tmp, phi_tmp 
             
          real*4      c_template(nchan) 
          real*4      pi 

C COUNTERS:
          integer*4   j,kk

C CONSTANTS:
          pi = 3.14159265358979

          do kk = 1,nchan ! number of channels
             phi_tmp = RM_in*L_sq(kk) + phase_in  
             ! Generate the inputs whose response 
             ! you wish: 
             ryt(kk) = W(kk)*cos(phi_tmp)
          enddo
          do j = 1,nrm
             do kk = 1,nchan ! number of channels
                phi_tmp = RM_samp(j)*L_sq(kk)
                c_template(kk) = cos(phi_tmp)
             enddo

             call dotproduct(ryt,c_template,rc_cor,nchan)
             rc_cor = rc_cor/real(nchan)
             ! Combine coherently to construct y(w)
             Q_beam(j) = rc_cor
          enddo
          return
          
          end


          subroutine compute_dirty_rmbeamU(L_sq,W,nchan,RM_in,phase_in, 
     -                        RM_samp, nrm, re_beam, im_beam,maxrm)
          
          implicit none

          !include '../INCLUDE/extract_rm.inc'

          integer*4   nchan, nrm
          integer*4   maxrm
          real*4      L_sq(*), RM_in, phase_in, W(*), RM_samp(*) 
          real*4      U_beam(maxrm) 
          real*4      ryt(nchan) 

          real*4      rs_cor
          real*4      ryw_tmp, phi_tmp 
             
          real*4      s_template(nchan)
          real*4      pi 

C COUNTERS:
          integer*4   j,kk

C CONSTANTS:
          pi = 3.14159265358979

          do kk = 1,nchan ! number of channels
             phi_tmp = RM_in*L_sq(kk) + phase_in  
             ! Generate the inputs whose response 
             ! you wish: 
             ryt(kk) = W(kk)*sin(phi_tmp)
          enddo
          do j = 1,nrm
             do kk = 1,nchan ! number of channels
                phi_tmp = RM_samp(j)*L_sq(kk)
                s_template(kk) = sin(phi_tmp)
             enddo

             call dotproduct(ryt,s_template,rs_cor,nchan)
             rs_cor = rs_cor/real(nchan)
             ! Combine coherently to construct y(w)
             U_beam(j) = rs_cor
          enddo
          return

          end
          
