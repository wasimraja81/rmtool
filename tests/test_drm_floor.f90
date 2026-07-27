program test_drm_floor
   !! Demonstrates WHY get_drm's own enforced floor is oversample>=2, not
   !! the bare 2-point Nyquist value of 1 -- the empirical evidence
   !! get_drm's own doc comment cites, made into a real, re-runnable
   !! regression check rather than a frozen number in a comment. Added
   !! after the user's own simpler redesign request: instead of letting a
   !! caller specify dRM directly, expose ONLY an oversampling factor
   !! (always applied on top of the mandatory Nyquist-type floor), so
   !! undersampling is architecturally impossible -- this test confirms
   !! the floor is set in the right place.
   !!
   !! Since get_drm now correctly REFUSES oversample<2 (tests/
   !! test_drm_floor_enforcement.f90 checks that refusal directly), the
   !! unsafe (oversample=1) case here is computed by hand, bypassing
   !! get_drm deliberately, purely to demonstrate what the enforcement is
   !! protecting against -- not a suggested usage pattern.
   !!
   !! Sky model: single Faraday-thin point source, RM=50 rad/m^2,
   !! amplitude=10 Jy/(rad/m^2), intrinsic angle chi0=0.3 rad (nonzero,
   !! same reasoning as tests/test_rmclean_lsqref_flex.f90's own sky
   !! model: a nonzero chi0 catches sign errors that "0 in, 0 out" would
   !! hide). Band-mean compute reference (lsq_ref_compute != 0), so any
   !! chi0 error from under-sampling actually shows up in the reported
   !! angle (see derotate_to_lsq_ref's own comment on why lsq_ref_
   !! compute=0 would hide this).
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp
   real(sp), parameter :: point_rm = 50.0_sp, point_amp = 10.0_sp, chi0_true = 0.3_sp
   real(sp), parameter :: pi_const = 3.14159265_sp
   integer, parameter :: nchan = 61
   real(sp), parameter :: rm_lo = -100.0_sp, rm_span = 300.0_sp

   real(sp) :: l_sq(nchan), q(nchan), u(nchan), lsq_ref_compute
   real(sp) :: chi0_at_oversample_1, chi0_at_oversample_2
   logical :: all_pass

   all_pass = .true.

   call band_channels(300.0_sp, 30.0_sp, nchan, l_sq)
   call get_lsq_ref_compute(l_sq, nchan, mode=lsq_ref_compute_centroid,&
   &lsq_ref_compute=lsq_ref_compute)
   call inject(l_sq, nchan, q, u)

   ! Unsafe: oversample=1 (bare Nyquist), computed BY HAND since get_drm
   ! itself now correctly refuses this -- see this file's own top comment.
   call run_case_raw_oversample(l_sq, nchan, lsq_ref_compute, q, u, 1.0_sp, chi0_at_oversample_1)
   ! Safe: oversample=2, get_drm's own enforced floor, via the real API.
   call run_case(l_sq, nchan, lsq_ref_compute, q, u, 2.0_sp, chi0_at_oversample_2)

   write(*,'(A)') '===================================================='
   write(*,'(A,F0.4,A,F0.4,A)') 'oversample=1 (bare Nyquist, computed by hand): chi0=',&
   &chi0_at_oversample_1, ' rad (true=', chi0_true, ' rad)'
   write(*,'(A,F0.4,A,F0.4,A)') 'oversample=2 (get_drm''s own enforced floor): chi0=',&
   &chi0_at_oversample_2, ' rad (true=', chi0_true, ' rad)'
   write(*,'(A)') '===================================================='

   call check(abs(chi0_at_oversample_1-chi0_true) > 0.1_sp,&
   &'oversample=1 (bare Nyquist) gives a WRONG chi0 (>0.1 rad off) -- this is'&
   &//' exactly what get_drm''s own floor=2 protects against', all_pass)
   call check(abs(chi0_at_oversample_2-chi0_true) <= 0.02_sp,&
   &'oversample=2 (get_drm''s own enforced floor) already recovers chi0 correctly', all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_drm_floor: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_drm_floor: one or more checks failed'
      stop 1
   endif

contains

   subroutine band_channels(nu_c, dnu, n, lsq_out)
      real(sp), intent(in) :: nu_c, dnu
      integer, intent(in) :: n
      real(sp), intent(out) :: lsq_out(n)
      real(sp) :: f_step, f_start, freq
      integer :: k
      f_step = dnu/real(n, sp)
      f_start = nu_c - 0.5_sp*dnu + 0.5_sp*f_step
      do k = 1, n
         freq = f_start + real(k-1, sp)*f_step
         lsq_out(k) = (c_light/freq)**2
      end do
   end subroutine band_channels

   subroutine inject(lsq_in, n, q_out, u_out)
      integer, intent(in) :: n
      real(sp), intent(in) :: lsq_in(n)
      real(sp), intent(out) :: q_out(n), u_out(n)
      real(sp) :: phase
      integer :: k
      do k = 1, n
         phase = 2.0_sp*(point_rm*lsq_in(k) + chi0_true)
         q_out(k) = point_amp*cos(phase)
         u_out(k) = point_amp*sin(phase)
      end do
   end subroutine inject

   subroutine run_case(lsq_in, n, lsq_ref, q_in, u_in, oversample, chi0_found)
      integer, intent(in) :: n
      real(sp), intent(in) :: lsq_in(n), lsq_ref, q_in(n), u_in(n), oversample
      real(sp), intent(out) :: chi0_found
      real(sp) :: drm

      call get_drm(lsq_in, n, lsq_ref, drm, oversample=oversample)
      call run_case_at_drm(lsq_in, n, lsq_ref, q_in, u_in, drm, chi0_found)
   end subroutine run_case

   subroutine run_case_raw_oversample(lsq_in, n, lsq_ref, q_in, u_in, oversample, chi0_found)
      !! Bypasses get_drm deliberately (it correctly refuses oversample<2)
      !! to demonstrate what its enforcement protects against -- NOT a
      !! pattern to copy elsewhere.
      integer, intent(in) :: n
      real(sp), intent(in) :: lsq_in(n), lsq_ref, q_in(n), u_in(n), oversample
      real(sp), intent(out) :: chi0_found
      real(dp), parameter :: pi = 3.14159265358979_dp
      real(dp) :: max_offset
      real(sp) :: drm

      max_offset = maxval(abs(real(lsq_in, dp) - real(lsq_ref, dp)))
      drm = real(pi/(2.0_dp*real(oversample, dp)*max_offset), sp)
      call run_case_at_drm(lsq_in, n, lsq_ref, q_in, u_in, drm, chi0_found)
   end subroutine run_case_raw_oversample

   subroutine run_case_at_drm(lsq_in, n, lsq_ref, q_in, u_in, drm, chi0_found)
      integer, intent(in) :: n
      real(sp), intent(in) :: lsq_in(n), lsq_ref, q_in(n), u_in(n), drm
      real(sp), intent(out) :: chi0_found

      integer :: nrm_l, j, ipeak, n_iter_used
      real(sp), allocatable :: rm(:), dirty_re(:), dirty_im(:)
      real(sp), allocatable :: comp_re(:), comp_im(:), resid_re(:), resid_im(:)
      real(sp), allocatable :: comp_rm_refined(:), comp_amp(:)
      real(sp) :: c_tmpl(n), s_tmpl(n)
      type(rmsf_table_t) :: table
      real(sp) :: rm_found

      nrm_l = nint(rm_span/drm) + 1
      allocate(rm(nrm_l), dirty_re(nrm_l), dirty_im(nrm_l))
      allocate(comp_re(nrm_l), comp_im(nrm_l), resid_re(nrm_l), resid_im(nrm_l))
      allocate(comp_rm_refined(nrm_l), comp_amp(nrm_l))
      do j = 1, nrm_l
         rm(j) = rm_lo + real(j-1, sp)*drm
      end do
      do j = 1, nrm_l
         c_tmpl = cos(2.0_sp*rm(j)*(lsq_in-lsq_ref))
         s_tmpl = -sin(2.0_sp*rm(j)*(lsq_in-lsq_ref))
         dirty_re(j) = (dot_product(q_in, c_tmpl) - dot_product(u_in, s_tmpl))/real(n, sp)
         dirty_im(j) = (dot_product(q_in, s_tmpl) + dot_product(u_in, c_tmpl))/real(n, sp)
      end do
      call build_rmsf_offset_table(lsq_in, n, lsq_ref, rm(nrm_l)-rm(1), drm, 20, table)
      call clean_complex(rm, nrm_l, dirty_re, dirty_im, table, 500, 0.1_sp, 1.0e-4_sp,&
      &comp_re, comp_im, resid_re, resid_im, n_iter_used, comp_rm_refined)
      comp_amp = sqrt(comp_re**2+comp_im**2)
      ipeak = 1
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) < abs(rm(ipeak)-point_rm)) ipeak = j
      end do
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) <= 10.0_sp .and. comp_amp(j) > comp_amp(ipeak)) ipeak = j
      end do
      rm_found = comp_rm_refined(ipeak)
      chi0_found = 0.5_sp*atan2(comp_im(ipeak), comp_re(ipeak)) - rm_found*lsq_ref
      chi0_found = chi0_found - nint(chi0_found/pi_const)*pi_const

      call destroy_rmsf_offset_table(table)
      deallocate(rm, dirty_re, dirty_im, comp_re, comp_im, resid_re, resid_im, comp_rm_refined, comp_amp)
   end subroutine run_case_at_drm

   subroutine check(cond, label, all_pass_io)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: label
      logical, intent(inout) :: all_pass_io
      if (cond) then
         write(*,'(A,A)') '  [PASS] ', label
      else
         write(*,'(A,A)') '  [FAIL] ', label
         all_pass_io = .false.
      endif
   end subroutine check

end program test_drm_floor
