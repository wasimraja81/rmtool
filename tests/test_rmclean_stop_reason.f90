program test_rmclean_stop_reason
   !! Validates clean_complex's THREE independent stopping criteria
   !! (planning/RMCLEAN_INTEGRATION_PLAN.md T8 -- replaced the old single
   !! overloaded `thresh` parameter, which silently treated an absolute
   !! flux value as a sigma-multiplier, comparing flux^2 against flux):
   !! niter cap, abs_flux_floor, and auto_nsigma (per-pixel full-spectrum
   !! IQR-based sigma, estimated once, not recomputed per iteration). n_iter_used
   !! alone cannot tell "hit the niter cap" apart from "a criterion
   !! happened to fire on the literal last allowed iteration" --
   !! stop_reason ('niter'/'abs_flux'/'auto_nsigma') disambiguates.
   !!
   !! Sky model: single Faraday-thin point source (RM=50 rad/m^2,
   !! amp=10, chi0=0.3 rad -- same convention as tests/test_drm_floor.f90),
   !! resolution-adequate grid (well above Gate 0's floor).
   !!
   !! Case A (niter-only, both optional criteria off): CLEAN must run
   !! every one of niter=50 iterations -- stop_reason=='niter',
   !! n_iter_used==niter exactly, no early exit possible by construction.
   !! Case B (abs_flux_floor, noiseless): stop_reason=='abs_flux',
   !! n_iter_used<niter, and the final residual's own peak amplitude is
   !! at/below the floor (the exact criterion being tested).
   !! Case C (auto_nsigma, WITH injected per-channel noise -- a
   !! noiseless spectrum's own IQR would still be ~0 (RMSF-sidelobe
   !! structure, not statistical noise), making auto_nsigma fire only
   !! at machine precision): stop_reason=='auto_nsigma', n_iter_used<niter.
   !!
   !! Also exercises the optional trace_peak_val/trace_rms_val/
   !! trace_flux_val outputs on Case B: confirms they're filled for
   !! iter=1..n_iter_used, that the peak trend is (weakly) monotonically
   !! decreasing and the flux trend weakly increasing (starting at
   !! exactly 0) -- the basic signature of a converging CLEAN -- and
   !! confirms the ORIGINAL call form (trace arguments omitted entirely,
   !! Case A/C) still works, since most call sites in rmclean_cubes.f90
   !! rely on that.
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp
   real(sp), parameter :: point_rm = 50.0_sp, point_amp = 10.0_sp, chi0_true = 0.3_sp
   integer, parameter :: nchan = 61
   integer, parameter :: nrm = 201
   real(sp), parameter :: rm_lo = -100.0_sp, rm_hi = 100.0_sp

   real(sp) :: l_sq(nchan), q(nchan), u(nchan), q_noisy(nchan), u_noisy(nchan)
   real(sp) :: lsq_ref
   real(sp) :: rm_samp(nrm), dirty_re(nrm), dirty_im(nrm)
   real(sp) :: dirty_re_noisy(nrm), dirty_im_noisy(nrm), resid_amp(nrm)
   real(sp) :: comp_re(nrm), comp_im(nrm), resid_re(nrm), resid_im(nrm)
   real(sp) :: comp_rm_refined(nrm)
   real(sp) :: trace_peak(500), trace_rms(500), trace_flux(500)
   type(rmsf_table_t) :: table
   integer :: n_iter_used, j
   character(len=16) :: stop_reason
   logical :: all_pass

   all_pass = .true.

   call band_channels(300.0_sp, 30.0_sp, nchan, l_sq)
   lsq_ref = sum(l_sq)/real(nchan, sp)
   call inject(l_sq, nchan, q, u)
   do j = 1, nrm
      rm_samp(j) = rm_lo + real(j-1, sp)*(rm_hi-rm_lo)/real(nrm-1, sp)
   end do
   call dirty_spectrum(l_sq, nchan, lsq_ref, q, u, rm_samp, nrm, dirty_re, dirty_im)
   call build_rmsf_offset_table(l_sq, nchan, lsq_ref, rm_samp(nrm)-rm_samp(1),&
   &rm_samp(2)-rm_samp(1), 20, table)

   ! --- Case A: niter-only, must run every iteration --------------------
   call clean_complex(l_sq, nchan, lsq_ref, rm_samp, nrm, dirty_re, dirty_im,&
   &table, 50, 0.1_sp, .false., 0.0_sp, .false., 0.0_sp,&
   &comp_re, comp_im, resid_re, resid_im, n_iter_used, stop_reason,&
   &comp_rm_refined)

   call check(trim(stop_reason) == 'niter',&
   &'Case A (niter-only, no abs_flux_floor/auto_nsigma):'&
   &//' stop_reason==niter', all_pass)
   call check(n_iter_used == 50,&
   &'Case A: n_iter_used==niter exactly (no early exit is possible'&
   &//' with both optional criteria off)', all_pass)

   ! --- Case B: abs_flux_floor, noiseless -------------------------------
   call clean_complex(l_sq, nchan, lsq_ref, rm_samp, nrm, dirty_re, dirty_im,&
   &table, 500, 0.1_sp, .true., 0.5_sp, .false., 0.0_sp,&
   &comp_re, comp_im, resid_re, resid_im, n_iter_used, stop_reason,&
   &comp_rm_refined, trace_peak_val=trace_peak, trace_rms_val=trace_rms,&
   &trace_flux_val=trace_flux)

   call check(trim(stop_reason) == 'abs_flux',&
   &'Case B (abs_flux_floor=0.5, noiseless): stop_reason==abs_flux',&
   &all_pass)
   call check(n_iter_used < 500,&
   &'Case B: stopped before exhausting niter (n_iter_used<500)', all_pass)
   resid_amp = sqrt(resid_re**2 + resid_im**2)
   call check(maxval(resid_amp) <= 0.5_sp + 0.05_sp,&
   &'Case B: final residual peak is at/near the abs_flux_floor (the'&
   &//' exact criterion just tested)', all_pass)
   call check(trace_peak(1) > 5.0_sp,&
   &'Case B: iteration-1 trace peak_val reflects the real source (>5)',&
   &all_pass)
   call check(is_weakly_decreasing(trace_peak, n_iter_used),&
   &'Case B: trace_peak_val is weakly monotonically decreasing (real'&
   &//' convergence trend)', all_pass)
   call check(trace_flux(1) == 0.0_sp,&
   &'Case B: trace_flux_val(1)==0 exactly (nothing cleaned before the'&
   &//' first iteration''s own subtraction)', all_pass)
   call check(is_weakly_increasing(trace_flux, n_iter_used),&
   &'Case B: trace_flux_val (cumulative cleaned flux) is weakly'&
   &//' monotonically increasing', all_pass)

   ! --- Case C: auto_nsigma, WITH injected per-channel noise ------------
   ! A noiseless spectrum's own IQR is dominated by near-zero/RMSF-
   ! sidelobe structure rather than genuine statistical noise, which
   ! would make auto_nsigma fire only once CLEAN has driven the residual
   ! to machine precision -- not a meaningful test of the criterion
   ! itself. Fixed-seed Gaussian-ish noise (Box-Muller from
   ! random_number) gives a genuine, reproducible noise sigma.
   call inject_noise(nchan, 0.3_sp, q, u, q_noisy, u_noisy)
   call dirty_spectrum(l_sq, nchan, lsq_ref, q_noisy, u_noisy, rm_samp, nrm,&
   &dirty_re_noisy, dirty_im_noisy)
   call clean_complex(l_sq, nchan, lsq_ref, rm_samp, nrm, dirty_re_noisy,&
   &dirty_im_noisy, table, 500, 0.1_sp, .false., 0.0_sp, .true.,&
   &5.0_sp, comp_re, comp_im, resid_re, resid_im, n_iter_used,&
   &stop_reason, comp_rm_refined)

   call check(trim(stop_reason) == 'auto_nsigma',&
   &'Case C (auto_nsigma=5.0, WITH injected noise): stop_reason=='&
   &//'auto_nsigma', all_pass)
   call check(n_iter_used < 500,&
   &'Case C: stopped before exhausting niter (n_iter_used<500)', all_pass)

   call destroy_rmsf_offset_table(table)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_rmclean_stop_reason: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_rmclean_stop_reason: one or more checks failed'
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

   subroutine inject_noise(n, sigma, q_in, u_in, q_out, u_out)
      !! Adds fixed-seed Gaussian-ish (Box-Muller) noise of the given
      !! per-channel sigma onto q_in/u_in -- reproducible run to run.
      integer, intent(in) :: n
      real(sp), intent(in) :: sigma, q_in(n), u_in(n)
      real(sp), intent(out) :: q_out(n), u_out(n)
      integer :: seed_size, k
      integer, allocatable :: seed_arr(:)
      real(sp) :: r1, r2

      call random_seed(size=seed_size)
      allocate(seed_arr(seed_size))
      seed_arr = 20250101
      call random_seed(put=seed_arr)
      do k = 1, n
         call random_number(r1)
         call random_number(r2)
         r1 = max(r1, 1.0e-7_sp)
         q_out(k) = q_in(k) + sigma*sqrt(-2.0_sp*log(r1))*cos(6.28318530718_sp*r2)
         call random_number(r1)
         call random_number(r2)
         r1 = max(r1, 1.0e-7_sp)
         u_out(k) = u_in(k) + sigma*sqrt(-2.0_sp*log(r1))*cos(6.28318530718_sp*r2)
      end do
      deallocate(seed_arr)
   end subroutine inject_noise

   subroutine dirty_spectrum(lsq_in, n, lsq_ref_in, q_in, u_in, rm_in, nrm_in,&
   &dre_out, dim_out)
      integer, intent(in) :: n, nrm_in
      real(sp), intent(in) :: lsq_in(n), lsq_ref_in, q_in(n), u_in(n), rm_in(nrm_in)
      real(sp), intent(out) :: dre_out(nrm_in), dim_out(nrm_in)
      real(sp) :: c_tmpl(n), s_tmpl(n)
      integer :: j
      do j = 1, nrm_in
         c_tmpl = cos(2.0_sp*rm_in(j)*(lsq_in-lsq_ref_in))
         s_tmpl = -sin(2.0_sp*rm_in(j)*(lsq_in-lsq_ref_in))
         dre_out(j) = (dot_product(q_in, c_tmpl) - dot_product(u_in, s_tmpl))/real(n, sp)
         dim_out(j) = (dot_product(q_in, s_tmpl) + dot_product(u_in, c_tmpl))/real(n, sp)
      end do
   end subroutine dirty_spectrum

   function is_weakly_decreasing(arr, n) result(ok)
      real(sp), intent(in) :: arr(:)
      integer, intent(in) :: n
      logical :: ok
      integer :: k
      ok = .true.
      do k = 2, n
         if (arr(k) > arr(k-1) + 1.0e-4_sp) then
            ok = .false.
            exit
         endif
      end do
   end function is_weakly_decreasing

   function is_weakly_increasing(arr, n) result(ok)
      real(sp), intent(in) :: arr(:)
      integer, intent(in) :: n
      logical :: ok
      integer :: k
      ok = .true.
      do k = 2, n
         if (arr(k) < arr(k-1) - 1.0e-4_sp) then
            ok = .false.
            exit
         endif
      end do
   end function is_weakly_increasing

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

end program test_rmclean_stop_reason
