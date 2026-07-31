program test_matched_filter_refine
   !! planning/RMCLEAN_INTEGRATION_PLAN.md tickets T3 + T3c: validates
   !! rmclean_mod's refine_peak_matched_filter -- both the original
   !! full-search path (T3) and the tiered fast-path/escalation design
   !! (T3c, from a design discussion with the user: a cheap
   !! fixed-location fit runs first every call, and the full search runs
   !! only when the fast fit's leftover misfit exceeds
   !! nsigma * noise_rms, a data-driven consistency cut).
   !!
   !! Scenario: identical to tests/test_drm_floor.f90's own P-band setup
   !! (nu_c=300 MHz, dnu=30 MHz, 61 channels), TWO point sources tested
   !! independently (RM=50/chi0=0.3, and RM=-27.3/chi0=-1.1, to rule out
   !! a result specific to one particular RM/phase combination),
   !! lsq_ref=0.0 (this project's historical convention -- the WORST
   !! case for the fast carrier, max_offset ~= the channels' own
   !! absolute lambda^2).
   !!
   !! Per (scenario x sampling density), THREE ways through the new
   !! routine are checked, each against truth:
   !!   FORCED SEARCH (nsigma=-1, threshold impossible to satisfy):
   !!     exercises the full-search tier alone -- T3's own original
   !!     validated numbers must keep holding here.
   !!   FORCED FAST (nsigma=huge): exercises the fast tier alone --
   !!     parabola location + closed-form anchored fit, no search. The
   !!     claim validated: this alone is already accurate in clean
   !!     single-component data, even at coarse resolution-only
   !!     sampling (the whole point of the tiered design's cheapness).
   !!   TIERED (nsigma=3, the production default): whichever path the
   !!     criterion actually picks must land on truth.
   !! Plus two dedicated tier-mechanics checks (same threshold, same
   !! noise level, differing ONLY in the data):
   !!   clean single component + moderate noise floor -> fast path must
   !!     be ACCEPTED (used_search=.false.);
   !!   two blended components + same noise floor -> the single-
   !!     component fast fit CANNOT explain the anchors, so the
   !!     escalation must FIRE (used_search=.true.).
   !! And a noisy single-component case at the default nsigma=3,
   !! confirming the default behaves sanely on data with a real noise
   !! floor (deterministic hand-rolled LCG noise, reproducible across
   !! compilers -- not compiler-seeded random_number).
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp
   real(sp), parameter :: pi_const = 3.14159265_sp
   integer, parameter :: nchan = 61
   real(sp), parameter :: lsq_ref = 0.0_sp

   real(sp) :: l_sq(nchan)
   real(sp) :: fwhm_rm, drm_os2, drm_os4
   logical :: all_pass

   all_pass = .true.

   call band_channels(300.0_sp, 30.0_sp, nchan, l_sq)
   call compute_rmsf_fwhm(l_sq, nchan, fwhm_rm)
   call get_drm(l_sq, nchan, lsq_ref, drm_os2, oversample=2.0_sp)
   call get_drm(l_sq, nchan, lsq_ref, drm_os4, oversample=4.0_sp)

   write(*,'(A)') '===================================================='
   write(*,'(A,F0.6,A)') 'fwhm = ', fwhm_rm, ' rad/m^2'
   write(*,'(A,F0.6,A)') 'get_drm(oversample=2, enforced floor) = ', drm_os2, ' rad/m^2'
   write(*,'(A,F0.6,A)') 'get_drm(oversample=4, default)        = ', drm_os4, ' rad/m^2'
   write(*,'(A)') '===================================================='

   write(*,'(A)') ''
   write(*,'(A)') '### Scenario 1: RM=50.0, chi0=0.3 ###'
   call run_all_cases(50.0_sp, 1.0_sp, 0.3_sp, drm_os2, drm_os4, fwhm_rm, all_pass)

   write(*,'(A)') ''
   write(*,'(A)') '### Scenario 2: RM=-27.3, chi0=-1.1 (rules out a result specific to scenario 1) ###'
   call run_all_cases(-27.3_sp, 1.0_sp, -1.1_sp, drm_os2, drm_os4, fwhm_rm, all_pass)

   write(*,'(A)') ''
   write(*,'(A)') '### Tier mechanics: clean accepts fast path, blended escalates (T3c) ###'
   call run_tier_mechanics(fwhm_rm, all_pass)

   write(*,'(A)') ''
   write(*,'(A)') '### Noisy single component at default nsigma=3 (T3c default check) ###'
   call run_noisy_case(fwhm_rm, all_pass)

   write(*,'(A)') ''
   write(*,'(A)') '===================================================='
   if (all_pass) then
      write(*,'(A)') '[PASS] test_matched_filter_refine: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_matched_filter_refine: one or more checks failed'
      stop 1
   endif

contains

   subroutine band_channels(nu_c, dnu, n, lsq_out)
      !! Identical construction to tests/test_drm_floor.f90's own --
      !! uniform in FREQUENCY (non-uniform in l_sq, ~1.34x spacing
      !! variation across this exact band, confirmed directly).
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

   subroutine build_dirty(q, u, rm, nrm_l, dirty_re, dirty_im)
      !! Direct-DFT dirty spectrum (no FFT), lsq_ref=0.0 -- the module's
      !! own convention throughout this test.
      real(sp), intent(in) :: q(nchan), u(nchan)
      integer, intent(in) :: nrm_l
      real(sp), intent(in) :: rm(nrm_l)
      real(sp), intent(out) :: dirty_re(nrm_l), dirty_im(nrm_l)
      real(sp) :: c_tmpl(nchan), s_tmpl(nchan)
      integer :: j
      do j = 1, nrm_l
         c_tmpl = cos(2.0_sp*rm(j)*(l_sq-lsq_ref))
         s_tmpl = -sin(2.0_sp*rm(j)*(l_sq-lsq_ref))
         dirty_re(j) = (dot_product(q, c_tmpl) - dot_product(u, s_tmpl))/real(nchan, sp)
         dirty_im(j) = (dot_product(q, s_tmpl) + dot_product(u, c_tmpl))/real(nchan, sp)
      end do
   end subroutine build_dirty

   subroutine run_all_cases(point_rm, point_amp, chi0_true, drm_os2_in,&
   &drm_os4_in, fwhm_rm_in, all_pass_io)
      real(sp), intent(in) :: point_rm, point_amp, chi0_true
      real(sp), intent(in) :: drm_os2_in, drm_os4_in, fwhm_rm_in
      logical, intent(inout) :: all_pass_io

      write(*,'(A)') '--- dRM = get_drm(oversample=4) [T3 regression floor] ---'
      call run_one_case(point_rm, point_amp, chi0_true, drm_os4_in, all_pass_io)
      write(*,'(A)') '--- dRM = get_drm(oversample=2) [T3 regression floor] ---'
      call run_one_case(point_rm, point_amp, chi0_true, drm_os2_in, all_pass_io)
      write(*,'(A)') '--- dRM = fwhm/4 [resolution-only -- old get_drm gate would REFUSE] ---'
      call run_one_case(point_rm, point_amp, chi0_true, fwhm_rm_in/4.0_sp, all_pass_io)
      write(*,'(A)') '--- dRM = fwhm/3 [3 samples/fwhm, coarsest resolution sampling still'//&
      &' accurate -- see run_tier_mechanics''s own comment for why 2 samples/fwhm is'//&
      &' excluded] ---'
      call run_one_case(point_rm, point_amp, chi0_true, fwhm_rm_in/3.0_sp, all_pass_io)
   end subroutine run_all_cases

   subroutine run_one_case(point_rm, point_amp, chi0_true, drm, all_pass_io)
      real(sp), intent(in) :: point_rm, point_amp, chi0_true, drm
      logical, intent(inout) :: all_pass_io

      real(sp), parameter :: rm_lo = -100.0_sp, rm_span = 300.0_sp
      integer :: nrm_l, j, imax
      real(sp), allocatable :: rm(:), dirty_re(:), dirty_im(:)
      real(sp) :: q(nchan), u(nchan), phase(nchan)
      real(sp) :: avg_abs_dummy, noise_rms
      real(sp) :: rm_new, chi0_new
      real(sp) :: peak_loc_new, re_new, im_new
      logical :: used_search

      nrm_l = nint(rm_span/drm) + 1
      allocate(rm(nrm_l), dirty_re(nrm_l), dirty_im(nrm_l))
      do j = 1, nrm_l
         rm(j) = rm_lo + real(j-1, sp)*drm
      end do

      phase = 2.0_sp*(point_rm*l_sq + chi0_true)
      q = point_amp*cos(phase)
      u = point_amp*sin(phase)
      call build_dirty(q, u, rm, nrm_l, dirty_re, dirty_im)

      call index_absmax(sqrt(dirty_re**2+dirty_im**2), nrm_l, imax, avg_abs_dummy)
      call rms_of_amp(dirty_re, dirty_im, nrm_l, noise_rms)

      write(*,'(A,F0.6,A,I0,A)') '  dRM=', drm, '  (nrm=', nrm_l, ')'

      ! FORCED SEARCH (nsigma=-1: threshold impossible, always escalate)
      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_rms, -1.0_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      rm_new = peak_loc_new
      chi0_new = wrap_chi0(0.5_sp*atan2(im_new, re_new) - rm_new*lsq_ref)
      write(*,'(A,F0.4,A,F0.4,A,L1)') '    SEARCH: RM=', rm_new, '  chi0=',&
      &chi0_new, ' rad  used_search=', used_search
      call check(used_search, '    [SEARCH] escalation forced (nsigma=-1)', all_pass_io)
      call check(abs(rm_new-point_rm) <= 0.5_sp*drm,&
      &'    [SEARCH] RM recovered within half a grid cell', all_pass_io)
      call check(abs(chi0_new-chi0_true) <= 0.05_sp,&
      &'    [SEARCH] chi0 recovered within 0.05 rad', all_pass_io)

      ! FORCED FAST (nsigma=huge: threshold always satisfied)
      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_rms, 1.0e30_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      rm_new = peak_loc_new
      chi0_new = wrap_chi0(0.5_sp*atan2(im_new, re_new) - rm_new*lsq_ref)
      write(*,'(A,F0.4,A,F0.4,A,L1)') '    FAST:   RM=', rm_new, '  chi0=',&
      &chi0_new, ' rad  used_search=', used_search
      call check(.not. used_search, '    [FAST] fast path forced (nsigma=huge)', all_pass_io)
      call check(abs(rm_new-point_rm) <= 0.5_sp*drm,&
      &'    [FAST] RM recovered within half a grid cell', all_pass_io)
      call check(abs(chi0_new-chi0_true) <= 0.05_sp,&
      &'    [FAST] chi0 recovered within 0.05 rad', all_pass_io)

      ! TIERED at the production default (nsigma=3)
      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_rms, 3.0_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      rm_new = peak_loc_new
      chi0_new = wrap_chi0(0.5_sp*atan2(im_new, re_new) - rm_new*lsq_ref)
      write(*,'(A,F0.4,A,F0.4,A,L1)') '    TIERED: RM=', rm_new, '  chi0=',&
      &chi0_new, ' rad  used_search=', used_search
      call check(abs(chi0_new-chi0_true) <= 0.05_sp,&
      &'    [TIERED nsigma=3] chi0 recovered within 0.05 rad', all_pass_io)

      write(*,'(A,F0.4,A,F0.4,A)') '    TRUE:   RM=', point_rm, '  chi0=',&
      &chi0_true, ' rad'

      deallocate(rm, dirty_re, dirty_im)
   end subroutine run_one_case

   subroutine run_tier_mechanics(fwhm_rm_in, all_pass_io)
      !! Same threshold, same claimed noise floor -- ONLY the data
      !! differs: a clean single component must be explained by the fast
      !! fixed-location single-component fit (accepted, no search); two
      !! blended components CANNOT be (the leftover misfit is real model
      !! error, far above the claimed noise), so the escalation must
      !! fire. Coarse resolution-only grid (fwhm/3), the tier's own
      !! hardest working regime.
      !!
      !! Why fwhm/3, not fwhm/2 (bug-fix follow-up, found once
      !! compute_rmsf_fwhm's own erroneous extra 0.5 factor was removed
      !! -- see that subroutine's own comment): at the TRUE fwhm/2 (2
      !! samples/fwhm), even a perfectly clean, noiseless single
      !! component has a fast-path fit residual of ~0.008 (pure parabola-
      !! location interpolation error at that coarse a grid, confirmed
      !! directly by instrumenting refine_peak_matched_filter's own
      !! leftover/threshold values) -- 2.7x this test's claimed
      !! noise_floor=0.001 threshold, so escalation fires even with no
      !! noise or model error at all, and the "clean single component:
      !! fast path ACCEPTED" assertion cannot hold there. At fwhm/3 the
      !! same residual drops to ~0.002, safely under threshold, while
      !! the blended-components case still escalates unambiguously
      !! (residual ~0.28, two orders of magnitude above threshold either
      !! way) -- fwhm/2 was only ever passing before this fix because
      !! compute_rmsf_fwhm's own bug made "fwhm/2" actually mean
      !! ~4 samples/fwhm, not the 2 its own label claimed.
      real(sp), intent(in) :: fwhm_rm_in
      logical, intent(inout) :: all_pass_io

      real(sp), parameter :: rm_lo = -100.0_sp, rm_span = 300.0_sp
      real(sp), parameter :: noise_floor = 0.001_sp   ! claimed noise, in dirty-map units
      integer :: nrm_l, j, imax
      real(sp), allocatable :: rm(:), dirty_re(:), dirty_im(:)
      real(sp) :: q(nchan), u(nchan), phase(nchan)
      real(sp) :: drm, avg_abs_dummy
      real(sp) :: peak_loc_new, re_new, im_new
      logical :: used_search

      drm = fwhm_rm_in/3.0_sp
      nrm_l = nint(rm_span/drm) + 1
      allocate(rm(nrm_l), dirty_re(nrm_l), dirty_im(nrm_l))
      do j = 1, nrm_l
         rm(j) = rm_lo + real(j-1, sp)*drm
      end do

      ! Clean single component
      phase = 2.0_sp*(50.0_sp*l_sq + 0.3_sp)
      q = cos(phase)
      u = sin(phase)
      call build_dirty(q, u, rm, nrm_l, dirty_re, dirty_im)
      call index_absmax(sqrt(dirty_re**2+dirty_im**2), nrm_l, imax, avg_abs_dummy)
      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_floor, 3.0_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      write(*,'(A,L1)') '  clean single component:  used_search=', used_search
      call check(.not. used_search,&
      &'  clean single component: fast path ACCEPTED (no escalation)', all_pass_io)

      ! Two blended components, ~0.8 fwhm apart, comparable amplitudes --
      ! a single-component model cannot explain the anchors here.
      phase = 2.0_sp*(50.0_sp*l_sq + 0.3_sp)
      q = cos(phase)
      u = sin(phase)
      phase = 2.0_sp*((50.0_sp+0.8_sp*fwhm_rm_in)*l_sq - 0.5_sp)
      q = q + 0.8_sp*cos(phase)
      u = u + 0.8_sp*sin(phase)
      call build_dirty(q, u, rm, nrm_l, dirty_re, dirty_im)
      call index_absmax(sqrt(dirty_re**2+dirty_im**2), nrm_l, imax, avg_abs_dummy)
      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_floor, 3.0_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      write(*,'(A,L1)') '  blended two components:  used_search=', used_search
      call check(used_search,&
      &'  blended two components: escalation FIRED (search used)', all_pass_io)

      deallocate(rm, dirty_re, dirty_im)
   end subroutine run_tier_mechanics

   subroutine run_noisy_case(fwhm_rm_in, all_pass_io)
      !! Single component + real per-channel noise (deterministic LCG,
      !! reproducible), default nsigma=3, data-driven noise estimate
      !! computed the same way clean_complex's own loop does
      !! (rms_about_mean of the residual amplitude). Confirms the
      !! default threshold neither breaks recovery nor forces the
      !! expensive path on ordinary noisy single-component data.
      real(sp), intent(in) :: fwhm_rm_in
      logical, intent(inout) :: all_pass_io

      real(sp), parameter :: rm_lo = -100.0_sp, rm_span = 300.0_sp
      real(sp), parameter :: point_rm = 50.0_sp, chi0_true = 0.3_sp
      real(sp), parameter :: noise_sigma = 0.02_sp
      integer :: nrm_l, j, imax, k
      real(sp), allocatable :: rm(:), dirty_re(:), dirty_im(:)
      real(sp) :: q(nchan), u(nchan), phase(nchan)
      real(sp) :: drm, avg_abs_dummy, noise_rms
      real(sp) :: rm_new, chi0_new
      real(sp) :: peak_loc_new, re_new, im_new
      logical :: used_search
      integer(dp) :: lcg

      drm = fwhm_rm_in/4.0_sp
      nrm_l = nint(rm_span/drm) + 1
      allocate(rm(nrm_l), dirty_re(nrm_l), dirty_im(nrm_l))
      do j = 1, nrm_l
         rm(j) = rm_lo + real(j-1, sp)*drm
      end do

      phase = 2.0_sp*(point_rm*l_sq + chi0_true)
      lcg = 20260728_dp
      do k = 1, nchan
         q(k) = cos(phase(k)) + noise_sigma*lcg_gauss(lcg)
         u(k) = sin(phase(k)) + noise_sigma*lcg_gauss(lcg)
      end do
      call build_dirty(q, u, rm, nrm_l, dirty_re, dirty_im)
      call index_absmax(sqrt(dirty_re**2+dirty_im**2), nrm_l, imax, avg_abs_dummy)
      call rms_of_amp(dirty_re, dirty_im, nrm_l, noise_rms)

      call refine_peak_matched_filter(l_sq, nchan, lsq_ref, rm, nrm_l,&
      &dirty_re, dirty_im, imax, drm, noise_rms, 3.0_sp,&
      &peak_loc_new, re_new, im_new, used_search)
      rm_new = peak_loc_new
      chi0_new = wrap_chi0(0.5_sp*atan2(im_new, re_new) - rm_new*lsq_ref)
      write(*,'(A,F0.4,A,F0.4,A,L1)') '  noisy (sigma=2% of amp): RM=', rm_new,&
      &'  chi0=', chi0_new, ' rad  used_search=', used_search
      write(*,'(A,F0.4,A,F0.4,A)') '  TRUE:                    RM=', point_rm,&
      &'  chi0=', chi0_true, ' rad'
      call check(abs(rm_new-point_rm) <= 0.5_sp*drm,&
      &'  noisy: RM recovered within half a grid cell at default nsigma=3', all_pass_io)
      call check(abs(chi0_new-chi0_true) <= 0.1_sp,&
      &'  noisy: chi0 recovered within 0.1 rad at default nsigma=3', all_pass_io)

      deallocate(rm, dirty_re, dirty_im)
   end subroutine run_noisy_case

   function lcg_gauss(state) result(g)
      !! Deterministic ~N(0,1) deviate: sum of 12 LCG uniforms minus 6
      !! (Irwin-Hall approximation) -- reproducible across compilers,
      !! unlike random_number's compiler-specific stream.
      integer(dp), intent(inout) :: state
      real(sp) :: g
      integer :: i
      real(dp) :: s
      s = 0.0_dp
      do i = 1, 12
         state = mod(state*6364136223846793005_dp + 1442695040888963407_dp,&
         &9223372036854775807_dp)
         s = s + real(abs(state), dp)/9223372036854775807.0_dp
      end do
      g = real(s - 6.0_dp, sp)
   end function lcg_gauss

   subroutine rms_of_amp(re_in, im_in, n, rms_out)
      !! Same data-driven noise proxy clean_complex's own loop uses
      !! (rms_about_mean of the amplitude) -- duplicated because
      !! rms_about_mean is private to rmclean_mod.
      integer, intent(in) :: n
      real(sp), intent(in) :: re_in(n), im_in(n)
      real(sp), intent(out) :: rms_out
      real(sp) :: amp(n), mean_val
      amp = sqrt(re_in**2 + im_in**2)
      mean_val = sum(amp)/real(n, sp)
      rms_out = sqrt(sum((amp - mean_val)**2)/real(n, sp))
   end subroutine rms_of_amp

   function wrap_chi0(x) result(w)
      real(sp), intent(in) :: x
      real(sp) :: w
      w = x - nint(x/pi_const)*pi_const
   end function wrap_chi0

   subroutine check(cond, label, all_pass_io)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: label
      logical, intent(inout) :: all_pass_io
      if (cond) then
         write(*,'(A,A)') '  [PASS]', label
      else
         write(*,'(A,A)') '  [FAIL]', label
         all_pass_io = .false.
      endif
   end subroutine check

end program test_matched_filter_refine
