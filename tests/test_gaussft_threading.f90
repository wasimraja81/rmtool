program test_gaussft_threading
   !! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): validates
   !! gaussft_mod's new in-plane threading (plan_convolution's own
   !! nthreads argument -- sfftw_plan_with_nthreads, plus the internal
   !! !$omp parallel do/workshare constructs added to convolve_to_beam)
   !! against the pre-T28 nthreads=1 behaviour, and against itself for
   !! run-to-run determinism.
   !!
   !! Three checks, all against a moderately large (401x401 -- big
   !! enough that OMP_NUM_THREADS>1 genuinely engages multiple threads
   !! across cores, unlike this project's tiny 32x32 fixtures) synthetic
   !! plane with several point sources, run through the SAME
   !! has_nan=.true. path T23 exercises (some pixels deliberately NaN'd)
   !! since that is the more complex of convolve_to_beam's two branches:
   !!
   !! 1. nthreads=4 matches nthreads=1 within a tolerance appropriate to
   !!    single precision (FFTW's own internal reduction order can
   !!    legitimately differ between thread counts -- floating-point
   !!    addition is not associative -- so this is NOT a bit-exact
   !!    check, unlike check 2 below).
   !! 2. nthreads=4 run TWICE gives BIT-IDENTICAL output both times --
   !!    confirms a given plan's own threaded decomposition is
   !!    deterministic across runs, not dependent on run-to-run OS
   !!    scheduling variance (asserted by FFTW's own documentation, not
   !!    otherwise verified anywhere in this project until now).
   !! 3. nthreads=4 output is not degenerate (has real, finite,
   !!    non-constant values away from the deliberately-NaN'd region) --
   !!    a cheap sanity check that the threaded path did not silently
   !!    zero out or corrupt the computation.
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use gaussft_mod
   implicit none

   integer, parameter :: n = 401
   real(dp), parameter :: dx = 5.0e-4_dp, dy = 5.0e-4_dp
   real(dp), parameter :: bmaj_in = 10.0_dp/3600.0_dp, bmin_in = 8.0_dp/3600.0_dp
   real(dp), parameter :: bpa_in = 15.0_dp
   real(dp), parameter :: bmaj = 20.0_dp/3600.0_dp, bmin = 16.0_dp/3600.0_dp
   real(dp), parameter :: bpa = 40.0_dp

   real(dp) :: image(n, n)
   real(dp) :: out_1thread(n, n), out_4thread_a(n, n), out_4thread_b(n, n)
   integer(kind=8) :: plan_fwd_1, plan_bwd_1, plan_fwd_4, plan_bwd_4
   integer :: nx_pad_1, ny_pad_1, nx_pad_4, ny_pad_4
   integer :: status_1, status_4a, status_4b
   real(dp) :: max_abs_diff, max_abs_diff_repeat
   logical :: all_pass
   integer :: ix, iy, n_nan_block, n_finite

   all_pass = .true.

   ! Several point sources plus a deliberately-NaN'd block (same shape
   ! of test T23's own regression uses -- a real plane in flight always
   ! has some NaN, so a threading test that never exercises the
   ! has_nan=.true. branch would miss the more complex code path.
   image = 0.0_dp
   image(80, 120) = 1.0_dp
   image(200, 200) = 0.6_dp
   image(310, 90) = 0.3_dp
   n_nan_block = 0
   do iy = 150, 170
      do ix = 150, 170
         image(ix, iy) = ieee_value_nan()
         n_nan_block = n_nan_block + 1
      enddo
   enddo

   call plan_convolution(n, n, plan_fwd_1, plan_bwd_1, nx_pad_1, ny_pad_1, 1)
   call convolve_to_beam(plan_fwd_1, plan_bwd_1, image, n, n, nx_pad_1,&
   &ny_pad_1, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa,&
   &out_1thread, status_1)
   call destroy_convolution_plan(plan_fwd_1, plan_bwd_1)

   call plan_convolution(n, n, plan_fwd_4, plan_bwd_4, nx_pad_4, ny_pad_4, 4)
   call convolve_to_beam(plan_fwd_4, plan_bwd_4, image, n, n, nx_pad_4,&
   &ny_pad_4, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa,&
   &out_4thread_a, status_4a)
   ! Second call against the SAME plan -- determinism check, not a
   ! fresh plan (re-planning is not what's being tested here).
   call convolve_to_beam(plan_fwd_4, plan_bwd_4, image, n, n, nx_pad_4,&
   &ny_pad_4, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa,&
   &out_4thread_b, status_4b)
   call destroy_convolution_plan(plan_fwd_4, plan_bwd_4)

   if (status_1.ne.0 .or. status_4a.ne.0 .or. status_4b.ne.0) then
      print '(A)', '[FAIL] convolve_to_beam returned nonzero status'
      all_pass = .false.
   endif

   ! Check 1: nthreads=4 matches nthreads=1 within a tolerance
   ! appropriate to single precision + legitimately different FP
   ! reduction order between thread counts (not bit-exact by design --
   ! see this program's own header comment). Compared only where BOTH
   ! sides are finite (the deliberately-NaN'd block should read NaN on
   ! both sides identically -- see check 3 below for that).
   max_abs_diff = 0.0_dp
   do iy = 1, n
      do ix = 1, n
         if (.not.ieee_is_nan(out_1thread(ix,iy)) .and.&
         &.not.ieee_is_nan(out_4thread_a(ix,iy))) then
            max_abs_diff = max(max_abs_diff,&
            &abs(out_1thread(ix,iy) - out_4thread_a(ix,iy)))
         endif
      enddo
   enddo
   if (max_abs_diff.lt.1.0e-5_dp) then
      print '(A,ES10.3,A)', '[PASS] nthreads=4 matches nthreads=1 (max|diff|=',&
      &max_abs_diff, ')'
   else
      print '(A,ES10.3)', '[FAIL] nthreads=4 diverges from nthreads=1, max|diff|=',&
      &max_abs_diff
      all_pass = .false.
   endif

   ! Check 2: nthreads=4 run twice against the SAME plan is
   ! bit-identical -- genuine determinism, not "close enough".
   max_abs_diff_repeat = maxval(abs(out_4thread_a - out_4thread_b),&
   &mask=(.not.ieee_is_nan(out_4thread_a)) .and. (.not.ieee_is_nan(out_4thread_b)))
   if (max_abs_diff_repeat.eq.0.0_dp) then
      print '(A)', '[PASS] nthreads=4 run twice is bit-identical (deterministic)'
   else
      print '(A,ES10.3)', '[FAIL] nthreads=4 run twice diverged, max|diff|=',&
      &max_abs_diff_repeat
      all_pass = .false.
   endif

   ! Check 3: not degenerate -- deep inside the NaN block stays NaN
   ! (threshold rejection still engaging under threading), and there
   ! is real, finite, non-constant data away from it.
   if (.not.ieee_is_nan(out_4thread_a(160, 160))) then
      print '(A)', '[FAIL] pixel deep inside the NaN block came back non-NaN under threading'
      all_pass = .false.
   else
      print '(A)', '[PASS] pixel deep inside the NaN block correctly remains NaN under threading'
   endif
   n_finite = count(.not.ieee_is_nan(out_4thread_a))
   if (n_finite.gt.(n*n)/2 .and. maxval(out_4thread_a, mask=.not.ieee_is_nan(out_4thread_a)).gt.0.0_dp) then
      print '(A,I0,A)', '[PASS] output is not degenerate (', n_finite, ' finite pixels, real positive peak present)'
   else
      print '(A)', '[FAIL] output looks degenerate under threading'
      all_pass = .false.
   endif

   if (all_pass) then
      print '(A)', '[PASS] test_gaussft_threading: all checks passed'
   else
      print '(A)', '[FAIL] test_gaussft_threading: one or more checks failed'
      stop 1
   endif

contains

   function ieee_value_nan() result(v)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: v
      v = ieee_value(1.0_dp, ieee_quiet_nan)
   end function ieee_value_nan

end program test_gaussft_threading
