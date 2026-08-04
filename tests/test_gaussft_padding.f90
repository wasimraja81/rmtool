program test_gaussft_padding
   !! Validates gaussft_mod's FFT-padding fix (next_fast_fft_size /
   !! plan_convolution's new nx_pad,ny_pad outputs / convolve_to_beam's
   !! new nx_pad,ny_pad arguments) -- added after a real ~46GB production
   !! run (4501x4501 = 7 x 643, 643 prime) took ~13min to convolve ONE
   !! cube, traced directly to FFTW falling back to a slow algorithm for
   !! that large prime factor (benchmarked: 2.76s/plane at n=4501 vs
   !! 1.26s/plane at the next 7-smooth size, 4608).
   !!
   !! Two checks:
   !! 1. next_fast_fft_size: a handful of known cases, including a size
   !!    that's ALREADY 7-smooth (must return itself unchanged, e.g. the
   !!    32x32 this project's own small test fixtures use, so existing
   !!    tests never actually exercised padding -- confirming why this
   !!    dedicated test exists).
   !! 2. Numerical equivalence: convolve_to_beam's own internal
   !!    zero-pad-then-crop, on a DELIBERATELY awkward image size (43,
   !!    prime -- guarantees next_fast_fft_size(43)=45 differs from 43,
   !!    so the padding path actually activates), must produce EXACTLY
   !!    the same result as manually zero-padding the same image to
   !!    45x45 BEFORE calling convolve_to_beam with nx=ny=nx_pad=ny_pad=45
   !!    (no auto-padding needed there, since the image already IS the
   !!    padded size) -- the two are mathematically the same operation,
   !!    so this proves the auto-padding code path does exactly what it
   !!    claims, not just "runs without crashing."
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use gaussft_mod
   implicit none

   integer :: n1, n2, n3, n4, n5
   logical :: all_pass

   all_pass = .true.

   ! --- Check 1: next_fast_fft_size ---
   n1 = next_fast_fft_size(4501)   ! 7 x 643 (643 prime) -> pads up
   n2 = next_fast_fft_size(32)     ! 2^5, already 7-smooth -> unchanged
   n3 = next_fast_fft_size(1)      ! degenerate edge case
   n4 = next_fast_fft_size(64)     ! 2^6, already 7-smooth -> unchanged
   n5 = next_fast_fft_size(43)     ! prime -> pads up

   call check_eq_int('next_fast_fft_size(4501) = 4536', n1, 4536, all_pass)
   call check_eq_int('next_fast_fft_size(32) = 32 (no-op, already fast)',&
   &n2, 32, all_pass)
   call check_eq_int('next_fast_fft_size(1) = 1', n3, 1, all_pass)
   call check_eq_int('next_fast_fft_size(64) = 64 (no-op, already fast)',&
   &n4, 64, all_pass)
   call check_eq_int('next_fast_fft_size(43) = 45', n5, 45, all_pass)

   ! --- Check 2: auto-pad-inside-convolve_to_beam matches manual padding ---
   call check_padding_equivalence(all_pass)

   if (all_pass) then
      print '(A)', '[PASS] test_gaussft_padding: all checks passed'
   else
      print '(A)', '[FAIL] test_gaussft_padding: one or more checks failed'
      stop 1
   endif

contains

   subroutine check_eq_int(label, got, want, ok)
      character(len=*), intent(in) :: label
      integer, intent(in) :: got, want
      logical, intent(inout) :: ok
      if (got.eq.want) then
         print '(A,A)', '[PASS] ', label
      else
         print '(A,A,A,I0,A,I0)', '[FAIL] ', label, ' -- got ', got,&
         &', want ', want
         ok = .false.
      endif
   end subroutine check_eq_int

   subroutine check_padding_equivalence(ok)
      logical, intent(inout) :: ok
      integer, parameter :: n = 43, n_pad_expected = 45
      real(dp), parameter :: dx = 5.0e-4_dp, dy = 5.0e-4_dp
      ! Native (source) PSF, and a larger target PSF -- genuine
      ! convolution, not a near-identity no-op.
      real(dp), parameter :: bmaj_in = 10.0_dp/3600.0_dp, bmin_in = 8.0_dp/3600.0_dp
      real(dp), parameter :: bpa_in = 15.0_dp
      real(dp), parameter :: bmaj = 20.0_dp/3600.0_dp, bmin = 16.0_dp/3600.0_dp
      real(dp), parameter :: bpa = 40.0_dp

      real(dp) :: image_native(n, n), image_padded_manual(n_pad_expected, n_pad_expected)
      real(dp) :: out_auto(n, n), out_manual_full(n_pad_expected, n_pad_expected)
      real(dp) :: out_identity(n, n)
      integer(kind=8) :: plan_fwd_native, plan_bwd_native
      integer(kind=8) :: plan_fwd_manual, plan_bwd_manual
      integer :: nx_pad_native, ny_pad_native, nx_pad_manual, ny_pad_manual
      integer :: status_auto, status_manual, status_identity
      real(dp) :: max_abs_diff, max_abs_diff_identity

      image_native = 0.0_dp
      ! An off-center point source (not the padding corner itself) --
      ! a symmetric center-of-array source wouldn't distinguish a
      ! padding-corner bug from a correct implementation as sharply.
      image_native(15, 20) = 1.0_dp

      image_padded_manual = 0.0_dp
      image_padded_manual(1:n, 1:n) = image_native

      call plan_convolution(n, n, plan_fwd_native, plan_bwd_native,&
      &nx_pad_native, ny_pad_native)
      call plan_convolution(n_pad_expected, n_pad_expected, plan_fwd_manual,&
      &plan_bwd_manual, nx_pad_manual, ny_pad_manual)

      call check_eq_int('plan_convolution(43,43) auto-pads to 45x45',&
      &nx_pad_native, n_pad_expected, ok)
      call check_eq_int('plan_convolution(45,45,...) is already fast (no pad)',&
      &nx_pad_manual, n_pad_expected, ok)

      ! Auto-padding path: image at its own true (43,43) extent,
      ! convolve_to_beam pads internally to (45,45) and crops back.
      call convolve_to_beam(plan_fwd_native, plan_bwd_native, image_native,&
      &n, n, nx_pad_native, ny_pad_native, dx, dy, bmaj_in, bmin_in, bpa_in,&
      &bmaj, bmin, bpa, out_auto, status_auto)

      ! Manual-padding reference: image ALREADY at the padded (45,45)
      ! extent (zero border added by hand above) -- no auto-padding
      ! needed (nx=ny=nx_pad=ny_pad=45), so this is a completely
      ! independent code path through the same routine.
      call convolve_to_beam(plan_fwd_manual, plan_bwd_manual,&
      &image_padded_manual, n_pad_expected, n_pad_expected, nx_pad_manual,&
      &ny_pad_manual, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa,&
      &out_manual_full, status_manual)

      if (status_auto.ne.0 .or. status_manual.ne.0) then
         print '(A)', '[FAIL] convolve_to_beam returned nonzero status'
         ok = .false.
         return
      endif

      max_abs_diff = maxval(abs(out_auto - out_manual_full(1:n, 1:n)))
      if (max_abs_diff.lt.1.0e-12_dp) then
         print '(A,A,ES10.3,A)', '[PASS] auto-pad-inside-convolve_to_beam',&
         &' matches manual pre-padding to the bit (max|diff|=', max_abs_diff, ')'
      else
         print '(A,ES10.3)', '[FAIL] auto-pad result diverges from manual'//&
         &' pre-padding reference, max|diff|=', max_abs_diff
         ok = .false.
      endif

      ! Independent physical sanity check, exact rather than approximate:
      ! target beam == native beam must be a no-op (G_ratio(u,v) == 1
      ! everywhere, so IFFT(FFT(image)*1) == image exactly, up to FFT
      ! round-off -- true regardless of image content or padding, unlike
      ! a flux-conservation check on a raw delta-function input, which
      ! this deconvolve/reconvolve operation was never designed for: it
      ! assumes its input already HAS the native beam's own shape baked
      ! in, not a bare unconvolved point source).
      call convolve_to_beam(plan_fwd_native, plan_bwd_native, image_native,&
      &n, n, nx_pad_native, ny_pad_native, dx, dy, bmaj_in, bmin_in, bpa_in,&
      &bmaj_in, bmin_in, bpa_in, out_identity, status_identity)

      call destroy_convolution_plan(plan_fwd_native, plan_bwd_native)
      call destroy_convolution_plan(plan_fwd_manual, plan_bwd_manual)

      if (status_identity.ne.0) then
         print '(A)', '[FAIL] identity-case convolve_to_beam returned nonzero status'
         ok = .false.
         return
      endif
      ! Tolerance loosened from the original 1.0e-9 (T27, docs/dev/
      ! MULTI_BAND_TOMOGRAPHY_PLAN.md): convolve_to_beam's FFT buffers
      ! are now single precision internally, so the achievable round-trip
      ! accuracy is bounded by float32 epsilon (~1.19e-7), not double
      ! precision -- 1.0e-9 is no longer achievable by design, not a
      ! regression. 1.0e-6 stays comfortably above the actually-observed
      ! error (~3.1e-8) while still tight enough to catch a real fault.
      max_abs_diff_identity = maxval(abs(out_identity - image_native))
      if (max_abs_diff_identity.lt.1.0e-6_dp) then
         print '(A,A,ES10.3,A)', '[PASS] target beam == native beam is a',&
         &' no-op (max|diff|=', max_abs_diff_identity, ')'
      else
         print '(A,ES10.3)', '[FAIL] target beam == native beam should be'//&
         &' a no-op, max|diff|=', max_abs_diff_identity
         ok = .false.
      endif
   end subroutine check_padding_equivalence

end program test_gaussft_padding
