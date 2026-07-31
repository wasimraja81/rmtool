module rmclean_mod
   !! RM-CLEAN core: Hogbom-style deconvolution of a complex Faraday
   !! dispersion function (dirty FDF) against its own RM-dependent dirty
   !! beam (RMSF). Pure computation only -- no FITS I/O -- mirroring
   !! gaussft_mod's own split between narrowly-scoped computation and a
   !! caller that owns I/O/config. Ported and modernized from the user's
   !! own thesis Fortran77 code (~/softwares/CURR_DEVEL/RM_CLEAN_TESTS/),
   !! not a raw port -- see planning/RMCLEAN_INTEGRATION_PLAN.md for the
   !! full design record (decisions confirmed with the user) this module
   !! implements, including that file's own "Choosing parameters" section
   !! -- start there for what to set lsq_ref_compute/lsq_ref_report/
   !! oversample/thresh to; this module has no config layer or defaults
   !! of its own yet (T2).
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   implicit none
   private
   public :: index_absmax, peak_interp_parabolic
   public :: plan_fourier_interp, destroy_fourier_interp_plan, fourier_interp_complex
   public :: rmsf_table_t, build_rmsf_offset_table, destroy_rmsf_offset_table
   public :: compute_dirty_rmbeam_direct, compute_dirty_rmbeam
   public :: rmsf_point_direct, refine_peak_matched_filter
   public :: clean_complex
   public :: compute_rmsf_fwhm, compute_rmsf_fwhm_multiband, restore_clean
   public :: get_drm, derotate_to_lsq_ref
   public :: get_lsq_ref_compute, get_lsq_ref_report
   public :: lsq_ref_compute_mid, lsq_ref_compute_intrinsic, lsq_ref_compute_centroid
   public :: lsq_ref_compute_min, lsq_ref_compute_max, lsq_ref_compute_fixed
   public :: lsq_ref_report_intrinsic, lsq_ref_report_centroid, lsq_ref_report_mid
   public :: lsq_ref_report_min, lsq_ref_report_max, lsq_ref_report_fixed

   ! FFTW3 constants (same convention as gaussft_mod's own -- declared
   ! directly rather than `include`d, fftw3.f is fixed-form F77 and
   ! cannot be included into a free-form .f90 file directly).
   integer, parameter :: fftw_forward = -1
   integer, parameter :: fftw_backward = 1
   integer, parameter :: fftw_estimate = 64

   ! get_lsq_ref_compute's own recognized modes -- the RECOMMENDED,
   ! default choice is lsq_ref_compute_mid (=1, listed first
   ! deliberately): it minimizes get_drm's own grid-cost bound (a
   ! genuine, provable property -- get_lsq_ref_compute's own comment)
   ! and, once comp_rm_refined/derotate_to_lsq_ref are used correctly,
   ! costs NOTHING in chi0 precision -- so there is no accuracy trade-off
   ! for choosing it, only a computational win. The rest exist so a
   ! caller can experiment or match some other external convention;
   ! none of them changes what get_drm subsequently enforces, since
   ! get_drm always derives its own bound from whichever
   ! lsq_ref_compute value is actually in use, regardless of mode.
   integer, parameter :: lsq_ref_compute_mid       = 1
   integer, parameter :: lsq_ref_compute_intrinsic = 2
   integer, parameter :: lsq_ref_compute_centroid  = 3
   integer, parameter :: lsq_ref_compute_min       = 4
   integer, parameter :: lsq_ref_compute_max       = 5
   integer, parameter :: lsq_ref_compute_fixed     = 6

   ! get_lsq_ref_report's own recognized modes -- lsq_ref_report_
   ! intrinsic (=1, listed first) is this project's own default
   ! reporting convention (matches the user's own thesis).
   integer, parameter :: lsq_ref_report_intrinsic = 1
   integer, parameter :: lsq_ref_report_centroid  = 2
   integer, parameter :: lsq_ref_report_min       = 3
   integer, parameter :: lsq_ref_report_max       = 4
   integer, parameter :: lsq_ref_report_mid       = 5
   integer, parameter :: lsq_ref_report_fixed     = 6

   ! Speed of light, Mm/s (matches rm_synthesis_mod.f90's own c_velocity
   ! exactly -- this module deliberately doesn't `use` that module, same
   ! standalone-computation-module precedent as gaussft_mod/commonbeam_mod,
   ! so the constant is duplicated here rather than imported).
   real(sp), parameter :: c_velocity = 299.792458_sp

   type :: rmsf_table_t
      !! Pre-BW-depol dirty-beam (RMSF) offset table -- planning/
      !! RMCLEAN_INTEGRATION_PLAN.md decisions 3/9. Derivation (worked
      !! through directly against compute_dirty_rmbeam.f's own dot-
      !! product formula, cross-checked by compute_dirty_rmbeam_direct
      !! below, not just asserted): with RM-independent per-channel
      !! weights (no BW-depol), the beam at trial RM=rm_samp(j) for a
      !! CLEAN component at (rm_in,phase_in) is
      !!   re_beam(j)+i*im_beam(j) = exp(2i*phase_in) * R(delta_j)
      !! with delta_j = rm_samp(j)-rm_in and
      !!   R(delta) = (1/nchan) * sum_k exp(-2i*delta*(L_sq(k)-lsq_ref_compute))
      !! -- i.e. the ENTIRE (rm_samp(j), rm_in) dependence collapses to
      !! one scalar offset delta, and phase_in enters only as a global
      !! post-multiply, not inside the sum at all. R(delta) is built
      !! once, on a grid finer than rm_samp's own spacing, spanning
      !! every offset a real (rm_in, rm_samp(j)) pair can produce;
      !! compute_dirty_rmbeam then costs one table lookup per (j) rather
      !! than an O(nchan) sum. NOT valid once BW-depol correction is
      !! reintroduced (planning doc T-future): per-channel weights
      !! become rm_in-dependent there, breaking the pure-offset
      !! structure this table relies on -- that will need its own
      !! design, not an extension of this one.
      integer :: n_fine = 0
      real(dp) :: delta_min = 0.0_dp, ddelta = 0.0_dp
      real(dp), allocatable :: re_fine(:), im_fine(:)
   end type rmsf_table_t

contains

   subroutine rms_about_mean(v, n, rms_val)
      !! Population RMS about the mean -- private, local copy of
      !! rm_synthesis_mod.f90's own compute_rms (added there for the
      !! eventual inline-into-rm_synthesis path, planning/
      !! RMCLEAN_INTEGRATION_PLAN.md decision 11's "later" goal). This
      !! module doesn't `use` rm_synthesis_mod at all (same standalone
      !! precedent as gaussft_mod/commonbeam_mod, see c_velocity's own
      !! comment above), so its own CLEAN loop needs a private copy for
      !! today's standalone-tool architecture -- a 2-line utility, not
      !! the kind of duplication planning doc decision 7 was about
      !! avoiding (compute_mean/dot_product_custom, both already
      !! nontrivial and already reused as-is, not touched here).
      integer, intent(in) :: n
      real(sp), intent(in) :: v(n)
      real(sp), intent(out) :: rms_val
      real(sp) :: mean_val
      mean_val = sum(v)/real(n, sp)
      rms_val = sqrt(sum((v-mean_val)**2)/real(n, sp))
   end subroutine rms_about_mean

   subroutine index_absmax(v, n, imax, avg_abs)
      !! Index of the largest |v(i)| and the mean of |v| over the whole
      !! array -- ported from index_absmax.f, rewritten with whole-array
      !! intrinsics (maxloc/sum) rather than a hand-rolled accumulator
      !! loop: at minimum as fast under -O3 vectorization, and reads as
      !! one expression each instead of a loop with a running-max check.
      integer, intent(in) :: n
      real(sp), intent(in) :: v(n)
      integer, intent(out) :: imax
      real(sp), intent(out) :: avg_abs

      imax = maxloc(abs(v), dim=1)
      avg_abs = sum(abs(v)) / real(n, sp)
   end subroutine index_absmax

   subroutine peak_interp_parabolic(re_in, im_in, n, i_center, peak_offset,&
   &re_peak, im_peak)
      !! Sub-pixel peak refinement by parabolic (quadratic) interpolation,
      !! extended per planning/RMCLEAN_INTEGRATION_PLAN.md decision 4:
      !! the ORIGINAL quad_interp.f only ever interpolated a magnitude
      !! array, leaving phase undefined at every mode except complex
      !! Fourier interpolation (a self-flagged incomplete state in that
      !! code, "Incorporate the scheme NOW!!!"). Here, the sub-pixel
      !! OFFSET is determined once from the log-magnitude of re_in/im_in
      !! (same numerically-stable log10 compression as the original, to
      !! tame dynamic range), then Re and Im are fit their OWN parabolas
      !! independently and evaluated at that SAME offset -- giving a
      !! genuine interpolated complex value (amplitude AND phase) instead
      !! of magnitude alone.
      !!
      !! re_in/im_in: n values (n>=3) centred on the sampled peak at
      !! index i_center (1-indexed) -- only the 3 samples immediately
      !! around i_center are used (i_center-1, i_center, i_center+1);
      !! n/the full arrays are accepted so a caller can pass a small
      !! subarray view without slicing, matching how rm_clean's own
      !! Hogbom loop already has a windowed subarray in hand around its
      !! current peak.
      !!
      !! peak_offset: the fitted vertex location, in fractional bins
      !! relative to i_center (typically in [-0.5,0.5]; the caller adds
      !! this to i_center to get the true fractional peak index -- kept
      !! separate rather than folded into an absolute index here, since
      !! that absolute-index/local-offset conflation is exactly the bug
      !! this port fixes, see below).
      !!
      !! Bug fixed vs. the original quad_interp.f: that code's own
      !! peak_val formula (`peak_val = beta - 0.25*(alpha-gama)*peak`)
      !! used `peak` -- by then reassigned to the ABSOLUTE fractional
      !! index (i_center + offset) -- where the correct closed-form
      !! vertex-value formula requires the LOCAL offset alone (derivable
      !! directly from the three-point parabola y=a(x-p)^2+b fit through
      !! (-1,alpha),(0,beta),(1,gamma): b = beta - 0.25*(alpha-gamma)*p,
      !! p being the local vertex offset, not i_center+p). Using the
      !! absolute index there is only numerically close to correct when
      !! i_center happens to be small (true in the original's own usage,
      !! always called on a 5-point window with i_center in [1,5]); it is
      !! not correct in general, and not carried forward here.
      integer, intent(in) :: n, i_center
      real(sp), intent(in) :: re_in(n), im_in(n)
      real(sp), intent(out) :: peak_offset, re_peak, im_peak

      real(sp) :: mag(3), alpha, beta, gama
      real(sp) :: alpha_re, beta_re, gama_re, alpha_im, beta_im, gama_im
      real(sp) :: abs_val(3)

      ! log10(0) is -Infinity -- not a numerical inconvenience to smooth
      ! over, a genuine domain violation (log of zero is mathematically
      ! undefined), so it cannot be guarded with an arbitrary small-
      ! number floor without silently fabricating curvature that isn't
      ! there. An exact zero here is a legitimate outcome of iterative
      ! floating-point subtraction (CLEAN can drive a residual bin to
      ! exactly 0.0 after enough iterations, confirmed directly: this
      ! was reached in practice, not a hypothetical). If any of the 3
      ! samples is exactly zero, there is no well-defined log-magnitude
      ! curvature to fit at all -- fall back to the sampled centre
      ! itself, the same principled fallback already used just below for
      ! the flat/linear (zero second-difference) case.
      abs_val = sqrt(re_in(i_center-1:i_center+1)**2 + im_in(i_center-1:i_center+1)**2)
      if (any(abs_val <= 0.0_sp)) then
         peak_offset = 0.0_sp
         re_peak = re_in(i_center)
         im_peak = im_in(i_center)
         return
      endif

      mag = log10(abs_val)
      alpha = mag(1)
      beta = mag(2)
      gama = mag(3)

      if (alpha + gama - 2.0_sp*beta == 0.0_sp) then
         ! Degenerate (perfectly flat or perfectly linear across the
         ! window) -- no curvature to fit a vertex to; the sampled centre
         ! is already the best estimate.
         peak_offset = 0.0_sp
      else
         peak_offset = 0.5_sp*(alpha-gama)/(alpha+gama-2.0_sp*beta)
      endif

      ! Evaluate Re's and Im's OWN quadratic fit (each through its own 3
      ! samples at x=-1,0,1) at x=peak_offset -- the exact quadratic
      ! through 3 equally-spaced points is y(x)=A*x^2+B*x+C with
      ! A=0.5*(alpha-2*beta+gama), B=0.5*(gama-alpha), C=beta (standard
      ! central-difference coefficients, exact for 3 points, not an
      ! approximation). This is deliberately NOT the vertex-value
      ! shortcut (b=beta-0.25*(alpha-gama)*p) quad_interp.f itself uses:
      ! that shortcut only gives the value AT THE CHANNEL'S OWN vertex,
      ! not at an externally supplied offset like peak_offset here
      ! (verified by direct derivation -- the two formulas coincide only
      ! when evaluating a channel at its own fitted peak).
      alpha_re = re_in(i_center-1); beta_re = re_in(i_center); gama_re = re_in(i_center+1)
      alpha_im = im_in(i_center-1); beta_im = im_in(i_center); gama_im = im_in(i_center+1)

      re_peak = 0.5_sp*(alpha_re - 2.0_sp*beta_re + gama_re)*peak_offset**2 +&
      &0.5_sp*(gama_re - alpha_re)*peak_offset + beta_re
      im_peak = 0.5_sp*(alpha_im - 2.0_sp*beta_im + gama_im)*peak_offset**2 +&
      &0.5_sp*(gama_im - alpha_im)*peak_offset + beta_im
   end subroutine peak_interp_parabolic

   subroutine plan_fourier_interp(npts, nout, plan_fwd, plan_bwd)
      !! Create the two FFTW plans (forward, size npts; backward, size
      !! nout) a fourier_interp_complex call for this (npts,nout) pair
      !! needs, once, to be reused across every subsequent call -- same
      !! plan-once/execute-many split as gaussft_mod's own
      !! plan_convolution, for the same reason (FFTW's planner is not
      !! thread-safe and must never run inside a parallel region; the
      !! "new-array execute" form used below is safe to call
      !! concurrently from multiple threads against a shared plan, each
      !! supplying its own arrays). MUST be called serially. Unlike the
      !! original fourier_interp.f (built on a home-grown radix-2-only
      !! FFT, forcing nout up to the next power of 2 above npts*ofac),
      !! FFTW handles arbitrary sizes directly -- nout can be exactly
      !! npts*ofac, no padding-then-truncating needed.
      integer, intent(in) :: npts, nout
      integer(kind=8), intent(out) :: plan_fwd, plan_bwd

      complex(dp), allocatable :: scratch(:)

      allocate(scratch(npts))
      call dfftw_plan_dft_1d(plan_fwd, npts, scratch, scratch, fftw_forward, fftw_estimate)
      deallocate(scratch)

      allocate(scratch(nout))
      call dfftw_plan_dft_1d(plan_bwd, nout, scratch, scratch, fftw_backward, fftw_estimate)
      deallocate(scratch)
   end subroutine plan_fourier_interp

   subroutine destroy_fourier_interp_plan(plan_fwd, plan_bwd)
      integer(kind=8), intent(inout) :: plan_fwd, plan_bwd
      call dfftw_destroy_plan(plan_fwd)
      call dfftw_destroy_plan(plan_bwd)
   end subroutine destroy_fourier_interp_plan

   subroutine fourier_interp_complex(plan_fwd, plan_bwd, re_in, im_in, npts,&
   &nout, re_out, im_out)
      !! Bandlimited (Fourier) interpolation of a complex signal from
      !! npts to nout equally-spaced samples over the same span: FFT the
      !! input, zero-pad the spectrum (insert nout-npts zeros at the
      !! Nyquist split, preserving low frequencies at both ends -- the
      !! standard "zero-pad in the frequency domain" upsampling
      !! construction), inverse FFT at size nout, scale by nout/npts.
      !! Exact for any signal that is itself bandlimited within the
      !! original npts-point Nyquist band (e.g. reproduces a pure
      !! sinusoid's true value at every interpolated point, not just an
      !! approximation) -- ported from fourier_interp.f, same
      !! zero-padding construction, FFTW instead of the original's
      !! home-grown fft1d/ifft1d (see plan_fourier_interp's own comment).
      !! plan_fwd/plan_bwd: from plan_fourier_interp, already sized for
      !! exactly this (npts,nout) pair -- not checked here, matching
      !! FFTW's own new-array-execute contract (passing plans for a
      !! different size is undefined behaviour).
      integer(kind=8), intent(in) :: plan_fwd, plan_bwd
      integer, intent(in) :: npts, nout
      real(sp), intent(in) :: re_in(npts), im_in(npts)
      real(sp), intent(out) :: re_out(nout), im_out(nout)

      complex(dp) :: cin(npts), cout(nout)
      integer :: n1, npad, i, j
      real(dp) :: norm

      cin = cmplx(real(re_in, dp), real(im_in, dp), dp)
      call dfftw_execute_dft(plan_fwd, cin, cin)

      cout = (0.0_dp, 0.0_dp)
      if (mod(npts, 2) == 0) then
         n1 = npts/2 + 1
      else
         n1 = (npts+1)/2
      endif
      npad = nout - npts

      cout(1:n1) = cin(1:n1)
      j = n1
      do i = n1+npad+1, nout
         j = j + 1
         cout(i) = cin(j)
      enddo

      call dfftw_execute_dft(plan_bwd, cout, cout)

      ! FFTW's dfftw_execute_dft is UNNORMALIZED in both directions
      ! (forward-then-backward returns N times the original signal, not
      ! the original signal) -- unlike the original fourier_interp.f's
      ! own home-grown ifft1d, which normalizes internally by 1/N as
      ! part of the inverse transform itself (the common convention for
      ! a hand-rolled FFT/IFFT pair). That's why the original code's own
      ! scale factor was nout/npts, not 1/npts: its own ifft1d already
      ! divided by nout internally, so nout/npts on top of that lands at
      ! the correct overall 1/npts. With FFTW doing no internal
      ! normalization at all, the only scale factor needed here is
      ! 1/npts directly.
      norm = 1.0_dp/real(npts, dp)
      re_out = real(real(cout, dp)*norm, sp)
      im_out = real(aimag(cout)*norm, sp)
   end subroutine fourier_interp_complex

   subroutine compute_dirty_rmbeam_direct(l_sq, nchan, rm_in, phase_in,&
   &cos_arr, sin_arr, nrm, maxrm, maxchan, lsq_ref_compute, re_beam, im_beam)
      !! Exact, O(nrm*nchan) reference implementation -- a direct port of
      !! compute_dirty_rmbeam.f's own dot-product formula (matched-filter
      !! correlation of a synthetic point-source signal at (rm_in,
      !! phase_in) against the SAME cos_arr/sin_arr templates
      !! rm_synthesis_mod.f90's own extract_general_setup already builds
      !! and reuses across the whole run -- planning doc decision 5, no
      !! duplicate template construction here). No BW-depol (that whole
      !! path is deferred, T-future -- see rmsf_table_t's own comment).
      !! Kept as the ground-truth compute_dirty_rmbeam is verified
      !! against, not meant for the CLEAN loop's own hot path.
      !!
      !! lsq_ref_compute: the phase-reference lambda-squared, EXPLICITLY supplied
      !! by the caller rather than computed internally as mean(l_sq) --
      !! this is a pure convention choice (the derivation in
      !! rmsf_table_t's own comment holds for ANY fixed reference point,
      !! not specifically the mean), and different callers legitimately
      !! want different choices: rm_synthesis_mod.f90's own
      !! extract_general_setup is UNCONDITIONALLY at lambda_sq=0 (its own
      !! phi_tmp=omega*t(kk) uses raw L_sq, no mean subtraction anywhere
      !! -- confirmed directly at rm_synthesis_mod.f90:675-791), matching
      !! the user's own thesis codebase exactly; this module's own
      !! get_lsq_ref_compute(mode=centroid) is the real example of a
      !! mean-style reference in this codebase -- a caller's CHOICE for
      !! compute cost (numerically smaller phase arguments), not
      !! something rm_synthesis's own dirty-map construction already
      !! does. Whatever the caller passes here MUST match whatever
      !! reference point was used to construct cos_arr/sin_arr AND the
      !! caller's own dirty-map/injected-data phase convention -- a
      !! mismatch here does not break the amplitude in a broad "noisy"
      !! way, it introduces a genuine phase distortion that differs
      !! per-component (verified directly: forcing two different
      !! reference conventions to disagree between an injected multi-
      !! component sky model and its own dirty-map construction gave
      !! wrong relative amplitudes between components, not just a
      !! rotated-but-otherwise-correct result).
      integer, intent(in) :: nchan, nrm, maxrm, maxchan
      real(sp), intent(in) :: l_sq(nchan), rm_in, phase_in
      real(sp), intent(in) :: cos_arr(maxrm, maxchan), sin_arr(maxrm, maxchan)
      real(sp), intent(in) :: lsq_ref_compute
      real(sp), intent(out) :: re_beam(nrm), im_beam(nrm)

      real(sp) :: phi_tmp
      real(sp) :: ryt(nchan), iyt(nchan)
      real(sp) :: rc_cor, rs_cor, ic_cor, is_cor
      integer :: j

      block
         integer :: kk
         do kk = 1, nchan
            phi_tmp = 2.0_sp*(rm_in*(l_sq(kk)-lsq_ref_compute) + phase_in)
            ryt(kk) = cos(phi_tmp)
            iyt(kk) = sin(phi_tmp)
         end do
      end block

      do j = 1, nrm
         rc_cor = dot_product(ryt, cos_arr(j,1:nchan)) / real(nchan, sp)
         rs_cor = dot_product(ryt, sin_arr(j,1:nchan)) / real(nchan, sp)
         ic_cor = dot_product(iyt, cos_arr(j,1:nchan)) / real(nchan, sp)
         is_cor = dot_product(iyt, sin_arr(j,1:nchan)) / real(nchan, sp)
         re_beam(j) = rc_cor - is_cor
         im_beam(j) = rs_cor + ic_cor
      end do
   end subroutine compute_dirty_rmbeam_direct

   subroutine build_rmsf_offset_table(l_sq, nchan, lsq_ref_compute, delta_span,&
   &native_ddelta, oversample, table)
      !! Build R(delta) once (see rmsf_table_t's own comment for the
      !! derivation) on a grid finer than the dirty map's own RM sampling
      !! by `oversample`, spanning [-delta_span,+delta_span] -- the
      !! caller passes delta_span = rm_samp(nrm)-rm_samp(1), the largest
      !! offset a (rm_in, rm_samp(j)) pair inside the search domain can
      !! ever produce, so every lookup this table will ever see falls
      !! inside its built range.
      !!
      !! lsq_ref_compute: MUST match whatever reference point the caller's own
      !! dirty-map/injected-data construction uses -- see
      !! compute_dirty_rmbeam_direct's own comment on why this is a real
      !! correctness requirement, not a cosmetic choice: R(delta) itself
      !! is only a pure function of delta for a FIXED, shared reference
      !! point; changing the reference multiplies R(delta) by an extra
      !! delta-dependent phase factor, exp(-2i*delta*(new_ref-old_ref)),
      !! not just a constant rotation.
      real(sp), intent(in) :: l_sq(nchan)
      integer, intent(in) :: nchan
      real(sp), intent(in) :: lsq_ref_compute
      real(sp), intent(in) :: delta_span, native_ddelta
      integer, intent(in) :: oversample
      type(rmsf_table_t), intent(out) :: table

      real(dp) :: lsq_ref_dp, delta_m, phase_k
      integer :: m

      lsq_ref_dp = real(lsq_ref_compute, dp)

      table%ddelta = real(native_ddelta, dp)/real(oversample, dp)
      table%n_fine = 2*ceiling(real(delta_span, dp)/table%ddelta) + 1
      table%delta_min = -0.5_dp*real(table%n_fine-1, dp)*table%ddelta
      allocate(table%re_fine(table%n_fine), table%im_fine(table%n_fine))

      !$omp parallel do default(shared) private(m, delta_m, phase_k)
      do m = 1, table%n_fine
         delta_m = table%delta_min + real(m-1, dp)*table%ddelta
         phase_k = -2.0_dp*delta_m
         table%re_fine(m) = sum(cos(phase_k*(real(l_sq, dp)-lsq_ref_dp))) / real(nchan, dp)
         table%im_fine(m) = sum(sin(phase_k*(real(l_sq, dp)-lsq_ref_dp))) / real(nchan, dp)
      end do
      !$omp end parallel do
   end subroutine build_rmsf_offset_table

   subroutine destroy_rmsf_offset_table(table)
      type(rmsf_table_t), intent(inout) :: table
      if (allocated(table%re_fine)) deallocate(table%re_fine)
      if (allocated(table%im_fine)) deallocate(table%im_fine)
      table%n_fine = 0
   end subroutine destroy_rmsf_offset_table

   subroutine compute_dirty_rmbeam(table, rm_samp, nrm, rm_in, phase_in,&
   &re_beam, im_beam)
      !! Fast path: one table lookup (linear interpolation between the
      !! two nearest fine-grid points -- the table is already
      !! oversampled well past rm_samp's own resolution, so linear is
      !! sufficient; no need for another Fourier/parabolic interpolation
      !! layer on top of an already-fine grid) per output point, plus one
      !! global phase multiply, instead of compute_dirty_rmbeam_direct's
      !! O(nchan) sum per point. See rmsf_table_t's own comment for why
      !! this is exact (not an approximation of the offset dependence
      !! itself) up to the table's own build resolution.
      type(rmsf_table_t), intent(in) :: table
      integer, intent(in) :: nrm
      real(sp), intent(in) :: rm_samp(nrm), rm_in, phase_in
      real(sp), intent(out) :: re_beam(nrm), im_beam(nrm)

      real(dp) :: cos_ph, sin_ph, delta_j, findex, frac, r_re, r_im
      integer :: j, idx

      cos_ph = cos(2.0_dp*real(phase_in, dp))
      sin_ph = sin(2.0_dp*real(phase_in, dp))

      do j = 1, nrm
         delta_j = real(rm_samp(j), dp) - real(rm_in, dp)
         findex = (delta_j - table%delta_min)/table%ddelta + 1.0_dp
         idx = int(findex)
         idx = max(1, min(table%n_fine-1, idx))
         frac = findex - real(idx, dp)
         r_re = table%re_fine(idx)*(1.0_dp-frac) + table%re_fine(idx+1)*frac
         r_im = table%im_fine(idx)*(1.0_dp-frac) + table%im_fine(idx+1)*frac
         re_beam(j) = real(r_re*cos_ph - r_im*sin_ph, sp)
         im_beam(j) = real(r_re*sin_ph + r_im*cos_ph, sp)
      end do
   end subroutine compute_dirty_rmbeam

   subroutine rmsf_point_direct(l_sq, nchan, lsq_ref_compute, delta, re_val, im_val)
      !! Lightweight, single-offset EXACT evaluation of R(delta) --
      !! O(nchan), no template arrays needed (unlike
      !! compute_dirty_rmbeam_direct above, which requires pre-built
      !! cos_arr/sin_arr sized (maxrm,maxchan) and always loops over nrm
      !! output points even for nrm=1). The same closed-form sum
      !! build_rmsf_offset_table's own fill loop already uses per
      !! fine-grid point, factored out here for reuse in a LOCAL search
      !! context (refine_peak_matched_filter below) where only a handful
      !! of offsets are needed per call, not a whole fine table.
      !! planning/RMCLEAN_INTEGRATION_PLAN.md ticket T3.
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan), lsq_ref_compute, delta
      real(sp), intent(out) :: re_val, im_val
      real(dp) :: phase_k

      phase_k = -2.0_dp*real(delta, dp)
      re_val = real(sum(cos(phase_k*(real(l_sq, dp)-real(lsq_ref_compute, dp))))/&
      &real(nchan, dp), sp)
      im_val = real(sum(sin(phase_k*(real(l_sq, dp)-real(lsq_ref_compute, dp))))/&
      &real(nchan, dp), sp)
   end subroutine rmsf_point_direct

   subroutine refine_peak_matched_filter(l_sq, nchan, lsq_ref_compute, rm_samp,&
   &nrm, resid_re, resid_im, imax, drm, noise_rms, nsigma, peak_loc,&
   &re_at_peak, im_at_peak, used_search)
      !! Replaces peak_interp_parabolic's OWN role for the CLEAN peak's
      !! complex value (planning/RMCLEAN_INTEGRATION_PLAN.md ticket T3,
      !! validated in isolation first via tests/test_matched_filter_
      !! refine.f90 before being ported here; TIERED fast-path design
      !! added in T3c, from a design discussion with the user).
      !! peak_interp_parabolic's OWN complex-value step fit a parabola
      !! directly through the raw, carrier-bearing stored Re/Im samples
      !! -- valid only if the OUTER grid already satisfied get_drm's own
      !! max_offset-based bound (~22-44x oversampling relative to fwhm
      !! at lsq_ref=0, thesis P-band numbers). This subroutine removes
      !! that dependency entirely: it never interpolates the stored
      !! samples at all, only evaluates the EXACTLY-known R(delta)
      !! (rmsf_point_direct above, computable from l_sq alone, no raw
      !! per-channel data needed) against the (up to 3) nearest STORED
      !! coarse samples ("anchors" -- the only real measurements in
      !! play; everything else is exact calculation).
      !!
      !! TIERED (T3c): two paths, cheap-first with a data-driven
      !! escalation criterion (the user's own design):
      !!
      !! FAST PATH (runs every call): take the peak LOCATION from the
      !! log-magnitude parabola (peak_interp_parabolic's own location
      !! step -- sound at resolution-level sampling, since magnitude
      !! carries no fast carrier; its own vertex formula guarantees
      !! |offset|<=0.5 whenever the centre bin is the discrete maximum,
      !! which index_absmax guarantees here). At that FIXED location,
      !! solve the one-complex-unknown amplitude A that best explains
      !! the anchors, closed-form (linear least squares / matched
      !! filtering):
      !!   A = sum_j(conj(R(delta_j))*data(j)) / sum_j(|R(delta_j)|^2)
      !! then measure the LEFTOVER misfit, sqrt(sum_j|data_j -
      !! A*R(delta_j)|^2 / nanchor), and compare it against
      !! nsigma*noise_rms (noise_rms: the caller's own data-driven noise
      !! estimate -- clean_complex passes its own per-iteration
      !! rms_about_mean of the residual amplitude). Within threshold:
      !! accept, done -- O(3*nchan) total, no search. This is NOT a
      !! division against a single sample (a rejected earlier design):
      !! a single-anchor division always "fits" perfectly (2 unknowns, 2
      !! equations) and so can never signal that the assumed location
      !! was wrong; fitting >=2 anchors at a FIXED location is
      !! over-determined, and its leftover misfit is a genuine
      !! self-consistency diagnostic -- small means the single-component
      !! model at the parabola's location really does explain the
      !! measurements, large means it doesn't (poor location, blended
      !! components, ...) and the answer should not be trusted.
      !!
      !! FULL SEARCH (escalation, only when the fast path's leftover
      !! misfit exceeds nsigma*noise_rms): try m_search+1 trial
      !! locations spanning [rm_samp(imax)-drm/2, rm_samp(imax)+drm/2]
      !! (the discrete maximum bin can be at most half a cell from the
      !! true peak, else a neighbour would have won), solving the same
      !! closed-form A at each and keeping the trial that maximizes the
      !! matched-filter statistic |sum_j(conj(R)*data)|^2/sum_j(|R|^2)
      !! (equivalently: minimizes the leftover misfit). No closed-form
      !! exists for the best location itself -- the fit statistic is a
      !! sum of per-channel sinusoids at different rates, a function
      !! class with no algebraic solution for its own extrema -- hence a
      !! bounded discrete search, not algebra.
      !!
      !! used_search reports which path produced the answer (fast=.false.)
      !! -- exposed so tests can verify each tier fires when it should,
      !! and so a caller can log escalation statistics.
      !!
      !! Search density (m_search), physically derived (the user's own
      !! requirement -- no arbitrary constants): the search samples the
      !! matched-filter statistic, whose fastest oscillation is TWICE
      !! the carrier rate (the statistic is a squared magnitude --
      !! squaring doubles frequency), carrier rate itself set by
      !! max_offset exactly as in get_drm's own outer-grid formula. So
      !! cycles_in_window = 2 * window_width * max_offset / pi
      !! (window_width=drm), sampled at samples_per_cycle=50 per
      !! statistic cycle -- an empirically validated margin (direct
      !! sweep, tests/test_matched_filter_refine.f90's scenarios: the
      !! failure boundary where the search locks onto a wrong local
      !! maximum sits at ~4 samples per statistic cycle; 50 is a
      !! checked ~12x margin, and the same per-trial spacing validated
      !! for T3's original full-window version, re-expressed for the
      !! halved window). m_floor=20 covers the near-zero-cycles case
      !! (lsq_ref close to max_offset: nothing oscillates, but the
      !! search still needs sub-cell resolution for the location) --
      !! 20 points across one cell, deliberately matching
      !! table_oversample's own default meaning elsewhere in this
      !! project ("how finely to resolve within one native grid cell"),
      !! not an independent invented constant.
      integer, intent(in) :: nchan, nrm, imax
      real(sp), intent(in) :: l_sq(nchan), lsq_ref_compute, rm_samp(nrm)
      real(sp), intent(in) :: resid_re(nrm), resid_im(nrm), drm
      real(sp), intent(in) :: noise_rms, nsigma
      real(sp), intent(out) :: peak_loc, re_at_peak, im_at_peak
      logical, intent(out) :: used_search

      integer, parameter :: samples_per_cycle = 50, m_floor = 20
      real(dp), parameter :: pi_dp = 3.14159265358979_dp
      integer :: j0, j1, nanchor, jj, kk, m_search
      integer :: anchor_idx(3)
      real(sp) :: trial_loc, max_offset, cycles_in_window
      real(sp) :: peak_offset, re_dummy, im_dummy, fast_loc
      real(dp) :: num_re, num_im, den, stat_val, best_stat
      real(dp) :: a_re, a_im, best_a_re, best_a_im, best_loc
      real(dp) :: sum_dd, leftover

      j0 = max(1, imax-1)
      j1 = min(nrm, imax+1)
      nanchor = j1-j0+1
      do jj = 1, nanchor
         anchor_idx(jj) = j0+jj-1
      end do
      sum_dd = 0.0_dp
      do jj = 1, nanchor
         sum_dd = sum_dd + real(resid_re(anchor_idx(jj)), dp)**2 +&
         &real(resid_im(anchor_idx(jj)), dp)**2
      end do

      ! --- FAST PATH: parabola location + closed-form amplitude ---
      if (imax > 1 .and. imax < nrm) then
         call peak_interp_parabolic(resid_re, resid_im, nrm, imax,&
         &peak_offset, re_dummy, im_dummy)
      else
         peak_offset = 0.0_sp
      endif
      fast_loc = rm_samp(imax) + peak_offset*drm

      call eval_trial(fast_loc, num_re, num_im, den, a_re, a_im, stat_val)
      if (den > 0.0_dp) then
         ! leftover misfit = sum|D|^2 - |sum(conj(R)D)|^2/sum|R|^2, the
         ! exact least-squares residual of the fit (never negative up to
         ! rounding -- clamped for the sqrt).
         leftover = max(0.0_dp, sum_dd - stat_val)
         if (real(sqrt(leftover/real(nanchor, dp)), sp) <=&
         &nsigma*max(noise_rms, 0.0_sp)) then
            peak_loc = fast_loc
            re_at_peak = real(a_re, sp)
            im_at_peak = real(a_im, sp)
            used_search = .false.
            return
         endif
      endif

      ! --- FULL SEARCH (escalation) ---
      used_search = .true.
      max_offset = maxval(abs(l_sq - lsq_ref_compute))
      cycles_in_window = real(2.0_dp*real(drm, dp)*real(max_offset, dp)/pi_dp, sp)
      m_search = max(m_floor, ceiling(cycles_in_window*real(samples_per_cycle, sp)))

      best_stat = -1.0_dp
      best_loc = real(rm_samp(imax), dp)
      best_a_re = real(resid_re(imax), dp)
      best_a_im = real(resid_im(imax), dp)

      do kk = 0, m_search
         trial_loc = rm_samp(imax) - 0.5_sp*drm + real(kk, sp)*drm/real(m_search, sp)
         call eval_trial(trial_loc, num_re, num_im, den, a_re, a_im, stat_val)
         if (den > 0.0_dp) then
            if (stat_val > best_stat) then
               best_stat = stat_val
               best_loc = real(trial_loc, dp)
               best_a_re = a_re
               best_a_im = a_im
            endif
         endif
      end do

      peak_loc = real(best_loc, sp)
      re_at_peak = real(best_a_re, sp)
      im_at_peak = real(best_a_im, sp)

   contains

      subroutine eval_trial(loc, num_re_o, num_im_o, den_o, a_re_o, a_im_o, stat_o)
         !! One trial location's closed-form fit against the anchors:
         !! the matched-filter numerator/denominator, the best-fit
         !! complex amplitude, and the fit statistic |num|^2/den (whose
         !! maximization over trials == leftover-misfit minimization).
         !! Shared by the fast path (one call, at the parabola's
         !! location) and the full search (one call per trial) so the
         !! two tiers cannot drift apart numerically.
         real(sp), intent(in) :: loc
         real(dp), intent(out) :: num_re_o, num_im_o, den_o, a_re_o, a_im_o, stat_o
         real(sp) :: re_r_l, im_r_l
         integer :: jj_l

         num_re_o = 0.0_dp
         num_im_o = 0.0_dp
         den_o = 0.0_dp
         do jj_l = 1, nanchor
            call rmsf_point_direct(l_sq, nchan, lsq_ref_compute,&
            &rm_samp(anchor_idx(jj_l))-loc, re_r_l, im_r_l)
            num_re_o = num_re_o + real(re_r_l, dp)*real(resid_re(anchor_idx(jj_l)), dp) +&
            &real(im_r_l, dp)*real(resid_im(anchor_idx(jj_l)), dp)
            num_im_o = num_im_o + real(re_r_l, dp)*real(resid_im(anchor_idx(jj_l)), dp) -&
            &real(im_r_l, dp)*real(resid_re(anchor_idx(jj_l)), dp)
            den_o = den_o + real(re_r_l, dp)**2 + real(im_r_l, dp)**2
         end do
         if (den_o > 0.0_dp) then
            a_re_o = num_re_o/den_o
            a_im_o = num_im_o/den_o
            stat_o = (num_re_o**2+num_im_o**2)/den_o
         else
            a_re_o = 0.0_dp
            a_im_o = 0.0_dp
            stat_o = -1.0_dp
         endif
      end subroutine eval_trial

   end subroutine refine_peak_matched_filter

   subroutine clean_complex(l_sq, nchan, lsq_ref_compute, rm_samp, nrm,&
   &dirty_re, dirty_im, table, niter, gain, thresh, comp_re, comp_im,&
   &resid_re, resid_im, n_iter_used, comp_rm_refined, nsigma_refine)
      !! Hogbom-style complex CLEAN on the dirty FDF (dirty_re/dirty_im),
      !! against its own RM-dependent dirty beam (rmsf_table_t). Ported
      !! from rm_clean.f, modernized, with planning/
      !! RMCLEAN_INTEGRATION_PLAN.md decision 1's restore-order fix
      !! baked into this subroutine's own OUTPUT CONTRACT: comp_re/
      !! comp_im (the pure delta-function clean-component map) and
      !! resid_re/resid_im (the final residual) are returned SEPARATELY,
      !! never summed here -- the original code summed them internally
      !! before its own restore step ever ran, meaning the Gaussian
      !! restoring-beam convolution operated on components+residual
      !! together. The caller (restore_clean) convolves ONLY comp_re/
      !! comp_im with the restoring beam and adds resid_re/resid_im
      !! afterward, unconvolved.
      !!
      !! Two further fixes vs. the original, found by direct derivation
      !! (not just carried over) while porting:
      !!
      !! 1) Component bookkeeping: the original accumulated
      !!    `frac*ResiQ(imax)`/`frac*ResiU(imax)` into its own
      !!    ClnFluxQ/ClnFluxU (frac=gain*peak_val already itself
      !!    proportional to ResiQ/ResiU's own magnitude at the peak --
      !!    dimensionally inconsistent with peak_val*ResiQ(imax) again,
      !!    and unrelated to `frac*re_beam(i)`/`frac*im_beam(i)`, what is
      !!    ACTUALLY subtracted from the residual every iteration).
      !!    Standard Hogbom convention: the complex value recorded for a
      !!    delta-function component is magnitude=frac at the peak's own
      !!    phase, i.e. frac*exp(i*phase_val) -- used here instead.
      !!
      !! 2) Phase convention into compute_dirty_rmbeam: compute_dirty_rmbeam's
      !!    own phase_in is DOUBLED internally (2*(RM_in*L_sq+phase_in),
      !!    matching extract_general_setup's own 2*RM*L_sq convention --
      !!    the standard P(L_sq)=p*exp(2i*chi0)*exp(2i*RM*L_sq) Faraday
      !!    relation, chi0 the UNDOUBLED intrinsic angle). phase_val here
      !!    is atan2(im,re) of the residual itself, i.e. already the
      !!    DOUBLED (2*chi0-convention) angle -- passing it straight into
      !!    compute_dirty_rmbeam as the original code does would double
      !!    it AGAIN (4*chi0 total). Verified directly against this
      !!    module's own compute_dirty_rmbeam (not just derived on paper):
      !!    with phase_in=phase_val/2, the beam's own value AT its exact
      !!    peak (delta=0) is exp(2i*phase_in)=exp(i*phase_val) -- i.e.
      !!    unit magnitude at exactly the residual's own phase, the
      !!    self-consistent choice for subtracting a real-scalar fraction
      !!    of the actual peak value. phase_val/2 is used here, not
      !!    phase_val.
      !!
      !! comp_rm_refined(nrm): per-bin FLUX-WEIGHTED sub-pixel RM location
      !! -- added at the user's own request, root-causing a real chi0-
      !! precision gap. Every iteration already computes a precise
      !! sub-pixel peak_loc (via refine_peak_matched_filter as of ticket
      !! T3 -- see that subroutine's own doc comment; originally
      !! peak_interp_parabolic, superseded because its OWN complex-value
      !! step required the outer grid to satisfy get_drm's demanding
      !! max_offset-based bound, not just resolution-level adequacy) --
      !! but that location information was previously DISCARDED once the
      !! component got filed into its integer grid bin imax, leaving
      !! comp_re/comp_im accurate in amplitude and phase but silently tied
      !! to whatever grid point imax happens to be, which can be biased up
      !! to half a grid cell (dRM/2) away from the true continuous
      !! location -- harmless for reading amplitude/phase off comp_re/
      !! comp_im directly (both already correctly flux-and-phase-weighted
      !! across iterations), but a real problem for anything that needs
      !! the RM value itself precisely, e.g. derotate_to_lsq_ref's own
      !! chi0 = 0.5*phase_val - RM_found*lsq_ref_compute (whose precision
      !! scales directly with RM_found's own). Confirmed empirically:
      !! reading chi0 off a coarse grid's own (perfectly clean, symmetric,
      !! but grid-centred) restored profile gave a ~0.22 rad error at a
      !! band-centroid lsq_ref_compute reference; using this flux-weighted
      !! comp_rm_refined instead, at the SAME coarse grid, gave chi0
      !! exactly matching the true value (0.0000 rad error) -- confirming
      !! the actual bottleneck was never insufficient global grid
      !! resolution or interpolation quality, but this discarded
      !! bookkeeping. Per-bin (not just for the single dominant component)
      !! so this stays correct for multi-component scenarios too (e.g. a
      !! point source AND a separate resolved feature, as in tests/
      !! thesis_scenario_rmclean.f90's own scenario). Bins that never
      !! receive any component flux fall back to their own rm_samp(j)
      !! (moot, since comp_re/comp_im are exactly zero there too, but
      !! keeps this array always well-defined, no NaN/garbage).
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan), lsq_ref_compute
      type(rmsf_table_t), intent(in) :: table
      integer, intent(in) :: nrm, niter
      real(sp), intent(in) :: rm_samp(nrm), dirty_re(nrm), dirty_im(nrm)
      real(sp), intent(in) :: gain, thresh
      real(sp), intent(out) :: comp_re(nrm), comp_im(nrm)
      real(sp), intent(out) :: resid_re(nrm), resid_im(nrm)
      integer, intent(out) :: n_iter_used
      real(sp), intent(out) :: comp_rm_refined(nrm)
      ! nsigma_refine: escalation threshold for refine_peak_matched_
      ! filter's own tiered design (T3c, see that subroutine's doc
      ! comment): the fast fixed-location fit is accepted when its
      ! leftover misfit is within nsigma_refine * (this iteration's own
      ! data-driven noise estimate, rms_val below); beyond that, the
      ! full local search runs instead. Optional; default 3.0 (a
      ! conventional n-sigma consistency cut, confirmed against tests/
      ! test_matched_filter_refine.f90's own noisy scenario rather than
      ! only asserted -- see that test).
      real(sp), intent(in), optional :: nsigma_refine

      real(sp) :: resid_re_snap(nrm), resid_im_snap(nrm)
      real(sp) :: resid_amp(nrm), re_beam(nrm), im_beam(nrm)
      real(sp) :: avg_abs, rms_val, peak_val, peak_loc, phase_val
      real(sp) :: re_at_peak, im_at_peak, frac, dRM, nsig
      logical :: used_search
      ! Accumulator state, kept in double precision internally -- see
      ! this subroutine's own top-of-file precision note: EVERY other
      ! numerically-delicate part of this module (the offset table,
      ! compute_dirty_rmbeam's own interpolation, restore_clean's
      ! convolution) already does its internal arithmetic in double
      ! precision, casting to real(sp) only at the public output
      ! boundary. This loop is the one place that didn't: resid_re/
      ! resid_im get REPEATEDLY subtracted from, up to niter times (often
      ! 1000+) -- unlike the peak-finding step below (which recomputes
      ! fresh from the CURRENT residual snapshot every iteration, so any
      ! single-precision rounding there is a bounded, non-compounding,
      ! per-iteration effect), each subtraction's rounding error
      ! compounds into every subsequent iteration. Confirmed empirically:
      ! a noise-free single-point-source case converged to a residual
      ! power of only ~1e-13 of the dirty map's own power after enough
      ! iterations to fully resolve the source -- suspiciously close to
      ! float32's own ~1.2e-7 relative-amplitude (squared: ~1.4e-14
      ! power) floor, meaning the OLD single-precision accumulator was
      ! the thing limiting how deep CLEAN could converge, not the
      ! algorithm or the data. Real (noisy) data has its own thermal
      ! noise floor almost always far above this, so this rarely mattered
      ! in practice -- but it is a real, identifiable precision
      ! limitation, now removed at negligible cost (a few extra KB per
      ! pixel, not per cube -- CLEAN runs one line of sight at a time).
      real(dp) :: resid_re_dp(nrm), resid_im_dp(nrm)
      real(dp) :: comp_re_dp(nrm), comp_im_dp(nrm)
      real(dp) :: rmloc_wsum(nrm), rmloc_wloc(nrm)
      integer :: iter, imax, j

      comp_re_dp = 0.0_dp
      comp_im_dp = 0.0_dp
      resid_re_dp = real(dirty_re, dp)
      resid_im_dp = real(dirty_im, dp)
      rmloc_wsum = 0.0_dp
      rmloc_wloc = 0.0_dp
      dRM = rm_samp(2) - rm_samp(1)
      nsig = 3.0_sp
      if (present(nsigma_refine)) nsig = nsigma_refine

      do iter = 1, niter
         n_iter_used = iter
         ! Fresh single-precision SNAPSHOT of the current (double-
         ! precision) residual, for peak-finding only -- recomputed from
         ! scratch every iteration, so this doesn't compound (see the
         ! precision note above).
         resid_re_snap = real(resid_re_dp, sp)
         resid_im_snap = real(resid_im_dp, sp)
         resid_amp = sqrt(resid_re_snap**2 + resid_im_snap**2)
         call rms_about_mean(resid_amp, nrm, rms_val)
         call index_absmax(resid_amp, nrm, imax, avg_abs)

         ! T3/T3c: refine_peak_matched_filter replaces peak_interp_
         ! parabolic's own role here (see this subroutine's own
         ! comp_rm_refined doc comment above, and refine_peak_matched_
         ! filter's own doc comment for the tiered fast-path/search
         ! design). It handles domain-edge imax internally (no separate
         ! edge-case branch needed here, unlike the old parabola, which
         ! could not fit through a missing neighbour). rms_val doubles as
         ! the data-driven noise estimate for the tier's own escalation
         ! criterion -- the same per-iteration quantity the stopping
         ! criterion below already uses, not a separate estimator.
         call refine_peak_matched_filter(l_sq, nchan, lsq_ref_compute,&
         &rm_samp, nrm, resid_re_snap, resid_im_snap, imax, dRM,&
         &rms_val, nsig, peak_loc, re_at_peak, im_at_peak, used_search)
         peak_val = sqrt(re_at_peak**2 + im_at_peak**2)
         phase_val = atan2(im_at_peak, re_at_peak)

         if ((peak_val - avg_abs) <= thresh*rms_val) exit

         call compute_dirty_rmbeam(table, rm_samp, nrm, peak_loc,&
         &0.5_sp*phase_val, re_beam, im_beam)

         frac = gain*peak_val
         comp_re_dp(imax) = comp_re_dp(imax) + real(frac*cos(phase_val), dp)
         comp_im_dp(imax) = comp_im_dp(imax) + real(frac*sin(phase_val), dp)
         rmloc_wsum(imax) = rmloc_wsum(imax) + real(frac, dp)
         rmloc_wloc(imax) = rmloc_wloc(imax) + real(frac, dp)*real(peak_loc, dp)

         resid_re_dp = resid_re_dp - real(frac, dp)*real(re_beam, dp)
         resid_im_dp = resid_im_dp - real(frac, dp)*real(im_beam, dp)
      end do

      comp_re = real(comp_re_dp, sp)
      comp_im = real(comp_im_dp, sp)
      resid_re = real(resid_re_dp, sp)
      resid_im = real(resid_im_dp, sp)

      do j = 1, nrm
         if (rmloc_wsum(j) > 0.0_dp) then
            comp_rm_refined(j) = real(rmloc_wloc(j)/rmloc_wsum(j), sp)
         else
            comp_rm_refined(j) = rm_samp(j)
         endif
      end do
   end subroutine clean_complex

   subroutine compute_rmsf_fwhm(l_sq, nchan, fwhm_rm)
      !! Theoretical restoring-beam FWHM in RM, from the lambda-squared
      !! span (edge channels extended by half a channel each, matching
      !! the original rm_restore.f's own edge-handling exactly). pi/
      !! lsq_span (planning/RMCLEAN_INTEGRATION_PLAN.md decision 8 --
      !! the user's deliberate choice, not the Brentjens & de Bruyn
      !! (2005) eq. 61 constant 2*sqrt(3)/lsq_span; not revisited here --
      !! pi and 2*sqrt(3) are close enough that this was never really the
      !! question). Uses this module's own c_velocity (299.792458,
      !! matching rm_synthesis_mod.f90's exactly) rather than the
      !! original code's rounded "300.0" -- a real precision improvement,
      !! not just a style change, since L_sq itself is built with the
      !! precise constant elsewhere in this project.
      !!
      !! Bug fixed (found by the user cross-checking a real ASKAP low-
      !! band restoring-beam value against the ~50 rad/m^2 expected from
      !! the RMSF, not by inspection): this used to compute 0.5*pi/
      !! lsq_span, an extra halving with no basis in this project's own
      !! documented pi/lsq_span intent above, nor in the original
      !! rm_restore.f's own definition of its FWHM_RM argument ("~fac/
      !! Lsq_span... fac is pi for RM and Lambda^2 kind of extraction",
      !! rm_restore.f's own top-of-file comment) -- the 0.5 belonged only
      !! inside THAT file's own separate sigma computation
      !! (sigma=0.5*(0.42466*FWHM_RM)), never in FWHM_RM itself. This
      !! project's own restore_clean applies the standard, un-doubled
      !! 0.42466 FWHM->sigma conversion to whatever this function
      !! returns, so the extra 0.5 here silently halved the actual
      !! restoring beam (and, via the same value, Gate 0's resolution
      !! criterion and the logged "Restoring beam FWHM") to half the true
      !! RMSF resolution. Confirmed via full-codebase search that nothing
      !! else compensates for this: fwhm_rm/compute_rmsf_fwhm(_multiband)
      !! are referenced nowhere outside rmclean.f90/rmclean_cubes.f90
      !! (rm_synthesis.f90/rm_synthesis_mod.f90 have zero occurrences of
      !! "fwhm" at all), and CLEAN's own component-finding
      !! (clean_complex, against compute_dirty_rmbeam's exact dirty RMSF)
      !! and grid-spacing (get_drm, a different bound entirely) never use
      !! this value -- only restore_clean's sigma, Gate 0's
      !! drm_required, and the informational log line do, so this was a
      !! pure restoring-beam-width/reporting bug, never a CLEAN-model
      !! bug.
      !!
      !! Separately, bug fixed vs. the original rm_restore.f (caught by
      !! tests/thesis_scenario_rmclean.f90, not by inspection): the
      !! original's own f1=c/sqrt(L_sq(nchan))/f2=c/sqrt(L_sq(1))
      !! indexing assumes a SPECIFIC l_sq ordering convention (l_sq(1)
      !! the smallest/highest-frequency channel) to end up with a
      !! positive lsq_span. That's the OPPOSITE of L_sq's own documented
      !! convention everywhere else in this project (rm_synthesis.f90's
      !! own comment: "L_sq for ALL channels... in descending lambda_sq
      !! order", i.e. l_sq(1) is the LARGEST/lowest-frequency channel) --
      !! with that convention (the one this module's own
      !! compute_dirty_rmbeam and every test in this project actually
      !! uses), the original's formula silently produces a NEGATIVE
      !! lsq_span (same magnitude, wrong sign), and therefore a negative
      !! FWHM_RM. A "span" is a magnitude, not an ordering-dependent
      !! signed quantity -- abs() here, rather than trying to track or
      !! document a specific required l_sq ordering convention this
      !! subroutine would otherwise silently depend on.
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan)
      real(sp), intent(out) :: fwhm_rm
      real(dp), parameter :: pi = 3.14159265358979_dp

      fwhm_rm = real(pi/padded_lsq_span(l_sq, nchan), sp)
   end subroutine compute_rmsf_fwhm

   subroutine compute_rmsf_fwhm_multiband(l_sq, nchan, band_offset, band_nz,&
   &n_bands, fwhm_rm)
      !! Multi-band restoring-beam FWHM -- NOT the same formula as
      !! compute_rmsf_fwhm above applied to the naively concatenated
      !! l_sq array. Caught empirically by tests/thesis_scenario_rmclean.f90:
      !! feeding a P-band+L-band concatenated array straight into
      !! compute_rmsf_fwhm gives a wildly wrong FWHM, because that
      !! subroutine's own span is (effectively) max(l_sq)-min(l_sq) across
      !! the WHOLE array -- for two widely separated, non-contiguous bands
      !! this is dominated by the large GAP between them, not by either
      !! band's own actual channel coverage. The thesis's own text (Sec
      !! 6.1.4) states the correct rule directly: "the effective
      !! RM-resolution is provided by the SUM of the lambda^2-spans at the
      !! two individual bands" -- verified numerically against Table 6.1
      !! itself (P alone: 0.2007 vs tabulated 0.201; L alone: 0.01255 vs
      !! tabulated 0.013; summed: 0.2133 vs tabulated P+L value 0.214 --
      !! all three matching to within the table's own rounding).
      !!
      !! l_sq(nchan): full concatenated multi-band array. band_offset/
      !! band_nz(n_bands): per-band segmentation, same convention as
      !! rm_synthesis.f90's own band_offset/band_nz (band k occupies
      !! l_sq(band_offset(k)+1 : band_offset(k)+band_nz(k))) -- the
      !! caller already has this from its own multi-band bookkeeping,
      !! not re-derived here by gap-detection heuristics on a flat array.
      integer, intent(in) :: nchan, n_bands
      real(sp), intent(in) :: l_sq(nchan)
      integer, intent(in) :: band_offset(n_bands), band_nz(n_bands)
      real(sp), intent(out) :: fwhm_rm

      real(dp) :: span_sum
      real(dp), parameter :: pi = 3.14159265358979_dp
      integer :: k, i0, i1

      span_sum = 0.0_dp
      do k = 1, n_bands
         i0 = band_offset(k) + 1
         i1 = band_offset(k) + band_nz(k)
         span_sum = span_sum + padded_lsq_span(l_sq(i0:i1), band_nz(k))
      end do
      ! pi/span_sum, not 0.5*pi/span_sum -- see compute_rmsf_fwhm's own
      ! comment for the full story on the erroneous extra 0.5 this used
      ! to carry (same bug, same fix, same rationale).
      fwhm_rm = real(pi/span_sum, sp)
   end subroutine compute_rmsf_fwhm_multiband

   function padded_lsq_span(l_sq, nchan) result(lsq_span)
      !! Shared by compute_rmsf_fwhm and compute_rmsf_fwhm_multiband: one
      !! (single, contiguous) band's own lambda-squared span, edge
      !! channels extended by half a channel each (matching the original
      !! rm_restore.f's own edge-handling), order-independent (abs() --
      !! see compute_rmsf_fwhm's own comment on the sign bug this fixes).
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan)
      real(dp) :: lsq_span
      real(dp) :: f1, f2, dfreq, lsq1, lsq2

      f1 = real(c_velocity, dp)/sqrt(real(l_sq(nchan), dp))
      f2 = real(c_velocity, dp)/sqrt(real(l_sq(1), dp))
      dfreq = (f2-f1)/real(nchan-1, dp)
      f1 = f1 - 0.5_dp*dfreq
      f2 = f2 + 0.5_dp*dfreq
      lsq2 = (real(c_velocity, dp)/f1)**2
      lsq1 = (real(c_velocity, dp)/f2)**2
      lsq_span = abs(lsq2 - lsq1)
   end function padded_lsq_span

   subroutine get_drm(l_sq, nchan, lsq_ref_compute, drm, oversample)
      !! The RM-grid spacing to actually use, GIVEN a chosen
      !! lsq_ref_compute -- a genuine sampling-theorem bound, not an
      !! empirically tuned constant, and (at the user's own explicit
      !! request) the ONLY way to influence dRM at all: there is no raw
      !! "requested dRM" input anywhere in this module's API, deliberately
      !! -- unlike lsq_ref_compute/lsq_ref_report (genuine free choices,
      !! no wrong answer), sampling below this bound is a CORRECTNESS
      !! failure, not a preference, so the only knob exposed is
      !! `oversample`, a multiplier ALWAYS applied on top of the mandatory
      !! floor -- it is architecturally impossible to ask this routine for
      !! an unsafe dRM.
      !!
      !! Distinct from compute_rmsf_fwhm/compute_rmsf_fwhm_multiband
      !! above, which answer a different question entirely:
      !!
      !! The dirty spectrum is P(phi) = (1/K) sum_k p_k
      !! exp(-2i*phi*(l_sq(k)-lsq_ref_compute)). Its MAGNITUDE envelope
      !! |P(phi)| is exactly independent of lsq_ref_compute (a reference
      !! change multiplies the whole sum by exp(-2i*phi*Delta_ref), unit
      !! modulus for every phi) and is set only by how the l_sq(k) are
      !! SPREAD OUT relative to each other -- the SPAN, compute_rmsf_
      !! fwhm's own quantity. This is the physically meaningful
      !! "resolution": how far apart two Faraday components must be to be
      !! distinguished, the analogue of a wave packet's GROUP
      !! velocity/envelope, which carries the actual information content
      !! and cannot be changed by a mere choice of coordinate origin.
      !!
      !! But CLEAN operates on Re(phi)/Im(phi) SEPARATELY (peak_interp_
      !! parabolic fits each independently; clean_complex subtracts a
      !! complex beam every iteration) -- it needs the actual complex
      !! VALUE, not just the magnitude, so it cannot simply discard
      !! lsq_ref_compute's effect. Each term k in the sum is its own
      !! complex sinusoid in phi with period pi/|l_sq(k)-lsq_ref_compute|
      !! -- the analogue of a carrier's PHASE velocity, which conveys no
      !! extra resolution/information by itself (the magnitude envelope
      !! doesn't care about it at all) but which the grid must still
      !! track faithfully if raw Re/Im values are what get sampled, fit,
      !! and subtracted. The fastest such term -- and therefore the
      !! binding constraint -- comes from whichever channel sits FARTHEST
      !! from lsq_ref_compute: max_k|l_sq(k)-lsq_ref_compute|. With
      !! lsq_ref_compute=0 (this project's thesis-matching convention)
      !! and a GHz-range band, that offset equals the channel's own
      !! (large) l_sq directly, far exceeding the (small) SPAN that
      !! governs resolution -- exactly why a grid adequate for
      !! compute_rmsf_fwhm's own resolution scale can still be wildly
      !! under-sampled for CLEAN's own Re/Im grid, a distinct failure mode
      !! from "not enough resolution" (confirmed empirically: CLEAN's
      !! residual diverged, 111 -> 8e16 over ~700 iterations, when this
      !! bound was violated -- see tests/thesis_scenario_rmclean.f90).
      !!
      !! oversample floor, root-caused empirically (not assumed): bare
      !! two-point Nyquist (oversample=1) does NOT diverge -- residual
      !! power stayed small in every scenario tested -- but it recovers
      !! the WRONG answer, not just an imprecise one: a single point
      !! source (RM=50, chi0=0.3 rad) came back as chi0=-0.025 rad at
      !! oversample=1 (tests/test_drm_floor.f90), a genuine error, not
      !! noise. oversample=2 already recovered chi0 to 0.3002 rad
      !! (matching oversample=15's own 0.3000 to within float32 rounding)
      !! -- the transition from "wrong" to "correct" happens somewhere in
      !! (1,2]. This is because peak_interp_parabolic's LOCAL 3-point
      !! quadratic fit needs several samples across a cycle to be a valid
      !! local approximation, not merely the bare two points needed to
      !! avoid ALIASING a global reconstruction (a different, weaker
      !! requirement -- stability alone does not imply correctness here).
      !! Enforced hard floor: oversample must be >= 2 (a margin above the
      !! observed transition point, checked on one scenario only -- not
      !! swept across band shapes/multi-component cases, so treat 2 as a
      !! floor with a safety margin already built in, not a knife-edge
      !! value). Default (if oversample is omitted): 4, giving further
      !! margin above the floor at a fraction of oversample=15's own grid
      !! cost (nrm=80 vs nrm=297 in the tested P-band scenario).
      !!
      !! UPDATE (planning/RMCLEAN_INTEGRATION_PLAN.md ticket T3): the
      !! root cause above -- peak_interp_parabolic's own local parabola
      !! fit through the raw, carrier-bearing samples -- has since been
      !! replaced (refine_peak_matched_filter, used by clean_complex as
      !! of T3). tests/test_drm_floor.f90 now confirms oversample=1 (and
      !! coarser) ALSO recovers correct chi0 through the full CLEAN loop
      !! -- this floor's own ORIGINAL justification (this comment, above)
      !! no longer describes what clean_complex actually needs. The
      !! enforcement below is UNCHANGED for now regardless (Gate 0's own
      !! conservative floor, planning ticket T3b, not yet revisited) --
      !! this is a documentation note about WHY, not a behaviour change.
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan), lsq_ref_compute
      real(sp), intent(out) :: drm
      real(sp), intent(in), optional :: oversample
      real(dp), parameter :: pi = 3.14159265358979_dp
      real(dp) :: max_offset, os

      os = 4.0_dp
      if (present(oversample)) os = real(oversample, dp)
      if (os < 2.0_dp) then
         write(*,'(A)') 'FATAL: get_drm: oversample must be >= 2 -- bare'//&
         &' Nyquist sampling (oversample=1) is confirmed empirically to'//&
         &' recover a WRONG (not just imprecise) chi0/RM; see this'//&
         &' routine''s own doc comment and tests/test_drm_floor.f90.'
         stop 1
      endif

      max_offset = maxval(abs(real(l_sq, dp) - real(lsq_ref_compute, dp)))
      drm = real(pi/(2.0_dp*os*max_offset), sp)
   end subroutine get_drm

   subroutine get_lsq_ref_compute(l_sq, nchan, mode, lsq_ref_compute, fixed_value)
      !! Computes lsq_ref_compute for the requested mode -- see the
      !! lsq_ref_compute_* constants declared at the top of this module
      !! for the full menu. mode=lsq_ref_compute_mid is the RECOMMENDED
      !! DEFAULT (the user's own explicit instruction): pick it unless you
      !! have a specific reason not to. The other modes exist so a caller
      !! can experiment or match some other external convention -- doing
      !! so costs nothing in SAFETY (get_drm enforces its own floor
      !! against whichever lsq_ref_compute value actually results,
      !! regardless of mode) but can cost real GRID SIZE (see below).
      !!
      !! Why mid is not just harmless but actively the smart choice --
      !! both halves of this claim verified, not asserted:
      !!
      !! 1) HARMLESS for accuracy: lsq_ref_compute is a pure bookkeeping
      !!    choice for the CLEAN computation, with NO effect on chi0
      !!    precision once clean_complex's own comp_rm_refined and
      !!    derotate_to_lsq_ref are used correctly (their own comments
      !!    derive and empirically confirm this: chi0 came back IDENTICAL,
      !!    to float32 rounding, whether lsq_ref_compute was 0 or the
      !!    band's own mean, tests/test_rmclean_lsqref_flex.f90). Choosing
      !!    mid therefore trades nothing away.
      !! 2) SMART for cost: get_drm's own bound is set by
      !!    max_k|l_sq(k)-lsq_ref_compute|, the SINGLE worst channel's
      !!    offset. For any fixed set of values, this maximum is
      !!    determined ENTIRELY by the two EXTREME values (min(l_sq),
      !!    max(l_sq)) -- every other channel, however many there are or
      !!    wherever they sit between the extremes, can never exceed
      !!    whichever extreme is farther from lsq_ref_compute, so it never
      !!    enters the max(). Minimizing max(|max(l_sq)-lsq_ref_compute|,
      !!    |lsq_ref_compute-min(l_sq)|) over choice of lsq_ref_compute is
      !!    solved exactly by making the two terms equal -- i.e.
      !!    lsq_ref_compute=(min(l_sq)+max(l_sq))/2, mid's own formula. A
      !!    classic 1D minimax/Chebyshev-centre result, not an
      !!    approximation. Channel COUNT (per-band or overall) does not
      !!    enter this at all -- confirmed numerically (tests/
      !!    test_optimal_lsq_ref.f90) against a P-band(61ch)+L-band(121ch)
      !!    combination, where mid (0.5806) needs a grid nrm=2253, against
      !!    nrm=3127 for the channel-count-weighted mean (0.3771, pulled
      !!    toward L-band's own values by its 2x channel count), nrm=4060
      !!    for centring on the lowest-frequency (P) band's own centroid
      !!    (1.0011 -- barely better than lsq_ref_compute=0's own 4748,
      !!    since doing so leaves the OTHER band almost as exposed as 0
      !!    did), and nrm=4748 for lsq_ref_compute=0 (intrinsic) itself --
      !!    mid is a genuine ~2x-or-more grid-size win here, for zero
      !!    accuracy cost.
      !!
      !! Channel count/SNR DOES matter for a different, unrelated
      !! question -- the actual achievable statistical precision -- but
      !! that is a property of the DATA, decoupled entirely from this
      !! choice (see point 1 above).
      integer, intent(in) :: nchan, mode
      real(sp), intent(in) :: l_sq(nchan)
      real(sp), intent(out) :: lsq_ref_compute
      real(sp), intent(in), optional :: fixed_value

      select case (mode)
      case (lsq_ref_compute_mid)
         lsq_ref_compute = 0.5_sp*(minval(l_sq) + maxval(l_sq))
      case (lsq_ref_compute_intrinsic)
         lsq_ref_compute = 0.0_sp
      case (lsq_ref_compute_centroid)
         lsq_ref_compute = sum(l_sq)/real(nchan, sp)
      case (lsq_ref_compute_min)
         lsq_ref_compute = minval(l_sq)
      case (lsq_ref_compute_max)
         lsq_ref_compute = maxval(l_sq)
      case (lsq_ref_compute_fixed)
         if (.not. present(fixed_value)) then
            write(*,'(A)') 'FATAL: get_lsq_ref_compute: mode=lsq_ref_compute_fixed'//&
            &' requires fixed_value to be supplied.'
            stop 1
         endif
         lsq_ref_compute = fixed_value
      case default
         write(*,'(A)') 'FATAL: get_lsq_ref_compute: unrecognized mode.'
         stop 1
      end select
   end subroutine get_lsq_ref_compute

   subroutine derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute, lsq_ref_report,&
   &re_in, im_in, re_out, im_out)
      !! Re-express a spectrum computed with internal reference
      !! lsq_ref_compute at an ARBITRARY reporting reference
      !! lsq_ref_report -- two DELIBERATELY separate knobs, addressing two
      !! unrelated questions (the user's own explicit ask, after finding
      !! the original lsq_ref-everywhere naming was starting to blur
      !! them): lsq_ref_compute is a pure COMPUTATIONAL COST lever (see
      !! get_drm/get_lsq_ref_compute -- it affects ONLY
      !! grid size, never precision, once comp_rm_refined is used
      !! correctly); lsq_ref_report is a REPORTING CONVENTION choice about
      !! where to QUOTE the intrinsic polarization angle, independent of
      !! how the computation was actually done. lsq_ref_report=0.0
      !! (mode=lsq_ref_report_intrinsic) recovers this project's own
      !! thesis convention (formerly this routine's own hardcoded target,
      !! when it was named derotate_to_lsq_zero); mode=lsq_ref_report_
      !! centroid recovers the Brentjens & de Bruyn (2005)-style
      !! convention, chosen there to minimize the STATISTICAL COVARIANCE
      !! between a jointly-fitted RM and chi0 in the presence of noise (a
      !! classic linear-regression centering result: centering the
      !! regressor at its own weighted mean decorrelates the fitted
      !! intercept and slope) -- a genuinely DIFFERENT optimization
      !! criterion from get_lsq_ref_compute's own mid mode (which
      !! minimizes get_drm's worst-case grid-stability bound instead).
      !! Neither is "more correct" than the other; they answer different
      !! questions, and this routine supports reporting at EITHER (or any
      !! other reference the caller wants) equally well.
      !!
      !! Derivation: the dirty/restored spectrum built at reference
      !! lsq_ref_compute is P_compute(phi) = P_report(phi) *
      !! exp(2i*phi*(lsq_ref_compute-lsq_ref_report)) (get_drm's
      !! own comment derives the general
      !! P_ref'(phi)=P_ref(phi)*exp(-2i*phi*(ref-ref')) relation between
      !! ANY two references; this is that relation with
      !! ref=lsq_ref_compute, ref'=lsq_ref_report). So P_report(phi) =
      !! P_compute(phi) * exp(-2i*phi*(lsq_ref_compute-lsq_ref_report))
      !! recovers the target convention exactly, at every phi
      !! independently -- a pointwise multiply, no re-synthesis needed.
      !! |P_report(phi)|=|P_compute(phi)| exactly (unit-modulus factor),
      !! so this changes ONLY phase, never amplitude -- a directly
      !! verifiable property (tests/test_rmclean_lsqref_flex.f90's own
      !! checks). lsq_ref_report=0.0 reduces exactly to this routine's own
      !! original (derotate_to_lsq_zero) behaviour.
      !!
      !! Practical use for the intrinsic polarization angle: apply this to
      !! clean_complex's own comp_re/comp_im (the pure component map, NOT
      !! the restored/convolved out_re/out_im -- see the precision
      !! guidance below for why), passing clean_complex's own
      !! comp_rm_refined(j) as rm_samp(j) (the actual continuous location,
      !! not the bare grid coordinate) for the bin j of interest, then
      !! read off chi0_report = 0.5*atan2(im_out(j), re_out(j)) -- the
      !! standard P(lsq_ref_report)=p*exp(2i*chi0_report) relation, now
      !! satisfied exactly regardless of which lsq_ref_compute the actual
      !! RM-CLEAN computation used internally. (chi0_report recovered
      !! only mod pi, the ordinary EVPA n*180-degree ambiguity inherent to
      !! any doubled-angle Faraday convention -- not introduced by this
      !! routine.) The restored map remains the right choice for visual
      !! inspection and for integrating an EXTENDED feature's total flux
      !! (tests/thesis_scenario_rmclean.f90's own window_flux) -- it is
      !! only for a single component's own precise RM/phase that
      !! comp_re/comp_im plus comp_rm_refined should be used instead of
      !! the restored map.
      !!
      !! IMPORTANT precision guidance, root-caused empirically (not just
      !! theorised): chi0_report = 0.5*phase_val -
      !! RM_found*(lsq_ref_compute-lsq_ref_report) (this routine's own
      !! derivation). When lsq_ref_compute EQUALS lsq_ref_report, ANY
      !! imprecision in RM_found multiplies against zero and cannot affect
      !! chi0_report at all -- the reason matching the two (e.g.
      !! lsq_ref_compute=lsq_ref_report=0, this routine's own original
      !! special case) is uniquely "free" for phase reporting, not a
      !! coincidence. Whenever the two DIFFER, the SAME RM_found
      !! imprecision gets multiplied by that difference before reaching
      !! chi0_report.
      !!
      !! USE clean_complex's OWN comp_rm_refined(j) for RM_found here, NOT
      !! rm_samp(j) and NOT a fresh peak_interp_parabolic pass over the
      !! restored map. First-principles root-cause (tests/
      !! test_rmclean_lsqref_flex.f90's own nonzero-chi0 case): a coarse,
      !! band-centroid-referenced run's own restored profile looked
      !! perfectly clean and symmetric around its peak, yet was centred
      !! exactly on the GRID POINT rm_samp(imax), ~0.22 rad/m^2 away from
      !! the true continuous RM -- no amount of re-interpolating that
      !! already-quantised output (parabolic or Fourier) can recover
      !! information that was already discarded. clean_complex computes a
      !! precise sub-pixel peak_loc every iteration already (used
      !! correctly for the beam subtraction itself) but used to throw it
      !! away once the component was filed into its integer grid bin --
      !! comp_rm_refined is exactly that discarded information, recovered
      !! (a flux-weighted average of peak_loc across every iteration that
      !! touched a given bin). Confirmed: using comp_rm_refined instead of
      !! rm_samp(imax) at the SAME coarse grid dropped the chi0 error from
      !! ~0.22 rad to 0.0000 rad (exact, to the precision checked) --
      !! confirming the bottleneck was never grid resolution or
      !! interpolation quality, but this specific discarded bookkeeping.
      integer, intent(in) :: nrm
      real(sp), intent(in) :: rm_samp(nrm), lsq_ref_compute, lsq_ref_report
      real(sp), intent(in) :: re_in(nrm), im_in(nrm)
      real(sp), intent(out) :: re_out(nrm), im_out(nrm)

      real(dp) :: ang, c, s, re_dp, im_dp, ref_diff
      integer :: j

      ref_diff = real(lsq_ref_compute, dp) - real(lsq_ref_report, dp)
      do j = 1, nrm
         ang = -2.0_dp*real(rm_samp(j), dp)*ref_diff
         c = cos(ang)
         s = sin(ang)
         re_dp = real(re_in(j), dp)
         im_dp = real(im_in(j), dp)
         re_out(j) = real(re_dp*c - im_dp*s, sp)
         im_out(j) = real(re_dp*s + im_dp*c, sp)
      end do
   end subroutine derotate_to_lsq_ref

   subroutine get_lsq_ref_report(l_sq, nchan, mode, lsq_ref_report, fixed_value)
      !! Computes lsq_ref_report for the requested mode -- see the
      !! lsq_ref_report_* constants declared at the top of this module.
      !! mode=lsq_ref_report_intrinsic (lambda_sq=0) is this project's own
      !! default reporting convention (the user's own thesis). mode=
      !! lsq_ref_report_centroid is the Brentjens & de Bruyn (2005)-style
      !! alternative, a genuinely DIFFERENT optimum from get_lsq_ref_
      !! compute's own mid mode, answering a different question (the
      !! user's own explicit ask to keep these separate, after asking
      !! directly why B&dB advocate the centroid rather than the mean,
      !! and whether that contradicts this module's own mid result -- it
      !! does not; they solve different problems).
      !!
      !! B&dB's own motivation is a classic linear-regression result: when
      !! jointly fitting chi(lambda^2) = chi0 + RM*lambda^2 to noisy
      !! multi-channel data, the estimated intercept (chi0) and slope (RM)
      !! are, in general, statistically CORRELATED -- exactly the same
      !! correlation any linear regression's intercept and slope have
      !! against a regressor not centred at its own mean. Centering
      !! lambda^2 at its own weighted mean before fitting/reporting
      !! removes that correlation entirely (a standard, provable property
      !! of least-squares regression, not specific to RM synthesis). This
      !! is UNRELATED to get_lsq_ref_compute's own mid criterion
      !! (minimizing get_drm's worst-case grid-stability bound, driven
      !! only by the two extreme channels, with no noise or covariance
      !! anywhere in it) -- the two answer different questions and are
      !! not expected to agree.
      !!
      !! mode=lsq_ref_report_centroid's own "weighted" mean -- CORRECTED
      !! after review: this project already has real per-channel
      !! weighting, at the rm_synthesis_mod.f90 level (MASK.CUBE.FITS/
      !! wts_gpu, a genuine per-PIXEL, per-CHANNEL 0/1 mask, confirmed by
      !! reading that module's own source rather than assumed) -- nowhere
      !! in this whole project, checked directly, is any GRADUATED (e.g.
      !! inverse-noise-variance) weight ever used; masking (full
      !! inclusion or full exclusion) is the entire weighting scheme this
      !! project has, and a 0/1 mask reduces exactly to "average over the
      !! surviving channels" -- precisely the plain arithmetic mean this
      !! mode already computes. rmclean_mod itself takes no mask argument
      !! directly (it only ever sees l_sq(nchan), same as every other
      !! routine here) -- masking is therefore the CALLER's
      !! responsibility: build l_sq per pixel from only that pixel's own
      !! unmasked channels (varying nchan pixel-to-pixel is fine) before
      !! calling this routine, exactly as any other rmclean_mod routine
      !! already expects. Given that, this mode is ALREADY correct for
      !! per-pixel masking today, not merely for a future addition. Only
      !! a genuinely graduated (non-0/1) weight, if this project ever
      !! adopts one, would require this branch to become a true weighted
      !! mean, Sum(w_k*l_sq(k))/Sum(w_k) -- flagged here so that addition
      !! doesn't silently leave this routine stale.
      !!
      !! mode=lsq_ref_report_fixed lets a caller request chi0 be quoted at
      !! any OTHER specific lambda_sq value they have their own reason to
      !! want (matching some external paper's own convention, say) --
      !! fixed_value is read only for this mode.
      integer, intent(in) :: nchan, mode
      real(sp), intent(in) :: l_sq(nchan)
      real(sp), intent(out) :: lsq_ref_report
      real(sp), intent(in), optional :: fixed_value

      select case (mode)
      case (lsq_ref_report_intrinsic)
         lsq_ref_report = 0.0_sp
      case (lsq_ref_report_centroid)
         lsq_ref_report = sum(l_sq)/real(nchan, sp)
      case (lsq_ref_report_min)
         lsq_ref_report = minval(l_sq)
      case (lsq_ref_report_max)
         lsq_ref_report = maxval(l_sq)
      case (lsq_ref_report_mid)
         lsq_ref_report = 0.5_sp*(minval(l_sq) + maxval(l_sq))
      case (lsq_ref_report_fixed)
         if (.not. present(fixed_value)) then
            write(*,'(A)') 'FATAL: get_lsq_ref_report: mode=lsq_ref_report_fixed'//&
            &' requires fixed_value to be supplied.'
            stop 1
         endif
         lsq_ref_report = fixed_value
      case default
         write(*,'(A)') 'FATAL: get_lsq_ref_report: unrecognized mode.'
         stop 1
      end select
   end subroutine get_lsq_ref_report

   subroutine restore_clean(rm_samp, nrm, comp_re, comp_im, resid_re,&
   &resid_im, fwhm_rm, plan_fwd, plan_bwd, out_re, out_im)
      !! Restore: convolve ONLY the pure clean-component map with a
      !! Gaussian restoring beam of the given FWHM, then add the
      !! UNCONVOLVED residual -- planning/RMCLEAN_INTEGRATION_PLAN.md
      !! decision 1, the critical fix vs. the original rm_restore.f
      !! (which convolved components+residual already summed, since
      !! rm_clean.f's own output had already merged them before restore
      !! ever ran -- clean_complex above returns them separately for
      !! exactly this reason).
      !!
      !! Implemented as a standard circular convolution via FFTW rather
      !! than porting the original's own fft1d/ifft1d + double-fftshift
      !! sequence (use FFTW per the user's instruction, and that
      !! particular shift/multiply/shift/unshift ordering is not
      !! straightforward to verify correct by inspection): the Gaussian
      !! kernel is built directly in WRAPPED form (index 1 = RM=0, index
      !! nrm/2+2 onward = the "negative lag" tail wrapped to the end of
      !! the array -- the standard layout circular convolution via FFT
      !! requires, equivalent to building the kernel centred in the
      !! array then ifftshift-ing it, just constructed directly). Both
      !! arrays are then FFT'd with NO shifting at all, multiplied,
      !! inverse FFT'd, and divided by nrm (FFTW's transforms are
      !! unnormalized -- see fourier_interp_complex's own comment on the
      !! same point). Verified against a delta-function input (see this
      !! module's own test suite): convolving a unit spike with this
      !! kernel reproduces the kernel itself, correctly centred at the
      !! spike's own location -- the defining property of a convolution
      !! identity element, and the simplest possible correctness check
      !! for a wrapped-kernel FFT convolution.
      !!
      !! plan_fwd/plan_bwd: from plan_fourier_interp(nrm, nrm, ...) --
      !! same size both directions, reused directly rather than a
      !! separate plan-builder for what is mathematically the same
      !! "plan a same-size forward+backward FFTW pair" operation.
      integer(kind=8), intent(in) :: plan_fwd, plan_bwd
      integer, intent(in) :: nrm
      real(sp), intent(in) :: rm_samp(nrm), comp_re(nrm), comp_im(nrm)
      real(sp), intent(in) :: resid_re(nrm), resid_im(nrm), fwhm_rm
      real(sp), intent(out) :: out_re(nrm), out_im(nrm)

      real(dp) :: dRM, sigma, two_sigma_sq, dist
      complex(dp) :: comp_c(nrm), kernel_c(nrm)
      integer :: i, k, half

      dRM = real(rm_samp(2) - rm_samp(1), dp)
      sigma = 0.42466_dp*real(fwhm_rm, dp)
      two_sigma_sq = 2.0_dp*sigma*sigma

      comp_c = cmplx(real(comp_re, dp), real(comp_im, dp), dp)

      half = nrm/2
      do i = 1, nrm
         k = i - 1
         if (k > half) k = k - nrm
         dist = real(k, dp)*dRM
         kernel_c(i) = cmplx(exp(-dist*dist/two_sigma_sq), 0.0_dp, dp)
      end do

      call dfftw_execute_dft(plan_fwd, comp_c, comp_c)
      call dfftw_execute_dft(plan_fwd, kernel_c, kernel_c)
      comp_c = comp_c * kernel_c
      call dfftw_execute_dft(plan_bwd, comp_c, comp_c)
      comp_c = comp_c / real(nrm, dp)

      out_re = real(real(comp_c, dp), sp) + resid_re
      out_im = real(aimag(comp_c), sp) + resid_im
   end subroutine restore_clean

end module rmclean_mod
