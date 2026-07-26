module rmclean_mod
   !! RM-CLEAN core: Hogbom-style deconvolution of a complex Faraday
   !! dispersion function (dirty FDF) against its own RM-dependent dirty
   !! beam (RMSF). Pure computation only -- no FITS I/O -- mirroring
   !! gaussft_mod's own split between narrowly-scoped computation and a
   !! caller that owns I/O/config. Ported and modernized from the user's
   !! own thesis Fortran77 code (~/softwares/CURR_DEVEL/RM_CLEAN_TESTS/),
   !! not a raw port -- see planning/RMCLEAN_INTEGRATION_PLAN.md for the
   !! full design record (decisions confirmed with the user) this module
   !! implements.
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   implicit none
   private
   public :: index_absmax, peak_interp_parabolic
   public :: plan_fourier_interp, destroy_fourier_interp_plan, fourier_interp_complex
   public :: rmsf_table_t, build_rmsf_offset_table, destroy_rmsf_offset_table
   public :: compute_dirty_rmbeam_direct, compute_dirty_rmbeam
   public :: clean_complex
   public :: compute_rmsf_fwhm, compute_rmsf_fwhm_multiband, restore_clean
   public :: required_drm_nyquist

   ! FFTW3 constants (same convention as gaussft_mod's own -- declared
   ! directly rather than `include`d, fftw3.f is fixed-form F77 and
   ! cannot be included into a free-form .f90 file directly).
   integer, parameter :: fftw_forward = -1
   integer, parameter :: fftw_backward = 1
   integer, parameter :: fftw_estimate = 64

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
      !!   R(delta) = (1/nchan) * sum_k exp(-2i*delta*(L_sq(k)-lsq_ref))
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
   &cos_arr, sin_arr, nrm, maxrm, maxchan, lsq_ref, re_beam, im_beam)
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
      !! lsq_ref: the phase-reference lambda-squared, EXPLICITLY supplied
      !! by the caller rather than computed internally as mean(l_sq) --
      !! this is a pure convention choice (the derivation in
      !! rmsf_table_t's own comment holds for ANY fixed reference point,
      !! not specifically the mean), and different callers legitimately
      !! want different choices: rm_synthesis_mod.f90's own
      !! extract_general_setup references to the mean (numerically
      !! smaller phase arguments across a whole run); the user's own
      !! thesis codebase references to lambda_sq=0 (no subtraction at
      !! all). Whatever the caller passes here MUST match whatever
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
      real(sp), intent(in) :: lsq_ref
      real(sp), intent(out) :: re_beam(nrm), im_beam(nrm)

      real(sp) :: phi_tmp
      real(sp) :: ryt(nchan), iyt(nchan)
      real(sp) :: rc_cor, rs_cor, ic_cor, is_cor
      integer :: j

      block
         integer :: kk
         do kk = 1, nchan
            phi_tmp = 2.0_sp*(rm_in*(l_sq(kk)-lsq_ref) + phase_in)
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

   subroutine build_rmsf_offset_table(l_sq, nchan, lsq_ref, delta_span,&
   &native_ddelta, oversample, table)
      !! Build R(delta) once (see rmsf_table_t's own comment for the
      !! derivation) on a grid finer than the dirty map's own RM sampling
      !! by `oversample`, spanning [-delta_span,+delta_span] -- the
      !! caller passes delta_span = rm_samp(nrm)-rm_samp(1), the largest
      !! offset a (rm_in, rm_samp(j)) pair inside the search domain can
      !! ever produce, so every lookup this table will ever see falls
      !! inside its built range.
      !!
      !! lsq_ref: MUST match whatever reference point the caller's own
      !! dirty-map/injected-data construction uses -- see
      !! compute_dirty_rmbeam_direct's own comment on why this is a real
      !! correctness requirement, not a cosmetic choice: R(delta) itself
      !! is only a pure function of delta for a FIXED, shared reference
      !! point; changing the reference multiplies R(delta) by an extra
      !! delta-dependent phase factor, exp(-2i*delta*(new_ref-old_ref)),
      !! not just a constant rotation.
      real(sp), intent(in) :: l_sq(nchan)
      integer, intent(in) :: nchan
      real(sp), intent(in) :: lsq_ref
      real(sp), intent(in) :: delta_span, native_ddelta
      integer, intent(in) :: oversample
      type(rmsf_table_t), intent(out) :: table

      real(dp) :: lsq_ref_dp, delta_m, phase_k
      integer :: m

      lsq_ref_dp = real(lsq_ref, dp)

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

   subroutine clean_complex(rm_samp, nrm, dirty_re, dirty_im, table, niter,&
   &gain, thresh, comp_re, comp_im, resid_re, resid_im, n_iter_used)
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
      type(rmsf_table_t), intent(in) :: table
      integer, intent(in) :: nrm, niter
      real(sp), intent(in) :: rm_samp(nrm), dirty_re(nrm), dirty_im(nrm)
      real(sp), intent(in) :: gain, thresh
      real(sp), intent(out) :: comp_re(nrm), comp_im(nrm)
      real(sp), intent(out) :: resid_re(nrm), resid_im(nrm)
      integer, intent(out) :: n_iter_used

      real(sp) :: resid_amp(nrm), re_beam(nrm), im_beam(nrm)
      real(sp) :: avg_abs, rms_val, peak_val, peak_loc, phase_val
      real(sp) :: peak_offset, re_at_peak, im_at_peak, frac, dRM
      integer :: iter, imax

      comp_re = 0.0_sp
      comp_im = 0.0_sp
      resid_re = dirty_re
      resid_im = dirty_im
      dRM = rm_samp(2) - rm_samp(1)

      do iter = 1, niter
         n_iter_used = iter
         resid_amp = sqrt(resid_re**2 + resid_im**2)
         call rms_about_mean(resid_amp, nrm, rms_val)
         call index_absmax(resid_amp, nrm, imax, avg_abs)

         if (imax > 1 .and. imax < nrm) then
            call peak_interp_parabolic(resid_re, resid_im, nrm, imax,&
            &peak_offset, re_at_peak, im_at_peak)
         else
            ! Sampled peak is on the domain edge -- no neighbour on one
            ! side to fit a parabola through; the sampled bin is already
            ! the best available estimate (same edge policy as the
            ! original quad_interp.f).
            peak_offset = 0.0_sp
            re_at_peak = resid_re(imax)
            im_at_peak = resid_im(imax)
         endif
         peak_loc = rm_samp(imax) + peak_offset*dRM
         peak_val = sqrt(re_at_peak**2 + im_at_peak**2)
         phase_val = atan2(im_at_peak, re_at_peak)

         if ((peak_val - avg_abs) <= thresh*rms_val) exit

         call compute_dirty_rmbeam(table, rm_samp, nrm, peak_loc,&
         &0.5_sp*phase_val, re_beam, im_beam)

         frac = gain*peak_val
         comp_re(imax) = comp_re(imax) + frac*cos(phase_val)
         comp_im(imax) = comp_im(imax) + frac*sin(phase_val)

         resid_re = resid_re - frac*re_beam
         resid_im = resid_im - frac*im_beam
      end do
   end subroutine clean_complex

   subroutine compute_rmsf_fwhm(l_sq, nchan, fwhm_rm)
      !! Theoretical restoring-beam FWHM in RM, from the lambda-squared
      !! span (edge channels extended by half a channel each, matching
      !! the original rm_restore.f's own edge-handling exactly). Kept as
      !! pi/lsq_span (planning/RMCLEAN_INTEGRATION_PLAN.md decision 8 --
      !! the user's deliberate choice, not the Brentjens & de Bruyn
      !! (2005) eq. 61 constant 2*sqrt(3)/lsq_span; not revisited here).
      !! Uses this module's own c_velocity (299.792458, matching
      !! rm_synthesis_mod.f90's exactly) rather than the original code's
      !! rounded "300.0" -- a real precision improvement, not just a
      !! style change, since L_sq itself is built with the precise
      !! constant elsewhere in this project.
      !!
      !! Bug fixed vs. the original rm_restore.f (caught by
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

      fwhm_rm = real(0.5_dp*pi/padded_lsq_span(l_sq, nchan), sp)
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
      fwhm_rm = real(0.5_dp*pi/span_sum, sp)
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

   subroutine required_drm_nyquist(l_sq, nchan, lsq_ref, oversample, drm_max)
      !! The RM-grid spacing an off-grid-peak CLEAN needs, GIVEN a chosen
      !! lsq_ref -- a genuine sampling-theorem bound, not an empirically
      !! tuned constant. Distinct from compute_rmsf_fwhm/
      !! compute_rmsf_fwhm_multiband above, which answer a different
      !! question entirely:
      !!
      !! The dirty spectrum is P(phi) = (1/K) sum_k p_k
      !! exp(-2i*phi*(l_sq(k)-lsq_ref)). Its MAGNITUDE envelope |P(phi)|
      !! is exactly independent of lsq_ref (a reference change multiplies
      !! the whole sum by exp(-2i*phi*Delta_ref), unit modulus for every
      !! phi) and is set only by how the l_sq(k) are SPREAD OUT relative
      !! to each other -- the SPAN, compute_rmsf_fwhm's own quantity. This
      !! is the physically meaningful "resolution": how far apart two
      !! Faraday components must be to be distinguished, the analogue of
      !! a wave packet's GROUP velocity/envelope, which carries the
      !! actual information content and cannot be changed by a mere
      !! choice of coordinate origin.
      !!
      !! But CLEAN operates on Re(phi)/Im(phi) SEPARATELY (peak_interp_
      !! parabolic fits each independently; clean_complex subtracts a
      !! complex beam every iteration) -- it needs the actual complex
      !! VALUE, not just the magnitude, so it cannot simply discard
      !! lsq_ref's effect. Each term k in the sum is its own complex
      !! sinusoid in phi with period pi/|l_sq(k)-lsq_ref| -- the
      !! analogue of a carrier's PHASE velocity, which conveys no extra
      !! resolution/information by itself (the magnitude envelope doesn't
      !! care about it at all) but which the grid must still track
      !! faithfully if raw Re/Im values are what get sampled, fit, and
      !! subtracted. The fastest such term -- and therefore the binding
      !! constraint -- comes from whichever channel sits FARTHEST from
      !! lsq_ref: max_k|l_sq(k)-lsq_ref|. With lsq_ref=0 (this project's
      !! thesis-matching convention) and a GHz-range band, that offset
      !! equals the channel's own (large) l_sq directly, far exceeding
      !! the (small) SPAN that governs resolution -- exactly why a grid
      !! adequate for compute_rmsf_fwhm's own resolution scale can still
      !! be wildly under-sampled for CLEAN's own Re/Im grid, a distinct
      !! failure mode from "not enough resolution" (confirmed empirically:
      !! CLEAN's residual diverged, 111 -> 8e16 over ~700 iterations, when
      !! this bound was violated -- see tests/thesis_scenario_rmclean.f90).
      !!
      !! Bare two-point Nyquist for that fastest term is
      !! drm <= pi/(2*max_offset). oversample (>1) tightens this by the
      !! given factor, since peak_interp_parabolic's LOCAL 3-point
      !! quadratic fit needs several samples across a cycle to be a valid
      !! local approximation, not merely the bare two points needed to
      !! avoid ALIASING a global reconstruction (a different, weaker,
      !! requirement) -- an oversample of ~10-15 was empirically confirmed
      !! sufficient (CLEAN's convergence trajectory matched the
      !! lsq_ref=mean baseline closely at oversample~14) and is a
      !! documented safety margin above the bare bound, not a second
      !! hidden ad-hoc constant.
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan), lsq_ref, oversample
      real(sp), intent(out) :: drm_max
      real(dp), parameter :: pi = 3.14159265358979_dp
      real(dp) :: max_offset

      max_offset = maxval(abs(real(l_sq, dp) - real(lsq_ref, dp)))
      drm_max = real(pi/(2.0_dp*real(oversample, dp)*max_offset), sp)
   end subroutine required_drm_nyquist

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
