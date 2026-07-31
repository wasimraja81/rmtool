program thesis_scenario_rmclean
   !! Reproduces Raja (2014) "Faraday Slicing Polarized Radio Sources"
   !! Chapter 6 Figures 6.1/6.2/6.3 (L-alone/P-alone/P+L-combined RM-CLEAN
   !! of a simulated line of sight), Table 6.1/6.2 exact, as a single
   !! line-of-sight test of rmclean_mod (planning/
   !! RMCLEAN_INTEGRATION_PLAN.md T1) -- no cube/FITS I/O needed, per the
   !! user's own suggestion ("single line of sight instead of a cube").
   !! Reuses the SAME P/L band channel counts as tests/
   !! make_thesis_scenario_cubes.py (61/121 channels) for consistency
   !! with that already-established thesis-scenario fixture, but does
   !! NOT reuse its F2/F3 addition (that's this project's own extension
   !! beyond the thesis's own demo, not part of Figures 6.1-6.3) and
   !! adds no noise (the thesis's own simulated spectra are themselves
   !! noise-free -- see the smooth, unscattered dirty profiles in
   !! Figures 6.1-6.3).
   !!
   !! Sky model (Table 6.2, exact): a Faraday-thin point source
   !! (RM=-100 rad/m^2, amplitude 15 Jy/(rad/m^2), PA=0) plus a
   !! Faraday-thick top-hat (RM 100-130 rad/m^2, total amplitude
   !! 5 Jy/(rad/m^2), PA=0) -- approximated as 61 closely-spaced thin
   !! slices spanning the top-hat (0.5 rad/m^2 spacing, far below any
   !! delta-RM achievable in this scenario), same discretisation
   !! convention as make_thesis_scenario_cubes.py's own.
   !!
   !! CLEAN parameters niter=1000/loop-gain=0.1 taken directly from the
   !! thesis figures' own panel labels. thresh=0.0 (not shown in the
   !! panels, only niter/loop-gain are) -- deliberately near-disables the
   !! peak-vs-mean early-stop test, letting niter alone govern
   !! convergence. Needed because that test's own "mean(|residual|)" term
   !! is a poor noise-floor proxy for a genuinely Faraday-thick (broad,
   !! non-spiky) profile: for the P+L case specifically, a moderate
   !! thresh (checked: 3.0) made the loop exit after a single iteration,
   !! since the broad hump's own mean amplitude sits close to its peak
   !! from the very first iteration -- confirmed a real early-termination
   !! effect (not a CLEAN correctness bug: clean_complex's own math is
   !! the same faithfully-ported formula either way) by disabling the
   !! stopping test and watching the manual iteration converge steadily.
   !! With thresh=0.0, P-alone/L-alone still run close to the full 1000
   !! (visually unchanged from a moderate thresh -- once the real
   !! structure is extracted, extra iterations at gain=0.1 only chip at
   !! already-small residual ripple), while P+L now converges naturally
   !! partway through (stops at 50) instead of after 1.
   !!
   !! Expected qualitative outcome (thesis text, Sec 6.1.2-6.1.4):
   !! L-band alone -- poor RM-resolution (delta-RM ~250, far coarser
   !! than the 30 rad/m^2 top-hat), extended feature not resolved.
   !! P-band alone -- fine RM-resolution (delta-RM ~15.6) but the
   !! extended feature is not usefully detected regardless (thesis's own
   !! wording: "cleaned tomograph shows no significant trace of the
   !! extended feature"), because P-band's own max RM-scale (~3.6) is
   !! far too small to properly represent structure spread this widely;
   !! the point source is recovered reliably. P+L combined -- "a
   !! significant fraction of the flux density for the extended feature"
   !! is recovered, and the point source remains reliable.
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp  ! Mm/s
   real(sp), parameter :: point_rm = -100.0_sp, point_amp = 15.0_sp
   real(sp), parameter :: thick_lo = 100.0_sp, thick_hi = 130.0_sp
   real(sp), parameter :: thick_amp_total = 5.0_sp
   ! thick_flux_true: the top-hat's own TOTAL integrated flux (height *
   ! width, Jy) -- what a real "total flux of an extended component"
   ! measurement should be compared against, analogous to integrating an
   ! extended source's surface brightness over its own solid angle in
   ! ordinary imaging. The point source's own "total" is just its
   ! amplitude directly (peak=total for an unresolved/delta-function
   ! component, same convention as measuring an unresolved source's flux
   ! straight off its peak pixel).
   real(sp), parameter :: thick_flux_true = thick_amp_total*(thick_hi-thick_lo)
   integer, parameter :: thick_n_slices = 61
   real(sp), parameter :: chi0 = 0.0_sp

   ! RM-grid sampling: with lsq_ref=0 (below), the dirty map's Re/Im
   ! components (not just |amplitude|) oscillate at a rate set by the
   ! CHANNELS' OWN ABSOLUTE lambda^2 (~1.1 for P-band), not the small
   ! spread around a mean (~0.2) -- a genuine Nyquist-type sampling
   ! requirement for peak_interp_parabolic's own Re/Im-separate quadratic
   ! fit to be valid at all (fitting a quadratic through 3 samples of a
   ! function that completes a full oscillation cycle within a handful of
   ! samples is meaningless, not just imprecise). dRM is therefore DERIVED
   ! below via rmclean_mod's own get_drm (not a hardcoded
   ! literal): confirmed directly that a coarser grid violating this bound
   ! (dRM=2, fine enough only for lsq_ref=mean's much slower oscillation)
   ! made CLEAN's own residual DIVERGE under lsq_ref=0 -- peak growing from
   ! 111 to 8e16 over 700 iterations before hitting floating-point
   ! overflow -- while a grid satisfying this bound converges cleanly
   ! (peak monotonically falling from 0.155 to 5.3e-4 over 1000
   ! iterations, matching the lsq_ref=mean case's own convergence
   ! trajectory closely). rm_drm_oversample=15 was empirically confirmed
   ! sufficient by that same convergence check -- comfortably above
   ! get_drm's own ENFORCED floor of 2 (not the bare 2-point Nyquist
   ! value of 1: tests/test_drm_floor.f90 found oversample=1 recovers a
   ! WRONG chi0, not just an imprecise one, since peak_interp_parabolic's
   ! local quadratic fit needs several samples per cycle to be valid, not
   ! just the two points that merely avoid aliasing a global
   ! reconstruction).
   real(sp), parameter :: rm_lo = -200.0_sp, rm_span = 450.0_sp
   real(sp), parameter :: rm_drm_oversample = 15.0_sp
   integer :: nrm
   real(sp) :: drm
   real(sp), allocatable :: rm_samp(:)
   real(sp) :: l_sq_p_worst(61), l_sq_l_worst(121), l_sq_worst(61+121)

   logical :: all_pass
   character(len=256) :: outdir

   if (command_argument_count() >= 1) then
      call get_command_argument(1, outdir)
   else
      outdir = '.'
   endif

   all_pass = .true.

   ! Derive the RM-grid spacing from the ACTUAL channel set that will be
   ! used (P-band's own l_sq dominates, being farthest from lsq_ref=0 --
   ! see get_drm's own doc comment -- but L-band is folded in
   ! too so this stays correct even if the band setup changes later),
   ! rather than hardcoding a literal value tuned to just this scenario.
   call band_channels(300.0_sp, 30.0_sp, 61, l_sq_p_worst)
   call band_channels(1200.0_sp, 120.0_sp, 121, l_sq_l_worst)
   l_sq_worst(1:61) = l_sq_p_worst
   l_sq_worst(62:182) = l_sq_l_worst
   call get_drm(l_sq_worst, 182, 0.0_sp, drm, oversample=rm_drm_oversample)
   nrm = nint(rm_span/drm) + 1
   allocate(rm_samp(nrm))
   call build_rm_axis(rm_samp, nrm, rm_lo, drm)
   write(*,'(A,F0.1,A,F0.4,A,I0,A)') 'Derived RM-grid spacing (Nyquist-type, oversample=',&
   &rm_drm_oversample, '): dRM=', drm, ' rad/m^2 (nrm=', nrm, ')'

   write(*,'(A)') '===================================================='
   write(*,'(A)') 'Input sky model (Table 6.2, exact):'
   write(*,'(A,F0.1,A,F0.1,A,F0.1,A)') '  Point source:  RM=', point_rm,&
   &' rad/m^2, amplitude=', point_amp, ' Jy/(rad/m^2), PA=0 (total flux=',&
   &point_amp, ' Jy, unresolved: peak=total)'
   write(*,'(A,F0.1,A,F0.1,A,F0.1,A,F0.1,A,F0.1,A,F0.1,A)') '  Faraday-thick top-hat: RM=[',&
   &thick_lo, ',', thick_hi, '] rad/m^2, height=', thick_amp_total,&
   &' Jy/(rad/m^2) (total flux = height * width = ', thick_amp_total, ' * ',&
   &thick_hi-thick_lo, ' = ', thick_flux_true, ' Jy)'
   write(*,'(A)') '===================================================='

   call run_band_case('L-alone (Fig 6-1)', 'fig6_1_Lalone', 1200.0_sp, 120.0_sp,&
   &121, rm_samp, nrm, trim(outdir), all_pass)
   call run_band_case('P-alone (Fig 6-2)', 'fig6_2_Palone', 300.0_sp, 30.0_sp,&
   &61, rm_samp, nrm, trim(outdir), all_pass)
   call run_combined_case(rm_samp, nrm, trim(outdir), all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] thesis_scenario_rmclean: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] thesis_scenario_rmclean: one or more checks failed'
      stop 1
   endif

contains

   subroutine build_rm_axis(rm, n, lo, d)
      integer, intent(in) :: n
      real(sp), intent(out) :: rm(n)
      real(sp), intent(in) :: lo, d
      integer :: j
      do j = 1, n
         rm(j) = lo + real(j-1, sp)*d
      end do
   end subroutine build_rm_axis

   subroutine sky_model_qu(l_sq, nchan, q, u)
      !! Table 6.2 sky model: linear superposition of the point source
      !! and the (finely-sliced) Faraday-thick top-hat, both PA=0.
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan)
      real(sp), intent(out) :: q(nchan), u(nchan)
      real(sp) :: slice_amp, slice_rm, phase
      integer :: k, s

      q = 0.0_sp
      u = 0.0_sp
      do k = 1, nchan
         phase = 2.0_sp*(point_rm*l_sq(k) + chi0)
         q(k) = q(k) + point_amp*cos(phase)
         u(k) = u(k) + point_amp*sin(phase)
      end do

      ! thick_amp_total is a DENSITY (Table 6.2's own units: Jy per unit
      ! RM), the top-hat's own HEIGHT -- not a total flux to divide
      ! across slices. Confirmed against the thesis figures themselves
      ! (Chapter 6, 3rd panel "Input Amplitude"): the top-hat is drawn
      ! reaching height 5 on the SAME axis the point source reaches 15,
      ! not some much smaller divided-by-width value. Each slice
      ! therefore carries amplitude = density * slice_width, giving
      ! total integrated flux = thick_amp_total * (thick_hi-thick_lo) =
      ! 150, not 5. (tests/make_thesis_scenario_cubes.py, this project's
      ! existing RM-synthesis-focused thesis fixture, divides by
      ! thick_n_slices instead -- i.e. treats "5" as already the total --
      ! flagged as a likely discrepancy there too, but not changed here:
      ! that script's own checks are relative/qualitative and unaffected
      ! by a uniform rescaling, and it's out of scope for this ticket.)
      slice_amp = thick_amp_total*(thick_hi-thick_lo)/real(thick_n_slices-1, sp)
      do s = 1, thick_n_slices
         slice_rm = thick_lo + (thick_hi-thick_lo)*real(s-1, sp)/real(thick_n_slices-1, sp)
         do k = 1, nchan
            phase = 2.0_sp*(slice_rm*l_sq(k) + chi0)
            q(k) = q(k) + slice_amp*cos(phase)
            u(k) = u(k) + slice_amp*sin(phase)
         end do
      end do
   end subroutine sky_model_qu

   subroutine band_channels(nu_c, dnu, nchan, l_sq)
      !! Table 6.1: nu_c/dnu define the band; channel centres uniformly
      !! spaced across it (same convention as make_thesis_scenario_cubes.py).
      real(sp), intent(in) :: nu_c, dnu
      integer, intent(in) :: nchan
      real(sp), intent(out) :: l_sq(nchan)
      real(sp) :: f_step, f_start, freq
      integer :: k
      f_step = dnu/real(nchan, sp)
      f_start = nu_c - 0.5_sp*dnu + 0.5_sp*f_step
      do k = 1, nchan
         freq = f_start + real(k-1, sp)*f_step
         l_sq(k) = (c_light/freq)**2
      end do
   end subroutine band_channels

   subroutine dirty_spectrum(l_sq, nchan, q, u, rm, nrm_l, dirty_re, dirty_im)
      !! Direct DFT dirty tomograph. Phase-referenced to lambda_sq=0 (NO
      !! subtraction of any mean/reference lambda^2), per the user's own
      !! explicit instruction, matching their own thesis codebase's
      !! convention exactly -- and ALSO matching rm_synthesis_mod.f90's
      !! own extract_general_setup convention, which is itself
      !! unconditionally at lambda_sq=0 (confirmed directly at
      !! rm_synthesis_mod.f90:675-791: its own phi_tmp=omega*t(kk) uses
      !! raw L_sq, no mean subtraction anywhere) -- this test's own
      !! convention is not a deliberate divergence from rm_synthesis's,
      !! it is the same convention. Must stay consistent
      !! with sky_model_qu's own convention (also raw l_sq, no
      !! subtraction) and with whatever lsq_ref is passed into
      !! build_rmsf_offset_table/compute_dirty_rmbeam_direct below (0.0
      !! throughout this file) -- a real correctness requirement, not a
      !! cosmetic choice: R(delta) itself is only a pure function of
      !! delta for a single, FIXED, shared reference point (see
      !! rmclean_mod's own compute_dirty_rmbeam_direct/
      !! build_rmsf_offset_table comments, which this test's own earlier
      !! empirical check confirmed directly -- mismatching references
      !! between the injected sky model and the dirty-map template gave
      !! wrong RELATIVE amplitudes between the point source and the
      !! Faraday-thick component, not just a benign rotation).
      integer, intent(in) :: nchan, nrm_l
      real(sp), intent(in) :: l_sq(nchan), q(nchan), u(nchan), rm(nrm_l)
      real(sp), intent(out) :: dirty_re(nrm_l), dirty_im(nrm_l)
      real(sp) :: c_tmpl(nchan), s_tmpl(nchan)
      integer :: j
      do j = 1, nrm_l
         c_tmpl = cos(2.0_sp*rm(j)*l_sq)
         s_tmpl = -sin(2.0_sp*rm(j)*l_sq)
         dirty_re(j) = (dot_product(q, c_tmpl) - dot_product(u, s_tmpl))/real(nchan, sp)
         dirty_im(j) = (dot_product(q, s_tmpl) + dot_product(u, c_tmpl))/real(nchan, sp)
      end do
   end subroutine dirty_spectrum

   subroutine input_model_amp(rm, n, amp)
      !! The TRUE sky model's own amplitude density vs RM, for the
      !! "Input Amplitude" panel -- drawn exactly as the thesis figures'
      !! own 3rd panel does: the point source as a single-bin spike of
      !! height point_amp, the top-hat as a flat plateau of height
      !! thick_amp_total across [thick_lo,thick_hi] (the TRUE model's
      !! own height, independent of how sky_model_qu discretises it into
      !! slices for the DFT -- those are a computational approximation,
      !! this is the thing being approximated).
      integer, intent(in) :: n
      real(sp), intent(in) :: rm(n)
      real(sp), intent(out) :: amp(n)
      integer :: j, ispike
      amp = 0.0_sp
      ispike = 1
      do j = 2, n
         if (abs(rm(j)-point_rm) < abs(rm(ispike)-point_rm)) ispike = j
      end do
      amp(ispike) = point_amp
      do j = 1, n
         if (rm(j) >= thick_lo .and. rm(j) <= thick_hi) amp(j) = amp(j) + thick_amp_total
      end do
   end subroutine input_model_amp

   subroutine write_profile_csv(outdir, slug, rm, n, dirty_re, dirty_im,&
   &clean_re, clean_im, input_amp)
      character(len=*), intent(in) :: outdir, slug
      integer, intent(in) :: n
      real(sp), intent(in) :: rm(n), dirty_re(n), dirty_im(n)
      real(sp), intent(in) :: clean_re(n), clean_im(n), input_amp(n)
      integer :: unit, j
      real(sp) :: dirty_amp, dirty_pha_deg, clean_amp, clean_pha_deg
      real(sp), parameter :: rad2deg = 180.0_sp/3.14159265358979_sp

      unit = 60
      open(unit, file=trim(outdir)//'/'//trim(slug)//'_profile.csv',&
      &status='unknown', action='write')
      write(unit,'(A)') 'RM,dirty_re,dirty_im,dirty_amp,dirty_phase_deg,'//&
      &'clean_amp,clean_phase_deg,input_amp'
      do j = 1, n
         dirty_amp = sqrt(dirty_re(j)**2 + dirty_im(j)**2)
         dirty_pha_deg = atan2(dirty_im(j), dirty_re(j))*rad2deg
         clean_amp = sqrt(clean_re(j)**2 + clean_im(j)**2)
         clean_pha_deg = atan2(clean_im(j), clean_re(j))*rad2deg
         write(unit,'(F0.4,7(",",F0.6))') rm(j), dirty_re(j), dirty_im(j),&
         &dirty_amp, dirty_pha_deg, clean_amp, clean_pha_deg, input_amp(j)
      end do
      close(unit)
   end subroutine write_profile_csv

   subroutine write_lsq_csv(outdir, slug, l_sq, nchan, band_offset_in, band_nz_in, n_bands_in)
      !! band_offset_in/band_nz_in/n_bands_in: optional, same convention as
      !! run_one's own. When present, a second CSV column records which
      !! band (1-based) each channel belongs to, so the plotting script
      !! can group channels by their TRUE band membership instead of
      !! guessing band boundaries from lambda^2-spacing statistics --
      !! guessing is fragile because a band's own channel-to-channel
      !! lambda^2 spacing (e.g. P-band's, ~0.003 here) can be orders of
      !! magnitude larger than another band's (e.g. L-band's, ~0.0001),
      !! so no single global "gap" threshold can separate real
      !! inter-band gaps from a coarser band's own ordinary channel
      !! spacing (caught empirically: the P-band was shattered into 61
      !! spurious single-point "segments" by exactly this ambiguity).
      character(len=*), intent(in) :: outdir, slug
      integer, intent(in) :: nchan
      real(sp), intent(in) :: l_sq(nchan)
      integer, intent(in), optional :: band_offset_in(:), band_nz_in(:), n_bands_in
      integer :: unit, k, b, band_id(nchan)
      unit = 61
      band_id = 1
      if (present(n_bands_in)) then
         do b = 1, n_bands_in
            band_id(band_offset_in(b)+1:band_offset_in(b)+band_nz_in(b)) = b
         end do
      endif
      open(unit, file=trim(outdir)//'/'//trim(slug)//'_lsq.csv',&
      &status='unknown', action='write')
      write(unit,'(A)') 'lambda_sq,band'
      do k = 1, nchan
         write(unit,'(F0.6,",",I0)') l_sq(k), band_id(k)
      end do
      close(unit)
   end subroutine write_lsq_csv

   function window_flux(rm, spec, n, lo, hi, drm_l) result(flux)
      !! Integrated flux (area under the curve, Riemann sum * bin width)
      !! within [lo,hi] -- the analogue of measuring an EXTENDED source's
      !! total flux by integrating its surface brightness over its own
      !! solid angle, rather than reading off a single peak pixel (which
      !! only makes sense for an unresolved/point source, where
      !! peak=total by convention). Compared against the true simulated
      !! top-hat's own total flux (thick_flux_true = height*width) to
      !! report a recovered FRACTION, not just a bare peak amplitude.
      integer, intent(in) :: n
      real(sp), intent(in) :: rm(n), spec(n), lo, hi, drm_l
      real(sp) :: flux
      integer :: j
      flux = 0.0_sp
      do j = 1, n
         if (rm(j) >= lo .and. rm(j) <= hi) flux = flux + spec(j)*drm_l
      end do
   end function window_flux

   subroutine run_one(label, slug, outdir, l_sq, nchan, rm, nrm_l, thick_flux,&
   &point_peak, point_rm_found, band_offset_in, band_nz_in, n_bands_in)
      !! band_offset_in/band_nz_in/n_bands_in: optional -- when present,
      !! FWHM_RM is computed via compute_rmsf_fwhm_multiband (sum of each
      !! band's own lambda^2-span) instead of compute_rmsf_fwhm's own
      !! single-band formula, which is wrong for a concatenated
      !! multi-band array (see compute_rmsf_fwhm_multiband's own comment
      !! -- caught by exactly this test).
      !! slug/outdir: write <outdir>/<slug>_profile.csv and _lsq.csv,
      !! for tests/plot_thesis_scenario_rmclean.py to render as PNGs in
      !! the thesis figures' own panel style.
      character(len=*), intent(in) :: label, slug, outdir
      integer, intent(in) :: nchan, nrm_l
      real(sp), intent(in) :: l_sq(nchan), rm(nrm_l)
      real(sp), intent(out) :: thick_flux, point_peak, point_rm_found
      integer, intent(in), optional :: band_offset_in(:), band_nz_in(:), n_bands_in

      real(sp) :: q(nchan), u(nchan)
      real(sp) :: dirty_re(nrm_l), dirty_im(nrm_l)
      real(sp) :: comp_re(nrm_l), comp_im(nrm_l), resid_re(nrm_l), resid_im(nrm_l)
      real(sp) :: comp_rm_refined(nrm_l)
      real(sp) :: out_re(nrm_l), out_im(nrm_l), out_amp(nrm_l), fwhm_rm
      real(sp) :: input_amp(nrm_l)
      type(rmsf_table_t) :: table
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: n_iter_used, ipeak, j
      character(len=16) :: stop_reason

      call sky_model_qu(l_sq, nchan, q, u)
      call dirty_spectrum(l_sq, nchan, q, u, rm, nrm_l, dirty_re, dirty_im)

      call build_rmsf_offset_table(l_sq, nchan, 0.0_sp, rm(nrm_l)-rm(1), rm(2)-rm(1), 20, table)
      if (present(n_bands_in)) then
         call compute_rmsf_fwhm_multiband(l_sq, nchan, band_offset_in,&
         &band_nz_in, n_bands_in, fwhm_rm)
      else
         call compute_rmsf_fwhm(l_sq, nchan, fwhm_rm)
      endif
      call plan_fourier_interp(nrm_l, nrm_l, plan_fwd, plan_bwd)

      ! niter-only (no abs_flux_floor/auto_nsigma): exact match for the
      ! old thresh=0.0 (letting niter alone govern convergence, this
      ! scenario's own documented choice for broad Faraday-thick
      ! features -- see planning/RMCLEAN_INTEGRATION_PLAN.md's own
      ! "Choosing parameters" section).
      call clean_complex(l_sq, nchan, 0.0_sp, rm, nrm_l, dirty_re, dirty_im,&
      &table, 1000, 0.1_sp, fwhm_rm, .false., 0.0_sp, .false., 0.0_sp,&
      &3.0_sp, comp_re, comp_im, resid_re, resid_im, n_iter_used,&
      &stop_reason, comp_rm_refined)
      call restore_clean(rm, nrm_l, comp_re, comp_im, resid_re, resid_im,&
      &fwhm_rm, plan_fwd, plan_bwd, out_re, out_im)
      out_amp = sqrt(out_re**2 + out_im**2)

      ! Targeted search around the KNOWN point-source RM, not a blind
      ! global max: with the amplitude-scaling fix above, the extended
      ! component's true total flux (150) exceeds the point source's
      ! (15), so a global max would (correctly, physically) often land
      ! in the extended feature's own window instead once a band is
      ! sensitive enough to recover much of it -- exactly the thesis's
      ! own point. Mirrors check_thesis_scenario.py's own targeted
      ! local_stat search rather than assuming the point source is
      ! always the single largest peak anywhere in the spectrum.
      ipeak = 1
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) < abs(rm(ipeak)-point_rm)) ipeak = j
      end do
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) <= 10.0_sp .and. out_amp(j) > out_amp(ipeak)) ipeak = j
      end do
      point_peak = out_amp(ipeak)
      ! comp_rm_refined(ipeak): the flux-weighted sub-pixel RM location
      ! clean_complex's own CLEAN loop already computed (rmclean_mod's own
      ! comment on this output) -- a real precision improvement over the
      ! bare grid coordinate rm(ipeak), not just a style change (see
      ! derotate_to_lsq_zero's own comment for the empirical evidence).
      point_rm_found = comp_rm_refined(ipeak)

      ! Standard "Jy/beam -> Jy" integrated-flux-density measurement
      ! (the same technique used to measure a resolved source's total
      ! flux from a real restored/CLEANed image): the restoring beam is
      ! normalised to unit PEAK height, not unit area (restore_clean's
      ! own kernel_c(i)=exp(...), height 1 at center -- this is exactly
      ! why a point source's PEAK alone already equals its total flux,
      ! independent of beam width: a delta function convolved with a
      ! height-1 kernel has peak = original flux exactly). Integrating
      ! the restored map directly (as a first attempt here did) instead
      ! double-counts the beam's own area on top of the true flux --
      ! caught empirically (L-alone came back at 1526% of the true
      ! flux, obviously unphysical) before being reported as a real
      ! result. Dividing the integrated restored amplitude by the
      ! beam's own area (sigma*sqrt(2*pi), continuous Gaussian area for
      ! a height-1 kernel) converts the "per-beam" integral back to a
      ! true flux, the same beam-area correction real interferometric
      ! flux measurements apply.
      block
         real(sp) :: sigma_rm, beam_area_rm
         sigma_rm = 0.42466_sp*fwhm_rm
         beam_area_rm = sigma_rm*sqrt(2.0_sp*3.14159265358979_sp)
         thick_flux = window_flux(rm, out_amp, nrm_l, thick_lo, thick_hi, rm(2)-rm(1))&
         &/beam_area_rm
      end block

      write(*,'(A)') '--------------------------------------------------'
      write(*,'(A,A)') 'Case: ', trim(label)
      write(*,'(A,I0,A,F0.3)') '  channels: ', nchan, '  FWHM_RM: ', fwhm_rm
      write(*,'(A,I0)') '  CLEAN iterations used: ', n_iter_used
      write(*,'(A,F0.3,A,F0.3,A,F0.3,A,F0.1,A)') '  Point source: restored peak at RM=',&
      &point_rm_found, ' amp=', point_peak, ' (true amp=', point_amp, ', fraction recovered=',&
      &100.0_sp*point_peak/point_amp, '%)'
      write(*,'(A,F0.3,A,F0.1,A,F0.1,A)') '  Thick component: restored integrated flux in [100,130]=',&
      &thick_flux, ' (true total flux=', thick_flux_true, ', fraction recovered=',&
      &100.0_sp*thick_flux/thick_flux_true, '%)'

      call input_model_amp(rm, nrm_l, input_amp)
      call write_profile_csv(outdir, slug, rm, nrm_l, dirty_re, dirty_im,&
      &out_re, out_im, input_amp)
      call write_lsq_csv(outdir, slug, l_sq, nchan, band_offset_in, band_nz_in, n_bands_in)

      call destroy_rmsf_offset_table(table)
      call destroy_fourier_interp_plan(plan_fwd, plan_bwd)
   end subroutine run_one

   subroutine run_band_case(label, slug, nu_c, dnu, nchan, rm, nrm_l, outdir, all_pass_io)
      character(len=*), intent(in) :: label, slug, outdir
      real(sp), intent(in) :: nu_c, dnu
      integer, intent(in) :: nchan, nrm_l
      real(sp), intent(in) :: rm(nrm_l)
      logical, intent(inout) :: all_pass_io
      real(sp), allocatable :: l_sq(:)
      real(sp) :: thick_flux, point_peak, point_rm_found

      allocate(l_sq(nchan))
      call band_channels(nu_c, dnu, nchan, l_sq)
      call run_one(label, slug, outdir, l_sq, nchan, rm, nrm_l, thick_flux,&
      &point_peak, point_rm_found)
      deallocate(l_sq)

      if (index(label, 'P-alone') > 0) then
         ! Strict, quantitative pass/fail -- reserved for the ISOLATED
         ! Faraday-simple (point) source only, per the user's own
         ! explicit direction: a simple source must be retrieved at its
         ! true RM with nearly equal power to the simulated input. 1.0
         ! rad/m^2 is a tight tolerance relative to this scenario's own
         ! resolution (FWHM_RM~7.8) -- not a loosened-until-it-passes
         ! number (actual recovered RM is within ~0.1); 15% is "nearly
         ! equal power" (actual recovered amplitude is within ~2.5%).
         call check(abs(point_rm_found-point_rm) <= 1.0_sp,&
         &'P-alone: isolated point source recovered within 1 rad/m^2 of true RM=-100', all_pass_io)
         call check(abs(point_peak-point_amp) <= 0.15_sp*point_amp,&
         &'P-alone: isolated point source power within 15% of true simulated input (15.0)', all_pass_io)
      endif
   end subroutine run_band_case

   subroutine run_combined_case(rm, nrm_l, outdir, all_pass_io)
      real(sp), intent(in) :: rm(nrm_l)
      integer, intent(in) :: nrm_l
      character(len=*), intent(in) :: outdir
      logical, intent(inout) :: all_pass_io
      real(sp), allocatable :: l_sq(:)
      real(sp), allocatable :: l_sq_p(:), l_sq_l(:)
      real(sp) :: thick_flux_pl, point_peak_pl, point_rm_found_pl
      real(sp) :: thick_flux_p, point_peak_p, point_rm_found_p
      integer :: n_p, n_l
      integer :: band_offset(2), band_nz(2)

      n_p = 61
      n_l = 121
      allocate(l_sq_p(n_p), l_sq_l(n_l), l_sq(n_p+n_l))
      call band_channels(300.0_sp, 30.0_sp, n_p, l_sq_p)
      call band_channels(1200.0_sp, 120.0_sp, n_l, l_sq_l)
      l_sq(1:n_p) = l_sq_p
      l_sq(n_p+1:n_p+n_l) = l_sq_l
      band_offset = (/0, n_p/)
      band_nz = (/n_p, n_l/)

      call run_one('P+L combined (Fig 6-3)', 'fig6_3_PLcombined', outdir,&
      &l_sq, n_p+n_l, rm, nrm_l, thick_flux_pl, point_peak_pl,&
      &point_rm_found_pl, band_offset, band_nz, 2)

      ! Re-run P-alone here (rather than threading its result through from
      ! run_band_case) to keep this comparison self-contained and exactly
      ! mirror check_thesis_scenario.py's own P-alone-vs-P+L relative
      ! assertion style for the thick component. Written under its own
      ! slug (not fig6_2_Palone) to avoid clobbering that case's own CSVs.
      call run_one('P-alone (re-run for comparison)', 'fig6_2_Palone_rerun',&
      &outdir, l_sq_p, n_p, rm, nrm_l, thick_flux_p, point_peak_p,&
      &point_rm_found_p)

      write(*,'(A)') '--------------------------------------------------'
      write(*,'(A,F0.1,A,F0.1,A)') 'Thick-component recovered flux fraction: P-alone=',&
      &100.0_sp*thick_flux_p/thick_flux_true, '%  P+L=',&
      &100.0_sp*thick_flux_pl/thick_flux_true, '%'

      ! Strict, quantitative pass/fail -- the isolated Faraday-simple
      ! point source, same tolerances as run_band_case's own (see that
      ! subroutine's own comment): must be retrieved at its true RM with
      ! nearly equal power.
      call check(abs(point_rm_found_pl-point_rm) <= 1.0_sp,&
      &'P+L: isolated point source recovered within 1 rad/m^2 of true RM=-100', all_pass_io)
      call check(abs(point_peak_pl-point_amp) <= 0.15_sp*point_amp,&
      &'P+L: isolated point source power within 15% of true simulated input (15.0)', all_pass_io)

      ! Qualitative markers for the EXTENDED (Faraday-thick) component --
      ! deliberately NOT a tight quantitative target (the user's own
      ! direction: do not over-restrict this, look for qualitative
      ! evidence of the thesis's own claim). Two generous, clearly-
      ! qualitative thresholds: combining bands must recover CLEARLY more
      ! than P-alone does (5x is a generous multiplier given the actual
      ! ratio here is ~27x, comfortably robust against run-to-run detail
      ! changes), and P+L's own recovered fraction must be a genuinely
      ! non-negligible ("significant") share of the true input, not
      ! merely "more than a negligible amount more than P-alone's own
      ! near-zero" -- 10% is a generous floor for "significant fraction",
      ! not a target tuned to the ~40% actually observed.
      call check(thick_flux_pl > 5.0_sp*thick_flux_p,&
      &'P+L thick-component recovery clearly exceeds P-alone (>5x) -- qualitative '//&
      &'marker for the thesis''s own "significant fraction of the flux density '//&
      &'for the extended feature ... only recovered when combined" result',&
      &all_pass_io)
      call check(thick_flux_pl/thick_flux_true > 0.10_sp,&
      &'P+L recovers a significant fraction (>10%) of the true simulated top-hat flux',&
      &all_pass_io)

      deallocate(l_sq_p, l_sq_l, l_sq)
   end subroutine run_combined_case

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

end program thesis_scenario_rmclean
