program test_rmclean_lsqref_flex
   !! Validates rmclean_mod capabilities added at the user's explicit
   !! request, now organized around the TWO DELIBERATELY SEPARATE axes the
   !! user asked to keep clean: (1) lsq_ref_compute, a pure COMPUTATIONAL
   !! COST lever -- lsq_ref_compute is a free choice for the CALLER
   !! (already true of compute_dirty_rmbeam_direct/build_rmsf_offset_
   !! table's own signatures), demonstrated here by running the SAME
   !! scenario at two very different compute references -- 0 (this
   !! project's thesis-matching convention, tests/thesis_scenario_
   !! rmclean.f90's own) and the band's own mean (this project's
   !! rm_synthesis_mod.f90's own extract_general_setup convention,
   !! numerically cheaper per suggest_drm); and (2)
   !! lsq_ref_report, a separate REPORTING CONVENTION choice -- via
   !! derotate_to_lsq_ref (generalized from this test's own earlier,
   !! lsq_zero-only version), letting a caller report the intrinsic
   !! polarization angle at EITHER lambda_sq=0 (this project's own
   !! thesis convention) OR at suggest_lsq_ref_report's own Brentjens &
   !! de Bruyn (2005)-style value (chosen there to decorrelate the
   !! estimated RM and chi0's own statistical uncertainty -- a different
   !! optimization criterion from suggest_lsq_ref_compute's own, not a
   !! competing answer to the same question), confirmed here to agree
   !! with an independently-derived closed-form expectation regardless of
   !! which lsq_ref_compute was used internally; and (3) clean_complex's
   !! own comp_rm_refined output, root-caused during this test's own
   !! development to be the actual fix for a real chi0-precision gap at
   !! nonzero lsq_ref_compute (a coarse grid's own restored profile is
   !! grid-quantised, up to dRM/2 away from the true continuous RM even
   !! though it looks perfectly clean and symmetric -- no amount of
   !! re-interpolating that already-quantised output can recover
   !! information already discarded; comp_rm_refined is exactly that
   !! discarded sub-pixel information, recovered). Confirms both
   !! lsq_ref_compute choices now recover chi0 to the SAME tight
   !! tolerance, at EITHER reporting reference, at their own (very
   !! different) grid costs -- the actual answer to "can we keep the
   !! RM-grid coarse and still trust chi0 at any chosen reporting
   !! convention": yes, provided comp_rm_refined is used, not a bare grid
   !! coordinate.
   !!
   !! Sky model: a single Faraday-thin point source, RM=50 rad/m^2,
   !! amplitude=10 Jy/(rad/m^2), intrinsic angle chi0=0.3 rad -- DELIBERATELY
   !! NONZERO (unlike thesis_scenario_rmclean.f90's own PA=0 sky model), so
   !! this test cannot pass merely because "0 in, 0 out" looks right
   !! regardless of a sign error in the derotation formula.
   !!
   !! Single synthetic band (P-band-like, Table 6.1's own nu_c=300/dnu=30/
   !! nchan=61) -- this test exercises the lsq_ref/derotation MACHINERY
   !! itself, not multi-band tomography, so one band is enough and keeps
   !! the comparison focused.
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp
   real(sp), parameter :: point_rm = 50.0_sp, point_amp = 10.0_sp
   real(sp), parameter :: chi0_true = 0.3_sp
   integer, parameter :: nchan = 61
   real(sp), parameter :: rm_drm_oversample = 15.0_sp
   real(sp), parameter :: rm_lo = -100.0_sp, rm_span = 300.0_sp
   real(sp), parameter :: pi_const = 3.14159265_sp

   real(sp) :: l_sq(nchan), q(nchan), u(nchan)
   real(sp) :: drm_zero, drm_centroid, lsq_ref_compute_mean, lsq_ref_report_val
   logical :: all_pass

   all_pass = .true.

   call band_channels(300.0_sp, 30.0_sp, nchan, l_sq)
   call sky_model_qu(l_sq, nchan, q, u)

   ! lsq_ref_compute_mean here is being used as a CHEAP COMPUTE reference
   ! (the "band-mean" case below) -- numerically identical to what
   ! suggest_lsq_ref_report computes (same arithmetic mean formula), but
   ! playing a conceptually DIFFERENT role: a compute-cost lever here,
   ! not a reporting convention. lsq_ref_report_val (below) is the
   ! SEPARATE, reporting-side use of that same formula, called via its
   ! own dedicated routine to exercise that public entry point directly.
   lsq_ref_compute_mean = sum(l_sq)/real(nchan, sp)
   call suggest_lsq_ref_report(l_sq, nchan, lsq_ref_report_val)

   call suggest_drm(l_sq, nchan, 0.0_sp, rm_drm_oversample, drm_zero)
   call suggest_drm(l_sq, nchan, lsq_ref_compute_mean, rm_drm_oversample, drm_centroid)

   write(*,'(A)') '===================================================='
   write(*,'(A,F0.1,A,F0.1,A,F0.2,A)') 'Sky model: point source RM=', point_rm,&
   &' amp=', point_amp, ' chi0=', chi0_true, ' rad'
   write(*,'(A,F0.6,A,F0.6,A)') 'suggest_drm: lsq_ref_compute=0 -> dRM<=',&
   &drm_zero, '   lsq_ref_compute=band-mean -> dRM<=', drm_centroid, ' rad/m^2'
   write(*,'(A,F0.1,A)') '  (band-mean compute reference is ', drm_centroid/drm_zero,&
   &'x coarser -- the computational motivation for choosing it)'
   write(*,'(A,F0.4)') 'suggest_lsq_ref_report (B&dB-style reporting reference): ',&
   &lsq_ref_report_val
   write(*,'(A)') '===================================================='

   call check(drm_centroid > 5.0_sp*drm_zero,&
   &'band-mean compute reference allows a >5x coarser RM grid than lsq_ref_compute=0', all_pass)

   ! chi0 acceptance tolerance: BOTH compute cases now use clean_complex's
   ! own comp_rm_refined (a flux-weighted sub-pixel RM location it already
   ! computes internally every iteration, previously discarded once filed
   ! into an integer grid bin) rather than a bare grid coordinate -- this
   ! decouples chi0 precision from dRM entirely (root-caused empirically:
   ! the earlier ~0.22 rad discrepancy at a coarse, nonzero-lsq_ref_compute
   ! grid was a bookkeeping loss, not a fundamental grid-resolution limit;
   ! using comp_rm_refined at that SAME coarse grid recovered chi0
   ! exactly). Both compute cases therefore get the same tight tolerance,
   ! at EITHER reporting reference.
   call run_case('lsq_ref_compute=0 (thesis convention)', 0.0_sp, drm_zero,&
   &l_sq, q, u, lsq_ref_report_val, 0.02_sp, all_pass)
   call run_case('lsq_ref_compute=band-mean (cheaper grid)', lsq_ref_compute_mean, drm_centroid,&
   &l_sq, q, u, lsq_ref_report_val, 0.02_sp, all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_rmclean_lsqref_flex: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_rmclean_lsqref_flex: one or more checks failed'
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

   subroutine sky_model_qu(lsq_in, n, q_out, u_out)
      !! Physical injection: ALWAYS raw l_sq, no reference at all -- a
      !! real source's own phase cannot depend on a numerical convention
      !! choice (see tests/thesis_scenario_rmclean.f90's own sky_model_qu
      !! comment on why this must never be lsq_ref_compute-dependent).
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
   end subroutine sky_model_qu

   subroutine dirty_spectrum_at_ref(lsq_in, n, q_in, u_in, rm, nrm_l,&
   &lsq_ref_compute, dirty_re, dirty_im)
      !! Same direct-DFT convention as tests/thesis_scenario_rmclean.f90's
      !! own dirty_spectrum, generalized with an explicit lsq_ref_compute
      !! (that file hardcodes 0.0; this test needs both).
      integer, intent(in) :: n, nrm_l
      real(sp), intent(in) :: lsq_in(n), q_in(n), u_in(n), rm(nrm_l), lsq_ref_compute
      real(sp), intent(out) :: dirty_re(nrm_l), dirty_im(nrm_l)
      real(sp) :: c_tmpl(n), s_tmpl(n)
      integer :: j
      do j = 1, nrm_l
         c_tmpl = cos(2.0_sp*rm(j)*(lsq_in-lsq_ref_compute))
         s_tmpl = -sin(2.0_sp*rm(j)*(lsq_in-lsq_ref_compute))
         dirty_re(j) = (dot_product(q_in, c_tmpl) - dot_product(u_in, s_tmpl))/real(n, sp)
         dirty_im(j) = (dot_product(q_in, s_tmpl) + dot_product(u_in, c_tmpl))/real(n, sp)
      end do
   end subroutine dirty_spectrum_at_ref

   subroutine run_case(label, lsq_ref_compute, drm, lsq_in, q_in, u_in,&
   &lsq_ref_report, chi0_tol, all_pass_io)
      character(len=*), intent(in) :: label
      real(sp), intent(in) :: lsq_ref_compute, drm
      real(sp), intent(in) :: lsq_in(nchan), q_in(nchan), u_in(nchan)
      real(sp), intent(in) :: lsq_ref_report, chi0_tol
      logical, intent(inout) :: all_pass_io

      integer :: nrm_l, j, ipeak, n_iter_used
      real(sp), allocatable :: rm(:), dirty_re(:), dirty_im(:)
      real(sp), allocatable :: comp_re(:), comp_im(:), resid_re(:), resid_im(:)
      real(sp), allocatable :: comp_rm_refined(:), comp_amp(:)
      real(sp), allocatable :: out_re(:), out_im(:), out_amp(:)
      real(sp), allocatable :: derot_re(:), derot_im(:)
      real(sp) :: fwhm_rm, amp_before, max_amp_err
      type(rmsf_table_t) :: table
      integer(kind=8) :: plan_fwd, plan_bwd
      real(sp) :: rm_found, chi0_at_zero, chi0_at_report, chi0_expected_at_report
      real(sp) :: comp_re_ref(1), comp_im_ref(1), rm_found_arr(1)
      real(sp) :: derot_re_report(1), derot_im_report(1)

      nrm_l = nint(rm_span/drm) + 1
      allocate(rm(nrm_l), dirty_re(nrm_l), dirty_im(nrm_l))
      allocate(comp_re(nrm_l), comp_im(nrm_l), resid_re(nrm_l), resid_im(nrm_l))
      allocate(comp_rm_refined(nrm_l), comp_amp(nrm_l))
      allocate(out_re(nrm_l), out_im(nrm_l), out_amp(nrm_l))
      allocate(derot_re(nrm_l), derot_im(nrm_l))

      do j = 1, nrm_l
         rm(j) = rm_lo + real(j-1, sp)*drm
      end do

      call dirty_spectrum_at_ref(lsq_in, nchan, q_in, u_in, rm, nrm_l,&
      &lsq_ref_compute, dirty_re, dirty_im)
      call build_rmsf_offset_table(lsq_in, nchan, lsq_ref_compute, rm(nrm_l)-rm(1), drm, 20, table)
      call compute_rmsf_fwhm(lsq_in, nchan, fwhm_rm)
      call plan_fourier_interp(nrm_l, nrm_l, plan_fwd, plan_bwd)

      call clean_complex(rm, nrm_l, dirty_re, dirty_im, table, 500, 0.1_sp,&
      &1.0e-4_sp, comp_re, comp_im, resid_re, resid_im, n_iter_used, comp_rm_refined)
      call restore_clean(rm, nrm_l, comp_re, comp_im, resid_re, resid_im,&
      &fwhm_rm, plan_fwd, plan_bwd, out_re, out_im)
      out_amp = sqrt(out_re**2 + out_im**2)
      comp_amp = sqrt(comp_re**2 + comp_im**2)

      ! Find the recovered COMPONENT (not the restored map's own peak):
      ! comp_rm_refined(j) is a flux-weighted average of the sub-pixel
      ! peak_loc clean_complex already computes every iteration for
      ! accurate beam subtraction -- previously discarded once filed into
      ! its integer grid bin j. Reading chi0 off the restored map's own
      ! (grid-quantised, up to dRM/2 biased) peak location was the actual
      ! bug behind the earlier ~0.22 rad discrepancy at a nonzero
      ! lsq_ref_compute -- not insufficient interpolation, a genuine
      ! bookkeeping loss now fixed at the source (rmclean_mod's own
      ! comp_rm_refined comment).
      ipeak = 1
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) < abs(rm(ipeak)-point_rm)) ipeak = j
      end do
      do j = 1, nrm_l
         if (abs(rm(j)-point_rm) <= 10.0_sp .and. comp_amp(j) > comp_amp(ipeak)) ipeak = j
      end do
      rm_found = comp_rm_refined(ipeak)

      ! Report chi0 at TWO different references via the SAME general
      ! derotate_to_lsq_ref, using the single recovered component
      ! (comp_re/comp_im at ipeak, with rm_found as its own location) --
      ! the array-of-length-1 calling form exercises the exact same code
      ! path a full-array call would, just for one bin.
      comp_re_ref(1) = comp_re(ipeak)
      comp_im_ref(1) = comp_im(ipeak)
      rm_found_arr(1) = rm_found

      call derotate_to_lsq_ref(rm_found_arr, 1, lsq_ref_compute, 0.0_sp,&
      &comp_re_ref, comp_im_ref, derot_re_report, derot_im_report)
      chi0_at_zero = 0.5_sp*atan2(derot_im_report(1), derot_re_report(1))
      chi0_at_zero = chi0_at_zero - nint(chi0_at_zero/pi_const)*pi_const

      call derotate_to_lsq_ref(rm_found_arr, 1, lsq_ref_compute, lsq_ref_report,&
      &comp_re_ref, comp_im_ref, derot_re_report, derot_im_report)
      chi0_at_report = 0.5_sp*atan2(derot_im_report(1), derot_re_report(1))
      chi0_at_report = chi0_at_report - nint(chi0_at_report/pi_const)*pi_const

      ! Independent, closed-form expectation for the report-reference
      ! case (NOT derived from derotate_to_lsq_ref itself, to actually
      ! test it rather than tautologically confirm it): chi0 at reference
      ! R = chi0_true (at reference 0) + RM_true*R, mod pi -- the same
      ! P_A(phi)=P_B(phi)*exp(2i*phi*(B-A))-style relation, evaluated at
      ! A=0, B=R, phi=RM_true.
      chi0_expected_at_report = chi0_true + point_rm*lsq_ref_report
      chi0_expected_at_report = chi0_expected_at_report&
      &- nint(chi0_expected_at_report/pi_const)*pi_const

      ! Separately, confirm derotate_to_lsq_ref's own array-wise identity
      ! (a pure per-grid-point phase rotation, |P| unchanged everywhere --
      ! a property of the transform itself, independent of any peak-
      ! finding precision question above). Applied to the restored map
      ! here specifically to exercise that array-wise path too (not the
      ! single-component path above) -- both are part of the public API
      ! and both need their own coverage.
      call derotate_to_lsq_ref(rm, nrm_l, lsq_ref_compute, 0.0_sp, out_re, out_im,&
      &derot_re, derot_im)

      ! |P| invariance under derotation: a pure phase rotation must not
      ! touch the amplitude, anywhere in the array, not just at the peak.
      amp_before = maxval(sqrt(out_re**2+out_im**2))
      max_amp_err = maxval(abs(sqrt(out_re**2+out_im**2) - sqrt(derot_re**2+derot_im**2)))

      write(*,'(A)') '--------------------------------------------------'
      write(*,'(A,A,A,I0,A,F0.4,A)') 'Case: ', trim(label), ' (nrm=', nrm_l, ', dRM=', drm, ')'
      write(*,'(A,F0.3,A,F0.3,A)') '  Recovered RM=', rm_found, ' (true=', point_rm, ')'
      write(*,'(A,F0.4,A,F0.4,A)') '  chi0 @ lambda_sq=0: ', chi0_at_zero,&
      &' rad (true=', chi0_true, ' rad)'
      write(*,'(A,F0.4,A,F0.4,A)') '  chi0 @ centroid: ', chi0_at_report,&
      &' rad (closed-form expected=', chi0_expected_at_report, ' rad)'
      write(*,'(A,ES10.3)') '  max |P| change under derotation (expect ~0): ', max_amp_err
      write(*,'(A,F0.4,A)') '  chi0 acceptance tolerance for this case: ', chi0_tol, ' rad'

      call check(abs(rm_found-point_rm) <= 1.0_sp,&
      &trim(label)//': point source recovered within 1 rad/m^2 of true RM', all_pass_io)
      call check(abs(chi0_at_zero-chi0_true) <= chi0_tol,&
      &trim(label)//': chi0 @ lambda_sq=0 within its derived tolerance of true intrinsic angle',&
      &all_pass_io)
      call check(abs(chi0_at_report-chi0_expected_at_report) <= chi0_tol,&
      &trim(label)//': chi0 @ centroid matches the independent closed-form expectation',&
      &all_pass_io)
      call check(max_amp_err <= 1.0e-4_sp*amp_before,&
      &trim(label)//': derotation changes |P| by < 1e-4 relative (phase-only rotation)', all_pass_io)

      call destroy_rmsf_offset_table(table)
      call destroy_fourier_interp_plan(plan_fwd, plan_bwd)
      deallocate(rm, dirty_re, dirty_im, comp_re, comp_im, resid_re, resid_im)
      deallocate(out_re, out_im, out_amp, derot_re, derot_im)
   end subroutine run_case

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

end program test_rmclean_lsqref_flex
