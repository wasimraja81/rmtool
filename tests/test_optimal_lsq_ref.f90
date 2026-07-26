program test_optimal_lsq_ref
   !! Validates rmclean_mod's optimal_lsq_ref_midpoint against the actual
   !! claim it's built on: minimizing required_drm_nyquist's own bound
   !! (i.e. the cheapest safe RM grid) is achieved by the midpoint of the
   !! channel set's own l_sq extent, NOT by a channel-count-weighted mean
   !! and NOT by centring on any single band's own centroid -- added after
   !! the user pushed back with two sharp, specific questions ("is the
   !! midpoint really optimal?", "does channel count not matter?"),
   !! answered here empirically rather than just asserted.
   !!
   !! Reuses the same P-band(61ch)/L-band(121ch) combination as tests/
   !! thesis_scenario_rmclean.f90 and tests/test_rmclean_lsqref_flex.f90 --
   !! deliberately IMBALANCED channel counts (2x), since a channel-count-
   !! weighted mean and the true midpoint-of-extent only coincide for a
   !! single, symmetric, evenly-sampled band; an imbalanced multi-band
   !! case is exactly where they diverge, and exactly where a wrong
   !! ("centre of the lowest-frequency band") intuition looks most
   !! plausible but is actually the worst of the candidates tested here.
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(sp), parameter :: c_light = 299.792458_sp
   integer, parameter :: n_p = 61, n_l = 121
   real(sp) :: l_sq_p(n_p), l_sq_l(n_l), l_sq_all(n_p+n_l)
   real(sp) :: lsq_mean_all, lsq_mean_p, lsq_mean_l, lsq_midpoint, lsq_expected_midpoint
   real(sp) :: drm_zero, drm_mean_all, drm_p_centre, drm_l_centre, drm_midpoint
   logical :: all_pass

   all_pass = .true.

   call band_channels(300.0_sp, 30.0_sp, n_p, l_sq_p)
   call band_channels(1200.0_sp, 120.0_sp, n_l, l_sq_l)
   l_sq_all(1:n_p) = l_sq_p
   l_sq_all(n_p+1:n_p+n_l) = l_sq_l

   lsq_mean_all = sum(l_sq_all)/real(n_p+n_l, sp)
   lsq_mean_p = sum(l_sq_p)/real(n_p, sp)
   lsq_mean_l = sum(l_sq_l)/real(n_l, sp)
   lsq_expected_midpoint = 0.5_sp*(minval(l_sq_all)+maxval(l_sq_all))

   call optimal_lsq_ref_midpoint(l_sq_all, n_p+n_l, lsq_midpoint)

   write(*,'(A)') '===================================================='
   write(*,'(A,F0.4,A,F0.4)') 'overall l_sq range: min=', minval(l_sq_all),&
   &' max=', maxval(l_sq_all)
   write(*,'(A,F0.4)') 'channel-count-weighted mean (182 channels): ', lsq_mean_all
   write(*,'(A,F0.4)') 'P-band-only centroid (61 channels): ', lsq_mean_p
   write(*,'(A,F0.4)') 'L-band-only centroid (121 channels): ', lsq_mean_l
   write(*,'(A,F0.4)') 'optimal_lsq_ref_midpoint result: ', lsq_midpoint
   write(*,'(A)') '===================================================='

   call check(abs(lsq_midpoint - lsq_expected_midpoint) < 1.0e-6_sp,&
   &'optimal_lsq_ref_midpoint returns exactly (min+max)/2', all_pass)

   call required_drm_nyquist(l_sq_all, n_p+n_l, 0.0_sp, 15.0_sp, drm_zero)
   call required_drm_nyquist(l_sq_all, n_p+n_l, lsq_mean_all, 15.0_sp, drm_mean_all)
   call required_drm_nyquist(l_sq_all, n_p+n_l, lsq_mean_p, 15.0_sp, drm_p_centre)
   call required_drm_nyquist(l_sq_all, n_p+n_l, lsq_mean_l, 15.0_sp, drm_l_centre)
   call required_drm_nyquist(l_sq_all, n_p+n_l, lsq_midpoint, 15.0_sp, drm_midpoint)

   write(*,'(A,F0.5,A,I0)') 'lsq_ref=0:                 dRM<=', drm_zero,&
   &'  nrm(span=450)=', nint(450.0_sp/drm_zero)+1
   write(*,'(A,F0.5,A,I0)') 'lsq_ref=weighted mean:     dRM<=', drm_mean_all,&
   &'  nrm(span=450)=', nint(450.0_sp/drm_mean_all)+1
   write(*,'(A,F0.5,A,I0)') 'lsq_ref=P-band centroid:   dRM<=', drm_p_centre,&
   &'  nrm(span=450)=', nint(450.0_sp/drm_p_centre)+1
   write(*,'(A,F0.5,A,I0)') 'lsq_ref=L-band centroid:   dRM<=', drm_l_centre,&
   &'  nrm(span=450)=', nint(450.0_sp/drm_l_centre)+1
   write(*,'(A,F0.5,A,I0)') 'lsq_ref=midpoint (optimal):dRM<=', drm_midpoint,&
   &'  nrm(span=450)=', nint(450.0_sp/drm_midpoint)+1

   ! The whole claim: the midpoint must allow a dRM at least as large
   ! (grid at least as cheap) as every other candidate tested, despite
   ! not being "the middle of the data" in any weighted-average sense.
   call check(drm_midpoint >= drm_zero,&
   &'midpoint allows dRM >= lsq_ref=0''s own bound', all_pass)
   call check(drm_midpoint >= drm_mean_all,&
   &'midpoint allows dRM >= the channel-count-weighted mean''s own bound', all_pass)
   call check(drm_midpoint >= drm_p_centre,&
   &'midpoint allows dRM >= centring on the lowest-frequency (P) band''s own bound', all_pass)
   call check(drm_midpoint >= drm_l_centre,&
   &'midpoint allows dRM >= centring on the L-band''s own bound', all_pass)

   ! And the specific, sharper claim: centring on the lowest-frequency
   ! band is a POOR choice, not a good one -- barely better than
   ! lsq_ref=0, since it merely relocates the worst-case offset to
   ! whichever OTHER band is now farthest from the reference.
   call check(drm_p_centre < 1.3_sp*drm_zero,&
   &'centring on P-band alone is NOT a meaningfully better choice than lsq_ref=0'&
   &//' (confirms channel count/band identity does not matter -- only the overall extremes do)',&
   &all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_optimal_lsq_ref: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_optimal_lsq_ref: one or more checks failed'
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

end program test_optimal_lsq_ref
