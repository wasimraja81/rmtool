program test_table_oversample_floor
   !! Validates table_oversample_floor (planning/RMCLEAN_INTEGRATION_
   !! PLAN.md ticket T21): the RMSF table's own Nyquist floor.
   !! Gate 0's min_samples_per_fwhm check only happens to guarantee
   !! this table condition at lsq_ref_compute_mode=mid (max_offset
   !! there is provably span/2, half of Gate 0's own span-based bound)
   !! -- every other mode has no such guarantee, confirmed here with a
   !! constructed case that fails Nyquist by >3x even at the default
   !! table_oversample=20 (docs/user/APP_REFERENCE.md's own
   !! table_oversample row has the full derivation this test checks).
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   real(dp), parameter :: pi_dp = 3.14159265358979_dp
   logical :: all_pass
   integer :: eff

   all_pass = .true.

   ! Case 1: mid mode, real WALLABY+EMU-scale numbers (span=0.097030,
   ! lsq_ref_compute=mid => max_offset=span/2=0.048515, drm=8.6207,
   ! from the real validation run's own .rm_axis.npy) -- Gate-0-
   ! guaranteed case, the floor must NOT raise the default 20.
   block
      real(sp) :: l_sq(2), lsq_ref_compute, drm
      l_sq = [0.043404_sp, 0.140434_sp]
      lsq_ref_compute = 0.5_sp*(l_sq(1)+l_sq(2))
      drm = 8.620689392089787_sp
      eff = table_oversample_floor(l_sq, 2, lsq_ref_compute, drm, 20)
      call check(eff == 20,&
      &'mid mode, real WALLABY+EMU scale: floor does not raise the'&
      &//' default 20 (Gate 0''s own guarantee already covers it)',&
      &all_pass)
   end block

   ! Case 2: zero mode, SAME real WALLABY+EMU numbers -- max_offset=
   ! max(l_sq)=0.140434 (>= span, no Gate-0 guarantee), but the real
   ! margin there is still wide enough that 20 remains sufficient.
   block
      real(sp) :: l_sq(2), drm
      l_sq = [0.043404_sp, 0.140434_sp]
      drm = 8.620689392089787_sp
      eff = table_oversample_floor(l_sq, 2, 0.0_sp, drm, 20)
      call check(eff == 20,&
      &'zero mode, real WALLABY+EMU scale: floor does not need to'&
      &//' raise 20 either here -- real margin, just smaller than mid''s',&
      &all_pass)
   end block

   ! Case 3: the constructed pathological case (span=1, lsq range
   ! [99,100] -- a real ~30 MHz low-frequency band, not an unphysical
   ! number) that breaks even table_oversample=20 at zero mode. The
   ! floor must raise it, and the raised value must genuinely satisfy
   ! >=4 samples/cycle (get_drm's own hard-enforced oversample>=2
   ! standard, not bare Nyquist).
   block
      real(sp) :: l_sq(2), drm
      real(dp) :: max_offset, samples_per_cycle
      l_sq = [99.0_sp, 100.0_sp]
      drm = 1.0_sp
      eff = table_oversample_floor(l_sq, 2, 0.0_sp, drm, 20)
      call check(eff > 20,&
      &'pathological zero-mode case (span=1, lsq=[99,100]): floor'&
      &//' raises table_oversample above the default 20', all_pass)
      max_offset = maxval(abs(real(l_sq, dp)))
      samples_per_cycle = (pi_dp/max_offset) / (real(drm, dp)/real(eff, dp))
      call check(samples_per_cycle >= 4.0_dp - 1.0e-6_dp,&
      &'the raised value genuinely satisfies >=4 samples/cycle'&
      &//' (get_drm''s own oversample>=2 standard)', all_pass)
   end block

   ! Case 4: the floor must never LOWER an already-generous user value.
   block
      real(sp) :: l_sq(2), drm
      l_sq = [0.043404_sp, 0.140434_sp]
      drm = 8.620689392089787_sp
      eff = table_oversample_floor(l_sq, 2, 0.0_sp, drm, 500)
      call check(eff == 500,&
      &'floor never lowers an already-generous user table_oversample',&
      &all_pass)
   end block

   if (all_pass) then
      write(*,'(A)') '[PASS] test_table_oversample_floor: all checks passed'
      stop 0
   else
      write(*,'(A)') '[FAIL] test_table_oversample_floor: one or more checks failed'
      stop 1
   endif

contains

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

end program test_table_oversample_floor
