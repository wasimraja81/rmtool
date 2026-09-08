program test_drm_floor_enforcement
   !! Confirms get_drm actually REFUSES an unsafe oversample (<2), not
   !! just that it produces a bad answer if you slip one past it (that's
   !! tests/test_drm_floor.f90's own job). This program is EXPECTED to
   !! be terminated by get_drm's own stop(1) -- tests/run_tests.sh's own
   !! section for this program checks for a NONZERO exit code as the
   !! PASS condition, the opposite of every other test in this suite,
   !! since successful enforcement means this program does not complete
   !! normally.
   use, intrinsic :: iso_fortran_env, only: sp => real32
   use rmclean_mod
   implicit none
   integer, parameter :: nchan = 61
   real(sp) :: l_sq(nchan), drm
   integer :: k

   do k = 1, nchan
      l_sq(k) = 1.0_sp + real(k, sp)*0.001_sp
   end do

   write(*,'(A)') 'Calling get_drm with oversample=1.5 (below the enforced'//&
   &' floor of 2) -- this program should NOT reach the line after this call.'
   call get_drm(l_sq, nchan, 0.0_sp, drm, oversample=1.5_sp)
   write(*,'(A)') 'FAIL: get_drm incorrectly allowed oversample=1.5 to proceed.'
   stop 1
end program test_drm_floor_enforcement
