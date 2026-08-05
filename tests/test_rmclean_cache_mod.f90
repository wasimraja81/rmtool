program test_rmclean_cache_mod
   !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): standalone unit tests
   !! for rmclean_cache_mod.f90 -- built as this ticket's increments
   !! land (fnv1a_hash first, then the pattern_registry_t machinery).
   !! Deliberately FITS-free and FFTW-free: every check here is pure
   !! in-memory logic, no mask files, no rmclean_mod dependency at all.
   use, intrinsic :: iso_fortran_env, only: int8, int64
   use rmclean_cache_mod, only: fnv1a_hash
   implicit none

   logical :: all_pass
   integer(kind=1) :: empty_bytes(0)
   integer(kind=1) :: seq5(5), seq7(7)
   integer(kind=8) :: h

   all_pass = .true.

   ! Known FNV-1a 64-bit test vectors, independently computed (Python
   ! reference implementation, not derived from this Fortran code) --
   ! see docs/dev/RMCLEAN_INTEGRATION_PLAN.md T14 for the derivation.
   ! Empty byte array: no loop iterations at all, hash must equal the
   ! offset basis constant itself -- also a direct sanity check that
   ! fnv_offset in the implementation IS the standard FNV-1a 64-bit
   ! offset basis (14695981039346656037 unsigned == this signed value).
   h = fnv1a_hash(empty_bytes, 0)
   call check(h.eq.-3750763034362895579_8, 'fnv1a_hash: empty byte array == offset basis', all_pass)

   seq5 = [integer(kind=1) :: 1, 2, 3, 4, 5]
   h = fnv1a_hash(seq5, 5)
   call check(h.eq.1109817072422714760_8, 'fnv1a_hash: [1,2,3,4,5] matches reference value', all_pass)

   ! A genuine 0/1 pattern, the actual shape fnv1a_hash is used for in
   ! practice (one pixel's own valid-channel mask row).
   seq7 = [integer(kind=1) :: 1, 0, 1, 1, 0, 0, 1]
   h = fnv1a_hash(seq7, 7)
   call check(h.eq.2228391952209758033_8, 'fnv1a_hash: [1,0,1,1,0,0,1] matches reference value', all_pass)

   ! Same bytes, called twice -- determinism (pure function, no hidden
   ! state), a basic sanity check worth having explicitly rather than
   ! assuming.
   call check(fnv1a_hash(seq7, 7).eq.fnv1a_hash(seq7, 7), 'fnv1a_hash: deterministic across repeated calls', all_pass)

   ! Different patterns must (for these specific short test vectors,
   ! verified directly, not just "hashes usually differ") produce
   ! different hashes -- catches a degenerate always-same-value bug.
   call check(fnv1a_hash(seq5, 5).ne.fnv1a_hash(seq7, 7), 'fnv1a_hash: distinct patterns hash differently', all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_rmclean_cache_mod: all checks passed'
   else
      write(*,'(A)') '[FAIL] test_rmclean_cache_mod: one or more checks failed'
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

end program test_rmclean_cache_mod
