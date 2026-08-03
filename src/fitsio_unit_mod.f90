module fitsio_unit_mod
   !! Thread-safe wrappers around CFITSIO's own FTGIOU/FTFIOU (get/free
   !! I/O unit) -- the CFITSIO-native equivalent of Fortran's own
   !! newunit=, needed because FTOPEN/FTINIT are CFITSIO library calls
   !! and don't support newunit= directly. Replaces every hardcoded/
   !! manually-offset unit number (200, 221-227, 300+thread_id,
   !! 1000+omp_get_thread_num(), ...) across convolve_cubes/
   !! reproject_cubes/match_cubes/rmclean_cubes.
   !!
   !! safe_ftopen/safe_ftinit/safe_ftclos bundle FTGIOU+FTOPEN (or
   !! FTINIT) and FTCLOS+FTFIOU into ONE critical section each --
   !! deliberately, not two separate smaller critical sections. Found
   !! the hard way (a real SIGSEGV, reproduced in isolation with a
   !! minimal 8-thread probe: concurrent FTGIOU+FTOPEN+FTGHSP+FTCLOS+
   !! FTFIOU cycles crash reliably within a few hundred iterations
   !! unless FTOPEN/FTCLOS themselves are inside the SAME lock as
   !! FTGIOU/FTFIOU, not just the unit-number bookkeeping alone):
   !! FTGIOU/FTFIOU's own internal free-unit bookkeeping evidently
   !! shares mutable global state with FTOPEN/FTCLOS's own open-file
   !! table, so serializing only the former still lets it race against
   !! another thread's concurrent FTOPEN/FTCLOS. The actual per-plane/
   !! per-tile data transfer calls (FTGSVE/FTPSSE/FTGPVB/FTGREC/...) stay
   !! OUTSIDE these critical sections -- only the structural open/close
   !! bookkeeping needs serializing, not the I/O itself, so
   !! io_read_threads/nwriters/per-thread WCS loading keep their
   !! real parallelism.
   implicit none
   private
   public :: safe_ftopen, safe_ftinit, safe_ftclos, free_fits_unit

contains

   subroutine free_fits_unit(unit)
      !! Releases a unit reserved by safe_ftopen/safe_ftinit when the
      !! FTOPEN/FTINIT call itself then failed (so there is no open file
      !! to FTCLOS -- safe_ftclos does not apply). Idempotent against a
      !! unit already freed (or never assigned).
      integer, intent(inout) :: unit
      integer :: status
      if (unit.le.0) return
      status = 0
      !$omp critical (fitsio_unit_mod_lock)
      call FTFIOU(unit, status)
      !$omp end critical (fitsio_unit_mod_lock)
      unit = -1
   end subroutine free_fits_unit

   subroutine safe_ftopen(unit, filename, rwmode, blocksize, status)
      integer, intent(out) :: unit
      character(len=*), intent(in) :: filename
      integer, intent(in) :: rwmode
      integer, intent(out) :: blocksize
      integer, intent(out) :: status
      status = 0
      !$omp critical (fitsio_unit_mod_lock)
      call FTGIOU(unit, status)
      call FTOPEN(unit, filename, rwmode, blocksize, status)
      !$omp end critical (fitsio_unit_mod_lock)
   end subroutine safe_ftopen

   subroutine safe_ftinit(unit, filename, blocksize, status)
      integer, intent(out) :: unit
      character(len=*), intent(in) :: filename
      integer, intent(in) :: blocksize
      integer, intent(out) :: status
      status = 0
      !$omp critical (fitsio_unit_mod_lock)
      call FTGIOU(unit, status)
      call FTINIT(unit, filename, blocksize, status)
      !$omp end critical (fitsio_unit_mod_lock)
   end subroutine safe_ftinit

   subroutine safe_ftclos(unit, status)
      !! Idempotent against a unit already closed/freed (or never
      !! assigned): left alone, so callers on an error path that already
      !! closed it can call this again without a double-close/free.
      integer, intent(inout) :: unit
      integer, intent(out) :: status
      status = 0
      if (unit.le.0) return
      !$omp critical (fitsio_unit_mod_lock)
      call FTCLOS(unit, status)
      call FTFIOU(unit, status)
      !$omp end critical (fitsio_unit_mod_lock)
      unit = -1
   end subroutine safe_ftclos

end module fitsio_unit_mod
