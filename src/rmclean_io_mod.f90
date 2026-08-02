module rmclean_io_mod
   !! T10: shared, independently-testable FITS I/O helpers for
   !! rmclean_cubes, split out of the program itself so they can be
   !! unit-tested in isolation (planning/RMCLEAN_INTEGRATION_PLAN.md
   !! ticket T10's own regression test needs to call read_mask_cube
   !! directly against a synthetic oversized fixture, which is not
   !! possible for a procedure nested inside a `program`). Intended as
   !! the home for T10b's own future read_mask_tile too, once the mask
   !! cube's whole-cube-resident read is unified with the AMP/PHA tiled
   !! scheme -- see that ticket's own text.
   use fitsio_unit_mod
   implicit none
   private
   public :: read_mask_cube

contains

   subroutine read_mask_cube(filename, nx_in, ny_in, nchan_in, mask_out, status)
      !! T10a: reads the WHOLE mask cube in one FTGPVBLL call (genuine
      !! 64-bit element count/offset) rather than the plain FTGPVB entry
      !! point, whose Fortran-visible firstelem/nelem are narrowed to a
      !! 32-bit LONG even though the underlying C ffgpvb() (and this same
      !! library's FTGSVE, used everywhere else in this codebase) is
      !! LONGLONG-safe throughout -- confirmed directly against this
      !! project's own bundled cfitsio-4.3.1 source (f77_wrap2.c's own
      !! FCALLSCSUB8 declarations for ffgpvb: the plain FTGPVB entry
      !! declares firstelem/nelem as LONG, a SEPARATE FTGPVBLL entry
      !! declares the identical two arguments LONGLONG -- both wrapping
      !! the same C function, so FTGPVBLL is not a slower/different
      !! code path, just the correctly-typed entry point). A real
      !! dataset can trivially exceed 2^31 total mask-cube elements
      !! (nx*ny*nchan) -- e.g. this project's own Jennifer ASKAP mask
      !! cube (4501x4501x288 ~= 5.83 billion) overflowed the plain
      !! FTGPVB's 32-bit nelem, silently truncating the read at ~76 of
      !! 288 channels for EVERY pixel in the image (see planning/
      !! RMCLEAN_INTEGRATION_PLAN.md ticket T10 for the full
      !! investigation/evidence). FTGPVBLL's firstelem/nelem below are
      !! genuine integer(kind=8), so this is correct at any dataset size,
      !! not just today's.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx_in, ny_in, nchan_in
      integer(kind=1), allocatable, intent(out) :: mask_out(:,:,:)
      integer, intent(out) :: status
      integer :: unit, blocksize, fitsstat
      integer(kind=8) :: n_elements_total
      logical :: anyflag

      status = 0
      allocate(mask_out(nx_in,ny_in,nchan_in))
      n_elements_total = int(nx_in,8)*int(ny_in,8)*int(nchan_in,8)
      fitsstat = 0
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif
      call FTGPVBLL(unit, 1, 1_8, n_elements_total, 0_1, mask_out, anyflag,&
      &fitsstat)
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read data from: ', trim(filename)
         status = -1
         return
      endif
   end subroutine read_mask_cube

end module rmclean_io_mod
