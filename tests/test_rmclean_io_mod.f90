program test_rmclean_io_mod
   !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): standalone unit test
   !! for rmclean_io_mod.f90's own next_tile_extent -- the tile-geometry
   !! stepper shared by the real per-tile compute loop (Pass 1) and the
   !! new pattern pre-scan (Pass 0). Links against fitsio_unit_mod/
   !! cfitsio (rmclean_io_mod's own module dependency) even though this
   !! specific function never calls a FITS routine -- next_tile_extent
   !! itself is pure, FITS-free logic.
   use rmclean_io_mod, only: next_tile_extent
   implicit none

   logical :: all_pass
   ! Hand-computed expected sequence for nx=10, ny=7, tile_ra=4,
   ! tile_dec=3 -- 9 tiles, RA-strips-first, matching the ORIGINAL
   ! nested `do while (iy_tile_beg.le.ny) / do while (ix_tile_beg.le.nx)`
   ! rule this function replaces (verified by hand, not derived from
   ! the function under test).
   integer, parameter :: n_expected = 9
   integer, parameter :: exp_ix(n_expected) = [1,5,9, 1,5,9, 1,5,9]
   integer, parameter :: exp_iy(n_expected) = [1,1,1, 4,4,4, 7,7,7]
   integer, parameter :: exp_tx(n_expected) = [4,4,2, 4,4,2, 4,4,2]
   integer, parameter :: exp_ty(n_expected) = [3,3,3, 3,3,3, 1,1,1]

   all_pass = .true.
   call walk_and_check(all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_rmclean_io_mod: all checks passed'
   else
      write(*,'(A)') '[FAIL] test_rmclean_io_mod: one or more checks failed'
      stop 1
   endif

contains

   subroutine walk_and_check(all_pass_io)
      !! Replays next_tile_extent's own full sequence for the 10x7/4x3
      !! grid above and checks both the tile COUNT (9, no more no less
      !! -- catches an off-by-one at the grid edge) and each tile's own
      !! origin+extent against the hand-computed expectation. Origin is
      !! captured BEFORE each call, since next_tile_extent both reports
      !! the CURRENT tile and advances (ix,iy) to the NEXT one in the
      !! same call.
      logical, intent(inout) :: all_pass_io
      integer :: ix, iy, tx, ty, ix_before, iy_before, k
      logical :: done
      logical :: sequence_ok, count_ok

      ix = 1
      iy = 1
      k = 0
      sequence_ok = .true.
      count_ok = .true.
      do
         ix_before = ix
         iy_before = iy
         call next_tile_extent(10, 7, 4, 3, ix, iy, tx, ty, done)
         if (done) exit
         k = k + 1
         if (k.gt.n_expected) then
            count_ok = .false.
            exit
         endif
         if (ix_before.ne.exp_ix(k) .or. iy_before.ne.exp_iy(k) .or.&
         &tx.ne.exp_tx(k) .or. ty.ne.exp_ty(k)) then
            sequence_ok = .false.
         endif
      end do
      if (k.ne.n_expected) count_ok = .false.
      call check(count_ok, 'next_tile_extent: 10x7 grid (tile 4x3) produces exactly the expected 9 tiles', all_pass_io)
      call check(sequence_ok, 'next_tile_extent: full tile sequence (origin+extent) matches hand-computed expectation', all_pass_io)
   end subroutine walk_and_check

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

end program test_rmclean_io_mod
