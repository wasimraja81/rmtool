module rmclean_io_mod
   !! T10: shared, independently-testable FITS I/O helpers for
   !! rmclean_cubes, split out of the program itself so they can be
   !! unit-tested in isolation. Originally also held read_mask_cube
   !! (T10a's fix for a 32-bit element-count overflow in a whole-cube
   !! FTGPVB read), retired once T10b replaced its only call site with
   !! read_mask_tile below -- FTGSVB (tiled subsection read) has no
   !! equivalent flat pre-multiplied element-count argument at all, so
   !! the overflow class T10a fixed cannot recur in this module; see
   !! planning/RMCLEAN_INTEGRATION_PLAN.md ticket T10 for the full
   !! investigation/evidence of both.
   use fitsio_unit_mod
   implicit none
   private
   public :: read_mask_tile, next_tile_extent

contains

   subroutine next_tile_extent(nx_in, ny_in, tile_ra_in, tile_dec_in,&
   &ix_tile_beg, iy_tile_beg, tx, ty, done)
      !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): given the CURRENT
      !! tile origin (ix_tile_beg, iy_tile_beg), returns that tile's own
      !! (tx, ty) extent and ADVANCES (ix_tile_beg, iy_tile_beg) in
      !! place to the NEXT tile's own origin (RA-strips-first). This is
      !! the exact tile-geometry rule the main program's own
      !! `do while (iy_tile_beg.le.ny) / do while (ix_tile_beg.le.nx)`
      !! nesting implements (rmclean_cubes.f90) -- flattened into one
      !! flat `do while (.not. done)` loop here specifically so it can
      !! be SHARED, verbatim, by both that real per-tile compute loop
      !! (Pass 1) and the new pattern pre-scan (Pass 0, T14) -- the two
      !! passes must never silently diverge on tile geometry, since
      !! Pass-0's own scan-position bookkeeping is only meaningful if
      !! Pass 1 visits pixels in EXACTLY the same order.
      !! done=.true. (tx=ty=0) once iy_tile_beg has already advanced
      !! past ny_in on entry -- the caller's own signal to stop.
      integer, intent(in) :: nx_in, ny_in, tile_ra_in, tile_dec_in
      integer, intent(inout) :: ix_tile_beg, iy_tile_beg
      integer, intent(out) :: tx, ty
      logical, intent(out) :: done

      if (iy_tile_beg.gt.ny_in) then
         done = .true.
         tx = 0
         ty = 0
         return
      endif
      done = .false.
      ty = min(tile_dec_in, ny_in-iy_tile_beg+1)
      tx = min(tile_ra_in, nx_in-ix_tile_beg+1)

      ix_tile_beg = ix_tile_beg + tx
      if (ix_tile_beg.gt.nx_in) then
         ix_tile_beg = 1
         iy_tile_beg = iy_tile_beg + ty
      endif
   end subroutine next_tile_extent

   subroutine split_range_rmclean_io(n_total, n_threads, base, rem)
      !! Same convention as rmclean_cubes.f90's own split_range_rmclean
      !! (itself matching rm_synthesis_mod.f90's split_channels_across_
      !! threads) -- kept as a small private duplicate here rather than
      !! shared via host association, since this module cannot host-
      !! associate into the program that calls it.
      integer, intent(in) :: n_total, n_threads
      integer, intent(out) :: base, rem
      base = n_total / n_threads
      rem = mod(n_total, n_threads)
   end subroutine split_range_rmclean_io

   subroutine read_mask_chunk(maskfile, nx_full, ny_full, nchan_full, ix0,&
   &iy0, tx, ty, chan_beg, chan_len, mask_out, status_par)
      !! One io_read_threads worker's own disjoint channel-chunk of one
      !! mask tile -- see read_mask_tile's own comment. FTGSVB (byte
      !! subsection read) is the mask-cube analogue of FTGSVE (used for
      !! AMP/PHA tile reads, read_amp_pha_chunk in rmclean_cubes.f90):
      !! both take per-axis extents (naxes/blc/trc), never a single
      !! pre-multiplied flat element count, so neither has the 32-bit
      !! overflow exposure T10a fixed in the old whole-cube FTGPVB call
      !! (confirmed directly for FTGSVE's own underlying C implementation,
      !! planning/RMCLEAN_INTEGRATION_PLAN.md ticket T10a's own evidence;
      !! FTGSVB shares the same underlying subsection-read machinery,
      !! parameterised only by output datatype).
      character(len=*), intent(in) :: maskfile
      integer, intent(in) :: nx_full, ny_full, nchan_full
      integer, intent(in) :: ix0, iy0, tx, ty, chan_beg, chan_len
      integer(kind=1), intent(out) :: mask_out(tx,ty,chan_len)
      integer, intent(inout) :: status_par
      integer :: unit, blocksize, fitsstat, group
      integer :: fpixel(3), lpixel(3), incs(3), naxes_full(3)
      logical :: anyflag

      group = 1
      naxes_full = (/ nx_full, ny_full, nchan_full /)
      fpixel = (/ ix0, iy0, chan_beg /)
      lpixel = (/ ix0+tx-1, iy0+ty-1, chan_beg+chan_len-1 /)
      incs = (/ 1, 1, 1 /)

      fitsstat = 0
      call safe_ftopen(unit, trim(maskfile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(maskfile)
         !$omp atomic write
         status_par = -1
         call free_fits_unit(unit)
         return
      endif
      call FTGSVB(unit, group, 3, naxes_full, fpixel, lpixel, incs, 0_1,&
      &mask_out, anyflag, fitsstat)
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read tile data from: ', trim(maskfile)
         !$omp atomic write
         status_par = -1
      endif
   end subroutine read_mask_chunk

   subroutine read_mask_tile(maskfile, nx_full, ny_full, nchan_full, ix0,&
   &iy0, tx, ty, io_read_threads_eff, mask_out, status)
      !! T10b: reads one tile's own mask subregion (full channel depth,
      !! [ix0:ix0+tx-1, iy0:iy0+ty-1, 1:nchan_full]) via FTGSVB, the same
      !! io_read_threads-parallel-chunk scheme read_amp_pha_tile already
      !! uses for AMP/PHA (rmclean_cubes.f90) -- replaces the old
      !! whole-cube-resident read_mask_cube in production use. Correct at
      !! any tile size (up to and including tx=nx_full,ty=ny_full, i.e.
      !! today's former whole-cube-resident behaviour, recovered for free
      !! when a tile that large fits the memory budget) and at any
      !! dataset size, since FTGSVB has no flat-element-count overflow
      !! exposure (see read_mask_chunk's own comment).
      character(len=*), intent(in) :: maskfile
      integer, intent(in) :: nx_full, ny_full, nchan_full, ix0, iy0, tx, ty
      integer, intent(in) :: io_read_threads_eff
      integer(kind=1), intent(out) :: mask_out(tx,ty,nchan_full)
      integer, intent(out) :: status
      integer :: base_chan, rem_chan, ith, chan_beg, chan_len
      integer :: status_par

      status = 0
      status_par = 0
      call split_range_rmclean_io(nchan_full, io_read_threads_eff, base_chan,&
      &rem_chan)

      !$omp parallel do schedule(static) default(shared)&
      !$omp& private(ith,chan_beg,chan_len) num_threads(io_read_threads_eff)
      do ith = 0, io_read_threads_eff-1
         if (ith.lt.rem_chan) then
            chan_len = base_chan + 1
            chan_beg = ith*(base_chan+1) + 1
         else
            chan_len = base_chan
            chan_beg = rem_chan*(base_chan+1) + (ith-rem_chan)*base_chan + 1
         endif
         if (chan_len.gt.0) then
            call read_mask_chunk(maskfile, nx_full, ny_full, nchan_full, ix0,&
            &iy0, tx, ty, chan_beg, chan_len,&
            &mask_out(:,:,chan_beg:chan_beg+chan_len-1), status_par)
         endif
      enddo
      !$omp end parallel do
      if (status_par.ne.0) status = -1
   end subroutine read_mask_tile

end module rmclean_io_mod
