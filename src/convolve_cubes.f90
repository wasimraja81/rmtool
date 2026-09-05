! convolve_cubes -- standalone common-resolution convolution tool, built on
! gaussft_mod (src/gaussft.f90, pure elliptical-Gaussian FFT-domain
! deconvolve/reconvolve computation) and commonbeam_mod (src/commonbeam.f90,
! smallest-common-beam geometry). This is the "main program" gaussft_mod's
! own header comment anticipated: FITS I/O, per-channel PSF bookkeeping,
! bad-channel handling, and user interaction, mirroring reproject_cubes.f90's
! own split between a narrowly-scoped computational module and a main
! program that drives it.
!
! Per-channel source PSFs: two routes, matching this project's own goal of
! not depending on external tools for basic portability (rmtool should be
! self-contained). Route 1: a CASA-style BEAMS binary table extension
! (CASAMBM=T in the primary header, EXTNAME='BEAMS', columns BMAJ/BMIN/BPA
! in arcsec/arcsec/deg, CHAN 0-indexed) -- confirmed against a real ASKAP
! cube (/data1/tmp/cutout-stokesQ.fits, used by that dataset's own
! full-image rmsynth cfg): CASAMBM=T, 288-row BEAMS table, one row per FREQ
! channel. Route 2:
! a plain ASCII text file, one line per channel: "channel bmaj_arcsec
! bmin_arcsec bpa_deg" (1-indexed channel, matching this project's existing
! bad-channel-file convention), '#'-prefixed or blank lines skipped -- for
! cubes with no BEAMS table, or for a user who wants to override/hand-supply
! per-channel PSFs without editing FITS binary tables. Whitespace- or
! comma-separated (or a mix), so a plain CSV export works unchanged --
! see read_beams_ascii's own comment for why no separate CSV parsing
! path is needed. Two ready-to-adapt examples, using real ASKAP
! per-channel beam values, so nobody has to reinvent this format from
! this description alone: cfg/example_beamLog.txt (aligned columns,
! human-readable at a glance) and cfg/example_beamLog.csv (comma-
! separated). A degenerate BEAMS-
! table entry (BMAJ effectively zero -- confirmed against the real cube:
! CASA writes ~1.18e-38, the smallest normal single-precision float, for a
! channel with no valid restoring beam) is treated as a bad channel, same as
! this project's existing global_badchan_file convention (verified: the
! real cube's 2 degenerate-beam channels, CHAN 160/177 0-indexed, are
! exactly the 2 channels -- 161/178 1-indexed -- already listed in
! cfg/askap_nan_channels.burdies) -- a bad channel's output plane is written
! as all-NaN, not convolved, same policy as rm_synthesis's own bad-channel
! handling.
!
! Target beam: by default, automatically derived via commonbeam_mod as the
! smallest common beam every GOOD channel of every input file can be
! deconvolved from -- pooled across ALL input files (one file for a
! single-band run, several for multi-band; multi-band support needs no
! extra machinery here: every channel of every input's own per-channel
! beam is simply added to the same pool before finding one common beam,
! then every file is convolved to that single shared target -- matching
! the physical requirement that the target resolution be identical
! across every channel, whether those channels come from one band or
! several). Can also be set explicitly (target_bmaj/target_bmin/
! target_bpa), which skips auto-derivation entirely -- gaussft_mod itself
! does no target-vs-source beam validity checking (a deliberate scope
! decision, see its own header comment: "We can of course have target beam
! less than PSF"), so an explicit target is the user's own call, unchecked.
! An auto-derived common beam can optionally be capped (max_common_bmaj) --
! if it comes out fatter than that, this program stops and refuses to
! proceed rather than silently convolving to a resolution the user never
! sanity-checked; this check does NOT apply to an explicit target (already
! an explicit user decision) or to commonbeam_mod itself (pure geometry, no
! policy -- the same reasoning that kept the earlier "is target >= source"
! check out of gaussft_mod also keeps this one out of commonbeam_mod: policy
! belongs here, in the main program, not in a computational module).
!
! BPA sky-to-pixel conversion: gaussft_mod's own bpa_in/bmaj_in convention
! is "the angle to rotate the pixel-frame (ix,iy) major axis by" (its own
! header comment: "this module does no coordinate-system reasoning of its
! own, it just rotates by the angle it's given"), whereas a FITS BMAJ/BMIN/
! BPA keyword (and this program's own BEAMS-table/ASCII readers) is the
! standard radio-astronomy convention: position angle measured in the SKY
! plane, from North, increasing through East. For an axis-aligned pixel
! grid (CDELT1/CDELT2 only, no CROTA/PC/CD rotation -- checked and refused
! loudly otherwise, see read_axis_info below), the local tangent-plane
! pixel-frame direction of North is sign(CDELT2)*(+iy), and of East is
! sign(CDELT1)*(+ix) (CDELT1 is standardly negative -- +ix is then WEST, so
! East is -ix; confirmed against the real cube's own header: CDELT1<0,
! CDELT2>0). A sky-frame unit vector at position angle theta has components
! (sin(theta) East, cos(theta) North); converting that into the pixel-frame
! unit vector gaussft_mod's own (cos(bpa_pixel), sin(bpa_pixel)) convention
! expects (see ellipse_edges' identical convention in commonbeam.f90) gives
! bpa_pixel = atan2(sign(CDELT2)*cos(theta), sign(CDELT1)*sin(theta)).
! Sanity-checked against 2 special cases for the real cube's own sign
! combination (CDELT1<0, CDELT2>0): theta=0 (pure North) -> bpa_pixel=90
! degrees -> gaussft's (cos90,sin90)=(0,1) = pure +iy = North, correct;
! theta=90 (pure East) -> bpa_pixel=180 degrees -> (cos180,sin180)=(-1,0) =
! pure -ix = East (since CDELT1<0 here), correct.
!
! Usage: convolve_cubes infiles=<file1>[,<file2>...] [outsuffix=<suffix>]
!    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]
!    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
!    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]
!    [npts=<n>] [khachiyan_tol=<tol>]
!    or: convolve_cubes --config <cfgfile>
!    or: convolve_cubes --help | -h
! Full usage text in print_usage below (shared by --help and the
! argument-error path, same convention as reproject_cubes.f90).
program convolve_cubes
   use, intrinsic :: iso_fortran_env, only: dp => real64, int8, int32, int64
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_ptr, c_funptr,&
   &c_null_ptr, c_funloc, c_loc, c_f_pointer
   use logging_mod
   use fitsio_unit_mod
   implicit none

   ! --- Logging & timing (planning-doc ticket, same convention as rm_
   ! synthesis's own log_level/timing_enabled/log_output_file, ported via
   ! the shared logging_mod rather than rm_synthesis_mod.f90 directly --
   ! see logging_mod.f90's own header comment for why). Stage names used
   ! below: startup, header, block_read, block_convolve, block_write,
   ! finalize.
   character(len=16) :: log_level
   logical :: timing_enabled
   character(len=272) :: log_output_file

   integer, parameter :: max_axes = 10
   integer, parameter :: max_inputs = 50
   integer, parameter :: max_channels = 20000
   ! Hard ceiling on elements-per-block (nx*ny*block_planes), independent
   ! of mem_frac_ram: CFITSIO's own Fortran wrapper (FTGSVE/FTPSSE) takes
   ! a default-INTEGER element count, so a single call describing more
   ! than 2^31-1 elements would silently wrap/misbehave on a real
   ! full-size image with a generous mem_frac_ram. 2e9 leaves comfortable
   ! headroom under 2,147,483,647 for ANY nx/ny, not just today's real
   ! 4501x4501 validation cube -- this clamp activates (shrinks
   ! block_planes below what mem_frac_ram alone would allow) only for
   ! genuinely huge images/generous mem_frac_ram combinations; ordinary
   ! runs never notice it.
   integer(kind=8), parameter :: max_elements_per_block = 2000000000_8

   ! --- io_overlap: background-thread block write (planning-doc ticket)
   ! --- see write_convolved_file's own comment for the single-writer-
   ! at-a-time design rationale. block_write_job_t is fully self-
   ! contained (unlike rmclean_cubes.f90's own tile_write_job_t, which
   ! leans on host association to program-level nx/ny/out_unit etc.)
   ! because here the block loop lives inside write_convolved_file, a
   ! subroutine called once per input file with its own LOCAL nx/ny/
   ! out_unit/naxes -- not program-level state a sibling procedure could
   ! reach by host association alone.
   type :: block_write_job_t
      ! file_path/datastart (not out_unit): raw Fortran stream I/O at a
      ! computed byte offset (T19, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.
      ! md), not a live CFITSIO handle -- see write_convolved_file's own
      ! comment for why. naxis/sky1/sky2/freq_axis/naxes are no longer
      ! needed here.
      character(len=1024) :: file_path = ' '
      integer(kind=8) :: datastart = 0_8
      integer :: chan_start = 0, chan_len = 0, nx = 0, ny = 0
      integer :: nwriters_eff = 1
      real, pointer :: data(:,:,:) => null()
   end type block_write_job_t
   ! target+save: must remain valid (untouched, undeallocated) from
   ! dispatch until block_write_join returns -- see
   ! block_write_dispatch_async's own comment. One instance reused
   ! across every input file (convolve_cubes processes files strictly
   ! one after another, never concurrently, so this is never shared
   ! across two in-flight files).
   type(block_write_job_t), target, save :: write_job
   integer(c_long) :: write_thread_id = 0
   logical :: write_pending = .false.
   ! Set by do_block_write on a write failure -- read back by
   ! write_convolved_file right after block_write_join returns (pthread_
   ! join's own synchronization makes this safe to read then, without
   ! needing an atomic: the background thread has already fully exited).
   logical :: write_failed = .false.

   interface
      function c_pthread_create(thread, attr, start_routine, arg)&
      &bind(C, name="pthread_create") result(rc)
         import :: c_int, c_long, c_ptr, c_funptr
         integer(c_long) :: thread
         type(c_ptr), value :: attr
         type(c_funptr), value :: start_routine
         type(c_ptr), value :: arg
         integer(c_int) :: rc
      end function c_pthread_create

      function c_pthread_join(thread, retval)&
      &bind(C, name="pthread_join") result(rc)
         import :: c_int, c_long, c_ptr
         integer(c_long), value :: thread
         type(c_ptr), value :: retval
         integer(c_int) :: rc
      end function c_pthread_join
   end interface

   character(len=512) :: infiles(max_inputs), beamfiles(max_inputs)
   integer :: n_inputs
   character(len=64) :: outsuffix
   ! badchan_file: one entry per infile, same comma-list convention as
   ! beamfiles -- each infile gets its own independent bad-channel list
   ! (empty entry = none for that infile). Genuinely per-band, unlike the
   ! single shared list this used to be (T16, planning/MULTI_BAND_
   ! TOMOGRAPHY_PLAN.md): a bad-channel index used to apply identically to
   ! every infile regardless of that infile's own channel numbering.
   character(len=512) :: badchan_file(max_inputs)
   logical :: have_target
   real(dp) :: target_bmaj, target_bmin, target_bpa
   real(dp) :: max_common_bmaj
   logical :: have_max_common_bmaj
   real :: mem_frac_ram
   ! io_overlap (default n): background-thread block write, overlapped
   ! with the NEXT block's read+convolve -- same scheme/key name as rm_
   ! synthesis/rmclean_cubes' own io_overlap (planning-doc ticket, added
   ! after a ~46GB end-to-end run measured convolve's own
   ! write I/O -- a real disk, 150MB/s measured -- as fully serial with
   ! compute, ~44s/block of dead time). See write_convolved_file's own
   ! comment for why this is a single-writer-at-a-time design (block_out
   ! is double-buffered so read+compute for the NEXT block can proceed
   ! immediately, but only one background write is ever in flight on
   ! out_unit at a time -- concurrent CFITSIO handle use from two threads
   ! is unsafe, same hazard rm_synthesis_mod.f90's own nwriters
   ! doc comment describes).
   logical :: io_overlap
   ! nwriters (default 1, T21 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! split one block's own chan_len planes into nwriters_eff disjoint
   ! chunks, each written concurrently via write_freq_block_raw -- same
   ! key name/clamp formula as rm_synthesis/rmclean_cubes' own nwriters
   ! (T6/T12): max(1, min(nwriters, omp_get_max_threads())), further
   ! bounded by chan_len itself inside write_convolved_file. Orthogonal
   ! to io_overlap -- io_overlap decides WHEN a block's write runs
   ! (inline vs. overlapped with the next block's read+compute);
   ! nwriters decides how many concurrent writers do that write once
   ! dispatched. See do_block_write's own comment for how they compose.
   integer :: nwriters
   ! dry_run (default n, T24 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! same key name as rm_synthesis' own dry_run -- checks the target
   ! output disk's rotational status and writes a suggested
   ! convolve_cubes_dryrun.cfg (io_overlap/nwriters) instead of
   ! processing any real file. See run_dry_run_check's own comment.
   logical :: dry_run
   integer :: npts
   real(dp) :: khachiyan_tol

   integer :: i, status

   ! Per-file axis/beam bookkeeping, gathered up front for every input
   ! before any convolution happens (needed to pool per-channel beams
   ! across ALL files before find_common_beam is called once).
   integer :: naxis_f(max_inputs), sky1_f(max_inputs), sky2_f(max_inputs)
   integer :: freq_axis_f(max_inputs), naxes_f(max_inputs, max_axes)
   real(dp) :: cdelt1_f(max_inputs), cdelt2_f(max_inputs)
   integer :: nfreq_f(max_inputs)
   real(dp), allocatable :: bmaj_f(:,:), bmin_f(:,:), bpa_f(:,:)
   logical, allocatable :: isbad_f(:,:)

   ! --- Skip-if-already-matched (planning-doc ticket, same logic as
   ! match_cubes.f90's own -- see there for the full rationale) ---
   logical :: needs_processing(max_inputs), out_exists

   real(dp), allocatable :: pool_bmaj(:), pool_bmin(:), pool_bpa(:)
   integer :: n_pool
   real(dp) :: common_bmaj, common_bmin, common_bpa
   integer :: badchan_list(max_channels), n_badchan

   call parse_args(status)
   if (status.ne.0) stop 1

   call init_logging(log_level, timing_enabled, log_output_file, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: cannot open log_output_file: ', trim(log_output_file)
      stop 1
   endif
   call log_message('info', 'startup', 'convolve_cubes run started')

   if (dry_run) then
      call run_dry_run_check(infiles(1), 'convolve_cubes')
      stop
   endif

   allocate(bmaj_f(max_inputs, max_channels), bmin_f(max_inputs, max_channels))
   allocate(bpa_f(max_inputs, max_channels), isbad_f(max_inputs, max_channels))

   do i = 1, n_inputs
      call read_axis_info(infiles(i), naxis_f(i), sky1_f(i), sky2_f(i),&
      &freq_axis_f(i), naxes_f(i,:), cdelt1_f(i), cdelt2_f(i), status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read axis info for: ', trim(infiles(i))
         stop 1
      endif
      nfreq_f(i) = naxes_f(i, freq_axis_f(i))
      if (nfreq_f(i).gt.max_channels) then
         write(*,*) 'ERROR: ', trim(infiles(i)), ' has ', nfreq_f(i),&
         &' channels, exceeding this program''s max_channels=', max_channels
         stop 1
      endif

      call read_beams(infiles(i), trim(beamfiles(i)), nfreq_f(i),&
      &bmaj_f(i,1:nfreq_f(i)), bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)),&
      &isbad_f(i,1:nfreq_f(i)), status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read per-channel beams for: ', trim(infiles(i))
         stop 1
      endif

      n_badchan = 0
      if (len_trim(badchan_file(i)).gt.0 .and. trim(badchan_file(i)).ne.'none') then
         call read_badchan_file(badchan_file(i), badchan_list, n_badchan, status)
         if (status.ne.0) stop 1
      endif
      call apply_badchan_list(badchan_list, n_badchan, nfreq_f(i), isbad_f(i,1:nfreq_f(i)))

      write(*,'(A,A,A,I0,A,I0,A)') 'Read ', trim(infiles(i)), ': ', nfreq_f(i),&
      &' channels, ', count(isbad_f(i,1:nfreq_f(i))), ' flagged bad'
   enddo

   if (have_target) then
      common_bmaj = target_bmaj
      common_bmin = target_bmin
      common_bpa = target_bpa
      write(*,'(A,F0.4,A,F0.4,A,F0.4)') 'Using explicit target beam: BMAJ=',&
      &common_bmaj, ' BMIN=', common_bmin, ' PA=', common_bpa
   else
      n_pool = 0
      do i = 1, n_inputs
         n_pool = n_pool + count(.not. isbad_f(i,1:nfreq_f(i)))
      enddo
      if (n_pool.lt.1) then
         write(*,*) 'ERROR: no good (non-bad) channels across any input file'
         stop 1
      endif
      allocate(pool_bmaj(n_pool), pool_bmin(n_pool), pool_bpa(n_pool))
      call pool_good_beams(n_inputs, nfreq_f, bmaj_f, bmin_f, bpa_f, isbad_f,&
      &max_inputs, max_channels, n_pool, pool_bmaj, pool_bmin, pool_bpa)

      call find_common_beam_wrap(n_pool, pool_bmaj, pool_bmin, pool_bpa,&
      &npts, khachiyan_tol, common_bmaj, common_bmin, common_bpa, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: could not find a common beam deconvolvable from',&
         &' every good input channel'
         stop 1
      endif
      write(*,'(A,I0,A,F0.4,A,F0.4,A,F0.4)') 'Derived common beam from ',&
      &n_pool, ' good channels: BMAJ=', common_bmaj, ' BMIN=', common_bmin,&
      &' PA=', common_bpa
      deallocate(pool_bmaj, pool_bmin, pool_bpa)

      if (have_max_common_bmaj .and. common_bmaj.gt.max_common_bmaj) then
         write(*,'(A,F0.4,A,F0.4,A)') 'ERROR: derived common beam BMAJ=',&
         &common_bmaj, ' arcsec exceeds max_common_bmaj=', max_common_bmaj,&
         &' arcsec -- refusing to proceed. Investigate why the required'//&
         &' common resolution is this coarse (e.g. an outlier per-channel'//&
         &' beam not already flagged bad), or raise max_common_bmaj if this'//&
         &' resolution is genuinely intended.'
         stop 1
      endif
   endif

   ! === Skip-if-already-matched pre-flight (planning-doc ticket) ===
   ! Same scheme as match_cubes.f90's own (see there for the full
   ! rationale): every safety check and skip decision for the WHOLE
   ! batch happens up front, before any file is processed, so a bad run
   ! (a stale output already on disk) fails fast. A pre-existing output
   ! path is always refused, never silently reused or overwritten,
   ! regardless of what this run's own skip decision would have been.
   do i = 1, n_inputs
      inquire(file=trim(strip_fits_ext(infiles(i)))//trim(outsuffix), exist=out_exists)
      if (out_exists) then
         write(*,*) 'ERROR: output path already exists, refusing to proceed'//&
         &' (stale output from a previous run? remove it first): ',&
         &trim(strip_fits_ext(infiles(i)))//trim(outsuffix)
         stop 1
      endif
      needs_processing(i) = .not. beam_matches_target(nfreq_f(i),&
      &bmaj_f(i,1:nfreq_f(i)), bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)),&
      &isbad_f(i,1:nfreq_f(i)), common_bmaj, common_bmin, common_bpa)
      if (.not. needs_processing(i)) then
         write(*,'(A,A,A)') 'SKIP: ', trim(infiles(i)),&
         &' already matches the target beam -- no output written, use it directly'
      endif
   enddo

   do i = 1, n_inputs
      if (.not. needs_processing(i)) cycle
      call write_convolved_file(infiles(i), trim(strip_fits_ext(infiles(i)))//trim(outsuffix),&
      &naxis_f(i), sky1_f(i), sky2_f(i), freq_axis_f(i), naxes_f(i,:),&
      &cdelt1_f(i), cdelt2_f(i), nfreq_f(i), bmaj_f(i,1:nfreq_f(i)),&
      &bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)), isbad_f(i,1:nfreq_f(i)),&
      &common_bmaj, common_bmin, common_bpa, mem_frac_ram, io_overlap,&
      &nwriters, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write convolved output for: ', trim(infiles(i))
         stop 1
      endif
      write(*,*) 'OK: wrote ', trim(strip_fits_ext(infiles(i)))//trim(outsuffix)
   enddo

   deallocate(bmaj_f, bmin_f, bpa_f, isbad_f)
   write(*,*) 'OK: all inputs convolved to common resolution.'
   call timer_report_summary()
   call log_message('info', 'finalize', 'convolve_cubes run completed')

contains

   subroutine pool_good_beams(n_inputs, nfreq, bmaj, bmin, bpa, isbad,&
   &dim1, dim2, n_pool, pool_bmaj, pool_bmin, pool_bpa)
      integer, intent(in) :: n_inputs, dim1, dim2, nfreq(dim1)
      real(dp), intent(in) :: bmaj(dim1,dim2), bmin(dim1,dim2), bpa(dim1,dim2)
      logical, intent(in) :: isbad(dim1,dim2)
      integer, intent(in) :: n_pool
      real(dp), intent(out) :: pool_bmaj(n_pool), pool_bmin(n_pool), pool_bpa(n_pool)
      integer :: ii, jj, k

      k = 0
      do ii = 1, n_inputs
         do jj = 1, nfreq(ii)
            if (.not. isbad(ii,jj)) then
               k = k + 1
               pool_bmaj(k) = bmaj(ii,jj)
               pool_bmin(k) = bmin(ii,jj)
               pool_bpa(k) = bpa(ii,jj)
            endif
         enddo
      enddo
   end subroutine pool_good_beams

   subroutine find_common_beam_wrap(n, bmaj, bmin, bpa, npts_in, tol_in,&
   &out_bmaj, out_bmin, out_bpa, status)
      !! Thin wrapper so the commonbeam_mod use-association only needs to
      !! appear once (real*8 <-> the module's own real64 are the same
      !! kind on every platform this project targets, so no conversion
      !! is needed, just an explicit interface boundary).
      use commonbeam_mod, only: find_common_beam
      integer, intent(in) :: n, npts_in
      real(dp), intent(in) :: bmaj(n), bmin(n), bpa(n), tol_in
      real(dp), intent(out) :: out_bmaj, out_bmin, out_bpa
      integer, intent(out) :: status

      call find_common_beam(n, bmaj, bmin, bpa, npts_in, tol_in,&
      &out_bmaj, out_bmin, out_bpa, status)
   end subroutine find_common_beam_wrap

   logical function beam_matches_target(nfreq, bmaj, bmin, bpa, isbad,&
   &common_bmaj, common_bmin, common_bpa) result(matches)
      !! Ticket: memory-in-plan-file "match_cubes: skip already-matched
      !! files". True only if EVERY good channel of this one input
      !! already has (bmaj,bmin,bpa) equal to the common target beam
      !! within a tight absolute tolerance -- see this ticket's own
      !! tolerance rationale (planning doc): "already-processed-
      !! identically", not "close enough to convolve negligibly". A bad
      !! (flagged) channel is skipped in this check the same way it's
      !! skipped everywhere else in this file (write_convolved_file's own
      !! all-NaN convention for bad channels) -- it carries no beam
      !! information either way.
      integer, intent(in) :: nfreq
      real(dp), intent(in) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: common_bmaj, common_bmin, common_bpa
      real(dp), parameter :: tol_beam = 1.0d-6
      integer :: k

      matches = .true.
      do k = 1, nfreq
         if (isbad(k)) cycle
         if (abs(bmaj(k)-common_bmaj).gt.tol_beam .or.&
         &abs(bmin(k)-common_bmin).gt.tol_beam .or.&
         &abs(bpa(k)-common_bpa).gt.tol_beam) then
            matches = .false.
            return
         endif
      enddo
   end function beam_matches_target

   function flag_from_value_convolve(val) result(flag)
      !! Same convention as rm_synthesis_mod.f90's own flag_from_value:
      !! first non-blank character '1'/'y'/'Y'/'t'/'T' -> true, anything
      !! else (including blank) -> false.
      character(len=*), intent(in) :: val
      logical :: flag
      character(len=64) :: t
      integer :: i

      t = adjustl(val)
      do i = 1, len_trim(t)
         t(i:i) = achar(iachar(t(i:i)) + merge(32, 0, t(i:i).ge.'A'.and.t(i:i).le.'Z'))
      enddo
      flag = .false.
      if (len_trim(t).eq.0) return
      if (t(1:1).eq.'1' .or. t(1:1).eq.'y' .or. t(1:1).eq.'t') flag = .true.
   end function flag_from_value_convolve

   function strip_fits_ext(filename) result(base)
      !! Output-name helper: infile with its trailing extension (whatever
      !! follows the last '.' in its basename -- .fits, .FITS, .FITSCUBE,
      !! ...) removed, so outsuffix lands where a human expects it -- e.g.
      !! "cutout-stokesQ.fits" + outsuffix "_CONV.FITS" gives "cutout-
      !! stokesQ_CONV.FITS", not the double-extension "cutout-
      !! stokesQ.fits_CONV.FITS" this used to produce by simply
      !! concatenating outsuffix onto the raw infile string. Generic by
      !! design (not hardcoded to ".fits") so any input extension, this
      !! project's own .FITSCUBE test fixtures included, is handled the
      !! same way. A filename with no '.' in its basename (only in a
      !! parent directory component) is returned unchanged.
      character(len=*), intent(in) :: filename
      character(len=len(filename)) :: base
      integer :: n, i, slash, dot

      base = filename
      n = len_trim(filename)
      if (n.lt.1) return

      slash = 0
      do i = n, 1, -1
         if (filename(i:i).eq.'/') then
            slash = i
            exit
         endif
      enddo

      dot = 0
      do i = n, slash+1, -1
         if (filename(i:i).eq.'.') then
            dot = i
            exit
         endif
      enddo

      if (dot.gt.slash+1) base = filename(1:dot-1)
   end function strip_fits_ext

   subroutine parse_args(status)
      integer, intent(out) :: status
      ! this_arg/cli_val and the raw_* CSV-list staging variables below
      ! hold a WHOLE comma-separated infiles=/beamfiles=/badchan_file=
      ! argument, not one path -- up to max_inputs (50) entries, each up
      ! to 512 chars (infiles(:)'s own per-entry length). 16384 gives
      ! generous headroom over the worst case (~25,600 chars); 512 here
      ! silently truncated a real run's own infiles= argument (found via
      ! scripts/run_pipeline.sh's symlink-redirection scheme, which
      ! lengthens every path enough to cross 512 for just 4 real files).
      character(len=16384) :: this_arg, cli_val
      character(len=512) :: cli_key, cfgfile
      character(len=16384) :: raw_infiles, raw_beamfiles, raw_badchan_file
      integer :: argc, iarg, ios
      logical :: has_kv, have_cfgfile, seen_infiles

      status = 0
      n_inputs = 0
      outsuffix = '_CONV.FITS'
      badchan_file = ' '
      have_target = .false.
      target_bmaj = 0.0d0
      target_bmin = 0.0d0
      target_bpa = 0.0d0
      have_max_common_bmaj = .false.
      max_common_bmaj = 0.0d0
      mem_frac_ram = 0.25
      io_overlap = .false.
      nwriters = 1
      dry_run = .false.
      log_level = 'info'
      timing_enabled = .false.
      log_output_file = ' '
      npts = 2000
      khachiyan_tol = 1.0d-5
      raw_infiles = ' '
      raw_beamfiles = ' '
      raw_badchan_file = ' '
      have_cfgfile = .false.
      seen_infiles = .false.

      argc = command_argument_count()
      if (argc.eq.0) then
         call print_usage()
         status = -1
         return
      endif

      iarg = 1
      do while (iarg.le.argc)
         call get_command_argument(iarg, this_arg)
         if (trim(this_arg).eq.'--help' .or. trim(this_arg).eq.'-h') then
            call print_usage()
            status = -1
            return
         else if (trim(this_arg).eq.'--config') then
            if (iarg.eq.argc) then
               write(*,*) 'ERROR: --config requires a file path argument'
               status = -1
               return
            endif
            call get_command_argument(iarg+1, cfgfile)
            have_cfgfile = .true.
            iarg = iarg + 2
         else
            call split_cli_kv(this_arg, cli_key, cli_val, has_kv)
            if (.not. has_kv) then
               write(*,*) 'ERROR: unrecognised argument "', trim(this_arg),&
               &'" -- expected key=value, --config <file>, or --help'
               status = -1
               return
            endif
            call apply_kv(trim(cli_key), trim(cli_val), raw_infiles,&
            &raw_beamfiles, raw_badchan_file, seen_infiles, status)
            if (status.ne.0) return
            iarg = iarg + 1
         endif
      enddo

      if (have_cfgfile) call read_cfg_file(cfgfile, raw_infiles, raw_beamfiles,&
      &raw_badchan_file, seen_infiles, status)
      if (status.ne.0) return

      if (.not. seen_infiles) then
         call print_usage()
         status = -1
         return
      endif

      n_inputs = cfg_csv_count(raw_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.max_inputs) then
         write(*,*) 'ERROR: infiles must list between 1 and ', max_inputs, ' files'
         status = -1
         return
      endif
      do i = 1, n_inputs
         call cfg_csv_get_item(raw_infiles, i, infiles(i))
         beamfiles(i) = 'auto'
      enddo
      if (len_trim(raw_beamfiles).gt.0) then
         if (cfg_csv_count(raw_beamfiles).ne.n_inputs) then
            write(*,*) 'ERROR: beamfiles must list exactly ', n_inputs,&
            &' entries (one per infile; use ''auto'' for a file''s own',&
            &' BEAMS table)'
            status = -1
            return
         endif
         do i = 1, n_inputs
            call cfg_csv_get_item(raw_beamfiles, i, beamfiles(i))
         enddo
      endif
      if (len_trim(raw_badchan_file).gt.0) then
         if (cfg_csv_count(raw_badchan_file).ne.n_inputs) then
            write(*,*) 'ERROR: badchan_file must list exactly ', n_inputs,&
            &' entries (one per infile; give an entry the literal value',&
            &' ''none'' -- e.g. ''none,file2.txt'' -- for a file with no',&
            &' manual bad-channel list)'
            status = -1
            return
         endif
         do i = 1, n_inputs
            call cfg_csv_get_item(raw_badchan_file, i, badchan_file(i))
         enddo
         ! T34 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): once badchan_file=
         ! is given at all, every position must be explicit -- a real path
         ! or the literal value 'none' -- matching rm_synthesis's own T33
         ! convention. Omitting badchan_file= entirely is untouched (this
         ! whole block is skipped, raw_badchan_file stays blank); this only
         ! rejects a blank POSITION inside an otherwise-given list, which
         ! was previously silently accepted as "nothing to list."
         do i = 1, n_inputs
            if (len_trim(badchan_file(i)).eq.0) then
               write(*,*) 'ERROR: badchan_file entry ', i, ' is blank --',&
               &' use the literal value ''none'' for a file with nothing',&
               &' to list'
               status = -1
               return
            endif
         enddo
      endif

      if (mem_frac_ram.le.0.0 .or. mem_frac_ram.gt.0.95) then
         write(*,*) 'ERROR: mem_frac_ram must be > 0 and <= 0.95, got ', mem_frac_ram
         status = -1
         return
      endif
      ios = 0
      if (npts.lt.12) then
         write(*,*) 'ERROR: npts must be at least 12, got ', npts
         status = -1
         return
      endif
   end subroutine parse_args

   subroutine apply_kv(key, val, raw_infiles, raw_beamfiles, raw_badchan_file,&
   &seen_infiles, status)
      character(len=*), intent(in) :: key, val
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles, raw_badchan_file
      logical, intent(inout) :: seen_infiles
      integer, intent(out) :: status
      integer :: ios

      status = 0
      select case (key)
      case ('infiles')
         raw_infiles = val
         seen_infiles = .true.
      case ('beamfiles')
         raw_beamfiles = val
      case ('outsuffix')
         outsuffix = val
      case ('badchan_file')
         raw_badchan_file = val
      case ('target_bmaj')
         read(val, *, iostat=ios) target_bmaj
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bmaj must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('target_bmin')
         read(val, *, iostat=ios) target_bmin
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bmin must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('target_bpa')
         read(val, *, iostat=ios) target_bpa
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bpa must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('max_common_bmaj')
         read(val, *, iostat=ios) max_common_bmaj
         if (ios.ne.0) then
            write(*,*) 'ERROR: max_common_bmaj must be a number'
            status = -1
            return
         endif
         have_max_common_bmaj = .true.
      case ('mem_frac_ram')
         read(val, *, iostat=ios) mem_frac_ram
         if (ios.ne.0) then
            write(*,*) 'ERROR: mem_frac_ram must be a number'
            status = -1
            return
         endif
      case ('io_overlap')
         io_overlap = flag_from_value_convolve(val)
      case ('nwriters')
         read(val, *, iostat=ios) nwriters
         if (ios.ne.0 .or. nwriters.lt.1) then
            write(*,*) 'ERROR: nwriters must be an integer >= 1'
            status = -1
            return
         endif
      case ('dry_run')
         dry_run = flag_from_value_convolve(val)
      case ('log_level')
         log_level = trim(val)
      case ('timing_enabled')
         timing_enabled = flag_from_value_logging(val)
      case ('log_output_file')
         log_output_file = trim(val)
      case ('npts')
         read(val, *, iostat=ios) npts
         if (ios.ne.0) then
            write(*,*) 'ERROR: npts must be an integer'
            status = -1
            return
         endif
      case ('khachiyan_tol')
         read(val, *, iostat=ios) khachiyan_tol
         if (ios.ne.0) then
            write(*,*) 'ERROR: khachiyan_tol must be a number'
            status = -1
            return
         endif
      case default
         write(*,*) 'ERROR: unrecognised key "', key, '"'
         status = -1
         return
      end select
   end subroutine apply_kv

   subroutine read_cfg_file(cfgfile, raw_infiles, raw_beamfiles, raw_badchan_file,&
   &seen_infiles, status)
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles, raw_badchan_file
      logical, intent(inout) :: seen_infiles
      integer, intent(out) :: status
      character(len=16384) :: line, val
      character(len=512) :: key
      integer :: unit_cfg, ios, line_no
      logical :: has_kv

      status = 0
      open(newunit=unit_cfg, file=trim(cfgfile), status='old', action='read', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open config file: ', trim(cfgfile)
         status = -1
         return
      endif
      line_no = 0
      do
         read(unit_cfg, '(A)', iostat=ios) line
         if (ios.ne.0) exit
         line_no = line_no + 1
         call cfg_split_key_value(line, key, val, has_kv)
         if (.not. has_kv) cycle
         call apply_kv(trim(key), trim(val), raw_infiles, raw_beamfiles,&
         &raw_badchan_file, seen_infiles, status)
         if (status.ne.0) then
            write(*,*) '  (at line ', line_no, ' in ', trim(cfgfile), ')'
            close(unit_cfg)
            return
         endif
      enddo
      close(unit_cfg)
   end subroutine read_cfg_file

   subroutine print_usage()
      write(*,'(A)') 'convolve_cubes -- convolve FITS cubes to a common angular resolution'
      write(*,'(A)') ''
      write(*,'(A)') 'Usage:'
      write(*,'(A)') '  convolve_cubes infiles=<file1>[,<file2>...] [outsuffix=<suffix>]'
      write(*,'(A)') '    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]'
      write(*,'(A)') '    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]'
      write(*,'(A)') '    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]'
      write(*,'(A)') '    [npts=<n>] [khachiyan_tol=<tol>] [io_overlap=y|n] [nwriters=<n>]'//&
      &' [dry_run=y|n]'
      write(*,'(A)') '  convolve_cubes --config <cfgfile>'
      write(*,'(A)') '  convolve_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'infiles: 1 or more FITS cubes (comma-separated), each with 2 sky axes'//&
      &' and exactly one other non-degenerate axis (FREQ, CTYPE starting'
      write(*,'(A)') '  ''FREQ''); every other axis (e.g. STOKES) must have extent 1 -- run'//&
      &' separate Stokes/etc slices as separate infiles.'
      write(*,'(A)') ''
      write(*,'(A)') 'outsuffix: appended to each infile''s own path for its output filename'//&
      &' (default _CONV.FITS).'
      write(*,'(A)') ''
      write(*,'(A)') 'beamfiles: per-channel source PSF for each infile, in the same order.'//&
      &' Each entry is either the literal word ''auto'' (read that infile''s own'
      write(*,'(A)') '  CASA-style BEAMS binary table extension) or a path to an ASCII text'//&
      &' file, one line per channel: "channel bmaj_arcsec bmin_arcsec bpa_deg"'
      write(*,'(A)') '  (1-indexed channel; ''#''-prefixed or blank lines skipped;'//&
      &' whitespace- or comma-separated, so a plain CSV export works too). Omit'
      write(*,'(A)') '  entirely to use ''auto'' for every infile. See'//&
      &' cfg/example_beamLog.txt and cfg/example_beamLog.csv for ready examples.'
      write(*,'(A)') ''
      write(*,'(A)') 'badchan_file: per-infile bad-channel list, in the same order as'//&
      &' infiles (comma-separated, one entry per infile -- a blank entry is a'
      write(*,'(A)') '  hard error; give it the literal value ''none'' instead, e.g.'//&
      &' ''none,file2.txt'', for an infile with no manual list of its own).'
      write(*,'(A)') '  Each file is the same one-integer-per-line, 1-indexed'//&
      &' convention as'
      write(*,'(A)') '  rm_synthesis''s own global_badchan_file. Independent of, and in'//&
      &' addition to, automatic bad-channel detection from each infile''s own'
      write(*,'(A)') '  BEAMS table (a degenerate near-zero beam entry) -- use this for'//&
      &' channels known bad for reasons a beam table alone would not capture'
      write(*,'(A)') '  (e.g. RFI). Either way, a bad channel is written as an all-NaN'//&
      &' plane, not convolved.'
      write(*,'(A)') ''
      write(*,'(A)') 'target_bmaj/target_bmin/target_bpa: explicit target beam (all three'//&
      &' required together) -- skips automatic common-beam derivation entirely.'
      write(*,'(A)') ''
      write(*,'(A)') 'max_common_bmaj: if the AUTO-derived common beam''s BMAJ exceeds this'//&
      &' (arcsec), refuse to proceed rather than silently convolve to an'
      write(*,'(A)') '  unexpectedly coarse resolution. Ignored when target_bmaj/etc is given'//&
      &' explicitly (already an explicit user decision).'
      write(*,'(A)') ''
      write(*,'(A)') 'mem_frac_ram (default 0.25): fraction of total system RAM budgeted for'//&
      &' one read/convolve/write block of planes at a time, same concept as'
      write(*,'(A)') '  reproject_cubes'' and rm_synthesis'' own mem_frac_ram.'
      write(*,'(A)') ''
      write(*,'(A)') 'npts (default 2000), khachiyan_tol (default 1e-5): passed straight to'//&
      &' commonbeam_mod''s find_common_beam -- boundary points sampled per beam and'
      write(*,'(A)') '  Khachiyan-algorithm convergence tolerance for the common-beam fit.'
      write(*,'(A)') ''
      write(*,'(A)') 'io_overlap (default n): y -- write each block on a background thread,'//&
      &' overlapped with the NEXT block''s own read+convolve, instead of blocking'
      write(*,'(A)') '  on the write before starting it. Only one background write is'//&
      &' ever in flight at a time.'
      write(*,'(A)') ''
      write(*,'(A)') 'nwriters (default 1): split one block''s own planes into this many'//&
      &' disjoint concurrent writers -- same key/clamp as rm_synthesis'' and'
      write(*,'(A)') '  rmclean_cubes'' own nwriters (max(1, min(nwriters, OMP thread count))).'//&
      &' Orthogonal to io_overlap: io_overlap decides WHEN a block''s write runs'
      write(*,'(A)') '  (inline vs. overlapped with the next block''s read+convolve);'//&
      &' nwriters decides how many concurrent writers do it once dispatched.'
      write(*,'(A)') ''
      write(*,'(A)') 'dry_run (default n): y -- check the target output disk''s own'//&
      &' rotational status and write a suggested convolve_cubes_dryrun.cfg'
      write(*,'(A)') '  (io_overlap/nwriters), instead of processing any real file.'//&
      &' Advisory only. Touches no data.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional keys (CLI or config, logging/timing):'
      write(*,'(A)') '  log_level       = error|warn|info|debug (default info)'
      write(*,'(A)') '  timing_enabled  = y|n -- print a stage timing summary (default n)'
      write(*,'(A)') '  log_output_file = path -- append log/timing output to this file'//&
      &' instead of stdout (default empty = stdout)'
   end subroutine print_usage

   subroutine cfg_split_key_value(raw_line, key, val, has_kv)
      character(len=*), intent(in) :: raw_line
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      character(len=len(raw_line)) :: line
      integer :: p1, p2, peq, pcut

      key = ' '
      val = ' '
      has_kv = .false.
      line = raw_line
      p1 = index(line, ';')
      p2 = index(line, '#')
      if (p1 > 0 .and. p2 > 0) then
         pcut = min(p1, p2)
      else if (p1 > 0) then
         pcut = p1
      else
         pcut = p2
      endif
      if (pcut > 0) line = line(1:pcut - 1)
      line = adjustl(line)
      if (len_trim(line) == 0) return
      peq = index(line, '=')
      if (peq <= 1) return
      key = adjustl(line(1:peq - 1))
      val = adjustl(line(peq + 1:))
      if (len_trim(key) == 0 .or. len_trim(val) == 0) return
      key = trim(key)
      val = trim(val)
      has_kv = .true.
   end subroutine cfg_split_key_value

   subroutine split_cli_kv(token, key, val, has_kv)
      character(len=*), intent(in) :: token
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      integer :: peq

      key = ' '
      val = ' '
      has_kv = .false.
      peq = index(token, '=')
      if (peq <= 1) return
      key = adjustl(token(1:peq - 1))
      val = adjustl(token(peq + 1:))
      if (len_trim(key) == 0 .or. len_trim(val) == 0) return
      key = trim(key)
      val = trim(val)
      has_kv = .true.
   end subroutine split_cli_kv

   function cfg_csv_count(str) result(n)
      character(len=*), intent(in) :: str
      integer :: n, ii

      n = 0
      if (len_trim(str) == 0) return
      n = 1
      do ii = 1, len_trim(str)
         if (str(ii:ii) == ',') n = n + 1
      enddo
   end function cfg_csv_count

   subroutine cfg_csv_get_item(str, idx, item)
      character(len=*), intent(in) :: str
      integer, intent(in) :: idx
      character(len=*), intent(out) :: item
      integer :: ii, cur, p0, n

      item = ' '
      n = len_trim(str)
      if (n == 0) return
      cur = 1
      p0 = 1
      do ii = 1, n
         if (str(ii:ii) == ',') then
            if (cur == idx) then
               item = adjustl(str(p0:ii - 1))
               return
            endif
            cur = cur + 1
            p0 = ii + 1
         endif
      enddo
      if (cur == idx) item = adjustl(str(p0:n))
   end subroutine cfg_csv_get_item

   subroutine read_badchan_file(filename, list, n, status)
      !! Same one-integer-per-line, 1-indexed convention as rm_synthesis's
      !! own global_badchan_file (see rm_synthesis.f90's own bad-channel
      !! read loop) -- confirmed against cfg/askap_nan_channels.burdies.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: list(:)
      integer, intent(out) :: n
      integer, intent(out) :: status
      integer :: unit_bc, ios

      status = 0
      n = 0
      open(newunit=unit_bc, file=trim(filename), status='old', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open badchan_file: ', trim(filename)
         status = -1
         return
      endif
      do
         if (n.ge.size(list)) then
            write(*,*) 'ERROR: too many entries in badchan_file (max ', size(list), ')'
            status = -1
            close(unit_bc)
            return
         endif
         n = n + 1
         read(unit_bc, *, iostat=ios) list(n)
         if (ios.ne.0) then
            n = n - 1
            exit
         endif
      enddo
      close(unit_bc)
   end subroutine read_badchan_file

   subroutine apply_badchan_list(list, n, nfreq, isbad)
      integer, intent(in) :: list(:), n, nfreq
      logical, intent(inout) :: isbad(nfreq)
      integer :: ii

      do ii = 1, n
         if (list(ii).ge.1 .and. list(ii).le.nfreq) isbad(list(ii)) = .true.
      enddo
   end subroutine apply_badchan_list

   subroutine read_axis_info(filename, naxis, sky1, sky2, freq_axis, naxes,&
   &cdelt1, cdelt2, status)
      !! Identify the 2 sky axes (CTYPE starting 'RA'/'DEC') and the FREQ
      !! axis (CTYPE starting 'FREQ') by inspecting CTYPEn directly via
      !! CFITSIO -- no AST needed here (unlike reproject_cubes, this tool
      !! never resamples/reprojects, so it only needs axis ROLES and pixel
      !! SCALE, not a full WCS Mapping). Every other axis must have extent
      !! 1 (a documented scope limit: run separate Stokes/etc slices as
      !! separate infiles, matching how rm_synthesis itself is invoked on
      !! single-Stokes cubes). Also refuses (loud error, not a silent
      !! mishandling) any CROTA or off-diagonal PC/CD rotation on the sky
      !! axes -- this program's own sky-to-pixel BPA conversion (see this
      !! file's own top-of-file comment) assumes an axis-aligned grid.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: naxis, sky1, sky2, freq_axis
      integer, intent(out) :: naxes(max_axes)
      real(dp), intent(out) :: cdelt1, cdelt2
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, k
      character(len=68) :: ctype, comment
      character(len=8) :: axstr

      status = 0
      fitsstat = 0
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.lt.2 .or. naxis.gt.max_axes) then
         write(*,*) 'ERROR: bad or missing NAXIS in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      sky1 = 0
      sky2 = 0
      freq_axis = 0
      naxes = 0
      do k = 1, naxis
         write(axstr,'(I0)') k
         fitsstat = 0
         call FTGKYJ(unit, 'NAXIS'//trim(axstr), naxes(k), comment, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: missing NAXIS', k, ' in: ', trim(filename)
            status = -1
            call safe_ftclos(unit, fitsstat)
            return
         endif
         fitsstat = 0
         call FTGKYS(unit, 'CTYPE'//trim(axstr), ctype, comment, fitsstat)
         if (fitsstat.ne.0) cycle
         ctype = adjustl(ctype)
         if (ctype(1:2).eq.'RA') then
            if (sky1.eq.0) then
               sky1 = k
            else
               write(*,*) 'ERROR: more than one RA-like axis in: ', trim(filename)
               status = -1
               call safe_ftclos(unit, fitsstat)
               return
            endif
         else if (ctype(1:3).eq.'DEC') then
            if (sky2.eq.0) then
               sky2 = k
            else
               write(*,*) 'ERROR: more than one DEC-like axis in: ', trim(filename)
               status = -1
               call safe_ftclos(unit, fitsstat)
               return
            endif
         else if (ctype(1:4).eq.'FREQ') then
            freq_axis = k
         endif
      enddo

      if (sky1.eq.0 .or. sky2.eq.0) then
         write(*,*) 'ERROR: could not identify RA/DEC sky axes in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif
      if (freq_axis.eq.0) then
         write(*,*) 'ERROR: could not identify a FREQ axis in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      do k = 1, naxis
         if (k.ne.sky1 .and. k.ne.sky2 .and. k.ne.freq_axis) then
            if (naxes(k).gt.1) then
               write(*,*) 'ERROR: ', trim(filename), ' axis ', k,&
               &' has extent ', naxes(k), ' > 1 -- only 2 sky axes plus one'//&
               &' FREQ axis are supported; run separate slices (e.g. per'//&
               &' Stokes) as separate infiles'
               status = -1
               call safe_ftclos(unit, fitsstat)
               return
            endif
         endif
      enddo

      write(axstr,'(I0)') sky1
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt1, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the RA axis in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the DEC axis in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      call check_no_rotation(unit, sky1, sky2, filename, status)
      if (status.ne.0) then
         call safe_ftclos(unit, fitsstat)
         return
      endif

      call safe_ftclos(unit, fitsstat)
   end subroutine read_axis_info

   subroutine check_no_rotation(unit, sky1, sky2, filename, status)
      !! Refuse (loudly) any CROTA or off-diagonal PC/CD rotation on the
      !! sky axes -- see this file's own top comment for why this
      !! program's sky-to-pixel BPA formula requires an axis-aligned grid.
      !! A diagonal PC (PCi_i=1, or absent -- FITS default) is fine and
      !! common (confirmed on the real ASKAP cube: PC1_1=PC2_2=1, no
      !! off-diagonal entries, no CROTA at all).
      integer, intent(in) :: unit, sky1, sky2
      character(len=*), intent(in) :: filename
      integer, intent(out) :: status
      integer :: fitsstat
      character(len=68) :: comment
      character(len=8) :: a1, a2
      real(dp) :: dval

      status = 0
      write(a1,'(I0)') sky1
      fitsstat = 0
      call FTGKYD(unit, 'CROTA'//trim(a1), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero CROTA',&
         &trim(a1), ' -- rotated sky grids are not supported by this'//&
         &' program''s sky-to-pixel BPA conversion; reproject onto an'//&
         &' axis-aligned grid first (see reproject_cubes)'
         status = -1
         return
      endif
      write(a2,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CROTA'//trim(a2), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero CROTA',&
         &trim(a2), ' -- rotated sky grids are not supported'
         status = -1
         return
      endif

      call check_offdiag(unit, 'PC', sky1, sky2, filename, status)
      if (status.ne.0) return
      call check_offdiag(unit, 'CD', sky1, sky2, filename, status)
   end subroutine check_no_rotation

   subroutine check_offdiag(unit, prefix, sky1, sky2, filename, status)
      integer, intent(in) :: unit, sky1, sky2
      character(len=*), intent(in) :: prefix, filename
      integer, intent(out) :: status
      integer :: fitsstat
      character(len=68) :: comment
      character(len=16) :: k12, k21
      character(len=8) :: a1, a2
      real(dp) :: dval

      status = 0
      write(a1,'(I0)') sky1
      write(a2,'(I0)') sky2
      k12 = trim(prefix)//trim(a1)//'_'//trim(a2)
      fitsstat = 0
      call FTGKYD(unit, trim(k12), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero ', trim(k12),&
         &' -- rotated/sheared sky grids are not supported'
         status = -1
         return
      endif
      k21 = trim(prefix)//trim(a2)//'_'//trim(a1)
      fitsstat = 0
      call FTGKYD(unit, trim(k21), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero ', trim(k21),&
         &' -- rotated/sheared sky grids are not supported'
         status = -1
         return
      endif
   end subroutine check_offdiag

   subroutine read_beams(filename, beamspec, nfreq, bmaj, bmin, bpa, isbad, status)
      character(len=*), intent(in) :: filename, beamspec
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      if (trim(beamspec).eq.'auto' .or. len_trim(beamspec).eq.0) then
         call read_beams_table(filename, nfreq, bmaj, bmin, bpa, isbad, status)
      else
         call read_beams_ascii(beamspec, nfreq, bmaj, bmin, bpa, isbad, status)
      endif
   end subroutine read_beams

   subroutine read_beams_table(filename, nfreq, bmaj, bmin, bpa, isbad, status)
      !! CASA-style BEAMS binary table extension (EXTNAME='BEAMS',
      !! columns BMAJ/BMIN/BPA in arcsec/arcsec/deg, CHAN 0-indexed) --
      !! see this file's own top comment for the real-cube verification.
      !! A degenerate row (BMAJ or BMIN < 1.0e-6 arcsec -- the real cube's
      !! own placeholder is ~1.18e-38, comfortably below any real beam) is
      !! flagged bad, not treated as a read error.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, nrows, k
      integer :: col_bmaj, col_bmin, col_bpa, col_chan
      real, allocatable :: rb_bmaj(:), rb_bmin(:), rb_bpa(:)
      integer, allocatable :: rb_chan(:)
      logical :: anyflag
      integer :: chan1

      status = 0
      fitsstat = 0
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      fitsstat = 0
      call FTMNHD(unit, -1, 'BEAMS', 0, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: no BEAMS binary table extension found in: ',&
         &trim(filename), ' -- pass an ASCII beamfile instead (see --help)'
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      fitsstat = 0
      call FTGNRW(unit, nrows, fitsstat)
      if (fitsstat.ne.0 .or. nrows.lt.1) then
         write(*,*) 'ERROR: could not read row count of BEAMS table in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      fitsstat = 0
      call FTGCNO(unit, .false., 'BMAJ', col_bmaj, fitsstat)
      call FTGCNO(unit, .false., 'BMIN', col_bmin, fitsstat)
      call FTGCNO(unit, .false., 'BPA', col_bpa, fitsstat)
      call FTGCNO(unit, .false., 'CHAN', col_chan, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: BEAMS table in ', trim(filename),&
         &' missing one of BMAJ/BMIN/BPA/CHAN columns'
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      allocate(rb_bmaj(nrows), rb_bmin(nrows), rb_bpa(nrows), rb_chan(nrows))
      fitsstat = 0
      call FTGCVE(unit, col_bmaj, 1, 1, nrows, 0.0, rb_bmaj, anyflag, fitsstat)
      call FTGCVE(unit, col_bmin, 1, 1, nrows, 0.0, rb_bmin, anyflag, fitsstat)
      call FTGCVE(unit, col_bpa, 1, 1, nrows, 0.0, rb_bpa, anyflag, fitsstat)
      call FTGCVJ(unit, col_chan, 1, 1, nrows, 0, rb_chan, anyflag, fitsstat)
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed reading BEAMS table columns in: ', trim(filename)
         status = -1
         deallocate(rb_bmaj, rb_bmin, rb_bpa, rb_chan)
         return
      endif

      isbad = .true.
      bmaj = 0.0d0
      bmin = 0.0d0
      bpa = 0.0d0
      do k = 1, nrows
         chan1 = rb_chan(k) + 1
         if (chan1.lt.1 .or. chan1.gt.nfreq) cycle
         bmaj(chan1) = real(rb_bmaj(k), kind=8)
         bmin(chan1) = real(rb_bmin(k), kind=8)
         bpa(chan1) = real(rb_bpa(k), kind=8)
         isbad(chan1) = (bmaj(chan1).lt.1.0d-6 .or. bmin(chan1).lt.1.0d-6)
      enddo
      deallocate(rb_bmaj, rb_bmin, rb_bpa, rb_chan)
   end subroutine read_beams_table

   subroutine read_beams_ascii(beamfile, nfreq, bmaj, bmin, bpa, isbad, status)
      !! One line per channel: "channel bmaj_arcsec bmin_arcsec bpa_deg"
      !! (1-indexed channel, '#'-prefixed/blank lines skipped -- see this
      !! file's own top comment). A channel never listed is flagged bad
      !! (no PSF known for it), same as a listed channel with BMAJ or BMIN
      !! set to 0 (or any value below 1e-6 arcsec) -- explicitly present
      !! but degenerate, same policy as a degenerate BEAMS-table row.
      !!
      !! Whitespace- or comma-separated, or a mix of both, on the same
      !! line -- no separate CSV code path needed, since Fortran's
      !! list-directed read (the `read(line, *, ...)` below) already
      !! treats commas and blanks as equivalent value separators
      !! (verified directly: "1, 14.0, 12.0, 70.0", "1 14.0 12.0 70.0"
      !! and "1,14.0,12.0,70.0" all parse identically). A conventional
      !! CSV header row (e.g. "channel,bmaj_arcsec,bmin_arcsec,bpa_deg")
      !! must still start with '#' to be skipped, like any other comment
      !! line here -- one simple, consistent rule rather than guessing
      !! whether an unmarked first line is a header. See
      !! cfg/example_beamLog.txt (aligned columns) and
      !! cfg/example_beamLog.csv (comma-separated) for ready-to-adapt
      !! examples using real ASKAP per-channel beam values, so a user
      !! never has to reinvent this format from the description alone.
      character(len=*), intent(in) :: beamfile
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      integer :: unit_bf, ios
      character(len=512) :: line
      integer :: ich
      real(dp) :: a, b, p

      status = 0
      isbad = .true.
      bmaj = 0.0d0
      bmin = 0.0d0
      bpa = 0.0d0

      open(newunit=unit_bf, file=trim(beamfile), status='old', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open ASCII beam file: ', trim(beamfile)
         status = -1
         return
      endif
      do
         read(unit_bf, '(A)', iostat=ios) line
         if (ios.ne.0) exit
         line = adjustl(line)
         if (len_trim(line).eq.0) cycle
         if (line(1:1).eq.'#') cycle
         read(line, *, iostat=ios) ich, a, b, p
         if (ios.ne.0) then
            write(*,*) 'ERROR: malformed line in ', trim(beamfile), ': ', trim(line)
            write(*,*) '  expected "channel bmaj_arcsec bmin_arcsec bpa_deg"',&
            &' (whitespace- or comma-separated) -- an unmarked header row',&
            &' must start with ''#'' to be skipped, see cfg/example_beamLog.csv'
            status = -1
            close(unit_bf)
            return
         endif
         if (ich.lt.1 .or. ich.gt.nfreq) then
            write(*,*) 'ERROR: channel ', ich, ' in ', trim(beamfile),&
            &' is out of range 1..', nfreq
            status = -1
            close(unit_bf)
            return
         endif
         bmaj(ich) = a
         bmin(ich) = b
         bpa(ich) = p
         isbad(ich) = (a.lt.1.0d-6 .or. b.lt.1.0d-6)
      enddo
      close(unit_bf)
   end subroutine read_beams_ascii

   subroutine sky_to_pixel_bpa(bpa_sky_deg, cdelt1, cdelt2, bpa_pixel_deg)
      !! See this file's own top comment for the full derivation.
      real(dp), intent(in) :: bpa_sky_deg, cdelt1, cdelt2
      real(dp), intent(out) :: bpa_pixel_deg
      real(dp), parameter :: pi = 3.14159265358979323846d0
      real(dp) :: theta, s1, s2

      theta = bpa_sky_deg*pi/180.0d0
      s1 = sign(1.0d0, cdelt1)
      s2 = sign(1.0d0, cdelt2)
      bpa_pixel_deg = atan2(s2*cos(theta), s1*sin(theta))*180.0d0/pi
   end subroutine sky_to_pixel_bpa

   subroutine write_convolved_file(infile, outfile, naxis, sky1, sky2,&
   &freq_axis, naxes, cdelt1, cdelt2, nfreq, bmaj_in, bmin_in, bpa_in,&
   &isbad, tgt_bmaj, tgt_bmin, tgt_bpa, mem_frac_ram, io_overlap_l,&
   &nwriters, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
      use gaussft_mod, only: plan_convolution, convolve_to_beam,&
      &destroy_convolution_plan, next_fast_fft_size
      character(len=*), intent(in) :: infile, outfile
      integer, intent(in) :: naxis, sky1, sky2, freq_axis, naxes(max_axes)
      real(dp), intent(in) :: cdelt1, cdelt2
      integer, intent(in) :: nfreq
      real(dp), intent(in) :: bmaj_in(nfreq), bmin_in(nfreq), bpa_in(nfreq)
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: tgt_bmaj, tgt_bmin, tgt_bpa
      real, intent(in) :: mem_frac_ram
      logical, intent(in) :: io_overlap_l
      integer, intent(in) :: nwriters
      integer, intent(out) :: status
      integer :: nwriters_eff

      integer :: nx, ny, nx_pad, ny_pad, in_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: naxes_out(max_axes)
      ! T19: this file's own pixel-data byte offset (FTGHAD, fetched
      ! once), for write_freq_block_raw -- see the comment right before
      ! its own use, below.
      integer(kind=8) :: datastart, headstart_dum, dataend_dum
      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes, block_planes64
      integer(kind=8) :: transient_bytes_per_plane, mem_safe_bytes_io
      integer :: block_planes, chan_start, chan_len, local_iplane
      integer :: omp_threads_cap
      integer :: cur_slot
      logical :: write_dispatched_ok
      real, allocatable, target :: block_in(:,:,:), block_out(:,:,:,:)
      real(dp) :: bpa_in_pixel(nfreq), tgt_bpa_pixel
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: status_par, ich, k
      real(dp) :: dx_deg, dy_deg
      real(dp) :: nanval
      ! Per-plane timing instrumentation (T28, docs/dev/MULTI_BAND_
      ! TOMOGRAPHY_PLAN.md): one 'thread_timing stage=convolve
      ! event=done' line per plane -- finer-grained than the old
      ! per-thread-per-block lines it replaces, since the per-plane loop
      ! is serial now (see below) and each plane's own convolve_to_beam
      ! call is internally multi-threaded, not this loop running
      ! multiple iterations concurrently. Still gated behind
      ! log_level=debug like every other debug-level log_message call.
      integer :: iblock
      real(dp) :: t_thread_start, t_thread_elapsed
      character(len=160) :: thread_msg
      ! ALLOCATABLE, not automatic/stack -- one instance per OpenMP thread
      ! (private() clause below), heap-allocated (first touch inside the
      ! parallel loop) rather than stack-allocated. For a real full-size
      ! image (e.g. 4501x4501) these are ~160MB EACH in real(dp); as
      ! automatic arrays they were silently allocated on each OMP worker
      ! thread's own stack, which is typically only a few MB (OMP_STACKSIZE
      ! default) -- a guaranteed stack-overflow SIGSEGV, invisible on this
      ! project's own tiny (32x32) test fixtures but fatal on real
      ! production-scale data (found via a ~46GB end-to-end
      ! verification run). Needed only because gaussft_mod's
      ! convolve_to_beam works in real(dp), while block_in/block_out
      ! (this file's own I/O buffers) are single precision, matching every
      ! other block-I/O buffer in this project (reproject_cubes.f90's own
      ! block_data_in/out).
      real(dp), allocatable :: plane_in(:,:), plane_out(:,:)
      real(dp) :: t_stage

      status = 0
      call timer_reset_file_stages()
      call log_message('info', 'convolve', 'starting: '//trim(infile))
      nx = naxes(sky1)
      ny = naxes(sky2)
      dx_deg = abs(cdelt1)
      dy_deg = abs(cdelt2)
      do ich = 1, nfreq
         call sky_to_pixel_bpa(bpa_in(ich), cdelt1, cdelt2, bpa_in_pixel(ich))
      enddo
      call sky_to_pixel_bpa(tgt_bpa, cdelt1, cdelt2, tgt_bpa_pixel)

      fitsstat = 0
      call safe_ftopen(in_unit, trim(infile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot reopen input for output header: ', trim(infile)
         status = -1
         call free_fits_unit(in_unit)
         return
      endif

      fitsstat = 0
      ! Plain filename (NOT '!'-prefixed): FTINIT fails if outfile already
      ! exists rather than silently clobbering it. Fixed bug, found while
      ! implementing match_cubes' own skip-if-already-matched feature
      ! (planning-doc ticket) and auditing every FTINIT call for the same
      ! pattern: this call previously used the '!'-prefix CLOBBER
      ! convention, silently deleting and overwriting a pre-existing
      ! output with no warning.
      call safe_ftinit(out_unit, trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create output file (already exists?): ',&
         &trim(outfile)
         status = -1
         call safe_ftclos(in_unit, fitsstat)
         call free_fits_unit(out_unit)
         return
      endif

      naxes_out(1:naxis) = naxes(1:naxis)
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis, naxes_out(1:naxis), 0, 1, extend, fitsstat)

      call copy_generic_header_convolve(in_unit, out_unit, status)
      call FTPKYD(out_unit, 'BMAJ', tgt_bmaj/3600.0d0, 13,&
      &'common-resolution major axis FWHM (deg)', fitsstat)
      call FTPKYD(out_unit, 'BMIN', tgt_bmin/3600.0d0, 13,&
      &'common-resolution minor axis FWHM (deg)', fitsstat)
      call FTPKYD(out_unit, 'BPA', tgt_bpa, 13,&
      &'common-resolution position angle (deg)', fitsstat)
      call FTPHIS(out_unit, 'convolve_cubes: convolved from '//trim(infile)//&
      &' to a common resolution', fitsstat)
      call safe_ftclos(in_unit, fitsstat)

      ! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): CFITSIO is not
      ! guaranteed thread-safe even across different file units/handles
      ! without a specific reentrant build -- confirmed the hard way, a
      ! real, timing-dependent hang on real WALLABY+EMU data (found via
      ! match_cubes.f90's own identical design): the main thread's own
      ! next-block read blocked concurrently with the background write
      ! thread's own write, both on what is almost certainly an internal
      ! CFITSIO lock. Fetch this file's own pixel-data byte offset ONCE,
      ! right here, then close out_unit immediately -- CFITSIO's job for
      ! this file is done. Every pixel write from here on goes through
      ! write_freq_block_raw (plain Fortran stream I/O, computed byte
      ! offsets), exactly mirroring rm_synthesis_mod.f90's own
      ! io_write_threads>1 design (write_rm_chunk_raw) -- CFITSIO is
      ! never touched concurrently because it is not touched AT ALL
      ! during the write loop.
      datastart = 0_8
      headstart_dum = 0_8
      dataend_dum = 0_8
      fitsstat = 0
      call ftghad(out_unit, headstart_dum, datastart, dataend_dum, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to get data-start offset (FTGHAD) for: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif
      call safe_ftclos(out_unit, fitsstat)

      block_planes64 = 0
      nx_pad = next_fast_fft_size(nx)
      ny_pad = next_fast_fft_size(ny)
      call get_mem_total_kb(mem_total_kb)
      ! *3 not *2: block_in (1x) + block_out's own 2 double-buffer slots
      ! (io_overlap's own read+compute/write overlap, see below) --
      ! budgeted whether or not io_overlap is actually on, so turning it
      ! on never silently blows past mem_frac_ram.
      bytes_per_plane = int(4,8) * int(nx,8) * int(ny,8) * 3_8
      mem_safe_bytes = int(real(mem_frac_ram,8) * real(mem_total_kb,8) * 1024.0d0, 8)

      ! T26 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): bytes_per_plane
      ! above only budgets the block I/O buffer, not convolve_to_beam's
      ! own NaN-aware-path transient working set (a native-size NaN
      ! mask, two padded-size complex FFT buffers alive at once, two
      ! native-size real buffers, on top of this caller's own
      ! native-size plane_in, real(dp), held for the whole call),
      ! entirely un-budgeted by the block buffer above.
      ! 16*(nx*ny) + 16*(nx_pad*ny_pad) bytes. The two 16s are NOT the
      ! same "16" -- each is its own sum (see match_cubes.f90's own
      ! identical comment for the full per-array breakdown: native term
      ! = nan_mask(4) + c_d-or-c_m(4, sp) + this caller's own
      ! plane_in(8, real(dp), deliberately not converted) = 16; padded
      ! term = g_final(8, sp) + cimg-or-cmsk(8, sp) = 16). NOT computed
      ! from gaussft_mod's own type declarations, so it needs
      ! re-deriving by hand if those kinds ever change again.
      !
      ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): reserved ONCE now,
      ! not omp_get_max_threads() copies of it -- convolve_to_beam is called
      ! serially, once per plane (see the per-plane loop below), with
      ! each call internally multi-threaded via plan_convolution's own
      ! nthreads (sfftw_plan_with_nthreads) rather than N callers each
      ! running their own independent copy of this working set
      ! concurrently. So this is a per-PLANE cost now, not a per-THREAD
      ! one -- T26's own decrementing search (find the largest safe
      ! thread count) no longer applies to THIS reservation; the
      ! variable is named accordingly (transient_bytes_per_plane, not
      ! _per_thread). omp_threads_cap here is a plain performance
      ! choice (how many threads help each plane's own compute), not a
      ! memory-safety cap.
      transient_bytes_per_plane = 16_8*int(nx,8)*int(ny,8) +&
      &16_8*int(nx_pad,8)*int(ny_pad,8)
      omp_threads_cap = omp_get_max_threads()
      mem_safe_bytes_io = mem_safe_bytes - transient_bytes_per_plane
      if (mem_safe_bytes_io.le.0_8) then
         write(*,'(A,F0.2,A)') 'WARNING: convolve working memory alone (',&
         &real(transient_bytes_per_plane,8)/(1024.0d0**3),&
         &' GB/plane) exceeds the mem_frac_ram budget -- block shrunk to'//&
         &' 1 plane, but even that may not fit; raise mem_frac_ram or add'//&
         &' more RAM.'
         mem_safe_bytes_io = bytes_per_plane
      endif
      block_planes64 = max(1_8, mem_safe_bytes_io / bytes_per_plane)
      block_planes64 = min(block_planes64,&
      &max_elements_per_block / max(1_8, int(nx,8)*int(ny,8)))
      block_planes64 = max(1_8, block_planes64)
      block_planes = int(min(block_planes64, int(nfreq,8)))
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ', nfreq,&
      &' plane(s), in blocks of up to ', block_planes, ' plane(s)'

      allocate(block_in(nx, ny, block_planes))
      allocate(block_out(nx, ny, block_planes, 0:1))

      ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): nthreads here means
      ! the plan itself now executes each transform using this many
      ! threads INTERNALLY -- see plan_convolution's own comment. The
      ! caller's per-plane loop below is SERIAL (no longer an !$omp
      ! parallel do across planes), so this genuinely is the compute
      ! parallelism now, not a duplicate of it.
      call plan_convolution(nx, ny, plan_fwd, plan_bwd, nx_pad, ny_pad, omp_threads_cap)
      if (nx_pad.ne.nx .or. ny_pad.ne.ny) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)') 'Convolution FFT padded from ',&
         &nx, 'x', ny, ' to ', nx_pad, 'x', ny_pad,&
         &' (next 7-smooth size -- avoids a large-prime-factor slowdown).'
      endif

      ! io_overlap: block_out is double-buffered (trailing 0:1 slot) so
      ! the NEXT block's read+compute can proceed into the OTHER slot
      ! immediately, while THIS block's write runs on a background
      ! pthread -- but only ONE write is ever in flight on out_unit at a
      ! time (joined right here, before dispatching the next one): two
      ! concurrent CFITSIO calls on the SAME handle from different
      ! threads is unsafe (same hazard rm_synthesis_mod.f90's own
      ! io_write_threads doc comment describes for read-write handles),
      ! and this program has no raw-stream-write bypass (unlike rm_
      ! synthesis/rmclean_cubes' own io_write_threads>1 path) to avoid
      ! it. This still captures the real win: write(N) overlaps
      ! read+compute(N+1), which is where the actual dead time was
      ! (measured directly: ~44s/block write, fully serial with compute,
      ! on that dataset's own storage).
      cur_slot = 0
      write_pending = .false.
      write_failed = .false.
      status_par = 0
      ! nwriters_eff: same clamp formula as rm_synthesis/rmclean_cubes'
      ! own nwriters (T7/T12) -- never changed at rename time, kept
      ! identical here (see do_block_write's own comment for why writer
      ! threads piggyback on cores genuinely idle at write-dispatch time
      ! rather than meaningfully contending with the OMP compute pool).
      nwriters_eff = max(1, min(nwriters, omp_get_max_threads()))
      chan_start = 1
      iblock = 0
      do while (chan_start.le.nfreq)
         chan_len = min(block_planes, nfreq - chan_start + 1)
         iblock = iblock + 1

         call timer_start(t_stage)
         call read_freq_block(infile, naxis, sky1, sky2, freq_axis,&
         &naxes, chan_start, chan_len, nx, ny, block_in(:,:,1:chan_len), status_par)
         call timer_stop('block_read', t_stage)
         if (status_par.ne.0) exit

         call timer_start(t_stage)
         ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): SERIAL over
         ! planes now (no !$omp parallel do here) -- each convolve_to_beam
         ! call is internally multi-threaded instead (plan_convolution's
         ! own nthreads, above), so the compute parallelism moved INSIDE
         ! this loop's own single call per iteration rather than being
         ! this loop running N iterations concurrently. See T28's own
         ! ticket for why (T26's per-thread memory multiplication).
         ! Per-plane (not per-thread-per-block) debug timing now, since
         ! there is no longer a "per-thread" axis to this loop at all --
         ! finer-grained than the old per-block-per-thread lines, not
         ! coarser.
         do local_iplane = 1, chan_len
            ich = chan_start + local_iplane - 1
            if (isbad(ich)) then
               nanval = ieee_value(1.0_dp, ieee_quiet_nan)
               block_out(:,:,local_iplane,cur_slot) = real(nanval)
            else
               t_thread_start = omp_get_wtime()
               if (.not. allocated(plane_in)) allocate(plane_in(nx,ny))
               if (.not. allocated(plane_out)) allocate(plane_out(nx,ny))
               plane_in = real(block_in(:,:,local_iplane), dp)
               call convolve_to_beam(plan_fwd, plan_bwd,&
               &plane_in, nx, ny, nx_pad, ny_pad,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel, plane_out, k)
               block_out(:,:,local_iplane,cur_slot) = real(plane_out)
               if (k.ne.0) status_par = -1
               t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
               write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)')&
               &'thread_timing stage=convolve event=done tid=0 block=',&
               &iblock,' plane=',ich,' nthreads=',omp_threads_cap,&
               &' dur_ms=',t_thread_elapsed
               call log_message('debug','tile_thread',trim(thread_msg))
            endif
            if (status_par.ne.0) exit
         enddo
         call timer_stop('block_convolve', t_stage)
         if (status_par.ne.0) exit

         call timer_start(t_stage)
         if (write_pending) then
            call block_write_join(write_thread_id)
            write_pending = .false.
            if (write_failed) then
               status_par = -1
               exit
            endif
         endif
         call timer_stop('block_write_join', t_stage)
         write_job%file_path = outfile
         write_job%datastart = datastart
         write_job%chan_start = chan_start
         write_job%chan_len = chan_len
         write_job%nwriters_eff = max(1, min(nwriters_eff, chan_len))
         write_job%nx = nx
         write_job%ny = ny
         write_job%data => block_out(:,:,1:chan_len,cur_slot)
         call timer_start(t_stage)
         if (io_overlap_l) then
            call block_write_dispatch_async(write_job, write_thread_id,&
            &write_dispatched_ok)
            write_pending = write_dispatched_ok
         else
            call do_block_write(write_job)
            if (write_failed) then
               status_par = -1
               exit
            endif
         endif
         call timer_stop('block_write', t_stage)

         chan_start = chan_start + chan_len
         cur_slot = 1 - cur_slot
      enddo

      call timer_start(t_stage)
      if (write_pending) then
         call block_write_join(write_thread_id)
         write_pending = .false.
         if (write_failed) status_par = -1
      endif
      call timer_stop('block_write_join', t_stage)

      call destroy_convolution_plan(plan_fwd, plan_bwd)
      deallocate(block_in, block_out)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to convolve/write one or more planes for: ', trim(infile)
         status = -1
         ! out_unit is NOT open here -- closed right after FTGHAD, above,
         ! before the block-write loop ever started.
         return
      endif

      ! CASAMBM/BEAMS: always attached (see write_beams_table's own
      ! comment for why) -- convolve_cubes always tracks a real
      ! per-channel good/bad split internally, even from a source with
      ! no BEAMS table of its own (an ASCII beam log, or a plain scalar
      ! BMAJ/BMIN/BPA header), so the scalar BMAJ/BMIN/BPA written above
      ! alone would misrepresent every bad/NaN channel as sharing the
      ! common target beam. Reopen out_unit now that every block write/
      ! join above has genuinely finished (block_write_join, right
      ! before this point, guarantees no writer thread is still active).
      fitsstat = 0
      call safe_ftopen(out_unit, trim(outfile), 1, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to reopen output for BEAMS table: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         return
      endif
      call FTPKYL(out_unit, 'CASAMBM', .true.,&
      &'Multiple beams per plane (see BEAMS ext)', fitsstat)
      call write_beams_table(out_unit, nfreq, isbad, tgt_bmaj, tgt_bmin,&
      &tgt_bpa, status_par)
      if (status_par.ne.0) then
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif

      call safe_ftclos(out_unit, fitsstat)
      call log_message('info', 'convolve', 'finished: '//trim(infile))
      call timer_report_file_summary(infile)
   end subroutine write_convolved_file

   subroutine do_block_write(job)
      !! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): raw stream write,
      !! not CFITSIO -- see write_freq_block_raw's own comment. Callable
      !! either inline (io_overlap=n) or as a pthread entry point's own
      !! payload (io_overlap=y); identical logic either way, so output is
      !! bit-for-bit the same regardless of which mode dispatched it (same
      !! invariant as rm_synthesis_mod.f90/rmclean_cubes.f90's own
      !! do_tile_write). Any write failure is reported but not propagated
      !! as a return status: once dispatched onto a background thread,
      !! there is no synchronous caller left to hand a status to -- same
      !! convention rm_synthesis's own do_tile_write uses.
      type(block_write_job_t), intent(inout) :: job
      integer :: status_local, wk, wbase, wrem, wlen_k, woff_k, wchan_start_k
      logical :: any_failed

      any_failed = .false.
      if (job%nwriters_eff.le.1) then
         status_local = 0
         call write_freq_block_raw(trim(job%file_path), job%datastart,&
         &job%nx, job%ny, job%chan_start, job%chan_len, job%data, status_local)
         if (status_local.ne.0) any_failed = .true.
      else
         ! nwriters_eff>1 (T21, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
         ! split this block's own chan_len channels into nwriters_eff
         ! disjoint, contiguous chunks -- same even-split scheme
         ! (base+remainder) as rm_synthesis_mod.f90's own do_tile_write.
         ! Each chunk is written by its own write_freq_block_raw call
         ! (its own byte range, its own slice of job%data), concurrently,
         ! via a plain nested !$omp parallel do -- reusing OMP's own
         ! thread-spawning rather than managing raw pthreads directly.
         ! This composes cleanly with io_overlap: when io_overlap=y, THIS
         ! whole call already runs on one background pthread (see
         ! block_write_thread_entry below), which itself becomes the
         ! "master" of this nested parallel region.
         wbase = job%chan_len / job%nwriters_eff
         wrem = mod(job%chan_len, job%nwriters_eff)
         !$omp parallel do num_threads(job%nwriters_eff) default(none)&
         !$omp& shared(job, wbase, wrem, any_failed)&
         !$omp& private(wk, wlen_k, woff_k, wchan_start_k, status_local)
         do wk = 0, job%nwriters_eff-1
            wlen_k = wbase + merge(1, 0, wk.lt.wrem)
            woff_k = wk*wbase + min(wk, wrem)
            wchan_start_k = job%chan_start + woff_k
            status_local = 0
            if (wlen_k.gt.0) then
               call write_freq_block_raw(trim(job%file_path), job%datastart,&
               &job%nx, job%ny, wchan_start_k, wlen_k,&
               &job%data(:,:,woff_k+1:woff_k+wlen_k), status_local)
               if (status_local.ne.0) then
                  !$omp atomic write
                  any_failed = .true.
               endif
            endif
         enddo
         !$omp end parallel do
      endif
      if (any_failed) then
         write(*,*) 'ERROR: background write failed for channels ',&
         &job%chan_start, '-', job%chan_start+job%chan_len-1
         write_failed = .true.
      endif
   end subroutine do_block_write

   function block_write_thread_entry(arg) bind(C) result(res)
      !! pthread start routine. Unpacks the opaque context pointer back
      !! into the block_write_job_t it was created from (same process,
      !! same build, so the round-trip through c_loc/c_f_pointer is safe
      !! even though block_write_job_t is not a bind(C) type).
      type(c_ptr), value :: arg
      type(c_ptr) :: res
      type(block_write_job_t), pointer :: job

      call c_f_pointer(arg, job)
      call do_block_write(job)
      res = c_null_ptr
   end function block_write_thread_entry

   subroutine block_write_dispatch_async(job, thread_id, dispatched)
      !! Launches do_block_write(job) on a background pthread. `job` must
      !! have the TARGET attribute at the call site (write_job's own
      !! declaration) and must remain valid -- untouched and
      !! undeallocated -- until block_write_join(thread_id) has returned.
      !! On pthread_create failure this runs the write synchronously right
      !! here instead (safe fallback: a write is never silently dropped),
      !! and reports dispatched=.false. so the caller knows there is
      !! nothing to join later.
      type(block_write_job_t), intent(inout), target :: job
      integer(c_long), intent(out) :: thread_id
      logical, intent(out) :: dispatched
      integer(c_int) :: rc

      rc = c_pthread_create(thread_id, c_null_ptr,&
      &c_funloc(block_write_thread_entry), c_loc(job))
      if (rc.ne.0_c_int) then
         write(*,*) 'WARNING: pthread_create failed for async block write;'//&
         &' running inline'
         call do_block_write(job)
         dispatched = .false.
      else
         dispatched = .true.
      endif
   end subroutine block_write_dispatch_async

   subroutine block_write_join(thread_id)
      integer(c_long), intent(in) :: thread_id
      integer(c_int) :: rc
      rc = c_pthread_join(thread_id, c_null_ptr)
   end subroutine block_write_join

   subroutine write_beams_table(unit, nfreq, isbad, tgt_bmaj, tgt_bmin,&
   &tgt_bpa, status)
      !! Appends a CASA-style BEAMS binary table extension (EXTNAME=
      !! 'BEAMS', columns BMAJ/BMIN/BPA/CHAN/POL -- matching the real
      !! ASKAP cube's own 5-column layout this file's top comment
      !! already documents: BMAJ/BMIN in arcsec, BPA in deg, CHAN
      !! 0-indexed, POL always 0 since convolve_cubes processes one
      !! Stokes product -- one file -- at a time) to the current
      !! (primary) HDU of unit, one row per channel. A GOOD channel gets
      !! the common target beam every good plane was actually convolved
      !! to; a BAD (skipped, all-NaN) channel gets the same degenerate
      !! sentinel (tiny(1.0), CASA's own ~1.18e-38 placeholder -- see
      !! read_beams_table's own comment) this program's own readers
      !! already treat as "no valid beam", so a downstream reader sees
      !! exactly which channels were actually convolved, rather than a
      !! single scalar BMAJ/BMIN/BPA that would claim every channel
      !! (including bad/NaN ones) shares the common beam.
      integer, intent(inout) :: unit
      integer, intent(in) :: nfreq
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: tgt_bmaj, tgt_bmin, tgt_bpa
      integer, intent(out) :: status

      character(len=8) :: ttype(5), tform(5), tunit_(5)
      real, allocatable :: col_bmaj(:), col_bmin(:), col_bpa(:)
      integer, allocatable :: col_chan(:), col_pol(:)
      integer :: fitsstat, ich, colnum

      ttype = (/'BMAJ    ', 'BMIN    ', 'BPA     ', 'CHAN    ', 'POL     '/)
      tform = (/'1E      ', '1E      ', '1E      ', '1J      ', '1J      '/)
      tunit_ = (/'arcsec  ', 'arcsec  ', 'deg     ', '        ', '        '/)

      allocate(col_bmaj(nfreq), col_bmin(nfreq), col_bpa(nfreq),&
      &col_chan(nfreq), col_pol(nfreq))
      do ich = 1, nfreq
         col_chan(ich) = ich - 1
         col_pol(ich) = 0
         if (isbad(ich)) then
            col_bmaj(ich) = tiny(1.0)
            col_bmin(ich) = tiny(1.0)
            col_bpa(ich) = 0.0
         else
            col_bmaj(ich) = real(tgt_bmaj)
            col_bmin(ich) = real(tgt_bmin)
            col_bpa(ich) = real(tgt_bpa)
         endif
      enddo

      fitsstat = 0
      call FTIBIN(unit, nfreq, 5, ttype, tform, tunit_, 'BEAMS', 0, fitsstat)
      call FTGCNO(unit, .false., 'BMAJ', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bmaj, fitsstat)
      call FTGCNO(unit, .false., 'BMIN', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bmin, fitsstat)
      call FTGCNO(unit, .false., 'BPA', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bpa, fitsstat)
      call FTGCNO(unit, .false., 'CHAN', colnum, fitsstat)
      call FTPCLJ(unit, colnum, 1, 1, nfreq, col_chan, fitsstat)
      call FTGCNO(unit, .false., 'POL', colnum, fitsstat)
      call FTPCLJ(unit, colnum, 1, 1, nfreq, col_pol, fitsstat)
      deallocate(col_bmaj, col_bmin, col_bpa, col_chan, col_pol)

      status = 0
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to write BEAMS binary table extension'
         status = -1
      endif
   end subroutine write_beams_table

   subroutine copy_generic_header_convolve(src_unit, dst_unit, status)
      !! No reprojection happens here (unlike reproject_cubes) -- input
      !! and output share IDENTICAL axis layout and numbering, so every
      !! header card copies through verbatim EXCEPT the structural
      !! keywords FTPHPR already wrote, and BMAJ/BMIN/BPA/CASAMBM (the
      !! caller overwrites BMAJ/BMIN/BPA with the new common beam right
      !! after calling this, and writes its own CASAMBM=T plus a freshly
      !! synthesized BEAMS table -- see write_beams_table -- once the
      !! per-channel good/bad split is known; a stale CASAMBM=T/BEAMS
      !! copied verbatim from the input here, describing the INPUT's own
      !! pre-convolution per-channel beams, would be actively misleading
      !! once every good channel has actually been convolved to one
      !! common target beam).
      integer, intent(in) :: src_unit, dst_unit
      integer, intent(inout) :: status
      integer :: nkeys, nmore, i, fitsstat
      character(len=80) :: card
      character(len=8) :: key

      if (status.ne.0) return
      fitsstat = 0
      call FTGHSP(src_unit, nkeys, nmore, fitsstat)
      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(src_unit, i, card, fitsstat)
         if (fitsstat.ne.0) cycle
         key = adjustl(card(1:8))
         select case (trim(key))
         case ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT',&
         &'END', 'BMAJ', 'BMIN', 'BPA', 'CASAMBM')
            cycle
         end select
         if (is_naxis_keyword(key)) cycle
         fitsstat = 0
         call FTPREC(dst_unit, card, fitsstat)
      enddo
   end subroutine copy_generic_header_convolve

   logical function is_naxis_keyword(key)
      character(len=8), intent(in) :: key
      integer :: i, klen
      is_naxis_keyword = .false.
      klen = len_trim(key)
      if (klen.le.5) return
      if (key(1:5).ne.'NAXIS') return
      do i = 6, klen
         if (key(i:i).lt.'0' .or. key(i:i).gt.'9') return
      enddo
      is_naxis_keyword = .true.
   end function is_naxis_keyword

   subroutine read_freq_block(filename, naxis, sky1, sky2, freq_axis,&
   &naxes, chan_start, chan_len, nx, ny, block_data, status)
      !! Single-varying-axis specialisation of reproject_cubes.f90's own
      !! read_one_block (see that subroutine's own comment for the full
      !! FTGSVE axis-order-permute reasoning this reuses directly): reads
      !! chan_len consecutive planes of the FREQ axis, full sky extent,
      !! every axis besides sky1/sky2/freq_axis already verified
      !! degenerate (extent 1) by read_axis_info, so there is no "other
      !! group" loop needed here at all, unlike the fully general
      !! N-non-sky-axis case reproject_cubes itself has to handle.
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: naxis, sky1, sky2, freq_axis, naxes(max_axes)
      integer, intent(in) :: chan_start, chan_len, nx, ny
      real, intent(out) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, group
      integer :: fpixels(max_axes), lpixels(max_axes), incs(max_axes)
      logical :: anyflg
      real :: badval
      integer :: rank_sky1, rank_sky2, rank_freq
      integer :: dims(3), idxvec(3), i, j, c
      logical :: natural_order
      real, allocatable :: natural_buf(:,:,:)

      if (status.ne.0) return

      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(sky1) = nx
      lpixels(sky2) = ny
      fpixels(freq_axis) = chan_start
      lpixels(freq_axis) = chan_start + chan_len - 1

      rank_sky1 = 1
      if (sky2.lt.sky1) rank_sky1 = rank_sky1 + 1
      if (freq_axis.lt.sky1) rank_sky1 = rank_sky1 + 1
      rank_sky2 = 1
      if (sky1.lt.sky2) rank_sky2 = rank_sky2 + 1
      if (freq_axis.lt.sky2) rank_sky2 = rank_sky2 + 1
      rank_freq = 1
      if (sky1.lt.freq_axis) rank_freq = rank_freq + 1
      if (sky2.lt.freq_axis) rank_freq = rank_freq + 1
      natural_order = (rank_sky1.eq.1 .and. rank_sky2.eq.2 .and. rank_freq.eq.3)

      fitsstat = 0
      blocksize = 1
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (natural_order) then
         call FTGSVE(unit, group, naxis, naxes(1:naxis), fpixels(1:naxis),&
         &lpixels(1:naxis), incs(1:naxis), badval, block_data, anyflg, fitsstat)
      else
         dims(rank_sky1) = nx
         dims(rank_sky2) = ny
         dims(rank_freq) = chan_len
         allocate(natural_buf(dims(1), dims(2), dims(3)))
         call FTGSVE(unit, group, naxis, naxes(1:naxis), fpixels(1:naxis),&
         &lpixels(1:naxis), incs(1:naxis), badval, natural_buf, anyflg, fitsstat)
         if (fitsstat.eq.0) then
            do c = 1, chan_len
               do j = 1, ny
                  do i = 1, nx
                     idxvec(rank_sky1) = i
                     idxvec(rank_sky2) = j
                     idxvec(rank_freq) = c
                     block_data(i,j,c) = natural_buf(idxvec(1), idxvec(2), idxvec(3))
                  enddo
               enddo
            enddo
         endif
         deallocate(natural_buf)
      endif
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') from ', trim(filename)
         status = -1
      endif
   end subroutine read_freq_block

   logical function host_is_big_endian_cv() result(is_be)
      !! Verbatim port of rm_synthesis_mod.f90's own host_is_big_endian.
      integer(int32) :: probe
      integer(int8) :: bytes(4)

      probe = 1_int32
      bytes = transfer(probe, bytes)
      is_be = (bytes(1) == 0_int8)
   end function host_is_big_endian_cv

   subroutine swap_bytes_r4_inplace_cv(buf, n)
      !! Verbatim port of rm_synthesis_mod.f90's own
      !! swap_bytes_r4_inplace.
      integer(kind=int64), intent(in) :: n
      real, intent(inout) :: buf(n)
      integer(kind=int64) :: i
      integer(int8) :: b(4), t

      do i = 1_int64, n
         b = transfer(buf(i), b)
         t = b(1); b(1) = b(4); b(4) = t
         t = b(2); b(2) = b(3); b(3) = t
         buf(i) = transfer(b, buf(i))
      end do
   end subroutine swap_bytes_r4_inplace_cv

   subroutine write_freq_block_raw(file_path, datastart, nx, ny, chan_start,&
   &chan_len, block_data, status)
      !! Raw Fortran stream I/O at a computed byte offset -- CFITSIO is
      !! never touched here at all (T19, docs/dev/MULTI_BAND_
      !! TOMOGRAPHY_PLAN.md -- see write_convolved_file's own comment
      !! for the full rationale: a real, timing-dependent hang, found
      !! via match_cubes.f90's identical design on real WALLABY+EMU
      !! data). Verbatim design port of rm_synthesis_mod.f90's own
      !! write_rm_chunk_raw. Output axis layout is IDENTICAL to the
      !! input's own here (no reprojection, unlike reproject_cubes'
      !! write_one_block) and this tool has no spatial sub-tiling at
      !! all -- so every channel plane is always exactly nx*ny
      !! contiguous elements, planes back-to-back in channel order
      !! (rm_synthesis_mod.f90's own "full-width" case, never its
      !! "partial-width" one).
      !!
      !! Only ever one writer in flight at a time by this whole file's
      !! own design (block_write_join always runs before the next
      !! dispatch) -- unlike rm_synthesis' own io_write_threads>1 (N
      !! genuinely concurrent writers), so no critical section is
      !! needed around the newunit= open here.
      character(len=*), intent(in) :: file_path
      integer(kind=8), intent(in) :: datastart
      integer, intent(in) :: nx, ny, chan_start, chan_len
      real, intent(in) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: u, ios, ip
      logical :: need_swap
      integer(kind=8) :: plane_stride_bytes, plane_elems, byte_pos
      real, allocatable :: plane_buf(:,:)

      if (status.ne.0) return

      need_swap = .not. host_is_big_endian_cv()
      plane_elems = int(nx,8) * int(ny,8)
      plane_stride_bytes = plane_elems * 4_8

      open(newunit=u, file=trim(file_path), access='stream',&
      &form='unformatted', status='old', action='write', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: write_freq_block_raw: failed to open ', trim(file_path)
         status = -1
         return
      endif

      allocate(plane_buf(nx,ny))
      do ip = 1, chan_len
         plane_buf = block_data(:,:,ip)
         if (need_swap) call swap_bytes_r4_inplace_cv(plane_buf, plane_elems)
         byte_pos = datastart + int(chan_start-1+ip-1,8)*plane_stride_bytes + 1_8
         write(u, pos=byte_pos, iostat=ios) plane_buf
         if (ios.ne.0) then
            write(*,*) 'ERROR: write_freq_block_raw: write failed (plane ',&
            &chan_start+ip-1, ') for ', trim(file_path)
            status = -1
            exit
         endif
      enddo
      deallocate(plane_buf)
      close(u)
   end subroutine write_freq_block_raw

   subroutine get_mem_total_kb(mem_total_kb)
      integer(kind=8), intent(out) :: mem_total_kb
      integer :: mem_unit, ios_mem
      character(len=128) :: mem_line
      integer(kind=8) :: mem_kb_tmp

      mem_total_kb = 0_8
      open(newunit=mem_unit, file='/proc/meminfo', status='old', iostat=ios_mem)
      if (ios_mem.eq.0) then
         do
            read(mem_unit, '(A)', iostat=ios_mem) mem_line
            if (ios_mem.ne.0) exit
            if (index(mem_line, 'MemTotal:').eq.1) then
               read(mem_line(index(mem_line,':')+1:), *, iostat=ios_mem) mem_kb_tmp
               if (ios_mem.eq.0) mem_total_kb = mem_kb_tmp
               exit
            endif
         enddo
         close(mem_unit)
      endif
      if (mem_total_kb.le.0_8) mem_total_kb = 4194304_8
   end subroutine get_mem_total_kb

   subroutine run_dry_run_check(target_path, tool_name)
      !! dry_run=y (T24, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): planning-
      !! only pass, mirroring rm_synthesis' own dry_run/tile_autotune.cfg
      !! convention -- touches no real data, just checks the target output
      !! disk's own rotational status (via /sys/block/<dev>/queue/
      !! rotational, resolved from target_path's own mount point) and
      !! writes a suggested <tool_name>_dryrun.cfg with recommended
      !! io_overlap/nwriters settings. Advisory only: this project's own
      !! nwriters clamp formula was deliberately left tied to
      !! omp_get_max_threads(), not disk type (see T7/T12's own
      !! discussion) -- this dry-run pass suggests a starting VALUE within
      !! that clamp, it does not change the clamp itself. See match_cubes.
      !! f90's own copy of this subroutine for the full rationale (adapt,
      !! don't share, per this project's own module convention).
      character(len=*), intent(in) :: target_path, tool_name

      character(len=600) :: shell_cmd
      character(len=64) :: capture_file
      character(len=16) :: rot_str
      integer :: capture_unit, ios_cap, exitstat_cmd, cfg_unit
      integer :: rot_val
      logical :: is_rotational, have_answer

      write(capture_file, '(A,I0,A)') '.dryrun_rotational_', getpid_wrap(), '.tmp'

      write(shell_cmd, '(A,A,A,A,A)')&
      &'dev=$(df --output=source ''', trim(target_path), ''' 2>/dev/null | tail -1); ',&
      &'base=$(basename "$dev"); ',&
      &'if echo "$base" | grep -Eq ''^nvme[0-9]+n[0-9]+p[0-9]+$''; then parent="${base%p*}"; '//&
      &'elif echo "$base" | grep -Eq ''^nvme[0-9]+n[0-9]+$''; then parent="$base"; '//&
      &'else parent=$(echo "$base" | sed -E ''s/[0-9]+$//''); fi; '//&
      &'cat "/sys/block/$parent/queue/rotational" 2>/dev/null > '//trim(capture_file)

      call execute_command_line(trim(shell_cmd), wait=.true., exitstat=exitstat_cmd)

      have_answer = .false.
      rot_val = -1
      open(newunit=capture_unit, file=trim(capture_file), status='old',&
      &action='read', iostat=ios_cap)
      if (ios_cap.eq.0) then
         read(capture_unit, '(A)', iostat=ios_cap) rot_str
         close(capture_unit, status='delete')
         if (ios_cap.eq.0) then
            read(rot_str, *, iostat=ios_cap) rot_val
            if (ios_cap.eq.0 .and. (rot_val.eq.0 .or. rot_val.eq.1)) have_answer = .true.
         endif
      else
         close(capture_unit, status='delete', iostat=ios_cap)
      endif

      is_rotational = (rot_val.eq.1)

      write(*,'(A)') 'dry_run=y: no files will be processed.'
      write(*,'(A,A)') 'Target path checked: ', trim(target_path)
      if (.not. have_answer) then
         write(*,'(A)') 'Could not determine the target disk''s rotational status'//&
         &' (df/sysfs check failed or path not found) -- no suggestion made.'
         return
      endif
      if (is_rotational) then
         write(*,'(A)') 'Target disk: spinning (rotational=1).'
         write(*,'(A)') 'Suggestion: io_overlap=n, nwriters=1 -- on a single spinning'//&
         &' disk, concurrent reads/writes contend for the same physical head;'
         write(*,'(A)') '  see docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md T24 for the full'//&
         &' reasoning (this is advisory, not enforced).'
      else
         write(*,'(A)') 'Target disk: non-rotational (SSD/NVMe, rotational=0).'
         write(*,'(A)') 'Suggestion: io_overlap=y, nwriters=2 -- concurrent I/O is'//&
         &' generally safe on SSD/NVMe; nwriters starts conservative (this machine''s'
         write(*,'(A)') '  own spare-core headroom outside OMP_NUM_THREADS), raise it'//&
         &' if your own machine has more idle cores available.'
      endif

      open(newunit=cfg_unit, file=trim(tool_name)//'_dryrun.cfg', status='replace',&
      &action='write', iostat=ios_cap)
      if (ios_cap.ne.0) then
         write(*,'(A)') 'WARNING: could not write suggested cfg file.'
         return
      endif
      write(cfg_unit, '(A)') '# Autogenerated by '//trim(tool_name)//' dry_run=y'
      write(cfg_unit, '(A)') '# Target path checked: '//trim(target_path)
      if (is_rotational) then
         write(cfg_unit, '(A)') '# Target disk: spinning (rotational=1)'
         write(cfg_unit, '(A)') 'io_overlap=n'
         write(cfg_unit, '(A)') 'nwriters=1'
      else
         write(cfg_unit, '(A)') '# Target disk: non-rotational (SSD/NVMe, rotational=0)'
         write(cfg_unit, '(A)') 'io_overlap=y'
         write(cfg_unit, '(A)') 'nwriters=2'
      endif
      write(cfg_unit, '(A)') '# Copy these KEY=VALUE lines into your own cfg, or pass'
      write(cfg_unit, '(A)') '# them directly on the command line.'
      close(cfg_unit)
      write(*,'(A)') 'Wrote suggested settings to '//trim(tool_name)//'_dryrun.cfg'
   end subroutine run_dry_run_check

   integer function getpid_wrap() result(pid)
      !! Wraps the getpid() C library call (via iso_c_binding) purely to
      !! make the dry_run capture-file name collision-proof if multiple
      !! dry_run invocations somehow overlap in the same directory --
      !! not needed for correctness of a single run, just cheap safety.
      use, intrinsic :: iso_c_binding, only: c_int
      interface
         function c_getpid() bind(C, name="getpid") result(r)
            import :: c_int
            integer(c_int) :: r
         end function c_getpid
      end interface
      pid = int(c_getpid())
   end function getpid_wrap

end program convolve_cubes
