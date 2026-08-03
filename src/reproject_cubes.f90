! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! FITS cubes -- two or more, not tied to the multi-band-tomography "band"
! concept specifically -- onto a single reference cube's grid using
! Starlink AST for WCS handling and resampling, so the existing
! rm_synthesis exact-match ingestion (T1) can consume genuinely misaligned
! bands unchanged.
!
! Current stage: N-input footprint-mode computation (intersection/union/
! reference), built on top of two previously-verified building blocks:
! extract_sky_mapping (automatic sky-axis detection, correct for any axis
! order/adjacency -- see its own comment for why the earlier consecutive-
! axis-only search was wrong) and astMapBox-based footprint bounds
! (verified against a genuine partial-overlap case, not just full
! overlap). Cross-file alignment composes each file's own pixel->sky
! Mapping with an astConvert between the two SkyFrame *objects* (not
! whole FrameSets -- astConvert's domain search does not recurse into a
! CmpFrame's internal components, confirmed against the actual SUN/211
! manual, so it cannot align two compound "STOKES-SKY-SPECTRUM" current
! frames directly).
!
! Footprint policy: zero overlap between any input and the running output
! grid is always a hard failure, regardless of mode (loudly refuse before
! compute, matching the rest of this project's philosophy) -- this is
! different from *partial* overlap, which is the normal, expected case
! for combining bands with different sky coverage and is not rejected.
!
! Also demonstrates the actual regridding step: reads planes and
! resamples them onto the final output grid via astResampleR (originally
! one plane at a time, via a since-removed resample_one_plane -- see the
! OpenMP/block paragraph below for its replacement). Caught a real bug on
! the first attempt: the documented
! 20-argument signature (this, ndim_in, lbnd_in, ubnd_in, in, in_var,
! interp, finterp, params, flags, tol, maxpix, badval, ndim_out,
! lbnd_out, ubnd_out, lbnd, ubnd, out, out_var) was missing "params"
! entirely, silently shifting every later argument by one position (a
! REAL array landing where an INTEGER "flags" scalar was expected, etc)
! -- segfaulted rather than erroring cleanly, so this needed bisecting
! with diagnostic prints rather than reading a status code. Verified
! after the fix: resampling the axis-swapped fixture onto the reference
! grid reproduces the exact known ground-truth value at a known pixel
! (read independently via Python beforehand), and the union-mode test's
! NaN count for uncovered output pixels matches the expected uncovered
! area exactly (16x32 = 512).
!
! Loops over every combination of an input's non-sky axes (channel,
! Stokes, or any other axis a cube happens to have), decoding a flat
! plane index via a mixed-radix counter so it works for any number of
! non-sky axes, not just the 2 (channel, Stokes) this project's own
! fixtures happen to have -- and reuses the same input->reference Mapping
! across every plane (it depends only on the 2 sky axes, never on which
! plane is being read). Verified across all 200 channels of the axis-
! swapped fixture: 3 spot-checked channels (1, 100, 200) all reproduce
! their known ground-truth values exactly, and the whole 200-channel run
! completes in well under a second.
!
! Actually writes an output FITS file now (write_reprojected_file):
! output axis layout puts the 2 sky axes at OUTPUT positions 1,2 (RA-
! fastest-on-disk, matching what rm_synthesis's own tile-read I/O already
! assumes -- confirmed against rm_synthesis_mod.f90's own auto-tile-
! planner comment), in whatever order the REFERENCE presents them; other
! (non-sky) axes keep INFILE's own values unchanged, at output positions
! 3... This guarantees "sky is axes 1,2" but only guarantees "axis 1 is
! literally RA" for a conventionally-ordered reference (RA before Dec) --
! a deliberate, documented scope limit, not a silent assumption. CRPIX
! for the sky axes is shifted for the output grid's own origin using the
! exact same formula rm_synthesis.f90 already uses for its own subimage
! CRPIX shift. Verified independently via Python/astropy, header and
! data both: reference mode reproduces the reference's own header
! exactly and all 3 spot-checked channels match known ground truth;
! intersection mode's CRPIX1 correctly shifts from 17.0 to 1.0 for a
! 16-pixel crop, CRPIX2 stays untouched (that axis wasn't cropped), and
! the pixel value at the shifted position matches independently-computed
! ground truth exactly.
!
! OpenMP-parallelised across planes, and reads/resamples/writes in
! memory-budgeted BLOCKS of planes (mem_frac_ram, same concept as
! rm_synthesis's own tile planner -- see get_mem_total_kb) rather than
! one plane at a time. Each block goes through three strictly separated
! phases: read (one thread, OMP `single`), resample (every thread, in
! parallel -- each thread builds its own private AST Mapping, since this
! Fortran AST binding has no astLock_/astUnlock_ to hand one Mapping
! between threads), write (one thread, `single`). Blocks mean a whole
! block's CFITSIO I/O happens on a single thread BY CONSTRUCTION, so
! there is nothing to lock -- an earlier plane-at-a-time version
! serialised every CFITSIO call behind an OMP critical section instead,
! which this replaced (see write_reprojected_file's own comment). A real
! bug surfaced getting here: FTGSVE fills its output array in ascending-
! axis-number order among the axes actually being read, and reading a
! whole block (not a single degenerate plane) exposed a case the old
! per-plane code never hit -- TEST_NONADJACENT.Q.FITSCUBE has FREQ on
! axis 1, *before* RA/DEC on axes 2 and 4, so a block read's fastest
! dimension is the block axis, not the sky axes; read_one_block now
! ranks the 3 relevant axes and only pays for a permute copy on
! non-conventional orderings like that one, not the common case.
! Verified byte-identical (header and data) against the pre-blocking,
! pre-OpenMP committed output across all fixtures, multiple block sizes
! (including a 3-plane block against 200 channels, exercising a
! non-exact remainder), and 25 repeated stress runs at default thread
! count with no failures.
!
! Usage: reproject_cubes mode=<intersection|union|reference> reffile=<reference_file> infiles=<input_file>[,<input_file>...]
!    or: reproject_cubes --config <cfgfile>
!    or: reproject_cubes --config <cfgfile> mode=<...> [reffile=<...>] [infiles=<...>]
!    or: reproject_cubes --help | -h
! No positional args: mode/reffile/infiles are always named key=value
! (no spaces around '='), on the command line or via --config, never
! inferred from argument order -- deliberate, to leave no room for a
! user mistake on which bare word means what. Each CLI key=value token
! overrides only that same key from --config (per-key precedence, not
! an atomic-group replacement -- unambiguous now that every CLI value
! names its own field). Config file is the same key=value style:
! mode=..., reffile=..., and infiles=file1,file2,file3 (comma-
! separated, same csv-list convention rm_synthesis's own multi-band
! config keys use). Full usage text is in print_usage below (shared by
! --help and the argument-error path).
program reproject_cubes
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_ptr, c_funptr,&
   &c_null_ptr, c_funloc, c_loc, c_f_pointer
   use logging_mod
   use fitsio_unit_mod
   implicit none

   ! --- Logging & timing (planning-doc ticket) -- see convolve_cubes.f90
   ! /logging_mod.f90's own header comments. Stage names used below:
   ! startup, block_read, block_resample, block_write, block_write_join,
   ! finalize.
   character(len=16) :: log_level
   logical :: timing_enabled
   character(len=272) :: log_output_file
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 (`*`-column comments) and cannot be `include`d
   ! into a free-form .f90 file directly (gfortran misparses its comments
   ! as code). Only the handful of symbols actually used are declared
   ! directly instead, matching AST_PAR's own declared types exactly
   ! (checked against /usr/include/AST_PAR).
   external :: ast_null
   integer, parameter :: ast__null = 0
   integer, parameter :: ast__szchr = 200
   integer, parameter :: ast__base = 0
   integer, parameter :: ast__current = -1
   integer, external :: ast_fitschan, ast_read, ast_geti
   integer, external :: ast_getmapping, ast_simplify, ast_getframe
   integer, external :: ast_pickaxes, ast_cmpmap, ast_convert
   integer, external :: ast_resampler
   integer, parameter :: ast__linear = 5
   logical, external :: ast_isaframeset, ast_isaskyframe
   character(len=ast__szchr), external :: ast_getc

   integer, parameter :: max_axes = 10
   integer, parameter :: max_inputs = 50
   ! Hard ceiling on elements-per-block -- see convolve_cubes.f90's own
   ! identical constant/comment (CFITSIO's Fortran wrapper takes a
   ! default-INTEGER element count; this clamp keeps
   ! max(nx_in*ny_in,nx_out*ny_out)*block_planes safely under 2^31-1 for
   ! ANY image size/mem_frac_ram combination).
   integer(kind=8), parameter :: max_elements_per_block = 2000000000_8

   ! --- io_overlap: background-thread block write (planning-doc ticket)
   ! --- same design as convolve_cubes.f90's own (see write_convolved_
   ! file's own comment there for the full single-writer-at-a-time
   ! rationale). Dispatched from inside write_reprojected_file's own
   ! `!$omp single` write phase (see there): the single thread kicks off
   ! the pthread and returns immediately, so the `single` construct's own
   ! implicit barrier only waits for the (fast) dispatch call, not the
   ! full write -- the NEXT block's read+resample then proceeds
   ! concurrently with the write running on the background thread.
   type :: block_write_job_t
      integer :: out_unit = 0
      integer :: naxes_in(max_axes) = 0, other_axes(max_axes) = 0
      integer :: other_idx(max_axes) = 0, n_other = 0
      integer :: chan_start = 0, chan_len = 0, nx = 0, ny = 0
      real, pointer :: data(:,:,:) => null()
   end type block_write_job_t
   type(block_write_job_t), target, save :: write_job
   integer(c_long) :: write_thread_id = 0
   logical :: write_pending = .false.
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

   character(len=16) :: mode
   character(len=512) :: reffile
   character(len=512) :: infiles(max_inputs)
   integer :: n_inputs, i

   character(len=512) :: cfgfile
   ! this_arg/cli_val/raw_cli_infiles hold a WHOLE comma-separated
   ! infiles= argument, not one path -- up to max_inputs (50) entries,
   ! each up to 512 chars (infiles(:)'s own per-entry length). 16384
   ! gives generous headroom over the worst case (~25,600 chars); 512
   ! here silently truncated a real run's own infiles= argument (found
   ! via scripts/run_pipeline.sh's symlink-redirection scheme, which
   ! lengthens every path enough to cross 512 for just 4 real files).
   character(len=16384) :: this_arg, cli_val
   character(len=512) :: cli_key
   character(len=16) :: cli_mode
   character(len=512) :: cli_reffile
   character(len=16384) :: raw_cli_infiles
   integer :: argc, iarg
   logical :: has_kv
   logical :: have_cfgfile
   logical :: seen_mode, seen_reffile, seen_infiles
   logical :: cli_seen_mode, cli_seen_reffile, cli_seen_infiles

   ! mem_frac_ram: optional (unlike mode/reffile/infiles, has a default
   ! and is never required) -- fraction of total system RAM (see
   ! get_mem_total_kb) budgeted for one read/resample/write block of
   ! planes at a time, same concept and default as rm_synthesis's own
   ! cfg%mem_frac_ram (rm_synthesis_mod.f90's plan_tile).
   real :: mem_frac_ram, cli_mem_frac_ram
   logical :: cli_seen_mem_frac_ram
   integer :: ios_mfr
   ! io_overlap (default n): background-thread block write, overlapped
   ! with the NEXT block's read+resample -- same scheme/key name as rm_
   ! synthesis/rmclean_cubes/convolve_cubes' own io_overlap (planning-doc
   ! ticket, added after the real ~46GB Jennifer end-to-end run measured
   ! this class of block write as fully serial with compute on real
   ! storage). Same CLI-then-cfg-then-override merge pattern as
   ! mem_frac_ram just above.
   logical :: io_overlap, cli_io_overlap
   logical :: cli_seen_io_overlap
   ! Logging & timing (planning-doc ticket) -- log_level/timing_enabled/
   ! log_output_file themselves declared near the top of this program
   ! (alongside the `use logging_mod`); only their CLI-override
   ! counterparts are declared here, same CLI-then-cfg-then-override
   ! merge pattern as mem_frac_ram/io_overlap above.
   character(len=16) :: cli_log_level
   logical :: cli_timing_enabled
   character(len=272) :: cli_log_output_file
   logical :: cli_seen_log_level, cli_seen_timing_enabled
   logical :: cli_seen_log_output_file

   integer :: wcs_ref, skymap_ref, skyframe_ref
   integer :: naxes_ref(max_axes), pixaxes_ref(2)
   integer :: wcs_in, skymap_in, skyframe_in
   integer :: naxes_in(max_axes), pixaxes_in(2)
   integer :: map_in2ref

   ! --- Skip-if-already-matched (planning-doc ticket, same logic as
   ! match_cubes.f90's own -- see there for the full rationale) ---
   logical :: needs_processing(max_inputs), out_exists
   integer :: nx_out_common, ny_out_common
   integer :: ref_skip_unit, cand_skip_unit, fitsstat_skip, blocksize_skip
   double precision :: lbnd_out(2), ubnd_out(2)
   double precision :: this_lbnd(2), this_ubnd(2)
   integer :: status

   ! Two ways to supply mode/reffile/infiles: explicit key=value CLI
   ! tokens (mode=..., reffile=..., infiles=...) or a --config key=value
   ! file (read_reproject_cfg below, same csv-list convention). No
   ! positional args -- deliberately: a bare "reproject_cubes
   ! intersection ref.fits a.fits b.fits" leaves the user to remember
   ! argument order from memory, which is exactly the kind of mistake
   ! named key=value args rule out. Both sources may be given together;
   ! each CLI key=value token overrides only its own field from --config
   ! (per-key, not atomic-group, precedence -- unambiguous now that every
   ! CLI value names the field it sets).
   have_cfgfile = .false.
   cli_seen_mode = .false.
   cli_seen_reffile = .false.
   cli_seen_infiles = .false.
   cli_seen_mem_frac_ram = .false.
   mem_frac_ram = 0.25
   cli_seen_io_overlap = .false.
   io_overlap = .false.
   cli_seen_log_level = .false.
   log_level = 'info'
   cli_seen_timing_enabled = .false.
   timing_enabled = .false.
   cli_seen_log_output_file = .false.
   log_output_file = ' '
   argc = command_argument_count()
   iarg = 1
   do while (iarg.le.argc)
      call get_command_argument(iarg, this_arg)
      if (trim(this_arg).eq.'--help' .or. trim(this_arg).eq.'-h') then
         call print_usage()
         stop
      else if (trim(this_arg).eq.'--config') then
         if (iarg.eq.argc) then
            write(*,*) 'ERROR: --config requires a file path argument'
            stop 1
         endif
         call get_command_argument(iarg+1, cfgfile)
         have_cfgfile = .true.
         iarg = iarg + 2
      else
         call split_cli_kv(this_arg, cli_key, cli_val, has_kv)
         if (.not. has_kv) then
            write(*,*) 'ERROR: unrecognised argument "', trim(this_arg),&
            &'" -- expected key=value (mode=..., reffile=..., infiles=...),',&
            &' --config <file>, or --help'
            stop 1
         endif
         select case (trim(cli_key))
         case ('mode')
            if (cli_seen_mode) then
               write(*,*) 'ERROR: mode given more than once on the command line'
               stop 1
            endif
            cli_mode = trim(cli_val)
            cli_seen_mode = .true.
         case ('reffile')
            if (cli_seen_reffile) then
               write(*,*) 'ERROR: reffile given more than once on the command line'
               stop 1
            endif
            cli_reffile = trim(cli_val)
            cli_seen_reffile = .true.
         case ('infiles')
            if (cli_seen_infiles) then
               write(*,*) 'ERROR: infiles given more than once on the command line'
               stop 1
            endif
            raw_cli_infiles = trim(cli_val)
            cli_seen_infiles = .true.
         case ('mem_frac_ram')
            if (cli_seen_mem_frac_ram) then
               write(*,*) 'ERROR: mem_frac_ram given more than once on the command line'
               stop 1
            endif
            read(cli_val, *, iostat=ios_mfr) cli_mem_frac_ram
            if (ios_mfr.ne.0) then
               write(*,*) 'ERROR: mem_frac_ram must be a number, got "',&
               &trim(cli_val), '"'
               stop 1
            endif
            cli_seen_mem_frac_ram = .true.
         case ('io_overlap')
            if (cli_seen_io_overlap) then
               write(*,*) 'ERROR: io_overlap given more than once on the command line'
               stop 1
            endif
            cli_io_overlap = flag_from_value_reproject(cli_val)
            cli_seen_io_overlap = .true.
         case ('log_level')
            if (cli_seen_log_level) then
               write(*,*) 'ERROR: log_level given more than once on the command line'
               stop 1
            endif
            cli_log_level = trim(cli_val)
            cli_seen_log_level = .true.
         case ('timing_enabled')
            if (cli_seen_timing_enabled) then
               write(*,*) 'ERROR: timing_enabled given more than once on the command line'
               stop 1
            endif
            cli_timing_enabled = flag_from_value_logging(cli_val)
            cli_seen_timing_enabled = .true.
         case ('log_output_file')
            if (cli_seen_log_output_file) then
               write(*,*) 'ERROR: log_output_file given more than once on the command line'
               stop 1
            endif
            cli_log_output_file = trim(cli_val)
            cli_seen_log_output_file = .true.
         case default
            write(*,*) 'ERROR: unrecognised key "', trim(cli_key), '" -- expected',&
            &' mode, reffile, infiles, mem_frac_ram, io_overlap, log_level,'//&
            &' timing_enabled, or log_output_file'
            stop 1
         end select
         iarg = iarg + 1
      endif
   enddo

   seen_mode = .false.
   seen_reffile = .false.
   seen_infiles = .false.
   if (have_cfgfile) then
      call read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs,&
      &mem_frac_ram, io_overlap, log_level, timing_enabled, log_output_file,&
      &status)
      if (status.ne.0) stop 1
      seen_mode = .true.
      seen_reffile = .true.
      seen_infiles = .true.
   endif

   if (cli_seen_mode) then
      mode = cli_mode
      seen_mode = .true.
   endif
   if (cli_seen_reffile) then
      reffile = cli_reffile
      seen_reffile = .true.
   endif
   if (cli_seen_infiles) then
      n_inputs = cfg_csv_count(raw_cli_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.max_inputs) then
         write(*,*) 'ERROR: infiles must list between 1 and ', max_inputs, ' files'
         stop 1
      endif
      do i = 1, n_inputs
         call cfg_csv_get_item(raw_cli_infiles, i, infiles(i))
      enddo
      seen_infiles = .true.
   endif
   if (cli_seen_mem_frac_ram) mem_frac_ram = cli_mem_frac_ram
   if (cli_seen_io_overlap) io_overlap = cli_io_overlap
   if (cli_seen_log_level) log_level = cli_log_level
   if (cli_seen_timing_enabled) timing_enabled = cli_timing_enabled
   if (cli_seen_log_output_file) log_output_file = cli_log_output_file

   if (.not. seen_mode .or. .not. seen_reffile .or. .not. seen_infiles) then
      call print_usage()
      stop 1
   endif

   if (mem_frac_ram.le.0.0 .or. mem_frac_ram.gt.0.95) then
      write(*,*) 'ERROR: mem_frac_ram must be > 0 and <= 0.95, got ', mem_frac_ram
      stop 1
   endif

   if (trim(mode).ne.'intersection' .and. trim(mode).ne.'union' .and.&
   &trim(mode).ne.'reference') then
      write(*,*) 'ERROR: mode must be intersection, union, or reference'
      stop 1
   endif

   call init_logging(log_level, timing_enabled, log_output_file, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: cannot open log_output_file: ', trim(log_output_file)
      stop 1
   endif
   call log_message('info', 'startup', 'reproject_cubes run started')

   status = 0
   call ast_begin(status)

   call load_wcs(reffile, wcs_ref, naxes_ref, status)
   call extract_sky_mapping(wcs_ref, skymap_ref, skyframe_ref, pixaxes_ref, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to load the reference file''s WCS'
      stop 1
   endif

   ! Output grid starts as the reference's own full extent; intersection
   ! shrinks it, union grows it, reference mode leaves it untouched (the
   ! loop below is skipped entirely for reference mode).
   lbnd_out(1) = 1.0d0
   lbnd_out(2) = 1.0d0
   ubnd_out(1) = real(naxes_ref(pixaxes_ref(1)), kind=8)
   ubnd_out(2) = real(naxes_ref(pixaxes_ref(2)), kind=8)
   write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Reference (', trim(reffile),&
   &') own extent: [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
   &lbnd_out(2), ',', ubnd_out(2), ']'

   if (trim(mode).ne.'reference') then
      do i = 1, n_inputs
         call load_wcs(infiles(i), wcs_in, naxes_in, status)
         call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to load input file: ', trim(infiles(i))
            stop 1
         endif

         call compose_pix2pix(skymap_in, skyframe_in, skymap_ref, skyframe_ref,&
         &map_in2ref, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to align input file to the reference: ',&
            &trim(infiles(i))
            stop 1
         endif

         call footprint_bounds(map_in2ref, naxes_in, pixaxes_in,&
         &this_lbnd, this_ubnd, status)
         write(*,'(A,A,A,F0.2,A,F0.2,A,F0.2,A,F0.2,A)') '  ', trim(infiles(i)),&
         &' footprint in reference space: [', this_lbnd(1), ',', this_ubnd(1),&
         &'] x [', this_lbnd(2), ',', this_ubnd(2), ']'

         ! Zero overlap with the running output grid is always a hard
         ! failure, regardless of mode -- a band sharing no sky at all
         ! with the rest is not legitimate partial coverage (that implies
         ! at least some shared sky), it is almost certainly the wrong
         ! file. Checked here, per-input, before folding it into the
         ! running bound, so the diagnostic names the specific offending
         ! file rather than a generic "empty result" at the end.
         if (this_ubnd(1).lt.lbnd_out(1) .or. this_lbnd(1).gt.ubnd_out(1) .or.&
         &this_ubnd(2).lt.lbnd_out(2) .or. this_lbnd(2).gt.ubnd_out(2)) then
            write(*,*) 'ERROR: zero sky overlap between the reference and: ',&
            &trim(infiles(i))
            write(*,*) 'Quitting now...'
            stop 1
         endif

         if (trim(mode).eq.'intersection') then
            lbnd_out = max(lbnd_out, this_lbnd)
            ubnd_out = min(ubnd_out, this_ubnd)
         else ! union
            lbnd_out = min(lbnd_out, this_lbnd)
            ubnd_out = max(ubnd_out, this_ubnd)
         endif

         call ast_annul(map_in2ref, status)
         call ast_annul(skymap_in, status)
         call ast_annul(skyframe_in, status)
         call ast_annul(wcs_in, status)
      enddo
   endif

   ! Round to the integer pixel grid: intersection wants the largest
   ! integer range fully CONTAINED within the real-valued bound (ceiling
   ! the lower edge, floor the upper edge); union wants the smallest
   ! integer range that fully CONTAINS it (floor/ceiling the other way).
   ! Reference mode's bound is already exactly integer (a file's own
   ! NAXIS), so this is a no-op for it either way.
   if (trim(mode).eq.'intersection') then
      lbnd_out = ceiling(lbnd_out)
      ubnd_out = floor(ubnd_out)
   else
      lbnd_out = floor(lbnd_out)
      ubnd_out = ceiling(ubnd_out)
   endif

   if (lbnd_out(1).gt.ubnd_out(1) .or. lbnd_out(2).gt.ubnd_out(2)) then
      write(*,*) 'ERROR: computed output grid is empty (', trim(mode), ' mode)'
      stop 1
   endif
   nx_out_common = nint(ubnd_out(1) - lbnd_out(1)) + 1
   ny_out_common = nint(ubnd_out(2) - lbnd_out(2)) + 1

   ! === Skip-if-already-matched pre-flight (planning-doc ticket) ===
   ! Same scheme as match_cubes.f90's own (see there for the full
   ! rationale): every safety check and skip decision for the WHOLE
   ! batch happens up front, before any file is processed, so a bad run
   ! (a stale output already on disk) fails fast. A pre-existing output
   ! path is always refused, never silently reused or overwritten,
   ! regardless of what this run's own skip decision would have been.
   do i = 1, n_inputs
      inquire(file=trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS', exist=out_exists)
      if (out_exists) then
         write(*,*) 'ERROR: output path already exists, refusing to proceed'//&
         &' (stale output from a previous run? remove it first): ',&
         &trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS'
         stop 1
      endif
   enddo

   fitsstat_skip = 0
   call safe_ftopen(ref_skip_unit, trim(reffile), 0, blocksize_skip, fitsstat_skip)
   if (fitsstat_skip.ne.0) then
      write(*,*) 'ERROR: cannot reopen reference file for the geometry check: ',&
      &trim(reffile)
      stop 1
   endif

   do i = 1, n_inputs
      call load_wcs(infiles(i), wcs_in, naxes_in, status)
      call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read input''s WCS for the geometry check: ',&
         &trim(infiles(i))
         stop 1
      endif
      fitsstat_skip = 0
      call safe_ftopen(cand_skip_unit, trim(infiles(i)), 0, blocksize_skip, fitsstat_skip)
      if (fitsstat_skip.ne.0) then
         write(*,*) 'ERROR: cannot reopen input for the geometry check: ',&
         &trim(infiles(i))
         stop 1
      endif
      needs_processing(i) = .not. sky_wcs_matches_target(cand_skip_unit,&
      &pixaxes_in(1), pixaxes_in(2), naxes_in(pixaxes_in(1)), naxes_in(pixaxes_in(2)),&
      &ref_skip_unit, pixaxes_ref(1), pixaxes_ref(2), lbnd_out(1)-1.0d0,&
      &lbnd_out(2)-1.0d0, nx_out_common, ny_out_common)
      call safe_ftclos(cand_skip_unit, fitsstat_skip)
      call ast_annul(skymap_in, status)
      call ast_annul(skyframe_in, status)
      call ast_annul(wcs_in, status)
      if (.not. needs_processing(i)) then
         write(*,'(A,A,A)') 'SKIP: ', trim(infiles(i)),&
         &' already matches the target grid -- no output written, use it directly'
      endif
   enddo
   fitsstat_skip = 0
   call safe_ftclos(ref_skip_unit, fitsstat_skip)

   ! --- Resample and write every plane of every input onto the final
   ! output grid, via astResampleR + FTPSSE ---
   ! Only naxes_in/pixaxes_in (plain per-file array-shape data, not an
   ! AST Object) are needed here -- write_reprojected_file builds its
   ! own private input->reference Mapping per OpenMP thread internally
   ! (see its own comment for why), so this loop no longer needs to
   ! derive map_in2ref itself the way it used to.
   do i = 1, n_inputs
      if (.not. needs_processing(i)) cycle
      call load_wcs(infiles(i), wcs_in, naxes_in, status)
      call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read input''s WCS for resampling: ',&
         &trim(infiles(i))
         stop 1
      endif
      ! Plain filename (NOT '!'-prefixed): write_reprojected_file's own
      ! FTINIT fails if the output already exists rather than silently
      ! clobbering it. Fixed bug, found while implementing match_cubes'
      ! own skip-if-already-matched feature (planning-doc ticket) and
      ! auditing every FTINIT call site for the same pattern: this call
      ! previously baked the '!'-prefix CLOBBER convention into the
      ! output filename itself, silently deleting and overwriting a
      ! pre-existing output with no warning.
      call write_reprojected_file(reffile, infiles(i),&
      &trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS', pixaxes_ref,&
      &naxes_in, pixaxes_in, lbnd_out, ubnd_out, mem_frac_ram, io_overlap,&
      &status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write reprojected output for: ',&
         &trim(infiles(i))
         stop 1
      endif
      call ast_annul(skymap_in, status)
      call ast_annul(skyframe_in, status)
      call ast_annul(wcs_in, status)
   enddo

   write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Final output grid (',&
   &trim(mode), ' mode): [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
   &lbnd_out(2), ',', ubnd_out(2), ']'

   call ast_annul(skymap_ref, status)
   call ast_annul(skyframe_ref, status)
   call ast_annul(wcs_ref, status)
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: footprint-mode output grid computed successfully.'
   call timer_report_summary()
   call log_message('info', 'finalize', 'reproject_cubes run completed')

contains

   subroutine write_reprojected_file(reffile, infile, outfile, pixaxes_ref,&
   &naxes_in, pixaxes_in, lbnd_out_d, ubnd_out_d, mem_frac_ram, io_overlap_l,&
   &status)
      !! Create outfile and write every reprojected plane of infile into
      !! it. Output axis layout: the 2 sky axes always occupy OUTPUT
      !! positions 1,2 (matching what rm_synthesis's own tile-read I/O
      !! already assumes -- RA-fastest-on-disk -- confirmed against
      !! rm_synthesis_mod.f90's own auto-tile-planner comment), in
      !! whatever order the REFERENCE file itself presents them (pixaxes_ref);
      !! this guarantees "sky is axes 1,2" but only guarantees "axis 1 is
      !! literally RA" if the reference is conventionally ordered (RA
      !! before Dec) -- a deliberate, documented scope limit, not a
      !! silent assumption: full semantic RA/Dec canonicalisation
      !! regardless of the reference's own convention is a follow-on
      !! refinement. The other (non-sky) axes keep their own values from
      !! INFILE unchanged (that band's own channel/Stokes definitions,
      !! untouched by reprojection), placed at output positions 3.. in
      !! their original relative order.
      !! Planes are read, resampled, and written in BLOCKS, not one at a
      !! time: block size comes from mem_frac_ram (get_mem_total_kb
      !! below mirrors rm_synthesis's own /proc/meminfo MemTotal read
      !! exactly, same concept as its plan_tile). Each block goes through
      !! three strictly separated phases -- read (one thread, via OMP
      !! `single`), resample (every thread, in parallel, one
      !! astResampleR call per plane in the block), write (one thread,
      !! via `single`) -- relying on the implicit barrier every `single`/
      !! `do` construct has by default. This replaced an earlier
      !! plane-at-a-time version that serialised every CFITSIO call
      !! behind an OMP critical section: with blocks, a whole block's I/O
      !! happens on a single thread BY CONSTRUCTION, so there is no
      !! concurrent CFITSIO access to guard against in the first place --
      !! nothing to lock, rather than a lock relied on to make concurrent
      !! access safe. Batching also amortises CFITSIO's per-call
      !! (FTOPEN/FTGSVE/FTPSSE) overhead across many planes instead of
      !! paying it every single one.
      use, intrinsic :: ieee_arithmetic
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
      character(len=*), intent(in) :: reffile, infile, outfile
      integer, intent(in) :: pixaxes_ref(2)
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      real, intent(in) :: mem_frac_ram
      logical, intent(in) :: io_overlap_l
      integer, intent(inout) :: status

      integer :: cur_slot
      logical :: write_dispatched_ok
      integer :: naxis, k, other_axes(max_axes), n_other
      integer :: other_idx(max_axes), remainder, radix
      integer :: n_planes, status_par, nthreads
      integer :: nx_out, ny_out, naxis_out, naxes_out(max_axes)
      integer :: nx_in, ny_in
      integer :: ref_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: beams_unit, beams_status, casambm_status, hdutype_dum
      logical :: casambm_val
      character(len=80) :: comment

      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes
      integer(kind=8) :: block_planes64
      integer :: block_planes, n_groups, igroup, axis1_extent
      integer :: chan_start, chan_len, local_iplane
      real, allocatable, target :: block_data_in(:,:,:), block_data_out(:,:,:,:)

      ! Per-OpenMP-thread private AST working set (see the parallel
      ! region below for why each thread builds its own).
      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref

      ! Per-plane resample working set (private, one astResampleR call
      ! per plane inside a block's parallel do).
      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      real :: badval
      double precision :: params_dummy(1)
      real(dp) :: t_stage
      ! Per-thread swim-lane instrumentation (planning-doc ticket) -- see
      ! convolve_cubes.f90's own write_convolved_file for the full
      ! rationale, and match_cubes.f90's own process_one_file_general
      ! for why iblock is PRIVATE here (persistent parallel region, every
      ! thread redundantly computes the same chan_start/chan_len
      ! progression).
      integer :: iblock, tid_local
      real(dp) :: t_thread_start, t_thread_elapsed
      character(len=160) :: thread_msg

      if (status.ne.0) return
      call log_message('info', 'reproject', 'starting: '//trim(infile))

      naxis = 0
      do k = 1, size(naxes_in)
         if (naxes_in(k).gt.0) naxis = k
      enddo
      n_other = 0
      do k = 1, naxis
         if (k.ne.pixaxes_in(1) .and. k.ne.pixaxes_in(2)) then
            n_other = n_other + 1
            other_axes(n_other) = k
         endif
      enddo

      nx_in = naxes_in(pixaxes_in(1))
      ny_in = naxes_in(pixaxes_in(2))
      nx_out = nint(ubnd_out_d(1) - lbnd_out_d(1)) + 1
      ny_out = nint(ubnd_out_d(2) - lbnd_out_d(2)) + 1
      naxis_out = 2 + n_other
      naxes_out(1) = nx_out
      naxes_out(2) = ny_out
      do k = 1, n_other
         naxes_out(2+k) = naxes_in(other_axes(k))
      enddo

      ! --- Create the output file and its primary header ---
      fitsstat = 0
      blocksize = 1
      call safe_ftinit(out_unit, trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to create output file: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         call free_fits_unit(out_unit)
         return
      endif
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis_out, naxes_out(1:naxis_out),&
      &0, 1, extend, fitsstat)

      ! Sky axes (output 1,2): WCS copied from the REFERENCE's own
      ! pixaxes_ref-numbered keywords, CRPIX shifted for the output
      ! grid's own origin (same formula rm_synthesis.f90 already uses for
      ! its own subimage CRPIX shift: (old_crpix - pixel_offset) + 1,
      ! stride 1 since this is a crop/grow, never a sub-sample).
      fitsstat = 0
      call safe_ftopen(ref_unit, trim(reffile), 0, blocksize, fitsstat)
      call copy_axis_keywords(ref_unit, pixaxes_ref(1), out_unit, 1,&
      &lbnd_out_d(1)-1.0d0, status)
      call copy_axis_keywords(ref_unit, pixaxes_ref(2), out_unit, 2,&
      &lbnd_out_d(2)-1.0d0, status)
      ! PCi_j/CDi_j (the 2x2-matrix alternative to CROTA, which
      ! copy_axis_keywords just above already handles per-axis): also
      ! from the reference, same reasoning as CROTA -- the output grid's
      ! orientation is the reference's own, unaffected by cropping/
      ! growing its extent.
      call copy_sky_rotation_matrix(ref_unit, pixaxes_ref(1),&
      &pixaxes_ref(2), out_unit, status)
      ! EQUINOX/RADESYS describe the coordinate SYSTEM, not a specific
      ! axis or the data -- the output is expressed in the REFERENCE's
      ! sky frame (compose_pix2pix converts every input's sky coordinates
      ! into it), so these follow reffile, not infile, same reasoning as
      ! the sky axes themselves just above.
      call copy_wcs_system_keywords(ref_unit, out_unit, status)
      call safe_ftclos(ref_unit, fitsstat)

      ! Other axes (output 3..): WCS copied from INFILE's own axis
      ! numbering unchanged (no CRPIX shift -- reprojection never touches
      ! these axes).
      fitsstat = 0
      call safe_ftopen(ref_unit, trim(infile), 0, blocksize, fitsstat)
      do k = 1, n_other
         call copy_axis_keywords(ref_unit, other_axes(k), out_unit,&
         &2+k, 0.0d0, status)
      enddo
      ! Everything else non-structural and non-axis-indexed (BUNIT,
      ! BMAJ/BMIN/BPA, OBJECT, TELESCOP, DATE-OBS, RESTFRQ, HISTORY,
      ! COMMENT, ...) describes the DATA, unchanged by reprojection --
      ! follows infile, matching the other axes just above. Previously
      ! dropped entirely (only per-axis WCS keywords were ever copied),
      ! silently losing units/beam/provenance from the output.
      call copy_generic_header(ref_unit, out_unit, status)
      call FTPHIS(out_unit, 'reproject_cubes: reprojected from '//&
      &trim(infile)//' onto the grid of '//trim(reffile), fitsstat)
      call safe_ftclos(ref_unit, fitsstat)

      ! CASAMBM/BEAMS: reprojection is a pure spatial resample of each
      ! plane -- it never touches the beam itself (BMAJ/BMIN/BPA are
      ! stored in sky degrees, not pixels), so a genuine per-channel
      ! BEAMS table on infile is still exactly correct for the
      ! reprojected output and should follow it unchanged. The scalar
      ! CASAMBM keyword itself already rode along verbatim as a raw
      ! header card inside copy_generic_header just above; that copy
      ! only ever touches the PRIMARY header, though, and can't reach
      ! the separate BEAMS extension HDU the keyword refers to -- this
      ! attaches that extension explicitly, matching the ftcopy pattern
      ! rm_synthesis.f90 already uses for its own CASAMBM/BEAMS
      ! passthrough.
      casambm_status = 0
      call ftgkyl(out_unit, 'CASAMBM', casambm_val, comment, casambm_status)
      if (casambm_status.eq.0 .and. casambm_val) then
         beams_status = 0
         call safe_ftopen(beams_unit, trim(infile), 0, blocksize, beams_status)
         call ftmnhd(beams_unit, -1, 'BEAMS', 0, beams_status)
         if (beams_status.eq.0) then
            status = 0
            call ftcopy(beams_unit, out_unit, 0, status)
            call ftmahd(out_unit, 1, hdutype_dum, status)
         else
            write(*,*) 'WARNING: CASAMBM=T but no BEAMS extension found in: ',&
            &trim(infile), ' -- output keeps CASAMBM=T with no BEAMS table.'
         endif
         call safe_ftclos(beams_unit, beams_status)
      endif

      ! --- Block size: mem_frac_ram fraction of total system RAM,
      ! divided by the bytes one plane's worth of input+output costs
      ! (same budgeting concept as rm_synthesis's plan_tile, applied to
      ! "planes in a block" instead of "output pixels in a tile"). Never
      ! larger than a single cycle of other_axes(1) (the fastest-varying
      ! non-sky axis) -- a block always covers one CONTIGUOUS range of
      ! that axis with every slower other axis held fixed (see the
      ! group/block loop below), so it can never usefully span past
      ! where that axis wraps.
      n_planes = 1
      do k = 1, n_other
         n_planes = n_planes * naxes_in(other_axes(k))
      enddo
      call get_mem_total_kb(mem_total_kb)
      ! block_data_out's own term is x2 (not x1) -- double-buffered for
      ! io_overlap (see below), budgeted whether or not it's actually on.
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &2_8*int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram,8) * real(mem_total_kb,8) *&
      &1024.0d0, 8)
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes64 = min(block_planes64, max_elements_per_block /&
      &max(1_8, max(int(nx_in,8)*int(ny_in,8), int(nx_out,8)*int(ny_out,8))))
      block_planes64 = max(1_8, block_planes64)
      block_planes = int(min(block_planes64, int(n_planes,8)))
      axis1_extent = 1
      if (n_other.ge.1) then
         axis1_extent = naxes_in(other_axes(1))
         block_planes = min(block_planes, axis1_extent)
      endif
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ',&
      &n_planes, ' plane(s), in blocks of up to ', block_planes, ' plane(s)'

      ! Parallelism only happens WITHIN a block (its plane range, via the
      ! `!$omp do` below) -- blocks themselves are strictly sequential,
      ! never overlapped/pipelined, so block_planes is a hard ceiling on
      ! how many threads can ever do useful work at once, not just a
      ! memory knob. A too-small mem_frac_ram therefore doesn't just add
      ! CFITSIO call overhead -- it can silently throw away most of the
      ! OpenMP speedup: measured 12.0s (block_planes=1, forced to 1
      ! thread) vs. 5.2s (one block covering the whole cube, 16 threads)
      ! on the same 1024x1024x300 synthetic cube, over 2x slower from
      ! losing parallelism, not from smaller reads/writes. Warn rather
      ! than silently eat that cost or override the user's own budget.
      if (block_planes.lt.omp_get_max_threads()) then
         write(*,'(A,I0,A,I0,A)') 'WARNING: mem_frac_ram limits blocks to ',&
         &block_planes, ' plane(s), below the ', omp_get_max_threads(),&
         &' threads available -- parallelism (not just I/O) is reduced;'//&
         &' raise mem_frac_ram for full speedup if memory allows.'
      endif

      allocate(block_data_in(nx_in, ny_in, block_planes))
      allocate(block_data_out(nx_out, ny_out, block_planes, 0:1))

      n_groups = 1
      do k = 2, n_other
         n_groups = n_groups * naxes_in(other_axes(k))
      enddo

      cur_slot = 0
      write_pending = .false.
      write_failed = .false.
      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      !$omp parallel num_threads(nthreads) default(none)&
      !$omp& shared(infile, reffile, naxes_in, pixaxes_in, other_axes,&
      !$omp& n_other, lbnd_out_d, ubnd_out_d, out_unit, status_par,&
      !$omp& nx_in, ny_in, nx_out, ny_out, block_planes, block_data_in,&
      !$omp& block_data_out, n_groups, axis1_extent, cur_slot,&
      !$omp& io_overlap_l, write_dispatched_ok, write_pending,&
      !$omp& write_thread_id, write_failed, write_job, t_stage)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, other_idx, remainder,&
      !$omp& radix, k, igroup, chan_start, chan_len, local_iplane, nbad,&
      !$omp& lbnd_in, ubnd_in, lbnd_o, ubnd_o, badval, params_dummy,&
      !$omp& iblock, tid_local, t_thread_start, t_thread_elapsed, thread_msg)

      ! Each thread builds its OWN private input->reference pixel
      ! Mapping from scratch (own ast_begin context, own load_wcs/
      ! extract_sky_mapping/compose_pix2pix calls) rather than sharing or
      ! cloning one Mapping built by the caller. SUN/211 Sec 4.12 ("AST
      ! Objects within Multi-threaded Applications") requires astLock/
      ! astUnlock to hand an AST Object from the thread that created it
      ! to another thread; this Fortran binding doesn't export
      ! ast_lock_/ast_unlock_ (checked libstarlink_ast.so.9's actual
      ! symbol table) -- only ast_copy_ is present, useless on its own
      ! without the paired lock handoff. Every thread's whole AST object
      ! graph is therefore self-created and never shared, squarely
      ! inside AST's documented thread-safe model without needing
      ! lock/unlock at all. CFITSIO unit numbers used for this (in
      ! load_wcs) come from fitsio_unit_mod's safe_ftopen (CFITSIO's own
      ! FTGIOU+FTOPEN, both inside one critical section), so each thread
      ! gets a genuinely distinct unit with no collision risk regardless of
      ! thread count (see load_wcs's own comment for the earlier version
      ! that hand-picked a unit offset and what it actually broke: a
      ! plain unit collision, not an AST/CFITSIO concurrency limit).
      iblock = 0
      t_status = 0
      call ast_begin(t_status)
      call load_wcs(reffile, t_wcs_ref, t_naxes_ref, t_status)
      call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      &t_pixaxes_ref, t_status)
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      if (t_status.ne.0) then
         !$omp atomic write
         status_par = -1
      endif

      ! Outer "group" loop: one iteration per combination of the SLOWER
      ! non-sky axes (other_axes(2:n_other)), decoded from igroup via the
      ! same mixed-radix scheme used elsewhere in this file. Every thread
      ! computes the same igroup/other_idx/chan_start/chan_len values
      ! redundantly (cheap, deterministic from shared inputs) so they all
      ! hit the same single/do/single sequence together, as OpenMP
      ! worksharing constructs require.
      do igroup = 1, n_groups
         remainder = igroup - 1
         do k = 2, n_other
            radix = naxes_in(other_axes(k))
            other_idx(k) = mod(remainder, radix) + 1
            remainder = remainder / radix
         enddo

         chan_start = 1
         do while (chan_start.le.axis1_extent)
            chan_len = min(block_planes, axis1_extent - chan_start + 1)

            !$omp single
            call timer_start(t_stage)
            call read_one_block(infile, naxes_in, pixaxes_in, other_axes,&
            &other_idx, n_other, chan_start, chan_len, nx_in, ny_in,&
            &block_data_in(:,:,1:chan_len), status_par)
            call timer_stop('block_read', t_stage)
            call timer_start(t_stage)
            !$omp end single

            iblock = iblock + 1
            tid_local = omp_get_thread_num()
            t_thread_start = omp_get_wtime()
            write(thread_msg,'(A,I0,A,I0,A,I0)') 'thread_timing stage=resample event=start tid=',&
            &tid_local,' block=',iblock,' unit_count=',chan_len
            call log_message('debug','tile_thread',trim(thread_msg))
            !$omp do schedule(static)
            do local_iplane = 1, chan_len
               if (t_status.eq.0 .and. status_par.eq.0) then
                  lbnd_in(1) = 1
                  lbnd_in(2) = 1
                  ubnd_in(1) = nx_in
                  ubnd_in(2) = ny_in
                  lbnd_o(1) = nint(lbnd_out_d(1))
                  lbnd_o(2) = nint(lbnd_out_d(2))
                  ubnd_o(1) = nint(ubnd_out_d(1))
                  ubnd_o(2) = nint(ubnd_out_d(2))
                  badval = ieee_value(badval, ieee_quiet_nan)
                  params_dummy(1) = 0.0d0
                  ! Same 20-argument astResampleR signature/pitfall as
                  ! before (see git history for the "missing params"
                  ! segfault this avoids); in_var/out_var unused, block
                  ! arrays reused as harmless placeholders.
                  nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                  &block_data_in(:,:,local_iplane),&
                  &block_data_in(:,:,local_iplane),&
                  &ast__linear, ast_null, params_dummy, 0, 0.0d0, 100, badval,&
                  &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o,&
                  &block_data_out(:,:,local_iplane,cur_slot),&
                  &block_data_out(:,:,local_iplane,cur_slot), t_status)
                  if (t_status.ne.0) then
                     !$omp atomic write
                     status_par = -1
                  endif
               endif
            enddo
            !$omp end do
            t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
            write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)') 'thread_timing stage=resample event=done tid=',&
            &tid_local,' block=',iblock,' unit_count=',chan_len,' dur_ms=',t_thread_elapsed
            call log_message('debug','tile_thread',trim(thread_msg))

            !$omp single
            call timer_stop('block_resample', t_stage)

            call timer_start(t_stage)
            if (write_pending) then
               call block_write_join(write_thread_id)
               write_pending = .false.
               if (write_failed) status_par = -1
            endif
            call timer_stop('block_write_join', t_stage)

            call timer_start(t_stage)
            if (status_par.eq.0) then
               write_job%out_unit = out_unit
               write_job%naxes_in(1:max_axes) = naxes_in(1:max_axes)
               write_job%other_axes(1:max_axes) = other_axes(1:max_axes)
               write_job%other_idx(1:max_axes) = other_idx(1:max_axes)
               write_job%n_other = n_other
               write_job%chan_start = chan_start
               write_job%chan_len = chan_len
               write_job%nx = nx_out
               write_job%ny = ny_out
               write_job%data => block_data_out(:,:,1:chan_len,cur_slot)
               if (io_overlap_l) then
                  call block_write_dispatch_async(write_job, write_thread_id,&
                  &write_dispatched_ok)
                  write_pending = write_dispatched_ok
               else
                  call do_block_write(write_job)
                  if (write_failed) status_par = -1
               endif
            endif
            call timer_stop('block_write', t_stage)
            cur_slot = 1 - cur_slot
            !$omp end single

            chan_start = chan_start + chan_len
         enddo
      enddo

      call ast_end(t_status)
      !$omp end parallel

      call timer_start(t_stage)
      if (write_pending) then
         call block_write_join(write_thread_id)
         write_pending = .false.
         if (write_failed) status_par = -1
      endif
      call timer_stop('block_write_join', t_stage)

      deallocate(block_data_in)
      deallocate(block_data_out)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to resample/write one or more planes for: ',&
         &trim(infile)
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif

      call safe_ftclos(out_unit, fitsstat)
      call log_message('info', 'reproject', 'finished: '//trim(infile))
   end subroutine write_reprojected_file

   subroutine copy_axis_keywords(src_unit, src_axis, dst_unit, dst_axis,&
   &crpix_shift, status)
      !! Copy CTYPE/CRVAL/CRPIX/CDELT/CUNIT/CROTA for src_axis (in
      !! src_unit's own header) to dst_axis (in dst_unit's header,
      !! already created via FTPHPR). CRPIX is additionally shifted by
      !! -crpix_shift (0 for a straight passthrough; the output grid's
      !! own pixel-1 offset, reference-pixel-numbered, for a cropped/
      !! grown sky axis) -- matches rm_synthesis.f90's own existing
      !! subimage CRPIX-shift formula exactly, generalised to any axis
      !! number via a constructed keyword string ("CRVAL"//axis, etc.)
      !! rather than a literal "1"/"2" suffix. CROTA needs no shift --
      !! cropping/growing the same grid changes CRPIX (a position within
      !! it), never its rotation. (PCi_j/CDi_j, the 2x2-matrix
      !! alternative to CROTA, are handled separately by
      !! copy_sky_rotation_matrix below -- they're indexed by an axis
      !! *pair*, not a single axis, so they don't fit this per-axis
      !! copy loop.)
      integer, intent(in) :: src_unit, src_axis, dst_unit, dst_axis
      double precision, intent(in) :: crpix_shift
      integer, intent(inout) :: status

      integer :: fitsstat
      character(len=8) :: axstr
      character(len=68) :: comment
      character(len=68) :: sval
      double precision :: dval

      if (status.ne.0) return

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYS(src_unit, 'CTYPE'//trim(axstr), sval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYS(dst_unit, 'CTYPE'//trim(axstr), trim(sval), ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CRVAL'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CRVAL'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CRPIX'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         dval = dval - crpix_shift
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CRPIX'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CDELT'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CDELT'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYS(src_unit, 'CUNIT'//trim(axstr), sval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYS(dst_unit, 'CUNIT'//trim(axstr), trim(sval), ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CROTA'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CROTA'//trim(axstr), dval, 13, ' ', fitsstat)
      endif
   end subroutine copy_axis_keywords

   subroutine copy_sky_rotation_matrix(src_unit, src_axis1, src_axis2,&
   &dst_unit, status)
      !! PCi_j/CDi_j: the 2x2-matrix alternative to CROTA (a file uses
      !! one convention or the other, essentially never both -- copying
      !! whichever is actually present, doing nothing if neither is, is
      !! the correct behaviour either way). Indexed by an axis PAIR
      !! (e.g. PC2_4), so this reads the 4 entries for the sky axis pair
      !! specifically (src_axis1,src_axis2 -- in the REFERENCE file's own
      !! numbering) and writes them at the output's fixed sky positions
      !! (1,2), the same axis-renumbering copy_axis_keywords already does
      !! for CTYPE/CRVAL/etc, just for a 2-axis-indexed keyword instead
      !! of a 1-axis-indexed one. This is intentionally NOT attempted for
      !! the generic (infile, non-sky-axis) header copy elsewhere
      !! (copy_generic_header skips PC/CD entirely there) -- a rotation
      !! entry for a channel or Stokes axis would be highly unusual, and
      !! blindly relocating it under the input's own arbitrary axis
      !! numbers, unlike this deliberate sky-axis-pair copy, risks
      !! landing it on the wrong axis.
      integer, intent(in) :: src_unit, src_axis1, src_axis2, dst_unit
      integer, intent(inout) :: status

      if (status.ne.0) return

      call copy_one_matrix_entry(src_unit, 'PC', src_axis1, src_axis1,&
      &dst_unit, 1, 1)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis1, src_axis2,&
      &dst_unit, 1, 2)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis2, src_axis1,&
      &dst_unit, 2, 1)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis2, src_axis2,&
      &dst_unit, 2, 2)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis1, src_axis1,&
      &dst_unit, 1, 1)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis1, src_axis2,&
      &dst_unit, 1, 2)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis2, src_axis1,&
      &dst_unit, 2, 1)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis2, src_axis2,&
      &dst_unit, 2, 2)
   end subroutine copy_sky_rotation_matrix

   subroutine copy_one_matrix_entry(su, prefix, sa, sb, du, da, db)
      !! Single-entry helper for copy_sky_rotation_matrix (a sibling
      !! internal subroutine, not nested inside it -- Fortran doesn't
      !! allow a second level of CONTAINS). Copies keyword
      !! "<prefix><sa>_<sb>" from unit su to "<prefix><da>_<db>" on unit
      !! du, only if present on su (absent is not an error -- most files
      !! use CROTA or nothing, not PC/CD, so 0-of-8 entries found here is
      !! the common case).
      integer, intent(in) :: su, sa, sb, du, da, db
      character(len=*), intent(in) :: prefix

      character(len=16) :: srckey, dstkey
      character(len=68) :: comment
      double precision :: dval
      integer :: fitsstat
      character(len=4) :: si, sj

      write(si,'(I0)') sa
      write(sj,'(I0)') sb
      srckey = trim(prefix)//trim(si)//'_'//trim(sj)
      fitsstat = 0
      call FTGKYD(su, trim(srckey), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(si,'(I0)') da
         write(sj,'(I0)') db
         dstkey = trim(prefix)//trim(si)//'_'//trim(sj)
         call FTPKYD(du, trim(dstkey), dval, 13, ' ', fitsstat)
      endif
   end subroutine copy_one_matrix_entry

   logical function sky_wcs_matches_target(cand_unit, cand_axis1, cand_axis2,&
   &nx_cand, ny_cand, ref_unit, ref_axis1, ref_axis2, crpix_shift1,&
   &crpix_shift2, nx_out, ny_out) result(matches)
      !! Ticket: memory-in-plan-file "match_cubes: skip already-matched
      !! files". Tight, "already-processed-identically" comparison of a
      !! candidate file's own sky-axis WCS (cand_unit, its own axis
      !! numbering) against what THIS run's own output grid will
      !! actually be: the reference's own CTYPE/CRVAL/CDELT/rotation
      !! UNCHANGED, CRPIX shifted by crpix_shift1/2 -- the exact same
      !! values copy_axis_keywords/copy_sky_rotation_matrix above would
      !! themselves write, so this function is deliberately checking "is
      !! the candidate already what write_reprojected_file would
      !! produce" rather than any more general astrometric equivalence.
      !! NOT a general "are these two grids close enough" check -- see
      !! this ticket's own tolerance rationale (planning doc): a false
      !! match here would silently misalign downstream RM synthesis, so
      !! every comparison is a tight absolute tolerance, and a file that
      !! is very close but outside it is correctly treated as NOT
      !! matching (processed as normal, only unnecessary compute at
      !! stake, never correctness).
      integer, intent(in) :: cand_unit, cand_axis1, cand_axis2, nx_cand, ny_cand
      integer, intent(in) :: ref_unit, ref_axis1, ref_axis2
      double precision, intent(in) :: crpix_shift1, crpix_shift2
      integer, intent(in) :: nx_out, ny_out

      double precision, parameter :: tol_val = 1.0d-9
      double precision, parameter :: tol_rot = 1.0d-9
      character(len=68) :: ctype_c1, ctype_c2, ctype_r1, ctype_r2
      double precision :: crval_c1, crval_c2, crval_r1, crval_r2
      double precision :: crpix_c1, crpix_c2, crpix_r1, crpix_r2
      double precision :: cdelt_c1, cdelt_c2, cdelt_r1, cdelt_r2
      double precision :: pc_c(2,2), pc_r(2,2), crota_c1, crota_c2, crota_r1, crota_r2
      logical :: have_pc_c, have_pc_r

      matches = .false.

      ! Pixel extent must match exactly (integers -- no tolerance question).
      if (nx_cand.ne.nx_out .or. ny_cand.ne.ny_out) return

      call get_axis_sval(cand_unit, 'CTYPE', cand_axis1, ctype_c1)
      call get_axis_sval(cand_unit, 'CTYPE', cand_axis2, ctype_c2)
      call get_axis_sval(ref_unit, 'CTYPE', ref_axis1, ctype_r1)
      call get_axis_sval(ref_unit, 'CTYPE', ref_axis2, ctype_r2)
      if (trim(ctype_c1).ne.trim(ctype_r1) .or. trim(ctype_c2).ne.trim(ctype_r2)) return

      call get_axis_dval(cand_unit, 'CRVAL', cand_axis1, 0.0d0, crval_c1)
      call get_axis_dval(cand_unit, 'CRVAL', cand_axis2, 0.0d0, crval_c2)
      call get_axis_dval(ref_unit, 'CRVAL', ref_axis1, 0.0d0, crval_r1)
      call get_axis_dval(ref_unit, 'CRVAL', ref_axis2, 0.0d0, crval_r2)
      if (abs(crval_c1-crval_r1).gt.tol_val .or. abs(crval_c2-crval_r2).gt.tol_val) return

      call get_axis_dval(cand_unit, 'CDELT', cand_axis1, 1.0d0, cdelt_c1)
      call get_axis_dval(cand_unit, 'CDELT', cand_axis2, 1.0d0, cdelt_c2)
      call get_axis_dval(ref_unit, 'CDELT', ref_axis1, 1.0d0, cdelt_r1)
      call get_axis_dval(ref_unit, 'CDELT', ref_axis2, 1.0d0, cdelt_r2)
      if (abs(cdelt_c1-cdelt_r1).gt.tol_val .or. abs(cdelt_c2-cdelt_r2).gt.tol_val) return

      ! CRPIX: the candidate's own value must equal the reference's own
      ! value MINUS the same shift write_reprojected_file itself applies
      ! (copy_axis_keywords' own formula) -- this is the one place the
      ! candidate and reference are expected to differ even on a genuine
      ! match, by construction.
      call get_axis_dval(cand_unit, 'CRPIX', cand_axis1, 1.0d0, crpix_c1)
      call get_axis_dval(cand_unit, 'CRPIX', cand_axis2, 1.0d0, crpix_c2)
      call get_axis_dval(ref_unit, 'CRPIX', ref_axis1, 1.0d0, crpix_r1)
      call get_axis_dval(ref_unit, 'CRPIX', ref_axis2, 1.0d0, crpix_r2)
      if (abs(crpix_c1-(crpix_r1-crpix_shift1)).gt.tol_val) return
      if (abs(crpix_c2-(crpix_r2-crpix_shift2)).gt.tol_val) return

      ! Rotation: PCi_j/CDi_j takes precedence over CROTA if EITHER file
      ! has any entry present (same "essentially never both" convention
      ! copy_sky_rotation_matrix's own comment documents) -- absent PC
      ! entries default to the FITS standard identity matrix, absent
      ! CROTA defaults to 0.
      call get_matrix_2x2(cand_unit, cand_axis1, cand_axis2, pc_c, have_pc_c)
      call get_matrix_2x2(ref_unit, ref_axis1, ref_axis2, pc_r, have_pc_r)
      if (have_pc_c .or. have_pc_r) then
         if (any(abs(pc_c-pc_r).gt.tol_rot)) return
      else
         call get_axis_dval(cand_unit, 'CROTA', cand_axis1, 0.0d0, crota_c1)
         call get_axis_dval(cand_unit, 'CROTA', cand_axis2, 0.0d0, crota_c2)
         call get_axis_dval(ref_unit, 'CROTA', ref_axis1, 0.0d0, crota_r1)
         call get_axis_dval(ref_unit, 'CROTA', ref_axis2, 0.0d0, crota_r2)
         if (abs(crota_c1-crota_r1).gt.tol_rot .or. abs(crota_c2-crota_r2).gt.tol_rot) return
      endif

      matches = .true.
   end function sky_wcs_matches_target

   subroutine get_axis_sval(unit, prefix, axis, val)
      !! Read string keyword "<prefix><axis>" (e.g. CTYPE2); empty string
      !! if absent -- used only for equality comparison, so an absent
      !! keyword on both sides of a comparison still compares equal.
      integer, intent(in) :: unit, axis
      character(len=*), intent(in) :: prefix
      character(len=*), intent(out) :: val
      character(len=8) :: axstr
      character(len=68) :: comment
      integer :: fitsstat
      val = ' '
      write(axstr,'(I0)') axis
      fitsstat = 0
      call FTGKYS(unit, trim(prefix)//trim(axstr), val, comment, fitsstat)
      if (fitsstat.ne.0) val = ' '
   end subroutine get_axis_sval

   subroutine get_axis_dval(unit, prefix, axis, default_val, val)
      !! Read double-precision keyword "<prefix><axis>"; default_val if
      !! absent (matching the FITS standard's own default for that
      !! keyword -- 0 for CRVAL/CROTA, 1 for CDELT/CRPIX -- passed in by
      !! the caller rather than hard-coded here).
      integer, intent(in) :: unit, axis
      character(len=*), intent(in) :: prefix
      double precision, intent(in) :: default_val
      double precision, intent(out) :: val
      character(len=8) :: axstr
      character(len=68) :: comment
      integer :: fitsstat
      write(axstr,'(I0)') axis
      fitsstat = 0
      call FTGKYD(unit, trim(prefix)//trim(axstr), val, comment, fitsstat)
      if (fitsstat.ne.0) val = default_val
   end subroutine get_axis_dval

   subroutine get_matrix_2x2(unit, axis1, axis2, m, have_any)
      !! Read the 2x2 PCi_j (falling back to CDi_j if PC is entirely
      !! absent -- same "one or the other, not both" convention as
      !! copy_sky_rotation_matrix) rotation/scale matrix for the axis
      !! pair (axis1,axis2). have_any=.false. (m set to the identity)
      !! when NEITHER convention has any entry present -- the FITS
      !! standard default.
      integer, intent(in) :: unit, axis1, axis2
      double precision, intent(out) :: m(2,2)
      logical, intent(out) :: have_any
      logical :: any_pc, any_cd
      double precision :: mcd(2,2)

      m = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2,2/))
      call get_matrix_entry(unit, 'PC', axis1, axis1, m(1,1), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis1, axis2, m(1,2), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis2, axis1, m(2,1), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis2, axis2, m(2,2), any_pc)

      mcd = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2,2/))
      any_cd = .false.
      call get_matrix_entry(unit, 'CD', axis1, axis1, mcd(1,1), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis1, axis2, mcd(1,2), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis2, axis1, mcd(2,1), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis2, axis2, mcd(2,2), any_cd)
      if (any_cd) m = mcd

      have_any = any_pc .or. any_cd
   end subroutine get_matrix_2x2

   subroutine get_matrix_entry(unit, prefix, a, b, val, found_any)
      integer, intent(in) :: unit, a, b
      character(len=*), intent(in) :: prefix
      double precision, intent(out) :: val
      logical, intent(out) :: found_any
      found_any = .false.
      call get_matrix_entry_track(unit, prefix, a, b, val, found_any)
   end subroutine get_matrix_entry

   subroutine get_matrix_entry_track(unit, prefix, a, b, val, found_any)
      !! Read "<prefix><a>_<b>" if present, updating val and OR-ing into
      !! found_any; leaves both untouched if absent (val keeps whatever
      !! default the caller pre-seeded it with, e.g. the identity).
      integer, intent(in) :: unit, a, b
      character(len=*), intent(in) :: prefix
      double precision, intent(inout) :: val
      logical, intent(inout) :: found_any
      character(len=16) :: key
      character(len=68) :: comment
      character(len=4) :: sa, sb
      integer :: fitsstat
      double precision :: dval
      write(sa,'(I0)') a
      write(sb,'(I0)') b
      key = trim(prefix)//trim(sa)//'_'//trim(sb)
      fitsstat = 0
      call FTGKYD(unit, trim(key), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         val = dval
         found_any = .true.
      endif
   end subroutine get_matrix_entry_track

   subroutine copy_wcs_system_keywords(src_unit, dst_unit, status)
      !! EQUINOX/RADESYS describe the sky coordinate SYSTEM (not a
      !! specific axis), copied verbatim if present -- absent from
      !! src_unit is not an error (many FITS files simply omit them,
      !! relying on standard defaults; nothing to propagate then).
      integer, intent(in) :: src_unit, dst_unit
      integer, intent(inout) :: status

      integer :: fitsstat
      character(len=68) :: comment, sval
      double precision :: dval

      if (status.ne.0) return

      fitsstat = 0
      call FTGKYD(src_unit, 'EQUINOX', dval, comment, fitsstat)
      if (fitsstat.eq.0) call FTPKYD(dst_unit, 'EQUINOX', dval, 13, ' ', fitsstat)

      fitsstat = 0
      call FTGKYS(src_unit, 'RADESYS', sval, comment, fitsstat)
      if (fitsstat.eq.0) call FTPKYS(dst_unit, 'RADESYS', trim(sval), ' ', fitsstat)
   end subroutine copy_wcs_system_keywords

   subroutine copy_generic_header(src_unit, dst_unit, status)
      !! Copy every header card from src_unit (always called with this
      !! opened on INFILE, for the "other" -- non-sky -- axes) to
      !! dst_unit EXCEPT: structural keywords FTPHPR already wrote
      !! (SIMPLE/BITPIX/NAXIS/EXTEND/PCOUNT/GCOUNT/END), and any
      !! axis-indexed keyword (NAXISn/CTYPEn/CRVALn/CRPIXn/CDELTn/
      !! CUNITn/CROTAn/PCn_m/CDn_m, for ANY n/m). CTYPE/CRVAL/CRPIX/
      !! CDELT/CUNIT/CROTA are handled per-axis elsewhere
      !! (copy_axis_keywords), PCi_j/CDi_j elsewhere too
      !! (copy_sky_rotation_matrix), and EQUINOX/RADESYS elsewhere too
      !! (copy_wcs_system_keywords) -- but ALL THREE of those, unlike
      !! this generic copy, are only ever called against the REFERENCE
      !! file for the 2 sky axes specifically, with correct axis
      !! renumbering built in. A rotation/WCS-system keyword found here,
      !! on one of INFILE's own non-sky axes, is skipped rather than
      !! blindly relocated under this file's own arbitrary axis
      !! numbering (which the output does not share) -- landing it on
      !! the wrong axis would be worse than leaving it absent, and a
      !! genuine rotation entry on a channel/Stokes axis would be
      !! extremely unusual in practice anyway. Everything else -- BUNIT,
      !! BMAJ/BMIN/BPA, OBJECT, TELESCOP, INSTRUME, DATE-OBS, RESTFRQ,
      !! HISTORY, COMMENT, and anything this project has never heard of
      !! -- is copied verbatim via a raw 80-column header record
      !! (FTGREC/FTPREC), not decoded/re-encoded through a type-specific
      !! FTGKYx/FTPKYx pair, so formatting/precision/comments survive
      !! exactly. This is what actually fixes losing
      !! units/beam/provenance from the output (previously ONLY the
      !! axis-indexed WCS keywords were ever copied at all).
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
         if (skip_generic_header_key(key)) cycle
         fitsstat = 0
         call FTPREC(dst_unit, card, fitsstat)
      enddo
   end subroutine copy_generic_header

   logical function skip_generic_header_key(key)
      !! True for any keyword copy_generic_header must NOT copy verbatim
      !! -- see that subroutine's own comment for why each category is
      !! excluded.
      character(len=8), intent(in) :: key

      skip_generic_header_key = .true.
      select case (trim(key))
      case ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT',&
      &'END', 'EQUINOX', 'RADESYS', 'WCSAXES', 'LONPOLE', 'LATPOLE')
         return
      end select
      if (is_indexed_keyword(key, 'NAXIS')) return
      if (is_indexed_keyword(key, 'CTYPE')) return
      if (is_indexed_keyword(key, 'CRVAL')) return
      if (is_indexed_keyword(key, 'CRPIX')) return
      if (is_indexed_keyword(key, 'CDELT')) return
      if (is_indexed_keyword(key, 'CUNIT')) return
      if (is_indexed_keyword(key, 'CROTA')) return
      if (is_indexed_keyword(key, 'PC')) return
      if (is_indexed_keyword(key, 'CD')) return
      skip_generic_header_key = .false.
   end function skip_generic_header_key

   logical function is_indexed_keyword(key, prefix)
      !! True if key is exactly prefix followed by one or two groups of
      !! digits separated by at most one underscore (e.g. prefix='CD':
      !! matches "CD1_2", "CD10_2"; prefix='NAXIS': matches "NAXIS1",
      !! "NAXIS12"; does not match "CDELT1" -- the "ELT1" tail after
      !! stripping "CD" is not all-digits-and-at-most-one-underscore).
      character(len=8), intent(in) :: key
      character(len=*), intent(in) :: prefix

      character(len=8) :: rest
      integer :: plen, klen, i
      logical :: seen_underscore

      is_indexed_keyword = .false.
      plen = len_trim(prefix)
      klen = len_trim(key)
      if (klen.le.plen) return
      if (key(1:plen).ne.prefix(1:plen)) return

      rest = key(plen+1:klen)
      seen_underscore = .false.
      do i = 1, len_trim(rest)
         if (rest(i:i).eq.'_') then
            if (seen_underscore) return
            seen_underscore = .true.
         else if (rest(i:i).lt.'0' .or. rest(i:i).gt.'9') then
            return
         endif
      enddo
      if (len_trim(rest).eq.0) return
      if (rest(1:1).eq.'_' .or. rest(len_trim(rest):len_trim(rest)).eq.'_') return
      is_indexed_keyword = .true.
   end function is_indexed_keyword

   subroutine read_one_block(filename, naxes_in, pixaxes_in, other_axes,&
   &other_idx, n_other, chan_start, chan_len, nx_in, ny_in, block_data_in,&
   &status)
      !! Read chan_len consecutive planes of filename in one CFITSIO
      !! call: sky axes full extent, other_axes(1) spanning
      !! [chan_start, chan_start+chan_len-1], every slower other axis
      !! (other_axes(2:n_other)) fixed at other_idx(2:n_other) -- the
      !! write-side mirror of write_one_block below, same fpixels/lpixels
      !! construction as the old per-plane resample_one_plane just with a
      !! range instead of a single index on other_axes(1). Always called
      !! from a single OpenMP thread (write_reprojected_file's `!$omp
      !! single` region), so a fixed CFITSIO unit number is safe here --
      !! no other thread can be touching CFITSIO at the same moment, by
      !! construction, not by locking convention.
      !!
      !! FTGSVE fills its output array in ascending-axis-number order
      !! among the non-degenerate (extent>1) axes of THIS read: the 2 sky
      !! axes (always non-degenerate) and other_axes(1) (non-degenerate
      !! whenever chan_len>1). The single-plane version this replaced
      !! never had to care which of the 3 was numerically smallest --
      !! every non-sky axis was degenerate then (extent exactly 1), and a
      !! degenerate axis contributes no stride regardless of where it
      !! sits in the axis order. That stopped being true here: whenever
      !! other_axes(1) itself is numerically BEFORE one or both sky axes
      !! (e.g. TEST_NONADJACENT.Q.FITSCUBE: FREQ on axis 1, RA/DEC on
      !! axes 2 and 4 -- FREQ reads fastest, not slowest), the natural
      !! read order stops matching the caller's fixed
      !! (nx_in,ny_in,chan_len) block_data_in layout, and reading
      !! straight into it silently scrambles the data (caught by
      !! comparing block output against the pre-batching serial output
      !! byte-for-byte -- TEST_NONADJACENT differed in 102397 of 102400
      !! elements). The output side never has this problem -- its axis
      !! layout is always canonicalised to sky-first (see
      !! write_one_block) regardless of the input's own numbering.
      !!
      !! Fixed generally: read into a buffer shaped in the ACTUAL natural
      !! order (via pairwise-comparison ranking of the 3 axis numbers),
      !! then copy into the caller's fixed layout explicitly. Common case
      !! (other_axes(1) numerically after both sky axes, i.e. a
      !! conventionally-ordered cube with channels/Stokes after RA/Dec)
      !! reads directly into block_data_in with no extra buffer or copy
      !! -- the permute path only costs anything on the non-conventional
      !! axis orderings that actually need it.
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      integer, intent(in) :: other_axes(:), other_idx(:), n_other
      integer, intent(in) :: chan_start, chan_len, nx_in, ny_in
      real, intent(out) :: block_data_in(:,:,:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, group, naxis, k
      integer :: fpixels(max_axes), lpixels(max_axes), incs(max_axes)
      logical :: anyflg
      real :: badval
      integer :: ax_sky1, ax_sky2, ax_block
      integer :: rank_sky1, rank_sky2, rank_block
      integer :: dims(3), idxvec(3), i, j, c
      logical :: natural_order
      real, allocatable :: natural_buf(:,:,:)

      if (status.ne.0) return

      naxis = 0
      do k = 1, max_axes
         if (naxes_in(k).gt.0) naxis = k
      enddo
      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(pixaxes_in(1)) = nx_in
      lpixels(pixaxes_in(2)) = ny_in
      if (n_other.ge.1) then
         fpixels(other_axes(1)) = chan_start
         lpixels(other_axes(1)) = chan_start + chan_len - 1
      endif
      do k = 2, n_other
         fpixels(other_axes(k)) = other_idx(k)
         lpixels(other_axes(k)) = other_idx(k)
      enddo

      ! Rank the 3 axes that matter (sky1, sky2, block) by ascending raw
      ! axis number -- pairwise comparison, works for any 3 distinct
      ! integers. natural_order (ranks already 1,2,3 in that order) is
      ! the fast path; anything else needs the permute buffer.
      ax_sky1 = pixaxes_in(1)
      ax_sky2 = pixaxes_in(2)
      ax_block = 0
      if (n_other.ge.1) ax_block = other_axes(1)
      rank_sky1 = 1
      if (ax_sky2.lt.ax_sky1) rank_sky1 = rank_sky1 + 1
      if (n_other.ge.1 .and. ax_block.lt.ax_sky1) rank_sky1 = rank_sky1 + 1
      rank_sky2 = 1
      if (ax_sky1.lt.ax_sky2) rank_sky2 = rank_sky2 + 1
      if (n_other.ge.1 .and. ax_block.lt.ax_sky2) rank_sky2 = rank_sky2 + 1
      rank_block = 3
      if (n_other.ge.1) then
         rank_block = 1
         if (ax_sky1.lt.ax_block) rank_block = rank_block + 1
         if (ax_sky2.lt.ax_block) rank_block = rank_block + 1
      endif
      natural_order = (rank_sky1.eq.1 .and. rank_sky2.eq.2 .and. rank_block.eq.3)

      fitsstat = 0
      blocksize = 1
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (natural_order) then
         call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
         &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
         &badval, block_data_in, anyflg, fitsstat)
      else
         dims(rank_sky1) = nx_in
         dims(rank_sky2) = ny_in
         dims(rank_block) = chan_len
         allocate(natural_buf(dims(1), dims(2), dims(3)))
         call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
         &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
         &badval, natural_buf, anyflg, fitsstat)
         if (fitsstat.eq.0) then
            do c = 1, chan_len
               do j = 1, ny_in
                  do i = 1, nx_in
                     idxvec(rank_sky1) = i
                     idxvec(rank_sky2) = j
                     idxvec(rank_block) = c
                     block_data_in(i,j,c) = natural_buf(idxvec(1), idxvec(2), idxvec(3))
                  enddo
               enddo
            enddo
         endif
         deallocate(natural_buf)
      endif
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start,&
         &'-', chan_start+chan_len-1, ') from ', trim(filename)
         call printerror(fitsstat)
         status = -1
      endif
   end subroutine read_one_block

   subroutine write_one_block(out_unit, naxes_in, other_axes, other_idx,&
   &n_other, chan_start, chan_len, nx_out, ny_out, block_data_out, status)
      !! Write chan_len consecutive resampled planes in one CFITSIO call.
      !! Output axis layout (see write_reprojected_file's own comment):
      !! sky always at output positions 1,2, full extent; other axes at
      !! output positions 3.. (2+k for other_axes(k)) -- other_axes(1)'s
      !! output slot spans the block's own chan_start:chan_start+
      !! chan_len-1 range, every slower other axis fixed at
      !! other_idx(2:n_other), same value just relocated to the output
      !! file's own axis numbering. Always called from a single OpenMP
      !! thread by construction -- see read_one_block's own comment on
      !! why that means no locking is needed here either.
      integer, intent(in) :: out_unit, n_other
      integer, intent(in) :: naxes_in(:), other_axes(:), other_idx(:)
      integer, intent(in) :: chan_start, chan_len, nx_out, ny_out
      real, intent(in) :: block_data_out(:,:,:)
      integer, intent(inout) :: status

      integer :: fpixels_wr(max_axes), lpixels_wr(max_axes)
      integer :: naxes_wr(max_axes), naxis_wr, k, fitsstat

      if (status.ne.0) return

      naxis_wr = 2 + n_other
      naxes_wr(1) = nx_out
      naxes_wr(2) = ny_out
      fpixels_wr(1) = 1
      fpixels_wr(2) = 1
      lpixels_wr(1) = nx_out
      lpixels_wr(2) = ny_out
      if (n_other.ge.1) then
         naxes_wr(3) = naxes_in(other_axes(1))
         fpixels_wr(3) = chan_start
         lpixels_wr(3) = chan_start + chan_len - 1
      endif
      do k = 2, n_other
         naxes_wr(2+k) = naxes_in(other_axes(k))
         fpixels_wr(2+k) = other_idx(k)
         lpixels_wr(2+k) = other_idx(k)
      enddo

      fitsstat = 0
      call FTPSSE(out_unit, 1, naxis_wr, naxes_wr(1:naxis_wr),&
      &fpixels_wr(1:naxis_wr), lpixels_wr(1:naxis_wr), block_data_out, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write block (planes ', chan_start,&
         &'-', chan_start+chan_len-1, ') to output'
      endif
   end subroutine write_one_block

   subroutine do_block_write(job)
      !! T4d-style: verbatim port of convolve_cubes.f90's own
      !! do_block_write -- see write_convolved_file's own comment there
      !! for the full single-writer-at-a-time rationale. Uses its own
      !! LOCAL status (not write_reprojected_file's own status_par --
      !! that variable is touched by the main thread's `!$omp do` region
      !! concurrently with a background write, so sharing it directly
      !! would race; write_failed is the safe, checked-after-join
      !! handoff instead).
      type(block_write_job_t), intent(inout) :: job
      integer :: status_local

      status_local = 0
      call write_one_block(job%out_unit, job%naxes_in(1:max_axes),&
      &job%other_axes(1:max_axes), job%other_idx(1:max_axes), job%n_other,&
      &job%chan_start, job%chan_len, job%nx, job%ny, job%data, status_local)
      if (status_local.ne.0) then
         write(*,*) 'ERROR: background write failed for planes ',&
         &job%chan_start, '-', job%chan_start+job%chan_len-1
         write_failed = .true.
      endif
   end subroutine do_block_write

   function block_write_thread_entry(arg) bind(C) result(res)
      type(c_ptr), value :: arg
      type(c_ptr) :: res
      type(block_write_job_t), pointer :: job

      call c_f_pointer(arg, job)
      call do_block_write(job)
      res = c_null_ptr
   end function block_write_thread_entry

   subroutine block_write_dispatch_async(job, thread_id, dispatched)
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

   subroutine compose_pix2pix(skymap_from, skyframe_from, skymap_to,&
   &skyframe_to, map_out, status)
      !! Compose pixel_from -> sky_from -> sky_to -> pixel_to, aligning
      !! sky_from/sky_to via a Frame-to-Frame astConvert between the two
      !! SkyFrame objects (handles any axis-order/equinox/system
      !! difference between them; see extract_sky_mapping's own comment
      !! for why this cannot be skipped -- composing the two pixel->sky
      !! Mappings directly, without this alignment step, silently
      !! produces wrong results whenever the two files' SkyFrames present
      !! their axes in a different order).
      integer, intent(in) :: skymap_from, skyframe_from
      integer, intent(in) :: skymap_to, skyframe_to
      integer, intent(out) :: map_out
      integer, intent(inout) :: status

      integer :: sky2sky

      map_out = ast__null
      if (status.ne.0) return

      sky2sky = ast_convert(skyframe_from, skyframe_to, ' ', status)
      if (status.ne.0 .or. sky2sky.eq.ast__null) then
         write(*,*) 'ERROR: failed to align two SkyFrames, status=', status
         status = -1
         return
      endif

      call ast_invert(skymap_to, status)
      map_out = ast_cmpmap(skymap_from, sky2sky, .true., ' ', status)
      map_out = ast_cmpmap(map_out, skymap_to, .true., ' ', status)
      call ast_invert(skymap_to, status)
      if (status.ne.0 .or. map_out.eq.ast__null) then
         write(*,*) 'ERROR: failed to compose the pixel->pixel Mapping,',&
         &' status=', status
         status = -1
      endif

      call ast_annul(sky2sky, status)
   end subroutine compose_pix2pix

   subroutine footprint_bounds(map_from_to, naxes_from, pixaxes_from,&
   &lbnd, ubnd, status)
      !! "from" file's full pixel extent (its own NAXIS on the 2 axes its
      !! sky Mapping depends on), expressed in "to" pixel space, via
      !! astMapBox (the true enclosing bound of each output coordinate --
      !! not a naive 4-corner check, which can underestimate the true
      !! extent for a non-axis-aligned Mapping). map_from_to must be the
      !! forward "from"->"to" pixel Mapping (as returned by
      !! compose_pix2pix with "from" as its first, "to" as its second
      !! argument pair).
      integer, intent(in) :: map_from_to
      integer, intent(in) :: naxes_from(:), pixaxes_from(2)
      double precision, intent(out) :: lbnd(2), ubnd(2)
      integer, intent(inout) :: status

      double precision :: lbnd_in(2), ubnd_in(2), xl(2), xu(2)

      lbnd = 0.0d0
      ubnd = 0.0d0
      if (status.ne.0) return

      lbnd_in(1) = 1.0d0
      lbnd_in(2) = 1.0d0
      ubnd_in(1) = real(naxes_from(pixaxes_from(1)), kind=8)
      ubnd_in(2) = real(naxes_from(pixaxes_from(2)), kind=8)

      call ast_mapbox(map_from_to, lbnd_in, ubnd_in, .true., 1,&
      &lbnd(1), ubnd(1), xl, xu, status)
      call ast_mapbox(map_from_to, lbnd_in, ubnd_in, .true., 2,&
      &lbnd(2), ubnd(2), xl, xu, status)
   end subroutine footprint_bounds

   subroutine extract_sky_mapping(wcs, skymap, skyframe, pixel_axes, status)
      !! Extract the pixel-grid -> sky (RA/Dec) Mapping from a WCS
      !! FrameSet, with the sky axes' positions in the (possibly compound
      !! Stokes+Sky+Spectrum) current frame detected automatically -- no
      !! assumption about which positions they occupy. Every axis PAIR
      !! (not just consecutive ones -- an earlier version assumed a
      !! SkyFrame's 2 axes must be consecutive within a CmpFrame, since
      !! decompose splits a CmpFrame into two axis-contiguous components;
      !! that assumption turned out to be wrong, confirmed by direct
      !! test: a file with RA/Dec on non-adjacent pixel axes 2 and 4 still
      !! has a genuine SkyFrame recoverable via ast_pickaxes(2,4), so
      !! whatever internal structure connects them is not simply "the two
      !! components of one decompose split") is probed with ast_pickaxes
      !! + ast_isaskyframe -- a genuine AST class check, not a guess -- to
      !! find which pair it is. Once known, astMapSplit selects a
      !! Mapping's INPUT axes, but the sky axes are on the OUTPUT side
      !! (current frame) of the pixel->compound Mapping -- so invert first
      !! (making sky axes selectable as inputs), simplify (helps AST
      !! recognise separability), split, then invert the result back to a
      !! forward pixel-subset -> sky Mapping. Also returns the isolated
      !! SkyFrame object itself (not just the Mapping) -- a SkyFrame's own
      !! axis order is NOT a fixed RA-then-Dec convention (it reflects
      !! whichever axis the header declared as longitude vs latitude
      !! first, confirmed by direct comparison: a file with CTYPE1=DEC
      !! presents (Dec,RA), not (RA,Dec)), so the caller needs the actual
      !! Frame object to align two files' sky axes correctly via
      !! astConvert rather than assuming a shared order.
      integer, intent(in) :: wcs
      integer, intent(out) :: skymap, skyframe
      integer, intent(out) :: pixel_axes(2)
      integer, intent(inout) :: status

      integer :: curframe, nout, i, j
      integer :: probe_axes(2), probe_frame, probe_map
      integer :: sky_axes_in(2), out_axes(4)
      integer :: fullmap, simplemap
      logical :: found_sky

      skymap = ast__null
      skyframe = ast__null
      pixel_axes = 0
      if (status.ne.0) return

      nout = ast_geti(wcs, 'Nout', status)
      curframe = ast_getframe(wcs, ast__current, status)
      found_sky = .false.
      outer: do i = 1, nout - 1
         do j = i + 1, nout
            probe_axes(1) = i
            probe_axes(2) = j
            probe_frame = ast_pickaxes(curframe, 2, probe_axes, probe_map, status)
            if (ast_isaskyframe(probe_frame, status)) then
               sky_axes_in = probe_axes
               found_sky = .true.
               skyframe = probe_frame
               call ast_annul(probe_map, status)
               exit outer
            endif
            call ast_annul(probe_frame, status)
            call ast_annul(probe_map, status)
         enddo
      enddo outer
      call ast_annul(curframe, status)

      if (.not. found_sky) then
         write(*,*) 'ERROR: no axis-pair SkyFrame found in the WCS'
         status = -1
         return
      endif

      fullmap = ast_getmapping(wcs, ast__base, ast__current, status)
      call ast_invert(fullmap, status)
      simplemap = ast_simplify(fullmap, status)
      call ast_mapsplit(simplemap, 2, sky_axes_in, out_axes, skymap, status)
      if (status.ne.0 .or. skymap.eq.ast__null) then
         write(*,*) 'ERROR: ast_mapsplit failed to isolate the sky Mapping,',&
         &' status=', status
         status = -1
         call ast_annul(fullmap, status)
         call ast_annul(simplemap, status)
         return
      endif
      call ast_invert(skymap, status)
      pixel_axes = out_axes(1:2)

      call ast_annul(fullmap, status)
      call ast_annul(simplemap, status)
   end subroutine extract_sky_mapping

   subroutine load_wcs(filename, wcs, naxes, status)
      !! Read filename's FITS header via CFITSIO, load it into an AST
      !! FitsChan, and return the WCS FrameSet recovered from it plus this
      !! file's own per-axis pixel-grid extent (NAXISn), needed later to
      !! bound each axis's own footprint for astMapBox.
      !!
      !! The CFITSIO unit comes from fitsio_unit_mod's safe_ftopen --
      !! called concurrently by every OpenMP thread during
      !! write_reprojected_file's per-thread setup, so each thread needs
      !! a genuinely distinct unit with no collision risk. Found the hard
      !! way, twice:
      !! (1) an earlier version left this unit a bare literal (11) with
      !! no per-thread distinction at all, so every thread calling
      !! load_wcs concurrently collided on the very same CFITSIO unit
      !! number, which is what actually crashed/deadlocked multi-threaded
      !! runs -- not, as first suspected, any inherent AST/CFITSIO
      !! thread-safety limitation (a battery of isolated reproducers
      !! matching the real call pattern -- concurrent per-thread AST
      !! FitsChan/FrameSet/Mapping construction, concurrent astResampleR
      !! on independent per-thread Mappings, concurrent CFITSIO
      !! reads+writes on correctly-separated unit numbers -- ran cleanly,
      !! 16 threads x hundreds of iterations x many repeats, once unit
      !! numbers actually stopped colliding). A later, hand-picked
      !! 1000+thread_num offset scheme fixed that.
      !! (2) replacing the hand-picked offset with plain FTGIOU-acquired
      !! units (critical-section-guarded, but with FTOPEN/FTCLOS
      !! themselves OUTSIDE that critical section) reintroduced the exact
      !! same class of crash via a different mechanism: FTGIOU/FTFIOU's
      !! own internal bookkeeping shares mutable global state with
      !! FTOPEN/FTCLOS's own open-file table, so two threads racing a
      !! critical-guarded FTGIOU against another thread's unguarded
      !! FTOPEN still corrupts that shared state. safe_ftopen/safe_ftclos
      !! (fitsio_unit_mod.f90) fix this by keeping FTGIOU+FTOPEN (and
      !! FTCLOS+FTFIOU) inside the SAME critical section -- verified with
      !! an isolated 8-thread reproducer that crashed reliably within a
      !! few hundred iterations without this, and ran clean with it.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: wcs
      integer, intent(out) :: naxes(:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, nkeys, nmore, i, fitschan
      character(len=80) :: card

      wcs = ast__null
      naxes = 0
      if (status.ne.0) return

      fitsstat = 0
      blocksize = 1
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to open FITS file: ', trim(filename)
         call printerror(fitsstat)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      call FTGHSP(unit, nkeys, nmore, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGHSP failed for ', trim(filename)
         call printerror(fitsstat)
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif

      call FTGISZ(unit, size(naxes), naxes, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGISZ failed for ', trim(filename)
         call printerror(fitsstat)
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      if (status.ne.0) then
         call safe_ftclos(unit, fitsstat)
         return
      endif

      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(unit, i, card, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: FTGREC failed at card ', i, ' for ', trim(filename)
            call printerror(fitsstat)
            call safe_ftclos(unit, fitsstat)
            status = -1
            return
         endif
         call ast_putfits(fitschan, card, .false., status)
      enddo
      ! Rewind the channel's internal read pointer before ast_read
      ! consumes it (Card is a 1-based cursor into the card list;
      ! ast_putfits above leaves it sitting past the last card written).
      call ast_seti(fitschan, 'Card', 1, status)

      call safe_ftclos(unit, fitsstat)

      wcs = ast_read(fitschan, status)
      if (status.ne.0 .or. wcs.eq.ast__null) then
         write(*,*) 'ERROR: ast_read failed to recover a WCS FrameSet for ',&
         &trim(filename)
         status = -1
         return
      endif
      if (.not. ast_isaframeset(wcs, status)) then
         write(*,*) 'ERROR: object read from FitsChan is not a FrameSet for ',&
         &trim(filename)
         status = -1
         return
      endif
      call ast_annul(fitschan, status)
   end subroutine load_wcs

   subroutine get_mem_total_kb(mem_total_kb)
      !! Total system RAM in kB, from /proc/meminfo's MemTotal -- same
      !! source and same reasoning as rm_synthesis.f90's own memory
      !! planner (search that file for "Tile planning for memory-
      !! efficient cube processing"): budgeting against TOTAL RAM rather
      !! than instantaneously-available RAM makes the chosen block size
      !! deterministic for a given cube/mem_frac_ram (reproducible across
      !! runs) instead of fluctuating with whatever else the machine is
      !! doing -- with the same caveat rm_synthesis documents: on a busy/
      !! shared node, a large mem_frac_ram can over-commit, since memory
      !! used by other jobs is not subtracted here. 4 GiB fallback if
      !! /proc/meminfo can't be read (e.g. non-Linux), matching
      !! rm_synthesis's own fallback constant exactly.
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

   subroutine print_usage()
      !! Shared by --help/-h and the argument-error path, so the two
      !! can't drift out of sync with each other.
      write(*,'(A)') 'reproject_cubes -- reproject FITS cubes onto a common grid'
      write(*,'(A)') ''
      write(*,'(A)') 'Usage:'
      write(*,'(A)') '  reproject_cubes mode=<intersection|union|reference>'//&
      &' reffile=<reference_file> infiles=<input_file>[,<input_file>...]'
      write(*,'(A)') '  reproject_cubes --config <cfgfile>'
      write(*,'(A)') '  reproject_cubes --config <cfgfile> mode=<...>'//&
      &' [reffile=<...>] [infiles=<...>]'
      write(*,'(A)') '  reproject_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'No positional args -- mode/reffile/infiles must each be given as'//&
      &' key=value (no spaces around ''='', e.g. mode=union), either directly on the'
      write(*,'(A)') 'command line or via --config. If both are given, each CLI key=value'//&
      &' overrides only that same key from --config; unset keys still come from --config.'
      write(*,'(A)') ''
      write(*,'(A)') 'Modes:'
      write(*,'(A)') '  reference    output grid is the reference file''s own extent'
      write(*,'(A)') '  intersection output grid shrinks to the overlap of all inputs'//&
      &' with the reference'
      write(*,'(A)') '  union        output grid grows to cover all inputs and the reference'
      write(*,'(A)') ''
      write(*,'(A)') 'Config file: key=value text file with three required keys:'
      write(*,'(A)') '  mode    = intersection | union | reference'
      write(*,'(A)') '  reffile = /path/to/reference.fits'
      write(*,'(A)') '  infiles = /path/band1.fits,/path/band2.fits,/path/band3.fits'
      write(*,'(A)') '(infiles is comma-separated, no spaces required; ''#'' or '';'''//&
      &' starts a comment.)'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default 0.25):'
      write(*,'(A)') '  mem_frac_ram = fraction (0,0.95] of total system RAM budgeted'//&
      &' for one read/resample/write block of planes at a time -- same concept'
      write(*,'(A)') '  as rm_synthesis''s own mem_frac_ram. Smaller = more, smaller'//&
      &' blocks (less peak memory, more CFITSIO calls); larger = fewer, bigger blocks.'
      write(*,'(A)') '  Threads only parallelise WITHIN one block, never across blocks'//&
      &' (blocks are processed strictly one after another) -- too small a'
      write(*,'(A)') '  mem_frac_ram therefore also throws away most of the OpenMP'//&
      &' speedup, not just increases I/O calls (a printed WARNING flags this).'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default n):'
      write(*,'(A)') '  io_overlap = y|n -- write each block on a background thread,'//&
      &' overlapped with the NEXT block''s own read+resample, instead of blocking'
      write(*,'(A)') '  on the write before starting it. Only one background write is'//&
      &' ever in flight at a time (concurrent CFITSIO handle use from two'
      write(*,'(A)') '  threads is unsafe) -- still overlaps the write with the'//&
      &' following block''s read+compute, which is where the dead time is.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional keys (CLI or config, logging/timing):'
      write(*,'(A)') '  log_level       = error|warn|info|debug (default info)'
      write(*,'(A)') '  timing_enabled  = y|n -- print a stage timing summary (default n)'
      write(*,'(A)') '  log_output_file = path -- append log/timing output to this file'//&
      &' instead of stdout (default empty = stdout)'
   end subroutine print_usage

   subroutine read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs,&
   &mem_frac_ram, io_overlap_l, log_level_l, timing_enabled_l,&
   &log_output_file_l, status)
      !! Parse a --config key=value file: three required keys, mode,
      !! reffile, and infiles (comma-separated, same csv-list convention
      !! rm_synthesis's own multi-band keys use), plus optional keys
      !! mem_frac_ram/io_overlap/log_level/timing_enabled/log_output_file
      !! (all intent(inout) -- left untouched, keeping whatever default
      !! the caller set, if the file doesn't mention them).
      !! Standalone re-implementation of rm_synthesis_mod's
      !! split_key_value/csv_count/csv_get_item (below) rather than a
      !! `use` dependency -- this tool is deliberately kept off the main
      !! rm_synthesis build graph (own binary, own dependency set, see
      !! the Makefile comment), and these are a handful of generic
      !! string-parsing lines each, unlikely to drift.
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(out) :: mode, reffile
      character(len=*), intent(out) :: infiles(:)
      integer, intent(out) :: n_inputs
      real, intent(inout) :: mem_frac_ram
      logical, intent(inout) :: io_overlap_l
      character(len=*), intent(inout) :: log_level_l
      logical, intent(inout) :: timing_enabled_l
      character(len=*), intent(inout) :: log_output_file_l
      integer, intent(out) :: status

      character(len=16384) :: line, val, raw_infiles
      character(len=512) :: key
      integer :: unit_cfg, ios, line_no, j, ios_mfr
      logical :: has_kv, seen_mode, seen_reffile, seen_infiles

      status = 0
      n_inputs = 0
      seen_mode = .false.
      seen_reffile = .false.
      seen_infiles = .false.
      raw_infiles = ' '

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
         select case (trim(key))
         case ('mode')
            mode = trim(val)
            seen_mode = .true.
         case ('reffile')
            reffile = trim(val)
            seen_reffile = .true.
         case ('infiles')
            raw_infiles = trim(val)
            seen_infiles = .true.
         case ('mem_frac_ram')
            read(val, *, iostat=ios_mfr) mem_frac_ram
            if (ios_mfr.ne.0) then
               write(*,*) 'ERROR: mem_frac_ram must be a number, at line ',&
               &line_no, ' in ', trim(cfgfile)
               status = -1
               close(unit_cfg)
               return
            endif
         case ('io_overlap')
            io_overlap_l = flag_from_value_reproject(val)
         case ('log_level')
            log_level_l = trim(val)
         case ('timing_enabled')
            timing_enabled_l = flag_from_value_logging(val)
         case ('log_output_file')
            log_output_file_l = trim(val)
         case default
            write(*,*) 'ERROR: unrecognised config key "', trim(key), '" at line ',&
            &line_no, ' in ', trim(cfgfile)
            status = -1
            close(unit_cfg)
            return
         end select
      enddo
      close(unit_cfg)

      if (.not. seen_mode .or. .not. seen_reffile .or. .not. seen_infiles) then
         write(*,*) 'ERROR: config file ', trim(cfgfile),&
         &' must set mode, reffile, and infiles'
         status = -1
         return
      endif

      n_inputs = cfg_csv_count(raw_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.size(infiles)) then
         write(*,*) 'ERROR: infiles in config must list between 1 and ',&
         &size(infiles), ' files'
         status = -1
         n_inputs = 0
         return
      endif
      do j = 1, n_inputs
         call cfg_csv_get_item(raw_infiles, j, infiles(j))
      enddo
   end subroutine read_reproject_cfg

   subroutine cfg_split_key_value(raw_line, key, val, has_kv)
      !! Same convention as rm_synthesis_mod's split_key_value: strips
      !! ';'/'#' comments, splits on the first '=', blank/comment-only
      !! lines and lines missing either side of '=' yield has_kv=.false.
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
      !! Split a single CLI argument of the form key=value. Deliberately
      !! no '#'/';' comment-stripping (unlike cfg_split_key_value) --
      !! this is one shell-split argv token, not a config-file line, and
      !! a file path could legitimately contain either character.
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

   function flag_from_value_reproject(val) result(flag)
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
   end function flag_from_value_reproject

   function strip_fits_ext(filename) result(base)
      !! Output-name helper: verbatim port of convolve_cubes.f90's own --
      !! see there for the full rationale (avoids the double-extension
      !! "name.fits_REPROJ.FITS" this used to produce). Strips whatever
      !! the trailing extension actually is, not just ".fits".
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

   function cfg_csv_count(str) result(n)
      !! Number of comma-separated items in str (1 if no comma present,
      !! 0 for a blank/empty string).
      character(len=*), intent(in) :: str
      integer :: n
      integer :: i

      n = 0
      if (len_trim(str) == 0) return
      n = 1
      do i = 1, len_trim(str)
         if (str(i:i) == ',') n = n + 1
      enddo
   end function cfg_csv_count

   subroutine cfg_csv_get_item(str, idx, item)
      !! Extract the idx-th (1-based) comma-separated item from str,
      !! trimmed of surrounding blanks.
      character(len=*), intent(in) :: str
      integer, intent(in) :: idx
      character(len=*), intent(out) :: item
      integer :: i, cur, p0, n

      item = ' '
      n = len_trim(str)
      if (n == 0) return

      cur = 1
      p0 = 1
      do i = 1, n
         if (str(i:i) == ',') then
            if (cur == idx) then
               item = adjustl(str(p0:i - 1))
               return
            endif
            cur = cur + 1
            p0 = i + 1
         endif
      enddo
      if (cur == idx) item = adjustl(str(p0:n))
   end subroutine cfg_csv_get_item

end program reproject_cubes
