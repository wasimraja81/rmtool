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
! cube (/data1/tmp/cutout-stokesQ.fits, used by cfg/rmsynth-jennifer.fullim
! .cfg): CASAMBM=T, 288-row BEAMS table, one row per FREQ channel. Route 2:
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
! deconvolved from -- pooled across ALL input files (multi-band support
! needs no extra machinery here: for multiple bands convolved to one common
! resolution, every channel of every band's own per-channel beam is simply
! added to the same pool before finding one common beam, then every file is
! convolved to that single shared target -- matching the physical
! requirement that the target resolution be identical across every channel
! of every band). Can also be set explicitly (target_bmaj/target_bmin/
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
!    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]
!    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
!    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]
!    [npts=<n>] [khachiyan_tol=<tol>]
!    or: convolve_cubes --config <cfgfile>
!    or: convolve_cubes --help | -h
! Full usage text in print_usage below (shared by --help and the
! argument-error path, same convention as reproject_cubes.f90).
program convolve_cubes
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none

   integer, parameter :: max_axes = 10
   integer, parameter :: max_inputs = 50
   integer, parameter :: max_channels = 20000

   character(len=512) :: infiles(max_inputs), beamfiles(max_inputs)
   integer :: n_inputs
   character(len=64) :: outsuffix
   character(len=512) :: badchan_file
   logical :: have_badchan_file
   logical :: have_target
   real(dp) :: target_bmaj, target_bmin, target_bpa
   real(dp) :: max_common_bmaj
   logical :: have_max_common_bmaj
   real :: mem_frac_ram
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

   real(dp), allocatable :: pool_bmaj(:), pool_bmin(:), pool_bpa(:)
   integer :: n_pool
   real(dp) :: common_bmaj, common_bmin, common_bpa
   integer :: badchan_list(max_channels), n_badchan

   call parse_args(status)
   if (status.ne.0) stop 1

   allocate(bmaj_f(max_inputs, max_channels), bmin_f(max_inputs, max_channels))
   allocate(bpa_f(max_inputs, max_channels), isbad_f(max_inputs, max_channels))

   n_badchan = 0
   if (have_badchan_file) then
      call read_badchan_file(badchan_file, badchan_list, n_badchan, status)
      if (status.ne.0) stop 1
   endif

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

   do i = 1, n_inputs
      call write_convolved_file(infiles(i), trim(infiles(i))//trim(outsuffix),&
      &naxis_f(i), sky1_f(i), sky2_f(i), freq_axis_f(i), naxes_f(i,:),&
      &cdelt1_f(i), cdelt2_f(i), nfreq_f(i), bmaj_f(i,1:nfreq_f(i)),&
      &bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)), isbad_f(i,1:nfreq_f(i)),&
      &common_bmaj, common_bmin, common_bpa, mem_frac_ram, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write convolved output for: ', trim(infiles(i))
         stop 1
      endif
      write(*,*) 'OK: wrote ', trim(infiles(i))//trim(outsuffix)
   enddo

   deallocate(bmaj_f, bmin_f, bpa_f, isbad_f)
   write(*,*) 'OK: all inputs convolved to common resolution.'

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

   subroutine parse_args(status)
      integer, intent(out) :: status
      character(len=512) :: this_arg, cli_key, cli_val, cfgfile
      character(len=512) :: raw_infiles, raw_beamfiles
      integer :: argc, iarg, ios
      logical :: has_kv, have_cfgfile, seen_infiles

      status = 0
      n_inputs = 0
      outsuffix = '_CONV.FITS'
      have_badchan_file = .false.
      badchan_file = ' '
      have_target = .false.
      target_bmaj = 0.0d0
      target_bmin = 0.0d0
      target_bpa = 0.0d0
      have_max_common_bmaj = .false.
      max_common_bmaj = 0.0d0
      mem_frac_ram = 0.25
      npts = 2000
      khachiyan_tol = 1.0d-5
      raw_infiles = ' '
      raw_beamfiles = ' '
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
            &raw_beamfiles, seen_infiles, status)
            if (status.ne.0) return
            iarg = iarg + 1
         endif
      enddo

      if (have_cfgfile) call read_cfg_file(cfgfile, raw_infiles, raw_beamfiles,&
      &seen_infiles, status)
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

   subroutine apply_kv(key, val, raw_infiles, raw_beamfiles, seen_infiles, status)
      character(len=*), intent(in) :: key, val
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles
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
         badchan_file = val
         have_badchan_file = .true.
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

   subroutine read_cfg_file(cfgfile, raw_infiles, raw_beamfiles, seen_infiles, status)
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles
      logical, intent(inout) :: seen_infiles
      integer, intent(out) :: status
      character(len=512) :: line, key, val
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
         &seen_infiles, status)
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
      write(*,'(A)') '    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]'
      write(*,'(A)') '    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]'
      write(*,'(A)') '    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]'
      write(*,'(A)') '    [npts=<n>] [khachiyan_tol=<tol>]'
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
      write(*,'(A)') 'badchan_file: same one-integer-per-line, 1-indexed convention as'//&
      &' rm_synthesis''s own global_badchan_file -- these channels, and any channel'
      write(*,'(A)') '  with a degenerate (near-zero) beam entry, are written as all-NaN'//&
      &' planes, not convolved.'
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
      unit = 200
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif

      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.lt.2 .or. naxis.gt.max_axes) then
         write(*,*) 'ERROR: bad or missing NAXIS in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
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
            call FTCLOS(unit, fitsstat)
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
               call FTCLOS(unit, fitsstat)
               return
            endif
         else if (ctype(1:3).eq.'DEC') then
            if (sky2.eq.0) then
               sky2 = k
            else
               write(*,*) 'ERROR: more than one DEC-like axis in: ', trim(filename)
               status = -1
               call FTCLOS(unit, fitsstat)
               return
            endif
         else if (ctype(1:4).eq.'FREQ') then
            freq_axis = k
         endif
      enddo

      if (sky1.eq.0 .or. sky2.eq.0) then
         write(*,*) 'ERROR: could not identify RA/DEC sky axes in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif
      if (freq_axis.eq.0) then
         write(*,*) 'ERROR: could not identify a FREQ axis in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
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
               call FTCLOS(unit, fitsstat)
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
         call FTCLOS(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the DEC axis in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      call check_no_rotation(unit, sky1, sky2, filename, status)
      if (status.ne.0) then
         call FTCLOS(unit, fitsstat)
         return
      endif

      call FTCLOS(unit, fitsstat)
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
      unit = 201
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif

      fitsstat = 0
      call FTMNHD(unit, -1, 'BEAMS', 0, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: no BEAMS binary table extension found in: ',&
         &trim(filename), ' -- pass an ASCII beamfile instead (see --help)'
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      fitsstat = 0
      call FTGNRW(unit, nrows, fitsstat)
      if (fitsstat.ne.0 .or. nrows.lt.1) then
         write(*,*) 'ERROR: could not read row count of BEAMS table in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
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
         call FTCLOS(unit, fitsstat)
         return
      endif

      allocate(rb_bmaj(nrows), rb_bmin(nrows), rb_bpa(nrows), rb_chan(nrows))
      fitsstat = 0
      call FTGCVE(unit, col_bmaj, 1, 1, nrows, 0.0, rb_bmaj, anyflag, fitsstat)
      call FTGCVE(unit, col_bmin, 1, 1, nrows, 0.0, rb_bmin, anyflag, fitsstat)
      call FTGCVE(unit, col_bpa, 1, 1, nrows, 0.0, rb_bpa, anyflag, fitsstat)
      call FTGCVJ(unit, col_chan, 1, 1, nrows, 0, rb_chan, anyflag, fitsstat)
      call FTCLOS(unit, fitsstat)
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
   &isbad, tgt_bmaj, tgt_bmin, tgt_bpa, mem_frac_ram, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads
      use gaussft_mod, only: plan_convolution, convolve_to_beam, destroy_convolution_plan
      character(len=*), intent(in) :: infile, outfile
      integer, intent(in) :: naxis, sky1, sky2, freq_axis, naxes(max_axes)
      real(dp), intent(in) :: cdelt1, cdelt2
      integer, intent(in) :: nfreq
      real(dp), intent(in) :: bmaj_in(nfreq), bmin_in(nfreq), bpa_in(nfreq)
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: tgt_bmaj, tgt_bmin, tgt_bpa
      real, intent(in) :: mem_frac_ram
      integer, intent(out) :: status

      integer :: nx, ny, in_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: naxes_out(max_axes)
      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes, block_planes64
      integer :: block_planes, chan_start, chan_len, local_iplane, nthreads
      real, allocatable :: block_in(:,:,:), block_out(:,:,:)
      real(dp) :: bpa_in_pixel(nfreq), tgt_bpa_pixel
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: status_par, ich, k
      real(dp) :: dx_deg, dy_deg
      real(dp) :: nanval
      ! Automatic (stack) arrays sized directly from the dummy arguments
      ! naxes(sky1)/naxes(sky2) -- NOT from the local nx/ny copies below,
      ! which aren't assigned until the executable part runs, too late to
      ! bound a specification-part automatic array. One instance per
      ! OpenMP thread (see the parallel region's own private() clause),
      ! needed only because gaussft_mod's convolve_to_beam works in
      ! real(dp), while block_in/block_out (this file's own I/O buffers)
      ! are single precision, matching every other block-I/O buffer in
      ! this project (reproject_cubes.f90's own block_data_in/out).
      real(dp) :: plane_in(naxes(sky1), naxes(sky2))
      real(dp) :: plane_out(naxes(sky1), naxes(sky2))

      status = 0
      nx = naxes(sky1)
      ny = naxes(sky2)
      dx_deg = abs(cdelt1)
      dy_deg = abs(cdelt2)
      do ich = 1, nfreq
         call sky_to_pixel_bpa(bpa_in(ich), cdelt1, cdelt2, bpa_in_pixel(ich))
      enddo
      call sky_to_pixel_bpa(tgt_bpa, cdelt1, cdelt2, tgt_bpa_pixel)

      fitsstat = 0
      in_unit = 210
      call FTOPEN(in_unit, trim(infile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot reopen input for output header: ', trim(infile)
         status = -1
         return
      endif

      out_unit = 211
      fitsstat = 0
      call FTINIT(out_unit, '!'//trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create output file: ', trim(outfile)
         status = -1
         call FTCLOS(in_unit, fitsstat)
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
      call FTCLOS(in_unit, fitsstat)

      block_planes64 = 0
      call get_mem_total_kb(mem_total_kb)
      bytes_per_plane = int(4,8) * int(nx,8) * int(ny,8) * 2_8
      mem_safe_bytes = int(real(mem_frac_ram,8) * real(mem_total_kb,8) * 1024.0d0, 8)
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes = int(min(block_planes64, int(nfreq,8)))
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ', nfreq,&
      &' plane(s), in blocks of up to ', block_planes, ' plane(s)'
      if (block_planes.lt.omp_get_max_threads()) then
         write(*,'(A,I0,A,I0,A)') 'WARNING: mem_frac_ram limits blocks to ',&
         &block_planes, ' plane(s), below the ', omp_get_max_threads(),&
         &' threads available -- parallelism is reduced; raise mem_frac_ram',&
         &' for full speedup if memory allows.'
      endif

      allocate(block_in(nx, ny, block_planes))
      allocate(block_out(nx, ny, block_planes))

      call plan_convolution(nx, ny, plan_fwd, plan_bwd)

      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      chan_start = 1
      do while (chan_start.le.nfreq)
         chan_len = min(block_planes, nfreq - chan_start + 1)

         call read_freq_block(infile, naxis, sky1, sky2, freq_axis,&
         &naxes, chan_start, chan_len, nx, ny, block_in(:,:,1:chan_len), status_par)
         if (status_par.ne.0) exit

         !$omp parallel do num_threads(nthreads) schedule(dynamic)&
         !$omp& default(none) shared(chan_len, nx, ny, block_in, block_out,&
         !$omp& isbad, chan_start, plan_fwd, plan_bwd, dx_deg, dy_deg,&
         !$omp& bmaj_in, bmin_in, bpa_in_pixel, tgt_bmaj, tgt_bmin,&
         !$omp& tgt_bpa_pixel, status_par)&
         !$omp& private(local_iplane, ich, k, nanval, plane_in, plane_out)
         do local_iplane = 1, chan_len
            ich = chan_start + local_iplane - 1
            if (isbad(ich)) then
               nanval = ieee_value(1.0_dp, ieee_quiet_nan)
               block_out(:,:,local_iplane) = real(nanval)
            else
               plane_in = real(block_in(:,:,local_iplane), dp)
               call convolve_to_beam(plan_fwd, plan_bwd,&
               &plane_in, nx, ny,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel, plane_out, k)
               block_out(:,:,local_iplane) = real(plane_out)
               if (k.ne.0) then
                  !$omp atomic write
                  status_par = -1
               endif
            endif
         enddo
         !$omp end parallel do
         if (status_par.ne.0) exit

         call write_freq_block(out_unit, naxis, sky1, sky2, freq_axis,&
         &naxes, chan_start, chan_len, nx, ny, block_out(:,:,1:chan_len), status_par)
         if (status_par.ne.0) exit

         chan_start = chan_start + chan_len
      enddo

      call destroy_convolution_plan(plan_fwd, plan_bwd)
      deallocate(block_in, block_out)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to convolve/write one or more planes for: ', trim(infile)
         status = -1
         call FTCLOS(out_unit, fitsstat)
         return
      endif

      call FTCLOS(out_unit, fitsstat)
   end subroutine write_convolved_file

   subroutine copy_generic_header_convolve(src_unit, dst_unit, status)
      !! No reprojection happens here (unlike reproject_cubes) -- input
      !! and output share IDENTICAL axis layout and numbering, so every
      !! header card copies through verbatim EXCEPT the structural
      !! keywords FTPHPR already wrote, and BMAJ/BMIN/BPA/CASAMBM (the
      !! caller overwrites BMAJ/BMIN/BPA with the new common beam right
      !! after calling this, and never writes CASAMBM at all: the output
      !! has one uniform beam, not a per-channel BEAMS table, and this
      !! program never creates a BEAMS extension in its output, so a
      !! stale CASAMBM=T left over from the input would be actively
      !! misleading to a downstream reader).
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
      unit = 5100
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
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
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') from ', trim(filename)
         status = -1
      endif
   end subroutine read_freq_block

   subroutine write_freq_block(out_unit, naxis, sky1, sky2, freq_axis,&
   &naxes, chan_start, chan_len, nx, ny, block_data, status)
      !! Output axis layout is IDENTICAL to the input's own (no
      !! reprojection, unlike reproject_cubes' write_one_block) -- so
      !! this writes directly at the SAME axis numbering/order as the
      !! input, no sky-first canonicalisation needed.
      integer, intent(in) :: out_unit, naxis, sky1, sky2, freq_axis, naxes(max_axes)
      integer, intent(in) :: chan_start, chan_len, nx, ny
      real, intent(in) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: fpixels_wr(max_axes), lpixels_wr(max_axes), fitsstat

      if (status.ne.0) return

      fpixels_wr(1:naxis) = 1
      lpixels_wr(1:naxis) = 1
      lpixels_wr(sky1) = nx
      lpixels_wr(sky2) = ny
      fpixels_wr(freq_axis) = chan_start
      lpixels_wr(freq_axis) = chan_start + chan_len - 1

      fitsstat = 0
      call FTPSSE(out_unit, 1, naxis, naxes(1:naxis), fpixels_wr(1:naxis),&
      &lpixels_wr(1:naxis), block_data, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') to output'
      endif
   end subroutine write_freq_block

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

end program convolve_cubes
