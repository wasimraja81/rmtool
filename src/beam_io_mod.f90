module beam_io_mod
   !! Per-channel beam (PSF) reading shared by convolve_cubes,
   !! match_cubes and the reprojection beam-sampling check in
   !! target_grid_mod. The three readers lived as identical copies in
   !! convolve_cubes.f90 and match_cubes.f90 until UR-21
   !! (docs/dev/BRANCH_STATUS_update-reproject.md), and had already
   !! drifted (only one printed the hint for a malformed beam log); one
   !! copy here keeps the tools reading beams the same way.
   !!
   !! Sources: a CASA-style BEAMS binary table (EXTNAME='BEAMS', columns
   !! BMAJ/BMIN/BPA in arcsec/arcsec/deg, CHAN 0-indexed), or a plain
   !! text beam log, one line per channel. A degenerate entry (BMAJ or
   !! BMIN below 1e-6 arcsec; CASA's own placeholder is ~1.18e-38) marks
   !! that channel bad rather than failing the read.
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use fitsio_unit_mod
   implicit none
   private
   public :: read_beams, read_beams_table, read_beams_ascii, infile_beam_extent

contains

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
      !! confirmed against a real ASKAP cube (see convolve_cubes.f90's
      !! top comment).
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
      !! (1-indexed channel, '#'-prefixed/blank lines skipped). A channel never listed is flagged bad
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

   subroutine infile_beam_extent(filename, have_beam, bmin_min_arcsec,&
   &bmaj_max_arcsec, status)
      !! The narrowest minor axis and the widest major axis of
      !! filename's beam, in arcsec, over every good channel: from its
      !! BEAMS table when the header says CASAMBM=T, otherwise from the
      !! header's BMAJ/BMIN (degrees). have_beam is .false. when neither
      !! is present or no channel has a usable beam -- the caller decides
      !! what an unknown beam means, this routine does not guess one.
      character(len=*), intent(in) :: filename
      logical, intent(out) :: have_beam
      real(dp), intent(out) :: bmin_min_arcsec, bmaj_max_arcsec
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, nrows, k
      integer :: col_bmaj, col_bmin
      real, allocatable :: rb_bmaj(:), rb_bmin(:)
      logical :: anyflag, casambm
      double precision :: hdr_bmaj, hdr_bmin
      character(len=80) :: comment

      have_beam = .false.
      bmin_min_arcsec = 0.0_dp
      bmaj_max_arcsec = 0.0_dp
      status = 0
      fitsstat = 0
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      casambm = .false.
      fitsstat = 0
      call FTGKYL(unit, 'CASAMBM', casambm, comment, fitsstat)
      if (fitsstat.ne.0) casambm = .false.

      if (casambm) then
         fitsstat = 0
         call FTMNHD(unit, -1, 'BEAMS', 0, fitsstat)
         if (fitsstat.eq.0) call FTGNRW(unit, nrows, fitsstat)
         if (fitsstat.eq.0) call FTGCNO(unit, .false., 'BMAJ', col_bmaj, fitsstat)
         if (fitsstat.eq.0) call FTGCNO(unit, .false., 'BMIN', col_bmin, fitsstat)
         if (fitsstat.eq.0 .and. nrows.ge.1) then
            allocate(rb_bmaj(nrows), rb_bmin(nrows))
            call FTGCVE(unit, col_bmaj, 1, 1, nrows, 0.0, rb_bmaj, anyflag, fitsstat)
            call FTGCVE(unit, col_bmin, 1, 1, nrows, 0.0, rb_bmin, anyflag, fitsstat)
            if (fitsstat.eq.0) then
               do k = 1, nrows
                  if (rb_bmaj(k).lt.1.0e-6 .or. rb_bmin(k).lt.1.0e-6) cycle
                  if (.not. have_beam) then
                     bmin_min_arcsec = real(rb_bmin(k), dp)
                     bmaj_max_arcsec = real(rb_bmaj(k), dp)
                     have_beam = .true.
                  else
                     bmin_min_arcsec = min(bmin_min_arcsec, real(rb_bmin(k), dp))
                     bmaj_max_arcsec = max(bmaj_max_arcsec, real(rb_bmaj(k), dp))
                  endif
               enddo
            endif
            deallocate(rb_bmaj, rb_bmin)
         endif
      else
         fitsstat = 0
         call FTGKYD(unit, 'BMAJ', hdr_bmaj, comment, fitsstat)
         call FTGKYD(unit, 'BMIN', hdr_bmin, comment, fitsstat)
         if (fitsstat.eq.0 .and. hdr_bmaj.gt.0.0d0 .and. hdr_bmin.gt.0.0d0) then
            bmaj_max_arcsec = real(hdr_bmaj, dp) * 3600.0_dp
            bmin_min_arcsec = real(hdr_bmin, dp) * 3600.0_dp
            have_beam = .true.
         endif
      endif
      fitsstat = 0
      call safe_ftclos(unit, fitsstat)
   end subroutine infile_beam_extent

end module beam_io_mod
