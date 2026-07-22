! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! FITS cubes -- two or more, not tied to the multi-band-tomography "band"
! concept specifically -- onto a single reference cube's grid using
! Starlink AST for WCS handling and resampling, so the existing
! rm_synthesis exact-match ingestion (T1) can consume genuinely misaligned
! bands unchanged.
!
! Current stage: proof-of-concept for extracting a clean pixel->sky
! Mapping from a single file's (possibly compound Stokes+Sky+Spectrum) WCS
! via astMapSplit, with the sky axes' positions detected automatically --
! no assumption about which positions they occupy. Earlier attempts to
! infer this by walking astDecompose's Frame-splitting order turned out
! unreliable (the compound Frame's own component order does not
! necessarily match the axis order the WCS actually presents -- confirmed
! by direct transform-vs-CRVAL comparison). astConvert/astFindFrame
! domain-search were also tried, but per the actual SUN/211 manual
! (installed via libstarlink-ast-doc), domain search only examines each
! FrameSet's own top-level registered frames (base + current) by their
! own Domain attribute -- it does not recurse into a CmpFrame's internal
! components, so a domainlist of 'SKY' can never match here (the current
! frame's own Domain is the whole compound label "STOKES-SKY-SPECTRUM").
! This stage instead uses a genuine AST class check (ast_isaskyframe) on
! every consecutive axis pair (guaranteed consecutive since CmpFrame
! concatenates component axis blocks contiguously) to find the sky axes'
! true positions before handing them to astMapSplit.
!
! Usage: reproject_cubes <fits_file>
program reproject_cubes
   implicit none
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
   integer, external :: ast_pickaxes
   logical, external :: ast_isaframeset, ast_isaskyframe
   character(len=ast__szchr), external :: ast_getc

   character(len=512) :: infile
   integer :: wcs, fullmap, simplemap, skymap
   integer :: nin, nout, status
   integer :: sky_axes_in(2), out_axes(4)
   double precision :: xin(1), yin(1), xout(1), yout(1)
   integer :: curframe, probe_axes(2), probe_frame, probe_map, i
   logical :: found_sky

   if (command_argument_count() < 1) then
      write(*,*) 'Usage: reproject_cubes <fits_file>'
      stop 1
   endif
   call get_command_argument(1, infile)

   status = 0
   call ast_begin(status)

   call load_wcs(infile, wcs, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to load WCS FrameSet'
      stop 1
   endif

   nin = ast_geti(wcs, 'Nin', status)
   nout = ast_geti(wcs, 'Nout', status)
   write(*,'(A,I0,A)') 'Nin=', nin, ' pixel axes'

   ! --- Locate the sky axes automatically, with no assumption about
   ! which positions they occupy ---
   ! A SkyFrame is always an indivisible 2-axis unit within a CmpFrame
   ! (component Frames are concatenated as contiguous axis blocks, never
   ! interleaved), so its 2 axes are guaranteed consecutive -- but WHICH
   ! consecutive pair depends on the file (confirmed unreliable to infer
   ! from astDecompose's own traversal order, see module comment above).
   ! Checking every consecutive pair via ast_pickaxes + ast_isaskyframe
   ! is a genuine AST class check, not a guess: it tries each candidate
   ! and asks AST directly "is this actually a SkyFrame", stopping at the
   ! first (and structurally only possible) match.
   curframe = ast_getframe(wcs, ast__current, status)
   found_sky = .false.
   do i = 1, nout - 1
      probe_axes(1) = i
      probe_axes(2) = i + 1
      probe_frame = ast_pickaxes(curframe, 2, probe_axes, probe_map, status)
      if (ast_isaskyframe(probe_frame, status)) then
         sky_axes_in = probe_axes
         found_sky = .true.
         call ast_annul(probe_frame, status)
         exit
      endif
      call ast_annul(probe_frame, status)
   enddo
   call ast_annul(curframe, status)

   if (.not. found_sky) then
      write(*,*) 'ERROR: no consecutive-axis-pair SkyFrame found in the WCS'
      stop 1
   endif
   write(*,'(A,I0,A,I0)') 'Sky axes detected automatically at positions: ',&
   &sky_axes_in(1), ' , ', sky_axes_in(2)

   ! --- Extract the pixel-grid -> sky Mapping ---
   ! astMapSplit selects a Mapping's INPUT axes, but the sky axes'
   ! positions (just determined above) are on the OUTPUT side (current
   ! frame) of the pixel->compound Mapping -- so invert first (making sky
   ! axes selectable as inputs), simplify (helps AST recognise
   ! separability), split, then invert the result back to a forward
   ! pixel-subset -> sky Mapping.
   fullmap = ast_getmapping(wcs, ast__base, ast__current, status)
   call ast_invert(fullmap, status)
   simplemap = ast_simplify(fullmap, status)
   call ast_mapsplit(simplemap, 2, sky_axes_in, out_axes, skymap, status)
   if (status.ne.0 .or. skymap.eq.ast__null) then
      write(*,*) 'ERROR: ast_mapsplit failed to isolate the sky Mapping,',&
      &' status=', status
      stop 1
   endif
   call ast_invert(skymap, status)

   write(*,*) 'Pixel -> sky Mapping extracted:'
   write(*,'(A,I0)') '  Nin  (pixel axes it depends on): ',&
   &ast_geti(skymap, 'Nin', status)
   write(*,'(A,I0)') '  Nout (sky axes)                : ',&
   &ast_geti(skymap, 'Nout', status)
   write(*,'(A,I0,A,I0)') '  Corresponding pixel axis position(s): ',&
   &out_axes(1), ' , ', out_axes(2)

   if (ast_geti(skymap,'Nin',status).eq.2 .and.&
   &ast_geti(skymap,'Nout',status).eq.2) then
      xin(1) = 5.0d0
      yin(1) = 9.0d0
      call ast_tran2(skymap, 1, xin, yin, .true., xout, yout, status)
      write(*,'(A,F0.6,A,F0.6,A)') '  Pixel (5,9) -> sky (rad): (',&
      &xout(1), ' , ', yout(1), ')'
      write(*,*) '  compare to full pixel(5,9,3,1) -> world result',&
      &' from the earlier full-map diagnostic (-3.141383, -0.000140)'
   endif

   call ast_annul(skymap, status)
   call ast_annul(simplemap, status)
   call ast_annul(fullmap, status)
   call ast_annul(wcs, status)
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: pixel -> sky Mapping extraction verified.'

contains

   subroutine load_wcs(filename, wcs, status)
      !! Read filename's FITS header via CFITSIO, load it into an AST
      !! FitsChan, and return the WCS FrameSet recovered from it.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: wcs
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, nkeys, nmore, i, fitschan
      character(len=80) :: card

      wcs = ast__null
      if (status.ne.0) return

      fitsstat = 0
      blocksize = 1
      unit = 11
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to open FITS file: ', trim(filename)
         call printerror(fitsstat)
         status = -1
         return
      endif

      call FTGHSP(unit, nkeys, nmore, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGHSP failed for ', trim(filename)
         call printerror(fitsstat)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      if (status.ne.0) then
         call FTCLOS(unit, fitsstat)
         return
      endif

      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(unit, i, card, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: FTGREC failed at card ', i, ' for ', trim(filename)
            call printerror(fitsstat)
            call FTCLOS(unit, fitsstat)
            status = -1
            return
         endif
         call ast_putfits(fitschan, card, .false., status)
      enddo
      ! Rewind the channel's internal read pointer before ast_read
      ! consumes it (Card is a 1-based cursor into the card list;
      ! ast_putfits above leaves it sitting past the last card written).
      call ast_seti(fitschan, 'Card', 1, status)

      call FTCLOS(unit, fitsstat)

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

end program reproject_cubes
