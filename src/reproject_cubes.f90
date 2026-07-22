! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! Q/U (and optionally I) FITS cubes -- two or more, not tied to the
! multi-band-tomography "band" concept specifically -- onto a single common
! RA/Dec pixel grid using Starlink AST for WCS handling and resampling, so
! the existing rm_synthesis exact-match ingestion (T1) can consume them
! unchanged.
!
! Current stage: minimal proof-of-concept exercising the CFITSIO -> AST
! FitsChan -> FrameSet pipeline on a single input file, to validate the
! Fortran calling convention end-to-end before building the full
! multi-cube reprojection logic (footprint modes, resampling, NaN-fill,
! output write).
!
! Usage: reproject_cubes <fits_file>
program reproject_cubes
   implicit none
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 (`*`-column comments) and cannot be `include`d
   ! into a free-form .f90 file directly (gfortran misparses its comments
   ! as code). Only a handful of symbols are needed here, so they are
   ! declared directly instead, matching AST_PAR's own declared types
   ! exactly (checked against /usr/include/AST_PAR): AST_NULL is a real
   ! exported external procedure (a null source/sink function placeholder
   ! for ast_fitschan), AST__NULL=0 is the "no object" handle value,
   ! AST__SZCHR=200 is the fixed return-string length ast_getc uses.
   external :: ast_null
   integer, parameter :: ast__null = 0
   integer, parameter :: ast__szchr = 200
   integer, parameter :: ast__current = -1
   integer, external :: ast_fitschan, ast_read, ast_geti
   integer, external :: ast_getframe
   logical, external :: ast_isaframeset, ast_isaskyframe, ast_isacmpframe
   character(len=ast__szchr), external :: ast_getc
   character(len=ast__szchr) :: cur_class
   integer :: curframe

   character(len=512) :: infile
   integer :: unit, blocksize, status, fitsstat
   integer :: nkeys, nmore, i
   character(len=80) :: card
   integer :: fitschan, wcs, skyframe
   integer :: nin, nout, sky_nout
   character(len=ast__szchr) :: domain

   ! Stack-based search for the SkyFrame component within (possibly
   ! nested) compound current frames -- see comment at the search site
   ! below for why this replaces the more obvious ast_findframe approach.
   integer :: frame_stack(20), nstack
   integer :: comp1, comp2
   logical :: series, inv1, inv2

   if (command_argument_count() < 1) then
      write(*,*) 'Usage: reproject_cubes <fits_file>'
      stop 1
   endif
   call get_command_argument(1, infile)

   ! --- Read the FITS header cards via CFITSIO ---
   fitsstat = 0
   blocksize = 1
   unit = 11
   call FTOPEN(unit, trim(infile), 0, blocksize, fitsstat)
   if (fitsstat.ne.0) then
      write(*,*) 'ERROR: failed to open FITS file: ', trim(infile)
      call printerror(fitsstat)
      stop 1
   endif

   call FTGHSP(unit, nkeys, nmore, fitsstat)
   if (fitsstat.ne.0) then
      write(*,*) 'ERROR: FTGHSP failed'
      call printerror(fitsstat)
      stop 1
   endif
   write(*,'(A,I0,A)') 'Read header: ', nkeys, ' cards'

   ! --- Build an AST FitsChan and load the cards into it ---
   status = 0
   call ast_begin(status)

   fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
   if (status.ne.0) then
      write(*,*) 'ERROR: ast_fitschan failed, status=', status
      stop 1
   endif

   do i = 1, nkeys
      fitsstat = 0
      call FTGREC(unit, i, card, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGREC failed at card ', i
         call printerror(fitsstat)
         stop 1
      endif
      call ast_putfits(fitschan, card, .false., status)
   enddo
   ! Rewind the channel's internal read pointer before ast_read consumes it
   ! (Card is a 1-based cursor into the card list; ast_putfits above leaves
   ! it sitting past the last card written).
   call ast_seti(fitschan, 'Card', 1, status)

   call FTCLOS(unit, fitsstat)

   if (status.ne.0) then
      write(*,*) 'ERROR: loading cards into FitsChan failed, status=', status
      stop 1
   endif

   ! --- Read the WCS FrameSet out of the FitsChan ---
   wcs = ast_read(fitschan, status)
   if (status.ne.0 .or. wcs.eq.ast__null) then
      write(*,*) 'ERROR: ast_read failed to recover a WCS FrameSet,',&
      &' status=', status
      stop 1
   endif

   if (.not. ast_isaframeset(wcs, status)) then
      write(*,*) 'ERROR: object read from FitsChan is not a FrameSet'
      stop 1
   endif

   nin = ast_geti(wcs, 'Nin', status)
   nout = ast_geti(wcs, 'Nout', status)
   domain = ast_getc(wcs, 'Domain', status)

   write(*,*) 'FrameSet read successfully:'
   write(*,'(A,I0)') '  Nin  (base/pixel axes) : ', nin
   write(*,'(A,I0)') '  Nout (current/wcs axes): ', nout
   write(*,'(A,A)')  '  Domain (current frame)  : ', trim(domain)

   curframe = ast_getframe(wcs, ast__current, status)
   cur_class = ast_getc(curframe, 'Class', status)
   write(*,'(A,A)') '  Class (current frame)   : ', trim(cur_class)

   ! --- Isolate just the 2 sky (RA/Dec) axes from the compound WCS ---
   ! Tried ast_findframe first (the textbook way to pull a SkyFrame out of
   ! a compound WCS) but it came back AST__NULL with status=0 (a genuine
   ! "no match", not an AST error) despite the current frame confirmed
   ! to be a CmpFrame -- likely an argument-marshaling detail of this
   ! Fortran binding not yet understood, not a dead end in principle.
   ! Falling back to the more primitive, more certain approach: CmpFrame
   ! combines two component frames axis-wise (astDecompose splits it back
   ! into those two, plus whether they're combined "in series" -- always
   ! .FALSE. for a CmpFrame, which concatenates axes in parallel, as
   ! opposed to CmpMap's sequential composition); search the (possibly
   ! nested, e.g. Stokes+(Sky+Spectrum)) tree of components with a small
   ! stack until a genuine SkyFrame leaf is found.
   nstack = 1
   frame_stack(1) = curframe
   skyframe = ast__null
   do while (nstack.gt.0 .and. skyframe.eq.ast__null)
      curframe = frame_stack(nstack)
      nstack = nstack - 1
      if (ast_isaskyframe(curframe, status)) then
         skyframe = curframe
      else if (ast_isacmpframe(curframe, status)) then
         call ast_decompose(curframe, comp1, comp2, series, inv1, inv2, status)
         nstack = nstack + 1
         frame_stack(nstack) = comp1
         nstack = nstack + 1
         frame_stack(nstack) = comp2
         call ast_annul(curframe, status)
      else
         call ast_annul(curframe, status)
      endif
   enddo

   if (status.ne.0 .or. skyframe.eq.ast__null) then
      write(*,*) 'ERROR: no SkyFrame component found in the WCS,',&
      &' status=', status
      stop 1
   endif

   sky_nout = ast_geti(skyframe, 'Naxes', status)
   write(*,*) 'Sky sub-frame isolated successfully:'
   write(*,'(A,I0)') '  Naxes (sky axes)        : ', sky_nout

   call ast_annul(skyframe, status)
   call ast_annul(wcs, status)
   call ast_annul(fitschan, status)
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: CFITSIO -> AST FitsChan -> FrameSet pipeline verified.'
end program reproject_cubes
