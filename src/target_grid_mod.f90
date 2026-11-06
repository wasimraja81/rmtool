module target_grid_mod
   !! Resolves target_system=/target_projection=/target_center=/
   !! target_lon=/target_lat=/target_cdelt_lon=/target_cdelt_lat=
   !! (docs/dev/REPROJECT_TARGET_FRAME_PLAN.md, docs/dev/
   !! BRANCH_STATUS_update-reproject.md tickets UR-7/UR-8) shared by
   !! reproject_cubes.f90 and match_cubes.f90's own independent
   !! reffile=-consuming main programs, so the two tools' own target-
   !! grid handling cannot diverge the way AST's SincSinc/SincCos
   !! drifted from each other -- same reasoning as ast_interp_mod.f90.
   !!
   !! `reffile=` supersedes everything about the output grid -- full
   !! stop. When given, none of the target_* keys below may also be
   !! given (hard error, not a silent ignore); the caller's own
   !! existing, unmodified `load_wcs(reffile,...)` + `extract_sky_
   !! mapping` remains the entire story, exactly as before any of this
   !! module existed. The only thing this module does in that case is
   !! reject the combination.
   !!
   !! When `reffile=` is absent, this module builds a fresh WCS
   !! template for the caller to use exactly as if it had come from a
   !! file, entirely IN MEMORY -- no interim file is written at any
   !! point (confirmed unnecessary, 2026-10-13: write_reprojected_
   !! file's own header-copy step only ever reads CTYPE/CRVAL/CRPIX/
   !! CDELT/PC/EQUINOX/RADESYS for the 2 sky axes from "reffile" --
   !! never any other metadata, which always comes from infile -- so
   !! there was never a real reason to round-trip through a file just
   !! to get AST to compute those 7 values correctly). Mechanism,
   !! verified end-to-end before any of this was written (scratch/
   !! build_check/test_noref_template.f90, test_shiftmap.f90):
   !!
   !! 1. Resolve the target field centre (lon0,lat0 in target_system):
   !!    `infiles_centroid` (default when target_center is omitted) --
   !!    3-D-vector-averaged centroid of each infile's own pixel-grid-
   !!    centre sky position (not naive lon/lat averaging, which breaks
   !!    near the poles or across the 0/360 wrap); or `user_defined`
   !!    via target_lon=/target_lat=, parsed with the same flexible
   !!    notation (decimal degrees, colon-sexagesimal, letter-suffixed
   !!    h/d/m/s sexagesimal) as scripts/casda_fetch.py's own --lon/
   !!    --lat.
   !! 2. Resolve CDELT magnitude for both axes: explicit target_cdelt_
   !!    lon=/target_cdelt_lat= (given together or not at all), or, when
   !!    both are omitted, the smallest |CDELT| across the infiles, taken
   !!    per axis (smallest_infile_cdelt) -- one infile gives its own
   !!    scale, several differing ones give the finest.
   !! 3. Build a TRIAL pixel->sky Mapping via a FitsChan populated
   !!    entirely in memory (ast_putfits+ast_read on hand-built card
   !!    text -- no CFITSIO file I/O at any point for this step),
   !!    CRPIX=(1,1), CRVAL=(lon0,lat0), the resolved CDELT, CTYPE for
   !!    the target system and projection (target_projection=, SIN when
   !!    omitted). AST builds the projection from the CTYPE suffix, so
   !!    every projection travels through this same in-memory path;
   !!    this was checked against astropy's reproject_interp for each
   !!    of TAN/SIN/CAR/AIT/HPX, on and off the equator, before the key
   !!    was enabled.
   !! 4. For each infile: reuse `footprint_bounds`/`compose_pix2pix`
   !!    (the SAME machinery `mode=union` already uses, verbatim
   !!    adapted here -- see those subroutines below) to find its own
   !!    footprint in trial-pixel space, accumulating the union exactly
   !!    as the existing mode=union loop does.
   !! 5. Round to the integer pixel grid, derive the final NAXIS and
   !!    CRPIX via the exact shift formula copy_axis_keywords already
   !!    uses elsewhere in this project (crpix_new = crpix_old -
   !!    crpix_shift).
   !! 6. Build the FINAL pixel->sky Mapping by composing an AST
   !!    ShiftMap (shift = 1 - final_crpix, confirmed by direct test to
   !!    reproduce the file-based ast_write result exactly) with the
   !!    trial Mapping -- no new file, no ast_wcsmap_ construction
   !!    needed; `ast_write`, when the caller is ready to materialise a
   !!    header, serializes this exactly as it would any other Mapping
   !!    (confirmed: CRVAL/CRPIX come out correct to floating-point
   !!    precision from this in-memory-built Mapping, same as from one
   !!    loaded from an existing file).
   !!
   !! resolve_target_grid itself returns PLAIN VALUES for the no-
   !! reffile case (target system name, centre lon/lat in degrees,
   !! CDELT magnitudes, final CRPIX) -- NOT live AST object handles.
   !! Confirmed the hard way, by direct multi-threaded test (scratch/
   !! build_check/test_ast_copy_threads.f90): AST enforces per-thread
   !! object OWNERSHIP on EVERY operation, including ast_copy itself,
   !! not only on later use of a copied object (reproject_cubes.f90's
   !! own write_reprojected_file has a detailed comment on a related,
   !! previously-discovered failure mode -- a DIFFERENT thread handing
   !! out copies; this module's own first attempt hit a stricter one:
   !! even a thread copying an object it does not itself own fails,
   !! "Invalid Object pointer ... owned by another thread"). So no AST
   !! object built here may ever be touched by a different thread than
   !! the one that built it -- every caller that needs the actual
   !! Mapping/SkyFrame (write_reprojected_file's own per-thread setup,
   !! and write_target_wcs_cards below) calls build_target_skymap
   !! (public) to build it FRESH, itself, from these plain values,
   !! exactly mirroring how the existing reffile-based path already
   !! has each thread build its own Mapping from scratch via its own
   !! load_wcs+extract_sky_mapping rather than sharing one.
   !!
   !! build_target_skymap (public): rebuilds the exact same final
   !! pixel->sky Mapping resolve_target_grid derived, from the plain
   !! values alone -- safe to call from any thread, independently,
   !! each call producing its own thread-owned objects.
   !!
   !! write_target_wcs_cards (public, called from each tool's own
   !! write_reprojected_file in place of the file-based copy_axis_
   !! keywords/copy_sky_rotation_matrix/copy_wcs_system_keywords calls,
   !! only when resolve_target_grid reports used_target_grid): calls
   !! build_target_skymap itself (so it is likewise safe regardless of
   !! which thread OpenMP's own `single` construct happens to run it
   !! on), `ast_write`s the result into a throwaway in-memory FitsChan
   !! (Encoding forced to 'FITS-WCS' -- without it AST silently falls
   !! back to its own Native encoding, found the hard way earlier this
   !! session), and dumps only the WCS-relevant cards directly onto the
   !! caller's already-open output unit via FTPREC -- safe because that
   !! unit, at the point this is called, has only the structural cards
   !! FTPHPR already wrote, nothing WCS-related yet, so there is no
   !! duplicate-card risk.
   use fitsio_unit_mod
   use beam_io_mod, only: infile_beam_extent
   implicit none
   private
   public :: resolve_target_grid, write_target_wcs_cards, build_target_skymap,&
   &check_output_grid_validity, check_overlap_and_size_grid, check_disk_space,&
   &output_cell_size_arcsec, check_beam_sampling, common_beam_sampling,&
   &fits_reference_cell_arcsec, cell_equal_rtol

   integer, parameter :: ast__null = 0
   integer, parameter :: ast__base = 0
   integer, parameter :: ast__current = -1
   integer, parameter :: max_axes_tg = 10
   ! AST's own invalid-point sentinel (/usr/include/AST_PAR), not a
   ! tolerance comparison -- ast_tran2 returns this exact value for
   ! either coordinate when a point has no valid mapping (UR-10,
   ! docs/dev/BRANCH_STATUS_update-reproject.md; confirmed this is
   ! AST's own behaviour, not assumed, by REPROJECT_TARGET_FRAME_
   ! PLAN.md section 4c's own HPX facet-seam verification).
   double precision, parameter :: ast__bad_tg = -1.79769313486231571d308
   ! AST__SZCHR in /usr/include/AST_PAR: the declared length of every
   ! string AST_GETC returns. This module does not include AST_PAR, so
   ! the value is mirrored here, as reproject_cubes.f90 already does.
   integer, parameter :: ast__szchr_tg = 200
   ! Tolerance, in pixels, for treating a footprint bound as lying on a
   ! whole pixel before it is rounded; see check_overlap_and_size_grid.
   double precision, parameter :: bound_snap_tol = 1.0d-6
   ! Beam sampling on the output grid (UR-21): a beam's minor axis must
   ! span at least beam_refuse_cells of the larger output cell side, or
   ! the reprojection is refused; below beam_warn_cells it is warned as
   ! barely Nyquist. Interpolation reproduces what the output grid can
   ! sample and aliases what it cannot. Checked against BMIN alone, so
   ! the beam's PA does not enter: a Gaussian is nowhere narrower than
   ! its minor axis.
   double precision, parameter :: beam_refuse_cells = 2.0d0
   double precision, parameter :: beam_warn_cells = 2.5d0
   ! Relative tolerance for two grids' pixel AREAS (not linear cell
   ! sizes -- see check_beam_sampling) to count as equal. Both are
   ! measured at their own reference pixel, where a projection adds no
   ! distortion. For TAN/SIN/CAR/AIT this would agree to floating-point
   ! round-off, but AST's own cell-size measurement (a finite one-pixel
   ! chord, not an infinitesimal derivative) picks up a real, bounded,
   ! projection-curvature term for HPX/AIT -- confirmed independently
   ! against astropy's own WCS (same finite-chord method) across every
   ! system/case combination, worst case 1.058e-8 relative (HPX,
   ! rotated). cell_equal_rtol sits ~100x above that confirmed ceiling,
   ! and ~5+ orders of magnitude below any coarsening (even 1%) this
   ! check exists to catch.
   double precision, parameter :: cell_equal_rtol = 1.0d-6
   ! Largest departure, in pixels, of a one-pixel step through the
   ! infile->trial mapping from (1,0)/(0,1) for the two grids to count
   ! as differing by a pure shift; see infile_lattice_offset. Identical
   ! geometries give round-off (order 1e-11), a rotation by angle a
   ! gives about a.
   double precision, parameter :: lattice_shift_tol = 1.0d-6

contains

   subroutine resolve_target_grid(seen_reffile, target_system, target_projection,&
   &target_center, target_lon, target_lat, target_cdelt_lon, target_cdelt_lat,&
   &target_hpx_nside, infiles, n_inputs, used_target_grid, grid_system,&
   &grid_projection, center_lon_deg, center_lat_deg,&
   &cdelt1_mag, cdelt2_mag, final_crpix, pixaxes_ref, naxes_ref,&
   &n_planes_per_infile, status)
      !! Entry point, called once by each tool's own main program,
      !! after CLI/cfg parsing and before the existing alignment
      !! pipeline runs. When seen_reffile is true: validates that no
      !! target_* key is also given (hard error if one is), sets
      !! used_target_grid=.false., and returns immediately -- the
      !! caller proceeds with its own existing load_wcs(reffile,...)
      !! and its own direct call to check_overlap_and_size_grid, below,
      !! unchanged; n_planes_per_infile is NOT populated in this case
      !! (the caller gets it from that direct call instead). When
      !! seen_reffile is false: builds the fresh template described in
      !! this module's own header comment, returns used_target_grid=
      !! .true. with pixaxes_ref=(1,2)/naxes_ref populated exactly as
      !! load_wcs+extract_sky_mapping would have populated them from a
      !! file, the PLAIN VALUES (center_lon_deg/center_lat_deg/
      !! cdelt1_mag/cdelt2_mag/final_crpix) needed to rebuild the
      !! Mapping/SkyFrame -- deliberately not live AST object
      !! handles, see this module's own header comment for why -- and
      !! n_planes_per_infile (UR-15), piggybacked from the overlap
      !! check this path already runs internally.
      !!
      !! grid_system/grid_projection are the system and projection the
      !! grid was built in, and are what callers must pass to
      !! build_target_skymap/write_target_wcs_cards -- not their own
      !! target_system/target_projection, which may be blank:
      !! target_projection= alone inherits the infiles' own system, and
      !! an omitted target_projection= means sin.
      logical, intent(in) :: seen_reffile
      character(len=*), intent(in) :: target_system, target_projection
      character(len=*), intent(in) :: target_center
      character(len=*), intent(in) :: target_lon, target_lat
      character(len=*), intent(in) :: target_cdelt_lon, target_cdelt_lat
      character(len=*), intent(in) :: target_hpx_nside
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      logical, intent(out) :: used_target_grid
      character(len=*), intent(out) :: grid_system, grid_projection
      double precision, intent(out) :: center_lon_deg, center_lat_deg
      double precision, intent(out) :: cdelt1_mag, cdelt2_mag
      double precision, intent(out) :: final_crpix(2)
      integer, intent(out) :: pixaxes_ref(2)
      integer, intent(out) :: naxes_ref(:)
      integer, intent(out) :: n_planes_per_infile(:)
      integer, intent(out) :: status

      logical :: any_target_key_given
      character(len=4) :: ctype_suffix

      status = 0
      used_target_grid = .false.
      grid_system = ' '
      grid_projection = ' '
      center_lon_deg = 0.0d0; center_lat_deg = 0.0d0
      cdelt1_mag = 0.0d0; cdelt2_mag = 0.0d0
      final_crpix = 0.0d0
      pixaxes_ref = 0
      naxes_ref = 0

      any_target_key_given = len_trim(target_system).ne.0 .or.&
      &len_trim(target_projection).ne.0 .or. len_trim(target_center).ne.0 .or.&
      &len_trim(target_lon).ne.0 .or. len_trim(target_lat).ne.0 .or.&
      &len_trim(target_cdelt_lon).ne.0 .or. len_trim(target_cdelt_lat).ne.0 .or.&
      &len_trim(target_hpx_nside).ne.0

      if (seen_reffile) then
         if (any_target_key_given) then
            write(*,*) 'ERROR: reffile= already fully determines the output'//&
            &' grid -- target_system=/target_projection=/target_center=/'//&
            &'target_lon=/target_lat=/target_cdelt_lon=/target_cdelt_lat=/'//&
            &'target_hpx_nside=',&
            &' are not allowed together with reffile='
            status = -1
         endif
         return
      endif

      ! No reffile from here on.
      if (len_trim(target_system).eq.0 .and. len_trim(target_projection).eq.0) then
         write(*,*) 'ERROR: with no reffile=, at least one of target_system=',&
         &' or target_projection= is required -- nothing requested'
         status = -1
         return
      endif

      ! Validates the value and normalises it to lower case, so the
      ! name stored in grid_projection is the one every later call
      ! maps to a CTYPE suffix.
      call projection_name_to_ctype(target_projection, grid_projection,&
      &ctype_suffix, status)
      if (status.ne.0) return

      if (len_trim(target_system).ne.0) then
         grid_system = target_system
      else
         call infiles_common_system(infiles, n_inputs, grid_system, status)
         if (status.ne.0) return
      endif

      call resolve_target_grid_noref(grid_system, grid_projection, target_center,&
      &target_lon, target_lat, target_cdelt_lon, target_cdelt_lat, target_hpx_nside,&
      &infiles, n_inputs, center_lon_deg, center_lat_deg, cdelt1_mag, cdelt2_mag,&
      &final_crpix, naxes_ref, n_planes_per_infile, status)
      if (status.ne.0) return
      pixaxes_ref(1) = 1
      pixaxes_ref(2) = 2
      used_target_grid = .true.
   end subroutine resolve_target_grid

   subroutine resolve_target_grid_noref(target_system, target_projection,&
   &target_center, target_lon,&
   &target_lat, target_cdelt_lon, target_cdelt_lat, target_hpx_nside, infiles,&
   &n_inputs, center_lon_deg, center_lat_deg, cdelt1_mag_out, cdelt2_mag_out,&
   &final_crpix_out, naxes_ref, n_planes_per_infile, status)
      !! All AST object construction in this subroutine (the trial WCS,
      !! the per-infile footprint check) is local and transient, used
      !! only to compute the plain numbers returned -- never crosses a
      !! thread boundary, so the per-thread-ownership restriction
      !! this module's own header comment describes does not apply
      !! here. It applies to build_target_skymap, below, which every
      !! CALLER of these returned numbers must use instead of sharing
      !! an object built here.
      character(len=*), intent(in) :: target_system, target_center
      character(len=*), intent(in) :: target_projection
      character(len=*), intent(in) :: target_lon, target_lat
      character(len=*), intent(in) :: target_cdelt_lon, target_cdelt_lat
      character(len=*), intent(in) :: target_hpx_nside
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      double precision, intent(out) :: center_lon_deg, center_lat_deg
      double precision, intent(out) :: cdelt1_mag_out, cdelt2_mag_out
      double precision, intent(out) :: final_crpix_out(2)
      integer, intent(out) :: naxes_ref(:)
      integer, intent(out) :: n_planes_per_infile(:)
      integer, intent(out) :: status

      character(len=16) :: ast_system
      character(len=16) :: projection_name
      character(len=4) :: ctype_suffix
      logical :: is_equatorial
      double precision :: lon0_deg, lat0_deg, cdelt1_mag, cdelt2_mag
      integer :: skymap_trial, skyframe_trial, pixaxes_trial(2)
      integer :: annul_status
      double precision :: lbnd(2), ubnd(2)
      double precision :: final_crpix(2)
      double precision :: lattice_d(2)
      integer :: final_naxis1, final_naxis2
      integer :: nside, hpx_naxis, i_hpx, k_hpx, naxis_i_hpx, n_other_i_hpx
      integer :: wcs_i_hpx, skymap_i_hpx, skyframe_i_hpx, pixaxes_i_hpx(2)
      integer :: naxes_i_hpx(max_axes_tg)
      double precision :: hpx_crpix, nside_parsed

      center_lon_deg = 0.0d0; center_lat_deg = 0.0d0
      cdelt1_mag_out = 0.0d0; cdelt2_mag_out = 0.0d0
      final_crpix_out = 0.0d0
      naxes_ref = 0
      status = 0

      call system_name_to_ast(target_system, ast_system, status)
      if (status.ne.0) return
      call projection_name_to_ctype(target_projection, projection_name,&
      &ctype_suffix, status)
      if (status.ne.0) return
      is_equatorial = trim(target_system).ne.'galactic'

      ! UR-23: a HEALPix-aligned HPX grid is fully determined by Nside
      ! alone -- CRVAL/CDELT/CRPIX/NAXIS all follow from it (see
      ! hpx_cell_deg/hpx_naxis_crpix), so none of target_center=/
      ! target_lon=/target_lat=/target_cdelt_lon=/target_cdelt_lat=
      ! has anything left to choose; allowing them would silently
      ! ignore whichever one the caller gave, which is worse than
      ! refusing. This whole branch returns before any of the
      ! generic centroid/CDELT/footprint-union logic below runs.
      if (trim(projection_name).eq.'hpx') then
         if (len_trim(target_center).ne.0 .or. len_trim(target_lon).ne.0 .or.&
         &len_trim(target_lat).ne.0 .or. len_trim(target_cdelt_lon).ne.0 .or.&
         &len_trim(target_cdelt_lat).ne.0) then
            write(*,*) 'ERROR: target_center=/target_lon=/target_lat=/'//&
            &'target_cdelt_lon=/target_cdelt_lat= are not allowed with'//&
            &' target_projection=hpx -- a HEALPix-aligned grid''s centre,'//&
            &' pixel scale and reference pixel all follow from'//&
            &' target_hpx_nside= (or its own default) alone'
            status = -1
            return
         endif
         if (len_trim(target_hpx_nside).ne.0) then
            call parse_plain_number(target_hpx_nside, nside_parsed, status)
            if (status.ne.0) return
            nside = nint(nside_parsed)
            if (.not.nside_is_power_of_two(nside)) then
               write(*,*) 'ERROR: target_hpx_nside= must be a power of 2,'//&
               &' got "', trim(target_hpx_nside), '"'
               status = -1
               return
            endif
         else
            call default_hpx_nside(infiles, n_inputs, nside, status)
            if (status.ne.0) return
         endif

         call hpx_naxis_crpix(nside, hpx_naxis, hpx_crpix)
         center_lon_deg = 0.0d0
         center_lat_deg = 0.0d0
         cdelt1_mag_out = hpx_cell_deg(nside)
         cdelt2_mag_out = cdelt1_mag_out
         final_crpix_out(1) = hpx_crpix
         final_crpix_out(2) = hpx_crpix
         naxes_ref(1) = hpx_naxis
         naxes_ref(2) = hpx_naxis

         n_planes_per_infile = 1
         do i_hpx = 1, n_inputs
            call load_wcs_plain(infiles(i_hpx), wcs_i_hpx, naxes_i_hpx, status)
            if (status.ne.0) return
            call extract_sky(wcs_i_hpx, skymap_i_hpx, skyframe_i_hpx,&
            &pixaxes_i_hpx, status)
            annul_status = 0
            call ast_annul(skymap_i_hpx, annul_status)
            call ast_annul(skyframe_i_hpx, annul_status)
            call ast_annul(wcs_i_hpx, annul_status)
            if (status.ne.0) return
            naxis_i_hpx = 0
            do k_hpx = 1, size(naxes_i_hpx)
               if (naxes_i_hpx(k_hpx).gt.0) naxis_i_hpx = k_hpx
            enddo
            n_other_i_hpx = 1
            do k_hpx = 1, naxis_i_hpx
               if (k_hpx.ne.pixaxes_i_hpx(1) .and. k_hpx.ne.pixaxes_i_hpx(2))&
               &n_other_i_hpx = n_other_i_hpx * naxes_i_hpx(k_hpx)
            enddo
            n_planes_per_infile(i_hpx) = n_other_i_hpx
         enddo

         write(*,'(A,A,A,I0,A,I0,A,I0,A)') 'target_system=', trim(target_system),&
         &' target_projection=hpx: built a HEALPix-aligned ', hpx_naxis, 'x',&
         &hpx_naxis, ' template, Nside=', nside,&
         &' -- no reffile given'
         return
      else if (len_trim(target_hpx_nside).ne.0) then
         write(*,*) 'ERROR: target_hpx_nside= is only allowed with',&
         &' target_projection=hpx'
         status = -1
         return
      endif

      if (len_trim(target_center).eq.0 .or. trim(target_center).eq.'infiles_centroid') then
         if (len_trim(target_lon).ne.0 .or. len_trim(target_lat).ne.0) then
            write(*,*) 'ERROR: target_lon=/target_lat= require',&
            &' target_center=user_defined'
            status = -1
            return
         endif
      else if (trim(target_center).eq.'user_defined') then
         if (len_trim(target_lon).eq.0 .or. len_trim(target_lat).eq.0) then
            write(*,*) 'ERROR: target_center=user_defined requires both',&
            &' target_lon= and target_lat='
            status = -1
            return
         endif
      else if (trim(target_center).eq.'reffile') then
         write(*,*) 'ERROR: target_center=reffile requires reffile= -- got',&
         &' no reffile'
         status = -1
         return
      else
         write(*,*) 'ERROR: target_center must be infiles_centroid or',&
         &' user_defined, got "', trim(target_center), '"'
         status = -1
         return
      endif

      if ((len_trim(target_cdelt_lon).ne.0) .neqv. (len_trim(target_cdelt_lat).ne.0)) then
         write(*,*) 'ERROR: target_cdelt_lon= and target_cdelt_lat= must be given',&
         &' together, or both omitted (then the smallest |CDELT| across the',&
         &' infiles is used, per axis)'
         status = -1
         return
      endif

      ! -- centre --
      if (trim(target_center).eq.'user_defined') then
         call parse_coord_value(target_lon, is_equatorial, lon0_deg, status)
         call parse_coord_value(target_lat, .false., lat0_deg, status)
         if (status.ne.0) return
      else
         block
            integer :: target_skyframe_tmp, wcs0, skymap0, skyframe0, pixaxes0(2)
            integer :: naxes0(max_axes_tg)
            call load_wcs_plain(infiles(1), wcs0, naxes0, status)
            if (status.ne.0) return
            call extract_sky(wcs0, skymap0, skyframe0, pixaxes0, status)
            if (status.ne.0) then
               annul_status = 0
               call ast_annul(wcs0, annul_status)
               return
            endif
            target_skyframe_tmp = ast_copy_ast(skyframe0, status)
            call ast_setc(target_skyframe_tmp, 'System', trim(ast_system), status)
            call compute_infiles_centroid(infiles, n_inputs, target_skyframe_tmp,&
            &lon0_deg, lat0_deg, status)
            annul_status = 0
            call ast_annul(wcs0, annul_status)
            call ast_annul(skymap0, annul_status)
            call ast_annul(skyframe0, annul_status)
            call ast_annul(target_skyframe_tmp, annul_status)
            if (status.ne.0) return
         end block
      endif

      ! -- CDELT --
      ! Explicit target_cdelt_lon=/target_cdelt_lat= win. Otherwise the
      ! default is the SMALLEST |CDELT| across every infile, taken per
      ! axis (UR-16, docs/dev/BRANCH_STATUS_update-reproject.md): the
      ! output grid then samples every input at least as finely as that
      ! input's own grid at the grid centre, so no input is coarsened
      ! there; away from the centre the projection changes the sky per
      ! pixel (sin by 1/cos of the angle from it). Resampling onto an
      ! equal-or-finer grid conserves flux and position (REPROJECT_FLUX_
      ! NOTE.md sections 1-4); the integral errors that note measures
      ! come from resampling much COARSER than an input (section 8).
      ! One rule covers every case: one infile gives that file's own
      ! scale, several equal ones give the shared scale, several
      ! differing ones give the finest.
      if (len_trim(target_cdelt_lon).ne.0) then
         call parse_plain_number(target_cdelt_lon, cdelt1_mag, status)
         call parse_plain_number(target_cdelt_lat, cdelt2_mag, status)
         if (status.ne.0) return
         cdelt1_mag = abs(cdelt1_mag)
         cdelt2_mag = abs(cdelt2_mag)
         if (cdelt1_mag.le.0.0d0 .or. cdelt2_mag.le.0.0d0) then
            write(*,*) 'ERROR: target_cdelt_lon=/target_cdelt_lat= must be nonzero'
            status = -1
            return
         endif
      else
         call smallest_infile_cdelt(infiles, n_inputs, cdelt1_mag, cdelt2_mag, status)
         if (status.ne.0) return
      endif

      ! -- trial WCS, entirely in memory --
      call build_trial_wcs(ast_system, ctype_suffix, lon0_deg, lat0_deg,&
      &cdelt1_mag, cdelt2_mag, skymap_trial, skyframe_trial, pixaxes_trial, status)
      if (status.ne.0) return

      ! -- pixel lattice (UR-22): the trial grid's reference pixel
      ! centre sits on the infiles' centroid, which for an even-sized
      ! infile is a pixel CORNER of that infile, so without this the
      ! output lattice is half a pixel off the input's and the union
      ! rounding below adds a row and a column. Moving the trial
      ! reference pixel by lattice_d puts the first infile's pixel
      ! centres on whole trial pixels, so the output keeps the infile's
      ! own pixel centres. Done only where the two grids differ by a
      ! pure shift; otherwise lattice_d is zero and the grid is as
      ! before. The centroid stays at CRVAL. --
      call infile_lattice_offset(infiles(1), skymap_trial, skyframe_trial,&
      &lattice_d, status)
      if (status.ne.0) return
      if (any(lattice_d.ne.0.0d0)) then
         annul_status = 0
         call ast_annul(skymap_trial, annul_status)
         call ast_annul(skyframe_trial, annul_status)
         call build_trial_wcs(ast_system, ctype_suffix, lon0_deg, lat0_deg,&
         &cdelt1_mag, cdelt2_mag, skymap_trial, skyframe_trial, pixaxes_trial,&
         &status, crpix_trial=1.0d0 + lattice_d)
         if (status.ne.0) return
      endif

      ! -- footprint union across every infile (UR-9/UR-15,
      ! docs/dev/BRANCH_STATUS_update-reproject.md): shared with the
      ! reffile-given union/intersection/reference paths via check_
      ! overlap_and_size_grid, below -- see that subroutine's own
      ! header for the zero-overlap hard-failure reasoning (a band
      ! sharing no sky at all with the rest is not legitimate partial
      ! coverage, it is almost always the wrong file) and the collect-
      ! every-offender-before-reporting behaviour. have_fixed_bound=
      ! .false. here: there is no reffile to seed a bound from, so the
      ! first infile establishes it instead (exempt from the check,
      ! same as before this was shared). --
      lbnd(1) = huge(1.0d0); lbnd(2) = huge(1.0d0)
      ubnd(1) = -huge(1.0d0); ubnd(2) = -huge(1.0d0)
      call check_overlap_and_size_grid(skymap_trial, skyframe_trial, infiles,&
      &n_inputs, 'union', .false., lbnd, ubnd, n_planes_per_infile, status)
      if (status.ne.0) return

      final_naxis1 = nint(ubnd(1) - lbnd(1)) + 1
      final_naxis2 = nint(ubnd(2) - lbnd(2)) + 1
      ! lbnd is in the shifted trial grid's pixels; build_target_skymap
      ! rebuilds the trial grid with CRPIX=(1,1), whose pixel
      ! coordinates are lattice_d smaller, so final_crpix carries it.
      final_crpix(1) = 1.0d0 - (lbnd(1) - 1.0d0) + lattice_d(1)
      final_crpix(2) = 1.0d0 - (lbnd(2) - 1.0d0) + lattice_d(2)

      center_lon_deg = lon0_deg
      center_lat_deg = lat0_deg
      cdelt1_mag_out = cdelt1_mag
      cdelt2_mag_out = cdelt2_mag
      final_crpix_out = final_crpix
      naxes_ref(1) = final_naxis1
      naxes_ref(2) = final_naxis2

      annul_status = 0
      call ast_annul(skymap_trial, annul_status)
      call ast_annul(skyframe_trial, annul_status)

      write(*,'(A,A,A,A,A,I0,A,I0,A,F0.6,A,F0.6,A)') 'target_system=',&
      &trim(target_system), ' target_projection=', trim(projection_name),&
      &': built a ', final_naxis1, 'x', final_naxis2,&
      &' template, centre (', lon0_deg, ',', lat0_deg, ') deg -- no reffile given'
   end subroutine resolve_target_grid_noref

   subroutine check_overlap_and_size_grid(skymap_target, skyframe_target, infiles,&
   &n_inputs, mode, have_fixed_bound, lbnd_out, ubnd_out, n_planes_per_infile, status)
      !! UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): the piece
      !! shared by every footprint-mode variant across both tools --
      !! compute each infile's own footprint bounds against the target
      !! grid's sky Mapping, check it for zero overlap, and (for
      !! union/intersection) fold it into the running grid bound.
      !! Collects every offending infile instead of stopping at the
      !! first (so a systematic problem, e.g. a whole band's Q/U pair,
      !! is told apart from one bad file), reports them all together,
      !! and aborts before the caller rounds/sizes anything if the
      !! list is non-empty.
      !!
      !! mode='intersection'/'union': lbnd_out/ubnd_out is the running
      !! bound -- every infile that is not an offender is folded into
      !! it (intersection shrinks, union grows); an offender is
      !! checked but never folded, so a later infile is still checked
      !! against the bound built from the valid ones before it (the
      !! "daisy-chain" case -- infile B overlapping an already-folded
      !! infile A rather than the original reffile directly -- still
      !! works). mode='reference': lbnd_out/ubnd_out is fixed by the
      !! caller (reffile's own extent) and is never folded -- every
      !! infile is checked against it, full stop.
      !!
      !! have_fixed_bound=.true.: lbnd_out/ubnd_out already holds a
      !! real bound (reffile's own extent) on entry, and every infile
      !! (i=1..n_inputs) is checked against it. have_fixed_bound=
      !! .false.: no reffile exists to seed a bound from -- the caller
      !! passes lbnd_out/ubnd_out pre-set to (+huge,+huge)/(-huge,
      !! -huge) and mode must be 'union'; infile 1 folds unconditionally
      !! (it trivially becomes the bound, since min/max against
      !! +-huge is a no-op) and is exempt from the check (there is
      !! nothing to compare it against yet), infiles 2..n_inputs are
      !! checked and folded as usual.
      !!
      !! n_planes_per_infile(i) is the product of infile(i)'s own non-
      !! sky axis extents (spectral/Stokes/etc.) -- piggybacked here at
      !! no extra cost (naxes_i is already loaded for the bounds
      !! check) so the caller's own disk-space pre-check needs no
      !! second pass over the infiles.
      integer, intent(in) :: skymap_target, skyframe_target
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      character(len=*), intent(in) :: mode
      logical, intent(in) :: have_fixed_bound
      double precision, intent(inout) :: lbnd_out(2), ubnd_out(2)
      integer, intent(out) :: n_planes_per_infile(:)
      integer, intent(out) :: status

      integer :: i, k, naxis_i, n_other_i
      integer :: wcs_i, skymap_i, skyframe_i, pixaxes_i(2)
      integer :: naxes_i(max_axes_tg), map_i2target, annul_status
      double precision :: this_lbnd(2), this_ubnd(2)
      character(len=512) :: offenders(size(infiles))
      integer :: n_offenders
      logical :: check_this_one, is_offender

      status = 0
      n_offenders = 0
      n_planes_per_infile = 1

      do i = 1, n_inputs
         call load_wcs_plain(infiles(i), wcs_i, naxes_i, status)
         if (status.ne.0) return
         call extract_sky(wcs_i, skymap_i, skyframe_i, pixaxes_i, status)
         if (status.ne.0) then
            annul_status = 0
            call ast_annul(wcs_i, annul_status)
            return
         endif

         naxis_i = 0
         do k = 1, size(naxes_i)
            if (naxes_i(k).gt.0) naxis_i = k
         enddo
         n_other_i = 1
         do k = 1, naxis_i
            if (k.ne.pixaxes_i(1) .and. k.ne.pixaxes_i(2))&
            &n_other_i = n_other_i * naxes_i(k)
         enddo
         n_planes_per_infile(i) = n_other_i

         call compose_pix2pix(skymap_i, skyframe_i, skymap_target,&
         &skyframe_target, map_i2target, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to align input file to the target grid: ',&
            &trim(infiles(i))
            annul_status = 0
            call ast_annul(skymap_i, annul_status)
            call ast_annul(skyframe_i, annul_status)
            call ast_annul(wcs_i, annul_status)
            return
         endif

         call footprint_bounds(map_i2target, naxes_i, pixaxes_i, this_lbnd,&
         &this_ubnd, status)
         write(*,'(A,A,A,F0.2,A,F0.2,A,F0.2,A,F0.2,A)') '  ', trim(infiles(i)),&
         &' footprint in target-grid space: [', this_lbnd(1), ',', this_ubnd(1),&
         &'] x [', this_lbnd(2), ',', this_ubnd(2), ']'

         check_this_one = have_fixed_bound .or. (i.gt.1)
         is_offender = .false.
         if (check_this_one) then
            if (this_ubnd(1).lt.lbnd_out(1) .or. this_lbnd(1).gt.ubnd_out(1) .or.&
            &this_ubnd(2).lt.lbnd_out(2) .or. this_lbnd(2).gt.ubnd_out(2)) then
               is_offender = .true.
               n_offenders = n_offenders + 1
               offenders(n_offenders) = infiles(i)
            endif
         endif

         if (.not. is_offender .and. trim(mode).ne.'reference') then
            if (trim(mode).eq.'intersection') then
               lbnd_out = max(lbnd_out, this_lbnd)
               ubnd_out = min(ubnd_out, this_ubnd)
            else
               lbnd_out = min(lbnd_out, this_lbnd)
               ubnd_out = max(ubnd_out, this_ubnd)
            endif
         endif

         annul_status = 0
         call ast_annul(map_i2target, annul_status)
         call ast_annul(skymap_i, annul_status)
         call ast_annul(skyframe_i, annul_status)
         call ast_annul(wcs_i, annul_status)
      enddo

      if (n_offenders.gt.0) then
         write(*,*) 'ERROR: zero sky overlap between the target grid and the',&
         &' following infile(s):'
         do i = 1, n_offenders
            write(*,*) '  ', trim(offenders(i))
         enddo
         write(*,*) 'Quitting now...'
         status = -1
         return
      endif

      ! The bounds come back from ast_mapbox with round-off: a footprint
      ! edge that falls on a whole pixel arrives as e.g.
      ! 16.000000000000835, and a bare ceiling() turns that into 17,
      ! adding a row or column (measured, UR-22). bound_snap_tol is far
      ! above that round-off and far below any edge that is not on a
      ! whole pixel.
      if (trim(mode).eq.'intersection') then
         lbnd_out = ceiling(lbnd_out - bound_snap_tol)
         ubnd_out = floor(ubnd_out + bound_snap_tol)
      else if (trim(mode).eq.'union') then
         lbnd_out = floor(lbnd_out + bound_snap_tol)
         ubnd_out = ceiling(ubnd_out - bound_snap_tol)
      endif
   end subroutine check_overlap_and_size_grid

   subroutine check_disk_space(output_path, n_inputs, n_planes_per_infile,&
   &nx_out, ny_out, min_free_pct_after_write, required_bytes, avail_bytes,&
   &total_bytes, status)
      !! UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): disk-space
      !! pre-check, run before any output file is created. Refuses when
      !! writing every output would leave less than min_free_pct_after_
      !! write percent of the output disk's TOTAL size free. Measured
      !! against total size, not against the space free right now: a
      !! buffer defined as a fraction of current free space shrinks as
      !! the disk fills, so it protects least when the disk is fullest.
      !!
      !! required_bytes sums 4 bytes/pixel * nx_out * ny_out * n_planes_
      !! per_infile(i) across every infile in n_inputs -- not narrowed
      !! to the infiles that survive each caller's own "skip-if-already-
      !! matched" pre-flight, so an over-estimate, never an under-
      !! estimate. Sizes come from `df --output=size,avail <output_
      !! path>`, captured via the same short-lived-dotfile-in-cwd
      !! pattern run_dry_run_check uses (never /tmp). All three byte
      !! counts are returned (avail_bytes=total_bytes=-1 if df itself
      !! failed, in which case the check is skipped with a WARNING) so
      !! the caller can also report them under dry_run=y.
      character(len=*), intent(in) :: output_path
      integer, intent(in) :: n_inputs
      integer, intent(in) :: n_planes_per_infile(:)
      integer, intent(in) :: nx_out, ny_out
      double precision, intent(in) :: min_free_pct_after_write
      integer(kind=8), intent(out) :: required_bytes, avail_bytes, total_bytes
      integer, intent(out) :: status

      character(len=2048) :: shell_cmd
      character(len=600) :: capture_file
      character(len=600) :: target_dir
      integer :: capture_unit, ios_cap, exitstat_cmd, i, slash_pos
      integer(kind=8) :: avail_kb, total_kb
      double precision :: remaining_bytes, min_free_bytes
      double precision, parameter :: gib = 1024.0d0**3

      status = 0
      avail_bytes = -1_8
      total_bytes = -1_8
      required_bytes = 0_8
      do i = 1, n_inputs
         required_bytes = required_bytes + int(4,8) * int(nx_out,8) *&
         &int(ny_out,8) * int(n_planes_per_infile(i),8)
      enddo

      ! Capture file lands next to output_path itself (guaranteed
      ! writable -- that's where the real output is about to go), not
      ! a bare CWD-relative name -- the latter broke whenever the
      ! running process doesn't own the current working directory
      ! (e.g. a non-root container's WORKDIR).
      slash_pos = index(trim(output_path), '/', back=.true.)
      if (slash_pos.gt.0) then
         target_dir = output_path(1:slash_pos)
      else
         target_dir = ''
      endif
      write(capture_file, '(A,A,I0,A)') trim(target_dir), '.diskcheck_', getpid_wrap(), '.tmp'
      write(shell_cmd, '(A,A,A,A)') 'df --output=size,avail ''', trim(output_path),&
      &''' 2>/dev/null | tail -1 > ', trim(capture_file)
      call execute_command_line(trim(shell_cmd), wait=.true., exitstat=exitstat_cmd)

      total_kb = -1_8
      avail_kb = -1_8
      open(newunit=capture_unit, file=trim(capture_file), status='old',&
      &action='read', iostat=ios_cap)
      if (ios_cap.eq.0) then
         read(capture_unit, *, iostat=ios_cap) total_kb, avail_kb
         close(capture_unit, status='delete')
         if (ios_cap.ne.0) then
            total_kb = -1_8
            avail_kb = -1_8
         endif
      else
         close(capture_unit, status='delete', iostat=ios_cap)
      endif

      if (avail_kb.lt.0_8 .or. total_kb.le.0_8) then
         write(*,*) 'WARNING: could not determine the size/free space of the disk holding ',&
         &trim(output_path), ' (df check failed) -- skipping the disk-space check'
         return
      endif

      avail_bytes = avail_kb * 1024_8
      total_bytes = total_kb * 1024_8
      remaining_bytes = dble(avail_bytes) - dble(required_bytes)
      min_free_bytes = (min_free_pct_after_write / 100.0d0) * dble(total_bytes)
      if (remaining_bytes.lt.min_free_bytes) then
         write(*,'(A)') 'ERROR: not enough disk space on the output disk.'
         write(*,'(A,F0.2,A)') '  output needs:            ', dble(required_bytes)/gib, ' GB'
         write(*,'(A,F0.2,A)') '  free now:                ', dble(avail_bytes)/gib, ' GB'
         write(*,'(A,F0.2,A)') '  would remain after write: ', remaining_bytes/gib, ' GB'
         write(*,'(A,F0.2,A,F0.2,A)') '  minimum to keep free:    ', min_free_bytes/gib,&
         &' GB (min_free_pct_after_write=', min_free_pct_after_write, '% of the disk)'
         write(*,'(A)') 'Free up space, reduce the number of infiles/output grid size,'//&
         &' or lower min_free_pct_after_write.'
         status = -1
      endif
   end subroutine check_disk_space

   integer function getpid_wrap() result(pid)
      !! Wraps the getpid() C library call (via iso_c_binding) to make
      !! the disk-space check's capture-file name collision-proof if
      !! multiple runs somehow overlap in the same directory -- not
      !! needed for correctness of a single run, just cheap safety.
      !! Same copy as each tool's own run_dry_run_check carries
      !! (this project's own "adapt, don't share" convention for a
      !! utility this small).
      use, intrinsic :: iso_c_binding, only: c_int
      interface
         function c_getpid() bind(C, name="getpid") result(r)
            import :: c_int
            integer(c_int) :: r
         end function c_getpid
      end interface
      pid = int(c_getpid())
   end function getpid_wrap

   subroutine build_target_skymap(target_system, target_projection,&
   &center_lon_deg, center_lat_deg,&
   &cdelt1_mag, cdelt2_mag, final_crpix, skymap, skyframe, pixel_axes, status)
      !! Rebuilds the exact final pixel->sky Mapping resolve_target_
      !! frame derived, from the plain values alone -- every object
      !! created here is owned by whichever thread calls this
      !! subroutine, so it is safe to call independently, in parallel,
      !! from every OpenMP thread (confirmed necessary, not just
      !! convenient, by direct multi-threaded test -- see this
      !! module's own header comment). Mirrors resolve_target_grid_
      !! noref's own trial-WCS-plus-shift construction exactly; see
      !! there for the full rationale.
      character(len=*), intent(in) :: target_system, target_projection
      double precision, intent(in) :: center_lon_deg, center_lat_deg
      double precision, intent(in) :: cdelt1_mag, cdelt2_mag
      double precision, intent(in) :: final_crpix(2)
      integer, intent(out) :: skymap, skyframe
      integer, intent(out) :: pixel_axes(2)
      integer, intent(inout) :: status

      integer, external :: ast_shiftmap, ast_cmpmap
      character(len=16) :: ast_system, projection_name
      character(len=4) :: ctype_suffix
      integer :: skymap_trial, skyframe_trial, pixaxes_trial(2)
      double precision :: shift(2)
      integer :: shiftmap_h
      double precision :: pc45(2,2)
      double precision, parameter :: deg45 = 0.78539816339744828d0

      skymap = ast__null
      skyframe = ast__null
      pixel_axes = 0
      if (status.ne.0) return

      call system_name_to_ast(target_system, ast_system, status)
      if (status.ne.0) return
      call projection_name_to_ctype(target_projection, projection_name,&
      &ctype_suffix, status)
      if (status.ne.0) return

      if (trim(projection_name).eq.'hpx') then
         ! UR-23: a HEALPix-aligned HPX grid's pixel axes are rotated
         ! 45 deg from lon/lat (verified against astropy_healpix); see
         ! build_trial_wcs's own comment on why this is a plain PC
         ! matrix, not an AST object composition.
         pc45(1,1) = cos(deg45); pc45(1,2) = -sin(deg45)
         pc45(2,1) = sin(deg45); pc45(2,2) = cos(deg45)
         call build_trial_wcs(ast_system, ctype_suffix, center_lon_deg,&
         &center_lat_deg, cdelt1_mag, cdelt2_mag, skymap_trial, skyframe_trial,&
         &pixaxes_trial, status, pc_matrix=pc45)
      else
         call build_trial_wcs(ast_system, ctype_suffix, center_lon_deg, center_lat_deg,&
         &cdelt1_mag, cdelt2_mag, skymap_trial, skyframe_trial, pixaxes_trial, status)
      endif
      if (status.ne.0) return

      shift(1) = 1.0d0 - final_crpix(1)
      shift(2) = 1.0d0 - final_crpix(2)
      shiftmap_h = ast_shiftmap(2, shift, ' ', status)
      skymap = ast_cmpmap(shiftmap_h, skymap_trial, .true., ' ', status)
      skyframe = skyframe_trial
      pixel_axes(1) = 1
      pixel_axes(2) = 2
   end subroutine build_target_skymap

   subroutine check_output_grid_validity(skymap, pixel_offset1, pixel_offset2,&
   &nx, ny, invalid_frac, status)
      !! UR-10 (docs/dev/BRANCH_STATUS_update-reproject.md): samples
      !! the output grid's own pixel->sky inverse mapping at every
      !! output pixel, once, before the per-plane resample pass --
      !! counts how many return AST's own invalid-point sentinel
      !! (ast__bad_tg), for reasons having nothing to do with missing
      !! input data coverage (a facet seam in -HPX, or any other
      !! projection's own undefined region -- confirmed this is
      !! exactly how AST signals it, REPROJECT_TARGET_FRAME_PLAN.md
      !! section 4c's own facet-seam verification). One allocation
      !! sized to the sky plane only (never multiplied by spectral
      !! channels) -- a single pass, done once per run, not per block.
      !!
      !! Callable for EITHER case (reffile-given or
      !! used_target_grid): pixel_offset1/2 is the shift from this
      !! output grid's own 1..nx/1..ny numbering to skymap's own
      !! native pixel numbering (lbnd_out-1 in each caller's own
      !! variable names; always exactly (0,0) for used_target_grid,
      !! since that case's own output grid IS its skymap's native
      !! numbering -- mode is forced to 'reference' there, so there is
      !! never a narrower/wider final grid than the one
      !! build_target_skymap already built).
      integer, intent(in) :: skymap
      double precision, intent(in) :: pixel_offset1, pixel_offset2
      integer, intent(in) :: nx, ny
      double precision, intent(out) :: invalid_frac
      integer, intent(inout) :: status

      double precision, allocatable :: xin(:), yin(:), xout(:), yout(:)
      integer(kind=8) :: npoint8, k
      integer :: npoint, ox, oy, n_invalid

      invalid_frac = 0.0d0
      if (status.ne.0) return

      npoint8 = int(nx, kind=8) * int(ny, kind=8)
      npoint = int(npoint8)
      allocate(xin(npoint), yin(npoint), xout(npoint), yout(npoint))

      k = 0
      do oy = 1, ny
         do ox = 1, nx
            k = k + 1
            xin(k) = dble(ox) + pixel_offset1
            yin(k) = dble(oy) + pixel_offset2
         enddo
      enddo

      call ast_tran2(skymap, npoint, xin, yin, .true., xout, yout, status)

      n_invalid = 0
      do k = 1, npoint8
         if (xout(k).eq.ast__bad_tg .or. yout(k).eq.ast__bad_tg) n_invalid = n_invalid + 1
      enddo
      invalid_frac = dble(n_invalid) / dble(npoint8)

      deallocate(xin, yin, xout, yout)
   end subroutine check_output_grid_validity

   subroutine write_target_wcs_cards(target_system, target_projection, center_lon_deg,&
   &center_lat_deg, cdelt1_mag, cdelt2_mag, final_crpix, dst_axis1,&
   &dst_axis2, out_unit, status)
      !! Called from write_reprojected_file in place of the file-based
      !! copy_axis_keywords/copy_sky_rotation_matrix/copy_wcs_system_
      !! keywords calls, only when resolve_target_grid reported
      !! used_target_grid. Takes the PLAIN VALUES resolve_target_grid
      !! returned, not AST object handles, and calls build_target_
      !! skymap itself -- so this is safe regardless of which thread
      !! OpenMP's own `single` construct happens to run it on (see this
      !! module's own header comment). dst_axis1/dst_axis2 are the
      !! OUTPUT file's own axis numbers for the 2 sky axes (today
      !! always 1,2 -- see reproject_cubes.f90's own comment on why sky
      !! axes are always output positions 1,2). out_unit must already
      !! have had FTPHPR called on it (structural cards only, no WCS
      !! cards yet) -- see this module's own header comment for why
      !! that makes the FTPREC-append below safe.
      character(len=*), intent(in) :: target_system, target_projection
      double precision, intent(in) :: center_lon_deg, center_lat_deg
      double precision, intent(in) :: cdelt1_mag, cdelt2_mag
      double precision, intent(in) :: final_crpix(2)
      integer, intent(in) :: dst_axis1, dst_axis2
      integer, intent(in) :: out_unit
      integer, intent(inout) :: status

      integer, external :: ast_frame, ast_frameset, ast_write
      integer, external :: ast_fitschan
      logical, external :: ast_findfits
      external :: ast_null

      integer :: pix_frame, new_fs, fitschan, nwritten, fitsstat
      integer :: skymap_ref, skyframe_ref, pixaxes_ref(2)
      character(len=80) :: card
      logical :: found

      if (status.ne.0) return

      call build_target_skymap(target_system, target_projection, center_lon_deg,&
      &center_lat_deg, cdelt1_mag, cdelt2_mag, final_crpix, skymap_ref,&
      &skyframe_ref, pixaxes_ref, status)
      if (status.ne.0) return

      pix_frame = ast_frame(2, ' ', status)
      new_fs = ast_frameset(pix_frame, ' ', status)
      call ast_addframe(new_fs, ast__base, skymap_ref, skyframe_ref, status)

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      call ast_setc(fitschan, 'Encoding', 'FITS-WCS', status)
      nwritten = ast_write(fitschan, new_fs, status)
      if (status.ne.0 .or. nwritten.lt.1) then
         write(*,*) 'ERROR: ast_write failed to build the target WCS header'
         status = -1
         return
      endif

      call ast_seti(fitschan, 'Card', 1, status)
      fitsstat = 0
      do
         found = ast_findfits(fitschan, '%f', card, .true., status)
         if (.not. found) exit
         if (is_wcs_card(card)) then
            call renumber_wcs_card(card, dst_axis1, dst_axis2)
            call FTPREC(out_unit, card, fitsstat)
            if (fitsstat.ne.0) then
               write(*,*) 'ERROR: FTPREC failed writing the target WCS header'
               call printerror(fitsstat)
               status = -1
               return
            endif
         endif
      enddo
   end subroutine write_target_wcs_cards

   logical function is_wcs_card(card) result(r)
      character(len=*), intent(in) :: card
      character(len=8) :: kw
      kw = card(1:8)
      r = kw(1:5).eq.'CTYPE' .or. kw(1:5).eq.'CRVAL' .or.&
      &kw(1:5).eq.'CRPIX' .or. kw(1:5).eq.'CDELT' .or.&
      &kw(1:5).eq.'CROTA' .or. kw(1:2).eq.'PC' .or.&
      &kw(1:7).eq.'RADESYS' .or. kw(1:7).eq.'EQUINOX' .or.&
      &kw(1:7).eq.'LONPOLE' .or. kw(1:7).eq.'LATPOLE' .or.&
      &kw(1:5).eq.'WCSAX'
   end function is_wcs_card

   subroutine renumber_wcs_card(card, dst_axis1, dst_axis2)
      !! ast_write always numbers a fresh 2-axis Frame's own axes 1,2 --
      !! renumber to the caller's own dst_axis1/dst_axis2 (today always
      !! 1,2 too, since sky axes are always output positions 1,2, but
      !! kept general and explicit rather than assumed).
      character(len=*), intent(inout) :: card
      integer, intent(in) :: dst_axis1, dst_axis2
      character(len=8) :: kw
      character(len=1) :: axchar
      if (dst_axis1.eq.1 .and. dst_axis2.eq.2) return
      kw = card(1:8)
      if (kw(1:5).eq.'CTYPE' .or. kw(1:5).eq.'CRVAL' .or. kw(1:5).eq.'CRPIX'&
      &.or. kw(1:5).eq.'CDELT' .or. kw(1:5).eq.'CROTA') then
         axchar = kw(6:6)
         if (axchar.eq.'1') then
            write(card(6:6), '(I1)') dst_axis1
         else if (axchar.eq.'2') then
            write(card(6:6), '(I1)') dst_axis2
         endif
      endif
   end subroutine renumber_wcs_card

   subroutine build_trial_wcs(ast_system, ctype_suffix, lon0_deg, lat0_deg,&
   &cdelt1_mag, cdelt2_mag, skymap, skyframe, pixel_axes, status, crpix_trial,&
   &pc_matrix)
      !! Entirely in-memory FitsChan, CRPIX=(1,1), standard sign
      !! convention (longitude axis negative, latitude axis positive).
      !! ctype_suffix is the projection code from
      !! projection_name_to_ctype (TAN/SIN/CAR/AIT/HPX); AST derives
      !! the projection from it on ast_read.
      character(len=*), intent(in) :: ast_system, ctype_suffix
      double precision, intent(in) :: lon0_deg, lat0_deg, cdelt1_mag, cdelt2_mag
      ! Reference pixel of the trial grid, (1,1) when absent. Only
      ! resolve_target_grid_noref moves it, to put the first infile's
      ! pixel centres on whole trial pixels (see infile_lattice_offset).
      double precision, intent(in), optional :: crpix_trial(2)
      ! PC1_1/PC1_2/PC2_1/PC2_2, identity when absent. UR-23's only
      ! caller: a HEALPix-aligned HPX grid needs its pixel axes rotated
      ! 45 deg from lon/lat (verified against astropy_healpix). Written
      ! as ordinary FITS cards, same as every other card here, so
      ! ast_read parses it exactly as it would any PC matrix in a real
      ! file -- no AST object composition needed, and ast_write (used
      ! later to emit the final header) serializes it back out the
      ! same way, transparently to every other caller.
      double precision, intent(in), optional :: pc_matrix(2,2)
      double precision :: crpix_use(2)
      integer, intent(out) :: skymap, skyframe
      integer, intent(out) :: pixel_axes(2)
      integer, intent(inout) :: status

      integer, external :: ast_fitschan, ast_read
      external :: ast_null

      integer :: fitschan, wcs
      character(len=4) :: lon_code, lat_code

      skymap = ast__null; skyframe = ast__null; pixel_axes = 0
      if (status.ne.0) return

      select case (trim(ast_system))
      case ('GALACTIC')
         lon_code = 'GLON'; lat_code = 'GLAT'
      case default
         lon_code = 'RA--'; lat_code = 'DEC-'
      end select

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      call ast_putfits(fitschan, 'SIMPLE  =                    T', .false., status)
      call ast_putfits(fitschan, 'BITPIX  =                  -32', .false., status)
      call ast_putfits(fitschan, 'NAXIS   =                    2', .false., status)
      call ast_putfits(fitschan, 'NAXIS1  =                   10', .false., status)
      call ast_putfits(fitschan, 'NAXIS2  =                   10', .false., status)
      call put_card_str(fitschan, 'CTYPE1', trim(lon_code)//'-'//trim(ctype_suffix), status)
      call put_card_str(fitschan, 'CTYPE2', trim(lat_code)//'-'//trim(ctype_suffix), status)
      ! RADESYS/EQUINOX: without these, AST defaults any equatorial
      ! system to ICRS regardless of ast_system, silently discarding
      ! FK5/FK4 (UR-7/UR-8 reopened 2026-10-24, docs/dev/
      ! BRANCH_STATUS_update-reproject.md). GALACTIC needs neither --
      ! its CTYPE prefix alone disambiguates it from every equatorial
      ! system, unlike RA---/DEC--, which is shared by ICRS/FK5/FK4.
      select case (trim(ast_system))
      case ('FK5')
         call put_card_str(fitschan, 'RADESYS', 'FK5', status)
         call put_card_dbl(fitschan, 'EQUINOX', 2000.0d0, status)
      case ('FK4')
         call put_card_str(fitschan, 'RADESYS', 'FK4', status)
         call put_card_dbl(fitschan, 'EQUINOX', 1950.0d0, status)
      case ('ICRS')
         call put_card_str(fitschan, 'RADESYS', 'ICRS', status)
      end select
      call put_card_dbl(fitschan, 'CRVAL1', lon0_deg, status)
      call put_card_dbl(fitschan, 'CRVAL2', lat0_deg, status)
      crpix_use = 1.0d0
      if (present(crpix_trial)) crpix_use = crpix_trial
      call put_card_dbl(fitschan, 'CRPIX1', crpix_use(1), status)
      call put_card_dbl(fitschan, 'CRPIX2', crpix_use(2), status)
      call put_card_dbl(fitschan, 'CDELT1', -cdelt1_mag, status)
      call put_card_dbl(fitschan, 'CDELT2', cdelt2_mag, status)
      if (present(pc_matrix)) then
         call put_card_dbl(fitschan, 'PC1_1', pc_matrix(1,1), status)
         call put_card_dbl(fitschan, 'PC1_2', pc_matrix(1,2), status)
         call put_card_dbl(fitschan, 'PC2_1', pc_matrix(2,1), status)
         call put_card_dbl(fitschan, 'PC2_2', pc_matrix(2,2), status)
      endif
      call ast_seti(fitschan, 'Card', 1, status)

      wcs = ast_read(fitschan, status)
      if (status.ne.0 .or. wcs.eq.ast__null) then
         write(*,*) 'ERROR: ast_read failed to build the trial WCS'
         status = -1
         return
      endif
      call extract_sky(wcs, skymap, skyframe, pixel_axes, status)
   end subroutine build_trial_wcs

   subroutine put_card_str(fc, key, val, status)
      integer, intent(in) :: fc
      character(len=*), intent(in) :: key, val
      integer, intent(inout) :: status
      character(len=80) :: c
      c = ' '
      c(1:len(key)) = key
      c(9:10) = '= '
      c(11:) = "'"//trim(val)//"'"
      call ast_putfits(fc, c, .false., status)
   end subroutine put_card_str

   subroutine put_card_dbl(fc, key, val, status)
      integer, intent(in) :: fc
      character(len=*), intent(in) :: key
      double precision, intent(in) :: val
      integer, intent(inout) :: status
      character(len=80) :: c
      character(len=20) :: vstr
      c = ' '
      c(1:len(key)) = key
      c(9:10) = '= '
      write(vstr, '(F20.15)') val
      c(11:30) = vstr
      call ast_putfits(fc, c, .false., status)
   end subroutine put_card_dbl

   integer function ast_copy_ast(obj, status) result(r)
      !! Thin wrapper so resolve_target_grid_noref's own block can
      !! call ast_copy without a top-level external declaration
      !! colliding with anything else in this module.
      integer, intent(in) :: obj
      integer, intent(inout) :: status
      integer, external :: ast_copy
      r = ast_copy(obj, status)
   end function ast_copy_ast

   subroutine system_name_to_ast(target_system, ast_system, status)
      !! target_system= vocabulary -- the exact same set/spelling as
      !! casda_fetch.py's own --frame, for one consistent set of frame
      !! names across the toolchain. AST's own System attribute strings
      !! checked against its skyframe.c source, not assumed.
      character(len=*), intent(in) :: target_system
      character(len=*), intent(out) :: ast_system
      integer, intent(inout) :: status

      select case (trim(target_system))
      case ('icrs')
         ast_system = 'ICRS'
      case ('fk5')
         ast_system = 'FK5'
      case ('fk4')
         ast_system = 'FK4'
      case ('galactic')
         ast_system = 'GALACTIC'
      case default
         write(*,*) 'ERROR: target_system must be icrs, fk5, fk4, or galactic,'//&
         &' got "', trim(target_system), '"'
         status = -1
      end select
   end subroutine system_name_to_ast

   subroutine projection_name_to_ctype(target_projection, projection_name,&
   &ctype_suffix, status)
      !! target_projection= vocabulary: a small named set rather than
      !! AST's whole projection catalogue -- the two pointed-field
      !! projections (tan, sin) and three for wide or all-sky fields
      !! (car, ait, hpx). Matched case-insensitively; blank means sin,
      !! the projection the no-reffile grid used before this key existed.
      !! projection_name returns the normalised lower-case name, which
      !! is what callers carry forward and pass back in.
      character(len=*), intent(in) :: target_projection
      character(len=*), intent(out) :: projection_name
      character(len=*), intent(out) :: ctype_suffix
      integer, intent(inout) :: status

      character(len=16) :: name_lc
      integer :: k, code

      projection_name = ' '
      ctype_suffix = ' '
      name_lc = adjustl(target_projection)
      do k = 1, len_trim(name_lc)
         code = iachar(name_lc(k:k))
         if (code.ge.iachar('A') .and. code.le.iachar('Z')) name_lc(k:k) = achar(code + 32)
      enddo
      if (len_trim(name_lc).eq.0) name_lc = 'sin'

      select case (trim(name_lc))
      case ('tan')
         ctype_suffix = 'TAN'
      case ('sin')
         ctype_suffix = 'SIN'
      case ('car')
         ctype_suffix = 'CAR'
      case ('ait')
         ctype_suffix = 'AIT'
      case ('hpx')
         ctype_suffix = 'HPX'
      case default
         write(*,*) 'ERROR: target_projection must be tan, sin, car, ait, or hpx,'//&
         &' got "', trim(target_projection), '"'
         status = -1
         return
      end select
      projection_name = name_lc
   end subroutine projection_name_to_ctype

   subroutine infiles_common_system(infiles, n_inputs, system_name, status)
      !! For target_projection= given without target_system=: the grid
      !! stays in the infiles' own sky system, so every infile must be
      !! in one system that target_system= could also name. Differing
      !! systems are an error rather than a pick, since choosing one
      !! would silently re-express the others; every infile's system is
      !! listed so the user can see which differ.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      character(len=*), intent(out) :: system_name
      integer, intent(inout) :: status

      character(len=ast__szchr_tg), external :: ast_getc
      character(len=16) :: this_name, ast_sys
      character(len=16), allocatable :: names(:)
      integer :: i, wcs, skymap, skyframe, pixaxes(2), annul_status
      integer :: naxes(max_axes_tg)
      logical :: all_same, all_supported

      system_name = ' '
      if (status.ne.0) return
      allocate(names(n_inputs))
      all_supported = .true.

      do i = 1, n_inputs
         call load_wcs_plain(infiles(i), wcs, naxes, status)
         if (status.ne.0) return
         call extract_sky(wcs, skymap, skyframe, pixaxes, status)
         if (status.ne.0) then
            annul_status = 0
            call ast_annul(wcs, annul_status)
            return
         endif
         ast_sys = ast_getc(skyframe, 'System', status)
         select case (trim(ast_sys))
         case ('ICRS')
            this_name = 'icrs'
         case ('FK5')
            this_name = 'fk5'
         case ('FK4')
            this_name = 'fk4'
         case ('GALACTIC')
            this_name = 'galactic'
         case default
            this_name = 'AST:'//trim(ast_sys)
            all_supported = .false.
         end select
         names(i) = this_name
         annul_status = 0
         call ast_annul(wcs, annul_status)
         call ast_annul(skymap, annul_status)
         call ast_annul(skyframe, annul_status)
         if (status.ne.0) return
      enddo

      all_same = all(names(1:n_inputs).eq.names(1))
      if (.not. all_supported .or. .not. all_same) then
         if (.not. all_same) then
            write(*,*) 'ERROR: target_projection= without target_system= keeps',&
            &' the infiles'' own sky system, but the infiles are in different',&
            &' systems -- give target_system= to choose one:'
         else
            write(*,*) 'ERROR: target_projection= without target_system= keeps',&
            &' the infiles'' own sky system, which must be icrs, fk5, fk4 or',&
            &' galactic -- give target_system= to choose one:'
         endif
         do i = 1, n_inputs
            write(*,'(A,A,A,A)') '  ', trim(infiles(i)), ': ', trim(names(i))
         enddo
         status = -1
         return
      endif
      system_name = names(1)
   end subroutine infiles_common_system

   subroutine load_wcs_plain(filename, wcs, naxes, status)
      !! Same pattern as reproject_cubes.f90's own load_wcs (FTGREC+
      !! ast_putfits+ast_seti('Card',1)+ast_read, fitschan annulled
      !! right after a successful ast_read).
      character(len=*), intent(in) :: filename
      integer, intent(out) :: wcs
      integer, intent(inout) :: naxes(:)
      integer, intent(inout) :: status

      integer, external :: ast_fitschan, ast_read
      logical, external :: ast_isaframeset
      external :: ast_null

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
      call FTGISZ(unit, size(naxes), naxes, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGHSP/FTGISZ failed for ', trim(filename)
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
   end subroutine load_wcs_plain

   subroutine extract_sky(wcs, skymap, skyframe, pixel_axes, status)
      !! Verbatim adaptation of reproject_cubes.f90's own
      !! extract_sky_mapping -- see there for the full rationale.
      integer, intent(in) :: wcs
      integer, intent(out) :: skymap, skyframe
      integer, intent(out) :: pixel_axes(2)
      integer, intent(inout) :: status

      integer, external :: ast_geti, ast_getframe, ast_pickaxes
      integer, external :: ast_getmapping, ast_simplify
      logical, external :: ast_isaskyframe

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
   end subroutine extract_sky

   subroutine output_cell_size_arcsec(skymap, skyframe, x_ref, y_ref,&
   &cell1_arcsec, cell2_arcsec, status)
      !! A grid's cell size on the sky, per pixel axis, measured through
      !! its own pixel->sky mapping between (x_ref, y_ref) and its
      !! neighbours one pixel along each axis, so a CD/PC matrix's
      !! scaling and rotation are included. Callers pass the grid's
      !! reference pixel: a projection adds no distortion there, so the
      !! result is the pixel size the header sets, comparable between
      !! grids of any projection. Measured elsewhere it would include
      !! the projection's growth away from its reference point (1/cos
      !! of the angle from it, for SIN), which is geometry, not a
      !! sampling choice.
      integer, intent(in) :: skymap, skyframe
      double precision, intent(in) :: x_ref, y_ref
      double precision, intent(out) :: cell1_arcsec, cell2_arcsec
      integer, intent(inout) :: status
      double precision, external :: ast_distance
      double precision :: xin(3), yin(3), xo(3), yo(3), p0(2), p1(2), p2(2)
      double precision, parameter :: rad2arcsec = 206264.80624709636d0

      cell1_arcsec = 0.0d0; cell2_arcsec = 0.0d0
      if (status.ne.0) return
      xin = (/ x_ref, x_ref + 1.0d0, x_ref /)
      yin = (/ y_ref, y_ref, y_ref + 1.0d0 /)
      call ast_tran2(skymap, 3, xin, yin, .true., xo, yo, status)
      if (status.ne.0) return
      if (any(xo.eq.ast__bad_tg) .or. any(yo.eq.ast__bad_tg)) then
         write(*,*) 'ERROR: the grid''s reference pixel has no valid sky position;',&
         &' cannot measure its cell size'
         status = -1
         return
      endif
      p0 = (/ xo(1), yo(1) /); p1 = (/ xo(2), yo(2) /); p2 = (/ xo(3), yo(3) /)
      cell1_arcsec = ast_distance(skyframe, p0, p1, status) * rad2arcsec
      cell2_arcsec = ast_distance(skyframe, p0, p2, status) * rad2arcsec
   end subroutine output_cell_size_arcsec

   subroutine fits_reference_cell_arcsec(filename, cell1_arcsec, cell2_arcsec, status)
      !! A FITS file's cell size on the sky, per pixel axis, measured at
      !! its reference pixel (CRPIX of its two sky axes) through its own
      !! WCS; see output_cell_size_arcsec for why the reference pixel.
      character(len=*), intent(in) :: filename
      double precision, intent(out) :: cell1_arcsec, cell2_arcsec
      integer, intent(inout) :: status
      integer :: wcs_i, skymap_i, skyframe_i, pixaxes_i(2), naxes_i(max_axes_tg)
      integer :: annul_status, unit_i, blocksize_i, fitsstat_i, k
      double precision :: crpix_i(2)
      character(len=16) :: axstr
      character(len=80) :: comment
      ! cfitsio's KEY_NO_EXIST
      integer, parameter :: key_no_exist = 202

      cell1_arcsec = 0.0d0; cell2_arcsec = 0.0d0
      if (status.ne.0) return
      call load_wcs_plain(filename, wcs_i, naxes_i, status)
      if (status.ne.0) return
      call extract_sky(wcs_i, skymap_i, skyframe_i, pixaxes_i, status)
      if (status.ne.0) then
         annul_status = 0
         call ast_annul(wcs_i, annul_status)
         return
      endif

      fitsstat_i = 0; blocksize_i = 1
      call safe_ftopen(unit_i, trim(filename), 0, blocksize_i, fitsstat_i)
      if (fitsstat_i.ne.0) then
         write(*,*) 'ERROR: failed to reopen ', trim(filename), ' for its CRPIX'
         call printerror(fitsstat_i)
         call free_fits_unit(unit_i)
         status = -1
      else
         do k = 1, 2
            write(axstr, '(I0)') pixaxes_i(k)
            call FTGKYD(unit_i, 'CRPIX'//trim(axstr), crpix_i(k), comment, fitsstat_i)
            ! An absent CRPIX defaults to 0.0 (FITS WCS Paper I), as AST reads it.
            if (fitsstat_i.eq.key_no_exist) then
               crpix_i(k) = 0.0d0
               fitsstat_i = 0
            endif
         enddo
         call safe_ftclos(unit_i, fitsstat_i)
         if (fitsstat_i.ne.0) then
            write(*,*) 'ERROR: failed to read CRPIX from ', trim(filename)
            call printerror(fitsstat_i)
            status = -1
         endif
      endif

      if (status.eq.0) then
         call output_cell_size_arcsec(skymap_i, skyframe_i, crpix_i(1), crpix_i(2),&
         &cell1_arcsec, cell2_arcsec, status)
      endif
      annul_status = 0
      call ast_annul(skymap_i, annul_status)
      call ast_annul(skyframe_i, annul_status)
      call ast_annul(wcs_i, annul_status)
   end subroutine fits_reference_cell_arcsec

   integer function common_beam_sampling(bmin_arcsec, cell_max_arcsec) result(verdict)
      !! 0: adequately sampled; 1: barely Nyquist (warn); 2: undersampled
      !! (refuse). See beam_refuse_cells/beam_warn_cells.
      double precision, intent(in) :: bmin_arcsec, cell_max_arcsec
      if (bmin_arcsec .lt. beam_refuse_cells*cell_max_arcsec) then
         verdict = 2
      else if (bmin_arcsec .lt. beam_warn_cells*cell_max_arcsec) then
         verdict = 1
      else
         verdict = 0
      endif
   end function common_beam_sampling

   subroutine check_beam_sampling(infiles, n_inputs, cell1_arcsec, cell2_arcsec,&
   &tool_and_grid_keys, status)
      !! UR-21: refuses reprojecting an infile onto a grid that is
      !! coarser than its own and cannot sample its beam, before any
      !! output is written, naming every offender at once. The rule
      !! applies only when the output cells are larger than the infile's
      !! own: a grid no coarser than the input's cannot alias anything
      !! the input holds, so there an infile that already undersamples
      !! its own beam is warned about, not refused, and an infile with
      !! no beam is accepted. On a coarser grid an infile with no beam is
      !! refused, since aliasing cannot then be ruled out.
      !!
      !! "No larger" (the no-beam case) compares each grid's own pixel
      !! AREA (cell1_arcsec*cell2_arcsec), not a per-axis linear size.
      !! A per-axis comparison needs the two grids' axes to correspond
      !! direction-for-direction, which does not hold in general (an
      !! arbitrary relative rotation, or different coordinate frames
      !! entirely, e.g. galactic vs ICRS); crossing one grid's worst
      !! axis against the other's best axis (an earlier version of this
      !! check) avoids that assumption but is unsound on its own: it
      !! flags a grid as "coarser" than an exact copy of itself whenever
      !! that grid's own two axes are not equal, which HPX and AIT's
      !! real, bounded, projection-native anisotropy (confirmed,
      !! independent of resolution: exactly 8/(3*pi) for HPX) always
      !! triggers, including on a null reprojection. Area needs no
      !! notion of axis correspondence at all: it is unaffected by
      !! either grid's own internal anisotropy or by any rotation
      !! between them, while still catching real coarsening (a
      !! deliberately 2x/10x coarser grid changes area by that much
      !! squared, dwarfing cell_equal_rtol).
      !! tool_and_grid_keys is the tool name followed by the run's own
      !! grid keys (reffile=/mode= or target_*=), quoted back in the
      !! advice so it can be rerun as printed.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      double precision, intent(in) :: cell1_arcsec, cell2_arcsec
      character(len=*), intent(in) :: tool_and_grid_keys
      integer, intent(inout) :: status

      double precision :: cell_max, bmin_i, bmaj_i, bmaj_widest, in1, in2
      double precision :: area_out, area_in
      logical :: coarser
      double precision :: target_arcsec
      logical :: have_beam
      integer :: i, verdict, n_refused, n_nobeam
      character(len=1024) :: infile_list

      if (status.ne.0) return
      cell_max = max(cell1_arcsec, cell2_arcsec)
      area_out = cell1_arcsec * cell2_arcsec
      n_refused = 0; n_nobeam = 0; bmaj_widest = 0.0d0
      infile_list = ' '

      do i = 1, n_inputs
         if (i.eq.1) then
            infile_list = trim(infiles(i))
         else
            infile_list = trim(infile_list)//','//trim(infiles(i))
         endif
         call infile_beam_extent(infiles(i), have_beam, bmin_i, bmaj_i, status)
         if (status.ne.0) return
         call fits_reference_cell_arcsec(infiles(i), in1, in2, status)
         if (status.ne.0) return
         area_in = in1 * in2
         coarser = area_out .gt. area_in*(1.0d0 + cell_equal_rtol)
         if (have_beam .and. .not. coarser) then
            bmaj_widest = max(bmaj_widest, bmaj_i)
            if (common_beam_sampling(bmin_i, max(in1, in2)).eq.2) then
               write(*,'(A,A,A,F0.3,A,F0.2,A)') 'WARNING: ', trim(infiles(i)),&
               &': its beam BMIN=', bmin_i, ' arcsec spans only ', bmin_i/max(in1, in2),&
               &' of its OWN cells, so it is already undersampled; the output grid'//&
               &' is no coarser, so reprojection adds no aliasing, but interpolated'//&
               &' values of undersampled structure are inaccurate.'
            endif
         else if (have_beam) then
            bmaj_widest = max(bmaj_widest, bmaj_i)
            verdict = common_beam_sampling(bmin_i, cell_max)
            if (verdict.eq.2) then
               n_refused = n_refused + 1
               write(*,'(A,A,A,F0.3,A,F0.2,A)') 'ERROR: ', trim(infiles(i)),&
               &': narrowest beam BMIN=', bmin_i, ' arcsec spans only ',&
               &bmin_i/cell_max, ' output cells; at least 2.0 are needed.'
            else if (verdict.eq.1) then
               write(*,'(A,A,A,F0.3,A,F0.2,A)') 'WARNING: ', trim(infiles(i)),&
               &': narrowest beam BMIN=', bmin_i, ' arcsec spans ', bmin_i/cell_max,&
               &' output cells -- barely Nyquist (2.5 or more is recommended).'
            endif
         else
            ! no beam: allowed only if the output is no coarser than this infile
            if (coarser) then
               n_nobeam = n_nobeam + 1
               write(*,'(A,A,A,A,F0.4,A,F0.4,A)') 'ERROR: ', trim(infiles(i)),&
               &' has no beam information (no BMAJ/BMIN, no BEAMS table) and the',&
               &' output cell area (effective size ', sqrt(area_out),&
               &' arcsec) is larger than its own (effective size ', sqrt(area_in),&
               &' arcsec), so aliasing cannot be ruled out.'
            else
               write(*,'(A,A,A)') 'INFO: ', trim(infiles(i)), ' has no beam'//&
               &' information; not checked, since the output grid is no coarser'//&
               &' than its own.'
            endif
         endif
      enddo

      if (n_refused.eq.0 .and. n_nobeam.eq.0) return

      write(*,'(A,F0.4,A,F0.4,A,F0.4,A)') '  Output cells: ', cell1_arcsec, ' x ',&
      &cell2_arcsec, ' arcsec. A beam must span at least 2.0 cells of the larger side'//&
      &' (BMIN >= ', beam_refuse_cells*cell_max, ' arcsec).'
      if (n_refused.gt.0) then
         ! The suggested beam clears the barely-Nyquist warning too, and is
         ! at least the widest input beam so every channel can be
         ! convolved to it.
         target_arcsec = max(beam_warn_cells*cell_max, bmaj_widest)
         target_arcsec = ceiling(target_arcsec*100.0d0)/100.0d0
         write(*,'(A)') '  Reprojecting onto this grid would alias structure it cannot'//&
         &' sample. Convolve to a wider beam first -- either in one step:'
         write(*,'(A,A,F0.2,A,F0.2,A,A,A,A)') '    match_cubes stages=both order=convolve_reproject',&
         &' target_bmaj=', target_arcsec, ' target_bmin=', target_arcsec,&
         &' target_bpa=0 infiles=', trim(infile_list), ' ', trim(tool_and_grid_keys)
         write(*,'(A)') '  or in two:'
         write(*,'(A,A,A,F0.2,A,F0.2,A)') '    convolve_cubes infiles=', trim(infile_list),&
         &' target_bmaj=', target_arcsec, ' target_bmin=', target_arcsec, ' target_bpa=0'
         write(*,'(A)') '    then rerun this reprojection on the *_CONV.FITS outputs.'
         write(*,'(A)') '  Both read each infile''s per-channel beams from its BEAMS table;'//&
         &' an infile without one also needs beamfiles=<beam log> (see --help).'
         write(*,'(A,F0.2,A)') '  (', target_arcsec, ' arcsec = 2.5 output cells, or the'//&
         &' widest input beam if larger; any beam of at least 2.0 cells is accepted.)'
      endif
      if (n_nobeam.gt.0) then
         write(*,'(A)') '  For an infile with no beam: add BMAJ/BMIN/BPA to its header,'//&
         &' convolve it to a known beam first, or use output cells no larger than'//&
         &' its own (target_cdelt_lon=/target_cdelt_lat=, or a reffile= with such cells).'
      endif
      write(*,'(A)') 'Quitting now...'
      status = -1
   end subroutine check_beam_sampling

   subroutine infile_lattice_offset(infile, skymap_trial, skyframe_trial,&
   &lattice_d, status)
      !! The shift, per axis and within [-0.5, 0.5], that moves the
      !! trial grid's reference pixel so infile's pixel centres land on
      !! whole trial pixels. Applied only when the infile->trial
      !! mapping is a pure shift (one-pixel steps map to one-pixel
      !! steps, see lattice_shift_tol): only then can every pixel centre
      !! coincide. A rotated, rescaled or differently projected grid has
      !! no lattice to keep, and shifting its anchor would only move its
      !! reference pixel off a whole pixel. Measured at a central pixel
      !! of infile, so the point lies inside the projection's valid
      !! region. Zero if the mapping is not a pure shift or the point
      !! has no valid mapping.
      character(len=*), intent(in) :: infile
      integer, intent(in) :: skymap_trial, skyframe_trial
      double precision, intent(out) :: lattice_d(2)
      integer, intent(inout) :: status

      integer :: wcs_i, skymap_i, skyframe_i, pixaxes_i(2), map_i2t
      integer :: naxes_i(max_axes_tg), annul_status
      double precision :: xin(3), yin(3), xout(3), yout(3)

      lattice_d = 0.0d0
      if (status.ne.0) return
      call load_wcs_plain(infile, wcs_i, naxes_i, status)
      if (status.ne.0) return
      call extract_sky(wcs_i, skymap_i, skyframe_i, pixaxes_i, status)
      if (status.ne.0) then
         annul_status = 0
         call ast_annul(wcs_i, annul_status)
         return
      endif
      call compose_pix2pix(skymap_i, skyframe_i, skymap_trial, skyframe_trial,&
      &map_i2t, status)
      if (status.eq.0) then
         xin(1) = dble((naxes_i(pixaxes_i(1)) + 1) / 2)
         yin(1) = dble((naxes_i(pixaxes_i(2)) + 1) / 2)
         xin(2) = xin(1) + 1.0d0; yin(2) = yin(1)
         xin(3) = xin(1);         yin(3) = yin(1) + 1.0d0
         call ast_tran2(map_i2t, 3, xin, yin, .true., xout, yout, status)
         if (status.eq.0 .and. all(xout.ne.ast__bad_tg) .and. all(yout.ne.ast__bad_tg)) then
            if (abs(xout(2) - xout(1) - 1.0d0).le.lattice_shift_tol .and.&
            &abs(yout(2) - yout(1)).le.lattice_shift_tol .and.&
            &abs(xout(3) - xout(1)).le.lattice_shift_tol .and.&
            &abs(yout(3) - yout(1) - 1.0d0).le.lattice_shift_tol) then
               lattice_d(1) = dble(nint(xout(1))) - xout(1)
               lattice_d(2) = dble(nint(yout(1))) - yout(1)
            endif
         endif
         annul_status = 0
         call ast_annul(map_i2t, annul_status)
      endif
      annul_status = 0
      call ast_annul(skymap_i, annul_status)
      call ast_annul(skyframe_i, annul_status)
      call ast_annul(wcs_i, annul_status)
   end subroutine infile_lattice_offset

   subroutine compose_pix2pix(skymap_from, skyframe_from, skymap_to,&
   &skyframe_to, map_out, status)
      !! Verbatim adaptation of reproject_cubes.f90's own compose_
      !! pix2pix -- see there for the full rationale. Not shared
      !! directly (that one is private to that program's own `contains`
      !! block, and match_cubes.f90 carries its own separate copy too
      !! -- same "adapt, don't share" convention for this generic
      !! geometric utility).
      integer, intent(in) :: skymap_from, skyframe_from
      integer, intent(in) :: skymap_to, skyframe_to
      integer, intent(out) :: map_out
      integer, intent(inout) :: status

      integer, external :: ast_convert, ast_cmpmap
      integer :: sky2sky

      map_out = ast__null
      if (status.ne.0) return

      sky2sky = ast_convert(skyframe_from, skyframe_to, ' ', status)
      if (status.ne.0 .or. sky2sky.eq.ast__null) then
         status = -1
         return
      endif

      call ast_invert(skymap_to, status)
      map_out = ast_cmpmap(skymap_from, sky2sky, .true., ' ', status)
      map_out = ast_cmpmap(map_out, skymap_to, .true., ' ', status)
      call ast_invert(skymap_to, status)
      if (status.ne.0 .or. map_out.eq.ast__null) then
         status = -1
      endif

      call ast_annul(sky2sky, status)
   end subroutine compose_pix2pix

   subroutine footprint_bounds(map_from_to, naxes_from, pixaxes_from,&
   &lbnd, ubnd, status)
      !! Verbatim adaptation of reproject_cubes.f90's own
      !! footprint_bounds -- see there for the full rationale.
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

   subroutine compute_infiles_centroid(infiles, n_inputs, target_skyframe,&
   &new_lon_deg, new_lat_deg, status)
      !! target_center=infiles_centroid (also the default when
      !! target_center is omitted, with no reffile): average, as 3-D
      !! unit vectors, each infile's own pixel-grid-centre sky position
      !! converted into the target system -- not naive lon/lat
      !! averaging, which breaks near the poles or across the 0/360
      !! wrap. ast_tran2 confirmed empirically to operate in RADIANS
      !! for a SkyFrame-derived Mapping (scratch/build_check/
      !! test_tran2_units.f90, cross-checked against astropy to 14
      !! significant digits).
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      integer, intent(in) :: target_skyframe
      double precision, intent(out) :: new_lon_deg, new_lat_deg
      integer, intent(inout) :: status

      integer, external :: ast_convert

      integer :: i, wcs_i, skymap_i, skyframe_i, pixaxes_i(2)
      integer :: naxes_i(max_axes_tg), sky2sky_i, annul_status
      double precision :: xin(1), yin(1), xout(1), yout(1)
      double precision :: lon_t(1), lat_t(1)
      double precision :: sumx, sumy, sumz, lon_r, lat_r, norm
      double precision, parameter :: pi = 3.14159265358979323846d0

      new_lon_deg = 0.0d0
      new_lat_deg = 0.0d0
      if (status.ne.0) return

      sumx = 0.0d0; sumy = 0.0d0; sumz = 0.0d0

      do i = 1, n_inputs
         call load_wcs_plain(infiles(i), wcs_i, naxes_i, status)
         if (status.ne.0) return

         call extract_sky(wcs_i, skymap_i, skyframe_i, pixaxes_i, status)
         if (status.ne.0) then
            annul_status = 0
            call ast_annul(wcs_i, annul_status)
            return
         endif

         xin(1) = (dble(naxes_i(pixaxes_i(1))) + 1.0d0) / 2.0d0
         yin(1) = (dble(naxes_i(pixaxes_i(2))) + 1.0d0) / 2.0d0
         call ast_tran2(skymap_i, 1, xin, yin, .true., xout, yout, status)

         sky2sky_i = ast_convert(skyframe_i, target_skyframe, ' ', status)
         if (status.ne.0 .or. sky2sky_i.eq.ast__null) then
            write(*,*) 'ERROR: no conversion path from ', trim(infiles(i)),&
            &' own sky system to the target system (target_center=',&
            &'infiles_centroid)'
            status = -1
            annul_status = 0
            call ast_annul(wcs_i, annul_status)
            return
         endif
         call ast_tran2(sky2sky_i, 1, xout, yout, .true., lon_t, lat_t, status)

         sumx = sumx + cos(lat_t(1)) * cos(lon_t(1))
         sumy = sumy + cos(lat_t(1)) * sin(lon_t(1))
         sumz = sumz + sin(lat_t(1))

         annul_status = 0
         call ast_annul(wcs_i, annul_status)
         call ast_annul(skymap_i, annul_status)
         call ast_annul(skyframe_i, annul_status)
         call ast_annul(sky2sky_i, annul_status)
      enddo

      norm = sqrt(sumx**2 + sumy**2 + sumz**2)
      if (norm.le.0.0d0) then
         write(*,*) 'ERROR: target_center=infiles_centroid: the infiles'' own',&
         &' centre points average to a degenerate (zero-length) vector'
         status = -1
         return
      endif

      lat_r = atan2(sumz / norm, sqrt((sumx / norm)**2 + (sumy / norm)**2))
      lon_r = atan2(sumy, sumx)
      if (lon_r.lt.0.0d0) lon_r = lon_r + 2.0d0 * pi

      new_lon_deg = lon_r * 180.0d0 / pi
      new_lat_deg = lat_r * 180.0d0 / pi
   end subroutine compute_infiles_centroid

   subroutine smallest_infile_cdelt(infiles, n_inputs, cdelt1_mag, cdelt2_mag, status,&
   &quiet)
      !! UR-16 (docs/dev/BRANCH_STATUS_update-reproject.md): default
      !! output pixel scale for the no-reffile path -- the smallest
      !! |CDELT| across every infile, taken per axis (longitude and
      !! latitude independently; pixels are square in practice, but
      !! nothing here assumes it). See resolve_target_grid_noref's own
      !! CDELT comment for why smallest rather than largest or a mean.
      !!
      !! Each infile's own sky axes are located first (extract_sky), so
      !! the right CDELTn is read even when the sky axes are not axes
      !! 1,2 -- the same reason the rest of this module never assumes
      !! an axis order. Reports when the infiles disagree, since the
      !! output grid is then finer than some of its own inputs.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      double precision, intent(out) :: cdelt1_mag, cdelt2_mag
      integer, intent(inout) :: status
      ! UR-23's default_hpx_nside reuses this same per-infile scan to
      ! find the finest input pixel, where "Output pixel scale: ..."
      ! would misreport an HPX run's actual output scale (Nside-
      ! derived, not this value) -- quiet=.true. suppresses it there.
      logical, intent(in), optional :: quiet

      integer :: i, wcs_i, skymap_i, skyframe_i, pixaxes_i(2)
      integer :: naxes_i(max_axes_tg), annul_status
      integer :: unit_i, blocksize_i, fitsstat_i
      character(len=8) :: axstr
      character(len=68) :: comment
      double precision :: cdelt_raw1, cdelt_raw2, c1, c2
      double precision :: max1_mag, max2_mag
      logical :: differ

      cdelt1_mag = 0.0d0
      cdelt2_mag = 0.0d0
      if (status.ne.0) return

      max1_mag = 0.0d0
      max2_mag = 0.0d0

      do i = 1, n_inputs
         call load_wcs_plain(infiles(i), wcs_i, naxes_i, status)
         if (status.ne.0) return
         call extract_sky(wcs_i, skymap_i, skyframe_i, pixaxes_i, status)
         annul_status = 0
         call ast_annul(skymap_i, annul_status)
         call ast_annul(skyframe_i, annul_status)
         call ast_annul(wcs_i, annul_status)
         if (status.ne.0) return

         fitsstat_i = 0; blocksize_i = 1
         call safe_ftopen(unit_i, trim(infiles(i)), 0, blocksize_i, fitsstat_i)
         if (fitsstat_i.ne.0) then
            write(*,*) 'ERROR: failed to reopen ', trim(infiles(i)),&
            &' for its own CDELT'
            call printerror(fitsstat_i)
            status = -1
            call free_fits_unit(unit_i)
            return
         endif
         write(axstr, '(I0)') pixaxes_i(1)
         call FTGKYD(unit_i, 'CDELT'//trim(axstr), cdelt_raw1, comment, fitsstat_i)
         write(axstr, '(I0)') pixaxes_i(2)
         call FTGKYD(unit_i, 'CDELT'//trim(axstr), cdelt_raw2, comment, fitsstat_i)
         call safe_ftclos(unit_i, fitsstat_i)
         if (fitsstat_i.ne.0) then
            write(*,*) 'ERROR: failed to read CDELT from ', trim(infiles(i))
            call printerror(fitsstat_i)
            status = -1
            return
         endif

         c1 = abs(cdelt_raw1)
         c2 = abs(cdelt_raw2)
         if (c1.le.0.0d0 .or. c2.le.0.0d0) then
            write(*,*) 'ERROR: ', trim(infiles(i)), ' has a zero CDELT on a sky axis'
            status = -1
            return
         endif
         if (i.eq.1) then
            cdelt1_mag = c1; max1_mag = c1
            cdelt2_mag = c2; max2_mag = c2
         else
            cdelt1_mag = min(cdelt1_mag, c1); max1_mag = max(max1_mag, c1)
            cdelt2_mag = min(cdelt2_mag, c2); max2_mag = max(max2_mag, c2)
         endif
      enddo

      ! Relative tolerance, not equality: these headers store CDELT as
      ! 32-bit values (a 0.001 deg/pixel grid reads back as
      ! -1.00000005E-03), so bit-equal scales differ in the last digits.
      differ = (max1_mag - cdelt1_mag).gt.1.0d-6*max(max1_mag, 1.0d-30) .or.&
      &(max2_mag - cdelt2_mag).gt.1.0d-6*max(max2_mag, 1.0d-30)

      if (present(quiet)) then
         if (quiet) return
      endif
      write(*,'(A,ES12.5,A,ES12.5,A,I0,A)') 'Output pixel scale: ', cdelt1_mag,&
      &' x ', cdelt2_mag, ' deg/pixel (smallest |CDELT| across ', n_inputs,&
      &' infile(s))'
      if (differ) then
         write(*,'(A,ES12.5,A,ES12.5,A)') '  the infiles differ in pixel scale'//&
         &' (largest is ', max1_mag, ' x ', max2_mag, ' deg/pixel) -- the'//&
         &' smallest is used so no input is coarsened at the grid centre; give'
         write(*,'(A)') '  target_cdelt_lon=/target_cdelt_lat= to choose a'//&
         &' different scale.'
      endif
   end subroutine smallest_infile_cdelt

   logical function nside_is_power_of_two(nside) result(ok)
      !! UR-23: target_hpx_nside= must be a power of 2 (HEALPix's own
      !! requirement -- each doubling quarters every pixel by
      !! quadrisecting it, which only tiles evenly starting from the
      !! base 12-pixel grid when Nside itself is a power of 2).
      integer, intent(in) :: nside
      ok = nside.ge.1 .and. iand(nside, nside - 1).eq.0
   end function nside_is_power_of_two

   double precision function hpx_cell_deg(nside) result(cell_deg)
      !! UR-23: a HEALPix-aligned HPX grid's pixel size, exact at any
      !! Nside (verified against astropy_healpix to double-precision
      !! round-off, Nside 1-1024 inclusive, equatorial and polar
      !! facets): |CDELT| = 90/(Nside*sqrt(2)) degrees, equal magnitude
      !! on both axes.
      integer, intent(in) :: nside
      cell_deg = 90.0d0 / (dble(nside) * sqrt(2.0d0))
   end function hpx_cell_deg

   subroutine hpx_naxis_crpix(nside, naxis, crpix)
      !! UR-23: the minimal bounding pixel grid that covers every
      !! HEALPix facet without truncation or duplication, and the
      !! reference pixel that puts a FITS integer pixel centre exactly
      !! on every true HEALPix pixel centre. Both verified empirically
      !! against astropy_healpix (round-trip error at floating-point
      !! round-off, Nside 1-1024): Nside=1 is a genuinely different,
      !! smaller case (CRPIX on a whole pixel); every Nside>=2 follows
      !! NAXIS=5*Nside+2 with CRPIX on a half pixel -- the "wrong"
      !! parity for a given Nside is not merely suboptimal, it is wrong
      !! by tens of degrees, so there is no intermediate rule to
      !! approximate.
      integer, intent(in) :: nside
      integer, intent(out) :: naxis
      double precision, intent(out) :: crpix
      if (nside.eq.1) then
         naxis = 6
         crpix = dble(naxis) / 2.0d0 + 1.0d0
      else
         naxis = 5 * nside + 2
         crpix = dble(naxis) / 2.0d0 + 0.5d0
      endif
   end subroutine hpx_naxis_crpix

   subroutine default_hpx_nside(infiles, n_inputs, nside, status)
      !! UR-23: Nside closest to the finest input pixel (by log2
      !! rounding to the nearest power of 2, reusing smallest_infile_
      !! cdelt's own per-infile scan -- the smaller of its two
      !! returned axis values, since an HPX cell is isotropic and the
      !! comparison needs one scalar). If that choice would coarsen
      !! some infile, it is accepted only when EVERY infile has a beam
      !! spanning >= 2.5 of the HPX cell (common_beam_sampling's own
      !! "adequately sampled" verdict, reusing UR-21's beam_warn_cells
      !! threshold rather than inventing a second one) -- never
      !! coarser at all if any infile has no beam. Otherwise Nside is
      !! doubled (halving the cell) until a safe choice is found, or
      !! until Nside reaches the verified upper bound of 1024.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      integer, intent(out) :: nside
      integer, intent(inout) :: status

      integer, parameter :: nside_min = 1, nside_max = 1024
      double precision :: cdelt1_mag, cdelt2_mag, finest_deg, cell_deg
      double precision :: bmin_i, bmaj_i
      logical :: have_beam, safe
      integer :: i, candidate

      nside = nside_min
      if (status.ne.0) return

      call smallest_infile_cdelt(infiles, n_inputs, cdelt1_mag, cdelt2_mag, status,&
      &quiet=.true.)
      if (status.ne.0) return
      finest_deg = min(cdelt1_mag, cdelt2_mag)

      candidate = 2**nint(log(90.0d0 / (sqrt(2.0d0) * finest_deg)) / log(2.0d0))
      candidate = max(nside_min, min(nside_max, candidate))

      do
         cell_deg = hpx_cell_deg(candidate)
         if (cell_deg.le.finest_deg*(1.0d0 + cell_equal_rtol)) exit
         safe = .true.
         do i = 1, n_inputs
            call infile_beam_extent(infiles(i), have_beam, bmin_i, bmaj_i, status)
            if (status.ne.0) return
            if (.not.have_beam) then
               safe = .false.
               exit
            endif
            if (common_beam_sampling(bmin_i, cell_deg*3600.0d0).ne.0) then
               safe = .false.
               exit
            endif
         enddo
         if (safe .or. candidate.ge.nside_max) exit
         candidate = candidate * 2
      enddo
      nside = candidate
   end subroutine default_hpx_nside

   subroutine parse_plain_number(value_str, val, status)
      !! target_cdelt_lon=/target_cdelt_lat=: a plain decimal number
      !! (degrees per pixel), no sexagesimal/notation handling -- a
      !! pixel scale is never written in hours or arcmin:arcsec
      !! colon-sexagesimal notation in practice, unlike a sky position.
      character(len=*), intent(in) :: value_str
      double precision, intent(out) :: val
      integer, intent(inout) :: status
      integer :: ios
      val = 0.0d0
      if (status.ne.0) return
      read(value_str, *, iostat=ios) val
      if (ios.ne.0) then
         write(*,*) 'ERROR: could not parse target_cdelt value "',&
         &trim(value_str), '"'
         status = -1
      endif
   end subroutine parse_plain_number

   subroutine parse_coord_value(value_str, hour_fallback, deg_value, status)
      !! Fortran equivalent of scripts/casda_fetch.py's own
      !! _parse_coord_axis: a letter-suffixed sexagesimal string
      !! ('12h34m56.7s', '-45d30m00s', '180.5d') is self-describing and
      !! parsed directly; a bare colon-sexagesimal string ('12:34:56.7')
      !! is ambiguous and resolved by hour_fallback (true for the
      !! longitude-like axis of an equatorial target_system -- the
      !! standard RA convention; false everywhere else, including
      !! every axis of a frame with no hour-angle convention of its own,
      !! e.g. Galactic); a bare decimal number is always degrees, on
      !! every axis, never hour angle -- same rules, same order of
      !! preference, as the Python original.
      character(len=*), intent(in) :: value_str
      logical, intent(in) :: hour_fallback
      double precision, intent(out) :: deg_value
      integer, intent(inout) :: status

      logical :: ok
      character(len=len(value_str)) :: work
      double precision :: sign_val, d1, d2, d3
      integer :: p1, p2, ios

      deg_value = 0.0d0
      if (status.ne.0) return

      call parse_sexagesimal_lettered(value_str, deg_value, ok)
      if (ok) return

      work = adjustl(trim(value_str))
      if (len_trim(work).eq.0) then
         write(*,*) 'ERROR: empty coordinate value'
         status = -1
         return
      endif

      sign_val = 1.0d0
      if (work(1:1).eq.'-') then
         sign_val = -1.0d0
         work = work(2:)
      else if (work(1:1).eq.'+') then
         work = work(2:)
      endif

      p1 = index(work, ':')
      if (p1.gt.0) then
         read(work(1:p1-1), *, iostat=ios) d1
         if (ios.ne.0) then
            write(*,*) 'ERROR: could not parse coordinate value "',&
            &trim(value_str), '"'
            status = -1
            return
         endif
         p2 = index(work(p1+1:), ':')
         d3 = 0.0d0
         if (p2.gt.0) then
            p2 = p1 + p2
            read(work(p1+1:p2-1), *, iostat=ios) d2
            if (ios.eq.0) read(work(p2+1:), *, iostat=ios) d3
         else
            read(work(p1+1:), *, iostat=ios) d2
         endif
         if (ios.ne.0) then
            write(*,*) 'ERROR: could not parse coordinate value "',&
            &trim(value_str), '"'
            status = -1
            return
         endif
         deg_value = sign_val * (d1 + d2 / 60.0d0 + d3 / 3600.0d0)
         if (hour_fallback) deg_value = deg_value * 15.0d0
         return
      endif

      read(work, *, iostat=ios) d1
      if (ios.ne.0) then
         write(*,*) 'ERROR: could not parse coordinate value "',&
         &trim(value_str), '"'
         status = -1
         return
      endif
      deg_value = sign_val * d1
   end subroutine parse_coord_value

   subroutine parse_sexagesimal_lettered(str, deg_value, ok)
      !! '12h34m56.7s' / '-45d30m00s' / '180.5d' / '12h' / '12h34m' --
      !! a leading optional sign, a primary number followed by h|d
      !! (hour- or degree-based), optionally a minutes number followed
      !! by m, optionally a seconds number followed by s. ok=.false.
      !! (not an error) when no h/d letter is found at all -- the
      !! caller then falls back to colon-sexagesimal/bare-decimal
      !! parsing, exactly mirroring Angle(value_str) raising UnitsError
      !! only for a string with no self-describing unit.
      character(len=*), intent(in) :: str
      double precision, intent(out) :: deg_value
      logical, intent(out) :: ok

      character(len=len(str)) :: work
      double precision :: sign_val, primary_val, min_val, sec_val
      integer :: pos_prim, pos_min, pos_sec, ios
      logical :: is_hour

      ok = .false.
      deg_value = 0.0d0

      work = adjustl(trim(str))
      if (len_trim(work).eq.0) return

      sign_val = 1.0d0
      if (work(1:1).eq.'-') then
         sign_val = -1.0d0
         work = work(2:)
      else if (work(1:1).eq.'+') then
         work = work(2:)
      endif

      pos_prim = scan(work, 'hHdD')
      if (pos_prim.eq.0) return
      is_hour = (work(pos_prim:pos_prim).eq.'h' .or. work(pos_prim:pos_prim).eq.'H')

      read(work(1:pos_prim-1), *, iostat=ios) primary_val
      if (ios.ne.0) return

      min_val = 0.0d0
      sec_val = 0.0d0
      pos_min = scan(work(pos_prim+1:), 'mM')
      if (pos_min.gt.0) then
         pos_min = pos_prim + pos_min
         read(work(pos_prim+1:pos_min-1), *, iostat=ios) min_val
         if (ios.ne.0) return
         pos_sec = scan(work(pos_min+1:), 'sS')
         if (pos_sec.gt.0) then
            pos_sec = pos_min + pos_sec
            read(work(pos_min+1:pos_sec-1), *, iostat=ios) sec_val
            if (ios.ne.0) return
         endif
      endif

      deg_value = sign_val * (primary_val + min_val / 60.0d0 + sec_val / 3600.0d0)
      if (is_hour) deg_value = deg_value * 15.0d0
      ok = .true.
   end subroutine parse_sexagesimal_lettered

end module target_grid_mod
