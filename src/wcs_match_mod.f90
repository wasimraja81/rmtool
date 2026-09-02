module wcs_match_mod
   !! Sky-axis WCS equivalence check shared by match_cubes (skip-if-
   !! already-matched), reproject_cubes (same), and rm_synthesis
   !! (multi-band cross-band validation). One implementation instead of
   !! three (T36, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) -- match_cubes
   !! and reproject_cubes each previously carried their own verbatim
   !! copy of this function.
   implicit none
   private
   public :: sky_wcs_matches_target

contains

   logical function sky_wcs_matches_target(cand_unit, cand_axis1, cand_axis2,&
   &nx_cand, ny_cand, ref_unit, ref_axis1, ref_axis2, crpix_shift1,&
   &crpix_shift2, nx_out, ny_out) result(matches)
      !! Tight, "already-processed-identically" comparison of a candidate
      !! file's own sky-axis WCS (cand_unit, its own axis numbering)
      !! against a reference file's own sky-axis WCS (ref_unit): CTYPE/
      !! CRVAL/CDELT/rotation/RADESYS/EQUINOX/LONPOLE/LATPOLE required
      !! equal, CRPIX required equal to the reference's own value minus
      !! crpix_shift1/2 (zero for a plain "must match exactly"
      !! comparison; nonzero when the candidate is checked against a
      !! footprint-grown/cropped output canvas -- match_cubes/
      !! reproject_cubes' own mode=union/intersection). NOT a general
      !! "are these two grids close enough" check -- a false match here
      !! would silently misalign downstream processing (reprojection
      !! skipped when it shouldn't be, or RM synthesis combining
      !! genuinely different sky positions), so every comparison is a
      !! tight absolute tolerance, and a file that is very close but
      !! outside it is correctly treated as NOT matching.
      integer, intent(in) :: cand_unit, cand_axis1, cand_axis2, nx_cand, ny_cand
      integer, intent(in) :: ref_unit, ref_axis1, ref_axis2
      double precision, intent(in) :: crpix_shift1, crpix_shift2
      integer, intent(in) :: nx_out, ny_out

      double precision, parameter :: tol_val = 1.0d-9
      double precision, parameter :: tol_rot = 1.0d-9
      double precision, parameter :: tol_frame = 1.0d-6
      character(len=68) :: ctype_c1, ctype_c2, ctype_r1, ctype_r2
      double precision :: crval_c1, crval_c2, crval_r1, crval_r2
      double precision :: crpix_c1, crpix_c2, crpix_r1, crpix_r2
      double precision :: cdelt_c1, cdelt_c2, cdelt_r1, cdelt_r2
      double precision :: pc_c(2,2), pc_r(2,2), crota_c1, crota_c2, crota_r1, crota_r2
      logical :: have_pc_c, have_pc_r
      character(len=8) :: radesys_c, radesys_r
      logical :: have_radesys_c, have_radesys_r, have_equinox_c, have_equinox_r
      double precision :: equinox_c, equinox_r, eq_eff_c, eq_eff_r
      logical :: have_lonpole_c, have_lonpole_r, have_latpole_c, have_latpole_r
      double precision :: lonpole_c, lonpole_r, latpole_c, latpole_r

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
      ! value minus the caller's own shift (zero for a plain match, the
      ! output-canvas offset for match_cubes/reproject_cubes' own
      ! mode=union/intersection).
      call get_axis_dval(cand_unit, 'CRPIX', cand_axis1, 1.0d0, crpix_c1)
      call get_axis_dval(cand_unit, 'CRPIX', cand_axis2, 1.0d0, crpix_c2)
      call get_axis_dval(ref_unit, 'CRPIX', ref_axis1, 1.0d0, crpix_r1)
      call get_axis_dval(ref_unit, 'CRPIX', ref_axis2, 1.0d0, crpix_r2)
      if (abs(crpix_c1-(crpix_r1-crpix_shift1)).gt.tol_val) return
      if (abs(crpix_c2-(crpix_r2-crpix_shift2)).gt.tol_val) return

      ! Rotation: PCi_j/CDi_j takes precedence over CROTA if either file
      ! has any entry present -- absent PC entries default to the FITS
      ! standard identity matrix, absent CROTA defaults to 0.
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

      ! RADESYS/EQUINOX: resolved through the FITS WCS standard's own
      ! default chain before comparing, not compared as raw keyword
      ! strings -- an absent RADESYS is not "unknown", it has a defined
      ! meaning. RADESYS absent: infer from EQUINOX (<1984 -> FK4,
      ! >=1984 -> FK5); both absent -> ICRS. EQUINOX only matters for
      ! FK4/FK5 (an equinox-referenced frame); absent there: infer from
      ! the resolved RADESYS (FK4 -> 1950.0, FK5 -> 2000.0).
      call get_global_sval(cand_unit, 'RADESYS', radesys_c, have_radesys_c)
      call get_global_sval(ref_unit, 'RADESYS', radesys_r, have_radesys_r)
      call get_global_dval(cand_unit, 'EQUINOX', equinox_c, have_equinox_c)
      call get_global_dval(ref_unit, 'EQUINOX', equinox_r, have_equinox_r)
      call resolve_radesys(radesys_c, have_radesys_c, equinox_c, have_equinox_c)
      call resolve_radesys(radesys_r, have_radesys_r, equinox_r, have_equinox_r)
      if (trim(radesys_c).ne.trim(radesys_r)) return
      if (trim(radesys_c).eq.'FK4' .or. trim(radesys_c).eq.'FK5') then
         call resolve_equinox(radesys_c, equinox_c, have_equinox_c, eq_eff_c)
         call resolve_equinox(radesys_r, equinox_r, have_equinox_r, eq_eff_r)
         if (abs(eq_eff_c-eq_eff_r).gt.tol_frame) return
      endif

      ! LONPOLE/LATPOLE: both default from CTYPE/CRVAL through a
      ! projection-dependent formula this check does not reimplement.
      ! Compared only when at least one side states a value explicitly:
      ! both absent is accepted (both use the same unevaluated default,
      ! since CTYPE/CRVAL are already required equal above); exactly one
      ! explicit is treated as a mismatch rather than guessed at.
      call get_global_dval(cand_unit, 'LONPOLE', lonpole_c, have_lonpole_c)
      call get_global_dval(ref_unit, 'LONPOLE', lonpole_r, have_lonpole_r)
      if (have_lonpole_c .or. have_lonpole_r) then
         if (.not.(have_lonpole_c .and. have_lonpole_r)) return
         if (abs(lonpole_c-lonpole_r).gt.tol_frame) return
      endif

      call get_global_dval(cand_unit, 'LATPOLE', latpole_c, have_latpole_c)
      call get_global_dval(ref_unit, 'LATPOLE', latpole_r, have_latpole_r)
      if (have_latpole_c .or. have_latpole_r) then
         if (.not.(have_latpole_c .and. have_latpole_r)) return
         if (abs(latpole_c-latpole_r).gt.tol_frame) return
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
      !! absent (the FITS standard's own default for that keyword -- 0
      !! for CRVAL/CROTA, 1 for CDELT/CRPIX -- passed in by the caller
      !! rather than hard-coded here).
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
      !! absent) rotation/scale matrix for the axis pair (axis1,axis2).
      !! have_any=.false. (m set to the identity) when neither convention
      !! has any entry present -- the FITS standard default.
      integer, intent(in) :: unit, axis1, axis2
      double precision, intent(out) :: m(2,2)
      logical, intent(out) :: have_any
      logical :: any_pc, any_cd
      double precision :: mcd(2,2)

      m = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2,2/))
      any_pc = .false.
      call get_matrix_entry_track(unit, 'PC', axis1, axis1, m(1,1), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis1, axis2, m(1,2), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis2, axis1, m(2,1), any_pc)
      call get_matrix_entry_track(unit, 'PC', axis2, axis2, m(2,2), any_pc)

      mcd = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2,2/))
      any_cd = .false.
      call get_matrix_entry_track(unit, 'CD', axis1, axis1, mcd(1,1), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis1, axis2, mcd(1,2), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis2, axis1, mcd(2,1), any_cd)
      call get_matrix_entry_track(unit, 'CD', axis2, axis2, mcd(2,2), any_cd)
      if (any_cd) m = mcd

      have_any = any_pc .or. any_cd
   end subroutine get_matrix_2x2

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

   subroutine get_global_sval(unit, keyword, val, have_val)
      !! Read a non-axis-suffixed string keyword (e.g. RADESYS).
      integer, intent(in) :: unit
      character(len=*), intent(in) :: keyword
      character(len=*), intent(out) :: val
      logical, intent(out) :: have_val
      character(len=68) :: comment
      integer :: fitsstat
      val = ' '
      fitsstat = 0
      call FTGKYS(unit, keyword, val, comment, fitsstat)
      have_val = (fitsstat.eq.0)
   end subroutine get_global_sval

   subroutine get_global_dval(unit, keyword, val, have_val)
      !! Read a non-axis-suffixed double keyword (e.g. EQUINOX/LONPOLE/
      !! LATPOLE).
      integer, intent(in) :: unit
      character(len=*), intent(in) :: keyword
      double precision, intent(out) :: val
      logical, intent(out) :: have_val
      character(len=68) :: comment
      integer :: fitsstat
      val = 0.0d0
      fitsstat = 0
      call FTGKYD(unit, keyword, val, comment, fitsstat)
      have_val = (fitsstat.eq.0)
   end subroutine get_global_dval

   subroutine resolve_radesys(radesys, have_radesys, equinox, have_equinox)
      !! In place: resolves radesys to its effective value through the
      !! FITS WCS standard's default chain (radesys given -> use it,
      !! upper-cased; else infer from equinox; else ICRS).
      character(len=*), intent(inout) :: radesys
      logical, intent(in) :: have_radesys
      double precision, intent(in) :: equinox
      logical, intent(in) :: have_equinox
      integer :: i
      if (have_radesys) then
         do i = 1, len_trim(radesys)
            if (radesys(i:i).ge.'a' .and. radesys(i:i).le.'z') then
               radesys(i:i) = achar(iachar(radesys(i:i)) - 32)
            endif
         enddo
         return
      endif
      if (have_equinox) then
         if (equinox.lt.1984.0d0) then
            radesys = 'FK4'
         else
            radesys = 'FK5'
         endif
      else
         radesys = 'ICRS'
      endif
   end subroutine resolve_radesys

   subroutine resolve_equinox(radesys, equinox, have_equinox, eq_eff)
      !! Effective EQUINOX for an already-resolved FK4/FK5 radesys: the
      !! keyword's own value if given, else that frame's own standard
      !! default epoch (FK4 -> B1950.0, FK5 -> J2000.0).
      character(len=*), intent(in) :: radesys
      double precision, intent(in) :: equinox
      logical, intent(in) :: have_equinox
      double precision, intent(out) :: eq_eff
      if (have_equinox) then
         eq_eff = equinox
      else if (trim(radesys).eq.'FK4') then
         eq_eff = 1950.0d0
      else
         eq_eff = 2000.0d0
      endif
   end subroutine resolve_equinox

end module wcs_match_mod
