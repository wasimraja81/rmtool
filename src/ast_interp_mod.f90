module ast_interp_mod
   !! Resolves the interp=/interp_params= keys (docs/dev/
   !! REPROJECT_INTERP_OPTIONS_PLAN.md) shared by reproject_cubes.f90
   !! and match_cubes.f90's four independent ast_resampler call sites
   !! (that file does not call into reproject_cubes.f90's code -- see
   !! the design doc's call-site inventory) into the AST scheme
   !! constant and params(:) array ast_resampler expects, in one place,
   !! so the four call sites cannot drift the way AST's own SincSinc/
   !! SincCos drifted from each other (REPROJECT_INTERP_OPTIONS_PLAN.md
   !! section 3).
   !!
   !! Scheme constants below are AST's own (checked against the
   !! installed /usr/include/AST_PAR, not assumed), kept PRIVATE except
   !! ast__ukern1 -- callers never need any other raw constant, only
   !! resolve_interp's resolved scheme_const output, so there is no
   !! reason to export a name that could collide with a caller's own
   !! local declaration of the same constant. ast__ukern1 is the one
   !! exception (section 6d): callers need it to resolve their own
   !! local finterp.
   !!
   !! params(:) semantics (traced directly against AST's own mapping.c
   !! dispatcher, not SUN/211 prose alone -- see the design doc section
   !! 5 table): params(1) is a neighbour/support-radius count shared by
   !! every 1-D-kernel scheme (<=0 means "let AST choose
   !! automatically"); sincsinc/sinccos/sincgauss/gauss/sombcos
   !! additionally read params(2) as their own shape value (a sinc
   !! width or a Gaussian FWHM). nearest/linear/blockave/sinc/somb do
   !! not use params(2).
   !!
   !! Tier 2 (a corrected Lanczos kernel via AST__UKERN1, design doc
   !! section 6): lanczos2/lanczos3/lanczos share one callback
   !! (lanczos_kernel, below this module -- a plain Fortran EXTERNAL
   !! subroutine, not bind(C); the finterp calling convention was
   !! verified empirically against AST's own Fortran glue, not assumed
   !! from the C-facing synopsis -- see section 6b and
   !! scratch/build_check/test_ukern1_abi.f90/test_ukern1_ptr.f90).
   !!
   !! resolve_interp deliberately does NOT return a finterp pointer --
   !! see section 6d. Passing a procedure(),pointer as a dummy argument
   !! into a contains-ed procedure that then shares it into its own
   !! !$omp parallel region segfaulted (scratch/build_check/
   !! test_ukern1_internal_proc_bug.f90 has the isolated reproducer
   !! and fix), with every real-WCS-Mapping/array-shape/thread-count
   !! alternative explanation ruled out first. Callers instead compare
   !! their own already-safe interp_scheme against the ast__ukern1
   !! constant below (exported public for exactly this) and resolve
   !! finterp as a LOCAL variable, immediately before their own
   !! parallel region, never as anyone's dummy argument.
   implicit none
   private
   public :: resolve_interp, ast__ukern1, interp_exact

   integer, parameter :: ast__ukern1    = 1
   ! Not an AST scheme: interp=exact is the flux-conserving overlap
   ! method of exact_overlap_mod; callers branch on this value instead
   ! of calling ast_resampler. Negative so it can never equal an AST__
   ! constant.
   integer, parameter :: interp_exact   = -1
   integer, parameter :: ast__nearest   = 4
   integer, parameter :: ast__linear    = 5
   integer, parameter :: ast__sinc      = 6
   integer, parameter :: ast__sincsinc  = 7
   integer, parameter :: ast__sinccos   = 8
   integer, parameter :: ast__sincgauss = 9
   integer, parameter :: ast__blockave  = 10
   integer, parameter :: ast__gauss     = 11
   integer, parameter :: ast__somb      = 12
   integer, parameter :: ast__sombcos   = 13

   !! Default kernel widths for the three schemes whose width AST
   !! cannot choose itself (UR-19, docs/dev/BRANCH_STATUS_update-
   !! reproject.md). In INPUT-GRID PIXELS, which is AST's own unit for
   !! these -- a Gaussian FWHM for gauss/sincgauss, a box half-width
   !! for blockave. Measured, not guessed, by
   !! tests/check_interp_width_defaults.py: each is the width with the
   !! lowest WORST-CASE error over 3 beam samplings x 4 output
   !! geometries, since a default cannot know the caller's geometry and
   !! the per-geometry optima disagree. gauss and blockave degrade
   !! monotonically as the width grows, so theirs is the narrowest
   !! width that still resamples; sincgauss's Gaussian tapers a sinc
   !! rather than doing the interpolating, so it runs the other way and
   !! plateaus wide. None of the three scales with the beam -- an
   !! earlier design deriving them from BMAJ and the pixel scale was
   !! dropped because the measurement found no beam-proportional
   !! optimum to derive.
   double precision, parameter :: default_gauss_fwhm_px         = 0.8d0
   double precision, parameter :: default_sincgauss_fwhm_px     = 8.0d0
   double precision, parameter :: default_blockave_halfwidth_px = 1.0d0

contains

   subroutine resolve_interp(name, params_str, scheme_const, params_arr,&
   &status)
      !! name: one of nearest/linear/sinc/sincsinc/sinccos/sincgauss/
      !! blockave/gauss/somb/sombcos/lanczos2/lanczos3/lanczos,
      !! case-insensitive (empty string defaults to lanczos3 -- the
      !! CLI/cfg layer's own default since 2026-10-13, UR-12's own
      !! accuracy benchmark; this routine does not duplicate that
      !! default, it only resolves whatever name it is given plus
      !! the empty-string case for direct callability from a test
      !! harness without going through full CLI parsing).
      !! params_str: comma-separated numbers, same csv style as this
      !! project's other list-valued keys (infiles=, badchan_file=,
      !! ...) -- semantics depend on name: 0/1/2 numbers, any unused
      !! slot meaning "let AST choose automatically", for every name
      !! except lanczos2/lanczos3 (interp_params rejected -- see
      !! resolve_lanczos) and lanczos (exactly 2 required).
      !! status: unchanged on success; set non-zero with an ERROR
      !! written to stdout on an unknown name or a malformed
      !! params_str, same convention as every other key this project
      !! validates (mode=, stages=, ...).
      !!
      !! Does NOT resolve finterp (see this module's own header
      !! comment, section 6d) -- callers compare scheme_const against
      !! this module's public ast__ukern1 and resolve finterp as a
      !! LOCAL variable themselves, immediately before their own
      !! !$omp parallel region.
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: params_str
      integer, intent(out) :: scheme_const
      double precision, intent(out) :: params_arr(2)
      integer, intent(inout) :: status

      character(len=16) :: key
      integer :: n_given, ios
      double precision :: p1, p2

      if (status.ne.0) return

      params_arr(1) = 0.0d0
      params_arr(2) = 0.0d0

      key = trim(adjustl(name))
      call lowercase_inplace(key)
      if (len_trim(key).eq.0) key = 'lanczos3'

      select case (trim(key))
      case ('nearest');   scheme_const = ast__nearest
      case ('linear');    scheme_const = ast__linear
      case ('sinc');      scheme_const = ast__sinc
      case ('sincsinc');  scheme_const = ast__sincsinc
      case ('sinccos');   scheme_const = ast__sinccos
      case ('sincgauss'); scheme_const = ast__sincgauss
      case ('blockave');  scheme_const = ast__blockave
      case ('gauss');     scheme_const = ast__gauss
      case ('somb');      scheme_const = ast__somb
      case ('sombcos');   scheme_const = ast__sombcos
      case ('lanczos2', 'lanczos3', 'lanczos'); scheme_const = ast__ukern1
      case ('exact')
         scheme_const = interp_exact
         if (count_csv_items(params_str).gt.0) then
            write(*,*) 'ERROR: interp=exact takes no interp_params= (it is'//&
            &' not a kernel; got "', trim(params_str), '")'
            status = -1
         endif
         return
      case default
         write(*,*) 'ERROR: unrecognised interp "', trim(name), '" -- expected',&
         &' nearest, linear, sinc, sincsinc, sinccos, sincgauss, blockave,'//&
         &' gauss, somb, sombcos, lanczos2, lanczos3, lanczos, or exact'
         status = -1
         return
      end select

      ! interp=sincsinc carries a known AST cutoff-constant mismatch
      ! (REPROJECT_INTERP_OPTIONS_PLAN.md section 3: truncates at 2/pi
      ! of the window's own peak, not at its zero) -- offered anyway
      ! (the point is an informed choice, not deciding on anyone's
      ! behalf), flagged here so the warning cannot be missed by using
      ! this routine outside the --help text that also documents it.
      if (trim(key).eq.'sincsinc') then
         write(*,*) 'WARNING: interp=sincsinc carries a known AST cutoff'//&
         &' discontinuity (see docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md'//&
         &' section 3) -- interp=sinccos is the clean alternative with'//&
         &' similar character.'
      endif

      if (trim(key).eq.'lanczos2' .or. trim(key).eq.'lanczos3' .or.&
      &trim(key).eq.'lanczos') then
         call resolve_lanczos(trim(key), params_str, params_arr, status)
         return
      endif

      n_given = count_csv_items(params_str)
      if (n_given.eq.0) then
         ! interp=gauss/sincgauss/blockave have no usable "auto"
         ! default the way sinc/sinccos/sincsinc/somb/sombcos's own
         ! neighbour count does -- confirmed against AST's own
         ! mapping.c directly (docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md,
         ! UR-14/UR-12): gauss/sincgauss clamp an omitted FWHM to a
         ! floor of 0.1 pixels (far too narrow -- the kernel underflows
         ! to zero at almost any sub-pixel offset), and blockave's own
         ! box half-width is not auto-derived from 0 at all. Until
         ! UR-19 this was a hard error, which left three of the twelve
         ! schemes unreachable unless the caller already knew a width
         ! to pass. The measured defaults above are supplied instead,
         ! and reported, so the omitted-params case resamples rather
         ! than refusing while the value in use stays on the record.
         if (trim(key).eq.'gauss') then
            params_arr(2) = default_gauss_fwhm_px
         else if (trim(key).eq.'sincgauss') then
            params_arr(2) = default_sincgauss_fwhm_px
         else if (trim(key).eq.'blockave') then
            params_arr(1) = default_blockave_halfwidth_px
         else
            return
         endif
         call report_default_params(trim(key), params_arr)
         return
      endif
      if (n_given.gt.2) then
         write(*,*) 'ERROR: interp_params takes at most 2 numbers, got ',&
         &n_given, ' in "', trim(params_str), '"'
         status = -1
         return
      endif

      call nth_csv_double(params_str, 1, p1, ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: interp_params'' first value is not a number: "',&
         &trim(params_str), '"'
         status = -1
         return
      endif
      params_arr(1) = p1

      if (n_given.eq.2) then
         call nth_csv_double(params_str, 2, p2, ios)
         if (ios.ne.0) then
            write(*,*) 'ERROR: interp_params'' second value is not a number: "',&
            &trim(params_str), '"'
            status = -1
            return
         endif
         params_arr(2) = p2
      endif
   end subroutine resolve_interp

   subroutine report_default_params(key, params_arr)
      !! Names the width this run is using when the caller omitted
      !! interp_params=. Printed rather than left silent because the
      !! width changes the resampled values, so a run whose numbers are
      !! later compared against another's needs the value it used to be
      !! recoverable from its own log.
      character(len=*), intent(in) :: key
      double precision, intent(in) :: params_arr(2)

      if (trim(key).eq.'blockave') then
         write(*,'(A,F6.2,A)') ' INFO: interp=blockave with no interp_params= --'//&
         &' using a box half-width of ', params_arr(1), ' pixels (default).'
      else
         write(*,'(A,A,A,F6.2,A)') ' INFO: interp=', trim(key), ' with no'//&
         &' interp_params= -- using a Gaussian FWHM of ', params_arr(2),&
         &' pixels (default).'
      endif
      write(*,*) '      Override with interp_params= (see --help); this'//&
      &' default is in pixels and does not scale with the beam.'
   end subroutine report_default_params

   subroutine resolve_lanczos(key, params_str, params_arr, status)
      !! lanczos2/lanczos3: a fixed at 2.0/3.0, support fixed at
      !! ceil(a) -- interp_params given at all is an error. Letting
      !! support drift below ceil(a) here would mean "lanczos3" is
      !! sometimes not Lanczos-3, with no offsetting benefit (anyone
      !! who wants that trade-off can ask for it explicitly via
      !! interp=lanczos instead) -- see REPROJECT_INTERP_OPTIONS_PLAN.md
      !! section 6c.
      !! lanczos: general form, requires interp_params=<support>,<a>
      !! (exactly 2 values, no default a -- none exists). support<=0
      !! means "use ceil(a)", the correct default (a smaller support
      !! silently drops contributions the kernel's own formula says
      !! are nonzero; a larger one only wastes compute asking the
      !! kernel for points it correctly evaluates to zero -- section
      !! 6b). An explicit smaller support is still accepted, as a
      !! deliberate speed/accuracy trade-off -- the same escape hatch
      !! sincsinc/sinccos/sincgauss already allow.
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: params_str
      double precision, intent(out) :: params_arr(2)
      integer, intent(inout) :: status

      integer :: n_given, ios
      double precision :: a, support, p1, p2

      if (status.ne.0) return

      n_given = count_csv_items(params_str)

      if (key.eq.'lanczos2' .or. key.eq.'lanczos3') then
         if (n_given.gt.0) then
            write(*,*) 'ERROR: interp_params is not accepted for interp="',&
            &trim(key), '" (its support radius is fixed) -- use'//&
            &' interp=lanczos with interp_params=<support>,<a> for an'//&
            &' explicit a/support combination'
            status = -1
            return
         endif
         if (key.eq.'lanczos2') then
            a = 2.0d0
         else
            a = 3.0d0
         endif
         support = ceiling(a)
      else
         if (n_given.ne.2) then
            write(*,*) 'ERROR: interp=lanczos requires interp_params=',&
            &'<support>,<a> (got ', n_given, ' value(s) in "',&
            &trim(params_str), '") -- support<=0 means use ceil(a)'
            status = -1
            return
         endif
         call nth_csv_double(params_str, 1, p1, ios)
         if (ios.ne.0) then
            write(*,*) 'ERROR: interp_params'' first value (support) is'//&
            &' not a number: "', trim(params_str), '"'
            status = -1
            return
         endif
         call nth_csv_double(params_str, 2, p2, ios)
         if (ios.ne.0) then
            write(*,*) 'ERROR: interp_params'' second value (a) is not'//&
            &' a number: "', trim(params_str), '"'
            status = -1
            return
         endif
         a = p2
         if (a.le.0.0d0) then
            write(*,*) 'ERROR: interp=lanczos requires a > 0, got ', a
            status = -1
            return
         endif
         if (p1.le.0.0d0) then
            support = ceiling(a)
         else
            support = p1
         endif
      endif

      params_arr(1) = support
      params_arr(2) = a
   end subroutine resolve_lanczos

   subroutine lowercase_inplace(s)
      character(len=*), intent(inout) :: s
      integer :: i, c
      do i = 1, len_trim(s)
         c = iachar(s(i:i))
         if (c.ge.iachar('A') .and. c.le.iachar('Z')) s(i:i) = achar(c + 32)
      enddo
   end subroutine lowercase_inplace

   function count_csv_items(str) result(n)
      !! Number of comma-separated items in str (0 for an empty/blank
      !! string). Small, self-contained double-array parser, separate
      !! from reproject_cubes.f90's/match_cubes.f90's own local
      !! cfg_csv_count/cfg_csv_get_item (each file has its own copy of
      !! those already -- see REPROJECT_INTERP_OPTIONS_PLAN.md's call-
      !! site inventory for that pre-existing duplication) -- those
      !! parse arbitrary-length filename lists; this parses at most 2
      !! numbers, a different and much narrower job not worth sharing
      !! code with the filename-list case.
      character(len=*), intent(in) :: str
      integer :: n, i
      character(len=len(str)) :: work
      n = 0
      work = adjustl(trim(str))
      if (len_trim(work).eq.0) return
      n = 1
      do i = 1, len_trim(work)
         if (work(i:i).eq.',') n = n + 1
      enddo
   end function count_csv_items

   subroutine nth_csv_double(str, idx, val, ios)
      character(len=*), intent(in) :: str
      integer, intent(in) :: idx
      double precision, intent(out) :: val
      integer, intent(out) :: ios
      character(len=len(str)) :: work, token
      integer :: item_no, start_pos, comma_pos

      work = adjustl(trim(str))
      item_no = 0
      start_pos = 1
      ios = -1
      do while (start_pos.le.len_trim(work))
         comma_pos = index(work(start_pos:), ',')
         if (comma_pos.eq.0) then
            token = adjustl(work(start_pos:len_trim(work)))
            item_no = item_no + 1
            if (item_no.eq.idx) then
               read(token, *, iostat=ios) val
            endif
            exit
         else
            token = adjustl(work(start_pos:start_pos + comma_pos - 2))
            item_no = item_no + 1
            if (item_no.eq.idx) then
               read(token, *, iostat=ios) val
               exit
            endif
            start_pos = start_pos + comma_pos
         endif
      enddo
   end subroutine nth_csv_double

end module ast_interp_mod

subroutine lanczos_kernel(offset, params, flags, value, status)
   !! AST__UKERN1 callback for interp=lanczos2/lanczos3/lanczos
   !! (REPROJECT_INTERP_OPTIONS_PLAN.md section 6). Deliberately a
   !! plain top-level Fortran subroutine, not a procedure inside
   !! ast_interp_mod -- the finterp calling convention this matches
   !! (SUN/211's 5-argument Fortran template, OFFSET/PARAMS/FLAGS/
   !! VALUE/STATUS) was verified empirically against a bare external
   !! subroutine (scratch/build_check/test_ukern1_abi.f90/
   !! test_ukern1_ptr.f90), not a module procedure; staying in that
   !! exact, proven form avoids introducing an unverified module-
   !! procedure-as-external-argument question on top of it. AST's own
   !! dispatch reads params(1) itself (the neighbour count, section
   !! 6b) -- never reaches here. params(2) is "a", read only by us;
   !! resolve_lanczos (ast_interp_mod, above) is what sets both.
   !!
   !! L(x) = sinc(x)*sinc(x/a) for |x|<a, 0 otherwise, where
   !! sinc(y) = sin(pi*y)/(pi*y) -- the standard Lanczos-a kernel,
   !! vanishing continuously at |x|=a (unlike AST's own SincSinc,
   !! section 3).
   implicit none
   double precision, intent(in) :: offset
   double precision, intent(in) :: params(*)
   integer, intent(in) :: flags
   double precision, intent(out) :: value
   integer, intent(inout) :: status

   double precision, parameter :: pi = 3.14159265358979323846d0
   double precision :: a, ao, po, poa

   a = params(2)
   ao = abs(offset)

   if (ao.ge.a) then
      value = 0.0d0
      return
   endif
   if (ao.lt.1.0d-12) then
      value = 1.0d0
      return
   endif

   po = pi * offset
   poa = po / a
   value = (sin(po) / po) * (sin(poa) / poa)
end subroutine lanczos_kernel
