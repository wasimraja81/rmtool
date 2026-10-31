program ast_icrs_to_frame
   !! Converts ICRS positions to another sky system with Starlink AST, the
   !! library reproject_cubes and match_cubes use for every sky conversion
   !! (scripts/compare_ast_astropy_frames.py compares the result with
   !! astropy). The target frame is built the way the tools build theirs:
   !! from FITS header cards (CTYPE, RADESYS, EQUINOX), so it carries the
   !! same attributes -- including any epoch AST assumes when the header
   !! gives none.
   !! Usage: ast_icrs_to_frame <points file> icrs|fk5|fk4|galactic
   !! Reads "ra dec" pairs in degrees, writes "lon lat" in degrees of the
   !! target system, one pair per line, to stdout.
   implicit none
   integer, external :: ast_skyframe, ast_convert, ast_fitschan, ast_read, ast_getframe
   integer, parameter :: ast__current = -1
   double precision, parameter :: deg2rad = 3.14159265358979323846d0 / 180.0d0
   integer :: status, f_icrs, f_tgt, fchan, wcs, cvt, n, ios, k, unit
   double precision, allocatable :: ra(:), dec(:), lon(:), lat(:)
   double precision :: x, y
   character(len=1024) :: infile
   character(len=16) :: system
   character(len=80) :: cards(16)
   integer :: ncard
   external :: ast_null

   call get_command_argument(1, infile)
   call get_command_argument(2, system)
   open(newunit=unit, file=trim(infile), status='old', action='read')
   n = 0
   do
      read(unit, *, iostat=ios) x, y
      if (ios.ne.0) exit
      n = n + 1
   enddo
   allocate(ra(n), dec(n), lon(n), lat(n))
   rewind(unit)
   do k = 1, n
      read(unit, *) ra(k), dec(k)
   enddo
   close(unit)
   ra = ra * deg2rad
   dec = dec * deg2rad

   ncard = 0
   call card("NAXIS   =                    2")
   call card("NAXIS1  =                   10")
   call card("NAXIS2  =                   10")
   select case (trim(system))
   case ('galactic')
      call card("CTYPE1  = 'GLON-CAR'")
      call card("CTYPE2  = 'GLAT-CAR'")
   case ('icrs')
      call card("CTYPE1  = 'RA---CAR'")
      call card("CTYPE2  = 'DEC--CAR'")
      call card("RADESYS = 'ICRS'")
   case ('fk5')
      call card("CTYPE1  = 'RA---CAR'")
      call card("CTYPE2  = 'DEC--CAR'")
      call card("RADESYS = 'FK5'")
      call card("EQUINOX =               2000.0")
   case ('fk4')
      call card("CTYPE1  = 'RA---CAR'")
      call card("CTYPE2  = 'DEC--CAR'")
      call card("RADESYS = 'FK4'")
      call card("EQUINOX =               1950.0")
   case ('fk5-epoch1950')
      ! confirming test: Equinox=2000 (so the FK5->FK4 step and the later
      ! precession step are both identity) but Epoch=1950 explicitly -- if
      ! AST's ICRS->FK5 result changes vs case 'fk5' (Equinox=Epoch=2000),
      ! that isolates the epoch-dependent HFK5Z/eraHfk5z spin term alone
      call card("CTYPE1  = 'RA---CAR'")
      call card("CTYPE2  = 'DEC--CAR'")
      call card("RADESYS = 'FK5'")
      call card("EQUINOX =               2000.0")
      call card("MJD-OBS =            33281.923")
   case ('fk5-1950')
      ! isolation diagnostic: FK5 label (IAU1976/Lieske precession, no
      ! E-terms, no Fricke equinox-bias correction) evaluated at B1950 --
      ! tests the precession step alone, separate from the FK4-specific
      ! equinox-bias+E-terms step
      call card("CTYPE1  = 'RA---CAR'")
      call card("CTYPE2  = 'DEC--CAR'")
      call card("RADESYS = 'FK5'")
      call card("EQUINOX =               1950.0")
   case default
      write(0, *) 'unknown system ', trim(system)
      stop 1
   end select
   call card("CRVAL1  =                  0.0")
   call card("CRVAL2  =                  0.0")
   call card("CRPIX1  =                  1.0")
   call card("CRPIX2  =                  1.0")
   call card("CDELT1  =                 -1.0")
   call card("CDELT2  =                  1.0")

   status = 0
   fchan = ast_fitschan(ast_null, ast_null, ' ', status)
   do k = 1, ncard
      call ast_putfits(fchan, cards(k), .false., status)
   enddo
   call ast_seti(fchan, 'Card', 1, status)
   wcs = ast_read(fchan, status)
   f_tgt = ast_getframe(wcs, ast__current, status)
   f_icrs = ast_skyframe('System=ICRS', status)
   cvt = ast_convert(f_icrs, f_tgt, ' ', status)
   call ast_tran2(cvt, n, ra, dec, .true., lon, lat, status)
   if (status.ne.0) then
      write(0, *) 'AST conversion failed, status ', status
      stop 1
   endif
   do k = 1, n
      write(*, '(ES24.16E3, 1X, ES24.16E3)') lon(k) / deg2rad, lat(k) / deg2rad
   enddo

contains

   subroutine card(text)
      character(len=*), intent(in) :: text
      ncard = ncard + 1
      cards(ncard) = text
   end subroutine card
end program ast_icrs_to_frame
