module gaussft_mod
   !! Convolve a single 2D image plane from one elliptical-Gaussian PSF to
   !! another, via FFT-domain deconvolve-then-reconvolve: multiply
   !! FT(image) by FT(target beam)/FT(source beam), inverse-transform.
   !! Pure computation only -- no FITS I/O, no NaN/bad-data handling, no
   !! per-channel BMAJ/BMIN/BPA bookkeeping; those are a caller's job
   !! (planned: a main program mirroring reproject_cubes' own split
   !! between user-facing I/O and this kind of narrowly-scoped
   !! computational module).
   !!
   !! Corrected from an earlier version (src/gaussft.f, the original
   !! Fortran77 prototype -- and its direct Python port, both upstream
   !! racs_tools/gaussft.py and the local mirror src/gaussft.py, which
   !! carry the identical formula): those computed the amplitude of each
   !! Gaussian's 2D Fourier transform as sqrt(2*pi*sigma_x*sigma_y). The
   !! correct closed-form amplitude, for the "ordinary frequency" FT
   !! convention this scheme uses throughout (u,v in cycles per unit
   !! length -- matching FFTW's own DFT frequency indexing, see
   !! build_fftfreq below), is 2*pi*sigma_x*sigma_y: a 2D Gaussian is
   !! separable, and its FT is the PRODUCT of two independent 1D
   !! transforms, each contributing its own sqrt(2*pi)*sigma factor --
   !! sqrt(2*pi)*sigma_x * sqrt(2*pi)*sigma_y = 2*pi*sigma_x*sigma_y, not
   !! sqrt(2*pi*sigma_x*sigma_y) (which is dimensionally wrong for a 2D
   !! integral besides -- sqrt(sigma_x*sigma_y) is a length, but a 2D
   !! integral has units of length^2). Verified several independent ways
   !! before this rewrite: closed-form derivation, direct numerical
   !! integration matching the true 2D Gaussian integral to 15
   !! significant figures, and round-tripping an actual point source
   !! through the full FFT/multiply/IFFT pipeline against MIRIAD's own
   !! gaufac-derived scaling factor (au2.gauss_factor, also in this
   !! repo, both landing on the same figure independently).
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: convolve_to_beam

   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: deg2rad = pi/180.0_dp
   ! sqrt(8*ln2): standard FWHM -> Gaussian-sigma conversion factor.
   real(dp), parameter :: fwhm2sigma = 2.0_dp*sqrt(2.0_dp*log(2.0_dp))

   ! FFTW3 constants (from /usr/include/fftw3.f) -- declared directly
   ! rather than `include`d: fftw3.f is fixed-form Fortran 77 (same
   ! issue as AST_PAR, see reproject_cubes.f90's own comment on it) and
   ! cannot be included into a free-form .f90 file directly.
   integer, parameter :: fftw_forward = -1
   integer, parameter :: fftw_backward = 1
   integer, parameter :: fftw_estimate = 64

contains

   subroutine convolve_to_beam(image, nx, ny, dx, dy,&
   &bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa, image_out, status)
      !! image(nx,ny): input plane. dx,dy: pixel scale, DEGREES (same
      !! convention as CDELT1/2 -- converted to radians internally,
      !! alongside the beam parameters below, so u,v end up in cycles
      !! per radian, matching sx/sy/sx_in/sy_in). bmaj_in/bmin_in/bpa_in:
      !! the plane's OWN (native/source) PSF, degrees, standard FITS
      !! BMAJ/BMIN/BPA convention (BPA measured the same way the input
      !! header defines it -- this module does no coordinate-system
      !! reasoning of its own, it just rotates by the angle it's given).
      !! bmaj/bmin/bpa: the TARGET PSF to convolve to, same convention.
      !! image_out(nx,ny): the convolved plane. status: 0 on success
      !! (reserved for future use -- e.g. FFTW failures -- this module
      !! does not itself judge whether bmaj/bmin/bpa is a sensible
      !! request relative to bmaj_in/bmin_in/bpa_in; that policy call
      !! belongs to the caller, not this computation).
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: image(nx, ny)
      real(dp), intent(in) :: dx, dy
      real(dp), intent(in) :: bmaj_in, bmin_in, bpa_in
      real(dp), intent(in) :: bmaj, bmin, bpa
      real(dp), intent(out) :: image_out(nx, ny)
      integer, intent(out) :: status

      real(dp) :: sx, sy, sx_in, sy_in
      real(dp) :: bpa_rad, bpa_in_rad, dx_rad, dy_rad
      real(dp) :: cos_bpa, sin_bpa, cos_bpa_in, sin_bpa_in
      real(dp) :: g_amp, dg_amp, g_ratio
      real(dp) :: ur, vr, ur_in, vr_in, g_arg, dg_arg
      real(dp), allocatable :: u(:), v(:)
      complex(dp), allocatable :: cimg(:,:), g_final(:,:)
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: ix, iy

      status = 0

      dx_rad = dx*deg2rad
      dy_rad = dy*deg2rad
      bpa_rad = bpa*deg2rad
      bpa_in_rad = bpa_in*deg2rad
      sx = (bmaj*deg2rad)/fwhm2sigma
      sy = (bmin*deg2rad)/fwhm2sigma
      sx_in = (bmaj_in*deg2rad)/fwhm2sigma
      sy_in = (bmin_in*deg2rad)/fwhm2sigma

      g_amp = 2.0_dp*pi*sx*sy
      dg_amp = 2.0_dp*pi*sx_in*sy_in
      g_ratio = g_amp/dg_amp
      cos_bpa = cos(bpa_rad)
      sin_bpa = sin(bpa_rad)
      cos_bpa_in = cos(bpa_in_rad)
      sin_bpa_in = sin(bpa_in_rad)

      allocate(u(nx), v(ny))
      call build_fftfreq(nx, dx_rad, u)
      call build_fftfreq(ny, dy_rad, v)

      allocate(g_final(nx, ny))
      do iy = 1, ny
         do ix = 1, nx
            ur = u(ix)*cos_bpa - v(iy)*sin_bpa
            vr = u(ix)*sin_bpa + v(iy)*cos_bpa
            g_arg = -2.0_dp*pi**2 * ((sx*ur)**2 + (sy*vr)**2)

            ur_in = u(ix)*cos_bpa_in - v(iy)*sin_bpa_in
            vr_in = u(ix)*sin_bpa_in + v(iy)*cos_bpa_in
            dg_arg = -2.0_dp*pi**2 * ((sx_in*ur_in)**2 + (sy_in*vr_in)**2)

            g_final(ix,iy) = g_ratio * exp(cmplx(g_arg - dg_arg, 0.0_dp, dp))
         enddo
      enddo
      deallocate(u, v)

      allocate(cimg(nx, ny))
      cimg = cmplx(image, 0.0_dp, dp)
      call dfftw_plan_dft_2d(plan_fwd, nx, ny, cimg, cimg, fftw_forward, fftw_estimate)
      call dfftw_execute_dft(plan_fwd, cimg, cimg)
      call dfftw_destroy_plan(plan_fwd)

      cimg = cimg*g_final
      deallocate(g_final)

      call dfftw_plan_dft_2d(plan_bwd, nx, ny, cimg, cimg, fftw_backward, fftw_estimate)
      call dfftw_execute_dft(plan_bwd, cimg, cimg)
      call dfftw_destroy_plan(plan_bwd)

      ! FFTW's transforms are unnormalised (forward then backward scales
      ! the result by nx*ny, same convention as numpy.fft.fft2/ifft2 --
      ! numpy just applies the 1/N inside ifft2 for you; FFTW leaves it
      ! to the caller, matching its own documented convention).
      image_out = real(cimg, dp) / real(nx*ny, dp)
      deallocate(cimg)
   end subroutine convolve_to_beam

   subroutine build_fftfreq(n, d, freq)
      !! Same frequency layout as numpy.fft.fftfreq(n, d=d): index k
      !! (0-based) maps to k/(n*d) for k < ceil(n/2), and (k-n)/(n*d)
      !! for k >= ceil(n/2) -- the standard DFT "wraparound" frequency
      !! ordering, matching how FFTW itself indexes its output. Built by
      !! hand since FFTW's own API has no fftfreq-equivalent helper.
      integer, intent(in) :: n
      real(dp), intent(in) :: d
      real(dp), intent(out) :: freq(n)
      integer :: k, half

      half = (n+1)/2
      do k = 0, half-1
         freq(k+1) = real(k, dp) / (real(n, dp)*d)
      enddo
      do k = half, n-1
         freq(k+1) = real(k-n, dp) / (real(n, dp)*d)
      enddo
   end subroutine build_fftfreq

end module gaussft_mod
