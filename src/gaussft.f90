module gaussft_mod
   !! Convolve a single 2D image plane from one elliptical-Gaussian PSF to
   !! another, via FFT-domain deconvolve-then-reconvolve: multiply
   !! FT(image) by FT(target beam)/FT(source beam), inverse-transform.
   !! Pure computation only -- no FITS I/O, no NaN/bad-data handling, no
   !! per-channel BMAJ/BMIN/BPA bookkeeping; those are a caller's job
   !! (planned: a main program mirroring reproject_cubes' own split
   !! between user-facing I/O and this kind of narrowly-scoped
   !! computational module). Multi-band needs nothing extra here --
   !! every plane already carries its own independent source PSF against
   !! one shared target PSF, whether it's the only plane being processed
   !! or one of many from several bands; "multi-band" is purely a matter
   !! of how many times a caller invokes convolve_to_beam and where it
   !! reads each call's bmaj_in/bmin_in/bpa_in from, not something this
   !! module needs to know about.
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
   !! significant figures, and round-tripping a point source through the
   !! full FFT/multiply/IFFT pipeline against MIRIAD's own gaufac-derived
   !! scaling factor (au2.gauss_factor, independently re-derived, both
   !! landing on the same figure).
   !!
   !! Split into plan/execute/destroy (rather than one self-contained
   !! call that plans its own FFTs, as the first version of this module
   !! did) so a caller can parallelise across planes with OpenMP: FFTW's
   !! planner functions are not thread-safe and must never run inside a
   !! parallel region, but a single plan, once created, is safe to
   !! EXECUTE concurrently from multiple threads via the "new-array
   !! execute" form (dfftw_execute_dft with explicit in/out arguments,
   !! as convolve_to_beam uses below) as long as each concurrent call
   !! supplies its own distinct arrays -- true here, since image/
   !! image_out/the internal work arrays are all local to each call.
   !! Verified directly (see this module's own test suite): 16 OpenMP
   !! threads calling convolve_to_beam concurrently against the SAME
   !! shared plan, each on its own image/beam pair, matches a serial
   !! run of the same 16 calls exactly.
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: plan_convolution, convolve_to_beam, destroy_convolution_plan,&
   &next_fast_fft_size

   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: deg2rad = pi/180.0_dp
   ! sqrt(8*ln2): standard FWHM -> Gaussian-sigma conversion factor.
   real(dp), parameter :: fwhm2sigma = 2.0_dp*sqrt(2.0_dp*log(2.0_dp))

   ! FFTW3 constants (from /usr/include/fftw3.f) -- declared directly
   ! rather than `include`d: fftw3.f is fixed-form Fortran 77 (same
   ! issue as AST_PAR in reproject_cubes.f90, see its own comment) and
   ! cannot be included into a free-form .f90 file directly.
   integer, parameter :: fftw_forward = -1
   integer, parameter :: fftw_backward = 1
   integer, parameter :: fftw_estimate = 64

contains

   subroutine plan_convolution(nx, ny, plan_fwd, plan_bwd, nx_pad, ny_pad)
      !! Create the FFTW plans for an nx-by-ny transform, once, to be
      !! reused by every subsequent convolve_to_beam call for planes of
      !! this same size (the common case -- every plane of a cube, and
      !! indeed of every band's cube in a multi-band run, shares one
      !! nx,ny). MUST be called serially, before any parallel region --
      !! FFTW's planner functions are not thread-safe. Pair with
      !! destroy_convolution_plan once every plane is done, also
      !! serially.
      !!
      !! nx_pad/ny_pad (out): the ACTUAL transform size the plans are for
      !! -- next_fast_fft_size(nx)/next_fast_fft_size(ny), never smaller
      !! than nx/ny. The plan (and every subsequent convolve_to_beam call
      !! using it) operates on this padded size, not the raw nx/ny --
      !! see next_fast_fft_size's own comment for why. The caller must
      !! pass nx_pad/ny_pad back into every convolve_to_beam call for
      !! this plan (an image smaller than the plan's own transform size
      !! is zero-padded internally, and the result cropped back).
      integer, intent(in) :: nx, ny
      integer(kind=8), intent(out) :: plan_fwd, plan_bwd
      integer, intent(out) :: nx_pad, ny_pad

      complex(dp), allocatable :: scratch(:,:)

      nx_pad = next_fast_fft_size(nx)
      ny_pad = next_fast_fft_size(ny)

      ! FFTW_ESTIMATE plans don't depend on the array CONTENTS, or even
      ! specifically on the memory used here, only the shape -- every
      ! actual convolve_to_beam call below uses the "new-array execute"
      ! form with its own arrays instead. FFTW's own documentation
      ! describes this as fully supported (correct regardless of
      ! alignment, though a non-ESTIMATE plan might not be as fast on a
      ! differently-aligned array than the one it was planned with --
      ! irrelevant for ESTIMATE, which never does alignment-specific
      ! optimisation in the first place).
      allocate(scratch(nx_pad, ny_pad))
      call dfftw_plan_dft_2d(plan_fwd, nx_pad, ny_pad, scratch, scratch, fftw_forward, fftw_estimate)
      call dfftw_plan_dft_2d(plan_bwd, nx_pad, ny_pad, scratch, scratch, fftw_backward, fftw_estimate)
      deallocate(scratch)
   end subroutine plan_convolution

   function next_fast_fft_size(n) result(m)
      !! Smallest m >= n whose only prime factors are 2, 3, 5, or 7 (a
      !! "7-smooth" number) -- FFTW has fast hard-coded codelets for
      !! these small factors (and their products), but falls back to a
      !! much slower generic algorithm (Bluestein's, effectively an
      !! extra embedded FFT of its own) for any size with a large prime
      !! factor. Real-world example that motivated this (not
      !! hard-coded, just the case that exposed the problem): a genuine
      !! ASKAP cube at nx=ny=4501 = 7 x 643, and 643 is prime -- a
      !! single 4501x4501 complex FFT took 2.76s (measured directly,
      !! FFTW_ESTIMATE) vs 1.26s for the next 7-smooth size up (4608 =
      !! 2^9 x 3^2), a >2x slowdown from one large prime factor alone.
      !! This function makes NO assumption about what n will be -- it
      !! searches upward from n for ANY input, so every image size gets
      !! the same treatment, not just this one real-data case.
      integer, intent(in) :: n
      integer :: m, r
      integer, parameter :: small_primes(4) = (/2, 3, 5, 7/)
      integer :: p

      m = max(n, 1)
      do
         r = m
         do p = 1, size(small_primes)
            do while (mod(r, small_primes(p)).eq.0)
               r = r/small_primes(p)
            enddo
         enddo
         if (r.eq.1) exit
         m = m + 1
      enddo
   end function next_fast_fft_size

   subroutine destroy_convolution_plan(plan_fwd, plan_bwd)
      integer(kind=8), intent(inout) :: plan_fwd, plan_bwd

      call dfftw_destroy_plan(plan_fwd)
      call dfftw_destroy_plan(plan_bwd)
   end subroutine destroy_convolution_plan

   subroutine convolve_to_beam(plan_fwd, plan_bwd, image, nx, ny, nx_pad,&
   &ny_pad, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa, image_out,&
   &status)
      !! plan_fwd/plan_bwd: from plan_convolution, already created,
      !! describing an nx_pad-by-ny_pad transform (not checked here --
      !! passing plans for a different size is undefined behaviour, same
      !! as FFTW's own new-array execute contract). image(nx,ny): input
      !! plane, at its own true extent -- NOT necessarily nx_pad,ny_pad
      !! (see plan_convolution/next_fast_fft_size: the FFT itself always
      !! runs at the padded size for speed, never at a raw size with a
      !! large prime factor). Zero-padded into the top-left corner of an
      !! (nx_pad,ny_pad) work array before the forward transform, and
      !! cropped back to the original (nx,ny) from the same corner after
      !! the inverse transform -- the standard "linear convolution via
      !! zero-padded FFT" technique: since the real image data no longer
      !! occupies the full periodic extent the DFT assumes, this also
      !! reduces the edge wraparound a same-size (unpadded) circular
      !! convolution would have, rather than introducing any new
      !! approximation. dx,dy: pixel scale, DEGREES (same convention as
      !! CDELT1/2 -- converted to radians internally, alongside the beam
      !! parameters below, so u,v end up in cycles per radian, matching
      !! sx/sy/sx_in/sy_in) -- unchanged by padding, since padding only
      !! extends the field at the SAME pixel scale, never resamples it.
      !! bmaj_in/bmin_in/bpa_in: THIS plane's own (native/source) PSF,
      !! degrees, standard FITS BMAJ/BMIN/BPA convention (BPA measured
      !! the same way the input header defines it -- this module does no
      !! coordinate-system reasoning of its own, it just rotates by the
      !! angle it's given). bmaj/bmin/bpa: the TARGET PSF to convolve
      !! to, same convention -- shared across every call for a
      !! common-resolution run, whatever plane or band each call's
      !! image/source PSF came from. image_out(nx,ny): the convolved
      !! plane, same extent as image. status: 0 on success (reserved for
      !! future use -- this module does not itself judge whether
      !! bmaj/bmin/bpa is a sensible request relative to
      !! bmaj_in/bmin_in/bpa_in; that policy call belongs to the caller,
      !! not this computation). Thread-safe: see this module's own
      !! header comment.
      integer(kind=8), intent(in) :: plan_fwd, plan_bwd
      integer, intent(in) :: nx, ny, nx_pad, ny_pad
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

      allocate(u(nx_pad), v(ny_pad))
      call build_fftfreq(nx_pad, dx_rad, u)
      call build_fftfreq(ny_pad, dy_rad, v)

      allocate(g_final(nx_pad, ny_pad))
      do iy = 1, ny_pad
         do ix = 1, nx_pad
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

      allocate(cimg(nx_pad, ny_pad))
      cimg = cmplx(0.0_dp, 0.0_dp, dp)
      cimg(1:nx, 1:ny) = cmplx(image, 0.0_dp, dp)
      call dfftw_execute_dft(plan_fwd, cimg, cimg)

      cimg = cimg*g_final
      deallocate(g_final)

      call dfftw_execute_dft(plan_bwd, cimg, cimg)

      ! FFTW's transforms are unnormalised (forward then backward scales
      ! the result by nx_pad*ny_pad, same convention as numpy.fft.fft2/
      ! ifft2 -- numpy just applies the 1/N inside ifft2 for you; FFTW
      ! leaves it to the caller, matching its own documented convention).
      ! Cropped back from the same top-left corner the image was placed
      ! at above.
      image_out = real(cimg(1:nx, 1:ny), dp) / real(nx_pad*ny_pad, dp)
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
