module gaussft_mod
   !! Convolve a single 2D image plane from one elliptical-Gaussian PSF to
   !! another, via FFT-domain deconvolve-then-reconvolve: multiply
   !! FT(image) by FT(target beam)/FT(source beam), inverse-transform.
   !! Pure computation only -- no FITS I/O, no per-channel BMAJ/BMIN/BPA
   !! bookkeeping; those are a caller's job (planned: a main program
   !! mirroring reproject_cubes' own split between user-facing I/O and
   !! this kind of narrowly-scoped computational module). NaN/bad-data
   !! handling IS done here (see convolve_to_beam's own comment): real
   !! ASKAP data is genuinely partially blanked (primary-beam edge
   !! NaNs), and a single NaN anywhere in a plane would otherwise poison
   !! the entire FFT-convolved output, since the FFT is a global
   !! operation over the whole plane. Multi-band needs nothing extra here --
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
   !! did): FFTW's planner functions are not thread-safe and must never
   !! run inside a parallel region, so plan_convolution/
   !! destroy_convolution_plan must always be called serially --
   !! unchanged since this module's very first version.
   !!
   !! Parallelisation strategy (T28, docs/dev/MULTI_BAND_TOMOGRAPHY_
   !! PLAN.md) -- CURRENT: convolve_to_beam is called SERIALLY, once per
   !! plane, by its caller; the plan itself (see plan_convolution's own
   !! nthreads argument) executes each transform using multiple threads
   !! INTERNALLY (sfftw_plan_with_nthreads), and every other elementwise
   !! step in convolve_to_beam is parallelised with its own !$omp
   !! construct. Chosen over the alternative (many caller threads each
   !! independently calling convolve_to_beam on their own plane) because
   !! the alternative's per-thread memory cost multiplies by thread
   !! count (T26 found this the hard way, on real ASKAP-scale data) --
   !! this module's own working arrays are large enough, at real
   !! resolution, that this dominates over near-linear-vs-sub-linear
   !! scaling differences between the two strategies.
   !!
   !! HISTORICAL (pre-T28, still technically supported, just no longer
   !! how this module's own callers use it): a single plan created with
   !! nthreads=1 is safe to EXECUTE CONCURRENTLY from multiple caller
   !! threads via the "new-array execute" form (sfftw_execute_dft with
   !! explicit in/out arguments) as long as each concurrent call
   !! supplies its own distinct arrays -- true here, since image/
   !! image_out/the internal work arrays are all local to each call.
   !! Verified directly (see this module's own test suite, at the time):
   !! 16 OpenMP threads calling convolve_to_beam concurrently against
   !! the SAME shared plan, each on its own image/beam pair, matched a
   !! serial run of the same 16 calls exactly. Do NOT combine the two
   !! strategies -- a plan created with nthreads>1, executed
   !! concurrently by nthreads>1 caller threads, oversubscribes threads
   !! by a factor of nthreads (each concurrent execute call would try to
   !! internally fan out to nthreads threads of its own).
   use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
   implicit none
   private
   public :: plan_convolution, convolve_to_beam, destroy_convolution_plan,&
   &next_fast_fft_size

   ! T27 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): the FFT buffers
   ! (plan_convolution's own scratch, and convolve_to_beam's cimg/cmsk/
   ! g_final/c_d/c_m below) are single precision (sp) -- FITS flux data
   ! is already written to disk as BITPIX=-32 (see write_convolved_file/
   ! write_matched_file's own FTPHPR calls), so there is no real
   ! dynamic-range/precision argument for carrying it at double
   ! precision through the FFT, and halving these (the largest, padded-
   ! size) arrays roughly halves convolve_to_beam's own peak per-plane
   ! memory. The one place this does NOT apply: the kernel exponent
   ! itself (g_arg/dg_arg below, a DIFFERENCE of two potentially large
   ! terms before exp()) stays real(dp) throughout its own computation,
   ! only cast down to sp at the point g_final is actually stored --
   ! single precision there risks catastrophic cancellation in the
   ! subtraction, which halving a bulk data array never does. image/
   ! image_out (convolve_to_beam's own dummy arguments) deliberately
   ! stay real(dp) -- the precision reduction is entirely internal to
   ! this module, so no caller (match_cubes.f90/convolve_cubes.f90)
   ! needs to change at all.

   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: deg2rad = pi/180.0_dp
   ! sqrt(8*ln2): standard FWHM -> Gaussian-sigma conversion factor.
   real(dp), parameter :: fwhm2sigma = 2.0_dp*sqrt(2.0_dp*log(2.0_dp))

   ! Minimum acceptable valid-data weight fraction for a NaN-containing
   ! plane's convolved output pixel (see convolve_to_beam's NaN-aware
   ! path below). A pixel's output is kept only if at least this
   ! fraction of the KERNEL'S OWN WEIGHT (not a raw NaN/valid pixel
   ! count -- the kernel is a Gaussian, so a NaN near its centre counts
   ! far more than one in its tail) came from valid input pixels;
   ! otherwise the output pixel is re-NaN'd rather than reported as a
   ! number derived mostly from fabricated (zero-filled) data. At the
   ! chosen 0.5, this is equivalent to: reject if more of the kernel's
   ! weight fell on NaN pixels than on valid ones.
   real(dp), parameter :: conv_nan_reject_frac = 0.5_dp

   ! FFTW3 constants (from /usr/include/fftw3.f) -- declared directly
   ! rather than `include`d: fftw3.f is fixed-form Fortran 77 (same
   ! issue as AST_PAR in reproject_cubes.f90, see its own comment) and
   ! cannot be included into a free-form .f90 file directly.
   integer, parameter :: fftw_forward = -1
   integer, parameter :: fftw_backward = 1
   integer, parameter :: fftw_estimate = 64

contains

   subroutine plan_convolution(nx, ny, plan_fwd, plan_bwd, nx_pad, ny_pad, nthreads)
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
      !!
      !! nthreads (T28, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): the
      !! resulting plan EXECUTES using this many threads INTERNALLY for
      !! every subsequent sfftw_execute_dft call against it
      !! (sfftw_plan_with_nthreads, applies to plans created immediately
      !! after the call, per FFTW's own documented global/thread-local
      !! planning-time state). This is a DIFFERENT parallelisation
      !! strategy than this module's original one (see the module
      !! header comment, now historical): originally, ONE single-
      !! threaded plan was executed CONCURRENTLY by many caller threads,
      !! each on their own plane. Now, ONE thread calls
      !! convolve_to_beam serially, once per plane, and the plan itself
      !! fans out across nthreads threads for each transform. These are
      !! MUTUALLY EXCLUSIVE strategies -- a caller must not do both at
      !! once (nthreads>1 concurrent callers each executing an
      !! nthreads>1-internally-threaded plan would oversubscribe
      !! threads^2-fold). Pass nthreads=1 for the old concurrent-callers
      !! usage pattern (behaves identically to before T28).
      integer, intent(in) :: nx, ny, nthreads
      integer(kind=8), intent(out) :: plan_fwd, plan_bwd
      integer, intent(out) :: nx_pad, ny_pad

      complex(sp), allocatable :: scratch(:,:)
      integer :: threads_ok

      nx_pad = next_fast_fft_size(nx)
      ny_pad = next_fast_fft_size(ny)

      ! sfftw_init_threads is safe to call more than once (FFTW's own
      ! documented behaviour -- subsequent calls are a no-op) so no
      ! separate "call once globally" bookkeeping is needed here; this
      ! is already called once per file (plan_convolution's own
      ! existing per-file granularity), not per plane.
      call sfftw_init_threads(threads_ok)
      if (threads_ok.eq.0) then
         write(*,'(A)') 'WARNING: sfftw_init_threads failed -- falling'//&
         &' back to single-threaded FFT execution (nthreads=1).'
         call sfftw_plan_with_nthreads(1)
      else
         call sfftw_plan_with_nthreads(max(1, nthreads))
      endif

      ! FFTW_ESTIMATE plans don't depend on the array CONTENTS, or even
      ! specifically on the memory used here, only the shape -- every
      ! actual convolve_to_beam call below uses the "new-array execute"
      ! form with its own arrays instead. FFTW's own documentation
      ! describes this as fully supported (correct regardless of
      ! alignment, though a non-ESTIMATE plan might not be as fast on a
      ! differently-aligned array than the one it was planned with --
      ! irrelevant for ESTIMATE, which never does alignment-specific
      ! optimisation in the first place). sfftw_* (single precision,
      ! T27) -- the plan's own precision must match whatever array it
      ! will later execute against (convolve_to_beam's cimg/cmsk, both
      ! complex(sp)); a plan created via dfftw_* cannot be executed
      ! against a complex(sp) array (undefined behaviour, not a clean
      ! error) so this scratch array is complex(sp) too, matching.
      allocate(scratch(nx_pad, ny_pad))
      call sfftw_plan_dft_2d(plan_fwd, nx_pad, ny_pad, scratch, scratch, fftw_forward, fftw_estimate)
      call sfftw_plan_dft_2d(plan_bwd, nx_pad, ny_pad, scratch, scratch, fftw_backward, fftw_estimate)
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

      call sfftw_destroy_plan(plan_fwd)
      call sfftw_destroy_plan(plan_bwd)
   end subroutine destroy_convolution_plan

   subroutine convolve_to_beam(plan_fwd, plan_bwd, image, nx, ny, nx_pad,&
   &ny_pad, dx, dy, bmaj_in, bmin_in, bpa_in, bmaj, bmin, bpa, image_out,&
   &status)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value,&
      &ieee_quiet_nan
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
      !!
      !! NaN handling: if image has no NaN pixels, the computation below
      !! is byte-for-byte identical to the original NaN-agnostic version
      !! (zero regression risk for the common case). If image DOES
      !! contain NaN (real ASKAP data routinely does, from primary-beam
      !! edge blanking), a plain FFT convolution would let a single NaN
      !! poison the entire plane -- the FFT is a global operation, so
      !! there is no such thing as a "locally" NaN output pixel.
      !! Instead: NaN pixels are zero-filled and convolved to give C_D;
      !! separately, a 0/1 validity mask is convolved through the SAME
      !! kernel and normalised by g_ratio (the kernel's own DC gain) to
      !! give C_M, the fraction of each output pixel's kernel weight
      !! that came from valid input (0..1, exactly 1 with no NaN
      !! nearby). image_out = C_D/C_M where C_M >= conv_nan_reject_frac,
      !! and NaN otherwise -- i.e. an output pixel is only reported if
      !! most of the kernel's weight, not just most of the surrounding
      !! pixel COUNT, came from real data. Matches the precedent set by
      !! RACS-tools' own convolve() (verified directly against its
      !! source), except RACS-tools does not renormalise the
      !! partially-contaminated case at all -- this implementation is
      !! more rigorous for the pixels it keeps rather than only
      !! blanket-rejecting fully-contaminated ones.
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
      real(dp) :: nanval
      real(dp), allocatable :: u(:), v(:)
      real(sp), allocatable :: c_d(:,:), c_m(:,:)
      logical, allocatable :: nan_mask(:,:)
      logical :: has_nan
      complex(sp), allocatable :: cimg(:,:), cmsk(:,:), g_final(:,:)
      integer :: ix, iy

      status = 0

      ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): convolve_to_beam
      ! is now called SERIALLY, once per plane, by its caller (the
      ! caller's own outer per-plane loop is no longer itself an
      ! !$omp parallel do -- see T28) -- so every elementwise step below
      ! is parallelised internally with its own !$omp construct, or it
      ! would run single-threaded and give up the throughput the old
      ! per-plane-per-thread design got for free. The FFT itself is
      ! parallelised separately, via plan_convolution's own
      ! sfftw_plan_with_nthreads (see there).
      allocate(nan_mask(nx, ny))
      !$omp workshare
      nan_mask = ieee_is_nan(image)
      !$omp end workshare
      has_nan = any(nan_mask)

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
      ! T28: private(ur,vr,g_arg,ur_in,vr_in,dg_arg) -- each thread's own
      ! per-(ix,iy) scratch, not shared state; u/v/cos_bpa/etc. stay
      ! shared (read-only within this loop).
      !$omp parallel do collapse(2) default(shared)&
      !$omp& private(ix,iy,ur,vr,g_arg,ur_in,vr_in,dg_arg)
      do iy = 1, ny_pad
         do ix = 1, nx_pad
            ur = u(ix)*cos_bpa - v(iy)*sin_bpa
            vr = u(ix)*sin_bpa + v(iy)*cos_bpa
            g_arg = -2.0_dp*pi**2 * ((sx*ur)**2 + (sy*vr)**2)

            ur_in = u(ix)*cos_bpa_in - v(iy)*sin_bpa_in
            vr_in = u(ix)*sin_bpa_in + v(iy)*cos_bpa_in
            dg_arg = -2.0_dp*pi**2 * ((sx_in*ur_in)**2 + (sy_in*vr_in)**2)

            ! Exponent computed and combined entirely in dp above (g_arg
            ! - dg_arg is a difference of two potentially large terms --
            ! single precision here would risk real cancellation error);
            ! only the finished value is cast down to sp (T27) at the
            ! point of storage.
            g_final(ix,iy) = cmplx(g_ratio * exp(cmplx(g_arg - dg_arg, 0.0_dp, dp)), kind=sp)
         enddo
      enddo
      !$omp end parallel do
      deallocate(u, v)

      allocate(cimg(nx_pad, ny_pad))
      !$omp workshare
      cimg = cmplx(0.0_sp, 0.0_sp, sp)
      !$omp end workshare
      if (has_nan) then
         !$omp workshare
         cimg(1:nx, 1:ny) = cmplx(merge(0.0_dp, image, nan_mask), 0.0_dp, sp)
         !$omp end workshare
      else
         !$omp workshare
         cimg(1:nx, 1:ny) = cmplx(image, 0.0_dp, sp)
         !$omp end workshare
      endif
      call sfftw_execute_dft(plan_fwd, cimg, cimg)

      !$omp workshare
      cimg = cimg*g_final
      !$omp end workshare

      ! FFTW's transforms are unnormalised (forward then backward scales
      ! the result by nx_pad*ny_pad, same convention as numpy.fft.fft2/
      ! ifft2 -- numpy just applies the 1/N inside ifft2 for you; FFTW
      ! leaves it to the caller, matching its own documented convention).
      ! Cropped back from the same top-left corner the image was placed
      ! at above. Division itself stays in dp (nx_pad*ny_pad cast to dp
      ! before dividing) even though cimg's own bits are only sp-precise
      ! -- preserves the exact same normalisation pathway as before this
      ! module went single precision (T27), rather than also changing
      ! how the division itself is computed.
      if (.not. has_nan) then
         deallocate(g_final)
         call sfftw_execute_dft(plan_bwd, cimg, cimg)
         !$omp workshare
         image_out = real(cimg(1:nx, 1:ny), dp) / real(nx_pad*ny_pad, dp)
         !$omp end workshare
         deallocate(cimg)
      else
         call sfftw_execute_dft(plan_bwd, cimg, cimg)
         allocate(c_d(nx, ny))
         !$omp workshare
         c_d = real(cimg(1:nx, 1:ny), dp) / real(nx_pad*ny_pad, dp)
         !$omp end workshare
         deallocate(cimg)

         ! Convolve the 0/1 validity mask through the SAME kernel as the
         ! data -- the fraction of each output pixel's kernel weight
         ! that came from valid input. Raw mask convolution maxes out at
         ! g_ratio (the kernel's own DC gain), not 1, when the mask is
         ! entirely 1 nearby -- see this subroutine's own header comment
         ! for why g_ratio isn't unit gain -- so divide it out once here
         ! to get a clean 0..1 fraction.
         allocate(cmsk(nx_pad, ny_pad))
         !$omp workshare
         cmsk = cmplx(0.0_sp, 0.0_sp, sp)
         cmsk(1:nx, 1:ny) = cmplx(merge(0.0_dp, 1.0_dp, nan_mask), 0.0_dp, sp)
         !$omp end workshare
         call sfftw_execute_dft(plan_fwd, cmsk, cmsk)
         !$omp workshare
         cmsk = cmsk*g_final
         !$omp end workshare
         deallocate(g_final)
         call sfftw_execute_dft(plan_bwd, cmsk, cmsk)

         allocate(c_m(nx, ny))
         !$omp workshare
         c_m = real(cmsk(1:nx, 1:ny), dp) / real(nx_pad*ny_pad, dp) / g_ratio
         !$omp end workshare
         deallocate(cmsk)

         nanval = ieee_value(1.0_dp, ieee_quiet_nan)
         !$omp workshare
         where (c_m >= conv_nan_reject_frac)
            image_out = c_d/c_m
         elsewhere
            image_out = nanval
         endwhere
         !$omp end workshare
         deallocate(c_d, c_m)
      endif
      deallocate(nan_mask)
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
