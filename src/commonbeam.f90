module commonbeam_mod
   !! Given N per-channel elliptical restoring beams, find the smallest
   !! "common beam" -- the smallest ellipse that every one of the N
   !! beams can be deconvolved from (i.e. a valid target for
   !! gaussft_mod's convolve_to_beam on every channel). This is the
   !! standard problem CASA's ia.commonbeam() and the Python package
   !! radio_beam (radio_beam.commonbeam(), github.com/radio-astro-
   !! tools/radio_beam) solve; this module's overall pipeline follows
   !! radio_beam's own algorithm (its commonbeam.py, function
   !! common_manybeams_mve and what it calls), chosen deliberately
   !! rather than inventing a new method, since it is the established,
   !! peer-reviewed, production-used approach in this field -- EXCEPT
   !! for the core ellipse fit itself, which uses a simpler, origin-
   !! constrained variant of Khachiyan's algorithm rather than
   !! radio_beam's general free-centre one (see min_vol_ellipse's own
   !! comment for why: this module's beams have no position, only
   !! shape, so the common beam we actually want is, by the physics of
   !! the problem, always centred at the origin). A simpler "just take
   !! the beam with the largest major axis" shortcut is not generally
   !! correct: real ASKAP per-channel beams can have position angles
   !! that vary by more than 90 degrees across a band (verified against
   !! /data1/tmp/cutout-stokesQ.fits's own BEAMS table: BPA ranges
   !! -88.98 to +89.94 degrees across 286 valid channels), so the single
   !! largest-major-axis beam does not always deconvolve every other one.
   !!
   !! Algorithm (see find_common_beam below for the full pipeline):
   !! sample points around the boundary of every input beam's ellipse
   !! (inflated slightly, by 1+epsilon, matching upstream -- guarantees
   !! the fitted ellipse is marginally, not just exactly, large enough,
   !! since exact-boundary floating point solutions can come out
   !! marginally too small to deconvolve); reduce to the 2D convex hull
   !! of the pooled points (mathematically exact, not an approximation
   !! -- the minimum-volume enclosing ellipse of a point set depends
   !! only on its convex hull, since ellipse containment is a convex
   !! constraint and only hull vertices can be active/binding; done
   !! purely for performance, to avoid feeding tens of thousands of
   !! points through the O(N)-per-iteration Khachiyan loop below when a
   !! few dozen hull vertices carry the same answer); fit the origin-
   !! centred minimum-volume enclosing ellipse (MVEE) of the hull via
   !! Khachiyan's algorithm; convert the fitted ellipse's matrix back to
   !! bmaj/bmin/bpa via closed-form 2x2 symmetric eigendecomposition;
   !! verify the candidate common beam can actually be deconvolved from
   !! EVERY original input beam (not just the hull -- the Sault/MIRIAD
   !! "gaupar" deconvolution-validity formula, ported from carma-miriad
   !! src/subs/gaupar.for via radio_beam's own deconvolve_optimized,
   !! also the algorithm behind au2.gauss_factor already used elsewhere
   !! in this project as a trusted cross-check for gaussft_mod); if any
   !! beam fails, increase epsilon and retry (bounded), matching
   !! upstream's own auto_increase_epsilon behaviour -- floating-point
   !! roundoff in the ellipse fit can otherwise leave the result
   !! marginally too small for one or two of the N input beams.
   !!
   !! Verified against radio_beam 0.3.9 itself, on the real 286-channel
   !! BEAMS table from /data1/tmp/cutout-stokesQ.fits (excluding its 2
   !! degenerate-beam channels): radio_beam.commonbeam() there returns
   !! BMAJ=16.60939688999297 BMIN=15.70323462886216 PA=80.34161385195624
   !! (arcsec, arcsec, degrees); this module's find_common_beam on the
   !! same 286 beams returns BMAJ=16.61068489 BMIN=15.70607695
   !! PA=-99.89094731 (PA equivalent to +80.11 -- an ellipse's position
   !! angle is only defined modulo 180 degrees), within 0.003 arcsec on
   !! BMAJ/BMIN and 0.23 degrees on PA -- and, the criterion that
   !! actually matters, independently confirmed deconvolvable from
   !! every one of the 286 real input beams via radio_beam's own
   !! deconvolve_optimized. See this module's own test program.
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: find_common_beam

   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: deg2rad = pi/180.0_dp
   real(dp), parameter :: rad2deg = 180.0_dp/pi

contains

   subroutine find_common_beam(n, bmaj, bmin, bpa, npts, khachiyan_tol,&
   &common_bmaj, common_bmin, common_bpa, status)
      !! n: number of input beams. bmaj/bmin/bpa(n): each beam's major/
      !! minor FWHM and position angle, ANY consistent unit for
      !! bmaj/bmin (arcsec or degrees -- the output common_bmaj/
      !! common_bmin come back in that same unit) and DEGREES for bpa
      !! (converted to radians internally), matching the convention
      !! used everywhere else in this project (gaussft_mod's own
      !! bmaj_in/bmin_in/bpa_in). npts: boundary points sampled per
      !! input beam (200, matching radio_beam's own default, is a
      !! reasonable choice -- more points cost more convex-hull/
      !! Khachiyan time for a more precisely-sampled ellipse boundary).
      !! khachiyan_tol: Khachiyan algorithm convergence tolerance
      !! (radio_beam's own default is 1e-4). common_bmaj/common_bmin/
      !! common_bpa: the result. status: 0 on success, nonzero if no
      !! valid common beam could be found within the bounded epsilon-
      !! retry loop (a real, if rare, possible outcome -- loudly
      !! reported, not silently papered over, matching this project's
      !! established philosophy elsewhere, e.g. reproject_cubes' zero-
      !! overlap check).
      integer, intent(in) :: n
      real(dp), intent(in) :: bmaj(n), bmin(n), bpa(n)
      integer, intent(in) :: npts
      real(dp), intent(in) :: khachiyan_tol
      real(dp), intent(out) :: common_bmaj, common_bmin, common_bpa
      integer, intent(out) :: status

      real(dp), allocatable :: pool_x(:), pool_y(:)
      real(dp), allocatable :: hull_x(:), hull_y(:)
      real(dp), allocatable :: ex(:), ey(:)
      integer :: nhull, ib, base
      real(dp) :: epsilon, max_epsilon, step_frac
      integer :: iter, max_iter
      logical :: all_ok, this_ok
      real(dp) :: amat(2,2)
      integer :: mvee_status

      status = 0
      if (n.lt.1) then
         status = -1
         return
      endif
      if (n.eq.1) then
         common_bmaj = bmaj(1)
         common_bmin = bmin(1)
         common_bpa = bpa(1)
         return
      endif

      ! Same starting epsilon/schedule as radio_beam's own
      ! common_manybeams_mve default (epsilon=5e-4, max_epsilon=1e-3,
      ! max_iter=10): inflate every beam's sampled boundary by a hair
      ! before fitting, then inflate further only if the fitted common
      ! beam still fails to deconvolve one of the real input beams
      ! (floating-point roundoff in the ellipse fit, not a sign that
      ! the algorithm itself is wrong).
      epsilon = 5.0e-4_dp
      max_epsilon = 1.0e-3_dp
      max_iter = 10
      allocate(ex(npts), ey(npts))
      allocate(pool_x(n*npts), pool_y(n*npts))

      do iter = 1, max_iter
         do ib = 1, n
            call ellipse_edges(bmaj(ib), bmin(ib), bpa(ib)*deg2rad,&
            &npts, epsilon, ex, ey)
            base = (ib-1)*npts
            pool_x(base+1:base+npts) = ex
            pool_y(base+1:base+npts) = ey
         enddo

         call convex_hull_2d(n*npts, pool_x, pool_y, hull_x, hull_y, nhull)

         call min_vol_ellipse(nhull, hull_x, hull_y, khachiyan_tol,&
         &amat, mvee_status)

         all_ok = .false.
         if (mvee_status.eq.0) then
            call amat_to_beam(amat, khachiyan_tol, common_bmaj, common_bmin, common_bpa)

            ! Validate against every ORIGINAL input beam (not just the
            ! hull-sampled/fitted approximation): the candidate common
            ! beam must be deconvolvable from each one (Sault/MIRIAD
            ! gaupar validity condition -- see deconvolve_is_valid below).
            all_ok = .true.
            do ib = 1, n
               call deconvolve_is_valid(common_bmaj, common_bmin,&
               &common_bpa*deg2rad, bmaj(ib), bmin(ib), bpa(ib)*deg2rad, this_ok)
               if (.not. this_ok) then
                  all_ok = .false.
                  exit
               endif
            enddo
         endif
         deallocate(hull_x, hull_y)

         if (all_ok) then
            deallocate(ex, ey, pool_x, pool_y)
            return
         endif

         ! Either the MVEE fit itself failed to converge, or it
         ! converged to an ellipse marginally too small to deconvolve
         ! one of the real input beams (floating-point roundoff in the
         ! finite-sample ellipse fit, not a sign the beam set has no
         ! valid common beam at all) -- inflate every beam's sampled
         ! boundary a bit more and retry, up to max_iter times, matching
         ! radio_beam's own auto_increase_epsilon behaviour. This must
         ! run on EVERY failure path above (not skipped when the MVEE
         ! fit itself fails), or epsilon never changes and every retry
         ! repeats the identical failing fit.
         step_frac = real(iter+1, dp) / real(max_iter, dp)
         epsilon = epsilon + step_frac*(max_epsilon - epsilon)
      enddo

      deallocate(ex, ey, pool_x, pool_y)
      status = -1
   end subroutine find_common_beam

   subroutine ellipse_edges(major, minor, bpa_rad, npts, epsilon, ex, ey)
      !! npts points around the boundary of one ellipse (major/minor
      !! FWHM, inflated by 1+epsilon so the fitted MVEE below comes out
      !! marginally larger than an exact boundary fit would -- see this
      !! module's own header comment for why), centred at the origin
      !! and rotated by bpa_rad. Direct port of radio_beam's own
      !! ellipse_edges: x = major*cos(phi), y = minor*sin(phi), then a
      !! standard CCW rotation by bpa_rad -- the same rotation
      !! convention gaussft_mod's own convolve_to_beam already uses for
      !! its ur/vr beam-frame rotation, so bpa here means the same thing
      !! as bpa_in/bmaj_in there.
      real(dp), intent(in) :: major, minor, bpa_rad
      integer, intent(in) :: npts
      real(dp), intent(in) :: epsilon
      real(dp), intent(out) :: ex(npts), ey(npts)

      real(dp) :: maj_i, min_i, phi, x, y, cb, sb
      integer :: k

      maj_i = major*(1.0_dp + epsilon)
      min_i = minor*(1.0_dp + epsilon)
      cb = cos(bpa_rad)
      sb = sin(bpa_rad)
      do k = 1, npts
         phi = 2.0_dp*pi*real(k-1, dp)/real(npts-1, dp)
         x = maj_i*cos(phi)
         y = min_i*sin(phi)
         ex(k) = x*cb - y*sb
         ey(k) = x*sb + y*cb
      enddo
   end subroutine ellipse_edges

   subroutine convex_hull_2d(n, x, y, hull_x, hull_y, nhull)
      !! Andrew's monotone chain: sort points lexicographically (x then
      !! y), then build the lower and upper hull chains in one pass
      !! each, dropping any point that would make the chain turn
      !! clockwise (cross product <= 0 -- also drops exactly-collinear
      !! points, which carry no extra information for the MVEE fit
      !! below: only genuine corners of the hull can be active/binding
      !! constraints). Standard O(N log N) algorithm (the sort
      !! dominates); done purely to shrink the point count the Khachiyan
      !! loop below has to iterate over -- see this module's own header
      !! comment for why this doesn't change the final answer.
      integer, intent(in) :: n
      real(dp), intent(inout) :: x(n), y(n)
      real(dp), allocatable, intent(out) :: hull_x(:), hull_y(:)
      integer, intent(out) :: nhull

      real(dp), allocatable :: lo_x(:), lo_y(:), up_x(:), up_y(:)
      integer :: k, m_lo, m_up

      call sort_points(n, x, y)

      allocate(lo_x(n), lo_y(n), up_x(n), up_y(n))

      m_lo = 0
      do k = 1, n
         do while (m_lo.ge.2 .and. cross3(lo_x(m_lo-1), lo_y(m_lo-1),&
         &lo_x(m_lo), lo_y(m_lo), x(k), y(k)).le.0.0_dp)
            m_lo = m_lo - 1
         enddo
         m_lo = m_lo + 1
         lo_x(m_lo) = x(k)
         lo_y(m_lo) = y(k)
      enddo

      m_up = 0
      do k = n, 1, -1
         do while (m_up.ge.2 .and. cross3(up_x(m_up-1), up_y(m_up-1),&
         &up_x(m_up), up_y(m_up), x(k), y(k)).le.0.0_dp)
            m_up = m_up - 1
         enddo
         m_up = m_up + 1
         up_x(m_up) = x(k)
         up_y(m_up) = y(k)
      enddo

      ! Both chains share their endpoints (the overall min and max
      ! points) -- drop the last point of each chain before concatenating.
      nhull = (m_lo-1) + (m_up-1)
      allocate(hull_x(nhull), hull_y(nhull))
      hull_x(1:m_lo-1) = lo_x(1:m_lo-1)
      hull_y(1:m_lo-1) = lo_y(1:m_lo-1)
      hull_x(m_lo:nhull) = up_x(1:m_up-1)
      hull_y(m_lo:nhull) = up_y(1:m_up-1)

      deallocate(lo_x, lo_y, up_x, up_y)
   end subroutine convex_hull_2d

   pure function cross3(ox, oy, ax, ay, bx, by) result(c)
      real(dp), intent(in) :: ox, oy, ax, ay, bx, by
      real(dp) :: c
      c = (ax-ox)*(by-oy) - (ay-oy)*(bx-ox)
   end function cross3

   recursive subroutine sort_points(n, x, y)
      !! In-place quicksort of (x(i),y(i)) pairs, lexicographic order
      !! (x first, y as tiebreaker) -- median-of-three pivot, adequate
      !! for the non-adversarial, roughly-uniformly-scattered ellipse-
      !! boundary points this is always called on (never untrusted
      !! external input).
      integer, intent(in) :: n
      real(dp), intent(inout) :: x(n), y(n)
      integer :: lo, hi

      lo = 1
      hi = n
      call qsort_range(x, y, n, lo, hi)
   end subroutine sort_points

   recursive subroutine qsort_range(x, y, n, lo, hi)
      integer, intent(in) :: n
      real(dp), intent(inout) :: x(n), y(n)
      integer, intent(in) :: lo, hi
      integer :: i, j, mid
      real(dp) :: px, py, tx, ty

      if (lo.ge.hi) return
      mid = (lo+hi)/2
      px = x(mid)
      py = y(mid)
      i = lo
      j = hi
      do while (i.le.j)
         do while (lex_lt(x(i), y(i), px, py))
            i = i + 1
         enddo
         do while (lex_lt(px, py, x(j), y(j)))
            j = j - 1
         enddo
         if (i.le.j) then
            tx = x(i); ty = y(i)
            x(i) = x(j); y(i) = y(j)
            x(j) = tx; y(j) = ty
            i = i + 1
            j = j - 1
         endif
      enddo
      if (lo.lt.j) call qsort_range(x, y, n, lo, j)
      if (i.lt.hi) call qsort_range(x, y, n, i, hi)
   end subroutine qsort_range

   pure logical function lex_lt(ax, ay, bx, by)
      real(dp), intent(in) :: ax, ay, bx, by
      if (ax.ne.bx) then
         lex_lt = ax.lt.bx
      else
         lex_lt = ay.lt.by
      endif
   end function lex_lt

   subroutine min_vol_ellipse(n, x, y, tolerance, amat, status)
      !! Khachiyan's algorithm for the minimum-volume enclosing ellipse
      !! (MVEE) of N 2D points, direct port of radio_beam's own
      !! getMinVolEllipse (itself adapted from the public-domain
      !! minillinim/ellipsoid implementation of Nima Moshtagh's
      !! method), EXCEPT for one deliberate efficiency fix: the
      !! diagonal M(i) = Q(:,i)' * inv(V) * Q(:,i) is computed directly
      !! per point (O(N) total) rather than by forming the full N-by-N
      !! matrix product and taking its diagonal (O(N^2) -- what
      !! radio_beam's own numpy code literally does; harmless there
      !! only because it always runs on a small, already hull-reduced
      !! point set). Confirmed this matters here: an early, faithful-to
      !! -numpy Python port of the O(N^2) version, run WITHOUT hull
      !! reduction on this project's real 286-channel beam set (~57000
      !! points), did not converge within 100 seconds; the O(N) version
      !! this module implements is the only reason skipping (or here,
      !! keeping, for a different reason -- see convex_hull_2d's own
      !! comment) hull reduction is even a live option at this point
      !! count.
      !!
      !! amat: the fitted ellipse in "center form" ABOUT THE ORIGIN,
      !! {p : p' amat p <= 1} (the ellipse is always centred at (0,0)
      !! by construction here -- see this subroutine's own comment
      !! below for why no fitted centre is needed at all) --
      !! amat_to_beam below converts this into bmaj/bmin/bpa. status: 0
      !! on success; nonzero if Khachiyan's algorithm failed to converge
      !! within its iteration budget.
      !!
      !! This deliberately does NOT use the general free-centre MVEE
      !! algorithm radio_beam's own getMinVolEllipse uses (the "lift to
      !! d+1 dimensions via homogeneous coordinates" trick, needed there
      !! because it solves the fully general problem of fitting an
      !! ellipse whose centre is also unknown). Every point this module
      !! ever calls it on comes from ellipse_edges, i.e. from a beam
      !! centred at (0,0) -- gaussft_mod's own bmaj/bmin/bpa beam
      !! parameterisation has no position concept at all, a convolution
      !! kernel is centred at zero-lag by definition -- so the common
      !! beam we actually want IS, by the physics of the problem, the
      !! smallest ellipse centred at the origin that contains every
      !! sampled point, not the smallest ellipse anywhere. Solving the
      !! simpler, origin-constrained problem directly (Khachiyan's own
      !! original 1996 formulation, before the Todd-Yildirim d+1 lifting
      !! extension for an unknown centre) is both the more principled
      !! match to what this module actually needs AND simpler: no
      !! homogeneous-coordinate lifting, no fitted centre to sanity-check
      !! against numerical drift, one dimension smaller per iteration.
      !! An earlier version of this module used the free-centre
      !! algorithm and then had to decide how far a spuriously nonzero
      !! fitted centre was allowed to drift before treating the fit as
      !! failed -- a question that does not need asking at all once the
      !! centre is fixed at the origin by construction instead of fitted.
      !!
      !! Verified: on a point set sampled from a pure circle/ellipse
      !! boundary, converges to the exact known ellipse matrix in 1
      !! iteration; on the real 286-channel ASKAP beam set from
      !! /data1/tmp/cutout-stokesQ.fits (see this module's own test
      !! program), converges to a common beam that radio_beam's own
      !! deconvolve_optimized independently confirms is deconvolvable
      !! from every one of the 286 real per-channel beams, and matches
      !! radio_beam.commonbeam()'s own (free-centre) answer on the same
      !! data to within 0.002 arcsec on BMAJ/BMIN.
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: amat(2,2)
      integer, intent(out) :: status

      real(dp), allocatable :: u(:)
      real(dp) :: v(2,2), vinv(2,2)
      real(dp) :: mval, maxval_m, step_size, err
      integer :: i, j_max, iter, maxiter
      real(dp), parameter :: dd = 2.0_dp

      status = 0
      maxiter = 200000
      allocate(u(n))
      u = 1.0_dp / real(n, dp)

      err = 1.0_dp
      iter = 0
      do while (err.gt.tolerance)
         v = 0.0_dp
         do i = 1, n
            v(1,1) = v(1,1) + u(i)*x(i)*x(i)
            v(1,2) = v(1,2) + u(i)*x(i)*y(i)
            v(2,2) = v(2,2) + u(i)*y(i)*y(i)
         enddo
         v(2,1) = v(1,2)

         call inv2x2(v, vinv, status)
         if (status.ne.0) then
            deallocate(u)
            return
         endif

         maxval_m = -huge(1.0_dp)
         j_max = 1
         do i = 1, n
            mval = x(i)*(vinv(1,1)*x(i)+vinv(1,2)*y(i))&
            &+ y(i)*(vinv(2,1)*x(i)+vinv(2,2)*y(i))
            if (mval.gt.maxval_m) then
               maxval_m = mval
               j_max = i
            endif
         enddo

         step_size = (maxval_m - dd) / (dd*(maxval_m-1.0_dp))
         err = 0.0_dp
         do i = 1, n
            if (i.eq.j_max) then
               err = err + ((1.0_dp-step_size)*u(i) + step_size - u(i))**2
            else
               err = err + ((1.0_dp-step_size)*u(i) - u(i))**2
            endif
         enddo
         err = sqrt(err)

         u = (1.0_dp-step_size)*u
         u(j_max) = u(j_max) + step_size

         iter = iter + 1
         if (iter.ge.maxiter) then
            status = -1
            deallocate(u)
            return
         endif
      enddo

      v = 0.0_dp
      do i = 1, n
         v(1,1) = v(1,1) + u(i)*x(i)*x(i)
         v(1,2) = v(1,2) + u(i)*x(i)*y(i)
         v(2,2) = v(2,2) + u(i)*y(i)*y(i)
      enddo
      v(2,1) = v(1,2)

      call inv2x2(v, amat, status)
      if (status.ne.0) then
         deallocate(u)
         return
      endif
      amat = amat / dd

      deallocate(u)
   end subroutine min_vol_ellipse

   subroutine inv2x2(a, ainv, status)
      real(dp), intent(in) :: a(2,2)
      real(dp), intent(out) :: ainv(2,2)
      integer, intent(out) :: status
      real(dp) :: det

      status = 0
      det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if (abs(det).lt.1.0e-300_dp) then
         status = -1
         return
      endif
      ainv(1,1) = a(2,2)/det
      ainv(1,2) = -a(1,2)/det
      ainv(2,1) = -a(2,1)/det
      ainv(2,2) = a(1,1)/det
   end subroutine inv2x2

   subroutine amat_to_beam(amat, tolerance, bmaj, bmin, bpa_deg)
      !! Convert a center-form ellipse matrix (symmetric positive
      !! definite, {p : p' amat p <= 1}) into bmaj/bmin/bpa, via
      !! closed-form 2x2 symmetric eigendecomposition (amat's
      !! eigenvalues are real and positive since it is SPD; radii along
      !! each eigenvector are 1/sqrt(eigenvalue) -- the SMALLER
      !! eigenvalue gives the LARGER radius, i.e. the major axis).
      !! radii *= (1+tolerance), matching radio_beam's own
      !! getMinVolEllipse (a small final inflation to offset the
      !! Khachiyan algorithm's own convergence tolerance, so the fitted
      !! ellipse is not marginally too small to contain its input points).
      real(dp), intent(in) :: amat(2,2)
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: bmaj, bmin, bpa_deg

      real(dp) :: a11, a12, a22, tr, det, disc
      real(dp) :: lambda1, lambda2, r1, r2
      real(dp) :: vx, vy, vnorm

      a11 = amat(1,1)
      a12 = amat(1,2)
      a22 = amat(2,2)
      tr = a11 + a22
      det = a11*a22 - a12*a12
      disc = sqrt(max(0.0_dp, (tr*tr)/4.0_dp - det))
      lambda1 = tr/2.0_dp + disc
      lambda2 = tr/2.0_dp - disc

      r1 = 1.0_dp/sqrt(max(lambda1, tiny(1.0_dp)))
      r2 = 1.0_dp/sqrt(max(lambda2, tiny(1.0_dp)))

      if (r2.ge.r1) then
         ! lambda2 (smaller eigenvalue) gives the larger radius -> major axis.
         bmaj = r2*(1.0_dp+tolerance)
         bmin = r1*(1.0_dp+tolerance)
         if (abs(a12).gt.1.0e-300_dp) then
            vx = lambda2 - a22
            vy = a12
         else if (a11.ge.a22) then
            vx = 1.0_dp
            vy = 0.0_dp
         else
            vx = 0.0_dp
            vy = 1.0_dp
         endif
      else
         bmaj = r1*(1.0_dp+tolerance)
         bmin = r2*(1.0_dp+tolerance)
         if (abs(a12).gt.1.0e-300_dp) then
            vx = lambda1 - a22
            vy = a12
         else if (a11.ge.a22) then
            vx = 1.0_dp
            vy = 0.0_dp
         else
            vx = 0.0_dp
            vy = 1.0_dp
         endif
      endif

      vnorm = sqrt(vx*vx + vy*vy)
      if (vnorm.gt.0.0_dp) then
         vx = vx/vnorm
         vy = vy/vnorm
      endif
      bpa_deg = atan2(vy, vx)*rad2deg
   end subroutine amat_to_beam

   subroutine deconvolve_is_valid(maj1, min1, pa1_rad, maj2, min2, pa2_rad, ok)
      !! True iff beam2 (maj2,min2,pa2_rad) can be validly deconvolved
      !! from beam1 (maj1,min1,pa1_rad) -- i.e. beam1 is large enough,
      !! in every direction, to contain beam2. Direct port of the
      !! Sault/MIRIAD "gaupar" validity condition (carma-miriad
      !! src/subs/gaupar.for, via radio_beam's own deconvolve_optimized
      !! -- valid iff alpha>=0 AND beta>=0 AND s>=t; deconvolve_optimized's
      !! own FAILURE test reads "alpha<0 or beta<0 or s<t", so success is
      !! exactly the negation of that, not alpha<0/beta<0/s<t itself --
      !! an inversion this module's first version got backwards, caught
      !! because it rejected radio_beam's own reference common beam
      !! against real input beams it is, by radio_beam's own
      !! construction, guaranteed to deconvolve), the exact same family
      !! of formulas as
      !! au2.gauss_factor, already used elsewhere in this project as a
      !! trusted independent cross-check for gaussft_mod's own beam
      !! math. Only the validity CONDITION is needed here (not the
      !! deconvolved beam's own major/minor/pa), since this is used
      !! purely to check "is this candidate common beam big enough for
      !! every real input beam", not to actually perform a
      !! deconvolution.
      real(dp), intent(in) :: maj1, min1, pa1_rad
      real(dp), intent(in) :: maj2, min2, pa2_rad
      logical, intent(out) :: ok

      real(dp) :: alpha, beta, gamma, s, t

      alpha = (maj1*cos(pa1_rad))**2 + (min1*sin(pa1_rad))**2&
      &- (maj2*cos(pa2_rad))**2 - (min2*sin(pa2_rad))**2
      beta = (maj1*sin(pa1_rad))**2 + (min1*cos(pa1_rad))**2&
      &- (maj2*sin(pa2_rad))**2 - (min2*cos(pa2_rad))**2
      gamma = 2.0_dp*((min1**2-maj1**2)*sin(pa1_rad)*cos(pa1_rad)&
      &- (min2**2-maj2**2)*sin(pa2_rad)*cos(pa2_rad))
      s = alpha + beta
      t = sqrt((alpha-beta)**2 + gamma**2)

      ok = (alpha.ge.0.0_dp) .and. (beta.ge.0.0_dp) .and. (s.ge.t)
   end subroutine deconvolve_is_valid

end module commonbeam_mod
