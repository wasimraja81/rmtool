module exact_overlap_mod
   !! Flux-conserving reprojection by spherical overlap area
   !! (interp=exact). Every pixel, input and output, is a spherical
   !! quadrilateral: its four corners mapped to the output's sky frame,
   !! joined by great-circle arcs. An output pixel's value is the
   !! overlap-area-weighted mean of the input pixels it overlaps, so each
   !! input pixel's area is shared out among the output pixels it
   !! overlaps and the sky integral of a surface brightness (Jy/beam) is
   !! carried over. Design, measured accuracy and the reasons for each
   !! safeguard below: docs/dev/REPROJECT_EXACT_PLAN.md.
   !!
   !! The weights depend only on the two grids, so they are built once per
   !! input file and applied to every plane. exact_build_rows calls AST
   !! and must be called by the thread that owns the Mappings it is given;
   !! each call fills only its own output rows, so threads can share the
   !! rows of one table. exact_apply is plain Fortran.
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   use, intrinsic :: ieee_arithmetic
   implicit none
   private
   public :: exact_weights_t, exact_init, exact_build_rows, exact_apply,&
   &exact_apply_rows, exact_free, exact_weights_bytes, exact_estimate_bytes, mib_str

   type :: exact_row_t
      ! Output pixel i of this row uses entries ptr(i) .. ptr(i+1)-1 of
      ! idx/w; an empty range means nothing overlaps it (output NaN).
      integer, allocatable :: ptr(:)
      integer(int64), allocatable :: idx(:)
      real(dp), allocatable :: w(:)
   end type exact_row_t

   type :: exact_weights_t
      integer :: nx_out = 0, ny_out = 0, nx_in = 0, ny_in = 0
      integer :: lbnd_o(2) = 0
      type(exact_row_t), allocatable :: rows(:)
   end type exact_weights_t

   double precision, parameter :: ast__bad_eo = -1.79769313486231571d308
   ! A point within this distance (radians) of a clipping great circle
   ! counts as on it. Shared corners agree to round-off; measured
   ! round-off widths are at most 3.4e-17 rad.
   double precision, parameter :: on_plane_tol = 1.0d-15
   ! Corners closer than this (radians) are one point: a CAR pixel with
   ! an edge on a pole has two corners there and is a triangle.
   double precision, parameter :: dup_tol = 1.0d-15
   ! Overlaps thinner than this (area / output pixel perimeter, radians)
   ! are round-off from edges the two grids share, not overlap. A width,
   ! not a fraction of pixel area: a sliver of fixed width is a larger
   ! fraction of a smaller pixel. About 30x the largest round-off width
   ! measured (3.4e-17 rad, 0.36" to 900" pixels, 0-30 deg from the
   ! reference point) and about 2e-10 arcsec.
   double precision, parameter :: sliver_width = 1.0d-15
   integer, parameter :: max_vert = 16
   ! A corner lying exactly on a projection boundary -- native longitude
   ! +-180 deg on a CAR grid's seam -- gets no sky position from AST, and a
   ! pixel without all four corners would be dropped (measured: a 360-deg
   ! CAR grid lost both seam columns, 0.78% of a plane's flux). Such a
   ! corner is retried this far (pixels) inside each of the four pixels
   ! sharing it; the pixel's area changes by about this fraction.
   double precision, parameter :: boundary_nudge = 1.0d-9
   ! Quadrant qd seen from a corner: direction (quad_dx, quad_dy).
   double precision, parameter :: quad_dx(4) = (/ 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
   double precision, parameter :: quad_dy(4) = (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)
   ! Pixel (i, j)'s corners in order: grid corner (i + cdi, j + cdj), from
   ! which the pixel lies in quadrant cq.
   integer, parameter :: cdi(4) = (/ 0, 1, 1, 0 /), cdj(4) = (/ 0, 0, 1, 1 /)
   integer, parameter :: cq(4) = (/ 1, 2, 3, 4 /)

contains

   subroutine exact_init(wts, nx_in, ny_in, lbnd_o, ubnd_o)
      type(exact_weights_t), intent(inout) :: wts
      integer, intent(in) :: nx_in, ny_in, lbnd_o(2), ubnd_o(2)
      call exact_free(wts)
      wts%nx_in = nx_in; wts%ny_in = ny_in
      wts%nx_out = ubnd_o(1) - lbnd_o(1) + 1
      wts%ny_out = ubnd_o(2) - lbnd_o(2) + 1
      wts%lbnd_o = lbnd_o
      allocate(wts%rows(wts%ny_out))
   end subroutine exact_init

   subroutine exact_build_rows(wts, map_in2out, skymap_out, map_in2sky, lon_axis,&
   &jr0, jr1, status)
      !! Builds the weights of output rows jr0..jr1 (1-based within the
      !! output grid). Only wts%rows(jr0:jr1) is written.
      !! map_in2out: input GRID pixel -> output pixel; only its inverse is
      !!   used, to find which input pixels an output pixel can touch.
      !! skymap_out: output pixel -> output sky (radians), for output
      !!   corners.
      !! map_in2sky: input GRID pixel -> output sky (radians), for input
      !!   corners. Input corners must not go through output pixel
      !!   coordinates: near a pole (a whole row of a CAR grid) or beyond a
      !!   projection's edge (outside AIT's ellipse) a sky point has no
      !!   usable output pixel position, and the pixels using it would be
      !!   dropped (measured: B CAR and B AIT, REPROJECT_EXACT_PLAN.md).
      !! lon_axis: which of the two sky axes is longitude (1 or 2).
      type(exact_weights_t), intent(inout) :: wts
      integer, intent(in) :: map_in2out, skymap_out, map_in2sky, lon_axis, jr0, jr1
      integer, intent(inout) :: status

      integer :: nxo, ncx, ncy, i, j, k, ia, jb, li, lj, n_q, n_p, cnt, total, qd, m, nb
      integer :: in_i_lo, in_i_hi, in_j_lo, in_j_hi, nci, ncj, a0, a1, b0, b1
      double precision, allocatable :: oc_vec(:,:,:), oc_xin(:,:), oc_yin(:,:)
      double precision, allocatable :: ic_vec(:,:,:)
      logical, allocatable :: oc_ok(:,:), ic_ok(:,:)
      ! Corners AST rejects (no sky position, or no input-pixel position):
      ! four alternatives each, boundary_nudge pixel inside each of the four
      ! pixels sharing the corner, quadrant qd as in quad_dx/quad_dy.
      double precision, allocatable :: oa_vec(:,:,:,:), oa_xin(:,:,:), oa_yin(:,:,:)
      double precision, allocatable :: ia_vec(:,:,:,:)
      logical, allocatable :: oa_ok(:,:,:), ia_ok(:,:,:)
      double precision, allocatable :: xa(:), ya(:), xb(:), yb(:), xs(:), ys(:)
      integer, allocatable :: bi(:), bj(:)
      double precision :: xmin, xmax, ymin, ymax, q(3, 4), p(3, 4), a, a_sum, w_min
      double precision :: cx(4), cy(4)
      logical :: okc
      integer(int64), allocatable :: t_idx(:)
      real(dp), allocatable :: t_w(:)

      if (status.ne.0) return
      nxo = wts%nx_out
      ncx = nxo + 1
      ncy = jr1 - jr0 + 2

      ! -- output corners of these rows: sky unit vectors, and input
      !    pixel coordinates for the candidate search --
      allocate(oc_vec(3, ncx, ncy), oc_xin(ncx, ncy), oc_yin(ncx, ncy), oc_ok(ncx, ncy))
      allocate(xa(ncx), ya(ncx), xb(ncx), yb(ncx))
      do j = 1, ncy
         do i = 1, ncx
            xa(i) = dble(wts%lbnd_o(1) + i - 1) - 0.5d0
            ya(i) = dble(wts%lbnd_o(2) + jr0 - 1 + j - 1) - 0.5d0
         enddo
         call ast_tran2(skymap_out, ncx, xa, ya, .true., xb, yb, status)
         if (status.ne.0) return
         do i = 1, ncx
            oc_ok(i, j) = xb(i).ne.ast__bad_eo .and. yb(i).ne.ast__bad_eo
            if (oc_ok(i, j)) oc_vec(:, i, j) = sky_vec(xb(i), yb(i), lon_axis)
         enddo
         call ast_tran2(map_in2out, ncx, xa, ya, .false., xb, yb, status)
         if (status.ne.0) return
         do i = 1, ncx
            oc_xin(i, j) = xb(i); oc_yin(i, j) = yb(i)
            if (xb(i).eq.ast__bad_eo .or. yb(i).eq.ast__bad_eo) oc_ok(i, j) = .false.
         enddo
      enddo
      deallocate(xa, ya, xb, yb)

      ! -- alternatives for rejected output corners --
      allocate(oa_vec(3, 4, ncx, ncy), oa_xin(4, ncx, ncy), oa_yin(4, ncx, ncy), oa_ok(4, ncx, ncy))
      oa_ok = .false.
      nb = count(.not. oc_ok)
      if (nb.gt.0) then
         allocate(bi(nb), bj(nb), xa(4 * nb), ya(4 * nb), xb(4 * nb), yb(4 * nb), xs(4 * nb), ys(4 * nb))
         m = 0
         do j = 1, ncy
            do i = 1, ncx
               if (oc_ok(i, j)) cycle
               m = m + 1; bi(m) = i; bj(m) = j
               do qd = 1, 4
                  xa(4 * (m - 1) + qd) = dble(wts%lbnd_o(1) + i - 1) - 0.5d0 + quad_dx(qd) * boundary_nudge
                  ya(4 * (m - 1) + qd) = dble(wts%lbnd_o(2) + jr0 - 1 + j - 1) - 0.5d0 + quad_dy(qd) * boundary_nudge
               enddo
            enddo
         enddo
         call ast_tran2(skymap_out, 4 * nb, xa, ya, .true., xs, ys, status)
         if (status.ne.0) return
         call ast_tran2(map_in2out, 4 * nb, xa, ya, .false., xb, yb, status)
         if (status.ne.0) return
         do m = 1, nb
            do qd = 1, 4
               k = 4 * (m - 1) + qd
               oa_ok(qd, bi(m), bj(m)) = xs(k).ne.ast__bad_eo .and. ys(k).ne.ast__bad_eo .and.&
               &xb(k).ne.ast__bad_eo .and. yb(k).ne.ast__bad_eo
               if (oa_ok(qd, bi(m), bj(m))) then
                  oa_vec(:, qd, bi(m), bj(m)) = sky_vec(xs(k), ys(k), lon_axis)
                  oa_xin(qd, bi(m), bj(m)) = xb(k); oa_yin(qd, bi(m), bj(m)) = yb(k)
               endif
            enddo
         enddo
         deallocate(bi, bj, xa, ya, xb, yb, xs, ys)
      endif

      ! -- input pixels these rows can touch: the corners' input-pixel
      !    bounding box, one pixel wider (a great-circle edge can bow past
      !    the straight line between two mapped corners) --
      xmin = huge(1.0d0); xmax = -huge(1.0d0); ymin = huge(1.0d0); ymax = -huge(1.0d0)
      do j = 1, ncy
         do i = 1, ncx
            if (oc_ok(i, j)) then
               xmin = min(xmin, oc_xin(i, j)); xmax = max(xmax, oc_xin(i, j))
               ymin = min(ymin, oc_yin(i, j)); ymax = max(ymax, oc_yin(i, j))
            endif
            do qd = 1, 4
               if (oa_ok(qd, i, j)) then
                  xmin = min(xmin, oa_xin(qd, i, j)); xmax = max(xmax, oa_xin(qd, i, j))
                  ymin = min(ymin, oa_yin(qd, i, j)); ymax = max(ymax, oa_yin(qd, i, j))
               endif
            enddo
         enddo
      enddo
      if (xmin.gt.xmax) then
         in_i_lo = 1; in_i_hi = 0; in_j_lo = 1; in_j_hi = 0
      else
         in_i_lo = max(1, floor(xmin + 0.5d0) - 1)
         in_i_hi = min(wts%nx_in, ceiling(xmax - 0.5d0) + 1)
         in_j_lo = max(1, floor(ymin + 0.5d0) - 1)
         in_j_hi = min(wts%ny_in, ceiling(ymax - 0.5d0) + 1)
      endif
      nci = max(0, in_i_hi - in_i_lo + 2)
      ncj = max(0, in_j_hi - in_j_lo + 2)
      allocate(ic_vec(3, max(1, nci), max(1, ncj)), ic_ok(max(1, nci), max(1, ncj)))
      allocate(ia_vec(3, 4, max(1, nci), max(1, ncj)), ia_ok(4, max(1, nci), max(1, ncj)))
      ic_ok = .false.
      ia_ok = .false.
      if (nci.gt.1 .and. ncj.gt.1) then
         allocate(xa(nci), ya(nci), xs(nci), ys(nci))
         do j = 1, ncj
            do i = 1, nci
               xa(i) = dble(in_i_lo + i - 1) - 0.5d0
               ya(i) = dble(in_j_lo + j - 1) - 0.5d0
            enddo
            call ast_tran2(map_in2sky, nci, xa, ya, .true., xs, ys, status)
            if (status.ne.0) return
            do i = 1, nci
               ic_ok(i, j) = xs(i).ne.ast__bad_eo .and. ys(i).ne.ast__bad_eo
               if (ic_ok(i, j)) ic_vec(:, i, j) = sky_vec(xs(i), ys(i), lon_axis)
            enddo
         enddo
         deallocate(xa, ya, xs, ys)
         ! alternatives for rejected input corners
         nb = count(.not. ic_ok(1:nci, 1:ncj))
         if (nb.gt.0) then
            allocate(bi(nb), bj(nb), xa(4 * nb), ya(4 * nb), xs(4 * nb), ys(4 * nb))
            m = 0
            do j = 1, ncj
               do i = 1, nci
                  if (ic_ok(i, j)) cycle
                  m = m + 1; bi(m) = i; bj(m) = j
                  do qd = 1, 4
                     xa(4 * (m - 1) + qd) = dble(in_i_lo + i - 1) - 0.5d0 + quad_dx(qd) * boundary_nudge
                     ya(4 * (m - 1) + qd) = dble(in_j_lo + j - 1) - 0.5d0 + quad_dy(qd) * boundary_nudge
                  enddo
               enddo
            enddo
            call ast_tran2(map_in2sky, 4 * nb, xa, ya, .true., xs, ys, status)
            if (status.ne.0) return
            do m = 1, nb
               do qd = 1, 4
                  k = 4 * (m - 1) + qd
                  ia_ok(qd, bi(m), bj(m)) = xs(k).ne.ast__bad_eo .and. ys(k).ne.ast__bad_eo
                  if (ia_ok(qd, bi(m), bj(m))) ia_vec(:, qd, bi(m), bj(m)) = sky_vec(xs(k), ys(k), lon_axis)
               enddo
            enddo
            deallocate(bi, bj, xa, ya, xs, ys)
         endif
      endif

      ! -- overlaps, row by row --
      do j = 1, jr1 - jr0 + 1
         ! upper bound on this row's entries: every candidate box
         total = 0
         do i = 1, nxo
            call candidate_box(i, j, a0, a1, b0, b1)
            total = total + max(0, a1 - a0 + 1) * max(0, b1 - b0 + 1)
         enddo
         allocate(t_idx(max(1, total)), t_w(max(1, total)))
         allocate(wts%rows(jr0 + j - 1)%ptr(nxo + 1))
         cnt = 0
         do i = 1, nxo
            wts%rows(jr0 + j - 1)%ptr(i) = cnt + 1
            call candidate_box(i, j, a0, a1, b0, b1)
            if (a1.lt.a0 .or. b1.lt.b0) cycle
            call out_quad(i, j, q, okc)
            if (.not. okc) cycle
            call prepare_polygon(q, n_q)
            if (n_q.lt.3) cycle
            w_min = sliver_width * perimeter(q, n_q)
            a_sum = 0.0d0
            k = cnt
            do jb = b0, b1
               lj = jb - in_j_lo + 1
               do ia = a0, a1
                  li = ia - in_i_lo + 1
                  call in_quad(li, lj, p, okc)
                  if (.not. okc) cycle
                  call prepare_polygon(p, n_p)
                  if (n_p.lt.3) cycle
                  a = overlap_area(p, n_p, q, n_q)
                  if (a.le.w_min) cycle
                  cnt = cnt + 1
                  t_idx(cnt) = int(jb - 1, int64) * int(wts%nx_in, int64) + int(ia, int64)
                  t_w(cnt) = a
                  a_sum = a_sum + a
               enddo
            enddo
            if (cnt.gt.k) t_w(k + 1:cnt) = t_w(k + 1:cnt) / a_sum
         enddo
         wts%rows(jr0 + j - 1)%ptr(nxo + 1) = cnt + 1
         allocate(wts%rows(jr0 + j - 1)%idx(max(1, cnt)), wts%rows(jr0 + j - 1)%w(max(1, cnt)))
         if (cnt.gt.0) then
            wts%rows(jr0 + j - 1)%idx(1:cnt) = t_idx(1:cnt)
            wts%rows(jr0 + j - 1)%w(1:cnt) = t_w(1:cnt)
         endif
         deallocate(t_idx, t_w)
      enddo

   contains

      ! Corner c of pixel (i, j) is grid corner (i + cdi(c), j + cdj(c));
      ! seen from that corner the pixel lies in quadrant cq(c).
      subroutine out_corner(i, j, c, v, x, y, ok)
         integer, intent(in) :: i, j, c
         double precision, intent(out) :: v(3), x, y
         logical, intent(out) :: ok
         integer :: gi, gj
         gi = i + cdi(c); gj = j + cdj(c)
         if (oc_ok(gi, gj)) then
            v = oc_vec(:, gi, gj); x = oc_xin(gi, gj); y = oc_yin(gi, gj); ok = .true.
         else
            ok = oa_ok(cq(c), gi, gj)
            if (ok) then
               v = oa_vec(:, cq(c), gi, gj); x = oa_xin(cq(c), gi, gj); y = oa_yin(cq(c), gi, gj)
            else
               v = 0.0d0; x = 0.0d0; y = 0.0d0
            endif
         endif
      end subroutine out_corner

      subroutine out_quad(i, j, quad, ok)
         integer, intent(in) :: i, j
         double precision, intent(out) :: quad(3, 4)
         logical, intent(out) :: ok
         integer :: c
         logical :: okk
         double precision :: x, y
         ok = .true.
         do c = 1, 4
            call out_corner(i, j, c, quad(:, c), x, y, okk)
            ok = ok .and. okk
         enddo
      end subroutine out_quad

      subroutine in_quad(li, lj, quad, ok)
         integer, intent(in) :: li, lj
         double precision, intent(out) :: quad(3, 4)
         logical, intent(out) :: ok
         integer :: c, gi, gj
         ok = .true.
         do c = 1, 4
            gi = li + cdi(c); gj = lj + cdj(c)
            if (ic_ok(gi, gj)) then
               quad(:, c) = ic_vec(:, gi, gj)
            else if (ia_ok(cq(c), gi, gj)) then
               quad(:, c) = ia_vec(:, cq(c), gi, gj)
            else
               ok = .false.
               return
            endif
         enddo
      end subroutine in_quad

      subroutine candidate_box(i, j, a0, a1, b0, b1)
         integer, intent(in) :: i, j
         integer, intent(out) :: a0, a1, b0, b1
         double precision :: v(3)
         integer :: c
         logical :: okk
         a0 = 1; a1 = 0; b0 = 1; b1 = 0
         do c = 1, 4
            call out_corner(i, j, c, v, cx(c), cy(c), okk)
            if (.not. okk) return
         enddo
         a0 = max(in_i_lo, floor(minval(cx) + 0.5d0) - 1)
         a1 = min(in_i_hi, ceiling(maxval(cx) - 0.5d0) + 1)
         b0 = max(in_j_lo, floor(minval(cy) + 0.5d0) - 1)
         b1 = min(in_j_hi, ceiling(maxval(cy) - 0.5d0) + 1)
      end subroutine candidate_box
   end subroutine exact_build_rows

   subroutine exact_apply(wts, plane_in, plane_out)
      !! plane_out = overlap-weighted mean of plane_in (whole planes, x
      !! fastest, nx_in*ny_in and nx_out*ny_out values); NaN where an
      !! overlapped input pixel is NaN or where nothing overlaps.
      type(exact_weights_t), intent(in) :: wts
      real, intent(in) :: plane_in(*)
      real, intent(inout) :: plane_out(*)
      call exact_apply_rows(wts, plane_in, plane_out, 1, wts%ny_out)
   end subroutine exact_apply

   subroutine exact_apply_rows(wts, plane_in, plane_out, r0, r1)
      !! As exact_apply, for output rows r0..r1 only (1-based within the
      !! output grid); the rest of plane_out is not touched, so threads
      !! that split a plane's rows can share one output plane.
      type(exact_weights_t), intent(in) :: wts
      real, intent(in) :: plane_in(*)
      real, intent(inout) :: plane_out(*)
      integer, intent(in) :: r0, r1
      integer :: i, j, e
      integer(int64) :: k
      real(dp) :: acc
      real :: nan
      nan = ieee_value(nan, ieee_quiet_nan)
      do j = r0, r1
         do i = 1, wts%nx_out
            k = int(j - 1, int64) * int(wts%nx_out, int64) + int(i, int64)
            if (wts%rows(j)%ptr(i + 1).le.wts%rows(j)%ptr(i)) then
               plane_out(k) = nan
               cycle
            endif
            acc = 0.0_dp
            do e = wts%rows(j)%ptr(i), wts%rows(j)%ptr(i + 1) - 1
               acc = acc + wts%rows(j)%w(e) * real(plane_in(wts%rows(j)%idx(e)), dp)
            enddo
            plane_out(k) = real(acc)
         enddo
      enddo
   end subroutine exact_apply_rows

   integer(int64) function exact_weights_bytes(wts) result(nbytes)
      type(exact_weights_t), intent(in) :: wts
      integer :: j
      nbytes = 0
      if (.not. allocated(wts%rows)) return
      do j = 1, size(wts%rows)
         if (allocated(wts%rows(j)%ptr)) nbytes = nbytes + 4_int64 * size(wts%rows(j)%ptr, kind=int64)
         if (allocated(wts%rows(j)%idx)) nbytes = nbytes + 16_int64 * size(wts%rows(j)%idx, kind=int64)
      enddo
   end function exact_weights_bytes

   integer(int64) function exact_estimate_bytes(nx_out, ny_out, out_over_in) result(nbytes)
      !! Estimated weight-table size before it is built: about
      !! (out_over_in + 1)^2 entries per output pixel (input pixels a
      !! general-position output pixel of that relative size overlaps),
      !! 16 bytes each, plus the row pointers.
      integer, intent(in) :: nx_out, ny_out
      double precision, intent(in) :: out_over_in
      nbytes = int(dble(nx_out) * dble(ny_out) * (out_over_in + 1.0d0)**2 * 16.0d0, int64)&
      & + 4_int64 * int(nx_out + 1, int64) * int(ny_out, int64)
   end function exact_estimate_bytes

   function mib_str(nbytes) result(str)
      !! nbytes as MiB with one decimal and a leading zero below 1 (the F0
      !! edit descriptor drops it).
      integer(int64), intent(in) :: nbytes
      character(len=24) :: str
      write(str, '(F24.1)') real(nbytes, dp) / 1048576.0_dp
      str = adjustl(str)
   end function mib_str

   subroutine exact_free(wts)
      type(exact_weights_t), intent(inout) :: wts
      if (allocated(wts%rows)) deallocate(wts%rows)
      wts%nx_out = 0; wts%ny_out = 0
   end subroutine exact_free

   pure function sky_vec(a1, a2, lon_axis) result(v)
      !! Unit vector of a sky position given as AST's two sky-axis values
      !! (radians), in whichever order the SkyFrame has them.
      double precision, intent(in) :: a1, a2
      integer, intent(in) :: lon_axis
      double precision :: v(3), lon, lat
      if (lon_axis.eq.2) then
         lon = a2; lat = a1
      else
         lon = a1; lat = a2
      endif
      v = (/ cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat) /)
   end function sky_vec

   pure function cross(a, b) result(c)
      double precision, intent(in) :: a(3), b(3)
      double precision :: c(3)
      c = (/ a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1) /)
   end function cross

   subroutine prepare_polygon(c, n)
      !! Merges consecutive corners closer than dup_tol (and a last
      !! corner equal to the first), then orders the remaining n corners
      !! counter-clockwise seen from outside the sphere. Orientation comes
      !! from the total signed area of the fan, not from the first three
      !! corners, which may coincide. n < 3 means degenerate.
      double precision, intent(inout) :: c(3, 4)
      integer, intent(out) :: n
      double precision :: tmp(3, 4), t
      integer :: k
      n = 0
      do k = 1, 4
         if (n.gt.0) then
            if (norm(c(:, k) - tmp(:, n)).le.dup_tol) cycle
         endif
         n = n + 1
         tmp(:, n) = c(:, k)
      enddo
      if (n.gt.1) then
         if (norm(tmp(:, 1) - tmp(:, n)).le.dup_tol) n = n - 1
      endif
      if (n.lt.3) then
         n = 0
         return
      endif
      t = 0.0d0
      do k = 2, n - 1
         t = t + dot_product(tmp(:, 1), cross(tmp(:, k) - tmp(:, 1), tmp(:, k + 1) - tmp(:, 1)))
      enddo
      if (t.eq.0.0d0) then
         n = 0
         return
      endif
      if (t.gt.0.0d0) then
         c(:, 1:n) = tmp(:, 1:n)
      else
         do k = 1, n
            c(:, k) = tmp(:, n - k + 1)
         enddo
      endif
   end subroutine prepare_polygon

   pure double precision function norm(v)
      double precision, intent(in) :: v(3)
      norm = sqrt(dot_product(v, v))
   end function norm

   double precision function perimeter(c, n) result(s)
      double precision, intent(in) :: c(3, 4)
      integer, intent(in) :: n
      integer :: k
      s = 0.0d0
      do k = 1, n
         s = s + norm(c(:, mod(k, n) + 1) - c(:, k))
      enddo
   end function perimeter

   double precision function overlap_area(subj, ns, clp, nc) result(a)
      !! Area (sr) of polygon subj clipped to polygon clp, both convex and
      !! counter-clockwise: Sutherland-Hodgman against each great circle
      !! of clp, then the spherical excess. Every product is taken on
      !! differences of nearby vertices, which floating point subtracts
      !! exactly, never on two nearly equal unit vectors: for arcsecond
      !! pixels the direct cross product of two corners loses about
      !! 1e-16 / (edge length) of its direction (measured: false overlaps
      !! 2.8e-11 rad wide without this, 3.4e-17 rad with it).
      double precision, intent(in) :: subj(3, 4), clp(3, 4)
      integer, intent(in) :: ns, nc
      double precision :: pin(3, max_vert), pout(3, max_vert), nrm(3), ps(3), pe(3), x(3), c0(3)
      double precision :: ds, de
      integer :: nin, nout, e, k
      nout = ns
      pout(:, 1:ns) = subj(:, 1:ns)
      do e = 1, nc
         if (nout.eq.0) exit
         c0 = clp(:, e)
         nrm = cross(c0, clp(:, mod(e, nc) + 1) - c0)
         nrm = nrm / norm(nrm)
         nin = nout
         pin(:, 1:nin) = pout(:, 1:nin)
         nout = 0
         do k = 1, nin
            ps = pin(:, merge(nin, k - 1, k.eq.1))
            pe = pin(:, k)
            ds = dot_product(nrm, ps - c0)
            de = dot_product(nrm, pe - c0)
            if (de.ge.-on_plane_tol) then
               if (ds.lt.-on_plane_tol) then
                  ! sign-correct crossing: dividing by (ds - de) keeps the
                  ! point on the arc between ps and pe, not its antipode
                  x = ps + (ds / (ds - de)) * (pe - ps)
                  nout = nout + 1; pout(:, nout) = x / norm(x)
               endif
               nout = nout + 1; pout(:, nout) = pe
            else if (ds.ge.-on_plane_tol) then
               x = ps + (ds / (ds - de)) * (pe - ps)
               nout = nout + 1; pout(:, nout) = x / norm(x)
            endif
            if (nout.ge.max_vert - 1) exit
         enddo
      enddo
      a = polygon_area(pout, nout)
   end function overlap_area

   double precision function polygon_area(p, n) result(a)
      !! Spherical excess of convex polygon p(:,1:n) as a fan from p(:,1):
      !! tan(E/2) = |a.(b x c)| / (1 + a.b + b.c + c.a) (Van Oosterom and
      !! Strackee), the triple product taken on differences from p(:,1).
      double precision, intent(in) :: p(:, :)
      integer, intent(in) :: n
      double precision :: t, d
      integer :: k
      a = 0.0d0
      if (n.lt.3) return
      do k = 2, n - 1
         t = dot_product(p(:, 1), cross(p(:, k) - p(:, 1), p(:, k + 1) - p(:, 1)))
         d = 1.0d0 + dot_product(p(:, 1), p(:, k)) + dot_product(p(:, k), p(:, k + 1)) +&
         &dot_product(p(:, k + 1), p(:, 1))
         a = a + 2.0d0 * atan2(t, d)
      enddo
      a = abs(a)
   end function polygon_area

end module exact_overlap_mod
