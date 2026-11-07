module rm_cuda_backend
  !! Optional CUDA backend for the RM-synthesis DFT (build with GPU_BACKEND=cuda). The sums
  !!   Re = sum_k w_k (Q_k cos_k - U_k sin_k),  Im = sum_k w_k (Q_k sin_k + U_k cos_k)
  !! (rm_synthesis's cos/sin templates, with the sign convention of tile_extract_gpu_rm_blocked) are one real matrix product
  !!   [Re | Im] (npix x 2 nrm) = [wQ | wU] (npix x 2 nz) * [[cos, sin], [-sin, cos]] (2 nz x 2 nrm)
  !! done by cuBLAS Sgemm on chunks of pixels. The data layout is the one prepare_gpu_data already produces (pixel fastest),
  !! and the result lands in the (pixel, RM) layout of p_tile_arr / phi_tile_arr. The per-pixel weights are folded into the
  !! left matrix, so flagged samples and flagged channels need no special case; the epilogue divides by the per-pixel weight sum.
  use iso_c_binding
  use omp_lib
  use cublas_iface
  implicit none
  private
  public :: rm_cuda_extract

  type(c_ptr), save :: blas_handle = c_null_ptr
  logical, save :: blas_ready = .false.
  ! Pixels per matrix product: bounds the temporary device buffers (about (2 nz + 2 nrm) * 4 bytes per pixel).
  integer, parameter :: pix_chunk = 262144

contains

  subroutine check_blas(status, what)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= CUBLAS_STATUS_SUCCESS) then
      write(*,'(A,A,A,I0)') ' ERROR: cuBLAS call ', what, ' failed with status ', status
      error stop 1
    end if
  end subroutine check_blas

  function dev_alloc(nbytes, dev) result(p)
    integer(c_size_t), intent(in) :: nbytes
    integer, intent(in) :: dev
    type(c_ptr) :: p
    p = omp_target_alloc(max(nbytes, 4_c_size_t), dev)
    if (.not. c_associated(p)) then
      write(*,'(A,F10.1,A)') ' ERROR: could not allocate device memory (MiB): ', real(nbytes) / 1048576.0, &
        '. Lower mem_frac_vram / gpu_vram_mib.'
      error stop 1
    end if
  end function dev_alloc

  subroutine h2d(dst, src, nbytes, dev)
    type(c_ptr), intent(in) :: dst, src
    integer(c_size_t), intent(in) :: nbytes
    integer, intent(in) :: dev
    integer(c_int) :: rc
    rc = omp_target_memcpy(dst, src, nbytes, 0_c_size_t, 0_c_size_t, dev, omp_get_initial_device())
    if (rc /= 0) error stop 'omp_target_memcpy host-to-device failed'
  end subroutine h2d

  subroutine rm_cuda_extract(specQ, specU, wts, mean_q, mean_u, wsum, cos_arr, sin_arr, &
                             nx_tile, ny_tile, nz_out, i_rm_block, nrm_block_now, &
                             rem_mean, output_mode, ap_angle_mode, p_tile_arr, phi_tile_arr)
    integer, intent(in) :: nx_tile, ny_tile, nz_out, i_rm_block, nrm_block_now, rem_mean, output_mode, ap_angle_mode
    real(c_float), intent(in), target, contiguous :: specQ(:,:), specU(:,:), wts(:,:)
    real(c_float), intent(in), target, optional :: mean_q(nx_tile*ny_tile), mean_u(nx_tile*ny_tile)
    real(c_float), intent(in), target :: wsum(nx_tile*ny_tile)
    real(c_float), intent(in) :: cos_arr(:,:), sin_arr(:,:)
    real(c_float), intent(inout), target, contiguous :: p_tile_arr(:), phi_tile_arr(:)

    integer(c_int64_t) :: npix, p0, pg, pidx
    integer :: nz, nrm, nc, cap, j, col, dev, ip, iz, ir
    integer(c_int) :: st, rc
    integer(c_size_t) :: nbytes_pix_nz, nbytes_out, offs
    real(c_float), allocatable, target :: bmat(:,:)
    type(c_ptr) :: d_q, d_u, d_w, d_ws, d_mq, d_mu, d_b, d_a, d_c, d_amp, d_pha
    real(c_float), pointer :: q(:,:), u(:,:), w(:,:), a(:,:), c(:,:), ws(:), mq(:), mu(:), amp(:), pha(:)
    real(c_float) :: one, zero, mqv, muv, wt, wsv, re, im, ang, nanv, zero_v

    one = 1.0; zero = 0.0; zero_v = 0.0
    npix = int(nx_tile, c_int64_t) * int(ny_tile, c_int64_t); nz = nz_out; nrm = nrm_block_now
    dev = omp_get_default_device()
    if (.not. blas_ready) then
      call check_blas(cublasCreate_v2(blas_handle), 'cublasCreate'); blas_ready = .true.
    end if

    ! Template matrix (2 nz x 2 nrm) for this RM block, from rm_synthesis's own cos/sin tables.
    allocate(bmat(2*nz, 2*nrm))
    do j = 1, nrm
      col = i_rm_block + j - 1
      bmat(1:nz, j)             =  cos_arr(1:nz, col)
      bmat(nz+1:2*nz, j)        = -sin_arr(1:nz, col)
      bmat(1:nz, nrm+j)         =  sin_arr(1:nz, col)
      bmat(nz+1:2*nz, nrm+j)    =  cos_arr(1:nz, col)
    end do

    nbytes_pix_nz = int(npix, c_size_t) * int(nz, c_size_t) * 4_c_size_t
    nbytes_out = int(npix, c_size_t) * int(nrm, c_size_t) * 4_c_size_t
    cap = int(min(int(pix_chunk, c_int64_t), npix))
    d_q = dev_alloc(nbytes_pix_nz, dev); d_u = dev_alloc(nbytes_pix_nz, dev); d_w = dev_alloc(nbytes_pix_nz, dev)
    d_ws = dev_alloc(int(npix, c_size_t) * 4_c_size_t, dev)
    d_b = dev_alloc(int(4*nz, c_size_t) * int(nrm, c_size_t) * 4_c_size_t, dev)
    d_a = dev_alloc(int(cap, c_size_t) * int(2*nz, c_size_t) * 4_c_size_t, dev)
    d_c = dev_alloc(int(cap, c_size_t) * int(2*nrm, c_size_t) * 4_c_size_t, dev)
    d_amp = dev_alloc(nbytes_out, dev); d_pha = dev_alloc(nbytes_out, dev)
    if (rem_mean > 0 .and. present(mean_q)) then
      d_mq = dev_alloc(int(npix, c_size_t) * 4_c_size_t, dev); d_mu = dev_alloc(int(npix, c_size_t) * 4_c_size_t, dev)
      call h2d(d_mq, c_loc(mean_q), int(npix, c_size_t) * 4_c_size_t, dev)
      call h2d(d_mu, c_loc(mean_u), int(npix, c_size_t) * 4_c_size_t, dev)
    else
      d_mq = dev_alloc(4_c_size_t, dev); d_mu = dev_alloc(4_c_size_t, dev)
    end if
    call h2d(d_q, c_loc(specQ), nbytes_pix_nz, dev); call h2d(d_u, c_loc(specU), nbytes_pix_nz, dev)
    call h2d(d_w, c_loc(wts), nbytes_pix_nz, dev); call h2d(d_ws, c_loc(wsum), int(npix, c_size_t) * 4_c_size_t, dev)
    call h2d(d_b, c_loc(bmat), int(4*nz, c_size_t) * int(nrm, c_size_t) * 4_c_size_t, dev)

    call c_f_pointer(d_q, q, [int(npix), nz]); call c_f_pointer(d_u, u, [int(npix), nz]); call c_f_pointer(d_w, w, [int(npix), nz])
    call c_f_pointer(d_ws, ws, [int(npix)]); call c_f_pointer(d_amp, amp, [int(npix) * nrm]); call c_f_pointer(d_pha, pha, [int(npix) * nrm])
    call c_f_pointer(d_mq, mq, [1]); call c_f_pointer(d_mu, mu, [1])
    if (rem_mean > 0 .and. present(mean_q)) then
      call c_f_pointer(d_mq, mq, [int(npix)]); call c_f_pointer(d_mu, mu, [int(npix)])
    end if

    do p0 = 1_c_int64_t, npix, int(pix_chunk, c_int64_t)
      nc = int(min(int(pix_chunk, c_int64_t), npix - p0 + 1_c_int64_t))
      call c_f_pointer(d_a, a, [nc, 2*nz]); call c_f_pointer(d_c, c, [nc, 2*nrm])

      ! left matrix [w (Q - meanQ) | w (U - meanU)] for this pixel chunk; pixel is the SIMD loop, as in the OpenMP kernel
      !$omp target teams distribute parallel do simd is_device_ptr(q, u, w, a, mq, mu) private(pg, iz, mqv, muv, wt)
      do ip = 1, nc
        pg = p0 + int(ip - 1, c_int64_t)
        mqv = 0.0; muv = 0.0
        if (rem_mean > 0) then
          mqv = mq(pg); muv = mu(pg)
        end if
        do iz = 1, nz
          wt = w(pg, iz)
          a(ip, iz) = wt * (q(pg, iz) - mqv)
          a(ip, nz + iz) = wt * (u(pg, iz) - muv)
        end do
      end do

      call check_blas(cublasSgemm_v2(blas_handle, CUBLAS_OP_N, CUBLAS_OP_N, nc, 2*nrm, 2*nz, one, d_a, nc, d_b, 2*nz, zero, d_c, nc), 'Sgemm')

      ! epilogue: divide by the per-pixel weight sum, then amplitude/angle (or real/imaginary), NaN where no channel is valid
      !$omp target teams distribute parallel do simd is_device_ptr(c, ws, amp, pha) private(pg, wsv, ir, pidx, re, im, ang, nanv)
      do ip = 1, nc
        pg = p0 + int(ip - 1, c_int64_t)
        wsv = ws(pg)
        do ir = 1, nrm
          pidx = pg + int(ir - 1, c_int64_t) * npix
          if (wsv > 0.0) then
            re = c(ip, ir) / wsv; im = c(ip, nrm + ir) / wsv
            if (output_mode == 1) then
              amp(pidx) = re; pha(pidx) = im
            else
              amp(pidx) = sqrt(re * re + im * im)
              ang = atan2(im, re)
              if (ap_angle_mode == 1) ang = 0.5 * ang
              pha(pidx) = ang
            end if
          else
            nanv = zero_v / zero_v
            amp(pidx) = nanv; pha(pidx) = nanv
          end if
        end do
      end do
    end do

    ! only the RM planes of this block come back (they are contiguous in the (pixel, RM) layout)
    offs = int(i_rm_block - 1, c_size_t) * int(npix, c_size_t) * 4_c_size_t
    rc = omp_target_memcpy(c_loc(p_tile_arr), d_amp, nbytes_out, offs, 0_c_size_t, omp_get_initial_device(), dev)
    if (rc /= 0) error stop 'omp_target_memcpy device-to-host (amplitude) failed'
    rc = omp_target_memcpy(c_loc(phi_tile_arr), d_pha, nbytes_out, offs, 0_c_size_t, omp_get_initial_device(), dev)
    if (rc /= 0) error stop 'omp_target_memcpy device-to-host (phase) failed'

    call omp_target_free(d_q, dev); call omp_target_free(d_u, dev); call omp_target_free(d_w, dev); call omp_target_free(d_ws, dev)
    call omp_target_free(d_b, dev); call omp_target_free(d_a, dev); call omp_target_free(d_c, dev)
    call omp_target_free(d_amp, dev); call omp_target_free(d_pha, dev); call omp_target_free(d_mq, dev); call omp_target_free(d_mu, dev)
  end subroutine rm_cuda_extract

end module rm_cuda_backend
