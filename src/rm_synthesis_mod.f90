module rm_synthesis_mod
  !! Modern Fortran module for RM-synthesis extraction routines
  !! Wraps legacy fixed-form subroutines with explicit interfaces
  !! Author: Wasim Raja (modernized 2026)
  
  use iso_fortran_env, only: sp => real32, dp => real64, int8, int16, int32, int64
  implicit none
  
  private
  public :: extract_general_setup, extract_general, extract_general_ri
  public :: extract_general_w, extract_general_ri_w
  public :: tile_extract_gpu
  public :: linspace, nchar
  public :: read_cfg_keyval
  public :: write_runtime_estimate
  public :: sp, dp, int32, int64
  
  ! Include file parameters for RM-synthesis
  integer, parameter :: max_axis = 100
  integer, parameter :: max_ra = 1024
  integer, parameter :: max_dec = 1024
  integer, parameter :: maxchan = 256
  integer, parameter :: max_pix = 134217728  ! 512 MB in real32
  integer, parameter :: maxofac = 16
  integer, parameter :: maxnt = maxchan * maxofac
  
  ! Physical constants
  ! Speed of light in units of 10^6 m/s (for freq[MHz] <-> lambda[m] conversion)
  real(sp), parameter :: c_velocity = 299.792458_sp
  
  public :: max_axis, max_ra, max_dec, maxchan, max_pix, maxofac, maxnt
  public :: c_velocity

contains

  subroutine extract_general_setup(t, npts, fac, beg_rm, end_rm, nout, nu, cos_arr, sin_arr, maxout, maxpts, use_auto_rm_range, ofac)
    !! Pre-compute sine and cosine templates for RM-extraction
    !! This avoids redundant trig calculations across multiple pixels
    !! use_auto_rm_range: 1=derive beg/end/nrm from data, 0=use user beg/end
    !! nout is final output depth and should be nrm * ofac
    !! See: extract_general_setup.f for original implementation
    
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, use_auto_rm_range
    integer(int32), intent(in) :: ofac
    real(sp), intent(in) :: t(*), fac, beg_rm, end_rm
    real(sp), intent(out) :: nu(*)
    real(sp), intent(out) :: cos_arr(maxpts, maxout), sin_arr(maxpts, maxout)
    
    real(sp) :: freq_MHz(npts), f1, f2, Lsq1, Lsq2, dfreq
    real(sp) :: t_span, d_nu, nu_span, omega, h_tmp, phi_tmp, beg_eff, end_eff
    real(sp) :: neg_span, pos_span, span_ratio
    integer(int32) :: i, j, kk, nneg, npos, zero_idx
    real(sp), parameter :: pi = 3.14159265358979, twopi = 6.28318530717959
    
    ! Generate temporal frequencies from L_sq data
    ! t is lambda_squared (wavelength in meters squared)
    ! freq_MHz is frequency in MHz
    ! Using c = 299.792458 × 10^6 m/s (speed of light)
    j = npts + 1
    do kk = 1, npts
      j = j - 1
      freq_MHz(j) = c_velocity / sqrt(t(kk))
    end do
    
    ! Calculate edge L_sq
    dfreq = (freq_MHz(npts) - freq_MHz(1)) / dble(npts - 1)
    f1 = freq_MHz(1) - 0.5_sp * dfreq
    f2 = freq_MHz(npts) + 0.5_sp * dfreq
    Lsq2 = (c_velocity / f1)**2
    Lsq1 = (c_velocity / f2)**2
    
    ! Relation between RM and wavelength-squared domains
    t_span = Lsq2 - Lsq1
    d_nu = fac / t_span
    nu_span = dble(npts) * d_nu
    
    ! Build RM limits for the final nout samples.
    if (use_auto_rm_range == 1) then
      beg_eff = -0.5_sp * real(npts - 1) * d_nu
      end_eff =  0.5_sp * real(npts - 1) * d_nu
    else
      beg_eff = beg_rm
      end_eff = end_rm
    end if

    if (nout <= 1) then
      nu(1) = beg_eff
    else if (nout >= 3 .and. beg_eff < 0.0_sp .and. end_eff > 0.0_sp) then
      neg_span = abs(beg_eff)
      pos_span = abs(end_eff)
      span_ratio = neg_span / (neg_span + pos_span)

      nneg = nint(real(nout - 1, kind=sp) * span_ratio)
      nneg = max(1_int32, min(nout - 2, nneg))
      npos = (nout - 1) - nneg

      h_tmp = min(neg_span / real(nneg, kind=sp), &
                  pos_span / real(npos, kind=sp))
      zero_idx = nneg + 1

      do i = 1, nout
        nu(i) = real(i - zero_idx, kind=sp) * h_tmp
      end do
    else
      h_tmp = (end_eff - beg_eff) / real(nout - 1)
      do i = 1, nout
        nu(i) = beg_eff + real(i - 1) * h_tmp
      end do
    end if
    
    ! Pre-compute cos and sin templates
    do i = 1, nout
      omega = 2.0_sp * nu(i)
      do kk = 1, npts
        phi_tmp = omega * t(kk)
        cos_arr(kk, i) = cos(phi_tmp)
        sin_arr(kk, i) = -sin(phi_tmp)
      end do
    end do
    
  end subroutine extract_general_setup

  subroutine extract_general(ryt_in, iyt_in, npts, nout, p_ex, phi_ex, &
                             cos_arr, sin_arr, maxout, maxpts, mean_rem)
    !! Extract RM power using pre-computed templates
    !! Uses only dot products (no trig recomputation)
    !! See: extract_general_v4.f for original implementation
    
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, mean_rem
    real(sp), intent(in) :: ryt_in(*), iyt_in(*)
    real(sp), intent(out) :: p_ex(*), phi_ex(*)
    real(sp), intent(in) :: cos_arr(maxpts, maxout), sin_arr(maxpts, maxout)
    
    real(sp) :: ryt(npts), iyt(npts)
    real(sp) :: rc_cor, ic_cor, rs_cor, is_cor, ryw_tmp, iyw_tmp
    integer(int32) :: i, kk
    
    ! Remove mean from Q and U if requested
    if (mean_rem > 0) then
      call compute_mean(ryt_in, npts, ryw_tmp)
      call compute_mean(iyt_in, npts, iyw_tmp)
      do i = 1, npts
        ryt(i) = ryt_in(i) - ryw_tmp
        iyt(i) = iyt_in(i) - iyw_tmp
      end do
    else
      do i = 1, npts
        ryt(i) = ryt_in(i)
        iyt(i) = iyt_in(i)
      end do
    end if
    
    ! Extract using pre-computed templates.
    ! One fused loop reduces memory traffic versus copying template vectors
    ! and calling 4 separate dot products per RM bin.
    !$omp parallel do default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
    !$omp shared(nout,npts,ryt,iyt,cos_arr,sin_arr,p_ex,phi_ex)
    do i = 1, nout
      rc_cor = 0.0_sp
      rs_cor = 0.0_sp
      ic_cor = 0.0_sp
      is_cor = 0.0_sp
      !$omp simd reduction(+:rc_cor,rs_cor,ic_cor,is_cor)
      do kk = 1, npts
        rc_cor = rc_cor + ryt(kk) * cos_arr(kk, i)
        rs_cor = rs_cor + ryt(kk) * sin_arr(kk, i)
        ic_cor = ic_cor + iyt(kk) * cos_arr(kk, i)
        is_cor = is_cor + iyt(kk) * sin_arr(kk, i)
      end do
      
      rc_cor = rc_cor / dble(npts)
      rs_cor = rs_cor / dble(npts)
      ic_cor = ic_cor / dble(npts)
      is_cor = is_cor / dble(npts)
      
      ! Combine coherently to construct y(omega)
      ryw_tmp = rc_cor - is_cor
      iyw_tmp = rs_cor + ic_cor
      p_ex(i) = sqrt(ryw_tmp**2 + iyw_tmp**2)
      phi_ex(i) = atan2(iyw_tmp, ryw_tmp)
    end do
    !$omp end parallel do
    
  end subroutine extract_general

  subroutine extract_general_ri(ryt_in, iyt_in, npts, nout, re_ex, im_ex, &
                                cos_arr, sin_arr, maxout, maxpts, mean_rem)
    !! Extract RM complex spectrum directly as REAL/IMAG outputs
    !! Avoids amplitude/phase conversion when RI mode is requested
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, mean_rem
    real(sp), intent(in) :: ryt_in(*), iyt_in(*)
    real(sp), intent(out) :: re_ex(*), im_ex(*)
    real(sp), intent(in) :: cos_arr(maxpts, maxout), sin_arr(maxpts, maxout)

    real(sp) :: ryt(npts), iyt(npts)
    real(sp) :: rc_cor, ic_cor, rs_cor, is_cor, ryw_tmp, iyw_tmp
    integer(int32) :: i, kk

    if (mean_rem > 0) then
      call compute_mean(ryt_in, npts, ryw_tmp)
      call compute_mean(iyt_in, npts, iyw_tmp)
      do i = 1, npts
        ryt(i) = ryt_in(i) - ryw_tmp
        iyt(i) = iyt_in(i) - iyw_tmp
      end do
    else
      do i = 1, npts
        ryt(i) = ryt_in(i)
        iyt(i) = iyt_in(i)
      end do
    end if

    !$omp parallel do default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
    !$omp shared(nout,npts,ryt,iyt,cos_arr,sin_arr,re_ex,im_ex)
    do i = 1, nout
      rc_cor = 0.0_sp
      rs_cor = 0.0_sp
      ic_cor = 0.0_sp
      is_cor = 0.0_sp
      !$omp simd reduction(+:rc_cor,rs_cor,ic_cor,is_cor)
      do kk = 1, npts
        rc_cor = rc_cor + ryt(kk) * cos_arr(kk, i)
        rs_cor = rs_cor + ryt(kk) * sin_arr(kk, i)
        ic_cor = ic_cor + iyt(kk) * cos_arr(kk, i)
        is_cor = is_cor + iyt(kk) * sin_arr(kk, i)
      end do

      rc_cor = rc_cor / dble(npts)
      rs_cor = rs_cor / dble(npts)
      ic_cor = ic_cor / dble(npts)
      is_cor = is_cor / dble(npts)

      ryw_tmp = rc_cor - is_cor
      iyw_tmp = rs_cor + ic_cor
      re_ex(i) = ryw_tmp
      im_ex(i) = iyw_tmp
    end do
    !$omp end parallel do

  end subroutine extract_general_ri

  subroutine extract_general_w(ryt_in, iyt_in, wts_in, npts, nout, p_ex, phi_ex, &
                               cos_arr, sin_arr, maxout, maxpts, mean_rem)
    !! Weighted RM extraction (AP mode).
    !! Channels with zero weight are ignored via weighted sums.
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, mean_rem
    real(sp), intent(in) :: ryt_in(*), iyt_in(*), wts_in(*)
    real(sp), intent(out) :: p_ex(*), phi_ex(*)
    real(sp), intent(in) :: cos_arr(maxpts, maxout), sin_arr(maxpts, maxout)

    real(sp) :: ryt(npts), iyt(npts), wts(npts)
    real(sp) :: rc_cor, ic_cor, rs_cor, is_cor, ryw_tmp, iyw_tmp
    real(sp) :: wsum, mean_q, mean_u
    integer(int32) :: i, kk

    do kk = 1, npts
      wts(kk) = max(0.0_sp, wts_in(kk))
    end do

    wsum = 0.0_sp
    !$omp simd reduction(+:wsum)
    do kk = 1, npts
      wsum = wsum + wts(kk)
    end do

    if (mean_rem > 0 .and. wsum > 0.0_sp) then
      mean_q = 0.0_sp
      mean_u = 0.0_sp
      !$omp simd reduction(+:mean_q,mean_u)
      do kk = 1, npts
        mean_q = mean_q + wts(kk) * ryt_in(kk)
        mean_u = mean_u + wts(kk) * iyt_in(kk)
      end do
      mean_q = mean_q / wsum
      mean_u = mean_u / wsum
      do i = 1, npts
        ryt(i) = ryt_in(i) - mean_q
        iyt(i) = iyt_in(i) - mean_u
      end do
    else
      do i = 1, npts
        ryt(i) = ryt_in(i)
        iyt(i) = iyt_in(i)
      end do
    end if

    !$omp parallel do default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
    !$omp shared(nout,npts,ryt,iyt,wts,wsum,cos_arr,sin_arr,p_ex,phi_ex)
    do i = 1, nout
      rc_cor = 0.0_sp
      rs_cor = 0.0_sp
      ic_cor = 0.0_sp
      is_cor = 0.0_sp
      !$omp simd reduction(+:rc_cor,rs_cor,ic_cor,is_cor)
      do kk = 1, npts
        rc_cor = rc_cor + wts(kk) * ryt(kk) * cos_arr(kk, i)
        rs_cor = rs_cor + wts(kk) * ryt(kk) * sin_arr(kk, i)
        ic_cor = ic_cor + wts(kk) * iyt(kk) * cos_arr(kk, i)
        is_cor = is_cor + wts(kk) * iyt(kk) * sin_arr(kk, i)
      end do

      if (wsum > 0.0_sp) then
        rc_cor = rc_cor / wsum
        rs_cor = rs_cor / wsum
        ic_cor = ic_cor / wsum
        is_cor = is_cor / wsum
      else
        rc_cor = 0.0_sp
        rs_cor = 0.0_sp
        ic_cor = 0.0_sp
        is_cor = 0.0_sp
      end if

      ryw_tmp = rc_cor - is_cor
      iyw_tmp = rs_cor + ic_cor
      p_ex(i) = sqrt(ryw_tmp**2 + iyw_tmp**2)
      phi_ex(i) = atan2(iyw_tmp, ryw_tmp)
    end do
    !$omp end parallel do

  end subroutine extract_general_w

  subroutine extract_general_ri_w(ryt_in, iyt_in, wts_in, npts, nout, re_ex, im_ex, &
                                  cos_arr, sin_arr, maxout, maxpts, mean_rem)
    !! Weighted RM extraction (RI mode).
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, mean_rem
    real(sp), intent(in) :: ryt_in(*), iyt_in(*), wts_in(*)
    real(sp), intent(out) :: re_ex(*), im_ex(*)
    real(sp), intent(in) :: cos_arr(maxpts, maxout), sin_arr(maxpts, maxout)

    real(sp) :: ryt(npts), iyt(npts), wts(npts)
    real(sp) :: rc_cor, ic_cor, rs_cor, is_cor, ryw_tmp, iyw_tmp
    real(sp) :: wsum, mean_q, mean_u
    integer(int32) :: i, kk

    do kk = 1, npts
      wts(kk) = max(0.0_sp, wts_in(kk))
    end do

    wsum = 0.0_sp
    !$omp simd reduction(+:wsum)
    do kk = 1, npts
      wsum = wsum + wts(kk)
    end do

    if (mean_rem > 0 .and. wsum > 0.0_sp) then
      mean_q = 0.0_sp
      mean_u = 0.0_sp
      !$omp simd reduction(+:mean_q,mean_u)
      do kk = 1, npts
        mean_q = mean_q + wts(kk) * ryt_in(kk)
        mean_u = mean_u + wts(kk) * iyt_in(kk)
      end do
      mean_q = mean_q / wsum
      mean_u = mean_u / wsum
      do i = 1, npts
        ryt(i) = ryt_in(i) - mean_q
        iyt(i) = iyt_in(i) - mean_u
      end do
    else
      do i = 1, npts
        ryt(i) = ryt_in(i)
        iyt(i) = iyt_in(i)
      end do
    end if

    !$omp parallel do default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
    !$omp shared(nout,npts,ryt,iyt,wts,wsum,cos_arr,sin_arr,re_ex,im_ex)
    do i = 1, nout
      rc_cor = 0.0_sp
      rs_cor = 0.0_sp
      ic_cor = 0.0_sp
      is_cor = 0.0_sp
      !$omp simd reduction(+:rc_cor,rs_cor,ic_cor,is_cor)
      do kk = 1, npts
        rc_cor = rc_cor + wts(kk) * ryt(kk) * cos_arr(kk, i)
        rs_cor = rs_cor + wts(kk) * ryt(kk) * sin_arr(kk, i)
        ic_cor = ic_cor + wts(kk) * iyt(kk) * cos_arr(kk, i)
        is_cor = is_cor + wts(kk) * iyt(kk) * sin_arr(kk, i)
      end do

      if (wsum > 0.0_sp) then
        rc_cor = rc_cor / wsum
        rs_cor = rs_cor / wsum
        ic_cor = ic_cor / wsum
        is_cor = is_cor / wsum
      else
        rc_cor = 0.0_sp
        rs_cor = 0.0_sp
        ic_cor = 0.0_sp
        is_cor = 0.0_sp
      end if

      ryw_tmp = rc_cor - is_cor
      iyw_tmp = rs_cor + ic_cor
      re_ex(i) = ryw_tmp
      im_ex(i) = iyw_tmp
    end do
    !$omp end parallel do

  end subroutine extract_general_ri_w

  subroutine tile_extract_gpu(specQ, specU, specMask, specI, flag_arr, cos_arr, sin_arr, &
                              Q_tile, U_tile, wts_tile, ngood_tile, p_tile_arr, phi_tile_arr, &
                              mask_tile_arr, nvalid_tile_arr, &
                              nx_tile, ny_tile, nz_out, n_chan_tmpl, nrm_out, &
                              use_gpu_actual, &
                              use_input_mask, nan_check_on, rem_mean, output_mode, ap_angle_mode)
    !! CPU/GPU extraction kernel over a flat pixel index.
    !!
    !! Weight-based channel handling:
    !!   - Globally-bad channels (flag_arr==0) are skipped — they have no template entry.
    !!   - Per-pixel bad channels (NaN or mask cube) get wts_tile=0 and contribute zero
    !!     to the DFT sums. The template dimension n_chan_tmpl is therefore FIXED for all
    !!     pixels and no dense packing is required.
    !!
    !! cos/sin templates are explicit-shape (n_chan_tmpl, nrm_out), matching the storage
    !! layout written by extract_general_setup (channel as leading dimension).
    !!
    !! ngood_tile(ipix) stores n_chan_tmpl (fixed) so the calling code can size line-cut
    !! writes; nvalid_tile_arr stores the per-pixel valid-channel count.
    !!
    !! specI is accepted but currently unused (reserved for future bias removal).
    implicit none

    integer(int32), intent(in) :: nx_tile, ny_tile, nz_out, n_chan_tmpl, nrm_out
    integer(int32), intent(in) :: rem_mean, output_mode, ap_angle_mode
    logical, intent(in) :: use_gpu_actual
    logical, intent(in) :: use_input_mask, nan_check_on

    real(sp), intent(in) :: specQ(nx_tile*ny_tile*nz_out)
    real(sp), intent(in) :: specU(nx_tile*ny_tile*nz_out)
    real(sp), intent(in) :: specMask(nx_tile*ny_tile*nz_out)
    real(sp), intent(in) :: specI(nx_tile*ny_tile*nz_out)
    integer(int32), intent(in) :: flag_arr(nz_out)
    real(sp), intent(in) :: cos_arr(n_chan_tmpl, nrm_out)
    real(sp), intent(in) :: sin_arr(n_chan_tmpl, nrm_out)

    real(sp), intent(inout) :: Q_tile(nx_tile*ny_tile*nz_out)
    real(sp), intent(inout) :: U_tile(nx_tile*ny_tile*nz_out)
    real(sp), intent(inout) :: wts_tile(nx_tile*ny_tile*nz_out)
    integer(int32), intent(inout) :: ngood_tile(nx_tile*ny_tile)
    real(sp), intent(inout) :: p_tile_arr(nx_tile*ny_tile*nrm_out)
    real(sp), intent(inout) :: phi_tile_arr(nx_tile*ny_tile*nrm_out)
    integer(int8), intent(inout) :: mask_tile_arr(nx_tile*ny_tile*nz_out)
    integer(int16), intent(inout) :: nvalid_tile_arr(nx_tile*ny_tile)

    integer(int32) :: ipix, ix_loc, iy_loc, pix_base
    integer(int32) :: i, iz, cnt2, kk, tmp_index
    integer(int32) :: nvalid_pix
    logical :: per_pix_valid
    real(sp) :: mask_val, q_now, u_now
    real(sp) :: rc_cor, rs_cor, ic_cor, is_cor
    real(sp) :: ryw_tmp, iyw_tmp, wsum, mean_q, mean_u
    real(sp) :: q_eff, u_eff
    real(sp) :: p_ex(nrm_out), phi_ex(nrm_out)

#ifdef USE_GPU
    !$omp target teams distribute parallel do if(use_gpu_actual)     &
    !$omp     map(to: specQ, specU, specMask, specI, flag_arr,        &
    !$omp             cos_arr, sin_arr,                               &
    !$omp             use_input_mask, nan_check_on,                   &
    !$omp             rem_mean, output_mode, ap_angle_mode,           &
    !$omp             nx_tile, ny_tile, nz_out, n_chan_tmpl, nrm_out) &
    !$omp     map(tofrom: Q_tile, U_tile, wts_tile, ngood_tile,       &
    !$omp                 p_tile_arr, phi_tile_arr,                   &
    !$omp                 mask_tile_arr, nvalid_tile_arr)              &
    !$omp     private(ipix, ix_loc, iy_loc, pix_base,                 &
    !$omp             cnt2, iz, kk, tmp_index,                        &
    !$omp             q_now, u_now, per_pix_valid, mask_val,          &
    !$omp             nvalid_pix, wsum, mean_q, mean_u,               &
    !$omp             rc_cor, rs_cor, ic_cor, is_cor,                 &
    !$omp             ryw_tmp, iyw_tmp, q_eff, u_eff,                 &
    !$omp             p_ex, phi_ex, i)
#else
    !$omp parallel do default(none)                                   &
    !$omp     private(ipix, ix_loc, iy_loc, pix_base,                 &
    !$omp             cnt2, iz, kk, tmp_index,                        &
    !$omp             q_now, u_now, per_pix_valid, mask_val,          &
    !$omp             nvalid_pix, wsum, mean_q, mean_u,               &
    !$omp             rc_cor, rs_cor, ic_cor, is_cor,                 &
    !$omp             ryw_tmp, iyw_tmp, q_eff, u_eff,                 &
    !$omp             p_ex, phi_ex, i)                                &
    !$omp     shared(nx_tile, ny_tile, nz_out, n_chan_tmpl, nrm_out,  &
    !$omp            specQ, specU, specMask, specI, flag_arr,          &
    !$omp            cos_arr, sin_arr, Q_tile, U_tile, wts_tile,      &
    !$omp            ngood_tile, p_tile_arr, phi_tile_arr,             &
    !$omp            mask_tile_arr, nvalid_tile_arr,                   &
    !$omp            use_input_mask, nan_check_on,                     &
    !$omp            rem_mean, output_mode, ap_angle_mode)
#endif
    do ipix = 1, nx_tile*ny_tile
      iy_loc = ((ipix - 1) / nx_tile) + 1
      ix_loc = ipix - (iy_loc - 1) * nx_tile
      pix_base = (ipix - 1) * nz_out

      ! --- Channel selection: fixed-length weight-based loop ---
      ! Globally-bad channels (flag_arr==0) are skipped; their mask entry is set to 0.
      ! Per-pixel bad channels (NaN/mask) are included in the tile but with wts=0.
      ! Template index kk advances only for globally-good channels — it always maps
      ! to the same lambda^2 as cos_arr(kk,i), regardless of per-pixel masking.
      kk = 0
      nvalid_pix = 0
      do cnt2 = 1, nz_out
        iz = nz_out - cnt2 + 1
        tmp_index = ix_loc + (iy_loc - 1) * nx_tile + (iz - 1) * nx_tile * ny_tile

        if (flag_arr(cnt2) == 1) then
          kk = kk + 1
          q_now = specQ(tmp_index)
          u_now = specU(tmp_index)

          per_pix_valid = .true.
          if (use_input_mask) then
            mask_val = specMask(tmp_index)
            if (mask_val <= 0.5_sp) per_pix_valid = .false.
          end if
          if (nan_check_on) then
            if (q_now /= q_now) per_pix_valid = .false.
            if (u_now /= u_now) per_pix_valid = .false.
          end if

          Q_tile(pix_base + kk) = q_now
          U_tile(pix_base + kk) = u_now
          wts_tile(pix_base + kk) = merge(1.0_sp, 0.0_sp, per_pix_valid)
          mask_tile_arr(tmp_index) = merge(1_int8, 0_int8, per_pix_valid)
          if (per_pix_valid) nvalid_pix = nvalid_pix + 1
        else
          mask_tile_arr(tmp_index) = 0_int8
        end if
      end do

      ! ngood_tile stores n_chan_tmpl (same for every pixel) for line-cut sizing.
      ! nvalid_tile_arr stores the per-pixel count of actually valid (wts>0) channels.
      ngood_tile(ipix) = n_chan_tmpl
      nvalid_tile_arr(ix_loc + (iy_loc - 1) * nx_tile) = int(nvalid_pix, kind=int16)

      ! Sum of weights — used for early exit and normalisation
      wsum = 0.0_sp
      do kk = 1, n_chan_tmpl
        wsum = wsum + wts_tile(pix_base + kk)
      end do

      if (wsum <= 0.0_sp) then
        do i = 1, nrm_out
          p_ex(i) = 0.0_sp
          phi_ex(i) = 0.0_sp
        end do
      else

#ifdef USE_GPU
        ! Hoist weighted mean once before the RM loop (branch-free inner reduction)
        mean_q = 0.0_sp
        mean_u = 0.0_sp
        if (rem_mean > 0) then
          do kk = 1, n_chan_tmpl
            mean_q = mean_q + wts_tile(pix_base + kk) * Q_tile(pix_base + kk)
            mean_u = mean_u + wts_tile(pix_base + kk) * U_tile(pix_base + kk)
          end do
          mean_q = mean_q / wsum
          mean_u = mean_u / wsum
        end if

        ! Inline weighted DFT (GPU path — no subroutine call inside target region)
        do i = 1, nrm_out
          rc_cor = 0.0_sp
          rs_cor = 0.0_sp
          ic_cor = 0.0_sp
          is_cor = 0.0_sp
          do kk = 1, n_chan_tmpl
            q_eff = Q_tile(pix_base + kk) - mean_q
            u_eff = U_tile(pix_base + kk) - mean_u
            rc_cor = rc_cor + wts_tile(pix_base + kk) * q_eff * cos_arr(kk, i)
            rs_cor = rs_cor + wts_tile(pix_base + kk) * q_eff * sin_arr(kk, i)
            ic_cor = ic_cor + wts_tile(pix_base + kk) * u_eff * cos_arr(kk, i)
            is_cor = is_cor + wts_tile(pix_base + kk) * u_eff * sin_arr(kk, i)
          end do
          rc_cor = rc_cor / wsum
          rs_cor = rs_cor / wsum
          ic_cor = ic_cor / wsum
          is_cor = is_cor / wsum

          ryw_tmp = rc_cor - is_cor
          iyw_tmp = rs_cor + ic_cor
          if (output_mode == 1) then
            p_ex(i) = ryw_tmp
            phi_ex(i) = iyw_tmp
          else
            p_ex(i) = sqrt(ryw_tmp**2 + iyw_tmp**2)
            phi_ex(i) = atan2(iyw_tmp, ryw_tmp)
            if (ap_angle_mode == 1) phi_ex(i) = 0.5_sp * phi_ex(i)
          end if
        end do
#else
        ! CPU path: delegate to existing weighted extraction subroutines.
        ! Mean removal is handled internally by extract_general_w / extract_general_ri_w.
        if (output_mode == 1) then
          call extract_general_ri_w(Q_tile(pix_base + 1), U_tile(pix_base + 1), &
                                    wts_tile(pix_base + 1), n_chan_tmpl, nrm_out, &
                                    p_ex, phi_ex, cos_arr, sin_arr, &
                                    nrm_out, n_chan_tmpl, rem_mean)
        else
          call extract_general_w(Q_tile(pix_base + 1), U_tile(pix_base + 1), &
                                 wts_tile(pix_base + 1), n_chan_tmpl, nrm_out, &
                                 p_ex, phi_ex, cos_arr, sin_arr, &
                                 nrm_out, n_chan_tmpl, rem_mean)
          if (ap_angle_mode == 1) then
            do i = 1, nrm_out
              phi_ex(i) = 0.5_sp * phi_ex(i)
            end do
          end if
        end if
#endif
      end if

      do i = 1, nrm_out
        tmp_index = ix_loc + (iy_loc - 1) * nx_tile + (i - 1) * nx_tile * ny_tile
        p_tile_arr(tmp_index) = p_ex(i)
        phi_tile_arr(tmp_index) = phi_ex(i)
      end do
    end do
#ifdef USE_GPU
    !$omp end target teams distribute parallel do
#else
    !$omp end parallel do
#endif

  end subroutine tile_extract_gpu

  subroutine compute_mean(arr, n, mean_val)
    !! Compute mean of array
    integer(int32), intent(in) :: n
    real(sp), intent(in) :: arr(n)
    real(sp), intent(out) :: mean_val
    integer(int32) :: i
    mean_val = 0.0_sp
    do i = 1, n
      mean_val = mean_val + arr(i)
    end do
    mean_val = mean_val / real(n, sp)
  end subroutine compute_mean

  subroutine dot_product_custom(a, b, result, n)
    !! Compute dot product of two vectors
    integer(int32), intent(in) :: n
    real(sp), intent(in) :: a(n), b(n)
    real(sp), intent(out) :: result
    integer(int32) :: i
    result = 0.0_sp
    do i = 1, n
      result = result + a(i) * b(i)
    end do
  end subroutine dot_product_custom

  subroutine linspace(base, limit, n, v)
    !! Generate linearly spaced vector from base to limit
    !! Generates n points including both endpoints
    !! If n=1, v(1) = limit. If base=limit, all elements = limit
    real(sp), intent(in) :: base, limit
    integer(int32), intent(inout) :: n
    real(sp), intent(out) :: v(*)
    integer(int32) :: i
    real(sp) :: h
    
    if (n < 1) then
      write(*, '(A)') '----------------- WARNING -------------------'
      write(*, '(A)') '----------- SUBROUTINE "LINSPACE"------------'
      write(*, '(A)') '    Wrong vector length, N changed to 100'
      write(*, '(A)') '---------------------------------------------'
      n = 100
    end if
    
    if (n == 1) then
      v(1) = limit
    else if (abs(base - limit) < tiny(1.0_sp)) then
      do i = 1, n
        v(i) = limit
      end do
    else
      h = (limit - base) / real(n - 1, sp)
      do i = 1, n
        v(i) = base + real(i - 1, sp) * h
      end do
    end if
  end subroutine linspace

  function nchar(string) result(ipos)
    !! Find length of string excluding trailing whitespace
    !! Returns position of last non-whitespace character
    !! Used for trimming: string(1:nchar(string))
    character(len=*), intent(in) :: string
    integer(int32) :: ipos
    character :: c, blank, tab, null_char
    integer(int32) :: i
    
    blank = ' '
    tab = char(9)
    null_char = char(0)
    
    ipos = 0
    i = len(string)
    do while (i > 0 .and. ipos == 0)
      c = string(i:i)
      if (c /= blank .and. c /= tab .and. c /= null_char) then
        ipos = i
      end if
      i = i - 1
    end do
  end function nchar

  subroutine read_cfg_keyval(cfgfile, path, infileQ, infileU, outfile, &
                             remove_badchan, badchan_file, subim, subim_parfile, &
                             subim_ra_blc, subim_ra_trc, subim_ra_inc, &
                             subim_dec_blc, subim_dec_trc, subim_dec_inc, &
                             subim_chan_blc, subim_chan_trc, subim_chan_inc, &
                             tile_ra, tile_dec, tile_mem_frac, tile_auto, dry_run, &
                             rem_mean, remove_qu_bias, resiQ, slopeQ, resiU, slopeU, &
                             path_I, infileI, ofac, fac, beg_rm, end_rm, nrm_out_par, &
                             use_auto_rm_range, output_mode, &
                             ap_angle_mode, mask_cube_file, &
                             mask_input_cube_file, &
                             mask_trust_mode, write_mask_output, &
                             write_nvalid_output, use_gpu, status)
    !! Read all runtime parameters from a single KEY=VALUE config file.
    implicit none
    character(len=*), intent(in) :: cfgfile
    character(len=*), intent(inout) :: path, infileQ, infileU, outfile
    character(len=*), intent(inout) :: badchan_file, subim_parfile, path_I, infileI
    character(len=*), intent(inout) :: mask_cube_file
    character(len=*), intent(inout) :: mask_input_cube_file
    character(len=*), intent(inout) :: mask_trust_mode
    logical, intent(inout) :: write_mask_output
    logical, intent(inout) :: write_nvalid_output
    logical, intent(inout) :: use_gpu
    logical, intent(inout) :: remove_badchan, subim, remove_qu_bias
    integer(int32), intent(inout) :: subim_ra_blc, subim_ra_trc, subim_ra_inc
    integer(int32), intent(inout) :: subim_dec_blc, subim_dec_trc, subim_dec_inc
    integer(int32), intent(inout) :: subim_chan_blc, subim_chan_trc, subim_chan_inc
    integer(int32), intent(inout) :: tile_ra, tile_dec
    integer(int32), intent(inout) :: rem_mean, ofac, nrm_out_par, use_auto_rm_range
    integer(int32), intent(inout) :: output_mode
    integer(int32), intent(inout) :: ap_angle_mode
    logical, intent(inout) :: tile_auto, dry_run
    real(sp), intent(inout) :: resiQ, slopeQ, resiU, slopeU, fac, beg_rm, end_rm
    real(sp), intent(inout) :: tile_mem_frac
    integer(int32), intent(out) :: status

    character(len=512) :: line, key, val, key_lc
    integer(int32) :: unit_cfg, ios, line_no, io_stat
    logical :: has_kv
    logical :: seen_path, seen_infileQ, seen_infileU, seen_outfile
    logical :: seen_remove_badchan, seen_badchan_file
    logical :: seen_subim, seen_subim_parfile
    logical :: seen_subim_ra_blc, seen_subim_ra_trc, seen_subim_ra_inc
    logical :: seen_subim_dec_blc, seen_subim_dec_trc, seen_subim_dec_inc
    logical :: seen_subim_chan_blc, seen_subim_chan_trc, seen_subim_chan_inc
    logical :: seen_tile_ra, seen_tile_dec, seen_tile_mem_frac
    logical :: seen_tile_auto, seen_dry_run
    logical :: seen_rem_mean, seen_remove_qu_bias
    logical :: seen_resiQ, seen_slopeQ, seen_resiU, seen_slopeU
    logical :: seen_path_I, seen_infileI
    logical :: seen_ofac, seen_fac, seen_beg_rm, seen_end_rm, seen_nrm_out
    logical :: seen_use_auto_rm_range
    logical :: seen_output_mode
    logical :: seen_ap_angle_mode
    logical :: seen_mask_cube_file, seen_mask_input_cube_file
    logical :: seen_mask_trust_mode
    logical :: seen_write_mask_output, seen_write_nvalid_output
    logical :: seen_use_gpu

    status = 0
    line_no = 0

    seen_path = .false.
    seen_infileQ = .false.
    seen_infileU = .false.
    seen_outfile = .false.
    seen_remove_badchan = .false.
    seen_badchan_file = .false.
    seen_subim = .false.
    seen_subim_parfile = .false.
    seen_subim_ra_blc = .false.
    seen_subim_ra_trc = .false.
    seen_subim_ra_inc = .false.
    seen_subim_dec_blc = .false.
    seen_subim_dec_trc = .false.
    seen_subim_dec_inc = .false.
    seen_subim_chan_blc = .false.
    seen_subim_chan_trc = .false.
    seen_subim_chan_inc = .false.
    seen_tile_ra = .false.
    seen_tile_dec = .false.
    seen_tile_mem_frac = .false.
    seen_tile_auto = .false.
    seen_dry_run = .false.
    seen_rem_mean = .false.
    seen_remove_qu_bias = .false.
    seen_resiQ = .false.
    seen_slopeQ = .false.
    seen_resiU = .false.
    seen_slopeU = .false.
    seen_path_I = .false.
    seen_infileI = .false.
    seen_ofac = .false.
    seen_fac = .false.
    seen_beg_rm = .false.
    seen_end_rm = .false.
    seen_nrm_out = .false.
    seen_use_auto_rm_range = .false.
    seen_output_mode = .false.
    seen_ap_angle_mode = .false.
    seen_mask_cube_file = .false.
    seen_mask_input_cube_file = .false.
    seen_mask_trust_mode = .false.
    seen_write_mask_output = .false.
    seen_write_nvalid_output = .false.
    seen_use_gpu = .false.

    ! Defaults can be overridden by the config.
    path = '../DATA/'
    infileQ = ' '
    infileU = ' '
    outfile = 'output'
    remove_badchan = .false.
    badchan_file = 'bad_channels.txt'
    subim = .false.
    subim_parfile = 'subimage.par'
    subim_ra_blc = 1
    subim_ra_trc = 0
    subim_ra_inc = 1
    subim_dec_blc = 1
    subim_dec_trc = 0
    subim_dec_inc = 1
    subim_chan_blc = 0
    subim_chan_trc = 0
    subim_chan_inc = 1
    tile_ra = 0
    tile_dec = 0
    tile_mem_frac = 0.25_sp
    tile_auto = .true.
    dry_run = .false.
    rem_mean = 0
    remove_qu_bias = .false.
    resiQ = 0.0_sp
    slopeQ = 0.0_sp
    resiU = 0.0_sp
    slopeU = 0.0_sp
    path_I = path
    infileI = ' '
    ofac = 4
    fac = 3.14159265358979_sp
    beg_rm = -50.0_sp
    end_rm = 50.0_sp
    nrm_out_par = 100
    use_auto_rm_range = 1
    output_mode = 0
    ap_angle_mode = 0
    mask_cube_file = ''
    mask_input_cube_file = ''
    mask_trust_mode = 'safe'
    write_mask_output = .true.
    write_nvalid_output = .true.
    use_gpu = .false.

    unit_cfg = 11
    open(unit_cfg, file=cfgfile, status='old', iostat=ios)
    if (ios /= 0) then
      status = ios
      return
    end if

    do
      read(unit_cfg, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line_no = line_no + 1

      call split_key_value(line, key, val, has_kv)
      if (.not. has_kv) cycle

      key_lc = trim(lower_ascii(key))

      select case (key_lc)
      case ('path')
        if (seen_path) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': path'
          status = -100
          close(unit_cfg)
          return
        end if
        seen_path = .true.
        path = trim(val)
      case ('infileq')
        if (seen_infileQ) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': infileQ'
          status = -101
          close(unit_cfg)
          return
        end if
        seen_infileQ = .true.
        infileQ = trim(val)
      case ('infileu')
        if (seen_infileU) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': infileU'
          status = -102
          close(unit_cfg)
          return
        end if
        seen_infileU = .true.
        infileU = trim(val)
      case ('outfile')
        if (seen_outfile) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': outfile'
          status = -103
          close(unit_cfg)
          return
        end if
        seen_outfile = .true.
        outfile = trim(val)
      case ('remove_badchan')
        if (seen_remove_badchan) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': remove_badchan'
          status = -104
          close(unit_cfg)
          return
        end if
        seen_remove_badchan = .true.
        remove_badchan = flag_from_value(val)
      case ('badchan_file', 'global_badchan_file')
        if (seen_badchan_file) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': global_badchan_file'
          status = -105
          close(unit_cfg)
          return
        end if
        seen_badchan_file = .true.
        badchan_file = trim(val)
      case ('subim')
        if (seen_subim) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim'
          status = -106
          close(unit_cfg)
          return
        end if
        seen_subim = .true.
        subim = flag_from_value(val)
      case ('subim_parfile')
        if (seen_subim_parfile) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_parfile'
          status = -107
          close(unit_cfg)
          return
        end if
        seen_subim_parfile = .true.
        subim_parfile = trim(val)
      case ('subim_ra_blc')
        if (seen_subim_ra_blc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_ra_blc'
          status = -171
          close(unit_cfg)
          return
        end if
        seen_subim_ra_blc = .true.
        read(val, *, iostat=io_stat) subim_ra_blc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_ra_blc at line ', line_no
          status = -171
          close(unit_cfg)
          return
        end if
        if (subim_ra_blc < 1) then
          write(*,*) 'Error: subim_ra_blc must be >= 1 at line ', line_no
          status = -171
          close(unit_cfg)
          return
        end if
      case ('subim_ra_trc')
        if (seen_subim_ra_trc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_ra_trc'
          status = -172
          close(unit_cfg)
          return
        end if
        seen_subim_ra_trc = .true.
        read(val, *, iostat=io_stat) subim_ra_trc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_ra_trc at line ', line_no
          status = -172
          close(unit_cfg)
          return
        end if
        if (subim_ra_trc > 0 .and. subim_ra_trc < subim_ra_blc) then
          write(*,*) 'Error: subim_ra_trc must be >= subim_ra_blc at line ', line_no
          status = -172
          close(unit_cfg)
          return
        end if
      case ('subim_ra_inc')
        if (seen_subim_ra_inc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_ra_inc'
          status = -173
          close(unit_cfg)
          return
        end if
        seen_subim_ra_inc = .true.
        read(val, *, iostat=io_stat) subim_ra_inc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_ra_inc at line ', line_no
          status = -173
          close(unit_cfg)
          return
        end if
        if (subim_ra_inc < 1) then
          write(*,*) 'Error: subim_ra_inc must be >= 1 at line ', line_no
          status = -173
          close(unit_cfg)
          return
        end if
      case ('subim_dec_blc')
        if (seen_subim_dec_blc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_dec_blc'
          status = -174
          close(unit_cfg)
          return
        end if
        seen_subim_dec_blc = .true.
        read(val, *, iostat=io_stat) subim_dec_blc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_dec_blc at line ', line_no
          status = -174
          close(unit_cfg)
          return
        end if
        if (subim_dec_blc < 1) then
          write(*,*) 'Error: subim_dec_blc must be >= 1 at line ', line_no
          status = -174
          close(unit_cfg)
          return
        end if
      case ('subim_dec_trc')
        if (seen_subim_dec_trc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_dec_trc'
          status = -175
          close(unit_cfg)
          return
        end if
        seen_subim_dec_trc = .true.
        read(val, *, iostat=io_stat) subim_dec_trc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_dec_trc at line ', line_no
          status = -175
          close(unit_cfg)
          return
        end if
        if (subim_dec_trc > 0 .and. subim_dec_trc < subim_dec_blc) then
          write(*,*) 'Error: subim_dec_trc must be >= subim_dec_blc at line ', line_no
          status = -175
          close(unit_cfg)
          return
        end if
      case ('subim_dec_inc')
        if (seen_subim_dec_inc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_dec_inc'
          status = -176
          close(unit_cfg)
          return
        end if
        seen_subim_dec_inc = .true.
        read(val, *, iostat=io_stat) subim_dec_inc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_dec_inc at line ', line_no
          status = -176
          close(unit_cfg)
          return
        end if
        if (subim_dec_inc < 1) then
          write(*,*) 'Error: subim_dec_inc must be >= 1 at line ', line_no
          status = -176
          close(unit_cfg)
          return
        end if
      case ('subim_chan_blc')
        if (seen_subim_chan_blc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_chan_blc'
          status = -177
          close(unit_cfg)
          return
        end if
        seen_subim_chan_blc = .true.
        read(val, *, iostat=io_stat) subim_chan_blc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_chan_blc at line ', line_no
          status = -177
          close(unit_cfg)
          return
        end if
        if (subim_chan_blc < 0) then
          write(*,*) 'Error: subim_chan_blc must be >= 0 at line ', line_no
          status = -177
          close(unit_cfg)
          return
        end if
      case ('subim_chan_trc')
        if (seen_subim_chan_trc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_chan_trc'
          status = -178
          close(unit_cfg)
          return
        end if
        seen_subim_chan_trc = .true.
        read(val, *, iostat=io_stat) subim_chan_trc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_chan_trc at line ', line_no
          status = -178
          close(unit_cfg)
          return
        end if
        if (subim_chan_trc > 0 .and. subim_chan_trc < subim_chan_blc) then
          write(*,*) 'Error: subim_chan_trc must be >= subim_chan_blc at line ', line_no
          status = -178
          close(unit_cfg)
          return
        end if
      case ('subim_chan_inc')
        if (seen_subim_chan_inc) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': subim_chan_inc'
          status = -179
          close(unit_cfg)
          return
        end if
        seen_subim_chan_inc = .true.
        read(val, *, iostat=io_stat) subim_chan_inc
        if (io_stat /= 0) then
          write(*,*) 'Error reading subim_chan_inc at line ', line_no
          status = -179
          close(unit_cfg)
          return
        end if
        if (subim_chan_inc < 1) then
          write(*,*) 'Error: subim_chan_inc must be >= 1 at line ', line_no
          status = -179
          close(unit_cfg)
          return
        end if
      case ('tile_ra')
        if (seen_tile_ra) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': tile_ra'
          status = -180
          close(unit_cfg)
          return
        end if
        seen_tile_ra = .true.
        read(val, *, iostat=io_stat) tile_ra
        if (io_stat /= 0) then
          write(*,*) 'Error reading tile_ra at line ', line_no
          status = -180
          close(unit_cfg)
          return
        end if
      case ('tile_dec')
        if (seen_tile_dec) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': tile_dec'
          status = -181
          close(unit_cfg)
          return
        end if
        seen_tile_dec = .true.
        read(val, *, iostat=io_stat) tile_dec
        if (io_stat /= 0) then
          write(*,*) 'Error reading tile_dec at line ', line_no
          status = -181
          close(unit_cfg)
          return
        end if
      case ('tile_mem_frac')
        if (seen_tile_mem_frac) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': tile_mem_frac'
          status = -182
          close(unit_cfg)
          return
        end if
        seen_tile_mem_frac = .true.
        read(val, *, iostat=io_stat) tile_mem_frac
        if (io_stat /= 0) then
          write(*,*) 'Error reading tile_mem_frac at line ', line_no
          status = -182
          close(unit_cfg)
          return
        end if
      case ('tile_auto')
        if (seen_tile_auto) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': tile_auto'
          status = -183
          close(unit_cfg)
          return
        end if
        seen_tile_auto = .true.
        tile_auto = flag_from_value(val)
      case ('dry_run')
        if (seen_dry_run) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': dry_run'
          status = -184
          close(unit_cfg)
          return
        end if
        seen_dry_run = .true.
        dry_run = flag_from_value(val)
      case ('rem_mean')
        if (seen_rem_mean) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': rem_mean'
          status = -108
          close(unit_cfg)
          return
        end if
        seen_rem_mean = .true.
        read(val, *, iostat=ios) rem_mean
        if (ios /= 0) then
          write(*,*) 'Invalid integer for rem_mean at cfg line ', line_no
          status = -109
          close(unit_cfg)
          return
        end if
      case ('remove_qu_bias')
        if (seen_remove_qu_bias) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': remove_qu_bias'
          status = -110
          close(unit_cfg)
          return
        end if
        seen_remove_qu_bias = .true.
        remove_qu_bias = flag_from_value(val)
      case ('resiq')
        if (seen_resiQ) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': resiQ'
          status = -111
          close(unit_cfg)
          return
        end if
        seen_resiQ = .true.
        read(val, *, iostat=ios) resiQ
        if (ios /= 0) then
          write(*,*) 'Invalid real for resiQ at cfg line ', line_no
          status = -112
          close(unit_cfg)
          return
        end if
      case ('slopeq')
        if (seen_slopeQ) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': slopeQ'
          status = -113
          close(unit_cfg)
          return
        end if
        seen_slopeQ = .true.
        read(val, *, iostat=ios) slopeQ
        if (ios /= 0) then
          write(*,*) 'Invalid real for slopeQ at cfg line ', line_no
          status = -114
          close(unit_cfg)
          return
        end if
      case ('resiu')
        if (seen_resiU) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': resiU'
          status = -115
          close(unit_cfg)
          return
        end if
        seen_resiU = .true.
        read(val, *, iostat=ios) resiU
        if (ios /= 0) then
          write(*,*) 'Invalid real for resiU at cfg line ', line_no
          status = -116
          close(unit_cfg)
          return
        end if
      case ('slopeu')
        if (seen_slopeU) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': slopeU'
          status = -117
          close(unit_cfg)
          return
        end if
        seen_slopeU = .true.
        read(val, *, iostat=ios) slopeU
        if (ios /= 0) then
          write(*,*) 'Invalid real for slopeU at cfg line ', line_no
          status = -118
          close(unit_cfg)
          return
        end if
      case ('path_i')
        if (seen_path_I) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': path_I'
          status = -119
          close(unit_cfg)
          return
        end if
        seen_path_I = .true.
        path_I = trim(val)
      case ('infilei')
        if (seen_infileI) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': infileI'
          status = -120
          close(unit_cfg)
          return
        end if
        seen_infileI = .true.
        infileI = trim(val)
      case ('ofac')
        if (seen_ofac) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': ofac'
          status = -121
          close(unit_cfg)
          return
        end if
        seen_ofac = .true.
        read(val, *, iostat=ios) ofac
        if (ios /= 0) then
          write(*,*) 'Invalid integer for ofac at cfg line ', line_no
          status = -122
          close(unit_cfg)
          return
        end if
      case ('fac')
        if (seen_fac) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': fac'
          status = -123
          close(unit_cfg)
          return
        end if
        seen_fac = .true.
        read(val, *, iostat=ios) fac
        if (ios /= 0) then
          write(*,*) 'Invalid real for fac at cfg line ', line_no
          status = -124
          close(unit_cfg)
          return
        end if
      case ('beg_rm')
        if (seen_beg_rm) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': beg_rm'
          status = -125
          close(unit_cfg)
          return
        end if
        seen_beg_rm = .true.
        read(val, *, iostat=ios) beg_rm
        if (ios /= 0) then
          write(*,*) 'Invalid real for beg_rm at cfg line ', line_no
          status = -126
          close(unit_cfg)
          return
        end if
      case ('end_rm', 'max_rm')
        if (seen_end_rm) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': end_rm/max_rm'
          status = -161
          close(unit_cfg)
          return
        end if
        seen_end_rm = .true.
        read(val, *, iostat=ios) end_rm
        if (ios /= 0) then
          write(*,*) 'Invalid real for end_rm at cfg line ', line_no
          status = -162
          close(unit_cfg)
          return
        end if
      case ('nrm', 'nrm_out')
        if (seen_nrm_out) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': nrm/nrm_out'
          status = -127
          close(unit_cfg)
          return
        end if
        seen_nrm_out = .true.
        read(val, *, iostat=ios) nrm_out_par
        if (ios /= 0) then
          write(*,*) 'Invalid integer for nrm at cfg line ', line_no
          status = -128
          close(unit_cfg)
          return
        end if
      case ('use_auto_rm_range')
        if (seen_use_auto_rm_range) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': use_auto_rm_range'
          status = -129
          close(unit_cfg)
          return
        end if
        seen_use_auto_rm_range = .true.
        read(val, *, iostat=ios) use_auto_rm_range
        if (ios /= 0) then
          write(*,*) 'Invalid integer for use_auto_rm_range at cfg line ', line_no
          status = -130
          close(unit_cfg)
          return
        end if
      case ('output_mode')
        if (seen_output_mode) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': output_mode'
          status = -168
          close(unit_cfg)
          return
        end if
        seen_output_mode = .true.
        select case (trim(lower_ascii(val)))
        case ('ap')
          output_mode = 0
        case ('ri')
          output_mode = 1
        case default
          write(*,*) 'Invalid output_mode at cfg line ', line_no
          write(*,*) 'Allowed values: ap, ri'
          status = -169
          close(unit_cfg)
          return
        end select
      case ('ap_angle_mode')
        if (seen_ap_angle_mode) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': ap_angle_mode'
          status = -159
          close(unit_cfg)
          return
        end if
        seen_ap_angle_mode = .true.
        select case (trim(lower_ascii(val)))
        case ('phase')
          ap_angle_mode = 0
        case ('pol')
          ap_angle_mode = 1
        case default
          write(*,*) 'Invalid ap_angle_mode at cfg line ', line_no
          write(*,*) 'Allowed values: phase, pol'
          status = -160
          close(unit_cfg)
          return
        end select
      case ('mask_cube_file')
        if (seen_mask_cube_file) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': mask_cube_file'
          status = -162
          close(unit_cfg)
          return
        end if
        seen_mask_cube_file = .true.
        mask_cube_file = trim(val)
      case ('mask_input_cube_file')
        if (seen_mask_input_cube_file) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': mask_input_cube_file'
          status = -163
          close(unit_cfg)
          return
        end if
        seen_mask_input_cube_file = .true.
        mask_input_cube_file = trim(val)
      case ('mask_trust_mode')
        if (seen_mask_trust_mode) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': mask_trust_mode'
          status = -164
          close(unit_cfg)
          return
        end if
        seen_mask_trust_mode = .true.
        select case (trim(lower_ascii(val)))
        case ('safe')
          mask_trust_mode = 'safe'
        case ('strict')
          mask_trust_mode = 'strict'
        case default
          write(*,*) 'Invalid mask_trust_mode at cfg line ', line_no
          write(*,*) 'Allowed values: safe, strict'
          status = -164
          close(unit_cfg)
          return
        end select
      case ('write_mask_output')
        if (seen_write_mask_output) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': write_mask_output'
          status = -187
          close(unit_cfg)
          return
        end if
        seen_write_mask_output = .true.
        write_mask_output = flag_from_value(val)
      case ('write_nvalid_output')
        if (seen_write_nvalid_output) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': write_nvalid_output'
          status = -188
          close(unit_cfg)
          return
        end if
        seen_write_nvalid_output = .true.
        write_nvalid_output = flag_from_value(val)
      case ('use_gpu', 'use_gpus')
        if (seen_use_gpu) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': use_gpu/use_gpus'
          status = -189
          close(unit_cfg)
          return
        end if
        seen_use_gpu = .true.
        use_gpu = flag_from_value(val)
      case default
        write(*,*) 'Unknown key in cfg at line ', line_no, ': ', trim(key)
        status = -131
        close(unit_cfg)
        return
      end select
    end do

    if (ios > 0) then
      write(*,*) 'Error while reading cfg file: ', trim(cfgfile)
      status = -132
      close(unit_cfg)
      return
    end if

    if (.not. seen_path) then
      write(*,*) 'Missing required cfg key: path'
      status = -133
    else if (.not. seen_infileQ) then
      write(*,*) 'Missing required cfg key: infileQ'
      status = -134
    else if (.not. seen_infileU) then
      write(*,*) 'Missing required cfg key: infileU'
      status = -135
    else if (.not. seen_outfile) then
      write(*,*) 'Missing required cfg key: outfile'
      status = -136
    else if (.not. seen_remove_badchan) then
      write(*,*) 'Missing required cfg key: remove_badchan'
      status = -137
    else if (.not. seen_badchan_file) then
      write(*,*) 'Missing required cfg key: global_badchan_file'
      status = -138
    else if (.not. seen_subim) then
      write(*,*) 'Missing required cfg key: subim'
      status = -139
    else if (.not. seen_rem_mean) then
      write(*,*) 'Missing required cfg key: rem_mean'
      status = -141
    else if (.not. seen_remove_qu_bias) then
      write(*,*) 'Missing required cfg key: remove_qu_bias'
      status = -142
    else if (.not. seen_resiQ) then
      write(*,*) 'Missing required cfg key: resiQ'
      status = -143
    else if (.not. seen_slopeQ) then
      write(*,*) 'Missing required cfg key: slopeQ'
      status = -144
    else if (.not. seen_resiU) then
      write(*,*) 'Missing required cfg key: resiU'
      status = -145
    else if (.not. seen_slopeU) then
      write(*,*) 'Missing required cfg key: slopeU'
      status = -146
    else if (.not. seen_ofac) then
      write(*,*) 'Missing required cfg key: ofac'
      status = -149
    else if (.not. seen_fac) then
      write(*,*) 'Missing required cfg key: fac'
      status = -150
    else if (.not. seen_use_auto_rm_range) then
      write(*,*) 'Missing required cfg key: use_auto_rm_range'
      status = -153
    end if

    if (status == 0 .and. use_auto_rm_range /= 0 .and. use_auto_rm_range /= 1) then
      write(*,*) 'Invalid use_auto_rm_range: expected 0 or 1'
      status = -154
    end if
    if (status == 0 .and. output_mode /= 0 .and. output_mode /= 1) then
      write(*,*) 'Invalid output_mode: expected ap or ri'
      status = -170
    end if
    if (status == 0 .and. ofac < 1) then
      write(*,*) 'Invalid ofac: expected >= 1'
      status = -155
    end if
    if (status == 0 .and. tile_ra < 0) then
      write(*,*) 'Invalid tile_ra: expected >= 0 (0 means auto)'
      status = -185
    end if
    if (status == 0 .and. tile_dec < 0) then
      write(*,*) 'Invalid tile_dec: expected >= 0 (0 means auto)'
      status = -186
    end if
    if (status == 0 .and. (tile_mem_frac <= 0.0_sp .or. tile_mem_frac > 0.95_sp)) then
      write(*,*) 'Invalid tile_mem_frac: expected 0 < tile_mem_frac <= 0.95'
      status = -187
    end if
    if (status == 0 .and. nrm_out_par < 1) then
      write(*,*) 'Invalid nrm: expected >= 1'
      status = -156
    end if
    if (status == 0 .and. use_auto_rm_range == 0) then
      if (.not. seen_beg_rm) then
        write(*,*) 'Missing required cfg key: beg_rm (needed for use_auto_rm_range=0)'
        status = -165
      else if (.not. seen_end_rm) then
        write(*,*) 'Missing required cfg key: end_rm (needed for use_auto_rm_range=0)'
        status = -166
      else if (.not. seen_nrm_out) then
        write(*,*) 'Missing required cfg key: nrm (needed for use_auto_rm_range=0)'
        status = -167
      else if (end_rm <= beg_rm) then
        write(*,*) 'Invalid end_rm: expected end_rm > beg_rm'
        status = -168
      end if
    end if

    ! I-cube is only needed when Q/U bias correction is enabled.
    if (status == 0 .and. remove_qu_bias) then
      if (.not. seen_path_I) then
        write(*,*) 'Missing required cfg key: path_I (needed for remove_qu_bias=1)'
        status = -157
      else if (.not. seen_infileI) then
        write(*,*) 'Missing required cfg key: infileI (needed for remove_qu_bias=1)'
        status = -158
      end if
    end if

    close(unit_cfg)
  end subroutine read_cfg_keyval

  subroutine write_runtime_estimate(report_file, npix_total, nchan_total, nchan_good, &
                                    nbad_chan, nrm_out, output_mode, tile_ra, tile_dec, &
                                    nx_out, ny_out, tile_bytes_est, tile_mem_frac, status)
    !! Write a dry-run runtime estimate table.
    implicit none
    character(len=*), intent(in) :: report_file
    integer(int64), intent(in) :: npix_total
    integer(int32), intent(in) :: nchan_total, nchan_good, nbad_chan, nrm_out, output_mode
    integer(int32), intent(in) :: tile_ra, tile_dec, nx_out, ny_out
    integer(int64), intent(in) :: tile_bytes_est
    real(sp), intent(in) :: tile_mem_frac
    integer(int32), intent(out) :: status

    integer(int32) :: unit_out, i
    integer(int64) :: tiles_x, tiles_y, total_tiles
    integer(int64) :: tile_bytes_local
    real(dp) :: pix_dp, nchan_dp, nrm_dp
    real(dp) :: flops_total, flops_kernel, flops_per_term, flops_per_rm
    real(dp) :: gflops_rates(5), hours_est(5), seconds_est(5)
    real(dp) :: mem_fracs(5), gpu_budgets(3)
    real(dp) :: tile_scale, tile_side, gpu_mem_mib, gpu_mem_bytes
    real(dp) :: bytes_per_tile_pixel, frac_use
    integer(int32) :: tile_side_i, gpu_side_i
    integer(int32) :: env_len, env_stat, io_stat
    character(len=128) :: env_gpu_mem
    character(len=16) :: mode_name

    status = 0
    unit_out = 97
    open(unit_out, file=report_file, status='replace', action='write', iostat=status)
    if (status /= 0) return

    pix_dp = real(npix_total, dp)
    nchan_dp = real(nchan_good, dp)
    nrm_dp = real(nrm_out, dp)

    ! Lower-bound arithmetic model for the current implementation.
    ! Dot products dominate, so we count 8 FLOPs per channel/RM term.
    flops_per_term = 8.0_dp
    if (output_mode == 1) then
      mode_name = 'RI'
      flops_per_rm = 4.0_dp
    else
      mode_name = 'AP'
      flops_per_rm = 12.0_dp
    end if

    flops_kernel = pix_dp * nchan_dp * nrm_dp * flops_per_term
    flops_total = flops_kernel + pix_dp * nrm_dp * flops_per_rm

    gflops_rates = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 12.0_dp]
    do i = 1, size(gflops_rates)
      seconds_est(i) = flops_total / (gflops_rates(i) * 1.0d9)
      hours_est(i) = seconds_est(i) / 3600.0_dp
    end do

    write(unit_out,'(A)') 'RM-synthesis dry-run runtime estimate'
    write(unit_out,'(A)') '-------------------------------------'
    write(unit_out,'(A,1X,I0)') 'Total pixels:', npix_total
    write(unit_out,'(A,1X,I0)') 'Total frequency channels in cube:', nchan_total
    write(unit_out,'(A,1X,I0)') 'Good frequency channels used:', nchan_good
    write(unit_out,'(A,1X,I0)') 'Explicit bad channels masked:', nbad_chan
    write(unit_out,'(A,1X,I0)') 'RM samples:', nrm_out
    write(unit_out,'(A,1X,A)') 'Output mode:', trim(mode_name)
    write(unit_out,'(A,1X,F8.3)') 'Tile memory fraction target:', tile_mem_frac
    write(unit_out,'(A,1X,I0,1X,A,1X,I0)') 'Tile size (x by y):', tile_ra, 'x', tile_dec
    write(unit_out,'(A,1X,ES16.6)') 'Tile memory (bytes):', real(tile_bytes_est,dp)
    tiles_x = (int(nx_out,kind=int64) + int(tile_ra,kind=int64) - 1_int64) / &
              int(tile_ra,kind=int64)
    tiles_y = (int(ny_out,kind=int64) + int(tile_dec,kind=int64) - 1_int64) / &
              int(tile_dec,kind=int64)
    total_tiles = tiles_x * tiles_y
    write(unit_out,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') 'Tiles across x/y/total:', &
      tiles_x, 'x', tiles_y, '=>', total_tiles
    write(unit_out,'(A)') ' '
    if (nchan_good < nchan_total) then
      write(unit_out,'(A)') 'Channel reduction comes from explicit masking or channel selection.'
    else
      write(unit_out,'(A)') 'All channels in the selected span are used.'
    end if
    write(unit_out,'(A)') ' '
    write(unit_out,'(A,1X,ES16.6)') 'Kernel FLOPs (dot products):', flops_kernel
    write(unit_out,'(A,1X,ES16.6)') 'Total estimated FLOPs:', flops_total
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'Estimated wall time at sustained throughput:'
    write(unit_out,'(A)') '   GFLOP/s        seconds          hours'
    do i = 1, size(gflops_rates)
      write(unit_out,'(F8.1,2X,F12.3,2X,F10.3)') gflops_rates(i), seconds_est(i), hours_est(i)
    end do

    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'User tiling advice from tile_mem_frac'
    write(unit_out,'(A)') '-------------------------------------'
    write(unit_out,'(A)') 'How to read this section:'
    write(unit_out,'(A)') '- mem_frac is your config tile_mem_frac target.'
    write(unit_out,'(A)') '- tile(x=y) is an equivalent square tile estimate for that mem_frac.'
    write(unit_out,'(A)') '- total tiles is how many tiles cover the full output image.'
    write(unit_out,'(A)') ' '
    write(unit_out,'(A,1X,F8.3)') 'Current tile_mem_frac:', tile_mem_frac
    write(unit_out,'(A,1X,I0,1X,A,1X,I0)') 'Current planner tile:', tile_ra, 'x', tile_dec
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') '   mem_frac   tile(x=y)   tile bytes     x-tiles   y-tiles   total'
    mem_fracs = [0.05_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp]
    do i = 1, size(mem_fracs)
      tile_scale = sqrt(mem_fracs(i)/max(1.0e-6_dp, real(tile_mem_frac,dp)))
      tile_side = tile_scale * sqrt(real(tile_ra*tile_dec,dp))
      tile_side_i = int(tile_side)
      if (tile_side_i < 16) tile_side_i = 16
      if (tile_side_i > nx_out) tile_side_i = nx_out
      if (tile_side_i > ny_out) tile_side_i = ny_out
      tile_bytes_local = int(tile_side_i,kind=int64) * int(tile_side_i,kind=int64) * &
        tile_bytes_est / max(1_int64, int(tile_ra,kind=int64) * int(tile_dec,kind=int64))
      tiles_x = (int(nx_out,kind=int64) + int(tile_side_i,kind=int64) - 1_int64) / int(tile_side_i,kind=int64)
      tiles_y = (int(ny_out,kind=int64) + int(tile_side_i,kind=int64) - 1_int64) / int(tile_side_i,kind=int64)
      total_tiles = tiles_x * tiles_y
      write(unit_out,'(F8.2,5X,I4,6X,ES11.4,2X,I3,6X,I3,6X,I5)') mem_fracs(i), tile_side_i, &
        real(tile_bytes_local,dp), int(tiles_x), int(tiles_y), int(total_tiles)
    end do
    write(unit_out,'(A)') 'Rule: tile area scales approximately with tile_mem_frac.'

    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'GPU tile advisory (square tiles)'
    write(unit_out,'(A)') '-------------------------------'
    gpu_mem_mib = 6144.0_dp
    env_gpu_mem = ' '
    env_len = 0
    env_stat = 0
    call get_environment_variable('RM_GPU_MEM_MIB', env_gpu_mem, env_len, env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_gpu_mem(1:env_len), *, iostat=io_stat) gpu_mem_mib
      if (io_stat /= 0 .or. gpu_mem_mib <= 0.0_dp) gpu_mem_mib = 6144.0_dp
    end if
    gpu_mem_bytes = gpu_mem_mib * 1024.0_dp * 1024.0_dp
    bytes_per_tile_pixel = real(tile_bytes_est,dp) / max(1.0_dp, real(tile_ra*tile_dec,dp))
    write(unit_out,'(A,1X,I0,1X,A)') 'Detected GPU memory:', int(gpu_mem_mib), 'MiB'
    write(unit_out,'(A,1X,I0)') 'Estimated bytes per tile pixel:', int(bytes_per_tile_pixel)
    write(unit_out,'(A)') 'How to read this section:'
    write(unit_out,'(A)') '- gpu_budget% is the share of on-board VRAM used by one tile.'
    write(unit_out,'(A)') '- tile(x=y) is the recommended square tile side in output pixels.'
    write(unit_out,'(A)') '- lower total means fewer tile launches and less host/device traffic.'
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'Candidate sizes at GPU memory budgets:'
    write(unit_out,'(A)') ' gpu_budget% tile(x=y)    tile bytes     x-tiles   y-tiles   total'
    gpu_budgets = [20.0_dp, 35.0_dp, 50.0_dp]
    do i = 1, size(gpu_budgets)
      tile_bytes_local = int((gpu_budgets(i)/100.0_dp) * gpu_mem_bytes, kind=int64)
      gpu_side_i = int(sqrt(real(tile_bytes_local,dp)/max(1.0_dp,bytes_per_tile_pixel)))
      if (gpu_side_i < 16) gpu_side_i = 16
      if (gpu_side_i > nx_out) gpu_side_i = nx_out
      if (gpu_side_i > ny_out) gpu_side_i = ny_out
      tile_bytes_local = int(gpu_side_i,kind=int64) * int(gpu_side_i,kind=int64) * &
        int(bytes_per_tile_pixel,kind=int64)
      tiles_x = (int(nx_out,kind=int64) + int(gpu_side_i,kind=int64) - 1_int64) / int(gpu_side_i,kind=int64)
      tiles_y = (int(ny_out,kind=int64) + int(gpu_side_i,kind=int64) - 1_int64) / int(gpu_side_i,kind=int64)
      total_tiles = tiles_x * tiles_y
      write(unit_out,'(F8.1,3X,I4,6X,ES11.4,2X,I3,6X,I3,6X,I5)') gpu_budgets(i), gpu_side_i, &
        real(tile_bytes_local,dp), int(tiles_x), int(tiles_y), int(total_tiles)
    end do
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'Optimal setup suggestion for this run:'
    frac_use = 100.0_dp * real(tile_bytes_est,dp) / max(1.0_dp,gpu_mem_bytes)
    write(unit_out,'(A,F6.2,A)') 'Current planner tile uses about ', frac_use, ' % of GPU VRAM per tile.'
    tile_bytes_local = int((35.0_dp/100.0_dp) * gpu_mem_bytes, kind=int64)
    gpu_side_i = int(sqrt(real(tile_bytes_local,dp)/max(1.0_dp,bytes_per_tile_pixel)))
    if (gpu_side_i < 16) gpu_side_i = 16
    if (gpu_side_i > nx_out) gpu_side_i = nx_out
    if (gpu_side_i > ny_out) gpu_side_i = ny_out
    write(unit_out,'(A,I0,A,I0,A)') 'Recommended GPU starting tile: ', gpu_side_i, ' x ', gpu_side_i, ' at 35.0 % VRAM budget.'
    write(unit_out,'(A)') 'Suggested cfg for a first GPU-oriented run:'
    write(unit_out,'(A)') '  tile_auto=n'
    write(unit_out,'(A,I0)') '  tile_ra=', gpu_side_i
    write(unit_out,'(A,I0)') '  tile_dec=', gpu_side_i
    write(unit_out,'(A)') 'For CPU-only runs, keeping tile_auto=y with your chosen tile_mem_frac is preferred.'

    close(unit_out)
  end subroutine write_runtime_estimate

  subroutine split_key_value(raw_line, key, val, has_kv)
    implicit none
    character(len=*), intent(in) :: raw_line
    character(len=*), intent(out) :: key, val
    logical, intent(out) :: has_kv
    character(len=len(raw_line)) :: line
    integer(int32) :: p1, p2, peq, pcut

    key = ' '
    val = ' '
    has_kv = .false.

    line = raw_line
    p1 = index(line, ';')
    p2 = index(line, '#')
    if (p1 > 0 .and. p2 > 0) then
      pcut = min(p1, p2)
    else if (p1 > 0) then
      pcut = p1
    else
      pcut = p2
    end if
    if (pcut > 0) line = line(1:pcut - 1)

    line = adjustl(line)
    if (len_trim(line) == 0) return

    peq = index(line, '=')
    if (peq <= 1) return

    key = adjustl(line(1:peq - 1))
    val = adjustl(line(peq + 1:))
    if (len_trim(key) == 0 .or. len_trim(val) == 0) return

    key = trim(key)
    val = trim(val)
    has_kv = .true.
  end subroutine split_key_value

  function lower_ascii(str) result(out)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer(int32) :: i, c

    out = str
    do i = 1, len(out)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        out(i:i) = achar(c + 32)
      end if
    end do
  end function lower_ascii

  logical function flag_from_value(val)
    implicit none
    character(len=*), intent(in) :: val
    character(len=64) :: t

    t = lower_ascii(adjustl(trim(val)))
    flag_from_value = .false.
    if (len_trim(t) == 0) return

    if (t(1:1) == '1' .or. t(1:1) == 'y' .or. t(1:1) == 't') then
      flag_from_value = .true.
    end if
  end function flag_from_value

end module rm_synthesis_mod
