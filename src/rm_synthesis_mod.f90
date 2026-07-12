module rm_synthesis_mod
  !! Modern Fortran module for RM-synthesis extraction routines
  !! Wraps legacy fixed-form subroutines with explicit interfaces
  !! Author: Wasim Raja (modernized 2026)
  
  use iso_fortran_env, only: sp => real32, dp => real64, int8, int16, int32, int64
#if defined(HOST_OMP) && (HOST_OMP == 1)
  use omp_lib, only: omp_get_wtime
#endif
  implicit none
  
  private
  public :: extract_general_setup, extract_general, extract_general_ri
  public :: extract_general_w, extract_general_ri_w
  public :: prepare_gpu_data, prepare_cpu_data, tile_extract_gpu_rm_blocked
  public :: cubestat_tail_quantile_maps
  public :: linspace, nchar
  public :: read_cfg_keyval
  public :: write_runtime_estimate
  public :: init_logging, log_message
  public :: timer_reset, timer_start, timer_stop, timer_add
  public :: timer_report_summary, wall_time_seconds
  public :: STAGE_TOTAL, STAGE_CFG_PARSE, STAGE_IO_INIT, STAGE_HEADER
  public :: STAGE_TILE_TOTAL, STAGE_TILE_READ, STAGE_TILE_MASK
  public :: STAGE_TILE_PREP, STAGE_TILE_COMPUTE, STAGE_TILE_CUBESTAT
  public :: STAGE_TILE_WRITE, STAGE_FINALIZE
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

  integer, parameter :: LOG_ERROR = 0
  integer, parameter :: LOG_WARN  = 1
  integer, parameter :: LOG_INFO  = 2
  integer, parameter :: LOG_DEBUG = 3

  integer, parameter :: STAGE_TOTAL         = 1
  integer, parameter :: STAGE_CFG_PARSE     = 2
  integer, parameter :: STAGE_IO_INIT       = 3
  integer, parameter :: STAGE_HEADER        = 4
  integer, parameter :: STAGE_TILE_TOTAL    = 5
  integer, parameter :: STAGE_TILE_READ     = 6
  integer, parameter :: STAGE_TILE_MASK     = 7
  integer, parameter :: STAGE_TILE_PREP     = 8
  integer, parameter :: STAGE_TILE_COMPUTE  = 9
  integer, parameter :: STAGE_TILE_CUBESTAT = 10
  integer, parameter :: STAGE_TILE_WRITE    = 11
  integer, parameter :: STAGE_FINALIZE      = 12
  integer, parameter :: MAX_STAGES          = 32

  logical, save :: logger_initialized = .false.
  logical, save :: logger_owns_unit = .false.
  logical, save :: timing_enabled_glob = .false.
  logical, save :: timing_tile_enabled_glob = .false.
  logical, save :: timing_io_enabled_glob = .false.
  integer, save :: logger_unit = 6
  integer, save :: logger_level = LOG_INFO
  real(dp), save :: stage_totals(MAX_STAGES) = 0.0_dp
  character(len=24), save :: stage_names(MAX_STAGES)

#if defined(HOST_OMP) && (HOST_OMP == 1)
  logical, parameter :: host_omp_enabled = .true.
#else
  logical, parameter :: host_omp_enabled = .false.
#endif
  
  public :: max_axis, max_ra, max_dec, maxchan, max_pix, maxofac, maxnt
  public :: c_velocity

contains

  subroutine init_stage_names()
    implicit none
    stage_names = ' '
    stage_names(STAGE_TOTAL) = 'total'
    stage_names(STAGE_CFG_PARSE) = 'cfg_parse'
    stage_names(STAGE_IO_INIT) = 'io_init'
    stage_names(STAGE_HEADER) = 'header_write'
    stage_names(STAGE_TILE_TOTAL) = 'tile_total'
    stage_names(STAGE_TILE_READ) = 'tile_read'
    stage_names(STAGE_TILE_MASK) = 'tile_mask'
    stage_names(STAGE_TILE_PREP) = 'tile_prep'
    stage_names(STAGE_TILE_COMPUTE) = 'tile_compute'
    stage_names(STAGE_TILE_CUBESTAT) = 'tile_cubestat'
    stage_names(STAGE_TILE_WRITE) = 'tile_write'
    stage_names(STAGE_FINALIZE) = 'finalize'
  end subroutine init_stage_names

  real(dp) function wall_time_seconds()
    implicit none
    integer(int64) :: clk_count, clk_rate
#if defined(HOST_OMP) && (HOST_OMP == 1)
    wall_time_seconds = omp_get_wtime()
#else
    call system_clock(clk_count, clk_rate)
    if (clk_rate > 0_int64) then
      wall_time_seconds = real(clk_count, dp) / real(clk_rate, dp)
    else
      wall_time_seconds = 0.0_dp
    end if
#endif
  end function wall_time_seconds

  integer function level_from_name(level_name)
    implicit none
    character(len=*), intent(in) :: level_name
    character(len=16) :: tmp
    tmp = trim(lower_ascii(level_name))
    select case (tmp)
    case ('error')
      level_from_name = LOG_ERROR
    case ('warn', 'warning')
      level_from_name = LOG_WARN
    case ('debug')
      level_from_name = LOG_DEBUG
    case default
      level_from_name = LOG_INFO
    end select
  end function level_from_name

  character(len=32) function iso_timestamp_local()
    implicit none
    integer :: vals(8)
    character(len=5) :: zone
    character(len=1) :: zsgn
    integer :: zhh, zmm

    call date_and_time(values=vals, zone=zone)
    zsgn = '+'
    if (zone(1:1) == '-') zsgn = '-'
    read(zone(2:3), '(I2)', err=10) zhh
    read(zone(4:5), '(I2)', err=10) zmm
    write(iso_timestamp_local, &
      '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,A1,I2.2,":",I2.2)') &
      vals(1), vals(2), vals(3), vals(5), vals(6), vals(7), zsgn, zhh, zmm
    return
10  continue
    write(iso_timestamp_local, &
      '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
      vals(1), vals(2), vals(3), vals(5), vals(6), vals(7)
  end function iso_timestamp_local

  subroutine init_logging(log_level_name, timing_enabled, timing_tile_enabled, &
                          timing_io_enabled, timing_output_file, status)
    implicit none
    character(len=*), intent(in) :: log_level_name
    logical, intent(in) :: timing_enabled, timing_tile_enabled
    logical, intent(in) :: timing_io_enabled
    character(len=*), intent(in) :: timing_output_file
    integer(int32), intent(out) :: status
    integer :: ios_local

    status = 0
    call init_stage_names()

    logger_level = level_from_name(log_level_name)
    timing_enabled_glob = timing_enabled
    timing_tile_enabled_glob = timing_tile_enabled
    timing_io_enabled_glob = timing_io_enabled

    if (logger_owns_unit) then
      close(logger_unit)
      logger_owns_unit = .false.
      logger_unit = 6
    end if

    if (nchar(timing_output_file) > 0) then
      logger_unit = 99
      open(logger_unit, file=trim(timing_output_file), status='unknown', &
           position='append', action='write', iostat=ios_local)
      if (ios_local /= 0) then
        status = ios_local
        logger_unit = 6
        return
      end if
      logger_owns_unit = .true.
    end if

    logger_initialized = .true.
  end subroutine init_logging

  subroutine log_message(level_name, stage_name, message)
    implicit none
    character(len=*), intent(in) :: level_name, stage_name, message
    integer :: msg_level
    character(len=32) :: ts

    if (.not. logger_initialized) return

    msg_level = level_from_name(level_name)
    if (msg_level > logger_level) return

    ts = iso_timestamp_local()
    write(logger_unit, '(A," [",A,"] [",A,"] ",A)') &
      trim(ts), trim(level_name), trim(stage_name), trim(message)
  end subroutine log_message

  subroutine timer_reset()
    implicit none
    call init_stage_names()
    stage_totals = 0.0_dp
  end subroutine timer_reset

  subroutine timer_start(t0)
    implicit none
    real(dp), intent(out) :: t0
    t0 = wall_time_seconds()
  end subroutine timer_start

  subroutine timer_add(stage_id, dt)
    implicit none
    integer(int32), intent(in) :: stage_id
    real(dp), intent(in) :: dt

    if (.not. timing_enabled_glob) return
    if (stage_id == STAGE_TILE_TOTAL .or. stage_id == STAGE_TILE_READ .or. &
        stage_id == STAGE_TILE_MASK .or. stage_id == STAGE_TILE_PREP .or. &
        stage_id == STAGE_TILE_COMPUTE .or. stage_id == STAGE_TILE_CUBESTAT .or. &
        stage_id == STAGE_TILE_WRITE) then
      if (.not. timing_tile_enabled_glob) return
    end if
    if ((stage_id == STAGE_IO_INIT .or. stage_id == STAGE_TILE_READ .or. &
         stage_id == STAGE_TILE_WRITE) .and. (.not. timing_io_enabled_glob)) return

    if (stage_id >= 1 .and. stage_id <= MAX_STAGES) then
      stage_totals(stage_id) = stage_totals(stage_id) + max(0.0_dp, dt)
    end if
  end subroutine timer_add

  subroutine timer_stop(stage_id, t0)
    implicit none
    integer(int32), intent(in) :: stage_id
    real(dp), intent(in) :: t0
    real(dp) :: dt
    dt = wall_time_seconds() - t0
    call timer_add(stage_id, dt)
  end subroutine timer_stop

  subroutine timer_report_summary()
    implicit none
    integer :: i
    real(dp) :: total_t, pct

    if (.not. timing_enabled_glob) return

    total_t = stage_totals(STAGE_TOTAL)
    if (total_t <= 0.0_dp) then
      total_t = 0.0_dp
      do i = 1, MAX_STAGES
        if (i /= STAGE_TOTAL) total_t = total_t + stage_totals(i)
      end do
    end if

    write(logger_unit, '(A)') ' '
    write(logger_unit, '(A)') 'Timing summary (seconds):'
    write(logger_unit, '(A)') 'stage                     sec         pct'
    do i = 1, MAX_STAGES
      if (len_trim(stage_names(i)) > 0 .and. stage_totals(i) > 0.0_dp) then
        pct = 0.0_dp
        if (total_t > 0.0_dp) pct = 100.0_dp * stage_totals(i) / total_t
        write(logger_unit, '(A24,1X,F12.3,1X,F8.2)') trim(stage_names(i)), &
          stage_totals(i), pct
      end if
    end do
    write(logger_unit, '(A)') ' '
  end subroutine timer_report_summary

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
    !$omp parallel do if(host_omp_enabled) default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
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

    !$omp parallel do if(host_omp_enabled) default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
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

    !$omp parallel do if(host_omp_enabled) default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
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

    !$omp parallel do if(host_omp_enabled) default(none) private(i,kk,rc_cor,rs_cor,ic_cor,is_cor,ryw_tmp,iyw_tmp) &
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

  subroutine prepare_gpu_data(specQ_flat, specU_flat, mask_tile, &
                              nx_tile, ny_tile, nz_out, &
                              specQ_gpu, specU_gpu, wts_gpu, &
                              rem_mean, mean_Q, mean_U, wsum_gpu)
    !! ========================================================================
    !! GPU Data Preparation: Reshape and Pack Arrays
    !! ========================================================================
    !!
    !! Purpose: Transform flat FITS arrays into full-size GPU-friendly layouts
    !! with unified masking applied.
    !!
    !! Input layout (from FITS):
    !!   - specQ_flat(nx_tile*ny_tile*nz_out) — flat 1D array
    !!   - specU_flat(nx_tile*ny_tile*nz_out) — flat 1D array
    !!   - mask_tile(nx_tile*ny_tile*nz_out) — unified mask: 0=bad, 1=good (integer*1)
    !!     Contains all masking: global bad channels, NaN/Inf, per-pixel mask
    !!   - Index formula: ix + (iy-1)*nx_tile + (iz-1)*nx_tile*ny_tile
    !!
    !! Output layout: specQ_gpu(npix, nz_out) — GPU-optimal (pixels fastest).
    !! Adjacent warp threads (adjacent ipix) access the same channel: stride-1
    !! across the warp = coalesced. Use prepare_cpu_data for CPU-only runs.
    !!
    !! No dense packing: all nz_out channels stored, bad channels have wts=0.
    !!
    !! wsum_gpu(npix) — per-pixel valid-channel count, always computed.
    !! Pixel-dependent masking (NaN/Inf, input mask) means wsum varies by pixel
    !! but is RM-independent. Precomputing here avoids nrm_out redundant passes
    !! inside tile_extract_gpu_rm_blocked.
    
    implicit none
    integer(int32), intent(in) :: nx_tile, ny_tile, nz_out
    integer(int32), intent(in) :: rem_mean
    
    real(sp), intent(in) :: specQ_flat(nx_tile*ny_tile*nz_out)
    real(sp), intent(in) :: specU_flat(nx_tile*ny_tile*nz_out)
    integer*1, intent(in) :: mask_tile(nx_tile*ny_tile*nz_out)
    
    real(sp), allocatable, intent(out) :: specQ_gpu(:,:)
    real(sp), allocatable, intent(out) :: specU_gpu(:,:)
    real(sp), allocatable, intent(out) :: wts_gpu(:,:)
    real(sp), allocatable, intent(out) :: mean_Q(:), mean_U(:)
    real(sp), allocatable, intent(out) :: wsum_gpu(:)
    
    integer(int32) :: npix, ipix, iz, src_idx
    real(sp) :: q_val, u_val
    real(sp) :: wsum, q_sum, u_sum
    
    npix = nx_tile * ny_tile
    
    ! GPU layout: (npix, nz_out) — pixels fastest for warp coalescing
    allocate(specQ_gpu(npix, nz_out))
    allocate(specU_gpu(npix, nz_out))
    allocate(wts_gpu(npix, nz_out))
    if (rem_mean > 0) then
      allocate(mean_Q(npix))
      allocate(mean_U(npix))
      mean_Q = 0.0_sp
      mean_U = 0.0_sp
    end if
    
    ! Load all channels using unified mask
    ! mask_tile already contains all masking info: global bad channels, NaN/Inf, per-pixel mask
    do iz = 1, nz_out
      do ipix = 1, npix
        src_idx = ipix + (iz - 1) * npix
        
        q_val = specQ_flat(src_idx)
        u_val = specU_flat(src_idx)
        
        ! Use unified mask: if mask_tile==0, channel is bad; else it's good
        ! Store all channels (bad ones have wts=0, will contribute zero to DFT)
        specQ_gpu(ipix, iz) = q_val
        specU_gpu(ipix, iz) = u_val
        wts_gpu(ipix, iz) = real(mask_tile(src_idx), sp)  ! 0.0 if bad, 1.0 if good
      end do
    end do
    
    ! Always compute per-pixel weight sums.
    ! wsum_gpu(ipix) is RM-independent but pixel-dependent when NaN/Inf or
    ! input-mask channels vary spatially. Precomputing once here saves
    ! nrm_out redundant accumulations per pixel in the GPU kernel.
    allocate(wsum_gpu(npix))
    !$omp parallel do if(host_omp_enabled) default(none) &
    !$omp     private(ipix, iz) &
    !$omp     shared(npix, nz_out, wts_gpu, wsum_gpu)
    do ipix = 1, npix
      wsum_gpu(ipix) = 0.0_sp
      do iz = 1, nz_out
        wsum_gpu(ipix) = wsum_gpu(ipix) + wts_gpu(ipix, iz)
      end do
    end do
    !$omp end parallel do

    ! Pre-compute per-pixel means if rem_mean > 0
    if (rem_mean > 0) then
      allocate(mean_Q(npix))
      allocate(mean_U(npix))
      !$omp parallel do if(host_omp_enabled) default(none) &
      !$omp     private(ipix, iz, q_sum, u_sum) &
      !$omp     shared(npix, nz_out, specQ_gpu, specU_gpu, wts_gpu, wsum_gpu, mean_Q, mean_U)
      do ipix = 1, npix
        q_sum = 0.0_sp
        u_sum = 0.0_sp
        do iz = 1, nz_out
          q_sum = q_sum + wts_gpu(ipix, iz) * specQ_gpu(ipix, iz)
          u_sum = u_sum + wts_gpu(ipix, iz) * specU_gpu(ipix, iz)
        end do
        if (wsum_gpu(ipix) > 0.0_sp) then
          mean_Q(ipix) = q_sum / wsum_gpu(ipix)
          mean_U(ipix) = u_sum / wsum_gpu(ipix)
        end if
      end do
      !$omp end parallel do
    end if

  end subroutine prepare_gpu_data

  subroutine prepare_cpu_data(specQ_flat, specU_flat, mask_tile, &
                              nx_tile, ny_tile, nz_out, &
                              specQ_cpu, specU_cpu, wts_cpu, &
                              rem_mean, mean_Q, mean_U, wsum_cpu)
    !! ========================================================================
    !! CPU Data Preparation: Reshape arrays with CPU-optimal memory layout
    !! ========================================================================
    !!
    !! Output layout: specQ_cpu(nz_out, npix) — channels fastest-varying.
    !! For the inner DFT loop (do iz=1,nz_out with ipix fixed),
    !! specQ_cpu(iz, ipix) accesses stride-1 memory. All nz_out channels
    !! for one pixel (e.g. 288*4=1152 bytes) fit in ~19 cache lines, loaded
    !! once and reused for all nrm_out RM bins.
    !!
    !! Contrast with prepare_gpu_data: (npix, nz_out) where the same loop
    !! has stride = npix*4B >> L3 cache, causing one DRAM miss per channel.
    
    implicit none
    integer(int32), intent(in) :: nx_tile, ny_tile, nz_out
    integer(int32), intent(in) :: rem_mean
    
    real(sp), intent(in) :: specQ_flat(nx_tile*ny_tile*nz_out)
    real(sp), intent(in) :: specU_flat(nx_tile*ny_tile*nz_out)
    integer*1, intent(in) :: mask_tile(nx_tile*ny_tile*nz_out)
    
    real(sp), allocatable, intent(out) :: specQ_cpu(:,:)
    real(sp), allocatable, intent(out) :: specU_cpu(:,:)
    real(sp), allocatable, intent(out) :: wts_cpu(:,:)
    real(sp), allocatable, intent(out) :: mean_Q(:), mean_U(:)
    real(sp), allocatable, intent(out) :: wsum_cpu(:)
    
    integer(int32) :: npix, ipix, iz, src_idx
    real(sp) :: q_sum, u_sum
    
    npix = nx_tile * ny_tile
    
    ! CPU layout: (nz_out, npix) — channels fastest for stride-1 inner loop
    allocate(specQ_cpu(nz_out, npix))
    allocate(specU_cpu(nz_out, npix))
    allocate(wts_cpu(nz_out, npix))
    
    ! Load all channels using unified mask
    do iz = 1, nz_out
      do ipix = 1, npix
        src_idx = ipix + (iz - 1) * npix
        specQ_cpu(iz, ipix) = specQ_flat(src_idx)
        specU_cpu(iz, ipix) = specU_flat(src_idx)
        wts_cpu(iz, ipix) = real(mask_tile(src_idx), sp)
      end do
    end do
    
    ! Per-pixel weight sums (RM-independent, precomputed once)
    allocate(wsum_cpu(npix))
    !$omp parallel do if(host_omp_enabled) default(none) &
    !$omp     private(ipix, iz) &
    !$omp     shared(npix, nz_out, wts_cpu, wsum_cpu)
    do ipix = 1, npix
      wsum_cpu(ipix) = 0.0_sp
      do iz = 1, nz_out
        wsum_cpu(ipix) = wsum_cpu(ipix) + wts_cpu(iz, ipix)
      end do
    end do
    !$omp end parallel do

    if (rem_mean > 0) then
      allocate(mean_Q(npix))
      allocate(mean_U(npix))
      !$omp parallel do if(host_omp_enabled) default(none) &
      !$omp     private(ipix, iz, q_sum, u_sum) &
      !$omp     shared(npix, nz_out, specQ_cpu, specU_cpu, wts_cpu, wsum_cpu, mean_Q, mean_U)
      do ipix = 1, npix
        q_sum = 0.0_sp
        u_sum = 0.0_sp
        do iz = 1, nz_out
          q_sum = q_sum + wts_cpu(iz, ipix) * specQ_cpu(iz, ipix)
          u_sum = u_sum + wts_cpu(iz, ipix) * specU_cpu(iz, ipix)
        end do
        if (wsum_cpu(ipix) > 0.0_sp) then
          mean_Q(ipix) = q_sum / wsum_cpu(ipix)
          mean_U(ipix) = u_sum / wsum_cpu(ipix)
        end if
      end do
      !$omp end parallel do
    end if

  end subroutine prepare_cpu_data

  subroutine tile_extract_gpu_rm_blocked(specQ_gpu, specU_gpu, wts_gpu, &
                                         mean_Q, mean_U, wsum_gpu, &
                                         cos_arr_gpu, sin_arr_gpu, &
                                         nx_tile, ny_tile, nz_out, &
                                         i_rm_block, nrm_block_now, nrm_out, &
                                         use_gpu_actual, rem_mean, output_mode, ap_angle_mode, &
                                         p_tile_arr, phi_tile_arr)
    !! ========================================================================
    !! GPU Kernel: RM-Block Tiled Extraction (Optimized)
    !! ========================================================================
    !!
    !! Purpose: Compute P(RM, pixel) and Phi(RM, pixel) using optimized GPU kernel
    !! with RM-block tiling strategy and full-size channel arrays.
    !!
    !! Data Flow:
    !!   1. Input: Pre-packed GPU arrays from prepare_gpu_data
    !!      - specQ_gpu(npix, nz_out) — FULL array with all channels
    !!      - wts_gpu(npix, nz_out) — weight mask (0 for bad, 1 for good channels)
    !!      - mean_Q(npix), mean_U(npix) — pre-computed if rem_mean > 0
    !!      - wsum_gpu(npix) — per-pixel valid-channel count (precomputed)
    !!      - cos_arr_gpu(nz_out, nrm_out) — FULL-SIZE templates
    !!   2. Process: RM-block loop (CPU), GPU parallel over pixels × RM_in_block
    !!      - Use collapse(2): parallelize (pixel, RM_in_block) pairs
    !!      - Inner loop (sequential): full DFT over all nz_out channels
    !!   3. Output: p_tile_arr(npix*nrm_out), phi_tile_arr(npix*nrm_out)
    !!
    !! Direct indexing (no channel mapping):
    !!   - Template: cos_arr_gpu(iz, i) where iz ∈ [1..nz_out]
    !!   - Data: specQ_gpu(ipix, iz) where iz ∈ [1..nz_out]
    !!   - Masking: wts_gpu(ipix, iz) handles both global and per-pixel bad channels

    implicit none
    integer(int32), intent(in) :: nx_tile, ny_tile, nz_out
    integer(int32), intent(in) :: i_rm_block, nrm_block_now, nrm_out
    integer(int32), intent(in) :: rem_mean, output_mode, ap_angle_mode
    logical, intent(in) :: use_gpu_actual
    
    real(sp), intent(in) :: specQ_gpu(:,:)
    real(sp), intent(in) :: specU_gpu(:,:)
    real(sp), intent(in) :: wts_gpu(:,:)
    real(sp), intent(in), optional :: mean_Q(nx_tile*ny_tile)
    real(sp), intent(in), optional :: mean_U(nx_tile*ny_tile)
    real(sp), intent(in) :: wsum_gpu(nx_tile*ny_tile)
    real(sp), intent(in) :: cos_arr_gpu(:, :)
    real(sp), intent(in) :: sin_arr_gpu(:, :)
    
    real(sp), intent(inout) :: p_tile_arr(:)
    real(sp), intent(inout) :: phi_tile_arr(:)
    
    integer(int32) :: ipix, npix, i_rm_local, i_rm_global, iz
    integer(int32) :: p_idx
    real(sp) :: rc_cor, rs_cor, ic_cor, is_cor, ryw_tmp, iyw_tmp
    real(sp) :: q_eff, u_eff, wt, mean_q_pix, mean_u_pix
    real(sp) :: zero_val = 0.0_sp  ! Used for runtime NaN generation (0.0/0.0)

    npix = nx_tile * ny_tile
    
#ifdef USE_GPU
    !$omp target teams distribute parallel do collapse(2) if(use_gpu_actual) &
    !$omp     map(to: specQ_gpu, specU_gpu, wts_gpu, &
    !$omp             mean_Q, mean_U, wsum_gpu, &
    !$omp             cos_arr_gpu, sin_arr_gpu, &
    !$omp             nx_tile, ny_tile, nz_out, &
    !$omp             i_rm_block, nrm_block_now, rem_mean, &
    !$omp             output_mode, ap_angle_mode, npix, zero_val) &
    !$omp     map(tofrom: p_tile_arr, phi_tile_arr) &
    !$omp     private(i_rm_local, i_rm_global, iz, &
    !$omp             rc_cor, rs_cor, ic_cor, is_cor, &
    !$omp             q_eff, u_eff, wt, ryw_tmp, iyw_tmp, &
    !$omp             mean_q_pix, mean_u_pix)
#else
    !$omp parallel do if(host_omp_enabled) collapse(2) schedule(dynamic,64) default(none) &
    !$omp     private(ipix, i_rm_local, i_rm_global, iz, p_idx, &
    !$omp             rc_cor, rs_cor, ic_cor, is_cor, &
    !$omp             q_eff, u_eff, wt, ryw_tmp, iyw_tmp, &
    !$omp             mean_q_pix, mean_u_pix) &
    !$omp     shared(npix, nz_out, nrm_block_now, &
    !$omp            specQ_gpu, specU_gpu, wts_gpu, &
    !$omp            mean_Q, mean_U, wsum_gpu, &
    !$omp            cos_arr_gpu, sin_arr_gpu, &
    !$omp            i_rm_block, rem_mean, output_mode, ap_angle_mode, &
    !$omp            p_tile_arr, phi_tile_arr, zero_val)
#endif
    do ipix = 1, npix
      do i_rm_local = 1, nrm_block_now
        i_rm_global = i_rm_block + i_rm_local - 1
        
        ! Initialize accumulators for this (pixel, RM) pair
        rc_cor = 0.0_sp
        rs_cor = 0.0_sp
        ic_cor = 0.0_sp
        is_cor = 0.0_sp
        
        ! Load per-pixel mean if needed
        if (rem_mean > 0 .and. present(mean_Q)) then
          mean_q_pix = mean_Q(ipix)
          mean_u_pix = mean_U(ipix)
        else
          mean_q_pix = 0.0_sp
          mean_u_pix = 0.0_sp
        end if
        
        ! Full DFT: sum over all channels for this RM bin
        ! Direct indexing: iz ∈ [1..nz_out] indexes both data and template
        ! wts=0 automatically skips bad channels
        do iz = 1, nz_out
#ifdef USE_GPU
          ! GPU layout: (npix, nz_out) — coalesced across warp
          q_eff = specQ_gpu(ipix, iz) - mean_q_pix
          u_eff = specU_gpu(ipix, iz) - mean_u_pix
          wt = wts_gpu(ipix, iz)
#else
          ! CPU layout: (nz_out, npix) — stride-1 channel access
          q_eff = specQ_gpu(iz, ipix) - mean_q_pix
          u_eff = specU_gpu(iz, ipix) - mean_u_pix
          wt = wts_gpu(iz, ipix)
#endif
          rc_cor = rc_cor + wt * q_eff * cos_arr_gpu(iz, i_rm_global)
          rs_cor = rs_cor + wt * q_eff * sin_arr_gpu(iz, i_rm_global)
          ic_cor = ic_cor + wt * u_eff * cos_arr_gpu(iz, i_rm_global)
          is_cor = is_cor + wt * u_eff * sin_arr_gpu(iz, i_rm_global)
        end do
        
        ! Normalize by precomputed per-pixel weight sum and compute output
        if (wsum_gpu(ipix) > 0.0_sp) then
          rc_cor = rc_cor / wsum_gpu(ipix)
          rs_cor = rs_cor / wsum_gpu(ipix)
          ic_cor = ic_cor / wsum_gpu(ipix)
          is_cor = is_cor / wsum_gpu(ipix)
          
          ryw_tmp = rc_cor - is_cor
          iyw_tmp = rs_cor + ic_cor
          
          p_idx = ipix + (i_rm_global - 1) * npix
          
          if (output_mode == 1) then
            ! Output real and imaginary parts
            p_tile_arr(p_idx) = ryw_tmp
            phi_tile_arr(p_idx) = iyw_tmp
          else
            ! Output polarized intensity and angle
            p_tile_arr(p_idx) = sqrt(ryw_tmp**2 + iyw_tmp**2)
            phi_tile_arr(p_idx) = atan2(iyw_tmp, ryw_tmp)
            if (ap_angle_mode == 1) then
              phi_tile_arr(p_idx) = 0.5_sp * phi_tile_arr(p_idx)
            end if
          end if
        else
          ! No valid data for this pixel: output NaN
          ! Using runtime 0.0/0.0 to generate IEEE NaN (portable across platforms)
          p_idx = ipix + (i_rm_global - 1) * npix
          p_tile_arr(p_idx) = zero_val / zero_val
          phi_tile_arr(p_idx) = zero_val / zero_val
        end if
      end do
    end do
#ifdef USE_GPU
    !$omp end target teams distribute parallel do
#else
    !$omp end parallel do
#endif

  end subroutine tile_extract_gpu_rm_blocked

  subroutine cubestat_tail_quantile_maps(p_tile_arr, phi_tile_arr, rm_axis, &
                                         nx_tile, ny_tile, nrm_out, &
                                         peak_map, rm_peak_map, &
                                         ang_peak_map, snr_map)
    !! Compute cubestat maps from tile RM profiles using tail-quantile sigma.
    !! Sigma definition per pixel: sigma = (q50 - q16) / 0.67449
    !! Fallback when q50<=q16: MAD-based robust sigma.
    implicit none
    integer(int32), intent(in) :: nx_tile, ny_tile, nrm_out
    real(sp), intent(in) :: p_tile_arr(nx_tile*ny_tile*nrm_out)
    real(sp), intent(in) :: phi_tile_arr(nx_tile*ny_tile*nrm_out)
    real(sp), intent(in) :: rm_axis(nrm_out)
    real(sp), intent(out) :: peak_map(nx_tile*ny_tile)
    real(sp), intent(out) :: rm_peak_map(nx_tile*ny_tile)
    real(sp), intent(out) :: ang_peak_map(nx_tile*ny_tile)
    real(sp), intent(out) :: snr_map(nx_tile*ny_tile)

    integer(int32) :: npix, ipix, irm, idx, nvalid
    integer(int32) :: idx_peak, i16, i50, i_mad
    real(sp) :: pval, pmax, sigma_noise, q16, q50, median_val, eps_sigma
    real(sp) :: vals(nrm_out), dev(nrm_out)
    real(sp) :: zero_val

    npix = nx_tile * ny_tile
    eps_sigma = 1.0e-12_sp
    zero_val = 0.0_sp

    !$omp parallel do if(host_omp_enabled) schedule(static) default(none) &
    !$omp   private(ipix, irm, idx, pval, nvalid, pmax, idx_peak, i16, i50, &
    !$omp           i_mad, q16, q50, sigma_noise, median_val, vals, dev) &
    !$omp   shared(npix, nrm_out, p_tile_arr, phi_tile_arr, rm_axis, eps_sigma, zero_val, &
    !$omp          peak_map, rm_peak_map, ang_peak_map, snr_map)
    do ipix = 1, npix
      pmax = -huge(1.0_sp)
      idx_peak = 0
      nvalid = 0

      do irm = 1, nrm_out
        idx = ipix + (irm - 1) * npix
        pval = p_tile_arr(idx)
        if (pval == pval) then
          nvalid = nvalid + 1
          vals(nvalid) = pval
          if (pval > pmax) then
            pmax = pval
            idx_peak = irm
          end if
        end if
      end do

      if (nvalid <= 0 .or. idx_peak <= 0) then
        peak_map(ipix) = zero_val / zero_val
        rm_peak_map(ipix) = zero_val / zero_val
        ang_peak_map(ipix) = zero_val / zero_val
        snr_map(ipix) = zero_val / zero_val
        cycle
      end if

      call sort_real_inplace(vals, nvalid)

      i16 = int(0.16_sp * real(nvalid - 1, sp)) + 1
      if (i16 < 1) i16 = 1
      if (i16 > nvalid) i16 = nvalid
      i50 = (nvalid + 1) / 2
      if (i50 < 1) i50 = 1
      if (i50 > nvalid) i50 = nvalid

      q16 = vals(i16)
      q50 = vals(i50)
      sigma_noise = (q50 - q16) / 0.67449_sp

      if (sigma_noise <= 0.0_sp) then
        median_val = vals(i50)
        do irm = 1, nvalid
          dev(irm) = abs(vals(irm) - median_val)
        end do
        call sort_real_inplace(dev, nvalid)
        i_mad = (nvalid + 1) / 2
        if (i_mad < 1) i_mad = 1
        if (i_mad > nvalid) i_mad = nvalid
        sigma_noise = 1.4826_sp * dev(i_mad)
      end if

      if (sigma_noise < eps_sigma) sigma_noise = eps_sigma

      idx = ipix + (idx_peak - 1) * npix
      peak_map(ipix) = pmax
      rm_peak_map(ipix) = rm_axis(idx_peak)
      ang_peak_map(ipix) = phi_tile_arr(idx)
      snr_map(ipix) = pmax / sigma_noise
    end do
    !$omp end parallel do

  contains

    subroutine sort_real_inplace(arr, n)
      real(sp), intent(inout) :: arr(n)
      integer(int32), intent(in) :: n
      integer(int32) :: i, j
      real(sp) :: key

      do i = 2, n
        key = arr(i)
        j = i - 1
        do while (j >= 1 .and. arr(j) > key)
          arr(j + 1) = arr(j)
          j = j - 1
        end do
        arr(j + 1) = key
      end do
    end subroutine sort_real_inplace

  end subroutine cubestat_tail_quantile_maps

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
                             tile_ra, tile_dec, mem_frac_ram, mem_frac_vram, &
                             gpu_vram_mib, tile_auto, dry_run, &
                             rem_mean, remove_qu_bias, resiQ, slopeQ, resiU, slopeU, &
                             path_I, infileI, ofac, fac, beg_rm, end_rm, nrm_out_par, &
                             use_auto_rm_range, output_mode, &
                             ap_angle_mode, mask_cube_file, &
                             mask_input_cube_file, &
                             mask_trust_mode, write_mask_output, &
                             write_nvalid_output, cubestat, use_gpu, io_overlap, &
                             log_level, timing_enabled, timing_tile_enabled, &
                             timing_io_enabled, timing_output_file, status)
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
    logical, intent(inout) :: cubestat
    logical, intent(inout) :: use_gpu
    logical, intent(inout) :: io_overlap
    character(len=*), intent(inout) :: log_level
    logical, intent(inout) :: timing_enabled
    logical, intent(inout) :: timing_tile_enabled
    logical, intent(inout) :: timing_io_enabled
    character(len=*), intent(inout) :: timing_output_file
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
    real(sp), intent(inout) :: mem_frac_ram, mem_frac_vram
    integer(int32), intent(inout) :: gpu_vram_mib
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
    logical :: seen_tile_ra, seen_tile_dec
    logical :: seen_mem_frac_ram, seen_mem_frac_vram, seen_gpu_vram_mib
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
    logical :: seen_cubestat
    logical :: seen_use_gpu
    logical :: seen_io_overlap
    logical :: seen_log_level
    logical :: seen_timing_enabled
    logical :: seen_timing_tile_enabled
    logical :: seen_timing_io_enabled
    logical :: seen_timing_output_file

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
    seen_mem_frac_ram = .false.
    seen_mem_frac_vram = .false.
    seen_gpu_vram_mib = .false.
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
    seen_cubestat = .false.
    seen_use_gpu = .false.
    seen_io_overlap = .false.
    seen_log_level = .false.
    seen_timing_enabled = .false.
    seen_timing_tile_enabled = .false.
    seen_timing_io_enabled = .false.
    seen_timing_output_file = .false.

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
    mem_frac_ram = 0.25_sp
    mem_frac_vram = 0.70_sp
    gpu_vram_mib = 0
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
    cubestat = .false.
    use_gpu = .false.
    io_overlap = .false.
    log_level = 'info'
    timing_enabled = .false.
    timing_tile_enabled = .false.
    timing_io_enabled = .false.
    timing_output_file = ''

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
      case ('mem_frac_ram')
        if (seen_mem_frac_ram) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': mem_frac_ram'
          status = -182
          close(unit_cfg)
          return
        end if
        seen_mem_frac_ram = .true.
        read(val, *, iostat=io_stat) mem_frac_ram
        if (io_stat /= 0) then
          write(*,*) 'Error reading mem_frac_ram at line ', line_no
          status = -182
          close(unit_cfg)
          return
        end if
      case ('mem_frac_vram')
        if (seen_mem_frac_vram) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': mem_frac_vram'
          status = -190
          close(unit_cfg)
          return
        end if
        seen_mem_frac_vram = .true.
        read(val, *, iostat=io_stat) mem_frac_vram
        if (io_stat /= 0) then
          write(*,*) 'Error reading mem_frac_vram at line ', line_no
          status = -190
          close(unit_cfg)
          return
        end if
      case ('gpu_vram_mib')
        if (seen_gpu_vram_mib) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': gpu_vram_mib'
          status = -191
          close(unit_cfg)
          return
        end if
        seen_gpu_vram_mib = .true.
        read(val, *, iostat=io_stat) gpu_vram_mib
        if (io_stat /= 0) then
          write(*,*) 'Error reading gpu_vram_mib at line ', line_no
          status = -191
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
      case ('cubestat')
        if (seen_cubestat) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': cubestat'
          status = -190
          close(unit_cfg)
          return
        end if
        seen_cubestat = .true.
        cubestat = flag_from_value(val)
      case ('use_gpu', 'use_gpus')
        if (seen_use_gpu) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': use_gpu/use_gpus'
          status = -189
          close(unit_cfg)
          return
        end if
        seen_use_gpu = .true.
        use_gpu = flag_from_value(val)
      case ('io_overlap')
        if (seen_io_overlap) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': io_overlap'
          status = -192
          close(unit_cfg)
          return
        end if
        seen_io_overlap = .true.
        io_overlap = flag_from_value(val)
      case ('log_level')
        if (seen_log_level) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': log_level'
          status = -193
          close(unit_cfg)
          return
        end if
        seen_log_level = .true.
        log_level = trim(lower_ascii(val))
      case ('timing_enabled')
        if (seen_timing_enabled) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': timing_enabled'
          status = -194
          close(unit_cfg)
          return
        end if
        seen_timing_enabled = .true.
        timing_enabled = flag_from_value(val)
      case ('timing_tile_enabled')
        if (seen_timing_tile_enabled) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': timing_tile_enabled'
          status = -195
          close(unit_cfg)
          return
        end if
        seen_timing_tile_enabled = .true.
        timing_tile_enabled = flag_from_value(val)
      case ('timing_io_enabled')
        if (seen_timing_io_enabled) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': timing_io_enabled'
          status = -196
          close(unit_cfg)
          return
        end if
        seen_timing_io_enabled = .true.
        timing_io_enabled = flag_from_value(val)
      case ('timing_output_file')
        if (seen_timing_output_file) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': timing_output_file'
          status = -197
          close(unit_cfg)
          return
        end if
        seen_timing_output_file = .true.
        timing_output_file = trim(val)
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
    if (status == 0 .and. (mem_frac_ram <= 0.0_sp .or. mem_frac_ram > 0.95_sp)) then
      write(*,*) 'Invalid mem_frac_ram: expected 0 < mem_frac_ram <= 0.95'
      status = -187
    end if
    if (status == 0 .and. (mem_frac_vram <= 0.0_sp .or. mem_frac_vram > 0.95_sp)) then
      write(*,*) 'Invalid mem_frac_vram: expected 0 < mem_frac_vram <= 0.95'
      status = -190
    end if
    if (status == 0 .and. gpu_vram_mib < 0) then
      write(*,*) 'Invalid gpu_vram_mib: expected >= 0 (0 means auto-detect)'
      status = -191
    end if
    if (status == 0) then
      if (trim(log_level) /= 'error' .and. trim(log_level) /= 'warn' .and. &
          trim(log_level) /= 'info' .and. trim(log_level) /= 'debug') then
        write(*,*) 'Invalid log_level: expected error|warn|info|debug'
        status = -198
      end if
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
                                    nx_out, ny_out, tile_bytes_est, mem_frac_ram, status)
    !! Write a dry-run runtime estimate table.
    implicit none
    character(len=*), intent(in) :: report_file
    integer(int64), intent(in) :: npix_total
    integer(int32), intent(in) :: nchan_total, nchan_good, nbad_chan, nrm_out, output_mode
    integer(int32), intent(in) :: tile_ra, tile_dec, nx_out, ny_out
    integer(int64), intent(in) :: tile_bytes_est
    real(sp), intent(in) :: mem_frac_ram
    integer(int32), intent(out) :: status

    integer(int32) :: unit_out, i
    integer(int64) :: tiles_x, tiles_y, total_tiles
    real(dp) :: pix_dp, nchan_dp, nrm_dp
    real(dp) :: flops_total, flops_kernel, flops_per_term, flops_per_rm
    real(dp) :: gflops_rates(5), hours_est(5), seconds_est(5)
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
    write(unit_out,'(A,1X,F8.3)') 'RAM memory fraction target (mem_frac_ram):', mem_frac_ram
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
    write(unit_out,'(A)') 'RAM read tiling (full-RA Dec strips)'
    write(unit_out,'(A)') '------------------------------------'
    write(unit_out,'(A)') 'The cube is read in full-RA Dec strips (RA is the'
    write(unit_out,'(A)') 'contiguous FITS axis), sized so one strip fits the'
    write(unit_out,'(A)') 'mem_frac_ram budget. These are the ACTUAL values for'
    write(unit_out,'(A)') 'this run:'
    write(unit_out,'(A,1X,F8.3)') '  mem_frac_ram:', mem_frac_ram
    write(unit_out,'(A,1X,I0,1X,A,1X,I0)') '  RAM strip (RA x Dec) px:', &
      tile_ra, 'x', tile_dec
    write(unit_out,'(A,1X,ES12.5,1X,A)') '  RAM strip size:', &
      real(tile_bytes_est,dp)/(1024.0_dp*1024.0_dp), 'MB'
    write(unit_out,'(A,1X,I0)') '  Dec strips to cover image:', int(tiles_y)
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'To use larger/smaller strips, change mem_frac_ram'
    write(unit_out,'(A)') '(fraction of total system RAM) or set tile_auto=n with'
    write(unit_out,'(A)') 'explicit tile_ra/tile_dec.'

    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'GPU memory advisory (two-level tiling)'
    write(unit_out,'(A)') '--------------------------------------'
    write(unit_out,'(A)') 'Under two-level tiling the read tile and the GPU'
    write(unit_out,'(A)') 'offload unit are DECOUPLED, so the read tile does'
    write(unit_out,'(A)') 'NOT need to fit in VRAM:'
    write(unit_out,'(A)') '  - tile_ra / tile_dec size the host RAM read block'
    write(unit_out,'(A)') '    (bigger is better for disk I/O). Control it with'
    write(unit_out,'(A)') '    mem_frac_ram (or set tile_auto=n, tile_ra/tile_dec).'
    write(unit_out,'(A)') '  - The GPU footprint is bounded SEPARATELY by the VRAM'
    write(unit_out,'(A)') '    sub-block, controlled by mem_frac_vram and gpu_vram_mib.'
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'If you hit a GPU out-of-memory (nvptx_alloc error):'
    write(unit_out,'(A)') '  - do NOT shrink tile_ra/tile_dec for VRAM;'
    write(unit_out,'(A)') '  - lower mem_frac_vram (e.g. 0.4), and/or'
    write(unit_out,'(A)') '  - set gpu_vram_mib to your card size in MiB'
    write(unit_out,'(A)') '    (cfg > GPU_MEM_MIB env > built-in default).'
    write(unit_out,'(A)') ' '
    write(unit_out,'(A)') 'See the "Two-level memory tiling (RAM -> VRAM)" section'
    write(unit_out,'(A)') 'appended below for the actual RAM-block and VRAM'
    write(unit_out,'(A)') 'sub-block sizes computed for this run.'

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
