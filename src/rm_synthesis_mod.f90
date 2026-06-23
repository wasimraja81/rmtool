module rm_synthesis_mod
  !! Modern Fortran module for RM-synthesis extraction routines
  !! Wraps legacy fixed-form subroutines with explicit interfaces
  !! Author: Wasim Raja (modernized 2026)
  
  use iso_fortran_env, only: sp => real32, dp => real64, int32, int64
  implicit none
  
  private
  public :: extract_general_setup, extract_general, extract_general_ri
  public :: linspace, nchar
  public :: read_cfg_keyval
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
    real(sp), intent(out) :: cos_arr(maxout, maxpts), sin_arr(maxout, maxpts)
    
    real(sp) :: freq_MHz(npts), f1, f2, Lsq1, Lsq2, dfreq
    real(sp) :: t_span, d_nu, nu_span, omega, h_tmp, phi_tmp, beg_eff, end_eff
    integer(int32) :: i, j, kk
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
        cos_arr(i, kk) = cos(phi_tmp)
        sin_arr(i, kk) = -sin(phi_tmp)
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
    real(sp), intent(in) :: cos_arr(maxout, maxpts), sin_arr(maxout, maxpts)
    
    real(sp) :: ryt(npts), iyt(npts), c_template(npts), s_template(npts)
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
    
    ! Extract using pre-computed templates
    do i = 1, nout
      do kk = 1, npts
        c_template(kk) = cos_arr(i, kk)
        s_template(kk) = sin_arr(i, kk)
      end do
      
      call dot_product_custom(ryt, c_template, rc_cor, npts)
      call dot_product_custom(ryt, s_template, rs_cor, npts)
      call dot_product_custom(iyt, c_template, ic_cor, npts)
      call dot_product_custom(iyt, s_template, is_cor, npts)
      
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
    
  end subroutine extract_general

  subroutine extract_general_ri(ryt_in, iyt_in, npts, nout, re_ex, im_ex, &
                                cos_arr, sin_arr, maxout, maxpts, mean_rem)
    !! Extract RM complex spectrum directly as REAL/IMAG outputs
    !! Avoids amplitude/phase conversion when RI mode is requested
    implicit none
    integer(int32), intent(in) :: npts, nout, maxout, maxpts, mean_rem
    real(sp), intent(in) :: ryt_in(*), iyt_in(*)
    real(sp), intent(out) :: re_ex(*), im_ex(*)
    real(sp), intent(in) :: cos_arr(maxout, maxpts), sin_arr(maxout, maxpts)

    real(sp) :: ryt(npts), iyt(npts), c_template(npts), s_template(npts)
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

    do i = 1, nout
      do kk = 1, npts
        c_template(kk) = cos_arr(i, kk)
        s_template(kk) = sin_arr(i, kk)
      end do

      call dot_product_custom(ryt, c_template, rc_cor, npts)
      call dot_product_custom(ryt, s_template, rs_cor, npts)
      call dot_product_custom(iyt, c_template, ic_cor, npts)
      call dot_product_custom(iyt, s_template, is_cor, npts)

      rc_cor = rc_cor / dble(npts)
      rs_cor = rs_cor / dble(npts)
      ic_cor = ic_cor / dble(npts)
      is_cor = is_cor / dble(npts)

      ryw_tmp = rc_cor - is_cor
      iyw_tmp = rs_cor + ic_cor
      re_ex(i) = ryw_tmp
      im_ex(i) = iyw_tmp
    end do

  end subroutine extract_general_ri

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
                             ap_angle_mode, status)
    !! Read all runtime parameters from a single KEY=VALUE config file.
    implicit none
    character(len=*), intent(in) :: cfgfile
    character(len=*), intent(inout) :: path, infileQ, infileU, outfile
    character(len=*), intent(inout) :: badchan_file, subim_parfile, path_I, infileI
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
      case ('badchan_file')
        if (seen_badchan_file) then
          write(*,*) 'Duplicate key in cfg at line ', line_no, ': badchan_file'
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
      write(*,*) 'Missing required cfg key: badchan_file'
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
