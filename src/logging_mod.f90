module logging_mod
   !! Lightweight, generic timestamped/leveled logger + named-stage timer
   !! accumulator, shared by convolve_cubes/reproject_cubes/match_cubes/
   !! rmclean_cubes -- a trimmed-down, tool-agnostic sibling of
   !! rm_synthesis_mod.f90's own init_logging/log_message/timer_* (not
   !! `use`d directly: rm_synthesis_mod.f90 is on a separate build graph,
   !! and its own logging is tightly coupled to its own fixed 15-entry
   !! STAGE_* enum and its own rmsynth_config_t). Same wire format as
   !! rm_synthesis's own log lines (ISO timestamp, [level] [stage]
   !! [tid=N] message) and the same "Timing summary (seconds)" report
   !! shape, so a user reading both tools' logs sees one consistent
   !! convention -- but stages here are named at the CALL SITE
   !! (register_stage/timer_stop take a plain string) rather than a
   !! pre-declared enum, since each of these 4 tools' own block/tile
   !! pipeline has a different shape from rm_synthesis's own tile
   !! pipeline (and from each other).
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   implicit none
   private
   public :: init_logging, log_message, log_note, level_from_name_logging
   public :: timer_start, timer_stop, timer_report_summary
   public :: timer_reset_file_stages, timer_report_file_summary
   public :: flag_from_value_logging

   integer, parameter :: LOG_ERROR = 0, LOG_WARN = 1, LOG_INFO = 2, LOG_DEBUG = 3
   integer, parameter, public :: max_stages_logging = 16

   logical, save :: logger_initialized = .false.
   logical, save :: logger_owns_unit = .false.
   logical, save :: timing_enabled_glob = .false.
   integer, save :: logger_unit = 6
   integer, save :: logger_level = LOG_INFO
   real(dp), save :: stage_totals(max_stages_logging) = 0.0_dp
   ! file_stage_totals: T22 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) --
   ! a SEPARATE accumulator from stage_totals above, always updated
   ! (never gated behind timing_enabled_glob) and reset per file via
   ! timer_reset_file_stages, so timer_report_file_summary can print a
   ! per-file INFO-level stage breakdown unconditionally, without
   ! changing stage_totals/timer_report_summary's own existing,
   ! timing_enabled-gated, whole-run behaviour at all.
   real(dp), save :: file_stage_totals(max_stages_logging) = 0.0_dp
   character(len=24), save :: stage_names(max_stages_logging) = ' '
   integer, save :: n_stages_registered = 0

contains

   real(dp) function wall_time_seconds()
      use omp_lib, only: omp_get_wtime
      wall_time_seconds = omp_get_wtime()
   end function wall_time_seconds

   function level_from_name_logging(level_name) result(lvl)
      character(len=*), intent(in) :: level_name
      integer :: lvl
      character(len=16) :: tmp
      integer :: i

      tmp = adjustl(level_name)
      do i = 1, len_trim(tmp)
         tmp(i:i) = achar(iachar(tmp(i:i)) +&
         &merge(32, 0, tmp(i:i).ge.'A'.and.tmp(i:i).le.'Z'))
      enddo
      select case (trim(tmp))
      case ('error')
         lvl = LOG_ERROR
      case ('warn', 'warning')
         lvl = LOG_WARN
      case ('debug')
         lvl = LOG_DEBUG
      case default
         lvl = LOG_INFO
      end select
   end function level_from_name_logging

   function flag_from_value_logging(val) result(flag)
      !! Same first-non-blank-char '1'/'y'/'t' convention every cfg
      !! boolean key in this project uses.
      character(len=*), intent(in) :: val
      logical :: flag
      character(len=64) :: t
      integer :: i

      t = adjustl(val)
      do i = 1, len_trim(t)
         t(i:i) = achar(iachar(t(i:i)) +&
         &merge(32, 0, t(i:i).ge.'A'.and.t(i:i).le.'Z'))
      enddo
      flag = .false.
      if (len_trim(t).eq.0) return
      if (t(1:1).eq.'1' .or. t(1:1).eq.'y' .or. t(1:1).eq.'t') flag = .true.
   end function flag_from_value_logging

   character(len=32) function iso_timestamp_local_logging()
      integer :: vals(8)
      call date_and_time(values=vals)
      write(iso_timestamp_local_logging,&
      &'(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,".",I3.3)')&
      &vals(1), vals(2), vals(3), vals(5), vals(6), vals(7), vals(8)
   end function iso_timestamp_local_logging

   integer function current_thread_id_logging()
      use omp_lib, only: omp_get_thread_num
      current_thread_id_logging = omp_get_thread_num()
   end function current_thread_id_logging

   subroutine init_logging(log_level_name, timing_enabled, log_output_file, status)
      character(len=*), intent(in) :: log_level_name
      logical, intent(in) :: timing_enabled
      character(len=*), intent(in) :: log_output_file
      integer, intent(out) :: status
      integer :: ios_local

      status = 0
      logger_level = level_from_name_logging(log_level_name)
      timing_enabled_glob = timing_enabled
      stage_totals = 0.0_dp
      n_stages_registered = 0

      if (logger_owns_unit) then
         close(logger_unit)
         logger_owns_unit = .false.
         logger_unit = 6
      endif

      if (len_trim(log_output_file).gt.0) then
         logger_unit = 98
         open(logger_unit, file=trim(log_output_file), status='unknown',&
         &position='append', action='write', iostat=ios_local)
         if (ios_local.ne.0) then
            status = ios_local
            logger_unit = 6
            return
         endif
         logger_owns_unit = .true.
      endif

      logger_initialized = .true.
   end subroutine init_logging

   subroutine log_message(level_name, stage_name, message)
      character(len=*), intent(in) :: level_name, stage_name, message
      integer :: msg_level, tid
      character(len=32) :: ts

      if (.not. logger_initialized) return
      msg_level = level_from_name_logging(level_name)
      if (msg_level.gt.logger_level) return

      ts = iso_timestamp_local_logging()
      tid = current_thread_id_logging()
      !$omp critical (logging_mod_write_lock)
      write(logger_unit, '(A," [",A,"] [",A,"] [tid=",I0,"] ",A)')&
      &trim(ts), trim(level_name), trim(stage_name), tid, trim(message)
      flush(logger_unit)
      !$omp end critical (logging_mod_write_lock)
   end subroutine log_message

   subroutine log_note(stage_name, note)
      character(len=*), intent(in) :: stage_name, note
      call log_message('debug', stage_name, note)
   end subroutine log_note

   integer function register_stage(name) result(stage_id)
      !! First call with a given name registers it; later calls with the
      !! same name return the same id -- callers just pass a plain string
      !! at each timer_stop call site, no separate enum to keep in sync.
      character(len=*), intent(in) :: name
      integer :: i

      do i = 1, n_stages_registered
         if (trim(stage_names(i)).eq.trim(name)) then
            stage_id = i
            return
         endif
      enddo
      if (n_stages_registered.ge.max_stages_logging) then
         stage_id = max_stages_logging
         return
      endif
      n_stages_registered = n_stages_registered + 1
      stage_names(n_stages_registered) = name
      stage_id = n_stages_registered
   end function register_stage

   subroutine timer_start(t0)
      real(dp), intent(out) :: t0
      t0 = wall_time_seconds()
   end subroutine timer_start

   subroutine timer_stop(stage_name, t0)
      !! Always accumulates into file_stage_totals (T22, docs/dev/
      !! MULTI_BAND_TOMOGRAPHY_PLAN.md) -- cheap (one wall-clock read +
      !! one atomic add) regardless of timing_enabled, so
      !! timer_report_file_summary always has real data to print.
      !! stage_totals (the OLD, whole-run accumulator feeding
      !! timer_report_summary) keeps its exact original,
      !! timing_enabled-gated behaviour, unchanged.
      character(len=*), intent(in) :: stage_name
      real(dp), intent(in) :: t0
      integer :: stage_id
      real(dp) :: dt

      dt = max(0.0_dp, wall_time_seconds() - t0)
      stage_id = register_stage(stage_name)
      !$omp atomic
      file_stage_totals(stage_id) = file_stage_totals(stage_id) + dt
      if (.not. timing_enabled_glob) return
      !$omp atomic
      stage_totals(stage_id) = stage_totals(stage_id) + dt
   end subroutine timer_stop

   subroutine timer_report_summary()
      integer :: i
      real(dp) :: total_t, pct
      character(len=160) :: line

      if (.not. timing_enabled_glob) return

      total_t = 0.0_dp
      do i = 1, n_stages_registered
         total_t = total_t + stage_totals(i)
      enddo

      call log_timing_line_logging('Timing summary (seconds):')
      call log_timing_line_logging('stage                     sec         pct')
      do i = 1, n_stages_registered
         pct = 0.0_dp
         if (total_t.gt.0.0_dp) pct = 100.0_dp * stage_totals(i) / total_t
         write(line, '(A24,1X,F12.3,1X,F8.2)') trim(stage_names(i)),&
         &stage_totals(i), pct
         call log_timing_line_logging(trim(line))
      enddo
   end subroutine timer_report_summary

   subroutine timer_reset_file_stages()
      !! T22 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): call once before
      !! starting a new file's own block loop, so the next
      !! timer_report_file_summary reflects just that file, not a
      !! running total across every file processed so far. Does not
      !! touch stage_totals/n_stages_registered/stage_names -- those
      !! keep accumulating across the whole run for
      !! timer_report_summary's own, separate, end-of-run report.
      file_stage_totals = 0.0_dp
   end subroutine timer_reset_file_stages

   subroutine timer_report_file_summary(file_label)
      !! T22 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): prints file_
      !! stage_totals unconditionally at INFO level (never gated behind
      !! timing_enabled, unlike timer_report_summary) -- answers "where
      !! did the time for THIS file go" without needing debug-level
      !! per-block/per-thread detail. No-op if nothing was ever timed
      !! for this file (e.g. a skipped/already-matched file).
      character(len=*), intent(in) :: file_label
      integer :: i
      real(dp) :: total_t, pct
      character(len=200) :: line

      total_t = 0.0_dp
      do i = 1, n_stages_registered
         total_t = total_t + file_stage_totals(i)
      enddo
      if (total_t.le.0.0_dp) return

      call log_timing_line_logging('Stage timing for '//trim(file_label)//' (seconds):')
      do i = 1, n_stages_registered
         pct = 0.0_dp
         if (total_t.gt.0.0_dp) pct = 100.0_dp * file_stage_totals(i) / total_t
         write(line, '(A24,1X,F12.3,1X,F8.2)') trim(stage_names(i)),&
         &file_stage_totals(i), pct
         call log_timing_line_logging(trim(line))
      enddo
   end subroutine timer_report_file_summary

   subroutine log_timing_line_logging(message)
      character(len=*), intent(in) :: message
      character(len=32) :: ts
      integer :: tid

      ts = iso_timestamp_local_logging()
      tid = current_thread_id_logging()
      !$omp critical (logging_mod_write_lock)
      write(logger_unit, '(A," [info] [timing] [tid=",I0,"] ",A)')&
      &trim(ts), tid, trim(message)
      flush(logger_unit)
      !$omp end critical (logging_mod_write_lock)
   end subroutine log_timing_line_logging

end module logging_mod
