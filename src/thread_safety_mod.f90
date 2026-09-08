module thread_safety_mod
   !! Guards against an OMP_NUM_THREADS request the host/container this
   !! program runs on cannot satisfy -- warns and clamps down instead
   !! of leaving OpenMP to silently oversubscribe (confirmed directly,
   !! outside this codebase: OMP_NUM_THREADS set above the available
   !! core count still runs with no error, every requested thread just
   !! time-sliced onto too few cores instead of one thread per core).
   !! Call check_and_clamp_omp_threads() once, at the very start of an
   !! OMP-enabled program's own main routine, before any work begins.
   !!
   !! Both /proc/cpuinfo and `lscpu` report the full HOST topology
   !! unconditionally, even inside a cpuset-restricted container
   !! (confirmed directly: --cpuset-cpus=0,1 still shows all 16 host
   !! CPUs to both) -- only a cgroup's own cpuset file reflects what
   !! this specific process can use. This module reads that file for
   !! the assigned CPU ID list, and uses `lscpu -p`
   !! only for topology STRUCTURE (which IDs share a physical core),
   !! since hardware layout itself doesn't change based on a cgroup
   !! restriction -- only which subset of it is assigned to this
   !! process does.
   implicit none
   private
   public :: check_and_clamp_omp_threads

contains

   subroutine check_and_clamp_omp_threads()
      !! No-op when compiled without -fopenmp (_OPENMP undefined) --
      !! confirmed directly: omp_set_num_threads/omp_get_max_threads
      !! are only defined in libgomp, which is not linked without
      !! -fopenmp, so calling them unconditionally would break the
      !! link for any program built without OpenMP. A serial build has
      !! no threads to clamp anyway, so skipping the whole check there
      !! is also the semantically correct behaviour, not just a build
      !! workaround.
#if defined(_OPENMP)
      use omp_lib, only: omp_set_num_threads
      use, intrinsic :: iso_fortran_env, only: error_unit
      character(len=64) :: val
      integer :: requested, stat_env, stat_read, ceiling
      logical :: avoid_hyperthreading

      ! Was a thread count requested at all? If OMP_NUM_THREADS isn't
      ! set, nothing was explicitly asked for -- whatever OpenMP's
      ! own unset-default resolves to is that runtime's own choice, not
      ! a user request to sanity-check against hardware.
      call get_environment_variable("OMP_NUM_THREADS", val, status=stat_env)
      if (stat_env /= 0) return

      read(val, *, iostat=stat_read) requested
      if (stat_read /= 0 .or. requested <= 1) return

      ! OMP_PLACES=cores is this project's own existing signal for "one
      ! thread per physical core" (set globally by scripts/
      ! run_pipeline.sh and by scratch/run_rmsynthesis_test.sh's own
      ! CPU-mode branch) -- read it back rather than inventing a
      ! separate flag for the same intent.
      call get_environment_variable("OMP_PLACES", val, status=stat_env)
      avoid_hyperthreading = (stat_env == 0 .and. trim(adjustl(val)) == "cores")

      if (avoid_hyperthreading) then
         ceiling = available_physical_cores()
      else
         ceiling = available_logical_cores()
      end if

      ! ceiling<=0 means detection itself failed (lscpu missing, no
      ! readable cgroup file, odd container setup, etc.) -- fail safe
      ! by skipping the check rather than clamping to a nonsense value.
      if (ceiling <= 0) return

      if (requested > ceiling) then
         if (avoid_hyperthreading) then
            write(error_unit, '(A,I0,A,I0,A,I0,A)') &
               'WARNING: OMP_NUM_THREADS=', requested, ' exceeds the ', &
               ceiling, ' physical cores available on this host/' // &
               'container (OMP_PLACES=cores asks for one thread per ' // &
               'physical core) -- clamping OMP_NUM_THREADS down to ', &
               ceiling, '.'
         else
            write(error_unit, '(A,I0,A,I0,A,I0,A)') &
               'WARNING: OMP_NUM_THREADS=', requested, ' exceeds the ', &
               ceiling, ' logical CPUs available on this host/' // &
               'container -- clamping OMP_NUM_THREADS down to ', &
               ceiling, '.'
         end if
         call omp_set_num_threads(ceiling)
      end if
#endif
   end subroutine check_and_clamp_omp_threads

   function available_logical_cores() result(n)
      !! Container-aware logical (SMT-inclusive) CPU count -- the
      !! length of this process's own assigned-CPU-ID list, not the
      !! host's raw total.
      integer :: n
      integer, allocatable :: ids(:)
      call assigned_cpu_ids(ids)
      n = size(ids)
   end function available_logical_cores

   function available_physical_cores() result(n)
      !! Distinct (core,socket) pairs among just the IDs this process
      !! can use -- correctly distinguishes "N CPUs that are SMT
      !! siblings of the same physical core" from "N CPUs that are
      !! different cores", using the host's own topology table
      !! restricted to just the assigned IDs.
      integer :: n
      integer, allocatable :: ids(:), topo_cpu(:), topo_core(:), topo_socket(:)
      integer, allocatable :: keys(:)
      integer :: i, j, n_topo, key, n_ids
      logical :: found

      n = 0
      call assigned_cpu_ids(ids)
      n_ids = size(ids)
      if (n_ids <= 0) return

      call host_topology(topo_cpu, topo_core, topo_socket, n_topo)
      if (n_topo <= 0) return

      allocate(keys(n_ids))
      do i = 1, n_ids
         found = .false.
         do j = 1, n_topo
            if (topo_cpu(j) == ids(i)) then
               key = topo_socket(j) * 100000 + topo_core(j)
               found = .true.
               exit
            end if
         end do
         if (.not. found) cycle   ! assigned ID missing from topology -- skip it, fail safe
         if (.not. any(keys(1:n) == key)) then
            n = n + 1
            keys(n) = key
         end if
      end do
   end function available_physical_cores

   subroutine assigned_cpu_ids(ids)
      !! This process's own usable logical CPU IDs -- read from the
      !! cgroup's own cpuset file (v2's cpuset.cpus.effective first,
      !! v1's cpuset.cpus as fallback), NOT from /proc/self/status's
      !! Cpus_allowed_list. Confirmed directly, the hard way: any
      !! binary linked against libgomp (i.e. every program in this
      !! project that uses OpenMP) has OMP_PLACES=cores narrow ITS OWN
      !! process affinity down to a single place at OpenMP's own
      !! startup -- before this subroutine's own first line runs, no
      !! matter how early it's called -- making Cpus_allowed_list
      !! measure libgomp's own side effect, not the host/container's
      !! own boundary. A cgroup's own cpuset file is a separate piece
      !! of kernel state describing the container/job's configured
      !! boundary, unaffected by what any one process inside it has
      !! bound its own affinity to. Trade-off accepted: this misses a
      !! bare `taskset`-only restriction with no container or cgroup
      !! involved at all (taskset uses the same per-process affinity
      !! mechanism libgomp's own narrowing does, so the two can't be
      !! told apart from inside the process) -- the
      !! container/Singularity/SLURM cases this project needs are all
      !! cgroup-backed.
      integer, allocatable, intent(out) :: ids(:)
      character(len=256) :: line
      integer :: unit, ios
      logical :: got_it

      got_it = .false.

      open(newunit=unit, file='/sys/fs/cgroup/cpuset.cpus.effective', &
           status='old', action='read', iostat=ios)
      if (ios == 0) then
         read(unit, '(A)', iostat=ios) line
         close(unit)
         if (ios == 0) then
            call parse_cpu_list(line, ids)
            got_it = .true.
         end if
      end if

      if (.not. got_it) then
         open(newunit=unit, file='/sys/fs/cgroup/cpuset/cpuset.cpus', &
              status='old', action='read', iostat=ios)
         if (ios == 0) then
            read(unit, '(A)', iostat=ios) line
            close(unit)
            if (ios == 0) then
               call parse_cpu_list(line, ids)
               got_it = .true.
            end if
         end if
      end if

      if (.not. got_it) then
         ! Neither cgroup file exists -- no cpuset restriction is in
         ! effect at all (bare host, no container). Fall back to every
         ! ID the host's own topology reports.
         block
            integer, allocatable :: core_ids(:), socket_ids(:)
            integer :: n
            call host_topology(ids, core_ids, socket_ids, n)
         end block
      end if
   end subroutine assigned_cpu_ids

   subroutine parse_cpu_list(str, ids)
      !! Parses a cgroup cpuset-style range list ("0-1", "0,2,4-6")
      !! into an explicit array of CPU IDs.
      character(len=*), intent(in) :: str
      integer, allocatable, intent(out) :: ids(:)
      character(len=len(str)) :: work
      integer :: i, lo, hi, ios, start_pos, comma_pos, dash_pos

      allocate(ids(0))
      work = adjustl(trim(str))
      start_pos = 1
      do while (start_pos <= len_trim(work))
         comma_pos = index(work(start_pos:), ',')
         if (comma_pos == 0) comma_pos = len_trim(work) - start_pos + 2
         block
            character(len=:), allocatable :: token
            token = trim(work(start_pos:start_pos + comma_pos - 2))
            dash_pos = index(token, '-')
            if (dash_pos > 0) then
               read(token(1:dash_pos - 1), *, iostat=ios) lo
               if (ios == 0) read(token(dash_pos + 1:), *, iostat=ios) hi
               if (ios == 0) then
                  do i = lo, hi
                     ids = [ids, i]
                  end do
               end if
            else
               read(token, *, iostat=ios) lo
               if (ios == 0) ids = [ids, lo]
            end if
         end block
         start_pos = start_pos + comma_pos
      end do
   end subroutine parse_cpu_list

   subroutine host_topology(cpu_ids, core_ids, socket_ids, n)
      !! Full HOST cpu/core/socket mapping via `lscpu -p` -- accurate
      !! for topology structure even inside a cpuset-restricted
      !! container, used only to look up which of this process's own
      !! assigned IDs (assigned_cpu_ids) share a physical core. Shells
      !! out rather than parsing /proc/cpuinfo directly -- lscpu's own
      !! -p format is a stable, purpose-built machine-readable output,
      !! not a key-value format meant for human eyes first. Capture
      !! file lands in the current working directory, never /tmp --
      !! same convention this project's own convolve_cubes.f90/
      !! match_cubes.f90/reproject_cubes.f90 already use for their own
      !! dry-run disk-type detection, since HPC compute nodes often
      !! restrict or omit /tmp entirely.
      integer, allocatable, intent(out) :: cpu_ids(:), core_ids(:), socket_ids(:)
      integer, intent(out) :: n
      character(len=600) :: shell_cmd
      character(len=64) :: capture_file
      character(len=256) :: line
      integer :: capture_unit, ios, exitstat_cmd
      integer :: p1, p2, v_cpu, v_core, v_socket, ios2

      n = 0
      allocate(cpu_ids(0), core_ids(0), socket_ids(0))

      write(capture_file, '(A,I0,A)') '.thread_safety_lscpu_', getpid_wrap(), '.tmp'
      shell_cmd = 'lscpu -p=CPU,CORE,SOCKET > ' // trim(capture_file) // ' 2>/dev/null'
      call execute_command_line(trim(shell_cmd), wait=.true., exitstat=exitstat_cmd)

      open(newunit=capture_unit, file=trim(capture_file), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         close(capture_unit, status='delete', iostat=ios)
         return
      end if

      do
         read(capture_unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         p1 = index(line, ',')
         if (p1 <= 0) cycle
         p2 = index(line(p1 + 1:), ',')
         if (p2 <= 0) cycle
         p2 = p1 + p2
         read(line(1:p1 - 1), *, iostat=ios2) v_cpu
         if (ios2 /= 0) cycle
         read(line(p1 + 1:p2 - 1), *, iostat=ios2) v_core
         if (ios2 /= 0) cycle
         read(line(p2 + 1:), *, iostat=ios2) v_socket
         if (ios2 /= 0) cycle
         cpu_ids = [cpu_ids, v_cpu]
         core_ids = [core_ids, v_core]
         socket_ids = [socket_ids, v_socket]
         n = n + 1
      end do
      close(capture_unit, status='delete')
   end subroutine host_topology

   integer function getpid_wrap() result(pid)
      !! Wraps the getpid() C library call (via iso_c_binding) purely
      !! to make the capture-file name collision-proof if multiple
      !! instances of this check somehow overlap in the same directory
      !! -- same implementation already used by convolve_cubes.f90/
      !! match_cubes.f90/reproject_cubes.f90 for the identical reason.
      use, intrinsic :: iso_c_binding, only: c_int
      interface
         function c_getpid() bind(C, name="getpid") result(r)
            import :: c_int
            integer(c_int) :: r
         end function c_getpid
      end interface
      pid = int(c_getpid())
   end function getpid_wrap

end module thread_safety_mod
