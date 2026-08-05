module rmclean_cache_mod
   !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): mask-pattern cache
   !! eviction support for rmclean_cubes -- split out from
   !! rmclean_cubes.f90's own program-internal cache logic into a real
   !! module specifically so this logic can be unit-tested standalone
   !! (a program's own internal procedures cannot be `use`d by a
   !! separate test program; only a module's own public procedures
   !! can). Pure Fortran -- no CFITSIO/OMP/config dependency.
   implicit none
   private
   public :: fnv1a_hash, linear_scan_extreme
   public :: pattern_registry_t, pattern_registry_entry_t
   public :: registry_init, registry_lookup_or_insert, registry_lookup
   public :: registry_advance, registry_next_occurrence

   type :: pattern_registry_entry_t
      !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): one distinct
      !! channel-validity pattern's own complete occurrence timeline
      !! across the WHOLE mask cube -- built once by a Pass-0 pre-scan,
      !! before real CLEAN compute starts, then consulted (never
      !! rebuilt) during the real run to answer Belady's own question:
      !! "when do I next need this pattern again?" occurrences(:) is
      !! ascending scan-position order; next_occ_ptr indexes the
      !! smallest not-yet-passed occurrence -- a property of the
      !! PATTERN's own timeline, continuously advancing through the
      !! real run regardless of whether this pattern is currently
      !! resident in the runtime cache (it can be evicted and later
      !! re-admitted; the timeline itself doesn't reset).
      integer(kind=8) :: hash = 0_8
      integer(kind=1), allocatable :: pattern(:)
      ! Deliberately integer(kind=8), not the plain (kind=4) integer
      ! used for bucket/entry COUNTS below -- scan positions/occurrence
      ! counts are a per-PIXEL running total across the whole cube, and
      ! this project has real history with exactly this class of
      ! integer-overflow bug (T10a) elsewhere. Cheap insurance (a few
      ! hundred extra MB of bookkeeping at real survey scale, not GB),
      ! not a demonstrated necessity at current scale (nx*ny for the
      ! real WALLABY+EMU cube is ~77.5M, ~27x under the int32 ceiling)
      ! -- confirmed explicitly with the user, not assumed.
      integer(kind=8), allocatable :: occurrences(:)
      integer(kind=8) :: n_occurrences = 0_8
      integer(kind=8) :: next_occ_ptr = 1_8
   end type pattern_registry_entry_t

   type :: pattern_registry_t
      !! Open-addressing hash table over pattern_registry_entry_t,
      !! same probe/hash scheme as rmclean_cubes.f90's own runtime
      !! cache -- but UNCAPPED (grows via doubling, no
      !! mask_pattern_cache_max-style ceiling); the runtime cache's own
      !! size limit is a memory-vs-reuse-benefit budget for TABLES
      !! (each ~72KiB), this registry only ever stores small integers.
      type(pattern_registry_entry_t), allocatable :: entries(:)
      integer, allocatable :: buckets(:) ! 0 = empty; else 1-based index into entries(:)
      integer :: n_entries = 0
      integer :: n_buckets = 0
      integer :: entries_capacity = 0
   end type pattern_registry_t

contains

   function linear_scan_extreme(values, n, find_max) result(idx)
      !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): 1-based index of
      !! the smallest (find_max=.false.) or largest (find_max=.true.)
      !! value in values(1:n); ties resolved to the LOWEST index
      !! (deterministic, not an arbitrary tie-break). Shared by both
      !! eviction policies: cache_eviction_policy='hitcount' scans
      !! hit_count (find_max=.false., evict the least-used); 'belady'
      !! scans each cached pattern's own next-occurrence position
      !! (find_max=.true., evict whichever is needed furthest in the
      !! future, or never again -- see rmclean_cubes.f90's own Belady
      !! eviction call site for why a plain linear scan, not a
      !! priority queue, is the deliberate choice here).
      integer(kind=8), intent(in) :: values(n)
      integer, intent(in) :: n
      logical, intent(in) :: find_max
      integer :: idx
      integer :: i
      idx = 1
      do i = 2, n
         if (find_max) then
            if (values(i).gt.values(idx)) idx = i
         else
            if (values(i).lt.values(idx)) idx = i
         endif
      end do
   end function linear_scan_extreme

   subroutine registry_init(reg)
      !! T14: allocates an empty registry with a small initial capacity
      !! -- deliberately small (16) so a unit test can cheaply exercise
      !! the growth path, not a real-scale pre-sizing hint. Real usage
      !! grows via doubling (registry_grow_entries), amortized O(1) per
      !! insert, so starting small costs a handful of extra
      !! reallocations at real survey scale, not a meaningful overhead.
      type(pattern_registry_t), intent(inout) :: reg
      integer, parameter :: initial_capacity = 16
      if (allocated(reg%entries)) deallocate(reg%entries)
      if (allocated(reg%buckets)) deallocate(reg%buckets)
      allocate(reg%entries(initial_capacity))
      reg%entries_capacity = initial_capacity
      reg%n_entries = 0
      reg%n_buckets = max(16, 4*initial_capacity)
      allocate(reg%buckets(0:reg%n_buckets-1))
      reg%buckets = 0
   end subroutine registry_init

   subroutine registry_rebuild_buckets(reg)
      !! T14: re-hashes every currently-live entry into a freshly
      !! (re)sized bucket array -- needed whenever entries_capacity
      !! changes, since bucket count is derived from it (same 4x
      !! convention as rmclean_cubes.f90's own runtime cache) and the
      !! hash-to-bucket modulo depends on bucket count. No tombstones
      !! here (unlike the runtime cache) -- the registry never evicts
      !! anything, only grows, so a plain from-scratch rebuild is
      !! always correct and never needs deletion handling.
      type(pattern_registry_t), intent(inout) :: reg
      integer(kind=8) :: probe
      integer :: i, tries, slot
      reg%n_buckets = max(16, 4*reg%entries_capacity)
      if (allocated(reg%buckets)) deallocate(reg%buckets)
      allocate(reg%buckets(0:reg%n_buckets-1))
      reg%buckets = 0
      do i = 1, reg%n_entries
         probe = modulo(reg%entries(i)%hash, int(reg%n_buckets, 8))
         do tries = 1, reg%n_buckets
            slot = int(probe)
            if (reg%buckets(slot).eq.0) exit
            probe = modulo(probe+1, int(reg%n_buckets, 8))
         enddo
         reg%buckets(slot) = i
      enddo
   end subroutine registry_rebuild_buckets

   subroutine registry_grow_entries(reg)
      !! T14: doubles entries(:) capacity (allocate-bigger + copy +
      !! move_alloc -- intrinsic assignment on a derived type with
      !! allocatable components deep-copies those components, so each
      !! entry's own pattern(:)/occurrences(:) survive the copy
      !! correctly, not just the fixed-size fields), then rebuilds the
      !! bucket table for the new capacity.
      type(pattern_registry_t), intent(inout) :: reg
      type(pattern_registry_entry_t), allocatable :: tmp(:)
      integer :: new_capacity
      new_capacity = reg%entries_capacity * 2
      allocate(tmp(new_capacity))
      tmp(1:reg%n_entries) = reg%entries(1:reg%n_entries)
      call move_alloc(tmp, reg%entries)
      reg%entries_capacity = new_capacity
      call registry_rebuild_buckets(reg)
   end subroutine registry_grow_entries

   subroutine registry_lookup_or_insert(reg, pattern, n, scan_pos, entry_id)
      !! T14: the Pass-0 pre-scan's own per-pixel call -- hash+compare
      !! (open addressing, same scheme as the runtime cache's own
      !! lookup); if this exact pattern hasn't been seen before, insert
      !! a new registry entry; either way, append scan_pos to that
      !! entry's own occurrence list (growing it on demand). Always
      !! succeeds (the registry itself is uncapped) -- entry_id is
      !! always valid on return, never a "not found" sentinel; this
      !! subroutine's whole job is to guarantee that.
      type(pattern_registry_t), intent(inout) :: reg
      integer(kind=1), intent(in) :: pattern(n)
      integer, intent(in) :: n
      integer(kind=8), intent(in) :: scan_pos
      integer, intent(out) :: entry_id
      integer(kind=8) :: h, probe
      integer :: tries, slot
      integer(kind=8), allocatable :: tmp_occ(:)
      integer(kind=8) :: new_occ_capacity

      h = fnv1a_hash(pattern, n)
      probe = modulo(h, int(reg%n_buckets, 8))
      entry_id = -1
      do tries = 1, reg%n_buckets
         slot = int(probe)
         if (reg%buckets(slot).eq.0) exit
         if (reg%entries(reg%buckets(slot))%hash.eq.h) then
            if (size(reg%entries(reg%buckets(slot))%pattern).eq.n) then
               if (all(reg%entries(reg%buckets(slot))%pattern.eq.pattern)) then
                  entry_id = reg%buckets(slot)
                  exit
               endif
            endif
         endif
         probe = modulo(probe+1, int(reg%n_buckets, 8))
      enddo

      if (entry_id.eq.-1) then
         if (reg%n_entries.ge.reg%entries_capacity) call registry_grow_entries(reg)
         reg%n_entries = reg%n_entries + 1
         entry_id = reg%n_entries
         reg%entries(entry_id)%hash = h
         if (allocated(reg%entries(entry_id)%pattern)) deallocate(reg%entries(entry_id)%pattern)
         allocate(reg%entries(entry_id)%pattern(n))
         reg%entries(entry_id)%pattern = pattern
         if (allocated(reg%entries(entry_id)%occurrences)) deallocate(reg%entries(entry_id)%occurrences)
         allocate(reg%entries(entry_id)%occurrences(4))
         reg%entries(entry_id)%n_occurrences = 0_8
         reg%entries(entry_id)%next_occ_ptr = 1_8

         ! Insert this new entry's own bucket slot -- re-probe from
         ! scratch (registry_grow_entries above may have rebuilt the
         ! whole bucket table, invalidating any earlier probe state).
         probe = modulo(h, int(reg%n_buckets, 8))
         do tries = 1, reg%n_buckets
            slot = int(probe)
            if (reg%buckets(slot).eq.0) exit
            probe = modulo(probe+1, int(reg%n_buckets, 8))
         enddo
         reg%buckets(slot) = entry_id
      endif

      if (reg%entries(entry_id)%n_occurrences.ge.size(reg%entries(entry_id)%occurrences, kind=8)) then
         new_occ_capacity = int(size(reg%entries(entry_id)%occurrences), 8) * 2_8
         allocate(tmp_occ(new_occ_capacity))
         tmp_occ(1:reg%entries(entry_id)%n_occurrences) =&
         &reg%entries(entry_id)%occurrences(1:reg%entries(entry_id)%n_occurrences)
         call move_alloc(tmp_occ, reg%entries(entry_id)%occurrences)
      endif
      reg%entries(entry_id)%n_occurrences = reg%entries(entry_id)%n_occurrences + 1_8
      reg%entries(entry_id)%occurrences(reg%entries(entry_id)%n_occurrences) = scan_pos
   end subroutine registry_lookup_or_insert

   subroutine registry_lookup(reg, pattern, n, entry_id)
      !! T14: pure READ-ONLY lookup into the registry -- unlike
      !! registry_lookup_or_insert, never inserts a new entry and never
      !! appends an occurrence. Used during Pass 1 (the real run) to
      !! find an ALREADY-registered pattern's own entry_id -- every
      !! pattern Pass 1 encounters was already recorded by Pass 0's own
      !! whole-cube pre-scan, so calling registry_lookup_or_insert again
      !! here would silently (and wrongly) append a duplicate,
      !! out-of-order occurrence to that pattern's own timeline.
      !! entry_id=0 means "not found" -- under correct operation this
      !! should never happen when cache_eviction_policy='belady' (Pass 0
      !! and Pass 1 visit the same pixels in the same order via the same
      !! next_tile_extent stepper), so callers should treat a 0 result
      !! as a genuine scan-order-mismatch bug, not a normal case to
      !! silently handle.
      type(pattern_registry_t), intent(in) :: reg
      integer(kind=1), intent(in) :: pattern(n)
      integer, intent(in) :: n
      integer, intent(out) :: entry_id
      integer(kind=8) :: h, probe
      integer :: tries, slot

      h = fnv1a_hash(pattern, n)
      probe = modulo(h, int(reg%n_buckets, 8))
      entry_id = 0
      do tries = 1, reg%n_buckets
         slot = int(probe)
         if (reg%buckets(slot).eq.0) return
         if (reg%entries(reg%buckets(slot))%hash.eq.h) then
            if (size(reg%entries(reg%buckets(slot))%pattern).eq.n) then
               if (all(reg%entries(reg%buckets(slot))%pattern.eq.pattern)) then
                  entry_id = reg%buckets(slot)
                  return
               endif
            endif
         endif
         probe = modulo(probe+1, int(reg%n_buckets, 8))
      enddo
   end subroutine registry_lookup

   subroutine registry_advance(reg, entry_id)
      !! T14: moves entry_id's own next_occ_ptr past the CURRENT scan
      !! position -- called once per pixel visit during the real run
      !! (Pass 1), for EVERY pattern encountered, hit or newly cached,
      !! since the occurrence timeline is a property of the pattern
      !! itself, not of whether it's currently resident in the runtime
      !! cache. A no-op past the end of occurrences(:) (next_occ_ptr
      !! left one past n_occurrences) -- registry_next_occurrence's own
      !! sentinel handles that case.
      type(pattern_registry_t), intent(inout) :: reg
      integer, intent(in) :: entry_id
      if (reg%entries(entry_id)%next_occ_ptr.le.reg%entries(entry_id)%n_occurrences) then
         reg%entries(entry_id)%next_occ_ptr = reg%entries(entry_id)%next_occ_ptr + 1_8
      endif
   end subroutine registry_advance

   function registry_next_occurrence(reg, entry_id) result(next_pos)
      !! T14: the scan-position this pattern will next be needed at,
      !! strictly after "now" (i.e. after the most recent
      !! registry_advance call for this entry) -- or huge(next_pos) if
      !! no future occurrence remains, a sentinel Belady's own farthest-
      !! next-use comparison (linear_scan_extreme, find_max=.true.)
      !! ranks above any real, finite position by construction.
      type(pattern_registry_t), intent(in) :: reg
      integer, intent(in) :: entry_id
      integer(kind=8) :: next_pos
      if (reg%entries(entry_id)%next_occ_ptr.gt.reg%entries(entry_id)%n_occurrences) then
         next_pos = huge(next_pos)
      else
         next_pos = reg%entries(entry_id)%occurrences(reg%entries(entry_id)%next_occ_ptr)
      endif
   end function registry_next_occurrence

   function fnv1a_hash(bytes, n) result(h)
      !! Standard 64-bit FNV-1a over a byte pattern (here: one pixel's
      !! own valid-channel mask row) -- used as the bucket key for
      !! rmclean_cubes.f90's own table_cache_entry_t, and (T14) for
      !! this module's own pattern_registry_t. Collisions are possible
      !! (any hash is), so every lookup ALSO does a full byte-for-byte
      !! pattern compare on a hash match -- never trusts the hash
      !! alone. Moved verbatim from rmclean_cubes.f90 (T14) -- no
      !! change to the algorithm itself.
      integer(kind=1), intent(in) :: bytes(n)
      integer, intent(in) :: n
      integer(kind=8) :: h
      integer(kind=8), parameter :: fnv_offset = -3750763034362895579_8
      integer(kind=8), parameter :: fnv_prime = 1099511628211_8
      integer :: i
      h = fnv_offset
      do i = 1, n
         h = ieor(h, int(bytes(i), 8))
         h = h * fnv_prime
      end do
   end function fnv1a_hash

end module rmclean_cache_mod
