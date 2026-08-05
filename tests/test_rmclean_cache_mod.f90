program test_rmclean_cache_mod
   !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): standalone unit tests
   !! for rmclean_cache_mod.f90 -- built as this ticket's increments
   !! land (fnv1a_hash first, then the pattern_registry_t machinery).
   !! Deliberately FITS-free and FFTW-free: every check here is pure
   !! in-memory logic, no mask files, no rmclean_mod dependency at all.
   use, intrinsic :: iso_fortran_env, only: int8, int64
   use rmclean_cache_mod, only: fnv1a_hash, linear_scan_extreme,&
   &pattern_registry_t, registry_init, registry_lookup_or_insert,&
   &registry_lookup, registry_advance, registry_next_occurrence
   implicit none

   logical :: all_pass
   integer(kind=1) :: empty_bytes(0)
   integer(kind=1) :: seq5(5), seq7(7)
   integer(kind=8) :: h

   all_pass = .true.

   ! Known FNV-1a 64-bit test vectors, independently computed (Python
   ! reference implementation, not derived from this Fortran code) --
   ! see docs/dev/RMCLEAN_INTEGRATION_PLAN.md T14 for the derivation.
   ! Empty byte array: no loop iterations at all, hash must equal the
   ! offset basis constant itself -- also a direct sanity check that
   ! fnv_offset in the implementation IS the standard FNV-1a 64-bit
   ! offset basis (14695981039346656037 unsigned == this signed value).
   h = fnv1a_hash(empty_bytes, 0)
   call check(h.eq.-3750763034362895579_8, 'fnv1a_hash: empty byte array == offset basis', all_pass)

   seq5 = [integer(kind=1) :: 1, 2, 3, 4, 5]
   h = fnv1a_hash(seq5, 5)
   call check(h.eq.1109817072422714760_8, 'fnv1a_hash: [1,2,3,4,5] matches reference value', all_pass)

   ! A genuine 0/1 pattern, the actual shape fnv1a_hash is used for in
   ! practice (one pixel's own valid-channel mask row).
   seq7 = [integer(kind=1) :: 1, 0, 1, 1, 0, 0, 1]
   h = fnv1a_hash(seq7, 7)
   call check(h.eq.2228391952209758033_8, 'fnv1a_hash: [1,0,1,1,0,0,1] matches reference value', all_pass)

   ! Same bytes, called twice -- determinism (pure function, no hidden
   ! state), a basic sanity check worth having explicitly rather than
   ! assuming.
   call check(fnv1a_hash(seq7, 7).eq.fnv1a_hash(seq7, 7), 'fnv1a_hash: deterministic across repeated calls', all_pass)

   ! Different patterns must (for these specific short test vectors,
   ! verified directly, not just "hashes usually differ") produce
   ! different hashes -- catches a degenerate always-same-value bug.
   call check(fnv1a_hash(seq5, 5).ne.fnv1a_hash(seq7, 7), 'fnv1a_hash: distinct patterns hash differently', all_pass)

   ! linear_scan_extreme: min/max index, deterministic lowest-index
   ! tie-break, and the huge()-sentinel-always-wins-find_max property
   ! Belady's own admission control (T14 increment 8) will depend on.
   call test_linear_scan_extreme(all_pass)

   ! pattern_registry_t: known-repeat occurrence tracking (increment 4),
   ! growth past initial capacity (increment 4), and advance/next-
   ! occurrence timeline tracking (increment 5).
   call test_registry_known_repeats(all_pass)
   call test_registry_growth(all_pass)
   call test_registry_advance_and_next_occurrence(all_pass)
   call test_registry_lookup(all_pass)

   if (all_pass) then
      write(*,'(A)') '[PASS] test_rmclean_cache_mod: all checks passed'
   else
      write(*,'(A)') '[FAIL] test_rmclean_cache_mod: one or more checks failed'
      stop 1
   endif

contains

   subroutine check(cond, label, all_pass_io)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: label
      logical, intent(inout) :: all_pass_io
      if (cond) then
         write(*,'(A,A)') '  [PASS] ', label
      else
         write(*,'(A,A)') '  [FAIL] ', label
         all_pass_io = .false.
      endif
   end subroutine check

   subroutine test_linear_scan_extreme(all_pass_io)
      logical, intent(inout) :: all_pass_io
      integer(kind=8) :: vals(5)
      integer(kind=8) :: with_sentinel(4)

      vals = [integer(kind=8) :: 5, 1, 9, 1, 3]
      call check(linear_scan_extreme(vals, 5, .false.).eq.2,&
      &'linear_scan_extreme: min index (lowest-index tie-break, value 1 at idx 2 and 4)', all_pass_io)
      call check(linear_scan_extreme(vals, 5, .true.).eq.3,&
      &'linear_scan_extreme: max index (value 9 at idx 3)', all_pass_io)

      with_sentinel = [integer(kind=8) :: 100_8, huge(0_8), 50_8, 999999_8]
      call check(linear_scan_extreme(with_sentinel, 4, .true.).eq.2,&
      &'linear_scan_extreme: huge() sentinel always wins find_max=.true., even against large finite values', all_pass_io)
   end subroutine test_linear_scan_extreme

   subroutine test_registry_known_repeats(all_pass_io)
      !! Sequence A,B,A,C,A,B at scan positions 1..6 -- 3 distinct
      !! patterns, known occurrence lists.
      logical, intent(inout) :: all_pass_io
      type(pattern_registry_t) :: reg
      integer(kind=1) :: pat_a(3), pat_b(3), pat_c(3)
      integer :: id_a, id_b, id_c, id_tmp

      pat_a = [integer(kind=1) :: 1, 0, 1]
      pat_b = [integer(kind=1) :: 0, 1, 0]
      pat_c = [integer(kind=1) :: 1, 1, 0]

      call registry_init(reg)
      call registry_lookup_or_insert(reg, pat_a, 3, 1_8, id_a)
      call registry_lookup_or_insert(reg, pat_b, 3, 2_8, id_b)
      call registry_lookup_or_insert(reg, pat_a, 3, 3_8, id_tmp)
      call check(id_tmp.eq.id_a, 'pattern_registry: repeat lookup of A returns the SAME entry id', all_pass_io)
      call registry_lookup_or_insert(reg, pat_c, 3, 4_8, id_c)
      call registry_lookup_or_insert(reg, pat_a, 3, 5_8, id_tmp)
      call registry_lookup_or_insert(reg, pat_b, 3, 6_8, id_tmp)

      call check(reg%n_entries.eq.3, 'pattern_registry: 3 distinct patterns from a 6-pixel sequence with repeats', all_pass_io)
      call check(id_a.ne.id_b .and. id_a.ne.id_c .and. id_b.ne.id_c,&
      &'pattern_registry: 3 distinct patterns got 3 distinct entry ids', all_pass_io)

      call check(reg%entries(id_a)%n_occurrences.eq.3, 'pattern_registry: A occurs 3 times', all_pass_io)
      call check(all(reg%entries(id_a)%occurrences(1:3).eq.[1_8,3_8,5_8]),&
      &'pattern_registry: A occurrence list is [1,3,5]', all_pass_io)
      call check(reg%entries(id_b)%n_occurrences.eq.2, 'pattern_registry: B occurs 2 times', all_pass_io)
      call check(all(reg%entries(id_b)%occurrences(1:2).eq.[2_8,6_8]),&
      &'pattern_registry: B occurrence list is [2,6]', all_pass_io)
      call check(reg%entries(id_c)%n_occurrences.eq.1, 'pattern_registry: C occurs 1 time', all_pass_io)
      call check(reg%entries(id_c)%occurrences(1).eq.4_8,&
      &'pattern_registry: C occurrence list is [4]', all_pass_io)
   end subroutine test_registry_known_repeats

   subroutine test_registry_growth(all_pass_io)
      !! 20 distinct single-byte patterns (0..19) -- exceeds the
      !! registry's own initial capacity (16), forcing at least one
      !! entries(:)/buckets(:) growth+rebuild. Every pattern must
      !! remain correctly retrievable (stable entry id, correct
      !! occurrence list) afterward.
      logical, intent(inout) :: all_pass_io
      type(pattern_registry_t) :: reg
      integer(kind=1) :: pat(1)
      integer :: ids(0:19), id_tmp, k
      logical :: ids_stable, occurrences_ok

      call registry_init(reg)
      do k = 0, 19
         pat(1) = int(k, 1)
         call registry_lookup_or_insert(reg, pat, 1, int(k, 8)+1_8, ids(k))
      end do
      call check(reg%n_entries.eq.20,&
      &'pattern_registry: 20 distinct patterns registered (forces growth past initial capacity 16)', all_pass_io)

      ! Re-lookup each pattern with a NEW scan position -- entry id
      ! must be stable (same as the original insert) despite the
      ! entries(:) array having been reallocated by growth in between.
      ids_stable = .true.
      occurrences_ok = .true.
      do k = 0, 19
         pat(1) = int(k, 1)
         call registry_lookup_or_insert(reg, pat, 1, int(k, 8)+100_8, id_tmp)
         if (id_tmp.ne.ids(k)) ids_stable = .false.
         if (reg%entries(ids(k))%n_occurrences.ne.2) occurrences_ok = .false.
         if (.not. all(reg%entries(ids(k))%occurrences(1:2).eq.[int(k,8)+1_8, int(k,8)+100_8])) occurrences_ok = .false.
      end do
      call check(ids_stable, 'pattern_registry: entry ids survive growth (re-lookup returns the same id)', all_pass_io)
      call check(occurrences_ok, 'pattern_registry: occurrence lists correct for all 20 patterns after growth', all_pass_io)
   end subroutine test_registry_growth

   subroutine test_registry_advance_and_next_occurrence(all_pass_io)
      !! Same A,B,A,C,A,B sequence as test_registry_known_repeats --
      !! walks registry_advance in real scan order and checks
      !! registry_next_occurrence at each step, including convergence
      !! to the huge() sentinel once a pattern's last occurrence has
      !! passed.
      logical, intent(inout) :: all_pass_io
      type(pattern_registry_t) :: reg
      integer(kind=1) :: pat_a(3), pat_b(3), pat_c(3)
      integer :: id_a, id_b, id_c, id_tmp

      pat_a = [integer(kind=1) :: 1, 0, 1]
      pat_b = [integer(kind=1) :: 0, 1, 0]
      pat_c = [integer(kind=1) :: 1, 1, 0]

      call registry_init(reg)
      call registry_lookup_or_insert(reg, pat_a, 3, 1_8, id_a)
      call registry_lookup_or_insert(reg, pat_b, 3, 2_8, id_b)
      call registry_lookup_or_insert(reg, pat_a, 3, 3_8, id_tmp)
      call registry_lookup_or_insert(reg, pat_c, 3, 4_8, id_c)
      call registry_lookup_or_insert(reg, pat_a, 3, 5_8, id_tmp)
      call registry_lookup_or_insert(reg, pat_b, 3, 6_8, id_tmp)

      ! Before any advance: next occurrence is each pattern's own FIRST
      ! recorded position.
      call check(registry_next_occurrence(reg, id_a).eq.1_8,&
      &'pattern_registry: before any advance, A''s next occurrence is 1', all_pass_io)
      call check(registry_next_occurrence(reg, id_b).eq.2_8,&
      &'pattern_registry: before any advance, B''s next occurrence is 2', all_pass_io)
      call check(registry_next_occurrence(reg, id_c).eq.4_8,&
      &'pattern_registry: before any advance, C''s next occurrence is 4', all_pass_io)

      ! Walk the real scan order: position 1=A, 2=B, 3=A, 4=C, 5=A, 6=B,
      ! advancing whichever pattern is "encountered" at each position,
      ! and checking each pattern's own next-occurrence right after.
      call registry_advance(reg, id_a) ! consumed position 1
      call check(registry_next_occurrence(reg, id_a).eq.3_8,&
      &'pattern_registry: after 1 advance, A''s next occurrence is 3', all_pass_io)

      call registry_advance(reg, id_b) ! consumed position 2
      call check(registry_next_occurrence(reg, id_b).eq.6_8,&
      &'pattern_registry: after 1 advance, B''s next occurrence is 6', all_pass_io)

      call registry_advance(reg, id_a) ! consumed position 3
      call check(registry_next_occurrence(reg, id_a).eq.5_8,&
      &'pattern_registry: after 2 advances, A''s next occurrence is 5', all_pass_io)

      call registry_advance(reg, id_c) ! consumed position 4 -- C's only occurrence
      call check(registry_next_occurrence(reg, id_c).eq.huge(0_8),&
      &'pattern_registry: after C''s only occurrence is consumed, next occurrence is the sentinel (never again)', all_pass_io)

      call registry_advance(reg, id_a) ! consumed position 5 -- A's last occurrence
      call check(registry_next_occurrence(reg, id_a).eq.huge(0_8),&
      &'pattern_registry: after A''s 3rd (last) occurrence is consumed, next occurrence is the sentinel', all_pass_io)

      call registry_advance(reg, id_b) ! consumed position 6 -- B's last occurrence
      call check(registry_next_occurrence(reg, id_b).eq.huge(0_8),&
      &'pattern_registry: after B''s 2nd (last) occurrence is consumed, next occurrence is the sentinel', all_pass_io)
   end subroutine test_registry_advance_and_next_occurrence

   subroutine test_registry_lookup(all_pass_io)
      !! T14 increment 7: registry_lookup is the pure READ-ONLY
      !! counterpart to registry_lookup_or_insert -- Pass 1 (the real
      !! run) uses it to find an already-registered pattern's entry_id
      !! without appending a spurious duplicate occurrence. Must (a)
      !! find exactly the same entry_id registry_lookup_or_insert
      !! would, (b) NEVER mutate n_occurrences/occurrences (unlike
      !! registry_lookup_or_insert), and (c) report entry_id=0 for a
      !! pattern that was never registered.
      logical, intent(inout) :: all_pass_io
      type(pattern_registry_t) :: reg
      integer(kind=1) :: pat_a(3), pat_b(3), pat_unseen(3)
      integer :: id_a, id_b, id_tmp, id_unseen
      integer(kind=8) :: n_occ_a_before

      pat_a = [integer(kind=1) :: 1, 0, 1]
      pat_b = [integer(kind=1) :: 0, 1, 0]
      pat_unseen = [integer(kind=1) :: 1, 1, 1]

      call registry_init(reg)
      call registry_lookup_or_insert(reg, pat_a, 3, 1_8, id_a)
      call registry_lookup_or_insert(reg, pat_b, 3, 2_8, id_b)
      call registry_lookup_or_insert(reg, pat_a, 3, 3_8, id_tmp)

      call registry_lookup(reg, pat_a, 3, id_tmp)
      call check(id_tmp.eq.id_a, 'registry_lookup: finds the same entry_id as registry_lookup_or_insert for a known pattern', all_pass_io)
      call registry_lookup(reg, pat_b, 3, id_tmp)
      call check(id_tmp.eq.id_b, 'registry_lookup: finds the correct entry_id for a second known pattern', all_pass_io)

      call registry_lookup(reg, pat_unseen, 3, id_unseen)
      call check(id_unseen.eq.0, 'registry_lookup: entry_id=0 for a pattern that was never registered', all_pass_io)

      n_occ_a_before = reg%entries(id_a)%n_occurrences
      call registry_lookup(reg, pat_a, 3, id_tmp)
      call registry_lookup(reg, pat_a, 3, id_tmp)
      call check(reg%entries(id_a)%n_occurrences.eq.n_occ_a_before,&
      &'registry_lookup: repeated calls never mutate n_occurrences (read-only, unlike registry_lookup_or_insert)', all_pass_io)
   end subroutine test_registry_lookup

end program test_rmclean_cache_mod
