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
