program test_mask_read_overflow
   !! planning/RMCLEAN_INTEGRATION_PLAN.md ticket T10a's own regression
   !! test: proves rmclean_io_mod's read_mask_cube correctly reads a
   !! mask cube whose TOTAL element count (nx*ny*nchan) exceeds 2^31 --
   !! the exact overflow that silently truncated the real Jennifer ASKAP
   !! mask cube (4501x4501x288 ~= 5.83 billion elements) at ~76 of 288
   !! channels for every pixel, before this ticket's fix (FTGPVB ->
   !! FTGPVBLL, 32-bit -> 64-bit firstelem/nelem).
   !!
   !! A genuine overflow test needs a fixture whose logical size exceeds
   !! 2^31 bytes (1 byte/element, no way around that for a byte-typed
   !! mask cube) -- unavoidably ~2.4GB here. Kept fast and lightweight
   !! via an OS-level SPARSE file: only 3 bytes are actually written
   !! (first element, a known-untouched middle element, and the very
   !! last element), leaving the rest as unwritten holes that read back
   !! as zero without consuming real disk blocks. The 2880-byte primary
   !! header is built BY HAND via raw stream I/O rather than through
   !! CFITSIO's own FTINIT/FTPHPR/FTCLOS -- confirmed directly that
   !! closing a freshly-FTPHPR'd file actually materializes (writes real
   !! disk blocks for) its entire declared data section, which would
   !! defeat the whole point of a lightweight fixture. CFITSIO is used
   !! only for the READ side (safe_ftopen + read_mask_cube's own
   !! FTGPVBLL), which is the actual thing under test and does not care
   !! how the file was created, only that it is standards-compliant.
   use rmclean_io_mod
   use fitsio_unit_mod
   implicit none

   integer, parameter :: nx = 2, ny = 2
   integer, parameter :: nchan = 600000000  ! nx*ny*nchan = 2.4e9 > 2^31
   character(len=*), parameter :: fname = 'tests/output/mask_overflow_fixture.fits'
   integer(kind=1), allocatable :: mask_out(:,:,:)
   integer :: status, stream_unit
   integer(kind=8), parameter :: header_bytes = 2880_8
   integer(kind=8) :: total_data_bytes
   logical :: ok
   integer(kind=1) :: marker
   character(len=80) :: card1, card2, card3, card4, card5, card6, card7
   character(len=2880) :: header_block
   character(len=20) :: valfield

   ok = .true.

   ! --- hand-build a minimal, standards-compliant primary header (7
   ! required cards + blank padding to exactly one 2880-byte block) ---
   write(valfield,'(A20)') 'T'
   card1 = 'SIMPLE  = '//valfield//' / conforms to FITS standard'
   write(valfield,'(I20)') 8
   card2 = 'BITPIX  = '//valfield//' / array data type'
   write(valfield,'(I20)') 3
   card3 = 'NAXIS   = '//valfield//' / number of array dimensions'
   write(valfield,'(I20)') nx
   card4 = 'NAXIS1  = '//valfield
   write(valfield,'(I20)') ny
   card5 = 'NAXIS2  = '//valfield
   write(valfield,'(I20)') nchan
   card6 = 'NAXIS3  = '//valfield
   card7 = 'END'
   header_block = card1//card2//card3//card4//card5//card6//card7
   ! (the remainder of header_block, columns 561:2880, is already
   ! blank-filled by the character(len=2880) assignment's own
   ! right-padding with spaces)

   open(newunit=stream_unit, file=fname, access='stream',&
   &form='unformatted', status='replace', action='write')
   write(stream_unit, pos=1) header_block
   close(stream_unit)

   total_data_bytes = int(nx,8)*int(ny,8)*int(nchan,8)

   ! --- sparse-write 3 marker bytes: first element, a middle element
   ! (left untouched -- must read back 0), and the very last element
   ! (only reachable if the full 2.4e9-element range is actually read) ---
   open(newunit=stream_unit, file=fname, access='stream',&
   &form='unformatted', status='old', action='readwrite')
   marker = 1_1
   write(stream_unit, pos=header_bytes+1) marker
   write(stream_unit, pos=header_bytes+total_data_bytes) marker
   close(stream_unit)

   ! --- the actual test: call the fixed read_mask_cube directly ---
   call read_mask_cube(fname, nx, ny, nchan, mask_out, status)

   if (status.ne.0) then
      print *, 'FAIL: read_mask_cube returned nonzero status:', status
      ok = .false.
   endif
   if (allocated(mask_out)) then
      if (mask_out(1,1,1).ne.1_1) then
         print *, 'FAIL: mask_out(1,1,1) =', mask_out(1,1,1),&
         &' expected 1 (first element)'
         ok = .false.
      endif
      if (mask_out(1,1,nchan/2).ne.0_1) then
         print *, 'FAIL: mask_out(1,1,nchan/2) =', mask_out(1,1,nchan/2),&
         &' expected 0 (untouched sparse hole)'
         ok = .false.
      endif
      if (mask_out(nx,ny,nchan).ne.1_1) then
         print *, 'FAIL: mask_out(nx,ny,nchan) =', mask_out(nx,ny,nchan),&
         &' expected 1 (the LAST of 2.4e9 elements -- reads 0 if the'//&
         &' 32-bit overflow bug is reintroduced, since a wrapped nelem'//&
         &' would never reach this far)'
         ok = .false.
      endif
   else
      print *, 'FAIL: mask_out never allocated'
      ok = .false.
   endif

   ! --- cleanup: this fixture is generated fresh every run, never
   ! committed (tests/output/ is gitignored), and should not linger --
   ! even on failure, since a stale multi-GB file would break the next
   ! run's own status='replace'/status='old' assumptions.
   open(newunit=stream_unit, file=fname, status='old')
   close(stream_unit, status='delete')

   if (ok) then
      print *, 'PASS: test_mask_read_overflow'
      stop 0
   else
      print *, 'FAIL: test_mask_read_overflow'
      stop 1
   endif
end program test_mask_read_overflow
