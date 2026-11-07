! Minimal Fortran interface to the cuBLAS routines used by the CUDA backend (C-binding, no CUDA Fortran needed).
module cublas_iface
  use iso_c_binding
  implicit none
  integer(c_int), parameter :: CUBLAS_OP_N = 0, CUBLAS_OP_T = 1
  integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS = 0
  interface
    integer(c_int) function cublasCreate_v2(handle) bind(C, name="cublasCreate_v2")
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: handle
    end function
    integer(c_int) function cublasDestroy_v2(handle) bind(C, name="cublasDestroy_v2")
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
    end function
    ! C = alpha * op(A) * op(B) + beta * C, column-major; alpha and beta are host scalars (default pointer mode)
    integer(c_int) function cublasSgemm_v2(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
        bind(C, name="cublasSgemm_v2")
      import :: c_int, c_ptr, c_float
      type(c_ptr), value :: handle
      integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
      real(c_float), intent(in) :: alpha, beta
      type(c_ptr), value :: a, b, c
    end function
  end interface
end module
