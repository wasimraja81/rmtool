chelp+
C Include file for fft2d.f routine.
C------------------------------------
C You are expected to include this 
C file at the beginning of your main
C program. Care has been taken to 
C prevent clash with other include 
C files by appending the maxdim
C parameter with an "_fft2d"
C
C You are also expected to include
C the following routines in your main
C program:
C
C   1) 'fft_general_lin.f'  
C   2) 'fort_lib.f'  
C   3) 'shift.f' 
C   4) 'fftshift.f' 
C
C NOTE: Currently, the maximum array 
C       that can be handled is of 
C       size 16384 X 8192.
C
C       In case you need to modify
C       this default size, change 
C       the maxdim_fft2d parameter
C       accordingly, and recompile 
C       your codes.
C
C    --wasim raja, 30 Dec, 2009
C------------------------------------
chelp-

       integer*4 maxdimx_fft2d,maxdimy_fft2d
       parameter(maxdimx_fft2d = 16384)
       parameter(maxdimy_fft2d = 8192)

