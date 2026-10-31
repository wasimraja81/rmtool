chelp+
! This subroutine calculates the inverse of any
! 2x2 matrix. 
! 
! It was designed to calculate the inverse of the 
! broken up components of the Mueller Matrix. 
! 
! This subroutine is specific to 2x2 matrices only.
! 
!          --wasim raja, rri, 19 Apr, 2010.
! 
chelp-

      subroutine inv_2x2(A,I)

      implicit none

      real*4 A(2,2), I(2,2)

      real*4 a11, a12, 
     -       a21, a22, 
     -       i11, i12, 
     -       i21, i22

      real*4 det


      a11 = A(1,1)
      a12 = A(1,2)

      a21 = A(2,1)
      a22 = A(2,2)

      det = 1.0/(a11*a22 - a12*a21)

      i11 = det*a22
      i12 = -det*a12

      i21 = -det*a21
      i22 = det*a11


      I(1,1) = i11
      I(1,2) = i12

      I(2,1) = i21
      I(2,2) = i22


      return
      end
