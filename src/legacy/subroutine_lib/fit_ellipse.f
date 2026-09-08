chelp+
      ! ------------------------------------------------
      ! This subroutine determines the best-fit ellipse 
      ! to a given set of X-Y scatter points. 
      !
      !  Inputs: 
      !          X -> The X values
      !          Y -> The Y values 
      !          N -> Number of points 
      !
      !  Output: 
      !          a    -> semi-axis(along x) 
      !          b    -> semi-axis(along y) 
      !          phi  -> tilt angle (radians) 
      !          X0   -> X-center of non-tilted ellipse
      !          Y0   -> Y-center of non-tilted ellipse
      !          X0_t -> X-center of tilted ellipse
      !          Y0_t -> Y-center of tilted ellipse
      !          amaj -> major axis length of the ellipse
      !          amin -> minor axis length of the ellipse  
      !          
      !                               --wr, 31 Dec, 2011
      !  Reference: fit_ellipse.m (from www, could not get 
      !                            the author name though)
      ! ------------------------------------------------
chelp-
* --------------------------------------------------------------------------------------------
* fit_ellipse - finds the best fit to an ellipse for the given set of points.
*
* Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
*
* Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
*           axis_handle - optional. a handle to an axis, at which the estimated ellipse 
*                         will be drawn along with it's axes
*
* Output:   ellipse_t - structure that defines the best fit to an ellipse
*                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
*                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
*                       phi         - orientation in radians of the ellipse (tilt)
*                       X0          - center at the X axis of the non-tilt ellipse
*                       Y0          - center at the Y axis of the non-tilt ellipse
*                       X0_in       - center at the X axis of the tilted ellipse
*                       Y0_in       - center at the Y axis of the tilted ellipse
*                       long_axis   - size of the long axis of the ellipse
*                       short_axis  - size of the short axis of the ellipse
*                       status      - status of detection of an ellipse
*
* Note:     if an ellipse was not detected (but a parabola or hyperbola), then
*           an empty structure is returned

* =====================================================================================
*                  Ellipse Fit using Least Squares criterion
* =====================================================================================
* We will try to fit the best ellipse to the given measurements. the mathematical
* representation of use will be the CONIC Equation of the Ellipse which is:
* 
*    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
*   
* The fit-estimation method of use is the Least Squares method (without any weights)
* The estimator is extracted from the following equations:
*
*    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
*
*    where:
*       A   - is the vector of parameters to be estimated (a,b,c,d,e)
*       x,y - is a single measurement
*
* We will define the cost function to be:
*
*   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
*            = (X*A+f_c)'*(X*A+f_c) 
*            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
*
*   where:
*       g_c(x_c,y_c;A) - vector function of ALL the measurements
*                        each element of g_c() is g(x,y;A)
*       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
*       f_c            - is actually defined as ones(length(f),1)*f
*
* Derivation of the Cost function with respect to the vector of parameters "A" yields:
*
*   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
*
* Which yields the estimator:
*
*       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
*       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*
* (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
*  
* NOW, all that is left to do is to extract the parameters from the Conic Equation.
* We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
*
*    Recall the conic representation of an ellipse:
* 
*       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
* 
* We will check if the ellipse has a tilt (=orientation). The orientation is present
* if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
* tilt of the ellipse.
*
* If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
* we will remove the tilt of the ellipse so as to remain with a conic representation of an 
* ellipse without a tilt, for which the math is more simple:
*
* Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
*
* We will remove the orientation using the following substitution:
*   
*   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
*   
*   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
*
*   where:      c = cos(phi)    ,   s = sin(phi)
*
*   and simplify...
*
*       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
*           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
*
*   The orientation is easily found by the condition of (B_new=0) which results in:
* 
*   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
*   
*   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
*   all the other constants A`,C`,D`,E` can be found.
*
*   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
*   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c 
*   C` = A*s^2 + B*c*s + C*c^2
*
* Next, we want the representation of the non-tilted ellipse to be as:
*
*       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
*
*       where:  (X0,Y0) is the center of the ellipse
*               a,b     are the ellipse "radiuses" (or sub-axis)
*
* Using a square completion method we will define:
*       
*       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
*
*       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
*                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
*
*       which yields the transformations:
*       
*           X0  =   -D`/(2*A`)
*           Y0  =   -E`/(2*C`)
*           a   =   sqrt( abs( F``/A` ) )
*           b   =   sqrt( abs( F``/C` ) )
*
* And finally we can define the remaining parameters:
*
*   long_axis   = 2 * max( a,b )
*   short_axis  = 2 * min( a,b )
*   Orientation = phi
*


      subroutine fit_ellipse(X,Y,N,a,b,phi,X0,Y0,X0_t,Y0_t,
     -                       amaj,amin)

      implicit none

      real*4          X(*), Y(*)
      real*4          a, b, phi, X0, Y0, X0_t, Y0_t, 
     -                amaj, amin 
      integer*4       N, i  

      ! Define local variables here: 
      real*4          xarr(N), yarr(N) 
      real*4          mean_x, mean_y, ax, ay 
      real*4          F(N,5),FT(5,N),P(5,5),PInv(5,5),S(1,5),AA(1,5)
      real*4          c, d, e
      real*4          pi, orientation_tol 
      real*4          cos_phi, sin_phi 
      real*4          err 

      pi = acos(-1.0)

      orientation_tol = 1.0e-3 


      call mean(X,N,mean_x)
      call mean(Y,N,mean_y)
!      write(*,*)"mu_x, mu_y: ",mean_x, mean_y 
      ! Store the mean removed data for better 
      ! accuracy in the inversion process. We'll 
      ! add the mean later: 
      do i = 1,N
         xarr(i) = X(i) - mean_x 
         yarr(i) = Y(i) - mean_y 
      enddo

      S(1,1) = 0.0
      S(1,2) = 0.0
      S(1,3) = 0.0
      S(1,4) = 0.0
      S(1,5) = 0.0

      do i = 1,N
         F(i,1) = xarr(i)*xarr(i) 
         F(i,2) = xarr(i)*yarr(i) 
         F(i,3) = yarr(i)*yarr(i) 
         F(i,4) = xarr(i) 
         F(i,5) = yarr(i) 

         S(1,1) = S(1,1) + F(i,1)
         S(1,2) = S(1,2) + F(i,2)
         S(1,3) = S(1,3) + F(i,3)
         S(1,4) = S(1,4) + F(i,4)
         S(1,5) = S(1,5) + F(i,5)
      enddo
      call transpose_2d(F,N,5,N,5,FT)
      call matmult_2d(FT,5,N,5,N,
     -                 F,N,5,N,5,P)

      call inverse(P,5,5,PInv)
      call matmult_2d(S,1,5,1,5,PInv,5,5,5,5,AA)

!      write(*,*)(S(1,i),i=1,5)
!      write(*,*)" "
!      write(*,*)(AA(1,i),i=1,5)

      ! The parameters of the conic section for the ellipse: 
      a = AA(1,1) 
      b = AA(1,2) 
      c = AA(1,3) 
      d = AA(1,4) 
      e = AA(1,5) 

      if(min(abs(b/a),abs(b/c)) .gt. orientation_tol)then
              phi = 0.5*atan(b/(c-a))
              cos_phi = cos(phi)
              sin_phi = sin(phi)

              a = AA(1,1)*cos_phi**2 - AA(1,2)*cos_phi*sin_phi + 
     -             AA(1,3)*sin_phi**2
              b = 0.0 
              c = AA(1,1)*sin_phi**2 + AA(1,2)*cos_phi*sin_phi + 
     -             AA(1,3)*cos_phi**2
              d = AA(1,4)*cos_phi - AA(1,5)*sin_phi 
              e = AA(1,4)*sin_phi + AA(1,5)*cos_phi 

              ! Correct the means as well 
              ax = mean_x  ! Only temporarily
              ay = mean_y  ! Only temporarily 

              mean_x = cos_phi*ax - sin_phi*ay
              mean_y = sin_phi*ax + cos_phi*ay
      else
              phi = 0.0 
              cos_phi = cos(phi)
              sin_phi = sin(phi)
      endif

      ! Check if the equation represents an ellipse: 
      err = a*c 
      if(err .eq. 0.0)then
              write(*,*)"WARNING!!!"
              write(*,*)"fit_ellipse could not locate an ellipse!! "
              write(*,*)"Instead, PARABOLA found!!"
      else if(err .lt. 0.0)then
              write(*,*)"WARNING!!!"
              write(*,*)"fit_ellipse could not locate an ellipse!! "
              write(*,*)"Instead, HYPERBOLA found!!"
      !else
              !write(*,*)"Successfully fitted an ellipse!"
      endif

      if(err .gt. 0 )then 
              ! Take care of negative coefficients: 
              if(a .lt. 0.0)then
                      a = -a
                      b = -b 
                      c = -c 
                      d = -d 
                      e = -e
              endif
              X0 = mean_x - 0.5*d/a 
              Y0 = mean_y - 0.5*e/c 

              ax = 1.0 + (d**2)/(4.0*a) + (e**2)/(4.0*c) 
              a = sqrt(ax/a)
              b = sqrt(ax/c)

              amaj = 2.0*max(a,b)
              amin = 2.0*min(a,b)

              ! Rotate the ellipse to its original tilt: 
              X0_t = cos_phi*X0 + sin_phi*Y0 
              Y0_t = -sin_phi*X0 + cos_phi*Y0 

!              write(*,*)"The parameters of the ellipse are: "
!              write(*,*)"a: ",a
!              write(*,*)"b: ",b
!              write(*,*)"tilt(deg): ",phi*180.0/pi 
!              write(*,*)"X0: ",X0
!              write(*,*)"Y0: ",Y0
!              write(*,*)"X0_t: ",X0_t
!              write(*,*)"Y0_t: ",Y0_t
!              write(*,*)"amaj: ",amaj 
!              write(*,*)"amin: ",amin 
      else
              write(*,*)"No ellipse could be fitted..."
              a = 0.0
              b = 0.0 
              phi = 0.0
              X0 = 0.0
              Y0 = 0.0 
              X0_t = 0.0 
              Y0_t = 0.0 
              amaj = 0.0
              amin = 0.0 
      endif




      return 
      end



