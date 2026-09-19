chelp+
! This program fits a set of sin/cos functions and polynomials
! of specified orders (order_h and order_p respectively) to a 
! specified section of the input data, and returns the best fit 
! array.
! Now nleft and lright refer to the bin numbers and not to no. 
! of bins on either side.
!
! The threshold (in units of rms noise) is supplied.
! When threshold > 0, we iteratively search for outliers and 
! exclude them from the data to be fitted.  -- 04 Sept, 2004
!      
! The fit_array on entry contains 1/0 flags indicating samples 
! to be excluded from the fit.  
!               
!
! author        Desh 
!
! date          06-Oct-1999
!
! reference     smooth baseline fitting, removal, sin/cos fn.s
!
chelp-



      subroutine poly_harm_fit(npoints,x_array,y_array,
     -                         nleft,nright,lr_exclude,
     -                         order_h,order_p,threshold,
     -                         silent,data_tag,
     -                         fit_array,fit_param)

      
      implicit none
      real*4    y_array(*),
     -          x_array(*),
     -          fit_array(*),
     -          fit_param(*),
     -          threshold
      integer*4 npoints,
     -          order_h,
     -          order_p,
     -          nleft,
     -          nright
      logical*4 lr_exclude,silent
      character*(*) data_tag


      integer*4     number, maxpts, maxpar
      parameter    (maxpts = 16416, maxpar = 128 )

      integer*4     i, j, k, m, nused, nused0, loglun/50/,
     -              certain_exclusions(maxpts), nchar, niter
      real*4        xcoord(maxpts), ycoord(maxpts),
     -              standard_deviation(maxpts), 
     -              param(maxpar),
     -              u(maxpts,maxpar),
     -              v(maxpar,maxpar), w(maxpar), 
     -              chipol, rms, x_offset, 
     -              minimum_chi, x_factor, harm_poly_function

      character*128  templine 
      external      harm_poly_curvefit


      if(.not.silent)then
c       open an old log file to append; if not found, open it afresh anyway
        i = -1
        open(unit=loglun,file='POLY_FIT.log',
     -       status='old',access='append',err=100)
        i = 1
100     if(i.lt.0)open(unit=loglun,file='POLY_FIT.log',status='unknown')
        write(loglun,*)' start of INPUT samples'
        write(loglun,*)(y_array(j),j=1,npoints)
        write(loglun,*)' end of INPUT samples'
      end if

c  Initialise some arrays/variables used in the fitting routines

      x_factor = (x_array(npoints) - x_array(1))
     -             *real(npoints)/real(npoints-1)/2.
      x_offset = (x_array(npoints) + x_array(1))/2.
     -             + x_factor/real(npoints)

      do i = 1,npoints
        standard_deviation(i) = 1.0
        if(fit_array(i).gt.0.0)then
          certain_exclusions(i) = 0
        else
          certain_exclusions(i) = 1
          standard_deviation(i) = 1000.0
        end if
      end do

        
      minimum_chi = 1.e22


      do j = 1,npoints
        ycoord(j) = y_array(j)
        xcoord(j) = (x_array(j)- x_offset)/x_factor
      end do

      number = 2*order_h + order_p + 1
      if (number .gt. maxpar) then
        if(.not.silent)then
          write (loglun,'(a,i2)') 'The maximum order allowed is ', 
     -                 maxpar
        end if
        order_h = (maxpar-order_p-1)/2
        number = maxpar
        if(.not.silent)then
          write (loglun,*)'Hence reducing the harmonic order to ',
     -                 order_h
        end if
      end if

      if(.not.silent)then
          write(loglun,'(a,a)')' Data_tag: ',data_tag(1:nchar(data_tag))
          write (loglun,*)'polynomial order : ',order_p,' (i.e. ',
     -                (order_p + 1),' terms)'
          write (loglun,*)'Harmonic order : ',order_h,' (i.e. ',
     -                2*order_h ,' terms)'

          templine(1:) = 'f(X) = c0'
          j = 10
          if(order_p.gt.0)then
            do i=1,order_p
               if(i.lt.10)then
                 write(templine(j:),'(a,i1,a,i1)')
     -                      ' + c',i,'*X**',i
                 j = j + 10       
               else
                 write(templine(j:),'(a,i2,a,i2)')
     -                      ' + c',i,'*X**',i
                 j = j + 12       
               end if
               if(j.gt.69)then
                 write(loglun,'(a)')templine(1:nchar(templine))
                 templine(1:) = ' '
                 j = 6
               end if
            end do
            if(nchar(templine).gt.6)then
              write(loglun,'(a)')templine(1:nchar(templine))
              templine(1:) = ' '
              j = 6
            end if
          end if
          if(order_h.gt.0)then
            do i=1,order_h
               k = order_p + 2*i - 1
               do m=1,2
                 if(k.lt.10)then
                   write(templine(j:),'(a,i1)')
     -                        ' + c',k
                   j = j + 5
                 else
                   write(templine(j:),'(a,i2)')
     -                        ' + c',k
                   j = j + 6
                 end if
                 if((i/2)*2.eq.i)then
                   if(m.eq.1)then
                     write(templine(j:),'(a,i2,a)')
     -                      '*cos(',i/2,'*pi*X)'
                   else
                     write(templine(j:),'(a,i2,a)')
     -                      '*sin(',i/2,'*pi*X)'
                   end if
                   j = j + 13       
                 else
                   if(m.eq.1)then
                     write(templine(j:),'(a,i2,a)')
     -                      '*cos(',i,'*pi*X/2)'
                   else
                     write(templine(j:),'(a,i2,a)')
     -                      '*sin(',i,'*pi*X/2)'
                   end if
                   j = j + 15       
                 end if
                 k = k + 1
               end do
               if(j.gt.61)then
                 write(loglun,'(a)')templine(1:nchar(templine))
                 templine(1:) = ' '
                 j = 6
               end if
            end do
            if(nchar(templine).gt.6)then
              write(loglun,'(a)')templine(1:nchar(templine))
              templine(1:) = ' '
              j = 6
            end if
          end if
          write(loglun,*)'where X = (freq_MHz - ',x_offset,')/',x_factor
      end if

      if(nleft.lt.nright)then
        if(lr_exclude) then
          do i = 1,nleft
            standard_deviation(i) = 1000.0
            certain_exclusions(i) = 1
          end do
          do i = nright,npoints
            standard_deviation(i) = 1000.0
            certain_exclusions(i) = 1
          end do
        else
          do i = nleft,nright
            standard_deviation(i) = 1000.0
            certain_exclusions(i) = 1
          end do
        end if
      end if
      if(threshold.gt.0.0)then ! we need to search for outliers and reduce their weights
        if(threshold.lt.3.0)threshold = 3.0
        if(.not.silent)then
          write (loglun,*)' '
          write (loglun,*)' We will iteratively exclude outliers '
          write (loglun,*)' that deviate by ',threshold,' sigmas'
        end if
      else
        threshold = 1000.0
      end if
      nused0 = 0
      do i=1,npoints
        if(certain_exclusions(i).eq.0)nused0 = nused0 + 1
      end do
      nused = nused0   ! to begin with
      if(.not.silent)then
        write (loglun,*)' Beginning with ',nused,' points to fit'
      end if

      niter = 0
      do while(.true.)

        niter = niter + 1
c       do the fit
        call svdfit_harm_poly
     -            (xcoord,ycoord,standard_deviation,npoints,
     -             param, number,order_h,order_p,
     -             u,v,w,maxpts,maxpar,chipol,
     -             harm_poly_curvefit)
        chipol = 0.0
        do i=1,npoints
          fit_array(i) = 0.0    ! just to be safe
          fit_array(i) = harm_poly_function( xcoord(i),
     -                         param, order_h, order_p)
          if(standard_deviation(i).le.2.)then
            rms = ycoord(i) - fit_array(i)
            chipol = chipol + rms*rms 
          end if
        end do
        rms = sqrt(chipol/real(nused))       !real(npoints))

        if(.not.silent)then
c  Give the user a numerical estimate of the fit
          write (loglun,*)(param(i),i=1,number)
          write (loglun,*)
          write (loglun,'(a,f11.2,a,f6.2,a,i6,a,i6,a)') 
     -        ' Chi sq. for harmonic fit = ',chipol,
     -        '; RMS  = ', rms,
     -        '; points used = ',nused,
     -        ' (out of ',npoints,')'
C ------------- following lines added by wasim --------------
          if(rms.ne.rms)then
                  write(*,*)'-------------------------------'
                  write(*,*)'          WARNING              '
                  write(*,*)'   RMS seems to be NAN...      '
                  write(*,*)'   One possibility may be      '
                  write(*,*)'You forgot to fill fit_array...'
                  write(*,*)'You may think about aborting...'
                  write(*,*)'-------------------------------'
          endif
C ------------- above lines added by wasim , 20 DEC 2008--------

          write (loglun,*)
        end if
c       load some useful info (to be returned)
        do i=1,number
          fit_param(i) = param(i)
        end do
        fit_param(number+1) = rms
        fit_param(number+2) = real(nused - number)    ! degrees of freedom
        fit_param(number+3) = x_offset
        fit_param(number+4) = x_factor
        

        j = nused0
        do i=1,npoints
          fit_array(i) = 0.0    ! just to be safe
          fit_array(i) = harm_poly_function( xcoord(i),
     -                         param, order_h, order_p)
          if(certain_exclusions(i).eq.0)then
            if(abs(ycoord(i)-fit_array(i)).gt.threshold*rms)then
              standard_deviation(i) = 1000.0
              j = j - 1
            else
              standard_deviation(i) = 1.0  
            end if
          else
            standard_deviation(i) = 1000.0
          end if
        end do
        if(j.lt.npoints/4)then ! we will be rejecting too many this way; reject none to reassess
          do i=1,npoints
            if(certain_exclusions(i).eq.0)then
              standard_deviation(i) = 1.0
            end if 
          end do
          nused = nused0
        else
          if(j.eq.nused)then
            if(.not.silent)then
              write(loglun,*)'==============================='
              write(loglun,*)'==============================='
              write(loglun,*)'   '
              close(unit=loglun)
            end if
            return
          else
            if(j.gt.nused.and.niter.gt.20)then
              if(.not.silent)then
                write(loglun,*)'==============================='
                write(loglun,*)'==============================='
                write(loglun,*)'   '
                close(unit=loglun)
              end if
              return
            else
              nused = j
            end if
          end if
        end if
        if(threshold.lt.0.0)then
            if(.not.silent)then
              write(loglun,*)'==============================='
              write(loglun,*)'==============================='
              write(loglun,*)'   '
              close(unit=loglun)
            end if
            return
        end if
      end do

      return
      end


c-------------------------------------------------------------------------------
c                          S U B R O U T I N E S
c-------------------------------------------------------------------------------
c
c
c name          harm_poly_function
c
c function      This function calculates the value of a harmonic of
c               order n, with coefficients param(i) at the point x
c
c call
      real*4 function harm_poly_function(x, param, n_h, n_p)
      implicit none
      integer*4         maxpar
      parameter       ( maxpar = 128 )

      real*4            x,              ! The point to evaluate the harmonic 
     -                                  ! at
     -                  param(maxpar)   ! The coefficients of the harmonic
      integer*4         n_h, n_p        ! The order of the harmonic and
                                        ! polynomial function respectively
                                        ! (so that n_h*2 + n_p + 1 = the
                                        ! number of elements in param)
c
c author        desh 
c 
c date          31-Aug-1995
c
c refe          curve fitting, sin/cos fn.s, polynomials
c- 

      integer*4         i, iparam
      real*4            pi_by_2
  
      pi_by_2  = acos(-1.)/2.
      

      harm_poly_function = param(1)        ! dc; zeroth order term
      iparam = 1
      
      if(n_p.gt.0)then
        do i = 1,n_p
          iparam = iparam + 1
          harm_poly_function = harm_poly_function 
     -       + param(iparam)* x**real(i)
        end do
      end if

      if(n_h.gt.0)then
        do i = 1,n_h    
c         first cos term
          iparam = iparam + 1
          harm_poly_function = harm_poly_function 
     -       + param(iparam)* cos(real(i)*pi_by_2*x)
c         and then the sin term
          iparam = iparam + 1
          harm_poly_function = harm_poly_function 
     -       + param(iparam)* sin(real(i)*pi_by_2*x)
        end do
      end if

      return
      end

c+
c
c name          harm_poly_curvefit
c
c function      This subroutine provides the basis functions 
c               (that are sin/cos and/or polynomials) for the least
c               squares harmonic fitting routine, to be used with the
c               routine SVDFIT_harm_poly.
c
c call
      subroutine harm_poly_curvefit
     -           (x,basis_functions,order_h,order_p)
      implicit none
      integer*4         maxpar
      parameter       ( maxpar = 128 )

      integer*4         order_h         ! The order of the harmonic function
                                        ! to be fitted 
      integer*4         order_p         ! The order of the polynomial function
                                        ! to be fitted
      real*4            x,              ! x-coordinate
     -                  basis_functions(maxpar)
                                        ! The basis functions 
                                        ! for the fit
c
c author       desh 
c
c date          31-Aug-1995
c
c refe          general, curve/baseline fitting
c-

      integer*4         i, iparam
      real*4            pi_by_2


      pi_by_2  = acos(-1.)/2.



      basis_functions(1) = 1.0
      iparam = 1 
      
      if(order_p.gt.0)then
        do i = 1,order_p
          iparam = iparam + 1
          basis_functions(iparam) = x**real(i)
        end do
      end if

      if(order_h.gt.0)then
        do i = 1,order_h
          iparam = iparam + 1
          basis_functions(iparam) = cos(real(i)*pi_by_2*x)
          iparam = iparam + 1
          basis_functions(iparam) = sin(real(i)*pi_by_2*x)
        end do
      end if
      

      return
      end
c=======================================================================



C********************************************************************
c                               SUBROUTINES
C********************************************************************

      SUBROUTINE SVDFIT_harm_poly(X,Y,SIG,NDATA,A,MA,
     -             ma_h,ma_p,U,V,W,MP,NP,CHISQ,FUNCS)

      INTEGER*4   NDATA,MA,NP,MP,NMAX,MMAX,ma_h,ma_p
      REAL*4      TMP,WMAX,TOL,THRESH,CHISQ,SUM
      PARAMETER   (NMAX=16416,MMAX=128,TOL=1.E-4)
      REAL*4      X(*),Y(*),SIG(*),A(MA),V(NP,NP),
     *            U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)


      DO 12 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA_H,MA_P)
c       write(*,*)'afunc= ',afunc(1),afunc(2)
        IF (SIG(I).LE.0.0d0) SIG(I) = 1.0d0
        TMP=1.d0/SIG(I)
        DO 11 J=1,MA
          U(I,J)=AFUNC(J)*TMP
11      CONTINUE
        B(I)=Y(I)*TMP
12    CONTINUE

      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.d0
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.0d0
14    CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.d0

      DO 16 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA_H,MA_P)
        SUM=0.d0
        DO 15 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
15      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
16    CONTINUE

      RETURN
      END

c====================================================================

      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)

      PARAMETER   (NMAX=128)
      INTEGER*4   M,N,MP,NP,L,I
      REAL*4      A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      REAL*4      ANORM,C,F,G,H,S,SCALE,X,Y,Z
      !REAL*4     PYTHAG
      G=0.d0
      SCALE=0.d0
      ANORM=0.d0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.d0
        S=0.d0
        SCALE=0.d0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.d0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.d0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                IF (H.NE.0.0d0) THEN
                  F=S/H
                ENDIF
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            if(H.NE.0.0d0) THEN
              DO 19 K=L,N
                RV1(K)=A(I,K)/H
19            CONTINUE
            ENDIF
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.d0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.d0) THEN
            IF (A(I,L).NE.0.d0) THEN
              DO 26 J=L,N
                V(J,I)=(A(I,J)/A(I,L))/G
26            CONTINUE
            ENDIF
            DO 29 J=L,N
              S=0.d0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.d0
            V(J,I)=0.d0
31        CONTINUE
        ENDIF
        V(I,I)=1.d0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.d0
33        CONTINUE
        ENDIF
        IF (G.NE.0.d0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.d0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              IF (A(I,I).NE.0.d0) THEN
                F=(S/A(I,I))*G
                DO 35 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
35              CONTINUE
              ENDIF
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.d0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.d0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.d0
          S=1.d0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
c	      H=PYTHAG(F,G)
              W(I)=H
              IF (H.NE.0.d0) THEN
                H=1.d0/H
              ENDIF
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.d0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          !IF (ITS.EQ.100) PAUSE 'No convergence in 100 iterations'
          ! Modified by wasim to replace obsolete PAUSE in f95 
          IF (ITS.EQ.100) THEN
                  write(*,*)'No convergence in 100 iterations'
                  write(*,*)"Press ENTER to continue"
                  read(*,*)
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          IF (H*Y.NE.0.0d0) THEN
            F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          ENDIF
          G=SQRT(F*F+1.d0)
c	  G=PYTHAG(F,1.d0)
          IF (X.NE.0.0d0) THEN
            F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          ENDIF
          C=1.d0
          S=1.d0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
c	    Z=PYTHAG(F,H)
            RV1(J)=Z
            IF (Z.NE.0.0d0) THEN
              C=F/Z
              S=H/Z
            ENDIF
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
c	    Z=PYTHAG(F,H)
            W(J)=Z
            IF (Z.NE.0.d0) THEN
              Z=1.d0/Z
              C=F*Z
              S=H*Z
            ENDIF 
            F= (C*G)+(S*Y) 
            X=-(S*G)+(C*Y) 
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.d0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END

c====================================================================

      FUNCTION PYTHAG(A,B)

      REAL*4   A,B,PYTHAG
      REAL*4   ABSA,ABSB

      ABSA=ABS(A)
      ABSB=ABS(B)
      IF(ABSA.GT.ABSB)THEN
         PYTHAG=ABSA*SQRT(1.d0+(ABSB/ABSA)**2)
      ELSE
         IF(ABSB.EQ.0.)THEN
            PYTHAG=0.D0
         ELSE
            PYTHAG=ABSB*SQRT(1.d0+(ABSA/ABSB)**2)
         ENDIF
      ENDIF

      RETURN
      END

C=====================================================================

      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)

      PARAMETER (NMAX=128)
      REAL*4    U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP)
      REAL*4    TMP(NMAX),S

      DO 12 J=1,N
        S=0.d0
        IF(W(J).NE.0.d0)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.d0
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END

