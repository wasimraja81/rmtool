c_test      real*8          this_jdtime,
c_test     -                this_ra,this_dec,
c_test     -                new_jdtime,
c_test     -                new_ra,new_dec
c_test      real*4          pi
c_test
c_test      write(*,*)' RA-DEC (hr,deg)?'
c_test      read(*,*)this_ra, this_dec
c_test      this_ra = this_ra*2.*pi/24.
c_test      this_dec = this_dec*2.*pi/360.
c_test      write(*,*)'this_epoch, new_epoch ?'
c_test      read(*,*)this_jdtime,new_jdtime
c_test      call precess(this_jdtime,new_jdtime,
c_test     -            this_ra,this_dec,new_ra,new_dec) 
c_test      stop
c_test      end


      subroutine precess(this_jdtime_in,new_jdtime_in, 
     -            this_ra,this_dec,new_ra,new_dec)
c     ra,dec is expected in radians, and use years for B1950 or J2000

      implicit none
      real*8          this_jdtime_in,this_jdtime,temp_jdtime,
     -                this_ra,this_dec,
     -                new_jdtime_in,new_jdtime,
     -                new_ra,new_dec,
     -                x
      character*20    rastrg,
     -                decstrg
      logical*4       convert_fk4fk5
cc      character*5     bjepoch
c      
c----------------------------------------------------------------------------
      integer*4          nchar, 
     -                   i, j, k, hh, mm, dd        !,as
      real*4             ss
      real*8             ra,dec, 
     -                   pi,twopi
      double precision   pmat(3,3),v1(3),v2(3),v3(3)
cc      character*80       tempstring
      character*1        isgn
cc      character*1        colon,dummy1,dummy2
c----------------------------------------------------------------------------

c_bypass_test    
c_bypass_test               new_ra = this_ra 
c_bypass_test               new_dec = this_dec 
c_bypass_test               return
c_bypass_test    


      this_jdtime = this_jdtime_in
      new_jdtime = new_jdtime_in

      pi = acos(-1.d0)
      twopi = 2.d0*pi

      if(this_jdtime.ne.1950.and.this_jdtime.ne.2000.)then
       if(this_jdtime.lt.100000.)then
         this_jdtime = this_jdtime + 2400000.d0
       end if
      end if
      if(new_jdtime.ne.1950.and.new_jdtime.ne.2000.)then
       if(new_jdtime.lt.100000.)then
         new_jdtime = new_jdtime + 2400000.d0
       end if
      end if
      convert_fk4fk5 = .false.
      if(this_jdtime.eq.1950.d0)then
        this_jdtime = 2433282.423d0
        convert_fk4fk5 = .true.
      else if(this_jdtime.eq.2000.d0)then
        this_jdtime = 2451545.0d0
      end if
      

      if(new_jdtime.eq.1950.d0)then
        new_jdtime = 2433282.423d0
      else if(new_jdtime.eq.2000.d0)then
        new_jdtime = 2451545.0d0
      end if
        

      ra = this_ra
      ss = this_ra*3600.*24./(2.*pi)
      hh = ss/3600.
      ss = ss - real(hh)*3600.
      mm = ss/60.
      ss = ss -real(mm)*60.
      write(rastrg,"(i2,':',i2,':',f5.2)")hh,mm,ss

      do i=1,11   !nchar(rastrg)
        if(rastrg(i:i).eq.' ') rastrg(i:i) = '0'
      end do

      isgn = '+'
      if(this_dec.lt.0.0)isgn = '-'

      dec = this_dec
      ss = abs(this_dec*360./(2.*pi)*3600.)
      
      dd = ss/3600.
      ss = ss - real(dd)*3600.
      mm = ss/60.
      ss = ss -real(mm)*60.
      write(decstrg,"(a,i2,':',i2,':',f5.2)")isgn,dd,mm,ss
      do i=1,12   !nchar(decstrg)
        if(decstrg(i:i).eq.' ') decstrg(i:i) = '0'
      end do

      write(6,'(a,4x,a,i4)') 
     -   rastrg(1:nchar(rastrg)),decstrg(1:nchar(decstrg)),
     -   convert_fk4fk5
c_test      write(6,5) rastrg(1:11),decstrg(1:12)
c_test5     format(
c_test     -       '               RA : ',a20,'   (1950.)'/,
c_test     -       '              Dec : ',a20,'   (1950.)' )
c_test      write(6,*)'  *+*+*+*+*+*+*+*+*+*+*+*+*+'
c
      v1(1) = dcos(ra)*dcos(dec)
      v1(2) = dsin(ra)*dcos(dec)
      v1(3) = dsin(dec)
      if(convert_fk4fk5)then
        call fk4fk5(v1,v2)
        temp_jdtime = 2451545.0d0
        if(temp_jdtime.eq.new_jdtime)then
          new_dec = dasin(v2(3))
          new_ra  = datan2(v2(2),v2(1))
          if (new_ra.lt.0.0d0) new_ra = new_ra + 2.d0*pi
          go to 2000
        end if
      else
        v2(1) = v1(1)
        v2(2) = v1(2)
        v2(3) = v1(3)
      end if
      call precsn(this_jdtime,new_jdtime,pmat)
      do 2 j = 1, 3
        x = 0.0d0
        do 1 k = 1, 3
          x = x+pmat(j, k)*v2(k)
    1   continue
        v3(j) = x
    2 continue
      dec = dasin(v3(3))
      ra  = datan2(v3(2),v3(1))
      if (ra.lt.0.0d0) ra = ra + 2.d0*pi

      new_ra = ra
      new_dec = dec
c
c_test      write(*,*)'  '
c_test      write(*,*)'Current Coordinates:'
c_test      call getradec(v3)
c_test      write(*,*)'   '
c
 
2000  ra = new_ra
      ss = new_ra*3600.*24./(2.*pi)
      hh = ss/3600.
      ss = ss - real(hh)*3600.
      mm = ss/60.
      ss = ss -real(mm)*60.
      write(rastrg,"(i2,':',i2,':',f5.2)")hh,mm,ss

      do i=1,11   !nchar(rastrg)
        if(rastrg(i:i).eq.' ') rastrg(i:i) = '0'
      end do

      isgn = '+'
      if(this_dec.lt.0.0)isgn = '-'

      dec = new_dec
      ss = abs(new_dec*360./(2.*pi)*3600.)
      
      dd = ss/3600.
      ss = ss - real(dd)*3600.
      mm = ss/60.
      ss = ss -real(mm)*60.
      write(decstrg,"(a,i2,':',i2,':',f5.2)")isgn,dd,mm,ss
      do i=1,12   !nchar(decstrg)
        if(decstrg(i:i).eq.' ') decstrg(i:i) = '0'
      end do

      write(6,'(a,4x,a,a,f15.3,2x,f15.3)') 
     -   rastrg(1:nchar(rastrg)),decstrg(1:nchar(decstrg)),
     -   ' after precess:',this_jdtime,new_jdtime

      return
      end      



cc      OPTIONS/NOF77

      subroutine fk4fk5(Vin,Vout)
      implicit none
      double precision Vin(3), Vout(3)
c+
c_name        fk4fk5
c
c_function   To convert a FK4 based position referred to the
c            mean equinox and equator of B1950.0 to positions
c            in the FK5 system referred to the equinox and
c            equator of J2000.0 in accordance with the
c            new (1976) IAU resolutions. Reference :
c            Aoki,S., et al, 1983. Astron.Astrophys.,128,p263.
c
c_call       call fk4fk5(vin,vout)
c_/vin       r*8(3)  FK4 position as 3-vector
c_/vout      r*8(3)  Returned FK5 position as 3-vector
c
c_author      D McConnell
c_created     11 September 1984
c_latest      23-OCT-1986
c
c_reference  Mt Pleasant, library
c_keyword    Ephemerides, coordinates
c-

      double precision       P11,P12,P13,P21,P22,P23,P31,P32,P33,
     -                       A1,A2,A3,w

      data                   P11 /  +0.99992 56782 D0/,
     -                       P12 /  -0.01118 20610 D0/,
     -                       P13 /  -0.00485 79477 D0/,

     -                       P21 /  +0.01118 20609 D0/,
     -                       P22 /  +0.99993 74784 D0/,
     -                       P23 /  -0.00002 71765 D0/,

     -                       P31 /  +0.00485 79479 D0/,
     -                       P32 /  -0.00002 71474 D0/,
     -                       P33 /  +0.99998 81997 D0/

      data                   A1 / -1.62557 D-06/,
     -                       A2 / -0.31919 D-06/,
     -                       A3 / -0.13843 D-06/

c  Allow for e-terms
      w  = VIN(1) * a1 + VIN(2) * a2 + VIN(3) * a3

      VIN(1) = VIN(1) - A1 + w * VIN(1)
      VIN(2) = VIN(2) - A2 + w * VIN(2)
      VIN(3) = VIN(3) - A3 + w * VIN(3)

c  Precess
      Vout(1) = P11 * VIN(1) + P12 * VIN(2) + P13 * VIN(3)
      Vout(2) = P21 * VIN(1) + P22 * VIN(2) + P23 * VIN(3)
      Vout(3) = P31 * VIN(1) + P32 * VIN(2) + P33 * VIN(3)

      return
      end



C
C
	SUBROUTINE PRE(DEQ1,DEQ2,DPREMA)
C
C**********************************************************************
C  CALCULATION OF THE MATRIX OF GENERAL PRECESSION FROM DEQ1 TO DEQ2.
C  THE PRECESSION ANGLES (DZETA,DZETT,DTHET) ARE COMPUTED FROM THE
C  CONSTANTS (DC1-DC9) CORRESPONDING TO THE DEFINITIONS IN THE
C  EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL EPHEREMIS (1961,P.30F).
C**********************************************************************
C
	IMPLICIT REAL*8 (D)
	DIMENSION DPREMA(3,3)
	DATA DCSAR/4.848136812D-6/,DC1900/1900.0D0/,DC1M2/0.01D0/,
     * DC1/2304.25D0/,DC2/1.396D0/,DC3/0.302D0/,DC4/0.018D0/,
     * DC5/0.791D0/,DC6/2004.683D0/,DC7/-0.853D0/,DC8/-0.426D0/,
     * DC9/-0.042D0/
	DT0=(DEQ1-DC1900)*DC1M2
	DT=(DEQ2-DEQ1)*DC1M2
	DTS=DT*DT
	DTC=DTS*DT
	DZETA=((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT=DZETA + DC5*DTS*DCSAR
	DTHET=((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSZETA=DSIN(DZETA)
	DCZETA=DCOS(DZETA)
	DSZETT=DSIN(DZETT)
	DCZETT=DCOS(DZETT)
	DSTHET=DSIN(DTHET)
	DCTHET=DCOS(DTHET)
	DA=DSZETA*DSZETT
	DB=DCZETA*DSZETT
	DC=DSZETA*DCZETT
	DD=DCZETA*DCZETT
	DPREMA(1,1)= DD*DCTHET - DA
	DPREMA(1,2)=-DC*DCTHET - DB
	DPREMA(1,3)=-DSTHET*DCZETT
	DPREMA(2,1)= DB*DCTHET + DC
	DPREMA(2,2)=-DA*DCTHET + DD
	DPREMA(2,3)=-DSTHET*DSZETT
	DPREMA(3,1)= DCZETA*DSTHET
	DPREMA(3,2)=-DSZETA*DSTHET
	DPREMA(3,3)= DCTHET
	RETURN
	END


      double precision function gmst(julday)
      implicit none
      double precision julday
c+
c_name   gmst
c_function  To compute the Greenwich Mean Sidereal Time
c           for the time given by JULDAY. This is low precision in that
c           dUT1 is assumed to be 0.
c_call    st0 = gmst(julday)
c_/julday r*8  Double precision gives errors of up to 3 microseconds.
c_author   D McConnell
c_created  1-AUG-1988
c_latest   1-AUG-1988
c-
      double precision jd0, time, tu

      jd0 = (julday - dmod((julday+0.5d0),1.0d0))
      time = julday - jd0
      tu = (jd0 - 2451545.0d0)/36525.0d0
      gmst = dmod(0.279057273d0 + tu*(100.0021390
     -                                     + tu*(1.077592d-6
     -                                      - tu*7.2d-11)),1.0d0)
      gmst = dmod((101.0d0+ gmst + time*1.00273 79093 50795),1.0d0)
      return
      end

      subroutine precsn(Js,Je,P)
c+
c_name     PRECSN 
c_function To calculate the precession matrix P for
c          dates AFTER 1984.0 (JD = 2445700.5)
c          Given the position of an object referred
c          to the equator and equinox of the epoch Js
c          its position referred to the equator and equinox
c          of epoch Je can be calculated as follows :
c 
c  express position as directcosine 3-vector, V1
c  then the corresponding vector V2 for epoch Js is
c 
c             V2 = P.V1
c
c_call precsn(Js, Je, P)
c
c_/Js       i*4      Julian day number of starting epoch
c_/Je       i*4      Julian day number of ending epoch
c_/P        r*8(3,3) Precession matrix
c
c_author    D. McConnell
c_created   9-Mar-1984
c_latest    ???
c_refe      Control system, COORD
c-

      implicit none
      double precision      Je,     Js,     P(3,3), T,      T2,
     -                      st,    st2,     st3,    zeta,   z,
     -                      theta, coszet,  sinzet, costhe, sinthe,
     -                      cosz,  sinz,    J2000,  julcen,
     -                      a1, a2, a3, a4, a5, a6,
     -                      b1, b2, b3, b4, b5, b6,
     -                      c1, c2, c3, c4, c5, c6
 
      parameter            (J2000  =   245 1545.0 D00,
     -                      JULCEN =   36525.0 D00)
 
      parameter            (A1     =   0.011 180 860 865 024 D0,
     -                      A2     =   0.000 006 770 713 945 D0,
     -                      A3     = - 0.000 000 000 673 891 D0,
     -                      A4     =   0.000 001 463 555 541 D0,
     -                      A5     = - 0.000 000 001 672 607 D0,
     -                      A6     =   0.000 000 087 256 766 D0,
 
     -                      B1     =   0.011 180 860 865 024 D0,
     -                      B2     =   0.000 006 770 713 945 D0,
     -                      B3     = - 0.000 000 000 673 891 D0,
     -                      B4     =   0.000 005 307 158 404 D0,
     -                      B5     =   0.000 000 000 319 977 D0,
     -                      B6     =   0.000 000 088 250 634 D0,
 
     -                      C1     =   0.009 717 173 455 170 D0,
     -                      C2     = - 0.000 004 136 915 141 D0,
     -                      C3     = - 0.000 000 001 052 046 D0,
     -                      C4     =   0.000 002 068 457 570 D0,
     -                      C5     =   0.000 000 001 052 046 D0,
     -                      C6     = - 0.000 000 202 812 107 D0)
 
 
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
      T  = (Js - J2000)/JULCEN
      st = (Je - Js)/JULCEN
      T2 = T * T
      st2 = st * st
      st3 = st2 * st
 
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c          Calculate the Equatorial precession parameters
c              (ref.   USNO Circular no. 163      1981,
c                    Lieske et al., Astron. & Astrophys., 58, 1 1977)
 
      zeta  =   (A1 + A2*T + A3*T2) * st
     -        + (A4 + A5*T)         * st2
     -        +  A6                 * st3
 
      z     =   (B1 + B2*T + B3*T2) * st
     -        + (B4 + B5*T)         * st2
     -        +  B6                 * st3
 
      theta =   (C1 + C2*T + C3*T2) * st
     -        - (C4 + C5*T)         * st2
     -        +  C6                 * st3
 
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c          Calculate the P matrix

      cos zet = dcos (zeta)
      sin zet = dsin (zeta)
      cos z   = dcos (z)
      sin z   = dsin (z)
      cos the = dcos (theta)
      sin the = dsin (theta)
 
      P(1,1) =  cos zet * cos z * cos the - sin zet * sin z
      P(1,2) = -sin zet * cos z * cos the - cos zet * sin z
      P(1,3) = -cos z * sin the
 
      P(2,1) =  cos zet * sin z * cos the + sin zet * cos z
      P(2,2) = -sin zet * sin z * cos the + cos zet * cos z
      P(2,3) = -sin z * sin the
 
      P(3,1) =  cos zet * sin the
      P(3,2) = -sin zet * sin the
      P(3,3) =  cos the
 
      return
      end
c=====================================================================


        subroutine getradec(v)
        real*8 v(3),ra,dec
        real*8 rah,ram,decd,decm
        real*8 ras,decs,fracs1,fracs2
        real*8 deg2rad
        integer*4 decsign

        deg2rad = dacos(-1.0d0)/180.d0

        decsign = 1
        dec = dasin(v(3))/deg2rad
        ra = datan2(v(2),v(1))/deg2rad
        if (ra.lt.0.d0) ra = ra+360.d0
        if (dec.lt.0.d0) decsign = -1
        ra = ra/15.d0
 
        rah = int(ra)
        ram = int((ra - rah)*60.d0)
        ras = (ra - rah - ram/60.d0)*3600.d0
        fracs1 = 100.*(ras - int(ras))
 
        dec = abs(dec)
        decd = int(dec)
        decm = int((dec-decd)*60.d0)
        decs = (dec - decd -decm/60.d0)*3600.d0
        decd = decd*decsign
        fracs2 = 100.*(decs - int(decs))
 
        write(*,"('RA: ',I3,':',I2.2,':',I2.2,'.',I2.2)")
     +        int(rah),int(ram),int(ras),int(fracs1)
        write(*,"('Dec:',I3,':',I2.2,':',I2.2,'.',I2.2)")
     +        int(decd),int(decm),int(decs),int(fracs2)
 
        return
        end
 

