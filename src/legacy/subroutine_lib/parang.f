chelp+
c     **********************************************************
c     subroutine parang(za,az,feedth,zach,xp)
c     **********************************************************
c     parang computes the parallactic angle for an offset feed.
c     All inputs are in degrees. (za,az) are source coordinates,
c     feedth is theta from the feeds file, zach is the carriage
c     house zenith angle, and xp is the parallactic angle.
c     -- Tony Phillips,  April 1989
c     **********************************************************
c
chelp-
      subroutine parang(za,az,feedth,zach,xp)


      real lat
      parameter (pi=3.141592653589)
c
c     Arecibo latitude.
      data lat/18.34/                
c
      lat=lat*(pi/180.)
      sx=-sin(zach)*sin(feedth)/sin(za)
      x=asin(sx)
      top=sin(az)
      bottom=cos(az)*cos(za)-tan(lat)*sin(za)
      xp=atan(top/bottom)
      xp=xp-x
c     write(*,15) xp
 15   format(' pa = ',e9.2)
      return
      end
c--------------------------------------------------------------------
