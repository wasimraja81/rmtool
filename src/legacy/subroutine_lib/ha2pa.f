      real Function ha2pa(alat_deg,dec_deg,ha_hrs)
      !alat = -42.805        ! latitude of the 14m dish in deg.
      !dec  = -45.33             ! declination of Vela psr
      implicit none

      real*8 x, tanpa, dec_deg,
     -       alat_deg, pi, deg2rad
      real*4 ha_hrs
      

      pi = acos(-1.0)
      deg2rad = pi/180.0

      x = ha_hrs/24.*360.

      tanpa = (-sin(x*deg2rad)*cos(alat_deg*deg2rad)) /
     -        (cos(x*deg2rad)*cos(alat_deg*deg2rad)*sin(dec_deg*deg2rad)
     -         -    cos(dec_deg*deg2rad)*sin(alat_deg*deg2rad))
      !ha2pa = atan(tanpa)/deg2rad
      ha2pa = atan2(tanpa,1.0)/deg2rad
      if(ha2pa.gt.180.0)then
              ha2pa = ha2pa - 180.0
      else if(ha2pa.lt.0.0)then
              ha2pa = ha2pa + 180.0
      endif
c      type*,ha,x,ha2pa
      return
      end

