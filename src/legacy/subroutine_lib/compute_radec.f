cc   compute apparent ra, dec from az, el, lst
cc   Zombeck 1982 p. 71
cc   translated from Jeff Hagen's routine radec.c

cc   az, el, ra, dec in degrees
cc   lst in hours


cc  az za in degrees lst in hours 

      subroutine compute_radec( az, za, lst, pra, pdec )

      implicit none
      real*8  az, za, lst
      real*8  pra, pdec
      real*8  sindec, sinaz, sinel, cosaz, cosel
      real*8  cosdec, dec, ha, cosha, sinlat, coslat
      real*8  el, deg2rad
      real*8  AOLAT/18.353806/         ! degrees north 

      deg2rad = dacos(-1.0d0)/180.0d0

      az = az + 180.
      if( az.gt.360.0 )then
       az = az - 360.0
      end if
      el = 90.0 - za
      sinel = sin(el*deg2rad)
      cosel = cos(el*deg2rad)
      sinaz = sin(az*deg2rad)
      cosaz = cos(az*deg2rad)
      sinlat = sin(AOLAT*deg2rad)
      coslat = cos(AOLAT*deg2rad)

      sindec = sinel * sinlat + cosel * coslat * cosaz

      if( sindec.ge.-1.0.and.sindec.le.1.0)then
        dec = asin( sindec )/deg2rad
      else
        dec = 90.0
      end if

      cosdec = cos(dec*deg2rad)
      if(cosdec.ne.0.0)then
        cosha = ( sinel * coslat - cosel * cosaz * sinlat)/cosdec 
      else
        cosha = 1.0
      end if
      if(cosha.ge.-1.0.and.cosha.le.1.0)then
        ha = acos( cosha )/(15.0*deg2rad)
      else
        ha = 0.0
      end if

cc    Resolve quadrant ambiguity

      if( az.lt.180.0 )then
        ha = -ha
      end if 
      pra =  lst - ha
      if(pra.lt.0.0 )then
        pra = pra + 24.0
      else if (pra.gt.24.0)then
        pra = pra - 24.0
      end if
      pdec = dec
      return
      end
