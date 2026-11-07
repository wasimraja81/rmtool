chelp+
*      ---------------------------------------------      
*      subroutine azza2lst(az,za,ra,site,
*                        ->dec,ha,lst)
*       In: az(deg), za(deg), ra(hrs), site_code
*      Out: dec(deg), ha(hrs), lst(hrs)
*       site-codes: 
*                    1 --> AO
*                    2 --> GMRT
*                    3 --> GBT
*                    4 --> SUBARU
*              default --> AO
*      
*      Code to compute the HA and the LST given a 
*      source's position in terms of Az, Za and RA.
*      
*      The subroutine azza2lst.f was tailormade for 
*      Arecibo Observatory. We use it for any other
*      observatory, by changing :
*      1) The latitude of the Observatory
*      2) The convention for defining the Azimuth 
*         viz., 0 or 180 (Arecibo Az is 180 degrees 
*         away from the azel dish)
*      

*      ---------------------------------------------      
chelp-
*                  -- wr, 06 APRIL, 2010
*
cc     Reference: Zombeck 1982 p. 71
cc     translated from Jeff Hagen's routine radec.c

cc     az za in degrees; lst,ra,ha in hours 


      subroutine azza2lst(az, za, ra, site_code, pdec, ha, lst)

      implicit none
      real*8  az, za, lst
      real*8  ra, pdec
      real*8  sindec, sinaz, sinel, cosaz, cosel
      real*8  cosdec, dec, ha, cosha, sinlat, coslat
      real*8  el, deg2rad
      integer*4  site_code   
      real*8 LAT
      real*8 AOLAT/18.353806/         ! degrees north 
      real*8 GMRTLAT/19.360/         ! degrees north 
      real*8 GBTLAT/38.434317/         ! degrees north 
      real*8 SUBARULAT/19.8255/         ! degrees north 

      deg2rad = dacos(-1.0d0)/180.0d0

      if(site_code.eq.1)then
              az = az + 180.
              LAT = AOLAT
      else if(site_code.eq.2)then
              az = az
              LAT = GMRTLAT
      else if(site_code.eq.3)then
              az = az
              LAT = GBTLAT
      else if(site_code.eq.4)then
              az = az
              LAT = SUBARULAT
      else
              !use default site (AO)
              LAT = AOLAT
              az = az + 180.0
      endif

      if( az.gt.360.0 )then
       az = az - 360.0
      end if

      el = 90.0 - za
      sinel = sin(el*deg2rad)
      cosel = cos(el*deg2rad)
      sinaz = sin(az*deg2rad)
      cosaz = cos(az*deg2rad)
      sinlat = sin(LAT*deg2rad)
      coslat = cos(LAT*deg2rad)

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
      if(ra.lt.0.0 )then
        ra = ra + 24.0
      else if (ra.gt.24.0)then
        ra = ra - 24.0
      end if
      lst =  ra + ha
      pdec = dec
      return
      end
