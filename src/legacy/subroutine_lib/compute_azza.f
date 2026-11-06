cc   returns the azimuth and za for a source at apparent ra, dec and lst
cc   Zombeck 1982 p. 71
cc   translated from Jeff Hagen's routine azza.c

cc   az, el, ra, dec in degrees
cc   lst in hours


cc  az za in degrees, lst in hours 

      subroutine compute_azza(ra, dec, lst, paz, pza) 
      implicit none

      real*8  ra, paz, pza, lst, az, el
      real*8  sindec, sinel, cosaz, cosel, cosazel
      real*8  cosdec, dec, ha, cosha, sinlat, coslat
      real*8  deg2rad
      real*8  AOLAT/18.353806/         ! degrees north 

      deg2rad = dacos(-1.0d0)/180.0d0


cc    Calculate hour angle and convert to degrees (use convention
cc    that hour angle ranges between +/- 12 hours)   

      ha = lst - ra
      if( ha.gt.12.0 )then
        ha = ha - 24.0
      else if( ha.lt.-12.0)then
        ha = ha + 24.0
      end if
      ha = ha*15.0d0   ! now in degrees

cc    Compute trig quantities needed; make angle calculations double
cc    precision to cut down on round-off-induced overflows near the
cc    zenith and the meridian.    

      sindec = sin(dec*deg2rad)
      cosdec = cos(dec*deg2rad)
      sinlat = sin(AOLAT*deg2rad)
      coslat = cos(AOLAT*deg2rad)
      cosha  = cos(ha*deg2rad)

cc    Compute elevation (trap for blow-ups near the zenith):  

      sinel = sindec * sinlat + cosdec * cosha * coslat

      if(sinel.gt.-1.0.and.sinel.le.1.0)then
        el = asin(sinel)/deg2rad
      else
        write(*,*)'elevation too high !!'
      end if

cc    Compute azimuth (trap for blow-ups near HA = 0)    

      cosel = cos(el*deg2rad)
      cosazel = sindec * coslat - cosdec * cosha * sinlat
      cosaz = (1.0/cosel) * cosazel

      if (cosaz.ge.-1.0.and.cosaz.le.1.0)then
        az = acos(cosaz)/deg2rad
      else 
        if(dec.le.AOLAT)then
          az = 180.0
        else
          az = 0.0
        end if
      end if

cc    Resolve quadrant ambiguity of azimuth:    

      if (ha.gt.0.0)then
        az = 360.0 - az
      end if
      pza = 90.0 - el

      az = az + 180.0   ! arecibo az is 180 away from azel dish 
      if( az.gt.360.0)then
        az = az - 360.0
      end if
      paz = az

cc       printf("az %f dec %f ha %f\n", az, dec, ha ); */

      return
      end

