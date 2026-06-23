chelp+
*      Code to demonstrate usage of the subroutine 
*      azza2lst.f to compute the HA and the LST  
*      given a source's position in terms of Az, 
*      Za and RA.
*      Az and Za are expected in Degrees, while RA 
*      should be supplied in Hours. 
*      
*      The subroutine azza2lst.f was tailormade for 
*      Arecibo Observatory. We use it for any other
*      observatory: changes were made at two places:
*      1) The latitude of the Observatory
*      2) The convention for defining the Azimuth 
*         viz., 0 or 180 (Arecibo Az is 180 degrees 
*         away from the azel dish)
chelp-
*                  -- wr, 06 April, 2010

      implicit none


      integer*4     maxspec, nspec, i
      integer*4     site_code
      ! site-codes: 
      ! 1 --> AO
      ! 2 --> GMRT
      ! 3 --> GBT
      ! 4 --> SUBARU
      ! default --> AO
      parameter     (maxspec = 655360)
      real*8        LAT
      real*8        ra_h, dec_d, az, za, ha, lst
      integer*4     ra_hh, ra_mm, dec_dd, dec_mm
      integer*4     nchar
      real*8        ra_ss, dec_ss, pi, this_jdtime,
     -              new_jdtime, ra_radians, dec_radians,
     -              ra_out_radians, dec_out_radians
      real*8        tsamp, lst_now 
      real*4        ha_now(maxspec), pa(maxspec)
      real*4        ha2pa
      character*120 xlabel, ylabel, title, templine
      real*8        AOLAT/18.353806/        ! degrees north 
      real*8        GMRTLAT/19.360/         ! degrees north 
      real*8        GBTLAT/38.434317/       ! degrees north 
      real*8        SUBARULAT/19.8255/      ! degrees north 

      real*4        ha_start, ha_end, ha_test(maxspec)
      real*8        dec_test


      pi = acos(-1.0)

      !--------------------------------
      ! PARAMETERS to be set up:
      this_jdtime = 2000.0d0
      new_jdtime = dble(2008.0) + dble(11.0/12.0)
      !new_jdtime = this_jdtime
      site_code = 1

      
      az = 292.787800d0
      za =   9.888200d0
      ra_hh = 6
      ra_mm = 59
      ra_ss = 48.124000

      dec_dd = 14
      dec_mm = 14
      dec_ss = 21.530000

      tsamp = 0.002
      nspec = 299999
      !--------------------------------
      ha_start = -6.0
      ha_end = 6.0
      !dec_test = 14.239d0
      dec_test = 64.6d0
      write(*,*)"start HA, end HA, test DEC, site-code"
      read(*,*)ha_start,ha_end,dec_test,site_code
      call linspace(ha_start,ha_end,nspec,ha_test)


      if(site_code.eq.1)then
              LAT = AOLAT
      else if(site_code.eq.2)then
              LAT = GMRTLAT
      else if(site_code.eq.3)then
              LAT = GBTLAT
      else if(site_code.eq.4)then
              LAT = SUBARULAT
      else
              !use default site (AO)
              LAT = AOLAT
      endif

!      ! TEST:
!      ha_now(1) = 6.0
!      dec_d = 43.0 
!      pa(1) = ha2pa(LAT,dec_d,ha_now(1))
!      write(*,*)"PA: ",pa(1)
!      stop
!      !END TEST:




      ra_h = dble(ra_hh) + dble(ra_mm)/60.0 + dble(ra_ss)/3600.0
      dec_d = dble(dec_dd) + dble(dec_mm)/60.0 + dble(dec_ss)/3600.0


      if(ra_h.lt.0.0 )then
        ra_h = ra_h + 24.0
      else if (ra_h.gt.24.0)then
        ra_h = ra_h - 24.0
      end if
      ra_radians = ra_h*15.0*pi/180.0
      dec_radians = dec_d*pi/180.0

      write(*,*)"------------------- "
      write(*,*)"Parameters for epoch: ",this_jdtime
      write(*,*)"RA_h: ",ra_h
      write(*,*)"  RA: ",ra_hh,':',ra_mm,':',ra_ss
      write(*,*)" DEC: ",dec_dd,':',dec_mm,':',dec_ss
      write(*,*)"DECD: ",dec_d
      write(*,*)" "
      write(*,*)"------------------- "

      ! Convert to ra-now (current epoch):
      call precess(this_jdtime,new_jdtime,
     -             ra_radians,dec_radians,
     -             ra_out_radians,dec_out_radians)

      ra_h = (ra_out_radians*180.0/pi)/15.0
      dec_d = (dec_out_radians*180.0/pi)

      if(ra_h.lt.0.0 )then
        ra_h = ra_h + 24.0
      else if (ra_h.gt.24.0)then
        ra_h = ra_h - 24.0
      end if
      ra_hh = int(ra_h)
      ra_mm = int((ra_h - ra_hh)*60.0)
      ra_ss = ((ra_h - ra_hh)*60.0 - ra_mm)*60.0

      dec_dd = int(dec_d)
      dec_mm = int((dec_d - dec_dd)*60.0)
      dec_ss = ((dec_d - dec_dd)*60.0 - dec_mm)*60.0

      write(*,*)"------------------------- "
      write(*,*)"Parameters for epoch: ",new_jdtime
      write(*,*)"  "
      write(*,*)"  Az: ",az
      write(*,*)"  Za: ",za
      write(*,*)"RA_h: ",ra_h
      write(*,*)"  RA: ",ra_hh,':',ra_mm,':',ra_ss
      write(*,*)" DEC: ",dec_dd,':',dec_mm,':',dec_ss
      write(*,*)"DECD: ",dec_d
      write(*,*)" "
      write(*,*)"------------------------- "


      call azza2lst(az,za,ra_h,site_code,dec_d,ha,lst)


      write(*,*)"------------------------- "
      write(*,*)"OUTPUTS of azza2lst:  "
      dec_dd = int(dec_d)
      dec_mm = int((dec_d - dec_dd)*60.0)
      dec_ss = ((dec_d - dec_dd)*60.0 - dec_mm)*60.0
      write(*,*)" "
      write(*,*)"  HA: ",ha
      write(*,*)"RA_h: ",ra_h
      write(*,*)" LST: ",lst
      write(*,*)" DEC: ",dec_dd,':',dec_mm,':',dec_ss
      write(*,*)"DECD: ",dec_d
      write(*,*)" "
      write(*,*)"------------------------- "

      ! Now we shall simulate the PA curve for the 
      ! duration of observation:
      
      lst_now = lst
      do i = 1,nspec
         lst_now = lst_now + tsamp*366.2422d0/365.2422d0/3600.0d0
         ha_now(i) = lst_now - ra_h
         pa(i) = ha2pa(LAT,dec_d,ha_now(i))
      enddo
      xlabel = 'Hour Angle'
      ylabel = 'PA '
      title = 'PA vs HA for B0656+14 data'
      call pgbeg(0,'/xs',1,1)
      call myplot1(ha_now,pa,nspec,xlabel,ylabel,title,2)

      ! if test plots:
      do i = 1,nspec
         pa(i) = ha2pa(LAT,dec_test,ha_test(i))
      enddo
      write(templine,'(f5.2)')dec_test
      title = 'PA vs HA, Dec: '//templine(1:nchar(templine))//'- site: '
      write(templine,'(i1)')site_code
      title = title(1:nchar(title))//templine(1:nchar(templine))
      call myplot1(ha_test,pa,nspec,xlabel,ylabel,title,2)

      call pgend

      end
      include 'azza2lst.f'
      include 'precess.f'
      include 'nchar.f'
      include 'myplot1.f'
      include 'fort_lib.f'
      include 'ha2pa.f'
