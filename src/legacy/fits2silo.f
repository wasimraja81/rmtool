chelp+ 
      !--------------------------------------------------------------
      ! This code was developed to export FITS images to SILO. 
      ! SILO is the native file-format for VisIt. We intend to 
      ! use VisIt for various kinds of visualization of our 
      ! data. 
      !
      ! Currently we intend to update the SILO file with the 
      ! following information: 
      ! 
      ! 1. spectral averaged Stokes-I     : <I>            [load] 
      ! 2. spectral averaged Stokes-Q     : <Q>            [load] 
      ! 3. spectral averaged Stokes-U     : <U>            [load] 
      ! 4. spectral averaged Stokes-V     : <V>            [load] 
      ! 5. LP @ RM = i                    : LP(RM=rm_val)  [load] 
      ! 6. PA @ RM = i                    : PA(RM=rm_val)  [load] 
      ! 7. dLP = LP/<I>                   : dLP(RM=rm_val) [compute] 
      ! 8. Q @ RM = i [ = LP x cos(PA) ]  : Q(RM=rm_val)   [compute] 
      ! 9. U @ RM = i [ = LP x sin(PA) ]  : U(RM=rm_val)   [compute] 
      !10. RA                             : RA_pix          
      !11. Dec                            : Dec_pix         
      !12. SNR-I = <I>/rmsI (spectral)    : SNR_I          [load] 
      !13. XRAY image                     : XRAY_MAP       [load] 
      ! 
      !                                  --wr, 09 Jul, 2012
      !--------------------------------------------------------------
chelp- 
      ! Last modification: wr, 16 Jul, 2012.
      ! Last modification: wr, 31 Jul, 2012.
      !
      !---------------------------------------------------------


      implicit none 

      integer*4         max_axes, maxdimx, maxdimy, maxunit 
      parameter         (max_axes=99,maxdimx = 1024, maxdimy = 1024,
     -                   maxunit = 99 )

      integer*4         nchar 
      integer*4         nused, nbad, use_filter 
      real*4            min_I_cutoff, max_I_cutoff 
      integer*4         ix, iy, i, k 
      character*220     cfgfile 
      character*220     infile_I, infile_Q, infile_U, infile_V 
      character*220     infile_lp, infile_pa, infile_snr_I, infile_X 
      character*220     infile_Irms, infile_Qrms, 
     -                  infile_Urms, infile_Vrms, infile_filter 
      character*220     outfile  
      character*220     path, out_path 
      character*1       junkchar 
      integer*4         cxpix, cypix, nxpix, nypix 
      !integer*4         cxpix1, cypix1, nxpix1, nypix1 
      integer*4         cxpix2, cypix2, nxpix2, nypix2 


      real*4            I_arr(maxdimx*maxdimy), 
     -                  Q_arr(maxdimx*maxdimy), 
     -                  U_arr(maxdimx*maxdimy), 
     -                  V_arr(maxdimx*maxdimy), 
     -                  X_arr(maxdimx*maxdimy), 
     -                  lp_arr(maxdimx*maxdimy), 
     -                  pa_arr(maxdimx*maxdimy), 
     -                  Irms_arr(maxdimx*maxdimy), 
     -                  Qrms_arr(maxdimx*maxdimy), 
     -                  Urms_arr(maxdimx*maxdimy), 
     -                  Vrms_arr(maxdimx*maxdimy), 
     -                  snr_I_arr(maxdimx*maxdimy), 
     -                  RA_arr(maxdimx*maxdimy), 
     -                  Dec_arr(maxdimx*maxdimy), 
     -                  x(maxdimx*maxdimy), 
     -                  y(maxdimx*maxdimy), 
     -                  z(maxdimx*maxdimy), 
     -                  indx(maxdimx*maxdimy), 
     -                  indy(maxdimx*maxdimy), 
     -                  image_arr(maxdimx,maxdimy), 
     -                  image_arr2(maxdimx,maxdimy), 
     -                  dLP_arr(maxdimx*maxdimy), 
     -                  Q_rm_arr(maxdimx*maxdimy), 
     -                  U_rm_arr(maxdimx*maxdimy)

      real*4            qmean, umean, qtmp, utmp, linpol, polang 
      integer*4         err, ierr, ndims, dbfile 

      integer*4         iunit
      integer*4         status  
      character         out_name*72, out_class*16 
      character         templine*220 
      character*6       min_I_str, max_I_str 
      integer*4         pixtype 

      integer*4    IstokesOpts, QstokesOpts, UstokesOpts, VstokesOpts, 
     -             LinpolOpts, PolangleOpts, snrOpts, dLpolOpts, 
     -             RaOpts, DecOpts, QrmOpts, UrmOpts, 
     -             xrayOpts, IrmsOpts, QrmsOpts, UrmsOpts, VrmsOpts

      character*16 I_UNITS, Q_UNITS, U_UNITS, V_UNITS, 
     -             LP_UNITS, PA_UNITS, SNR_UNITS, dLP_UNITS, 
     -             RA_UNITS, Dec_UNITS, Q_rm_UNITS, U_rm_UNITS, 
     -             xray_UNITS, IRMS_UNITS, QRMS_UNITS, URMS_UNITS, 
     -             VRMS_UNITS 

      integer*4    lenI_UNITS, lenQ_UNITS, lenU_UNITS, lenV_UNITS, 
     -             lenLP_UNITS, lenPA_UNITS, lenSNR_UNITS, lendLP_UNITS,
     -             lenRA_UNITS, lenDec_UNITS, lenQ_rm_UNITS,
     -             lenU_rm_UNITS, lenxray_UNITS, lenIRMS_UNITS, 
     -             lenQRMS_UNITS, lenURMS_UNITS, lenVRMS_UNITS 
      integer*4    len_outfile, len_mesh_name 
      character*16 mesh_name 

      include '/home/wasim/SILO/SILO_ROOT/include/silo.inc'
      
      !--------------------------------------
      ! Some input parameters: 
      if(iargc().ne.1)then 
              write(*,*)"Usage: "
              write(*,*)"You need to use a config file: "
              write(*,*)"    comb_fits_image <config file> "
              write(*,*)" "
              stop 
      else
              call getarg(1,cfgfile) 
              cfgfile = '../CONFIG/'//cfgfile(1:nchar(cfgfile))
      endif


      call get_lun(iunit) 
      open(iunit,file=cfgfile,status='old',err=101)
      goto 102
101   write(*,*)"Error opening file: ",cfgfile(1:nchar(cfgfile))
      write(*,*)"Quitting now..."
      stop 

102   continue 
      read(iunit,*)junkchar   ! comment line 
      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      path = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_I = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_Q = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_U = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_V = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_lp = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_pa = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_Irms = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_Qrms = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_Urms = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_Vrms = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_snr_I = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      infile_X = templine(1:nchar(templine)) 

      !---------------------------------------------------
      read(iunit,*)use_filter 

      read(iunit,*)min_I_cutoff, max_I_cutoff 

      write(templine,'(f7.3)')min_I_cutoff
      ! Count the number of preceding SPACE characters: 
      i = 0 
      do k = 1,nchar(templine)
         if (templine(k:k).eq.' ')then 
                 i = i + 1 
         else
                 goto 1100
         endif
      enddo
1100  continue 
      min_I_str = templine(i+1:nchar(templine))

      write(templine,'(f7.3)')max_I_cutoff
      ! Count the number of preceding SPACE characters: 
      i = 0 
      do k = 1,nchar(templine)
         if (templine(k:k).eq.' ')then 
                 i = i + 1 
         else
                 goto 1101
         endif
      enddo
1101  continue 
      max_I_str = templine(i+1:nchar(templine))

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      out_path = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      out_name = templine(1:nchar(templine)) 

      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 

      ! Overwrite outclass if specified in .cfg file: 
      if(nchar(templine).gt.0)then 
              out_class = templine(1:nchar(templine)) 
      endif


      read(iunit,'(a)')templine 
      templine = templine(1:index(templine,';')-1) 
      mesh_name = templine(1:nchar(templine)) 
      len_mesh_name = nchar(mesh_name) 

      read(iunit,*)pixtype 
      close(iunit) 




      infile_I = path(1:nchar(path))//infile_I(1:nchar(infile_I))
      infile_Q = path(1:nchar(path))//infile_Q(1:nchar(infile_Q))
      infile_U = path(1:nchar(path))//infile_U(1:nchar(infile_U))
      infile_V = path(1:nchar(path))//infile_V(1:nchar(infile_V))
      infile_X = path(1:nchar(path))//infile_X(1:nchar(infile_X))
      infile_LP = path(1:nchar(path))//infile_LP(1:nchar(infile_LP))
      infile_PA = path(1:nchar(path))//infile_PA(1:nchar(infile_PA))
      infile_Irms = path(1:nchar(path))//
     -              infile_Irms(1:nchar(infile_Irms))
      infile_Qrms = path(1:nchar(path))//
     -              infile_Qrms(1:nchar(infile_Qrms))
      infile_Urms = path(1:nchar(path))//
     -              infile_Urms(1:nchar(infile_Urms))
      infile_Vrms = path(1:nchar(path))//
     -              infile_Vrms(1:nchar(infile_Vrms))
      infile_snr_I = path(1:nchar(path))//
     -               infile_snr_I(1:nchar(infile_snr_I))

!      write(*,*)"infile_I: ",infile_I(1:nchar(infile_I))
!      write(*,*)"infile_Q: ",infile_Q(1:nchar(infile_Q))
!      write(*,*)" "

      if (pixtype.eq.0)then
              outfile = out_path(1:nchar(out_path))// 
     -                  out_name(1:nchar(out_name))//'.'//
     -                  out_class(1:nchar(out_class))//
     -                  '.SILO'
      else if (pixtype.lt.0)then
              outfile = out_path(1:nchar(out_path))// 
     -                   out_name(1:nchar(out_name))//'.'//
     -                   out_class(1:nchar(out_class))//
     -                   '.NOISE_REGION.Icutoff_'//
     -                   min_I_str(1:nchar(min_I_str))//
     -                   '.SILO'
      else if(pixtype.gt.0.and.pixtype.lt.2)then
              outfile = out_path(1:nchar(out_path))// 
     -                   out_name(1:nchar(out_name))//'.'//
     -                   out_class(1:nchar(out_class))//
     -                   '.SIGNIFICANT_I_REGION.Icutoff_'//
     -                   min_I_str(1:nchar(min_I_str))//
     -                   '.SILO'
      else if(pixtype.ge.2)then
              outfile = out_path(1:nchar(out_path))// 
     -                   out_name(1:nchar(out_name))//'.'//
     -                   out_class(1:nchar(out_class))//
     -                   '.I_REGION_RANGE_'//
     -                   min_I_str(1:nchar(min_I_str))//'_to_'//
     -                   max_I_str(1:nchar(max_I_str))//
     -                   '.SILO'
      else
              ! Default -- ALL pixels will be used: 
              write(*,*)"Invalid pixtype..."
              write(*,*)"Using ALL pixels..."
              outfile = out_path(1:nchar(out_path))// 
     -                  out_name(1:nchar(out_name))//'.'//
     -                  out_class(1:nchar(out_class))//
     -                  '.SILO'
      endif


      !write(*,*)"min_I_str: ",min_I_str(1:nchar(min_I_str))
      !write(*,*)"max_I_str: ",max_I_str(1:nchar(max_I_str))
      write(*,*)"outfile: ",outfile(1:nchar(outfile))
      write(*,*)" "
      !stop
      
      
   
      !-------------------------------------- 
      ! We wish to load the entire images into 
      ! memory (Ensure that the images are NOT 
      ! very large): 
      !=======================================
      ! Freeze parameters to match criterion for 
      ! ENTIRE image reading: 
      cxpix = 0 
      cypix = 0 

      cxpix2 = 0 
      cypix2 = 0 
      
      nxpix = 0 
      nypix = 0 
      !nxpix1 = 0 
      !nypix1 = 0 
      nxpix2= 0 
      nypix2= 0 
      !=======================================


      !----------------------------------------------------
      ! READ the necessary FITS files to be used for filtering:  
      if (use_filter .ne. 13)then
         if (use_filter .eq.1)then 
              infile_filter = infile_I(1:nchar(infile_I))
         else if(use_filter .eq. 2)then
              infile_filter = infile_Q(1:nchar(infile_Q))
         else if(use_filter .eq. 3)then
              infile_filter = infile_U(1:nchar(infile_U))
         else if(use_filter .eq. 4)then
              infile_filter = infile_V(1:nchar(infile_V))
         else if(use_filter .eq. 5)then
              infile_filter = infile_LP(1:nchar(infile_LP))
         else if(use_filter .eq. 6)then
              infile_filter = infile_PA(1:nchar(infile_PA))
         else if(use_filter .eq. 7)then
              infile_filter = infile_Irms(1:nchar(infile_Irms))
         else if(use_filter .eq. 8)then
              infile_filter = infile_Qrms(1:nchar(infile_Qrms))
         else if(use_filter .eq. 9)then
              infile_filter = infile_Urms(1:nchar(infile_Urms))
         else if(use_filter .eq. 10)then
              infile_filter = infile_Vrms(1:nchar(infile_Vrms))
         else if(use_filter .eq. 11)then
              infile_filter = infile_snr_I(1:nchar(infile_snr_I))
         else if(use_filter .eq. 12)then
              infile_filter = infile_X(1:nchar(infile_X))
         endif
         call load_fits_image(infile_filter, cxpix2,cypix2,nxpix,nypix,
     -                  image_arr, maxdimx, maxdimy,status)
          if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_filter(1:nchar(infile_filter))
              write(*,*)"Quitting now..."
              stop 
          endif
      else 
          ! Compute dLP first and then use the filter: 
          call load_fits_image(infile_I, cxpix2,cypix2,nxpix,nypix,
     -                  image_arr, maxdimx, maxdimy,status)
          if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_I(1:nchar(infile_I))
              write(*,*)"Quitting now..."
              stop 
          endif
          call load_fits_image(infile_LP, cxpix2,cypix2,nxpix,nypix,
     -                  image_arr2, maxdimx, maxdimy,status)
          if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_LP(1:nchar(infile_LP))
              write(*,*)"Quitting now..."
              stop 
          endif
          do ix = 1,nxpix
             do iy = 1,nypix
                image_arr(ix,iy) = image_arr2(ix,iy)/image_arr(ix,iy) 
             enddo
          enddo
      endif
      nused = nxpix * nypix 
      if(pixtype.eq.0)then
              ! Use ALL pixels
              k = 0 
              i = 0 
              do ix = 1,nxpix
                 do iy = 1,nypix
                    k = k + 1 
                    indx(k) = ix 
                    indy(k) = iy 
                 enddo
              enddo
      else if(pixtype.gt.0.and.pixtype.lt.2)then
              ! Use pixels for which I(pix) > min_I_cutoff 
              k = 0 
              do ix = 1,nxpix
                 do iy = 1,nypix
                    if(image_arr(ix,iy).gt.min_I_cutoff)then 
                            k = k + 1 
                            indx(k) = ix 
                            indy(k) = iy 
                    endif
                 enddo
              enddo
              nused = k 
      else if(pixtype.lt.0)then
              ! Use pixels for which I(pix) <= min_I_cutoff 
              k = 0 
              do ix = 1,nxpix
                 do iy = 1,nypix
                    if(image_arr(ix,iy).le.min_I_cutoff)then 
                            k = k + 1 
                            indx(k) = ix 
                            indy(k) = iy 
                    endif
                 enddo
              enddo
              nused = k 
      else if(pixtype.ge.2)then
              ! Use pixels for which min_I_cutoff<= I(pix) <= max_I_cutoff 
              k = 0 
              do ix = 1,nxpix
                 do iy = 1,nypix
                    if(image_arr(ix,iy).ge.min_I_cutoff.and.
     -                 image_arr(ix,iy).le.max_I_cutoff)then 
                            k = k + 1 
                            indx(k) = ix 
                            indy(k) = iy 
                    endif
                 enddo
              enddo
              nused = k 
      endif
      if(k .lt. 30)then
              write(*,*)"Not enough data found!"
              write(*,*)"Quitting now..."
              stop
      endif
      ! End of data filtering based on input criterion
      ! --------------------------------------------------- 

      !----------------------------------------------------
      call load_fits_image(infile_I, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_I(1:nchar(infile_I))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif
      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         I_arr(k) = image_arr(ix,iy)
      enddo

      !----------------------------------------------------
      call load_fits_image(infile_Q, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_Q(1:nchar(infile_Q))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif
      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         Q_arr(k) = image_arr(ix,iy)
      enddo

      !----------------------------------------------------
      call load_fits_image(infile_U, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_U(1:nchar(infile_U))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif
      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         U_arr(k) = image_arr(ix,iy)
      enddo


      !----------------------------------------------------
      call load_fits_image(infile_V, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_V(1:nchar(infile_V))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         V_arr(k) = image_arr(ix,iy)
         ! Use mod(V) -- to assess leakage:  
         !V_arr(k) = abs(image_arr(ix,iy))
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_snr_I, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_snr_I(1:nchar(infile_snr_I))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         snr_I_arr(k) = image_arr(ix,iy)
      enddo

      call load_fits_image(infile_lp, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_lp(1:nchar(infile_lp))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         LP_arr(k) = image_arr(ix,iy) 
      enddo
   
      !----------------------------------------------------
      call load_fits_image(infile_pa, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_pa(1:nchar(infile_pa))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         PA_arr(k) = image_arr(ix,iy) 
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_Irms, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_Irms(1:nchar(infile_Irms))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         Irms_arr(k) = image_arr(ix,iy) 
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_Qrms, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_Qrms(1:nchar(infile_Qrms))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         Qrms_arr(k) = image_arr(ix,iy) 
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_Urms, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_Urms(1:nchar(infile_Urms))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         Urms_arr(k) = image_arr(ix,iy) 
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_Vrms, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_Vrms(1:nchar(infile_Vrms))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 
         Vrms_arr(k) = image_arr(ix,iy) 
      enddo
      !----------------------------------------------------
      call load_fits_image(infile_X, cxpix2,cypix2,nxpix2,nypix2,
     -                  image_arr, maxdimx, maxdimy,status)
      if(status.ne.0)then
              write(*,*)"Error opening file: ",
     -                        infile_X(1:nchar(infile_X))
              write(*,*)"Quitting now..."
              stop 
      endif
      ! Bare minimum check for dimensional mismatch: 
      if (nxpix.ne.nxpix2.or.nypix.ne.nypix2)then
              write(*,*)"Image dimensions do not match!"
              write(*,*)"Quitting now..."
              stop
      endif

      k = 0 
      nbad = 0 
      do k = 1,nused
         ix = indx(k) 
         iy = indy(k) 

         X_arr(k) = image_arr(ix,iy)
         ! Use scale if Xray is too high          
         !X_arr(k) = image_arr(ix,iy)*1.0e-3

         if(X_arr(k) .lt. 0.0)then
                 X_arr(k) = 0.0 
         endif
         ! Some more variables: 
         RA_arr(k) = real(ix) 
         Dec_arr(k) = real(iy) 

         ! RM Domain varables: 
         dLP_arr(k) = LP_arr(k) / I_arr(k) 
         !if(dLP_arr(k).ne.dLP_arr(k).or.LP_arr(k).gt.I_arr(k))then
         if(dLP_arr(k).ne.dLP_arr(k).or.LP_arr(k).gt.abs(I_arr(k)))then
                 dLP_arr(k) = 1.2 ! Dump noisy point to a high bin 
                 !dLP_arr(k) = 0.0 ! Dump noisy point 
                 nbad = nbad + 1 
         endif 
         Q_rm_arr(k) = LP_arr(k) * cos(PA_arr(k)) 
         U_rm_arr(k) = LP_arr(k) * sin(PA_arr(k)) 

         ! Prepare the mesh variables: 
         x(k) = real(k) 
         y(k) = real(k) 
         z(k) = real(k) 

      enddo
      ! Compute the mean of Q and U arrays: 
      ! [We may like to write out the mean-removed polarization vector]
      call mean (Q_rm_arr,nused,qmean) 
      call mean (U_rm_arr,nused,umean) 

      write(*,*)"Number of BAD dLpol pix: ",nbad 

      !----------------------------------------------------
      ! End of reading necessary INPUT files
      !----------------------------------------------------
      write(*,*)"Input files have "
      write(*,*)"    naxes(1): ", nxpix 
      write(*,*)"    naxes(2): ", nypix
      write(*,*)"Total pixels: ", nxpix * nypix 
      write(*,*)" "
      write(*,*)"Output file has: "
      write(*,*)"Total pixels: ", nused  
      write(*,*)"----------------------------"

      !--------------------------------------------- 
              
      ! Create the SILO FILE : 
      len_outfile = nchar(outfile)
      ierr = dbcreate(outfile, len_outfile,
     -                 DB_CLOBBER, DB_LOCAL, 
     -                 DB_F77NULL, 0, DB_HDF5, dbfile) 

      if(dbfile.eq.-1) then
              write (*,*) 'Could not create Silo file!\n'
              stop
      endif
      ! Generate the MESH for SILO: 
      ndims = 3 

      err = DBPutpm (dbfile, mesh_name, len_mesh_name, 
     -               ndims, x, y, z, nused, DB_FLOAT, 
     -               DB_F77NULL, ierr)

      !---------------------------------------------
      ! Write the values : 

      !---------------------------------------------
      ! 1. Stokes-I 
      I_units = 'Jy/bm'
      lenI_UNITS = nchar(I_units)

      err = dbmkoptlist(1,IstokesOpts)
      err = dbaddcopt(IstokesOpts,DBOPT_UNITS,I_UNITS,lenI_UNITS)

      err = dbputpv1 (dbfile, "specavg_I", 9,mesh_name,len_mesh_name, 
     -                I_arr, nused, DB_FLOAT, IstokesOpts, ierr)
      !---------------------------------------------
      ! 2. Stokes-Q 
      Q_units = 'Jy/bm'
      lenQ_UNITS = nchar(Q_units)

      err = dbmkoptlist(1,QstokesOpts)
      err = dbaddcopt(QstokesOpts,DBOPT_UNITS,Q_UNITS,lenQ_UNITS)

      err = dbputpv1 (dbfile, "specavg_Q", 9,mesh_name,len_mesh_name, 
     -                Q_arr, nused, DB_FLOAT, QstokesOpts, ierr)

      !---------------------------------------------
      ! 3. Stokes-U 
      U_units = 'Jy/bm'
      lenU_UNITS = nchar(U_units)

      err = dbmkoptlist(1,UstokesOpts)
      err = dbaddcopt(UstokesOpts,DBOPT_UNITS,U_UNITS,lenU_UNITS)

      err = dbputpv1 (dbfile, "specavg_U", 9,mesh_name,len_mesh_name, 
     -                U_arr, nused, DB_FLOAT, UstokesOpts, ierr)

      !---------------------------------------------
      ! 4. Stokes-V 
      V_units = 'Jy/bm'
      lenV_UNITS = nchar(V_units)

      err = dbmkoptlist(1,VstokesOpts)
      err = dbaddcopt(VstokesOpts,DBOPT_UNITS,V_UNITS,lenV_UNITS)

      err = dbputpv1 (dbfile, "specavg_V", 9,mesh_name,len_mesh_name, 
     -                V_arr, nused, DB_FLOAT, VstokesOpts, ierr) 

      !---------------------------------------------
      ! 5. LinPol 
      LP_units = 'Jy/bm'
      lenLP_UNITS = nchar(LP_units)

      err = dbmkoptlist(1,LinpolOpts)
      err = dbaddcopt(LinpolOpts,DBOPT_UNITS,LP_UNITS,lenLP_UNITS)

      err = dbputpv1 (dbfile, "LinPol",6,mesh_name,len_mesh_name, 
     -                LP_arr, nused, DB_FLOAT, LinpolOpts, ierr)

      !---------------------------------------------
      ! 6. PolAngle 
      PA_units = 'radians'
      lenPA_UNITS = nchar(PA_units)

      err = dbmkoptlist(1,PolangleOpts)
      err = dbaddcopt(PolangleOpts,DBOPT_UNITS,PA_UNITS,lenPA_UNITS)

      err = dbputpv1 (dbfile, "PolAngle",8,mesh_name,len_mesh_name, 
     -                PA_arr, nused, DB_FLOAT, PolangleOpts, ierr)

      !---------------------------------------------
      ! 7. dLinPol
      dLP_units = 'fraction'
      lendLP_UNITS = nchar(dLP_units)

      err = dbmkoptlist(1,dLPolOpts)
      err = dbaddcopt(dLPolOpts,DBOPT_UNITS,dLP_UNITS,lendLP_UNITS)

      err = dbputpv1 (dbfile, "dLinpol",7,mesh_name,len_mesh_name, 
     -                dLP_arr, nused, DB_FLOAT, dLpolOpts, ierr)

      !---------------------------------------------
      ! 8. Q(RM) 
      Q_rm_units = 'Jy/bm'
      lenQ_rm_UNITS = nchar(Q_rm_units)

      err = dbmkoptlist(1,QrmOpts)
      err = dbaddcopt(QrmOpts,DBOPT_UNITS,Q_rm_UNITS,lenQ_rm_UNITS)

      err = dbputpv1 (dbfile, "Q_RM",4,mesh_name,len_mesh_name, 
     -                Q_rm_arr, nused, DB_FLOAT, QrmOpts, ierr)

      !---------------------------------------------
      ! 9. U(RM) 
      U_rm_units = 'Jy/bm'
      lenU_rm_UNITS = nchar(U_rm_units)

      err = dbmkoptlist(1,UrmOpts)
      err = dbaddcopt(UrmOpts,DBOPT_UNITS,U_rm_UNITS,lenU_rm_UNITS)

      err = dbputpv1 (dbfile, "U_RM",4,mesh_name,len_mesh_name, 
     -                U_rm_arr, nused, DB_FLOAT, UrmOpts, ierr)

      !---------------------------------------------
      ! 10. SNR_I  
      SNR_units = 'ratio'
      lenSNR_UNITS = nchar(SNR_units)

      err = dbmkoptlist(1,snrOpts)
      err = dbaddcopt(snrOpts,DBOPT_UNITS,SNR_UNITS,lenSNR_UNITS)

      err = dbputpv1 (dbfile, "SNR_I",5,mesh_name,len_mesh_name, 
     -                snr_I_arr, nused, DB_FLOAT, snrOpts, ierr)

      !---------------------------------------------
      ! 11. RA  
      RA_units = 'PIXELS'
      lenRA_UNITS = nchar(RA_units)

      err = dbmkoptlist(1,RaOpts)
      err = dbaddcopt(RaOpts,DBOPT_UNITS,RA_UNITS,lenRA_UNITS)

      err = dbputpv1 (dbfile, "RA",2,mesh_name,len_mesh_name, 
     -                RA_arr, nused, DB_FLOAT, RaOpts, ierr) 

      !---------------------------------------------
      ! 12. DEC  
      Dec_units = 'PIXELS'
      lenDec_UNITS = nchar(Dec_units)

      err = dbmkoptlist(1,DecOpts)
      err = dbaddcopt(DecOpts,DBOPT_UNITS,Dec_UNITS,lenDec_UNITS)

      err = dbputpv1 (dbfile, "Dec",3,mesh_name,len_mesh_name, 
     -                Dec_arr, nused, DB_FLOAT, DecOpts, ierr) 



      !---------------------------------------------
      ! 13. X-RAY   
      XRAY_UNITS = 'counts'
      lenxray_UNITS = nchar(XRAY_UNITS)

      err = dbmkoptlist(1,xrayOpts)
      err = dbaddcopt(xrayOpts,DBOPT_UNITS,XRAY_UNITS,lenxray_UNITS)

      err = dbputpv1 (dbfile, "X_RAY",5,mesh_name,len_mesh_name, 
     -                X_arr, nused, DB_FLOAT, xrayOpts, ierr)


      !---------------------------------------------
      ! 14. RMS-I   
      IRMS_UNITS = 'Jy/Bm'
      lenIRMS_UNITS = nchar(IRMS_UNITS)

      err = dbmkoptlist(1,IrmsOpts)
      err = dbaddcopt(IrmsOpts,DBOPT_UNITS,IRMS_UNITS,lenIRMS_UNITS)

      err = dbputpv1 (dbfile, "RMS_I",5,mesh_name,len_mesh_name, 
     -                Irms_arr, nused, DB_FLOAT, IrmsOpts, ierr)

      !---------------------------------------------
      ! 15. RMS-Q   
      QRMS_UNITS = 'Jy/Bm'
      lenQRMS_UNITS = nchar(QRMS_UNITS)

      err = dbmkoptlist(1,QrmsOpts)
      err = dbaddcopt(QrmsOpts,DBOPT_UNITS,QRMS_UNITS,lenQRMS_UNITS)

      err = dbputpv1 (dbfile, "RMS_Q",5,mesh_name,len_mesh_name, 
     -                Qrms_arr, nused, DB_FLOAT, QrmsOpts, ierr)

      !---------------------------------------------
      ! 16. RMS-U   
      URMS_UNITS = 'Jy/Bm'
      lenURMS_UNITS = nchar(URMS_UNITS)

      err = dbmkoptlist(1,UrmsOpts)
      err = dbaddcopt(UrmsOpts,DBOPT_UNITS,URMS_UNITS,lenURMS_UNITS)

      err = dbputpv1 (dbfile, "RMS_U",5,mesh_name,len_mesh_name, 
     -                Urms_arr, nused, DB_FLOAT, UrmsOpts, ierr)

      !---------------------------------------------
      ! 17. RMS-V   
      VRMS_UNITS = 'Jy/Bm'
      lenVRMS_UNITS = nchar(VRMS_UNITS)

      err = dbmkoptlist(1,VrmsOpts)
      err = dbaddcopt(VrmsOpts,DBOPT_UNITS,VRMS_UNITS,lenVRMS_UNITS)

      err = dbputpv1 (dbfile, "RMS_V",5,mesh_name,len_mesh_name, 
     -                Vrms_arr, nused, DB_FLOAT, VrmsOpts, ierr)
      !---------------------------------------------
      ierr = DBClose(dbfile)

      !------------------------------------
      ! TODO: 
      ! Make provision to BLANK outliers: 
      ! 
      !------------------------------------
      ! Write to a TEXT file some data for quick plot:
      templine = outfile(1:index(outfile,'.SILO')-1)//'.TXT'
      
      open(21,file=templine,status='unknown')
      write(21,*)"# RA     Dec         I         LP        dLP     PA 
     -      XRAY       V          Irms        Qrms       Urms      
     -Vrms"
      do k = 1,nused
         write(21,fmt=201)RA_arr(k), Dec_arr(k), I_arr(k), 
     -              LP_arr(k), dLP_arr(k), PA_arr(k),X_arr(k),
     -              V_arr(k),Irms_arr(k), Qrms_arr(k), Urms_arr(k), 
     -              Vrms_arr(k)   
      enddo
201   format(F6.1,2x,F6.1,1x,F11.4,1x,F11.6,2x,F6.4,1x,F7.4,1x,
     -       F11.6,F11.6,1x,F11.6,1x,F11.6,1x,F11.6,1x,F11.6)

      ! ============================================================
      ! TEST: A temporary file: 
      !------------------------------------
      ! Write to a TEXT file some data for quick plot:
      templine = outfile(1:index(outfile,'.SILO')-1)//'.2.TXT'
      
      open(22,file=templine,status='unknown')
      write(22,*)"# RA     Dec         I         Q        U       V 
     -      LP       dLP          PA        XRay       LP & PA (mean rem
     -oved)"
      do k = 1,nused
!         write(22,fmt=202)RA_arr(k), Dec_arr(k), I_arr(k), 
!     -              Q_arr(k), U_arr(k), V_arr(k),LP_arr(k),
!     -              dLP_arr(k),PA_arr(k), X_arr(k)   
      ! Remove mean from the linpol vector: 
      

      qtmp = Q_rm_arr(k) - qmean 
      utmp = U_rm_arr(k) - umean 

      linpol = sqrt(qtmp**2 + utmp**2)  
      polang = atan2(utmp,qtmp) 

         write(22,fmt=202)RA_arr(k), Dec_arr(k), I_arr(k), 
     -              Q_rm_arr(k), U_rm_arr(k), V_arr(k),LP_arr(k),
     -              dLP_arr(k),PA_arr(k), X_arr(k), linpol, polang  
      enddo
202   format(F6.1,2x,F6.1,1x,F11.4,1x,F11.6,2x,F6.4,1x,F7.4,1x,
     -       F11.6,1x,F11.6,1x,F11.6,1x,F11.6,F11.6,1x,F11.6,1x)
   
      ! ============================================================

      close(21) 
      close(22) 
      write(*,*)" "
      write(*,*)" OUTFILES WRITTEN IN: ",out_path(1:nchar(out_path))

      end


      include '/usr/lib/subroutine_lib/fort_lib.f'
      include '/usr/lib/subroutine_lib/nchar.f'
      include '/usr/lib/subroutine_lib/printerror.f'
      !include '/usr/lib/subroutine_lib/load_fits_image.f'
      include 'load_fits_image.f'
