
      subroutine my_randn(npts,mu,sigma,OutArr)
      implicit none


      integer*4 npts
      integer*4 seed
      integer*4 iseed
      integer*4 ii
      real*4 gasdev
      real*4 OutArr(*)
      real*4 mu, sigma


      call get_seed_for_rand(seed)
      iseed = -seed 
      !write(*,*)"Current seed: ",iseed
      do ii = 1,npts
         OutArr(ii) = mu + sigma*gasdev(iseed)
      enddo

      return


      end
      include 'gasdev.f'
      include 'ran1.for'
      include 'ran2.for'

