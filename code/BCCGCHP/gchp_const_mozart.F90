!xiaolu
!Purpose: set constant number for MOZART
!XIAOLU,2017/06/01
!From mo_const_mozart.F
      module gchp_const_mozart

      use pmgrid,         only:  plon, plat   !add by zf 2008.03.11

      implicit none

      save

      real, parameter ::  gravit = 9.80616    ! m/s
      real, parameter ::  rgrav  = 1./gravit
      real, parameter ::  dayspy = 365.       ! days per year
      real, parameter ::  rearth = 6.37122e6  ! radius earth (m)

      real ::  pi                 ! radians
      real ::  twopi              ! 2*pi (radians)
      real ::  pid2               ! pi/2 (radians)
      real ::  r2d                ! radians to degrees
      real ::  d2r                ! degrees to radians 
      real ::  lat25 = 0.         ! 25 latitude (radians)
      real ::  lat40 = 0.         ! 40 latitude (radians)
      real ::  lat45 = 0.         ! 45 latitude (radians)
      real ::  lat59 = 0.         ! 59 latitude (radians)
      real ::  lat60 = 0.         ! 60 latitude (radians)
      real ::  lat70 = 0.         ! 70 latitude (radians)
      real ::  phi(plat)          ! latitudes (radians)
      real ::  lam(plon)          ! longitudes( radians )
      real ::  sinlam(plon)       ! sine of longitudes
      real ::  coslam(plon)       ! cose of longitudes

      integer ::  ktop       !add by zf 2008.07.11
      CONTAINS

      subroutine GC_CONSTANTS_INTI()

      implicit none

     pi     = 4.*ATAN( 1. )
!      pi     = z'400921fd54442d18'
      twopi  = 2.*pi
      pid2   = .5*pi
      d2r    = pi/180.
      r2d    = 180./pi

      lat25 = d2r * 25.
      lat40 = d2r * 40.
      lat45 = d2r * 45.
      lat59 = d2r * 59.
      lat60 = d2r * 60.
      lat70 = d2r * 70.

      end subroutine GC_CONSTANTS_INTI

      end module gchp_const_mozart
