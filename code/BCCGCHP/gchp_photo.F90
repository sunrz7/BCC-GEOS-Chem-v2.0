!xiaolu
!XIAOLU,2017/06/01
!From mo_photo.F, but abandon those not in need
!
      module gchp_photo
!----------------------------------------------------------------------
!       ... photolysis interp table and related arrays
!----------------------------------------------------------------------

      implicit none

      private
!      public :: prate_inti, photo, set_ub_col, setcol, diurnal_geom, sundis
      public ::diurnal_geom

!      character(len=256) , public  ::  bndprate

      save


      contains


      subroutine diurnal_geom( longitude, latitude, time_of_year, polar_night, polar_day, &
                               sunon, sunoff, loc_angle, zen_angle )
!------------------------------------------------------------------
!       ... diurnal geometry factors
!------------------------------------------------------------------

      use gchp_const_mozart, only : pi, twopi, pid2, dayspy, d2r

      implicit none

!------------------------------------------------------------------
!       ... dummy arguments
!------------------------------------------------------------------
      real,    intent(in)  ::     latitude           ! latitude 
      real,    intent(in)  ::     longitude          ! latitude 
      real, intent(in)     ::     time_of_year       ! time of year
      real, intent(out)    ::     sunon           ! sunrise angle in radians
      real, intent(out)    ::     sunoff          ! sunset angle in radians
      real, intent(out)    ::     zen_angle       ! solar zenith angle
      real, intent(out)    ::     loc_angle       ! "local" time angle
      logical, intent(out) ::     polar_day       ! continuous daylight flag
      logical, intent(out) ::     polar_night     ! continuous night flag

!------------------------------------------------------------------
!        ... local variables
!------------------------------------------------------------------
      real    ::  dec_max
      real    ::  declination
      real    ::  doy_loc            ! day of year
      real    ::  tod                ! time of day
      real    ::  sin_dec, cos_dec   ! sin, cos declination
      real    ::  cosphi             ! cos latitude
      real    ::  sinphi             ! sin latitude

      dec_max     = 23.45 * d2r
      sinphi      = sin( latitude )
      cosphi      = cos( latitude )
      polar_day   = .false.
      polar_night = .false.
!------------------------------------------------------------------
!        note: this formula assumes a 365 day year !
!------------------------------------------------------------------
      doy_loc     = aint( time_of_year )
      declination = dec_max * cos((doy_loc - 172.)*twopi/dayspy)
!------------------------------------------------------------------
!        determine if in polar day or night
!        if not in polar day or night then
!        calculate terminator longitudes
!------------------------------------------------------------------
      if( abs(latitude) >= (pid2 - abs(declination)) ) then
         if( sign(1.,declination) == sign(1.,latitude) ) then
            polar_day = .true.
            sunoff    = 2.*twopi
            sunon     = -twopi
         else
            polar_night  = .true.
            zen_angle = -1.0
            return
         end if
      else
         sunoff = acos( -tan(declination)*tan(latitude) )
         sunon  = twopi - sunoff
      end if

      sin_dec = sin( declination )
      cos_dec = cos( declination )
!------------------------------------------------------------------
!       ... compute base for zenith angle
!------------------------------------------------------------------
      tod = (time_of_year - doy_loc) + .5
!-------------------------------------------------------------------
!        note: longitude 0 (greenwich) at 0:00 hrs
!              maps to local angle = pi
!-------------------------------------------------------------------
      loc_angle = (tod + longitude/2.8125*180./pi/128.)*twopi
      loc_angle = mod( loc_angle,twopi )

      if( polar_day ) then
         zen_angle = acos( sinphi*sin_dec + cosphi*cos_dec*cos(loc_angle) )
      else
         if( loc_angle <= sunoff .or. loc_angle >= sunon ) then
            zen_angle = acos( sinphi*sin_dec + cosphi*cos_dec*cos(loc_angle) )
         else
            zen_angle = -1.
         endif
      end if

      end subroutine diurnal_geom

      end module gchp_photo
