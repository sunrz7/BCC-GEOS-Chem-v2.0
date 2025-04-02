!xiaolu
!---------------------------------------
!from MOZART2 mo_calendar.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/04/28
!---------------------------------------

      module GCHP_CALENDAR

      implicit none

      save

      character(len=16) :: &
        type = '365             '  ! calendar type.  Currently '365' or 'gregorian'

      CONTAINS

      real function DIFFDAT( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
! 	... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Input arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2

      if( type(1:3)  ==  '365' ) then
         DIFFDAT = DIFFDAT365( dat1, sec1, dat2, sec2 )
      else if ( type(1:9)  ==  'gregorian' ) then
         DIFFDAT = DIFFDATGRG( dat1, sec1, dat2, sec2 )
      end if

      end function DIFFDAT

!-------------------------------------------------------------------------

      real function DIFFDAT365( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
!       ... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           N.B. Assume 1 year = 365 days.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2


!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: lastd, lasts, firstd, firsts, &
                 yr2, mo2, dy2, doy2, yr1, mo1, dy1, doy1, ndays
      real :: sign, days

!-----------------------------------------------------------------------
!       ... Dates equal
!-----------------------------------------------------------------------
      if( dat1 == dat2 .and. sec1 == sec2 ) then
         DIFFDAT365 = 0.
         return
      end if

!-----------------------------------------------------------------------
!       ... Which date is later?
!-----------------------------------------------------------------------
      if( dat2 > dat1 ) then
         sign = 1.
         lastd  = dat2
         lasts  = sec2
         firstd = dat1
         firsts = sec1
      else if( dat2 < dat1 ) then
         sign   = -1.
         lastd  = dat1
         lasts  = sec1
         firstd = dat2
         firsts = sec2
      else
         if( sec2 > sec1 ) then
            sign   = 1.
            lastd  = dat2
            lasts  = sec2
            firstd = dat1
            firsts = sec1
         else
            sign   = -1.
            lastd  = dat1
            lasts  = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

!-----------------------------------------------------------------------
!       ... Compute number of days between lastd and firstd
!-----------------------------------------------------------------------
      yr2  = lastd / 10000
      mo2  = MOD( lastd, 10000 ) / 100
      dy2  = MOD( lastd, 100 )
      doy2 = DOY( mo2, dy2 )

      yr1  = firstd / 10000
      mo1  = MOD( firstd, 10000 ) / 100
      dy1  = MOD( firstd, 100 )
      doy1 = DOY( mo1, dy1 )

      ndays = 365*(yr2 - yr1) + doy2 - doy1

!-----------------------------------------------------------------------
!       ... Adjust for remaining seconds
!-----------------------------------------------------------------------
      days = REAL( ndays ) + REAL( lasts - firsts )/86400.

!-----------------------------------------------------------------------
!       ... Adjust sign
!-----------------------------------------------------------------------
      DIFFDAT365 = sign * days

      end function DIFFDAT365

!-------------------------------------------------------------------------

      real function DIFFDATGRG( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
!       ... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           N.B. Assume Gregorian calendar.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2


!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: lastd, lasts, firstd, firsts, ndays
      real    :: sign, days

!-----------------------------------------------------------------------
!       ... Dates equal?
!-----------------------------------------------------------------------
      if( dat1 == dat2 .and. sec1 == sec2 ) then
         DIFFDATGRG = 0.
         return
      end if

!-----------------------------------------------------------------------
!       ... Which date is later?
!-----------------------------------------------------------------------
      if( dat2 > dat1 ) then
         sign   = 1.
         lastd  = dat2
         lasts  = sec2
         firstd = dat1
         firsts = sec1
      else if( dat2 < dat1 ) then
         sign   = -1.
         lastd  = dat1
         lasts  = sec1
         firstd = dat2
         firsts = sec2
      else
         if( sec2 > sec1 ) then
            sign   = 1.
            lastd  = dat2
            lasts  = sec2
            firstd = dat1
            firsts = sec1
         else
            sign   = -1.
            lastd  = dat1
            lasts  = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

!-----------------------------------------------------------------------
!       ... Compute number of days between lastd and firstd
!-----------------------------------------------------------------------
      ndays = GREG2JDAY( lastd ) - GREG2JDAY( firstd )

!-----------------------------------------------------------------------
!       ... Adjust for remaining seconds
!-----------------------------------------------------------------------
      days = REAL( ndays ) + REAL( lasts - firsts )/86400.

!-----------------------------------------------------------------------
!       ... Adjust sign
!-----------------------------------------------------------------------
      DIFFDATGRG = sign * days

      end function DIFFDATGRG

!-----------------------------------------------------------------------

      integer function DOY( mon, cday )
!-----------------------------------------------------------------------
!       ... Compute day of year ignoring leap years.
!           Returns values in the range [1,365].
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: mon, cday

      integer, save :: jdcon(12) = &
            (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

      doy = jdcon(mon) + cday

      end function DOY

!-----------------------------------------------------------------------

      integer function GREG2JDAY( date )
!-----------------------------------------------------------------------
!       ... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: date

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: yy, mm, dd
      integer :: ap, mp
      integer :: y, d, n, g

!-----------------------------------------------------------------------
!       ... Extract year, month, and day from date
!-----------------------------------------------------------------------
      yy = date / 10000
      mm = MOD( ABS(date),10000 ) / 100
      dd = MOD( ABS(date),100 )

!-----------------------------------------------------------------------
!       ... Modify year and month numbers
!-----------------------------------------------------------------------
      ap = yy - (12 - mm)/10
      mp = MOD( mm-3,12 )
      if( mp < 0 ) then
         mp = mp + 12
      end if

!-----------------------------------------------------------------------
!       ... Julian day
!-----------------------------------------------------------------------
      y = INT( 365.25*( ap + 4712 ) )
      d = INT( 30.6*mp + .5 )
      n = y + d + dd  + 59
      g = INT( .75*INT( ap/100 + 49 ) ) - 38
      GREG2JDAY = n - g

      end function GREG2JDAY

      real function CALDAYR( idate, isec )
!-----------------------------------------------------------------------
!       ... Calendar day with fractional part.  Returns values in the
!           range [1., 366.)
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: idate, isec

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: mon, day

      mon = MOD( idate,10000 ) / 100
      day = MOD( idate,100 )

      CALDAYR = DOY( mon, day ) + REAL( isec )/86400.

      end function CALDAYR

      end module GCHP_CALENDAR
