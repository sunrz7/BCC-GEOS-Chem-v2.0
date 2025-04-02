!xiaolu
!---------------------------------------
!from MOZART2 mo_charutl.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/04/28
!---------------------------------------

      module GCHP_CHARUTL

      private :: TOKLEN, LASTND

      CONTAINS

      integer function TOKLEN( cs )
!-----------------------------------------------------------------------
! 	... Return token length, i.e., starting from begining of string
!           return length of first token as delimited by either a blank
!           or a null byte.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: cs       !  Input character string

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: m, n

      m = 0
      do n = 1,LEN(cs)
         if( cs(n:n) == ' ' .or. cs(n:n) == char(0) ) then
            exit
         else
            m = m + 1
         end if
      end do
      TOKLEN = m

      end function TOKLEN

      integer function GLC( cs )
!-----------------------------------------------------------------------
! 	... Position of last significant character in string. 
!           Here significant means non-blank or non-null.
!           Return values:
!               > 0  => position of last significant character
!               = 0  => no significant characters in string
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: cs       !  Input character string

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: l, n

      l = LEN( cs )
      if( l == 0 ) then
         GLC = 0
         return
      end if

      do n = l,1,-1
         if( cs(n:n) /= ' ' .and. cs(n:n) /= CHAR(0) ) then
            exit
         end if
      end do
      GLC = n

      end function GLC

      integer function LASTND( cs )
!-----------------------------------------------------------------------
! 	... Position of last non-digit in the first input token.
! 	    Return values:
!     	    > 0  => position of last non-digit
!     	    = 0  => token is all digits (or empty)
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: cs       !  Input character string

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: n, nn, digit

      n = GLC( cs )
      if( n == 0 ) then     ! empty string
         LASTND = 0
         return
      end if

      do nn = n,1,-1
         digit = ICHAR( cs(nn:nn) ) - ICHAR('0')
         if( digit < 0 .or. digit > 9 ) then
            LASTND = nn
            return
         end if
      end do

      LASTND = 0    ! all characters are digits

      end function LASTND

      integer function LASTSL( cs )
!-----------------------------------------------------------------------
! 	... Position of last slash (/) in the first input token.
!           Return values:
!               > 0  => position of last slash
!               = 0  => slash not found in token
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: cs       !  Input character string

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: l, n

      l = TOKLEN( cs )
      if( l == 0 ) then
         LASTSL = 0
         return
      end if

      do n = l,1,-1 
         if( cs(n:n) == '/' ) then
            exit
         end if
      end do
      LASTSL = n

      end function LASTSL

      function I2CHAR( kint, kilen )
!-----------------------------------------------------------------------
! 	... Convert integer kint to character type, left-justify and return number
!           of characters in i2char in kilen
!           Return value:
!           The ascii representation of kint
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        kint                                  ! integer to be converted
      integer, intent(out) :: &
        kilen                                 ! length of the ascii representation of kint

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: &
        i
      character(len=20) :: &
        chtem

!-----------------------------------------------------------------------
! 	... Function declarations
!-----------------------------------------------------------------------
      character(len=20) :: I2CHAR

      write(chtem,'(i20)') kint
      i = 20
      do
         if( chtem(i:i) /= ' ' ) then
            i = i - 1
            if( i < 1 ) then
               I2CHAR = chtem
               exit
            end if
            cycle
         else
            I2CHAR = chtem(i+1:20)
	    exit
         end if
      end do
      kilen = GLC( I2CHAR )

      end function I2CHAR

      subroutine UPCASE( lstring, outstring )
!-----------------------------------------------------------------------
! 	... Convert character string lstring to upper case
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy argument
!-----------------------------------------------------------------------
      character(len=*), intent(in)  :: lstring
      character(len=*), intent(out) :: outstring

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: i

      outstring = ' '
      do i = 1,GLC( lstring )
         outstring(i:i) = lstring(i:i)
         if( ICHAR( lstring(i:i) ) >= 97 .and. ICHAR( lstring(i:i) ) <= 122 ) then
             outstring(i:i) = CHAR(ICHAR(lstring(i:i)) - 32)
         end if
      end do

      end subroutine UPCASE

      integer function INCSTR( s, inc )
!-----------------------------------------------------------------------
! 	... Increment a string whose ending characters are digits.
!           The incremented integer must be in the range [0 - (10**n)-1]
!           where n is the number of trailing digits.
!           Return values:
!
!            0 success
!           -1 error: no trailing digits in string
!           -2 error: incremented integer is out of range
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy variables
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        inc                                       ! value to increment string (may be negative)
      character(len=*), intent(inout) :: &
        s                                         ! string with trailing digits


!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: &
        i, &                          ! index
        lstr, &                       ! number of significant characters in string
        lnd, &                        ! position of last non-digit
        ndigit, &                     ! number of trailing digits
        ival, &                       ! integer value of trailing digits
        pow, &                        ! power of ten
        digit                         ! integer value of a single digit

      lstr   = GLC( s )
      lnd    = LASTND( s )
      ndigit = lstr - lnd

      if( ndigit == 0 ) then
         INCSTR = -1
         return
      end if

!-----------------------------------------------------------------------
!     	... Calculate integer corresponding to trailing digits.
!-----------------------------------------------------------------------
      ival = 0
      pow  = 0
      do i = lstr,lnd+1,-1
         digit = ICHAR(s(i:i)) - ICHAR('0')
         ival  = ival + digit * 10**pow
         pow   = pow + 1
      end do

!-----------------------------------------------------------------------
!     	... Increment the integer
!-----------------------------------------------------------------------
      ival = ival + inc
      if( ival < 0 .or. ival > 10**ndigit-1 ) then
         INCSTR = -2
         return
      end if

!-----------------------------------------------------------------------
!     	... Overwrite trailing digits
!-----------------------------------------------------------------------
      pow = ndigit
      do i = lnd+1,lstr
         digit  = MOD( ival,10**pow ) / 10**(pow-1)
         s(i:i) = CHAR( ICHAR('0') + digit )
         pow    = pow - 1
      end do

      INCSTR = 0

      end function INCSTR

      end module GCHP_CHARUTL
