#include <misc.h>
#include <params.h>

      module gc_grid
!---------------------------------------------------------------------
! ... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none
      save
      integer, parameter :: &
                pcnst = 220+1, &  ! 63 +1, & ! number of advected constituents including cloud water
                pcnstm1 = 220     ! 63, & ! number of advected constituents excluding cloud water
      integer, parameter :: & 
                plev = PLEV ,   & ! number of vertical levels
                plevp = plev+1, & ! plev plus 1
                plevm = plev-1, & ! plev minus 1
                plon = PLON,    & ! number of longitudes
                plat = PLAT       ! number of latitudes
      integer, parameter :: &
                pnats = 0 ! number of non-advected trace species
      integer :: nodes ! mpi task count

!     integer :: plonl ! longitude tile dimension

      integer :: pplon ! longitude tile count
      integer :: plnplv ! plonl * plev
      end module gc_grid
