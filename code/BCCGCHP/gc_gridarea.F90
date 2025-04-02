#include <misc.h>
#include <params.h>

module gc_gridarea
!----------------------------------------------------------------------- 
! 
! Purpose: Module to calculate the area of the atmosphere grid
!          for use by GEOS-Chem 
! 
! Author: M. Long
!
! Adapted from ccsm_msg.F90.
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8             ! atmospheric model precision
  use pmgrid,       only: plat, plon, beglat, endlat, plond, masterproc, iam      ! Model grid
  use ppgrid,       only: pcols, pver, begchunk, endchunk  ! Physics grid
  use phys_grid,    only: scatter_field_to_chunk
  use shr_kind_mod, only: SHR_KIND_IN                    ! defines CCSM real & integer kinds

!------------------wangln 20050609  apply CPL5--------------------------
 
  use abortutils, only: endrun
  use rgrid, only: nlon                                 ! Reduced grid
!
  implicit none

!--------------------------------------------------------------------------
! Public interface and information
!--------------------------------------------------------------------------
  private         ! Make the default access private, explicitly declare public
  public :: gc_calcarea

  contains

!===============================================================================

  subroutine gc_calcarea(area_m2_chunk,lchnk)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
!-----------------------------------------------------------------------

    use infnan
    use commap, only: latdeg, londeg
    use dycore, only: dycore_is
    use time_manager, only: get_nstep, get_step_size
!    use CMN_GCTM_MOD, only: Re
    
    integer, optional, intent(in) :: lchnk
    real(r8), intent(out) :: area_m2_chunk (pcols)  ! Area in m^2  for each local chunk

!--------------------------Local Variables------------------------------
    integer lat, lon, i, j, n     ! loop indices
    integer nstep                 ! current time step
    integer msgpday               ! number of send/recv msgs per day
    integer sizebuf               ! size of buffer for sending grid data to coupler
    integer startpoint            ! starting value for grid numbering scheme
    integer(SHR_KIND_IN) ::  mask(plon,plat)       ! Mask of valid data
    real(r8),allocatable :: sbuf(:,:)  ! array for holding grid data to be sent to coupler
    real(r8) dtime                ! timestep size [s]
    real(r8) area(plon,plat)      ! Area in radians squared for each grid point
    real(r8) area_m2(plon,plat)   ! Area in m^2 for each gridbox
    real(r8) area_cm2(plon,plat)  ! Area in cm^2 for each gridbox
    real(r8) :: tmp(pcols,begchunk:endchunk)
    real(r8) clondeg(plon,plat)   ! Longitude grid
    real(r8) clatdeg(plon,plat)   ! latitude grid as 2 dimensional array
    real(r8) ns_vert(4,plon,plat) ! latitude grid vertices
    real(r8) ew_vert(4,plon,plat) ! longitude grid vertices
    real(r8) del_theta            ! difference in latitude at a grid point
    real(r8) del_phi              ! difference in longitude at a grid point
    real(r8) pi                   ! mathmatical constant 3.1415...
    real(r8) degtorad             ! convert degrees to radians
    real(r8) :: spval = 1.e30
    REAL*8, PARAMETER :: Re     =   6.375d6
!-----------------------------------------------------------------------

       if (masterproc) then

       nstep = get_nstep()
       dtime = get_step_size()

! Constants
!
       pi       = acos(-1.)
       degtorad  = pi / 180.0
!
! Mask for which cells are active and inactive and 2D latitude grid
!
       mask(:,:)    = 0        ! Initialize mask so that cells are inactive
       clatdeg(:,:) = spval
       clondeg(:,:) = spval
       do lat = 1, plat
         mask(1:nlon(lat),lat)    = 1     ! Active cells
         clatdeg(1:nlon(lat),lat) = latdeg(lat) ! Put latitude in 2D array
         clondeg(1:nlon(lat),lat) = londeg(1:nlon(lat),lat)
       end do
!
! Send vertices of each grid point
! Verticies are ordered as follows: 
! 1=lower left, 2 = upper left, 3 = upper right, 4 = lower right
!
       ns_vert(:,:,:) = spval
       ew_vert(:,:,:) = spval
!
! Longitude vertices
!
       do lat = 1, plat
         ew_vert(1,1,lat)             = (londeg(1,lat) - 360.0 + londeg(nlon(lat),lat))*0.5
         ew_vert(1,2:nlon(lat),lat)   = (londeg(1:nlon(lat)-1,lat) + &
                                         londeg(2:nlon(lat),lat))*0.5
         ew_vert(2,:nlon(lat),lat)    = ew_vert(1,:nlon(lat),lat)  ! Copy lowleft corner to upleft
         ew_vert(3,:nlon(lat)-1,lat)  = ew_vert(1,2:nlon(lat),lat)
         ew_vert(3,nlon(lat),lat)     = (londeg(nlon(lat),lat) + (360.0 + londeg(1,lat)))*0.5
         ew_vert(4,:nlon(lat),lat)    = ew_vert(3,:nlon(lat),lat)  ! Copy lowright corner to upright
       end do
!
! Latitude
!
       if ( dycore_is('LR') )then
         ns_vert(1,:nlon(1),1)         = -90.0 + (latdeg(1) - latdeg(2))*0.5
         ns_vert(2,:nlon(plat),plat)   =  90.0 + (latdeg(plat) - latdeg(plat-1))*0.5
       else
         ns_vert(1,:nlon(1),1)         = -90.0
         ns_vert(2,:nlon(plat),plat)   =  90.0
       end if
       ns_vert(4,:nlon(1),1)       = ns_vert(1,nlon(1),1)        ! Copy lower left to lower right
       ns_vert(3,:nlon(plat),plat) = ns_vert(2,nlon(plat),plat)  ! Copy up left to up right
       do lat = 2, plat
         ns_vert(1,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat-1) )*0.5
         ns_vert(4,:nlon(lat),lat) = ns_vert(1,:nlon(lat),lat)
       end do
       do lat = 1, plat-1
         ns_vert(2,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat+1) )*0.5
         ns_vert(3,:nlon(lat),lat) = ns_vert(2,:nlon(lat),lat)
       end do
!
! Get area of grid cells (as radians squared)
!
       area(:,:) = 0.0
       do lat = 1, plat
         do lon = 1, nlon(lat)
           del_phi = sin( ns_vert(2,lon,lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat) = del_theta*del_phi
         end do
       end do
!
! If grid has a pole point (as in Lin-Rood dynamics
!
      if ( dycore_is('LR') )then
         lat = 1
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi = -sin( latdeg(lat)*degtorad ) + sin( ns_vert(2,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
         lat = plat
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi =  sin( latdeg(lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
       end if
       if ( abs(sum(area) - 4.0*pi) > 1.e-12 )then
         write (6,*) ' GC_GRIDAREA: sum of areas on globe does not = 4*pi'
         write (6,*) ' sum of areas = ', sum(area)
         call endrun
       end if

! NOTE:  Numbering scheme is: West to East and South to North
! starting at south pole.  Should be the same as what's used
! in SCRIP


       area_m2  = Re*Re*area(:,:)   ! Radians^2 to  m^2
       area_cm2 = area_m2*1.d4      ! Radians^2 to cm^2

    endif  ! end of if-masterproc

!   write(*,*),'check area_m2 at gc_gridarea.F90',area_m2


    call scatter_field_to_chunk(1,1,1,plon,area_m2,tmp(1,begchunk) )
! write(*,*),'xiaolu check size of tmp geoschem/BCCGCHP/gc_gridarea.F90',pcols,endchunk-begchunk+1

! write(*,*),'check area_m2 at gc_gridarea.F90',area_m2
    if (masterproc) then
!       do i=1,pcols
!          write(*,'(a,i3,4e11.1)'),'<>tmp(pcols): ',i,tmp(i,begchunk:endchunk)
!       enddo
!       do i=1,plat
!          write(*,'(a,i3,e11.1)'),'<>area(plat): ',i,area_m2(1,i)
!       enddo
    endif
! write(*,*),'xiaolu check tmp geoschem/BCCGCHP/gc_gridarea.F90',shape(tmp)
 

    if (present(lchnk)) then
       area_m2_chunk  = tmp(:,lchnk)
    else
       area_m2_chunk  = tmp(:,1)
    endif

    return
  end subroutine gc_calcarea

!===============================================================================

end module gc_gridarea
