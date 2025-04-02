#include <misc.h>
#include <params.h>

      module gchp_surf
      use mpishorthand      , only: mpicom, mpiint, mpir8

!--------------------------------------------------------------------------
!	... surface data
!--------------------------------------------------------------------------

      implicit none

      character(len=256), public  ::  bndsurf  ! full pathname for Mozart surface variables dataset 

      save

      real, allocatable, dimension(:,:)   :: desertbdy 
      real, allocatable, dimension(:,:)   :: oceanbdy
      real, allocatable, dimension(:,:,:) :: seaicebdy
      real, allocatable, dimension(:,:)   :: ocean
      real, allocatable, dimension(:,:,:) :: seaice

      contains

      subroutine surf_inti( plonl, platl )
!-----------------------------------------------------------------------
! 	... read surface data database
!-----------------------------------------------------------------------

      use gchp_netcdf
      use gchp_const_mozart,  only : d2r
!add by zf 2008.06.02
      use ppgrid,           only: pcols, begchunk, endchunk
      use pmgrid,           only: masterproc,plond
      use ioFileMod,        only: getfil
      use abortutils,       only: endrun
      use error_messages,   only: handle_ncerr
      use phys_grid,        only:  scatter_field_to_chunk

!end zf

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl
      integer, intent(in) :: platl

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      real, allocatable, dimension(:,:) :: desert 

      integer :: ncid
      integer :: nlon
      integer :: nlat
      integer :: nmonth
      integer :: ierr
      integer :: m
      integer :: vid, dimid_lon, dimid_lat, dimid_month
      integer, dimension(3) :: start, count

      real, allocatable :: lat(:)
      real, allocatable :: lon(:)

      character(len=256) :: locfn    ! netcdf local filename to open   add by zf

!-----------------------------------------------------------------------
!       ... open netcdf file
!-----------------------------------------------------------------------

      allocate( desertbdy(pcols,begchunk:endchunk), stat=ierr)
      allocate( oceanbdy(pcols,begchunk:endchunk), stat=ierr)
      allocate( seaicebdy(pcols,begchunk:endchunk,12), stat=ierr)

!add by zf 2008.05.28
    if (masterproc) then
      call getfil (bndsurf, locfn)
      call wrap_open(locfn,0,ncid)
!end zf
!-----------------------------------------------------------------------
!       ... check number of months
!-----------------------------------------------------------------------
      call handle_ncerr( nf_inq_dimid( ncid, 'month', dimid_month ), &
                         'surf_inti: failed to find dimension month' )
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_month, nmonth ), &
                         'surf_inti: failed to get length of dimension month' )
      if( nmonth /= 12 ) then
         write(*,*) 'surf_inti: error! nmonth = ',nmonth,', expecting 12'
         call endrun
      end if
!-----------------------------------------------------------------------
!       ... get latitude and longitude
!-----------------------------------------------------------------------
      call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid_lat ), &
                         'surf_inti: failed to find dimension lat' )
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_lat, nlat ), &
                         'surf_inti: failed to get length of dimension lat' )
      allocate( lat(nlat), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'surf_inti: lat allocation error = ',ierr
         call endrun
      end if
      call handle_ncerr( nf_inq_varid( ncid, 'lat', vid ), &
                         'surf_inti: failed to find variable lat' )
      call handle_ncerr( nf_get_var_double( ncid, vid, lat ), &
                         'surf_inti: failed to read variable lat' )
      lat(:nlat) = lat(:nlat) * d2r

      call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid_lon ), &
                         'surf_inti: failed to find dimension lon' )
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_lon, nlon ), &
                         'surf_inti: failed to get length of dimension lon' )
      endif   ! master
!----------
#if (defined SPMD)
      call mpibcast( nlon, 1, mpiint, 0, mpicom )
      call mpibcast( nlat, 1, mpiint, 0, mpicom )
#endif
      allocate( desert(nlon,nlat), &
                ocean (nlon,nlat), &
                seaice(nlon,nlat,12), stat=ierr )
      desert(:,:)   = 0.0
      ocean(:,:)    = 0.0
      seaice(:,:,:) = 0.0

!-------------     
    if (masterproc) then
      allocate( lon(nlon), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'surf_inti: lon allocation error = ',ierr
         call endrun
      end if
      call handle_ncerr( nf_inq_varid( ncid, 'lon', vid ), &
                         'surf_inti: failed to find variable lon' )
      call handle_ncerr( nf_get_var_double( ncid, vid, lon ), &
                         'surf_inti: failed to read variable lon' )
      lon(:nlon) = lon(:nlon) * d2r

!-----------------------------------------------------------------------
! 	... read the surface data
!-----------------------------------------------------------------------
      start(:) = (/ 1, 1, 1 /)
      count(:) = (/ nlon, nlat, nmonth /)

      call handle_ncerr( nf_inq_varid( ncid, 'desert', vid ), 'surf_inti: getting desert id' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start(1:2), count(1:2), desert ), &
                         'surf_inti: getting desert' )

      call handle_ncerr( nf_inq_varid( ncid, 'ocean', vid ), 'surf_inti: getting ocean id' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start(1:2), count(1:2), ocean ), &
                         'surf_inti: getting ocean' )

      call handle_ncerr( nf_inq_varid( ncid, 'seaice', vid ), 'surf_inti: getting seaice id' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, seaice ), &
                         'surf_inti: getting seaice' )

      call handle_ncerr( nf_close(ncid), 'surf_inti: closing netcdf file' )

!-----------------------------------------------------------------------
!	... due to possible input data errors limit all
!           variables to a the range [0,1]
!lwh 10/00 -- not really necessary when using the netcdf input file.
!             the data in this file are already limited to [0,1].
!-----------------------------------------------------------------------
      desert(:,:)   = min( 1., max( 0., desert(:,:) ) )
      ocean(:,:)    = min( 1., max( 0., ocean(:,:) ) )
      seaice(:,:,:) = min( 1., max( 0., seaice(:,:,:) ) )

      call scatter_field_to_chunk(1,1,1,nlon,desert,desertbdy)
      call scatter_field_to_chunk(1,1,1,nlon,ocean,oceanbdy)
      call scatter_field_to_chunk(1,1,12,nlon,seaice,seaicebdy)
#if (defined SPMD)
      else
      call scatter_field_to_chunk(1,1,1,nlon,desert,desertbdy)
      call scatter_field_to_chunk(1,1,1,nlon,ocean,oceanbdy)
      call scatter_field_to_chunk(1,1,12,nlon,seaice,seaicebdy)
#endif
      end if !(masterproc)

      DEALLOCATE(desert)

      end subroutine surf_inti

      end module gchp_surf
