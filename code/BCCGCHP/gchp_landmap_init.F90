#include <misc.h>
      MODULE GCHP_Landmap_Init
!---------------------------------------------------------------------
!       Initiate BCC-AVIM for GCHP dry deposition,read landtype fraction
!       from input file
!History:
! XIAOLU,2017/08/20: initial version; mimic from dvel_inti
!---------------------------------------------------------------------
      implicit none

      character(len=256), public  ::  file_land     !full pathname for BCC-AVIM landtype file

      save

      private
      public  :: drydep_landmap_init
     ! public  :: GC_INITRUN_FORDRYDEP
      real, allocatable,public::ltfrac (:,:,:)

        CONTAINS
     
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !TEST ONLY
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !  SUBROUTINE GC_INITRUN_FORDRYDEP
     !   !-------------
     !   !XiaoLu, 2018/07
     !   !try to output landfield, it has to be done here after intht( )
     !   !--------------
     !   use gc_grid,         only : plon, plat, plev, plevp
     !   implicit none
     !   call drydep_landmap_init(plat,plon)

     !  ENDSUBROUTINE GC_INITRUN_FORDRYDEP
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        SUBROUTINE  drydep_landmap_init(platl,plonl)
      use gchp_netcdf
      use gc_grid,         only : plon, plat

      use ppgrid,           only: pcols, begchunk, endchunk
      use pmgrid,           only: masterproc
      use ioFileMod,        only: getfil
      use abortutils,       only: endrun
      use error_messages,   only: handle_ncerr
      use phys_grid,        only: scatter_field_to_chunk
      use history,          only: outfld
#if ( defined SPMD )
      use mpishorthand,         only: mpicom, mpiint
#endif

      implicit none

!---------------------------------------------------------------------------
!       ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: platl
      integer, intent(in) :: plonl
!---------------------------------------------------------------------------
!       ... Local variables
!---------------------------------------------------------------------------
      integer :: ncid, vid_ltfrac, nlat, nlon, ndims,ntype
      integer :: dimid_lat, dimid_lon, dimid_type
      integer :: ierr,ii,c
      integer :: count(3), start(3)

      character(len=256) :: locfn    ! netcdf local filename to open
      real, allocatable :: ltfrac_in(:,:,:)
!---------------------------------------------------------------------------
!       ... Open NetCDF file
!---------------------------------------------------------------------------
      if( masterproc ) then

      call getfil (file_land, locfn)
      call wrap_open(locfn,0,ncid)
!---------------------------------------------------------------------------
!       ... Get variable ID for dep vel array
!---------------------------------------------------------------------------
      call handle_ncerr( NF_INQ_VARID( ncid, 'LTFRAC', vid_ltfrac ),&
        'drydep_landmap_init: LTFRAC not found in input file' )
!---------------------------------------------------------------------------
!       ... Inquire about dimensions
!---------------------------------------------------------------------------
      call handle_ncerr( NF_INQ_DIMID( ncid, 'lon', dimid_lon ), &
        'drydep_landmap_init: getting lon dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lon, nlon ), &
        'drydep_landmap_init: getting nlon' )

      call handle_ncerr( NF_INQ_DIMID( ncid, 'lat', dimid_lat ), &
        'drydep_landmap_init: getting lat dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lat, nlat ), &
        'drydep_landmap_init: getting nlat' )

      call handle_ncerr( NF_INQ_DIMID( ncid, 'IOLSON', dimid_type ), &
        'drydep_landmap_init: getting landtype dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_type, ntype ),&
        'drydep_landmap_init: getting ntype' )

     write(*,*) 'BCC-AVIM land type: dimensions (nlon,nlat,ntype) =', &
         nlon,nlat,ntype

        end if !masterproc
!---------------------------------------------------------------------------
!       ... Check dimensions of dvel variable. Must be (lon, lat,
!       ntype).
!---------------------------------------------------------------------------
!   IGNORE


!---------------------------------------------------------------------------
!       ... Allocate and read data from file
!---------------------------------------------------------------------------
#if (defined SPMD )
            !write(*,*)'do we in SPMD?'.yes
            call mpibcast( nlon, 1, mpiint, 0, mpicom )
            call mpibcast( nlat, 1, mpiint, 0, mpicom )
            call mpibcast( ntype, 1, mpiint, 0, mpicom )
#endif


      allocate( ltfrac_in(nlon,nlat,ntype), stat=ierr )
         if( ierr /= 0 ) then
         write(*,*) 'drydep_landmap_init:: Failed to allocate ltfrac_in'
         call endrun
         end if


      if( masterproc ) then

      start = (/ 1, 1, 1 /)
      count = (/ nlon, nlat, ntype /)
      call handle_ncerr( NF_GET_VARA_DOUBLE( ncid, vid_ltfrac, start, &
      count, ltfrac_in ), &
                         'drydep_landmap_init: getting ltfrac_in' )

      end if !masterproc
!---------------------------------------------------------------------------
!       ... Regrid deposition velocities
!---------------------------------------------------------------------------
     allocate( ltfrac(pcols,begchunk:endchunk,ntype), stat=ierr)
      if( ierr /= 0 ) then
         write(*,*) 'landmap_inti: Failed to allocate landmap'
         call endrun
      end if

     call scatter_field_to_chunk(1,1,ntype,nlon,ltfrac_in(:,:,:),ltfrac(:,:,:))

     !do c=begchunk, endchunk
     !    call outfld('OceanLTF',ltfrac(:, c,2),pcols   ,c     )
     !end do
     !OK, we now make sure that oceanLTF is correct here.
      
        if( masterproc ) then
        deallocate( ltfrac_in )
        end if

        END SUBROUTINE drydep_landmap_init

      END MODULE GCHP_Landmap_Init
