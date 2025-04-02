!--------------------------------------
!Purpose:
!from MOZART2 mo_tracinp.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/13
!---------------------------------------


      module GCHP_TRACINP
!-----------------------------------------------------------------------
! 	... Tracer initial conditions
!-----------------------------------------------------------------------

      use gc_grid,         only : plon, plat, pcnstm1
      use pmgrid,          only : plev    !add by zf 2008.07.02
      use gchp_netcdf

      implicit none

      save

      public  :: initrac

      real , public                :: ic_vmr_int(plon,plat,plev,pcnstm1)   
      character(len=256) , public  ::  bndic

      integer :: ncid            ! netcdf file id

      CONTAINS

      subroutine initrac
!-----------------------------------------------------------------------
! 	... Initialize /tracinp/, and tracer arrays.   The tracer fields
!           must come directly from a time sample that exists on the input
!           file (though not necessarily the first time sample).
!-----------------------------------------------------------------------


!add by zf 2008.06.13
      use pmgrid,         only: masterproc
      use ioFileMod,      only: getfil
      use abortutils,     only: endrun
      use error_messages, only: handle_ncerr
!end zf

      implicit none

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: astat,dimid, varid,i, j, k, m

      character(len=256) :: locfn    ! netcdf local filename to open   add by zf
!-----------------------------------------------------------------------
!     	... Open netcdf ic file
!-----------------------------------------------------------------------

    if (masterproc) then
      call getfil (bndic, locfn)
      call wrap_open(locfn,0,ncid)

      do i = len_trim(bndic), 1, -1
         if(bndic(i:i).eq.'/') exit
      enddo

!-----------------------------------------------------------------------
!   	... Allocate memory for ic and ic surface pressure
!-----------------------------------------------------------------------
      call GETTRAC( ic_vmr_int )

!--------------------------------------------------------------------------
!  	... Close file
!--------------------------------------------------------------------------
      call HANDLE_NCERR( NF_CLOSE( ncid ), 'INITRAC: Failed to close file ' // bndic  )

      end if          !  masterproc

      end subroutine INITRAC

      subroutine GETTRAC( ic_vmr )
!-----------------------------------------------------------------------
! 	... Read constituent data from netcdf intial condition file
!-----------------------------------------------------------------------

      use netcdf
      use gchp_tracname, only : tracnam
      use gc_grid,     only : pcnstm1
      use pmgrid,      only : masterproc     !add by zf 2008.06.13
      use error_messages, only: handle_ncerr !add xiaolu,2017/06/13

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      real, intent(out) :: ic_vmr(plon,plat,plev,pcnstm1)     ! ic initial concentrations

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer          ::  j, k, m              ! indicies
      integer          ::  varid                ! variable id
      integer          ::  ncret                ! netcdf return code
      integer          ::  start3(3), count3(3)
      character(len=8) :: fldname

    if (masterproc) then
      start3(:) = (/ 1,1,1 /)
      count3(:) = (/ plon,plat,plev /)

!-----------------------------------------------------------------------
! 	... Get transported variables
!-----------------------------------------------------------------------
      do m = 1,pcnstm1
	    fldname = tracnam(m)

            ncret = NF_INQ_VARID( ncid, TRIM( fldname ), varid )
            if( ncret /= NF_NOERR ) then
                write(*,*) 'GETTRAC: Failed to get variable id for ' // TRIM( fldname ) &
                          // ', will initialize with 10^-38'
                ic_vmr(:,:,:,m) = 1.e-38
            else
                call HANDLE_NCERR( NF_GET_VARA_DOUBLE( ncid, varid, start3, count3, ic_vmr(1,1,1,m) ), &
                            'GETTRAC: Failed to get variable ' // TRIM( fldname ) )
                write(*,*) 'GETTRAC: Loaded variable ' // TRIM( fldname ) // '.'
            end if

      end do

      end if          !  masterproc

      end subroutine GETTRAC

      end module GCHP_TRACINP
