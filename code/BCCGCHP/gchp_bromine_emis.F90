#include <misc.h>
#include <params.h>
!----------------------------------
!This is a program to read bromine emissions (CHBr3,CH2Br2)in GEOS-Chem chemistry
!Data comes from Bromocarb_Liang2010.nc and is converted by Xiao Lu
!Code is mimiced from the prescribed_em_d2_anthro.F90
!History:
!Xiao Lu, 2019/02: Initial version
!xiaolu,2019/02
!----------------------------------
#ifdef GEOSCHEM
module gchp_bromine_emis

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use pmgrid,         only: plon, plat, masterproc
   use ppgrid,         only: begchunk, endchunk, pcols, pver
   use phys_grid,      only: scatter_field_to_chunk
   use time_manager,   only: get_curr_date, get_curr_calday
   use abortutils,     only: endrun
   use prescribed_em_d2_anthro, only: emis2dcyc

#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

   public ini_bromine_emis
   public get_bromine_emis
!---
!Private data
!---

    integer varid_chbr3                           ! netcdf id for variables
    integer varid_ch2br2

    real(r8), private :: cdaym     ! dataset calendar day previous month
    real(r8), private :: cdayp     ! dataset calendar day next month

    integer  :: np1         ! current forward time index 
    integer  :: ncid   ! netcdf ID
    integer  :: timesiz     ! size of time dimension

     real(r8), private, allocatable ::  flux_chbr3(:,:)
     real(r8), private, allocatable ::  flux_ch2br2(:,:)
!=================================
contains
!=================================

   subroutine ini_bromine_emis

   use iofilemod,          only : getfil

!---
! Local workspace
!---
   integer londimid                        ! netcdf id for longitude dimension
   integer latdimid                        ! netcdf id for latitude dimension
   integer timeid                          ! netcdf id for time variable
   integer dimids(3)                       ! variable shape
   integer cnt3(3)                         ! array of counts for each dimension
   integer strt3(3)                        ! array of starting indices

   integer i, j, k, lat, n                    ! longitude, level, latitude, time indices
   integer  :: yr, mon, day                ! components of a date
   integer  :: ncsec                       ! current time of day [seconds]
   integer  :: istat
   real(r8) :: calday                      ! current calendar day
   real(r8) :: caldayloc                     ! calendar day (includes yr if no cycling)

   character(len=256) :: locfn, locfn_albedo        ! netcdf local filename to open
   character(len=256) :: filename , filename_albedo
   character(len=25)  :: species_name
   character(len=50)  :: units

   integer latsiz, lonsiz, idx
   real(r8), allocatable  ::  tmp_chbr3(:,:),tmp_ch2br2(:,:)

!---
! Allocate data
!---
   allocate( flux_chbr3(pcols,begchunk:endchunk), stat=istat )
   allocate( flux_ch2br2(pcols,begchunk:endchunk), stat=istat )

!
! SPMD: Master does all the work. Sends needed info to slaves
!
   units(:) = ' '

   if (masterproc) then

      !---
      !read nc data
      !---
    filename='/BIGDATA1/pku_atmos_lzhang_1/xpye/BCC-GEOS-Chem/data_process/geos-chem-emissions/t159_from_t42_Bromine_emission_forBCC.nc'
      call getfil (filename, locfn, 0)
      call wrap_open (locfn, 0, ncid)

!
! Use year information only if not cycling ozone dataset
!
      calday = get_curr_calday()
      call get_curr_date(yr, mon, day, ncsec)

!-------zf 2016.06.12
      if (emis2dcyc) then
        caldayloc = calday
      else
        caldayloc = calday + yr*365.
      end if

!---
! Get and check dimension info
!---
      CALL WRAP_INQ_DIMID( ncid, 'lon', londimid   )
      CALL WRAP_INQ_DIMID( ncid, 'lat', latdimid   )
      CALL WRAP_INQ_DIMID( ncid, 'time', timeid  )

      CALL WRAP_INQ_DIMLEN( ncid, londimid, lonsiz   )
      CALL WRAP_INQ_DIMLEN( ncid, latdimid, latsiz   )
      CALL WRAP_INQ_DIMLEN( ncid, timeid, timesiz   )

      CALL WRAP_INQ_VARID( ncid, 'CHBr3', varid_chbr3   )
      CALL WRAP_INQ_VARDIMID (ncid, varid_chbr3, dimids)

      CALL WRAP_INQ_VARID( ncid, 'CH2Br2', varid_ch2br2   )
      CALL WRAP_INQ_VARDIMID (ncid, varid_ch2br2, dimids)

      if (dimids(1) /= londimid .and. dimids(2) /= latdimid) then
         write(6,*)'ini_bromine_emis: Data must be ordered lon, lat, time'
         call endrun
      end if
      if (lonsiz .ne. plon .or. latsiz .ne. plat ) then
         call endrun ('ini_bromine_emis: The lonsiz and latsiz are not the same as the model resolution')
      end if

      CALL WRAP_INQ_VARID( ncid, 'time', timeid   )

   end if   ! masterproc

#if (defined SPMD )
   call mpibcast (timesiz, 1, mpiint, 0, mpicom)
#endif

!---
! Dynamically allocated memory for module data
!---
      allocate( tmp_chbr3(plon,plat), stat=istat )
      allocate( tmp_ch2br2(plon,plat), stat=istat )
      !allocate( time_emis(timesiz), stat=istat )


   if (masterproc) then
      strt3(1) = 1
      strt3(2) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = 1

      flux_chbr3=0.
      flux_ch2br2=0.
      strt3(3) = mon

           call wrap_inq_varid( ncid, 'CHBr3', varid_chbr3)
           call wrap_get_vara_realx (ncid, varid_chbr3, strt3, cnt3, tmp_chbr3)

           call wrap_inq_varid( ncid, 'CH2Br2', varid_ch2br2)
           call wrap_get_vara_realx (ncid, varid_ch2br2, strt3, cnt3, tmp_ch2br2)
   endif !masterproc

      call scatter_field_to_chunk(1,1,1,plon,tmp_chbr3,flux_chbr3(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plon,tmp_ch2br2,flux_ch2br2(1,begchunk))

   if (masterproc) then
       ! Close netcdf file
       call wrap_close(ncid)
   end if

#if (defined SPMD )
   call mpibcast (np1, 1, mpiint, 0, mpicom)
   call mpibcast (cdaym, 1, mpir8, 0, mpicom)
   call mpibcast (cdayp, 1, mpir8, 0, mpicom)
#endif

   deallocate(tmp_chbr3)
   deallocate(tmp_ch2br2)

end subroutine ini_bromine_emis

!******************************************************************
 subroutine get_bromine_emis(lchnk, ncol, gc_chbr3_input,gc_ch2br2_input)

!-----------------------------------------------------------------------
!
! Purpose: Return slice of bromine emission data
!
!-----------------------------------------------------------------------

  implicit none

!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns
  real(r8), intent(out)  :: gc_chbr3_input(pcols)  ! GC
  real(r8), intent(out)  :: gc_ch2br2_input(pcols)  ! GC
!
! Local variables.
!
   integer :: i, k, n                 ! longitude, level indices
   integer :: L, L1, L2
!-----------------------------------------------------------------------

      do i = 1, ncol
         gc_chbr3_input(i) = flux_chbr3(i,lchnk)
         gc_ch2br2_input(i) = flux_ch2br2(i,lchnk)
      end do

   return
 end subroutine  get_bromine_emis

!******************************************************************

end module gchp_bromine_emis
#endif


