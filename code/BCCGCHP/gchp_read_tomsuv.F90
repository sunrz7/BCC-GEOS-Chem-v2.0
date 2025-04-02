#include <misc.h>
#include <params.h>
!----------------------------------
!This is a program to read TOMS-UV data to constrain the J-values in GEOS-Chem chemistry
!Data come from GEOS-Chem,and prepared by Xiao Lu
!Code is mimiced from the prescribed_em_d2_anthro.F90
!History:
!Xiao Lu, 2019/01: Initial version 
!Xiao Lu, 2019/01: Also read UV-Albedo in this routine
!xiaolu,2019/01
!-----------------------------------
#ifdef GEOSCHEM
module gchp_read_tomsuv
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

   public ini_tomsuv 
!   public int_tomsuv
   public get_tomsuv
   public get_uvalbedo
!---
!Public data
!---
   

!---
!Private data
!---

    integer varid                           ! netcdf id for variables
    integer varid_albedo                           ! netcdf id for variables

    real(r8), private :: cdaym     ! dataset calendar day previous month
    real(r8), private :: cdayp     ! dataset calendar day next month

    integer  :: np1         ! current forward time index of ozone dataset
    integer  :: ncid_toms   ! netcdf ID for TOMS
    integer  :: ncid_albedo   ! netcdf ID for UV-albedo
    integer  :: timesiz     ! size of time dimension on ozone dataset

    real(r8), private, allocatable :: gc_toms(:, :, :)!toms,toms1,toms2,dtoms1,dtoms2
    real(r8), private, allocatable :: gc_albedo(:,:) !for UValbedo

    real(r8), private, allocatable :: time_emis2(:)   ! Date on dataset (days since 0000 yr)

!=================================
contains
!=================================

   subroutine ini_tomsuv

   use iofilemod,          only : getfil

!---
! Local workspace
!---
   integer londimid                        ! netcdf id for longitude dimension
   integer latdimid                        ! netcdf id for latitude dimension
   integer timeid                          ! netcdf id for time variable
   integer dimids(3)                       ! variable shape
   integer dimids_albedo(3)                       ! variable shape
   integer cnt3(3)                         ! array of counts for each dimension
   integer strt3(3)                        ! array of starting indices
   integer cnt3_albedo(3)                         ! array of counts for each dimension
   integer strt3_albedo(3)                        ! array of starting indices

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
   real(r8), allocatable  ::  tmp_toms(:,:),tmp_albedo(:,:)
!---
! Allocate data
!---
   allocate( gc_toms(pcols,begchunk:endchunk,5), stat=istat )
   allocate( gc_albedo(pcols,begchunk:endchunk), stat=istat )

!
! SPMD: Master does all the work. Sends needed info to slaves
!
   units(:) = ' '

   if (masterproc) then

      !---
      !read nc data
      !---
    filename='/BIGDATA1/pku_atmos_lzhang_1/xpye/BCC-GEOS-Chem/data_process/TOMSUV/t159_from_t42_TOMS_O3col_for_BCCGC.nc'
      call getfil (filename, locfn, 0)
      call wrap_open (locfn, 0, ncid_toms)

    filename_albedo='/BIGDATA1/pku_atmos_lzhang_1/xpye/BCC-GEOS-Chem/data_process/TOMSUV/t159_from_t42_uvalbedo_for_BCCGC.nc'
      call getfil (filename_albedo, locfn_albedo, 0)
      call wrap_open (locfn_albedo, 0, ncid_albedo)

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
      CALL WRAP_INQ_DIMID( ncid_toms, 'lon', londimid   )
      CALL WRAP_INQ_DIMID( ncid_toms, 'lat', latdimid   )
      CALL WRAP_INQ_DIMID( ncid_toms, 'time', timeid  )

      CALL WRAP_INQ_DIMLEN( ncid_toms, londimid, lonsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_toms, latdimid, latsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_toms, timeid, timesiz   )

      CALL WRAP_INQ_VARID( ncid_toms, 'TOMS', varid   )
      CALL WRAP_INQ_VARDIMID (ncid_toms, varid, dimids)

      CALL WRAP_INQ_VARID( ncid_albedo, 'UVALBD', varid_albedo   )
      CALL WRAP_INQ_VARDIMID (ncid_albedo, varid_albedo, dimids_albedo)

      if (dimids(1) /= londimid .and. dimids(2) /= latdimid) then
         write(6,*)'ini_tomsuv: Data must be ordered lon, lat, time'
         call endrun
      end if
      if (lonsiz .ne. plon .or. latsiz .ne. plat ) then
         call endrun ('ini_tomsuv: The lonsiz and latsiz are not the same as the model resolution')
      end if

      CALL WRAP_INQ_VARID( ncid_toms, 'time', timeid   )

     ! CALL WRAP_GET_ATT_TEXT(ncid_toms, timeid, 'units', units )
     ! if ( trim(units) /= 'days since 1750-01-01 00:00:00' ) then
     !    call endrun ('em2ini_IPCC_anthro:  The time units are not days since 1750-01-01')
     ! end if

   end if   ! masterproc

#if (defined SPMD )
   call mpibcast (timesiz, 1, mpiint, 0, mpicom)
#endif


!---
! Dynamically allocated memory for module data
!---
      allocate( tmp_toms(plon,plat), stat=istat )
      allocate( tmp_albedo(plon,plat), stat=istat )
      allocate( time_emis2(timesiz), stat=istat )                   

   if (masterproc) then

!---
! Retrieve entire date and sec variables.
!---
      CALL WRAP_INQ_VARID( ncid_toms, 'time', timeid   )
      CALL WRAP_GET_VAR_REALX (ncid_toms,timeid,time_emis2)
      strt3(1) = 1
      strt3(2) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = 1
   endif   ! materproc
  
!---
! read data 
!---
!===============================================================
!the following for TOMS-UV
!===============================================================

   if (masterproc) then

      gc_toms=0. 
      gc_albedo=0.

     !
     ! Normal interpolation between consecutive time slices.
     !


      do n=1,timesiz-1
         np1 = n+1 
          cdaym = time_emis2(n) - 99*365.   !zf2017.06
          cdayp = time_emis2(n+1) - 99*365.   !zf2017.06

!---------zf 2016.06.13
         if (.not.emis2dcyc ) then
           cdaym = cdaym + 365.0*1849.0   !1850.0   zf 2017.06
           cdayp = cdayp + 365.0*1849.0   !1850.0   zf 2017.06
         end if

!---------------

         if (caldayloc >= cdaym .and. caldayloc < cdayp) then
!!           write(*,*)'xiaolu check toms day',caldayloc,cdaym,cdayp
           strt3(3) = n
!!               strt3(3) = 1
           call wrap_inq_varid( ncid_toms, 'TOMS', varid)
           call wrap_get_vara_realx (ncid_toms, varid, strt3, cnt3, tmp_toms)

!!           write(*,*)'xiaolu check toms data',minval(tmp_toms),maxval(tmp_toms)

          goto 10
         end if

     end do

      write(6,*)'em2ini_IPCC_anthro: Failed to find dates bracketing , year, month, day',  yr, mon, day
      call endrun
!------------------------------------------------
10    continue

   endif   ! masterproc

   call scatter_field_to_chunk(1,1,1,plon,tmp_toms,gc_toms(1,begchunk,1))

!===============================================================
!the following for albedo
!===============================================================
  if (masterproc) then
      strt3_albedo(1) = 1
      strt3_albedo(2) = 1
      cnt3_albedo(1)  = lonsiz
      cnt3_albedo(2)  = latsiz
      cnt3_albedo(3)  = 1

      gc_albedo=0.
      strt3_albedo(3) = mon

!      write(*,*)'xiaolu check mon',mon
           call wrap_inq_varid( ncid_albedo, 'UVALBD', varid_albedo)
           call wrap_get_vara_realx (ncid_albedo, varid_albedo, strt3_albedo, cnt3_albedo, tmp_albedo)

!      write(*,*)'xiaolu check tmp_albedo',maxval(tmp_albedo)

   endif   ! masterproc

   call scatter_field_to_chunk(1,1,1,plon,tmp_albedo,gc_albedo(1,begchunk))

  
   !xiaolu,2019/01/26
   if (masterproc) then 
       ! Close netcdf file
       call wrap_close(ncid_albedo)
       call wrap_close(ncid_toms)
   end if

#if (defined SPMD )
   call mpibcast (np1, 1, mpiint, 0, mpicom)
   call mpibcast (time_emis2, timesiz, mpir8, 0, mpicom )
   call mpibcast (cdaym, 1, mpir8, 0, mpicom)
   call mpibcast (cdayp, 1, mpir8, 0, mpicom)
#endif

   deallocate(tmp_toms)
   deallocate(tmp_albedo)

end subroutine ini_tomsuv

!*******************************************************************
!*******************************************************************
 subroutine get_tomsuv( lchnk, ncol, gc_toms_input)

!-----------------------------------------------------------------------
!
! Purpose: Return slice of TOMS-UV data
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns
  real(r8), intent(out)  :: gc_toms_input(pcols,5)  ! GCtoms
!
! Local variables.
!
   integer :: i, k, n                 ! longitude, level indices
   integer :: L, L1, L2
!-----------------------------------------------------------------------

!!!   do n = 1, 5
      n=1
      do i = 1, ncol
         gc_toms_input(i,n) = gc_toms(i,lchnk,n)
      end do
!!!   enddo

   return
 end subroutine  get_tomsuv

!******************************************************************

!*******************************************************************
 subroutine get_uvalbedo( lchnk, ncol, gc_albedo_input)

!-----------------------------------------------------------------------
!
! Purpose: Return slice of TOMS-UV data
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns
  real(r8), intent(out)  :: gc_albedo_input(pcols)  ! GCtoms
!
! Local variables.
!
   integer :: i, k, n                 ! longitude, level indices
   integer :: L, L1, L2
!-----------------------------------------------------------------------

      do i = 1, ncol
         gc_albedo_input(i) = gc_albedo(i,lchnk)
      end do

   return
 end subroutine  get_uvalbedo

!******************************************************************


end module gchp_read_tomsuv
#endif
