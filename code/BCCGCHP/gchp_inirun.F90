!---------------------------------------
!xiaolu
!Purpose: GEOS-Chem-MOZART initialization
!from MOZART2 mo_inirun.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/08
!---------------------------------------

#include <misc.h>
#include <params.h>

#ifdef MOZART2

      module gchp_inirun

      implicit none

      save 

      contains

      subroutine gc_inirun

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid,          only : masterproc
      use gc_grid,         only : plon, plat, plev, plevp

      !use gchp_const_mozart, only : gc_constants_inti, ktop
      use mo_const_mozart, only : mo_constants_inti, ktop
      use gchp_chemini,      only : gc_chemini
!      use gchp_hook,         only : gc_hook_inti !xiaolu,2017/10/08,lightning
      use mo_hook,         only : moz_hook_inti
      !use gchp_drydep,       only : dvel_inti
      use gchp_surface,      only : inisflx

      !xiaolu add use gchp_sim_chm
!      use gchp_sim_chm,      only:sim_chm_data
      use mo_sim_chm,       only :sim_chm_data    !zf test
!      use hycoef
      !xiaolu add GCHP_Landmap_Init,2017/08
      !use GCHP_Landmap_Init
      use mo_ub_vals,    only : UB_INTI!sunrz add 2025
      implicit none

#include <comhyb.h>

      real, parameter :: p_limit = 10.e2     ! findsp pressure limit (Pa)
      real, parameter :: rair   = 287.04 
      real, parameter :: gravit = 9.80616
      real, parameter :: cpair  = 1004.64
      real, parameter :: latvap = 2.5104e06
      real, parameter :: rhoh2o = 1.e3

      logical   ::         diconvccm
      logical   ::         arconvccm

      real      ::         ref_pmid(plev)
      integer   ::         k

!---------------------------------------------------------------------
 
      arconvccm = .false.
      diconvccm = .true.      
!      call gc_constants_inti
       call mo_constants_inti      
!---------------------------------------------------------------------
!       ... Diagnostics
!---------------------------------------------------------------------
!      pdiags%negtrc  = negtrc_diagprnt
!      pdiags%imp_slv = imp_slv_diagprnt
!      pdiags%adv     = adv_diagprnt
!---------------------------------------------------------------------
!     	... Initialize surface flux calculation.
!---------------------------------------------------------------------
      !write(*,*) "In gchp_inirun.F90, xpye, before call inisflx"
      call inisflx( rair, gravit )
!---------------------------------------------------------------------
!       ... Initialize Hacks convective mass flux adjustment parameterization.
!---------------------------------------------------------------------
!zftest      if( arconvccm .or. diconvccm ) then
!zftest         call mfinti( hypi, rair, cpair, gravit, latvap, rhoh2o )
!zftest         write(*,*) 'inirun: finished MFINTI'
!zftest      end if
!---------------------------------------------------------------------
!       ... Initialize tropospheric chemistry calculation
!---------------------------------------------------------------------
!      call gc_chemini( plon, plat )
      call ub_inti( plon, plat ) !sunrz add 2025
!---------------------------------------------------------------------
!     	... Initialize deposition velocities
!--------------------------------------------------------------------- 
!xiaolu,2017/08, USE GEOS-Chem DRY DEP instead of BCC DRY DEP
!      call dvel_inti( plat, plon )
!sunrz comment it 
!      call drydep_landmap_init(plat,plon)
!---------------------------------------------------------------------
!       ... Find level where etamids are all > 10hPa
!---------------------------------------------------------------------
      !xiaolu,2017/10/08,lightning
!      call gc_hook_inti(1.2)
      call moz_hook_inti(1.2)
      do k = 1, plev
        ref_pmid(k) = (hyam(k) + hybm(k))*1.e5
      end do

      ktop = 0
      if( ref_pmid(1) < p_limit ) then
         do k = 1,plev
            if( ref_pmid(k) < p_limit ) then
              ktop = k
            end if
         end do
      end if
     
      !xiaolu:call sim_chm_data here (in gchp_sim_chm.F90)
      call sim_chm_data

      end subroutine gc_inirun

      end module gchp_inirun
#endif

