!xiaolu
!---------------------------------------
!Purpose: GEOS-Chem-MOZART initialization
!from MOZART2 mo_chemini.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/08
!add using dust_initialization
!---------------------------------------

#include <misc.h>
#include <params.h>

#ifdef MOZART2

      module gchp_chemini

      private
      public :: gc_chemini

      contains

      subroutine gc_chemini( plonl, platl )
!-----------------------------------------------------------------------
! 	... Chemistry module intialization
!-----------------------------------------------------------------------

!      use mo_ub_vals,    only : UB_INTI
      use gchp_surf,       only : SURF_INTI
!      use mo_airplane,   only : airpl_src
!      use mo_photo,      only : PRATE_INTI
!      use mo_chem_utls,  only : chem_utls_inti
      use gchp_srf_emis,   only : srf_emis_inti
!      use mo_sethet,     only : sethet_inti
!      use mo_setext,     only : setext_inti
!      use mo_usrrxt,     only : usrrxt_inti
!      use mo_grp_ratios, only : set_grp_ratios_inti
!      use mo_sulf,       only : sulf_inti
!      use mo_chem_mods,  only : grpcnt, clscnt1, clscnt4, clscnt5
      use gc_grid,       only : pcnstm1
!      use MO_EXP_SOL,    only : EXP_SLV_INTI
!      use MO_IMP_SOL,    only : IMP_SLV_INTI

      use gchp_tracinp,         only : bndic
!      use bcc_solar_data,     only : solar_data_init
!      use bcc_photo,          only : photo_inti
      use gchp_bcc_dust_intr,      only : dust_initialize
!      use bcc_aerosols_intr,  only : mz_aero_initialize
!      use bcc_setsox,         only : sox_inti
!      use bcc_solar_parms,    only : solar_parms_init

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl, platl

    integer :: i
    character(len=256) :: emisspath
    character(len=256) :: solar_data_file
    character(len=256) :: euvacdat_file
    character(len=256) :: photon_file
    character(len=256) :: electron_file
    character(len=256) :: solar_parms_file
    character(len=256) :: xs_coef_file
    character(len=256) :: xs_short_file
    character(len=256) :: xs_long_file
    character(len=256) :: rsf_file
    character(len=256) :: exo_coldens_file
    character(len=256) :: soil_erod_file

!-----------------------------------------------------------------------
! 	... Initialize photorate module
!-----------------------------------------------------------------------
!      call prate_inti

!-----------------------------------------------------------------------
! 	... Read time-independent airplane emissions
!-----------------------------------------------------------------------
!      call airpl_src(  plonl, platl )

!-----------------------------------------------------------------------
! 	... Initialize the chem utils module
!-----------------------------------------------------------------------
!      call chem_utls_inti

!-----------------------------------------------------------------------
! 	... Read time-dependent surface flux dataset
!-----------------------------------------------------------------------
!      call srf_emis_inti( plonl, platl, 1 )

!-----------------------------------------------------------------------
! 	... Intialize the het rates module
!-----------------------------------------------------------------------
 !     call sethet_inti

!-----------------------------------------------------------------------
! 	... Intialize the ext frcing module
!-----------------------------------------------------------------------
!      call setext_inti

!-----------------------------------------------------------------------
! 	... Intialize the rxt rate constant module
!-----------------------------------------------------------------------
!      call usrrxt_inti

!!-----------------------------------------------------------------------
!! 	... Intialize the grp ratios module
!!-----------------------------------------------------------------------
!      call set_grp_ratios_inti

!-----------------------------------------------------------------------
! 	... Read time-dependent surface variables dataset
!-----------------------------------------------------------------------
      call surf_inti( plonl, platl )

!-----------------------------------------------------------------------
! 	... Read time-dependent upper boundary values
!-----------------------------------------------------------------------

!     call ub_inti( plonl, platl )

!!-----------------------------------------------------------------------
!! 	... Read time-dependent sulfate dataset
!!	    NOTE : This is now a netcdf dataset
!!-----------------------------------------------------------------------
    
!      call sulf_inti( plonl, platl )

!      if( clscnt1 > 0 ) then
!!-----------------------------------------------------------------------
!!	... Initialize the explicit solver
!!-----------------------------------------------------------------------
!         call exp_slv_inti
!      end if
!      if( clscnt4 > 0 ) then
!!-----------------------------------------------------------------------
!!	... Initialize the implicit solver
!!-----------------------------------------------------------------------
!         call imp_slv_inti
!      end if

!!----------------------------------------------
!      do i = len_trim(bndic), 1, -1
!         if(bndic(i:i).eq.'/') exit
!      enddo
!      emisspath = bndic(1:i)
   
!       solar_parms_file = trim(emisspath)//'non-resolution-data/proxy_solar_Solomon_Richmond_1845-2008_daily_noleap_c110526.nc' 
!       solar_data_file= trim(emisspath)//'non-resolution-data/spectral_irradiance_Lean_1610-2009_ann_c100405.nc'   
!       xs_long_file   = trim(emisspath)//'non-resolution-data/temp_prs_GT200nm_jpl06_c080930.nc'
!       xs_coef_file   = trim(emisspath)//'non-resolution-data/effxstex.txt'
!       xs_short_file  = trim(emisspath)//'non-resolution-data/xs_short_jpl06_c080930.nc'
!       rsf_file       = trim(emisspath)//'non-resolution-data/RSF_GT200nm_v3.0_c080811.nc'
!       euvacdat_file  = trim(emisspath)//'non-resolution-data/euvac_v1.2.dat'
!       photon_file    = trim(emisspath)//'non-resolution-data/photon.dat'
!       electron_file  = trim(emisspath)//'non-resolution-data/electron.dat'
!       exo_coldens_file = trim(emisspath)//'resolution-data/exo_coldens.nc'

!       call solar_parms_init(solar_parms_file)
!       call solar_data_init(solar_data_file)

!       call sox_inti
!       call photo_inti( solar_data_file, &
!                        xs_coef_file,  xs_short_file, xs_long_file, rsf_file,         &
!                        euvacdat_file, photon_file,  electron_file, exo_coldens_file  )
!
!!---------------------------------------------------------------------
!!  ... initialize dust erodibility
!!--------------------------------------------------------------------
!!    soil_erod_file = trim(emisspath)//'dst_64x128_t42_c20130311.nc'
!!    soil_erod_file = trim(emisspath)//'dst_160x320_t106_c20140910.nc'

    !uncomment,xiaolu,2017/06/16
!    soil_erod_file = trim(emisspath)//'dst_t159.nc'
!    call dust_initialize(soil_erod_file)

!    call mz_aero_initialize

!--------------------------------------------------

      end subroutine gc_chemini

      end module gchp_chemini
#endif

