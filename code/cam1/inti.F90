#include <misc.h>
#include <params.h>

subroutine inti ()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set constants and call initialization procedures for time independent
! physics routines
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,             only: plev, plevp          ! Needed for hypm passed to vd_inti
!------ zf 2008.03
#ifdef MOZART2
!   use chemistry,          only: trace_gas, chem_init, mozart2_init
!   use mo_read_sim_chm,    only: read_sim_chm
!   use prescribed_ch4,     only: ch4ini           !zf 2014.12.25
    use gigc_chunk_mod,     only: gigc_chunk_init
    use gchp_states
    use state_grid_mod, only: init_state_grid ! sunrz 2023/11 for V14.1.1
    USE bccgchp_vars_mod
    use physics_types,   only: physics_state
    use time_manager,    only: get_curr_calday,dtime,start_ymd,start_tod,stop_ymd,stop_tod
    use ioFileMod,       only: getfil
    use filenames,       only: bndtvoxid
    use acbnd,           only: acbndini
    use pmgrid,          only: masterproc,iam,plon,plat
    use ppgrid,          only: pcols,begchunk,endchunk
    use chemistry,       only: trace_gas,chem_init
    use phys_grid, only:get_lon_all_p,get_lat_all_p,get_rlon_all_p,get_rlat_all_p,clon_p,clat_p,beglat,endlat,nlon_p
    use gc_gridarea,     only:gc_calcarea
    !----------------------xiaolu,2017/06/13:add for emission
    use chemistry,         only: gchp_init
    use hemco_interface,  only: HCOI_Chunk_Init !sunrz 2024/7
#else
   use chemistry,          only: trace_gas, chem_init
#endif
#ifdef BCCCHEM
   use chemistry,          only: bccchem_init
#endif
#ifdef WACCM
   use mo_chemini,         only: chemini
   use iondrag,            only: iondrag_init
#endif
#ifdef GEOSCHEM
   use chemistry,          only: geoschem_init
#endif

#ifdef DO_COSP
   use cospsimulator_main,    only: cospsimulator_intr_init
#endif
   use tracers,            only: tracers_flag, tracers_init
   use aerosol_intr,       only: prognostic_aerosol_initialize
   use gw_drag,            only: gw_inti
   use vertical_diffusion, only: vd_inti
   use rayleigh_friction,  only: do_rayleigh, rayleigh_friction_init !add by luyx20170119
   use moistconvection,    only: mfinti
#ifdef MICROP-MG
   use cloud_fraction,     only: microp_cldfrc_init
   use microp_driver,      only: microp_driver_init
   use conv_water,         only: conv_water_init
#else
   use cloud_fraction,     only: cldfrc_init
#endif
   use cldcond,            only: cldcond_init
   use param_cldoptics,    only: param_cldoptics_init
   use wu_conv,            only: wu_convi
   use shr_const_mod,      only: shr_const_zvir, shr_const_cpwv, shr_const_rwv
   use physconst,          only: rair, cpair, cpwv, gravit, stebol, epsilo, tmelt, &
                                 latvap, latice, rh2o, zvir, cpvir, rhoh2o, pstd,  &
                                 karman, rhodair
#ifdef MAC_SP
   use prescribed_aerosols,     only: aerosol_initialize
   use prescribed_aerosols_cdnc,only: aerosol_initialize_cdnc
#else
   use prescribed_aerosols,     only: aerosol_initialize
#endif
   use aer_optics,         only: aer_optics_initialize
   use check_energy,       only: check_energy_init

#if (defined BCCCHEM) || (defined WACCM)
   use ioFileMod,         only: getfil
   use acbnd,             only: acbndini
   use time_manager,      only: get_curr_calday
   use filenames,         only: bndtvoxid

   use prescribed_em_d3,  only: em3ini_IPCC
   use prescribed_em_d2,  only: em2ini_IPCC
   use prescribed_em_d2_anthro,  only: em2ini_IPCC_anthro    ! zf2016.12.08
   use prescribed_emis_3d,  only: em3d_ini_IPCC              ! zf2017.01.04
   use prescribed_fco2_data,  only: fco2ini_IPCC
   use prescribed_fco2_anthro,  only: fco2ini_IPCC_anthro
#endif

#ifdef CO2
   use prescribed_fco2_data,  only: fco2ini_IPCC
   use prescribed_fco2_anthro, only: fco2ini_IPCC_anthro
   use prescribed_emis_3d,     only: em3d_ini_IPCC  
#endif
!----------------
! wutw 2014.3.3
   use dewpoint, only : dewp_init
!-------------------------
! Jing X.W.
#ifdef BCC_RAD
 use BCC_RAD_init_mod
#endif
!--------------
#ifdef MAC_SP
   use mo_simple_plumes, only: sp_setup
#endif

   implicit none

#include <comctl.h>
#include <comhyb.h>
!-------------------------------------------------
   character*4 ncnam(5)                           ! names of oxidants
   character(len=256) :: oxid_file                ! netcdf local filename for oxidants
   real(r8) :: calday  

#ifdef WACCM
  character(len=256) :: efield_lflux_file
  character(len=256) :: efield_hflux_file
  character(len=256) :: efield_wei96_file
#endif
 #if (defined MOZART2)
    integer  :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
     REAL(4) ,allocatable   :: lonCtr(:,:) ! Lon centers [radians]
     REAL(4) ,allocatable   :: latCtr(:,:) ! Lat centers [radians]
     REAL(4),allocatable           :: lonEdge(:,:)   ! Lon ctrs [deg] from ESMF
     REAL(4),allocatable           :: latEdge(:,:)   ! Lat ctrs [deg] from ESMF
    integer ::I_LO
    integer ::J_LO
    integer ::I_HI
    integer ::J_HI
    integer ::i,j,II   ! loop index
    integer ::nlat
    integer ::iam_bcc
 
    real(r8), dimension(:), allocatable :: area_all
    REAL    ::gctime
    type(physics_state) :: state
    INTEGER  :: RC            ! Success or failure
 #endif
!----------------------------------------------------------
!
! Initialize radiation data for aerosol forcing calculation
! Initialize aerosol fields from files
!
   call aer_optics_initialize()

#ifdef MAC_SP
   call aerosol_initialize()
   call aerosol_initialize_cdnc()
#else
   call aerosol_initialize()
#endif   

   call prognostic_aerosol_initialize()

!------------------------------------------------------------
! read cospsimulator_nl.txt
!
#ifdef DO_COSP
   call cospsimulator_intr_init
#endif   
!
!-----------------------------------------------------------------------
!
! Initialize physconst variables
! In adiabatic case, set zvir and cpvir explicitly to zero instead of 
! computing as (rh2o/rair - 1.) and (cpwv/cpair - 1.) respectively, in order 
! to guarantee an identical zero.
!
   if (adiabatic .or. ideal_phys) then
      rh2o  = rair
      zvir  = 0.
      cpwv  = cpair
      cpvir = 0.
   else
      rh2o  = shr_const_rwv
      zvir  = shr_const_zvir
      cpwv  = shr_const_cpwv
      cpvir = cpwv/cpair - 1.
   end if
!
! Call time independent initialization routines for parameterizations.
!
#ifdef MOZART2
!   select case (nsrest)
!   case (1)
!   call read_sim_chm
!   end select

!   call mozart2_init           !zf 2008.03
!   call fco2ini_IPCC
!   call fco2ini_IPCC_anthro

   calday = get_curr_calday()
!
!   call em3ini_IPCC
!   call em2ini_IPCC
!   call em2ini_IPCC_anthro    !zf 2016.12.08
!   call em3d_ini_IPCC           !zf 2017.01.04
!   call ch4ini           !zf20141225

    ncol  = 16
    allocate (lonCtr(ncol,1))
    allocate (latCtr(ncol,1))
    allocate (lonEdge(ncol+1,1))
    allocate (latEdge(ncol,2))
    gctime=dtime
    call gchp_init
     ! Initialize fields of the Grid State object
     CALL Init_State_Grid( Input_Opt, State_Grid, RC )
 
     ! Pass grid information obtained from Extract_ to State_Grid
     State_Grid%NX          = ncol            ! # lons   on this PET
     State_Grid%NY          = 1            ! # lats   on this PET
     State_Grid%NZ          = plev            ! # levels on this PET
     State_Grid%GlobalNX    = plon      ! # lons   in global grid
     State_Grid%GlobalNY    = plat      ! # lats   in global grid
     State_Grid%NativeNZ    = plev      ! # levels in global grid
     State_Grid%XMinOffset  = 1             ! X offset from global grid
     State_Grid%XMaxOffset  = State_Grid%NX ! X offset from global grid
     State_Grid%YMinOffset  = 1             ! Y offset from global grid
     State_Grid%YMaxOffset  = State_Grid%NY ! Y offset from global grid
     State_Grid%MaxTropLev  = 40            ! # trop. levels
     State_Grid%MaxStratLev = plev          ! # strat. levels
 
 
     call Set_Input_Opt( masterproc, Input_Opt, RC )
     CALL GIGC_Chunk_Init( nymdB     = start_ymd,      & ! YYYYMMDD @ start of run
                           nhmsB     = start_tod,      & ! hhmmss   @ start of run
                           nymdE     = stop_ymd,      & ! YYYYMMDD @ end of run
                           nhmsE     = stop_tod,      & ! hhmmss   @ end of run
                           tsChem    = gctime,     & ! Chemical timestep [s]
                           tsDyn     = gctime,      & ! Dynamic  timestep [s]
                           tsRad     = gctime,      & ! RRTMG    timestep [s]
                           lonCtr    = lonCtr,     & ! Lon centers [radians]
                           latCtr    = latCtr,     & ! Lat centers [radians]
                           lonEdge   = lonEdge,    &
                           latEdge   = latEdge,    &
                           Input_Opt = Input_Opt,  & ! Input Options obj
                           State_Chm = State_Chm,  & ! Chemistry State obj
                           State_Diag= State_Diag, & ! Diagnostics State obj
                           State_Grid= State_Grid, & ! Grid State obj
                           State_Met = State_Met,  & ! Meteorology State obj
                           HcoConfig = HcoConfig,  & ! HEMCO config obj
                           RC        = RC        )
 
 
     call HCOI_Chunk_Init()
#endif

#ifdef CO2
   call fco2ini_IPCC
   call fco2ini_IPCC_anthro
   call em3d_ini_IPCC
#endif
!----------------------------
! wutw 2014.0303  for S2S

    call dewp_init

!-------------------------

#ifdef BCCCHEM
    call em3ini_IPCC
    call em2ini_IPCC
    call fco2ini_IPCC
    call bccchem_init
#endif

#ifdef WACCM
    call em3ini_IPCC
    call em2ini_IPCC
!   call fco2ini_IPCC

    call chemini

#endif

!--------------------------------------------------------
#if (defined GEOSCHEM)

    calday = get_curr_calday()
    ncnam(1) = 'OH'
    ncnam(2) = 'HO2'

    call getfil(bndtvoxid, oxid_file, 0)
    call acbndini( oxid_file, calday, 2, ncnam )
    write (6,*) ' opened oxidant file ',trim(oxid_file)

   call em3ini_IPCC
   call em2ini_IPCC
#endif

!-------------------------------------------------

   if (trace_gas) call chem_init
   if (tracers_flag) call tracers_init
   if (tracers_flag) call tracers_init
#ifdef BCC_DRAG
   call gw_inti (cpair   ,cpwv    ,gravit  ,rair    ,hypi   , &
                 hypm    )
#else
   call gw_inti (cpair   ,cpwv    ,gravit  ,rair    ,hypi    )
#endif
   if (do_rayleigh) then
      call rayleigh_friction_init() ! add by luyx20170119
   endif 
   call vd_inti (cpair   ,cpwv    ,gravit  ,rair    ,zvir   , &
                 hypm    ,karman    )
   call tsinti  (tmelt   ,latvap  ,rair    ,stebol  ,latice  )
   call radini  (gravit  ,cpair   ,epsilo  ,stebol  ,pstd*10.0 )
   call esinti  (epsilo  ,latvap  ,latice  ,rh2o    ,cpair  , &
                 tmelt   )
   call mfinti  (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )

#ifdef MICROP-MG
   call microp_cldfrc_init
   call microp_driver_init
   call conv_water_init
#else
   call cldfrc_init
#endif
   call wu_convi( tmelt, epsilo, latvap, cpair )
   call cldinti ()

   call cldcond_init
   call param_cldoptics_init
   call check_energy_init

#ifdef BCC_RAD
   call BCC_RAD_init
#endif

#ifdef MAC_SP
   call sp_setup
#endif

   return
end subroutine inti
