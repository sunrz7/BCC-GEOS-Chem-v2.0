#include <misc.h>
#include <params.h>

subroutine tphysbc (ztodt,   pblht,   kbfs ,   tpert,   ts,      sst,     &
                    qpert,   precl,   precc,   precsl,  precsc,  &
                    asdir,   asdif,   aldir,   aldif,   snowh,   &
                    qrs,     qrl,     flwds,   fsns,    fsnt,    &
                    flns,    flnt,    lwup,    srfrad,  sols,    &
                    soll,    solsd,   solld,   state,   tend,    &
                    pbuf,    prcsnw,  fsds ,   landm,   landfrac,&
                    ocnfrac, icefrac, fhs  ,                     &
                    cosp_cnt,                                    &
#ifdef NSTEP_MOIST
                    cnt,     cnb,     cmfdqr,  nevapr,  prain   ,&
                    cmfmc,                                       &
                    dt_moist, dq_moist,    du_moist, dv_moist, &
                    q0_conv,  t0_conv,     date0_conv,        &
#endif
#ifdef MOZART2
                    chemflx, dq_chem,                          &
                    flux_ISOP, flux_ACET,                 &
                    flux_C3H6, flux_C2H4,                 &
                    flux_OC2,  flux_C10H16,               &
                    flux_CO2,  flux_N2O, flux_DMS,        &
                    co2mmr_latmean,shf,lhf,u10,v10,taux,tauy,eflx,lands,lais,area,   &
#endif
#ifdef CO2
                    co2mmr_latmean,                       &
#endif
#ifdef HAMOCC
                    dmsvmr_a,  n2ovmr_a,                  &
#endif
                    co2vmr_a, dyn_tend)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics BEFORE coupling to land, sea, and ice models.
! 
! Method: 
! Call physics subroutines and compute the following:
!     o cloud calculations (cloud fraction, emissivity, etc.)
!     o radiation calculations
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------


   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use param_cldoptics, only: param_cldoptics_calc
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, &
                              dynamics_tend, physics_ptend_init
   use diagnostics,     only: diag_dynvar
   use history,         only: outfld
   use physconst,       only: gravit, cpair, tmelt, cappa, zvir, rair, rga
   use radheat,         only: radheat_net
   use constituents,    only: pcnst, pnats, ppcnst, qmin, iaero

   use dycore,          only: get_resolution

   use constituents,    only: cnst_get_ind
   use time_manager,    only: is_first_step, get_nstep, get_curr_calday
   use check_energy,    only: check_energy_chng, check_energy_fix
   use check_energy,    only: check_tracers_data, check_tracers_init
   use dycore,          only: dycore_is
#ifdef DO_COSP
   use cospsimulator_main, only: cospsimulator_intr_run, cosp_nradsteps, first_run_cosp
#endif
   use abortutils,      only: endrun
   use ghg_surfvals,    only: chem_surfvals_get       ! for carbon cycle, wzz, 2008.11.20
   use ghg_surfvals_2d,   only: ghg_surfvals_get_co2vmr_2d
   use wu_cdnc,        only: get_cdnc
   use physconst,      only: mwdry, mwch4, mwn2o, mwco2, mwf11, mwf12

#ifdef MOZART2
    use pmgrid,             only: masterproc
    use wv_saturation, only: aqsat !to calculate relative humidity
    use gchp_states
    use chemistry,       only: ncnst
    use gc_gridarea, only:gc_calcarea
    use hco_bcc_convert_state_mod, only: State_CAM_AREAM2  
    use mo_ub_vals,    only : set_ub_vals ! sunrz 2025
#endif

#ifdef CO2
   use prescribed_emis_3d, only: getem3_IPCC6, nem3_IPCC6
   use bcc_co2_em3_prod,   only: co2_em3_prod
#endif

#ifdef BCCCHEM
   use bcc_chem_mods,          only: adv_mass
   use bcc_chem_utls,          only : get_spc_ndx
#endif

#ifdef WACCM
   use chem_mods,          only: adv_mass
   use m_spc_id,           only: id_CO2
#endif

#ifdef GEOSCHEM
   use acbnd,            only: acbndget
   use prescribed_em_d3, only: getem3_IPCC, nem3
   use prescribed_em_d2, only: getem2_IPCC, nem2, varname2d
   use prescribed_em_d2, only: idxCO, idxCH4, idxNO, idxC2H4
   use wv_saturation,    only: aqsat
   use chemistry,        only: ncnst, ixchm
   use gc_chemdr,        only: do_gc_chem, molar_mass
#endif

#if (defined BCCCHEM) 
   use bcc_chemdr,  only : bcc_chemdr_ub
#endif

   implicit none

#include <comctl.h>

   real(r8), parameter :: rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
   real(r8), parameter :: rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
   real(r8), parameter :: rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
   real(r8), parameter :: rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
   real(r8), parameter :: rmwco2 = mwco2/mwdry      ! ratio of molecular weights of co2 to dry air
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: ts(pcols)                      ! surface temperature
   real(r8), intent(in) :: sst(pcols)                     ! sea surface temperature
   real(r8), intent(in) :: kbfs(pcols)                     ! surface turbance flux
   real(r8), intent(inout) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(inout) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess
   real(r8), intent(in) :: asdir(pcols)                  ! Albedo: shortwave, direct
   real(r8), intent(in) :: asdif(pcols)                  ! Albedo: shortwave, diffuse
   real(r8), intent(in) :: aldir(pcols)                  ! Albedo: longwave, direct
   real(r8), intent(in) :: aldif(pcols)                  ! Albedo: longwave, diffuse
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: qrs(pcols,pver)            ! Shortwave heating rate
   real(r8), intent(inout) :: qrl(pcols,pver)            ! Longwave  heating rate
   real(r8), intent(inout) :: flwds(pcols)               ! Surface longwave down flux
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: lwup(pcols)                    ! Surface longwave up flux
   real(r8), intent(out) :: srfrad(pcols)                 ! Net surface radiative flux (watts/m**2)
   real(r8), intent(inout) :: sols(pcols)                 ! Direct beam solar rad. onto srf (sw)
   real(r8), intent(inout) :: soll(pcols)                 ! Direct beam solar rad. onto srf (lw)
   real(r8), intent(inout) :: solsd(pcols)                ! Diffuse solar radiation onto srf (sw)
   real(r8), intent(inout) :: solld(pcols)                ! Diffuse solar radiation onto srf (lw)
   real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
   real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
   real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
   real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
   real(r8), intent(out) :: prcsnw(pcols)                 ! snowfall rate (precsl + precsc)
   real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction
   real(r8), intent(in) :: ocnfrac(pcols)                 ! land fraction
   real(r8), intent(in) :: icefrac(pcols)                 ! land fraction
   real(r8), intent(in) :: fhs(pcols)    ! soil moisture content, wzz, 2008.3.12

   integer, intent(inout) :: cosp_cnt

#ifdef NSTEP_MOIST
   real(r8), intent(inout) :: cnt(pcols)                        ! Top level of convective activity
   real(r8), intent(inout) :: cnb(pcols)                        ! Lowest level of convective activity
   real(r8), intent(inout) :: cmfdqr(pcols,pver)
   real(r8), intent(inout) :: nevapr(pcols,pver)
   real(r8), intent(inout) :: prain(pcols,pver)
   real(r8), intent(inout) :: cmfmc(pcols,pverp)
   real(r8), intent(inout) :: dq_moist(pcols,pver,ppcnst)
   real(r8), intent(inout) :: du_moist(pcols,pver)
   real(r8), intent(inout) :: dv_moist(pcols,pver)
   real(r8), intent(inout) :: dt_moist(pcols,pver)
   real(r8), intent(inout) :: q0_conv(pcols,pver)           ! dq/dt due to moist processes
   real(r8), intent(inout) :: t0_conv(pcols,pver)            ! dt/dt due to moist processes
   real(r8) date0_conv 
#else
   real(r8) :: cnt(pcols)                        ! Top level of convective activity
   real(r8) :: cnb(pcols)                        ! Lowest level of convective activity
   real(r8) :: cmfdqr(pcols,pver)
   real(r8) :: nevapr(pcols,pver)
   real(r8) :: prain(pcols,pver)
   real(r8) :: cmfmc(pcols,pverp)
   real(r8) :: dq_moist(pcols,pver,ppcnst)
   real(r8) :: du_moist(pcols,pver)
   real(r8) :: dv_moist(pcols,pver)
   real(r8) :: dt_moist(pcols,pver)
   real(r8) :: q0_conv(pcols,pver)           ! dq/dt due to moist processes
   real(r8) :: t0_conv(pcols,pver)
   real(r8) date0_conv
#endif
#ifdef MOZART2
   real(r8),intent(inout) :: chemflx(pcols, ppcnst)
   real(r8),intent(inout) :: dq_chem(pcols,pver,ppcnst)
   real(r8),intent(in)    :: flux_ISOP(pcols)
   real(r8),intent(in)    :: flux_ACET(pcols)
   real(r8),intent(in)    :: flux_C3H6(pcols)
   real(r8),intent(in)    :: flux_C2H4(pcols)
   real(r8),intent(in)    :: flux_OC2(pcols)
   real(r8),intent(in)    :: flux_C10H16(pcols)

   real(r8),intent(in)    :: flux_CO2(pcols)
   real(r8),intent(in)    :: flux_N2O(pcols)
   real(r8),intent(in)    :: flux_DMS(pcols)
   real(r8),intent(in)    :: co2mmr_latmean(pcols,pver)
#endif
  
#ifdef CO2
   real(r8),intent(in)    :: co2mmr_latmean(pcols,pver) 
   real(r8)               :: em3_IPCC6(pcols, pver, nem3_IPCC6)
   real(r8)               :: emis_co2_vmr_d3(pcols,pver)
#endif
 
   real(r8)             :: co2vmr_a(pcols)          ! regional distribution of atmospheric CO2 (VMR)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!
   type(dynamics_tend), intent(inout) :: dyn_tend

#if (defined MOZART2)
   real(r8)               :: mbar(pcols,pver)      ! mean wet atmospheric mass ( amu )

   integer                :: n, m
   integer                :: co2_ndx
   integer                :: idx
   real(r8) :: troppre(pcols)
    real(r8) :: flxprec (pcols,pver+1)
    real(r8) :: flxsnow (pcols,pver+1)
    real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
    real(r8) :: gc_frln(pcols)
    real(r8) :: gc_tskin(pcols)
    real(r8) :: gc_sst(pcols)
    real(r8) :: gc_ualb(pcols)
    real(r8) :: gc_cosmi(pcols)
    real(r8) :: gc_lwi(pcols)
    real(r8) :: gc_tropp(pcols)
    real(r8) :: gc_rh(pcols,pver) !relatibe humidity
    real(r8) :: gc_es(pcols,pver) !saturation vapor pressure
    real(r8) :: gc_qs(pcols,pver) !Saturation specific humidity
    real(r8) :: papp1(pcols,pver)!pressure center
    real(r8) :: paphp1(pcols,pver+1) !pressure edge
    real(r8) :: dpp  (pcols,pver) !pressure thick
    real(r8) :: psp1 (pcols)
    real(r8), dimension(:), allocatable :: area_all
     ! First call?
    LOGICAL, SAVE                  :: FIRST = .TRUE.
 
    INTEGER,SAVE                :: GC2BCC(ncnst)
    CHARACTER(128)               :: GCName,wetdepname
    logical                     :: found
    INTEGER                     :: II,JJ,NN,ND
    real(r8)   ::  tauxcl(pcols,0:pver) ! water cloud extinction optical depth
    real(r8)   ::  tauxci(pcols,0:pver) ! ice cloud extinction optical depth
    real(r8)   ::  airmass(pcols,pver)
    real(r8)   ::  zdelt(pcols,pver)
    real(r8)   ::  tmp_data(pcols,pver)
    real(r8)   ::  cldtau(pcols,pver)
 
    real(r8) ::  gcpficu(pcols,pver+1)
    real(r8) ::  gcpflcu(pcols,pver+1)
    real(r8) :: gcconprain(pcols,pver)
    real(r8) :: gcconevapr(pcols,pver)
    real(r8) :: gc_du(pcols,pver)
    real(r8) :: gccmfmc(pcols,pver+1)

     ! sunrz 2024/1
     INTEGER   :: J,             T
     INTEGER   :: typeCounter, maxFracInd(1), sumIUSE
     CHARACTER(LEN=255) :: ErrMsg, ThisLoc
     !sunrz end
    !sunrz 2024/2
    real(r8), intent(inout) :: shf(pcols)          ! Sensible heat flux (w/m^2)
    real(r8), intent(inout) :: lhf(pcols)
    real(r8),intent(in)    :: u10(pcols)               ! zf 2015-05-19
    real(r8),intent(in)    :: v10(pcols)               ! zf 2015-05-19
    real(r8), intent(in) :: taux(pcols)            ! X surface stress (zonal)
    real(r8), intent(in) :: tauy(pcols)            ! Y surface stress(meridional)
    real(r8)                    :: th(pcols,pver)!potential temperature 
    real(r8)                    :: wrk(pcols),rrho(pcols),khfs(pcols),kqfs(pcols),ustar(pcols)
    real(r8)                    :: thvsrf(pcols)
    real(r8)                    :: dryflx(pcols,ppcnst)
     INTEGER                     :: ICOL,ITYPE,MAXTYPE
     real(r8)                    :: maxfrac,thisfrac
     real      :: eflx(pcols,pver,pcnst)
     real      :: lands(pcols,73)
     real      :: lais(pcols,73)
     real      :: area(pcols)
     real      :: area_test(pcols)
     real      :: to3(pcols)
     real      :: airdens(pcols,pver)
     real      :: bxhght(pcols,pver)
     real      :: ad(pcols,pver)
     real      :: olson(pcols)
     real      :: xxlai(pcols)
#endif
#ifdef HAMOCC
   integer              :: dms_idx, n2o_idx
   real(r8)             :: dmsvmr_a(pcols)       !  atmospheric DMS (VMR)
   real(r8)             :: n2ovmr_a(pcols)       !  atmospheric N2O (VMR)
#endif

#if (defined CO2)
   integer                :: idx
#endif
!
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies


   integer  :: nstep                         ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   integer  :: ltrop(pcols)                  ! 
   real(r8) lll(pcols)
   real(r8) CDNC(pcols,pver)         
   real(r8) cldst(pcols,pver)                 ! stratus cloud cover over ocean
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                      !       "     low  cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) prect(pcols)                      ! total (conv+large scale) precip rate
   real(r8) rtdt                              ! 1./ztodt
   real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
!                                             !    maximally overlapped region.
!                                             !    0->pmxrgn(i,1) is range of pressure for
!                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                             !    2nd region, etc
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns
   logical  do_moist,    do_chem

   integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
   integer  i,k                               ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.

   integer :: idox, idco2, idch4, idn2o

!---------------------------------
   real(r8), dimension(pcols,pverp) :: sh_flxprc      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), dimension(pcols,pverp) :: sh_flxsnw      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), dimension(pcols,pver) :: sh_cldliq
   real(r8), dimension(pcols,pver) :: sh_cldice

   real(r8), dimension(pcols,pverp) :: dp_flxprc      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), dimension(pcols,pverp) :: dp_flxsnw      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), dimension(pcols,pver) :: dp_cldliq
   real(r8), dimension(pcols,pver) :: dp_cldice

   real(r8), dimension(pcols,pverp) :: ls_flxprc      ! ls stratiform flux of precip at interfaces (kg/m2/s)
   real(r8), dimension(pcols,pverp) :: ls_flxsnw      ! ls stratiform flux of snow   at interfaces (kg/m2/s)

   real(r8), dimension(pcols,pver) :: ls_reffrain
   real(r8), dimension(pcols,pver) :: ls_reffsnow
   real(r8), dimension(pcols,pver) :: cv_reffliq
   real(r8), dimension(pcols,pver) :: cv_reffice

!----------------------------------------
                                           
   real(r8) rel(pcols,pver)                   ! Liquid cloud particle effective radius
   real(r8) rei(pcols,pver)                   ! Ice effective drop size (microns)
   real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
   real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
!


! physics buffer fields to compute tendencies for cloud condensation package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: qcwat, tcwat, lcwat, cld
   real(r8), pointer, dimension(:,:) :: concld   ! wzz for convective cloud cover
   real(r8), pointer, dimension(:,:) :: cstcld   ! wzz for layered cloud cover
   real(r8), pointer, dimension(:,:) :: pr_cn   ! for convective precipitation for aero-input, wzz, 2008.2.29
   real(r8), pointer, dimension(:,:) :: pr_cs   ! for large-scale precipitation for aero-input, wzz, 2008.2.29


! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
!                                          
! Used for OUTFLD only                     
!                                          
   real(r8), pointer, dimension(:,:,:) :: fracis         ! fraction of transported species that are insoluble
   real(r8) timestep(pcols)
!
!     Variables for doing deep convective transport outside of wu_convr
!
   real(r8) :: mc(pcols,pver)
!
   real(r8) du2(pcols,pver)
   real(r8) conicw(pcols,pver)

! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: flx_heat(pcols)
   logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps
   type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes

   real(r8) co2vmr_obs(pcols,pver)

   real(r8) o3vmr(pcols,pver)
   real(r8) co2vmr(pcols,pver)
   real(r8) ch4vmr(pcols,pver)
   real(r8) n2ovmr(pcols,pver)
   real(r8) cfc11vmr(pcols,pver)
   real(r8) cfc12vmr(pcols,pver)

   real(r8) ch4(pcols,pver)
   real(r8) n2o(pcols,pver)
   real(r8) cfc11(pcols,pver)
   real(r8) cfc12(pcols,pver)
   real(r8) coef, pt, pb

   real(r8) co2mmr_column_mean(pcols)
   real(r8) ptot(pcols)

!   real(r8) :: o3mmr(pcols,pver)
   real(r8) :: co2mmr(pcols,pver)
   real(r8) :: ch4mmr(pcols,pver)
   real(r8) :: n2ommr(pcols,pver)

   real(r8) :: mmr(pcols,pver,220) ! xported species ( mmr )
   real     ::  zmid(pcols,pver)      ! midpoint geopotential in km
   real     :: o3mmr(pcols,pver)   ! sunrz 2025
!--------------------------------------------------------------------
!
   zero = 0.
!
   lchnk = state%lchnk
   ncol  = state%ncol
!  if ( lchnk.eq.249) then
!   write(0,*)'wzz: aero 0 in b:', state%q(2,pver,15)
!  endif


   rtdt = 1./ztodt


   nstep = get_nstep()
   calday = get_curr_calday()


!
! Output NSTEP for debugging
!
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)


!  dry surface pressure
   call outfld ('PSDRY',  state%psdry, pcols, lchnk)


! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QCWAT')
   qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('TCWAT')
   tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('LCWAT')
   lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 
   
   
   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim) 
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('CONCLD')
   concld  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('CSTCLD')
   cstcld  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('CONPRE')
   pr_cn   => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('CSTPRE')
   pr_cs   => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('FRACS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1:ppcnst)  ! wzz, 2008.8.28

!------------------------------------------------
!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.


   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure
!
! Make sure that input tracers are all positive (probably unnecessary)
!
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              ppcnst,qmin  ,state%q )

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

!----------
#ifdef GEOSCHEM
   tend_ciw(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
   tend_clw(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
   tend_w(:ncol,:pver)   = state%q(:ncol,:pver,1)
#endif
!------------

   if ( is_first_step() ) then    ! wzz, 2008.8.29
   fracis (:ncol,:,1:ppcnst) = 1.
   endif

! compute mass integrals of input tracers state
   call check_tracers_init(state, tracerint)


!===================================================
! Global mean total energy fixer
!===================================================
   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
   qini(:ncol,:pver) = state%q(:ncol,:pver,1)


   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )
   call outfld('ENEK',  state%ke,     pcols, lchnk   )
!
!===================================================
! Dry adjustment
!===================================================


! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy


   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)


!#ifdef DRYADJ
!  call dadadj_new (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
!                   ptend%s, ptend%q(1,1,1))
!#else

   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
!#endif
!-------------------------------


   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call physics_update (state, tend, ptend, ztodt)

!===================================================

   do_moist= (nstep.eq.0 .or. nstep.eq.1) .or. (mod(nstep-1, nstep_moist).eq.0 .and. nstep.ne.1)

   if( do_moist ) then

       du_moist(:ncol,:pver)         = state%u(:ncol,:pver)
       dv_moist(:ncol,:pver)         = state%v(:ncol,:pver)
       dt_moist(:ncol,:pver)         = state%t(:ncol,:pver)
       dq_moist(:ncol,:pver,:ppcnst) = state%q(:ncol,:pver,:ppcnst)


       call tphysbc_moist (ztodt,   ts,      sst,     snowh,   landm,   &
                    landfrac,ocnfrac, icefrac, fhs,              &
                    pblht,   kbfs   ,tpert,   qpert,   precl,   &
                    precc,   precsl,  precsc,  prcsnw,  fsds,    &
                    cmfmc,   cldst,   cnt,     cnb,     &
                    dlf,     mc,      prain,   nevapr,  cmfdqr, &
                    du2,     conicw,  ixcldice,ixcldliq,tracerint, &
                    state,   tend,    pbuf,    ptend,   dyn_tend,  &
                    fracis,  qcwat,   tcwat,   lcwat,   cld,       &
                    q0_conv, t0_conv, date0_conv,                 &
                    concld,  cstcld,  pr_cn,   pr_cs , &
                    dp_flxprc,  dp_flxsnw,  dp_cldliq,  dp_cldice,  &
                    sh_flxprc,  sh_flxsnw,  sh_cldliq,  sh_cldice,  &
                    ls_flxprc,  ls_flxsnw,                          &
                    ls_reffrain, ls_reffsnow, cv_reffliq, cv_reffice,icimr,icwmr,gcconprain,gcconevapr, gc_du )

       du_moist(:ncol,:pver)         = (state%u(:ncol,:pver) - du_moist(:ncol,:pver))/ztodt
       dv_moist(:ncol,:pver)         = (state%v(:ncol,:pver) - dv_moist(:ncol,:pver))/ztodt
       dt_moist(:ncol,:pver)         = (state%t(:ncol,:pver) - dt_moist(:ncol,:pver))/ztodt
       dq_moist(:ncol,:pver,:ppcnst) = (state%q(:ncol,:pver,:ppcnst) - dq_moist(:ncol,:pver,:ppcnst)) &
                                             /ztodt

!----------------------------------    
   else
       ptend%lu          = .TRUE.
       ptend%lv          = .TRUE.
       ptend%ls          = .TRUE.
       ptend%lq(:ppcnst) = .TRUE.
       ptend%u(:ncol,:pver)         = du_moist(:ncol,:pver)
       ptend%v(:ncol,:pver)         = dv_moist(:ncol,:pver)
       ptend%s(:ncol,:pver)         = dt_moist(:ncol,:pver) * cpair
       ptend%q(:ncol,:pver,:ppcnst) = dq_moist(:ncol,:pver,:ppcnst)
       call physics_update (state, tend, ptend, ztodt)

   endif

!==========================================================   

   call find_tropopaus( lchnk, ncol, state%zm, state%pmid, state%t, ltrop )

   call diag_dynvar (lchnk, ncol, state)
!====================================================

#ifdef MOZART2
   lll=0
   do i=1, ncol
      troppre(i) =  state%pmid(i,ltrop(i))
      lll(i)     = ltrop(i)
   enddo
   call outfld('TROPP   ',lll ,pcols, lchnk)
     !===============================
     !STEP1:STORE/CALCULATE VARIABLES NEEDED FOR GEOS-CHEM
     !=============================== 
     gc_frln=landfrac !GEOS-Chem FRCLND=FRLAND
     gc_sst=ts
     gc_tskin=ts
     call get_rlat_all_p(lchnk, ncol, clat)
     call get_rlon_all_p(lchnk, ncol, clon)
     call zenith (calday, clat, clon, coszrs, ncol)
     gc_cosmi=coszrs
     gc_ualb=asdir
     !---------------
     !CALCULATE GRID AREA
     !---------------
!#ifdef GEOSCHEM
        IF(FIRST) THEN
!           call gc_calcarea(area_m2_chunk=area_all,lchnk=state%lchnk)
        ! State_Met%AREA_M2(1,:,1)=area_all
        Where(area .LT. 0)
           area=1e-30
        Endwhere
!write(6,*) 'sunrz check area=',size(area_all),size(State_Met%AREA_M2)
        State_Met%AREA_M2(:,1)=area(:ncol)!State_CAM_AREAM2(:ncol,lchnk)!AREA_M2!area_all(:ncol)
! in V14.1.1 area is 2D
        State_Grid%AREA_M2(:,1)=area(:ncol)!State_CAM_AREAM2(:ncol,lchnk)!AREA_M2!area_all(:ncol)
!sunrz 2023/12
        !----------------
        ENDIF
        FIRST=.FALSE.
!#ifdef GEOSCHEM
     DO i=1,pcols

        IF (landfrac(i) > ocnfrac(i).AND. landfrac(i) > icefrac(i)) THEN
        gc_lwi(i) = 1.0 !LAND
        ELSE IF (landfrac(i) > ocnfrac(i)) THEN
        gc_lwi(i) = 2.0 !ICE
        ELSE IF (icefrac(i) > ocnfrac(i)) THEN
        gc_lwi(i) = 2.0 !ICE
        ELSE
        gc_lwi(i) = 0.0 !ocean

        ENDIF
     ENDDO
     !CALCULATE PRESSURE and T-based TROPOPAUSE(gc_tropp)
     !can also use state%pmid
     !---------------------------------
     psp1(:ncol) = state%ps(:ncol)
     call plevs0(ncol, pcols, pver, psp1, paphp1, papp1, dpp)
     do i=1,pcols
     !gc_tropp(i)=papp1(i,ltrop(i)+1)*0.01 !Pa->hPa
     gc_tropp(i)=state%pmid(i,ltrop(i)+1)*0.01 !Pa->hPa
     enddo

     !---------------------------------
     !CALCULATE RELATIVE HUMIDITY
     !---------------------------------
     call aqsat (state%t    ,state%pmid  ,gc_es   ,gc_qs   ,pcols   , &
         ncol ,pver  ,1       ,pver    )
       gc_rh(:ncol,:) = state%q(:ncol,:,1)/gc_qs(:ncol,:)*100.

     call param_cldoptics_calc(state, cld, landfrac, landm,icefrac, &
                                cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh, cldtau)
     call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
     !STEP2: match BCC-GCHP variables meteorology here
     !========================================
     !********************************************
     !**************** 2-D ***********************
     !BC: Before coupling, pass value here
     !AC: After  coupling, pass value tphysac.F90
        ! CONV_DEPTH        Convective cloud depth [m]!sunrz 2024/1
          DO i=1,ncol
          State_Met%CONV_DEPTH(i,1)          = state%zm(i,cnt(i))
          ENDDO

          State_Met%FLASH_DENS(:,1)          = 0*State_Met%CONV_DEPTH(:,1) !no lightning yet
        ! FRSEAICE          Surface sea ice fraction [1] | Sea ice fraction
          State_Met%FRSEAICE(:,1)            = icefrac(:ncol)

        ! Surface albedo, unit:1, AC, but pass value here
          State_Met%ALBD(:,1)                = asdir(:ncol)             ! 1
        !----
        ! Cloud fraction, 1 , BC
          State_Met%CLDFRC(:,1)              = cltot(:ncol)

        ! Land fraction, 1 , BC
          State_Met%FRCLND(:,1)              = landfrac(:ncol)
        !----
        ! Land fraction, 1 , BC
          State_Met%FRLAND(:,1)              = landfrac(:ncol)

        !----
        ! Ocean fraction, 1 , BC
          State_Met%FROCEAN(:,1)             = ocnfrac(:ncol)
        !----
        ! Flags for land type, 1, BC
!sunrz          State_Met%LWI(1,:)                 = gc_LWI! 1

        !----
        !Planetary boundary layer height,m,BC
          State_Met%PBLH(:,1)                = pblht(:ncol)

        !----
        !Convective precipitation, mm day-1, BC
          State_Met%PRECCON(:,1)             = precc(:ncol)* 86400d0*1000! m/s->mm/day
        !Total precipitation,mm day-1, BC : Large scale + Convective
          State_Met%PRECTOT(:,1)             = precl(:ncol)* 86400d0*1000+precc(:ncol)*86400d0*1000   ! m/s->mm/day
        !--------
        !fraction of convection,xiaolu, 2018/09
          DO i=1,ncol
             IF (State_Met%PRECTOT(i,1) .GT. 0.0) THEN
             State_Met%CNV_FRC(i,1)  = State_Met%PRECCON(i,1)/State_Met%PRECTOT(i,1)
             ENDIF
          ENDDO
        !----
        !Tropopause, hPa, BC
          State_Met%TROPP(:,1)               = gc_tropp(:ncol)! hPa

          State_Met%PS1_WET(:,1)             = state%ps(:ncol)*0.01 !hPa

        !----
        !Surface temperature, K ,BC
          State_Met%TS(:,1)                  = ts(:ncol) ! K

        !Surface temperature, K ,BC
          State_Met%TSKIN(:,1)               = gc_tskin(:ncol) ! K

        !----
        !Sea surface temperature, K,
!sunrz          State_Met%SST(1,:)                 = gc_sst! K

        !----
        !cosine of zenith angle,1, BC
          State_Met%SUNCOS(:,1)              = gc_cosmi(:ncol) ! unitless

        !----
        !cosine of zenith angle,1,BC
          State_Met%SUNCOSmid(:,1)           = gc_cosmi(:ncol) ! 1
     !****************************************
     !************3-D (pcol, pver+1) *********
     !****************************************
        !----
        !Upward moist convective mass flux, kg m-2 s-2, not use
        !xiaolu,warning:current unit is ! kg/m2/s
        !devide the time step: ztodt
        !Unit should be kg/m2/s
          State_Met%CMFMC(:,1,:)   = cmfmc(:ncol,pver+1:1:-1)        ! kg/m2/s

!        Where(State_Met%CMFMC .LT. 0)
!           State_Met%CMFMC=0
!        Endwhere
!        Where(State_Met%CMFMC .GT. 0.2)
!           State_Met%CMFMC=0.2
!        Endwhere


          State_Met%PFILSAN(:,1,:) = ls_flxsnw(:ncol,pver+1:1:-1)+sh_flxsnw(:ncol,pver+1:1:-1)  ! kg/m2/s
        !----
        !xiaolu,2018/12/30
        !3d flux of ice convective precipitation,kg m-2 s-1, BC
        gcpficu=dp_flxsnw
        State_Met%PFICU(:,1,:) = dp_flxsnw(:ncol,pver+1:1:-1) ! kg/m2/s 

        !----
        !3d flux of liquid non-convective precipitation,kg m-2 s-1, BC
          State_Met%PFLLSAN(:,1,:) = sh_flxprc(:ncol,pver+1:1:-1)+ls_flxprc(:ncol,pver+1:1:-1)-State_Met%PFILSAN(:,1,:)!kg/m2/s

        State_Met%PFLCU(:,1,:) = dp_flxprc(:ncol,pver+1:1:-1)-State_Met%PFICU(:,1,:) ! kg/m2/s

        !Pressure edge for each grid, NOTE: PLE is reversed in the vertical!
          State_Met%PEDGE(:,1,:)  = state%pint(:ncol,pver+1:1:-1)*0.01 ! Pa -> hPa
          State_Met%PMID(:,1,:)   = state%pmid(:ncol,pver:1:-1)*0.01
          State_Met%CLDF(:,1,:)                = cld(:ncol,pver:1:-1)
        ! Detrainment cloud mass flux kg/m2/s 
        ! warning: present unit is 1/mb
        ! Set detrainment to 0, following wutw's suggestion,2019/01/14
           !State_Met%DTRAIN(:,1,:)=du2(:ncol,pver:1:-1)*100/9.8
           State_Met%DTRAIN(:,1,:)=gc_du(:ncol,pver:1:-1)
          !write(*,*)'xiaolu check Dtrain',minval(du2),maxval(du2)

        !----
        ! Convective rainwater source kg kg-1 s-1
          Where (gcconprain .LT. 0)
            gcconprain=0.0
          ENDWHERE

          State_Met%DQRCU(:,1,:)=gcconprain(:ncol,pver:1:-1)
        !  write(*,*)'xiaolu check DQRCU',minval(gcconprain),maxval(gcconprain)

        !----
        !Large-scale rainwater source,kg kg-1 s-1,BC
          airmass(:,:) = state%pmid(:,:)/(287.04*state%t(:,:))   ! kg/m3
          State_Met%DQRLSAN(:,1,:)             = prain(:ncol,pver:1:-1)

        !----
        !Cloud ice water mixing ratio
          State_Met%QI(:,1,:)                  = icimr(:ncol,pver:1:-1) ! kg/kg
          State_Met%QL(:,1,:)                  = icwmr(:ncol,pver:1:-1) ! kg/kg
        !----
        !Relative humidity
          State_Met%RH(:,1,:)                  = gc_rh(:ncol,pver:1:-1)!1 -> %

        !----
        !Specific humidity
          State_Met%SPHU(:,1,:)                = state%q(:ncol,pver:1:-1,1) * 1d3 ! kg/kg ->g/kg

        !----
        !Temperature
          State_Met%T(:,1,:)                   = state%t(:ncol,pver:1:-1) ! K

          State_Met%REEVAPCN(:,1,:)            =   gcconevapr(:ncol,pver:1:-1) ! kg/kg/s
          State_Met%REEVAPLS(:,1,:)            = nevapr(:ncol,pver:1:-1) ! kg/kg/s

           State_Met%OPTD(:,1,:)             = cldtau(:ncol,pver:1:-1)!*cld(:,pver:1:-1)

        !------------------------------------
        !xiaolu: besides listed Include_GCHP.H,the following vars also from GCHP
        !------------------------------------
        !----2D
        State_Met%PS1_DRY(:,1)=state%psdry(:ncol)*0.01

        !-----3D
        State_Met%DELP_DRY(:,1,:) = state%pdeldry(:ncol,pver:1:-1)*0.01
         do i=1,ncol
            rrho(i)   = rair*state%t(i,pver)/state%pmid(i,pver)
            ustar(i)  = max(sqrt(sqrt(taux(i)**2 + tauy(i)**2)*rrho(i)),0.01)
         !   khfs(i)   = shf(i)*rrho(i)/cpair
         enddo
         !  do i=1,ncol
         !     kqfs(i)= cflx(i,1)*rrho(i)
         !  end do
           th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)
        ! Compute Obukhov length virtual temperature flux and various arrays for
        ! use later:
       ! do i=1,ncol
       !     kbfs(i)      = khfs(i) + 0.61*th(i,pver)*kqfs(i)
       !     thvsrf(i)    = th(i,pver)*(1.0 + 0.61*state%q(i,pver,1))
       !     obklen(i)    = -thvsrf(i)*ustar(i)**3/(gravit*karman*(kbfs(i)
       !     +sign(1.e-10_r8,kbfs(i))))
       ! end do
        ! Sensible heat flux from turbulence, w m-2, AC ,use for calculate
        ! obklen,so not use here as obklen is passed directly from bcc
          State_Met%EFLUX(:,1)               = lhf(:ncol)
          State_Met%HFLUX(:,1)               = shf(:ncol)

          State_Met%U10M(:,1)                = u10(:ncol) ! m/s

        !----
        !Friction velocity,m s-1 AC
          State_Met%USTAR(:,1)               = ustar(:ncol)!*landfrac(:ncol)+0.02*(1-landfrac(:ncol)) ! m/s

        !----
        !10m V-wind, m s-1, AC
          State_Met%V10M(:,1)                =  v10(:ncol) ! m/s

        !----
        !Roughness length , m , AC, read from landtype
          State_Met%Z0(:,1)                  = 0.1                      !sunrz readd m  
        State_Met%PHIS(:,1)=state%phis(:ncol)
        State_Met%PARDR(:,1)=sols(:ncol)
        State_Met%PARDF(:,1)=solsd(:ncol)

        State_Met%SWGDN(:,1)=fsds(:ncol)
        State_Met%FRLAKE(:,1)=0
        State_Met%FRLANDIC(:,1)=0
        State_Met%FRSNO(:,1)=0
        State_Met%GWETROOT(:,1)=fhs(:ncol)
        State_Met%GWETTOP(:,1)=fhs(:ncol)
        State_Met%PRECANV(:,1)=0
        State_Met%PRECLSC(:,1)=precl(:ncol)* 86400d0*1000
        State_Met%PS2_WET(:,1)=State_Met%PS1_WET(:,1)
        State_Met%PS2_DRY(:,1)=State_Met%PS1_DRY(:,1)
        State_Met%SLP(:,1)= State_Met%PS1_WET(:,1)
        State_Met%SNODP(:,1)=snowh(:ncol)
        State_Met%SNOMAS(:,1)=snowh(:ncol)*916.7
        State_Met%UVALBEDO(:,1) = asdir(:ncol)
        DO i=1,ncol
        State_Met%Z0(i,1) = State_Met%FRLAND(i,1) * 0.01  &
             + State_Met%FRSEAICE(i,1)  * 0.04 &
             + State_Met%FROCEAN(i,1)  * 0.00001
        ENDDO
        !In-cloud ice cloud optical thickness (visible band),BC,no need
        State_Met%TAUCLI(:,1,:)              = tauxci(:ncol,pver:1:-1)     ! 1

        !----
        !In-cloud liquid cloud optical thickness (visible band),BC
        State_Met%TAUCLW(:,1,:)              = tauxcl(:ncol,pver:1:-1)     ! 1

        State_Met%CLDTOPS(:,1)=0 !not yet for lightning
        State_Met%U(:,1,:)=state%u(:ncol,pver:1:-1)
        State_Met%V(:,1,:)=state%v(:ncol,pver:1:-1)
        State_Met%OMEGA(:,1,:)=state%omega(:ncol,pver:1:-1)
    do i=1,ncol
        State_Met%TO3(i,1) = 0.0
        do k=1,pver
        State_Met%TO3(i,1) = State_Met%TO3(i,1)&
                   + State_Chm%Species(180)%Conc(i,1,k) &
                      * State_Met%AIRDEN(i,1,k)   &
                      * State_Met%BXHEIGHT(i,1,k)
        enddo
        State_Met%TO3(i,1) = State_Met%TO3(i,1) * ( 6.02e23 / 48 ) / 1e+1 / 2.69e+16
        to3(i)=State_Met%TO3(i,1)
        airdens(i,:)=State_Met%AIRDEN(i,1,pver:1:-1)
        bxhght(i,:)=State_Met%BXHEIGHT(i,1,pver:1:-1)
        ad(i,:)=State_Met%AD(i,1,pver:1:-1)
        area_test(i)=State_Met%AREA_M2(i,1)
    enddo

!#endif
 !=====================================
 !STEP 3: match chemical speices
 !=====================================
!     DO II=1,220
!       State_Chm%Species(II)%Conc(:,1,:) = state%q(:ncol,pver:1:-1,II+4)
!       State_Chm%Species(II)%Conc(:,1,:)=State_Chm%Species(II)%Conc(:,1,:)+eflx(:ncol,pver:1:-1,II+4)/100*9.80665/State_Met%DELP_DRY(:,1,:)*dtime
!     ENDDO
        State_Met%LandTypeFrac(:,1,:) = lands(:ncol,:)
        State_Met%XLAI_NATIVE(:,1,:)  = lais(:ncol,:)


DO II=1,ncol
olson(II)=0
xxlai(II)=0
DO i=1,73
olson(II)=olson(II)+lands(II,i)*i
xxlai(II)=xxlai(II)+lais(II,i)*1
ENDDO
ENDDO

      do i = 1,pver
         zmid(:,i) = 1.e-3 * state%zm(:,i)
      enddo

      
      DO i=1,220
      mmr(:,:,i)=state%q(:,:,i+3)
      ENDDO
      !o3vmr(:,:)=mmr(:,:,180)*48/28.966
      call set_ub_vals( ztodt, lchnk, mmr, state%pmid, zmid, &
                        state%t, calday, pcols , ltrop, o3vmr )

      !write(*,*) 'sunrz check ori=', state%q(:,1:2,183)
      !write(*,*) 'sunrz check ub=', o3vmr(:,1:2)
      state%q(:,:,183)=mmr(:,:,180)
      
       call tphysbc_geoschem (ztodt,state,calday,State_Met,&
            State_Chm,Input_Opt,State_Grid,State_Diag,ts,landfrac,ocnfrac,icefrac,&
            flux_ISOP, flux_ACET,  flux_C3H6,                &
            flux_C2H4, flux_OC2,   flux_C10H16,              &
            flux_CO2,  flux_N2O,   flux_DMS,                 &
            chemflx,cnt,cnb,                                 &
            GC2BCC,eflx)

    call outfld('T_3D', state%t, pcols, lchnk)
    call outfld('P_3D', state%pmid, pcols, lchnk)
    call outfld('SUNCOS',gc_cosmi, pcols, lchnk)
    call outfld('TO3',to3,pcols,lchnk)
    call outfld('AREA', area_test,pcols,lchnk)
    call outfld('REEVAPLS', nevapr,pcols,lchnk)
    call outfld('DQRLSAN', prain,pcols,lchnk)
    call outfld('RH', gc_rh,pcols,lchnk)
    call outfld('PFILSAN', ls_flxsnw+sh_flxsnw+dp_flxsnw,pcols,lchnk)
    call outfld('PFLLSAN', ls_flxprc+sh_flxprc+dp_flxprc-ls_flxsnw-sh_flxsnw-dp_flxsnw,pcols,lchnk)
    call outfld('DELP',state%pdeldry,pcols,lchnk)
    call outfld('PEDGE',state%pint,pcols,lchnk)
    call outfld('ncol',area_test*0+ncol,pcols,lchnk)
    call outfld('lchnk',area_test*0+lchnk,pcols,lchnk)
    call outfld('olson',olson,pcols,lchnk)
    call outfld('eflx',eflx(:,pver:1:-1,178),pcols,lchnk)
    call outfld('eflx2',eflx(:,pver:1:-1,179),pcols,lchnk)
    call outfld('USTAR2',ustar,pcols,lchnk)
    call outfld('Z0',state_met%z0(:,1),pcols,lchnk)
    call outfld('XLAI',xxlai,pcols,lchnk)
    call outfld('DTRAIN',gc_du,pcols,lchnk)
    call outfld('DQRCU',gcconprain,pcols,lchnk)
    call outfld('REEVAPCN',gcconevapr,pcols,lchnk)
    call outfld('GCCMFMC',cmfmc,pcols,lchnk)
    call outfld('HCO_178',eflx(:,pver:1:-1,178),pcols,lchnk)
    call outfld('HCO_59',eflx(:,pver:1:-1,59),pcols,lchnk)
    call outfld('HCO_137',eflx(:,pver:1:-1,137),pcols,lchnk)
    call outfld('HCO_209',eflx(:,pver:1:-1,209),pcols,lchnk)
    call outfld('HCO_212',eflx(:,pver:1:-1,212),pcols,lchnk)
    call outfld('HCO_62',eflx(:,pver:1:-1,62),pcols,lchnk)
    call outfld('HCO_63',eflx(:,pver:1:-1,63),pcols,lchnk)
    call outfld('HCO_64',eflx(:,pver:1:-1,64),pcols,lchnk)
    call outfld('HCO_65',eflx(:,pver:1:-1,65),pcols,lchnk)
    call outfld('HCO_14',eflx(:,pver:1:-1,14),pcols,lchnk)
    call outfld('HCO_15',eflx(:,pver:1:-1,15),pcols,lchnk)
    call outfld('HCO_185',eflx(:,pver:1:-1,185),pcols,lchnk)
    call outfld('HCO_186',eflx(:,pver:1:-1,186),pcols,lchnk)
#endif
 
#ifdef CO2

   em3_IPCC6(:,:,:) = 0.0_r8
   call cnst_get_ind('CO2', idco2)

  call getem3_IPCC6( lchnk, ncol, state%zi, em3_IPCC6 )

  call co2_em3_prod ( ncol, ztodt, state%q(1,1,idco2),  state%pmid, state%t,   &
                       em3_IPCC6, nem3_IPCC6, emis_co2_vmr_d3)
#endif

!===================================================
! Radiation computations
!===================================================
!
! Cosine solar zenith angle for current time step
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call zenith (calday, clat, clon, coszrs, ncol)


!====================================================================================
#ifdef MOZART2
!   call cnst_get_ind('O3', idox)
!   do i=1, ncol
!      do k=1, pver
!         o3vmr(i,k) = state%q(i,k,idox) * 28.966 / 47.998  ! mmr==> vmr
!      enddo 
!   enddo
   call radozn(lchnk   ,ncol    ,state%pmid    ,o3vmr   )
#else
   call radozn(lchnk   ,ncol    ,state%pmid    ,o3vmr   )
#endif
   call outfld('O3VMR   ',o3vmr ,pcols, lchnk)

!--------------
#ifdef MOZART2

#ifdef GHG_2D
   call ghg_surfvals_get_co2vmr_2d(lchnk, ncol, co2vmr_obs)
#else
   co2vmr_obs(:,:) = chem_surfvals_get('CO2VMR')
#endif

!   call cnst_get_ind('CO2', idco2)

!   do i=1, ncol
!      co2mmr_column_mean(i) = 0.0
!      ptot(i)               = 0.0
!      do k=1, pver
!         co2mmr_column_mean(i) = co2mmr_column_mean(i) + state%q(i,k,idco2) * state%pdeldry(i,k)
!         ptot(i)               = ptot(i) + state%pdeldry(i,k)
!      enddo
!      co2mmr_column_mean(i) = co2mmr_column_mean(i) / ptot(i)
!   enddo

!   do i=1, ncol
!      do k=1, pver
!          co2vmr(i,k) = co2vmr_obs(i,k)
!          state%q(i,k,idco2) = co2vmr(i,k) * 44.0 / 28.966
 !        co2vmr(i,k) = co2mmr_latmean(i,k) * 28.966 / 44.0    ! mmr = > vmr [ppm]
 !        co2vmr(i,k) = state%q(i,k,idco2)  * 28.966 / 44.0    ! mmr = > vmr [ppm]
 !        co2vmr(i,k) =  co2mmr_column_mean(i) * 28.966 / 44.0
 !     enddo
 !  enddo

!------
#elif (defined CO2)

#ifdef GHG_2D
   call ghg_surfvals_get_co2vmr_2d(lchnk, ncol, co2vmr_obs)
#else
   co2vmr_obs(:,:) = chem_surfvals_get('CO2VMR')
#endif

   call cnst_get_ind('CO2', idco2)

   do i=1, ncol
      co2mmr_column_mean(i) = 0.0
      ptot(i)               = 0.0 
      do k=1, pver
         co2mmr_column_mean(i) = co2mmr_column_mean(i) + state%q(i,k,idco2) * state%pdeldry(i,k)
         ptot(i)               = ptot(i) + state%pdeldry(i,k)
      enddo
      co2mmr_column_mean(i) = co2mmr_column_mean(i) / ptot(i)
   enddo

   do i=1, ncol
      do k=1, pver
 !        co2vmr(i,k) = co2vmr_obs(i,k) 
 !        co2vmr(i,k) = co2mmr_latmean(i,k) * 28.966 / 44.0    ! mmr = > vmr [ppm]
 !        co2vmr(i,k) = state%q(i,k,idco2)  * 28.966 / 44.0    ! mmr = > vmr [ppm]
          co2vmr(i,k) =  co2mmr_column_mean(i) * 28.966 / 44.0  
      enddo
   enddo
  
!----
#else

#ifdef GHG_2D
   call ghg_surfvals_get_co2vmr_2d(lchnk, ncol, co2vmr)
#else
   co2vmr(:,:) = chem_surfvals_get('CO2VMR')  ! wzz, for carbon cycle, 2008.11.20
#endif

#endif
!---

!------------------------------------------------

! get mass mixing ratios of CH4, N2O, CFC11
!
   call trcmix(lchnk         ,ncol    , &
               state%pmid    ,n2o     ,ch4     ,                     &
               cfc11         ,cfc12   )
   do i=1, ncol
      do k=1, pver
         n2ovmr(i,k) = n2o(i,k)/rmwn2o      ! mmr->vmr
         ch4vmr(i,k) = ch4(i,k)/rmwch4      ! mmr->vmr
         cfc11vmr(i,k) = cfc11(i,k)/rmwf11  ! mmr->vmr
         cfc12vmr(i,k) = cfc12(i,k)/rmwf12  ! mmr->vmr
      enddo
   enddo

#ifdef MOZART2
!--- update mixing ratios of CH4, N2O

   call  cnst_get_ind('CH4', idch4)
   call  cnst_get_ind('N2O', idn2o)

   do i=1, ncol
      do k=1, pver
         state%q(i,k,idn2o) = n2ovmr(i,k)*rmwn2o
         state%q(i,k,idch4) = ch4vmr(i,k)*rmwch4
      enddo
   enddo
#endif

    call outfld('CO2VMR', co2vmr, pcols, lchnk)
    call outfld('CH4VMR', ch4vmr, pcols, lchnk)
    call outfld('N2OVMR', n2ovmr, pcols, lchnk)
    call outfld('CFC11', cfc11vmr, pcols, lchnk)
    call outfld('CFC12', cfc12vmr, pcols, lchnk)
!===============================================================
   if (dosw .or. dolw) then

! Compute cloud water/ice paths and optical properties for input to radiation
      call t_startf('cldoptics')

      call param_cldoptics_calc(state, cld, landfrac, landm,icefrac, &
                                cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh, cldtau)
      call t_stopf('cldoptics')

!--------------------------------------------------------------------------------------
! Complete radiation calculations
!
#ifdef BCC_RAD
!
! zhanghua scheme
!
        call t_startf('radctl_zh')
        call radctl_zh (lchnk, ncol, lwup, coszrs, calday, &
             clat, clon, cld, rel, rei, state, &
             aldir, aldif, asdir, asdif, ocnfrac, landfrac, icefrac, &
             qrs, qrl, flwds, fsns, fsnt, flns, flnt, sols, soll, solsd, solld, fsds )

        call t_stopf ('radctl_zh')

       ! qrs(:,:)=0.0_r8
       ! qrl(:,:)=0.0_r8
       ! flwds(:)=0.0_r8
       ! fsns(:)=0.0_r8
       ! fsnt(:)=0.0_r8
       ! flns(:)=0.0_r8
       ! flnt(:)=0.0_r8
       ! sols(:)=0.0_r8
       ! soll(:)=0.0_r8
       ! solsd(:)=0.0_r8
       ! solld(:)=0.0_r8
       ! fsds(:)=0.0_r8

#else
      call t_startf ('radctl')

#ifdef GEOSCHEM
      call radctl (lchnk, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cld, cicewp, cliqwp, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, state%zm, state, fsds, tauoptd , &
                   o3vmr,  co2vmr,  ch4vmr, n2ovmr, cfc11vmr, cfc12vmr )
#else
      call radctl (lchnk, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cld, cicewp, cliqwp, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, state%zm, state, fsds, &
                   o3vmr,  co2vmr,  ch4vmr, n2ovmr, cfc11vmr, cfc12vmr )
#endif
      call t_stopf ('radctl')
#endif
!------------------------------------------------------------------
!
! Cloud cover diagnostics
! radctl can change pmxrgn and nmxrgn so cldsav needs to follow 
! radctl.
!
      call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
!
! Dump cloud field information to history tape buffer (diagnostics)
!
      call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
      call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
      call outfld('CLDMED  ',clmed  ,pcols,lchnk)
      call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
      call outfld('CLOUD   ',cld    ,pcols,lchnk)

#ifdef DO_COSP

     co2mmr(:,:) = co2vmr(:,:) *44./28.966
     o3mmr(:,:)  = o3vmr(:,:)  *47.998/28.966
     ch4mmr(:,:) = ch4vmr(:,:) /rmwch4
     n2ommr(:,:) = n2ovmr(:,:) /rmwn2o

     ! advance counter for this timestep (chunk dimension required for thread safety)
       cosp_cnt = cosp_cnt + 1

     ! if counter is the same as cosp_nradsteps, run cosp and reset counter

       call cospsimulator_intr_run( lchnk, ncol, state%q(1,1,1), state%pint, state%zi, &
                       state%phis, state%zm, state%lat, state%lon, state%t,              &
                       state%pmid, state%q(1,1,ixcldliq), state%q(1,1,ixcldice),         &
                       state%ps,   state%u,               state%v,                       &
                       landfrac,  coszrs, o3mmr, co2mmr, ch4mmr, n2ommr, ts,     &
                       cld,       concld, rel, rei, emis, cliqwp, cicewp,                &
                       dp_flxprc, dp_flxsnw, sh_flxprc, sh_flxsnw, ls_flxprc, ls_flxsnw, &
                       dp_cldliq, dp_cldice, sh_cldliq, sh_cldice,                       &
                       ls_reffrain, ls_reffsnow, cv_reffliq, cv_reffice )
                    ! optional      cld_swtau_in, snow_tau_in, snow_emis_in )

       cosp_cnt = 0 

#endif

   else


! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k =1 , pver
            do i = 1, ncol
               qrs(i,k) = qrs(i,k)/state%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state%pdel(i,k)
            end do
         end do
      end if
         
   end if

!
! Compute net flux
! Since fsns, fsnt, flns, and flnt are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   do i=1,ncol
      tend%flx_net(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do
!
! Compute net radiative heating
!
   call radheat_net (state, ptend, qrl, qrs)
!
! Add radiation tendencies to cummulative model tendencies and update profiles
!
   call physics_update(state, tend, ptend, ztodt)
!  if ( lchnk.eq.249) then
!   write(0,*)'wzz: aero 31 in b:', state%q(2,pver,15)
!  endif


! check energy integrals
   call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, tend%flx_net)
!
! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   srfrad(:ncol) = fsns(:ncol) + flwds(:ncol)
   call outfld('SRFRAD  ',srfrad,pcols,lchnk)
!
! Save atmospheric fields to force surface models
!
#if (defined MOZART2) 
!   call cnst_get_ind('CO2', idx)
   
!   do i=1, ncol
!      co2vmr_a(i) =  state%q(i,pver,idx) * 28.9644 / 44.0
	              ! kg/kg to m3/m3,  mean molecular weight of dry air: 28.9644 (g/mol)2
!   enddo

#ifdef HAMOCC
   call cnst_get_ind('DMS', dms_idx)
   call cnst_get_ind('N2O', n2o_idx)
   do i=1, ncol
      dmsvmr_a(i) =  state%q(i,pver,dms_idx) * 28.9644 / 62.132
                      ! kg/kg to m3/m3,  mean molecular weight of dry air: 28.9644 (g/mol)2
      n2ovmr_a(i) =  state%q(i,pver,n2o_idx) * 28.9644 / 44.0128
                      ! kg/kg to m3/m3,  mean molecular weight of dry air: 28.9644 (g/mol)2
   enddo
#endif

#elif (defined CO2)
   call cnst_get_ind('CO2', idx)
   do i=1, ncol
      co2vmr_a(i) =  state%q(i,pver,idx) * 28.9644 / 44.0
	                              ! kg/kg to m3/m3,  mean molecular weight of dry air: 28.9644 (g/mol)2
   enddo

#ifdef HAMOCC
   do i=1, ncol
      dmsvmr_a(i) =  0.0
      n2ovmr_a(i) =  n2ovmr(i,pver)
   enddo
#endif

#else
!---------------------
#ifdef GHG_2D
   call ghg_surfvals_get_co2vmr_2d(lchnk, ncol, co2vmr)
   do i=1, ncol
      co2vmr_a(i) = co2vmr(i,pver)   ! global mean CO2 (vmr)
   enddo
#else
   co2vmr_a(:) = chem_surfvals_get('CO2VMR')
#endif
!--------------------
#endif

   call srfxfer (lchnk, ncol, state%ps, state%u(1,pver), state%v(1,pver),    &
                 state%t(1,pver), & 
                 co2vmr_a,        &
		 state%q(1,pver,1), state%exner(1,pver), state%zm(1,pver), &
                 state%pmid, state%rpdel(1,pver))


!-------------------------------------------------------------------
! start to call GEOSCHEM
!-------------------------------------------------------------------
#ifdef GEOSCHEM


!   call acbndget( 'H2O2', ncol, lchnk, h2o23d )
!   call acbndget( 'O3'  , ncol, lchnk, o3 )
!   call acbndget( 'HO2' , ncol, lchnk, ho2im )        ! 1.52e8
!   call acbndget( 'NO3' , ncol, lchnk, no3   )        ! 1.64e7


!><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   ! Add emissions to tracer tendency array
   ! Emissions MUST BE IN CORRECT UNITS!!!

!
!  Get HO2: concentration    units: molec/cm3
!
   convert(:ncol,:) = (1.d0/6.022d23)*8.314d0*(state%T(:ncol,:)/state%pmid(:ncol,:))*1.d6 &
                      *33.001d0/28.98d0
   ! 1.d6 => cm3/m3; 17.001 = MW(OH); 28.98 = MW(dry air); 8.314 = Universal Gas Constant
   call cnst_get_ind('HO2', idx)
   molar_mass(idx) = 33.01
   call acbndget( 'HO2'  , ncol, lchnk, oh(:ncol,:))  !  To directly read from prescribed data file
   state%q(:ncol,1:pver,idx) = oh(:ncol,1:pver)*convert(:ncol,:) ! Converted to kg/kg
   call outfld('HO2',  state%q(:,:,idx)  ,pcols   ,lchnk   ) ! Output here before chemistry
!
!
!  Get OH: concentration    units: molec/cm3
!
   convert(:ncol,:) = (1.d0/6.022d23)*8.314d0*(state%T(:ncol,:)/state%pmid(:ncol,:))*1.d6 &
                      *17.001d0/28.98d0
   ! 1.d6 => cm3/m3; 17.001 = MW(OH); 28.98 = MW(dry air); 8.314 = Universal Gas Constant
   call cnst_get_ind('OH', idx)
   molar_mass(idx) = 17.01
   call acbndget( 'OH'  , ncol, lchnk, oh(:ncol,:))  !  To directly read from prescribed data file
   state%q(:ncol,1:pver,idx) = oh(:ncol,1:pver)*convert(:ncol,:) ! Converted to kg/kg
   call outfld('OH',  state%q(:,:,idx)  ,pcols   ,lchnk   ) ! Output here before chemistry
!
!
! Get surface emission (em3data and em2data)
!
   call  getem3_IPCC( lchnk, ncol, state%zm, em3data )
   call  getem2_IPCC( lchnk, ncol, em2data )
!
!  Get CO  & convert units
!
   call cnst_get_ind('CO', idx)


   molar_mass(idx)  = 28.0103988600000    ! CO
   srfflx(:ncol)  = amufac * molar_mass(idx) *  em2data(:ncol,idxCO)
                            ! units: molecules/cm2/s --> surface emission flux ( kg/m^2/s )
   ptend%q(:ncol,pver,idx) = srfflx(:ncol)*gravit*state%rpdel(:ncol,pver)
                            ! --> kg/kg/s
   ptend%lq(idx)           = .TRUE.
!  Get CO  & convert units
!
!   call cnst_get_ind('NOx', idx)
!
!   molar_mass(idx)  = 46.0    ! NO & NO2
!   srfflx(:ncol)  = amufac * molar_mass(idx) *  em2data(:ncol,idxNO)
!                            ! units: molecules/cm2/s --> surface emission flux ( kg/m^2/s )
!   ptend%q(:ncol,pver,idx) = srfflx(:ncol)*gravit*state%rpdel(:ncol,pver)
!                            ! --> kg/kg/s
!   ptend%lq(idx)           = .TRUE.
!
!  Get CH4  & convert units
!
   call cnst_get_ind('CH4', idx)


   molar_mass(idx)  = 16.0405998200000     ! CH4
   srfflx(:ncol)  = amufac * molar_mass(idx) *  em2data(:ncol,idxCH4)
                             ! units: molecules/cm2/s --> surface emission flux ( kg/m^2/s )
   ptend%q(:ncol,pver,idx) = srfflx(:ncol)*gravit*state%rpdel(:ncol,pver)


   ptend%lq(idx)           = .TRUE.
!  Update Physics Tendency Arrays
!


   !QUICK FIX
   call cnst_get_ind('CO2', idx)
   molar_mass(idx) = 44.00


   call cnst_get_ind('CH2O', idx)
   molar_mass(idx) = 30.02
   !---------


   call physics_update(state, tend, ptend, ztodt)


!><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   do i = 1, ncol
      if ( fsds(i) .lt. 1.0e-3) then
         albedo(i) = 1.
      else
         albedo(i) =  1.0 - (fsns(i)/fsds(i))
      end if
      uv_albedo(i) = albedo(i)
   end do


   ! pressure of tropopause (Pa)
   do i = 1, ncol
      coslat(i) = cos(clat(i))
      ptrop(i) = 250.0e2 - 150.0e2*coslat(i)**2.0
   end do


   do i = 1, ncol
      if (landfrac(i) .gt. 0.5) then
        z0(i) = 5.e-2
        ustar(i) = 1.0
      else if (ts(i) .lt. (273.16 - 1.8) .and. landfrac(i) .lt. 0.1) then
        z0(i) = 5.e-4
        ustar(i) = 0.2
      else
        z0(i) = 5.e-3
        ustar(i) = 0.7
      end if


      u_10m(i) = state%u(1,pver)*log(10.0/z0(i))/log((state%zi(i,pver)-state%zi(i,pver+1))*0.5/z0(i))
      v_10m(i) = state%v(1,pver)*log(10.0/z0(i))/log((state%zi(i,pver)-state%zi(i,pver+1))*0.5/z0(i))
      gwettop(i) = 0.30
   end do


   do i = 1, ncol
      do k = 1, pver
        rho(i,k) = state%pmid(i,k) / (rair * state%t(i,k))
        dz(i,k)  = state%zi(i,k)-state%zi(i,k+1)


        tend_ciw(i,k) = (state%q(i,k,ixcldice)-tend_ciw(i,k))/ztodt
        tend_clw(i,k) = (state%q(i,k,ixcldliq)-tend_clw(i,k))/ztodt
        tend_w(i,k) = (state%q(i,k,1)-tend_w(i,k))/ztodt
      end do
   end do


   call aqsat(state%t(1,1)   , state%pmid    ,es      ,qs      ,pcols   , &
               ncol    ,pver    ,1       ,pver    )
   do i=1, ncol
      do k=1, pver
         relhum(i,k) = qs(i,k)/qs(i,k)
      enddo
   enddo


   !!!call cnst_get_ind('CO',coidx) ! CAM index yields 4
   !!!coidx=coidx-3    ! 4-->1 (CAM --> MOZART2)


!   call cnst_get_ind('CH4',idx)
!   write(*,*) 'PRE PHYS CH4 SUM : ',  (state%q(1,26,idx))


   call do_gc_chem( state%q(:ncol,:,ixchm:) , state%q(:,:,1) , state%t  ,     &
                    state%pmid,          state%pdel, state%pint,                   &
                    albedo      , cltot   , landfrac  , pblht          ,     &
                    precl+precc , precc   , ptrop     , trefht         ,     &
                    sst         , u_10m   , v_10m     , ustar          ,     &
                    z0          , cld     , mc        , du2            ,     &
                    cicewp      , cliqwp  , relhum    , shflx          ,     &
                    rho         , coslat  , tauoptd(:,1:pver)          ,     &
                    sols        , solsd   , fsds      , uv_albedo      ,     &
                    dz          , tend_ciw, tend_clw  , tend_w         ,     &
                    tend_dynq + tend_w,     ts        , gwettop        ,     &
                    abs(state%lat(2)-state%lat(1)), abs(state%lon(2)-state%lon(1)), &
                    state%lat, state%lon, state%ncol, state%zm, lchnk)


!    call cnst_get_ind('CO2',idx)
!    write(*,*) 'PHYS CO2 SUM : ',  (state%q(:ncol,:,idx))
!    call cnst_get_ind('CO',idx)
!    write(*,*) 'PHYS CO SUM : ',   (state%q(:ncol,:,idx))
    call cnst_get_ind('CH4',idx)
!    write(*,*) 'PHYS CH4 SUM : ',  (state%q(:ncol,:,idx))
!    call cnst_get_ind('CH2O',idx)
!    write(*,*) 'PHYS CH2O SUM : ', (state%q(:ncol,:,idx))
!    call cnst_get_ind('OH',idx)
!    write(*,*) 'PHYS TEST OH SUM : ',  sum(state%q(:ncol,:,idx))


! IS CO2 BEING OVERWRITTEN ELSEWHERE IN CAM?
! THE REASULT FROM GEOS_CHEM IS A POSITIVE VALUE
! BUT ONLY ZEROS ARE WRITTEN TO THE HISTORY FILES
! M.LONG, Sept. 12, 2012
   call cnst_get_ind('CO2',idx)
   call outfld('CO2', state%q(:,:,idx)  ,pcols   ,lchnk   )
   call cnst_get_ind('CO',idx)
   call outfld('CO',  state%q(:,:,idx)  ,pcols   ,lchnk   )
   call cnst_get_ind('CH4',idx)
   call outfld('CH4', state%q(:,:,idx)  ,pcols   ,lchnk   )
   call cnst_get_ind('CH2O',idx)
   call outfld('CH2O',state%q(:,:,idx)  ,pcols   ,lchnk   )
!   call cnst_get_ind('NOx',idx)
!   call outfld('NOx',state%q(:,:,idx)  ,pcols   ,lchnk   )


#endif


!---------------------------------------------------------------------------------------
! Save history variables. These should move to the appropriate parameterization interface
!---------------------------------------------------------------------------------------


   call outfld('PRECL   ',precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',precc   ,pcols   ,lchnk       )
   call outfld('PRECSL  ',precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',precsc  ,pcols   ,lchnk       )
   
   prect(:ncol) = precc(:ncol) + precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )


#if ( defined COUP_CSM )
   call outfld('PRECLav ',precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',precc   ,pcols   ,lchnk   )
#endif


#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ',prect   ,pcols   ,lchnk       )
#endif
!     
! Compute heating rate for dtheta/dt 
!
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )


! convert radiative heating rates to Q*dp for energy conservation
   if (conserve_energy) then
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if


   return
end subroutine tphysbc
