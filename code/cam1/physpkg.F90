#include <misc.h>
#include <params.h>

subroutine physpkg(phys_state, gw, ztodt, phys_tend, pbuf, dyn_tend)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics
! 
! Method: 
! COUP_CSM and must be checked in order to invoke the proper calling
! sequence for running the CSM model
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, masterproc
   use tracers,      only: set_state_pdry
   use ppgrid,       only: pcols, pver
   use buffer,       only: pblht, kbfs, tpert, qpert, qrs, qrl, cosp_cnt
#ifdef UWMT
   use buffer,       only: kvm, kvh
#endif
#ifdef NSTEP_MOIST
   use buffer,       only: dq_moist, dt_moist, du_moist, dv_moist, cnt, cnb
   use buffer,       only: cmfdqr, nevapr, prain, cmfmc, q0_conv, t0_conv, date0_conv
#endif
#ifdef NSTEP_GW
   use buffer,       only: utgw, vtgw, stgw, qtgw
#endif
#ifdef MOZART2
   use buffer,       only: chemflx
   use buffer,       only: dq_chem
   use constituents, only: cnst_get_ind
#endif

   use check_energy, only: check_energy_gmean
   use comsrf
   use comsrfdiag
#ifdef COUP_CSM
   use ccsm_msg,   only: ccsmave, dorecv, dosend, ccsmsnd, ccsmrcv
#else
   use atm_lndMod, only: atmlnd_drv
#endif
#ifdef SPMD
   use mpishorthand
#endif
   use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_allocate, pbuf_deallocate, &
                             pbuf_update_tim_idx
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physics_types,  only: physics_state, physics_tend, dynamics_tend
   use diagnostics,    only: diag_surf
   use time_manager,   only: get_nstep, is_first_step, is_first_restart_step, &
                             is_end_curr_month, get_curr_date, is_end_curr_day
   use physconst,      only: stebol
   use dycore,         only: dycore_is
   use constituents,   only: ppcnst,pcnst, iaero     !liuqx 20071227
   use ghg_surfvals,   only: chem_surfvals_get       ! for carbon cycle, wzz, 2008.11.20

#ifdef OBS-FORCE
   use time_manager,     only: obs_ncep, obs_fnl, obs_rain
   use atm_force_ncep,   only: tempint
   use rain_force_ncep,  only: rainint
   use atm_force_fnl,    only: tempint_fnl
   use atm_force_t639,   only: tempint_t639
#endif

#if (!defined COUP_CSM)
   use ice_constants, only: TfrezK
#endif
   use history,       only: outfld

#if (!defined COUP_CSM)
#if (!defined COUP_SOM)
   use sst_data, only: sstint
   use ice_data, only: iceint
#endif
#endif


#ifdef MOZART2
!   use mo_sulf,         only : read_sulf
!   use mo_hook,         only : moz_radon
!   use chemistry,       only : ixchm
!   use bcc_timestep_init, only: timestep_init
    use constituents,        only: cnst_get_ind
    use time_manager, only: get_curr_calday
    use hemco_interface, only: HCOI_Chunk_Run, hemco_sunrz, land_sunrz, lai_sunrz
    use hco_extra,   only: AREA_M2
    use hco_bcc_convert_state_mod, only:State_CAM_AREAM2
    use gchp_hook,         only : gc_radon
#endif

#ifdef BCCCHEM
   use chemistry,             only: ixchm
   use bcc_timestep_init,     only: bccchem_timestep_init
   use bcc_lightning,         only: lightning_no_prod
#endif

#ifdef WACCM
   use mo_lightning,           only: lightning_no_prod
   use chem_timestep_init,     only: timestep_init
   use efield,                 only: get_efield
   use mo_radheat,             only: radheat_timestep_init
#endif

#ifdef GEOSCHEM
   use chemistry,       only : ixchm
#endif

#ifdef MOZART2
   use time_manager, only: get_curr_calday
   use acbnd,        only: acbndint

   use prescribed_em_d3, only: em3int_IPCC
   use prescribed_em_d2, only: em2int_IPCC
   use prescribed_em_d2_anthro, only: em2int_IPCC_anthro   !zf2016.12.08
   use prescribed_emis_3d, only: em3d_int_IPCC          !zf2017.01.04
   use prescribed_fco2_data, only: fco2int_IPCC
   use prescribed_fco2_anthro, only: fco2int_IPCC_anthro

   use prescribed_ch4, only : ch4int   !zf

   use constituents, only: cnst_get_ind
   use pmgrid,       only: plon, plat
   use phys_grid,    only: get_ncols_p, scatter_field_to_chunk, gather_chunk_to_field
#endif

#ifdef CO2
   use prescribed_fco2_data, only: fco2int_IPCC
   use prescribed_fco2_anthro, only: fco2int_IPCC_anthro
   use prescribed_emis_3d,     only: em3d_int_IPCC  

   use constituents, only: cnst_get_ind
   use pmgrid,       only: plon, plat
   use phys_grid,    only: get_ncols_p, scatter_field_to_chunk, gather_chunk_to_field

#endif

#if (defined BCCCHEM) || (defined WACCM)
   use prescribed_fco2_data, only: fco2int_IPCC
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf

   type(dynamics_tend), intent(inout), dimension(begchunk:endchunk) :: dyn_tend
!
!---------------------------Local workspace-----------------------------
!
#ifdef MOZART2
   real(r8)  cldtop(pcols,begchunk:endchunk)
   real(r8)  cldbot(pcols,begchunk:endchunk)
   real(r8)  co2_lndflux(pcols,begchunk:endchunk)
   real(r8)  fco2(pcols,begchunk:endchunk)
   real(r8)  co2emis(pcols,begchunk:endchunk)
   real(r8)  co2emis_anthro(pcols,begchunk:endchunk)
   real(r8)  co2emis_fos(pcols,begchunk:endchunk)
   integer   idx_CO2

   real(r8) :: co2mmr_field      (plon,plat,pver)
   real(r8) :: co2mmr_lat        (plat,pver)
   real(r8) :: co2mmr_lat_allgrid(plon,plat,pver)
   real(r8) :: co2mmr_mozart        (pcols, begchunk:endchunk)
   real(r8) :: co2mmr_mozart_latmean(pcols,pver,begchunk:endchunk)

    real(r8),  allocatable  :: tflx(:,:,:,:)
    real(r8),  allocatable  :: eflx(:,:,:,:)
    real(r8)  eflxtobc(pcols,pver,begchunk:endchunk,pcnst)
    real(r8) :: temp(plon,pver,plat,pcnst)
    real(r8) :: temp2(plon,plat,73)
    real(r8) :: temp3(plon,plat,73)
    real(r8) :: lands(pcols,begchunk:endchunk,73)
    real(r8) :: lais(pcols,begchunk:endchunk,73)
#endif

#ifdef CO2
   real(r8)  fco2(pcols,begchunk:endchunk)
   real(r8)  co2emis(pcols,begchunk:endchunk)
   real(r8)  co2emis_anthro(pcols,begchunk:endchunk)
   real(r8)  co2emis_fos(pcols,begchunk:endchunk)

   real(r8) :: co2mmr_field      (plon,plat,pver)
   real(r8) :: co2mmr_lat        (plat,pver)
   real(r8) :: co2mmr_lat_allgrid(plon,plat,pver)
   real(r8) :: co2mmr_mozart        (pcols, begchunk:endchunk)
   real(r8) :: co2mmr_mozart_latmean(pcols,pver,begchunk:endchunk)

   integer idx_co2

#endif
   real    :: calday
   integer :: i,m,lat,c,lchnk,j,k                 ! indices
   integer :: lats(pcols)                       ! array of latitude indices
   integer :: lons(pcols)                       ! array of longitude indices
   integer :: ncol                              ! number of columns
   integer :: nstep                             ! current timestep number
   integer :: ncdate                            ! current date in integer format [yyyymmdd]
   integer :: ncsec                             ! current time of day [seconds]
   integer :: yr, mon, day                      ! year, month, and day components of a date

   real(r8) fsds(pcols,begchunk:endchunk)        ! Surface solar down flux
   real(r8) :: tmp(pcols,begchunk:endchunk)

!-----------------------------------------------------------------------

   call t_startf ('physpkg_st')
   nstep = get_nstep()

   call pbuf_allocate('physpkg')

! Compute total energy of input state and previous output state
   call t_startf ('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf, ztodt, nstep)
   call t_stopf ('chk_en_gmean')

#ifdef OBS-FORCE
!----------------------------------------------------------------
! New added variable
!----------------------------------------------------------------
   if (obs_ncep) then     ! jwh 2010.9.10
         call t_startf ('tempint')
         call tempint  (phys_state)
         call t_stopf  ('tempint')
   else if (obs_fnl) then   !jwh 2013
         call t_startf ('tempint_fnl')
         call tempint_fnl  (phys_state)
         call t_stopf  ('tempint_fnl')
   else   ! use the T639 data
         call t_startf ('tempint_t639')
         call tempint_t639  (phys_state)
         call t_stopf  ('tempint_t639')
   endif
!
   if (obs_rain) then     ! jwh 2010.9.10
         call t_startf ('rainint')
         call rainint  (phys_state)
         call t_stopf  ('rainint')
   end if
#endif

!===================================================================
#if (defined MOZART2)

   calday = get_curr_calday()
!   call em3int_IPCC
!   call em2int_IPCC
!   call em2int_IPCC_anthro        !zf2016.12.08
!   call em3d_int_IPCC               !zf2017.01.04

!   call timestep_init( phys_state )
!   call fco2int_IPCC( co2emis , fco2 )
!   call fco2int_IPCC_anthro( co2emis_anthro , co2emis_fos )
!   co2emis(:,:) = co2emis_anthro(:, :)

!   do c=begchunk, endchunk
!      call outfld('CO2E_ANT', co2emis_anthro(:, c),pcols   ,c     )
!      call outfld('CO2E_FOS', co2emis_fos(:, c),pcols   ,c     )
!   end do

!   call read_sulf
!   call moz_radon

!   call ch4int   !zf
    allocate(tflx(plon,pver,plat,pcnst))
    allocate(eflx(pcols,pver,begchunk:endchunk,pcnst))
 
    call HCOI_Chunk_Run( phys_state, phase=2)
 
    if (masterproc) then
    do i=1,pver
    tflx(:,i,:,:)=hemco_sunrz(:,:,i,:)
    enddo
    do i=1,plon/2
    temp(i,:,:,:)=tflx(i+plon/2,:,:,:)
    temp(i+plon/2,:,:,:)=tflx(i,:,:,:)
    temp2(i,:,:)=land_sunrz(i+plon/2,:,:)
    temp2(i+plon/2,:,:)=land_sunrz(i,:,:)
    temp3(i,:,:)=lai_sunrz(i+plon/2,:,:)
    temp3(i+plon/2,:,:)=lai_sunrz(i,:,:)
    enddo
    tflx=temp
    land_sunrz=temp2
    lai_sunrz=temp3
    endif
 do i=1,pcnst
    call scatter_field_to_chunk(1,pver,1,plon,tflx(:,:,:,i),eflx(1,1,begchunk,i))
 enddo
    eflxtobc=eflx
 do i=1,73
    call scatter_field_to_chunk(1,1,1,plon,land_sunrz(:,:,i),lands(1,begchunk,i))
    call scatter_field_to_chunk(1,1,1,plon,lai_sunrz(:,:,i),lais(1,begchunk,i))
 enddo
 
     IF ( ALLOCATED( hemco_sunrz ) ) DEALLOCATE( hemco_sunrz )
     IF ( ALLOCATED( land_sunrz  ) ) DEALLOCATE( land_sunrz )
     IF ( ALLOCATED( lai_sunrz  ) ) DEALLOCATE( lai_sunrz )
#endif

#ifdef CO2
   call fco2int_IPCC( co2emis , fco2 )
   call fco2int_IPCC_anthro( co2emis_anthro , co2emis_fos )
   co2emis(:,:) = co2emis_anthro(:, :)

   call em3d_int_IPCC

   do c=begchunk, endchunk
      call outfld('CO2E_ANT', co2emis_anthro(:, c),pcols   ,c     )
      call outfld('CO2E_FOS', co2emis_fos(:, c),pcols   ,c     )
   end do
#endif

!#if (defined MOZART2) || (defined CO2) 
!-------------------------------------
!
!   call cnst_get_ind('CO2', idx_co2)
!   do k=1, pver
!     do c= begchunk,endchunk
!        ncol = get_ncols_p(c)
!        co2mmr_mozart(:ncol,c) = phys_state(c)%q(:ncol,k,idx_co2)
!     enddo
!     call gather_chunk_to_field(1,1,1,plon,co2mmr_mozart,co2mmr_field(1,1,k))   
!   enddo
!   if( masterproc ) then 
!      do k=1, pver
!         do lat=1, plat
!           co2mmr_lat(lat,k) = 0.0
!           do i=1, plon
!             co2mmr_lat(lat,k) = co2mmr_lat(lat,k) + co2mmr_field(i,lat,k) /float(plon)
!           enddo
!           do i=1, plon
!             co2mmr_lat_allgrid(i,lat,k) = co2mmr_lat(lat,k)
!           enddo
!         enddo
!      enddo
!   endif
!   call scatter_field_to_chunk(1, pver, 1, plon, co2mmr_lat_allgrid, co2mmr_mozart_latmean)
!
!#endif

#ifdef BCCHEM
!
! Time interpolate for chemistry, if appropriate
!
   call em3int_IPCC
   call em2int_IPCC
   call bccchem_timestep_init( phys_state )
   call fco2int_IPCC( co2emis , fco2 )
#endif

#ifdef WACCM

  ! Time interpolate for chemistry.
  !
   call em3int_IPCC
   call em2int_IPCC
!  call fco2int_IPCC( co2emis , fco2 )

   call timestep_init

!! Upper atmosphere radiative processes
!
!   call radheat_timestep_init

  ! Compute the electric field
   call get_efield

#endif

!=============================================================
! Advance time information
!-----------------------------------------------------------------------

   call advnce()
   call t_stopf ('physpkg_st')
!
!  set the state vector dry quantities
!

   if ( .not. dycore_is('LR') )then
      ! for LR, this is done in d_p_coupling since dynamics is called first

!$OMP PARALLEL DO PRIVATE (C,NCOL)

      do c=begchunk, endchunk
         call set_state_pdry( phys_state(c)) 
      end do

   endif ! .not. dycore_is('LR')

#ifdef TRACER_CHECK
   call gavglook ('before tphysbc DRY', phys_state, gw)
#endif


!-----------------------------------------------------------------------
! Tendency physics before flux coupler invokation
!-----------------------------------------------------------------------
!

#if (defined BFB_CAM_SCAM_IOP )
   do c=begchunk, endchunk
      call outfld('Tg',srfflx_state2d(c)%ts,pcols   ,c     )
   end do
#endif
   call t_startf ('bc_physics')

!$OMP PARALLEL DO PRIVATE (C,NCOL)

   do c=begchunk, endchunk
      call t_startf ('tphysbc')

#ifdef TDK_CU
      call tphysbc_tdk (ztodt, pblht(1,c), tpert(1,c),                    &
	              srfflx_state2d(c)%ts, srfflx_state2d(c)%sst,        &
                      qpert(1,1,c), surface_state2d(c)%precl,             &
	   	      surface_state2d(c)%precc, surface_state2d(c)%precsl,&
                      surface_state2d(c)%precsc,                          &
                      srfflx_state2d(c)%asdir, srfflx_state2d(c)%asdif,   &
                      srfflx_state2d(c)%aldir, srfflx_state2d(c)%aldif,   &
                      snowhland(1,c),                                     &
                      qrs(1,1,c), qrl(1,1,c), surface_state2d(c)%flwds,   &
                      fsns(1,c), fsnt(1,c),                               &
                      flns(1,c),    flnt(1,c), srfflx_state2d(c)%lwup,    &
                      surface_state2d(c)%srfrad, surface_state2d(c)%sols, &
                      surface_state2d(c)%soll, surface_state2d(c)%solsd,  &
                      surface_state2d(c)%solld,                           &
                      phys_state(c), phys_tend(c),                        &
	              pbuf, prcsnw(1,c), fsds(1,c), landm(1,c), landfrac(1,c), &
	              ocnfrac(1,c),icefrac(1,c)  ,                             &
                      dyn_tend(c), srfflx_state2d(c)%cflx           )  ! wtw 200508 add the last line
#else
       
      call tphysbc (ztodt, pblht(1,c), kbfs(1,c)   , tpert(1,c),                      &
                    srfflx_state2d(c)%ts, srfflx_state2d(c)%sst,        &
                    qpert(1,1,c), surface_state2d(c)%precl,             &
                    surface_state2d(c)%precc, surface_state2d(c)%precsl,&
                    surface_state2d(c)%precsc,                          &
                    srfflx_state2d(c)%asdir, srfflx_state2d(c)%asdif,   &
                    srfflx_state2d(c)%aldir, srfflx_state2d(c)%aldif,   &
                    snowhland(1,c),                                     &
                    qrs(1,1,c), qrl(1,1,c), surface_state2d(c)%flwds,   &
                    fsns(1,c), fsnt(1,c),                               &
                    flns(1,c),    flnt(1,c), srfflx_state2d(c)%lwup,    &
                    surface_state2d(c)%srfrad, surface_state2d(c)%sols, &
                    surface_state2d(c)%soll, surface_state2d(c)%solsd,  &
                    surface_state2d(c)%solld,                           &
                    phys_state(c), phys_tend(c),                        &
                    pbuf, prcsnw(1,c), fsds(1,c), landm(1,c), landfrac(1,c), &
                    ocnfrac(1,c),icefrac(1,c),srfflx_parm2d(c)%sw1,     &
                    cosp_cnt(c),                                        &
#ifdef NSTEP_MOIST                  
                    cnt(1,c), cnb(1,c), cmfdqr(1,1,c),  nevapr(1,1,c),  &
                    prain(1,1,c), cmfmc(1,1,c),                         &
                    dt_moist(1,1,c), dq_moist(1,1,1,c), du_moist(1,1,c), dv_moist(1,1,c), &
                    q0_conv(1,1,c),  t0_conv(1,1,c) , date0_conv(1,c), &
#endif 
#ifdef MOZART2
                    chemflx(1,1,c), dq_chem(1,1,1,c),                   &
                    srfflx_state2d(c)%flux_ISOP, srfflx_state2d(c)%flux_ACET,     &
                    srfflx_state2d(c)%flux_C3H6, srfflx_state2d(c)%flux_C2H4,     &
                    srfflx_state2d(c)%flux_OC2,  srfflx_state2d(c)%flux_C10H16,   &
                    srfflx_state2d(c)%flux_CO2,  srfflx_state2d(c)%flux_N2O,      &
                    srfflx_state2d(c)%flux_DMS,                                   &
                    co2mmr_mozart_latmean(1,1,c),srfflx_state2d(c)%shf, srfflx_state2d(c)%lhf,&
                    srfflx_state2d(c)%u_10m,srfflx_state2d(c)%v_10m,srfflx_state2d(c)%wsx,srfflx_state2d(c)%wsy,& 
                    eflxtobc(:,:,c,:), lands(:,c,:),lais(:,c,:),State_CAM_AREAM2(:,c),&
#endif
#ifdef CO2
                    co2mmr_mozart_latmean(1,1,c),                                 &
#endif
#ifdef HAMOCC
                    surface_state2d(c)%dmsvmr,   surface_state2d(c)%n2ovmr,       &
#endif
                    surface_state2d(c)%co2, dyn_tend(c) )


#endif
      call t_stopf ('tphysbc')

      if (dosw .or. dolw) then
	call output_flns_fsns_fluxes(surface_state2d(c),c)
      end if	

!-----------------------------------------------------

#if ( ! defined COUP_CSM )
!
! zero surface fluxes at beginning of each time step.  Land Ocean and Ice
! processes will will write into process specific flux variables
! at the end of the time step these separate fluxes will be combined over the
! entire grid
!
      call srfflx_state_reset (srfflx_state2d(c))
#endif

   end do
   call t_stopf ('bc_physics')

!----------------------------------------------------------------------
! Set lightning production of NO
!-------------------------------------------------------------
#ifdef BCCCHEM
   call lightning_no_prod( phys_state, pbuf, landfrac, ocnfrac, &
                           cldtop,    cldbot )
#endif

#ifdef WACCM
   call lightning_no_prod( phys_state, pbuf , cldtop,    cldbot )
#endif

!---------------------------------------------------------
#ifdef TRACER_CHECK
   call gavglook ('between DRY', phys_state, gw)
#endif

#if ( ! defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities - no flux coupler
!-----------------------------------------------------------------------
!
   if (.not. aqua_planet) then
!
! Call land model driving routine
!
#ifdef TIMING_BARRIERS
      call t_startf ('sync_tphysbc_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_tphysbc_lnd')
#endif
      call t_startf ('atmlnd_drv')

#if ( defined SCAM )
       if (landfrac(1,begchunk).gt.0) &
#endif

      call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0,&
                      mvelpp,surface_state2d,srfflx_parm2d)

      call t_stopf ('atmlnd_drv')

#ifdef TIMING_BARRIERS
      call t_startf ('sync_after_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_after_lnd')
#endif
!
!   if ( masterproc ) then
!   write(0,*)'wzz: cflx bc_physics0=', srfflx_state2d(c)%cflx(2,43), srfflx_state2d(c)%cflx(2,1)
!   endif
! save off albedos and longwave for som offline vars
!
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            if (landfrac(i,c) > 0.) then
               asdirlnd(i,c) = srfflx_parm2d(c)%asdir(i)
               asdiflnd(i,c) = srfflx_parm2d(c)%asdif(i)
               aldirlnd(i,c) = srfflx_parm2d(c)%aldir(i)
               aldiflnd(i,c) = srfflx_parm2d(c)%aldif(i)
               lwuplnd(i,c)  = srfflx_parm2d(c)%lwup(i)
            else
               asdirlnd(i,c) = 0. 
               asdiflnd(i,c) = 0. 
               aldirlnd(i,c) = 0. 
               aldiflnd(i,c) = 0. 
               lwuplnd(i,c)  = 0. 
            end if
         end do
!
!output shf/lhf fluxes for land model
!
         call output_shf_lhf_fluxes(srfflx_parm2d(c), c, ncol, landfrac(1,c), 'LND')
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)
      end do
   end if                    ! end of not aqua_planet if block

#if (defined COUP_SOM)
!
! Set ocean surface quantities - ocn model internal to atm
!
   if (is_end_curr_day ()) then
      call print_coverage ('icefrac', ' million km^2', icefrac, 1.d-12)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*sicthk(i,c)
         end do
      end do
      call print_coverage ('icevol ', ' 10^13m^3', tmp, 1.d-13)

      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*snowhice(i,c)
         end do
      end do
      call print_coverage ('snowvol', ' 10^13m^3', tmp, 1.d-13)
   end if

   call t_startf ('somint')
   call somint ()
   call t_stopf ('somint')

   call t_startf ('somoce')
   call somoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('somoce')

#else

   call t_startf ('sstint')
   call sstint ()
   call t_stopf ('sstint')
!
! iceint may change ocean fraction, so call it before camoce
!
   call t_startf ('iceint')
   call iceint ()
   call t_stopf ('iceint')

   call t_startf ('camoce')
   call camoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('camoce')
#endif
!
! Set ice surface quantities - icn model internal to atm
!
   call t_startf('camice')
   call camice (surface_state2d, srfflx_parm2d)
   call t_stopf('camice')
!
! output shf/lhf fluxes for ice/ocn/som_offline 
!
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do i=1,ncol
         if(icefrac(i,c) > 0.) then
            tsice_rad(i,c) = sqrt(sqrt(srfflx_parm2d(c)%lwup(i)/stebol))
         else
            tsice_rad(i,c) = TfrezK
         endif
      end do
      call output_shf_lhf_fluxes (srfflx_parm2d(c), c, ncol, icefrac(1,c), 'ICE')
      call output_shf_lhf_fluxes (srfflx_parm2d_ocn(c), c, ncol, ocnfrac(1,c), 'OCN')
      call output_shfoi_lhfoi_fluxes (srfflx_parm2d_ocn(c), srfflx_parm2d(c), c)

!JR SOM case: Have to wait to call update routine till after both ocean and ice have
!JR operated, since the fractions can change internal to the parameterization
      do i = 1, ncol
         srfflx_state2d(c)%sst(i) = srfflx_parm2d_ocn(c)%ts(i)
      enddo
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d_ocn(c), ocnfrac(1,c), ncol)
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), icefrac(1,c), ncol)
   end do
#endif

#if ( defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities using csm flux coupler
!-----------------------------------------------------------------------
!
! If send data to flux coupler only on radiation time steps:
!
   if (flxave) then
!
! Average the precipitation input to lsm between radiation calls.
!
      call ccsmave(iradsw, nstep, dosw)
!
! Use solar radiation flag to determine data exchange steps 
! with flux coupler. This processes are not independent since 
! instantaneous radiative fluxes are passed, valid over the 
! interval to the next radiation calculation. The same 
! considerations apply to the long and shortwave fluxes, so 
! the intervals must be the same. Data is received from the 
! coupler one step after it is sent.
!
      if (nstep == 0) then
         dorecv = .true.
         dosend = .true.
      else if (nstep == 1) then
         dorecv = .false.
         dosend = .false.
      else if ( (nstep == 2) .and. (iradsw == 1) ) then
         dorecv = .true.
         dosend = dosw
      else
         dorecv = dosend
         dosend = dosw
      end if
   endif
!
! If send data to flux coupler on every time step
!
   if (.not. flxave) then
      if (nstep /= 1) then
         dorecv = .true.
         dosend = .true.
      else 
         dorecv = .false.
         dosend = .false.
      endif
   endif
!
! Send/recv data to/from the csm flux coupler.
!
   if (dosend) call ccsmsnd ( )
   if (dorecv) call ccsmrcv ( )
#endif
!
!-----------------------------------------------------------------------
! Tendency physics after coupler 
! Not necessary at terminal timestep.
!-----------------------------------------------------------------------
!
   call t_startf ('ac_physics')

!$OMP PARALLEL DO PRIVATE (C, NCOL)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)

!
! surface diagnostics for history files
!
      call diag_surf (srfflx_state2d(c), surface_state2d(c), icefrac(1,c), ocnfrac(1,c), landfrac(1,c), &
                      sicthk(1,c), snowhland(1,c), snowhice(1,c), tsice(1,c), trefmxav(1,c), &
                      trefmnav(1,c) )

      call t_startf ('tphysac')

      call tphysac (ztodt, pblht(1,c), kbfs(1,c), qpert(1,1,c), tpert(1,c), srfflx_state2d(c)%shf,        &
             srfflx_state2d(c)%wsx,srfflx_state2d(c)%wsy,srfflx_state2d(c)%cflx,sgh(1,c),srfflx_state2d(c)%lhf, &
             landfrac(1,c), snowhland(1,c),srfflx_state2d(c)%tref,surface_state2d(c)%precc,surface_state2d(c)%precl,    &
             surface_state2d(c)%precsc, surface_state2d(c)%precsl, phys_state(c), phys_tend(c), pbuf, &
             ocnfrac(1,c), fsds(1,c), icefrac(1,c), fv(1,c), ram1(1,c), srfflx_state2d(c)%ts,  &
#ifdef MOZART2
             chemflx(1,1,c),                  &
             srfflx_state2d(c)%flx_mss_vrt_dst01, srfflx_state2d(c)%flx_mss_vrt_dst02, &
             srfflx_state2d(c)%flx_mss_vrt_dst03, srfflx_state2d(c)%flx_mss_vrt_dst04, &
             srfflx_state2d(c)%carbon, co2emis(1,c) , fco2(1,c),  &
             srfflx_state2d(c)%vlc_trb_SO4, srfflx_state2d(c)%vlc_trb_NH4, &
             srfflx_state2d(c)%vlc_trb_CB2, srfflx_state2d(c)%vlc_trb_SO2, &
             srfflx_state2d(c)%u_10m, srfflx_state2d(c)%v_10m,          &    !zf 2015-05-19
             srfflx_state2d(c)%lnd_ram1,  srfflx_state2d(c)%lnd_fv,        &
#endif
#ifdef CO2
             srfflx_state2d(c)%carbon, co2emis(1,c) , fco2(1,c),  &
#endif
#ifdef NSTEP_GW
             utgw(1,1,c), vtgw(1,1,c), stgw(1,1,c), qtgw(1,1,1,c), &
#endif
#ifdef UWMT
             kvm(1,1,c), kvh(1,1,c), &
#endif
             srfflx_parm2d(c)%sw1)
      call t_stopf ('tphysac')
!    write(0,*)'wzz: c1=',c

   end do                    ! Chunk loop

   call t_stopf('ac_physics')

#ifdef TRACER_CHECK
   call gavglook ('after tphysac FV:WET)', phys_state, gw )
#endif

   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()

end subroutine physpkg
