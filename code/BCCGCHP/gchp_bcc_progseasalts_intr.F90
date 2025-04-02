!xiaolu
!---------------------------------------
!Purpose: DUST AEROSOL DRY DEPOSITION and EMISSION
!from bcc_progseasalts_intr.F90

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/15
!---------------------------------------

#include <misc.h>
#include <params.h>

module gchp_bcc_progseasalts_intr

  !---------------------------------------------------------------------------------
  ! Module to interface the aerosol parameterizations with CAM
  ! written by PJR (extensively modified from chemistry module)
  ! prognostic sea salt module taken from dust module by nmm 11/12/03
  !---------------------------------------------------------------------------------

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use pmgrid,                only: masterproc
  use ppgrid,                only: pcols, pver,pverp
  use physconst,             only: mwdry, mwh2o,gravit,rair
  use constituents,          only: ppcnst, cnst_add, cnst_name, cnst_get_ind
  use gchp_bcc_drydep_mod,        only: calcram !mo->gchp
  use gchp_bcc_dust_sediment_mod, only: dust_sediment_tend !mo->gchp
  use gchp_bcc_chem_mods,         only: nsst !mo->gchp
  use gchp_chem_utls,          only: get_spc_ndx

  implicit none

  private          ! Make default type private to the module

  save

  ! 
  !  character(len=8), dimension(ncnst), parameter :: & ! constituent names
  !       cnst_names = (/'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04'/)

  !
  ! Public interfaces
  !
!  public bcc_progseasalts_wet                                  ! interface to wet deposition
  public bcc_progseasalts_emis_intr                            ! interface to emission
  public progseasalts_drydep_intr                          ! interface to tendency computation

!  xiaolu,2017/06/16 
!  real(r8) :: stk_crc(nsst)=(/ 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8 /) 
   real(r8) :: stk_crc(nsst)=(/ 1.0_r8, 1.0_r8/)

  ![frc] Correction to Stokes settling velocity--we are assuming they are close enough 1.0 to be
  !       set to one--unlikely to cause considerable error
  !
  real(r8),parameter:: dns_aer_sst = &
                       2200.0_r8     ![kg m-3] Aerosol density

  !(diameter not radius as in tie)
!--------zftest-------------------------------------------------------------------------------
!zftest  real(r8) :: sst_source(nsst)  &                                 ! source factors for each bin 
!zftest              =(/4.77e-15_r8, 5.19e-14_r8, 1.22e-13_r8, 6.91e-14_r8/)

  !xiaolu,2017/06/16
!  real(r8) :: sst_source(nsst)  &      ! source factors for each bin 
!     =(/4.77e-15_r8*0.65, 5.19e-14_r8*0.65, 1.22e-13_r8*0.71, 6.91e-14_r8*0.76/)

!  real(r8) :: sst_source(nsst)  &      ! source factors for each bin
!      =(/4.77e-15_r8*0.65, 5.19e-14_r8*0.65/)

!xiaolu, 2018/08/26:change sst_source to be consistent with Jaegle et al. (2011) ACP
!xiaolu do sensitivity test for seasalt
real(r8) :: sst_source(nsst)  &      ! source factors for each bin
      =(/4.75e-15_r8*0.65*0.9, 5.19e-14_r8*0.65*40*0.3*0.9*0.5*0.6/)

!             =(/4.77e-15_r8*0.6, 5.19e-14_r8*0.56, 1.22e-13_r8*0.46, 6.91e-14_r8*0.7/)
!--------zftest-------------------------------------------------------------------------------

  !xiaolu,2017/06/16
!  real(r8) :: smt_vwr(nsst) & 
 !             =(/0.52e-6_r8,2.38e-6_r8,4.86e-6_r8,15.14e-6_r8/) ![m] Mass-weighted mean diameter resolved
 real(r8) :: smt_vwr(nsst) &
              =(/0.25e-6_r8,4.5e-6_r8/) ![m] Mass-weighted mean diameter resolved


contains


!xiaolu,comment wet depostion

    subroutine progseasalts_drydep_intr (lchnk, ncol, q, pdel, t, pmid, lq, dq,  &
	                        dt, lat, clat, &
                            fsds, obklen, ts, ustar, pblh, month, landfrac, &
                            icefrac, ocnfrac,fvin,ram1in)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to dry deposition and sedimentation of progseasalts
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: Natalie Mahowald and Phil Rasch
    ! 
    !-----------------------------------------------------------------------
    use history,               only: outfld
    use phys_grid,             only: get_lat_all_p
    use constituents,          only: cnst_name
    use gchp_bcc_drydep_mod,        only: d3ddflux  
    use gchp_bcc_dust_sediment_mod, only: dust_sediment_tend

    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
	integer, intent(in)              :: lchnk
	integer, intent(in)              :: ncol
    real(r8),            intent(in)  :: dt             ! time step

    logical, intent(inout) :: lq(ppcnst)
    real(r8), intent(inout) :: dq(pcols,pver,ppcnst)

    real(r8), intent(in) :: q(pcols,pver,ppcnst)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)

    integer, intent(in)  :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                 ! latitude 
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel--used over oceans and sea ice.
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    integer, intent(in)  :: month
    real(r8), intent(in) :: fvin(pcols)         ! for dry dep velocities from land model for progseasalts
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model for progseasalts

    !
    ! Local variables
    !
    integer :: m                                  ! tracer index
    integer :: mm                                  ! tracer index
    real(r8) :: tvs(pcols,pver)
    real(r8) :: dvel(pcols)            ! deposition velocity
    real(r8) :: sflx(pcols)            ! deposition flux
    real(r8) :: vlc_dry(pcols,pver,nsst)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,nsst)            ! dep velocity
    real(r8)::  vlc_trb(pcols,nsst)            ! dep velocity
    real(r8)::  dep_trb(pcols)       !kg/m2/s
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_dry_tend(pcols,pver)       !kg/kg/s (total of grav and trb)
    real(r8) :: obuf(1)
    real(r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real(r8) :: pvprogseasalts(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: tsflx(pcols)
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice
    integer  :: ix_progseasalts

    integer :: i,k
    real(r8) :: oro(pcols)
    !
    !-----------------------------------------------------------------------

    tvs(:ncol,:) = t(:ncol,:)   !*(1+q(:ncol,k)
    rho(:ncol,:)=  pmid(:ncol,:)/(rair*t(:ncol,:))
    ! calculate oro--need to run match

    !   Dry deposition of Progseasalts Aerosols
    !   #################################
    !  we get the ram1,fv from the land model as ram1in,fvin,, but need to calculate it over oceans and ice.  
    !  better if we got thse from the ocean and ice model
    !  for friction velocity, we use ustar (from taux and tauy), except over land, where we use fv from the land model.

    ! copy fv,ram1 values from land model to variables to modify in calcram
    !    fv = fvin    
    !    ram1 = ram1in

    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
                 ustar,ram1in,ram1,t(:,pver),pmid(:,pver),&
                 pdel(:,pver),fvin,fv)

    ! this is the same as in the dust model--perhaps there is someway to 
    !calculate them up higher in the model

    call ProgseasaltsDryDep(ncol,t(:,:),pmid(:,:),q(:,:,:),ram1,fv,vlc_dry,vlc_trb,vlc_grv)

    !xiaolu change SSLT01 to SALA
    call cnst_get_ind('SALA', ix_progseasalts)

    tsflx(:)=0._r8
    do m=1,nsst
       mm = ix_progseasalts + m-1

       ! use pvprogseasalts instead (means making the top level 0)

       pvprogseasalts(:ncol,1)=0._r8
       pvprogseasalts(:ncol,2:pverp) = vlc_dry(:ncol,:,m)

       sflx = 0.
       call d3ddflux(ncol, vlc_dry(:,:,m), q(:,:,mm),pmid,pdel, tvs,sflx,dq(:,:,mm),dt)

       sflx(:) = 0.
       do k=1, pver
          do i=1, ncol
             sflx(i) = sflx(i) + dq(i,k,mm)*pdel(i,k)/gravit
          enddo
       enddo 

       ! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i,m)/vlc_dry(i,pver,m)
          dep_grv(i)=sflx(i)*vlc_grv(i,pver,m)/vlc_dry(i,pver,m)
       enddo
       tsflx(:ncol)=tsflx(:ncol)+sflx(:ncol)

!        write(*,*),'xiaolu check sea salt dry dep',sflx(:ncol)


!-------------------------2015-05-29---------------------------------
        
        if (m == 1) call outfld('DRY_SLT1', sflx, pcols, lchnk)
!        if (m == 2) call outfld('DRY_SLT2', sflx, pcols, lchnk)
!        if (m == 3) call outfld('DRY_SLT3', sflx, pcols, lchnk)
!        if (m == 4) call outfld('DRY_SLT4', sflx, pcols, lchnk)

!        if (m == 1) call outfld('GRV_SLT1', dep_grv, pcols, lchnk)
!        if (m == 2) call outfld('GRV_SLT2', dep_grv, pcols, lchnk)
!        if (m == 3) call outfld('GRV_SLT3', dep_grv, pcols, lchnk)
!        if (m == 4) call outfld('GRV_SLT4', dep_grv, pcols, lchnk)
!--------------------------------------------------------------------

       ! set flags for tendencies (water and 4 ghg's)
       lq(mm) = .TRUE.
    end do
    return
    end subroutine progseasalts_drydep_intr

!=================================================================================================
    subroutine bcc_progseasalts_emis_intr (lchnk, ncol, u, v, zm, pdel, dt,sflx,ocnfrc, icefrc, u10, v10)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to emission of all progseasaltss.
    ! Derives from Xue Xie Tie's seasalts in MOZART, which derives from 
    !   Chin et al., 2001 (JAS) and Gong et al., 1997 (JGR-102,3805-3818).
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: Phil Rasch and Natalie Mahowald
    !
    ! Derives from Martensson et al. (2003) (JGR-108, 4297,doi:10.1029/2002JD002263)
    ! valid from 20nm to ~2500nm dry diameter (based on lab experiment with artificial sea water)
    !
    ! currently we recommend that it is combined with
    ! the parameterisation by Monahan et al. (1986) for
    ! dry diameters > 2-3 um even if it lacks
    ! temperature dependence (despite that Bowyer et
    ! al. (1990) found a similar dependency in the tank
    ! from Monahan et al. (1986))
    !
    !-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday
    use time_manager,  only: is_perpetual

    !   use progseasalts, only: progseasaltssf
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    integer, intent(in)              :: lchnk
    integer, intent(in)              :: ncol
    real(r8), intent(in)             :: u10(pcols), v10(pcols)     !zf 2015-05-19
    real(r8), intent(in)             :: pdel(pcols,pver)
    real(r8),            intent(in)  :: dt             ! time step
    real(r8),            intent(inout) :: sflx(pcols, ppcnst)
    real(r8):: ocnfrc(pcols)
    real(r8):: icefrc(pcols)
    real(r8):: u(pcols, pver)
    real(r8):: v(pcols, pver)
    real(r8):: zm(pcols, pver)

    integer  lat(pcols)                  ! latitude index 
    integer  lon(pcols)                  ! longitude index
    integer i
    integer m,mm

    real(r8):: tot_sflx(pcols,2)                ! accumulate over all bins for output ;2 SSLT

    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: u10cubed(pcols)
    !    real (r8), parameter :: z0=0.5  ! m roughness length over oceans--from Tie.
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model
    !

    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: ix_progseasalts

    real(r8) :: emis_sslt(pcols)

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

!    tot_sflx(:)=0._r8

!zf    u10cubed(:ncol)=sqrt(u(:ncol,pver)**2+v(:ncol,pver)**2)

    ! move the winds to 10m high from the midpoint of the gridbox:
    ! follows Tie and Seinfeld and Pandis, p.859 with math.

!zf    u10cubed(:ncol)=u10cubed(:ncol)*log(10._r8/z0)/log(zm(:ncol,pver)/z0)
    u10cubed(:ncol)=sqrt(u10(:ncol)**2+v10(:ncol)**2)    !zf 2015-05-19

    ! we need them to the 3.41 power, according to Gong et al., 1997:
    u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8


!    call cnst_get_ind('SSLT01', ix_progseasalts)
    call cnst_get_ind('SALA', ix_progseasalts)
!    call cnst_get_ind('SALC', ix_progseasalts)

    tot_sflx(:,:)=0._r8 
    do m=1,nsst
       mm = ix_progseasalts + m -1

       do i=1, ncol
          
          if( icefrc(i) .lt. 0.95) then

              !sflx already includes drydep flux of SSLT
              sflx(i,mm)  = sflx(i,mm) + sst_source(m)* u10cubed(i)*max((ocnfrc(i)-icefrc(i)),0.0)
              tot_sflx(i,m) = tot_sflx(i,m) + sst_source(m)* u10cubed(i)*max((ocnfrc(i)-icefrc(i)),0.0)
          endif
       enddo
       ! this is being done inside of the vertical diffusion automatically
       !         ptend%lq(m) = .true. ! tendencies for all progseasalts on
       !         ptend%q(:ncol,pver,mm) = cflx(:ncol,m)*gravit/state%pdel(:ncol,pver)

!----------------------------------------------------------------------
       if(m == 1) call outfld('EMS_SLT1', tot_sflx(:, m), pcols, lchnk)
       if(m == 2) call outfld('EMS_SLT2', tot_sflx(:, m), pcols, lchnk)
!       if(m == 3) call outfld('EMS_SLT3', sflx(:, mm), pcols, lchnk)
!       if(m == 4) call outfld('EMS_SLT4', sflx(:, mm), pcols, lchnk)
!----------------------------------------------------------------------

    enddo

    return
  end subroutine bcc_progseasalts_emis_intr

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine ProgseasaltsDryDep(c)
  !
  ! !INTERFACE:
  !
  subroutine ProgseasaltsDryDep(ncol,t,pmid,q,ram1,fv,vlc_dry,vlc_trb,vlc_grv)
    !
    ! !DESCRIPTION: 
    !
    ! Dry deposition for seasalts: modified from dust dry deposition following
    ! Xue Xie Tie's MOZART seasalts by NMM.  Sam Levis did hte first dust dry 
    ! deposition in CLM/CCSM2 from Charlie Zender's codes, modified by NMM for CAM.
    ! Sam's Notes in CLM:
    ! Determine Turbulent dry deposition for progseasalts. Calculate the turbulent 
    ! component of progseasalts dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the progseasalts dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, progseasaltsini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
    ! Note that because sea salts' radius changes with humidity we cannot 
    ! precalculate slip coefficients, etc. in the initialization subroutine, but
    ! we have to calculate them at the time--how slow???
    !
    ! !USES
    !
    use physconst,     only: rair,pi,boltz
    use wv_saturation, only: aqsat
    use history,       only: outfld

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8) :: q(pcols,pver,ppcnst)       !atm temperature (K)
    real(r8) :: t(pcols,pver)       !atm temperature (K)
    real(r8) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8) :: rho(pcols,pver)     !atm density (kg/m**3)
    real(r8) :: fv(pcols)           !friction velocity (m/s)
    real(r8) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: vlc_trb(pcols,nsst)  !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,nsst)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,nsst)  !dry deposn velocity (m/s)
    integer, intent(in) :: ncol
    !
    ! !REVISION HISTORY
    ! Created by Sam Levis
    ! Modified for CAM by Natalie Mahowald
    ! Modified for Seasalts by NMM
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k,ix,lchnk          !indices
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(pcols,pver,nsst) ![frc] Slip correction factor
    real(r8) :: rss_lmn(nsst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp ,r          !temporary 
    real(r8) :: wetdia(pcols,pver,nsst)        ! wet diameter of seasalts
    real(r8) :: RH(pcols,pver),es(pcols,pver),qs(pcols,pver)  ! for wet radius calculation

    ! constants

!wtw    real(r8),parameter::shm_nbr_xpn_lnd=-1.5            ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over ccean
    real(r8),parameter:: c1=0.7674_r8, c2=3.0790_r8, c3=2.57e-11_r8,c4=-1.424_r8  ! wet radius calculation constants


    ! needs fv and ram1 passed in from lnd model

    !------------------------------------------------------------------------


    call aqsat(t,pmid,es,qs,pcols,ncol,pver,1,pver)

    RH(:ncol,:)=q(:ncol,:,1)/qs(:ncol,:)
    RH(:ncol,:)=max(0.01_r8,min(0.99_r8,RH(:ncol,:)))
    ! set stokes correction to 1.0 for now not a bad assumption for our size range)
    do m=1,nsst
       stk_crc(m)=1._r8
    enddo
    do k=1,pver
       do i=1,ncol
          rho(i,k)=pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho(i,k) ![m2 s-1] Kinematic viscosity of air

          do m = 1, nsst
             r=smt_vwr(m)/2.0_r8
             wetdia(i,k,m)=((r**3+c1*r**c2/(c3*r**c4-log(RH(i,k))))**(1._r8/3._r8))*2.0_r8
             slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*wetdia(i,k,m)/(2.0_r8*mfp_atm(i,k)))) / &
                  wetdia(i,k,m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * wetdia(i,k,m) * wetdia(i,k,m) * dns_aer_sst * &
                  gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m)=vlc_grv(i,k,m)
          end do

       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do m = 1, nsst
       do i=1,ncol
          r=smt_vwr(m)/2.0_r8
          wetdia(i,k,m)=((r**3+c1*r**c2/(c3*r**c4-log(RH(i,k))))**(1._r8/3._r8))*2.0_r8

          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
               (3.0_r8*pi*vsc_dyn_atm(i,k)*wetdia(i,k,m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
          shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          !           if(ocnfrac.gt.0.5) shm_nbr_xpn=shm_nbr_xpn_ocn
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965

          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do !ncol
    end do
    return
    end subroutine ProgseasaltsDryDep

end module gchp_bcc_progseasalts_intr
