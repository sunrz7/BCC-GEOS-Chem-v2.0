!xiaolu
!---------------------------------------
!Purpose: AEROSOL DEPOSITION
!from other_aerosols_drydep.F90

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/15
!---------------------------------------

#include <misc.h>
#include <params.h>

#if (defined MOZART2) 
    subroutine gchp_other_aerosols_drydep &
        ( lchnk, ncol,  q, pdel, t, pmid, lq, dq, dtime, &
                                 obklen, ustar, pblh, &
                                 landfrac, icefrac, ocnfrac,fvin,ram1in , vlc_trb_SO2)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to dry deposition
    ! 
    ! Author: WU Tongwen 
    ! 
    !-----------------------------------------------------------------------
    use shr_kind_mod,          only: r8 => shr_kind_r8
    use constituents,          only: ppcnst, cnst_get_ind
    use ppgrid,                only: pcols, pver
    use gchp_bcc_drydep_mod,        only: d3ddflux, calcram
    use physconst,             only: rair,pi, gravit, boltz
    use history,               only: outfld

    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------
    ! 	... parameters for log-normal distribution by number
    ! references:
    !   Chin et al., JAS, 59, 461, 2003
    !   Liao et al., JGR, 108(D1), 4001, 2003
    !   Martin et al., JGR, 108(D3), 4097, 2003
    !-----------------------------------------------------------------

!---
! wtw  refer to CAM5
!
    real(r8), parameter :: rm_sulf  = 0.05e-6           ! mean radius of sulfate particles (m) (wutw)
    real(r8), parameter :: sd_sulf  = 2.03_r8           ! standard deviation of radius for sulfate (Chin)
    real(r8), parameter :: rho_sulf = 1.77e3_r8         ! density of sulfate aerosols (kg/m3) (wutw) 
    real(r8), parameter :: dispersion_sulf = 2.03_r8    ! geometric standard deviation of
                                                        ! aerosol size distribution

    real(r8), parameter :: rm_orgc  = 0.03e-6           ! mean radius of organic carbon particles (m) (wutw)
    real(r8), parameter :: sd_orgc  = 2.24_r8           ! standard deviation of radius for OC (Chin)
    real(r8), parameter :: rho_orgc = 1.8e3_r8          ! density of OC aerosols (kg/m3) (Chin)
    real(r8), parameter :: dispersion_orgc = 2.24_r8    ! geometric standard deviation of 
                                                        ! aerosol size distribution

    real(r8), parameter :: rm_bc    = 0.03e-6           ! mean radius of soot/BC particles (m) (Chin)
    real(r8), parameter :: sd_bc    = 2.00_r8           ! standard deviation of radius for BC (Chin)
    real(r8), parameter :: rho_bc   = 1.7e3_r8          ! density of BC aerosols (kg/m3) (Chin)
    real(r8), parameter :: dispersion_bc = 2._r8       ! geometric standard deviation of    
                                                        ! aerosol size distribution

    !
    ! idxSUL, ixdSSLT, idxdst1, idxdst2, idxdst3, idxdst4, idxOCPHO, idxBCPHO, idxOCPHI, idxBCPHI, idxBG, idxVOLC
    !  real(r8) :: dispersion_aer (naer_all)        ! in,geometric standard deviation of aerosol size distribution
    !  real(r8) :: density_aer    (naer_all)        ! in,density of aerosol (kg/m3)
    !  data  dispersion_aer /2.03, 1.59138, 1.9, 1.9, 1.9, 1.9, 2.24, 2, 2.24, 2, 2, 2.03/
    !  data  density_aer    /1170., 2200., 2600., 2600., 2600., 2600., 1800., 1000., &
    !                        2600., 1000., 1770., 1770./
    !
    ! Arguments:
    !
    integer, intent(in) :: lchnk                              ! chunk identifier
    integer, intent(in) :: ncol                               ! number of atmospheric columns
    real(r8),            intent(in)  :: dtime                 ! time step

    logical, intent(inout) :: lq(ppcnst)                   !ptend%lq
    real(r8), intent(inout) :: dq(pcols,pver,ppcnst)       !ptend%q

    real(r8), intent(in) :: q(pcols,pver,ppcnst)           !state%q
    real(r8), intent(in) :: t(pcols,pver)                  !state%t
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)

    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel--used over oceans and sea ice.
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    real(r8), intent(in) :: fvin(pcols)        ! for dry dep velocities from land model for dust
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model for dust
    real(r8),intent(in)  :: vlc_trb_SO2(pcols)       ! dry deposite velcity (m/s)

    !
    ! Local variables
    !
    integer :: m                                  ! tracer index
    real(r8) :: tvs(pcols,pver)
    real(r8) :: sflx(pcols)                    ! deposition flux
    real(r8) :: vlc_dry(pcols,pver)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver)            ! dep velocity
    real(r8)::  vlc_trb(pcols)                 ! dep velocity
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice

    real(r8) :: vsc_dyn_atm(pcols,pver)

	integer so4_idx, oc1_idx, oc2_idx, cb1_idx, cb2_idx, so2_idx
	real(r8) :: radius, dispersion, density
	real(r8) :: rho        (pcols,pver)
	real(r8) :: mfp_atm    (pcols,pver)
	real(r8) :: vsc_knm_atm(pcols,pver)
	real(r8) :: slp_crc    (pcols,pver)

        integer, parameter :: moment = 3 ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    real(r8) :: radius_moment
    real(r8) :: lnsig, sig_part
    real(r8) :: stk_nbr, dff_aer, shm_nbr, shm_nbr_xpn, tmp,  rss_lmn, rss_trb 
    integer :: i,k
    real(r8) :: dep_trb(pcols), dep_grv(pcols) 

    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8          ![frc] shm_nbr_xpn over land
!org    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
!    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over land
    !
    !-----------------------------------------------------------------------

    tvs(:ncol,:) = t(:ncol,:) 
	
    !--------------------------------------------------------------------------
    ! copy fv,ram1 values from land model to variables to modify in calcram
    !    fv = fvin    
    !    ram1 = ram1in

    !-----------------------------------------------------------------------
    !Calc aerodynamic resistance over oceans and sea ice (comes in from land
    !model) ,xiaolu
    !-----------------------------------------------------------------------


    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
                 ustar,ram1in, ram1, t(:,pver),pmid(:,pver),&
                 pdel(:,pver),fvin,fv)
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
!xiaolu change to match GEOSCHEM
!    call cnst_get_ind('SO4', so4_idx)
!    call cnst_get_ind('OC1', oc1_idx)
!    call cnst_get_ind('OC2', oc2_idx)
!    call cnst_get_ind('CB1', cb1_idx)
!    call cnst_get_ind('CB2', cb2_idx)
!    call cnst_get_ind('SO2', so2_idx)

    call cnst_get_ind('SO4', so4_idx)
    call cnst_get_ind('OCPO', oc1_idx)
    call cnst_get_ind('OCPI', oc2_idx)
    call cnst_get_ind('BCPO', cb1_idx)
    call cnst_get_ind('BCPI', cb2_idx)
    call cnst_get_ind('SO2', so2_idx)

    tvs(:ncol,:) = t(:ncol,:)   !*(1+q(:ncol,k)
	
    do m=1, ppcnst 

   	  vlc_dry(:,:) = 0.0
	  vlc_grv(:,:) = 0.0

      if( m == so4_idx ) then     
          radius     = rm_sulf 
          sig_part   = sd_sulf                     
          density    = rho_sulf    
	  goto 111
		                  !
      elseif ( m == oc1_idx .or.  m == oc2_idx )  then
          radius     = rm_orgc 
          sig_part   = sd_orgc                     
          density    = rho_orgc
	  goto 111

      elseif ( m == cb1_idx .or. m == cb2_idx ) then 
          radius     = rm_bc 
          sig_part   = sd_bc                     
          density    = rho_bc 
          goto 111

      elseif ( m == so2_idx ) then
          goto 333

      else
	  goto 222
      endif
!-------------------
!-------------------
!For OC,BC,sulfate particles
111   continue
      do i=1,ncol
         do k=1,pver

          lnsig = log(sig_part)

! use a maximum radius of 50 microns when calculating deposition velocity
          radius_moment = radius*   &  
                          exp((float(moment)-1.5)*lnsig*lnsig)
          dispersion = exp(2.*lnsig*lnsig)

          rho(i,k) =  pmid(i,k)/(rair*t(i,k))

          !==================================================
          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties

          !-----------
          !Despription: 
          !following is to calculate QUASI-LAMINAR RESISTANCE for PARTICLE
          !based on SeP, Chapter 19, dry deposition
          !xiaolu,2017/04/25
          !
          !PART1: CALCULATE PARTICLE SETTLING VELOCITY
          !vsc_dyn_atm:dynamic viscosity
          !vsc_knm_atm:kinematic viscosity
          !mfp_atm:mean free pass of air
          !slp_crc:slip correction factor for Stoke's Law, SeP Chap. 9
          !vlc_grv:particle settling velocity

          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
                             (t(i,k)+120.0_r8)                ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k)     = 2.0_r8 * vsc_dyn_atm(i,k) / &                 ![m] SeP97 p. 455
                             (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))

          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho(i,k)       ![m2 s-1] Kinematic viscosity of air
          slp_crc(i,k)     = 1.0_r8 + mfp_atm(i,k) * &
                             (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment/(mfp_atm(i,k)))) / &
                             radius_moment             ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(i,k) = (4.0_r8/18.0_r8) * radius_moment*radius_moment*density* &
                         gravit*slp_crc(i,k) / vsc_dyn_atm(i,k)                 ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(i,k) = vlc_grv(i,k) * dispersion

          vlc_dry(i,k) = vlc_grv(i,k)
       enddo
     enddo

        !----------------
        !Description:
        !PART2: CALCULATE rb, refer to SeP, Chap. 19
        !xiaolu,2017/04/25
        !
        !stk_nbr:Stokes number to calculate Eim(colletcion efficiency of
        !impaction)
        !dff_aer:Brownian diffusivity of particle
        !shm_nbr:Schmidt number
        !shm_nbr_xpn:[1/3,2/3], Eb=shm_nbr**shm_nbr_xpn 
        
 
     k=pver  ! only look at bottom level for next part
     do i=1,ncol

          stk_nbr = vlc_grv(i,k) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k) / &                 ![m2 s-1]
                      (3.0_r8*pi*vsc_dyn_atm(i,k)* radius_moment )           !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
                            !          shm_nbr_xpn=shm_nbr_xpn_ocn
                            !          if(landfrac(i) .gt. 0.5_r8 ) shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          shm_nbr_xpn = shm_nbr_xpn_lnd 
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn = 1.0_r8 / (tmp*fv(i))                                ![s m-1] SeP97 p.972,965

          rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i,k) ![s m-1]
          vlc_trb(i) = 1.0_r8 / rss_trb                                 ![m s-1]
          vlc_dry(i,k) = vlc_trb(i)  +vlc_grv(i,k)

       end do !end ncols loop
  
       sflx(:) = 0.0
       call d3ddflux(ncol, vlc_dry(:,:), q(:,:,m), pmid, pdel, tvs, sflx, dq(:,:,m),dtime)

!-------------------------2015-06-01---------------------------------------------
       sflx(:) = 0.
       do k=1, pver
          do i=1, ncol
             sflx(i) = sflx(i) + dq(i,k,m)*pdel(i,k)/gravit
          enddo
       enddo

       ! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i)/vlc_dry(i,pver)                 !   kg /m^s/sec
          dep_grv(i)=sflx(i)*vlc_grv(i,pver)/vlc_dry(i,pver)            !   kg /m^s/sec
       enddo

!       if( m == so4_idx ) call outfld('DRY_SO4 ', sflx, pcols, lchnk)
       if( m == oc1_idx ) call outfld('DRY_OC1 ', sflx, pcols, lchnk)
!       if( m == oc2_idx ) call outfld('DRY_OC2 ', sflx, pcols, lchnk)
!       if( m == cb1_idx ) call outfld('DRY_CB1 ', sflx, pcols, lchnk)
!       if( m == cb2_idx ) call outfld('DRY_CB2 ', sflx, pcols, lchnk)
        
!       if( m == so4_idx ) call outfld('GRV_SO4 ', dep_grv, pcols, lchnk)
!       if( m == oc1_idx ) call outfld('GRV_OC1 ', dep_grv, pcols, lchnk)
!       if( m == oc2_idx ) call outfld('GRV_OC2 ', dep_grv, pcols, lchnk)
!       if( m == cb1_idx ) call outfld('GRV_CB1 ', dep_grv, pcols, lchnk)
!       if( m == cb2_idx ) call outfld('GRV_CB2 ', dep_grv, pcols, lchnk)
!-------------------------------------------------------------------------------

       lq(m) = .TRUE.

       goto 222
!-------------------------------------------------------------------------------------
 333   continue
!
! SO2 dry deposition
!

     vlc_dry = 0.

     k=pver  ! only look at bottom level for next part
     do i=1,ncol
        ! -- ocean dry deposition velocity
          rss_trb = ram1(i) ![s m-1]
          vlc_trb(i) = 1.0_r8 / rss_trb                                 ![m s-1]
          vlc_trb(i) = vlc_trb_SO2(i) + vlc_trb(i) * (1.-landfrac(i))

          vlc_dry(i,pver) = vlc_trb(i)*0.2 
     enddo 
     sflx(:) = 0.0
     call d3ddflux(ncol, vlc_dry(:,:), q(:,:,m), pmid, pdel, tvs, sflx, dq(:,:,m),dtime)

     dep_trb(:) = 0.
     do k=1, pver
        do i=1, ncol
           dep_trb(i) = dep_trb(i) + dq(i,k,m)*pdel(i,k)/gravit
        enddo
     enddo
!     call outfld('DRY_SO2 ', dep_trb, pcols, lchnk)

     lq(m) = .TRUE.

222 continue 
    end do
    return
    end subroutine gchp_other_aerosols_drydep
#endif
