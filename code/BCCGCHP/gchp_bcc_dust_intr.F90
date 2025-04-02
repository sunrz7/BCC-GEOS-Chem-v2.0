!xiaolu
!---------------------------------------
!Purpose: DUST AEROSOL DRY DEPOSITION and EMISSION
!from bcc_dust_intr.F90 

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/15
!---------------------------------------

#include <misc.h> 
#include <params.h>
module gchp_bcc_dust_intr

  !---------------------------------------------------------------------------------
  ! Module to interface the aerosol parameterizations with CAM
  ! written by PJR (extensively modified from chemistry module)
  !---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8, cl => shr_kind_cl
  use pmgrid,        only: masterproc
  use ppgrid,        only: pcols, pver,pverp
  use physconst,     only: mwdry, mwh2o,gravit,rair
  use constituents,  only: ppcnst, cnst_add, cnst_name, cnst_get_ind
  use abortutils,    only: endrun
  use gchp_bcc_chem_mods, only: ndst ! mo_->gchp_
  use ppgrid,        only: begchunk, endchunk
  use gchp_chem_utls,  only: get_spc_ndx!mo_ -> gchp_
  use phys_grid,     only: scatter_field_to_chunk
  use mpishorthand,  only: mpicom, mpiint, mpir8

  implicit none

  private          ! Make default type private to the module

  save

  real(r8), allocatable ::  soil_erodibility(:,:)     ! soil erodibility factor

  integer, private  :: ncid

  integer, parameter:: dst_src_nbr =3
  integer, parameter:: sz_nbr =200

!  character(len=8), dimension(ncnst), parameter :: & ! constituent names
!       cnst_names = (/'DST01', 'DST02', 'DST03', 'DST04'/)

  !
  ! Public interfaces
  !
  public dust_initialize                           ! initialize (history) variables
!  public bcc_dust_wet                                  ! interface to wet deposition
  public bcc_dust_emis_intr                            ! interface to emission
  public dust_drydep_intr                          ! interface to tendency computation

!----------------------------

  real(r8) ::stk_crc(ndst)             ![frc] Correction to Stokes settling velocity
  real(r8) ::dns_aer                                    ![kg m-3] Aerosol density
  real(r8) ::ovr_src_snk_mss(dst_src_nbr,ndst)  
  real(r8) ::dmt_vwr(ndst)            ![m] Mass-weighted mean diameter resolved

  real(r8) ::tmp1
  
! Namelist variables
  real(r8)      :: dust_emis_fact = -1.e36   ! tuning parameter for dust emissions

!===============================================================================
contains
!===============================================================================

subroutine dust_initialize(soil_erod_file)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: initialize parameterization of dust chemistry
   !          (declare history variables)
   ! 
   !-----------------------------------------------------------------------

   use ioFileMod,        only : getfil
   use phys_grid,        only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
   use gchp_const_mozart,  only : pi, d2r

   character(len=*), intent(in)  :: soil_erod_file

   ! local variables
   integer :: did, vid, nlat, nlon
   integer :: m
   integer :: c 

   real(r8), allocatable :: erodibility_in(:,:)  ! temporary input array
   real(r8), allocatable :: dst_lons(:)
   real(r8), allocatable :: dst_lats(:)

   character(len=256)    :: locfn
   !-----------------------------------------------------------------------

   ! use Sam's dust initialize subroutine:  call equivalent here:

   call Dustini()

   ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
   ! read in soil erodibility factors, similar to Zender's boundary conditions
!
! wtw  soil_erod: atm/cam/dst/dst_64x128_c090203.nc
!
   ! Get file name.  

   if( masterproc ) then

      call getfil(soil_erod_file, locfn)
      call wrap_open(locfn, 0, ncid)

      call wrap_inq_dimid( ncid, 'lon', did )
      call wrap_inq_dimlen( ncid, did, nlon )

      call wrap_inq_dimid( ncid, 'lat', did )
      call wrap_inq_dimlen( ncid, did, nlat )

   endif

#if (defined SPMD )
   call mpibcast (nlon, 1, mpiint, 0, mpicom)
   call mpibcast (nlat, 1, mpiint, 0, mpicom)
#endif

   allocate(dst_lons(nlon))
   allocate(dst_lats(nlat))
   allocate(erodibility_in(nlon,nlat))
   allocate(soil_erodibility(pcols,begchunk:endchunk))

   if( masterproc ) then
       call wrap_inq_varid( ncid, 'lon', vid )
       call wrap_get_var_realx( ncid, vid, dst_lons  )

       call wrap_inq_varid( ncid, 'lat', vid )
       call wrap_get_var_realx( ncid, vid, dst_lats  )
   end if

   if( masterproc ) then
       call wrap_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
       call wrap_get_var_realx( ncid, vid, erodibility_in )
   endif

   call scatter_field_to_chunk(1,1,1,nlon, erodibility_in, soil_erodibility(1,begchunk))

   deallocate(erodibility_in)
   
    return
    end subroutine dust_initialize

!===============================================================================

    subroutine dust_drydep_intr (lchnk, ncol,  q, pdel, t, pmid, lq, dq, dt, lat, clat, &
                                 fsds, obklen, ts, ustar, pblh, &
		  	         month, landfrac, icefrac, ocnfrac,fvin,ram1in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to dry deposition and sedimentation of dust
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
    use gchp_bcc_drydep_mod,        only: d3ddflux, calcram  !mo_ -> gchp_
    use gchp_bcc_dust_sediment_mod, only: dust_sediment_tend !mo_ -> gchp_
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    integer, intent(in) :: lchnk                              ! chunk identifier
    integer, intent(in) :: ncol                               ! number of atmospheric columns
    real(r8),            intent(in)  :: dt             ! time step

    logical, intent(inout) :: lq(ppcnst)
    real(r8), intent(inout) :: dq(pcols,pver,ppcnst)

    real(r8), intent(in) :: q(pcols,pver,ppcnst)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)

    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
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
    real(r8), intent(in) :: fvin(pcols)        ! for dry dep velocities from land model for dust
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model for dust


    ! Local variables
    !
    integer :: m                                  ! tracer index
    integer :: mm                                  ! tracer index
    real(r8) :: tvs(pcols,pver)
    real(r8) :: dvel(pcols)            ! deposition velocity
    real(r8) :: sflx(pcols)            ! deposition flux
    real(r8) :: vlc_dry(pcols,pver,ndst)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,ndst)            ! dep velocity
    real(r8)::  vlc_trb(pcols,ndst)            ! dep velocity
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_trb(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_dry_tend(pcols,pver)       !kg/kg/s (total of grav and trb)
    real(r8) :: obuf(1)
    real(r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real(r8) :: pvdust(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: tsflx(pcols)
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice

    integer  :: i,k
    real(r8) :: oro(pcols)
	integer  :: dust_idx1
    !
    !-----------------------------------------------------------------------

    tvs(:ncol,:) = t(:ncol,:)     !*(1+q(:ncol,k)
    rho(:ncol,:)=  pmid(:ncol,:)/(rair*t(:ncol,:))
    tsflx(:)=0._r8
    ! calculate oro--need to run match


    do i=1,ncol
       oro(i)=1._r8
       if(icefrac(i)>0.5_r8) oro(i)=2._r8
       if(ocnfrac(i)>0.5_r8) oro(i)=0._r8
    enddo

    !   Dry deposition of Dust Aerosols
    !  we get the ram1,fv from the land model as ram1in,fvin, but need to calculate it over oceans and ice.  
    !  better if we got thse from the ocean and ice model
    !  for friction velocity, we use ustar (from taux and tauy), except over land, where we use fv from the land model.

    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
                 ustar,ram1in,ram1,t(:,pver),pmid(:,pver),&
                 pdel(:,pver),fvin,fv)

    call DustDryDep(ncol,t(:,:),pmid(:,:),ram1,fv,vlc_dry,vlc_trb,vlc_grv,landfrac)
  
!--------------------------------
!xiaolu,change DST01->DST1 to match geoschem
!--------------------------------
!    call cnst_get_ind('DST01', dust_idx1)
    call cnst_get_ind('DST1', dust_idx1)   
 
    do m=1,ndst
       mm = dust_idx1 + m - 1

       ! use pvdust instead (means making the top level 0)
       pvdust(:ncol,1)=0._r8
       pvdust(:ncol,2:pverp) = vlc_dry(:ncol,:,m)

       call d3ddflux(ncol, vlc_dry(:,:,m), q(:,:,mm),pmid,pdel, tvs,sflx, dq(:,:,mm),dt)

       ! apportion dry deposition into turb and gravitational settling for tapes

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

!-----------------------2015-05-29-------------------------------------

       if( m == 1) call outfld( 'DRY_DST1', sflx, pcols, lchnk)
!       if( m == 2) call outfld( 'DRY_DST2', sflx, pcols, lchnk)
!       if( m == 3) call outfld( 'DRY_DST3', sflx, pcols, lchnk)
!       if( m == 4) call outfld( 'DRY_DST4', sflx, pcols, lchnk)

!       if( m == 1) call outfld( 'GRV_DST1', dep_grv, pcols, lchnk)
!       if( m == 2) call outfld( 'GRV_DST2', dep_grv, pcols, lchnk)
!       if( m == 3) call outfld( 'GRV_DST3', dep_grv, pcols, lchnk)
!       if( m == 4) call outfld( 'GRV_DST4', dep_grv, pcols, lchnk)
!---------------------------------------------------------------------

       lq(mm) = .TRUE.
    end do
    return
    end subroutine dust_drydep_intr


    subroutine bcc_dust_emis_intr (lchnk, ncol, dt, sflx, landfrac, pdel, &
                     flx_mss_vrt_dst01, flx_mss_vrt_dst02, flx_mss_vrt_dst03, flx_mss_vrt_dst04)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model (need #define BGC) and
    ! the soil erodibility factor is applied here.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: Phil Rasch and Natalie Mahowald
    !
    ! 
    !-----------------------------------------------------------------------
    use history,       only: outfld
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !

    integer, intent(in)                :: lchnk
    integer, intent(in)                :: ncol
    real(r8),            intent(in)    :: dt             ! time step
    real(r8),            intent(inout) :: sflx(pcols, ppcnst)   !  ! surface emissions ( kg/m^2/s )
    real(r8),            intent(in)    :: landfrac(pcols)
    real(r8),            intent(in)    :: pdel(pcols,pver)
    real(r8), intent(in)  :: flx_mss_vrt_dst01(pcols)    !surface dust emission (kg/m**2/s) [ + = to atm]    
    real(r8), intent(in)  :: flx_mss_vrt_dst02(pcols)    !surface dust emission (kg/m**2/s) [ + = to atm]    
    real(r8), intent(in)  :: flx_mss_vrt_dst03(pcols)    !surface dust emission (kg/m**2/s) [ + = to atm]    
    real(r8), intent(in)  :: flx_mss_vrt_dst04(pcols)    !surface dust emission (kg/m**2/s) [ + = to atm]    


    integer  :: i
    integer  :: m
    real(r8) :: soil_erod_tmp(pcols)
    !
    real(r8) :: fact  ! tuning factor for dust emissions
    real(r8) :: to_lats(pcols), to_lons(pcols)
    real(r8) :: cflx(pcols, ppcnst)
    integer  :: dust_idx1
   
    real(r8) :: emis_dust(pcols)
!--------------------------------------------------------------------
    !--------------------------------------------------
    !xiaolu,change DST01 to DST1 to match GEOSCHEM
    !--------------------------------------------------
    !call cnst_get_ind('DST01', dust_idx1)
     call cnst_get_ind('DST1', dust_idx1)

    do i = 1, ncol
       soil_erod_tmp(i) = soil_erodibility( i,lchnk )

       ! change test to from 0.1 to 0.001 after
       ! discussion with N. Mahowald April 8 2009
       if(soil_erod_tmp(i) .lt. 0.001_r8) soil_erod_tmp(i)=0._r8
    end do

    !cflx(:ncol,dust_idx1)  =flx_mss_vrt_dst01(:ncol)*0.038_r8/0.032456_r8
    !cflx(:ncol,dust_idx1+1)=flx_mss_vrt_dst02(:ncol)*0.11_r8/0.174216_r8
    !cflx(:ncol,dust_idx1+2)=flx_mss_vrt_dst03(:ncol)*0.17_r8/0.4085517_r8
    !cflx(:ncol,dust_idx1+3)=flx_mss_vrt_dst04(:ncol)*0.67_r8/0.384811_r8

    !xiaolu note 2018/08:
    !as GEOS-Chem and BCC-AVIM has different bin description of dust, we scale a factor in the emission
    !to reconcile

    cflx(:ncol,dust_idx1)  =flx_mss_vrt_dst01(:ncol)*0.038_r8/0.032456_r8 *1.5*0.9
    cflx(:ncol,dust_idx1+1)=flx_mss_vrt_dst02(:ncol)*0.11_r8/0.174216_r8 *2.5
    cflx(:ncol,dust_idx1+2)=flx_mss_vrt_dst03(:ncol)*0.17_r8/0.4085517_r8 *2*0.9
    cflx(:ncol,dust_idx1+3)=flx_mss_vrt_dst04(:ncol)*0.67_r8/0.384811_r8 *0.25*0.9*0.95



    ! tuning factor (fact) for dust emissions
    fact = 1.0 
!---------------------

    do m=dust_idx1,dust_idx1+3
       ! multiply by soil erodibility factor
       do i=1, ncol
          if( landfrac(i) .lt. 0.001) then
	      cflx(i,m) = 0.0
	  else
    !wtw      cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*1.15_r8
    !zf       cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*1.30_r8
!----------------
            if( m .eq. dust_idx1 ) then
              cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*0.24_r8
            else if ( m .eq. dust_idx1 + 1) then
              cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*0.23_r8
            else if ( m .eq. dust_idx1 + 2) then
              cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*0.22_r8
            else
              cflx(i,m) =cflx(i,m)*soil_erod_tmp(i)/fact*0.18_r8
            end if
!----------------
	  endif
          sflx(i,m) =sflx(i,m) + cflx(i,m)
       enddo
    enddo

!  write(*,*),'xiaolu check dust emissions:',cflx

!--------------------- 2015-05-29---------------------------
    call outfld('EMS_DST1',cflx(:, dust_idx1), pcols, lchnk)
    call outfld('EMS_DST2',cflx(:, dust_idx1+1), pcols, lchnk)
    call outfld('EMS_DST3',cflx(:, dust_idx1+2), pcols, lchnk)
    call outfld('EMS_DST4',cflx(:, dust_idx1+3), pcols, lchnk)
!----------------------2015-05-29---------------------------
    return
  end subroutine bcc_dust_emis_intr

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine Dustini()
  !
  ! !INTERFACE:
  !
  subroutine Dustini()
    !
    ! !DESCRIPTION: 
    !
    ! Compute source efficiency factor from topography
    ! Initialize other variables used in subroutine Dust:
    ! ovr_src_snk_mss(m,n) and tmp1.
    ! Define particle diameter and density needed by atm model
    ! as well as by dry dep model
    ! Source: Paul Ginoux (for source efficiency factor)
    ! Modifications by C. Zender and later by S. Levis
    ! Rest of subroutine from C. Zender's dust model
    !
    ! !USES
    !
    use physconst,    only: pi,rair
    !
    ! !ARGUMENTS:
    !
    implicit none
    !
    ! !REVISION HISTORY
    ! Created by Samual Levis
    ! Revised for CAM by Natalie Mahowald
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !Local Variables
    integer  :: ci,m,n                  !indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ![frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ![frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ![frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 !Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 !Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ![m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ![m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ![m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ![m] Width of size bin
    real(r8) :: slp_crc(ndst)           ![frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ![m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ![m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ![m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ![frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ![frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     !temporary 
    real(r8) :: ln_gsd                  ![frc] ln(gsd)
    real(r8) :: gsd_anl                 ![frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ![m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ![m] Number median particle diameter
    real(r8) :: lgn_dst                 !Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ![frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ![frc] Current relative accuracy
    real(r8) :: itr_idx                 ![idx] Counting index
    real(r8) :: dns_mdp                 ![kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ![m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ![kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ![m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            !Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      !Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ![m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ![m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ![m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ![m] Size Bin widths

    ! constants
    real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ![m] Mass median diameter
         (/ 0.832e-6_r8 , 4.82e-6_r8 , 19.38e-6_r8 /)        !BSM96 p. 73 Table 2
    real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ![frc] Geometric std deviation
         (/ 2.10_r8     ,  1.90_r8   , 1.60_r8     /)        !BSM96 p. 73 Table 2
    real(r8) :: mss_frc_src(dst_src_nbr) =    &     ![frc] Mass fraction 
         (/ 0.036_r8, 0.957_r8, 0.007_r8 /)                  !BSM96 p. 73 Table 2

    real(r8) :: dmt_grd(5) =                  &     ![m] Particle diameter grid
         (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6_r8    ![m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0_r8         ![kg m-3] Density of optimal saltation particles


    ! declare erf intrinsic function
    real(r8) :: dum     !dummy variable for erf test

    real(r8) derf
    !------------------------------------------------------------------------

    ! Sanity check on erf: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

    dum = 1.0_r8
    if (abs(0.8427_r8-ERF(dum))/0.8427_r8>0.001_r8) then
       write(*,*) 'erf(1.0) = ',ERF(dum)
       call endrun ('Dustini: Error function error')
    endif
    dum = 0.0_r8
    if (ERF(dum) /= 0.0_r8) then
       write(*,*) 'erf(0.0) = ',ERF(dum)
       call endrun ('Dustini: Error function error')
    endif

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
       sqrt2lngsdi = sqrt(2.0_r8) * log(gsd_anl_src(m))
       do n = 1, ndst
          lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
          lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
          ovr_src_snk_frc = 0.5_r8 * (ERF(lndmaxjovrdmdni/sqrt2lngsdi) - &
               ERF(lndminjovrdmdni/sqrt2lngsdi))
          ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
       enddo
    enddo

!-------------------------
! wtw
!
    ! The following code from subroutine wnd_frc_thr_slt_get was placed 
    ! here because tmp1 needs to be defined just once

    ryn_nbr_frc_thr_prx_opt = 0.38_r8 + 1331.0_r8 * (100.0_r8*dmt_slt_opt)**1.56_r8

    if (ryn_nbr_frc_thr_prx_opt < 0.03_r8) then
       write(6,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
       call endrun
    else if (ryn_nbr_frc_thr_prx_opt < 10.0_r8) then
       ryn_nbr_frc_thr_opt_fnc = -1.0_r8 + 1.928_r8 * (ryn_nbr_frc_thr_prx_opt**0.0922_r8)
       ryn_nbr_frc_thr_opt_fnc = 0.1291_r8 * 0.1291_r8 / ryn_nbr_frc_thr_opt_fnc
    else
       ryn_nbr_frc_thr_opt_fnc = 1.0_r8 - 0.0858_r8 * exp(-0.0617_r8*(ryn_nbr_frc_thr_prx_opt-10.0_r8))
       ryn_nbr_frc_thr_opt_fnc = 0.120_r8 * 0.120_r8 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
    end if

    icf_fct = 1.0_r8 + 6.0e-07_r8 / (dns_slt * gravit * (dmt_slt_opt**2.5_r8))
    dns_fct = dns_slt * gravit * dmt_slt_opt
    tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

! wtw
!------------------------

    !delete portions that have to do only with the sink

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here ndst = dst_nbr

    ! Override automatic grid with preset grid if available

    if (ndst == 4) then
       do n = 1, ndst
          dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
          dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
          dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
          dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
       end do
    else
       call endrun('Dustini error: ndst must equal to 4 with current code')  
    endif

    ! Bin physical properties
    gsd_anl = 2.0_r8      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e+3_r8   ! [kg m-3] Aerosol density
    ! Set a fundamental statistic for each bin
    dmt_vma = 2.524e-6_r8 ! [m] Mass median diameter analytic She84 p.75 Table1
    dmt_vma = 3.5e-6_r8
    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)
    dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]
    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl
    do n = 1, ndst
       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_r8/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do
       lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*pi))
       dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved
       do m = 1, sz_nbr
          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)
          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               pi / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               pi / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
       end do
       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do
    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)
    eps_max = 1.0e-4_r8
    dns_mdp = 100000._r8 / (295.0_r8*rair) ![kg m-3] const prs_mdp & tpt_vrt
    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_r8 * ((295.0_r8/273.0_r8)**1.5_r8) * 393.0_r8 / &
         (295.0_r8+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000._r8*sqrt(8.0_r8/(pi*rair*295.0_r8)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, ndst
       slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
            (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
            gravit * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, ndst

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]
       do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_r8) then
             cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                  (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                  9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                  (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_r8) then
             cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
          else
             write(*,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
             call endrun ('Dustini error: Reynolds number too large in stk_crc_get()')
          endif

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)
          vlc_grv(m) = sqrt(4.0_r8 * gravit * dmt_vwr(m) * slp_crc(m) * dns_aer / &
               (3.0_r8*cff_drg_grv(m)*dns_mdp))   
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          endif
          if (itr_idx > 20) then
             write(*,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
             goto 100                                        !to next iteration
          endif
          itr_idx = itr_idx + 1

       end do                                                !end while
100    continue   !Label to jump to when iteration does not converge
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities
    do m = 1, ndst
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

    return
    end subroutine Dustini

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine DustDryDep(c)
  !
  ! !INTERFACE:
  !
  subroutine DustDryDep(ncol,t,pmid,ram1,fv,vlc_dry,vlc_trb,vlc_grv,landfrac)
    !
    ! !DESCRIPTION: 
    !
    ! Determine Turbulent dry deposition for dust. Calculate the turbulent 
    ! component of dust dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the dust dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, dustini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
    !
    ! !USES
    !
    use physconst,     only: rair,pi,boltz
    !
    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8) :: t(pcols,pver)       !atm temperature (K)
    real(r8) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8),intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8),intent(in) :: ram1(pcols)           
                                               !    real(r8) :: ramxo1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: vlc_trb(pcols,ndst)  !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,ndst)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,ndst)  !dry deposn velocity (m/s)
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    integer, intent(in) :: ncol
    !
    ! !REVISION HISTORY
    ! Created by Sam Levis
    ! Modified for CAM by Natalie Mahowald
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k          !indices
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(pcols,pver,ndst) ![frc] Slip correction factor
    real(r8) :: rss_lmn(ndst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp           !temporary 

    ! constants

!wtw    real(r8),parameter::shm_nbr_xpn_lnd=-1.5       ! ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over land
    ! needs fv and ram1 passed in from lnd model

    !------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol
          rho=pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          do m = 1, ndst
             slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm(i,k)))) / &
                  dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                  gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m) = vlc_grv(i,k,m)
          end do
       enddo
    enddo

    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       do m = 1, ndst
          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
               (3.0_r8*pi*vsc_dyn_atm(i,k)*dmt_vwr(m)) !SeP97 p.474
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
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965
       end do

       ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)
       do m = 1, ndst
          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do

    end do !end ncols loop

    return
  end subroutine DustDryDep

end module gchp_bcc_dust_intr
