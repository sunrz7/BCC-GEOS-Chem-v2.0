!-----------------------------------------------------------------------
! Purpose:
! CALL GEOS-Chem in BCC
!
! REVISION HISTORY:
! Xiao Lu,2017/01/11: initial version
! Xiao Lu,2017/06/01: call gas deposition
! Xiao Lu,2017/06/07: add flux variables (sflx,flux_NNN)
! Xiao Lu,2017/06/13: add gas emission
! Xiao Lu,2017/10/05: change gigc_chunk_run to gigc_chunk_run_bc
! Xiao Lu,2017/10/08: add lightning NO emission
!-----------------------------------------------------------------------

#include <misc.h>
#include <params.h>

#if (defined MOZART2)

subroutine tphysbc_geoschem(ztodt,state,calday,State_Met,&
           State_Chm,Input_Opt,State_Grid,State_Diag,ts,landfrac,ocnfrac,icefrac,&
           flux_ISOP, flux_ACET,  flux_C3H6,                &
           flux_C2H4, flux_OC2,   flux_C10H16,              &
           flux_CO2,  flux_N2O,   flux_DMS,                 &
           sflx,cnt,cnb,                                    &
           GC2BCC,eflx)

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p,get_lat_all_p,get_lon_all_p
   use phys_buffer,     only: pbuf_size_max, pbuf_fld
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update
   use history,         only: outfld
   use time_manager,    only: get_step_size, dtime, get_nstep
   use constituents,    only: pcnst, pnats, ppcnst, qmin, iaero
   use check_energy,    only: check_energy_chng

#ifdef MOZART2
   use pmgrid,             only: masterproc,iam,plev,plon,plat
   use gigc_chunk_mod,     only: gigc_chunk_run
   USE input_opt_mod,           ONLY : optinput
   USE state_chm_mod,           ONLY : chmstate
   USE state_met_mod,           ONLY : metstate
   USE state_grid_mod,          ONLY : grdState!sunrz
   USE state_diag_mod           !sunrz
   USE hco_state_gc_mod,   ONLY : hcostate !sunrz 2024/3
   use time_manager,    only:get_curr_calday, get_curr_date
   USE cmn_size_mod
   USE gc_grid_mod,        ONLY :  setgridfromctr,setgridfromctredges
!sunrz add 2023/11
   USE roundoff_mod
   USE physconstants
!-------------
!FOR DRY DEPOSITION
!-------------
   USE gchp_surface, ONLY:gc_sflxdr
   USE chemistry,    ONLY: ncnst,ixchm
   USE Species_Mod,          ONLY : Species,  SpcConc
    use Unitconv_Mod,        only : Convert_Spc_Units
!------------
!FOR EMISSION
!------------

!   USE prescribed_em_d3,  ONLY: getem3_IPCC, nem3
!   USE prescribed_em_d2,  ONLY: getem2_IPCC, nem2, varname2d
!   USE prescribed_em_d2,  ONLY: idxCO, idxCH4, idxNO, idxC2H4
!   USE prescribed_em_d3,   ONLY: getem3_IPCC, nem3
!   USE prescribed_emis_3d, ONLY: getem3_IPCC6, nem3_IPCC6
!   USE prescribed_em_d2,   ONLY: getem2_IPCC, nem2, varname2d
!   USE prescribed_em_d2_anthro, ONLY: getem2_IPCC_anthro, nem2_anthro
!   USE prescribed_em_d2, ONLY: idxCO, idxCH4, idxNO, idxC2H4

   !3-D emission
!   USE gchp_bcc_em3_prod, ONLY:bcc_em3_prod

   !2019/01/14
!   USE gchp_read_tomsuv,   ONLY: get_tomsuv
   !2019/01/25
!   USE gchp_read_tomsuv,   ONLY: get_uvalbedo

   !2019/02/17 
!   USE gchp_bromine_emis,  ONLY: get_bromine_emis
!------------
!FOR LIGHTNING
!------------
!        USE gchp_hook,        ONLY: gc_hook !2017/10/08
        use mo_hook,         only: moz_hook    !2008.07.10
        use commap,          only: latdeg, londeg
        use shr_const_mod,  only: SHR_CONST_PI, SHR_CONST_REARTH!sunrz 2024/2
#endif

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   type(physics_state), intent(inout) :: state
   
#ifdef MOZART2
    TYPE(OptInput),intent(inout) :: Input_Opt   ! Input Options object
    TYPE(ChmState),intent(inout) :: State_Chm   ! Chemistry State object
    TYPE(MetState),intent(inout) :: State_Met   ! Meteorology State object
!sunrz add 2023/11
    TYPE(GrdState), INTENT(inout)    :: State_Grid     ! Grid State object
    TYPE(DgnState), INTENT(inout)    :: State_Diag
    !-------
    !FOR DRY DEPOSITION/EMISSION
    !-------
   REAL(r8), intent(in) :: ts(pcols)       ! surface temperature
   REAL(r8), intent(in) :: landfrac(pcols)       ! land fraction
    !XIAOLU,2017/06/06
   REAL(r8), intent(inout) :: sflx(pcols, ppcnst) !surface flux?
    real(r8),intent(in)    :: flux_ISOP(pcols)
    real(r8),intent(in)    :: flux_ACET(pcols)
    real(r8),intent(in)    :: flux_C3H6(pcols)
    real(r8),intent(in)    :: flux_C2H4(pcols)
    real(r8),intent(in)    :: flux_OC2(pcols)
    real(r8),intent(in)    :: flux_C10H16(pcols)
    real(r8),intent(in)    :: flux_CO2(pcols)
    real(r8),intent(in)    :: flux_N2O(pcols)
    real(r8),intent(in)    :: flux_DMS(pcols)
    !-------
    !FOR LIGHTNING NO EMISSION
    !-------
    REAL(r8), intent(in) :: ocnfrac(pcols)
    REAL(r8), intent(in) :: icefrac(pcols)
    real(r8),intent(in)  :: cnt(pcols)    ! Top level ofconvective activity
    real(r8),intent(in)  :: cnb(pcols)     ! Lowest level of convective activity

    !-------
    !for 3-D emission
    !------
    Integer,  intent(in) ::GC2BCC(ncnst)
#endif

!
!---------------------------Local workspace-----------------------------
!
!   type(physics_ptend),intent(inout)  :: ptend                  ! indivdual parameterization tendencies

   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)

   real(r8),intent(in)  :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   !real(r8) :: zero(pcols)                    ! array of zeros

   real(r8) :: prod_no(pcols,pver)
   real(r8) :: emilnox(pcols,pver)
#   if (defined MOZART2)
   integer ::IM
   integer ::JM
   integer ::LM
   integer ::i,j,l,n,II,NA   ! loop index
   integer ::nlat
   integer ::iam_bcc
!sunrz 2023/11
   integer ::ierr,my_id 
   integer :: siipar,sjjpar,sllpar

   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
   INTEGER         :: nymd        ! YYYY/MM/DD @ current time
   INTEGER         :: nhms        ! hh:mm:ss   @ current time
   INTEGER         :: dayOfYr     ! UTC day of year
   INTEGER         :: hour        ! UTC hour
   INTEGER         :: minute      ! UTC minute
   INTEGER         :: second      ! UTC second
   REAL*4          :: utc         ! UTC time [hrs]
   REAL*4          :: hElapsed    ! Elapsed hours


    INTEGER        :: Phase       ! Run phase (1 or 2)
    LOGICAL         :: IsChemTime  ! Time for chemistry?
    LOGICAL         :: IsRadTime   ! sunrz add 2023/11

   integer :: ngcols
   real(r8), dimension(:), allocatable :: clat_d
   real(r8), dimension(:), allocatable :: clon_d
   real(r8), dimension(:), allocatable :: area_d
   real(r8), dimension(:), allocatable ::area_d_out
   real(r8), dimension(:), allocatable :: area_all
   integer :: hdim1_d, hdim2_d


 ! First call?
    LOGICAL, SAVE                  :: FIRST = .TRUE.

    REAL*4,allocatable           :: lonCtr(:,:)   ! Lon ctrs [deg] from ESMF
    REAL*4,allocatable           :: latCtr(:,:)   ! Lat ctrs [deg] from ESMF
    REAL*4,allocatable           :: lonEdge(:,:)   ! Lon ctrs [deg] from ESMF
    REAL*4,allocatable           :: latEdge(:,:)   ! Lat ctrs [deg] from ESMF
        real(r8), parameter :: rad2deg = 180._r8/SHR_CONST_PI
        real(r8) :: up,dn !sunrz 2024/3
    INTEGER  :: RC            ! Success or failure
!-------------------------
!VARIABLES FOR EMISSION
!-------------------------
!        real(r8) :: em2data (pcols, nem2)
!        real(r8) :: em3data (pcols, pver, nem3)
!        real(r8) :: em2data_anthro (pcols, nem2_anthro)
!        real(r8) :: em3_IPCC6(pcols, pver, nem3_IPCC6)

!--------------------------
!VARIABLES FOR TOMS-UV OZONE DATA
!--------------------------
!        real(r8) :: gc_toms_input(pcols, 5) !toms,toms1,toms2,dtoms1,dtoms2
!--------------------------
!VARIABLES FOR UVALBEDO DATA
!--------------------------
!        real(r8) :: gc_uvalbedo_input(pcols)
        real(r8) :: aa1,aa2
        REAL(KIND( 0.0_4 )),        POINTER :: Ptr2D(:,:)
        REAL(KIND( 0.0_4 )),        POINTER :: Ptr3D(:,:,:)
        integer    :: ERR, FLAG
!--------------------------
!VARIABLES FOR BROMINE EMISSION
!--------------------------
!        real(r8) :: gc_chbr3_input(pcols)
!        real(r8) :: gc_ch2br2_input(pcols)

    INTEGER                 :: ND
    REAL(fp), TARGET        :: dflx(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Chm%nAdvect                       )
    TYPE(Species),  POINTER :: ThisSpc
    CHARACTER(LEN=63)      :: OrigUnit
    real      :: airdens(pcols,pver)
    real      :: bxhght(pcols,pver)
    real      :: ad(pcols,pver)
    real      :: eflx(pcols,pver,pcnst)
!--------------------------
!VARIABLES FOR BROMINE EMISSION
!--------------------------
!        real(r8) :: gc_chbr3_input(pcols)
!        real(r8) :: gc_ch2br2_input(pcols)
    real(r8) :: h2o2jv(pcols,plev) !sunrz 2024/3
    real(r8) :: hno3jv(pcols,plev)
    real(r8) :: no2jv(pcols,plev)
    real(r8) :: o3dv(pcols)
    real(r8) :: o3dflx(pcols)
    real(r8) :: no2dv(pcols)
    real(r8) :: no2dflx(pcols)
    real(r8) :: nodv(pcols)
    real(r8) :: nodflx(pcols)
    real(r8) :: hno3dv(pcols)
    real(r8) :: hcldv(pcols)
    real(r8) :: hno3dflx(pcols)
    real(r8) :: so4dv(pcols)
    real(r8) :: nh4dv(pcols)
    real(r8) :: hno3wd(pcols,plev)
    real(r8) :: no2wd(pcols,plev)
    real(r8) :: nowd(pcols,plev)
    real(r8) :: hclwd(pcols,plev)
    real(r8) :: hclwd2(pcols,plev)
    real(r8) :: gcpox(pcols,plev)
    real(r8) :: gclox(pcols,plev)
    real(r8) :: megan01(pcols,1)
    real(r8) :: megan02(pcols,1)
    real(r8) :: megan03(pcols,1)
    real(r8) :: megan04(pcols,1)
    real(r8) :: soil01(pcols,1)
    real(r8) :: soil02(pcols,1)
    real(r8) :: sunrz01(pcols,plev)
    real(r8) :: sunrz02(pcols,plev)
    real(r8) :: sunrz03(pcols,plev)
    real(r8) :: sunrz04(pcols,plev)
    real(r8) :: aod(pcols,plev)
#  endif

!-----------------------------------------------------------------------

  !zero = 0.
  lchnk = state%lchnk
  ncol  = state%ncol


!     write(*,*)'xiaolu',shape(state_chm%Species)
!     write(*,*)'xiaolu check state_chm #68-1',maxval(state_chm%Species(1,:,:,68))
!     write(*,*)'xiaolu check state_chm #69-1',maxval(state_chm%Species(1,:,:,81))

#if (defined MOZART2)
        !--------------------------
        !return LON/LAT from BCC
        !--------------------------
   call get_lon_all_p (lchnk, ncol, lon )
   call get_lat_all_p (lchnk, ncol, lat )
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

!write(6,*) 'clat=',clat,'clon=',clon
!write(6,*) 'lat=',lat,'lon=',lon,shape(lat),shape(lon)
!stop
        !------------------------
        !Set GEOS-Chem GRID DEMENSION
        !------------------------
!        IM=1
!        JM=ncol
!        LM=plev

! sunrz 2023/11 use local variables because IIPAR is not used in V14.1.1
!siipar=1
!sjjpar=pcols
!sllpar=plev

        allocate (lonCtr(ncol,1))
        allocate (latCtr(ncol,1))
!        lonCtr(1,:)=clon!*3.1415926535897932384626433/180!sunrz 2024/2
!        latCtr(1,:)=clat!*3.1415926535897932384626433/180!sunrz 2024/2
!        do i = 1,ncol
!        if (clon(i)>180/rad2deg) then
!        clon(i)=clon(i)-360/rad2deg
!        endif
!        enddo  
        lonCtr(:,1)=clon(:ncol)!*3.1415926535897932384626433/180!sunrz 2024/2
 !       if (masterproc) then 
 !       write(6,*) 'sunrz check lon=',(clon-180/rad2deg)*rad2deg, 'number=',lchnk,clat
        
 !       endif
 !       lonCtr(:,1)=clon
        latCtr(:,1)=clat(:ncol)!*3.1415926535897932384626433/180!sunrz 2024/2
!if (masterproc) then 
!write(6,*) 'sunrz check latedge=',latEdge(1,1)*rad2deg,'l=',lchnk
!write(6,*) 'sunrz check lonedge=',lonEdge(1,1)*rad2deg,'l=',lchnk
!write(6,*) 'sunrz check lat=',latCtr*rad2deg,'latedge=','l=',lchnk
!write(6,*) 'sunrz check lon=',lonCtr*rad2deg,'lonedge=','l=',lchnk
!endif
!        allocate (lonCtr(2,pcols/2))
!        allocate (latCtr(2,pcols/2))
!    ENDDO

    !-----------------------
    !create GEOS-Chem GRID
    !-----------------------
!    CALL SetGridFromCtr(masterproc, 1, pcols, lonCtr, latCtr, RC )
    CALL SetGridFromCtr(Input_Opt, State_Grid, lonCtr, latCtr, RC )!sunrz comment 2024/2
!   CALL SetGridFromCtrEdges( Input_Opt, State_Grid, lonCtr, latCtr, lonEdge,latEdge, RC )!sunrz 2024/2
!    ENDIF

        !----------------------
        !GET TIME INFORMATION
        !---------------------
        call get_curr_date(yr, mon, day, tod)
        hour = tod/3600
        minute = (tod-hour*3600) / 60
        second = (tod-hour*3600-minute*60) !/ 60!sunrz 2024/3

        nymd=yr*10000+mon*100+day                  !20010101
        nhms=hour*10000+minute*100+second          ! 4623
 
        utc=hour+minute/60+second/3600             !2
        hElapsed = get_nstep()*dtime/3600          !2
        dayOfYr= calday  !have problem if not start 1.1 and over 1 yr sunrz 2024/3
!write(6,*) 'sunrz check calday=',calday
!if (masterproc) then 
!write(6,*) 'sunrz check time',yr,mon,day,tod,nymd,nhms,dayOfYr,hour,minute,second,utc,hElapsed
!endif
        Phase=-1
        !Phase=2 ! hhan, 20210925
        IsChemTime=.True.
        IsRadTime = .False.!sunrz
        !WRITE(6,*) 'sunrz check CO in tphysbc_geoschem','max=',maxval(state%q(:,pver:1:-1,8))
!-----------------------------------
! Get surface emission (em3data and em2data)
! update,2018/03,xiaolu
!-----------------------------------
!   em3data(:,:,:) = 0.0_r8
!   em2data(:,:)   = 0.0
!   em2data_anthro(:,:) = 0.0
!   em3_IPCC6(:,:,:) = 0.0_r8
!   call gc_hook(ncol, prod_no, cnt, cnb, state%zm, state%zi, state%t, &
!                landfrac, ocnfrac, icefrac, pcols, clat, lchnk )
   call moz_hook(ncol, prod_no, cnt, cnb, state%zm, state%zi, state%t, &
                 landfrac, ocnfrac, icefrac, pcols, clat, lchnk )

   !
   ! out: em2data
   !
!!   call  getem2_IPCC( lchnk, ncol, em2data )
!!   call  getem2_IPCC_anthro( lchnk, ncol, em2data_anthro )
  
   !xiaolu,2019/01 ,read in TOMSUV here ; and UV-albedo
!   call get_tomsuv( lchnk, ncol, gc_toms_input)
!   State_Met%TO3(:,1) = gc_toms_input(:,1)*0+300

!   call get_uvalbedo( lchnk, ncol, gc_uvalbedo_input)
!   State_Met%UVALBEDO(:,1) = gc_uvalbedo_input(:)
   !xiaolu,2019/02, read in Bromine emission
   !call get_bromine_emis(lchnk, ncol, gc_chbr3_input,gc_ch2br2_input)

   ! write(*,*)'xiaolu check tomsuv',minval(gc_toms_input),maxval(gc_toms_input)

   !
   ! out: em3data
       !State_Chm%Species(1,:,:,i) = state%q(:,pver:1:-1,i+4)
       !ENDDO
     DO II=1,220
       State_Chm%Species(II)%Conc(:,1,:) = state%q(:ncol,pver:1:-1,II+3)
     ENDDO
!write(*,*) 'sunrz check before so4=',maxval(State_Chm%Species(213)%Conc),minval(State_Chm%Species(213)%Conc)
!                         second,    utc,       hElapsed,  Input_Opt,  &
!                         State_Chm, State_Met, Phase,     IsChemTime, &
!                         RC                                            )
    State_Grid%MaxTropLev  = State_Grid%NZ ! # trop. levels below
    State_Grid%MaxStratLev = State_Grid%NZ ! # strat. levels below

     call  GIGC_Chunk_Run   (   masterproc,                                  &
                             nymd,       nhms,       yr,           mon,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime, IsRadTime,              &
                             RC )


       dflx    =  0.0
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=====================================================================
       ! Apply dry deposition frequencies
       ! These are the frequencies calculated in drydep_mod.F90
       ! The HEMCO drydep frequencies (from air-sea exchange and
       ! PARANOX) were already added above.
       !
       ! NOTES:
       ! (1) Loops over only the drydep species
       ! (2) If drydep is turned off, nDryDep=0 and the loop won't execute
       ! (3) Tagged species are included in this loop. via species database
       !=====================================================================
       DO ND = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          IF ( N <= 0 ) CYCLE

          ! Point to the corresponding Species Database entry
          ThisSpc => State_Chm%SpcData(N)%Info

          ! only use the lowest model layer for calculating drydep fluxes
          ! given that spc is in v/v
          dflx(I,J,N) = dflx(I,J,N) + State_Chm%DryDepFreq(I,J,ND) &
                        * State_Chm%Species(N)%Conc(I,J,1) !      &
    !                             /  ( 28.9644 / ThisSpc%MW_g )
       !   stop
          ! Free species database pointer
          ThisSpc => NULL()
       ENDDO
       ! Convert DFLX from 1/s to kg/m2/s
       dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1)                        &
                                 / State_Grid%Area_M2(I,J)                    
    ENDDO
    ENDDO

     do i=1,pcols
     emilnox(i,:)=prod_no(i,:)*state_met%bxheight(i,1,pver:1:-1)*30*1e-3*1e6/6.02e23
     enddo
     call outfld('EMILNOX ',emilnox, pcols, lchnk)
      !DO i=1,2
      !WRITE(6,*) 'sunrz after chem check','max=',maxval(state%q(:,:,i+4)),'min=',minval(state%q(:,:,i+4)),'number:',i,shape(state%q(:,:,i+4))
      !ENDDO
!   State%q(:,:,5:220+4)=State_Chm%Species(1,:,pver:1:-1,1:220)
!     DO II=1,220
!        state%q(:ncol,pver:1:-1,II+4) = State_Chm%Species(II)%Conc(:,1,:)
!     ENDDO

!     DO II=1,220
!        state%q(:ncol,pver,II+4)=State_Chm%Species(II)%Conc(:,1,1)-dflx(:,1,II)/100*9.80665/State_Met%DELP_DRY(:,1,1)*dtime
!     ENDDO
 !    DO II=1,220
 !       state%q(:,pver,II+4) =State_Chm%Species(II)%Conc(:,1,1)-dflx(:,1,II)/100*9.80665/State_Met%DELP_DRY(:,1,1)*dtime
        !state%q(:,pver:1:-1,II+4)=state%q(:,pver:1:-1,II+4)+eflx(:,:,II+4)/100*9.80665/State_Met%DELP_DRY(:,1,:)*dtime

 !    ENDDO
    eflx(:,:,178)= eflx(:,:,178)+emilnox(:,pver:1:-1)
    CALL Convert_Spc_Units( Input_Opt,            State_Chm,   State_Grid,   &
                            State_Met,           'kg/m2',     RC,            &
                            OrigUnit   )


     DO II=1,220
        State_Chm%Species(II)%Conc(:,1,1)=State_Chm%Species(II)%Conc(:,1,1)-dflx(:,1,II)*dtime
        State_Chm%Species(II)%Conc(:,1,:)=State_Chm%Species(II)%Conc(:,1,:)+eflx(:ncol,:,II+3)*dtime
     ENDDO

    CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                            OrigUnit,  RC )

     DO II=1,220
        state%q(:ncol,pver:1:-1,II+3) =State_Chm%Species(II)%Conc(:,1,:)
     ENDDO

     DEALLOCATE(latctr)
   !  DEALLOCATE(latedge)
     DEALLOCATE(lonctr)
   !  DEALLOCATE(lonedge)

    do ii=1,ncol
        airdens(ii,:)=State_Met%AIRDEN(ii,1,pver:1:-1)
        bxhght(ii,:)=State_Met%BXHEIGHT(ii,1,pver:1:-1)
        ad(ii,:)=State_Met%AD(ii,1,pver:1:-1)
    enddo

    call outfld('AIRDENS',airdens,pcols,lchnk)
    call outfld('BXHGHT', bxhght,pcols,lchnk)
    call outfld('AD', ad,pcols,lchnk)

     sunrz01=State_Diag%OHconcAfterChem(:,1,pver:1:-1)
     call outfld('BCCOH',sunrz01,pcols,lchnk)
     o3dv=State_Diag%DryDepVel(:,1,128)
     o3dflx=State_Diag%DryDep(:,1,128)
     no2dv=State_Diag%DryDepVel(:,1,125)
     no2dflx=State_Diag%DryDep(:,1,125)
     hno3dv=State_Diag%DryDepVel(:,1,49)
     hcldv=State_Diag%DryDepVel(:,1,43)
     hno3dflx=State_Diag%DryDep(:,1,49)
     so4dv=State_Diag%DryDepVel(:,1,155)
     nh4dv=State_Diag%DryDepVel(:,1,122)
     call outfld('O3DV',o3dv,pcols,lchnk)
     call outfld('O3DFLX',dflx(:,1,180),pcols,lchnk)
     call outfld('NO2DV',no2dv,pcols,lchnk)
     call outfld('NO2DFLX',no2dflx,pcols,lchnk)
     call outfld('NO2DFLX2',dflx(:,1,176),pcols,lchnk)
     call outfld('HClDV',hno3dv,pcols,lchnk)
     call outfld('HClDFLX',dflx(:,1,83),pcols,lchnk)
!     call outfld('HNO3DV',hcldv,pcols,lchnk)
!     call outfld('HNO3DFLX',hno3dflx,pcols,lchnk)
!     call outfld('NH4DV',nh4dv,pcols,lchnk)
!     call outfld('SO4DV',so4dv,pcols,lchnk)
     hno3wd=State_Diag%WetLossLs(:,1,pver:1:-1,41)
!     call outfld('HNO3WD',hno3wd,pcols,lchnk)
     hclwd=State_Diag%WetLossLs(:,1,pver:1:-1,35)
     call outfld('HClWD',hclwd,pcols,lchnk)
     hclwd2=State_Diag%WetLossConv(:,1,pver:1:-1,35)
     call outfld('HClWD2',hclwd2,pcols,lchnk)
     gcpox=State_Diag%Prod(:,1,pver:1:-1,1)
     gclox=State_Diag%Loss(:,1,pver:1:-1,1)
     call outfld('GCPOx',gcpox,pcols,lchnk)
     call outfld('GCLOx',gclox,pcols,lchnk)
     aod=State_Diag%AODHygWL1(:,1,pver:1:-1,1)+State_Diag%AODHygWL1(:,1,pver:1:-1,2)+State_Diag%AODHygWL1(:,1,pver:1:-1,3)+State_Diag%AODHygWL1(:,1,pver:1:-1,4)+State_Diag%AODHygWL1(:,1,pver:1:-1,5)
     call outfld('ODAER',aod,pcols,lchnk)
#endif

   return
end subroutine tphysbc_geoschem
#endif
