!xiaolu
!---------------------------
!Purpose:
!ROUTINE FOR GEOS-CHEM DRY DEPOSTION
!FROM MO_SURFACE IN MOZART2

!REVISION HISTORY:
!
!XIAO LU,2017/04/27:ONLY MAINTAIN CALL GAS DRY DEPOSITION PART
!XIAO LU,2017/06/01
!---------------------------

      module gchp_surface

      implicit none

      public    ::    inisflx, gc_sflxdr 
      save

      real :: rair, gravit

      contains

!--------------------------
!SUBROUTINE inisflx:
!set the physical constants for surface flux routines
!obtain from MOZART
!XIAOLU,2017/06/01
!---------------------------
      subroutine inisflx( xrair, xgravit )

      implicit none

!----------------------------------------------------------------------------
! 	... dummy arguments:
!----------------------------------------------------------------------------

      real, intent(in) :: &
        xrair, &    ! gas constant for dry air
        xgravit     ! gravitational acceleration

!----------------------------------------------------------------------------
!  	... set the physical constants for surface flux routines
!----------------------------------------------------------------------------
      rair   = xrair
      gravit = xgravit

      end subroutine inisflx
!===========================
!===========================


!==========================
!SUBROUTINE gc_sflxdr
!obtain from sflxdr from MOZART
!XIAOLU,2017/06/01
!==========================
      subroutine gc_sflxdr( delt, calday  ,as  ,pmid  ,pdel    ,ts    ,   &
                         q       ,t   ,plonl ,  cflx,  &
                         lchnk   ,ncol,landfrac , &
                         flux_ISOP      , flux_ACET      ,  flux_C3H6,      &
                         flux_C2H4      , flux_OC2       ,  flux_C10H16,    &
                         flux_CO2       , flux_N2O       ,  flux_DMS,   &
                         nem2_anthro    , em2data_anthro ,                  &
                          nem2    ,em2data, &
                         gc_chbr3_input,gc_ch2br2_input)

!-----------------------------------------------------------------------
! 	... set surface fluxes or make adjustment to field due to
!           surface flux. dry deposition routines are called from
!           here because they work by setting surface fluxes.
!-----------------------------------------------------------------------

      use gc_grid,         only : plev, pcnstm1, plon !cam1/
        !check pcnstm1 in gc_grid
      use gchp_const_mozart, only : pi, rearth !geoschem/BCCGCHP/
      use gchp_sim_chm,      only : latwtsbdy  
!      use gchp_srf_emis,     only : srf_emis   
      use gchp_drydep,       only : gc_drydep
      use gchp_photo,        only : diurnal_geom !geoschem/BCCGCHP/
      use phys_grid,       only : get_rlat_all_p, get_rlon_all_p       !2008.08.11
      use history,         only : outfld
      use gchp_chem_utls,    only : get_spc_ndx

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl
      integer, intent(in) :: lchnk
      integer, intent(in) :: ncol 

      real, intent(in) :: delt
      real, intent(in) :: calday                    ! julian day + fraction (greenwich)
      real, intent(in) :: pmid(plonl,plev)          ! pressure at layer midpoints
      real, intent(in) :: pdel(plonl,plev)          ! 
      real, intent(in) :: ts(plonl)                 ! surface temperature
      real, intent(in) :: t(plonl,plev)             ! temperature
      real, intent(in) :: q(plonl,plev)             ! specific humidity 
      real, intent(in) :: landfrac(plonl)                ! landfrac 

      real,intent(in)    :: flux_ISOP(plonl)
      real,intent(in)    :: flux_ACET(plonl)
      real,intent(in)    :: flux_C3H6(plonl)
      real,intent(in)    :: flux_C2H4(plonl)
      real,intent(in)    :: flux_OC2(plonl)
      real,intent(in)    :: flux_C10H16(plonl)

      real,intent(in)    :: flux_CO2(plonl)
      real,intent(in)    :: flux_N2O(plonl)
      real,intent(in)    :: flux_DMS(plonl)

       !xiaolu,add for emission 2017/06/13
      integer, intent(in) :: nem2
      real, intent(in)    :: em2data(plonl,nem2)
       !update ,2018/03 
      integer, intent(in) :: nem2_anthro
      real, intent(in)    :: em2data_anthro(plonl,nem2_anthro)

       !xiaolu,add bromine emission
      real, intent(in)    :: gc_chbr3_input(plonl)
      real, intent(in)    :: gc_ch2br2_input(plonl)
!----
      real, intent(in)    :: as(plonl,plev,pcnstm1)   ! advected species
      real, intent(inout) :: cflx(plonl,pcnstm1)        ! surface flux for advected species

      integer ::  ioro(plonl)                       ! orography 
!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: dms_ndx, so2_ndx, cb1_ndx, cb2_ndx, oc1_ndx, oc2_ndx
      integer :: ch3coch3_ndx, c2h4_ndx, c2h5oh_ndx, c2h6_ndx, c3h6_ndx
      integer :: c3h8_ndx, ch2o_ndx, ch3cho_ndx, ch3oh_ndx, n2o_ndx
      integer :: co_ndx, h2_ndx, isop_ndx, no_ndx, c4h10_ndx, ch4_ndx
      integer :: nh3_ndx, c10h16_ndx
      integer :: c2h2_ndx,c6h6_ndx,c7h8_ndx,c8h10_ndx
      integer :: chbr3_ndx,ch2br2_ndx
      integer :: i, m, n, base, file
      real    :: sunon(plonl)         ! sunrise longitude (radians)
      real    :: sunoff(plonl)        ! sunset longitude (radians)
      real    :: zen_angle(plonl)     ! zenith angle as function of longitude (radians)
      real    :: loc_angle(plonl)     ! time angle as function of longitude (radians)
      real    :: depvel(plonl,pcnstm1)  ! deposition velocity ( cm/s )

      real    :: dflx(plonl,pcnstm1)    ! deposition flux ( kg/m^2/s )
      real    :: sflx(plonl,pcnstm1) 

      real    :: dep_flx(plonl)       ! deposition flux ( kg/m^2/s )
      real    :: clon(plonl),clat(plonl)
      real    :: tv(plonl)            ! virtual temperature in surface layer (K)
      real    :: hsa_fac
      real    :: sflx2vmr(plonl)

      real    :: sflx_c5h12(plonl)
      real    :: sflx_c6h14(plonl)

      logical :: polar_night(plonl)          ! wrk flag for diurnal_geom
      logical :: polar_day(plonl)            ! wrk flag for diurnal_geom
   
      LOGICAL, SAVE                  :: FIRSTDRYDEP = .TRUE.

      character(len=32) :: fldname
      character(len=3)  :: num

      do i = 1, plonl
        if(landfrac(i) .gt. 0.5) then
          ioro(i) = 1
        else
          ioro(i) = 0
        endif
      enddo


      do m = 1,pcnstm1
        sflx(:,m) = 0.
        dflx(:,m) = 0.
      end do

!----------------------------------------------------------------------
!virtual temperature
!----------------------------------------------------------------------
      call virtem(ncol, plonl, 1, t(1,plev), q(1,plev), 0.61, tv)  !2008.08.20
!-----------------------------------------------------------------------
!	... diurnal geometry
!-----------------------------------------------------------------------
      call get_rlat_all_p(lchnk, ncol, clat)    !add by zf 2008.08.11
      call get_rlon_all_p(lchnk, ncol, clon)    !add by zf 2008.08.11

      do i = 1, plonl
        call diurnal_geom( clon(i), clat(i), calday, polar_night(i), polar_day(i), &
                           sunon(i), sunoff(i), loc_angle(i), zen_angle(i) )
      enddo
!-----------------------------------------------------------------------
!	... the surface emissions
! XIAOLU,2017/06/06:do emission later 
! XIAOLU,2017/06/13:TURN ON emission module,add nem2, em2data
!-----------------------------------------------------------------------

!      call srf_emis( lchnk, clat, calday, sflx, ioro, &
!                     loc_angle, polar_night, polar_day, sunon, sunoff, plonl , &
!                     flux_ISOP      , flux_ACET      ,  flux_C3H6,             &
!                     flux_C2H4      , flux_OC2       ,  flux_C10H16,           &
!                     flux_CO2       , flux_N2O       ,  flux_DMS,       &
!                     nem2_anthro    , em2data_anthro ,                  &
!                     nem2, em2data,sflx_c5h12, sflx_c6h14)
!-------------------2015-05-29------------------------------------------
      dms_ndx       = get_spc_ndx( 'DMS' )
      so2_ndx       = get_spc_ndx( 'SO2' )
      !xiaolu, change CB,OC
      !cb1_ndx       = get_spc_ndx( 'CB1' )
      !cb2_ndx       = get_spc_ndx( 'CB2' )
      !oc1_ndx       = get_spc_ndx( 'OC1' )
      !oc2_ndx       = get_spc_ndx( 'OC2' )

      cb1_ndx       = get_spc_ndx( 'BCPO' )
      cb2_ndx       = get_spc_ndx( 'BCPI' )
      oc1_ndx       = get_spc_ndx( 'OCPO' )
      oc2_ndx       = get_spc_ndx( 'OCPI' )


      !ch3coch3_ndx  = get_spc_ndx( 'CH3COCH3' )
      ch3coch3_ndx  = get_spc_ndx( 'ACET')
      c2h4_ndx      = get_spc_ndx( 'C2H4' )
      !c2h5oh_ndx    = get_spc_ndx( 'C2H5OH' )
      c2h5oh_ndx    = get_spc_ndx( 'EOH' )
      c2h6_ndx      = get_spc_ndx( 'C2H6' )
      !c3h6_ndx      = get_spc_ndx( 'C3H6' )
      c3h6_ndx      = get_spc_ndx( 'PRPE' )
      c3h8_ndx      = get_spc_ndx( 'C3H8' )
      ch2o_ndx      = get_spc_ndx( 'CH2O' )
      !ch3cho_ndx    = get_spc_ndx( 'CH3CHO' )
      ch3cho_ndx    = get_spc_ndx( 'ALD2' )
      ch3oh_ndx     = get_spc_ndx( 'CH3OH' )
      n2o_ndx       = get_spc_ndx( 'N2O' )
      co_ndx        = get_spc_ndx( 'CO' )
      h2_ndx        = get_spc_ndx( 'H2' )
      isop_ndx      = get_spc_ndx( 'ISOP' )
      no_ndx        = get_spc_ndx( 'NO' )
      !c4h10_ndx     = get_spc_ndx( 'C4H10' )
      c4h10_ndx       = get_spc_ndx( 'ALK4' )
      ch4_ndx       = get_spc_ndx( 'CH4' )
      nh3_ndx       = get_spc_ndx( 'NH3' )
      c10h16_ndx    = get_spc_ndx( 'C10H16' )

      !-------------------
      !xiaolu,2019/02
      !add bromine emission
      !-------------------
      chbr3_ndx     = get_spc_ndx( 'CHBr3')
      ch2br2_ndx    = get_spc_ndx( 'CH2Br2')

      IF (chbr3_ndx > 1) THEN
        sflx(:,chbr3_ndx)=gc_chbr3_input
        call outfld('EM_CHBr3',sflx(:,chbr3_ndx), plonl, lchnk)
      ENDIF

      IF (ch2br2_ndx > 1) THEN
        sflx(:,ch2br2_ndx)=gc_ch2br2_input
        call outfld('E_CH2Br2',sflx(:,ch2br2_ndx), plonl, lchnk)
      ENDIF

      !---------------------
      !xiaolu,2019/04
      !add VOCs
      !--------------------
      c2h2_ndx       = get_spc_ndx('C2H2')
      c6h6_ndx       = get_spc_ndx('BENZ')
      c7h8_ndx       = get_spc_ndx('TOLU')
      c8h10_ndx      = get_spc_ndx('XYLE')


      !write(*,*)'xiaolu check index emis',c2h2_ndx,c6h6_ndx,c7h8_ndx,c8h10_ndx
      IF (c6h6_ndx > 1 .and. c7h8_ndx>1 .and. c8h10_ndx > 1) THEN 
      !call outfld('EMS_C2H2',sflx(:,c2h2_ndx),plonl, lchnk)
      call outfld('EMS_BENZ',sflx(:,c6h6_ndx),plonl, lchnk)
      call outfld('EMS_TOLU',sflx(:,c7h8_ndx),plonl, lchnk)
      call outfld('EMS_XYLE',sflx(:,c8h10_ndx),plonl, lchnk)
      call outfld('EM_C5H12',sflx_c5h12,plonl, lchnk)
      call outfld('EM_C6H14',sflx_c6h14,plonl, lchnk)
      ENDIF

      call outfld('EMIS_DMS',sflx(:,dms_ndx), plonl, lchnk)
      call outfld('EMIS_SO2',sflx(:,so2_ndx), plonl, lchnk)
      call outfld('EMIS_CB1',sflx(:,cb1_ndx), plonl, lchnk)
      call outfld('EMIS_CB2',sflx(:,cb2_ndx), plonl, lchnk)
      call outfld('EMIS_OC1',sflx(:,oc1_ndx), plonl, lchnk)
      call outfld('EMIS_OC2',sflx(:,oc2_ndx), plonl, lchnk)

      call outfld('EMS_ACET',sflx(:,ch3coch3_ndx), plonl, lchnk)
      call outfld('EMS_C2H4',sflx(:,c2h4_ndx), plonl, lchnk)
      call outfld('EMS_EOH', sflx(:,c2h5oh_ndx), plonl, lchnk)
      call outfld('EMS_C2H6',sflx(:,c2h6_ndx), plonl, lchnk)
      call outfld('EMS_C3H6',sflx(:,c3h6_ndx), plonl, lchnk)
      call outfld('EMS_C3H8',sflx(:,c3h8_ndx), plonl, lchnk)
      call outfld('EMS_CH2O',sflx(:,ch2o_ndx), plonl, lchnk)
      call outfld('EMS_ALD2',sflx(:,ch3cho_ndx), plonl, lchnk)
      call outfld('EMS_CH4O',sflx(:,ch3oh_ndx), plonl, lchnk)
      call outfld('EMS_N2O ',sflx(:,n2o_ndx), plonl, lchnk)
      call outfld('EMS_CO  ',sflx(:,co_ndx), plonl, lchnk)
     ! WRITE(6,*) 'sunrz check sflx(:,co_ndx):', maxval(sflx(:,co_ndx))
      call outfld('EMS_H2  ',sflx(:,h2_ndx), plonl, lchnk)
      call outfld('EMS_ISOP',sflx(:,isop_ndx), plonl, lchnk)
      call outfld('EMS_NO  ',sflx(:,no_ndx), plonl, lchnk)
      call outfld('EMS_ALK4',sflx(:,c4h10_ndx), plonl, lchnk)
      call outfld('EMS_CH4 ',sflx(:,ch4_ndx), plonl, lchnk)
      call outfld('EMS_NH3 ',sflx(:,nh3_ndx), plonl, lchnk)
      call outfld('E_C10H16',sflx(:,c10h16_ndx), plonl, lchnk)

      call outfld('ISOP_flx',flux_ISOP, plonl, lchnk)
      call outfld('ACET_flx',flux_ACET, plonl, lchnk)
      call outfld('C3H6_flx',flux_C3H6, plonl, lchnk)
      call outfld('C2H4_flx',flux_C2H4, plonl, lchnk)
      call outfld('OC2_flux',flux_OC2 , plonl, lchnk)
      call outfld('C10H16_f',flux_C10H16, plonl, lchnk)

      call outfld('CO2_flux',flux_CO2 , plonl, lchnk)
      call outfld('N2O_flux',flux_N2O , plonl, lchnk)
      call outfld('DMS_flux',flux_DMS , plonl, lchnk)
!-----------------------------------------------------------------------
!XIAOLU,2017/06
!MAIN ROUTINE TO CALCULATE DRY DEPOSITION
!	... the dry deposition
!-----------------------------------------------------------------------
!      call gc_drydep( lchnk, calday, ts, zen_angle, &
!                   depvel, dflx, rair, as, pmid(1,plev), tv, plonl )

!-----------------------------------------------------------------------
!       ... dry deposition velocity to history files
!-----------------------------------------------------------------------
!      call outfld('OXDV    ',depvel(:,1)  ,plonl,lchnk)
!      call outfld('NO2DV   ',depvel(:,5)  ,plonl,lchnk)
!      call outfld('HNO3DV  ',depvel(:,7)  ,plonl,lchnk)
!      call outfld('CH4DV   ',depvel(:,10)  ,plonl,lchnk)
!      call outfld('CH3OOHDV',depvel(:,12)  ,plonl,lchnk)
!      call outfld('CH2ODV  ',depvel(:,13)  ,plonl,lchnk)
       call outfld('CODV    ',depvel(:,co_ndx)  ,plonl,lchnk)
!      call outfld('H2O2DV  ',depvel(:,17)  ,plonl,lchnk)
!      call outfld('POOHDV  ',depvel(:,22)  ,plonl,lchnk)
!      call outfld('CH3CODV ',depvel(:,24)  ,plonl,lchnk)
!      call outfld('PANDV   ',depvel(:,25)  ,plonl,lchnk)
!      call outfld('MPANDV  ',depvel(:,30)  ,plonl,lchnk)
!      call outfld('C2H5ODV ',depvel(:,38)  ,plonl,lchnk)
!      call outfld('ONITDV  ',depvel(:,26)  ,plonl,lchnk)
!      call outfld('C3H7HDV ',depvel(:,42)  ,plonl,lchnk)
!      call outfld('ROOHDV  ',depvel(:,44)  ,plonl,lchnk)
!      call outfld('CH3CnDV ',depvel(:,53)  ,plonl,lchnk)
!      call outfld('CH3C4DV ',depvel(:,43)  ,plonl,lchnk)
!      call outfld('PbDV    ',depvel(:,55)  ,plonl,lchnk)
!      call outfld('O3INERDV',depvel(:,63)  ,plonl,lchnk)
!      call outfld('O3SDV   ',depvel(:,62)  ,plonl,lchnk)
!      call outfld('H2DV    ',depvel(:,61)  ,plonl,lchnk)
!      call outfld('ONITRDV ',depvel(:,57)  ,plonl,lchnk)
!      call outfld('MACRHDV ',depvel(:,35)  ,plonl,lchnk)
!      call outfld('XOOHDV  ',depvel(:,59)  ,plonl,lchnk)
!      call outfld('ISOPHDV ',depvel(:,60)  ,plonl,lchnk)
!      call outfld('CH3CHODV',depvel(:,21)  ,plonl,lchnk)
!      call outfld('NODV    ',depvel(:,4)  ,plonl,lchnk)
!      call outfld('HO2NO2DV',depvel(:,8)  ,plonl,lchnk)
!      call outfld('GLYALDDV',depvel(:,47)  ,plonl,lchnk)
!      call outfld('HYACDV  ',depvel(:,48)  ,plonl,lchnk)
!      call outfld('CH3OHDV ',depvel(:,45)  ,plonl,lchnk)
!      call outfld('C2H5OHDV',depvel(:,46)  ,plonl,lchnk)
!      call outfld('HYDRDDV ',depvel(:,51)  ,plonl,lchnk)

!-----------------------------------------------------------------------
!       ... dry deposition flux to history files
!-----------------------------------------------------------------------
      hsa_fac = 2*pi*rearth*rearth/real(plon)

      !call sim_chm_data
      !  write(*,*),'xiaolu check latwtsbdy',shape(latwtsbdy)
      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,1)
!      call outfld('OXDFLX  ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,62)
!      call outfld('O3SDFLX ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,co_ndx)
      call outfld('CODFLX  ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,4)
!      call outfld('NODFLX  ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,5)
!      call outfld('NO2DFL  ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,7)
!      call outfld('HNO3DFLX',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,25)
!      call outfld('PANDFLX ',dep_flx  ,plonl,lchnk)

      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,30)
!      call outfld('MPANDFLX',dep_flx  ,plonl,lchnk)

!      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,10)
!      call outfld('CH4DFLX ',dep_flx  ,plonl,lchnk)

!      dep_flx(:) = hsa_fac*latwtsbdy(:,lchnk)*delt*dflx(:,61)
!      call outfld('H2DFLX  ',dep_flx  ,plonl,lchnk)

!----------------------------------------------------------------------
!	... form surface flux
!-----------------------------------------------------------------------
      do m = 1,pcnstm1
!comment dlfx as we do not use BCC dep
!        cflx(:,m) = cflx(:,m) + (sflx(:,m) - dflx(:,m))
         cflx(:,m) = cflx(:,m) + sflx(:,m)
      end do

!-------------------------------
!      do m = 1, pcnstm1
!         as(:,plev,m) = as(:,plev,m) + sflx(:ncol,m) /pdel(:ncol,plev)*gravit * delt
!      end do
!--------------------------------------
      end subroutine gc_sflxdr

      end module gchp_surface
