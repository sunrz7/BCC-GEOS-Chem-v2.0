!xiaolu
!---------------------------------------
!from MOZART2 mo_srf_emis.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/07
!---------------------------------------

#include <misc.h>
#include <params.h>

      module gchp_srf_emis
!---------------------------------------------------------------
! 	... Interpolation surface emissions
!---------------------------------------------------------------

      use pmgrid,         only: plon, plat, masterproc

      implicit none

      private
      public  :: srf_emis_inti, nppbdy, baseflux, landflux

      character(len=256), public  :: bndacet    ! full pathname for Mozart
      character(len=256), public  :: bndc2h4    ! full pathname for Mozart
      character(len=256), public  :: bndc3h6    ! full pathname for Mozart
      character(len=256), public  :: bndc2h6    ! full pathname for Mozart
      character(len=256), public  :: bndc3h8    ! full pathname for Mozart
      character(len=256), public  :: bndch2o    ! full pathname for Mozart
      character(len=256), public  :: bndch3oh   ! full pathname for Mozart
      character(len=256), public  :: bndch4     ! full pathname for Mozart
      character(len=256), public  :: bndco      ! full pathname for Mozart
      character(len=256), public  :: bndh2      ! full pathname for Mozart
      character(len=256), public  :: bndisoprene! full pathname for Mozart
      character(len=256), public  :: bndn2o     ! full pathname for Mozart
      character(len=256), public  :: bndnmv     ! full pathname for Mozart
      character(len=256), public  :: bndnox     ! full pathname for Mozart
      character(len=256), public  :: bndnpp     ! full pathname for Mozart
      character(len=256), public  :: bndterpense! full pathname for Mozart

      save

      real, allocatable, dimension(:,:,:) :: &
               noxembdy, ch4embdy, ch2oembdy, coembdy, c3h6embdy, &
               isoembdy, terpembdy, c2h4embdy, c2h6embdy, c4h10embdy, &
               n2oembdy, c3h8embdy, acetembdy, h2embdy, ch3ohembdy, &
               nppbdy

      real ::  baseflux = 0., &
               landflux = 0.

      real    :: sf(plat)
      real    :: factor
      integer :: no_ndx, ch4_ndx, co_ndx, c3h6_ndx, isop_ndx, c10h16_ndx, &
                 c2h4_ndx, c2h6_ndx, c3h8_ndx, ch3coch3_ndx, c4h10_ndx, &
		 n2o_ndx, h2_ndx, ch3oh_ndx, Rn_ndx, ch2o_ndx, &
                 c2h2_ndx,c6h6_ndx,c7h8_ndx,c8h10_ndx                

      integer :: dms_ndx, so2_ndx, cb1_ndx, cb2_ndx, oc1_ndx, oc2_ndx
      integer :: c2h5oh_ndx, ch3cho_ndx,  hcn_ndx, nh3_ndx ,ch3cocho_ndx 
 
      logical :: no_em, ch4_em, co_em, c3h6_em, isop_em, c10h16_em, &
                 c2h4_em, c2h6_em, c3h8_em, ch3coch3_em, c4h10_em, &
		 n2o_em, h2_em, ch3oh_em, Rn_em, ch2o_em

      contains

      subroutine srf_emis_inti( plonl, platl, pplon )
!-----------------------------------------------------------------------
! 	... Read the primary emission flux datasets 
!-----------------------------------------------------------------------

      use gchp_chem_utls,     only : get_spc_ndx, has_srfems
      use gchp_sim_chm,       only : latwts    !zf test
      use gc_chem_mods,     only : adv_mass
      use gchp_const_mozart,  only : d2r, pi, rearth

!add by zf 2008.03.13
      use ppgrid,         only:  pcols, begchunk, endchunk   !2008.08.11
      use phys_grid,      only:  scatter_field_to_chunk      !2008.08.11
      use ioFileMod,      only: getfil
      use abortutils,     only: endrun
!end zf
      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl
      integer, intent(in) :: platl
      integer, intent(in) :: pplon

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      real, allocatable, dimension(:,:,:) :: &
               noxem, ch4em, ch2oem, coem, c3h6em, &
               isoem, terpem, c2h4em, c2h6em, c4h10em, &
               n2oem, c3h8em, acetem, h2em, ch3ohem, &
               npp

      real, parameter :: amufac = 1.65979e-23         ! 1.e4* kg / amu
      integer :: astat
      integer :: i, j, ncid                              ! Indices
      real    :: mw
      real    :: seq
      character(len=10) :: species
      character(len=256) :: locfn    ! netcdf local filename to open   add by zf 2008.03.13

      write(*,*) "In gchp_srf_emis.F90, subroutine srf_emis_inti, xpye"
!-----------------------------------------------------------------------
! 	... Setup the indicies
!-----------------------------------------------------------------------
      no_ndx       = get_spc_ndx( 'NO' )
      ch4_ndx      = get_spc_ndx( 'CH4' )
      ch2o_ndx     = get_spc_ndx( 'CH2O' )
      co_ndx       = get_spc_ndx( 'CO' )
!xiaolu change c3h6 to PRPE
!      c3h6_ndx     = get_spc_ndx( 'C3H6' )
      c3h6_ndx     = get_spc_ndx( 'PRPE' )
      isop_ndx     = get_spc_ndx( 'ISOP' )
      c10h16_ndx   = get_spc_ndx( 'C10H16' )
      c2h4_ndx     = get_spc_ndx( 'C2H4' )
      c2h6_ndx     = get_spc_ndx( 'C2H6' )
      c3h8_ndx     = get_spc_ndx( 'C3H8' )
!xiaolu change CH3COCH3 to ACET
!      ch3coch3_ndx = get_spc_ndx( 'CH3COCH3' )
      ch3coch3_ndx = get_spc_ndx( 'ACET' )
!xiaolu change C4H10 to ALK4
      c4h10_ndx    = get_spc_ndx( 'ALK4' )
      write(*,*)'ALK4',c4h10_ndx
      n2o_ndx      = get_spc_ndx( 'N2O' )
      h2_ndx       = get_spc_ndx( 'H2' )
      ch3oh_ndx    = get_spc_ndx( 'CH3OH' )
      Rn_ndx       = get_spc_ndx( 'Rn' )

      dms_ndx       = get_spc_ndx( 'DMS' )
      so2_ndx       = get_spc_ndx( 'SO2' )

      !xiaolu change OC/BC to match CB/OC
!      cb1_ndx       = get_spc_ndx( 'CB1' )
!      cb2_ndx       = get_spc_ndx( 'CB2' )
!      oc1_ndx       = get_spc_ndx( 'OC1' )
!      oc2_ndx       = get_spc_ndx( 'OC2' )

      cb1_ndx       = get_spc_ndx( 'BCPO' )
      cb2_ndx       = get_spc_ndx( 'BCPI' )
      oc1_ndx       = get_spc_ndx( 'OCPI' )
      oc2_ndx       = get_spc_ndx( 'OCPO' )
!xiaolu,change C2H5OH to EOH
!     c2h5oh_ndx    = get_spc_ndx( 'C2H5OH' )
      c2h5oh_ndx    = get_spc_ndx('EOH')
!xiaolu,change ch3cho to ALD2
!      ch3cho_ndx    = get_spc_ndx( 'CH3CHO' )
      ch3cho_ndx    = get_spc_ndx( 'ALD2' )
      hcn_ndx       = get_spc_ndx( 'HCN' )
      nh3_ndx       = get_spc_ndx( 'NH3' )

!xiaolu,2019/04
      c2h2_ndx       = get_spc_ndx('C2H2')
      c6h6_ndx       = get_spc_ndx('BENZ')
      c7h8_ndx       = get_spc_ndx('TOLU')
      c8h10_ndx      = get_spc_ndx('XYLE')

!      write(*,*)'xiaolu check c6h6', c6h6_ndx
      no_em       = has_srfems( 'NO' )
      ch4_em      = has_srfems( 'CH4' )
      ch2o_em     = has_srfems( 'CH2O' )
      co_em       = has_srfems( 'CO' )
!     c3h6_em     = has_srfems( 'C3H6' )
      c3h6_em     = has_srfems( 'PRPE' )
      isop_em     = has_srfems( 'ISOP' )
      c10h16_em   = has_srfems( 'C10H16' )
      c2h4_em     = has_srfems( 'C2H4' )
      c2h6_em     = has_srfems( 'C2H6' )
      c3h8_em     = has_srfems( 'C3H8' )
!     ch3coch3_em = has_srfems( 'CH3COCH3' )
!     c4h10_em    = has_srfems( 'C4H10' )
!xiaolu change tracer name
      ch3coch3_em = has_srfems( 'ACET' )
      c4h10_em    = has_srfems( 'ALK4' )

      n2o_em      = has_srfems( 'N2O' )
      h2_em       = has_srfems( 'H2' )
      ch3oh_em    = has_srfems( 'CH3OH' )
      Rn_em       = has_srfems( 'Rn' )

      if (masterproc) then
      write(*,*) 'srf_emis_inti : diagnostic'
      write(*,'(10i5)') no_ndx, ch4_ndx, co_ndx, c3h6_ndx, isop_ndx, c10h16_ndx, &
                 c2h4_ndx, c2h6_ndx, c3h8_ndx, ch3coch3_ndx, c4h10_ndx, &
		 n2o_ndx, h2_ndx, ch3oh_ndx, Rn_ndx, ch2o_ndx
      write(*,*)  no_em, ch4_em, co_em, c3h6_em, isop_em, c10h16_em, &
                 c2h4_em, c2h6_em, c3h8_em, ch3coch3_em, c4h10_em, &
		 n2o_em, h2_em, ch3oh_em, Rn_em, ch2o_em
      endif

!-----------------------------------------------------------------------
! 	... allocate working arrays
!-----------------------------------------------------------------------
      allocate( noxem(plonl,platl,12), ch4em(plonl,platl,12), &
                ch2oem(plonl,platl,12), coem(plonl,platl,12), &
                c3h6em(plonl,platl,12), isoem(plonl,platl,12), &
                terpem(plonl,platl,12), c2h4em(plonl,platl,12), &
                c2h6em(plonl,platl,12), c4h10em(plonl,platl,12), &
                n2oem(plonl,platl,12), c3h8em(plonl,platl,12), &
                acetem(plonl,platl,12), h2em(plonl,platl,12), &
                ch3ohem(plonl,platl,12), npp(plonl,platl,12), &
		stat=astat )
      if( astat/= 0 ) then
	 write(*,*) 'srf_emis_inti: failed to allocate wrking arrays; error = ',astat
	 call endrun
      end if

!-----------------------------------------------------------------------
! 	... Diagnostics setup
!-----------------------------------------------------------------------
      seq = 2.*pi*1.e4*rearth**2/REAL(plon)
      sf(:plat) = seq*latwts(:plat)
      factor = 86400. * 365. &   ! sec / year
             / 6.022e23 &        ! molec / mole 
             * 1.e-12            ! Tg / g

      if (masterproc) then
!********************* NOT USE?***************************************
!-----------------------------------------------------------------------
! 	... NO
!-----------------------------------------------------------------------
      if( no_em ) then
         species  = 'NO'
         mw       = 14.00674                      ! g N / mole
         call getfil (bndnox, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgN/y', ncid, noxem )
         noxem(:,:,:)   = amufac * adv_mass(no_ndx) * noxem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... CO
!-----------------------------------------------------------------------
      if( co_em ) then
         species  = 'CO'
         mw       = adv_mass(co_ndx)               ! g / mole
         write(*,*)"CO mole mass",mw
         write(*,*)"all mass",adv_mass
         call getfil (bndco, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, coem )
         coem(:,:,:) = amufac * adv_mass(co_ndx) * coem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... CH4
!-----------------------------------------------------------------------
      if( ch4_em ) then
         species  = 'CH4'
         mw       = adv_mass(ch4_ndx)               ! g / mole
         call getfil (bndch4, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, ch4em )
         ch4em(:,:,:)   = amufac * adv_mass(ch4_ndx) * ch4em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... C2H6
!-----------------------------------------------------------------------
      if( c2h6_em ) then
         species  = 'C2H6'
         mw       = 2.*12.011                       ! g C / mole
         call getfil (bndc2h6, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, c2h6em )
         c2h6em(:,:,:) = amufac * adv_mass(c2h6_ndx) * c2h6em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... C3H8
!-----------------------------------------------------------------------
      if( c3h8_em ) then
         species  = 'C3H8'
         mw       = 3.*12.011                       ! g C / mole
         call getfil (bndc3h8, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, c3h8em )
         c3h8em(:,:,:) = amufac * adv_mass(c3h8_ndx) * c3h8em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... C2H4
!-----------------------------------------------------------------------
      if( c2h4_em ) then
         species  = 'C2H4'
         mw       = 2.*12.011                       ! g C / mole
         call getfil (bndc2h4, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, c2h4em )
         c2h4em(:,:,:) = amufac * adv_mass(c2h4_ndx) * c2h4em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... C3H6
!-----------------------------------------------------------------------
      if( c3h6_em ) then
         species  = 'C3H6'
         mw       = 3.*12.011                       ! g C / mole
         call getfil (bndc3h6, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, c3h6em )
         c3h6em(:,:,:) = amufac * adv_mass(c3h6_ndx) * c3h6em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... ACET
!-----------------------------------------------------------------------
      if( ch3coch3_em ) then
         species  = 'ACET'
         mw       = adv_mass(ch3coch3_ndx)           ! g / mole
         call getfil (bndacet, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, acetem )
         acetem(:,:,:) = amufac * adv_mass(ch3coch3_ndx) * acetem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... CH2O
!-----------------------------------------------------------------------
      if( ch2o_em ) then
         species  = 'CH2O'
         mw       = adv_mass(ch2o_ndx)               ! g / mole
         call getfil (bndch2o, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, ch2oem )
         ch2oem(:,:,:)  = amufac * adv_mass(ch2o_ndx) * ch2oem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... C4H10
!-----------------------------------------------------------------------
      if( c4h10_em ) then
         species  = 'C4H10'
         mw       = 4.*12.011                        ! g C / mole
         call getfil (bndnmv, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, c4h10em )
         c4h10em(:,:,:) = amufac * adv_mass(c4h10_ndx) * c4h10em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... ISO
!-----------------------------------------------------------------------
      if( isop_em ) then
         species  = 'ISOP'
         mw       = 5.*12.011                      ! g C / mole
         call getfil (bndisoprene, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, isoem )
         isoem(:,:,:) = amufac * adv_mass(isop_ndx) * isoem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... TERP
!-----------------------------------------------------------------------
      if( c10h16_em ) then
         species  = 'TERP'
         mw       = 10.*12.011                     ! g C / mole
         call getfil (bndterpense, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'TgC/y', ncid, terpem )
         terpem(:,:,:) = amufac * adv_mass(c10h16_ndx) * terpem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... N2O
!-----------------------------------------------------------------------
      if( n2o_em ) then
         species  = 'N2O'
         mw       = adv_mass(n2o_ndx)               ! g / mole
         call getfil (bndn2o , locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, n2oem )
         n2oem(:,:,:) = amufac * adv_mass(n2o_ndx) * n2oem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... H2
!-----------------------------------------------------------------------
      if( h2_em ) then
         species  = 'H2'
         mw       = adv_mass(h2_ndx)               ! g / mole
         call getfil (bndh2, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, h2em )
         h2em(:,:,:) = amufac * adv_mass(h2_ndx) * h2em(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... CH3OH
!-----------------------------------------------------------------------
      if( ch3oh_em ) then
         species  = 'CH3OH'
         mw       = adv_mass(ch3oh_ndx)               ! g / mole
         call getfil (bndch3oh, locfn)
         call wrap_open(locfn,0,ncid)
         call read_emis( plonl, platl, pplon, species, mw, 'Tg/y', ncid, ch3ohem )
         ch3ohem(:,:,:) = amufac * adv_mass(ch3oh_ndx) * ch3ohem(:,:,:)
      end if

!-----------------------------------------------------------------------
! 	... NPP
!-----------------------------------------------------------------------
      species  = 'NPP'
      mw       = 1.               ! g / mole
      call getfil (bndnpp, locfn)
      call wrap_open(locfn,0,ncid)
      call read_emis( plonl, platl, pplon, species, mw, '', ncid, npp )

      endif

      allocate( noxembdy(pcols,begchunk:endchunk,12), ch4embdy(pcols,begchunk:endchunk,12), &
                ch2oembdy(pcols,begchunk:endchunk,12), coembdy(pcols,begchunk:endchunk,12), &
                c3h6embdy(pcols,begchunk:endchunk,12), isoembdy(pcols,begchunk:endchunk,12), &
                terpembdy(pcols,begchunk:endchunk,12), c2h4embdy(pcols,begchunk:endchunk,12), &
                c2h6embdy(pcols,begchunk:endchunk,12), c4h10embdy(pcols,begchunk:endchunk,12), &
                n2oembdy(pcols,begchunk:endchunk,12), c3h8embdy(pcols,begchunk:endchunk,12), &
                acetembdy(pcols,begchunk:endchunk,12), h2embdy(pcols,begchunk:endchunk,12), &
                ch3ohembdy(pcols,begchunk:endchunk,12), nppbdy(pcols,begchunk:endchunk,12), &
		stat=astat )
      if( astat/= 0 ) then
	 write(*,*) 'srf_emis_inti: failed to allocate wrking arrays; error = ',astat
	 call endrun
      end if

      call scatter_field_to_chunk(1,1,12,plonl,noxem,noxembdy)
      call scatter_field_to_chunk(1,1,12,plonl,ch4em,ch4embdy)
      call scatter_field_to_chunk(1,1,12,plonl,ch2oem,ch2oembdy)
      call scatter_field_to_chunk(1,1,12,plonl,coem,coembdy)
      call scatter_field_to_chunk(1,1,12,plonl,c3h6em,c3h6embdy)
      call scatter_field_to_chunk(1,1,12,plonl,isoem,isoembdy)
      call scatter_field_to_chunk(1,1,12,plonl,terpem,terpembdy)
      call scatter_field_to_chunk(1,1,12,plonl,c2h4em,c2h4embdy)
      call scatter_field_to_chunk(1,1,12,plonl,c2h6em,c2h6embdy)
      call scatter_field_to_chunk(1,1,12,plonl,c4h10em,c4h10embdy)
      call scatter_field_to_chunk(1,1,12,plonl,n2oem,n2oembdy)
      call scatter_field_to_chunk(1,1,12,plonl,c3h8em,c3h8embdy)
      call scatter_field_to_chunk(1,1,12,plonl,acetem,acetembdy)
      call scatter_field_to_chunk(1,1,12,plonl,h2em,h2embdy)
      call scatter_field_to_chunk(1,1,12,plonl,ch3ohem,ch3ohembdy)
      call scatter_field_to_chunk(1,1,12,plonl,npp,nppbdy)

      end subroutine srf_emis_inti

      subroutine read_emis( plonl, platl, pplon, species, mw, units, ncid, emis )
!-----------------------------------------------------------------------
!       ... Read surface emissions from NetCDF file
!-----------------------------------------------------------------------

      use gchp_netcdf

      use gchp_const_mozart,  only : phi, lam, d2r

      use abortutils,     only: endrun
      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl, platl, pplon
      integer, intent(in) :: ncid
      real, intent(in)    :: mw
      real, dimension(plonl,platl,12), intent(out) :: emis
      character(len=*), intent(in) :: species, units

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: j, jl, ju, m               ! Indices
      integer :: jlim(2)
      integer :: ierr
      integer :: vid, nlat, nlon, nmonth, nvars, ndims, attlen
      integer :: dimid_lat, dimid_lon, dimid_month
      integer, allocatable :: dimid(:)
      integer :: gndx
      real, allocatable :: lat(:), lon(:)
      real :: total, total_wrk, scale_factor
      real :: wrk(plonl,platl,12)
      character(len=NF_MAX_NAME) :: varname
      character(len=100) :: var_longname

!-----------------------------------------------------------------------
!       ... Get grid dimensions from file
!-----------------------------------------------------------------------
      call handle_ncerr( NF_INQ_DIMID( ncid, 'lat', dimid_lat ), &
                         'read_emis: Failed to find dimension lat for species ' // trim(species) )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lat, nlat ), &
                         'read_emis: Failed to get length of dimension lat for species ' // &
                         trim(species) )
      allocate( lat(nlat), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'read_emis: lat allocation error = ',ierr
         call endrun
      end if
      call handle_ncerr( NF_INQ_VARID( ncid, 'lat', vid ), &
                         'read_emis: Failed to find variable lat for species ' // trim(species) )
      call handle_ncerr( NF_GET_VAR_DOUBLE( ncid, vid, lat ), &
                         'read_emis: Failed to read variable lat for species ' // trim(species) )
      lat(:nlat) = lat(:nlat) * d2r
 
      call handle_ncerr( NF_INQ_DIMID( ncid, 'lon', dimid_lon ), &
                         'read_emis: Failed to find dimension lon for species ' // trim(species) )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lon, nlon ), &
                         'read_emis: Failed to find dimension lon for species ' // trim(species) )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lon, nlon ), &
                         'read_emis: Failed to get length of dimension lon for species ' // &
                         trim(species) )
      allocate( lon(nlon), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'read_emis: lon allocation error = ',ierr
         call endrun
      end if
      call handle_ncerr( NF_INQ_VARID( ncid, 'lon', vid ), &
                         'read_emis: Failed to find variable lon for species ' // trim(species) )
      call handle_ncerr( NF_GET_VAR_DOUBLE( ncid, vid, lon ), &
                         'read_emis: Failed to read variable lon for species ' // trim(species) )
      lon(:nlon) = lon(:nlon) * d2r
 
      call handle_ncerr( NF_INQ_DIMID( ncid, 'month', dimid_month ), &
                         'read_emis: Failed to find dimension month for species ' // trim(species) )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_month, nmonth ), &
                         'read_emis: Failed to get length of dimension month for species ' // &
                         trim(species) )
      if( nmonth /= 12 ) then
         write(*,*) 'read_emis: Error, number of months = ',nmonth,' for species ' // trim(species)
         call endrun
      end if

      call handle_ncerr( NF_INQ_NVARS( ncid, nvars ), &
                         'read_emis: Failed to get number of variables for species ' // &
                         trim(species) )

!-----------------------------------------------------------------------
!       ... Set total emissions to zero
!-----------------------------------------------------------------------
      emis(:,:,:) = 0.
!-----------------------------------------------------------------------
!       ... Loop over file variables
!-----------------------------------------------------------------------
      write(*,*) 'Annual average ',trim(species),' emissions ',trim(units)
      do vid = 1,nvars
         ierr = NF_INQ_VARNAME( ncid, vid, varname )
         if( ierr /= 0 ) then
            write(*,*) 'read_emis: Failed to get name of variable # ',vid, &
                '; species=' // trim(species)
            call endrun
         end if
         call handle_ncerr( NF_INQ_VARNDIMS( ncid, vid, ndims ), &
                            'read_emis: Failed to get number of dimensions for ' // &
                            'variable ' // trim(varname) // ', species=' // trim(species) )
         if( ndims < 3 ) then
            cycle
         elseif( ndims > 3 ) then
            write(*,*) 'read_emis: Skipping variable ', trim(varname),', ndims = ',ndims, &
                       ' , species=',trim(species)
            cycle
         end if
!-----------------------------------------------------------------------
!       ... Check order of dimensions. Must be (lon, lat, month).
!-----------------------------------------------------------------------
         allocate( dimid( ndims ), stat=ierr )
         if( ierr /= 0 ) then
            write(*,*) 'read_emis: dimid allocation error = ',ierr
            call endrun
         end if
         call handle_ncerr( NF_INQ_VARDIMID( ncid, vid, dimid ), &
                            'read_emis: Failed to get dimension IDs for variable ' // &
                            trim(varname)  // ', species=' // trim(species) )
         if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. dimid(3) /= dimid_month ) then
            write(*,*) 'read_emis: Dimensions in wrong order for variable ',trim(varname)
            write(*,*) '...      Expecting (lon, lat, month)'
            call endrun
         end if
         deallocate( dimid, stat=ierr )
         if( ierr /= 0 ) then
            write(*,*) 'read_emis: Failed to deallocate dimid, ierr = ',ierr
            call endrun
         end if
!-----------------------------------------------------------------------
!       ... Check for long_name of variable
!-----------------------------------------------------------------------
         var_longname = ' '
         ierr = NF_GET_ATT_TEXT( ncid, vid, 'long_name', var_longname )
         if( ierr /= 0 ) then
            var_longname = trim( varname )
         end if
!-----------------------------------------------------------------------
!       ... Read data from this variable
!-----------------------------------------------------------------------
         call handle_ncerr( NF_GET_VARA_DOUBLE( ncid, vid, &
                                                (/ 1, 1, 1/), &                 ! start
                                                (/ nlon, nlat, 12 /), &  ! count
                                                wrk), &
                            'read_emis: Failed to read variable ' // trim( varname ) // &
                            ', species=' // trim(species) )
!-----------------------------------------------------------------------
!       ... Check for scale_factor
!-----------------------------------------------------------------------
         ierr = NF_GET_ATT_DOUBLE( ncid, vid, 'scale_factor', scale_factor )
         if( ierr /= NF_NOERR  ) then
            scale_factor = 1.
         end if

         total = 0.
         do m = 1,12
           do j = 1,plat
             total = total + sf(j) * SUM( wrk( :plon, j, m ) )
           enddo
           emis (:,:,m) = emis(:,:,m) + wrk(:,:,m)
         enddo
!-----------------------------------------------------------------------
!       ... Get global emission from this source
!-----------------------------------------------------------------------
         total = total * factor * mw / 12. ! convert from molec/s to Tg/y, divide by 12 (months)
 
            write(*,'(2a10,1x,a50,1x,f10.3)') trim(species), trim(varname), trim(var_longname), total
      end do  ! (vid)
!-----------------------------------------------------------------------
!       ... Close NetCDF file
!-----------------------------------------------------------------------
      call handle_ncerr( NF_CLOSE( ncid ), &
                         'read_emis: Failed to close NetCDF file, species=' // trim(species) )

!-----------------------------------------------------------------------
!       ... Get total global emission for this species
!-----------------------------------------------------------------------
      total = 0.
      do j = 1,platl
         total = total + sf(j) * SUM( emis( :plonl, j, :12 ) )
      end do
      total = total * factor * mw / 12. ! convert from molec/s to Tg/y, divide by 12 (months)
 
      write(*,'(2a10,1x,a50,1x,f10.3)') trim(species), 'TOTAL', 'TOTAL', total

      end subroutine read_emis


      end module gchp_srf_emis

