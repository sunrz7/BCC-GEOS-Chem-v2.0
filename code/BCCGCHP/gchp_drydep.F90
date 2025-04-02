!xiaolu
!---------------------------------------
!Purpose: GAS DEPOSITION
!from MOZART2 mo_drydep.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/04/28
!---------------------------------------

#include <misc.h>
      module gchp_drydep
!---------------------------------------------------------------------
!       ... Dry deposition velocity input data and code for netcdf input
!---------------------------------------------------------------------

      use gc_grid, only : pcnstm1

      implicit none
 
      character(len=256), public  ::  bnddvel     !full pathname for Mozart dry deposition velocity

      save

      private
      public  :: dvel_inti , gc_drydep

      real, allocatable :: dvel(:,:,:,:), &   ! depvel array interpolated to model grid
                           dvel_interp(:,:,:) ! depvel array interpolated to grid and time
      real, allocatable :: sea_ice(:,:), &    ! sea ice interpolated to correct time
                           npp_interp(:,:)    ! NPP interpolated to correct time
      integer :: map(pcnstm1)                     ! indices for drydep species
      integer :: nspecies                       ! number of depvel species in input file
      integer :: pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
		 h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
		 ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, &
		 c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
		 no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
		 hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
		 xooh_ndx, ch3cho_ndx, isopooh_ndx
      integer :: o3_tab_ndx = -1, h2o2_tab_ndx = -1, &
                 ch3ooh_tab_ndx = -1, co_tab_ndx = -1, &
                 ch3cho_tab_ndx = -1
      logical :: pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
		 h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
		 ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, ch3cho_dd, c2h5oh_dd, &
		 c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
		 glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd
      logical :: o3_in_tab = .false., h2o2_in_tab = .false., ch3ooh_in_tab = .false., &
                 co_in_tab = .false., ch3cho_in_tab = .false.


      contains

      subroutine dvel_inti( platl, plonl )
!---------------------------------------------------------------------------
!       ... Initialize module, depvel arrays, and a few other variables.
!           The depvel fields will be linearly interpolated to the correct time
!---------------------------------------------------------------------------

      use gchp_netcdf
      use gchp_const_mozart, only : d2r
      use gc_grid,         only : plon, plat
      use gchp_tracname,     only : tracnam
      use gchp_chem_utls,    only : get_spc_ndx, has_drydep
      use gchp_charutl,      only : glc

!add by zf 2008.08.13
      use ppgrid,           only: pcols, begchunk, endchunk
      use pmgrid,           only: masterproc
      use ioFileMod,        only: getfil
      use abortutils,       only: endrun
      use error_messages,   only: handle_ncerr
      use phys_grid,        only: scatter_field_to_chunk

#if ( defined SPMD )
  use mpishorthand,         only: mpicom, mpiint
#endif

!end zf
      implicit none

!---------------------------------------------------------------------------
!       ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: platl
      integer, intent(in) :: plonl

!---------------------------------------------------------------------------
!       ... Local variables
!---------------------------------------------------------------------------
      integer :: ncid, vid, vid_dvel, nlat, nlon, nmonth, ndims
      integer :: dimid_lat, dimid_lon, dimid_species, dimid_time
      integer :: dimid(4), count(4), start(4)
      integer :: jlim(2), jl, ju
      integer :: gndx, m, ispecies, nchar, ierr, i
      integer :: k,j    !zf
      real :: scale_factor
      real, dimension(plon,platl) :: wrk2d
      real, allocatable :: dvel_lats(:), dvel_lons(:)
      real, allocatable :: dvel_in(:,:,:,:)                          ! input depvel array
      character(len=50) :: units
      character(len=80) :: filename, lpath, rpath
      character(len=20), allocatable :: species_names(:)             ! names of depvel species
      logical :: found

      character(len=256) :: locfn    ! netcdf local filename to open   2008.08.13
!---------------------------------------------------------------------------
!       ... Setup species maps
!---------------------------------------------------------------------------
      pan_ndx  = get_spc_ndx( 'PAN')
      mpan_ndx = get_spc_ndx( 'MPAN')
      no2_ndx  = get_spc_ndx( 'NO2')
      hno3_ndx = get_spc_ndx( 'HNO3')
      co_ndx   = get_spc_ndx( 'CO')
      o3_ndx   = get_spc_ndx( 'O3')
      if( o3_ndx < 1 ) then
         o3_ndx = get_spc_ndx( 'OX')
      end if
      h2o2_ndx     = get_spc_ndx( 'H2O2')
      onit_ndx     = get_spc_ndx( 'ONIT')
      onitr_ndx    = get_spc_ndx( 'ONITR')
      ch4_ndx      = get_spc_ndx( 'CH4')
      ch2o_ndx     = get_spc_ndx( 'CH2O')
      ch3ooh_ndx   = get_spc_ndx( 'CH3OOH')
      ch3cho_ndx   = get_spc_ndx( 'CH3CHO')
      ch3cocho_ndx = get_spc_ndx( 'CH3COCHO')
      pooh_ndx     = get_spc_ndx( 'POOH')
      ch3coooh_ndx = get_spc_ndx( 'CH3COOOH')
      c2h5ooh_ndx  = get_spc_ndx( 'C2H5OOH')
      c3h7ooh_ndx  = get_spc_ndx( 'C3H7OOH')
      rooh_ndx     = get_spc_ndx( 'ROOH')
      ch3coch3_ndx = get_spc_ndx( 'CH3COCH3')
      no_ndx       = get_spc_ndx( 'NO')
      ho2no2_ndx   = get_spc_ndx( 'HO2NO2')
      glyald_ndx   = get_spc_ndx( 'GLYALD')
      hyac_ndx     = get_spc_ndx( 'HYAC')
      ch3oh_ndx    = get_spc_ndx( 'CH3OH')
      c2h5oh_ndx   = get_spc_ndx( 'C2H5OH')
      macrooh_ndx  = get_spc_ndx( 'MACROOH')
      isopooh_ndx  = get_spc_ndx( 'ISOPOOH')
      xooh_ndx     = get_spc_ndx( 'XOOH')
      hydrald_ndx  = get_spc_ndx( 'HYDRALD')
      h2_ndx       = get_spc_ndx( 'H2')
      Pb_ndx       = get_spc_ndx( 'Pb')
      o3s_ndx      = get_spc_ndx( 'O3S')
      o3inert_ndx  = get_spc_ndx( 'O3INERT')

      pan_dd  = has_drydep( 'PAN')
      mpan_dd = has_drydep( 'MPAN')
      no2_dd  = has_drydep( 'NO2')
      hno3_dd = has_drydep( 'HNO3')
      co_dd   = has_drydep( 'CO')
      o3_dd   = has_drydep( 'O3')
      if( .not. o3_dd ) then
         o3_dd = has_drydep( 'OX')
      end if
      h2o2_dd     = has_drydep( 'H2O2')
      onit_dd     = has_drydep( 'ONIT')
      onitr_dd    = has_drydep( 'ONITR')
      ch4_dd      = has_drydep( 'CH4')
      ch2o_dd     = has_drydep( 'CH2O')
      ch3ooh_dd   = has_drydep( 'CH3OOH')
      ch3cho_dd   = has_drydep( 'CH3CHO')
      c2h5oh_dd   = has_drydep( 'C2H5OH')
      ch3cocho_dd = has_drydep( 'CH3COCHO')
      pooh_dd     = has_drydep( 'POOH')
      ch3coooh_dd = has_drydep( 'CH3COOOH')
      c2h5ooh_dd  = has_drydep( 'C2H5OOH')
      c3h7ooh_dd  = has_drydep( 'C3H7OOH')
      rooh_dd     = has_drydep( 'ROOH')
      ch3coch3_dd = has_drydep( 'CH3COCH3')
      glyald_dd   = has_drydep( 'GLYALD')
      hyac_dd     = has_drydep( 'HYAC')
      ch3oh_dd    = has_drydep( 'CH3OH')
      macrooh_dd  = has_drydep( 'MACROOH')
      isopooh_dd  = has_drydep( 'ISOPOOH')
      xooh_dd     = has_drydep( 'XOOH')
      hydrald_dd  = has_drydep( 'HYDRALD')
      h2_dd       = has_drydep( 'H2')
      Pb_dd       = has_drydep( 'Pb')
      o3s_dd      = has_drydep( 'O3S')
      o3inert_dd  = has_drydep( 'O3INERT')

    if (masterproc) then
      write(*,*) 'xiaolu check dvel_inti: diagnostics'
      write(*,'(10i5)') pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
		 h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
		 ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, &
		 c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
		 no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
		 hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
		 xooh_ndx, ch3cho_ndx, isopooh_ndx
      write(*,*) pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
		 h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
		 ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, ch3cho_dd, c2h5oh_dd, &
		 c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
		 glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd
     end if   !masterproc

!---------------------------------------------------------------------------
!       ... Open NetCDF file
!---------------------------------------------------------------------------
      if( masterproc ) then

      call getfil (bnddvel, locfn)
      call wrap_open(locfn,0,ncid)

!---------------------------------------------------------------------------
!       ... Get variable ID for dep vel array
!---------------------------------------------------------------------------
      call handle_ncerr( NF_INQ_VARID( ncid, 'dvel', vid_dvel ), 'dvel_inti: dvel not found in depvel input file' )

!---------------------------------------------------------------------------
!       ... Inquire about dimensions
!---------------------------------------------------------------------------
      call handle_ncerr( NF_INQ_DIMID( ncid, 'lon', dimid_lon ), 'dvel_inti: getting lon dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lon, nlon ), 'dvel_inti: getting nlon' )
      call handle_ncerr( NF_INQ_DIMID( ncid, 'lat', dimid_lat ), 'dvel_inti: getting lat dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_lat, nlat ), 'dvel_inti: getting nlat' )
      call handle_ncerr( NF_INQ_DIMID( ncid, 'species', dimid_species ), 'dvel_inti: getting species dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_species, nspecies ), 'dvel_inti: getting nspecies' )
      call handle_ncerr( NF_INQ_DIMID( ncid, 'time', dimid_time ), 'dvel_inti: getting time dimension ID' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid_time, nmonth ), 'dvel_inti: getting nmonth' )
      write(*,*) 'dvel_inti: dimensions (nlon,nlat,nspecies,nmonth) = ',nlon,nlat,nspecies,nmonth

!---------------------------------------------------------------------------
!       ... Check dimensions of dvel variable. Must be (lon, lat, species, month).
!---------------------------------------------------------------------------
      call handle_ncerr( NF_INQ_VARNDIMS( ncid, vid_dvel, ndims ), &
                         'dvel_inti: Failed to get ndims for dvel' )
      if( ndims /= 4 ) then
         write(*,*) 'dvel_inti: dvel has ',ndims,' dimensions. Expecting 4.'
         call endrun
      end if
      call handle_ncerr( NF_INQ_VARDIMID( ncid, vid_dvel, dimid ), &
                         'dvel_inti: Failed to get dimension IDs for dvel' )
      if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. &
          dimid(3) /= dimid_species .or. dimid(4) /= dimid_time ) then
         write(*,*) 'dvel_inti: Dimensions in wrong order for dvel'
         write(*,*) '...      Expecting (lon, lat, species, month)'
         call endrun
      end if

!---------------------------------------------------------------------------
!       ... Allocate depvel lats, lons and read
!---------------------------------------------------------------------------
      allocate( dvel_lats(nlat), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate dvel_lats vector'
         call endrun
      end if
      allocate( dvel_lons(nlon), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate dvel_lons vector'
         call endrun
      end if
      call handle_ncerr( NF_INQ_VARID( ncid, 'lat', vid ), 'dvel_inti: lat dim not in dvel file' )
      call handle_ncerr( NF_GET_VAR_DOUBLE( ncid, vid, dvel_lats ), 'dvel_inti: getting dvel lats' )
      call handle_ncerr( NF_INQ_VARID( ncid, 'lon', vid ), 'dvel_inti: lon dim not in dvel file' )
      call handle_ncerr( NF_GET_VAR_DOUBLE( ncid, vid, dvel_lons ), 'dvel_inti: getting dvel lons' )

!---------------------------------------------------------------------------
!       ... Set the transform from inputs lats to simulation lats
!---------------------------------------------------------------------------

      dvel_lats(:nlat) = d2r * dvel_lats(:nlat)
      dvel_lons(:nlon) = d2r * dvel_lons(:nlon)

      deallocate( dvel_lats, dvel_lons )

      end if             !masterproc

!---------------------------------------------------------------------------
!     	... Allocate dvel and read data from file
!---------------------------------------------------------------------------
#if (defined SPMD )
    call mpibcast( nlon, 1, mpiint, 0, mpicom )
    call mpibcast( nlat, 1, mpiint, 0, mpicom )
    call mpibcast( nmonth, 1, mpiint, 0, mpicom )
    call mpibcast( nspecies, 1, mpiint, 0, mpicom )
#endif
      allocate( dvel_in(nlon,nlat,nspecies,nmonth), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate dvel_in'
         call endrun
      end if

      if( masterproc ) then

      start = (/ 1, 1, 1, 1 /)
      count = (/ nlon, nlat, nspecies, nmonth /)
      call handle_ncerr( NF_GET_VARA_DOUBLE( ncid, vid_dvel, start, count, dvel_in ), &
                         'dvel_inti: getting dvel' )

!---------------------------------------------------------------------------
!     	... Check units of deposition velocity. If necessary, convert to cm/s.
!---------------------------------------------------------------------------
      units(:) = ' '
      call handle_ncerr( NF_GET_ATT_TEXT( ncid, vid_dvel, 'units', units  ), 'dvel_inti: dvel units not found' )
      if( units(:glc(units)) == 'm/s' ) then
         scale_factor = 100.
      elseif( units(:glc(units)) == 'cm/s' ) then
         scale_factor = 1.
      else
         scale_factor = 1.
      end if

      dvel_in(:,:,:,:) = dvel_in(:,:,:,:) * scale_factor

      end if             !masterproc
!---------------------------------------------------------------------------
!     	... Regrid deposition velocities
!---------------------------------------------------------------------------
      allocate( dvel(pcols,begchunk:endchunk,nspecies,nmonth), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate dvel'
         call endrun
      end if
      do i = 1, nmonth
        call scatter_field_to_chunk(1,1,nspecies,nlon,dvel_in(:,:,:,i),dvel(:,:,:,i))
      enddo

!---------------------------------------------------------------------------
!     	... Read in species names and determine mapping to tracer numbers
!---------------------------------------------------------------------------
      allocate( species_names(nspecies), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: species_names allocation error = ',ierr
         call endrun
      end if

      map(:) = 0

      if( masterproc ) then
      deallocate( dvel_in )

      call handle_ncerr( NF_INQ_VARID( ncid, 'species_name', vid ), &
                         'dvel_inti: Getting species_name id' )
      call handle_ncerr( NF_INQ_VARNDIMS( ncid, vid, ndims ), &
                         'dvel_inti: Getting number of dimensions for species_name' )
      call handle_ncerr( NF_INQ_VARDIMID( ncid, vid, dimid ), &
                         'dvel_inti: Getting dimensions for species_name' )
      call handle_ncerr( NF_INQ_DIMLEN( ncid, dimid(1), nchar ), &
                         'dvel_inti: Getting dimension length' )
      do ispecies = 1,nspecies
         start(:2) = (/ 1, ispecies /)
         count(:2) = (/ nchar, 1 /)
         species_names(ispecies)(:) = ' '
         call handle_ncerr( NF_GET_VARA_TEXT( ncid, vid, start(:2), count(:2), species_names(ispecies) ), &
                            'dvel_inti: Getting species names' )
    
         if( species_names(ispecies) == 'O3' ) then
	    o3_in_tab  = .true.
	    o3_tab_ndx = ispecies
         else if( species_names(ispecies) == 'H2O2' ) then
	    h2o2_in_tab  = .true.
	    h2o2_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CH3OOH' ) then
	    ch3ooh_in_tab  = .true.
	    ch3ooh_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CO' ) then
	    co_in_tab  = .true.
	    co_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CH3CHO' ) then
	    ch3cho_in_tab  = .true.
	    ch3cho_tab_ndx = ispecies
	 end if
         found = .false.
         do m = 1,pcnstm1
            if( species_names(ispecies) == tracnam(m) .or. &
                (species_names(ispecies) == 'O3' .and. tracnam(m) == 'OX') .or. &
                (species_names(ispecies) == 'HNO4' .and. tracnam(m) == 'HO2NO2') ) then
               map(m) = ispecies
               found = .true.
!  write(*,*),'xiaolu check map speices',map(m),species_names(ispecies),tracnam(m)
!  xiaolu:matched,2017/06/16
               exit
            end if
         end do
         if( .not. found ) then
            write(*,*) 'dvel_inti: Warning! DVEL species ',trim(species_names(ispecies)),' not found'
         endif
      end do
      deallocate( species_names )

      end  if    !masterproc

#if (defined SPMD )
    do m = 1, pcnstm1
      call mpibcast( map(m), 1, mpiint, 0, mpicom )
    end do
    
    call mpibcast( o3_tab_ndx, 1, mpiint, 0, mpicom )
    call mpibcast( h2o2_tab_ndx, 1, mpiint, 0, mpicom )
    call mpibcast( ch3ooh_tab_ndx, 1, mpiint, 0, mpicom )
    call mpibcast( co_tab_ndx, 1, mpiint, 0, mpicom )
    call mpibcast( ch3cho_tab_ndx, 1, mpiint, 0, mpicom )

    call mpibcast( o3_in_tab, 1, mpiint, 0, mpicom )
    call mpibcast( h2o2_in_tab, 1, mpiint, 0, mpicom )
    call mpibcast( ch3ooh_in_tab, 1, mpiint, 0, mpicom )
    call mpibcast( co_in_tab, 1, mpiint, 0, mpicom )
    call mpibcast( ch3cho_in_tab, 1, mpiint, 0, mpicom )
#endif
!---------------------------------------------------------------------------
!     	... Allocate dvel_interp array
!---------------------------------------------------------------------------
      allocate( dvel_interp(pcols,begchunk:endchunk,nspecies),stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate dvel_interp; error = ',ierr
         call endrun
      end if
      allocate( npp_interp(pcols,begchunk:endchunk),stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate npp_interp; error = ',ierr
         call endrun
      end if
      allocate( sea_ice(pcols,begchunk:endchunk),stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'dvel_inti: Failed to allocate sea_ice; error = ',ierr
         call endrun
      end if

      end subroutine dvel_inti

      subroutine interpdvel( plonl, calday, lat )
!---------------------------------------------------------------------------
!       ... Interpolate the fields whose values are required at the
!           begining of a timestep.
!---------------------------------------------------------------------------

      use gchp_calendar, only : caldayr 
      use gchp_surf,     only : seaicebdy 
      use gchp_srf_emis, only : nppbdy

      implicit none

!---------------------------------------------------------------------------
!       ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: plonl        ! latitude tile dimension
      integer, intent(in) :: lat          ! local latitude index
      real, intent(in) :: calday          ! Interpolate the input data to calday

!---------------------------------------------------------------------------
!       ... Local variables
!---------------------------------------------------------------------------
      integer  ::  m, last, next, i
      integer  ::  dates(12) = (/ 116, 214, 316, 415,  516,  615, &
                                  716, 816, 915, 1016, 1115, 1216 /)
      real :: calday_loc, last_days, next_days
      real, save ::  days(12)
      logical, save  ::  entered = .false.

      if( .not. entered ) then
         do m = 1,12
            days(m) = caldayr( dates(m), 0 )
         end do
         entered = .true.
      end if

      if( calday < days(1) ) then
         next = 1
         last = 12
      else if( calday >= days(12) ) then
         next = 1
         last = 12
      else
         do m = 11,1,-1
            if( calday >= days(m) ) then
               exit
            end if
         end do
         last = m
         next = m + 1
      end if

      last_days  = days( last )
      next_days  = days( next )
      calday_loc = calday

      if( next_days < last_days ) then
         next_days = next_days + 365.
      end if
      if( calday_loc < last_days ) then
         calday_loc = calday_loc + 365.
      end if

      do m = 1,nspecies
         call intp2d( last_days, next_days, calday_loc, &
                      dvel(1:plonl,lat,m,last), &
                      dvel(1:plonl,lat,m,next), &
                      dvel_interp(1:plonl,lat,m), plonl )
      end do


      call intp2d( last_days, next_days, calday_loc, &
                   seaicebdy(1:plonl,lat,last), &
                   seaicebdy(1:plonl,lat,next), &
                   sea_ice(1:plonl,lat),plonl )

      !xiaolu error, seaicebdy->nppbdy,may introduce error
      call intp2d( last_days, next_days, calday_loc, &
                   nppbdy(1:plonl,lat,last), &
                   nppbdy(1:plonl,lat,next), &
                   npp_interp(1:plonl,lat),plonl )
    !    call intp2d( last_days, next_days, calday_loc, &
    !               seaicebdy(1:plonl,lat,last), &
    !               seaicebdy(1:plonl,lat,next), &
    !               npp_interp(1:plonl,lat),plonl )


      end subroutine interpdvel

      subroutine intp2d( t1, t2, tint, f1, f2, fint, plonl )
!-----------------------------------------------------------------------
!       ... Linearly interpolate between f1(t1) and f2(t2) to fint(tint).
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)  :: plonl                                   

      real, intent(in)  :: t1            ! time level of f1
      real, intent(in)  :: t2            ! time level of f2
      real, intent(in)  :: tint          ! interpolant time
      real, intent(in)  :: f1(plonl)     ! field at time t1
      real, intent(in)  :: f2(plonl)     ! field at time t2

      real, intent(out) :: fint(plonl)   ! field at time tint

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      real :: factor

      factor = (tint - t1)/(t2 - t1)
      fint(:) = f1(:) + (f2(:) - f1(:))*factor

      end subroutine intp2d

      subroutine gc_drydep( lat, calday, tsurf, zen_angle, &
                         depvel, dflx, rair, q, p, tv, plonl )
!--------------------------------------------------------
!       ... Form the deposition velocities for this
!           latitude slice
!--------------------------------------------------------

      use gchp_chem_utls, only : get_spc_ndx, has_drydep
      use gchp_surf,      only : oceanbdy
      use gc_grid,      only : pcnstm1, plev
      use pmgrid,           only: masterproc

      implicit none

!--------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------
      integer, intent(in) ::  plonl               ! longitude tile dimension
      integer, intent(in) ::  lat                 ! incoming latitude index
      real, intent(in)    ::  rair                ! gas constant of dry air in J/deg/kg
      real, intent(in)    ::  q(plonl,plev,pcnstm1) ! tracer mmr (kg/kg)
      real, intent(in)    ::  p(plonl)            ! midpoint pressure in surface layer (Pa)
      real, intent(in)    ::  calday              ! time of year in days
      real, intent(in)    ::  tsurf(plonl)        ! surface temperature (K)
      real, intent(in)    ::  zen_angle(plonl)    ! zenith angle (radians)
      real, intent(in)    ::  tv(plonl)           ! virtual temperature in surface layer (K)
      real, intent(inout) ::  dflx(plonl,pcnstm1)   ! flux due to dry deposition (kg/m^2/sec)
      real, intent(out)   ::  depvel(plonl,pcnstm1) ! deposition vel (cm/s)

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: m, spc_ndx, tmp_ndx, i
      real    :: vel(plonl)
      real    :: glace(plonl)
      real    :: open_ocean(plonl)
      real    :: temp_fac(plonl)
      real    :: wrk(plonl)
      real    :: tmp(plonl)
      real    :: o3_tab_dvel(plonl)

!-----------------------------------------------------------------------
!       ... Note the factor 1.e-2 in the wrk array calculation is
!           to transform the incoming dep vel from cm/s to m/s
!-----------------------------------------------------------------------
      wrk(:) =  1.e-2 * p(:) / (rair * tv(:))

!--------------------------------------------------------
!       ... Initialize all deposition velocities to zero
!--------------------------------------------------------
      depvel(:,:) = 0.

!--------------------------------------------------------
!       ... Time interpolate primary depvel array
!           (also seaice and npp)
!--------------------------------------------------------
      call interpdvel( plonl, calday, lat )

     if( o3_in_tab ) then
        o3_tab_dvel(:) = dvel_interp(:,lat,o3_tab_ndx)
     end if
!--------------------------------------------------------
!       ... Set deposition velocities
!--------------------------------------------------------
      do m = 1,pcnstm1
         if( map(m) /= 0 ) then
            depvel(:,m) = dvel_interp(:,lat,map(m))
            dflx(:,m)   = wrk(:) * depvel(:,m) * q(:,plev,m)
!  write(*,*),'xiaolu check in drydep',m,map(m),dflx(:,m)
!  xiaolu, dflx are updated:2017/06/16
         end if
      end do
!>>>>????>>>>>
!xiaolu: the following are not used?
!>>>>>>>>>>>>>
!--------------------------------------------------------
!       ... Set some variables needed for some dvel calculations
!--------------------------------------------------------

      temp_fac(:)   = min( 1., max( 0., (tsurf(:) - 268.) / 5. ) )
      open_ocean(:) = max( 0.,oceanbdy(:,lat) - sea_ice(:,lat) )
      glace(:)      = sea_ice(:,lat) + (1. - oceanbdy(:,lat)) * (1. - temp_fac(:))
      glace(:)      = min( 1.,glace(:) )

!--------------------------------------------------------
!       ... Set pan & mpan
!--------------------------------------------------------
      if( o3_in_tab ) then
         tmp(:) = o3_tab_dvel(:) / 3.
      else
         tmp(:) = 0.
      end if
      if( pan_dd ) then
         if( map(pan_ndx) == 0 ) then
            depvel(:,pan_ndx) = tmp(:)
            dflx(:,pan_ndx)   = wrk(:) * tmp(:) * q(:,plev,pan_ndx)
         end if
      end if
      if( mpan_dd ) then
         if( map(mpan_ndx) == 0 ) then
            depvel(:,mpan_ndx) = tmp(:)
            dflx(:,mpan_ndx)   = wrk(:) * tmp(:) * q(:,plev,mpan_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set no2 dvel
!--------------------------------------------------------
      if( no2_dd ) then
         if( map(no2_ndx) == 0 .and. o3_in_tab ) then
            depvel(:,no2_ndx) = (.6*o3_tab_dvel(:) + .055*oceanbdy(:,lat)) * .9
            dflx(:,no2_ndx)   = wrk(:) * depvel(:,no2_ndx) * q(:,plev,no2_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set hno3 dvel
!--------------------------------------------------------
      tmp(:) = (2. - open_ocean(:)) * (1. - glace(:)) + .05 * glace(:)
      if( hno3_dd ) then
         if( map(hno3_ndx) == 0 ) then
            depvel(:,hno3_ndx) = tmp(:)
            dflx(:,hno3_ndx)   = wrk(:) * tmp(:) * q(:,plev,hno3_ndx)
         else
            tmp(:) = depvel(:,hno3_ndx)
         end if
      end if
      if( onitr_dd ) then
         if( map(onitr_ndx) == 0 ) then
            depvel(:,onitr_ndx) = tmp(:)
            dflx(:,onitr_ndx)   = wrk(:) * tmp(:) * q(:,plev,onitr_ndx)
         end if
      end if
      if( isopooh_dd ) then
         if( map(isopooh_ndx) == 0 ) then
            depvel(:,isopooh_ndx) = tmp(:)
            dflx(:,isopooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,isopooh_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set h2o2 dvel
!--------------------------------------------------------
      if( h2o2_dd ) then
         if( .not. h2o2_in_tab ) then
            if( o3_in_tab ) then
               tmp(:) = .05*glace(:) + oceanbdy(:,lat) - sea_ice(:,lat) &
                                     + (1. - (glace(:) + oceanbdy(:,lat)) + sea_ice(:,lat)) &
                         *max( 1.,1./(.5 + 1./(6.*o3_tab_dvel(:))) )
            else
               tmp(:) = 0.
            end if
         else
            tmp(:) = dvel_interp(:,lat,h2o2_tab_ndx)
         end if
         if( map(h2o2_ndx) == 0 ) then
            depvel(:,h2o2_ndx) = tmp(:)
            dflx(:,h2o2_ndx)   = wrk(:) * tmp(:) * q(:,plev,h2o2_ndx)
         end if
      end if
!--------------------------------------------------------
!       ... Set onit
!--------------------------------------------------------
      if( onit_dd ) then
         if( map(onit_ndx) == 0 ) then
            depvel(:,onit_ndx) = tmp(:)
            dflx(:,onit_ndx)   = wrk(:) * tmp(:) * q(:,plev,onit_ndx)
         end if
      end if
      if( ch3cocho_dd ) then
         if( map(ch3cocho_ndx) == 0 ) then
            depvel(:,ch3cocho_ndx) = tmp(:)
            dflx(:,ch3cocho_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch3cocho_ndx)
         end if
      end if
      if( ch3ooh_in_tab ) then
         tmp(:) = dvel_interp(:,lat,ch3ooh_tab_ndx)
      else
         tmp(:) = .5 * tmp(:)
      end if
      if( ch3ooh_dd ) then
         if( map(ch3ooh_ndx) == 0 ) then
            depvel(:,ch3ooh_ndx) = tmp(:)
            dflx(:,ch3ooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch3ooh_ndx)
         end if
      end if
      if( pooh_dd ) then
         if( map(pooh_ndx) == 0 ) then
            depvel(:,pooh_ndx) = tmp(:)
            dflx(:,pooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,pooh_ndx)
         end if
      end if
      if( ch3coooh_dd ) then
         if( map(ch3coooh_ndx) == 0 ) then
            depvel(:,ch3coooh_ndx) = tmp(:)
            dflx(:,ch3coooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch3coooh_ndx)
         end if
      end if
      if( c2h5ooh_dd ) then
         if( map(c2h5ooh_ndx) == 0 ) then
            depvel(:,c2h5ooh_ndx) = tmp(:)
            dflx(:,c2h5ooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,c2h5ooh_ndx)
         end if
      end if
      if( c3h7ooh_dd ) then
         if( map(c3h7ooh_ndx) == 0 ) then
            depvel(:,c3h7ooh_ndx) = tmp(:)
            dflx(:,c3h7ooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,c3h7ooh_ndx)
         end if
      end if
      if( rooh_dd ) then
         if( map(rooh_ndx) == 0 ) then
            depvel(:,rooh_ndx) = tmp(:)
            dflx(:,rooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,rooh_ndx)
         end if
      end if
      if( macrooh_dd ) then
         if( map(macrooh_ndx) == 0 ) then
            depvel(:,macrooh_ndx) = tmp(:)
            dflx(:,macrooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,macrooh_ndx)
         end if
      end if
      if( xooh_dd ) then
         if( map(xooh_ndx) == 0 ) then
            depvel(:,xooh_ndx) = tmp(:)
            dflx(:,xooh_ndx)   = wrk(:) * tmp(:) * q(:,plev,xooh_ndx)
         end if
      end if
      if( ch3oh_dd ) then
         if( map(ch3oh_ndx) == 0 ) then
            depvel(:,ch3oh_ndx) = tmp(:)
            dflx(:,ch3oh_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch3oh_ndx)
         end if
      end if
      if( c2h5oh_dd ) then
         if( map(c2h5oh_ndx) == 0 ) then
            depvel(:,c2h5oh_ndx) = tmp(:)
            dflx(:,c2h5oh_ndx)   = wrk(:) * tmp(:) * q(:,plev,c2h5oh_ndx)
         end if
      end if

      if( o3_in_tab ) then
         tmp(:) = o3_tab_dvel(:)
      else
         tmp(:) = 0.
      end if
      if( ch2o_dd ) then
         if( map(ch2o_ndx) == 0 ) then
            depvel(:,ch2o_ndx) = tmp(:)
            dflx(:,ch2o_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch2o_ndx)
         end if
      end if
      if( hydrald_dd ) then
         if( map(hydrald_ndx) == 0 ) then
            depvel(:,hydrald_ndx) = tmp(:)
            dflx(:,hydrald_ndx)   = wrk(:) * tmp(:) * q(:,plev,hydrald_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set co and related species dep vel
!--------------------------------------------------------
      if( co_in_tab ) then
         tmp(:) = dvel_interp(:,lat,co_tab_ndx)
      else
         tmp(:) = .3333 * npp_interp(:,lat) * ( 1. - glace(:) )
      end if
      if( co_dd ) then
         if( map(co_ndx) == 0 ) then
            depvel(:,co_ndx) = tmp(:)
            dflx(:,co_ndx)   = wrk(:) * tmp(:) * q(:,plev,co_ndx)
         end if
      end if
      if( ch3coch3_dd ) then
         if( map(ch3coch3_ndx) == 0 ) then
            depvel(:,ch3coch3_ndx) = tmp(:)
            dflx(:,ch3coch3_ndx)   = wrk(:) * tmp(:) * q(:,plev,ch3coch3_ndx)
         end if
      end if
      if( hyac_dd ) then
         if( map(hyac_ndx) == 0 ) then
            depvel(:,hyac_ndx) = tmp(:)
            dflx(:,hyac_ndx)   = wrk(:) * tmp(:) * q(:,plev,hyac_ndx)
         end if
      end if
      if( h2_dd ) then
         if( map(h2_ndx) == 0 ) then
            depvel(:,h2_ndx) = tmp(:) * 1.5                ! Hough(1991)
            dflx(:,h2_ndx)   = wrk(:) * depvel(:,h2_ndx) * q(:,plev,h2_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set glyald
!--------------------------------------------------------
      if( glyald_dd ) then
         if( map(glyald_ndx) == 0 ) then
            if( ch3cho_dd ) then
               depvel(:,glyald_ndx) = depvel(:,ch3cho_ndx)
            else if( ch3cho_in_tab ) then
               depvel(:,glyald_ndx) = dvel_interp(:,lat,ch3cho_tab_ndx)
            else
               depvel(:,glyald_ndx) = 0.
            end if
            dflx(:,glyald_ndx)   = wrk(:) * depvel(:,glyald_ndx) * q(:,plev,glyald_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set Lead deposition
!--------------------------------------------------------
      if( Pb_dd ) then
         if( map(Pb_ndx) == 0 ) then
            depvel(:,Pb_ndx) = oceanbdy(:,lat)  * .05 + (1. - oceanbdy(:,lat)) * .2
            dflx(:,Pb_ndx)   = wrk(:) * depvel(:,Pb_ndx) * q(:,plev,Pb_ndx)
         end if
      end if

!--------------------------------------------------------
!       ... Set diurnal dependence for OX dvel
!--------------------------------------------------------
      if( o3_dd .or. o3s_dd .or. o3inert_dd ) then
         if( o3_dd .or. o3_in_tab ) then
            if( o3_dd ) then
               tmp(:) = max( 1.,sqrt( (depvel(:,o3_ndx) - .2)**3/.27 + 4.*depvel(:,o3_ndx) + .67 ) )
               vel(:) = depvel(:,o3_ndx)
            else if( o3_in_tab ) then
               tmp(:) = max( 1.,sqrt( (o3_tab_dvel(:) - .2)**3/.27 + 4.*o3_tab_dvel(:) + .67 ) )
               vel(:) = o3_tab_dvel(:)
            end if
            where( zen_angle(:) <= 0. )
               vel(:) = vel(:) / tmp(:)
            elsewhere
               vel(:) = vel(:) * tmp(:)
            endwhere
         else
            vel(:) = 0.
         end if
         if( o3_dd ) then
            depvel(:,o3_ndx) = vel(:)
            dflx(:,o3_ndx)   = wrk(:) * vel(:) * q(:,plev,o3_ndx)
         end if
!--------------------------------------------------------
!       ... Set stratospheric O3 deposition
!--------------------------------------------------------
         if( o3s_dd ) then
            depvel(:,o3s_ndx) = vel(:)
            dflx(:,o3s_ndx)   = wrk(:) * vel(:) * q(:,plev,o3s_ndx)
         end if
         if( o3inert_dd ) then
            depvel(:,o3inert_ndx) = vel(:)
            dflx(:,o3inert_ndx)   = wrk(:) * vel(:) * q(:,plev,o3inert_ndx)
         end if
      end if

      end subroutine gc_drydep

      end module gchp_drydep
