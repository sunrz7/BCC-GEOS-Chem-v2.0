!xiaolu
!---------------------------------------
!Purpose: Lightning NO emission
!from MOZART2 mo_hook.F90

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/10/08
!---------------------------------------

#include <misc.h>
#include <params.h>

      module gchp_hook

      use shr_kind_mod, only: r8 => shr_kind_r8
#if ( defined SPMD )
  use mpishorthand      , only: mpicom, mpir8
#endif

      use gc_grid, only : plev

      implicit none

      private
      save

      public  :: gc_hook_inti,gc_hook, gc_radon, factor

  !   real(r8), allocatable :: prod_no(:,:,:)
  !   real(r8), allocatable :: prod_no_col(:,:)
  !   real(r8), allocatable :: flash_freq(:,:)

      real(r8) :: csrf
!----------------------------------------------------------------------
!	set global lightning nox scaling factor
!----------------------------------------------------------------------
      real(r8) :: factor = 1.             ! user-controlled scaling factor to achieve arbitrary no prod.
      !xiaolu, change to 4 types following Ott et al. (2010) @JGR
      real(r8) :: vdist(16,3)             ! vertical distribution of lightning

      contains

      subroutine gc_hook_inti( lght_no_prd_factor )
!----------------------------------------------------------------------
!       ... initialize the chemistry "hook" routine
!----------------------------------------------------------------------

      use ppgrid,          only : pcols, begchunk, endchunk   !add by zf 2008.11.17
      use gchp_const_mozart, only : twopi, rearth
      use gc_grid,      only : plong => plon
      use abortutils,   only: endrun    !add by zf 2008.07.10

      implicit none

!----------------------------------------------------------------------
!	... dummy args
!----------------------------------------------------------------------
      real(r8), intent(in), optional :: lght_no_prd_factor        ! lightning no production factor

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer :: astat

      if( present( lght_no_prd_factor ) ) then
         if( lght_no_prd_factor /= 1. ) then
            factor = lght_no_prd_factor
         end if
      end if
      csrf = twopi*rearth*rearth/real(plong)    ! rearth in m
!----------------------------------------------------------------------
!       ... vdist(kk,itype) = % of lightning nox between (kk-1) and (kk)
!           km for profile itype
!----------------------------------------------------------------------

     !=================================================================
     !Set lightning NO vertical distribution, xiaolu, 2019/03
      !xiaolu: now change the vertical distribution of Pickering et al. (1998) to Ott et al. (2010)
      vdist(:,1) = (/ 20.1, 2.3, 0.8, 1.5, 3.4, 5.3, 3.6, 3.8, &       ! midlat cont
                       5.4, 6.6, 8.3, 9.6,12.8,10.0, 6.2, 0.3 /)
      vdist(:,2) = (/  5.8, 2.9, 2.6, 2.4, 2.2, 2.1, 2.3, 6.1, &       ! trop marine
                      16.5,14.1,13.7,12.8,12.5, 2.8, 0.9, 0.3 /)
      vdist(:,3) = (/  8.2, 1.9, 2.1, 1.6, 1.1, 1.6, 3.0, 5.8, &       ! trop cont
                       7.6, 9.6,10.5,12.3,11.8,12.5, 8.1, 2.3 /)

      !=================================================================

   !   allocate( prod_no(pcols,plev,begchunk:endchunk),stat=astat )
   !   if( astat /= 0 ) then
!	 write(*,*) 'moz_hook_inti: failed to allocate prod_no; error = ',astat
!	 call endrun
!      end if
   !   allocate( prod_no_col(pcols,begchunk:endchunk),stat=astat )
   !   if( astat /= 0 ) then
   !	 write(*,*) 'moz_hook_inti: failed to allocate prod_no_col; error = ',astat
!	 call endrun
!      end if
!      allocate( flash_freq(pcols,begchunk:endchunk),stat=astat )
!      if( astat /= 0 ) then
!	 write(*,*) 'moz_hook_inti: failed to allocate flash_freq; error = ',astat
!	 call endrun
!      end if

      end subroutine gc_hook_inti

      subroutine gc_hook( ncol, prod_no, cldtop, cldbot, zm, zint, t,     &
			   landfrac, ocnfrac, icefrac, plonl, clat, lchnk )
!----------------------------------------------------------------------
!	... general purpose chemistry "hook" routine.
!           update deposition velocity and sulfur input fields,
!           and calculate lightning nox source & rn emissions.
!----------------------------------------------------------------------

      use pmgrid,          only:  masterproc
      use gchp_const_mozart, only : dayspy, lat25, lat40,lat60, lat70
      use gc_chem_mods,    only : nadv_mass, adv_mass
      use gc_grid,         only : plat, plev, plevp, plong => plon
      use gchp_sim_chm,      only : latwtsbdy
      use history,         only: outfld

      implicit none

!----------------------------------------------------------------------
!	... dummy args
!----------------------------------------------------------------------
      integer, intent(in) :: plonl, lchnk, ncol
      real(r8), intent(in) :: clat(plonl)                          ! latitudes
      real(r8), intent(in) :: cldtop(plonl)     ! cloud top level index
      real(r8), intent(in) :: cldbot(plonl)     ! cloud bottom level index
      real(r8), intent(in) :: landfrac(plonl)        ! orography "flag"
      real(r8), intent(in) :: ocnfrac(plonl)        ! orography "flag"
      real(r8), intent(in) :: icefrac(plonl)        ! orography "flag"
      real(r8), intent(in) :: zm(plonl,plev)    ! geopot height above surface at midpoints (m)
      real(r8), intent(in) :: zint(plonl,plevp) ! geopot height above surface at interfaces (m)
      real(r8), intent(in) :: t(plonl,plev)     ! temperature
      real(r8), intent(out) :: prod_no(plonl,plev)     ! temperature

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer, parameter :: land = 1, ocean = 0
      real(r8), parameter    :: secpyr = dayspy * 8.64e4

      integer :: i, istat, j, jgbl, &
                 cldtind, &         ! level index for cloud top
                 cldbind, &         ! level index for cloud base > 273k
                 surf_type, &
                 node, &
                 file               ! file index
      integer :: k, kk, zlow_ind, zhigh_ind, itype
      real(r8)    :: glob_flashfreq, &  ! global flash frequency [s-1]
                 glob_noprod, &     ! global rate of no production [as tgn/yr]
                 frac_sum, &        !
                 wrk                ! work variable
      real(r8)       :: zlow, zhigh, zlow_scal, zhigh_scal, fraction
      real(r8), dimension( plonl ) :: &
                 oro,  &
                 dchgzone, &        ! depth of discharge zone [km]
                 cldhgt, &          ! cloud top height [km]
                 cgic, &            ! cloud-ground/intracloud discharge ratio
                 flash_energy, &    ! energy of flashes per second
                 glob_prod_no_col   ! global no production rate for diagnostics
     real(r8) :: prod_no_col(plonl)
     real(r8) :: flash_freq(plonl) 
!----------------------------------------------------------------------
! 	... parameters to determine cg/ic ratio [price and rind, 1993]
!----------------------------------------------------------------------
      real(r8), parameter  :: ca = .021, cb = -.648, cc = 7.49, cd = -36.54, ce = 64.09

!----------------------------------------------------------------------
!	lightning no production : initialize ...
!----------------------------------------------------------------------
      prod_no(:,:) = 0.
      prod_no_col(:)  = 0.
      flash_freq(:)   = 0.

      cldhgt(:)     = 0.
      dchgzone(:)   = 0.
      cgic(:)       = 0.
      glob_prod_no_col(:) = 0.
      flash_energy(:) = 0.

      do i = 1, ncol   !plonl
!-----------test--------------------------
        if( landfrac(i) .gt. 0.5 ) then
          oro(i) = 1.
        else
          oro(i) = 0.
        endif
!-----------test--------------------------
!zf        if( abs(max(landfrac(i),ocnfrac(i),icefrac(i)) - landfrac(i)) .lt. 1.e-8) then
!zf          oro(i) = 1.
!zf        else if( abs(max(landfrac(i),ocnfrac(i),icefrac(i)) - ocnfrac(i)) .lt. 1.e-8) then
!zf          oro(i) = 0.
!zf        else
!zf          oro(i) = 2.
!zf        endif
      enddo
!--------------------------------------------------------------------------------
!	... estimate flash frequency and resulting no emissions
!           [price, penner, prather, 1997 (jgr)]
!    lightning only occurs in convective clouds with a discharge zone, i.e.
!    an altitude range where liquid water, ice crystals, and graupel coexist.
!    we test this by examining the temperature at the cloud base.
!    it is assumed that only one thunderstorm occurs per grid box, and its
!    flash frequency is determined by the maximum cloud top height (not the
!    depth of the discharge zone). this is somewhat speculative but yields
!    reasonable results.
!
!       the cg/ic ratio is determined by an empirical formula from price and
!    rind [1993]. the average energy of a cg flash is estimated as 6.7e9 j,
!    and the average energy of a ic flash is assumed to be 1/10 of that value.
!       the no production rate is assumed proportional to the discharge energy
!    with 1e17 n atoms per j. the total number of n atoms is then distributed
!    over the complete column of grid boxes.
!--------------------------------------------------------------------------------
            do i = 1,ncol    !plonl
!--------------------------------------------------------------------------------
! 	... find cloud top and bottom level above 273k
!--------------------------------------------------------------------------------
	       cldtind = nint( cldtop(i) )
               cldbind = nint( cldbot(i) )
               do
                  if( cldbind <= cldtind .or. t(i,cldbind) < 273. ) then
		     exit
		  end if
                  cldbind = cldbind - 1
               end do
	       if( cldtind < plev .and. cldtind > 0 .and. cldtind < cldbind ) then
!--------------------------------------------------------------------------------
!       ... compute cloud top height and depth of charging zone
!--------------------------------------------------------------------------------
	          cldhgt(i) = 1.e-3*max( 0.,zint(i,cldtind) )
                  dchgzone(i) = cldhgt(i)-1.e-3*zm(i,cldbind)
!--------------------------------------------------------------------------------
!       ... compute flash frequency for given cloud top height
!           (flashes storm^-1 min^-1)
!--------------------------------------------------------------------------------
          if( nint( oro(i) ) == land ) then
	             flash_freq(i) = 3.44e-5 * cldhgt(i)**4.9 
	  else
             flash_freq(i) = 6.40e-4 * cldhgt(i)**1.7
          end if
!--------------------------------------------------------------------------------
!       ... compute cg/ic ratio
!           cgic = proportion of cg flashes (=pg from ppp paper)
!--------------------------------------------------------------------------------
                  cgic(i) = 1./((((ca*dchgzone(i) + cb)*dchgzone(i) + cc) &
                                      *dchgzone(i) + cd)*dchgzone(i) + ce)
                  if( dchgzone(i) < 5.5 ) then
		     cgic(i) = 0.
		  end if
                  if( dchgzone(i) > 14. ) then
		     cgic(i) = .02
		  end if
!--------------------------------------------------------------------------------
!       ... compute flash energy (cg*6.7e9 + ic*6.7e8)
!           and convert to total energy per second
!--------------------------------------------------------------------------------
                  flash_energy(i) = cgic(i)*6.7e9 + (1. - cgic(i))*6.7e8
                  flash_energy(i) = flash_energy(i)*flash_freq(i)/60.

!--------------------------------------------------------------------------------
! 	... compute number of n atoms produced per second
!           and convert to n atoms per second per cm2 and apply fudge factor
!--------------------------------------------------------------------------------
               prod_no_col(i) = 1.e17*flash_energy(i) &
                                     /(1.e4*csrf*latwtsbdy(i,lchnk)) * factor
!--------------------------------------------------------------------------------
! 	... compute global no production rate in tgn/yr:
!           tgn per second: * 14.00674 * 1.65979e-24 * 1.e-12
!             nb: 1.65979e-24 = 1/avo
!           tgn per year: * secpyr
!--------------------------------------------------------------------------------
               glob_prod_no_col(i) = 1.e17*flash_energy(i) &
                                        * 14.00674 * 1.65979e-24 * 1.e-12 * secpyr * factor
	       end if
            end do

!--------------------------------------------------------------------------------
! 	... accumulate global total, convert to flashes per second
!--------------------------------------------------------------------------------
      glob_flashfreq = sum( flash_freq(:) )/60.

!--------------------------------------------------------------------------------
! 	... accumulate global no production rate
!--------------------------------------------------------------------------------
      glob_noprod = 0.
      do i = 1,ncol    !plonl
        glob_noprod = glob_noprod + glob_prod_no_col(i)
      end do

      if( masterproc ) then
!         write(*,*) ' '
!         write(*,'('' global flash freq (/s), lightning nox (tgn/y) = '',2f10.4)') &
!                     glob_flashfreq, glob_noprod
!         write(*,*) 'moz_hook : global flash freq (/s), lightning nox (tgn/y) = ', glob_flashfreq, glob_noprod
      end if

      if( glob_noprod > 0. ) then
!--------------------------------------------------------------------------------
!	... distribute production up to cloud top [pickering et al., 1998 (jgr)]
!--------------------------------------------------------------------------------
	       do i = 1,ncol    !plonl
	          cldtind = nint( cldtop(i) )
 	          if( prod_no_col(i) > 0. ) then
	             if( cldhgt(i) > 0. ) then
                        if( abs( clat(i) ) > lat25 ) then         !midlatitude continental
                           itype = 1
                        else if ( nint( oro(i) ) == land ) then
                           itype = 3                              ! tropical continental
                        else
                           itype = 2                              ! topical marine
                        end if
                        frac_sum = 0.
                        do k = cldtind,plev
                           zlow       = zint(i,k+1) * 1.e-3   ! lower interface height (km)
                           zlow_scal  = zlow * 16./cldhgt(i)  ! scale to 16 km convection height
                           zlow_ind   = max( 1,int(zlow_scal)+1 )  ! lowest vdist index to include in layer
                           zhigh      = zint(i,k) * 1.e-3     ! upper interface height (km)
                           zhigh_scal = zhigh * 16./cldhgt(i) ! height (km) scaled to 16km convection height
                           zhigh_ind  = max( 1,min( 16,int(zhigh_scal)+1 ) )  ! highest vdist index to include in layer

                           !---------
                           !xiaolu check vertical distribution of LNO
                           !---------
                           !write(*,*)'xiaolu LNO',k,cldhgt(i),zlow,zhigh,zhigh_scal
                           !write(*,*)'xiaolu LNO',zlow_ind,zhigh_ind
                           do kk = zlow_ind,zhigh_ind
                              fraction = min( zhigh_scal,real(kk) ) &         ! fraction of vdist in this model layer
                                         - max( zlow_scal,real(kk-1) )
                              fraction = max( 0., min( 1.,fraction ) )
                              frac_sum = frac_sum + fraction*vdist(kk,itype)
                              prod_no(i,k) = prod_no(i,k) &         ! sum the fraction of column nox in layer k
                                             + fraction*vdist(kk,itype)*.01
                           end do
			   prod_no(i,k) = prod_no_col(i) * prod_no(i,k) & ! multiply fraction by column amount
                                               / (1.e5*(zhigh - zlow))                   ! and convert to atom n cm^-3 s^-1
                        end do
	             end if
	          end if
	       end do
      end if

      if (masterproc) then
        !write(*,*)'Online lightning NO production'
      endif

!-----------------------------------------------------------------------------------
!  ... output lightning no production to history file
!-----------------------------------------------------------------------------------

      call outfld('LNO_PROD',glob_prod_no_col, plonl, lchnk)
      call outfld('FLASHFRQ',flash_freq(:), plonl, lchnk)
      call outfld('CLDHGT  ',cldhgt, plonl, lchnk)
!     call outfld('DCHGZONE',dchgzone, plonl, lchnk)
      call outfld('DCHGZONE',prod_no_col, plonl, lchnk)
      call outfld('CGIC    ',cgic, plonl, lchnk)

      end subroutine gc_hook

      subroutine gc_radon()

      use pmgrid,          only:  masterproc
      use gchp_const_mozart, only : dayspy, lat60, lat70
      use gc_grid,      only : plat, plon
      use gchp_srf_emis,  only : baseflux, landflux
      use gchp_sim_chm,   only : latwts
      use gchp_surf,      only : ocean, seaice
      use commap,       only : clat
      implicit none

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter    :: secpyr = dayspy * 8.64e4

      integer    :: i,j
      integer    :: surf_type

      real(r8)       :: sflux1, sflux2, sumoc, sumln
      real(r8)       :: total_flux

!--------------------------------------------------------------------------------
!	... form radon surface emission factors
!	    note : radon global emission = 15 kg/yr
!--------------------------------------------------------------------------------
      if( masterproc ) then
      total_flux = 15. / secpyr                                ! kg/s
      sflux1 = 0.
      sflux2 = 0.

      do j =1 ,plat
        if( clat(j) < lat70 .and. clat(j)  > -lat60 ) then
	  sumoc = 0.
	  sumln = 0.
          do i = 1,plon
              if( (ocean(i,j)+seaice(i,j,3)).gt.0.5 ) then
	      sumoc = sumoc + 1.
	    else
	      if( clat(j) >=  lat60 ) then
	        sumln = sumln + .5
              else
	        sumln = sumln + 1.
	      end if
	    end if
	  end do
	  sflux1 = sflux1 + csrf * latwts(j) * sumoc * baseflux
	  sflux2 = sflux2 + csrf * latwts(j) * sumln
        end if
      enddo

      landflux = (total_flux - sflux1) / sflux2

      end if   !  (masterproc) 

#if (defined SPMD )
     call mpibcast( landflux, 1, mpir8, 0, mpicom )
#endif

      end subroutine gc_radon

      end module gchp_hook
