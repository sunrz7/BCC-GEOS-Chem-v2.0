!xiaolu
!-------------------
!XIAOLU,2017/06/01
!From mo_chem_utls
!------------------
      module gchp_chem_utls
 
      private
!------------
!MOZART
!      public :: adjh2o, inti_mr_xform, mmr2vmr, vmr2mmr, negtrc, &
!                get_spc_ndx, get_het_ndx, get_extfrc_ndx, &
!                has_drydep, has_srfems, get_rxt_ndx, get_grp_ndx, &
!                get_grp_mem_ndx, chem_utls_inti
!-----------
       public :: get_spc_ndx,has_drydep,has_srfems

      save
    
      integer   ::  ox_ndx, o3_ndx, o1d_ndx, o_ndx
      logical   ::  do_ox

      contains

!      subroutine chem_utls_inti
!!-----------------------------------------------------------------------
!!     ... Initialize the chem utils module
!!-----------------------------------------------------------------------

!      implicit none

!      ox_ndx = get_spc_ndx( 'OX' )
!      if( ox_ndx > 0 ) then
!         o3_ndx  = get_grp_mem_ndx( 'O3' )
!         o1d_ndx = get_grp_mem_ndx( 'O1D' )
!         o_ndx   = get_grp_mem_ndx( 'O' )
!         do_ox   = o3_ndx > 0 .and. o1d_ndx > 0 .and. o_ndx > 0
!      else
!         o3_ndx  = 1
!         o1d_ndx = 1
!         o_ndx   = 1
!         do_ox = .false.
!      end if

!      end subroutine chem_utls_inti

!      subroutine adjh2o( h2o, sh, mbar, vmr, plonl )
!!-----------------------------------------------------------------------
!!     ... transform water vapor from mass to volumetric mixing ratio
!!-----------------------------------------------------------------------

!      use gc_grid,  only : pcnstm1  !xiaolu,2017/06/01
!      use pmgrid,   only : plev    !add by zf 2008.07.02

!      implicit none

!!-----------------------------------------------------------------------
!!       ... dummy arguments
!!-----------------------------------------------------------------------
!      integer, intent(in) :: plonl
!      real, dimension(plonl,plev,pcnstm1), intent(in) :: &
!                               vmr                    ! xported species vmr
!      real, dimension(plonl,plev), intent(in) :: &
!                                sh                     ! specific humidity ( mmr )
!      real, dimension(plonl,plev), intent(in)  :: &
!                                mbar                   ! atmos mean mass
!      real, dimension(plonl,plev), intent(out) :: &
!                                 h2o                   ! water vapor vmr

!!-----------------------------------------------------------------------
!!       ... local variables
!!-----------------------------------------------------------------------
!      real, parameter :: mh2o = 1. /18.01528

!      integer ::   k, ndx_ch4
!      real    ::   t_value(plonl)

!-----------------------------------------------------------------------
!       ... limit dyn files water vapor
!-----------------------------------------------------------------------
!      ndx_ch4 = get_spc_ndx( 'CH4' )
!      if( ndx_ch4 > 0 ) then
!         do k = 1,plev
!            h2o(:,k)   = mbar(:,k) * sh(:plonl,k) * mh2o
!            t_value(:) = 6.e-6 - 2.*vmr(:,k,ndx_ch4)
!            where( t_value(:) > h2o(:,k) )
!               h2o(:,k) = t_value(:)
!            endwhere
!         end do
!      end if

!      end subroutine adjh2o

!      subroutine inti_mr_xform( sh, mbar, plonl )
!!-----------------------------------------------------------------
!!       ... initialize mean atmospheric "wet" mass
!!-----------------------------------------------------------------

!      use pmgrid,   only : plev    !add by zf 2008.07.02

!      implicit none

!!-----------------------------------------------------------------
!!       ... dummy args
!!-----------------------------------------------------------------
!      integer, intent(in) :: plonl
!      real, intent(in)  :: sh(plonl,plev)     ! specific humidity (kg/kg)
!      real, intent(out) :: mbar(plonl,plev)   ! mean wet atm mass ( amu )

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
!      real, parameter :: dry_mass = 28.966    ! amu
!      real, parameter :: mfac = 1. / .622

!      integer :: k

!      do k = 1,plev
!         mbar(:,k) = dry_mass
!      end do

!      end subroutine inti_mr_xform

!================
!FROM MOZART, NO NEED:XIAOLU,2017/06/07
!================
!      subroutine mmr2vmr( vmr, mmr, mbar, plonl )
!!-----------------------------------------------------------------
!!       ... xfrom from mass to volume mixing ratio
!!-----------------------------------------------------------------

!      use gc_chem_mods, only : adv_mass !xiaolu,2017/06/01
!      use gc_grid,      only : pcnstm1  !xiaolu,2017/06/01
!      use pmgrid,       only : plev

!      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
!      integer, intent(in) :: plonl
!      real, intent(in)  :: mbar(plonl,plev)
!      real, intent(in)  :: mmr(plonl,plev,pcnstm1)
!      real, intent(out) :: vmr(plonl,plev,pcnstm1)

!!-----------------------------------------------------------------
!!       ... local variables
!!-----------------------------------------------------------------
!      integer :: k, m

!      do m = 1,pcnstm1
!         if( adv_mass(m) /= 0. ) then
!            do k = 1,plev
!               vmr(:,k,m) = mbar(:,k) * mmr(:,k,m) / adv_mass(m)
!            end do
!         end if
!      end do

!      end subroutine mmr2vmr

!      subroutine vmr2mmr( vmr, mmr, nas, grp_ratios, mbar, plonl )
!!-----------------------------------------------------------------
!!       ... xfrom from mass to volume mixing ratio
!!-----------------------------------------------------------------

!      use gc_chem_mods, only : adv_mass, nadv_mass, grpcnt
!      use gc_grid,      only : pcnstm1
!      use pmgrid,       only : plev    !add by zf 2008.07.02

!      implicit none

!!-----------------------------------------------------------------
!!       ... dummy args
!!-----------------------------------------------------------------
!      integer, intent(in) :: plonl
!      real, intent(in)  :: mbar(plonl,plev)
!      real, intent(in)  :: vmr(plonl,plev,pcnstm1)
!      real, intent(out) :: mmr(plonl,plev,pcnstm1)
!      real, intent(in)  :: grp_ratios(:,:,:)
!      real, intent(out) :: nas(:,:,:)

!!-----------------------------------------------------------------
!!       ... local variables
!!-----------------------------------------------------------------
!      integer :: k, m
!      real    :: grp_mass(plonl)            ! weighted group mass

!!-----------------------------------------------------------------
!!       ... the non-group species
!!-----------------------------------------------------------------
!      do m = 1,pcnstm1
!         if( adv_mass(m) /= 0. ) then
!            do k = 1,plev
!               mmr(:,k,m) = adv_mass(m) * vmr(:,k,m) / mbar(:,k)
!            end do
!         end if
!      end do
!!-----------------------------------------------------------------
!!       ... the "group" species
!!-----------------------------------------------------------------
!      if( do_ox ) then
!         do k = 1,plev
!            grp_mass(:)     = grp_ratios(:,k,o3_ndx) * nadv_mass(o3_ndx) &
!                              + grp_ratios(:,k,o_ndx) * nadv_mass(o_ndx) &
!                              + grp_ratios(:,k,o1d_ndx) * nadv_mass(o1d_ndx)
!            mmr(:,k,ox_ndx)  = grp_mass(:) * vmr(:,k,ox_ndx) / mbar(:,k)
!            grp_mass(:)     = mmr(:,k,ox_ndx) / grp_mass(:)
!            nas(:,k,o3_ndx)  = nadv_mass(o3_ndx) * grp_ratios(:,k,o3_ndx) * grp_mass(:)
!            nas(:,k,o_ndx)   = nadv_mass(o_ndx) * grp_ratios(:,k,o_ndx) * grp_mass(:)
!            nas(:,k,o1d_ndx) = nadv_mass(o1d_ndx) * grp_ratios(:,k,o1d_ndx) * grp_mass(:)
!         end do
!      end if

!      end subroutine vmr2mmr
!==============================================================
!==============================================================
      integer function get_spc_ndx( spc_name )
!-----------------------------------------------------------------------
!     ... return overall species index associated with spc_name
!-----------------------------------------------------------------------

      use gc_chem_mods,  only : pcnstm1
      use gchp_tracname,   only : tracnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_spc_ndx = -1
      do m = 1,pcnstm1
         if( trim( spc_name ) == trim( tracnam(m) ) ) then
            get_spc_ndx = m
            exit
         end if
      end do

      end function get_spc_ndx

!comment XIAOLU,2017/06/07
!      integer function get_grp_ndx( grp_name )
!!-----------------------------------------------------------------------
!!     ... return group index associated with spc_name
!!-----------------------------------------------------------------------

!      use mo_chem_mods,  only : ngrp, grp_lst

!      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
!      character(len=*), intent(in) :: grp_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
!      integer :: m

!      get_grp_ndx = -1
!      do m = 1,ngrp
!         if( trim( grp_name ) == trim( grp_lst(m) ) ) then
!            get_grp_ndx = m
!            exit
!         end if
!      end do

!      end function get_grp_ndx

!      integer function get_grp_mem_ndx( mem_name )
!!-----------------------------------------------------------------------
!!     ... return group member index associated with spc_name
!!-----------------------------------------------------------------------

!      use mo_chem_mods,  only : grpcnt
!      use mo_tracname,    only : natsnam

!      implicit none

!!-----------------------------------------------------------------------
!!     ... dummy arguments
!!-----------------------------------------------------------------------
!      character(len=*), intent(in) :: mem_name

!!-----------------------------------------------------------------------
!!     ... local variables
!!-----------------------------------------------------------------------
!      integer :: m

!      get_grp_mem_ndx = -1
!      if( grpcnt > 0 ) then
!         do m = 1,max(1,grpcnt)
!            if( trim( mem_name ) == trim( natsnam(m) ) ) then
!               get_grp_mem_ndx = m
!               exit
!            end if
!         end do
!      end if

!      end function get_grp_mem_ndx

!      integer function get_het_ndx( het_name )
!!-----------------------------------------------------------------------
!!     ... return overall het process index associated with spc_name
!!-----------------------------------------------------------------------

!      use mo_chem_mods,  only : hetcnt, het_lst

!      implicit none

!!-----------------------------------------------------------------------
!!     ... dummy arguments
!!-----------------------------------------------------------------------
!      character(len=*), intent(in) :: het_name

!!-----------------------------------------------------------------------
!!     ... local variables
!!-----------------------------------------------------------------------
!      integer :: m

!      get_het_ndx = -1
!      do m = 1,max(1,hetcnt)
!         if( trim( het_name ) == trim( het_lst(m) ) ) then
!            get_het_ndx = m
!            exit
!         end if
!      end do

!      end function get_het_ndx

!      integer function get_extfrc_ndx( frc_name )
!!-----------------------------------------------------------------------
!!     ... return overall external frcing index associated with spc_name
!!-----------------------------------------------------------------------

!      use mo_chem_mods,  only : extcnt, extfrc_lst

!      implicit none

!!-----------------------------------------------------------------------
!!     ... dummy arguments
!!-----------------------------------------------------------------------
!      character(len=*), intent(in) :: frc_name

!!-----------------------------------------------------------------------
!!     ... local variables
!!-----------------------------------------------------------------------
!      integer :: m

!      get_extfrc_ndx = -1
!      if( extcnt > 0 ) then
!         do m = 1,max(1,extcnt)
!            if( trim( frc_name ) == trim( extfrc_lst(m) ) ) then
!               get_extfrc_ndx = m
!               exit
!            end if
!         end do
!      end if

!      end function get_extfrc_ndx

!      integer function get_rxt_ndx( rxt_alias )
!!-----------------------------------------------------------------------
!!     ... return overall external frcing index associated with spc_name
!!-----------------------------------------------------------------------

!     use mo_chem_mods,  only : rxt_alias_cnt, rxt_alias_lst, rxt_alias_map

!      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
!      character(len=*), intent(in) :: rxt_alias

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
!      integer :: m

!      get_rxt_ndx = -1
!      do m = 1,rxt_alias_cnt
!         if( trim( rxt_alias ) == trim( rxt_alias_lst(m) ) ) then
!            get_rxt_ndx = rxt_alias_map(m)
!            exit
!         end if
!      end do

!     end function get_rxt_ndx

      logical function has_drydep( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species dry deposition
!-----------------------------------------------------------------------

      use gc_chem_mods,  only : drydep_cnt, drydep_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_drydep = .false.
      do m = 1,drydep_cnt
         if( trim( spc_name ) == trim( drydep_lst(m) ) ) then
            has_drydep = .true.
            exit
         end if
      end do

      end function has_drydep

      logical function has_srfems( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species surface emission
!-----------------------------------------------------------------------

      use gc_chem_mods,     only : srfems_cnt, srfems_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name
!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_srfems = .false.
      do m = 1,srfems_cnt
         if( trim( spc_name ) == trim( srfems_lst(m) ) ) then
            has_srfems = .true.
            exit
         end if
      end do

      end function has_srfems

      end module gchp_chem_utls
