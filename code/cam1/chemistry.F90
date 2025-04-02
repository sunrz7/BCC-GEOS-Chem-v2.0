#include <misc.h>
#include <params.h>

module chemistry

!---------------------------------------------------------------------------------
! Module to parameterized greenhouse gas chemical loss frequencies from 
! Portmann and Solomon
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat, plev, plevp, plon, masterproc
  use ppgrid,       only: pcols, pver
  use physconst,    only: mwdry, mwch4, mwn2o, mwf11, mwf12, mwh2o, gravit
#ifdef MOZART2
  use constituents, only: ppcnst, cnst_add, cnst_name, advected   !zf 2010.03.15
#else
  use constituents, only: ppcnst, cnst_add, cnst_name, advected
#endif
#ifdef BCCCHEM
  use bcc_chem_mods, only : gas_pcnst
#endif
#ifdef WACCM
  use chem_mods,     only : gas_pcnst
#endif
#ifdef GEOSCHEM
  use TRACER_MOD,   only: TRACER_NAME
  use constituents, only: nonadvec, pnats, pcnst
#endif
#ifdef GHG_2D
  use ghg_surfvals_2d, only:  ghg_surfvals_get_ghg_2d
#endif
  use ghg_surfvals, only: ch4vmr, n2ovmr, f11vmr, f12vmr
  use abortutils,   only: endrun
  use timeinterp,   only: getfactors

  implicit none

  private
  save
!
! Public interfaces
!
#ifdef MOZART2
  public geoschem_register                ! register mozart2
  public gchp_init                    ! initialize (history) variables
  public geoschem_implements_cnst         ! returns true if consituent is implemented by this package  --zf 2008.06.26
#endif

#ifdef CO2
  public co2_register                ! register co2
#endif

#ifdef BCCCHEM
  public bcc_register
  public bccchem_init
  public mozart2_implements_cnst         ! returns true if consituent is implemented by this package  --zf 2008.06.26
#endif
#ifdef WACCM
  public waccm_chem_register
  public mozart2_implements_cnst         ! returns true if consituent is implemented by this package  --zf 2008.06.26
#endif

#ifdef GEOSCHEM
  public gc_register
  public geoschem_init
  public geoschem_implements_cnst
  public geoschem_init_cnst                        ! initialize mixing ratios if not read from initial file
#endif
  public :: ncnst
  public chem_register                             ! register consituents
  public chem_implements_cnst                      ! returns true if consituent is implemented by this package
  public chem_init_cnst                            ! initialize mixing ratios if not read from initial file
  public chem_init                                 ! initialize (history) variables
  public chem_timestep_init                        ! time interpolate chemical loss frequencies
  public chem_timestep_tend                        ! interface to tendency computation

! Namelist variable
  logical, public :: trace_gas=.false.             ! set true to activate this module

  integer, public :: ixchm

  integer, parameter :: ptrlon=01                  ! number of longitudes in input dataset
  integer, parameter :: ptrlat=36                  ! number of latitudes in input dataset
  integer, parameter :: ptrlev=56                  ! number of levels in input dataset
  integer, parameter :: ptrtim=12                  ! number of times(months) in input dataset

! Ratios of molecular weights
  real(r8), parameter :: rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
  real(r8), parameter :: rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
  real(r8), parameter :: rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
  real(r8), parameter :: rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
  real(r8), parameter :: rh2och4= mwh2o/mwch4      ! ratio of molecular weight h2o   to ch4

! Time/space interpolation of loss frequencies
  real(r8) :: tch4i  (plat,plev,ptrtim) ! input data  ch4   loss rate interp. in lat and lev
  real(r8) :: tn2oi  (plat,plev,ptrtim) ! input data  n2o   loss rate interp. in lat and lev
  real(r8) :: tcfc11i(plat,plev,ptrtim) ! input data  cfc11 loss rate interp. in lat and lev
  real(r8) :: tcfc12i(plat,plev,ptrtim) ! input data  cfc12 loss rate interp. in lat and lev
  real(r8) :: tch4m  (plat,plev,2)      ! input data  ch4   loss rate interp. in lat and lev
  real(r8) :: tn2om  (plat,plev,2)      ! input data  n2o   loss rate interp. in lat and lev
  real(r8) :: tcfc11m(plat,plev,2)      ! input data  cfc11 loss rate interp. in lat and lev
  real(r8) :: tcfc12m(plat,plev,2)      ! input data  cfc12 loss rate interp. in lat and lev
  real(r8) :: tch4   (plat,plev)        ! instantaneous ch4   loss rate 
  real(r8) :: tn2o   (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: tcfc11 (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: tcfc12 (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: cdaytrm                   ! calendar day for previous month data
  real(r8) :: cdaytrp                   ! calendar day for next month data

  integer :: np                               ! array index for previous month tracer data
  integer :: nm                               ! array index for next month tracer data
  integer :: np1                              ! current forward time index of tracer dataset
  integer :: date_tr(ptrtim)                  ! date on tracer dataset (YYYYMMDD)
  integer :: sec_tr(ptrtim)                   ! seconds of date on tracer dataset (0-86399)
  real(r8) :: ch4vmr_2d   (pcols,plev)
  real(r8) :: n2ovmr_2d   (pcols,plev)
  real(r8) :: f11vmr_2d   (pcols,plev)
  real(r8) :: f12vmr_2d   (pcols,plev)

! dummy values for specific heats at constant pressure
  real(r8), parameter:: cpch4 = 666.
  real(r8), parameter:: cpn2o = 666.
  real(r8), parameter:: cpf11 = 666.
  real(r8), parameter:: cpf12 = 666.

!-------zf 2008.2.02
#ifdef MOZART2
  integer, parameter :: ncnst=220
  character(len=8), dimension(ncnst) :: cnst_names  ! constituent names
   data  cnst_names/'ACET','ACTA','AERI','ALD2','ALK4','AONITA','AROMP4','AROMP5',& !8
                  'ATOOH','BALD','BCPI','BCPO','BENZ','BENZP','Br','Br2',&!8
                  'BrCl','BrNO2','BrNO3','BrO','BrSALA','BrSALC','BZCO3H','BZPAN',&!8
                  'C2H2','C2H4','C2H6','C3H8','CCl4','CFC11','CFC113','CFC114',&!8
                  'CFC115','CFC12','CH2Br2','CH2Cl2','CH2I2','CH2IBr','CH2ICl','CH2O',&!8
                  'CH3Br','CH3CCl3','CH3Cl','CH3I','CH4','CHBr3','CHCl3','Cl',&!8
                  'Cl2','Cl2O2','ClNO2','ClNO3','ClO','ClOO','CLOCK','CO',&!8
                  'CSL','DMS','DST1','DST2','DST3','DST4','EOH','ETHLN',&!8
                  'ETHN','ETHP','ETNO3','ETP', & !4
                  'GLYC','GLYX','H1211','H1301','H2402','H2O',& !6
                  'H2O2','HAC','HBr','HC5A','HCFC123','HCFC141b','HCFC142b','HCFC22',&!75-82
                  'HCl','HCOOH','HI','HMHP','HMML','HMS','HNO2','HNO3', &!83-90
                  'HNO4','HOBr','HOCl','HOI','HONIT','HPALD1','HPALD2','HPALD3', &!91-98
                  'HPALD4','HPETHNL','I','I2','I2O2','I2O3','I2O4','IBr', &!99-106
                  'ICHE','ICl','ICN','ICPDH','IDC','IDCHP','IDHDP','IDHPE',&!107-114
                  'IDN','IEPOXA','IEPOXB','IEPOXD','IHN1','IHN2','IHN3','IHN4',&!115-122
                  'INDIOL','INO','INPB','INPD','IO','IONITA','IONO','IONO2',&!123-130
                  'IPRNO3','ISALA','ISALC','ISOP','ITCN','ITHN','LIMO','LVOC',&!131-138
                  'LVOCOA','MACR','MACR1OOH','MAP','MCRDH','MCRENOL','MCRHN','MCRHNB',&!139-146
                  'MCRHP','MCT','MEK','MENO3','MGLY','MOH','MONITA','MONITS',&!147-154
                  'MONITU','MP','MPAN','MPN','MSA','MTPA','MTPO','MVK',& !155-162
                  'MVKDH','MVKHC','MVKHCB','MVKHP','MVKN','MVKPC','N2O','N2O5',& !163-170
                  'NH3','NH4','NIT','NITs','NO','NO2','NO3','NPHEN',& !171-178
                  'NPRNO3','O3','OClO','OCPI','OCPO','OCS','OIO','PAN',&!179-186
                  'pFe','PHEN','PIP','PP','PPN','PROPNN','PRPE','PRPN',& !187-194     
                  'PYAC','R4N2','R4P','RA3P','RB3P','RCHO','RIPA','RIPB',&!195-202
                  'RIPC','RIPD','RP','SALA','SALAAL','SALACL','SALC','SALCAL',&!203-210
                  'SALCCL','SO2','SO4','SO4s','SOAGX','SOAIE','SOAP','SOAS',& !211-218
                  'TOLU','XYLE'/   !220
  real(r8), parameter :: cpche(ncnst) = 666.
  real(r8), parameter :: mwche(ncnst) = 9999.
!------------------------

#elif (defined CO2)

  integer, parameter :: ncnst=1
  character(len=8), dimension(ncnst) :: cnst_names  ! constituent names
  data cnst_names /'CO2     '/ 
  real(r8), parameter :: cpche(ncnst) = 666.
  real(r8), parameter :: mwche(ncnst) = 9999.

#elif (defined BCCCHEM)
  integer, parameter :: ncnst= 75                    ! 59 + 16

  character(len=8), dimension(ncnst) :: cnst_names ! constituent names, be set by bcc_register

!  character(len=8), dimension(ncnst), parameter :: & ! constituent names
!     cnst_names = (/ 'O3      ','O       ','O1D     ','O2      ','O2_1S   ', &
!                     'O2_1D   ','N2O     ','N       ','NO      ','NO2     ', &
!                     'NO3     ','HNO3    ','HO2NO2  ','N2O5    ','CH4     ', &
!                     'CH3O2   ','CH3OOH  ','CH2O    ','CO      ','H2      ', &
!                     'H       ','OH      ','HO2     ','H2O2    ','CLY     ', &
!                     'BRY     ','CL      ','CL2     ','CLO     ','OCLO    ', &
!                     'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','BRCL    ', &
!                     'BR      ','BRO     ','HBR     ','HOBR    ','BRONO2  ', &
!                     'CH3CL   ','CH3BR   ','CFC11   ','CFC12   ','CFC113  ', &
!                     'HCFC22  ','CCL4    ','CH3CCL3 ','CF3BR   ','CF2CLBR ', &
!                     'CO2     ','N2p     ','O2p     ','Np      ','Op      ', &
!                     'NOp     ','e       ','N2D     ','H2O     ',            &
!                     'SO2     ','SO4     ','DMS     ','OC1     ','OC2     ', &
!                     'CB1     ','CB2     ','SSLT01  ','SSLT02  ','SSLT03  ', &
!                     'SSLT04  ','DST01   ','DST02   ','DST03   ','DST04   ' /)
!

  real(r8), parameter :: cpche(ncnst) = 666.
  real(r8), parameter :: mwche(ncnst) = 9999.

#elif (defined WACCM)
  integer, parameter :: ncnst= 57

  character(len=8), dimension(ncnst) :: cnst_names ! constituent names, be set by bcc_register

!  character(len=8), dimension(ncnst), parameter :: & ! constituent names
!     cnst_names = (/ 'O3      ','O       ','O1D     ','O2      ','O2_1S   ', &
!                     'O2_1D   ','N2O     ','N       ','NO      ','NO2     ', &
!                     'NO3     ','HNO3    ','HO2NO2  ','N2O5    ','CH4     ', &
!                     'CH3O2   ','CH3OOH  ','CH2O    ','CO      ','H2      ', &
!                     'H       ','OH      ','HO2     ','H2O2    ',            &
!                     'CL      ','CL2     ','CLO     ','OCLO    ', &
!                     'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','BRCL    ', &
!                     'BR      ','BRO     ','HBR     ','HOBR    ','BRONO2  ', &
!                     'CH3CL   ','CH3BR   ','CFC11   ','CFC12   ','CFC113  ', &
!                     'HCFC22  ','CCL4    ','CH3CCL3 ','CF3BR   ','CF2CLBR ', &
!                     'CO2     ','N2p     ','O2p     ','Np      ','Op      ', &
!                     'NOp     ','e       ','N2D     ','H2O     '

  real(r8), parameter :: cpche(ncnst) = 666.
  real(r8), parameter :: mwche(ncnst) = 9999.

#else
!----
#ifdef GEOSCHEM
#ifdef TAGGED_CO2
  integer, parameter :: ncnst=54
#else

  integer, parameter :: ncnst=4 ! TEMPORARY - MSLSept09,2012
!
!  integer, parameter :: ncnst=45 ! GEOS-Chem with FULLCHEM - MSLSept17,2012
!                                 ! including CH4 as an advected tracer in
!                                 ! GC's input.geos___.rc file. OH, CO2 not advected
#endif
  character(len=8), dimension(ncnst) :: & ! constituent names
                cnst_names != TRACER_NAME!(/'CO   '/)

  real(r8), parameter :: cpche(ncnst) = 666.
  real(r8), parameter :: mwche(ncnst) = 9999.
!------
#else
  integer, parameter :: ncnst=4                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)
#endif
!----------
#endif
!----------------------------------------
  character(len=8) :: srcnam(ncnst)                  ! names of source/sink tendencies
  integer :: ixghg                                   ! index of 1st constituent (N2O)


contains

!------ zf 2008.2.02
#ifdef MOZART2
  subroutine geoschem_register
!-----------------------------------------------------------------------
!
! Purpose: register advected constituents for mozart2
!
! Author: ZF
!-----------------------------------------------------------------------
    use gc_chem_mods,     only : adv_mass, nadv_mass
!    use bcc_chem_mods,    only : bcc_chem_mods_inti, imozart

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    integer :: i
    real(r8), parameter :: zero  = 0._r8
real(r8)            :: mw(220) !sunrz 2024/8

     mw=(/58.09,60.06,126.90,44.06,58.12,189.12,68.08,98.10,              &
                  90.09,106.12,12.01,12.01,78.12,110.11,79.90,159.80,     &
                  115.45,125.91,141.91,95.90,79.90,79.90,138.12,183.12,   &
                  26.05,28.05,30.08,44.11,153.82,137.37,187.38,170.92,    &
                  154.47,120.91,173.83,84.93,267.84,220.84,176.38,30.03,  &
                  94.94,133.35,50.45,141.94,16.0,252.73,119.35,35.45,     &
                  70.90,102.91,81.45,97.45,51.45,67.45,1.0,28.0,          &
                  108.14,62.13,29.0,29.0,29.0,29.0,46.07,105.06,          &
                  107.07,78.07,91.08,62.08,                               &
                  60.06,58.04,165.36,148.91,259.82,18.02,                 &
                  34.02,74.08,80.91,100.13,152.93,116.94,100.50,86.47,    &
                  36.45,46.03,127.91,64.05,102.10,111.10,47.01,63.01,     &
                  79.01,96.91,52.45,143.89,215.0,116.13,116.13,116.13,    &
                  116.13,76.06,126.90,253.80,285.80,301.80,317.80,206.90, &
                  116.13,162.45,145.13,150.15,98.11,148.13,168.17,150.15, &
                  192.15,106.14,106.14,106.14,147.15,147.15,147.15,147.15,&
                  102.0,156.91,163.15,163.15,142.90,14.01,172.91,188.91,  &
                  105.11,126.90,126.90,68.13,195.15,197.17,136.26,154.19, &
                  154.19,70.10,102.10,76.06,104.12,86.10,149.11,149.11,   &
                  120.12,124.0,72.11,77.05,72.07,32.05,14.01,215.28,      &
                  215.28,48.05,147.10,93.05,96.10,136.26,136.26,70.09,    &
                  105.13,102.10,102.10,120.12,149.12,118.10,44.02,108.02, &
                  17.04,18.05,62.01,31.4,30.01,46.01,62.01,139.11,        &
                  105.11,48.00,67.45,12.01,12.01,60.07,158.90,121.06,     &
                  55.85,94.11,186.28,92.11,135.08,119.08,42.09,137.11,    &
                  88.07,119.10,90.14,76.11,76.11,58.09,118.15,118.15,     &
                  118.15,118.15,90.09,31.4,31.4,35.45,31.4,31.4,          &
                  35.45,64.04,96.06,31.4,58.04,118.15,150.0,150.0,        &
                  92.15,106.18/)
!-----------------------------------------------------------------------
! Set names of diffused variable tendencies and declare them as history variables

    do i = 1, ncnst
!wtw    call cnst_add(cnst_names(i),advected,mwche(i),cpche(i),zero,m,mixtype='dry')
    call cnst_add(cnst_names(i),advected,mw(i),cpche(i),zero,m)
    if (i==1) ixchm = m
    enddo


!--------------------------
!    call bcc_chem_mods_inti

  end subroutine geoschem_register
#endif

#ifdef CO2
  subroutine co2_register
!-----------------------------------------------------------------------
!
! Purpose: register advected constituents for co2
!
! Author: ZF
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    integer :: i
    real(r8), parameter :: zero  = 0._r8
    real(r8) :: adv_mass
!-----------------------------------------------------------------------
! Set names of diffused variable tendencies and declare them as history variables

    adv_mass = 44.0    !only include CO2

    do i = 1, ncnst
      call cnst_add(cnst_names(i),advected,adv_mass, cpche(i),zero,m)
      if (i==1) ixchm = m
    enddo

  end subroutine co2_register
#endif


#ifdef BCCCHEM

  subroutine bcc_register
!-----------------------------------------------------------------------
!
! Purpose: register advected constituents and physics buffer fields
!
!-----------------------------------------------------------------------
    use constituents,        only : cnst_fixed_ubc, ppcnst
    use bcc_chem_mods,       only : gas_pcnst
    use bcc_chem_mods,       only : set_sim_dat
    use bcc_chem_mods,       only : solsym, adv_mass
    use bcc_chem_mods,       only : map2chm, imozart
    use bcc_exbdrift,        only : exbdrift_register
    use bcc_cfc11star,       only : register_cfc11star
    use bcc_chem_mods,       only : drydep_cnt, drydep_list, drydep_mapping
    use bcc_chem_mods,       only : gas_wetdep_cnt, gas_wetdep_list, gas_wetdep_mapping
    use bcc_chem_mods,       only : aer_wetdep_cnt, aer_wetdep_list, aer_wetdep_mapping
!
!wtw    use short_lived_species, only : slvd_index, short_lived_map=>map, register_short_lived_species
!
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
    integer :: m, n                            ! tracer index
    integer :: i, L
    real(r8), parameter :: zero  = 0._r8
        real(r8)            :: qminc

!-----------------------------------------------------------------------
! Set the simulation chemistry variables
!-----------------------------------------------------------------------

    if( gas_pcnst .ne. ncnst) then
            write(*,*) 'STOP in bcc_register: ncnst =\= gas_pcnst'
        call endrun
    endif

    cnst_names(1:ncnst) = solsym(1:gas_pcnst)

    call set_sim_dat
!----
    qminc = 1.e-36_r8
    do i = 1, ncnst
       if ( trim(cnst_names(i)) == 'O3' ) then
          qminc          = 1.e-12_r8
       else if ( trim(cnst_names(i)) == 'CH4' ) then
          qminc          = 1.e-12_r8
       else if ( trim(cnst_names(i)) == 'N2O' .or. trim(cnst_names(i)) == 'CO2' ) then
          qminc = 1.e-15_r8
       else if ( trim(cnst_names(i)) == 'CFC11' .or. trim(cnst_names(i)) == 'CFC12' ) then
          qminc = 1.e-20_r8
       else if ( trim(cnst_names(i)) == 'O2_1S' .or. trim(cnst_names(i)) == 'O2_1D' ) then
          qminc = 1.e-20_r8
       endif
!--
       call cnst_add(cnst_names(i),advected,adv_mass(i),cpche(i),qminc,m,mixtype='wet')

       if( i == 1) then
           ixchm   = m
           imozart = m
       end if
!----
!wtw
       map2chm(m) = i
!---
    enddo

! wtw
!
    drydep_mapping(:) = 0
    do i=1, drydep_cnt
       do L = 1, gas_pcnst
          if(  trim( drydep_list(i) ) == trim( solsym(L) ) ) then
             drydep_mapping(i)  = L
             exit
          end if
       end do
        enddo

    gas_wetdep_mapping(:) = 0
    do i=1, gas_wetdep_cnt
       do L = 1, gas_pcnst
          if(  trim( gas_wetdep_list(i) ) == trim( solsym(L) ) ) then
             gas_wetdep_mapping(i)  = L
             exit
          end if
       end do
    enddo

    aer_wetdep_mapping(:) = 0
    do i=1, aer_wetdep_cnt
       do L = 1, gas_pcnst
          if(  trim( aer_wetdep_list(i) ) == trim( solsym(L) ) ) then
             aer_wetdep_mapping(i)  = L
             exit
          end if
       end do
    enddo

!----------------
    ! BAB: 2004-09-01 kludge to define a fixed ubc for water vapor
    !      required because water vapor is not declared by chemistry but
    !      has a fixed ubc only if chemistry is running.
    !-----------------------------------------------------------------------

    cnst_fixed_ubc(:) = .false.

    do n = 1,ppcnst
       m = map2chm(n)
       if( trim(solsym(m)) == 'CH4' .or.   trim(solsym(m)) ==  'CO2'  .or. &
           trim(solsym(m)) == 'O2'  .or.   trim(solsym(m)) ==  'N'    .or. &
           trim(solsym(m)) == 'NO'  .or.   trim(solsym(m)) ==  'CO'   .or. &
           trim(solsym(m)) == 'H'   .or.   trim(solsym(m)) ==  'H2'   .or. &
           trim(solsym(m)) == 'O'   .or.   n == 1 ) then

           cnst_fixed_ubc(n)   = .true.
       endif
    enddo

    call register_cfc11star()

    !
    ! Initialize e and b fields
        !
    call exbdrift_register()

        return
    end subroutine bcc_register

#endif
!==============================================================================

#ifdef WACCM
subroutine waccm_chem_register
!-----------------------------------------------------------------------
!
! Purpose: register advected constituents and physics buffer fields
!
!-----------------------------------------------------------------------

    use m_spc_id
    use chem_mods,       only : gas_pcnst, adv_mass
    use mo_tracname,     only : solsym
    use mo_sim_dat,      only : set_sim_dat
    use chem_mods,       only : imozart, map2chm
    use constituents,    only : cnst_fixed_ubc

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: m, n, i                            ! tracer index
    real(r8) :: qmin                                ! min value
    logical  :: ic_from_cam2                        ! wrk variable for initial cond input
    logical  :: has_fixed_ubc                       ! wrk variable for upper bndy cond
    character(len=128) :: lng_name                  ! variable long name

    real(r8) :: qminc
!-----------------------------------------------------------------------
! Set the simulation chemistry variables
!-----------------------------------------------------------------------

    call set_sim_dat

    if( gas_pcnst .ne. ncnst) then
            write(*,*) 'STOP in bcc_register: ncnst =\= gas_pcnst'
            call endrun
    endif

    cnst_names(1:ncnst) = solsym(1:gas_pcnst)
!-----------------------------------------------------------------------
! BAB: 2004-09-01 kludge to define a fixed ubc for water vapor
!      required because water vapor is not declared by chemistry but
!      has a fixed ubc only if chemistry is running.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Set names of diffused variable tendencies and declare them as history variables
!-----------------------------------------------------------------------
    qminc = 1.e-36_r8
    do i = 1, ncnst

       if ( trim(cnst_names(i)) == 'O3' ) then
          qminc          = 1.e-12_r8
       else if ( trim(cnst_names(i)) == 'CH4' ) then
          qminc          = 1.e-12_r8
       else if ( trim(cnst_names(i)) == 'N2O' .or. trim(cnst_names(i)) == 'CO2' ) then
          qminc = 1.e-15_r8
       else if ( trim(cnst_names(i)) == 'CFC11' .or. trim(cnst_names(i)) == 'CFC12' ) then
          qminc = 1.e-20_r8  
       else if ( trim(cnst_names(i)) == 'O2_1S' .or. trim(cnst_names(i)) == 'O2_1D' ) then
          qminc = 1.e-20_r8  
       endif
!-- 
!      call cnst_add(cnst_names(i),advected,mwche(i),cpche(i),qminc,m,mixtype='wet')
       call cnst_add(cnst_names(i),advected,adv_mass(i),cpche(i),qminc,m,mixtype='wet')
    
       if( i == 1) then      
           ixchm   = m
           imozart = m
       end if
       map2chm(m) = i
    end do
    write(*,*) 'imozart = ',imozart
    write(*,*) 'map2chm = ',map2chm(:)

    cnst_fixed_ubc(:) = .false.

    do n = 1,ppcnst
       m = map2chm(n)
!       if( trim(solsym(m)) .eq. 'CH4' .or.   trim(solsym(m)) .eq. 'CO2'  .or. &
!           trim(solsym(m)) .eq. 'O2'  .or.   trim(solsym(m)) .eq. 'N'    .or. &
!           trim(solsym(m)) .eq. 'NO'  .or.   trim(solsym(m)) .eq. 'CO'   .or. &
!           trim(solsym(m)) .eq. 'H'   .or.   trim(solsym(m)) .eq. 'H2'   .or. &
!           trim(solsym(m)) .eq. 'O'   .or.   n == 1 ) then

       if( trim(solsym(m)) .eq. 'CO2'  .or. &
           trim(solsym(m)) .eq. 'O2'  .or.   trim(solsym(m)) .eq. 'N'    .or. &
           trim(solsym(m)) .eq. 'NO'  .or.   trim(solsym(m)) .eq. 'CO'   .or. &
           trim(solsym(m)) .eq. 'H'   .or.   trim(solsym(m)) .eq. 'H2'   .or. &
           trim(solsym(m)) .eq. 'O'   .or.   n == 1 ) then

           cnst_fixed_ubc(n)   = .true.
       endif
    enddo

  end subroutine waccm_chem_register

#endif
!===============================================================================
#ifdef GEOSCHEM
  subroutine gc_register
!-----------------------------------------------------------------------
!
! Purpose: register advected constituents for GEOSCHEM
!
! Author: Michael
!-----------------------------------------------------------------------

    use ppgrid,             only : pcols, pver
    use gc_initrun_mod,     only : gc_getopts, gc_init_dimensions
    use history,            only: addfld, add_default, phys_decomp
!    use gc_environment_mod, only : tracer_index

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    integer :: i
    real(r8), parameter :: zero  = 0._r8
!-----------------------------------------------------------------------
! Set names of diffused variable tendencies and declare them as history variables

    call gc_init_dimensions(pcols, pver) ! Allocate GEOS-Chem arrays
    ! Is there a better location for this? It needs to be done through
    ! 'call initindx'... prior to 'call inti' so that GEOS-Chem
    ! index arrays are allocated

!    if ( .not. ALLOCATED( TRACER_INDEX ) ) ALLOCATE( TRACER_INDEX(ncnst) )
    call gc_getopts ! Read "input.geos"
                    ! Currently executed for all nodes. It would be more
                    ! efficient to read @ masterproc and scatter just the

#ifdef TAGGED_CO2
       cnst_names = (/ 'CO2     ','CO2ff   ','CO2oc   ','CO2bal  ', &
            'CO2bb   ','CO2bf   ','CO2nte  ','CO2se   ','CO2av   ', &
            'CO2ch   ','CO2corr ','CO2bg   ','CO2xx   ','CO2_1   ', &
            'CO2_2   ','CO2_3   ','CO2_4   ','CO2_5   ','CO2_6   ', &
            'CO2_7   ','CO2_8   ','CO2_9   ','CO2_10  ','CO2_11  ', &
            'CO2_12  ','CO2_13  ','CO2_14  ','CO2_15  ','CO2_16  ', &
            'CO2_17  ','CO2_18  ','CO2_19  ','CO2_20  ','CO2_21  ', &
            'CO2_22  ','CO2_23  ','CO2_24  ','CO2_25  ','CO2_26  ', &
            'CO2_27  ','CO2_28  ','CO2_29  ','CO2_30  ','CO2_31  ', &
            'CO2_32  ','CO2_33  ','CO2_34  ','CO2_35  ','CO2_36  ', &
            'CO2_37  ','CO2_38  ','CO2_39  ','CO2se2  ','CO2av2  '/)
#else
       write(*,*) 'TRACER NAMES: ', TRACER_NAME, shape(TRACER_NAME)
       cnst_names = TRACER_NAME
#endif
    do i = 1, ncnst
!wtw      call cnst_add(cnst_names(i),advected,mwche(i),cpche(i),zero,m,mixtype='dry')
       call cnst_add(cnst_names(i),advected,adv_mass(i),cpche(i),zero,m)

!       call addfld (cnst_names(i),'kg/kg ',pver, 'A',trim(cnst_name(i)),phys_decomp)
!       call add_default (cnst_name(i), 1, ' ')

       if ( masterproc) then
          write(*,*) 'Registered ', trim(cnst_names(i))
       endif
       if (i==1) ixchm = m
    enddo
    i = ncnst+pnats
!wtw   call cnst_add('OH', nonadvec,17.01,cpche(ncnst),zero,m,mixtype='dry')
    call cnst_add('OH', nonadvec,17.01,cpche(ncnst),zero,m)
!wtw   call cnst_add('HO2',nonadvec,33.01,cpche(ncnst),zero,m,mixtype='dry')
    call cnst_add('HO2',nonadvec,33.01,cpche(ncnst),zero,m)

  end subroutine gc_register
#endif

!===============================================================================
  subroutine chem_register
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents for parameterized greenhouse gas chemistry
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------

! Set names of diffused variable tendencies and declare them as history variables

    call cnst_add(cnst_names(1), advected, mwn2o, cpn2o, 0._r8, ixghg, longname='Nitrous Oxide')
    call cnst_add(cnst_names(2), advected, mwch4, cpch4, 0._r8, m, longname='Methane')
    call cnst_add(cnst_names(3), advected, mwf11, cpf11, 0._r8, m)
    call cnst_add(cnst_names(4), advected, mwf12, cpf12, 0._r8, m)

    return
  end subroutine chem_register

!===============================================================================

  function chem_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: chem_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     chem_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           chem_implements_cnst = .true.
           return
        end if
     end do
  end function chem_implements_cnst

!===============================================================================

#if (defined MOZART2 ) || (defined BCCCHEM) || (defined WACCM)

  function geoschem_implements_cnst(name)
!-----------------------------------------------------------------------
!
! Purpose: return true if specified constituent is implemented by this package
!
! Author: ZF
!
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: geoschem_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     geoschem_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           geoschem_implements_cnst = .true.
           return
        end if
     end do
  end function geoschem_implements_cnst

#endif

#ifdef GEOSCHEM

  function geoschem_implements_cnst(name)
!-----------------------------------------------------------------------
!
! Purpose: return true if specified constituent is implemented by this package
!
! Author:
!
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: geoschem_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     geoschem_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           geoschem_implements_cnst = .true.
           return
        end if
     end do
  end function geoschem_implements_cnst

#endif

  subroutine chem_init
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterized greenhouse gas chemistry
!          (declare history variables)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    use history,    only: addfld, add_default, phys_decomp

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------

! Set names of diffused variable tendencies and declare them as history variables
    do m = 1, ncnst
       srcnam(m) = trim(cnst_name(ixghg-1+m)) // 'SRC'
       call addfld (srcnam(m),'kg/kg/s ',pver, 'A',trim(cnst_name(ixghg-1+m))//' source/sink',phys_decomp)
       call add_default (srcnam(m), 1, ' ')
    end do

    call chem_init_loss

    return
  end subroutine chem_init

!===============================================================================

#ifdef MOZART2
  subroutine gchp_init

! Purpose: initialize parameterized mozart2
!          (declare history variables)
!-----------------------------------------------------------------------
    use history,     only: addfld, add_default, phys_decomp
    use gchp_inirun,   only: gc_inirun

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------

! Set names of diffused variable tendencies and declare them as history variables
!    do m = 1, ncnst
!       srcnam(m) = trim(cnst_name(ixchm-1+m)) // 'SRC'
!       call addfld (srcnam(m),'kg/kg/s ',pver, 'A',trim(cnst_name(ixchm-1+m))//' source/sink',phys_decomp)
!       call add_default (srcnam(m), 1, ' ')
!    end do

    call gc_inirun    ! zf 2008.06.24

    return
  end subroutine gchp_init
#endif

#ifdef GEOSCHEM
      subroutine geoschem_init

        use gc_initrun_mod,   only: gc_initrun!, gc_getopts
        use gc_interface_mod, only: gc_met, gc_state

        implicit none

!        call gc_getopts
        call gc_initrun(gc_met, gc_state, pcols, pver)
      end subroutine geoschem_init
#endif


#ifdef BCCCHEM
  subroutine bccchem_init
!-----------------------------------------------------------------------
!
! Purpose: initialize parameterized greenhouse gas chemistry
!          (declare history variables)
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: NCAR CMS
!
!-----------------------------------------------------------------------
    use history,           only: addfld, add_default, phys_decomp
    use bcc_inirun,        only: inirun
        use bcc_cfc11star,     only: init_cfc11star

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------
!
!wtw    call bcc_constants_inti
!
!-------------------------------------------

    call inirun

!--------------------------------------------------
!wtw     ghg_chem = .false.
!wtw     if ( ghg_chem ) then
!wtw        call ghg_chem_init(phys_state, bndtvg, h2orates)
!wtw     endif

!wtw     call lin_strat_chem_inti(phys_state)
!wtw     call chlorine_loading_init( chlorine_loading_file, &
!wtw                                 type = chlorine_loading_type, &
!wtw                                 ymd = chlorine_loading_fixed_ymd, &
!wtw                                 tod = chlorine_loading_fixed_tod )

     call init_cfc11star()

!----------------------------------

    return
  end subroutine bccchem_init
#endif

!======================================================================
  subroutine chem_timestep_tend (state, ptend, cflx, dt, fh2o)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_lat_all_p
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

    type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
    real(r8), intent(inout) :: cflx(pcols,ppcnst)      ! Surface constituent flux (kg/m^2/s)
    real(r8), intent(out) :: fh2o(pcols)          ! H2O flux required to balance the H2O source from methane chemistry (kg/m^2/s)
!
! Local variables
!
    integer :: i, k, m                            ! indices
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: lat(pcols)                         ! latitude index for S->N storage
!
!-----------------------------------------------------------------------
! Initialize output tendency structure
    call physics_ptend_init(ptend)

    ioff  = ixghg - 1
    lchnk = state%lchnk
    ncol  = state%ncol

! get latitude indices
    call get_lat_all_p(lchnk, ncol, lat)

! compute tendencies and surface fluxes
    call ghg_chem ( lchnk, ncol, lat,                                               &
         state%q(:,:,1), state%q(:,:,ixghg:ixghg+3),                              &
         ptend%q(:,:,1), ptend%q(:,:,ixghg:ixghg+3), cflx(:,ixghg:ixghg+3), dt)

! set flags for tracer tendencies (water and 4 ghg's)
    ptend%lq(1)             = .TRUE.
    ptend%lq(ioff+1:ioff+4) = .TRUE.
!
! record tendencies on history files
    do m = 1, 4
       call outfld (srcnam(m),ptend%q(:,:,ioff+m),pcols,lchnk)
    end do

! Compute water vapor flux required to make the conservation checker happy
    fh2o = 0
    do k = 1, pver
       do i = 1, ncol
          fh2o(i) = fh2o(i) + ptend%q(i,k,1)*state%pdel(i,k)/gravit
       end do
    end do

    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine ghg_chem (lchnk, ncol, lat, qh2o, qghg, dqh2o, dqghg, fghg, dt)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Apply the interpolated chemical loss rates from the input data to
! N2O, CH4, CFC11 and CFC12. Set the surface values to a constant.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: lat(pcols)             ! latitude index for S->N storage

    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: qh2o(pcols,pver)      ! mass mixing ratios of water vapor
    real(r8), intent(in) :: qghg(pcols,pver,4)    ! mass mixing ratios of greenhouse gases

    real(r8), intent(out) :: dqh2o(pcols,pver)    ! tendency of mass mixing ratios (water)
    real(r8), intent(out) :: dqghg(pcols,pver,4)  ! tendency of mass mixing ratios (ghg's)
    real(r8), intent(out) :: fghg(pcols,4)        ! Surface constituent flux (kg/m^2/s)
!
! Local variables
!
    integer i,k                                   ! loop indexes
    real(r8) xch4                                 ! new methane mass mixing ratio
    real(r8) xn2o                                 ! new nitrous oxide mass mixing ratio
    real(r8) xcfc11                               ! new cfc11 mass mixing ratio
    real(r8) xcfc12                               ! new cfc12 mass mixing ratio
!
!-----------------------------------------------------------------------
!
! Apply chemical rate coefficient using time split implicit method. The
! turn the new value back into a tendency. NOTE that water tendency is
! twice methane tendency. Water is specific humidity which is in mass
! mixing ratio units. Note that
!  o 1   => indx of n2o
!  o 2 => indx of ch4
!  o 3 => indx of cfc11
!  o 4 => indx of cfc12
!
    do k=1,pver-2
       do i=1,ncol
          xn2o         = qghg(i,k,1) / (1. + tn2o  (lat(i),k) * dt)
          xch4         = qghg(i,k,2) / (1. + tch4  (lat(i),k) * dt)
          xcfc11       = qghg(i,k,3) / (1. + tcfc11(lat(i),k) * dt)
          xcfc12       = qghg(i,k,4) / (1. + tcfc12(lat(i),k) * dt)

          dqghg(i,k,1) =(xn2o   - qghg(i,k,1)) / dt
          dqghg(i,k,2) =(xch4   - qghg(i,k,2)) / dt
          dqghg(i,k,3) =(xcfc11 - qghg(i,k,3)) / dt
          dqghg(i,k,4) =(xcfc12 - qghg(i,k,4)) / dt

          dqh2o(i,k)   = -2. * rh2och4 * dqghg(i,k,2)
       end do
    end do
!
! Set the "surface" tendencies (bottom 2 levels) to maintain specified
! tropospheric concentrations.
!
#ifdef GHG_2D
    call ghg_surfvals_get_ghg_2d(lchnk,ch4vmr_2d,n2ovmr_2d,f11vmr_2d,f12vmr_2d)
    do k = pver-1, pver
       do i=1,ncol
          dqghg(i,k,1) =((rmwn2o*n2ovmr_2d(i,k)) - qghg(i,k,1)) / dt
          dqghg(i,k,2) =((rmwch4*ch4vmr_2d(i,k)) - qghg(i,k,2)) / dt
          dqghg(i,k,3) =((rmwf11*f11vmr_2d(i,k)) - qghg(i,k,3)) / dt
          dqghg(i,k,4) =((rmwf12*f12vmr_2d(i,k)) - qghg(i,k,4)) / dt
          dqh2o(i,k)   = 0.
       end do
    end do
#else
    do k = pver-1, pver
       do i=1,ncol
          dqghg(i,k,1) =((rmwn2o*n2ovmr) - qghg(i,k,1)) / dt
          dqghg(i,k,2) =((rmwch4*ch4vmr) - qghg(i,k,2)) / dt
          dqghg(i,k,3) =((rmwf11*f11vmr) - qghg(i,k,3)) / dt
          dqghg(i,k,4) =((rmwf12*f12vmr) - qghg(i,k,4)) / dt
          dqh2o(i,k)   = 0.
       end do
    end do
#endif
!
! For now set all tracer fluxes to 0
!
    do i=1,ncol
       fghg(i,1) = 0.
       fghg(i,2) = 0.
       fghg(i,3) = 0.
       fghg(i,4) = 0.
    end do

    return
  end subroutine ghg_chem

!===============================================================================
  subroutine chem_init_loss
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Do initial read of time-variant chemical loss rate frequency dataset, containing
! loss rates as a function of latitude and pressure.  Determine the two
! consecutive months between which the current date lies.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
!-----------------------------------------------------------------------
    use ioFileMod
    use commap
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use filenames, only: bndtvg
#if ( defined SPMD )
    use mpishorthand
#endif

    implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>
    include 'netcdf.inc'

!
! Local variables
!
    integer dateid            ! netcdf id for date variable
    integer secid             ! netcdf id for seconds variable
    integer lonid             ! netcdf id for longitude variable
    integer latid             ! netcdf id for latitude variable
    integer levid             ! netcdf id for level variable
    integer timid             ! netcdf id for time variable
    integer tch4id            ! netcdf id for ch4   loss rate
    integer tn2oid            ! netcdf id for n2o   loss rate
    integer tcfc11id          ! netcdf id for cfc11 loss rate
    integer tcfc12id          ! netcdf id for cfc12 loss rate
    integer lonsiz            ! size of longitude dimension on tracer dataset
    integer levsiz            ! size of level dimension on tracer dataset
    integer latsiz            ! size of latitude dimension on tracer dataset
    integer timsiz            ! size of time dimension on tracer dataset
    integer j,n,k,nt          ! indices
    integer ki,ko,ji,jo       ! indices
    integer :: yr, mon, day   ! components of a date
    integer :: ncdate         ! current date in integer format [yyyymmdd]
    integer :: ncsec          ! current time of day [seconds]

    real(r8) :: calday        ! current calendar day
    real(r8) lato(plat)       ! cam model latitudes (degrees)
    real(r8) zo(plev)         ! cam model heights (m)
    real(r8) lati(ptrlat)     ! input data latitudes (degrees)
    real(r8) pin(ptrlev)      ! input data pressure values (mbars)
    real(r8) zi(ptrlev)       ! input data heights (m)

    real(r8) xch4i  (ptrlon,ptrlev,ptrlat,ptrtim) ! input ch4   loss rate coeff
    real(r8) xn2oi  (ptrlon,ptrlev,ptrlat,ptrtim) ! input n2o   loss rate coeff
    real(r8) xcfc11i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc11 loss rate coeff
    real(r8) xcfc12i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc12 loss rate coeff

    real(r8) xch4(ptrlat, ptrlev)   ! input ch4   loss rate coeff indices changed
    real(r8) xn2o(ptrlat, ptrlev)   ! input n2o   loss rate coeff indices changed
    real(r8) xcfc11(ptrlat, ptrlev) ! input cfc11 loss rate coeff indices changed
    real(r8) xcfc12(ptrlat, ptrlev) ! input cfc12 loss rate coeff indices changed

    real(r8) xch4lv(ptrlat, plev)   ! input ch4   loss rate coeff interp to cam levels
    real(r8) xn2olv(ptrlat, plev)   ! input n2o   loss rate coeff interp to cam levels
    real(r8) xcfc11lv(ptrlat, plev) ! input cfc11 loss rate coeff interp to cam levels
    real(r8) xcfc12lv(ptrlat, plev) ! input cfc12 loss rate coeff interp to cam levels

    character(len=256) :: locfn    ! netcdf local filename to open 
!
!-----------------------------------------------------------------------
!
! Initialize
!
    nm = 1
    np = 2
!
! SPMD: Master reads dataset and does spatial interpolation on all time samples.
!       Spatial interpolents are broadcast.  All subsequent time interpolation is
!       done in every process.
!
    if (masterproc) then
       write(6,*)'CHEM_INIT_LOSS: greenhouse loss rates from: ', trim(bndtvg)
       call getfil(bndtvg, locfn)
       call wrap_open(locfn, 0, ncid_trc)
       write(6,*)'CHEM_INIT_LOSS: NCOPN returns id ',ncid_trc,' for file ',trim(locfn)
!
!------------------------------------------------------------------------
! Read tracer data
!------------------------------------------------------------------------
!
! Get dimension info
!
       call wrap_inq_dimid(ncid_trc, 'lat' , latid)
       call wrap_inq_dimid(ncid_trc, 'lev' , levid)
       call wrap_inq_dimid(ncid_trc, 'lon' , lonid)
       call wrap_inq_dimid(ncid_trc, 'time', timid)

       call wrap_inq_dimlen(ncid_trc, lonid, lonsiz)
       call wrap_inq_dimlen(ncid_trc, levid, levsiz)
       call wrap_inq_dimlen(ncid_trc, latid, latsiz)
       call wrap_inq_dimlen(ncid_trc, timid, timsiz)
!
! Check dimension info
!
       if (ptrlon/=1) then
          call endrun ('CHEM_INIT_LOSS: longitude dependence not implemented')
       endif
       if (lonsiz /= ptrlon) then
          write(6,*)'CHEM_INIT_LOSS: lonsiz=',lonsiz,' must = ptrlon=',ptrlon
          call endrun
       end if
       if (levsiz /= ptrlev) then
          write(6,*)'CHEM_INIT_LOSS: levsiz=',levsiz,' must = ptrlev=',ptrlev
          call endrun
       end if
       if (latsiz /= ptrlat) then
          write(6,*)'CHEM_INIT_LOSS: latsiz=',latsiz,' must = ptrlat=',ptrlat
          call endrun
       end if
       if (timsiz /= ptrtim) then
          write(6,*)'CHEM_INIT_LOSS: timsiz=',timsiz,' must = ptrtim=',ptrtim
          call endrun
       end if
!
! Determine necessary dimension and variable id's
!
       call wrap_inq_varid(ncid_trc, 'lat'    , latid)
       call wrap_inq_varid(ncid_trc, 'lev'    , levid)
       call wrap_inq_varid(ncid_trc, 'date'   , dateid)
       call wrap_inq_varid(ncid_trc, 'datesec', secid)
       call wrap_inq_varid(ncid_trc, 'TCH4'   , tch4id)
       call wrap_inq_varid(ncid_trc, 'TN2O'   , tn2oid)
       call wrap_inq_varid(ncid_trc, 'TCFC11' , tcfc11id)
       call wrap_inq_varid(ncid_trc, 'TCFC12' , tcfc12id)
!
! Obtain entire date and sec variables. Assume that will always
! cycle over 12 month data.
!
       call wrap_get_var_int(ncid_trc, dateid, date_tr)
       call wrap_get_var_int(ncid_trc, secid , sec_tr)

       if (mod(date_tr(1),10000)/100 /= 1) then
          call endrun ('(CHEM_INIT_LOSS): error when cycling data: 1st month must be 1')
       end if
       if (mod(date_tr(ptrtim),10000)/100 /= 12) then
          call endrun ('(CHEM_INIT_LOSS): error when cycling data: last month must be 12')
       end if
!
! Obtain input data latitude and level arrays.
!
       call wrap_get_var_realx(ncid_trc, latid, lati)
       call wrap_get_var_realx(ncid_trc, levid, pin )
!
! Convert input pressure levels to height (m).
! First convert from millibars to pascals.
!
       do k=1,ptrlev
          pin(k) = pin(k)*100.
          zi(k) = 7.0e3 * log (1.0e5 / pin(k))
       end do
!
! Convert approximate cam pressure levels to height (m).
!
       do k=1,plev
          zo (k) = 7.0e3 * log (1.0e5 / hypm(k))
       end do
!
! Convert cam model latitudes to degrees.
! Input model latitudes already in degrees.
!
       do j=1,plat
          lato(j) = clat(j)*45./atan(1.)
       end do
!
! Obtain all time samples of tracer data.
!
       call wrap_get_var_realx(ncid_trc, tch4id  , xch4i  )
       call wrap_get_var_realx(ncid_trc, tn2oid  , xn2oi  )
       call wrap_get_var_realx(ncid_trc, tcfc11id, xcfc11i)
       call wrap_get_var_realx(ncid_trc, tcfc12id, xcfc12i)
!
! Close netcdf file
!
       call wrap_close(ncid_trc)
!
!------------------------------------------------------------------------
! Interpolate tracer data to model grid
!------------------------------------------------------------------------
!
! Loop over all input times.
!
       do nt = 1, ptrtim
!
! Remove longitude and time index and switch level and latitude indices
! for the loss coefficients.
!
          do j=1,ptrlat
             do k=1,ptrlev
                xch4  (j,k) = xch4i  (1,k,j,nt)
                xn2o  (j,k) = xn2oi  (1,k,j,nt)
                xcfc11(j,k) = xcfc11i(1,k,j,nt)
                xcfc12(j,k) = xcfc12i(1,k,j,nt)
             end do
          end do
!
! Interpolate input data to model levels.
! If the CAM level is outside the range of the input data (this
! can happen only in troposphere) put zero for every latitude.
! Otherwise determine the input data levels bounding the current
! CAM level and interpolate.
!
          do ko=1,plev
             if (zo(ko) < zi(ptrlev)) then
                do j=1,ptrlat
                   xch4lv  (j,ko) = 0.0
                   xn2olv  (j,ko) = 0.0
                   xcfc11lv(j,ko) = 0.0
                   xcfc12lv(j,ko) = 0.0
                end do
                goto 50
             end if
             do ki=1,ptrlev-1
                if (zo(ko) < zi(ki) .and. zo(ko) >= zi(ki+1)) then
                   do j=1,ptrlat
                      xch4lv(j,ko) = xch4(j,ki) + (xch4(j,ki+1) - xch4(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xn2olv(j,ko) = xn2o(j,ki) + (xn2o(j,ki+1) - xn2o(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xcfc11lv(j,ko) = xcfc11(j,ki) + (xcfc11(j,ki+1) - xcfc11(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xcfc12lv(j,ko) = xcfc12(j,ki) + (xcfc12(j,ki+1) - xcfc12(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                   end do
                   goto 50
                endif
             end do
             call endrun ('(CHEM_INIT_LOSS): Error in vertical interpolation')
50           continue
          end do
!
! Interpolate input data to model latitudes.
! Determine the input data latitudes bounding the current CAM latitude and
! interpolate. Use last value from input data if the cam latitude is
! outside the range of the input data latitudes.
!
          do jo=1,plat
             if (lato(jo) <= lati(1)) then
                do k = 1, plev
                   tch4i(jo,k,nt)   = xch4lv(1,k)
                   tn2oi(jo,k,nt)   = xn2olv(1,k)
                   tcfc11i(jo,k,nt) = xcfc11lv(1,k)
                   tcfc12i(jo,k,nt) = xcfc12lv(1,k)
                end do
             else if (lato(jo) >= lati(ptrlat)) then
                do k = 1, plev
                   tch4i(jo,k,nt)   = xch4lv(ptrlat,k)
                   tn2oi(jo,k,nt)   = xn2olv(ptrlat,k)
                   tcfc11i(jo,k,nt) = xcfc11lv(ptrlat,k)
                   tcfc12i(jo,k,nt) = xcfc12lv(ptrlat,k)
                end do
             else
                do ji=1,ptrlat-1
                   if ( (lato(jo) > lati(ji)) .and. (lato(jo) <= lati(ji+1))) then
                      do k=1,plev
                         tch4i(jo,k,nt) = xch4lv(ji,k) + (xch4lv(ji+1,k) - xch4lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tn2oi(jo,k,nt) = xn2olv(ji,k) + (xn2olv(ji+1,k) - xn2olv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tcfc11i(jo,k,nt) = xcfc11lv(ji,k) + (xcfc11lv(ji+1,k) - xcfc11lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tcfc12i(jo,k,nt) = xcfc12lv(ji,k) + (xcfc12lv(ji+1,k) - xcfc12lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                      end do
                      goto 90
                   endif
                end do
             end if
             write (6,*)'(CHEM_INIT_LOSS): Error in horizontal interpolation'
90           continue
          end do

       end do                 ! end loop over time samples

    endif                     ! end of masterproc
#if ( defined SPMD )
    call mpibcast (date_tr, ptrtim, mpiint, 0, mpicom)
    call mpibcast (sec_tr , ptrtim, mpiint, 0, mpicom)
    call mpibcast (tch4i  , plev*plat*ptrtim, mpir8, 0, mpicom)
    call mpibcast (tn2oi  , plev*plat*ptrtim, mpir8, 0, mpicom)
    call mpibcast (tcfc11i, plev*plat*ptrtim, mpir8, 0, mpicom)
    call mpibcast (tcfc12i, plev*plat*ptrtim, mpir8, 0, mpicom)
#endif

!
! Initial time interpolation between December and January
!
    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    n = 12
    np1 = 1
    call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
    call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
    if (calday <= cdaytrp .or. calday > cdaytrm) then
       do j=1,plat
          do k=1,plev
             tch4m  (j,k,nm) = tch4i  (j,k,n)
             tn2om  (j,k,nm) = tn2oi  (j,k,n)
             tcfc11m(j,k,nm) = tcfc11i(j,k,n)
             tcfc12m(j,k,nm) = tcfc12i(j,k,n)
             tch4m  (j,k,np) = tch4i  (j,k,np1)
             tn2om  (j,k,np) = tn2oi  (j,k,np1)
             tcfc11m(j,k,np) = tcfc11i(j,k,np1)
             tcfc12m(j,k,np) = tcfc12i(j,k,np1)
          end do
       end do
       goto 10
    end if
!
! Initial normal interpolation between consecutive time slices.
!
    do n=1,ptrtim-1
       np1 = n + 1
       call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
       call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
       if (calday > cdaytrm .and. calday <= cdaytrp) then
          do j=1,plat
             do k=1,plev
                tch4m  (j,k,nm) = tch4i  (j,k,n)
                tn2om  (j,k,nm) = tn2oi  (j,k,n)
                tcfc11m(j,k,nm) = tcfc11i(j,k,n)
                tcfc12m(j,k,nm) = tcfc12i(j,k,n)
                tch4m  (j,k,np) = tch4i  (j,k,np1)
                tn2om  (j,k,np) = tn2oi  (j,k,np1)
                tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                tcfc12m(j,k,np) = tcfc12i(j,k,np1)
             end do
          end do
          goto 10
       end if
    end do
    write(6,*)'CHEM_INIT_LOSS: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
    call endrun
!
! Data positioned correctly
!
10  continue


    return
  end subroutine chem_init_loss

!===============================================================================
  subroutine chem_timestep_init
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Time interpolate chemical loss rates to current time, updating the
! bounding month arrays with new data if necessary
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------

    use commap
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual

    implicit none

#include <comctl.h>
!-----------------------------------------------------------------------
!
! Local workspace
!
    integer ntmp               ! temporary
    integer j,k                ! indices
    integer :: yr, mon, day    ! components of a date
    integer :: ncdate          ! current date in integer format [yyyymmdd]
    integer :: ncsec           ! current time of day [seconds]
    real(r8) :: calday         ! current calendar day
    real(r8) fact1, fact2      ! time interpolation factors
!-----------------------------------------------------------------------
!
! If model time is past current forward timeslice, obtain the next
! timeslice for time interpolation.  Messy logic is for
! interpolation between December and January (np1.eq.1).
!
    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day
    
    if (calday > cdaytrp .and.  .not. (np1 == 1 .and. calday > cdaytrm)) then
       np1 = mod(np1,12) + 1
       if (np1 > ptrtim) then
          call endrun ('CHEMINT: Attempt to access bad month')
       end if
       cdaytrm = cdaytrp
       call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
       if (np1 == 1 .or. calday <= cdaytrp) then
          ntmp = nm
          nm   = np
          np   = ntmp
          do j=1,plat
             do k=1,plev
                tch4m  (j,k,np) = tch4i  (j,k,np1)
                tn2om  (j,k,np) = tn2oi  (j,k,np1)
                tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                tcfc12m(j,k,np) = tcfc12i(j,k,np1)
             end do
          end do
       else
          write(6,*)'CHEMINT: Input data for date',date_tr(np1), &
             ' sec ',sec_tr(np1), 'does not exceed model date', &
             ncdate,' sec ',ncsec,' Stopping.'
          call endrun
       end if
    end if
!
! Determine time interpolation factors.  1st arg says we are cycling yearly data
!
    call getfactors (.true., np1, cdaytrm, cdaytrp, calday, &
                     fact1, fact2, 'SULFINT:')
!
! Do time interpolation
!
    do j=1,plat
       do k=1,plev
          tch4(j,k)   = tch4m  (j,k,nm)*fact1 + tch4m  (j,k,np)*fact2
          tn2o(j,k)   = tn2om  (j,k,nm)*fact1 + tn2om  (j,k,np)*fact2
          tcfc11(j,k) = tcfc11m(j,k,nm)*fact1 + tcfc11m(j,k,np)*fact2
          tcfc12(j,k) = tcfc12m(j,k,nm)*fact1 + tcfc12m(j,k,np)*fact2
       end do
    end do

    return
  end subroutine chem_timestep_init

!===============================================================================
  subroutine chem_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial mass mixing ratios of CH4, N2O, CFC11 and CFC12.
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name
    real(r8), intent(out) :: q(plon,plev,plat)   !  mass mixing ratio
!-----------------------------------------------------------------------

#ifdef GHG_2D
    if ( name == 'N2O' ) then
       q = rmwn2o * n2ovmr_2d(1,1)          !xinxg change to (1,1), because currently not used
       return
    else if ( name == 'CH4' ) then
       q = rmwch4 * ch4vmr_2d(1,1)
       return
    else if ( name == 'CFC11' ) then
       q = rmwf11 * f11vmr_2d(1,1)
       return
    else if ( name == 'CFC12' ) then
       q = rmwf12 * f12vmr_2d(1,1)
       return
    end if

#else
    if ( name == 'N2O' ) then
       q = rmwn2o * n2ovmr
       return
    else if ( name == 'CH4' ) then
       q = rmwch4 * ch4vmr
       return
    else if ( name == 'CFC11' ) then
       q = rmwf11 * f11vmr
       return
    else if ( name == 'CFC12' ) then
       q = rmwf12 * f12vmr
       return
    end if
#endif

  end subroutine chem_init_cnst

#ifdef GEOSCHEM

  subroutine geoschem_init_cnst(name, q)
!-----------------------------------------------------------------------
!
! Purpose:
! Set initial mass mixing ratios of CO, ... for tracers in GEOSCHEM
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name

    real(r8), intent(out) :: q(plon,plev,plat)   !  mass mixing ratio
!-----------------------------------------------------------------------

    if ( trim(name) == 'NULL' ) then
       q = 0.d0
       return
    else if ( name == 'CH4' ) then
       q = 6.31d-7 ! 1140 ppb converted to kg/kg based on IPCC 1950 values
       return
    else if ( name == 'CO2' ) then
       q = 0.d0 !4.71d-4 ! 310 ppm converted to kg/kg based on 1950 IPCC values
       return
    else if ( name == 'NOx' ) then
       q = 1.d-35
       return
    else
       q = 1.d-20
    end if

  end subroutine geoschem_init_cnst
#endif

end module chemistry
