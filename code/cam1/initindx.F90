#include <misc.h>
subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef TDK_CU
  use constituents, only: pcnst, ppcnst, cnst_add, advected, cnst_chk_dim, cnst_name, &
                          nonadvec, ixcldw
#else
  use constituents, only: pcnst, ppcnst, cnst_add, advected, cnst_chk_dim, cnst_name
#endif
  use phys_buffer,  only: pbuf_init
!------ zf 2008.03
#ifdef MOZART2
  use chemistry,    only: trace_gas, chem_register, geoschem_register
#else
  use chemistry,    only: trace_gas, chem_register
#endif
#ifdef BCCCHEM  
  use chemistry,    only: bcc_register
#endif

#ifdef CO2  
  use chemistry,    only: co2_register
#endif

#ifdef WACCM  
  use chemistry,    only: waccm_chem_register
  use exbdrift,     only: exbdrift_register
#endif

#ifdef GEOSCHEM
  use chemistry,    only: gc_register
#endif
  use cldcond,      only: cldcond_register
  use physconst,    only: mwdry, cpair, mwh2o, cph2o
  use tracers, only: tracers_register
!  use constituents, only: dcconnam, sflxnam, hadvnam, vadvnam, fixcnam, 
  use constituents, only: dcconnam, sflxnam, tendnam, tottnam
  use check_energy, only: check_energy_register
  use aerosol_intr, only: aerosol_register_cnst
#ifdef BCC_DRAG
  use gw_drag,            only: gw_drag_register
#endif
#if ( defined BFB_CAM_SCAM_IOP )
  use iop
#endif
  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index 

#ifdef TDK_CU
  logical, parameter :: cldw_adv=.false.  ! true => cloud water is treated as advected tracer
#endif
!-----------------------------------------------------------------------

! Initialize physics buffer
  call pbuf_init()

! Register water vapor.
! ***** N.B. ***** This must be the first call to cnst_add so that
!                  water vapor is constituent 1.
  call cnst_add('Q', advected, mwh2o, cph2o, 1.E-12_r8, mm, &
                longname='Specific humidity', readiv=.true.)
#ifdef TDK_CU

! Register cloud water and determine index (either advected or non-adv).
  if (cldw_adv) then
     call cnst_add('CWAT', advected, mwdry, cpair, 0._r8, ixcldw, &
                    longname='Total Grid box averaged Condensate Amount (liquid + ice)')
  else
     call cnst_add('CWAT', nonadvec, mwdry, cpair, 0._r8, ixcldw, &
                    longname='Total Grid box averaged Condensate Amount (liquid + ice)')
  endif

#endif

!
! Register cloud water
  call cldcond_register()
!
! Register chemical constituents
  if (trace_gas) then
     call chem_register()
  endif
!
! register aerosols
  call aerosol_register_cnst()
!
! Register mozart2 constituents
#ifdef MOZART2
  call geoschem_register                 ! zf 2008.03
#endif

#ifdef CO2
  call co2_register                 ! zf 2008.03
#endif

#ifdef BCCCHEM
  call bcc_register
#endif

#ifdef WACCM
  call waccm_chem_register()

  ! Initialize e and b fields
  call exbdrift_register()
#endif

! register gravity wave drag
#ifdef BCC_DRAG
  call gw_drag_register()
#endif

! Register geoschem constituents
#ifdef GEOSCHEM
  call gc_register    
#endif

! Register advected test tracers and determine starting index
  call tracers_register()

!
! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()
!
! Set default names for non-water advected and non-advected tracers
! Set names of advected and non-advected tracer diagnostics
!
  do m=1,ppcnst
     dcconnam(m) = 'DC'//cnst_name(m)
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do
  do m=1,pcnst
!     hadvnam(m)  = 'HA'//cnst_name(m)
!     vadvnam(m)  = 'VA'//cnst_name(m)
!     fixcnam(m)  = 'DF'//cnst_name(m)
     tendnam(m)  = 'TE'//cnst_name(m)
     tottnam(m)  = 'TA'//cnst_name(m)
  end do

#if ( defined BFB_CAM_SCAM_IOP )
  do m=1,pcnst
     alphanam(m) = 'AFIX'//cnst_name(m)
     alphanam(m)=to_lower(alphanam(m))
     dqfxnam(m) = 'DQFX'//cnst_name(m)
     dqfxnam(m) = to_lower(dqfxnam(m))
  end do
#endif
  call check_energy_register()

end subroutine initindx
