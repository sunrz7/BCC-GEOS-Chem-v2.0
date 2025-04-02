MODULE BCCGCHP_VARS_MOD
!
! !USES:
!
!  USE GIGC_Type_Mod                                  ! Derived type defs
  USE input_opt_mod                                  ! Input Options obj
  USE state_chm_mod                                  ! Chemistry State obj
  USE state_met_mod                                  ! Meteorology State obj
  USE species_mod,   ONLY : species

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!

!
! !PRIVATE MEMBER FUNCTIONS:
!

!
! !PRIVATE TYPES:
!
  ! Objects for GEOS-Chem
  TYPE(OptInput)                   :: Input_Opt      ! Input Options
  TYPE(MetState)                   :: State_Met      ! Meteorology state
  TYPE(ChmState)                   :: State_Chm      ! Chemistry state
  TYPE(Species), POINTER           :: ThisSpc => NULL()

  ! Scalars

  ! List here GEOS-Chem tracer names and corresponding names to be assigned
  ! to the AERO bundle (if GC is the AERO provider). The names in the AERO
  ! bundle must be the names that are expected by the irradiation component:
  ! - OCphobic, OCphilic, BCphobic, and BCphilic for hydrophobic and hydrophilic
  !   organic and black carbon, respectively
  ! - SO4 for SO4
  ! - du001 - du005 for the following five dust bins (see DU_GridComp.rc in
  !   GOCART):
  !   radius_lower: 0.1 1.0 1.8 3.0 6.0
  !   radius_upper: 1.0 1.8 3.0 6.0 10.0
  !
  !   The GEOS-Chem dust bins are: 
  !   Reff: 0.7 1.4 2.4 4.5
  !   Those become simply mapped onto the GOCART dust bins 1-4 (du001 ... du004).
  !
  ! - ss001-ss005 for the following five sea salt aerosol bins (see SS_GridComp.rc
  !   in GOCART):
  !   radius_lower: 0.03 0.1 0.5 1.5 5.0
  !   radius_upper: 0.1  0.5 1.5 5.0 10.0
  !
  !   The GEOS-Chem sea salt aerosols are (SALA and SALC):
  !   radius_lower: 0.01 0.5
  !   radius_upper: 0.5  8.0
  !   SALA becomes mapped onto ss001 and ss002, and SALC onto ss003, ss004, ss005. 
  !   For now, we assume uniform size distribution within the GEOS-Chem bins, i.e.
  !   the GEOS-Chem size bins are evenly split into the GOCART bins. The fractions can
  !   be specified below.
  !   At some point, we may revisit these fractions (at least take into account the
  !   log-normal behavior of the aerosol distribution)
  INTEGER, PARAMETER           :: NumAERO = 11
  CHARACTER(LEN=16)   :: GcNames(NumAero) = &
                                  (/ 'DST1',     'DST2',     'DST3',     'DST4',     &
                                     'SALA',     'SALC',     'BCPO',     'BCPI',     &
                                     'OCPO',     'OCPI',     'SO4 '                   /)

  CHARACTER(LEN=16)   :: AeroNames(NumAero) = &
                                  (/ 'du001   ', 'du002   ', 'du003   ', 'du004   ', &
                                     'ss001   ', 'ss003   ', 'BCphobic', 'BCphilic', &
                                     'OCphobic', 'OCphilic', 'SO4     '               /)

  ! Fraction of SALA in ss001 and ss002, respectively
  CHARACTER(LEN=16)   :: SALAnames(2) = (/ 'ss001', 'ss002' /)
  REAL, PARAMETER              :: SALAsplit(2) = (/  0.2,     0.8    /)
  ! Fraction of SALC in ss003, ss004, and ss005.
  CHARACTER(LEN=16)   :: SALCnames(3) = (/ 'ss003', 'ss004' , 'ss005' /)
  REAL, PARAMETER              :: SALCsplit(3) = (/  0.13,    0.47,     0.4    /) 

  ! Prefix of the tracer and species names in the internal state. Those have to match
  ! the prefixes given in GEOSCHEMchem_Registry.rc. 
  CHARACTER(LEN=4), PARAMETER  :: SPFX = 'SPC_'
 
CONTAINS
END MODULE BCCGCHP_VARS_MOD
 
