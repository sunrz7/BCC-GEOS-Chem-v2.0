!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! sunrz 2023/11
! !MODULE: gigc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  run, and finalize methods for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Chunk_Mod
!
! !USES:
!      
  use gchp_utils
  use errcode_mod
  use precision_mod
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
!xiaolu
!  PUBLIC :: GIGC_Chunk_Run
  PUBLIC :: GIGC_Chunk_Run
!  PUBLIC :: GIGC_Chunk_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
!  PRIVATE :: GIGC_Cap_Tropopause_Prs
!  PRIVATE :: GIGC_Revert_Units
!  PRIVATE :: GIGC_Assert_Units
!
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_init
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_INIT is the ESMF init method for
!  the Grid-Independent GEOS-Chem (aka "GIGC").  This routine calls lower-
!  level routines to allocate arrays and read input files.
!\\
!\\
! !INTERFACE:

!
  SUBROUTINE GIGC_Chunk_Init( nymdB,         nhmsB,      nymdE,           &
                              nhmsE,         tsChem,     tsDyn,           &
                              tsRad,         lonCtr,     latCtr,       lonEdge,latEdge,   &
                              Input_Opt,     State_Chm,  State_Diag,      &
                              State_Grid,    State_Met,  HcoConfig,       &
                               RC )
!
! !USES:
!
    use chemistry_mod,           only : init_chemistry
    use emissions_mod,           only : emissions_init
    use gc_environment_mod
    use gc_grid_mod,             only : setgridfromctr,setgridfromctredges
    !use gchp_historyexports_mod, only : historyconfigobj
    use hco_types_mod,           only : configobj
    use input_mod,               only : read_input_file
    use input_opt_mod,           only : optinput, set_input_opt
    use linear_chem_mod,         only : init_linear_chem
    use linoz_mod,               only : linoz_read
    use physconstants,           only : pi_180
    use pressure_mod,            only : init_pressure
    use roundoff_mod,            only : roundoff
    use state_chm_mod,           only : chmstate, ind_
    use state_diag_mod,          only : dgnstate
    use state_grid_mod,          only : grdstate, init_state_grid
    use state_met_mod,           only : metstate
    use time_mod,                only : set_timesteps
    use ucx_mod,                 only : init_ucx
    use unitconv_mod,            only : convert_spc_units
    ! for conc initialization
    use species_mod,             only : species
    ! sunrz for diag 2024/3
    use diaglist_mod,       only: init_diaglist, print_diaglist
    use taggeddiaglist_mod, only: init_taggeddiaglist, print_taggeddiaglist
    use diaglist_mod, only: dgnlist
    use taggeddiaglist_mod, only: taggeddgnlist
    ! sunrz end
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsRad       ! Chemistry timestep [s]
    REAL(selected_real_kind(3,25)), INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(selected_real_kind(3,25)), INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
    REAL(selected_real_kind(3,25)), INTENT(IN)     :: lonEdge(:,:)! Lonedges[radians]
    REAL(selected_real_kind(3,25)), INTENT(IN)     :: latEdge(:,:)! Latedges[radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag     ! Diag State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    TYPE(ConfigObj),     POINTER       :: HcoConfig      ! HEMCO config obj
    !TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config obj
    ! sunrz for diag 2024/3
    type(DgnList)                              :: Global_DiagList
    type(TaggedDgnList)                        :: Global_TaggedDiagList
    ! sunrz end 
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: I, J, L, STATUS
    CHARACTER(LEN=256)             :: I_am
    TYPE(Species),      POINTER    :: ThisSpc     ! Species pointer for init bg values
    INTEGER                        :: IND         ! Current species index
    !=======================================================================
    ! GIGC_CHUNK_INIT begins here 
    !=======================================================================

    ! Error trap
    I_am = 'GIGC_CHUNK_INIT (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    ! Update Input_Opt with timing fields
    ! We will skip defining these in READ_INPUT_FILE
    Input_Opt%NYMDb   = nymdB           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSb   = nhmsB           ! hhmmss   @ end   of simulation
    Input_Opt%NYMDe   = nymdE           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSe   = nhmsE           ! hhmmss   @ end   of simulation
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_RAD  = INT( tsRad  )   ! RRTMG     timestep [sec]
    !Input_Opt%Max_AdvectSpc = 600
    ! Read geoschem_config.yml at very beginning of simulation on every CPU
    CALL Read_Input_File( Input_Opt, State_Grid, RC )

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt, State_Grid, RC )

    ! Set maximum number of levels in the chemistry grid
    State_Grid%MaxChemLev  = State_Grid%MaxStratLev

    ! In the ESMF/MPI environment, we can get the total overhead ozone
    ! either from the met fields (GCHPsa) or from the Import State (GEOS-5)
    Input_Opt%USE_O3_FROM_MET = .TRUE.

    ! Read LINOZ climatology
    IF ( Input_Opt%LLINOZ ) THEN
       CALL Linoz_Read( Input_Opt, RC )
    ENDIF
    ! sunrz for diag 2024/3
        CALL Init_DiagList(Input_Opt%AmIRoot, "HISTORY.rc",Global_DiagList, RC)
        if(RC /= GC_SUCCESS) then
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Init_DiagList)"
            return
        endif

        CALL Print_DiagList(Input_Opt%AmIRoot, Global_DiagList, RC)

        ! Initialize the TaggedDiag_List
        CALL Init_TaggedDiagList(Input_Opt%amIroot, Global_DiagList,  &
                                 Global_TaggedDiagList,   RC         )
        IF ( RC /= GC_SUCCESS ) THEN
           write(6,*) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Init_TaggedDiagList)"
           return
        ENDIF

        CALL Print_TaggedDiagList(Input_Opt%amIRoot, Global_TaggedDiagList, RC)
    ! sunrz end
    ! Allocate all lat/lon arrays
    CALL GC_Allocate_All( Input_Opt, State_Grid, RC )

    ! Set grid based on passed mid-points
    !CALL SetGridFromCtr( Input_Opt, State_Grid, lonCtr, latCtr, RC )
    CALL SetGridFromCtrEdges( Input_Opt, State_Grid, lonCtr, latCtr, lonEdge, latEdge, RC )!sunrz 2024/2
    ! Set GEOS-Chem timesteps on all CPUs
    CALL Set_Timesteps( Input_Opt  = Input_Opt,                              &
                        Chemistry  = Input_Opt%TS_CHEM,                      &
                        Convection = Input_Opt%TS_CONV,                      &
                        Dynamics   = Input_Opt%TS_DYN,                       &
                        Emission   = Input_Opt%TS_EMIS,                      &
                        Radiation  = Input_Opt%TS_RAD,                       &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,                  &
                                          Input_Opt%TS_CONV ),               &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( Global_DiagList, Global_TaggedDiagList, Input_Opt, &
                           State_Chm, State_Diag, State_Grid, State_Met, RC )

    CALL GC_Init_Extra( Global_DiagList, Input_Opt,    &
                        State_Chm, State_Diag, State_Grid, RC )
    ! Set initial species units to internal state units, the same
    ! units as the restart file values. Note that species concentrations
    ! are all still zero at this point since internal state values are not
    ! copied to State_Chm%Species%Conc until Run (post-initialization).
    State_Chm%Spc_Units = 'kg/kg dry'
    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, RC )
    ENDIF

    ! Initialize HEMCO
!    CALL EMISSIONS_INIT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
!                         HcoConfig=HcoConfig )

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    ! Initialize UCX routines
    CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )
!be cautious sunrz
    IF ( Input_Opt%LINEAR_CHEM ) THEN
       CALL INIT_LINEAR_CHEM( Input_Opt, State_Chm, State_Met, State_Grid, RC )
    ENDIF


    RC = GC_Success

    do I = 1, State_Chm%nSpecies
      ThisSpc => State_Chm%SpcData(I)%Info

      if(trim(ThisSpc%Name) == '') cycle
      IND = IND_(trim(ThisSpc%Name))
      if(IND < 0) cycle

      ! Initialize using background values from species database.
      call SET_BACKGROUND_CONC( Input_Opt%AmIRoot, ThisSpc, &
                               State_Chm, State_Met, State_Grid, &
                               Input_Opt, IND, RC)

      ThisSpc => NULL()
    enddo
    write(6,*) "Chunk_Init: Set background concentrations"

  END SUBROUTINE GIGC_Chunk_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_run
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_RUN is the ESMF run method for
!  the Grid-Independent GEOS-Chem (aka "GIGC").  This routine is the driver
!  for the following operations:
!
! \begin{itemize}
! \item Dry deposition
! \item Chemistry
! \end{itemize}
!
! !INTERFACE:
!

        !-------------------------------------
        !xiaolu,2017/10/05: 
        ! 1.change GIGC_CHUNK_RUN to GIGC_CHUNK_RUN_BC: chemistry and wet dep.
        ! 2.add GIGC_CHUNK_RUN_AC: dry dep.
        !-------------------------------------
!  SUBROUTINE GIGC_Chunk_Run_BC( am_I_Root,  IM,    JM,        LM,         &
!                             nymd,      nhms,      year,      month,      &
!                             day,       dayOfYr,   hour,      minute,     &
!                             second,    utc,       hElapsed,  Input_Opt,  &
!                             State_Chm, State_Met, Phase,     IsChemTime, &
!                             RC                                            )
!

  SUBROUTINE GIGC_Chunk_Run   ( am_I_Root,                                   &
                             nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime, IsRadTime,              &
                             RC )
!
! !USES:
!
    ! GEOS-Chem state objects
    use input_opt_mod,      only : optinput
    use state_chm_mod,      only : chmstate
    use state_diag_mod
    use state_grid_mod,     only : grdstate
    use state_met_mod,      only : metstate

    ! geos-chem components
    use chemistry_mod,      only : do_chemistry, recompute_od
    use convection_mod,     only : do_convection
    use drydep_mod,         only : do_drydep
    use emissions_mod,      only : emissions_run
    use mixing_mod,         only : do_tend, do_mixing
    use wetscav_mod,        only : setup_wetscav, do_wetdep

    ! hemco components (eventually moved to a separate gridcomp?)
    use hco_state_gc_mod,   only : hcostate, extstate
    use hco_interface_common, only : sethcotime
    use hco_interface_gc_mod, only : compute_sflx_for_vdiff

    ! specialized subroutines
    use calc_met_mod,       only : airqnt
    use calc_met_mod,       only : set_dry_surface_pressure
    use calc_met_mod,       only : set_clock_tracer
    use calc_met_mod,       only : gchp_cap_tropopause_prs
    use set_global_ch4_mod, only : set_ch4
    use modis_lai_mod,      only : compute_xlai,get_xlainative_from_hemco
    use pbl_mix_mod,        only : compute_pbl_height
    use pressure_mod,       only : set_floating_pressures
    use toms_mod,           only : compute_overhead_o3
    use ucx_mod,            only : set_h2o_trac
    use vdiff_mod,          only : max_pblht_for_vdiff

    ! utilities
    use errcode_mod
    use hco_error_mod
    use pressure_mod,       only : accept_external_pedge
    use state_chm_mod,      only : ind_
    use time_mod,           only : accept_external_date_time
    use unitconv_mod,       only : convert_spc_units, print_global_species_kg

    ! diagnostics
    use diagnostics_mod,    only : zero_diagnostics_startoftimestep
    use diagnostics_mod,    only : set_diagnostics_endoftimestep
    use aerosol_mod,        only : set_aermass_diagnostic

    use species_mod,   only : species
    use pmgrid,             only: masterproc!sunrz 2023/12
!xiaolu, add olson map initialization
!    use gchp_landmap_mod     ! computes ireg, iland, iuse from olson map
    use olson_landmap_mod,  only : compute_olson_landmap, init_landtypefrac!sunrz 2024/2
!    use hco_state_gc_mod,   only : hcostate !sunrz 2024/3
    use pmgrid,        only: masterproc,iam
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    INTEGER,        INTENT(IN)    :: Phase       ! Run phase (-1, 1 or 2)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry? 
    LOGICAL,        INTENT(IN)    :: IsRadTime   ! Time for RRTMG? 
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on root CPU?
!    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref to this GridComp
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code


    REAL*8                         :: DT
    CHARACTER(LEN=255)             :: I_am, OrigUnit
    INTEGER                        :: STATUS, HCO_PHASE, RST

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend
    LOGICAL                        :: DoTurb
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep
    LOGICAL                        :: DoRad

    ! First call?
    LOGICAL, SAVE                  :: FIRST    = .TRUE.
    LOGICAL, SAVE                  :: FIRST_RT = .TRUE. ! RRTMG

    ! # of times this routine has been called. Only temporary for printing 
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings
    LOGICAL                        :: SetStratH2O

    ! For RRTMG
    INTEGER                        :: N,i

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                  :: scaleMR = .FALSE.

    ! Debug variables
    INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1
    !=======================================================================
    ! GIGC_CHUNK_RUN begins here 
    !=======================================================================

    ! Error trap
    I_am = 'GIGC_CHUNK_RUN (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

if (maxval(State_Chm%Species(213)%Conc)>1) then
write(*,*) 'sunrz check101 before so4=',maxval(State_Chm%Species(213)%Conc),maxval(State_Chm%Species(180)%Conc)
endif
    ! Populate grid box areas from gc_grid_mod.F on first call. We have to do 
    ! this in the run phase and not in initialize because it seems like the
    ! AREA pointer (imported from superdynamics) is only properly filled in
    ! run.
!    IF ( FIRST ) THEN
!       AREA_M2 = State_Met%AREA_M2
!    ENDIF

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation 
    !
    ! Here, we use the following operator sequence:
    ! 
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2 
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2     
    ! 
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the input.geos file is enabled. If the physics component already
    ! covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    !=======================================================================

    ! By default, do processes as defined in input.geos. DoTend is defined
    ! below. 
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
! sunrz 2023/11 in BC force drydep to false
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime    ! chemistry time step
    !DoEmis   = IsChemTime                         ! chemistry time step
    DoEmis   = .false.                            ! sunrz 2023/12
    DoTurb   = .false.                            ! dynamic time step
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
!    DoChem   = .false.
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 
    DoRad    = Input_Opt%LRAD  .AND. IsRadTime    ! radiation time step 
    ! Only do selected processes for given phases: 
    ! Phase 1: disable turbulence, chemistry and wet deposition. 
    IF ( Phase == 1 ) THEN
       DoTurb   = .FALSE.
       DoChem   = .FALSE.
       DoWetDep = .FALSE.

    ! Phase 2: disable convection, drydep and emissions. 
    ELSEIF ( Phase == 2 ) THEN
       DoConv   = .FALSE.
       DoDryDep = .FALSE.
       DoEmis   = .FALSE. 
    ENDIF

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    DoTend = .false.!( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB
    ! DoTend = .FALSE. !sunrz 2023/12
    ! testing only
    IF ( am_I_Root .and. NCALLS < 10 ) THEN 
       write(*,*) 'GEOS-Chem phase ', Phase, ':'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
    ENDIF
    ! Zero out certain State_Diag arrays. This should not be done in a phase 2
    ! call since this can erase diagnostics filled during phase 1 (e.g., drydep)
    ! (ckeller, 1/21/2022).
    CALL Zero_Diagnostics_StartOfTimestep( Input_Opt, State_Diag, RC )!sunrz add for diag 2024/3
    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = nymd,       &  
                                    value_NHMS     = nhms,       &  
                                    value_YEAR     = year,       &  
                                    value_MONTH    = month,      &  
                                    value_DAY      = day,        &  
                                    value_DAYOFYR  = dayOfYr,    &  
                                    value_HOUR     = hour,       &  
                                    value_MINUTE   = minute,     &  
                                    value_HELAPSED = hElapsed,   & 
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Set HEMCO time
!    CALL SetHcoTime ( HcoState,   ExtState,   year,    month,   day,   &
!                      dayOfYr,    hour,       minute,  second,  DoEmis,  RC )

    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge( State_Met  = State_Met,   &
                                State_Grid = State_Grid,  &
                                RC         = RC          )
    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    ! and compute avg surface pressures near polar caps
!---------------
!xiaolu comment,201703,as it is passed from BCC
!    CALL SET_DRY_SURFACE_PRESSURE( State_Met, 1 )
!---------------

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    ! and average as polar caps
!-----------------
!xiaolu comment,2017/03,as it is passed from BCC
!    CALL SET_DRY_SURFACE_PRESSURE( State_Met, 2 )
!----------------

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
!    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
!    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )
    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ! Define airmass and related quantities
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .false.)!sunrz 2023/12
    ! Initialize/reset wetdep after air quantities computed
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    ENDIF
    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GCHP_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )

    ! Call PBL quantities. Those are always needed
    ! xiaolu
    ! XIAOLU comment COMPUTE_PBL_HEIGHT as PBLH is from BCC dynamics
    ! CALL COMPUTE_PBL_HEIGHT( State_Met )
    ! Update clock tracer if relevant
    IF (  IND_('CLOCK','A') > 0 ) THEN
       CALL Set_Clock_Tracer( State_Chm, State_Grid )
    ENDIF
    ! Call PBL quantities. Those are always needed
! sunrz 2023/11
!    CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )
    !write(6,*) 'sunrz check so4=',maxval(State_Chm%Species(213)%Conc),minval(State_Chm%Species(213)%Conc)
    !stop
    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
       IF ( Input_Opt%LSETH2O ) THEN
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, &
                          State_Grid,  State_Met, RC )

      ! Only force strat once if using UCX
      ! FIXME: hplin this needs to be handled in another way for nested domains
      ! --
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF
    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    !HCO_PHASE = 1
    !CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
    !                    State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1/-1                              !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection (in v/v)
    ! 
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in geoschem_config.yml. This should only be done if convection is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do convection now' 

       CALL DO_CONVECTION ( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
 
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF  
    !=======================================================================
    ! 2. Dry deposition.
    !
    ! This calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       ! testing only
       if(am_I_Root.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif
       ! Make sure tracers are in kg
!       CALL ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met, State_Chm, RC )

      !----------------------------
      !For BCC, we need to initialize olsen map here
      !xiaolu,2017/08/18
      !----------------------------
!       CALL Compute_BCCAVIM_Landmap_GCHP( am_I_Root, State_Met, RC)
       ! Calculate drydep rates 
!       CALL Do_DryDep   ( am_I_Root = am_I_Root,            & ! Root CPU?
!                          Input_Opt = Input_Opt,            & ! Input Options
!                          State_Chm = State_Chm,            & ! Chemistry State
!                          State_Met = State_Met,            & ! Met State
!                          RC        = RC                   )  ! Success?

       ! Revert units
!       Call GIGC_Revert_Units( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       ! Do dry deposition
    CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )!sunrz 2024/2
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )

       CALL Do_DryDep ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF ! Do drydep
    !=======================================================================
    ! 3. Emissions (HEMCO). HEMCO must be called on first time step to make
    ! sure that the HEMCO data lists are all properly set up. 
    !=======================================================================
    IF ( DoEmis ) THEN

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do emissions now'

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Emissions done!'
    ENDIF
    !CALL Init_LandTypeFrac( Input_Opt, State_Grid, State_Met, RC ) !sunrz return to olson
    !CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )!sunrz 2024/2
    ! Initialize the State_Met%XLAI_NATIVE field from HEMCO
    !CALL Get_XlaiNative_from_HEMCO( Input_Opt, State_Grid, State_Met, RC )
    ! Calculate MODIS leaf area indexes needed for dry deposition
    !CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )
    ! Check that units are correct
    !ASSERT_(GIGC_Assert_Units(am_I_Root, State_Chm))

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry 
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    ! Subroutine DO_TEND operates in mass units, e.g. the tracers must be 
    ! in [kg].
    ! SDE 2016-04-05: Input units should now be v/v dry. 
    !=======================================================================
    IF ( DoTend ) THEN 
  
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'
 
       ! Make sure tracers are in v/v
!       CALL ConvertSpc_KgKgDry_to_VVDry( am_I_Root, State_Chm, RC )

       ! Get emission time step [s]. 
       !ASSERT_(ASSOCIATED(HcoState))
!       DT = HcoState%TS_EMIS 
       DT = Input_Opt%TS_EMIS
       ! Apply tendencies over entire PBL. Use emission time step.
!       CALL DO_TEND ( am_I_Root, Input_Opt, State_Met, State_Chm, .FALSE., RC, DT=DT )
       !ASSERT_(RC==GC_SUCCESS)

       ! Revert units
!       Call GIGC_Revert_Units( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       !ASSERT_(RC==GC_SUCCESS)
       ! Apply tendencies over entire PBL. Use emission time step.
!       CALL DO_TEND( Input_Opt, State_Chm, State_Diag, &
!                     State_Grid, State_Met, .FALSE., RC, DT=DT )

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT 
 
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 

    ENDIF ! Tendencies 

    ! Check that units are correct
    !ASSERT_(GIGC_Assert_Units(am_I_Root, State_Chm))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2                                !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence (v/v)
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in input.geos. This should only be done if turbulence is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    ! Subroutine DO_MIXING operates in mixing ratios, e.g. the tracers must 
    ! be in [v/v]. 
    !=======================================================================
!      DO i=1,220
!      WRITE(6,*) 'sunrz before doturb check','max=',maxval(State_Chm%Species(i)%Conc(1,:,:)),'min=',minval(State_Chm%Species(i)%Conc(1,:,:)),'number:',i,size(State_Chm%Species(i)%Conc(1,:,:))
!      ENDDO
    IF ( DoTurb ) THEN

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do turbulence now'
       ! Only do the following for the non-local PBL mixing
       IF ( Input_Opt%LNLPBL ) THEN

          ! Once the initial met fields have been read in, we need to find
          ! the maximum PBL level for the non-local mixing algorithm.
          ! This only has to be done once. (bmy, 5/28/20)
          !
          ! hplin: this might need change for nested domains
          ! but it appears to be just a 0-D quantity, so leave it here.
          ! do not check for FIRST, because the first call due to off-centered
          ! clock
          ! only runs emissions. Check for NCALLS instead as a hack (hplin,
          ! 4/25/22)
          !
          ! Note: NCALLS needs to account for MAX_DOM * 2 at the very least,
          ! because GIGC_Chunk_Run will be called empty for the first MAX_DOM
          ! times.
          ! Some hours of debugging were wasted here. Bug reported by axzhang
          ! from SUSTech.
          ! (hplin, 2/16/23)
          IF ( NCALLS < 17 ) THEN
             CALL Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
          ENDIF

          ! Compute the surface flux for the non-local mixing,
          ! (which means getting emissions & drydep from HEMCO)
          ! and store it in State_Chm%Surface_Flux
          CALL Compute_Sflx_For_Vdiff( Input_Opt,  State_Chm, State_Diag,    &
                                       State_Grid, State_Met, RC            )
       ENDIF

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step.
       CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag,                    &
                        State_Grid, State_Met, RC                           )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
       CALL SET_CH4 ( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
    ENDIF

    ! Check that units are correct
    !ASSERT_(GIGC_Assert_Units(am_I_Root, State_Chm))
!      DO i=1,220
!      WRITE(6,*) 'sunrz before dochem check','max=',maxval(State_Chm%Species(i)%Conc(1,:,:)),'min=',minval(State_Chm%Species(i)%Conc(1,:,:)),'number:',i,size(State_Chm%Species(i)%Conc(1,:,:))
!      ENDDO
    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do chemistry now'
       ! Calculate TOMS O3 overhead. For now, always use it from the
       ! Met field. State_Met%TO3 is imported from PCHEM.
       ! (ckeller, 10/21/2014).
       CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, &
                                    .TRUE., State_Met%TO3, RC )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( .FALSE., Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF

       ! Do chemistry
       CALL Do_Chemistry( Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Chemistry done!'
    ENDIF
    ! Check that units are correct
    !ASSERT_(GIGC_Assert_Units(am_I_Root, State_Chm))

    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do wetdep now'
       ! Do wet deposition
       CALL DO_WETDEP( Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
       !ASSERT_(RC==GC_SUCCESS)
       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
       ! write(*,*) 'after do wetdeposition, 20210925'
    ENDIF

    ! Check that units are correct
    !ASSERT_(GIGC_Assert_Units(am_I_Root, State_Chm))

    ! SDE 2017-01-06: Archive the specific humidity to allow tracers to be
    ! modified currectly after SPHU is updated
    State_Met%SPHU_prev = State_Met%SPHU
    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
    IF ( DoChem ) THEN
       CALL RECOMPUTE_OD ( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
    ENDIF
    !=======================================================================
    ! Diagnostics 
    ! In an ESMF environment, all diagnostics are passed to the MAPL 
    ! HISTORY component every time step. Thus, we can call the diagnostics 
    ! at the end of the call sequence even though we don't know yet what the
    ! next step will be (e.g. we do not know yet if this is the last time 
    ! step of this month, etc.).
    !=======================================================================
!    CALL Diagnostics_Write ( am_I_Root, Input_Opt, State_Chm, .FALSE., RC ) 
    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                    State_Grid, State_Met, RC )
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'
    !=======================================================================
    ! Clean up
    !=======================================================================

    ! Make sure tracers leave routine in v/v dry
! xiaolu
!    CALL ConvertSpc_KgKgDry_to_VVDry( am_I_Root, State_Chm, RC )
    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )

    !=======================================================================
    ! Clean up
    !=======================================================================
    ! testing only
    IF ( NCALLS < 10 ) NCALLS = NCALLS + 1

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS
  END SUBROUTINE GIGC_Chunk_Run
!*****************************************************************************



!*****************************************************************************
!SUBROUTINE GIGC_Chunk_Run_AC 
!xiaolu,2017/10/05
!*****************************************************************************




!BOP
  SUBROUTINE GCHP_PRINT_MET(I, J, L,         &
       Input_Opt, State_Grid, State_Met, LOC, RC )

    !
    ! !USES:
    !
    use state_met_mod,        only : metstate
    use input_opt_mod,        only : optinput
    use state_grid_mod,       only : grdstate

    !
    ! !INPUT PARAMETERS:
    !
    INTEGER,          INTENT(IN)    :: I         ! Grid cell lat index
    INTEGER,          INTENT(IN)    :: J         ! Grid cell lon index
    INTEGER,          INTENT(IN)    :: L         ! Grid cell lev index
    CHARACTER(LEN=*), INTENT(IN)    :: LOC       ! Call location string
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met ! Meteorology State object
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !

    !
    ! !OUTPUT PARAMETERS:
    !
    INTEGER,          INTENT(OUT)   :: RC        ! Success or failure?!
    ! !REMARKS:
    !
    ! !REVISION HISTORY:
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !
    CHARACTER(LEN=255) :: ErrorMsg, ThisLoc


    !=========================================================================
    ! GCHP_PRINT_MET begins here!
    !=========================================================================

    ErrorMsg  = ''
    ThisLoc   = ' -> at GCHP_Print_Met (in module ' // &
         'Interfaces/GCHP/gchp_chunk_mod.F)'

    ! Assume success
    RC = GC_SUCCESS

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) TRIM( LOC )
       WRITE( 6, 113 ) State_Grid%YMid(I,J), State_Grid%XMid(I,J)
    ENDIF
100 FORMAT( /, '%%%%% GCHP_PRINT_MET at ', a )
113 FORMAT( 'Lat: ', f5.1, '   Lon: ', f5.1 )
    ! Write formatted output
    IF ( Input_Opt%amIRoot ) THEN
       ! 2-D Fields
       WRITE( 6, 114 ) 'PBLH',     State_Met%PBLH(I,J),     I, J
       WRITE( 6, 114 ) 'PSC2_WET', State_Met%PSC2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PSC2_DRY', State_Met%PSC2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS1_WET',  State_Met%PS1_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS1_DRY',  State_Met%PS1_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS2_WET',  State_Met%PS2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS2_DRY',  State_Met%PS2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'TS',       State_Met%TS(I,J),       I, J
       WRITE( 6, 114 ) 'U10M',     State_Met%U10M(I,J),     I, J
       ! 3-D Fields
       WRITE( 6, 115 ) 'CLDF',     State_Met%CLDF(I,J,L),      I, J, L
       WRITE( 6, 115 ) 'OMEGA',    State_Met%OMEGA(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'PEDGE',    State_Met%PEDGE(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'T',        State_Met%T(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'U',        State_Met%U(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'V',        State_Met%V(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'AD',       State_Met%AD(I,J,L),        I, J, L
       WRITE( 6, 115 ) 'PREVSPHU', State_Met%SPHU_PREV(I,J,L), I, J, L
       WRITE( 6, 115 ) 'SPHU',     State_Met%SPHU(I,J,L),      I, J, L
       ! terminator
       WRITE( 6, 120 )
    ENDIF
114 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J  = ',2I4 )
115 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J,L= ',3I4 )
120 FORMAT( / )


  END SUBROUTINE GCHP_PRINT_MET
!EOC
END MODULE GIGC_Chunk_Mod
