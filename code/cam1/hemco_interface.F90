!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)
!                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemco_interface
!
! !DESCRIPTION: Module HEMCO\_INTERFACE is the HEMCO-CESM interface module.
!               CESM operates on chunks thus the interface is called HCOI\_Chunk.
!               Internally it uses a gridded component to interact with the CAM
!               physics grid; these functions are internal and called HCO\_GC...
!\\
!\\
! !INTERFACE:
!
module hemco_interface
!
! !USES:

    ! Species information 
    use constituents,             only: pcnst       ! # of species
    use constituents,             only: cnst_name   ! species names
    use constituents,             only: cnst_mw     ! advected mass
    ! Grid
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs
    use hco_extra,                only: xmid, xedge, ymid, yedge, yedge_r, ysin, area_m2, ap, bp
    ! Time
    use time_manager,             only: get_curr_date,get_curr_calday
    use time_manager,             only: get_step_size, dtime, get_nstep

    use hco_extra,                only: hco_grid_init
    ! HEMCO types
    use hco_error_mod,            only: hp          ! HEMCO precision
    use hco_error_mod,            only: sp          ! HEMCO single precision used for Ptrs
    use hco_error_mod,            only: hco_success, hco_fail, hco_version
    use hco_state_mod,            only: hco_state
    use hcox_state_mod,           only: ext_state
    use hco_types_mod,            only: configobj

    use shr_kind_mod,             only: r8 => shr_kind_r8


    use hco_bcc_convert_state_mod,only: hcoi_allocate_all, cam_getbefore_hcoi, cam_regridset_hcoi
    implicit none
    private
!
! !PRIVATE MEMBER FUNCTIONS:
    public  :: HCOI_Chunk_Init
    public  :: HCOI_Chunk_Run

! !PRIVATE TYPES:
!
    character(len=256)               :: HcoRoot             ! HEMCO data root path
    character(len=256)               :: HcoConfigFile       ! HEMCO configuration file path
    character(len=256)               :: HcoDiagnFile        ! HEMCO diagnostics config file path
    type(ConfigObj), pointer         :: HcoConfig => NULL()

    type(HCO_State), pointer, public :: HcoState  => NULL()
    type(Ext_State), pointer, public :: ExtState  => NULL()

    ! HEMCO internal grid parameters
    integer                          :: HcoGridIM           ! # of lons
    integer                          :: HcoGridJM           ! # of lats

    ! HEMCO configuration parameters that are set by namelist in CESM
    integer                          :: HcoFixYY            ! if > 0, force 'Emission year'

    ! Last execution times for the HEMCO component. We are assuming that time
    ! flows unidirectionally (and forwards, for now). (hplin, 3/30/20)
    integer                          :: last_HCO_day, last_HCO_second

    ! Meteorological fields used by HEMCO to be regridded to the HEMCO grid
    ! (hplin, 3/31/20)
    ! We have to store the fields because the regridding can only take place
    ! within the GridComp.
    ! Fields are allocated after the internal grid is initialized (so my_* are
    ! avail)
    !
    ! Moved to hco_cam_convert_state_mod, 12/16/20
    real(r8), allocatable, public        :: hemco_sunrz (:,:,:,:)
    real(r8), allocatable, public        :: land_sunrz (:,:,:)
    real(r8), allocatable, public        :: lai_sunrz (:,:,:)
contains
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Chunk_Init
!
! !DESCRIPTION: HCOI\_Chunk\_Init is the initialization method for the CAM
!  interface to HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Chunk_Init()
!
! !USES:
!
        use pmgrid,            only: plon, plat, masterproc
        ! HEMCO Initialization
        use hco_config_mod,   only: config_readfile, configinit
        use hco_driver_mod,   only: hco_init
        use hco_error_mod,    only: hco_logfile_open
        use hco_logfile_mod,  only: hco_spec2log
        use hco_state_mod,    only: hcostate_init
        use hco_types_mod,    only: configobj
        use hco_types_mod,    only: listcont
        use hco_vertgrid_mod, only: hco_vertgrid_define

        ! For some overriding work in HEMCO configuration
        use hco_extlist_mod,  only: corenr
        use hco_types_mod,    only: opt, ext

        ! HEMCO extensions initialization
        use hcox_driver_mod,  only: hcox_init
!
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Init'
        integer                      :: HMRC                 ! HEMCO return code
        integer                      :: RC
        integer                      :: N                    ! Loop idx

        ! Gridded component properties.
        ! Note that while edyn_grid_comp initializes the GridComp directly using
        ! cam_instance's inst_name, I think this may cause namespace clashing.
        ! So I'll be prefixing this with hco_ just incase.
        character(len=32)            :: HCO_GC_InstName = ''
!        type(ESMF_VM)                :: hco_esmf_vm

        ! MPI stuff
        integer                      :: localPET, PETcount
        integer, allocatable         :: PETlist(:)           ! PETs for each instance of the physics grid

        ! HEMCO properties
        integer                      :: nHcoSpc

        ! Temporary string for species naming in exports
        character(len=128)           :: exportName, exportDesc
        character(len=128)           :: exportNameTmp

        ! Timing properties
        integer                      :: year, month, day, tod
        integer                      :: hour, minute, second, dt
        integer                      :: prev_day, prev_s, now_day, now_s
        integer                      :: stepsize_tmp

        ! Temporaries
        logical                      :: IsExtfrc3DEmis

        ! Temporaries for HEMCO cycling
        logical                      :: OptFound
        type(Ext), pointer           :: ThisExt
        type(Opt), pointer           :: ThisOpt

        !-----------------------------------------------------------------------


        if(masterproc) then
            write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            write(6,*) "HEMCO: Harmonized Emissions Component"
            write(6,*) "https://doi.org/10.5194/gmd-14-5487-2021 (Lin et al., 2021)"
            write(6,*) "HEMCO_CESM interface version 1.3.0"
            write(6,*) "You are using HEMCO version ", ADJUSTL(HCO_VERSION)
            !write(6,*) "ROOT: ", HcoRoot
            !write(6,*) "Config File: ", HcoConfigFile
            !write(6,*) "Diagn File: ", HcoDiagnFile
            write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        endif

        ! Assume success
        HMRC = HCO_SUCCESS


        HcoGridIM=plon
        HcoGridJM=plat
        ! Initialize the HEMCO intermediate grid
        call HCO_Grid_Init (IM_in = HcoGridIM, JM_in = HcoGridJM,  RC=RC)
        if(masterproc) then
            write(6,*) "> Initialized HEMCO Grid environment successfully!"
            write(6,*) "> Global Dimensions: ", HcoGridIM, HcoGridJM
        !    write(6,*) "> my_IM, my_JM, LM, my_CE", my_IM, my_JM, LM, my_CE
        endif

        !-----------------------------------------------------------------------
        ! Allocate HEMCO meteorological objects
        ! We are allocating globally for the whole HEMCO component here. This
        ! may
        ! clash if my_CE changes (multiple CAM instances). To be verified.
        ! Should be a easy fix regardless, simply allocate and dealloc in the
        ! run
        ! (hplin, 3/31/20)
        !-----------------------------------------------------------------------
        call HCOI_Allocate_All()
        if(masterproc) then
            write(6,*) "> Allocated HEMCO temporary met fields"
        endif
        ! We are using pcnst here, which is # of constituents.
        nHcoSpc             = pcnst          ! # of hco species?

        call ConfigInit(HcoConfig, HMRC, nModelSpecies=nHcoSpc)

        HcoConfig%amIRoot   = masterproc
        ! HcoConfig%amIRoot   = .true. ! for debug only so verbosity is higher

        HcoConfig%MetField  = 'MERRA2'
        HcoConfig%GridRes   = ''

        !-----------------------------------------------------------------------
        ! Retrieve the species list and register exports
        !-----------------------------------------------------------------------
        ! Below we directly use nHcoSpc which corresponds to the number of
        ! constituents
        ! (nHcoSpc = pcnst). Only constituents may be advected.
        HcoConfig%nModelSpc = nHcoSpc
        HcoConfig%nModelAdv = nHcoSpc            ! # of adv spc?

        do N = 1, nHcoSpc
            HcoConfig%ModelSpc(N)%ModID   = N ! model id

            HcoConfig%ModelSpc(N)%SpcName = trim(cnst_name(N)) ! only constituents can be emitted

            !----------------------------------------------
            ! Register export properties.
            !----------------------------------------------
            ! History output (this will be moved to hco_cam_exports soon
            ! hopefully)
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(N)%SpcName)

            !if(masterproc) write(iulog,*) "Exported exportName " //
            !trim(exportName) // " to history"

            ! Physics buffer
            ! Note that _AddField will prepend HCO_, so do not add it here
            !
            ! Update hplin 1/13/23: Verify if part of extfrc_lst. If yes,
            ! then allocate as 3-D. Otherwise, this field can be allocated
            ! as 2-D. This scan is somewhat intensive as it uses get_extfrc_ndx
            ! which loops through extcnt in extfrc_lst.
            !IsExtfrc3DEmis = (get_extfrc_ndx(trim(HcoConfig%ModelSpc(N)%SpcName)) .gt. 0)

            if(masterproc) write(6,*) "Setting up HEMCO exportName " // trim(exportName)
        enddo

        HcoConfigFile = 'HEMCO_Config.rc'
        HcoDiagnFile  = 'HEMCO_Diagn.rc'
        HcoRoot       = '/g12/expert3/sunrz/Extdata/HEMCO'      
 
        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 1, HMRC, IsDryRun=.false.)

        ! Open the log file
        if(masterproc) then
            call HCO_LOGFILE_OPEN(HcoConfig%Err, RC=HMRC)
        endif

        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 2, HMRC, IsDryRun=.false.)

        call HcoState_Init(HcoState, HcoConfig, nHcoSpc, HMRC)

        stepsize_tmp = get_step_size()
        if(masterproc) write(6,*) "> HEMCO_BCC: Step size is ", stepsize_tmp

        HcoState%TS_EMIS = stepsize_tmp * 1.0
        HcoState%TS_CHEM = stepsize_tmp * 1.0
        HcoState%TS_DYN  = stepsize_tmp * 1.0

        ! Not a MAPL simulation. isESMF is deceiving.
        HcoState%Options%isESMF = .false.

        ! Deposition length scale. Used for computing dry deposition frequencies
        ! over the entire PBL or the first model layer. Hardcoded for now,
        ! should load Input_Opt%PBL_DRYDEP from GEOS-Chem-CESM (hplin, 3/29/20)
        ! !FIXME
        HcoState%Options%PBL_DRYDEP = .false.

        ! Don't support DryRun option (for now)
        HcoState%Options%IsDryRun = .false.

        !-----------------------------------------------------------------------
        ! Register HEMCO species information (HEMCO state object)
        !-----------------------------------------------------------------------
        do N = 1, nHcoSpc
            HcoState%Spc(N)%ModID         = N               ! model id
            HcoState%Spc(N)%SpcName       = trim(cnst_name(N)) ! species name
            HcoState%Spc(N)%MW_g          = cnst_mw(N)     ! mol. weight [g/mol]

            ! !!! We don't set Henry's law coefficients in HEMCO_CESM !!!
            ! they are mostly used in HCOX_SeaFlux_Mod, but HCOX are unsupported
            ! (for now)
            ! (hplin, 3/29/20)

            ! If is CESM-GC, then set Henry's law constants for SeaFlux
            ! KLUDGE by hplin: 1/3/21
            ! For defined species, hard code the Henry* values for now so we can
            ! work
            ! with SeaFlux. This fragmentation will cause issues down the road,
            ! FIXME
            if(HcoState%Spc(N)%SpcName .eq. "CH3I") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 0.20265_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 3.6e+3_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "DMS") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 4.80e-1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 3100.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "ACET") then
                ! 101.325_r8. using new henry constants
                HcoState%Spc(N)%HenryK0  = 2.74e+1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5500.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "MOH") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 2.03e+2_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5600.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "ALD2") then
                ! 101.325_r8. using new henry constants
                HcoState%Spc(N)%HenryK0  = 1.30e-01_r8 * 101.325_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5900.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "MENO3") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 1.1e+1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 4700.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "ETNO3") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 1.6_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5400.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif
            ! Write to log too
            if(masterproc) then
                write(6,*) ">> Spc", N, " = ", cnst_name(N), "MW_g", cnst_mw(N)
            endif
        enddo

        HcoState%NX = plon
        HcoState%NY = plat
        HcoState%NZ = pver

        ! Pass Ap, Bp values, units [Pa], [unitless]
        ! later remove masterproc
        call HCO_VertGrid_Define(HcoState%Config,                &
                                 zGrid = HcoState%Grid%zGrid,    &
                                 nz    = HcoState%NZ,            &
                                 Ap    = Ap,                     &
                                 Bp    = Bp,                     &
                                 RC    = HMRC)

        ! Point to grid variables
        HcoState%Grid%XMID%Val         => XMid   
        HcoState%Grid%YMID%Val         => YMid   
        HcoState%Grid%XEdge%Val        => XEdge  
        HcoState%Grid%YEdge%Val        => YEdge  
        HcoState%Grid%YSin%Val         => YSin  
        HcoState%Grid%AREA_M2%Val      => AREA_M2

        ! This might not be sufficient to update the configuration within the
        ! extlist -1.
        HcoConfig%ROOT = HcoRoot

        ! DiagnFile property. This is part of the extension options and we
        ! replicate
        ! some of the functionality in order to override it here. Maybe this
        ! work could
        ! be done upstream.
        OptFound = .false.
        ThisOpt => NULL()
        ThisExt => NULL()
        ThisExt => HcoConfig%ExtList
        do while(associated(ThisExt))
            if(ThisExt%ExtNr /= CoreNr) then    ! Looking for the core extension.
                ThisExt => ThisExt%NextExt
                cycle
            endif

            ThisOpt => ThisExt%Opts
            do while(associated(ThisOpt))
                if(trim(ThisOpt%OptName) == "DiagnFile") then
                    ThisOpt%OptValue = HcoDiagnFile
                    OptFound = .true.
                    exit
                endif

                ThisOpt => ThisOpt%NextOpt
            enddo

            if(OptFound) then
                ThisExt => NULL()
            else
                ThisExt => ThisExt%NextExt
            endif
        enddo

        ThisOpt => NULL()
        ThisExt => NULL()

        call HCO_Init(HcoState, HMRC)
        if(masterproc) write(6,*) "> HEMCO initialized successfully!"


        call HCOX_Init(HcoState, ExtState, HMRC)
        if(masterproc) write(6,*) "> HEMCO extensions initialized successfully!"

     end subroutine HCOI_Chunk_Init

!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Chunk_Run
!
! !DESCRIPTION: HCOI\_Chunk\_Run is the run method for the CAM interface to
! HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Chunk_Run(phys_state, phase)
!
! !USES:
!
        use pmgrid,         only: plon, plat, masterproc
        ! HEMCO
        use hco_interface_common,      only: gethcoval, gethcodiagn
        use hco_clock_mod,             only: hcoclock_set, hcoclock_get
        use hco_clock_mod,             only: hcoclock_emissionsdone
        use hco_diagn_mod,             only: hcodiagn_autoupdate
        use hco_driver_mod,            only: hco_run
        use hco_emislist_mod,          only: hco_getptr
        use hco_calc_mod,              only: hco_evalfld
        use hco_fluxarr_mod,           only: hco_fluxarrreset
        use hco_geotools_mod,          only: hco_calcvertgrid, hco_setpblm

        use hco_state_mod,             only: hco_gethcoid
        USE hco_interface_common,      ONLY : sethcotime ! sunrz 2024/8
        ! HEMCO Extensions
        use hcox_driver_mod,           only: hcox_run
        use physics_types,  only: physics_state

        ! Necessary imported properties for physics calculations
        use hco_bcc_convert_state_mod, only: State_HCO_PSFC, State_HCO_TK, State_HCO_PBLH, State_CAM_chmDMS
        use hco_bcc_convert_state_mod, only: State_CAM_chmACET, State_CAM_chmALD2, State_CAM_chmMOH, State_CAM_chmMENO3, State_CAM_DELP_DRYs,State_CAM_chmETNO3
        use phys_grid,    only: scatter_field_to_chunk
!
! !INPUT PARAMETERS:
!



!
! !REVISION HISTORY:
!  06 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Run'

        ! Parameters
        real(r8), parameter                   :: G0_100 = 100.e+0_r8 /9.80665e+0_r8

        integer                               :: I, J, K
        integer                               :: HI, HJ, HL
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
        integer, intent(in)                :: phase               ! 1, 2
        ! Temporaries for exports
        real                                  :: TMP
        logical                               :: FND
        character(len=128)                    :: exportName, exportNameTmp
        integer                               :: spcID, N

        ! Temporaries used for export
        real                                  :: exportFldHco(plon,plat, 1:pver)
        real                                  :: exportFldCAM(pcols, pver, begchunk:endchunk)

        ! Temporaries used for export (2-D data)
        real                                  :: exportFldHco2(plon,plat)
        real                                  :: exportFldCAM2(1:pcols)

        ! For grabbing data from HEMCO Ptrs (uses HEMCO single-precision)
        real(sp), pointer                     :: Ptr2D(:,:)
        real(sp), pointer                     :: Ptr3D(:,:,:)

        ! Timing properties
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
        ! HEMCO vertical grid property pointers
        ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
        ! make PEDGE be a pointer and allocate/deallocate it on each call.
        real(hp), pointer            :: BXHEIGHT(:,:,:)    ! Grid box height [m ]
        real(hp), pointer            :: PEDGE   (:,:,:)    ! Pressure @ lvl edges [Pa]
        real(hp), pointer            :: ZSFC    (:,:  )    ! Surface geopotential [m ]


        real(hp), pointer            :: PBLM    (:,:  )    ! PBL height [m ]
        real(hp), pointer            :: PSFC    (:,:  )    ! Surface pressure [Pa]
        real(hp), pointer            :: TK      (:,:,:)    ! Temperature [K ]

        ! HEMCO return code
        integer                      :: HMRC
        logical, save                :: FIRST = .True.
        logical                      :: doExport = .False.
        integer, save                :: nCalls = 0

        ! Assume success
        HMRC = HCO_SUCCESS

        nCalls = nCalls + 1

        call CAM_GetBefore_HCOI(phys_state, phase, HcoState, ExtState)
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
        dayOfYr= get_curr_calday()  !sunrz 2024/3


        if(masterproc) then
            write(6,*) "HEMCO_BCC: Running HCOI_Chunk_Run"
            write(6,*) "HEMCO_BCC: Updating HEMCO clock to set", yr, mon, day, hour, minute, second, tod, nymd, nhms, utc, hElapsed, dayOfYr
        endif


        call HCOClock_Set(HcoState, yr, mon, day,  &
                          hour, minute, second, IsEmisTime=.true., RC=HMRC)

!        CALL HcoClock_Get( HcoState%Clock, IsEmisTime=.true., RC=HMRC )

        call HCO_FluxArrReset(HcoState, HMRC)

        !-----------------------------------------------------------------------
        ! Regrid necessary meteorological quantities (Phase 1)
        ! Computes the absolute minimum (PSFC and TK) necessary for HEMCO
        ! to define its grid.
        !-----------------------------------------------------------------------
        call CAM_RegridSet_HCOI(HcoState, ExtState, Phase=1)

        if(masterproc .and. nCalls < 10) then
            write(6,*) "HEMCO_BCC: Finished regridding CAM met fields to HEMCO (1)"
        endif

        !-----------------------------------------------------------------------
        ! Get grid properties
        !-----------------------------------------------------------------------
        ! Calculate HEMCO vertical grid properties, e.g. PEDGE,
        ! then PHIS, BXHEIGHT, T, ..., from GridEdge_Set
        !
        ! The below conversions mostly borrowed from tfritz's CESM2-GC.
        ! HEMCO CalcVertGrid can approximate all quantities. We provide them
        ! with
        ! PSFC (surface pressure), TK (temperature)
        ! in the form of allocated pointers. The rest can be inferred from Ap,
        ! Bp

        if (masterproc) then
        PSFC => State_HCO_PSFC
        TK   => State_HCO_TK


        call HCO_CalcVertGrid(HcoState, PSFC, ZSFC, TK, BXHEIGHT, PEDGE, HMRC)

        ! Pass boundary layer height to HEMCO (PBLm = PBL mixing height) [m]
        call HCO_SetPBLm(HcoState, PBLM=State_HCO_PBLH, &
                         DefVal=1000.0_hp, & ! default value
                         RC=HMRC)

        endif
        !-----------------------------------------------------------------------
        ! Regrid necessary meteorological quantities (Phase 2)
        ! Has to be below grid properties because vertical grid needs to be
        ! defined
        ! for quantities to be computed!
        !-----------------------------------------------------------------------
        call CAM_RegridSet_HCOI(HcoState, ExtState, Phase=2)
        if(masterproc .and. nCalls < 10) then
            write(6,*) "HEMCO_BCC: Finished regridding CAM met fields to HEMCO (2)"
        endif
        !-----------------------------------------------------------------------
        ! Set HEMCO options
        !-----------------------------------------------------------------------
        ! Range of species and emission categories.
        ! Set Extension number ExtNr to 0, indicating that the core
        ! module shall be executed.
        HcoState%Options%SpcMin = 1
        HcoState%Options%SpcMax = -1
        HcoState%Options%CatMin = 1
        HcoState%Options%CatMax = -1
        HcoState%Options%ExtNr  = 0

        ! Use temporary array?
        HcoState%Options%FillBuffer = .FALSE.

        !-----------------------------------------------------------------------
        ! Run HEMCO!
        !----------------------------------------------------------------------- 
        if (masterproc) then
        call HCO_Run( HcoState, 1, HMRC, IsEndStep=.false. )
        if(masterproc .and. nCalls < 10) write(6,*) "HEMCO_BCC: HCO_Run Phase 1"

        call HCO_Run( HcoState, 2, HMRC, IsEndStep=.false. )
        if(masterproc .and. nCalls < 10) write(6,*) "HEMCO_BCC: HCO_Run Phase 2"

        call HCOX_Run(HcoState, ExtState, HMRC)
        if(masterproc .and. nCalls < 10) write(6,*) "HEMCO_BCC: HCOX_Run"

        call HcoDiagn_AutoUpdate(HcoState, HMRC)
        call HcoClock_EmissionsDone(HcoState%Clock, HMRC)
!        endif !masterproc

        allocate(hemco_sunrz (plon, plat, pver, pcnst))
        hemco_sunrz(:,:,:,:) =0.0_r8    
   
        allocate(land_sunrz (plon, plat, 73))
        allocate(lai_sunrz (plon, plat, 73))
        ! For each species...
        do spcID = 1, HcoConfig%nModelSpc

            ! TODO: Eventually convert aerosol number emissions from mass fluxes
            ! directly rather than using scale factors for num_ax (1/12/23,
            ! hplin)

            ! Build history / pbuf field name (HCO_NO, HCO_CO, etc.)
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(spcID)%SpcName)
            !doExport   = (FIRST .or. associated(HcoState%Spc(spcID)%Emis%Val))
            ! if(masterproc) write(iulog,*) "HEMCO_CESM: Begin exporting " //
            ! trim(exportName)

            ! Get HEMCO emissions flux [kg/m2/s].
            ! For performance optimization ... tap into HEMCO structure directly
            ! (ugly ugly)
            ! No need to flip vertical here. The regridder will do it for us
            if(associated(HcoState%Spc(spcID)%Emis%Val) ) then
                    exportFldHco(:,:,:)   = 0.0_r8
                    exportFldCAM(:,:,:)   = 0.0_r8

                    ! Retrieve flux from HEMCO...
                    exportFldHco = HcoState%Spc(spcID)%Emis%Val
!                    call scatter_field_to_chunk(1,pver,1,plon,exportFldHco,exportFldCAM)
                ! No emission value. No need to run regrid, instead populate
                ! with zeros as needed
                ! Why not populate at top, you ask? Because zeroing out arrays
                ! is expensive, and we do not want to be doing twice the work.
            else
                    exportFldHco(:,:,:)   = 0.0_r8
            endif
        hemco_sunrz(:,:,:,spcID) = exportFldHco(:,:,:)
        enddo
        write(6,*) 'HEMCO done!'

! extra landtype and lai
            do N = 0, 72
                ! Assume success
                HMRC = HCO_SUCCESS

                ! LANDTYPExx
                exportFldCAM2(:)   = 0.0_r8

                ! Grab the pointer if available
                write(exportNameTmp, '(a,i2.2)') 'LANDTYPE', N
                call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
                if(HMRC == HCO_SUCCESS .and. FND) then
                    exportFldHco2(:,:) = Ptr2D(:,:) ! Have to promote precision
                    land_sunrz(:,:,N+1) = exportFldHco2(:,:)
!                    call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
                endif

                Ptr2D => NULL()

                ! XLAIxx
                exportFldCAM2(:)   = 0.0_r8

                ! Grab the pointer if available
                write(exportNameTmp, '(a,i2.2)') 'XLAI', N
                call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
                if(HMRC == HCO_SUCCESS .and. FND) then
                    exportFldHco2(:,:) = Ptr2D(:,:) ! Have to promote precision
                    lai_sunrz(:,:,N+1) = exportFldHco2(:,:)
!                    call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
                endif

                Ptr2D => NULL()
            enddo

        endif !masterproc
!write(6,*) 'sunrz check=',maxval(hemco_sunrz(:,:,:,4)),minval(hemco_sunrz(:,:,:,4))
    end subroutine HCOI_Chunk_Run
! EOC
end module hemco_interface
