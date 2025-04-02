!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)
!                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_cam_convert_state_mod
!
! !DESCRIPTION: Module HCO\_CAM\_Convert\_State\_Mod handles state
! conversion
!  between the CAM meteorological fields and the HEMCO state.
!
!\\
!\\
! !INTERFACE:
!
module hco_bcc_convert_state_mod

    ! HEMCO types
    use hco_error_mod,            only: sp, hp
    use hco_state_mod,            only: HCO_State
    use hcox_state_mod,           only: Ext_State

    ! MPI status in CESM
    use pmgrid,               only: masterproc,plat,plon
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs
    use shr_kind_mod,             only: r8 => shr_kind_r8
    implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public       :: HCOI_Allocate_All
    public       :: CAM_GetBefore_HCOI
    public       :: CAM_RegridSet_HCOI


! !PRIVATE TYPES:
!
    ! Flag for supported features
    logical                          :: feat_JValues

    ! Indices for reactions (rxt) and pbuf fields
    integer                          :: index_rxt_jno2, index_rxt_joh,index_JNO2, index_JOH

    ! On the CAM grid (state%psetcols, pver) (LM, my_CE)
    ! Arrays are flipped in order (k, i) for the regridder
    real(r8), pointer, public        :: State_CAM_t(:,:,:)
    real(r8), pointer, public        :: State_CAM_ps(:,:)
    real(r8), pointer, public        :: State_CAM_psdry(:,:)
    real(r8), pointer, public        :: State_CAM_pblh(:,:)

    real(r8), pointer, public        :: State_CAM_TS(:,:)
    real(r8), pointer, public        :: State_CAM_SST(:,:)
    real(r8), pointer, public        :: State_CAM_U10M(:,:)
    real(r8), pointer, public        :: State_CAM_V10M(:,:)
    real(r8), pointer, public        :: State_CAM_ALBD(:,:)
    real(r8), pointer, public        :: State_CAM_USTAR(:,:)

    real(r8), pointer, public        :: State_CAM_CSZA(:,:)

    real(r8), pointer, public        :: State_CAM_AREAM2(:,:)
    real(r8), pointer, public        :: State_CAM_AIRs  (:,:)
    real(r8), pointer, public        :: State_CAM_DELP_DRYs(:,:)

    ! Land fractions (converted to Olson) from CAM - 1D
    real(r8), pointer, public        :: State_CAM_FRLAND   (:,:)
    real(r8), pointer, public        :: State_CAM_FROCEAN  (:,:)
    real(r8), pointer, public        :: State_CAM_FRSEAICE (:,:)

    ! Chem Constituents on CAM grid
    real(r8), pointer, public        :: State_CAM_chmO3 (:,:,:)
    real(r8), pointer, public        :: State_CAM_chmNO (:,:,:)
    real(r8), pointer, public        :: State_CAM_chmNO2(:,:,:)
    real(r8), pointer, public        :: State_CAM_chmHNO3(:,:,:)

    ! For surface dep calculation, only surface needs to be copied
    real(r8), pointer, public        :: State_CAM_chmDMS(:,:)
    real(r8), pointer, public        :: State_CAM_chmACET(:,:)
    real(r8), pointer, public        :: State_CAM_chmALD2(:,:)
    real(r8), pointer, public        :: State_CAM_chmMENO3(:,:)
    real(r8), pointer, public        :: State_CAM_chmETNO3(:,:)
    real(r8), pointer, public        :: State_CAM_chmMOH(:,:)

    ! J-values from chemistry (2-D only, on surface)
    real(r8), pointer, public        :: State_CAM_JNO2(:,:)
    real(r8), pointer, public        :: State_CAM_JOH (:,:)

    ! Q at 2m [kg H2O/kg air]
    real(r8), pointer, public        :: State_CAM_QV2M(:,:)

    !------------------------------------------------------------------
    ! On the HEMCO grid (my_IM, my_JM, LM) or possibly LM+1
    ! HEMCO grid are set as POINTERs so it satisfies HEMCO which wants to point
    real(r8), pointer, public        :: State_GC_PSC2_DRY(:,:)  ! Dry pressure from PSC2_DRY
    real(r8), pointer, public        :: State_GC_DELP_DRY(:,:,:)! Delta dry pressure

    real(r8), pointer, public        :: State_HCO_AIR (:,:,:)   ! GC_AD, mass of dry air in grid box [kg]
    real(r8), pointer, public        :: State_HCO_AIRVOL(:,:,:) ! GC_AIRVOL, volume of grid box [m^3]

    real(r8), pointer, public        :: State_HCO_TK  (:,:,:)
    real(r8), pointer, public        :: State_HCO_PSFC(:,:)   ! Wet?
    real(r8), pointer, public        :: State_HCO_PBLH(:,:)   ! PBLH [m]

    real(r8), pointer, public        :: State_HCO_TS(:,:)
    real(r8), pointer, public        :: State_HCO_TSKIN(:,:)
    real(r8), pointer, public        :: State_HCO_U10M(:,:)
    real(r8), pointer, public        :: State_HCO_V10M(:,:)
    real(r8), pointer, public        :: State_HCO_ALBD(:,:)
    real(r8), pointer, public        :: State_HCO_USTAR(:,:)
    real(r8), pointer, public        :: State_HCO_F_OF_PBL(:,:,:)

    real(r8), pointer, public        :: State_HCO_CSZA(:,:)

    real(r8), pointer, public        :: State_HCO_FRLAND  (:,:)
    real(r8), pointer, public        :: State_HCO_FRLANDIC(:,:)
    real(r8), pointer, public        :: State_HCO_FROCEAN (:,:)
    real(r8), pointer, public        :: State_HCO_FRSEAICE(:,:)
    real(r8), pointer, public        :: State_HCO_FRLAKE  (:,:)

    ! Chem Constituents on HEMCO grid
    real(r8), pointer, public        :: State_HCO_chmO3 (:,:,:)
    real(r8), pointer, public        :: State_HCO_chmNO (:,:,:)
    real(r8), pointer, public        :: State_HCO_chmNO2(:,:,:)
    real(r8), pointer, public        :: State_HCO_chmHNO3(:,:,:)

    ! J-values from chemistry (passed in reverse direction)
    real(r8), pointer, public        :: State_HCO_JNO2(:,:)
    real(r8), pointer, public        :: State_HCO_JOH (:,:)

    real(r8), pointer, public        :: State_HCO_QV2M(:,:)

    ! constituent indices:
    integer                          :: id_O3, id_NO, id_NO2, id_HNO3
    integer                          :: id_H2O, id_Q
    integer                          :: id_DMS, id_ACET, id_ALD2, id_MENO3, id_ETNO3, id_MOH

contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Allocate_All
!
! !DESCRIPTION: HCOI\_Allocate\_All allocates temporary met fields for use in
!  HEMCO and performs initialization of pbuf fields for interfacing with chem.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Allocate_All()

    use hco_extra,                only: AREA_M2
    use phys_grid,    only: scatter_field_to_chunk
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Allocate_All'
        integer :: my_CE, my_IM, my_JM, LM
        integer                      :: RC

        my_CE=pcols
        my_IM=plon
        my_JM=plat
        LM=pver

        ! PBL height [m]
        ! Comes from pbuf
        allocate(State_CAM_pblh(my_CE, begchunk:endchunk), stat=RC)
        if (masterproc) then
        allocate(State_HCO_PBLH(my_IM, my_JM), stat=RC)
        endif
        ! Grid box area [m2]
        ! Used for calculation of deposition fluxes directly on CAM grid
        allocate(State_CAM_AREAM2(my_CE, begchunk:endchunk), stat=RC)

        ! Surface grid box weight [kg]
        allocate(State_CAM_AIRs(my_CE, begchunk:endchunk), stat=RC)

        ! Surface delta grid box pressure differential [hPa]
        allocate(State_CAM_DELP_DRYs(my_CE, begchunk:endchunk), stat=RC)

        ! Surface pressure (wet) [Pa]
        allocate(State_CAM_ps(my_CE, begchunk:endchunk), stat=RC)
        
        if (masterproc) then
        allocate(State_HCO_PSFC(my_IM, my_JM), stat=RC)
        endif

        ! Surface pressure (dry) [hPa]
        allocate(State_CAM_psdry(my_CE, begchunk:endchunk), stat=RC)
        if (masterproc) then
        allocate(State_GC_PSC2_DRY(my_IM, my_JM), stat=RC)
        endif
        ! T
        allocate(State_CAM_t (my_CE, LM, begchunk:endchunk), stat=RC)
        if (masterproc) then
        allocate(State_HCO_TK(my_IM, my_JM, LM), stat=RC)
        endif

        ! For HEMCO extensions...
        allocate(State_CAM_TS   (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_SST  (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_U10M (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_V10M (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_ALBD (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_USTAR(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_CSZA (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_FRLAND(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_FROCEAN(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_FRSEAICE(my_CE, begchunk:endchunk), stat=RC)

        ! QV2M
        allocate(State_CAM_QV2M(my_CE, begchunk:endchunk), stat=RC)
        if (masterproc) then
        allocate(State_HCO_QV2M(my_IM, my_JM), stat=RC)
        endif
        ! Constituents
        allocate(State_CAM_chmO3(my_CE, LM, begchunk:endchunk), stat=RC)

        allocate(State_CAM_chmNO(my_CE, LM, begchunk:endchunk), stat=RC)

        allocate(State_CAM_chmNO2(my_CE, LM, begchunk:endchunk), stat=RC)

        allocate(State_CAM_chmHNO3(my_CE, LM, begchunk:endchunk), stat=RC)

        ! Constituents for deposition flux, copy sfc only
        allocate(State_CAM_chmDMS(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_chmACET(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_chmALD2(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_chmMENO3(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_chmETNO3(my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_chmMOH(my_CE, begchunk:endchunk), stat=RC)

        ! J-values
        allocate(State_CAM_JNO2 (my_CE, begchunk:endchunk), stat=RC)
        allocate(State_CAM_JOH  (my_CE, begchunk:endchunk), stat=RC)

        ! On HEMCO grid
        if (masterproc) then
        allocate(State_HCO_AIR(my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_TS(my_IM, my_JM), stat=RC)
        allocate(State_HCO_TSKIN(my_IM, my_JM), stat=RC)
        allocate(State_HCO_U10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_V10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_ALBD(my_IM, my_JM), stat=RC)
        allocate(State_HCO_CSZA(my_IM, my_JM), stat=RC)
        allocate(State_HCO_USTAR(my_IM, my_JM), stat=RC)
        allocate(State_HCO_F_OF_PBL(my_IM, my_JM, LM), stat=RC)
        allocate(State_GC_DELP_DRY(my_IM, my_JM, LM), stat=RC)

        allocate(State_HCO_chmO3 (my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmNO (my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmNO2(my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmHNO3(my_IM, my_JM, LM), stat=RC)

        allocate(State_HCO_JOH (my_IM, my_JM), stat=RC)
        allocate(State_HCO_JNO2(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLAND(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLANDIC(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FROCEAN(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRSEAICE(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLAKE(my_IM, my_JM), stat=RC)
        endif
        ! Clear values
        State_HCO_AIR(:,:,:) = 0.0_r8
        State_HCO_PBLH(:,:) = 0.0_r8
        State_HCO_PSFC(:,:) = 0.0_r8
        State_GC_PSC2_DRY(:,:) = 0.0_r8
        State_HCO_TK(:,:,:) = 0.0_r8
        State_HCO_TS(:,:) = 0.0_r8
        State_HCO_TSKIN(:,:) = 0.0_r8
        State_HCO_U10M(:,:) = 0.0_r8
        State_HCO_V10M(:,:) = 0.0_r8
        State_HCO_ALBD(:,:) = 0.0_r8
        State_HCO_USTAR(:,:) = 0.0_r8
        State_HCO_CSZA(:,:) = 0.0_r8
        State_HCO_F_OF_PBL(:,:,:) = 0.0_r8
        State_GC_DELP_DRY(:,:,:) = 0.0_r8

        State_HCO_chmO3(:,:,:) = 0.0_r8
        State_HCO_chmNO(:,:,:) = 0.0_r8
        State_HCO_chmNO2(:,:,:) = 0.0_r8
        State_HCO_chmHNO3(:,:,:) = 0.0_r8

        State_HCO_QV2M(:,:) = 0.0_r8

        State_CAM_chmDMS(:,:) = 0.0_r8
        State_CAM_chmACET(:,:) = 0.0_r8
        State_CAM_chmALD2(:,:) = 0.0_r8
        State_CAM_chmMENO3(:,:) = 0.0_r8
        State_CAM_chmETNO3(:,:) = 0.0_r8
        State_CAM_chmMOH(:,:) = 0.0_r8

        State_CAM_JNO2(:,:) = 0.0_r8
        State_CAM_JOH (:,:) = 0.0_r8
        State_CAM_QV2M(:,:) = 0.0_r8
        State_CAM_USTAR(:,:) = 0.0_r8

        State_HCO_JNO2(:,:) = 0.0_r8
        State_HCO_JOH(:,:) = 0.0_r8

        ! Unsupported fields in CESM are directly assigned to zero.
        ! These, if ever available, will be updated along with chemistry.F90
        ! under src/chemistry/geoschem. (hplin, 1/20/23)
        State_HCO_FRLANDIC(:,:) = 0.0_r8
        State_HCO_FRLAKE(:,:) = 0.0_r8
        ! Populate persistent values
        !call HCO_Grid_HCO2CAM_2D(AREA_M2, State_CAM_AREAM2)
        call scatter_field_to_chunk (1,1,1,plon,AREA_M2,State_CAM_AREAM2)
        !if (masterproc) then
        !write(6,*) 'sunrz check area=',State_CAM_AREAM2,shape(State_CAM_AREAM2)
        !endif
    end subroutine HCOI_Allocate_All

! !INTERFACE:
!
    subroutine CAM_GetBefore_HCOI(phys_state, phase, HcoState, ExtState)
!
! !USES:
!
        ! Pbuf wrappers by hplin
!        use hco_cam_exports,only: HCO_Export_Pbuf_QueryField

        ! Type descriptors
!        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physconst,       only: gravit, cpair, tmelt, cappa, zvir, rair, rga
        ! Physics grid
        use phys_grid,      only: get_ncols_p
        use phys_grid,      only: get_rlon_all_p, get_rlat_all_p

        use ppgrid,         only: pcols                   ! max. # of columns in chunk

        ! CAM physics buffer (some fields are here and some are in phys state)
!        use physics_buffer, only: pbuf_get_chunk, pbuf_get_field
!        use physics_buffer, only: pbuf_get_index

        ! Time description and zenith angle data
!        use orbit,          only: zenith
        use time_manager,   only: get_curr_calday

        ! Constituent information to retrieve from physics state
        use constituents,   only: cnst_get_ind
        use buffer,       only: pblht
        use comsrf!,       only: srfflx_state2d
        ! Output and mpi
!        use cam_logfile,    only: iulog
!        use spmd_utils,     only: masterproc, mpicom, masterprocid, iam

! !INPUT PARAMETERS:
!
!        type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
!        type(physics_buffer_desc), pointer :: pbuf2d(:,:)
        integer, intent(in)                :: phase               ! 1, 2

        type(HCO_State), pointer           :: HcoState
        type(Ext_State), pointer           :: ExtState
!
! !REMARKS:
!  Note that arrays in CAM format stored in the HEMCO interface are (k, i) idxd
!  for the regridder
!
!  Also needs HEMCO and HEMCO extensions state information in order to check
!  whether we actually need to populate required meteorology fields
!
!  Note for CSZA:
!  - We do not have access to the state here, so we need to get geo data and
!  zenith
!    using a different method, by looping through all chunks.
!
! !REVISION HISTORY:
!  16 Dec 2020 - H.P. Lin    - Initial version
!  04 Feb 2021 - H.P. Lin    - Add CSZA calculation with geographical data
!EOP
!------------------------------------------------------------------------------
!BOC
!

        character(len=*), parameter  :: subname = 'CAM_GetBefore_HCOI'
        integer                      :: RC                   ! ESMF return code

        integer                      :: I, J, K, lchnk            ! Loop idx
        integer                      :: ncol



        ! pbuf indices:
        integer                      :: index_pblh, index_JNO2, index_JOH

        ! Temporary geographical indices, allocated to max size (pcols)
        ! need to use actual column # ncol = get_ncols_p to fill to my_CE, which
        ! is exact
        real(r8)                     :: lchnk_rlats(1:pcols), lchnk_rlons(1:pcols)
        real(r8)                     :: lchnk_zenith(1:pcols)

        ! Current calday for SZA
        real(r8)                     :: calday

        ! is this first timestep? skip reading pbuf from certain data if so
        logical, save                :: FIRST = .true.
        integer                      :: LM
        !----------------------------------------------------
        ! Assume success
        !RC = ESMF_SUCCESS
         
        LM=pver
  
        if(masterproc .and. FIRST) then
            write(6,*) "> CAM_GetBefore_HCOI entering"
        endif

        ! Get calday for cosza (current time, not midpoint of dt)
        calday = get_curr_calday()

        ! TODO: Move constituent index calculation (which only needs to be
        ! initialized once)
        ! to a place where it only runs once ...
        ! Setup constituent indices so their concentrations can be retrieved
        ! from state%q (MMR)
        call cnst_get_ind('O3', id_O3)
        call cnst_get_ind('NO', id_NO)
        call cnst_get_ind('NO2', id_NO2)
        call cnst_get_ind('HNO3', id_HNO3)

        ! Get constitutent index for specific humidity
        call cnst_get_ind('Q', id_Q)
        ! call cnst_get_ind('H2O', id_H2O)
        ! id_H2O not used for now, and also not present in CAM-chem. hplin,
        ! 9/9/22

        ! Retrieve optional - for deposition - constituent IDs
        call cnst_get_ind('DMS', id_DMS, abort=.False.)
        call cnst_get_ind('MENO3', id_MENO3, abort=.False.)
        call cnst_get_ind('ETNO3', id_ETNO3, abort=.False.)
        call cnst_get_ind('ACET', id_ACET, abort=.False.)

        ! Phase 1: Store the fields in hemco_interface (copy)
!        I = 0
        do lchnk = begchunk, endchunk    ! loop over all chunks in the physics grid
            ncol = get_ncols_p(lchnk)    ! columns per chunk
            call get_rlat_all_p(lchnk, ncol, lchnk_rlats)
            call get_rlon_all_p(lchnk, ncol, lchnk_rlons)

            ! Compute zenith for chunk
            ! (could also do it all at once but it would require a separate
            ! buffer to store ll...)
            ! FIXME hplin: this slicing might be a little inefficient
            call zenith(calday, lchnk_rlats(1:ncol), lchnk_rlons(1:ncol),lchnk_zenith(1:ncol), ncol)
            I = 0
            do J = 1, ncol               ! loop over columns in the chunk
                I = I + 1                ! advance one column
                ! 3-D Fields
                do K = 1, LM             !        chunk    col, lev
                    State_CAM_t(I,K,lchnk) = phys_state(lchnk)%t(J,K)

                    ! FIXME: hplin, check if indices actually exist!!
                    ! Chemical concentrations should be in MMR [kg/kg air]
                    State_CAM_chmO3(I,K,lchnk) = phys_state(lchnk)%q(J,K,id_O3)
                    State_CAM_chmNO(I,K,lchnk) = phys_state(lchnk)%q(J,K,id_NO)
                    State_CAM_chmNO2(I,K,lchnk) = phys_state(lchnk)%q(J,K,id_NO2)
                    State_CAM_chmHNO3(I,K,lchnk) = phys_state(lchnk)%q(J,K,id_HNO3)
                enddo

                ! Verify and retrieve surface fluxes for deposition, if
                ! available (hplin, 5/7/21)
                if(id_DMS  > 0) State_CAM_chmDMS (I,lchnk)   = phys_state(lchnk)%q(J,LM,id_DMS )
                if(id_MENO3 > 0) State_CAM_chmMENO3(I,lchnk) = phys_state(lchnk)%q(J,LM,id_MENO3)
                if(id_ETNO3 > 0) State_CAM_chmETNO3(I,lchnk) = phys_state(lchnk)%q(J,LM,id_ETNO3)
                if(id_ACET > 0) State_CAM_chmACET(I,lchnk)   = phys_state(lchnk)%q(J,LM,id_ACET)
                if(id_ALD2 > 0) State_CAM_chmALD2(I,lchnk)   = phys_state(lchnk)%q(J,LM,id_ALD2)
                if(id_MOH  > 0) State_CAM_chmMOH (I,lchnk)   = phys_state(lchnk)%q(J,LM,id_MOH )

                !----------------------------------------
                ! 2-D Fields
                !----------------------------------------

                ! DEBUG: Write to CSZA as a test for latitude to make sure we
                ! are doing correctly

                ! QV2M [kg H2O/kg air] (MMR -- this is converted to VMR by
                ! *MWdry/18 in SeaSalt)
                ! at 2M, roughly surface ~ LM (1 is TOA)
                ! (hplin, 8/10/22)
                State_CAM_QV2M(I,lchnk) = phys_state(lchnk)%q(J,LM,id_Q)

                ! PBLH [m]
                State_CAM_pblh(I,lchnk) = pblht(J,lchnk)

                ! COSZA Cosine of zenith angle [1]
                State_CAM_CSZA(I,lchnk) = lchnk_zenith(J)

                ! USTAR Friction velocity [m/s]
                State_CAM_USTAR(I,lchnk) = max(sqrt(sqrt(srfflx_state2d(lchnk)%wsx(J)**2 + srfflx_state2d(lchnk)%wsy(J)**2)*rair*phys_state(lchnk)%t(J,LM)/phys_state(lchnk)%pmid(J,LM)),0.01)

                ! Sea level pressure [Pa] (note difference in units!!)
                State_CAM_ps(I,lchnk) = phys_state(lchnk)%ps(J)
               
                ! Dry pressure [hPa] (Pa -> hPa, x0.01)
                State_CAM_psdry(I,lchnk) = phys_state(lchnk)%psdry(J) * 0.01_r8

                ! Surface temperature [K]
                if(ExtState%T2M%DoUse) then
                    State_CAM_TS(I,lchnk) = srfflx_state2d(lchnk)%ts(J)
                endif

                ! Sea surface temperature [K]
                ! DO NOT use TS - the definition in CESM-GC is wrong and will
                ! give wrong SST values
                ! which are too high, and the ocean being too warm will cause
                ! huge issues with Iodine emis.
                if(ExtState%TSKIN%DoUse) then
                    State_CAM_SST(I,lchnk) = srfflx_state2d(lchnk)%sst(J)
                endif

                ! 10M E/W and N/S wind speed [m/s] (fixme: use pver?)
                if(ExtState%U10M%DoUse) then
                    State_CAM_U10M(I,lchnk) = phys_state(lchnk)%U(J, pver)
                    State_CAM_V10M(I,lchnk) = phys_state(lchnk)%V(J, pver)
                endif

                ! Visible surface albedo [1]
                if(ExtState%ALBD%DoUse) then
                    State_CAM_ALBD(I,lchnk) = srfflx_state2d(lchnk)%asdir(J)
                endif

                ! Converted-to-Olson land fractions [1] (hplin, 8/10/22)
                if(ExtState%FRLAND%DoUse) then
                    State_CAM_FRLAND(I,lchnk) = landfrac(J,lchnk)
                endif

                ! FRLANDIC unsupported

                if(ExtState%FROCEAN%DoUse) then
                    State_CAM_FROCEAN(I,lchnk) = ocnFrac(J,lchnk) + iceFrac(J,lchnk)
                endif

                if(ExtState%FRSEAICE%DoUse) then
                    State_CAM_FRSEAICE(I,lchnk) = iceFrac(J,lchnk)
                endif

                ! FRLAKE unsupported
                ! FRSNO unsupported

            enddo
        enddo
        ! The below are stubs and are not regridded or processed for now due to
        ! lack of data
        if(FIRST) then
            FIRST = .false.

            if(masterproc) then
                write(6,*) "> CAM_GetBefore_HCOI finished"!,State_CAM_USTAR
            endif
        endif
    end subroutine CAM_GetBefore_HCOI

!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CAM_RegridSet_HCOI
!
! !DESCRIPTION: CAM\_GetBefore\_HCOI populates the internal copy of CAM state
!  within the conversion module to prepare for regridding within the gridcomp.
!\\
!\\
! !INTERFACE:
!
    subroutine CAM_RegridSet_HCOI(HcoState, ExtState, Phase)
!
! !USES:
!
        ! Type descriptors
        use physics_types,  only: physics_state

        ! Physics grid
        use phys_grid,      only: get_ncols_p
        use phys_grid,      only: scatter_field_to_chunk, gather_chunk_to_field

        ! HEMCO output container state
        USE hcox_state_mod, only: ExtDat_Set

        ! Vertical grid specification
        use hco_extra,  only: Ap, Bp, AREA_M2

        ! HEMCO utilities
        use hco_geotools_mod, only: HCO_GetSUNCOS
!
! !INPUT PARAMETERS:
        type(HCO_State), pointer           :: HcoState
        type(Ext_State), pointer           :: ExtState
        integer                            :: Phase

! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'CAM_GetBefore_HCOI'
        integer                      :: RC                   ! ESMF return code

        logical, save                :: FIRST = .true.
        integer, save                :: nCalls = 0

        integer                      :: I, J, L              ! Loop index

        ! Temporary quantities needed for PBL computation
        real(r8)                     :: BLTOP, BLTHIK, DELP
        integer                      :: LTOP

        ! Physical constants from physconstants.F90 (from GEOS-Chem)
        ! TODO: Verify consistency with CAM model! (hplin, 3/3/21)
        real(r8), parameter          :: G0_100 = 100.e+0_r8 / 9.80665e+0_r8
        real(r8), parameter          :: SCALE_HEIGHT = 7600.0_r8

        nCalls = nCalls + 1

        if(masterproc .and. nCalls < 10) then
            write(6,*) "> CAM_RegridSet_HCOI entering", phase
        endif

        !-----------------------------------------------------------------------
        ! Regrid necessary physics quantities from the CAM grid to the HEMCO
        ! grid
        ! Phase 1: Regrid
        !-----------------------------------------------------------------------
        if(Phase == 1) then
            call gather_chunk_to_field(1,1,1,plon,State_CAM_ps,State_HCO_PSFC)
            call gather_chunk_to_field(1,1,1,plon,State_CAM_psdry,State_GC_PSC2_DRY)
            call gather_chunk_to_field(1,1,1,plon,State_CAM_pblh,State_HCO_PBLH)
            call gather_chunk_to_field(1,pver,1,plon,State_CAM_t,State_HCO_TK)
            !call HCO_Grid_CAM2HCO_2D(State_CAM_ps,     State_HCO_PSFC   )
            !call HCO_Grid_CAM2HCO_2D(State_CAM_psdry,  State_GC_PSC2_DRY)
            !call HCO_Grid_CAM2HCO_2D(State_CAM_pblh,   State_HCO_PBLH   )
            !call HCO_Grid_CAM2HCO_3D(State_CAM_t,      State_HCO_TK     )

            if(masterproc .and. nCalls < 10) then
                write(6,*) "> CAM_RegridSet_HCOI exiting phase 1"
            endif


            return
        endif

        ! Below only Phase 2...
        !-----------------------------------------------------------------------
        ! Compute air quantities (hplin, 2/4/21)
        !-----------------------------------------------------------------------
        if (masterproc) then
        ! Unified loop 1 (LJI)
        do L = 1, pver
        do J = 1, plat
        do I = 1, plon
            ! Calculate DELP_DRY (from pressure_mod)
            State_GC_DELP_DRY(I,J,L) = (Ap(L)   + (Bp(L)   *State_GC_PSC2_DRY(I,J))) - &
                                       (Ap(L+1) + (Bp(L+1) *State_GC_PSC2_DRY(I,J)))

            ! Calculate AD (AIR mass)
            ! Note that AREA_M2 are GLOBAL indices so you need to perform
            ! offsetting!!
            !
            ! DELP_DRY is in [hPa]. G0_100 is 100/g, converts to [Pa], and
            ! divides by [m/s2] (accel to kg)
            State_HCO_AIR(I,J,L) = State_GC_DELP_DRY(I,J,L) * G0_100 * AREA_M2(I,J)
        enddo
        enddo
        enddo
        endif
        ! Populate CAM information equivalent
        call scatter_field_to_chunk (1,1,1,plon,State_GC_DELP_DRY,State_CAM_DELP_DRYs)
        call scatter_field_to_chunk (1,1,1,plon,State_HCO_AIR,State_CAM_AIRs)
        !call HCO_Grid_HCO2CAM_2D(State_GC_DELP_DRY(:,:,1), State_CAM_DELP_DRYs)
        !call HCO_Grid_HCO2CAM_2D(State_HCO_AIR(:,:,1), State_CAM_AIRs)

        !-----------------------------------------------------------------------
        ! Surface temperature [K] - use both for T2M and TSKIN for now according
        ! to CESM-GC,
        ! hplin 12/21/2020
        if(ExtState%T2M%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_TS, State_HCO_TS    )
            call gather_chunk_to_field (1,1,1,plon,State_CAM_TS, State_HCO_TS)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%T2M,  'T2M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_TS)
            endif
        endif

        ! Sea surface temperature [K] (hplin 3/20/23)
        if(ExtState%TSKIN%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_SST, State_HCO_TSKIN)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_SST, State_HCO_TSKIN)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_TSKIN)
            endif
        endif

        ! 10M E/W and N/S wind speed [m/s] (fixme: use pver?)
        if(ExtState%U10M%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_U10M, State_HCO_U10M)
            !call HCO_Grid_CAM2HCO_2D(State_CAM_V10M, State_HCO_V10M)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_U10M, State_HCO_U10M)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_V10M, State_HCO_V10M)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%U10M,  'U10M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_U10M)

            call ExtDat_Set(HcoState, ExtState%V10M,  'V10M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_V10M)
            endif
        endif

        ! Cos of Zenith Angle [1]
        if(ExtState%SUNCOS%DoUse) then
            ! call HCO_Grid_CAM2HCO_2D(State_CAM_CSZA, State_HCO_CSZA)
            if (masterproc) then
            ! Use native CSZA from HEMCO for consistency?
            call HCO_GetSUNCOS(HcoState, State_HCO_CSZA, 0, RC)

            call ExtDat_Set(HcoState, ExtState%SUNCOS,'SUNCOS_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_CSZA)
            endif
        endif

        ! Visible surface albedo [1]
        if(ExtState%ALBD%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_ALBD, State_HCO_ALBD)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_ALBD, State_HCO_ALBD)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%ALBD,  'ALBD_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_ALBD)
            endif
        endif

        ! Friction velocity
        if(ExtState%USTAR%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_USTAR, State_HCO_USTAR)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_USTAR, State_HCO_USTAR)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%USTAR, 'USTAR_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_USTAR)
            endif
        endif

        ! Air mass [kg]
        if(ExtState%AIR%DoUse) then
            ! This is computed above using GC routines for air quantities, so
            ! it does not necessitate a regrid from CAM.
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%AIR,   'AIRMASS_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_AIR)
            endif
        endif


        ! Constituents [MMR]
        if(ExtState%O3%DoUse) then
            !call HCO_Grid_CAM2HCO_3D(State_CAM_chmO3, State_HCO_chmO3)
            call gather_chunk_to_field (1,pver,1,plon,State_CAM_chmO3, State_HCO_chmO3)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%O3,   'HEMCO_O3_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmO3)
            endif
        endif

        if(ExtState%NO2%DoUse) then
            !call HCO_Grid_CAM2HCO_3D(State_CAM_chmNO2, State_HCO_chmNO2)
            call gather_chunk_to_field (1,pver,1,plon,State_CAM_chmNO2, State_HCO_chmNO2)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%NO2,   'HEMCO_NO2_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmNO2)
            endif
        endif

        if(ExtState%NO%DoUse) then
            !call HCO_Grid_CAM2HCO_3D(State_CAM_chmNO, State_HCO_chmNO)
            call gather_chunk_to_field (1,pver,1,plon,State_CAM_chmNO, State_HCO_chmNO)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%NO,   'HEMCO_NO_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmNO)
            endif
        endif

        ! MMR of H2O at 2m [kg H2O/kg air] (2-D only)
        if(ExtState%QV2M%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_QV2M, State_HCO_QV2M)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_QV2M, State_HCO_QV2M)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%QV2M,  'QV2M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_QV2M)
            endif
        endif

if (masterproc) then
        do J = 1, plat
            do I = 1, plon
                ! use barometric law for pressure at PBL top
                BLTOP = HcoState%Grid%PEDGE%Val(I,J,1) * EXP(-State_HCO_PBLH(I,J)/SCALE_HEIGHT)

                ! PBL thickness [hPa]
                BLTHIK = HcoState%Grid%PEDGE%Val(I,J,1) - BLTOP

                ! Now, find the PBL top level
                do L = 1, pver
                    if(BLTOP > HcoState%Grid%PEDGE%Val(I,J,L+1)) then
                        LTOP = L
                        exit
                    endif
                enddo

                ! Find the fraction of grid box (I,J,L) within the PBL
                do L = 1, pver
                    DELP = HcoState%Grid%PEDGE%Val(I,J,L) - HcoState%Grid%PEDGE%Val(I,J,L+1)
                    ! ...again, PEDGE goes up to LM+1

                    if(L < LTOP) then
                        ! grid cell lies completely below the PBL top
                        State_HCO_F_OF_PBL(I,J,L) = DELP / BLTHIK
                    else if(L == LTOP) then
                        ! grid cell straddles PBL top
                        State_HCO_F_OF_PBL(I,J,L) = (HcoState%Grid%PEDGE%Val(I,J,L) - BLTOP) / BLTHIK
                    else
                        ! grid cells lies completely above the PBL top
                        State_HCO_F_OF_PBL(I,J,L) = 0.0
                    endif

                    ! write(6,*) "I,J/L", I, J, L, ": F_OF_PBL",
                    ! State_HCO_F_OF_PBL(I,J,L)
                enddo
            enddo
        enddo
endif

        if(ExtState%FRAC_OF_PBL%DoUse) then
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_F_OF_PBL)
            endif
        endif

        ! J-values - if supported
        if(ExtState%JOH%DoUse .or. ExtState%JNO2%DoUse) then
                !call HCO_Grid_CAM2HCO_2D(State_CAM_JOH, State_HCO_JOH)
            if (masterproc) then
                call ExtDat_Set(HcoState, ExtState%JOH,  'JOH_FOR_EMIS', &
                                RC,       FIRST,         State_HCO_JOH)

                !call HCO_Grid_CAM2HCO_2D(State_CAM_JNO2, State_HCO_JNO2)

                call ExtDat_Set(HcoState, ExtState%JNO2, 'JNO2_FOR_EMIS', &
                                RC,       FIRST,         State_HCO_JNO2)
            endif
        endif

        if(ExtState%FRLAND%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_FRLAND, State_HCO_FRLAND)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_FRLAND, State_HCO_FRLAND)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FRLAND,       'FRLAND_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLAND)
            endif
        endif

        if(ExtState%FRLANDIC%DoUse) then
            ! Unsupported - State_HCO_FRLANDIC is always zero
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FRLANDIC, 'FRLANDIC_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLANDIC)
            endif
        endif

        if(ExtState%FROCEAN%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_FROCEAN, State_HCO_FROCEAN)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_FROCEAN, State_HCO_FROCEAN)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FROCEAN,      'FROCEAN_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FROCEAN)
            endif
        endif

        if(ExtState%FRSEAICE%DoUse) then
            !call HCO_Grid_CAM2HCO_2D(State_CAM_FRSEAICE, State_HCO_FRSEAICE)
            call gather_chunk_to_field (1,1,1,plon,State_CAM_FRSEAICE, State_HCO_FRSEAICE)
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FRSEAICE, 'FRSEAICE_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRSEAICE)
            endif
        endif

        if(ExtState%FRLAKE%DoUse) then
            ! Unsupported - State_HCO_FRLAKE is always zero
            if (masterproc) then
            call ExtDat_Set(HcoState, ExtState%FRLAKE,       'FRLAKE_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLAKE)
            endif
        endif

        if(FIRST) then
            FIRST = .false.
        endif

        if(masterproc .and. nCalls < 10) then
            write(6,*) "> CAM_RegridSet_HCOI finished"
        endif




    end subroutine CAM_RegridSet_HCOI


!EOC
end module hco_bcc_convert_state_mod
