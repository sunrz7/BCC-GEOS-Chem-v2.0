
        MODULE gchp_states
         use input_opt_mod                                  ! input options obj
         use state_chm_mod                                  ! chemistry stateobj
         use state_met_mod
         use state_diag_mod        !sunrz 2023/11 for v14.1.1
         use state_grid_mod        !sunrz 2023/11 
         use hco_types_mod
         use species_mod,   only : species

            TYPE(OptInput)      :: Input_Opt   ! Input Options object
            TYPE(ChmState),SAVE :: State_Chm
            TYPE(DgnState),SAVE :: State_Diag
            TYPE(GrdState),SAVE :: State_Grid
            TYPE(MetState),SAVE :: State_Met
            TYPE(ConfigObj),pointer :: HcoConfig
            TYPE(Species),  POINTER :: ThisSpc => NULL()

        END MODULE  gchp_states

 


