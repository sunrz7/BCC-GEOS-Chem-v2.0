!--------
!xiaolu,2016/12/11:mimic mo_chem_mods
!xiaolu,2018/03/26:need to add srfems_list
!-------- 

      module gc_chem_mods
!--------------------------------------------------------------
! ... basic chemistry array parameters
!--------------------------------------------------------------
      use gc_grid,       only : pcnstm1
      implicit none

      save

!      integer, parameter :: hetcnt = 26, & ! number of heterogeneous processes
!                            phtcnt = 33, & ! number of photo processes
!                            rxntot = 168, & ! number of total reactions
!                            gascnt = 135, & ! number of gas phase reactions
       integer, parameter ::         nfs = 4, & ! number of "fixed" species
!                            relcnt = 0, & ! number of relationship species
                             grpcnt = 3 ! number of group members
!                            imp_nzcnt = 587, & ! number of non-zero implicit matrix entries
!                            rod_nzcnt = 0, & ! number of non-zero rodas matrix entries
!                            extcnt = 3, & ! number of species with external forcing
!                            clscnt1 = 8, & ! number of species in explicit class
!                            clscnt2 = 0, & ! number of species in hov class
!                            clscnt3 = 0, & ! number of species in ebi class
!                            clscnt4 = 55, & ! number of species in implicit class

!                            clscnt5 = 0, & ! number of species in rodas class
!                            indexm = 1, & ! index of total atm density in invariant array
!                            ncol_abs = 2, & ! number of column densities
!                            indexh2o = 4, & ! index of water vapor density
!                            clsze = 4 ! loop length for implicit chemistry

      integer :: ngrp = 0
      integer :: drydep_cnt = 0
!      integer :: srfems_cnt = 0
!      integer :: rxt_alias_cnt = 0
!      integer, allocatable :: grp_mem_cnt(:)
!      integer, allocatable :: rxt_alias_map(:)
      character(len=8), allocatable :: drydep_lst(:)
!      character(len=8), allocatable :: srfems_lst(:)
      character(len=8) :: inv_lst(max(1,nfs))
      real :: adv_mass(pcnstm1)
      real :: nadv_mass(grpcnt)

!xiaolu add the list for model emission,2018/03/26
       integer, parameter :: srfems_cnt = 21
       character(len=8), dimension(srfems_cnt), parameter :: & ! constituent names
       srfems_lst    = (/'NO      ','CO      ','ALK4    ','ISOP    ','ACET    ', &
                         'MEK     ','ALD2    ','PRPE    ','CH2O    ','C2H6    ', &
                         'DMS     ','SO2     ','SO4     ','NH3     ','BCPI    ', &
                         'OCPI    ','BCPO    ','OCPO    ','CH2Br2  ','CH3Br   ', &
                         'C3H8    '/)



      end module gc_chem_mods
