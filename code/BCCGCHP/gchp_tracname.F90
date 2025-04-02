!XIAOLU,2017/06/01
!From mo_tracname.F90
      module gchp_tracname
!-----------------------------------------------------------
! 	... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use gc_grid,      only : pcnstm1
      use gc_chem_mods, only : grpcnt

      implicit none

      save

      character(len=8) :: tracnam(pcnstm1)          ! species names
      !-----------------
      !XIAOLU:define tracnam here, the same as cnst_name in
      !chemistry.F90

          data tracnam/'ACET','ACTA','AERI','ALD2','ALK4','AONITA', &
                 'AROMP4','AROMP5',& !8
                 'ATOOH','BALD','BCPI','BCPO','BENZ','BENZP','Br', &
                 'Br2',&!8
                 'BrCl','BrNO2','BrNO3','BrO','BrSALA','BrSALC', & 
                 'BZCO3H','BZPAN',&!8
                 'C2H2','C2H4','C2H6','C3H8','CCl4','CFC11','CFC113', & 
                 'CFC114',&!8
                 'CFC115','CFC12','CH2Br2','CH2Cl2','CH2I2','CH2IBr', & 
                 'CH2ICl','CH2O',&!8
                 'CH3Br','CH3CCl3','CH3Cl','CH3I','CH4','CHBr3', &
                 'CHCl3','Cl',&!8
                 'Cl2','Cl2O2','ClNO2','ClNO3','ClO','ClOO','CLOCK', &
                 'CO',&!8
                 'CSL','DMS','DST1','DST2','DST3','DST4','EOH','ETHLN',&!8
                 'ETHN','ETHP','ETNO3','ETP', & !4
                 'GLYC','GLYX','H1211','H1301','H2402','H2O',& !6
                 'H2O2','HAC','HBr','HC5A','HCFC123','HCFC141b', &
                 'HCFC142b','HCFC22',&!75-82
                 'HCl','HCOOH','HI','HMHP','HMML','HMS','HNO2','HNO3', &!83-90
                 'HNO4','HOBr','HOCl','HOI','HONIT','HPALD1','HPALD2', &
                 'HPALD3', &!91-98
                 'HPALD4','HPETHNL','I','I2','I2O2','I2O3','I2O4', &
                 'IBr',& !99-106
                 'ICHE','ICl','ICN','ICPDH','IDC','IDCHP','IDHDP', &
                 'IDHPE',&!107-114
                 'IDN','IEPOXA','IEPOXB','IEPOXD','IHN1','IHN2', &
                 'IHN3','IHN4',&!115-122
                 'INDIOL','INO','INPB','INPD','IO','IONITA','IONO', &
                 'IONO2',&!123-130
                 'IPRNO3','ISALA','ISALC','ISOP','ITCN','ITHN', &
                 'LIMO','LVOC',&!131-138
                 'LVOCOA','MACR','MACR1OOH','MAP','MCRDH','MCRENOL', &
                 'MCRHN','MCRHNB',&!139-146
                 'MCRHP','MCT','MEK','MENO3','MGLY','MOH','MONITA', &
                 'MONITS',&!147-154
                 'MONITU','MP','MPAN','MPN','MSA','MTPA','MTPO','MVK',& !155-162
                 'MVKDH','MVKHC','MVKHCB','MVKHP','MVKN','MVKPC',&
                 'N2O','N2O5',& !163-170
                 'NH3','NH4','NIT','NITs','NO','NO2','NO3','NPHEN',& !171-178
                 'NPRNO3','O3','OClO','OCPI','OCPO','OCS','OIO','PAN',& !179-186
                 'pFe','PHEN','PIP','PP','PPN','PROPNN','PRPE','PRPN',& !187-194
                 'PYAC','R4N2','R4P','RA3P','RB3P','RCHO','RIPA', &
                 'RIPB',&!195-202
                 'RIPC','RIPD','RP','SALA','SALAAL','SALACL','SALC',& 
                 'SALCAL',&!203-210
                 'SALCCL','SO2','SO4','SO4s','SOAGX','SOAIE','SOAP', &
                 'SOAS',&!211-218
                 'TOLU','XYLE'/   !220


      character(len=8) :: natsnam(max(1,grpcnt))  ! names of non-advected trace species



      end module gchp_tracname
