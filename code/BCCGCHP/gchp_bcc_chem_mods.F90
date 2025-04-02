!xiaolu
!---------------------------------------
!Purpose: BCC AEROSOL CHEM MODS
!from bcc_chem_mods.F90

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/15
!---------------------------------------


#include <misc.h>
#include <params.h>

      module gchp_bcc_chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod, only : r8 => shr_kind_r8
      use constituents, only : ppcnst
      implicit none

      public :: get_rxt_ndx_bcc, rebin, flt_date, get_inv_ndx
      public :: bcc_chem_mods_inti

      save
      integer, parameter :: & 
                            hetcnt = 59,         & ! = gas_pcnst,  number of heterogeneous processes
                            phtcnt = 73,         & ! number of photolysis reactions
                            nabscol = 2            ! number of absorbing column densities
      integer, parameter :: ndst  = 4, &
                            nsst  = 2 !4->2, xiao
      integer :: imozart

      character(len=8), dimension(ndst), parameter :: &
                       dust_names  = (/'DST1', 'DST2', 'DST3', 'DST4'/)

      character(len=8), dimension(nsst), parameter :: &
               progseasalts_names = (/'SALA', 'SALC'/)


       integer, parameter :: aer_wetdep_cnt = 3
       integer            :: aer_wetdep_mapping (aer_wetdep_cnt)

       character(len=8), dimension(aer_wetdep_cnt), parameter :: & ! constituent names
          aer_wetdep_list    = (/'CB2     ','OC2     ','SO4     '/)

   integer, parameter :: rxt_tag_cnt = 142   ! 142 in WACCM

   character(len=16), dimension(rxt_tag_cnt), parameter :: &  
                    rxt_tag_lst = (/ 'jo2_a           ', 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', &
                                     'jn2o            ', 'jno             ', 'jno_i           ', 'jno2            ', &
                                     'jn2o5_a         ', 'jn2o5_b         ', 'jhno3           ', 'jno3_a          ', &
                                     'jno3_b          ', 'jho2no2_a       ', 'jho2no2_b       ', 'jch3ooh         ', &
                                     'jch2o_a         ', 'jch2o_b         ', 'jh2o_a          ', 'jh2o_b          ', &
                                     'jh2o_c          ', 'jh2o2           ', 'jcl2            ', 'jclo            ', &
                                     'joclo           ', 'jcl2o2          ', 'jhocl           ', 'jhcl            ', &
                                     'jclono2_a       ', 'jclono2_b       ', 'jbrcl           ', 'jbro            ', &
                                     'jhobr           ', 'jbrono2_a       ', 'jbrono2_b       ', 'jch3cl          ', &
                                     'jccl4           ', 'jch3ccl3        ', 'jcfcl3          ', 'jcf2cl2         ', &
                                     'jcfc113         ', 'jhcfc22         ', 'jch3br          ', 'jcf3br          ', &
                                     'jcf2clbr        ', 'jco2            ', 'jch4_a          ', 'jch4_b          ', &
                                     'jeuv_1          ', 'jeuv_2          ', 'jeuv_3          ', 'jeuv_4          ', &
                                     'jeuv_5          ', 'jeuv_6          ', 'jeuv_7          ', 'jeuv_8          ', &
                                     'jeuv_9          ', 'jeuv_10         ', 'jeuv_11         ', 'jeuv_12         ', &
                                     'jeuv_13         ', 'jeuv_14         ', 'jeuv_15         ', 'jeuv_16         ', &
                                     'jeuv_17         ', 'jeuv_18         ', 'jeuv_19         ', 'jeuv_20         ', &
                                     'jeuv_21         ', 'jeuv_22         ', 'jeuv_23         ', 'jeuv_24         ', &
                                     'jeuv_25         ', 'usr_O_O2        ', 'cph1            ', 'usr_O_O         ', &
                                     'cph18           ', 'cph19           ', 'cph20           ', 'cph21           ', &
                                     'ag2             ', 'cph22           ', 'cph23           ', 'cph24           ', &
                                     'ag1             ', 'cph17           ', 'cph16           ', 'cph29           ', &
                                     'cph25           ', 'cph26           ', 'cph27           ', 'cph28           ', &
                                     'cph8            ', 'cph12           ', 'cph13           ', 'tag_NO2_NO3     ', &
                                     'usr_N2O5_M      ', 'usr_HNO3_OH     ', 'tag_NO2_HO2     ', 'usr_HO2NO2_M    ', &
                                     'usr_CO_OH_b     ', 'cph5            ', 'cph7            ', 'cph15           ', &
                                     'cph3            ', 'cph11           ', 'cph14           ', 'cph4            ', &
                                     'cph9            ', 'usr_HO2_HO2     ', 'tag_CLO_CLO     ', 'usr_CL2O2_M     ', &
                                     'het1            ', 'het2            ', 'het3            ', 'het4            ', &
                                     'het5            ', 'het6            ', 'het7            ', 'het8            ', &
                                     'het9            ', 'het10           ', 'het11           ', 'het12           ', &
                                     'het13           ', 'het14           ', 'het15           ', 'het16           ', &
                                     'het17           ', 'ion1            ', 'ion2            ', 'ion3            ', &
                                     'ion4            ', 'ion5            ', 'ion6            ', 'ion7            ', &
                                     'ion8            ', 'ion9            ', 'ion11           ', 'elec1           ', &
                                     'elec2           ', 'elec3           ' /)

  integer, dimension(rxt_tag_cnt), parameter :: &  
      rxt_tag_map                  = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                       31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                       51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                       61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                       71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
                                       82, 83, 84, 85, 86, 87, 88, 89, 108, 109, &
                                      110, 111, 114, 115, 116, 119, 120, 122, 127, 129, &
                                      138, 139, 140, 142, 144, 145, 146, 151, 152, 153, &
                                      171, 172, 202, 203, 204, 205, 206, 207, 208, 209, &
                                      210, 211, 212, 213, 214, 215, 216, 217, 218, 219, &
                                      220, 221, 222, 223, 224, 225, 226, 227, 229, 230, &
                                      231, 232 /)

      character(len=16)  :: pht_alias_lst (phtcnt,2)
      real(r8)           :: pht_alias_mult(phtcnt,2)

contains

      subroutine bcc_chem_mods_inti
      implicit none

      pht_alias_lst(1:phtcnt,1) = &
                               (/ 'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     '/)
      pht_alias_lst(1:phtcnt,2) = &
                               (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     '/)
      pht_alias_mult(1:phtcnt,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(1:phtcnt,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)

     end subroutine bcc_chem_mods_inti


 integer function get_rxt_ndx_bcc( rxt_tag )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: rxt_tag

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_rxt_ndx_bcc = -1
    do m = 1,rxt_tag_cnt
       if( trim( rxt_tag ) == trim( rxt_tag_lst(m) ) ) then
          get_rxt_ndx_bcc = rxt_tag_map(m)
          exit
       end if
    end do

  end function get_rxt_ndx_bcc


  integer function get_inv_ndx( invariant )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use gc_chem_mods,  only : nfs, inv_lst

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: invariant

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_inv_ndx = -1
    do m = 1,nfs
       if( trim( invariant ) == trim( inv_lst(m) ) ) then
          get_inv_ndx = m
          exit
       end if
    end do

  end function get_inv_ndx

  real(r8) function flt_date( ncdate, ncsec )
    !-----------------------------------------------------------------------
    ! Purpose: Convert date and seconds of day to floating point days since
    !          0001/01/01
    !-----------------------------------------------------------------------
    use time_manager, only : timemgr_datediff
    implicit none

    !-----------------------------------------------------------------------
    !   ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)   :: ncdate      ! Current date as yyyymmdd
    integer, intent(in)   :: ncsec       ! Seconds of day for current date

    integer :: refymd = 00010101
    integer :: reftod = 0

    call timemgr_datediff(refymd, reftod, ncdate, ncsec, flt_date)

  end function flt_date

  subroutine rebin( nsrc, ntrg, src_x, trg_x, src, trg )
    !---------------------------------------------------------------
    !   ... rebin src to trg
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !   ... dummy arguments
    !---------------------------------------------------------------
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8), intent(in)      :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)      :: trg_x(ntrg+1)         ! target coordinates
    real(r8), intent(in)      :: src(nsrc)             ! source array
    real(r8), intent(out)     :: trg(ntrg)             ! target array

    !---------------------------------------------------------------
    !   ... local variables
    !---------------------------------------------------------------
    integer  :: i, l
    integer  :: si, si1
    integer  :: sil, siu
    real(r8)     :: y
    real(r8)     :: sl, su
    real(r8)     :: tl, tu

    !---------------------------------------------------------------
    !   ... check interval overlap
    !---------------------------------------------------------------
    !     if( trg_x(1) < src_x(1) .or. trg_x(ntrg+1) > src_x(nsrc+1) ) then
    !        write(iulog,*) 'rebin: target grid is outside source grid'
    !        write(iulog,*) '       target grid from ',trg_x(1),' to ',trg_x(ntrg+1)
    !        write(iulog,*) '       source grid from ',src_x(1),' to ',src_x(nsrc+1)
    !        call endrun
    !     end if

    do i = 1,ntrg
       tl = trg_x(i)
       if( tl < src_x(nsrc+1) ) then
          do sil = 1,nsrc+1
             if( tl <= src_x(sil) ) then
                exit
             end if
          end do
          tu = trg_x(i+1)
          do siu = 1,nsrc+1
             if( tu <= src_x(siu) ) then
                exit
             end if
          end do
          y   = 0._r8
          sil = max( sil,2 )
          siu = min( siu,nsrc+1 )
          do si = sil,siu
             si1 = si - 1
             sl  = max( tl,src_x(si1) )
             su  = min( tu,src_x(si) )
             y   = y + (su - sl)*src(si1)
          end do
          trg(i) = y/(trg_x(i+1) - trg_x(i))
       else
          trg(i) = 0._r8
       end if
    end do

  end subroutine rebin

end module gchp_bcc_chem_mods
