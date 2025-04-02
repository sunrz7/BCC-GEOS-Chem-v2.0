#include <misc.h>
#include <params.h>
!----------
!Author: Xiao LU
!Purpose: add BCC 3-D emission (prescribed and lightning NO) to tracer concentration
!Initial version: 2018/03/10, mimic from BCC codes from Tongwen Wu
!----------
#if (defined GEOSCHEM)

MODULE gchp_bcc_em3_prod

	implicit none

	private
      public :: bcc_em3_prod


      contains
      subroutine bcc_em3_prod (  pdel,ncol,lchnk, ztodt, mmr,  pmid, tfld,                 &
				 em3data,nem3,                                             &
                                 em3_IPCC6, nem3_IPCC6,                                    &
                                 prod_no)

!--------------------------------------------------------------------
      use ppgrid,             only : pver
      use shr_kind_mod,       only : r8 => shr_kind_r8
      use ppgrid,             only : pcols
      use constituents,       only : cnst_get_ind
      use gc_grid,      only : pcnstm1
      use prescribed_em_d3,    only : idx3_NO2, idx3_SO2, idx3_CB2, idx3_SO4
      use prescribed_emis_3d,  only : idxCB1_IPCC6, idxCB2_IPCC6, idxOC1_IPCC6, idxOC2_IPCC6,   &
                                      idxSO2_IPCC6, idxNH3_IPCC6, idxCO2_IPCC6,                 &
				      idxCO_IPCC6,  idxNO_IPCC6
      use gchp_chem_utls, only : get_spc_ndx
      use physconst,     only: gravit
      use history,         only : outfld
      implicit none
!--------------------------------------------------------------------
! ... Dummy args
!--------------------------------------------------------------------
      integer,  intent(in)    :: nem3,nem3_IPCC6
      integer,  intent(in)    :: ncol
      integer,  intent(in)    :: lchnk
      real(r8), intent(in)    :: pdel(pcols,pver)
      real(r8), intent(in)    :: ztodt 
      real(r8), intent(in)    :: pmid(pcols,pver)
      real(r8), intent(in)    :: tfld(pcols,pver)
      real(r8), intent(in)    :: em3_IPCC6(pcols, pver, nem3_IPCC6)   ! mol/cm3/s
      real(r8), intent(in)    :: em3data(pcols, pver, nem3)   ! mol/cm3/s
      real(r8), intent(inout) :: mmr(pcols,pver,pcnstm1)

      !add lightning NO emission here?
      real(r8), intent(in)    :: prod_no(pcols,pver)
  
!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer :: ipcc_em3_advmass(nem3_IPCC6)
     ! data ipcc_em3_advmass /12.0,12.0,12.0,12.0,28.0,&
     !                     30.0,64.0,17.0,44.0,16.0/

      integer  :: k, i

      real(r8), parameter ::  boltz = 1.38044e-16      ! erg/k

      real(r8) :: m(pcols, pver)  ! invariant densities (molecules/cm^3)
      real(r8) :: vmr(pcols,pver)


	 !emission in vmr

      real  ::  emis_so2_vmr_d3(pcols,pver)
      real ::  emis_cb2_vmr_d3(pcols,pver)
      real ::  emis_so4_vmr_d3(pcols,pver)
      real ::  emis_co2_vmr_d3(pcols,pver)

      real ::  emis_cb1_vmr_d3(pcols,pver)
      real ::  emis_oc1_vmr_d3(pcols,pver)
      real ::  emis_oc2_vmr_d3(pcols,pver)
      real ::  emis_nh3_vmr_d3(pcols,pver)

      real ::  emis_lno_vmr_d3(pcols,pver)
      real ::  emis_no_vmr_d3(pcols,pver)
      real ::  emis_co_vmr_d3(pcols,pver)

      real ::  mbar(pcols,pver)
      real ::  sflux(pcols)
      integer  :: so2_ndx, so4_ndx, cb2_ndx, co2_ndx
      integer  :: oc1_ndx, oc2_ndx, cb1_ndx, nh3_ndx
      integer  :: co_ndx,no_ndx


!===================================================================================


ipcc_em3_advmass(1) = 12.0
ipcc_em3_advmass(2) = 12.0
ipcc_em3_advmass(3) = 12.0
ipcc_em3_advmass(4) = 12.0
ipcc_em3_advmass(5) = 28.0
ipcc_em3_advmass(6) = 30.0
ipcc_em3_advmass(7) = 64.0
ipcc_em3_advmass(8) = 17.0
ipcc_em3_advmass(9) = 44.0
ipcc_em3_advmass(10) = 16.0

	!------------
	!initialize values
	!------------
     so2_ndx = get_spc_ndx( 'SO2' )
     so4_ndx = get_spc_ndx( 'SO4' )
     cb2_ndx = get_spc_ndx( 'CB2' )
     co2_ndx = get_spc_ndx( 'CO2' )

     cb1_ndx = get_spc_ndx( 'CB1' )
     oc1_ndx = get_spc_ndx( 'OC1' )
     oc2_ndx = get_spc_ndx( 'OC2' )
     nh3_ndx = get_spc_ndx( 'NH3' )

      !different name for GEOS-Chem VS emission file
      cb1_ndx       = get_spc_ndx( 'BCPO' )
      cb2_ndx       = get_spc_ndx( 'BCPI' )
      oc1_ndx       = get_spc_ndx( 'OCPO' )
      oc2_ndx       = get_spc_ndx( 'OCPI' )

     co_ndx = get_spc_ndx( 'CO' )
     no_ndx = get_spc_ndx( 'NO' )


     emis_so2_vmr_d3(:,:) = 0.0
     emis_cb2_vmr_d3(:,:) = 0.0
     emis_so4_vmr_d3(:,:) = 0.0
     emis_cb1_vmr_d3(:,:) = 0.0
     emis_oc1_vmr_d3(:,:) = 0.0
     emis_oc2_vmr_d3(:,:) = 0.0
     emis_nh3_vmr_d3(:,:) = 0.0
     emis_co2_vmr_d3(:,:) = 0.0

     emis_lno_vmr_d3(:,:) = 0.0
     emis_no_vmr_d3(:,:) = 0.0
     emis_co_vmr_d3(:,:) = 0.0


mbar=28.966
!----------------------
!add emission to state%q (mmr)
!  for each species,
!1) mmr-> vmr
!2) add emission to vmr
!3) vmr-> mmr
!----------------------

      do k=1, pver
         m(:ncol,k) = 10. * pmid(:ncol,k) / (boltz*tfld(:ncol,k))

	 !-----
         !IPCC6 emissions
         !-----
         !SO2
         if( idxSO2_IPCC6 > 0 ) then
	     vmr(:,:) = mmr(:,:,so2_ndx) * 28.966 / ipcc_em3_advmass(idxSO2_IPCC6) 
             emis_so2_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxSO2_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_so2_vmr_d3(:ncol,k) 
             mmr(:,:,so2_ndx) = vmr(:,:) * ipcc_em3_advmass(idxSO2_IPCC6)/ 28.966
         endif

         !-----
         !CB1
         if( idxCB2_IPCC6 > 0 ) then
             vmr(:,:) = mmr(:,:,cb1_ndx) * 28.966 / ipcc_em3_advmass(idxCB1_IPCC6)
             emis_cb1_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxCB1_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_cb1_vmr_d3(:ncol,k)
             mmr(:,:,cb1_ndx) = vmr(:,:) * ipcc_em3_advmass(idxCB1_IPCC6)/ 28.966
         endif


         !-----
         !CB2
         if( idxCB2_IPCC6 > 0 ) then
             vmr(:,:) = mmr(:,:,cb2_ndx) * 28.966 / ipcc_em3_advmass(idxCB2_IPCC6)
             emis_cb2_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxCB2_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_cb2_vmr_d3(:ncol,k)
             mmr(:,:,cb2_ndx) = vmr(:,:) * ipcc_em3_advmass(idxCB2_IPCC6)/ 28.966
         endif


         !-----
         !OC1
         if( idxOC1_IPCC6 > 0 ) then
             vmr(:,:) = mmr(:,:,oc1_ndx) * 28.966 / ipcc_em3_advmass(idxOC1_IPCC6)
             emis_oc1_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxOC1_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_oc1_vmr_d3(:ncol,k)
             mmr(:,:,oc1_ndx) = vmr(:,:) * ipcc_em3_advmass(idxOC1_IPCC6)/ 28.966
         endif


         !-----
         !OC2
         if( idxOC2_IPCC6 > 0 ) then
             vmr(:,:) = mmr(:,:,oc2_ndx) * 28.966 / ipcc_em3_advmass(idxOC2_IPCC6)
             emis_oc2_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxOC2_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_oc2_vmr_d3(:ncol,k)
             mmr(:,:,oc2_ndx) = vmr(:,:) * ipcc_em3_advmass(idxOC2_IPCC6)/ 28.966
         endif


         !-----
         !NH3
         if( idxNH3_IPCC6 > 0 ) then
             vmr(:,:) = mmr(:,:,nh3_ndx) * 28.966 / ipcc_em3_advmass(idxNH3_IPCC6)
             emis_nh3_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxNH3_IPCC6) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_nh3_vmr_d3(:ncol,k)
             mmr(:,:,nh3_ndx) = vmr(:,:) * ipcc_em3_advmass(idxNH3_IPCC6)/ 28.966
         endif

         !-----
         !CO2



         !=============================
         !??????????
         !=============================
            !---------------
	    !lighting NO
            !write(*,*)'xiaolu check no_ndx in bccemd3',no_ndx
             vmr(:,:) = mmr(:,:,no_ndx) * 28.966 / ipcc_em3_advmass(idxNO_IPCC6)
             emis_lno_vmr_d3(:ncol,k) = prod_no(:ncol,k) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_lno_vmr_d3(:ncol,k)
             mmr(:,:,no_ndx) = vmr(:,:) * ipcc_em3_advmass(idxNO_IPCC6)/ 28.966
             

            !---------------
            if( idxNO_IPCC6 > 0 ) then
                vmr(:,:) = mmr(:,:,no_ndx) * 28.966 / ipcc_em3_advmass(idxNO_IPCC6)
                emis_no_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxNO_IPCC6) /m(:ncol,k) *ztodt
                !vmr(:ncol,k) = vmr(:ncol,k) + emis_no_vmr_d3(:ncol,k)*30.0/46.0 !NOx->NO
                mmr(:,:,no_ndx) = vmr(:,:) * ipcc_em3_advmass(idxNO_IPCC6)/ 28.966
            endif

            !---------------
            if( idxCO_IPCC6 > 0 ) then
                vmr(:,:) = mmr(:,:,co_ndx) * 28.966 / ipcc_em3_advmass(idxCO_IPCC6)
                emis_co_vmr_d3(:ncol,k) = em3_IPCC6(:ncol,k,idxCO_IPCC6) /m(:ncol,k) *ztodt
                vmr(:ncol,k) = vmr(:ncol,k) + emis_co_vmr_d3(:ncol,k)
                mmr(:,:,co_ndx) = vmr(:,:) * ipcc_em3_advmass(idxCO_IPCC6)/ 28.966
            endif
         !=============================

         !-----
         !USE 'old' emissions if defined
         !-----
         !SO2
         if( idx3_SO2 > 0 ) then
             vmr(:,:) = mmr(:,:,so2_ndx) * 28.966 / ipcc_em3_advmass(idxSO2_IPCC6)
             emis_so2_vmr_d3(:ncol,k) = em3data(:ncol,k,idx3_SO2) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_so2_vmr_d3(:ncol,k)
             mmr(:,:,so2_ndx) = vmr(:,:) * ipcc_em3_advmass(idxSO2_IPCC6)/ 28.966
         endif

         !-----
         !CB2
         if( idx3_CB2 > 0 ) then
             vmr(:,:) = mmr(:,:,cb2_ndx) * 28.966 / ipcc_em3_advmass(idxCB2_IPCC6)
             emis_cb2_vmr_d3(:ncol,k) = em3data(:ncol,k,idx3_CB2) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_cb2_vmr_d3(:ncol,k)
             mmr(:,:,cb2_ndx) = vmr(:,:) * ipcc_em3_advmass(idxCB2_IPCC6)/ 28.966
         endif

         !-----
         !SO4
         if( idx3_SO4 > 0 ) then
             vmr(:,:) = mmr(:,:,so4_ndx) * 28.966 /96.0 
             emis_so4_vmr_d3(:ncol,k) = em3data(:ncol,k,idx3_SO4) /m(:ncol,k) *ztodt
             vmr(:ncol,k) = vmr(:ncol,k) + emis_so4_vmr_d3(:ncol,k)
             mmr(:,:,so4_ndx) = vmr(:,:) * 96.0/ 28.966
         endif

 
      enddo

!--------------------
!output 3-D emissions
!--------------------
      !-----
      !SO2
      sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols 
            sflux(i) = sflux(i) + emis_so2_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxSO2_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_SO2 ' ,sflux , pcols, lchnk)

      !-----
      !CO
      sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols
            sflux(i) = sflux(i) + emis_co_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxCO_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_CO ' ,sflux , pcols, lchnk)

      !-----
      !NH3
      sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols
            sflux(i) = sflux(i) + emis_nh3_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxNH3_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_NH3' ,sflux , pcols, lchnk)


      !-----
      !NO
      sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols
            sflux(i) = sflux(i) + emis_no_vmr_d3(i,k)*30.0/46.0  &
             /ztodt  * ipcc_em3_advmass(idxNO_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_NO ' ,sflux , pcols, lchnk)


      !-----
      !CB2
      sflux(:) = 0.0
      do k=1, pver
         do i=1, pcols
            sflux(i) = sflux(i) + emis_cb2_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxCB2_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_CB2 ' ,sflux , pcols, lchnk)

      !-----
      !CB1
      sflux(:) = 0.0
      do k=1, pver
         do i=1, pcols
            sflux(i) = sflux(i) + emis_cb1_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxCB1_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_CB1 ' ,sflux , pcols, lchnk)

      !-----
      !OC2
      sflux(:) = 0.0
      do k=1, pver
         do i=1, pcols
            sflux(i) = sflux(i) + emis_oc2_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxOC2_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_OC2 ' ,sflux , pcols, lchnk)

      !-----
      !OC1
      sflux(:) = 0.0
      do k=1, pver
         do i=1, pcols
            sflux(i) = sflux(i) + emis_oc1_vmr_d3(i,k)/ztodt  * ipcc_em3_advmass(idxOC1_IPCC6)/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_OC1 ' ,sflux , pcols, lchnk)


      !------
      !SO4
      sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols 
            sflux(i) = sflux(i) + emis_so4_vmr_d3(i,k)/ztodt  * 96.0/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_SO4 ' ,sflux , pcols, lchnk)

      !------
      !lightning NO
      !-------
        sflux(:) = 0.0
      do k=1, pver
         do i=1,pcols
            sflux(i) = sflux(i) + emis_lno_vmr_d3(i,k)/ztodt  * 30.0/mbar(i,k) &
                              *pdel(i,k)/gravit   ! kg/m2/s
         enddo
      enddo
      call outfld('EM3D_LNO ' ,sflux , pcols, lchnk)


      end subroutine

END MODULE gchp_bcc_em3_prod
#endif
