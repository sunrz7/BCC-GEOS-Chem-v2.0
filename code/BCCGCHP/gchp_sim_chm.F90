!xiaolu
!---------------------------------------
!from MOZART2 mo_sim_chm.f

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/07
!---------------------------------------

      module gchp_sim_chm

      use pmgrid,         only: masterproc, plat, plon
      use units,          only: getunit
      use phys_grid,      only: scatter_field_to_chunk   !2009.04.20
      use ppgrid,         only: begchunk, endchunk, pcols       !2009.04.20

      implicit none

      private    !zf 2008.04.18
      save

      public   ::    sim_chm_data

      real ,public  :: latwts(plat)                       ! latitude weights 
      real ,allocatable ,public  :: latwtsbdy(:,:)       ! 2009.04.20

      contains
!-----------------------------------------------------------------
      subroutine sim_chm_data()
      use commap,    only: w

      implicit none

      integer   ::   iu, astat, i, j
      real      ::   latwtst(plon,plat)
!----------------------------------

      allocate( latwtsbdy(pcols,begchunk:endchunk) )
      
!      iu = getunit ()
!      open(unit = iu, file='latwts.dat',status='old')
!      read(iu,*) latwts
!      close(iu)

      latwts(:) = w(:)
 
      do i = 1, plon
      do j = 1, plat
        latwtst(i,j) = latwts(j)
      end do
      end do

      call scatter_field_to_chunk(1,1,1,plon,latwtst,latwtsbdy)   !2009.04.20

      end subroutine sim_chm_data

      end module gchp_sim_chm
