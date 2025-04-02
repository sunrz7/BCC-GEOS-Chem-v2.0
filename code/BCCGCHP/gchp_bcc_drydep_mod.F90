!xiaolu
!---------------------------------------
!Purpose: AEROSOL DRY DEPOSITION
!from bcc_drydep_mod.F90

!AUTHOR:
!XIAO LU

!REVISION HISTORY:
!FIRST VERSION,2017/06/15
!---------------------------------------
#include <misc.h>
#include <params.h>
module gchp_bcc_drydep_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use physconst,    only: gravit, rair

  public :: ddflux, d3ddflux, calcram

contains

      subroutine ddflux( ncol, vg, q, p, tv, flux )

! Compute surface flux due to dry deposition processes.


      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vg(pcols)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    p(pcols)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    tv(pcols)  ! midpoint virtual temperature in surface layer (K)

! Output arguments:

      real(r8), intent(out) ::    flux(pcols) ! flux due to dry deposition in kg/m^s/sec

! Local variables:

      integer i

      do i = 1, ncol
         flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
      end do

      return
      end subroutine ddflux

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine d3ddflux
!
! !INTERFACE:
!
   subroutine  d3ddflux ( ncol, vlc_dry, q,pmid,pdel, tv, dep_dry,dep_dry_tend,dt)
! Description:
!Do 3d- settling deposition calculations following Zender's dust codes, Dec 02.
!
! Author: Natalie Mahowald
!
      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vlc_dry(pcols,pver)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols,pver)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    pmid(pcols,pver)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    pdel(pcols,pver)   ! delta pressure across level (Pa)
      real(r8), intent(in) ::    tv(pcols,pver)  ! midpoint virtual temperature in surface layer (K)
    real(r8),            intent(in)  :: dt             ! time step

! Output arguments:

      real(r8), intent(out)   ::    dep_dry(pcols) ! flux due to dry deposition in kg /m^s/sec
      real(r8), intent(inout) ::    dep_dry_tend(pcols,pver) ! flux due to dry deposition in kg /m^s/sec

! Local variables:

      real(r8) :: flux(pcols,0:pver)  ! downward flux at each level:  kg/m2/s 
      integer i,k
      do i=1,ncol
         flux(i,0)=0._r8
      enddo
      do k=1,pver
         do i = 1, ncol
            flux(i,k) = -min(vlc_dry(i,k) * q(i,k) * pmid(i,k) /(tv(i,k) * rair), &
                      q(i,k)*pdel(i,k)/gravit/dt)
            dep_dry_tend(i,k)=  dep_dry_tend(i,k) + (flux(i,k)-flux(i,k-1))/pdel(i,k)*gravit  !kg/kg/s

         end do
      enddo
! surface flux:
      do i=1,ncol
         dep_dry(i)= flux(i,pver)
      enddo
      return
      end subroutine d3ddflux

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Calcram
!
! !INTERFACE:
!

      subroutine  calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
           ustar,ram1in,ram1,t,pmid,&
           pdel,fvin,fv)
        !
        ! !DESCRIPTION: 
        !  
        ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
        ! from Seinfeld and Pandis, p.963.
        !  
        ! Author: Natalie Mahowald
        !
        implicit none
        integer, intent(in) :: ncol
        real(r8),intent(in) :: ram1in(pcols)        !aerodynamical resistance (s/m)
        real(r8),intent(in) :: fvin(pcols)                ! sfc frc vel from land
        real(r8),intent(out) :: ram1(pcols)         !aerodynamical resistance (s/m)
        real(r8),intent(out) :: fv(pcols)                 ! sfc frc vel from land
        real(r8), intent(in) :: obklen(pcols)                 ! obklen
        real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
        real(r8), intent(in) :: landfrac(pcols)               ! land fraction
        real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
        real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
        real(r8), intent(in) :: t(pcols)       !atm temperature (K)
        real(r8), intent(in) :: pmid(pcols)    !atm pressure (Pa)
        real(r8), intent(in) :: pdel(pcols)    !atm pressure (Pa)
        real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
        real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
        real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

        ! local variables
        real(r8) :: z,psi,psi0,nu,nu0,temp,ram
        integer :: i

        !-------------------------------
        !Description: based on Seinfeld and Pandis,Chapter 19,dry deposition
        !xiaolu,2017/04/25
        !1)z,roughness length
        !
        !2)psi,psi0:dimensionless height = z/L, L is Obk length,
        !psi~[-1,1],psi0:ocean
        !
        !3)temp=z/z(ocean), or z(seaice)
        !
        !4)u0=(1-15*ps0)**0.25; u1=(1-15*ps1)**0.25

        do i=1,ncol
           z=pdel(i)*rair*t(i)/pmid(i)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995

           if(obklen(i).eq.0) then
              psi=0._r8
              psi0=0._r8
           else
              psi=min(max(z/obklen(i),-1.0_r8),1.0_r8)
              psi0=min(max(zzocen/obklen(i),-1.0_r8),1.0_r8)
           endif

           !ocean or sea ice? xiaolu
           temp=z/zzocen
           if(icefrac(i) > 0.5_r8) then 
              if(obklen(i).gt.0) then 
                 psi0=min(max(zzsice/obklen(i),-1.0_r8),1.0_r8)
              else
                 psi0=0.0_r8
              endif
              temp=z/zzsice
	   endif

           !ram for different stability regime based on psi value. xiaolu
           if(psi> 0._r8) then !stable
              ram=1/xkar/ustar(i)*(log(temp)+4.7_r8*(psi-psi0))
           else!unstable
              nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
              nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
              if(ustar(i).ne.0._r8) then
                 ram=1/xkar/ustar(i)*(log(temp) &
                      +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                      +2.0_r8*(atan(nu)-atan(nu0)))
              else
	         ram=0._r8
              endif
           endif


           if(landfrac(i) < 0.000000001_r8) then
              fv(i)=ustar(i)
              ram1(i)=ram
           else
              fv(i)=fvin(i) + ustar(i)*(1.0-landfrac(i))   ! wtw
              ram1(i)=ram1in(i) + ram*(1.0-landfrac(i))    ! wtw 
           endif
           !          write(iulog,*) i,pdel(i),t(i),pmid(i),gravit,obklen(i),psi,psi0,icefrac(i),nu,nu0,ram,ustar(i),&
           !             log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)),2.0*(atan(nu)-atan(nu0))

        enddo

        ! fvitt -- fv == 0 causes a floating point exception in 
        ! dry dep of sea salts and dust
        where ( fv(:ncol) == 0._r8 ) 
           fv(:ncol) = 1.e-12_r8
        endwhere

        return
      end subroutine calcram


!##############################################################################
end module gchp_bcc_drydep_mod
