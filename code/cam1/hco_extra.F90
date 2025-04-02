#include <misc.h>
#include <params.h>

module hco_extra
    use shr_kind_mod,             only: r8 => shr_kind_r8
    implicit none
    private
    save
    public        :: HCO_Grid_Init

! !PUBLIC TYPES:
!
    ! Global grid parameters.
    ! this may be refactored into some other structure that isn't global later.
    ! for now this will do (hplin, 2/11/20) -- and I am confident this will
    ! be the way for at least a few more years, because who touches working
    ! code? ;)
    integer, public, protected :: IM                 ! # of lons
    integer, public, protected :: JM                 ! # of lats
    integer, public, protected :: LM                 ! # of levs

    ! Computed parameters for compatibility with GEOS-Chem
    real(r8), public, protected:: DX                 ! Delta X           [deg long]
    real(r8), public, protected:: DY                 ! Delta X           [deg lat]

    ! Horizontal Coordinates
    real(r8), public, pointer  ::                  &
                                  XMid (:,:),      & ! Longitude centers [deg]
                                  XEdge(:,:),      & ! Longitude edges   [deg]
                                  YMid (:,:),      & ! Latitude  centers [deg]
                                  YEdge(:,:),      & ! Latitude  edges   [deg]
                                  YEdge_R(:,:),    & ! Latitude  edges R [rad]
                                  YSin (:,:)         ! SIN( lat edges )  [1]

    real(r8), public, pointer  ::                  &
                                  AREA_M2(:,:),    & ! Area of grid box [m^2]
                                  Ap     (:),      & ! "hyai" Hybrid-sigma Ap value [Pa]
                                  Bp     (:)         ! "hybi" Hybrid-sigma Bp value [Pa]

contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Grid_Init
!
! !DESCRIPTION: Subroutine HCO\_Grid\_Init initializes the HEMCO-CAM interface
!  grid descriptions and MPI distribution.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Grid_Init( IM_in, JM_in, RC )


        use pmgrid         
        use ppgrid,         only: pver                     ! # of levs
        ! Physical constants
        use shr_const_mod,      only: pi => shr_const_pi
        use shr_const_mod,      only: Re => shr_const_rearth
        !use hycoef
        implicit none

#include <comhyb.h>

! !INPUT PARAMETERS:
!
        integer, intent(in)         :: IM_in, JM_in            ! # lon, lat, lev global
        integer, intent(inout)      :: RC                      ! Return code
!
! !REMARKS:

! !LOCAL VARIABLES:
!
        character(len=*), parameter :: subname = 'HCO_Grid_Init'
        integer                     :: I, J, L, N
        real(r8)                    :: SIN_N, SIN_S, PI_180

        ! MPI stuff
        integer                     :: color
        integer                     :: lons_per_task, lons_overflow
        integer                     :: lats_per_task, lats_overflow
        integer                     :: lon_beg, lon_end, lat_beg, lat_end, task_cnt


        ! Some physical constants...
        PI_180 = pi / 180.0_r8

        ! Accept external dimensions.
        IM    = IM_in
        JM    = JM_in


        !-----------------------------------------------------------------------
        ! Compute vertical grid parameters
        !-----------------------------------------------------------------------
        ! Can be directly retrieved from hyai, hybi
        ! Although they need to be flipped to be passed from CAM (from tfritz)
        !
        ! Note: In GEOS-Chem, Ap, Bp are defined in hPa, 1
        !       but in HEMCO, they are defined as    Pa, 1 (see
        !       hcoi_gc_main_mod.F90 :2813)
        ! So you have to be especially wary of the units.

        ! For now, use the CAM vertical grid verbatim
        LM = pver

        ! Ap, Bp has LM+1 edges for LM levels
        allocate(Ap(LM + 1), STAT=RC)          ! LM levels, LM+1 edges
        allocate(Bp(LM + 1), STAT=RC)

        ! Allocate PEDGE information

        ! G-C def: Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
        ! CAM def: Pifce(    L) = hyai(k)*ps0 + [ hybi(k) * ps ]
        !   w.r.t. ps0 = base state srfc prs; ps = ref srfc prs.
        !
        ! Note that the vertical has to be flipped and this will need to be done
        ! everywhere else within HEMCO_CESM, too.
        do L = 1, (LM+1)
            Ap(L) = hyai(LM+2-L) * ps0
            Bp(L) = hybi(LM+2-L)
        enddo

        !-----------------------------------------------------------------------
        ! Compute horizontal grid parameters
        !-----------------------------------------------------------------------
        ! Notes: long range (i) goes from -180.0_r8 to +180.0_r8
        !        lat  range (j) goes from - 90.0_r8 to + 90.0_r8

        allocate(XMid (IM,   JM  ), STAT=RC)
        allocate(XEdge(IM+1, JM  ), STAT=RC)
        allocate(YMid (IM,   JM  ), STAT=RC)
        allocate(YEdge(IM,   JM+1), STAT=RC)
        allocate(YEdge_R(IM, JM+1), STAT=RC)
        allocate(YSin (IM,   JM+1), STAT=RC)
        allocate(AREA_M2(IM, JM  ), STAT=RC)

        ! Compute DX, DY (lon, lat)
        DX = 360.0_r8 / real(IM, r8)
        DY = 180.0_r8 / real((JM - 1), r8)

        ! Loop over horizontal grid
        ! Note: Might require special handling at poles. FIXME. (hplin, 2/11/20)
        do J = 1, JM
        do I = 1, IM
            ! Longitude centers [deg]
            XMid(I, J) = (DX * (I-1)) - 180.0_r8

            ! Latitude centers [deg]
            YMid(I, J) = (DY * (J-1)) -  90.0_r8

            ! Note half-sized polar boxes for global grid, multiply DY by 1/4 at
            ! poles
            if(J == 1) then
                YMid(I, 1)  = -90.0_r8 + (0.25_r8 * DY)
            endif
            if(J == JM) then
                YMid(I, JM) =  90.0_r8 - (0.25_r8 * DY)
            endif
            ! Edges [deg] (or called corners in CAM ionos speak)
            XEdge(I, J) = XMid(I, J) - DX * 0.5_r8
            YEdge(I, J) = YMid(I, J) - DY * 0.5_r8
            YEdge_R(I, J) = (PI_180 * YEdge(I, J))
            YSin (I, J) = SIN( YEdge_R(I, J) ) ! Needed for MAP_A2A regridding

            ! Compute the LAST edges
            if(I == IM) then
                XEdge(I+1,J) = XEdge(I, J) + DX
            endif

            ! Enforce half-sized polar boxes where northern edge of grid boxes
            ! along the SOUTH POLE to be -90 deg lat.
            if(J == 1) then
                YEdge(I, 1) = -90.0_r8
            endif

            if(J == JM) then
                ! Northern edge of grid boxes along the north pole to be +90 deg
                ! lat
                YEdge(I,J+1) = 90.0_r8

                ! Adjust for second-to-last lat edge
                YEdge(I,J  ) = YEdge(I,J+1) - (DY * 0.5_r8)
                YEdge_R(I,J) = YEdge(I,J) * PI_180
                YSin(I, J)   = SIN( YEdge_R(I, J) )

                ! Last latitude edge [radians]
                YEdge_R(I,J+1) = YEdge(I,J+1) * PI_180
                YSin(I,J+1)    = SIN( YEdge_R(I,J+1) )
            endif
        enddo
        enddo

        ! Compute grid box areas after everything is populated...
        do J = 1, JM
        do I = 1, IM
            ! Sine of latitudes at N and S edges of grid box (I,J)
           SIN_N = SIN( YEdge_R(I,J+1) )
           SIN_S = SIN( YEdge_R(I,J  ) )

           ! Grid box surface areas [m2]
           AREA_M2(I,J) = ( DX * PI_180 ) * ( Re**2 ) * (SIN_N - SIN_S)
        enddo
        enddo

        ! Output debug information on the global grid information
        ! Copied from gc_grid_mod.F90 and pressure_mod.F
        if(masterproc) then
            write( 6, '(a)' )
            write( 6, '(''%%%%%%%%%%%%%%% HEMCO GRID %%%%%%%%%%%%%%%'')' )
            write( 6, '(a)' )
            write( 6, *) 'DX', DX, 'DY', DY
            write( 6, '(''Grid box longitude centers [degrees]: '')' )
            write( 6, * ) size(XMid, 1), size(XMid, 2)
            write( 6, '(8(f8.3,1x))' ) ( XMid(I,1), I=1,IM )
            write( 6, '(a)' )
            write( 6, '(''Grid box longitude edges [degrees]: '')' )
            write( 6, * ) size(XEdge, 1), size(XEdge, 2)
            write( 6, '(8(f8.3,1x))' ) ( XEdge(I,1), I=1,IM+1 )
            write( 6, '(a)' )
            write( 6, '(''Grid box latitude centers [degrees]: '')' )
            write( 6, * ) size(YMid, 1), size(YMid, 2)
            write( 6, '(8(f8.3,1x))' ) ( YMid(1,J), J=1,JM )
            write( 6, '(a)' )
            write( 6, '(''Grid box latitude edges [degrees]: '')' )
            write( 6, * ) size(YEdge, 1), size(YEdge, 2)
            write( 6, '(8(f8.3,1x))' ) ( YEdge(1,J), J=1,JM+1 )
            write( 6, '(a)' )
            write( 6, '(''SIN( grid box latitude edges )'')' )
            write( 6, '(8(f8.3,1x))' ) ( YSin(1,J), J=1,JM+1 )

            write( 6, '(a)'   ) REPEAT( '=', 79 )
            write( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
            write( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) AP(1:LM+1)
            write( 6, '(a)'   )
            write( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) BP(1:LM+1)
            write( 6, '(a)'   ) REPEAT( '=', 79 )
        endif

    end subroutine HCO_Grid_Init
end module hco_extra
