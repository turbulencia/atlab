#include "tlab_error.h"

module Radiation
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: MAX_PARS, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Error
    use TLab_Memory, only: TLab_Allocate_Real
    use FDM, only: fdm_Int0
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_AirWater, only: inb_scal_e, inb_scal_qt, inb_scal_ql, inb_scal_T
    use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_InPlace
    use TLab_CacheBlock
    use OPR_Radiation
    implicit none
    private

    public :: Radiation_Initialize
    public :: Radiation_AddInfrared_Z

    ! -------------------------------------------------------------------
    integer, parameter :: ncomps_max = 10               ! maximum number of radiatively active components
    integer, parameter :: nbands_max = 10               ! maximum number of spectral bands
    type :: infrared_dt
        sequence
        integer type
        logical active(MAX_VARS), lpadding(3)           ! evolution equations affected by this term
        integer ncomps                                  ! number of radiatively active components
        integer nbands                                  ! number of spectral bands
        real(wp) :: kappa(ncomps_max, nbands_max)       ! mass absorption coefficients for each radiatively active component
        real(wp) beta(3, nbands_max)                    ! polynomial coefficients for band emission functions; assuming second-order polynomial
        real(wp) bcs_t(nbands_max)                      ! boundary conditions at the top, typically downward fluxes
        real(wp) bcs_b(nbands_max)                      ! boundary conditions at the bottom, upward fluxes or surface emissivity
    end type infrared_dt

    type(infrared_dt), public, protected :: infraredProps    ! Radiation parameters
    ! type(visible_dt), public, protected :: visibleProps    ! Radiation parameters

    integer, parameter, public :: TYPE_RAD_NONE = 0
    integer, parameter :: TYPE_IR_GRAY_LIQUID = 1
    integer, parameter :: TYPE_IR_GRAY = 2
    integer, parameter :: TYPE_IR_BAND = 3

    real(wp), allocatable, target :: bcs_ht(:)          ! flux boundary condition at the top of the domain
    real(wp), allocatable, target :: bcs_hb(:)          ! flux boundary condition at the bottom of the domain

    ! -----------------------------------------------------------------------
    ! Cache blocking information; tuned to levante in a subdomain 48x32x768
    ! integer, parameter :: groupSizeX = 1024, groupSizeY = 512, groupSizeZ = 32
    integer, parameter :: groupSizeZ = CacheBlockSizeZ

contains
    !########################################################################
    !########################################################################
    subroutine Radiation_Initialize(inifile)
        use TLab_Memory, only: imax, jmax
        use TLab_ScanFiles
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer(wi) idummy, ib, ic
        real(wp) :: dummy(MAX_PARS)

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Infrared'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#AbsorptionComponent#=<values>')
        call TLab_Write_ASCII(bakfile, '#BoundaryConditionsKmax=<values>')
        call TLab_Write_ASCII(bakfile, '#BoundaryConditionsKmin=<values>')
        call TLab_Write_ASCII(bakfile, '#BetaCoefficient#=<values>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; infraredProps%type = TYPE_RAD_NONE
        else if (trim(adjustl(sRes)) == 'grayliquid') then; infraredProps%type = TYPE_IR_GRAY_LIQUID
        else if (trim(adjustl(sRes)) == 'gray') then; infraredProps%type = TYPE_IR_GRAY
        else if (trim(adjustl(sRes)) == 'band') then; infraredProps%type = TYPE_IR_BAND
        else
            call TLab_Error(trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        infraredProps%active = .false.
        if (infraredProps%type /= TYPE_RAD_NONE) then
            infraredProps%bcs_t(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'BoundaryConditionsKmax', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, infraredProps%bcs_t)
            infraredProps%nbands = idummy

            infraredProps%bcs_b(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'BoundaryConditionsKmin', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, infraredProps%bcs_b)
            if (idummy == 1) then
                ! epsilon = infraredProps%bcs_b(idummy)
            else if (idummy /= infraredProps%nbands) then
                call TLab_Error(__FILE__//'. Error in '//trim(adjustl(block))//'.AbsorptionComponent.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            infraredProps%kappa(:, :) = 0.0_wp
            do ic = 1, ncomps_max
                write (sRes, *) ic
                call ScanFile_Char(bakfile, inifile, block, 'AbsorptionComponent'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= infraredProps%nbands) then
                        call TLab_Error(__FILE__//'. Error in '//trim(adjustl(block))//'.AbsorptionComponent.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    infraredProps%kappa(ic, 1:idummy) = dummy(1:idummy)
                else
                    exit
                end if
            end do
            infraredProps%ncomps = ic - 1           ! correct for the increment in the loop

            infraredProps%beta(:, :) = 0.0_wp
            infraredProps%beta(1:3, 1) = [2.6774e-1_wp, -1.3344e-3_wp, 1.8017e-6_wp] ! default coefficients for band 1 according to Jeevanjee, 2023 for vapor bands
            infraredProps%beta(1:3, 2) = [-2.2993e-2_wp, 8.7439e-5_wp, 1.4744e-7_wp] ! default coefficients for band 2
            do ic = 1, 3                                                        ! so far, only 3 coefficients, 2. order polynomial
                write (sRes, *) ic
                call ScanFile_Char(bakfile, inifile, block, 'BetaCoefficient'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= infraredProps%nbands - 1) then                      ! The coefficients need to add to 1, that is why we only read nbands -1
                        call TLab_Error(__FILE__//'. Error in '//trim(adjustl(block))//'.BetaCoefficient.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    infraredProps%beta(ic, 1:idummy) = dummy(1:idummy)
                end if
            end do
            infraredProps%beta(1:3, infraredProps%nbands) = [1.0_wp, 0.0_wp, 0.0_wp]                ! last band from equation sum beta_i = 1
            do ib = 1, infraredProps%nbands - 1
                infraredProps%beta(1:3, infraredProps%nbands) = infraredProps%beta(1:3, infraredProps%nbands) - infraredProps%beta(1:3, ib)
            end do

            ! -------------------------------------------------------------------
            select case (imixture)
            case (MIXT_TYPE_AIRWATER)
                infraredProps%active(inb_scal_e) = .true.               ! energy

                ! first radiatively active scalar is liquid
                ! second radiatively active scalar is vapor
                ! third radiatively active scalar is a homogeneous field, e.g., CO2

            case default
                call TLab_Error(trim(adjustl(eStr))//'Infrared only derived for airwater mixture.')
                call TLab_Stop(DNS_ERROR_OPTION)

            end select

        end if

        ! !########################################################################
        ! bakfile = trim(adjustl(inifile))//'.bak'
        ! block = 'Visible'

        !########################################################################
        ! Local allocation
        allocate (bcs_ht(imax*jmax), bcs_hb(imax*jmax))!, t_ht(imax*jmax))

        return
    end subroutine Radiation_Initialize

    !########################################################################
    !########################################################################
    subroutine Radiation_AddInfrared_Z(localProps, nx, ny, nz, s, rhs, b, tmp1, tmp2, flux_down, flux_up)
        use TLab_Arrays, only: wrk3d
        use TLab_Constants, only: BCS_MAX, BCS_MIN, BCS_BOTH
        type(infrared_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(inout) :: rhs(nx*ny*nz)            ! array to which heating rate is added
        real(wp), intent(inout) :: b(nx*ny*nz)              ! emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! flux down, local aux arrays
        real(wp), intent(inout) :: tmp2(nx*ny*nz)           ! flux up, local aux arrays
        real(wp), intent(out), optional :: flux_down(nx*ny*nz)
        real(wp), intent(out), optional :: flux_up(nx*ny*nz)

        ! -----------------------------------------------------------------------
        integer iband, ibc, nlines
        integer(wi) ib

        !########################################################################

#define a_source(i) sCache(i,1)
#define qv(i) sCache(i,2)
#define ql(i) sCache(i,3)
#define T(i) sCache(i,4)

        do ib = 1, nx*ny, groupSizeZ
            nlines = min(groupSizeZ, nx*ny - ib + 1)

            select case (localProps%type)
            case (TYPE_IR_GRAY_LIQUID)
                ! memory arrangement
                call tlab_cache_reduce_z(s(ib, inb_scal_ql), nlines, nx*ny, nz, a_source(1))

                ! absorption coefficient
                a_source(1:nlines*nz) = localProps%kappa(1, 1)*a_source(1:nlines*nz)
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_InPlace(nlines, 1, nz, rbackground, a_source(1))
                end if

                ! solve radiative transfer equation, a_source overwritten with heating rate
                ibc = BCS_MAX                                   ! Flux given at the top
                bcs_ht = localProps%bcs_t(1)                    ! downward flux at domain top
                if (localProps%bcs_b(1) > 0.0_wp) then
                    ibc = BCS_BOTH                              ! Flux also given at the bottom
                    bcs_hb = localProps%bcs_b(1)                ! upward flux at domain bottom
                end if
                if (present(flux_up)) then
                    call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, a_source(1), wrk3d, bcs_hb, bcs_ht, ibc, tmp1, tmp2)
                    ! memory arrangement
                    call tlab_cache_spread_add_z(a_source(1), nlines, nx*ny, nz, rhs(ib))
                    call tlab_cache_spread_z(tmp1, nlines, nx*ny, nz, flux_down(ib))
                    call tlab_cache_spread_z(tmp2, nlines, nx*ny, nz, flux_up(ib))
                else
                    call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, a_source(1), wrk3d, bcs_hb, bcs_ht, ibc)
                    ! memory arrangement
                    call tlab_cache_spread_add_z(a_source(1), nlines, nx*ny, nz, rhs(ib))
                end if

            case (TYPE_IR_BAND, TYPE_IR_GRAY)
                ! memory arrangement
                call tlab_cache_reduce_z(s(ib, inb_scal_ql), nlines, nx*ny, nz, ql(1))
                call tlab_cache_reduce_z(s(ib, inb_scal_qt), nlines, nx*ny, nz, qv(1))
                qv(1:nlines*nz) = qv(1:nlines*nz) - ql(1:nlines*nz)
                call tlab_cache_reduce_z(s(ib, inb_scal_T), nlines, nx*ny, nz, T(1))

                do iband = 1, localProps%nbands
                    ! emission function, Stefan-Boltzmann law
                    b(1:nlines*nz) = sigma*T(1:nlines*nz)**4* &
                                   (localProps%beta(1, iband) + T(1:nlines*nz)*(localProps%beta(2, iband) + T(1:nlines*nz)*localProps%beta(3, iband)))

                    ! absorption coefficient
                    a_source(1:nlines*nz) = localProps%kappa(1, iband)*ql(1:nlines*nz) + &
                                            localProps%kappa(2, iband)*qv(1:nlines*nz) + &
                                            localProps%kappa(3, iband)
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nlines, 1, nz, rbackground, a_source(1))
                    end if

                    ! solve radiative transfer equation, a_source overwritten with heating rate
                    bcs_ht(1:nlines) = localProps%bcs_t(iband)   ! downward flux at domain top
                    epsilon = infraredProps%bcs_b(1)            ! surface emmisivity at the bottom of the domain
                    call IR_RTE1_Global(nlines, nz, fdm_Int0, a_source(1), b, wrk3d, bcs_hb, bcs_ht, tmp1, tmp2)
                    ! call IR_RTE1_Local()
                    ! call IR_RTE1_Incremental()

                    ! memory arrangement
                    call tlab_cache_spread_add_z(a_source(1), nlines, nx*ny, nz, rhs(ib))
                    if (present(flux_up)) then
                        if (iband == 1) then
                            call tlab_cache_spread_z(tmp1, nlines, nx*ny, nz, flux_down(ib))
                            call tlab_cache_spread_z(tmp2, nlines, nx*ny, nz, flux_up(ib))
                        else
                            call tlab_cache_spread_add_z(tmp1, nlines, nx*ny, nz, flux_down(ib))
                            call tlab_cache_spread_add_z(tmp2, nlines, nx*ny, nz, flux_up(ib))
                        end if
                    end if

                end do

            end select

        end do

#undef a_source
#undef qv
#undef ql
#undef T

    end subroutine Radiation_AddInfrared_Z

end module
