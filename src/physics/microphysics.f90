#include "tlab_error.h"

module Microphysics
    use TLab_Constants, only: wp, wi, pi_wp, MAX_VARS, MAX_PARS
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use TLab_Memory, only: inb_scal_array, TLab_Allocate_Real
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Warning, TLab_Error
    use TLab_CacheBlock
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_Anelastic, only: dsmooth
    use Thermo_AirWater, only: inb_scal_ql, inb_scal_qt, inb_scal_e, inb_scal_T
    use OPR_Partial, only: OPR_Partial_Z, OPR_P1
    ! use OPR_ODES
    implicit none
    private

    public :: Microphysics_Initialize
    public :: Microphysics_AddSedimentation_Z
    public :: Microphysics_Evaporation
    public :: Microphysics_Evaporation_Impl

    ! -------------------------------------------------------------------
    type microphysics_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical :: active(MAX_VARS) = .false.   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type microphysics_dt

    type(microphysics_dt), public, protected :: sedimentationProps          ! Microphysics parameters
    integer, parameter :: TYPE_SED_NONE = 0
    integer, parameter :: TYPE_SED_AIRWATER = 1
    integer, parameter :: TYPE_SED_AIRWATERSIMPLIFIED = 2

    type(microphysics_dt), public, protected :: evaporationProps            ! Microphysics parameters
    integer, parameter :: TYPE_EVA_NONE = 0
    integer, parameter, public :: TYPE_EVA_EQUILIBRIUM = 1
    integer, parameter, public :: TYPE_EVA_EXPONENTIAL = 2
    integer, parameter, public :: TYPE_EVA_EXPONENTIAL_IMPL = 3
    integer, parameter, public :: TYPE_EVA_QSCONST = 4
    integer, parameter, public :: TYPE_EVA_QSCALC = 5
    integer, parameter, public :: TYPE_EVA_QSCALC_IMPL = 6
    integer, parameter, public :: TYPE_EVA_QSCALC_SEMIIMPL = 7

    real(wp) :: settling                         ! sedimentation effects
    real(wp) :: damkohler                        ! reaction
    ! real(wp) :: stokes                           ! particle inertial effects

    real(wp), allocatable, target :: tmp_eva(:,:,:)  ! 3D temporary arrays for evaporation routine

    ! -----------------------------------------------------------------------
    ! Cache blocking information; tuned to levante in a subdomain 48x32x768
    ! integer, parameter :: groupSizeX = 1024, groupSizeY = 512, groupSizeZ = 32
    integer, parameter :: groupSizeZ = CacheBlockSizeZ

contains
!########################################################################
!########################################################################
    subroutine Microphysics_Initialize(inifile)
        use TLab_ScanFiles
        use TLab_Memory, only: imax, jmax, kmax
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy
        integer(wi) :: inb_tmp_eva = 0

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        block = 'Evaporation'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Damkohler=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; evaporationProps%type = TYPE_EVA_NONE
        elseif (trim(adjustl(sRes)) == 'equilibrium') then; evaporationProps%type = TYPE_EVA_EQUILIBRIUM
        elseif (trim(adjustl(sRes)) == 'exponential') then; evaporationProps%type = TYPE_EVA_EXPONENTIAL
        elseif (trim(adjustl(sRes)) == 'exponential-impl') then; evaporationProps%type = TYPE_EVA_EXPONENTIAL_IMPL
        elseif (trim(adjustl(sRes)) == 'qsconst') then; evaporationProps%type = TYPE_EVA_QSCONST
        elseif (trim(adjustl(sRes)) == 'qscalc') then; evaporationProps%type = TYPE_EVA_QSCALC
        elseif (trim(adjustl(sRes)) == 'qscalc-impl') then; evaporationProps%type = TYPE_EVA_QSCALC_IMPL
        elseif (trim(adjustl(sRes)) == 'qscalc-semiimpl') then; evaporationProps%type = TYPE_EVA_QSCALC_SEMIIMPL
        else
            call TLab_Error(trim(adjustl(eStr))//'Error in entry Type.'//trim(adjustl(sRes)))
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Real(bakfile, inifile, block, 'Damkohler', '0.0', damkohler)

        if (evaporationProps%type /= TYPE_EVA_NONE) then
            evaporationProps%parameters(:) = 1.0_wp         ! default values
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, evaporationProps%parameters)
            end if
            call ScanFile_Real(bakfile, inifile, block, 'Exponent', '0.0', evaporationProps%auxiliar(1))
            call ScanFile_Real(bakfile, inifile, block, 'MaxNR', '1.E-12', evaporationProps%auxiliar(2))

        end if

        ! -------------------------------------------------------------------
        ! Initialize
        if (imixture == MIXT_TYPE_AIRWATER) then

            select case (evaporationProps%type)

            case (TYPE_EVA_EQUILIBRIUM)
                inb_scal_array = inb_scal_array + 1         ! space for ql as diagnostic variable.
!                use Thermo_Base, only: dsmooth, NEWTONRAPHSON_ERROR
                dsmooth = evaporationProps%parameters(1)    ! Smooth factor for thermodynamic partial derivatives

            case (TYPE_EVA_EXPONENTIAL,TYPE_EVA_EXPONENTIAL_IMPL,TYPE_EVA_QSCONST,TYPE_EVA_QSCALC,TYPE_EVA_QSCALC_IMPL,TYPE_EVA_QSCALC_SEMIIMPL)
                evaporationProps%active(inb_scal_ql) = .true.           ! Only ql scalars is affected
                if (damkohler <= 0.0_wp) then
                    call TLab_Error(__FILE__//'. Damkohler number must be nonzero if non-equilibrium evaporation is retained.')
                    call TLab_Stop(DNS_ERROR_OPTION)
                end if
                dsmooth = evaporationProps%parameters(3)    ! Smooth factor for thermodynamic partial derivatives

            end select

            if (inb_scal_array <= 3) then
                call TLab_Error(trim(adjustl(eStr))//'AirWater mixture needs at least 4 scalar arrays.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            ! Local allocation
            select case (evaporationProps%type)

            case (TYPE_EVA_QSCALC_IMPL,TYPE_EVA_QSCALC_SEMIIMPL)
                inb_tmp_eva = 2 ! Additional memory space
                call TLab_Allocate_Real(__FILE__, tmp_eva, [imax*jmax, kmax, inb_tmp_eva], 'tmp-eva')

            end select

        end if

        ! ###################################################################
        ! Read
        block = 'Sedimentation'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Settling=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')
        call TLab_Write_ASCII(bakfile, '#Exponent=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; sedimentationProps%type = TYPE_SED_NONE
        elseif (trim(adjustl(sRes)) == 'airwater') then; sedimentationProps%type = TYPE_SED_AIRWATER
        elseif (trim(adjustl(sRes)) == 'airwatersimplified') then; sedimentationProps%type = TYPE_SED_AIRWATERSIMPLIFIED
        else
            call TLab_Error(trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Real(bakfile, inifile, block, 'Settling', '0.0', settling)

        if (sedimentationProps%type /= TYPE_SED_NONE) then
            sedimentationProps%parameters(:) = 1.0_wp           ! default values
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, sedimentationProps%parameters)
            end if

            if (imixture == MIXT_TYPE_AIRWATER) then
                sedimentationProps%active(:) = .true.           ! All scalars are affected
                call ScanFile_Real(bakfile, inifile, block, 'Exponent', '0.0', sedimentationProps%auxiliar(1))
            end if

        end if

        ! -------------------------------------------------------------------
        ! Initialize
        sedimentationProps%scalar = inb_scal_array      ! By default, sedimentation is caused by last scalar
        if (imixture == MIXT_TYPE_AIRWATER) then
            sedimentationProps%scalar = inb_scal_ql     ! sedimentation is caused by liquid
        end if

        if (sedimentationProps%type /= TYPE_SED_NONE) then
            if (settling > 0.0_wp) then
                sedimentationProps%parameters = sedimentationProps%parameters*settling ! adding the settling number in the parameter definitions
            else
                call TLab_Error(__FILE__//'. Settling number must be nonzero if sedimentation is retained.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Microphysics_Initialize

    subroutine Microphysics_AddSedimentation_Z(locProps, nx, ny, nz, is, s, rhs, tmp1, flux)
        use TLab_Arrays, only: wrk3d
        use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_StaticL
        use FDM, only: fdm_der1_Z
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(inout) :: rhs(nx*ny*nz)              ! array to which heating rate is added
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) dummy, exponent
        integer ib, nlines

        !########################################################################
        exponent = locProps%auxiliar(1)
        dummy = 1.0_wp + exponent
        if (present(flux)) flux = 0.0_wp

#define h(i) sCache(i,1)
#define qt(i) sCache(i,2)
#define ql(i) sCache(i,3)
#define T(i) sCache(i,4)
#define rhoql(i) wrk3d(i)

        do ib = 1, nx*ny, groupSizeZ
            nlines = min(groupSizeZ, nx*ny - ib + 1)

            select case (locProps%type)
            case (TYPE_SED_AIRWATER)
                ! memory arrangement
                call tlab_cache_reduce_z(s(ib, inb_scal_ql), nlines, nx*ny, nz, ql(1))
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_OutPlace(nlines, 1, nz, rbackground, ql(1), rhoql(1))
                else
                    rhoql(1:nlines*nz) = ql(1:nlines*nz)
                end if

                if (is == inb_scal_ql) then
                    if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(1.0_wp - ql(1:nlines*nz))*(rhoql(1:nlines*nz)**dummy)
                    else
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(1.0_wp - ql(1:nlines*nz))*rhoql(1:nlines*nz)
                    end if

                else if (is == inb_scal_qt) then
                    ! memory arrangement
                    call tlab_cache_reduce_z(s(ib, inb_scal_qt), nlines, nx*ny, nz, qt(1))

                    if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(1.0_wp - qt(1:nlines*nz))*(rhoql(1:nlines*nz)**dummy)
                    else
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(1.0_wp - qt(1:nlines*nz))*rhoql(1:nlines*nz)
                    end if

                else if (is == inb_scal_e) then
                    ! memory arrangement
                    call tlab_cache_reduce_z(s(ib, inb_scal_e), nlines, nx*ny, nz, h(1))
                    call tlab_cache_reduce_z(s(ib, inb_scal_qt), nlines, nx*ny, nz, qt(1))
                    call tlab_cache_reduce_z(s(ib, inb_scal_T), nlines, nx*ny, nz, T(1))

                    call Thermo_Anelastic_StaticL(nlines, 1, nz, T(1), tmp1)
                    if (exponent > 0.0_wp) then
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(tmp1(1:nlines*nz) - h(1:nlines*nz))*(rhoql(1:nlines*nz)**dummy)
                    else
                        tmp1(1:nlines*nz) = locProps%parameters(is)*(tmp1(1:nlines*nz) - h(1:nlines*nz))*rhoql(1:nlines*nz)
                    end if

                end if

                call fdm_der1_Z%compute(nlines, tmp1, wrk3d)

                call tlab_cache_spread_add_z(wrk3d, nlines, nx*ny, nz, rhs(ib))
                if (present(flux)) then
                    call tlab_cache_spread_subtract_z(tmp1, nlines, nx*ny, nz, flux(ib))
                end if

            case (TYPE_SED_AIRWATERSIMPLIFIED)
                ! to be done
                ! call OPR_Partial_Z(OPR_P1, nx, ny, nz, rhoql, tmp1)
                ! if (exponent > 0.0_wp) tmp1 = tmp1*(rhoql**exponent)
                ! source = locProps%parameters(is)*dummy*tmp1

                ! if (present(flux)) flux = -locProps%parameters(is)*(rhoql**dummy)

            end select

        end do

        return
    end subroutine Microphysics_AddSedimentation_Z

!########################################################################
!########################################################################
    subroutine Microphysics_Evaporation(locProps, nx, ny, nz, is, s, source, qsat)
        use Thermo_AirWater, only: inb_scal_T, rd_ov_rv
        use Thermo_Base, only: THERMO_PSAT, NPSAT
        use Thermo_Anelastic, only: pbackground, rbackground, Thermo_Anelastic_Weight_OutPlace
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: s(nx*ny, nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny, nz)
        real(wp), intent(inout) :: qsat(nx*ny, nz)

        target s, source

        integer(wi) :: k, ipsat

        !real(wp) :: qsat(nx*ny, nz)

        real(wp) :: dummy, mod_exponent

        real(wp), pointer :: s_active(:, :) => null()

        character(len=512) :: istr, tstr

        source = 0.0_wp

        mod_exponent = locProps%auxiliar(1)

        select case (locProps%type)
        case (TYPE_EVA_EXPONENTIAL)
            if (is == inb_scal_ql) then
                ! -------------------------------------------------------------------
                ! Calculate condensation rate as Da * q_l
                ! -------------------------------------------------------------------
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:, :, inb_scal_ql), source)
                    s_active => source
                else
                    s_active => s(:, :, inb_scal_ql)
                end if
                source = damkohler*locProps%parameters(1)*s_active

            else
                write (tstr, "(I2)") inb_scal_ql
                write (istr, "(I2)") is
                call TLab_Warning('Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

        case (TYPE_EVA_QSCONST, TYPE_EVA_QSCALC)
            if (is == inb_scal_ql) then

                select case (locProps%type)
                case (TYPE_EVA_QSCONST)
                    ! -------------------------------------------------------------------
                    ! Use a constant saturation specific humidity q_sat
                    ! -------------------------------------------------------------------
                    qsat = locProps%parameters(2)

                case (TYPE_EVA_QSCALC)
                    ! -------------------------------------------------------------------
                    ! Calculate saturation specific humidity q_sat(T,p)
                    ! -------------------------------------------------------------------
#define P_LOC pbackground
#define T_LOC s(:,:,inb_scal_T)
#define psat qsat
#define qstmp qsat
                    psat(:, :) = THERMO_PSAT(NPSAT)
                    do ipsat = NPSAT - 1, 1, -1
                        psat(:, :) = psat(:, :)*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    do k = 1, nz
                        qstmp(:,k) = rd_ov_rv/(P_LOC(k)/psat(:,k) - 1.0_wp)
                    end do
                    qsat(:,:) = qstmp(:,:)/(1.0_wp + qstmp(:,:))
#undef P_LOC
#undef T_LOC
#undef psat
#undef qstmp
                end select

                ! -------------------------------------------------------------------
                ! Calculate condensation rate as Da * ( [ 1 - (1-q_t)/(1-q_l) ]/q_s - 1 ) * q_l^(1+alpha)
                ! -------------------------------------------------------------------
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:, :, inb_scal_ql), source)
                    s_active => source
                else
                    s_active => s(:, :, inb_scal_ql)
                end if
                if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    dummy = 1.0_wp + mod_exponent
                    source = damkohler*locProps%parameters(1) * ( ( 1.0_wp - (1.0_wp-s(:,:,inb_scal_qt))/(1.0_wp-s(:,:,inb_scal_ql)))/qsat - 1.0_wp ) * sign(1.0_wp,s_active)*(abs(s_active))**dummy 
                else
                   source = damkohler*locProps%parameters(1)*((1.0_wp - (1.0_wp - s(:, :, inb_scal_qt))/(1.0_wp - s(:, :, inb_scal_ql)))/qsat - 1.0_wp)*s_active
                end if

            else
                write (tstr, "(I2)") inb_scal_ql
                write (istr, "(I2)") is
                call TLab_Warning('WARNING: Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

        end select

        nullify (s_active)

    end subroutine Microphysics_Evaporation

!########################################################################
!########################################################################
    subroutine Microphysics_Evaporation_Impl(locProps, nx, ny, nz, is, s, dtrkm, tmp1, tmp2, s_EXPL, s_active) !tmp3, tmp4)
        use Thermo_AirWater, only: inb_scal_T, rd_ov_rv, Cd, Cdv, Cvl, Lv0
        use Thermo_Base, only: THERMO_PSAT, NPSAT
        use Thermo_Anelastic, only: pbackground, rbackground, ribackground, epbackground
        use Thermo_Anelastic, only: Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_Weight_InPlace
        use Thermo_Anelastic, only: NEWTONRAPHSON_ERROR
        use TLab_Arrays, only: wrk2d
#ifdef USE_MPI
        use mpi_f08
        use TLab_WorkFlow, only: ims_err
        real(wp) diff_glob
        integer(wi) it_glob
#endif
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: dtrkm
        real(wp), intent(inout) :: s(nx*ny, nz, inb_scal_array)
        !real(wp), intent(inout), dimension(nx*ny,nz) :: tmp1, tmp2, tmp3, tmp4
        real(wp), intent(inout), dimension(nx*ny,nz) :: tmp1, tmp2, s_EXPL, s_active
        !target tmp3, tmp4


        integer(wi) :: k, ipsat, it !, ij

        !real(wp), dimension(nx*ny,nz) :: s_active
        !real(wp), pointer :: s_active(:,:) => null()
        !real(wp), pointer :: s_active(:,:) => tmp4(1:nx*ny,1:nz)
        !!!real(wp), dimension(nx*ny,nz) :: s_active_to_modexp
        
        !real(wp), dimension(nx*ny,nz) :: s_EXPL
        !real(wp), pointer :: s_EXPL(:,:) => null()
        
        !real(wp), dimension(nx*ny,nz) :: qsat
        !real(wp), dimension(nx*ny,nz) :: diqsdql
        !real(wp), dimension(nx*ny,nz) :: dTdql
            
        real(wp), pointer :: dTdql(:)   => null()
        !real(wp), pointer :: dTdql(:)   => wrk2d(1:nx*ny,1)
        real(wp), pointer :: qsat(:,:)  => null()
        !real(wp), pointer :: qsat(1:nx*ny,1:nz)    => tmp_eva(1:nx*ny*nz,1)
        real(wp), pointer :: diqsdql(:,:) => null()
        !real(wp), pointer :: diqsdql(1:nx*ny,1:nz) => tmp_eva(1:nx*ny*nz,2)

        real(wp) :: mod_exponent, diff, eps


        !target s, tmp2

        character(len=512) :: istr, tstr

        !s_active => tmp4
        !s_EXPL => tmp3

        mod_exponent = locProps%auxiliar(1)

        eps = locProps%auxiliar(2)
        !eps=1.e-12_wp

        select case (locProps%type)
        case (TYPE_EVA_EXPONENTIAL_IMPL)
            if (is == inb_scal_ql) then
                s_EXPL(:,:) = s(:,:,inb_scal_ql)
                it = 0
                diff = 1.0_wp
                do while (diff > eps .and. it <= 10)
                    ! -------------------------------------------------------------------
                    ! Calculate condensation rate as Da * q_l
                    ! -------------------------------------------------------------------
                    ! --- Derivative df(q_l)/dq_l for NRmethod ---
                    tmp1 = 1.0_wp - dtrkm*damkohler*locProps%parameters(1)

                    ! --- Source term delta f(q_l) for NRmethod ---
                    tmp2 = s(:,:,inb_scal_ql) - s_EXPL - dtrkm*damkohler*locProps%parameters(1)*s(:,:,inb_scal_ql)

                    ! -------------------------------------------------------------------
                    ! Newton-Raphson iteration
                    ! -------------------------------------------------------------------
                    s(:,:,inb_scal_ql) = s(:,:,inb_scal_ql) - tmp2/tmp1
                    s(:,:,inb_scal_ql) = min(max(s(:,:,inb_scal_ql), 0.0_wp), s(:,:,inb_scal_qt))

#ifdef _DEBUG
                    diff = maxval(abs(tmp2/tmp1/s(:, :, inb_scal_ql)), mask=merge(.True.,.False., abs(tmp1*s(:,:,inb_scal_ql)) > 0.0_wp))
#else
                    diff = maxval(abs(tmp2/tmp1/s(:, :, inb_scal_ql)))
#endif

#ifdef USE_MPI
                    call MPI_ALLREDUCE(diff, diff_glob, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
                    diff = diff_glob
#endif

                    it = it + 1
                end do
                NEWTONRAPHSON_ERROR = diff

            else
                write (tstr, "(I2)") inb_scal_ql
                write (istr, "(I2)") is
                call TLab_Warning('WARNING: Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

        case (TYPE_EVA_QSCALC_IMPL, TYPE_EVA_QSCALC_SEMIIMPL)
            dTdql   => wrk2d(1:nx*ny,1)
            qsat    => tmp_eva(:,:,1)
            diqsdql => tmp_eva(:,:,2)

            if (is == inb_scal_ql) then

                if (locProps%type == TYPE_EVA_QSCALC_SEMIIMPL) then
                    ! -------------------------------------------------------------------
                    ! Update Temperature
                    ! -------------------------------------------------------------------
                    call TLab_Diagnostic(nx, ny, nz, s)

                    ! -------------------------------------------------------------------
                    ! Update saturation specific humidity q_sat(T,p)
                    ! -------------------------------------------------------------------
#define P_LOC pbackground
#define T_LOC s(:,:,inb_scal_T)
#define psat qsat
#define qstmp qsat
                    psat(:, :) = THERMO_PSAT(NPSAT)
                    do ipsat = NPSAT - 1, 1, -1
                        psat(:, :) = psat(:, :)*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    do k = 1, nz
                        qstmp(:, k) = rd_ov_rv/(P_LOC(k)/psat(:, k) - 1.0_wp)
                    end do
                    qsat(:, :) = qstmp(:, :)/(1.0_wp + qstmp(:, :))
#undef P_LOC
#undef T_LOC
#undef psat
#undef qstmp
                end if

                diqsdql(:, :) = 0.0_wp
                s_EXPL(:,:) = s(:,:,inb_scal_ql)
                it = 0
                diff = 1.0_wp
                do while (diff > eps .and. it <= 10)
                !do while (it<=10)
                    if (locProps%type == TYPE_EVA_QSCALC_IMPL) then
                        ! -------------------------------------------------------------------
                        ! Update Temperature
                        ! -------------------------------------------------------------------
                        call TLab_Diagnostic(nx, ny, nz, s)

                        ! -------------------------------------------------------------------
                        ! Update saturation specific humidity q_sat(T,p)
                        ! -------------------------------------------------------------------
#define P_LOC pbackground
#define T_LOC s(:,:,inb_scal_T)
#define psat qsat
#define qstmp qsat

#define E_LOC epbackground
#define dpsdT diqsdql
                        psat(:, :) = THERMO_PSAT(NPSAT)
                        dpsdT(:, :) = 0.0_wp
                        do ipsat = NPSAT - 1, 1, -1
                            psat(:, :) = psat(:, :)*T_LOC + THERMO_PSAT(ipsat)
                            dpsdT(:, :) = dpsdT(:, :)*T_LOC + THERMO_PSAT(ipsat + 1)*ipsat
                        end do
                        do k = 1, nz
                            dTdql(:) = ( Lv0*(Cd+Cdv*s(:,k,inb_scal_qt)) - Cvl*(s(:,k,inb_scal_e)-E_LOC(k)) )/(( Cd + Cdv*s(:,k,inb_scal_qt) + Cvl*s(:,k,inb_scal_ql) )**2)
                            diqsdql(:,k)= -P_LOC(k)/(psat(:,k)**2)/rd_ov_rv * dpsdT(:,k) * dTdql(:)
                            qstmp(:, k) = rd_ov_rv/(P_LOC(k)/psat(:, k) - 1.0_wp)
                        end do
                        qsat(:, :) = qstmp(:, :)/(1.0_wp + qstmp(:, :))
#undef P_LOC
#undef T_LOC
#undef psat
#undef qstmp

#undef E_LOC
#undef dpsdT
                    end if

                    ! -------------------------------------------------------------------
                    ! Calculate condensation rate as Da * ( [ 1 - (1-q_t)/(1-q_l) ]/q_s - 1 ) * q_l^(1+alpha)
                    ! -------------------------------------------------------------------
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        !call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:,:,inb_scal_ql), tmp2)
                        !s_active => tmp2
                        call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:, :, inb_scal_ql), s_active)
                    else
                        !s_active => s(:,:,inb_scal_ql)
                        s_active = s(:, :, inb_scal_ql)
                    end if

                    s_active = max(s_active, 0.0_wp)

                    ! --- Derivative df(q_l)/dq_l for NRmethod ---
                    ! ( rescaled by q_l^alpha to avoid division by zero )
                    tmp1 = ((1.0_wp - (1.0_wp - s(:, :, inb_scal_qt))/(1.0_wp - s(:, :, inb_scal_ql)))/qsat - 1.0_wp)*(1.0_wp + mod_exponent)
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, tmp1)
                    end if
                    tmp1 = tmp1  + ( -(1.0_wp-s(:,:,inb_scal_qt))/(1.0_wp-s(:,:,inb_scal_ql))**2/qsat + ( 1.0_wp - (1.0_wp-s(:,:,inb_scal_qt))/(1.0_wp-s(:,:,inb_scal_ql)))*diqsdql ) * s_active
                    tmp1 = tmp1*damkohler*locProps%parameters(1)
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, ribackground, tmp1)
                    end if
                    if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        !tmp1 = sign(1.0_wp,s_active)*(abs(s_active))**(-mod_exponent) - dtrkm*tmp1
                        tmp1 = s_active**(-mod_exponent) - dtrkm*tmp1
                    else
                        tmp1 = 1.0_wp - dtrkm*tmp1
                    end if

                    ! --- Source term delta f(q_l) for NRmethod ---
                    ! ( rescaled by q_l^alpha to avoid division by zero )
                    tmp2 = damkohler*locProps%parameters(1)*((1.0_wp - (1.0_wp - s(:, :, inb_scal_qt))/(1.0_wp - s(:, :, inb_scal_ql)))/qsat - 1.0_wp)*s_active
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, ribackground, tmp2)
                    end if
                    if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        !tmp2 = sign(1.0_wp,s_active)*(abs(s_active))**(-mod_exponent) * (s(:,:,inb_scal_ql) - s_EXPL) - dtrkm*tmp2
                        tmp2 = s_active**(-mod_exponent) * (s(:,:,inb_scal_ql) - s_EXPL) - dtrkm*tmp2
                    else
                        tmp2 = s(:,:,inb_scal_ql) - s_EXPL - dtrkm*tmp2
                    end if

                    ! -------------------------------------------------------------------
                    ! Newton-Raphson iteration
                    ! -------------------------------------------------------------------
                    s(:,:,inb_scal_ql) = s(:,:,inb_scal_ql) - tmp2/tmp1
                    s(:,:,inb_scal_ql) = min(max(s(:,:,inb_scal_ql), 0.0_wp), s(:,:,inb_scal_qt))

#ifdef _DEBUG
                    diff = maxval(abs(tmp2/tmp1/s(:, :, inb_scal_ql)), mask=merge(.True.,.False., abs(tmp1*s(:,:,inb_scal_ql)) > 0.0_wp))
#else
                    diff = maxval(abs(tmp2/tmp1/s(:, :, inb_scal_ql)))
#endif

#ifdef USE_MPI
                    call MPI_ALLREDUCE(diff, diff_glob, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
                    diff = diff_glob
#endif

                    it = it + 1
                end do
                NEWTONRAPHSON_ERROR = diff

            else
                write (tstr, "(I2)") inb_scal_ql
                write (istr, "(I2)") is
                call TLab_Warning('WARNING: Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

            nullify (qsat,diqsdql,dTdql)

        end select

        !nullify (s_active)
        !nullify (s_EXPL,s_active)

    end subroutine Microphysics_Evaporation_Impl

end module
