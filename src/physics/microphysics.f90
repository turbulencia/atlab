#include "tlab_error.h"

module Microphysics
    use TLab_Constants, only: wp, wi, pi_wp, efile, MAX_VARS, MAX_PARS, wfile
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use TLab_Memory, only: inb_scal_array
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_Anelastic, only: dsmooth
    use Thermo_AirWater, only: inb_scal_ql, inb_scal_qt, inb_scal_e
    use OPR_Partial, only: OPR_Partial_Z, OPR_P1
    ! use OPR_ODES
    implicit none
    private

    public :: Microphysics_Initialize
    public :: Microphysics_Sedimentation_Z
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
    integer, parameter, public :: TYPE_EVA_QSCONST = 2
    integer, parameter, public :: TYPE_EVA_QSCALC = 3
    integer, parameter, public :: TYPE_EVA_QSCALC_IMPL = 4
    integer, parameter, public :: TYPE_EVA_QSCALC_SEMIIMPL = 5

    real(wp) :: settling                         ! sedimentation effects
    real(wp) :: damkohler                        ! reaction
    ! real(wp) :: stokes                           ! particle inertial effects

contains
!########################################################################
!########################################################################
    subroutine Microphysics_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy

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
        elseif (trim(adjustl(sRes)) == 'qsconst') then; evaporationProps%type = TYPE_EVA_QSCONST
        elseif (trim(adjustl(sRes)) == 'qscalc')  then; evaporationProps%type = TYPE_EVA_QSCALC
        elseif (trim(adjustl(sRes)) == 'qscalc-impl')  then; evaporationProps%type = TYPE_EVA_QSCALC_IMPL
        elseif (trim(adjustl(sRes)) == 'qscalc-semiimpl')  then; evaporationProps%type = TYPE_EVA_QSCALC_SEMIIMPL
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.'//trim(adjustl(sRes)))
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

        end if

        ! -------------------------------------------------------------------
        ! Initialize
        if (imixture == MIXT_TYPE_AIRWATER) then

            select case (evaporationProps%type)

            case(TYPE_EVA_EQUILIBRIUM)
                inb_scal_array = inb_scal_array + 1         ! space for ql as diagnostic variable.
!                use Thermo_Base, only: dsmooth, NEWTONRAPHSON_ERROR
                dsmooth = evaporationProps%parameters(1)    ! Smooth factor for thermodynamic partial derivatives

            case (TYPE_EVA_QSCONST,TYPE_EVA_QSCALC,TYPE_EVA_QSCALC_IMPL,TYPE_EVA_QSCALC_SEMIIMPL)
                evaporationProps%active(inb_scal_ql) = .true.           ! Only ql scalars is affected
                if (damkohler <= 0.0_wp) then
                    call TLab_Write_ASCII(efile, __FILE__//'. Damkohler number must be nonzero if non-equilibrium evaporation is retained.')
                    call TLab_Stop(DNS_ERROR_OPTION)
                end if
                dsmooth = evaporationProps%parameters(3)    ! Smooth factor for thermodynamic partial derivatives

            end select

            if (inb_scal_array <= 3) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'AirWater mixture needs at least 4 scalar arrays.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

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
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
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
                call TLab_Write_ASCII(efile, __FILE__//'. Settling number must be nonzero if sedimentation is retained.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Microphysics_Initialize

!########################################################################
!########################################################################
    subroutine Microphysics_Sedimentation_Z(locProps, nx, ny, nz, is, s, source, tmp1, flux)
        use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_StaticL
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        target s, source

        ! -----------------------------------------------------------------------
        real(wp) dummy, exponent

        real(wp), pointer :: s_active(:) => null()

        !########################################################################
        exponent = locProps%auxiliar(1)
        dummy = 1.0_wp + exponent

        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:, locProps%scalar(is)), source)
            s_active => source
        else
            s_active => s(:, locProps%scalar(is))
        end if

        select case (locProps%type)
        case (TYPE_SED_AIRWATER)
            if (any([inb_scal_qt, inb_scal_ql] == is)) then
                if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*s_active
                end if

            else if (is == inb_scal_e) then
                call Thermo_Anelastic_StaticL(nx, ny, nz, s, tmp1)
                if (exponent > 0.0_wp) then
                    tmp1 = locProps%parameters(is)*tmp1*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*tmp1*s_active
                end if

            end if

            ! select case (is)
            ! case (2, 3)         ! q_t, q_l
            !     if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
            !         tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*(s_active**dummy)
            !     else
            !         tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*s_active
            !     end if

            ! case default        ! energy variables
            !     call Thermo_Anelastic_StaticL(nx, ny, nz, s, tmp1)
            !     if (exponent > 0.0_wp) then
            !         tmp1 = locProps%parameters(is)*tmp1*(s_active**dummy)
            !     else
            !         tmp1 = locProps%parameters(is)*tmp1*s_active
            !     end if

            ! end select

            call OPR_Partial_Z(OPR_P1, nx, ny, nz, tmp1, source)
            if (present(flux)) flux = -tmp1

        case (TYPE_SED_AIRWATERSIMPLIFIED)
            ! if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
            !     tmp1 = locProps%parameters(is)*(s_active**dummy)
            ! else
            !     tmp1 = locProps%parameters(is)*s_active
            ! end if

            ! call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g, tmp1, source)
            ! if (present(flux)) flux = -tmp1

            ! the previous formulation yields oscillations at sharp gradients
            call OPR_Partial_Z(OPR_P1, nx, ny, nz, s_active, tmp1)
            if (exponent > 0.0_wp) tmp1 = tmp1*(s_active**exponent)
            source = locProps%parameters(is)*dummy*tmp1

            if (present(flux)) flux = -locProps%parameters(is)*(s_active**dummy)

        end select

        nullify (s_active)

    end subroutine Microphysics_Sedimentation_Z

!########################################################################
!########################################################################
    subroutine Microphysics_Evaporation(locProps, nx, ny, nz, is, s, source)
        use Thermo_AirWater, only: inb_scal_T, rd_ov_rv
        use Thermo_Base, only: THERMO_PSAT, NPSAT
        use Thermo_Anelastic, only: pbackground, rbackground, Thermo_Anelastic_Weight_OutPlace
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: s(nx*ny, nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny, nz)

        target s, source

        integer(wi) :: k, ipsat
        
        real(wp) :: qsat(nx*ny, nz)

        real(wp) :: dummy, mod_exponent

        real(wp), pointer :: s_active(:,:) => null()

        character(len=512) :: istr,tstr

        source = 0.0_wp

        mod_exponent = locProps%auxiliar(1)

        select case (locProps%type)
        case (TYPE_EVA_QSCONST,TYPE_EVA_QSCALC)
            if (is .eq. inb_scal_ql) then
                
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
#define tmp1 qsat
                    psat(:,:) = THERMO_PSAT(NPSAT)
                    do ipsat = NPSAT - 1, 1, -1
                        psat(:,:) = psat(:,:)*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    do k = 1, nz
                        tmp1(:,k) = rd_ov_rv/(P_LOC(k)/psat(:,k) - 1.0_wp)
                    end do
                    qsat(:,:) = tmp1(:,:)/(1.0_wp + tmp1(:,:))
#undef P_LOC
#undef T_LOC
#undef psat
#undef tmp1
                end select

                ! -------------------------------------------------------------------
                ! Calculate condensation rate as 1/Sq * ( (q_t-q_l)/q_s - 1 ) * q_l^(1+alpha)
                ! -------------------------------------------------------------------
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:,:,inb_scal_ql), source)
                    s_active => source
                else
                    s_active => s(:,:,inb_scal_ql)
                end if
                if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    dummy = 1.0_wp + mod_exponent
                    source = damkohler*locProps%parameters(1) * ( (s(:,:,inb_scal_qt)-s(:,:,inb_scal_ql))/qsat - 1.0_wp ) * sign(1.0_wp,s_active)*(abs(s_active))**dummy 
                else
                    source = damkohler*locProps%parameters(1) * ( (s(:,:,inb_scal_qt)-s(:,:,inb_scal_ql))/qsat - 1.0_wp ) * s_active
                end if

            else
                write(tstr,"(I2)") inb_scal_ql
                write(istr,"(I2)") is
                call TLab_Write_ASCII(wfile, 'WARNING: Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

        end select

        nullify (s_active)

    end subroutine Microphysics_Evaporation

!########################################################################
!########################################################################
    subroutine Microphysics_Evaporation_Impl(locProps, nx, ny, nz, is, s, dte)
        use Thermo_AirWater, only: inb_scal_T, rd_ov_rv, Cd, Cdv, Cvl, Lv0
        use Thermo_Base, only: THERMO_PSAT, NPSAT
        use Thermo_Anelastic, only: pbackground, rbackground, ribackground, epbackground
        use Thermo_Anelastic, only: Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_Weight_InPlace
        use Thermo_Anelastic, only: NEWTONRAPHSON_ERROR
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_err
        real(wp) diff_glob
#endif
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in)    :: dte
        real(wp), intent(inout) :: s(nx*ny, nz, inb_scal_array)
        

        integer(wi) :: k, ipsat, it
        
        real(wp), dimension(nx*ny,nz) :: s_RHS, tmp1, tmp2, qsat, diqsdql, dTdql

        real(wp) :: mod_exponent, diff, eps

        real(wp), dimension(nx*ny,nz) :: s_active
        !real(wp), pointer :: s_active(:,:) => null()

        !target s, tmp2

        character(len=512) :: istr,tstr


        
        
        mod_exponent = locProps%auxiliar(1)

        select case (locProps%type)
        case (TYPE_EVA_QSCALC_IMPL,TYPE_EVA_QSCALC_SEMIIMPL)
            if (is .eq. inb_scal_ql) then

                if ( locProps%type == TYPE_EVA_QSCALC_SEMIIMPL ) then
                    ! -------------------------------------------------------------------
                    ! Update Temperature
                    ! -------------------------------------------------------------------
                    call TLab_Diagnostic(nx,ny,nz,s)

                    ! -------------------------------------------------------------------
                    ! Update saturation specific humidity q_sat(T,p)
                    ! -------------------------------------------------------------------
#define P_LOC pbackground
#define T_LOC s(:,:,inb_scal_T)
#define psat qsat
#define qstmp qsat
                    psat(:,:) = THERMO_PSAT(NPSAT)
                    do ipsat = NPSAT - 1, 1, -1
                        psat(:,:) = psat(:,:)*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    do k = 1, nz
                        qstmp(:,k) = rd_ov_rv/(P_LOC(k)/psat(:,k) - 1.0_wp)
                    end do
                    qsat(:,:) = qstmp(:,:)/(1.0_wp + qstmp(:,:))
#undef P_LOC
#undef T_LOC
#undef psat
#undef qstmp
                end if

                diqsdql(:,:)=0.0_wp
                s_RHS(:,:) = s(:,:,inb_scal_ql)
                it=0
                eps=1.e-12_wp
                diff=1.0_wp
                do while (diff>eps .and. it<=10)
                !do while (it<=10)
                    if ( locProps%type == TYPE_EVA_QSCALC_IMPL ) then
                        ! -------------------------------------------------------------------
                        ! Update Temperature
                        ! -------------------------------------------------------------------
                        call TLab_Diagnostic(nx,ny,nz,s)

                        ! -------------------------------------------------------------------
                        ! Update saturation specific humidity q_sat(T,p)
                        ! -------------------------------------------------------------------
#define P_LOC pbackground
#define T_LOC s(:,:,inb_scal_T)
#define psat qsat
#define qstmp qsat

#define E_LOC epbackground
#define dpsdT diqsdql
                        psat(:,:) = THERMO_PSAT(NPSAT)
                        dpsdT(:,:) = 0.0_wp
                        do ipsat = NPSAT - 1, 1, -1
                            psat(:,:) = psat(:,:)*T_LOC + THERMO_PSAT(ipsat)
                            dpsdT(:,:) = dpsdT(:,:)*T_LOC + THERMO_PSAT(ipsat+1)*ipsat
                        end do
                        do k = 1, nz
                            dTdql(:,k) = ( Lv0*(Cd+Cdv*s(:,k,inb_scal_qt)) - Cvl*(s(:,k,inb_scal_e)-E_LOC(k)) )/(( Cd + Cdv*s(:,k,inb_scal_qt) + Cvl*s(:,k,inb_scal_ql) )**2)
                            diqsdql(:,k)= -P_LOC(k)/(psat(:,k)**2)/rd_ov_rv * dpsdT(:,k) * dTdql(:,k)
                            qstmp(:,k) = rd_ov_rv/(P_LOC(k)/psat(:,k) - 1.0_wp)
                        end do
                        qsat(:,:) = qstmp(:,:)/(1.0_wp + qstmp(:,:))
#undef P_LOC
#undef T_LOC
#undef psat
#undef qstmp

#undef E_LOC
#undef dpsdT
                    end if

                    ! -------------------------------------------------------------------
                    ! Calculate condensation rate as 1/Sq * ( (q_t-q_l)/q_s - 1 ) * q_l^(1+alpha)
                    ! -------------------------------------------------------------------
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        !call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:,:,inb_scal_ql), tmp2)
                        !s_active => tmp2
                        call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:,:,inb_scal_ql), s_active)
                    else
                        !s_active => s(:,:,inb_scal_ql)
                        s_active = s(:,:,inb_scal_ql)
                    end if

                    s_active = max(s_active, 0.0_wp)

                    ! --- Derivative df(q_l)/dq_l for NRmethod --- 
                    ! ( rescaled by q_l^alpha to avoid division by zero )
                    tmp1 = ( (s(:,:,inb_scal_qt)-s(:,:,inb_scal_ql))/qsat - 1.0_wp ) * (1.0_wp + mod_exponent)
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, tmp1)
                    end if
                    tmp1 = tmp1  + ( -1.0_wp/qsat + (s(:,:,inb_scal_qt)-s(:,:,inb_scal_ql))*diqsdql ) * s_active
                    tmp1 = tmp1  * damkohler*locProps%parameters(1)
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, ribackground, tmp1)
                    end if
                    if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        !tmp1 = sign(1.0_wp,s_active)*(abs(s_active))**(-mod_exponent) - dte*tmp1
                        tmp1 = s_active**(-mod_exponent) - dte*tmp1
                    else
                        tmp1 = 1.0_wp - dte*tmp1
                    end if

                    ! --- Source term delta f(q_l) for NRmethod ---
                    ! ( rescaled by q_l^alpha to avoid division by zero )
                    tmp2 = damkohler*locProps%parameters(1) * ( (s(:,:,inb_scal_qt)-s(:,:,inb_scal_ql))/qsat - 1.0_wp ) * s_active
                    if (nse_eqns == DNS_EQNS_ANELASTIC) then
                        call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, ribackground, tmp2)
                    end if
                    if (mod_exponent /= 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                        !tmp2 = sign(1.0_wp,s_active)*(abs(s_active))**(-mod_exponent) * (s(:,:,inb_scal_ql) - s_RHS) - dte*tmp2
                        tmp2 = s_active**(-mod_exponent) * (s(:,:,inb_scal_ql) - s_RHS) - dte*tmp2
                    else
                        tmp2 = s(:,:,inb_scal_ql) - s_RHS - dte*tmp2
                    end if

                    ! -------------------------------------------------------------------
                    ! Newton-Raphson iteration
                    ! -------------------------------------------------------------------
                    s(:,:,inb_scal_ql) = s(:,:,inb_scal_ql) - tmp2/tmp1

                    diff = MAXVAL( ABS( tmp2/tmp1 / s(:,:,inb_scal_ql) ) )
#ifdef USE_MPI
                    call MPI_ALLREDUCE(diff, diff_glob, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
                    diff = diff_glob
#endif

                    
                    it = it + 1
                end do
                NEWTONRAPHSON_ERROR = diff

            else
                write(tstr,"(I2)") inb_scal_ql
                write(istr,"(I2)") is
                call TLab_Write_ASCII(wfile, 'WARNING: Unexpectedly Scalar'//trim(istr)//' (other than q_l=Scalar'//trim(tstr)//') set to active for NON-EQUILIBRIUM EVAPORATION')

            end if

        end select

        !nullify (s_active)

    end subroutine Microphysics_Evaporation_Impl

end module
