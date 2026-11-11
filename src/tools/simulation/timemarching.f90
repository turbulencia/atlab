#include "tlab_error.h"

module TimeMarching
    use TLab_Constants, only: wp, wi, big_wp
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_Time, only: rtime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Log, TLab_Error
#ifdef USE_MPI
    use mpi_f08
    use TLab_WorkFlow, only: ims_err
#endif
    use NavierStokes, only: nse_eqns, DNS_EQNS_COMPRESSIBLE, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use NavierStokes, only: visc, schmidt, prandtl
    use RungeKutta

    implicit none
    private

    public :: TMarch_Initialize
    public :: TMarch_Advance_Step
    public :: TMarch_Compute_Step

    public :: hq        ! needed in statistics and planes; to be removed from public

    ! -----------------------------------------------------------------------
    class(rungekutta_lowstorage_dt), allocatable, protected :: TMarchScheme

    real(wp), public :: dtime                       ! time step
    logical, public :: use_variable_timestep = .true.

    logical :: remove_divergence = .true.           ! Remove residual divergence every time step

    real(wp) :: cfla, cfld                          ! CFL numbers
    real(wp) schmidtfactor, dx2i
    type :: ds_dt
        real(wp), allocatable :: one_ov_ds1(:)
        real(wp), allocatable :: one_ov_ds2(:)
    end type
    type(ds_dt) :: ds(3)

    real(wp), allocatable, target :: hq(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, target :: hs(:, :)       ! Right-hand sides Eulerian fields
    real(wp), pointer :: pxy_hq(:, :, :) => null()
    real(wp), pointer :: pxy_hs(:, :, :) => null()

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_Initialize(inifile)
        use TLab_ScanFiles
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Grid, only: globalGrid, x, y, z, xSubgrid, ySubgrid, zSubgrid
        use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, lstr
        character(len=128) eStr
        character(len=512) sRes
        integer ig
        integer(wi) i, j, k

        real(wp) dummy

        ! ###################################################################
        ! read
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Scheme=<RungeKuttaExplicit3/RungeKuttaExplicit4>')
        call TLab_Write_ASCII(bakfile, '#TimeStep=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxDiffusiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxReactiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#RemoveDivergence=<none/remove>')

        call ScanFile_Char(bakfile, inifile, block, 'Scheme', 'dummy', sRes)
        select case (trim(adjustl(sRes)))
        case ('rungekuttaexplicit3')
            allocate (rk3_dt :: TMarchScheme)
            lstr = '0.6'; 
        case ('rungekuttaexplicit4')
            allocate (rk45_dt :: TMarchScheme)
            lstr = '1.2'; 
        case default
            call TLab_Error(trim(adjustl(eStr))//'Wrong Scheme option.')
            call TLab_Stop(DNS_ERROR_RKORDER)
        end select

        ! Default cfla value set in lstr while reading Scheme
        call ScanFile_Real(bakfile, inifile, block, 'MaxCFL', trim(adjustl(lstr)), cfla)
        write (lstr, *) 0.25_wp*cfla ! Default value for diffusive CFL
        call ScanFile_Real(bakfile, inifile, block, 'MaxDiffusiveCFL', trim(adjustl(lstr)), cfld)

        call ScanFile_Char(bakfile, inifile, block, 'TimeStep', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            use_variable_timestep = .false.
            read (sRes, *) dtime

            call ScanFile_Char(bakfile, inifile, block, 'MaxCFL', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call TLab_Error(trim(adjustl(eStr))//'Cannot impose both time step and max CFL.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end if

        call ScanFile_Char(bakfile, inifile, block, 'RemoveDivergence', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') then; remove_divergence = .false.
        else if (trim(adjustl(sRes)) == 'yes') then; remove_divergence = .true.
        else
            call TLab_Error(trim(adjustl(eStr))//'Wrong RemoveDivergence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        ! Initialize scheme
        call TMarchScheme%initialize()

        ! ###################################################################
        ! Memory management
        call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
        call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

        pxy_hq(1:imax*jmax, 1:kmax, 1:inb_flow) => hq(1:imax*jmax*kmax*inb_flow, 1)
        pxy_hs(1:imax*jmax, 1:kmax, 1:inb_scal) => hs(1:imax*jmax*kmax*inb_scal, 1)

        ! ###################################################################
        ! maximum diffusivities for TMarch_Compute_Step
        schmidtfactor = 1.0_wp
        dummy = 1.0_wp/prandtl
        schmidtfactor = max(schmidtfactor, dummy)
        dummy = 1.0_wp/minval(schmidt(1:inb_scal))
        schmidtfactor = max(schmidtfactor, dummy)

        ! ###################################################################
        ! Calculate inverse of Jacobian (inverse of grid spacing)
        do ig = 1, 3
            allocate (ds(ig)%one_ov_ds1(globalGrid%axes(ig)%size))
            allocate (ds(ig)%one_ov_ds2(globalGrid%axes(ig)%size))

            if (globalGrid%axes(ig)%uniform) then
                if (globalGrid%axes(ig)%size > 1) then
                    ds(ig)%one_ov_ds1(:) = 1.0_wp/(globalGrid%axes(ig)%nodes(2) - globalGrid%axes(ig)%nodes(1))
                else
                    ds(ig)%one_ov_ds1(:) = 1.0_wp   ! 2d case
                end if
            else
                ! we use one_ov_ds2 as auxiliary array here
#define aux(i) ds(ig)%one_ov_ds2(i)
                aux(:) = [(real(i - 1, wp), i=1, globalGrid%axes(ig)%size)]
                select case (ig)
                case (1)
                    call fdm_der1_X%compute(1, aux(:), ds(ig)%one_ov_ds1(:))
                case (2)
                    call fdm_der1_Y%compute(1, aux(:), ds(ig)%one_ov_ds1(:))
                case (3)
                    call fdm_der1_Z%compute(1, aux(:), ds(ig)%one_ov_ds1(:))
                end select
            end if
#undef aux

            ds(ig)%one_ov_ds2(:) = ds(ig)%one_ov_ds1(:)*ds(ig)%one_ov_ds1(:)

        end do

        ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TMarch_Compute_Step
        dx2i = 0.0_wp
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    dummy = 0.0_wp
                    if (x%size > 1) dummy = dummy + ds(1)%one_ov_ds2(i + xSubgrid%offset)
                    if (y%size > 1) dummy = dummy + ds(2)%one_ov_ds2(j + ySubgrid%offset)
                    if (z%size > 1) dummy = dummy + ds(3)%one_ov_ds2(k + zSubgrid%offset)
                    dx2i = max(dx2i, dummy)
                end do
            end do
        end do

#ifdef USE_MPI
        call MPI_ALLREDUCE(dx2i, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        dx2i = dummy
#endif

        call TMarch_Compute_Step()

        return
    end subroutine TMarch_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_Advance_Step()
        use TLab_Arrays, only: q, s, txc
        use TLab_Sources, only: TLab_Sources_Flow, TLab_Sources_Scal, TLab_Sources_Scal_Implicit
        use TLab_Pointers_2D, only: pxy_q, pxy_s
        use DNS_Control, only: DNS_Limit_Bounds
        use Thermo_Anelastic, only: ribackground
        use Buffer
        use BoundaryConditions

        ! -------------------------------------------------------------------
#ifdef USE_PROFILE
        integer(wi) t_srt, t_end, t_dif, idummy, PROC_CYCLES, MAX_CYCLES
        character*256 time_string
#endif

        !########################################################################
        hq(:, :) = 0.0_wp                                   ! Initialize the accumulation of rhs terms
        hs(:, :) = 0.0_wp

        TMarchScheme%time = rtime
        TMarchScheme%substep = 1

        do while (TMarchScheme%substep <= TMarchScheme%num_substep)     ! Loop over the sub-stages
#ifdef USE_PROFILE
            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
#endif

            ! Explicit part
            call TLab_Sources_Flow(q, s, TMarchScheme%time, hq, txc(:, 1))
            call TLab_Sources_Scal(s, hs, TMarchScheme%time, txc(:, 1), txc(:, 2), txc(:, 3))

            call Buffer_Nudge_Flow(q, hq)
            call Buffer_Nudge_Scal(s, hs)

            call NavierStokes_PerVolume(hq, hs, &
                                        TMarchScheme%coef_a(TMarchScheme%substep)*dtime, &
                                        remove_divergence)
            call Bcs_Explicit_Flow(hq)
            call Bcs_Explicit_Scal(hs)

            select case (nse_eqns)
            case (DNS_EQNS_BOUSSINESQ)
                call TMarchScheme%AdvanceSubstep_Boussinesq(pxy_q, pxy_hq, dtime)
                call TMarchScheme%AdvanceSubstep_Boussinesq(pxy_s, pxy_hs, dtime)

            case (DNS_EQNS_ANELASTIC)
                call TMarchScheme%AdvanceSubstep_Anelastic(pxy_q, pxy_hq, dtime, ribackground)
                call TMarchScheme%AdvanceSubstep_Anelastic(pxy_s, pxy_hs, dtime, ribackground)

            case (DNS_EQNS_COMPRESSIBLE)

            end select

            ! Implicit part
            call TLab_Sources_Scal_Implicit(TMarchScheme%coef_t(TMarchScheme%substep)*dtime, s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
            ! Update of boundary condition still missing

            ! Limiters
            call DNS_Limit_Bounds()

            ! Diagnostics
            call TLab_Diagnostic(imax, jmax, kmax, s)

            TMarchScheme%time = TMarchScheme%time + TMarchScheme%coef_t(TMarchScheme%substep)*dtime
            TMarchScheme%substep = TMarchScheme%substep + 1

            ! -------------------------------------------------------------------
            ! Profiling data
#ifdef USE_PROFILE
            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)
            idummy = t_end - t_srt

#ifdef USE_MPI
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
            if (mpiGrid%rank == 0) then
         write (time_string, 999) mpiGrid%num_processors, xMpi%num_processors, zMpi%num_processors, rkm_substep, t_dif/1.0_wp/PROC_CYCLES/mpiGrid%num_processors
999             format(I5.5, ' (xMpi%num_processors X zMpi%num_processors:', I4.4, 'x', I4.4, 1x, ') RK-Substep', I1, ':', E13.5, 's')
                call TLab_Log(time_string)
            end if
#else
            t_dif = idummy
            write (time_string, 999) rkm_substep, t_dif/1.0_wp/PROC_CYCLES/mpiGrid%num_processors
999         format('RK-Substep', I1, ':', E13.5, 's')
            call TLab_Log(time_string)
#endif

#endif
        end do

        return
    end subroutine TMarch_Advance_Step

    !########################################################################
    !#
    !# Determine the variable time step.
    !# For constant time step, this routine
    !# calculates CFL and diffustion numbers for log files
    !#
    !# The diffusion number is fixed in terms of the CFL, which is the input.
    !# This depends on the scheme used. From Lele (1992), page 32, we have that, if
    !# the sixth order tridiagonal scheme is used, then the maximum CFL number
    !# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994
    !# used here we have 3.36/1.989, about 1.69.
    !# This holds for periodic case. A safety margin leads to the common value of 1.2.
    !#
    !# If second order finite different operator is used, then the maximum
    !# diffusion number is 2.9/6.857, about 0.42.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
    !# If the extension by Lamballais et al is used, then the maximum
    !# diffusion number is 2.9/pi^2, about 0.29.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
    !#
    !# If twice the first order finite difference operator is used, then the
    !# maximum diffusion number is 2.9/1.989^2, about 0.73.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
    !#
    !########################################################################
    subroutine TMarch_Compute_Step()
        use TLab_Grid, only: y, xSubgrid, ySubgrid
        use DNS_Control, only: logs_data, logs_dtime
        use TLab_Pointers_3D, only: u, v, w, p_wrk3d

        ! -------------------------------------------------------------------
        integer(wi) i, j, k
        integer(wi) ipmax, j_glo
        real(wp) pmax(3), dtc, dtd
#ifdef USE_MPI
        real(wp) pmax_aux(3)
#endif

        ! ###################################################################
        dtc = big_wp    ! So that the minimum non-zero determines dt at the end
        dtd = big_wp

        ipmax = 0       ! Initialize counter of time constraints

        ! ###################################################################
        ! CFL number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            if (y%size > 1) then
                do k = 1, kmax
                    do j = 1, jmax
                        j_glo = j + ySubgrid%offset
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + xSubgrid%offset) &
                                               + abs(v(i, j, k))*ds(2)%one_ov_ds1(j_glo) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            else    ! do I need this?
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + xSubgrid%offset) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            end if

        end select

        pmax(1) = maxval(p_wrk3d)

        ! ###################################################################
        ! Diffusion number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of \mu*(1/dx^2 + 1/dy^2 + 1/dz^2)
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            pmax(2) = schmidtfactor*visc*dx2i

        end select

        ! ###################################################################
        ! Final operations
        ! ###################################################################
#ifdef USE_MPI
        call MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

        if (use_variable_timestep) then
            if (pmax(1) > 0.0_wp) dtc = cfla/pmax(1) ! Set time step for the given CFL number
            if (pmax(2) > 0.0_wp) dtd = cfld/pmax(2) ! Set time step for the given diffusion number

            dtime = min(dtc, big_wp)
            ! select case (rkm_mode)
            ! case (RKM_EXP3, RKM_EXP4)       ! Explicit diffusion
            dtime = min(dtd, dtime)

            ! end select

        end if

        ! Real CFL and diffusion numbers being used, for the logfile
        logs_dtime = dtime
        logs_data(2) = dtime*pmax(1)
        logs_data(3) = dtime*pmax(2)

        return

    end subroutine TMarch_Compute_Step

end module TimeMarching
