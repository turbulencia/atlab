#include "tlab_error.h"

!########################################################################
!#
!# Runge-Kutta explicit 3th order from Williamson 1980
!# Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
!# Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
!#
!########################################################################
module TimeMarching

    use TLab_Constants, only: efile, wfile, wp, wi, big_wp
    use TLab_WorkFlow, only: flow_on, scal_on
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    ! use PARTICLE_VARS
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
#endif
    use TLab_Grid, only: x, y, z
    use NavierStokes, only: nse_eqns, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use NavierStokes, only: nse_advection, EQNS_CONVECTIVE, EQNS_DIVERGENCE, EQNS_SKEWSYMMETRIC
    use NavierStokes, only: visc, schmidt, prandtl
    use DNS_Arrays
    use BoundaryConditions
    use Buffer

    implicit none
    private

    public :: TMarch_Initialize
    public :: TMarch_RungeKutta
    public :: TMarch_Courant

    real(wp), public :: dtime                       ! time step
    real(wp), public :: dte                         ! time step of each substep
    logical, public :: remove_divergence            ! Remove residual divergence every time step
    logical, public :: use_variable_timestep = .true.

    ! -------------------------------------------------------------------
    ! type :: tmarch_dt
    !     sequence
    !     integer type
    !     integer nb_stages                           ! number of stages
    !     real(wp) kdt(5), kco(4), ktime(5)           ! explicit scheme coefficients
    !     real(wp) kex(3), kim(3)                     ! implicit scheme coefficients
    !     procedure(tmarch_interface), pointer, nopass :: tmarch_scheme
    ! end type tmarch_dt

    integer(wi) :: rkm_mode                     ! Type of Runge-Kutta scheme
    integer, parameter :: RKM_EXP3 = 3
    integer, parameter :: RKM_EXP4 = 4
    integer, parameter :: RKM_IMP3_DIFFUSION = 5
    integer, parameter :: RKM_IMP3_SOURCE = 6
    integer, parameter :: RKM_IMP3_DIFFSOURCE = 7

    integer(wi) :: rkm_endstep                  ! number of substeps
    integer(wi) :: rkm_substep                  ! substep counter

    real(wp) :: cfla, cfld, cflr                ! CFL numbers
    real(wp) etime                              ! time at each substep

    real(wp) kdt(5), kco(4), ktime(5)           ! explicit scheme coefficients
    real(wp) kex(3), kim(3)                     ! implicit scheme coefficients

    real(wp) schmidtfactor, dx2i
    integer(wi) i, j, k, kdsp, jdsp, idsp, is
    real(wp) dummy

    type :: ds_dt
        real(wp), allocatable :: one_ov_ds1(:)
        real(wp), allocatable :: one_ov_ds2(:)
    end type
    type(ds_dt) :: ds(3)

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_Initialize(inifile)
        use TLab_Memory, only: TLab_Allocate_Real
        use FDM, only: g

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, lstr
        character(len=128) eStr
        character(len=512) sRes
        integer ig

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Scheme=<RungeKuttaExplicit3/RungeKuttaExplicit4/RungeKuttaDiffusion3>')
        call TLab_Write_ASCII(bakfile, '#TimeStep=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxDiffusiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxReactiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#RemoveDivergence=<none/remove>')

        call ScanFile_Char(bakfile, inifile, block, 'Scheme', 'dummy', sRes)
        if (trim(adjustl(sRes)) == 'rungekuttaexplicit3') then; rkm_mode = RKM_EXP3; lstr = '0.6'; 
        elseif (trim(adjustl(sRes)) == 'rungekuttaexplicit4') then; rkm_mode = RKM_EXP4; lstr = '1.2'; 
        elseif (trim(adjustl(sRes)) == 'rungekuttadiffusion3') then; rkm_mode = RKM_IMP3_DIFFUSION; lstr = '0.6'; 
            !  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttasource3'    ) THEN; rkm_mode = RKM_IMP3_SOURCE;
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Scheme option.')
            call TLab_Stop(DNS_ERROR_RKORDER)
        end if

        ! Default cfla value set in lstr while reading Scheme
        call ScanFile_Real(bakfile, inifile, block, 'MaxCFL', trim(adjustl(lstr)), cfla)
        write (lstr, *) 0.25_wp*cfla ! Default value for diffusive CFL
        call ScanFile_Real(bakfile, inifile, block, 'MaxDiffusiveCFL', trim(adjustl(lstr)), cfld)
        write (lstr, *) 0.5_wp*cfla ! Default value for reactive CFL
        call ScanFile_Real(bakfile, inifile, block, 'MaxReactiveCFL', trim(adjustl(lstr)), cflr)

        call ScanFile_Char(bakfile, inifile, block, 'TimeStep', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            use_variable_timestep = .false.
            read (sRes, *) dtime

            call ScanFile_Char(bakfile, inifile, block, 'MaxCFL', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Cannot impose both time step and max CFL.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end if

        call ScanFile_Char(bakfile, inifile, block, 'RemoveDivergence', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') then; remove_divergence = .false.
        else if (trim(adjustl(sRes)) == 'yes') then; remove_divergence = .true.
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong RemoveDivergence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! -------------------------------------------------------------------
        ! Consistency check
        if (rkm_mode == RKM_IMP3_DIFFUSION) then
            do is = 1, inb_scal
                if (BcsScalKmin%type(is) == DNS_BCS_Neumann .or. &
                    BcsScalKmax%type(is) == DNS_BCS_Neumann) then
                    write (sRes, *) is; sRes = trim(adjustl(eStr))//'Scalar'//trim(adjustl(sRes))//'. Finite flux BC not implemented for SEMI-IMPLICITE DIFFUSION'
                    call TLab_Write_ASCII(wfile, trim(adjustl(sRes)))
                    write (sRes, *) is; sRes = trim(adjustl(eStr))//'Scalar'//trim(adjustl(sRes))//'. Setting fluxes at boundary to zero'
                    call TLab_Write_ASCII(wfile, trim(adjustl(sRes)))
                end if
            end do

        end if

        ! ###################################################################
        ! RK coefficients
        select case (rkm_mode)
        case (RKM_EXP3)             ! Runge-Kutta explicit 3th order from Williamson 1980
            rkm_endstep = 3

            kdt(1:3) = [1.0_wp/3.0_wp, 15.0_wp/16.0_wp, 8.0_wp/15.0_wp]
            ktime(1:3) = [0.0_wp, 1.0_wp/3.0_wp, 3.0_wp/4.0_wp]
            kco(1:2) = [-5.0_wp/9.0_wp, -153.0_wp/128.0_wp]

        case (RKM_EXP4)             ! Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
            rkm_endstep = 5

            kdt(1) = 1432997174477.0_wp/9575080441755.0_wp
            kdt(2) = 5161836677717.0_wp/13612068292357.0_wp
            kdt(3) = 1720146321549.0_wp/2090206949498.0_wp
            kdt(4) = 3134564353537.0_wp/4481467310338.0_wp
            kdt(5) = 2277821191437.0_wp/14882151754819.0_wp

            ktime(1) = 0.0_wp
            ktime(2) = kdt(1)
            ktime(3) = 2526269341429.0_wp/6820363962896.0_wp
            ktime(4) = 2006345519317.0_wp/3224310063776.0_wp
            ktime(5) = 2802321613138.0_wp/2924317926251.0_wp

            kco(1) = -567301805773.0_wp/1357537059087.0_wp
            kco(2) = -2404267990393.0_wp/2016746695238.0_wp
            kco(3) = -3550918686646.0_wp/2091501179385.0_wp
            kco(4) = -1275806237668.0_wp/842570457699.0_wp

        case (RKM_IMP3_DIFFUSION)   ! Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
            rkm_endstep = 3

            kdt(1:3) = [8.0_wp/15.0_wp, 5.0_wp/12.0_wp, 3.0_wp/4.0_wp]

            kim(1:3) = [111.0_wp/256.0_wp, 1.0_wp/2.0_wp, 2.0_wp/9.0_wp]
            kex(1:3) = [145.0_wp/256.0_wp, -9.0_wp/50.0_wp, 2.0_wp/9.0_wp]
            kco(1:3) = [0.0_wp, -17.0_wp/25.0_wp, -5.0_wp/9.0_wp]
            ! TO DO - calculate ktime from coefficients  ktime
            ktime(1:3) = [0.0_wp, 0.0_wp, 0.0_wp]

            ! Coefficients from Spalart, Moser, Rogers (1991)
            ! kim = beta/gamma
            ! kex = alpha/gamma
            ! kco = zeta/gamma
            !
            ! alpha = [ 29./96.,   -3./40,    1./6. ]
            ! beta  = [ 37./160.,   5./24.,   1./6. ]
            ! gamma = [  8./15.,    5./12.,   3./4. ]
            ! zeta  = [  0.,      -17./60.,  -5./12.]

        end select

        ! ###################################################################
        ! Memory management
        call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
        call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

        p_hq(1:imax, 1:jmax, 1:kmax, 1:inb_flow) => hq(1:imax*jmax*kmax*inb_flow, 1)
        p_hs(1:imax, 1:jmax, 1:kmax, 1:inb_scal) => hs(1:imax*jmax*kmax*inb_scal, 1)

        ! ###################################################################
        ! maximum diffusivities for TMarch_Courant
        schmidtfactor = 1.0_wp
        dummy = 1.0_wp/prandtl
        schmidtfactor = max(schmidtfactor, dummy)
        dummy = 1.0_wp/minval(schmidt(1:inb_scal))
        schmidtfactor = max(schmidtfactor, dummy)

        ! ###################################################################
        do ig = 1, 3
            allocate (ds(ig)%one_ov_ds1(g(ig)%size))
            ds(ig)%one_ov_ds1(:) = 1.0_wp/g(ig)%jac(:, 1)

            allocate (ds(ig)%one_ov_ds2(g(ig)%size))
            ds(ig)%one_ov_ds2(:) = ds(ig)%one_ov_ds1(:)*ds(ig)%one_ov_ds1(:)

        end do

        ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TMarch_Courant
#ifdef USE_MPI
        idsp = ims_offset_i
        jdsp = ims_offset_j
        kdsp = ims_offset_k
#else
        idsp = 0
        jdsp = 0
        kdsp = 0
#endif

        dx2i = 0.0_wp
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    dummy = 0.0_wp
                    if (x%size > 1) dummy = dummy + ds(1)%one_ov_ds2(i + idsp)
                    if (y%size > 1) dummy = dummy + ds(2)%one_ov_ds2(j + jdsp)
                    if (z%size > 1) dummy = dummy + ds(3)%one_ov_ds2(k + kdsp)
                    dx2i = max(dx2i, dummy)
                end do
            end do
        end do

#ifdef USE_MPI
        call MPI_ALLREDUCE(dx2i, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        dx2i = dummy
#endif

        call TMarch_Courant()

        return
    end subroutine TMarch_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_RungeKutta()
        use TLab_Arrays
        use DNS_LOCAL
        use DNS_Control
        use DNS_Arrays

        ! -------------------------------------------------------------------
        real(wp) alpha

#ifdef USE_PROFILE
        integer(wi) t_srt, t_end, t_dif, idummy, PROC_CYCLES, MAX_CYCLES
        character*256 time_string
#endif

        !########################################################################
        ! -------------------------------------------------------------------
        ! Initialize arrays to zero for the explcit low-storage algorithm
        ! -------------------------------------------------------------------
        if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
            if (flow_on) hq = 0.0_wp
            if (scal_on) hs = 0.0_wp
            ! if (part%type /= PART_TYPE_NONE) l_hq = 0.0_wp
        end if
        !########################################################################
        ! Loop over the sub-stages
        !########################################################################
        do rkm_substep = 1, rkm_endstep

            ! -------------------------------------------------------------------
            ! Update transported (or prognostic) variables q and s
            ! -------------------------------------------------------------------
            dte = dtime*kdt(rkm_substep)
            etime = rtime + dtime*ktime(rkm_substep)

#ifdef USE_PROFILE
            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
#endif
            ! I could define procedure pointers to handle this...
            select case (nse_eqns)
            case (DNS_EQNS_BOUSSINESQ)
                if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                    call TMarch_Substep_Boussinesq_Explicit()
                else if (rkm_mode == RKM_IMP3_DIFFUSION) then
                    call TMarch_Substep_Boussinesq_Implicit()
                end if

            case (DNS_EQNS_ANELASTIC)
                if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                    call TMarch_Substep_Anelastic_Explicit()
                    ! else if (rkm_mode == RKM_IMP3_DIFFUSION) then
                    !     call TMarch_Substep_Boussinesq_Implicit()
                end if

            case (DNS_EQNS_COMPRESSIBLE)

            end select

            call DNS_BOUNDS_LIMIT()
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call TLab_Diagnostic(imax, jmax, kmax, s)
            end if
            if (int(logs_data(1)) /= 0) return ! Error detected

            ! -------------------------------------------------------------------
            ! Update RHS hq and hs in the explicit low-storage algorithm
            ! -------------------------------------------------------------------
            if ((rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) .and. &
                rkm_substep < rkm_endstep) then

                alpha = kco(rkm_substep)

                if (flow_on) then
                    do is = 1, inb_flow
#ifdef USE_BLAS
                        call DSCAL(isize_field, alpha, hq(:, is), 1)
#else
                        hq(:, is) = alpha*hq(:, is)
#endif
                    end do
                end if

                if (scal_on) then
                    do is = 1, inb_scal
#ifdef USE_BLAS
                        call DSCAL(isize_field, alpha, hs(:, is), 1)
#else
                        hs(:, is) = alpha*hs(:, is)
#endif
                    end do
                end if

            end if

            ! -------------------------------------------------------------------
            ! Profiling data
            ! -------------------------------------------------------------------
#ifdef USE_PROFILE
            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)
            idummy = t_end - t_srt

#ifdef USE_MPI
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
            if (ims_pro == 0) then
                write (time_string, 999) ims_npro, ims_npro_i, ims_npro_k, rkm_substep, t_dif/1.0_wp/PROC_CYCLES/ims_npro
999             format(I5.5, ' (ims_npro_i X ims_npro_k:', I4.4, 'x', I4.4, 1x, ') RK-Substep', I1, ':', E13.5, 's')
                call TLab_Write_ASCII(lfile, time_string)
            end if
#else
            t_dif = idummy
            write (time_string, 999) rkm_substep, t_dif/1.0_wp/PROC_CYCLES/ims_npro
999         format('RK-Substep', I1, ':', E13.5, 's')
            call TLab_Write_ASCII(lfile, time_string)
#endif

#endif
        end do

        return
    end subroutine TMarch_RungeKutta

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
    subroutine TMarch_Courant()
        use DNS_Control, only: logs_data, logs_dtime
        use TLab_Pointers_3D, only: u, v, w, p_wrk3d

        ! -------------------------------------------------------------------
        integer(wi) ipmax, j_glo
        real(wp) pmax(3), dtc, dtd
#ifdef USE_MPI
        real(wp) pmax_aux(3)
#endif

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i
        jdsp = ims_offset_j
#else
        idsp = 0
        jdsp = 0
#endif

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
                        j_glo = j + jdsp
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + idsp) &
                                               + abs(v(i, j, k))*ds(2)%one_ov_ds1(j_glo) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            else    ! do I need this?
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + idsp) &
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
            select case (rkm_mode)
            case (RKM_EXP3, RKM_EXP4)       ! Explicit diffusion
                dtime = min(dtd, dtime)

            end select

        end if

        ! Real CFL and diffusion numbers being used, for the logfile
        logs_dtime = dtime
        logs_data(2) = dtime*pmax(1)
        logs_data(3) = dtime*pmax(2)

        return

    end subroutine TMarch_Courant

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Boussinesq_Explicit()
        use TLab_Arrays, only: q, s, txc
        use DNS_Arrays, only: hq, hs
        use TLab_Sources

        ! #######################################################################
        ! Accumulate RHS terms
        call TLab_Sources_Flow(q, s, hq, txc(:, 1))
        call TLab_Sources_Scal(s, hs, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))

        if (bufferType == BUFFER_TYPE_NUDGE) call Buffer_Nudge()

        call NSE_Boussinesq()

        ! #######################################################################
        ! Perform the time stepping
        do is = 1, inb_flow
            q(:, is) = q(:, is) + dte*hq(:, is)
        end do

        do is = 1, inb_scal
            s(:, is) = s(:, is) + dte*hs(:, is)
        end do

        return
    end subroutine TMarch_Substep_Boussinesq_Explicit

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Anelastic_Explicit()
        use TLab_Arrays, only: q, s, txc
        use DNS_Arrays, only: hq, hs
        use TLab_Sources
        use Microphysics, only: Microphysics_Evaporation_Impl, evaporationProps
        use Thermo_AirWater, only: inb_scal_ql

        ! #######################################################################
        ! Accumulate RHS terms
        call TLab_Sources_Flow(q, s, hq, txc(:, 1))
        call TLab_Sources_Scal(s, hs, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))

        if (bufferType == BUFFER_TYPE_NUDGE) call Buffer_Nudge()

        call NSE_Anelastic()

        ! #######################################################################
        ! Perform the time stepping
        do is = 1, inb_flow
            q(:, is) = q(:, is) + dte*hq(:, is)
        end do

        do is = 1, inb_scal
            s(:, is) = s(:, is) + dte*hs(:, is)
        end do

        ! #######################################################################
        ! Iterate implicit non-evaporation
        call Microphysics_Evaporation_Impl(evaporationProps,imax,jmax,kmax,inb_scal_ql,s,dte)

        ! #######################################################################
        !call TLab_Diagnostic(imax, jmax, kmax, s) !RH: tbd after re-enforcing scalar limits -> moved to TMarch_RungeKutta

        return
    end subroutine TMarch_Substep_Anelastic_Explicit

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Boussinesq_Implicit()

        ! call RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(kex(rkm_substep), kim(rkm_substep), kco(rkm_substep))

        ! ! ! pressure-correction algorithm; to be checked
        ! ! CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(&
        ! !      kex,kim,kco,  &
        ! !      q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
        ! !      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), txc(1,8))

        return
    end subroutine TMarch_Substep_Boussinesq_Implicit

    !     !########################################################################
    !     !########################################################################
    !     subroutine TMarch_SUBSTEP_PARTICLE()
    !         use DNS_Arrays, only: l_hq
    !         use PARTICLE_VARS
    !         use PARTICLE_ARRAYS

    !         ! -------------------------------------------------------------------
    !         integer(wi) is

    ! #ifdef USE_MPI
    !         integer(wi) nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
    ! #endif
    !         real(wp) x_right, z_right
    !         real(wp) y_right

    !         !#####################################################################
    !         call RHS_PART_1()

    !         !#######################################################################
    !         ! Update particle properties
    !         !#######################################################################
    !         do is = 1, inb_part
    !             l_q(1:l_g%np, is) = l_q(1:l_g%np, is) + dte*l_hq(1:l_g%np, is)
    !         end do

    !         !#####################################################################
    !         ! Boundary control to see if particles leave processor
    !         !#####################################################################
    ! #ifdef USE_MPI

    !         ! -------------------------------------------------------------------
    !         ! Particle sorting for Send/Recv X-Direction
    !         ! -------------------------------------------------------------------
    !         if (ims_npro_i > 1) then
    !             call PARTICLE_MPI_SORT(1, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    !             if (ims_pro_i == 0) then !Take care of periodic boundary conditions west
    !                 if (nzone_west /= 0) then
    !                     l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) = &
    !                         l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) + x%scale
    !                 end if
    !             end if

    !             if (ims_pro_i == (ims_npro_i - 1)) then !Take care of periodic boundary conditions east
    !                 if (nzone_east /= 0) then
    !                     l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) = &
    !                         l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) - x%scale
    !                 end if
    !             end if

    !             call PARTICLE_MPI_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, l_q, l_hq, l_g%tags, l_g%np)

    !         else
    !             x_right = x%nodes(1) + x%scale

    !             do i = 1, l_g%np
    !                 if (l_q(i, 1) > x_right) then
    !                     l_q(i, 1) = l_q(i, 1) - x%scale

    !                 elseif (l_q(i, 1) < x%nodes(1)) then
    !                     l_q(i, 1) = l_q(i, 1) + x%scale

    !                 end if

    !             end do

    !         end if

    !         ! -------------------------------------------------------------------
    !         ! Particle sorting for Send/Recv Z-Direction
    !         ! -------------------------------------------------------------------
    !         if (ims_npro_k > 1) then
    !             call PARTICLE_MPI_SORT(3, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    !             if (ims_pro_k == 0) then !Take care of periodic boundary conditions south
    !                 if (nzone_south /= 0) then
    !                     l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) = &
    !                         l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) + z%scale
    !                 end if
    !             end if

    !             if (ims_pro_k == (ims_npro_k - 1)) then !Take care of periodic boundary conditions north
    !                 if (nzone_north /= 0) then
    !                     l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) = &
    !                         l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) - z%scale
    !                 end if
    !             end if

    !             call PARTICLE_MPI_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, l_q, l_hq, l_g%tags, l_g%np)

    !         else
    !             call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    !             z_right = z%nodes(1) + z%scale

    !             do i = 1, l_g%np
    !                 if (l_q(i, 3) > z_right) then
    !                     l_q(i, 3) = l_q(i, 3) - z%scale

    !                 elseif (l_q(i, 3) < z%nodes(1)) then
    !                     l_q(i, 3) = l_q(i, 3) + z%scale

    !                 end if

    !             end do

    !         end if

    ! #else
    !         !#######################################################################
    !         ! Serial; would it be faster to use MOD?
    !         x_right = x%nodes(1) + x%scale
    !         z_right = z%nodes(1) + z%scale

    !         do i = 1, l_g%np
    !             if (l_q(i, 1) > x_right) then
    !                 l_q(i, 1) = l_q(i, 1) - x%scale

    !             elseif (l_q(i, 1) < x%nodes(1)) then
    !                 l_q(i, 1) = l_q(i, 1) + x%scale

    !             end if

    !             if (l_q(i, 3) > z_right) then
    !                 l_q(i, 3) = l_q(i, 3) - z%scale

    !             elseif (l_q(i, 3) < z%nodes(1)) then
    !                 l_q(i, 3) = l_q(i, 3) + z%scale

    !             end if

    !         end do

    ! #endif

    !         y_right = y%nodes(1) + y%scale
    !         select case (part_bcs)
    !         case (PART_BCS_SPECULAR)
    !             do i = 1, l_g%np
    !                 if (l_q(i, 2) > y_right) then
    !                     l_q(i, 2) = 2*y_right - l_q(i, 2)
    !                     l_q(i, 5) = -l_q(i, 5)
    !                 elseif (l_q(i, 2) < y%nodes(1)) then
    !                     l_q(i, 2) = 2*y%nodes(1) - l_q(i, 2)
    !                     l_q(i, 5) = -l_q(i, 5)
    !                 end if
    !             end do

    !         case (PART_BCS_STICK)
    !             do i = 1, l_g%np
    !                 if (l_q(i, 2) > y_right) then
    !                     l_q(i, 2) = y_right
    !                 elseif (l_q(i, 2) < y%nodes(1)) then
    !                     l_q(i, 2) = y%nodes(1)
    !                 end if
    !             end do

    !         end select

    !         !#######################################################################
    !         ! Recalculating closest node below in Y direction
    !         !#######################################################################
    !         call LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, y%size, y%nodes)

    !         return
    !     end subroutine TMarch_SUBSTEP_PARTICLE

end module TimeMarching
