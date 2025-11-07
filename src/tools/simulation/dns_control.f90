#include "tlab_error.h"

module DNS_Control
    use TLab_Constants, only: MAX_VARS, wp, wi, sp, big_wp
    use TLab_Constants, only: efile, lfile
    use TLab_Memory, only: inb_scal
    use TLab_WorkFlow, only: scal_on, flow_on, TLab_Stop, TLab_Write_ASCII
    use TLab_Time
    use NavierStokes
    ! use TLab_Constants, only: MAX_PATH_LENGTH
    ! use TIME
    use DNS_LOCAL, only: wall_time, start_clock
    implicit none

    character(len=*), parameter :: ofile = 'dns.out'    ! data logger filename
    ! character(len=*), parameter :: ofile_base = 'dns.out'    ! data logger filename
    character(len=*), parameter :: vfile = 'dns.obs'    ! insitu obs. logger filename
    ! character(len=*), parameter :: vfile_base = 'dns.obs'    ! insitu obs. logger filename
    ! character(len=MAX_PATH_LENGTH) :: ofile, vfile
    ! character(len=MAX_PATH_LENGTH) :: logger_path
    real(wp) :: logs_dtime, logs_data(20)       ! information (time, time step, cfls, dilatation...)

    real(wp) :: obs_data(20)        ! information (custom variables / insitu measurements ...)
    integer :: dns_obs_log

    type bounds_dt
        sequence
        logical active
        real(wp) min
        real(wp) max
    end type bounds_dt

    type(bounds_dt) bound_p             ! limit pressure in compressible flows
    type(bounds_dt) bound_r             ! limit density in compressible flows
    type(bounds_dt) bound_s(MAX_VARS)    ! limit scalars
    type(bounds_dt) bound_d             ! control dilatation in incompressible/anelastic flows

    ! Obs log-file type
    integer, parameter :: OBS_TYPE_NONE = 0
    integer, parameter :: OBS_TYPE_EKMAN = 1

contains
    !########################################################################
    !########################################################################
    subroutine DNS_Control_Initialize(inifile)
        use TLab_Memory, only: inb_flow, inb_scal
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer idummy
        real(wp) dummy(inb_flow + inb_scal + 1)
        integer inb_scal_local1

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#ObsLog=<None/Ekman>')
        call ScanFile_Char(bakfile, inifile, block, 'ObsLog', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; dns_obs_log = OBS_TYPE_NONE
        else if (trim(adjustl(sRes)) == 'ekman') then; dns_obs_log = OBS_TYPE_EKMAN
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'. Incorrect type ObsLog.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        ! ###################################################################
        block = 'Control'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[Control]')
        call TLab_Write_ASCII(bakfile, '#FlowLimit=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#MinPressure=<pressure>')
        call TLab_Write_ASCII(bakfile, '#MaxPressure=<pressure>')
        call TLab_Write_ASCII(bakfile, '#MinDensity=<density>')
        call TLab_Write_ASCII(bakfile, '#MaxDensity=<density>')
        call TLab_Write_ASCII(bakfile, '#ScalLimit=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#MinScalar=<scalar>')
        call TLab_Write_ASCII(bakfile, '#MaxScalar=<scalar>')

        bound_r%active = .false.
        bound_p%active = .false.
        call ScanFile_Char(bakfile, inifile, 'Control', 'FlowLimit', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') then
            bound_r%active = .true.
            bound_p%active = .true.
        end if
        ! Final check in last section of the routine
        call ScanFile_Real(bakfile, inifile, 'Control', 'MinPressure', '-1.0', bound_p%min)
        call ScanFile_Real(bakfile, inifile, 'Control', 'MaxPressure', '-1.0', bound_p%max)
        call ScanFile_Real(bakfile, inifile, 'Control', 'MinDensity', '-1.0', bound_r%min)
        call ScanFile_Real(bakfile, inifile, 'Control', 'MaxDensity', '-1.0', bound_r%max)

        bound_d%active = .false.
        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            bound_d%active = .true.
            bound_d%max = big_wp ! default
            call ScanFile_Char(bakfile, inifile, 'Control', 'MaxDilatation', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = 1
                call LIST_REAL(sRes, idummy, dummy)
                bound_d%max = dummy(1)
            end if
        end if

        bound_s(:)%active = .false.
        call ScanFile_Char(bakfile, inifile, 'Control', 'ScalLimit', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') then
            bound_s(:)%active = .true.
        end if

        bound_s(:)%min = 0.0_wp; inb_scal_local1 = MAX_VARS
        if (any(bound_s(:)%active)) then
            call ScanFile_Char(bakfile, inifile, 'Control', 'MinScalar', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call LIST_REAL(sRes, inb_scal_local1, dummy)
                if (inb_scal_local1 /= inb_scal) then ! Consistency check
                    call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'MinScalar size does not match inb_scal.')
                    call TLab_Stop(DNS_ERROR_OPTION)
                else
                    bound_s(1:inb_scal)%min = dummy(1:inb_scal)
                end if
            end if
        end if

        bound_s%max = 1.0_wp; inb_scal_local1 = MAX_VARS
        if (any(bound_s(:)%active)) then
            call ScanFile_Char(bakfile, inifile, 'Control', 'MaxScalar', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call LIST_REAL(sRes, inb_scal_local1, dummy)
                if (inb_scal_local1 /= inb_scal) then ! Consistency check
                    call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'MaxScalar size does not match inb_scal.')
                    call TLab_Stop(DNS_ERROR_OPTION)
                else
                    bound_s(1:inb_scal)%max = dummy(1:inb_scal)
                end if
            end if
        end if

        ! ###################################################################
        ! ###################################################################
        ! if (bound_p%min < 0.0_wp) bound_p%min = pbg%mean*1.0e-6_wp
        ! if (bound_p%max < 0.0_wp) bound_p%max = pbg%mean*1.0e6_wp
        ! if (bound_r%min < 0.0_wp) bound_r%min = rbg%mean*1.0e-6_wp
        ! if (bound_r%max < 0.0_wp) bound_r%max = rbg%mean*1.0e6_wp

        logs_data(1) = 0; obs_data(1) = 0 ! Status
        call DNS_BOUNDS_CONTROL()
        call DNS_OBS_CONTROL()
        ! call DNS_BOUNDS_LIMIT() !RH: schould be obsolete here

        ! call DNS_LOGS_PATH_INITIALIZE()
        call DNS_LOGS_INITIALIZE()
        call DNS_LOGS()
        if (dns_obs_log /= OBS_TYPE_NONE) then
            call DNS_OBS_INITIALIZE()
            call DNS_OBS()
        end if

        return
    end subroutine DNS_Control_Initialize

    !########################################################################
    ! Limit min/max values of fields, if required
    !########################################################################
    subroutine DNS_BOUNDS_LIMIT()
        use TLab_Memory, only: inb_scal
        use TLab_Arrays

        ! -------------------------------------------------------------------
        integer(wi) is

        ! ###################################################################
        do is = 1, inb_scal
            if (bound_s(is)%active) then
                s(:, is) = min(max(s(:, is), bound_s(is)%min), bound_s(is)%max)
            end if
        end do

        if (bound_r%active) then
            q(:, 5) = min(max(q(:, 5), bound_r%min), bound_r%max)
        end if

        if (bound_p%active) then
            q(:, 6) = min(max(q(:, 6), bound_p%min), bound_p%max)
        end if

        return
    end subroutine DNS_BOUNDS_LIMIT

    !########################################################################
    !########################################################################
    subroutine DNS_BOUNDS_CONTROL()
        use TLab_Constants, only: efile, lfile, fmt_r
        use NavierStokes, only: nse_eqns, DNS_EQNS_COMPRESSIBLE, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Arrays
        use TLab_WorkFlow, only: TLab_Write_ASCII
        use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_OutPlace
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
        use TLabMPI_VARS, only: ims_time_min, ims_err
#endif
        use FI_VECTORCALCULUS

        ! -------------------------------------------------------------------
        integer(wi) idummy(3)
        real(wp) dummy
#ifdef USE_MPI
#else
        integer wall_time_loc, int_dummy
#endif
        character*128 line
        character*32 str

        ! Pointers to existing allocated space
        real(wp), dimension(:, :, :), pointer :: loc_max

        ! ###################################################################
        ! Check wall time bounds - maximum runtime
#ifdef USE_MPI
        wall_time = MPI_WTIME() - ims_time_min
        call MPI_BCast(wall_time, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#else
        ! call ETIME(tdummy, wall_time_loc)
        call system_clock(wall_time_loc, int_dummy)
        wall_time = real(wall_time_loc - start_clock)/int_dummy
#endif
        ! ###################################################################
        ! Compressible flow
        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_COMPRESSIBLE)

#define p_min_loc logs_data(5)
#define p_max_loc logs_data(6)
#define r_min_loc logs_data(7)
#define r_max_loc logs_data(8)

            ! Check density
            call MINMAX(imax, jmax, kmax, q(:, 5), r_min_loc, r_max_loc)
            if (r_min_loc < bound_r%min .or. r_max_loc > bound_r%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Density out of bounds.')
                logs_data(1) = DNS_ERROR_NEGDENS
            end if

            ! Check pressure
            call MINMAX(imax, jmax, kmax, q(:, 6), p_min_loc, p_max_loc)
            if (p_min_loc < bound_p%min .or. p_max_loc > bound_p%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Pressure out of bounds.')
                logs_data(1) = DNS_ERROR_NEGPRESS
            end if

        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, q(1, 1), txc(1, 3))
                call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, q(1, 2), txc(1, 4))
                call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, q(1, 3), txc(1, 5))
                call FI_INVARIANT_P(imax, jmax, kmax, txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 1), txc(1, 2))
            else
                call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2))
            end if

#define d_max_loc logs_data(11)
#define d_min_loc logs_data(10)

            call MINMAX(imax, jmax, kmax, txc(1, 1), d_max_loc, d_min_loc)
            d_min_loc = -d_min_loc; d_max_loc = -d_max_loc

            if (max(abs(d_min_loc), abs(d_min_loc)) > bound_d%max) then
                call TLab_Write_ASCII(efile, 'DNS_CONTROL. Dilatation out of bounds.')
                logs_data(1) = DNS_ERROR_DILATATION

                ! Locating the points where the maximum dilatation occurs
                wrk3d = -txc(:, 1)
                loc_max(1:imax, 1:jmax, 1:kmax) => wrk3d(1:imax*jmax*kmax)

                dummy = maxval(wrk3d)
                if (abs(dummy) > bound_d%max) then
                    idummy = maxloc(loc_max)
                    write (str, fmt_r) dummy; line = 'Maximum dilatation '//trim(adjustl(str))
#ifdef USE_MPI
                    idummy(1) = idummy(1) + ims_offset_i
                    idummy(3) = idummy(3) + ims_offset_k
#endif
                    write (str, *) idummy(1); line = trim(adjustl(line))//' at grid node '//trim(adjustl(str))
                    write (str, *) idummy(2); line = trim(adjustl(line))//':'//trim(adjustl(str))
                    write (str, *) idummy(3); line = trim(adjustl(line))//':'//trim(adjustl(str))//'.'
                    call TLab_Write_ASCII(lfile, line, .true.)
                end if

                dummy = minval(wrk3d)
                if (abs(dummy) > bound_d%max) then
                    idummy = minloc(loc_max)
                    write (str, fmt_r) dummy; line = 'Minimum dilatation '//trim(adjustl(str))
#ifdef USE_MPI
                    idummy(1) = idummy(1) + ims_offset_i
                    idummy(3) = idummy(3) + ims_offset_k
#endif
                    write (str, *) idummy(1); line = trim(adjustl(line))//' at grid node '//trim(adjustl(str))
                    write (str, *) idummy(2); line = trim(adjustl(line))//':'//trim(adjustl(str))
                    write (str, *) idummy(3); line = trim(adjustl(line))//':'//trim(adjustl(str))//'.'
                    call TLab_Write_ASCII(lfile, line, .true.)
                end if

            end if

        end select

        return

    end subroutine DNS_BOUNDS_CONTROL
    !########################################################################
    !########################################################################
    subroutine DNS_OBS_CONTROL()
        use TLab_Memory, only: imax, jmax, kmax, inb_scal
        use TLab_WorkFlow, only: scal_on
        use TLab_Grid, only: z
        use FI_VORTICITY_EQN, only: FI_VORTICITY
        use TLab_Arrays
        use Averages, only: AVG_IK_V
        use Integration, only: Int_Simpson

        integer(wi) :: ip, is

        ! -------------------------------------------------------------------

#define ubulk    obs_data(2)
#define wbulk    obs_data(3)
#define uy1      obs_data(4)
#define wy1      obs_data(5)
#define alpha_1  obs_data(6)
#define alpha_ny obs_data(7)
#define int_ent  obs_data(8)

        ip = 8

        select case (dns_obs_log)

        case (OBS_TYPE_EKMAN)
            ! ubulk, wbulk
            call AVG_IK_V(imax, jmax, kmax, q(1, 1), wrk1d(:, 1), wrk1d(:, 2))
            call AVG_IK_V(imax, jmax, kmax, q(1, 3), wrk1d(:, 3), wrk1d(:, 4))
            ubulk = (1.0_wp/z%scale)*Int_Simpson(wrk1d(1:kmax, 1), z%nodes(1:kmax))
            wbulk = (1.0_wp/z%scale)*Int_Simpson(wrk1d(1:kmax, 3), z%nodes(1:kmax))

            ! dudy(1), dwdy(1) approximation
            uy1 = wrk1d(2, 1)/z%nodes(2)
            wy1 = wrk1d(2, 3)/z%nodes(2)

            ! turning angles (in degrees)
            alpha_1 = ATAN2D(wy1, uy1)
            alpha_ny = ATAN2D(wrk1d(kmax, 3), wrk1d(kmax, 1))

            ! integrated entstrophy
            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
            call AVG_IK_V(imax, jmax, kmax, txc(1, 1), wrk1d(:, 1), wrk1d(:, 2))
            int_ent = (1.0_wp/z%scale)*Int_Simpson(wrk1d(1:kmax, 1), z%nodes(1:kmax))

            if (scal_on) then
                do is = 1, inb_scal
                    call AVG_IK_V(imax, jmax, kmax, s(1, is), wrk1d(:, 1), wrk1d(:, 2))
                    obs_data(ip + is) = (wrk1d(2, 1) - wrk1d(1, 1))/z%nodes(2)
                end do
            end if

        end select

        return

    end subroutine DNS_OBS_CONTROL
    !########################################################################
    !# Initialize path to write dns.out & tlab.logs
    !########################################################################
    ! subroutine DNS_LOGS_PATH_INITIALIZE()

    !     integer env_status, path_len

    !     call get_environment_variable("DNS_LOGGER_PATH", logger_path, path_len, env_status, .true.)

    !     select case (env_status)
    !     case (-1)
    !         call TLab_Write_ASCII(efile, "DNS_START. The environment variable  is too long and cannot be handled in the foreseen array.")
    !         call TLab_Stop(DNS_ERROR_OPTION)
    !     case (0)
    !         if (.not. logger_path(path_len:path_len) == '/') then
    !             logger_path = trim(adjustl(logger_path))//'/'
    !         end if

    !     case (1:)
    !         logger_path = trim(adjustl(''))

    !     end select

    ! end subroutine DNS_LOGS_PATH_INITIALIZE

    !########################################################################
    ! Create headers or dns.out file
    !
    !# logs_data01 State (>0 if error)
    !#
    !# logs_data02 Maximum CFL number
    !# logs_data03 Maximum diffusion number
    !# logs_data04 Maximum source number
    !#
    !# logs_data05 Minimum pressure
    !# logs_data06 Maximum pressure
    !# logs_data07 Minimum density
    !# logs_data08 Maximum density
    !# logs_data09 NEWTONRAPHSON_ERROR
    !#
    !# logs_data10 Minimum dilatation
    !# logs_data11 Maximum dilatation
    !########################################################################
    subroutine DNS_LOGS_INITIALIZE()
        use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
        use Microphysics

        integer ip
        character(len=256) line1

        ! ofile = trim(adjustl(logger_path))//trim(adjustl(ofile_base))
        ! ofile = trim(adjustl(ofile))

        line1 = '#'; ip = 1
        line1 = line1(1:ip)//' '//' Itn.'; ip = ip + 1 + 7
        line1 = line1(1:ip)//' '//' time'; ip = ip + 1 + 13
        line1 = line1(1:ip)//' '//' dt'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' CFL#'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' D#'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' visc'; ip = ip + 1 + 10

        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            line1 = line1(1:ip)//' '//' DilMin'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' DilMax'; ip = ip + 1 + 13

        case (DNS_EQNS_COMPRESSIBLE)
            line1 = line1(1:ip)//' '//' PMin'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' PMax'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' RMin'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' RMax'; ip = ip + 1 + 10

        end select

        if (imixture == MIXT_TYPE_AIRWATER .and. (evaporationProps%type == TYPE_EVA_EQUILIBRIUM .or. evaporationProps%type == TYPE_EVA_QSCALC_IMPL .or. evaporationProps%type == TYPE_EVA_QSCALC_SEMIIMPL)) then
            line1 = line1(1:ip)//' '//' NewtonRs'; ip = ip + 1 + 10
        end if

        line1 = line1(1:ip - 1)//'#'
        call TLab_Write_ASCII(ofile, repeat('#', len_trim(line1)))
        call TLab_Write_ASCII(ofile, trim(adjustl(line1)))
        call TLab_Write_ASCII(ofile, repeat('#', len_trim(line1)))

    end subroutine DNS_LOGS_INITIALIZE

    !########################################################################

    subroutine DNS_LOGS()
        use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
        use Thermo_Anelastic, only: NEWTONRAPHSON_ERROR
        use Microphysics
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_err
        real(wp) dummy
#endif

        integer ip
        character(len=256) line1, line2

        write (line1, 100) int(logs_data(1)), itime, rtime, logs_dtime, (logs_data(ip), ip=2, 3), visc
100     format((1x, I1), (1x, I7), (1x, E13.6), 4(1x, E10.3))

        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            write (line2, 200) logs_data(10), logs_data(11)
200         format(2(1x, E13.6))
            line1 = trim(line1)//trim(line2)

        case (DNS_EQNS_COMPRESSIBLE)
            write (line2, 300) (logs_data(ip), ip=5, 8)
300         format(4(1x, E10.3))
            line1 = trim(line1)//trim(line2)

        end select

        if (imixture == MIXT_TYPE_AIRWATER .and. (evaporationProps%type == TYPE_EVA_EQUILIBRIUM .or. evaporationProps%type == TYPE_EVA_QSCALC_IMPL .or. evaporationProps%type == TYPE_EVA_QSCALC_SEMIIMPL)) then
#ifdef USE_MPI
            call MPI_ALLREDUCE(NEWTONRAPHSON_ERROR, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
            NEWTONRAPHSON_ERROR = dummy
#endif
            write (line2, 400) NEWTONRAPHSON_ERROR
400         format(1(1x, E10.3))
            line1 = trim(line1)//trim(line2)
        end if

        call TLab_Write_ASCII(ofile, trim(adjustl(line1)))

    end subroutine DNS_LOGS

    !########################################################################
    !# Create headers or dns.obs file
    !########################################################################
    subroutine DNS_OBS_INITIALIZE()

        implicit none

        integer(wi) :: ip, is
        character(len=256) :: line1, str

        ! vfile = trim(adjustl(logger_path))//trim(adjustl(vfile_base))
        ! vfile = trim(adjustl(vfile))

        line1 = '#'; ip = 1
        line1 = line1(1:ip)//' '//' Itn.'; ip = ip + 1 + 7
        line1 = line1(1:ip)//' '//' time'; ip = ip + 1 + 13

        select case (dns_obs_log)
        case (OBS_TYPE_EKMAN)
            line1 = line1(1:ip)//' '//' u_bulk'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' w_bulk'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' u_y(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' w_y(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' alpha(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' alpha(ny)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' entstrophy'; ip = ip + 1 + 13
            if (scal_on) then
                do is = 1, inb_scal
                    write (str, *) is
                    line1 = line1(1:ip)//' '//' s'//trim(adjustl(str))//'_y(1)'; ip = ip + 1 + 13
                end do
            end if
        end select

        line1 = line1(1:ip - 1)//'#'
        call TLab_Write_ASCII(vfile, repeat('#', len_trim(line1)))
        call TLab_Write_ASCII(vfile, trim(adjustl(line1)))
        call TLab_Write_ASCII(vfile, repeat('#', len_trim(line1)))

    end subroutine DNS_OBS_INITIALIZE

    !########################################################################

    subroutine DNS_OBS()

        implicit none

        integer(wi) :: ip, is
        character(len=256) :: line1, line2

        write (line1, 100) int(obs_data(1)), itime, rtime
100     format((1x, I1), (1x, I7), (1x, E13.6))

        select case (dns_obs_log)
        case (OBS_TYPE_EKMAN)
            write (line2, 200) (obs_data(ip), ip=2, 8)
200         format(7(1x, E13.6))
            line1 = trim(line1)//trim(line2)
            if (scal_on) then
                do is = 1, inb_scal
                    write (line2, 300) obs_data(8 + is)
300                 format(1x, E13.6)
                    line1 = trim(line1)//trim(line2)
                end do
            end if
        end select

        call TLab_Write_ASCII(vfile, trim(adjustl(line1)))

    end subroutine DNS_OBS

end module DNS_Control
