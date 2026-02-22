#include "tlab_error.h"

module TLabMPI_PROCS
    use mpi_f08
    use TLab_Constants, only: wp, dp, sp, wi, lfile, efile
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS
    implicit none
    private

    public :: TLabMPI_Initialize
    public :: TLabMPI_Halos_X, TLabMPI_Halos_Y
    public :: TLabMPI_Panic

contains
    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) dims(2), coord(2)
        logical period(2), remain_dims(2), reorder

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes!, line
        character*64 lstr
        integer(wi) nsize_total

        ! #######################################################################
        ! Domain decomposition in parallel mode
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Grid'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Imax(*)=<imax_proc>')
        call TLab_Write_ASCII(bakfile, '#Jmax(*)=<jmax_proc>')

        if (ims_npro > 1) then
            nsize_total = jmax
            write (lstr, *) nsize_total
            call ScanFile_Int(bakfile, inifile, block, 'Jmax(*)', trim(adjustl(lstr)), jmax)
            if (jmax > 0 .and. mod(nsize_total, jmax) == 0) then
                ims_npro_j = nsize_total/jmax
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Input jmax incorrect')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if

            nsize_total = imax
            write (lstr, *) nsize_total
            call ScanFile_Int(bakfile, inifile, block, 'Imax(*)', trim(adjustl(lstr)), imax)
            if (imax > 0 .and. mod(nsize_total, imax) == 0) then
                ims_npro_i = nsize_total/imax
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Input imax incorrect')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if

            ims_npro_k = 1

            ! consistency check
            if (ims_npro_i*ims_npro_j == ims_npro) then
                write (lstr, *) ims_npro_i; write (sRes, *) ims_npro_j
                lstr = trim(adjustl(lstr))//'x'//trim(adjustl(sRes))
                call TLab_Write_ASCII(lfile, 'Initializing domain partition '//trim(adjustl(lstr)))
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Inconsistency in total number of PEs')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if

        end if

        ! #######################################################################
        call TLab_Write_ASCII(lfile, 'Creating MPI communicators.')

        ! the first index in the grid corresponds to j, the second to i
        dims(1) = ims_npro_j; dims(2) = ims_npro_i; period = .true.; reorder = .false.
        ! dims(1) = ims_npro_i; dims(2) = ims_npro_j; period = .true.; reorder = .false.
        call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xy, ims_err)

        call MPI_CART_COORDS(ims_comm_xy, ims_pro, 2, coord, ims_err)
        ims_pro_j = coord(1); ims_pro_i = coord(2)      ! starting at 0
        ! ims_pro_j = coord(2); ims_pro_i = coord(1)
        !
        ! equivalent to:
        ! ims_pro_i = mod(ims_pro, ims_npro_i) ! Starting at 0
        ! ims_pro_j = ims_pro/ims_npro_i  ! Starting at 0
        ! to revert them:

        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_CART_SUB(ims_comm_xy, remain_dims, ims_comm_x, ims_err)
        ! call MPI_CART_SUB(ims_comm_xy, remain_dims, ims_comm_y, ims_err)

        remain_dims(1) = .true.; remain_dims(2) = .false.
        call MPI_CART_SUB(ims_comm_xy, remain_dims, ims_comm_y, ims_err)
        ! call MPI_CART_SUB(ims_comm_xy, remain_dims, ims_comm_x, ims_err)

        ! ip = ims_pro
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_x, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along X', id
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_y, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along Z', id

        ims_offset_i = ims_pro_i*imax       ! local offset in grid points
        ims_offset_j = ims_pro_j*jmax
        ims_offset_k = 0

        ! -----------------------------------------------------------------------
        ! Control of MPI type
        select case (wp)
        case (dp)
            TLAB_MPI_REAL_TYPE = MPI_REAL8
        case (sp)
            TLAB_MPI_REAL_TYPE = MPI_REAL4
        end select

        return
    end subroutine TLabMPI_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Panic(location, mpi_error_code)
        character(len=*), intent(in) :: location
        integer, intent(in) :: mpi_error_code

        !##############################
        character error_string*1024
        integer error_local, error_len

        call MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
        call TLab_Write_ASCII(efile, 'MPI-ERROR: Source file'//trim(adjustl(LOCATION)), .true.)
        call TLab_Write_ASCII(efile, error_string, .true.)

        call TLab_Stop(mpi_error_code)
        ! Not supposed to return from this subroutine

    end subroutine TLabMPI_Panic

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_X(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(:)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(:)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i)
        source = mod(ims_pro_i + 1, ims_npro_i)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          ims_comm_x, MPI_STATUS_IGNORE, ims_err)

        ! pass to following processor
        dest = mod(ims_pro_i + 1, ims_npro_i)
        source = mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          ims_comm_x, MPI_STATUS_IGNORE, ims_err)

        return
    end subroutine TLabMPI_Halos_X

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_Y(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(:)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(:)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(ims_pro_j - 1 + ims_npro_j, ims_npro_j)
        source = mod(ims_pro_j + 1, ims_npro_j)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          ims_comm_y, MPI_STATUS_IGNORE, ims_err)

        ! pass to following processor
        dest = mod(ims_pro_j + 1, ims_npro_j)
        source = mod(ims_pro_j - 1 + ims_npro_j, ims_npro_j)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          ims_comm_y, MPI_STATUS_IGNORE, ims_err)

        return
    end subroutine TLabMPI_Halos_Y

end module TLabMPI_PROCS
