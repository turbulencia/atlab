module TLabMPI_VARS
    use TLab_Constants, only: wp, wi
    use mpi_f08, only: MPI_Comm, MPI_Datatype
    implicit none

    type(MPI_Comm) :: ims_comm_xy                               ! Plane communicators
    type(MPI_Comm) :: ims_comm_x, ims_comm_y                    ! Directional communicators
    type(MPI_Comm) :: ims_comm_xy_aux
    type(MPI_Comm) :: ims_comm_x_aux, ims_comm_y_aux

    integer :: ims_pro                                          ! Local PE number in global communicator
    integer :: ims_pro_i, ims_pro_j, ims_pro_k                  ! Local PE number in directional communicators
    integer :: ims_npro                                         ! Total number of PE
    integer :: ims_npro_i, ims_npro_j, ims_npro_k               ! Number of PEs in directional communicators
    integer :: ims_offset_i, ims_offset_j, ims_offset_k         ! Grid offsets

    integer :: ims_err, ims_tag

    real(wp) :: ims_time_min, ims_time_max, ims_time_trans      ! Profiling
    integer(wi) :: ims_bcs_imax, ims_bcs_jmax

    type(MPI_Datatype) :: TLAB_MPI_REAL_TYPE                    ! MPI Type control

    type :: mpi_grid_dt
        type(MPI_Comm) :: comm
        integer :: num_processors
        integer :: rank
    end type

end module TLabMPI_VARS
