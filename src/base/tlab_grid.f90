#include "tlab_error.h"

module TLab_Grid
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: TLab_Grid_Read
    public :: TLab_Grid_Write

    public :: grid_dt
    public :: grid, x, y, z
    public :: subGrid, subX, subY, subZ

    ! -----------------------------------------------------------------------
    type :: grid_dt
        ! sequence
        character*8 name
        integer(wi) size
        logical :: uniform = .false.
        logical :: periodic = .false.
        real(wp) scale
        real(wp), allocatable :: nodes(:)
    end type

    type(grid_dt), target :: grid(3)
    type(grid_dt), pointer :: x => grid(1), y => grid(2), z => grid(3)

    type :: subgrid_dt
        type(grid_dt), pointer :: parent
        integer :: offset
        integer :: size
    contains
        procedure :: initialize => sub_grid_initialize
    end type

    type(subgrid_dt), target :: subGrid(3)
    type(subgrid_dt), pointer :: subX => subGrid(1), subY => subGrid(2), subZ => subGrid(3)

contains
    !########################################################################
    !########################################################################
    subroutine sub_grid_initialize(self, grid)
        class(subgrid_dt), intent(inout) :: self
        type(grid_dt), intent(in), target :: grid

        self%parent => grid
        self%offset = 0
        self%size = grid%size

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine TLab_Grid_Read(name, x, y, z, sizes)
        character*(*) name
        type(grid_dt), intent(inout) :: x, y, z
        integer(wi), intent(in), optional :: sizes(3)

        ! -----------------------------------------------------------------------
        character*(32) line

        ! #######################################################################
        open (50, file=name, status='old', form='unformatted')
        rewind (50)

        ! -----------------------------------------------------------------------
        read (50) x%size, y%size, z%size

        if (present(sizes)) then        ! check
            if (any([x%size, y%size, z%size] /= sizes)) then
                close (50)
                write (line, 100) x%size, y%size, z%size
                call TLab_Write_ASCII(efile, __FILE__//'. Dimensions ('//trim(line)//') unmatched.')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if
        end if

        read (50) x%scale, y%scale, z%scale

        if (allocated(x%nodes)) deallocate (x%nodes)
        if (allocated(y%nodes)) deallocate (y%nodes)
        if (allocated(z%nodes)) deallocate (z%nodes)
        allocate (x%nodes(x%size), y%nodes(y%size), z%nodes(z%size))

        read (50) x%nodes(:)
        read (50) y%nodes(:)
        read (50) z%nodes(:)

        ! -----------------------------------------------------------------------
        close (50)

        return

100     format(I5, ',', I5, ',', I5)

    end subroutine TLab_Grid_Read

!########################################################################
!########################################################################
    subroutine TLab_Grid_Write(name, x, y, z)
        character*(*) name
        type(grid_dt), intent(in) :: x, y, z

        !########################################################################
        open (unit=51, file=name, form='unformatted', status='unknown')

        write (51) x%size, y%size, z%size
        write (51) x%scale, y%scale, z%scale

        write (51) x%nodes(1:x%size)
        write (51) y%nodes(1:y%size)
        write (51) z%nodes(1:z%size)

        close (51)

        return
    end subroutine TLab_Grid_Write

end module TLab_Grid
