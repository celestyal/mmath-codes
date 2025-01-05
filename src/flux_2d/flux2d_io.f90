module flux2d_io

    use flux2d_types

    implicit none

contains
    subroutine read_grid(g, id)
        character(len=*), intent(in) :: id
        type(grid), intent(out) :: g

        namelist /grid/ g
        open(unit=1, file=("exports/" // trim(adjustl(id)) // "/parameters"))
        read(unit=1, nml=grid)
        close(unit=1)
    end subroutine read_grid


    subroutine read_debug(d, id)
        character(len=*), intent(in) :: id
        type(debug), intent(out) :: d

        namelist /debug/ d

        open(unit=1, file=("exports/" // trim(adjustl(id)) // "/parameters"))
        read(unit=1, nml=debug)
        close(unit=1)
    end subroutine read_debug


    subroutine read_profile(p, id)
        character(len=*), intent(in) :: id
        type(profile), intent(out) :: p

        namelist /profile/ p

        open(unit=1, file=("exports/" // trim(adjustl(id)) // "/parameters"))
        read(unit=1, nml=profile)
        close(unit=1)
    end subroutine read_profile


    subroutine export_data(id, data, filename, mode)
        real, intent(in) :: data(:,:)
        character(len=*), intent(in):: filename, mode
        character(512), intent(in) :: id

        open(unit=2, file=("exports/" // trim(adjustl(id)) // "/" // filename), action="write", position=mode, form="unformatted")
        write(2) data
        close(unit=2)
    end subroutine export_data


    subroutine save_parameters(id, d,g,p)
        type(debug), intent(in) :: d
        type(grid), intent(in) :: g
        type(profile), intent(in) :: p
        character(512), intent(in) :: id
        character(len=512) :: filepath

        namelist /debug/ d
        namelist /grid/ g
        namelist /profile/ p

        filepath = "exports/" // trim(adjustl(id)) // "/parameters"

        open(unit=1, file=trim(filepath), action="write", position="asis", form="formatted")
        write(unit=1, nml=debug)
        close(unit=1)

        open(unit=1, file=trim(filepath), action="write", position="append", form="formatted")
        write(unit=1, nml=grid)
        close(unit=1)

        open(unit=1, file=trim(filepath), action="write", position="append", form="formatted")
        write(unit=1, nml=profile)
        close(unit=1)
    end subroutine save_parameters


    ! write the information needed for idl axis and timestamps
    subroutine export_for_idl(id, g)
        character(512), intent(in) :: id
        type(grid), intent(in) :: g

        open(unit=2, file=("exports/" // trim(adjustl(id)) // "/simulation"), action="write", position="asis", form="formatted")
        write(2, *) g%nx
        write(2, *) g%ny
        write(2, *) g%xmax
        write(2, *) g%xmin
        write(2, *) g%ymax
        write(2, *) g%ymin
        close(unit=2)
    end subroutine export_for_idl


    ! append current timestamp to the system file, so that saves can be timestamped.
    subroutine save_timestamp(id, time)
        character(512), intent(in) :: id
        real, intent(in) :: time

        open(unit=2, file=("exports/" // trim(adjustl(id)) // "/simulation"), action="write", position="append", form="formatted")
        write(2, *) time
        close(unit=2)
    end subroutine save_timestamp


    function read_id() result(id)
        character(512) :: id
        logical :: exists

        inquire(file=".simulation_id", exist=exists)
        if (.not. exists) stop " The file './.simulation_id' must exist and contain a non-empty new or existing simulation id string"

        open(unit=2, file=".simulation_id", action="read", status="unknown")
        read(2, *), id
        close(2)
    end function read_id
end module flux2d_io
