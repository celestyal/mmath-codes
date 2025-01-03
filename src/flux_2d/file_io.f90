module file_io

    use types

    implicit none

contains
    subroutine read_grid(g, identifier)
        character(len=*), intent(in) :: identifier
        type(grid), intent(out) :: g

        namelist /grid/ g
        open(unit=1, file=("exports/" // trim(adjustl(identifier)) // "/parameters"))
        read(unit=1, nml=grid)
        close(unit=1)
    end subroutine read_grid


    subroutine read_debug(d, identifier)
        character(len=*), intent(in) :: identifier
        type(debug), intent(out) :: d

        namelist /debug/ d

        open(unit=1, file=("exports/" // trim(adjustl(identifier)) // "/parameters"))
        read(unit=1, nml=debug)
        close(unit=1)
    end subroutine read_debug


    subroutine read_profile(p, identifier)
        character(len=*), intent(in) :: identifier
        type(profile), intent(out) :: p

        namelist /profile/ p

        open(unit=1, file=("exports/" // trim(adjustl(identifier)) // "/parameters"))
        read(unit=1, nml=profile)
        close(unit=1)
    end subroutine read_profile


    subroutine export_data(identifier, data, filename, mode)
        real, intent(in) :: data(:,:)
        character(len=*), intent(in):: filename, mode
        character(512), intent(in) :: identifier

        open(unit=2, file=("exports/" // trim(adjustl(identifier)) // "/" // filename), action="write", position=mode, form="unformatted")
        write(2) data
        close(unit=2)
    end subroutine export_data


    subroutine save_parameters(identifier, d,g,p)
        type(debug), intent(in) :: d
        type(grid), intent(in) :: g
        type(profile), intent(in) :: p
        character(512), intent(in) :: identifier
        character(len=512) :: filepath

        namelist /debug/ d
        namelist /grid/ g
        namelist /profile/ p

        filepath = "exports/" // trim(adjustl(identifier)) // "/parameters"

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
    subroutine export_for_idl(identifier, g)
        character(512), intent(in) :: identifier
        type(grid), intent(in) :: g

        open(unit=2, file=("exports/" // trim(adjustl(identifier)) // "/simulation"), action="write", position="asis", form="formatted")
        write(2, *) g%nx
        write(2, *) g%ny
        write(2, *) g%xmax
        write(2, *) g%xmin
        write(2, *) g%ymax
        write(2, *) g%ymin
        close(unit=2)
    end subroutine export_for_idl


    ! append current timestamp to the system file, so that saves can be timestamped.
    subroutine save_timestamp(identifier, time)
        character(512), intent(in) :: identifier
        real, intent(in) :: time

        open(unit=2, file=("exports/" // trim(adjustl(identifier)) // "/simulation"), action="write", position="append", form="formatted")
        write(2, *) time
        close(unit=2)
    end subroutine save_timestamp


    function read_identifier() result(identifier)
        character(512) :: identifier
        logical :: exists

        inquire(file=".simulation_id", exist=exists)
        if (.not. exists) stop " The file './.simulation_id' must exist and contain a non-empty new or existing simulation identifier string"

        open(unit=2, file=".simulation_id", action="read", status="unknown")
        read(2, *), identifier
        close(2)
    end function read_identifier
end module file_io
