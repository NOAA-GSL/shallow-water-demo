module NetCDF_Utils

use netcdf

contains

  !------------------------------------------------------------------
  ! nc_check
  !
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  subroutine nc_check(istatus)

    integer, intent (in) :: istatus

    character(len=512) :: error_msg

    ! if no error, nothing to do here.  we are done.
    if( istatus == nf90_noerr) return

    error_msg = nf90_strerror(istatus)

    print *,error_msg
    stop 1

  end subroutine nc_check

end module NetCDF_Utils
