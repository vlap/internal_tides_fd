module err_manager

     use precisions, only: dp

     implicit none

        integer, parameter               :: nf90_NoErr = 0

        contains
!********************************************************************
!********************************************************************
!     Error Management
!********************************************************************
!********************************************************************

     subroutine handle_nc_err(err_code, section)

     integer, intent(in)                    ::     err_code
     character (len = *), intent(in) :: section

     if (err_code /= nf90_NoErr) then
           write(*,*) 'Something went wrong in section ', section
    end if

     end subroutine

!********************************************************************

     subroutine handle_av_err(message, weight1, weight2)

     character (len = *), intent(in) :: message
     real(dp)     ::     weight1, weight2

           print *,  message, ", using weights: ", weight1, weight2

     end subroutine

end

