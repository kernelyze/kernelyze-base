! taylor_series_rankn_kernel_of_params.f90
!
! Copyright (c) 2017 by Kernelyze LLC
! Author: Thomas A. Knox
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Affero General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Affero General Public License for more details.
! 
! You should have received a copy of the GNU Affero General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! created on: 2017-04-20 (based on earlier taylor_series_rankn_kernel)
! updated on: 2017-04-22 (remove unnecessary use of LAPACK module)
!
! This module contains a constructor subroutine
! that builds a rank-$n$ two-variable Taylor
! series approximating kernel using parameter
! matrices that factor the matrix of scaled partial
! derivatives.
  
module taylor_series_rankn_kernel_of_params_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use factorial_mod, only : factorial
  use kernel_mod, only : kernel
  use kernel_nmonom_mod, only : kernel_rankn
  use monomial_mod, only : monomial
  
  implicit none

  contains
  
  subroutine taylor_series_rankn_kernel_of_params( &
      x_center, &
      y_center, &
      v_matrix, &
      w_matrix, &
      rankn_kernel, &
      err_stat, &
      err_msg)
    ! Arguments
    real(wp), intent(in)  :: x_center
    real(wp), intent(in)  :: y_center
    real(wp), intent(in)  :: v_matrix( : , : )
    real(wp), intent(in)  :: w_matrix( : , : )
    class(kernel_rankn), allocatable, intent(out) :: rankn_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = &
        'taylor_series_rankn_kernel_of_params: '
    ! Local variables
    type(monomial)  :: x_compfuncs( size( v_matrix , 1 ) )
    type(monomial)  :: y_compfuncs( size( w_matrix , 1 ) )
    integer         :: j
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Be certain that V and W conform so that $V W^{T}$ makes sense
    if (size(v_matrix, 2) /= size(w_matrix, 2)) then
      ! Set the error message and error status if they were provided
      if (present(err_msg)) then
        err_msg = proc_name &
            // 'The V and W matrices must have the same number of columns.'
      end if
      if (present(err_stat)) then
        err_stat = -101
      end if
      return
    end if
    ! Allocate the rankn_kernel
    allocate(rankn_kernel, stat = local_err_stat, errmsg = local_err_msg)
    ! Handle any allocation problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    ! Set up the components of the rankn_kernel
    ! Allocate and build the component functions
    do j = 1, size(x_compfuncs)
      call x_compfuncs(j)%set_power(j - 1)
      call x_compfuncs(j)%set_center(x_center)
    end do
    call rankn_kernel%set_f_array(x_compfuncs, local_err_stat, local_err_msg)
    ! Handle any problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    do j = 1, size(y_compfuncs)
      call y_compfuncs(j)%set_power(j - 1)
      call y_compfuncs(j)%set_center(y_center)
    end do
    call rankn_kernel%set_g_array(y_compfuncs, local_err_stat, local_err_msg)
    ! Handle any problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    call rankn_kernel%set_v_matrix(v_matrix)
    call rankn_kernel%set_w_matrix(w_matrix)
  end subroutine taylor_series_rankn_kernel_of_params

end module taylor_series_rankn_kernel_of_params_mod