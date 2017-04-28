! taylor_series_rankn_kernel.f90
!
! Copyright (c) 2016, 2017 by Kernelyze LLC
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
! created on: 2016-02-08
! updated on: 2016-02-08
! updated on: 2016-02-10 (added proper centering of monomials)
! updated on: 2016-02-13 (fixed indexing off-by-one for partials
!             and i-vs-j confusion of eigenval indexing in 
!             V matrix construction)
! updated on: 2016-07-05 (Accounted for kernel_rankn change that includes
!             a W matrix to permit kernel_rankn objects to accommodate
!             numerically-optimal rank-$n$ approximations to totally
!             positive kernels that are asymmetric -- note that the
!             procedure taylor_series_rankn_kernel only applies to *symmetric*
!             nondegenerate totally positive kernels, this change is
!             to give greater flexibility to kernel_rankn objects in
!             *other* contexts)
! updated on: 2016-09-30 (added check for negative eigenvalues to catch
!             errors)
! updated on: 2016-10-02 (switched from symmetric eigendecomp to singular
!             value decomposition to accommodate asymmetric kernels)
! updated on: 2016-04-20 (refactored into numerical procedure and
!             rank-$n$ kernel constructor)
! updated on: 2016-04-27 (more flexibility for finite differences)
!
! This module contains a constructor subroutine
! that builds a rank-$n$ two-variable Taylor
! series approximating kernel, given a kernel to 
! approximate and a permitted approximation rank, $n$.
  
module taylor_series_rankn_kernel_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use kernel_mod, only : kernel
  use kernel_nmonom_mod, only : kernel_rankn
  use taylor_series_rankn_params_mod, only : taylor_series_rankn_params
  use taylor_series_rankn_kernel_of_params_mod, only : &
      taylor_series_rankn_kernel_of_params
  
  implicit none

  contains
  
  subroutine taylor_series_rankn_kernel( &
      kernel_to_approx, &
      x_center, &
      y_center, &
      rank, &
      rankn_kernel, &
      err_stat, &
      err_msg, &
      fin_diff_delta)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    real(wp), intent(in)          :: x_center
    real(wp), intent(in)          :: y_center
    integer, intent(in)           :: rank
    class(kernel_rankn), allocatable, intent(out) :: rankn_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Optional finite-diff parameter
    real(wp), intent(in), optional  :: fin_diff_delta
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'taylor_series_rankn_kernel: '
    ! Local variables
    real(wp), dimension(rank , rank)  :: v_matrix
    real(wp), dimension(rank , rank)  :: w_matrix
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Obtain the V and W matrices (the factorization of the
    ! matrix of scaled partial derivatives of the kernel)
    if (present(fin_diff_delta)) then
      call taylor_series_rankn_params( &
          kernel_to_approx = kernel_to_approx, &
          x_center = x_center, &
          y_center = y_center, &
          rank = rank, &
          v_matrix = v_matrix, &
          w_matrix = w_matrix, &
          err_stat = local_err_stat, &
          err_msg = local_err_msg, &
          fin_diff_delta = fin_diff_delta)
    else
      call taylor_series_rankn_params( &
          kernel_to_approx = kernel_to_approx, &
          x_center = x_center, &
          y_center = y_center, &
          rank = rank, &
          v_matrix = v_matrix, &
          w_matrix = w_matrix, &
          err_stat = local_err_stat, &
          err_msg = local_err_msg)
    end if
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
    ! Using the V and W matrices, build the rank-$n$
    ! approximating kernel
    call taylor_series_rankn_kernel_of_params( &
        x_center = x_center, &
        y_center = y_center, &
        v_matrix = v_matrix, &
        w_matrix = w_matrix, &
        rankn_kernel = rankn_kernel, &
        err_stat = local_err_stat, &
        err_msg = local_err_msg)
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
  end subroutine taylor_series_rankn_kernel

end module taylor_series_rankn_kernel_mod
