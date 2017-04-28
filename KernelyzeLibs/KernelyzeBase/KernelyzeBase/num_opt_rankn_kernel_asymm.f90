! num_opt_rankn_kernel_asymm.f90
!
! Copyright (c) 2016 by Kernelyze LLC
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
! created on: 2016-07-05
! updated on: 2016-07-05
!
! This module contains a constructor subroutine
! that builds a numerically-optimal rank-$n$
! approximating kernel (optimality being in the
! supremum norm).
!
! Unlike num_opt_rankn_kernel, this procedure can handle 
! (nondegenerate) totally positive kernels that may be
! asymmetric (or symmetric).
  
module num_opt_rankn_kernel_asymm_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use kernel_mod, only : kernel
  use kernel_nsect_mod, only : kernel_rankn
  use num_opt_rankn_params_mod, only : num_opt_rankn_params
  use num_opt_rankn_kernel_of_params_mod, only : &
      num_opt_rankn_kernel_of_params

  implicit none

  contains
  
  subroutine num_opt_rankn_kernel_asymm( &
      kernel_to_approx, &
      rank, &
      max_iter, &
      toler, &
      rankn_kernel, &
      borsuk_lb, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    integer, intent(in)           :: rank
    integer, intent(in)           :: max_iter
    real(wp), intent(in)          :: toler
    class(kernel_rankn), &
         allocatable, intent(out) :: rankn_kernel
    real(wp), intent(out)         :: borsuk_lb
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'num_opt_rankn_kernel_asymm: '
    ! Local variables
    real(wp)  :: rho_vec(rank + 1)
    real(wp)  :: gamma_vec(rank + 1)
    real(wp)  :: v_mat(rank + 1, rank)
    real(wp)  :: w_mat(rank + 1, rank)
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Compute the numerically-optimal parameters
    call num_opt_rankn_params( &
        kernel_to_approx = kernel_to_approx, &
        rank = rank, &
        max_iter = max_iter, &
        toler = toler, &
        rho_vec = rho_vec, &
        gamma_vec = gamma_vec, &
        v_mat = v_mat, &
        w_mat = w_mat, &
        borsuk_lb = borsuk_lb, &
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
    ! Use the numerically-optimal parameters to build the
    ! rank-$n$ numerically-optimal approximating kernel
    call num_opt_rankn_kernel_of_params( &
        kernel_to_approx = kernel_to_approx, &
        rho_vec = rho_vec, &
        gamma_vec = gamma_vec, &
        v_mat = v_mat, &
        w_mat = w_mat, &
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
  end subroutine num_opt_rankn_kernel_asymm

end module num_opt_rankn_kernel_asymm_mod