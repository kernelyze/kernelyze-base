! num_opt_rankn_kernel_of_params.f90
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
! updated on: 2016-11-05 (Added code to account for analytical integral
!             calculations that may be required later [this uses
!             calculate_integral_kernel_svar_impl])
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
!
! This module contains a constructor subroutine
! that builds a numerically-optimal rank-$n$
! approximating kernel (optimality being in the
! supremum norm) given parameters for the kernel,
! vectors rho and gamma and matrices V and W.  Such
! parameters could, for example, come from num_opt_rankn_params.
!
! Unlike num_opt_rankn_kernel, this procedure can handle 
! (nondegenerate) totally positive kernels that may be
! asymmetric (or symmetric).
  
module num_opt_rankn_kernel_of_params_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use kernel_mod, only : kernel
  use kernel_nsect_mod, only : kernel_rankn
  use compfunc_kernel_mod, only : compfunc_kernel
  
  implicit none

  contains
  
  subroutine num_opt_rankn_kernel_of_params( &
      kernel_to_approx, &
      rho_vec, &
      gamma_vec, &
      v_mat, &
      w_mat, &
      rankn_kernel, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    real(wp), intent(in)          :: rho_vec(:)
    real(wp), intent(in)          :: gamma_vec(size(rho_vec))
    real(wp), intent(in)          :: v_mat(size(rho_vec), size(rho_vec) - 1)
    real(wp), intent(in)          :: w_mat(size(rho_vec), size(rho_vec) - 1)
    class(kernel_rankn), allocatable, intent(out) :: rankn_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'num_opt_rankn_kernel_of_params: '
    ! Local variables
    type(compfunc_kernel)               :: compfuncs_x(size(rho_vec))
    type(compfunc_kernel)               :: compfuncs_y(size(rho_vec))
    integer                             :: j
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
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
    ! Allocate and build the component functions that are functions of x
    do j = 1, size(compfuncs_x)
      call compfuncs_x(j)%set_kernel(kernel_to_approx, &
          local_err_stat, local_err_msg)
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
      call compfuncs_x(j)%set_eval_pt(rho_vec(j))
      call compfuncs_x(j)%set_x_is_fixed(.false.)
    end do
    ! Allocate and build the component functions that are functions of y
    do j = 1, size(compfuncs_y)
      call compfuncs_y(j)%set_kernel(kernel_to_approx, &
          local_err_stat, local_err_msg)
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
      call compfuncs_y(j)%set_eval_pt(gamma_vec(j))
      call compfuncs_y(j)%set_x_is_fixed(.true.)
    end do
    call rankn_kernel%set_f_array(compfuncs_x, local_err_stat, local_err_msg)
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
    call rankn_kernel%set_g_array(compfuncs_y, local_err_stat, local_err_msg)
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
    call rankn_kernel%set_v_matrix(v_mat)
    call rankn_kernel%set_w_matrix(w_mat)
  end subroutine num_opt_rankn_kernel_of_params

end module num_opt_rankn_kernel_of_params_mod