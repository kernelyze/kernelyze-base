! num_opt_rankn_params.f90
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
! created on: 2016-07-04 
! updated on: 2016-07-04
! updated on: 2016-07-05 (continue work toward first working draft)
!
! This module contains a procedure that computes
! all of the parameters that define a numerically-
! optimal rank-$n$ approximation to a kernel that
! is nondegenerate totally positive on its (rectangular)
! domain.  Unlike num_opt_rankn_kernel, the output here
! is composed of the parameters needed to evaluate the
! numerically-optimal rank-$n$ approximation; no
! approximate-kernel object is built.
!
! The kernel need ***not*** be symmetric, unlike the
! procedure num_opt_rankn_kernel.  It does, of course,
! need to be (nondegenerate) totally positive.
  
module num_opt_rankn_params_mod

  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use kernel_mod, only : kernel
  use worst_rho_and_gamma_mod, only : worst_rho_and_gamma
  use all_num_opt_mats_asymm_mod , only : all_num_opt_mats_asymm

  implicit none

  contains
  
  subroutine num_opt_rankn_params( &
      kernel_to_approx, &
      rank, &
      max_iter, &
      toler, &
      rho_vec, &
      gamma_vec, &
      v_mat, &
      w_mat, &
      borsuk_lb, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    integer, intent(in)           :: rank
    integer, intent(in)           :: max_iter
    real(wp), intent(in)          :: toler
    real(wp), intent(out)         :: rho_vec(rank + 1)
    real(wp), intent(out)         :: gamma_vec(rank + 1)
    real(wp), intent(out)         :: v_mat(rank + 1, rank)
    real(wp), intent(out)         :: w_mat(rank + 1, rank)
    real(wp), intent(out)         :: borsuk_lb
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'num_opt_rankn_params: '
    ! Local variables
    real(wp), dimension(rank + 1 , rank + 1)    :: a_mat
    real(wp), dimension(rank + 1 , rank + 1)    :: b_mat
    real(wp), dimension(rank + 1)               :: coeffs_rho
    real(wp), dimension(rank + 1)               :: coeffs_gamma
    real(wp), dimension(rank + 1)               :: nodes_rho
    real(wp), dimension(rank + 1)               :: nodes_gamma
    real(wp), dimension(rank + 1)               :: errors_at_nodes_rho
    real(wp), dimension(rank + 1)               :: errors_at_nodes_gamma
    real(wp), dimension(rank + 1)               :: b_coeffs_rho
    real(wp), dimension(rank + 1)               :: b_coeffs_gamma
    real(wp), dimension(rank + 1)               :: b_nodes_rho
    real(wp), dimension(rank + 1)               :: b_nodes_gamma
    real(wp), dimension(rank + 1)               :: b_errors_at_nodes_rho
    real(wp), dimension(rank + 1)               :: b_errors_at_nodes_gamma
    real(wp), dimension(rank)                   :: d_vec
    real(wp)                                    :: discrep, ker_discrep
    integer                                     :: num_iter, ker_num_iter
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Get the worst rho and gamma vectors for the kernel
    call worst_rho_and_gamma( &
        kernel_obj = kernel_to_approx, &
        num_terms = rank + 1, &
        tolerance = toler, &
        max_iter = max_iter, &
        rho_vec = rho_vec, &
        gamma_vec = gamma_vec, &
        coeffs_rho = coeffs_rho, &
        coeffs_gamma = coeffs_gamma, &
        a_matrix_rho = a_mat, &
        a_matrix_gamma = b_mat, &
        nodes_rho = nodes_rho, &
        nodes_gamma = nodes_gamma, &
        errors_at_nodes_rho = errors_at_nodes_rho, &
        errors_at_nodes_gamma = errors_at_nodes_gamma, &
        b_coeffs_rho = b_coeffs_rho, &
        b_coeffs_gamma = b_coeffs_gamma, &
        b_nodes_rho = b_nodes_rho, &
        b_nodes_gamma = b_nodes_gamma, &
        b_errors_at_nodes_rho = b_errors_at_nodes_rho, &
        b_errors_at_nodes_gamma = b_errors_at_nodes_gamma, &
        discrepancy = discrep, &
        num_iter = num_iter, &
        ker_discrep = ker_discrep, &
        ker_num_iter = ker_num_iter, &
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
    ! Set the intent(out) borsuk lower bound scalar:
    borsuk_lb = min( &
        minval( abs(b_errors_at_nodes_rho) ), &
        minval( abs( b_errors_at_nodes_gamma ) ) )
    ! Obtain the matrices used to form the numerically-optimal approximation
    call all_num_opt_mats_asymm( &
        kernel_obj = kernel_to_approx, &
        rho_vec = rho_vec, &
        gamma_vec = gamma_vec, &
        coeffs_rho = coeffs_rho, &
        coeffs_gamma = coeffs_gamma, &
        v_mat = v_mat, &
        w_mat = w_mat, &
        a_mat = a_mat, & 
        b_mat = b_mat, &
        d_vec = d_vec, &
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
  end subroutine num_opt_rankn_params

end module num_opt_rankn_params_mod