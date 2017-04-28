! worst_rho_and_gamma_for_c.f90
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
! created on: 2017-04-14
!
! A module that gives C types and a binding for
! Fortran functionality to find the worst rho 
! and gamma arrays for approximation of a given
! kernel with a rank-$n$ kernel.

module worst_rho_and_gamma_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use worst_rho_and_gamma_mod, only : worst_rho_and_gamma
use kernel_proc_ptr_mod, only : kernel_proc_ptr
use fort_c_interop_mod, only : copy_string
use, intrinsic :: iso_c_binding

implicit none

! Here is the abstract interface for a pure
! function that takes two C doubles and returns
! a C double.
abstract interface
    pure function bivar_func(x, y) result(res)
      import c_double
      ! Argument
      real(c_double), intent(in)  :: x
      real(c_double), intent(in)  :: y
      ! Result
      real(c_double)              :: res
    end function bivar_func
  end interface
contains
 
! The subroutine worst_rho_and_gamma_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine worst_rho_and_gamma_for_c( &
    c_ker_func, &
    x_lb, &
    x_ub, &
    y_lb, &
    y_ub, &
    num_terms, &
    tolerance, &
    max_iter, &
    rho_vec, &
    gamma_vec, &
    coeffs_rho, &
    coeffs_gamma, &
    a_matrix_rho, &
    a_matrix_gamma, &
    nodes_rho, &
    nodes_gamma, &
    errors_at_nodes_rho, &
    errors_at_nodes_gamma, &
    b_coeffs_rho, &
    b_coeffs_gamma, &
    b_nodes_rho, &
    b_nodes_gamma, &
    b_errors_at_nodes_rho, &
    b_errors_at_nodes_gamma, &
    discrepancy, &
    num_iter, &
    ker_discrep, &
    ker_num_iter, &
    c_err_stat, &
    c_err_msg ) bind(C, name='worst_rho_and_gamma')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  real(c_double), intent(in)        :: x_lb
  real(c_double), intent(in)        :: x_ub
  real(c_double), intent(in)        :: y_lb
  real(c_double), intent(in)        :: y_ub
  integer(c_int), intent(in)        :: num_terms
  real(c_double), intent(in)        :: tolerance
  integer(c_int), intent(in)        :: max_iter
  real(c_double), intent(out)       :: rho_vec(num_terms)
  real(c_double), intent(out)       :: gamma_vec(num_terms)
  real(c_double), intent(out)       :: coeffs_rho(num_terms)
  real(c_double), intent(out)       :: coeffs_gamma(num_terms)
  real(c_double), intent(out)       :: a_matrix_rho(num_terms, num_terms)
  real(c_double), intent(out)       :: a_matrix_gamma(num_terms, num_terms)
  real(c_double), intent(out)       :: nodes_rho(num_terms)
  real(c_double), intent(out)       :: nodes_gamma(num_terms)
  real(c_double), intent(out)       :: errors_at_nodes_rho(num_terms)
  real(c_double), intent(out)       :: errors_at_nodes_gamma(num_terms)
  real(c_double), intent(out)       :: b_coeffs_rho(num_terms)
  real(c_double), intent(out)       :: b_coeffs_gamma(num_terms)
  real(c_double), intent(out)       :: b_nodes_rho(num_terms)
  real(c_double), intent(out)       :: b_nodes_gamma(num_terms)
  real(c_double), intent(out)       :: b_errors_at_nodes_rho(num_terms)
  real(c_double), intent(out)       :: b_errors_at_nodes_gamma(num_terms)
  real(c_double), intent(out)       :: discrepancy
  integer(c_int), intent(out)       :: num_iter
  real(c_double), intent(out)       :: ker_discrep
  integer(c_int), intent(out)       :: ker_num_iter
  ! The 'intent out' error variables
  integer(c_int), intent(out)       :: c_err_stat
  ! How to handle the error message string?
  character(c_char), intent(out)    :: c_err_msg(*)
  ! Local variables
  procedure(bivar_func), pointer    :: f_ker_func
  type(kernel_proc_ptr)             :: f_kernel
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(c_ker_func, f_ker_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  call f_kernel%set_x_lb(x_lb)
  call f_kernel%set_x_ub(x_ub)
  call f_kernel%set_y_lb(y_lb)
  call f_kernel%set_y_ub(y_ub)
  ! Call the Fortran worst_rho_and_gamma subroutine
  call worst_rho_and_gamma( &
      f_kernel, &
      num_terms, & 
      tolerance, &
      max_iter, &
      rho_vec, &
      gamma_vec, &
      coeffs_rho, &
      coeffs_gamma, &
      a_matrix_rho, &
      a_matrix_gamma, &
      nodes_rho, &
      nodes_gamma, &
      errors_at_nodes_rho, &
      errors_at_nodes_gamma, &
      b_coeffs_rho, &
      b_coeffs_gamma, &
      b_nodes_rho, &
      b_nodes_gamma, &
      b_errors_at_nodes_rho, &
      b_errors_at_nodes_gamma, &
      discrepancy, &
      num_iter, &
      ker_discrep, &
      ker_num_iter, &
      c_err_stat, &
      err_msg)
  ! Handle any errors
  if (c_err_stat /= 0) then
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
end subroutine worst_rho_and_gamma_for_c
  
end module worst_rho_and_gamma_for_c_mod  