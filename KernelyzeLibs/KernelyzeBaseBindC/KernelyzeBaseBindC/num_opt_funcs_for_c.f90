! num_opt_funcs_for_c.f90
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
! created on: 2016-07-06
! updated on: 2016-07-06
! updated on: 2016-07-07 (continuing work toward a first version)
! updated on: 2016-07-14 (gaussian_sigma_fixed now has a setter for
!             its now-private sigma)
! updated on: 2016-08-11 (to accommodate changes to gaussian_pdf)
! updated on: 2017-04-14 (added num_opt_rankn_params_ptr, which uses
!             a C function pointer to a bivariate function [kernel]
!             to build a Fortran kernel which is then passed as
!             an argument to the Fortran procedure num_opt_rankn_params)
! updated on: 2017-04-20 (comment correction)
! updated on: 2017-04-24 (remove map-maker and other unused segments)
!
! A module that gives C types and a binding for functionality involving
! numerically-optimal approximations to (nondegenerate) totally positive
! kernels.

module num_opt_funcs_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use fort_c_interop_mod, only : copy_string
use kernel_proc_ptr_mod, only : kernel_proc_ptr
use kernel_nsect_mod, only : kernel_rankn
use num_opt_rankn_params_mod, only : num_opt_rankn_params
use num_opt_rankn_kernel_of_params_mod, only : &
    num_opt_rankn_kernel_of_params
use integral_nsect_disc_mod, only : integral_rankn_disc
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

subroutine num_opt_rankn_params_for_c( &
    c_ker_func, &
    c_x_lb, &
    c_x_ub, &
    c_y_lb, &
    c_y_ub, &
    c_approx_rank, &
    c_max_iter, &
    c_toler, &
    c_rho_vec, &
    c_gamma_vec, &
    c_v_mat, &
    c_w_mat, &
    c_borsuk_lb, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='num_opt_rankn_params')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  real(c_double), intent(in)        :: c_x_lb
  real(c_double), intent(in)        :: c_x_ub
  real(c_double), intent(in)        :: c_y_lb
  real(c_double), intent(in)        :: c_y_ub
  integer(c_int), intent(in)        :: c_approx_rank
  integer(c_int), intent(in)        :: c_max_iter
  real(c_double), intent(in)        :: c_toler
  real(c_double), intent(out)       :: c_rho_vec(c_approx_rank + 1)
  real(c_double), intent(out)       :: c_gamma_vec(c_approx_rank + 1)
  real(c_double), intent(out)       :: c_v_mat( &
      c_approx_rank + 1, c_approx_rank)
  real(c_double), intent(out)       :: c_w_mat( &
      c_approx_rank + 1, c_approx_rank)
  real(c_double), intent(out)       :: c_borsuk_lb
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
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
  call f_kernel%set_x_lb(c_x_lb)
  call f_kernel%set_x_ub(c_x_ub)
  call f_kernel%set_y_lb(c_y_lb)
  call f_kernel%set_y_ub(c_y_ub)
  ! Obtain the numerically-optimal parameters for this kernel
  call num_opt_rankn_params( &
      kernel_to_approx = f_kernel, &
      rank = c_approx_rank, &
      max_iter = c_max_iter, &
      toler = c_toler, &
      rho_vec = c_rho_vec, &
      gamma_vec = c_gamma_vec, &
      v_mat = c_v_mat, &
      w_mat = c_w_mat, &
      borsuk_lb = c_borsuk_lb, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
end subroutine num_opt_rankn_params_for_c  
    
subroutine num_opt_rankn_eval_for_c( &
    c_ker_func, &
    c_approx_rank, &
    c_rho_vec, &
    c_gamma_vec, &
    c_v_mat, &
    c_w_mat, &
    c_num_eval_pts, &
    c_eval_x_pts, &
    c_eval_y_pts, &
    c_results_of_eval, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='num_opt_rankn_eval')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_rho_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_gamma_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_v_mat( &
      c_approx_rank + 1, c_approx_rank)
  real(c_double), intent(in)        :: c_w_mat( &
      c_approx_rank + 1, c_approx_rank)
  integer(c_int), intent(in)        :: c_num_eval_pts
  real(c_double), intent(in)        :: c_eval_x_pts(c_num_eval_pts)
  real(c_double), intent(in)        :: c_eval_y_pts(c_num_eval_pts)
  real(c_double), intent(out)       :: c_results_of_eval(c_num_eval_pts)
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
  ! Local variables
  procedure(bivar_func), pointer    :: f_ker_func
  type(kernel_proc_ptr)             :: f_kernel
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  integer                           :: i
  class(kernel_rankn), allocatable  :: rankn_kernel
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(c_ker_func, f_ker_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  ! Build the rank-$n$ kernel to evaluate
  call num_opt_rankn_kernel_of_params( &
      kernel_to_approx = f_kernel, &
      rho_vec = c_rho_vec, &
      gamma_vec = c_gamma_vec, &
      v_mat = c_v_mat, &
      w_mat = c_w_mat, &
      rankn_kernel = rankn_kernel, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
  if (c_err_stat /= 0) then
    return
  end if
  ! Evaluate the rank-$n$ kernel at the given points
  do i = 1, c_num_eval_pts
    c_results_of_eval(i) = rankn_kernel%eval(c_eval_x_pts(i), c_eval_y_pts(i))
  end do
end subroutine num_opt_rankn_eval_for_c    

subroutine num_opt_rankn_integral_eval_for_c( &
    c_ker_func, &
    c_approx_rank, &
    c_rho_vec, &
    c_gamma_vec, &
    c_v_mat, &
    c_w_mat, &
    c_integral_is_over_x, &
    c_num_integral_pts, &
    c_integral_pts, &
    c_integral_wts, &
    c_num_eval_pts, &
    c_eval_pts, &
    c_results_of_eval, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='num_opt_rankn_integral_eval')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_rho_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_gamma_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_v_mat( &
      c_approx_rank + 1, c_approx_rank)
  real(c_double), intent(in)        :: c_w_mat( &
      c_approx_rank + 1, c_approx_rank)
  logical(c_bool), intent(in)       :: c_integral_is_over_x
  integer(c_int), intent(in)        :: c_num_integral_pts
  real(c_double), intent(in)        :: c_integral_pts(c_num_integral_pts)
  real(c_double), intent(in)        :: c_integral_wts(c_num_integral_pts)
  integer(c_int), intent(in)        :: c_num_eval_pts
  real(c_double), intent(in)        :: c_eval_pts(c_num_eval_pts)
  real(c_double), intent(out)       :: c_results_of_eval(c_num_eval_pts)
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
  ! Local variables
  procedure(bivar_func), pointer          :: f_ker_func
  type(kernel_proc_ptr)                   :: f_kernel
  character(len=err_msg_len)              :: err_msg
  integer(c_int)                          :: c_str_len
  integer                                 :: i
  class(kernel_rankn), allocatable        :: rankn_kernel
  class(integral_rankn_disc), allocatable :: rankn_integral
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(c_ker_func, f_ker_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  ! Build the rank-$n$ kernel to evaluate
  call num_opt_rankn_kernel_of_params( &
      kernel_to_approx = f_kernel, &
      rho_vec = c_rho_vec, &
      gamma_vec = c_gamma_vec, &
      v_mat = c_v_mat, &
      w_mat = c_w_mat, &
      rankn_kernel = rankn_kernel, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
  if (c_err_stat /= 0) then
    return
  end if
  ! Compute the integral (a discrete sum, in this case)
  call rankn_kernel%kernel_integral( &
      sum_is_over_x = logical(c_integral_is_over_x), &
      weights = c_integral_wts, &
      eval_pts = c_integral_pts, &
      integral = rankn_integral, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
  if (c_err_stat /= 0) then
    return
  end if
  ! Evaluate the rank-$n$ kernel integral at the given points
  do i = 1, c_num_eval_pts
    c_results_of_eval(i) = rankn_integral%eval(c_eval_pts(i))
  end do
    end subroutine num_opt_rankn_integral_eval_for_c    
    
subroutine num_opt_rankn_integral_coeffs_for_c( &
    c_ker_func, &
    c_approx_rank, &
    c_rho_vec, &
    c_gamma_vec, &
    c_v_mat, &
    c_w_mat, &
    c_integral_is_over_x, &
    c_num_integral_pts, &
    c_integral_pts, &
    c_integral_wts, &
    c_num_coeffs, &
    c_coeffs, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='num_opt_rankn_integral_coeffs')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_rho_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_gamma_vec(c_approx_rank + 1)
  real(c_double), intent(in)        :: c_v_mat( &
      c_approx_rank + 1, c_approx_rank)
  real(c_double), intent(in)        :: c_w_mat( &
      c_approx_rank + 1, c_approx_rank)
  logical(c_bool), intent(in)       :: c_integral_is_over_x
  integer(c_int), intent(in)        :: c_num_integral_pts
  real(c_double), intent(in)        :: c_integral_pts(c_num_integral_pts)
  real(c_double), intent(in)        :: c_integral_wts(c_num_integral_pts)
  integer(c_int), intent(in)        :: c_num_coeffs
  real(c_double), intent(out)       :: c_coeffs(c_num_coeffs)
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
  ! Local variables
  procedure(bivar_func), pointer          :: f_ker_func
  type(kernel_proc_ptr)                   :: f_kernel
  character(len=err_msg_len)              :: err_msg
  integer(c_int)                          :: c_str_len
  integer                                 :: i
  class(kernel_rankn), allocatable        :: rankn_kernel
  class(integral_rankn_disc), allocatable :: rankn_integral
  real(wp), allocatable                   :: coeff_vec(:)
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(c_ker_func, f_ker_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  ! Build the rank-$n$ kernel to evaluate
  call num_opt_rankn_kernel_of_params( &
      kernel_to_approx = f_kernel, &
      rho_vec = c_rho_vec, &
      gamma_vec = c_gamma_vec, &
      v_mat = c_v_mat, &
      w_mat = c_w_mat, &
      rankn_kernel = rankn_kernel, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
  if (c_err_stat /= 0) then
    return
  end if
  ! Compute the integral (a discrete sum, in this case)
  call rankn_kernel%kernel_integral( &
      sum_is_over_x = logical(c_integral_is_over_x), &
      weights = c_integral_wts, &
      eval_pts = c_integral_pts, &
      integral = rankn_integral, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
  if (c_err_stat /= 0) then
    return
  end if
  ! Obtain the coefficient vector
  coeff_vec = rankn_integral%get_coeff_vec()
  ! Check to be sure that the claimed size is accurate
  if (size(coeff_vec) /= c_num_coeffs) then
    c_err_stat = 6
    err_msg = "num_opt_rankn_integral_coeffs: incorrect number of coeffs given"
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
  ! Copy the coefficients into the output array
  do i = 1, c_num_coeffs
    c_coeffs(i) = coeff_vec(i)
  end do
end subroutine num_opt_rankn_integral_coeffs_for_c

end module num_opt_funcs_for_c_mod
