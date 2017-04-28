! taylor_funcs_for_c.f90
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
!
! created on: 2017-04-20
! updated on: 2017-04-22 (added ability to pass ptr to deriv func)
! updated on: 2017-04-23 (corrected comment)
! updated on: 2017-04-27 (use increased finite-diff step size flexibility)
!
! A module that gives C types and a binding for functionality involving
! truncated-Taylor-series approximations to kernels.

module taylor_funcs_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use fort_c_interop_mod, only : copy_string
use kernel_proc_ptr_mod, only : kernel_proc_ptr
use kernel_nmonom_mod, only : kernel_rankn
use taylor_series_rankn_params_mod, only : taylor_series_rankn_params
use taylor_series_rankn_kernel_of_params_mod, only : &
    taylor_series_rankn_kernel_of_params
use integral_nmonom_disc_mod, only : integral_rankn_disc
use, intrinsic :: iso_c_binding

implicit none

! Here are the abstract interfaces for a pure
! function that takes two C doubles and returns
! a C double (evaluation) and for a pure 
! function that takes two C integers and two
! C doubles and returns the partial derivative
! of the order specified by the two integers
! evaluated at the two doubles.
abstract interface
    pure function bivar_func(x, y) result(res)
      import c_double
      ! Argument
      real(c_double), intent(in)  :: x
      real(c_double), intent(in)  :: y
      ! Result
      real(c_double)              :: res
    end function bivar_func
    ! mth deriv w. r. t. x of nth deriv w. r. t. y,
    ! so m == 0 means just a y partial and n == 0
    ! means just an x partial, and m == n == 0 would
    ! be just the value.
    pure function bivar_mnder(m, n , x, y) result(deriv)
      ! Bring in the working precision
      import wp
      ! Arguments
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      real(wp), intent(in)  :: x
      real(wp), intent(in)  :: y
      ! Result
      real(wp)              :: deriv
    end function bivar_mnder
end interface

contains

subroutine taylor_rankn_params_for_c( &
    c_ker_func, &
    c_approx_rank, &
    c_x_center, &
    c_y_center, &
    c_v_matrix, &
    c_w_matrix, &
    c_err_stat, &
    c_err_msg, &
    c_fin_diff_delta) &
    bind(C, name='taylor_rankn_params')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_x_center
  real(c_double), intent(in)        :: c_y_center
  real(c_double), intent(in)        :: c_fin_diff_delta
  real(c_double), intent(out)       :: c_v_matrix(c_approx_rank, c_approx_rank)
  real(c_double), intent(out)       :: c_w_matrix(c_approx_rank, c_approx_rank)
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
  ! Obtain the numerically-optimal parameters for this kernel
  call taylor_series_rankn_params( &
      kernel_to_approx = f_kernel, &
      x_center = c_x_center, &
      y_center = c_y_center, &
      rank = c_approx_rank, &
      v_matrix = c_v_matrix, &
      w_matrix = c_w_matrix, &
      err_stat = c_err_stat, &
      err_msg = err_msg, &
      fin_diff_delta = c_fin_diff_delta)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
end subroutine taylor_rankn_params_for_c
    
subroutine taylor_rankn_params_fdf_for_c( &
    c_ker_func, &
    c_deriv_func, &
    c_approx_rank, &
    c_x_center, &
    c_y_center, &
    c_v_matrix, &
    c_w_matrix, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='taylor_rankn_params_fdf')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  type(c_funptr), intent(in), value :: c_deriv_func
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_x_center
  real(c_double), intent(in)        :: c_y_center
  real(c_double), intent(out)       :: c_v_matrix(c_approx_rank, c_approx_rank)
  real(c_double), intent(out)       :: c_w_matrix(c_approx_rank, c_approx_rank)
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
  ! Local variables
  procedure(bivar_func), pointer    :: f_ker_func
  procedure(bivar_mnder), pointer   :: f_deriv_func
  type(kernel_proc_ptr)             :: f_kernel
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  ! Body
  ! First convert the C procedure pointers to Fortran procedure pointers
  call c_f_procpointer(c_ker_func, f_ker_func)
  call c_f_procpointer(c_deriv_func, f_deriv_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  call f_kernel%set_deriv_ptr(f_deriv_func)
  ! Obtain the numerically-optimal parameters for this kernel
  call taylor_series_rankn_params( &
      kernel_to_approx = f_kernel, &
      x_center = c_x_center, &
      y_center = c_y_center, &
      rank = c_approx_rank, &
      v_matrix = c_v_matrix, &
      w_matrix = c_w_matrix, &
      err_stat = c_err_stat, &
      err_msg = err_msg)
  ! Copy the Fortran error message string to the output C string
  c_str_len = err_msg_len
  call copy_string(c_err_msg, err_msg, c_str_len, 1)
end subroutine taylor_rankn_params_fdf_for_c    
    
subroutine taylor_rankn_eval_for_c( &
    c_approx_rank, &
    c_x_center, &
    c_y_center, &
    c_v_matrix, &
    c_w_matrix, &
    c_num_eval_pts, &
    c_eval_x_pts, &
    c_eval_y_pts, &
    c_results_of_eval, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='taylor_rankn_eval')
  ! Arguments
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_x_center
  real(c_double), intent(in)        :: c_y_center
  real(c_double), intent(in)        :: c_v_matrix(c_approx_rank, c_approx_rank)
  real(c_double), intent(in)        :: c_w_matrix(c_approx_rank, c_approx_rank)
  integer(c_int), intent(in)        :: c_num_eval_pts
  real(c_double), intent(in)        :: c_eval_x_pts(c_num_eval_pts)
  real(c_double), intent(in)        :: c_eval_y_pts(c_num_eval_pts)
  real(c_double), intent(out)       :: c_results_of_eval(c_num_eval_pts)
  ! Arguments for error-handling
  integer(c_int), intent(out)         :: c_err_stat
  character(kind=c_char), intent(out) :: c_err_msg(*)
  ! Local variables
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  integer                           :: i
  class(kernel_rankn), allocatable  :: rankn_kernel
  ! Body
  ! Build the rank-$n$ kernel to evaluate
  call taylor_series_rankn_kernel_of_params( &
      x_center = c_x_center, &
      y_center = c_y_center, &
      v_matrix = c_v_matrix, &
      w_matrix = c_w_matrix, &
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
end subroutine taylor_rankn_eval_for_c
    
subroutine taylor_rankn_integral_eval_for_c( &
    c_approx_rank, &
    c_x_center, &
    c_y_center, &
    c_v_matrix, &
    c_w_matrix, &
    c_integral_is_over_x, &
    c_num_integral_pts, &
    c_integral_pts, &
    c_integral_wts, &
    c_num_eval_pts, &
    c_eval_pts, &
    c_results_of_eval, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='taylor_rankn_integral_eval')
  ! Arguments
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_x_center
  real(c_double), intent(in)        :: c_y_center
  real(c_double), intent(in)        :: c_v_matrix(c_approx_rank, c_approx_rank)
  real(c_double), intent(in)        :: c_w_matrix(c_approx_rank, c_approx_rank)
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
  character(len=err_msg_len)                :: err_msg
  integer(c_int)                            :: c_str_len
  integer                                   :: i
  class(kernel_rankn), allocatable          :: rankn_kernel
  class(integral_rankn_disc), allocatable   :: rankn_integral
  ! Body
  ! Build the rank-$n$ kernel to evaluate
  call taylor_series_rankn_kernel_of_params( &
      x_center = c_x_center, &
      y_center = c_y_center, &
      v_matrix = c_v_matrix, &
      w_matrix = c_w_matrix, &
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
end subroutine taylor_rankn_integral_eval_for_c
    
subroutine taylor_rankn_integral_coeffs_for_c( &
    c_approx_rank, &
    c_x_center, &
    c_y_center, &
    c_v_matrix, &
    c_w_matrix, &
    c_integral_is_over_x, &
    c_num_integral_pts, &
    c_integral_pts, &
    c_integral_wts, &
    c_num_coeffs, &
    c_coeffs, &
    c_err_stat, &
    c_err_msg) &
    bind(C, name='taylor_rankn_integral_coeffs')
  ! Arguments
  integer(c_int), intent(in)        :: c_approx_rank
  real(c_double), intent(in)        :: c_x_center
  real(c_double), intent(in)        :: c_y_center
  real(c_double), intent(in)        :: c_v_matrix(c_approx_rank, c_approx_rank)
  real(c_double), intent(in)        :: c_w_matrix(c_approx_rank, c_approx_rank)
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
  character(len=err_msg_len)                :: err_msg
  integer(c_int)                            :: c_str_len
  integer                                   :: i
  class(kernel_rankn), allocatable          :: rankn_kernel
  class(integral_rankn_disc), allocatable   :: rankn_integral
  real(wp), allocatable                     :: coeff_vec(:)
  ! Body
  ! Build the rank-$n$ kernel to evaluate
  call taylor_series_rankn_kernel_of_params( &
      x_center = c_x_center, &
      y_center = c_y_center, &
      v_matrix = c_v_matrix, &
      w_matrix = c_w_matrix, &
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
    err_msg = "taylor_rankn_integral_coeffs: incorrect number of coeffs given"
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
  ! Copy the coefficients into the output array
  do i = 1, c_num_coeffs
    c_coeffs(i) = coeff_vec(i)
  end do
end subroutine taylor_rankn_integral_coeffs_for_c
    
end module taylor_funcs_for_c_mod
