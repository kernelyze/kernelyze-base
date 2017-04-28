! borsuk_lower_bound_for_c.f90
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
! Fortran functionality to find Borsuk lower
! bounds on the uniform error of approximation
! of a given kernel with a rank-$n$ kernel.

module borsuk_lower_bound_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use borsuk_lower_bound_mod, only : borsuk_lower_bound
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
 
! The subroutine borsuk_lower_bound_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine borsuk_lower_bound_for_c( &
    c_ker_func, &
    size_of_rho, &
    rho_vec, &
    is_over_x, &
    tolerance, &
    max_iter, &
    n_grid_pts, &
    grid, &
    coeffs, &
    nodes, &
    errors_at_nodes, &
    discrepancy, &
    num_iter, &
    c_err_stat, &
    c_err_msg ) bind(C, name='borsuk_lower_bound')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  integer(c_int), intent(in)        :: size_of_rho
  real(c_double), intent(in)        :: rho_vec(size_of_rho)
  logical(c_bool), intent(in)       :: is_over_x
  real(c_double), intent(in)        :: tolerance
  integer(c_int), intent(in)        :: max_iter
  integer(c_int), intent(in)        :: n_grid_pts
  real(c_double), intent(in)        :: grid(n_grid_pts)
  real(c_double), intent(out)       :: coeffs(size_of_rho)
  real(c_double), intent(out)       :: nodes(size_of_rho)
  real(c_double), intent(out)       :: errors_at_nodes(size_of_rho)
  real(c_double), intent(out)       :: discrepancy
  integer(c_int), intent(out)       :: num_iter
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
  ! Call the Fortran borsuk_lower_bound subroutine
  call borsuk_lower_bound( &
      f_kernel, &
      rho_vec, &
      logical(is_over_x), &
      tolerance, &
      max_iter, &
      grid, &
      coeffs, &
      nodes, &
      errors_at_nodes, &
      discrepancy, &
      num_iter, &
      c_err_stat, &
      err_msg )
  ! Handle any errors
  if (c_err_stat /= 0) then
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
end subroutine borsuk_lower_bound_for_c
  
end module borsuk_lower_bound_for_c_mod  