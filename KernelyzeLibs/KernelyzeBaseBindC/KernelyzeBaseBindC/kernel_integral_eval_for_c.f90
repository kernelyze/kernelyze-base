! kernel_integral_eval_for_c.f90
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
! created on: 2017-04-23
!
! A module that gives C types and a binding for functionality involving
! evaluation of discrete integrals (weighted sums) of kernels.

module kernel_integral_eval_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use fort_c_interop_mod, only : copy_string
use kernel_proc_ptr_mod, only : kernel_proc_ptr
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

! The subroutine kernel_integral_eval_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.    
subroutine kernel_integral_eval_for_c( &
    c_ker_func, &
    c_integral_is_over_x, &
    c_num_integral_pts, &
    c_integral_pts, &
    c_integral_wts, &
    c_num_eval_pts, &
    c_eval_pts, &
    c_results_of_eval) &
    bind(C, name='kernel_integral_eval')
  ! Arguments
  type(c_funptr), intent(in), value :: c_ker_func
  logical(c_bool), intent(in)       :: c_integral_is_over_x
  integer(c_int), intent(in)        :: c_num_integral_pts
  real(c_double), intent(in)        :: c_integral_pts(c_num_integral_pts)
  real(c_double), intent(in)        :: c_integral_wts(c_num_integral_pts)
  integer(c_int), intent(in)        :: c_num_eval_pts
  real(c_double), intent(in)        :: c_eval_pts(c_num_eval_pts)
  real(c_double), intent(out)       :: c_results_of_eval(c_num_eval_pts)
  ! Local variables
  procedure(bivar_func), pointer    :: f_ker_func
  type(kernel_proc_ptr)             :: f_kernel
  integer                           :: i
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(c_ker_func, f_ker_func)
  ! Now build a singlevar_proc_ptr object
  call f_kernel%set_proc_ptr(f_ker_func)
  ! Evaluate the rank-$n$ kernel integral at the given points
  if (c_integral_is_over_x) then
    do i = 1, c_num_eval_pts
      c_results_of_eval(i) = &
          dot_product( c_integral_wts, &
                       f_kernel%eval_elt( c_integral_pts , c_eval_pts(i) ) )
    end do
  else ! In this case, the integral is over y
    do i = 1, c_num_eval_pts
      c_results_of_eval(i) = &
          dot_product( c_integral_wts, &
                       f_kernel%eval_elt( c_eval_pts(i) , c_integral_pts ) )
    end do
  end if
end subroutine kernel_integral_eval_for_c
    
end module kernel_integral_eval_for_c_mod    