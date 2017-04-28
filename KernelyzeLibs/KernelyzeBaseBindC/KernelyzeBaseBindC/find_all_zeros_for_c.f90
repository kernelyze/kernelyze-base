! find_all_zeros_for_c.f90
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
! Fortran find-all-zeros functionality.

module find_all_zeros_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use find_all_zeros_mod, only : find_all_zeros
use singlevar_proc_ptr_mod, only : singlevar_proc_ptr
use fort_c_interop_mod, only : copy_string
use, intrinsic :: iso_c_binding
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none

! Here is the abstract interface for a pure
! function that takes a C double and returns
! a C double.
abstract interface
  pure function func(x)
    import :: c_double
    real(c_double), intent(in)  :: x
    real(c_double)              :: func
  end function func
end interface

  contains
 
! The subroutine find_all_zeros_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
! Because result_zeros is of unknown size, I pass its allocated
! size as the first argument (max_num_zeros) and then the
! number of elements actually used to store zeros in the
! second (output) argument (num_zeros).  The remaining elements
! are filled with quiet NaNs as a precaution.
subroutine find_all_zeros_for_c( &
    max_num_zeros, &
    num_zeros, &
    result_zeros, &
    f, &
    n_grid_pts, &
    grid, &
    toler, &
    c_err_stat, &
    c_err_msg ) bind(C, name='find_all_zeros')
  ! Arguments
  integer(c_int), intent(in)        :: max_num_zeros
  integer(c_int), intent(out)       :: num_zeros
  real(c_double), intent(out)       :: result_zeros(max_num_zeros)
  type(c_funptr), intent(in), value :: f
  integer(c_int), intent(in)        :: n_grid_pts
  real(c_double), intent(in)        :: grid(n_grid_pts)
  real(c_double), intent(in)        :: toler
  ! The 'intent out' error variables
  integer(c_int), intent(out)       :: c_err_stat
  ! How to handle the error message string?
  character(c_char), intent(out)    :: c_err_msg(*)
  ! Local variables
  integer                           :: i
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  real(wp), allocatable             :: res_zeros(:)
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Call the Fortran find_all_zeros subroutine
  call find_all_zeros( &
      res_zeros, &
      svar, &
      grid, &
      toler, &
      c_err_stat, &
      err_msg )
  ! Handle any errors
  if (c_err_stat /= 0) then
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
  ! Check the size of the res_zeros array; if it is greater
  ! than max_num_zeros then res_zeros is larger than
  ! result_zeros and I must throw an error
  num_zeros = size(res_zeros)
  if (num_zeros > max_num_zeros) then
    c_err_stat = -1
    err_msg = 'number of zeros found greater than max number given.'
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
  ! Copy res_zeros into the result_zeros array
  do i = 1, num_zeros
    result_zeros(i) = res_zeros(i)
  end do
  ! As a safety measure, fill the remaining elements of
  ! result_zeros (that are not occupied by actual zeros
  ! that were found) with quiet NaNs.  Of course, this
  ! loop may have zero iterations if num_zeros == max_num_zeros.
  do i = num_zeros + 1, size(result_zeros)
    result_zeros(i) = ieee_value(1E0_wp, ieee_quiet_nan)
  end do
end subroutine find_all_zeros_for_c
  
end module find_all_zeros_for_c_mod  