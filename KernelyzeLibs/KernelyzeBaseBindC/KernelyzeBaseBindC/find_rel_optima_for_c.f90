! find_rel_optima_for_c.f90
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
! Fortran find-relative-optima functionality.

module find_rel_optima_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use find_rel_optima_mod, only : find_rel_optima
use singlevar_proc_ptr_mod, only : singlevar_proc_ptr
use fort_c_interop_mod, only : copy_string
use, intrinsic :: iso_c_binding

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
 
! The subroutine find_rel_optima_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine find_rel_optima_for_c( &
    result_optima, &
    f, &
    n_grid_pts, &
    grid, &
    toler, &
    c_err_stat, &
    c_err_msg ) bind(C, name='find_rel_optima')
  ! Arguments
  ! result_optima holds both location and value for each relative
  ! optimum, so it has 2 elements for each interval implied by
  ! the grid points.
  integer(c_int), intent(in)        :: n_grid_pts
  real(c_double), intent(out)       :: result_optima(n_grid_pts - 1, 2)
  type(c_funptr), intent(in), value :: f
  real(c_double), intent(in)        :: grid(n_grid_pts)
  real(c_double), intent(in)        :: toler
  ! The 'intent out' error variables
  integer(c_int), intent(out)       :: c_err_stat
  ! How to handle the error message string?
  character(c_char), intent(out)    :: c_err_msg(*)
  ! Local variables
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  real(wp), allocatable             :: res_optima(: , :)
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Call the Fortran find_rel_optima subroutine
  call find_rel_optima( &
      res_optima, &
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
  ! Check the size of the res_optima array; if it is greater
  ! than max_num_optima then res_optima is larger than
  ! result_optima and I must throw an error
  if ( size(res_optima, 1) /= (n_grid_pts - 1) ) then
    c_err_stat = -1
    err_msg = 'number of optima must equal number of grid points minus 1.'
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
  ! Copy res_optima into the result_optima array
  result_optima = res_optima
end subroutine find_rel_optima_for_c
  
end module find_rel_optima_for_c_mod  