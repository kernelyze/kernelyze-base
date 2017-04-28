! brent_for_c.f90
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
! A module that gives C types and a binding for Brent
! root-finding and one-dimensional optimization.

module brent_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use brent_mod, only : zeroin, fmin, brentZero, brentmin
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
  
! The function zeroin_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
function zeroin_for_c( &
    a, &
    b, &
    f, &
    tol ) result(z) bind(C, name='zeroin')
  ! Arguments
  real(c_double), intent(in)        :: a
  real(c_double), intent(in)        :: b
  type(c_funptr), intent(in), value :: f
  real(c_double), intent(in)        :: tol
  ! Result
  real(c_double)                    :: z
  ! Local variables
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Finally, call the Fortran zeroin function
  z = zeroin(a, b, svar, tol)
end function zeroin_for_c
   
! The subroutine brent_zero_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine brent_zero_for_c( &
    a, &
    b, &
    f, &
    tol, &
    max_iter, &
    neval, &
    err_code, &
    x_zero, &
    f_zero ) bind(C, name='brent_zero')
  ! Arguments
  real(c_double), intent(in)        :: a
  real(c_double), intent(in)        :: b
  type(c_funptr), intent(in), value :: f
  real(c_double), intent(in)        :: tol
  integer(c_int), intent(in)        :: max_iter
  integer(c_int), intent(out)       :: neval
  integer(c_int), intent(out)       :: err_code
  real(c_double), intent(out)       :: x_zero
  real(c_double), intent(out)       :: f_zero
  ! Local variables
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Call the Fortran brentZero subroutine
  call brentZero(a, b, svar, tol, max_iter, neval, err_code, x_zero, f_zero)
end subroutine brent_zero_for_c
    
! The subroutine fmin_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine fmin_for_c( &
    a, &
    b, &
    f, &
    tol, &
    xopt, &
    fopt, &
    c_err_stat, &
    c_err_msg ) bind(C, name='f_min')
  ! Arguments
  real(c_double), intent(in)        :: a
  real(c_double), intent(in)        :: b
  type(c_funptr), intent(in), value :: f
  real(c_double), intent(in)        :: tol
  real(c_double), intent(out)       :: xopt
  real(c_double), intent(out)       :: fopt
  ! The 'intent out' error variables
  integer(c_int), intent(out)       :: c_err_stat
  ! How to handle the error message string?
  character(c_char), intent(out)    :: c_err_msg(*)
  ! Local variables
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  character(len=err_msg_len)        :: err_msg
  integer(c_int)                    :: c_str_len
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Call the Fortran fmin subroutine
  call fmin(a, b, svar, tol, xopt, fopt, c_err_stat, err_msg)
  ! Handle any errors
  if (c_err_stat /= 0) then
    ! Copy the Fortran error message string to the output C string
    c_str_len = err_msg_len
    call copy_string(c_err_msg, err_msg, c_str_len, 1)
    return
  end if
end subroutine fmin_for_c
    
! The subroutine brent_min_for_c cannot be declared as pure,
! even though it has no side effects, because the intrinsic
! conversion procedure c_f_procpointer is not declared pure.
subroutine brent_min_for_c( &
    a, &
    b, &
    f, &
    tol, &
    max_iter, &
    neval, &
    err_code, &
    x_min, &
    f_min ) bind(C, name='brent_min')
  ! Arguments
  real(c_double), intent(in)        :: a
  real(c_double), intent(in)        :: b
  type(c_funptr), intent(in), value :: f
  real(c_double), intent(in)        :: tol
  integer(c_int), intent(in)        :: max_iter
  integer(c_int), intent(out)       :: neval
  integer(c_int), intent(out)       :: err_code
  real(c_double), intent(out)       :: x_min
  real(c_double), intent(out)       :: f_min
  ! Local variables
  procedure(func), pointer          :: f_fortptr
  type(singlevar_proc_ptr)          :: svar
  ! Body
  ! First convert the C procedure pointer to a Fortran procedure pointer
  call c_f_procpointer(f, f_fortptr)
  ! Now build a singlevar_proc_ptr object
  call svar%set_proc_ptr(f_fortptr)
  ! Call the Fortran brentZero subroutine
  call brentmin(a, b, svar, tol, max_iter, neval, err_code, x_min, f_min)
end subroutine brent_min_for_c
    
end module brent_for_c_mod