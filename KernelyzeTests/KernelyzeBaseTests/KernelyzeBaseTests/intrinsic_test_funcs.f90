! intrinsic_test_funcs.f90
!
! Copyright (c) 2015 by Kernelyze LLC
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
! created on: 2015-04-25
! updated on: 2015-04-25
!
! A module of example intrinsic functions to use in testing functionality.
  
module intrinsic_test_funcs_mod

  use set_precision, only : wp
  use singlevar_proc_ptr_mod, only : singlevar_proc_ptr

  implicit none
  
  abstract interface
    pure function univar_func(x) result(y)
      ! Bring in the working precision
      import wp
      ! Argument
      real(wp), intent(in)  :: x
      ! Result
      real(wp)              :: y
    end function univar_func
  end interface

  type(singlevar_proc_ptr), protected :: sin_func
  type(singlevar_proc_ptr), protected :: exp_func
  type(singlevar_proc_ptr), protected :: log_func

  contains
  
  subroutine init_test_funcs()
    procedure(univar_func), pointer :: sin_ptr, exp_ptr, log_ptr
    sin_ptr => dsin
    exp_ptr => dexp
    log_ptr => dlog
    call sin_func%set_proc_ptr(sin_ptr)
    call exp_func%set_proc_ptr(exp_ptr)
    call log_func%set_proc_ptr(log_ptr)
  end subroutine init_test_funcs

end module intrinsic_test_funcs_mod