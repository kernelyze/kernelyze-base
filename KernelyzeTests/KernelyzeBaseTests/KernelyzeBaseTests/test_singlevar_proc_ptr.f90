! test_singlevar_proc_ptr.f90
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
! created on: 2016-08-29
! updated on: 2016-08-30 (added tests for find_all_zeros
!             and find_rel_optima)
! updated on: 2016-08-31 (fixes to test for find_rel_optima)
!
! A module of unit tests for singlevar_proc_ptrs.
  
module test_singlevar_proc_ptr_mod

use set_precision, only : wp
use constants_mod, only : pi, err_msg_len
use chebyshev_points_mod, only : chebyshev_points
use singlevar_proc_ptr_mod, only : singlevar_proc_ptr
use find_all_zeros_mod, only : find_all_zeros
use find_rel_optima_mod, only : find_rel_optima

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
  pure function univar_nder(n , x) result(deriv)
    ! Bring in the working precision
    import wp
    ! Arguments
    integer, intent(in)   :: n
    real(wp), intent(in)  :: x
    ! Result
    real(wp)              :: deriv
  end function univar_nder
end interface

contains

function test_singlevar_proc_ptr(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality to use ' &
      // 'singlevar_funcs based on procedure pointers: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. &
      test_singlevar_proc_ptr_sin(unit_to_use)
  test_pass = test_pass .and. &
      test_svar_pptr_find_all_zeros(unit_to_use)
  test_pass = test_pass .and. &
      test_svar_pptr_find_rel_optima(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_singlevar_proc_ptr

pure function sinfunc(x) result(y)
  ! Arguments
  real(wp), intent(in)  :: x
  ! Result
  real(wp)              :: y
  ! Body
  y = sin(x)
end function sinfunc

function test_singlevar_proc_ptr_sin(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  integer                           :: i
  real(wp)                          :: ptr_eval, sin_eval
  type(singlevar_proc_ptr)          :: svar_pptr
  procedure(univar_func), pointer   :: sin_ptr
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test use of sin in a singlevar_proc_ptr: '
  sin_ptr => sinfunc
  call svar_pptr%set_proc_ptr(sin_ptr)
  do i = 1, 10
    ptr_eval = svar_pptr%eval( real(i, wp) / 1E1_wp )
    sin_eval = sin( real(i, wp) / 1E1_wp )
    write(unit_to_use , * ) 'Evals from singlevar_proc_ptr and sin: ', &
        ptr_eval, ' ', sin_eval
    if (abs(ptr_eval - sin_eval) > epsilon(1E0_wp)) then
      test_pass = .false.
      write(unit_to_use , *) 'Test failed, singlevar_proc_ptr and ', &
          'sin should give same result but did not.'
    end if
  end do
  write(unit_to_use , *) 'Did the singlevar_proc_ptr test for the ' &
      // 'sin function pass? ', test_pass
end function test_singlevar_proc_ptr_sin

function test_svar_pptr_find_all_zeros(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  integer                           :: i
  real(wp), allocatable             :: zeros_found(:)
  real(wp)                          :: grid(200)
  type(singlevar_proc_ptr)          :: svar_pptr
  procedure(univar_func), pointer   :: sin_ptr
  integer                           :: local_err_stat
  character(len=err_msg_len)        :: local_err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test use of a singlevar_proc_ptr ', &
      'in find_all_zeros.'
  sin_ptr => sinfunc
  call svar_pptr%set_proc_ptr(sin_ptr)
  grid = 1E1_wp * chebyshev_points(200) ! From -10 to 10
  call find_all_zeros( &
      zeros_found, &
      svar_pptr, &
      grid, &
      1E-15_wp, & ! toler
      local_err_stat, &
      local_err_msg )
  if (local_err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use , *) 'Test failed with message ' // local_err_msg
  end if
  if (size(zeros_found) /= 7) then
    test_pass = .false.
    write(unit_to_use , *) 'Test failed, should have found 7 zeros ', &
        'but found ', size(zeros_found), ' instead'
  else
    do i = -3, 3
      if (abs(zeros_found(i+4) - pi * real(i , wp)) > &
          1E1_wp * epsilon(1E0_wp)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed, zero number ', i+4, ' should ', &
            'be ', pi * real(i, wp), ' but is ', zeros_found(i+4), ' instead'
      end if
    end do
  end if
  write(unit_to_use , *) 'Zeros found were: ', zeros_found
  write(unit_to_use , *) 'Did the singlevar_proc_ptr test of ' &
      // 'find_all_zeros for the sin function pass? ', test_pass
end function test_svar_pptr_find_all_zeros

function test_svar_pptr_find_rel_optima(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  integer                           :: i
  real(wp), allocatable             :: zeros_found(:)
  real(wp), allocatable             :: optima_found(: , :)
  real(wp)                          :: grid(200)
  type(singlevar_proc_ptr)          :: svar_pptr
  procedure(univar_func), pointer   :: sin_ptr
  integer                           :: local_err_stat
  character(len=err_msg_len)        :: local_err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test use of a singlevar_proc_ptr ', &
      'in find_rel_optima.'
  sin_ptr => sinfunc
  call svar_pptr%set_proc_ptr(sin_ptr)
  grid = 1E1_wp * chebyshev_points(200) ! From -10 to 10
  call find_all_zeros( &
      zeros_found, &
      svar_pptr, &
      grid, &
      1E-15_wp, & ! toler
      local_err_stat, &
      local_err_msg )
  if (local_err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use , *) 'Test failed at zero-finding with message ' &
        // local_err_msg
  end if
  call find_rel_optima( &
      optima_found, &
      svar_pptr, &
      zeros_found, &
      1E-15_wp, & ! toler
      local_err_stat, &
      local_err_msg )
  if (local_err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use , *) 'Test failed with message ' // local_err_msg
  end if
  if (size(optima_found , 1) /= 6) then
    test_pass = .false.
    write(unit_to_use , *) 'Test failed, should have found 6 relative ', &
        'optima but found ', size(optima_found , 1), ' instead'
  else
    do i = 1, 6
      if (abs(svar_pptr%eval(optima_found(i , 1)) &
          - svar_pptr%eval(pi * real(2 * i - 7, wp) / 2E0_wp)) > &
          1E1_wp * epsilon(1E0_wp)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed, relative optimum number ', i, &
            ' should be ', pi * real(2 * i - 7, wp) / 2E0_wp, ' with val ', &
            svar_pptr%eval(pi * real(2 * i - 7, wp) / 2E0_wp), ' but is ', &
            optima_found(i , 1), ' with val ', &
            svar_pptr%eval(optima_found(i , 1)), ' instead'
      end if
    end do
  end if
  write(unit_to_use , *) 'Optima found were located at: ', optima_found(: , 1)
  write(unit_to_use , *) 'Optima values found were: ', optima_found(: , 2)
  write(unit_to_use , *) 'Did the singlevar_proc_ptr test of ' &
      // 'find_rel_optima for the sin function pass? ', test_pass
end function test_svar_pptr_find_rel_optima

end module test_singlevar_proc_ptr_mod