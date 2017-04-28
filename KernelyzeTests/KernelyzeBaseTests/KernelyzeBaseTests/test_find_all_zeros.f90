! test_find_all_zeros.f90
!
! Copyright (c) 2015, 2016 by Kernelyze LLC
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
! created on: 2015-02-16
! updated on: 2015-04-22
! updated on: 2016-04-25 (to use singlevar_funcs)
!
! A module of unit tests for the find_all_zeros module.
  
module test_find_all_zeros_mod

use set_precision, only : wp
use intrinsic_test_funcs_mod, only : sin_func, exp_func, log_func
use find_all_zeros_mod, only  : find_all_zeros

implicit none

contains
  
function test_find_all_zeros(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Test find_all_zeros
  write(unit_to_use , *) 'Test of find-all-zeros functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_find_all_zeros_sin(unit_to_use)
  test_pass = test_pass .and. test_find_all_zeros_log(unit_to_use)
  test_pass = test_pass .and. test_find_all_zeros_exp(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_find_all_zeros

function test_find_all_zeros_sin(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(21)                  :: grid_vec
  real(wp) , dimension(:) , allocatable     :: res_zeros
  integer                                   :: j
  ! Body
  test_pass = .true.
  do j = 1,21
    grid_vec(j) = -10.0E0_wp + DBLE(j - 1) 
  end do
  call find_all_zeros(res_zeros, sin_func, grid_vec, 1.0E-15_wp)
  write(unit_to_use , *) 'The zeros of sin(x) from -10.0 to 10.0 are ' &
      // '(should be seven of them): '
  do j = 1,size(res_zeros)
    write(unit_to_use , *) res_zeros(j)
  end do
  test_pass = test_pass .and. (size(res_zeros) == 7)
end function test_find_all_zeros_sin

function test_find_all_zeros_log(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(21)                  :: grid_vec
  real(wp) , dimension(:) , allocatable     :: res_zeros
  integer                                   :: j
  ! Body
  test_pass = .true.
  do j = 1,21
    grid_vec(j) = 0.5_wp + 0.5_wp * DBLE(j - 1) 
  end do
  call find_all_zeros(res_zeros, log_func, grid_vec, 1.0E-15_wp)
  write(unit_to_use , *) 'The zeros of log(x) from 0.5 to 10.5 are ' &
      // '(should be one of them): '
  do j = 1,size(res_zeros)
    write(unit_to_use , *) res_zeros(j)
  end do
  test_pass = test_pass .and. (size(res_zeros) == 1)
end function test_find_all_zeros_log

function test_find_all_zeros_exp(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(21)                  :: grid_vec
  real(wp) , dimension(:) , allocatable     :: res_zeros
  integer                                   :: j
  ! Body
  test_pass = .true.
  do j = 1,21
    grid_vec(j) = -10.0E0_wp + DBLE(j - 1) 
  end do
  call find_all_zeros(res_zeros, exp_func, grid_vec, 1.0E-15_wp)
  ! Of course, the exponential function has no zeros on the 
  ! real line.
  write(unit_to_use , *) 'The zeros of exp(x) from -10.0 to 10.0 are ' &
      // '(should be none): '
  do j = 1,size(res_zeros)
    write(unit_to_use , *) res_zeros(j)
  end do
  test_pass = test_pass .and. (size(res_zeros) == 0)
end function test_find_all_zeros_exp

end module test_find_all_zeros_mod