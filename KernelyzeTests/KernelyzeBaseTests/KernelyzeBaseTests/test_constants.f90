! test_constants.f90
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
! created on: 2015-02-23
! updated on: 2015-02-23
!
! A module of unit tests for the constants module.
  
module test_constants_mod

use set_precision, only : wp
use constants_mod, only : pi, exp_of_one, alloc_errmsg_len

implicit none

contains
  
function test_constants(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Test of constants
  write(unit_to_use , *) 'Test of constants'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_constants_pi(unit_to_use)
  test_pass = test_pass .and. test_constants_exp(unit_to_use)
  test_pass = test_pass .and. test_constants_alloc_errmsg_len(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_constants

function test_constants_pi(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variable
  real(wp)                  :: value_at_pi
  ! Body
  test_pass = .true.
  
  write(unit_to_use , *) 'Test the constant representing pi: '
  
  value_at_pi = sin(pi)
  write(unit_to_use , *) 'sin(pi) should be zero and is: ', value_at_pi
  test_pass = test_pass .and. (abs(value_at_pi) < epsilon(1.0_wp))
  
  value_at_pi = cos(pi)
  write(unit_to_use , *) 'cos(pi) should be -1 and is: ', value_at_pi
  test_pass = test_pass .and. (abs(value_at_pi + 1.0_wp) < epsilon(1.0_wp))

end function test_constants_pi

function test_constants_exp(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variable
  real(wp)                  :: value_at_exp
  ! Body
  test_pass = .true.
  
  write(unit_to_use , *) 'Test the constant representing exp(1): '
  
  value_at_exp = log(exp_of_one)
  write(unit_to_use , *) 'log(exp(1)) should be one and is: ', value_at_exp
  test_pass = test_pass .and. (abs(value_at_exp - 1.0_wp) < epsilon(1.0_wp))
  
end function test_constants_exp

function test_constants_alloc_errmsg_len(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Body
  test_pass = .true.
  
  write(unit_to_use , *) 'Test the constant representing ' // &
      'the length of any allocation error message: '
  
  write(unit_to_use , *) 'alloc_errmsg_len should be 255 and is: ', & 
      alloc_errmsg_len
  test_pass = test_pass .and. (alloc_errmsg_len == 255)
  
end function test_constants_alloc_errmsg_len

end module test_constants_mod