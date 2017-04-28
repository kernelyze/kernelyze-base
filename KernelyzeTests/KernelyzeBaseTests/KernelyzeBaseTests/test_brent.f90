! test_brent.f90
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
! A module of unit tests for the brent module.
  
module test_brent_mod

use set_precision, only : wp
use brent_mod, only  : zeroin, fmin
use constants_mod, only : pi
use intrinsic_test_funcs_mod, only : sin_func
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

implicit none

contains
  
function test_brent(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  write(unit_to_use , *) 'Test of Brent-related functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_brent_nice(unit_to_use)
  test_pass = test_pass .and. test_brent_no_zero(unit_to_use)
  test_pass = test_pass .and. test_brent_nice_min(unit_to_use)
  test_pass = test_pass .and. test_brent_large_intvl_min(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line

end function test_brent

function test_brent_nice(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  ! Local variables
  real(wp)  :: a, b, tol, zero_x
  
  ! Test zeroin
  write(unit_to_use , *) 'Check a nice problem: '
  a = -2.5E-1_wp
  b = 2.5E-1_wp
  tol = 1.0E-15_wp
  zero_x = zeroin(a, b, sin_func, tol)
  write(unit_to_use , *) &
    'The zero of sin(x) between -0.25 and 0.25, should be 0: ', zero_x
  test_pass = .not. ieee_is_nan(zero_x)
end function test_brent_nice

function test_brent_no_zero(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  ! Local variables
  real(wp)  :: a, b, tol, zero_x
  
  ! Test zeroin
  write(unit_to_use , *) &
    'Check a problem without a zero: '
  a = 1E-3_wp
  b = 2E-3_wp
  tol = 1.0E-15_wp
  zero_x = zeroin(a, b, sin_func, tol)
  write(unit_to_use , *) &
    'sin(x) has no zero between 0.001 and 0.002, ', &
    'this should be NaN: ', zero_x
  ! This test is passed if zero_x is NaN
  test_pass = ieee_is_nan(zero_x)
end function test_brent_no_zero

function test_brent_nice_min(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  ! Local variables
  real(wp)  :: a, b, tol, min_x, min_f
  
  ! Test fmin
  write(unit_to_use , *) &
    'Find a minimum over a nice interval: '
  a = -pi
  b = 0E0_wp
  tol = 1.0E-15_wp
  call fmin(a, b, sin_func, tol, min_x, min_f)
  write(unit_to_use , *) 'The minimum of sin(x) between ', &
    '-pi and zero is: ', min_f
  write(unit_to_use , *) 'This minimum is achieved at: ', min_x
  test_pass = .not.( &
      ieee_is_nan(min_f) .or. ieee_is_nan(min_x))
end function test_brent_nice_min

function test_brent_large_intvl_min(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  ! Local variables
  real(wp)  :: a, b, tol, min_x, min_f
  
  ! Test fmin
  write(unit_to_use , *) &
    'Find a minimum over a large interval: '
  a = 0E0_wp
  b = 1E2_wp ! 100.0
  tol = 1.0E-15_wp
  call fmin(a, b, sin_func, tol, min_x, min_f)
  write(unit_to_use , *) 'The minimum of sin(x) between ', &
    '0 and 100 is: ', min_f
  write(unit_to_use , *) 'This minimum is achieved at: ', min_x
  test_pass = .not.( &
      ieee_is_nan(min_f) .or. ieee_is_nan(min_x))
end function test_brent_large_intvl_min

end module test_brent_mod