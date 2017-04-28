! test_falling_factorial.f90
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
! created on: 2016-02-11
! updated on: 2016-02-11
!
! A module of unit tests to check the computation of
! falling factorials.
  
module test_falling_factorial_mod

use set_precision, only : wp
use falling_factorial_mod, only : falling_factorial

implicit none

contains
  
function test_falling_factorial(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  write(unit_to_use , *) 'Test of falling factorial calculations'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_falling_factorial_small_n(unit_to_use)
  test_pass = test_pass .and. test_falling_factorial_large_n(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line

end function test_falling_factorial

function test_falling_factorial_small_n(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  real(wp)                    :: testfac
  real(wp)                    :: expected_factorial
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Check 3rd falling factorial of 2.3: '
  testfac = falling_factorial(2.3E0_wp, 3)
  expected_factorial = 2.3E0_wp * 1.3E0_wp * 0.3E0_wp
  write(unit_to_use , *) 'Falling factorial: ', &
      testfac, ' expected factorial: ', expected_factorial
  if (abs(testfac - expected_factorial) > 1E1_wp * epsilon(testfac)) then
    test_pass = .false.
    write(unit_to_use, *) 'Falling factorial test failed, difference of: ', &
        testfac - expected_factorial
  end if
end function test_falling_factorial_small_n

function test_falling_factorial_large_n(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  integer                     :: i
  real(wp)                    :: testfac
  real(wp)                    :: expected_factorial
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Check the 10th falling factorial of 11.1: '
  testfac = falling_factorial(11.1E0_wp, 10)
  expected_factorial = 11.1E0_wp
  do i = 1, 9
    expected_factorial = expected_factorial * (11.1E0_wp - real(i, wp))
  end do
  write(unit_to_use , *) 'Falling factorial: ', &
      testfac, ' expected factorial: ', expected_factorial
  if (abs(testfac - expected_factorial) > 1E1_wp * epsilon(testfac)) then
    test_pass = .false.
    write(unit_to_use, *) 'Falling factorial test failed, difference of: ', &
        testfac - expected_factorial
  end if
end function test_falling_factorial_large_n

end module test_falling_factorial_mod