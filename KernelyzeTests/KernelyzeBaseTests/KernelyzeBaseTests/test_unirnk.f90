! test_unirnk.f90
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
! created on: 2015-02-16
! updated on: 2015-04-22
!
! A module of unit tests for the unirnk module.
  
module test_unirnk_mod

use set_precision, only : wp
use m_unirnk, only  : unirnk
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none

contains
  
function test_unirnk(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Test of unirnk
  write(unit_to_use , *) 'Test of unirnk functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_unirnk_identical(unit_to_use)
  test_pass = test_pass .and. test_unirnk_ordered(unit_to_use)
  test_pass = test_pass .and. test_unirnk_reverse_ordered(unit_to_use)
  test_pass = test_pass .and. test_unirnk_with_nans(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_unirnk

function test_unirnk_identical(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(200)  :: test_value_array
  integer, dimension(200)   :: test_order_array
  integer                   :: num_unique
  ! Test of unirnk
  test_pass = .true.
  
  write(unit_to_use , *) 'First test on 200 identical elements: '
  test_value_array = 1.0E0_wp
  call unirnk(test_value_array, test_order_array, num_unique)
  write(unit_to_use , *) 'The number of unique entries, ', &
    'should be one: ', num_unique
  test_pass = test_pass .and. (num_unique == 1)
end function test_unirnk_identical

function test_unirnk_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(200)  :: test_value_array
  integer, dimension(200)   :: test_order_array
  integer                   :: num_unique, j
  ! Test of unirnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test on ordered array: '
  do j = 1,200
    test_value_array(j) = real(j, wp)
  end do
  call unirnk(test_value_array, test_order_array, num_unique)
  write(unit_to_use , *) 'The number of unique entries, ', &
    'should be 200: ', num_unique
  test_pass = test_pass .and. (num_unique == 200)
  write(unit_to_use , *) 'Should go from 1 to 200 in ascending ', &
    'order: '
  do j = 1,num_unique
    if (j > 1) then
      if (test_order_array(j) <= test_order_array(j-1)) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
end function test_unirnk_ordered

function test_unirnk_reverse_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(200)  :: test_value_array
  integer, dimension(200)   :: test_order_array
  integer                   :: num_unique, j
  ! Test of unirnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test on reverse-ordered array: '
  do j = 1,200
    test_value_array(j) = real(201 - j, wp)
  end do
  call unirnk(test_value_array, test_order_array, num_unique)
  write(unit_to_use , *) 'The number of unique entries, ', &
    'should be 200: ', num_unique
  test_pass = test_pass .and. (num_unique == 200)
  write(unit_to_use , *) 'Should go from 200 to 1 in descending ', &
    'order: '
  do j = 1,num_unique
    if (j > 1) then
      if (test_order_array(j) >= test_order_array(j-1)) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
end function test_unirnk_reverse_ordered

function test_unirnk_with_nans(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(200)  :: test_value_array
  integer, dimension(200)   :: test_order_array
  integer                   :: num_unique, j
  ! Test of unirnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test handling of quiet NaNs: '
  ! Should set all odd entries to a quiet NaN.
  ! The result is that 100 entries are quiet NaN, while
  ! the remaining 100 entries are distinct numbers.
  do j = 1,200
    test_value_array(j) = real(201 - j, wp)
  end do
  do j = 1, 200, 2
    test_value_array(j) = ieee_value(test_value_array(j), ieee_quiet_nan)
  end do
  call unirnk(test_value_array, test_order_array, num_unique)
  write(unit_to_use , *) 'The number of unique entries, ', &
    'should be 200: ', num_unique
  test_pass = test_pass .and. (num_unique == 200)
  write(unit_to_use , *) 'Should be evens from 200 to 2 in descending ', &
    'order (the non-NaNs), then the odds from 1 to 199 in ascending ', &
    'order (the NaNs): '
  do j = 1,num_unique
    if (j <= 100) then
      if (test_order_array(j) /= (202 - (2 * j)) ) then
        test_pass = .false.
      end if
    else
      if (test_order_array(j) /= ((2 * j) - 201) ) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
end function test_unirnk_with_nans

end module test_unirnk_mod