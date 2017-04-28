! test_mrgrnk.f90
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
! created on: 2015-03-14
! updated on: 2015-04-22
!
! A module of unit tests for the mrgrnk module.
  
module test_mrgrnk_mod

use set_precision, only : wp
use mrgrnk_mod, only  : mrgrnk
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none

contains
  
function test_mrgrnk(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Test of unirnk
  write(unit_to_use , *) 'Test of mrgrnk functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_mrgrnk_identical(unit_to_use)
  test_pass = test_pass .and. test_mrgrnk_ordered(unit_to_use)
  test_pass = test_pass .and. test_mrgrnk_reverse_ordered(unit_to_use)
  test_pass = test_pass .and. test_mrgrnk_with_nans(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_mrgrnk

function test_mrgrnk_identical(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(51)  :: test_value_array
  integer, dimension(51)   :: test_order_array
  integer                   :: j
  ! Test of mrgrnk
  test_pass = .true.
  
  write(unit_to_use , *) 'First test on 51 identical elements: '
  test_value_array = 1.0E0_wp
  call mrgrnk(test_value_array, test_order_array)
  do j = 2,size(test_order_array)
    test_pass = test_pass .and. &
        (test_value_array(test_order_array(j)) &
        >= test_value_array(test_order_array(j-1)))
  end do
  if (test_pass) then
    write(unit_to_use , *) 'Test on identical elements passed.'
  else 
    write(unit_to_use , *) 'Test on identical elements failed.'
  end if
end function test_mrgrnk_identical

function test_mrgrnk_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(20)   :: test_value_array
  integer, dimension(20)    :: test_order_array
  integer                   :: j
  ! Test of mrgrnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test on ordered array: '
  do j = 1,size(test_value_array)
    test_value_array(j) = real(j, wp)
  end do
  call mrgrnk(test_value_array, test_order_array)
  write(unit_to_use , *) 'Should go from 1 to 20 in ascending ', &
    'order: '
  do j = 1,size(test_order_array)
    if (j > 1) then
      if (test_order_array(j) <= test_order_array(j-1)) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
  if (test_pass) then
    write(unit_to_use , *) 'Test on ordered elements passed.'
  else 
    write(unit_to_use , *) 'Test on ordered elements failed.'
  end if
end function test_mrgrnk_ordered

function test_mrgrnk_reverse_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(27)   :: test_value_array
  integer, dimension(27)    :: test_order_array
  integer                   :: j
  ! Test of mrgrnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test on reverse-ordered array: '
  do j = 1,size(test_value_array)
    test_value_array(j) = real(201 - j, wp)
  end do
  call mrgrnk(test_value_array, test_order_array)
  write(unit_to_use , *) 'Should go from 27 to 1 in descending ', &
    'order: '
  do j = 1,size(test_order_array)
    if (j > 1) then
      if (test_order_array(j) >= test_order_array(j-1)) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
  if (test_pass) then
    write(unit_to_use , *) 'Test on reverse-ordered elements passed.'
  else 
    write(unit_to_use , *) 'Test on reverse-ordered elements failed.'
  end if
end function test_mrgrnk_reverse_ordered

function test_mrgrnk_with_nans(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(19)   :: test_value_array
  integer, dimension(19)    :: test_order_array
  integer                   :: j
  ! Test of mrgrnk
  test_pass = .true.
  
  write(unit_to_use , *) 'Test handling of quiet NaNs: '
  ! Set all odd entries to a quiet NaN
  do j = 1,size(test_value_array)
    test_value_array(j) = real(201 - j, wp)
  end do
  do j = 1, size(test_value_array), 2
    test_value_array(j) = ieee_value(test_value_array(j), ieee_quiet_nan)
  end do
  call mrgrnk(test_value_array, test_order_array)
  write(unit_to_use , *) 'Should be evens from 18 to 2 in descending ', &
    'order (the non-NaNs), then the odds from 1 to 19 in ascending ', &
    'order (the NaNs): '
  do j = 1,size(test_order_array)
    if (j <= 9) then
      if (test_order_array(j) /= (20 - (2 * j)) ) then
        test_pass = .false.
      end if
    else
      if (test_order_array(j) /= ((2 * j) - 19) ) then
        test_pass = .false.
      end if
    end if
    write(unit_to_use , *) test_order_array(j)
  end do
  if (test_pass) then
    write(unit_to_use , *) 'Test of NaN handling passed.'
  else 
    write(unit_to_use , *) 'Test of NaN handling failed.'
  end if
end function test_mrgrnk_with_nans

end module test_mrgrnk_mod