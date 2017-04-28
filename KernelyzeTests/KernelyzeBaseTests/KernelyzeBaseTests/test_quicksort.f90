! test_quicksort.f90
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
! created on: 2015-03-08
! updated on: 2015-04-22
!
! A module of unit tests for the quicksort module.
  
module test_quicksort_mod

use set_precision, only : wp
use quicksort , only : qsort

implicit none

contains

function test_quicksort(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of quicksort functionality: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_quicksort_already_ordered(unit_to_use)
  test_pass = test_pass .and. test_quicksort_reverse_ordered(unit_to_use)
  test_pass = test_pass .and. test_quicksort_one_out_of_order(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_quicksort

function test_quicksort_already_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)             :: unit_to_use
  ! Return value
  logical                         :: test_pass
  ! Local variables
  real(wp) , dimension(1000)      :: test_array
  real(wp)                        :: temp_vbl
  integer                         :: j
  ! Test the quicksort functionality
  do j = 1 , 1000
    ! The constants are not critical here, they just
    ! exist to be sure that I am working with reals of
    ! kind wp rather than integers
    test_array(j) = 1.2_wp * real(j, wp) + 0.3456_wp 
  end do
  call qsort(test_array, 1, 1000, temp_vbl)
  test_pass = .true.
  do j = 2, 1000
    if (test_array(j) < test_array(j - 1)) then
      test_pass = .false.
    end if
  end do
  write(unit_to_use , *) 'Did the already-ordered test pass? ', test_pass
end function test_quicksort_already_ordered

function test_quicksort_reverse_ordered(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)             :: unit_to_use
  ! Return value
  logical                         :: test_pass
  ! Local variables
  real(wp) , dimension(5000)      :: test_array
  real(wp)                        :: temp_vbl
  integer                         :: j
  ! Test the quicksort functionality
  do j = 1 , 5000
    ! The constants are not critical here, they just
    ! exist to be sure that I am working with reals of
    ! kind wp rather than integers
    test_array(j) = 2552_wp - (2.3_wp * real(j, wp) + 0.4567_wp)
  end do
  call qsort(test_array, 1, 5000, temp_vbl)
  test_pass = .true.
  do j = 2, 5000
    if (test_array(j) < test_array(j - 1)) then
      test_pass = .false.
    end if
  end do
  write(unit_to_use , *) 'Did the reverse-ordered test pass? ', test_pass
end function test_quicksort_reverse_ordered

function test_quicksort_one_out_of_order(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)             :: unit_to_use
  ! Return value
  logical                         :: test_pass
  ! Local variables
  real(wp) , dimension(4321)      :: test_array
  real(wp)                        :: temp_vbl
  integer                         :: j
  ! Test the quicksort functionality
  ! Fill the array in an ordered fashion
  do j = 1 , 4321
    ! The constants are not critical here, they just
    ! exist to be sure that I am working with reals of
    ! kind wp rather than integers
    test_array(j) = -1432_wp + (2.3_wp * real(j, wp) + 0.4567_wp)
  end do
  ! Now make the last element of the test array the smallest
  test_array(4321) = -1E5_wp
  call qsort(test_array, 1, 4321, temp_vbl)
  test_pass = .true.
  do j = 2, 4321
    if (test_array(j) < test_array(j - 1)) then
      test_pass = .false.
    end if
  end do
  write(unit_to_use , *) 'Did the one-out-of-order test pass? ', test_pass
end function test_quicksort_one_out_of_order

end module test_quicksort_mod