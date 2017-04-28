! test_binomial_coeffs.f90
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
! created on: 2016-02-06
! updated on: 2016-02-06
!
! A module of unit tests to check the computation of
! binomial coefficients.
  
module test_binomial_coeffs_mod

use set_precision, only : wp
use binomial_coeffs_mod, only : binomial_coeffs

implicit none

contains
  
function test_binomial_coeffs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  write(unit_to_use , *) 'Test of binomial coefficient calculations'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_binomial_coeffs_small_n(unit_to_use)
  test_pass = test_pass .and. test_binomial_coeffs_large_n(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line

end function test_binomial_coeffs

function test_binomial_coeffs_small_n(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  integer                     :: i
  real(wp)                    :: bin_coeffs(6)
  real(wp)                    :: expected_coeffs(6)
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Check the n == 5 binomial coefficients: '
  bin_coeffs = binomial_coeffs(5)
  expected_coeffs(1) = 1E0_wp
  expected_coeffs(2) = 5E0_wp
  expected_coeffs(3) = 10E0_wp
  expected_coeffs(4) = 10E0_wp
  expected_coeffs(5) = 5E0_wp
  expected_coeffs(6) = 1E0_wp
  do i = 1, size(bin_coeffs)
    write(unit_to_use , *) 'Binomial coeff: ', &
        bin_coeffs(i), ' expected coeff: ', expected_coeffs(i)
    if (abs(bin_coeffs(i) - expected_coeffs(i)) > epsilon(bin_coeffs(i))) then
      test_pass = .false.
      write(unit_to_use, *) 'Binomial coeff test failed, difference of: ', &
          bin_coeffs(i) - expected_coeffs(i)
    end if
  end do
end function test_binomial_coeffs_small_n

function test_binomial_coeffs_large_n(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  integer                     :: i
  real(wp)                    :: prev_bin_coeffs(20)
  real(wp)                    :: bin_coeffs(21)
  real(wp)                    :: expected_coeffs(21)
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Check the n == 20 binomial coefficients: '
  prev_bin_coeffs = binomial_coeffs(19)
  bin_coeffs = binomial_coeffs(20)
  ! Use the recurrence relation (a la Pascal's triangle) to get
  ! the expected coefficients:
  expected_coeffs(1) = 1E0_wp
  expected_coeffs(size(expected_coeffs)) = 1E0_wp
  do i = 2, size(expected_coeffs) - 1
    expected_coeffs(i) = prev_bin_coeffs(i) + prev_bin_coeffs(i - 1)
  end do
  do i = 1, size(bin_coeffs)
    write(unit_to_use , *) 'Binomial coeff: ', &
        bin_coeffs(i), ' expected coeff: ', expected_coeffs(i)
    if (abs(bin_coeffs(i) - expected_coeffs(i)) > epsilon(bin_coeffs(i))) then
      test_pass = .false.
      write(unit_to_use, *) 'Binomial coeff test failed, difference of: ', &
          bin_coeffs(i) - expected_coeffs(i)
    end if
  end do
end function test_binomial_coeffs_large_n

end module test_binomial_coeffs_mod