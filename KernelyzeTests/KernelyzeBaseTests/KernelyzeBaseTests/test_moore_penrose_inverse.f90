! test_moore_penrose_inverse.f90
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
! created on: 2015-03-03
! updated on: 2015-04-22
! updated on: 2016-04-25 (to use kernel objects)
!
! A module of unit tests for the moore_penrose_inverse module.
  
module test_moore_penrose_inverse_mod

use set_precision, only : wp
use kernel_test_funcs_mod , only : exp_prod_kernel
use moore_penrose_inverse_mod, only  : moore_penrose_inverse

implicit none

contains
  
function test_moore_penrose_inverse(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  ! Test moore_penrose_inverse
  write(unit_to_use , *) 'Test of Moore-Penrose inverse functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_moore_penrose_inverse_identity(unit_to_use)
  test_pass = test_pass .and. test_moore_penrose_inverse_ones(unit_to_use)
  test_pass = test_pass .and. test_moore_penrose_inverse_exp_prod(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_moore_penrose_inverse

function test_moore_penrose_inverse_identity(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(5 , 5)             :: identity_mat, inverse_result
  real(wp) , dimension(5 , 5)             :: should_be_identity
  real(wp)                                :: toler
  integer                                 :: i
  ! Body
  test_pass = .true.
  identity_mat = 0E0_wp
  do i = 1, size(identity_mat , 1)
    identity_mat(i , i) = 1E0_wp
  end do
  toler = 20 * epsilon(1E0_wp)
  inverse_result = moore_penrose_inverse(identity_mat , toler)
  ! Output the test results
  write(unit_to_use , *) 'The Moore-Penrose inverse of the ' &
      // 'identity matrix is: '
  ! Use an unlimited format item to control output
  do i = 1, size(inverse_result, 1)
    write(unit_to_use , '(*(2x,F10.8))') inverse_result(i , : )
  end do
  ! Reset the identity_mat matrix; it was changed in moore_penrose_inverse
  identity_mat = 0E0_wp
  do i = 1, size(identity_mat , 1)
    identity_mat(i , i) = 1E0_wp
  end do
  should_be_identity = matmul(identity_mat , inverse_result)
  write(unit_to_use , *) 'The product of the identity matrix and ' &
      // 'its Moore-Penrose inverse is: '
  ! Use an unlimited format item to control output
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_mat - should_be_identity)) <= toler)
end function test_moore_penrose_inverse_identity

function test_moore_penrose_inverse_ones(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(6 , 6)             :: ones_mat, inverse_result
  real(wp) , dimension(6 , 6)             :: identity_mat, should_be_ones
  real(wp)                                :: toler
  integer                                 :: i
  ! Body
  test_pass = .true.
  ones_mat = 1E0_wp
  identity_mat = 0E0_wp
  do i = 1, size(identity_mat , 1)
    identity_mat(i , i) = 1E0_wp
  end do
  toler = 20 * epsilon(1E0_wp)
  inverse_result = moore_penrose_inverse(ones_mat , toler)
  ! Output the test results
  write(unit_to_use , *) 'The Moore-Penrose inverse of the ' &
      // 'ones matrix is: '
  ! Use an unlimited format item to control output
  do i = 1, size(inverse_result, 1)
    write(unit_to_use , '(*(2x,F10.8))') inverse_result(i , : )
  end do
  ! Reset the ones_mat matrix, since it was changed in moore_penrose_inverse
  ones_mat = 1E0_wp
  should_be_ones = matmul(matmul(ones_mat , inverse_result) , ones_mat)
  write(unit_to_use , *) 'The following should be all ones: '
  ! Use an unlimited format item to control output
  do i = 1, size(should_be_ones, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_ones(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the ones matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(ones_mat - should_be_ones)) <= toler)
end function test_moore_penrose_inverse_ones

function test_moore_penrose_inverse_exp_prod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(10 , 10)           :: exp_prod_mat, inverse_result
  real(wp) , dimension(10 , 10)           :: should_be_exp_prod, identity_mat
  integer                                 :: i, j
  ! Body
  test_pass = .true.
  identity_mat = 0E0_wp
  do i = 1, size(identity_mat , 1)
    identity_mat(i , i) = 1E0_wp
  end do
  do j = 1, size(exp_prod_mat, 2)
    do i = 1, size(exp_prod_mat, 1)
      exp_prod_mat(i, j) = exp_prod_kernel%eval( &
          real(i, wp) / 1E1_wp, real(j, wp) / 1E1_wp)
    end do
  end do
  ! Output the test exponential product matrix
  write(unit_to_use , *) 'The test exponential product matrix is: '
  ! Use an unlimited format item to control output
  do i = 1, size(inverse_result, 1)
    write(unit_to_use , '(*(2x,F16.14))') exp_prod_mat(i , : )
  end do
  ! Compute the Moore-Penrose generalized inverse
  inverse_result = moore_penrose_inverse(exp_prod_mat)
  ! Output the test results
  write(unit_to_use , *) 'The Moore-Penrose inverse of the ' &
      // 'test exponential product matrix is: '
  ! Use an unlimited format item to control output
  do i = 1, size(inverse_result, 1)
    write(unit_to_use , '(*(2x,F16.2))') inverse_result(i , : )
  end do
  ! Reset the exp_prod_mat matrix; it was changed in moore_penrose_inverse
  do j = 1, size(exp_prod_mat, 2)
    do i = 1, size(exp_prod_mat, 1)
      exp_prod_mat(i, j) = exp_prod_kernel%eval( &
          real(i, wp) / 1E1_wp, real(j, wp) / 1E1_wp)
    end do
  end do
  should_be_exp_prod = matmul(matmul(exp_prod_mat , inverse_result) , & 
      exp_prod_mat)
  write(unit_to_use , *) 'This should be the test exponential product ' &
      // 'matrix, but only approximately: '
  ! Use an unlimited format item to control output
  do i = 1, size(should_be_exp_prod, 1)
    write(unit_to_use , '(*(2x,F16.14))') should_be_exp_prod(i , : )
  end do
end function test_moore_penrose_inverse_exp_prod

end module test_moore_penrose_inverse_mod