! test_compute_a_matrix.f90
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
! created on: 2015-03-08
! updated on: 2015-04-22
! updated on: 2016-04-25 (to use kernel objects)
!
! A module of unit tests for the compute_a_matrix module.
  
module test_compute_a_matrix_mod

use set_precision, only : wp
use chebyshev_points_mod, only : chebyshev_points
use kernel_test_funcs_mod, only : &
    exp_prod_kernel, gaussian_kernel, cauchy_kernel
use compute_a_matrix_mod , only : compute_a_matrix

implicit none

contains

function test_compute_a_matrix(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of compute-A-matrix functionality: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_compute_a_matrix_exp_prod(unit_to_use)
  test_pass = test_pass .and. test_compute_a_matrix_gaussian(unit_to_use)
  test_pass = test_pass .and. test_compute_a_matrix_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_compute_a_matrix

function test_compute_a_matrix_exp_prod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  real(wp), dimension(9 , 9)  :: test_a_matrix, kernel_matrix
  real(wp), dimension(9 , 9)  :: identity_matrix, should_be_identity
  real(wp), dimension(9)      :: test_rho
  real(wp), parameter         :: toler = 1E-6_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A matrix computation for the ' &
      // 'exponential product kernel: '
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = chebyshev_points(size(test_rho))
  test_a_matrix = compute_a_matrix(exp_prod_kernel , test_rho)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = exp_prod_kernel%eval(test_rho(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_a_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product A * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A matrix computation for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_compute_a_matrix_exp_prod

function test_compute_a_matrix_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  real(wp), dimension(6 , 6)  :: test_a_matrix, kernel_matrix
  real(wp), dimension(6 , 6)  :: identity_matrix, should_be_identity
  real(wp), dimension(6)      :: test_rho
  real(wp), parameter         :: toler = 1E-6_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A matrix computation for the ' &
      // 'Gaussian kernel: '
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = chebyshev_points(size(test_rho))
  test_a_matrix = compute_a_matrix(gaussian_kernel , test_rho)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = gaussian_kernel%eval(test_rho(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_a_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product A * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A matrix computation for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_compute_a_matrix_gaussian

function test_compute_a_matrix_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)           :: unit_to_use
  ! Return value
  logical                       :: test_pass
  ! Local variables
  real(wp), dimension(8 , 8)  :: test_a_matrix, kernel_matrix
  real(wp), dimension(8 , 8)  :: identity_matrix, should_be_identity
  real(wp), dimension(8)      :: test_rho
  ! Note: greater margin for error given here than in other
  ! two kernel tests of compute_a_matrix
  real(wp), parameter         :: toler = 1E-4_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A matrix computation for the ' &
      // 'Cauchy kernel: '
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = chebyshev_points(size(test_rho))
  test_a_matrix = compute_a_matrix(cauchy_kernel , test_rho)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = cauchy_kernel%eval(test_rho(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_a_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product A * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A matrix computation for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_compute_a_matrix_cauchy

end module test_compute_a_matrix_mod