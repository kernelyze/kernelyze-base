! test_compute_a_and_b_matrices.f90
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
! created on: 2016-07-02
! updated on: 2016-07-02
! updated on: 2016-07-03 (continued toward a first working version)
!
! A module of unit tests for the compute_a_and_b_matrices module.
  
module test_compute_a_and_b_matrices_mod

use set_precision, only : wp
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use compute_a_and_b_matrices_mod , only : compute_a_and_b_matrices

implicit none

contains

function test_compute_a_and_b_matrices(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of compute-A-and-B-matrices functionality: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_a_and_b_matrices_exp_prod(unit_to_use)
  test_pass = test_pass .and. test_a_and_b_matrices_gaussian(unit_to_use)
  test_pass = test_pass .and. test_a_and_b_matrices_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_compute_a_and_b_matrices

function test_a_and_b_matrices_exp_prod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_expprod)        :: exp_prod_kernel
  real(wp), dimension(9 , 9)  :: test_a_matrix, test_b_matrix, kernel_matrix
  real(wp), dimension(9 , 9)  :: identity_matrix, should_be_identity
  real(wp), dimension(9)      :: test_rho, test_gamma
  real(wp), parameter         :: toler = 1E-6_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A and B matrix computations for the ' &
      // 'exponential product kernel: '
  call exp_prod_kernel%set_c(1E0_wp)
  call exp_prod_kernel%set_x_lb(0E0_wp)
  call exp_prod_kernel%set_x_ub(5E0_wp)
  call exp_prod_kernel%set_y_lb(-1E0_wp)
  call exp_prod_kernel%set_y_ub(0.5E0_wp)
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = exp_prod_kernel%cheb_pts(.false., size(test_rho))
  test_gamma = exp_prod_kernel%cheb_pts(.true., size(test_gamma))
  call compute_a_and_b_matrices( &
      exp_prod_kernel, &
      test_rho, &
      test_gamma, &
      test_a_matrix, &
      test_b_matrix)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = exp_prod_kernel%eval(test_gamma(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_b_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product B * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A and B matrix computations for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_a_and_b_matrices_exp_prod

function test_a_and_b_matrices_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_gaussian)       :: gaussian_kernel
  real(wp), dimension(6 , 6)  :: test_a_matrix, test_b_matrix, kernel_matrix
  real(wp), dimension(6 , 6)  :: identity_matrix, should_be_identity
  real(wp), dimension(6)      :: test_rho, test_gamma
  real(wp), parameter         :: toler = 1E-6_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A and B matrix computations for the ' &
      // 'Gaussian kernel: '
  call gaussian_kernel%set_c(1E0_wp)
  call gaussian_kernel%set_x_lb(-4E0_wp)
  call gaussian_kernel%set_x_ub(-1E0_wp)
  call gaussian_kernel%set_y_lb(-2E0_wp)
  call gaussian_kernel%set_y_ub(2E0_wp)
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = gaussian_kernel%cheb_pts(.false., size(test_rho))
  test_gamma = gaussian_kernel%cheb_pts(.true., size(test_gamma))
  call compute_a_and_b_matrices( &
      gaussian_kernel, &
      test_rho, &
      test_gamma, &
      test_a_matrix, &
      test_b_matrix)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = gaussian_kernel%eval(test_gamma(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_b_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product B * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A and B matrix computations for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_a_and_b_matrices_gaussian

function test_a_and_b_matrices_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)           :: unit_to_use
  ! Return value
  logical                       :: test_pass
  ! Local variables
  type(kernel_cauchy)         :: cauchy_kernel
  real(wp), dimension(8 , 8)  :: test_a_matrix, test_b_matrix, kernel_matrix
  real(wp), dimension(8 , 8)  :: identity_matrix, should_be_identity
  real(wp), dimension(8)      :: test_rho, test_gamma
  ! Note: greater margin for error given here than in other
  ! two kernel tests of compute_a_matrix
  real(wp), parameter         :: toler = 1E-4_wp
  integer                     :: i, j
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test the A and B matrix computations for the ' &
      // 'Cauchy kernel: '
  call cauchy_kernel%set_c(3E0_wp)
  call cauchy_kernel%set_alpha(1E0_wp)
  call cauchy_kernel%set_x_lb(0E0_wp)
  call cauchy_kernel%set_x_ub(5E0_wp)
  call cauchy_kernel%set_y_lb(-0.5E0_wp)
  call cauchy_kernel%set_y_ub(2E0_wp)
  identity_matrix = 0E0_wp
  do j = 1, size(identity_matrix, 1)
    identity_matrix(j , j) = 1E0_wp
  end do
  test_rho = cauchy_kernel%cheb_pts(.false., size(test_rho))
  test_gamma = cauchy_kernel%cheb_pts(.true., size(test_gamma))
  call compute_a_and_b_matrices( &
      cauchy_kernel, &
      test_rho, &
      test_gamma, &
      test_a_matrix, &
      test_b_matrix)
  do j = 1, size(kernel_matrix , 2)
    do i = 1, size(kernel_matrix , 1)
      kernel_matrix(i , j) = cauchy_kernel%eval(test_gamma(i) , test_rho(j))
    end do
  end do
  should_be_identity = matmul(test_b_matrix , matmul(kernel_matrix , &
      transpose(test_a_matrix)))
  ! Now output the matrix that should be the identity, then
  ! check it against the identity matrix to be sure that each element
  ! is close to the corresponding element of the identity
  ! Use an unlimited format item to control output
  write(unit_to_use , *) 'The product B * kernel_mat * A^T should ' &
      // 'be the identity and is: '
  do i = 1, size(should_be_identity, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_identity(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! the identity matrix's corresponding element.
  test_pass = test_pass .and. &
      (maxval(abs(identity_matrix - should_be_identity)) <= toler)
  write(unit_to_use , *) 'Did the A and B matrix computations for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_a_and_b_matrices_cauchy

end module test_compute_a_and_b_matrices_mod