! test_all_num_opt_mats_asymm.f90
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
! created on: 2016-07-03 
! updated on: 2016-07-03
! updated on: 2016-11-11 (changed tolerance in Cauchy kernel test as
!             results changed marginally with an update to 
!             Intel Fortran 2016 version 4)
!
! A module of unit tests for the all_num_opt_mats_asymm module.
  
module test_all_num_opt_mats_asymm_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use worst_rho_and_gamma_mod, only : worst_rho_and_gamma
use all_num_opt_mats_asymm_mod , only : all_num_opt_mats_asymm

implicit none

contains

function test_all_num_opt_mats_asymm(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of all-num-optimal-matrices-asymmetric ' &
      // 'functionality: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. &
      test_all_num_opt_mats_asymm_exp_prod(unit_to_use)
  test_pass = test_pass .and. &
      test_all_num_opt_mats_asymm_gaussian(unit_to_use)
  test_pass = test_pass .and. &
      test_all_num_opt_mats_asymm_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_all_num_opt_mats_asymm

function test_all_num_opt_mats_asymm_exp_prod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_expprod)          :: exp_prod_kernel
  real(wp), dimension(9 , 9)    :: test_a_matrix, test_b_matrix
  real(wp), dimension(9 , 9)    :: identity_minus_rank_one
  real(wp), dimension(9 , 8)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(9 , 9)    :: c_matrix, prod_matrix, should_be_zeros
  real(wp), dimension(9)        :: test_rho, test_gamma
  real(wp), dimension(9)        :: test_coeffs_rho, test_coeffs_gamma
  real(wp), dimension(9)        :: nodes_rho, errors_at_nodes_rho
  real(wp), dimension(9)        :: b_coeffs_rho, b_nodes_rho
  real(wp), dimension(9)        :: b_errors_at_nodes_rho
  real(wp), dimension(9)        :: nodes_gamma, errors_at_nodes_gamma
  real(wp), dimension(9)        :: b_coeffs_gamma, b_nodes_gamma
  real(wp), dimension(9)        :: b_errors_at_nodes_gamma
  real(wp), dimension(8)        :: should_be_ones
  real(wp), parameter           :: toler = 1E-8_wp
  real(wp)                      :: discrep, ker_discrep
  real(wp)                      :: b_lowerbd
  integer                       :: i, j, num_iter, ker_num_iter
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp), dimension(201, 201) :: error_grid
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test all-num-opt-mat-asymm for the ' &
      // 'exponential product kernel: '
  call exp_prod_kernel%set_c(1E0_wp)
  call exp_prod_kernel%set_x_lb(0E0_wp)
  call exp_prod_kernel%set_x_ub(5E0_wp)
  call exp_prod_kernel%set_y_lb(-1E0_wp)
  call exp_prod_kernel%set_y_ub(0.5E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call worst_rho_and_gamma( &
      kernel_obj = exp_prod_kernel, &
      num_terms = size(test_coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      a_matrix_rho = test_a_matrix, &
      a_matrix_gamma = test_b_matrix, &
      nodes_rho = nodes_rho, &
      nodes_gamma = nodes_gamma, &
      errors_at_nodes_rho = errors_at_nodes_rho, &
      errors_at_nodes_gamma = errors_at_nodes_gamma, &
      b_coeffs_rho = b_coeffs_rho, &
      b_coeffs_gamma = b_coeffs_gamma, &
      b_nodes_rho = b_nodes_rho, &
      b_nodes_gamma = b_nodes_gamma, &
      b_errors_at_nodes_rho = b_errors_at_nodes_rho, &
      b_errors_at_nodes_gamma = b_errors_at_nodes_gamma, &
      discrepancy = discrep, &
      num_iter = num_iter, &
      ker_discrep = ker_discrep, &
      ker_num_iter = ker_num_iter, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Obtain the matrices used to form the numerically-optimal approximation
  call all_num_opt_mats_asymm( &
      kernel_obj = exp_prod_kernel, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      a_mat = test_a_matrix, & 
      b_mat = test_b_matrix, &
      d_vec = should_be_ones, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Output the vector that should be ones but note that numerical
  ! issues may cause elements of the vector to differ from one.
  ! Since the V and W matrices account for the possibility of such
  ! a departure from expectations, such a situation should not
  ! create any problems downstream.
  write(unit_to_use , *) 'The following vector of singular values ' &
      // 'should be all ones, but may be different due to ' &
      // 'numerical challenges: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,F10.8))') should_be_ones
  ! Compute I - coeffs * coeffs'
  do j = 1, size(identity_minus_rank_one, 2)
    do i = 1, size(identity_minus_rank_one, 1)
      identity_minus_rank_one(i , j) = &
          - test_coeffs_rho(i) * test_coeffs_gamma(j)
    end do
    identity_minus_rank_one(j , j) = identity_minus_rank_one(j , j) + 1E0_wp
  end do
  ! The computed V and A matrices should satisfy
  ! A' * (I - coeffs_rho * coeffs_gamma') * B == V * W'
  c_matrix = matmul(transpose(test_a_matrix) , &
      matmul(identity_minus_rank_one , test_b_matrix))
  prod_matrix = matmul(test_v_matrix, transpose(test_w_matrix))
  should_be_zeros =  c_matrix - prod_matrix
  ! Now output the matrix that should be zero, then
  ! check it against the zero matrix to be sure that each element
  ! is close to the corresponding element of the zero matrix.
  ! Use an unlimited format item to control output
  write(unit_to_use , *) &
      'The A^T * (I - coeffs_rho * coeffs_gamma^T) * B - V * W^T ' &
      // 'matrix should be zero and is: '
  do i = 1, size(should_be_zeros, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_zeros(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! zero.
  test_pass = test_pass .and. &
      (maxval(abs(should_be_zeros)) <= toler)
  ! Now check the results using the Borsuk lower bound:
  b_lowerbd = min( &
      minval( abs( b_errors_at_nodes_rho ) ), &
      minval( abs( b_errors_at_nodes_gamma ) ) )
  test_grid_x = exp_prod_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = exp_prod_kernel%cheb_pts(.false., size(test_grid_y))
  do j = 1, 201
    do i = 1, 201
      error_grid(i , j) = exp(test_grid_x(i) * test_grid_y(j)) - &
          dot_product( matmul( exp(test_rho * test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( exp(test_gamma * test_grid_y(j)) , &
                  test_w_matrix ) )
    end do
  end do
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write(unit_to_use , *) 'Borsuk errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_rho
  write(unit_to_use , *) 'Borsuk errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_gamma
  write(unit_to_use , *) 'Kernel Remez errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_rho
  write(unit_to_use , *) 'Kernel Remez errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_gamma
  write(unit_to_use , *) 'Prod of x and y kernel Remez errors at nodes: '
  write(unit_to_use , '(*(2x,e20.13))') &
      errors_at_nodes_rho * errors_at_nodes_gamma
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Largest error on test grid is: ', &
      maxval( abs(error_grid) )
  write(unit_to_use , *) 'Did the all-num-opt-mat computation for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_all_num_opt_mats_asymm_exp_prod

function test_all_num_opt_mats_asymm_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_gaussian)         :: gaussian_kernel
  real(wp), dimension(6 , 6)    :: test_a_matrix, test_b_matrix
  real(wp), dimension(6 , 6)    :: identity_minus_rank_one
  real(wp), dimension(6 , 5)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(6 , 6)    :: c_matrix, prod_matrix, should_be_zeros
  real(wp), dimension(6)        :: test_rho, test_gamma
  real(wp), dimension(6)        :: test_coeffs_rho, test_coeffs_gamma
  real(wp), dimension(6)        :: nodes_rho, errors_at_nodes_rho
  real(wp), dimension(6)        :: b_coeffs_rho, b_nodes_rho
  real(wp), dimension(6)        :: b_errors_at_nodes_rho
  real(wp), dimension(6)        :: nodes_gamma, errors_at_nodes_gamma
  real(wp), dimension(6)        :: b_coeffs_gamma, b_nodes_gamma
  real(wp), dimension(6)        :: b_errors_at_nodes_gamma
  real(wp), dimension(5)        :: should_be_ones
  real(wp), parameter           :: toler = 1E-8_wp
  real(wp)                      :: discrep, ker_discrep
  real(wp)                      :: b_lowerbd
  integer                       :: i, j, num_iter, ker_num_iter
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp), dimension(201, 201) :: error_grid
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test all-num-opt-mat-asymm for the ' &
      // 'Gaussian kernel: '
  call gaussian_kernel%set_c(1E0_wp)
  call gaussian_kernel%set_x_lb(-4E0_wp)
  call gaussian_kernel%set_x_ub(-1E0_wp)
  call gaussian_kernel%set_y_lb(-2E0_wp)
  call gaussian_kernel%set_y_ub(2E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call worst_rho_and_gamma( &
      kernel_obj = gaussian_kernel, &
      num_terms = size(test_coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      a_matrix_rho = test_a_matrix, &
      a_matrix_gamma = test_b_matrix, &
      nodes_rho = nodes_rho, &
      nodes_gamma = nodes_gamma, &
      errors_at_nodes_rho = errors_at_nodes_rho, &
      errors_at_nodes_gamma = errors_at_nodes_gamma, &
      b_coeffs_rho = b_coeffs_rho, &
      b_coeffs_gamma = b_coeffs_gamma, &
      b_nodes_rho = b_nodes_rho, &
      b_nodes_gamma = b_nodes_gamma, &
      b_errors_at_nodes_rho = b_errors_at_nodes_rho, &
      b_errors_at_nodes_gamma = b_errors_at_nodes_gamma, &
      discrepancy = discrep, &
      num_iter = num_iter, &
      ker_discrep = ker_discrep, &
      ker_num_iter = ker_num_iter, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Obtain the matrices used to form the numerically-optimal approximation
  call all_num_opt_mats_asymm( &
      kernel_obj = gaussian_kernel, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      a_mat = test_a_matrix, & 
      b_mat = test_b_matrix, &
      d_vec = should_be_ones, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Output the vector that should be ones but note that numerical
  ! issues may cause elements of the vector to differ from one.
  ! Since the V and W matrices account for the possibility of such
  ! a departure from expectations, such a situation should not
  ! create any problems downstream.
  write(unit_to_use , *) 'The following vector of singular values ' &
      // 'should be all ones, but may be different due to ' &
      // 'numerical challenges: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,F10.8))') should_be_ones
  ! Compute I - coeffs * coeffs'
  do j = 1, size(identity_minus_rank_one, 2)
    do i = 1, size(identity_minus_rank_one, 1)
      identity_minus_rank_one(i , j) = &
          - test_coeffs_rho(i) * test_coeffs_gamma(j)
    end do
    identity_minus_rank_one(j , j) = identity_minus_rank_one(j , j) + 1E0_wp
  end do
  ! The computed V and A matrices should satisfy
  ! A' * (I - coeffs_rho * coeffs_gamma') * B == V * W'
  c_matrix = matmul(transpose(test_a_matrix) , &
      matmul(identity_minus_rank_one , test_b_matrix))
  prod_matrix = matmul(test_v_matrix, transpose(test_w_matrix))
  should_be_zeros =  c_matrix - prod_matrix
  ! Now output the matrix that should be zero, then
  ! check it against the zero matrix to be sure that each element
  ! is close to the corresponding element of the zero matrix.
  ! Use an unlimited format item to control output
  write(unit_to_use , *) &
      'The A^T * (I - coeffs_rho * coeffs_gamma^T) * B - V * W^T ' &
      // 'matrix should be zero and is: '
  do i = 1, size(should_be_zeros, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_zeros(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! zero.
  test_pass = test_pass .and. &
      (maxval(abs(should_be_zeros)) <= toler)
  ! Now check the results using the Borsuk lower bound:
  b_lowerbd = min( &
      minval( abs( b_errors_at_nodes_rho ) ), &
      minval( abs( b_errors_at_nodes_gamma ) ) )
  test_grid_x = gaussian_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = gaussian_kernel%cheb_pts(.false., size(test_grid_y))
  do j = 1, 201
    do i = 1, 201
      error_grid(i , j) = &
          gaussian_kernel%eval(test_grid_x(i) , test_grid_y(j)) - &
          dot_product( &
              matmul( gaussian_kernel%eval_elt(test_rho , test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( gaussian_kernel%eval_elt(test_gamma , test_grid_y(j)) , &
                  test_w_matrix ) )
    end do
  end do
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write(unit_to_use , *) 'Borsuk errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_rho
  write(unit_to_use , *) 'Borsuk errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_gamma
  write(unit_to_use , *) 'Kernel Remez errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_rho
  write(unit_to_use , *) 'Kernel Remez errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_gamma
  write(unit_to_use , *) 'Prod of x and y kernel Remez errors at nodes: '
  write(unit_to_use , '(*(2x,e20.13))') &
      errors_at_nodes_rho * errors_at_nodes_gamma
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Largest error on test grid is: ', &
      maxval( abs(error_grid) )
  write(unit_to_use , *) 'Did the all-num-opt-mat computation for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_all_num_opt_mats_asymm_gaussian

function test_all_num_opt_mats_asymm_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_cauchy)           :: cauchy_kernel
  real(wp), dimension(8 , 8)    :: test_a_matrix, test_b_matrix
  real(wp), dimension(8 , 8)    :: identity_minus_rank_one
  real(wp), dimension(8 , 7)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(8 , 8)    :: c_matrix, prod_matrix, should_be_zeros
  real(wp), dimension(8)        :: test_rho, test_gamma
  real(wp), dimension(8)        :: test_coeffs_rho, test_coeffs_gamma
  real(wp), dimension(8)        :: nodes_rho, errors_at_nodes_rho
  real(wp), dimension(8)        :: b_coeffs_rho, b_nodes_rho
  real(wp), dimension(8)        :: b_errors_at_nodes_rho
  real(wp), dimension(8)        :: nodes_gamma, errors_at_nodes_gamma
  real(wp), dimension(8)        :: b_coeffs_gamma, b_nodes_gamma
  real(wp), dimension(8)        :: b_errors_at_nodes_gamma
  real(wp), dimension(7)        :: should_be_ones
  real(wp), parameter           :: toler = 4E-7_wp
  real(wp)                      :: discrep, ker_discrep
  real(wp)                      :: b_lowerbd
  integer                       :: i, j, num_iter, ker_num_iter
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp), dimension(201, 201) :: error_grid
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test all-num-opt-mat-asymm for the ' &
      // 'Cauchy kernel: '
  call cauchy_kernel%set_c(3E0_wp)
  call cauchy_kernel%set_alpha(1E0_wp)
  call cauchy_kernel%set_x_lb(0E0_wp)
  call cauchy_kernel%set_x_ub(5E0_wp)
  call cauchy_kernel%set_y_lb(-0.5E0_wp)
  call cauchy_kernel%set_y_ub(2E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call worst_rho_and_gamma( &
      kernel_obj = cauchy_kernel, &
      num_terms = size(test_coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      a_matrix_rho = test_a_matrix, &
      a_matrix_gamma = test_b_matrix, &
      nodes_rho = nodes_rho, &
      nodes_gamma = nodes_gamma, &
      errors_at_nodes_rho = errors_at_nodes_rho, &
      errors_at_nodes_gamma = errors_at_nodes_gamma, &
      b_coeffs_rho = b_coeffs_rho, &
      b_coeffs_gamma = b_coeffs_gamma, &
      b_nodes_rho = b_nodes_rho, &
      b_nodes_gamma = b_nodes_gamma, &
      b_errors_at_nodes_rho = b_errors_at_nodes_rho, &
      b_errors_at_nodes_gamma = b_errors_at_nodes_gamma, &
      discrepancy = discrep, &
      num_iter = num_iter, &
      ker_discrep = ker_discrep, &
      ker_num_iter = ker_num_iter, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Obtain the matrices used to form the numerically-optimal approximation
  call all_num_opt_mats_asymm( &
      kernel_obj = cauchy_kernel, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      coeffs_rho = test_coeffs_rho, &
      coeffs_gamma = test_coeffs_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      a_mat = test_a_matrix, & 
      b_mat = test_b_matrix, &
      d_vec = should_be_ones, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Output the vector that should be ones but note that numerical
  ! issues may cause elements of the vector to differ from one.
  ! Since the V and W matrices account for the possibility of such
  ! a departure from expectations, such a situation should not
  ! create any problems downstream.
  write(unit_to_use , *) 'The following vector of singular values ' &
      // 'should be all ones, but may be different due to ' &
      // 'numerical challenges: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,F10.8))') should_be_ones
  ! Compute I - coeffs * coeffs'
  do j = 1, size(identity_minus_rank_one, 2)
    do i = 1, size(identity_minus_rank_one, 1)
      identity_minus_rank_one(i , j) = &
          - test_coeffs_rho(i) * test_coeffs_gamma(j)
    end do
    identity_minus_rank_one(j , j) = identity_minus_rank_one(j , j) + 1E0_wp
  end do
  ! The computed V and A matrices should satisfy
  ! A' * (I - coeffs_rho * coeffs_gamma') * B == V * W'
  c_matrix = matmul(transpose(test_a_matrix) , &
      matmul(identity_minus_rank_one , test_b_matrix))
  prod_matrix = matmul(test_v_matrix, transpose(test_w_matrix))
  should_be_zeros =  c_matrix - prod_matrix
  ! Now output the matrix that should be zero, then
  ! check it against the zero matrix to be sure that each element
  ! is close to the corresponding element of the zero matrix.
  ! Use an unlimited format item to control output
  write(unit_to_use , *) &
      'The A^T * (I - coeffs_rho * coeffs_gamma^T) * B - V * W^T ' &
      // 'matrix should be zero and is: '
  do i = 1, size(should_be_zeros, 1)
    write(unit_to_use , '(*(2x,F10.8))') should_be_zeros(i , : )
  end do
  ! No element of the product should differ by more than "toler" from
  ! zero.
  test_pass = test_pass .and. &
      (maxval(abs(should_be_zeros)) <= toler)
  ! Now check the results using the Borsuk lower bound:
  b_lowerbd = min( &
      minval( abs( b_errors_at_nodes_rho ) ), &
      minval( abs( b_errors_at_nodes_gamma ) ) )
  test_grid_x = cauchy_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = cauchy_kernel%cheb_pts(.false., size(test_grid_y))
  do j = 1, 201
    do i = 1, 201
      error_grid(i , j) = &
          cauchy_kernel%eval(test_grid_x(i) , test_grid_y(j)) - &
          dot_product( &
              matmul( cauchy_kernel%eval_elt(test_rho , test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( cauchy_kernel%eval_elt(test_gamma , test_grid_y(j)) , &
                  test_w_matrix ) )
    end do
  end do
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write(unit_to_use , *) 'Borsuk errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_rho
  write(unit_to_use , *) 'Borsuk errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') b_errors_at_nodes_gamma
  write(unit_to_use , *) 'Kernel Remez errors at nodes over x: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_rho
  write(unit_to_use , *) 'Kernel Remez errors at nodes over y: '
  write(unit_to_use , '(*(2x,e20.13))') errors_at_nodes_gamma
  write(unit_to_use , *) 'Prod of x and y kernel Remez errors at nodes: '
  write(unit_to_use , '(*(2x,e20.13))') &
      errors_at_nodes_rho * errors_at_nodes_gamma
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Largest error on test grid is: ', &
      maxval( abs(error_grid) )
  write(unit_to_use , *) 'Did the all-num-opt-mat computation for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_all_num_opt_mats_asymm_cauchy

end module test_all_num_opt_mats_asymm_mod