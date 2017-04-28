! test_worst_rho_and_gamma.f90
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
!
! A module of unit tests for the worst_rho_and_gamma module.
  
module test_worst_rho_and_gamma_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use worst_rho_and_gamma_mod, only : worst_rho_and_gamma

implicit none

contains
  
function test_worst_rho_and_gamma(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality finding the rho ' // &
      'and gamma that maximize the Borsuk lower bound for a given kernel ' // &
      '(which need not be symmetric): '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_worst_rho_and_gamma_expprod(unit_to_use)
  test_pass = test_pass .and. test_worst_rho_and_gamma_gaussian(unit_to_use)
  test_pass = test_pass .and. test_worst_rho_and_gamma_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_worst_rho_and_gamma

function test_worst_rho_and_gamma_expprod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_expprod)        :: exp_prod_kernel
  real(wp)                    :: discrepancy, ker_discrep
  real(wp) , dimension(9)     :: rho_vec, gamma_vec
  real(wp) , dimension(9)     :: coeffs_rho, coeffs_gamma
  real(wp) , dimension(9)     :: nodes_rho, nodes_gamma
  real(wp) , dimension(9)     :: errors_at_nodes_rho, errors_at_nodes_gamma
  real(wp) , dimension(9)     :: b_coeffs_rho, b_coeffs_gamma
  real(wp) , dimension(9)     :: b_nodes_rho, b_nodes_gamma
  real(wp) , dimension(9)     :: b_errors_at_nodes_rho
  real(wp) , dimension(9)     :: b_errors_at_nodes_gamma
  real(wp) , dimension(9,9)   :: a_matrix_rho, a_matrix_gamma
  integer                     :: j , num_iter, ker_num_iter
  integer                     :: err_stat
  character(len=err_msg_len)  :: err_msg
  ! Body
  test_pass = .true.
  call exp_prod_kernel%set_c(1E0_wp)
  call exp_prod_kernel%set_x_lb(0E0_wp)
  call exp_prod_kernel%set_x_ub(5E0_wp)
  call exp_prod_kernel%set_y_lb(-1E0_wp)
  call exp_prod_kernel%set_y_ub(0.5E0_wp)
  ! Test the worst-rho-vector functionality
  call worst_rho_and_gamma( &
      kernel_obj = exp_prod_kernel, &
      num_terms = size(coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = rho_vec, &
      gamma_vec = gamma_vec, &
      coeffs_rho = coeffs_rho, &
      coeffs_gamma = coeffs_gamma, &
      a_matrix_rho = a_matrix_rho, &
      a_matrix_gamma = a_matrix_gamma, &
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
      discrepancy = discrepancy, &
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
  write(unit_to_use, *) 'Worst rho and gamma calculation, ', &
      'exponential product kernel'
  write(unit_to_use, *) 'x range from ', exp_prod_kernel%get_x_lb(), &
      ' to ', exp_prod_kernel%get_x_ub()
  write(unit_to_use, *) 'y range from ', exp_prod_kernel%get_y_lb(), &
      ' to ', exp_prod_kernel%get_y_ub()
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'kernel discrepancy: ', ker_discrep
  write(unit_to_use, *) 'number of iterations: ', num_iter
  write(unit_to_use, *) 'kernel number of iterations: ', ker_num_iter
  write(unit_to_use, *) 'rho vector: ', rho_vec
  write(unit_to_use, *) 'b_nodes_rho: ', b_nodes_rho
  write(unit_to_use, *) 'b_errors_at_nodes_rho: ', b_errors_at_nodes_rho
  write(unit_to_use, *) 'gamma vector: ', gamma_vec
  write(unit_to_use, *) 'b_nodes_gamma: ', b_nodes_gamma
  write(unit_to_use, *) 'b_errors_at_nodes_gamma: ', b_errors_at_nodes_gamma
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' rho: ' , rho_vec(j) , &
        ' coeff: ' , coeffs_rho(j) , ' nodes: ', nodes_rho(j) , &
        ' error at node: ' , errors_at_nodes_rho(j)
  end do
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' gamma: ' , gamma_vec(j) , &
        ' coeff: ' , coeffs_gamma(j) , ' nodes: ', nodes_gamma(j) , &
        ' error at node: ' , errors_at_nodes_gamma(j)
  end do
end function test_worst_rho_and_gamma_expprod

function test_worst_rho_and_gamma_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_gaussian)       :: gaussian_kernel
  real(wp)                    :: discrepancy, ker_discrep
  real(wp) , dimension(6)     :: rho_vec, gamma_vec
  real(wp) , dimension(6)     :: coeffs_rho, coeffs_gamma
  real(wp) , dimension(6)     :: nodes_rho, nodes_gamma
  real(wp) , dimension(6)     :: errors_at_nodes_rho, errors_at_nodes_gamma
  real(wp) , dimension(6)     :: b_coeffs_rho, b_coeffs_gamma
  real(wp) , dimension(6)     :: b_nodes_rho, b_nodes_gamma
  real(wp) , dimension(6)     :: b_errors_at_nodes_rho
  real(wp) , dimension(6)     :: b_errors_at_nodes_gamma
  real(wp) , dimension(6,6)   :: a_matrix_rho, a_matrix_gamma
  integer                     :: j , num_iter, ker_num_iter
  integer                     :: err_stat
  character(len=err_msg_len)  :: err_msg
  ! Body
  test_pass = .true.
  call gaussian_kernel%set_c(1E0_wp)
  call gaussian_kernel%set_x_lb(-4E0_wp)
  call gaussian_kernel%set_x_ub(-1E0_wp)
  call gaussian_kernel%set_y_lb(-2E0_wp)
  call gaussian_kernel%set_y_ub(2E0_wp)
  ! Test the worst-rho-vector functionality
  call worst_rho_and_gamma( &
      kernel_obj = gaussian_kernel, &
      num_terms = size(coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = rho_vec, &
      gamma_vec = gamma_vec, &
      coeffs_rho = coeffs_rho, &
      coeffs_gamma = coeffs_gamma, &
      a_matrix_rho = a_matrix_rho, &
      a_matrix_gamma = a_matrix_gamma, &
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
      discrepancy = discrepancy, &
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
  write(unit_to_use, *) 'Worst rho and gamma calculation, ', &
      'Gaussian kernel'
  write(unit_to_use, *) 'x range from ', gaussian_kernel%get_x_lb(), &
      ' to ', gaussian_kernel%get_x_ub()
  write(unit_to_use, *) 'y range from ', gaussian_kernel%get_y_lb(), &
      ' to ', gaussian_kernel%get_y_ub()
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'kernel discrepancy: ', ker_discrep
  write(unit_to_use, *) 'number of iterations: ', num_iter
  write(unit_to_use, *) 'kernel number of iterations: ', ker_num_iter
  write(unit_to_use, *) 'rho vector: ', rho_vec
  write(unit_to_use, *) 'b_nodes_rho: ', b_nodes_rho
  write(unit_to_use, *) 'b_errors_at_nodes_rho: ', b_errors_at_nodes_rho
  write(unit_to_use, *) 'gamma vector: ', gamma_vec
  write(unit_to_use, *) 'b_nodes_gamma: ', b_nodes_gamma
  write(unit_to_use, *) 'b_errors_at_nodes_gamma: ', b_errors_at_nodes_gamma
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' rho: ' , rho_vec(j) , &
        ' coeff: ' , coeffs_rho(j) , ' nodes: ', nodes_rho(j) , &
        ' error at node: ' , errors_at_nodes_rho(j)
  end do
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' gamma: ' , gamma_vec(j) , &
        ' coeff: ' , coeffs_gamma(j) , ' nodes: ', nodes_gamma(j) , &
        ' error at node: ' , errors_at_nodes_gamma(j)
  end do
end function test_worst_rho_and_gamma_gaussian

function test_worst_rho_and_gamma_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_cauchy)         :: cauchy_kernel
  real(wp)                    :: discrepancy, ker_discrep
  real(wp) , dimension(8)     :: rho_vec, gamma_vec
  real(wp) , dimension(8)     :: coeffs_rho, coeffs_gamma
  real(wp) , dimension(8)     :: nodes_rho, nodes_gamma
  real(wp) , dimension(8)     :: errors_at_nodes_rho, errors_at_nodes_gamma
  real(wp) , dimension(8)     :: b_coeffs_rho, b_coeffs_gamma
  real(wp) , dimension(8)     :: b_nodes_rho, b_nodes_gamma
  real(wp) , dimension(8)     :: b_errors_at_nodes_rho
  real(wp) , dimension(8)     :: b_errors_at_nodes_gamma
  real(wp) , dimension(8,8)   :: a_matrix_rho, a_matrix_gamma
  integer                     :: j , num_iter, ker_num_iter
  integer                     :: err_stat
  character(len=err_msg_len)  :: err_msg
  ! Body
  test_pass = .true.
  call cauchy_kernel%set_c(3E0_wp)
  call cauchy_kernel%set_alpha(1E0_wp)
  call cauchy_kernel%set_x_lb(0E0_wp)
  call cauchy_kernel%set_x_ub(5E0_wp)
  call cauchy_kernel%set_y_lb(-0.5E0_wp)
  call cauchy_kernel%set_y_ub(2E0_wp)
  ! Test the worst-rho-vector functionality
  call worst_rho_and_gamma( &
      kernel_obj = cauchy_kernel, &
      num_terms = size(coeffs_rho), &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      rho_vec = rho_vec, &
      gamma_vec = gamma_vec, &
      coeffs_rho = coeffs_rho, &
      coeffs_gamma = coeffs_gamma, &
      a_matrix_rho = a_matrix_rho, &
      a_matrix_gamma = a_matrix_gamma, &
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
      discrepancy = discrepancy, &
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
  write(unit_to_use, *) 'Worst rho and gamma calculation, ', &
      'Cauchy kernel'
  write(unit_to_use, *) 'x range from ', cauchy_kernel%get_x_lb(), &
      ' to ', cauchy_kernel%get_x_ub()
  write(unit_to_use, *) 'y range from ', cauchy_kernel%get_y_lb(), &
      ' to ', cauchy_kernel%get_y_ub()
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'kernel discrepancy: ', ker_discrep
  write(unit_to_use, *) 'number of iterations: ', num_iter
  write(unit_to_use, *) 'kernel number of iterations: ', ker_num_iter
  write(unit_to_use, *) 'rho vector: ', rho_vec
  write(unit_to_use, *) 'b_nodes_rho: ', b_nodes_rho
  write(unit_to_use, *) 'b_errors_at_nodes_rho: ', b_errors_at_nodes_rho
  write(unit_to_use, *) 'gamma vector: ', gamma_vec
  write(unit_to_use, *) 'b_nodes_gamma: ', b_nodes_gamma
  write(unit_to_use, *) 'b_errors_at_nodes_gamma: ', b_errors_at_nodes_gamma
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' rho: ' , rho_vec(j) , &
        ' coeff: ' , coeffs_rho(j) , ' nodes: ', nodes_rho(j) , &
        ' error at node: ' , errors_at_nodes_rho(j)
  end do
  do j = 1, size(rho_vec)
    write(unit_to_use, *) 'index: ', j , ' gamma: ' , gamma_vec(j) , &
        ' coeff: ' , coeffs_gamma(j) , ' nodes: ', nodes_gamma(j) , &
        ' error at node: ' , errors_at_nodes_gamma(j)
  end do
end function test_worst_rho_and_gamma_cauchy

end module test_worst_rho_and_gamma_mod