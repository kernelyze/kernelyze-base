! test_borsuk_lower_bound.f90
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
! created on: 2015-02-15
! updated on: 2015-04-22
! updated on: 2016-07-01 (Accommodate interface changes
!             made to allow for asymmetry and general
!             rectangular domains)
!
! A module of unit tests for the borsuk_lower_bound module.
  
module test_borsuk_lower_bound_mod

use set_precision, only : wp
use chebyshev_points_mod, only : chebyshev_points
use kernel_test_funcs_mod, only : &
    exp_prod_kernel, gaussian_kernel, cauchy_kernel
use borsuk_lower_bound_mod, only : borsuk_lower_bound

implicit none

contains
  
function test_borsuk_lower_bound(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of Borsuk-lower-bound functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_borsuk_lower_bound_kernel_one(unit_to_use)
  test_pass = test_pass .and. test_borsuk_lower_bound_kernel_two(unit_to_use)
  test_pass = test_pass .and. test_borsuk_lower_bound_kernel_three(unit_to_use)
  test_pass = test_pass .and. test_borsuk_lower_bound_duplicates(unit_to_use)
  test_pass = test_pass .and. test_borsuk_lower_bound_tiny(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_borsuk_lower_bound

function test_borsuk_lower_bound_kernel_one(unit_to_use) result(test_pass)

  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  ! Local variables
  real(wp)                                  :: discrepancy
  real(wp) , dimension(40)                  :: cheb_pts
  real(wp) , dimension(5)                   :: rho_vec, coeffs, nodes, errors_at_nodes
  integer                                   :: j , num_iter
  ! Body
  test_pass = .true.
  ! Get a convenient set of Chebyshev points
  cheb_pts = chebyshev_points(40)
  ! Test the Borsuk lower bound functionality
  rho_vec = chebyshev_points(5)
  call borsuk_lower_bound( &
      kernel_obj = exp_prod_kernel, &
      rho_vec = rho_vec, &
      is_over_x = .true., &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      grid = cheb_pts, &
      coeffs = coeffs, &
      nodes = nodes, &
      errors_at_nodes = errors_at_nodes, &
      discrepancy = discrepancy, &
      num_iter = num_iter)
  write(unit_to_use, *) 'Borsuk lower bound calculation'
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'number of iterations: ', num_iter
  do j = 1,5
    write(unit_to_use, *) 'index: ', j , ' coeff: ' , coeffs(j) , ' nodes: ', nodes(j) , ' error at node: ' , errors_at_nodes(j)
  end do
end function test_borsuk_lower_bound_kernel_one

function test_borsuk_lower_bound_kernel_two(unit_to_use) result(test_pass)

  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  ! Local variables
  real(wp)                                  :: discrepancy
  real(wp) , dimension(40)                  :: cheb_pts
  real(wp) , dimension(5)                   :: rho_vec, coeffs, nodes, errors_at_nodes
  integer                                   :: j , num_iter
  ! Body
  test_pass = .true.
  ! Get a convenient set of Chebyshev points
  cheb_pts = chebyshev_points(40)
  ! Test the Borsuk lower bound functionality on a different kernel
  rho_vec = chebyshev_points(5)
  call borsuk_lower_bound( &
      kernel_obj = gaussian_kernel, &
      rho_vec = rho_vec, &
      is_over_x = .true., &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      grid = cheb_pts, &
      coeffs = coeffs, &
      nodes = nodes, &
      errors_at_nodes = errors_at_nodes, &
      discrepancy = discrepancy, &
      num_iter = num_iter)
  write(unit_to_use, *) 'Borsuk lower bound calculation for Gaussian kernel'
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'number of iterations: ', num_iter
  do j = 1,5
    write(unit_to_use, *) 'index: ', j , ' coeff: ' , coeffs(j) , ' nodes: ', nodes(j) , ' error at node: ' , errors_at_nodes(j)
  end do
end function test_borsuk_lower_bound_kernel_two

function test_borsuk_lower_bound_kernel_three(unit_to_use) result(test_pass)

  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  ! Local variables
  real(wp)                                  :: discrepancy
  real(wp) , dimension(40)                  :: cheb_pts
  real(wp) , dimension(5)                   :: rho_vec, coeffs, nodes, errors_at_nodes
  integer                                   :: j , num_iter
  ! Body
  test_pass = .true.
  ! Get a convenient set of Chebyshev points
  cheb_pts = chebyshev_points(40)
  ! Test the Borsuk lower bound functionality on a different kernel
  rho_vec = chebyshev_points(5)
  call borsuk_lower_bound( &
      kernel_obj = cauchy_kernel, &
      rho_vec = rho_vec, &
      is_over_x = .true., &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      grid = cheb_pts, &
      coeffs = coeffs, &
      nodes = nodes, &
      errors_at_nodes = errors_at_nodes, &
      discrepancy = discrepancy, &
      num_iter = num_iter)
  write(unit_to_use, *) 'Borsuk lower bound calculation for Cauchy kernel'
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'number of iterations: ', num_iter
  do j = 1,5
    write(unit_to_use, *) 'index: ', j , ' coeff: ' , coeffs(j) , ' nodes: ', nodes(j) , ' error at node: ' , errors_at_nodes(j)
  end do
end function test_borsuk_lower_bound_kernel_three

function test_borsuk_lower_bound_duplicates(unit_to_use) result(test_pass)

  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  ! Local variables
  real(wp)                                  :: discrepancy
  real(wp) , dimension(40)                  :: cheb_pts
  real(wp) , dimension(5)                   :: rho_vec, coeffs, nodes, errors_at_nodes
  integer                                   :: j , num_iter
  ! Body
  test_pass = .true.
  ! Get a convenient set of Chebyshev points
  cheb_pts = chebyshev_points(40)
  rho_vec = chebyshev_points(5)
  ! Test the Borsuk lower bound functionality when there is a duplicate in the rho vector
  rho_vec(1) = -1.0E0_wp ! This should be redundant
  rho_vec(2) = -1.0E0_wp
  call borsuk_lower_bound( &
      kernel_obj = exp_prod_kernel, &
      rho_vec = rho_vec, &
      is_over_x = .true., &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      grid = cheb_pts, &
      coeffs = coeffs, &
      nodes = nodes, &
      errors_at_nodes = errors_at_nodes, &
      discrepancy = discrepancy, &
      num_iter = num_iter)
  write(unit_to_use, *) 'Borsuk lower bound calculation with duplicate entries in rho_vec(1) and rho_vec(2)'
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'number of iterations: ', num_iter
  do j = 1,5
    write(unit_to_use, *) 'index: ', j , ' coeff: ' , coeffs(j) , ' nodes: ', nodes(j) , ' error at node: ' , errors_at_nodes(j)
  end do
end function test_borsuk_lower_bound_duplicates

function test_borsuk_lower_bound_tiny(unit_to_use) result(test_pass)

  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  ! Local variables
  real(wp)                                  :: discrepancy
  real(wp) , dimension(15)                  :: rho_vec_large, coeffs_large, nodes_large, errors_at_nodes_large
  integer                                   :: j , num_iter
  ! Body
  test_pass = .true.
  ! Now test the Borsuk lower bound for a problem with error close to rounding error
  rho_vec_large = chebyshev_points(15)
  call borsuk_lower_bound( &
      kernel_obj = exp_prod_kernel, &
      rho_vec = rho_vec_large, &
      is_over_x = .true., &
      tolerance = 1.0E-15_wp, &
      max_iter = 100, &
      grid = chebyshev_points(200), &
      coeffs = coeffs_large, &
      nodes = nodes_large, &
      errors_at_nodes = errors_at_nodes_large, &
      discrepancy = discrepancy, &
      num_iter = num_iter)
  write(unit_to_use, *) 'Borsuk lower bound calculation close to rounding error'
  write(unit_to_use, *) 'discrepancy: ', discrepancy
  write(unit_to_use, *) 'number of iterations: ', num_iter
  do j = 1,15
    write(unit_to_use, *) 'index: ', j , ' coeff: ' , coeffs_large(j) , &
        ' nodes: ', nodes_large(j) , ' error at node: ' , errors_at_nodes_large(j)
  end do
end function test_borsuk_lower_bound_tiny

end module test_borsuk_lower_bound_mod