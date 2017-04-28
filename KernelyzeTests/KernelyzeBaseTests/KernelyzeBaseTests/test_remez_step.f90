! test_remez_step.f90
!
! Copyright (c) 2015, 2016, 2017 by Kernelyze LLC
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
! created on: 2015-02-25
! updated on: 2015-04-22
! updated on: 2016-04-25 (to use kernel objects)
! updated on: 2016-10-03 (to account for rescaling of kernel_gaussian)
! updated on: 2017-04-22 (switch to use of tk_simple_lapack module)
!
! A module of unit tests for the remez_step module.
  
module test_remez_step_mod

use set_precision, only : wp
use constants_mod, only : inv_sqrt2pi
use remez_step_mod, only  : remez_step
use kernel_test_funcs_mod, only : gaussian_kernel, exp_prod_kernel
use chebyshev_points_mod, only : chebyshev_points
use tk_simple_lapack_mod, only : tk_gesv
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none

contains
  
function test_remez_step(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Test of unirnk
  write(unit_to_use , *) 'Test of Remez-like step functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_remez_step_simple(unit_to_use)
  test_pass = test_pass .and. test_remez_step_many_pts(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_remez_step

function test_remez_step_simple(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp), dimension(1)    :: test_rho_vec
  real(wp), dimension(1)    :: test_coeffs
  real(wp), dimension(200)  :: test_grid
  real(wp), dimension(1)    :: test_nodes
  real(wp), dimension(1)    :: test_vals_at_nodes
  real(wp)                  :: test_toler
  ! Test of (very simple) Remez-like step
  test_pass = .true.
  
  test_rho_vec(1) = 0E0_wp
  test_coeffs(1) = 1.0E0_wp
  test_grid = chebyshev_points(200)
  test_toler = 1.0E-15_wp
  
  write(unit_to_use , *) 'Test of Remez-like step in very simple case: '
  ! Note that the optional argument "a_matrix_in" is *not*
  ! supplied to remez_step below, so a_matrix within
  ! remez_step is set equal to the identity matrix.
  call remez_step(kernel_obj = gaussian_kernel, &
      rho_vec = test_rho_vec, &
      coeffs = test_coeffs, &
      tolerance = test_toler, &
      grid = test_grid, &
      nodes = test_nodes, &
      values_at_nodes = test_vals_at_nodes, &
      is_over_x = .true.)
  write(unit_to_use , *) 'Node should be zero and is: ' , test_nodes(1)
  write(unit_to_use , *) 'Value at node should be ', - inv_sqrt2pi, &
      ' and is: ', test_vals_at_nodes(1)
  test_pass = test_pass &
      .and. (abs(test_nodes(1)) < sqrt(test_toler)) &
      .and. (abs(test_vals_at_nodes(1) + inv_sqrt2pi) < test_toler)
end function test_remez_step_simple

function test_remez_step_many_pts(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Result
  logical                     :: test_pass
  ! Local variables
  real(wp), dimension(20)     :: test_rho_vec
  real(wp), dimension(20)     :: target_vec
  real(wp), dimension(20)     :: test_coeffs
  real(wp), dimension(200)    :: test_grid
  real(wp), dimension(20)     :: test_nodes
  real(wp), dimension(20)     :: test_vals_at_nodes
  real(wp), dimension(20,20)  :: func_mat, a_mat
  integer, dimension(20)      :: ipiv
  real(wp)                    :: test_toler
  integer                     :: i, j, lapack_info
  ! Test of (very simple) Remez-like step
  test_pass = .true.
  
  test_rho_vec = chebyshev_points(20)
  test_grid = chebyshev_points(200)
  test_toler = 1.0E-15_wp
  
  ! A decent initial guess as to the nodes where the max abs
  ! error is attained is the Chebyshev points
  test_nodes = chebyshev_points(size(test_rho_vec))
  do j = 1 , size(test_nodes)
    do i = 1 , size(test_rho_vec)
      func_mat(i , j) = exp_prod_kernel%eval(test_rho_vec(i) , test_nodes(j))
    end do
  end do
  ! Alternating array: -1, 1, -1, 1, etc.
  target_vec = (/ ( (-1.0E0_wp)**j, j = 1, size(test_rho_vec) ) /)
  test_coeffs = target_vec ! This will be overwritten momentarily
  a_mat = func_mat ! This will be overwritten momentarily
  ! The subroutine tk_gesv overwrites a_mat using the L and U
  ! factors in the LU decomposition; it also overwrites
  ! coeffs with the solution of the linear system
  call tk_gesv( a_mat , test_coeffs , ipiv , lapack_info)
  test_pass = test_pass .and. (lapack_info == 0)
  ! Normalize the coeffs so that their abs values sum to
  ! one.
  test_coeffs = test_coeffs / sum(abs(test_coeffs))
  
  write(unit_to_use , *) 'Test of Remez-like step with many points: '
  ! Note that the optional argument "a_matrix_in" is *not*
  ! supplied to remez_step below, so a_matrix within
  ! remez_step is set equal to the identity matrix.
  call remez_step(kernel_obj = exp_prod_kernel, &
      rho_vec = test_rho_vec, &
      coeffs = test_coeffs, &
      tolerance = test_toler, &
      grid = test_grid, &
      nodes = test_nodes, &
      values_at_nodes = test_vals_at_nodes, &
      is_over_x = .true.)
  write(unit_to_use , *) 'Coeffs: '
  do j = 1,size(test_coeffs)
    write(unit_to_use , *) test_coeffs(j)
  end do
  write(unit_to_use , *) 'Nodes: '
  do j = 1,size(test_nodes)
    write(unit_to_use , *) test_nodes(j)
  end do
  write(unit_to_use , *) 'Values at nodes: '
  do j = 1,size(test_vals_at_nodes)
    write(unit_to_use , *) test_vals_at_nodes(j)
  end do
  ! I expect the difference of this function from zero to be
  ! below the given tolerance
  test_pass = test_pass &
      .and. (maxval(abs(test_vals_at_nodes)) < test_toler)
end function test_remez_step_many_pts

end module test_remez_step_mod