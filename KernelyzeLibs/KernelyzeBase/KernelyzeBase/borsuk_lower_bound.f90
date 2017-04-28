! borsuk_lower_bound.f90
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
! created on: 2015-02-07
! updated on: 2015-04-22
! updated on: 2016-04-24 (Avoid upstream use of internal procedures 
!             for thread safety)
! updated on: 2016-07-01 (Accommodate rectangular domains and
!             asymmetry)
!
! A module containing the subroutine borsuk_lower_bound.
! Given a set of points in rho_vec, this subroutine
! minimizes the maximum deviation from zero 
! (if over x)
! $max_{\left[x_L, x_U \right]} \left|coeffs' * vec_of_kernel(x) \right|$ or
! (if over y)
! $max_{\left[y_L, y_U \right]} \left|coeffs' * vec_of_kernel(y) \right|$
! The term 'vec_of_kernel' is a vector of kernel evaluations, 
! $K \left( x, \rho_j \right)$ if over x or 
! $K \left( \gamma_j ,  y \right)$ if over y.
! The 'coeffs' term is a vector of 
! coefficients to be solved for (with the constraint that they must
! sum to one in absolute value) to minimize the maximum deviation
! from zero.
!
! NOTE that the kernel $K \left( x , y \right) $ must be
! nondegenerate totally positive on the square
! $\left[x_L , x_U \right] \times \left[ y_L , y_U \right]$
! to use the subroutine borsuk_lower_bound.
  
module borsuk_lower_bound_mod
  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use remez_variant_mod, only : remez_variant
  use kernel_mod, only : kernel
  implicit none
  contains
  subroutine borsuk_lower_bound( &
      kernel_obj, &
      rho_vec, &
      is_over_x, &
      tolerance, &
      max_iter, &
      grid, &
      coeffs, &
      nodes, &
      errors_at_nodes, &
      discrepancy, &
      num_iter, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)                         :: kernel_obj
    real(wp), intent(in), dimension(:)                :: rho_vec
    logical, intent(in)                               :: is_over_x
    real(wp), intent(in)                              :: tolerance
    integer, intent(in)                               :: max_iter
    real(wp), optional, intent(in), dimension(:)      :: grid
    real(wp), intent(out), dimension(size(rho_vec))   :: coeffs
    real(wp), intent(out), dimension(size(rho_vec))   :: nodes
    real(wp), intent(out), dimension(size(rho_vec))   :: errors_at_nodes
    real(wp), intent(out)                             :: discrepancy
    integer, intent(out)                              :: num_iter
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Continue appropriately
    if (present(err_stat)) then
      if (present(err_msg)) then
        call remez_variant( &
            is_borsuk = .true., &
            kernel_obj = kernel_obj, &
            is_over_x = is_over_x, &
            rho_vec = rho_vec, &
            tolerance = tolerance, &
            max_iter = max_iter, &
            grid = grid, &
            coeffs = coeffs, &
            nodes = nodes, &
            errors_at_nodes = errors_at_nodes, &
            discrepancy = discrepancy, &
            num_iter = num_iter, &
            err_stat = err_stat, &
            err_msg = err_msg)
      else
        call remez_variant( &
            is_borsuk = .true., &
            kernel_obj = kernel_obj, &
            is_over_x = is_over_x, &
            rho_vec = rho_vec, &
            tolerance = tolerance, &
            max_iter = max_iter, &
            grid = grid, &
            coeffs = coeffs, &
            nodes = nodes, &
            errors_at_nodes = errors_at_nodes, &
            discrepancy = discrepancy, &
            num_iter = num_iter, &
            err_stat = err_stat)
      end if
    else
      if (present(err_msg)) then
        call remez_variant( &
            is_borsuk = .true., &
            kernel_obj = kernel_obj, &
            is_over_x = is_over_x, &
            rho_vec = rho_vec, &
            tolerance = tolerance, &
            max_iter = max_iter, &
            grid = grid, &
            coeffs = coeffs, &
            nodes = nodes, &
            errors_at_nodes = errors_at_nodes, &
            discrepancy = discrepancy, &
            num_iter = num_iter, &
            err_msg = err_msg)
      else
        call remez_variant( &
            is_borsuk = .true., &
            kernel_obj = kernel_obj, &
            is_over_x = is_over_x, &
            rho_vec = rho_vec, &
            tolerance = tolerance, &
            max_iter = max_iter, &
            grid = grid, &
            coeffs = coeffs, &
            nodes = nodes, &
            errors_at_nodes = errors_at_nodes, &
            discrepancy = discrepancy, &
            num_iter = num_iter)
      end if
    end if
  end subroutine borsuk_lower_bound
end module borsuk_lower_bound_mod