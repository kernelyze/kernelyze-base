! eigen_rankn_kernel.f90
!
! Copyright (c) 2016, 2017 by Kernelyze LLC
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
! created on: 2016-12-04
! updated on: 2017-01-01 (work toward initial compilation)
! updated on: 2017-01-02 (debugging)
! updated on: 2017-04-20 (refactored to make the parameter
!             calculation separate from the approximate-
!             kernel allocation and build)
!
! This module contains a constructor subroutine
! that builds a rank-$n$ truncated eigenfunction
! series approximation, given a kernel to 
! approximate and a permitted approximation rank, $n$.
! Since asymmetric kernels are permitted, this is
! perhaps more precisely termed a singular function
! approximation (a continuous analog to the discrete
! singular value decomposition for a matrix). The
! technique used is Nystrom's method, familar from
! the integral equation literature.  The subtlety
! here, which is alluded to for the symmetric case
! with eigenfunctions and eigenvalues in Delves & Mohamed (1985)
! (Computational Methods for Integral Equations) on
! page 152ff., is that the weights used to approximate
! the integrals involved (via quadrature) must be 
! accounted for.  Here I require non-negative weights,
! which allows the SVD to be applied to
! $ M \equiv \sqrt{W_x} K \sqrt{W_y} $, where $W_x$ is the
! diagonal matrix with diagonal elements $ w_x ( x_i )$
! and $W_y$ is the diagonal matrix with diagonal
! elements $ w_y ( y_j ) $ (so that the square roots
! are unambiguous), while $K = K \left( x_i , y_j \right)$
! is the discretization of the kernel.  The resulting 
! singular vectors have orthogonality that is to be 
! interpreted with respect to the weight vectors
! given.  To recover an approximation to the kernel,
! I first transform the left and right singular
! vectors $ U $ and $ V $ in the SVD of $ M $ as
! defined above: $ M = U D V^{T} $:
! $ \tilde{ U } = \sqrt{ W_x^{-1} } U $ and
! $ \tilde{ V } = \sqrt{ W_y^{-1} } V $.
! Then I follow the Nystrom method and form approximate
! integrals with respect to the kernel:
! $ \tilde{ g_j } \left( y \right) =
! \frac{ 1 }{ d_j } \sum_{i = 1}^{n} K \left( x_i , y \right) 
! w_x ( x_i ) \tilde{ u_{i j} }$.
  
module eigen_rankn_kernel_mod

  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use kernel_mod, only : kernel
  use kernel_neigen_mod, only : kernel_rankn
  use eigen_rankn_params_mod, only : eigen_rankn_params
  use eigen_rankn_kernel_of_params_mod, only : eigen_rankn_kernel_of_params

  implicit none

  contains
  
  subroutine eigen_rankn_kernel( &
      kernel_to_approx, &
      rank, &
      x_points, &
      y_points, &
      x_weights, &
      y_weights, &
      rankn_kernel, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    integer, intent(in)           :: rank
    real(wp), intent(in)          :: x_points( : )
    real(wp), intent(in)          :: y_points( : )
    real(wp), intent(in)          :: x_weights( size(x_points) )
    real(wp), intent(in)          :: y_weights( size(y_points) )
    class(kernel_rankn), allocatable, intent(out) :: rankn_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'eigen_rankn_kernel: '
    ! Intermediate outputs -- the parameters of the
    ! rank-$n$ truncated singular function series
    ! approximation
    real(wp)  :: x_func_wts( size(y_points), rank )
    real(wp)  :: y_func_wts( size(x_points), rank )
    real(wp)  :: trunc_sing_vals( rank )
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Obtain the rank-$n$ parameters
    call eigen_rankn_params( &
        kernel_to_approx = kernel_to_approx, &
        rank = rank, &
        x_points = x_points, &
        y_points = y_points, &
        x_weights = x_weights, &
        y_weights = y_weights, &
        x_func_wts = x_func_wts, &
        y_func_wts = y_func_wts, &
        trunc_sing_vals = trunc_sing_vals, &
        err_stat = local_err_stat, &
        err_msg = local_err_msg)
    ! Handle any error
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    ! Build the approximation using the rank-$n$ parameters
    call eigen_rankn_kernel_of_params( &
        kernel_to_approx = kernel_to_approx, &
        x_points = x_points, &
        y_points = y_points, &
        trunc_sing_vals = trunc_sing_vals, &
        x_func_wts = x_func_wts, &
        y_func_wts = y_func_wts, &
        rankn_kernel = rankn_kernel, &
        err_stat = local_err_stat, &
        err_msg = local_err_msg)
    ! Handle any error
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
  end subroutine eigen_rankn_kernel

end module eigen_rankn_kernel_mod