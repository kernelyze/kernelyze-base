! eigen_rankn_params.f90
!
! Copyright (c) 2017 by Kernelyze LLC
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
! created on: 2017-04-20 (based on earlier eigen_rankn_kernel)
! updated on: 2017-04-21 (correcting x-y in singular vectors)
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! This module contains a numerical subroutine
! that calculates a rank-$n$ truncated eigenfunction
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
  
module eigen_rankn_params_mod

  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use kernel_mod, only : kernel
  use tk_simple_lapack_mod, only : tk_gesvd ! singular value decomposition

  implicit none

  contains
  
  pure subroutine eigen_rankn_params( &
      kernel_to_approx, &
      rank, &
      x_points, &
      y_points, &
      x_weights, &
      y_weights, &
      x_func_wts, &
      y_func_wts, &
      trunc_sing_vals, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in) :: kernel_to_approx
    integer, intent(in)       :: rank
    real(wp), intent(in)      :: x_points( : )
    real(wp), intent(in)      :: y_points( : )
    real(wp), intent(in)      :: x_weights( size(x_points) )
    real(wp), intent(in)      :: y_weights( size(y_points) )
    ! Note that the Nystrom method constructs the functions of
    ! x as sums over y values and the functions of y as sums
    ! over x values, hence the dimensions of x_func_wts and
    ! y_func_wts below.
    real(wp), intent(out)     :: x_func_wts( size(y_points), rank )
    real(wp), intent(out)     :: y_func_wts( size(x_points), rank )
    real(wp), intent(out)     :: trunc_sing_vals( rank )
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'eigen_rankn_params: '
    ! Local variables
    integer                                                 :: i, j
    integer                                                 :: info
    ! Outputs of the singular value decomposition
    real(wp)  :: singular_vals( min( size(x_points) , size(y_points) ) )
    real(wp)  :: left_singular_vecs( size(x_points) , size(x_points) )
    real(wp)  :: right_singular_vecs_t( size(y_points) , size(y_points) )
    real(wp)  :: wwvec( min( size(x_points) , size(y_points) ) - 1)
    ! Coefficient matrix
    real(wp)  :: coeff_mat( size(x_points) , size(y_points) )
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Check that the rank is <= min( size(x_points), size(y_points) )
    if (rank > size(x_points) .or. rank > size(y_points) ) then
      if (present(err_stat)) then
        err_stat = -1
      end if
      if (present(err_msg)) then
        err_msg = proc_name // 'rank must be <= min size of points'
      end if
      return
    end if
    ! Check that all weights are non-negative:
    if (any(x_weights < 0E0_wp)) then
      if (present(err_stat)) then
        err_stat = -2
      end if
      if (present(err_msg)) then
        err_msg = proc_name // 'all x weights must be >= 0'
      end if
      return
    end if
    if (any(y_weights < 0E0_wp)) then
      if (present(err_stat)) then
        err_stat = -3
      end if
      if (present(err_msg)) then
        err_msg = proc_name // 'all y weights must be >= 0'
      end if
      return
    end if
    ! Discretize the kernel
    do j = 1, size(y_points)
      do i = 1, size(x_points)
        coeff_mat(i , j) = &
            sqrt( x_weights(i) ) * &
            sqrt( y_weights(j) ) * &
            kernel_to_approx%eval( x_points(i) , y_points(j) )
      end do
    end do
    ! Call the interface tk_gesvd to the LAPACK routine dgesvd to
    ! compute the singular value decomposition of the matrix
    ! coeff_mat.
    call tk_gesvd(  a = coeff_mat, &
                    s = singular_vals, &
                    u = left_singular_vecs, &
                    vt = right_singular_vecs_t, &
                    ww = wwvec, &
                    info = info)
    if (info /= 0) then
      ! This means that the execution of the SVD
      ! routine did not go well.
      ! Set the error message and error status if they were provided
      if (present(err_msg)) then
        err_msg = proc_name &
            // 'Problem in computing SVD'
      end if
      if (present(err_stat)) then
        err_stat = info
      end if
      return
    end if
    ! Now use singular values and singular vectors to approximate
    ! the singular functions (eigenfunctions, in the symmetric case)
    ! of the kernel.  First get the square root of the weights out
    ! of the singular vectors:
    do i = 1, size(left_singular_vecs, 1)
      left_singular_vecs(i , :) = &
          left_singular_vecs(i , :) / sqrt( x_weights( i ) )
    end do
    do i = 1, size(right_singular_vecs_t, 2)
      right_singular_vecs_t(: , i) = &
          right_singular_vecs_t(: , i) / sqrt( y_weights( i ) )
    end do
    ! Compute the weights for the functions of x and for the
    ! functions of y that will be part of the truncated
    ! singular function approximation to the kernel.
    do j = 1, rank
      ! Note that the functions of x are sums over y points
      ! and use y weights; this is a feature of the 
      ! Nystrom method.
      do i = 1, size(y_points)
        x_func_wts(i, j) = &
            y_weights(i) * right_singular_vecs_t(j, i) / singular_vals(j)
      end do
    end do
    do j = 1, rank
      ! Note that the functions of y are sums over x points
      ! and use x weights; this is a feature of the
      ! Nystrom method.
      do i = 1, size(x_points)
        y_func_wts(i, j) = &
            x_weights(i) * left_singular_vecs(i, j) / singular_vals(j)
      end do
    end do
    ! Fill the "short" array of singular values
    trunc_sing_vals = singular_vals(1:rank)
  end subroutine eigen_rankn_params

end module eigen_rankn_params_mod
