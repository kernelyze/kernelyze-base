! taylor_series_rankn_params.f90
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
! created on: 2017-04-20 (based on earlier taylor_series_rankn_kernel)
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! This module contains a numerical subroutine that builds
! the left and right coefficient matrices of a rank-$n$
! two-variable Taylor series approximating the given kernel.
! The Taylor series can be expressed as:
! $ \left( 1 x x^2 \cdots x^{n-1} \right) V 
! W^T \left( 1 y y^2 \cdots y^{n-1} \right)^T$,
! and the V and W matrices are computed here.
! Why not just compute the overall coefficient
! matrix $V W^T$?  Because this approach fits 
! downstream interfaces and because it implicitly
! makes a (reasonable) choice about the form
! of the functions $ f_i \left( x \right) $ and
! $ g_i \left( y \right) $ in the equivalent expression
! $ \sum_i^n f_i \left( x \right) g_i \left( y \right) $
! for the approximation to the kernel.
  
module taylor_series_rankn_params_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use factorial_mod, only : factorial
  use kernel_mod, only : kernel
  use tk_simple_lapack_mod, only : tk_gesvd ! singular value decomposition

  implicit none

  contains
  
  pure subroutine taylor_series_rankn_params( &
      kernel_to_approx, &
      x_center, &
      y_center, &
      rank, &
      v_matrix, &
      w_matrix, &
      err_stat, &
      err_msg, &
      fin_diff_delta)
    ! Arguments
    class(kernel), intent(in) :: kernel_to_approx    
    real(wp), intent(in)      :: x_center
    real(wp), intent(in)      :: y_center
    integer, intent(in)       :: rank
    real(wp), intent(out)     :: v_matrix(rank, rank)
    real(wp), intent(out)     :: w_matrix(rank, rank)
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Optional finite-diff parameter
    real(wp), intent(in), optional  :: fin_diff_delta
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'taylor_series_rankn_params: '
    integer                           :: i, j
    integer                           :: info
    ! Outputs of the singular value decomposition
    real(wp)  :: singular_vals(rank)
    real(wp)  :: left_singular_vecs(rank, rank)
    real(wp)  :: right_singular_vecs_t(rank, rank)
    real(wp)  :: wwvec(rank - 1)
    ! Coefficient matrix
    real(wp)  :: coeff_mat(rank, rank)
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Fill the coefficient matrix with appropriate partial
    ! derivatives of the kernel:
    do j = 1, rank
      do i = 1, rank
        if (present(fin_diff_delta)) then
          coeff_mat(i , j) = kernel_to_approx%mth_nth_partial( &
              i - 1 , j - 1 , x_center , y_center , fin_diff_delta )
        else
          coeff_mat(i , j) = kernel_to_approx%mth_nth_partial( &
              i - 1 , j - 1 , x_center , y_center )
        end if
      end do
    end do
    ! Now scale the (i, j) element by (1 / ((i-1)! * (j-1)!)) to get
    ! the Taylor series coefficient matrix:
    do j = 1, rank
      do i = 1, rank
        coeff_mat(i , j) = coeff_mat(i , j) / &
            ( factorial(i - 1) * factorial(j - 1) )
      end do
    end do
    ! Call the interface gesvd to the LAPACK routines ?gesvd to
    ! compute the singular value decomposition of the matrix
    ! coeff_mat.
    call tk_gesvd( a = coeff_mat, &
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
    ! Use the singular values and singular vectors to build
    ! the V and W matrices
    do i = 1, rank
      do j = 1, rank
        v_matrix(i , j) = left_singular_vecs(i , j) &
            * sqrt( singular_vals(j) )
        ! The right_singular_vec_t is the transpose of the matrix of right
        ! singular vectors (it is V^{T} in the U * S * V^{T} expression)
        w_matrix(i , j) = right_singular_vecs_t(j , i) &
            * sqrt( singular_vals(j) )
      end do
    end do
  end subroutine taylor_series_rankn_params

end module taylor_series_rankn_params_mod