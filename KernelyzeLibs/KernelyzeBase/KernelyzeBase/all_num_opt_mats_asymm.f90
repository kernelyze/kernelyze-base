! all_num_opt_mats_asymm.f90
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
! created on: 2016-07-03
! updated on: 2016-07-03
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! A module containing a subroutine that computes
! the A and B matrices for a given kernel, \rho vector,
! and \gamma vector using compute_a_and_b_matrices and 
! also computes the key matrices V and W in constructing 
! a numerically-optimal approximation (in a sup-norm 
! sense) to the given kernel; unlike the module all_num_opt_mats_mod,
! this module's procedure works for asymmetric kernels.  
! Finally, a vector D is returned which should be a vector 
! of ones; it is returned only to permit output checking 
! for numerical issues.  V, W are matrices of dimension 
! (n + 1) x n where A, B are (n + 1) x (n + 1).  They 
! satisfy
! V * W' = A' * (I - coeffs_rho * coeffs_gamma') * B.
! To see a proof of the claim that d_vec should
! be a vector of ones, in the symmetric case, look to:
! Golub (1973), "Some Modified Eigenvalue Problems,"
! Section 5.  This asymmetric case is not formally covered
! there, but the logic is similar.
!
! Even though the vector D should be a vector of ones, I
! use it in building the V and W matrices because in 
! numerically-challenging cases D may be numerically
! different from a vector of ones, but the desired
! matrices are still both computable and very useful.

module all_num_opt_mats_asymm_mod
  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use mrgrnk_mod, only : mrgrnk
  use compute_a_and_b_matrices_mod, only : compute_a_and_b_matrices
  use kernel_mod, only : kernel
  use tk_simple_lapack_mod, only : tk_gesvd
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  pure subroutine all_num_opt_mats_asymm(& 
      kernel_obj, &
      rho_vec, &
      gamma_vec, &
      coeffs_rho, &
      coeffs_gamma, &
      v_mat, &
      w_mat, &   
      a_mat, &
      b_mat, &
      d_vec, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)                                 :: kernel_obj
    real(wp), dimension(:), intent(in)                        :: rho_vec
    real(wp), dimension(size(rho_vec)), intent(in)            :: gamma_vec
    real(wp), dimension(size(rho_vec)), intent(in)            :: coeffs_rho
    real(wp), dimension(size(rho_vec)), intent(in)            :: coeffs_gamma
    real(wp), dimension(size(rho_vec), size(rho_vec)-1), &
        intent(out)                                           :: v_mat
    real(wp), dimension(size(rho_vec), size(rho_vec)-1), &
        intent(out)                                           :: w_mat
    real(wp), dimension(size(rho_vec), size(rho_vec)), &
        intent(out)                                           :: a_mat
    real(wp), dimension(size(rho_vec), size(rho_vec)), &
        intent(out)                                           :: b_mat
    real(wp), dimension(size(rho_vec)-1), intent(out)         :: d_vec
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                            :: err_stat
    character(len=err_msg_len), intent(out), optional         :: err_msg
    ! Local variables
    integer                                                   :: info, i, j
    integer, dimension(size(rho_vec))                         :: sort_inds
    real(wp), dimension(size(rho_vec))                        :: svals
    real(wp), dimension(size(rho_vec), size(rho_vec))         :: left_svecs
    real(wp), dimension(size(rho_vec), size(rho_vec))         :: right_svecs_t
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Helpful local variables for the numerical work here
    real(wp) :: coeff_mat(size(rho_vec) , size(rho_vec))
    real(wp) :: wwvec(size(rho_vec) - 1)
    ! Local variable to label any problems:
    character(*), parameter         :: proc_name = &
        "Computation of all matrices needed for numerically-optimal " &
        // "approximation: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Obtain the A and B matrices
    call compute_a_and_b_matrices( &
        kernel_obj, &
        rho_vec, &
        gamma_vec, &
        a_mat, &
        b_mat)
    ! Compute  the I - coeffs_rho * coeffs_gamma' matrix
    do j = 1 , size(rho_vec)
      do i = 1 , size(rho_vec)
        coeff_mat(i , j) = - coeffs_rho(i) * coeffs_gamma(j)
      end do
      coeff_mat(j , j) = coeff_mat(j , j) + 1E0_wp
    end do
    ! Call the interface tk_gesvd to the LAPACK SVD routine for 
    ! the matrix coeff_mat.
    call tk_gesvd(  a = coeff_mat, &
                    s = svals, &
                    u = left_svecs, &
                    vt = right_svecs_t, &
                    ww = wwvec, &
                    info = info)
    if (info /= 0) then
      ! This means that the execution of the SVD routine did not go well.
      ! Set the error message and error status if they were provided
      if (present(err_msg)) then
        err_msg = proc_name &
            // 'Problem in computing SVD'
      end if
      if (present(err_stat)) then
        err_stat = info
      end if
      a_mat = ieee_value(a_mat, ieee_quiet_nan)
      b_mat = ieee_value(b_mat, ieee_quiet_nan)
      v_mat = ieee_value(v_mat, ieee_quiet_nan)
      w_mat = ieee_value(w_mat, ieee_quiet_nan)
      d_vec = ieee_value(d_vec, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Now sort the singular values to find the smallest one.
    ! In exact arithmetic, I would get one singular value of
    ! precisely zero and the rest equal to one.
    call mrgrnk(svals, sort_inds, local_err_stat, local_err_msg)
    if (local_err_stat /= 0) then
      ! This means that the mergesort did not go well.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      a_mat = ieee_value(a_mat, ieee_quiet_nan)
      b_mat = ieee_value(b_mat, ieee_quiet_nan)
      v_mat = ieee_value(v_mat, ieee_quiet_nan)
      w_mat = ieee_value(w_mat, ieee_quiet_nan)
      d_vec = ieee_value(d_vec, ieee_quiet_nan)
      ! Return
      return
    end if
    d_vec = svals(sort_inds(2:size(rho_vec)))
    ! Scale by the singular values, even though they should
    ! be trivial (one zero, the rest all equal to one):
    do i = 1, size(rho_vec)
      do j = 1, size(rho_vec)
        ! The max(0 , ...) should not be needed, it is just a safeguard
        left_svecs(   i , j) = left_svecs(   i , j) &
            * sqrt( max( 0E0_wp , svals(j) ) )
        right_svecs_t(i , j) = right_svecs_t(i , j) &
            * sqrt( max( 0E0_wp , svals(i) ) )
      end do
    end do
    ! ***Note*** that at this point the "svecs" are not
    ! really singular vectors any more: in the presence of
    ! numerical noise, the rescaling by sqrt singular values
    ! performed above is helpful, but certainly does not
    ! preserve orthogonality.
    v_mat = matmul( transpose( a_mat ), &
        left_svecs( : , sort_inds( 2 : size(rho_vec) ) ) )
    w_mat = matmul( transpose( b_mat ), &
        transpose( right_svecs_t( sort_inds( 2 : size(rho_vec) ) , : ) ) )
  end subroutine all_num_opt_mats_asymm
end module all_num_opt_mats_asymm_mod