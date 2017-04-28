! compute_a_and_b_matrices.f90
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
! created on: 2016-07-02
! updated on: 2016-07-02
! updated on: 2016-07-03 (added some temporary checks, then deleted them;
!             this helped to find a silly loop-bound error, which I fixed)
! updated on: 2016-09-10 (avoid an auto-parallelization that causes
!             performance issues due to nested parallelization)
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! A module containing a subroutine that computes
! the A and B matrices for a given kernel, $\rho$ vector,
! and $\gamma$ vector using the SVD provided by LAPACK.
! The A and B matrices have the essential property that
! B * kermat * A' is approximately the identity, where
! kermat(i, j) = kernel($\gamma_i$, $\rho_j$).

module compute_a_and_b_matrices_mod
  use set_precision, only : wp
  use kernel_mod, only : kernel
  use tk_simple_lapack_mod, only : tk_gesvd
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  pure subroutine compute_a_and_b_matrices( &
      kernel_obj, &
      rho_vec, &
      gamma_vec, &
      a_matrix, &
      b_matrix)
    ! Arguments
    class(kernel), intent(in) :: kernel_obj
    real(wp), intent(in)      :: rho_vec(:)
    real(wp), intent(in)      :: gamma_vec(size(rho_vec))
    ! Return value
    real(wp), intent(out)     :: a_matrix(size(rho_vec), size(rho_vec))
    real(wp), intent(out)     :: b_matrix(size(rho_vec), size(rho_vec))
    ! Local variables
    integer   :: info, i, j
    real(wp)  :: singular_vals(size(rho_vec))
    real(wp)  :: left_singular_vecs(size(rho_vec), size(rho_vec))
    real(wp)  :: right_singular_vecs_t(size(rho_vec), size(rho_vec))
    real(wp)  :: wwvec(size(rho_vec) - 1)
    real(wp)  :: kermat(size(rho_vec) , size(rho_vec))
    ! Local variable to label any problems:
    character(*), parameter :: proc_name = &
        "Computation of an A matrix and a B matrix: "
    ! Body
    ! Fill in the matrix k(i,j) = kernel( gamma_vec(i), rho_vec(j) ):
    do j = 1 , size(rho_vec)
      do i = 1 , size(rho_vec)
        kermat(i , j) = kernel_obj%eval(gamma_vec(i) , rho_vec(j))
      end do
    end do
    ! Call the tk_simple_lapack routine tk_gesvd to compute the SVD of kermat.
    call tk_gesvd(  a = kermat, &
                    s = singular_vals, &
                    u = left_singular_vecs, &
                    vt = right_singular_vecs_t, &
                    ww = wwvec, &
                    info = info)
    if (info /= 0) then
      ! This means that the execution of the SVD routine did not
      ! go well.
      a_matrix = ieee_value(a_matrix, ieee_quiet_nan)
      b_matrix = ieee_value(b_matrix, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Replacing any singular values that are less than epsilon
    ! with epsilon does not interfere with the property that
    ! B * kermat * A' is approximately the identity.
    where (singular_vals < epsilon(1E0_wp)) singular_vals = epsilon(1E0_wp)
    ! Make A = diag(1/sqrt(singular_vals)) * right_singular_vecs'
    !      B = diag(1/sqrt(singular_vals)) * left_singular_vals'
    ! Because this function may be used heavily within threads during
    ! multithreaded execution, and because this do loop would otherwise
    ! be auto-parallelized, be sure that it is *not* auto-parallelized
    ! (experience and VTune indicate that allowing auto-parallelization
    ! can cause a substantive performance problem due to
    ! __kmp_fork_barrier).
    !$DIR NOPARALLEL
    do j = 1 , size(rho_vec)
      do i = 1 , size(rho_vec)
        ! The right_singular_vec_t is already the transpose of the matrix of
        ! right singular vectors (it is V^{T} in the U * S * V^{T} expression)
        a_matrix(i, j) = right_singular_vecs_t(i , j) &
            * sqrt(1E0_wp / singular_vals(i))
        b_matrix(i, j) = left_singular_vecs(j , i) &
            * sqrt(1E0_wp / singular_vals(i))
      end do
    end do
  end subroutine compute_a_and_b_matrices
end module compute_a_and_b_matrices_mod