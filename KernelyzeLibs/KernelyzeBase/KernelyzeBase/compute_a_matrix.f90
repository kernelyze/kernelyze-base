! compute_a_matrix.f90
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
! created on: 2015-03-07
! updated on: 2015-04-19
! updated on: 2016-04-24 (Avoid upstream use of internal procedures 
!             for thread safety)
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! A module containing a subroutine that computes
! the A matrix for a given kernel and \rho vector
! using the symmetric eigendecomposition provided
! by LAPACK.

module compute_a_matrix_mod
  use set_precision, only : wp
  use kernel_mod, only : kernel
  ! Symmetric eigendecomposition, packed storage
  use tk_simple_lapack_mod, only : tk_spev
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  pure function compute_a_matrix(kernel_obj, rho_vec) result(a_matrix)
    ! Arguments
    class(kernel), intent(in)                                 :: kernel_obj
    real(wp) , dimension(:), intent(in)                       :: rho_vec
    ! Return value
    real(wp) , dimension(size(rho_vec) , size(rho_vec))       :: a_matrix
    ! Local variables
    integer                                                   :: info, i, j
    real(wp) , dimension(size(rho_vec))                       :: eigvals
    real(wp) , dimension(size(rho_vec) , size(rho_vec))       :: eigvecs
    ! The dimension of kermat arises as a result of the packed
    ! storage of the upper-triangular part of the full matrix.
    real(wp) , &
        dimension((size(rho_vec) * (size(rho_vec) + 1)) / 2)  :: kermat
    ! Local variable to label any problems:
    character(*), parameter         :: proc_name = &
        "Computation of an A matrix: "
    ! Body
    ! Fill in the packed-form version of the upper triangle of
    ! the matrix kernel( rho_vec * rho_vec'):
    do j = 1 , size(rho_vec)
      do i = 1 , j
        ! Observe that j * (j - 1) is always even, so integer division
        ! of that quantity by 2 does not inadvertently round down.
        kermat(i + ((j * (j -1)) / 2)) = kernel_obj%eval(rho_vec(i) , rho_vec(j))
      end do
    end do
    ! Call the interface tk_spev to the LAPACK routine dspev to
    ! compute the symmetric eigendecomposition of the packed-storage
    ! matrix kermat.
    call tk_spev(  ap = kermat , &
                w = eigvals , &
                uplo = 'U' , &  ! The upper triangle is packed
                z = eigvecs, &
                info = info)
    if (info /= 0) then
      ! This means that the execution of the eigendecomposition
      ! routine did not go well.
      a_matrix = ieee_value(a_matrix, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Replacing any eigenvalues that are less than epsilon
    ! with epsilon does not interfere with the property that
    ! A * kernel(rho_vec * rho_vec') * A' is approximately
    ! the identity.
    where (eigvals < epsilon(1E0_wp)) eigvals = epsilon(1E0_wp)
    ! Make A = diag(1/sqrt(eigvals)) * eigvecs'
    do j = 1 , size(rho_vec)
      do i = 1 , size(rho_vec)
        ! a_matrix(i, j) = eigvecs(j , i) / sqrt(eigvals(i))
        a_matrix(i, j) = eigvecs(j , i) * sqrt(1E0_wp / eigvals(i))
      end do
    end do
  end function compute_a_matrix
end module compute_a_matrix_mod