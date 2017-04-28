! moore_penrose_inverse.f90
!
! Copyright (c) 2015, 2017 by Kernelyze LLC
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
! created on: 2015-03-02
! updated on: 2015-04-19
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! A module containing a subroutine that computes
! the Moore-Penrose generalized inverse. The singular
! value decomposition from LAPACK is used to perform
! the computation.
  
module moore_penrose_inverse_mod
  use set_precision, only : wp
  use tk_simple_lapack_mod, only : tk_gesvd
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  pure function moore_penrose_inverse(amat, thresh) result(mpinv)
    ! Arguments
    real(wp) , dimension(: , :) , intent(in)                        :: amat
    real(wp) , intent(in) , optional                                :: thresh
    ! Return value
    real(wp) , dimension(size(amat , 2) , size(amat , 1))           :: mpinv
    ! Local variables
    integer                                                         :: j
    integer                                                         :: info
    integer                                                         :: num_s
    real(wp)                                                        :: toler
    real(wp) , dimension(size(amat , 1) , size(amat , 2))           :: amatcopy
    real(wp) , dimension(min(size(amat , 1) , size(amat , 2)))      :: svec
    real(wp) , dimension(size(amat , 1) , size(amat , 1))           :: umat
    real(wp) , dimension(size(amat , 2) , size(amat , 2))           :: vmat_t
    real(wp) , dimension(min(size(amat , 1) , size(amat , 2)) - 1)  :: wwvec
    ! Body
    ! gesvd below needs an argument that it can overwrite, so
    ! make a copy of amat to preserve the ability to declare it
    ! intent(in) and thus maintain the purity of this function
    amatcopy = amat
    ! Call the interface tk_gesvd to the LAPACK routine dgesvd to
    ! compute the singular value decomposition of the matrix
    ! amat.
    call tk_gesvd( a = amatcopy , s = svec , u = umat , vt = vmat_t , &
        ww = wwvec , info = info)
    if (info /= 0) then
      ! This means that the execution of the SVD
      ! routine did not go well.
      mpinv = ieee_value(mpinv, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Set the threshhold for numerically negligible singular
    ! values of the matrix amat.  If the tolerance was provided,
    ! use it; otherwise, use a reasonable default value.
    if (present(thresh)) then
      toler = thresh
    else
      toler = real(max(size(amat , 1) , size(amat , 2)), wp) &
          * spacing(maxval(abs(svec)))
    end if
    ! Now build a truncated SVD in which all singular
    ! values that are below the threshhold "thresh" are
    ! omitted and the rest are inverted.  This is the
    ! Moore-Penrose generalized inverse.  Recall that
    ! the singular values stored in svec are, by
    ! definition, non-negative.
    num_s = count(svec > toler)
    do j = 1,num_s
      vmat_t(j,:) = vmat_t(j,:) / svec(j)
    end do
    mpinv = transpose(matmul(umat( : , 1:num_s), vmat_t( 1:num_s , : )))
  end function moore_penrose_inverse
end module moore_penrose_inverse_mod