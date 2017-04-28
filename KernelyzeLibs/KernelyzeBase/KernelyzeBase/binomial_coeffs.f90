! binomial_coeffs.f90
!
! Copyright (c) 2016 by Kernelyze LLC
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
! created on: 2016-02-07
! updated on: 2016-02-07
!
! This module contains functionality to compute
! binomial coefficients.  Why make them reals?
! They are integers, but for very large n they
! become huge and a floating-point approximation
! is an efficient way to handle them.  Also, my
! uses for them typically involve multiplying them
! by reals after generating them, so this avoids
! an explicit type conversion.
  
module binomial_coeffs_mod

use set_precision, only : wp

implicit none

private

public :: binomial_coeffs

contains
  
! This function computes and returns a row of
! Pascal's Triangle; that is, the result returned
! is an array of all (n + 1) binomial coefficients
! of the form (n choose k) for some k >= 0 and <= n.
pure function binomial_coeffs(n) result(coeffs)
  ! Arguments
  integer, intent(in) :: n
  ! Result
  real(wp)            :: coeffs(n + 1)
  ! Local variables
  integer             :: i, j
  real(wp)            :: coeffs_old(n + 1)
  ! Body
  coeffs = 0E0_wp
  coeffs(1) = 1E0_wp
  do i = 2, n + 1
    coeffs_old = coeffs
    do j = 2, i
      coeffs(j) = coeffs_old(j) + coeffs_old(j - 1)
    end do
  end do
end function binomial_coeffs

end module binomial_coeffs_mod