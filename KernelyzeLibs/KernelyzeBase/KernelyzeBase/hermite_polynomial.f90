! hermite_polynomial.f90
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
! created on: 2016-02-08
! updated on: 2016-02-09
! updated on: 2016-02-11 (fixed H_1, was off by factor of 2)
!
! This module provides the Hermite polynomials.
! The Hermite polynomials here are the ``physicist's''
! version.
  
module hermite_polynomial_mod

  use set_precision, only : wp
  
  implicit none

  private

  public :: hermite_polynomial
  public :: hermite_polynomial_vec

contains
  
  pure function hermite_polynomial( n , x ) result(hpoly)
    ! Arguments
    integer, intent(in)   :: n ! n should be >= 0
    real(wp), intent(in)  :: x
    ! Result
    real(wp)              :: hpoly
    ! Local variables
    real(wp)              :: polys_up_to_n( n + 1 )
    ! Body
    polys_up_to_n = hermite_polynomial_vec( n , x )
    ! The last index in polys_up_to_n is the n^th index
    ! but is the (n + 1)^st element, since the Hermite
    ! polynomials are typically indexed starting with 0.
    hpoly = polys_up_to_n( n + 1 )
  end function hermite_polynomial
  
  pure function hermite_polynomial_vec( n , x ) result(hpolys_up_to_n)
    ! Arguments
    integer, intent(in)   :: n ! n should be >= 0
    real(wp), intent(in)  :: x
    ! Result
    real(wp)              :: hpolys_up_to_n( n + 1 )
    ! Local variables
    integer               :: i
    ! Body
    !
    ! Use the recurrence relation for Hermite polynomials.
    ! The result array needs to have (n + 1) elements in 
    ! order to provide indices up to and including n, since
    ! the Hermite polynomials are typically indexed starting
    ! with 0.
    !
    ! The recurrence relation is: 
    ! H_{n+1} \left( x \right) = 
    ! 2 x H_{n} \left( x \right) - 2 n H_{n-1} \left( x \right)
    ! and H_{0} \left( x \right) = 1, H_{1} \left( x \right) = x
    ! are the initial conditions that begin the recursion
    hpolys_up_to_n(1) = 1E0_wp
    if (n > 0) then
      hpolys_up_to_n(2) = 2E0_wp * x
    end if
    do i = 2, n
      hpolys_up_to_n(i + 1) = &
          2E0_wp * x * hpolys_up_to_n(i) &
          - 2E0_wp * real(i - 1, wp) * hpolys_up_to_n(i - 1)
    end do
  end function hermite_polynomial_vec

end module hermite_polynomial_mod