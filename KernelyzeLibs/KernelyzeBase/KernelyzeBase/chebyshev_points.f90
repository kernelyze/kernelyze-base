! chebyshev_points.f90
!
! Copyright (c) 2015 by Kernelyze LLC
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
! updated on: 2015-04-19
!
! A simple utility module containing a function to
! return a set of $n$ Chebyshev points (that is,
! the points at which the $n^{th}$ Chebyshev
! polynomial achieves its maximal absolute value)
! in the closed interval $\left[-1, 1 \right]$.
  
module chebyshev_points_mod
  use set_precision, only : wp
  use constants_mod, only : pi
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  pure function chebyshev_points(n) result(c_points)
    integer, intent(in)     :: n
    real(wp), dimension(n)  :: c_points
    real(wp), dimension(n)  :: grid
    integer                 :: j
    if (n < 2) then
      ! Set the results to NaN
      c_points = ieee_value(c_points, ieee_quiet_nan)
      ! Return
      return
    end if
    do j = 1,n
      grid(j) = -(n-1) + 2.0E0_wp * (j - 1)
    end do
    grid = pi * grid / ( 2.0E0_wp * (n - 1) )
    c_points = sin(grid)
  end function chebyshev_points
end module chebyshev_points_mod