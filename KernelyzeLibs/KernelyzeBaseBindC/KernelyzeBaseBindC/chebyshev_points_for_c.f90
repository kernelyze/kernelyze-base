! chebyshev_points_for_c.f90
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
! created on: 2017-04-14
!
! A module that gives C types and a binding for
! Fortran functionality that gets a set of
! Chebyshev points (that is, the points
! at which the $n^{th}$ Chebyshev polynomial
! achieves its maximal absolute value) in the
! closed interval $\left[-1, 1 \right]$.

module chebyshev_points_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use chebyshev_points_mod, only : chebyshev_points
use, intrinsic :: iso_c_binding

implicit none

contains
  
pure subroutine chebyshev_points_for_c(n, pts) &
    bind(C, name='chebyshev_points')
  ! Arguments
  integer(c_int), intent(in)        :: n
  real(c_double), intent(out)       :: pts(n)
  ! Body
  ! Call the Fortran chebyshev_points function
  pts = chebyshev_points(n)
end subroutine chebyshev_points_for_c
    
end module chebyshev_points_for_c_mod