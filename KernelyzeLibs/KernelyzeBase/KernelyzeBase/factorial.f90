! factorial.f90
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
! updated on: 2016-02-08
!
! This module contains functionality to compute
! factorials.
  
module factorial_mod

use set_precision, only : wp

implicit none

private

public :: factorial

contains
  
! This function computes the factorial:
! n * (n - 1) * (n - 2) * . . . * 2 * 1
pure function factorial(n) result(fact)
  ! Arguments
  integer, intent(in)   :: n
  ! Result
  real(wp)              :: fact
  ! Local variables
  integer               :: j
  ! Body
  fact = 1E0_wp
  do j = 1, n
    fact = fact * real(j , wp)
  end do
end function factorial

end module factorial_mod