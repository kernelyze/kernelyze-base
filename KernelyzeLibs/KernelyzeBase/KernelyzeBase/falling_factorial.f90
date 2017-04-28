! falling_factorial.f90
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
! falling factorials.
  
module falling_factorial_mod

use set_precision, only : wp

implicit none

private

public :: falling_factorial

contains
  
! This function computes the falling factorial:
! x * (x - 1) * (x - 2) * . . . * (x - i + 1)
pure function falling_factorial(x, i) result(fact)
  ! Arguments
  real(wp), intent(in)  :: x
  integer, intent(in)   :: i
  ! Result
  real(wp)              :: fact
  ! Local variables
  integer               :: j
  ! Body
  ! If i is 0, then the result is 1 (this is useful
  ! in, e. g., using the falling factorial when
  ! taking the i^th derivative of x ** n)
  if (i == 0) then
    fact = 1E0_wp
    return
  end if
  ! If I am here, then i /= 0
  fact = x
  do j = 1, i - 1
    fact = fact * ( x - real(j , wp) )
  end do
end function falling_factorial

end module falling_factorial_mod