! itoa.f90
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
! created on: 2016-04-19
! updated on: 2016-04-19
! updated on: 2016-04-20 (streamlined)
! updated on: 2016-04-21 (fixed spacing)
!
! This module contains a simple pure convenience function
! that, given an integer, returns a string.

module itoa_mod

implicit none

private

public :: itoa

contains
  
pure function itoa(i) result(a)
  ! Arguments
  integer, intent(in)           :: i
  ! Result
  character(len=:), allocatable :: a
  ! Body
  allocate(character(len=(range(i) + 2)) :: a)
  write(a, '(i0)') i
  a = trim(a)
end function itoa
  
end module itoa_mod