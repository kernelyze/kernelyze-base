! matmul_for_c.f90
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
!
! created on: 2017-04-23
!
! A module that gives C types and a binding for an interface to
! the matmul intrinsic of Fortran.
  
module matmul_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use, intrinsic :: iso_c_binding

implicit none

  contains
  
! Compute A B = D, where A is m by n, B is n by p, and D is thus
! m by p.
pure subroutine matmul_for_c( &
    c_m, & ! The number of rows of A (and rows of D)
    c_n, & ! The number of columns of A (and rows of B)
    c_p, & ! The number of columns of B (and columns of D)
    c_a, & ! The A matrix
    c_b, & ! The B matrix
    c_d) & ! The D matrix
    bind(C, name='matmul')
  ! Arguments
  integer(c_int), intent(in)    :: c_m
  integer(c_int), intent(in)    :: c_n
  integer(c_int), intent(in)    :: c_p
  real(c_double), intent(in)    :: c_a(c_m, c_n)
  real(c_double), intent(in)    :: c_b(c_n, c_p)
  real(c_double), intent(out)   :: c_d(c_m, c_p)
  ! Body
  c_d = matmul(c_a, c_b)
end subroutine matmul_for_c

end module matmul_for_c_mod  