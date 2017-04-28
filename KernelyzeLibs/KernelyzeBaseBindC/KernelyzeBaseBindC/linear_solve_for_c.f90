! linear_solve_for_c.f90
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
! the LAPACK subroutine DGESV.
  
module linear_solve_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use tk_simple_lapack_mod, only : tk_gesv
use, intrinsic :: iso_c_binding

implicit none

contains
  
! Solve A X = B, where A is n by n and B is n by nrhs  
pure subroutine linear_solve_for_c( &
    c_n, & ! The number of rows (== number of columns) of the A matrix
    c_a, & ! The A matrix
    c_nrhs, & ! The number of columns of the B matrix
    c_b, & ! The B matrix
    c_x, & ! The solved X matrix
    c_info ) &
    bind(C, name='linear_solve')
  ! Arguments
  integer(c_int), intent(in)    :: c_n
  real(c_double), intent(in)    :: c_a(c_n, c_n)
  integer(c_int), intent(in)    :: c_nrhs
  real(c_double), intent(in)    :: c_b(c_n, c_nrhs)
  real(c_double), intent(out)   :: c_x(c_n, c_nrhs)
  integer(c_int), intent(out)   :: c_info
  ! Local variables
  real(wp)  :: temp_a(c_n, c_n)
  integer   :: ipiv(c_n)
  ! Body
  ! On entry, c_x will be the B matrix; on exit, it will be the X matrix
  c_x = c_b 
  ! The A matrix gets overwritten in the LAPACK call, so copy the
  ! intent(in) argument
  temp_a = c_a
  call tk_gesv(temp_a, c_x, ipiv, c_info)
end subroutine linear_solve_for_c

end module linear_solve_for_c_mod  
  