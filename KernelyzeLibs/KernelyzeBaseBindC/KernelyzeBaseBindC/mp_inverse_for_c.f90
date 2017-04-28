! mp_inverse_for_c.f90
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
! created on: 2015-05-31
! updated on: 2015-06-03
!
! A module that gives C types and a binding for moore_penrose_inverse.

module mp_inverse_for_c_mod

use set_precision, only : wp
use moore_penrose_inverse_mod, only : moore_penrose_inverse
use, intrinsic :: iso_c_binding

implicit none

contains
  
pure subroutine mp_inverse_for_c(to_return, to_invert, nrows, ncols, toler) &
    bind(C, name='mp_inverse')

  ! Arguments
  ! Note the use of the "value" descriptor to be sure that
  ! values rather than references can be passed in from C for
  ! the nrows, ncols, and toler parameters.
  integer(c_int), intent(in), value                       :: nrows
  integer(c_int), intent(in), value                       :: ncols
  real(c_double), intent(out), dimension(nrows, ncols)    :: to_return
  real(c_double), intent(in), dimension(nrows, ncols)     :: to_invert
  real(c_double), intent(in), value                       :: toler
  ! Body
  ! Note the kind conversion, which may be necessary if Fortran
  ! working precision (wp) is not the same as C double precision
  ! (c_double).
  to_return = real( &
      moore_penrose_inverse(real(to_invert , wp), real(toler , wp)), &
      c_double)
end subroutine mp_inverse_for_c

end module mp_inverse_for_c_mod