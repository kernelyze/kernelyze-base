! test_set_precision.f90
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
! created on: 2015-02-25
! updated on: 2015-02-25
!
! A module of unit tests for the set_precision module.
  
module test_set_precision_mod

use set_precision, only : wp

implicit none

contains
  
function test_set_precision(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Test of constants
  write(unit_to_use , *) 'Test of set_precision'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_set_precision_wp(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_set_precision

function test_set_precision_wp(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp)                  :: vbl_wp_precision
  integer                   :: kind_of_vbl
  ! Body
  test_pass = .true.
  
  write(unit_to_use , *) 'Test the precision of a real(wp) variable: '
  
  vbl_wp_precision = 1.0_wp
  kind_of_vbl = kind(vbl_wp_precision)
  write(unit_to_use , *) 'The kind type parameter of a real(wp) variable: ' , &
      kind_of_vbl
  ! The test passes if wp represents double precision; of course, if 
  ! I ever changed this to single or quadruple precision this test would fail,
  ! but 1) I would want to be alerted to "double-check" such a change and 2)
  ! it would be very easy to modify this test in the event of such a
  ! change.
  test_pass = test_pass .and. (kind_of_vbl == kind(1.0D0))

end function test_set_precision_wp

end module test_set_precision_mod