! test_elemental.f90
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
! created on: 2015-04-04
! updated on: 2015-04-04
!
! A module testing simple applications of elemental procedures.
! This is embarrassingly simple but clarifies potentially
! interesting applications.
  
module test_elemental_mod

  use set_precision, only : wp
  
  implicit none
  
  contains
  
  real(wp) elemental function elemental_example(x, y)
    ! Arguments
    real(wp), intent(in)  :: x
    real(wp), intent(in)  :: y
    ! Body
    elemental_example = x + y
  end function elemental_example

  logical function test_elemental(unit_to_use)
    ! Arguments
    integer, intent(in)     :: unit_to_use
    ! Local variables
    real(wp), dimension(10) :: a_arr, b_arr, res_arr
    real(wp)                :: a, b, res
    integer                 :: j
    ! Body
    write(unit_to_use , *) 'Test of simple elemental function: '
    a = 2E0_wp
    b = -1E0_wp
    do j = 1, size(a_arr)
      a_arr(j) = - real(j, wp) + 5E0_wp
      b_arr(j) = real(j, wp)
    end do
    write(unit_to_use , *) 'The a variable is: ', a
    write(unit_to_use , *) 'The b variable is: ', b
    write(unit_to_use , *) 'The a array variable is: '
    write(unit_to_use , '(*(2x,F10.7))') a_arr
    write(unit_to_use , *) 'The b array variable is: '
    write(unit_to_use , '(*(2x,F10.7))') b_arr
    res = elemental_example(a , b)
    write(unit_to_use , *) 'An elemental add of a to b gives: ', res
    res_arr = elemental_example(a_arr , b)
    write(unit_to_use , *) 'An elemental add of a_arr to b gives: '
    write(unit_to_use , '(*(2x,F10.7))') res_arr
    res_arr = elemental_example(a , b_arr)
    write(unit_to_use , *) 'An elemental add of a to b_arr gives: '
    write(unit_to_use , '(*(2x,F10.7))') res_arr
    res_arr = elemental_example(a_arr , b_arr)
    write(unit_to_use , *) 'An elemental add of a_arr to b_arr gives: '
    write(unit_to_use , '(*(2x,F10.7))') res_arr
    test_elemental = .true.
  end function test_elemental

end module test_elemental_mod