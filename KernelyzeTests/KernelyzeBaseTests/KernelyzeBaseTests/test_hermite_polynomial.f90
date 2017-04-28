! test_hermite_polynomial.f90
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
! created on: 2016-02-11
! updated on: 2016-02-11
!
! A module testing the hermite_polynomial derived type.

module test_hermite_polynomial_mod

  use set_precision, only : wp
  use hermite_polynomial_mod, only : hermite_polynomial
  
  implicit none
  
  contains
  
  function test_hermite_polynomial(unit_to_use) result(test_pass)
    ! Arguments
    integer, intent(in)                       :: unit_to_use
    ! Return value
    logical                                   :: test_pass

    write(unit_to_use , *) 'Test of functionality to evaluate ' &
        // 'Hermite polynomials (physicists version): '
  
    test_pass = .true.
  
    test_pass = test_pass .and. test_hermite_index_two(unit_to_use)
    test_pass = test_pass .and. test_hermite_low_index(unit_to_use)
    test_pass = test_pass .and. test_hermite_high_index(unit_to_use)
  
    write(unit_to_use , *) ' ' ! Output a blank line
  
  end function test_hermite_polynomial
  
  function test_hermite_index_two(unit_to_use) result(test_pass)
    ! Arguments
    integer, intent(in)       :: unit_to_use
    ! Result
    logical                   :: test_pass
    ! Local variables
    real(wp), dimension(10)   :: test_arr, func_arr
    real(wp)                  :: x
    integer                   :: j
    ! Body
    write(unit_to_use , *) 'Test of n == 2 Hermite polynomial: '
    do j = 1, size(test_arr)
      test_arr(j) = (real(j, wp) - real(size(test_arr) + 1, wp) / 2E0_wp) &
          / ( real(size(test_arr) + 1, wp) / 2E0_wp )
    end do
    write(unit_to_use , *) 'The test array is: '
    write(unit_to_use , '(*(2x,F9.6))') test_arr
    do j = 1, size(test_arr)
      func_arr(j) = hermite_polynomial(2, test_arr(j))
    end do
    write(unit_to_use , *) 'The result array is: '
    write(unit_to_use , '(*(2x,F9.6))') func_arr
    test_pass = .true.
    ! The index-4 Hermite polynomial (physicist's version) is:
    ! 4 x**2 - 2
    do j = 1, size(test_arr)
      x = test_arr(j)
      if ( abs(func_arr(j) &
                - ( 4E0_wp * (x**2) - 2E0_wp ) ) &
              > 10E0_wp * epsilon(func_arr(j)) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed, computed Hermite polynomial: ', &
            func_arr(j), ' but expected: ', &
            4E0_wp * (x**2) - 2E0_wp
      end if
    end do
  end function test_hermite_index_two

  function test_hermite_low_index(unit_to_use) result(test_pass)
    ! Arguments
    integer, intent(in)       :: unit_to_use
    ! Result
    logical                   :: test_pass
    ! Local variables
    real(wp), dimension(10)   :: test_arr, func_arr
    real(wp)                  :: x
    integer                   :: j
    ! Body
    write(unit_to_use , *) 'Test of n == 4 Hermite polynomial: '
    do j = 1, size(test_arr)
      test_arr(j) = (real(j, wp) - real(size(test_arr) + 1, wp) / 2E0_wp) &
          / ( real(size(test_arr) + 1, wp) / 2E0_wp )
    end do
    write(unit_to_use , *) 'The test array is: '
    write(unit_to_use , '(*(2x,F9.6))') test_arr
    do j = 1, size(test_arr)
      func_arr(j) = hermite_polynomial(4, test_arr(j))
    end do
    write(unit_to_use , *) 'The result array is: '
    write(unit_to_use , '(*(2x,F15.6))') func_arr
    test_pass = .true.
    ! The index-4 Hermite polynomial (physicist's version) is:
    ! 16 * x**4 - 48 x**2 + 12
    do j = 1, size(test_arr)
      x = test_arr(j)
      if ( abs(func_arr(j) &
                - ( 16E0_wp * (x**4) - 48E0_wp * (x**2) + 12E0_wp ) ) &
              > 100E0_wp * epsilon(func_arr(j)) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed, computed Hermite polynomial: ', &
            func_arr(j), ' but expected: ', &
            16E0_wp * (x**4) - 48E0_wp * (x**2) + 12E0_wp
      end if
    end do
  end function test_hermite_low_index
  
  function test_hermite_high_index(unit_to_use) result(test_pass)
    ! Arguments
    integer, intent(in)       :: unit_to_use
    ! Result
    logical                   :: test_pass
    ! Local variables
    real(wp), dimension(10)   :: test_arr, func_arr
    real(wp)                  :: x
    integer                   :: j
    ! Body
    write(unit_to_use , *) 'Test of n == 10 Hermite polynomial: '
    do j = 1, size(test_arr)
      test_arr(j) = (real(j, wp) - real(size(test_arr) + 1, wp) / 2E0_wp) &
          / ( real(size(test_arr) + 1, wp) / 2E0_wp )
    end do
    write(unit_to_use , *) 'The test array is: '
    write(unit_to_use , '(*(2x,F9.6))') test_arr
    do j = 1, size(test_arr)
      func_arr(j) = hermite_polynomial(10, test_arr(j))
    end do
    write(unit_to_use , *) 'The result array is: '
    write(unit_to_use , '(*(2x,F15.6))') func_arr
    test_pass = .true.
    ! The index-10 Hermite polynomial (physicist's version) is:
    ! 1024 * x**10 - 23040 * x**8 + 161280 * x**6 
    ! - 403200 * x**4 + 302400 * x**2 - 30240
    do j = 1, size(test_arr)
      x = test_arr(j)
      if ( abs(func_arr(j) &
                - (   1024E0_wp   * (x**10) &
                    - 23040E0_wp  * (x**8) &
                    + 161280E0_wp * (x**6) &
                    - 403200E0_wp * (x**4) &
                    + 302400E0_wp * (x**2) &
                    - 30240E0_wp ) ) &
              > 1E5_wp * epsilon(func_arr(j)) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed, computed Hermite polynomial: ', &
            func_arr(j), ' but expected: ', &
              1024E0_wp   * (x**10) &
            - 23040E0_wp  * (x**8) &
            + 161280E0_wp * (x**6) &
            - 403200E0_wp * (x**4) &
            + 302400E0_wp * (x**2) &
            - 30240E0_wp
      end if
    end do
  end function test_hermite_high_index

end module test_hermite_polynomial_mod
