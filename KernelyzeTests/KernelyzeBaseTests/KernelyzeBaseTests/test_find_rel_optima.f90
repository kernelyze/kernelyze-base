! test_find_rel_optima.f90
!
! Copyright (c) 2015, 2016 by Kernelyze LLC
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
! created on: 2015-02-16
! updated on: 2015-04-22
! updated on: 2016-04-25 (to use singlevar_funcs)
!
! A module of unit tests for the find_rel_optima module.
  
module test_find_rel_optima_mod

use set_precision, only : wp
use intrinsic_test_funcs_mod, only : sin_func, exp_func
use find_all_zeros_mod, only  : find_all_zeros
use find_rel_optima_mod, only  : find_rel_optima

implicit none

contains
  
function test_find_rel_optima(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Test search for relative optima
  write(unit_to_use , *) 'Test of find-relative-optima functionality'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_find_rel_optima_sin(unit_to_use)
  test_pass = test_pass .and. test_find_rel_optima_exp(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
end function test_find_rel_optima

function test_find_rel_optima_sin(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(21)                  :: grid_vec
  real(wp) , dimension(:) , allocatable     :: res_zeros
  real(wp) , dimension(: , :), allocatable  :: res_rel_optima
  integer                                   :: j
  ! Body
  test_pass = .true.
  
  ! Get zeros between which to find relative optima
  do j = 1,21
    grid_vec(j) = -10.0E0_wp + DBLE(j - 1) 
  end do
  call find_all_zeros(res_zeros, sin_func, grid_vec, 1.0E-15_wp)
  ! Now test find_rel_optima, using the zeros computed above
  call find_rel_optima(res_rel_optima, sin_func, res_zeros, 1.0E-15_wp)
  write(unit_to_use, *) 'The relative optima of sin(x) from -10.0 ' &
      // 'to 10.0 are (should be six of them): '
  do j = 1,size(res_rel_optima, 1)
    write(unit_to_use, *) 'x: ', res_rel_optima(j,1) , ' and -abs(value): ', res_rel_optima(j, 2)
  end do
  test_pass = test_pass .and. (size(res_rel_optima, 1) == 6)
end function test_find_rel_optima_sin

function test_find_rel_optima_exp(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  ! Local variables
  real(wp) , dimension(2)                   :: grid_vec
  real(wp) , dimension(: , :), allocatable  :: res_rel_optima
  integer                                   :: j
  ! Body
  test_pass = .true.
  
  ! Get zeros between which to find relative optima
  grid_vec(1) = -1.0E0_wp
  grid_vec(2) = 1.0E0_wp
  ! Now test find_rel_optima, using the zeros computed above
  call find_rel_optima(res_rel_optima, exp_func, grid_vec, 1.0E-15_wp)
  write(unit_to_use, *) 'The relative optima of exp(x) from -1.0 ' &
      // 'to 1.0 are (should be one of them): '
  do j = 1,size(res_rel_optima, 1)
    write(unit_to_use, *) 'x: ', res_rel_optima(j,1) , ' and -abs(value): ', res_rel_optima(j, 2)
  end do
  test_pass = test_pass .and. (size(res_rel_optima, 1) == 1)
end function test_find_rel_optima_exp

end module test_find_rel_optima_mod