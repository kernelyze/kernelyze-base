! test_chebyshev_points.f90
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
! created on: 2015-02-16
! updated on: 2015-04-22
!
! A module of unit tests for the chebyshev_points module.
  
module test_chebyshev_points_mod

use set_precision, only : wp
use chebyshev_points_mod, only : chebyshev_points
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

implicit none

contains

function test_chebyshev_points(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  write(unit_to_use , *) 'Test of Chebyshev-points functionality'
  
  test_pass = .true.
  
  test_pass = test_chebyshev_points_small(unit_to_use)
  test_pass = test_pass .and. test_chebyshev_points_large(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line

end function test_chebyshev_points

function test_chebyshev_points_small(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)       :: unit_to_use
  ! Result
  logical                   :: test_pass
  
  ! Local variables
  real(wp) , dimension(10)  :: cheb_pts
  integer                   :: j
  
  ! Check on the first 10 Chebyshev points
  cheb_pts = chebyshev_points(10)
  write(unit_to_use, *) 'The first 10 Chebyshev points are: '
  do j = 1,size(cheb_pts)
    write(unit_to_use, *) cheb_pts(j)
  end do
  test_pass = .not. (any(ieee_is_nan(cheb_pts)))
end function test_chebyshev_points_small

function test_chebyshev_points_large(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Result
  logical                     :: test_pass
  
  ! Local variables
  real(wp) , dimension(200)   :: cheb_pts
  integer                     :: j
  
  ! Check on the first 200 Chebyshev points
  cheb_pts = chebyshev_points(200)
  write(unit_to_use, *) 'The first 200 Chebyshev points are: '
  do j = 1,size(cheb_pts)
    write(unit_to_use, *) cheb_pts(j)
  end do
  test_pass = .not. (any(ieee_is_nan(cheb_pts)))
end function test_chebyshev_points_large

end module test_chebyshev_points_mod