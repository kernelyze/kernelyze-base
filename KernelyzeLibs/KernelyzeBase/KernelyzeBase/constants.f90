! constants.f90
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
! created on: 2015-02-07
! updated on: 2015-04-19
! updated on: 2016-10-03 (added sqrtpi and inv_sqrtpi)
! updated on: 2016-10-06 (added twopi)
!
! A simple utility module containing frequently-used constants
  
module constants_mod
  use set_precision, only : wp
  implicit none
  private ! Make everything in this module private by default
  public pi, twopi, sqrt2, twice_sqrt2, sqrtpi, inv_sqrtpi, sqrt2pi, &
      inv_sqrt2pi, exp_of_one, alloc_errmsg_len, description_len, err_msg_len
  
  real(wp), parameter :: pi               = 4.0E0_wp * atan(1.0E0_wp)
  real(wp), parameter :: twopi            = 2E0_wp * pi
  real(wp), parameter :: sqrt2            = sqrt(2E0_wp)
  real(wp), parameter :: twice_sqrt2      = 2E0_wp * sqrt2
  real(wp), parameter :: sqrtpi           = sqrt( pi )
  real(wp), parameter :: inv_sqrtpi       = 1E0_wp / sqrtpi
  real(wp), parameter :: sqrt2pi          = sqrt(twopi)
  real(wp), parameter :: inv_sqrt2pi      = 1E0_wp / sqrt2pi
  real(wp), parameter :: exp_of_one       = exp(1.0E0_wp)
  integer, parameter  :: alloc_errmsg_len = 255
  integer, parameter  :: description_len  = 255
  integer, parameter  :: err_msg_len      = 255
  
end module constants_mod