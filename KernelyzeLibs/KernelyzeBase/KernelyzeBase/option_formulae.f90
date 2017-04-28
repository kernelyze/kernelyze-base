! option_formulae.f90
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
! created on: 2016-08-19
! updated on: 2016-08-21
!
! This module contains option-pricing formulae.
  
module option_formulae_mod

use set_precision, only : wp
use normal_distn, only : std_normal_pdf, std_normal_cdf
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none

private

public :: black_formula
public :: normopt_formula
! Separate elemental versions (as in, e. g., kernel derived type)
public :: black_formula_elt
public :: normopt_formula_elt

contains
  
pure function black_formula( &
    is_call, & ! .true. for call, .false. for put
    strike, &  
    forward, & 
    vol, &     ! vol of option, including time: sigma * sqrt(T)
    disc_fac ) result(price)
  ! Arguments
  logical, intent(in)   :: is_call
  real(wp), intent(in)  :: strike
  real(wp), intent(in)  :: forward
  real(wp), intent(in)  :: vol
  real(wp), intent(in)  :: disc_fac
  ! Result
  real(wp)              :: price
  ! Local variables
  real(wp)  :: d1
  real(wp)  :: d2
  ! Body
  ! If vol is negative by more than epsilon, this is an error: return NaN
  if (vol < - epsilon(vol)) then
    price = ieee_value(price, ieee_quiet_nan)
    return
  end if
  ! If vol is less than epsilon, return discounted intrinsic
  if (vol < epsilon(vol)) then
    if (is_call) then
      price = disc_fac * max( 0E0_wp , forward - strike )
    else
      price = disc_fac * max( 0E0_wp , strike - forward )
    end if
    return
  end if
  ! At this point, I know that vol is greater than epsilon
  d1 = ( log(forward / strike) / vol ) + ( 5E-1_wp * vol )
  d2 = d1 - vol
  if (is_call) then
    price = forward * std_normal_cdf( d1 ) - strike * std_normal_cdf( d2 )
  else
    price = strike * std_normal_cdf( - d2 ) - forward * std_normal_cdf( - d1 )
  end if
  price = disc_fac * price
end function black_formula
    
elemental function black_formula_elt( &
    is_call, & ! .true. for call, .false. for put
    strike, &  
    forward, & 
    vol, &     ! vol of option, including time: sigma * sqrt(T)
    disc_fac ) result(price)
  ! Arguments
  logical, intent(in)   :: is_call
  real(wp), intent(in)  :: strike
  real(wp), intent(in)  :: forward
  real(wp), intent(in)  :: vol
  real(wp), intent(in)  :: disc_fac
  ! Result
  real(wp)              :: price
  ! Body
  price = black_formula(is_call, strike, forward, vol, disc_fac)
end function black_formula_elt
    
pure function normopt_formula( &
    is_call, & ! .true. for call, .false. for put
    strike, &  
    forward, & 
    vol, &     ! vol of option, including time: sigma * sqrt(T)
    disc_fac ) result(price)
  ! Arguments
  logical, intent(in)   :: is_call
  real(wp), intent(in)  :: strike
  real(wp), intent(in)  :: forward
  real(wp), intent(in)  :: vol
  real(wp), intent(in)  :: disc_fac
  ! Result
  real(wp)              :: price
  ! Local variables
  real(wp)  :: d
  real(wp)  :: h
  ! Body
  ! If vol is negative by more than epsilon, this is an error: return NaN
  if (vol < - epsilon(vol)) then
    price = ieee_value(price, ieee_quiet_nan)
    return
  end if
  if (is_call) then
    d = forward - strike
  else
    d = strike - forward
  end if
  ! If vol is less than epsilon, return discounted intrinsic
  if (vol < epsilon(vol)) then
    price = disc_fac * max( 0E0_wp , d )
    return
  end if
  ! At this point, I know that vol is greater than epsilon
  h = d / vol
  price = vol * std_normal_pdf( h ) + d * std_normal_cdf( h )
  price = disc_fac * price
end function normopt_formula
    
elemental function normopt_formula_elt( &
    is_call, & ! .true. for call, .false. for put
    strike, &  
    forward, & 
    vol, &     ! vol of option, including time: sigma * sqrt(T)
    disc_fac ) result(price)
  ! Arguments
  logical, intent(in)   :: is_call
  real(wp), intent(in)  :: strike
  real(wp), intent(in)  :: forward
  real(wp), intent(in)  :: vol
  real(wp), intent(in)  :: disc_fac
  ! Result
  real(wp)              :: price
  ! Body
  price = normopt_formula(is_call, strike, forward, vol, disc_fac)
end function normopt_formula_elt

end module option_formulae_mod