! option_formulae_for_c.f90
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
! A module that gives C types and a binding for functionality involving
! option-pricing formulae.
  
module option_formulae_for_c_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use option_formulae_mod, only : black_formula, normopt_formula
use fort_c_interop_mod, only : copy_string
use, intrinsic :: iso_c_binding

implicit none

contains
  
pure function black_formula_for_c( &
    c_is_call, &
    c_strike, &
    c_forward, &
    c_vol, &
    c_disc_fac ) &
    bind(C, name='black_formula') &
    result(price)
  ! Arguments
  logical(c_bool), intent(in) :: c_is_call
  real(c_double), intent(in)  :: c_strike
  real(c_double), intent(in)  :: c_forward
  real(c_double), intent(in)  :: c_vol
  real(c_double), intent(in)  :: c_disc_fac
  ! Result
  real(c_double)              :: price
  ! Body
  price = black_formula( &
      logical(c_is_call), &
      c_strike, &
      c_forward, &
      c_vol, &
      c_disc_fac )
end function black_formula_for_c
    
pure function normopt_formula_for_c( &
    c_is_call, &
    c_strike, &
    c_forward, &
    c_vol, &
    c_disc_fac ) &
    bind(C, name='normopt_formula') &
    result(price)
  ! Arguments
  logical(c_bool), intent(in) :: c_is_call
  real(c_double), intent(in)  :: c_strike
  real(c_double), intent(in)  :: c_forward
  real(c_double), intent(in)  :: c_vol
  real(c_double), intent(in)  :: c_disc_fac
  ! Result
  real(c_double)              :: price
  ! Body
  price = normopt_formula( &
      logical(c_is_call), &
      c_strike, &
      c_forward, &
      c_vol, &
      c_disc_fac )
end function normopt_formula_for_c    
  
end module option_formulae_for_c_mod