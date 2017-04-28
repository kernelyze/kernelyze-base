! kernel_blackscholes.f90
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
! created on: 2015-04-03
! updated on: 2015-07-03
! updated on: 2016-08-26 (use black_formula)
! updated on: 2016-09-27 (added is_call flag)
! updated on: 2016-09-30 (use inputs directly rather than exponentiating,
!             make default sigma value more reasonable)
! updated on: 2016-10-02 (revert to exponentiating inputs -- this makes
!             sense when considering scenario shifts, where parallel
!             shifts in log forwards are much more plausible than
!             parallel shifts in forwards themselves)
!
! This module contains a derived type that extends the 
! kernel_summable abstract derived type and makes it concrete:
! the extension here is to a Black-Scholes kernel.
  
module kernel_blackscholes_mod

  use set_precision, only : wp
  use kernel_summable_mod, only : kernel_summable
  use option_formulae_mod, only : black_formula

  implicit none

  private

  type, public, extends(kernel_summable) :: kernel_blackscholes
    ! The key constant "sigma" in the definition of
    ! the Black-Scholes kernel
    real(wp), private :: sigma = 1E-1_wp ! Default value
    ! Is this a call?
    logical, private  :: is_call = .true. ! Default value
  contains
    ! The evaluation procedure for a blackscholes
    ! product kernel
    procedure, pass :: eval => eval_blackscholes_kernel
    ! Get the parameter "alpha"
    procedure, pass :: get_sigma => get_sigma_param
    ! Set the parameter "alpha"
    procedure, pass :: set_sigma => set_sigma_param
    ! Get the is_call flag
    procedure, pass :: get_is_call
    ! Set the is_call flag
    procedure, pass :: set_is_call
  end type kernel_blackscholes
  
  contains
  
  pure recursive function eval_blackscholes_kernel(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_blackscholes), intent(in)  :: this
    real(wp), intent(in)                    :: x
    real(wp), intent(in)                    :: y
    ! Function result
    real(wp)                                :: eval_res
    ! Body
    eval_res = black_formula( &
        is_call = this%is_call, &
        strike = exp( x ), &
        forward = exp( y ), &
        vol = this%sigma, &
        disc_fac = 1E0_wp )
  end function eval_blackscholes_kernel
  
  pure function get_sigma_param(this) result(sigma_param)
    ! Argument
    class(kernel_blackscholes), intent(in)  :: this
    ! Function result
    real(wp)                                :: sigma_param
    ! Body
    sigma_param = this%sigma
  end function get_sigma_param
  
  pure subroutine set_sigma_param(this, sigma_param)
    ! Arguments
    class(kernel_blackscholes), intent(inout) :: this
    real(wp), intent(in)                      :: sigma_param
    ! Body
    this%sigma = sigma_param
  end subroutine set_sigma_param
  
  pure function get_is_call(this) result(call_flag)
    ! Argument
    class(kernel_blackscholes), intent(in)  :: this
    ! Function result
    logical                                 :: call_flag
    ! Body
    call_flag = this%is_call
  end function get_is_call
  
  pure subroutine set_is_call(this, call_flag)
    ! Arguments
    class(kernel_blackscholes), intent(inout) :: this
    logical, intent(in)                       :: call_flag
    ! Body
    this%is_call = call_flag
  end subroutine set_is_call

end module kernel_blackscholes_mod