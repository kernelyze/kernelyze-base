! monomial.f90
!
! Copyright (c) 2016, 2017 by Kernelyze LLC
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
! created on: 2016-02-07
! updated on: 2016-02-07
! updated on: 2016-02-08 (eliminated "use" of unnecessary modules)
! updated on: 2016-02-10 (added recentering)
! updated on: 2016-02-19 (fixed get_center type-bound procedure)
! updated on: 2016-11-05 (conform to new 
!             calculate_integral_kernel_svar; this allows
!             smooth use in template_kernelrankn and sets
!             the stage for potential future functionality)
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
! updated on: 2017-04-23 (remove use statements for
!             modules or parts of modules that are
!             not actually used)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! This module contains a derived type extending
! the singlevar_func abstract derived type. The
! derived type here provides a monomial, that is,
! a function of the form x**n for some positive
! integer n.
  
module monomial_mod

  use set_precision, only : wp
  use singlevar_func_mod, only : singlevar_func
  use falling_factorial_mod, only : falling_factorial

  implicit none

  private

  type, public, extends(singlevar_func) :: monomial
    ! The power of the monomial; that is, n in the
    ! functional form x**n
    integer, private    :: mon_power
    ! The centering of the monomial, that is, c in the
    ! functional form (x - c)**n
    real(wp), private   :: center
  contains
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of singlevar_func
    procedure, pass :: eval
    ! Override the singlevar_func derivative procedures
    procedure, pass :: first_deriv
    procedure, pass :: second_deriv
    procedure, pass :: nth_deriv
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => monomial_assign
    ! Inspector function for the monomial power
    procedure, pass :: get_power
    ! Setter procedure for the monomial power
    procedure, pass :: set_power
    ! Inspector function for the monomial center
    procedure, pass :: get_center
    ! Setter procedure for the monomial center
    procedure, pass :: set_center
  end type monomial

  contains
  
  ! Implement the deferred procedure for non-elemental
  ! (but still pure) evaluation.
  pure recursive function eval(this, arg_pt) result(eval_res)
    ! Arguments
    class(monomial), intent(in) :: this
    real(wp), intent(in)        :: arg_pt
    ! Function result
    real(wp)                    :: eval_res
    ! Body
    ! I must guard against 0**0 issues; if 0**0 is
    ! encountered, the result should be 1E0_wp
    if (abs(arg_pt - this%center) < epsilon(arg_pt) &
        .and. this%mon_power == 0) then
      eval_res = 1E0_wp
    else
      eval_res = (arg_pt - this%center) ** this%mon_power
    end if
  end function eval
  ! First derivative
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(monomial), intent(in) :: this
    real(wp), intent(in)        :: x
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                    :: deriv
    ! Body
    ! I must guard against 0**0 issues; if 0**0 is
    ! encountered, the result should be 1E0_wp
    if (abs(x - this%center) < epsilon(x) &
        .and. this%mon_power - 1 == 0) then
      deriv = real(this%mon_power, wp)
    else
      deriv = real(this%mon_power, wp) * &
          ((x - this%center) ** (this%mon_power - 1) )
    end if
  end function first_deriv
  ! Second derivative
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(monomial), intent(in) :: this
    real(wp), intent(in)        :: x
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                    :: deriv
    ! Body
    ! I must guard against 0**0 issues; if 0**0 is
    ! encountered, the result should be 1E0_wp
    if (abs(x - this%center) < epsilon(x) &
        .and. this%mon_power - 2 == 0) then
      deriv = real( this%mon_power * (this%mon_power - 1) , wp )
    else
      deriv = real( this%mon_power * (this%mon_power - 1) , wp ) * &
          ((x - this%center) ** (this%mon_power - 2) )
    end if    
  end function second_deriv
  ! nth derivative, using the kernel being sectioned
  pure function nth_deriv(this, n, x, d) result(deriv)
    ! Arguments
    class(monomial), intent(in) :: this
    integer, intent(in)         :: n
    real(wp), intent(in)        :: x
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                    :: deriv
    ! Body
    ! I must guard against 0**0 issues; if 0**0 is
    ! encountered, the result should be 1E0_wp
    if (abs(x - this%center) < epsilon(x) &
        .and. this%mon_power - n == 0) then
      deriv = falling_factorial( real( this%mon_power , wp ) , n )
    else
      deriv = falling_factorial( real( this%mon_power , wp ) , n ) * &
          ((x - this%center) ** (this%mon_power - n) )
    end if
  end function nth_deriv
  ! Assignment
  subroutine monomial_assign(left, right)
    ! Arguments
    class(monomial), intent(out)      :: left
    class(singlevar_func), intent(in) :: right
    ! Body
    select type (right)
    type is (monomial)
      left%mon_power = right%mon_power
      left%center    = right%center
    end select
  end subroutine monomial_assign
  ! Inspector for the monomial power
  pure function get_power(this) result(curr_power)
    ! Argument
    class(monomial), intent(in)  :: this
    ! Function result
    integer                      :: curr_power
    ! Body
    curr_power = this%mon_power
  end function get_power
  ! Setter for the monomial power
  pure subroutine set_power(this, new_power)
    ! Arguments
    class(monomial), intent(inout)  :: this
    integer, intent(in)             :: new_power
    ! Body
    this%mon_power = new_power
  end subroutine set_power
  ! Inspector for the monomial center
  pure function get_center(this) result(curr_center)
    ! Argument
    class(monomial), intent(in)  :: this
    ! Function result
    real(wp)                     :: curr_center
    ! Body
    curr_center = this%center
  end function get_center
  ! Setter for the monomial center
  pure subroutine set_center(this, new_center)
    ! Arguments
    class(monomial), intent(inout)  :: this
    real(wp), intent(in)            :: new_center
    ! Body
    this%center = new_center
  end subroutine set_center

end module monomial_mod
