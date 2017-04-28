! kernel_cauchy.f90
!
! Copyright (c) 2015, 2016, 2017 by Kernelyze LLC
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
! updated on: 2016-01-09
! updated on: 2016-02-07 (added higher-order analytical partials)
! updated on: 2017-04-23 (remove map-maker)
! updated on: 2017-04-27 (interfaces adjusted for greater
!             flexibility for other kernels without analytical partials)
!
! This module contains a derived type that extends the 
! kernel_map_maker abstract derived type and makes it concrete:
! the extension here is to a (generalized) Cauchy kernel
! of the form $ \frac{1}{ \left( \left( c + x + y \right)^{\alpha} \right) }$.
  
module kernel_cauchy_mod

  use set_precision, only : wp
  use falling_factorial_mod, only : falling_factorial
  use kernel_summable_mod, only : kernel_summable
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  private

  type, public, extends(kernel_summable) :: kernel_cauchy
    ! The first key constant "c" in the definition of the
    ! Cauchy kernel
    real(wp), private :: c = 3E0_wp ! Default value
    ! The second key constant "alpha" in the definition of
    ! the Cauchy kernel
    real(wp), private :: alpha = 1E0_wp ! Default value
  contains
    ! The evaluation procedure for a Cauchy
    ! product kernel
    procedure, pass :: eval => eval_cauchy_kernel
    ! First and second partial derivatives
    procedure, pass :: first_partial_dx
    procedure, pass :: first_partial_dy
    procedure, pass :: second_partial_dx_dx
    procedure, pass :: second_partial_dy_dy
    procedure, pass :: second_partial_dx_dy
    ! Higher-order partial derivatives
    procedure, pass :: mth_nth_partial
    ! Get the parameter "c"
    procedure, pass :: get_c => get_c_param
    ! Set the parameter "c"
    procedure, pass :: set_c => set_c_param
    ! Get the parameter "alpha"
    procedure, pass :: get_alpha => get_alpha_param
    ! Set the parameter "alpha"
    procedure, pass :: set_alpha => set_alpha_param
  end type kernel_cauchy
  
  contains
  
  pure recursive function eval_cauchy_kernel(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_cauchy), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    ! Function result
    real(wp)                            :: eval_res
    ! Local variable
    real(wp)                            :: z
    ! Body
    z = this%c + x + y
    if (z > epsilon(z)) then
      eval_res = 1E0_wp / (z**this%alpha)
    else
      eval_res = ieee_value(eval_res, ieee_quiet_nan)
    end if
  end function eval_cauchy_kernel
  
  elemental function first_partial_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = ( -this%alpha ) * (this%c + x + y)**( -(this%alpha + 1) )
  end function first_partial_dx
  
  elemental function first_partial_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%first_partial_dx( x , y )
  end function first_partial_dy
  
  elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = ( this%alpha ) * ( this%alpha + 1 ) * &
        (this%c + x + y)**( -(this%alpha + 2) )
  end function second_partial_dx_dx
  
  elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%second_partial_dx_dx( x , y )
  end function second_partial_dy_dy
  
  elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%second_partial_dx_dx( x , y )
  end function second_partial_dx_dy
  
  pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
    ! Arguments
    class(kernel_cauchy), intent(in) :: this
    integer, intent(in)       :: m
    integer, intent(in)       :: n
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = falling_factorial(-this%alpha, m + n) * &
        (this%c + x + y)**( -(this%alpha + real(m + n, wp) ) )
  end function mth_nth_partial
  
  pure function get_c_param(this) result(c_param)
    ! Argument
    class(kernel_cauchy), intent(in)  :: this
    ! Function result
    real(wp)                          :: c_param
    ! Body
    c_param = this%c
  end function get_c_param
  
  pure subroutine set_c_param(this, c_param)
    ! Arguments
    class(kernel_cauchy), intent(inout)   :: this
    real(wp), intent(in)                  :: c_param
    ! Body
    this%c = c_param
  end subroutine set_c_param
  
  pure function get_alpha_param(this) result(alpha_param)
    ! Argument
    class(kernel_cauchy), intent(in)  :: this
    ! Function result
    real(wp)                          :: alpha_param
    ! Body
    alpha_param = this%alpha
  end function get_alpha_param
  
  pure subroutine set_alpha_param(this, alpha_param)
    ! Arguments
    class(kernel_cauchy), intent(inout)   :: this
    real(wp), intent(in)                  :: alpha_param
    ! Body
    this%alpha = alpha_param
  end subroutine set_alpha_param

end module kernel_cauchy_mod