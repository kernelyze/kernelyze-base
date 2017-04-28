! kernel_expprod.f90
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
! updated on: 2016-01-07
! updated on: 2016-02-06 (added higher-order analytical partials)
! updated on: 2017-04-23 (remove map-maker)
! updated on: 2017-04-27 (interfaces adjusted for greater
!             flexibility for other kernels without analytical partials)
!
! This module contains a derived type that extends the 
! kernel_map_maker abstract derived type and makes it concrete:
! the extension here is to an exponential product kernel
! of the form $ \exp{ \left( c x y \right) }$.
  
module kernel_expprod_mod

  use set_precision, only : wp
  use binomial_coeffs_mod, only : binomial_coeffs
  use falling_factorial_mod, only : falling_factorial
  use kernel_summable_mod, only : kernel_summable

  implicit none

  private

  type, public, extends(kernel_summable) :: kernel_expprod
    ! The key constant "c" in the definition of the
    ! exponential product kernel
    real(wp), private :: c = 1E0_wp ! Default value
  contains
    ! The evaluation procedure for an exponential
    ! product kernel
    procedure, pass :: eval => eval_expprod_kernel
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
  end type kernel_expprod
  
  contains
  
  pure recursive function eval_expprod_kernel(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)              :: x
    real(wp), intent(in)              :: y
    ! Function result
    real(wp)                          :: eval_res
    ! Body
    eval_res = exp(this%c * x * y)
  end function eval_expprod_kernel
  
  elemental function first_partial_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%c * y * exp(this%c * x * y)
  end function first_partial_dx
  
  elemental function first_partial_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%c * x * exp(this%c * x * y)
  end function first_partial_dy
  
  elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%c * this%c * y * y * exp(this%c * x * y)
  end function second_partial_dx_dx
  
  elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%c * this%c * x * x * exp(this%c * x * y)
  end function second_partial_dy_dy
  
  elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%c * (1E0_wp + this%c * x * y ) * exp(this%c * x * y)
  end function second_partial_dx_dy
  
  pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
    ! Arguments
    class(kernel_expprod), intent(in) :: this
    integer, intent(in)       :: m
    integer, intent(in)       :: n
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Local variables
    integer                   :: i
    real(wp)                  :: bin_coeffs(n + 1)
    real(wp)                  :: xterm
    real(wp)                  :: yterm
    real(wp)                  :: cmfactor
    ! Body
    ! I consider first taking the m^th derivative
    ! with respect to x to get:
    ! c^m y^m exp( c x y ).
    ! Then use the higher-order product rule to
    ! differentiate this with respect to y
    ! n times.
    bin_coeffs = binomial_coeffs(n)
    deriv = 0E0_wp
    do i = 1, min(m + 1, n + 1)
      ! It is important to guard against problems
      ! that arise when the ill-defined expression 0**0
      ! would arise in the calculation below.  If
      ! 0**0 would arise, I want the corresponding
      ! quantity to be treated as 1E0_wp.
      if (abs(this%c * x) < epsilon(x) .and. n - i + 1 == 0) then
        xterm = 1E0_wp
      else
        xterm = (this%c * x) ** (n - i + 1)
      end if
      if (abs(y) < epsilon(y) .and. m - i + 1 == 0) then
        yterm = 1E0_wp
      else
        yterm = y ** (m - i + 1)
      end if
      deriv = deriv + &
          bin_coeffs(i) * & 
          ( xterm * &
              falling_factorial(real(m, wp) , i - 1) * &
              yterm ) * &
          exp( this%c * x * y )
    end do
    ! Account for the c^m factor
    if (abs(this%c) < epsilon(this%c) .and. m == 0) then
      cmfactor = 1E0_wp
    else
      cmfactor = this%c ** m
    end if
    deriv = cmfactor * deriv
  end function mth_nth_partial
  
  pure function get_c_param(this) result(c_param)
    ! Argument
    class(kernel_expprod), intent(in) :: this
    ! Function result
    real(wp)                          :: c_param
    ! Body
    c_param = this%c
  end function get_c_param
  
  pure subroutine set_c_param(this, c_param)
    ! Arguments
    class(kernel_expprod), intent(inout)  :: this
    real(wp), intent(in)                  :: c_param
    ! Body
    this%c = c_param
  end subroutine set_c_param

end module kernel_expprod_mod
