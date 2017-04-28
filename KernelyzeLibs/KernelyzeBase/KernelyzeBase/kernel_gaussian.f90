! kernel_gaussian.f90
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
! updated on: 2016-01-08
! updated on: 2016-01-10 (fixed factor of two in partials)
! updated on: 2016-02-09 (added higher-order partials)
! updated on: 2016-10-03 (added scaling to represent p. d. f.)
! updated on: 2016-11-12 (override the eval_integral_svar
!             type-bound procedure of the kernel parent type
!             to provide analytical integrals for certain
!             types of singlevar_funcs)
! updated on: 2017-04-23 (remove map-maker and unnecessary
!             use statements)
! updated on: 2017-04-27 (interfaces adjusted for greater
!             flexibility for other kernels without analytical partials)
!
! This module contains a derived type that extends the 
! kernel_map_maker abstract derived type and makes it concrete: 
! the extension here is to a Gaussian kernel
! of the form $ \frac{ 1 }{ \pi } \sqrt{c} 
! \exp{ \left( - c \left( x - y \right)^2 \right) } $. Because
! of the normalization, this is a probability density function
! over x for fixed y (or over y for fixed x), the normal density
! with variance $ \frac{ 1 }{ 2 c }$.
  
module kernel_gaussian_mod

  use set_precision, only : wp
  use constants_mod, only : inv_sqrtpi
  use hermite_polynomial_mod, only : hermite_polynomial
  use kernel_summable_mod, only : kernel_summable

  implicit none

  private

  type, public, extends(kernel_summable) :: kernel_gaussian
    ! The key constant "c" in the definition of the
    ! Gaussian kernel
    real(wp), private :: c = 1E0_wp ! Default value
  contains
    ! The evaluation procedure for a Gaussian
    ! product kernel
    procedure, pass :: eval => eval_gaussian_kernel
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
  end type kernel_gaussian
  
  contains
  
  pure recursive function eval_gaussian_kernel(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_gaussian), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    ! Function result
    real(wp)                            :: eval_res
    eval_res = inv_sqrtpi * sqrt(this%c) * exp(-this%c * (x - y) * (x - y))
  end function eval_gaussian_kernel
  
  elemental function first_partial_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = -2E0_wp * inv_sqrtpi * sqrt(this%c) * &
        this%c * (x - y) * exp(-this%c * (x - y) * (x - y))
  end function first_partial_dx
  
  elemental function first_partial_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = 2E0_wp * inv_sqrtpi * sqrt(this%c) * &
        this%c * (x - y) * exp(-this%c * (x - y) * (x - y))
  end function first_partial_dy
  
  elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = 2E0_wp * inv_sqrtpi * sqrt(this%c) * &
        this%c * (2E0_wp * this%c * (x - y) * (x - y) - 1E0_wp) &
        * exp(-this%c * (x - y) * (x - y))
  end function second_partial_dx_dx
  
  elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = this%second_partial_dx_dx(x , y)
  end function second_partial_dy_dy
  
  elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    deriv = -this%second_partial_dx_dx(x , y)
  end function second_partial_dx_dy
  
  pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
    ! Arguments
    class(kernel_gaussian), intent(in)  :: this
    integer, intent(in)                 :: m
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Local variables
    real(wp)                            :: hermite_poly
    ! Body
    ! This expression is a consequence of the Rodrigues formula for the
    ! Hermite polynomials (where I use the ``physicist's'' version),
    ! H_{n} \left( x \right) = 
    ! (-1)^{n} \exp{ \left( x^{2} \right)} 
    ! \frac{ d^{n} }{ d x^{n} } \exp{ \left( x^{2} \right)}
    ! along with the scaling formula for Hermite polynomials, as in the
    ! connectiion between the physicist's and probabilist's versions
    ! of the Hermite polynomials, to take care of the c factor.  
    ! Note that the (-1)^{n} becomes (-1)^{m + n} but then gets modified 
    ! due to the chain rule to become (-1)^{m}.
    hermite_poly = hermite_polynomial( m + n , sqrt(this%c) * ( x - y ) )
    deriv = inv_sqrtpi * sqrt(this%c) * &
        ( (-1E0_wp) ** m ) * ( sqrt(this%c) ** (m + n) ) &
        * hermite_poly * exp(-this%c * (x - y) * (x - y))
  end function mth_nth_partial
  
  pure function get_c_param(this) result(c_param)
    ! Argument
    class(kernel_gaussian), intent(in) :: this
    ! Function result
    real(wp)                          :: c_param
    ! Body
    c_param = this%c
  end function get_c_param
  
  pure subroutine set_c_param(this, c_param)
    ! Arguments
    class(kernel_gaussian), intent(inout)  :: this
    real(wp), intent(in)                  :: c_param
    ! Body
    this%c = c_param
  end subroutine set_c_param

end module kernel_gaussian_mod