! kernel.f90
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
! updated on: 2016-01-10 (fix to finite-diff for 2nd deriv w. r. t. y)
! updated on: 2016-01-10 (first partial finite diff -> centered diff)
! updated on: 2016-02-07 (added (m, n) partial and finite diffs)
! updated on: 2016-07-01 (let kernels know their [rectangular] domains, 
!             add ability to generate a grid of Chebyshev points over
!             either the x or the y interval defining the domain)
! updated on: 2016-07-05 (small change to chebpts type bound procedure
!             in effort to avoid threading problems)
! updated on: 2016-09-29 (fixed stupidity in default mth_nth_partial
!             implementation -- not noticed earlier because the
!             default implementation was overridden in the frequently-
!             used kernels until now)
! updated on: 2016-11-12 (add ability to evaluate an integral of
!             a kernel times a given singlevar_func)
! updated on: 2016-04-27 (more flexibility for finite differences)
!
! This module contains an abstract derived type representing
! a kernel (a function of two scalar variables):
! $K \left( x , y \right)$.
  
module kernel_mod

  use set_precision, only : wp
  use binomial_coeffs_mod, only : binomial_coeffs
  use chebyshev_points_mod, only : chebyshev_points
  use singlevar_func_mod, only : singlevar_func

  implicit none

  private
  
  ! Increment for finite-difference calculations
  real(wp), parameter, private  :: default_delta = 1E-4_wp
  
  ! Number of points to use in quadrature
  integer, parameter, private   :: n_quad_pts = 200

  type, public, abstract :: kernel
    ! The endpoints of the intervals defining the
    ! rectangular domain of the kernel; these default
    ! to create the standard square [-1 , 1] \cross [-1 , 1].
    ! Note that these are defaults and do *not* implicitly
    ! save these variables.
    real(wp), private :: x_lb = -1E0_wp
    real(wp), private :: x_ub = 1E0_wp
    real(wp), private :: y_lb = -1E0_wp
    real(wp), private :: y_ub = 1E0_wp
  contains
    ! Elemental evaluation procedure for both convenience
    ! and parallelization efficiency.
    procedure, pass                             :: eval_elt
    ! Non-elemental evaluation procedure; this is pure
    ! but not elemental, since the standard states that
    ! an elemental procedure cannot be passed as an
    ! actual argument and there are cases (e. g., for
    ! Borsuk lower bounds or kernel Remez) in which I may 
    ! wish to pass a kernel's evaluation procedure as an 
    ! actual argument.
    procedure(kernel_func_eval), pass, deferred :: eval
    ! Overridable procedures to compute first and second partial
    ! derivatives.  The implementation here is finite-difference-based
    ! but overrides could provide analytical derivatives in types
    ! that extend this abstract type (or, if an extending type is
    ! non-differentiable, could set derivatives to a quiet NaN).
    procedure, pass                             :: first_partial_dx
    procedure, pass                             :: first_partial_dy
    procedure, pass                             :: second_partial_dx_dx
    procedure, pass                             :: second_partial_dy_dy
    procedure, pass                             :: second_partial_dx_dy
    ! Procedures that are explicitly finite-difference approximations
    ! to first and second derivatives
    procedure, pass                             :: finite_diff_dx
    procedure, pass                             :: finite_diff_dy
    procedure, pass                             :: finite_diff_dx_dx
    procedure, pass                             :: finite_diff_dy_dy
    procedure, pass                             :: finite_diff_dx_dy
    ! Higher-order partial derivatives; if m == 0 this is just a
    ! partial w. r. t. y and if n == 0 this is just a partial w. r. t.
    ! x.  As for first_partial and second_partial routines, if an
    ! extending type is non-differentiable it could set derivatives
    ! to a quiet NaN.  The implementation here uses the specific first 
    ! or second partial if m + n <= 2, and otherwise is finite-difference-
    ! based, but types that extend this type could provide analytical
    ! derivatives for m + n > 2.
    procedure, pass                             :: mth_nth_partial
    ! Higher-order finite-difference approximations.  The implementation
    ! does *not* call the specific first- or second-order finite-difference
    ! procedures if m + n <= 2.
    procedure, pass                             :: mth_nth_finite_diff
    ! Procedures that package up the first and second partial
    ! derivatives in ways that users may find convenient
    procedure, pass                             :: gradient
    procedure, pass                             :: hessian
    ! Evaluate the integral of this kernel times a given singlevar_func;
    ! overridable procedure, the implementation here is quadrature
    ! but extending types may override to provide analytical
    ! integrals.
    procedure, pass                             :: eval_integral_svar
    ! A quadrature (numerical integral) of the product of this kernel
    ! and a singlevar_func
    procedure, pass                             :: eval_quadrature_svar
    ! Inspectors and mutators for the interval endpoints defining the
    ! rectangular domain of the kernel
    procedure, pass                             :: get_x_lb
    procedure, pass                             :: get_x_ub
    procedure, pass                             :: get_y_lb
    procedure, pass                             :: get_y_ub
    procedure, pass                             :: set_x_lb
    procedure, pass                             :: set_x_ub
    procedure, pass                             :: set_y_lb
    procedure, pass                             :: set_y_ub
    ! Generate a set of Chebyshev points on either the x interval 
    ! of the domain or the y interval of the domain
    procedure, pass                             :: cheb_pts
  end type kernel

  abstract interface
    ! The key property of a kernel: it can be evaluated
    ! at a pair of scalars.
    pure recursive function kernel_func_eval(this, x, y) result(eval_res)
      ! Bring in the working precision
      import wp
      ! Bring in the abstract derived type declared here
      import kernel
      ! Arguments
      class(kernel), intent(in) :: this
      real(wp), intent(in)      :: x
      real(wp), intent(in)      :: y
      ! Function result
      real(wp)                  :: eval_res
    end function kernel_func_eval
  end interface

  contains
  
  elemental function eval_elt(this, x, y) result(eval_res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    ! Function result
    real(wp)                  :: eval_res
    ! Body -- just invoke the pure but
    ! not elemental procedure.
    eval_res = this%eval(x , y)
  end function eval_elt
  
  elemental function first_partial_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dx(x, y, d)
    else 
      deriv = this%finite_diff_dx(x, y)
    end if
  end function first_partial_dx
  
  elemental function first_partial_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dy(x, y, d)
    else
      deriv = this%finite_diff_dy(x, y)
    end if
  end function first_partial_dy
  
  elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dx_dx(x, y, d)
    else
      deriv = this%finite_diff_dx_dx(x, y)
    end if
  end function second_partial_dx_dx
  
  elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dy_dy(x, y, d)
    else
      deriv = this%finite_diff_dy_dy(x, y)
    end if
  end function second_partial_dy_dy
  
  elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dx_dy(x, y, d)
    else
      deriv = this%finite_diff_dx_dy(x, y)
    end if
  end function second_partial_dx_dy
  
  ! Why pure and not elemental?  In overrides provided by
  ! extending types, I will typically need to use m and / or
  ! n to dimension arrays, etc.
  pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
    ! Arguments
    class(kernel), intent(in) :: this
    integer, intent(in)       :: m
    integer, intent(in)       :: n
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: deriv
    ! Local variables
    real(wp)                  :: locald
    ! Body
    !
    ! If this is actually just a first or second derivative,
    ! call the appropriate specific function.  Otherwise,
    ! call the corresponding finite-difference routine.
    !
    ! If a finite-difference change was given, use it;
    ! otherwise, use the default delta.
    if (present(d)) then
      locald = d
    else
      locald = default_delta
    endif
    
    if ( m + n <= 2 ) then
      if (m == 0 .and. n == 0 ) then
        deriv = this%eval( x , y )
      else if ( m == 0 .and. n == 1 ) then
        deriv = this%first_partial_dy( x , y , locald )
      else if ( m == 1 .and. n == 0 ) then
        deriv = this%first_partial_dx( x , y , locald )
      else if ( m == 0 .and. n == 2 ) then
        deriv = this%second_partial_dy_dy( x , y , locald )
      else if ( m == 2 .and. n == 0 ) then
        deriv = this%second_partial_dx_dx( x , y , locald )
      else ! in this case, ( m == 1 .and. n == 1 )
        deriv = this%second_partial_dx_dy( x , y , locald )
      end if
    else 
      deriv = this%mth_nth_finite_diff( m, n, x, y, locald )
    end if
  end function mth_nth_partial
  
  elemental function finite_diff_dx(this, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    res = (this%eval( x + delta , y ) - this%eval( x - delta , y )) &
          / (2E0_wp * delta)
  end function finite_diff_dx
  
  elemental function finite_diff_dy(this, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    res = (this%eval( x , y + delta ) - this%eval( x , y - delta )) &
          / (2E0_wp * delta)
  end function finite_diff_dy
  
  elemental function finite_diff_dx_dx(this, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    res = (this%eval( x + delta , y ) &
         + this%eval( x - delta , y ) &
         - 2E0_wp * this%eval( x , y )) &
         / (delta * delta)
  end function finite_diff_dx_dx
  
  elemental function finite_diff_dy_dy(this, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    res = (this%eval( x , y + delta ) &
         + this%eval( x , y - delta ) &
         - 2E0_wp * this%eval( x , y )) &
         / (delta * delta)
  end function finite_diff_dy_dy
  
  elemental function finite_diff_dx_dy(this, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    res = (this%eval( x + delta , y + delta ) &
         - this%eval( x + delta , y - delta ) &
         - this%eval( x - delta , y + delta ) &
         + this%eval( x - delta , y - delta )) &
         / (4E0_wp * delta * delta)
  end function finite_diff_dx_dy
  
  pure function mth_nth_finite_diff(this, m, n, x, y, d) result(res)
    ! Arguments
    class(kernel), intent(in) :: this
    integer, intent(in)       :: m
    integer, intent(in)       :: n
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: res
    ! Local variables
    integer                   :: i
    integer                   :: j
    real(wp)                  :: denom
    real(wp)                  :: binoms_for_x(m + 1)
    real(wp)                  :: binoms_for_y(n + 1)
    real(wp)                  :: fun_lattice(m + 1, n + 1)
    real(wp)                  :: coeff_lattice(m + 1, n + 1)
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    denom = 1E0_wp
    do i = 1, m + n
      denom = denom * delta
    end do
    binoms_for_x = binomial_coeffs(m)
    binoms_for_y = binomial_coeffs(n)
    do i = 1, m + 1
      do j = 1, n + 1
        fun_lattice(i, j) = this%eval( &
            x + ( (real(m, wp) / 2E0_wp) - real(i - 1, wp)) * delta, &
            y + ( (real(n, wp) / 2E0_wp) - real(j - 1, wp)) * delta)
        coeff_lattice(i, j) = binoms_for_x(i) * binoms_for_y(j)
        if (mod(i + j, 2) == 1) then
          coeff_lattice(i, j) = -coeff_lattice(i, j)
        end if
      end do
    end do
    res = sum(sum(fun_lattice * coeff_lattice , 1))
    res = res / denom
  end function mth_nth_finite_diff
  
  pure function gradient(this, x, y, d) result(grad)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: grad(2)
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    grad(1) = this%first_partial_dx( x , y , delta )
    grad(2) = this%first_partial_dy( x , y , delta )
  end function gradient
  
  pure function hessian(this, x, y, d) result(hess)
    ! Arguments
    class(kernel), intent(in) :: this
    real(wp), intent(in)      :: x
    real(wp), intent(in)      :: y
    real(wp), intent(in), optional  :: d
    ! Function result
    real(wp)                  :: hess( 2 , 2 )
    ! Local variables
    real(wp)                  :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    hess( 1 , 1 ) = this%second_partial_dx_dx( x , y , delta )
    hess( 1 , 2 ) = this%second_partial_dx_dy( x , y , delta )
    hess( 2 , 1 ) = hess( 1 , 2 )
    hess( 2 , 2 ) = this%second_partial_dy_dy( x , y , delta )
  end function hessian
  
  pure recursive function eval_integral_svar( &
      this, &
      int_is_over_x, &
      arg_pt, &
      svar_to_use ) result(res)
    ! Arguments
    class(kernel), intent(in)         :: this
    logical, intent(in)               :: int_is_over_x
    real(wp), intent(in)              :: arg_pt
    class(singlevar_func), intent(in) :: svar_to_use
    ! Result
    real(wp)                          :: res
    ! Body
    res = this%eval_quadrature_svar( &
        int_is_over_x, &
        arg_pt, &
        svar_to_use )
  end function eval_integral_svar
      
  pure recursive function eval_quadrature_svar( &
      this, &
      int_is_over_x, &
      arg_pt, &
      svar_to_use ) result(res)
    ! Arguments
    class(kernel), intent(in)         :: this
    logical, intent(in)               :: int_is_over_x
    real(wp), intent(in)              :: arg_pt
    class(singlevar_func), intent(in) :: svar_to_use
    ! Result
    real(wp)                          :: res
    ! Local variables
    integer                           :: i
    real(wp)                          :: lower_bd
    real(wp)                          :: upper_bd
    real(wp)                          :: quadpts(n_quad_pts)
    ! Body
    if (int_is_over_x) then
      lower_bd = this%x_lb
      upper_bd = this%x_ub
    else
      lower_bd = this%y_lb
      upper_bd = this%y_ub
    end if
    do i = 1, n_quad_pts
      quadpts(i) = lower_bd + &
          ( real(i - 1, wp) / real(n_quad_pts - 1, wp) ) * &
          ( upper_bd - lower_bd )
    end do
    ! Use elemental evaluation -- this is the most basic
    ! quadrature imaginable
    if (int_is_over_x) then
      res = sum( this%eval_elt(quadpts , arg_pt) * &
                 svar_to_use%eval_elt(quadpts) ) / &
          real(n_quad_pts, wp)
    else
      res = sum( this%eval_elt(arg_pt , quadpts) * &
                svar_to_use%eval_elt(quadpts) ) / &
          real(n_quad_pts, wp)
    end if
  end function eval_quadrature_svar
  
  pure function get_x_lb(this) result(pt)
    ! Arguments
    class(kernel), intent(in) :: this
    ! Function result
    real(wp)                  :: pt
    ! Body
    pt = this%x_lb
  end function get_x_lb
  
  pure function get_x_ub(this) result(pt)
    ! Arguments
    class(kernel), intent(in) :: this
    ! Function result
    real(wp)                  :: pt
    ! Body
    pt = this%x_ub
  end function get_x_ub
  
  pure function get_y_lb(this) result(pt)
    ! Arguments
    class(kernel), intent(in) :: this
    ! Function result
    real(wp)                  :: pt
    ! Body
    pt = this%y_lb
  end function get_y_lb
  
  pure function get_y_ub(this) result(pt)
    ! Arguments
    class(kernel), intent(in) :: this
    ! Function result
    real(wp)                  :: pt
    ! Body
    pt = this%y_ub
  end function get_y_ub
  
  pure subroutine set_x_lb(this, newpt)
    ! Arguments
    class(kernel), intent(inout)  :: this
    real(wp), intent(in)          :: newpt
    ! Body
    this%x_lb = newpt
  end subroutine set_x_lb
  
  pure subroutine set_x_ub(this, newpt)
    ! Arguments
    class(kernel), intent(inout)  :: this
    real(wp), intent(in)          :: newpt
    ! Body
    this%x_ub = newpt
  end subroutine set_x_ub
  
  pure subroutine set_y_lb(this, newpt)
    ! Arguments
    class(kernel), intent(inout)  :: this
    real(wp), intent(in)          :: newpt
    ! Body
    this%y_lb = newpt
  end subroutine set_y_lb
  
  pure subroutine set_y_ub(this, newpt)
    ! Arguments
    class(kernel), intent(inout)  :: this
    real(wp), intent(in)          :: newpt
    ! Body
    this%y_ub = newpt
  end subroutine set_y_ub
  
  pure function cheb_pts(this, over_x, n) result(pts)
    ! Arguments
    class(kernel), intent(in) :: this
    logical, intent(in)       :: over_x
    integer, intent(in)       :: n
    ! Function result
    real(wp)                  :: pts(n)
    ! Body
    pts = chebyshev_points(n)
    ! Now map these points from [-1 , 1] to the x or y
    ! interval
    if (over_x) then
      ! Note that the 1/2 in the slope and intercept of
      ! the linear mapping has been factored out (in front)
      pts = 0.5E0_wp * (  (this%x_ub - this%x_lb) * pts &
                        + (this%x_ub + this%x_lb) )
    else
      ! Note that the 1/2 in the slope and intercept of
      ! the linear mapping has been factored out (in front)
      pts = 0.5E0_wp * (  (this%y_ub - this%y_lb) * pts &
                        + (this%y_ub + this%y_lb) )
    end if  
  end function cheb_pts

end module kernel_mod