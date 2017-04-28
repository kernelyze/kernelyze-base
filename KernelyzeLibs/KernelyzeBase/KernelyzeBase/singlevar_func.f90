! singlevar_func.f90
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
! created on: 2015-04-04
! updated on: 2015-06-16
! added assignment on: 2015-11-12
! updated on: 2016-01-09 (derivatives, finite-diff)
! updated on: 2016-01-10 (centered difference for finite diff)
! updated on: 2016-02-07 (higher-order partials and finite diffs)
! updated on: 2016-02-11 (fixed sign of binom coeff in higher-order fin diff)
! updated on: 2016-08-23 (added integrate_svar, x_lb, x_ub, and Chebyshev pts)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! This module contains an abstract derived type representing
! a univariate function (suitable to be, e. g., a component in 
! one of the $n$ function pairs used in a rank-$n$ kernel).
  
module singlevar_func_mod

  use set_precision, only : wp
  use binomial_coeffs_mod, only : binomial_coeffs
  use chebyshev_points_mod, only : chebyshev_points

  implicit none

  private
  
  ! Increment for finite-difference calculations
  real(wp), parameter, private  :: default_delta = 1E-4_wp
  
  ! Number of points to use in quadrature
  integer, parameter, private   :: n_quad_pts = 200

  type, public, abstract :: singlevar_func
    real(wp), private :: x_lb
    real(wp), private :: x_ub
  contains
    ! Elemental evaluation procedure for both convenience
    ! and parallelization efficiency.
    procedure, pass                             :: eval_elt
    ! Non-elemental evaluation procedure; this is pure
    ! but not elemental, since the standard states that
    ! an elemental procedure cannot be passed as an
    ! actual argument and there are cases (e. g., for
    ! root-finding or maximization) in which I may wish
    ! to pass a singlevar_func's evaluation procedure
    ! as an actual argument.
    procedure(nonelt_func_eval), pass, deferred :: eval
    ! Overridable procedures to compute first and second derivatives
    ! The implementation here is finite-difference-based but
    ! overrides could provide analytical derivatives in types
    ! that extend this abstract type (or, if an extending type is
    ! non-differentiable, could set derivatives to a quiet NaN).
    procedure, pass                             :: first_deriv
    procedure, pass                             :: second_deriv
    ! Similarly an overridable procedure to compute the nth derivative.
    procedure, pass                             :: nth_deriv
    ! Procedures that are explicitly finite-difference approximations
    ! to first and second derivatives
    procedure, pass                             :: finite_diff_dx
    procedure, pass                             :: finite_diff_dx_dx
    ! Similarly, a procedure to compute the finite-difference
    ! approximation to the nth derivative
    procedure, pass                             :: nth_finite_diff
    ! Integrate the product of this and another singlevar_func;
    ! overridable procedure, the implementation here is quadrature
    ! but extending types may override to provide analytical
    ! integrals.
    procedure, pass                             :: integrate_svar
    ! A quadrature (numerical integral) of the product of this
    ! and another singlevar_func
    procedure, pass                             :: quadrature_svar
    ! Inspector and mutator for the x_lb and x_ub private components
    procedure, pass                             :: get_x_lb
    procedure, pass                             :: get_x_ub
    procedure, pass                             :: set_x_lb
    procedure, pass                             :: set_x_ub
    ! Generate a set of Chebyshev points on the x interval 
    ! of the domain
    procedure, pass                             :: cheb_pts
    ! A procedure and its binding to assignment
    procedure(singlevar_assignment), deferred   :: singlevar_assign
    generic, public :: assignment(=) => singlevar_assign
  end type singlevar_func

  abstract interface
    ! The key property of a component function: it can be evaluated
    ! at a scalar.  This evaluation is pure but not elemental,
    ! see comment above.
    pure recursive function nonelt_func_eval(this, arg_pt) result(eval_res)
      ! Bring in the working precision
      import wp
      ! Bring in the abstract derived type declared here
      import singlevar_func
      ! Arguments
      class(singlevar_func), intent(in) :: this
      real(wp), intent(in)              :: arg_pt
      ! Function result
      real(wp)                          :: eval_res
    end function nonelt_func_eval
    ! This is the interface for the deferred assignment operation
    subroutine singlevar_assignment(left, right)
      ! Bring in the abstract derived type declared here
      import singlevar_func
      ! Arguments
      class(singlevar_func), intent(out)  :: left
      class(singlevar_func), intent(in)   :: right
    end subroutine singlevar_assignment
  end interface

  contains
  
  elemental function eval_elt(this, x) result(eval_res)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    real(wp), intent(in)              :: x
    ! Function result
    real(wp)                          :: eval_res
    ! Body -- just invoke the pure but
    ! not elemental procedure.
    eval_res = this%eval(x)
  end function eval_elt
  
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dx( x , d )    
    else
      deriv = this%finite_diff_dx( x )    
    end if
  end function first_deriv
  
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: deriv
    ! Body
    if (present(d)) then
      deriv = this%finite_diff_dx_dx( x , d )    
    else
      deriv = this%finite_diff_dx_dx( x )    
    end if
  end function second_deriv
  
  pure function nth_deriv(this, n, x, d) result(deriv)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    integer, intent(in)               :: n
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: deriv
    ! Body
    if (present(d)) then
      deriv = this%nth_finite_diff( n , x , d )
    else
      deriv = this%nth_finite_diff( n , x )
    end if
  end function nth_deriv
  
  elemental function finite_diff_dx(this, x, d) result(fin_diff)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: fin_diff
    ! Local variables
    real(wp)                          :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    fin_diff = (this%eval( x + delta ) - this%eval( x - delta )) &
              / (2E0_wp * delta)   
  end function finite_diff_dx
  
  elemental function finite_diff_dx_dx(this, x, d) result(fin_diff)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: fin_diff
    ! Local variables
    real(wp)                          :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    fin_diff = (this%eval( x + delta ) &
         + this%eval( x - delta ) &
         - 2E0_wp * this%eval( x )) &
         / (delta * delta)   
  end function finite_diff_dx_dx
  
  pure function nth_finite_diff(this, n, x, d) result(fin_diff)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    integer, intent(in)               :: n
    real(wp), intent(in)              :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                          :: fin_diff
    ! Local variables
    integer                   :: i
    real(wp)                  :: denom
    real(wp)                  :: binoms(n + 1)
    real(wp)                  :: fun_lattice(n + 1)
    real(wp)                  :: coeff_lattice(n + 1)
    ! Local variables
    real(wp)                          :: delta
    ! Body
    if (present(d)) then
      delta = d
    else
      delta = default_delta
    end if
    denom = 1E0_wp
    do i = 1, n
      denom = denom * delta
    end do
    binoms = binomial_coeffs(n)
    do i = 1, n + 1
      fun_lattice(i) = this%eval( &
          x + ( (real(n, wp) / 2E0_wp) - real(i - 1, wp)) * delta)
      coeff_lattice(i) = binoms(i)
      if (mod(i, 2) == 0) then
        coeff_lattice(i) = -coeff_lattice(i)
      end if
    end do
    fin_diff = sum(fun_lattice * coeff_lattice)
    fin_diff = fin_diff / denom
  end function nth_finite_diff
  
  pure function integrate_svar(this, svar) result(res)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    class(singlevar_func), intent(in) :: svar
    ! Result
    real(wp)                          :: res
    ! Body
    ! Default is to use quadrature -- extending
    ! types may override to provide analytical
    ! integrals.
    res = quadrature_svar(this, svar)
  end function integrate_svar
  
  pure function quadrature_svar(this, svar) result(res)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    class(singlevar_func), intent(in) :: svar
    ! Result
    real(wp)                          :: res
    ! Local variables
    integer                           :: i
    real(wp)                          :: quadpts(n_quad_pts)
    ! Body
    do i = 1, n_quad_pts
      quadpts(i) = this%x_lb + &
          ( real(i - 1, wp) / real(n_quad_pts - 1, wp) ) * &
          ( this%x_ub - this%x_lb )
    end do
    ! Use elemental evaluation -- this is the most basic
    ! quadrature imaginable
    res = sum( this%eval_elt(quadpts) * svar%eval_elt(quadpts) ) / &
        real(n_quad_pts, wp)
  end function quadrature_svar
  
  pure function cheb_pts(this, n_pts) result(cpts)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    integer, intent(in)               :: n_pts
    ! Result
    real(wp)                          :: cpts(n_pts)
    ! Body
    cpts = chebyshev_points(n_pts)
    ! Note that the 1/2 in the slope and intercept of
    ! the linear mapping has been factored out (in front)
    cpts = 0.5E0_wp * (  (this%x_ub - this%x_lb) * cpts &
                       + (this%x_ub + this%x_lb) )
  end function cheb_pts
  
  pure function get_x_lb(this) result(curr_x_lb)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    ! Result
    real(wp)                          :: curr_x_lb
    ! Body
    curr_x_lb = this%x_lb
  end function get_x_lb
  
  pure function get_x_ub(this) result(curr_x_ub)
    ! Arguments
    class(singlevar_func), intent(in) :: this
    ! Result
    real(wp)                          :: curr_x_ub
    ! Body
    curr_x_ub = this%x_ub
  end function get_x_ub
  
  pure subroutine set_x_lb(this, new_x_lb)
    ! Arguments
    class(singlevar_func), intent(inout)  :: this
    real(wp), intent(in)                  :: new_x_lb
    ! Body
    this%x_lb = new_x_lb
  end subroutine set_x_lb
  
  pure subroutine set_x_ub(this, new_x_ub)
    ! Arguments
    class(singlevar_func), intent(inout)  :: this
    real(wp), intent(in)                  :: new_x_ub
    ! Body
    this%x_ub = new_x_ub
  end subroutine set_x_ub

end module singlevar_func_mod
