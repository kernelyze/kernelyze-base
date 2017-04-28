! singlevar_proc_ptr.f90
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
! created on: 2016-04-25
! updated on: 2016-04-25
! updated on: 2016-08-10
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! Sometimes I want to pass something defined
! as a plain function (e. g., an intrinsic
! function like sin, exp, or log) as an argument
! to a procedure that takes a singlevar_func.
! This module provides a derived type that
! extends singlevar_func by implementing
! evaluation using a procedure pointer and
! (optionally) derivatives using another
! procedure pointer, enabling this use.
  
module singlevar_proc_ptr_mod

  use set_precision, only : wp
  use singlevar_func_mod, only : singlevar_func
  use, intrinsic  :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  type, public, extends(singlevar_func) :: singlevar_proc_ptr
    procedure(univar_func), pointer, private, nopass  :: my_funptr => null()
    procedure(univar_nder), pointer, private, nopass  :: my_derivptr => null()
  contains
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of singlevar_func
    procedure, pass :: eval
    ! For now, this is evaluated using the derivative procedure pointer
    ! if that procedure pointer is associated; otherwise, it just uses
    ! the finite difference.  Could modify to have a separate first
    ! derivative procedure pointer to add flexibility at the cost of
    ! more code.
    procedure, pass :: first_deriv
    ! For now, this is evaluated using the derivative procedure pointer
    ! if that procedure pointer is associated; otherwise, it just uses
    ! the 2nd finite difference.  Could modify to have a separate 2nd
    ! derivative procedure pointer to add flexibility at the cost of
    ! more code.
    procedure, pass :: second_deriv
    ! If a derivative procedure pointer is given, use it in the nth_deriv
    ! procedure; if not, just use the nth finite difference
    procedure, pass :: nth_deriv
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => singlevar_proc_ptr_assign
    ! Set the procedure pointer
    procedure, pass :: set_proc_ptr
    ! Get the procedure pointer
    procedure, pass :: get_proc_ptr
    ! Set the derivative procedure pointer
    procedure, pass :: set_deriv_ptr
    ! Get the derivative procedure pointer
    procedure, pass :: get_deriv_ptr
  end type singlevar_proc_ptr

  abstract interface
    pure function univar_func(x) result(y)
      ! Bring in the working precision
      import wp
      ! Argument
      real(wp), intent(in)  :: x
      ! Result
      real(wp)              :: y
    end function univar_func
    pure function univar_nder(n , x) result(deriv)
      ! Bring in the working precision
      import wp
      ! Arguments
      integer, intent(in)   :: n
      real(wp), intent(in)  :: x
      ! Result
      real(wp)              :: deriv
    end function univar_nder
  end interface

  contains
  
  ! Implement the deferred procedure for non-elemental
  ! (but still pure) evaluation of singlevar_func for the
  ! singlevar_proc_ptr derived type.
  pure recursive function eval(this, arg_pt) result(eval_res)
    ! Arguments
    class(singlevar_proc_ptr), intent(in) :: this
    real(wp), intent(in)                  :: arg_pt
    ! Function result
    real(wp)                              :: eval_res
    ! Body
    if (associated(this%my_funptr)) then
      eval_res = this%my_funptr(arg_pt)
    else
      eval_res = ieee_value(eval_res, ieee_quiet_nan)
    end if
  end function eval
  
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(singlevar_proc_ptr), intent(in) :: this
    real(wp), intent(in)                  :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                              :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr( 1 , x )
    else
      if (present(d)) then
        deriv = this%finite_diff_dx( x , d )    
      else
        deriv = this%finite_diff_dx( x )    
      end if
    end if
  end function first_deriv
  
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(singlevar_proc_ptr), intent(in) :: this
    real(wp), intent(in)                  :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                              :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr( 2 , x )
    else
      if (present(d)) then
        deriv = this%finite_diff_dx( x , d )    
      else
        deriv = this%finite_diff_dx( x )    
      end if
    end if    
  end function second_deriv
  
  pure function nth_deriv(this, n, x, d) result(deriv)
    ! Arguments
    class(singlevar_proc_ptr), intent(in) :: this
    integer, intent(in)                   :: n
    real(wp), intent(in)                  :: x
    real(wp), intent(in), optional    :: d
    ! Function result
    real(wp)                              :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr( n , x )
    else
      if (present(d)) then
        deriv = this%nth_finite_diff( n , x )
      else
        deriv = this%nth_finite_diff( n , x , d )
      end if
    end if
  end function nth_deriv
  
  ! Implement the deferred assignment procedure for
  ! singlevar_proc_ptr objects
  subroutine singlevar_proc_ptr_assign(left, right)
    ! Arguments
    class(singlevar_proc_ptr), intent(out)  :: left
    class(singlevar_func), intent(in)       :: right
    ! Body
    select type (right)
    type is (singlevar_proc_ptr)
      left%my_funptr    => right%my_funptr
      left%my_derivptr  => right%my_derivptr
    end select
  end subroutine singlevar_proc_ptr_assign
  
  ! Set the procedure pointer
  subroutine set_proc_ptr(this, new_funptr)
    ! Arguments
    class(singlevar_proc_ptr), intent(inout)  :: this
    procedure(univar_func), pointer           :: new_funptr
    ! Body
    this%my_funptr => new_funptr
  end subroutine set_proc_ptr
  
  ! Get the procedure pointer
  subroutine get_proc_ptr(this, curr_funptr)
    ! Argument
    class(singlevar_proc_ptr), intent(inout)  :: this
    ! Result
    procedure(univar_func), pointer           :: curr_funptr
    ! Body
    curr_funptr => this%my_funptr
  end subroutine get_proc_ptr
  
  ! Set the derivative procedure pointer
  subroutine set_deriv_ptr(this, new_derptr)
    ! Arguments
    class(singlevar_proc_ptr), intent(inout)  :: this
    procedure(univar_nder), pointer           :: new_derptr
    ! Body
    this%my_derivptr => new_derptr
  end subroutine set_deriv_ptr
  
  ! Get the derivative procedure pointer
  subroutine get_deriv_ptr(this, curr_derptr)
    ! Argument
    class(singlevar_proc_ptr), intent(inout) :: this
    ! Result
    procedure(univar_nder), pointer          :: curr_derptr
    ! Body
    curr_derptr => this%my_derivptr
  end subroutine get_deriv_ptr

end module singlevar_proc_ptr_mod
