! kernel_proc_ptr.f90
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
! created on: 2016-08-10
! updated on: 2017-04-27 (interfaces adjusted for greater flexibility
!             for kernels without analytical partials)
!
! Sometimes I want to pass something defined
! as a plain bivariate function as an argument
! to a procedure that takes a kernel.
! This module provides a derived type that
! extends kernel by implementing
! evaluation using a procedure pointer and
! (optionally) derivatives using another
! procedure pointer, enabling this use.
  
module kernel_proc_ptr_mod

  use set_precision, only : wp
  use kernel_summable_mod, only : kernel_summable
  use, intrinsic  :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  type, public, extends(kernel_summable) :: kernel_proc_ptr
    procedure(bivar_func), pointer, private, nopass   :: my_funptr => null()
    procedure(bivar_mnder), pointer, private, nopass  :: my_derivptr => null()
  contains
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of kernel
    procedure, pass :: eval
    ! For now, these are evaluated using the derivative procedure pointer
    ! if that procedure pointer is associated; otherwise, it just uses
    ! the finite difference.  Could modify to have a separate derivative
    ! procedure pointer for each partial to add flexibility at the cost of
    ! more code.
    procedure, pass :: first_partial_dx
    procedure, pass :: first_partial_dy
    procedure, pass :: second_partial_dx_dx
    procedure, pass :: second_partial_dy_dy
    procedure, pass :: second_partial_dx_dy
    ! If a derivative procedure pointer is given, use it in the nth_deriv
    ! procedure; if not, just use the nth finite difference
    procedure, pass :: mth_nth_partial
    ! Implement an assignment procedure
    procedure       :: kernel_proc_ptr_assign
    ! Set the procedure pointer
    procedure, pass :: set_proc_ptr
    ! Get the procedure pointer
    procedure, pass :: get_proc_ptr
    ! Set the derivative procedure pointer
    procedure, pass :: set_deriv_ptr
    ! Get the derivative procedure pointer
    procedure, pass :: get_deriv_ptr
  end type kernel_proc_ptr

  abstract interface
    pure function bivar_func(x, y) result(res)
      ! Bring in the working precision
      import wp
      ! Argument
      real(wp), intent(in)  :: x
      real(wp), intent(in)  :: y
      ! Result
      real(wp)              :: res
    end function bivar_func
    ! mth deriv w. r. t. x of nth deriv w. r. t. y,
    ! so m == 0 means just a y partial and n == 0
    ! means just an x partial, and m == n == 0 would
    ! be just the value.
    pure function bivar_mnder(m, n , x, y) result(deriv)
      ! Bring in the working precision
      import wp
      ! Arguments
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      real(wp), intent(in)  :: x
      real(wp), intent(in)  :: y
      ! Result
      real(wp)              :: deriv
    end function bivar_mnder
  end interface

  contains
  
  pure recursive function eval(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    ! Function result
    real(wp)                            :: eval_res
    ! Body
    if (associated(this%my_funptr)) then
      eval_res = this%my_funptr(x, y)
    else
      eval_res = ieee_value(eval_res, ieee_quiet_nan)
    end if
  end function eval
  
  elemental function first_partial_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(1, 0, x, y)
    else
      if (present(d)) then
        deriv = this%finite_diff_dx(x, y, d)
      else
        deriv = this%finite_diff_dx(x, y)
      end if
    end if
  end function first_partial_dx
  
  elemental function first_partial_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(0, 1, x, y)
    else
      if (present(d)) then
        deriv = this%finite_diff_dy(x, y, d)
      else
        deriv = this%finite_diff_dy(x, y)
      end if
    end if
  end function first_partial_dy
  
  elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(2, 0, x, y)
    else
      if (present(d)) then
        deriv = this%finite_diff_dx_dx(x, y, d)
      else
        deriv = this%finite_diff_dx_dx(x, y)
      end if
    end if
  end function second_partial_dx_dx
  
  elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(0, 2, x, y)
    else
      if (present(d)) then
        deriv = this%finite_diff_dy_dy(x, y, d)
      else
        deriv = this%finite_diff_dy_dy(x, y)
      end if
    end if
  end function second_partial_dy_dy
  
  elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(1, 1, x, y)
    else
      if (present(d)) then
        deriv = this%finite_diff_dx_dy(x, y, d)
      else
        deriv = this%finite_diff_dx_dy(x, y)
      end if
    end if
  end function second_partial_dx_dy
  
  pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
    ! Arguments
    class(kernel_proc_ptr), intent(in)  :: this
    integer, intent(in)                 :: m
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    real(wp), intent(in), optional :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (associated(this%my_derivptr)) then
      deriv = this%my_derivptr(m, n, x, y)
    else
      if (present(d)) then
        deriv = this%mth_nth_finite_diff(m, n, x, y, d)
      else
        deriv = this%mth_nth_finite_diff(m, n, x, y)
      end if
    end if
  end function mth_nth_partial
  
  ! Implement an assignment procedure for kernel_proc_ptr objects
  subroutine kernel_proc_ptr_assign(left, right)
    ! Arguments
    class(kernel_proc_ptr), intent(out) :: left
    class(kernel_summable), intent(in)  :: right
    ! Body
    select type (right)
    type is (kernel_proc_ptr)
      left%my_funptr    => right%my_funptr
      left%my_derivptr  => right%my_derivptr
    end select
  end subroutine kernel_proc_ptr_assign
  
  ! Set the procedure pointer
  subroutine set_proc_ptr(this, new_funptr)
    ! Arguments
    class(kernel_proc_ptr), intent(inout) :: this
    procedure(bivar_func), pointer        :: new_funptr
    ! Body
    this%my_funptr => new_funptr
  end subroutine set_proc_ptr
  
  ! Get the procedure pointer
  subroutine get_proc_ptr(this, this_funptr)
    ! Argument
    class(kernel_proc_ptr), intent(inout) :: this
    ! Result
    procedure(bivar_func), pointer        :: this_funptr
    ! Body
    this_funptr => this%my_funptr
  end subroutine get_proc_ptr
  
  ! Set the derivative procedure pointer
  subroutine set_deriv_ptr(this, new_derptr)
    ! Arguments
    class(kernel_proc_ptr), intent(inout) :: this
    procedure(bivar_mnder), pointer       :: new_derptr
    ! Body
    this%my_derivptr => new_derptr
  end subroutine set_deriv_ptr
  
  ! Get the derivative procedure pointer
  subroutine get_deriv_ptr(this, this_derptr)
    ! Argument
    class(kernel_proc_ptr), intent(inout) :: this
    ! Result
    procedure(bivar_mnder), pointer       :: this_derptr
    ! Body
    this_derptr => this%my_derivptr
  end subroutine get_deriv_ptr

end module kernel_proc_ptr_mod
