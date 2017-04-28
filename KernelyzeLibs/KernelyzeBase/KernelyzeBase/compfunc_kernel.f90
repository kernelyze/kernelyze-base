! compfunc_kernel.f90
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
! created on: 2015-04-05
! updated on: 2015-06-16
! added assignment on : 2015-11-12
! updated on: 2016-01-09 (derivatives)
! updated on: 2016-02-07 (higher-order derivatives)
! updated on: 2016-08-22 (integral of a compfunc_kernel
!             times a given singlevar_func)
! updated on: 2016-10-27 (conform to new calculate_integral_kernel_svar)
! updated on: 2016-10-30 (rename getter and setter for integral calcer)
! updated on: 2016-11-05 (improve spacing)
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
! updated on: 2016-11-12 (replace calculate_integral_kernel_svar with
!             a type-bound procedure on kernel)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! This module contains a derived type extending
! the singlevar_func abstract derived type. The
! derived type here provides the function obtained
! by fixing one of the arguments to a kernel.  This
! is sometimes referred to as a "section" of a
! kernel
  
module compfunc_kernel_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use singlevar_func_mod, only : singlevar_func
  use kernel_mod, only : kernel
  
  implicit none

  private

  type, public, extends(singlevar_func) :: compfunc_kernel
    ! If true, form the kernel section by fixing the
    ! $x$ argument of the kernel; otherwise, form the
    ! kernel section by fixing the $y$ argument of the
    ! kernel
    logical, private                    :: x_is_fixed
    ! The point at which one of the arguments of the
    ! kernel is evaluated to obtain this compfunc_kernel
    real(wp), private                   :: eval_pt
    ! The kernel whose section is taken.  It must be
    ! allocatable to permit polymorphic behavior.
    class(kernel), allocatable, private :: kernel_sectioned
  contains
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of singlevar_func
    procedure, pass :: eval
    ! Override the singlevar_func derivative procedures
    procedure, pass :: first_deriv
    procedure, pass :: second_deriv
    procedure, pass :: nth_deriv
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => compfunc_assign
    ! Inspector function for x_is_fixed
    procedure, pass :: get_x_is_fixed
    ! Setter procedure for x_is_fixed
    procedure, pass :: set_x_is_fixed
    ! Inspector for eval_pt
    procedure, pass :: get_eval_pt
    ! Setter for eval_pt
    procedure, pass :: set_eval_pt
    ! Inspector for the kernel being sectioned
    procedure, pass :: get_kernel
    ! Setter for the kernel being sectioned
    procedure, pass :: set_kernel
    procedure, pass :: integrate_svar
  end type compfunc_kernel

  contains
  
  ! Implement the deferred procedure for non-elemental
  ! (but still pure) evaluation.
  pure recursive function eval(this, arg_pt) result(eval_res)
    ! Arguments
    class(compfunc_kernel), intent(in)  :: this
    real(wp), intent(in)                :: arg_pt
    ! Function result
    real(wp)                            :: eval_res
    ! Body
    if (this%x_is_fixed) then
      eval_res = this%kernel_sectioned%eval(this%eval_pt , arg_pt)
    else
      eval_res = this%kernel_sectioned%eval(arg_pt , this%eval_pt)
    end if
  end function eval
  ! First derivative, using the kernel being sectioned
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(compfunc_kernel), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in), optional      :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (present(d)) then
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%first_partial_dy(this%eval_pt , x , d)
      else
        deriv = this%kernel_sectioned%first_partial_dx(x , this%eval_pt , d)
      end if
    else
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%first_partial_dy(this%eval_pt , x)
      else
        deriv = this%kernel_sectioned%first_partial_dx(x , this%eval_pt)
      end if
    end if
  end function first_deriv
  ! Second derivative, using the kernel being sectioned
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(compfunc_kernel), intent(in)  :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in), optional      :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (present(d)) then
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%second_partial_dy_dy( &
            this%eval_pt , x , d)
      else
        deriv = this%kernel_sectioned%second_partial_dx_dx( &
            x , this%eval_pt , d)
      end if
    else
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%second_partial_dy_dy( &
            this%eval_pt , x)
      else
        deriv = this%kernel_sectioned%second_partial_dx_dx( &
            x , this%eval_pt)
      end if
    end if
  end function second_deriv
  ! nth derivative, using the kernel being sectioned
  pure function nth_deriv(this, n, x, d) result(deriv)
    ! Arguments
    class(compfunc_kernel), intent(in)  :: this
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: x
    real(wp), intent(in), optional      :: d
    ! Function result
    real(wp)                            :: deriv
    ! Body
    if (present(d)) then
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%mth_nth_partial( &
            0 , n , this%eval_pt , x , d)
      else
        deriv = this%kernel_sectioned%mth_nth_partial( &
            n , 0 , x , this%eval_pt , d)
      end if    
    else
      if (this%x_is_fixed) then
        deriv = this%kernel_sectioned%mth_nth_partial( &
            0 , n , this%eval_pt , x)
      else
        deriv = this%kernel_sectioned%mth_nth_partial( &
            n , 0 , x , this%eval_pt)
      end if    
    end if
  end function nth_deriv
  ! Assignment
  subroutine compfunc_assign(left, right)
    ! Arguments
    class(compfunc_kernel), intent(out) :: left
    class(singlevar_func), intent(in)   :: right
    ! Body
    select type (right)
    type is (compfunc_kernel)
      left%x_is_fixed = right%x_is_fixed
      left%eval_pt    = right%eval_pt
      if (allocated(left%kernel_sectioned)) then
        deallocate(left%kernel_sectioned)
      end if
      allocate(left%kernel_sectioned, source = right%kernel_sectioned)
    end select
  end subroutine compfunc_assign
  ! Inspector for whether x is the fixed variable
  ! (if not, then y is the fixed variable)
  pure function get_x_is_fixed(this) result(is_x_fixed)
    ! Argument
    class(compfunc_kernel), intent(in)  :: this
    ! Function result
    logical                             :: is_x_fixed
    ! Body
    is_x_fixed = this%x_is_fixed
  end function get_x_is_fixed
  ! Setter for whether x is the fixed variable
  pure subroutine set_x_is_fixed(this, fix_x)
    ! Arguments
    class(compfunc_kernel), intent(inout)   :: this
    logical, intent(in)                     :: fix_x
    ! Body
    this%x_is_fixed = fix_x
  end subroutine set_x_is_fixed
  ! Inspector for eval point
  pure function get_eval_pt(this) result(eval_point)
    ! Argument
    class(compfunc_kernel), intent(in)  :: this
    ! Function result
    real(wp)                            :: eval_point
    ! Body
    eval_point = this%eval_pt
  end function get_eval_pt
  ! Setter for eval point
  pure subroutine set_eval_pt(this, eval_point)
    ! Arguments
    class(compfunc_kernel), intent(inout)   :: this
    real(wp), intent(in)                    :: eval_point
    ! Body
    this%eval_pt = eval_point
  end subroutine set_eval_pt
  ! When I get the kernel to section, I make a 
  ! copy of it and return that to be sure that the
  ! object returned has the proper dynamic type.
  subroutine get_kernel(this, kernel_used, err_stat, err_msg)
    ! Arguments
    class(compfunc_kernel), intent(in)      :: this
    class(kernel), allocatable, intent(out) :: kernel_used
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "get_kernel procedure of " &
        // "the compfunc_kernel derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Allocate kernel_used, copying into it the correct
    ! dynamic type and type component values from
    ! the kernel_sectioned component of this.
    allocate(kernel_used, source = this%kernel_sectioned, &
        stat = alloc_stat, errmsg = alloc_err_msg)
    ! Handle any allocation problem
    if (alloc_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      return
    end if
  end subroutine get_kernel
  ! Setting the kernel to integrate involves an allocation and copy, 
  ! which preserves: 1) safety and 2) the dynamic type of the input
  ! kernel object.
  subroutine set_kernel(this, new_kernel, err_stat, err_msg)
    ! Arguments
    class(compfunc_kernel), intent(inout) :: this
    class(kernel), intent(in)             :: new_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "set_kernel procedure of " &
        // "the compfunc_kernel derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! If the kernel being sectioned is currently allocated,
    ! deallocate it
    if (allocated(this%kernel_sectioned)) then
      deallocate(this%kernel_sectioned, &
          stat = alloc_stat, errmsg = alloc_err_msg)
      ! Handle any allocation problem
      if (alloc_stat /= 0) then
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = alloc_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // alloc_err_msg
        end if
        return
      end if
    end if
    ! Allocate the kernel being sectioned, using the
    ! dynamic type and type components of the new_kernel
    allocate(this%kernel_sectioned, source = new_kernel, &
        stat = alloc_stat, errmsg = alloc_err_msg)
    ! Handle any allocation problem
    if (alloc_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      return
    end if
  end subroutine set_kernel
  
  pure function integrate_svar(this, svar) result(res)
    ! Arguments
    class(compfunc_kernel), intent(in)          :: this
    class(singlevar_func), intent(in)           :: svar
    ! Result
    real(wp)                            :: res
    ! Body
    res = this%kernel_sectioned%eval_integral_svar( &
        .not. this%x_is_fixed, &
        this%eval_pt, &
        svar )
  end function integrate_svar

end module compfunc_kernel_mod