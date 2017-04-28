! integral_kernel_measure.f90
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
! created on: 2015-03-28
! updated on: 2015-06-17
! updated on: 2016-01-09 (derivatives)
! updated on: 2016-02-07 (higher-order derivatives)
! updated on: 2016-08-18 (add support for integrals 
!             using singlevar_funcs)
! updated on: 2016-08-21 (correct typo, add use statement
!             for calculate_integral_kernel_svar)
! updated on: 2016-08-22 (added flag in call to
!             calculate_integral_kernel_svar)
! updated on: 2016-10-27 (conform to new calculate_integral_kernel_svar)
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
! updated on: 2016-11-12 (replace calculate_integral_kernel_svar with
!             a type-bound procedure on kernel)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! This module contains an abstract derived type
! that provides an interface for objects that
! represent either the function of $y$ that results
! from integrating the kernel $K \left( x , y \right)$
! with respect to some signed measure, whether
! discrete or continuous, over $x$ or the function
! of $x$ that results from similar integration
! with respect to some signed measure over $y$.
  
module integral_kernel_measure_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, description_len, err_msg_len
  use kernel_mod, only : kernel
  use singlevar_func_mod, only : singlevar_func
  
  implicit none

  private

  type, public, abstract, extends(singlevar_func) :: integral_kernel_measure
      ! If true, then the measure is over $x$
      ! (otherwise, the measure is over $y$).
      logical, private                        :: integral_is_over_x
      ! A string description of this integral over
      ! one of the variables of a kernel (with respect
      ! to some signed measure)
      character(len=description_len), private :: description
      ! The kernel to integrate; it is allocatable
      ! because that allows it to preserve dynamic
      ! type (to exhibit useful polymorphism) with
      ! respect to child derived types of kernel
      class(kernel), allocatable, private     :: kernel_to_use
    contains
      ! Get the string description of this integral of a kernel
      procedure, pass                             :: get_description
      ! Set the string description of this integral of a kernel
      procedure, pass                             :: set_description
      ! A simple instpector for which variable to integrate over;
      ! returns .true. for integration over $x$, and .false. for
      ! integration over $y$.
      procedure, pass                             :: get_integral_vbl
      ! A simple setter for which variable to integrate over;
      ! set to .true. for integration over $x$, and to .false.
      ! for integration over $y$.
      procedure, pass                             :: set_integral_vbl
      ! The key deferred procedure to evaluate this integral
      procedure(eval_int_kernel), pass, deferred  :: eval
      ! The function to evaluate the underlying kernel 
      procedure, pass                             :: eval_kernel
      ! Evaluate partial derivatives of the underlying kernel
      ! (never need the 2nd cross partials, so it is not provided)
      procedure, pass   :: eval_dkernel_dx
      procedure, pass   :: eval_dkernel_dy
      procedure, pass   :: eval_d2kernel_dx_dx
      procedure, pass   :: eval_d2kernel_dy_dy
      procedure, pass   :: eval_mth_nth_partial_kernel
      ! Compute the integral of the kernel and a singlevar_func
      procedure, pass   :: integrate_kernel_singlevar
      ! Return the kernel object for the kernel to integrate
      procedure, pass                             :: get_kernel
      ! Set the kernel to integrate
      procedure, pass                             :: set_kernel
  end type integral_kernel_measure
   
  abstract interface  
    ! The interface of the function that evaluates
    ! the integral of the kernel with respect to the measure
    ! at a particular point.
    pure function eval_int_kernel(this, arg_pt) result(eval_res)
      ! Bring in the working precision
      import wp
      ! Bring in the abstract derived type defined here
      import integral_kernel_measure
      ! Arguments
      class(integral_kernel_measure), intent(in)  :: this
      real(wp), intent(in)                        :: arg_pt
      ! Function result
      real(wp)                                    :: eval_res
    end function eval_int_kernel
  end interface

contains
  
  ! A simple inspector to get the string description
  ! of this integral of a kernel
  pure function get_description(this) result(desc)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    ! Function result
    character(len=description_len)              :: desc
    ! Body
    desc = this%description
  end function get_description
  ! A simple setter to change the string description
  ! of this integral of a kernel
  pure subroutine set_description(this, new_desc)
    ! Arguments
    class(integral_kernel_measure), intent(inout)   :: this
    ! Note: allow arbitrary-length character
    ! argument.  Longer strings will be truncated
    ! and shorter strings will be padded out, 
    ! as per the Fortran standard for character
    ! assignment.
    character(len=*), intent(in)                    :: new_desc
    ! Body
    this%description = new_desc
  end subroutine set_description
  ! A simple inspector to determine which variable
  ! is being integrated over; returns .true. if
  ! integration is over $x$ and .false. if integration
  ! is over $y$.
  pure logical function get_integral_vbl(this)
    class(integral_kernel_measure), intent(in) :: this
    get_integral_vbl = this%integral_is_over_x
  end function get_integral_vbl
  ! A simple setter to choose which variable to integrate
  ! over; pass in .true. to set integration to be over $x$
  ! and .false. to set integration to be over $y$.
  pure subroutine set_integral_vbl(this, new_integral_dir)
    class(integral_kernel_measure), intent(inout) :: this
    logical, intent(in)                           :: new_integral_dir
    this%integral_is_over_x = new_integral_dir
  end subroutine set_integral_vbl
  ! The pure recursive function that simply evaluates the kernel itself.
  pure recursive function eval_kernel(this, eval_at_x, eval_at_y) result(kernel_res)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: eval_at_x
    real(wp), intent(in)                        :: eval_at_y
    ! Function result
    real(wp)                                    :: kernel_res
    ! Body
    kernel_res = this%kernel_to_use%eval(eval_at_x , eval_at_y)
  end function eval_kernel
  ! First partial of kernel w. r. t. x
  elemental function eval_dkernel_dx( this, x, y, d) result(deriv)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: x
    real(wp), intent(in)                        :: y
    real(wp), intent(in), optional              :: d
    ! Function result
    real(wp)                                    :: deriv
    ! Body
    if (present(d)) then
      deriv = this%kernel_to_use%first_partial_dx(x , y , d)
    else
      deriv = this%kernel_to_use%first_partial_dx(x , y)
    end if
  end function eval_dkernel_dx
  ! First partial of kernel w. r. t. y
  elemental function eval_dkernel_dy( this, x, y, d) result(deriv)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: x
    real(wp), intent(in)                        :: y
    real(wp), intent(in), optional              :: d
    ! Function result
    real(wp)                                    :: deriv
    ! Body
    if (present(d)) then
      deriv = this%kernel_to_use%first_partial_dy(x , y , d)
    else
      deriv = this%kernel_to_use%first_partial_dy(x , y)
    end if
  end function eval_dkernel_dy
  ! Second partial of kernel w. r. t. x
  elemental function eval_d2kernel_dx_dx( this, x, y, d) result(deriv)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: x
    real(wp), intent(in)                        :: y
    real(wp), intent(in), optional              :: d
    ! Function result
    real(wp)                                    :: deriv
    ! Body
    if (present(d)) then
      deriv = this%kernel_to_use%second_partial_dx_dx(x , y , d)
    else
      deriv = this%kernel_to_use%second_partial_dx_dx(x , y)
    end if
  end function eval_d2kernel_dx_dx
  ! Second partial of kernel w. r. t. y
  elemental function eval_d2kernel_dy_dy( this, x, y, d) result(deriv)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: x
    real(wp), intent(in)                        :: y
    real(wp), intent(in), optional              :: d
    ! Function result
    real(wp)                                    :: deriv
    ! Body
    if (present(d)) then
      deriv = this%kernel_to_use%second_partial_dy_dy(x , y , d)
    else
      deriv = this%kernel_to_use%second_partial_dy_dy(x , y)
    end if
  end function eval_d2kernel_dy_dy
  ! The general higher-order mixed partial derivative of the kernel
  pure function eval_mth_nth_partial_kernel( this, m, n, x, y, d ) &
      result(deriv)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    integer, intent(in)                         :: m
    integer, intent(in)                         :: n
    real(wp), intent(in)                        :: x
    real(wp), intent(in)                        :: y
    real(wp), intent(in), optional              :: d
    ! Function result
    real(wp)                                    :: deriv
    ! Body
    if (present(d)) then
      deriv = this%kernel_to_use%mth_nth_partial(m, n, x, y , d)
    else
      deriv = this%kernel_to_use%mth_nth_partial(m, n, x, y)
    end if
  end function eval_mth_nth_partial_kernel
  ! Integrate the kernel and a singlevar_func
  pure function integrate_kernel_singlevar( this, arg_pt, svar ) &
      result(res)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    real(wp), intent(in)                        :: arg_pt
    class(singlevar_func), intent(in)           :: svar
    ! Result
    real(wp)                                    :: res
    ! Body
    res = this%kernel_to_use%eval_integral_svar( &
        this%integral_is_over_x, &
        arg_pt, &
        svar )
  end function integrate_kernel_singlevar
  ! When I get the kernel to integrate, I make a 
  ! copy of it and return that to be sure that the
  ! object returned has the proper dynamic type.
  subroutine get_kernel(this, kernel_used, err_stat, err_msg)
    ! Arguments
    class(integral_kernel_measure), intent(in)  :: this
    class(kernel), allocatable, intent(out)     :: kernel_used
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "get_kernel procedure of " &
        // "the integral_kernel_measure derived type: "
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
    ! the kernel_to_use component of this.
    allocate(kernel_used, source = this%kernel_to_use, &
        stat = alloc_stat, errmsg = alloc_err_msg)
    ! Handle any allocation problem
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
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
    class(integral_kernel_measure), intent(inout) :: this
    class(kernel), intent(in)                     :: new_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "set_kernel procedure of " &
        // "the integral_kernel_measure derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! If the kernel to integrate is currently allocated,
    ! deallocate it
    if (allocated(this%kernel_to_use)) then
      deallocate(this%kernel_to_use, &
          stat = alloc_stat, errmsg = alloc_err_msg)
      ! Handle any allocation problem
      if (alloc_stat /= 0) then
        ! If I am here, there was an allocation problem.
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
    ! Allocate the kernel to integrate, using the
    ! dynamic type and type components of the new_kernel
    allocate(this%kernel_to_use, source = new_kernel, &
        stat = alloc_stat, errmsg = alloc_err_msg)
    ! Handle any allocation problem
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
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

end module integral_kernel_measure_mod