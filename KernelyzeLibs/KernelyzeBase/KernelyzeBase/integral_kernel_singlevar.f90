! integral_kernel_singlevar.f90
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
! created on: 2016-08-17
! updated on: 2016-08-17 (renamed and recoded for generality)
! updated on: 2016-08-21 (small fixes)
! updated on: 2016-10-30 (conform to new calculate_integral_kernel_svar)
! updated on: 2016-11-05 (try using calculate_integral_kernel_svar_impl)
! updated on: 2016-11-06 (note that the try on 2016-11-05 did not work,
!             I reverted it)
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
! updated on: 2017-04-23 (purity, elementality, and recursive
!             nature of eval)
!
! This module contains a derived type extending the abstract 
! derived type integral_kernel_measure. This extension restricts 
! to the case of a "measure" that is a singlevar_func.
  
module integral_kernel_singlevar_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use singlevar_func_mod, only : singlevar_func
  use kernel_mod, only : kernel
  use integral_kernel_measure_mod, only : integral_kernel_measure
  
  implicit none

  private

  type, public, extends(integral_kernel_measure) :: integral_kernel_singlevar
    class(singlevar_func), allocatable, private           :: svar_func
  contains
    ! Implement the deferred procedure of the
    ! parent abstract derived type
    procedure, pass :: eval => eval_kernel_svar
    ! How would I override the derivative procedures
    ! of singlevar_func?  I would need to build a new
    ! derived type extending kernel (say, kernel_deriv)
    ! that evaluates to a specified (mth, nth) deriv
    ! of a kernel given arguments (x, y) and then use
    ! this as the kernel to integrate against the singlevar.
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => kernel_svar_assign
    ! Inspectors and mutators
    procedure, pass :: get_singlevar
    procedure, pass :: set_singlevar
  end type integral_kernel_singlevar

contains
  
  pure recursive function eval_kernel_svar(this, arg_pt) result(eval_res)
    ! Arguments
    class(integral_kernel_singlevar), intent(in)  :: this
    real(wp), intent(in)                          :: arg_pt
    ! Function result
    real(wp)                                      :: eval_res
    ! Body
    eval_res = this%integrate_kernel_singlevar( &
        arg_pt, &
        this%svar_func )
  end function eval_kernel_svar
  
  ! Implement the deferred assignment operation
  subroutine kernel_svar_assign(left, right)
    ! Arguments
    class(integral_kernel_singlevar), intent(out) :: left
    class(singlevar_func), intent(in)             :: right
    ! Local variables
    class(kernel), allocatable                    :: kernel_used
    class(singlevar_func), allocatable            :: singlevar_used
    ! Body
    select type (right)
    type is (integral_kernel_singlevar)
      call right%get_singlevar(singlevar_used)
      call left%set_singlevar(singlevar_used)
      call left%set_description(right%get_description())
      call right%get_kernel(kernel_used)
      call left%set_kernel(kernel_used)
    end select
  end subroutine kernel_svar_assign
  
  ! Inspectors and mutators
  
  subroutine get_singlevar(this, curr_svar, err_stat, err_msg)
    ! Arguments
    class(integral_kernel_singlevar), intent(in)      :: this
    class(singlevar_func), allocatable, intent(out)   :: curr_svar
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "get_singlevar procedure of " &
        // "the integral_kernel_singlevar derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Allocate curr_svar, copying into it the correct
    ! dynamic type and type component values from
    ! the svar_func component of this.
    allocate(curr_svar, source = this%svar_func, &
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
  end subroutine get_singlevar
  
  ! Setting the singlevar_func to integrate involves an allocation and 
  ! copy, which preserves: 1) safety and 2) the dynamic type of the input
  ! kernel object.
  subroutine set_singlevar(this, new_svar, err_stat, err_msg)
    ! Arguments
    class(integral_kernel_singlevar), intent(inout) :: this
    class(singlevar_func), allocatable, intent(in)  :: new_svar
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle allocation problems:
    integer                         :: alloc_stat
    character(len=alloc_errmsg_len) :: alloc_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "set_singlevar procedure of " &
        // "the integral_kernel_singlevar derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! If the singlevar_func to integrate is currently allocated,
    ! deallocate it
    if (allocated(this%svar_func)) then
      deallocate(this%svar_func, &
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
    ! Allocate the singlevar_func to integrate, using the
    ! dynamic type and type components of the new_kernel
    allocate(this%svar_func, source = new_svar, &
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
  end subroutine set_singlevar

end module integral_kernel_singlevar_mod