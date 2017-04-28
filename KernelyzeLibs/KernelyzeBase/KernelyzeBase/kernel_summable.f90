! kernel_summable.f90
!
! Copyright (c) 2015, 2017 by Kernelyze LLC
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
! created on: 2015-06-19
! updated on: 2015-06-19
! updated on: 2016-09-02 (added kernel_integral_svar type-bound procedure)
! updated on: 2017-04-23 (remove unnecessary use statements)
!
! This module contains an abstract derived type representing
! a kernel (a function of two scalar variables)
! $K \left( x , y \right)$ that is also summable, so that
! it can form integrals over either $x$ or $y$ with respect
! to discrete signed distributions (so that the integration
! is a summation).
  
module kernel_summable_mod

  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use singlevar_func_mod, only : singlevar_func
  use kernel_mod, only : kernel
  use integral_kernel_disc_mod, only : integral_kernel_disc
  use integral_kernel_singlevar_mod, only : integral_kernel_singlevar

  implicit none

  private

  type, public, abstract, extends(kernel) :: kernel_summable
  contains
    ! Construct the integral of this kernel with respect to
    ! a discrete signed measure over either $x$ or $y$
    procedure, pass :: kernel_integral
    ! Build an object that integrates the kernel
    ! over either $x$ or $y$ using a singlevar_func.
    procedure, pass :: kernel_integral_svar
  end type kernel_summable

  contains
  
  ! Build an object that integrates the kernel
  ! over either $x$ or $y$ using a distribution
  ! with discrete support (so that the actual operation
  ! is a sum).
  subroutine kernel_integral( &
      this, &
      sum_is_over_x, &
      weights, &
      eval_pts, &
      integral, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel_summable), intent(in)                      :: this
    logical, intent(in)                                     :: sum_is_over_x
    real(wp), dimension(:), intent(in)                      :: weights
    real(wp), dimension(size(weights)), intent(in)          :: eval_pts
    class(integral_kernel_disc), allocatable, intent(out)   :: integral
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "kernel_integral procedure of " &
        // "the kernel_summable derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Allocate and set up the integral
    allocate(integral, &
        stat = local_err_stat, errmsg = local_err_msg)
    ! Handle any allocation problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    call integral%set_integral_vbl(sum_is_over_x)
    call integral%set_kernel(this , local_err_stat , local_err_msg)
    ! Handle any problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    call integral%set_weights(weights)
    call integral%set_eval_pts(eval_pts)
  end subroutine kernel_integral
      
  ! Build an object that integrates the kernel
  ! over either $x$ or $y$ using a singlevar_func.
  subroutine kernel_integral_svar( &
      this, &
      int_is_over_x, &
      svar, &
      integral, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel_summable), intent(in)                :: this
    logical, intent(in)                               :: int_is_over_x
    class(singlevar_func), allocatable, intent(in)    :: svar
    class(integral_kernel_singlevar), allocatable, &
                                          intent(out) :: integral
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = "kernel_integral_svar " &
        // "procedure of the kernel_summable derived type: "
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Allocate the integral object
    allocate(integral, &
        stat = local_err_stat, errmsg = local_err_msg)
    ! Handle any allocation problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    ! Set the values of the new object that I don't need to compute here
    call integral%set_integral_vbl(int_is_over_x)
    call integral%set_description('Integral vs svar from summable kernel.')
    call integral%set_kernel(this, local_err_stat, local_err_msg)
    ! Handle any problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
    call integral%set_singlevar(svar, local_err_stat, local_err_msg)
    ! Handle any problem
    if (local_err_stat /= 0) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      return
    end if
  end subroutine kernel_integral_svar      

end module kernel_summable_mod