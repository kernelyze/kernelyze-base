! find_rel_optima.f90
!
! Copyright (C) 2015, 2016 by Kernelyze LLC
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
! created on: 2015-02-07
! updated on: 2015-04-19
! updated on: 2016-04-24 (Avoid use of internal procedures 
!             for thread safety, both here and upstream)
! updated on: 2016-04-25 (fixed typo)
!
! Code to find relative optima between grid points.

module find_rel_optima_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use singlevar_func_mod, only : singlevar_func
  use brent_mod, only : fmin
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  
  implicit none
  
  ! Extend singlevar_func to obtain a derived type
  ! whose eval type-bound procedure returns the negative
  ! of the absolute value of another singlevar_func
  type, private, extends(singlevar_func)  :: neg_abs_func_obj
    ! The function whose negative abs val is to be 
    ! returned
    class(singlevar_func), private, allocatable :: base_func
  contains
    ! Initialize objects of type neg_abs_func_obj
    procedure, pass :: init_neg_abs
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of singlevar_func
    procedure, pass :: eval
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => neg_abs_assign 
  end type neg_abs_func_obj

  contains
  
  ! Implement the deferred procedure for non-elemental
  ! (but still pure) evaluation of singlevar_func for the
  ! neg_abs_func_obj derived type.
  pure recursive function eval(this, arg_pt) result(eval_res)
    ! Arguments
    class(neg_abs_func_obj), intent(in) :: this
    real(wp), intent(in)                :: arg_pt
    ! Function result
    real(wp)                            :: eval_res
    ! Body
    eval_res = -abs(this%base_func%eval(arg_pt))
  end function eval
  
  ! Implement the deferred assignment procedure for
  ! neg_abs_func_obj objects
  subroutine neg_abs_assign(left, right)
    ! Arguments
    class(neg_abs_func_obj), intent(out)  :: left
    class(singlevar_func), intent(in)     :: right
    ! Body
    select type (right)
    type is (neg_abs_func_obj)
      if (allocated(left%base_func)) then
        deallocate(left%base_func)
      end if
      allocate(left%base_func, source = right%base_func)
    end select
  end subroutine neg_abs_assign
  
  ! Initialize a neg_abs_func_obj object
  subroutine init_neg_abs( &
      this, &
      func_to_use)
    ! Arguments
    class(neg_abs_func_obj), intent(out)  :: this
    class(singlevar_func), intent(in)     :: func_to_use
    ! Body
    if (allocated(this%base_func)) then
      deallocate(this%base_func)
    end if
    allocate(this%base_func, source = func_to_use)
  end subroutine init_neg_abs
  
  ! Subroutine to compute relative optima
  subroutine find_rel_optima( &
      result_optima, &
      func_obj, &
      grid, &
      toler, &
      err_stat, &
      err_msg)
    ! Arguments
    real(wp), intent(out), dimension(:, :), allocatable :: result_optima
    class(singlevar_func), intent(in)                   :: func_obj
    real(wp), intent(in), dimension(:)                  :: grid
    real(wp), intent(in)                                :: toler
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables
    type(neg_abs_func_obj)  :: neg_abs_func
    real(wp)                :: xopt, fopt
    integer                 :: j, n_grid_m1
    ! Local variables to handle any allocation problems:
    character(*), parameter         :: proc_name = "Find relative optima: "
    character(len=alloc_errmsg_len) :: alloc_err_msg
    integer                         :: alloc_stat
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! Set up the result grid
    n_grid_m1 = size(grid,1) - 1
    allocate (result_optima(n_grid_m1, 2), stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Since the allocation problem was with result_optima, the intent(out)
      ! variable, I cannot set result_optima to NaN as I would like to.
      ! Simply return.
      return
    end if
    ! Check characteristics of the arguments
    if (n_grid_m1 < 1) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = -1
      end if
      if (present(err_msg)) then
        err_msg = proc_name // &
            'The grid must have at least two elements'
      end if
      ! Set result_optima to NaN
      result_optima = ieee_value(result_optima, ieee_quiet_nan)
      ! Return
      return
    end if
    if (toler <= 0E0_wp) then
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = -1
      end if
      if (present(err_msg)) then
        err_msg = proc_name // &
            'Tolerance must be a positive number'
      end if
      ! Set result_optima to NaN
      result_optima = ieee_value(result_optima, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Initialize neg_abs_func
    call neg_abs_func%init_neg_abs(func_obj)
    do j = 1,n_grid_m1
      ! Minimizing the negative of the absolute value of the function
      ! func will find the point giving the largest absolute value
      ! of func on the (closed) interval from grid(j) to grid(j + 1)
      call fmin(grid(j), grid(j + 1), neg_abs_func, toler, &
          xopt, fopt, local_err_stat, local_err_msg)
      if (local_err_stat /= 0) then
        ! This means that finding the minimum did not go well.
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = local_err_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // local_err_msg
        end if
        ! Set result_optima to NaN
        result_optima = ieee_value(result_optima, ieee_quiet_nan)
        ! Return
        return
      end if
      if (neg_abs_func%eval(grid(j)) < fopt) then
        ! Getting here indicates numerical issues in the fmin call
        ! above, since in this case the left endpoint of the 
        ! interval delivers a value of func that is more extreme
        ! in absolute value than the purportedly most extreme
        ! value on the interval, fopt
        xopt = grid(j)
        fopt = neg_abs_func%eval(grid(j))
      else if (neg_abs_func%eval(grid(j + 1)) < fopt) then
        ! Getting here indicates numerical issues in the fmin call
        ! above, since in this case the right endpoint of the 
        ! interval delivers a value of func that is more extreme
        ! in absolute value than the purportedly most extreme
        ! value on the interval, fopt
        xopt = grid(j + 1)
        fopt = neg_abs_func%eval(grid(j + 1))
      end if
      result_optima(j, 1) = xopt
      result_optima(j, 2) = fopt
    end do
  end subroutine find_rel_optima
end module find_rel_optima_mod