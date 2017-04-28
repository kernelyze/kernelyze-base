! find_all_zeros.f90
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
! created on: 2015-02-01
! updated on: 2015-04-22
! updated on: 2016-04-24 (Avoid upstream use of internal procedures 
!             for thread safety)
!
! Code to find zeros between grid points.

module find_all_zeros_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use brent_mod, only : zeroin
  use m_unirnk, only : unirnk
  use singlevar_func_mod, only : singlevar_func
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  
  implicit none
  
  contains

  pure subroutine find_all_zeros( &
      result_zeros, &
      func_obj, &
      grid, &
      toler, &
      err_stat, &
      err_msg)
    ! Arguments
    real(wp), intent(out), dimension(:), allocatable  :: result_zeros
    class(singlevar_func), intent(in)                 :: func_obj
    real(wp), intent(in), dimension(:)                :: grid
    real(wp), intent(in)                              :: toler
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables
    integer                             :: j, n_grid_m1, zeros_j, nunique
    real(wp)                            :: toler_sq, fun_at_j, fun_at_j_plus1, prod_fun_evals
    real(wp), dimension(:), allocatable :: temp
    integer, dimension(:), allocatable  :: indices_of_unique
    ! Local variables to handle any allocation problems:
    character(*), parameter         :: proc_name = "Find all zeros: "
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
    ! Set up the proper grid size and allocate the result array
    n_grid_m1 = size(grid,1) - 1
    allocate (result_zeros(n_grid_m1), stat = alloc_stat , errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Since the allocation problem was with result_zeros, the intent(out)
      ! variable, I cannot set result_zeros to NaN as I would like to.
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
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
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
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
      ! Return
      return
    end if
    ! I use the squared tolerance repeatedly within the following loop, so precompute it
    toler_sq = toler**2
    ! Initialize the index into the results array
    zeros_j = 1
    do j = 1,n_grid_m1
      fun_at_j = func_obj%eval(grid(j))
      fun_at_j_plus1 = func_obj%eval(grid(j+1))
      prod_fun_evals = fun_at_j * fun_at_j_plus1
      if (prod_fun_evals <= 0E0_wp) then
        result_zeros(zeros_j) = zeroin(grid(j), grid(j+1), func_obj, toler)
        zeros_j = zeros_j + 1
      else if (prod_fun_evals <= toler_sq) then
        if (abs(fun_at_j) <= abs(fun_at_j_plus1)) then
          result_zeros(zeros_j) = grid(j)
        else
          result_zeros(zeros_j) = grid(j+1)
        end if
        zeros_j = zeros_j + 1
      end if
    end do
    ! Now resize the result array appropriately
    allocate (indices_of_unique(zeros_j-1), stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
      ! Return
      return
    end if
    call unirnk(result_zeros(1:(zeros_j-1)), indices_of_unique, nunique, &
        local_err_stat, local_err_msg)
    if (local_err_stat /= 0) then
      ! This means that unique ranking did not go well.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
      ! Return
      return
    end if
    allocate (temp(nunique), stat = alloc_stat , errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
      ! Return
      return
    end if
    temp = result_zeros(indices_of_unique(1:nunique))
    deallocate (result_zeros, stat = alloc_stat , errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set result_zeros to NaN
      result_zeros = ieee_value(result_zeros, ieee_quiet_nan)
      ! Return
      return
    end if
    call move_alloc(temp, result_zeros)
  end subroutine find_all_zeros
end module find_all_zeros_mod