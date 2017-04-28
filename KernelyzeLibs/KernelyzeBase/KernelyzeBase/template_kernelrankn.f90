! template_kernelrankn.f90
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
! created on: 2015-07-07
! updated on: 2015-11-12
! updated on: 2016-01-16 (to add derivatives)
! updated on: 2016-02-11 (to add higher-order derivatives)
! updated on: 2016-07-05 (to accommodate asymmetry)
! updated on: 2016-08-22 (added kernel_integral_svar)
! updated on: 2016-11-05 (Added code to account for analytical integral
!             calculations [this uses calculate_integral_kernel_svar_impl];
!             also change to avoid code duplication in kernel_integral and
!             kernel_integal_svar)
! updated on: 2016-11-06 (revert to old submodule approach for
!             calculate_integral_kernel_svar)
! updated on: 2017-04-27 (interfaces adjusted for greater flexibility
!             for kernels without analytical partials)
!
! This is an include file that helps to achieve
! "poor man's templating" for rank-$n$
! kernels.  The natural template parameters are
! the type of the single-variable functions
! that are used in building the rank-$n$ kernel
! and the resulting type of integrals against
! the kernel.

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len, err_msg_len
use kernel_mod, only : kernel
use singlevar_func_mod, only : singlevar_func

implicit none

private

! Cannot parameterize this derived type, because the point
! of using rank as a length parameter would be to dimension
! the arrays f_array and g_array appropriately, but both
! of these arrays must be of abstract derived type -- they
! need to exhibit polymorphism, so they must be allocatable
! (or pointer -- I implement as allocatable).
type, public, extends(kernel) :: kernel_rankn
! An array of component functions to account for the $n$ functions
! of $x$.
class(singlevar_type), dimension(:), allocatable, private :: f_array
! An array of component functions to account for the $n$ functions
! of $y$.
class(singlevar_type), dimension(:), allocatable, private :: g_array
! A coefficient matrix to use in conjunction with the f_array
! component functions
real(wp), dimension(: , :), allocatable, private          :: v_matrix
! A coefficient matrix to use in conjunction with the g_array
! component functions
real(wp), dimension(: , :), allocatable, private          :: w_matrix
contains
! The evaluation procedure for a rank-$n$ kernel, this
! implements the deferred evaluation procedure of the
! parent abstract derived type kernel.
procedure, pass :: eval => eval_rankn_kernel
! First and second partial derivatives
procedure, pass                             :: first_partial_dx
procedure, pass                             :: first_partial_dy
procedure, pass                             :: second_partial_dx_dx
procedure, pass                             :: second_partial_dy_dy
procedure, pass                             :: second_partial_dx_dy
procedure, pass                             :: mth_nth_partial
! An inspector to obtain the rank of the kernel
procedure, pass :: get_rank
! An inspector to obtain the number of component functions.
! This number may be different from the rank if the V
! matrix is not square
procedure, pass :: get_num_component_funcs
! Get the number of f functions
procedure, pass :: get_size_f_array
! Get the number of g functions
procedure, pass :: get_size_g_array
! Get a single function of the f_array by its index,
! as a singlevar_func (not a singlevar_type)
procedure, pass :: get_f_func
! An inspector for the f_array component
procedure, pass :: get_f_array
! A setter for the f_array component
procedure, pass :: set_f_array
! Get a single function of the g_array by its index,
! as a singlevar_func (not a singlevar_type)
procedure, pass :: get_g_func
! An inspector for the g_array component
procedure, pass :: get_g_array
! A setter for the g_array component
procedure, pass :: set_g_array
! An inspector for the v_matrix component
procedure, pass :: get_v_matrix
! A setter for the v_matrix component
procedure, pass :: set_v_matrix
! An inspector for the w_matrix component
procedure, pass :: get_w_matrix
! A setter for the w_matrix component
procedure, pass :: set_w_matrix
! Function to evaluate the f_array functions at points
procedure, pass :: eval_f_array
! Function to evaluate the g_array functions at points
procedure, pass :: eval_g_array
! Construct the integral of this kernel with respect to
! a discrete signed measure over either $x$ or $y$.  Because
! this is a rank-$n$ kernel, the integral should be built
! differently than for a full-rank (that is, infinite-rank)
! kernel
procedure, pass :: kernel_integral
! Build an object that integrates the kernel
! over either $x$ or $y$ using a singlevar_func.
procedure, pass :: kernel_integral_svar
end type kernel_rankn
  
contains
  
! The key evaluation function
pure recursive function eval_rankn_kernel(this, x, y) result(eval_res)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  ! Function result
  real(wp)                        :: eval_res
  ! Body
  ! Use the fact that the eval procedure of singlevar_type is
  ! elemental.
  eval_res = dot_product( &
      matmul( this%f_array%eval_elt(x) , this%v_matrix ) , &
      matmul( this%g_array%eval_elt(y) , this%w_matrix ) )
end function eval_rankn_kernel
! Partial derivative functions:

elemental function first_partial_dx(this, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Body
  ! Use the fact that the eval and first_deriv procedures of 
  ! singlevar_type are elemental.
  if (present(d)) then
    deriv = dot_product( &
        matmul( this%f_array%first_deriv(x, d) , this%v_matrix ) , &
        matmul( this%g_array%eval_elt(y) , this%w_matrix ) )
  else
    deriv = dot_product( &
        matmul( this%f_array%first_deriv(x) , this%v_matrix ) , &
        matmul( this%g_array%eval_elt(y) , this%w_matrix ) )
  end if
end function first_partial_dx
  
elemental function first_partial_dy(this, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Body
  ! Use the fact that the eval and first_deriv procedures of 
  ! singlevar_type are elemental.
  if (present(d)) then
    deriv = dot_product( &
        matmul( this%f_array%eval_elt(x) , this%v_matrix ) , &
        matmul( this%g_array%first_deriv(y, d) , this%w_matrix ) )
  else
    deriv = dot_product( &
        matmul( this%f_array%eval_elt(x) , this%v_matrix ) , &
        matmul( this%g_array%first_deriv(y) , this%w_matrix ) )
  end if
end function first_partial_dy
  
elemental function second_partial_dx_dx(this, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Body
  ! Use the fact that the eval and second_deriv procedures of 
  ! singlevar_type are elemental.
  if (present(d)) then
    deriv = dot_product( &
        matmul( this%f_array%second_deriv(x, d) , this%v_matrix ) , &
        matmul( this%g_array%eval_elt(y) , this%w_matrix ) )
  else
    deriv = dot_product( &
        matmul( this%f_array%second_deriv(x) , this%v_matrix ) , &
        matmul( this%g_array%eval_elt(y) , this%w_matrix ) )
  end if
end function second_partial_dx_dx
  
elemental function second_partial_dy_dy(this, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Body
  ! Use the fact that the eval and second_deriv procedures of 
  ! singlevar_type are elemental.
  if (present(d)) then
    deriv = dot_product( &
        matmul( this%f_array%eval_elt(x) , this%v_matrix ) , &
        matmul( this%g_array%second_deriv(y, d) , this%w_matrix ) )
  else
    deriv = dot_product( &
        matmul( this%f_array%eval_elt(x) , this%v_matrix ) , &
        matmul( this%g_array%second_deriv(y) , this%w_matrix ) )
  end if
end function second_partial_dy_dy
  
elemental function second_partial_dx_dy(this, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Body
  ! Use the fact that the first_deriv procedure of 
  ! singlevar_type is elemental.
  if (present(d)) then
    deriv = dot_product( &
        matmul( this%f_array%first_deriv(x, d) , this%v_matrix ) , &
        matmul( this%g_array%first_deriv(y, d) , this%w_matrix ) )
  else
    deriv = dot_product( &
        matmul( this%f_array%first_deriv(x) , this%v_matrix ) , &
        matmul( this%g_array%first_deriv(y) , this%w_matrix ) )
  end if
end function second_partial_dx_dy

pure function mth_nth_partial(this, m, n, x, y, d) result(deriv)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  integer, intent(in)             :: m
  integer, intent(in)             :: n
  real(wp), intent(in)            :: x
  real(wp), intent(in)            :: y
  real(wp), intent(in), optional :: d
  ! Function result
  real(wp)                        :: deriv
  ! Local variables
  integer                         :: i
  real(wp)                        :: f_deriv(size(this%f_array))
  real(wp)                        :: g_deriv(size(this%g_array))
  ! Body
  ! Use the nth_deriv procedure of singlevar_type.
  if (present(d)) then
    do i = 1, size(f_deriv)
      f_deriv(i) = this%f_array(i)%nth_deriv(m, x, d)
    end do
    do i = 1, size(g_deriv)
      g_deriv(i) = this%g_array(i)%nth_deriv(n, y, d)
    end do
    deriv = dot_product( &
        matmul( f_deriv , this%v_matrix ) , &
        matmul( g_deriv , this%w_matrix ) )
  else
    do i = 1, size(f_deriv)
      f_deriv(i) = this%f_array(i)%nth_deriv(m, x)
    end do
    do i = 1, size(g_deriv)
      g_deriv(i) = this%g_array(i)%nth_deriv(n, y)
    end do
    deriv = dot_product( &
        matmul( f_deriv , this%v_matrix ) , &
        matmul( g_deriv , this%w_matrix ) )
  end if
end function mth_nth_partial

! An inspector to return the rank of the kernel
pure function get_rank(this) result(rank_res)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  ! Function result
  integer                         :: rank_res
  ! Body
  ! Should really only be necessary to include either f or g
  ! in the below, but including both seems safer.
  rank_res = min( &
      size(this%f_array), &
      size(this%g_array), &
      size(this%v_matrix, 1), &
      size(this%v_matrix, 2), &
      size(this%w_matrix, 1), &
      size(this%w_matrix, 2))
end function get_rank
! An inspector to return the number of component functions
pure function get_num_component_funcs(this) result(num_res)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  ! Function result
  integer                         :: num_res
  ! Body
  ! Should really only be necessary to include either f
  ! or g below, but taking the maximum seems safer.
  num_res = max( &
      size(this%f_array), &
      size(this%g_array) )
end function get_num_component_funcs
! An inspector to return the size of the f_array
pure function get_size_f_array(this) result(num_res)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  ! Function result
  integer                         :: num_res
  ! Body
  num_res = size(this%f_array)
end function get_size_f_array
! An inspector to return the size of the g_array
pure function get_size_g_array(this) result(num_res)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  ! Function result
  integer                         :: num_res
  ! Body
  num_res = size(this%g_array)
end function get_size_g_array
! Get a single function of the f_array by its index,
! as a singlevar_func (not a singlevar_type)
subroutine get_f_func(this, res_func, func_index, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  integer, intent(in)             :: func_index
  class(singlevar_func), &
        allocatable, intent(out)  :: res_func
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "get_f_func procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! The below select type statement seems totally unnecessary . . .
  ! however, if I do not use it (if I attempt to use
  ! this%f_array(func_index) directly as the source of the
  ! allocation for res_func), I run into what appears to be a
  ! compiler bug: instead of using the stated index of the
  ! f_array as the source, the compiler uses the *first*
  ! index of f_array as the source.
  select type(temp_ind_single => this%f_array(func_index))
  class is (singlevar_type)
    allocate( &
      res_func, &
      source = temp_ind_single, &
      stat = alloc_stat, &
      errmsg = alloc_err_msg)
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
  end select
end subroutine get_f_func
! An inspector for the f_array component
subroutine get_f_array(this, curr_f_array, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(in)         :: this
  class(singlevar_type), dimension(:), &
      allocatable, intent(out)            :: curr_f_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "get_f_array procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! Allocate the curr_f_array, copying into it the correct
  ! dynamic type and type component values from
  ! the f_array component of this.
  allocate(curr_f_array, source = this%f_array, &
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
end subroutine get_f_array
! A setter for the f_array component
subroutine set_f_array(this, new_f_array, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(inout)              :: this
  class(singlevar_type), dimension(:), intent(in) :: new_f_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "set_f_array procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! If the f_array is currently allocated, deallocate it
  if (allocated(this%f_array)) then
    deallocate(this%f_array, &
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
  ! Allocate the f_array, using the dynamic type and type components 
  ! of the new_f_array
  allocate(this%f_array, source = new_f_array, &
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
end subroutine set_f_array
! Get a single function of the g_array by its index,
! as a singlevar_func (not a singlevar_type)
subroutine get_g_func(this, res_func, func_index, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(in) :: this
  integer, intent(in)             :: func_index
  class(singlevar_func), &
        allocatable, intent(out)  :: res_func
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "get_g_func procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! The below select type statement seems totally unnecessary . . .
  ! however, if I do not use it (if I attempt to use
  ! this%g_array(func_index) directly as the source of the
  ! allocation for res_func), I run into what appears to be a
  ! compiler bug: instead of using the stated index of the
  ! g_array as the source, the compiler uses the *first*
  ! index of g_array as the source.
  select type(temp_ind_single => this%g_array(func_index))
  class is (singlevar_type)
    allocate( &
      res_func, &
      source = temp_ind_single, &
      stat = alloc_stat, &
      errmsg = alloc_err_msg)
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
  end select
end subroutine get_g_func
! An inspector for the g_array component
subroutine get_g_array(this, curr_g_array, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(in)         :: this
  class(singlevar_type), dimension(:), &
      allocatable, intent(out)            :: curr_g_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "get_g_array procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! Allocate the curr_g_array, copying into it the correct
  ! dynamic type and type component values from
  ! the g_array component of this.
  allocate(curr_g_array, source = this%g_array, &
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
end subroutine get_g_array
! A setter for the g_array component
subroutine set_g_array(this, new_g_array, err_stat, err_msg)
  ! Arguments
  class(kernel_rankn), intent(inout)              :: this
  class(singlevar_type), dimension(:), intent(in) :: new_g_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "set_g_array procedure of " &
      // "the kernel_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! If the g_array is currently allocated, deallocate it
  if (allocated(this%g_array)) then
    deallocate(this%g_array, &
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
  ! Allocate the g_array, using the dynamic type and type components 
  ! of the new_g_array
  allocate(this%g_array, source = new_g_array, &
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
end subroutine set_g_array
! An inspector for the v_matrix component
pure function get_v_matrix(this) result(curr_v_matrix)
  ! Arguments
  class(kernel_rankn), intent(in)         :: this
  ! Function result
  real(wp), dimension(: , :), allocatable :: curr_v_matrix
  ! Body
  curr_v_matrix = this%v_matrix
end function get_v_matrix
! A setter for the v_matrix component
pure subroutine set_v_matrix(this, new_v_matrix)
  ! Arguments
  class(kernel_rankn), intent(inout)              :: this
  real(wp), dimension(: , :), intent(in)          :: new_v_matrix
  ! Body
  this%v_matrix = new_v_matrix
end subroutine set_v_matrix
! An inspector for the w_matrix component
pure function get_w_matrix(this) result(curr_w_matrix)
  ! Arguments
  class(kernel_rankn), intent(in)         :: this
  ! Function result
  real(wp), dimension(: , :), allocatable :: curr_w_matrix
  ! Body
  curr_w_matrix = this%w_matrix
end function get_w_matrix
! A setter for the w_matrix component
pure subroutine set_w_matrix(this, new_w_matrix)
  ! Arguments
  class(kernel_rankn), intent(inout)              :: this
  real(wp), dimension(: , :), intent(in)          :: new_w_matrix
  ! Body
  this%w_matrix = new_w_matrix
end subroutine set_w_matrix
! Function to evaluate the f_array functions at points
pure function eval_f_array(this , eval_at) result(eval_res)
  ! Arguments
  class(kernel_rankn), intent(in)     :: this
  real(wp), dimension(:), intent(in)  :: eval_at
  ! Function result
  real(wp), dimension(size(eval_at), size(this%f_array))  :: eval_res
  ! Local variables
  integer :: j
  ! Body
  ! Use the fact that the eval procedure of singlevar_type is elemental
  do j = 1 , size(this%f_array)
    eval_res(: , j) = this%f_array(j)%eval_elt(eval_at)
  end do
end function eval_f_array
! Function to evaluate the g_array functions at points
pure function eval_g_array(this , eval_at) result(eval_res)
  ! Arguments
  class(kernel_rankn), intent(in)     :: this
  real(wp), dimension(:), intent(in)  :: eval_at
  ! Function result
  real(wp), dimension(size(eval_at), size(this%g_array))  :: eval_res
  ! Local variables
  integer :: j
  ! Body
  ! Use the fact that the eval procedure of singlevar_type is elemental
  do j = 1 , size(this%g_array)
    eval_res(: , j) = this%g_array(j)%eval_elt(eval_at)
  end do
end function eval_g_array
    
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
  class(kernel_rankn), intent(in)                   :: this
  logical, intent(in)                               :: sum_is_over_x
  real(wp), dimension(:), intent(in)                :: weights
  real(wp), dimension(size(weights)), intent(in)    :: eval_pts
  class(integral_type), allocatable, &
                                        intent(out) :: integral
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle information on any local errors
  integer                         :: local_err_stat
  character(len=err_msg_len)      :: local_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "kernel_integral procedure of " &
      // "the kernel_rankn derived type: "
  ! Local variables
  real(wp), dimension(this%get_num_component_funcs()) :: temp_coeffs
  class(singlevar_type), dimension(:), allocatable    :: temp_compfunc_array
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! Check weights and eval_pts inputs
  if (size(weights) /= size(eval_pts)) then
    ! If err_stat and / or err_msg were provided, fill them in
    if (present(err_stat)) then
      err_stat = -1
    end if
    if (present(err_msg)) then
      err_msg = proc_name // 'eval_pts and weights must have equal sizes'
    end if
    return
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
  call integral%set_integral_vbl(sum_is_over_x)
  call integral%set_description('Integral built from rank-n kernel.')
  call integral%set_eval_pts(eval_pts)
  call integral%set_weights(weights)
  if (sum_is_over_x) then
    call integral%set_v_matrix(this%get_w_matrix())
  else
    call integral%set_v_matrix(this%get_v_matrix())
  end if
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
  ! Now to compute the coefficient vector
  ! First put the proper array of singlevar_types
  ! into comp_array
  if (sum_is_over_x) then
    call this%get_g_array( &
        temp_compfunc_array, &
        local_err_stat, &
        local_err_msg)
  else
    call this%get_f_array( &
        temp_compfunc_array, &
        local_err_stat, &
        local_err_msg)
  end if
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
  call integral%set_comp_array( &
      temp_compfunc_array, &
      local_err_stat, &
      local_err_msg)
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
  ! Now compute the temp_coeffs:
  if (sum_is_over_x) then
    temp_coeffs = matmul( weights , this%eval_f_array(eval_pts) )
  else
    temp_coeffs = matmul( weights , this%eval_g_array(eval_pts) )
  end if
  ! Account for the V or W matrix of the rank-$n$ kernel.
  ! This computation to is carefully ordered to improve
  ! accuracy in high-precision cases (large-rank approximations),
  ! since in those cases floating-point arithmetic means that
  ! order matters: operations that are associative in exact
  ! arithmetic are not associative in floating-point arithmetic.
  if (sum_is_over_x) then
    call integral%set_coeff_vec( &
        matmul( temp_coeffs , this%get_v_matrix() ) )
  else
    call integral%set_coeff_vec( &
        matmul( transpose( this%get_w_matrix() ) , &
            temp_coeffs ))
  end if
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
  class(kernel_rankn), intent(in)                   :: this
  logical, intent(in)                               :: int_is_over_x
  class(singlevar_func), allocatable, intent(in)    :: svar
  class(integral_sv_type), allocatable, &
                                        intent(out) :: integral
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle information on any local errors
  integer                         :: local_err_stat
  character(len=err_msg_len)      :: local_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "kernel_integral_svar " &
      // "procedure of the kernel_rankn derived type: "
  ! Local variables
  integer                                             :: i
  real(wp), dimension(this%get_num_component_funcs()) :: temp_coeffs
  class(singlevar_type), dimension(:), allocatable    :: temp_compfunc_array
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
  call integral%set_description('Integral vs svar from rank-n kernel.')
  if (int_is_over_x) then
    call integral%set_v_matrix(this%get_w_matrix())
  else
    call integral%set_v_matrix(this%get_v_matrix())
  end if
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
  ! Now to compute the coefficient vector
  ! First put the proper array of singlevar_types
  ! into comp_array
  if (int_is_over_x) then
    call this%get_g_array( &
        temp_compfunc_array, &
        local_err_stat, &
        local_err_msg)
  else
    call this%get_f_array( &
        temp_compfunc_array, &
        local_err_stat, &
        local_err_msg)
  end if
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
  call integral%set_comp_array( &
      temp_compfunc_array, &
      local_err_stat, &
      local_err_msg)
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
  ! Now compute the temp_coeffs:
  if (int_is_over_x) then
    do i = 1, size(this%f_array)
      temp_coeffs(i) = this%f_array(i)%integrate_svar( svar )
    end do
  else
    do i = 1, size(this%g_array)
      temp_coeffs(i) = this%g_array(i)%integrate_svar( svar )
    end do
  end if
  ! Account for the V or W matrix of the rank-$n$ kernel.
  ! This computation to is carefully ordered to improve
  ! accuracy in high-precision cases (large-rank approximations),
  ! since in those cases floating-point arithmetic means that
  ! order matters: operations that are associative in exact
  ! arithmetic are not associative in floating-point arithmetic.
  if (int_is_over_x) then
    call integral%set_coeff_vec( &
        matmul( temp_coeffs , this%get_v_matrix() ) )
  else
    call integral%set_coeff_vec( &
        matmul( transpose( this%get_w_matrix() ) , &
            temp_coeffs ))
  end if
end subroutine kernel_integral_svar    