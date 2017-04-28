! template_integralrankn.f90
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
! updated on: 2015-07-18
! assignment added on : 2015-11-12 and then removed
! updated on: 2016-02-07 (move get_effective_eval_pts)
! updated on: 2016-08-22 (adjust for greater flexibility in
!             the parent_type being extended)
! updated on: 2016-08-29 (expose interface for eval given
!             values of the compfunc_kernel array)
! updated on: 2017-04-23 (remove use statements for
!             modules or parts of modules that are
!             not actually used)
!
! This is an include file that helps to achieve
! "poor man's templating" for rank-$n$
! integrals.  The natural template parameter is 
! the type of the single-variable functions
! that are used in building the rank-$n$
! integral.

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len, err_msg_len

implicit none

private

! The integral is represented by 
! dot_product(coeff_vec , matmul( comp_array( vbl ) , v_matrix), 
! where vbl is either $x$ (if integration was over $y$) or $y$ 
! (if integration was over $x$), and comp_array is either comp_array 
! (if integration was over $y$) or g_array (if integration was 
! over $x$).
type, public, extends(parent_type) :: integral_rankn
  ! The coefficient vector
  real(wp), dimension(:), allocatable, private    :: coeff_vec
  ! The component functions; treating these functions as a
  ! vector, first multiply the vector by the matrix v_matrix
  ! and then dot-product the result with the coeff_vec to
  ! evaluate the integral_rankn.  Why not collapse everything
  ! to a single dot-product?  Order matters in floating-point
  ! arithmetic (operations that are associative in exact
  ! arithmetic are not associative in floating-point
  ! arithmetic), and keeping this order of operations can
  ! matter in high-precision cases, as when the kernel is a
  ! rank-$n$ approximation where $n$ is rather large.
  class(singlevar_type), dimension(:), &
                            allocatable, private  :: comp_array
  ! A coefficient matrix
  real(wp), dimension(:, :), allocatable, private :: v_matrix
contains
  ! Implement the eval procedure
  procedure, pass :: eval_nonelt => eval_integral_rankn
  ! Evaluate given a vector representing the compfunc_kernel values
  procedure, pass :: eval_given_vals
  ! An inspector for the coefficient vector
  procedure, pass :: get_coeff_vec
  ! A setter for the coefficient vector
  procedure, pass :: set_coeff_vec
  ! An inspector for the component functions
  procedure, pass :: get_comp_array
  ! A setter for the component functions
  procedure, pass :: set_comp_array
  ! An inspector for the v_matrix component
  procedure, pass :: get_v_matrix
  ! A setter for the v_matrix component
  procedure, pass :: set_v_matrix
  ! Effective eval_pts: because of reduced rank, there may
  ! be fewer effective eval points than explicit eval pts
  procedure, pass :: get_effective_eval_pts
  ! Effective weights: because of reduced rank, there may
  ! be fewer effective weights than explicit weights
  procedure, pass :: get_effective_weights
end type integral_rankn

contains
  
pure function eval_integral_rankn(this, arg_pt) result(eval_res)
  ! Arguments
  class(integral_rankn), intent(in) :: this
  real(wp), intent(in)              :: arg_pt
  ! Function result
  real(wp)                          :: eval_res
  ! Local variable to label any problems:
  character(*), parameter :: proc_name = &
      'Evaluation of the integral of a rank-n kernel using ' // &
      'a discrete signed measure: '
  ! Body
  eval_res = this%eval_given_vals( this%comp_array%eval_elt(arg_pt) )
end function eval_integral_rankn
! Given values for the compfunc_kernels, evaluate the integral
pure function eval_given_vals(this, given_vals) result(eval_res)
  ! Arguments
  class(integral_rankn), intent(in) :: this
  real(wp), intent(in)              :: given_vals(size(this%comp_array))
  ! Function result
  real(wp)                          :: eval_res
  ! Local variable to label any problems:
  character(*), parameter :: proc_name = &
      'Evaluation of the integral of a rank-n kernel ' // &
      'given values of the compfunc_kernel array: '
  ! Body
  if (this%get_integral_vbl()) then
    eval_res = dot_product( &
        matmul( transpose(this%v_matrix) , &
            given_vals ) , & ! this%comp_array%eval(arg_pt) ) , &
        this%coeff_vec )
  else
    eval_res = dot_product( &
        matmul( given_vals, & ! this%comp_array%eval(arg_pt) , &
            this%v_matrix ) , &
        this%coeff_vec )
  end if
end function eval_given_vals
! An inspector for the coefficient vector
pure function get_coeff_vec(this) result(curr_coeffs)
  ! Arguments
  class(integral_rankn), intent(in)   :: this
  ! Function result
  real(wp), dimension(:), allocatable :: curr_coeffs
  ! Body
  curr_coeffs = this%coeff_vec
end function get_coeff_vec
! A setter for the coefficient vector
pure subroutine set_coeff_vec(this, new_coeffs)
  ! Arguments
  class(integral_rankn), intent(inout)  :: this
  real(wp), dimension(:), intent(in)    :: new_coeffs
  ! Body
  this%coeff_vec = new_coeffs
end subroutine set_coeff_vec
! An inspector for the comp_array component
subroutine get_comp_array(this, curr_comp_array, err_stat, err_msg)
  ! Arguments
  class(integral_rankn), intent(in)       :: this
  class(singlevar_type), dimension(:), &
      allocatable, intent(out)            :: curr_comp_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "get_comp_array procedure of " &
      // "the integral_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! Allocate the curr_comp_array, copying into it the correct
  ! dynamic type and type component values from
  ! the comp_array component of this.
  allocate(curr_comp_array, source = this%comp_array, &
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
end subroutine get_comp_array
! A setter for the comp_array component
subroutine set_comp_array(this, new_comp_array, err_stat, err_msg)
  ! Arguments
  class(integral_rankn), intent(inout)              :: this
  class(singlevar_type), dimension(:), intent(in)   :: new_comp_array
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  ! Local variables to handle allocation problems:
  integer                         :: alloc_stat
  character(len=alloc_errmsg_len) :: alloc_err_msg
  ! Local parameter to label any error messages:
  character(*), parameter :: proc_name = "set_comp_array procedure of " &
      // "the integral_rankn derived type: "
  ! Body
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  ! If the comp_array is currently allocated, deallocate it
  if (allocated(this%comp_array)) then
    deallocate(this%comp_array, &
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
  ! Allocate the comp_array, using the dynamic type and type components 
  ! of the new_comp_array
  allocate(this%comp_array, source = new_comp_array, &
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
end subroutine set_comp_array
! An inspector for the v_matrix component
pure function get_v_matrix(this) result(curr_v_matrix)
  ! Arguments
  class(integral_rankn), intent(in)       :: this
  ! Function result
  real(wp), dimension(: , :), allocatable :: curr_v_matrix
  ! Body
  curr_v_matrix = this%v_matrix
end function get_v_matrix
! A setter for the v_matrix component
pure subroutine set_v_matrix(this, new_v_matrix)
  ! Arguments
  class(integral_rankn), intent(inout)    :: this
  real(wp), dimension(: , :), intent(in)  :: new_v_matrix
  ! Body
  this%v_matrix = new_v_matrix
end subroutine set_v_matrix
! Effective weights
pure function get_effective_weights(this) result(curr_eff_weights)
  ! Arguments
  class(integral_rankn), intent(in)   :: this
  ! Result
  real(wp), dimension(:), allocatable :: curr_eff_weights
  ! Body
  curr_eff_weights = matmul( this%get_v_matrix(), this%get_coeff_vec() )
end function get_effective_weights