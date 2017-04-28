! integral_kernel_disc.f90
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
! created on: 2015-04-02
! updated on: 2015-07-10
! assignment added on : 2015-11-12
! updated on: 2016-01-09 (derivatives)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! This module contains a derived type extending the abstract 
! derived type integral_kernel_measure. This extension restricts 
! to the case of a discrete signed measure.
  
module integral_kernel_disc_mod

  use set_precision, only : wp
  use singlevar_func_mod, only : singlevar_func
  use kernel_mod, only : kernel
  use integral_kernel_measure_mod, only : integral_kernel_measure

  implicit none

  private

  type, public, extends(integral_kernel_measure) :: integral_kernel_disc
      ! The support points of the discrete signed measure
      real(wp), dimension(:), allocatable, private  :: eval_pts
      ! The weights of the discrete signed measure
      real(wp), dimension(:), allocatable, private  :: weights
  contains
      ! Implement the deferred procedure of the
      ! parent abstract derived type
      procedure, pass :: eval => eval_integral_disc
      ! Override the singlevar_func derivative procedures
      procedure, pass :: first_deriv
      procedure, pass :: second_deriv
      procedure, pass :: nth_deriv
      ! Implement the deferred assignment procedure
      procedure       :: singlevar_assign => integral_assign
      ! Basic inspector for the eval_pts
      procedure, pass :: get_eval_pts
      ! Basic setter for the eval_pts
      procedure, pass :: set_eval_pts
      ! Basic inspector for the weights
      procedure, pass :: get_weights
      ! Basic setter for the weights
      procedure, pass :: set_weights
      ! Effective eval_pts: here just the eval_pts, 
      ! but derived types may override
      procedure, pass :: get_effective_eval_pts
      ! Effective weights: here just the eval_pts, 
      ! but derived types may override
      procedure, pass :: get_effective_weights
  end type integral_kernel_disc

contains
  
  pure recursive function eval_integral_disc(this, arg_pt) result(eval_res)
    ! Arguments
    class(integral_kernel_disc), intent(in)     :: this
    real(wp), intent(in)                        :: arg_pt
    ! Function result
    real(wp)                                    :: eval_res
    ! Local variables
    integer   :: i
    real(wp)  :: kernel_eval_res(size(this%eval_pts))
    ! Body
    ! Use elemental kernel evaluation; less code and
    ! ease of parallelization.
    if (this%get_integral_vbl()) then
      do i = 1, size(this%eval_pts)
        kernel_eval_res(i) = this%eval_kernel(this%eval_pts(i), arg_pt)
      end do
      eval_res = dot_product( this%weights , kernel_eval_res )
    else
      do i = 1, size(this%eval_pts)
        kernel_eval_res(i) = this%eval_kernel(arg_pt, this%eval_pts(i))
      end do
      eval_res = dot_product( this%weights , kernel_eval_res )
    end if
  end function eval_integral_disc
  
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    real(wp), intent(in)                    :: x
    real(wp), intent(in), optional          :: d
    ! Function result
    real(wp)                                :: deriv
    ! Local variables
    integer                                 :: i
    ! Body
    deriv = 0E0_wp
    if (present(d)) then
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_dkernel_dy(this%eval_pts(i) , x , d)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_dkernel_dx(x , this%eval_pts(i) , d)
        end do
      end if
    else
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_dkernel_dy(this%eval_pts(i) , x)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_dkernel_dx(x , this%eval_pts(i))
        end do
      end if
    end if
  end function first_deriv
  
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    real(wp), intent(in)                    :: x
    real(wp), intent(in), optional          :: d
    ! Function result
    real(wp)                                :: deriv
    ! Local variables
    integer                                 :: i
    ! Body
    deriv = 0E0_wp
    if (present(d)) then
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_d2kernel_dy_dy(this%eval_pts(i) , x , d)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_d2kernel_dx_dx(x , this%eval_pts(i) , d)
        end do
      end if
    else
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_d2kernel_dy_dy(this%eval_pts(i) , x)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_d2kernel_dx_dx(x , this%eval_pts(i))
        end do
      end if
    end if
  end function second_deriv
  
  pure function nth_deriv(this, n, x, d) result(deriv)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    integer, intent(in)                     :: n
    real(wp), intent(in)                    :: x
    real(wp), intent(in), optional          :: d
    ! Function result
    real(wp)                                :: deriv
    ! Local variables
    integer                                 :: i
    ! Body
    deriv = 0E0_wp
    if (present(d)) then
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_mth_nth_partial_kernel( &
              0, n, this%eval_pts(i) , x , d)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_mth_nth_partial_kernel( &
              n, 0, x , this%eval_pts(i) , d)
        end do
      end if
    else
      if (this%get_integral_vbl()) then
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_mth_nth_partial_kernel(0, n, this%eval_pts(i) , x)
        end do
      else
        do i = 1, size(this%weights)
          deriv = deriv + this%weights(i) &
              * this%eval_mth_nth_partial_kernel(n, 0, x , this%eval_pts(i))
        end do
      end if
    end if
  end function nth_deriv
  
  ! Implement the deferred assignment operation
  subroutine integral_assign(left, right)
    ! Arguments
    class(integral_kernel_disc), intent(out)  :: left
    class(singlevar_func), intent(in)         :: right
    ! Local variables
    class(kernel), allocatable                :: kernel_used
    ! Body
    select type (right)
    type is (integral_kernel_disc)
      left%eval_pts = right%eval_pts
      left%weights  = right%weights
      call left%set_integral_vbl(right%get_integral_vbl())
      call left%set_description(right%get_description())
      call right%get_kernel(kernel_used)
      call left%set_kernel(kernel_used)
    end select
  end subroutine
  
  ! A simple inspector for the evaluation points
  pure function get_eval_pts(this) result(curr_eval_pts)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    ! Result
    real(wp), dimension(:), allocatable     :: curr_eval_pts
    ! Body
    curr_eval_pts = this%eval_pts
  end function get_eval_pts
  
  ! A simple setter for the eval_pts
  pure subroutine set_eval_pts(this, new_eval_pts)
    ! Arguments
    class(integral_kernel_disc), intent(inout)  :: this
    real(wp), dimension(:), intent(in)          :: new_eval_pts
    ! Body
    this%eval_pts = new_eval_pts
  end subroutine set_eval_pts
  
  ! A simple inspector for the weights
  pure function get_weights(this) result(curr_weights)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    ! Result
    real(wp), dimension(:), allocatable     :: curr_weights
    ! Body
    curr_weights = this%weights
  end function get_weights
  
  ! A simple setter for the weights
  pure subroutine set_weights(this, new_weights)
    ! Arguments
    class(integral_kernel_disc), intent(inout)  :: this
    real(wp), dimension(:), intent(in)          :: new_weights
    ! Body
    this%weights = new_weights
  end subroutine set_weights
  
  ! Effective eval_pts: here just the eval_pts, 
  ! but derived types may override
  pure function get_effective_eval_pts(this) result(curr_eff_eval_pts)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    ! Result
    real(wp), dimension(:), allocatable     :: curr_eff_eval_pts
    ! Body
    curr_eff_eval_pts = this%get_eval_pts()
  end function get_effective_eval_pts
  
  ! Effective weights: here just the eval_pts, 
  ! but derived types may override
  pure function get_effective_weights(this) result(curr_eff_weights)
    ! Arguments
    class(integral_kernel_disc), intent(in) :: this
    ! Result
    real(wp), dimension(:), allocatable     :: curr_eff_weights
    ! Body
    curr_eff_weights = this%get_weights()
  end function get_effective_weights

end module integral_kernel_disc_mod