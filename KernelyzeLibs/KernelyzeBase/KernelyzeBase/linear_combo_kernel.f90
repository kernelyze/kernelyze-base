! linear_combo_kernel.f90
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
! created on: 2016-04-24
! updated on: 2016-04-24
! updated on: 2016-06-30 (Added switch to allow
!             combos that are functions of y as
!             well as combos that are functions
!             of x [the latter remains the default])
! updated on: 2016-07-01 (Continue the 2016-06-30
!             work, no longer permit the choice
!             of variable to default)
! updated on: 2017-04-27 (added flexibility to finite-difference step size)
!
! A module containing the linear_combo_kernel
! derived type, which extends singlevar_func
! in a way that is useful in several kernel-approximation
! computations.
  
module linear_combo_kernel_mod
  
  use set_precision, only : wp
  use kernel_mod, only : kernel
  use singlevar_func_mod, only : singlevar_func
  
  implicit none
  
  ! Extend singlevar_func to obtain a derived type
  ! that can capture the quantity
  ! coeffs' * (a_matrix * kernel(x, rho_vec))
  type, public, extends(singlevar_func) :: linear_combo_kernel
    ! If true, the linear combo is a function of x; else, it
    ! is a function of y
    logical, private                    :: is_func_of_x
    ! The rho vector used to compute the combo
    real(wp), private, allocatable      :: combo_rho_vec(:)
    ! The A matrix used to compute the combo
    real(wp), private, allocatable      :: combo_a_matrix(: , :)
    ! The coefficients used to compute the combo
    real(wp), private, allocatable      :: combo_coeffs(:)
    ! The kernel object used to compute the combo
    class(kernel), private, allocatable :: combo_kernel
  contains
    ! Initialize the combo data members (rho vector and A matrix)
    ! (note that I do not implement the usual inspectors and setters
    ! because this is an implementation object -- it should be used
    ! only in remez_step, and it is marked private here)
    procedure, pass :: init_combo
    ! Implement the deferred non-elemental (but still pure)
    ! evaluation procedure of singlevar_func
    procedure, pass :: eval
    ! Derivatives
    procedure, pass :: first_deriv
    procedure, pass :: second_deriv
    ! Could provide an nth_deriv, but the base kernel procedure
    ! is not elemental (just pure), so the implementation is
    ! a bit more cumbersome and it is not currently useful.
    ! procedure, pass :: nth_deriv
    ! Implement the deferred assignment procedure
    procedure       :: singlevar_assign => combo_assign 
  end type linear_combo_kernel
  
  contains
  
  ! Implement the deferred procedure for non-elemental
  ! (but still pure) evaluation of singlevar_func for the
  ! combo derived type.
  pure recursive function eval(this, arg_pt) result(eval_res)
    ! Arguments
    class(linear_combo_kernel), intent(in)  :: this
    real(wp), intent(in)              :: arg_pt
    ! Function result
    real(wp)                          :: eval_res
    ! Body
    if (this%is_func_of_x) then
      eval_res = dot_product( &
          this%combo_coeffs, &
         matmul( this%combo_a_matrix , &
                 this%combo_kernel%eval_elt( arg_pt , this%combo_rho_vec ) ) )
    else
      eval_res = dot_product( &
          this%combo_coeffs, &
         matmul( this%combo_a_matrix , &
                 this%combo_kernel%eval_elt( this%combo_rho_vec , arg_pt ) ) )
    end if
  end function eval
  
  ! Note that this "x" argument could actually be applied
  ! to either the "x" or "y" argument positions of the kernel.
  ! The need to name this argument "x" is enforced by the Intel
  ! compiler ("The dummy arguments of an overriding and overridden
  ! binding that correspond by position must have the same names.")
  elemental function first_deriv(this, x, d) result(deriv)
    ! Arguments
    class(linear_combo_kernel), intent(in)  :: this
    real(wp), intent(in)                    :: x
    real(wp), intent(in), optional          :: d
    ! Function result
    real(wp)                                :: deriv
    ! Body
    if (present(d)) then
      if (this%is_func_of_x) then
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%first_partial_dy( &
                      x , this%combo_rho_vec , d) ) )
      else
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%first_partial_dy( &
                      this%combo_rho_vec , x , d) ) )
      end if
    else
      if (this%is_func_of_x) then
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%first_partial_dy( &
                      x , this%combo_rho_vec ) ) )
      else
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%first_partial_dy( &
                      this%combo_rho_vec , x ) ) )
      end if
    end if
  end function first_deriv
  
  ! Note that this "x" argument could actually be applied
  ! to either the "x" or "y" argument positions of the kernel.
  ! The need to name this argument "x" is enforced by the Intel
  ! compiler ("The dummy arguments of an overriding and overridden
  ! binding that correspond by position must have the same names.")
  elemental function second_deriv(this, x, d) result(deriv)
    ! Arguments
    class(linear_combo_kernel), intent(in)  :: this
    real(wp), intent(in)                    :: x
    real(wp), intent(in), optional          :: d
    ! Function result
    real(wp)                                :: deriv
    ! Body
    if (present(d)) then
      if (this%is_func_of_x) then
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%second_partial_dy_dy( &
                      x , this%combo_rho_vec , d ) ) )
      else
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%second_partial_dy_dy( &
                      this%combo_rho_vec , x , d ) ) )
      end if
    else
      if (this%is_func_of_x) then
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%second_partial_dy_dy( &
                      x , this%combo_rho_vec ) ) )
      else
        deriv = dot_product( &
            this%combo_coeffs, &
           matmul( this%combo_a_matrix , &
                   this%combo_kernel%second_partial_dy_dy( &
                      this%combo_rho_vec , x ) ) )
      end if
    end if
  end function second_deriv
  
  ! Implement the deferred assignment procedure for combo objects
  subroutine combo_assign(left, right)
    ! Arguments
    class(linear_combo_kernel), intent(out) :: left
    class(singlevar_func), intent(in) :: right
    ! Body
    select type (right)
    type is (linear_combo_kernel)
      left%is_func_of_x   = right%is_func_of_x
      left%combo_rho_vec  = right%combo_rho_vec
      left%combo_a_matrix = right%combo_a_matrix
      left%combo_coeffs   = right%combo_coeffs
      if (allocated(left%combo_kernel)) then
        deallocate(left%combo_kernel)
      end if
      allocate(left%combo_kernel, source = right%combo_kernel)
    end select
  end subroutine combo_assign
  
  ! Initialize a combo object
  subroutine init_combo( &
      this, &
      kernel_to_use, &
      rho_vec_to_use, &
      a_matrix_to_use, &
      coeffs_to_use, &
      func_of_x)
    ! Arguments
    class(linear_combo_kernel), intent(out) :: this
    class(kernel), intent(in)               :: kernel_to_use
    real(wp), intent(in)                    :: rho_vec_to_use(:)
    real(wp), intent(in)                    :: a_matrix_to_use(: , :)
    real(wp), intent(in)                    :: coeffs_to_use(:)
    logical, intent(in)                     :: func_of_x
    ! Body
    this%is_func_of_x = func_of_x
    this%combo_rho_vec  = rho_vec_to_use
    this%combo_a_matrix = a_matrix_to_use
    this%combo_coeffs   = coeffs_to_use
    if (allocated(this%combo_kernel)) then
      deallocate(this%combo_kernel)
    end if
    allocate(this%combo_kernel, source = kernel_to_use)
  end subroutine init_combo

end module linear_combo_kernel_mod