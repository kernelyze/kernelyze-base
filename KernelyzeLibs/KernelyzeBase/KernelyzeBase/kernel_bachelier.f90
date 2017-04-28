! kernel_bachelier.f90
!
! Copyright (c) 2015, 2016 by Kernelyze LLC
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
! created on: 2015-04-03
! updated on: 2015-07-03
! updated on: 2016-08-26 (use normopt_formula)
! updated on: 2016-09-27 (added is_call flag)
!
! This module contains a derived type that extends the 
! kernel_summable abstract derived type and makes it concrete:
! the extension here is to a Bachelier kernel.
  
module kernel_bachelier_mod

  use set_precision, only : wp
  use kernel_summable_mod, only : kernel_summable
  use option_formulae_mod, only : normopt_formula

  implicit none

  private

  type, public, extends(kernel_summable) :: kernel_bachelier
    ! The key constant "sigma" in the definition of
    ! the Bachelier kernel
    real(wp), private :: sigma = 1E-2_wp ! Default value
    ! Is this a call?
    logical, private  :: is_call = .true. ! Default value
  contains
    ! The evaluation procedure for a bachelier
    ! product kernel
    procedure, pass :: eval => eval_bachelier_kernel
    ! Get the parameter "alpha"
    procedure, pass :: get_sigma => get_sigma_param
    ! Set the parameter "alpha"
    procedure, pass :: set_sigma => set_sigma_param
    ! Get the is_call flag
    procedure, pass :: get_is_call
    ! Set the is_call flag
    procedure, pass :: set_is_call
  end type kernel_bachelier
  
  contains
  
  pure recursive function eval_bachelier_kernel(this, x, y) result(eval_res)
    ! Arguments
    class(kernel_bachelier), intent(in) :: this
    real(wp), intent(in)                :: x
    real(wp), intent(in)                :: y
    ! Function result
    real(wp)                            :: eval_res
    ! Body
    eval_res = normopt_formula( &
        is_call = .true., &
        strike = x, &
        forward = y, &
        vol = this%sigma, &
        disc_fac = 1E0_wp )
  end function eval_bachelier_kernel
  
  pure function get_sigma_param(this) result(sigma_param)
    ! Argument
    class(kernel_bachelier), intent(in) :: this
    ! Function result
    real(wp)                            :: sigma_param
    ! Body
    sigma_param = this%sigma
  end function get_sigma_param
  
  pure subroutine set_sigma_param(this, sigma_param)
    ! Arguments
    class(kernel_bachelier), intent(inout)  :: this
    real(wp), intent(in)                    :: sigma_param
    ! Body
    this%sigma = sigma_param
  end subroutine set_sigma_param
  
  pure function get_is_call(this) result(call_flag)
    ! Argument
    class(kernel_bachelier), intent(in) :: this
    ! Function result
    logical                             :: call_flag
    ! Body
    call_flag = this%is_call
  end function get_is_call
  
  pure subroutine set_is_call(this, call_flag)
    ! Arguments
    class(kernel_bachelier), intent(inout)  :: this
    logical, intent(in)                     :: call_flag
    ! Body
    this%is_call = call_flag
  end subroutine set_is_call

end module kernel_bachelier_mod