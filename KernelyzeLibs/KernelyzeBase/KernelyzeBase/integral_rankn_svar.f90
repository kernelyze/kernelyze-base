! integral_rankn_svar.f90
!
! Copyright (c) 2016 by Kernelyze LLC
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
! created on: 2016-10-07
!
! This module contains a derived type extending the abstract 
! derived type integral_kernel_singlevar. This extension restricts 
! to the case of a rank-$n$ integral.
  
module integral_rankn_svar_mod

  use set_precision, only : wp
  use integral_kernel_singlevar_mod, only : integral_kernel_singlevar
  
  implicit none

  private

  type, public, abstract, &
      extends(integral_kernel_singlevar) :: integral_rankn_svar
  contains
    ! An inspector for the coefficient vector
    procedure(getcoeff), pass, deferred :: get_coeff_vec
    ! A setter for the coefficient vector
    procedure(setcoeff), pass, deferred :: set_coeff_vec
    ! An inspector for the v_matrix component
    procedure(getvmat), pass, deferred  :: get_v_matrix
    ! A setter for the v_matrix component
    procedure(setvmat), pass, deferred  :: set_v_matrix
  end type integral_rankn_svar
  
  abstract interface
    ! An inspector for the coefficient vector
    pure function getcoeff(this) result(curr_coeffs)
      import wp
      import integral_rankn_svar
      ! Arguments
      class(integral_rankn_svar), intent(in)  :: this
      ! Function result
      real(wp), dimension(:), allocatable     :: curr_coeffs
    end function getcoeff
    ! A setter for the coefficient vector
    pure subroutine setcoeff(this, new_coeffs)
      import wp
      import integral_rankn_svar
      ! Arguments
      class(integral_rankn_svar), intent(inout) :: this
      real(wp), dimension(:), intent(in)        :: new_coeffs
    end subroutine setcoeff
    ! An inspector for the v_matrix component
    pure function getvmat(this) result(curr_v_matrix)
      import wp
      import integral_rankn_svar
      ! Arguments
      class(integral_rankn_svar), intent(in)  :: this
      ! Function result
      real(wp), dimension(: , :), allocatable :: curr_v_matrix
    end function getvmat
    ! A setter for the v_matrix component
    pure subroutine setvmat(this, new_v_matrix)
      import wp
      import integral_rankn_svar
      ! Arguments
      class(integral_rankn_svar), intent(inout) :: this
      real(wp), dimension(: , :), intent(in)    :: new_v_matrix
    end subroutine setvmat
  end interface

contains
  
end module integral_rankn_svar_mod