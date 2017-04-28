! integral_nsect_disc.f90
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
! created on: 2015-07-07
! updated on: 2015-07-07
! updated on: 2016-02-07 (moved get_effective_eval_points)
! updated on: 2016-08-22 (adjust for greater flexibility in
!             template_integralrankn.f90)
!
! This module contains a derived type extending the derived
! type integral_kernel_disc. This extension restricts 
! to the case of a rank-$n$ kernel, and even more, to
! the case in which the functions making up the rank-$n$
! kernel are themselves kernel sections (so, of type
! compfunc_kernel).
  
module integral_nsect_disc_mod

  use compfunc_kernel_mod, only : singlevar_type => compfunc_kernel
  use integral_kernel_disc_mod, only : parent_type => integral_kernel_disc
  
  ! Use the preprocessor to substitute appropriately
  ! for the name of the derived type.
#define integral_rankn integral_rankn_disc
  
  ! It is important to use the preprocessor's inclusion
  ! directive rather than a Fortran include statement
  ! because use of a Fortran include statement would not
  ! lead to the preprocessor define directives getting
  ! applied to the included file.  I need the define
  ! directives to be applied to this file, of course.
#include "template_integralrankn.f90"
  
  ! Effective eval_pts
  pure function get_effective_eval_pts(this) result(curr_eff_eval_pts)
    ! Arguments
    class(integral_rankn_disc), intent(in)  :: this
    ! Result
    real(wp), dimension(:), allocatable     :: curr_eff_eval_pts
    ! Local variables
    integer :: i
    ! Body
    curr_eff_eval_pts = this%v_matrix(:,1) ! A hack, automates allocation
    do i = 1, size(this%comp_array)
      curr_eff_eval_pts(i) = this%comp_array(i)%get_eval_pt()
    end do
  end function get_effective_eval_pts

end module integral_nsect_disc_mod