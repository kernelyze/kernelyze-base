! integral_neigen_disc.f90
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
! created on: 2016-12-04 (simple instatiation, just
!             minor changes and updates relative to
!             the code of integral_nsect_disc on which
!             it is based)
!
! This module contains a derived type extending the derived
! type integral_kernel_disc. This extension restricts 
! to the case of a rank-$n$ kernel, and even more, to
! the case in which the functions making up the rank-$n$
! kernel are linear combinations of  kernel sections 
! (so, of type integral_kernel_disc).  The name of
! the module comes from the fact that this type of
! rank-$n$ kernel is interesting primarily because it
! can be used to model a truncated eigenfunction
! decomposition of a kernel, where the eigenfunctions
! are found numerically using Nystrom's method, see
! kernel_neigen.f90.
  
module integral_neigen_disc_mod

  use integral_kernel_disc_mod, only : singlevar_type => integral_kernel_disc
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
    integer               :: i, j, eval_pts_ind, num_pts
    real(wp), allocatable :: temp_array( : )
    ! Body
    num_pts = 0
    do i = 1, size(this%comp_array)
      num_pts = num_pts + size(this%comp_array(i)%get_eval_pts())
    end do
    allocate(curr_eff_eval_pts(num_pts))
    eval_pts_ind = 1
    do i = 1, size(this%comp_array)
      temp_array = this%comp_array(i)%get_eval_pts()
      do j = 1, size(temp_array)
        curr_eff_eval_pts(eval_pts_ind) = temp_array(j)
        eval_pts_ind = eval_pts_ind + 1
      end do
    end do
  end function get_effective_eval_pts

end module integral_neigen_disc_mod