! eigen_rankn_kernel_of_params.f90
!
! Copyright (c) 2017 by Kernelyze LLC
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
! created on: 2017-04-20 (based on earlier eigen_rankn_kernel)
! updated on: 2017-04-21 (fixed x-y confusion in a compfunc)
!
! This module contains a constructor subroutine
! that builds a rank-$n$ truncated eigenfunction
! series approximation, given parameters already
! computed for that purpose.  See eigen_rankn_params.f90
! for a numerical subroutine that computes the
! relevant parameters.
  
module eigen_rankn_kernel_of_params_mod

  use set_precision, only : wp
  use constants_mod, only : err_msg_len
  use kernel_mod, only : kernel
  use kernel_neigen_mod, only : kernel_rankn
  use integral_kernel_disc_mod, only : integral_kernel_disc

  implicit none

  contains
  
  subroutine eigen_rankn_kernel_of_params( &
      kernel_to_approx, &
      x_points, &
      y_points, &
      trunc_sing_vals, &
      x_func_wts, &
      y_func_wts, &
      rankn_kernel, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)     :: kernel_to_approx
    real(wp), intent(in)          :: x_points( : )
    real(wp), intent(in)          :: y_points( : )
    real(wp), intent(in)          :: trunc_sing_vals( : )
    real(wp), intent(in)          :: x_func_wts( &
        size(y_points) , size(trunc_sing_vals) )
    real(wp), intent(in)          :: y_func_wts( &
        size(y_points) , size(trunc_sing_vals) )
    class(kernel_rankn), allocatable, intent(out) :: rankn_kernel
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables to handle information on any local errors
    integer                         :: local_err_stat
    character(len=err_msg_len)      :: local_err_msg
    ! Local parameter to label any error messages:
    character(*), parameter :: proc_name = 'eigen_rankn_kernel_of_params: '
    ! Local variables
    type(integral_kernel_disc)  :: x_compfuncs( size(trunc_sing_vals) )
    type(integral_kernel_disc)  :: y_compfuncs( size(trunc_sing_vals) )
    real(wp)                    :: diag_matrix( &
        size(trunc_sing_vals), size(trunc_sing_vals) )
    integer                     :: i, j
    integer                     :: rank
    ! Body
    ! Initialize the optional error arguments if present
    if (present(err_stat)) then
      err_stat = 0
    end if
    if (present(err_msg)) then
      err_msg = ''
    end if
    ! The rank is the number of singular values provided
    rank = size( trunc_sing_vals )
    ! Check that the rank is <= min( size(x_points), size(y_points) )
    if (rank > size(x_points) .or. rank > size(y_points) ) then
      if (present(err_stat)) then
        err_stat = -1
      end if
      if (present(err_msg)) then
        err_msg = proc_name // 'rank must be <= min size of points'
      end if
      return
    end if
    ! Allocate the rankn_kernel
    allocate(rankn_kernel, stat = local_err_stat, errmsg = local_err_msg)
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
    ! Set up the components of the rankn_kernel
    ! Allocate and build the component functions
    do j = 1, size(x_compfuncs)
      call x_compfuncs(j)%set_kernel( &
          kernel_to_approx, &
          local_err_stat, &
          local_err_msg)
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
      call x_compfuncs(j)%set_integral_vbl( .false. ) ! Integral over y
      call x_compfuncs(j)%set_eval_pts( y_points )
      ! Note that the functions of x are weighted sums over y points;
      ! this is a feature of the Nystrom method.
      call x_compfuncs(j)%set_weights( x_func_wts( : , j ) )
    end do
    call rankn_kernel%set_f_array(x_compfuncs, local_err_stat, local_err_msg)
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
    do j = 1, size(y_compfuncs)
      call y_compfuncs(j)%set_kernel( &
          kernel_to_approx, &
          local_err_stat, &
          local_err_msg)
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
      call y_compfuncs(j)%set_integral_vbl( .true. ) ! Integral over x
      call y_compfuncs(j)%set_eval_pts( x_points )
      ! Note that the functions of y are weighted sums over x points;
      ! this is a feature of the Nystrom method.
      call y_compfuncs(j)%set_weights( y_func_wts( : , j ) )
    end do
    call rankn_kernel%set_g_array(y_compfuncs, local_err_stat, local_err_msg)
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
    diag_matrix = 0E0_wp
    do i = 1, rank
      diag_matrix( i , i ) = sqrt( trunc_sing_vals( i ) )  
    end do
    call rankn_kernel%set_v_matrix(diag_matrix)
    call rankn_kernel%set_w_matrix(diag_matrix)
  end subroutine eigen_rankn_kernel_of_params

end module eigen_rankn_kernel_of_params_mod