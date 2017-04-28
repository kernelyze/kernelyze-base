! remez_step.f90
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
! created on: 2015-02-08
! updated on: 2015-04-22
! updated on: 2016-04-24 (Avoid use of internal procedures 
!             for thread safety, both here and upstream)
! updated on: 2016-06-30 (Add optional left and right endpoints
!             for the intervals on which optimization occurs)
! updated on: 2016-07-01 (Modifications to make generalization
!             to asymmetric cases easier; pack interval endpoints
!             into the kernel [this is the kernel domain, after
!             all] and make the is_over_x flag non-optional)
! updated on: 2017-04-22 (use tk_simple_lapack)
!
! A module containing the subroutine remez_step.
!
! NOTE that the kernel $K \left( x , y \right) $ must be
! nondegenerate totally positive on the rectangle
! $\left[a , b \right] \times \left[ c , d \right]$
! to use the subroutine remez_step.  The interval endpoints
! default to a = c = -1 and b = d = 1 if not provided
! (they are optional arguments).
  
module remez_step_mod
  
  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use find_all_zeros_mod, only : find_all_zeros
  use find_rel_optima_mod, only : find_rel_optima
  use kernel_mod, only : kernel
  use linear_combo_kernel_mod, only : linear_combo_kernel
  use tk_simple_lapack_mod, only : tk_gesvd
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  
  implicit none
  
  contains
  
  ! The subroutine to calculate a Remez step
  subroutine remez_step( &
      kernel_obj, &
      rho_vec, &
      coeffs, &
      tolerance, &
      a_matrix_in, &
      grid, &
      nodes, &
      values_at_nodes, &
      is_over_x, &
      err_stat, &
      err_msg)
    ! Arguments
    class(kernel), intent(in)                       :: kernel_obj
    real(wp), intent(in), dimension(:)              :: rho_vec
    real(wp), intent(in), dimension(size(rho_vec))  :: coeffs
    real(wp), intent(in)                            :: tolerance
    real(wp), optional, intent(in), &
        dimension(size(rho_vec) , size(rho_vec))    :: a_matrix_in
    real(wp), intent(in), dimension(:)              :: grid
    real(wp), intent(out), dimension(size(rho_vec)) :: nodes
    real(wp), intent(out), dimension(size(rho_vec)) :: values_at_nodes
    ! If true, the Remez step is over the x variable; otherwise, it
    ! is over the y variable ('is over' means that the sup in the
    ! Remez step is taken over that variable).
    logical, intent(in)                             :: is_over_x
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables
    type(linear_combo_kernel)                       :: deviation
    real(wp), dimension(:), allocatable             :: deviation_zeros
    real(wp), dimension(:), allocatable             :: finer_grid
    real(wp), dimension(:), allocatable             :: optima_grid
    real(wp), dimension(:,:), allocatable           :: result_optima
    real(wp), dimension(size(grid))                 :: check_grid_points
    real(wp), &
          dimension(size(rho_vec) , size(rho_vec))  :: a_matrix, a_matrix_copy
    real(wp)                                        :: condition_number
    integer                                         :: j, curr_num_grid_pts
    ! Local variables to handle the SVD of a_matrix
    integer                                         :: info
    real(wp), dimension(size(rho_vec))              :: svec
    real(wp), &
        dimension(size(rho_vec) , size(rho_vec))    :: umat
    real(wp), &
        dimension(size(rho_vec) , size(rho_vec))    :: vmat_t
    real(wp), dimension(size(rho_vec) - 1)          :: wwvec
    ! Local variables to handle any allocation problems:
    character(*), parameter         :: proc_name = "Remez step: "
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
    ! If an a_matrix_in was supplied, use it
    if (present(a_matrix_in)) then
      a_matrix = a_matrix_in
    ! If no a_matrix_in was given, then set
    ! a_matrix to be the identity matrix
    else
      a_matrix = 0E0_wp
      do j = 1 , size(rho_vec)
        a_matrix(j , j) = 1E0_wp
      end do
    end if
    ! Call the interface tk_gesvd to the LAPACK routine dgesvd to
    ! compute the singular value decomposition of the matrix
    ! a_matrix (pass a copy to avoid interfering with a_matrix).
    ! This is not necessary if a_matrix_in was not supplied,
    ! since in that case the condition_number value will be 1
    ! (since a_matrix will be set to the identity).
    if (present(a_matrix_in)) then
      a_matrix_copy = a_matrix
      call tk_gesvd(  a = a_matrix_copy , &
                      s = svec , &
                      u = umat , &
                      vt = vmat_t , &
                      ww = wwvec , &
                      info = info)
      if (info /= 0) then
        ! This means that the execution of the SVD
        ! routine did not go well.
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = info
        end if
        if (present(err_msg)) then
          err_msg = proc_name // 'Problem in computing SVD'
        end if
        nodes = ieee_value(nodes, ieee_quiet_nan)
        values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
        ! Return
        return
      end if
      ! The absolute value below is used to ensure that numerical
      ! error does not cause the condition number to be -\infty,
      ! since that would cause checks below to fail 
      ! inappropriately.
      condition_number = abs( maxval(svec) / minval(svec) )
    else
      condition_number = 1E0_wp
    end if
    ! Initialize the deviation object -- I only need to do this once,
    ! since the kernel, rho vector, A matrix, and coefficients are all
    ! constant through the remainder of the execution.
    call deviation%init_combo( &
        kernel_obj, &
        rho_vec, &
        a_matrix, &
        coeffs, &
        is_over_x)
    ! Now find all zeros
    call find_all_zeros(deviation_zeros, deviation, grid, &
        tolerance, local_err_stat, local_err_msg)
    if (local_err_stat /= 0) then
      ! This means that finding all zeros did not go well.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      nodes = ieee_value(nodes, ieee_quiet_nan)
      values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
      ! Return
      return
    end if
    ! Then execute numerical checks
    do j = 1 , size(grid)
      check_grid_points(j) = deviation%eval(grid(j))
    end do
    curr_num_grid_pts = size(grid)
    do while (size(deviation_zeros) < (size(rho_vec) - 1))
      curr_num_grid_pts = 2 * curr_num_grid_pts
      if (curr_num_grid_pts > 16000) then
        if (present(err_stat)) then
          err_stat = -1
        end if
        if (present(err_msg)) then
          err_msg = proc_name // 'The search grid has grown too large.'
        end if
        nodes = ieee_value(nodes, ieee_quiet_nan)
        values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
        return
      end if
      if (allocated(finer_grid)) then
        deallocate(finer_grid, stat = alloc_stat, errmsg = alloc_err_msg)
        if (alloc_stat /= 0) then
          ! The deallocation did not complete succesfully
          ! If err_stat and / or err_msg were provided, fill them in
          if (present(err_stat)) then
            err_stat = alloc_stat
          end if
          if (present(err_msg)) then
            err_msg = proc_name // alloc_err_msg
          end if
          ! Set the intent(out) variables to IEEE quiet NaNs
          nodes = ieee_value(nodes, ieee_quiet_nan)
          values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
          ! Return
          return
        end if
      end if
      allocate(finer_grid(curr_num_grid_pts), stat = alloc_stat, errmsg = alloc_err_msg)
      if (alloc_stat /= 0) then
        ! The allocation did not complete succesfully
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = alloc_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // alloc_err_msg
        end if
        ! Set the intent(out) variables to IEEE quiet NaNs
        nodes = ieee_value(nodes, ieee_quiet_nan)
        values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
        ! Return
        return
      end if
      finer_grid = kernel_obj%cheb_pts(is_over_x, curr_num_grid_pts)
      ! Find all zeros
      call find_all_zeros(deviation_zeros, deviation, finer_grid, &
          tolerance, local_err_stat, local_err_msg)
      if (local_err_stat /= 0) then
        ! This means that finding all zeros did not go well.
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = local_err_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // local_err_msg
        end if
        nodes = ieee_value(nodes, ieee_quiet_nan)
        values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
        ! Return
        return
      end if
    end do
    if (size(deviation_zeros) >= size(rho_vec)) then
      if (maxval(abs(check_grid_points)) >= &
          (3E0_wp * condition_number * &
              kernel_obj%eval(1E0_wp, 1E0_wp) * epsilon(1.0E0_wp))) then
        ! Here use a mask to check if the (larger than expected) number of
        ! zeros is a result of the entire function being within tolerance 
        ! of zero.  If so, try searching between the elements of rho_vec.
        ! If even that does not generate a satisfactory result, generate an
        ! error and return
        if (count(abs(check_grid_points) < tolerance) >= size(rho_vec)) then
          ! Find all zeros
          call find_all_zeros(deviation_zeros, deviation, rho_vec, &
              tolerance, local_err_stat, local_err_msg)
          if (local_err_stat /= 0) then
            ! This means that finding all zeros did not go well.
            ! If err_stat and / or err_msg were provided, fill them in
            if (present(err_stat)) then
              err_stat = local_err_stat
            end if
            if (present(err_msg)) then
              err_msg = proc_name // local_err_msg
            end if
            nodes = ieee_value(nodes, ieee_quiet_nan)
            values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
            ! Return
            return
          end if
          if (size(deviation_zeros) /= (size(rho_vec) - 1)) then
            ! This would be very unexpected, and is an error.
            ! If err_stat and / or err_msg were provided, fill them in
            if (present(err_stat)) then
              err_stat = -1
            end if
            if (present(err_msg)) then
              err_msg = proc_name // 'It should not be possible to find this many zeros'
            end if
            nodes = ieee_value(nodes, ieee_quiet_nan)
            values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
            return
          end if
        else
          ! There are more zeros than expected, and they are *not*
          ! the result of the entire function being within tolerance of
          ! zero.  This is an error.
          ! If err_stat and / or err_msg were provided, fill them in
          if (present(err_stat)) then
            err_stat = -1
          end if
          if (present(err_msg)) then
            err_msg = proc_name // 'There are too many zeros, and ' // &
            'not because the entire function is almost zero'
          end if
          nodes = ieee_value(nodes, ieee_quiet_nan)
          values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
          return
        end if
      end if
    end if
    ! Even if the number of zeros is *not* too large, I will
    ! return if all of the grid points deliver a value that is
    ! effectively zero.
    if (maxval(abs(check_grid_points)) < &
        (3E0_wp * condition_number &
        * kernel_obj%eval(1E0_wp, 1E0_wp) * epsilon(1.0E0_wp))) then
      nodes = rho_vec
      do j = 1 , size(rho_vec)
        values_at_nodes(j) = deviation%eval(rho_vec(j))
      end do
      return
    end if
    ! Finally, find all relative optima between the zeros
    ! The grid between whose points I seek optima is formed by
    ! supplementing the deviation_zeros array with a -1.0 before
    ! all other elements and a 1.0 after all other elements
    allocate(optima_grid(size(deviation_zeros) + 2), &
        stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! The allocation did not complete succesfully
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set the intent(out) variables to IEEE quiet NaNs
      nodes = ieee_value(nodes, ieee_quiet_nan)
      values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
      ! Return
      return
    end if
    optima_grid(2:(size(deviation_zeros) + 1)) = deviation_zeros
    if (is_over_x) then
      optima_grid(1)                    = kernel_obj%get_x_lb()
      optima_grid(ubound(optima_grid))  = kernel_obj%get_x_ub()
    else
      optima_grid(1)                    = kernel_obj%get_y_lb()
      optima_grid(ubound(optima_grid))  = kernel_obj%get_y_ub()
    end if
    ! Find the relative optima
    call find_rel_optima(result_optima, deviation, optima_grid, &
        tolerance, local_err_stat, local_err_msg)
    if (local_err_stat /= 0) then
      ! This means that finding relative optima did not go well.
      ! If err_stat and / or err_msg were provided, fill them in
      if (present(err_stat)) then
        err_stat = local_err_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // local_err_msg
      end if
      nodes = ieee_value(nodes, ieee_quiet_nan)
      values_at_nodes = ieee_value(values_at_nodes, ieee_quiet_nan)
      ! Return
      return
    end if
    nodes = result_optima(:,1)
    values_at_nodes = result_optima(:,2)
  end subroutine remez_step
end module remez_step_mod