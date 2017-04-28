! remez_variant.f90
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
! created on: 2015-03-07
! updated on: 2015-04-23
! updated on: 2016-04-24 (Avoid upstream use of internal procedures 
!             for thread safety)
! updated on: 2016-07-01 (Allow for different rectangular domains
!             and operation over either x or y [for asymmetry])
! updated on: 2016-07-03 (Correctly handle the A matrix when 
!             vec_for_dot_prod is given, requiring it to be
!             passed in when I am working with a potentially
!             asymmetric kernel; also, no need to fail if
!             the dot product of coeffs with vec_for_dot_prod
!             is negative, the scaling will be fine in that
!             case)
! updated on: 2017-04-22 (use tk_simple_lapack, correct comment)
!
! A module containing the subroutine remez_variant.
! This is the key functionality underlying both 
! borsuk_lower_bound and kernel_remez.
! Given a set of points in rho_vec, this subroutine
! minimizes the maximum deviation from zero (if over x)
! $max_{\left[x_L, x_U \right]} 
! \left|coeffs' * func_matrix * vec_of_kernel(x) \right|$ or
! (if over y)
! $max_{\left[y_L, y_U \right]} 
! \left|coeffs' * func_matrix * vec_of_kernel(y) \right|$
! The term 'vec_of_kernel' is a vector of kernel evaluations, 
! $K \left( x, \rho_j \right)$ if over x or 
! $K \left( \gamma_j ,  y \right)$ if over y.
! The term 'func_matrix' is the identity for the Borsuk lower
! bound or a matrix based on the kernel and the \rho vector
! for the kernel Remez method.
! The 'coeffs' term is a vector of coefficients to be solved for 
! (with the constraint that they must have: unit L^2 norm 
! for the kernel Remez method, unit L^1 norm for the Borsuk 
! lower bound, and unit inner product with a fixed vector for the
! asymmetric kernel Remez method) to minimize the maximum 
! deviation from zero.
!
! NOTE that the kernel $K \left( x , y \right) $ must be
! strictly sign regular on the rectangle
! $\left[x_L , x_U \right] \times \left[ y_L , y_U \right]$
! to use the subroutine borsuk_lower_bound.
  
module remez_variant_mod
  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len
  use chebyshev_points_mod, only : chebyshev_points
  use compute_a_matrix_mod, only : compute_a_matrix
  use moore_penrose_inverse_mod, only : moore_penrose_inverse
  use remez_step_mod, only : remez_step
  use kernel_mod, only : kernel
  use tk_simple_lapack_mod, only : tk_gesv
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  contains
  subroutine remez_variant( &
      is_borsuk, &
      kernel_obj, &
      is_over_x, &
      rho_vec, &
      tolerance, &
      max_iter, &
      vec_for_dot_prod, &
      key_matrix, &
      grid, &
      coeffs, &
      nodes, &
      errors_at_nodes, &
      discrepancy, &
      num_iter, &
      err_stat, &
      err_msg)
    ! Arguments
    ! The first argument determines whether to use Borsuk lower bound
    ! or kernel Remez values.
    logical, intent(in)                               :: is_borsuk
    class(kernel), intent(in)                         :: kernel_obj
    logical, intent(in)                               :: is_over_x
    real(wp), intent(in), dimension(:)                :: rho_vec
    real(wp), intent(in)                              :: tolerance
    integer, intent(in)                               :: max_iter
    real(wp), optional, intent(in), &
      dimension(size(rho_vec))                        :: vec_for_dot_prod
    real(wp), optional, intent(in), &
      dimension(size(rho_vec) , size(rho_vec))        :: key_matrix
    real(wp), optional, intent(in), dimension(:)      :: grid
    real(wp), intent(out), dimension(size(rho_vec))   :: coeffs
    real(wp), intent(out), dimension(size(rho_vec))   :: nodes
    real(wp), intent(out), dimension(size(rho_vec))   :: errors_at_nodes
    real(wp), intent(out)                             :: discrepancy
    integer, intent(out)                              :: num_iter
    ! Optional arguments to pass back information on any errors
    integer, intent(out), optional                    :: err_stat
    character(len=err_msg_len), intent(out), optional :: err_msg
    ! Local variables
    real(wp), dimension(:), allocatable               :: grid_to_use
    real(wp), dimension(size(rho_vec))                :: target_vec
    real(wp), dimension(size(rho_vec))                :: new_nodes
    real(wp), dimension(size(rho_vec))                :: new_coeffs
    real(wp), dimension(size(rho_vec))                :: new_errors_at_nodes
    integer, dimension(size(rho_vec))                 :: ipiv
    real(wp), dimension(size(rho_vec),size(rho_vec))  :: func_mat, func_mat_copy
    real(wp), dimension(size(rho_vec),size(rho_vec))  :: a_matrix
    real(wp)                                          :: old_max
    logical                                           :: rho_has_duplicates
    integer                                           :: i, j, k, lapack_info
    ! Local variables to handle any allocation problems:
    character(*), parameter         :: proc_name = "Remez variant: "
    character(len=alloc_errmsg_len) :: alloc_err_msg
    integer                         :: alloc_stat
    ! Local variables to handle any local errors
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
    ! Initialize the number of iterations to zero
    num_iter = 0
    ! If this is a call for a Borsuk lower bound, the
    ! a_matrix is the identity.  Otherwise, it is a call
    ! for a kernel Remez method and the a_matrix must be
    ! computed.
    if (is_borsuk) then
      ! Set the a_matrix to the identity matrix
      a_matrix = 0E0_wp
      do j = 1 , size(rho_vec)
        a_matrix(j, j) = 1E0_wp
      end do
      ! There should be no vec_for_dot_prod if this is a
      ! Borsuk lower bound call -- the vec_for_dot_prod is
      ! used (if present) to normalize coefficients in a
      ! kernel Remez call.
      if (present(vec_for_dot_prod)) then
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = -11
        end if
        if (present(err_msg)) then
          err_msg = proc_name // 'vector for dot product should not ' &
              // 'be passed for a Borsuk lower bound call.'
        end if
        ! Set all intent(out) variables to IEEE quiet NaNs
        coeffs = ieee_value(coeffs, ieee_quiet_nan)
        nodes = ieee_value(nodes, ieee_quiet_nan)
        errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
        discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
        ! Return
        return
      end if
    else
      ! Calculate the a_matrix
      if (present(vec_for_dot_prod)) then    
        if (.not. present(key_matrix)) then
          ! If err_stat and / or err_msg were provided, fill them in
          if (present(err_stat)) then
            err_stat = -11
          end if
          if (present(err_msg)) then
            err_msg = proc_name // 'if vec_for_dot_prod is given, then ' &
                // 'key_matrix must also be provided.'
          end if
          ! Set all intent(out) variables to IEEE quiet NaNs
          coeffs = ieee_value(coeffs, ieee_quiet_nan)
          nodes = ieee_value(nodes, ieee_quiet_nan)
          errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
          discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
          ! Return
          return
        end if
        a_matrix = key_matrix
      else
        a_matrix = compute_a_matrix(kernel_obj , rho_vec)
      end if
    end if
    ! Handle the grid, since this is an optional argument;
    ! using 200 Chebyshev points will work well across a broad
    ! range of use cases.
    if (present(grid)) then
      allocate(grid_to_use(size(grid)), stat = alloc_stat, errmsg = alloc_err_msg)
      ! What to do if there is an allocation error:
      if (alloc_stat /= 0) then
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = alloc_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // alloc_err_msg
        end if
        ! Set all intent(out) variables to IEEE quiet NaNs
        coeffs = ieee_value(coeffs, ieee_quiet_nan)
        nodes = ieee_value(nodes, ieee_quiet_nan)
        errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
        discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
        ! Return
        return
      end if
      grid_to_use = grid
    else
      allocate(grid_to_use(200), stat = alloc_stat, errmsg = alloc_err_msg)
      ! What to do if there is an allocation error:
      if (alloc_stat /= 0) then
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = alloc_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // alloc_err_msg
        end if
        ! Set all intent(out) variables to IEEE quiet NaNs
        coeffs = ieee_value(coeffs, ieee_quiet_nan)
        nodes = ieee_value(nodes, ieee_quiet_nan)
        errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
        discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
        ! Return
        return
      end if
      grid_to_use = kernel_obj%cheb_pts(is_over_x, 200)      
    end if
    ! If there are duplicates in the rho vector, handle them.
    ! I don't use the unirnk functionality here, because
    ! if there are duplicates I actually need to find a
    ! pair of them in order to return properly.
    rho_has_duplicates = .false.
    outerduploop: do j = 1,size(rho_vec)
      do k = j + 1, size(rho_vec)
        if (rho_vec(j) == rho_vec(k)) then
          rho_has_duplicates = .true.
          exit outerduploop
        end if
      end do
    end do outerduploop
    ! At this point, if I exited from the above loop then
    ! rho_has_duplicates will be .true., and j and k
    ! will be the indices of two duplicate entries of
    ! rho_vec.  Respond accordingly.
    ! If there are any duplicate entries in the rho
    ! vector, then it is possible to achieve a zero
    ! maximum abs error by setting the coefficient of one
    ! copy of the duplicated entry to the negative of
    ! the coefficient on the other copy of the duplicated
    ! entry (with all other coefficients zero).
    if (rho_has_duplicates) then
      coeffs = 0E0_wp  ! Set all entries of coeff to zero
      coeffs(j) = -0.5E0_wp
      coeffs(k) = 0.5E0_wp
      errors_at_nodes = 0E0_wp ! Set all errors at nodes to zero
      nodes = rho_vec ! The max abs error of zero occurs everywhere
      discrepancy = 0E0_wp
      return
    end if
    ! A decent initial guess as to the nodes where the max abs
    ! error is attained is the Chebyshev points
    nodes = kernel_obj%cheb_pts(is_over_x, size(rho_vec))
    ! Fill in the func_mat
    do j = 1 , size(rho_vec)
      do i = 1 , size(rho_vec)
        if (is_over_x) then
          func_mat(i , j) = kernel_obj%eval(nodes(i) , rho_vec(j))
        else
          func_mat(i , j) = kernel_obj%eval(rho_vec(j) , nodes(i))
        end if
      end do
    end do
    ! Alternating array: -1, 1, -1, 1, etc.
    target_vec = (/ ( (-1.0E0_wp)**j, j = 1, size(rho_vec) ) /)
    ! For a Borsuk lower bound call, I can solve the linear
    ! system directly.  For a kernel Remez method call, I need
    ! to be quite careful in handling the linear system.
    if (is_borsuk) then
      func_mat_copy = func_mat ! This will be overwritten momentarily
      coeffs = target_vec ! This will be overwritten momentarily
      ! The subroutine tk_gesv overwrites func_mat_copy using the L and U
      ! factors in the LU decomposition; it also overwrites
      ! coeffs with the solution of the linear system
      call tk_gesv( func_mat_copy , coeffs , ipiv , lapack_info)
      if ( lapack_info /= 0) then
        ! This means that the execution of the linear system
        ! solver did not go well.
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = lapack_info
        end if
        if (present(err_msg)) then
          err_msg = proc_name // 'Problem in solving linear system'
        end if
        coeffs = ieee_value(coeffs, ieee_quiet_nan)
        nodes = ieee_value(nodes, ieee_quiet_nan)
        errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
        discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
        ! Return
        return
      end if
      ! Normalize the coeffs to have unit L^1 norm.
      coeffs = coeffs / sum(abs(coeffs))
    else
      ! If this is not a Borsuk lower bound call, then
      ! the a_matrix will generally not be the identity
      ! and I need to multiply func_mat by its transpose:
      func_mat = matmul(func_mat , transpose(a_matrix))
      func_mat_copy = func_mat ! This will be overwritten momentarily
      ! Using the Moore-Penrose generalized inverse is very
      ! helpful in handling the possible near-singularity of the
      ! matrix func_mat.
      coeffs = matmul(moore_penrose_inverse(func_mat_copy) , target_vec)
      ! Normalize the coeffs to have unit L^2 norm (if vec_for_dot_prod
      ! is not provided) or to have unit dot product with vec_for_dot_prod
      if (present(vec_for_dot_prod)) then
        coeffs = coeffs / dot_product(coeffs , vec_for_dot_prod)
      else
        coeffs = coeffs / sqrt(dot_product(coeffs , coeffs))
      end if
    end if
    errors_at_nodes = matmul(func_mat,coeffs)
    discrepancy = tolerance + 1.0E0_wp
    ! Note: num_iter is correctly initialized to zero above
    ! Iterate to find the optimum
    do while (discrepancy > tolerance .and. num_iter <= max_iter)
      old_max = maxval(abs(errors_at_nodes))
      ! Use the current coefficients to find the max abs values of
      ! the errors. Note that the "new_errors_at_nodes" will, after
      ! this call, hold the attained node values under the 
      ! *current* coefficients; this will be overwritten
      ! below by the node values under the *new* coefficients.
      call remez_step( &
          kernel_obj = kernel_obj , &
          rho_vec = rho_vec , &
          coeffs = coeffs , &
          tolerance = tolerance , &
          a_matrix_in = a_matrix, & ! This will just be the identity if Borsuk
          grid = grid_to_use , &
          nodes = new_nodes , &
          values_at_nodes = new_errors_at_nodes, &
          is_over_x = is_over_x, &
          err_stat = local_err_stat, &
          err_msg = local_err_msg)
      if (local_err_stat /= 0) then
        ! This means that the Remez step did not go well.
        ! If err_stat and / or err_msg were provided, fill them in
        if (present(err_stat)) then
          err_stat = local_err_stat
        end if
        if (present(err_msg)) then
          err_msg = proc_name // local_err_msg
        end if
        ! Set all intent(out) variables to IEEE quiet NaNs
        coeffs = ieee_value(coeffs, ieee_quiet_nan)
        nodes = ieee_value(nodes, ieee_quiet_nan)
        errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
        discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
        ! Return
        return
      end if
      ! Fill in the func_mat
      do j = 1 , size(rho_vec)
        do i = 1 , size(rho_vec)
          if (is_over_x) then
            func_mat(i , j) = kernel_obj%eval(new_nodes(i) , rho_vec(j))
          else
            func_mat(i , j) = kernel_obj%eval(rho_vec(j) , new_nodes(i))
          end if
        end do
      end do
      ! For a Borsuk lower bound call, I can solve the linear
      ! system directly.  For a kernel Remez method call, I need
      ! to be quite careful in handling the linear system.
      if (is_borsuk) then
        new_coeffs = target_vec ! This will be overwritten momentarily
        func_mat_copy = func_mat ! This will be overwritten momentarily
        ! The subroutine tk_gesv overwrites func_mat_copy using the L and U
        ! factors in the LU decomposition; it also overwrites
        ! new_coeffs with the solution of the linear system
        call tk_gesv( func_mat_copy , new_coeffs , ipiv , lapack_info)
        if ( lapack_info /= 0) then
          ! This means that the execution of the linear system
          ! solver did not go well.
          ! If err_stat and / or err_msg were provided, fill them in
          if (present(err_stat)) then
            err_stat = lapack_info
          end if
          if (present(err_msg)) then
            err_msg = proc_name // 'Problem in solving linear system'
          end if
          ! Set all intent(out) variables to IEEE quiet NaNs
          coeffs = ieee_value(coeffs, ieee_quiet_nan)
          nodes = ieee_value(nodes, ieee_quiet_nan)
          errors_at_nodes = ieee_value(errors_at_nodes, ieee_quiet_nan)
          discrepancy = ieee_value(discrepancy, ieee_quiet_nan)
          ! Return
          return
        end if
        ! Normalize the coeffs to have unit L^1 norm.
        new_coeffs = new_coeffs / sum(abs(new_coeffs))
      else
        ! If this is not a Borsuk lower bound call, then
        ! the a_matrix will generally not be the identity
        ! and I need to multiply func_mat by its transpose:
        func_mat = matmul(func_mat , transpose(a_matrix))
        func_mat_copy = func_mat ! This will be overwritten momentarily
        ! Using the Moore-Penrose generalized inverse is very
        ! helpful in handling the possible near-singularity of the
        ! matrix func_mat.
        new_coeffs = matmul(moore_penrose_inverse(func_mat_copy) , target_vec)
        ! Normalize the coeffs to have unit L^2 norm (if vec_for_dot_prod
        ! is not provided) or to have unit dot product with vec_for_dot_prod
        if (present(vec_for_dot_prod)) then
          new_coeffs = new_coeffs &
              / dot_product(new_coeffs , vec_for_dot_prod)
        else
          new_coeffs = new_coeffs / sqrt(dot_product(new_coeffs , new_coeffs))
        end if
      end if
      ! Calculate the errors at the new nodes
      new_errors_at_nodes = matmul(func_mat,new_coeffs)
      discrepancy = maxval(abs(new_errors_at_nodes)) - old_max
      if (discrepancy < 0E0_wp) then
        discrepancy = 0E0_wp
        exit
      end if
      num_iter = num_iter + 1
      nodes = new_nodes
      coeffs = new_coeffs
      errors_at_nodes = new_errors_at_nodes
    end do
  end subroutine remez_variant
end module remez_variant_mod