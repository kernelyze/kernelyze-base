! test_num_opt_rankn_params.f90
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
! created on: 2016-07-05 
! updated on: 2016-07-05
!
! A module of unit tests for the num_opt_rankn_params module.
! This is notably different from num_opt_rankn_kernel because:
!
! 1) The kernel to approximate does ***not*** need to be symmetric.
! 2) The output is composed of the numerical parameters needed to
!    evaluate the numerically-optimal rank-$n$ approximation; no
!    approximate-kernel object is built.
  
module test_num_opt_rankn_params_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use num_opt_rankn_params_mod, only : num_opt_rankn_params

implicit none

contains

function test_num_opt_rankn_params(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of all num-opt-rankn-params ' &
      // 'functionality: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_num_opt_rankn_params_exp_prod(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_params_gaussian(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_params_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_num_opt_rankn_params

function test_num_opt_rankn_params_exp_prod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_expprod)          :: exp_prod_kernel
  real(wp), dimension(9 , 8)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(9)        :: test_rho, test_gamma
  real(wp), parameter           :: num_toler = 1E-14_wp
  real(wp)                      :: b_lowerbd
  integer                       :: i, j
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp)                      :: rankn_eval, kernel_eval
  real(wp)                      :: max_abs_err
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt-rankn-params for the ' &
      // 'exponential product kernel: '
  call exp_prod_kernel%set_c(1E0_wp)
  call exp_prod_kernel%set_x_lb(0E0_wp)
  call exp_prod_kernel%set_x_ub(5E0_wp)
  call exp_prod_kernel%set_y_lb(-1E0_wp)
  call exp_prod_kernel%set_y_ub(0.5E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call num_opt_rankn_params( &
      kernel_to_approx = exp_prod_kernel, &
      rank = size(test_rho) - 1, &
      max_iter = 100, &
      toler = 1E-15_wp, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      borsuk_lb = b_lowerbd, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the rho vector
  write(unit_to_use , *) 'The computed rho vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_rho
  ! Output the gamma vector
  write(unit_to_use , *) 'The computed gamma vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_gamma
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Now check the results using the Borsuk lower bound:
  test_grid_x = exp_prod_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = exp_prod_kernel%cheb_pts(.false., size(test_grid_y))
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = 1, size(test_grid_x)
    do j = 1, size(test_grid_y)
      rankn_eval = dot_product( &
              matmul( exp_prod_kernel%eval_elt(test_rho , test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( exp_prod_kernel%eval_elt(test_gamma , test_grid_y(j)) , &
                  test_w_matrix ) )
      kernel_eval = exp_prod_kernel%eval(test_grid_x(i), test_grid_y(j))
      max_abs_err = max( max_abs_err , abs(rankn_eval - kernel_eval) )
      if (abs(rankn_eval - kernel_eval) > (b_lowerbd + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            b_lowerbd + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - kernel_eval)
      end if
    end do
  end do
  ! Output key summary results:
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write( unit_to_use , * ) 'Largest abs val err on test grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_num_opt_rankn_params_exp_prod

function test_num_opt_rankn_params_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_gaussian)         :: gaussian_kernel
  real(wp), dimension(6 , 5)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(6)        :: test_rho, test_gamma
  real(wp), parameter           :: num_toler = 1E-14_wp
  real(wp)                      :: b_lowerbd
  integer                       :: i, j
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp)                      :: rankn_eval, kernel_eval
  real(wp)                      :: max_abs_err
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt-rankn-params for the ' &
      // 'Gaussian kernel: '
  call gaussian_kernel%set_c(1E0_wp)
  call gaussian_kernel%set_x_lb(-4E0_wp)
  call gaussian_kernel%set_x_ub(-1E0_wp)
  call gaussian_kernel%set_y_lb(-2E0_wp)
  call gaussian_kernel%set_y_ub(2E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call num_opt_rankn_params( &
      kernel_to_approx = gaussian_kernel, &
      rank = size(test_rho) - 1, &
      max_iter = 100, &
      toler = 1E-15_wp, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      borsuk_lb = b_lowerbd, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the rho vector
  write(unit_to_use , *) 'The computed rho vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_rho
  ! Output the gamma vector
  write(unit_to_use , *) 'The computed gamma vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_gamma
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Now check the results using the Borsuk lower bound:
  test_grid_x = gaussian_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = gaussian_kernel%cheb_pts(.false., size(test_grid_y))
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = 1, size(test_grid_x)
    do j = 1, size(test_grid_y)
      rankn_eval = dot_product( &
              matmul( gaussian_kernel%eval_elt(test_rho , test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( gaussian_kernel%eval_elt(test_gamma , test_grid_y(j)) , &
                  test_w_matrix ) )
      kernel_eval = gaussian_kernel%eval(test_grid_x(i), test_grid_y(j))
      max_abs_err = max( max_abs_err , abs(rankn_eval - kernel_eval) )
      if (abs(rankn_eval - kernel_eval) > (b_lowerbd + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            b_lowerbd + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - kernel_eval)
      end if
    end do
  end do
  ! Output key summary results:
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write( unit_to_use , * ) 'Largest abs val err on test grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_num_opt_rankn_params_gaussian

function test_num_opt_rankn_params_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Local variables
  type(kernel_cauchy)           :: cauchy_kernel
  real(wp), dimension(8 , 7)    :: test_v_matrix, test_w_matrix
  real(wp), dimension(8)        :: test_rho, test_gamma
  real(wp), parameter           :: num_toler = 1E-14_wp
  real(wp)                      :: b_lowerbd
  integer                       :: i, j
  real(wp), dimension(201)      :: test_grid_x, test_grid_y
  real(wp)                      :: rankn_eval, kernel_eval
  real(wp)                      :: max_abs_err
  integer                       :: err_stat
  character(len=err_msg_len)    :: err_msg
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt-rankn-params for the ' &
      // 'Cauchy kernel: '
  call cauchy_kernel%set_c(3E0_wp)
  call cauchy_kernel%set_alpha(1E0_wp)
  call cauchy_kernel%set_x_lb(0E0_wp)
  call cauchy_kernel%set_x_ub(5E0_wp)
  call cauchy_kernel%set_y_lb(-0.5E0_wp)
  call cauchy_kernel%set_y_ub(2E0_wp)
  ! Get the worst rho and gamma vectors for the kernel
  call num_opt_rankn_params( &
      kernel_to_approx = cauchy_kernel, &
      rank = size(test_rho) - 1, &
      max_iter = 100, &
      toler = 1E-15_wp, &
      rho_vec = test_rho, &
      gamma_vec = test_gamma, &
      v_mat = test_v_matrix, &
      w_mat = test_w_matrix, &
      borsuk_lb = b_lowerbd, &
      err_stat = err_stat, &
      err_msg = err_msg)
  if (err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use, *) 'Test failed with message ', err_msg
    return
  end if
  ! Output the rho vector
  write(unit_to_use , *) 'The computed rho vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_rho
  ! Output the gamma vector
  write(unit_to_use , *) 'The computed gamma vector is: '
  ! Use an unlimited format to control output
  write(unit_to_use , '(*(2x,E15.8))') test_gamma
  ! Output the V matrix
  write(unit_to_use , *) 'The computed V matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_v_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_v_matrix(i , :)
  end do
  ! Output the W matrix
  write(unit_to_use , *) 'The computed W matrix is: '
  ! Use an unlimited format to control output
  do i = 1, size(test_w_matrix, 1)
    write(unit_to_use , '(*(2x,E15.8))') test_w_matrix(i , :)
  end do
  ! Now check the results using the Borsuk lower bound:
  test_grid_x = cauchy_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = cauchy_kernel%cheb_pts(.false., size(test_grid_y))
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = 1, size(test_grid_x)
    do j = 1, size(test_grid_y)
      rankn_eval = dot_product( &
              matmul( cauchy_kernel%eval_elt(test_rho , test_grid_x(i)) , &
                  test_v_matrix ) , &
              matmul( cauchy_kernel%eval_elt(test_gamma , test_grid_y(j)) , &
                  test_w_matrix ) )
      kernel_eval = cauchy_kernel%eval(test_grid_x(i), test_grid_y(j))
      max_abs_err = max( max_abs_err , abs(rankn_eval - kernel_eval) )
      if (abs(rankn_eval - kernel_eval) > (b_lowerbd + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            b_lowerbd + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - kernel_eval)
      end if
    end do
  end do
  ! Output key summary results:
  write(unit_to_use , *) 'Test grid over x runs from ', &
      minval(test_grid_x), ' to ', maxval(test_grid_x)
  write(unit_to_use , *) 'Test grid over y runs from ', &
      minval(test_grid_y), ' to ', maxval(test_grid_y)
  write(unit_to_use , *) 'Borsuk lower bound is: ' , b_lowerbd
  write( unit_to_use , * ) 'Largest abs val err on test grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_num_opt_rankn_params_cauchy

end module test_num_opt_rankn_params_mod