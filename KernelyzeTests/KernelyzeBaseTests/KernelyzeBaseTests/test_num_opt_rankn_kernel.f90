! test_num_opt_rankn_kernel.f90
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
! created on: 2015-04-06
! updated on: 2016-01-02
! updated on: 2016-09-27 (added kernel_blackscholes and kernel_bachelier)
! updated on: 2016-09-30 (conform to new argument convention for 
!             kernel_blackscholes: arguments to the kernel are now
!             strike and forward rather than ln(strike) and ln(forward))
! updated on: 2016-10-02 (reverted to previous convention for
!             kernel_blackscholes inputs -- this makes sense for
!             scenario shifts, e. g., parallel shifts in ln(forward)
!             are plausible, while parallel shifts in the forward
!             curve itself are harder to justify economically)
! updated on: 2017-04-23 (remove the symmetric-only num_opt_rankn_kernel
!             tests)
!
! A module of unit tests for the num_opt_rankn_kernel module and
! the num_opt_rankn_kernel_assym module.
  
module test_num_opt_rankn_kernel_mod

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len, err_msg_len
use kernel_mod, only : kernel
use kernel_nsect_mod, only : kernel_rankn
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use kernel_bachelier_mod, only : kernel_bachelier
use kernel_blackscholes_mod, only : kernel_blackscholes
use num_opt_rankn_kernel_asymm_mod, only : num_opt_rankn_kernel_asymm

implicit none

contains

function test_num_opt_rankn_kernel(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality to build ' &
      // 'numerically-optimal rank-n kernels: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_num_opt_rankn_kernel_expprod(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_kernel_gaussian(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_kernel_cauchy(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_kernel_bachelier(unit_to_use)
  test_pass = test_pass .and. test_num_opt_rankn_kernel_black(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_num_opt_rankn_kernel

function test_num_opt_rankn_kernel_expprod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter               :: num_toler = 2E-12_wp
  type(kernel_expprod)              :: expprod_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: borsuk_lb
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, expprod_eval, max_abs_err
  real(wp), dimension( : , : ), &
      allocatable                   :: v_matrix
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt rank-n kernel construction for the ' &
      // 'exponential product kernel: '
  ! Set the "c" parameter of the exponential product kernel
  call expprod_kernel%set_c(1E0_wp)
  ! Call the num_opt_rankn_kernel subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! expprod_kernel.
  call num_opt_rankn_kernel_asymm( &
      expprod_kernel, & ! The kernel to approximate
      8, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the Borsuk lower bound for this approximation
  write(unit_to_use , *) 'The Borsuk lower bound on the supremum of ' &
      // 'the approximation error is: ', borsuk_lb
  ! Get the V matrix of the rank-$n$ approximation and output it:
  v_matrix = rankn_kernel%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n kernel is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 10E0_wp
      test_y = real(j, wp) / 10E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      expprod_eval = expprod_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'rank-n eval: ', rankn_eval, &
          ' kernel eval: ', expprod_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - expprod_eval) )
      if (abs(rankn_eval - expprod_eval) > (borsuk_lb + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            borsuk_lb + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - expprod_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_num_opt_rankn_kernel_expprod

function test_num_opt_rankn_kernel_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter               :: num_toler = 2E-13_wp
  type(kernel_gaussian)             :: gaussian_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: borsuk_lb
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, gaussian_eval, max_abs_err
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt rank-n kernel construction for the ' &
      // 'Gaussian kernel: '
  ! Set the "c" parameter of the Gaussian kernel
  call gaussian_kernel%set_c(1E0_wp)
  ! Call the num_opt_rankn_kernel subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! gaussian_kernel.
  call num_opt_rankn_kernel_asymm( &
      gaussian_kernel, & ! The kernel to approximate
      8, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the Borsuk lower bound for this approximation
  write(unit_to_use , *) 'The Borsuk lower bound on the supremum of ' &
      // 'the approximation error is: ', borsuk_lb
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 10E0_wp
      test_y = real(j, wp) / 10E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      gaussian_eval = gaussian_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'rank-n eval: ', rankn_eval, &
          ' kernel eval: ', gaussian_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - gaussian_eval) )
      if (abs(rankn_eval - gaussian_eval) > (borsuk_lb + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            borsuk_lb + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - gaussian_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_num_opt_rankn_kernel_gaussian

function test_num_opt_rankn_kernel_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter               :: num_toler = 3E-14_wp
  type(kernel_cauchy)               :: cauchy_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: borsuk_lb
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, cauchy_eval, max_abs_err
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt rank-n kernel construction for the ' &
      // 'Cauchy kernel: '
  ! Set the "c" parameter of the Cauchy kernel
  call cauchy_kernel%set_c(3E0_wp)
  ! Set the "alpha" parameter of the Cauchy kernel
  call cauchy_kernel%set_alpha(1E0_wp)
  ! Call the num_opt_rankn_kernel subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! cauchy_kernel.
  call num_opt_rankn_kernel_asymm( &
      cauchy_kernel, & ! The kernel to approximate
      5, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the Borsuk lower bound for this approximation
  write(unit_to_use , *) 'The Borsuk lower bound on the supremum of ' &
      // 'the approximation error is: ', borsuk_lb
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 10E0_wp
      test_y = real(j, wp) / 10E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      cauchy_eval = cauchy_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'rank-n eval: ', rankn_eval, &
          ' kernel eval: ', cauchy_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - cauchy_eval) )
      if (abs(rankn_eval - cauchy_eval) > (borsuk_lb + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            borsuk_lb + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - cauchy_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_num_opt_rankn_kernel_cauchy

function test_num_opt_rankn_kernel_bachelier(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter               :: num_toler = 3E-14_wp
  type(kernel_bachelier)            :: bachelier_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: err_stat
  character(len=err_msg_len)        :: err_msg
  integer                           :: i, j
  real(wp)                          :: borsuk_lb
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, bachelier_eval, max_abs_err
  real(wp)                          :: test_grid_x(11), test_grid_y(11)
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt rank-n kernel construction for the ' &
      // 'Bachelier kernel: '
  ! Set the parameters of the Bachelier kernel
  call bachelier_kernel%set_is_call(.true.)
  call bachelier_kernel%set_sigma(1E-2_wp) ! 100 bps of vol
  call bachelier_kernel%set_x_lb(1E-3_wp) ! 10 bp min strike
  call bachelier_kernel%set_x_ub(3E-2_wp) ! 300 bp max strike
  call bachelier_kernel%set_y_lb(-1E-2_wp) ! -100 bp min forward
  call bachelier_kernel%set_y_ub(5E-2_wp) ! 500 bp max forward
  ! Call the num_opt_rankn_kernel_asymm subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! bachelier_kernel.
  call num_opt_rankn_kernel_asymm( &
      bachelier_kernel, & ! The kernel to approximate
      3, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      err_stat, & ! The error status flag
      err_msg) ! The error message (if any)
  ! Output the Borsuk lower bound for this approximation
  write(unit_to_use , *) 'The Borsuk lower bound on the supremum of ' &
      // 'the approximation error is: ', borsuk_lb
  ! Check some evaluation points
  test_grid_x = bachelier_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = bachelier_kernel%cheb_pts(.false., size(test_grid_y))
  max_abs_err = 0E0_wp
  do i = 1, size(test_grid_x)
    do j = 1, size(test_grid_y)
      test_x = test_grid_x(i)
      test_y = test_grid_y(j)
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      bachelier_eval = bachelier_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'rank-n eval: ', rankn_eval, &
          ' kernel eval: ', bachelier_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - bachelier_eval) )
      if (abs(rankn_eval - bachelier_eval) > (borsuk_lb + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            borsuk_lb + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - bachelier_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Bachelier kernel pass? ', test_pass
end function test_num_opt_rankn_kernel_bachelier

function test_num_opt_rankn_kernel_black(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter               :: num_toler = 3E-14_wp
  type(kernel_blackscholes)         :: black_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: err_stat
  character(len=err_msg_len)        :: err_msg
  integer                           :: i, j
  real(wp)                          :: borsuk_lb
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, black_eval, max_abs_err
  real(wp)                          :: test_grid_x(11), test_grid_y(11)
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test num-opt rank-n kernel construction for the ' &
      // 'Black-Scholes kernel: '
  ! Set the parameters of the Bachelier kernel
  call black_kernel%set_is_call(.true.)
  call black_kernel%set_sigma(1E-1_wp) ! 10% vol
  call black_kernel%set_x_lb(1E-2_wp) ! min log strike
  call black_kernel%set_x_ub(3E-1_wp) ! max log strike
  call black_kernel%set_y_lb(-1E-1_wp) ! min log forward
  call black_kernel%set_y_ub(5E-1_wp) ! max log forward
  ! Call the num_opt_rankn_kernel_asymm subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! black_kernel.
  call num_opt_rankn_kernel_asymm( &
      black_kernel, & ! The kernel to approximate
      3, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      err_stat, & ! The error status flag
      err_msg) ! The error message (if any)
  ! Output the Borsuk lower bound for this approximation
  write(unit_to_use , *) 'The Borsuk lower bound on the supremum of ' &
      // 'the approximation error is: ', borsuk_lb
  ! Check some evaluation points
  test_grid_x = black_kernel%cheb_pts(.true., size(test_grid_x))
  test_grid_y = black_kernel%cheb_pts(.false., size(test_grid_y))
  max_abs_err = 0E0_wp
  do i = 1, size(test_grid_x)
    do j = 1, size(test_grid_y)
      test_x = test_grid_x(i)
      test_y = test_grid_y(j)
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      black_eval = black_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'rank-n eval: ', rankn_eval, &
          ' kernel eval: ', black_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - black_eval) )
      if (abs(rankn_eval - black_eval) > (borsuk_lb + num_toler)) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            borsuk_lb + num_toler
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - black_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the num-opt rank-n computation for the ' &
      // 'Black-Scholes kernel pass? ', test_pass
end function test_num_opt_rankn_kernel_black

end module test_num_opt_rankn_kernel_mod