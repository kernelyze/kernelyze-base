! test_integral_rankn_disc.f90
!
! Copyright (c) 2015, 2017 by Kernelyze LLC
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
! created on: 2015-04-10
! updated on: 2015-07-07
! updated on: 2017-04-23 (replace symmetric-only num_opt_rankn_kernel
!             with num_opt_rankn_kernel_asymm)
!
! A module of unit tests for the integral_rankn_disc module.
  
module test_integral_rankn_disc_mod

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len
use kernel_mod, only : kernel
use kernel_nsect_mod, only : kernel_rankn
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use num_opt_rankn_kernel_asymm_mod, only : num_opt_rankn_kernel_asymm
use integral_nsect_disc_mod, only : integral_rankn_disc

implicit none

contains

function test_integral_rankn_disc(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality to build ' &
      // 'integrals of rank-n kernels: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_integral_rankn_disc_expprod(unit_to_use)
  test_pass = test_pass .and. test_integral_rankn_disc_gaussian(unit_to_use)
  test_pass = test_pass .and. test_integral_rankn_disc_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_integral_rankn_disc

function test_integral_rankn_disc_expprod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter                   :: num_toler = 2E-13_wp
  type(kernel_expprod)                  :: expprod_kernel
  class(kernel_rankn), allocatable      :: rankn_kernel
  class(integral_rankn_disc), &
      allocatable                       :: rankn_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: borsuk_lb, test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  real(wp), dimension( : , : ), &
      allocatable                       :: v_matrix
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_rankn_disc construction for the ' &
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
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the rank-$n$ kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_rankn_disc object
  call rankn_kernel%kernel_integral( &
      .true., & ! Is the integral over $x$?  (If not, it is over $y$.)
      weights, &  ! The weights in the discrete signed measure
      eval_pts, & ! The support points of the discrete signed measure
      rankn_integral, & ! The integral object to build
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the coefficient vector
  write(unit_to_use , *) 'The coefficient vector is: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_coeff_vec()
  ! Output the V matrix
  v_matrix = rankn_integral%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n integral is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_weights()
  ! Now try evaluating the integral of the rank-$n$ kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = rankn_integral%eval(test_pts(i))
    reference_res = dot_product( &
        expprod_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' rank-n: ', test_res , ' ref: ' , &
        reference_res
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < (borsuk_lb + num_toler) )
  end do
  write(unit_to_use , *) 'Did the integral of rank-n computation for the ' &
      // 'approx. of the exponential product kernel pass? ', test_pass
end function test_integral_rankn_disc_expprod

function test_integral_rankn_disc_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter                   :: num_toler = 2E-13_wp
  type(kernel_gaussian)                 :: gaussian_kernel
  class(kernel_rankn), allocatable      :: rankn_kernel
  class(integral_rankn_disc), &
      allocatable                       :: rankn_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: borsuk_lb, test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  real(wp), dimension( : , : ), &
      allocatable                       :: v_matrix
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_rankn_disc construction for the ' &
      // 'Gaussian kernel: '
  ! Set the "c" parameter of the Gaussian product kernel
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
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the rank-$n$ kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_rankn_disc object
  call rankn_kernel%kernel_integral( &
      .true., & ! Is the integral over $x$?  (If not, it is over $y$.)
      weights, &  ! The weights in the discrete signed measure
      eval_pts, & ! The support points of the discrete signed measure
      rankn_integral, & ! The integral object to build
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the coefficient vector
  write(unit_to_use , *) 'The coefficient vector is: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_coeff_vec()
  ! Output the V matrix
  v_matrix = rankn_integral%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n integral is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_weights()
  ! Now try evaluating the integral of the rank-$n$ kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = rankn_integral%eval(test_pts(i))
    reference_res = dot_product( &
        gaussian_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' rank-n: ', test_res , ' ref: ' , &
        reference_res
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < (borsuk_lb + num_toler) )
  end do
  write(unit_to_use , *) 'Did the integral of rank-n computation for the ' &
      // 'approx. of the Gaussian kernel pass? ', test_pass
end function test_integral_rankn_disc_gaussian

function test_integral_rankn_disc_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  ! num_toler is the amount by which any error between the
  ! rank-n approximation and the original kernel is allowed
  ! to exceed (in absolute value) the computed Borsuk lower
  ! bound without failing the test.
  real(wp), parameter                   :: num_toler = 2E-13_wp
  type(kernel_cauchy)                   :: cauchy_kernel
  class(kernel_rankn), allocatable      :: rankn_kernel
  class(integral_rankn_disc), &
      allocatable                       :: rankn_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: borsuk_lb, test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  real(wp), dimension( : , : ), &
      allocatable                       :: v_matrix
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_rankn_disc construction for the ' &
      // 'Cauchy kernel: '
  ! Set the "c" parameter of the Cauchy kernel
  call cauchy_kernel%set_c(3E0_wp)
  ! Set the "alpha" parameter of the Cauchy kernel
  call cauchy_kernel%set_alpha(1E0_wp)
  ! Call the num_opt_rankn_kernel subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! Cauchy_kernel.
  call num_opt_rankn_kernel_asymm( &
      cauchy_kernel, & ! The kernel to approximate
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
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the rank-$n$ kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_rankn_disc object
  call rankn_kernel%kernel_integral( &
      .true., & ! Is the integral over $x$?  (If not, it is over $y$.)
      weights, &  ! The weights in the discrete signed measure
      eval_pts, & ! The support points of the discrete signed measure
      rankn_integral, & ! The integral object to build
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Output the coefficient vector
  write(unit_to_use , *) 'The coefficient vector is: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_coeff_vec()
  ! Output the V matrix
  v_matrix = rankn_integral%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n integral is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') rankn_integral%get_weights()
  ! Now try evaluating the integral of the rank-$n$ kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = rankn_integral%eval(test_pts(i))
    reference_res = dot_product( &
        cauchy_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' rank-n: ', test_res , ' ref: ' , &
        reference_res
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < (borsuk_lb + num_toler) )
  end do
  write(unit_to_use , *) 'Did the integral of rank-n computation for the ' &
      // 'approx. of the Cauchy kernel pass? ', test_pass
end function test_integral_rankn_disc_cauchy

end module test_integral_rankn_disc_mod