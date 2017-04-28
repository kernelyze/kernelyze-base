! test_integral_kernel_disc.f90
!
! Copyright (c) 2015 by Kernelyze LLC
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
! created on: 2015-04-16
! updated on: 2015-07-07
!
! A module of unit tests for the integral_kernel_disc module.
  
module test_integral_kernel_disc_mod

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len, description_len
use kernel_nsect_mod, only : kernel_rankn
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use integral_kernel_disc_mod, only : integral_kernel_disc

implicit none

contains

function test_integral_kernel_disc(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality to build ' &
      // 'integrals of kernels w. r. t. discrete signed measures: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_integral_kernel_disc_expprod(unit_to_use)
  test_pass = test_pass .and. test_integral_kernel_disc_gaussian(unit_to_use)
  test_pass = test_pass .and. test_integral_kernel_disc_cauchy(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_integral_kernel_disc

function test_integral_kernel_disc_expprod(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  type(kernel_expprod)                  :: expprod_kernel
  type(integral_kernel_disc)            :: kernel_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  character(len=description_len)        :: descrip_str
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_kernel_disc construction for the ' &
      // 'exponential product kernel: '
  ! Set the "c" parameter of the exponential product kernel
  call expprod_kernel%set_c(1E0_wp)
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_kernel_disc object
  descrip_str = 'Test integral of kernel'
  call kernel_integral%set_description(descrip_str) ! Description string
  call kernel_integral%set_integral_vbl(.true.) ! The integral is over $x$
  ! Set the kernel
  call kernel_integral%set_kernel(expprod_kernel, alloc_stat, alloc_err_msg)
  call kernel_integral%set_eval_pts(eval_pts) ! Set the evaluation points
  call kernel_integral%set_weights(weights) ! Set the weights
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_weights()
  ! Now try evaluating the integral of the kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = kernel_integral%eval(test_pts(i))
    reference_res = dot_product( &
        expprod_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' kernel: ', test_res , ' ref: ' , &
        reference_res
    ! Any difference should be only rounding error
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < ( 10E0_wp * epsilon(test_res) ) )
  end do
  write(unit_to_use , *) 'Did the integral of ' &
      // 'the exponential product kernel pass? ', test_pass
end function test_integral_kernel_disc_expprod

function test_integral_kernel_disc_gaussian(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  type(kernel_gaussian)                 :: gaussian_kernel
  type(integral_kernel_disc)            :: kernel_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  character(len=description_len)        :: descrip_str
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_kernel_disc construction for the ' &
      // 'Gaussian kernel: '
  ! Set the "c" parameter of the Gaussian kernel
  call gaussian_kernel%set_c(1E0_wp)
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_kernel_disc object
  descrip_str = 'Test integral of kernel'
  call kernel_integral%set_description(descrip_str) ! Description string
  call kernel_integral%set_integral_vbl(.true.) ! The integral is over $x$
  ! Set the kernel
  call kernel_integral%set_kernel(gaussian_kernel, alloc_stat, alloc_err_msg)
  call kernel_integral%set_eval_pts(eval_pts) ! Set the evaluation points
  call kernel_integral%set_weights(weights) ! Set the weights
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_weights()
  ! Now try evaluating the integral of the kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = kernel_integral%eval(test_pts(i))
    reference_res = dot_product( &
        gaussian_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' kernel: ', test_res , ' ref: ' , &
        reference_res
    ! Any difference should be only rounding error
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < ( 10E0_wp * epsilon(test_res) ) )
  end do
  write(unit_to_use , *) 'Did the integral of ' &
      // 'the Gaussian kernel pass? ', test_pass
end function test_integral_kernel_disc_gaussian

function test_integral_kernel_disc_cauchy(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                   :: unit_to_use
  ! Return value
  logical                               :: test_pass
  ! Local variables
  type(kernel_cauchy)                   :: cauchy_kernel
  type(integral_kernel_disc)            :: kernel_integral
  integer                               :: alloc_stat
  character(len=alloc_errmsg_len)       :: alloc_err_msg
  character(len=description_len)        :: descrip_str
  integer                               :: i, eval_sz, test_sz
  real(wp)                              :: test_res, reference_res
  real(wp), dimension(21)               :: test_pts
  real(wp), dimension(11)               :: eval_pts, weights 
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test integral_kernel_disc construction for the ' &
      // 'Cauchy kernel: '
  ! Set the "c" parameter of the Cauchy kernel
  call cauchy_kernel%set_c(3E0_wp)
  ! Set the "alpha" parameter of the Cauchy kernel
  call cauchy_kernel%set_alpha(1E0_wp)
  ! Set the weights and support points of the discrete signed measure
  ! to use in integrating the kernel
  eval_sz = size(eval_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  eval_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(eval_sz - 1, wp) ) ), &
      i = 1, eval_sz ) /)
  ! These are simply equal weights
  weights = (/ ( ( 1E0_wp / real(eval_sz , wp) ) , i = 1 , eval_sz ) /)
  ! Build the integral_kernel_disc object
  descrip_str = 'Test integral of kernel'
  call kernel_integral%set_description(descrip_str) ! Description string
  call kernel_integral%set_integral_vbl(.true.) ! The integral is over $x$
  ! Set the kernel
  call kernel_integral%set_kernel(cauchy_kernel, alloc_stat, alloc_err_msg)
  call kernel_integral%set_eval_pts(eval_pts) ! Set the evaluation points
  call kernel_integral%set_weights(weights) ! Set the weights
  ! Output the eval points
  write(unit_to_use , *) 'The evaluation points of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_eval_pts()
  ! Output the weights
  write(unit_to_use , *) 'The weights of the integral are: '
  write(unit_to_use , '(*(2x,e15.8))') kernel_integral%get_weights()
  ! Now try evaluating the integral of the kernel at a number
  ! of points:
  test_sz = size(test_pts)
  ! Note that the implied do loop below produces evenly-spaced points
  ! from -1E0_wp to 1E0_wp.
  test_pts = (/ &
      ( -1E0_wp + ( real(i - 1, wp) * ( 2E0_wp / real(test_sz - 1, wp) ) ), &
      i = 1, test_sz ) /)
  do i = 1, test_sz
    test_res = kernel_integral%eval(test_pts(i))
    reference_res = dot_product( &
        cauchy_kernel%eval_elt(test_pts(i) , eval_pts) , weights)
    write(unit_to_use , *) test_pts(i) , &
        ' kernel: ', test_res , ' ref: ' , &
        reference_res
    ! Any difference should be only rounding error
    test_pass = test_pass .and. &
      ( abs(test_res - reference_res) < ( 10E0_wp * epsilon(test_res) ) )
  end do
  write(unit_to_use , *) 'Did the integral of ' &
      // 'the Cauchy kernel pass? ', test_pass
end function test_integral_kernel_disc_cauchy

end module test_integral_kernel_disc_mod