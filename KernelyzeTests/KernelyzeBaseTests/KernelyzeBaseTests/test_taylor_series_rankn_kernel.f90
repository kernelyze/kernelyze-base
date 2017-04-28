! test_taylor_series_rankn_kernel.f90
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
! created on: 2016-02-13
! updated on: 2016-02-13
! updated on: 2016-02-14 (added output of V matrices and Taylor-series
!             matrices)
! updated on: 2016-09-30 (added Bachelier and Black-Scholes kernels)
! updated on: 2016-10-02 (necessary modifications for Bachelier kernel test)
!
! A module of unit tests for the taylor_series_rankn_kernel module.
  
module test_taylor_series_rankn_kernel_mod

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len
use kernel_mod, only : kernel
use kernel_nmonom_mod, only : kernel_rankn
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use kernel_bachelier_mod, only : kernel_bachelier
use kernel_blackscholes_mod, only : kernel_blackscholes
use taylor_series_rankn_kernel_mod , only : taylor_series_rankn_kernel

implicit none

contains

function test_taylor_series_rankn_kernel(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass

  write(unit_to_use , *) 'Test of functionality to build ' &
      // 'rank-n Taylor-series kernels: '
  
  test_pass = .true.
  
  test_pass = test_pass .and. &
      test_taylor_series_rankn_kernel_expprod(unit_to_use)
  test_pass = test_pass .and. &
      test_taylor_series_rankn_kernel_gaussian(unit_to_use)
  test_pass = test_pass .and. &
      test_taylor_series_rankn_kernel_cauchy(unit_to_use)
  test_pass = test_pass .and. &
      test_taylor_series_rankn_kernel_bachelier(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line
  
end function test_taylor_series_rankn_kernel

function test_taylor_series_rankn_kernel_expprod(unit_to_use) &
    result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  real(wp)                          :: expected_error
  type(kernel_expprod)              :: expprod_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, expprod_eval, max_abs_err
  real(wp), allocatable             :: v_matrix( : , : )
  real(wp), allocatable             :: w_matrix( : , : )
  real(wp), allocatable             :: taylor_matrix( : , : )
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test 4-term Taylor series for the ' &
      // 'exponential product kernel: '
  ! I know that the maximal error of the Taylor series occurs at
  ! (1, 1) and (-1, -1); both result in the maximal error given
  ! below for the 4-term (1, x, x**2, x**3) Taylor series.
  expected_error = ( exp(1E0_wp) - &
      (1E0_wp + 1E0_wp + 0.5E0_wp + (1E0_wp / 6E0_wp)) ) + &
      1E-14_wp ! The last term, just to handle rounding error
  ! Set the "c" parameter of the exponential product kernel
  call expprod_kernel%set_c(1E0_wp)
  ! Call the taylor_series_rankn_kernel subroutine to compute the
  ! Taylor series for the expprod_kernel.
  call taylor_series_rankn_kernel( &
      expprod_kernel, & ! The kernel to approximate
      0E0_wp, & ! The x center for the Taylor series
      0E0_wp, & ! The y center for the Taylor series
      4, & ! The approximation rank
      rankn_kernel, & ! The rank-n approximation kernel
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Get the V matrix of the rank-$n$ approximation and output it:
  v_matrix = rankn_kernel%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n kernel is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Get the W matrix of the rank-$n$ approximation and output it:
  w_matrix = rankn_kernel%get_w_matrix()  
  write(unit_to_use , *) 'The W matrix used in the rank-n kernel is: '
  do i = 1, size(w_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') w_matrix(i, :)
  end do
  ! Get and output the Taylor-series matrix
  taylor_matrix = matmul(v_matrix , transpose(w_matrix) )
  write(unit_to_use , *) 'The Taylor-series matrix is: '
  do i = 1, size(taylor_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') taylor_matrix(i, :)
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
      write(unit_to_use , *) 'Taylor-series eval: ', rankn_eval, &
          ' kernel eval: ', expprod_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - expprod_eval) )
      if (abs(rankn_eval - expprod_eval) > expected_error) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            expected_error
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - expprod_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the Taylor-series computation for the ' &
      // 'exponential product kernel pass? ', test_pass
end function test_taylor_series_rankn_kernel_expprod

function test_taylor_series_rankn_kernel_gaussian(unit_to_use) &
    result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  real(wp)                          :: expected_error
  type(kernel_gaussian)             :: gaussian_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, gaussian_eval, max_abs_err
  real(wp), allocatable             :: v_matrix( : , : )
  real(wp), allocatable             :: w_matrix( : , : )
  real(wp), allocatable             :: taylor_matrix( : , : )
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test 18-term Taylor series for the ' &
      // 'Gaussian kernel: '
  ! I know that the maximal error of the Taylor series occurs at
  ! (1, -1) and (-1, 1), but the precise form is complex due to
  ! the cross-partials of the Gaussian kernel
  expected_error = 1E-5_wp ! A reasonable guess
  ! Set the "c" parameter of the Gaussian kernel
  call gaussian_kernel%set_c(1E0_wp)
  ! Call the taylor_series_rankn_kernel subroutine to compute the
  ! Taylor series for the gaussian_kernel.
  call taylor_series_rankn_kernel( &
      gaussian_kernel, & ! The kernel to approximate
      0E0_wp, & ! The x center for the Taylor series
      0E0_wp, & ! The y center for the Taylor series
      18, & ! The approximation rank
      rankn_kernel, & ! The rank-n approximation kernel
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Get the V matrix of the rank-$n$ approximation and output it:
  v_matrix = rankn_kernel%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n kernel is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Get the W matrix of the rank-$n$ approximation and output it:
  w_matrix = rankn_kernel%get_w_matrix()  
  write(unit_to_use , *) 'The W matrix used in the rank-n kernel is: '
  do i = 1, size(w_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') w_matrix(i, :)
  end do
  ! Get and output the Taylor-series matrix
  taylor_matrix = matmul(v_matrix , transpose(w_matrix) )
  write(unit_to_use , *) 'The Taylor-series matrix is: '
  do i = 1, size(taylor_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') taylor_matrix(i, :)
  end do
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 10E0_wp
      test_y = real(j, wp) / 10E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      gaussian_eval = gaussian_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'Taylor-series eval: ', rankn_eval, &
          ' kernel eval: ', gaussian_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - gaussian_eval) )
      if (abs(rankn_eval - gaussian_eval) > expected_error) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            expected_error
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - gaussian_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the Taylor-series computation for the ' &
      // 'Gaussian kernel pass? ', test_pass
end function test_taylor_series_rankn_kernel_gaussian

function test_taylor_series_rankn_kernel_cauchy(unit_to_use) &
    result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  real(wp)                          :: expected_error
  type(kernel_cauchy)               :: cauchy_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, cauchy_eval, max_abs_err
  real(wp), allocatable             :: v_matrix( : , : )
  real(wp), allocatable             :: w_matrix( : , : )
  real(wp), allocatable             :: taylor_matrix( : , : )
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test 12-term Taylor series for the ' &
      // 'Cauchy kernel: '
  ! A reasonable guess
  expected_error = 1E-3_wp
  ! Set the "c" parameter of the Cauchy kernel
  call cauchy_kernel%set_c(3E0_wp)
  ! Set the "alpha" parameter of the Cauchy kernel
  call cauchy_kernel%set_alpha(1E0_wp)
  ! Call the taylor_series_rankn_kernel subroutine to compute the
  ! Taylor series for the cauchy_kernel.
  call taylor_series_rankn_kernel( &
      cauchy_kernel, & ! The kernel to approximate
      0E0_wp, & ! The x center for the Taylor series
      0E0_wp, & ! The y center for the Taylor series
      12, & ! The approximation rank
      rankn_kernel, & ! The rank-n approximation kernel
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Get the V matrix of the rank-$n$ approximation and output it:
  v_matrix = rankn_kernel%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n kernel is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Get the W matrix of the rank-$n$ approximation and output it:
  w_matrix = rankn_kernel%get_w_matrix()  
  write(unit_to_use , *) 'The W matrix used in the rank-n kernel is: '
  do i = 1, size(w_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') w_matrix(i, :)
  end do
  ! Get and output the Taylor-series matrix
  taylor_matrix = matmul(v_matrix , transpose(w_matrix) )
  write(unit_to_use , *) 'The Taylor-series matrix is: '
  do i = 1, size(taylor_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') taylor_matrix(i, :)
  end do
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 10E0_wp
      test_y = real(j, wp) / 10E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      cauchy_eval = cauchy_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'Taylor-series eval: ', rankn_eval, &
          ' kernel eval: ', cauchy_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - cauchy_eval) )
      if (abs(rankn_eval - cauchy_eval) > expected_error) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            expected_error
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - cauchy_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the Taylor-series computation for the ' &
      // 'Cauchy kernel pass? ', test_pass
end function test_taylor_series_rankn_kernel_cauchy
    
function test_taylor_series_rankn_kernel_bachelier(unit_to_use) &
    result(test_pass)
  ! Arguments
  integer, intent(in)               :: unit_to_use
  ! Return value
  logical                           :: test_pass
  ! Local variables
  real(wp)                          :: expected_error
  type(kernel_bachelier)            :: bachelier_kernel
  class(kernel_rankn), allocatable  :: rankn_kernel
  integer                           :: alloc_stat
  character(len=alloc_errmsg_len)   :: alloc_err_msg
  integer                           :: i, j
  real(wp)                          :: test_x, test_y
  real(wp)                          :: rankn_eval, bachelier_eval, max_abs_err
  real(wp), allocatable             :: v_matrix( : , : )
  real(wp), allocatable             :: w_matrix( : , : )
  real(wp), allocatable             :: taylor_matrix( : , : )
  ! Body
  test_pass = .true.
  write(unit_to_use , *) 'Test 3-term Taylor series for the ' &
      // 'Bachelier kernel: '
  ! A reasonable guess
  expected_error = 1E-3_wp
  ! Set the parameters of the Bachelier kernel
  call bachelier_kernel%set_is_call(.true.)
  call bachelier_kernel%set_sigma(1E-2_wp)
  ! Call the taylor_series_rankn_kernel subroutine to compute the
  ! Taylor series for the bachelier_kernel.
  call taylor_series_rankn_kernel( &
      bachelier_kernel, & ! The kernel to approximate
      0E0_wp, & ! The x center for the Taylor series
      0E0_wp, & ! The y center for the Taylor series
      3, & ! The approximation rank
      rankn_kernel, & ! The rank-n approximation kernel
      alloc_stat, & ! The allocation status flag
      alloc_err_msg) ! The allocation error message (if any)
  ! Get the V matrix of the rank-$n$ approximation and output it:
  v_matrix = rankn_kernel%get_v_matrix()
  write(unit_to_use , *) 'The V matrix used in the rank-n kernel is: '
  do i = 1, size(v_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') v_matrix(i, :)
  end do
  ! Get the W matrix of the rank-$n$ approximation and output it:
  w_matrix = rankn_kernel%get_w_matrix()
  write(unit_to_use , *) 'The W matrix used in the rank-n kernel is: '
  do i = 1, size(w_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') w_matrix(i, :)
  end do
  ! Get and output the Taylor-series matrix
  taylor_matrix = matmul(v_matrix , transpose(w_matrix) )
  write(unit_to_use , *) 'The Taylor-series matrix is: '
  do i = 1, size(taylor_matrix, 1)
    ! Use an unlimited format to control output
    write(unit_to_use , '(*(2x,e15.8))') taylor_matrix(i, :)
  end do
  ! Check some evaluation points
  max_abs_err = 0E0_wp
  do i = -10, 10, 2
    do j = -10, 10, 2
      test_x = real(i, wp) / 1000E0_wp
      test_y = real(j, wp) / 1000E0_wp
      rankn_eval = rankn_kernel%eval(test_x, test_y)
      bachelier_eval = bachelier_kernel%eval(test_x, test_y)
      write(unit_to_use , *) 'For i of ', i, ' and j of ', j, ':'
      write(unit_to_use , *) 'Taylor-series eval: ', rankn_eval, &
          ' kernel eval: ', bachelier_eval
      max_abs_err = max( max_abs_err , abs(rankn_eval - bachelier_eval) )
      if (abs(rankn_eval - bachelier_eval) > expected_error) then
        test_pass = .false.
        write(unit_to_use , *) 'Test failed at i of ', i, ' and j of ', j
        write(unit_to_use , *) 'Permitted abs val of error: ', &
            expected_error
        write(unit_to_use , *) 'Actual abs val of error: ', &
            abs(rankn_eval - bachelier_eval)
      end if
    end do
  end do
  ! Output the largest absolute value of approximation error that was observed:
  write( unit_to_use , * ) 'Largest abs val err on grid: ' , max_abs_err
  write(unit_to_use , *) 'Did the Taylor-series computation for the ' &
      // 'Bachelier kernel pass? ', test_pass
end function test_taylor_series_rankn_kernel_bachelier 

end module test_taylor_series_rankn_kernel_mod
