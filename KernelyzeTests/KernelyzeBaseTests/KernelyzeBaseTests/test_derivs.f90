! test_derivs.f90
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
! created on: 2016-01-10
! updated on: 2016-01-10
! updated on: 2016-01-16 (added rank-n kernel deriv tests)
! updated on: 2016-01-17 (finished rank-n kernel deriv tests)
! updated on: 2016-02-11 (added higher-order derivatives)
! updated on: 2016-02-12 (added further higher-order derivatives)
! updated on: 2016-02-22 (loosen tolerance for 3rd-order derivatives
!             to accommodate alternate machines)
! updated on: 2016-11-22 (further loosen 3rd-deriv tolerance to
!             accommodate NUC)
!
! A module of unit tests to check the analytical
! derivatives provided by some kernels and functions.
  
module test_derivs_mod

use set_precision, only : wp
use constants_mod, only : err_msg_len, description_len
use kernel_expprod_mod, only : kernel_expprod
use kernel_gaussian_mod, only : kernel_gaussian
use kernel_cauchy_mod, only : kernel_cauchy
use kernel_nsect_mod, only : kernel_rankn
use compfunc_kernel_mod, only : compfunc_kernel
use integral_kernel_disc_mod, only : integral_kernel_disc
use num_opt_rankn_kernel_asymm_mod, only : num_opt_rankn_kernel_asymm

implicit none

contains
  
function test_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)                       :: unit_to_use
  ! Return value
  logical                                   :: test_pass
  
  write(unit_to_use , *) 'Test of analytical kernel and function derivatives'
  
  test_pass = .true.
  
  test_pass = test_pass .and. test_expprod_derivs(unit_to_use)
  test_pass = test_pass .and. test_gaussian_derivs(unit_to_use)
  test_pass = test_pass .and. test_cauchy_derivs(unit_to_use)
  test_pass = test_pass .and. test_rankn_derivs(unit_to_use)
  test_pass = test_pass .and. test_compfunc_derivs(unit_to_use)
  test_pass = test_pass .and. test_integral_derivs(unit_to_use)
  
  write(unit_to_use , *) ' ' ! Output a blank line

end function test_derivs

function test_expprod_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(kernel_expprod)        :: expprod_kernel
  real(wp)                    :: x
  real(wp)                    :: y
  integer                     :: i
  integer                     :: j
  real(wp)                    :: dkernel_dx
  real(wp)                    :: dkernel_dy
  real(wp)                    :: d2kernel_dx_dx
  real(wp)                    :: d2kernel_dy_dy
  real(wp)                    :: d2kernel_dx_dy
  real(wp)                    :: d3kernel_dx_dx_dy
  real(wp)                    :: d3kernel_dx_dy_dy
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dy
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dy_dy
  real(wp)                    :: findiff_dx_dy
  real(wp)                    :: findiff_dx_dx_dy
  real(wp)                    :: findiff_dx_dy_dy
  ! Body
  test_pass = .true.
  ! Set up the kernel
  call expprod_kernel%set_c(1E0_wp)
  write(unit_to_use , *) 'Check the derivs of an exponential product kernel: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    do j = 1, 21
      y = real(j - 11, wp) / 10E0_wp
      ! Analytical derivatives:
      dkernel_dx        = expprod_kernel%first_partial_dx(x , y)
      dkernel_dy        = expprod_kernel%first_partial_dy(x , y)
      d2kernel_dx_dx    = expprod_kernel%second_partial_dx_dx(x , y)
      d2kernel_dy_dy    = expprod_kernel%second_partial_dy_dy(x , y)
      d2kernel_dx_dy    = expprod_kernel%second_partial_dx_dy(x , y)
      d3kernel_dx_dx_dy = expprod_kernel%mth_nth_partial(2, 1, x, y)
      d3kernel_dx_dy_dy = expprod_kernel%mth_nth_partial(1, 2, x, y)
      ! Finite differences:
      findiff_dx        = expprod_kernel%finite_diff_dx(x , y)
      findiff_dy        = expprod_kernel%finite_diff_dy(x , y)
      findiff_dx_dx     = expprod_kernel%finite_diff_dx_dx(x , y)
      findiff_dy_dy     = expprod_kernel%finite_diff_dy_dy(x , y)
      findiff_dx_dy     = expprod_kernel%finite_diff_dx_dy(x , y)
      findiff_dx_dx_dy  = expprod_kernel%mth_nth_finite_diff(2, 1, x, y)
      findiff_dx_dy_dy  = expprod_kernel%mth_nth_finite_diff(1, 2, x, y)
      write(unit_to_use , *) 'Arguments x = ', x, ' and y = ', y
      write(unit_to_use , *) 'Partial w. r. t. x, analytical ', &
          'and finite diff: ', dkernel_dx, ' ', findiff_dx
      write(unit_to_use , *) 'Partial w. r. t. y, analytical ', &
          'and finite diff: ', dkernel_dy, ' ', findiff_dy
      write(unit_to_use , *) '2nd partial w. r. t. x, analytical ', &
          'and finite diff: ', d2kernel_dx_dx, ' ', findiff_dx_dx
      write(unit_to_use , *) '2nd partial w. r. t. y, analytical ', &
          'and finite diff: ', d2kernel_dy_dy, ' ', findiff_dy_dy
      write(unit_to_use , *) 'Cross partial, analytical ', &
          'and finite diff: ', d2kernel_dx_dy, ' ', findiff_dx_dy
      write(unit_to_use , *) '3rd partial, dx dx dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dx_dy, ' ', findiff_dx_dx_dy
      write(unit_to_use , *) '3rd partial, dx dy dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dy_dy, ' ', findiff_dx_dy_dy
      if ( (abs(dkernel_dx        - findiff_dx)       > toler) .or. &
           (abs(dkernel_dy        - findiff_dy)       > toler) .or. &
           (abs(d2kernel_dx_dx    - findiff_dx_dx)    > toler) .or. &
           (abs(d2kernel_dy_dy    - findiff_dy_dy)    > toler) .or. &
           (abs(d2kernel_dx_dy    - findiff_dx_dy)    > toler) .or. &
           (abs(d3kernel_dx_dx_dy - findiff_dx_dx_dy) > 1.5E1_wp * toler) .or. &
           (abs(d3kernel_dx_dy_dy - findiff_dx_dy_dy) > 1.5E1_wp * toler) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Expprod test failed: analytical and ', &
            'finite-difference derivatives differ too much for x = ', &
            x, ' and y = ', y
      end if
    end do
  end do
end function test_expprod_derivs

function test_gaussian_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(kernel_gaussian)       :: gaussian_kernel
  real(wp)                    :: x
  real(wp)                    :: y
  integer                     :: i
  integer                     :: j
  real(wp)                    :: dkernel_dx
  real(wp)                    :: dkernel_dy
  real(wp)                    :: d2kernel_dx_dx
  real(wp)                    :: d2kernel_dy_dy
  real(wp)                    :: d2kernel_dx_dy
  real(wp)                    :: d3kernel_dx_dx_dy
  real(wp)                    :: d3kernel_dx_dy_dy
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dy
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dy_dy
  real(wp)                    :: findiff_dx_dy
  real(wp)                    :: findiff_dx_dx_dy
  real(wp)                    :: findiff_dx_dy_dy
  ! Body
  test_pass = .true.
  ! Set up the kernel
  call gaussian_kernel%set_c(1E0_wp)
  write(unit_to_use , *) 'Check the derivs of a Gaussian kernel: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    do j = 1, 21
      y = real(j - 11, wp) / 10E0_wp
      ! Analytical derivatives:
      dkernel_dx        = gaussian_kernel%first_partial_dx(x , y)
      dkernel_dy        = gaussian_kernel%first_partial_dy(x , y)
      d2kernel_dx_dx    = gaussian_kernel%second_partial_dx_dx(x , y)
      d2kernel_dy_dy    = gaussian_kernel%second_partial_dy_dy(x , y)
      d2kernel_dx_dy    = gaussian_kernel%second_partial_dx_dy(x , y)
      d3kernel_dx_dx_dy = gaussian_kernel%mth_nth_partial(2, 1, x, y)
      d3kernel_dx_dy_dy = gaussian_kernel%mth_nth_partial(1, 2, x, y)
      ! Finite differences:
      findiff_dx        = gaussian_kernel%finite_diff_dx(x , y)
      findiff_dy        = gaussian_kernel%finite_diff_dy(x , y)
      findiff_dx_dx     = gaussian_kernel%finite_diff_dx_dx(x , y)
      findiff_dy_dy     = gaussian_kernel%finite_diff_dy_dy(x , y)
      findiff_dx_dy     = gaussian_kernel%finite_diff_dx_dy(x , y)
      findiff_dx_dx_dy  = gaussian_kernel%mth_nth_finite_diff(2, 1, x, y)
      findiff_dx_dy_dy  = gaussian_kernel%mth_nth_finite_diff(1, 2, x, y)
      write(unit_to_use , *) 'Arguments x = ', x, ' and y = ', y
      write(unit_to_use , *) 'Partial w. r. t. x, analytical ', &
          'and finite diff: ', dkernel_dx, ' ', findiff_dx
      write(unit_to_use , *) 'Partial w. r. t. y, analytical ', &
          'and finite diff: ', dkernel_dy, ' ', findiff_dy
      write(unit_to_use , *) '2nd partial w. r. t. x, analytical ', &
          'and finite diff: ', d2kernel_dx_dx, ' ', findiff_dx_dx
      write(unit_to_use , *) '2nd partial w. r. t. y, analytical ', &
          'and finite diff: ', d2kernel_dy_dy, ' ', findiff_dy_dy
      write(unit_to_use , *) 'Cross partial, analytical ', &
          'and finite diff: ', d2kernel_dx_dy, ' ', findiff_dx_dy
      write(unit_to_use , *) '3rd partial, dx dx dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dx_dy, ' ', findiff_dx_dx_dy
      write(unit_to_use , *) '3rd partial, dx dy dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dy_dy, ' ', findiff_dx_dy_dy
      if ( (abs(dkernel_dx        - findiff_dx)       > toler) .or. &
           (abs(dkernel_dy        - findiff_dy)       > toler) .or. &
           (abs(d2kernel_dx_dx    - findiff_dx_dx)    > toler) .or. &
           (abs(d2kernel_dy_dy    - findiff_dy_dy)    > toler) .or. &
           (abs(d2kernel_dx_dy    - findiff_dx_dy)    > toler) .or. &
           (abs(d3kernel_dx_dx_dy - findiff_dx_dx_dy) > 1.5E1_wp * toler) .or. &
           (abs(d3kernel_dx_dy_dy - findiff_dx_dy_dy) > 1.5E1_wp * toler) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Gaussian test failed: analytical and ', &
            'finite-difference derivatives differ too much for x = ', &
            x, ' and y = ', y
      end if
    end do
  end do
end function test_gaussian_derivs

function test_cauchy_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(kernel_cauchy)         :: cauchy_kernel
  real(wp)                    :: x
  real(wp)                    :: y
  integer                     :: i
  integer                     :: j
  real(wp)                    :: dkernel_dx
  real(wp)                    :: dkernel_dy
  real(wp)                    :: d2kernel_dx_dx
  real(wp)                    :: d2kernel_dy_dy
  real(wp)                    :: d2kernel_dx_dy
  real(wp)                    :: d3kernel_dx_dx_dy
  real(wp)                    :: d3kernel_dx_dy_dy
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dy
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dy_dy
  real(wp)                    :: findiff_dx_dy
  real(wp)                    :: findiff_dx_dx_dy
  real(wp)                    :: findiff_dx_dy_dy
  ! Body
  test_pass = .true.
  ! Set up the kernel
  call cauchy_kernel%set_c(3E0_wp)
  call cauchy_kernel%set_alpha(1E0_wp)
  write(unit_to_use , *) 'Check the derivs of a Cauchy kernel: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    do j = 1, 21
      y = real(j - 11, wp) / 10E0_wp
      ! Analytical derivatives:
      dkernel_dx        = cauchy_kernel%first_partial_dx(x , y)
      dkernel_dy        = cauchy_kernel%first_partial_dy(x , y)
      d2kernel_dx_dx    = cauchy_kernel%second_partial_dx_dx(x , y)
      d2kernel_dy_dy    = cauchy_kernel%second_partial_dy_dy(x , y)
      d2kernel_dx_dy    = cauchy_kernel%second_partial_dx_dy(x , y)
      d3kernel_dx_dx_dy = cauchy_kernel%mth_nth_partial(2, 1, x, y)
      d3kernel_dx_dy_dy = cauchy_kernel%mth_nth_partial(1, 2, x, y)
      ! Finite differences:
      findiff_dx        = cauchy_kernel%finite_diff_dx(x , y)
      findiff_dy        = cauchy_kernel%finite_diff_dy(x , y)
      findiff_dx_dx     = cauchy_kernel%finite_diff_dx_dx(x , y)
      findiff_dy_dy     = cauchy_kernel%finite_diff_dy_dy(x , y)
      findiff_dx_dy     = cauchy_kernel%finite_diff_dx_dy(x , y)
      findiff_dx_dx_dy  = cauchy_kernel%mth_nth_finite_diff(2, 1, x, y)
      findiff_dx_dy_dy  = cauchy_kernel%mth_nth_finite_diff(1, 2, x, y)
      write(unit_to_use , *) 'Arguments x = ', x, ' and y = ', y
      write(unit_to_use , *) 'Partial w. r. t. x, analytical ', &
          'and finite diff: ', dkernel_dx, ' ', findiff_dx
      write(unit_to_use , *) 'Partial w. r. t. y, analytical ', &
          'and finite diff: ', dkernel_dy, ' ', findiff_dy
      write(unit_to_use , *) '2nd partial w. r. t. x, analytical ', &
          'and finite diff: ', d2kernel_dx_dx, ' ', findiff_dx_dx
      write(unit_to_use , *) '2nd partial w. r. t. y, analytical ', &
          'and finite diff: ', d2kernel_dy_dy, ' ', findiff_dy_dy
      write(unit_to_use , *) 'Cross partial, analytical ', &
          'and finite diff: ', d2kernel_dx_dy, ' ', findiff_dx_dy
      write(unit_to_use , *) '3rd partial, dx dx dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dx_dy, ' ', findiff_dx_dx_dy
      write(unit_to_use , *) '3rd partial, dx dy dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dy_dy, ' ', findiff_dx_dy_dy
      if ( (abs(dkernel_dx        - findiff_dx)       > toler) .or. &
           (abs(dkernel_dy        - findiff_dy)       > toler) .or. &
           (abs(d2kernel_dx_dx    - findiff_dx_dx)    > toler) .or. &
           (abs(d2kernel_dy_dy    - findiff_dy_dy)    > toler) .or. &
           (abs(d2kernel_dx_dy    - findiff_dx_dy)    > toler) .or. &
           (abs(d3kernel_dx_dx_dy - findiff_dx_dx_dy) > 1.5E1_wp * toler) .or. &
           (abs(d3kernel_dx_dy_dy - findiff_dx_dy_dy) > 1.5E1_wp * toler) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Cauchy test failed: analytical and ', &
            'finite-difference derivatives differ too much for x = ', &
            x, ' and y = ', y
      end if
    end do
  end do
end function test_cauchy_derivs

function test_rankn_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(kernel_expprod)        :: expprod_kernel
  class(kernel_rankn), &
                 allocatable  :: rankn_kernel
  real(wp)                    :: borsuk_lb
  integer                     :: local_err_stat
  character(len=err_msg_len)  :: local_err_msg
  real(wp)                    :: x
  real(wp)                    :: y
  integer                     :: i
  integer                     :: j
  real(wp)                    :: dkernel_dx
  real(wp)                    :: dkernel_dy
  real(wp)                    :: d2kernel_dx_dx
  real(wp)                    :: d2kernel_dy_dy
  real(wp)                    :: d2kernel_dx_dy
  real(wp)                    :: d3kernel_dx_dx_dy
  real(wp)                    :: d3kernel_dx_dy_dy
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dy
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dy_dy
  real(wp)                    :: findiff_dx_dy
  real(wp)                    :: findiff_dx_dx_dy
  real(wp)                    :: findiff_dx_dy_dy
  ! Body
  test_pass = .true.
  ! Set up the kernel
  call expprod_kernel%set_c(1E0_wp)
  ! Call the num_opt_rankn_kernel subroutine to compute the
  ! numerically-optimal rank-$n$ kernel approximation to
  ! expprod_kernel.
  call num_opt_rankn_kernel_asymm( &
      expprod_kernel, & ! The kernel to approximate
      4, & ! The approximation rank
      100, &  ! The maximum number of iterations
      1e-15_wp, & ! The numerical tolerance
      rankn_kernel, & ! The rank-n approximation kernel
      borsuk_lb, &  ! The Borsuk lower bound for the sup of the error
      local_err_stat, & ! The allocation status flag
      local_err_msg) ! The allocation error message (if any)
  write(unit_to_use , *) 'Check the derivs of a rank-n kernel: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    do j = 1, 21
      y = real(j - 11, wp) / 10E0_wp
      ! Analytical derivatives:
      dkernel_dx        = rankn_kernel%first_partial_dx(x , y)
      dkernel_dy        = rankn_kernel%first_partial_dy(x , y)
      d2kernel_dx_dx    = rankn_kernel%second_partial_dx_dx(x , y)
      d2kernel_dy_dy    = rankn_kernel%second_partial_dy_dy(x , y)
      d2kernel_dx_dy    = rankn_kernel%second_partial_dx_dy(x , y)
      d3kernel_dx_dx_dy = rankn_kernel%mth_nth_partial(2, 1, x, y)
      d3kernel_dx_dy_dy = rankn_kernel%mth_nth_partial(1, 2, x, y)
      ! Finite differences:
      findiff_dx        = rankn_kernel%finite_diff_dx(x , y)
      findiff_dy        = rankn_kernel%finite_diff_dy(x , y)
      findiff_dx_dx     = rankn_kernel%finite_diff_dx_dx(x , y)
      findiff_dy_dy     = rankn_kernel%finite_diff_dy_dy(x , y)
      findiff_dx_dy     = rankn_kernel%finite_diff_dx_dy(x , y)
      findiff_dx_dx_dy  = rankn_kernel%mth_nth_finite_diff(2, 1, x, y)
      findiff_dx_dy_dy  = rankn_kernel%mth_nth_finite_diff(1, 2, x, y)
      write(unit_to_use , *) 'Arguments x = ', x, ' and y = ', y
      write(unit_to_use , *) 'Partial w. r. t. x, analytical ', &
          'and finite diff: ', dkernel_dx, ' ', findiff_dx
      write(unit_to_use , *) 'Partial w. r. t. y, analytical ', &
          'and finite diff: ', dkernel_dy, ' ', findiff_dy
      write(unit_to_use , *) '2nd partial w. r. t. x, analytical ', &
          'and finite diff: ', d2kernel_dx_dx, ' ', findiff_dx_dx
      write(unit_to_use , *) '2nd partial w. r. t. y, analytical ', &
          'and finite diff: ', d2kernel_dy_dy, ' ', findiff_dy_dy
      write(unit_to_use , *) 'Cross partial, analytical ', &
          'and finite diff: ', d2kernel_dx_dy, ' ', findiff_dx_dy
      write(unit_to_use , *) '3rd partial, dx dx dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dx_dy, ' ', findiff_dx_dx_dy
      write(unit_to_use , *) '3rd partial, dx dy dy, analytical ', &
          'and finite diff: ', d3kernel_dx_dy_dy, ' ', findiff_dx_dy_dy
      if ( (abs(dkernel_dx        - findiff_dx)       > toler) .or. &
           (abs(dkernel_dy        - findiff_dy)       > toler) .or. &
           (abs(d2kernel_dx_dx    - findiff_dx_dx)    > toler) .or. &
           (abs(d2kernel_dy_dy    - findiff_dy_dy)    > toler) .or. &
           (abs(d2kernel_dx_dy    - findiff_dx_dy)    > toler) .or. &
           (abs(d3kernel_dx_dx_dy - findiff_dx_dx_dy) > 3E1_wp * toler) .or. &
           (abs(d3kernel_dx_dy_dy - findiff_dx_dy_dy) > 3E1_wp * toler) ) then
        test_pass = .false.
        write(unit_to_use , *) 'Rank-n test failed: analytical and ', &
            'finite-difference derivatives differ too much for x = ', &
            x, ' and y = ', y
      end if
    end do
  end do
end function test_rankn_derivs

function test_compfunc_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(compfunc_kernel)       :: testfunc
  type(kernel_expprod)        :: expprod_kernel
  real(wp)                    :: x
  integer                     :: i
  real(wp)                    :: dfunc_dx
  real(wp)                    :: d2func_dx_dx
  real(wp)                    :: d3func_dx_dx_dx
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dx_dx_dx
  integer                     :: local_err_stat
  character(len=err_msg_len)  :: local_err_msg
  ! Body
  test_pass = .true.
  ! Set up the kernel
  call expprod_kernel%set_c(1E0_wp)
  call testfunc%set_kernel(expprod_kernel, local_err_stat, local_err_msg)
  if (local_err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use , *) 'Failed to set testfunc kernel, message: ' &
        // local_err_msg
    return
  end if
  call testfunc%set_eval_pt(5E-1_wp)
  write(unit_to_use , *) 'Check the derivs of a kernel section: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    dfunc_dx         = testfunc%first_deriv( x )
    d2func_dx_dx     = testfunc%second_deriv( x )
    d3func_dx_dx_dx  = testfunc%nth_deriv( 3 , x )
    findiff_dx       = testfunc%finite_diff_dx( x )
    findiff_dx_dx    = testfunc%finite_diff_dx_dx( x )
    findiff_dx_dx_dx = testfunc%nth_finite_diff( 3 , x )
    write(unit_to_use , *) 'Argument x = ', x
    write(unit_to_use , *) 'First derivative, analytical ', &
        'and finite diff: ', dfunc_dx, ' ', findiff_dx
    write(unit_to_use , *) 'Second derivative, analytical ', &
        'and finite diff: ', d2func_dx_dx, ' ', findiff_dx_dx
    write(unit_to_use , *) 'Third derivative, analytical ', &
        'and finite diff: ', d3func_dx_dx_dx, ' ', findiff_dx_dx_dx
    if ( (abs(dfunc_dx        - findiff_dx)       > toler) .or. &
         (abs(d2func_dx_dx    - findiff_dx_dx)    > toler) .or. &
         (abs(d3func_dx_dx_dx - findiff_dx_dx_dx) > 1.5E1_wp * toler) ) then
      test_pass = .false.
      write(unit_to_use , *) 'Compfunc test failed: analytical and ', &
          'finite-difference derivatives differ too much for x = ', x
    end if
  end do
end function test_compfunc_derivs

function test_integral_derivs(unit_to_use) result(test_pass)
  ! Arguments
  integer, intent(in)         :: unit_to_use
  ! Return value
  logical                     :: test_pass
  ! Parameter for error tolerance
  real(wp), parameter         :: toler = 1E-4_wp
  ! Local variables
  type(integral_kernel_disc)  :: kernel_integral
  type(kernel_gaussian)       :: gaussian_kernel
  real(wp)                    :: x
  integer                     :: i
  integer                     :: eval_sz
  real(wp)                    :: dfunc_dx
  real(wp)                    :: d2func_dx_dx
  real(wp)                    :: d3func_dx_dx_dx
  real(wp)                    :: findiff_dx
  real(wp)                    :: findiff_dx_dx
  real(wp)                    :: findiff_dx_dx_dx
  real(wp)                    :: eval_pts(11)
  real(wp)                    :: weights(11)
  integer                     :: local_err_stat
  character(len=err_msg_len)  :: local_err_msg
  character( &
        len=description_len)  :: descrip_str
  ! Body
  test_pass = .true.
  ! Set up the kernel
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
  call kernel_integral%set_kernel( &
      gaussian_kernel, local_err_stat, local_err_msg)
  if (local_err_stat /= 0) then
    test_pass = .false.
    write(unit_to_use , *) 'Failed to set kernel of integral, message: ' &
        // local_err_msg
    return
  end if
  call kernel_integral%set_eval_pts(eval_pts) ! Set the evaluation points
  call kernel_integral%set_weights(weights) ! Set the weights
  write(unit_to_use , *) 'Check the derivs of an integral of a kernel: '
  do i = 1, 21
    x = real(i - 11, wp) / 10E0_wp
    dfunc_dx         = kernel_integral%first_deriv( x )
    d2func_dx_dx     = kernel_integral%second_deriv( x )
    d3func_dx_dx_dx  = kernel_integral%nth_deriv( 3 , x )
    findiff_dx       = kernel_integral%finite_diff_dx( x )
    findiff_dx_dx    = kernel_integral%finite_diff_dx_dx( x )
    findiff_dx_dx_dx = kernel_integral%nth_finite_diff( 3 , x )
    write(unit_to_use , *) 'Argument x = ', x
    write(unit_to_use , *) 'First derivative, analytical ', &
        'and finite diff: ', dfunc_dx, ' ', findiff_dx
    write(unit_to_use , *) 'Second derivative, analytical ', &
        'and finite diff: ', d2func_dx_dx, ' ', findiff_dx_dx
    write(unit_to_use , *) 'Third derivative, analytical ', &
        'and finite diff: ', d3func_dx_dx_dx, ' ', findiff_dx_dx_dx
    if ( (abs(dfunc_dx        - findiff_dx)       > toler) .or. &
         (abs(d2func_dx_dx    - findiff_dx_dx)    > toler) .or. &
         (abs(d3func_dx_dx_dx - findiff_dx_dx_dx) > 1.5E1_wp * toler) ) then
      test_pass = .false.
      write(unit_to_use , *) 'Integral test failed: analytical and ', &
          'finite-difference derivatives differ too much for x = ', x
    end if
  end do
end function test_integral_derivs

end module test_derivs_mod
