!  fort_approx.f90 
!
!  Copyright (C) 2015, 2016, 2017 by Kernelyze LLC
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
!  created on: 2015-02-01
!  updated on: 2016-01-18
!  updated on: 2016-02-06
!  updated on: 2016-02-11
!  updated on: 2016-02-19
!  updated on: 2016-02-23
!  updated on: 2016-03-02 (temporary check of effects on
!              test_build_and_eval_num_opt)
!  updated on: 2016-04-08 (added tk_threads tests)
!  updated on: 2016-04-12 (added tk_threads maps tests)
!  updated on: 2016-04-24 (removed test_fmin_twovar)
!  updated on: 2016-04-25 (initialize kernel_test_funcs module variables
!              and initialize intrinsic_test_funcs module variables)
!  updated on: 2016-05-06 (added testing of the 
!              tk_threads_multivar_kernel_integral_store and of the 
!              discrete_random_numbers functionality)
!  updated on: 2016-05-12 (added gradient and Hessian tests)
!  updated on: 2016-05-19 (test thread-safe stores with variable
!              length strings as keys)
!  updated on: 2016-07-01 (rename worst_rho_for_kernel
!             -> worst_rho_for_kernel_symm)
!  updated on: 2016-07-03 (add tests for compute_a_and_b_matrices,
!              for worst_rho_and_gamma, and for all_num_opt_mats_asymm)
!  updated on: 2016-07-05 (add tests for num_opt_rankn_params)
!  updated on: 2016-07-13 (add gaussian_pdf tests)
!  updated on: 2016-08-23 (add num_opt_rankn_trivar tests)
!  updated on: 2016-08-29 (add trivar_integral_rankn tests and
!              singlevar_proc_ptr tests)
!  updated on: 2016-09-09 (added tk_threads_trivar_store tests)
!  updated on: 2016-10-04 (added test_option_port_approx tests)
!  updated on: 2016-11-24 (added kernel_integral_option tests)
!  updated on: 2017-01-01 (added singular function series tests)
!  updated on: 2017-04-23 (subset to only Base tests)
!
!  FUNCTIONS:
!  FortApprox - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: fort_approx
!
!  PURPOSE:  Runs entire suite of tests.
!
!****************************************************************************

program fort_approx
  
  use intrinsic_test_funcs_mod
  use kernel_test_funcs_mod
  use test_set_precision_mod, only : test_set_precision
  use test_constants_mod, only : test_constants
  use test_elemental_mod, only : test_elemental
  use test_unirnk_mod, only : test_unirnk
  use test_mrgrnk_mod, only : test_mrgrnk
  use test_quicksort_mod, only : test_quicksort
  use test_binomial_coeffs_mod, only : test_binomial_coeffs
  use test_factorial_mod, only : test_factorial
  use test_falling_factorial_mod, only : test_falling_factorial
  use test_hermite_polynomial_mod, only : test_hermite_polynomial
  use test_singlevar_proc_ptr_mod, only : test_singlevar_proc_ptr
  use test_brent_mod, only : test_brent
  use test_moore_penrose_inverse_mod, only : test_moore_penrose_inverse
  use test_compute_a_matrix_mod, only : test_compute_a_matrix
  use test_compute_a_and_b_matrices_mod, only : test_compute_a_and_b_matrices
  use test_find_all_zeros_mod, only : test_find_all_zeros
  use test_find_rel_optima_mod, only : test_find_rel_optima
  use test_chebyshev_points_mod, only : test_chebyshev_points
  use test_remez_step_mod, only : test_remez_step
  use test_borsuk_lower_bound_mod, only : test_borsuk_lower_bound
  use test_worst_rho_and_gamma_mod, only : test_worst_rho_and_gamma
  use test_all_num_opt_mats_asymm_mod, only : test_all_num_opt_mats_asymm
  use test_num_opt_rankn_kernel_mod, only : test_num_opt_rankn_kernel
  use test_num_opt_rankn_params_mod, only : test_num_opt_rankn_params
  use test_taylor_series_rankn_kernel_mod, &
      only : test_taylor_series_rankn_kernel
  use test_eigen_rankn_kernel_mod, only : test_eigen_rankn_kernel
  use test_integral_kernel_disc_mod, only : test_integral_kernel_disc
  use test_integral_rankn_disc_mod, only : test_integral_rankn_disc
  use test_derivs_mod, only : test_derivs
  use, intrinsic :: iso_fortran_env, only : output_unit
  
  implicit none
  
  ! True if all tests pass
  logical :: test_pass
  ! Unit to use for file output
  integer :: test_unit
  
  test_pass = .true.
  
  ! Open an output file for write
  open(file = 'test_output.txt', form = 'formatted', newunit = test_unit)
  
  ! Initialize the intrinsic_test_funcs function objects:
  call init_test_funcs()
  ! Initialize the kernel_test_funcs kernels:
  call init_test_kernels()
    
  ! Run the unit tests for the precision that is set for "wp" 
  ! (working precision)
  test_pass = test_pass .and. test_set_precision(test_unit)
  write(test_unit , *) &
      'All set_precision tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for constants used
  test_pass = test_pass .and. test_constants(test_unit)
  write(test_unit , *) &
      'All constants tests completed; did all pass so far?  ', test_pass
  ! Run the tests for elemental function behavior
  test_pass = test_pass .and. test_elemental(test_unit)
  write(test_unit , *) &
      'All elemental tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the unique ranking functionality
  test_pass = test_pass .and. test_unirnk(test_unit)
  write(test_unit , *) &
      'All unirnk tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the merge-sort-based ranking functionality
  test_pass = test_pass .and. test_mrgrnk(test_unit)
  write(test_unit , *) &
      'All mrgrnk tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the quicksort functionality
  test_pass = test_pass .and. test_quicksort(test_unit)
  write(test_unit , *) &
      'All quicksort tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the binomial-coefficients functionality
  test_pass = test_pass .and. test_binomial_coeffs(test_unit)
  write(test_unit , *) &
      'All binomial-coeff. tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the factorial functionality
  test_pass = test_pass .and. test_factorial(test_unit)
  write(test_unit , *) &
      'All factorial tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the falling factorial functionality
  test_pass = test_pass .and. test_falling_factorial(test_unit)
  write(test_unit , *) &
      'All falling fact. tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the Hermite polynomial functionality
  test_pass = test_pass .and. test_hermite_polynomial(test_unit)
  write(test_unit , *) &
      'All Hermite polynomial tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the singlevar_proc_ptr functionality
  test_pass = test_pass .and. test_singlevar_proc_ptr(test_unit)
  write(test_unit , *) &
      'All singlevar_proc_ptr tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the Brent zero-finding and 
  ! minimization functionality
  test_pass = test_pass .and. test_brent(test_unit)
  write(test_unit , *) &
      'All Brent tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the Moore-Penrose generalized inverse
  test_pass = test_pass .and. test_moore_penrose_inverse(test_unit)
  write(test_unit , *) &
      'All Moore-Penrose tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the computation of the A matrix
  test_pass = test_pass .and. test_compute_a_matrix(test_unit)
  write(test_unit , *) &
      'All compute-A-matrix tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the computation of the A matrix
  test_pass = test_pass .and. test_compute_a_and_b_matrices(test_unit)
  write(test_unit , *) &
      'All compute-A-and-B-matrices tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality to find all zeros
  test_pass = test_pass .and. test_find_all_zeros(test_unit)
  write(test_unit , *) &
      'All find-all-zeros tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the functionality to find relative optima
  test_pass = test_pass .and. test_find_rel_optima(test_unit)
  write(test_unit , *) &
      'All find-rel-optima tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the Chebyshev points functionality
  test_pass = test_pass .and. test_chebyshev_points(test_unit)
  write(test_unit , *) &
      'All chebyshev-points tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the Remez-like step functionality
  test_pass = test_pass .and. test_remez_step(test_unit)
  write(test_unit , *) &
      'All Remez-like-step tests completed; did all pass so far?  ', test_pass
  ! Run the unit tests for the Borsuk lower bound functionality
  test_pass = test_pass .and. test_borsuk_lower_bound(test_unit)
  write(test_unit , *) &
      'All borsuk-lower-bound tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the worst-rho-and-gamma functionality
  test_pass = test_pass .and. test_worst_rho_and_gamma(test_unit)
  write(test_unit , *) &
      'All worst-rho-and-gamma tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality that finds all matrices 
  ! used in constructing the numerically optimal approximation to a given 
  ! (totally positive but potentially asymmetric) kernel
  test_pass = test_pass .and. test_all_num_opt_mats_asymm(test_unit)
  write(test_unit , *) &
      'All all-num-opt-mats-asymm tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality that constructs and evaluates the
  ! rank-$n$ numerically optimal approximating kernel to a given (totally 
  ! positive) kernel
  test_pass = test_pass .and. test_num_opt_rankn_kernel(test_unit)
  write(test_unit , *) &
      'All num-opt-rankn-kernel tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality that finds the numerical
  ! parameters needed to evaluate the rank-$n$ numerically optimal 
  ! approximating kernel for a given (totally positive) kernel; unlike
  ! num_opt_rankn_kernel, the kernel to approximation does ***not***
  ! need to be symmetric, and the output of the procedure being tested
  ! is not a rank-$n$ kernel object -- it is only the vectors and
  ! matrices that mathematically define such a rank-$n$ kernel.
  test_pass = test_pass .and. test_num_opt_rankn_params(test_unit)
  write(test_unit , *) &
      'All num-opt-rankn-params tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality that constructs and evaluates the
  ! rank-$n$ Taylor series approximating a given kernel
  test_pass = test_pass .and. test_taylor_series_rankn_kernel(test_unit)
  write(test_unit , *) &
      'All Taylor-series kernel tests completed; did all pass so far? ', &
      test_pass
  ! Run the unit tests for the functionality that constructs and evaluates the
  ! rank-$n$ singular function series approximating a given kernel
  test_pass = test_pass .and. test_eigen_rankn_kernel(test_unit)
  write(test_unit , *) &
      'All singular function series kernel tests completed; ', &
      'did all pass so far? ', &
      test_pass
  ! Run the unit tests for the functionality that builds and evaluates the 
  ! "integral" (really a sum) of a kernel with respect to a discrete signed 
  ! measure
  test_pass = test_pass .and. test_integral_kernel_disc(test_unit)
  write(test_unit , *) &
      'All integral-kernel-disc tests completed; did all pass so far?  ', &
      test_pass
  ! Run the unit tests for the functionality that builds and evaluates the 
  ! "integral" (really a sum) of a rank-$n$ kernel with respect to a discrete 
  ! signed measure
  test_pass = test_pass .and. test_integral_rankn_disc(test_unit)
  write(test_unit , *) &
      'All integral-rankn-disc tests completed; did all pass so far?  ', &
      test_pass
  ! Run the tests for analytical derivative accuracy
  test_pass = test_pass .and. test_derivs(test_unit)
  write(test_unit , *) &
      'All analytical derivative tests completed; did all pass so far?  ', &
      test_pass
  
  ! Notify of completion and output test_pass logical
  write(test_unit , *) 'All tests completed; did all pass?  ', test_pass
  
  ! Output to the command line
  write(output_unit , *) 'All tests completed; did all pass? ', test_pass
  write(output_unit , *) 'Full output written to test_output.txt.'
  
end program fort_approx