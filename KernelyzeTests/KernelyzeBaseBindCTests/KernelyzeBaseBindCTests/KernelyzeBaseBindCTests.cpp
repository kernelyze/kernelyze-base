// KernelyzeBaseBindCTests.cpp : Defines the entry point for the console application.
//
// Copyright (c) 2015, 2016, 2017 by Kernelyze LLC
// Author: Thomas A. Knox
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// created:   2015-06-02
// updated:   2015-06-03
// updated:   2016-02-28 (added build_eval_num_opt functionality)
// updated:   2016-05-14 (added thread-safe multivar store functionality)
// updated:   2016-05-15 (better approach to thread-safe multivar store tests)
// updated:	  2016-07-07 (added num_opt_funcs testing)
// updated:   2016-09-12 (added trivar store testing)
// updated:   2017-04-14 (added testing for Brent zero-finding,
//            Brent minimization, find-all-zeros, find-rel-optima,
//            Borsuk lower bounds, and worst rho and gamma arrays)
// updated:   2017-04-21 (added Taylor-series and singular-function
//            approximation tests using function pointers)
// updated:   2017-04-23 (added option-formulae, kernel-integral,
//            matmul, and linear-solver tests)
// updated:   2017-04-24 (rename to "Base" rather than "Core" to reflect
//            different solution)

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "mp_inverse.h"
#include "num_opt_funcs_test.h"
#include "brent_test.h"
#include "find_all_zeros_test.h"
#include "find_rel_optima_test.h"
#include "borsuk_lower_bound_test.h"
#include "worst_rho_and_gamma_test.h"
#include "taylor_approx_test.h"
#include "eigen_approx_test.h"
#include "kernel_integral_eval_test.h"
#include "option_formulae_test.h"
#include "linear_solve_test.h"
#include "matmul_test.h"

int _tmain(int argc, _TCHAR* argv[])
{
  int matsize = 3;
  double *to_invert = (double *)malloc(matsize * matsize * sizeof(double));
  double *to_return = (double *)malloc(matsize * matsize * sizeof(double));
  // Fill in the arrays
  for (int i = 0; i < matsize; ++i) {
    for (int j = 0; j < matsize; ++j) {
      // Initialize the to_invert array to the identity
      to_invert[matsize * j + i] = (i == j) ? double(1.0) : double(0.1);
      // Initialize the to_return array to -99; this is just
      // a convenient check to be sure that this array does
      // actually get overwritten as desired, I don't really
      // need to initialize it at all.
      to_return[matsize * j + i] = double(-99.0);
    }
  }
  // Make the call
  mp_inverse(to_return, to_invert, matsize, matsize, 1e-14);
  // Output first the matrix to invert, then the results
  printf("The matrix to invert is:\n");
  for (int i = 0; i < matsize; ++i) {
    for (int j = 0; j < matsize; ++j) {
      printf("%f ", to_invert[matsize * i + j]);
    }
    printf("\n");
  }
  printf("The result is:\n");
  for (int i = 0; i < matsize; ++i) {
    for (int j = 0; j < matsize; ++j) {
      printf("%f ", to_return[matsize * i + j]);
    }
    printf("\n");
  }
  free(to_invert);
  free(to_return);

  brent_test();

  find_all_zeros_test();

  find_rel_optima_test();

  matmul_test();

  linear_solve_test();

  option_formulae_test();

  kernel_integral_eval_test();

  borsuk_lower_bound_test();

  worst_rho_and_gamma_test();

  taylor_approx_test();

  eigen_approx_test();

  num_opt_funcs_test();

  return 0;
}

