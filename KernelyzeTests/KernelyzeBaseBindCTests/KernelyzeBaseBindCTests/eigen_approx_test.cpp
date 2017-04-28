// Copyright (c) 2017 by Kernelyze LLC
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
// created:	2017-04-21
// updated: 2017-04-23 (added integral eval and coeffs)

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "test_function_ptrs.h"
#include "eigen_approx_test.h"

void eigen_approx_test() {
	// Index variables
	int i, j;
	// Arguments not requiring dynamic allocation
	int err_stat[1] = { 0 };
	char err_msg[256];
	
	static const int approx_rank = 3;
	static const int num_eigen_x_pts = 50;
	static const int num_eigen_y_pts = 50;

	double x_points[num_eigen_x_pts];
	double y_points[num_eigen_y_pts];
	double x_weights[num_eigen_x_pts];
	double y_weights[num_eigen_y_pts];

	double x_func_wts[num_eigen_y_pts * approx_rank];
	double y_func_wts[num_eigen_x_pts * approx_rank];
	double trunc_sing_vals[approx_rank];

	// Fill x and y points and weights
	// Use the square [0, 2] \times [-1, 1] for the domain
	for (i = 0; i < num_eigen_x_pts; ++i) {
		x_points[i] = 2.0 * ((double)i) / ((double)(num_eigen_x_pts - 1));
		x_weights[i] = 2.0 / ((double)num_eigen_x_pts);
	}
	for (i = 0; i < num_eigen_y_pts; ++i) {
		y_points[i] = -1.0 + 2.0 * ((double)i) / ((double)(num_eigen_y_pts - 1));;
		y_weights[i] = 2.0 / ((double)num_eigen_y_pts);
	}

	// Call into the C binding of the Fortran procedure to get the
	// numerically-optimal rank-$n$ parameters for this kernel
	eigen_rankn_params(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_weights,
		y_weights,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		err_stat,
		err_msg);

	printf("Singular func approx params from func ptr to exp prod kernel:\n");
	printf("The x func wts matrix is:\n");
	for (i = 0; i < num_eigen_x_pts; ++i) {
		for (j = 0; j < approx_rank; ++j) {
			printf("%f ", x_func_wts[num_eigen_x_pts * j + i]);
		}
		printf("\n ");
	}
	printf("The y func wts matrix is:\n");
	for (i = 0; i < num_eigen_y_pts; ++i) {
		for (j = 0; j < approx_rank; ++j) {
			printf("%f ", y_func_wts[num_eigen_y_pts * j + i]);
		}
		printf("\n ");
	}

	// Next test evaluating the rank-$n$ Taylor-series
	// approximation using the parameters produced above
	// from the function pointer
	static const int eval_size = 100;
	double test_eval_x_pts[eval_size];
	double test_eval_y_pts[eval_size];
	double results_eval[eval_size];
	double x_lb[1] = { 0.0 };
	double x_ub[1] = { 2.0 };
	double y_lb[1] = { -1.0 };
	double y_ub[1] = { 1.0 };

	for (i = 0; i < eval_size; ++i) {
		test_eval_x_pts[i] = *x_lb +
			(((double)i) / ((double)(eval_size - 1)))
			* (*x_ub - *x_lb);
		test_eval_y_pts[i] = *y_lb +
			(((double)i) / ((double)(eval_size - 1)))
			* (*y_ub - *y_lb);
	}

	eigen_rankn_eval(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		&eval_size,
		test_eval_x_pts,
		test_eval_y_pts,
		results_eval,
		err_stat,
		err_msg);

	printf("Results from singular-function approximation via function pointer of exp-prod kernel.\n");
	for (i = 0; i < eval_size; ++i) {
		double exact_result = expprod(&test_eval_x_pts[i], &test_eval_y_pts[i]);
		printf("At (x, y) of (%f, %f) singular-func-approx from params is %f and exact is %f\n",
			test_eval_x_pts[i], test_eval_y_pts[i], results_eval[i], exact_result);
	}

	static const int num_integral_pts = 21;
	double integral_wts[num_integral_pts];
	double integral_pts[num_integral_pts];
	// For use in making the rank-$n$ kernel eval give intermediate results
	// in the test calculation
	double dummy_eval_pts[num_integral_pts];
	double dummy_eval_res[num_integral_pts];

	static const int num_eval_pts = 11;
	double eval_pts[num_eval_pts];
	double eval_results[num_eval_pts];
	double eval_results_by_hand[num_eval_pts];

	bool is_over_x = true;

	// Evenly-spaced points from -1.0 to 1.0, with
	// equal weights
	for (int i = 0; i < num_integral_pts; ++i) {
		integral_wts[i] = 1.0 / ((double)num_integral_pts);
		integral_pts[i] = -1.0 + 2.0 * ((double)i) / ((double)(num_integral_pts - 1));
	}

	// Test at similarly evenly-spaced points
	for (int i = 0; i < num_eval_pts; ++i) {
		eval_pts[i] = -1.0 + 2.0 * ((double)i) / ((double)(num_eval_pts - 1));
	}

	// Compute the integral (really a sum) using the library:
	eigen_rankn_integral_eval(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		&is_over_x,
		&num_integral_pts,
		integral_pts,
		integral_wts,
		&num_eval_pts,
		eval_pts,
		eval_results,
		err_stat,
		err_msg);

	for (int i = 0; i < num_eval_pts; ++i) {
		eval_results_by_hand[i] = 0.0;
		for (int j = 0; j < num_integral_pts; ++j) {
			dummy_eval_pts[j] = eval_pts[i];
		}
		eigen_rankn_eval(
			&expprod,
			&approx_rank,
			&num_eigen_x_pts,
			&num_eigen_y_pts,
			x_points,
			y_points,
			x_func_wts,
			y_func_wts,
			trunc_sing_vals,
			&num_integral_pts,
			integral_pts,
			dummy_eval_pts,
			dummy_eval_res,
			err_stat,
			err_msg);
		for (int j = 0; j < num_integral_pts; ++j) {
			eval_results_by_hand[i] +=
				integral_wts[j] * dummy_eval_res[j];
		}
		printf("Singular func. rank-n integral (sum) from library at %f is %f; by rank-n kernel eval it is %f\n",
			eval_pts[i], eval_results[i], eval_results_by_hand[i]);
	}

	// Obtain the coefficients of the rank-$n$ integral
	static const int num_coeffs = 3; // must equal approx_rank value
	double integral_coeffs[num_coeffs];
	eigen_rankn_integral_coeffs(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		&is_over_x,
		&num_integral_pts,
		integral_pts,
		integral_wts,
		&num_coeffs,
		integral_coeffs,
		err_stat,
		err_msg);

	for (int i = 0; i < num_coeffs; ++i) {
		printf("Coefficient %d of the singular func. rank-n integral is %f\n",
			i, integral_coeffs[i]);
	}

	// Check the linearity of the rank-$n$ integral coeffs
	static const int num_integral2_pts = 4;
	double integral2_wts[num_integral2_pts];
	double integral2_pts[num_integral2_pts];

	// The line below only works if integral and integral2 have disjoint sets of
	// points; obviously, the linearity of the coefficients is more general than this.
	static const int num_integral3_pts = num_integral_pts + num_integral2_pts;
	double integral3_wts[num_integral3_pts];
	double integral3_pts[num_integral3_pts];

	// For integral2, evenly-spaced points from -2.0 to -1.5, with negative weights
	for (int i = 0; i < num_integral2_pts; ++i) {
		integral2_wts[i] = -1.0 / ((double)num_integral_pts);
		integral2_pts[i] = -2.0 + 0.5 * ((double)i) / ((double)(num_integral_pts - 1));
	}

	// Now take the integral that is -0.5 * integral + 2.0 * integral2
	for (int i = 0; i < num_integral3_pts; ++i) {
		if (i < num_integral2_pts) {
			integral3_wts[i] = 2.0 * integral2_wts[i];
			integral3_pts[i] = integral2_pts[i];
		}
		else {
			int tempi = i - num_integral2_pts;
			integral3_wts[i] = -0.5 * integral_wts[tempi];
			integral3_pts[i] = integral_pts[tempi];
		}
	}

	// Compute the coefficients of integral2 and of integral3 (which
	// is -0.5 * integral + 2.0 * integral2)
	double integral2_coeffs[num_coeffs];
	double integral3_coeffs[num_coeffs];
	eigen_rankn_integral_coeffs(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		&is_over_x,
		&num_integral2_pts,
		integral2_pts,
		integral2_wts,
		&num_coeffs,
		integral2_coeffs,
		err_stat,
		err_msg);
	eigen_rankn_integral_coeffs(
		&expprod,
		&approx_rank,
		&num_eigen_x_pts,
		&num_eigen_y_pts,
		x_points,
		y_points,
		x_func_wts,
		y_func_wts,
		trunc_sing_vals,
		&is_over_x,
		&num_integral3_pts,
		integral3_pts,
		integral3_wts,
		&num_coeffs,
		integral3_coeffs,
		err_stat,
		err_msg);

	// Check that integral3_coeffs = integral_coeffs + integral2_coeffs
	// This is a test of the property that rank-$n$ integral coefficients
	// should be linear in the underlying integration measures, where
	// multiplying a weighted sum (discrete integration measure) by a
	// scalar means multiplying each of its weights by that measure
	// and adding two weighted sums (discrete integration measures) means
	// taking the union of their sets of support points and attaching to
	// each support point the sum of the weights that point gets in each
	// of the weighted sums to be combined (if a point exists in only one
	// of the weighted sums, its weight in the other is zero).
	printf("Test of the linearity of singular func. rank-n integral coefficients:\n");
	for (int i = 0; i < num_coeffs; ++i) {
		printf("Coeff %d of integral3 is %f, should be %f\n",
			i,
			integral3_coeffs[i],
			-0.5 * integral_coeffs[i] + 2.0 * integral2_coeffs[i]);
	}

}
