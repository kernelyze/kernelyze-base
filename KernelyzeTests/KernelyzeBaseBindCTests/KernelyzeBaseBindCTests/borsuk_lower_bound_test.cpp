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
// created:	2017-04-14
// updated: 2017-04-21 (put function for function pointer
//          test in separate header)

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "test_function_ptrs.h"
#include "borsuk_lower_bound_test.h"

void borsuk_lower_bound_test() {

	int n_grid_pts = 40;
	double *grid_pts;
	grid_pts = (double *)malloc(sizeof(double) * n_grid_pts);
	chebyshev_points(&n_grid_pts, grid_pts);

	int size_of_rho = 5;
	double* rho_vec;
	rho_vec = (double *)malloc(sizeof(double) * size_of_rho);
	chebyshev_points(&size_of_rho, rho_vec);

	bool is_over_x = true;
	double tolerance = 1e-15;
	int max_iter = 100;

	double *coeffs			= (double *)malloc(sizeof(double) * size_of_rho);
	double *nodes			= (double *)malloc(sizeof(double) * size_of_rho);
	double *errors_at_nodes = (double *)malloc(sizeof(double) * size_of_rho);

	double discrepancy;
	int num_iter;

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	borsuk_lower_bound(
		&expprod,
		&size_of_rho,
		rho_vec,
		&is_over_x,
		&tolerance,
		&max_iter,
		&n_grid_pts,
		grid_pts,
		coeffs,
		nodes,
		errors_at_nodes,
		&discrepancy,
		&num_iter,
		&c_err_stat,
		c_err_msg);

	printf("Borsuk lower bound, rank 4 approx (5-elt rho) of exp prod:\n");
	printf("Discrepancy: %f\n", discrepancy);
	for (int i = 0; i < size_of_rho; ++i) {
		printf("Index %d: coeff: %f node: %f error at node: %f\n",
			i, coeffs[i], nodes[i], errors_at_nodes[i]);
	}

	free(grid_pts);
	free(rho_vec);
	free(coeffs);
	free(nodes);
	free(errors_at_nodes);
}