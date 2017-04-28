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
// updated: 2017-04-21 (put function for function ptr
//          into separate header)

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "test_function_ptrs.h"
#include "worst_rho_and_gamma_test.h"

void worst_rho_and_gamma_test() {

	int num_terms = 9;
	double tolerance = 1e-15;
	int max_iter = 100;

	double x_lb =  0.0;
	double x_ub =  5.0;
	double y_lb = -1.0;
	double y_ub =  0.5;

	double *rho_vec = (double *)malloc(sizeof(double) * num_terms);
	double *gamma_vec = (double *)malloc(sizeof(double) * num_terms);
	
	double *coeffs_rho = (double *)malloc(sizeof(double) * num_terms);
	double *coeffs_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *nodes_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *errors_at_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *errors_at_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);

	double* a_matrix_rho = (double *)malloc(
		sizeof(double) * num_terms * num_terms);
	double* a_matrix_gamma = (double *)malloc(
		sizeof(double) * num_terms * num_terms);

	double *b_coeffs_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_coeffs_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *b_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *b_errors_at_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_errors_at_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);

	double discrepancy;
	int num_iter;

	double ker_discrep;
	int ker_num_iter;

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	worst_rho_and_gamma(
		&expprod,
		&x_lb,
		&x_ub,
		&y_lb,
		&y_ub,
		&num_terms,
		&tolerance,
		&max_iter,
		rho_vec,
		gamma_vec,
		coeffs_rho,
		coeffs_gamma,
		a_matrix_rho,
		a_matrix_gamma,
		nodes_rho,
		nodes_gamma,
		errors_at_nodes_rho,
		errors_at_nodes_gamma,
		b_coeffs_rho,
		b_coeffs_gamma,
		b_nodes_rho,
		b_nodes_gamma,
		b_errors_at_nodes_rho,
		b_errors_at_nodes_gamma,
		&discrepancy,
		&num_iter,
		&ker_discrep,
		&ker_num_iter,
		&c_err_stat,
		c_err_msg);

	if (c_err_stat != 0) {
		printf("Test failed with message: %s\n", c_err_msg);
	}
	printf("Worst rho and gamma, rank 4 approx (5-elt rho) of exp prod:\n");
	printf("x range from %f to %f\n", x_lb, x_ub);
	printf("y range from %f to %f\n", y_lb, y_ub);
	printf("discrepancy: %12.10e\n", discrepancy);
	printf("kernel discrepancy: %12.10e\n", ker_discrep);
	printf("number of iterations: %d\n", num_iter);
	printf("kernel number of iterations: %d\n", ker_num_iter);
	
	printf("rho primary outputs\n");
	for (int i = 0; i < num_terms; ++i) {
		printf("Index %d: rho: %12.10e coeff: %12.10e node: %12.10e error at node: %12.10e\n",
			i, 
			rho_vec[i], 
			coeffs_rho[i], 
			nodes_rho[i], 
			errors_at_nodes_rho[i]);
	}

	printf("gamma primary outputs\n");
	for (int i = 0; i < num_terms; ++i) {
		printf("Index %d: gamma: %12.10e coeff: %12.10e node: %12.10e error at node: %12.10e\n",
			i,
			gamma_vec[i],
			coeffs_gamma[i],
			nodes_gamma[i],
			errors_at_nodes_gamma[i]);
	}

	printf("rho Borsuk outputs\n");
	for (int i = 0; i < num_terms; ++i) {
		printf("Index %d: b_coeff: %12.10e b_node: %12.10e b_error at node: %12.10e\n",
			i,
			b_coeffs_rho[i],
			b_nodes_rho[i],
			b_errors_at_nodes_rho[i]);
	}

	printf("gamma Borsuk outputs\n");
	for (int i = 0; i < num_terms; ++i) {
		printf("Index %d: b_coeff: %12.10e b_node: %12.10e b_error at node: %12.10e\n",
			i,
			b_coeffs_gamma[i],
			b_nodes_gamma[i],
			b_errors_at_nodes_gamma[i]);
	}

	free(rho_vec);
	free(gamma_vec);
	free(coeffs_rho);
	free(coeffs_gamma);
	free(nodes_rho);
	free(nodes_gamma);
	free(errors_at_nodes_rho);
	free(errors_at_nodes_gamma);
	free(a_matrix_rho);
	free(a_matrix_gamma);
	free(b_coeffs_rho);
	free(b_coeffs_gamma);
	free(b_nodes_rho);
	free(b_nodes_gamma);
	free(b_errors_at_nodes_rho);
	free(b_errors_at_nodes_gamma);

}