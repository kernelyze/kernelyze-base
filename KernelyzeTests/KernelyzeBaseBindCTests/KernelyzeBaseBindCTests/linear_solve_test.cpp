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
// created:	2017-04-23

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "matmul_test.h"
#include "linear_solve_test.h"

void linear_solve_test() {

	static const int n = 15;
	static const int nrhs = 4;

	double amat[n * n];
	double bmat[n * nrhs];
	double xmat[n * nrhs];
	double test_should_be_bmat[n * nrhs];

	// Efficiency is not critical in this small test,
	// hence no stride optimizations
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			amat[j * n + i] = 1.0 / (1.0 + (double)( (i - j) * (i - j) ));
		}
		for (int j = 0; j < nrhs; ++j) {
			bmat[j * n + i] = sin((double)(i + j));
		}
	}

	int info;

	linear_solve(&n, amat, &nrhs, bmat, xmat, &info);

	matmul(&n, &n, &nrhs, amat, xmat, test_should_be_bmat);

	for (int i = 0; i < n; ++i) {
		printf("Row %d difference in A * X_solved and B, should be zero:\n", i);
		for (int j = 0; j < nrhs; ++j) {
			printf("%f ", test_should_be_bmat[j * n + i] - bmat[j * n + i]);
		}
		printf("\n");
	}

}