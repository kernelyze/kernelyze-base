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

void matmul_test() {

	static const int m = 15;
	static const int n = 4;
	static const int p = 3;

	double amat[m * n];
	double bmat[n * p];
	double xmat[m * p];
	double test_should_be_xmat[m * p];

	// Efficiency is not critical in this small test,
	// hence no stride optimizations
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			amat[j * m + i] = exp((double)((i - j)) / 20.0);
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < p; ++j) {
			bmat[j * n + i] = sin((double)(i + j));
		}
	}

	matmul(&m, &n, &p, amat, bmat, xmat);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < p; ++j) {
			test_should_be_xmat[j * m + i] = 0.0;
			for (int k = 0; k < n; ++k) {
				test_should_be_xmat[j * m + i] +=
					amat[k * m + i] * bmat[j * n + k];
			}
		}
	}

	for (int i = 0; i < m; ++i) {
		printf("Row %d difference in A * B from library and by hand, should be zero:\n", i);
		for (int j = 0; j < p; ++j) {
			printf("%f ", test_should_be_xmat[j * m + i] - xmat[j * m + i]);
		}
		printf("\n");
	}

}