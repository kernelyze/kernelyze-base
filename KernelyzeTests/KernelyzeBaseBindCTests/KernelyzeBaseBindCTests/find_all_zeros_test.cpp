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
// updated: 2017-04-14

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "find_all_zeros_test.h"

double sin_func(double *x) {
	return sin(*x);
}

void find_all_zeros_test() {
	
	int n_grid_pts = 21;
	double *grid_pts;
	grid_pts = (double *)malloc(sizeof(double) * n_grid_pts);
	for (int i = 0; i < n_grid_pts; ++i) {
		grid_pts[i] = -10.0 + (double)i;
	}
	double toler = 1e-15;
	int max_num_zeros = 7;
	int num_zeros;
	double *result_zeros;
	result_zeros = (double *)malloc(sizeof(double) * max_num_zeros);
	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	find_all_zeros(
		&max_num_zeros,
		&num_zeros,
		result_zeros,
		&sin_func,
		&n_grid_pts,
		grid_pts,
		&toler,
		&c_err_stat,
		c_err_msg);

	printf("finding all zeros of sin(x) from -10.0 to 10.0:\n");
	printf("%d zeros found, error code: %d, error message: %s\n", 
		num_zeros, c_err_stat, c_err_msg);
	for (int i = 0; i < num_zeros; ++i) {
		printf("zero %d: %f\n", i, result_zeros[i]);
	}

	free(grid_pts);
	free(result_zeros);
}