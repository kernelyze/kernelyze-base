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
// updated: 2017-04-23 (correct printed message)

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "find_rel_optima_test.h"

double sine_func(double *x) {
	return sin(*x);
}

void find_rel_optima_test() {

	int n_grid_pts = 21;
	double *grid_pts;
	grid_pts = (double *)malloc(sizeof(double) * n_grid_pts);
	for (int i = 0; i < n_grid_pts; ++i) {
		grid_pts[i] = -10.0 + (double)i;
	}
	double toler = 1e-15;
	double *result_optima;
	result_optima = (double *)malloc(
		sizeof(double) * 2 * (n_grid_pts - 1));
	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	find_rel_optima(
		result_optima,
		&sine_func,
		&n_grid_pts,
		grid_pts,
		&toler,
		&c_err_stat,
		c_err_msg);

	printf("finding rel optima of sin(x) from -10.0 to 10.0:\n");
	printf("%d optima expected, error code: %d, error message: %s\n",
		n_grid_pts - 1, c_err_stat, c_err_msg);
	for (int i = 0; i < n_grid_pts - 1; ++i) {
		// Keep in mind column-major ordering in the array
		// result_optima, which is (n_grid_pts - 1) rows by
		// 2 columns.
		printf("opt number %d: x_opt %f, -abs(f_opt) %f\n", 
			i, result_optima[i], result_optima[i + n_grid_pts - 1]);
	}

	free(grid_pts);
	free(result_optima);
}