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
// created:   2017-04-14
// updated:   2017-04-14

#ifndef BORSUK_LOWER_BOUND_FOR_C
#define BORSUK_LOWER_BOUND_FOR_C

extern "C" void chebyshev_points(
	int* n,
	double* pts);

extern "C" void borsuk_lower_bound(
	double (*f)(double *x, double *y),
	int* size_of_rho,
	double* rho_vec, 
	bool* is_over_x,
	double* tolerance,
	int* max_iter,
	int* n_grid_pts,
	double* grid,
	double* coeffs,
	double* nodes,
	double* errors_at_nodes,
	double* discrepancy,
	int* num_iter,
	int* err_stat,
	char* err_msg);

void borsuk_lower_bound_test();

#endif // BORSUK_LOWER_BOUND_FOR_C

