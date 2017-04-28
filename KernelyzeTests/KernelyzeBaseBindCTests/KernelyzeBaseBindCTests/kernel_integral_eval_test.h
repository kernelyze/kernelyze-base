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

#ifndef KERNEL_INTEGRAL_EVAL_C
#define KERNEL_INTEGRAL_EVAL_C

extern "C" void kernel_integral_eval(
	double(*f)(double *x, double *y),
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_eval_pts,
	const double* c_eval_pts,
	double* c_results_of_eval);

void kernel_integral_eval_test();

#endif // KERNEL_INTEGRAL_EVAL_C