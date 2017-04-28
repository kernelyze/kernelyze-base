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

#ifndef TAYLOR_FUNCS_C
#define TAYLOR_FUNCS_C

extern "C" void taylor_rankn_params(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const double* c_x_center,
	const double* c_y_center,
	double* c_v_mat,
	double* c_w_mat,
	int* c_err_stat,
	char* c_err_msg,
	double* fin_diff_step);

extern "C" void taylor_rankn_eval(
	const int* c_approx_rank,
	const double* c_x_center,
	const double* c_y_center,
	const double* c_v_mat,
	const double* c_w_mat,
	const int* c_num_eval_pts,
	const double* c_eval_x_pts,
	const double* c_eval_y_pts,
	double* c_results_of_eval,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void taylor_rankn_integral_eval(
	const int* c_approx_rank,
	const double* c_x_center,
	const double* c_y_center,
	const double* c_v_mat,
	const double* c_w_mat,
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_eval_pts,
	const double* c_eval_pts,
	double* c_results_of_eval,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void taylor_rankn_integral_coeffs(
	const int* c_approx_rank,
	const double* c_x_center,
	const double* c_y_center,
	const double* c_v_mat,
	const double* c_w_mat,
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_coeffs,
	const double* c_coeffs,
	int* c_err_stat,
	char* c_err_msg);

void taylor_approx_test();

#endif // TAYLOR_FUNCS_C
