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

#ifndef EIGEN_FUNCS_C
#define EIGEN_FUNCS_C

extern "C" void eigen_rankn_params(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const int* c_num_x_pts,
	const int* c_num_y_pts,
	const double* c_x_points,
	const double* c_y_points,
	const double* c_x_weights,
	const double* c_y_weights,
	double* c_x_func_wts,
	double* c_y_func_wts,
	double* c_trunc_sing_vals,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void eigen_rankn_eval(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const int* c_num_x_pts,
	const int* c_num_y_pts,
	const double* c_x_points,
	const double* c_y_points,
	double* c_x_func_wts,
	double* c_y_func_wts,
	double* c_trunc_sing_vals,
	const int* c_num_eval_pts,
	const double* c_eval_x_pts,
	const double* c_eval_y_pts,
	double* c_results_of_eval,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void eigen_rankn_integral_eval(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const int* c_num_x_pts,
	const int* c_num_y_pts,
	const double* c_x_points,
	const double* c_y_points,
	double* c_x_func_wts,
	double* c_y_func_wts,
	double* c_trunc_sing_vals,
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_eval_pts,
	const double* c_eval_pts,
	double* c_results_of_eval,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void eigen_rankn_integral_coeffs(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const int* c_num_x_pts,
	const int* c_num_y_pts,
	const double* c_x_points,
	const double* c_y_points,
	double* c_x_func_wts,
	double* c_y_func_wts,
	double* c_trunc_sing_vals,
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_coeffs,
	const double* c_coeffs,
	int* c_err_stat,
	char* c_err_msg);

void eigen_approx_test();

#endif // EIGEN_FUNCS_C
