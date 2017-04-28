// Copyright (c) 2016, 2017 by Kernelyze LLC
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
// created:	2016-07-07
// updated: 2016-07-07
// updated: 2016-07-14 (added tests for gaussian_expprod 
//          and gaussian_sigma_fixed functionality)
// updated: 2017-04-14 (added function-pointer version
//          of num_opt_rankn_params, named
//          num_opt_rankn_params_ptr)
// updated: 2017-04-21 (added num_opt_rankn_eval)
// updated: 2017-04-24 (added new functionality and deleted deprecated)

#ifndef NUM_OPT_FUNCS_C
#define NUM_OPT_FUNCS_C

extern "C" void num_opt_rankn_params(
	double(*f)(double *x, double *y),
	const double* c_x_lb,
	const double* c_x_ub,
	const double* c_y_lb,
	const double* c_y_ub,
	const int* c_approx_rank,
	const int* c_max_iter,
	const double* c_toler,
	double* c_rho_vec,
	double* c_gamma_vec,
	double* c_v_mat,
	double* c_w_mat,
	double* c_borsuk_lb,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void num_opt_rankn_eval(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const double* c_rho_vec,
	const double* c_gamma_vec,
	const double* c_v_mat,
	const double* c_w_mat,
	const int* c_num_eval_pts,
	const double* c_eval_x_pts,
	const double* c_eval_y_pts,
	double* c_results_of_eval,
	int* c_err_stat,
	char* c_err_msg);

extern "C" void num_opt_rankn_integral_eval(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const double* c_rho_vec,
	const double* c_gamma_vec,
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

extern "C" void num_opt_rankn_integral_coeffs(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const double* c_rho_vec,
	const double* c_gamma_vec,
	const double* c_v_mat,
	const double* c_w_mat,
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_coeffs,
	double* c_coeffs,
	int* c_err_stat,
	char* c_err_msg);

void num_opt_funcs_test();

#endif // NUM_OPT_FUNCS_C
