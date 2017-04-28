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
// created:   2017-04-26 (based on earlier separate files)

#ifndef KERNLCAML_STUBS_H
#define KERNLCAML_STUBS_H

#ifndef CAML_NAME_SPACE
#define CAML_NAME_SPACE
#endif

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <caml/mlvalues.h>
#include <caml/callback.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/bigarray.h>

// Approach: give declaration of C interface and then declaration
// of OCaml primitive stub

// Zeroin
double zeroin(
	double* a,
	double* b,
	double(*f)(double *x),
	double* toler);

CAMLprim value caml_zeroin(
	value vf,
	value va,
	value vb,
	value vtol);

// Brent zero
void brent_zero(
	double* a,
	double* b,
	double(*f)(double *x),
	double* toler,
	int* max_iter,
	int* neval,
	int* err_code,
	double* x_zero,
	double* f_zero);

CAMLprim value caml_brent_zero(
	value vf,
	value va,
	value vb,
	value vtol,
	value vmax_iter);

// f_min
void f_min(
	double* a,
	double* b,
	double(*f)(double *x),
	double* tol,
	double* x_opt,
	double* f_opt,
	int* err_stat,
	char* err_msg);

CAMLprim value caml_f_min(
	value vf,
	value va,
	value vb,
	value vtol);

// brent_min
void brent_min(
	double* a,
	double* b,
	double(*f)(double *x),
	double* toler,
	int* max_iter,
	int* neval,
	int* err_code,
	double* x_min,
	double* f_min);

CAMLprim value caml_brent_min(
	value vf,
	value va,
	value vb,
	value vtol,
	value vmax_iter);

// find_all_zeros
void find_all_zeros(
	int* max_num_zeros,
	int* num_zeros,
	double* result_zeros,
	double(*f)(double *x),
	int* n_grid_pts,
	double* grid,
	double* toler,
	int* err_stat,
	char* err_msg);

CAMLprim value caml_find_all_zeros(
	value vf,
	value vgrid,
	value vtol);

// find_rel_optima
void find_rel_optima(
	double* result_optima,
	double(*f)(double *x),
	int* n_grid_pts,
	double* grid,
	double* toler,
	int* err_stat,
	char* err_msg);

CAMLprim value caml_find_rel_optima(
	value vf,
	value vgrid,
	value vtol);

// chebyshev_points
void chebyshev_points(
	int* n,
	double* pts);

CAMLprim value caml_cheb_pts(
	value vn);

// matmul
void matmul(
	const int* c_m,
	const int* c_n,
	const int* c_p,
	const double* c_a,
	const double* c_b,
	double* c_d);

CAMLprim value caml_matmul(
	value va,
	value vb);

// linear_solve
void linear_solve(
	const int* c_n,
	const double* c_a,
	const int* c_nrhs,
	const double* c_b,
	double* c_x,
	int* c_info);

CAMLprim value caml_linear_solve(
	value va,
	value vb);

// Moore-Penrose inverse
void mp_inverse(
	double* to_return, 
	double* to_invert, 
	int nrows, 
	int ncols, 
	double toler);

CAMLprim value caml_mp_inverse(
	value va,
	value vtol);

// Black (1976) option pricing formula
double black_formula(
	const bool* c_is_call,
	const double* c_strike,
	const double* c_forward,
	const double* c_vol,
	const double* c_disc_fac);

CAMLprim value caml_black_formula(
	value vcall,
	value vstrike,
	value vforward,
	value vvol,
	value vdiscfac);

// Bachelier (or normal-distribution) option pricing formula
double normopt_formula(
	const bool* c_is_call,
	const double* c_strike,
	const double* c_forward,
	const double* c_vol,
	const double* c_disc_fac);

CAMLprim value caml_normopt_formula(
	value vcall,
	value vstrike,
	value vforward,
	value vvol,
	value vdiscfac);

// Borsuk lower bound
void borsuk_lower_bound(
	double(*f)(double *x, double *y),
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

// Two OCaml stubs because there are more than
// 5 arguments for borsuk_lower_bound's OCaml
// version.
CAMLprim value caml_borsuk_lower_bound_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_borsuk_lower_bound_native(
	value vf,
	value vrho_vec,
	value voverx,
	value vtol,
	value vmax_iter,
	value vgrid);

// Find rho and gamma that make the Borsuk lower bound large
void worst_rho_and_gamma(
	double(*f)(double *x, double *y),
	double* x_lb,
	double* x_ub,
	double* y_lb,
	double* y_ub,
	int* num_terms,
	double* tolerance,
	int* max_iter,
	double* rho_vec,
	double* gamma_vec,
	double* coeffs_rho,
	double* coeffs_gamma,
	double* a_matrix_rho,
	double* a_matrix_gamma,
	double* nodes_rho,
	double* nodes_gamma,
	double* errors_at_nodes_rho,
	double* errors_at_nodes_gamma,
	double* b_coeffs_rho,
	double* b_coeffs_gamma,
	double* b_nodes_rho,
	double* b_nodes_gamma,
	double* b_errors_at_nodes_rho,
	double* b_errors_at_nodes_gamma,
	double* discrepancy,
	int* num_iter,
	double* ker_discrep,
	int* ker_num_iter,
	int* err_stat,
	char* err_msg);

// Two OCaml stubs because there are more than
// 5 arguments for worst_rho_and_gamma's OCaml
// version.
CAMLprim value caml_worst_rho_and_gamma_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_worst_rho_and_gamma_native(
	value vf,
	value vx_lb,
	value vx_ub,
	value vy_lb,
	value vy_ub,
	value vnumterms,
	value vtol,
	value vmax_iter);

// Compute the integral (really a discrete sum) over one
// variable of a two-variable function (a "kernel")
void kernel_integral_eval(
	double(*f)(double *x, double *y),
	const bool* c_integral_is_over_x,
	const int* c_num_integral_pts,
	const double* c_integral_pts,
	const double* c_integral_wts,
	const int* c_num_eval_pts,
	const double* c_eval_pts,
	double* c_results_of_eval);

CAMLprim value caml_kernel_integral_eval(
	value vf,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts);

// Compute the parameters of a Taylor-series approximation
// to a two-variable function
void taylor_rankn_params(
	double(*f)(double *x, double *y),
	const int* c_approx_rank,
	const double* c_x_center,
	const double* c_y_center,
	double* c_v_mat,
	double* c_w_mat,
	int* c_err_stat,
	char* c_err_msg,
	double* fin_diff_step);

CAMLprim value caml_taylor_rankn_params(
	value vf,
	value vrank,
	value vx_center,
	value vy_center,
	value vfin_diff_step);

// Evaluate the Taylor-series approximation to a two-
// variable function
void taylor_rankn_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for taylor_rankn_eval's OCaml
// version.
CAMLprim value caml_taylor_rankn_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_taylor_rankn_eval_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value veval_x_pts,
	value veval_y_pts);

// Evaluate the integral (really a weighted sum) over 
// one variable of a Taylor-series approximation to a
// two-variable function
void taylor_rankn_integral_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for taylor_rankn_integral_eval's
// OCaml version.
CAMLprim value caml_taylor_rankn_integral_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_taylor_rankn_integral_eval_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts);

// If the Taylor-series approximation to a function
// of two variables has rank $n$, then any integral
// (or sum, here) over one of the two variables can
// be fully described by $n$ coefficients.
void taylor_rankn_integral_coeffs(
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

// Two OCaml stubs because there are more than
// 5 arguments for taylor_rankn_integral_coeff's
// OCaml version.
CAMLprim value caml_taylor_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_taylor_rankn_integral_coeffs_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts);

// Compute the parameters of an approximation
// to a two-variable function that uses
// a truncated series of singular functions
// (the continuous analogs of singular vectors
// in the singular value decomposition).
void eigen_rankn_params(
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

// Two OCaml stubs because there are more than
// 5 arguments for eigen_rankn_params'
// OCaml version.
CAMLprim value caml_eigen_rankn_params_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_eigen_rankn_params_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_wts,
	value vy_wts);

// Evaluate the singular-function approximation to a two-
// variable function, given pre-computed parameters.
void eigen_rankn_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for eigen_rankn_eval's
// OCaml version.
CAMLprim value caml_eigen_rankn_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_eigen_rankn_eval_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_funcwts,
	value vy_funcwts,
	value vtrunc_svals,
	value veval_x_pts,
	value veval_y_pts);

// Evaluate the integral(really a weighted sum) over
// one variable of a singular-function approximation to a
// two-variable function
void eigen_rankn_integral_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for eigen_rankn_integral_eval's
// OCaml version.
CAMLprim value caml_eigen_rankn_integral_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_eigen_rankn_integral_eval_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_funcwts,
	value vy_funcwts,
	value vtrunc_svals,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts);

// If the singular-function approximation to a function
// of two variables has rank $n$, then any integral
// (or sum, here) over one of the two variables can
// be fully described by $n$ coefficients.
void eigen_rankn_integral_coeffs(
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

// Two OCaml stubs because there are more than
// 5 arguments for eigen_rankn_integral_coeffs'
// OCaml version.
CAMLprim value caml_eigen_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_eigen_rankn_integral_coeffs_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_funcwts,
	value vy_funcwts,
	value vtrunc_svals,
	value voverx,
	value vintegral_pts,
	value vintegral_wts);

// Compute the parameters of a numerically-optimal
// rank-$n$ approximation to a two-variable function,
// where approximation quality is assessed in the
// supremum norm over a given rectangle.
void num_opt_rankn_params(
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

// Two OCaml stubs because there are more than
// 5 arguments for num_opt_rankn_params' OCaml
// version.
CAMLprim value caml_num_opt_rankn_params_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_num_opt_rankn_params_native(
	value vf,
	value vx_lb,
	value vx_ub,
	value vy_lb,
	value vy_ub,
	value vrank,
	value vtol,
	value vmax_iter);

// Evaluate the numerically-optimal rank-$n$
// approximation to a two-variable function.
void num_opt_rankn_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for num_opt_rankn_eval's OCaml
// version.
CAMLprim value caml_num_opt_rankn_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_num_opt_rankn_eval_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value veval_x_pts,
	value veval_y_pts);

// Evaluate the integral (here, a sum) over one
// variable of a numerically-optimal rank-$n$
// approximation to a two-variable function.
void num_opt_rankn_integral_eval(
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

// Two OCaml stubs because there are more than
// 5 arguments for num_opt_rankn_integral_eval's OCaml
// version.
CAMLprim value caml_num_opt_rankn_integral_eval_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_num_opt_rankn_integral_eval_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts);

// If the numerically-optimal approximation to a function
// of two variables has rank $n$, then any integral
// (or sum, here) over one of the two variables can
// be fully described by $n$ coefficients.
void num_opt_rankn_integral_coeffs(
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

// Two OCaml stubs because there are more than
// 5 arguments for num_opt_rankn_integral_coeffs' OCaml
// version.
CAMLprim value caml_num_opt_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn);

CAMLprim value caml_num_opt_rankn_integral_coeffs_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts);

#endif // KERNLCAML_STUBS_H
