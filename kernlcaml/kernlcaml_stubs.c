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

#include "kernlcaml_stubs.h"
#include <malloc.h>

// Brent zero-finder:

// First create a static value and evaluation routine
// for the primitive that will operate on an OCaml
// closure.

static value zeroin_closure;

double eval_zeroin_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(zeroin_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

// Now define the primitive itself.

CAMLprim value caml_zeroin(
	value vf,
	value va,
	value vb,
	value vtol)
{
	CAMLparam4(vf, va, vb, vtol);
	CAMLlocal1(vzero);

	// Set the static closure 
	zeroin_closure = vf;
	caml_register_global_root(&zeroin_closure);
	// Get raw doubles from the input params
	double a = Double_val(va);
	double b = Double_val(vb);
	double tol = Double_val(vtol);

	double zerof = zeroin(
		&a,
		&b,
		&eval_zeroin_closure,
		&tol);

	vzero = caml_copy_double(zerof);

	caml_remove_global_root(&zeroin_closure);
	CAMLreturn(vzero);
}

// brent_zero

static value brent_zero_closure;

double eval_brent_zero_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(brent_zero_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

CAMLprim value caml_brent_zero(
	value vf,
	value va,
	value vb,
	value vtol,
	value vmax_iter)
{
	CAMLparam5(vf, va, vb, vtol, vmax_iter);
	CAMLlocal5(vxzero, vfzero, vniter, verr_code, vres);

	// Set the static closure 
	brent_zero_closure = vf;
	caml_register_global_root(&brent_zero_closure);

	// Get raw doubles from the input params
	double a = Double_val(va);
	double b = Double_val(vb);
	double tol = Double_val(vtol);
	int max_iter = Int_val(vmax_iter);

	int niter;
	int err_code;
	double xzero;
	double fzero;

	brent_zero(
		&a,
		&b,
		&eval_brent_zero_closure,
		&tol,
		&max_iter,
		&niter,
		&err_code,
		&xzero,
		&fzero);

	// Copy results into OCaml values
	vxzero = caml_copy_double(xzero);
	vfzero = caml_copy_double(fzero);
	vniter = Val_int(niter);
	verr_code = Val_int(err_code);

	// Allocate an OCaml tuple
	vres = caml_alloc(4, 0);

	// Fill the tuple
	Store_field(vres, 0, vxzero);
	Store_field(vres, 1, vfzero);
	Store_field(vres, 2, vniter);
	Store_field(vres, 3, verr_code);

	caml_remove_global_root(&brent_zero_closure);

	CAMLreturn(vres);

}

// f_min functionality

static value f_min_closure;

double eval_f_min_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(f_min_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

CAMLprim value caml_f_min(
	value vf,
	value va,
	value vb,
	value vtol) 
{

	CAMLparam4(vf, va, vb, vtol);
	CAMLlocal4(vxopt, vfopt, verr_code, vres);

	// Set the static closure 
	f_min_closure = vf;
	caml_register_global_root(&f_min_closure);

	// Get raw doubles from the input params
	double a = Double_val(va);
	double b = Double_val(vb);
	double tol = Double_val(vtol);
	
	double xopt;
	double fopt;
	int err_code = 0;
	char err_msg[256];
	err_msg[0] = '\0';

	f_min(
		&a,
		&b,
		&eval_f_min_closure,
		&tol,
		&xopt,
		&fopt,
		&err_code,
		err_msg);

	// Copy results to return into OCaml values
	vxopt = caml_copy_double(xopt);
	vfopt = caml_copy_double(fopt);
	verr_code = Val_int(err_code);

	// Allocate an OCaml tuple
	vres = caml_alloc(3, 0);

	// Fill the tuple
	Store_field(vres, 0, vxopt);
	Store_field(vres, 1, vfopt);
	Store_field(vres, 2, verr_code);

	caml_remove_global_root(&f_min_closure);

	CAMLreturn(vres);

}

// Brent minimization (like f_min, but with control over max_iter)
static value brent_min_closure;

double eval_brent_min_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(brent_min_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

CAMLprim value caml_brent_min(
	value vf,
	value va,
	value vb,
	value vtol,
	value vmax_iter)
{

	CAMLparam5(vf, va, vb, vtol, vmax_iter);
	CAMLlocal5(vxopt, vfopt, vniter, verr_code, vres);

	// Set the static closure 
	brent_min_closure = vf;
	caml_register_global_root(&brent_min_closure);

	// Get raw doubles from the input params
	double a = Double_val(va);
	double b = Double_val(vb);
	double tol = Double_val(vtol);
	int max_iter = Int_val(vmax_iter);

	double xopt;
	double fopt;
	int niter;
	int err_code = 0;
	
	brent_min(
		&a,
		&b,
		&eval_brent_min_closure,
		&tol,
		&max_iter,
		&niter,
		&err_code,
		&xopt,
		&fopt);

	// Copy results to return into OCaml values
	vxopt = caml_copy_double(xopt);
	vfopt = caml_copy_double(fopt);
	vniter = Val_int(niter);
	verr_code = Val_int(err_code);

	// Allocate an OCaml tuple
	vres = caml_alloc(4, 0);

	// Fill the tuple
	Store_field(vres, 0, vxopt);
	Store_field(vres, 1, vfopt);
	Store_field(vres, 2, vniter);
	Store_field(vres, 3, verr_code);

	caml_remove_global_root(&brent_min_closure);

	CAMLreturn(vres);

}

static value find_all_zeros_closure;

double eval_findallz_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(find_all_zeros_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

CAMLprim value caml_find_all_zeros(
	value vf,
	value vgrid,
	value vtol)
{
	CAMLparam3(vf, vgrid, vtol);
	CAMLlocal1(vzeros);

	// Set the static closure 
	find_all_zeros_closure = vf;

	caml_register_global_root(&find_all_zeros_closure);

	double *grid_pts = Caml_ba_data_val(vgrid);
	int n_grid_pts = Caml_ba_array_val(vgrid)->dim[0];

	double toler = Double_val(vtol);

	int max_num_zeros = n_grid_pts;
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
		&eval_findallz_closure,
		&n_grid_pts,
		grid_pts,
		&toler,
		&c_err_stat,
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = num_zeros;
	vzeros = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);
	// Now copy from result_zeros -- this avoids the need to
	// manage memory on the C side, as I can free result_zeros
	// after the copy
	double *temp_ba_data = Caml_ba_data_val(vzeros);
	int i;
	for (i = 0; i < num_zeros; ++i) {
		temp_ba_data[i] = result_zeros[i];
	}

	free(result_zeros);

	caml_remove_global_root(&find_all_zeros_closure);

	CAMLreturn(vzeros);
}

// Find the set of relative optima on the intervals defined
// by the given grid
static value find_rel_optima_closure;

double eval_find_rel_optima_closure(double *x) {
	CAMLparam0();
	CAMLlocal2(vx, vy);
	vx = caml_copy_double(*x);
	vy = caml_callback(find_rel_optima_closure, vx);
	CAMLreturnT(double, Double_val(vy));
}

CAMLprim value caml_find_rel_optima(
	value vf,
	value vgrid,
	value vtol)
{
	CAMLparam3(vf, vgrid, vtol);
	CAMLlocal1(voptima);

	// Set the static closure 
	find_rel_optima_closure = vf;
	caml_register_global_root(&find_rel_optima_closure);

	double *grid_pts = Caml_ba_data_val(vgrid);
	int n_grid_pts = Caml_ba_array_val(vgrid)->dim[0];

	double toler = Double_val(vtol);

	int num_optima = n_grid_pts - 1;
	double *result_optima;
	result_optima = (double *)malloc(sizeof(double) * 2 * num_optima);
	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	find_rel_optima(
		result_optima,
		&eval_find_rel_optima_closure,
		&n_grid_pts,
		grid_pts,
		&toler,
		&c_err_stat,
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[2];
	ba_dim[0] = num_optima;
	ba_dim[1] = 2;
	voptima = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);
	// Now copy from result_optima -- this avoids the need to
	// manage memory on the C side, as I can free result_optima
	// after the copy
	double *temp_ba_data = Caml_ba_data_val(voptima);
	int i;
	for (i = 0; i < num_optima; ++i) {
		temp_ba_data[i] = result_optima[i];
		temp_ba_data[num_optima + i] = result_optima[num_optima + i];
	}

	free(result_optima);

	caml_remove_global_root(&find_rel_optima_closure);

	CAMLreturn(voptima);
}

// Get n Chebyshev points
CAMLprim value caml_cheb_pts(
	value vn)
{

	CAMLparam1(vn);
	CAMLlocal1(vcheb);

	int n = Int_val(vn);

	double *result_pts;
	result_pts = (double *)malloc(sizeof(double) * n);

	chebyshev_points(
		&n,
		result_pts);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n;
	vcheb = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);
	// Now copy from result_pts -- this avoids the need to
	// manage memory on the C side, as I can free result_pts
	// after the copy
	double *temp_ba_data = Caml_ba_data_val(vcheb);
	int i;
	for (i = 0; i < n; ++i) {
		temp_ba_data[i] = result_pts[i];
	}

	free(result_pts);

	CAMLreturn(vcheb);

}

// Matrix multiplication, provided by Fortran
CAMLprim value caml_matmul(
	value va,
	value vb)
{
	CAMLparam2(va, vb);
	CAMLlocal1(vc);

	double *amat = Caml_ba_data_val(va);
	int m = Caml_ba_array_val(va)->dim[0];
	int n = Caml_ba_array_val(va)->dim[1];

	double *bmat = Caml_ba_data_val(vb);
	int p = Caml_ba_array_val(vb)->dim[1];

	double *cmat;
	cmat = (double *)malloc(sizeof(double) * m * p);

	matmul(
		&m,
		&n,
		&p,
		amat,
		bmat,
		cmat);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[2];
	ba_dim[0] = m;
	ba_dim[1] = p;
	vc = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);
	
	double *temp_ba_data = Caml_ba_data_val(vc);
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < p; ++j) {
			temp_ba_data[m*j + i] = cmat[m*j + i];
		}
	}

	free(cmat);

	CAMLreturn(vc);
}

// Solution of a general linear system, via LAPACK
CAMLprim value caml_linear_solve(
	value va,
	value vb)
{
	CAMLparam2(va, vb);
	CAMLlocal1(vx);

	double *amat = Caml_ba_data_val(va);
	int m = Caml_ba_array_val(va)->dim[0];
	
	double *bmat = Caml_ba_data_val(vb);
	int p = Caml_ba_array_val(vb)->dim[1];

	double *xmat;
	xmat = (double *)malloc(sizeof(double) * m * p);

	int info;

	linear_solve(
		&m,
		amat,
		&p,
		bmat,
		xmat,
		&info);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[2];
	ba_dim[0] = m;
	ba_dim[1] = p;
	vx = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vx);
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < p; ++j) {
			temp_ba_data[m*j + i] = xmat[m*j + i];
		}
	}

	free(xmat);

	CAMLreturn(vx);
}

// Moore-Penrose inverse
CAMLprim value caml_mp_inverse(
	value va,
	value vtol)
{
	CAMLparam2(va, vtol);
	CAMLlocal1(vx);

	double *amat = Caml_ba_data_val(va);
	int m = Caml_ba_array_val(va)->dim[0];
	int n = Caml_ba_array_val(va)->dim[1];

	double toler = Double_val(vtol);

	double *xmat;
	xmat = (double *)malloc(sizeof(double) * m * n);

	int info;

	mp_inverse(
		xmat,
		amat,
		m,
		n,
		toler);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[2];
	ba_dim[0] = m;
	ba_dim[1] = n;
	vx = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vx);
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			temp_ba_data[m*j + i] = xmat[m*j + i];
		}
	}

	free(xmat);

	CAMLreturn(vx);
}

// Black (1976) option pricing formula
CAMLprim value caml_black_formula(
	value vcall,
	value vstrike,
	value vforward,
	value vvol,
	value vdiscfac)
{
	CAMLparam5(vcall, vstrike, vforward, vvol, vdiscfac);
	CAMLlocal1(vres);

	bool is_call = Bool_val(vcall);
	double c_strike = Double_val(vstrike);
	double c_forward = Double_val(vforward);
	double c_vol = Double_val(vvol);
	double c_disc_fac = Double_val(vdiscfac);

	double price = black_formula(
		&is_call,
		&c_strike,
		&c_forward,
		&c_vol,
		&c_disc_fac);

	vres = caml_copy_double(price);

	CAMLreturn(vres);
}

// Bachelier (or normal-distribution) option pricing formula
CAMLprim value caml_normopt_formula(
	value vcall,
	value vstrike,
	value vforward,
	value vvol,
	value vdiscfac)
{
	CAMLparam5(vcall, vstrike, vforward, vvol, vdiscfac);
	CAMLlocal1(vres);

	bool is_call = Bool_val(vcall);
	double c_strike = Double_val(vstrike);
	double c_forward = Double_val(vforward);
	double c_vol = Double_val(vvol);
	double c_disc_fac = Double_val(vdiscfac);

	double price = normopt_formula(
		&is_call,
		&c_strike,
		&c_forward,
		&c_vol,
		&c_disc_fac);

	vres = caml_copy_double(price);

	CAMLreturn(vres);
}

// Closure and evaluator for Borsuk lower bound
static value borsuk_closure;

double eval_borsuk_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(borsuk_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Borsuk lower bound, bytecode version
CAMLprim value caml_borsuk_lower_bound_bytecode(
	value* argv,
	int argn) 
{
	return caml_borsuk_lower_bound_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5]);
}

// Borsuk lower bound, native-code version
CAMLprim value caml_borsuk_lower_bound_native(
	value vf,
	value vrho_vec,
	value voverx,
	value vtol,
	value vmax_iter,
	value vgrid)
{
	CAMLparam5(vf, vrho_vec, voverx, vtol, vmax_iter);
	CAMLxparam1(vgrid);
	CAMLlocal5(varray, vdiscrep, vniter, verrstat, vres);

	// Set the static closure 
	borsuk_closure = vf;
	caml_register_global_root(&borsuk_closure);

	double *grid_pts = Caml_ba_data_val(vgrid);
	int n_grid_pts = Caml_ba_array_val(vgrid)->dim[0];

	double *rho_vec = Caml_ba_data_val(vrho_vec);
	int size_of_rho = Caml_ba_array_val(vrho_vec)->dim[0];

	bool overx = Bool_val(voverx);
	double toler = Double_val(vtol);
	int max_iter = Int_val(vmax_iter);

	double* coeffs = (double *)malloc(sizeof(double) * size_of_rho);
	double* nodes = (double *)malloc(sizeof(double) * size_of_rho);
	double* errors_at_nodes = (double *)malloc(sizeof(double) * size_of_rho);

	double discrep;
	int niter;
	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	borsuk_lower_bound(
		&eval_borsuk_closure,
		&size_of_rho,
		rho_vec,
		&overx,
		&toler,
		&max_iter,
		&n_grid_pts,
		grid_pts,
		coeffs,
		nodes,
		errors_at_nodes,
		&discrep,
		&niter,
		&c_err_stat,
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[2];
	ba_dim[0] = size_of_rho;
	ba_dim[1] = 3; // Coeffs in first column, then nodes, then errors at nodes
	varray = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);
	
	double *temp_ba_data = Caml_ba_data_val(varray);
	int i;
	for (i = 0; i < size_of_rho; ++i) {
		temp_ba_data[i] = coeffs[i];
		temp_ba_data[size_of_rho + i] = nodes[i];
		temp_ba_data[2 * size_of_rho + i] = errors_at_nodes[i];
	}

	free(coeffs);
	free(nodes);
	free(errors_at_nodes);

	// Copy the other outputs to OCaml values
	vdiscrep = caml_copy_double(discrep);
	vniter = Val_int(niter);
	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(4, 0);

	// Fill the tuple
	Store_field(vres, 0, varray);
	Store_field(vres, 1, vdiscrep);
	Store_field(vres, 2, vniter);
	Store_field(vres, 3, verrstat);

	caml_remove_global_root(&borsuk_closure);

	CAMLreturn(vres);

}

// Closure and evaluator for worst rho and gamma
static value worst_rho_gamma_closure;

double eval_worst_rho_gamma_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(worst_rho_gamma_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Worst rho and gamma, bytecode version
CAMLprim value caml_worst_rho_and_gamma_bytecode(
	value* argv,
	int argn)
{
	return caml_worst_rho_and_gamma_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7]);
}

// Worst rho and gamma, native code version
CAMLprim value caml_worst_rho_and_gamma_native(
	value vf,
	value vx_lb,
	value vx_ub,
	value vy_lb,
	value vy_ub,
	value vnumterms,
	value vtol,
	value vmax_iter)
{
	CAMLparam5(vf, vx_lb, vx_ub, vy_lb, vy_ub);
	CAMLxparam3(vnumterms, vtol, vmax_iter);
	CAMLlocal5(vrho_res, vgamma_res, vamat_rho, vamat_gamma, vres);
	CAMLlocal5(vdiscrep, vniter, vkerdiscrep, vkerniter, verrstat);

	// Set the static closure 
	worst_rho_gamma_closure = vf;
	caml_register_global_root(&worst_rho_gamma_closure);

	double x_lb = Double_val(vx_lb);
	double x_ub = Double_val(vx_ub);
	double y_lb = Double_val(vy_lb);
	double y_ub = Double_val(vy_ub);
	int num_terms = Int_val(vnumterms);
	double toler = Double_val(vtol);
	int max_iter = Int_val(vmax_iter);

	double *rho_vec = (double *)malloc(sizeof(double) * num_terms);
	double *gamma_vec = (double *)malloc(sizeof(double) * num_terms);

	double *coeffs_rho = (double *)malloc(sizeof(double) * num_terms);
	double *coeffs_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *nodes_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *errors_at_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *errors_at_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);

	double* a_matrix_rho = (double *)malloc(
		sizeof(double) * num_terms * num_terms);
	double* a_matrix_gamma = (double *)malloc(
		sizeof(double) * num_terms * num_terms);

	double *b_coeffs_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_coeffs_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *b_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);
	double *b_errors_at_nodes_rho = (double *)malloc(sizeof(double) * num_terms);
	double *b_errors_at_nodes_gamma = (double *)malloc(sizeof(double) * num_terms);

	double discrepancy;
	int num_iter;

	double ker_discrep;
	int ker_num_iter;

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	worst_rho_and_gamma(
		&eval_worst_rho_gamma_closure,
		&x_lb,
		&x_ub,
		&y_lb,
		&y_ub,
		&num_terms,
		&toler,
		&max_iter,
		rho_vec,
		gamma_vec,
		coeffs_rho,
		coeffs_gamma,
		a_matrix_rho,
		a_matrix_gamma,
		nodes_rho,
		nodes_gamma,
		errors_at_nodes_rho,
		errors_at_nodes_gamma,
		b_coeffs_rho,
		b_coeffs_gamma,
		b_nodes_rho,
		b_nodes_gamma,
		b_errors_at_nodes_rho,
		b_errors_at_nodes_gamma,
		&discrepancy,
		&num_iter,
		&ker_discrep,
		&ker_num_iter,
		&c_err_stat,
		c_err_msg);

	/**********************  Rho Results       **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat rho_ba_dim[2];
	rho_ba_dim[0] = num_terms;
	// Rho in first column, then coeffs, then nodes, then errors at nodes, 
	// then Borsuk coeffs, then Borsuk nodes, then Borsuk errors at nodes
	rho_ba_dim[1] = 7; 
	vrho_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		rho_ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vrho_res);
	int i;
	for (i = 0; i < num_terms; ++i) {
		temp_ba_data[i] = rho_vec[i];
		temp_ba_data[num_terms + i] = coeffs_rho[i];
		temp_ba_data[2 * num_terms + i] = nodes_rho[i];
		temp_ba_data[3 * num_terms + i] = errors_at_nodes_rho[i];
		temp_ba_data[4 * num_terms + i] = b_coeffs_rho[i];
		temp_ba_data[5 * num_terms + i] = b_nodes_rho[i];
		temp_ba_data[6 * num_terms + i] = b_errors_at_nodes_rho[i];
	}

	/**********************  Gamma Results     **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat gamma_ba_dim[2];
	gamma_ba_dim[0] = num_terms;
	// Gamma in first column, then coeffs, then nodes, then errors at nodes, 
	// then Borsuk coeffs, then Borsuk nodes, then Borsuk errors at nodes
	gamma_ba_dim[1] = 7;
	vgamma_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		gamma_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vgamma_res);
	for (i = 0; i < num_terms; ++i) {
		temp_ba_data[i] = gamma_vec[i];
		temp_ba_data[num_terms + i] = coeffs_gamma[i];
		temp_ba_data[2 * num_terms + i] = nodes_gamma[i];
		temp_ba_data[3 * num_terms + i] = errors_at_nodes_gamma[i];
		temp_ba_data[4 * num_terms + i] = b_coeffs_gamma[i];
		temp_ba_data[5 * num_terms + i] = b_nodes_gamma[i];
		temp_ba_data[6 * num_terms + i] = b_errors_at_nodes_gamma[i];
	}

	/**********************  Amat Rho Results  **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat amat_rho_ba_dim[2];
	amat_rho_ba_dim[0] = num_terms;
	amat_rho_ba_dim[1] = num_terms;
	vamat_rho = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		amat_rho_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vamat_rho);
	int j;
	for (i = 0; i < num_terms; ++i) {
		for (j = 0; j < num_terms; ++j) {
			temp_ba_data[num_terms * j + i] = a_matrix_rho[num_terms * j + i];
		}
	}

	/********************** Amat Gamma Results **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat amat_gamma_ba_dim[2];
	amat_gamma_ba_dim[0] = num_terms;
	amat_gamma_ba_dim[1] = num_terms;
	vamat_gamma = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		amat_gamma_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vamat_gamma);
	for (i = 0; i < num_terms; ++i) {
		for (j = 0; j < num_terms; ++j) {
			temp_ba_data[num_terms * j + i] = a_matrix_gamma[num_terms * j + i];
		}
	}

	free(rho_vec);
	free(gamma_vec);
	free(coeffs_rho);
	free(coeffs_gamma);
	free(nodes_rho);
	free(nodes_gamma);
	free(errors_at_nodes_rho);
	free(errors_at_nodes_gamma);
	free(a_matrix_rho);
	free(a_matrix_gamma);
	free(b_coeffs_rho);
	free(b_coeffs_gamma);
	free(b_nodes_rho);
	free(b_nodes_gamma);
	free(b_errors_at_nodes_rho);
	free(b_errors_at_nodes_gamma);

	// Copy the other outputs to OCaml values
	vdiscrep = caml_copy_double(discrepancy);
	vkerdiscrep = caml_copy_double(ker_discrep);
	vniter = Val_int(num_iter);
	vkerniter = Val_int(ker_num_iter);
	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(9, 0);

	// Fill the tuple
	Store_field(vres, 0, vrho_res);
	Store_field(vres, 1, vgamma_res);
	Store_field(vres, 2, vdiscrep);
	Store_field(vres, 3, vniter);
	Store_field(vres, 4, vkerdiscrep);
	Store_field(vres, 5, vkerniter);
	Store_field(vres, 6, verrstat);
	Store_field(vres, 7, vamat_rho);
	Store_field(vres, 8, vamat_gamma);

	caml_remove_global_root(&worst_rho_gamma_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for kernel_integral_eval
static value kernel_integral_eval_closure;

double eval_kernel_integral_eval_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(kernel_integral_eval_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_kernel_integral_eval(
	value vf,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts)
{

	CAMLparam5(vf, voverx, vintegral_pts, vintegral_wts, veval_pts);
	CAMLlocal1(vres);

	// Set the static closure 
	kernel_integral_eval_closure = vf;
	caml_register_global_root(&kernel_integral_eval_closure);

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double *eval_pts = Caml_ba_data_val(veval_pts);
	int n_eval_pts = Caml_ba_array_val(veval_pts)->dim[0];

	bool overx = Bool_val(voverx);

	double* results = (double *)malloc(sizeof(double) * n_eval_pts);

	kernel_integral_eval(
		&eval_kernel_integral_eval_closure, 
		&overx, 
		&n_integral_pts, 
		integral_pts, 
		integral_wts, 
		&n_eval_pts, 
		eval_pts, 
		results);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_pts;
	vres = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vres);
	int i;
	for (i = 0; i < n_eval_pts; ++i) {
		temp_ba_data[i] = results[i];
	}

	free(results);

	caml_remove_global_root(&kernel_integral_eval_closure);

	CAMLreturn(vres);

}

// Closure and evaluator for taylor_rankn_params
static value taylor_rankn_params_closure;

double eval_taylor_rankn_params_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(taylor_rankn_params_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_taylor_rankn_params(
	value vf,
	value vrank,
	value vx_center,
	value vy_center,
	value vfin_diff_step)
{
	CAMLparam5(vf, vrank, vx_center, vy_center, vfin_diff_step);
	CAMLlocal4(v_vmat, v_wmat, verrstat, vres);

	// Set the static closure 
	taylor_rankn_params_closure = vf;
	caml_register_global_root(&taylor_rankn_params_closure);

	int rank = Int_val(vrank);
	double xcenter = Double_val(vx_center);
	double ycenter = Double_val(vy_center);
	double fdiff_step = Double_val(vfin_diff_step);

	double* vmat = (double *)malloc(sizeof(double) * rank * rank);
	double* wmat = (double *)malloc(sizeof(double) * rank * rank);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	taylor_rankn_params(
		&eval_taylor_rankn_params_closure, 
		&rank, 
		&xcenter, 
		&ycenter, 
		vmat, 
		wmat, 
		&c_err_stat, 
		c_err_msg,
		&fdiff_step);

	/********************* V matrix              ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat v_ba_dim[2];
	v_ba_dim[0] = rank;
	v_ba_dim[1] = rank;
	v_vmat = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		v_ba_dim);

	double *temp_ba_data = Caml_ba_data_val(v_vmat);
	int i, j;
	for (i = 0; i < rank; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[rank * j + i] = vmat[rank * j + i];
		}
	}

	/********************* W matrix              ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat w_ba_dim[2];
	w_ba_dim[0] = rank;
	w_ba_dim[1] = rank;
	v_wmat = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		w_ba_dim);

	// Reuse the pointer from above
	temp_ba_data = Caml_ba_data_val(v_wmat);
	for (i = 0; i < rank; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[rank * j + i] = wmat[rank * j + i];
		}
	}

	free(vmat);
	free(wmat);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(3, 0);

	// Fill the tuple
	Store_field(vres, 0, v_vmat);
	Store_field(vres, 1, v_wmat);
	Store_field(vres, 2, verrstat);

	caml_remove_global_root(&taylor_rankn_params_closure);

	CAMLreturn(vres);

}

// Taylor rank-$n$ eval, bytecode version
CAMLprim value caml_taylor_rankn_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_taylor_rankn_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6]);
}

// Taylor rank-$n$ eval, native code version
CAMLprim value caml_taylor_rankn_eval_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value veval_x_pts,
	value veval_y_pts)
{
	CAMLparam5(vrank, vx_center, vy_center, v_vmat, v_wmat);
	CAMLxparam2(veval_x_pts, veval_y_pts);
	CAMLlocal3(v_eval_res, verrstat, vres);

	int rank = Int_val(vrank);
	double xcenter = Double_val(vx_center);
	double ycenter = Double_val(vy_center);

	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *eval_x_pts = Caml_ba_data_val(veval_x_pts);
	int n_eval_x_pts = Caml_ba_array_val(veval_x_pts)->dim[0];

	double *eval_y_pts = Caml_ba_data_val(veval_y_pts);
	int n_eval_y_pts = Caml_ba_array_val(veval_y_pts)->dim[0];

	double* results = (double*)malloc(sizeof(double) * n_eval_x_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	taylor_rankn_eval(
		&rank, 
		&xcenter, 
		&ycenter, 
		vmat, 
		wmat, 
		&n_eval_x_pts, 
		eval_x_pts,
		eval_y_pts, 
		results, 
		&c_err_stat, 
		c_err_msg);
		
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_x_pts;
	v_eval_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(v_eval_res);
	int i;
	for (i = 0; i < n_eval_x_pts; ++i) {
		temp_ba_data[i] = results[i];
	}

	free(results);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, v_eval_res);
	Store_field(vres, 1, verrstat);

	CAMLreturn(vres);
}

// Evaluation of a rank-$n$ integral of a Taylor
// series approximation, bytecode version
CAMLprim value caml_taylor_rankn_integral_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_taylor_rankn_integral_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7],
		argv[8]);
}

// Evaluation of a rank-$n$ integral of a Taylor
// series approximation, native code version
CAMLprim value caml_taylor_rankn_integral_eval_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts)
{
	CAMLparam5(vrank, vx_center, vy_center, v_vmat, v_wmat);
	CAMLxparam4(voverx, vintegral_pts, vintegral_wts, veval_pts);
	CAMLlocal3(v_eval_res, verrstat, vres);

	int rank = Int_val(vrank);
	double xcenter = Double_val(vx_center);
	double ycenter = Double_val(vy_center);
	bool overx = Bool_val(voverx);

	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double *eval_pts = Caml_ba_data_val(veval_pts);
	int n_eval_pts = Caml_ba_array_val(veval_pts)->dim[0];

	double* results = (double*)malloc(sizeof(double) * n_eval_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	taylor_rankn_integral_eval(
		&rank,
		&xcenter,
		&ycenter,
		vmat,
		wmat,
		&overx,
		&n_integral_pts,
		integral_pts,
		integral_wts,
		&n_eval_pts,
		eval_pts,
		results,
		&c_err_stat,
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_pts;
	v_eval_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(v_eval_res);
	int i;
	for (i = 0; i < n_eval_pts; ++i) {
		temp_ba_data[i] = results[i];
	}

	free(results);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, v_eval_res);
	Store_field(vres, 1, verrstat);

	CAMLreturn(vres);
}

// Coefficients of a rank-$n$ integral of a Taylor
// series approximation, bytecode version
CAMLprim value caml_taylor_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn)
{
	return caml_taylor_rankn_integral_coeffs_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7]);
}

// Coefficients of a rank-$n$ integral of a Taylor
// series approximation, native code version
CAMLprim value caml_taylor_rankn_integral_coeffs_native(
	value vrank,
	value vx_center,
	value vy_center,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts)
{
	CAMLparam5(vrank, vx_center, vy_center, v_vmat, v_wmat);
	CAMLxparam3(voverx, vintegral_pts, vintegral_wts);
	CAMLlocal3(v_coeff_res, verrstat, vres);

	int rank = Int_val(vrank);
	double xcenter = Double_val(vx_center);
	double ycenter = Double_val(vy_center);
	bool overx = Bool_val(voverx);

	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double* coeffs = (double*)malloc(sizeof(double) * rank);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	taylor_rankn_integral_coeffs(
		&rank,
		&xcenter,
		&ycenter,
		vmat,
		wmat,
		&overx,
		&n_integral_pts,
		integral_pts,
		integral_wts,
		&rank,
		coeffs,
		&c_err_stat,
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = rank;
	v_coeff_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(v_coeff_res);
	int i;
	for (i = 0; i < rank; ++i) {
		temp_ba_data[i] = coeffs[i];
	}

	free(coeffs);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, v_coeff_res);
	Store_field(vres, 1, verrstat);

	CAMLreturn(vres);
}

// Closure and evaluator for eigen_rankn_params
static value eigen_rankn_params_closure;

double eval_eigen_rankn_params_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(eigen_rankn_params_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Parameters of a singular function series
// approximation, bytecode version
CAMLprim value caml_eigen_rankn_params_bytecode(
	value* argv,
	int argn)
{
	return caml_eigen_rankn_params_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5]);
}

// Parameters of a singular function series
// approximation, native code version
CAMLprim value caml_eigen_rankn_params_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_wts,
	value vy_wts)
{
	CAMLparam5(vf, vrank, vx_pts, vy_pts, vx_wts);
	CAMLxparam1(vy_wts);
	CAMLlocal5(vx_funcwts, vy_funcwts, vtrunc_svals, verrstat, vres);

	// Set the static closure 
	eigen_rankn_params_closure = vf;
	caml_register_global_root(&eigen_rankn_params_closure);

	int rank = Int_val(vrank);

	double *x_pts = Caml_ba_data_val(vx_pts);
	int n_x_pts = Caml_ba_array_val(vx_pts)->dim[0];

	double *y_pts = Caml_ba_data_val(vy_pts);
	int n_y_pts = Caml_ba_array_val(vy_pts)->dim[0];

	double *x_wts = Caml_ba_data_val(vx_wts);
	int n_x_wts = Caml_ba_array_val(vx_wts)->dim[0];

	double *y_wts = Caml_ba_data_val(vy_wts);
	int n_y_wts = Caml_ba_array_val(vy_wts)->dim[0];

	// The sizing below is not a typo; in the Nystrom approach,
	// the number of points in the X function weighting is equal to
	// the number of original Y points, and vice versa.
	double* x_funcwts = (double *)malloc(sizeof(double) * n_y_pts * rank);
	double* y_funcwts = (double *)malloc(sizeof(double) * n_x_pts * rank);
	double* trunc_svals = (double *)malloc(sizeof(double) * rank);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	eigen_rankn_params(
		&eval_eigen_rankn_params_closure, 
		&rank, 
		&n_x_pts, 
		&n_y_pts, 
		x_pts, 
		y_pts, 
		x_wts, 
		y_wts, 
		x_funcwts, 
		y_funcwts, 
		trunc_svals, 
		&c_err_stat, 
		c_err_msg);

	/********************* X Funcwts             ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat x_ba_dim[2];
	x_ba_dim[0] = n_y_pts;
	x_ba_dim[1] = rank;
	vx_funcwts = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		x_ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vx_funcwts);
	int i, j;
	for (i = 0; i < n_y_pts; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[n_y_pts * j + i] = x_funcwts[n_y_pts * j + i];
		}
	}

	/********************* Y Funcwts             ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat y_ba_dim[2];
	y_ba_dim[0] = n_x_pts;
	y_ba_dim[1] = rank;
	vy_funcwts = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		y_ba_dim);

	// Reuse the pointer from above
	temp_ba_data = Caml_ba_data_val(vy_funcwts);
	for (i = 0; i < n_x_pts; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[n_x_pts * j + i] = y_funcwts[n_x_pts * j + i];
		}
	}

	/********************* Truncated Svals       ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat t_ba_dim[1];
	t_ba_dim[0] = rank;
	vtrunc_svals = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		t_ba_dim);

	// Reuse the pointer from above
	temp_ba_data = Caml_ba_data_val(vtrunc_svals);
	for (i = 0; i < rank; ++i) {
		temp_ba_data[i] = trunc_svals[i];
	}

	free(x_funcwts);
	free(y_funcwts);
	free(trunc_svals);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(4, 0);

	// Fill the tuple
	Store_field(vres, 0, vx_funcwts);
	Store_field(vres, 1, vy_funcwts);
	Store_field(vres, 2, vtrunc_svals);
	Store_field(vres, 3, verrstat);

	caml_remove_global_root(&eigen_rankn_params_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for eigen_rankn_eval
static value eigen_rankn_eval_closure;

double eval_eigen_rankn_eval_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(eigen_rankn_eval_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Evaluation of a rank-$n$ singular function approximation,
// bytecode version
CAMLprim value caml_eigen_rankn_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_eigen_rankn_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7],
		argv[8]);
}

// Evaluation of a rank-$n$ singular function approximation,
// native code version
CAMLprim value caml_eigen_rankn_eval_native(
	value vf,
	value vrank,
	value vx_pts,
	value vy_pts,
	value vx_funcwts,
	value vy_funcwts,
	value vtrunc_svals,
	value veval_x_pts,
	value veval_y_pts)
{
	CAMLparam5(vf, vrank, vx_pts, vy_pts, vx_funcwts);
	CAMLxparam4(vy_funcwts, vtrunc_svals, veval_x_pts, veval_y_pts);
	CAMLlocal3(vevalres, verrstat, vres);

	// Set the static closure 
	eigen_rankn_eval_closure = vf;
	caml_register_global_root(&eigen_rankn_eval_closure);

	int rank = Int_val(vrank);

	double *x_pts = Caml_ba_data_val(vx_pts);
	int n_x_pts = Caml_ba_array_val(vx_pts)->dim[0];

	double *y_pts = Caml_ba_data_val(vy_pts);
	int n_y_pts = Caml_ba_array_val(vy_pts)->dim[0];

	double *x_funcwts = Caml_ba_data_val(vx_funcwts);
	int n_rows_x = Caml_ba_array_val(vx_funcwts)->dim[0];
	int n_cols_x = Caml_ba_array_val(vx_funcwts)->dim[1];

	double *y_funcwts = Caml_ba_data_val(vy_funcwts);
	int n_rows_y = Caml_ba_array_val(vy_funcwts)->dim[0];
	int n_cols_y = Caml_ba_array_val(vy_funcwts)->dim[1];

	double *trunc_svals = Caml_ba_data_val(vtrunc_svals);
	int n_svals = Caml_ba_array_val(vtrunc_svals)->dim[0];

	double *eval_x_pts = Caml_ba_data_val(veval_x_pts);
	int n_eval_x_pts = Caml_ba_array_val(veval_x_pts)->dim[0];

	double *eval_y_pts = Caml_ba_data_val(veval_y_pts);
	int n_eval_y_pts = Caml_ba_array_val(veval_y_pts)->dim[0];

	double* evalres = (double *)malloc(sizeof(double) * n_eval_x_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	eigen_rankn_eval(
		&eval_eigen_rankn_eval_closure,
		&rank,
		&n_x_pts,
		&n_y_pts,
		x_pts,
		y_pts,
		x_funcwts,
		y_funcwts,
		trunc_svals,
		&n_eval_x_pts,
		eval_x_pts,
		eval_y_pts,
		evalres,
		&c_err_stat,
		c_err_msg);

	/********************* Eval result           ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_x_pts;
	vevalres = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vevalres);
	int i, j;
	for (i = 0; i < n_eval_x_pts; ++i) {
		temp_ba_data[i] = evalres[i];
	}

	free(evalres);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, vevalres);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&eigen_rankn_eval_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for eigen_rankn_integral_eval
static value eigen_rankn_integral_eval_closure;

double eval_eigen_rankn_integral_eval_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(eigen_rankn_integral_eval_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Evaluation of an integral (really, a weighted sum) of
// a rank-$n$ singular function approximation, bytecode version
CAMLprim value caml_eigen_rankn_integral_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_eigen_rankn_integral_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7],
		argv[8],
		argv[9],
		argv[10]);
}

// Evaluation of an integral (really, a weighted sum) of
// a rank-$n$ singular function approximation, native code version
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
	value veval_pts)
{
	CAMLparam5(vf, vrank, vx_pts, vy_pts, vx_funcwts);
	CAMLxparam5(vy_funcwts, vtrunc_svals, voverx, vintegral_pts, vintegral_wts);
	CAMLxparam1(veval_pts);
	CAMLlocal3(vevalres, verrstat, vres);

	// Set the static closure 
	eigen_rankn_integral_eval_closure = vf;
	caml_register_global_root(&eigen_rankn_integral_eval_closure);

	int rank = Int_val(vrank);
	bool overx = Bool_val(voverx);

	double *x_pts = Caml_ba_data_val(vx_pts);
	int n_x_pts = Caml_ba_array_val(vx_pts)->dim[0];

	double *y_pts = Caml_ba_data_val(vy_pts);
	int n_y_pts = Caml_ba_array_val(vy_pts)->dim[0];

	double *x_funcwts = Caml_ba_data_val(vx_funcwts);
	int n_rows_x = Caml_ba_array_val(vx_funcwts)->dim[0];
	int n_cols_x = Caml_ba_array_val(vx_funcwts)->dim[1];

	double *y_funcwts = Caml_ba_data_val(vy_funcwts);
	int n_rows_y = Caml_ba_array_val(vy_funcwts)->dim[0];
	int n_cols_y = Caml_ba_array_val(vy_funcwts)->dim[1];

	double *trunc_svals = Caml_ba_data_val(vtrunc_svals);
	int n_svals = Caml_ba_array_val(vtrunc_svals)->dim[0];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double *eval_pts = Caml_ba_data_val(veval_pts);
	int n_eval_pts = Caml_ba_array_val(veval_pts)->dim[0];

	double* evalres = (double *)malloc(sizeof(double) * n_eval_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	eigen_rankn_integral_eval(
		&eval_eigen_rankn_integral_eval_closure,
		&rank,
		&n_x_pts,
		&n_y_pts,
		x_pts,
		y_pts,
		x_funcwts,
		y_funcwts,
		trunc_svals,
		&overx,
		&n_integral_pts,
		integral_pts,
		integral_wts,
		&n_eval_pts,
		eval_pts,
		evalres,
		&c_err_stat,
		c_err_msg);

	/********************* Eval result           ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_pts;
	vevalres = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vevalres);
	int i, j;
	for (i = 0; i < n_eval_pts; ++i) {
		temp_ba_data[i] = evalres[i];
	}

	free(evalres);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, vevalres);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&eigen_rankn_integral_eval_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for eigen_rankn_integral_coeffs
static value eigen_rankn_integral_coeffs_closure;

double eval_eigen_rankn_integral_coeffs_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(eigen_rankn_integral_coeffs_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

// Coefficients of an integral (really, a weighted sum) of
// a rank-$n$ singular function approximation, bytecode version
CAMLprim value caml_eigen_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn)
{
	return caml_eigen_rankn_integral_coeffs_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7],
		argv[8],
		argv[9]);
}

// Coefficients of an integral (really, a weighted sum) of
// a rank-$n$ singular function approximation, native code version
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
	value vintegral_wts)
{
	CAMLparam5(vf, vrank, vx_pts, vy_pts, vx_funcwts);
	CAMLxparam5(vy_funcwts, vtrunc_svals, voverx, vintegral_pts, vintegral_wts);
	CAMLlocal3(vcoeffs, verrstat, vres);

	// Set the static closure 
	eigen_rankn_integral_coeffs_closure = vf;
	caml_register_global_root(&eigen_rankn_integral_coeffs_closure);

	int rank = Int_val(vrank);
	bool overx = Bool_val(voverx);

	double *x_pts = Caml_ba_data_val(vx_pts);
	int n_x_pts = Caml_ba_array_val(vx_pts)->dim[0];

	double *y_pts = Caml_ba_data_val(vy_pts);
	int n_y_pts = Caml_ba_array_val(vy_pts)->dim[0];

	double *x_funcwts = Caml_ba_data_val(vx_funcwts);
	int n_rows_x = Caml_ba_array_val(vx_funcwts)->dim[0];
	int n_cols_x = Caml_ba_array_val(vx_funcwts)->dim[1];

	double *y_funcwts = Caml_ba_data_val(vy_funcwts);
	int n_rows_y = Caml_ba_array_val(vy_funcwts)->dim[0];
	int n_cols_y = Caml_ba_array_val(vy_funcwts)->dim[1];

	double *trunc_svals = Caml_ba_data_val(vtrunc_svals);
	int n_svals = Caml_ba_array_val(vtrunc_svals)->dim[0];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double* coeffs = (double *)malloc(sizeof(double) * rank);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	eigen_rankn_integral_coeffs(
		&eval_eigen_rankn_integral_eval_closure,
		&rank,
		&n_x_pts,
		&n_y_pts,
		x_pts,
		y_pts,
		x_funcwts,
		y_funcwts,
		trunc_svals,
		&overx,
		&n_integral_pts,
		integral_pts,
		integral_wts,
		&rank,
		coeffs,
		&c_err_stat,
		c_err_msg);

	/********************* Coefficients result   ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = rank;
	vcoeffs = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vcoeffs);
	int i, j;
	for (i = 0; i < rank; ++i) {
		temp_ba_data[i] = coeffs[i];
	}

	free(coeffs);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, vcoeffs);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&eigen_rankn_integral_coeffs_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for num_opt_rankn_params
static value num_opt_rankn_params_closure;

double eval_num_opt_rankn_params_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(num_opt_rankn_params_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_num_opt_rankn_params_bytecode(
	value* argv,
	int argn)
{
	return caml_num_opt_rankn_params_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7]);
}

CAMLprim value caml_num_opt_rankn_params_native(
	value vf,
	value vx_lb,
	value vx_ub,
	value vy_lb,
	value vy_ub,
	value vrank,
	value vtol,
	value vmax_iter)
{
	CAMLparam5(vf, vx_lb, vx_ub, vy_lb, vy_ub);
	CAMLxparam3(vrank, vtol, vmax_iter);
	CAMLlocal5(vrho_res, vgamma_res, vvmat, vwmat, vres);
	CAMLlocal2(vborsuk_lb, verrstat);

	// Set the static closure 
	num_opt_rankn_params_closure = vf;
	caml_register_global_root(&num_opt_rankn_params_closure);

	double x_lb = Double_val(vx_lb);
	double x_ub = Double_val(vx_ub);
	double y_lb = Double_val(vy_lb);
	double y_ub = Double_val(vy_ub);
	int rank = Int_val(vrank);
	double toler = Double_val(vtol);
	int max_iter = Int_val(vmax_iter);

	double *rho_vec = (double *)malloc(sizeof(double) * (rank + 1));
	double *gamma_vec = (double *)malloc(sizeof(double) * (rank + 1));

	double* v_matrix = (double *)malloc(sizeof(double) * (rank + 1) * rank);
	double* w_matrix = (double *)malloc(sizeof(double) * (rank + 1) * rank);

	double borsuk_lb;

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	num_opt_rankn_params(
		&eval_num_opt_rankn_params_closure, 
		&x_lb, 
		&x_ub, 
		&y_lb, 
		&y_ub, 
		&rank, 
		&max_iter, 
		&toler, 
		rho_vec, 
		gamma_vec, 
		v_matrix, 
		w_matrix, 
		&borsuk_lb, 
		&c_err_stat, 
		c_err_msg);

	/**********************  Rho Results       **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat rho_ba_dim[1];
	rho_ba_dim[0] = rank + 1;
	vrho_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		rho_ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vrho_res);
	int i;
	for (i = 0; i < rank + 1; ++i) {
		temp_ba_data[i] = rho_vec[i];
	}

	/**********************  Gamma Results     **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat gamma_ba_dim[1];
	gamma_ba_dim[0] = rank + 1;
	vgamma_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		gamma_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vgamma_res);
	for (i = 0; i < rank + 1; ++i) {
		temp_ba_data[i] = gamma_vec[i];
	}

	/**********************  V matrix Results  **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat v_ba_dim[2];
	v_ba_dim[0] = rank + 1;
	v_ba_dim[1] = rank;
	vvmat = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		v_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vvmat);
	int j;
	for (i = 0; i < rank + 1; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[(rank + 1) * j + i] = v_matrix[(rank + 1) * j + i];
		}
	}

	/**********************  W matrix Results  **************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat w_ba_dim[2];
	w_ba_dim[0] = rank + 1;
	w_ba_dim[1] = rank;
	vwmat = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		2, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		w_ba_dim);

	// Reusing pointer from above
	temp_ba_data = Caml_ba_data_val(vwmat);
	for (i = 0; i < rank + 1; ++i) {
		for (j = 0; j < rank; ++j) {
			temp_ba_data[(rank + 1) * j + i] = w_matrix[(rank + 1) * j + i];
		}
	}

	free(rho_vec);
	free(gamma_vec);
	free(v_matrix);
	free(w_matrix);

	// Copy the other outputs to OCaml values
	vborsuk_lb = caml_copy_double(borsuk_lb);
	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(6, 0);

	// Fill the tuple
	Store_field(vres, 0, vrho_res);
	Store_field(vres, 1, vgamma_res);
	Store_field(vres, 2, vvmat);
	Store_field(vres, 3, vwmat);
	Store_field(vres, 4, vborsuk_lb);
	Store_field(vres, 5, verrstat);

	caml_remove_global_root(&num_opt_rankn_params_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for num_opt_rankn_eval
static value num_opt_rankn_eval_closure;

double eval_num_opt_rankn_eval_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(num_opt_rankn_eval_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_num_opt_rankn_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_num_opt_rankn_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6]);
}

CAMLprim value caml_num_opt_rankn_eval_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value veval_x_pts,
	value veval_y_pts)
{
	CAMLparam5(vf, vrho_vec, vgamma_vec, v_vmat, v_wmat);
	CAMLxparam2(veval_x_pts, veval_y_pts);
	CAMLlocal3(v_eval_res, verrstat, vres);

	// Set the static closure 
	num_opt_rankn_eval_closure = vf;
	caml_register_global_root(&num_opt_rankn_eval_closure);

	double *rho_vec = Caml_ba_data_val(vrho_vec);
	int size_of_rho = Caml_ba_array_val(vrho_vec)->dim[0];

	double *gamma_vec = Caml_ba_data_val(vgamma_vec);
	int size_of_gamma = Caml_ba_array_val(vgamma_vec)->dim[0];
	
	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *eval_x_pts = Caml_ba_data_val(veval_x_pts);
	int n_eval_x_pts = Caml_ba_array_val(veval_x_pts)->dim[0];

	double *eval_y_pts = Caml_ba_data_val(veval_y_pts);
	int n_eval_y_pts = Caml_ba_array_val(veval_y_pts)->dim[0];

	double* results = (double*)malloc(sizeof(double) * n_eval_x_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	num_opt_rankn_eval(
		&eval_num_opt_rankn_eval_closure, 
		&n_cols_v, 
		rho_vec, 
		gamma_vec, 
		vmat, 
		wmat, 
		&n_eval_x_pts, 
		eval_x_pts, 
		eval_y_pts, 
		results, 
		&c_err_stat, 
		c_err_msg);

	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_x_pts;
	v_eval_res = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(v_eval_res);
	int i;
	for (i = 0; i < n_eval_x_pts; ++i) {
		temp_ba_data[i] = results[i];
	}

	free(results);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, v_eval_res);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&num_opt_rankn_eval_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for num_opt_rankn_integral_eval
static value num_opt_rankn_integral_eval_closure;

double eval_num_opt_rankn_integral_eval_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(num_opt_rankn_integral_eval_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_num_opt_rankn_integral_eval_bytecode(
	value* argv,
	int argn)
{
	return caml_num_opt_rankn_integral_eval_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7],
		argv[8]);
}

CAMLprim value caml_num_opt_rankn_integral_eval_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts,
	value veval_pts)
{
	CAMLparam5(vf, vrho_vec, vgamma_vec, v_vmat, v_wmat);
	CAMLxparam4(voverx, vintegral_pts, vintegral_wts, veval_pts);
	CAMLlocal3(vevalres, verrstat, vres);

	// Set the static closure 
	num_opt_rankn_integral_eval_closure = vf;
	caml_register_global_root(&num_opt_rankn_integral_eval_closure);

	bool overx = Bool_val(voverx);

	double *rho_vec = Caml_ba_data_val(vrho_vec);
	int size_of_rho = Caml_ba_array_val(vrho_vec)->dim[0];

	double *gamma_vec = Caml_ba_data_val(vgamma_vec);
	int size_of_gamma = Caml_ba_array_val(vgamma_vec)->dim[0];

	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double *eval_pts = Caml_ba_data_val(veval_pts);
	int n_eval_pts = Caml_ba_array_val(veval_pts)->dim[0];

	double* evalres = (double *)malloc(sizeof(double) * n_eval_pts);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	num_opt_rankn_integral_eval(
		&eval_num_opt_rankn_integral_eval_closure, 
		&n_cols_v, 
		rho_vec, 
		gamma_vec, 
		vmat, 
		wmat, 
		&overx, 
		&n_integral_pts, 
		integral_pts, 
		integral_wts, 
		&n_eval_pts, 
		eval_pts, 
		evalres, 
		&c_err_stat, 
		c_err_msg);

	/********************* Eval result           ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_eval_pts;
	vevalres = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vevalres);
	int i, j;
	for (i = 0; i < n_eval_pts; ++i) {
		temp_ba_data[i] = evalres[i];
	}

	free(evalres);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, vevalres);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&num_opt_rankn_integral_eval_closure);

	CAMLreturn(vres);
}

// Closure and evaluator for num_opt_rankn_integral_coeffs
static value num_opt_rankn_integral_coeffs_closure;

double eval_num_opt_rankn_integral_coeffs_closure(double *x, double *y) {
	CAMLparam0();
	CAMLlocal3(vx, vy, veval);
	vx = caml_copy_double(*x);
	vy = caml_copy_double(*y);
	veval = caml_callback2(num_opt_rankn_integral_coeffs_closure, vx, vy);
	CAMLreturnT(double, Double_val(veval));
}

CAMLprim value caml_num_opt_rankn_integral_coeffs_bytecode(
	value* argv,
	int argn)
{
	return caml_num_opt_rankn_integral_coeffs_native(
		argv[0],
		argv[1],
		argv[2],
		argv[3],
		argv[4],
		argv[5],
		argv[6],
		argv[7]);
}

CAMLprim value caml_num_opt_rankn_integral_coeffs_native(
	value vf,
	value vrho_vec,
	value vgamma_vec,
	value v_vmat,
	value v_wmat,
	value voverx,
	value vintegral_pts,
	value vintegral_wts)
{
	CAMLparam5(vf, vrho_vec, vgamma_vec, v_vmat, v_wmat);
	CAMLxparam3(voverx, vintegral_pts, vintegral_wts);
	CAMLlocal3(vcoeffs, verrstat, vres);

	// Set the static closure 
	num_opt_rankn_integral_coeffs_closure = vf;
	caml_register_global_root(&num_opt_rankn_integral_coeffs_closure);

	bool overx = Bool_val(voverx);

	double *rho_vec = Caml_ba_data_val(vrho_vec);
	int size_of_rho = Caml_ba_array_val(vrho_vec)->dim[0];

	double *gamma_vec = Caml_ba_data_val(vgamma_vec);
	int size_of_gamma = Caml_ba_array_val(vgamma_vec)->dim[0];

	double *vmat = Caml_ba_data_val(v_vmat);
	int n_rows_v = Caml_ba_array_val(v_vmat)->dim[0];
	int n_cols_v = Caml_ba_array_val(v_vmat)->dim[1];

	double *wmat = Caml_ba_data_val(v_wmat);
	int n_rows_w = Caml_ba_array_val(v_wmat)->dim[0];
	int n_cols_w = Caml_ba_array_val(v_wmat)->dim[1];

	double *integral_pts = Caml_ba_data_val(vintegral_pts);
	int n_integral_pts = Caml_ba_array_val(vintegral_pts)->dim[0];

	double *integral_wts = Caml_ba_data_val(vintegral_wts);
	int n_integral_wts = Caml_ba_array_val(vintegral_wts)->dim[0];

	double* coeffs = (double *)malloc(sizeof(double) * n_cols_v);

	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	num_opt_rankn_integral_coeffs(
		&eval_num_opt_rankn_integral_coeffs_closure,
		&n_cols_v,
		rho_vec,
		gamma_vec,
		vmat,
		wmat,
		&overx,
		&n_integral_pts,
		integral_pts,
		integral_wts,
		&n_cols_v,
		coeffs,
		&c_err_stat,
		c_err_msg);

	/********************* Coefficients          ************************/
	// Only inform OCaml of the memory that actually needs to be preserved
	intnat ba_dim[1];
	ba_dim[0] = n_cols_v;
	vcoeffs = caml_ba_alloc(
		CAML_BA_FLOAT64 | CAML_BA_FORTRAN_LAYOUT,
		1, // Number of dimensions
		NULL, // This causes OCaml to allocate memory
		ba_dim);

	double *temp_ba_data = Caml_ba_data_val(vcoeffs);
	int i, j;
	for (i = 0; i < n_cols_v; ++i) {
		temp_ba_data[i] = coeffs[i];
	}

	free(coeffs);

	verrstat = Val_int(c_err_stat);

	// Allocate an OCaml tuple
	vres = caml_alloc(2, 0);

	// Fill the tuple
	Store_field(vres, 0, vcoeffs);
	Store_field(vres, 1, verrstat);

	caml_remove_global_root(&num_opt_rankn_integral_coeffs_closure);

	CAMLreturn(vres);
}