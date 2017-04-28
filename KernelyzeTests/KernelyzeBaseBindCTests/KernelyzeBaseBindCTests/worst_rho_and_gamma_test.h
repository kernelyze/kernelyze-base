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

#ifndef WORST_RHO_AND_GAMMA_FOR_C
#define WORST_RHO_AND_GAMMA_FOR_C

extern "C" void worst_rho_and_gamma(
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

void worst_rho_and_gamma_test();

#endif // WORST_RHO_AND_GAMMA_FOR_C
