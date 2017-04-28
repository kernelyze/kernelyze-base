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

#ifndef BRENT_FOR_C
#define BRENT_FOR_C

extern "C" double zeroin(
	double* a, 
	double* b, 
	double (*f)(double *x), 
	double* toler);

extern "C" void brent_zero(
	double* a, 
	double* b, 
	double(*f)(double *x), 
	double* toler, 
	int* max_iter, 
	int* neval, 
	int* err_code, 
	double* x_zero, 
	double* f_zero);

extern "C" void f_min(
	double* a,
	double* b,
	double(*f)(double *x),
	double* tol,
	double* x_opt,
	double* f_opt,
	int* err_stat,
	char* err_msg);

extern "C" void brent_min(
	double* a,
	double* b,
	double(*f)(double *x),
	double* toler,
	int* max_iter,
	int* neval,
	int* err_code,
	double* x_min,
	double* f_min);

void brent_test();

#endif // BRENT_FOR_C