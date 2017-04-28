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
// updated: 2017-04-14

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "brent_test.h"

double sinfunc(double *x) {
	return sin(*x);
}

double quadratic_func(double* x) {
	return (*x) * ((*x) - 1.0);
}

void brent_test() {
	double a = -1.0;
	double b = 1.0;
	double toler = 1e-15;

	double testeval = sinfunc(&b);

	printf("The sin(x) evaluation at 1.0 is: %f\n", testeval);

	double zero_of_sin = zeroin(&a, &b, &sinfunc, &toler);

	printf("zeroin of sin(x) on [-1.0, 1.0] should be 0, is: %f\n", 
		zero_of_sin);

	int max_iter = 100;
	int neval;
	int err_code;
	double x_zero;
	double f_zero;

	brent_zero(
		&a, 
		&b, 
		&sinfunc, 
		&toler, 
		&max_iter, 
		&neval, 
		&err_code, 
		&x_zero, 
		&f_zero);

	printf("brent_zero of sin(x), x in [-1, 1]:\n");
	printf("root %f, value at root %f, %d iterations, %d error code\n",
		x_zero, f_zero, neval, err_code);

	testeval = quadratic_func(&a);
	printf("The x * (x - 1.0) evaluation at -1.0 is %f\n", testeval);

	double x_opt;
	double f_opt;
	int c_err_stat;
	char c_err_msg[256];
	c_err_msg[0] = '\0';

	f_min(
		&a,
		&b,
		&quadratic_func,
		&toler,
		&x_opt,
		&f_opt,
		&c_err_stat,
		c_err_msg);

	printf("f_min of x * (x - 1.0), x in [-1, 1]:\n");
	printf("min %f, value at min %f, err code %d, err msg %s\n",
		x_opt, f_opt, c_err_stat, c_err_msg);

	brent_min(
		&a,
		&b,
		&quadratic_func,
		&toler,
		&max_iter,
		&neval,
		&err_code,
		&x_opt,
		&f_opt);

	printf("brent_min of x * (x - 1.0), x in [-1, 1]:\n");
	printf("min %f, value at min %f, %d iterations, %d error code\n",
		x_opt, f_opt, neval, err_code);
}