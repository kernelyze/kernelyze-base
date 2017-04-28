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
// updated: 2017-04-23 (added Cauchy, Black, and Bachelier kernels)

#include <math.h>
#include "test_function_ptrs.h"
#include "option_formulae_test.h"

double expprod(double *x, double *y) {
	return exp((*x) * (*y));
}

double gaussian(double* x, double* y) {
	static const double pi = 4.0 * atan(1.0);
	return (1.0 / sqrt(2 * pi)) * exp(-0.5 * ((*x) - (*y)) * ((*x) - (*y)));
}

double cauchy(double* x, double* y) {
	return ( 1.0 / ( (*x) - (*y) ) );
}

// The x argument is the forward, the y argument is the strike
double bachelier(double* x, double* y) {
	// Some reasonable test parameters for a call;
	// for the normal (Bachelier) option-pricing
	// formula, the vol is in absolute terms, so 
	// the "1.0" vol below is in accord with pricing
	// if the standard deviation of the underlier
	// over the life of the option ($\sigma \sqrt{T}$)
	// is 1.0.  In pricing for rates, for example,
	// this would be a 100 bp vol (and, of course,
	// to actually compute the price one would need
	// an annuity DV01 for a typical rates option).
	// The discount factor is 1.0, so this example
	// has no discounting.
	static const bool is_call = true;
	static const double vol = 1.0;
	static const double discfac = 1.0;
	return normopt_formula(&is_call, y, x, &vol, &discfac);
}

// The x argument is the forward, the y argument is the strike
double black(double* x, double* y) {
	// Some reasonable test parameters for a call;
	// the vol of 0.5 corresponds to a 50% vol, 
	// which is high but not totally outlandish.
	// The discount factor is 1.0, so this example
	// has no discounting.
	static const bool is_call = true;
	static const double vol = 0.5;
	static const double discfac = 1.0;
	return black_formula(&is_call, y, x, &vol, &discfac);
}
