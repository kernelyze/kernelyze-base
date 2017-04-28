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

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "option_formulae_test.h"

void option_formulae_test() {
	
	static const double pi = 4.0 * atan(1.0);

	// Under the Bachelier (normal) model, an ATMF straddle 
	// should be worth vol * sqrt(2 / pi) * discfac
	bool is_call = true;
	double strike = 42.0;
	double forward = 42.0;
	double vol = 2.0;
	double disc_fac = 1.0;

	double atmf_straddle = normopt_formula(
		&is_call, &strike, &forward, &vol, &disc_fac);
	is_call = false;
	atmf_straddle += normopt_formula(
		&is_call, &strike, &forward, &vol, &disc_fac);
	
	double expected_atmf_straddle = vol * sqrt(2.0 / pi);

	printf("Normal-distribution ATMF straddle: %f; expected: %f\n",
		atmf_straddle, expected_atmf_straddle);

	// Check the delta of a Black-model option that has d1 == 0.0
	is_call = true;
	vol = 0.5;
	strike = 100.0;
	forward = strike * exp(-0.5 * vol * vol);
	disc_fac = 1.0;
	double bump_amt = 0.0001;
	// double atmf_call    = black_formula(&is_call, &strike, &forward, &vol, &disc_fac);
	forward = forward + bump_amt;
	double atmf_call_up = black_formula(&is_call, &strike, &forward, &vol, &disc_fac);
	forward = forward - 2.0 * bump_amt;
	double atmf_call_dn = black_formula(&is_call, &strike, &forward, &vol, &disc_fac);
	double bump_delta = (atmf_call_up - atmf_call_dn) / ( 2.0 * bump_amt );

	printf("Black model call with d1 == 0.0; delta should be 0.5, bump delta is %f\n",
		bump_delta);
}
