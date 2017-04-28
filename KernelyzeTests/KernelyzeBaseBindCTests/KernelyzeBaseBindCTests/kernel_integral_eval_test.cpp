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
#include "test_function_ptrs.h"
#include "kernel_integral_eval_test.h"

void kernel_integral_eval_test() {
	static const int num_integral_pts = 21;
	double integral_wts[num_integral_pts];
	double integral_pts[num_integral_pts];

	static const int num_eval_pts = 11;
	double eval_pts[num_eval_pts];
	double eval_results[num_eval_pts];
	double eval_results_by_hand[num_eval_pts];

	bool is_over_x = true;

	// Evenly-spaced points from -1.0 to 1.0, with
	// equal weights
	for (int i = 0; i < num_integral_pts; ++i) {
		integral_wts[i] = 1.0 / ((double)num_integral_pts);
		integral_pts[i] = -1.0 + 2.0 * ((double)i) / ((double)(num_integral_pts-1));
	}

	// Test at similarly evenly-spaced points
	for (int i = 0; i < num_eval_pts; ++i) {
		eval_pts[i] = -1.0 + 2.0 * ((double)i) / ((double)(num_eval_pts - 1));
	}

	// Compute the integral (really a sum) using the library:
	kernel_integral_eval(
		&expprod, 
		&is_over_x, 
		&num_integral_pts, 
		integral_pts, 
		integral_wts, 
		&num_eval_pts, 
		eval_pts, 
		eval_results);

	for (int i = 0; i < num_eval_pts; ++i) {
		eval_results_by_hand[i] = 0.0;
		for (int j = 0; j < num_integral_pts; ++j) {
			eval_results_by_hand[i] +=
				integral_wts[j] * expprod(&integral_pts[j], &eval_pts[i]);
		}
		printf("Kernel discrete integral (sum) from library at %f is %f; by hand it is %f\n",
			eval_pts[i], eval_results[i], eval_results_by_hand[i]);
	}
}
