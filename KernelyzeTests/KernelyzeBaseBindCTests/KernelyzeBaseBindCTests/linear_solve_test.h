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

#ifndef LINEAR_SOLVE_C
#define LINEAR_SOLVE_C

extern "C" void linear_solve(
	const int* c_n,
	const double* c_a,
	const int* c_nrhs,
	const double* c_b,
	double* c_x,
	int* c_info);

void linear_solve_test();

#endif // LINEAR_SOLVE_C
