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

#ifndef FIND_REL_OPTIMA_FOR_C
#define FIND_REL_OPTIMA_FOR_C

extern "C" void find_rel_optima(
	double* result_optima,
	double(*f)(double *x),
	int* n_grid_pts,
	double* grid,
	double* toler,
	int* err_stat,
	char* err_msg);

void find_rel_optima_test();

#endif // FIND_REL_OPTIMA_FOR_C
