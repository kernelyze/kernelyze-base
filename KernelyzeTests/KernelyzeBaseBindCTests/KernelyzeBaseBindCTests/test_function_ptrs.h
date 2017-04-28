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

#ifndef TEST_FUN_PTRS
#define TEST_FUN_PTRS

double expprod(double *x, double *y);

double gaussian(double* x, double* y);

double cauchy(double* x, double* y);

double bachelier(double* x, double* y);

double black(double* x, double* y);

#endif // TEST_FUN_PTRS
