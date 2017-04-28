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

#ifndef OPTION_FORMULAE_C
#define OPTION_FORMULAE_C

extern "C" double black_formula(
	const bool* c_is_call,
	const double* c_strike,
	const double* c_forward,
	const double* c_vol,
	const double* c_disc_fac);

extern "C" double normopt_formula(
	const bool* c_is_call,
	const double* c_strike,
	const double* c_forward,
	const double* c_vol,
	const double* c_disc_fac);

void option_formulae_test();

#endif // OPTION_FORMULAE_C