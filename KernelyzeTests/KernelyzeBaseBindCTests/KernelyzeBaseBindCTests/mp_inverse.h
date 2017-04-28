// Copyright (c) 2015 by Kernelyze LLC
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
// created:   2015-06-03
// updated:   2015-06-03

#ifndef MP_INVERSE_FOR_C
#define MP_INVERSE_FOR_C

extern "C" void mp_inverse(double* to_return, double* to_invert, int nrows, int ncols, double toler);

#endif // MP_INVERSE_FOR_C