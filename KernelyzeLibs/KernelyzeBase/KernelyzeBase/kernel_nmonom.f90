! kernel_nmonom.f90
!
! Copyright (c) 2016 by Kernelyze LLC
! Author: Thomas A. Knox
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Affero General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Affero General Public License for more details.
! 
! You should have received a copy of the GNU Affero General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! created on: 2016-02-07
! updated on: 2016-02-07
! updated on: 2016-08-22 (added support for integrals
!             against singlevar_funcs)
! updated on: 2016-10-07 (conform to new name for
!             integral_nmonom_svar)
!
! This module contains a derived type that represents
! a rank-$n$ kernel whose $n$ component functions of $x$
! and $n$ component functions of $y$ are all monomials.
  
module kernel_nmonom_mod

  use monomial_mod, only : singlevar_type => monomial
  use integral_nmonom_disc_mod, only : integral_type => integral_rankn_disc
  use integral_nmonom_svar_mod, only : integral_sv_type => integral_nmonom_svar

  include "template_kernelrankn.f90"

end module kernel_nmonom_mod