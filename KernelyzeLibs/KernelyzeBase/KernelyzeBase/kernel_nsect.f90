! kernel_nsect.f90
!
! Copyright (c) 2015, 2016 by Kernelyze LLC
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
! created on: 2015-07-07
! updated on: 2015-07-07
! updated on: 2016-10-07 (conform to new name for
!             integral_nsect_svar)
!
! This module contains a derived type that represents
! a rank-$n$ kernel whose $n$ component functions of $x$
! and $n$ component functions of $y$ are all sections
! of a kernel, that is, are compfunc_kernel instances.
  
module kernel_nsect_mod

  use compfunc_kernel_mod, only : singlevar_type => compfunc_kernel
  use integral_nsect_disc_mod, only : integral_type => integral_rankn_disc
  use integral_nsect_svar_mod, only : integral_sv_type => integral_nsect_svar

  include "template_kernelrankn.f90"

end module kernel_nsect_mod