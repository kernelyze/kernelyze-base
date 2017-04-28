! kernel_neigen.f90
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
! created on: 2016-12-01 (am)
! updated on: 2016-12-04 (simple updates)
!
! This module contains a derived type that represents
! a rank-$n$ kernel whose $n$ component functions of $x$
! and $n$ component functions of $y$ are all eigenfunctions,
! but are represented by integral_kernel_disc types because
! these are the result of a Nystrom approach to finding
! the eigenfunctions numerically: sample the kernel at
! a vector of x points and a vector of y points, then
! SVD the resulting matrix, then compute
! f_i \left( x \right) = \sum_{j=1}^{n} K \left( x , y_j \right) g_j^{singvec} 
  
module kernel_neigen_mod

  use integral_kernel_disc_mod, only : singlevar_type => integral_kernel_disc
  use integral_neigen_disc_mod, only : integral_type => integral_rankn_disc
  use integral_neigen_svar_mod, only : integral_sv_type => integral_neigen_svar

  include "template_kernelrankn.f90"

end module kernel_neigen_mod