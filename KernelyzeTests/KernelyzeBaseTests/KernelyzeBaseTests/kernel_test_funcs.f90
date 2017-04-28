! kernel_test_funcs.f90
!
! Copyright (c) 2015 by Kernelyze LLC
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
! created on: 2015-04-25
! updated on: 2015-04-25
!
! A module of example kernels to use in testing functionality.
  
module kernel_test_funcs_mod
  
  use set_precision, only : wp
  use kernel_expprod_mod, only : kernel_expprod
  use kernel_gaussian_mod, only : kernel_gaussian
  use kernel_cauchy_mod, only : kernel_cauchy
  
  implicit none
  
  type(kernel_expprod), protected  :: exp_prod_kernel
  type(kernel_gaussian), protected :: gaussian_kernel
  type(kernel_cauchy), protected   :: cauchy_kernel
  
  contains
  
  subroutine init_test_kernels()
    call exp_prod_kernel%set_c(1E0_wp)
    call gaussian_kernel%set_c(0.5E0_wp)
    call cauchy_kernel%set_c(3E0_wp)
    call cauchy_kernel%set_alpha(1E0_wp)
  end subroutine init_test_kernels
  
end module kernel_test_funcs_mod