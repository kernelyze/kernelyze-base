! normal_distn.f90
!
! Copyright (c) 2016, 2017 by Kernelyze LLC
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
! created on: 2016-08-21
! updated on: 2017-04-23 (remove use statements for
!             modules or parts of modules that are
!             not actually used)
!
! This module contains functions to evaluate the
! probability density function and cumulative distribution
! function of the univariate normal (Gaussian) distribution.
  
module normal_distn

use set_precision, only : wp
use constants_mod, only : inv_sqrt2pi, sqrt2

implicit none

private

public  :: normal_pdf
public  :: normal_cdf
! Elemental versions
public  :: normal_pdf_elt
public  :: normal_cdf_elt
! Standard normal (mean of zero, standard deviation of one)
public  :: std_normal_pdf
public  :: std_normal_cdf
! Elemental versions for standard normal
public  :: std_normal_pdf_elt
public  :: std_normal_cdf_elt

contains
  
pure function normal_pdf(x, mu, sigma) result(pdf)
  ! Arguments
  real(wp), intent(in)  :: x
  real(wp), intent(in)  :: mu
  real(wp), intent(in)  :: sigma
  ! Result
  real(wp)              :: pdf
  ! Body
  pdf = ( inv_sqrt2pi / sigma ) &
        * exp( - (1E0_wp / (2E0_wp * sigma * sigma)) * (x - mu) * (x - mu) )
end function normal_pdf

elemental function normal_pdf_elt(x, mu, sigma) result(pdf)
  ! Arguments
  real(wp), intent(in)  :: x
  real(wp), intent(in)  :: mu
  real(wp), intent(in)  :: sigma
  ! Result
  real(wp)              :: pdf
  ! Body
  pdf = normal_pdf(x, mu, sigma)
end function normal_pdf_elt

pure function normal_cdf(x, mu, sigma) result(cdf)
  ! Arguments
  real(wp), intent(in)  :: x
  real(wp), intent(in)  :: mu
  real(wp), intent(in)  :: sigma
  ! Result
  real(wp)              :: cdf
  ! Body
  ! Use the Fortran intrinsic function erfc and the relation
  ! standard_normal_cdf(z) = 0.5 * erfc( - z / sqrt(2) )
  cdf = 5E-1_wp * erfc( - ( x - mu ) / (sigma * sqrt2) )
end function normal_cdf

elemental function normal_cdf_elt(x, mu, sigma) result(cdf)
  ! Arguments
  real(wp), intent(in)  :: x
  real(wp), intent(in)  :: mu
  real(wp), intent(in)  :: sigma
  ! Result
  real(wp)              :: cdf
  ! Body
  cdf = normal_cdf(x, mu, sigma)
end function normal_cdf_elt

pure function std_normal_pdf(x) result(pdf)
  ! Arguments
  real(wp), intent(in)  :: x
  ! Result
  real(wp)              :: pdf
  ! Body
  ! Could improve performance a bit by implementing this directly,
  ! but that would create code that is much harder to maintain.
  pdf = normal_pdf( x , 0E0_wp , 1E0_wp )
end function std_normal_pdf

elemental function std_normal_pdf_elt(x) result(pdf)
  ! Arguments
  real(wp), intent(in)  :: x
  ! Result
  real(wp)              :: pdf
  ! Body
  pdf = std_normal_pdf( x )
end function std_normal_pdf_elt

pure function std_normal_cdf(x) result(cdf)
  ! Arguments
  real(wp), intent(in)  :: x
  ! Result
  real(wp)              :: cdf
  ! Body
  ! Could improve performance a bit by implementing this directly,
  ! but that would create code that is much harder to maintain.
  cdf = normal_cdf( x , 0E0_wp , 1E0_wp )
end function std_normal_cdf

elemental function std_normal_cdf_elt(x) result(cdf)
  ! Arguments
  real(wp), intent(in)  :: x
  ! Result
  real(wp)              :: cdf
  ! Body
  cdf = std_normal_cdf( x )
end function std_normal_cdf_elt

end module normal_distn