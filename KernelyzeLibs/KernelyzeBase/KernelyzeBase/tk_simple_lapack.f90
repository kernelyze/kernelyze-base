! tk_simple_lapack.f90
!
! Copyright (c) 2017 by Kernelyze LLC
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
! created on: 2017-04-22
!
! This module contains simple interfaces for the LAPACK
! routines dgesv, dgbsv, dgesvd, dspev, dsytrf, and dsytrs.
  
module tk_simple_lapack_mod

  use set_precision, only : wp
  use constants_mod, only : alloc_errmsg_len, err_msg_len

  implicit none
  
  interface
    ! The dgbsv LAPACK pure subroutine
    pure subroutine dgbsv( &
        n, &
        kl, &
        ku, &
        nrhs, &
        ab, &
        ldab, &
        ipiv, &
        b, &
        ldb, &
        info)
      import wp
      integer, intent(inout)          :: n, kl, ku, nrhs, ldab, ldb, info
      real(wp), intent(inout)         :: ab(ldab, *)
      real(wp), intent(inout)         :: b(ldb, *)
      integer, intent(inout)          :: ipiv(*)
    end subroutine dgbsv
    ! The dgesv LAPACK pure subroutine
    pure subroutine dgesv( &
        n, &
        nrhs, &
        a, &
        lda, &
        ipiv, &
        b, &
        ldb, &
        info)
      import wp
      integer, intent(inout)          :: n, nrhs, lda, ldb, info
      real(wp), intent(inout)         :: a(lda, *)
      real(wp), intent(inout)         :: b(ldb, *)
      integer, intent(inout)          :: ipiv(*)
    end subroutine dgesv  
    ! The dgesvd LAPACK pure subroutine
    pure subroutine dgesvd( &
        jobu, &
        jobvt, &
        m, &
        n, &
        a, &
        lda, &
        s, &
        u, &
        ldu, &
        vt, &
        ldvt, &
        work, &
        lwork, &
        info)
      import wp
      character(len=1), intent(inout) :: jobu
      character(len=1), intent(inout) :: jobvt
      integer, intent(inout)          :: m, n, lda, ldu, ldvt, lwork, info
      real(wp), intent(inout)         :: a(lda, *)
      real(wp), intent(inout)         :: s(*)
      real(wp), intent(inout)         :: u(ldu, *)
      real(wp), intent(inout)         :: vt(ldvt, *)
      real(wp), intent(inout)         :: work(*)
    end subroutine dgesvd
    ! The dspev LAPACK pure subroutine
    pure subroutine dspev( &
        jobz, &
        uplo, &
        n, &
        ap, &
        w, &
        z, &
        ldz, &
        work, &
        info)
      import wp
      character(len=1), intent(inout) :: jobz
      character(len=1), intent(inout) :: uplo
      integer, intent(inout)          :: n, ldz, info
      real(wp), intent(inout)         :: ap(*)
      real(wp), intent(inout)         :: w(*)
      real(wp), intent(inout)         :: z(ldz,*)
      real(wp), intent(inout)         :: work(*)
    end subroutine dspev
    ! The dsytrf LAPACK pure subroutine
    pure subroutine dsytrf( &
        uplo, &
        n, &
        a, &
        lda, &
        ipiv, &
        work, &
        lwork, &
        info)
      import wp
      character(len=1), intent(inout) :: uplo
      integer, intent(inout)          :: n, lda, lwork, info
      real(wp), intent(inout)         :: a(lda, *)
      integer, intent(inout)          :: ipiv(*)
      real(wp), intent(inout)         :: work(*)
    end subroutine dsytrf  
    ! The dsytrf LAPACK pure subroutine
    pure subroutine dsytrs( &
        uplo, &
        n, &
        nrhs, &
        a, &
        lda, &
        ipiv, &
        b, &
        ldb, &
        info)
      import wp
      character(len=1), intent(in)    :: uplo
      integer, intent(inout)          :: n, nrhs, lda, ldb, info
      real(wp), intent(in)            :: a(lda, *)
      real(wp), intent(inout)         :: b(ldb, *)
      integer, intent(in)             :: ipiv(*)
    end subroutine dsytrs  
  end interface
  
  ! Interface to allow tk_gbsv to represent
  ! a call to solve with a single right-hand side
  ! vector or a set of right-hand side vectors
  ! that are packed into a $n \times nrhs$
  ! matrix.
  interface tk_gbsv
    module procedure tk_gbsv1, tk_gbsv2
  end interface tk_gbsv
  
  ! Interface to allow tk_gesv to represent
  ! a call to solve with a single right-hand side
  ! vector or a set of right-hand side vectors
  ! that are packed into a $n \times nrhs$
  ! matrix.
  interface tk_gesv
    module procedure tk_gesv1, tk_gesv2
  end interface tk_gesv
  
contains
  
pure subroutine tk_gbsv1( &
    ab, &
    b, &
    kl, &
    ipiv, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: ab(:, :)
  real(wp), intent(inout)       :: b(:)
  integer, intent(in)           :: kl
  integer, intent(inout)        :: ipiv(size(ab, 2))
  integer, intent(out)          :: info
  ! Local variables
  integer :: n, nrhs, ldab, ldb, ku, local_kl
  ! Body
  n = size(ab, 2)
  ldab = size(ab, 1)
  ldb = size(b)
  nrhs = 1
  local_kl = kl
  ku = ldab - (2 * kl + 1)
  ! Now call the LAPACK pure subroutine
  call dgbsv( &
      n, &
      local_kl, &
      ku, &
      nrhs, &
      ab, &
      ldab, &
      ipiv, &
      b, &
      ldb, &
      info)
end subroutine tk_gbsv1
    
pure subroutine tk_gbsv2( &
    ab, &
    b, &
    kl, &
    ipiv, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: ab(:, :)
  real(wp), intent(inout)       :: b(:, :)
  integer, intent(in)           :: kl
  integer, intent(inout)        :: ipiv(size(ab, 2))
  integer, intent(out)          :: info
  ! Local variables
  integer :: n, nrhs, ldab, ldb, ku, local_kl
  ! Body
  n = size(ab, 2)
  ldab = size(ab, 1)
  ldb = size(b, 1)
  nrhs = size(b, 2)
  local_kl = kl
  ku = ldab - (2 * kl + 1)
  ! Now call the LAPACK pure subroutine
  call dgbsv( &
      n, &
      local_kl, &
      ku, &
      nrhs, &
      ab, &
      ldab, &
      ipiv, &
      b, &
      ldb, &
      info)
end subroutine tk_gbsv2    
  
pure subroutine tk_gesv1( &
    a, &
    b, &
    ipiv, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: a(:, :)
  real(wp), intent(inout)       :: b(:)
  integer, intent(inout)        :: ipiv(size(a, 1))
  integer, intent(out)          :: info
  ! Local variables
  integer :: n, nrhs, lda, ldb
  ! Body
  n = size(a, 1)
  lda = size(a, 1)
  ldb = size(b)
  nrhs = 1
  ! Now call the LAPACK pure subroutine
  call dgesv( &
      n, &
      nrhs, &
      a, &
      lda, &
      ipiv, &
      b, &
      ldb, &
      info)
end subroutine tk_gesv1
    
pure subroutine tk_gesv2( &
    a, &
    b, &
    ipiv, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: a(:, :)
  real(wp), intent(inout)       :: b(:, :)
  integer, intent(inout)        :: ipiv(size(a, 1))
  integer, intent(out)          :: info
  ! Local variables
  integer :: n, nrhs, lda, ldb
  ! Body
  n = size(a, 1)
  lda = size(a, 1)
  ldb = size(b, 1)
  nrhs = size(b, 2)
  ! Now call the LAPACK pure subroutine
  call dgesv( &
      n, &
      nrhs, &
      a, &
      lda, &
      ipiv, &
      b, &
      ldb, &
      info)
end subroutine tk_gesv2    
  
pure subroutine tk_gesvd( &
    a, &
    s, &
    u, &
    vt, &
    ww, &
    info )
  ! Arguments
  real(wp), intent(inout)       :: a( :, : )
  real(wp), intent(out)         :: s( min( size(a, 1), size(a, 2) ) )
  real(wp), intent(out)         :: u( size(a, 1), size(a, 1) )
  real(wp), intent(out)         :: vt( size(a, 2), size(a, 2) )
  real(wp), intent(out)         :: ww( size(s) - 1 )
  integer, intent(out)          :: info
  ! Local variables
  real(wp), allocatable :: work( : )
  character(len=1)      :: jobu, jobvt
  integer               :: m, n, lda, ldu, ldvt, lwork
  integer               :: i
  ! Body
  m = size(a, 1)
  n = size(a, 2)
  ! Use m * n or quadruple the minimum size, whichever is larger
  lwork = max(m * n, 4 * max(3*min(m, n)+max(m, n), 5*min(m,n)) )
  lda = m
  ldu = m
  ldvt = n
  jobu = 'A'
  jobvt = 'A'
  ! Allocate the work array
  allocate( work(lwork), stat = info )
  if (info /= 0) then
    return
  end if
  ! Now call the LAPACK pure subroutine
  call dgesvd( &
      jobu, &
      jobvt, &
      m, &
      n, &
      a, &
      lda, &
      s, &
      u, &
      ldu, &
      vt, &
      ldvt, &
      work, &
      lwork, &
      info)
  do i = 2, min(m, n)
    ww(i-1) = work(i)
  end do
end subroutine tk_gesvd  

pure subroutine tk_spev( &
    ap, &
    w, &
    uplo, &
    z, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: ap(:)
  real(wp), intent(out)         :: w(:)
  character(len=1), intent(in)  :: uplo
  real(wp), intent(out)         :: z(size(w),size(w))
  integer, intent(out)          :: info
  ! Local variables
  real(wp)  :: work( 4 * size(w) )
  character(len=1)  :: jobz
  character(len=1)  :: local_uplo
  integer           :: n, ldz
  ! Body
  jobz = 'V'
  local_uplo = uplo
  n = size(w)
  ldz = size(z, 1)
  ! Now call the LAPACK pure subroutine
  call dspev( &
    jobz = jobz, &
    uplo = local_uplo, &
    n = n, &
    ap = ap, &
    w = w, &
    z = z, &
    ldz = ldz, &
    work = work, &
    info = info )
end subroutine tk_spev
    
pure subroutine tk_sytrf( &
    a, &
    uplo, &
    ipiv, &
    info)
  ! Arguments
  real(wp), intent(inout)       :: a(:, :)
  character(len=1), intent(in)  :: uplo
  integer, intent(out)          :: ipiv(size(a, 1))
  integer, intent(out)          :: info
  ! Local variables
  real(wp)  :: work( 64 * size(a, 1) )
  character(len=1)  :: local_uplo
  integer           :: n, lda, lwork
  ! Body
  local_uplo = uplo
  n = size(a, 1)
  lda = size(a, 1)
  lwork = size(work)
  ! Now call the LAPACK pure subroutine
  call dsytrf( &
      local_uplo, &
      n, &
      a, &
      lda, &
      ipiv, &
      work, &
      lwork, &
      info)
end subroutine tk_sytrf
    
pure subroutine tk_sytrs( &
    a, &
    b, &
    ipiv, &
    uplo, &
    info)
  ! Arguments
  real(wp), intent(in)          :: a(:, :)
  real(wp), intent(inout)       :: b(:)
  character(len=1), intent(in)  :: uplo
  integer, intent(in)           :: ipiv(size(a, 1))
  integer, intent(out)          :: info
  ! Local variables
  integer :: n, nrhs, lda, ldb
  ! Body
  n = size(a, 1)
  lda = size(a, 1)
  ldb = size(b)
  nrhs = 1
  ! Now call the LAPACK pure subroutine
  call dsytrs( &
      uplo, &
      n, &
      nrhs, &
      a, &
      lda, &
      ipiv, &
      b, &
      ldb, &
      info)
end subroutine tk_sytrs 
    
end module tk_simple_lapack_mod