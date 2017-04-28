! unirnk.f90
!
! Based on Michel Olagnon's orderpack.f90 (see
! www.fortran-2000.com, "A WEB site on the way 
! to more public domain utilities")
!
! Modifications are Copyright (c) 2015 by Kernelyze LLC
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
! created on: 2015-02-21
! updated on: 2015-04-23
!
! A module containing unique ranking functionality,
! created using Michel Olagnon's orderpack.f90

module m_unirnk

use set_precision, only : wp
use constants_mod, only : alloc_errmsg_len, err_msg_len
use , intrinsic :: ieee_arithmetic, only : ieee_is_nan

implicit none

public :: unirnk
private :: i_unirnk, d_unirnk, d_unirnk_no_nan
private :: i_nearless, d_nearless, nearless
interface unirnk
  module procedure d_unirnk, i_unirnk
end interface unirnk
interface nearless
  module procedure d_nearless, i_nearless
end interface nearless

contains
  
!   TAK: Treat NaN as a value; since NaN is not equal to
!   anything, even itself, put all NaNs at the end of the
!   array of unique values, which includes each NaN as
!   its own unique value.
pure subroutine d_unirnk( &
    xvalt, &
    irngt, &
    nuni, &
    err_stat, &
    err_msg)

  ! Arguments
  real(wp), dimension (:), intent (in)  :: xvalt
  integer, dimension (:), intent (out)  :: irngt
  integer, intent (out)                 :: nuni
  ! Optional arguments to handle errors if desired
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  
  ! Local variables
  real(wp) , dimension(:) , allocatable ::  xvalt_no_nan
  integer , dimension(:) , allocatable  ::  ind_xvalt_nans, ind_xvalt_non_nans
  integer                               ::  num_nans, j, nan_k, nonan_k
  
  ! Local variables to handle any allocation problems:
  character(*), parameter         :: proc_name = "Unique ranking for floating-point arrays: "
  character(len=alloc_errmsg_len) :: alloc_err_msg
  integer                         :: alloc_stat
  
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  
  ! Extract any NaNs in the input array (xvalt)
  num_nans = count(ieee_is_nan(xvalt))
  if (num_nans > 0) then
    
    ! Allocate the helper arrays
    allocate(xvalt_no_nan(size(xvalt) - num_nans), stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If it was passed, fill in the err_msg appropriately
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set irngt to -1
      irngt = -1
      ! Return
      return
    end if
    allocate(ind_xvalt_nans(num_nans), stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If it was passed, fill in the err_msg appropriately
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set irngt to -1
      irngt = -1
      ! Return
      return
    end if
    allocate(ind_xvalt_non_nans(size(xvalt_no_nan)), stat = alloc_stat, errmsg = alloc_err_msg)
    if (alloc_stat /= 0) then
      ! If I am here, there was an allocation problem.
      ! If it was passed, fill in the err_msg appropriately
      if (present(err_stat)) then
        err_stat = alloc_stat
      end if
      if (present(err_msg)) then
        err_msg = proc_name // alloc_err_msg
      end if
      ! Set irngt to -1
      irngt = -1
      ! Return
      return
    end if
    
    ! Find the indices of NaN and non-NaN elements of the xvalt array
    nan_k = 1
    nonan_k = 1
    do j = 1,size(xvalt)
      if (ieee_is_nan(xvalt(j))) then
        ind_xvalt_nans(nan_k) = j
        nan_k = nan_k + 1
      else
        ind_xvalt_non_nans(nonan_k) = j
        nonan_k = nonan_k + 1
      end if
    end do
    ! Put all non-NaN values of xvalt into
    ! xvalt_no_nan.  This is easily done with
    ! a pack statement, but since I need the
    ! indices of these elements in the original
    ! xvalt array anyway, it's even easier to
    ! use those indices directly.
    xvalt_no_nan = xvalt(ind_xvalt_non_nans)
    ! Now call the routine that assumes no NaNs are
    ! present in the array, using the "NaNs excluded"
    ! array (xvalt_no_nan)
    call d_unirnk_no_nan(xvalt_no_nan , irngt , nuni)
    ! Use a concise vector-of-indices statement
    ! to map the indices produced (which are indices
    ! in the non-NaN array) back to indices in the 
    ! original array.
    irngt(1:nuni) = ind_xvalt_non_nans(irngt(1:nuni))
    ! Account for the NaN values
    irngt((nuni + 1):) = ind_xvalt_nans(1:min(size(ind_xvalt_nans), &
        (size(irngt) - nuni)))

    ! Update the number of unique entries to
    ! account for the number of NaNs
    nuni = nuni + num_nans
  else
    ! Call the routine that assumes no NaNs are
    ! present in the array, since none are
    call d_unirnk_no_nan(xvalt , irngt , nuni)
  end if

end subroutine d_unirnk

pure subroutine d_unirnk_no_nan (xvalt, irngt, nuni)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      real (kind=wp), dimension (:), intent (in) :: xvalt
      integer, dimension (:), intent (out) :: irngt
      integer, intent (out) :: nuni
! __________________________________________________________
      integer, dimension (size(irngt)) :: jwrkt
      integer :: lmtna, lmtnc, irng, irng1, irng2
      integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
      real (kind=wp) :: xtst, xvala, xvalb
!
!
      nval = min (size(xvalt), size(irngt))
      nuni = nval
!
      select case (nval)
      case (:0)
         return
      case (1)
         irngt (1) = 1
         return
      case default
         continue
      end select
!
!  Fill-in the index array, creating ordered couples
!
      do iind = 2, nval, 2
         if (xvalt(iind-1) < xvalt(iind)) then
            irngt (iind-1) = iind - 1
            irngt (iind) = iind
         else
            irngt (iind-1) = iind
            irngt (iind) = iind - 1
         end if
      end do
      if (modulo(nval, 2) /= 0) then
         irngt (nval) = nval
      end if
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      lmtna = 2
      lmtnc = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      do
         if (nval <= 4) exit
!
!   Loop on merges of A and B into C
!
         do iwrkd = 0, nval - 1, 4
            if ((iwrkd+4) > nval) then
               if ((iwrkd+2) >= nval) exit
!
!   1 2 3
!
               if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
!
!   1 3 2
!
               if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                  irng2 = irngt (iwrkd+2)
                  irngt (iwrkd+2) = irngt (iwrkd+3)
                  irngt (iwrkd+3) = irng2
!
!   3 1 2
!
               else
                  irng1 = irngt (iwrkd+1)
                  irngt (iwrkd+1) = irngt (iwrkd+3)
                  irngt (iwrkd+3) = irngt (iwrkd+2)
                  irngt (iwrkd+2) = irng1
               end if
               exit
            end if
!
!   1 2 3 4
!
            if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
!
!   1 3 x x
!
            if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
               irng2 = irngt (iwrkd+2)
               irngt (iwrkd+2) = irngt (iwrkd+3)
               if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
!   1 3 2 4
                  irngt (iwrkd+3) = irng2
               else
!   1 3 4 2
                  irngt (iwrkd+3) = irngt (iwrkd+4)
                  irngt (iwrkd+4) = irng2
               end if
!
!   3 x x x
!
            else
               irng1 = irngt (iwrkd+1)
               irng2 = irngt (iwrkd+2)
               irngt (iwrkd+1) = irngt (iwrkd+3)
               if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
                  irngt (iwrkd+2) = irng1
                  if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
!   3 1 2 4
                     irngt (iwrkd+3) = irng2
                  else
!   3 1 4 2
                     irngt (iwrkd+3) = irngt (iwrkd+4)
                     irngt (iwrkd+4) = irng2
                  end if
               else
!   3 4 1 2
                  irngt (iwrkd+2) = irngt (iwrkd+4)
                  irngt (iwrkd+3) = irng1
                  irngt (iwrkd+4) = irng2
               end if
            end if
         end do
!
!  The Cs become As and Bs
!
         lmtna = 4
         exit
      end do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      do
         if (2*lmtna >= nval) exit
         iwrkf = 0
         lmtnc = 2 * lmtnc
!
!   Loop on merges of A and B into C
!
         do
            iwrk = iwrkf
            iwrkd = iwrkf + 1
            jinda = iwrkf + lmtna
            iwrkf = iwrkf + lmtnc
            if (iwrkf >= nval) then
               if (jinda >= nval) exit
               iwrkf = nval
            end if
            iinda = 1
            iindb = jinda + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            jwrkt (1:lmtna) = irngt (iwrkd:jinda)
            xvala = xvalt (jwrkt(iinda))
            xvalb = xvalt (irngt(iindb))
!
            do
               iwrk = iwrk + 1
!
!  We still have unprocessed values in both A and B
!
               if (xvala > xvalb) then
                  irngt (iwrk) = irngt (iindb)
                  iindb = iindb + 1
                  if (iindb > iwrkf) then
!  Only A still with unprocessed values
                     irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                     exit
                  end if
                  xvalb = xvalt (irngt(iindb))
               else
                  irngt (iwrk) = jwrkt (iinda)
                  iinda = iinda + 1
                  if (iinda > lmtna) exit! Only B still with unprocessed values
                  xvala = xvalt (jwrkt(iinda))
               end if
!
            end do
         end do
!
!  The Cs become As and Bs
!
         lmtna = 2 * lmtna
      end do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      iinda = 1
      iindb = lmtna + 1
      nuni = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      jwrkt (1:lmtna) = irngt (1:lmtna)
      if (iindb <= nval) then
        xtst = nearless (min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
      else
        xtst = nearless (xvalt(jwrkt(1)))
      endif
      do iwrk = 1, nval
!
!  We still have unprocessed values in both A and B
!
         if (iinda <= lmtna) then
            if (iindb <= nval) then
               if (xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
                  irng = irngt (iindb)
                  iindb = iindb + 1
               else
                  irng = jwrkt (iinda)
                  iinda = iinda + 1
               end if
            else
!
!  Only A still with unprocessed values
!
               irng = jwrkt (iinda)
               iinda = iinda + 1
            end if
         else
!
!  Only B still with unprocessed values
!
            irng = irngt (iwrk)
         end if
         if (xvalt(irng) > xtst) then
            xtst = xvalt (irng)
            nuni = nuni + 1
            irngt (nuni) = irng
         end if
!
      end do
!
      return
!
end subroutine d_unirnk_no_nan

pure subroutine i_unirnk (xvalt, irngt, nuni)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      integer, dimension (:), intent (in) :: xvalt
      integer, dimension (:), intent (out) :: irngt
      integer, intent (out) :: nuni
! __________________________________________________________
      integer, dimension (size(irngt)) :: jwrkt
      integer :: lmtna, lmtnc, irng, irng1, irng2
      integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
      integer :: xtst, xvala, xvalb
!
!
      nval = min (size(xvalt), size(irngt))
      nuni = nval
!
      select case (nval)
      case (:0)
         return
      case (1)
         irngt (1) = 1
         return
      case default
         continue
      end select
!
!  Fill-in the index array, creating ordered couples
!
      do iind = 2, nval, 2
         if (xvalt(iind-1) < xvalt(iind)) then
            irngt (iind-1) = iind - 1
            irngt (iind) = iind
         else
            irngt (iind-1) = iind
            irngt (iind) = iind - 1
         end if
      end do
      if (modulo(nval, 2) /= 0) then
         irngt (nval) = nval
      end if
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      lmtna = 2
      lmtnc = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      do
         if (nval <= 4) exit
!
!   Loop on merges of A and B into C
!
         do iwrkd = 0, nval - 1, 4
            if ((iwrkd+4) > nval) then
               if ((iwrkd+2) >= nval) exit
!
!   1 2 3
!
               if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
!
!   1 3 2
!
               if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                  irng2 = irngt (iwrkd+2)
                  irngt (iwrkd+2) = irngt (iwrkd+3)
                  irngt (iwrkd+3) = irng2
!
!   3 1 2
!
               else
                  irng1 = irngt (iwrkd+1)
                  irngt (iwrkd+1) = irngt (iwrkd+3)
                  irngt (iwrkd+3) = irngt (iwrkd+2)
                  irngt (iwrkd+2) = irng1
               end if
               exit
            end if
!
!   1 2 3 4
!
            if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
!
!   1 3 x x
!
            if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
               irng2 = irngt (iwrkd+2)
               irngt (iwrkd+2) = irngt (iwrkd+3)
               if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
!   1 3 2 4
                  irngt (iwrkd+3) = irng2
               else
!   1 3 4 2
                  irngt (iwrkd+3) = irngt (iwrkd+4)
                  irngt (iwrkd+4) = irng2
               end if
!
!   3 x x x
!
            else
               irng1 = irngt (iwrkd+1)
               irng2 = irngt (iwrkd+2)
               irngt (iwrkd+1) = irngt (iwrkd+3)
               if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
                  irngt (iwrkd+2) = irng1
                  if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
!   3 1 2 4
                     irngt (iwrkd+3) = irng2
                  else
!   3 1 4 2
                     irngt (iwrkd+3) = irngt (iwrkd+4)
                     irngt (iwrkd+4) = irng2
                  end if
               else
!   3 4 1 2
                  irngt (iwrkd+2) = irngt (iwrkd+4)
                  irngt (iwrkd+3) = irng1
                  irngt (iwrkd+4) = irng2
               end if
            end if
         end do
!
!  The Cs become As and Bs
!
         lmtna = 4
         exit
      end do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      do
         if (2*lmtna >= nval) exit
         iwrkf = 0
         lmtnc = 2 * lmtnc
!
!   Loop on merges of A and B into C
!
         do
            iwrk = iwrkf
            iwrkd = iwrkf + 1
            jinda = iwrkf + lmtna
            iwrkf = iwrkf + lmtnc
            if (iwrkf >= nval) then
               if (jinda >= nval) exit
               iwrkf = nval
            end if
            iinda = 1
            iindb = jinda + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            jwrkt (1:lmtna) = irngt (iwrkd:jinda)
            xvala = xvalt (jwrkt(iinda))
            xvalb = xvalt (irngt(iindb))
!
            do
               iwrk = iwrk + 1
!
!  We still have unprocessed values in both A and B
!
               if (xvala > xvalb) then
                  irngt (iwrk) = irngt (iindb)
                  iindb = iindb + 1
                  if (iindb > iwrkf) then
!  Only A still with unprocessed values
                     irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                     exit
                  end if
                  xvalb = xvalt (irngt(iindb))
               else
                  irngt (iwrk) = jwrkt (iinda)
                  iinda = iinda + 1
                  if (iinda > lmtna) exit! Only B still with unprocessed values
                  xvala = xvalt (jwrkt(iinda))
               end if
!
            end do
         end do
!
!  The Cs become As and Bs
!
         lmtna = 2 * lmtna
      end do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      iinda = 1
      iindb = lmtna + 1
      nuni = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      jwrkt (1:lmtna) = irngt (1:lmtna)
      if (iindb <= nval) then
        xtst = nearless (min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
      else
        xtst = nearless (xvalt(jwrkt(1)))
      endif
      do iwrk = 1, nval
!
!  We still have unprocessed values in both A and B
!
         if (iinda <= lmtna) then
            if (iindb <= nval) then
               if (xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
                  irng = irngt (iindb)
                  iindb = iindb + 1
               else
                  irng = jwrkt (iinda)
                  iinda = iinda + 1
               end if
            else
!
!  Only A still with unprocessed values
!
               irng = jwrkt (iinda)
               iinda = iinda + 1
            end if
         else
!
!  Only B still with unprocessed values
!
            irng = irngt (iwrk)
         end if
         if (xvalt(irng) > xtst) then
            xtst = xvalt (irng)
            nuni = nuni + 1
            irngt (nuni) = irng
         end if
!
      end do
!
      return
!
end subroutine i_unirnk

pure function d_nearless (xval) result (d_nl)
!  Nearest value less than given value
! __________________________________________________________
      real (kind=wp), intent (in) :: xval
      real (kind=wp) :: d_nl
! __________________________________________________________
      d_nl = nearest (xval, -1.0_wp)
      return
!
end function d_nearless

pure function i_nearless (xval) result (i_nl)
!  Nearest value less than given value
! __________________________________________________________
      integer, intent (in) :: xval
      integer :: i_nl
! __________________________________________________________
      i_nl = xval - 1
      return
!
end function i_nearless

end module m_unirnk