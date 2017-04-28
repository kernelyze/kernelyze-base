! From "Numerical Computing with Modern Fortran"
! by Richard J. Hanson and Tim Hopkins,
! SIAM OT134 (2013).  "SIAM grants a royalty-free license to copy
! and distribute the software code posted on the book's
! supplemental Web page, provided the source is acknowledged."
! (Statement above Library of Congress Cataloging-in-Publication
!  Data.)  The authors also state that (page x of the Introduction):
! "there is no restriction on the use of the software for any
! purpose.  The only requirement is that, if use is made of
! our codes, our book and SIAM are referenced in any derived
! work."  I gratefully acknowledge the book (cited above)
! and SIAM.
!
! This is not a module because it is accessed via
! an "include" statement in the quicksort module.
  
    subroutine exchange(i,j)
! Exchange the contents of the ith and jth elements
    INTEGER, INTENT(IN) :: i,j
    REAL(wp) :: t
    t = a(i)
    a(i) = a(j)
    a(j) = t
    END SUBROUTINE exchange
  
    LOGICAL FUNCTION compare(i,j)
! Compare the contents of the ith and jth elements
! This determines the final sorting order.
! Code for ascending order.
    INTEGER, INTENT(IN) :: i,j
    compare = a(i) < a(j)
    END FUNCTION compare
  
    subroutine compex(i,j)
    INTEGER, INTENT(IN) :: i,j
    IF(compare(j,i)) CALL exchange(i,j)
    END SUBROUTINE compex
  
    subroutine moveValue(i,j)
! Overwrite the contents of the jth element
! with the contents of the ith element
    INTEGER, INTENT(IN) :: i,j
    a(i) = a(j)
    END SUBROUTINE moveValue
  
! The next three subprograms are used to store, 
! compare against and restore a particular element.
    LOGICAL FUNCTION compareValue(j)
    INTEGER, INTENT(IN) :: j
    compareValue = savedVal < a(j)
    END FUNCTION compareValue
  
    subroutine saveValue(i)
    INTEGER, INTENT(IN) :: i
    savedVal = a(i)
    END SUBROUTINE saveValue
  
    subroutine restoreValue(i)
    INTEGER, INTENT(IN) :: i
    a(i) = savedVal
    END SUBROUTINE restoreValue
