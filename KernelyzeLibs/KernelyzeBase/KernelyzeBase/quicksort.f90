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
  
MODULE quicksort
USE set_precision, ONLY : wp

CONTAINS

  subroutine qsort(a, left, right, savedVal)
  INTEGER, INTENT(IN) :: left, right
  REAL(wp), INTENT(INOUT) :: a(left:right)
  REAL(wp), INTENT(INOUT) :: savedVal

  CALL quicksort(left, right)
  CALL insertion

  CONTAINS
    
    subroutine insertion
    INTEGER :: i, j

    DO i = left+1, right
      CALL compex(left, i)
    END DO

    DO i = left+2, right
      j = i
      CALL saveValue(i)
      DO WHILE(compareValue(j-1))
        CALL moveValue(j, j-1)
        j = j-1
      END DO

      CALL restoreValue(j)
    END DO

    END SUBROUTINE insertion

    RECURSIVE SUBROUTINE quicksort(left, right)
    INTEGER, INTENT(IN) :: left, right

    INTEGER, PARAMETER:: switchsorts=10
    INTEGER :: i

    IF((right-left) > switchsorts) THEN
      CALL exchange((right+left)/2, (right-1))
      CALL compex(left, right-1)
      CALL compex(right, left)
      CALL compex(right-1, right)
      i = partition(left+1, right-1)
      CALL quicksort(left, i-1)
      CALL quicksort(i+1, right)
    END IF

    END SUBROUTINE quicksort

    function partition(left, right) RESULT(i)
    INTEGER, INTENT(IN) :: left, right
    INTEGER :: i, j

    i = left - 1
    j = right
    CALL saveValue(right)

    DO
      DO
        i = i+1
        IF(i>right) EXIT
        IF(compareValue(i)) EXIT
      END DO

      DO
        j = j-1
        IF(.NOT.compareValue(j) .OR. j==left) EXIT
      END DO

      IF(i>= j) EXIT
      CALL exchange(i,j)

    END DO

    CALL exchange(i, right)

    END FUNCTION partition

    INCLUDE 'qsfuns.f90'
  
  END SUBROUTINE qsort

END MODULE quicksort
