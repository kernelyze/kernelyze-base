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
! Module set_precision provides the kind type parameter needed
! to define the precision of a complete package along
! with values for all commonly used precisions
    MODULE set_precision
      use, intrinsic  :: iso_c_binding, only : c_char
! ..
! .. Intrinsic Functions ..
      INTRINSIC KIND
! .. Parameters ..
! Define the standard precisions
! For IEEE standard arithmetic we could also use
!     INTEGER, PARAMETER :: skind = SELECTED_REAL_KIND(p=6, r=37)
!     INTEGER, PARAMETER :: dkind = SELECTED_REAL_KIND(p=15, r=307)
      INTEGER, PARAMETER :: skind = KIND(0.0E0)
      INTEGER, PARAMETER :: dkind = KIND(0.0D0)
! The next statement is required to run the codes
! associated with Chapter 10 on IEEE arithmetic.
! This is non-standard and may not be available.
! Different compilers set this value in different ways;
! see comments below for advice.
!     INTEGER, PARAMETER:: qkind = ...
! Set the precision for the whole package
      INTEGER, PARAMETER :: wp = dkind
! To change the default package precision to single precision change
! the parameter assignment to wp above to
!     INTEGER, PARAMETER :: wp = skind
! and recompile the complete package.

!-----------------------------------------------------------
! For the non-standard quadruple precision:
! IBM and Intel compilers recognize KIND(0.0Q0) 
!     INTEGER, PARAMETER:: qkind = KIND(0.0Q0)
!-----------------------------------------------------------
! If you are using the NAG compiler then you can
! make use of the NAG-supplied f90_kind module
!     USE, INTRINSIC :: f90_kind
!     INTEGER, PARAMETER:: skind = single ! single precision
!     INTEGER, PARAMETER:: dkind = double ! double precision
!     INTEGER, PARAMETER:: qkind = quad   ! quad   precision
! but note that this module is not part of the standard
! and is thus likely to be non-portable.
!-----------------------------------------------------------
!
! TAK: character kind
      integer, parameter  :: chark = c_char

    END MODULE set_precision
