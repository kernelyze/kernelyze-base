! From "Numerical Computing with Modern Fortran"
! by Richard J. Hanson and Tim Hopkins,
! SIAM OT134 (2013).
! Modified by Thomas A. Knox (placed in module, reformatted,
! ragged array functionality made into a converter rather than
! a pass-and-print subroutine, further comments, renames)
!
! c_f_string modified version of approach in at least two posts by IanH
! 
! Modifications Copyright (c) 2016, 2017 by Kernelyze LLC
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
! For the modifications:
! created on: 2016-02-25
! updated on: 2016-02-25
! updated on: 2016-02-26 (removed rename in binding to C for ragged_array,
!             made it an actual conversion routine)
! updated on: 2016-02-27 (continued work on conversion for ragged arrays
!             in three dimensions)
! updated on: 2016-05-12 (added intent(in) c-to-fortran string converter)
! updated on: 2016-05-13 (fixes to c_str_to_f_str)
! updated on: 2016-05-16 (fixed ridiculous bug in c_str_to_f_str)
! updated on: 2016-05-19 (changed c_str_to_f_str to use variable-length
!             Fortran allocatable character variables)
! updated on: 2016-05-23 (better handling for conversion of string
!             arrays)
! updated on: 2017-04-24 (remove ragged array functionality)

module fort_c_interop_mod

implicit none

contains
  
  pure subroutine copy_string(cstring, fstring, clen, ftoc) bind(c)
    ! Use statements
    use, intrinsic :: iso_c_binding, only : c_char, c_int, c_null_char
    ! Arguments
    character(c_char), intent(inout)  :: cstring(*), fstring(*)
    integer(c_int), intent(inout)     :: clen
    integer(c_int), intent(in)        :: ftoc
    ! Body
    ! If ftoc is set to 1 then convert a Fortran string of length clen
    ! stored in fstring into a C string stored in cstring. This contains
    ! the same characters but is terminated by the null character.

    ! If ftoc is set to 0 then convert a C string stored in cstring
    ! into a Fortran string stored in fstring. This contains the same
    ! characters but is NOT terminated by the null character. clen is
    ! set to the clength of the copied string

    if (ftoc == 1) then
      cstring(1:clen) = fstring(1:clen)
      cstring(clen+1) = c_null_char
    else
      clen = 0
      ! The do while construct has been deprecated :-(
      do 
        if (cstring(clen+1) == c_null_char) exit ! for a while!
        clen = clen+1
        fstring(clen) = cstring(clen)
      end do
    end if
  end subroutine copy_string
  
  subroutine c_str_to_f_str(cstring, fstring)
    ! Use statements
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    ! Arguments
    character(c_char), intent(in)   :: cstring(*)
    character(kind=c_char,len=:), &
        allocatable, intent(out)    :: fstring
    ! Local variables
    integer(c_int)  :: clen
    integer(c_int)  :: i
    ! Body
    ! Convert a C string stored in cstring
    ! into a Fortran string stored in fstring. This contains the same
    ! characters but is NOT terminated by the null character. clen is
    ! set to the clength of the copied string
    clen = 0
    do
      if (cstring(clen+1) == c_null_char) exit ! while condition
      clen = clen+1
    end do
    allocate(character(kind=c_char,len=clen)::fstring)
    do i = 1, clen
      fstring(i:i) = cstring(i)
    end do
  end subroutine c_str_to_f_str
  
  ! Written as a modification of IanH posts:
  ! Copy a null terminated C string (specified via a non-null c_ptr) to an 
  ! allocatable deferred length default character variable.
  subroutine c_f_string(c_string, f_string)
    ! Use statement
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_f_pointer
    ! Arguments
    type(c_ptr), intent(in)                                 :: c_string
    character(kind=c_char, len=:), allocatable, intent(out) :: f_string
    ! Local variables
    ! Array for accessing string pointed at by C pointer 
    character(kind=c_char, len=1), pointer :: string_ptr(:)
    ! String index
    integer :: i
    interface
      ! Steal std C library function rather than writing our own.
      function strlen(s) bind(c, name='strlen')
        use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
        implicit none
        ! Argument
        type(c_ptr), intent(in), value :: s
        ! Result
        integer(c_size_t) :: strlen
      end function strlen
    end interface
    ! Body
    ! Map C pointer to fortran character array
    call c_f_pointer(c_string, string_ptr, [strlen(c_string)])
    ! Allocate fortran character variable to the c string's length
    allocate(character(size(string_ptr)) :: f_string)
    ! Copy across (with possible kind conversion) characters
    do i = 1, size(string_ptr)
      f_string(i:i) = string_ptr(i)
    end do
  end subroutine c_f_string
  
  ! This is not ideal; there's a real problem in returning
  ! variable-length strings *from* Fortran *to* C.  If memory
  ! is allocated for them on the Fortran side, as would seem
  ! natural, then it must be persistent -- this is obviously
  ! quite bad (think SAVE or module variables), and leads
  ! almost inevitably to memory leaks (how does Fortran know
  ! when to deallocate?).  The approach below uses the assumption
  ! that the C string has been allocated (with c_str_size) on
  ! the C side, where c_str_size includes the space to be used
  ! by the trailing null character.  A more sophisticated approach
  ! would demand use of a descriptor on the C side a la Fortran
  ! 2015, but using arrays of strings would still be troublesome
  ! because such descriptors would need to be applied to elements
  ! of a structure that was (in general) itself allocatable . . ..
  subroutine f_c_string(f_string, c_str_size, c_string, c_size_too_small)
    ! Use statement
    use, intrinsic :: iso_c_binding, only : &
        c_bool, c_int, c_char, c_null_char, c_ptr, c_f_pointer
    ! Arguments
    character(kind=c_char, len=*), intent(in) :: f_string
    integer(c_int), intent(in)                :: c_str_size
    type(c_ptr), intent(out)                  :: c_string
    logical(c_bool), intent(out)              :: c_size_too_small
    ! Local variables
    ! Array for accessing string pointed at by C pointer 
    character(kind=c_char, len=1), pointer :: string_ptr(:)
    ! String index
    integer :: i
    ! Body
    if (len(f_string) > c_str_size - 1) then
      c_size_too_small = .true.
    else
      c_size_too_small = .false.
    end if
    ! Map C pointer to fortran character array
    call c_f_pointer(c_string, string_ptr, [c_str_size])
    ! Copy over
    do i = 1, min(len(f_string), c_str_size - 1)
      string_ptr(i) = f_string(i:i)
    end do
    ! Remember to null-terminate the C string
    string_ptr(min(len(f_string) + 1, c_str_size)) = c_null_char
  end subroutine f_c_string
  
end module fort_c_interop_mod
