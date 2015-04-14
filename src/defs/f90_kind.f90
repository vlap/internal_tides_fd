! Copyright (c) 1994 Unicomp, Inc.
!
! Developed at Unicomp, Inc.
!
! Permission to use, copy, modify, and distribute this
! software is freely granted, provided that this notice 
! is preserved.

! NAGWare F90

module f90_kind

   implicit none

   intrinsic kind
   private kind

!  For each intrinsic data type, the number of kinds and the kind numbers:

   integer, parameter :: NUMBER_OF_INTEGER_KINDS = 3
   integer, parameter :: INTEGER_KINDS (NUMBER_OF_INTEGER_KINDS) = &
                               (/ 1, 2, 3 /)

   integer, parameter :: NUMBER_OF_LOGICAL_KINDS = 2
   integer, parameter :: LOGICAL_KINDS (NUMBER_OF_LOGICAL_KINDS) = &
                               (/ 1, 2 /)

   integer, parameter :: NUMBER_OF_REAL_KINDS = 2
   integer, parameter :: REAL_KINDS (NUMBER_OF_REAL_KINDS) = &
                               (/ 1, 2 /)

   integer, parameter :: NUMBER_OF_CHARACTER_KINDS = 1
   integer, parameter :: CHARACTER_KINDS (NUMBER_OF_CHARACTER_KINDS) = &
                               (/ 1 /)

!  The default kinds

   integer, parameter :: DEFAULT_INTEGER_KIND = kind (0)
   integer, parameter :: DEFAULT_LOGICAL_KIND = kind (.false.)
   integer, parameter :: single = kind (0.0)
   integer, parameter :: double = kind (0.0d0)
   integer, parameter :: character = kind ("A")


end module f90_kind
