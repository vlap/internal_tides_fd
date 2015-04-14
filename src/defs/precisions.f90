    MODULE precisions

!      MODULE of precision parameters to be USEd by free-format source.

!      WP = 'working precision'
!      HP = 'higher precision'
!      dp = 'double precision'
!      zp = 'double precision' (for complex)
!      I4B = 4-byte INTEGER

!      Example: for a single-precision implementation
!        RP = KIND(0.0E0)
!        WP = KIND(0.0E0)
!        HP = KIND(0.0D0)
!        I4B = SELECTED_INT_KIND(9)

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..

       integer, parameter              :: si = selected_int_kind(8) ! ShortInteger, equivalent to integer(kind=4), max integer = 2147483647, can index up to 2Gb
       integer, parameter              :: li = selected_int_kind(10) ! LongInteger, equivalent to integer(kind=8), max integer = 9223372036854775807, can index ALOT

       integer, parameter              :: sp = kind(1.0)
       integer, parameter              :: dp = kind(1.d0)
       integer, parameter              :: qp = selected_real_kind(2*precision(1.0_dp))

!------------- GRID, BASIC VECTORS AND MATRICES ARE GENERATED AND SAVED IN DOUBLE PRECISION ------------------!
!-------------- THIS ALLOWES TO LINK TO PREVIOUSLY GENERATED GRIDS REGARDLESS OF PRECISION -------------------!
       integer, parameter              :: wp = dp
       integer, parameter              :: cwp = dp


       real(sp), parameter             :: eps = 1e-10                   ! epsilon for rounding error

!      .. Intrinsic Functions ..
       INTRINSIC                          KIND, SELECTED_INT_KIND

    END MODULE precisions
