! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
MODULE  Extended  ! precision specification for real computations

    ! This requests the processor to use a real implementation 'stnd'
    ! which provides at least 20 decimal digits of precision and an
    ! exponent range of at least 10 ^ +- 80.  It is expected that this
    ! precision may not be available on all machines.
    ! In July 2002, we found this available on
    !     SUN (Solaris) with f95
    !     IBM (AIX)     with xlf90
    !     DEC alpha     with f90

    IMPLICIT NONE
    Integer, PUBLIC, PARAMETER :: stnd = Selected_Real_Kind ( 20, 80 )
    !-------------------------

    ! A few computations are preferably done in higher precision 'extd'. The
    ! numbers chosen here should be such that the underlying hardware will
    ! select a higher precision for kind 'extd' than for kind 'stnd', if
    ! this is feasible.  If a higher precision is not readily available,
    ! the same values may be used as are given above for 'stnd'. It is
    ! anticipated that on many machines, such an even higher precision may
    ! not be available.

    Integer, PUBLIC, PARAMETER :: extd = Selected_Real_Kind ( 30, 80 )
    !-------------------------

end Module Extended
MODULE  Integers  ! precision specification for integer computations

    ! This is provided for those machines where using short integers
    ! has some advantage. The range here is from -999 to +999. Note that
    ! 999 = 10^3 -1. No harm will be done if short integers are made
    ! the same as long integers.

    IMPLICIT NONE
    Integer, PUBLIC, PARAMETER :: short = Selected_Int_Kind ( 3 )
    !-------------------------

    ! The range here is at least from  -9 999 999  to +9 999 999. Note
    ! that 10^7 - 1 = 9,999,999. This may limit the largest possible value
    ! of the dimension  n  of problems which can be solved. If n is to
    ! be larger, the 7 should be replaced with a number k so that
    ! n is considerably less than 10^k.

    Integer, PUBLIC, PARAMETER :: long  = Selected_Int_Kind ( 7 )
    !-------------------------

end Module Integers
MODULE  Low  ! precision specification for real computations

    ! This requests the processor to use a real implementation 'stnd'
    ! which provides at least 6 decimal digits of precision and an
    ! exponent range of at least 10 ^ +- 35.  This would be suitable for
    ! low accuracy computations.  It is expected that this precision will
    ! be available on all machines.

    IMPLICIT NONE
    Integer, PUBLIC, PARAMETER :: stnd = Selected_Real_Kind ( 6, 35 )
    !-------------------------

    ! A few computations are preferably done in higher precision 'extd'. The
    ! numbers chosen here should be such that the underlying hardware will
    ! select a higher precision for kind 'extd' than for kind 'stnd', if
    ! this is feasible.  If a higher precision is not readily available,
    ! the same values may be used as are given above for 'stnd'. It is
    ! anticipated that on most machines this higher precision will also
    ! be available.

    Integer, PUBLIC, PARAMETER :: extd = Selected_Real_Kind ( 12, 35 )
    !-------------------------

end Module Low
MODULE  Normal  ! precision specification for real computations

    ! This requests the processor to use a real implementation 'stnd'
    ! which provides at least 12 decimal digits of precision and an
    ! exponent range of at least 10 ^ +- 50.  It is expected that this
    ! precision will be available on all machines.

    IMPLICIT NONE
    Integer, PUBLIC, PARAMETER :: stnd = Selected_Real_Kind ( 12, 50 )
    !-------------------------

    ! A few computations are preferably done in higher precision 'extd'. The
    ! numbers chosen here should be such that the underlying hardware will
    ! select a higher precision for kind 'extd' than for kind 'stnd', if
    ! this is feasible.  If a higher precision is not readily available,
    ! the same values may be used as are given above for 'stnd'. It is
    ! anticipated that on many machines this higher precision may
    ! not be available.

   !Integer, PUBLIC, PARAMETER :: extd = Selected_Real_Kind ( 20, 50 ) ! preferred
    Integer, PUBLIC, PARAMETER :: extd = Selected_Real_Kind ( 12, 50 ) ! NAG f95
    !-------------------------

end Module Normal
MODULE Precision_Model

    ! This provides a convenient way of selecting the precision
    ! required for a computation. By simply ensuring that a leading '!'
    ! appears on all but exactly one of the following USE statements,
    ! and then recompiling all routines, the precision of an entire
    ! computation can be altered.

    !    USE Low
         USE Normal
    !    USE Extended

         USE Integers

   ! This is the original F90 code
     !   PRIVATE
     !   PUBLIC :: stnd, extd, short, long
   ! This is the F code
         PUBLIC

end Module Precision_Model
