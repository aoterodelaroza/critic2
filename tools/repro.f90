! Minimal reproducer for a gfortran code-generation bug that crashes the
! critic2 GUI on Windows. See tools/check-gfortran-bug.md for context.
!
! The bug: passing an array section of a POLYMORPHIC object's allocatable
! component (md%r(:,i)) as the argument of a TYPE-BOUND dispatch call
! (c%c2x(...)) makes the compiler emit a bogus array descriptor for the
! component, so the access faults (spurious out-of-bounds / segfault). This is
! the exact shape of dynamics@proc.f90:157 (keep_in_cell).
!
! Build and run with the compiler you intend to use:
!   gfortran -g -O0 -fcheck=all repro.f90 -o repro && ./repro          # native
!   x86_64-w64-mingw32-gfortran -g -O0 -fcheck=all repro.f90 -o repro.exe   # cross
!
! An AFFECTED compiler aborts with:
!   Fortran runtime error: Index '1' of dimension 2 of array 'md%r'
!   outside of expected range (0:0)
! A GOOD compiler prints: COMPILER OK

module m
  implicit none
  type t
     real*8, allocatable :: r(:,:)   ! like mdrun%r (positions)
     integer :: nat = 0
  end type t
  type cr
     real*8 :: m(3,3) = reshape((/2d0,0d0,0d0, 0d0,2d0,0d0, 0d0,0d0,2d0/),(/3,3/))
   contains
     procedure :: c2x => cr_c2x      ! type-bound, like crystal%c2x
  end type cr
contains
  function cr_c2x(this,x) result(y)
    class(cr), intent(in) :: this
    real*8, intent(in) :: x(3)
    real*8 :: y(3)
    y = matmul(this%m,x)
  end function cr_c2x

  subroutine trigger(md,c)
    class(t),  intent(inout) :: md   ! polymorphic, like class(mdrun)
    class(cr), intent(inout) :: c    ! polymorphic, like class(crystal)
    integer :: i
    real*8 :: xf(3)
    do i = 1, md%nat
       xf = c%c2x(md%r(:,i))         ! <-- the miscompiled pattern
       if (any(abs(xf - 2d0) > 1d-9)) stop "wrong result"
    end do
  end subroutine trigger
end module m

program repro
  use m
  implicit none
  type(t)  :: a
  type(cr) :: c
  a%nat = 2
  allocate(a%r(3,2))
  a%r = 1d0
  call trigger(a,c)                  ! aborts here on an affected compiler
  print '(a)', "COMPILER OK (not affected by the polymorphic-dispatch codegen bug)"
end program repro
