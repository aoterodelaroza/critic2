! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>. 
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Keast library: cubatures of a tetrahedron
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!  Modified:
!    02 July 2008
!  Author:
!    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!  Downloaded from John Burkardt's page, Florida State University:
!    http://people.sc.fsu.edu/~jburkardt/f_src/keast/keast.html

! module keast
! COMP_NEXT computes the compositions of the integer N into K parts.
! I4_MODP returns the nonnegative remainder of I4 division.
! I4_WRAP forces an I4 to lie between given limits by wrapping.
! KEAST_DEGREE returns the degree of a Keast rule for the tetrahedron.
! KEAST_ORDER_NUM returns the order of a Keast rule for the tetrahedron.
! KEAST_RULE returns the points and weights of a Keast rule.
! KEAST_RULE_NUM returns the number of Keast rules for the tetrahedron.
! KEAST_SUBORDER returns the suborders for a Keast rule.
! KEAST_SUBORDER_NUM returns the number of suborders for a Keast rule.
! KEAST_SUBRULE returns a compressed Keast rule.
! MONOMIAL_VALUE evaluates a monomial.
! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
! TETRAHEDRON_REFERENCE_TO_PHYSICAL maps T4 reference points to physical points.
! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
! TIMESTAMP prints the current YMDHMS date as a time stamp.
module keast
  implicit none

  private
  public :: keast_order_num
  public :: keast_rule
  public :: keast_rule_num

contains
  
  subroutine comp_next ( n, k, a, more, h, t )

    !*****************************************************************************80
    !
    !! COMP_NEXT computes the compositions of the integer N into K parts.
    !
    !  Discussion:
    !
    !    A composition of the integer N into K parts is an ordered sequence
    !    of K nonnegative integers which sum to N.  The compositions (1,2,1)
    !    and (1,1,2) are considered to be distinct.
    !
    !    The routine computes one composition on each call until there are no more.
    !    For instance, one composition of 6 into 3 parts is
    !    3+2+1, another would be 6+0+0.
    !
    !    On the first call to this routine, set MORE = FALSE.  The routine
    !    will compute the first element in the sequence of compositions, and
    !    return it, as well as setting MORE = TRUE.  If more compositions
    !    are desired, call again, and again.  Each time, the routine will
    !    return with a new composition.
    !
    !    However, when the LAST composition in the sequence is computed 
    !    and returned, the routine will reset MORE to FALSE, signaling that
    !    the end of the sequence has been reached.
    !
    !    This routine originally used a SAVE statement to maintain the
    !    variables H and T.  I have decided that it is safer
    !    to pass these variables as arguments, even though the user should
    !    never alter them.  This allows this routine to safely shuffle
    !    between several ongoing calculations.
    !
    !
    !    There are 28 compositions of 6 into three parts.  This routine will
    !    produce those compositions in the following order:
    !
    !     I         A
    !     -     ---------
    !     1     6   0   0
    !     2     5   1   0
    !     3     4   2   0
    !     4     3   3   0
    !     5     2   4   0
    !     6     1   5   0
    !     7     0   6   0
    !     8     5   0   1
    !     9     4   1   1
    !    10     3   2   1
    !    11     2   3   1
    !    12     1   4   1
    !    13     0   5   1
    !    14     4   0   2
    !    15     3   1   2
    !    16     2   2   2
    !    17     1   3   2
    !    18     0   4   2
    !    19     3   0   3
    !    20     2   1   3
    !    21     1   2   3
    !    22     0   3   3
    !    23     2   0   4
    !    24     1   1   4
    !    25     0   2   4
    !    26     1   0   5
    !    27     0   1   5
    !    28     0   0   6
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    02 July 2008
    !
    !  Author:
    !
    !    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    Albert Nijenhuis, Herbert Wilf,
    !    Combinatorial Algorithms for Computers and Calculators,
    !    Second Edition,
    !    Academic Press, 1978,
    !    ISBN: 0-12-519260-6,
    !    LC: QA164.N54.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
    !
    !    Input, integer ( kind = 4 ) K, the number of parts in the composition.
    !
    !    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
    !
    !    Input/output, logical MORE, set by the user to start the computation,
    !    and by the routine to terminate it.
    !
    !    Input/output, integer ( kind = 4 )  H, T, two internal parameters needed for the
    !    computation.  The user should allocate space for these in the calling
    !    program, include them in the calling sequence, but never alter them!
    !
    implicit none

    integer ( kind = 4 ) k

    integer ( kind = 4 ) a(k)
    integer ( kind = 4 ) h
    logical more
    integer ( kind = 4 ) n
    integer ( kind = 4 ) t
    !
    !  The first computation.
    !
    if ( .not. more ) then

       t = n
       h = 0
       a(1) = n
       a(2:k) = 0
       !
       !  The next computation.
       !
    else

       if ( 1 < t ) then
          h = 0
       end if

       h = h + 1
       t = a(h)
       a(h) = 0
       a(1) = t - 1
       a(h+1) = a(h+1) + 1

    end if
    !
    !  This is the last element of the sequence if all the
    !  items are in the last slot.
    !
    more = ( a(k) /= n )

    return
  end subroutine comp_next
  subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    logical ( kind = 4 ) lopen

    iunit = 0

    do i = 1, 99

       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          inquire ( unit = i, opened = lopen, iostat = ios )

          if ( ios == 0 ) then
             if ( .not. lopen ) then
                iunit = i
                return
             end if
          end if

       end if

    end do

    return
  end subroutine get_unit
  function i4_modp ( i, j )

    !*****************************************************************************80
    !
    !! I4_MODP returns the nonnegative remainder of I4 division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !  Example:
    !
    !        I     J     MOD I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number to be divided.
    !
    !    Input, integer ( kind = 4 ) J, the number that divides I.
    !
    !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) j
    integer ( kind = 4 ) value

    if ( j == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'I4_MODP - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
       stop
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
       value = value + abs ( j )
    end if

    i4_modp = value

    return
  end function i4_modp
  function i4_wrap ( ival, ilo, ihi )

    !*****************************************************************************80
    !
    !! I4_WRAP forces an I4 to lie between given limits by wrapping.
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  Value
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) IVAL, an integer value.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
    !
    !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
    !
    implicit none

    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) ival
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    integer ( kind = 4 ) value
    integer ( kind = 4 ) wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
       value = jlo
    else
       value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value

    return
  end function i4_wrap
  subroutine keast_degree ( rule, degree )

    !*****************************************************************************80
    !
    !! KEAST_DEGREE returns the degree of a Keast rule for the tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 June 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) DEGREE, the degree of the rule.
    !
    implicit none

    integer ( kind = 4 ) degree
    integer ( kind = 4 ) rule

    if ( rule == 1 ) then
       degree = 0
    else if ( rule == 2 ) then
       degree = 1
    else if ( rule == 3 ) then
       degree = 2
    else if ( rule == 4 ) then
       degree = 3 
    else if ( rule == 5 ) then
       degree = 4
    else if ( rule == 6 ) then
       degree = 4 
    else if ( rule == 7 ) then
       degree = 5 
    else if ( rule == 8 ) then
       degree = 6 
    else if ( rule == 9 ) then
       degree = 7 
    else if ( rule == 10 ) then
       degree = 8
    else
       degree = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KEAST_DEGREE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
       stop
    end if

    return
  end subroutine keast_degree
  subroutine keast_order_num ( rule, order_num )

    !*****************************************************************************80
    !
    !! KEAST_ORDER_NUM returns the order of a Keast rule for the tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 December 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) ORDER_NUM, the order of the rule.
    !
    implicit none

    integer ( kind = 4 ) order_num
    integer ( kind = 4 ) rule
    integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
    integer ( kind = 4 ) suborder_num

    call keast_suborder_num ( rule, suborder_num )

    allocate ( suborder(1:suborder_num) )

    call keast_suborder ( rule, suborder_num, suborder )

    order_num = sum ( suborder(1:suborder_num) )

    deallocate ( suborder )

    return
  end subroutine keast_order_num
  subroutine keast_rule ( rule, order_num, xyz, w )

    !*****************************************************************************80
    !
    !! KEAST_RULE returns the points and weights of a Keast rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 June 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) ORDER_NUM, the order (number of points) of the rule.
    !
    !    Output, real ( kind = 8 ) XYZ(3,ORDER_NUM), the points of the rule.
    !
    !    Output, real ( kind = 8 ) W(ORDER_NUM), the weights of the rule.
    !
    implicit none

    integer ( kind = 4 ) order_num

    integer ( kind = 4 ) k
    integer ( kind = 4 ) o
    integer ( kind = 4 ) rule
    integer ( kind = 4 ) s
    integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
    integer ( kind = 4 ) suborder_num
    real    ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
    real    ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyzz
    real    ( kind = 8 ) w(order_num)
    real    ( kind = 8 ) xyz(3,order_num)
    !
    !  Get the suborder information.
    !
    call keast_suborder_num ( rule, suborder_num )

    allocate ( suborder(suborder_num) )
    allocate ( suborder_xyzz(4,suborder_num) )
    allocate ( suborder_w(suborder_num) )

    call keast_suborder ( rule, suborder_num, suborder )

    call keast_subrule ( rule, suborder_num, suborder_xyzz, suborder_w )
    !
    !  Expand the suborder information to a full order rule.
    !
    o = 0

    do s = 1, suborder_num

       if ( suborder(s) == 1 ) then

          o = o + 1
          xyz(1:3,o) = suborder_xyzz(1:3,s)
          w(o) = suborder_w(s)
          !
          !  For SUBORDER = 4, we list the coordinates of the generator as
          !
          !    A,B,B,B
          !
          !  and we generate
          !
          !    A, B, B = (1,2,3)
          !    B, B, B = (2,3,4)
          !    B, B, A = (3,4,1)
          !    B, A, B = (4,1,2)
          !
       else if ( suborder(s) == 4 ) then

          do k = 1, 4
             o = o + 1
             xyz(1,o) = suborder_xyzz ( i4_wrap(k,  1,4), s )
             xyz(2,o) = suborder_xyzz ( i4_wrap(k+1,1,4), s )
             xyz(3,o) = suborder_xyzz ( i4_wrap(k+2,1,4), s )
             w(o) = suborder_w(s)
          end do
          !
          !  For SUBORDER = 6, we list the coordinates of the generator as
          !
          !    A,A,B,B
          !
          !  and we generate
          !
          !    B, A, A = (4,1,2)
          !    A, B, A = (1,4,2)
          !    A, A, B = (1,2,4)
          !
          !    A, B, B = (1,3,4)
          !    B, A, B = (4,2,3)
          !    B, B, A = (4,3,1)
          !
       else if ( suborder(s) == 6 ) then

          do k = 1, 3
             o = o + 1
             xyz(1:3,o) = suborder_xyzz ( 1, s )
             xyz(k,o)   = suborder_xyzz ( 3, s )
             w(o) = suborder_w(s)
          end do

          do k = 1, 3
             o = o + 1
             xyz(1:3,o) = suborder_xyzz ( 3, s )
             xyz(k,o)   = suborder_xyzz ( 1, s )
             w(o) = suborder_w(s)
          end do
          !
          !  For SUBORDER = 12, we list the coordinates of the generator as
          !
          !    A,A,B,C
          !
          !  and we generate
          !
          !    B, A, A
          !    A, B, A
          !    A, A, B
          !
          !    C, A, A
          !    A, C, A
          !    A, A, C
          !   
          !    A, B, C
          !    B, C, A
          !    C, A, B
          !    A, C, B
          !    C, B, A
          !    B, A, C
          !
       else if ( suborder(s) == 12 ) then

          do k = 1, 3
             o = o + 1
             xyz(1:3,o) = suborder_xyzz ( 1, s )
             xyz(k,o)   = suborder_xyzz ( 3, s )
             w(o) = suborder_w(s)
          end do

          do k = 1, 3
             o = o + 1
             xyz(1:3,o) = suborder_xyzz ( 1, s )
             xyz(k,o)   = suborder_xyzz ( 4, s )
             w(o) = suborder_w(s)
          end do

          do k = 1, 3
             o = o + 1
             xyz(1,o) = suborder_xyzz ( i4_wrap(k+1,2,4), s )
             xyz(2,o) = suborder_xyzz ( i4_wrap(k+2,2,4), s )
             xyz(3,o) = suborder_xyzz ( i4_wrap(k+3,2,4), s )
             w(o) = suborder_w(s)
          end do

          do k = 1, 3
             o = o + 1
             xyz(1,o) = suborder_xyzz ( i4_wrap(k+1,2,4), s )
             xyz(2,o) = suborder_xyzz ( i4_wrap(k+3,2,4), s )
             xyz(3,o) = suborder_xyzz ( i4_wrap(k+2,2,4), s )
             w(o) = suborder_w(s)
          end do


       else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'KEAST_RULE - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s) 
          write ( *, '(a,i8)' ) '  RULE =    ', rule
          write ( *, '(a,i8)' ) '  ORDER_NUM = ', order_num
          stop

       end if

    end do

    deallocate ( suborder )
    deallocate ( suborder_xyzz )
    deallocate ( suborder_w )

    return
  end subroutine keast_rule
  subroutine keast_rule_num ( rule_num )

    !*****************************************************************************80
    !
    !! KEAST_RULE_NUM returns the number of Keast rules for the tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 December 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) RULE_NUM, the number of rules.
    !
    implicit none

    integer ( kind = 4 ) rule_num

    rule_num = 10

    return
  end subroutine keast_rule_num
  subroutine keast_suborder ( rule, suborder_num, suborder )

    !*****************************************************************************80
    !
    !! KEAST_SUBORDER returns the suborders for a Keast rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 June 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders.
    !
    !    Output, integer ( kind = 4 ) SUBORDER(SUBORDER_NUM), the suborders.
    !
    implicit none

    integer ( kind = 4 ) suborder_num

    integer ( kind = 4 ) rule
    integer ( kind = 4 ) suborder(suborder_num)

    if ( rule == 1 ) then
       suborder(1:suborder_num) = (/ 1 /)
    else if ( rule == 2 ) then
       suborder(1:suborder_num) = (/ 4 /)
    else if ( rule == 3 ) then
       suborder(1:suborder_num) = (/ 1, 4 /)
    else if ( rule == 4 ) then
       suborder(1:suborder_num) = (/ 4, 6 /)
    else if ( rule == 5 ) then
       suborder(1:suborder_num) = (/ 1, 4, 6 /)
    else if ( rule == 6 ) then
       suborder(1:suborder_num) = (/ 6, 4, 4 /)
    else if ( rule == 7 ) then
       suborder(1:suborder_num) = (/ 1, 4, 4, 6 /)
    else if ( rule == 8 ) then
       suborder(1:suborder_num) = (/ 4, 4, 4, 12 /)
    else if ( rule == 9 ) then
       suborder(1:suborder_num) = (/ 1, 4, 4, 4, 6, 12 /)
    else if ( rule == 10 ) then
       suborder(1:suborder_num) = (/ 1, 4, 4, 6, 6, 12, 12 /)
    else
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KEAST_SUBORDER - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
       stop
    end if

    return
  end subroutine keast_suborder
  subroutine keast_suborder_num ( rule, suborder_num )

    !*****************************************************************************80
    !
    !! KEAST_SUBORDER_NUM returns the number of suborders for a Keast rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 June 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders.
    !
    implicit none

    integer ( kind = 4 ) rule
    integer ( kind = 4 ) suborder_num

    if ( rule == 1 ) then
       suborder_num = 1
    else if ( rule == 2 ) then
       suborder_num = 1
    else if ( rule == 3 ) then
       suborder_num = 2
    else if ( rule == 4 ) then
       suborder_num = 2
    else if ( rule == 5 ) then
       suborder_num = 3
    else if ( rule == 6 ) then
       suborder_num = 3
    else if ( rule == 7 ) then
       suborder_num = 4
    else if ( rule == 8 ) then
       suborder_num = 4
    else if ( rule == 9 ) then
       suborder_num = 6
    else if ( rule == 10 ) then
       suborder_num = 7
    else
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KEAST_SUBORDER_NUM - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
       stop
    end if

    return
  end subroutine keast_suborder_num
  subroutine keast_subrule ( rule, suborder_num, suborder_xyzz, suborder_w )

    !*****************************************************************************80
    !
    !! KEAST_SUBRULE returns a compressed Keast rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 June 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Patrick Keast,
    !    Moderate Degree Tetrahedral Quadrature Formulas,
    !    Computer Methods in Applied Mechanics and Engineering,
    !    Volume 55, Number 3, May 1986, pages 339-348.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders 
    !    of the rule.
    !
    !    Output, real ( kind = 8 ) SUBORDER_XYZZ(4,SUBORDER_NUM),
    !    the barycentric coordinates of the abscissas.
    !
    !    Output, real ( kind = 8 ) SUBORDER_W(SUBORDER_NUM), the
    !    suborder weights.
    !
    implicit none

    integer ( kind = 4 ) suborder_num

    integer ( kind = 4 ) :: i4_4 = 4
    integer ( kind = 4 ) rule
    real    ( kind = 8 ) suborder_w(suborder_num)
    real    ( kind = 8 ) suborder_xyzz(4,suborder_num)

    if ( rule == 1 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,  0.250000000000000000D+00, &
          0.250000000000000000D+00,  0.250000000000000000D+00  &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.166666666666666667D+00 &
          /)

    else if ( rule == 2 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.585410196624968500D+00,  0.138196601125010500D+00, &
          0.138196601125010500D+00,  0.138196601125010500D+00  &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.0416666666666666667D+00 &
          /)

    else if ( rule == 3 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,  0.250000000000000000D+00, &
          0.250000000000000000D+00,  0.250000000000000000D+00, &
          0.500000000000000000D+00,  0.166666666666666667D+00, &
          0.166666666666666667D+00,  0.166666666666666667D+00 &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          -0.133333333333333333D+00, &
          0.075000000000000000D+00 &
          /)

    else if ( rule == 4 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.568430584196844400D+00,  0.143856471934385200D+00, &
          0.143856471934385200D+00,  0.143856471934385200D+00, &
          0.500000000000000000D+00,  0.500000000000000000D+00, &
          0.000000000000000000D+00,  0.000000000000000000D+00  &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.0362941783134009000D+00,  &
          0.00358165890217718333D+00  &
          /)

    else if ( rule == 5 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,   0.250000000000000000D+00,  &
          0.250000000000000000D+00,   0.250000000000000000D+00,  &
          0.785714285714285714D+00,   0.0714285714285714285D+00, &
          0.0714285714285714285D+00,  0.0714285714285714285D+00, &
          0.399403576166799219D+00,   0.399403576166799219D+00,  &
          0.100596423833200785D+00,   0.100596423833200785D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          -0.0131555555555555556D+00,  &
          0.00762222222222222222D+00, &
          0.0248888888888888889D+00   &
          /)

    else if ( rule == 6 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.500000000000000000D+00,   0.500000000000000000D+00,  &
          0.000000000000000000D+00,   0.000000000000000000D+00, &
          0.698419704324386603D+00,   0.100526765225204467D+00,  &
          0.100526765225204467D+00,   0.100526765225204467D+00, &
          0.0568813795204234229D+00,  0.314372873493192195D+00,  &
          0.314372873493192195D+00,   0.314372873493192195D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.00317460317460317450D+00,  &
          0.0147649707904967828D+00, &
          0.0221397911142651221D+00 &
          /)

    else if ( rule == 7 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,  0.250000000000000000D+00,  &
          0.250000000000000000D+00,  0.250000000000000000D+00,  &
          0.00000000000000000D+00,   0.333333333333333333D+00,  &
          0.333333333333333333D+00,  0.333333333333333333D+00,  &
          0.727272727272727273D+00,  0.0909090909090909091D+00, &
          0.0909090909090909091D+00, 0.0909090909090909091D+00, &
          0.0665501535736642813D+00, 0.0665501535736642813D+00, &
          0.433449846426335728D+00,  0.433449846426335728D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.0302836780970891856D+00,  &
          0.00602678571428571597D+00, &
          0.0116452490860289742D+00,  &
          0.0109491415613864534D+00   &
          /)

    else if ( rule == 8 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.356191386222544953D+00,  0.214602871259151684D+00,  &
          0.214602871259151684D+00,  0.214602871259151684D+00,  &
          0.877978124396165982D+00,  0.0406739585346113397D+00, &
          0.0406739585346113397D+00, 0.0406739585346113397D+00, &
          0.0329863295731730594D+00, 0.322337890142275646D+00,  &
          0.322337890142275646D+00,  0.322337890142275646D+00,  &
          0.0636610018750175299D+00, 0.0636610018750175299D+00, &
          0.269672331458315867D+00,  0.603005664791649076D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.00665379170969464506D+00, &
          0.00167953517588677620D+00, &
          0.00922619692394239843D+00, &
          0.00803571428571428248D+00  &
          /)

    else if ( rule == 9 ) then 

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,    0.250000000000000000D+00,  &
          0.250000000000000000D+00,    0.250000000000000000D+00,  &
          0.765360423009044044D+00,    0.0782131923303186549D+00, &
          0.0782131923303186549D+00,   0.0782131923303186549D+00, &
          0.634470350008286765D+00,    0.121843216663904411D+00,  &
          0.121843216663904411D+00,    0.121843216663904411D+00,  &
          0.00238250666073834549D+00,  0.332539164446420554D+00,  &
          0.332539164446420554D+00,    0.332539164446420554D+00,  &
          0.500000000000000000D+00,    0.500000000000000000D+00,  &
          0.00000000000000000D+00,     0.00000000000000000D+00,   &
          0.100000000000000000D+00,    0.100000000000000000D+00,  &
          0.200000000000000000D+00,    0.600000000000000000D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          0.0182642234661087939D+00,   &
          0.0105999415244141609D+00,   &
          -0.0625177401143299494D+00,   &
          0.00489142526307353653D+00,  &
          0.000970017636684296702D+00, &
          0.0275573192239850917D+00    &
          /)

    else if ( rule == 10 ) then

       suborder_xyzz(1:4,1:suborder_num) = reshape ( (/ &
          0.250000000000000000D+00,  0.250000000000000000D+00,  &
          0.250000000000000000D+00,  0.250000000000000000D+00,  &
          0.617587190300082967D+00,  0.127470936566639015D+00,  &
          0.127470936566639015D+00,  0.127470936566639015D+00,  &
          0.903763508822103123D+00,  0.0320788303926322960D+00, &
          0.0320788303926322960D+00, 0.0320788303926322960D+00, &
          0.0497770956432810185D+00, 0.0497770956432810185D+00, &
          0.450222904356718978D+00,  0.450222904356718978D+00,  &
          0.183730447398549945D+00,  0.183730447398549945D+00,  &
          0.316269552601450060D+00,  0.316269552601450060D+00,  &
          0.231901089397150906D+00,  0.231901089397150906D+00,  &
          0.0229177878448171174D+00, 0.513280033360881072D+00,  &
          0.0379700484718286102D+00, 0.0379700484718286102D+00, &
          0.730313427807538396D+00,  0.193746475248804382D+00   &
          /), (/ i4_4, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
          -0.0393270066412926145D+00,   &
          0.00408131605934270525D+00,  &
          0.000658086773304341943D+00, &
          0.00438425882512284693D+00,  &
          0.0138300638425098166D+00,   &
          0.00424043742468372453D+00,  &
          0.00223873973961420164D+00   &
          /)

    else

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KEAST_SUBRULE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if
    !
    !  Renormalize the weights so that they sum to 1.
    !
    suborder_w(1:suborder_num) = 6.0D+00 * suborder_w(1:suborder_num)

    return
  end subroutine keast_subrule
  subroutine monomial_value ( dim_num, point_num, x, expon, value )

    !*****************************************************************************80
    !
    !! MONOMIAL_VALUE evaluates a monomial.
    !
    !  Discussion:
    !
    !    This routine evaluates a monomial of the form
    !
    !      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
    !
    !    where the exponents are nonnegative integers.  Note that
    !    if the combination 0^0 is encountered, it should be treated
    !    as 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 May 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of points at which the
    !    monomial is to be evaluated.
    !
    !    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the point coordinates.
    !
    !    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
    !
    !    Output, real ( kind = 8 ) VALUE(POINT_NUM), the value of the monomial.
    !
    implicit none

    integer ( kind = 4 ) dim_num
    integer ( kind = 4 ) point_num

    integer ( kind = 4 ) dim
    integer ( kind = 4 ) expon(dim_num)
    real    ( kind = 8 ) value(point_num)
    real    ( kind = 8 ) x(dim_num,point_num)

    value(1:point_num) = 1.0D+00

    do dim = 1, dim_num
       if ( 0 /= expon(dim) ) then
          value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**expon(dim)
       end if
    end do

    return
  end subroutine monomial_value
  function r8mat_det_4d ( a )

    !*****************************************************************************80
    !
    !! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is a two dimensional matrix of double precision real values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
    !
    !    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
    !
    implicit none

    real    ( kind = 8 ) a(4,4)
    real    ( kind = 8 ) r8mat_det_4d

    r8mat_det_4d = &
       a(1,1) * ( &
       a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
       - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
       + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
       a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
       - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
       + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
       a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
       - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
       + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
       a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
       - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
       + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

    return
  end function r8mat_det_4d
  subroutine tetrahedron_reference_to_physical ( t, n, ref, phy )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_REFERENCE_TO_PHYSICAL maps T4 reference points to physical points.
    !
    !  Discussion:
    !
    !    Given the vertices of an order 4 physical tetrahedron and a point
    !    (R,S,T) in the reference tetrahedron, the routine computes the value
    !    of the corresponding image point (X,Y,Z) in physical space.
    !
    !    This routine will also be correct for an order 10 tetrahedron,
    !    if the mapping between reference and physical space
    !    is linear.  This implies, in particular, that the sides of the
    !    image tetrahedron are straight, the faces are flat, and
    !    the "midside" nodes in the physical tetrahedron are
    !    halfway along the edges of the physical tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 December 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(3,4), the coordinates of the vertices.
    !    The vertices are assumed to be the images of (0,0,0), (1,0,0),
    !    (0,1,0) and (0,0,1) respectively.
    !
    !    Input, integer ( kind = 4 ) N, the number of points to transform.
    !
    !    Input, real ( kind = 8 ) REF(3,N), points in the reference tetrahedron.
    !
    !    Output, real ( kind = 8 ) PHY(3,N), corresponding points in the
    !    physical tetrahedron.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) i
    real    ( kind = 8 ) phy(3,n)
    real    ( kind = 8 ) ref(3,n)
    real    ( kind = 8 ) t(3,4)

    do i = 1, 3
       phy(i,1:n) = t(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) - ref(3,1:n) ) &
          + t(i,2) *             ref(1,1:n)                             &
          + t(i,3) *                          ref(2,1:n)                &
          + t(i,4) *                                       ref(3,1:n)
    end do

    return
  end subroutine tetrahedron_reference_to_physical
  subroutine tetrahedron_volume ( tetra, volume )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
    !
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 3

    real    ( kind = 8 ) a(4,4)
    real    ( kind = 8 ) tetra(dim_num,4)
    real    ( kind = 8 ) volume

    a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
    a(4,1:4) = 1.0D+00

    volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

    return
  end subroutine tetrahedron_volume
  subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
       ampm = 'AM'
    else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
       else
          ampm = 'PM'
       end if
    else
       h = h - 12
       if ( h < 12 ) then
          ampm = 'PM'
       else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
             ampm = 'Midnight'
          else
             ampm = 'AM'
          end if
       end if
    end if

    write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
       d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
  end subroutine timestamp

end module keast
