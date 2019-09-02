! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Standalone library for I/O and error handling. Heavily modified
! from previous code (pi). 
module tools_io
  implicit none

  public
  public :: string
  private :: string_int4
  private :: string_int8
  private :: string_real8
  private :: string_real4
  private :: string_char
  private :: string_logical
  public :: equal
  public :: equali
  public :: getline
  public :: getline_raw
  public :: nameguess
  public :: zatguess
  public :: upper
  public :: lower
  public :: getword
  public :: lgetword
  public :: isletter
  public :: isdigit
  public :: isinteger
  public :: isreal
  public :: isexpression
  public :: isexpression_or_word
  public :: isassignment
  public :: ioinit
  public :: fopen_read
  public :: fopen_write
  public :: fopen_append
  public :: fopen_scratch
  public :: fclose
  public :: falloc
  public :: fdealloc
  public :: ferror
  public :: tictac

  ! Justification parameters
  integer, parameter :: ioj_left = -1
  integer, parameter :: ioj_right = -2
  integer, parameter :: ioj_center = -3
  integer, parameter :: ioj_atpoint = -4

  ! General string-building routine
  interface string
     module procedure string_int4
     module procedure string_int8
     module procedure string_real8
     module procedure string_real4
     module procedure string_char
     module procedure string_logical
  end interface string

  ! character parameters
  character*(1), parameter :: tab = char(9) !< tab char (ascii 9)

  ! main input and output, logical unit allocation
  integer :: uin !< main input lu
  integer :: uout !< main output lu
  integer :: ucopy !< logical unit where the copy of the input is written
  character(len=:), allocatable :: filepath !< relative path to find related files
  logical :: interactive !< is this an interactive sesion?

  ! error system
  integer, parameter :: faterr = -1 !< fatal error flag
  integer, parameter :: warning = 1 !< warning flag
  integer, parameter :: noerr = 0 !< info flag
  integer :: nwarns = 0 !< Number of warnings
  integer :: ncomms = 0 !< Number of comments
  
  interface
     module subroutine stdargs(optv,ghome,uroot)
       character(len=:), allocatable, intent(out) :: optv
       character(len=:), allocatable, intent(out) :: ghome
       character(len=:), allocatable, intent(out) :: uroot
     end subroutine stdargs
     module function string_int4(a,length,justify,pad0,padspace) result(s)
       character(len=:), allocatable :: s
       integer*4, intent(in) :: a
       integer, intent(in), optional :: length
       integer, intent(in), optional :: justify
       logical, intent(in), optional :: pad0
       integer, intent(in), optional :: padspace
     end function string_int4
     module function string_int8(a,length,justify,pad0,padspace) result(s)
       character(len=:), allocatable :: s
       integer*8, intent(in) :: a
       integer, intent(in), optional :: length
       integer, intent(in), optional :: justify
       logical, intent(in), optional :: pad0
       integer, intent(in), optional :: padspace
     end function string_int8
     module function string_real8(a,desc,length,decimal,justify,padspace) result(s)
       character(len=:), allocatable :: s
       real*8, intent(in) :: a
       character*1, intent(in) :: desc
       integer, intent(in), optional :: length
       integer, intent(in), optional :: decimal
       integer, intent(in), optional :: justify
       integer, intent(in), optional :: padspace
     end function string_real8
     module function string_real4(a,desc,length,decimal,justify,padspace) result(s)
       character(len=:), allocatable :: s
       real*4, intent(in) :: a
       character*1, intent(in) :: desc
       integer, intent(in), optional :: length
       integer, intent(in), optional :: decimal
       integer, intent(in), optional :: justify
       integer, intent(in), optional :: padspace
     end function string_real4
     module function string_char(a,length,justify,padspace) result(s)
       character(len=:), allocatable :: s
       character(len=*), intent(in) :: a
       integer, intent(in), optional :: length
       integer, intent(in), optional :: justify
       integer, intent(in), optional :: padspace
     end function string_char
     module function string_logical(a) result(s)
       character(len=:), allocatable :: s
       logical, intent(in) :: a
     end function string_logical
     module function equal(s,t)
       character*(*), intent(in) :: s
       character*(*), intent(in) :: t
       logical :: equal
     end function equal
     module function equali(s,t)
       character*(*), intent(in) :: s
       character*(*), intent(in) :: t
       logical :: equali
     end function equali
     module function getline(u,oline,eofstop,ucopy)
       character(len=:), allocatable, intent(out) :: oline
       integer, intent(in) :: u
       logical, intent(in), optional :: eofstop
       integer, intent(in), optional :: ucopy
       logical :: getline
     end function getline
     module function getline_raw(u,line,eofstop) result(ok)
       integer, intent(in) :: u
       character(len=:), allocatable, intent(out) :: line
       logical, intent(in), optional :: eofstop
       logical :: ok
     end function getline_raw
     module function zatguess(atname)
       use param, only: maxzat0
       integer :: zatguess
       character*(*), intent(in) :: atname
     end function zatguess
     module function nameguess (zat,nounderscore)
       use param, only: maxzat0
       integer, intent(in) :: zat
       logical, intent(in), optional :: nounderscore
       character*(2) :: nameguess
     end function nameguess
     module function upper(a)
       character*(*), intent(in) :: a
       character*(len(a)) :: upper
     end function upper
     module function lower(a)
       character*(*), intent(in) :: a
       character*(len(a)) :: lower
     end function lower
     module function deblank(a,toup,todn)
       character*(*), intent(in) :: a
       logical, intent(in), optional :: toup
       logical, intent(in), optional :: todn
       character*(len(a)) :: deblank
     end function deblank
     module function stripchar(a,ch)
       character*(*), intent(in) :: a
       character, intent(in) :: ch
       character*(len(a)) :: stripchar
     end function stripchar
     module function getword(line,lp) result(word)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
       character(len=:), allocatable :: word
     end function getword
     module function lgetword (line,lp) result(word)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
       character(len=:), allocatable :: word
     endfunction lgetword
     module function isletter(a)
       character*(1), intent(in) :: a
       logical :: isletter
     end function isletter
     module function isdigit(a)
       character*(1), intent(in) :: a
       logical :: isdigit
     end function isdigit
     module function isinteger(ival,line,lp0)
       logical :: isinteger
       integer, intent(inout) :: ival
       character*(*), intent(in) :: line
       integer, intent(inout), optional :: lp0
     end function isinteger
     module function isreal(rval, line, lp0)
       real*8, intent(out) :: rval
       character*(*), intent(in) :: line
       integer, intent(inout), optional :: lp0
       logical :: isreal
     end function isreal
     module function isexpression(word, line, lp0)
       character(len=:), allocatable, intent(out) :: word
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp0
       logical :: isexpression
     end function isexpression
     module function isexpression_or_word(word, line, lp0)
       character(len=:), allocatable, intent(out) :: word
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp0
       logical :: isexpression_or_word
     end function isexpression_or_word
     module function isassignment(var, expr, line)
       use param, only: constlist, funclist
       character*(*), intent(in) :: line
       character(len=:), allocatable, intent(out) :: var
       character(len=:), allocatable, intent(out) :: expr
       logical :: isassignment
     end function isassignment
     module subroutine ioinit ()
     end subroutine ioinit
     module function fopen_read(file,form,abspath0,errstop) result(lu)
       character*(*), intent(in) :: file
       character*(*), intent(in), optional :: form
       logical, intent(in), optional :: abspath0
       logical, intent(in), optional :: errstop
       integer :: lu
     end function fopen_read
     module function fopen_write(file,form,abspath0,errstop) result(lu)
       character*(*), intent(in) :: file
       character*(*), intent(in), optional :: form
       logical, intent(in), optional :: abspath0
       logical, intent(in), optional :: errstop
       integer :: lu
     end function fopen_write
     module function fopen_append(file,form,abspath0,errstop) result(lu)
       character*(*), intent(in) :: file
       character*(*), intent(in), optional :: form
       logical, intent(in), optional :: abspath0
       logical, intent(in), optional :: errstop
       integer :: lu
     end function fopen_append
     module function fopen_scratch(form,errstop) result(lu)
       character*(*), intent(in), optional :: form
       logical, intent(in), optional :: errstop
       integer :: lu
     end function fopen_scratch
     module subroutine fclose(lu)
       integer, intent(in) :: lu
     end subroutine fclose
     module function falloc()
       integer :: falloc
     end function falloc
     module subroutine fdealloc(u)
       integer, intent(in) :: u
     end subroutine fdealloc
     module subroutine ferror(routine,message,errortype,inputline,syntax)
       character*(*), intent(in) :: routine
       character*(*), intent(in) :: message
       integer, intent(in) :: errortype
       character*(*), intent(in), optional :: inputline
       logical, intent(in), optional :: syntax
     end subroutine ferror
     module subroutine tictac(mesg)
       character*(*), intent(in) :: mesg
     end subroutine tictac
     module subroutine start_clock()
     end subroutine start_clock
     module subroutine print_clock()
     end subroutine print_clock
     module subroutine sectotime(s,x0,s0)
       real, intent(in) :: s
       integer, intent(out), optional :: x0(4)
       character(len=:), allocatable, optional :: s0
     end subroutine sectotime
  end interface

end module tools_io
