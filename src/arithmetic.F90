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

! Evaluation of arithmetic expressions. This entire module is
! THREAD-SAFE (or should be, at least).
module arithmetic
#ifdef HAVE_LIBXC
  use xc_f90_types_m, only: xc_f90_pointer_t
#endif
  implicit none

  public :: eval
  public :: fields_in_eval
  public :: tokenize
  private :: fieldeval
  private :: isnumber
  private :: isoperator
  private :: isfunction
  private :: isidentifier
  public :: iprec
  public :: iassoc
  private :: pop
  public :: istype
  private :: die
  public :: setvariable
  public :: isvariable
  public :: clearvariable
  public :: clearallvariables
  public :: listvariables
  private :: chemfunction

  private
  integer, parameter, public :: fun_openpar  = 1  !< open parenthesis
  integer, parameter, public :: fun_closepar = 2  !< close parenthesis
  integer, parameter, public :: fun_uplus    = 3  !< unary +
  integer, parameter, public :: fun_uminus   = 4  !< unary -
  integer, parameter, public :: fun_abs      = 5  !< abs(.)
  integer, parameter, public :: fun_exp      = 6  !< exp(.)
  integer, parameter, public :: fun_sqrt     = 7  !< sqrt(.)
  integer, parameter, public :: fun_floor    = 8  !< floor(.)
  integer, parameter, public :: fun_ceiling  = 9  !< ceiling(.)
  integer, parameter, public :: fun_round    = 10 !< round(.)
  integer, parameter, public :: fun_log      = 11 !< log(.)
  integer, parameter, public :: fun_log10    = 12 !< log10(.)
  integer, parameter, public :: fun_sin      = 13 !< sin(.)
  integer, parameter, public :: fun_asin     = 14 !< asin(.)
  integer, parameter, public :: fun_cos      = 15 !< cos(.)
  integer, parameter, public :: fun_acos     = 16 !< acos(.)
  integer, parameter, public :: fun_tan      = 17 !< tan(.)
  integer, parameter, public :: fun_atan     = 18 !< atan(.)
  integer, parameter, public :: fun_atan2    = 19 !< atan2(.)
  integer, parameter, public :: fun_sinh     = 20 !< sinh(.)
  integer, parameter, public :: fun_cosh     = 21 !< cosh(.)
  integer, parameter, public :: fun_erf      = 22 !< erf(.)
  integer, parameter, public :: fun_erfc     = 23 !< erfc(.)
  integer, parameter, public :: fun_min      = 24 !< min(.)
  integer, parameter, public :: fun_max      = 25 !< max(.)
  integer, parameter, public :: fun_power    = 26 !< **
  integer, parameter, public :: fun_leq      = 27 !< <=
  integer, parameter, public :: fun_geq      = 28 !< >=
  integer, parameter, public :: fun_equal    = 29 !< ==
  integer, parameter, public :: fun_neq      = 30 !< !=
  integer, parameter, public :: fun_and      = 31 !< &&
  integer, parameter, public :: fun_or       = 32 !< ||
  integer, parameter, public :: fun_plus     = 33 !< binary +
  integer, parameter, public :: fun_minus    = 34 !< binary -
  integer, parameter, public :: fun_prod     = 35 !< *
  integer, parameter, public :: fun_div      = 36 !< /
  integer, parameter, public :: fun_modulo   = 37 !< %
  integer, parameter, public :: fun_great    = 38 !< >
  integer, parameter, public :: fun_less     = 39 !< <
  integer, parameter, public :: fun_xc       = 40 !< xc(.,...)  (X)
  integer, parameter, public :: fun_gtf      = 41 !< Thomas-Fermi kinetic energy density
  integer, parameter, public :: fun_vtf      = 42 !< Potential energy density calcd using gtf and the local virial theorem
  integer, parameter, public :: fun_htf      = 43 !< Total energy density calcd using gtf and the local virial theorem
  integer, parameter, public :: fun_gtf_kir  = 44 !< Thomas-Fermi ked with Kirzhnits gradient correction
  integer, parameter, public :: fun_vtf_kir  = 45 !< Potential energy density calcd using gtf_kir and the local virial theorem
  integer, parameter, public :: fun_htf_kir  = 46 !< Total energy density calcd using gtf_kir and the local virial theorem
  integer, parameter, public :: fun_gkin     = 47 !< Kinetic energy density, G-version (grho * grho)
  integer, parameter, public :: fun_kkin     = 48 !< Kinetic energy density, K-version (rho * laprho)
  integer, parameter, public :: fun_l        = 49 !< Lagrangian density (-1/4 * laprho)
  integer, parameter, public :: fun_elf      = 50 !< Electron localization function (ELF)
  integer, parameter, public :: fun_vir      = 51 !< Electronic potential energy density, virial field
  integer, parameter, public :: fun_he       = 52 !< Energy density, G + V
  integer, parameter, public :: fun_lol      = 53 !< Localized-orbital locator (LOL)
  integer, parameter, public :: fun_lol_kir  = 54 !< Localized-orbital locator (LOL), with Kirzhnits G

#ifdef HAVE_LIBXC
  integer, parameter :: maxfun = 600
  type libxc_functional
     logical :: init = .false.
     integer :: family ! LDA, GGA, etc.
     integer :: id     ! identifier
     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional
  type(libxc_functional) :: ifun(maxfun)
#endif

  ! token type
  type token
     integer :: type = 0
     real*8 :: fval = 0d0
     integer :: ival = 0
     character(len=:), allocatable :: sval
     character*10 :: fder = ""
  end type token
  public :: token
  integer, parameter, public :: token_undef = 0
  integer, parameter, public :: token_num = 1
  integer, parameter, public :: token_fun = 2
  integer, parameter, public :: token_op = 3
  integer, parameter, public :: token_lpar = 4
  integer, parameter, public :: token_rpar = 5
  integer, parameter, public :: token_comma = 6
  integer, parameter, public :: token_field = 7
  
  ! module procedure interfaces
  interface
     recursive module function eval(expr,hardfail,iok,x0,sptr,fh,fcheck,feval,periodic)
       use iso_c_binding, only: c_ptr
       use hashmod, only: hash
       real*8 :: eval
       character(*), intent(in) :: expr
       logical, intent(in) :: hardfail
       logical, intent(out) :: iok
       real*8, intent(in), optional :: x0(3)
       type(c_ptr), intent(in), optional :: sptr
       type(hash), intent(in), optional :: fh
       optional :: fcheck, feval
       logical, intent(in), optional :: periodic
       interface
          function fcheck(sptr,id,iout)
            use iso_c_binding, only: c_ptr
            logical :: fcheck
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(out), optional :: iout
          end function fcheck
          function feval(sptr,id,nder,fder,x0,periodic)
            use iso_c_binding, only: c_ptr
            use types, only: scalar_value
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
       end interface
     end function eval
     module subroutine fields_in_eval(expr,fh,n,idlist)
       use hashmod, only: hash
       character(*), intent(in) :: expr
       type(hash), intent(in) :: fh
       integer, intent(out) :: n
       character*255, allocatable, intent(inout) :: idlist(:)
     end subroutine fields_in_eval
     module function tokenize(expr,ntok,toklist,lpexit,fh) 
       use hashmod, only: hash
       logical :: tokenize
       character(*), intent(in) :: expr
       integer, intent(out) :: ntok
       type(token), intent(inout), allocatable :: toklist(:)
       integer, intent(inout) :: lpexit
       type(hash), intent(in), optional :: fh
     end function tokenize
     recursive module function fieldeval(fid,fder,x0,sptr,fcheck,feval,periodic)
       use iso_c_binding, only: c_ptr
       real*8 :: fieldeval
       character*(*), intent(in) :: fid
       character*(*), intent(in) :: fder
       real*8, intent(in), optional :: x0(3)
       type(c_ptr), intent(in), optional :: sptr
       optional :: fcheck, feval
       logical, intent(in), optional :: periodic
       interface
          function fcheck(sptr,id,iout)
            use iso_c_binding, only: c_ptr
            logical :: fcheck
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(out), optional :: iout
          end function fcheck
          function feval(sptr,id,nder,fder,x0,periodic)
            use iso_c_binding, only: c_ptr
            use types, only: scalar_value
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
       end interface
     end function fieldeval
     module function isnumber(rval,expr,lp)
       logical :: isnumber
       character*(*), intent(in) :: expr
       integer, intent(inout) :: lp
       real*8, intent(out) :: rval
     end function isnumber
     module function isoperator(c,expr,lp)
       logical :: isoperator
       character*(*), intent(in) :: expr
       integer, intent(inout) :: lp
       integer, intent(out) :: c
     end function isoperator
     module function isfunction(c,expr,lp,wasop)
       logical :: isfunction
       character*(*), intent(in) :: expr
       integer, intent(inout) :: lp
       integer, intent(out) :: c
       logical, intent(in) :: wasop
     end function isfunction
     module function isidentifier(id,expr,lp,fder)
       logical :: isidentifier
       character(len=:), allocatable, intent(out) :: id
       character*(*), intent(in) :: expr
       integer, intent(inout) :: lp
       character*(*), intent(out), optional :: fder
     end function isidentifier
     module function iprec(c)
       integer :: iprec
       integer, intent(in) :: c
     end function iprec
     module function iassoc(c)
       integer :: iassoc
       integer, intent(in) :: c
     end function iassoc
     module subroutine pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,fail)
       use iso_c_binding, only: c_ptr
       real*8, intent(inout) :: q(:)
       integer, intent(inout) :: s(:)
       integer, intent(inout) :: nq, ns
       real*8, intent(in), optional :: x0(3)
       type(c_ptr), intent(in), optional :: sptr
       optional :: fcheck, feval
       logical, intent(in), optional :: periodic
       logical, intent(out) :: fail
       interface
          function fcheck(sptr,id,iout)
            use iso_c_binding, only: c_ptr
            logical :: fcheck
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(out), optional :: iout
          end function fcheck
          function feval(sptr,id,nder,fder,x0,periodic)
            use iso_c_binding, only: c_ptr
            use types, only: scalar_value
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
       end interface
     end subroutine pop
     module function istype(c,type)
       integer, intent(in) :: c
       character*(*), intent(in) :: type
       logical :: istype
     endfunction istype
     module subroutine die(msg,msg2)
       character*(*), intent(in) :: msg
       character*(*), intent(in), optional :: msg2
     end subroutine die
     module subroutine setvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(in) :: ival
     end subroutine setvariable
     module function isvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(out) :: ival
       logical :: isvariable
     end function isvariable
     module subroutine clearvariable(ikey)
       character*(*), intent(in) :: ikey
     end subroutine clearvariable
     module subroutine clearallvariables()
     end subroutine clearallvariables
     module subroutine listvariables()
     end subroutine listvariables
     module function chemfunction(c,sia,x0,sptr,feval,periodic) result(q)
       use iso_c_binding, only: c_ptr
       integer, intent(in) :: c
       character*(*), intent(in) :: sia
       real*8, intent(in) :: x0(3)
       real*8 :: q
       logical, intent(in), optional :: periodic
       type(c_ptr), intent(in), optional :: sptr
       interface
          function feval(sptr,id,nder,fder,x0,periodic)
            use types, only: scalar_value
            use iso_c_binding, only: c_ptr
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
       end interface
     end function chemfunction
     module function isspecialfield(fid)
       character*(*), intent(in) :: fid
       logical :: isspecialfield
     end function isspecialfield
  end interface

contains

end module arithmetic
