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
  use hashmod, only: hash
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
     module subroutine fields_in_eval(expr,fh,n,idlist)
       character(*), intent(in) :: expr
       type(hash), intent(in) :: fh
       integer, intent(out) :: n
       character*255, allocatable, intent(inout) :: idlist(:)
     end subroutine fields_in_eval
     module function tokenize(expr,ntok,toklist,lpexit,fh) 
       logical :: tokenize
       character(*), intent(in) :: expr
       integer, intent(out) :: ntok
       type(token), intent(inout), allocatable :: toklist(:)
       integer, intent(inout) :: lpexit
       type(hash), intent(in), optional :: fh
     end function tokenize
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
     module function isspecialfield(fid)
       character*(*), intent(in) :: fid
       logical :: isspecialfield
     end function isspecialfield
  end interface

contains

  !> Evaluate an arithmetic expression expr. If the expression
  !> contains fields ($), use x0 as the evaluation point
  !> (Cartesian). If hardfail is true, stop with error if the
  !> expression fails to evaluate. Otherwise, return the exit status
  !> in iok (tue=success). If periodic is present and false, evaluate
  !> the expression at x0 considering the field as
  !> non-periodic. Otherwise, evaluate it as if in a periodic system.
  !> This routine is thread-safe.
  recursive function eval(expr,hardfail,iok,x0,sptr,fh,fcheck,feval,periodic)
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
       !> Check that the id is a grid and is a sane field
       function fcheck(sptr,id,iout)
         use iso_c_binding, only: c_ptr
         logical :: fcheck
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(out), optional :: iout
       end function fcheck
       !> Evaluate the field at a point
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

    integer :: i, ntok, lp
    integer :: c, s(100)
    logical :: again, ok, ifail
    real*8 :: q(100)
    integer :: nq, ns
    type(token), allocatable :: toklist(:)

    ! tokenize the expression in input
    iok = .false.
    eval = 0d0
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,fh)
    if (.not.ok) then
       call dofail(expr(1:lp-1) // " -- " // expr(lp:))
       return
    end if

    ! initialize
    nq = 0
    ns = 0
    q(1) = 0d0

    ! run over tokens
    do i = 1, ntok
       if (toklist(i)%type == token_num) then
          ! a number
          nq = nq + 1
          q(nq) = toklist(i)%fval
       elseif (toklist(i)%type == token_fun) then
          ! a function
          ns = ns + 1
          s(ns) = toklist(i)%ival
       elseif (toklist(i)%type == token_op) then
          ! a binary operator
          c = toklist(i)%ival
          again = .true.
          do while (again)
             again = .false.
             if (ns > 0) then
                if (iprec(c) < iprec(s(ns)) .or. iassoc(c)==-1 .and. iprec(c)<=iprec(s(ns))) then
                   call pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,ifail)
                   if (ifail) then
                      call dofail()
                      return
                   endif
                   again = .true.
                end if
             end if
          end do
          ns = ns + 1
          s(ns) = c
       elseif (toklist(i)%type == token_lpar) then
          ! left parenthesis
          ns = ns + 1
          s(ns) = fun_openpar
       elseif (toklist(i)%type == token_rpar) then
          ! right parenthesis
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,ifail)
              if (ifail) then
                 call dofail()
                 return
              end if
           end do
           if (ns == 0) then
              call dofail('mismatched parentheses')
              return
           end if
           ns = ns - 1
           ! if the top of the stack is a function, pop it
           if (ns > 0) then
              c = s(ns)
              if (istype(c,'function')) then
                 call pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,ifail)
                 if (ifail) then
                    call dofail()
                    return
                 end if
              end if
           end if
        elseif (toklist(i)%type == token_comma) then
           ! a comma
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,ifail)
              if (ifail) then
                 call dofail()
                 return
              endif
           end do
           if (s(ns) /= fun_openpar) then
              call dofail('mismatched parentheses')
              return
           end if
        elseif (toklist(i)%type == token_field) then
           ! a field
           nq = nq + 1
           if (present(x0).and.present(fcheck).and.present(feval)) then
              q(nq) = fieldeval(toklist(i)%sval,toklist(i)%fder,x0,sptr,fcheck,feval,periodic)
           else
              call dofail()
              return
           end if
        else
           call dofail()
           return
        end if
    end do

    ! unwind the stack
    do while (ns > 0)
       call pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,ifail)
       if (ifail) then
          call dofail()
          return
       endif
    end do
    iok = .true.
    eval = q(1)
    return

  contains
    subroutine dofail(errmsg)
      character*(*), intent(in), optional :: errmsg
      
      iok = .false.
      if (hardfail) then
         if (present(errmsg)) then
            call die("Error evaluating expression. ",errmsg)
         else
            call die("Error evaluating expression. ",expr)
         end if
      else
         return
      endif
    end subroutine dofail

  end function eval

  !> Evaluate field with identifier fid and field flag fder at 
  !> point x0. sptr is the system C pointer. fcheck checks whether
  !> the field is sane. feval is the evaluation function. If periodic
  !> is true, evaluate the field under pbc.
  recursive function fieldeval(fid,fder,x0,sptr,fcheck,feval,periodic)
    use tools_io, only: string, isinteger, lower
    use types, only: scalar_value
    use iso_c_binding, only: c_ptr
    real*8 :: fieldeval
    character*(*), intent(in) :: fid
    character*(*), intent(in) :: fder
    real*8, intent(in), optional :: x0(3) !< position
    type(c_ptr), intent(in), optional :: sptr
    optional :: fcheck, feval
    logical, intent(in), optional :: periodic

    interface
       !> Check that the id is a grid and is a sane field
       function fcheck(sptr,id,iout)
         use iso_c_binding, only: c_ptr
         logical :: fcheck
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(out), optional :: iout
       end function fcheck
       !> Evaluate the field at a point
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

    integer :: nder, ider, lp
    type(scalar_value) :: res
    logical :: ok
    character*10 :: fderl

    fieldeval = 0d0
    if (present(x0).and.present(feval).and.present(fcheck)) then
       if (.not.fcheck(sptr,fid).and..not.isspecialfield(fid)) &
          call die('wrong field in expression: ' // string(fid))
       fderl = lower(fder)
       select case (trim(fderl))
       case ("","v","c")
          nder = 0
       case ("x","y","z","g")
          nder = 1
       case ("xx","xy","xz","yx","yy","yz","zx","zy","zz","l","lv","lc")
          nder = 2
       case default
          ! let feval interpret fder - field-specific (e.g. a MO)
          nder = -1
       end select
       res = feval(sptr,fid,nder,fder,x0,periodic)

       select case (trim(fderl))
       case ("")
          fieldeval = res%f
       case ("v")
          fieldeval = res%fval
       case ("c")
          fieldeval = res%f - res%fval
       case ("x")
          fieldeval = res%gf(1)
       case ("y")
          fieldeval = res%gf(2)
       case ("z")
          fieldeval = res%gf(3)
       case ("xx")
          fieldeval = res%hf(1,1)
       case ("xy")
          fieldeval = res%hf(1,2)
       case ("xz")
          fieldeval = res%hf(1,3)
       case ("yy")
          fieldeval = res%hf(2,2)
       case ("yz")
          fieldeval = res%hf(2,3)
       case ("zz")
          fieldeval = res%hf(3,3)
       case ("l")
          fieldeval = res%del2f
       case ("lv")
          fieldeval = res%del2fval
       case ("lc")
          fieldeval = res%del2f - res%del2fval
       case ("g")
          fieldeval = res%gfmod
       case default
          fieldeval = res%fspc
       end select
    else
       call die('evaluating field ' // string(fid) // ' without point')
    end if

  end function fieldeval

  !> Pop from the stack and operate on the queue.  This routine is
  !> thread-safe.
  subroutine pop(q,nq,s,ns,x0,sptr,fcheck,feval,periodic,fail)
    use tools_io, only: string
    use iso_c_binding, only: c_ptr
#ifdef HAVE_LIBXC
    use xc_f90_types_m
    use libxc_funcs_m
    use xc_f90_lib_m
#endif

    real*8, intent(inout) :: q(:)
    integer, intent(inout) :: s(:)
    integer, intent(inout) :: nq, ns
    real*8, intent(in), optional :: x0(3)
    type(c_ptr), intent(in), optional :: sptr
    optional :: fcheck, feval
    logical, intent(in), optional :: periodic
    logical, intent(out) :: fail

    interface
       !> Check that the id is a grid and is a sane field
       function fcheck(sptr,id,iout)
         use iso_c_binding, only: c_ptr
         logical :: fcheck
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(out), optional :: iout
       end function fcheck
       !> Evaluate the field at a point
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

    integer :: ia
    integer :: c
    real*8 :: a, b
    character*8 :: sia
#ifdef HAVE_LIBXC
    real*8 :: rho, grho, lapl, tau, zk
#endif

    ! pop from the stack
    fail = .false.
    if (ns == 0) call die('error in expression')
    c = s(ns)
    ns = ns - 1

    ! And apply to the queue
    if (c == fun_xc) then
       ! Functional from the xc library
#ifdef HAVE_LIBXC
       ia = tointeger(q(nq))

       if (.not.ifun(ia)%init) then
          ifun(ia)%id = ia
          ifun(ia)%family = xc_f90_family_from_id(ifun(ia)%id)
          call xc_f90_func_init(ifun(ia)%conf,ifun(ia)%info,ifun(ia)%id,XC_UNPOLARIZED)
          ifun(ia)%init = .true.
       endif

       if (ifun(ia)%family == XC_FAMILY_LDA .and. nq <= 1.or.&
          ifun(ia)%family == XC_FAMILY_GGA .and. nq <= 2.or.&
          ifun(ia)%family == XC_FAMILY_MGGA .and. nq <= 4) then
          call die("insufficient argument list in xc")
       endif

       select case(ifun(ia)%family)
       case (XC_FAMILY_LDA)
          rho = max(q(nq-1),1d-14)
          call xc_f90_lda_exc(ifun(ia)%conf, 1, rho, zk)
          nq = nq - 1
       case (XC_FAMILY_GGA)
          rho = max(q(nq-2),1d-14)
          grho = q(nq-1)*q(nq-1)
          call xc_f90_gga_exc(ifun(ia)%conf, 1, rho, grho, zk)
          nq = nq - 2
       case (XC_FAMILY_MGGA)
          rho = max(q(nq-4),1d-14)
          grho = q(nq-3)*q(nq-3)
          lapl = q(nq-2)
          tau = 2 * q(nq-1)
          call xc_f90_mgga_exc(ifun(ia)%conf, 1, rho, grho, lapl, tau, zk)
          nq = nq - 4
       end select
       q(nq) = zk * rho
#else
       call die('(/"!! ERROR !! critic2 was not compiled with libxc support !!"/)')
#endif
    elseif (istype(c,'binary')) then
       ! a binary operator or function
       if (nq < 2) call die('error in expression')
       a = q(nq-1)
       b = q(nq)
       nq = nq - 1
       select case(c)
       case (fun_plus)
          q(nq) = a + b
       case (fun_minus)
          q(nq) = a - b
       case (fun_prod)
          q(nq) = a * b
       case (fun_div)
          q(nq) = a / b
       case (fun_power)
          q(nq) = a ** b
       case (fun_modulo)
          q(nq) = modulo(a,b)
       case (fun_atan2)
          q(nq) = atan2(a,b)
       case (fun_min)
          q(nq) = min(a,b)
       case (fun_max)
          q(nq) = max(a,b)
       case (fun_great)
          if (a > b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_less)
          if (a < b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_leq)
          if (a <= b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_geq)
          if (a >= b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_equal)
          if (a == b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_neq)
          if (a /= b) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_and)
          if (.not.(a == 0d0).and..not.(b == 0d0)) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       case (fun_or)
          if (.not.(a == 0d0).or..not.(b == 0d0)) then
             q(nq) = 1d0
          else
             q(nq) = 0d0
          endif
       end select
    elseif (istype(c,'unary')) then
       ! a unary operator or function
       if (nq < 1) call die('error in expression')
       select case(c)
       case (fun_uplus)
          q(nq) = +q(nq)
       case (fun_uminus)
          q(nq) = -q(nq)
       case (fun_abs)
          q(nq) = abs(q(nq))
       case (fun_exp)
          q(nq) = exp(q(nq))
       case (fun_sqrt)
          q(nq) = sqrt(q(nq))
       case (fun_floor)
          q(nq) = floor(q(nq))
       case (fun_ceiling)
          q(nq) = ceiling(q(nq))
       case (fun_round)
          q(nq) = nint(q(nq))
       case (fun_log)
          q(nq) = log(q(nq))
       case (fun_log10)
          q(nq) = log10(q(nq))
       case (fun_sin)
          q(nq) = sin(q(nq))
       case (fun_asin)
          q(nq) = asin(q(nq))
       case (fun_cos)
          q(nq) = cos(q(nq))
       case (fun_acos)
          q(nq) = acos(q(nq))
       case (fun_tan)
          q(nq) = tan(q(nq))
       case (fun_atan)
          q(nq) = atan(q(nq))
       case (fun_sinh)
          q(nq) = sinh(q(nq))
       case (fun_cosh)
          q(nq) = cosh(q(nq))
       case (fun_erf)
          q(nq) = erf(q(nq))
       case (fun_erfc)
          q(nq) = erfc(q(nq))
       end select
    elseif (istype(c,'chemfunction')) then
       ! We need a point and the evaluator
       if (.not.present(x0).or..not.present(feval).or..not.present(fcheck)) then
          fail = .true.
          return
       endif
    
       ! Also an integer field identifier as the first argument
       ia = tointeger(q(nq))
       write (sia,'(I8)') ia
       sia = adjustl(sia)
       if (.not.fcheck(sptr,sia)) &
          call die('wrong field ' // string(sia))
    
       ! Use the library of chemical functions
       q(nq) = chemfunction(c,sia,x0,sptr,feval,periodic)
    else
       call die('error in expression')
    end if
  contains
    function tointeger(a) result(ia)
      real*8 :: a
      integer :: ia

      if (a - nint(a) > 1d-13) call die("An integer was required in this expression.")
      ia = nint(a)
    endfunction tointeger
  end subroutine pop

  !> Calculate a chemical function for a given field.  This routine is
  !> thread-safe.
  function chemfunction(c,sia,x0,sptr,feval,periodic) result(q)
    use types, only: scalar_value
    use param, only: pi
    use iso_c_binding, only: c_ptr
    integer, intent(in) :: c
    character*(*), intent(in) :: sia
    real*8, intent(in) :: x0(3)
    real*8 :: q
    logical, intent(in), optional :: periodic
    type(c_ptr), intent(in), optional :: sptr
  
    interface
       !> Evaluate the field at a point
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
  
    type(scalar_value) :: res
    real*8 :: f0, ds, ds0, g, g0
  
    ! a constant
    real*8 :: ctf = 2.8712340001881911d0 ! Thomas-Fermi k.e.d. constant, 3/10 * (3*pi^2)^(2/3)
  
    q = 0d0
    select case(c)
    case (fun_gtf)
       ! Thomas-Fermi kinetic energy density for the uniform electron gas
       ! See Yang and Parr, Density-Functional Theory of Atoms and Molecules
       res = feval(sptr,sia,0,"",x0,periodic)
       q = ctf * res%f**(5d0/3d0)
    case (fun_vtf)
       ! Potential energy density calculated using fun_gtf and the local
       ! virial theorem (2g(r) + v(r) = 1/4*lap(r)).
       res = feval(sptr,sia,2,"",x0,periodic)
       q = ctf * res%f**(5d0/3d0)
       q = 0.25d0 * res%del2f - 2 * q
    case (fun_htf)
       ! Total energy density calculated using fun_gtf and the local
       ! virial theorem (h(r) + v(r) = 1/4*lap(r)).
       res = feval(sptr,sia,2,"",x0,periodic)
       q = ctf * res%f**(5d0/3d0)
       q = 0.25d0 * res%del2f - q
    case (fun_gtf_kir)
       ! Thomas-Fermi kinetic energy density with Kirzhnits'
       ! semiclassical gradient correction for inhomogeinities. See:
       !   Kirzhnits, D. (1957). Sov. Phys. JETP, 5, 64-72.
       !   Kirzhnits, D. Field Theoretical Methods in Many-body Systeins (Pergamon, New York, 1967)
       ! and also:
       !   Zhurova and Tsirelson, Acta Cryst. B (2002) 58, 567-575.
       ! for more references and its use applied to experimental electron
       ! densities.
       res = feval(sptr,sia,2,"",x0,periodic)
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
    case (fun_vtf_kir)
       ! Potential energy density calculated using fun_gtf_kir and the
       ! local virial theorem (2g(r) + v(r) = 1/4*lap(r)).
       res = feval(sptr,sia,2,"",x0,periodic)
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = 0.25d0 * res%del2f - 2 * q
    case (fun_htf_kir)
       ! Total energy density calculated using fun_gtf_kir and the
       ! local virial theorem (h(r) + v(r) = 1/4*lap(r)).
       res = feval(sptr,sia,2,"",x0,periodic)
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = 0.25d0 * res%del2f - q
    case (fun_gkin)
       ! G-kinetic energy density (sum grho * grho)
       res = feval(sptr,sia,1,"",x0,periodic)
       q = res%gkin
    case (fun_kkin)
       ! K-kinetic energy density (sum rho * laprho)
       res = feval(sptr,sia,2,"",x0,periodic)
       q = res%gkin - 0.25d0 * res%del2f
    case (fun_l)
       ! Lagrangian density (-1/4 * lap)
       res = feval(sptr,sia,2,"",x0,periodic)
       q = - 0.25d0 * res%del2f
    case (fun_elf)
       ! Electron localization function
       ! Becke and Edgecombe J. Chem. Phys. (1990) 92, 5397-5403
       res = feval(sptr,sia,1,"",x0,periodic)
       if (res%f < 1d-30) then
          q = 0d0
       else
          f0 = max(res%f,1d-30)
          ds = res%gkin  - 1d0/8d0 * res%gfmod**2 / f0
          ds0 = ctf * f0**(5d0/3d0)
          q = ds / ds0
          q = 1d0 / (1d0 + q**2)
       end if
    case (fun_vir)
       ! Electronic potential energy density (virial field)
       ! Keith et al. Int. J. Quantum Chem. (1996) 57, 183-198.
       res = feval(sptr,sia,2,"",x0,periodic)
       q = res%vir
    case (fun_he)
       ! Energy density, fun_vir + fun_gkin
       !   Keith et al. Int. J. Quantum Chem. (1996) 57, 183-198.
       res = feval(sptr,sia,2,"",x0,periodic)
       q = res%vir + res%gkin
    case (fun_lol)
       ! Localized-orbital locator
       !   Schmider and Becke, J. Mol. Struct. (Theochem) (2000) 527, 51-61
       !   Schmider and Becke, J. Chem. Phys. (2002) 116, 3184-3193.
       res = feval(sptr,sia,2,"",x0,periodic)
       q = ctf * res%f**(5d0/3d0) / max(res%gkin,1d-30)
       q = q / (1d0+q)
    case (fun_lol_kir)
       ! Localized-orbital locator using Kirzhnits k.e.d.
       !   Tsirelson and Stash, Acta Cryst. (2002) B58, 780.
       res = feval(sptr,sia,2,"",x0,periodic)
       f0 = max(res%f,1d-30)
       g0 = ctf * f0**(5d0/3d0) 
       g = g0 + 1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = g0 / g
       q = q / (1d0+q)
    end select
  
  end function chemfunction

end module arithmetic
