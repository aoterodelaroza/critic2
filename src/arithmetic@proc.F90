! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
submodule (arithmetic) proc
#ifdef HAVE_LIBXC
  use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_info_t
#endif
  implicit none

  !xx! private procedures
  ! function tokenize(expr,ntok,toklist,lpexit,syl,errmsg)
  ! function loaded_fields(syl) result(str)
  ! function arith_usage(c,arg) result(str)
  ! function arith_code(c) result(idx)
  ! function arith_nargs(idx) result(n)
  ! function check_nargs(c,na) result(errmsg)
  ! function iprec(c)
  ! function iassoc(c)
  ! function istype(c,type,iwantarg)
  ! function isstructvar(fid,c,fder)
  ! function arith_find(name,kind) result(idx)
  ! recursive function fieldeval(fid,fder,errmsg,x0,sptr,periodic)
  ! function structvareval(svar,x0,syl,periodic) result(q)
  ! function isnumber(rval,expr,lp)
  ! function isoperator(c,expr,lp)
  ! function isfunction(c,expr,lp,wasop)
  ! function isidentifier(id,expr,lp,fder)
  ! subroutine pop(q,nq,s,ns,x0,sptr,periodic,errmsg,nargs)
  ! subroutine pop_grid(q,nq,s,ns,fail)
  ! function chemfunction(c,sia,x0,args,sptr,periodic,errmsg) result(q)

  ! libxc functional
#ifdef HAVE_LIBXC
  integer, parameter :: maxfun = 800
  type libxc_functional
     logical :: init = .false.
     integer :: family ! LDA, GGA, etc.
     integer :: id     ! identifier
     type(xc_f03_func_t) :: conf ! the pointer used to call the library
     type(xc_f03_func_info_t) :: info ! information about the functional
  end type libxc_functional
  type(libxc_functional) :: ifun(maxfun)
#endif

  integer, parameter :: token_undef = 0
  integer, parameter :: token_num = 1
  integer, parameter :: token_fun = 2
  integer, parameter :: token_op = 3
  integer, parameter :: token_lpar = 4
  integer, parameter :: token_rpar = 5
  integer, parameter :: token_comma = 6
  integer, parameter :: token_field = 7
  integer, parameter :: token_structvar = 8

contains

  !> Evaluate an arithmetic expression expr. If the expression
  !> contains fields ($), use x0 as the evaluation point
  !> (Cartesian). If hardfail is true, stop with error if the
  !> expression fails to evaluate. Otherwise, return the exit status
  !> in iok (tue=success). If periodic is present and false, evaluate
  !> the expression at x0 considering the field as
  !> non-periodic. Otherwise, evaluate it as if in a periodic system.
  !> sptr C pointer to the calling system. toklistin = use the token
  !> list provided by the user instead of the expression. This routine
  !> is thread-safe.
  recursive module function eval(expr,errmsg,x0,sptr,periodic,toklistin)
    use ieee_arithmetic, only: ieee_is_finite
    use ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_usual
    use iso_c_binding, only: c_ptr
    real*8 :: eval
    character(*), intent(in) :: expr
    character(len=:), allocatable, intent(inout) :: errmsg
    real*8, intent(in), optional :: x0(3)
    type(c_ptr), intent(in), optional :: sptr
    logical, intent(in), optional :: periodic
    type(token), intent(in), optional, target :: toklistin(:)

    logical :: usual_halt(3)

    ! Suspend FPE halting while the expression is evaluated: the guarded
    ! operations report an error instead, and the remaining ones (power,
    ! modulo, overflows, ...) must not abort the program on a
    ! user-provided expression. The result is checked for finiteness.
    call ieee_get_halting_mode(ieee_usual,usual_halt)
    call ieee_set_halting_mode(ieee_usual,.false.)
    eval = eval_core(expr,errmsg,x0,sptr,periodic,toklistin)
    call ieee_set_halting_mode(ieee_usual,usual_halt)
    if (len_trim(errmsg) == 0) then
       if (.not.ieee_is_finite(eval)) &
          errmsg = 'expression evaluates to a non-finite value (Inf/NaN): ' // trim(expr)
    end if

  end function eval

  recursive function eval_core(expr,errmsg,x0,sptr,periodic,toklistin) result(eval)
    use systemmod, only: system
    use iso_c_binding, only: c_ptr, c_f_pointer
    use hashmod, only: hash
    real*8 :: eval
    character(*), intent(in) :: expr
    character(len=:), allocatable, intent(inout) :: errmsg
    real*8, intent(in), optional :: x0(3)
    type(c_ptr), intent(in), optional :: sptr
    logical, intent(in), optional :: periodic
    type(token), intent(in), optional, target :: toklistin(:)

    integer :: i, ntok, lp
    integer :: c, s(100), sb(100), nb
    logical :: again, ok
    real*8 :: q(100)
    integer :: nq, ns
    type(token), allocatable, target :: toklist(:)
    type(token), pointer :: tl(:)
    type(system), pointer :: syl

    errmsg = ""

    ! recover the system pointer
    syl => null()
    eval = 0d0
    if (present(sptr)) then
       call c_f_pointer(sptr,syl)
       if (.not.syl%isinit) then
          errmsg = 'system not initialized'
          return
       end if
    end if

    ! tokenize the expression, or point at the token list provided in input
    if (present(toklistin)) then
       tl => toklistin
       ntok = size(toklistin,1)
    else
       lp = 1
       ok = tokenize(expr,ntok,toklist,lp,syl,errmsg)
       if (.not.ok) return
       tl => toklist
    end if

    ! initialize
    nq = 0
    ns = 0
    q(1) = 0d0

    ! run over tokens
    do i = 1, ntok
       if (tl(i)%type == token_num) then
          ! a number
          nq = nq + 1
          q(nq) = tl(i)%fval
       elseif (tl(i)%type == token_fun) then
          ! a function
          ns = ns + 1
          s(ns) = tl(i)%ival
       elseif (tl(i)%type == token_op) then
          ! a binary operator
          c = tl(i)%ival
          again = .true.
          do while (again)
             again = .false.
             if (ns > 0) then
                if (iprec(c) < iprec(s(ns)) .or. iassoc(c)==-1 .and. iprec(c)<=iprec(s(ns))) then
                   call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
                   if (len_trim(errmsg) > 0) return
                   again = .true.
                end if
             end if
          end do
          ns = ns + 1
          s(ns) = c
       elseif (tl(i)%type == token_lpar) then
          ! left parenthesis; if it opens the arguments of a function,
          ! remember the queue length to count them at the closing one
          ns = ns + 1
          s(ns) = fun_openpar
          sb(ns) = -1
          if (i > 1) then
             if (tl(i-1)%type == token_fun) sb(ns) = nq
          end if
       elseif (tl(i)%type == token_rpar) then
          ! right parenthesis
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
              if (len_trim(errmsg) > 0) return
           end do
           if (ns == 0) then
              errmsg = "mismatched parentheses: ')' without '('"
              return
           end if
           nb = sb(ns)
           ns = ns - 1
           ! if these were the arguments of a function, check their
           ! number and apply the function
           if (nb >= 0) then
              errmsg = check_nargs(s(ns),nq-nb)
              if (len_trim(errmsg) > 0) return
              call pop(q,nq,s,ns,x0,syl,periodic,errmsg,nq-nb)
              if (len_trim(errmsg) > 0) return
           end if
        elseif (tl(i)%type == token_comma) then
           ! a comma
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
              if (len_trim(errmsg) > 0) return
           end do
           if (ns == 0) then
              errmsg = "comma outside the arguments of a function"
              return
           end if
        elseif (tl(i)%type == token_field) then
           ! a field
           nq = nq + 1
           q(nq) = fieldeval(tl(i)%sval,tl(i)%fder,errmsg,x0,syl,periodic)
           if (len_trim(errmsg) > 0) return
        elseif (tl(i)%type == token_structvar) then
           ! a structural variable
           nq = nq + 1
           q(nq) = structvareval(tl(i)%ival,tl(i)%fder,errmsg,x0,syl,periodic)
           if (len_trim(errmsg) > 0) return
        else
           errmsg = 'syntax error'
           return
        end if
    end do

    ! unwind the stack
    do while (ns > 0)
       call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
       if (len_trim(errmsg) > 0) return
    end do
    if (nq > 1) then
       errmsg = "missing operator between two values, or too many arguments in a function"
       return
    elseif (nq == 0 .and. ntok > 0) then
       errmsg = "missing value"
       return
    end if
    eval = q(1)

  end function eval_core

  !> Evaluate an arithmetic expression (expr) on a grid with
  !> dimensions n. sptr = C pointer to the system. Returns the grid
  !> of values in f (and iok=.true.) if success, or zero and iok=.false.
  !> if failed. A failure gives no message: the callers then
  !> evaluate the expression point by point, which reports the
  !> errors.
  module subroutine eval_grid(n,expr,sptr,f,iok)
    use ieee_arithmetic, only: ieee_is_finite
    use ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_usual
    use systemmod, only: system
    use types, only: realloc
    use tools_io, only: uout, string
    use iso_c_binding, only: c_ptr, c_f_pointer
    integer, intent(in) :: n(3)
    character(*), intent(in) :: expr
    type(c_ptr), intent(in) :: sptr
    real*8, intent(out) :: f(n(1),n(2),n(3))
    logical, intent(out) :: iok

    integer :: i, ntok, lp, idx
    integer :: c, s(100), sb(100), nb
    logical :: again, ok, ifail
    integer :: nq, ns
    real*8, allocatable :: q(:,:,:,:)
    type(token), allocatable :: toklist(:)
    type(system), pointer :: syl
    logical :: usual_halt(3)
    integer :: nbad
    character(len=:), allocatable :: errmsg

    ! suspend FPE halting while the expression is evaluated (see eval)
    call ieee_get_halting_mode(ieee_usual,usual_halt)
    call ieee_set_halting_mode(ieee_usual,.false.)

    ! recover the system pointer
    syl => null()
    call c_f_pointer(sptr,syl)
    if (.not.syl%isinit) goto 999

    ! tokenize the expression in input
    iok = .true.
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,syl,errmsg)
    if (.not.ok) goto 999

    ! initialize
    nq = 0
    ns = 0
    allocate(q(n(1),n(2),n(3),1))

    ! The grid version of the arithmetic evaluator does not support
    ! certain types of operators (xc, chemfunction) or structural
    ! variables. Check that the grid fields are all correct
    do i = 1, ntok
       if (toklist(i)%ival == fun_xc) goto 999
       if (toklist(i)%type == token_structvar) goto 999
       if (istype(toklist(i)%ival,'chemfunction')) goto 999
       if (toklist(i)%type == token_field) then
          ifail = .not.syl%goodfield(key=toklist(i)%sval,n=n)
          ifail = ifail .or. (trim(toklist(i)%fder) /= "" .and. trim(toklist(i)%fder)/="v")
          if (ifail) goto 999
       end if
    end do

    ! run over tokens
    do i = 1, ntok
       if (toklist(i)%type == token_num) then
          ! a number
          nq = nq + 1
          if (nq > size(q,4)) call realloc(q,n(1),n(2),n(3),nq)
          q(:,:,:,nq) = toklist(i)%fval
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
                   call pop_grid(q,nq,s,ns,ifail)
                   if (ifail) goto 999
                   again = .true.
                end if
             end if
          end do
          ns = ns + 1
          s(ns) = c
       elseif (toklist(i)%type == token_lpar) then
          ! left parenthesis (see eval_core)
          ns = ns + 1
          s(ns) = fun_openpar
          sb(ns) = -1
          if (i > 1) then
             if (toklist(i-1)%type == token_fun) sb(ns) = nq
          end if
       elseif (toklist(i)%type == token_rpar) then
          ! right parenthesis
          do while (ns > 0)
             if (s(ns) == fun_openpar) exit
             call pop_grid(q,nq,s,ns,ifail)
             if (ifail) goto 999
          end do
          if (ns == 0) goto 999
          nb = sb(ns)
          ns = ns - 1
          ! if these were the arguments of a function, check their
          ! number and apply the function
          if (nb >= 0) then
             if (len_trim(check_nargs(s(ns),nq-nb)) > 0) goto 999
             call pop_grid(q,nq,s,ns,ifail)
             if (ifail) goto 999
          end if
       elseif (toklist(i)%type == token_comma) then
          ! a comma
          do while (ns > 0)
             if (s(ns) == fun_openpar) exit
             call pop_grid(q,nq,s,ns,ifail)
             if (ifail) goto 999
          end do
          if (ns == 0) goto 999
       elseif (toklist(i)%type == token_field) then
          ! a field (already made all the checks at the parsing stage)
          nq = nq + 1
          if (nq > size(q,4)) call realloc(q,n(1),n(2),n(3),nq)
          idx = syl%fieldname_to_idx(toklist(i)%sval)
          q(:,:,:,nq) = syl%f(idx)%grid%f
       else
          goto 999
       end if
    end do

    ! unwind the stack
    do while (ns > 0)
       call pop_grid(q,nq,s,ns,ifail)
       if (ifail) goto 999
    end do
    if (nq /= 1) goto 999
    f = q(:,:,:,1)
    if (allocated(q)) deallocate(q)
    call ieee_set_halting_mode(ieee_usual,usual_halt)

    ! warn about non-finite values in the result
    nbad = count(.not.ieee_is_finite(f))
    if (nbad > 0) &
       write (uout,'("!! Warning !! expression generated ",A," non-finite values (Inf/NaN) out of ",A,&
          &" grid points: ",A)') string(nbad), string(n(1)*n(2)*n(3)), trim(expr)

    return
999 continue
    call ieee_set_halting_mode(ieee_usual,usual_halt)
    iok = .false.
    if (allocated(q)) deallocate(q)
    f = 0d0
    return

  end subroutine eval_grid

  !> Parse an expression (expr) and return the n field ids in it
  !> (idlist).  sptr = C pointer to the calling system. If chemargs,
  !> include the fields given as arguments of chemical functions
  !> (gtf(rho), gtf(2)). This routine is thread-safe. If error,
  !> return non-zero-length string in errmsg.
  module subroutine fields_in_eval(expr,errmsg,n,idlist,sptr,chemargs)
    use systemmod, only: system
    use tools_io, only: string
    use types, only: realloc
    use param, only: mlen
    use iso_c_binding, only: c_ptr, c_f_pointer
    character(*), intent(in) :: expr
    character(len=:), allocatable, intent(inout) :: errmsg
    integer, intent(out) :: n
    character(len=mlen), allocatable, intent(inout) :: idlist(:)
    type(c_ptr), intent(in) :: sptr
    logical, intent(in), optional :: chemargs

    integer :: lp, i
    logical :: dochem, keep
    logical :: ok
    integer :: ntok
    type(token), allocatable :: toklist(:)
    type(system), pointer :: syl

    errmsg = ""
    syl => null()
    ! recover the system pointer
    call c_f_pointer(sptr,syl)
    if (.not.syl%isinit) then
       errmsg = "system not initialized"
       return
    end if

    ! allocate space for the field ids
    if (allocated(idlist)) deallocate(idlist)
    allocate(idlist(10))

    ! tokenize the expression
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,syl,errmsg)
    if (.not. ok) return

    ! return the ids of the fields in an array
    dochem = .false.
    if (present(chemargs)) dochem = chemargs
    n = 0
    do i = 1, ntok
       keep = toklist(i)%type == token_field
       if (dochem .and. toklist(i)%type == token_num) keep = len_trim(toklist(i)%sval) > 0
       if (keep) then
          n = n + 1
          if (n > size(idlist)) call realloc(idlist,2*n)
          idlist(n) = trim(toklist(i)%sval)
       end if
    end do
    call realloc(idlist,n)

  end subroutine fields_in_eval

  !> Set the value of a variable.
  module subroutine setvariable(ikey,ival)
    use param, only: vh

    character*(*), intent(in) :: ikey
    real*8, intent(in) :: ival

    call vh%put(trim(ikey),ival)

  end subroutine setvariable

  !> Determine if a variable is defined This routine is thread-safe.
  module function isvariable(ikey,ival)
    use param, only: vh

    character*(*), intent(in) :: ikey
    real*8, intent(out) :: ival
    logical :: isvariable

    isvariable = vh%iskey(trim(ikey))
    if (isvariable) then
       ival = vh%get(trim(ikey),ival)
    end if

  end function isvariable

  !> Clear a variable.
  module subroutine clearvariable(ikey)
    use param, only: vh

    character*(*), intent(in) :: ikey

    call vh%delkey(trim(ikey))

  end subroutine clearvariable

  !> Clear all the variables by freeing the hash.
  module subroutine clearallvariables()
    use param, only: vh

    call vh%free()

  end subroutine clearallvariables

  !> List of the variables in the internal database.  This routine is
  !> thread-safe.
  module subroutine listvariables()
    use tools_io, only: uout, string, ioj_right
    use param, only: vh

    integer :: i, nkeys, idum
    character(len=:), allocatable :: key, typx, val
    character*2 :: typ
    real*4 :: rdum
    real*8 :: ddum
    character*1 :: sdum

    nkeys = vh%keys()
    write (uout,'("* LIST of variables (",A,")")') string(nkeys)
    do i = 1, nkeys
       key = vh%getkey(i)
       typ = vh%type(key)
       val = ""
       select case (typ)
       case ("i_")
          typx = "integer"
          val = string(vh%get(key,idum))
       case ("r_")
          typx = "real"
          val = string(vh%get(key,rdum),'G')
       case ("d_")
          typx = "double"
          val = string(vh%get(key,ddum),'G')
       case ("s_")
          typx = "string"
          val = string(vh%get(key,sdum))
       end select
       write (uout,'(A,". ",A," = ",A)') string(i,3,ioj_right), &
          string(key), string(val)
    end do
    write (uout,*)

  end subroutine listvariables

  !> List of the variables in the internal database.  This routine is
  !> thread-safe.
  module subroutine listlibxc(doref,doname,doflags)
#ifdef HAVE_LIBXC
    use tools_io, only: uout, string, ioj_left
    use xc_f03_lib_m
#endif
    logical, intent(in) :: doref
    logical, intent(in) :: doname
    logical, intent(in) :: doflags

#ifdef HAVE_LIBXC
     integer :: i, nfun, id, ifam, ikind, iflags, mlen, iref
     integer, allocatable :: idlist(:)
     character(len=256) :: name
     character(len=2048) :: longname
     type(xc_f03_func_t) :: p
     type(xc_f03_func_reference_t) :: ref
     type(xc_f03_func_info_t) :: info
     character(len=:), allocatable :: str, aux

     ! header and some libxc info
     write (uout,'("* LIST of libxc functionals")')
     nfun = xc_f03_number_of_functionals()
     allocate(idlist(nfun))
     call xc_f03_available_functional_numbers(idlist(1)) ! hack
     mlen = xc_f03_maximum_name_length()

     ! gather and write info about the functionals in this libxc
     do i = 1, nfun
        id = idlist(i)
        name = xc_f03_functional_get_name(id)
        call xc_f03_func_init(p, id, XC_UNPOLARIZED)
        info = xc_f03_func_get_info(p)

        ifam = xc_f03_func_info_get_family(info)
        ikind = xc_f03_func_info_get_kind(info)
        iflags = xc_f03_func_info_get_flags(info)
        longname = xc_f03_func_info_get_name(info)

        str = string(id,5,ioj_left)
        str = str // " " // string(name,mlen+1,ioj_left)

       select case (ikind)
       case(XC_EXCHANGE)
          str = str // " " // "x "
       case(XC_CORRELATION)
          str = str // " " // "c "
       case(XC_EXCHANGE_CORRELATION)
          str = str // " " // "xc"
       case(XC_KINETIC)
          str = str // " " // "k "
       end select

       select case (ifam)
       case(XC_FAMILY_UNKNOWN)
          str = str // " " // " unknown"
       case(XC_FAMILY_NONE)
          str = str // " " // "  none  "
       case(XC_FAMILY_LDA)
          str = str // " " // "  LDA   "
       case(XC_FAMILY_GGA)
          str = str // " " // "  GGA   "
       case(XC_FAMILY_MGGA)
          str = str // " " // "meta-GGA"
       case(XC_FAMILY_LCA)
          str = str // " " // "  LCA   "
       case(XC_FAMILY_OEP)
          str = str // " " // "  OEP   "
       case(XC_FAMILY_HYB_GGA)
          str = str // " " // " hyb-GGA"
       case(XC_FAMILY_HYB_MGGA)
          str = str // " " // "hyb-mGGA"
       end select

       if (doflags) then
          aux = "-"
          if (iand(iflags,XC_FLAGS_HAVE_EXC) /= 0) aux = aux // "E-"
          if (iand(iflags,XC_FLAGS_HAVE_VXC) /= 0) aux = aux // "V-"
          if (iand(iflags,XC_FLAGS_HAVE_FXC) /= 0) aux = aux // "F-"
          if (iand(iflags,XC_FLAGS_HAVE_KXC) /= 0) aux = aux // "K-"
          if (iand(iflags,XC_FLAGS_HAVE_LXC) /= 0) aux = aux // "L-"
          if (iand(iflags,XC_FLAGS_1D) /= 0) aux = aux // "1D-"
          if (iand(iflags,XC_FLAGS_2D) /= 0) aux = aux // "2D-"
          if (iand(iflags,XC_FLAGS_3D) /= 0) aux = aux // "3D-"
          if (iand(iflags,XC_FLAGS_HYB_CAM) /= 0) aux = aux // "CAM-"
          if (iand(iflags,XC_FLAGS_HYB_CAMY) /= 0) aux = aux // "CAMY-"
          if (iand(iflags,XC_FLAGS_VV10) /= 0) aux = aux // "VV10-"
          if (iand(iflags,XC_FLAGS_HYB_LC) /= 0) aux = aux // "LC-"
          if (iand(iflags,XC_FLAGS_HYB_LCY) /= 0) aux = aux // "LCY-"
          if (iand(iflags,XC_FLAGS_STABLE) /= 0) aux = aux // "STB-"
          if (iand(iflags,XC_FLAGS_DEVELOPMENT) /= 0) aux = aux // "DEV-"
          if (iand(iflags,XC_FLAGS_NEEDS_LAPLACIAN) /= 0) aux = aux // "LAP-"
          str = str // " " // string(aux,25,ioj_left)
       end if

       if (doname) then
          str = str // " " // trim(longname)
       end if

       if (doref) then
          iref = 0
          str = str // " ["
          do while (iref >= 0)
             ref = xc_f03_func_info_get_references(info,iref)
             if (iref < 0) exit
             longname = xc_f03_func_reference_get_ref(ref)
             if (iref == 1) then
                str = str // trim(longname)
             else
                str = str // "|" // trim(longname)
             end if
             if (iref == XC_MAX_REFERENCES-1) exit
          end do
          str = str // "]"
       end if

       call xc_f03_func_end(p)
       write (uout,'(A)') trim(str)
     end do
     write (uout,*)
#endif

  end subroutine listlibxc

  !> Tokenize the expresion expr and return a list of tokens in
  !> toklist. sptr C pointer to the calling system.
  module subroutine pretokenize(expr,toklist,errmsg,sptr)
    use systemmod, only: system
    use iso_c_binding, only: c_ptr, c_f_pointer
    character(*), intent(in) :: expr
    type(token), intent(inout), allocatable :: toklist(:)
    character(len=:), allocatable, intent(inout) :: errmsg
    type(c_ptr), intent(in), optional :: sptr

    type(system), pointer :: syl
    integer :: ntok, lp
    logical :: ok
    type(token), allocatable :: aux(:)

    ! initialize
    errmsg = ""

    ! recover the system pointer
    syl => null()
    if (present(sptr)) then
       call c_f_pointer(sptr,syl)
       if (.not.syl%isinit) then
          errmsg = 'system not initialized'
          return
       end if
    end if

    ! tokenize the expression
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,syl,errmsg)
    if (.not.ok) then
       errmsg = errmsg // ' in expression: ' // trim(expr)
       return
    end if

    ! reallocate
    if (ntok /= size(toklist,1)) then
       allocate(aux(ntok))
       aux = toklist(1:ntok)
       call move_alloc(aux,toklist)
    end if

  end subroutine pretokenize

  !xx! private procedures

  !> Given an expression in string expr starting at lpexit, parse all
  !> tokens for the arithmetic evaluation. Return the tokens in
  !> toklist and the number of tokens in ntok, and advance the string
  !> pointer lpexit. syl = calling system. This routine is
  !> thread-safe.
  function tokenize(expr,ntok,toklist,lpexit,syl,errmsg)
    use systemmod, only: system
    use tools_io, only: isinteger, string
    use param, only: vh, tab
    logical :: tokenize
    character(*), intent(in) :: expr
    integer, intent(out) :: ntok
    type(token), intent(inout), allocatable :: toklist(:)
    integer, intent(inout) :: lpexit
    type(system), intent(inout), optional :: syl
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: lp, ll, lpt
    character(len=:), allocatable :: str
    character*10 :: fder
    logical :: ok, wasop, inchem, iskey
    real*8 :: a
    integer :: c, npar, id, lp2, cchem

    errmsg = ""
    lpt = lpexit

    ! initialize
    tokenize = .true.
    inchem = .false.
    npar = 0
    cchem = 0

    ! allocate the token list
    if (allocated(toklist)) deallocate(toklist)
    ntok = 0
    allocate(toklist(10))

    ! length of the expression and initialization
    lp = lpexit
    ll = len_trim(expr)

    ! parse
    wasop = .true.
    main: do while (lp <= ll)
       lpexit = lp
       ! skip blanks and quotes
       do while (expr(lp:lp) == ' '.or.expr(lp:lp) == tab.or.expr(lp:lp)==new_line('a').or.&
          expr(lp:lp) == '"'.or.expr(lp:lp) == "'")
          lp = lp + 1
          if (lp > ll) exit main
       enddo
       lpt = lp
       if (isnumber(a,expr,lp)) then
          ! a number (without sign); if it is the field of a chemical
          ! function, keep it in sval for fields_in_eval, as addfieldnum
          if (inchem .and. ntok >= 2) then
             if (toklist(ntok)%type == token_lpar .and. toklist(ntok-1)%type == token_fun .and.&
                toklist(ntok-1)%ival == cchem) then
                call addtok(token_num,fval=a,sval=string(nint(a)))
             else
                call addtok(token_num,fval=a)
             end if
          else
             call addtok(token_num,fval=a)
          end if
          wasop = .false.
       elseif (isfunction(c,expr,lp,wasop)) then
          ! a function
          call addtok(token_fun,ival=c)
          wasop = .false.
          if (istype(c,'chemfunction')) then
             inchem = .true.
             npar = 0
             cchem = c
          end if
       elseif (isoperator(c,expr,lp)) then
          ! a binary operator
          call addtok(token_op,ival=c)
          wasop = .true.
       elseif (expr(lp:lp) == '(') then
          ! a left parenthesis
          call addtok(token_lpar)
          wasop = .true.
          lp = lp + 1
          npar = npar + 1
       elseif (expr(lp:lp) == ')') then
          ! a right parenthesis
          if (ntok > 0) then
             if (toklist(ntok)%type == token_comma) then
                errmsg = "empty argument before ')'"
                goto 999
             end if
          end if
          call addtok(token_rpar)
          wasop = .false.
          lp = lp + 1
          npar = npar - 1
          if (npar <= 0 .and. inchem) inchem = .false.
       elseif (expr(lp:lp) == ',') then
          ! a comma
          if (ntok > 0) then
             if (toklist(ntok)%type == token_comma .or. toklist(ntok)%type == token_lpar) then
                errmsg = "empty argument before ','"
                goto 999
             end if
          end if
          call addtok(token_comma)
          wasop = .true.
          lp = lp + 1
       elseif (expr(lp:lp) == '@') then
          ! a structural variable
          lp = lp + 1
          ok = isidentifier(str,expr,lp,fder)
          if (.not.ok) then
             errmsg = "a structural variable name must follow @, e.g. @dnuc"
             goto 999
          end if
          ok = isstructvar(str,c,fder)
          if (.not.ok) then
             if (arith_find(str,arithkind_svar) > 0) then
                errmsg = "@" // trim(str) // " does not take a selector (:" // trim(fder) // ")"
             else
                errmsg = "unknown structural variable: @" // trim(str)
             end if
             goto 999
          end if
          call addtok(token_structvar,ival=c,fder=fder)
          wasop=.false.

       elseif (expr(lp:lp) == '$') then
          ! a field read the field identifier and the : identifier
          lp = lp + 1
          ok = isidentifier(str,expr,lp,fder)
          if (.not.ok) then
             errmsg = "a field number or name must follow $, e.g. $1"
             goto 999
          end if

          ! an integer or a known key
          lp2 = 1
          ok = isinteger(id,str,lp2)
          iskey = isfieldkey(str)
          if (.not.ok .and. .not.iskey) then
             errmsg = "unknown field: $" // trim(str) // loaded_fields(syl)
             goto 999
          end if
          if (.not.inchem) then
             ! normal interpretation
             call addtok(token_field,sval=str,fder=fder)
          elseif (iskey) then
             ! inside a chemical function
             call addfieldnum(str)
          else
             if (syl%goodfield(key=str)) then
                errmsg = arith_usage(cchem,"") // " takes a field number or name without $, e.g. " //&
                   arith_usage(cchem,trim(str))
             else
                errmsg = "unknown field in " // arith_usage(cchem,"") // ": " // trim(str) // loaded_fields(syl)
             end if
             goto 999
          end if
          wasop = .false.
       elseif (isidentifier(str,expr,lp)) then
          if (.not.inchem) then
             ! a constant or a variable
             ok = vh%iskey(trim(str))
             if (.not.ok) then
                if (isfieldkey(str)) then
                   errmsg = "'" // trim(str) // "' is a field: write $" // trim(str)
                elseif (arith_find(str,arithkind_svar) > 0) then
                   errmsg = "'" // trim(str) // "' is a structural variable: write @" // trim(str)
                else
                   errmsg = "unknown name: " // trim(str)
                end if
                goto 999
             end if
             call addtok(token_num,fval=vh%get(trim(str),a))
          else
             ! inside a chemical function -> field identifier
             if (.not.isfieldkey(str)) then
                errmsg = "unknown field in " // arith_usage(cchem,"") // ": " // trim(str) // loaded_fields(syl)
                goto 999
             end if
             call addfieldnum(str)
          end if
          wasop = .false.
       elseif (expr(lp:lp) == '=') then
          errmsg = "use == to compare two values"
          goto 999
       else
          errmsg = "unexpected character: " // expr(lp:lp)
          goto 999
       end if
    end do main

    lpexit = lp
    tokenize = .true.
    return

999 continue ! could not tokenize
    errmsg = errmsg // " (column " // string(lpt) // ")"
    tokenize = .false.
    return
  contains
    subroutine addtok(type,ival,fval,sval,fder)
      integer, intent(in) :: type
      integer, intent(in), optional :: ival
      real*8, intent(in), optional :: fval
      character*(*), intent(in), optional :: sval
      character*(*), intent(in), optional :: fder

      type(token), allocatable :: auxtoklist(:)

      ntok = ntok + 1
      if (ntok > size(toklist)) then
         allocate(auxtoklist(2*ntok))
         auxtoklist(1:size(toklist)) = toklist
         call move_alloc(auxtoklist,toklist)
      end if
      toklist(ntok)%type = type
      if (present(ival)) toklist(ntok)%ival = ival
      if (present(fval)) toklist(ntok)%fval = fval
      if (present(sval)) then
         toklist(ntok)%sval = sval
      else
         toklist(ntok)%sval = ""
      end if
      if (present(fder)) toklist(ntok)%fder = fder

    end subroutine addtok

    !> Is name a field name (key) in the system?
    function isfieldkey(name)
      character*(*), intent(in) :: name
      logical :: isfieldkey

      isfieldkey = .false.
      if (present(syl)) isfieldkey = syl%fh%iskey(trim(name))

    end function isfieldkey

    !> Add the id of the field with this name as a number (the
    !> argument of a chemical function). The name is kept in sval
    !> for fields_in_eval.
    subroutine addfieldnum(name)
      character*(*), intent(in) :: name

      integer :: idf

      idf = syl%fh%get(trim(name),1)
      call addtok(token_num,fval=real(idf,8),sval=trim(name))

    end subroutine addfieldnum
  end function tokenize

  !> Return the list of fields loaded in system syl, to append to
  !> error messages ("; loaded fields: $0, $1"), or an empty string
  !> if syl is not present. This routine is thread-safe.
  function loaded_fields(syl) result(str)
    use systemmod, only: system
    use tools_io, only: string
    type(system), intent(in), optional :: syl
    character(len=:), allocatable :: str

    integer :: i

    str = ""
    if (.not.present(syl)) return
    if (allocated(syl%f)) then
       do i = 0, syl%nf
          if (.not.syl%f(i)%isinit) cycle
          if (len(str) > 0) str = str // ", "
          str = str // "$" // string(i)
       end do
    end if
    if (len(str) == 0) then
       str = "; no fields are loaded"
    else
       str = "; loaded fields: " // str
    end if

  end function loaded_fields

  !> Return the name of the function with code c (fun_*) as it is
  !> written in an expression, with its arguments: gtf(id), or with
  !> argument arg if present: gtf(1), gtf(). This routine is
  !> thread-safe.
  function arith_usage(c,arg) result(str)
    integer, intent(in) :: c
    character*(*), intent(in), optional :: arg
    character(len=:), allocatable :: str

    integer :: i

    str = "the function"
    do i = 1, size(arith_names,1)
       if (arith_names(i)%kind == arithkind_svar) cycle
       if (arith_names(i)%code /= c) cycle
       if (present(arg)) then
          str = trim(arith_names(i)%name) // "(" // arg // ")"
       else
          str = trim(arith_names(i)%name) // trim(arith_names(i)%args)
       end if
       return
    end do

  end function arith_usage

  !> Return the index in the table of names (arith_names) of the
  !> function with code c (fun_*), or zero if c is not a function
  !> (operators, parentheses). This routine is thread-safe.
  function arith_code(c) result(idx)
    integer, intent(in) :: c
    integer :: idx

    do idx = 1, size(arith_names,1)
       if (arith_names(idx)%kind == arithkind_svar) cycle
       if (arith_names(idx)%code == c) return
    end do
    idx = 0

  end function arith_code

  !> Return the number of arguments of the function with index idx
  !> in the table of names (arith_names), from the argument list in
  !> the table: (x) = 1, (x,y) = 2,... Returns -1 for xc(), which
  !> takes a number that depends on the functional. This routine is
  !> thread-safe.
  function arith_nargs(idx) result(n)
    integer, intent(in) :: idx
    integer :: n

    integer :: i

    if (arith_names(idx)%code == fun_xc) then
       n = -1
       return
    end if
    n = 1
    do i = 1, len_trim(arith_names(idx)%args)
       if (arith_names(idx)%args(i:i) == ',') n = n + 1
    end do

  end function arith_nargs

  !> Check that the function with code c (fun_*) was given na
  !> arguments, and return an error message if it was not, or an
  !> empty string if it was. xc() needs at least two here; its exact
  !> number of arguments depends on the functional and is checked in
  !> pop. This routine is thread-safe.
  function check_nargs(c,na) result(errmsg)
    integer, intent(in) :: c, na
    character(len=:), allocatable :: errmsg

    integer :: idx, want

    errmsg = ""
    idx = arith_code(c)
    if (idx == 0) return
    want = arith_nargs(idx)
    if (want < 0) then
       if (na < 2) errmsg = "missing arguments in " // arith_usage(c)
    elseif (na > want) then
       errmsg = "too many arguments in " // arith_usage(c)
    elseif (na < want) then
       if (arith_names(idx)%kind == arithkind_chem .and. want == 1) then
          errmsg = "missing field in " // arith_usage(c) // ", e.g. " // arith_usage(c,"1")
       elseif (want == 1) then
          errmsg = "missing argument in " // arith_usage(c)
       else
          errmsg = "missing arguments in " // arith_usage(c)
       end if
    end if

  end function check_nargs

  !> Return an operator precedence. This routine is thread-safe.
  function iprec(c)
    integer :: iprec
    integer, intent(in) :: c

    iprec = 0
    if (c == fun_or) then ! or
       iprec = 1
    else if (c == fun_and) then ! and
       iprec = 2
    else if (c == fun_neq .or. c == fun_equal .or. c == fun_leq .or. &
       c == fun_geq .or. c == fun_less .or. c == fun_great) then ! relational
       iprec = 3
    else if (c == fun_plus .or. c == fun_minus) then ! add, sub
       iprec = 4
    elseif (c == fun_prod .or. c == fun_div .or. c == fun_modulo) then ! mult, div
       iprec = 5
    elseif (c == fun_power) then ! exp
       iprec = 6
    elseif (c == fun_uminus .or. c == fun_uplus) then ! unary plus and minus
       iprec = 7
    elseif (c == fun_openpar .or. c == fun_closepar) then ! parentheses
       iprec = 0
    end if

  end function iprec

  !> Return 1 for right-associative operator and -1 for
  !> left-associative.  This routine is thread-safe.
  function iassoc(c)
    integer :: iassoc
    integer, intent(in) :: c

    iassoc = 0
    if (c == fun_plus .or. c == fun_minus .or. c == fun_prod .or. &
       c == fun_div .or. c == fun_modulo .or. c == fun_or .or. &
       c == fun_and .or. c == fun_equal .or. c == fun_leq .or. &
       c == fun_geq .or. c == fun_less .or. c == fun_great) then
       iassoc = -1
    elseif (c == fun_power .or. c == fun_openpar .or. c == fun_closepar .or.&
       c == fun_uplus .or. c == fun_uminus) then
       iassoc = 1
    end if

  end function iassoc

  !> Return the type of stack element (operator, function, unary,
  !> binary). If iwantarg is present and c is a chemical function,
  !> return the number of arguments the chemical function wants;
  !> otherwise, return zero. The functions and their arguments come
  !> from the table of names (arith_names); xc() is not a chemical
  !> function here because it takes field values, not a field id.
  !> This routine is thread-safe.
  function istype(c,type,iwantarg)
    integer, intent(in) :: c
    character*(*), intent(in) :: type
    integer, intent(inout), optional :: iwantarg
    logical :: istype

    integer :: idx

    if (present(iwantarg)) iwantarg = 0
    istype = .false.

    ! operators and signs
    select case (type)
    case ('operator','binary')
       istype = &
          c == fun_power .or. c == fun_leq .or. c == fun_geq .or.&
          c == fun_equal .or. c == fun_neq .or. c == fun_and .or.&
          c == fun_or .or. c == fun_plus .or. c == fun_minus .or.&
          c == fun_prod .or. c == fun_div .or. c == fun_modulo .or.&
          c == fun_great .or. c == fun_less
    case ('unary')
       istype = c == fun_uplus .or. c == fun_uminus
    end select
    if (istype .or. type == 'operator') return

    ! functions, from the table
    idx = arith_code(c)
    if (idx == 0) return
    select case (type)
    case ('function')
       istype = .true.
    case ('chemfunction')
       istype = arith_names(idx)%kind == arithkind_chem .and. c /= fun_xc
       if (istype .and. present(iwantarg)) iwantarg = arith_nargs(idx)
    case ('binary')
       istype = arith_names(idx)%kind == arithkind_math .and. arith_nargs(idx) == 2
    case ('unary')
       istype = arith_names(idx)%kind == arithkind_math .and. arith_nargs(idx) == 1
    end select

  endfunction istype

  !> If this identifier corresponds to a structural variable, return
  !> .true. and the integer ID for the variable in c. Otherwise,
  !> return .false.. This routine is thread-safe.
  function isstructvar(fid,c,fder)
    character*(*), intent(in) :: fid
    integer, intent(out) :: c
    character*(*), intent(in) :: fder
    logical :: isstructvar

    integer :: i

    isstructvar = .false.
    c = 0
    i = arith_find(fid,arithkind_svar)
    if (i == 0) return
    c = arith_names(i)%code

    ! only some take the atom selector (:n)
    isstructvar = (len_trim(fder) == 0 .or. arith_names(i)%args == ":n")

  end function isstructvar

  !> Find name in the table of names (arith_names) of the given kind
  !> (arithkind_*; zero = any function, math or chemical),
  !> case-insensitive. Returns its index in the table,
  !> or zero if not found. This routine is thread-safe.
  function arith_find(name,kind) result(idx)
    use tools_io, only: lower
    character*(*), intent(in) :: name
    integer, intent(in) :: kind
    integer :: idx

    character(len=len(arith_names(1)%name)) :: lname
    integer :: n

    ! can run once per name token and evaluation point: keep it short
    idx = 0
    n = len_trim(name)
    if (n > len(lname)) return
    lname = lower(name(1:n))
    do idx = 1, size(arith_names,1)
       if (kind == 0) then
          if (arith_names(idx)%kind == arithkind_svar) cycle
       elseif (arith_names(idx)%kind /= kind) then
          cycle
       end if
       if (lname == arith_names(idx)%name) return
    end do
    idx = 0

  end function arith_find

  !> Evaluate field with identifier fid and field flag fder at point
  !> x0 and return the value. syl = calling system. If periodic is
  !> true, evaluate the field under pbc. Returns non-empty errmsg if
  !> the evaluation failed.
  recursive function fieldeval(fid,fder,errmsg,x0,syl,periodic)
    use systemmod, only: system
    use tools_math, only: det3sym
    use tools_io, only: isinteger, lower
    use types, only: scalar_value, field_evaluation_avail, fieldeval_category_f,&
       fieldeval_category_spin, fieldeval_category_fder1, fieldeval_category_fder2,&
       fieldeval_category_mo, id_mo_homo, id_mo_lumo, id_mo_ahomo, id_mo_alumo,&
       id_mo_bhomo, id_mo_blumo, id_mo_a, id_mo_b, id_mo_id
    character*(*), intent(in) :: fid
    character*(*), intent(in) :: fder
    character(len=:), allocatable, intent(inout) :: errmsg
    real*8, intent(in), optional :: x0(3)
    type(system), intent(inout), optional :: syl
    logical, intent(in), optional :: periodic
    real*8 :: fieldeval

    integer :: idx
    type(scalar_value) :: res
    character*10 :: fderl
    type(field_evaluation_avail) :: request
    logical :: ok

    errmsg = ""

    ! recover the system pointer
    if (present(syl)) then
       if (.not.syl%isinit) then
          errmsg = fnam() // ': system not initialized'
          return
       end if
    end if

    fieldeval = 0d0
    if (present(x0).and.present(syl)) then
       if (syl%goodfield(key=fid,idout=idx)) then
          ! request from the field
          call request%clear()
          fderl = lower(adjustl(fder))
          select case (fderl)
          case ("","v","c")
             request%avail(fieldeval_category_f) = .true.
          case("up","dn","sp")
             request%avail(fieldeval_category_f) = .true.
             request%avail(fieldeval_category_spin) = .true.
          case ("x","y","z","g")
             request%avail(fieldeval_category_fder1) = .true.
          case ("xx","xy","xz","yx","yy","yz","zx","zy","zz","l","lv","lc","h1","h2","h3","hd","s","r")
             request%avail(fieldeval_category_fder2) = .true.
          case default
             ! this must be a molecular orbital
             request%avail(fieldeval_category_mo) = .true.
             ok = .true.
             if (fderl == "homo") then
                request%moini = id_mo_homo
             elseif (fderl == "lumo") then
                request%moini = id_mo_lumo
             elseif (fderl == "ahomo") then
                request%moini = id_mo_ahomo
             elseif (fderl == "alumo") then
                request%moini = id_mo_alumo
             elseif (fderl == "bhomo") then
                request%moini = id_mo_bhomo
             elseif (fderl == "blumo") then
                request%moini = id_mo_blumo
             elseif (fderl(1:1) == "a") then
                request%moini = id_mo_a
                ok = isinteger(request%moend,fderl(2:))
             elseif (fderl(1:1) == "b") then
                request%moini = id_mo_b
                ok = isinteger(request%moend,fderl(2:))
             else
                request%moini = id_mo_id
                ok = isinteger(request%moend,fderl)
             end if
             if (.not.ok) then
                errmsg = fnam() // ': unknown field modifier'
                return
             end if
          end select

          ! run the field evaluation
          call syl%f(idx)%grd(x0,request,res,periodic=periodic)
          if (.not.res%valid) then
             errmsg = fnam() // ': invalid evaluation result'
             fieldeval = 0d0
             return
          end if
          if (.not.res%satisfied) then
             errmsg = fnam() //&
                ': cannot compute the property with this field type'
             fieldeval = 0d0
             return
          end if

          ! Unpack the result
          select case (trim(fderl))
          case ("")
             fieldeval = res%f
          case ("v")
             fieldeval = res%fval
          case ("c")
             fieldeval = res%f - res%fval
          case ("up")
             fieldeval = res%fspin(1)
          case ("dn")
             fieldeval = res%fspin(2)
          case ("sp")
             fieldeval = res%fspin(1) - res%fspin(2)
          case ("x")
             fieldeval = res%gf(1)
          case ("y")
             fieldeval = res%gf(2)
          case ("z")
             fieldeval = res%gf(3)
          case ("xx")
             fieldeval = res%hf(1,1)
          case ("xy","yx")
             fieldeval = res%hf(1,2)
          case ("xz","zx")
             fieldeval = res%hf(1,3)
          case ("yy")
             fieldeval = res%hf(2,2)
          case ("yz","zy")
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
          case ("h1")
             fieldeval = res%hfeval(1)
          case ("h2")
             fieldeval = res%hfeval(2)
          case ("h3")
             fieldeval = res%hfeval(3)
          case ("hd")
             fieldeval = det3sym(res%hf)
          case ("r")
             fieldeval = res%r
          case ("s")
             fieldeval = res%s
          case default
             fieldeval = res%fspc
          end select
       else
          errmsg = fnam() // ': unknown field' // loaded_fields(syl)
          return
       end if
    else
       errmsg = fnam() // ': no point for field evaluation'
       return
    end if

  contains
    !> The field and modifier as written in the expression ($1:x).
    function fnam() result(str)
      character(len=:), allocatable :: str

      str = "$" // trim(adjustl(fid))
      if (len_trim(fder) > 0) str = str // ":" // trim(adjustl(fder))

    end function fnam
  end function fieldeval

  !> Evaluate a structural variable at point x0 using the crystal/molecular
  !> structure in syl. periodic=.true. if the system is assumed periodic.
  function structvareval(svar,fder,errmsg,x0,syl,periodic) result(q)
    use systemmod, only: system
    use grid1mod, only: agrid
    use global, only: dunit0, iunit
    use tools_io, only: string, isinteger, lower
    use param, only: icrd_cart
    real*8 :: q
    integer, intent(in) :: svar
    character*(*), intent(in) :: fder
    character(len=:), allocatable, intent(inout) :: errmsg
    real*8, intent(in), optional :: x0(3)
    type(system), intent(inout), optional :: syl
    logical, intent(in), optional :: periodic

    integer :: nid, nid0, lvec(3), iz, lp
    real*8 :: dist, x(3), rrho1, rrho2
    logical :: ok

    errmsg = ""

    ! recover the system pointer
    if (present(syl)) then
       if (.not.syl%isinit) then
          errmsg = 'evaluating structural variable but system not initialized'
          return
       end if
       if (.not.syl%c%isinit) then
          errmsg = 'evaluating structural variable but structure not initialized'
          return
       end if
    end if
    if (.not.present(x0)) then
       errmsg = 'evaluating structural variable without point'
       return
    end if

    select case (svar)
    case(svar_xc)
       q = x0(1)
    case(svar_yc)
       q = x0(2)
    case(svar_zc)
       q = x0(3)
    case(svar_xm)
       q = (x0(1) + syl%c%molx0(1)) * dunit0(iunit)
    case(svar_ym)
       q = (x0(2) + syl%c%molx0(2)) * dunit0(iunit)
    case(svar_zm)
       q = (x0(3) + syl%c%molx0(3)) * dunit0(iunit)
    case(svar_x,svar_y,svar_z)
       if (syl%c%ismolecule) then
          x = (x0 + syl%c%molx0) * dunit0(iunit)
       else
          x = syl%c%c2x(x0)
       end if
       if (svar == svar_x) then
          q = x(1)
       elseif (svar == svar_y) then
          q = x(2)
       else
          q = x(3)
       end if
    case(svar_xx,svar_yx,svar_zx)
       x = syl%c%c2x(x0)
       if (svar == svar_xx) then
          q = x(1)
       elseif (svar == svar_yx) then
          q = x(2)
       else
          q = x(3)
       end if
    case(svar_xxr,svar_yxr,svar_zxr)
       x = syl%c%c2xr(x0)
       if (svar == svar_xxr) then
          q = x(1)
       elseif (svar == svar_yxr) then
          q = x(2)
       else
          q = x(3)
       end if
    case(svar_dnuc,svar_xnucx,svar_ynucx,svar_znucx,svar_xnucc,svar_ynucc,svar_znucc,&
         svar_idnuc,svar_nidnuc,svar_rho0nuc,svar_spcnuc,svar_zatnuc)
       q = 0d0
       if (len_trim(fder) > 0) then
          lp = 1
          ok = isinteger(nid,fder,lp)
          if (.not.ok) then
             errmsg = 'wrong selector in rho0nuc structural variable'
             return
          end if
          if (nid < 1 .or. nid > syl%c%ncel) then
             errmsg = 'atom ID in rho0nuc structural variable out of range'
             return
          end if
          call syl%c%nearest_atom(x0,icrd_cart,nid0,dist,lvec=lvec,id0=nid)
          if (nid /= nid0) return ! fixme: the atom was too far
       else
          call syl%c%nearest_atom(x0,icrd_cart,nid,dist,lvec=lvec)
       end if

       if (nid == 0) then
          ! fixme: the atom was too far
          q = 0d0
       else if (svar == svar_dnuc) then
          q = dist * dunit0(iunit)
       elseif (svar == svar_xnucx) then
          q = syl%c%atcel(nid)%x(1)
       elseif (svar == svar_ynucx) then
          q = syl%c%atcel(nid)%x(2)
       elseif (svar == svar_znucx) then
          q = syl%c%atcel(nid)%x(3)
       elseif (svar == svar_xnucc) then
          q = (syl%c%atcel(nid)%r(1) + syl%c%molx0(1)) * dunit0(iunit)
       elseif (svar == svar_ynucc) then
          q = (syl%c%atcel(nid)%r(2) + syl%c%molx0(2)) * dunit0(iunit)
       elseif (svar == svar_znucc) then
          q = (syl%c%atcel(nid)%r(3) + syl%c%molx0(3)) * dunit0(iunit)
       elseif (svar == svar_idnuc) then
          q = nid
       elseif (svar == svar_nidnuc) then
          q = syl%c%atcel(nid)%idx
       elseif (svar == svar_rho0nuc) then
          iz = syl%c%spc(syl%c%atcel(nid)%is)%z
          q = 0d0
          if (iz > 0) then
             if (agrid(iz)%isinit) then
                if (dist <= agrid(iz)%rmax) then
                   call agrid(iz)%interp(dist,q,rrho1,rrho2)
                end if
             end if
          end if
       elseif (svar == svar_spcnuc) then
          q = syl%c%atcel(nid)%is
       elseif (svar == svar_zatnuc) then
          q = syl%c%spc(syl%c%atcel(nid)%is)%z
       end if
    case default
       errmsg = 'evaluating structural variable: unknown variable'
       return
    end select

  end function structvareval

  !> Read an unsigned number or return false and leave lp unchanged
  !> This routine is thread-safe.
  function isnumber(rval,expr,lp)
    use tools_io, only: isdigit
    logical :: isnumber
    character*(*), intent(in) :: expr !< Input string
    integer, intent(inout) :: lp !< Pointer to current position on string
    real*8, intent(out) :: rval !< Real value read

    integer :: i, ll

    ll = len_trim(expr)
    isnumber = .false.
    rval = 0d0
    if (lp < 1 .or. lp > len(expr)) return
    i = lp
    do while (isdigit(expr(i:i)))
       i = i + 1
       if (i > ll) goto 999
    enddo
    if (expr(i:i) == '.') then
       i = i + 1
       if (i > ll) goto 999
    end if
    do while (isdigit(expr(i:i)))
       i = i + 1
       if (i > ll) goto 999
    enddo
    if (i == lp) return

    ! exponent
    if (expr(i:i)=='e' .or. expr(i:i)=='E' .or. expr(i:i)=='d' .or. expr(i:i)=='D'.or.&
       expr(i:i)=='q' .or. expr(i:i)=='Q') then
       i = i + 1
       if (i > ll) return
       if (expr(i:i)=='-' .or. expr(i:i)=='+') then
          i = i + 1
          if (i > ll) return
       end if
       if (.not.isdigit(expr(i:i))) return
       do while (isdigit(expr(i:i)))
          i = i + 1
          if (i > ll) goto 999
       enddo
    end if

    if (lp == i-1) then
       if (expr(lp:i-1) == ".") return
    end if

999 continue
    isnumber=.true.
    read(expr(lp:i-1),*) rval
    lp = i

  end function isnumber

  !> Read a binary operator or return false and leave lp unchanged
  !> This routine is thread-safe.
  function isoperator(c,expr,lp)
    use param, only: ampersand
    logical :: isoperator
    character*(*), intent(in) :: expr
    integer, intent(inout) :: lp
    integer, intent(out) :: c

    logical :: ok1, ok2
    integer :: ll

    ll = len_trim(expr)
    isoperator = .false.
    if (lp > len(expr)) return
    ok1 = expr(lp:lp) == '+'.or.expr(lp:lp) == '-'.or.expr(lp:lp) == '*'.or.&
          expr(lp:lp) == '/'.or.expr(lp:lp) == '^'.or.expr(lp:lp) == '%'.or.&
          expr(lp:lp) == '<'.or.expr(lp:lp) == '>'
    if (lp < ll) then
       ok2 = expr(lp:lp+1) == '**'.or.expr(lp:lp+1) == '=='.or.&
             expr(lp:lp+1) == '<='.or.expr(lp:lp+1) == '>='.or.&
             expr(lp:lp+1) == '!='.or.expr(lp:lp+1) == '||'.or.&
             expr(lp:lp+1) == ampersand//ampersand
    else
       ok2 = .false.
    endif

    if (ok1 .or. ok2) then
       if (ok2) then
          if (expr(lp:lp+1) == "**") then
             c = fun_power
             lp = lp + 1
          elseif (expr(lp:lp+1) == "<=") then
             c = fun_leq
             lp = lp + 1
          elseif (expr(lp:lp+1) == ">=") then
             c = fun_geq
             lp = lp + 1
          elseif (expr(lp:lp+1) == "==") then
             c = fun_equal
             lp = lp + 1
          elseif (expr(lp:lp+1) == "!=") then
             c = fun_neq
             lp = lp + 1
          elseif (expr(lp:lp+1) == ampersand//ampersand) then
             c = fun_and
             lp = lp + 1
          elseif (expr(lp:lp+1) == "||") then
             c = fun_or
             lp = lp + 1
          else
             c = fun_unknown
          end if
       else
          if (expr(lp:lp) == "+") then
             c = fun_plus
          elseif (expr(lp:lp) == "-") then
             c = fun_minus
          elseif (expr(lp:lp) == "*") then
             c = fun_prod
          elseif (expr(lp:lp) == "/") then
             c = fun_div
          elseif (expr(lp:lp) == "%") then
             c = fun_modulo
          elseif (expr(lp:lp) == ">") then
             c = fun_great
          elseif (expr(lp:lp) == "<") then
             c = fun_less
          elseif (expr(lp:lp) == "^") then
             c = fun_power
          else
             c = fun_unknown
          end if
       endif
       lp = lp + 1
       isoperator = .true.
    end if

  end function isoperator

  !> Read a unary operator (function) or return false and leave lp
  !> unchanged This routine is thread-safe.
  function isfunction(c,expr,lp,wasop)
    logical :: isfunction
    character*(*), intent(in) :: expr
    integer, intent(inout) :: lp
    integer, intent(out) :: c
    logical, intent(in) :: wasop

    integer :: lpo, ll, i
    character*(len(expr)) :: word

    isfunction = .false.
    lpo = lp
    c = 0
    ll = len_trim(expr)
    if (lp > len(expr)) return
    if (expr(lp:lp) == '+' .and. wasop) then
       ! unary +
       c = fun_uplus
       lp = lp + 1
    elseif (expr(lp:lp) == '-' .and. wasop) then
       ! unary -
       c = fun_uminus
       lp = lp + 1
    elseif (expr(lp:lp) >= 'a' .and. expr(lp:lp)<='z' .or. &
       expr(lp:lp) >= 'A' .and. expr(lp:lp)<='Z') then
       ! function
       word = ""
       do while (expr(lp:lp) >= 'a' .and. expr(lp:lp)<='z' .or. &
          expr(lp:lp) >= 'A' .and. expr(lp:lp)<='Z' .or. &
          expr(lp:lp) >= '0' .and. expr(lp:lp)<='9'.or. &
          expr(lp:lp) == '_')
          word = trim(word) // expr(lp:lp)
          lp = lp + 1
          if (lp > ll) exit
       end do
       i = arith_find(word,0)
       if (i == 0) then
          lp = lpo
          return
       end if
       c = arith_names(i)%code
    else
       lp = lpo
       return
    end if

    isfunction = .true.

  end function isfunction

  !> Read an identifier (field, variable, constant) or return false
  !> and leave lp unchanged. If fder is present, try to find the field
  !> derivative selector after the identifier (:xx).  This routine is
  !> thread-safe.
  function isidentifier(id,expr,lp,fder)
    use tools_io, only: isletter, isdigit
    logical :: isidentifier
    character(len=:), allocatable, intent(out) :: id
    character*(*), intent(in) :: expr
    integer, intent(inout) :: lp
    character*(*), intent(out), optional :: fder

    character*(len(expr)) :: word
    integer :: lpo, i, ll

    isidentifier = .false.
    lpo = lp
    ll = len_trim(expr)
    if (lp < 1 .or. lp > len(expr)) return
    if (isletter(expr(lp:lp)) .or. isdigit(expr(lp:lp)) .or. expr(lp:lp) == '_') then
       ! read the identifier
       word = ""
       do while (isletter(expr(lp:lp)) .or. isdigit(expr(lp:lp)) .or. expr(lp:lp) == '_')
          word = trim(word) // expr(lp:lp)
          lp = lp + 1
          if (lp > ll) exit
       end do
       id = word
    else
       lp = lpo
       return
    end if
    isidentifier = .true.

    if (present(fder)) then
       fder = ""
       if (lp <= ll) then
          if (expr(lp:lp) == ":") then
             lp = lp + 1
             i = lp
             do while (i <= ll)
                if (.not.(isletter(expr(i:i)).or.isdigit(expr(i:i)))) exit
                i = i + 1
             enddo
             fder = expr(lp:i-1)
             lp = i
          end if
       end if
    endif

  end function isidentifier

  !> Pop from the stack and operate on the queue. q is the stack of
  !> values (with nq elements). s is the stack of operations (with ns
  !> elements). If the expression contains fields ($), use x0 as the
  !> evaluation point (Cartesian). syl = calling system.  If periodic
  !> is present and false, evaluate the expression at x0 considering
  !> the field as non-periodic. Return non-zero errmsg if an error was
  !> found. This routine is thread-safe.
  subroutine pop(q,nq,s,ns,x0,syl,periodic,errmsg,nargs)
    use systemmod, only: system
#ifdef HAVE_LIBXC
    use iso_c_binding, only: c_size_t
    use xc_f03_lib_m
    use tools_io, only: string
#endif
    real*8, intent(inout) :: q(:)
    integer, intent(inout) :: s(:)
    integer, intent(inout) :: nq, ns
    real*8, intent(in), optional :: x0(3)
    type(system), intent(inout), optional :: syl
    logical, intent(in), optional :: periodic
    character(len=:), allocatable, intent(inout) :: errmsg
    integer, intent(in), optional :: nargs

    integer :: ia, iwantarg
    integer :: c, i
#ifdef HAVE_LIBXC
    integer :: ierr
#endif
    real*8 :: a, b
    real*8, allocatable :: args(:)
    character*8 :: sia
#ifdef HAVE_LIBXC
    real*8 :: rho(1), grho(1), lapl(1), tau(1), zk(1)
#endif

    ! recover the system pointer
    if (present(syl)) then
       if (.not.syl%isinit) then
          errmsg = 'system not initialized'
          return
       end if
    end if

    ! pop from the stack
    errmsg = ""
    if (ns == 0) then
       errmsg = 'syntax error'
       return
    end if
    c = s(ns)
    ns = ns - 1

    ! And apply to the queue
    if (c == fun_xc) then
       ! Functional from the xc library
#ifdef HAVE_LIBXC
       if (nq < 1) then
          errmsg = "missing arguments in " // arith_usage(c)
          return
       end if
       ia = nint(q(nq))
       if (abs(q(nq) - ia) > 1d-13) then
          errmsg = "need an integer for the functional ID"
          return
       end if
       if (ia < 1 .or. ia > maxfun) then
          errmsg = "libxc functional ID out of range: " // string(ia)
          return
       end if

       ! dirty trick for syncing the threads
       if (.not.ifun(ia)%init) then
          !$omp critical (ifuninit)
          if (.not.ifun(ia)%init) then
             ! an id libxc does not know leaves the functional unusable
             call xc_f03_func_init(ifun(ia)%conf, ia, XC_UNPOLARIZED, ierr)
             if (ierr == 0) then
                ifun(ia)%id = ia
                ifun(ia)%info = xc_f03_func_get_info(ifun(ia)%conf)
                ifun(ia)%family = xc_f03_func_info_get_family(ifun(ia)%info)
                ifun(ia)%init = .true.
             end if
          end if
          !$omp end critical (ifuninit)
          if (.not.ifun(ia)%init) then
             errmsg = "unknown libxc functional ID: " // string(ia)
             return
          end if
       endif

       ! the number of arguments between the parentheses, if known;
       ! otherwise, the values in the queue
       if (.not.present(nargs)) then
          if (ifun(ia)%family == XC_FAMILY_LDA .and. nq <= 1.or.&
             ifun(ia)%family == XC_FAMILY_GGA .and. nq <= 2.or.&
             ifun(ia)%family == XC_FAMILY_MGGA .and. nq <= 4) then
             errmsg = "missing arguments in " // arith_usage(c)
             return
          endif
       else
          if (ifun(ia)%family == XC_FAMILY_LDA .and. nargs /= 2) then
             errmsg = "libxc functional " // string(ia) // " is an LDA: write xc(rho," // string(ia) // ")"
          elseif (ifun(ia)%family == XC_FAMILY_GGA .and. nargs /= 3) then
             errmsg = "libxc functional " // string(ia) // " is a GGA: write xc(rho,grad," // string(ia) // ")"
          elseif (ifun(ia)%family == XC_FAMILY_MGGA .and. nargs /= 5) then
             errmsg = "libxc functional " // string(ia) // " is a meta-GGA: write xc(rho,grad,lapl,tau," //&
                string(ia) // ")"
          end if
          if (len_trim(errmsg) > 0) return
       end if

       select case(ifun(ia)%family)
       case (XC_FAMILY_LDA)
          rho = q(nq-1)
          call xc_f03_lda_exc(ifun(ia)%conf, 1_c_size_t, rho, zk)
          nq = nq - 1
       case (XC_FAMILY_GGA)
          rho = q(nq-2)
          grho = q(nq-1)*q(nq-1)
          call xc_f03_gga_exc(ifun(ia)%conf, 1_c_size_t, rho, grho, zk)
          nq = nq - 2
       case (XC_FAMILY_MGGA)
          rho = q(nq-4)
          grho = q(nq-3)*q(nq-3)
          lapl = q(nq-2)
          tau = q(nq-1)
          call xc_f03_mgga_exc(ifun(ia)%conf, 1_c_size_t, rho, grho, lapl, tau, zk)
          nq = nq - 4
       end select
       q(nq) = zk(1) * rho(1)
#else
       errmsg = "critic2 was not compiled with libxc support"
       return
#endif
    elseif (istype(c,'binary')) then
       ! a binary operator or function
       if (nq < 2) then
          if (istype(c,'function')) then
             errmsg = "missing arguments in " // arith_usage(c)
          else
             errmsg = "missing value next to an operator"
          end if
          return
       end if
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
          if (b == 0d0) then
             errmsg = 'divided by zero'
             return
          end if
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
       case default
          errmsg = "function not implemented: " // arith_usage(c)
          return
       end select
    elseif (istype(c,'unary')) then
       ! a unary operator or function
       if (nq < 1) then
          if (istype(c,'function')) then
             errmsg = "missing argument in " // arith_usage(c)
          else
             errmsg = "missing value after a sign"
          end if
          return
       end if
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
          if (q(nq) <= 0d0) then
             errmsg = 'logarithm of a negative number'
             return
          end if
          q(nq) = log(q(nq))
       case (fun_log10)
          if (q(nq) <= 0d0) then
             errmsg = 'logarithm of a negative number'
             return
          end if
          q(nq) = log10(q(nq))
       case (fun_sin)
          q(nq) = sin(q(nq))
       case (fun_asin)
          if (q(nq) > 1d0 .or. q(nq) < -1d0) then
             errmsg = 'asin argument out of range'
             return
          end if
          q(nq) = asin(q(nq))
       case (fun_cos)
          q(nq) = cos(q(nq))
       case (fun_acos)
          if (q(nq) > 1d0 .or. q(nq) < -1d0) then
             errmsg = 'acos argument out of range'
             return
          end if
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
       case default
          errmsg = "function not implemented: " // arith_usage(c)
          return
       end select
    elseif (istype(c,'chemfunction',iwantarg)) then
       ! We need a point and the evaluator
       if (.not.present(x0).or..not.present(syl)) then
          errmsg = "a point is needed for evaluating chemical function"
          return
       endif

       ! Consume the extra arguments for this function
       if (nq < max(iwantarg,1)) then
          if (iwantarg > 1) then
             errmsg = "missing arguments in " // arith_usage(c)
          else
             errmsg = "missing field in " // arith_usage(c) // ", e.g. " // arith_usage(c,"1")
          end if
          return
       end if
       if (iwantarg > 1) then
          allocate(args(iwantarg-1))
          do i = 1, iwantarg-1
             args(i) = q(nq - (iwantarg-1) + i)
          end do
          nq = nq - (iwantarg - 1)
       end if

       ! An integer field identifier as the first argument
       ia = nint(q(nq))
       if (abs(ia - q(nq)) > 1d-13) then
          errmsg = arith_usage(c,"") // " takes a field number or name, e.g. " // arith_usage(c,"1")
          return
       end if
       write (sia,'(I8)') ia
       sia = adjustl(sia)
       if (.not.syl%goodfield(key=sia)) then
          errmsg = "unknown field in " // arith_usage(c,"") // ": " // trim(sia) // loaded_fields(syl)
          return
       end if

       ! Use the library of chemical functions
       q(nq) = chemfunction(c,sia,x0,args,syl,periodic,errmsg)
       if (len_trim(errmsg) > 0) return
       if (iwantarg > 1) &
          deallocate(args)
    elseif (c == fun_openpar) then
       errmsg = "mismatched parentheses: '(' without ')'"
       return
    else
       errmsg = 'error in expression'
       return
    end if
  end subroutine pop

  !> Pop from the stack and operate on the queue, grid version. q is
  !> the stack of values (with nq elements). s is the stack of
  !> operations (with ns elements). Return fail=.true. if an error was
  !> found. This routine is thread-safe.
  subroutine pop_grid(q,nq,s,ns,fail)
    real*8, allocatable, intent(inout) :: q(:,:,:,:)
    integer, intent(inout) :: s(:)
    integer, intent(inout) :: nq, ns
    logical, intent(out) :: fail

    integer :: c

    ! pop from the stack
    fail = .false.
    if (ns == 0) goto 999
    c = s(ns)
    ns = ns - 1

    ! And apply to the queue
    if (c == fun_xc) then
       goto 999
    elseif (istype(c,'binary')) then
       ! a binary operator or function
       if (nq < 2) goto 999
       nq = nq - 1
       select case(c)
       case (fun_plus)
          q(:,:,:,nq) = q(:,:,:,nq) + q(:,:,:,nq+1)
       case (fun_minus)
          q(:,:,:,nq) = q(:,:,:,nq) - q(:,:,:,nq+1)
       case (fun_prod)
          q(:,:,:,nq) = q(:,:,:,nq) * q(:,:,:,nq+1)
       case (fun_div)
          q(:,:,:,nq) = q(:,:,:,nq) / q(:,:,:,nq+1)
       case (fun_power)
          q(:,:,:,nq) = q(:,:,:,nq) ** q(:,:,:,nq+1)
       case (fun_modulo)
          q(:,:,:,nq) = modulo(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_atan2)
          q(:,:,:,nq) = atan2(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_min)
          q(:,:,:,nq) = min(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_max)
          q(:,:,:,nq) = max(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_great)
          where (q(:,:,:,nq) > q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_less)
          where (q(:,:,:,nq) < q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_leq)
          where (q(:,:,:,nq) <= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_geq)
          where (q(:,:,:,nq) >= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_equal)
          where (q(:,:,:,nq) == q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_neq)
          where (q(:,:,:,nq) /= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_and)
          where (.not.(q(:,:,:,nq)==0d0).and..not.(q(:,:,:,nq+1)==0d0))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_or)
          where (.not.(q(:,:,:,nq)==0d0).or..not.(q(:,:,:,nq+1)==0d0))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case default
          goto 999
       end select
    elseif (istype(c,'unary')) then
       ! a unary operator or function
       if (nq < 1) goto 999
       select case(c)
       case (fun_uplus)
          q(:,:,:,nq) = +q(:,:,:,nq)
       case (fun_uminus)
          q(:,:,:,nq) = -q(:,:,:,nq)
       case (fun_abs)
          q(:,:,:,nq) = abs(q(:,:,:,nq))
       case (fun_exp)
          q(:,:,:,nq) = exp(q(:,:,:,nq))
       case (fun_sqrt)
          q(:,:,:,nq) = sqrt(q(:,:,:,nq))
       case (fun_floor)
          q(:,:,:,nq) = floor(q(:,:,:,nq))
       case (fun_ceiling)
          q(:,:,:,nq) = ceiling(q(:,:,:,nq))
       case (fun_round)
          q(:,:,:,nq) = nint(q(:,:,:,nq))
       case (fun_log)
          q(:,:,:,nq) = log(q(:,:,:,nq))
       case (fun_log10)
          q(:,:,:,nq) = log10(q(:,:,:,nq))
       case (fun_sin)
          q(:,:,:,nq) = sin(q(:,:,:,nq))
       case (fun_asin)
          q(:,:,:,nq) = asin(q(:,:,:,nq))
       case (fun_cos)
          q(:,:,:,nq) = cos(q(:,:,:,nq))
       case (fun_acos)
          q(:,:,:,nq) = acos(q(:,:,:,nq))
       case (fun_tan)
          q(:,:,:,nq) = tan(q(:,:,:,nq))
       case (fun_atan)
          q(:,:,:,nq) = atan(q(:,:,:,nq))
       case (fun_sinh)
          q(:,:,:,nq) = sinh(q(:,:,:,nq))
       case (fun_cosh)
          q(:,:,:,nq) = cosh(q(:,:,:,nq))
       case (fun_erf)
          q(:,:,:,nq) = erf(q(:,:,:,nq))
       case (fun_erfc)
          q(:,:,:,nq) = erfc(q(:,:,:,nq))
       case default
          goto 999
       end select
    else
       goto 999
    end if

    return
999 continue
    fail = .true.

  end subroutine pop_grid

  !> Calculate a chemical function for a given field.  This routine is
  !> thread-safe.
  function chemfunction(c,sia,x0,args,syl,periodic,errmsg) result(q)
    use systemmod, only: system
    use global, only: dunit0, iunit
    use tools_math, only: bhole
    use tools_io, only: string
    use types, only: scalar_value, field_evaluation_avail, fieldeval_category_f, fieldeval_category_fder1,&
       fieldeval_category_fder2, fieldeval_category_gkin, fieldeval_category_stress,&
       fieldeval_category_spin, fieldeval_category_mep,&
       fieldeval_category_uslater, fieldeval_category_xhole, fieldeval_category_nheff
    use param, only: pi
    integer, intent(in) :: c
    character*(*), intent(in) :: sia
    real*8, intent(in) :: x0(3)
    real*8, intent(in), allocatable :: args(:)
    real*8 :: q
    logical, intent(in), optional :: periodic
    type(system), intent(inout), optional :: syl
    character(len=:), allocatable, intent(inout) :: errmsg

    type(scalar_value) :: res
    real*8 :: f0, ds, ds0, g, g0, dsigs, quads, tau, drhos2, rhos, laps
    real*8 :: br_b, br_alf, br_a, raux(3)
    integer :: idx
    logical :: dohole, use1, use2
    type(field_evaluation_avail) :: request

    ! a constant
    real*8, parameter :: ctf = 2.8712340001881911d0 ! Thomas-Fermi k.e.d. constant, 3/10 * (3*pi^2)^(2/3)

    q = 0d0
    if (.not.present(syl)) return
    if (.not.syl%isinit) return

    idx = syl%fieldname_to_idx(sia)
    if (.not.syl%goodfield(id=idx)) then
       errmsg = "Invalid field or incorrect number of arguments in expression."
       return
    end if

    call request%clear()
    request%avail(fieldeval_category_f) = .true.
    if (c == fun_rdg) then
       request%avail(fieldeval_category_fder1) = .true.
    elseif (c == fun_vtf .or. c == fun_htf .or. c == fun_gtf_kir .or. c == fun_vtf_kir .or.&
       c == fun_htf_kir .or. c == fun_l .or. c == fun_lol_kir) then
       request%avail(fieldeval_category_fder1) = .true.
       request%avail(fieldeval_category_fder2) = .true.
    elseif (c == fun_gkin .or. c == fun_lol) then
       request%avail(fieldeval_category_gkin) = .true.
    elseif (c == fun_elf) then
       request%avail(fieldeval_category_fder1) = .true.
       request%avail(fieldeval_category_gkin) = .true.
    elseif (c == fun_kkin .or. c == fun_brhole_a1 .or. c == fun_brhole_a2 .or. c == fun_brhole_a .or.&
       c == fun_brhole_b1 .or. c == fun_brhole_b2 .or. c == fun_brhole_b .or. c == fun_brhole_alf1 .or.&
       c == fun_brhole_alf2 .or. c == fun_brhole_alf .or. c == fun_xhcurv1 .or. c == fun_xhcurv2 .or.&
       c == fun_xhcurv .or. c == fun_dsigs1 .or. c == fun_dsigs2 .or. c == fun_dsigs) then
       request%avail(fieldeval_category_fder1) = .true.
       request%avail(fieldeval_category_fder2) = .true.
       request%avail(fieldeval_category_gkin) = .true.
    elseif (c == fun_vir) then
       request%avail(fieldeval_category_stress) = .true.
    elseif (c == fun_he) then
       request%avail(fieldeval_category_gkin) = .true.
       request%avail(fieldeval_category_stress) = .true.
    elseif (c == fun_mep) then
       request%avail(fieldeval_category_mep) = .true.
    elseif (c == fun_uslater) then
       request%avail(fieldeval_category_uslater) = .true.
    elseif (c == fun_nheff) then
       request%avail(fieldeval_category_nheff) = .true.
    elseif (c == fun_xhole) then
       request%avail(fieldeval_category_xhole) = .true.
       if (.not.allocated(args)) then
          errmsg = "xhole requires arguments for the reference point."
          return
       end if
       if (size(args,1) /= 3) then
          errmsg = "xhole requires three arguments for the reference point."
          return
       end if
       if (syl%c%ismolecule) then
          request%xref = args / dunit0(iunit) - syl%c%molx0
       else
          request%xref = syl%c%x2c(args)
       end if
    end if

    ! evaluate the field
    call syl%f(idx)%grd(x0,request,res,periodic)
    if (.not.res%valid) then
       errmsg = "Arithmetic expression evaluated outside the field domain"
       return
    end if
    if (.not.res%satisfied) then
       errmsg = "The requested arithmetic expression cannot be evaluated with this field"
       return
    end if

    select case(c)
    case (fun_gtf)
       ! Thomas-Fermi kinetic energy density for the uniform electron gas
       ! See Yang and Parr, Density-Functional Theory of Atoms and Molecules
       q = ctf * res%f**(5d0/3d0)
    case (fun_vtf)
       ! Potential energy density calculated using fun_gtf and the local
       ! virial theorem (2g(r) + v(r) = 1/4*lap(r)).
       q = ctf * res%f**(5d0/3d0)
       q = 0.25d0 * res%del2f - 2 * q
    case (fun_htf)
       ! Total energy density calculated using fun_gtf and the local
       ! virial theorem (h(r) + v(r) = 1/4*lap(r)).
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
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
    case (fun_vtf_kir)
       ! Potential energy density calculated using fun_gtf_kir and the
       ! local virial theorem (2g(r) + v(r) = 1/4*lap(r)).
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = 0.25d0 * res%del2f - 2 * q
    case (fun_htf_kir)
       ! Total energy density calculated using fun_gtf_kir and the
       ! local virial theorem (h(r) + v(r) = 1/4*lap(r)).
       f0 = max(res%f,1d-30)
       q = ctf * f0**(5d0/3d0) + &
          1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = 0.25d0 * res%del2f - q
    case (fun_gkin)
       ! G-kinetic energy density (sum grho * grho)
       q = res%gkin
    case (fun_kkin)
       ! K-kinetic energy density (sum rho * laprho)
       q = res%gkin - 0.25d0 * res%del2f
    case (fun_l)
       ! Lagrangian density (-1/4 * lap)
       q = - 0.25d0 * res%del2f
    case (fun_elf)
       ! Electron localization function
       ! Becke and Edgecombe J. Chem. Phys. (1990) 92, 5397-5403
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
       q = res%vir
    case (fun_he)
       ! Energy density, fun_vir + fun_gkin
       !   Keith et al. Int. J. Quantum Chem. (1996) 57, 183-198.
       q = res%vir + res%gkin
    case (fun_lol)
       ! Localized-orbital locator
       !   Schmider and Becke, J. Mol. Struct. (Theochem) (2000) 527, 51-61
       !   Schmider and Becke, J. Chem. Phys. (2002) 116, 3184-3193.
       q = ctf * res%f**(5d0/3d0) / max(res%gkin,1d-30)
       q = q / (1d0+q)
    case (fun_lol_kir)
       ! Localized-orbital locator using Kirzhnits k.e.d.
       !   Tsirelson and Stash, Acta Cryst. (2002) B58, 780.
       f0 = max(res%f,1d-30)
       g0 = ctf * f0**(5d0/3d0)
       g = g0 + 1/72d0 * res%gfmod**2 / f0 + 1d0/6d0 * res%del2f
       q = g0 / g
       q = q / (1d0+q)
    case (fun_rdg)
       ! Reduced density gradient
       ! s = |gradrho| / (2 * (3*pi^2)^(1/3) * rho^(4/3))
       if (res%f < 1d-30) then
          q = 0d0
       else
          f0 = max(res%f,1d-30)
          q = res%gfmod / (2d0 * (3d0*pi*pi)**(1d0/3d0) * f0**(4d0/3d0))
       end if
    case (fun_brhole_a1,fun_brhole_a2,fun_brhole_a,fun_brhole_b1,fun_brhole_b2,fun_brhole_b,&
          fun_brhole_alf1,fun_brhole_alf2,fun_brhole_alf,fun_xhcurv1,fun_xhcurv2,fun_xhcurv,&
          fun_dsigs1,fun_dsigs2,fun_dsigs)
       ! - brhole: Becke-Roussel (BR) hole parameters. The spherically
       ! averaged exchange hole in the BR model is an exponential
       ! A*exp(-alpha * r) at a distance b from the reference point.
       !  A.D. Becke and M.R. Roussel, Phys. Rev. A 39 (1989) 3761
       ! - xhcurv: curvature of the spherically averaged exchange
       ! hole at the reference point. (Q_sigma)
       ! - dsigs: leading coefficient of the same-spin pair density
       ! (D_sigma).
       dohole = (c==fun_brhole_a1).or.(c==fun_brhole_a2).or.(c==fun_brhole_a).or.&
          (c==fun_brhole_b1).or.(c==fun_brhole_b2).or.(c==fun_brhole_b).or.&
          (c==fun_brhole_alf1).or.(c==fun_brhole_alf2).or.(c==fun_brhole_alf)
       use1 = (c==fun_brhole_a1).or.(c==fun_brhole_b1).or.(c==fun_brhole_alf1).or.(c==fun_xhcurv1).or.(c==fun_dsigs1)
       use2 = (c==fun_brhole_a2).or.(c==fun_brhole_b2).or.(c==fun_brhole_alf2).or.(c==fun_xhcurv2).or.(c==fun_dsigs2)

       if (res%ev%avail(fieldeval_category_spin)) then
          if (use1) then
             call assign_bhole_variables(res%fspinval(1),res%lapspin(1),res%gkinspin(1),res%gfmodspinval(1),.false.)
             if (dohole) call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
          elseif (use2) then
             call assign_bhole_variables(res%fspinval(2),res%lapspinval(2),res%gkinspin(2),res%gfmodspinval(2),.false.)
             if (dohole) call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
          else
             call assign_bhole_variables(res%fspinval(1),res%lapspinval(1),res%gkinspin(1),res%gfmodspinval(1),.false.)
             if (dohole) then
                call bhole(rhos,quads,1d0,raux(1),raux(2),raux(3))
             else
                raux(1) = dsigs
                raux(2) = quads
             end if
             call assign_bhole_variables(res%fspinval(2),res%lapspinval(2),res%gkinspin(2),res%gfmodspinval(2),.false.)
             if (dohole) then
                call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
                br_b = 0.5d0 * (raux(1) + br_b)
                br_alf = 0.5d0 * (raux(2) + br_alf)
                br_a = 0.5d0 * (raux(3) + br_a)
             else
                dsigs = 0.5d0 * (dsigs + raux(1))
                quads = 0.5d0 * (quads + raux(2))
             end if
          end if
       else
          call assign_bhole_variables(res%f,res%del2f,res%gkin,res%gfmod,.true.)
          if (dohole) call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
       end if
       if (c == fun_brhole_a1 .or. c == fun_brhole_a2 .or. c == fun_brhole_a) then
          q = br_a
       elseif (c == fun_brhole_b1 .or. c == fun_brhole_b2 .or. c == fun_brhole_b) then
          q = br_b
       elseif (c == fun_brhole_alf1 .or. c == fun_brhole_alf2 .or. c == fun_brhole_alf) then
          q = br_alf
       elseif (c == fun_xhcurv1 .or. c == fun_xhcurv2 .or. c == fun_xhcurv) then
          q = quads
       elseif (c == fun_dsigs1 .or. c == fun_dsigs2 .or. c == fun_dsigs) then
          q = dsigs
       end if
    case (fun_mep,fun_uslater,fun_nheff,fun_xhole)
       q = res%fspc
    case default
       errmsg = "function not implemented: " // arith_usage(c)
    end select

  contains
    subroutine assign_bhole_variables(rhos_,laps_,tau_,gfmod_,dohalf)
      real*8, intent(in) :: rhos_, laps_, tau_, gfmod_
      logical, intent(in) :: dohalf

      if (dohalf) then
         rhos = 0.5d0 * rhos_
         laps = 0.5d0 * laps_
         tau = tau_
         drhos2 = (0.5d0 * gfmod_)
         drhos2 = drhos2 * drhos2
      else
         rhos = rhos_
         laps = laps_
         tau = 0.5d0 * tau_
         drhos2 = gfmod_ * gfmod_
      end if
      dsigs = tau - 0.25d0 * drhos2 / max(rhos,1d-30)
      quads = (laps - 2d0 * dsigs) / 6d0

    end subroutine assign_bhole_variables
  end function chemfunction

end submodule proc
