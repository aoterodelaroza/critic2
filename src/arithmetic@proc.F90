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
  ! function tokenize(expr,ntok,toklist,lpexit,syl)
  ! function iprec(c)
  ! function iassoc(c)
  ! function istype(c,type,iwantarg)
  ! function isspecialfield(fid)
  ! function isstructvar(fid,c,fder)
  ! recursive function fieldeval(fid,fder,errmsg,x0,sptr,periodic)
  ! function structvareval(svar,x0,syl,periodic) result(q)
  ! function isnumber(rval,expr,lp)
  ! function isoperator(c,expr,lp)
  ! function isfunction(c,expr,lp,wasop)
  ! function isidentifier(id,expr,lp,fder)
  ! subroutine pop(q,nq,s,ns,x0,sptr,periodic,errmsg)
  ! subroutine pop_grid(q,nq,s,ns,fail)
  ! function chemfunction(c,sia,x0,args,sptr,periodic,errmsg) result(q)
  ! function specialfieldeval(fid,syl,x0) result(res)

  ! enum for operators and functions
  integer, parameter :: fun_unknown     = 0  !< unknown operator
  integer, parameter :: fun_openpar     = 1  !< open parenthesis
  integer, parameter :: fun_closepar    = 2  !< close parenthesis
  integer, parameter :: fun_uplus       = 3  !< unary +
  integer, parameter :: fun_uminus      = 4  !< unary -
  integer, parameter :: fun_abs         = 5  !< abs(.)
  integer, parameter :: fun_exp         = 6  !< exp(.)
  integer, parameter :: fun_sqrt        = 7  !< sqrt(.)
  integer, parameter :: fun_floor       = 8  !< floor(.)
  integer, parameter :: fun_ceiling     = 9  !< ceiling(.)
  integer, parameter :: fun_round       = 10 !< round(.)
  integer, parameter :: fun_log         = 11 !< log(.)
  integer, parameter :: fun_log10       = 12 !< log10(.)
  integer, parameter :: fun_sin         = 13 !< sin(.)
  integer, parameter :: fun_asin        = 14 !< asin(.)
  integer, parameter :: fun_cos         = 15 !< cos(.)
  integer, parameter :: fun_acos        = 16 !< acos(.)
  integer, parameter :: fun_tan         = 17 !< tan(.)
  integer, parameter :: fun_atan        = 18 !< atan(.)
  integer, parameter :: fun_atan2       = 19 !< atan2(.)
  integer, parameter :: fun_sinh        = 20 !< sinh(.)
  integer, parameter :: fun_cosh        = 21 !< cosh(.)
  integer, parameter :: fun_erf         = 22 !< erf(.)
  integer, parameter :: fun_erfc        = 23 !< erfc(.)
  integer, parameter :: fun_min         = 24 !< min(.)
  integer, parameter :: fun_max         = 25 !< max(.)
  integer, parameter :: fun_power       = 26 !< **
  integer, parameter :: fun_leq         = 27 !< <=
  integer, parameter :: fun_geq         = 28 !< >=
  integer, parameter :: fun_equal       = 29 !< ==
  integer, parameter :: fun_neq         = 30 !< !=
  integer, parameter :: fun_and         = 31 !< &&
  integer, parameter :: fun_or          = 32 !< ||
  integer, parameter :: fun_plus        = 33 !< binary +
  integer, parameter :: fun_minus       = 34 !< binary -
  integer, parameter :: fun_prod        = 35 !< *
  integer, parameter :: fun_div         = 36 !< /
  integer, parameter :: fun_modulo      = 37 !< %
  integer, parameter :: fun_great       = 38 !< >
  integer, parameter :: fun_less        = 39 !< <
  integer, parameter :: fun_xc          = 40 !< xc(.,...)  (X)
  integer, parameter :: fun_gtf         = 41 !< Thomas-Fermi kinetic energy density
  integer, parameter :: fun_vtf         = 42 !< Potential energy density calcd using gtf and the local virial theorem
  integer, parameter :: fun_htf         = 43 !< Total energy density calcd using gtf and the local virial theorem
  integer, parameter :: fun_gtf_kir     = 44 !< Thomas-Fermi ked with Kirzhnits gradient correction
  integer, parameter :: fun_vtf_kir     = 45 !< Potential energy density calcd using gtf_kir and the local virial theorem
  integer, parameter :: fun_htf_kir     = 46 !< Total energy density calcd using gtf_kir and the local virial theorem
  integer, parameter :: fun_gkin        = 47 !< Kinetic energy density, G-version (grho * grho)
  integer, parameter :: fun_kkin        = 48 !< Kinetic energy density, K-version (rho * laprho)
  integer, parameter :: fun_l           = 49 !< Lagrangian density (-1/4 * laprho)
  integer, parameter :: fun_elf         = 50 !< Electron localization function (ELF)
  integer, parameter :: fun_vir         = 51 !< Electronic potential energy density, virial field
  integer, parameter :: fun_he          = 52 !< Energy density, G + V
  integer, parameter :: fun_lol         = 53 !< Localized-orbital locator (LOL)
  integer, parameter :: fun_lol_kir     = 54 !< Localized-orbital locator (LOL), with Kirzhnits G
  integer, parameter :: fun_brhole_a1   = 55 !< BR hole, A prefactor, spin up
  integer, parameter :: fun_brhole_a2   = 56 !< BR hole, A prefactor, spin down
  integer, parameter :: fun_brhole_a    = 57 !< BR hole, A prefactor, spin average
  integer, parameter :: fun_brhole_alf1 = 58 !< BR hole, alpha exponent, spin up
  integer, parameter :: fun_brhole_alf2 = 59 !< BR hole, alpha exponent, spin down
  integer, parameter :: fun_brhole_alf  = 60 !< BR hole, alpha exponent, spin average
  integer, parameter :: fun_brhole_b1   = 61 !< BR hole, hole distance (b), spin up
  integer, parameter :: fun_brhole_b2   = 62 !< BR hole, hole distance (b), spin down
  integer, parameter :: fun_brhole_b    = 63 !< BR hole, hole distance (b), spin average
  integer, parameter :: fun_xhcurv1     = 64 !< exchange-hole curvature, spin up
  integer, parameter :: fun_xhcurv2     = 65 !< exchange-hole curvature, spin down
  integer, parameter :: fun_xhcurv      = 66 !< exchange-hole curvature, spin average
  integer, parameter :: fun_dsigs1      = 67 !< leading coefficient same-spin pair density, spin up
  integer, parameter :: fun_dsigs2      = 68 !< leading coefficient same-spin pair density, spin down
  integer, parameter :: fun_dsigs       = 69 !< leading coefficient same-spin pair density, spin avg.
  integer, parameter :: fun_mep         = 70 !< Molecular electrostatic potential
  integer, parameter :: fun_uslater     = 71 !< Slater potential
  integer, parameter :: fun_nheff       = 72 !< Effective exchange hole normalization
  integer, parameter :: fun_xhole       = 73 !< exchange-hole
  integer, parameter :: fun_rdg         = 74 !< reduced density gradient (RDG)

  ! enum for structural variables
  integer, parameter :: svar_dnuc    = 1  !< Distance to the closest nucleus
  integer, parameter :: svar_xnucx   = 2  !< x of the closest nucleus (crystallographic)
  integer, parameter :: svar_ynucx   = 3  !< y of the closest nucleus (crystallographic)
  integer, parameter :: svar_znucx   = 4  !< z of the closest nucleus (crystallographic)
  integer, parameter :: svar_xnucc   = 5  !< x of the closest nucleus (Cartesian)
  integer, parameter :: svar_ynucc   = 6  !< y of the closest nucleus (Cartesian)
  integer, parameter :: svar_znucc   = 7  !< z of the closest nucleus (Cartesian)
  integer, parameter :: svar_x       = 8  !< x of the evaluation point (default units)
  integer, parameter :: svar_y       = 9  !< y of the evaluation point (default units)
  integer, parameter :: svar_z       = 10 !< z of the evaluation point (default units)
  integer, parameter :: svar_xx      = 11 !< x of the evaluation point (crystallographic)
  integer, parameter :: svar_yx      = 12 !< y of the evaluation point (crystallographic)
  integer, parameter :: svar_zx      = 13 !< z of the evaluation point (crystallographic)
  integer, parameter :: svar_xc      = 14 !< x of the evaluation point (Cartesian)
  integer, parameter :: svar_yc      = 15 !< y of the evaluation point (Cartesian)
  integer, parameter :: svar_zc      = 16 !< z of the evaluation point (Cartesian)
  integer, parameter :: svar_xm      = 17 !< x of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_ym      = 18 !< y of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_zm      = 19 !< z of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_xxr     = 20 !< x of the evaluation point (reduced cryst.)
  integer, parameter :: svar_yxr     = 21 !< y of the evaluation point (reduced cryst.)
  integer, parameter :: svar_zxr     = 22 !< z of the evaluation point (reduced cryst.)
  integer, parameter :: svar_idnuc   = 23 !< complete-list id of the closest nucleus
  integer, parameter :: svar_nidnuc  = 24 !< non-equivalent-list id of the closest nucleus
  integer, parameter :: svar_rho0nuc = 25 !< atomic density contribution from the closest nucleus
  integer, parameter :: svar_spcnuc  = 26 !< species id of the closest nucleus
  integer, parameter :: svar_zatnuc  = 27 !< atomic number of the closest nucleus

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
    use systemmod, only: system
    use iso_c_binding, only: c_ptr, c_f_pointer
    use hashmod, only: hash
    real*8 :: eval
    character(*), intent(in) :: expr
    character(len=:), allocatable, intent(inout) :: errmsg
    real*8, intent(in), optional :: x0(3)
    type(c_ptr), intent(in), optional :: sptr
    logical, intent(in), optional :: periodic
    type(token), intent(in), optional :: toklistin(:)

    integer :: i, ntok, lp
    integer :: c, s(100)
    logical :: again, ok
    real*8 :: q(100)
    integer :: nq, ns
    type(token), allocatable :: toklist(:)
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

    ! tokenize the expression or use the tokens in input
    if (present(toklistin)) then
       toklist = toklistin
       ntok = size(toklistin,1)
    else
       lp = 1
       ok = tokenize(expr,ntok,toklist,lp,syl)
       if (.not.ok) then
          errmsg = 'syntax error evaluating: ' // trim(expr)
          return
       end if
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
                   call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
                   if (len_trim(errmsg) > 0) return
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
              call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
              if (len_trim(errmsg) > 0) return
           end do
           if (ns == 0) then
              errmsg = 'mismatched parentheses'
              return
           end if
           ns = ns - 1
           ! if the top of the stack is a function, pop it
           if (ns > 0) then
              c = s(ns)
              if (istype(c,'function')) then
                 call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
                 if (len_trim(errmsg) > 0) return
              end if
           end if
        elseif (toklist(i)%type == token_comma) then
           ! a comma
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop(q,nq,s,ns,x0,syl,periodic,errmsg)
              if (len_trim(errmsg) > 0) return
           end do
           if (s(ns) /= fun_openpar) then
              errmsg = 'mismatched parentheses'
              return
           end if
        elseif (toklist(i)%type == token_field) then
           ! a field
           nq = nq + 1
           q(nq) = fieldeval(toklist(i)%sval,toklist(i)%fder,errmsg,x0,syl,periodic)
           if (len_trim(errmsg) > 0) return
        elseif (toklist(i)%type == token_structvar) then
           ! a structural variable
           nq = nq + 1
           q(nq) = structvareval(toklist(i)%ival,toklist(i)%fder,errmsg,x0,syl,periodic)
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
    eval = q(1)

  end function eval

  !> Evaluate an arithmetic expression (expr) on a grid with
  !> dimensions n. sptr = C pointer to the system. Returns the grid
  !> of values in f (and iok=.true.) if success, or zero and iok=.false.
  !> if failed.
  module subroutine eval_grid(n,expr,sptr,f,iok)
    use systemmod, only: system
    use types, only: realloc
    use iso_c_binding, only: c_ptr, c_f_pointer
    integer, intent(in) :: n(3)
    character(*), intent(in) :: expr
    type(c_ptr), intent(in) :: sptr
    real*8, intent(out) :: f(n(1),n(2),n(3))
    logical, intent(out) :: iok

    integer :: i, ntok, lp, idx
    integer :: c, s(100)
    logical :: again, ok, ifail
    integer :: nq, ns
    real*8, allocatable :: q(:,:,:,:)
    type(token), allocatable :: toklist(:)
    type(system), pointer :: syl

    ! recover the system pointer
    syl => null()
    call c_f_pointer(sptr,syl)
    if (.not.syl%isinit) goto 999

    ! tokenize the expression in input
    iok = .true.
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,syl)
    if (.not.ok) then
       goto 999
       return
    end if

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
          ! left parenthesis
          ns = ns + 1
          s(ns) = fun_openpar
       elseif (toklist(i)%type == token_rpar) then
          ! right parenthesis
          do while (ns > 0)
             if (s(ns) == fun_openpar) exit
             call pop_grid(q,nq,s,ns,ifail)
             if (ifail) goto 999
          end do
          if (ns == 0) goto 999
          ns = ns - 1
          ! if the top of the stack is a function, pop it
          if (ns > 0) then
             c = s(ns)
             if (istype(c,'function')) then
                call pop_grid(q,nq,s,ns,ifail)
                if (ifail) goto 999
             end if
          end if
       elseif (toklist(i)%type == token_comma) then
          ! a comma
          do while (ns > 0)
             if (s(ns) == fun_openpar) exit
             call pop_grid(q,nq,s,ns,ifail)
             if (ifail) goto 999
          end do
          if (s(ns) /= fun_openpar) goto 999
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
    f = q(:,:,:,1)
    if (allocated(q)) deallocate(q)

    return
999 continue
    iok = .false.
    if (allocated(q)) deallocate(q)
    f = 0d0
    return

  end subroutine eval_grid

  !> Parse an expression (expr) and return the n field ids in it
  !> (idlist).  sptr = C pointer to the calling system. This routine
  !> is thread-safe. If error, return non-zero-length string in
  !> errmsg.
  module subroutine fields_in_eval(expr,errmsg,n,idlist,sptr)
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

    integer :: lp, i
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
    ok = tokenize(expr,ntok,toklist,lp,syl)
    if (.not. ok) then
       errmsg = "syntax error"
       return
    end if

    ! return the ids of the fields in an array
    n = 0
    do i = 1, ntok
       if (toklist(i)%type == token_field .and..not.isspecialfield(toklist(i)%sval)) then
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
    ok = tokenize(expr,ntok,toklist,lp,syl)
    if (.not.ok) then
       errmsg = 'syntax error in expression: ' // trim(expr)
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
  function tokenize(expr,ntok,toklist,lpexit,syl)
    use systemmod, only: system
    use tools_io, only: lower, isinteger
    use param, only: vh, tab
    logical :: tokenize
    character(*), intent(in) :: expr
    integer, intent(out) :: ntok
    type(token), intent(inout), allocatable :: toklist(:)
    integer, intent(inout) :: lpexit
    type(system), intent(inout), optional :: syl

    integer :: lp, ll
    character(len=:), allocatable :: str
    character*10 :: fder
    logical :: ok, wasop, inchem
    real*8 :: a
    integer :: c, npar, id, lp2

    ! recover the system pointer
    if (present(syl)) then
       if (.not.syl%isinit) goto 999
    end if

    ! initialize
    tokenize = .true.
    inchem = .false.
    npar = 0

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
       if (isnumber(a,expr,lp)) then
          ! a number (without sign)
          call addtok(token_num,fval=a)
          wasop = .false.
       elseif (isfunction(c,expr,lp,wasop)) then
          ! a function
          call addtok(token_fun,ival=c)
          wasop = .false.
          if (istype(c,'chemfunction')) then
             inchem = .true.
             npar = 0
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
          call addtok(token_rpar)
          wasop = .false.
          lp = lp + 1
          npar = npar - 1
          if (npar <= 0 .and. inchem) inchem = .false.
       elseif (expr(lp:lp) == ',') then
          ! a comma
          call addtok(token_comma)
          wasop = .true.
          lp = lp + 1
       elseif (expr(lp:lp) == '@') then
          ! a structural variable
          lp = lp + 1
          ok = isidentifier(str,expr,lp,fder)
          if (.not.ok) goto 999
          ok = isstructvar(str,c,fder)
          if (.not.ok) goto 999
          call addtok(token_structvar,ival=c,fder=fder)
          wasop=.false.

       elseif (expr(lp:lp) == '$') then
          ! a field read the field identifier and the : identifier
          lp = lp + 1
          ok = isidentifier(str,expr,lp,fder)
          if (.not.ok) goto 999
          ok = .false.

          ! an integer or a known key
          lp2 = 1
          ok = isinteger(id,str,lp2)
          if (present(syl)) ok = ok .or. syl%fh%iskey(trim(str))

          if (.not.ok) then
             ! a special field
             str = lower(str)
             ok = isspecialfield(trim(str))
             if (.not.ok) goto 999
             fder = ""
          end if
          if (.not.inchem) then
             ! normal interpretation
             call addtok(token_field,sval=str,fder=fder)
          else
             ! inside a chemical function
             if (.not.present(syl)) goto 999
             ok = syl%fh%iskey(trim(str))
             if (.not.ok) goto 999
             id = syl%fh%get(trim(str),1)
             call addtok(token_num,fval=real(id,8))
          end if
          wasop = .false.
       elseif (isidentifier(str,expr,lp)) then
          if (.not.inchem) then
             ! a constant or a variable
             ok = vh%iskey(trim(str))
             if (.not.ok) goto 999
             call addtok(token_num,fval=vh%get(trim(str),a))
          else
             ! inside a chemical function -> field identifier
             if (.not.present(syl)) goto 999
             ok = syl%fh%iskey(trim(str))
             if (.not.ok) goto 999
             id = syl%fh%get(trim(str),1)
             call addtok(token_num,fval=real(id,8))
          end if
          wasop = .false.
       else
          goto 999
       end if
    end do main

    lpexit = lp
    tokenize = .true.
    return

999 continue ! could not tokenize
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
  end function tokenize

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
  !> otherwise, return zero. This routine is thread-safe.
  function istype(c,type,iwantarg)
    integer, intent(in) :: c
    character*(*), intent(in) :: type
    integer, intent(inout), optional :: iwantarg
    logical :: istype

    if (present(iwantarg)) iwantarg = 0
    istype = .false.
    if (type == 'function') then
       istype = &
          c == fun_abs .or. c == fun_exp .or. c == fun_sqrt .or. &
          c == fun_floor .or. c == fun_ceiling .or. c == fun_round .or. &
          c == fun_log .or. c == fun_log10 .or. c == fun_sin .or. &
          c == fun_asin .or. c == fun_cos .or. c == fun_acos .or. &
          c == fun_tan .or. c == fun_atan .or. c == fun_atan2 .or. &
          c == fun_sinh .or. c == fun_cosh .or. c == fun_erf .or. &
          c == fun_erfc .or. c == fun_min .or. c == fun_max .or. &
          c == fun_xc .or. c == fun_gtf .or. c == fun_vtf .or. c == fun_htf .or. &
          c == fun_gtf_kir .or. c == fun_vtf_kir .or. c == fun_htf_kir .or.&
          c == fun_gkin .or. c == fun_kkin .or. c == fun_l .or. c == fun_elf .or.&
          c == fun_vir .or. c == fun_he .or. c == fun_lol .or. c == fun_lol_kir .or.&
          c == fun_brhole_a1 .or. c == fun_brhole_a2 .or. c == fun_brhole_a .or. &
          c == fun_brhole_alf1 .or. c == fun_brhole_alf2 .or. c == fun_brhole_alf .or. &
          c == fun_brhole_b1 .or. c == fun_brhole_b2 .or. c == fun_brhole_b .or. &
          c == fun_xhcurv1 .or. c == fun_xhcurv2 .or. c == fun_xhcurv .or.&
          c == fun_dsigs1 .or. c == fun_dsigs2 .or. c == fun_dsigs .or.&
          c == fun_mep .or. c == fun_uslater .or. c == fun_nheff .or. c == fun_xhole .or.&
          c == fun_rdg
    elseif (type == 'chemfunction') then
       istype = &
          c == fun_gtf .or. c == fun_vtf .or. c == fun_htf .or. &
          c == fun_gtf_kir .or. c == fun_vtf_kir .or. c == fun_htf_kir .or.&
          c == fun_gkin .or. c == fun_kkin .or. c == fun_l .or. c == fun_elf .or.&
          c == fun_vir .or. c == fun_he .or. c == fun_lol .or. c == fun_lol_kir .or.&
          c == fun_brhole_a1 .or. c == fun_brhole_a2 .or. c == fun_brhole_a .or. &
          c == fun_brhole_alf1 .or. c == fun_brhole_alf2 .or. c == fun_brhole_alf .or. &
          c == fun_brhole_b1 .or. c == fun_brhole_b2 .or. c == fun_brhole_b .or.&
          c == fun_xhcurv1 .or. c == fun_xhcurv2 .or. c == fun_xhcurv .or.&
          c == fun_dsigs1 .or. c == fun_dsigs2 .or. c == fun_dsigs .or.&
          c == fun_mep .or. c == fun_uslater .or. c == fun_nheff .or. c == fun_xhole .or.&
          c == fun_rdg
       if (present(iwantarg)) then
          if (c == fun_xhole) then
             iwantarg = 4
          else
             iwantarg = 1
          end if
       end if

    elseif (type == 'operator') then
       istype = &
          c == fun_power .or. c == fun_leq .or. c == fun_geq .or.&
          c == fun_equal .or. c == fun_neq .or. c == fun_and .or.&
          c == fun_or .or. c == fun_plus .or. c == fun_minus .or.&
          c == fun_prod .or. c == fun_div .or. c == fun_modulo .or.&
          c == fun_great .or. c == fun_less
    elseif (type == 'binary') then
       istype = &
          c == fun_atan2 .or. c == fun_min .or. c == fun_max .or.&
          c == fun_power .or. c == fun_leq .or. c == fun_geq .or.&
          c == fun_equal .or. c == fun_neq .or. c == fun_and .or.&
          c == fun_or .or. c == fun_plus .or. c == fun_minus .or.&
          c == fun_prod .or. c == fun_div .or. c == fun_modulo .or.&
          c == fun_great .or. c == fun_less
    elseif (type == 'unary') then
       istype = &
          c == fun_uplus .or. c == fun_uminus .or. c == fun_abs .or.&
          c == fun_exp .or. c == fun_sqrt .or. c == fun_floor .or.&
          c == fun_ceiling .or. c == fun_round .or. c == fun_log .or.&
          c == fun_log10 .or. c == fun_sin .or. c == fun_asin .or.&
          c == fun_cos .or. c == fun_acos .or. c == fun_tan .or.&
          c == fun_atan .or. c == fun_sinh .or. c == fun_cosh .or.&
          c == fun_erf .or. c == fun_erfc
    endif

  endfunction istype

  !> Does this identifier correspond to a special field. This routine
  !> is thread-safe.
  function isspecialfield(fid)
    use tools_io, only: lower
    character*(*), intent(in) :: fid
    logical :: isspecialfield

    isspecialfield = (trim(fid) == "ewald")

  end function isspecialfield

  !> If this identifier corresponds to a structural variable, return
  !> .true. and the integer ID for the variable in c. Otherwise,
  !> return .false.. This routine is thread-safe.
  function isstructvar(fid,c,fder)
    use tools_io, only: lower
    character*(*), intent(in) :: fid
    integer, intent(out) :: c
    character*(*), intent(in) :: fder
    logical :: isstructvar
    logical :: fderallow

    isstructvar = .true.
    select case (trim(lower(fid)))
    case("dnuc")
       c = svar_dnuc
    case("xnucx")
       c = svar_xnucx
    case("ynucx")
       c = svar_ynucx
    case("znucx")
       c = svar_znucx
    case("xnucc")
       c = svar_xnucc
    case("ynucc")
       c = svar_ynucc
    case("znucc")
       c = svar_znucc
    case("x")
       c = svar_x
    case("y")
       c = svar_y
    case("z")
       c = svar_z
    case("xx")
       c = svar_xx
    case("yx")
       c = svar_yx
    case("zx")
       c = svar_zx
    case("xc")
       c = svar_xc
    case("yc")
       c = svar_yc
    case("zc")
       c = svar_zc
    case("xm")
       c = svar_xm
    case("ym")
       c = svar_ym
    case("zm")
       c = svar_zm
    case("xxr")
       c = svar_xxr
    case("yxr")
       c = svar_yxr
    case("zxr")
       c = svar_zxr
    case("idnuc")
       c = svar_idnuc
    case("nidnuc")
       c = svar_nidnuc
    case("rho0nuc")
       c = svar_rho0nuc
    case("spcnuc")
       c = svar_spcnuc
    case("zatnuc")
       c = svar_zatnuc
    case default
       isstructvar = .false.
    end select

    fderallow = (c==svar_dnuc).or.(c==svar_xnucx).or.(c==svar_ynucx).or.(c==svar_znucx).or.&
       (c==svar_xnucc).or.(c==svar_ynucc).or.(c==svar_znucc).or.(c==svar_rho0nuc)

    if (.not.fderallow.and.len_trim(fder) > 0) isstructvar = .false.

  end function isstructvar

  !> Evaluate field with identifier fid and field flag fder at point
  !> x0 and return the value. syl = calling system. If periodic is
  !> true, evaluate the field under pbc. Returns non-empty errmsg if
  !> the evaluation failed.
  recursive function fieldeval(fid,fder,errmsg,x0,syl,periodic)
    use systemmod, only: system
    use tools_math, only: det3sym
    use tools_io, only: string, isinteger, lower
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

    integer :: nder, idx
    type(scalar_value) :: res
    character*10 :: fderl
    type(field_evaluation_avail) :: request
    logical :: ok

    errmsg = ""

    ! recover the system pointer
    if (present(syl)) then
       if (.not.syl%isinit) then
          errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - system not initialized'
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
                if (.not.ok) then
                   errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - wrong field modifier'
                   return
                end if
             elseif (fderl(1:1) == "b") then
                request%moini = id_mo_b
                ok = isinteger(request%moend,fderl(2:))
                if (.not.ok) then
                   errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - wrong field modifier'
                   return
                end if
             else
                request%moini = id_mo_id
                ok = isinteger(request%moend,fderl)
                if (.not.ok) then
                   errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - wrong field modifier'
                   return
                end if
             end if
          end select

          ! run the field evaluation
          call syl%f(idx)%grd(x0,request,res,periodic=periodic)
          if (.not.res%valid) then
             errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - invalid evaluation result'
             fieldeval = 0d0
             return
          end if
          if (.not.res%satisfied) then
             errmsg = 'field ' // string(fid) // ':' // string(fder) //&
                ' - cannot compute the property with this field type'
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
       else if (isspecialfield(fid)) then
          fieldeval = specialfieldeval(fid,syl,x0)
       else
          errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - unknown field'
          return
       end if
    else
       errmsg = 'field ' // string(fid) // ':' // string(fder) // ' - no point for field evaluation'
       return
    end if

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
    use tools_io, only: lower

    logical :: isfunction
    character*(*), intent(in) :: expr
    integer, intent(inout) :: lp
    integer, intent(out) :: c
    logical, intent(in) :: wasop

    integer :: lpo, ll
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
       select case (trim(lower(word)))
       case ("abs")
          c = fun_abs
       case ("exp")
          c = fun_exp
       case ("sqrt")
          c = fun_sqrt
       case ("floor")
          c = fun_floor
       case ("ceil")
          c = fun_ceiling
       case ("ceiling")
          c = fun_ceiling
       case ("round")
          c = fun_round
       case ("log")
          c = fun_log
       case ("log10")
          c = fun_log10
       case ("sin")
          c = fun_sin
       case ("asin")
          c = fun_asin
       case ("cos")
          c = fun_cos
       case ("acos")
          c = fun_acos
       case ("tan")
          c = fun_tan
       case ("atan")
          c = fun_atan
       case ("atan2")
          c = fun_atan2
       case ("sinh")
          c = fun_sinh
       case ("cosh")
          c = fun_cosh
       case ("erf")
          c = fun_erf
       case ("erfc")
          c = fun_erfc
       case ("min")
          c = fun_min
       case ("max")
          c = fun_max
       case ("xc")
          c = fun_xc
       case ("gtf")
          c = fun_gtf
       case ("vtf")
          c = fun_vtf
       case ("htf")
          c = fun_htf
       case ("gtf_kir")
          c = fun_gtf_kir
       case ("vtf_kir")
          c = fun_vtf_kir
       case ("htf_kir")
          c = fun_htf_kir
       case ("gkin")
          c = fun_gkin
       case ("kkin")
          c = fun_kkin
       case ("lag")
          c = fun_l
       case ("elf")
          c = fun_elf
       case ("vir")
          c = fun_vir
       case ("he")
          c = fun_he
       case ("lol")
          c = fun_lol
       case ("lol_kir")
          c = fun_lol_kir
       case ("rdg")
          c = fun_rdg
       case ("brhole_a1")
          c = fun_brhole_a1
       case ("brhole_a2")
          c = fun_brhole_a2
       case ("brhole_a")
          c = fun_brhole_a
       case ("brhole_alf1")
          c = fun_brhole_alf1
       case ("brhole_alf2")
          c = fun_brhole_alf2
       case ("brhole_alf")
          c = fun_brhole_alf
       case ("brhole_b1")
          c = fun_brhole_b1
       case ("brhole_b2")
          c = fun_brhole_b2
       case ("brhole_b")
          c = fun_brhole_b
       case ("xhcurv1")
          c = fun_xhcurv1
       case ("xhcurv2")
          c = fun_xhcurv2
       case ("xhcurv")
          c = fun_xhcurv
       case ("dsigs1")
          c = fun_dsigs1
       case ("dsigs2")
          c = fun_dsigs2
       case ("dsigs")
          c = fun_dsigs
       case ("mep")
          c = fun_mep
       case ("uslater")
          c = fun_uslater
       case ("xhole")
          c = fun_xhole
       case ("nheff")
          c = fun_nheff
       case default
          lp = lpo
          return
       end select
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
    if (expr(lp:lp) >= 'a' .and. expr(lp:lp)<='z' .or. &
       expr(lp:lp) >= 'A' .and. expr(lp:lp)<='Z' .or. &
       expr(lp:lp) >= '0' .and. expr(lp:lp)<='9' .or. &
       expr(lp:lp) == '_') then
       ! read the identifier
       word = ""
       do while (expr(lp:lp) >= 'a' .and. expr(lp:lp)<='z' .or. &
          expr(lp:lp) >= 'A' .and. expr(lp:lp)<='Z' .or. &
          expr(lp:lp) >= '0' .and. expr(lp:lp)<='9' .or. &
          expr(lp:lp) == '_')
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
             do while (isletter(expr(i:i)).or.isdigit(expr(i:i)))
                i = i + 1
                if (i > ll) exit
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
  subroutine pop(q,nq,s,ns,x0,syl,periodic,errmsg)
    use systemmod, only: system
    use tools_io, only: string
#ifdef HAVE_LIBXC
    use iso_c_binding, only: c_size_t
    use xc_f03_lib_m
#endif
    real*8, intent(inout) :: q(:)
    integer, intent(inout) :: s(:)
    integer, intent(inout) :: nq, ns
    real*8, intent(in), optional :: x0(3)
    type(system), intent(inout), optional :: syl
    logical, intent(in), optional :: periodic
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: ia, iwantarg
    integer :: c, i
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
       ia = nint(q(nq))
       if (abs(q(nq) - ia) > 1d-13) then
          errmsg = "need an integer for the functional ID"
          return
       end if

       ! dirty trick for syncing the threads
       if (.not.ifun(ia)%init) then
          !$omp critical (ifuninit)
          if (.not.ifun(ia)%init) then
             ifun(ia)%id = ia
             call xc_f03_func_init(ifun(ia)%conf, ifun(ia)%id, XC_UNPOLARIZED)
             ifun(ia)%info = xc_f03_func_get_info(ifun(ia)%conf)
             ifun(ia)%family = xc_f03_func_info_get_family(ifun(ia)%info)
             ifun(ia)%init = .true.
          end if
          !$omp end critical (ifuninit)
       endif

       if (ifun(ia)%family == XC_FAMILY_LDA .and. nq <= 1.or.&
          ifun(ia)%family == XC_FAMILY_GGA .and. nq <= 2.or.&
          ifun(ia)%family == XC_FAMILY_MGGA .and. nq <= 4) then
          errmsg = "insufficient argument list in xc"
          return
       endif

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
          errmsg = 'syntax error'
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
       end select
    elseif (istype(c,'unary')) then
       ! a unary operator or function
       if (nq < 1) then
          errmsg = 'syntax error'
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
       end select
    elseif (istype(c,'chemfunction',iwantarg)) then
       ! We need a point and the evaluator
       if (.not.present(x0).or..not.present(syl)) then
          errmsg = "a point is needed for evaluating chemical function"
          return
       endif

       ! Consume the extra arguments for this function
       if (iwantarg > 1) then
          if (nq < iwantarg) then
             errmsg = 'syntax error'
             return
          end if
          allocate(args(iwantarg-1))
          do i = 1, iwantarg-1
             args(i) = q(nq - (iwantarg-1) + i)
          end do
          nq = nq - (iwantarg - 1)
       end if

       ! An integer field identifier as the first argument
       ia = nint(q(nq))
       if (abs(ia - q(nq)) > 1d-13) then
          errmsg = "an integer field is required as argument"
          return
       end if
       write (sia,'(I8)') ia
       sia = adjustl(sia)
       if (.not.syl%goodfield(key=sia)) then
          errmsg = 'wrong field: ' // string(sia)
          return
       end if

       ! Use the library of chemical functions
       q(nq) = chemfunction(c,sia,x0,args,syl,periodic,errmsg)
       if (len_trim(errmsg) > 0) return
       if (iwantarg > 1) &
          deallocate(args)
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
    use fieldmod, only: type_wfn
    use global, only: dunit0, iunit
    use tools_math, only: bhole
    use tools_io, only: string
    use types, only: scalar_value, field_evaluation_avail, fieldeval_category_f, fieldeval_category_fder1,&
       fieldeval_category_fder2, fieldeval_category_gkin, fieldeval_category_stress,&
       fieldeval_category_spin, fieldeval_category_mo, fieldeval_category_mep,&
       fieldeval_category_uslater, fieldeval_category_xhole, fieldeval_category_nheff
    use param, only: pi
    integer, intent(in) :: c
    character*(*), intent(in) :: sia
    real*8, intent(in) :: x0(3)
    real*8, intent(in), allocatable :: args(:)
    real*8 :: q, xref(3)
    logical, intent(in), optional :: periodic
    type(system), intent(inout), optional :: syl
    character(len=:), allocatable, intent(inout) :: errmsg

    type(scalar_value) :: res
    real*8 :: f0, ds, ds0, g, g0, dsigs, quads, tau, drhos2, rhos, laps
    real*8 :: br_b, br_alf, br_a, raux(3), ux
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

  !> Evaluate a special field (id string fid) at point x0 (cryst. coords.)
  !> using system syl
  function specialfieldeval(fid,syl,x0) result(res)
    use systemmod, only: system
    character*(*), intent(in) :: fid
    type(system), intent(inout) :: syl
    real*8, intent(in) :: x0(3)
    real*8 :: res

    real*8 :: xp(3)
    real*8 :: rcut, hcut, eta, qsum
    integer :: lrmax(3), lhmax(3)

    if (trim(fid) == "ewald") then
       ! Ewald potential at point x0
       xp = syl%c%c2x(x0)
       call syl%c%calculate_ewald_cutoffs(rcut,hcut,eta,qsum,lrmax,lhmax)
       res = syl%c%ewald_pot(xp,rcut,hcut,eta,qsum,lrmax,lhmax)
    end if

  end function specialfieldeval

end submodule proc
