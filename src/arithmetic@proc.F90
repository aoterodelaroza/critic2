! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (arithmetic) proc
  implicit none

contains
  
  !> Return field ids from the evaluation of an expression.  This
  !> routine is thread-safe.
  module subroutine fields_in_eval(expr,fh,n,idlist)
    use hashmod, only: hash
    use tools_io, only: string
    use types, only: realloc
    character(*), intent(in) :: expr
    type(hash), intent(in) :: fh
    integer, intent(out) :: n
    character*255, allocatable, intent(inout) :: idlist(:)

    integer :: lp, i
    logical :: ok
    integer :: ntok
    type(token), allocatable :: toklist(:)

    ! allocate space for the field ids
    if (allocated(idlist)) deallocate(idlist)
    allocate(idlist(10))

    ! tokenize the expression
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,fh)
    if (.not. ok) &
       call die("error evaluating expression: " // string(expr))

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

  !> Given an expression in string expr starting at lpexit, parse all
  !> tokens for the arithmetic evaluation. Return the tokens in
  !> toklist and the number of tokens in ntok, and advance the string
  !> pointer lpexit.  This routine is thread-safe.
  module function tokenize(expr,ntok,toklist,lpexit,fh) 
    use hashmod, only: hash
    use tools_io, only: lower, isinteger
    use param, only: vh
    logical :: tokenize
    character(*), intent(in) :: expr
    integer, intent(out) :: ntok
    type(token), intent(inout), allocatable :: toklist(:)
    integer, intent(inout) :: lpexit
    type(hash), intent(in), optional :: fh

    integer :: lp, ll
    character(len=:), allocatable :: str
    character*10 :: fder
    logical :: ok, wasop, inchem
    real*8 :: a
    integer :: c, npar, id, lp2

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
       do while (expr(lp:lp) == ' '.or.expr(lp:lp) == '\t'.or.expr(lp:lp)=='\n'.or.&
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
       elseif (expr(lp:lp) == '$') then
          ! a field read the field identifier and the : identifier
          lp = lp + 1
          ok = isidentifier(str,expr,lp,fder)
          if (.not.ok) goto 999
          ok = .false.

          ! an integer or a known key
          lp2 = 1
          ok = isinteger(id,str,lp2)
          if (present(fh)) ok = ok .or. fh%iskey(trim(str))

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
             if (.not.present(fh)) goto 999
             ok = fh%iskey(trim(str))
             if (.not.ok) goto 999
             id = fh%get(trim(str),1)
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
             if (.not.present(fh)) goto 999
             ok = fh%iskey(trim(str))
             if (.not.ok) goto 999
             id = fh%get(trim(str),1)
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

  !> Read an unsigned number or return false and leave lp unchanged
  !> This routine is thread-safe.
  module function isnumber(rval,expr,lp)
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

999 continue
    isnumber=.true.
    read(expr(lp:i-1),*) rval
    lp = i

  end function isnumber

  !> Read a binary operator or return false and leave lp unchanged
  !> This routine is thread-safe.
  module function isoperator(c,expr,lp)
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
             call die("isoperator","unknown operator")
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
             call die("isoperator","unknown operator")
          end if
       endif
       lp = lp + 1
       isoperator = .true.
    end if

  end function isoperator

  !> Read a unary operator (function) or return false and leave lp
  !> unchanged This routine is thread-safe.
  module function isfunction(c,expr,lp,wasop)
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
  module function isidentifier(id,expr,lp,fder)
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

  !> Return an operator precedence. This routine is thread-safe.
  module function iprec(c)
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
  module function iassoc(c)
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
  !> binary).  This routine is thread-safe.
  module function istype(c,type)
    integer, intent(in) :: c
    character*(*), intent(in) :: type
    logical :: istype

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
          c == fun_vir .or. c == fun_he .or. c == fun_lol .or. c == fun_lol_kir
    elseif (type == 'chemfunction') then
       istype = &
          c == fun_gtf .or. c == fun_vtf .or. c == fun_htf .or. &
          c == fun_gtf_kir .or. c == fun_vtf_kir .or. c == fun_htf_kir .or.&
          c == fun_gkin .or. c == fun_kkin .or. c == fun_l .or. c == fun_elf .or.&
          c == fun_vir .or. c == fun_he .or. c == fun_lol .or. c == fun_lol_kir
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

  module subroutine die(msg,msg2)
    use tools_io, only: ferror, faterr
    character*(*), intent(in) :: msg
    character*(*), intent(in), optional :: msg2
    !$omp critical (error)
    call ferror('eval',msg,faterr,msg2)
    !$omp end critical (error)
  end subroutine die

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

  !> Does this identifier correspond to a special field. This routine
  !> is thread-safe.
  module function isspecialfield(fid)
    use tools_io, only: lower
    character*(*), intent(in) :: fid
    logical :: isspecialfield

    ! xxxx model1r and model1v are a hack
    isspecialfield = (trim(fid) == "ewald" .or. trim(fid) == "model1r" .or.&
       trim(fid) == "model1v") 
    
  end function isspecialfield

end submodule proc
