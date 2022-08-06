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

submodule (tools_io) proc
  implicit none

  ! string parameters
  integer, parameter :: maxlen = 255 !< maximum length of a number

  ! clock
  integer*8 :: stime0
  real :: ctime0

  ! logical unit allocation
  logical :: alloc(0:100) !< allocation flag array

contains

  !> Initialize the LU array.
  module subroutine lualloc_init()
    use iso_fortran_env, only: error_unit, input_unit, output_unit

    alloc = .false.
    alloc(error_unit) = .true.
    alloc(input_unit) = .true.
    alloc(output_unit) = .true.

  end subroutine lualloc_init

  !> Initialize the LU array, set default I/O values and interpret the
  !> command-line options. On output, optv contains the dash-options
  !> passed to critic, ghome is the path passed with -r and uroot is
  !> the root for the run.
  module subroutine stdargs(optv,ghome,uroot)
    use iso_fortran_env, only: input_unit, output_unit
    use param, only: dirsep
    character(len=:), allocatable, intent(out) :: optv
    character(len=:), allocatable, intent(out) :: ghome
    character(len=:), allocatable, intent(out) :: uroot

    integer, parameter :: arglen = 1024

    integer :: n
    integer :: argc, idx
    character(len=arglen) :: argv, aux
    logical :: local

    ! initialize the alloc array
    call lualloc_init()

    ! default values
    local = .false.
    optv=""
    ghome=""
    uin = input_unit
    uout = output_unit
    interactive = .true.
    usegui = .false.
    ucopy = uout
    uroot = "stdin"
    filepath = "."

    ! check if the -g argument is present -> usegui
    argc = command_argument_count()
    do n = 1, argc
       call getarg(n,argv)
       argv = adjustl(argv)
       if (len(argv) >= 2) then
          if (argv(1:2) == "-g") then
             usegui = .true.
             exit
          end if
       end if
    end do

    ! process arguments
    if (argc > 0) then
       n=0
       do while (n < argc)
          n=n+1
          call getarg(n,argv)
          argv = adjustl(argv)

          if (argv(1:1) == "-" .and. len(argv) > 1) then
             if (trim(argv(2:2)) == "r") then
                n = n + 1
                call getarg(n,argv)
                ghome = trim(adjustl(argv))
             elseif (trim(argv(2:2)) == "l") then
                local = .true.
             else
                optv = trim(optv) // argv(2:2)
             endif
          elseif (uin == input_unit .and..not.usegui) then
             uroot = trim(argv)
             if (index(uroot,dirsep) > 0) then
                filepath = uroot(1:index(uroot,dirsep,.true.)-1)
                if (len_trim(filepath) == 0 .or. local) filepath = "."
                ! xx(note1)xx
                ! gfortran-5 and gfortran-6 have problems with the much simpler:
                !   uroot = uroot(index(uroot,dirsep,.true.)+1:)
                ! where they bail out with internal compiler error. Strangely,
                ! gfortran-4 works perfectly. -- A.
                idx = index(uroot,dirsep,.true.)
                aux = uroot(idx+1:)
                uroot = trim(aux)
             end if
             if (index(uroot,'.') > 0) then
                ! see xx(note1)xx, tools_io.f90
                idx = index(uroot,'.',.true.)
                aux = uroot(1:idx-1)
                uroot = trim(aux)
             end if
             uroot = trim(adjustl(uroot))
             uin = fopen_read(argv,abspath0=.true.)
             interactive = .false.
          elseif (uout == output_unit .and..not.usegui) then
             uout = fopen_write(argv)
          endif
       end do
    end if

    ! the gui imposes a few restrictions
    if (usegui) then
       interactive = .true.
       ucopy = -1
       uroot = "gui"
    end if

  end subroutine stdargs

  !> Build a string from an integer*4. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the minimal representation for the integer, then the
  !> smallest possible length is used instead. Justify selects the
  !> justification of the integer in the field (left, center, and
  !> right). pad0, pad the full length with zeros. padspace, pad with
  !> this number of spaces. If padspace > 0, pad in the left, if < 0,
  !> in the right.
  module function string_int4(a,length,justify,pad0,padspace) result(s)
    character(len=:), allocatable :: s
    integer*4, intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    logical, intent(in), optional :: pad0
    integer, intent(in), optional :: padspace

    logical :: ispad
    character*(maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ialen0

    if (present(pad0)) then
       ispad = pad0
    else
       ispad = .false.
    end if
    ialen = 0
    if (present(length)) then
       if (length > 0) then
          ialen = length
       end if
    end if
    write (aux,'(I0)') a
    ialen0 = len_trim(aux)
    if (ialen == 0) then
       ialen = ialen0
    else
       ialen = max(ialen,ialen0)
       if (ispad) then
          write (fmt,'("(I0.",I0,")")') ialen
          write (aux,fmt) a
       else
          write (aux,'(I0)') a
       end if
    end if
    aux = adjustl(aux)
    s = aux(1:ialen)

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       elseif (justify == ioj_center) then
          ialen = len_trim(s)
          if (len(s) /= ialen) then
             i = (len(s) - ialen)/2 + 1
             s(i:i+ialen) = trim(s)
             s(1:i-1) = ""
             s(i+ialen+1:) = ""
          end if
       endif
    end if

    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          ! see xx(note1)xx, tools_io.f90
          aux = s(1:len(s)-ipad)
          s(ipad+1:) = trim(adjustl(aux))
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          ! see xx(note1)xx, tools_io.f90
          aux = s(ipad+1:)
          s(1:ialen) = trim(adjustl(aux))
          s(ialen+1:) = ""
       endif
    endif

  end function string_int4

  !> Build a string from an integer*8. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the minimal representation for the integer, then the
  !> smallest possible length is used instead. Justify selects the
  !> justification of the integer in the field (left, center, and
  !> right). pad0, pad the full length with zeros. padspace, pad with
  !> this number of spaces. If padspace > 0, pad in the left, if < 0,
  !> in the right.
  module function string_int8(a,length,justify,pad0,padspace) result(s)
    character(len=:), allocatable :: s
    integer*8, intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    logical, intent(in), optional :: pad0
    integer, intent(in), optional :: padspace

    logical :: ispad
    character*(maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ialen0

    if (present(pad0)) then
       ispad = pad0
    else
       ispad = .false.
    end if
    ialen = 0
    if (present(length)) then
       if (length > 0) then
          ialen = length
       end if
    end if
    write (aux,'(I0)') a
    ialen0 = len_trim(aux)
    if (ialen == 0) then
       ialen = ialen0
    else
       ialen = max(ialen,ialen0)
       if (ispad) then
          write (fmt,'("(I0.",I0,")")') ialen
          write (aux,fmt) a
       else
          write (aux,'(I0)') a
       end if
    end if
    aux = adjustl(aux)
    s = aux(1:ialen)

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       elseif (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:i+ialen) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       endif
    end if

    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          ! see xx(note1)xx, tools_io.f90
          aux = s(1:len(s)-ipad)
          s(ipad+1:) = trim(adjustl(aux))
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          ! see xx(note1)xx, tools_io.f90
          aux = s(ipad+1:)
          s(1:ialen) = trim(adjustl(aux))
          s(ialen+1:) = ""
       endif
    endif

  end function string_int8

  !> Build a string from a real*8 number. desc is the descriptor used
  !> (f, g, e, d, etc.). length is the length of the output string. If
  !> the length given by the user is not enough to contain the minimal
  !> representation for the number, then the smallest possible length is
  !> used instead. decimal is the number of decimal places (default:
  !> 7). Justify selects the justification of the integer in the field
  !> (left, center, and right). If a postive value is given to justify,
  !> then align the point to that character. padspace, pad with this
  !> number of spaces. If padspace > 0, pad in the left, if < 0, in the
  !> right.
  module function string_real8(a,desc,length,decimal,justify,padspace) result(s)
    character(len=:), allocatable :: s
    real*8, intent(in) :: a
    character*1, intent(in) :: desc
    integer, intent(in), optional :: length
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace

    integer, parameter :: decimal_default = 7

    character*(maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ipad0, idx, dec
    logical :: islen, first

    islen = .false.
    if (present(length)) islen = (length > 0)
    if (present(decimal)) then
       dec = decimal
    else
       dec = decimal_default
    endif

    ! We'll need two passes through here if length is present but
    ! decimal is not. The first is to measure the length with the
    ! default number of decimal places. The second sets the number of
    ! decimals to match the requested length.
    first = .true.
10  continue
    if (desc == 'e' .or. desc == 'E' .or. desc == 'd' .or. desc == 'D') then
       write (fmt,'("(1p,",A1,I0,".",I0,")")') desc, maxlen, dec
    else
       write (fmt,'("(",A1,I0,".",I0,")")') desc, maxlen, dec
    end if
    write (aux,fmt) a
    aux = adjustl(aux)
    ialen = len_trim(aux)
    if (islen .and..not.present(decimal)) then
       if (ialen > length .and. first) then
          ! Calculate the number of decimals to get that length, then
          ! loop back and write the actual string.
          dec = max(decimal_default - (ialen-length),1)
          first = .false.
          goto 10
       end if
    end if

    if (islen) ialen = max(length,ialen)
    s = aux(1:ialen)

    ipad0 = 0
    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       elseif (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:i+ialen) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       elseif (justify > 0) then
          idx = index(s,".")
          if (idx > 0) then
             ipad0 = max(justify - idx,0)
          endif
       endif
    end if

    if (present(padspace)) then
       if (ipad0 == 0) then
          ipad0 = padspace
       else
          ipad0 = ipad0 + padspace
       end if
    endif
    if (ipad0 > 0) then
       do ipad = ipad0,0,-1
          if (s(len(s)-ipad+1:) == "") exit
       end do
       s(ipad+1:) = s(1:len(s)-ipad)
       s(1:ipad) = ""
    elseif (ipad0 < 0) then
       do ipad = -ipad0,0,-1
          if (s(1:ipad) == "") exit
       end do
       ialen = len(s)-ipad
       s(1:ialen) = s(ipad+1:)
       s(ialen+1:) = ""
    endif

  end function string_real8

  !> Use the real*8 routine for real*4.
  module function string_real4(a,desc,length,decimal,justify,padspace) result(s)
    character(len=:), allocatable :: s
    real*4, intent(in) :: a
    character*1, intent(in) :: desc
    integer, intent(in), optional :: length
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace

    s = string_real8(real(a,8),desc,length,decimal,justify,padspace)

  end function string_real4

  !> Build a string from another string. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the original string, then it is ignored. Justify selects
  !> the justification of the integer in the field (left, center, and
  !> right). padspace, pad with this number of spaces. If padspace > 0,
  !> pad in the left, if < 0, in the right.
  module function string_char(a,length,justify,padspace) result(s)
    character(len=:), allocatable :: s
    character(len=*), intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace

    integer :: i, ialen, ipad, i2len

    i2len = len_trim(adjustl(a))
    ialen = i2len
    if (present(length)) then
       if (length > ialen) then
          ialen = length
       end if
    end if
    allocate(character(ialen)::s)
    s(1:i2len) = trim(adjustl(a))
    s(i2len+1:ialen) = ""

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       elseif (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:min(i+ialen,len(s))) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       endif
    end if

    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          s(ipad+1:) = s(1:len(s)-ipad)
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          s(1:ialen) = s(ipad+1:)
          s(ialen+1:) = ""
       endif
    endif

  end function string_char

  !> Build a string from a logical value. Returns "Yes"
  !> if input is true or "No" if it is false.
  module function string_logical(a) result(s)
    character(len=:), allocatable :: s
    logical, intent(in) :: a

    if (a) then
       s = "Yes"
    else
       s = "No"
    end if

  end function string_logical

  !> Compare two strings for equality.
  module function equal(s,t)
    character*(*), intent(in) :: s !< First string
    character*(*), intent(in) :: t !< Second string
    logical :: equal

    equal = string(s) == string(t)

  end function equal

  !> Compare two strings for equality (case insensitive).
  module function equali(s,t)
    character*(*), intent(in) :: s !< First string
    character*(*), intent(in) :: t !< Second string
    logical :: equali

    equali = lower(string(s)) == lower(string(t))

  end function equali

  !> Read a line from logical unit u, and return true if read was
  !> successful. Continuation with "\", skip blank lines and comments.
  !> If eofstop is true, raise error on EOF. If ucopy (integer) exists
  !> and is positive, write a copy of the output line to that logical
  !> unit, preceded by a prefix. If nprompt, show this number in the
  !> prompt.
  !> If lstr, read the line from the instr global variable and ignore
  !> the logical unit u. If lstr is true, eofstop is always .false.
  !> ucopy is not used (-1).
  module function getline(u,lstr,oline,eofstop,ucopy,nprompt)
    character(len=:), allocatable, intent(out) :: oline
    integer, intent(in) :: u
    logical, intent(in) :: lstr
    logical, intent(in), optional :: eofstop
    integer, intent(in), optional :: ucopy
    integer, intent(in), optional :: nprompt
    logical :: getline

    integer :: i, lenu
    character(len=:), allocatable :: line
    logical :: notfirst, ok
    character*(1) :: inquote

    character*(2), parameter :: prfx = '%%'

    getline = .false.
    oline = ""
    notfirst = .false.
    do while (.true.)
       ! read the line
       if (.not.lstr) then
          ok = getline_raw(u,line,nprompt=nprompt)
       else
          instrptr = instrptr + 1
          ok = (instrptr <= size(instr,1))
          if (ok) then
             line = instr(instrptr)%s
          else
             line = ""
          end if
       end if
       ! exit if eof
       if (.not.ok) then
          if (present(eofstop).and..not.lstr) then
             if (eofstop) call ferror("getline","unexpected end of file",faterr)
          end if
          return
       end if

       ! remove tabs
       do i = 1, len(line)
          if (line(i:i) == tab) line(i:i) = " "
       end do

       ! remove trailing blanks
       line = adjustl(line)
       lenu = len_trim(line)

       ! remove comments outside quotes
       inquote = ""
       do i = 1, lenu
          if (line(i:i) == "'".or.line(i:i) == '"') then
             if (len_trim(inquote) == 0) then
                inquote = line(i:i)
                cycle
             elseif (line(i:i) == inquote) then
                inquote = ""
                cycle
             endif
          endif
          if (len_trim(inquote) > 0) cycle
          if (line(i:i) == "#") then
             line = line(1:i-1)
             exit
          endif
       end do

       ! skip blank lines
       lenu = len_trim(line)
       if (lenu == 0) then
          if (notfirst) then
             exit
          else
             cycle
          end if
       end if

       ! continuation
#ifdef __PGI
       if (line(lenu:lenu) /= "\\") then
#else
       if (line(lenu:lenu) /= "\") then ! Hey, emacs, turn the lights on! " Thanks.
#endif
          oline = trim(oline) // " " // trim(line)
          oline = adjustl(oline)
          exit
       end if
       oline = trim(oline) // " " // line(1:lenu-1)
       notfirst = .true.
    end do
    if (present(ucopy).and..not.lstr) then
       if (ucopy >= 0) &
          write (ucopy,'(A,X,A/)') prfx, trim(oline)
    endif
    getline = .true.

  end function getline

  !> Read a line from logical unit u, and return true if read was
  !> successful. Don't process it to remove comments, etc. If eofstop
  !> is true, raise error on EOF. If compiled with readline support
  !> and u is the standard input, use readline to fetch the line. In
  !> that case, eofstop is ignored.
  module function getline_raw(u,line,eofstop,nprompt) result(ok)
    use iso_fortran_env, only: input_unit
#ifdef HAVE_READLINE
    use, intrinsic :: iso_c_binding, only: c_ptr, c_associated, c_size_t
    use c_interface_module, only: f_c_string_dup, c_f_string, c_strlen_safe, c_free
#endif
    integer, intent(in) :: u
    character(len=:), allocatable, intent(out) :: line
    logical, intent(in), optional :: eofstop
    integer, intent(in), optional :: nprompt
    logical :: ok

#ifdef HAVE_READLINE
    interface
       function readline(prompt) bind(c)
         import c_ptr
         type(c_ptr), value :: prompt
         type(c_ptr) :: readline
       end function readline
       subroutine add_history(line) bind(c)
         import c_ptr
         type(c_ptr), value :: line
       end subroutine add_history
    end interface

    type(c_ptr) :: line_c, prompt_c
    integer(c_size_t) :: ll
#endif
    character(len=:), allocatable :: aux
    integer :: ios, nread, i

    integer, parameter :: blocksize = 128

    ! use readline if stdin and readline is available
#ifdef HAVE_READLINE
    if (u == input_unit) then
       ! use readline to fetch a new line from the user
       if (present(nprompt)) then
          prompt_c = f_c_string_dup("critic2:"//string(nprompt)//"> ")
       else
          prompt_c = f_c_string_dup("critic2> ")
       end if
       line_c = readline(prompt_c)
       call c_free(prompt_c)

       ok = c_associated(line_c)
       if (ok) then
          ! transform to fortran
          ll = c_strlen_safe(line_c)
          allocate(character(len=ll) :: line)
          call c_f_string(line_c,line)

          ! save the line in the history
          if (ll > 0) call add_history(line_c)

          ! clean up memory and exit
          call c_free(line_c)
       else
          line = ""
       end if
       return
    end if
#endif

    ! write the prompt
    if (u == input_unit) then
       if (present(nprompt)) then
          write (*,'(A)',advance='no') "critic2:"//string(nprompt)//"> "
       else
          write (*,'(A)',advance='no') "critic2> "
       end if
    end if

    ! normal line reading
    allocate(character(len=blocksize)::aux)
    do i = 1,blocksize
       aux(i:i) = " "
    end do
    line = ""
    ok = .true.
    do while(.true.)
       read(u,'(A)',advance="no",iostat=ios,size=nread) aux
       line = line(1:len(line)) // aux
       ok = .not.is_iostat_end(ios)
       if (is_iostat_eor(ios) .or. is_iostat_end(ios)) exit
    end do
    if (.not.ok.and.present(eofstop)) then
       if (eofstop) call ferror("getline_raw","unexpected end of file",faterr)
    end if

    ! insert a blank line if eof was received in stdin
    if (u == input_unit .and..not.ok) write (*,*)

  end function getline_raw

  !> Guess the atomic number from the atomic name. A value of -1 is
  !> returned if the atom is not recognized.  A valid atomic name is
  !> made of one or two letters corresponding to a legal atomic symbol
  !> followed by any other additional characters.
  module function zatguess(atname)
    use param, only: maxzat0
    integer :: zatguess !< Ouptut atomic number
    character*(*), intent(in) :: atname !< Input atomic symbol (case insensitive)

    character*(len(atname)) :: aux
    character*(2) :: atsymbol
    integer :: i
    character*(2), parameter :: an(1:maxzat0) =(/&
       'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE',& ! 1-10
       'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA',& ! 11-20
       'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',& ! 21-30
       'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR',& ! 31-40
       'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',& ! 41-50
       'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',& ! 51-60
       'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',& ! 61-70
       'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',& ! 71-80
       'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',& ! 81-90
       'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM',& ! 91-100
       'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS',& ! 101-110
       'RG', 'CN', 'NH', 'FL', 'MC', 'LV', 'TS', 'OG',& ! 111-118
       'XN', 'XB', 'XR', 'XC', 'XZ'& ! 119-123
       /)

    aux = adjustl(atname)
    atsymbol = aux(1:min(2,len(aux)))
    atsymbol = upper(atsymbol)
    if (.not.isletter(atsymbol(2:2))) atsymbol(2:2) = ' '

1   continue
    ! the normal elements
    do i = 1, size(an)
       if (atsymbol == an(i)) then
          zatguess = i
          return
       endif
    enddo

    ! deuterium, ghost atoms, etc.
    zatguess = -1
    if (atsymbol == 'D ') zatguess = 1
    if (atsymbol == 'GH') zatguess = 0
    if (atsymbol == 'XX') zatguess = 0
    if (atsymbol == 'BQ') zatguess = 0

    ! If the second character was junk, give it another go with just the first letter
    if (zatguess == -1 .and. atsymbol(2:2) /= ' ') then
       atsymbol(2:2) = ' '
       goto 1
    end if

  end function zatguess

  !> The reverse of zatguess. Given the atomic number, return the
  !> chemical symbol of the atom, or XX if it is not in range.
  !> If nounderscore is true, use blanks instead of underscores to pad
  !> the symbol.
  module function nameguess(zat,nounderscore)
    use param, only: maxzat0
    integer, intent(in) :: zat !< Input atomic number
    logical, intent(in), optional :: nounderscore !< Use blanks instead of underscore to fill
    character*(2) :: nameguess !< Output atomic symbol

    character*(2), parameter :: an(1:maxzat0) =(/&
       'H_', 'He', 'Li', 'Be', 'B_', 'C_', 'N_', 'O_', 'F_', 'Ne',& ! 1-10
       'Na', 'Mg', 'Al', 'Si', 'P_', 'S_', 'Cl', 'Ar', 'K_', 'Ca',& ! 11-20
       'Sc', 'Ti', 'V_', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',& ! 21-30
       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y_', 'Zr',& ! 31-40
       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',& ! 41-50
       'Sb', 'Te', 'I_', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',& ! 51-60
       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',& ! 61-70
       'Lu', 'Hf', 'Ta', 'W_', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',& ! 71-80
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',& ! 81-90
       'Pa', 'U_', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',& ! 91-100
       'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',& ! 101-110
       'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',& ! 111-118
       'Xn', 'Xb', 'Xr', 'Xc', 'Xz'& ! 119-123
       /)

    if (zat == 0) then
       nameguess = 'Bq'
    elseif (zat > 0 .and. zat <= maxzat0) then
       nameguess = an(zat)
    elseif (zat == 201) then
       nameguess = 'Xx'
    end if
    if (present(nounderscore)) then
       if (nounderscore) then
          if (nameguess(2:2) == "_") nameguess(2:2) = " "
       end if
    end if

  end function nameguess

  !> Convert a string to uppercase, portable.
  module function upper(a)
    character*(*), intent(in) :: a
    character*(len(a)) :: upper

    character(*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, idx

    upper = a
    do i = 1,len(a)
       idx = index(lo,a(i:i))
       if (idx > 0) upper(i:i) = up(idx:idx)
    end do

  end function upper

  !> Convert a string to lowercase, portable.
  module function lower(a)
    character*(*), intent(in) :: a
    character*(len(a)) :: lower

    character(*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, idx

    lower = a
    do i = 1,len(a)
       idx = index(up,a(i:i))
       if (idx > 0) lower(i:i) = lo(idx:idx)
    end do

  end function lower

  !> Remove blanks from a string. Optionally, convert
  !> to uppercase or lowercase.
  module function deblank(a,toup,todn)
    character*(*), intent(in) :: a
    logical, intent(in), optional :: toup
    logical, intent(in), optional :: todn
    character*(len(a)) :: deblank

    character(*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, j, idx

    deblank = a
    j = 0
    do i = 1, len(a)
       if (a(i:i) == ' ') cycle
       j = j + 1
       idx = 0
       if (present(toup)) then
          if (toup) then
             idx = index(lo,a(i:i))
             if (idx > 0) deblank(j:j) = up(idx:idx)
          end if
       end if
       if (present(todn)) then
          if (todn) then
             idx = index(up,a(i:i))
             if (idx > 0) deblank(j:j) = lo(idx:idx)
          end if
       end if
       if (idx == 0) deblank(j:j) = a(i:i)
    end do
    if (j < len(a)) deblank(j+1:len(a)) = ""

  end function deblank

  !> Strip a character from a string. Return the stripped string.
  module function stripchar(a,ch)
    character*(*), intent(in) :: a
    character, intent(in) :: ch
    character*(len(a)) :: stripchar

    integer :: i, j

    stripchar = a
    if (index(a,ch) == 0) return

    j = 0
    do i = 1, len(a)
       if (a(i:i) == ch) cycle
       j = j + 1
       stripchar(j:j) = a(i:i)
    end do
    if (j < len(a)) stripchar(j+1:len(a)) = ""

  end function stripchar

  !> Get next word from line at lp. If successful, increase lp. A word
  !> is defined as any sequence on nonblanks.
  module function getword(line,lp) result(word)
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp !< Pointer to position on input line, updated after reading.
    character(len=:), allocatable :: word !< Word read

    integer :: len0, wp, i

    len0 = len(line)
    word = ""
    if (lp > len0) return

    i = lp
    do while (line(i:i) == " ")
       i = i+1
       if (i > len0) return
    enddo

    wp = i
    do while (i<len0 .and. line(i:i)/=" ")
       i = i+1
    enddo
    word = line(wp:i)
    lp = i+1

  end function getword

  ! Same as getword, but convert to lowercase
  module function lgetword (line,lp) result(word)
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp !< Pointer to position on input line, updated after reading.
    character(len=:), allocatable :: word

    word = getword(line,lp)
    if (len_trim(word) > 0) word = lower(word)

  endfunction lgetword

  !> True if the input is a letter.
  module function isletter(a)
    character*(1), intent(in) :: a
    logical :: isletter

    character(*), parameter :: letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    isletter = (index(letters,a) > 0)

  end function isletter

  !> True if the input is a digit.
  module function isdigit(a)
    character*(1), intent(in) :: a
    logical :: isdigit

    character(*), parameter :: digits = '0123456789'

    isdigit = (index(digits,a) > 0)

  end function isdigit

  !> Get an integer value from input text and return it in
  !> ival. Advance the pointer (lp0) if provided. If lp0 is not
  !> present, use 1. If a valid integer is not found, returns .false.
  module function isinteger(ival,line,lp0)
    logical :: isinteger
    integer, intent(inout) :: ival !< Integer value read, if not a valid integer is found, it is not modified
    character*(*), intent(in) :: line !< Input string
    integer, intent(inout), optional :: lp0 !< Pointer to current position on string

    integer :: i, len0, lp
    logical :: readit

    lp = 1
    if (present(lp0)) lp = lp0

    isinteger = .false.
    len0 = len(line)
    if (lp > len0) goto 999
    i = lp
    do while (line(i:i)==" " .or. line(i:i)==tab)
       i=i+1
       if (i > len0) goto 999
    enddo
    if (line(i:i) == '+' .or. line(i:i) == '-') i=i+1
    if (i > len0) goto 999
    if (isdigit(line(i:i))) then
       do while (isdigit(line(i:i)))
          if (i >= len0) exit
          i=i+1
       enddo
       readit = (i == len0)
       if (.not.readit) readit = (line(i:i) == " ")
       if (readit) then
          read (line(lp:i),*) ival
          lp = i+1
          isinteger=.true.
       endif
    endif

999 continue
    if (present(lp0)) lp0 = lp

  end function isinteger

  !> Read a list of integers from a line, starting at pointer lp0. Returns
  !> the number of integers read (n), the integers (ix), and updates the
  !> position of the pointer.
  module subroutine readintegers(n,ix,line,lp0)
    use types, only: realloc
    integer, intent(out) :: n
    integer, intent(inout), allocatable :: ix(:)
    character*(*), intent(in) :: line
    integer, intent(inout), optional :: lp0

    integer :: lp

    if (present(lp0)) then
       lp = lp0
    else
       lp = 1
    end if
    if (allocated(ix)) deallocate(ix)
    allocate(ix(1))

    n = 0
    do while(isinteger(ix(n+1),line,lp))
       n = n + 1
       if (n >= size(ix,1)) call realloc(ix,2*n)
    end do
    if (n > 0) then
       call realloc(ix,n)
    else
       deallocate(ix)
    end if
    if (present(lp0)) lp0 = lp

  end subroutine readintegers

  !> Get a real number from line and sets rval to it. If a valid real
  !> number is not found, isreal returns .false. and rval is not changed.
  module function isreal(rval, line, lp0)
    real*8, intent(out) :: rval !< Real value read
    character*(*), intent(in) :: line !< Input string
    integer, intent(inout), optional :: lp0 !< Pointer to current position on string
    logical :: isreal

    integer :: len0, lp
    logical :: matched

    character(*), parameter :: exponents = 'eEdDqQ+-'
    character(*), parameter :: lexponents = 'eEdDqQ'

    integer :: i

    ! initialize
    lp = 1
    if (present(lp0)) lp = lp0

    isreal = .false.
    len0 = len(line)
    if (lp > len0) return

    ! trailing blanks
    i = lp
    do while (line(i:i)==" " .or. line(i:i)==tab)
       i=i+1
       if (i > len0) return
    enddo

    ! sign
    matched = .false.
    if (line(i:i) == '+' .or. line(i:i) == '-') i = i + 1
    if (i > len0) return

    ! first sequence of digits
    if (isdigit(line(i:i))) then
       matched = .true.
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       enddo
    end if

    ! decimal point
    if (line(i:i) == '.') then
       i=i+1
       if (i > len0) goto 10
    end if

    ! second sequence of digits
    if (isdigit(line(i:i))) then
       matched = .true.
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       enddo
    end if

    ! that did not look of a real number at all!
    if (.not.matched) return

    ! get optional exponent
    if (index(exponents,line(i:i)) > 0) then
       i = i+1
       if (i > len0) return

       if (line(i:i)=='+' .or. line(i:i)=='-') then
          if (index(lexponents,line(i-1:i-1)) > 0) then
             i = i+1
             if (i > len0) return
          else
             return
          endif
       end if

       if (.not.isdigit(line(i:i))) return
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       enddo

       if (i < len0 .and..not.line(i:i) == " ") return
    endif

10  continue
    read (line(lp:i-1),*) rval
    if (present(lp0)) lp0 = i
    isreal = .true.

  end function isreal

  !> Get a expression from the line. Expressions start with " or '.
  !> If the read fails, then lp0 is unchanged, and isexpression is
  !> false.
  module function isexpression(word, line, lp0)
    character(len=:), allocatable, intent(out) :: word !< Expression read
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    logical :: isexpression !< True if the read was successful

    integer :: l, lp, lpi
    character*1 :: sclose

    isexpression = .false.
    lp = lp0
    l = len(line)
    if (lp > l) return
    do while (line(lp:lp)==" ")   !skip blanks
       lp = lp + 1
       if (lp > l) return
    enddo

    if (line(lp:lp) /= '"' .and. line(lp:lp) /= "'") return
    sclose = line(lp:lp)
    lp = lp + 1
    lpi = lp
    if (lp > l) return
    do while (line(lp:lp) /= sclose)
       lp = lp + 1
       if (lp > l) return
    end do
    if (line(lp:lp) /= sclose) return

    isexpression = .true.
    word = line(lpi:lp-1)
    lp0 = lp+1

  end function isexpression

  !> Get a expression from line. Expressions start with " or '. If no
  !> expression is available, then get the next word.
  module function isexpression_or_word(word, line, lp0)
    character(len=:), allocatable, intent(out) :: word !< Word/expression read
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    logical :: isexpression_or_word !< True if the read was successful

    character(len=:), allocatable :: word0

    ! is expression...
    isexpression_or_word = isexpression(word,line,lp0)
    if (isexpression_or_word) return

    ! ... or word?
    word0 = getword(line,lp0)
    if (len_trim(word0) > 0) then
       word = word0
       isexpression_or_word = .true.
    else
       isexpression_or_word = .false.
    end if

  end function isexpression_or_word

  !> Read an assignment from the input line
  module function isassignment(var, expr, line)
    use param, only: constlist, funclist
    character*(*), intent(in) :: line !< Input line
    character(len=:), allocatable, intent(out) :: var !< Output variable name (left hand side)
    character(len=:), allocatable, intent(out) :: expr !< Output expression (right hand side)
    logical :: isassignment !< True if the read was successful

    integer :: idx, i
    logical :: nok

    ! get the = in the middle
    var = ""
    expr = ""
    isassignment = .false.
    idx = index(line,'=')
    if (idx == 0) return

    ! assign variable name and expression
    var = trim(adjustl(line(1:idx-1)))
    expr = trim(adjustl(line(idx+1:)))

    ! check that there are no more =
    if (index(var,'=') /= 0) return
    if (index(expr,'=') /= 0) return

    ! check that the variable name is only one word
    if (index(trim(var),' ') /= 0) return

    ! and that it starts with a letter
    if (.not.isletter(var(1:1))) return

    ! and also that it doesn't have any funny symbols
    do i = 1, len_trim(var)
       nok = .not.isletter(var(i:i))
       nok = nok .and. .not.isdigit(var(i:i))
       nok = nok .and. (var(i:i) /= '_')
       if (nok) return
    end do

    ! check that the variable name is not a constant
    do i = 1, size(constlist)
       if (trim(constlist(i)) == trim(lower(var))) return
    end do

    ! check that the variable name is not a function
    do i = 1, size(funclist)
       if (trim(funclist(i)) == trim(lower(var))) return
    end do

    isassignment = .true.
  end function isassignment

  !> Initialize file system and connect standard units.
  !> Use this before falloc and fdealloc.
  module subroutine history_init()
#ifdef HAVE_READLINE
    use iso_c_binding, only: c_int, c_ptr
    use c_interface_module, only: f_c_string_dup, c_free
    interface
       function read_history(file) bind(c)
         import c_int, c_ptr
         integer(c_int) :: read_history
         type(c_ptr), value :: file
       end function read_history
    end interface

    character(len=1024) :: home
    type(c_ptr) :: file_c
    integer(c_int) :: istat

    ! This probably works only for linux, but it doesn't crash if the
    ! read fails
    if (interactive) then
       call get_environment_variable("HOME",home,status=istat)
       if (istat == 0) then
          file_c = f_c_string_dup(trim(home) // "/.critic2_history")
          istat = read_history(file_c)
          call c_free(file_c)
       end if
    end if
#endif
  end subroutine history_init

  !> If readline is used and the run is interactive, write the history file.
  module subroutine history_end()
#ifdef HAVE_READLINE
    use iso_c_binding, only: c_int, c_ptr
    use c_interface_module, only: f_c_string_dup, c_free
    interface
       function write_history(file) bind(c)
         import c_int, c_ptr
         integer(c_int) :: write_history
         type(c_ptr), value :: file
       end function write_history
    end interface

    character(len=1024) :: home
    type(c_ptr) :: file_c
    integer(c_int) :: istat

    ! This probably works only for linux, but it doesn't crash if the
    ! write fails
    if (interactive) then
       call get_environment_variable("HOME",home,status=istat)
       if (istat == 0) then
          file_c = f_c_string_dup(trim(home) // "/.critic2_history")
          istat = write_history(file_c)
          call c_free(file_c)
       end if
    end if
#endif
  end subroutine history_end

  !> Open a file for reading. The argument form controls the
  !> formatting, and is passed directly to open(). If abspath is
  !> present, file in input is as an absolute path. If errstop is
  !> present, stop the execution if the file could not be opened;
  !> otherwise, return a negative lu. Default is errstop = .true.
  module function fopen_read(file,form,abspath0,errstop) result(lu)
    use param, only: dirsep
    character*(*), intent(in) :: file
    character*(*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    logical, intent(in), optional :: errstop
    integer :: lu

    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath, errstop_

    errstop_ = .true.
    if (present(errstop)) errstop_ = errstop
    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='old',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='old',iostat=ios)
    endif
    if (ios /= 0) then
       lu = -1
       if (errstop_) &
          call ferror("fopen_read","error opening file: "//string(file),faterr)
    end if

  end function fopen_read

  !> Open a file for reading. The argument form controls the
  !> formatting, and is passed directly to open(). If abspath is
  !> present, file in input is as an absolute path.
  module function fopen_write(file,form,abspath0,errstop) result(lu)
    use param, only: dirsep
    character*(*), intent(in) :: file
    character*(*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    logical, intent(in), optional :: errstop
    integer :: lu

    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath, errstop_

    errstop_ = .true.
    if (present(errstop)) errstop_ = errstop
    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='unknown',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='unknown',iostat=ios)
    end if
    if (ios /= 0) then
       lu = -1
       if (errstop_) &
          call ferror("fopen_write","error opening file: "//string(file),faterr)
    end if

  end function fopen_write

  !> Open a file for appending
  module function fopen_append(file,form,abspath0,errstop) result(lu)
    use param, only: dirsep
    character*(*), intent(in) :: file
    character*(*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    logical, intent(in), optional :: errstop
    integer :: lu

    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath, errstop_

    errstop_ = .true.
    if (present(errstop)) errstop_ = errstop
    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='old',access='append',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='old',access='append',iostat=ios)
    end if
    if (ios /= 0) then
       lu = -1
       if (errstop_) &
          call ferror("fopen_append","error opening file: "//string(file),faterr)
    end if

  end function fopen_append

  !> Open a scratch file for writing
  module function fopen_scratch(form,errstop) result(lu)
    character*(*), intent(in), optional :: form
    logical, intent(in), optional :: errstop
    integer :: lu

    integer :: ios
    logical :: errstop_

    errstop_ = .true.
    if (present(errstop)) errstop_ = errstop
    lu = falloc()
    if (present(form)) then
       if (lower(form) == "unformatted") then
          open(unit=lu,status='scratch',form=form,access="stream",iostat=ios)
       else
          open(unit=lu,status='scratch',form=form,iostat=ios)
       end if
    else
       open(unit=lu,status='scratch',form="unformatted",access="stream",iostat=ios)
    end if
    if (ios /= 0) then
       lu = -1
       if (errstop_) &
          call ferror("fopen_write","error opening scratch file",faterr)
    end if

  end function fopen_scratch

  !> Close a file
  module subroutine fclose(lu)
    integer, intent(in) :: lu

    close(lu)
    call fdealloc(lu)

  end subroutine fclose

  !> Allocate a logical unit
  module function falloc()
    integer :: falloc

    do falloc = 1, size(alloc)
       if (.not.alloc(falloc)) exit
    end do
    if (alloc(falloc)) &
       call ferror("falloc","could not allocate logical unit",faterr)
    alloc(falloc) = .true.

  end function falloc

  !> Deallocate a logical unit
  module subroutine fdealloc(u)
    integer, intent(in) :: u !< LU to deallocate

    if (u < 0 .or. u > size(alloc)) &
       call ferror("fdealloc","could not deallocate logical unit",faterr)
    alloc(u)=.false.

  end subroutine fdealloc

  !> Send an error message 'message' to stdout, coming from routine
  !> 'routine'. errortype is the error code (see mod_param.f90).
  module subroutine ferror(routine,message,errortype,inputline,syntax)
    character*(*), intent(in) :: routine !< Routine calling the error
    character*(*), intent(in) :: message !< The message
    integer, intent(in) :: errortype !< Fatal, warning or info
    character*(*), intent(in), optional :: inputline !< Previous-line input message
    logical, intent(in), optional :: syntax !< Is this a syntax error? (do not stop if interactive)

    character*(20) :: chtype

    ! message styles shamelessly copied from abinit.
    if (errortype.eq.faterr) then
       chtype='ERROR'
    else if (errortype.eq.warning) then
       chtype='WARNING'
    else if (errortype.eq.noerr) then
       chtype='COMMENT'
    else
       chtype='UNKNOWN'
    endif
    !$omp critical (IO)
    if (present(inputline)) &
       write (uout,'("!error! ",A)') trim(inputline)
    write (uout,100) trim(chtype),trim(routine),trim(message)
    !$omp end critical (IO)
    if (errortype.eq.faterr) then
       !$omp critical (IO)
       write (uout,100) trim(chtype),trim(routine), trim(message)
       !$omp end critical (IO)
       if (present(syntax)) then
          if (syntax .and. interactive) return
       end if
       stop 1
    else if(errortype.eq.warning) then
       nwarns = nwarns + 1
    else if(errortype.eq.noerr) then
       ncomms = ncomms + 1
    endif

100 format (A,"(",A,"): ",A)

  end subroutine ferror

  !> Write a timestamp to the output, with an optional prefix message.
  module subroutine tictac(mesg)
    character*(*), intent(in) :: mesg !< Prefix message
    integer :: values(8)

    character(len=:), allocatable :: output

    call date_and_time(values=values)

    output = string(mesg)
    output = output // "--" // string(values(1))
    output = output // "/" // string(values(2))
    output = output // "/" // string(values(3))
    output = output // ", " // string(values(5),2,pad0=.true.)
    output = output // ":" // string(values(6),2,pad0=.true.)
    output = output // ":" // string(values(7),2,pad0=.true.)
    output = output // "." // string(values(8),3,pad0=.true.)
    write (uout,'(A)') output

  end subroutine tictac

  !> Start the clock
  module subroutine start_clock()

    call system_clock(stime0)
    call cpu_time(ctime0)

  end subroutine start_clock

  !> Print the elapsed time to output
  module subroutine print_clock()

    integer*8 :: irate, imax
    real :: ctime
    integer*8 :: stime
    character(len=:), allocatable :: sout

    call system_clock(count_rate=irate,count_max=imax)
    call cpu_time(ctime)
    call system_clock(stime)

    call sectotime(real(stime-stime0)/real(irate),s0=sout)
    write (uout,'("Elapsed wall time: ",A)') string(sout)
    call sectotime(ctime-ctime0,s0=sout)
    write (uout,'("Elapsed CPU time: ",A)') string(sout)

  end subroutine print_clock

  !> Seconds to days/hours/minutes/seconds format. x(1): seconds,
  !> x(2): minutes, x(3): hours, x(4): days.
  module subroutine sectotime(s,x0,s0)
    real, intent(in) :: s
    integer, intent(out), optional :: x0(4)
    character(len=:), allocatable, optional :: s0

    integer :: i
    integer, parameter :: f(3) = (/60, 60, 24/)
    integer :: x(4)
    character(len=:), allocatable :: saux

    x = 0
    x(1) = floor(s)
    do i = 1, 3
       x(i+1) = x(i) / f(i)
       x(i) = x(i) - f(i) * x(i+1)
    end do

    if (present(x0)) x0 = x
    if (present(s0)) then
       s0 = ""
       if (x(4) > 0) then
          saux = trim(s0) // " " // string(x(4)) // "d"
          s0 = saux
       end if
       if (x(3) > 0) then
          saux = trim(s0) // " " // string(x(3)) // "h"
          s0 = saux
       end if
       if (x(2) > 0) then
          saux = trim(s0) // " " // string(x(2)) // "m"
          s0 = saux
       end if
       saux = trim(s0) // " " // string(x(1)) // "s"
       s0 = saux
    end if

  end subroutine sectotime

end submodule proc
