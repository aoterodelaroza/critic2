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

submodule (global) proc
  implicit none

contains

  !> Initialize basic variables at the beginning of the run.
  !> Also sets default values by calling global_set_defaults.
  module subroutine global_init(ghome,datadir)
    use tools_io, only: string, ferror, warning, uout
    use param, only: dirsep
    character*(*) :: ghome, datadir
    integer :: isenv
    logical :: lchk
    character(len=:), allocatable :: cifstr, msgr1, msgr2, msg1, msg2, msg3
    integer, parameter :: maxlenpath = 1024

    cifstr = dirsep // "cif" // dirsep // "cif_core.dic"

    ! read the -r option
    msgr1 = ""
    msgr2 = ""
    if (len_trim(ghome) > 0) then
       critic_home = string(ghome)
       inquire(file=trim(critic_home) // cifstr,exist=lchk)
       if (lchk) goto 99
       msgr1 = "(!) 0. Not found (-r option): " // trim(critic_home) // cifstr

       critic_home = string(ghome) // dirsep // "dat"
       inquire(file=trim(critic_home) // cifstr,exist=lchk)
       if (lchk) goto 99
       msgr2 = "(!) 0. Not found (-r option): " // trim(critic_home) // cifstr
    endif

    ! read env variable CRITIC_HOME
    if (allocated(critic_home)) deallocate(critic_home)
    allocate(character(len=maxlenpath)::critic_home)
    call get_environment_variable("CRITIC_HOME",critic_home,status=isenv)
    if (isenv ==0) then
       critic_home = trim(critic_home) // dirsep // "dat"
       inquire(file=trim(critic_home) // cifstr,exist=lchk)
       if (lchk) goto 99
       msg1 = "(!) 1. Not found (CRITIC_HOME): " // trim(critic_home) // cifstr
    else
       msg1 = "(!) 1. CRITIC_HOME environment variable not set"
    end if

    ! then the install path
    critic_home = trim(adjustl(datadir))
    inquire(file=trim(critic_home) // cifstr,exist=lchk)
    if (lchk) goto 99
    msg2 = "(!) 2. Not found (install path): " // trim(critic_home) // cifstr

    ! then the current directory
    critic_home = "."
    inquire(file=trim(critic_home) // cifstr,exist=lchk)
    if (lchk) goto 99
    msg3 = "(!) 3. Not found (pwd): " // trim(critic_home) // cifstr

    ! argh!
    call ferror("grda_init","Could not find data files.",warning)
    if (len_trim(msgr1) > 0 .and. len_trim(msgr2) > 0) &
       write (uout,'(A/A)') msgr1, msgr2
    write (uout,'(A/A/A)') msg1, msg2, msg3
    write (uout,'("(!) The cif dict file, the density files, and the structure library")')
    write (uout,'("(!) will not be available.")')
    critic_home = "."

99  continue

    ! library file path
    clib_file = trim(critic_home) // dirsep // "lib" // dirsep // "crystal.dat"
    mlib_file = trim(critic_home) // dirsep // "lib" // dirsep // "molecule.dat"

    ! set all default values
    call global_set_defaults()

  end subroutine global_init

  !> Set the default values for all the global variables
  module subroutine global_set_defaults()
    use meshmod, only: mesh_type_franchini, mesh_level_good

    ! global flags
    precisecube = .true.

    ! bond factor
    bondfactor = 1.4d0

    ! symmetry
    doguess = -1
    symprec = 1d-2

    ! units
    iunit = iunit_bohr
    iunit_isdef = .true.

    ! navigation
    NAV_stepper = NAV_stepper_bs
    NAV_step = 0.3d0
    NAV_maxerr = 1d-4
    NAV_gradeps = 1d-7
    prunedist = 0.1d0

    ! integration
    INT_radquad_type = INT_gauleg
    INT_radquad_nr = 50
    INT_radquad_abserr = 1d0
    INT_radquad_relerr = 1d-12
    INT_radquad_errprop = 3
    INT_radquad_errprop_default = .true.
    INT_iasprec = 1d-5

    ! molecular mesh
    mesh_type = mesh_type_franchini
    mesh_level = mesh_level_good

  end subroutine global_set_defaults

  !> Print out the initial critic2 banner
  module subroutine initial_banner()
    use tools_io, only: uout

    write (uout,'("                  _   _     _          ___     ")')
    write (uout,'("                 (_) | |   (_)        |__ \    ")')
    write (uout,'("     ___   _ __   _  | |_   _    ___     ) |   ")')
    write (uout,'("    / __| | ''__| | | | __| | |  / __|   / /    ")')
    write (uout,'("   | (__  | |    | | | |_  | | | (__   / /_    ")')
    write (uout,'("    \___| |_|    |_|  \__| |_|  \___| |____|   ")')
    write (uout,*)
    write (uout,'("* CRITIC2: analysis of real-space scalar fields in solids and molecules.")')
    write (uout,'("  (c) 1996-2019 A. Otero-de-la-Roza, A. Martin-Pendas, V. Lua~na")')
    write (uout,'("  Distributed under GNU GPL v.3 (see COPYING for details)")')
    write (uout,'("  Bugs, requests, and rants: aoterodelaroza@gmail.com ")')
    write (uout,'("  Website: https://aoterodelaroza.github.io/critic2/")')
    write (uout,'("  If you find this software useful, please cite:")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 185 (2014) 1007-1018.")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 180 (2009) 157-166.")')
    write (uout,*)

  end subroutine initial_banner

  !> Print out the help message
  module subroutine help_me()
    use tools_io, only: uout

    write (uout,'("Critic2 requires a sequence of keywords to operate. Please read the manual")')
    write (uout,'("distributed with the program (https://aoterodelaroza.github.io/critic2/).")')
    write (uout,'("The command-line syntax is:")')
    write (uout,'("")')
    write (uout,'("     critic2 [-q] [-h] [-r /path/to/critic2] [-l] [inputfile [outputfile]] ")')
    write (uout,'("")')
    write (uout,'("If the input or the output file are not present, stdin and stdout are used ")')
    write (uout,'("instead. The additonal options are:")')
    write (uout,'("")')
    write (uout,'("   -h: print this message.")')
    write (uout,'("   -r: have critic2 find its data in /path/to/critic2 or /path/to/critic2/dat ")')
    write (uout,'("   -l: critic2 always generates files in the current working directory")')
    write (uout,'("   -q: quiet mode.")')
    write (uout,'("")')

  end subroutine help_me

  !> Print out the compilation details and the hardwired paths
  module subroutine config_write()
    use spglib, only: spg_get_major_version, spg_get_minor_version, spg_get_micro_version
#ifdef HAVE_LIBXC
    use xc_f90_lib_m
#endif
    use config, only: getstring, istring_package, istring_atarget,&
       istring_adate, istring_datadir, istring_version
    use param, only: dirsep
    use tools_io, only: uout, string

    logical :: lchk
    integer :: iver(3)

    write (uout,'("+ ",A," (development), version ",A,"")') getstring(istring_package), getstring(istring_version)
    write (uout,'("         host: ",A)') getstring(istring_atarget)
    write (uout,'("         date: ",A)') getstring(istring_adate)
    write (uout,'(" compiled dat: ",A)') getstring(istring_datadir)
    write (uout,'("      datadir: ",A)') trim(critic_home)
    inquire(file=trim(critic_home) // dirsep // "cif" // dirsep // "cif_core.dic",exist=lchk)
    write (uout,'("     dic file: ",A)') trim(critic_home) // dirsep // "cif" // dirsep // "cif_core.dic"
    write (uout,'("...was found?: ",L)') lchk

    iver(1) = spg_get_major_version()
    iver(2) = spg_get_minor_version()
    iver(3) = spg_get_micro_version()
    write (uout,'("       spglib: ",A,".",A,".",A)') string(iver(1)), string(iver(2)), string(iver(3))

#ifdef HAVE_LIBXC
    call xc_f90_version(iver(1),iver(2),iver(3))
    write (uout,'("        libxc: ",A,".",A,".",A)') string(iver(1)), string(iver(2)), string(iver(3))
#else
    write (uout,'("        libxc: <unavailable>")')
#endif
    write (uout,*)

  end subroutine config_write

  !> Parse the command line and set a global variable
  module subroutine critic_setvariables(line,lp)
    use meshmod, only: mesh_type_becke, mesh_type_franchini, mesh_level_small,&
       mesh_level_normal, mesh_level_good, mesh_level_vgood, mesh_level_amazing
    use arithmetic, only: eval, setvariable
    use tools_io, only: lgetword, getword, equal, isinteger, isreal, ferror, &
       faterr, string, uout, isassignment, getword, zatguess
    use param, only: maxzat0, atmcov, atmvdw
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    character(len=:), allocatable :: word, var
    logical :: ok
    real*8 :: rdum
    integer :: idum, lp2, iz
    logical :: iok, iscov

    word = lgetword(line,lp)
    if (equal(word,'bondfactor')) then
       ok = isreal(bondfactor,line,lp)
       if (.not.ok) &
          call ferror('critic_setvariables','Wrong bondfactor',faterr,line,syntax=.true.)
       bondfactor = min(bondfactor,2.0d0)
       call check_no_extra_word(ok)
    else if (equal(word,'ode_mode')) then
       do while (.true.)
          word = lgetword(line,lp)
          if (equal(word,'method')) then
             word = lgetword(line,lp)
             if (equal(word,'euler')) then
                NAV_stepper = NAV_stepper_euler
             elseif (equal(word,'heun')) then
                NAV_stepper = NAV_stepper_heun
             else if (equal(word,'bs')) then
                NAV_stepper = NAV_stepper_bs
             else if (equal(word,'rkck')) then
                NAV_stepper = NAV_stepper_rkck
             else if (equal(word,'dp')) then
                NAV_stepper = NAV_stepper_dp
             else
                call ferror('critic_setvariables','Wrong ODE_MODE stepper',faterr,line,syntax=.true.)
             end if
          else if (equal(word,'maxstep')) then
             ok = isreal(NAV_step,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/MAXSTEP',faterr,line,syntax=.true.)
             NAV_step = NAV_step / dunit0(iunit)
          else if (equal(word,'maxerr')) then
             ok = isreal(NAV_maxerr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/MAXERR',faterr,line,syntax=.true.)
          else if (equal(word,'gradeps')) then
             ok = isreal(NAV_gradeps,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/GRADEPS',faterr,line,syntax=.true.)
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in ODE_MODE',faterr,line,syntax=.true.)
          else
             exit
          end if
       end do
    else if (equal(word,'prune_distance')) then
       ok = isreal(prunedist,line,lp)
       if (.not.ok) call ferror('critic_setvariables','Wrong PRUNE_DISTANCE',faterr,line,syntax=.true.)
       prunedist = prunedist / dunit0(iunit)
    else if (equal (word,'int_radial')) then
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'type')) then
             word = lgetword(line,lp)
             if (equal(word,'gauleg')) then
                INT_radquad_type = INT_gauleg
             else if (equal(word,'qags')) then
                INT_radquad_type = INT_qags
             else if (equal(word,'qng')) then
                INT_radquad_type = INT_qng
             else if (equal(word,'qag')) then
                INT_radquad_type = INT_qag
             else
                call ferror('critic_setvariables','Wrong INT_RADIAL type',faterr,line,syntax=.true.)
             end if
          elseif (equal(word,'nr')) then
             ok = isinteger(INT_radquad_nr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL nr',faterr,line,syntax=.true.)
          elseif (equal(word,'abserr')) then
             ok = isreal(INT_radquad_abserr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL abserr',faterr,line,syntax=.true.)
          elseif (equal(word,'relerr')) then
             ok = isreal(INT_radquad_relerr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL relerr',faterr,line,syntax=.true.)
          elseif (equal(word,'errprop')) then
             ok = isinteger(INT_radquad_errprop,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL errprop',faterr,line,syntax=.true.)
             INT_radquad_errprop_default = .false.
          elseif (equal(word,'prec')) then
             ok = isreal(INT_iasprec,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL prec',faterr,line,syntax=.true.)
             INT_iasprec = INT_iasprec / dunit0(iunit)
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in INT_RADIAL',faterr,line,syntax=.true.)
          else
             exit
          end if
       end do
    else if (equal (word,'meshtype')) then
       ! meshtype {becke|franchini} [small|normal|good|verygood|amazing]
       word = lgetword(line,lp)
       if (equal(word,'becke')) then
          mesh_type = mesh_type_becke
       elseif (equal(word,'franchini')) then
          mesh_type = mesh_type_franchini
       else
          call ferror('critic_setvariables','Unknown keyword in MESHTYPE',faterr,line,syntax=.true.)
       end if

       lp2 = lp
       word = lgetword(line,lp)
       if (equal(word,'small')) then
          MESH_level = mesh_level_small
       else if (equal(word,'normal')) then
          MESH_level = mesh_level_normal
       else if (equal(word,'good')) then
          MESH_level = mesh_level_good
       else if (equal(word,'verygood')) then
          MESH_level = mesh_level_vgood
       else if (equal(word,'amazing')) then
          MESH_level = mesh_level_amazing
       else
          lp = lp2
       end if
       call check_no_extra_word(ok)
    elseif (equal(word,'library')) then
       word = lgetword(line,lp)
       if (equal(word,'crystal')) then
          clib_file = string(line(lp:))
       elseif (equal(word,'molecule')) then
          mlib_file = string(line(lp:))
       else
          call ferror('critic_setvariables','Unknown keyword in LIBRARY',faterr,line,syntax=.true.)
       end if
    elseif (equal(word,'units')) then
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'bohr').or.equal(word,'au').or.equal(word,'a.u.')) then
             iunit = iunit_bohr
             iunit_isdef = .false.
          elseif (equal(word,'ang').or.equal(word,'angstrom')) then
             iunit = iunit_ang
             iunit_isdef = .false.
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in UNITS',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do
    elseif (equal(word,'standardcube')) then
       precisecube = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'precisecube')) then
       precisecube = .true.
       call check_no_extra_word(ok)
    elseif (equal(word,'radii')) then
       word = lgetword(line,lp)
       if (equal(word,'cov')) then
          iscov = .true.
       elseif (equal(word,'vdw')) then
          iscov = .false.
       else
          call list_radii()
          return
       end if

       do while(.true.)
          lp2 = lp
          iz = 0
          ok = isinteger(iz,line,lp)
          if (.not.ok) then
             word = lgetword(line,lp)
             if (len_trim(word) == 0) then
                lp = lp2
                exit
             else
                iz = zatguess(word)
             end if
          end if

          if (iz < 1 .or. iz > maxzat0) then
             call ferror('critic2','Syntax error or wrong expression in RADII',faterr,line,syntax=.true.)
             return
          end if
          if (eval_next_real(rdum,line,lp)) then
             rdum = rdum / dunit0(iunit)
          else
             call ferror('critic2','Syntax error or wrong expression in RADII',faterr,line,syntax=.true.)
             return
          end if

          if (iscov) then
             atmcov(iz) = rdum
          else
             atmvdw(iz) = rdum
          end if
       end do
    elseif (isassignment(var,word,line)) then
       rdum = eval(word,.false.,iok)
       if (.not.iok) then
          call ferror('critic2','Syntax error or wrong assignment',faterr,line,syntax=.true.)
          return
       end if
       call setvariable(trim(var),rdum)
    else
       word = string(line)
       rdum = eval(word,.false.,iok)
       if (.not.iok) then
          call ferror('critic2','Syntax error or wrong expression',faterr,line,syntax=.true.)
          return
       end if
       if (.not.quiet) then
          write (uout,'(1p,G22.14/)') rdum
       else
          write (uout,'(1p,G22.14)') rdum
       endif
    end if
    
  contains
    subroutine check_no_extra_word(ok)
      character(len=:), allocatable :: aux2
      logical :: ok
      aux2 = getword(line,lp)
      ok = .true.
      if (len_trim(aux2) > 0) then
         call ferror('critic_setvariables','Unknown extra keyword',faterr,line,syntax=.true.)
         ok = .false.
      end if
    end subroutine check_no_extra_word
  end subroutine critic_setvariables

  !> Clear the value of a variable
  module subroutine critic_clearvariable(line)
    use arithmetic, only: clearallvariables, clearvariable
    use tools_io, only: getword, lower, equal
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word
    integer :: lp

    lp = 1
    do while(.true.)
       word = getword(line,lp)
       if (equal(lower(word),"all")) then
          call clearallvariables()
       elseif (len_trim(word) > 0) then
          call clearvariable(word)
       else
          exit
       endif
    end do

  end subroutine critic_clearvariable

  !> Evaluate next expression of word in line and return a real number.
  module function eval_next_real(res,line,lp0)
    use arithmetic, only: eval
    use tools_io, only: isexpression_or_word, string

    logical :: eval_next_real
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    real*8, intent(out) :: res
    
    integer :: lp
    character(len=:), allocatable :: word 

    res = 0d0
    lp = lp0
    word = ""
    eval_next_real = isexpression_or_word(word,line,lp)
    if (eval_next_real) then
       res = eval(string(word),.false.,eval_next_real)
       if (eval_next_real) lp0 = lp
    endif
          
  end function eval_next_real

  !> Evaluate next expression of word in line and return an integer.
  module function eval_next_int(res,line,lp0)
    use arithmetic
    use tools_io

    logical :: eval_next_int
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    integer, intent(out) :: res
    
    character(len=:), allocatable :: word 
    integer :: lp
    real*8 :: rdum

    real*8, parameter :: eps = 1d-20

    res = 0
    lp = lp0
    word = ""
    eval_next_int = isexpression_or_word(word,line,lp)
    if (eval_next_int) then
       rdum = eval(word,.false.,eval_next_int)
       if (abs(rdum - nint(rdum)) > eps) then
          eval_next_int = .false.
          return
       endif
       if (eval_next_int) then
          lp0 = lp
          res = nint(rdum)
       endif
    endif
          
  end function eval_next_int

  !> Write to standard output the list of atomic radii
  module subroutine list_radii()
    use tools_io, only: uout, string, nameguess
    use global, only: dunit0, iunit, iunitname0
    use param, only: atmcov, atmvdw
    integer :: i
    character*(2) :: name

    write (uout,'("* List of atomic radii (per atomic number)")')
    write (uout,'("# All radii in ",A)') iunitname0(iunit)
    write (uout,'("# Z at  rcov  rvdw")')
    do i = 1, maxzat0
       name = nameguess(i,.true.)
       write (uout,'(999(X,A))') string(i,length=3), name, &
          string(atmcov(i) * dunit0(iunit),'f',length=5,decimal=2), &
          string(atmvdw(i) * dunit0(iunit),'f',length=5,decimal=2)
    end do
    write (uout,*)

  end subroutine list_radii

end submodule proc
