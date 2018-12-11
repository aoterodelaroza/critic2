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

    character(len=:), allocatable :: cifstr, msg1, msg2, msg3
    integer, parameter :: maxlenpath = 1024

    cifstr = dirsep // "cif" // dirsep // "cif_core.dic"

    ! read the -r option
    if (len_trim(ghome) > 0) then
       critic_home = string(ghome) // dirsep // "dat"
       inquire(file=trim(critic_home) // cifstr,exist=lchk)
       if (lchk) goto 99
       write (uout,'("  File not found: ",A)') trim(critic_home) // cifstr
    endif

    ! read env variable CRITIC_HOME
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

    precisecube = .true.

    ! bond factor
    bondfactor = 1.4d0

    ! symmetry
    doguess = -1
    symprec = 1d-4

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
    MESH_type = 3

    ! qtree
    ! gradient_mode set in critic->set_reference
    gradient_mode = 1
    qtree_ode_mode = 8
    stepsize = 0.1d0
    ode_abserr = -1d0
    mpstep = 0
    qtreefac = 2d0
    integ_mode = 0
    integ_scheme = 1
    keastnum = 5
    cub_abs = 1d0
    cub_rel = 1d-6
    cub_mpts = 1000
    prop_mode = 1
    sphfactor = 0d0
    sphintfactor = 1d0
    plot_mode = 0
    docontacts = .false.
    ws_origin = (/0d0, 0d0, 0d0/)
    ws_scale = -1d0
    killext = .true.
    plotsticks = .true.
    autosph = 2
    checkbeta = .false.
    minl = 4
    qtree_presplit = 0
    color_allocate = -1
    setsph_lvl = 6
    vcutoff = 0d0

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
    write (uout,'("  (c) 1996-2015 A. Otero-de-la-Roza, A. Martin-Pendas, V. Lua~na")')
    write (uout,'("  Distributed under GNU GPL v.3 (see COPYING for details)")')
    write (uout,'("  Bugs, requests, and rants: aoterodelaroza@gmail.com ")')
    write (uout,'("  If you find this software useful, please cite:")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 185 (2014) 1007-1018.")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 180 (2009) 157-166.")')
    write (uout,*)

  end subroutine initial_banner

  !> Print out the help message
  module subroutine help_me()
    use tools_io, only: uout

    write (uout,'("Critic2 requires a sequence of keywords to operate. Please read the manual")')
    write (uout,'("distributed with the program (dat/user-guide.txt) and the syntax file")')
    write (uout,'("(doc/syntax.txt). The command-line syntax is:")')
    write (uout,'("")')
    write (uout,'("     critic2 [-q] [-h] [-r /path/to/critic2] [inputfile [outputfile]] ")')
    write (uout,'("")')
    write (uout,'("If the input or the output file are not present, stdin and stdout are used ")')
    write (uout,'("instead. The additonal options are:")')
    write (uout,'("")')
    write (uout,'("   -r: tell critic2 that its data is in /path/to/critic2/dat ")')
    write (uout,'("   -h: print this message.")')
    write (uout,'("   -q: quiet mode.")')
    write (uout,'("")')

  end subroutine help_me

  !> Print out the compilation details and the hardwired paths
  module subroutine config_write()
    use config, only: getstring, istring_package,&
       istring_f77, istring_fc, istring_cc, istring_fflags,&
       istring_fcflags, istring_cflags, istring_ldflags,&
       istring_atarget, istring_adate, istring_enabledebug,&
       istring_datadir, istring_version
    use param, only: dirsep
    use tools_io, only: uout
    logical :: lchk

    write (uout,'("+ ",A," (development), commit ",A,"")') getstring(istring_package), getstring(istring_version)
    write (uout,'(" compile host: ",A)') getstring(istring_atarget)
    write (uout,'(" compile date: ",A)') getstring(istring_adate)
    write (uout,'("    using f77: ",A," ",A)') getstring(istring_f77), getstring(istring_fflags)
    write (uout,'("          f90: ",A," ",A)') getstring(istring_fc), getstring(istring_fcflags)
    write (uout,'("            c: ",A," ",A)') getstring(istring_cc), getstring(istring_cflags)
    write (uout,'("      ldflags: ",A)') getstring(istring_ldflags)
    write (uout,'("       debug?: ",A)') getstring(istring_enabledebug)
    write (uout,'(" compiled dat: ",A)') getstring(istring_datadir)
    write (uout,'("      datadir: ",A)') trim(critic_home)
    inquire(file=trim(critic_home) // dirsep // "cif" // dirsep // "cif_core.dic",exist=lchk)
    write (uout,'("     dic file: ",A)') trim(critic_home) // dirsep // "cif" // dirsep // "cif_core.dic"
    write (uout,'("...was found?: ",L)') lchk
    write (uout,*)

  end subroutine config_write

  !> Parse the command line and set a global variable
  module subroutine critic_setvariables(line,lp)
    use arithmetic, only: eval, setvariable
    use tools_io, only: lgetword, getword, equal, isinteger, isreal, ferror, &
       faterr, string, uout, isassignment, getword, zatguess
    use param, only: maxzat0, atmcov
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    character(len=:), allocatable :: word, var
    logical :: ok
    real*8 :: rdum
    integer :: idum, lp2
    logical :: iok

    word = lgetword(line,lp)
    if (equal(word,'nosymm') .or. equal(word,'nosym')) then
       doguess = 0
       call check_no_extra_word(ok)
    elseif (equal(word,'symm').or.equal(word,'sym')) then
       ok = isinteger(doguess,line,lp)
       if (.not.ok) doguess = 1
       call check_no_extra_word(ok)
    elseif (equal(word,'symprec')) then
       ok = isreal(symprec,line,lp)
       if (.not.ok) &
          call ferror('critic_setvariables','Wrong symprec',faterr,line,syntax=.true.)
       call check_no_extra_word(ok)
    elseif (equal(word,'bondfactor')) then
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
       word = lgetword(line,lp)
       if (equal(word,'becke')) then
          MESH_type = 0
       elseif (equal(word,'franchini')) then
          ok = isinteger(MESH_type,line,lp)
          if (MESH_type < 0 .or. MESH_type > 5) &
             call ferror('critic_setvariables','Invalid Franchini quality in MESHTYPE',faterr,line,syntax=.true.)
       else
          call ferror('critic_setvariables','Unknown keyword in MESHTYPE',faterr,line,syntax=.true.)
       end if
       call check_no_extra_word(ok)
       
    else if (equal (word,'gradient_mode')) then
       ok = isinteger(gradient_mode,line,lp)
       if (.not.ok) then
          call ferror('critic_setvariables','Wrong gradient_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    else if (equal (word,'qtree_ode_mode')) then
       ok = isinteger(qtree_ode_mode,line,lp)
       if (.not.ok) then
          call ferror('critic_setvariables','Wrong qtree_ode_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    else if (equal (word,'stepsize')) then
       ok = isreal(stepsize,line,lp)
       if (.not.ok) then
          call ferror('critic_setvariables','Wrong stepsize',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
       stepsize = stepsize / dunit0(iunit)
    elseif (equal(word,'ode_abserr')) then
       ok = isreal(ode_abserr,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong ode_abserr',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'mpstep')) then
       ok = isinteger(mpstep,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong mpstep',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'qtreefac')) then
       ok = isreal(qtreefac,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong qtreefac',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'integ_mode')) then
       ok = isinteger(idum,line,lp)
       if (idum > 0) then
          ok = ok .and. isinteger(integ_mode(idum),line,lp)
       else
          ok = ok .and. isinteger(idum,line,lp)
          integ_mode = idum
       end if
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong integ_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'integ_scheme')) then
       ok = isinteger(integ_scheme,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong integ_scheme',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'keastnum')) then
       ok = isinteger(keastnum,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong keastnum',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'prop_mode')) then
       ok = isinteger(prop_mode,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong prop_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'sphintfactor')) then
       ok = isinteger(idum,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong sphintfactor',faterr,line,syntax=.true.)
       else
          ok = isreal(rdum,line,lp)
          if (.not. ok .or. idum == 0d0) then
             sphintfactor = rdum
          else
             sphintfactor(idum) = rdum
          end if
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'cub_abs')) then
       ok = isreal(cub_abs,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong cub_abs',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'cub_rel')) then
       ok = isreal(cub_rel,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong cub_rel',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'cub_mpts')) then
       ok = isinteger(cub_mpts,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong cub_mpts',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'plot_mode')) then
       ok = isinteger(plot_mode,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong plot_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'docontacts')) then
       docontacts = .true.
       call check_no_extra_word(ok)
    elseif (equal(word,'nocontacts')) then
       docontacts = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'ws_origin')) then
       ok = isreal(ws_origin(1),line,lp)
       ok = ok .and. isreal(ws_origin(2),line,lp)
       ok = ok .and. isreal(ws_origin(3),line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong ws_origin',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'ws_scale')) then
       ok = isreal(ws_scale,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong ws_scale',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'killext')) then
       killext = .true.
       call check_no_extra_word(ok)
    elseif (equal(word,'nokillext')) then
       killext = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'autosph')) then
       ok = isinteger(autosph,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong autosph',faterr,line,syntax=.true.)
       else if (autosph /= 1 .and. autosph /= 2) then
          call ferror('critic_setvariables','autosph must be 1 or 2',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'checkbeta')) then
       checkbeta = .true.
       call check_no_extra_word(ok)
    elseif (equal(word,'nocheckbeta')) then
       checkbeta = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'qtree_minl')) then
       ok = isinteger(minl,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong gradient_mode',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'plotsticks')) then
       plotsticks = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'noplotsticks')) then
       plotsticks = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'color_allocate')) then
       ok = isinteger(color_allocate,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong color_allocate',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
    elseif (equal(word,'setsph_lvl')) then
       ok = isinteger(color_allocate,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong setsph_lvl',faterr,line,syntax=.true.)
       else if (setsph_lvl > 7 .or. setsph_lvl < minl) then
          call ferror('critic->qtree','setsph_lvl must be > minl and < 8',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       endif
    elseif (equal(word,'vcutoff')) then
       ok = isreal(vcutoff,line,lp)
       if (.not. ok) then
          call ferror('critic_setvariables','Wrong vcutoff',faterr,line,syntax=.true.)
       else
          call check_no_extra_word(ok)
       end if
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
       do while(.true.)
          lp2 = lp
          ok = isinteger(idum,line,lp)
          if (.not.ok) then
             word = lgetword(line,lp)
             if (len_trim(word) == 0) then
                lp = lp2
                exit
             else
                idum = zatguess(word)
                if (idum < 1 .or. idum > maxzat0) then
                   call ferror('critic2','Syntax error or wrong expression',faterr,line,syntax=.true.)
                   return
                end if
             end if
          end if
          if (eval_next_real(rdum,line,lp)) then
             atmcov(idum) = rdum / dunit0(iunit)
          else
             call ferror('critic2','Syntax error or wrong expression',faterr,line,syntax=.true.)
             return
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

end submodule proc
