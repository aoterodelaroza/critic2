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

! Global variables, parameters and basic initialization subroutines.
! Contains global variables and parameters, and initialization
! procedures run at the beginning of the execution.
module global
  use param, only: bohrtoa, maxzat0
  implicit none

  public

  interface eval_next
     module procedure eval_next_real
     module procedure eval_next_int
  end interface eval_next
  private :: eval_next_real, eval_next_int

  ! *** dimension parameters ***
  integer, parameter :: mneq = 200 !< maximum number of non-equivalent atoms

  ! Distance cutoffs
  real*8, parameter :: atomeps = 1d-2 !< Minimum distance to consider two atoms different
  real*8, parameter :: atomeps2 = atomeps*atomeps 

  ! Environment variables
  character(len=:), allocatable :: critic_home !< CRITIC_HOME directory (with trailing "/")
  character(len=:), allocatable :: clib_file !< The path to the crystal library file
  character(len=:), allocatable :: mlib_file !< The path to the molecular library file

  ! Cutoff radius for 1d-12 atomic densities (max of r_LDA, r_PBE). Up
  ! to Pu (93).  VAlues are bohr. This distance is usually enough to
  ! converge any density or orbital (dftb+, pi) contribution from an
  ! atom.
  real*8, parameter :: cutrad(maxzat0) = (/&
     2.149886192475d+01, 1.169139170668d+01, 3.430831385801d+01,& !001-003 H,  He, Li
     2.502075396007d+01, 2.801001395722d+01, 2.167675592180d+01,& !004-006 Be, B,  C 
     1.749805708313d+01, 1.465173060207d+01, 1.263885024136d+01,& !007-009 N,  O,  F 
     1.110521599057d+01, 3.523728402162d+01, 2.763528367271d+01,& !010-012 Ne, Na, Mg
     3.395549316507d+01, 2.847261278601d+01, 2.487715217494d+01,& !013-015 Al, Si, P 
     2.222930087269d+01, 2.022231415676d+01, 1.857150175607d+01,& !016,018 S,  Cl, Ar
     3.884428729523d+01, 3.144587224767d+01, 2.970981796151d+01,& !019-021 K,  Ca, Sc
     2.864438811442d+01, 2.784088946336d+01, 2.925194799711d+01,& !022-024 Ti, V,  Cr
     2.660566532177d+01, 2.609916866690d+01, 2.558901439004d+01,& !025-027 Mn, Fe, Co
     2.517359152887d+01, 2.691554955610d+01, 2.435659411320d+01,& !028-030 Ni, Cu, Zn
     3.467478603212d+01, 2.914443825602d+01, 2.572575006996d+01,& !031-033 Ga, Ge, As
     2.323452863278d+01, 2.134146595122d+01, 1.981582897591d+01,& !034-036 Se, Br, Kr
     3.976877622180d+01, 3.266858263171d+01, 3.027851405458d+01,& !037-039 Rb, Sr, Y 
     2.899491720657d+01, 2.967865003580d+01, 2.914637014504d+01,& !040-042 Zr, Nb, Mo
     2.697201600611d+01, 2.844039136970d+01, 2.814409350112d+01,& !043-045 Tc, Ru, Rh
     1.659926809140d+01, 2.771163603049d+01, 2.519886588880d+01,& !046-048 Pd, Ag, Cd
     3.538116802480d+01, 3.026454800251d+01, 2.695514982633d+01,& !049-051 In, Sn, Sb
     2.460024202780d+01, 2.277601677390d+01, 2.130658017554d+01,& !052-054 Te, I,  Xe
     4.137458546886d+01, 3.442036204804d+01, 3.242561614450d+01,& !055-057 Cs, Ba, La
     3.212250868201d+01, 3.306457792690d+01, 3.284026775197d+01,& !058-060 Ce, Pr, Nd
     3.262654222620d+01, 3.242292112974d+01, 3.222895504883d+01,& !061-063 Pm, Sm, Eu
     3.128431696800d+01, 3.180465714311d+01, 3.157435555252d+01,& !064-066 Gd, Tb, Dy
     3.135291924508d+01, 3.120231512615d+01, 3.099709515455d+01,& !067-069 Ho, Er, Tm
     3.079969409503d+01, 3.160515129459d+01, 2.709458469010d+01,& !070-072 Yb, Lu, Hf
     2.614193052742d+01, 2.548104664032d+01, 2.489113924347d+01,& !073-075 Ta, W,  Re
     2.441668377017d+01, 2.405143298004d+01, 2.466268008529d+01,& !076-078 Os, Ir, Pt
     2.439924398342d+01, 2.305709117567d+01, 3.643576493190d+01,& !079-081 Au, Hg, Tl
     3.110226831614d+01, 2.780342993946d+01, 2.541102668192d+01,& !082-084 Pb, Bi, Po
     2.360240806573d+01, 2.210165966674d+01, 4.053200388132d+01,& !085-087 At, Rn, Fr
     3.407838067822d+01, 3.585071927373d+01, 3.175945034367d+01,& !088-090 Ra, Ac, Th
     3.478340806986d+01, 3.489038964505d+01, 3.514212336660d+01,& !091-093 Pa, U,  Np
     3.120895952111d+01, 37d0, 37d0,&                             !094-096 Pu, Am, Cm
     37d0, 37d0, 37d0,&                                           !097-099 Bk, Cf, Es
     37d0, 37d0, 37d0,&                                           !100-102 Fm, Md, No
     37d0, 37d0, 37d0,&                                           !103-105 Lr, Rf, Db
     37d0, 37d0, 37d0,&                                           !106-108 Sg, Bh, Hs
     37d0, 37d0, 37d0,&                                           !109-111 Mt
     37d0, 37d0, 37d0,&                                           !112-114
     37d0, 37d0, 37d0,&                                           !115-117
     37d0,&                                                       !118
     0d0, 0d0, 0d0, 0d0, 0d0&                                     !119-123
     /) !< Cutoff radius for 1d-12 atomic densities (max of r_LDA, r_PBE)

  ! global control for critic
  character(len=:), allocatable :: fileroot !< file prefix
  logical :: quiet
  logical :: precisecube

  ! units
  integer :: iunit
  logical :: iunit_isdef
  character(len=:), allocatable :: iunitname
  real*8 :: dunit

  integer, parameter :: iunit_bohr = 1
  integer, parameter :: iunit_ang = 2
  character*4, parameter :: iunitname0(2) = (/"bohr","ang_"/)
  real*8, parameter :: dunit0(2) = (/1d0,bohrtoa/)

  ! default border for a molecular unit cell
  real*8, parameter :: rborder_def = 10d0 / bohrtoa

  ! guess and symmetry option (-1 = only for small systems, 0 = no, 1 = cen, 2 = full)
  integer :: doguess

  ! A crystal is considered small if it has less than this number of
  ! atoms in the unit cell.
  integer, parameter :: crsmall = 500

  ! reference scalar field
  integer :: refden

  ! navigation options
  real*8 :: NAV_step !< gradient path step length (bohr)
  real*8 :: NAV_gradeps !< gradient path gradient mod termination threshold
  real*8 :: NAV_maxerr !< maximum error in gradient path tracing (bohr)
  integer :: NAV_stepper  !< the stepper in gp tracing
  integer, parameter :: NAV_stepper_euler = 1 !< Euler stepper (1 eval), poor-man's adaptive step
  integer, parameter :: NAV_stepper_rkck  = 2 !< Runge-Kutta-Cash-Karp embedded 4(5)-order, local extrapolation (6 eval), with error estimate
  integer, parameter :: NAV_stepper_dp    = 3 !< Dormand-prince embedded 4(5)-order, local extrapolation (7 eval), with error estimate
  integer, parameter :: NAV_stepper_bs    = 4 !< Bogacki-Shampine embedded 2(3)-order method, (5-1=4 eval, fsal), with error estimate
  integer, parameter :: NAV_stepper_heun  = 5 !< Heun stepper (2 eval), poor-man's adaptive step
  real*8 :: prunedist

  ! critical points
  real*8 :: CP_hdegen = 1d-8 !< a CP is degenerate if any Hessian element is less than this value

  ! radial integration
  integer :: INT_radquad_type !< type of radial integration
  integer :: INT_radquad_nr !< number of nodes in radial integration
  real*8 :: INT_radquad_abserr !< req'd absolute error in radial integration
  real*8 :: INT_radquad_relerr !< req'd relative error in radial integration
  integer :: INT_radquad_errprop !< adaptive integration driver property
  logical :: INT_radquad_errprop_default !< has the property been set by the user?
  real*8 :: INT_iasprec !< precision of the IAS for radial integration
  integer, parameter :: INT_gauleg = 1 !< gauss legendre
  integer, parameter :: INT_qags = 2 !< quadpack qags
  integer, parameter :: INT_qng = 3 !< quadpack qng
  integer, parameter :: INT_qag = 4 !< quadpack qag
  integer, parameter :: INT_lebedev = 5 !< lebedev

  ! mesh type and quality for molecular integrations
  ! 0 = Becke, >0 = Franchini with 1 (small), 2 (normal), 3 (good),
  ! 4(very good), 5 (excellent).
  integer :: MESH_type !< type and quality of mesh for molecular integration

  ! qtree
  ! note: ws_origin also used in auto
  integer :: gradient_mode   
  integer :: qtree_ode_mode    
  real*8 :: stepsize
  real*8 :: ode_abserr
  integer :: mpstep 
  real*8  :: qtreefac
  integer :: integ_mode(20)
  integer :: integ_scheme
  integer :: keastnum
  integer :: prop_mode 
  real*8  :: sphfactor(mneq), sphintfactor(mneq)
  real*8  :: cub_abs, cub_rel
  integer :: cub_mpts
  integer :: plot_mode
  logical :: docontacts
  real*8 :: ws_origin(3), ws_scale
  logical :: ws_use
  integer :: autosph
  logical :: killext, checkbeta
  integer :: minl 
  logical :: plotsticks
  integer :: qtree_presplit
  integer :: color_allocate
  integer :: setsph_lvl
  real*8  :: vcutoff

contains

  !> Initialize basic variables at the beginning of the run.
  !> Also sets default values by calling global_set_defaults.
  subroutine global_init(ghome,datadir)
    use tools_io, only: string, ferror, warning, uout
    use param, only: dirsep

    character*(*) :: ghome, datadir

    integer :: isenv
    logical :: lchk

    character(len=:), allocatable :: msg1, msg2, msg3, msg4
    integer, parameter :: maxlenpath = 1024

    ! read the -r option
    if (len_trim(ghome) > 0) then
       critic_home = string(ghome) // dirsep // "dat"
       inquire(file=trim(critic_home)// dirsep // "cif_core.dic",exist=lchk)
       if (lchk) goto 99
       write (uout,'("  File not found: ",A)') trim(critic_home)// dirsep // "cif_core.dic"
    endif

    ! read env variable CRITIC_HOME
    allocate(character(len=maxlenpath)::critic_home)
    call get_environment_variable("CRITIC_HOME",critic_home,status=isenv)
    if (isenv ==0) then
       critic_home = trim(critic_home) // dirsep // "dat"
       inquire(file=trim(critic_home)// dirsep // "cif_core.dic",exist=lchk)
       if (lchk) goto 99
       msg1 = "(!) 1. Not found (CRITIC_HOME): " // trim(critic_home)// dirsep // "cif_core.dic"
    else
       msg1 = "(!) 1. CRITIC_HOME environment variable not set"
    end if

    ! then the install path
    critic_home = trim(adjustl(datadir))
    inquire(file=trim(critic_home)// dirsep // "cif_core.dic",exist=lchk)
    if (lchk) goto 99
    msg2 = "(!) 2. Not found (install path): " // trim(critic_home)// dirsep // "cif_core.dic"

    ! then the current directory
    critic_home = "."
    inquire(file=trim(critic_home)// dirsep // "cif_core.dic",exist=lchk)
    if (lchk) goto 99
    msg3 = "(!) 3. Not found (pwd): " // trim(critic_home)// dirsep // "cif_core.dic"

    ! argh!
    call ferror("grda_init","Could not find data files.",warning)
    write (uout,'(A/A/A)') msg1, msg2, msg3
    write (uout,'("(!) The cif dict file, the density files, and the structure library")')
    write (uout,'("(!) will not be available.")')
    critic_home = "."

99  continue

    ! library file path
    clib_file = trim(critic_home) // dirsep // "crystal.dat"
    mlib_file = trim(critic_home) // dirsep // "molecule.dat"

    ! set all default values
    call global_set_defaults()

  end subroutine global_init

  !> Set the default values for all the global variables
  subroutine global_set_defaults()

    doguess = -1
    refden = 0
    precisecube = .true.

    ! units
    iunit = iunit_bohr
    iunitname = trim(iunitname0(iunit))
    dunit = dunit0(iunit)
    iunit_isdef = .true.

    ! navigation
    NAV_stepper = NAV_stepper_bs
    NAV_step = 0.1d0
    NAV_maxerr = 1d-3
    NAV_gradeps = 1d-9
    prunedist = 0.1d0

    ! integration
    INT_radquad_type = INT_gauleg
    INT_radquad_nr = 50
    INT_radquad_abserr = 1d0
    INT_radquad_relerr = 1d-12
    INT_radquad_errprop = 2
    INT_radquad_errprop_default = .true.
    INT_iasprec = 1d-5

    ! molecular mesh
    MESH_type = 2

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
    ws_use = .true.
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
  subroutine initial_banner()
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
  subroutine help_me()
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
  subroutine config_write(package,version,atarget,adate,f77,fflags,fc,&
     fcflags,ldflags,enable_debug,datadir)
    use param, only: dirsep
    use tools_io, only: uout
    character*(*), intent(in) :: package, version, atarget, adate
    character*(*), intent(in) :: f77, fflags, fc, fcflags
    character*(*), intent(in) :: ldflags, enable_debug, datadir
    logical :: lchk

    write (uout,'("+ ",A,", commit ",A,"")') package, version
    write (uout,'(" compile host: ",A)') atarget
    write (uout,'(" compile date: ",A)') adate
    write (uout,'("    using f77: ",A," ",A)') f77, fflags
    write (uout,'("          f90: ",A," ",A)') fc, fcflags
    write (uout,'("      ldflags: ",A)') ldflags
    write (uout,'("       debug?: ",A)') enable_debug
    write (uout,'(" compiled dat: ",A)') datadir
    write (uout,'("      datadir: ",A)') trim(critic_home)
    inquire(file=trim(critic_home)// dirsep // "cif_core.dic",exist=lchk)
    write (uout,'("     dic file: ",A)') trim(critic_home)// dirsep // "cif_core.dic"
    write (uout,'("...was found?: ",L)') lchk
    write (uout,*)

  end subroutine config_write

  !> Parse the command line and set a global variable
  subroutine critic_setvariables(line,lp)
    use arithmetic, only: eval, setvariable
    use tools_io, only: lgetword, getword, equal, isinteger, isreal, ferror, &
       faterr, string, uout, isassignment, getword
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    character(len=:), allocatable :: word, var
    logical :: ok
    real*8 :: rdum
    integer :: idum
    logical :: iok

    word = lgetword(line,lp)
    if (equal(word,'nosymm') .or. equal(word,'nosym')) then
       doguess = 0
       call check_no_extra_word(ok)
    elseif (equal(word,'symm')) then
       ok = isinteger(doguess,line,lp)
       if (.not.ok) doguess = 2
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
    elseif (equal(word,'nows')) then
       ws_use = .false.
       call check_no_extra_word(ok)
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
       dunit = dunit0(iunit)
    elseif (equal(word,'standardcube')) then
       precisecube = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'precisecube')) then
       precisecube = .true.
       call check_no_extra_word(ok)
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
  subroutine critic_clearvariable(line)
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
  function eval_next_real(res,line,lp0)
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
  function eval_next_int(res,line,lp0)
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

end module global
