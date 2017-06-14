! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
! Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
! Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
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

!> Interface for the GUI.
module gui_interface
  use crystalmod, only: crystal
  use c_interface_module
  implicit none

  private

  public :: critic2_initialize
  private :: critic2_initialize_library
  public :: critic2_end
  public :: open_structure
  public :: open_structure_from_library
  public :: new_structure
  public :: preview_structure
  public :: accept_previewed_structure
  public :: reject_previewed_structure
  private :: initialize_structure
  public :: call_auto
  public :: update_scene
  public :: clear_scene
  public :: get_text_info
  public :: get_seed_from_current_structure
  public :: set_library_file
  private :: stick_from_endpoints
  private :: denewline
  private :: save_state
  private :: load_state
  private :: clear_state
  public :: run_critic2_command

  ! constants
  character*1, parameter :: nl = new_line('a')

  ! C-interoperable stick type
  type, bind(c) :: c_stick
     real(c_float) :: r1(3) !< position of the first end (bohr) 
     real(c_float) :: r2(3) !< position of the second end (bohr)
     real(c_float) :: rmid(3) !< position of the bond center (bohr)
     real(c_float) :: length !< length of the bond (bohr)
     real(c_float) :: thick !< thickness (bohr)
     real(c_float) :: rgb(3) !< color
     real(c_float) :: rot(4,4) !< rotation matrix for the stick
  end type c_stick

  ! C-interoperable ball type
  type, bind(c) :: c_ball
     real(c_float) :: r(3) !< center position (bohr) 
     real(c_float) :: rad !< ball radius (bohr) 
     real(c_float) :: rgb(3) !< Color
  end type c_ball

  ! C-interoperable atom type
  type, bind(c) :: c_atom
     character(kind=c_char,len=1) :: name(11) !< Atomic name
     integer(c_int) :: z !< Atomic number
     type(c_ball) :: b !< ball representation of this atom
  end type c_atom

  ! C-interoperable bond type
  type, bind(c) :: c_bond
     integer(c_int) :: i1 !< index for the first atom (C-style)
     integer(c_int) :: i2 !< index for the second atom (C-style)
     type(c_stick) :: s !< stick representation of this bond
  end type c_bond

  ! C-interoperable critical point
  type, bind(c) :: c_critp
     integer(c_int) :: type !< type (-3,-1,1,3)
     character(kind=c_char,len=1) :: name(11) !< Name
     type(c_ball) :: b !< ball representation of this CP
  end type c_critp

  ! C-interoperable crystal seed (for the GUI input)
  type, bind(c) :: c_crystalseed
     integer(c_int) :: type ! 0 = molecule, 1 = crystal
     integer(c_int) :: achoice ! 0 = aa, 1 = rr
     character(kind=c_char,len=1) :: straa(255)
     character(kind=c_char,len=1) :: strbb(255)
     character(kind=c_char,len=1) :: strrr(255*5)
     character(kind=c_char,len=1) :: strat(255*1024)
     character(kind=c_char,len=1) :: strspg(255)
     logical(c_bool) :: molcubic
     real(c_float) :: molborder
     integer(c_int) :: borunits ! 0 = bohr, 1 = angstrom
     integer(c_int) :: aaunits ! 0 = bohr, 1 = angstrom
     integer(c_int) :: rrunits ! 0 = bohr, 1 = angstrom
     integer(c_int) :: atunits ! 0 = bohr, 1 = angstrom, 2 = fractional
     integer(c_int) :: errcode
     character(kind=c_char,len=1) :: errmsg(255)
  end type c_crystalseed

  !! Interface to the GUI code !!
  ! global flags
  logical(c_bool), bind(c) :: isinit
  logical(c_bool), bind(c) :: ispreview
  logical(c_bool), bind(c) :: ismolecule

  ! number of atoms in the scene
  integer(c_int), bind(c) :: nat

  ! allocatable pointer for atom information (fortran side)
  type(c_atom), pointer :: at_f(:)

  ! C pointer to atom information
  type(c_ptr), bind(c) :: at 
  
  ! number of bonds in the scene
  integer(c_int), bind(c) :: nbond

  ! allocatable pointer for bond information (fortran side)
  type(c_bond), pointer :: bond_f(:)

  ! C pointer to bond information
  type(c_ptr), bind(c) :: bond
  
  ! number of atoms in the scene
  integer(c_int), bind(c) :: ncritp

  ! allocatable pointer for atom information (fortran side)
  type(c_critp), pointer :: critp_f(:)

  ! C pointer to atom information
  type(c_ptr), bind(c) :: critp
  
  ! unit cell and lattice vectors
  real(c_float), bind(c) :: cell_x0(3)
  real(c_float), bind(c) :: cell_lat(3,3)
  integer(c_int), bind(c) :: cell_nstick
  type(c_stick), bind(c) :: cell_s(12)

  ! bounding box limits and center
  real(c_float), bind(c) :: box_xmin(3)
  real(c_float), bind(c) :: box_xmax(3)
  real(c_float), bind(c) :: box_xcm(3)
  real(c_float), bind(c) :: box_xmaxlen
  real(c_float), bind(c) :: box_xmaxclen

  ! structures in the library
  integer(c_int), bind(c) :: nlib_crys
  integer(c_int), bind(c) :: nlib_mol
  character(kind=c_char,len=1), pointer :: lib_crys_f(:,:)
  type(c_ptr), bind(c) :: lib_crys
  character(kind=c_char,len=1), pointer :: lib_mol_f(:,:)
  type(c_ptr), bind(c) :: lib_mol

  !! Save state (for GUI previews) !!
  logical :: issaved = .false.
  logical(c_bool) :: isinit_
  logical(c_bool) :: ismolecule_
  integer(c_int) :: nat_
  type(c_atom), pointer :: at_f_(:)
  integer(c_int) :: nbond_
  type(c_bond), pointer :: bond_f_(:)
  integer(c_int) :: ncritp_
  type(c_critp), pointer :: critp_f_(:)
  real(c_float) :: cell_x0_(3)
  real(c_float) :: cell_lat_(3,3)
  integer(c_int) :: cell_nstick_
  type(c_stick) :: cell_s_(12)
  real(c_float) :: box_xmin_(3)
  real(c_float) :: box_xmax_(3)
  real(c_float) :: box_xcm_(3)
  real(c_float) :: box_xmaxlen_
  real(c_float) :: box_xmaxclen_
  type(crystal) :: cr_

contains

  !> Initialize critic2.
  subroutine critic2_initialize() bind(c)
    use graphics, only: graphics_init
    use spgs, only: spgs_init
    use fields, only: fields_init, fields_end
    use crystalmod, only: cr
    use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
       fcflags, cc, cflags, ldflags, enable_debug, package
    use global, only: global_init, fileroot, config_write, initial_banner
    use tools_io, only: stdargs, ioinit, ucopy, uout, start_clock, &
       tictac
    use param, only: param_init

    character(len=:), allocatable :: optv, ghome

    ! initialize parameters
    call start_clock()
    call param_init()

    ! input/output, arguments (tools_io)
    call ioinit()
    call stdargs(optv,ghome,fileroot)

    ! set default values and initialize the rest of the modules
    call global_init(ghome,datadir)
    call cr%init()
    call fields_init()
    call spgs_init()
    call graphics_init()

    ! banner and compilation info; do not copy input
    call initial_banner()
    call config_write(package,version,atarget,adate,f77,fflags,fc,&
       fcflags,cc,cflags,ldflags,enable_debug,datadir)
    call tictac('CRITIC2')
    write (uout,*)
    ucopy = -1

    ! initialize libraries
    call critic2_initialize_library()

    ! clear the scene
    call clear_scene(logical(.false.,c_bool))

  end subroutine critic2_initialize

  subroutine critic2_initialize_library()
    use global, only: clib_file, mlib_file
    use tools_io, only: fopen_read, fclose, lgetword, getline, equal,&
       getword
    character(len=:), allocatable :: word, line
    logical :: lchk
    integer :: lu, lp, n

    ! entries in the crystal library
    nlib_crys = 0
    inquire(file=clib_file,exist=lchk)
    if (lchk) then
       ! number of entries
       lu = fopen_read(clib_file,abspath0=.true.)
       do while (getline(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,'structure')) then
             word = getword(line,lp)
             nlib_crys = nlib_crys + 1
          endif
       end do
       rewind(lu)

       ! read the entries
       if (associated(lib_crys_f)) deallocate(lib_crys_f)
       allocate(lib_crys_f(255,nlib_crys))
       n = 0
       do while (getline(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,'structure')) then
             word = getword(line,lp)
             n = n + 1
             call f_c_string(word,lib_crys_f(:,n))
          endif
       end do
       lib_crys = c_loc(lib_crys_f)
       call fclose(lu)
    endif

    ! entries in the molecule library
    nlib_mol = 0
    inquire(file=mlib_file,exist=lchk)
    if (lchk) then
       ! number of entries
       lu = fopen_read(mlib_file,abspath0=.true.)
       do while (getline(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,'structure')) then
             word = getword(line,lp)
             nlib_mol = nlib_mol + 1
          endif
       end do
       rewind(lu)

       ! read the entries
       if (associated(lib_mol_f)) deallocate(lib_mol_f)
       allocate(lib_mol_f(255,nlib_mol))
       n = 0
       do while (getline(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,'structure')) then
             word = getword(line,lp)
             n = n + 1
             call f_c_string(word,lib_mol_f(:,n))
          endif
       end do
       lib_mol = c_loc(lib_mol_f)
       call fclose(lu)
    endif

  end subroutine critic2_initialize_library

  !> End of the critic2 run
  subroutine critic2_end() bind(c)
    use fields, only: fields_end
    use crystalmod, only: cr
    use pi_private, only: pi_end
    use wfn_private, only: wfn_end
    use varbas, only: varbas_end
    use grd_atomic, only: grda_end
    use tools_io, only: uout, nwarns, ncomms, print_clock, tictac, string
    
    call pi_end()
    call wfn_end()
    call grda_end()
    call varbas_end()
    call cr%end()
    call fields_end()
    
    write (uout,'("CRITIC2 ended succesfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
       string(nwarns), string(ncomms)
    call print_clock()
    call tictac('CRITIC2')
  end subroutine critic2_end

  ! Read a new molecule/crystal from an external file
  subroutine open_structure(filename0, ismolecule) bind(c)
    use struct_drivers, only: struct_crystal_input
    use crystalmod, only: cr

    type(c_ptr), intent(in) :: filename0
    integer(c_int), value :: ismolecule
    character(len=:), allocatable :: filename

    ! transform to fortran string
    filename = c_string_value(filename0)

    ! read the external file name 
    call struct_crystal_input(cr, filename, ismolecule, .true.)

    ! initialize and update the scene
    if (cr%isinit) then
       ispreview = .false.
       call initialize_structure()
       call update_scene()
    end if

  end subroutine open_structure

  ! Read a new molecule/crystal from an external file
  subroutine open_structure_from_library(nstr, ismolecule) bind(c)
    use crystalmod, only: cr
    use crystalseedmod, only: crystalseed
    use tools_io, only: ferror, faterr

    integer(c_int), value :: nstr
    integer(c_int), value :: ismolecule
    character*255 :: name
    type(crystalseed) :: seed
    logical :: oksyn

    ! transform to fortran string
    if (ismolecule == 0) then
       call c_f_string(lib_crys_f(:,nstr),name)
    else
       call c_f_string(lib_mol_f(:,nstr),name)
    end if

    ! read the external file name 
    call seed%struct_read_library(name, ismolecule==1, oksyn)
    if (.not.oksyn) &
       call ferror("open_structure_from_library","name not found",faterr)

    ! build the crystal structure
    call cr%struct_new(seed,.false.)

    ! initialize and update the scene
    if (cr%isinit) then
       ispreview = .false.
       call initialize_structure()
       call update_scene()
    end if

  end subroutine open_structure_from_library

  ! Read a new molecule/crystal from an external file
  function new_structure(useed,preview) bind(c)
    use crystalmod, only: cr
    use crystalseedmod, only: crystalseed
    use global, only: eval_next
    use tools_io, only: getword, zatguess
    use tools_math, only: matinv, det
    use types, only: realloc
    use param, only: bohrtoa
    type(c_crystalseed), intent(inout) :: useed
    logical(c_bool), value :: preview
    integer(c_int) :: new_structure

    type(crystalseed) :: fseed
    logical :: ok
    integer :: lp, i, j
    character(len=255*1024) :: aux
    real*8 :: xat(3), r(3,3), dd
    character(len=:), allocatable :: atsym

    ! set default return value
    new_structure = 0
    useed%errcode = 0
    call f_c_string("",useed%errmsg,254)

    ! build a crystal seed from the information passed from the GUI
    fseed%isused = .true.
    fseed%file = "<gui input>"
    if (useed%type == 1) then
       ! a crystal
       fseed%ismolecule = .false.
       if (useed%achoice == 0) then
          fseed%useabr = 1
       else
          fseed%useabr = 2
       end if
    elseif  (useed%type == 0) then
       ! a molecule
       fseed%ismolecule = .true.
       fseed%cubic = useed%molcubic
       fseed%border = real(useed%molborder,8)
       if (useed%borunits == 1) fseed%border = fseed%border / bohrtoa
       fseed%useabr = 0
    else
       useed%errcode = 1
       call f_c_string("You must choose between molecule and crystal.",useed%errmsg,254)
       return
    end if

    ! parse the cell information
    ok = .true.
    lp = 1
    if (fseed%useabr == 1) then
       fseed%aa = 0d0
       call c_f_string(useed%straa,aux)
       do i = 1, 3
          ok = ok .and. eval_next(fseed%aa(i),aux,lp)
       end do
       if (.not.ok .or. any(fseed%aa < 1d-10)) then
          useed%errcode = 2
          call f_c_string("Incorrect cell lengths.",useed%errmsg,254)
          return
       end if
       if (useed%aaunits == 1) fseed%aa = fseed%aa / bohrtoa

       fseed%bb = 0d0
       call c_f_string(useed%strbb,aux)
       lp = 1
       do i = 1, 3
          ok = ok .and. eval_next(fseed%bb(i),aux,lp)
       end do
       if (.not.ok .or. any(fseed%bb < 1d-10)) then
          useed%errcode = 3
          call f_c_string("Incorrect cell angles.",useed%errmsg,254)
          return
       end if
    elseif (fseed%useabr == 2) then
       call c_f_string(useed%strrr,aux)
       call denewline(aux)
       do i = 1, 3
          do j = 1, 3
             ok = ok .and. eval_next(fseed%crys2car(j,i),aux,lp)
          end do
       end do
       if (.not.ok) then
          useed%errcode = 4
          call f_c_string("Incorrect lattice vectors.",useed%errmsg,254)
          return
       end if
       if (abs(det(fseed%crys2car)) < 1d-6) then
          useed%errcode = 5
          call f_c_string("Cell has zero volume.",useed%errmsg,254)
          return
       end if
       if (useed%rrunits == 1) fseed%crys2car = fseed%crys2car / bohrtoa
    end if

    ! parse the atoms
    call c_f_string(useed%strat,aux)
    call denewline(aux)
    lp = 1
    fseed%nat = 0
    allocate(fseed%x(3,10),fseed%name(10))
    do while (.true.)
       atsym = getword(aux,lp)
       if (len_trim(atsym) == 0) exit

       ok = (zatguess(atsym) >= 0)
       if (.not.ok) then
          useed%errcode = 6
          call f_c_string("Unknown atomic symbol: " // trim(atsym) // ".",useed%errmsg,254)
          return
       end if

       ok = eval_next(xat(1),aux,lp)
       ok = ok .and. eval_next(xat(2),aux,lp)
       ok = ok .and. eval_next(xat(3),aux,lp)
       if (.not.ok) exit
       fseed%nat = fseed%nat + 1
       if (fseed%nat > size(fseed%x,2)) then
          call realloc(fseed%x,3,2*fseed%nat)
          call realloc(fseed%name,2*fseed%nat)
       end if
       fseed%x(:,fseed%nat) = xat
       fseed%name(fseed%nat) = atsym
    end do
    call realloc(fseed%x,3,fseed%nat)
    call realloc(fseed%name,fseed%nat)
    fseed%usezname = 2
    if (fseed%nat == 0) then
       useed%errcode = 7
       call f_c_string("No atoms in the cell.",useed%errmsg,254)
       return
    end if

    ! atomic unit conversion
    if (useed%atunits == 0 .and. .not.fseed%ismolecule) then
       ! bohr
       r = matinv(fseed%crys2car)
       do i = 1, fseed%nat
          fseed%x(:,i) = matmul(r,fseed%x(:,i))
       end do
    elseif (useed%atunits == 1) then
       ! angstrom
       fseed%x = fseed%x / bohrtoa
       if (.not.fseed%ismolecule) then
          r = matinv(fseed%crys2car)
          do i = 1, fseed%nat
             fseed%x(:,i) = matmul(r,fseed%x(:,i))
          end do
       end if
    end if

    ! space group
    call c_f_string(useed%strspg,aux)
    if (len_trim(aux) > 0) then
       call fseed%spgs_wrap(aux,.false.)
       if (fseed%havesym == 0) then
          useed%errcode = 8
          call f_c_string("Incorrect space group.",useed%errmsg,254)
          return
       end if
    end if
    if (.not.fseed%ismolecule .and. fseed%havesym == 0) then
       fseed%findsym = 1
    end if

    ! build the new crystal structure
    call cr%struct_new(fseed,.false.)
    if (fseed%nat == 0) then
       useed%errcode = 9
       call f_c_string("Could not initialize the crystal.",useed%errmsg,254)
       return
    end if

    ! initialize and update the scene
    if (cr%isinit) then
       ispreview = preview
       if (.not.preview) then
          call initialize_structure()
       else
          call cr%struct_fill(.true.,-1,.false.,.true.,.false.)
       end if
       call update_scene()
       new_structure = 1
    end if

  end function new_structure

  ! Save the data for the current structure and load an alternative
  ! structure that may be accepted or not.
  function preview_structure(useed) bind(c)
    use crystalseedmod, only: crystalseed
    use crystalmod, only: cr
    type(c_crystalseed), intent(inout) :: useed
    integer(c_int) :: preview_structure

    logical :: clearlater
    integer :: i

    if (.not.issaved) then
       call save_state()
       clearlater = .true.
    else
       clearlater = .false.
    end if
    preview_structure = new_structure(useed,logical(.true.,c_bool))
    if (preview_structure == 0 .and. clearlater) &
       call clear_state()

  end function preview_structure

  ! accept the previewed structure 
  subroutine accept_previewed_structure() bind(c)

    ! clear the save state and properly initialize the structure
    call clear_state()
    call initialize_structure()
    ispreview = .false.
       
  end subroutine accept_previewed_structure

  ! reject the previewed structure 
  subroutine reject_previewed_structure() bind(c)
    
    call load_state()
    call clear_state()
    ispreview = .false.
    call update_scene()

  end subroutine reject_previewed_structure

  ! The cr variable contains a new structure. Use this routine to
  ! initialize the rest of the modules to prepare it for usage.
  subroutine initialize_structure()
    use fields, only: nprops, integ_prop, f, type_grid, itype_fval, itype_lapval,&
       fields_integrable_report
    use grd_atomic, only: grda_init
    use crystalmod, only: cr
    use tools_io, only: uout, string
    use global, only: refden, gradient_mode, INT_radquad_errprop_default, INT_radquad_errprop
    use autocp, only: init_cplist

    if (.not.cr%isinit) return

    ! fill environments, asterisms, nearest neighbors
    call cr%struct_fill(.true.,-1,.false.,.true.,.false.)
    ! print some information about the structure
    call cr%struct_report()
    ! initialize the radial densities
    call grda_init(.true.,.true.,.false.)
    ! header and change refden
    write (uout,'("* Field number ",A," is now REFERENCE."/)') string(0)
    refden = 0
    call init_cplist(.true.)
    ! define second integrable property as the valence charge.
    nprops = max(2,nprops)
    integ_prop(2)%used = .true.
    integ_prop(2)%itype = itype_fval
    integ_prop(2)%fid = 0
    integ_prop(2)%prop_name = "Pop"
    ! define third integrable property as the valence laplacian.
    nprops = max(3,nprops)
    integ_prop(3)%used = .true.
    integ_prop(3)%itype = itype_lapval
    integ_prop(3)%fid = 0
    integ_prop(3)%prop_name = "Lap"
    ! reset defaults for qtree
    if (f(refden)%type == type_grid) then
       gradient_mode = 1
       if (INT_radquad_errprop_default) INT_radquad_errprop = 2
    else
       gradient_mode = 2
       if (INT_radquad_errprop_default) INT_radquad_errprop = 3
    end if

    ! report
    call fields_integrable_report()

  end subroutine initialize_structure

  ! Calculate critical points for the current field
  subroutine call_auto() bind (c)
    use crystalmod, only: cr
    use autocp, only: init_cplist, autocritic

    if (.not.cr%isinit) return
    call autocritic("")
    call update_scene()

  end subroutine call_auto

  ! Update the data in the module variables - makes it available
  ! to the C++ code
  subroutine update_scene() bind(c)
    use fragmentmod, only: fragment_merge_array
    use varbas, only: ncpcel, cpcel
    use crystalmod, only: cr
    use tools_math, only: norm, cross
    use types, only: fragment
    use param, only: atmcov, jmlcol, maxzat
    
    integer :: i, j, k, iz, nup, ix(3)
    real*8 :: x1(3), x2(3), xcm(3)
    integer :: nmol, idxi, idxj
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)

    real*8, parameter :: bondthickness = 0.05d0
    real*8, parameter :: bondcolor(3) = (/0.d0,0.d0,0.d0/)
    real*8, parameter :: cpradius = 0.2d0
    real*8, parameter :: cellcolor(3) = (/1.0d0,0.7d0,0.0d0/)
    real*8, parameter :: cellxcolor(3) = (/1.0d0,0.0d0,0.0d0/)
    real*8, parameter :: cellycolor(3) = (/0.0d0,1.0d0,0.0d0/)
    real*8, parameter :: cellzcolor(3) = (/0.0d0,0.0d0,1.0d0/)
    real*8, parameter :: cellthick = 0.05d0

    if (.not.cr%isinit) return

    ! global flags
    isinit = cr%isinit
    ismolecule = cr%ismolecule

    ! Calculate the bounding box
    box_xmin = 1e30
    box_xmax = -1e30
    xcm = 0.
    if (.not.ispreview) then
       do i = 1, ncpcel
          x1 = cpcel(i)%r + cr%molx0
          box_xmin = min(x1,box_xmin)
          box_xmax = max(x1,box_xmax)
          xcm = xcm + x1
       end do
       xcm = xcm / ncpcel
    else
       do i = 1, cr%ncel
          x1 = cr%atcel(i)%r + cr%molx0
          box_xmin = min(x1,box_xmin)
          box_xmax = max(x1,box_xmax)
          xcm = xcm + x1
       end do
       xcm = xcm / cr%ncel
    end if
    box_xcm = xcm
    box_xmaxlen = max(maxval(box_xmax - box_xmin),1d0)

    ! molecules -> whole molecule. crystals -> border and
    ! molmotif if it is molecular crystal.
    ix = 1
    if (cr%ismolecule) then
       fr = cr%listatoms_cells(ix,.false.)
    else
       fr = cr%listatoms_cells(ix,.true.)
       if (all(cr%moldiscrete(1:cr%nmol))) then
          call cr%listmolecules(fr,nmol,fr0,isdiscrete)
          fr = fragment_merge_array(fr0)
          deallocate(fr0,isdiscrete)
       end if
    end if

    ! Allocate space for atoms
    if (associated(at_f)) deallocate(at_f)
    allocate(at_f(max(fr%nat,1)))

    ! For now, just go ahead and represent the whole cell/molecule
    nat = fr%nat
    do i = 1, nat
       iz = cr%at(fr%at(i)%idx)%z
       at_f(i)%z = iz
       at_f(i)%b%r = fr%at(i)%r + cr%molx0 - xcm
       if (atmcov(iz) > 1) then
          at_f(i)%b%rad = atmcov(iz)
       else
          at_f(i)%b%rad = 2d0*atmcov(iz)
       end if
       at_f(i)%b%rgb = real(jmlcol(:,iz),4) / 255.
       call f_c_string(cr%at(fr%at(i)%idx)%name,at_f(i)%name,11)
    end do
    at = c_loc(at_f)

    ! Allocate space for bonds. For now, only bonds between atoms in
    ! the main cell.
    nbond = 0
    do i = 1, fr%nat
       idxi = fr%at(i)%cidx
       do k = 1, cr%nstar(idxi)%ncon
          do j = i+1, fr%nat
             idxj = fr%at(j)%cidx
             if (cr%nstar(idxi)%idcon(k) == idxj .and. all(fr%at(i)%lvec+cr%nstar(idxi)%lcon(:,k) == fr%at(j)%lvec)) &
                nbond = nbond + 1
          end do
       end do
    end do
    if (associated(bond_f)) deallocate(bond_f)
    allocate(bond_f(max(nbond,1)))
    nbond = 0
    do i = 1, fr%nat
       idxi = fr%at(i)%cidx
       do k = 1, cr%nstar(idxi)%ncon
          do j = i+1, fr%nat
             idxj = fr%at(j)%cidx
             if (cr%nstar(idxi)%idcon(k) == idxj .and. all(fr%at(i)%lvec+cr%nstar(idxi)%lcon(:,k) == fr%at(j)%lvec)) then
                nbond = nbond + 1
                bond_f(nbond)%i1 = i-1
                bond_f(nbond)%i2 = j-1
                x1 = fr%at(i)%r + cr%molx0 - xcm
                x2 = fr%at(j)%r + cr%molx0 - xcm
                bond_f(nbond)%s = stick_from_endpoints(x1,x2,bondthickness,bondcolor)
             end if
          end do
       end do
    end do
    bond = c_loc(bond_f)

    if (.not.ispreview) then
       ! Allocate space for critical points
       ncritp = ncpcel - cr%ncel
       if (associated(critp_f)) deallocate(critp_f)
       allocate(critp_f(max(ncritp,1)))
       
       ! For now, just go ahead and represent the whole cell
       j = 0
       do i = cr%ncel+1, ncpcel
          j = j + 1
          iz = maxzat + 1 + cpcel(i)%typind
          critp_f(j)%type = cpcel(i)%typ
          critp_f(j)%b%r = cpcel(i)%r + cr%molx0 - xcm
          critp_f(j)%b%rgb = real(jmlcol(:,iz),4) / 255.
          critp_f(j)%b%rad = cpradius
          call f_c_string(cpcel(i)%name,critp_f(j)%name,11)
       end do
       critp = c_loc(critp_f)
    else
       critp = C_NULL_PTR
    end if

    ! unit cell and lattice vectors
    cell_x0 = cr%molx0
    do i = 1, 3
       cell_lat(:,i) = cr%crys2car(:,i)
    end do
    x1 = cr%molx0 - xcm
    x2 = x1 + cr%crys2car(:,1)+cr%crys2car(:,2)+cr%crys2car(:,3)
    cell_nstick = 12;
    cell_s(1) = stick_from_endpoints(x1,x1+cr%crys2car(:,1),cellthick,cellxcolor)
    cell_s(2) = stick_from_endpoints(x1,x1+cr%crys2car(:,2),cellthick,cellycolor)
    cell_s(3) = stick_from_endpoints(x1,x1+cr%crys2car(:,3),cellthick,cellzcolor)
    cell_s(4) = stick_from_endpoints(x1+cr%crys2car(:,1),x1+cr%crys2car(:,1)+cr%crys2car(:,2),cellthick,cellcolor)
    cell_s(5) = stick_from_endpoints(x1+cr%crys2car(:,1),x1+cr%crys2car(:,1)+cr%crys2car(:,3),cellthick,cellcolor)
    cell_s(6) = stick_from_endpoints(x1+cr%crys2car(:,2),x1+cr%crys2car(:,2)+cr%crys2car(:,1),cellthick,cellcolor)
    cell_s(7) = stick_from_endpoints(x1+cr%crys2car(:,2),x1+cr%crys2car(:,2)+cr%crys2car(:,3),cellthick,cellcolor)
    cell_s(8) = stick_from_endpoints(x1+cr%crys2car(:,3),x1+cr%crys2car(:,3)+cr%crys2car(:,1),cellthick,cellcolor)
    cell_s(9) = stick_from_endpoints(x1+cr%crys2car(:,3),x1+cr%crys2car(:,3)+cr%crys2car(:,2),cellthick,cellcolor)
    cell_s(10) = stick_from_endpoints(x1+cr%crys2car(:,1)+cr%crys2car(:,2),x2,cellthick,cellcolor)
    cell_s(11) = stick_from_endpoints(x1+cr%crys2car(:,1)+cr%crys2car(:,3),x2,cellthick,cellcolor)
    cell_s(12) = stick_from_endpoints(x1+cr%crys2car(:,2)+cr%crys2car(:,3),x2,cellthick,cellcolor)
    box_xmaxclen = norm(x2-x1)

  end subroutine update_scene

  ! Clear the scene and (optionally) unload the crystal structure
  subroutine clear_scene(unload) bind(c)
    use crystalmod, only: cr
    logical(c_bool), intent(in), value :: unload

    ! global flags
    isinit = .false.
    ismolecule = .false.

    ! atoms
    if (associated(at_f)) deallocate(at_f)
    nat = 0
    at = C_NULL_PTR
    
    ! bonds
    if (associated(bond_f)) deallocate(bond_f)
    nbond = 0
    bond = C_NULL_PTR
    
    ! critical points
    if (associated(critp_f)) deallocate(critp_f)
    ncritp = 0
    critp = C_NULL_PTR
    
    ! cell
    cell_nstick = 0

    ! bounding box
    box_xmin = 0.
    box_xmax = 0.
    box_xcm = 0. 
    box_xmaxlen = 0. 

    ! unload the crystal structure
    if (unload) &
       call cr%end()

  end subroutine clear_scene

  function get_text_info(imode) bind(c) result(txt)
    use crystalmod, only: cr, pointgroup_info, laue_string, holo_string
    use fragmentmod, only: fragment_cmass
    use tools_io, only: string, ioj_center, ioj_left
    use param, only: bohrtoa, maxzat
    integer(c_int), intent(in), value :: imode
    type(c_ptr) :: txt
    
    character*1, parameter :: nl = new_line('a')

    integer :: i
    character(len=:), allocatable :: txto, aux
    character(len=3) :: schpg
    integer :: holo, laue, nelec, idx
    character*255 :: strout(cr%neqv)
    real*8 :: xcm(3), x(3)

    txt = c_null_ptr
    select case (imode)
    case(0) ! unselected
    case(1) ! General (unavailable)
    case(2) ! General (molecule)
       nelec = 0
       do i = 1, cr%nneq
          if (cr%at(i)%z >= maxzat) cycle
          nelec = nelec + cr%at(i)%z * cr%at(i)%mult
       end do
       txto = &
          "From: " // string(cr%file) // nl //&
          "Number of atoms: "  // string(cr%ncel) // nl //&
          "Number of electrons: " // string(nelec) // nl 
       txt = f_c_string_dup(txto)
    case(3) ! General (crystal) 
       nelec = 0
       do i = 1, cr%nneq
          if (cr%at(i)%z >= maxzat) cycle
          nelec = nelec + cr%at(i)%z * cr%at(i)%mult
       end do
       txto = &
 "From: " // string(cr%file) // nl //&
 "Lattice parameters (bohr): " // string(cr%aa(1),'f',decimal=4) // " " // string(cr%aa(2),'f',decimal=4) // " " // string(cr%aa(3),'f',decimal=4) // nl //&
 "Lattice parameters (ang): " // string(cr%aa(1)*bohrtoa,'f',decimal=4) // " " // string(cr%aa(2)*bohrtoa,'f',decimal=4) // " " // string(cr%aa(3)*bohrtoa,'f',decimal=4) // nl //&
 "Lattice angles: " // string(cr%bb(1),'f',decimal=3) // " " // string(cr%bb(2),'f',decimal=3) // " " // string(cr%bb(3),'f',decimal=3) // nl //&
 "Atoms in the unit cell: "  // string(cr%ncel) // nl //&
 "Atoms in the asymmetric unit: " // string(cr%nneq) // nl //&
 "Z: " // string(cr%ncel/cr%nneq) // nl //&
 "Atoms in the environment: " // string(cr%nenv) // nl //&
 "Number of electrons: " // string(nelec) // nl //&
 "Cell volume (bohr^3): " // string(cr%omega,'f',decimal=3) // nl //&
 "Cell volume (ang^3): " // string(cr%omega * bohrtoa**3,'f',decimal=3) // nl 

       if (len_trim(cr%spg%choice) > 0) then
          aux = txto // "Space group (Hermann-Mauguin): " // string(cr%spg%international_symbol) // " (number " // &
             string(cr%spg%spacegroup_number) // ", setting " // string(cr%spg%choice) // ")" // nl
       else
          aux = txto // "Space group (Hermann-Mauguin): " // string(cr%spg%international_symbol) // " (number " //  &
             string(cr%spg%spacegroup_number) // ")" // nl
       end if
       call pointgroup_info(cr%spg%pointgroup_symbol,schpg,holo,laue)
       txto = aux // "Space group (Hall): " // string(cr%spg%hall_symbol) // " (number " // &
          string(cr%spg%hall_number) // ")" // nl // &
          "Point group (Hermann-Mauguin): " // string(cr%spg%pointgroup_symbol) // nl // &
          "Point group (Schoenflies): " // string(schpg) // nl // &
          "Holohedry: " // string(holo_string(holo)) // nl // &
          "Laue class: " // string(laue_string(laue)) // nl // &
          "Orthogonal? " // string(cr%isortho)
       txt = f_c_string_dup(txto)
    case(4) ! Asymmetric unit (crystal)
       txto = "# nat      x             y             z          name  mult  Z  " // nl
       do i = 1, cr%nneq
          aux = txto // "  " // string(i,3,ioj_center) // string(cr%at(i)%x(1),'f',length=14,decimal=10,justify=3) //&
             string(cr%at(i)%x(2),'f',length=14,decimal=10,justify=3) // string(cr%at(i)%x(3),'f',length=14,decimal=10,justify=3) //&
             string(cr%at(i)%name,10,ioj_center) // string(cr%at(i)%mult,4,ioj_center) // string(cr%at(i)%z,4,ioj_center) // nl
          txto = aux
       end do
       txt = f_c_string_dup(txto)
    case(5) ! Unit cell (crystal)
       txto = "# nat        position (cryst. coords.)            name    Z" // nl
       do i = 1, cr%ncel
          idx = cr%atcel(i)%idx
          aux = txto // "  " // string(i,3,ioj_center) // string(cr%atcel(i)%x(1),'f',length=14,decimal=10,justify=3) //&
             string(cr%atcel(i)%x(2),'f',length=14,decimal=10,justify=3) // string(cr%atcel(i)%x(3),'f',length=14,decimal=10,justify=3) //&
             string(cr%at(idx)%name,10,ioj_center) // string(cr%at(idx)%z,4,ioj_center) // nl
          txto = aux
       end do
       txt = f_c_string_dup(txto)
    case(6) ! Symmetry (crystal) 
       if (len_trim(cr%spg%choice) > 0) then
          aux = "Space group (Hermann-Mauguin): " // string(cr%spg%international_symbol) // " (number " // &
             string(cr%spg%spacegroup_number) // ", setting " // string(cr%spg%choice) // ")" // nl
       else
          aux = "Space group (Hermann-Mauguin): " // string(cr%spg%international_symbol) // " (number " //  &
             string(cr%spg%spacegroup_number) // ")" // nl
       end if
       call pointgroup_info(cr%spg%pointgroup_symbol,schpg,holo,laue)
       txto = aux // "Space group (Hall): " // string(cr%spg%hall_symbol) // " (number " // &
          string(cr%spg%hall_number) // ")" // nl // &
          "Point group (Hermann-Mauguin): " // string(cr%spg%pointgroup_symbol) // nl // &
          "Point group (Schoenflies): " // string(schpg) // nl // &
          "Holohedry: " // string(holo_string(holo)) // nl // &
          "Laue class: " // string(laue_string(laue)) // nl
       aux = txto // "List of symmetry operations (" // string(cr%neqv) // "):" // nl
       txto = aux

       call cr%struct_report_symxyz(strout)
       do i = 1, cr%neqv
          aux = txto // "   " // string(i) // ": " // string(strout(i)) // nl
          txto = aux
       enddo
       aux = txto // "List of centering vectors (" // string(cr%ncv) // "):" // nl
       txto = aux
       do i = 1, cr%ncv
          aux = txto // "   " // string(i) // ": " // string(cr%cen(1,i),'f',length=8,decimal=5) // " " //&
             string(cr%cen(2,i),'f',length=8,decimal=5) // " " // string(cr%cen(3,i),'f',length=8,decimal=5) // nl
          txto = aux
       enddo
       txt = f_c_string_dup(txto)
    case(7) ! Fragments (crystal)
       txto = "# Id nat            Center of mass       Discrete? " // nl
       do i = 1, cr%nmol
          xcm = cr%c2x(fragment_cmass(cr%mol(i)))
          aux = txto // "  " // string(i,3,ioj_left) // " " // string(cr%mol(i)%nat,4,ioj_left) // " " //&
             string(xcm(1),'f',10,6,3) // " "  // string(xcm(2),'f',10,6,3) // " "  // &
             string(xcm(3),'f',10,6,3) // " " // string(cr%moldiscrete(i)) // nl
          txto = aux
       end do
       txt = f_c_string_dup(txto)
    case(8) ! Atoms (molecule)
       txto = "# nat           position (angstrom)               name    Z" // nl
       do i = 1, cr%ncel
          idx = cr%atcel(i)%idx
          x = (cr%atcel(i)%r + cr%molx0) * bohrtoa
          aux = txto // "  " // string(i,3,ioj_center) // string(x(1),'f',length=14,decimal=10,justify=3) //&
             string(x(2),'f',length=14,decimal=10,justify=3) // string(x(3),'f',length=14,decimal=10,justify=3) //&
             string(cr%at(idx)%name,10,ioj_center) // string(cr%at(idx)%z,4,ioj_center) // nl
          txto = aux
       end do
       txt = f_c_string_dup(txto)
    case(9) ! Fragments (molecule)
       txto = "# Id nat            Center of mass       Discrete? " // nl
       do i = 1, cr%nmol
          xcm = (fragment_cmass(cr%mol(i))+cr%molx0) * bohrtoa
          aux = txto // "  " // string(i,3,ioj_left) // " " // string(cr%mol(i)%nat,4,ioj_left) // " " //&
             string(xcm(1),'f',10,6,3) // " "  // string(xcm(2),'f',10,6,3) // " "  // &
             string(xcm(3),'f',10,6,3) // " " // string(cr%moldiscrete(i)) // nl
          txto = aux
       end do
       txt = f_c_string_dup(txto)
    end select

  end function get_text_info

  function get_seed_from_current_structure() result(useed) bind(c)
    use crystalmod, only: cr
    use global, only: rborder_def
    use tools_io, only: string
    use param, only: bohrtoa
    type(c_crystalseed) :: useed

    character(len=:), allocatable :: aux, rr, aa, bb, at
    integer :: i, j

    if (.not.cr%isinit) then
       useed%type = -1
       useed%errcode = 0
       return
    end if
    if (cr%ismolecule) then
       useed%type = 0
    else
       useed%type = 1
    end if
    useed%achoice = 1
    
    aa = ""
    bb = ""
    rr = ""
    do i = 1, 3
       aux = trim(aa) // " " // string(cr%aa(i),'f',15,10)
       aa = aux
       aux = trim(bb) // " " // string(cr%bb(i),'f',10,5)
       bb = aux
       do j = 1, 3
          aux = trim(rr) // " " // string(cr%crys2car(i,j),'f',10,5)
          rr = aux
       end do
       aux = rr // nl
       rr = aux
    end do
    call f_c_string(aa,useed%straa)
    call f_c_string(bb,useed%strbb)
    call f_c_string(rr,useed%strrr)
    useed%aaunits = 0
    useed%rrunits = 0
    
    at = ""
    do i = 1, cr%ncel
       aux = trim(at) // string(cr%at(cr%atcel(i)%idx)%name,2)
       at = aux
       do j = 1, 3
          if (cr%ismolecule) then
             aux = trim(at) // " " // string(cr%atcel(i)%r(j)*bohrtoa,'f',15,10)
          else
             aux = trim(at) // " " // string(cr%atcel(i)%x(j),'f',15,10)
          end if
          at = aux
       end do
       aux = trim(at) // nl
       at = aux
    end do
    call f_c_string(at,useed%strat)
    if (cr%ismolecule) then
       useed%atunits = 1
    else
       useed%atunits = 2
    end if

    aux = ""
    call f_c_string(aux,useed%strspg)
    aux = ""
    call f_c_string(aux,useed%errmsg)

    useed%molcubic = .false.
    useed%errcode = 0
    useed%molborder = rborder_def * bohrtoa
    useed%borunits = 1

  end function get_seed_from_current_structure

  subroutine set_library_file(filename0,type) bind(c)
    use global, only: clib_file, mlib_file
    type(c_ptr), intent(in) :: filename0
    integer(c_int), value :: type
    
    character(len=:), allocatable :: filename

    ! transform to fortran string
    filename = c_string_value(filename0)
    if (type == 1) then
       clib_file = filename
    elseif (type == 2) then
       mlib_file = filename
    end if
    call critic2_initialize_library()

  end subroutine set_library_file

  ! Build a c_stick from the two endpoints, thickness, and rgb
  function stick_from_endpoints(x1,x2,thick,rgb) result(stick)
    use tools_math, only: norm, cross
    real*8, intent(in) :: x1(3), x2(3), thick, rgb(3)
    type(c_stick) :: stick
    
    real*8 :: xd(3), normd, uz(3), u(3), nu, ca, sa

    ! points
    stick%r1 = x1
    stick%r2 = x2
    stick%rmid = 0.5d0 * (x1+x2)
    
    ! thickness and rgb
    stick%thick = thick
    stick%rgb = rgb

    ! length
    xd = x1 - x2
    normd = norm(xd)
    xd = xd / normd
    stick%length = normd
    
    ! rotation matrix
    uz = (/0d0, 0d0, 1d0/)
    u = cross(uz,xd)
    nu = norm(u)
    if (nu < 1d-8) then
       if (xd(3) > 0d0) then
          ca = 1d0
          sa = 0d0
          u = (/1d0,0d0,0d0/)
       else
          ca = -1d0
          sa = 0d0
          u = (/1d0,0d0,0d0/)
       end if
    else
       u = u / nu
       ca = xd(3)
       sa = sqrt(1d0 - ca**2)
    end if

    stick%rot(1,1) = ca + u(1)*u(1)*(1d0 - ca)
    stick%rot(2,1) = u(1)*u(2)*(1d0 - ca) - u(3)*sa
    stick%rot(3,1) = u(1)*u(3)*(1d0 - ca) + u(2)*sa
    stick%rot(4,1) = 0d0

    stick%rot(1,2) = u(1)*u(2)*(1 - ca) + u(3)*sa
    stick%rot(2,2) = ca + u(2)*u(2)*(1d0 - ca)
    stick%rot(3,2) = u(2)*u(3)*(1d0 - ca) - u(1)*sa
    stick%rot(4,2) = 0d0

    stick%rot(1,3) = u(1)*u(3)*(1d0 - ca) - u(2)*sa
    stick%rot(2,3) = u(3)*u(2)*(1d0 - ca) + u(1)*sa
    stick%rot(3,3) = ca + u(3)*u(3)*(1d0 - ca)
    stick%rot(4,3) = 0d0

    stick%rot(1,4) = 0d0
    stick%rot(2,4) = 0d0
    stick%rot(3,4) = 0d0
    stick%rot(4,4) = 1d0
             
  end function stick_from_endpoints

  subroutine denewline(str)
    character*(*), intent(inout) :: str
    
    integer :: i

    do i = 1, len_trim(str)
       if (str(i:i) == nl) str(i:i) = " "
    end do

  end subroutine denewline

  ! put the current scene into the save state
  subroutine save_state()
    use crystalmod, only: cr

    isinit_ = isinit
    ismolecule_ = ismolecule
    nat_ = nat
    if (associated(at_f_)) deallocate(at_f_)
    if (associated(at_f)) then
       allocate(at_f_(size(at_f,1)))
       at_f_ = at_f
    end if
    nbond_ = nbond
    if (associated(bond_f_)) deallocate(bond_f_)
    if (associated(bond_f)) then
       allocate(bond_f_(size(bond_f,1)))
       bond_f_ = bond_f
    end if
    ncritp_ = ncritp
    if (associated(critp_f_)) deallocate(critp_f_)
    if (associated(critp_f)) then
       allocate(critp_f_(size(critp_f,1)))
       critp_f_ = critp_f
    end if
    cell_x0_ = cell_x0
    cell_lat_ = cell_lat
    cell_nstick_ = cell_nstick
    cell_s_ = cell_s
    box_xmin_ = box_xmin
    box_xmax_ = box_xmax
    box_xcm_ = box_xcm
    box_xmaxlen_ = box_xmaxlen
    box_xmaxclen_ = box_xmaxclen
    cr_ = cr
    issaved = .true.

  end subroutine save_state

  ! load the save state, if available
  subroutine load_state()
    use crystalmod, only: cr

    if (issaved) then
       isinit = isinit_
       ismolecule = ismolecule_
       nat = nat_
       if (associated(at_f)) deallocate(at_f)
       if (associated(at_f_)) then
          allocate(at_f(size(at_f_,1)))
          at_f = at_f_
       end if
       nbond = nbond_
       if (associated(bond_f)) deallocate(bond_f)
       if (associated(bond_f_)) then
          allocate(bond_f(size(bond_f_,1)))
          bond_f = bond_f_
       end if
       ncritp = ncritp_
       if (associated(critp_f)) deallocate(critp_f)
       if (associated(critp_f_)) then
          allocate(critp_f(size(critp_f_,1)))
          critp_f = critp_f_
       end if
       cell_x0 = cell_x0_
       cell_lat = cell_lat_
       cell_nstick = cell_nstick_
       cell_s = cell_s_
       box_xmin = box_xmin_
       box_xmax = box_xmax_
       box_xcm = box_xcm_
       box_xmaxlen = box_xmaxlen_
       box_xmaxclen = box_xmaxclen_
       cr = cr_
    end if

  end subroutine load_state

  ! clear the save state
  subroutine clear_state()

    if (associated(at_f_)) deallocate(at_f_)
    if (associated(bond_f_)) deallocate(bond_f_)
    if (associated(critp_f_)) deallocate(critp_f_)
    issaved = .false.

  end subroutine clear_state
  
  subroutine run_critic2_command(command0) bind(c)
    use tricks, only: trick
    use stm, only: stm_driver
    use xdm, only: xdm_driver
    use hirshfeld, only: hirsh_props_grid
    use qtree, only: qtree_integration, qtree_setsphfactor
    use bisect, only: basinplot, bundleplot, sphereintegrals, integrals
    use integration, only: intgrid_driver
    use flux, only: fluxprint
    use autocp, only: init_cplist, autocritic, cpreport
    use nci, only: nciplot
    use rhoplot, only: rhoplot_point, rhoplot_line, rhoplot_plane, rhoplot_cube,&
       rhoplot_grdvec
    use fields, only: fieldname_to_idx, goodfield, f, fused, type_grid, nprops,&
       integ_prop, itype_fval, itype_lapval, fields_init, fields_end, &
       fields_load, fields_unload, setfield, fieldinfo, benchmark, &
       fields_integrable, fields_pointprop, testrmt, listfields, listfieldalias,&
       fields_integrable_report
    use varbas, only: varbas_end, varbas_identify
    use grd_atomic, only: grda_init, grda_end
    use struct_drivers, only: struct_crystal_input, struct_newcell, struct_molcell,&
       struct_clearsym, struct_charges, struct_write, struct_powder, struct_rdf,&
       struct_compare, struct_environ, struct_packing, struct_atomlabel
    use crystalmod, only: cr
    use wfn_private, only: wfn_end
    use pi_private, only: pi_end
    use spgs, only: spgs_init
    use global, only: fileroot, quiet, refden, eval_next, gradient_mode,&
       int_radquad_errprop_default, int_radquad_errprop, global_init,&
       initial_banner, help_me, config_write, critic_clearvariable,&
       critic_setvariables,global_set_defaults, iunit, iunit_isdef,&
       iunit_ang, iunit_bohr
    use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
       fcflags, cc, cflags, ldflags, enable_debug, package
    use graphics, only: graphics_init
    use arithmetic, only: listvariables
    use tools_io, only: uout, ucopy, uin, getline, lgetword, equal, faterr,&
       ferror, getword, string, nwarns, ncomms, ioinit, stdargs, tictac, &
       start_clock, print_clock
    use param, only: param_init
    type(c_ptr), intent(in), value :: command0

    character(len=:), allocatable :: command, subline, word
    integer :: level, plevel
    integer :: lp, lpold
    integer :: i, id, nn, ismoli
    logical :: ll1, ok
    real*8 :: rdum

    ! transform to fortran string
    command = c_string_value(command0)

    ! parse
    lp = 1
    word = lgetword(command,lp)
    subline = command(lp:)
    ! crystal
    if (equal(word,'crystal') .or. equal(word,'molecule')) then
       ! there is a previous crystal structure, clean up...
       if (cr%isinit) call clean_structure()

       if (equal(word,'crystal')) then
          ismoli = 0
       else
          ismoli = 1
       end if

       ! change default output units
       if (iunit_isdef) then
          if (equal(word,'molecule')) then
             iunit = iunit_ang
          else
             iunit = iunit_bohr
          end if
       end if

       ! read the crystal enviornment
       call struct_crystal_input(cr,subline,ismoli,.true.)

       if (cr%isinit) then
          ispreview = .false.
          call initialize_structure()
          call update_scene()
       end if

    elseif (equal(word,'newcell')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before newcell',faterr,command,syntax=.true.)
          return
       end if
       call struct_newcell(subline)
       if (cr%isinit) then
          ispreview = .false.
          call initialize_structure()
          call update_scene()
       end if

       ! molcell
    elseif (equal(word,'molcell')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before molcell',faterr,command,syntax=.true.)
          return
       end if
       call struct_molcell(subline)
       call update_scene()

       ! clearsym/clearsymm
    elseif (equal(word,'clearsym') .or. equal(word,'clearsymm')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before clearsym',faterr,command,syntax=.true.)
          return
       end if
       call struct_clearsym() 
       call check_no_extra_word(ok)
       if (.not.ok) return

       ! q/qat, zpsp, nocore
    elseif (equal(word,'q') .or. equal(word,'qat') &
       .or. equal(word,'zpsp') .or. equal(word,'nocore')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before q/qat/zpsp/nocore',faterr,command,syntax=.true.)
          return
       end if
       call struct_charges(command,ok)
       if (ok) then
          ll1 = equal(word,'zpsp') .or. equal(word,'nocore')
          call grda_init(ll1,.not.ll1,.true.)
       end if

       ! atomlabel
    elseif (equal(word,'atomlabel')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before atomlabel',faterr,command,syntax=.true.)
          return
       end if
       call struct_atomlabel(cr,subline)

       ! write
    elseif (equal(word,'write')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before write',faterr,command,syntax=.true.)
          return
       end if
       call struct_write(cr,subline)

       ! load
    elseif (equal(word,'load')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before load',faterr,command,syntax=.true.)
          return
       end if
       call fields_load(subline,id,ok)
       if (ok) then
          if (refden == 0) call set_reference(id)
       else
          call fields_unload(id)
       end if

       ! unload
    elseif (equal(word,'unload')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before unload',faterr,command,syntax=.true.)
          return
       end if
       lpold = lp
       word = getword(command,lp)
       nn = fieldname_to_idx(word)
       if (nn >= 0) then
          if (.not.goodfield(nn)) then
             call ferror('critic2','wrong field in UNLOAD',faterr,command,syntax=.true.)
             return
          end if
          if (nn == 0) then
             call ferror('critic2','can not unload the promolecular density',faterr,command,syntax=.true.)
             return
          end if
          call fields_unload(nn)
       else
          lp = lpold
          word = lgetword(command,lp)
          if (equal(word,'all')) then
             do i = 1, ubound(f,1)
                if (fused(i)) call fields_unload(i)
             end do
          else
             call ferror('critic2','Unknown keyword in UNLOAD',faterr,command,syntax=.true.)
          endif
       end if

       ! powder
    elseif (equal(word,'powder')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before powder',faterr,command,syntax=.true.)
          return
       end if
       call struct_powder(command(lp:),cr)

       ! rdf
    elseif (equal(word,'rdf')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before rdf',faterr,command,syntax=.true.)
          return
       end if
       call struct_rdf(command(lp:),cr)

       ! compare
    elseif (equal(word,'compare')) then
       call struct_compare(command(lp:))

       ! setfield
    elseif (equal(word,'setfield')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before setfield',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','wrong field in setfield',faterr,command,syntax=.true.)
          return
       end if
       call setfield(f(id),id,command(lp:),ok)
       call fieldinfo(id,.false.,.true.)

       ! reference
    elseif (equal(word,'reference')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before reference',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) then
          call ferror('critic2','unknown field in REFERENCE',faterr,command,syntax=.true.)
          return
       end if
       if (.not.goodfield(id)) then
          call ferror('critic2','REFERENCE: field is not allocated',faterr,command,syntax=.true.)
          return
       end if
       call set_reference(id)
       call check_no_extra_word(ok)

       ! point
    elseif (equal(word,'point')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before point',faterr,command,syntax=.true.)
          return
       end if
       call rhoplot_point(subline)

       ! line
    elseif (equal(word,'line')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before line',faterr,command,syntax=.true.)
          return
       end if
       call rhoplot_line(subline)

       ! plane
    elseif (equal(word,'plane')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before plane',faterr,command,syntax=.true.)
          return
       end if
       call rhoplot_plane(subline)

       ! cube
    elseif (equal(word,'cube')) then
       if (.not. cr%isinit) then 
          call ferror('critic2','need crystal before cube',faterr,command,syntax=.true.)
          return
       end if
       call rhoplot_cube(subline)

       ! grdvec
    elseif (equal(word,'grdvec')) then
       if (.not. cr%isinit) then
          call ferror('critic','need crystal before grdvec',faterr,command,syntax=.true.)
          return
       end if
       call rhoplot_grdvec()

       ! nciplot
    elseif (equal(word,'nciplot')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before nciplot',faterr,command,syntax=.true.)
          return
       end if
       call nciplot()

       ! benchmark 
    elseif (equal(word,'benchmark')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before benchmark',faterr,command,syntax=.true.)
          return
       end if
       ok = eval_next(nn,command,lp)
       if (.not. ok) nn = 10000
       call check_no_extra_word(ok)
       if (.not.ok) return
       call benchmark(nn)

       ! auto
    elseif (equal(word,'auto')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before auto',faterr,command,syntax=.true.)
          return
       end if
       call autocritic(subline)

       ! cpreport
    elseif (equal(word,'cpreport')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before cpreport',faterr,command,syntax=.true.)
          return
       end if
       call cpreport(subline)

       ! fluxprint
    elseif (equal(word,'fluxprint')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before fluxprint',faterr,command,syntax=.true.)
          return
       end if
       call fluxprint()

       ! integrable
    elseif (equal(word,'integrable')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before integrable',faterr,command,syntax=.true.)
          return
       end if
       call fields_integrable(subline)

       ! pointprop
    elseif (equal(word,'pointprop')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before pointprop',faterr,command,syntax=.true.)
          return
       end if
       call fields_pointprop(subline)

       ! basinplot
    elseif (equal(word,'basinplot')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before basinplot',faterr,command,syntax=.true.)
          return
       end if
       call basinplot(subline)

       ! bundleplot
    elseif (equal(word,'bundleplot')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before bundleplot',faterr,command,syntax=.true.)
          return
       end if
       call bundleplot(subline)

       ! sphereintegrals
    elseif (equal(word,'sphereintegrals')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before sphereintegrals',faterr,command,syntax=.true.)
          return
       end if
       call sphereintegrals(subline)

       ! integrals
    elseif (equal(word,'integrals')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before integrals',faterr,command,syntax=.true.)
          return
       end if
       call integrals(subline)

       ! qtree 
    elseif (equal(word,'qtree')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before qtree',faterr,command,syntax=.true.)
          return
       end if
       ok = eval_next(level,command,lp)
       if (.not.ok) level = 6
       ok = eval_next(plevel,command,lp)
       if (.not.ok) plevel = 0
       call check_no_extra_word(ok)
       if (.not.ok) return
       call qtree_integration(level,plevel)

       ! yt
    elseif (equal(word,'yt')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before yt',faterr,command,syntax=.true.)
          return
       end if
       call intgrid_driver(command)

       ! bader
    elseif (equal(word,'bader')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before bader',faterr,command,syntax=.true.)
          return
       end if
       call intgrid_driver(command)

       ! xdm
    elseif (equal(word,'xdm')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before xdm',faterr,command,syntax=.true.)
          return
       end if
       call xdm_driver(subline)

       ! stm
    elseif (equal(word,'stm')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before stm',faterr,command,syntax=.true.)
          return
       end if
       call stm_driver(command(lp:))

       ! sphfactor
    elseif (equal(word,'sphfactor')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before sphfactor',faterr,command,syntax=.true.)
          return
       end if
       call qtree_setsphfactor(subline)

       ! root
    elseif (equal(word,'root')) then
       if (len_trim(command(lp:)) < 1) then
          call ferror('critic2','need a string for root',faterr,command,syntax=.true.)
          return
       end if
       fileroot = command(lp:)

       ! hirshfeld
    elseif (equal(word,'hirshfeld')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before hirshfeld',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       call hirsh_props_grid()

       ! ewald
    elseif (equal(word,'ewald')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before ewald',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       if (cr%ismolecule) then
          call ferror("critic2","EWALD can not be used with molecules",faterr)
          return
       end if

       rdum = cr%ewald_energy()
       write (uout,'("* Ewald electrostatic energy (Hartree) = ",A/)') &
          string(rdum,'e',decimal=12)

       ! environ
    elseif (equal(word,'environ')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before environ',faterr,command,syntax=.true.)
          return
       end if
       call struct_environ(command(lp:))

       ! packing
    elseif (equal(word,'packing')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before packing',faterr,command,syntax=.true.)
          return
       end if
       call struct_packing(command(lp:))

       ! identify
    elseif (equal(word,'identify')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before identify',faterr,command,syntax=.true.)
          return
       end if
       call varbas_identify(command,lp)

       ! sum
    elseif (equal(word,'sum')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before sum',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','SUM: field is not allocated',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       write (uout,'("SUM(",A,") = ",A/)') string(id), string(sum(f(id)%f),'e',decimal=12)

       ! min
    elseif (equal(word,'min')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before min',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','MIN: field is not allocated',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       write (uout,'("MIN(",A,") = ",A/)') string(id), string(minval(f(id)%f),'e',decimal=12)

       ! max
    elseif (equal(word,'max')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before max',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','MAX: field is not allocated',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       write (uout,'("MAX(",A,") = ",A/)') string(id), string(maxval(f(id)%f),'e',decimal=12)

       ! mean
    elseif (equal(word,'mean')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before mean',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','MEAN: field is not allocated',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       write (uout,'("MEAN(",A,") = ",A/)') string(id),&
          string(sum(f(id)%f) / (f(id)%n(1)*f(id)%n(2)*f(id)%n(3)),'e',decimal=12)

       ! count
    elseif (equal(word,'count')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before count',faterr,command,syntax=.true.)
          return
       end if
       word = getword(command,lp)
       id = fieldname_to_idx(word)
       if (id < 0) id = refden
       if (.not.goodfield(id)) then
          call ferror('critic2','COUNT: field is not allocated',faterr,command,syntax=.true.)
          return
       end if

       ok = eval_next(rdum,command,lp)
       if (.not.ok) rdum = 0d0

       if (f(id)%type /= type_grid) then
          call ferror("count","count can only be used with grids",faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       write (uout,'("COUNT(",A," > ",A,") = ",A/)') string(id), &
          
          string(rdum,'e',decimal=7), string(count(f(id)%f > rdum))

       ! testrmt
    elseif (equal(word,'testrmt')) then
       if (.not. cr%isinit) then
          call ferror('critic2','need crystal before testrmt',faterr,command,syntax=.true.)
          return
       end if
       call check_no_extra_word(ok)
       if (.not.ok) return
       call testrmt(refden,3)

       ! clear
    elseif (equal(word,'clear')) then
       call critic_clearvariable(command(lp:))

       ! list
    elseif (equal(word,'list')) then
       call check_no_extra_word(ok)
       if (.not.ok) return
       call listvariables()
       call listfields()
       call listfieldalias()

       ! run/system
    elseif (equal(word,'run') .or. equal(word,'system')) then
       call system(command(lp:))

       ! echo
    elseif (equal(word,'echo')) then
       write (uout,'(A)') string(command(lp:))

       ! trick
    elseif (equal(word,'trick')) then
       call trick(command(lp:))

    else
       lp = 1
       call critic_setvariables(command, lp)
    endif
    
  contains
    !> Set field number id as reference
    subroutine set_reference(id)

      integer, intent(in) :: id

      ! header and change refden
      write (uout,'("* Field number ",A," is now REFERENCE."/)') string(id)

      ! initialize CP list, defer the calculation of nuclei properties to the report
      if (refden /= id .or. id == 0) then
         refden = id
         call init_cplist(.true.)
      end if

      ! define second integrable property as the valence charge.
      nprops = max(2,nprops)
      integ_prop(2)%used = .true.
      integ_prop(2)%itype = itype_fval
      integ_prop(2)%fid = id
      integ_prop(2)%prop_name = "Pop"

      ! define third integrable property as the valence laplacian.
      nprops = max(3,nprops)
      integ_prop(3)%used = .true.
      integ_prop(3)%itype = itype_lapval
      integ_prop(3)%fid = id
      integ_prop(3)%prop_name = "Lap"

      ! reset defaults for qtree
      if (f(refden)%type == type_grid) then
         gradient_mode = 1
         if (INT_radquad_errprop_default) INT_radquad_errprop = 2
      else
         gradient_mode = 2
         if (INT_radquad_errprop_default) INT_radquad_errprop = 3
      end if

      ! report
      call fields_integrable_report()

    end subroutine set_reference

    !> Cleans the structure and associated data
    subroutine clean_structure()
      ! ...the structural data
      call cr%end()
      call cr%init()
      ! ...the fields associated to the previous structure
      call fields_end()
      call fields_init()
      ! ...the loaded radial atomic and core densities
      call grda_end()
      ! ...the CP list
      call varbas_end()
    end subroutine clean_structure

    subroutine check_no_extra_word(ok)
      character(len=:), allocatable :: aux2
      logical, intent(out) :: ok
      aux2 = getword(command,lp)
      ok = .true.
      if (len_trim(aux2) > 0) then
         call ferror('critic','Unknown extra keyword',faterr,command,syntax=.true.)
         ok = .false.
      end if
    end subroutine check_no_extra_word

  end subroutine run_critic2_command

end module gui_interface
