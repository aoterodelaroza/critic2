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
  use struct_basic, only: crystal
  use c_interface_module
  implicit none

  private

  public :: critic2_initialize
  public :: critic2_end
  public :: open_structure
  public :: new_structure
  public :: preview_structure
  public :: accept_previewed_structure
  public :: reject_previewed_structure
  private :: initialize_structure
  public :: call_auto
  public :: update_scene
  public :: clear_scene
  public :: get_text_info
  private :: stick_from_endpoints
  private :: denewline
  private :: save_state
  private :: load_state
  private :: clear_state

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
    use struct_basic, only: cr
    use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
       fcflags, cc, cflags, ldflags, enable_debug, package
    use global, only: global_init, fileroot, config_write, initial_banner
    use tools_io, only: stdargs, ioinit, ucopy, uout, start_clock, &
       tictac
    use param, only: param_init
    character(len=:), allocatable :: optv
    character(len=:), allocatable :: ghome

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

    ! clear the scene
    call clear_scene(logical(.false.,c_bool))

  end subroutine critic2_initialize

  !> End of the critic2 run
  subroutine critic2_end() bind(c)
    use fields, only: fields_end
    use struct_basic, only: cr
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
    use struct, only: struct_crystal_input
    use struct_basic, only: cr

    type(c_ptr), intent(in) :: filename0
    integer(c_int), value :: ismolecule
    character(len=:), allocatable :: filename

    integer :: isformat

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
  function new_structure(useed,preview) bind(c)
    use struct_basic, only: crystalseed, spgs_wrap, cr
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
       call spgs_wrap(fseed,aux,.false.)
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
    use struct_basic, only: crystalseed, cr
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
    use struct_basic, only: cr
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
    use struct_basic, only: cr
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
    use struct_basic, only: cr
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
    use struct_basic, only: cr
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
    use struct_basic, only: cr, pointgroup_info, laue_string, holo_string
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
    
    character*1, parameter :: nl = new_line('a')
    integer :: i

    do i = 1, len_trim(str)
       if (str(i:i) == nl) str(i:i) = " "
    end do

  end subroutine denewline

  ! put the current scene into the save state
  subroutine save_state()
    use struct_basic, only: cr

    if (.not.isinit) return
    isinit_ = isinit
    ismolecule_ = ismolecule
    nat_ = nat
    if (associated(at_f_)) deallocate(at_f_)
    allocate(at_f_(size(at_f,1)))
    at_f_ = at_f
    nbond_ = nbond
    if (associated(bond_f_)) deallocate(bond_f_)
    allocate(bond_f_(size(bond_f,1)))
    bond_f_ = bond_f
    ncritp_ = ncritp
    if (associated(critp_f_)) deallocate(critp_f_)
    allocate(critp_f_(size(critp_f,1)))
    critp_f_ = critp_f
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
    use struct_basic, only: cr

    if (issaved) then
       isinit = isinit_
       ismolecule = ismolecule_
       nat = nat_
       if (associated(at_f)) deallocate(at_f)
       allocate(at_f(size(at_f_,1)))
       at_f = at_f_
       nbond = nbond_
       if (associated(bond_f)) deallocate(bond_f)
       allocate(bond_f(size(bond_f_,1)))
       bond_f = bond_f_
       ncritp = ncritp_
       if (associated(critp_f)) deallocate(critp_f)
       allocate(critp_f(size(critp_f_,1)))
       critp_f = critp_f_
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

end module gui_interface
