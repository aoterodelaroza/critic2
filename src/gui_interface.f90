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
  use c_interface_module
  implicit none

  private

  public :: critic2_initialize
  public :: critic2_end
  public :: call_structure
  public :: call_auto
  public :: update_scene
  public :: clear_scene
  public :: get_text_info

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

  ! global flags
  logical(c_bool), bind(c) :: isinit
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
  subroutine call_structure(filename0, ismolecule) bind(c)
    use fields, only: nprops, integ_prop, f, type_grid, itype_fval, itype_lapval,&
       fields_integrable_report
    use grd_atomic, only: grda_init
    use struct, only: struct_crystal_input
    use struct_basic, only: cr
    use tools_io, only: uout, string
    use global, only: refden, gradient_mode, INT_radquad_errprop_default, INT_radquad_errprop
    use autocp, only: init_cplist

    type(c_ptr), intent(in) :: filename0
    integer(c_int), value :: ismolecule
    character(len=:), allocatable :: filename

    integer :: isformat

    ! transform to fortran string
    filename = c_string_value(filename0)

    call struct_crystal_input(cr, filename, ismolecule, .true.)
    if (cr%isinit) then
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
    else
       call cr%init()
    end if

    ! update the scene
    call update_scene()

  end subroutine call_structure

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
    use varbas, only: ncpcel, cpcel
    use struct_basic, only: cr
    use tools_math, only: norm, cross
    use param, only: atmcov, jmlcol, maxzat
    
    integer :: i, j, iz
    real*8 :: x1(3), x2(3), xcm(3)

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
    do i = 1, ncpcel
       x1 = cpcel(i)%r + cr%molx0
       box_xmin = min(x1,box_xmin)
       box_xmax = max(x1,box_xmax)
       xcm = xcm + x1
    end do
    xcm = xcm / ncpcel
    box_xcm = xcm
    box_xmaxlen = maxval(box_xmax - box_xmin)

    ! Allocate space for atoms
    if (associated(at_f)) deallocate(at_f)
    allocate(at_f(max(cr%ncel,1)))

    ! For now, just go ahead and represent the whole cell/molecule
    nat = cr%ncel
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       at_f(i)%z = iz
       at_f(i)%b%r = cr%atcel(i)%r + cr%molx0 - xcm
       if (atmcov(iz) > 1) then
          at_f(i)%b%rad = atmcov(iz)
       else
          at_f(i)%b%rad = 2d0*atmcov(iz)
       end if
       at_f(i)%b%rgb = real(jmlcol(:,iz),4) / 255.
       call f_c_string(cr%at(cr%atcel(i)%idx)%name,at_f(i)%name,11)
    end do
    at = c_loc(at_f)

    ! Allocate space for bonds. For now, only bonds between atoms in
    ! the main cell.
    nbond = 0
    do i = 1, cr%ncel
       do j = 1, cr%nstar(i)%ncon
          if (cr%nstar(i)%idcon(j) > i .and. all(cr%nstar(i)%lcon(:,j) == 0)) then
             nbond = nbond + 1
          end if
       end do
    end do
    if (associated(bond_f)) deallocate(bond_f)
    allocate(bond_f(max(nbond,1)))
    nbond = 0
    do i = 1, cr%ncel
       do j = 1, cr%nstar(i)%ncon
          if (cr%nstar(i)%idcon(j) > i .and. all(cr%nstar(i)%lcon(:,j) == 0)) then
             nbond = nbond + 1
             bond_f(nbond)%i1 = i-1
             bond_f(nbond)%i2 = cr%nstar(i)%idcon(j)-1
             x1 = cr%atcel(i)%r + cr%molx0 - xcm
             x2 = cr%atcel(cr%nstar(i)%idcon(j))%r + cr%molx0 - xcm
             bond_f(nbond)%s = stick_from_endpoints(x1,x2,bondthickness,bondcolor)
          end if
       end do
    end do
    bond = c_loc(bond_f)

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

end module gui_interface
