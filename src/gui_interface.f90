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

  ! C-interoperable atom type
  type, bind(c) :: c_atom
     character(kind=c_char,len=1) :: name(11) !< Atomic name
     integer(c_int) :: z !< Atomic number
     real(c_float) :: r(3) !< Atomic position (Cartesian, bohr)
  end type c_atom

  ! C-interoperable bond type
  type, bind(c) :: c_bond
     integer(c_int) :: i1 !< index for the first atom (C-style)
     integer(c_int) :: i2 !< index for the second atom (C-style)
     real(c_float) :: r1(3) !< position of the first atom (bohr)
     real(c_float) :: r2(3) !< position of the second atom (bohr)
     real(c_float) :: rmid(3) !< position of the bond center (bohr)
     real(c_float) :: length !< length of the bond (bohr)
     real(c_float) :: rot(4,4) !< rotation matrix
  end type c_bond

  ! C-interoperable critical point
  type, bind(c) :: c_critp
     real(c_float) :: r(3) !< position (bohr)
     integer(c_int) :: type !< type (-3,-1,1,3)
     character(kind=c_char,len=1) :: name(11) !< Name
  end type c_critp

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
  
  ! bounding box limits and center
  real(c_float), bind(c) :: xmin(3)
  real(c_float), bind(c) :: xmax(3)
  real(c_float), bind(c) :: xcm(3)
  real(c_float), bind(c) :: xmaxlen

  public :: get_positions
  public :: get_atom_position
  public :: get_num_atoms
  public :: num_of_bonds
  public :: get_atom_bond
  public :: num_of_crit_points
  public :: get_cp_pos_type

  ! scene blueprint (also, see escher's representation):
  ! 
  ! critical point
  ! position
  ! type
  ! typename
  ! selected
  !
  ! cell
  ! molecular cell

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
    call clear_scene()

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
  subroutine call_structure(filename0, nc, isMolecule) bind(c)
    use fields, only: nprops, integ_prop, f, type_grid, itype_fval, itype_lapval,&
       fields_integrable_report
    use grd_atomic, only: grda_init
    use struct, only: struct_crystal_input
    use struct_basic, only: cr
    use tools_io, only: uout, string
    use global, only: refden, gradient_mode, INT_radquad_errprop_default, INT_radquad_errprop
    use autocp, only: init_cplist

    character (kind=c_char, len=1), dimension (*), intent (in) :: filename0
    integer (kind=c_int), value :: nc
    integer (kind=c_int), value :: isMolecule

    character(len=:), allocatable :: filename

    filename = str_c_to_f(filename0, nc)

    call struct_crystal_input(cr, filename, isMolecule == 1, .true.)
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

    if (cr%isinit) then
       call autocritic("")
    end if

    call update_scene()

  end subroutine call_auto

  ! Update the data in the module variables - makes it available
  ! to the C++ code
  subroutine update_scene() bind(c)
    use varbas, only: ncpcel, cpcel
    use struct_basic, only: cr
    use tools_math, only: norm, cross
    
    integer :: i, j, npts
    real*8 :: x1(3), x2(3), xd(3), uz(3), normd, u(3), nu
    real*8 :: ca, sa, rot(4,4)

    ! Prepare for bounding box calculation
    xmin = 1e30
    xmax = -1e30
    xcm = 0.
    xmaxlen = 0.
    npts = 0

    ! Allocate space for atoms
    if (associated(at_f)) deallocate(at_f)
    allocate(at_f(max(cr%ncel,1)))

    ! For now, just go ahead and represent the whole cell/molecule
    nat = cr%ncel
    do i = 1, cr%ncel
       at_f(i)%z = cr%at(cr%atcel(i)%idx)%z
       x1 = cr%atcel(i)%r + cr%molx0
       at_f(i)%r = x1
       call f_c_string(cr%at(cr%atcel(i)%idx)%name,at_f(i)%name,11)
       xmin = min(x1,xmin)
       xmax = max(x1,xmax)
       xcm = xcm + x1
       npts = npts + 1
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
             x1 = cr%atcel(i)%r + cr%molx0
             x2 = cr%atcel(cr%nstar(i)%idcon(j))%r + cr%molx0
             bond_f(nbond)%r1 = x1
             bond_f(nbond)%r2 = x2
             bond_f(nbond)%rmid = 0.5d0 * (x1+x2)
             xd = x1 - x2
             normd = norm(xd)
             bond_f(nbond)%length = normd

             ! calculate the axis and angle to rotate the bond to the z-axis
             xd = xd / normd
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

             ! rotation matrix for this transformation
             rot(1,1) = ca + u(1)*u(1)*(1d0 - ca)
             rot(2,1) = u(1)*u(2)*(1d0 - ca) - u(3)*sa
             rot(3,1) = u(1)*u(3)*(1d0 - ca) + u(2)*sa
             rot(4,1) = 0d0

             rot(1,2) = u(1)*u(2)*(1 - ca) + u(3)*sa
             rot(2,2) = ca + u(2)*u(2)*(1d0 - ca)
             rot(3,2) = u(2)*u(3)*(1d0 - ca) - u(1)*sa
             rot(4,2) = 0d0

             rot(1,3) = u(1)*u(3)*(1d0 - ca) - u(2)*sa
             rot(2,3) = u(3)*u(2)*(1d0 - ca) + u(1)*sa
             rot(3,3) = ca + u(3)*u(3)*(1d0 - ca)
             rot(4,3) = 0d0

             rot(1,4) = 0d0
             rot(2,4) = 0d0
             rot(3,4) = 0d0
             rot(4,4) = 1d0

             ! save the rotation matrix 
             bond_f(nbond)%rot = rot
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
       x1 = cpcel(i)%r + cr%molx0
       critp_f(j)%r = x1
       critp_f(j)%type = cpcel(i)%typ
       call f_c_string(cpcel(i)%name,critp_f(j)%name,11)
       xmin = min(x1,xmin)
       xmax = max(x1,xmax)
       xcm = xcm + x1
       npts = npts + 1
    end do
    critp = c_loc(critp_f)

    ! wrap up center of the scene and dimensions
    xcm = xcm / npts
    xmaxlen = maxval(xmax - xmin)

  end subroutine update_scene

  ! Clear the scene
  subroutine clear_scene() bind(c)

    ! atoms
    if (associated(at_f)) deallocate(at_f)
    nat = 0
    at = C_NULL_PTR
    
    ! bounding box
    xmin = 0.
    xmax = 0.
    xcm = 0. 

  end subroutine clear_scene

  function str_c_to_f(strc, nchar) result(strf)
    integer(kind=C_INT), intent(in), value :: nchar
    character(kind=C_CHAR,len=1), intent(in) :: strc(nchar)
    character(len=:), allocatable :: strf

    integer :: i

    allocate(character(len=nchar) :: strf)
    do i = 1, nchar
      strf(i:i) = strc(i)
    end do

  endfunction str_c_to_f

  subroutine get_positions(n,z,x) bind(c,name="get_positions")
    use struct_basic, only: cr
    integer(c_int), intent(out) :: n
    type(c_ptr), intent(out) :: z
    type(c_ptr), intent(out) :: x
    integer(c_int), allocatable, target :: iz(:)
    real(c_double), allocatable, target :: ix(:,:)

    integer :: i, j

    n = cr%ncel

    allocate(iz(cr%ncel))
    do i = 1, cr%ncel
       iz(i) = int(cr%at(cr%atcel(i)%idx)%z,C_INT)
    end do

    allocate(ix(cr%ncel,3))
    do i = 1, cr%ncel
       do j = 1, 3
          ix(i,j) = real(cr%atcel(i)%r(j),c_double)
       end do
    end do

    x = c_loc(ix)
    z = c_loc(iz)

    deallocate(ix)
    deallocate(iz)

  end subroutine get_positions

  subroutine get_num_atoms(n) bind (c, name="get_num_atoms")
    use struct_basic, only: cr

    integer(c_int), intent(out) :: n

    n = cr%ncel

  end subroutine get_num_atoms

  subroutine get_atom_position(index, atomicN, x, y, z) bind (c, name="get_atom_position")
    use struct_basic, only: cr
    use types, only: celatom
    use global, only: dunit0, iunit

    integer (kind=c_int), value :: index
    integer(c_int), intent(out) :: atomicN
    real(c_double), intent(out) :: x
    real(c_double), intent(out) :: y
    real(c_double), intent(out) :: z

    atomicN = int(cr%at(cr%atcel(index)%idx)%z, c_int)

    x = real((cr%atcel(index)%r(1)+cr%molx0(1))*dunit0(iunit),c_double)
    y = real((cr%atcel(index)%r(2)+cr%molx0(2))*dunit0(iunit),c_double)
    z = real((cr%atcel(index)%r(3)+cr%molx0(3))*dunit0(iunit),c_double)

  end subroutine get_atom_position

  subroutine num_of_bonds(n_atom, nstarNum) bind (c, name="num_of_bonds")
    use struct_basic, only: cr
    integer (kind=c_int), value :: n_atom
    integer(c_int), intent(out) :: nstarNum

    call cr%find_asterisms()

    nstarNum = cr%nstar(n_atom)%ncon

  end subroutine num_of_bonds

  subroutine get_atom_bond(n_atom, nstarIdx, connected_atom, neighCrystal) bind (c, name="get_atom_bond")
    use struct_basic, only: cr
    integer (kind=c_int), value :: n_atom
    integer (kind=c_int), value :: nstarIdx
    integer(c_int), intent(out) :: connected_atom
    logical (kind=c_bool), intent(out) :: neighCrystal

    integer :: lconTrans(3)

    connected_atom = cr%nstar(n_atom)%idcon(nstarIdx)

    if (.NOT. cr%ismolecule) then
      lconTrans = cr%nstar(n_atom)%lcon(:, nstarIdx)

      if (lconTrans(1) /= 0 .OR. lconTrans(2) /= 0 .OR. lconTrans(3) /= 0) then
        neighCrystal = .true.
      end if
    end if

  end subroutine get_atom_bond

  subroutine num_of_crit_points(n_critp) bind (c, name="num_of_crit_points")
    use varbas, only: ncpcel
    integer(c_int), intent(out) :: n_critp

    n_critp = ncpcel

  end subroutine num_of_crit_points

  subroutine get_cp_pos_type(cpIdx, type, x, y, z) bind (c, name="get_cp_pos_type")
    use struct_basic, only: cr
    use global, only: dunit0, iunit
    use varbas, only: cpcel

    integer (kind=c_int), value :: cpIdx
    integer(c_int), intent(out) :: type
    real(c_double), intent(out) :: x
    real(c_double), intent(out) :: y
    real(c_double), intent(out) :: z

    x = (cpcel(cpIdx)%r(1) + cr%molx0(1))*dunit0(iunit)
    y = (cpcel(cpIdx)%r(2) + cr%molx0(2))*dunit0(iunit)
    z = (cpcel(cpIdx)%r(3) + cr%molx0(3))*dunit0(iunit)

    type = cpcel(cpIdx)%typ

  end subroutine get_cp_pos_type

end module gui_interface
