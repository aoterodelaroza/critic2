! Copyright (c) 2015 Alberto Otero de la Roza
! <alberto@fluor.quimica.uniovi.es>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! The struct_read_qein and qe_latgen routines in this module were
! adapted from
! Quantum ESPRESSO, version 4.3.2.
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!> Read the crystal structure from files in several formats
module struct_readers
  implicit none
  private

  public :: parse_crystal_env
  public :: parse_molecule_env
  public :: struct_read_library
  public :: struct_read_cif
  public :: struct_read_cube
  public :: struct_read_wien
  public :: struct_read_vasp
  public :: struct_read_potcar
  public :: struct_read_abinit
  public :: struct_read_elk
  public :: struct_read_mol
  public :: struct_read_qeout
  public :: struct_read_qein
  public :: struct_read_siesta
  private :: qe_latgen
  private :: spgs_wrap
  private :: fill_molecule
  
contains

  !> Parse a crystal environment
  subroutine parse_crystal_env(c,lu)
    use struct_basic
    use global
    use arithmetic
    use tools_math
    use tools_io
    use types
    use param

    type(crystal), intent(inout) :: c !< Crystal
    integer, intent(in) :: lu !< Logical unit for input

    character(len=:), allocatable :: word, aux, aexp, line
    character*255, allocatable :: sline(:)
    integer :: i, j, k, lp, nsline, idx, luout, iat, lp2
    real*8 :: gmat(3,3), rmat(3,3), scal, ascal, x(3), xn(3)
    logical :: ok, goodcell, goodspg, useit

    character*(1), parameter :: ico(3) = (/"x","y","z"/)
    logical :: icodef(3), iok
    real*8 :: icoval(3)

    goodcell = .false.
    goodspg = .false.
    c%nneq = 0
    nsline = 0
    if (lu == uin) then
       luout = ucopy
    else
       luout = -1
    endif
    do while (getline(lu,line,ucopy=luout))
       lp = 1
       word = lgetword(line,lp)
       if (equal (word,'cell')) then
          ! cell <a> <b> <c> <alpha> <beta> <gamma>
          goodcell = .true.
          ok = eval_next(c%aa(1),line,lp)
          ok = ok .and. eval_next(c%aa(2),line,lp)
          ok = ok .and. eval_next(c%aa(3),line,lp)
          ok = ok .and. eval_next(c%bb(1),line,lp)
          ok = ok .and. eval_next(c%bb(2),line,lp)
          ok = ok .and. eval_next(c%bb(3),line,lp)
          if (.not.ok) then
             call ferror("parse_crystal_env","Wrong CELL syntax",faterr,line)
          endif
          word = lgetword(line,lp)
          if (equal(word,'angstrom') .or.equal(word,'ang')) then
             c%aa = c%aa / bohrtoa
          elseif (.not.equal(word,'bohr').and..not.equal(word,'au').and.len_trim(word) > 0) then
             call ferror('parse_crystal_input','Unknown extra keyword in CELL',faterr,line)
          endif
          call c%set_cryscar()

          ! cartesian <scale> .. endcartesian
       else if (equal (word,'cartesian')) then
          ok = eval_next(scal,line,lp)
          if (.not.ok) scal = 1d0
          ascal = 1d0/dunit
          aux = getword(line,lp)
          if (len_trim(aux) > 0) call ferror('parse_crystal_input','Unknown extra keyword in CARTESIAN',faterr,line)

          i = 0
          rmat = 0d0
          do while(.true.)
             lp = 1
             ok = getline(lu,line,ucopy=luout)
             word = lgetword(line,lp)
             if (equal(word,'angstrom') .or.equal(word,'ang')) then
                ! angstrom/ang
                ascal = 1d0/bohrtoa
             else if (equal(word,'bohr') .or.equal(word,'au')) then
                ! bohr/au
                ascal = 1d0
             else if (equal(word,'end').or.equal(word,'endcartesian')) then
                ! end/endcartesian
                aux = getword(line,lp)
                if (len_trim(aux) > 0) call ferror('parse_crystal_input','Unknown extra keyword in CARTESIAN',faterr,line)
                exit
             else
                ! matrix row
                i = i + 1
                if (i > 3) then
                   ok = .false.
                else
                   lp = 1
                   ok = ok .and. eval_next(rmat(i,1),line,lp)
                   ok = ok .and. eval_next(rmat(i,2),line,lp)
                   ok = ok .and. eval_next(rmat(i,3),line,lp)
                end if
             end if
             if (.not.ok) call ferror('parse_crystal_env','Bad CARTESIAN environment',faterr,line)
             aux = getword(line,lp)
             if (len_trim(aux) > 0) call ferror('parse_crystal_input','Unknown extra keyword in CARTESIAN',faterr,line)
          end do
          gmat = matmul(rmat,transpose(rmat)) * scal**2 * ascal**2
          do i = 1, 3
             c%aa(i) = sqrt(gmat(i,i))
          end do
          c%bb(1) = acos(gmat(2,3)/c%aa(2)/c%aa(3)) * 180d0 / pi
          c%bb(2) = acos(gmat(1,3)/c%aa(1)/c%aa(3)) * 180d0 / pi
          c%bb(3) = acos(gmat(1,2)/c%aa(1)/c%aa(2)) * 180d0 / pi
          c%crys2car = transpose(rmat) * scal * ascal
          c%car2crys = matinv(c%crys2car)
          goodcell = .true.

       else if (equal(word,'spg').or.equal(word,'spgr')) then
          ! spg <spg>
          useit = equal(word,'spgr')
          word = line(lp:)
          call spgs_wrap(c,word,useit)
          goodspg = .true.

       else if (equal(word,'symm')) then
          ! symm <line>
          if (.not.allocated(sline)) allocate(sline(10))
          nsline = nsline + 1
          if (nsline > size(sline)) call realloc(sline,2*size(sline))
          sline(nsline) = lower(adjustl(line(lp:)))
          
       else if (equal(word,'endcrystal') .or. equal(word,'end')) then
          ! endcrystal/end
          exit
       else
          ! keyword not found, must be an atom. The syntax:
          !    neq <x> <y> <z> <atom> ...
          !    <atom> <x> <y> <z> ...
          !    <atnumber> <x> <y> <z> ...
          ! are acceptable
          if (.not.equal(word,'neq')) then
             c%nneq = c%nneq+1
             if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))

             ! try to read four fields from the input
             lp2 = 1
             ok = isinteger(iat,line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(1),line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp2)
             if (.not.ok) then
                ! then it must be <atom> <x> <y> <z>
                ok = eval_next(c%at(c%nneq)%x(1),line,lp)
                ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp)
                ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp)
                if (.not.ok) &
                   call ferror("parse_crystal_env","Wrong atomic input syntax",faterr,line)
                c%at(c%nneq)%name = string(word)
             else
                lp = lp2
                c%at(c%nneq)%name = nameguess(iat,.true.)
             end if
          else
             c%nneq = c%nneq+1
             if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))
             ok = eval_next(c%at(c%nneq)%x(1),line,lp)
             ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp)
             ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp)
             if (.not.ok) &
                call ferror("parse_crystal_env","Wrong NEQ syntax",faterr,line)
             c%at(c%nneq)%name = getword(line,lp)
             c%at(c%nneq)%name = string(c%at(c%nneq)%name)
          end if

          c%at(c%nneq)%z = zatguess(c%at(c%nneq)%name)
          if (c%at(c%nneq)%z < 0) &
             call ferror('parse_crystal_input','Unknown atomic symbol in NEQ',faterr,line)
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'zpsp')) then
                ok = eval_next(c%at(c%nneq)%zpsp,line,lp)
                if (.not.ok) call ferror('parse_crystal_env','Wrong ZPSP in neq',faterr,line)
             else if (equal(word,'q')) then
                ok = eval_next(c%at(c%nneq)%qat,line,lp)
                if (.not.ok) call ferror('parse_crystal_env','Wrong Q in neq',faterr,line)
             else if (equal(word,'ang') .or. equal(word,'angstrom')) then
                if (.not.goodcell) &
                   call ferror('struct_parse_input','Need cell parameters before cartesian neq',faterr,line)
                c%at(c%nneq)%x = c%c2x(c%at(c%nneq)%x / bohrtoa)
             else if (equal(word,'bohr') .or. equal(word,'au')) then
                if (.not.goodcell) &
                   call ferror('struct_parse_input','Need cell parameters before cartesian neq',faterr,line)
                c%at(c%nneq)%x = c%c2x(c%at(c%nneq)%x)
             else if (len_trim(word) > 0) then
                call ferror('parse_crystal_input','Unknown keyword in NEQ',faterr,line)
             else
                exit
             end if
          end do
       end if
    end do
    aux = getword(line,lp)
    if (len_trim(aux) > 0) call ferror('parse_crystal_input','Unknown extra keyword in ENDCRYSTAL',faterr,line)
    if (c%nneq == 0) call ferror('parse_crystal_input','No atoms in input',faterr)

    ! symm transformation
    if (nsline > 0 .and. allocated(sline)) then
       ! save the old x,y,z variables if they are defined
       do k = 1, 3
          icodef(k) = isvariable(ico(k),icoval(k))
       end do

       do i = 1, nsline ! run over symm lines
          do j = 1, c%nneq ! run over atoms
             line = trim(adjustl(sline(i)))
             xn = c%at(j)%x - floor(c%at(j)%x)
             ! push the atom coordinates into x,y,z variables
             do k = 1, 3 
                call setvariable(ico(k),c%at(j)%x(k))
             end do

             ! parse the three fields in the arithmetic expression 
             do k = 1, 3
                if (k < 3) then
                   idx = index(line,",")
                   if (idx == 0) &
                      call ferror('parse_crystal_env','error reading symmetry operation',faterr,line)
                   aexp = line(1:idx-1)
                   aux = adjustl(line(idx+1:))
                   line = aux
                else
                   aexp = line
                end if
                x(k) = eval(aexp,.true.,iok)
             end do
             x = x - floor(x)

             ! check if this atom already exists
             ok = .true.
             do k = 1, c%nneq
                if (all(abs(x - xn) < 1d-5)) then
                   ok = .false.
                   exit
                endif
             end do
             
             ! add this atom to the list
             if (ok) then
                c%nneq = c%nneq+1
                if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))
                c%at(c%nneq) = c%at(j)
                c%at(c%nneq)%x = x
             endif
          end do
       end do
       deallocate(sline)

       ! re-set the previous values of x, y, z
       do k = 1, 3
          if (icodef(k)) call setvariable(ico(k),icoval(k))
       end do
    end if

    if (.not.goodspg) call c%guessspg(doguess,.false.)

  end subroutine parse_crystal_env

  !> Parse a molecule environment
  subroutine parse_molecule_env(c,lu)
    use struct_basic
    use global
    use arithmetic
    use tools_math
    use tools_io
    use types
    use param

    type(crystal), intent(inout) :: c !< Crystal
    integer, intent(in) :: lu !< Logical unit for input

    character(len=:), allocatable :: word, aux, line
    integer :: i, j, k, lp, lp2, luout, iat
    real*8 :: rborder
    logical :: ok, docube

    docube = .false.
    rborder = rborder_def 
    c%nneq = 0
    if (lu == uin) then
       luout = ucopy
    else
       luout = -1
    endif
    do while (getline(lu,line,ucopy=luout))
       lp = 1
       word = lgetword(line,lp)

       if (equal(word,'cube').or.equal(word,'cubic')) then
          ! cube
          docube = .true.
          word = lgetword(line,lp)
          call check_no_extra_word()

       else if (equal(word,'border')) then
          ! border [border.r]
          ok = eval_next(rborder,line,lp)
          if (.not.ok) &
             call ferror('parse_molecule_input','Wrong syntax in BORDER',faterr,line)
          call check_no_extra_word()

       else if (equal(word,'endmolecule') .or. equal(word,'end')) then
          ! endmolecule/end
          exit
       else
          ! keyword not found, must be an atom. The syntax:
          !    neq <x> <y> <z> <atom> ...
          !    <atom> <x> <y> <z> ...
          !    <atnumber> <x> <y> <z> ...
          ! are acceptable
          if (.not.equal(word,'neq')) then
             c%nneq = c%nneq+1
             if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))

             ! try to read four fields from the input
             lp2 = 1
             ok = isinteger(iat,line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(1),line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp2)
             ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp2)
             if (.not.ok) then
                ! then it must be <atom> <x> <y> <z>
                ok = eval_next(c%at(c%nneq)%x(1),line,lp)
                ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp)
                ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp)
                if (.not.ok) &
                   call ferror("parse_molecule_env","Wrong atomic input syntax",faterr,line)
                c%at(c%nneq)%name = string(word)
             else
                lp = lp2
                c%at(c%nneq)%name = nameguess(iat,.true.)
             endif
          else
             c%nneq = c%nneq+1
             if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))
             ok = eval_next(c%at(c%nneq)%x(1),line,lp)
             ok = ok .and. eval_next(c%at(c%nneq)%x(2),line,lp)
             ok = ok .and. eval_next(c%at(c%nneq)%x(3),line,lp)
             if (.not.ok) &
                call ferror("parse_molecule_env","Wrong NEQ syntax",faterr,line)
             c%at(c%nneq)%name = getword(line,lp)
             c%at(c%nneq)%name = string(c%at(c%nneq)%name)
          endif

          c%at(c%nneq)%x = c%at(c%nneq)%x / dunit
          c%at(c%nneq)%z = zatguess(c%at(c%nneq)%name)
          if (c%at(c%nneq)%z < 0) &
             call ferror('parse_molecule_input','Unknown atomic symbol or incorrect syntax',faterr,line)
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'zpsp')) then
                ok = eval_next(c%at(c%nneq)%zpsp,line,lp)
                if (.not.ok) call ferror('parse_molecule_env','Wrong ZPSP',faterr,line)
             else if (equal(word,'q')) then
                ok = eval_next(c%at(c%nneq)%qat,line,lp)
                if (.not.ok) call ferror('parse_molecule_env','Wrong Q',faterr,line)
             else if (equal(word,'ang') .or. equal(word,'angstrom')) then
                c%at(c%nneq)%x = c%at(c%nneq)%x / bohrtoa
             else if (equal(word,'bohr') .or. equal(word,'au')) then
                c%at(c%nneq)%x = c%at(c%nneq)%x
             else if (len_trim(word) > 0) then
                call ferror('parse_molecule_input','Unknown extra keyword in atomic input',faterr,line)
             else
                exit
             end if
          end do
       endif
    end do
    aux = getword(line,lp)
    if (len_trim(aux) > 0) call ferror('parse_molecule_input','Unknown extra keyword in ENDMOLECULE',faterr,line)
    if (c%nneq == 0) call ferror('parse_molecule_input','No atoms in input',faterr)

    ! fill the missing information
    call fill_molecule(c,rborder,docube)
    call c%set_cryscar()
    call c%guessspg(doguess,.false.)

  contains
    subroutine check_no_extra_word()
      character(len=:), allocatable :: aux2
      aux2 = getword(line,lp)
      if (len_trim(aux2) > 0) &
         call ferror('critic','Unknown extra keyword',faterr,line)
    end subroutine check_no_extra_word
  
  end subroutine parse_molecule_env

  !> Read a structure from the critic2 structure library
  subroutine struct_read_library(c,line,mol)
    use global
    use tools_io
    use struct_basic

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: line !< Library entry
    logical, intent(in) :: mol !< Is this a molecule?

    character(len=:), allocatable :: word, l2, stru, aux, libfile
    logical :: lchk, found, ok
    integer :: lu, lp, lpo

    ! read the structure
    lpo = 1
    stru = lgetword(line,lpo)
    if (len_trim(stru) < 1) &
       call ferror("struct_read_library","structure label missing in CRYSTAL/MOLECULE LIBRARY",faterr,line)

    if (mol) then
       libfile = mlib_file
    else
       libfile = clib_file
    endif

    ! open the library file
    inquire(file=libfile,exist=lchk)
    if (.not.lchk) then
       write (uout,'("(!) Library file:"/8X,A)') trim(libfile)
       call ferror("struct_read_library","library file not found!",faterr)
    endif
    lu = fopen_read(libfile,abspath=.true.)

    ! find the block
    found = .false.
    main: do while (getline(lu,l2))
       lp = 1
       word = lgetword(l2,lp)
       if (equal(word,'structure')) then
          do while(len_trim(word) > 0)
             word = lgetword(l2,lp)
             if (equal(word,stru)) then
                found = .true.
                exit main
             endif
          end do
       endif
    end do main
    if (.not.found) then
       write (uout,'("(!) Structure not found in file:"/8X,A)') trim(libfile)
       call ferror("struct_read_library","structure not found in library!",faterr)
    end if

    ! read the crystal/molecule environment inside
    ok = getline(lu,l2)
    if (mol) then
       call parse_molecule_env(c,lu)
    else
       call parse_crystal_env(c,lu)
    endif
    call fclose(lu)

    ! make sure there's no more input
    aux = getword(line,lpo)
    if (len_trim(aux) > 0) call ferror('struct_read_library','Unknown extra keyword in CRYSTAL/MOLECULE LIBRARY',faterr,line)

  end subroutine struct_read_library

  !> Read the structure from a CIF file (uses ciftbx)
  subroutine struct_read_cif(c,file,dblock,verbose,mol)
    use struct_basic
    use arithmetic
    use global
    use tools_io
    use types
    use param

    include 'ciftbx/ciftbx.cmv'
    include 'ciftbx/ciftbx.cmf'

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: file !< Input file name
    character*(*), intent(in) :: dblock !< Data block
    logical, intent(in) :: mol !< Is this a molecule? 
    logical, intent(in) :: verbose !< Write information to the output

    character(len=1024) :: dictfile, sym, tok
    character*30 :: atname, spg
    real*8 :: x(3)
    real*8 :: sigx, rot0(3,4), xo, yo, zo
    logical :: fl, fl1, fl2, found, ok, ix, iy, iz, iok
    integer :: i, j, ludum, luscr, idx

    character*(1), parameter :: ico(3) = (/"x","y","z"/)

    ludum = falloc()
    luscr = falloc()
    fl = init_(ludum, uout, luscr, uout)

    ! header
    if (verbose) then
       write (uout,'("* CIFtbx initialized")')
       write (uout,'("  ... by Sydney R. Hall and Herbert J. Bernstein (see src/ciftbx/ for details).")')
       write (uout,'("  ",A)') trim(adjustl(tbxver_))
       write (uout,*)
    end if

    ! open dictionary
    dictfile = trim(adjustl(critic_home)) // dirsep // 'cif_core.dic'
    if (verbose) &
       write (uout,'("+ Using dictionary: ",A/)') trim(adjustl(dictfile))
    fl = dict_(dictfile,'valid')
    if (.not.fl) &
       call ferror('struct_read_cif','Dictionary file (cif_core.dic) not found. Check CRITIC_HOME',faterr)

    ! open cif file
    fl = ocif_(file)
    if (.not.fl) &
       call ferror('struct_read_cif','CIF file not found',faterr,file)

    ! read data block
    if (verbose) &
       write (uout,'("+ Reading first data block")')
    fl = data_(dblock)
    if (.not.fl) &
       call ferror('struct_read_cif','incorrect named data block',faterr,file)
    if (verbose) &
       write (uout,*)

    ! read cell dimensions
    fl = numd_('_cell_length_a',c%aa(1),sigx)
    fl = fl .and. numd_('_cell_length_b',c%aa(2),sigx)
    fl = fl .and. numd_('_cell_length_c',c%aa(3),sigx)
    if (.not.fl) &
       call ferror('struct_read_cif','error reading cell lengths',faterr,file)
    c%aa = c%aa / bohrtoa

    ! read cell angles
    fl = numd_('_cell_angle_alpha',c%bb(1),sigx)
    fl = fl .and. numd_('_cell_angle_beta',c%bb(2),sigx)
    fl = fl .and. numd_('_cell_angle_gamma',c%bb(3),sigx)
    if (.not.fl) &
       call ferror('struct_read_cif','error reading cell angles',faterr,file)

    ! read atomic positions
    c%nneq = 1
    do while(.true.)
       if (c%nneq > size(c%at)) call realloc(c%at,2*c%nneq)
       fl = char_('_atom_site_label',atname)
       if (.not.fl) &
          fl = char_('_atom_site_type_symbol',atname)
       c%at(c%nneq)%z = zatguess(atname)
       if (c%at(c%nneq)%z < 0) &
          call ferror('struct_read_cif','Unrecognized atomic symbol',faterr,file)
       c%at(c%nneq)%name = atname(1:10)

       fl = fl .and. numd_('_atom_site_fract_x',x(1),sigx)
       fl = fl .and. numd_('_atom_site_fract_y',x(2),sigx)
       fl = fl .and. numd_('_atom_site_fract_z',x(3),sigx)
       c%at(c%nneq)%x = x
       if (.not.fl) &
          call ferror('struct_read_cif','error reading atomic positions',faterr,file)
       if (.not.loop_) exit
       c%nneq = c%nneq + 1
    end do

    ! save the old value of x, y, and z variables
    ix = isvariable("x",xo)
    iy = isvariable("y",yo)
    iz = isvariable("z",zo)

    ! use the symmetry information from _symmetry_equiv_pos_as_xyz
    found = .false.
    fl1 = .false.
    fl2 = .false.
    c%neqv = 0
    c%ncv = 1
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%cen(:,1) = 0d0
    do while(.true.)
       if (.not.found) then
          fl1 = char_('_symmetry_equiv_pos_as_xyz',sym)
          if (.not.fl1) fl2 = char_('_space_group_symop_operation_xyz',sym)
          if (.not.(fl1.or.fl2)) exit
          found = .true.
       else
          if (fl1) fl1 = char_('_symmetry_equiv_pos_as_xyz',sym)
          if (fl2) fl2 = char_('_space_group_symop_operation_xyz',sym)
       endif

       ! do stuff with sym
       if (.not.(fl1.or.fl2)) &
          call ferror('struct_read_cif','error reading symmetry xyz elements',faterr,file)

       ! process the three symmetry elements
       rot0 = 0d0
       sym = trim(adjustl(sym)) // ","
       do i = 1, 3
          ! extract the next token
          idx = index(sym,",")
          if (idx == 0) &
             call ferror('struct_read_cif','error reading symmetry operation',faterr,sym)
          tok = sym(1:idx-1)
          sym = sym(idx+1:)

          ! the translation component
          do j = 1, 3
             call setvariable(ico(j),0d0)
          end do
          rot0(i,4) = eval(tok,.true.,iok)

          ! the x-, y-, z- components
          do j = 1, 3
             call setvariable(ico(j),1d0)
             rot0(i,j) = eval(tok,.true.,iok) - rot0(i,4)
             call setvariable(ico(j),0d0)
          enddo
       enddo

       ! now we have a rot0
       if (all(abs(eyet - rot0) < 1d-12)) then
          ! the identity
          c%neqv = c%neqv + 1
          c%rotm(:,:,c%neqv) = rot0
       elseif (all(abs(eye - rot0(1:3,1:3)) < 1d-12)) then
          ! a non-zero pure translation
          ! check if I have it already
          ok = .true.
          do i = 1, c%ncv
             if (all(abs(rot0(:,4) - c%cen(:,i)) < 1d-12)) then
                ok = .false.
                exit
             endif
          end do
          if (ok) then
             c%ncv = c%ncv + 1
             if (c%ncv > size(c%cen,2)) call realloc(c%cen,3,2*c%ncv)
             c%cen(:,c%ncv) = rot0(:,4)
          endif
       else
          ! a rotation, with some pure translation in it
          ! check if I have this rotation matrix already
          ok = .true.
          do i = 1, c%neqv
             if (all(abs(c%rotm(1:3,1:3,i) - rot0(1:3,1:3)) < 1d-12)) then
                ok = .false.
                exit
             endif
          end do
          if (ok) then
             c%neqv = c%neqv + 1
             c%rotm(:,:,c%neqv) = rot0
          endif
       endif
       ! exit the loop
       if (.not.loop_) exit
    end do
    call realloc(c%cen,3,c%ncv)

    ! let struct_fill handle this
    c%lcent = 0

    ! restore the old values of x, y, and z
    if (ix) call setvariable("x",xo)
    if (iy) call setvariable("y",yo)
    if (iz) call setvariable("z",zo)

    ! read and process spg information
    if (.not.found) then
       ! the "official" Hermann-Mauginn symbol from the dictionary: many cif files don't have one
       fl = char_('_symmetry_space_group_name_H-M',spg)

       ! the "alternative" symbol... the core dictionary says I shouldn't be using this
       if (.not.fl) fl = char_('_space_group_name_H-M_alt',spg)

       ! oh, well, that's that...
       if (.not.fl) &
          call ferror('struct_read_cif','error reading symmetry',faterr,file)

       ! call spgs and hope for the best
       call spgs_wrap(c,spg,.false.)
    endif
    c%havesym = 2

    ! clean up
    call purge_()
    call fdealloc(ludum)
    call fdealloc(luscr)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_cif

  !> Read the structure from a gaussian cube file
  subroutine struct_read_cube(c,file,verbose,mol)
    use struct_basic
    use global
    use types
    use tools_io
    use tools_math
    use param

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: verbose !< Verbose?
    logical, intent(in) :: mol !< Is this a molecule?

    integer :: lu
    integer :: i, nstep(3), nn, iz
    real*8 :: x0(3), rmat(3,3), g(3,3), rdum, rx(3)

    lu = fopen_read(file)

    ! ignore the title lines
    read (lu,*)
    read (lu,*)

    ! number of atoms and unit cell
    read (lu,*) c%nneq, x0
    if (verbose) then
       write (uout,'("* Structure read from cube file: ",A)') trim(file)
       write (uout,'("  Number of atoms: ",I4)') c%nneq
       write (uout,'("  Cube origin (x0): ",3(F13.6,X))') x0
    endif

    do i = 1, 3
       read (lu,*) nstep(i), rmat(:,i)
       if (verbose) then
          write (uout,'("  Axis (nstep,dir): ",I5,3(F13.6,X))') nstep(i), rmat(:,i)
       endif
       rmat(:,i) = rmat(:,i) * nstep(i)
    end do

    g = matmul(transpose(rmat),rmat)
    c%crys2car = rmat
    c%car2crys = matinv(rmat)
    do i = 1, 3
       c%aa(i) = sqrt(g(i,i))
    end do
    c%bb(1) = acos(g(2,3) / c%aa(2) / c%aa(3)) * 180d0 / pi
    c%bb(2) = acos(g(1,3) / c%aa(1) / c%aa(3)) * 180d0 / pi
    c%bb(3) = acos(g(1,2) / c%aa(1) / c%aa(2)) * 180d0 / pi
    rmat = matinv(transpose(rmat))

    ! Atomic positions.
    if (c%nneq > size(c%at)) call realloc(c%at,c%nneq)
    nn = c%nneq
    c%nneq = 0
    do i = 1, nn
       read (lu,*) iz, rdum, rx
       if (iz > 0) then
          c%nneq = c%nneq + 1
          c%at(c%nneq)%z = iz
          rx = matmul(rx - x0,rmat)
          c%at(c%nneq)%x = rx - floor(rx)
       endif
    end do
    if (c%nneq /= nn) call realloc(c%at,c%nneq)

    ! no symmetry
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,1:3,1) = eye
    c%ncv = 0
    c%lcent = 1

    ! initialize atoms
    c%at(1:c%nneq)%zpsp = -1
    c%at(1:c%nneq)%qat = 0

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c,x0)

    call fclose(lu)

  end subroutine struct_read_cube

  !> Read the crystal structure from a WIEN2k STRUCT file.
  !> Code adapted from the WIEN2k distribution.
  subroutine struct_read_wien(c,file,readall,mol)
    use struct_basic
    use tools_io
    use param
    use types

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: file !< struct file
    logical, intent(in) :: readall !< if true, read the atomic positions into c%at%x
    logical, intent(in) :: mol !< is this a molecule?

    integer :: lut
    integer :: i, j, i1, i2, j1, iat, istart
    real*8 :: mat(3,3), rnot, rmt, pos(3), tau(3), znuc
    integer :: multw, iatnr, iz(3,3), jatom, mu, jri
    character*4 :: lattic, cform
    character*80 :: titel
    character*10 :: aname

    lut = fopen_read(file)

    READ(lut,102) TITEL
    READ(lut,103) LATTIC, c%nneq, cform
102 FORMAT(A80)
103 FORMAT(A4,23X,I3,1x,a4,/,4X,4X) ! new

    c%ncv = 0
    IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
       c%lcent=1
    ELSE IF(LATTIC(1:1).EQ.'F') THEN
       c%lcent=6
    ELSE IF(LATTIC(1:1).EQ.'B') THEN
       c%lcent=5
    ELSE IF(LATTIC(1:1).EQ.'H') THEN
       c%lcent=1
    ELSE IF(LATTIC(1:1).EQ.'R') THEN
       c%lcent=7
    ELSE IF(LATTIC(1:3).EQ.'CXY') THEN
       c%lcent=4
    ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN
       c%lcent=2
    ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN
       c%lcent=3
    ELSE
       STOP 'LATTIC NOT DEFINED'
    END IF

    READ(lut,100) c%aa(1:3), c%bb(1:3)
100 FORMAT(6F10.5)
    if(c%bb(3) == 0.d0) c%bb(3)=90.d0

    iat = 0
    DO JATOM=1,c%nneq
       iat = iat + 1
       if (iat > size(c%at)) call realloc(c%at,2*size(c%at))
       READ(lut,1012) iatnr,c%at(iat)%x(1:3),MULTW

       istart = iat
       if (readall) then
          DO MU=1,MULTW-1
             iat = iat + 1
             if (iat > size(c%at)) call realloc(c%at,2*size(c%at))
             READ(lut,1013) iatnr, c%at(iat)%x(1:3)
          end DO
       else
          DO MU=1,MULTW-1
             READ(lut,1013) iatnr, pos
          end DO
       end if

       READ(lut,113) ANAME,JRI,RNOT,RMT,Znuc
       c%at(istart:iat)%z = nint(Znuc)
       if (any(c%at(istart:iat)%z <= 0) .or. any(c%at(istart:iat)%z > 118)) then
          call ferror('struct_read_wien','Unrecognized type of atom',faterr,file)
       end if

       aname = adjustl(aname)
       c%at(istart:iat)%name = aname(1:10)

       READ(lut,1051) ((mat(I1,J1),I1=1,3),J1=1,3)
    end DO
113 FORMAT(A10,5X,I5,5X,F10.5,5X,F10.5,5X,F5.2)
1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/15X,I2) ! new
1013 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7) ! new
1051 FORMAT(20X,3F10.8)
    c%nneq = iat

    !.read number of symmetry operations, sym. operations
    READ(lut,114) c%neqv
114 FORMAT(I4)

    do i=1,c%neqv
       read(lut,115) ((iz(i1,i2),i1=1,3),tau(i2),i2=1,3)
       do j=1,3
          c%rotm(:,j,i)=dble(iz(j,:))
       enddo
       c%rotm(:,4,i)=tau
    end do
    if (c%neqv > 0) c%havesym = 2

115 FORMAT(3(3I2,F10.5,/))

    ! initialize charge and pseudopotential charge
    c%at(1:c%nneq)%zpsp = -1
    c%at(1:c%nneq)%qat = 0
    call realloc(c%at,c%nneq)

    ! clean up
    call fclose(lut)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_wien

  !> Read everything except the grid from a VASP POSCAR, etc. file
  subroutine struct_read_vasp(c,filename,ntypat,ztypat,mol)
    use struct_basic
    use types
    use tools_io
    use tools_math
    use param

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: filename !< Input file name
    integer, intent(inout) :: ntypat !< Number of atom types
    character*5, intent(inout) :: ztypat(100) !< Atomic numbers for the types
    logical, intent(in) :: mol !< Is this a molecule?

    integer :: lu, lp, nn, typ
    character(len=:), allocatable :: word, line
    logical :: ok, iscar

    integer :: i, j
    real*8 :: scalex, scaley, scalez, scale
    real*8 :: rprim(3,3), gprim(3,3)
    real*8 :: omegaa

    ! open
    lu = fopen_read(filename)

    ! read the title and the scale line
    ok = getline_raw(lu,line,.true.)
    ok = getline_raw(lu,line,.true.)
    lp = 1
    ok = isreal(scalex,line,lp)
    ok = ok .and. isreal(scaley,line,lp)
    if (.not.ok) then
       scale = scalex
       scalex = 1d0
       scaley = 1d0
       scalez = 1d0
    else
       ok = isreal(scalez,line,lp)
       scale = 1d0
    end if

    ! read the cell vectors and calculate the metric tensor
    do i = 1, 3
       read (lu,*) rprim(1,i), rprim(2,i), rprim(3,i)
    end do
    if (scale < 0d0) then
       gprim = matmul(transpose(rprim),rprim)
       omegaa = sqrt(detsym(gprim))
       ! adjust the lengths to give the volume
       scale = abs(scale) / abs(omegaa)**(1d0/3d0)
    end if
    rprim(1,:) = rprim(1,:) * scalex * scale
    rprim(2,:) = rprim(2,:) * scaley * scale
    rprim(3,:) = rprim(3,:) * scalez * scale
    rprim = rprim / bohrtoa
    gprim = matmul(transpose(rprim),rprim)
    omegaa = sqrt(detsym(gprim))
    if (omegaa < 0d0) call ferror('readvasp','negative cell volume',faterr,filename)
    c%crys2car = rprim
    c%car2crys = matinv(rprim)

    ! For versions >= 5.2, a line indicating the atom types appears here
    ok = getline_raw(lu,line,.true.)
    lp = 1
    word = getword(line,lp)
    if (zatguess(word) >= 0) then
       ! An atom name has been read -> read the rest of the line
       ntypat = 0
       do while (zatguess(word) >= 0)
          ntypat = ntypat + 1
          ztypat(ntypat) = string(word) // " "
          word = getword(line,lp)
       end do
       ok = getline_raw(lu,line,.true.)
    else
       if (ztypat(1) == "") then
          call ferror('struct_read_vasp','Atom types are required for VASP < 5.2 inputs',faterr,filename)
       end if
    end if

    ! read number of atoms of each type
    lp = 1
    c%nneq = 0
    do i = 1, ntypat
       ok = isinteger(nn,line,lp)
       if (.not.ok) call ferror('struct_read_vasp','Too many atom types in CRYSTAL',faterr,line)

       typ = zatguess(ztypat(i))
       do j = c%nneq+1, c%nneq+nn
          if (j > size(c%at)) call realloc(c%at,2*size(c%at))
          c%at(j)%z = typ
          c%at(j)%name = nameguess(typ)
       end do
       c%nneq = c%nneq + nn
    end do

    ! Read atomic positions (cryst. coords.)
    read(lu,*) line
    line = adjustl(line)
    if (line(1:1) == 's' .or. line(1:1) == 'S') then
       read(lu,*) line
       line = adjustl(line)
    endif
    iscar = .false.
    if (line(1:1) == 'd' .or. line(1:1) == 'D') then
       iscar = .false.
    elseif (line(1:1) == 'c' .or. line(1:1) == 'C' .or. line(1:1) == 'k' .or. line(1:1) == 'K') then
       iscar = .true.
    endif
    do i=1,c%nneq
       read(lu,*) c%at(i)%x
       if (iscar) &
          c%at(i)%x = matmul(c%car2crys,c%at(i)%x / bohrtoa)
    enddo

    call fclose(lu)

    ! fill in the required values for critic
    do i = 1, 3
       c%aa(i) = sqrt(gprim(i,i))
    end do
    c%bb(1)=acos(gprim(2,3)/sqrt(gprim(2,2)*gprim(3,3)))/pi*180d0
    c%bb(2)=acos(gprim(1,3)/sqrt(gprim(1,1)*gprim(3,3)))/pi*180d0
    c%bb(3)=acos(gprim(1,2)/sqrt(gprim(1,1)*gprim(2,2)))/pi*180d0

    ! initialize charge and pseudopotential charge
    c%at(1:c%nneq)%zpsp = -1
    c%at(1:c%nneq)%qat = 0

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_vasp

  !> Read everything except the grid from a VASP POSCAR, etc. file
  subroutine struct_read_potcar(filename,ntyp,ztyp)
    use tools_io
    use param
    use types

    character*(*), intent(in) :: filename !< Input file name
    integer, intent(out) :: ntyp !< Number of atom types
    character*5, intent(out) :: ztyp(100) !< Atomic numbers for the types

    integer :: lu, lp
    character(len=:), allocatable :: aux1, aatom, line
    logical :: ok

    ntyp = 0

    ! open
    lu = fopen_read(filename)

    ! read the atoms
    do while (getline_raw(lu,line))
       lp = 1
       aux1 = getword(line,lp)
       aatom = getword(line,lp)
       ntyp = ntyp + 1
       ztyp(ntyp) = aatom(1:2)
       line = ""
       do while (.not. (trim(adjustl(line)) == 'End of Dataset'))
          ok = getline_raw(lu,line,.true.)
       end do
    end do

    ! close
    call fclose(lu)

  end subroutine struct_read_potcar

  !> Read the structure from an abinit DEN file (and similar files: ELF, LDEN, etc.)
  subroutine struct_read_abinit(c,file,mol)
    use struct_basic
    use types
    use tools_math
    use tools_io
    use abinit_private
    use param

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?

    integer :: lu, fform0
    type(hdr_type) :: hdr
    integer :: i
    real*8 :: rmat(3,3), g(3,3)

    lu = fopen_read(file,"unformatted")

    ! read the header of the DEN file
    call hdr_io(fform0,hdr,1,lu)

    ! cell parameters
    rmat = hdr%rprimd(:,:)
    g = matmul(transpose(rmat),rmat)
    c%crys2car = rmat
    c%car2crys = matinv(rmat)
    do i = 1, 3
       c%aa(i) = sqrt(g(i,i))
    end do
    c%bb(1) = acos(g(2,3) / c%aa(2) / c%aa(3)) * 180d0 / pi
    c%bb(2) = acos(g(1,3) / c%aa(1) / c%aa(3)) * 180d0 / pi
    c%bb(3) = acos(g(1,2) / c%aa(1) / c%aa(2)) * 180d0 / pi

    ! abinit always uses primitive cells, even when it does not (chkprim = -1)
    c%ncv = 0
    c%lcent = 1

    ! atoms
    c%nneq = hdr%natom
    if (c%nneq > size(c%at)) call realloc(c%at,c%nneq)
    do i = 1, c%nneq
       c%at(i)%x = hdr%xred(:,i)
       c%at(i)%z = nint(hdr%znucltypat(hdr%typat(i)))
       write (c%at(i)%name,'(A2,I2.2)') nameguess(c%at(i)%z),i
    end do

    ! symmetry and transformations
    c%neqv = hdr%nsym
    do i = 1, hdr%nsym
       c%rotm(:,1:3,i) = real(hdr%symrel(:,:,i),8)
       c%rotm(:,4,i) = hdr%tnons(:,i)
    end do
    if (hdr%nsym > 0) c%havesym = 2

    ! charges and pseudopotential charges
    if (hdr%ntypat /= hdr%npsp) call ferror('struct_read_abinit','Can not handle ntypat/=npsp (?)',faterr,file)
    do i = 1, c%nneq
       ! c%at(i)%zpsp = nint(hdr%zionpsp(hdr%typat(i))) ! changed this behavior for consistency
       c%at(i)%zpsp = -1
       c%at(i)%qat = 0
    end do

    ! clean up
    call fclose(lu)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_abinit

  ! The following code has been adapted from the elk distribution, version 1.3.2
  ! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU General Public License.
  subroutine struct_read_elk(c,filename,mol)
    use struct_basic
    use tools_io
    use tools_math
    use types
    use param

    type(crystal), intent(inout) :: c !< crystal
    character*(*), intent(in) :: filename !< input filename
    logical, intent(in) :: mol !< is this a molecule?

    character(len=:), allocatable :: line, atname
    integer :: lu, i, zat, j, lp
    real*8 :: g(3,3)
    integer :: natoms
    integer :: nspecies
    logical :: ok

    lu = fopen_read(filename)

    ! ignore the 'scale' stuff
    do i = 1, 14
       read(lu,*)
    end do

    read(lu,'(3G18.10)') c%crys2car(:,1)
    read(lu,'(3G18.10)') c%crys2car(:,2)
    read(lu,'(3G18.10)') c%crys2car(:,3)
    g = matmul(transpose(c%crys2car),c%crys2car)
    c%car2crys = matinv(c%crys2car)

    do i = 1, 3
       c%aa(i) = sqrt(g(i,i))
    end do
    c%bb(1) = acos(g(2,3)/c%aa(2)/c%aa(3)) * 180 / pi
    c%bb(2) = acos(g(1,3)/c%aa(1)/c%aa(3)) * 180 / pi
    c%bb(3) = acos(g(1,2)/c%aa(1)/c%aa(2)) * 180 / pi

    ok = getline_raw(lu,line,.true.)
    ok = getline_raw(lu,line,.true.)
    if (equal(line,'molecule')) call ferror('read_elk_geometry','Isolated molecules not supported',faterr,line)

    read(lu,'(I4)') nspecies
    do i = 1, nspecies
       ok = getline_raw(lu,line,.true.)
       lp = 1
       atname = getword(line,lp)
       do j = 1, len(atname)
          if (atname(j:j) == "'") atname(j:j) = " "
          if (atname(j:j) == '"') atname(j:j) = " "
       end do
       zat = zatguess(atname)
       if (zat == -1) call ferror('read_elk_geometry','Species file name must start with an atomic symbol',faterr,filename)
       read(lu,*) natoms
       do j = 1, natoms
          c%nneq = c%nneq + 1
          if (c%nneq > size(c%at)) call realloc(c%at,2*size(c%at))
          read(lu,*) c%at(c%nneq)%x(:)
          c%at(c%nneq)%r = c%x2c(c%at(c%nneq)%x)
          c%at(c%nneq)%z = zat
          c%at(c%nneq)%zpsp = -1
          c%at(c%nneq)%qat = 0
          c%at(c%nneq)%name = trim(adjustl(atname))
       end do
    end do

    ! initialize charge and pseudopotential charge
    c%at(1:c%nneq)%zpsp = -1
    c%at(1:c%nneq)%qat = 0

    ! clean up
    call fclose(lu)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_elk

  !> Read the structure from an xyz/wfn/wfx file
  subroutine struct_read_mol(c,file,fmt,rborder,docube)
    use struct_basic
    use wfn_private
    use grid1_tools
    use types
    use tools_io
    use param

    type(crystal), intent(inout) :: c !< crystal
    character*(*), intent(in) :: file !< Input file name
    character*3, intent(in) :: fmt !< wfn/wfx/xyz
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic

    integer, allocatable :: iz(:)
    real*8, allocatable :: x(:,:)
    integer :: lu, i, j, i1, i2
    real*8 :: zreal
    character(len=:), allocatable :: line
    character*4 :: atsym
    character*4 :: orbtyp

    ! open the file
    lu = fopen_read(file)

    if (equal(fmt,'xyz')) then
       ! xyz
       read (lu,*) c%nneq
       read (lu,*)
       call realloc(c%at,c%nneq)
       do i = 1, c%nneq
          read (lu,*) atsym, c%at(i)%x
          c%at(i)%z = zatguess(atsym)
          if (c%at(i)%z == 0) &
             call ferror('struct_read_mol','could not determine atomic number',faterr)
          c%at(i)%name = trim(adjustl(atsym))
          c%at(i)%x = c%at(i)%x / bohrtoa
       end do
    elseif (equal(fmt,'wfn')) then
       ! wfn
       read (lu,*)
       read (lu,101) orbtyp, i1, i2, c%nneq
       if (c%nneq <= 0) &
          call ferror('struct_read_mol','wrong number of atoms',faterr)
       call realloc(c%at,c%nneq)
       do i = 1, c%nneq
          read(lu,106) atsym, c%at(i)%x, zreal
          c%at(i)%z = zatguess(atsym)
          c%at(i)%name = trim(adjustl(atsym))
       end do
101    format (4X,A4,10X,3(I5,15X))
106    format(2X,A2,20X,3F12.8,10X,F5.1)
    elseif (equal(fmt,'wfx')) then
       ! wfx
       do while (getline_raw(lu,line))
          if (line(1:1) == "<" .and. line(2:2) /= "/") then
             if (trim(line) == "<Number of Nuclei>") then
                read (lu,*) c%nneq
             endif
          endif
       enddo

       rewind(lu)
       allocate(x(3,c%nneq),iz(c%nneq))
       if (c%nneq == 0) &
          call ferror("struct_read_mol","Number of Nuclei tag not found",faterr)
       call realloc(c%at,c%nneq)
       do while (.true.)
          read(lu,'(A)',end=20) line
          line = adjustl(line)
          if (line(1:1) == "<" .and. line(2:2) /= "/") then
             if (trim(line) == "<Atomic Numbers>") then
                iz = wfx_read_integers(lu,c%nneq)
                do i = 1, c%nneq
                   c%at(i)%z = iz(i)
                end do
             elseif (trim(line) == "<Nuclear Cartesian Coordinates>") then
                x = reshape(wfx_read_reals1(lu,3*c%nneq),shape(x))
                do i = 1, c%nneq
                   c%at(i)%x = x(:,i)
                end do
             endif
          endif
       enddo
20     continue
       deallocate(x,iz)
    end if

    ! close the file
    call fclose(lu)

    ! fill the rest of the info
    call fill_molecule(c,rborder,docube)

    ! Do note use atomic charges or pseudos
    do i = 1, c%nneq
       c%at(i)%zpsp = -1
       c%at(i)%qat = 0
    end do

  end subroutine struct_read_mol

  !> Read the structure from a quantum espresso output
  subroutine struct_read_qeout(c,file,mol)
    use struct_basic
    use tools_io
    use tools_math
    use param
    use types

    type(crystal), intent(inout) :: c !< crystal
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?

    integer :: lu, nstrucs, is0, ideq, i, k
    character(len=:), allocatable :: line
    integer :: ibrav, nat, ntyp, id, idum
    real*8 :: alat, r(3,3), qaux, g(3,3)
    logical :: ok
    character*(10), allocatable :: attyp(:), atn(:)
    integer, allocatable :: zpsptyp(:)
    real*8, allocatable :: x(:,:)

    lu = fopen_read(file)

    ! first pass: read the number of structures
    nstrucs = 0
    do while (getline_raw(lu,line))
       if (index(line,"Title:") > 0) then
          nstrucs = nstrucs + 1
       end if
    end do

    ! which of them? -> last one
    is0 = nstrucs

    ! rewind and read the correct structure
    rewind(lu)
    nstrucs = 0
    a: do while (getline_raw(lu,line))
       if (index(line,"Title:") > 0) then
          nstrucs = nstrucs + 1
       end if
       if (nstrucs == is0) then
          ! read header
          do while (getline_raw(lu,line))
             ideq = index(line,"=") + 1
             if (index(line,"bravais-lattice index") > 0) then
                ok = isinteger(ibrav,line,ideq)
             elseif (index(line,"lattice parameter (alat)") > 0) then
                ok = isreal(alat,line,ideq)
             elseif (index(line,"number of atoms/cell") > 0) then
                ok = isinteger(nat,line,ideq)
             elseif (index(line,"number of atomic types") > 0) then
                ok = isinteger(ntyp,line,ideq)
             elseif (index(line,"crystal axes:") > 0) then
                do i = 1, 3
                   ok = getline_raw(lu,line,.true.)
                   ideq = index(line,"(",.true.) + 1
                   ok = isreal(r(i,1),line,ideq)
                   ok = ok.and.isreal(r(i,2),line,ideq)
                   ok = ok.and.isreal(r(i,3),line,ideq)
                end do
                exit
             end if
          end do

          ! allocate types
          allocate(attyp(ntyp),zpsptyp(ntyp))
          allocate(atn(nat), x(3,nat))
          zpsptyp = 0

          ! atomic species
          do while (getline_raw(lu,line))
             if (index(line,"atomic species   valence    mass     pseudopotential")>0) then
                do i = 1, ntyp
                   ok = getline_raw(lu,line,.true.)
                   read (line,*) attyp(i), qaux
                   zpsptyp(i) = nint(qaux)
                end do
                exit
             end if
          end do

          ! coordinates
          do while (getline_raw(lu,line))
             if (index(line,"Cartesian axes")>0) then
                ok = getline_raw(lu,line,.true.)
                ok = getline_raw(lu,line,.true.)
                do i = 1, nat
                   ok = getline_raw(lu,line,.true.)
                   read(line,*) idum, atn(i)
                   line = line(index(line,"(",.true.)+1:)
                   read(line,*) x(:,i)
                end do
             end if
          end do

          ! last structure, done
          exit a
       end if
    end do a

    ! fill struct quantities
    ! cell
    r = r * alat
    g = matmul(r,transpose(r))
    do i = 1, 3
       c%aa(i) = sqrt(g(i,i))
    end do
    c%bb(1) = acos(g(2,3)/c%aa(2)/c%aa(3)) * 180d0 / pi
    c%bb(2) = acos(g(1,3)/c%aa(1)/c%aa(3)) * 180d0 / pi
    c%bb(3) = acos(g(1,2)/c%aa(1)/c%aa(2)) * 180d0 / pi
    c%crys2car = transpose(r)
    c%car2crys = matinv(c%crys2car)
    ! atoms
    c%nneq = nat
    call realloc(c%at,c%nneq)
    do i = 1, nat
       c%at(i)%x = c%c2x(x(:,i) * alat)
       c%at(i)%x = c%at(i)%x - floor(c%at(i)%x)

       ! identify type
       id = 0
       do k = 1, ntyp
          if (trim(atn(i)) == trim(attyp(k))) then
             id = k
             exit
          end if
       end do
       if (id == 0) call ferror('struct_read_qeout','atom type not found',faterr)
       c%at(i)%zpsp = zpsptyp(id)
       c%at(i)%qat = 0
       c%at(i)%name = attyp(id)
       c%at(i)%z = zatguess(c%at(i)%name)
       if (c%at(i)%z == 0) &
          call ferror('struct_read_qeout','could not determine atomic number',faterr)
    end do

    ! deallocate
    deallocate(attyp,zpsptyp,atn,x)

    ! no symmetry
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,1:3,1) = eye
    c%ncv = 0
    c%lcent = 1

    ! close the shop
    call fclose(lu)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_qeout

  !> Read the structure from a quantum espresso input
  subroutine struct_read_qein(c0,file,mol)
    ! This subroutine has been adapted from parts of the Quantum
    ! ESPRESSO code, version 4.3.2.  
    ! Copyright (C) 2002-2009 Quantum ESPRESSO group
    ! This file is distributed under the terms of the
    ! GNU General Public License. See the file `License'
    ! in the root directory of the present distribution,
    ! or http://www.gnu.org/copyleft/gpl.txt .
    use struct_basic
    use tools_io
    use tools_math
    use param
    use types

    type(crystal), intent(inout) :: c0 !< crystal
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?

    integer, parameter :: dp = selected_real_kind(14,200)
    integer, parameter :: ntypx = 10
    integer, parameter :: nsx = ntypx
    integer, parameter :: nspinx = 2
    integer, parameter :: lmaxx = 3
    integer, parameter :: lqmax = 2*lmaxx+1

    ! from QE
    ! namelist control
    character(len=80) :: title, calculation, verbosity, restart_mode,&
       disk_io
    integer :: nstep, iprint, isave, ndr, ndw, gdir, nppstr, nberrycyc, &
       printwfc
    logical :: tstress, tprnfor, tefield, tefield2, lelfield, dipfield, &
       lberry, wf_collect, saverho, tabps, lkpoint_dir, use_wannier, &
       lecrpa
    real*8 :: dt, refg, max_seconds, ekin_conv_thr, etot_conv_thr, &
       forc_conv_thr
    character(len=256) :: outdir, prefix, pseudo_dir, wfcdir, vdw_table_name
    namelist /control/ title, calculation, verbosity, restart_mode,  &
       nstep, iprint, isave, tstress, tprnfor, dt, ndr, ndw, outdir,   &
       prefix, wfcdir, max_seconds, ekin_conv_thr, etot_conv_thr,      &
       forc_conv_thr, pseudo_dir, disk_io, tefield, dipfield, lberry,  &
       gdir, nppstr, wf_collect, printwfc, lelfield, nberrycyc, refg,  &
       tefield2, saverho, tabps, lkpoint_dir, use_wannier, lecrpa,     &
       vdw_table_name
    
    ! namelist system
    integer :: ibrav = 14
    real*8 :: celldm(6) = 0.d0
    real*8 :: a, b, c, cosab, cosac, cosbc
    integer :: nat = 0
    integer :: ntyp = 0
    real*8 :: tot_charge, tot_magnetization, ecutwfc, ecutrho, degauss, &
       ecfixed, qcutz, q2sigma, starting_magnetization(nsx), &
       starting_ns_eigenvalue(lqmax,nspinx,nsx), hubbard_u(nsx), &
       hubbard_alpha(nsx), a_pen(10,nspinx), sigma_pen(10), alpha_pen(10), &
       emaxpos, eopreg, eamp, lambda, fixed_magnetization(3), angle1(nsx), &
       angle2(nsx), b_field(3), sic_epsilon, sic_alpha, london_s6, london_rcut, &
       xdm_a1, xdm_a2, ts_sr, esm_efield, esm_w
    integer :: nbnd, nr1, nr2, nr3, nr1s, nr2s, nr3s, nr1b, nr2b, nr3b, &
       nspin, edir, report, xdm_usehigh, esm_nfit, esm_debug_gpmax
    character(len=80) :: occupations, smearing, input_dft, u_projection_type, &
       constrained_magnetization, sic, assume_isolated
    logical :: nosym, noinv, nosym_evc, force_symmorphic, lda_plus_u, la2f, &
       step_pen, noncolin, lspinorb, starting_spin_angle, no_t_rev, force_pairing, &
       spline_ps, one_atom_occupations, london, xdm, xdm_onlyc, xdm_fixc6, &
       xdm_usec9, ts, ts_onlyc, esm_debug
    character(len=3) :: esm_bc
    namelist /system/ ibrav, celldm, a, b, c, cosab, cosac, cosbc, nat, &
       ntyp, nbnd, ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s,  &
       nr3s, nr1b, nr2b, nr3b, nosym, nosym_evc, noinv,                 &
       force_symmorphic, starting_magnetization,                        &
       occupations, degauss, nspin, ecfixed,              &
       qcutz, q2sigma, lda_plus_u, hubbard_u, hubbard_alpha,            &
       edir, emaxpos, eopreg, eamp, smearing, starting_ns_eigenvalue,   &
       u_projection_type, input_dft, la2f, assume_isolated,             &
       noncolin, lspinorb, starting_spin_angle, lambda, angle1, angle2, &
       report,              &
       constrained_magnetization, b_field, fixed_magnetization,         &
       sic, sic_epsilon, force_pairing, sic_alpha,                      &
       tot_charge, tot_magnetization,                                   &
       spline_ps, one_atom_occupations, london, london_s6, london_rcut, &
       xdm, xdm_onlyc, xdm_fixc6, xdm_usec9, xdm_usehigh, xdm_a1,       &
       xdm_a2, ts, ts_onlyc, ts_sr,                                     &
       step_pen, a_pen, sigma_pen, alpha_pen, no_t_rev,                 &
       esm_bc, esm_efield, esm_w, esm_nfit, esm_debug, esm_debug_gpmax

    ! namelist electrons
    real*8 :: emass, emass_cutoff, ortho_eps, electron_damping, ekincw, fnosee, &
       ampre, grease, diis_hcut, diis_wthr, diis_delt, diis_fthr, diis_temp, &
       diis_achmix, diis_g0chmix, diis_g1chmix, diis_rothr, diis_ethr, mixing_beta,&
       diago_thr_init, conv_thr, lambda_cold, fermi_energy, rotmass, occmass,&
       occupation_damping, rotation_damping, etresh, passop, efield, efield_cart(3),&
       efield2
    character(len=80) :: orthogonalization, electron_dynamics, electron_velocities,&
       electron_temperature, startingwfc, mixing_mode, diagonalization, startingpot,&
       rotation_dynamics, occupation_dynamics
    integer :: ortho_max, electron_maxstep, diis_size, diis_nreset, diis_maxstep, &
       diis_nchmix, diis_nrot(3), mixing_ndim, diago_cg_maxiter, diago_david_ndim, &
       mixing_fixed_ns, n_inner, niter_cold_restart, maxiter, niter_cg_restart, &
       epol, epol2
    logical :: diis_rot, diis_chguess, diago_full_acc, tcg, real_space, tqr,&
       occupation_constraints
    namelist /electrons/ emass, emass_cutoff, orthogonalization, &
       electron_maxstep, ortho_eps, ortho_max, electron_dynamics,   &
       electron_damping, electron_velocities, electron_temperature, &
       ekincw, fnosee, ampre, grease,                               &
       diis_size, diis_nreset, diis_hcut,                           &
       diis_wthr, diis_delt, diis_maxstep, diis_rot, diis_fthr,     &
       diis_temp, diis_achmix, diis_g0chmix, diis_g1chmix,          &
       diis_nchmix, diis_nrot, diis_rothr, diis_ethr, diis_chguess, &
       mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns,      &
       tqr, diago_cg_maxiter, diago_david_ndim, diagonalization ,   &
       startingpot, startingwfc , conv_thr,                         &
       diago_thr_init, n_inner, fermi_energy, rotmass, occmass,     &
       rotation_damping, occupation_damping, rotation_dynamics,     &
       occupation_dynamics, tcg, maxiter, etresh, passop, epol,     &
       efield, epol2, efield2, diago_full_acc,                      &
       occupation_constraints, niter_cg_restart,                    &
       niter_cold_restart, lambda_cold, efield_cart, real_space
    
    ! namelist ions
    character(len=80) :: phase_space, ion_dynamics, ion_positions, ion_velocities,&
       ion_temperature, pot_extrapolation, wfc_extrapolation
    integer, parameter :: nhclm   = 4
    integer, parameter :: max_nconstr = 100
    integer :: nhpcl, nhptyp, nhgrp(nsx), ndega, ion_nstepe, ion_maxstep, nraise,&
       bfgs_ndim, fe_nstep, sw_nstep, eq_nstep
    real*8 :: ion_radius(nsx), ion_damping, tempw, fnosep(nhclm), tolp, fnhscl(nsx),&
       amprp(nsx), greasp, upscale, delta_t, trust_radius_max, trust_radius_min,&
       trust_radius_ini, w_1, w_2, sic_rloc, g_amplitude, fe_step(max_nconstr)
    logical :: tranp(nsx), refold_pos, remove_rigid_rot
    namelist /ions/ phase_space, ion_dynamics, ion_radius, ion_damping,  &
       ion_positions, ion_velocities, ion_temperature,      &
       tempw, fnosep, nhgrp, fnhscl, nhpcl, nhptyp, ndega, tranp,   &
       amprp, greasp, tolp, ion_nstepe, ion_maxstep,        &
       refold_pos, upscale, delta_t, pot_extrapolation,     &
       wfc_extrapolation, nraise, remove_rigid_rot,         &
       trust_radius_max, trust_radius_min,                  &
       trust_radius_ini, w_1, w_2, bfgs_ndim, sic_rloc,     &
       fe_step, fe_nstep, sw_nstep, eq_nstep, g_amplitude

    ! namelist cell
    character(len=80) :: cell_parameters, cell_dynamics, cell_velocities, &
       cell_temperature, cell_dofree
    real(dp) :: press, wmass, temph, fnoseh, greash, cell_factor, cell_damping,&
       press_conv_thr
    integer :: cell_nstepe = 1
    namelist /cell/ cell_parameters, cell_dynamics, cell_velocities, &
       press, wmass, cell_temperature, temph, fnoseh,   &
       cell_dofree, greash, cell_factor, cell_nstepe,   &
       cell_damping, press_conv_thr

    ! local to this routine
    integer :: lu, ios, lp, i, j
    character(len=:), allocatable :: line, word
    logical :: havecell
    real*8 :: r(3,3), g(3,3)
    integer :: iunit
    integer, parameter :: icrystal = 1
    integer, parameter :: ibohr = 2
    integer, parameter :: iang = 3
    integer, parameter :: ialat = 4

    ! open
    lu = fopen_read(file)

    ! read the namelists
    read(lu,control,iostat=ios)
    if (ios/=0) call ferror("struct_read_qein","wrong namelist control",faterr)
    read(lu,system,iostat=ios)
    if (ios/=0) call ferror("struct_read_qein","wrong namelist system",faterr)
    read(lu,electrons,iostat=ios)
    if (ios/=0) call ferror("struct_read_qein","wrong namelist electrons",faterr)
    if (trim(calculation)=='relax'.or.trim(calculation)=='md'.or.&
        trim(calculation)=='vc-relax'.or.trim(calculation)=='vc-md'.or.&
        trim(calculation)=='cp'.or.trim(calculation)=='vc-cp'.or.&
        trim(calculation)=='smd'.or.trim(calculation)=='cp-wf') then
       read(lu,ions,iostat=ios)
       if (ios/=0) call ferror("struct_read_qein","wrong namelist ions",faterr)
    endif
    if (trim(calculation)=='vc-relax'.or.trim(calculation)=='vc-md'.or.&
        trim(calculation)=='vc-cp') then
       read(lu,cell,iostat=ios)
       if (ios/=0) call ferror("struct_read_qein","wrong namelist ions",faterr)
    end if
    
    ! allocate space for atoms
    c0%nneq = nat
    call realloc(c0%at,c0%nneq)

    ! read the cards
    havecell = .false.
    do while (getline_raw(lu,line))
       line = upper(line)
       lp = 1
       word = getword(line,lp)
       if (equal(word,'ATOMIC_POSITIONS')) then
          word = getword(line,lp)
          if (equal(word,"CRYSTAL")) then
             iunit = icrystal
          elseif (equal(word,"BOHR")) then
             iunit = ibohr
          elseif (equal(word,"ANGSTROM")) then
             iunit = iang
          elseif (equal(word,"ALAT")) then
             iunit = ialat
          else
             iunit = ialat
          end if
          do i = 1, nat
             read (lu,*) c0%at(i)%name, c0%at(i)%x
          end do
       elseif (equal(word,'CELL_PARAMETERS')) then
          havecell = .true.
          do i = 1, 3
             read (lu,*) (r(i,j),j=1,3)
          end do
       endif
    end do

    ! figure it out
    if (ibrav == 0) then
       if (celldm(1) /= 0.D0) r = r * celldm(1)
    else
       call qe_latgen(ibrav,celldm,r(1,:),r(2,:),r(3,:))
    endif

    ! fill the cell metrics
    g = matmul(r,transpose(r))
    do i = 1, 3
       c0%aa(i) = sqrt(g(i,i))
    end do
    c0%bb(1) = acos(g(2,3)/c0%aa(2)/c0%aa(3)) * 180d0 / pi
    c0%bb(2) = acos(g(1,3)/c0%aa(1)/c0%aa(3)) * 180d0 / pi
    c0%bb(3) = acos(g(1,2)/c0%aa(1)/c0%aa(2)) * 180d0 / pi
    c0%crys2car = transpose(r)
    c0%car2crys = matinv(c0%crys2car)

    ! do the atom stuff
    do i = 1, nat
       if (iunit == ialat) then
          c0%at(i)%x = c0%c2x(c0%at(i)%x * celldm(1))
       elseif (iunit == ibohr) then
          c0%at(i)%x = c0%c2x(c0%at(i)%x)
       elseif (iunit == iang) then
          c0%at(i)%x = c0%c2x(c0%at(i)%x / bohrtoa)
       endif
       c0%at(i)%x = c0%at(i)%x - floor(c0%at(i)%x)
       c0%at(i)%zpsp = -1
       c0%at(i)%qat = 0
       c0%at(i)%z = zatguess(c0%at(i)%name)
       if (c0%at(i)%z == 0) &
          call ferror('struct_read_qeout','could not determine atomic number',faterr)
    end do

    ! no symmetry
    c0%havesym = 0
    c0%neqv = 1
    c0%rotm = 0d0
    c0%rotm(:,1:3,1) = eye
    c0%ncv = 0
    c0%lcent = 1

    ! close
    call fclose(lu)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c0)

  end subroutine struct_read_qein

  !> Read the structure from a siesta STRUCT_OUT input
  subroutine struct_read_siesta(c,file,mol)
    use struct_basic
    use types
    use tools_io
    use tools_math
    use param

    type(crystal), intent(inout) :: c !< Crystal
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?

    integer :: lu
    real*8 :: r(3,3), g(3,3)
    integer :: i, iz, idum

    ! open
    lu = fopen_read(file)

    ! the lattice vectors
    do i = 1, 3
       read (lu,*) r(i,:)
    end do
    r = r / bohrtoa

    ! the atoms
    read (lu,*) c%nneq
    call realloc(c%at,c%nneq)
    do i = 1, c%nneq
       read (lu,*) idum, iz, c%at(i)%x
       c%at(i)%z = iz
       c%at(i)%name = nameguess(iz)
       c%at(i)%zpsp = -1
       c%at(i)%qat = 0
    end do

    ! fill the cell metrics
    g = matmul(r,transpose(r))
    do i = 1, 3
       c%aa(i) = sqrt(g(i,i))
    end do
    c%bb(1) = acos(g(2,3)/c%aa(2)/c%aa(3)) * 180d0 / pi
    c%bb(2) = acos(g(1,3)/c%aa(1)/c%aa(3)) * 180d0 / pi
    c%bb(3) = acos(g(1,2)/c%aa(1)/c%aa(2)) * 180d0 / pi
    c%crys2car = transpose(r)
    c%car2crys = matinv(c%crys2car)

    ! no symmetry
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,1:3,1) = eye
    c%ncv = 0
    c%lcent = 1

    ! close
    call fclose(lu)

    ! if this is a molecule, set up the origin and the molecular cell
    if (mol) call fill_molecule_given_cell(c)

  end subroutine struct_read_siesta

  !> From QE, generate the lattice from the ibrav
  subroutine qe_latgen(ibrav,celldm,a1,a2,a3)
    ! This subroutine has been adapted from parts of the Quantum
    ! ESPRESSO code, version 4.3.2.  
    ! Copyright (C) 2002-2009 Quantum ESPRESSO group
    ! This file is distributed under the terms of the
    ! GNU General Public License. See the file `License'
    ! in the root directory of the present distribution,
    ! or http://www.gnu.org/copyleft/gpl.txt .
    use tools_io
    !-----------------------------------------------------------------------
    !     sets up the crystallographic vectors a1, a2, and a3.
    !
    !     ibrav is the structure index:
    !       1  cubic P (sc)                8  orthorhombic P
    !       2  cubic F (fcc)               9  one face centered orthorhombic
    !       3  cubic I (bcc)              10  all face centered orthorhombic
    !       4  hexagonal and trigonal P   11  body centered orthorhombic
    !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
    !       6  tetragonal P (st)          13  one face centered monoclinic
    !       7  tetragonal I (bct)         14  triclinic P
    !     Also accepted:
    !       0  "free" structure          -12  monoclinic P (unique axis: b)
    !      -5  trigonal R, threefold axis along (111) 
    !
    !     NOTA BENE: all axis sets are right-handed
    !     Boxes for US PPs do not work properly with left-handed axis
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    integer, intent(in) :: ibrav
    real(DP), intent(inout) :: celldm(6)
    real(DP), intent(out) :: a1(3), a2(3), a3(3)

    real(DP), parameter:: sr2 = 1.414213562373d0, sr3 = 1.732050807569d0
    integer :: ir
    real(DP) :: term, cbya, term1, term2, singam, sen

    if (celldm(1) <= 0.d0) call ferror('qe_latgen','wrong celldm(1)',faterr)

    ! index of bravais lattice supplied
    if (ibrav == 1) then
       ! simple cubic lattice
       a1(1)=celldm(1)
       a2(2)=celldm(1)
       a3(3)=celldm(1)
       !
    else if (ibrav == 2) then
       ! fcc lattice
       term=celldm(1)/2.d0
       a1(1)=-term
       a1(3)=term
       a2(2)=term
       a2(3)=term
       a3(1)=-term
       a3(2)=term
       !
    else if (ibrav == 3) then
       ! bcc lattice
       term=celldm(1)/2.d0
       do ir=1,3
          a1(ir)=term
          a2(ir)=term
          a3(ir)=term
       end do
       a2(1)=-term
       a3(1)=-term
       a3(2)=-term
       !
    else if (ibrav == 4) then
       ! hexagonal lattice
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       cbya=celldm(3)
       a1(1)=celldm(1)
       a2(1)=-celldm(1)/2.d0
       a2(2)=celldm(1)*sr3/2.d0
       a3(3)=celldm(1)*cbya
       !
    else if (ibrav == 5) then
       ! trigonal lattice, threefold axis along c (001)
       if (celldm (4) <= -0.5d0 .or. celldm (4) >= 1) &
          call ferror('qe_latgen','wrong celldm(4)',faterr)
       !
       term1=sqrt(1.d0+2.d0*celldm(4))
       term2=sqrt(1.d0-celldm(4))
       a2(2)=sr2*celldm(1)*term2/sr3
       a2(3)=celldm(1)*term1/sr3
       a1(1)=celldm(1)*term2/sr2
       a1(2)=-a1(1)/sr3
       a1(3)= a2(3)
       a3(1)=-a1(1)
       a3(2)= a1(2)
       a3(3)= a2(3)
       !
    else if (ibrav ==-5) then
       ! trigonal lattice, threefold axis along (111)
       if (celldm (4) <= -0.5d0 .or. celldm (4) >= 1) &
          call ferror('qe_latgen','wrong celldm(4)',faterr)
       !
       term1 = sqrt(1.0_dp + 2.0_dp*celldm(4))
       term2 = sqrt(1.0_dp - celldm(4))
       a1(1) = celldm(1)*(term1-2.0_dp*term2)/3.0_dp
       a1(2) = celldm(1)*(term1+term2)/3.0_dp
       a1(3) = a1(2)
       a2(1) = a1(3)
       a2(2) = a1(1)
       a2(3) = a1(2)
       a3(1) = a1(2)
       a3(2) = a1(3)
       a3(3) = a1(1)
    else if (ibrav == 6) then
       ! tetragonal lattice
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       cbya=celldm(3)
       a1(1)=celldm(1)
       a2(2)=celldm(1)
       a3(3)=celldm(1)*cbya
       !
    else if (ibrav == 7) then
       ! body centered tetragonal lattice
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       cbya=celldm(3)
       a2(1)=celldm(1)/2.d0
       a2(2)=a2(1)
       a2(3)=cbya*celldm(1)/2.d0
       a1(1)= a2(1)
       a1(2)=-a2(1)
       a1(3)= a2(3)
       a3(1)=-a2(1)
       a3(2)=-a2(1)
       a3(3)= a2(3)
       !
    else if (ibrav == 8) then
       ! Simple orthorhombic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       a1(1)=celldm(1)
       a2(2)=celldm(1)*celldm(2)
       a3(3)=celldm(1)*celldm(3)
       !
    else if (ibrav == 9) then
       ! One face centered orthorhombic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       a1(1) = 0.5d0 * celldm(1)
       a1(2) = a1(1) * celldm(2)
       a2(1) = - a1(1)
       a2(2) = a1(2)
       a3(3) = celldm(1) * celldm(3)
       !
    else if (ibrav == 10) then
       ! All face centered orthorhombic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       a2(1) = 0.5d0 * celldm(1)
       a2(2) = a2(1) * celldm(2)
       a1(1) = a2(1)
       a1(3) = a2(1) * celldm(3)
       a3(2) = a2(1) * celldm(2)
       a3(3) = a1(3)
       !
    else if (ibrav == 11) then
       ! Body centered orthorhombic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       !
       a1(1) = 0.5d0 * celldm(1)
       a1(2) = a1(1) * celldm(2)
       a1(3) = a1(1) * celldm(3)
       a2(1) = - a1(1)
       a2(2) = a1(2)
       a2(3) = a1(3)
       a3(1) = - a1(1)
       a3(2) = - a1(2)
       a3(3) = a1(3)
       !
    else if (ibrav == 12) then
       ! Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       if (abs(celldm(4))>=1.d0) call ferror('qe_latgen','wrong celldm(4)',faterr)
       !
       sen=sqrt(1.d0-celldm(4)**2)
       a1(1)=celldm(1)
       a2(1)=celldm(1)*celldm(2)*celldm(4)
       a2(2)=celldm(1)*celldm(2)*sen
       a3(3)=celldm(1)*celldm(3)
       !
    else if (ibrav ==-12) then
       ! Simple monoclinic lattice, unique axis: b (more common)
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       if (abs(celldm(5))>=1.d0) call ferror('qe_latgen','wrong celldm(5)',faterr)
       !
       sen=sqrt(1.d0-celldm(5)**2)
       a1(1)=celldm(1)
       a2(2)=celldm(1)*celldm(2)
       a3(1)=celldm(1)*celldm(3)*celldm(5)
       a3(3)=celldm(1)*celldm(3)*sen
       !
    else if (ibrav == 13) then
       ! One face centered monoclinic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       if (abs(celldm(4))>=1.d0) call ferror('qe_latgen','wrong celldm(4)',faterr)
       !
       sen = sqrt( 1.d0 - celldm(4) ** 2 )
       a1(1) = 0.5d0 * celldm(1) 
       a1(3) =-a1(1) * celldm(3)
       a2(1) = celldm(1) * celldm(2) * celldm(4)
       a2(2) = celldm(1) * celldm(2) * sen
       a3(1) = a1(1)
       a3(3) =-a1(3)
       !
    else if (ibrav == 14) then
       ! Triclinic lattice
       if (celldm (2) <= 0.d0) call ferror('qe_latgen','wrong celldm(2)',faterr)
       if (celldm (3) <= 0.d0) call ferror('qe_latgen','wrong celldm(3)',faterr)
       if (abs(celldm(4))>=1.d0) call ferror('qe_latgen','wrong celldm(4)',faterr)
       if (abs(celldm(5))>=1.d0) call ferror('qe_latgen','wrong celldm(5)',faterr)
       if (abs(celldm(6))>=1.d0) call ferror('qe_latgen','wrong celldm(6)',faterr)
       !
       singam=sqrt(1.d0-celldm(6)**2)
       term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)-celldm(4)**2-celldm(5)**2-celldm(6)**2)
       if (term < 0.d0) call ferror('qe_latgen','celldm do not make sense, check your data',faterr)
       term= sqrt(term/(1.d0-celldm(6)**2))
       a1(1)=celldm(1)
       a2(1)=celldm(1)*celldm(2)*celldm(6)
       a2(2)=celldm(1)*celldm(2)*singam
       a3(1)=celldm(1)*celldm(3)*celldm(5)
       a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
       a3(3)=celldm(1)*celldm(3)*term
       !
    else
       call ferror('qe_latgen','nonexistent bravais lattice',faterr)
    end if

  end subroutine qe_latgen

  !> Wrapper to the spgs module
  subroutine spgs_wrap(c,spg,usespgr)
    use struct_basic
    use spgs
    use global
    use types
    use tools_io
    use param

    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: spg
    logical, intent(in) :: usespgr

    call spgs_driver(spg,usespgr)
    c%lcent = spgs_lcent
    c%ncv = spgs_ncv
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,c%ncv))
    c%cen(:,1:c%ncv) = real(spgs_cen(:,1:c%ncv),8) / 12d0
    c%neqv = spgs_n
    c%rotm = real(spgs_m,8)
    c%rotm(:,4,:) = c%rotm(:,4,:) / 12d0
    c%havesym = 2

  end subroutine spgs_wrap

  !> Using the atomic position (bohr) in the field cr%x directly read
  !> from a molecular description (a xyz file, for example), center
  !> the molecule, set c%aa and c%bb by applying a border to it, set
  !> the c%ismolecule, the c%molx0, and the c%molborder, and
  !> deactivate the symmetry. Used in two cases: when the molecule is
  !> read from an external file or when given in the critic2 input.
  !> If docube is true, make a cubic cell instead of a parallelepiped.
  subroutine fill_molecule(c,rborder,docube)
    use struct_basic
    use types
    use param

    type(crystal), intent(inout) :: c
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    
    real*8 :: xmin(3), xmax(3), xcm(3)
    integer :: i, j

    ! minimum, maximum, and center-of-mass
    xmax = -1d40
    xmin =  1d40
    xcm = 0d0
    do i = 1, c%nneq
       do j = 1, 3
          xmax(j) = max(c%at(i)%x(j)+rborder,xmax(j))
          xmin(j) = min(c%at(i)%x(j)-rborder,xmin(j))
       end do
       xcm = xcm + c%at(i)%x
    end do
    xcm = xcm / c%nneq
    if (docube) then
       xmin = minval(xmin)
       xmax = maxval(xmax)
    end if

    ! write the positions and the rest of the information
    c%aa = xmax - xmin
    c%bb = 90d0
    do i = 1, c%nneq
       c%at(i)%x = (c%at(i)%x-xcm+0.5d0*c%aa) / c%aa
    end do

    ! Keep the (1/2,1/2,1/2) translation we applied; this is a moleucle
    c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa + xcm 

    ! Set up the molecular cell. c%molborder is in fractional coordinates
    ! and gives the position of the molecular cell in each axis. By default,
    ! choose the molecular cell as the minimal encompassing cell for the molecule
    ! plus 80% of the border or 2 bohr, whichever is larger. The molecular cell 
    ! can not exceed the actual unit cell
    c%molborder = max(rborder - max(2d0,0.8d0 * rborder),0d0) / (xmax - xmin)

    ! no symmetry for now
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,1:3,1) = eye
    c%ncv = 0
    c%lcent = 1
    
  end subroutine fill_molecule

  !> Fills the remaining information for a molecule loaded from a crystal
  !> cell description (e.g. a cube file). Sets the origin to the given x0
  !> or (if not present) to the center of the cell. Then, sets up the
  !> molecular cell. Do not allow molecules in non-orthogonal cells.
  subroutine fill_molecule_given_cell(c,x0)
    use struct_basic
    use tools_io
    type(crystal), intent(inout) :: c
    real*8, intent(in), optional :: x0(3)

    integer :: i, j
    real*8 :: xmin(3)

    if (any(abs(c%bb - 90d0) > 1d-5)) &
       call ferror('fill_molecule_given_cell','Can not use MOLECULE with a non-orthogonal cell',faterr)

    ! save the origin
    if (present(x0)) then
       c%molx0 = x0
    else
       c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa 
    endif

    ! Set up the molecular cell. c%molborder is in fractional coordinates
    ! and gives the position of the molecular cell in each axis. By default,
    ! choose the molecular cell as the minimal encompassing cell for the molecule
    ! plus 80% of the border or 2 bohr, whichever is larger. The molecular cell 
    ! can not exceed the actual unit cell
    xmin =  1d40
    do i = 1, c%nneq
       do j = 1, 3
          xmin(j) = min(c%at(i)%x(j),xmin(j))
          xmin(j) = min(1d0-max(c%at(i)%x(j),1d0-xmin(j)),xmin(j))
       end do
    end do

    c%molborder = max(xmin - max(0.8d0 * xmin,2d0/c%aa),0d0)

  end subroutine fill_molecule_given_cell

end module struct_readers
