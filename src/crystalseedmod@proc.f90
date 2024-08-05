! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Crystal seed class, contains the structure file readers.
submodule (crystalseedmod) proc
  implicit none

  !xx! private subroutines
  ! subroutine read_all_cif(file,mol,errmsg,nseed,mseed,seed0,dblock,ti)
  ! subroutine read_all_mol2(file,errmsg,nseed,mseed,seed0,name,ti)
  ! subroutine read_all_qeout(nseed,seed,file,mol,istruct,errmsg,ti)
  ! subroutine read_all_xyz(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_log(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_aimsout(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_castep_geom(nseed,seed,file,errmsg,ti)
  ! subroutine read_pdb_geometry(file,nat,x,z,name,errmsg,ti)
  ! subroutine read_zmat_geometry(file,nat,x,z,name,errmsg,ti)
  ! function which_out_format(file,ti)
  ! subroutine which_in_format(file,isformat,ti)
  ! function string_to_symop(str)

contains

  !> Deallocate arrays in the seed
  module subroutine seed_end(seed)
    use param, only: isformat_unknown
    class(crystalseed), intent(inout) :: seed !< Crystal seed output

    seed%isused = .false.
    seed%file = ""
    seed%name = ""
    seed%isformat = isformat_unknown
    seed%nat = 0
    if (allocated(seed%x)) deallocate(seed%x)
    if (allocated(seed%is)) deallocate(seed%is)
    seed%nspc = 0
    if (allocated(seed%spc)) deallocate(seed%spc)
    seed%useabr = 0
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.
    seed%neqv = 0
    seed%ncv = 0
    if (allocated(seed%cen)) deallocate(seed%cen)
    if (allocated(seed%rotm)) deallocate(seed%rotm)
    seed%ismolecule = .false.
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%energy = huge(1d0)
    seed%pressure = huge(1d0)

  end subroutine seed_end

  !> Parse a crystal environment
  module subroutine parse_crystal_env(seed,lu,oksyn)
    use spglib, only: spg_get_hall_number_from_symbol, spg_get_symmetry_from_database
    use global, only: eval_next, dunit0, iunit, iunit_isdef, iunit_bohr
    use tools_math, only: matinv
    use tools_io, only: uin, getline, ucopy, lgetword, equal, ferror, faterr,&
       getword, lower, isinteger, string, nameguess, zatguess, equali
    use param, only: bohrtoa, isformat_from_input
    use types, only: realloc

    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    integer, intent(in) :: lu !< Logical unit for input
    logical, intent(out) :: oksyn !< Was there a syntax error?

    character(len=:), allocatable :: word, aux, line, name, errmsg
    character*255, allocatable :: sline(:)
    integer :: i, j, k, lp, nsline, luout, iat, lp2, iunit0, it
    integer :: hnum, ier
    real*8 :: rmat(3,3), scal, ascal, x(3), xdif(3)
    logical :: ok, goodspg
    logical :: isset
    real*8 :: rot(3,4)

    call seed%end()
    if (iunit_isdef) then
       iunit0 = iunit_bohr
    else
       iunit0 = iunit
    end if

    oksyn = .false.
    goodspg = .false.
    seed%nat = 0
    seed%nspc = 0
    seed%useabr = 0
    nsline = 0
    if (lu == uin) then
       luout = ucopy
    else
       luout = -1
    endif
    allocate(seed%x(3,10),seed%is(10),seed%spc(2))
    do while (getline(lu,line,ucopy=luout))
       lp = 1
       word = lgetword(line,lp)
       if (equal (word,'cell')) then
          ! cell <a> <b> <c> <alpha> <beta> <gamma>
          ok = eval_next(seed%aa(1),line,lp)
          ok = ok .and. eval_next(seed%aa(2),line,lp)
          ok = ok .and. eval_next(seed%aa(3),line,lp)
          ok = ok .and. eval_next(seed%bb(1),line,lp)
          ok = ok .and. eval_next(seed%bb(2),line,lp)
          ok = ok .and. eval_next(seed%bb(3),line,lp)
          if (.not.ok) then
             call ferror("parse_crystal_env","Wrong cell parameter (CELL) syntax",faterr,line,syntax=.true.)
             return
          endif
          isset = .false.
          word = lgetword(line,lp)
          if (equal(word,'angstrom').or.equal(word,'ang')) then
             isset = .true.
             seed%aa = seed%aa / bohrtoa
          elseif (equal(word,'bohr').or.equal(word,'au')) then
             isset = .true.
          elseif (len_trim(word) > 0) then
             call ferror('parse_crystal_env','Unknown extra keyword in cell parameters (CELL)',faterr,line,syntax=.true.)
             return
          endif
          if (.not.isset) then
             seed%aa = seed%aa / dunit0(iunit0)
          end if
          seed%useabr = 1

          ! cartesian <scale> .. endcartesian
       else if (equal (word,'cartesian')) then
          ok = eval_next(scal,line,lp)
          if (.not.ok) scal = 1d0
          ascal = 1d0/dunit0(iunit)
          aux = getword(line,lp)
          if (len_trim(aux) > 0) then
             call ferror('parse_crystal_env','Unknown extra keyword in lattice vectors (CARTESIAN)',faterr,line,syntax=.true.)
             return
          end if

          i = 0
          rmat = 0d0
          isset = .false.
          do while(.true.)
             lp = 1
             ok = getline(lu,line,ucopy=luout)
             word = lgetword(line,lp)
             if (equal(word,'angstrom') .or.equal(word,'ang')) then
                ! angstrom/ang
                ascal = 1d0/bohrtoa
                isset = .true.
             else if (equal(word,'bohr') .or.equal(word,'au')) then
                ! bohr/au
                ascal = 1d0
                isset = .true.
             else if (equal(word,'end').or.equal(word,'endcartesian')) then
                ! end/endcartesian
                aux = getword(line,lp)
                if (len_trim(aux) > 0) then
                   call ferror('parse_crystal_env','Unknown extra keyword in lattice vectors (CARTESIAN)',faterr,line,syntax=.true.)
                   return
                end if
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
             if (.not.ok) then
                call ferror('parse_crystal_env','Bad lattice vector (CARTESIAN) environment',faterr,line,syntax=.true.)
                return
             end if
             aux = getword(line,lp)
             if (len_trim(aux) > 0) then
                call ferror('parse_crystal_env','Unknown extra keyword in lattice vectors (CARTESIAN)',faterr,line,syntax=.true.)
                return
             end if
          end do
          if (.not.isset) then
             ascal = 1d0 / dunit0(iunit0)
          end if
          seed%m_x2c = transpose(rmat) * scal * ascal
          rmat = seed%m_x2c
          call matinv(rmat,3,ier)
          if (ier /= 0) then
             call ferror('parse_crystal_env','Invalid lattice vectors (matrix not invertible)',faterr)
             return
          end if
          seed%useabr = 2

       else if (equal(word,'spg')) then
          hnum = spg_get_hall_number_from_symbol(line(lp:))
          call spg_get_symmetry_from_database(hnum,seed%neqv,seed%ncv,seed%rotm,seed%cen)
          if (seed%neqv <= 0) then
             call ferror('parse_crystal_env','Error interpreting space group symbol',faterr,line,syntax=.true.)
             return
          end if
          seed%havesym = 1
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
          seed%nat = seed%nat + 1
          if (seed%nat > size(seed%x,2)) then
             call realloc(seed%x,3,2*seed%nat)
             call realloc(seed%is,2*seed%nat)
          end if

          if (.not.equal(word,'neq')) then
             ! try to read four fields from the input
             lp2 = 1
             ok = isinteger(iat,line,lp2)
             ok = ok .and. eval_next(seed%x(1,seed%nat),line,lp2)
             ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp2)
             ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp2)
             if (.not.ok) then
                ! then it must be <atom> <x> <y> <z>
                ok = eval_next(seed%x(1,seed%nat),line,lp)
                ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp)
                ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp)
                if (.not.ok) then
                   call ferror("parse_crystal_env","Wrong atomic input syntax",faterr,line,syntax=.true.)
                   return
                end if
                name = string(word)
             else
                lp = lp2
                name = nameguess(iat,.true.)
             end if
          else
             ok = eval_next(seed%x(1,seed%nat),line,lp)
             ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp)
             ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp)
             if (.not.ok) then
                call ferror("parse_crystal_env","Wrong atom input syntax",faterr,line,syntax=.true.)
                return
             end if
             name = trim(getword(line,lp))
          end if

          it = 0
          do i = 1, seed%nspc
             if (equali(seed%spc(i)%name,name)) then
                it = i
                exit
             end if
          end do
          if (it == 0) then
             seed%nspc = seed%nspc + 1
             if (seed%nspc > size(seed%spc,1)) &
                call realloc(seed%spc,2*seed%nspc)
             it = seed%nspc
             seed%spc(it)%name = name
             seed%spc(it)%z = zatguess(name)
             if (seed%spc(it)%z < 0) then
                call ferror('parse_crystal_env','Unknown atomic symbol in atom input',faterr,line,syntax=.true.)
                return
             end if
             seed%spc(it)%qat = 0d0
          end if
          seed%is(seed%nat) = it

          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'ang') .or. equal(word,'angstrom')) then
                if (seed%useabr /= 2) then
                   call ferror('parse_crystal_env','Need lattice vectors (CARTESIAN) for angstrom coordinates',&
                      faterr,line,syntax=.true.)
                   return
                end if
                seed%x(:,seed%nat) = matmul(rmat,seed%x(:,seed%nat) / bohrtoa)
             else if (equal(word,'bohr') .or. equal(word,'au')) then
                if (seed%useabr /= 2) then
                   call ferror('parse_crystal_env','Need lattice vectors (CARTESIAN) for bohr coordinates',&
                      faterr,line,syntax=.true.)
                   return
                end if
                seed%x(:,seed%nat) = matmul(rmat,seed%x(:,seed%nat))
             else if (len_trim(word) > 0) then
                call ferror('parse_crystal_env','Unknown keyword in atom input',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do
       end if
    end do
    aux = getword(line,lp)
    if (len_trim(aux) > 0) then
       call ferror('parse_crystal_env','Unknown extra keyword in ENDCRYSTAL',faterr,line,syntax=.true.)
       return
    end if
    if (seed%nat == 0) then
       call ferror('parse_crystal_env','No atoms in input',faterr,syntax=.true.)
       return
    end if
    if (seed%useabr == 0) then
       call ferror('parse_crystal_env','No cell information given',faterr,syntax=.true.)
       return
    end if

    ! symm transformation
    if (nsline > 0 .and. allocated(sline)) then
       do i = 1, nsline ! run over symm lines
          rot = string_to_symop(sline(i),errmsg)
          if (len_trim(errmsg) > 0) then
             call ferror('parse_crystal_env','Error parsing SYMM:' // errmsg,faterr,syntax=.true.)
             return
          end if

          do j = 1, seed%nat ! run over atoms
             x = matmul(rot(1:3,1:3),seed%x(:,j)) + rot(:,4)
             x = x - floor(x)

             ! check if this atom already exists
             ok = .true.
             do k = 1, seed%nat
                xdif = x - seed%x(:,k)
                xdif = xdif - nint(xdif)
                if (all(abs(xdif) < 1d-5)) then
                   ok = .false.
                   exit
                endif
             end do

             ! add this atom to the list
             if (ok) then
                seed%nat = seed%nat + 1
                if (seed%nat > size(seed%x,2)) then
                   call realloc(seed%x,3,2*seed%nat)
                   call realloc(seed%is,2*seed%nat)
                end if
                seed%is(seed%nat) = seed%is(j)
                seed%x(:,seed%nat) = x
             endif
          end do
       end do
       deallocate(sline)
    end if
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%spc,seed%nspc)

    ! symmetry
    if (goodspg) then
       seed%havesym = 1
       seed%findsym = 0
       seed%checkrepeats = .false.
    else
       seed%havesym = 0
       seed%findsym = -1
    end if
    oksyn = .true.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = .false.
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = "<input>"
    seed%name = "<input>"
    seed%isformat = isformat_from_input

  end subroutine parse_crystal_env

  !> Parse a molecule environment
  module subroutine parse_molecule_env(seed,lu,oksyn)
    use global, only: rborder_def, eval_next, dunit0, iunit, iunit_ang, iunit_isdef
    use tools_io, only: uin, ucopy, getline, lgetword, equal, ferror, faterr,&
       string, isinteger, nameguess, getword, zatguess, equali
    use param, only: bohrtoa, isformat_from_input
    use types, only: realloc

    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    integer, intent(in) :: lu !< Logical unit for input
    logical, intent(out) :: oksyn !< Was there a syntax error?

    character(len=:), allocatable :: word, aux, line, name
    integer :: lp, lp2, luout, iat, iunit0, it, i
    real*8 :: rborder
    logical :: ok, docube, isset

    call seed%end()
    if (iunit_isdef) then
       iunit0 = iunit_ang
    else
       iunit0 = iunit
    end if

    ok = .false.
    docube = .false.
    rborder = rborder_def
    seed%nat = 0
    seed%nspc = 0
    allocate(seed%x(3,10),seed%is(10),seed%spc(2))
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
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'border')) then
          ! border [border.r]
          ok = eval_next(rborder,line,lp)
          if (.not.ok) then
             call ferror('parse_molecule_input','Wrong syntax in BORDER',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'endmolecule') .or. equal(word,'end')) then
          ! endmolecule/end
          exit
       else
          ! keyword not found, must be an atom. The syntax:
          !    neq <x> <y> <z> <atom> ...
          !    <atom> <x> <y> <z> ...
          !    <atnumber> <x> <y> <z> ...
          ! are acceptable
          seed%nat = seed%nat + 1
          if (seed%nat > size(seed%x,2)) then
             call realloc(seed%x,3,2*seed%nat)
             call realloc(seed%is,2*seed%nat)
          end if

          if (.not.equal(word,'neq')) then
             ! try to read four fields from the input
             lp2 = 1
             ok = isinteger(iat,line,lp2)
             ok = ok .and. eval_next(seed%x(1,seed%nat),line,lp2)
             ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp2)
             ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp2)
             if (.not.ok) then
                ! then it must be <atom> <x> <y> <z>
                ok = eval_next(seed%x(1,seed%nat),line,lp)
                ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp)
                ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp)
                if (.not.ok) then
                   call ferror("parse_molecule_env","Wrong atomic input syntax",faterr,line,syntax=.true.)
                   return
                end if
                name = string(word)
             else
                lp = lp2
                name = nameguess(iat,.true.)
             endif
          else
             ok = eval_next(seed%x(1,seed%nat),line,lp)
             ok = ok .and. eval_next(seed%x(2,seed%nat),line,lp)
             ok = ok .and. eval_next(seed%x(3,seed%nat),line,lp)
             if (.not.ok) then
                call ferror("parse_molecule_env","Wrong atom input syntax",faterr,line,syntax=.true.)
                return
             end if
             name = trim(getword(line,lp))
          endif

          it = 0
          do i = 1, seed%nspc
             if (equali(seed%spc(i)%name,name)) then
                it = i
                exit
             end if
          end do
          if (it == 0) then
             seed%nspc = seed%nspc + 1
             if (seed%nspc > size(seed%spc,1)) &
                call realloc(seed%spc,2*seed%nspc)
             it = seed%nspc
             seed%spc(it)%name = name
             seed%spc(it)%z = zatguess(name)
             if (seed%spc(it)%z < 0) then
                call ferror('parse_molecule_env','Unknown atomic symbol in atom input',faterr,line,syntax=.true.)
                return
             end if
             seed%spc(it)%qat = 0d0
          end if
          seed%is(seed%nat) = it

          isset = .false.
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'ang') .or. equal(word,'angstrom')) then
                isset = .true.
                seed%x(:,seed%nat) = seed%x(:,seed%nat) / bohrtoa
             else if (equal(word,'bohr') .or. equal(word,'au')) then
                isset = .true.
             else if (len_trim(word) > 0) then
                call ferror('parse_molecule_input','Unknown extra keyword in atomic input',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do
          if (.not.isset) then
             seed%x(:,seed%nat) = seed%x(:,seed%nat) / dunit0(iunit0)
          end if
       endif
    end do
    aux = getword(line,lp)
    if (len_trim(aux) > 0) then
       call ferror('parse_molecule_input','Unknown extra keyword in ENDMOLECULE',faterr,line,syntax=.true.)
       return
    end if
    if (seed%nat == 0) then
       call ferror('parse_molecule_input','No atoms in input',faterr,syntax=.true.)
       return
    end if
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%spc,seed%nspc)
    oksyn = .true.
    seed%useabr = 0

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = .true.
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = "<input>"
    seed%name = "<input>"
    seed%isformat = isformat_from_input

  contains
    function check_no_extra_word()
      character(len=:), allocatable :: aux2
      logical :: check_no_extra_word
      aux2 = getword(line,lp)
      if (len_trim(aux2) > 0) then
         call ferror('critic','Unknown extra keyword',faterr,line,syntax=.true.)
         check_no_extra_word = .false.
      else
         check_no_extra_word = .true.
      end if
    end function check_no_extra_word

  end subroutine parse_molecule_env

  !> Read a structure from the critic2 structure library and return
  !> the seed. line is the structure ID. oksyn is true if there were
  !> no syntax errors. If mol is present and true, use the default
  !> molecular library. If mol is present and false, use the default
  !> crystal library. If file is present, use this file as the
  !> library.
  module subroutine read_library(seed,line,oksyn,mol,file,ti)
    use global, only: mlib_file, clib_file
    use tools_io, only: lgetword, ferror, faterr, uout, fopen_read, getline,&
       equal, getword, fclose
    use param, only: isformat_from_library
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: line
    logical, intent(out) :: oksyn
    logical, intent(in), optional :: mol
    character*(*), intent(in), optional :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: word, l2, stru, aux, libfile
    logical :: lchk, found, ok
    integer :: lu, lp, lpo

    ! read the structure
    call seed%end()
    oksyn = .false.
    lpo = 1
    stru = lgetword(line,lpo)
    if (len_trim(stru) < 1) then
       call ferror("read_library","structure label missing in CRYSTAL/MOLECULE LIBRARY",faterr,line,syntax=.true.)
       return
    endif

    libfile = clib_file
    if (present(mol)) then
       if (mol) then
          libfile = mlib_file
       else
          libfile = clib_file
       endif
    elseif (present(file)) then
       libfile = file
    end if

    ! open the library file
    inquire(file=libfile,exist=lchk)
    if (.not.lchk) then
       write (uout,'("(!) Library file:"/"        ",A)') trim(libfile)
       call ferror("read_library","library file not found!",faterr,syntax=.true.)
       return
    endif
    lu = fopen_read(libfile,abspath0=.true.,ti=ti)

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
       write (uout,'("(!) Structure not found in file:"/"        ",A)') trim(libfile)
       call ferror("read_library","structure not found in library!",faterr,syntax=.true.)
       call fclose(lu)
       return
    end if

    ! get the crystal/molecule keyword
    ok = getline(lu,l2)
    if (.not.ok) then
       call fclose(lu)
       call ferror("read_library","error parsing the crystal/molecule environment",faterr,syntax=.true.)
       return
    end if
    lp = 1
    word = lgetword(l2,lp)
    if (equal(word,"crystal")) then
       call seed%parse_crystal_env(lu,ok)
    elseif (equal(word,"molecule")) then
       call seed%parse_molecule_env(lu,ok)
    else
       call fclose(lu)
       call ferror("read_library","error parsing the crystal/molecule keyword",faterr,syntax=.true.)
       return
    endif
    seed%file = trim(line) // " (library)"
    seed%name = trim(line) // " (library)"
    seed%isformat = isformat_from_library

    call fclose(lu)
    if (.not.ok) return

    ! make sure there's no more input
    aux = getword(line,lpo)
    if (len_trim(aux) > 0) then
       call ferror('read_library','Unknown extra keyword in CRYSTAL/MOLECULE LIBRARY',faterr,line,syntax=.true.)
       return
    end if
    oksyn = .true.

  end subroutine read_library

  !> Create a crystal seed from a molecular fragment. If order_by_cidx0
  !> is present and .true., order the resulting seed in increasing value of
  !> cidx from the fragment.
  module subroutine from_fragment(seed,fr,order_by_cidx0)
    use global, only: rborder_def
    use tools, only: qcksort
    use tools_io, only: ferror, faterr
    use param, only: isformat_derived
    class(crystalseed), intent(inout) :: seed
    type(fragment), intent(in) :: fr
    logical, intent(in), optional :: order_by_cidx0

    integer :: i
    logical :: order_by_cidx
    integer, allocatable :: iord(:), midx(:)

    call seed%end()
    if (fr%nat==0) &
       call ferror('from_fragment','fragment has zero atoms',faterr)
    if (.not.fr%discrete) &
       call ferror('from_fragment','cannot handle non-discrete fragments',faterr)
    order_by_cidx = .false.
    if (present(order_by_cidx0)) order_by_cidx = order_by_cidx0

    ! copy species
    seed%nspc = fr%nspc
    if (allocated(seed%spc)) deallocate(seed%spc)
    allocate(seed%spc(seed%nspc))
    seed%spc = fr%spc

    ! determine the final order
    allocate(iord(fr%nat))
    do i = 1, fr%nat
       iord(i) = i
    end do
    if (order_by_cidx) then
       allocate(midx(fr%nat))
       do i = 1, fr%nat
          midx(i) = fr%at(i)%cidx
       end do
       call qcksort(midx,iord,1,fr%nat)
       deallocate(midx)
    end if

    ! copy atoms
    seed%nat = fr%nat
    if (allocated(seed%x)) deallocate(seed%x)
    if (allocated(seed%is)) deallocate(seed%is)
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    do i = 1, seed%nat
       seed%x(:,i) = fr%at(iord(i))%r
       seed%is(i) = fr%at(iord(i))%is
    end do
    deallocate(iord)

    ! rest of the info
    seed%useabr = 0
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.
    seed%isused = .true.
    seed%ismolecule = .true.
    seed%cubic = .false.
    seed%border = rborder_def
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = ""
    seed%name = ""
    seed%isformat = isformat_derived

  end subroutine from_fragment

  !> Remove all hydrogens from the seed
  module subroutine strip_hydrogens(seed)
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed

    integer :: i, nnat

    nnat = 0
    do i = 1, seed%nat
       if (seed%spc(seed%is(i))%z /= 1) then
          nnat = nnat + 1
          seed%x(:,nnat) = seed%x(:,i)
          seed%is(nnat) = seed%is(i)
       end if
    end do
    seed%nat = nnat
    call realloc(seed%x,3,nnat)
    call realloc(seed%is,nnat)

  end subroutine strip_hydrogens

  !> Detect the file format of a file from the extension and read a
  !> crystal seed from it. If mol0 == 1, force a molecule. If mol0 ==
  !> 0, force a crystal. If mol0 == -1, let the format detection
  !> decide. If the read was successful, return an empty error
  !> message
  module subroutine read_any_file(seed,file,mol0,errmsg,ti)
    use global, only: rborder_def
    use param, only: isformat_cif, isformat_shelx, isformat_f21,&
       isformat_cube, isformat_bincube, isformat_struct, isformat_abinit,&
       isformat_elk, isformat_fploout,&
       isformat_qein, isformat_qeout, isformat_crystal, isformat_xyz,&
       isformat_wfn, isformat_wfx, isformat_fchk, isformat_molden,&
       isformat_gaussian, isformat_siesta, isformat_xsf, isformat_gen,&
       isformat_vasp, isformat_pwc, isformat_axsf, isformat_dat,&
       isformat_pgout, isformat_orca, isformat_dmain, isformat_aimsin,&
       isformat_aimsout, isformat_tinkerfrac, isformat_gjf, isformat_zmat,&
       isformat_castepcell, isformat_castepgeom, isformat_mol2,&
       isformat_pdb
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    integer, intent(in) :: mol0
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: isformat
    logical :: mol, hastypes

    ! detect the format for this file
    errmsg = ""
    call struct_detect_format(file,isformat,ti=ti)

    ! is this a crystal or a molecule?
    if (mol0 == 1) then
       mol = .true.
    elseif (mol0 == 0) then
       mol = .false.
    elseif (mol0 == -1) then
       call struct_detect_ismol(file,isformat,mol,ti=ti)
    else
       errmsg = "read_any_file: unknown mol0"
       return
    end if

    ! build the seed
    if (isformat == isformat_cif) then
       call seed%read_cif(file," ",mol,errmsg,ti=ti)

    elseif (isformat == isformat_mol2) then
       call seed%read_mol2(file,rborder_def,.false.," ",errmsg,ti=ti)

    elseif (isformat == isformat_shelx) then
       call seed%read_shelx(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_f21) then
       call seed%read_f21(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_cube) then
       call seed%read_cube(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_bincube) then
       call seed%read_bincube(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_struct) then
       call seed%read_wien(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_vasp) then
       call seed%read_vasp(file,mol,hastypes,errmsg,ti=ti)
       if (len_trim(errmsg) == 0 .and..not.hastypes) then
          errmsg = "Atom types not found in VASP file"
          return
       end if

    elseif (isformat == isformat_abinit) then
       call seed%read_abinit(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_elk) then
       call seed%read_elk(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_qeout) then
       call seed%read_qeout(file,mol,0,errmsg,ti=ti)

    elseif (isformat == isformat_crystal) then
       call seed%read_crystalout(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_fploout) then
       call seed%read_fploout(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_qein) then
       call seed%read_qein(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_xyz.or.isformat == isformat_wfn.or.&
       isformat == isformat_wfx.or.isformat == isformat_fchk.or.&
       isformat == isformat_molden.or.isformat == isformat_gaussian.or.&
       isformat == isformat_dat.or.isformat == isformat_pgout.or.&
       isformat == isformat_orca.or.isformat == isformat_gjf.or.&
       isformat == isformat_zmat.or.isformat == isformat_pdb) then
       call seed%read_mol(file,isformat,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_siesta) then
       call seed%read_siesta(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_castepcell) then
       call seed%read_castep_cell(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_castepgeom) then
       call seed%read_castep_geom(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_xsf) then
       call seed%read_xsf(file,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_gen) then
       call seed%read_dftbp(file,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_pwc) then
       call seed%read_pwc(file,mol,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_axsf) then
       call seed%read_axsf(file,1,0d0,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_dmain) then
       call seed%read_dmain(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_aimsin) then
       call seed%read_aimsin(file,mol,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_aimsout) then
       call seed%read_aimsout(file,mol,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_tinkerfrac) then
       call seed%read_tinkerfrac(file,mol,errmsg,ti=ti)

    else
       errmsg = "unrecognized file format"
       return
    end if

  end subroutine read_any_file

  !> Read a structure seed from a CIF file. If dblock is empty, return
  !> the first structure in the cif, otherwise return the seed with
  !> data name equal to dblock. If mol, force a molecule/crystal
  !> system.  If error, return non-empty errmsg.
  module subroutine read_cif(seed,file,dblock,mol,errmsg,ti)
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    character*(*), intent(in) :: dblock !< Data block
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    call seed%end()
    call read_all_cif(file,mol,errmsg,seed0=seed,dblock=dblock,ti=ti)

  end subroutine read_cif

  !> Read a structure seed from a TRIPOS mol2 file. If name is empty,
  !> return the first molecule in the mol2 file, otherwise return the
  !> molecule with that name. If mol, force a molecule/crystal system.
  !> If error, return non-empty errmsg.
  module subroutine read_mol2(seed,file,rborder,docube,name,errmsg,ti)
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character*(*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    call seed%end()
    call read_all_mol2(file,errmsg,seed0=seed,name=name,ti=ti)
    seed%border = rborder
    seed%cubic = docube

  end subroutine read_mol2

  !> Read the structure from a res or ins (shelx) file
  module subroutine read_shelx(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, isreal, isinteger,&
       lower, zatguess, fclose
    use param, only: eyet, eye, bohrtoa, isformat_shelx
    use types, only: realloc
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, ilat
    logical :: ok, iscent, havecell, found
    character(len=:), allocatable :: word, line, aux
    real*8 :: raux, rot0(3,4)
    integer :: i, j, n
    integer :: iz
    integer :: lncv
    real*8, allocatable :: lcen(:,:)

    real*8, parameter :: eps = 1d-5

    ! file and seed name
    call seed%end()
    seed%file = file
    seed%name = file
    seed%isformat = isformat_shelx
    errmsg = ""

    ! initialize symmetry
    iscent = .false.
    seed%ncv = 1
    allocate(seed%cen(3,4))
    seed%cen = 0d0
    seed%neqv = 1
    allocate(seed%rotm(3,4,48))
    seed%rotm = 0d0
    seed%rotm(:,:,seed%neqv) = eyet
    havecell = .false.
    seed%nat = 0
    if (.not.allocated(seed%x)) allocate(seed%x(3,10))
    if (.not.allocated(seed%is)) allocate(seed%is(10))

    ! centering vectors may come in symm. If that happens,
    ! replicate the atoms and let LATT determine the global
    ! centering vectors
    lncv = 0
    allocate(lcen(3,1))
    lcen = 0d0

    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    do while (.true.)
       ok = getline_local()
       if (.not.ok) exit
       lp = 1
       word = lgetword(line,lp)
       if (len_trim(word) > 4) word = word(1:4)
       if (equal(word,"cell")) then
          ! read the cell parameters from the cell card
          ok = isreal(raux,line,lp)
          ok = ok .and. isreal(seed%aa(1),line,lp)
          ok = ok .and. isreal(seed%aa(2),line,lp)
          ok = ok .and. isreal(seed%aa(3),line,lp)
          ok = ok .and. isreal(seed%bb(1),line,lp)
          ok = ok .and. isreal(seed%bb(2),line,lp)
          ok = ok .and. isreal(seed%bb(3),line,lp)
          if (.not.ok) then
             errmsg = "Error reading CELL card."
             goto 999
          end if
          seed%aa = seed%aa / bohrtoa
          havecell = .true.
       elseif (equal(word,"latt")) then
          ! read the centering vectors from the latt card
          ok = isinteger(ilat,line,lp)
          if (.not.ok) then
             errmsg = "Error reading LATT card."
             goto 999
          end if
          select case(abs(ilat))
          case(1)
             ! P
             seed%ncv=1
          case(2)
             ! I
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
             seed%cen(3,2)=0.5d0
          case(3)
             ! R obverse
             seed%ncv=3
             seed%cen(:,2) = (/2d0,1d0,1d0/) / 3d0
             seed%cen(:,3) = (/1d0,2d0,2d0/) / 3d0
          case(4)
             ! F
             seed%ncv=4
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
             seed%cen(2,3)=0.5d0
             seed%cen(3,3)=0.5d0
             seed%cen(1,4)=0.5d0
             seed%cen(3,4)=0.5d0
          case(5)
             ! A
             seed%ncv=2
             seed%cen(2,2)=0.5d0
             seed%cen(3,2)=0.5d0
          case(6)
             ! B
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(3,2)=0.5d0
          case(7)
             ! C
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
          case default
             errmsg = "Unknown LATT value."
             goto 999
          end select
          iscent = (ilat > 0)
       elseif (equal(word,"symm")) then
          ! symmetry operations from the symm card
          aux = lower(line(lp:)) // ","
          rot0 = string_to_symop(aux,errmsg)
          if (len_trim(errmsg) > 0) goto 999

          if (all(abs(eye - rot0(1:3,1:3)) < 1d-12)) then
             ! a non-zero pure translation or the identity
             if (all(abs(rot0(:,4)) < 1d-12)) then
                ! ignore the identity
             else
                ! must be a pure translation
                lncv = lncv + 1
                if (lncv > size(lcen,2)) &
                   call realloc(lcen,3,2*lncv)
                lcen(:,lncv) = rot0(:,4)
             endif
          else
             ! a rotation, with some pure translation in it
             ! check if I have this rotation matrix already
             ok = .true.
             do i = 1, seed%neqv
                if (all(abs(seed%rotm(:,:,i) - rot0(:,:)) < 1d-12)) then
                   ok = .false.
                   exit
                endif
             end do
             if (ok) then
                seed%neqv = seed%neqv + 1
                seed%rotm(:,:,seed%neqv) = rot0
             else
                errmsg = "Found repeated entry in SYMM."
                goto 999
             endif
          endif

       elseif (equal(word,"sfac")) then
          ! atomic types from the sfac card
          seed%nspc = 0
          allocate(seed%spc(2))
          do while (.true.)
             word = lgetword(line,lp)
             iz = zatguess(word)
             if (iz <= 0 .or. len_trim(word) < 1) exit
             seed%nspc = seed%nspc + 1
             if (seed%nspc > size(seed%spc,1)) call realloc(seed%spc,2*seed%nspc)
             seed%spc(seed%nspc)%z = iz
             seed%spc(seed%nspc)%name = trim(word)
          end do
          call realloc(seed%spc,seed%nspc)
       elseif (equal(word,"unit")) then
          ! ignore the unit card... some res files don't have it,
          ! and we can count the atoms in the list anyway

          ! ignore all the following cards
       elseif (equal(word,"abin").or.equal(word,"acta").or.equal(word,"afix").or.&
          equal(word,"anis").or.equal(word,"ansc").or.equal(word,"ansr").or.&
          equal(word,"basf").or.equal(word,"bind").or.equal(word,"bloc").or.&
          equal(word,"bond").or.equal(word,"bump").or.equal(word,"cgls").or.&
          equal(word,"chiv").or.equal(word,"conf").or.equal(word,"conn").or.&
          equal(word,"damp").or.equal(word,"dang").or.equal(word,"defs").or.&
          equal(word,"delu").or.equal(word,"dfix").or.equal(word,"disp").or.&
          equal(word,"eadp").or.equal(word,"eqiv").or.equal(word,"exti").or.&
          equal(word,"exyz").or.equal(word,"flat").or.equal(word,"fmap").or.&
          equal(word,"free").or.equal(word,"fvar").or.equal(word,"grid").or.&
          equal(word,"hfix").or.equal(word,"hklf").or.equal(word,"hope").or.&
          equal(word,"htab").or.&
          equal(word,"isor").or.equal(word,"laue").or.equal(word,"list").or.&
          equal(word,"l.s.").or.equal(word,"merg").or.equal(word,"mole").or.&
          equal(word,"more").or.&
          equal(word,"move").or.equal(word,"mpla").or.equal(word,"ncsy").or.&
          equal(word,"neut").or.equal(word,"omit").or.equal(word,"part").or.&
          equal(word,"plan").or.equal(word,"prig").or.equal(word,"rem").or.&
          equal(word,"resi").or.equal(word,"rigu").or.equal(word,"rtab").or.&
          equal(word,"sadi").or.equal(word,"same").or.equal(word,"shel").or.&
          equal(word,"simu").or.equal(word,"size").or.equal(word,"spec").or.&
          equal(word,"stir").or.equal(word,"sump").or.equal(word,"swat").or.&
          equal(word,"temp").or.equal(word,"time").or.&
          equal(word,"titl").or.equal(word,"twin").or.&
          equal(word,"twst").or.equal(word,"wght").or.equal(word,"wigl").or.&
          equal(word,"wpdb").or.equal(word,"xnpd").or.equal(word,"zerr")) then
          cycle

          ! also ignore the frag...fend blocks
       elseif (equal(word,"frag")) then
          do while (.true.)
             ok = getline_local()
             if (.not.ok) then
                errmsg = "Unexpected end of file inside frag block."
                goto 999
             end if
             lp = 1
             word = lgetword(line,lp)
             if (equal(word,"fend")) exit
          end do
       elseif (equal(word,"end")) then
          ! end of the input
          exit
       else
          ! maybe this is an atom, but if we can not tell, it could be a new or very old keyword
          iz = zatguess(word)
          if (iz < 0) cycle

          ! check if this is an atom
          seed%nat = seed%nat + 1
          if (seed%nat > size(seed%is)) then
             call realloc(seed%x,3,2*seed%nat)
             call realloc(seed%is,2*seed%nat)
          end if
          ok = isinteger(iz,line,lp)
          ok = ok .and. isreal(seed%x(1,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(2,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(3,seed%nat),line,lp)
          if (.not.ok) then
             seed%nat = seed%nat - 1
             continue
          end if
          if (iz < 1 .or. iz > seed%nspc) then
             errmsg = "Atom type not found in SFAC list."
             goto 999
          end if
          seed%is(seed%nat) = iz
       end if
    end do

    if (seed%nspc == 0) then
       errmsg = "No SFAC information (atomic types) found."
       goto 999
    end if
    if (seed%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    if (.not.havecell) then
       errmsg = "No cell found."
       goto 999
    end if
    seed%useabr = 1 ! use aa and bb

    if (iscent) then
       ! add the -1 operations
       n = seed%neqv
       do i = 1, n
          seed%rotm(1:3,1:3,n+i) = -seed%rotm(1:3,1:3,i)
          seed%rotm(:,4,n+i) = seed%rotm(:,4,i)
       end do
       seed%neqv = 2*n
    end if

    ! Replicate the atoms using the local centering vectors passed in
    ! SYMM, if there are any. This is contrary to the SHELX res format
    ! specs, but some programs do it anyway.
    if (lncv > 0) then
       do i = 1, lncv
          found = .false.
          do j = 1, seed%ncv
             if (all(lcen(:,i) - seed%cen(:,j) < eps)) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) then
             seed%ncv = seed%ncv + 1
             if (seed%ncv > size(seed%cen,2)) &
                call realloc(seed%cen,3,2*seed%ncv)
             seed%cen(:,seed%ncv) = lcen(:,i)
          end if
       end do
    end if
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%cen,3,seed%ncv)

    ! use the symmetry in this file
    seed%havesym = 1
    seed%checkrepeats = .true.
    seed%findsym = -1
    call realloc(seed%rotm,3,4,seed%neqv)
    call realloc(seed%cen,3,seed%ncv)

999 continue
    call fclose(lu)

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0

  contains
    function getline_local() result(ok_)
      logical :: ok_
      integer :: idx
      character(len=:), allocatable :: aux

      ok_ = getline_raw(lu,line,.false.)
      if (.not.ok_) return
      idx = index(line,"=",.true.)
      do while (idx == len_trim(line))
         ok_ = getline_raw(lu,aux,.false.)
         line = line(1:idx-1) // trim(aux)
         if (.not.ok_) return
         idx = index(line,"=",.true.)
      end do
    end function getline_local
  end subroutine read_shelx

  !> Read the structure from a fort.21 from neighcrys
  module subroutine read_f21(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, nameguess
    use param, only: bohrtoa, maxzat0, mlen, isformat_f21
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    character(len=:), allocatable :: line
    character(len=mlen) :: dum
    character*15 :: str
    character*1 :: latt
    logical :: ok
    integer :: isused(maxzat0), idum, iz, i
    real*8 :: rdum

    ! open the file for reading
    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) goto 999
    seed%file = file
    seed%name = file
    seed%isformat = isformat_f21

    ! read the lattice constants and the atomic coordinates
    isused = 0
    seed%useabr = 1
    seed%nat = 0
    seed%nspc = 0
    seed%neqv = 0
    seed%ncv = 0
    allocate(seed%spc(10),seed%is(10),seed%rotm(3,4,48),seed%cen(3,4))
    seed%cen = 0d0
    seed%rotm = 0d0
    do while(getline_raw(lu,line))
       if (index(line,"LATTICE CENTRING TYPE") > 0) then
          read (line,*,err=999,end=999) dum, dum, dum, latt
          if (latt == "P") then
             seed%ncv=1
          elseif (latt == "I") then
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
             seed%cen(3,2)=0.5d0
          elseif (latt == "R") then
             seed%ncv=3
             seed%cen(:,2) = (/2d0,1d0,1d0/) / 3d0
             seed%cen(:,3) = (/1d0,2d0,2d0/) / 3d0
          elseif (latt == "F") then
             seed%ncv=4
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
             seed%cen(2,3)=0.5d0
             seed%cen(3,3)=0.5d0
             seed%cen(1,4)=0.5d0
             seed%cen(3,4)=0.5d0
          elseif (latt == "A") then
             seed%ncv=2
             seed%cen(2,2)=0.5d0
             seed%cen(3,2)=0.5d0
          elseif (latt == "B") then
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(3,2)=0.5d0
          elseif (latt == "C") then
             seed%ncv=2
             seed%cen(1,2)=0.5d0
             seed%cen(2,2)=0.5d0
          else
             errmsg = "unknown lattice centring type"
             goto 999
          end if
       elseif (index(line,"ROTATION MATRICES") > 0) then
          ok = getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          if (.not.ok) goto 999

          do while (.true.)
             ok = getline_raw(lu,line)
             ok = ok.and.getline_raw(lu,line)
             if (.not.ok) goto 999
             if (len_trim(line) == 0) exit

             seed%neqv = seed%neqv + 1
             read(line,*,err=999,end=999) (rdum,i=1,4), seed%rotm(1,:,seed%neqv)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read(line,*,err=999,end=999) (rdum,i=1,4), seed%rotm(2,:,seed%neqv)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read(line,*,err=999,end=999) (rdum,i=1,4), seed%rotm(3,:,seed%neqv)
          end do
       elseif (index(line,"LATTICE CONSTANTS") > 0) then
          read (lu,*,err=999,end=999) dum, dum, seed%aa(1), dum, dum, seed%aa(2), dum, dum, seed%aa(3)
          read (lu,*,err=999,end=999) dum, dum, seed%bb(1), dum, dum, seed%bb(2), dum, dum, seed%bb(3)
          seed%aa = seed%aa / bohrtoa

       elseif (index(line,"Inequivalent basis atoms") > 0) then

          ! read the first block, with the atomic numbers and indices
          ok = getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          if (.not.ok) goto 999
          do while(getline_raw(lu,line))
             if (len_trim(line) == 0) exit
             ! a new atom
             seed%nat = seed%nat + 1
             read (line,*,err=999,end=999) idum, iz, str

             ! assign atomic species by atomic number
             if (isused(iz) == 0) then
                seed%nspc = seed%nspc + 1
                if (seed%nspc > size(seed%spc,1)) &
                   call realloc(seed%spc,2*seed%nspc)
                seed%spc(seed%nspc)%z = iz
                seed%spc(seed%nspc)%name = nameguess(iz,.true.)
                isused(iz) = seed%nspc
             end if

             ! assign atom type
             if (seed%nat > size(seed%is,1)) &
                call realloc(seed%is,2*seed%nat)
             seed%is(seed%nat) = isused(iz)
          end do

          ! read the second block, with the atomic coordinates
          allocate(seed%x(3,seed%nat))
          ok = getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          do i = 1, seed%nat
             ok = getline_raw(lu,line)
             ok = ok.and.getline_raw(lu,line)
             if (.not.ok) goto 999
             read (line,*,err=999,end=999) seed%x(:,i)
          end do
       end if
    end do
    if (seed%nat == 0 .or. seed%nspc == 0 .or. seed%neqv == 0 .or. seed%ncv == 0) goto 999
    call realloc(seed%is,seed%nat)
    call realloc(seed%spc,seed%nspc)
    call realloc(seed%rotm,3,4,seed%neqv)
    call realloc(seed%cen,3,seed%ncv)

    ! close the file and clean up
    call fclose(lu)

    ! have symmetry, but recalculate
    seed%havesym = 1
    seed%checkrepeats = .false.
    seed%findsym = -1

    ! molecule
    seed%ismolecule = mol
    seed%havex0 = .true.
    seed%molx0 = 0d0

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = .false.
    seed%border = 0d0

    errmsg = ""
999 continue
    if (lu > 0) call fclose(lu)

  end subroutine read_f21

  !> Read the structure from a gaussian cube file
  module subroutine read_cube(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, nameguess, getline_raw
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: isformat_cube
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    integer :: i, j, nstep(3), nn, iz, it, ier
    real*8 :: x0(3), rmat(3,3), rdum, rx(3), rxt(3)
    logical :: ismo, ok
    character(len=:), allocatable :: line

    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! the name of the seed is the first line
    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    seed%file = file
    seed%name = file
    seed%isformat = isformat_cube

    ! ignore the title lines
    read (lu,*,err=999,end=999)

    ! number of atoms and unit cell
    read (lu,*,err=999,end=999) seed%nat, x0
    ismo = (seed%nat < 0)
    seed%nat = abs(seed%nat)

    do i = 1, 3
       read (lu,*,err=999,end=999) nstep(i), rmat(:,i)
       rmat(:,i) = rmat(:,i) * nstep(i)
    end do

    seed%m_x2c = rmat
    rmat = transpose(rmat)
    call matinv(rmat,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if

    seed%useabr = 2

    ! Atomic positions.
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    allocate(seed%spc(2))
    nn = seed%nat
    seed%nat = 0
    do i = 1, nn
       read (lu,*,err=999,end=999) iz, rdum, rx
       if (iz > 0) then
          seed%nat = seed%nat + 1
          !! intel compiler errors with
          ! rx = matmul(rx - x0,rmat)
          rxt = rx - x0
          rx = matmul(rxt,rmat)
          seed%x(:,seed%nat) = rx - floor(rx)
          it = 0
          do j = 1, seed%nspc
             if (seed%spc(j)%z == iz) then
                it = j
                exit
             end if
          end do
          if (it == 0) then
             seed%nspc = seed%nspc + 1
             if (seed%nspc > size(seed%spc,1)) &
                call realloc(seed%spc,2*seed%nspc)
             seed%spc(seed%nspc)%z = iz
             seed%spc(seed%nspc)%name = nameguess(iz)
             it = seed%nspc
          end if
          seed%is(seed%nat) = it
       endif
    end do
    if (seed%nat /= nn) then
       call realloc(seed%x,3,seed%nat)
       call realloc(seed%is,seed%nat)
    end if
    call realloc(seed%spc,seed%nspc)

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%checkrepeats = .false.
    seed%findsym = -1

    ! molecule
    seed%ismolecule = mol
    seed%havex0 = .true.
    seed%molx0 = x0

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = .false.
    seed%border = 0d0

  end subroutine read_cube

  !> Read the structure from a binary cube file
  module subroutine read_bincube(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, nameguess, getline_raw
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: isformat_bincube
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ier
    integer :: i, j, nstep(3), nn, iz, it
    real*8 :: x0(3), rmat(3,3), rdum, rx(3), rxt(3)

    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,form="unformatted",ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! number of atoms and unit cell
    read (lu,err=999,end=999) seed%nat, x0

    read (lu,err=999,end=999) nstep, rmat
    do i = 1, 3
       rmat(:,i) = rmat(:,i) * nstep(i)
    end do

    seed%m_x2c = rmat
    rmat = transpose(rmat)
    call matinv(rmat,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if
    seed%useabr = 2

    ! Atomic positions.
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    allocate(seed%spc(2))
    nn = seed%nat
    seed%nat = 0
    do i = 1, nn
       read (lu,err=999,end=999) iz, rdum, rx
       if (iz > 0) then
          seed%nat = seed%nat + 1
          !! intel compiler errors with
          ! rx = matmul(rx - x0,rmat)
          rxt = rx - x0
          rx = matmul(rxt,rmat)
          seed%x(:,seed%nat) = rx - floor(rx)
          it = 0
          do j = 1, seed%nspc
             if (seed%spc(j)%z == iz) then
                it = j
                exit
             end if
          end do
          if (it == 0) then
             seed%nspc = seed%nspc + 1
             if (seed%nspc > size(seed%spc,1)) &
                call realloc(seed%spc,2*seed%nspc)
             seed%spc(seed%nspc)%z = iz
             seed%spc(seed%nspc)%name = nameguess(iz)
             it = seed%nspc
          end if
          seed%is(seed%nat) = it
       endif
    end do
    if (seed%nat /= nn) then
       call realloc(seed%x,3,seed%nat)
       call realloc(seed%is,seed%nat)
    end if
    call realloc(seed%spc,seed%nspc)

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%checkrepeats = .false.
    seed%findsym = -1

    ! molecule
    seed%ismolecule = mol
    seed%havex0 = .true.
    seed%molx0 = x0

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = .false.
    seed%border = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_bincube

  end subroutine read_bincube

  !> Read the crystal structure from a WIEN2k STRUCT file.
  !> Code adapted from the WIEN2k distribution.
  module subroutine read_wien(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, zatguess, fclose, equal, equali
    use types, only: realloc
    use param, only: pi, isformat_struct
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< struct file
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lut
    integer :: i, j, i1, i2, j1, iat, iat0, istart, it
    real*8 :: mat(3,3), rnot, rmt, pos(3), tau(3), znuc
    integer :: multw, iatnr, iz(3,3), jatom, mu, jri
    character*4 :: lattic, cform
    character*80 :: titel
    character*10 :: aname
    logical :: readall
    real*8 :: ahex, chex

    ! seed file
    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    seed%file = file
    seed%isformat = isformat_struct

    ! first pass to see whether we have symmetry or not
    lut = fopen_read(file,ti=ti)
    if (lut < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    READ(lut,102,err=999,end=999) TITEL
    READ(lut,103,err=999,end=999) LATTIC, seed%nat, cform
    READ(lut,100,err=999,end=999) seed%aa(1:3), seed%bb(1:3)
    DO JATOM=1,seed%nat
       READ(lut,1012,err=999,end=999) iatnr,pos,MULTW
       DO MU=1,MULTW-1
          READ(lut,1013,err=999,end=999) iatnr, pos
       end DO
       READ(lut,113,err=999,end=999) ANAME,JRI,RNOT,RMT,Znuc
       READ(lut,1051,err=999,end=999) ((mat(I1,J1),I1=1,3),J1=1,3)
    end DO
    READ(lut,114,err=999,end=999) seed%neqv

    readall = (seed%neqv <= 0)

    ! second pass -> actually process the information
    rewind(lut)
    READ(lut,102,err=999,end=999) TITEL
    seed%name = file

    READ(lut,103,err=999,end=999) LATTIC, seed%nat, cform
102 FORMAT(A80)
103 FORMAT(A4,23X,I3,1x,a4,/,4X,4X) ! new

    seed%ncv = 1
    allocate(seed%cen(3,4))
    seed%cen = 0d0
    IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
       seed%ncv=1
    ELSE IF(LATTIC(1:1).EQ.'F') THEN
       seed%ncv=4
       seed%cen(1,2)=0.5d0
       seed%cen(2,2)=0.5d0
       seed%cen(2,3)=0.5d0
       seed%cen(3,3)=0.5d0
       seed%cen(1,4)=0.5d0
       seed%cen(3,4)=0.5d0
    ELSE IF(LATTIC(1:1).EQ.'B') THEN
       seed%ncv=2
       seed%cen(1,2)=0.5d0
       seed%cen(2,2)=0.5d0
       seed%cen(3,2)=0.5d0
    ELSE IF(LATTIC(1:1).EQ.'H') THEN
       seed%ncv=1
    ELSE IF(LATTIC(1:1).EQ.'R') THEN
       seed%ncv=1
    ELSE IF(LATTIC(1:3).EQ.'CXY') THEN
       seed%ncv=2
       seed%cen(1,2)=0.5d0
       seed%cen(2,2)=0.5d0
    ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN
       seed%ncv=2
       seed%cen(2,2)=0.5d0
       seed%cen(3,2)=0.5d0
    ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN
       seed%ncv=2
       seed%cen(1,2)=0.5d0
       seed%cen(3,2)=0.5d0
    ELSE
       errmsg = "Unknown lattice."
       goto 999
    END IF

    READ(lut,100,err=999,end=999) seed%aa(1:3), seed%bb(1:3)

    if (LATTIC(1:1) == 'R') then
       ahex = seed%aa(1)
       chex = seed%aa(3)
       seed%aa = sqrt((chex/3d0)**2 + ahex**2/3d0)
       seed%bb = 180d0 * (1d0 - 2d0 * acos(ahex/2d0/seed%aa(1)) / pi)
    endif

100 FORMAT(6F10.5)
    if(seed%bb(3) == 0.d0) seed%bb(3)=90.d0
    seed%useabr = 1

    seed%nspc = 0
    allocate(seed%spc(2))
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    iat = 0
    DO JATOM=1,seed%nat
       iat0 = iat
       iat = iat + 1
       if (iat > size(seed%is)) then
          call realloc(seed%x,3,2*iat)
          call realloc(seed%is,2*iat)
       end if
       READ(lut,1012,err=999,end=999) iatnr,seed%x(:,iat),MULTW

       istart = iat
       if (readall) then
          DO MU=1,MULTW-1
             iat = iat + 1
             if (iat > size(seed%is)) then
                call realloc(seed%x,3,2*iat)
                call realloc(seed%is,2*iat)
             end if
             READ(lut,1013,err=999,end=999) iatnr, seed%x(:,iat)
          end DO
       else
          DO MU=1,MULTW-1
             READ(lut,1013,err=999,end=999) iatnr, pos
          end DO
       end if

       READ(lut,113,err=999,end=999) ANAME,JRI,RNOT,RMT,Znuc
       aname = adjustl(aname)
       it = 0
       do i = 1, seed%nspc
          if (equali(aname,seed%spc(i)%name)) then
             it = i
             exit
          end if
       end do
       if (it == 0) then
          seed%nspc = seed%nspc + 1
          if (seed%nspc > size(seed%spc,1)) &
             call realloc(seed%spc,2*seed%nspc)
          seed%spc(seed%nspc)%name = aname
          seed%spc(seed%nspc)%z = zatguess(aname)
          it = seed%nspc
       end if
       do i = iat0+1, iat
          seed%is(i) = it
       end do
       READ(lut,1051,err=999,end=999) ((mat(I1,J1),I1=1,3),J1=1,3)
    end DO
113 FORMAT(A10,5X,I5,5X,F10.5,5X,F10.5,5X,F5.2)
1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/15X,I2) ! new
1013 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7) ! new
1051 FORMAT(20X,3F10.8)
    seed%nat = iat
    call realloc(seed%x,3,iat)
    call realloc(seed%is,iat)
    call realloc(seed%spc,seed%nspc)

    !.read number of symmetry operations, sym. operations
    READ(lut,114,err=999,end=999) seed%neqv
114 FORMAT(I4)

    if (seed%neqv > 0) then
       allocate(seed%rotm(3,4,seed%neqv))
       do i=1, seed%neqv
          read(lut,115,err=999,end=999) ((iz(i1,i2),i1=1,3),tau(i2),i2=1,3)
          do j=1,3
             seed%rotm(:,j,i)=dble(iz(j,:))
          enddo
          seed%rotm(:,4,i)=tau
       end do
    end if

115 FORMAT(3(3I2,F10.5,/))

    ! symmetry
    if (seed%neqv > 0) then
       seed%havesym = 1
       seed%findsym = 0
    else
       seed%havesym = 0
       seed%findsym = -1
    end if
    seed%checkrepeats = .false.

    errmsg = ""
999 continue

    ! clean up
    call fclose(lut)

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0

  end subroutine read_wien

  !> Read everything except the grid from a VASP POSCAR, etc. file
  !> If hastypes is present, it is equal to .true. if the file
  !> could be read successfully or .false. if the atomic types
  !> are missing.
  module subroutine read_vasp(seed,file,mol,hastypes,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, getline_raw, isreal, &
       getword, zatguess, string, isinteger, nameguess, fclose
    use tools_math, only: det3sym, matinv
    use param, only: bohrtoa, isformat_vasp
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    logical, intent(out) :: hastypes
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nn
    character(len=:), allocatable :: word, line
    logical :: ok, iscar

    integer :: i, j, ier
    real*8 :: scalex, scaley, scalez, scale
    real*8 :: rprim(3,3), gprim(3,3)
    real*8 :: omegaa

    ! open
    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    hastypes = .true.
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

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
       omegaa = sqrt(det3sym(gprim))
       ! adjust the lengths to give the volume
       scale = (abs(scale) / abs(omegaa))**(1d0/3d0)
    end if
    rprim(1,:) = rprim(1,:) * scalex * scale
    rprim(2,:) = rprim(2,:) * scaley * scale
    rprim(3,:) = rprim(3,:) * scalez * scale
    rprim = rprim / bohrtoa
    gprim = matmul(transpose(rprim),rprim)
    omegaa = sqrt(det3sym(gprim))
    if (omegaa < 0d0) then
       errmsg = "Negative cell volume."
       goto 999
    end if
    seed%m_x2c = rprim
    call matinv(rprim,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if
    seed%useabr = 2

    ! For versions >= 5.2, a line indicating the atom types appears here
    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    lp = 1
    word = getword(line,lp)
    if (zatguess(word) >= 0) then
       ! An atom name has been read -> read the rest of the line
       seed%nspc = 0
       if (allocated(seed%spc)) deallocate(seed%spc)
       allocate(seed%spc(2))
       do while (zatguess(word) >= 0)
          seed%nspc = seed%nspc + 1
          if (seed%nspc > size(seed%spc,1)) &
             call realloc(seed%spc,2*seed%nspc)
          seed%spc(seed%nspc)%name = word
          seed%spc(seed%nspc)%z = zatguess(word)
          word = getword(line,lp)
       end do
       call realloc(seed%spc,seed%nspc)
       ok = getline_raw(lu,line,.true.)
       if (.not.ok) goto 999
       hastypes = .true.
    else
       errmsg = ""
       hastypes = .false.
    end if

    ! read number of atoms of each type
    lp = 1
    seed%nat = 0
    i = 0
    allocate(seed%is(10))
    do while(isinteger(nn,line,lp))
       i = i + 1
       do j = seed%nat+1, seed%nat+nn
          if (j > size(seed%is)) &
             call realloc(seed%is,2*(seed%nat+nn))
          seed%is(j) = i
       end do
       seed%nat = seed%nat + nn
    end do
    if (hastypes) then
       if (i /= seed%nspc) then
          errmsg = "Too many atom types"
          goto 999
       end if
    end if
    allocate(seed%x(3,seed%nat))
    call realloc(seed%is,seed%nat)

    ! check there are no more atoms in this line
    nn = -1
    ok = isinteger(nn,line,lp)
    if (ok .and. nn /= -1) then
       errmsg = "Too few atom types"
       goto 999
    end if

    ! Read atomic positions (cryst. coords.)
    read(lu,*,err=999,end=999) line
    line = adjustl(line)
    if (line(1:1) == 's' .or. line(1:1) == 'S') then
       read(lu,*,err=999,end=999) line
       line = adjustl(line)
    endif
    iscar = .false.
    if (line(1:1) == 'd' .or. line(1:1) == 'D') then
       iscar = .false.
    elseif (line(1:1) == 'c' .or. line(1:1) == 'C' .or. line(1:1) == 'k' .or. line(1:1) == 'K') then
       iscar = .true.
    endif
    do i = 1, seed%nat
       read(lu,*,err=999,end=999) seed%x(:,i)
       if (iscar) &
          seed%x(:,i) = matmul(rprim,seed%x(:,i) / bohrtoa)
    enddo

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_vasp

  end subroutine read_vasp

  !> Read the structure from an abinit DEN file (and similar files: ELF, LDEN, etc.)
  module subroutine read_abinit(seed,file,mol,errmsg,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, nameguess, fclose
    use abinit_private, only: hdr_type, hdr_io
    use types, only: realloc
    use param, only: isformat_abinit
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, fform0
    type(hdr_type) :: hdr
    integer :: i, iz
    real*8 :: rmat(3,3)

    call seed%end()
    errmsg = ""
    lu = fopen_read(file,"unformatted",errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! read the header of the DEN file
    call hdr_io(fform0,hdr,1,lu,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! cell parameters
    rmat = hdr%rprimd(:,:)
    seed%m_x2c = rmat
    seed%useabr = 2

    ! types
    seed%nspc = hdr%ntypat
    allocate(seed%spc(seed%nspc))
    do i = 1, seed%nspc
       iz = nint(hdr%znucltypat(i))
       seed%spc(i)%z = iz
       seed%spc(i)%name = nameguess(iz)
    end do

    ! atoms
    seed%nat = hdr%natom
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    do i = 1, seed%nat
       seed%x(:,i) = hdr%xred(:,i)
       seed%is(i) = hdr%typat(i)
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! abinit has symmetry in hdr%nsym/hdr%symrel, but there is no
    ! distinction between pure centering and rotation operations, and
    ! the user may not want any symmetry - let critic2 guess.
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_abinit

  end subroutine read_abinit

  ! The following code has been adapted from the elk distribution, version 1.3.2
  ! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU General Public License.
  module subroutine read_elk(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, equal, getword,&
       zatguess, nameguess, fclose, string
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: isformat_elk
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< input filename
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, atname
    integer :: lu, i, zat, j, lp, idx
    integer :: natoms
    logical :: ok

    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! ignore the 'scale' stuff
    do i = 1, 14
       read(lu,*,err=999,end=999)
    end do

    read(lu,'(3G18.10)',err=999,end=999) seed%m_x2c(:,1)
    read(lu,'(3G18.10)',err=999,end=999) seed%m_x2c(:,2)
    read(lu,'(3G18.10)',err=999,end=999) seed%m_x2c(:,3)
    seed%useabr = 2

    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    ok = getline_raw(lu,line,.true.)
    if (.not.ok) goto 999
    if (equal(line,'molecule')) then
       errmsg = "Isolated molecules not supported."
       goto 999
    end if

    seed%nat = 0
    allocate(seed%x(3,10),seed%is(10))
    read(lu,'(I4)',err=999,end=999) seed%nspc
    allocate(seed%spc(seed%nspc))
    do i = 1, seed%nspc
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
       lp = 1
       atname = getword(line,lp)
       do j = 1, len(atname)
          if (atname(j:j) == "'") atname(j:j) = " "
          if (atname(j:j) == '"') atname(j:j) = " "
       end do
       zat = zatguess(atname)
       if (zat == -1) then
          errmsg = "Species file name must start with an atomic symbol"
          goto 999
       end if
       seed%spc(i)%z = zat

       idx = index(atname,".in",.true.)
       if (idx > 1) then
          seed%spc(i)%name = trim(atname(1:idx-1))
       else
          seed%spc(i)%name = trim(atname)
       end if

       read(lu,*,err=999,end=999) natoms
       do j = 1, natoms
          seed%nat = seed%nat + 1
          if (seed%nat > size(seed%x,2)) then
             call realloc(seed%x,3,2*seed%nat)
             call realloc(seed%is,2*seed%nat)
          end if
          read(lu,*,err=999,end=999) seed%x(:,seed%nat)
          seed%is(seed%nat) = i
       end do
    end do
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_elk

  end subroutine read_elk

  !> Read the structure from a molecule file
  module subroutine read_mol(seed,file,fmt,rborder,docube,errmsg,ti)
    use wfn_private, only: wfn_read_xyz_geometry, wfn_read_wfn_geometry, &
       wfn_read_wfx_geometry, wfn_read_fchk_geometry, wfn_read_molden_geometry,&
       wfn_read_log_geometry, wfn_read_dat_geometry, wfn_read_pgout_geometry,&
       wfn_read_orca_geometry, wfn_read_gjf_geometry
    use param, only: isformat_xyz, isformat_wfn, isformat_wfx,&
       isformat_fchk, isformat_molden, isformat_gaussian, isformat_dat,&
       isformat_pgout, isformat_orca, isformat_gjf, isformat_zmat, isformat_pdb
    use tools_io, only: equali
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: fmt !< wfn/wfx/xyz
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer, allocatable :: z(:)
    character*(10), allocatable :: name(:) !< Atomic names
    integer :: i, j, it

    call seed%end()
    errmsg = ""
    if (fmt == isformat_xyz) then
       ! xyz
       call wfn_read_xyz_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_gjf) then
       ! xyz
       call wfn_read_gjf_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_wfn) then
       ! wfn
       call wfn_read_wfn_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_wfx) then
       ! wfx
       call wfn_read_wfx_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_fchk) then
       ! fchk
       call wfn_read_fchk_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_molden) then
       ! molden (psi4)
       call wfn_read_molden_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_gaussian) then
       ! Gaussian output file
       call wfn_read_log_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_dat) then
       ! psi4 output file
       call wfn_read_dat_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_pgout) then
       ! postg output file
       call wfn_read_pgout_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_orca) then
       ! orca output file
       call wfn_read_orca_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_pdb) then
       call read_pdb_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    elseif (fmt == isformat_zmat) then
       call read_zmat_geometry(file,seed%nat,seed%x,z,name,errmsg,ti=ti)
    end if
    seed%useabr = 0
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.
    if (len_trim(errmsg) > 0) goto 999

    seed%nspc = 0
    allocate(seed%is(seed%nat),seed%spc(2))
    do i = 1, seed%nat
       it = 0
       do j = 1, seed%nspc
          if (equali(seed%spc(j)%name,name(i))) then
             it = j
             exit
          end if
       end do
       if (it == 0) then
          seed%nspc = seed%nspc + 1
          if (seed%nspc > size(seed%spc,1)) &
             call realloc(seed%spc,2*seed%nspc)
          seed%spc(seed%nspc)%name = name(i)
          seed%spc(seed%nspc)%z = z(i)
          it = seed%nspc
       end if
       seed%is(i) = it
    end do
    if (seed%nspc == 0) then
       errmsg = "No atomic species found."
       goto 999
    end if
    if (seed%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    call realloc(seed%spc,seed%nspc)

    errmsg = ""
999 continue

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = .true.
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = fmt

  end subroutine read_mol

  !> Read the structure from a quantum espresso output (file) and
  !> return it as a crystal seed. If mol, the structure is assumed to
  !> be a molecule (currently, the only effect is that its value is
  !> passed to the %ismolecule field). If istruct is zero, read the
  !> last geometry; otherwise, read geometry number istruct. If an
  !> error condition is found, return the error message in errmsg
  !> (zero-length string if no error).
  module subroutine read_qeout(seed,file,mol,istruct,errmsg,ti)
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    integer, intent(in) :: istruct !< structure number
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed), allocatable :: seedaux(:)
    integer :: nseed

    call seed%end()
    call read_all_qeout(nseed,seedaux,file,mol,istruct,errmsg,ti=ti)
    if (allocated(seedaux) .and. nseed >= 1) then
       seed = seedaux(1)
    end if

  end subroutine read_qeout

  !> Read the structure from a quantum espresso input
  module subroutine read_qein(seed,file,mol,errmsg,ti)
    ! This subroutine has been adapted from parts of the Quantum
    ! ESPRESSO code, version 4.3.2 and later.
    ! Copyright (C) 2002-2009 Quantum ESPRESSO group
    ! This file is distributed under the terms of the
    ! GNU General Public License. See the file `License'
    ! in the root directory of the present distribution,
    ! or http://www.gnu.org/copyleft/gpl.txt .
    use qe_private, only: qe_latgen, sup_spacegroup
    use tools_io, only: fopen_read, getline_raw, lower, getword,&
       equal, zatguess, fclose
    use tools_math, only: matinv
    use param, only: bohrtoa, isformat_qein
    use types, only: realloc
    integer, parameter :: dp = selected_real_kind(14,200)
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: nattot
    real(dp), allocatable :: tautot(:,:)
    integer, allocatable :: ityptot(:)

    integer, parameter :: ntypx = 10
    integer, parameter :: nsx = ntypx
    integer, parameter :: nspinx = 2
    integer, parameter :: lmaxx = 3
    integer, parameter :: lqmax = 2*lmaxx+1

    !!! Up to date with quantum espresso 6.3. More recent versions may
    !!! need additional keywords.

    ! from QE
    ! namelist control
    character(len=80) :: title, calculation, verbosity, restart_mode,&
       disk_io, memory
    character(len=10) :: point_label_type
    character(len=256) :: input_xml_schema_file
    integer :: nstep, iprint, isave, ndr, ndw, gdir, nppstr, nberrycyc, &
       printwfc
    logical :: tstress, tprnfor, tefield, tefield2, lelfield, dipfield, &
       lberry, wf_collect, saverho, tabps, lkpoint_dir, use_wannier, &
       lecrpa, tqmmm, lorbm, lfcpopt, lfcpdyn, gate
    real*8 :: dt, refg, max_seconds, ekin_conv_thr, etot_conv_thr, &
       forc_conv_thr
    character(len=256) :: outdir, prefix, pseudo_dir, wfcdir, vdw_table_name
    namelist /control/ title, calculation, verbosity, restart_mode,  &
       nstep, iprint, isave, tstress, tprnfor, dt, ndr, ndw, outdir,   &
       prefix, wfcdir, max_seconds, ekin_conv_thr, etot_conv_thr,      &
       forc_conv_thr, pseudo_dir, disk_io, tefield, dipfield, lberry,  &
       gdir, nppstr, wf_collect, printwfc, lelfield, nberrycyc, refg,  &
       tefield2, saverho, tabps, lkpoint_dir, use_wannier, lecrpa,     &
       vdw_table_name, tqmmm, lorbm, memory, point_label_type,         &
       lfcpopt, lfcpdyn, input_xml_schema_file, gate

    ! namelist system
    integer :: ibrav
    real*8 :: celldm(6)
    real*8 :: a, b, c, cosab, cosac, cosbc
    integer :: nat
    integer :: ntyp
    integer :: origin_choice
    integer :: space_group
    logical :: rhombohedral
    logical :: uniqueb
    real*8 :: tot_charge, tot_magnetization, ecutwfc, ecutrho, degauss, &
       ecfixed, qcutz, q2sigma, starting_magnetization(nsx), &
       starting_ns_eigenvalue(lqmax,nspinx,nsx), hubbard_u(nsx), &
       hubbard_alpha(nsx), a_pen(10,nspinx), sigma_pen(10), alpha_pen(10), &
       emaxpos, eopreg, eamp, lambda, fixed_magnetization(3), angle1(nsx), &
       angle2(nsx), b_field(3), sic_epsilon, sic_alpha, london_s6, london_rcut, &
       xdm_a1, xdm_a2, ts_sr, esm_efield, esm_w, &
       block_1, block_2, block_height, ecutfock, ecutvcut, esm_a, esm_zb, exx_fraction, &
       fcp_mass, fcp_mdiis_step, fcp_mu, fcp_relax_crit, fcp_relax_step, fcp_tempw, &
       hubbard_beta(nsx), hubbard_j0(nsx), hubbard_j(3,nsx), localization_thr, london_c6(nsx), &
       london_rvdw(nsx), ref_alat, scdmden, scdmgrd, screening_parameter, starting_charge(nsx), &
       ts_vdw_econv_thr, yukawa, zgate
    integer :: nbnd, nr1, nr2, nr3, nr1s, nr2s, nr3s, nr1b, nr2b, nr3b, &
       nspin, edir, report, xdm_usehigh, esm_nfit, esm_debug_gpmax, &
       dftd3_version, fcp_mdiis_size, lda_plus_u_kind, n_proj, nqx1, &
       nqx2, nqx3
    character(len=80) :: occupations, smearing, input_dft, u_projection_type, &
       constrained_magnetization, sic, assume_isolated
    logical :: nosym, noinv, nosym_evc, force_symmorphic, lda_plus_u, la2f, &
       step_pen, noncolin, lspinorb, starting_spin_angle, no_t_rev, force_pairing, &
       spline_ps, one_atom_occupations, london, xdm, xdm_onlyc, xdm_fixc6, &
       xdm_usec9, ts, ts_onlyc, esm_debug, &
       ace, block, dftd3_threebody, lforcet, relaxz, scdm, ts_vdw, ts_vdw_isolated, &
       use_all_frac, x_gamma_extrapolation
    character(len=3) :: esm_bc
    character(len=80) :: exxdiv_treatment, vdw_corr
    character(len=8) :: fcp_relax

    namelist /system/ ibrav, celldm, a, b, c, cosab, cosac, cosbc, nat,     &
       ntyp, nbnd, ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s, nr3s, nr1b, &
       nr2b, nr3b, nosym, nosym_evc, noinv, force_symmorphic, starting_magnetization, &
       occupations, degauss, nspin, ecfixed, qcutz, q2sigma, lda_plus_u, &
       hubbard_u, hubbard_alpha, edir, emaxpos, eopreg, eamp, smearing, &
       starting_ns_eigenvalue, u_projection_type, input_dft, la2f, assume_isolated, &
       noncolin, lspinorb, starting_spin_angle, lambda, angle1, angle2, report, &
       constrained_magnetization, b_field, fixed_magnetization, sic, sic_epsilon, &
       force_pairing, sic_alpha, tot_charge, tot_magnetization, spline_ps, &
       one_atom_occupations, london, london_s6, london_rcut, xdm, xdm_onlyc, &
       xdm_fixc6, xdm_usec9, xdm_usehigh, xdm_a1, xdm_a2, ts, ts_onlyc, ts_sr, &
       step_pen, a_pen, sigma_pen, alpha_pen, no_t_rev, esm_bc, esm_efield, &
       esm_w, esm_nfit, esm_debug, esm_debug_gpmax, use_all_frac, starting_charge, &
       lda_plus_u_kind, hubbard_j, hubbard_j0, hubbard_beta, nqx1, nqx2, nqx3, &
       ecutfock, localization_thr, scdm, ace, scdmden, scdmgrd, n_proj, exxdiv_treatment, &
       x_gamma_extrapolation, yukawa, ecutvcut, exx_fraction, screening_parameter, &
       ref_alat, lforcet, vdw_corr, london_c6, london_rvdw, dftd3_version, dftd3_threebody, &
       ts_vdw, ts_vdw_isolated, ts_vdw_econv_thr, esm_a, esm_zb, fcp_mu, fcp_mass, &
       fcp_tempw, fcp_relax, fcp_relax_step, fcp_relax_crit, fcp_mdiis_size, &
       fcp_mdiis_step, space_group, uniqueb, origin_choice, rhombohedral, &
       zgate, relaxz, block, block_1, block_2, block_height

    ! namelist electrons
    real*8 :: emass, emass_cutoff, ortho_eps, electron_damping, ekincw, fnosee, &
       ampre, grease, diis_hcut, diis_wthr, diis_delt, diis_fthr, diis_temp, &
       diis_achmix, diis_g0chmix, diis_g1chmix, diis_rothr, diis_ethr, mixing_beta,&
       diago_thr_init, conv_thr, lambda_cold, fermi_energy, rotmass, occmass,&
       occupation_damping, rotation_damping, etresh, passop, efield, efield_cart(3),&
       efield2, emass_emin, emass_cutoff_emin, electron_damping_emin, dt_emin
    character(len=80) :: orthogonalization, electron_dynamics, electron_velocities,&
       electron_temperature, startingwfc, mixing_mode, diagonalization, startingpot,&
       rotation_dynamics, occupation_dynamics, efield_phase
    integer :: ortho_max, electron_maxstep, diis_size, diis_nreset, diis_maxstep, &
       diis_nchmix, diis_nrot(3), mixing_ndim, diago_cg_maxiter, diago_david_ndim, &
       mixing_fixed_ns, n_inner, niter_cold_restart, maxiter, niter_cg_restart, &
       epol, epol2
    logical :: diis_rot, diis_chguess, diago_full_acc, tcg, real_space, tqr,&
       occupation_constraints, &
       scf_must_converge, tq_smoothing, tbeta_smoothing, adaptive_thr, tcpbo
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
       niter_cold_restart, lambda_cold, efield_cart, real_space,    &
       scf_must_converge, tq_smoothing, tbeta_smoothing, adaptive_thr, &
       tcpbo,emass_emin, emass_cutoff_emin, electron_damping_emin,  &
       dt_emin, efield_phase

    ! namelist ions
    character(len=80) :: phase_space, ion_dynamics, ion_positions, ion_velocities,&
       ion_temperature, pot_extrapolation, wfc_extrapolation
    integer, parameter :: nhclm   = 4
    integer, parameter :: max_nconstr = 100
    integer :: nhpcl, nhptyp, nhgrp(nsx), ndega, ion_nstepe, ion_maxstep, nraise,&
       bfgs_ndim, fe_nstep, sw_nstep, eq_nstep, n_muller, np_muller
    real*8 :: ion_radius(nsx), ion_damping, tempw, fnosep(nhclm), tolp, fnhscl(nsx),&
       amprp(nsx), greasp, upscale, delta_t, trust_radius_max, trust_radius_min,&
       trust_radius_ini, w_1, w_2, sic_rloc, g_amplitude, fe_step(max_nconstr)
    logical :: tranp(nsx), refold_pos, remove_rigid_rot, l_mplathe, l_exit_muller

    namelist /ions/ phase_space, ion_dynamics, ion_radius, ion_damping,  &
       ion_positions, ion_velocities, ion_temperature,      &
       tempw, fnosep, nhgrp, fnhscl, nhpcl, nhptyp, ndega, tranp,   &
       amprp, greasp, tolp, ion_nstepe, ion_maxstep,        &
       refold_pos, upscale, delta_t, pot_extrapolation,     &
       wfc_extrapolation, nraise, remove_rigid_rot,         &
       trust_radius_max, trust_radius_min,                  &
       trust_radius_ini, w_1, w_2, bfgs_ndim, sic_rloc,     &
       fe_step, fe_nstep, sw_nstep, eq_nstep, g_amplitude, &
       l_mplathe, n_muller, np_muller, l_exit_muller

    ! namelist cell
    character(len=80) :: cell_parameters, cell_dynamics, cell_velocities, &
       cell_temperature, cell_dofree
    real(dp) :: press, wmass, temph, fnoseh, greash, cell_factor, cell_damping,&
       press_conv_thr
    integer :: cell_nstepe
    namelist /cell/ cell_parameters, cell_dynamics, cell_velocities, &
       press, wmass, cell_temperature, temph, fnoseh,   &
       cell_dofree, greash, cell_factor, cell_nstepe,   &
       cell_damping, press_conv_thr

    ! local to this routine
    integer :: lu, ios, lp, i, j, ier
    character(len=:), allocatable :: line, word
    character*10 :: atm
    real*8 :: r(3,3)
    integer :: iunit, cunit, ibrav_sg
    integer, parameter :: icrystal = 1
    integer, parameter :: ibohr = 2
    integer, parameter :: iang = 3
    integer, parameter :: ialat = 4
    integer, parameter :: icrystalsg = 5
    integer, allocatable :: rd_if_pos(:,:)
    real*8, allocatable :: rd_for(:,:)

    ! initialize
    call seed%end()
    celldm = 0d0
    ibrav = -1
    nat = 0
    ntyp = 0
    origin_choice = 1
    space_group = 0
    rhombohedral = .false.
    uniqueb = .false.
    cell_nstepe = 1

    ! open
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    r = 0d0
    calculation = ""

    ! read the namelists
    read(lu,control,iostat=ios)
    if (ios /= 0) then
       errmsg = "Wrong namelist control."
       goto 999
    end if
    read(lu,system,iostat=ios)
    if (ios/=0) then
       errmsg = "Wrong namelist system."
       goto 999
    end if
    read(lu,electrons,iostat=ios)
    if (ios/=0) then
       errmsg = "Wrong namelist electrons."
       goto 999
    end if
    if (trim(calculation)=='relax'.or.trim(calculation)=='md'.or.&
       trim(calculation)=='vc-relax'.or.trim(calculation)=='vc-md'.or.&
       trim(calculation)=='cp'.or.trim(calculation)=='vc-cp'.or.&
       trim(calculation)=='smd'.or.trim(calculation)=='cp-wf') then
       read(lu,ions,iostat=ios)
       if (ios/=0) then
          errmsg = "Wrong namelist ions."
          goto 999
       end if
    endif
    if (trim(calculation)=='vc-relax'.or.trim(calculation)=='vc-md'.or.&
       trim(calculation)=='vc-cp') then
       read(lu,cell,iostat=ios)
       if (ios/=0) then
          errmsg = "Wrong namelist ions."
          goto 999
       end if
    end if

    ! allocate space for atoms
    seed%nat = nat
    seed%nspc = ntyp
    allocate(seed%x(3,nat),seed%is(nat),seed%spc(ntyp))

    ! read the cards
    iunit = icrystal
    cunit = ialat
    do while (getline_raw(lu,line))
       line = lower(line)
       lp = 1
       word = getword(line,lp)
       if (equal(word,'atomic_species')) then
          do i = 1, ntyp
             read (lu,*,iostat=ios) seed%spc(i)%name
             if (ios/=0) then
                errmsg = "Error reading atomic species."
                goto 999
             end if
             seed%spc(i)%z = zatguess(seed%spc(i)%name)
          end do

       else if (equal(word,'atomic_positions')) then
          word = getword(line,lp)
          if (index(word,"crystal_sg") > 0) then
             iunit = icrystalsg
          elseif (index(word,"crystal") > 0) then
             iunit = icrystal
          elseif (index(word,"bohr") > 0) then
             iunit = ibohr
          elseif (index(word,"angstrom") > 0) then
             iunit = iang
          elseif (index(word,"alat") > 0) then
             iunit = ialat
          else
             iunit = ialat
          end if
          do i = 1, nat
             read (lu,*,iostat=ios) atm, seed%x(:,i)
             if (ios/=0) then
                errmsg = "Error reading atomic positions."
                goto 999
             end if
             seed%is(i) = 0
             do j = 1, seed%nspc
                if (equal(seed%spc(j)%name,atm)) then
                   seed%is(i) = j
                   exit
                end if
             end do
             if (seed%is(i) == 0) then
                errmsg = "Could not find atomic species "//trim(atm)//"."
                goto 999
             end if
          end do
       elseif (equal(word,'cell_parameters')) then
          word = getword(line,lp)
          cunit = ialat
          if (index(word,"bohr") > 0) then
             cunit = ibohr
          elseif (index(word,"angstrom") > 0) then
             cunit = iang
          elseif (index(word,"alat") > 0) then
             cunit = ialat
          elseif (len_trim(word) == 0) then
             cunit = ialat
          else
             cunit = ibohr
          end if
          do i = 1, 3
             read (lu,*,iostat=ios) (r(i,j),j=1,3)
             if (ios/=0) then
                errmsg = "Error reading cell parameters."
                goto 999
             end if
          end do
       endif
    end do

    ! figure it out
    if (space_group /= 0) then
       ! space group number was given
       if (iunit /= icrystalsg) then
          errmsg = "space_group requires crystal_sg atomic coordinates"
          goto 999
       end if

       allocate(rd_if_pos(3,nat),rd_for(3,nat))
       rd_if_pos = 1
       rd_for = 0d0
       CALL sup_spacegroup(seed%x,seed%is,rd_for,rd_if_pos,space_group,&
          nat,uniqueb,rhombohedral,origin_choice,ibrav_sg,nattot,tautot,&
          ityptot)
       deallocate(rd_if_pos, rd_for)

       if (ibrav/=-1 .and. ibrav /= ibrav_sg) then
          errmsg = "input ibrav not compatible with space_group number"
          goto 999
       end if
       ibrav = ibrav_sg

       ! reallocate using the new info
       seed%nat = nattot
       if (allocated(seed%x)) deallocate(seed%x)
       if (allocated(seed%is)) deallocate(seed%is)
       allocate(seed%x(3,nattot),seed%is(nattot))
       do i = 1, nattot
          seed%x(:,i) = tautot(:,i)
          seed%is(i) = ityptot(i)
       end do

       ! calculate the new r
       call qe_latgen(ibrav,celldm,r(:,1),r(:,2),r(:,3),errmsg)
       if (len_trim(errmsg) > 0) goto 999

    else if (ibrav == 0) then
       ! ibrav = 0 and CELL_PARAMETERS
       if (cunit == ialat) then
          if (celldm(1) /= 0.D0) r = r * celldm(1)
       elseif (cunit == iang) then
          r = r / bohrtoa
       end if
       r = transpose(r)
    else if (ibrav == -1) then
       ! neither ibrav nor space_group, this is an error
       errmsg = "no ibrav or space_group found"
       goto 999
    else
       ! a non-zero ibrav
       call qe_latgen(ibrav,celldm,r(:,1),r(:,2),r(:,3),errmsg)
       if (len_trim(errmsg) > 0) goto 999
    endif

    ! fill the cell metrics
    seed%m_x2c = r
    call matinv(r,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if
    seed%useabr = 2

    ! do the atom stuff
    do i = 1, nat
       if (iunit == ialat) then
          seed%x(:,i) = matmul(r,seed%x(:,i) * celldm(1))
       elseif (iunit == ibohr) then
          seed%x(:,i) = matmul(r,seed%x(:,i))
       elseif (iunit == iang) then
          seed%x(:,i) = matmul(r,seed%x(:,i) / bohrtoa)
       endif
       seed%x(:,i) = seed%x(:,i) - floor(seed%x(:,i))
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_qein

  end subroutine read_qein

  !> Read the structure from a crystal output
  module subroutine read_crystalout(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, isinteger, isreal,&
       zatguess, fclose, equali
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: bohrtoa, isformat_crystal
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, ier
    character(len=:), allocatable :: line
    integer :: idum, iz, lp
    real*8 :: r(3,3), x(3)
    logical :: ok, iscrystal
    character*(10) :: ats

    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    errmsg = "Error reading file: " // trim(file)
    r = 0d0
    iscrystal = .false.
    allocate(seed%x(3,10),seed%is(10),seed%spc(2))
    seed%nat = 0
    seed%nspc = 0
    ! rewind and read the correct structure
    rewind(lu)
    do while (getline_raw(lu,line))
       if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) then
          cycle
       elseif (index(line,"CRYSTAL CALCULATION") > 0) then
          iscrystal = .true.
       elseif (index(line,"3D - CRYSTAL") > 0) then
          iscrystal = .true.
       elseif (index(line,"MOLECULAR CALCULATION") > 0) then
          iscrystal = .false.
       elseif (index(line,"DIRECT LATTICE VECTORS CARTESIAN COMPONENTS") > 0) then
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
          if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
          do i = 1, 3
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
             lp = 1
             ok = isreal(r(i,1),line,lp)
             ok = ok.and.isreal(r(i,2),line,lp)
             ok = ok.and.isreal(r(i,3),line,lp)
             if (.not.ok) then
                errmsg = "Wrong lattice vectors."
                goto 999
             end if
          end do
          r = r / bohrtoa
       elseif (index(line,"CARTESIAN COORDINATES - PRIMITIVE CELL") > 0) then
          do i = 1, 3
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
          end do
          line = ""
          seed%nat = 0
          do while (.true.)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
             if (len_trim(line) < 1) exit
             if (index(line,"ERROR") > 0) exit
             if (index(line,"WARNING") > 0) exit
             seed%nat = seed%nat + 1
             if (seed%nat > size(seed%x,2)) then
                call realloc(seed%x,3,2*seed%nat)
                call realloc(seed%is,2*seed%nat)
             end if
             read (line,*,err=999,end=999) idum, iz, ats, x
             seed%x(:,seed%nat) = x / bohrtoa
             seed%is(seed%nat) = 0
             do j = 1, seed%nspc
                if (equali(trim(ats),seed%spc(j)%name)) then
                   seed%is(seed%nat) = j
                   exit
                end if
             end do
             if (seed%is(seed%nat) == 0) then
                seed%nspc = seed%nspc + 1
                if (seed%nspc > size(seed%spc,1)) &
                   call realloc(seed%spc,2*seed%nspc)
                seed%spc(seed%nspc)%name = trim(ats)
                seed%spc(seed%nspc)%z = zatguess(ats)
                seed%is(seed%nat) = seed%nspc
             end if
          end do
       end if
    end do
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%spc,seed%nspc)

    if (.not.iscrystal) then
       errmsg = "Only CRYSTAL calculations supported (no MOLECULE, SLAB or POLYMER)."
       goto 999
    end if
    if (all(r == 0d0)) then
       errmsg = "Could not find lattice vectors."
       goto 999
    end if

    ! cell
    seed%m_x2c = transpose(r)
    r = seed%m_x2c
    call matinv(r,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if
    seed%useabr = 2

    ! atoms
    do i = 1, seed%nat
       seed%x(:,i) = matmul(r,seed%x(:,i))
       seed%x(:,i) = seed%x(:,i) - floor(seed%x(:,i))
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_crystal

  end subroutine read_crystalout

  !> Read the structure from an FPLO output
  module subroutine read_fploout(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, isinteger, isreal,&
       zatguess, nameguess, fclose
    use tools_math, only: matinv
    use param, only: maxzat, isformat_fploout
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, ll, ier
    character(len=:), allocatable :: line
    integer :: iz, lp, idum1, idum2, idum3
    logical :: ok, iscell, isatoms
    character*(10) :: ats
    integer :: usez(maxzat)
    real*8 :: r(3,3)

    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    errmsg = "Error reading file: " // trim(file)
    iscell = .false.
    isatoms = .false.
    seed%nat = 0
    seed%nspc = 0
    ! rewind and read the correct structure
    do while (getline_raw(lu,line))
       ll = len(line)
       if (ll == 15) then
          if (line == "lattice vectors") then
             do i = 1, 3
                lp = 1
                ok = getline_raw(lu,line)
                line = line(12:)
                ok = ok .and. isreal(r(i,1),line,lp)
                ok = ok .and. isreal(r(i,2),line,lp)
                ok = ok .and. isreal(r(i,3),line,lp)
             end do
             if (.not.ok) then
                errmsg = "Error reading lattice vectors."
                goto 999
             end if
             iscell = .true.
          end if
       elseif (ll == 27) then
          if (line(18:27) == "Atom sites") then
             ! number of sites
             ok = getline_raw(lu,line)
             ok = ok .and. getline_raw(lu,line)
             ok = ok .and. isinteger(seed%nat,line(18:))
             do i = 1, 3
                ok = ok .and. getline_raw(lu,line)
             end do
             if (.not.ok) then
                errmsg = "Error reading number of sites."
                goto 999
             end if

             ! sites
             seed%nspc = 0
             if (allocated(seed%x)) deallocate(seed%x)
             if (allocated(seed%is)) deallocate(seed%is)
             allocate(seed%x(3,seed%nat),seed%is(seed%nat))
             usez = 0
             do i = 1, seed%nat
                ok = getline_raw(lu,line)
                read (line,*) idum1, ats, idum2, idum3, seed%x(:,i)
                iz = zatguess(ats)
                if (iz < 1 .or. iz > maxzat) then
                   errmsg = "Unknown atomic species: " // ats // "."
                   goto 999
                end if

                if (usez(iz) == 0) then
                   seed%nspc = seed%nspc + 1
                   usez(iz) = seed%nspc
                   seed%is(i) = seed%nspc
                else
                   seed%is(i) = usez(iz)
                end if
             end do

             ! species
             if (allocated(seed%spc)) deallocate(seed%spc)
             allocate(seed%spc(seed%nspc))
             do iz = 1, maxzat
                if (usez(iz) > 0) then
                   i = usez(iz)
                   seed%spc(i)%name = nameguess(iz,.true.)
                   seed%spc(i)%z = iz
                   seed%spc(i)%qat = 0d0
                end if
             end do

             ! all done
             isatoms = (seed%nat > 0) .and. (seed%nspc > 0)
          end if
       end if
    end do
    if (.not.iscell) then
       errmsg = "No lattice parameters found."
       goto 999
    end if
    if (.not.isatoms) then
       errmsg = "No atoms found."
       goto 999
    end if

    ! cell
    seed%m_x2c = transpose(r)
    r = seed%m_x2c
    call matinv(r,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting lattice parameter matrix"
       goto 999
    end if
    seed%useabr = 2

    ! atoms
    do i = 1, seed%nat
       seed%x(:,i) = matmul(r,seed%x(:,i))
       seed%x(:,i) = seed%x(:,i) - floor(seed%x(:,i))
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_fploout

  end subroutine read_fploout

  !> Read the structure from a siesta OUT input
  module subroutine read_siesta(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, nameguess, fclose
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: bohrtoa, isformat_siesta
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, idum

    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! the lattice vectors
    do i = 1, 3
       read (lu,*,err=999,end=999) r(i,:)
    end do
    r = r / bohrtoa

    ! the atoms
    seed%nspc = 0
    read (lu,*,err=999,end=999) seed%nat
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%spc(2))
    do i = 1, seed%nat
       read (lu,*,err=999,end=999) seed%is(i), idum, seed%x(:,i)
       if (idum > size(seed%spc,1)) &
          call realloc(seed%spc,2*idum)
       seed%nspc = max(seed%nspc,seed%is(i))
       seed%spc(seed%is(i))%z = idum
       seed%spc(seed%is(i))%name = nameguess(idum)
    end do
    call realloc(seed%spc,seed%nspc)

    ! fill the cell metrics
    seed%m_x2c = transpose(r)
    seed%useabr = 2

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_siesta

  end subroutine read_siesta

  !> Read the structure from a CASTEP cell file
  module subroutine read_castep_cell(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, lgetword, getline_raw,&
       equal, getword, isinteger, zatguess, isreal, lower
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: bohrtoa, bohrtom, bohrtocm, bohrtonm, isformat_castepcell
    use hashmod, only: hash
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, lword
    integer :: lu, lp, idum, ier, i
    logical :: ok, iscart
    real*8 :: rconv, m(3,3)
    type(hash) :: usen

    call usen%init()
    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    iscart = .false.
    seed%useabr = 0
    do while(get_next_line())
       lp = 1
       word = lgetword(line,lp)
       if (index(word,"%block") == 1) then
          word = lgetword(line,lp)
          if (equal(word,"lattice_abc")) then
             !! cell parameters
             ! optional unit
             if (.not.get_optional_unit(rconv)) goto 999

             ! read the rest
             read(line,*,err=999,end=999) seed%aa
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) seed%bb

             ! conversion
             seed%useabr = 1
             seed%aa = seed%aa * rconv
          elseif (equal(word,"lattice_cart")) then
             !! lattice vectors
             ! optional unit
             if (.not.get_optional_unit(rconv)) goto 999

             ! read the rest
             read(line,*,err=999,end=999) seed%m_x2c(:,1)
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) seed%m_x2c(:,2)
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) seed%m_x2c(:,3)
             if (.not.get_next_line()) goto 999

             ! conversion
             seed%useabr = 2
             seed%m_x2c = seed%m_x2c * rconv

          elseif (equal(word,"positions_frac").or.equal(word,"positions_abs")) then
             iscart = equal(word,"positions_abs")
             if (iscart) then
                if (.not.get_optional_unit(rconv)) goto 999
             else
                if (.not.get_next_line()) goto 999
             end if

             !! positions in fractional/Cartesian coordinates
             ok = .true.
             seed%nat = 0
             seed%nspc = 0
             allocate(seed%x(3,10),seed%is(10),seed%spc(2))
             do while(ok)
                if (line(1:1) == "%") exit
                seed%nat = seed%nat + 1
                if (seed%nat > size(seed%x,2)) then
                   call realloc(seed%x,3,2*seed%nat)
                   call realloc(seed%is,2*seed%nat)
                end if

                ! read the atomic symbol or number
                lp = 1
                word = getword(line,lp)
                lword = lower(word)
                if (usen%iskey(lword)) then
                   seed%is(seed%nat) = usen%get(lword,1)
                else
                   seed%nspc = seed%nspc + 1
                   if (seed%nspc > size(seed%spc,1)) &
                      call realloc(seed%spc,2*seed%nspc)
                   seed%spc(seed%nspc)%name = trim(word)
                   if (isinteger(idum,word)) then
                      seed%spc(seed%nspc)%z = idum
                   else
                      seed%spc(seed%nspc)%z = zatguess(word)
                      if (seed%spc(seed%nspc)%z < 0) then
                         errmsg = "Unknown atomic symbol: " // word
                         goto 999
                      end if
                   end if
                   seed%spc(seed%nspc)%qat = 0d0
                   call usen%put(lword,seed%nspc)
                   seed%is(seed%nat) = seed%nspc
                end if

                ! read the atomic coordinates
                ok = isreal(seed%x(1,seed%nat),line,lp)
                ok = ok .and. isreal(seed%x(2,seed%nat),line,lp)
                ok = ok .and. isreal(seed%x(3,seed%nat),line,lp)
                if (.not.ok) goto 999

                ok = get_next_line()
             end do
             if (seed%nat == 0 .or. seed%nspc == 0) then
                errmsg = "Error reading atoms"
                goto 999
             end if
             call realloc(seed%x,3,seed%nat)
             call realloc(seed%is,seed%nat)
             call realloc(seed%spc,seed%nspc)

             ! conversion factor
             if (iscart) seed%x = seed%x * rconv
          end if
       end if
    end do
    ! consistency checks
    if (seed%useabr == 0) then
       errmsg = "No lattice block found"
       goto 999
    end if
    if (seed%nat == 0 .or. seed%nspc == 0) then
       errmsg = "No atoms found"
       goto 999
    end if
    if (iscart .and. seed%useabr /= 2) then
       errmsg = "Atomic Cartesian (absolute) coordinates require lattice_cart"
       goto 999
    end if

    ! transform to fractional coordinates
    if (iscart) then
       m = seed%m_x2c
       call matinv(m,3,ier)
       if (ier /= 0) then
          errmsg = "error inverting lattice vector matrix"
          goto 999
       end if
       do i = 1, seed%nat
          seed%x(:,i) = matmul(m,seed%x(:,i))
       end do
    end if

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_castepcell

  contains
    function get_optional_unit(rconv)
      real*8, intent(out) :: rconv
      logical :: get_optional_unit

      logical :: ok, readnew
      integer :: lp

      get_optional_unit = .false.
      rconv = 1d0 / bohrtoa
      ok = get_next_line()
      if (.not.ok) return
      lp = 1
      word = lgetword(line,lp)
      readnew = .true.
      if (equal(word,"ang")) then
         rconv = 1d0 / bohrtoa
      elseif (equal(word,"bohr") .or. equal(word,"a0")) then
         rconv = 1d0
      elseif (equal(word,"m")) then
         rconv = 1d0 / bohrtom
      elseif (equal(word,"cm")) then
         rconv = 1d0 / bohrtocm
      elseif (equal(word,"nm")) then
         rconv = 1d0 / bohrtonm
      else
         readnew = .false.
      end if
      if (readnew) then
         ok = getline_raw(lu,line)
         if (.not.ok) return
      end if
      get_optional_unit = .true.

    end function get_optional_unit

    function get_next_line()
      logical :: get_next_line

      do while(getline_raw(lu,line))
         line = trim(adjustl(line))
         if (skipline(line)) cycle
         get_next_line = .true.
         return
      end do
      get_next_line = .false.

    end function get_next_line

    function skipline(line)
      use tools_io, only: lower
      character*(*), intent(in) :: line
      logical :: skipline

      integer :: alen

      alen = len(line)
      skipline = .true.
      if (alen == 0) return
      if (line(1:1) == "#" .or. line(1:1) == "!" .or. line(1:1) == ";") return
      if (alen >= 7) then
         if (lower(line(1:7)) == "comment") return
      end if
      skipline = .false.
    end function skipline

  end subroutine read_castep_cell

  !> Read the structure from a CASTEP geom file
  module subroutine read_castep_geom(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword,&
       getword, lower, isinteger, isreal, zatguess
    use tools_math, only: matinv
    use hashmod, only: hash
    use param, only: isformat_castepgeom
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, lword
    integer :: lu, ll, i, lp, idum, is, ier
    integer :: nlast, n, nat, nspc
    logical :: ok
    type(hash) :: usen
    logical, allocatable :: usespc(:)
    real*8 :: m(3,3)

    call usen%init()
    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! first pass, read the number of lines and the number of atoms
    nlast = 0
    n = 0
    nat = 0
    nspc = 0
    do while(getline_raw(lu,line))
       n = n + 1
       line = trim(adjustl(line))
       ll = len(line)
       if (line(ll-4:ll) == "<-- E") then
          nlast = n
          nat = 0
       elseif (line(ll-4:ll) == "<-- R") then
          nat = nat + 1
          lp = 1
          word = lgetword(line,lp)
          if (.not.usen%iskey(word)) then
             nspc = nspc + 1
             call usen%put(word,nspc)
          end if
       end if
    end do
    if (nat == 0) goto 999
    if (nspc == 0) goto 999

    ! second pass, actual read
    seed%useabr = 2
    seed%nat = nat
    seed%nspc = nspc
    allocate(seed%x(3,nat),seed%is(nat),seed%spc(nspc),usespc(nspc))
    usespc = .false.
    rewind(lu)
    do i = 1, nlast-1
       read (lu,*,err=999,end=999)
    end do
    nat = 0
    do while(getline_raw(lu,line))
       line = trim(adjustl(line))
       ll = len(line)
       if (line(ll-4:ll) == "<-- E") then
          read(line,*,err=999,end=999) seed%energy
       elseif (line(ll-4:ll) == "<-- h") then
          read(line,*,err=999,end=999) seed%m_x2c(:,1)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) seed%m_x2c(:,2)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) seed%m_x2c(:,3)
       elseif (line(ll-4:ll) == "<-- R") then
          nat = nat + 1
          lp = 1
          word = getword(line,lp)
          lword = lower(word)
          ok = isinteger(idum,line,lp)
          ok = ok .and. isreal(seed%x(1,nat),line,lp)
          ok = ok .and. isreal(seed%x(2,nat),line,lp)
          ok = ok .and. isreal(seed%x(3,nat),line,lp)
          if (.not.ok) goto 999

          is = usen%get(lword,1)
          seed%is(nat) = is
          if (.not.usespc(is)) then
             seed%spc(is)%name = trim(word)
             seed%spc(is)%qat = 0d0
             if (isinteger(idum,word)) then
                seed%spc(is)%z = idum
             else
                seed%spc(is)%z = zatguess(word)
                if (seed%spc(is)%z < 0) then
                   errmsg = "Unknown atomic symbol: " // word
                   goto 999
                end if
             end if
             usespc(is) = .true.
          end if
       end if
    end do

    ! transform to fractional coordinates
    m = seed%m_x2c
    call matinv(m,3,ier)
    if (ier /= 0) then
       errmsg = "error inverting lattice vector matrix"
       goto 999
    end if
    do i = 1, seed%nat
       seed%x(:,i) = matmul(m,seed%x(:,i))
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_castepgeom

  end subroutine read_castep_geom

  !> Read the structure from a DMACRYS input file (dmain)
  module subroutine read_dmain(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, nameguess, fclose, getline_raw, lgetword, isreal,&
       getword, zatguess, isinteger, lower
    use tools_math, only: matinv
    use param, only: bohrtoa, maxzat, isformat_dmain
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer, allocatable :: usedz(:)
    character(len=:), allocatable :: line, word, lword
    character*2 :: atsym
    logical :: havelatt, haverlscal, ok
    integer :: lu, lp, iz, lvl, nrec, ier
    real*8 :: r(3,3), rlscal
    integer :: i

    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! loop over the lines
    haverlscal = .false.
    havelatt = .false.
    rlscal = 1d0
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       lp = 1
       word = lgetword(line,lp)

       ! read the lattice vectors
       if (word == "latt" .and..not.havelatt) then
          do i = 1, 3
             read (lu,*,err=999,end=999) r(i,:)
          end do
          havelatt = .true.
       elseif (word == "cuto") then
          ! lattice vector scaling factor
          ok = isreal(rlscal,line,lp)
          haverlscal = .true.
       elseif (word == "basi") then
          seed%nspc = 0
          seed%nat = 0
          allocate(seed%x(3,10),seed%is(10),seed%spc(2),usedz(maxzat))
          usedz = 0
          do while(.true.)
             ! get the atom
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             lp = 1
             word = getword(line,lp)
             lword = lower(word)
             if (lword == "ends") exit
             atsym = word(1:2)
             iz = zatguess(atsym)

             ! write down the species
             if (usedz(iz) == 0) then
                seed%nspc = seed%nspc + 1
                if (seed%nspc > size(seed%spc,1)) &
                   call realloc(seed%spc,2*seed%nspc)
                seed%spc(seed%nspc)%z = iz
                seed%spc(seed%nspc)%name = atsym
                usedz(iz) = seed%nspc
             end if

             ! write down the atom
             seed%nat = seed%nat + 1
             if (seed%nat > size(seed%x,2)) then
                call realloc(seed%x,3,2*seed%nat)
                call realloc(seed%is,2*seed%nat)
             end if
             ok = isreal(seed%x(1,seed%nat),line,lp)
             ok = ok .and. isreal(seed%x(2,seed%nat),line,lp)
             ok = ok .and. isreal(seed%x(3,seed%nat),line,lp)
             if (.not.ok) goto 999
             seed%is(seed%nat) = usedz(iz)

             ! process the block for this atom
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             lp = 1
             word = lgetword(line,lp)
             if (word == "level") then
                ok = isinteger(lvl,line,lp)
                nrec = 0
                do i = 0, lvl
                   nrec = nrec + (2*i+1 - 1)/7 + 1
                end do
             else
                nrec = 0
             end if

             ! skip the extra lines
             ok = .true.
             do i = 1, nrec
                ok = ok .and. getline_raw(lu,line)
             end do
             if (.not.ok) goto 999
          end do
          call realloc(seed%spc,seed%nspc)
          call realloc(seed%x,3,seed%nat)
          exit
       elseif (word == "stop") then
          exit
       end if
    end do

    ! lattice vectors
    seed%useabr = 2
    seed%m_x2c = transpose(r) * rlscal / bohrtoa

    ! transform atoms
    r = transpose(r)
    call matinv(r,3,ier)
    if (ier /= 0) goto 999
    do i = 1, seed%nat
       seed%x(:,i) = matmul(r,seed%x(:,i))
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_dmain

  end subroutine read_dmain

  !> Read the structure from a file in DFTB+ gen format.
  module subroutine read_dftbp(seed,file,rborder,docube,errmsg,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, getline, lower, equal, &
       getword, zatguess, nameguess, fclose
    use param, only: bohrtoa, isformat_gen
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ier
    real*8 :: r(3,3)
    integer :: i, iz, idum, lp
    logical :: ok, molout
    character*1 :: isfrac
    character(len=:), allocatable :: line, word

    ! open
    call seed%end()
    molout = .false.
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! number of atoms and type of coordinates
    ok = getline(lu,line)
    if (.not.ok) goto 999
    read (line,*,err=999,end=999) seed%nat, isfrac
    isfrac = lower(isfrac)
    if (.not.(equal(isfrac,"f").or.equal(isfrac,"c").or.equal(isfrac,"s"))) then
       errmsg = 'Wrong coordinate selector.'
       goto 999
    end if
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))

    ! atom types
    seed%nspc = 0
    allocate(seed%spc(2))
    ok = getline(lu,line)
    if (.not.ok) goto 999
    lp = 1
    word = getword(line,lp)
    iz = zatguess(word)
    do while (iz >= 0)
       seed%nspc = seed%nspc + 1
       if (seed%nspc > size(seed%spc,1)) &
          call realloc(seed%spc,2*seed%nspc)
       seed%spc(seed%nspc)%z = iz
       seed%spc(seed%nspc)%name = nameguess(iz)
       word = getword(line,lp)
       iz = zatguess(word)
    end do
    if (seed%nspc == 0) then
       errmsg = 'No atomic types found.'
       goto 999
    end if
    call realloc(seed%spc,seed%nspc)

    ! read atomic positions
    do i = 1, seed%nat
       ok = getline(lu,line)
       if (.not.ok) goto 999
       read (line,*,err=999,end=999) idum, seed%is(i), seed%x(:,i)
       if (isfrac /= "f") &
          seed%x(:,i) = seed%x(:,i) / bohrtoa
    end do

    ! read lattice vectors, if they exist
    ok = getline(lu,line)
    if (ok) then
       do i = 1, 3
          ok = getline(lu,line,.true.)
          read (line,*) r(i,:)
       end do
       r = r / bohrtoa

       ! fill the cell metrics
       seed%m_x2c = transpose(r)
       r = seed%m_x2c
       call matinv(r,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting matrix"
          goto 999
       end if

       if (isfrac == "c") then
          errmsg = 'Lattice plus C not supported.'
          goto 999
       elseif (isfrac == "s") then
          do i = 1, seed%nat
             seed%x(:,i) = matmul(r,seed%x(:,i))
          end do
       end if
       seed%useabr = 2
       molout = .false.
    else
       ! molecule and no lattice -> set up the origin and the molecular cell
       if (isfrac == "f" .or. isfrac == "s") then
          errmsg = 'S or F coordinates but no lattice vectors.'
          goto 999
       end if
       seed%useabr = 0
       molout = .true.
    end if

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = molout
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_gen

  end subroutine read_dftbp

  !> Read the structure from an xsf file.
  module subroutine read_xsf(seed,file,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, nameguess, equal,&
       zatguess, isinteger, getword, isreal, lower, string
    use tools_math, only: matinv
    use param, only: bohrtoa, isformat_xsf
    use types, only: realloc
    use hashmod, only: hash
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, name
    character*10 :: atn, latn
    integer :: lu, lp, i, j, iz, it, ier
    real*8 :: r(3,3), x(3)
    logical :: ok, ismol, xread
    type(hash) :: usen

    ! open
    call seed%end()
    ismol = .false.
    errmsg = ""
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    xread = .false.
    errmsg = "Error reading file: " // trim(file)
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"primvec")) then
          do i = 1, 3
             read (lu,*,err=999,end=999) r(i,:)
          end do
          r = r / bohrtoa
          ismol = .false.
       elseif (equal(word,"primcoord").and..not.xread) then
          read (lu,*,err=999,end=999) seed%nat
          allocate(seed%x(3,seed%nat),seed%is(seed%nat))
          seed%nspc = 0
          allocate(seed%spc(2))
          do i = 1, seed%nat
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             lp = 1
             ok = isinteger(iz,line,lp)
             if (ok) then
                ! Z x y z
                name = nameguess(iz,.true.)
             else
                word = getword(line,lp)
                name = trim(adjustl(word))
                iz = zatguess(name)
             end if
             ok = isreal(seed%x(1,i),line,lp)
             ok = ok.and.isreal(seed%x(2,i),line,lp)
             ok = ok.and.isreal(seed%x(3,i),line,lp)
             if (.not.ok) then
                errmsg = 'Wrong atomic position.'
                goto 999
             end if
             seed%x(:,i) = seed%x(:,i) / bohrtoa

             it = 0
             do j = 1, seed%nspc
                if (seed%spc(j)%z == iz) then
                   it = j
                   exit
                end if
             end do
             if (it == 0) then
                seed%nspc = seed%nspc + 1
                if (seed%nspc > size(seed%spc,1)) &
                   call realloc(seed%spc,2*seed%nspc)
                seed%spc(seed%nspc)%z = iz
                seed%spc(seed%nspc)%name = name
                it = seed%nspc
             end if
             seed%is(i) = it
          end do
          ismol = .false.
          xread = .true.
       elseif (equal(word,"atoms").and..not.xread) then
          ismol = .true.
          call usen%init()
          seed%nat = 0
          seed%nspc = 0
          allocate(seed%x(3,10),seed%spc(5),seed%is(10))
          do while (getline_raw(lu,line))
             if (len_trim(line) == 0) exit

             read (line,*,err=999,end=999) atn, x
             seed%nat = seed%nat + 1
             if (seed%nat > size(seed%x,2)) then
                call realloc(seed%x,3,2*seed%nat)
                call realloc(seed%is,2*seed%nat)
             end if
             seed%x(:,seed%nat) = x / bohrtoa

             ok = isinteger(iz,atn)
             if (.not.ok) then
                iz = zatguess(atn)
             end if
             if (iz < 0) then
                errmsg = "Unknown atomic symbol: "//trim(atn)//"."
                goto 999
             end if

             latn = lower(atn)
             if (usen%iskey(latn)) then
                seed%is(seed%nat) = usen%get(latn,1)
             else
                seed%nspc = seed%nspc + 1
                if (seed%nspc > size(seed%spc,1)) &
                   call realloc(seed%spc,2*seed%nspc)
                seed%spc(seed%nspc)%name = trim(atn)
                seed%spc(seed%nspc)%z = iz
                seed%spc(seed%nspc)%qat = 0d0
                call usen%put(latn,seed%nspc)
                seed%is(seed%nat) = seed%nspc
             end if
          end do
          call realloc(seed%x,3,seed%nat)
          call realloc(seed%is,seed%nat)
          call realloc(seed%spc,seed%nspc)
          xread = .true.
       end if
    end do
    call realloc(seed%spc,seed%nspc)
    if (seed%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    if (seed%nspc == 0) then
       errmsg = "No atomic species found."
       goto 999
    end if

    if (.not.ismol) then
       ! fill the cell metrics
       seed%m_x2c = transpose(r)
       r = seed%m_x2c
       call matinv(r,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting matrix"
          goto 999
       end if
       seed%useabr = 2

       ! convert atoms to crystallographic
       do i = 1, seed%nat
          seed%x(:,i) = matmul(r,seed%x(:,i))
       end do
    else
       seed%useabr = 0
    end if

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = ismol
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_xsf

  end subroutine read_xsf

  !> Read the structure from a pwc file.
  module subroutine read_pwc(seed,file,mol,errmsg,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, fclose, zatguess
    use param, only: isformat_pwc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    integer :: version, i, ier
    character*3, allocatable :: atm(:)
    real*8 :: r(3,3), alat

    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,form="unformatted",ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! header
    read (lu,err=999,end=999) version
    if (version < 2) then
       errmsg = "This pwc file is too old. Please update your QE and regenerate it."
       goto 999
    end if

    read (lu,err=999,end=999) seed%nspc, seed%nat, alat

    ! species
    allocate(atm(seed%nspc),seed%spc(seed%nspc))
    read (lu,err=999,end=999) atm
    do i = 1, seed%nspc
       seed%spc(i)%name = trim(atm(i))
       seed%spc(i)%z = zatguess(seed%spc(i)%name)
    end do
    deallocate(atm)

    ! read the rest
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    read (lu,err=999,end=999) seed%is
    read (lu,err=999,end=999) seed%x
    read (lu,err=999,end=999) seed%m_x2c
    seed%m_x2c = seed%m_x2c * alat

    ! convert to crystallographic
    r = seed%m_x2c
    call matinv(r,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if
    do i = 1, seed%nat
       seed%x(:,i) = matmul(r,seed%x(:,i)) * alat
    end do
    seed%useabr = 2

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_pwc

  end subroutine read_pwc

  !> Read the structure from an axsf file (xcrysden). Read the
  !> coordinates from PRIMCOORD block and nudge them using the
  !> eigenvector on the same block scaled by the value of xnudge
  !> (bohr).
  module subroutine read_axsf(seed,file,nread0,xnudge,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isinteger, &
       string, getword, isreal, nameguess, zatguess
    use tools_math, only: matinv
    use param, only: bohrtoa, isformat_axsf
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: nread0
    real*8, intent(in) :: xnudge
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, name
    integer :: lu, lp, iprim, i, it, iz, j, ier
    real*8 :: r(3,3), x(3)
    logical :: ok, ismol, didreadr, didreadx

    ! open
    call seed%end()
    ismol = .false.
    errmsg = ""
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    errmsg = "Error reading file: " // trim(file)
    didreadr = .false.
    didreadx = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"primvec")) then
          didreadr = .true.
          do i = 1, 3
             read (lu,*,err=999,end=999) r(i,:)
          end do
          r = r / bohrtoa
          ismol = .false.
       elseif (equal(word,"primcoord")) then
          ok = isinteger(iprim,line,lp)
          if (iprim == nread0) then
             didreadx = .true.
             read (lu,*,err=999,end=999) seed%nat
             allocate(seed%x(3,seed%nat),seed%is(seed%nat))
             seed%nspc = 0
             allocate(seed%spc(2))
             do i = 1, seed%nat
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999

                ! read the atomic coordinates
                lp = 1
                ok = isinteger(iz,line,lp)
                if (ok) then
                   ! Z x y z
                   name = nameguess(iz,.true.)
                else
                   word = getword(line,lp)
                   name = trim(adjustl(word))
                   iz = zatguess(name)
                end if
                ok = isreal(seed%x(1,i),line,lp)
                ok = ok.and.isreal(seed%x(2,i),line,lp)
                ok = ok.and.isreal(seed%x(3,i),line,lp)
                if (.not.ok) then
                   errmsg = 'Wrong atomic position.'
                   goto 999
                end if
                seed%x(:,i) = seed%x(:,i) / bohrtoa

                ! file this species if it is a new species
                it = 0
                do j = 1, seed%nspc
                   if (seed%spc(j)%z == iz) then
                      it = j
                      exit
                   end if
                end do
                if (it == 0) then
                   seed%nspc = seed%nspc + 1
                   if (seed%nspc > size(seed%spc,1)) &
                      call realloc(seed%spc,2*seed%nspc)
                   seed%spc(seed%nspc)%z = iz
                   seed%spc(seed%nspc)%name = name
                   it = seed%nspc
                end if
                seed%is(i) = it

                ! read the displacement vector and apply the nudge
                ok = isreal(x(1),line,lp)
                ok = ok.and.isreal(x(2),line,lp)
                ok = ok.and.isreal(x(3),line,lp)
                if (.not.ok) then
                   errmsg = 'Wrong displacement vector.'
                   goto 999
                end if
                seed%x(:,i) = seed%x(:,i) + xnudge * x
             end do
          end if
       end if
    end do
    if (.not.didreadr) then
       errmsg = "Could not find PRIMVEC block "
       goto 999
    end if
    if (.not.didreadx) then
       errmsg = "Could not find PRIMCOORD block number " // string(nread0)
       goto 999
    end if
    call realloc(seed%spc,seed%nspc)
    if (seed%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    if (seed%nspc == 0) then
       errmsg = "No atomic species found."
       goto 999
    end if

    if (.not.ismol) then
       ! fill the cell metrics
       seed%m_x2c = transpose(r)
       r = seed%m_x2c
       call matinv(r,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting matrix"
          goto 999
       end if
       seed%useabr = 2

       ! convert atoms to crystallographic
       do i = 1, seed%nat
          seed%x(:,i) = matmul(r,seed%x(:,i))
       end do
    else
       seed%useabr = 0
    end if

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = ismol
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_axsf

  end subroutine read_axsf

  !> Read an FHIaims input file.
  module subroutine read_aimsin(seed,file,mol,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isreal, &
       getword, zatguess
    use tools_math, only: matinv
    use types, only: realloc
    use hashmod, only: hash
    use param, only: bohrtoa, isformat_aimsin
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nlat, i, idx, ier
    logical :: is_file_mol, ok
    character(len=:), allocatable :: line, word
    real*8 :: rlat(3,3)
    logical, allocatable :: isfrac(:)
    type(hash) :: usespc

    ! open
    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)
    call usespc%init()

    ! allocate atoms
    allocate(isfrac(10),seed%x(3,10),seed%is(10))

    ! collect the information
    is_file_mol = .true.
    nlat = 0
    seed%nspc = 0
    seed%nat = 0
    do while (getline_raw(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (len_trim(word) == 0) cycle
       if (word(1:1) == "#") cycle

       if (equal(word,'atom').or.equal(word,'atom_frac')) then
          seed%nat = seed%nat + 1
          if (seed%nat > size(isfrac,1)) then
             call realloc(isfrac,2*seed%nat)
             call realloc(seed%x,3,2*seed%nat)
             call realloc(seed%is,2*seed%nat)
          end if
          ok = isreal(seed%x(1,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(2,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(3,seed%nat),line,lp)
          if (.not.ok) goto 999
          isfrac(seed%nat) = equal(word,'atom_frac')

          word = getword(line,lp)
          if (.not.usespc%iskey(word)) then
             seed%nspc = seed%nspc + 1
             call usespc%put(word,seed%nspc)
             seed%is(seed%nat) = seed%nspc
          else
             seed%is(seed%nat) = usespc%get(word,1)
          end if

       elseif (equal(word,'lattice_vector')) then
          is_file_mol = .false.
          nlat = nlat + 1
          ok = isreal(rlat(1,nlat),line,lp)
          ok = ok .and. isreal(rlat(2,nlat),line,lp)
          ok = ok .and. isreal(rlat(3,nlat),line,lp)
          if (.not.ok) then
             errmsg = "error reading lattice vectors from input file"
             goto 999
          end if
       end if
    end do
    if (.not.is_file_mol.and. nlat /= 3) then
       errmsg = "lattice_vector found but wrong number of lattice vectors"
       goto 999
    end if
    call realloc(isfrac,seed%nat)
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)

    ! fill the species array
    allocate(seed%spc(seed%nspc))
    do i = 1, seed%nspc
       word = usespc%getkey(i)
       idx = usespc%get(word,1)
       seed%spc(idx)%name = word
       seed%spc(idx)%z = zatguess(word)
       if (seed%spc(idx)%z <= 0) then
          errmsg = "unknown atom type: " // word
          goto 999
       end if
    end do

    ! handle the molecule/crystal expectation/contents of the file
    if (is_file_mol) then
       seed%ismolecule = .true.
       seed%useabr = 0
       seed%m_x2c = 0d0
    else
       if (mol) then
          errmsg = "tried to load a crystal as a molecule; use the CRYSTAL keyword"
          goto 999
       end if
       seed%ismolecule = .false.
       seed%useabr = 2
       rlat = rlat / bohrtoa
       seed%m_x2c = rlat
       call matinv(rlat,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting lattice vector matrix"
          goto 999
       end if
    end if

    ! convert the atomic coordinates
    do i = 1, seed%nat
       if (.not.isfrac(i)) then
          seed%x(:,i) = seed%x(:,i) / bohrtoa
          if (.not.is_file_mol) seed%x(:,i) = matmul(rlat,seed%x(:,i))
       end if
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_aimsin

  end subroutine read_aimsin

  !> Read an FHI aims output file.
  module subroutine read_aimsout(seed,file,mol,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isreal, &
       getword, zatguess
    use tools_math, only: matinv
    use types, only: realloc
    use hashmod, only: hash
    use param, only: bohrtoa, hartoev, eva3togpa, isformat_aimsout
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nlat, i, j, idx, ier, nupdate, iup, iat
    logical :: is_file_mol, ok, isfinal
    character*1 :: cdum
    character*10 :: dum1, dum2, splbl
    character(len=:), allocatable :: line, word
    real*8 :: rlat(3,3)
    logical, allocatable :: isfrac(:)
    type(hash) :: usespc

    ! open
    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)
    call usespc%init()

    ! allocate atoms
    allocate(isfrac(10),seed%x(3,10),seed%is(10))

    ! advance until we are at the "Input geometry" line
    ok = .false.
    do while (getline_raw(lu,line))
       if (line == "  Input geometry:") then
          ok = .true.
          exit
       end if
    end do
    if (.not.ok) then
       errmsg = "Failed to locate the Input geometry block in the FHIaims output file"
       goto 999
    end if

    ! get the cell geometry
    ok = getline_raw(lu,line)
    if (.not.ok) goto 999
    if (index(line,"Unit cell:") > 0) then
       do i = 1, 3
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
          read (line,*,err=999,end=999) cdum, (rlat(j,i),j=1,3)
       end do
       is_file_mol = .false.
    elseif (index(line,"No unit cell requested.") > 0) then
       is_file_mol = .true.
    else
       goto 999
    end if

    ! get the atomic positions
    seed%nspc = 0
    seed%nat = 0
    ok = getline_raw(lu,line)
    ok = ok .and. getline_raw(lu,line)
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) exit

       seed%nat = seed%nat + 1
       if (seed%nat > size(isfrac,1)) then
          call realloc(isfrac,2*seed%nat)
          call realloc(seed%x,3,2*seed%nat)
          call realloc(seed%is,2*seed%nat)
       end if
       read(line,*,err=999,end=999) cdum, dum1, dum2, splbl, (seed%x(j,seed%nat),j=1,3)
       isfrac(seed%nat) = .false.

       word = trim(adjustl(splbl))
       if (.not.usespc%iskey(word)) then
          seed%nspc = seed%nspc + 1
          call usespc%put(word,seed%nspc)
          seed%is(seed%nat) = seed%nspc
       else
          seed%is(seed%nat) = usespc%get(word,1)
       end if
    end do
    call realloc(isfrac,seed%nat)
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)

    ! fill the species array
    allocate(seed%spc(seed%nspc))
    do i = 1, seed%nspc
       word = usespc%getkey(i)
       idx = usespc%get(word,1)
       seed%spc(idx)%name = word
       seed%spc(idx)%z = zatguess(word)
       if (seed%spc(idx)%z <= 0) then
          errmsg = "unknown atom type: " // word
          goto 999
       end if
    end do

    ! read the rest of the file and search for "Updated atomic structure"
    ! or "Final atomic structure" blocks
    nupdate = 0
    isfinal = .false.
    do while (getline_raw(lu,line))
       if (index(line,'| Total energy uncorrected') > 0) then
          idx = index(line,':')
          ok = isreal(seed%energy,line(idx+1:))
          if (ok) then
             seed%energy = seed%energy / hartoev
          else
             seed%energy = huge(1d0)
          end if
       elseif (index(line,'|  Pressure') > 0) then
          idx = index(line,':')
          ok = isreal(seed%pressure,line(idx+1:))
          if (ok) then
             seed%pressure = seed%pressure * eva3togpa
          else
             seed%pressure = huge(1d0)
          end if
       elseif (trim(line) == "  Updated atomic structure:") then
          nupdate = nupdate + 1
       elseif (trim(line) == "  Final atomic structure:") then
          isfinal = .true.
          exit
       end if
    end do

    ! if the "final atomic structure" is not found but "updated atomic structure"
    ! was, this must be an aborted run; read the last geometry
    if (.not.isfinal.and.nupdate > 0) then
       rewind(lu)
       iup = 0
       do while (getline_raw(lu,line))
          if (trim(line) == "  Updated atomic structure:") then
             iup = iup + 1
             if (iup == nupdate) exit
          end if
       end do
    end if

    ! read the last geometry block
    if (isfinal.or.nupdate > 0) then
       nlat = 0
       iat = 0
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       do while (getline_raw(lu,line))
          if (line == "  Fractional coordinates:" .or. line(1:4) == "----") exit
          if (len_trim(line) == 0) cycle
          lp = 1
          word = lgetword(line,lp)
          if (word == "lattice_vector") then
             nlat = nlat + 1
             ok = isreal(rlat(1,nlat),line,lp)
             ok = ok .and. isreal(rlat(2,nlat),line,lp)
             ok = ok .and. isreal(rlat(3,nlat),line,lp)
             if (.not.ok) goto 999
          elseif (word == "atom") then
             iat = iat + 1
             ok = isreal(seed%x(1,iat),line,lp)
             ok = ok .and. isreal(seed%x(2,iat),line,lp)
             ok = ok .and. isreal(seed%x(3,iat),line,lp)
             isfrac(iat) = .false.
          end if
       end do
    end if

    ! handle the molecule/crystal expectation/contents of the file
    if (is_file_mol) then
       seed%ismolecule = .true.
       seed%useabr = 0
       seed%m_x2c = 0d0
    else
       if (mol) then
          errmsg = "tried to load a crystal as a molecule; use the CRYSTAL keyword"
          goto 999
       end if
       seed%ismolecule = .false.
       seed%useabr = 2
       rlat = rlat / bohrtoa
       seed%m_x2c = rlat
       call matinv(rlat,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting lattice vector matrix"
          goto 999
       end if
    end if

    ! convert the atomic coordinates
    do i = 1, seed%nat
       if (.not.isfrac(i)) then
          seed%x(:,i) = seed%x(:,i) / bohrtoa
          if (.not.is_file_mol) seed%x(:,i) = matmul(rlat,seed%x(:,i))
       end if
    end do

    ! read the last energy
    do while (getline_raw(lu,line))
       if (index(line,'| Total energy uncorrected') > 0) then
          idx = index(line,':')
          ok = isreal(seed%energy,line(idx+1:))
          if (ok) then
             seed%energy = seed%energy / hartoev
          else
             seed%energy = huge(1d0)
          end if
       end if
    end do

    errmsg = ""
999 continue
    call fclose(lu)

    ! symmetry
    seed%havesym = 0
    seed%findsym = -1
    seed%checkrepeats = .false.

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = file
    seed%name = file
    seed%isformat = isformat_aimsout

  end subroutine read_aimsout

  !> Read the structure from a TINKER frac file (must have cell
  !> parameters in the second line).
  module subroutine read_tinkerfrac(seed,file,mol,errmsg,ti)
    use tools_io, only: getline_raw, fopen_read, fclose, isinteger, isreal, zatguess, nameguess
    use param, only: bohrtoa, maxzat, isformat_tinkerfrac
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, i, iz, idum
    logical :: ok
    character(len=:), allocatable :: line
    character*30 :: atsym
    integer, allocatable :: imap(:)

    ! open the file
    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    seed%file = file
    seed%isformat = isformat_tinkerfrac

    ! line 1: number of atoms and name of the seed
    lp = 1
    ok = getline_raw(lu,line,.false.)
    ok = ok .and. isinteger(seed%nat,line,lp)
    if (.not.ok) goto 999
    seed%name = seed%file

    ! line 2: cell parameters
    read (lu,*,err=999,end=999) seed%aa, seed%bb
    seed%aa = seed%aa / bohrtoa
    seed%useabr = 1

    ! atoms
    seed%nspc = 0
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    allocate(seed%spc(2),imap(maxzat))
    imap = 0
    do i = 1, seed%nat
       read(lu,*,err=999,end=999) idum, atsym, seed%x(:,i)
       iz = zatguess(atsym)
       if (iz <= 0 .or. iz > maxzat) then
          errmsg = "Unknown atomic symbol."
          goto 999
       endif
       if (imap(iz) == 0) then
          seed%nspc = seed%nspc + 1
          if (seed%nspc > size(seed%spc,1)) call realloc(seed%spc,2*seed%nspc)
          seed%spc(seed%nspc)%name = nameguess(iz,.true.)
          seed%spc(seed%nspc)%z = iz
          seed%spc(seed%nspc)%qat = 0d0
          imap(iz) = seed%nspc
       endif
       seed%is(i) = imap(iz)
    end do
    call realloc(seed%spc,seed%nspc)
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    deallocate(imap)

    errmsg = ""
999 continue
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%checkrepeats = .false.
    seed%findsym = -1

    ! molecule
    seed%ismolecule = mol
    seed%havex0 = .true.
    seed%molx0 = 0d0

    ! rest of the seed information
    seed%isused = .true.
    seed%cubic = .false.
    seed%border = 0d0

  end subroutine read_tinkerfrac

  !> Adapt the size of an allocatable 1D type(crystalseed) array
  module subroutine realloc_crystalseed(a,nnew)

    type(crystalseed), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(crystalseed), allocatable :: temp(:)
    integer :: l1, u1

    if (.not.allocated(a)) return
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_crystalseed

  !> Detect the format for the structure-containing file. Normally,
  !> this works by detecting the extension, but the file may be
  !> opened and searched if ambiguity is present. The format and
  !> whether the file contains a molecule or crysatl is returned.
  !> If alsofield is present, then return .true. if the file also
  !> contains a scalar field.
  module subroutine struct_detect_format(file,isformat,alsofield,ti)
    use param, only: isformat_unknown, isformat_cif, isformat_shelx,&
       isformat_f21, isformat_xyz, isformat_gjf,&
       isformat_cube, isformat_bincube, isformat_struct, isformat_abinit, isformat_elk,&
       isformat_wfn, isformat_wfx, isformat_fchk, isformat_molden,&
       isformat_gaussian, isformat_siesta, isformat_xsf, isformat_gen,&
       isformat_vasp, isformat_pwc, isformat_axsf, isformat_dat, isformat_pgout,&
       isformat_dmain, isformat_aimsin, isformat_aimsout, isformat_tinkerfrac,&
       isformat_castepcell, isformat_castepgeom, isformat_qein, isformat_qeout,&
       isformat_mol2, isformat_pdb, isformat_zmat
    use tools_io, only: equal, fopen_read, fclose, lower, getline,&
       getline_raw, equali
    use param, only: dirsep
    character*(*), intent(in) :: file
    integer, intent(out) :: isformat
    logical, intent(out), optional :: alsofield
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: basename, wextdot, wextdot2, wext_, line
    logical :: isvasp, alsofield_
    integer :: lu, nat, ios, idx
    character*1 :: isfrac

    if (present(alsofield)) alsofield = .false.
    alsofield_ = .false.
    basename = file(index(file,dirsep,.true.)+1:)
    wext_ = basename(index(basename,'_',.true.)+1:)

    idx = index(basename,'.',.true.)
    wextdot = basename(idx+1:)
    if (idx > 0) then
       wextdot2 = basename(index(basename(1:idx-1),'.',.true.)+1:)
    else
       wextdot2 = ""
    end if

    isvasp = (index(basename,'CONTCAR') > 0) .or. &
       (index(basename,'CHGCAR') > 0) .or. (index(basename,'CHG') > 0).or.&
       (index(basename,'ELFCAR') > 0) .or. (index(basename,'AECCAR0') > 0).or.&
       (index(basename,'AECCAR1') > 0) .or. (index(basename,'AECCAR2') > 0) .or.&
       (index(basename,'POSCAR') > 0)

    if (equal(lower(wextdot),'cif')) then
       isformat = isformat_cif
    elseif (equal(wextdot,'pwc')) then
       isformat = isformat_pwc
       alsofield_ = .true.
    elseif (equal(wextdot,'res').or.equal(wextdot,'ins').or.equal(wextdot,'16')) then
       isformat = isformat_shelx
    elseif (equal(wextdot,'21')) then
       isformat = isformat_f21
    elseif (equal(wextdot,'cube')) then
       isformat = isformat_cube
       alsofield_ = .true.
    elseif (equal(wextdot,'bincube')) then
       isformat = isformat_bincube
       alsofield_ = .true.
    elseif (equal(wextdot,'struct')) then
       isformat = isformat_struct
    elseif (equal(wextdot,'DEN').or.equal(wext_,'DEN').or.equal(wextdot,'ELF').or.equal(wext_,'ELF').or.&
       equal(wextdot,'POT').or.equal(wext_,'POT').or.equal(wextdot,'VHA').or.equal(wext_,'VHA').or.&
       equal(wextdot,'VHXC').or.equal(wext_,'VHXC').or.equal(wextdot,'VXC').or.equal(wext_,'VXC').or.&
       equal(wextdot,'GDEN1').or.equal(wext_,'GDEN1').or.equal(wextdot,'GDEN2').or.equal(wext_,'GDEN2').or.&
       equal(wextdot,'GDEN3').or.equal(wext_,'GDEN3').or.equal(wextdot,'LDEN').or.equal(wext_,'LDEN').or.&
       equal(wextdot,'KDEN').or.equal(wext_,'KDEN').or.equal(wextdot,'PAWDEN').or.equal(wext_,'PAWDEN').or.&
       equal(wextdot,'VCLMB').or.equal(wext_,'VCLMB').or.equal(wextdot,'VPSP').or.equal(wext_,'VPSP')) then
       isformat = isformat_abinit
       alsofield_ = .true.
    elseif (equal(wextdot,'OUT')) then
       isformat = isformat_elk
    elseif (equal(wextdot,'out')) then
       call which_out_format(file,isformat,ti=ti)
    elseif (equal(wextdot,'own')) then
       call which_out_format(file,isformat,ti=ti)
       if (isformat /= isformat_aimsout) goto 999
    elseif (equal(wextdot,'in')) then
       call which_in_format(file,isformat,ti=ti)
    elseif (equal(wextdot2,'in.next_step')) then
       call which_in_format(file,isformat,ti=ti)
       if (isformat /= isformat_aimsin) goto 999
    elseif (equal(wextdot,'pwi')) then
       isformat = isformat_qein
    elseif (equal(wextdot,'pwo')) then
       isformat = isformat_qeout
    elseif (equal(wextdot,'xyz')) then
       isformat = isformat_xyz
    elseif (equal(wextdot,'gjf').or.equal(wextdot,'com')) then
       isformat = isformat_gjf
    elseif (equal(wextdot,'zmat')) then
       isformat = isformat_zmat
    elseif (equal(wextdot,'pgout')) then
       isformat = isformat_pgout
    elseif (equal(wextdot,'wfn')) then
       isformat = isformat_wfn
       alsofield_ = .true.
    elseif (equal(wextdot,'wfx')) then
       isformat = isformat_wfx
       alsofield_ = .true.
    elseif (equal(wextdot,'log')) then
       isformat = isformat_gaussian
       alsofield_ = .false.
    elseif (equal(wextdot,'fchk')) then
       isformat = isformat_fchk
       alsofield_ = .true.
    elseif (equal(wextdot,'molden') .or. equal(wextdot2,'molden.input')) then
       isformat = isformat_molden
       alsofield_ = .true.
    elseif (equal(wextdot,'STRUCT_OUT').or.equal(wextdot,'STRUCT_IN')) then
       isformat = isformat_siesta
    elseif (equal(wextdot,'cell')) then
       isformat = isformat_castepcell
    elseif (equal(wextdot,'geom')) then
       isformat = isformat_castepgeom
    elseif (equal(wextdot,'dmain')) then
       isformat = isformat_dmain
    elseif (equal(wextdot,'xsf')) then
       isformat = isformat_xsf
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) goto 999
       do while (getline(lu,line))
          if (len_trim(line) > 0) exit
       end do
       if (present(alsofield)) then
          do while (getline(lu,line))
             if (equali(line,"begin_block_datagrid_3d")) then
                alsofield_ = .true.
                exit
             end if
          end do
       end if
       call fclose(lu)
    elseif (equal(wextdot,'gen')) then
       isformat = isformat_gen

       ! determine whether it is a molecule or crystal
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) goto 999
       do while (getline_raw(lu,line))
          if (len_trim(line) > 0) exit
       end do
       read (line,*,iostat=ios) nat, isfrac
       if (ios /= 0) goto 999
       isfrac = lower(isfrac)
       call fclose(lu)
    elseif (equal(wextdot,'axsf')) then
       isformat = isformat_axsf
    elseif (equal(wextdot,'dat')) then
       isformat = isformat_dat
    elseif (equal(wextdot,'frac')) then
       isformat = isformat_tinkerfrac
    elseif (equal(lower(wextdot),'vasp')) then
       isformat = isformat_vasp
       alsofield_ = .false.
    elseif (isvasp) then
       isformat = isformat_vasp
       alsofield_ = (index(basename,'CHGCAR') > 0) .or. (index(basename,'CHG') > 0) .or. &
          (index(basename,'ELFCAR') > 0) .or. (index(basename,'AECCAR0') > 0) .or. &
          (index(basename,'AECCAR1') > 0) .or. (index(basename,'AECCAR2') > 0)
    elseif (equal(wextdot,'mol2')) then
       isformat = isformat_mol2
    elseif (equal(wextdot,'pdb')) then
       isformat = isformat_pdb
    else
       goto 999
    endif
    if (present(alsofield)) alsofield = alsofield_

    if (isformat == 0) goto 999

    return
999 continue
    isformat = isformat_unknown

  end subroutine struct_detect_format

  !> Detect the format for a file containing molecular or
  !> crystal vibrations.
  module subroutine vibrations_detect_format(file,ivformat,ti)
    use param, only: isformat_unknown, isformat_v_matdynmodes, dirsep
    use tools_io, only: equal, lower
    character*(*), intent(in) :: file
    integer, intent(out) :: ivformat
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: basename, wextdot, wextdot2, wext_
    integer :: idx

    basename = file(index(file,dirsep,.true.)+1:)
    wext_ = basename(index(basename,'_',.true.)+1:)
    idx = index(basename,'.',.true.)
    wextdot = basename(idx+1:)
    if (idx > 0) then
       wextdot2 = basename(index(basename(1:idx-1),'.',.true.)+1:)
    else
       wextdot2 = ""
    end if

    if (equal(lower(wextdot),'modes')) then
       ivformat = isformat_v_matdynmodes
    else
       ivformat = isformat_unknown
    endif

  end subroutine vibrations_detect_format

  !> Detect whether a file with format isformat contains a molecule
  !> (ismol=.true.)  or a crystal (.false.)
  module subroutine struct_detect_ismol(file,isformat,ismol,ti)
    use tools_io, only: fopen_read, fclose, getline, equali,&
       getline_raw, lgetword, equal, lower
    use param, only: isformat_cif, isformat_shelx, isformat_f21,&
       isformat_cube, isformat_bincube, isformat_struct, isformat_abinit, isformat_elk,&
       isformat_qein, isformat_qeout, isformat_crystal, isformat_fploout,&
       isformat_xyz, isformat_gjf, isformat_wfn,&
       isformat_wfx, isformat_fchk, isformat_molden, isformat_gaussian, isformat_siesta,&
       isformat_xsf, isformat_gen, isformat_vasp, isformat_pwc, isformat_axsf,&
       isformat_dat, isformat_pgout, isformat_orca, isformat_dmain, isformat_aimsin,&
       isformat_aimsout, isformat_tinkerfrac, isformat_castepcell, isformat_castepgeom,&
       isformat_mol2, isformat_pdb, isformat_zmat
    character*(*), intent(in) :: file
    integer, intent(in) :: isformat
    logical, intent(out) :: ismol
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word
    integer :: lu, ios, lp, nat
    character*1 :: isfrac
    logical :: ok

    ismol = .false.
    select case (isformat)
    case (isformat_cif,isformat_pwc,isformat_shelx,isformat_f21,&
       isformat_cube,isformat_bincube,isformat_struct,isformat_abinit,&
       isformat_elk,isformat_siesta,isformat_dmain,isformat_vasp,&
       isformat_axsf,isformat_tinkerfrac,isformat_qein,isformat_qeout,&
       isformat_crystal,isformat_fploout,isformat_castepcell,isformat_castepgeom)
       ismol = .false.

    case (isformat_xyz,isformat_gjf,isformat_pgout,isformat_wfn,isformat_wfx,&
       isformat_gaussian,isformat_fchk,isformat_molden,isformat_dat,&
       isformat_orca,isformat_mol2,isformat_pdb,isformat_zmat)
       ismol = .true.

    case(isformat_aimsout)
       lu = fopen_read(file,errstop=.false.,ti=ti)
       ismol = .false.
       do while(getline_raw(lu,line))
          if (adjustl(trim(line)) == "Input geometry:") then
             ok = getline_raw(lu,line)
             if (.not.ok) then
                return
             else if (index(line,"Unit cell:") > 0) then
                ismol = .false.
             elseif (index(line,"No unit cell requested.") > 0) then
                ismol = .true.
             endif
             exit
          end if
       end do
       call fclose(lu)

    case (isformat_aimsin)
       ismol = .true.
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) return
       do while(getline_raw(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,"lattice_vector")) then
             ismol = .false.
             exit
          end if
       end do
       call fclose(lu)

    case (isformat_gen)
       ismol = .false.
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) return
       do while (getline_raw(lu,line))
          if (len_trim(line) > 0) exit
       end do
       read (line,*,iostat=ios) nat, isfrac
       if (ios /= 0) return
       isfrac = lower(isfrac)
       if (equal(isfrac,"c")) then
          ismol = .true.
       else
          ismol = .false.
       end if
       call fclose(lu)

    case (isformat_xsf)
       ismol = .false.
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) return
       do while (getline(lu,line))
          if (len_trim(line) > 0) exit
       end do
       if (equali(line,"atoms")) then
          ismol = .true.
       else
          ismol = .false.
       end if
       call fclose(lu)

    case default
       ismol = .false.
    end select

  end subroutine struct_detect_ismol

  !> Read the species into the seed from a VASP POTCAR file.
  module subroutine read_potcar(seed,file,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, getword, fclose, zatguess
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp
    character(len=:), allocatable :: aux1, aatom, line
    logical :: ok

    errmsg = ""
    seed%nspc = 0
    if (allocated(seed%spc)) deallocate(seed%spc)
    allocate(seed%spc(2))

    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening POTCAR file."
       return
    end if

    ! read the atoms
    do while (getline_raw(lu,line))
       lp = 1
       aux1 = getword(line,lp)
       aatom = getword(line,lp)
       seed%nspc = seed%nspc + 1
       if (seed%nspc > size(seed%spc,1)) &
          call realloc(seed%spc,2*seed%nspc)

       seed%spc(seed%nspc)%name = aatom
       seed%spc(seed%nspc)%z = zatguess(aatom)
       line = ""
       do while (.not. (trim(adjustl(line)) == 'End of Dataset'))
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) then
             errmsg = "Unexpected termination of POTCAR file."
             call fclose(lu)
             return
          end if
       end do
    end do
    call realloc(seed%spc,seed%nspc)

    ! close
    call fclose(lu)

  end subroutine read_potcar

  !> Read all seeds from a file. If iafield is present, then return
  !> the seed number for which the file can be read as a field (or 0
  !> if none). If mol0 == 1, force reading molecules, if mol0 == 0,
  !> force reading crystals, if mol0 == -1, let the routine figure it
  !> out. Returns the number of seeds read (nseed), the seeds themselves
  !> (seed) and collapse=.true. if the seeds are the steps of a
  !> geometry optimization.
  module subroutine read_seeds_from_file(file,mol0,isformat0,readlastonly,&
     nseed,seed,collapse,errmsg,iafield,ti)
    use global, only: rborder_def, doguess
    use tools_io, only: getword, equali, fopen_read, fclose
    use param, only: isformat_cube, isformat_bincube, isformat_xyz, isformat_wfn,&
       isformat_wfx, isformat_fchk, isformat_molden, isformat_gaussian, isformat_gjf,&
       isformat_zmat, isformat_abinit, isformat_cif, isformat_pwc, isformat_fploout,&
       isformat_crystal, isformat_elk, isformat_gen, isformat_qein, isformat_qeout,&
       isformat_shelx, isformat_siesta, isformat_struct, isformat_vasp, isformat_axsf,&
       isformat_xsf, isformat_castepcell, isformat_castepgeom,&
       isformat_dat, isformat_f21, isformat_unknown, isformat_pgout, isformat_orca,&
       isformat_dmain, isformat_aimsin, isformat_aimsout, isformat_tinkerfrac,&
       isformat_mol2, isformat_pdb
    character*(*), intent(in) :: file
    integer, intent(in) :: mol0
    integer, intent(in) :: isformat0
    logical, intent(in) :: readlastonly
    integer, intent(out) :: nseed
    type(crystalseed), allocatable, intent(inout) :: seed(:)
    logical, intent(out) :: collapse
    character(len=:), allocatable, intent(out) :: errmsg
    integer, intent(out), optional :: iafield
    type(thread_info), intent(in), optional :: ti

    integer :: isformat, is, mol0_, i, lu
    logical :: mol, alsofield, ok

    errmsg = ""
    alsofield = .false.
    mol0_ = mol0
    isformat = isformat0
    nseed = 0
    if (allocated(seed)) deallocate(seed)

    ! check if the file exists and is openable
    inquire(file=file,exist=ok)
    if (.not.ok) then
       errmsg = "File does not exist."
       goto 999
    end if
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       goto 999
    end if
    call fclose(lu)

    call struct_detect_format(file,is,alsofield,ti=ti)
    if (isformat == isformat_unknown) isformat = is
    if (isformat == isformat_unknown) then
       errmsg = "Unknown file format/extension in file."
       goto 999
    end if
    if (mol0_ == 1) then
       mol = .true.
    elseif (mol0_ == 0) then
       mol = .false.
    elseif (mol0_ == -1) then
       call struct_detect_ismol(file,isformat,mol,ti=ti)
    end if

    ! by default, we expect one seed only
    nseed = 1
    allocate(seed(1))

    ! read all available seeds in the file
    if (isformat == isformat_cif) then
       call read_all_cif(file,mol,errmsg,nseed=nseed,mseed=seed,ti=ti)
    elseif (isformat == isformat_mol2) then
       call read_all_mol2(file,errmsg,nseed=nseed,mseed=seed,ti=ti)
    elseif (isformat == isformat_pwc) then
       call seed(1)%read_pwc(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_shelx) then
       call seed(1)%read_shelx(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_f21) then
       call seed(1)%read_f21(file,mol,errmsg,ti=ti)
    else if (isformat == isformat_cube) then
       call seed(1)%read_cube(file,mol,errmsg,ti=ti)
    else if (isformat == isformat_bincube) then
       call seed(1)%read_bincube(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_struct) then
       call seed(1)%read_wien(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_vasp) then
       call read_all_vasp(nseed,seed,file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_abinit) then
       call seed(1)%read_abinit(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_elk) then
       call seed(1)%read_elk(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_qeout) then
       call read_all_qeout(nseed,seed,file,mol,-1,errmsg,ti=ti)
    elseif (isformat == isformat_crystal) then
       call seed(1)%read_crystalout(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_fploout) then
       call seed(1)%read_fploout(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_qein) then
       call seed(1)%read_qein(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_xyz) then
       call read_all_xyz(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_gaussian) then
       call read_all_log(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_wfn .or. isformat == isformat_wfx.or.&
       isformat == isformat_fchk.or.isformat == isformat_molden.or.&
       isformat == isformat_dat.or.isformat == isformat_pgout.or.&
       isformat == isformat_orca.or.isformat == isformat_gjf.or.&
       isformat == isformat_pdb.or.isformat == isformat_zmat) then
       call seed(1)%read_mol(file,isformat,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_siesta) then
       call seed(1)%read_siesta(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_castepcell) then
       call seed(1)%read_castep_cell(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_castepgeom) then
       call read_all_castep_geom(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_dmain) then
       call seed(1)%read_dmain(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_aimsin) then
       call seed(1)%read_aimsin(file,mol,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_aimsout) then
       call read_all_aimsout(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_tinkerfrac) then
       call seed(1)%read_tinkerfrac(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_axsf) then
       call seed(1)%read_axsf(file,1,0d0,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_xsf) then
       call seed(1)%read_xsf(file,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_gen) then
       call seed(1)%read_dftbp(file,rborder_def,.false.,errmsg,ti=ti)
    end if
    if (mol0 /= -1) &
       seed(1)%ismolecule = mol

999 continue
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

    ! handle the readlastonly option
    if (readlastonly .and. nseed > 1) then
       seed(1) = seed(nseed)
       nseed = 1
       call realloc_crystalseed(seed,1)
    end if

    ! handle the doguess option
    do i = 1, nseed
       if (.not.seed(i)%ismolecule) then
          if (doguess == 0) then
             seed(i)%havesym = 0
             seed(i)%findsym = 0
             seed(i)%checkrepeats = .false.
          elseif (doguess == 1 .and. seed(i)%havesym == 0) then
             seed(i)%findsym = 1
          else
             seed(i)%findsym = -1
          end if
       end if
    end do

    ! output iafield
    if (present(iafield)) then
       if (alsofield) then
          iafield = nseed
       else
          iafield = 0
       end if
    end if

    ! output collapse
    collapse = ((isformat == isformat_qeout .or. isformat == isformat_gaussian .or.&
       isformat == isformat_aimsout .or. isformat == isformat_castepgeom).and.nseed > 1)

  end subroutine read_seeds_from_file

  !> Report the contents of a crystalseed (debug only)
  module subroutine report(seed)
    use tools_io, only: uout, string
    class(crystalseed), intent(inout) :: seed

    integer :: i, j

    write (uout,'("isused = ",A)') string(seed%isused)
    write (uout,'("file = ",A)') string(seed%file)
    write (uout,'("name = ",A)') string(seed%name)
    write (uout,*)

    write (uout,'("## species (",A,")")') string(seed%nspc)
    do i = 1, seed%nspc
       write (uout,'(99(A,X))') string(seed%spc(i)%z), string(seed%name)
    end do
    write (uout,*)

    write (uout,'("## atomic positions (",A,")")') string(seed%nat)
    do i = 1, seed%nat
       write (uout,'(99(A,X))') string(seed%is(i)), &
          (string(seed%x(j,i),'f',decimal=10),j=1,3)
    end do
    write (uout,*)

    write (uout,'("## cell")')
    write (uout,'("useabr = ",A)') string(seed%useabr)
    if (seed%useabr == 1) then
       write (uout,'("aa = ",3(A,X))') (string(seed%aa(i),'f',decimal=8),i=1,3)
       write (uout,'("bb = ",3(A,X))') (string(seed%bb(i),'f',decimal=5),i=1,3)
    else
       write (uout,'("x2c = ",3(A,X))') (string(seed%m_x2c(1,i),'f',decimal=8),i=1,3)
       write (uout,'("x2c = ",3(A,X))') (string(seed%m_x2c(2,i),'f',decimal=8),i=1,3)
       write (uout,'("x2c = ",3(A,X))') (string(seed%m_x2c(3,i),'f',decimal=8),i=1,3)
    end if
    write (uout,*)

    write (uout,'("## symmetry")')
    write (uout,'("havesym = ",A)') string(seed%havesym)
    write (uout,'("findsym = ",A)') string(seed%findsym)
    write (uout,'("checkrepeats = ",A)') string(seed%checkrepeats)
    write (uout,'("neqv = ",A)') string(seed%neqv)
    write (uout,'("ncv = ",A)') string(seed%ncv)
    write (uout,*)

    write (uout,'("## extra fields")')
    write (uout,'("ismolecule = ",A)') string(seed%ismolecule)
    write (uout,'("cubic = ",A)') string(seed%cubic)
    write (uout,'("border = ",A)') string(seed%border,'f',decimal=4)
    write (uout,'("havex0 = ",A)') string(seed%havex0)
    write (uout,'("molx0 = ",3(A,X))') (string(seed%molx0(i),'f',decimal=4),i=1,3)
    write (uout,'("energy = ",A)') string(seed%energy,'e',decimal=10)
    write (uout,'("pressure = ",A)') string(seed%pressure,'e',decimal=2)

  end subroutine report

  !> Define the assignment operator for the crystal seed class.
  module subroutine assign_crystalseed(to,from)
    class(crystalseed), intent(out) :: to
    type(crystalseed), intent(in) :: from

    to%isused = from%isused
    to%file = from%file
    to%name = from%name
    to%isformat = from%isformat
    to%nat = from%nat
    if (allocated(from%x)) then
       to%x = from%x
    else
       if (allocated(to%x)) deallocate(to%x)
    end if
    if (allocated(from%is)) then
       to%is = from%is
    else
       if (allocated(to%is)) deallocate(to%is)
    end if
    to%nspc = from%nspc
    if (allocated(from%spc)) then
       to%spc = from%spc
    else
       if (allocated(to%spc)) deallocate(to%spc)
    end if
    to%useabr = from%useabr
    to%aa = from%aa
    to%bb = from%bb
    to%m_x2c = from%m_x2c
    to%havesym = from%havesym
    to%findsym = from%findsym
    to%checkrepeats = from%checkrepeats
    to%neqv = from%neqv
    to%ncv = from%ncv
    if (allocated(from%cen)) then
       to%cen = from%cen
    else
       if (allocated(to%cen)) deallocate(to%cen)
    end if
    if (allocated(from%rotm)) then
       to%rotm = from%rotm
    else
       if (allocated(to%rotm)) deallocate(to%rotm)
    end if
    to%ismolecule = from%ismolecule
    to%cubic = from%cubic
    to%border = from%border
    to%havex0 = from%havex0
    to%molx0 = from%molx0

  end subroutine assign_crystalseed

  !> Read the alat from a Quantum ESPRESSO output file (file) and
  !> return it in alat. If error, return non-zero errmsg. ti = thread
  !> info.
  module subroutine read_alat_from_qeout(file,alat,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, isreal, fclose
    character*(*), intent(in) :: file
    real*8, intent(out) :: alat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    integer :: lu, idx
    character(len=:), allocatable :: line

    ! open
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! read the alat
    alat = -1d0
    do while (getline_raw(lu,line))
       if (index(line,"lattice parameter (alat)") > 0) then
          idx = index(line,"=")
          ok = isreal(alat,line(idx+1:))
       end if
    end do
    if (alat < 0d0) then
       errmsg = "Error reading alat from file: " // trim(file)
       alat = 0d0
    end if
    call fclose(lu)

  end subroutine read_alat_from_qeout

  !xx! private subroutines

  !> Read one or all structures from a VASP POSCAR-style file
  !> (filename file) and return the corresponding crystal seeds in
  !> seed. If mol=.true., interpret the structure as a molecule
  !> (currently, this only sets the %ismolecule field). If an error
  !> condition is found, return the error message in errmsg
  !> (zero-length string if no error).
  subroutine read_all_vasp(nseed,seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, isreal, getword, zatguess,&
       isinteger
    use types, only: realloc
    use tools_math, only: det3sym, matinv
    use param, only: bohrtoa, dirsep, isformat_vasp
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, ier, nn
    integer :: i, j
    character(len=:), allocatable :: line, word, ofile, path
    logical :: ok, iscar
    real*8 :: scalex, scaley, scalez, scale
    real*8 :: rprim(3,3), gprim(3,3)
    real*8 :: omegaa

    ! open
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! initialize the seeds
    nseed = 0
    if (allocated(seed)) deallocate(seed)
    allocate(seed(10))

    do while (.true.)
       ! new seed
       nseed = nseed + 1
       if (nseed > size(seed,1)) call realloc_crystalseed(seed,2*nseed)

       ! generic error message
       errmsg = "Error reading file: " // trim(file)

       ! read the title and the scale line
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       if (len_trim(line) == 0) exit ! USPEX does not insert a blank line after
       lp = 1
       word = getword(line,lp)
       if (len_trim(word) > 0) then
          seed(nseed)%name = trim(file) // "|" // trim(word)
       else
          seed(nseed)%name = trim(file)
       end if

       ok = getline_raw(lu,line)
       if (.not.ok) exit
       lp = 1
       ok = isreal(scalex,line,lp)
       if (.not.ok) exit
       ok = ok .and. isreal(scaley,line,lp)
       if (.not.ok) then
          scale = scalex
          scalex = 1d0
          scaley = 1d0
          scalez = 1d0
       else
          ok = isreal(scalez,line,lp)
          if (.not.ok) exit
          scale = 1d0
       end if

       ! read the cell vectors and calculate the metric tensor
       do i = 1, 3
          read (lu,*,err=998,end=998) rprim(1,i), rprim(2,i), rprim(3,i)
       end do

       if (scale < 0d0) then
          gprim = matmul(transpose(rprim),rprim)
          omegaa = sqrt(det3sym(gprim))
          ! adjust the lengths to give the volume
          scale = (abs(scale) / abs(omegaa))**(1d0/3d0)
       end if
       rprim(1,:) = rprim(1,:) * scalex * scale
       rprim(2,:) = rprim(2,:) * scaley * scale
       rprim(3,:) = rprim(3,:) * scalez * scale
       rprim = rprim / bohrtoa
       gprim = matmul(transpose(rprim),rprim)
       omegaa = sqrt(det3sym(gprim))
       if (omegaa < 0d0) then
          errmsg = "Negative cell volume."
          exit
       end if
       seed(nseed)%m_x2c = rprim
       call matinv(rprim,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting matrix"
          exit
       end if
       seed(nseed)%useabr = 2

       ! For versions >= 5.2, a line indicating the atom types appears here
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       lp = 1
       word = getword(line,lp)
       if (zatguess(word) >= 0) then
          ! An atom name has been read -> read the rest of the line
          seed(nseed)%nspc = 0
          if (allocated(seed(nseed)%spc)) deallocate(seed(nseed)%spc)
          allocate(seed(nseed)%spc(2))
          do while (zatguess(word) >= 0)
             seed(nseed)%nspc = seed(nseed)%nspc + 1
             if (seed(nseed)%nspc > size(seed(nseed)%spc,1)) &
                call realloc(seed(nseed)%spc,2*seed(nseed)%nspc)
             seed(nseed)%spc(seed(nseed)%nspc)%name = word
             seed(nseed)%spc(seed(nseed)%nspc)%z = zatguess(word)
             word = getword(line,lp)
          end do
          call realloc(seed(nseed)%spc,seed(nseed)%nspc)
          ok = getline_raw(lu,line)
          if (.not.ok) exit
       else
          ! see if we can locate a POTCAR in the same path
          path = file(1:index(file,dirsep,.true.))
          if (len_trim(path) < 1) &
             path = "."
          ofile = trim(path) // "/POTCAR"
          inquire(file=ofile,exist=ok)
          if (.not.ok) then
             errmsg = "Atom types not found in POSCAR and no POTCAR in the same directory"
             exit
          end if

          call seed(nseed)%read_potcar(ofile,errmsg,ti=ti)
          if (len_trim(errmsg) > 0) exit
          if (seed(nseed)%nspc <= 0) then
             errmsg = "No atoms found in POTCAR."
             exit
          end if
       end if

       ! read number of atoms of each type
       lp = 1
       seed(nseed)%nat = 0
       allocate(seed(nseed)%is(10))
       do i = 1, seed(nseed)%nspc
          ok = isinteger(nn,line,lp)
          if (.not.ok) then
             errmsg = "Too many atom types"
             exit
          end if
          do j = seed(nseed)%nat+1, seed(nseed)%nat+nn
             if (j > size(seed(nseed)%is)) &
                call realloc(seed(nseed)%is,2*(seed(nseed)%nat+nn))
             seed(nseed)%is(j) = i
          end do
          seed(nseed)%nat = seed(nseed)%nat + nn
       end do
       allocate(seed(nseed)%x(3,seed(nseed)%nat))
       call realloc(seed(nseed)%is,seed(nseed)%nat)

       ! check there are no more atoms in this line
       nn = -1
       ok = isinteger(nn,line,lp)
       if (ok .and. nn /= -1) then
          errmsg = "Too few atom types"
          exit
       end if

       ! Read atomic positions (cryst. coords.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = adjustl(line)
       if (line(1:1) == 's' .or. line(1:1) == 'S') then
          ok = getline_raw(lu,line)
          if (.not.ok) exit
          line = adjustl(line)
       endif
       iscar = .false.
       if (line(1:1) == 'd' .or. line(1:1) == 'D') then
          iscar = .false.
       elseif (line(1:1) == 'c' .or. line(1:1) == 'C' .or. line(1:1) == 'k' .or. line(1:1) == 'K') then
          iscar = .true.
       endif
       do i = 1, seed(nseed)%nat
          read(lu,*,err=998,end=998) seed(nseed)%x(:,i)
          if (iscar) &
             seed(nseed)%x(:,i) = matmul(rprim,seed(nseed)%x(:,i) / bohrtoa)
       enddo
    end do
998 continue

    ! last seed is not valid, make sure we have at least one
    nseed = nseed - 1
    if (nseed == 0) then
       call fclose(lu)
       return
    else
       errmsg = ""
    end if
    call realloc_crystalseed(seed,nseed)

    call fclose(lu)

    ! fill the rest of the information
    do i = 1, nseed
       ! symmetry
       seed(i)%havesym = 0
       seed(i)%findsym = -1
       seed(i)%checkrepeats = .false.

       ! rest of the seed information
       seed(i)%isused = .true.
       seed(i)%ismolecule = mol
       seed(i)%cubic = .false.
       seed(i)%border = 0d0
       seed(i)%havex0 = .false.
       seed(i)%molx0 = 0d0
       seed(i)%file = file
       seed(i)%isformat = isformat_vasp
    end do

  end subroutine read_all_vasp

  !> Read multiple structure seeds from a CIF file. If mol, force a
  !> molecule/crystal system.  If error, return non-empty errmsg.
  !> If nseed and mseed are present, read all seeds into the mseed
  !> array and return the number of seeds in nseed. If seed0 is
  !> present, return the data block dblock. If dblock is not present
  !> or is empty, return the first data block in seed0.
  subroutine read_all_cif(file,mol,errmsg,nseed,mseed,seed0,dblock,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lower, isexpression_or_word,&
       isreal, zatguess
    use types, only: realloc
    use param, only: tab, isformat_cif
    integer, intent(out), optional :: nseed
    type(crystalseed), intent(inout), allocatable, optional :: mseed(:)
    type(crystalseed), intent(inout), optional :: seed0
    character*(*), intent(in), optional :: dblock
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    real*8 :: rdum, aa(3), bb(3)
    integer :: lu, idx, leng, lp, nhead, i
    character(len=:), allocatable :: word, line, blockname, spg
    logical :: indata, inloopheader, inloop, ldum, ok
    logical :: havefields(10) ! 1-3 = abc, 4-6 = angles, 7=symbol, 8-10=xyz
    integer :: looptype ! 1 = symmetry, 2 = atoms
    integer :: cols(5) ! 1=symbol, 2=x, 3=y, 4=z, 5=symop
    integer :: nat, nop
    real*8, allocatable :: xat(:,:)
    integer, allocatable :: zat(:)
    character*256, allocatable :: ops(:)
    type(crystalseed) :: seed

    ! consistency check
    if (.not.(present(nseed).and.present(mseed)).and..not.present(seed0)) then
       errmsg = "Error reading cif file (incorrect use of the interface)"
       return
    end if

    ! initialize
    errmsg = ""
    if (present(nseed).and.present(mseed)) then
       nseed = 0
       if (allocated(mseed)) deallocate(mseed)
       allocate(mseed(10))
    end if

    ! initialize atom list
    nat = 0
    nop = 0
    allocate(xat(3,20),zat(20))
    allocate(ops(48))

    ! run over lines in the file
    indata = .false.
    inloopheader = .false.
    inloop = .false.
    spg = ""
    havefields = .false.
    lu = fopen_read(file,ti=ti)
    main: do while (getline_raw(lu,line))
       line = adjustl(line)

       ! replace tabs with blanks
       do while (.true.)
          idx = index(line,tab)
          if (idx == 0) exit
          line = line(1:idx-1) // " " // line(idx+1:)
       end do

       ! skip blank lines and comments
       leng = len_trim(line)
       if (leng == 0) cycle
       if (line(1:1) == "#") cycle

       ! data_
       if (leng >= 5) then
          if (lower(line(1:5)) == "data_") then
             ! finalize the previous seeed
             if (indata) then
                call fill_seed(seed)
                if (present(nseed).and.present(mseed)) then
                   nseed = nseed + 1
                   if (nseed > size(mseed,1)) call realloc_crystalseed(mseed,2*nseed)
                   mseed(nseed) = seed
                elseif (present(seed0)) then
                   seed0 = seed
                   return
                else
                   errmsg = "Error reading cif file (incorrect use of the interface)"
                   return
                end if
             end if
             spg = ""
             indata = .false.
             blockname = trim(line(6:))
             havefields = .false.
             inloopheader = .false.
             inloop = .false.
             nat = 0
             nop = 0

             ! reset the flags and set the block name
             ok = .not.present(dblock)
             if (.not.ok) &
                ok = (len_trim(dblock) == 0 .or. blockname == dblock)
             if (ok) indata = .true.
             cycle
          end if
       end if
       if (.not.indata) cycle

       ! loop_
       if (leng >= 5) then
          if (lower(line(1:5)) == "loop_") then
             inloop = .false.
             inloopheader = .true.
             looptype = 0
             cols = 0
             nhead = 0
             cycle
          end if
       end if

       ! A _ keyword
       lp = 1
       if (line(1:1) == "_") then
          inloop = .false.
          ldum = isexpression_or_word(word,line,lp)
          word = lower(word)
          if (.not.inloopheader) then
             ! one of the variables not in a loop
             if (word == "_cell_length_a") then
                aa(1) = read_field_float()
                havefields(1) = .true.
             elseif (word == "_cell_length_b") then
                aa(2) = read_field_float()
                havefields(2) = .true.
             elseif (word == "_cell_length_c") then
                aa(3) = read_field_float()
                havefields(3) = .true.
             elseif (word == "_cell_angle_alpha") then
                bb(1) = read_field_float()
                havefields(4) = .true.
             elseif (word == "_cell_angle_beta") then
                bb(2) = read_field_float()
                havefields(5) = .true.
             elseif (word == "_cell_angle_gamma") then
                bb(3) = read_field_float()
                havefields(6) = .true.
             elseif (word == "_symmetry_space_group_name_h-m") then
                spg = read_field_string()
             elseif (word == "_space_group_name_h-m_alt") then
                if (len_trim(spg) == 0) spg = read_field_string()
             end if
             if (len_trim(errmsg) /= 0) return
          else
             nhead = nhead + 1
             if (word == "_atom_site_type_symbol" .or. word == "_atom_site_label" .or.&
                word == "_atom_site_fract_x" .or. word == "_atom_site_fract_y" .or.&
                word == "_atom_site_fract_z") then
                ! an atom loop
                looptype = 2
                nat = 0
                if (word == "_atom_site_type_symbol") then
                   cols(1) = nhead
                elseif (word == "_atom_site_label") then
                   if (cols(1) == 0) cols(1) = nhead
                elseif (word == "_atom_site_fract_x") then
                   cols(2) = nhead
                elseif (word == "_atom_site_fract_y") then
                   cols(3) = nhead
                elseif (word == "_atom_site_fract_z") then
                   cols(4) = nhead
                end if
             elseif (word == "_symmetry_equiv_pos_as_xyz" .or. &
                word == "_space_group_symop_operation_xyz") then
                ! a symmetry loop
                looptype = 1
                nop = 0
                if (word == "_symmetry_equiv_pos_as_xyz") then
                   cols(5) = nhead
                elseif (word == "_space_group_symop_operation_xyz") then
                   if (cols(5) == 0) cols(5) = nhead
                end if
             end if
          end if
          cycle
       else
          if (inloopheader) inloop = .true.
          inloopheader = .false.
       end if

       ! we must be in a loop, read or skip
       if (inloop) then
          if (looptype == 2) then
             ! the atoms loop
             ! add one atom
             nat = nat + 1
             if (nat > size(zat,1)) then
                call realloc(xat,3,2*nat)
                call realloc(zat,2*nat)
             end if

             ! read the fields
             lp = 1
             do i = 1, maxval(cols(1:4))
                ok = isexpression_or_word(word,line,lp)
                if (ok) then
                   if (cols(2)==i .or. cols(3)==i .or. cols(4)==i) then
                      if (ok) ok = isreal(rdum,word)
                      if (ok) then
                         if (cols(2) == i) then
                            xat(1,nat) = rdum
                            havefields(8) = .true.
                         elseif (cols(3) == i) then
                            xat(2,nat) = rdum
                            havefields(9) = .true.
                         elseif (cols(4) == i) then
                            xat(3,nat) = rdum
                            havefields(10) = .true.
                         end if
                      end if
                   elseif (cols(1)==i) then
                      zat(nat) = zatguess(word)
                      ok = (zat(nat) /= 0)
                      havefields(7) = .true.
                   end if
                end if

                if (.not.ok) then
                   errmsg = "Error reading cif file (atom loop)"
                   return
                end if
             end do
          elseif (looptype == 1) then
             ! the symmetry operations loop
             ! add one operation
             nop = nop + 1
             if (nop > size(ops,1)) call realloc(ops,2*nop)

             ! read the fields
             lp = 1
             do i = 1, cols(5)
                ok = isexpression_or_word(word,line,lp)
                if (ok.and.cols(5)==i) ops(nop) = word
                if (.not.ok) then
                   errmsg = "Error reading cif file (symop loop)"
                   return
                end if
             end do
          end if
          cycle
       end if

       ! skip the ; blocks associated with an unknown field
       if (line(1:1) == ";") then
          do while (getline_raw(lu,line))
             line = adjustl(line)
             if (len_trim(line) == 0) cycle
             if (line(1:1) == ";") cycle main
          end do
       end if
    end do main
    call fclose(lu)

    ! fill the last seed
    call fill_seed(seed)
    if (present(nseed).and.present(mseed)) then
       nseed = nseed + 1
       call realloc_crystalseed(mseed,nseed)
       mseed(nseed) = seed
    elseif (present(seed0)) then
       seed0 = seed
    end if

  contains
    function read_field_float()
      use tools_io, only: isreal
      real*8 :: read_field_float

      read_field_float = 0d0
      if (len_trim(line(lp:)) == 0) then
         ! read the next line if empty
         ok = getline_raw(lu,line)
         if (.not.ok) then
            errmsg = "Error reading cif file (float read error)"
            return
         end if
         line = adjustl(line)
         lp = 1
      end if
      if (.not.isreal(read_field_float,line,lp)) &
         errmsg = "Error reading cif file (float read error)"

    end function read_field_float

    function read_field_string()
      use tools_io, only: isexpression_or_word
      use param, only: newline
      character(len=:), allocatable :: read_field_string

      ! read the next line if empty
      read_field_string = ""
      if (len_trim(line(lp:)) == 0) then
         ok = getline_raw(lu,line)
         if (.not.ok) then
            errmsg = "Error reading cif file (float read error)"
            return
         end if
         line = adjustl(line)
         lp = 1
      end if

      if (line(1:1) /= ";") then
         if (isexpression_or_word(read_field_string,line,lp)) return
         errmsg = "Error reading cif file (float read error)"
      else
         ! read a ; block
         do while (getline_raw(lu,line))
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == ";") then
               errmsg = ""
               return
            end if
            read_field_string = read_field_string // trim(line) // newline
         end do
      end if

    end function read_field_string

    subroutine fill_seed(seed)
      use spglib, only: spg_get_hall_number_from_symbol, spg_get_symmetry_from_database
      use types, only: realloc
      use tools_io, only: nameguess
      use param, only: bohrtoa, maxzat, eyet, eye
      type(crystalseed), intent(out) :: seed

      integer :: zspc(maxzat), hnum
      integer :: i, j
      real*8 :: rot0(3,4)
      logical :: ok, foundsym

      ! check we have all the info
      if (.not.all(havefields)) then
         errmsg = "Error reading cif file (missing fields or data block not found)"
         return
      end if

      ! initialize seed
      call seed%end()
      seed%file = trim(file)
      seed%name = trim(file) // "|" // trim(blockname)
      seed%isformat = isformat_cif
      seed%useabr = 1
      seed%aa = aa / bohrtoa
      seed%bb = bb

      ! fill the species
      zspc = 0
      do i = 1, nat
         zspc(zat(i)) = maxzat+1
      end do
      if (allocated(seed%spc)) deallocate(seed%spc)
      allocate(seed%spc(count(zspc > 0)))
      seed%nspc = 0
      do i = 1, maxzat
         if (zspc(i) > 0) then
            seed%nspc = seed%nspc + 1
            seed%spc(seed%nspc)%z = i
            seed%spc(seed%nspc)%name = nameguess(i,.true.)
            zspc(i) = seed%nspc
         end if
      end do

      ! fill the atoms
      seed%nat = nat
      if (allocated(seed%x)) deallocate(seed%x)
      if (allocated(seed%is)) deallocate(seed%is)
      allocate(seed%x(3,nat),seed%is(nat))
      do i = 1, nat
         seed%is(i) = zspc(zat(i))
         seed%x(:,i) = xat(:,i)
      end do

      ! pre-allocate the symmetry
      seed%neqv = 0
      seed%ncv = 1
      if (allocated(seed%rotm)) deallocate(seed%rotm)
      if (allocated(seed%cen)) deallocate(seed%cen)
      allocate(seed%rotm(3,4,48),seed%cen(3,4))
      seed%cen = 0d0

      ! fill the symmetry with the symops from the cif file
      foundsym = .false.
      if (nop > 0) then
         foundsym = .true.
         ! from the symmetry operations in the cif file
         do j = 1, nop
            rot0 = string_to_symop(ops(j),errmsg)
            if (len_trim(errmsg) > 0) return

            if (all(abs(eyet - rot0) < 1d-12)) then
               ! the identity
               seed%neqv = seed%neqv + 1
               if (seed%neqv > size(seed%rotm,3)) &
                  call realloc(seed%rotm,3,4,2*seed%neqv)
               seed%rotm(:,:,seed%neqv) = rot0
            elseif (all(abs(eye - rot0(1:3,1:3)) < 1d-12)) then
               ! a non-zero pure translation
               ! check if I have it already
               ok = .true.
               do i = 1, seed%ncv
                  if (all(abs(rot0(:,4) - seed%cen(:,i)) < 1d-12)) then
                     ok = .false.
                     exit
                  endif
               end do
               if (ok) then
                  seed%ncv = seed%ncv + 1
                  if (seed%ncv > size(seed%cen,2)) call realloc(seed%cen,3,2*seed%ncv)
                  seed%cen(:,seed%ncv) = rot0(:,4)
               endif
            else
               ! a rotation, with some pure translation in it
               ! check if I have this rotation matrix already
               ok = .true.
               do i = 1, seed%neqv
                  if (all(abs(seed%rotm(1:3,1:3,i) - rot0(1:3,1:3)) < 1d-12)) then
                     ok = .false.
                     exit
                  endif
               end do
               if (ok) then
                  seed%neqv = seed%neqv + 1
                  seed%rotm(:,:,seed%neqv) = rot0
               endif
            endif
         end do
      end if

      ! get the symmetry from the space group label
      if (.not.foundsym) then
         hnum = spg_get_hall_number_from_symbol(spg)
         if (hnum > 0) &
            call spg_get_symmetry_from_database(hnum,seed%neqv,seed%ncv,seed%rotm,seed%cen)
         foundsym = (seed%neqv > 0) .and. (seed%ncv > 0)
      end if

      if (.not.foundsym) then
         ! could not find symmetry in the cif file, calculate it later
         seed%neqv = 1
         seed%ncv = 1
         if (allocated(seed%rotm)) deallocate(seed%rotm)
         if (allocated(seed%cen)) deallocate(seed%cen)
         allocate(seed%rotm(3,4,1),seed%cen(3,1))
         seed%rotm(:,:,1) = eyet
         seed%cen(:,1) = 0d0
         seed%havesym = 0
         seed%checkrepeats = .false.
         seed%findsym = -1
      else
         ! we have symmetry, check for repeats
         seed%havesym = 1
         seed%checkrepeats = .true.
         seed%findsym = 0
      end if

      ! reallocate
      if (seed%neqv > 0) call realloc(seed%rotm,3,4,seed%neqv)
      if (seed%ncv > 0) call realloc(seed%cen,3,seed%ncv)

      ! rest of the seed information
      seed%isused = .true.
      seed%ismolecule = mol
      seed%cubic = .false.
      seed%border = 0d0
      seed%havex0 = .false.
      seed%molx0 = 0d0

    end subroutine fill_seed

  end subroutine read_all_cif

  !> Read multiple structure seeds from a mol2 file. If mol, force a
  !> molecule/crystal system.  If error, return non-empty errmsg.  If
  !> nseed and mseed are present, read all seeds into the mseed array
  !> and return the number of seeds in nseed. If seed0 is present,
  !> return the molecule with name name. If name is not present or is
  !> empty, return the first molecule in seed0.
  subroutine read_all_mol2(file,errmsg,nseed,mseed,seed0,name,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, getline, lower, isexpression_or_word,&
       isreal, zatguess, equal, isinteger, zatguess
    use types, only: realloc
    use param, only: maxzat, bohrtoa, isformat_mol2
    integer, intent(out), optional :: nseed
    type(crystalseed), intent(inout), allocatable, optional :: mseed(:)
    type(crystalseed), intent(inout), optional :: seed0
    character*(*), intent(in), optional :: name
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    logical :: ok, indata, havename
    character(len=:), allocatable :: line, molname
    integer :: nat, i, lp, idum
    integer :: zat
    character*10 :: attyp, atname
    type(crystalseed) :: seed
    integer, allocatable :: usedz(:)

    ! consistency check
    if (.not.(present(nseed).and.present(mseed)).and..not.present(seed0)) then
       errmsg = "Error reading mol2 file (incorrect use of the interface)"
       return
    end if

    ! initialize
    errmsg = ""
    if (present(nseed).and.present(mseed)) then
       nseed = 0
       if (allocated(mseed)) deallocate(mseed)
       allocate(mseed(10))
    end if
    havename = present(name)
    if (havename) havename = (len_trim(name) > 0)
    nat = 0
    indata = .false.

    ! main loop
    lu = fopen_read(file,ti=ti)
    main: do while (getline(lu,line))

       ! scan for ATOM (if in data-reading mode)
       if (indata) then
          if (line(1:1) == "@") then
             if (equal(line,"@<TRIPOS>ATOM")) then
                if (nat <= 0) then
                   errmsg = "Zero atoms to be read at the beginning of ATOM block"
                   goto 999
                end if

                ! fill the seed with information from this ATOM block
                call seed%end()
                seed%nat = nat
                allocate(seed%x(3,nat),seed%is(nat),usedz(maxzat))
                usedz = 0
                seed%nspc = 0
                allocate(seed%spc(10))

                do i = 1, nat
                   ok = getline(lu,line)
                   lp = 1
                   read(line,*,err=999,end=999) idum, atname, seed%x(:,i), attyp
                   zat = zatguess(attyp)
                   if (zat < 1 .or. zat > maxzat) then
                      errmsg = "error reading atom type: " // trim(attyp)
                      goto 999
                   end if
                   if (usedz(zat) == 0) then
                      seed%nspc = seed%nspc + 1
                      if (seed%nspc > size(seed%spc,1)) call realloc(seed%spc,2*seed%nspc)
                      seed%spc(seed%nspc)%name = trim(attyp)
                      seed%spc(seed%nspc)%z = zat
                      usedz(zat) = seed%nspc
                   end if
                   seed%is(i) = usedz(zat)
                end do
                call realloc(seed%spc,seed%nspc)
                deallocate(usedz)

                ! mol2 coordinates in angstrom
                seed%x = seed%x / bohrtoa

                ! rest of the seed information
                seed%useabr = 0
                seed%havesym = 0
                seed%findsym = -1
                seed%checkrepeats = .false.
                seed%isused = .true.
                seed%ismolecule = .true.
                seed%cubic = .false.
                seed%border = rborder_def
                seed%havex0 = .false.
                seed%molx0 = 0d0
                seed%file = file
                seed%isformat = isformat_mol2
                if (len_trim(molname) > 0) then
                   seed%name = trim(file) // "|" // trim(molname)
                else
                   seed%name = file
                end if
             end if
          end if
       end if

       ! scan for a MOLECULE
       if (line(1:1) == "@") then
          if (equal(line,"@<TRIPOS>MOLECULE")) then
             ! fill the last seed, if in data-reading mode
             if (indata) then
                if (present(nseed).and.present(mseed)) then
                   nseed = nseed + 1
                   call realloc_crystalseed(mseed,nseed)
                   mseed(nseed) = seed
                elseif (present(seed0)) then
                   seed0 = seed
                   exit ! read the first structure only
                end if
                indata = .false.
             end if

             ! read the MOLECULE block
             ! molecule name
             ok = getline(lu,line)
             if (.not.ok) then
                errmsg = "Error reading molecule name in MOLECULE specification"
                goto 999
             end if
             if (havename) then
                indata = equal(line,name)
             else
                indata = .true.
             end if

             if (indata) then
                ! save the molecule name
                molname = trim(adjustl(line))

                ! number of atoms
                ok = getline(lu,line)
                ok = ok .and. isinteger(nat,line)
                if (.not.ok) then
                   errmsg = "Error reading number of atoms in MOLECULE specification"
                   goto 999
                end if
             end if
          end if
       end if

    end do main
    call fclose(lu)

    ! fill the last seed, if in data-reading mode
    if (indata) then
       if (present(nseed).and.present(mseed)) then
          nseed = nseed + 1
          call realloc_crystalseed(mseed,nseed)
          mseed(nseed) = seed
       elseif (present(seed0)) then
          seed0 = seed
       end if
    end if

    ! check if a seed has been read
    if (present(nseed).and.present(mseed)) then
       if (nseed == 0) then
          errmsg = "No structures found in the mol2 file"
          goto 999
       end if
    elseif (present(seed0)) then
       if (.not.seed0%isused) then
          if (havename) then
             errmsg = "No structures with name " // trim(name) // " found in the mol2 file"
             goto 999
          else
             errmsg = "No structures found in the mol2 file"
             goto 999
          end if
       end if
    end if

    return
999 continue
    call fclose(lu)
    if (present(nseed)) nseed = 0
    if (present(mseed)) deallocate(mseed)

  end subroutine read_all_mol2

  !> Read one or all structures from a QE output (filename file) and
  !> return the corresponding crystal seeds in seed. If istruct < 0,
  !> read all seeds and return the number of seeds read in nseed. If
  !> istruct = 0, return a single seed for the last structure. If
  !> istruct > 0, return that particular structure. If mol=.true.,
  !> interpret the structure as a molecule (currently, this only sets
  !> the %ismolecule field). If an error condition is found, return
  !> the error message in errmsg (zero-length string if no error).
  subroutine read_all_qeout(nseed,seed,file,mol,istruct,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, isinteger, isreal,&
       zatguess, fclose, equali, string
    use tools_math, only: matinv
    use param, only: bohrtoa, isformat_qeout
    use types, only: realloc, species
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: istruct !< ID of the structure
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ideq, i, j, is0, ier
    character(len=:), allocatable :: line, str
    character*10 :: atn, sdum
    character*40 :: sene
    integer :: idum, npad, idx
    real*8 :: alat, r(3,3), qaux, rfac, cfac, rdum
    logical :: ok, tox
    ! interim copy of seed info
    integer :: nat, nspc, iuse
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    type(species), allocatable :: spc(:) !< Species
    real*8 :: m_x2c(3,3)
    logical :: hasx, hasis, hasspc, hasr

    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! first pass: read the number of structures
    nseed = 0
    do while (getline_raw(lu,line))
       if (index(line,"!") == 1) then
          nseed = nseed + 1
       end if
    end do
    if (nseed == 0) then
       errmsg = "No valid structures found."
       goto 999
    end if
    if (allocated(seed)) deallocate(seed)

    is0 = 0
    if (istruct >= 0) then
       allocate(seed(1))
       seed(1)%nspc = 0
       seed(1)%nat = 0
    else
       allocate(seed(nseed))
       do i = 1, nseed
          seed(i)%nspc = 0
          seed(i)%nat = 0
       end do
    end if
    alat = 1d0
    npad = ceiling(log10(nseed-1+0.1d0))

    ! rewind and read all the structures
    rewind(lu)
    errmsg = "Error reading file: " // trim(file)
    nat = 0
    nspc = 0
    tox = .false.
    hasx = .false.
    hasis = .false.
    hasspc = .false.
    hasr = .false.
    do while (getline_raw(lu,line))
       ideq = index(line,"=") + 1

       ! Count the structures
       if (index(line,"lattice parameter (alat)") > 0) then
          ok = isreal(alat,line,ideq)

       elseif (index(line,"number of atoms/cell") > 0) then
          ok = isinteger(nat,line,ideq)
          if (allocated(x)) deallocate(x)
          if (allocated(is)) deallocate(is)
          allocate(x(3,nat),is(nat))

       elseif (index(line,"number of atomic types") > 0) then
          ok = isinteger(nspc,line,ideq)
          if (allocated(spc)) deallocate(spc)
          allocate(spc(nspc))

       elseif (index(line,"atomic species   valence    mass     pseudopotential")>0) then
          do i = 1, nspc
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read (line,*,err=999,end=999) spc(i)%name, qaux
             spc(i)%z = zatguess(spc(i)%name)
             if (spc(i)%z < 0) then
                errmsg = "Unknown atomic symbol: "//trim(spc(i)%name)//"."
                goto 999
             end if
          end do
          hasspc = .true.

       elseif (index(line,"crystal axes:") > 0) then
          do i = 1, 3
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             ideq = index(line,"(",.true.) + 1
             ok = isreal(r(i,1),line,ideq)
             ok = ok.and.isreal(r(i,2),line,ideq)
             ok = ok.and.isreal(r(i,3),line,ideq)
             if (.not.ok) goto 999
          end do
          r = r * alat ! alat comes before crystal axes
          m_x2c = transpose(r)
          tox = .false.
          hasr = .true.

       elseif (index(line,"Cartesian axes")>0) then
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
          is = 0
          do i = 1, nat
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read(line,*,err=999,end=999) idum, atn
             line = line(index(line,"(",.true.)+1:)
             read(line,*,err=999,end=999) x(:,i)
             do j = 1, nspc
                if (equali(spc(j)%name,atn)) then
                   is(i) = j
                   exit
                end if
             end do
             if (is(i) == 0) then
                errmsg = "Unknown atom type: "//atn
                goto 999
             end if
          end do
          ! this is Cartesian in alat units
          x = x * alat
          tox = .true.
          hasx = .true.
          hasis = .true.

       elseif (line(1:15) == "CELL_PARAMETERS") then
          cfac = 1d0
          if (index(line,"angstrom") > 0) then
             cfac = 1d0 / bohrtoa
          elseif (index(line,"alat") > 0) then
             cfac = alat
          elseif (index(line,"bohr") > 0) then
             cfac = 1d0
          end if
          do i = 1, 3
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             ideq = 1
             ok = isreal(r(i,1),line,ideq)
             ok = ok.and.isreal(r(i,2),line,ideq)
             ok = ok.and.isreal(r(i,3),line,ideq)
             if (.not.ok) goto 999
          end do
          r = r * cfac
          m_x2c = transpose(r)
          hasr = .true.

       elseif (line(1:16) == "ATOMIC_POSITIONS") then
          rfac = 1d0
          if (index(line,"angstrom") > 0) then
             tox = .true.
             rfac = 1d0 / bohrtoa
          elseif (index(line,"alat") > 0) then
             tox = .true.
             rfac = alat
          elseif (index(line,"bohr") > 0) then
             tox = .true.
             rfac = 1d0
          elseif (index(line,"crystal") > 0) then
             tox = .false.
             rfac = 1d0
          end if
          is = 0
          do i = 1, nat
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read(line,*,err=999,end=999) atn, x(:,i)
             do j = 1, nspc
                if (equali(spc(j)%name,atn)) then
                   is(i) = j
                   exit
                end if
             end do
             if (is(i) == 0) then
                errmsg = "Unknown atom type: "//atn
                goto 999
             end if
          end do
          x = x * rfac
          hasx = .true.
       else if (index(line,"!") == 1) then
          if (.not.hasx .or. nat == 0) then
             errmsg = "Missing atomic positions."
             goto 999
          end if
          if (.not.hasis) then
             errmsg = "Missing atomic types."
             goto 999
          end if
          if (.not.hasspc .or. nspc == 0) then
             errmsg = "Missing atomic species."
             goto 999
          end if
          if (.not.hasr) then
             errmsg = "Missing cell dimensions."
             goto 999
          end if
          is0 = is0 + 1
          hasx = .false.

          ! decide whether we want to keep this structure in a seed
          iuse = 0
          if (istruct < 0) then
             iuse = is0
          elseif (istruct == 0 .and. is0 == nseed) then
             iuse = 1
          elseif (istruct == is0) then
             iuse = 1
          end if

          ! keep the seed
          if (iuse > 0) then
             seed(iuse)%nat = nat
             seed(iuse)%nspc = nspc
             seed(iuse)%spc = spc
             seed(iuse)%x = x
             seed(iuse)%is = is
             seed(iuse)%m_x2c = m_x2c

             seed(iuse)%useabr = 2
             r = seed(iuse)%m_x2c
             call matinv(r,3,ier)
             if (ier /= 0) then
                errmsg = "Error inverting matrix"
                goto 999
             end if

             do i = 1, seed(iuse)%nat
                if (tox) then
                   seed(iuse)%x(:,i) = matmul(r,seed(iuse)%x(:,i))
                end if
                seed(iuse)%x(:,i) = seed(iuse)%x(:,i) - floor(seed(iuse)%x(:,i))
             end do

             seed(iuse)%havesym = 0
             seed(iuse)%checkrepeats = .false.
             seed(iuse)%findsym = -1
             seed(iuse)%isused = .true.
             seed(iuse)%ismolecule = mol
             seed(iuse)%cubic = .false.
             seed(iuse)%border = 0d0
             seed(iuse)%havex0 = .false.
             seed(iuse)%molx0 = 0d0
             seed(iuse)%file = file
             seed(iuse)%isformat = isformat_qeout

             read (line,*,err=999,end=999) sdum, sdum, sdum, sdum, sene
             read (sene,*,err=999,end=999) rdum
             if (istruct < 0) then
                if (is0 == nseed) then
                   seed(iuse)%name = trim(file) // "|(fin) (" //&
                      trim(adjustl(string(rdum,'f',20,8))) // " Ry)"
                   seed(iuse)%energy = rdum / 2d0
                else
                   str = string(iuse,npad,pad0=.true.)
                   str = string(str,length=max(5,len(str)))
                   seed(iuse)%name = trim(file) // "|" // str // " (" //&
                      trim(adjustl(string(rdum,'f',decimal=8))) // " Ry)"
                   seed(iuse)%energy = rdum / 2d0
                end if
             else
                seed(iuse)%name = file
             end if
          end if
       else if (iuse > 0 .and. index(line,"total   stress") > 0) then
          ! add the pressure to the last seed, if available
          idx = index(line,'=')
          ok = isreal(seed(iuse)%pressure,line(idx+1:))
          if (ok) then
             seed(iuse)%pressure = seed(iuse)%pressure / 10d0 ! kbar -> GPa
          else
             seed(iuse)%pressure = huge(1d0)
          end if
       end if
    end do

    if (istruct >= 0) nseed = 1

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  end subroutine read_all_qeout

  !> Read all structures from an xyz file. Returns all seeds.
  subroutine read_all_xyz(nseed,seed,file,errmsg,ti)
    use global, only: rborder_def
    use hashmod, only: hash
    use tools_io, only: fopen_read, fclose, getline_raw, lower, zatguess,&
       isinteger, isreal, string, nameguess
    use types, only: realloc
    use param, only: maxzat, bohrtoa, isformat_xyz
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, nat, i, iz, lp
    logical :: ok
    real*8 :: edum
    real*8, allocatable :: energy(:)
    character(len=:), allocatable :: line, latn
    character*10 :: atn
    type(hash) :: usen

    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    errmsg = "Error reading file: " // trim(file)
    nseed = 0
    allocate(energy(10))
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) cycle
       read (line,*,err=999,end=999) nat
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999

       ! interpret the title line, if possible
       lp = 1
       ok = isreal(edum,line,lp)
       if (.not.ok .or. len_trim(line(lp:)) > 0) then ! a single numerical field -> energy
          edum = 0d0
       end if

       do i = 1, nat
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
       end do
       nseed = nseed + 1
       if (nseed > size(energy,1)) call realloc(energy,2*nseed)
       energy(nseed) = edum
    end do

    if (allocated(seed)) deallocate (seed)
    allocate(seed(nseed))
    rewind(lu)
    nseed = 0
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) cycle
       nseed = nseed + 1
       call usen%init()
       read (line,*,err=999,end=999) nat
       seed(nseed)%nat = nat

       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       seed(nseed)%file = file
       seed(nseed)%name = seed(nseed)%file
       seed(nseed)%isformat = isformat_xyz

       seed(nseed)%nspc = 0
       allocate(seed(nseed)%x(3,nat),seed(nseed)%is(nat),seed(nseed)%spc(10))
       do i = 1, nat
          read (lu,*,err=999,end=999) atn, seed(nseed)%x(:,i)

          ok = isinteger(iz,atn)
          if (ok) then
             if (iz < 0 .or. iz > maxzat) then
                errmsg = "Invalid atomic number: "//string(iz)//"."
                goto 999
             end if
             atn = nameguess(iz,.true.)
          else
             iz = zatguess(atn)
             if (iz < 0) then
                errmsg = "Unknown atomic symbol: "//trim(atn)//"."
                goto 999
             end if
          end if

          latn = lower(trim(atn))
          if (usen%iskey(latn)) then
             seed(nseed)%is(i) = usen%get(latn,1)
          else
             seed(nseed)%nspc = seed(nseed)%nspc + 1
             if (seed(nseed)%nspc > size(seed(nseed)%spc,1)) &
                call realloc(seed(nseed)%spc,2*seed(nseed)%nspc)
             seed(nseed)%spc(seed(nseed)%nspc)%name = trim(atn)
             seed(nseed)%spc(seed(nseed)%nspc)%z = iz
             seed(nseed)%spc(seed(nseed)%nspc)%qat = 0d0
             call usen%put(latn,seed(nseed)%nspc)
             seed(nseed)%is(i) = seed(nseed)%nspc
          end if
       end do
       call realloc(seed(nseed)%spc,seed(nseed)%nspc)
       seed(nseed)%x = seed(nseed)%x / bohrtoa
       seed(nseed)%useabr = 0
       seed(nseed)%havesym = 0
       seed(nseed)%checkrepeats = .false.
       seed(nseed)%findsym = -1
       seed(nseed)%isused = .true.
       seed(nseed)%ismolecule = .true.
       seed(nseed)%cubic = .false.
       seed(nseed)%border = rborder_def
       seed(nseed)%havex0 = .false.
       seed(nseed)%molx0 = 0d0
       seed(nseed)%energy = energy(nseed)
    end do
    deallocate(energy)

    if (nseed > 1) then
       do i = 1, nseed
          seed(i)%name = trim(seed(i)%name) // "|" // string(i)
       end do
    end if

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  end subroutine read_all_xyz

  !> Read all structures from a Gaussian output (log) file. Returns
  !> all crystal seeds.
  subroutine read_all_log(nseed,seed,file,errmsg,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, getline_raw, nameguess, string
    use types, only: species
    use param, only: maxzat, bohrtoa, isformat_gaussian
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, str
    integer :: lu, nat, idum, iz, nspc, i, npad
    integer :: usez(0:maxzat), idx, in
    logical :: ok, laste, lastinputor
    type(species), allocatable :: spc(:)
    real*8 :: energy
    real*8, allocatable :: esave(:)

    errmsg = ""

    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! count the number of seeds, atoms, and build the species
    energy = huge(1d0)
    nat = 0
    nseed = 0
    lastinputor = .false.
    do while (getline_raw(lu,line))
       ok = (index(line,"Input orientation:") > 0)
       if (ok) then
          lastinputor = .true.
       elseif (.not.lastinputor) then
          ok = (index(line,"Standard orientation:") > 0)
          if (ok) lastinputor = .false.
       end if

       if (ok) then
          nseed = nseed + 1

          if (nat == 0) then
             usez = 0
             ok = getline_raw(lu,line)
             ok = ok .and. getline_raw(lu,line)
             ok = ok .and. getline_raw(lu,line)
             ok = ok .and. getline_raw(lu,line)
             if (.not.ok) goto 999
             do while (.true.)
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999
                if (index(line,"---------") > 0) exit
                nat = nat + 1
                read(line,*,err=999,end=999) idum, iz
                usez(iz) = 1
             end do
          end if
       end if
    end do
    if (nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if

    ! build the species
    nspc = count(usez > 0)
    if (nspc == 0) then
       errmsg = "No species found."
       goto 999
    end if
    allocate(spc(nspc))
    nspc = 0
    do i = 0, maxzat
       if (usez(i) > 0) then
          nspc = nspc + 1
          spc(nspc)%z = i
          spc(nspc)%name = nameguess(i,.true.)
          spc(nspc)%qat = 0d0
          usez(i) = nspc
       end if
    end do

    if (allocated(seed)) deallocate(seed)
    allocate(seed(nseed),esave(nseed))
    esave = huge(1d0)
    rewind(lu)
    in = 0
    lastinputor = .false.
    do while (getline_raw(lu,line))
       ok = (index(line,"Input orientation:") > 0)
       if (ok) then
          lastinputor = .true.
       elseif (.not.lastinputor) then
          ok = (index(line,"Standard orientation:") > 0)
          if (ok) lastinputor = .false.
       end if

       if (ok) then
          in = in + 1
          seed(in)%nat = nat
          allocate(seed(in)%x(3,nat),seed(in)%is(nat))

          ok = getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          if (.not.ok) goto 999

          do i = 1, nat
             read (lu,*,err=999,end=999) idum, iz, idum, seed(in)%x(:,i)
             seed(in)%is(i) = usez(iz)
          end do

          seed(in)%x = seed(in)%x / bohrtoa
          seed(in)%isused = .true.
          seed(in)%file = file
          seed(in)%isformat = isformat_gaussian
          seed(in)%name = file
          seed(in)%nspc = nspc
          seed(in)%spc = spc
          seed(in)%useabr = 0
          seed(in)%havesym = 0
          seed(in)%checkrepeats = .false.
          seed(in)%findsym = -1
          seed(in)%isused = .true.
          seed(in)%ismolecule = .true.
          seed(in)%cubic = .false.
          seed(in)%border = rborder_def
          seed(in)%havex0 = .false.
          seed(in)%molx0 = 0d0
          laste = .false.
       elseif (index(line,"SCF Done") > 0) then
          idx = index(line,"=")
          if (idx > 0) then
             line = line(idx+1:)
             read (line,*) energy
             esave(in) = energy
          end if
          laste = .true.
       end if
    end do
    if (.not.laste.and.nseed > 1) then
       nseed = nseed - 1
       call realloc_crystalseed(seed,nseed)
    end if

    if (nseed > 1) then
       npad = ceiling(log10(nseed-1+0.1d0))
       do i = 1, nseed-1
          str = string(i,npad,pad0=.true.)
          str = string(str,length=max(5,len(str)))
          seed(i)%name = trim(file) // "|" // str // " (" //&
             trim(adjustl(string(esave(i),'f',decimal=9))) // " Ha)"
          seed(i)%energy = esave(i)
       end do
       seed(nseed)%name = trim(file) // "|(fin) (" //&
          trim(adjustl(string(energy,'f',decimal=9))) // " Ha)"
       seed(nseed)%energy = energy
    else
       seed(in)%name = trim(file)
       if (seed(in)%energy /= huge(1d0)) seed(in)%energy = energy
    end if

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  end subroutine read_all_log

  !> Read an FHI aims output file.
  subroutine read_all_aimsout(nseed,seed,file,errmsg,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isreal, &
       getword, zatguess, string
    use tools_math, only: matinv
    use types, only: realloc
    use hashmod, only: hash
    use param, only: bohrtoa, hartoev, eva3togpa, isformat_aimsout
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nlat, i, j, idx, ier, iat, npad
    logical :: is_file_mol, ok, isfinal
    character*1 :: cdum
    character*10 :: dum1, dum2, splbl
    character(len=:), allocatable :: line, word, str
    real*8 :: rlat(3,3)
    logical, allocatable :: isfrac(:)
    type(hash) :: usespc

    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)
    call usespc%init()

    ! allocate initial seed
    nseed = 1
    if (allocated(seed)) deallocate(seed)
    allocate(seed(10))

    ! allocate atoms
    allocate(isfrac(10),seed(1)%x(3,10),seed(1)%is(10))

    ! advance until we are at the "Input geometry" line
    ok = .false.
    do while (getline_raw(lu,line))
       if (line == "  Input geometry:") then
          ok = .true.
          exit
       end if
    end do
    if (.not.ok) then
       errmsg = "Failed to locate the Input geometry block in the FHIaims output file"
       goto 999
    end if

    ! get the cell geometry
    ok = getline_raw(lu,line)
    if (.not.ok) goto 999
    if (index(line,"Unit cell:") > 0) then
       do i = 1, 3
          ok = getline_raw(lu,line)
          if (.not.ok) goto 999
          read (line,*,err=999,end=999) cdum, (rlat(j,i),j=1,3)
       end do
       is_file_mol = .false.
    elseif (index(line,"No unit cell requested.") > 0) then
       is_file_mol = .true.
    else
       goto 999
    end if

    ! get the atomic positions
    seed(1)%nspc = 0
    seed(1)%nat = 0
    ok = getline_raw(lu,line)
    ok = ok .and. getline_raw(lu,line)
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) exit

       seed(1)%nat = seed(1)%nat + 1
       if (seed(1)%nat > size(isfrac,1)) then
          call realloc(isfrac,2*seed(1)%nat)
          call realloc(seed(1)%x,3,2*seed(1)%nat)
          call realloc(seed(1)%is,2*seed(1)%nat)
       end if
       read(line,*,err=999,end=999) cdum, dum1, dum2, splbl, (seed(1)%x(j,seed(1)%nat),j=1,3)
       isfrac(seed(1)%nat) = .false.

       word = trim(adjustl(splbl))
       if (.not.usespc%iskey(word)) then
          seed(1)%nspc = seed(1)%nspc + 1
          call usespc%put(word,seed(1)%nspc)
          seed(1)%is(seed(1)%nat) = seed(1)%nspc
       else
          seed(1)%is(seed(1)%nat) = usespc%get(word,1)
       end if
    end do
    call realloc(isfrac,seed(1)%nat)
    call realloc(seed(1)%x,3,seed(1)%nat)
    call realloc(seed(1)%is,seed(1)%nat)

    ! fill the species array
    allocate(seed(1)%spc(seed(1)%nspc))
    do i = 1, seed(1)%nspc
       word = usespc%getkey(i)
       idx = usespc%get(word,1)
       seed(1)%spc(idx)%name = word
       seed(1)%spc(idx)%z = zatguess(word)
       if (seed(1)%spc(idx)%z <= 0) then
          errmsg = "unknown atom type: " // word
          goto 999
       end if
    end do

    ! ismolecule and cell
    seed(1)%ismolecule = is_file_mol
    if (is_file_mol) then
       seed(1)%useabr = 0
       seed(1)%m_x2c = 0d0
    else
       seed(1)%useabr = 2
       rlat = rlat / bohrtoa
       seed(1)%m_x2c = rlat
       call matinv(rlat,3,ier)
       if (ier /= 0) then
          errmsg = "Error inverting lattice vector matrix"
          goto 999
       end if
    end if

    ! convert the atomic coordinates
    do i = 1, seed(1)%nat
       if (.not.isfrac(i)) then
          seed(1)%x(:,i) = seed(1)%x(:,i) / bohrtoa
          if (.not.is_file_mol) seed(1)%x(:,i) = matmul(rlat,seed(1)%x(:,i))
       end if
    end do

    ! read the rest of the file
    main: do while (.true.)
       ! search for "Updated atomic structure" or "Final atomic structure" blocks
       ! read energy and pressure as we go
       do while (.true.)
          ok = getline_raw(lu,line)
          if (.not.ok) exit main
          if (index(line,'| Total energy uncorrected') > 0) then
             idx = index(line,':')
             ok = isreal(seed(nseed)%energy,line(idx+1:))
             if (.not.ok) seed(nseed)%energy = huge(1d0)
          elseif (index(line,'|  Pressure') > 0) then
             idx = index(line,':')
             ok = isreal(seed(nseed)%pressure,line(idx+1:))
             if (ok) then
                seed(nseed)%pressure = seed(nseed)%pressure * eva3togpa
             else
                seed(nseed)%pressure = huge(1d0)
             end if
          elseif (trim(line) == "  Updated atomic structure:") then
             nseed = nseed + 1
             isfinal = .false.
             exit
          elseif (trim(line) == "  Final atomic structure:") then
             nseed = nseed + 1
             isfinal = .true.
             exit
          end if
       end do

       ! We are about to read a new seed, make space for it and initialize
       if (nseed > size(seed,1)) call realloc_crystalseed(seed,2*nseed)
       seed(nseed)%nspc = seed(1)%nspc
       seed(nseed)%spc = seed(1)%spc
       seed(nseed)%nat = seed(1)%nat
       seed(nseed)%is = seed(1)%is
       allocate(seed(nseed)%x(3,seed(1)%nat))

       ! read the geometry block
       nlat = 0
       iat = 0
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       do while (getline_raw(lu,line))
          if (line == "  Fractional coordinates:" .or. line(1:4) == "----") exit
          if (len_trim(line) == 0) cycle
          lp = 1
          word = lgetword(line,lp)
          if (word == "lattice_vector") then
             nlat = nlat + 1
             ok = isreal(rlat(1,nlat),line,lp)
             ok = ok .and. isreal(rlat(2,nlat),line,lp)
             ok = ok .and. isreal(rlat(3,nlat),line,lp)
             if (.not.ok) goto 999
          elseif (word == "atom") then
             iat = iat + 1
             ok = isreal(seed(nseed)%x(1,iat),line,lp)
             ok = ok .and. isreal(seed(nseed)%x(2,iat),line,lp)
             ok = ok .and. isreal(seed(nseed)%x(3,iat),line,lp)
             isfrac(iat) = .false.
          end if
       end do

       ! handle the flags that are the same as the first seed
       seed(nseed)%ismolecule = seed(1)%ismolecule
       seed(nseed)%useabr = seed(1)%useabr

       ! lattice vectors
       if (is_file_mol) then
          seed(nseed)%m_x2c = 0d0
       else
          rlat = rlat / bohrtoa
          seed(nseed)%m_x2c = rlat
          call matinv(rlat,3,ier)
          if (ier /= 0) then
             errmsg = "Error inverting lattice vector matrix"
             goto 999
          end if
       end if

       ! convert the atomic coordinates
       do i = 1, seed(nseed)%nat
          if (.not.isfrac(i)) then
             seed(nseed)%x(:,i) = seed(nseed)%x(:,i) / bohrtoa
             if (.not.is_file_mol) seed(nseed)%x(:,i) = matmul(rlat,seed(nseed)%x(:,i))
          end if
       end do

       ! if this is the final structure and we do not have an energy, take it from last step
       if (isfinal .and. seed(nseed)%energy == huge(1d0) .and. nseed > 1)&
          seed(nseed)%energy = seed(nseed-1)%energy
    end do main
    call fclose(lu)
    call realloc_crystalseed(seed,nseed)

    npad = ceiling(log10(nseed-1+0.1d0))
    do i = 1, nseed
       ! symmetry
       seed(i)%havesym = 0
       seed(i)%findsym = -1
       seed(i)%checkrepeats = .false.

       ! rest of the seed information
       seed(i)%isused = .true.
       seed(i)%cubic = .false.
       seed(i)%border = rborder_def
       seed(i)%havex0 = .false.
       seed(i)%molx0 = 0d0
       seed(i)%file = file
       seed(i)%name = file
       seed(i)%isformat = isformat_aimsout

       ! name and energy conversion
       if (i == nseed) then
          if (seed(i)%energy /= huge(1d0)) then
             seed(i)%name = trim(file) // "|(fin) (" //&
                trim(adjustl(string(seed(i)%energy,'f',decimal=8))) // " eV)"
          else
             seed(i)%name = trim(file) // "|(fin)"
          end if
       else
          str = string(i,npad,pad0=.true.)
          str = string(str,length=max(5,len(str)))
          if (seed(i)%energy /= huge(1d0)) then
             seed(i)%name = trim(file) // "|" // str // " (" //&
                trim(adjustl(string(seed(i)%energy,'f',decimal=8))) // " eV)"
          else
             seed(i)%name = trim(file) // "|" // str
          end if
       end if
       if (seed(i)%energy /= huge(1d0)) &
          seed(i)%energy = seed(i)%energy / hartoev
    end do

    ! fin
    errmsg = ""
    return
999 continue ! error condition
    call fclose(lu)
    nseed = 0
    deallocate(seed)

  end subroutine read_all_aimsout

  !> Read all seeds from a CASTEP geom file.
  subroutine read_all_castep_geom(nseed,seed,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword,&
       getword, lower, isinteger, isreal, zatguess, string
    use tools_math, only: matinv
    use hashmod, only: hash
    use param, only: hartoev, isformat_castepgeom
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, lword, str
    integer :: lu, ll, i, j, lp, idum, is, ier
    integer :: nat, nspc, npad
    logical :: ok
    type(hash) :: usen
    logical, allocatable :: usespc(:)
    real*8 :: m(3,3)

    call usen%init()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! first pass, read the number of lines and the number of atoms
    nseed = 0
    nat = 0
    nspc = 0
    do while(getline_raw(lu,line))
       line = trim(adjustl(line))
       ll = len(line)
       if (line(ll-4:ll) == "<-- E") then
          nseed = nseed + 1
       elseif (nseed == 1) then
          if (line(ll-4:ll) == "<-- R") then
             nat = nat + 1
             lp = 1
             word = lgetword(line,lp)
             if (.not.usen%iskey(word)) then
                nspc = nspc + 1
                call usen%put(word,nspc)
             end if
          end if
       end if
    end do
    if (nseed == 0) goto 999

    ! allocate seeds
    if (allocated(seed)) deallocate(seed)
    allocate(seed(nseed))
    do i = 1, nseed
       allocate(seed(i)%x(3,nat),seed(i)%is(nat),seed(i)%spc(nspc))
       seed(i)%useabr = 2
       seed(i)%nat = nat
       seed(i)%nspc = nspc
    end do

    ! second pass, actual read
    allocate(usespc(nspc))
    usespc = .false.
    rewind(lu)
    nseed = 0
    do while(getline_raw(lu,line))
       line = trim(adjustl(line))
       ll = len(line)
       if (line(ll-4:ll) == "<-- E") then
          nseed = nseed + 1
          nat = 0
          read(line,*,err=999,end=999) seed(nseed)%energy
       elseif (line(ll-4:ll) == "<-- h") then
          read(line,*,err=999,end=999) seed(nseed)%m_x2c(:,1)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) seed(nseed)%m_x2c(:,2)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) seed(nseed)%m_x2c(:,3)
       elseif (line(ll-4:ll) == "<-- R") then
          nat = nat + 1
          lp = 1
          word = getword(line,lp)
          lword = lower(word)
          ok = isinteger(idum,line,lp)
          ok = ok .and. isreal(seed(nseed)%x(1,nat),line,lp)
          ok = ok .and. isreal(seed(nseed)%x(2,nat),line,lp)
          ok = ok .and. isreal(seed(nseed)%x(3,nat),line,lp)
          if (.not.ok) goto 999

          is = usen%get(lword,1)
          seed(nseed)%is(nat) = is
          if (.not.usespc(is)) then
             seed(1)%spc(is)%name = trim(word)
             seed(1)%spc(is)%qat = 0d0
             if (isinteger(idum,word)) then
                seed(1)%spc(is)%z = idum
             else
                seed(1)%spc(is)%z = zatguess(word)
                if (seed(1)%spc(is)%z < 0) then
                   errmsg = "Unknown atomic symbol: " // word
                   goto 999
                end if
             end if
             usespc(is) = .true.
          end if
       end if
    end do

    ! copy the species
    do i = 2, nseed
       seed(i)%spc = seed(1)%spc
    end do

    ! transform to fractional coordinates
    do i = 1, nseed
       m = seed(i)%m_x2c
       call matinv(m,3,ier)
       if (ier /= 0) then
          errmsg = "error inverting lattice vector matrix"
          goto 999
       end if
       do j = 1, seed(i)%nat
          seed(i)%x(:,j) = matmul(m,seed(i)%x(:,j))
       end do
    end do

    ! rest of the seed info
    npad = ceiling(log10(nseed-1+0.1d0))
    do i = 1, nseed
       ! no symmetry
       seed(i)%havesym = 0
       seed(i)%findsym = -1
       seed(i)%checkrepeats = .false.

       ! rest of the seed information
       seed(i)%isused = .true.
       seed(i)%ismolecule = .false.
       seed(i)%cubic = .false.
       seed(i)%border = 0d0
       seed(i)%havex0 = .false.
       seed(i)%molx0 = 0d0
       seed(i)%file = file
       seed(i)%isformat = isformat_castepgeom
       if (i == nseed) then
          seed(i)%name = trim(file) // "|(fin) (" //&
             trim(adjustl(string(seed(i)%energy*hartoev,'f',20,8))) // " eV)"
       else
          str = string(i,npad,pad0=.true.)
          str = string(str,length=max(5,len(str)))
          seed(i)%name = trim(file) // "|" // str // " (" //&
             trim(adjustl(string(seed(i)%energy*hartoev,'f',20,8))) // " eV)"
       end if
    end do

    errmsg = ""
999 continue
    call fclose(lu)

  end subroutine read_all_castep_geom

  !> Read a molecular geometry from a pdb file. Returns the number of
  !> atoms (n), atomic positions (x, in bohr), atomic numbers (z), and
  !> atomic names (name). If error, return a non-empty errmsg.
  subroutine read_pdb_geometry(file,n,x,z,name,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, zatguess
    use types, only: realloc
    use param, only: bohrtoa
    character*(*), intent(in) :: file
    integer, intent(out) :: n
    real*8, allocatable, intent(inout) :: x(:,:)
    integer, allocatable, intent(inout) :: z(:)
    character*(10), allocatable, intent(inout) :: name(:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    character(len=:), allocatable :: line

    errmsg = ""
    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    n = 0
    allocate(x(3,10),z(10),name(10))
    main: do while (getline_raw(lu,line))
       if (len(line) > 6) then
          if (line(1:4) == "ATOM" .or. line(1:6) == "HETATM") then
             n = n + 1
             if (n > size(z,1)) then
                call realloc(x,3,2*n)
                call realloc(z,2*n)
                call realloc(name,2*n)
             end if

             z(n) = zatguess(line(77:78))
             if (z(n) <= 0) goto 999
             name(n) = line(13:16)
             read (line(31:38),*,err=999,end=999) x(1,n)
             read (line(39:46),*,err=999,end=999) x(2,n)
             read (line(47:54),*,err=999,end=999) x(3,n)
          end if
       end if
    end do main

    if (n == 0) then
       errmsg = "No atoms found."
       goto 999
    else
       call realloc(x,3,n)
       call realloc(z,n)
       call realloc(name,n)
       x = x / bohrtoa
    endif

    errmsg = ""
999 continue
    call fclose(lu)
  end subroutine read_pdb_geometry

  !> Read a molecular geometry from a zmat file. Returns the number of
  !> atoms (n), atomic positions (x, in bohr), atomic numbers (z), and
  !> atomic names (name). If error, return a non-empty errmsg.
  subroutine read_zmat_geometry(file,n,x,z,name,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, zatguess,&
       isinteger, getword, isreal
    use types, only: realloc
    use param, only: bohrtoa
    character*(*), intent(in) :: file
    integer, intent(out) :: n
    real*8, allocatable, intent(inout) :: x(:,:)
    integer, allocatable, intent(inout) :: z(:)
    character*(10), allocatable, intent(inout) :: name(:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, idum
    character(len=:), allocatable :: line, word
    logical :: first, ok
    real*8 :: dist, ang, dih, xaux(3)
    integer :: iat(3)

    errmsg = ""
    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    first = .true.
    n = 0
    allocate(x(3,10),z(10),name(10))
    main: do while (getline_raw(lu,line))
       ! skip the first line if it contains two integers (q and mult)
       if (first) then
          lp = 1
          ok = isinteger(idum,line,lp)
          ok = ok .and. isinteger(idum,line,lp)
          if (ok) cycle
       end if
       if (len_trim(line) == 0) cycle

       ! new atom
       n = n + 1
       if (n > size(z,1)) then
          call realloc(x,3,2*n)
          call realloc(z,2*n)
          call realloc(name,2*n)
       end if

       ! read symbol
       lp = 1
       word = getword(line,lp)
       name(n) = word
       z(n) = zatguess(word)
       x(:,n) = 0d0
       if (z(n) <= 0) goto 999

       ! first atom: stop here
       if (n == 1) cycle

       ! read distance
       ok = isinteger(iat(1),line,lp)
       ok = ok .and. isreal(dist,line,lp)
       if (.not.ok) goto 999
       if (iat(1) < 1 .or. iat(1) >= n) goto 999

       ! second atom: stop here
       if (n == 2) then
          x(3,n) = dist
          cycle
       end if

       ! read angle
       ok = isinteger(iat(2),line,lp)
       ok = ok .and. isreal(ang,line,lp)
       if (.not.ok) goto 999
       if (iat(2) < 1 .or. iat(2) >= n) goto 999

       ! third atom: stop here
       if (n == 3) then
          xaux = (/1d0,0d0,0d0/)
          dih = 0d0
          x(:,n) = zmat_step(x(:,iat(1)),x(:,iat(2)),xaux,dist,ang,dih)
          cycle
       end if

       ! read dihedral
       ok = isinteger(iat(3),line,lp)
       ok = ok .and. isreal(dih,line,lp)
       if (.not.ok) goto 999
       if (iat(3) < 1 .or. iat(3) >= n) goto 999

       ! calculate the position
       x(:,n) = zmat_step(x(:,iat(1)),x(:,iat(2)),x(:,iat(3)),dist,ang,dih)
    end do main

    if (n == 0) then
       errmsg = "No atoms found."
       goto 999
    else
       call realloc(x,3,n)
       call realloc(z,n)
       call realloc(name,n)
       x = x / bohrtoa
    endif

    errmsg = ""
999 continue
    call fclose(lu)
  contains
    ! Calculate the coordinates of atom number 4 given the coordinates
    ! of atoms 1, 2 and 3, and the distance 4-1, angle 4-1-2 and dihedral
    ! 4-1-2-3.
    function zmat_step(x0_,x1_,x2_,d_,ang_,dieh_)
      use tools_math, only: cross, matinv, tosphere
      use param, only: pi
      real*8, intent(in) :: x0_(3), x1_(3), x2_(3), d_, ang_, dieh_
      real*8 :: zmat_step(3)

      real*8 :: x0(3), x1(3), x2(3), x2p(3), d, ang, dieh
      real*8 :: crot(3,3), xaux(3)
      real*8 :: asph(2), r2, rf, phf, thf, xf(3)

      ! copy the variables for work space
      x0 = x0_
      x1 = x1_
      x2 = x2_
      d = d_
      ang = ang_
      dieh = dieh_

      ! convert to radians
      ang = ang * pi / 180d0
      dieh = dieh * pi / 180d0

      ! subtract the origin
      x1 = x1 - x0
      x2 = x2 - x0

      ! x1-x0 is aligned to z
      crot(:,3) = x1 / norm2(x1)
      xaux = (/0d0,0d0,1d0/)
      crot(:,1) = cross(x1,xaux)
      if (abs(norm2(crot(:,1))) < 1d-12) then
         xaux = (/0d0,1d0,0d0/)
         crot(:,1) = cross(x1,xaux)
      end if
      crot(:,1) = crot(:,1) / norm2(crot(:,1))
      crot(:,2) = cross(crot(:,3),crot(:,1))

      ! transform x2 b transforming to spherical coordinates and back
      x2p = matmul(x2,crot)
      call tosphere(x2p,r2,asph)
      rf = d
      phf = 0.5d0 * pi - ang
      thf = asph(2) + dieh
      xf(1) = rf * cos(phf) * cos(thf)
      xf(2) = rf * cos(phf) * sin(thf)
      xf(3) = rf * sin(phf)
      call matinv(crot,3)
      zmat_step = matmul(xf,crot) + x0

    end function zmat_step

  end subroutine read_zmat_geometry

  !> Determine whether a given output file (.scf.out or .out) comes
  !> from a crystal, quantum espresso, or orca calculation.
  subroutine which_out_format(file,isformat,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, equal, lower, lgetword
    use param, only: isformat_qeout, isformat_crystal, isformat_fploout,&
       isformat_orca, isformat_aimsout
    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: isformat
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    character(len=:), allocatable :: line
    logical :: ok

    isformat = 0
    lu = fopen_read(file,ti=ti)
    if (lu < 0) return
    line = ""
    do while(getline_raw(lu,line))
       if (index(line,"--- An Ab Initio, DFT and Semiempirical electronic structure package ---") > 0) then
          isformat = isformat_orca
          exit
       else if (index(line,"EEEEEEEEEE STARTING  DATE") > 0) then
          isformat = isformat_crystal
          exit
       else if (index(line,"Program PWSCF") > 0) then
          isformat = isformat_qeout
          exit
       elseif (index(line,"Invoking FHI-aims ...") > 0) then
          isformat = isformat_aimsout
          exit
       elseif (index(line,"FULL-POTENTIAL LOCAL-ORBITAL MINIMUM BASIS BANDSTRUCTURE CODE") > 0) then
          isformat = isformat_fploout
          exit
       end if
    end do

    ! determine whether the aims output contains a molecule or a crystal
    if (isformat == isformat_aimsout) then
       isformat = 0
       do while(getline_raw(lu,line))
          if (adjustl(trim(line)) == "Input geometry:") then
             isformat = isformat_aimsout
             ok = getline_raw(lu,line)
             if (.not.ok) then
                isformat = 0
             else if (index(line,"Unit cell:") > 0) then
             elseif (index(line,"No unit cell requested.") > 0) then
             endif
             exit
          end if
       end do
    end if
    call fclose(lu)

  end subroutine which_out_format

  !> Determine whether a given input file (.in) comes from an FHIaims
  !> or a quantum espresso calculation.
  subroutine which_in_format(file,isformat,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, equal, lgetword, lower
    use param, only: isformat_qein, isformat_aimsin
    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: isformat
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp
    character(len=:), allocatable :: line, word, lline

    ! parse the input file looking for mandatory keywords
    isformat = 0
    lu = fopen_read(file,ti=ti)
    if (lu < 0) return
    line = ""
    do while(getline_raw(lu,line))
       if (len_trim(line) == 0) cycle
       lp = 1
       lline = lower(line)
       word = lgetword(line,lp)
       if (word(1:1) == "#") cycle
       if ((index(lline,"&control") > 0) .or. (index(lline,"&system") > 0)) then
          isformat = isformat_qein
          exit
       end if
       if (equal(word,"atom").or.equal(word,"atom_frac").or.equal(word,"lattice_vector")) then
          isformat = isformat_aimsin
          exit
       end if
    end do

    ! If this is an FHI-aims input, determine whether it is a molecule or a crystal
    if (isformat == isformat_aimsin) then
       rewind(lu)
       do while(getline_raw(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,"lattice_vector")) then
             exit
          end if
       end do
    end if

    ! clean up
    call fclose(lu)

  end subroutine which_in_format

  !> Convert a cif-file-style string (x,y,z) to a symmetry operation.
  function string_to_symop(str,errmsg) result(rot0)
    use arithmetic, only: eval
    use tools_io, only: lower
    character*(*), intent(in) :: str
    character(len=:), allocatable, intent(out) :: errmsg
    real*8 :: rot0(3,4)

    character(len=:), allocatable :: sym, tok, tok0
    integer :: i, j, idx

    errmsg = ""
    rot0 = 0d0
    sym = str
    sym = trim(adjustl(lower(sym))) // ","
    do i = 1, 3
       ! extract the next token
       idx = index(sym,",")
       if (idx == 0) then
          errmsg = "Error reading symmetry operation."
          return
       end if
       tok0 = sym(1:idx-1)
       sym = sym(idx+1:)

       ! the translation component
       tok = tok0
       call replace(tok,'x',0)
       call replace(tok,'y',0)
       call replace(tok,'z',0)
       rot0(i,4) = eval(tok,errmsg)
       if (len_trim(errmsg) > 0) then
          errmsg = "Error reading symmetry operation: " // trim(errmsg)
          return
       end if

       ! the other components
       do j = 1, 3
          tok = tok0
          if (j == 1) then
             call replace(tok,'x',1)
             call replace(tok,'y',0)
             call replace(tok,'z',0)
          elseif (j == 2) then
             call replace(tok,'x',0)
             call replace(tok,'y',1)
             call replace(tok,'z',0)
          elseif (j == 3) then
             call replace(tok,'x',0)
             call replace(tok,'y',0)
             call replace(tok,'z',1)
          end if
          rot0(i,j) = eval(tok,errmsg) - rot0(i,4)
          if (len_trim(errmsg) > 0) then
             errmsg = "Error reading symmetry operation: " // trim(errmsg)
             return
          end if
       enddo
    enddo

  contains
    ! Replace character ch in string tok with integer ival.
    subroutine replace(tok,ch,ival)
      use tools_io, only: string
      character(len=:), allocatable, intent(inout) :: tok
      character*1, intent(in) :: ch
      integer, intent(in) :: ival

      integer :: idx

      do while (.true.)
         idx = index(tok,ch)
         if (idx == 0) return
         tok = tok(1:idx-1) // "(" // string(ival) // ")" // tok(idx+1:)
      end do

    end subroutine replace

  end function string_to_symop

end submodule proc
