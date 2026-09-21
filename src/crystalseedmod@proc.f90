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

  ! Interim record for one structure, filled by the file readers and
  ! converted into a crystalseed by rawseed_to_seed.
  integer, parameter :: spc_label = 1 !< species = first appearance of the (case-sensitive) label
  integer, parameter :: spc_label_nocase = 2 !< species = first appearance of the label, case-insensitive
  integer, parameter :: spc_z = 3 !< species = first appearance of the atomic number
  integer, parameter :: spc_zsorted = 4 !< species = atomic numbers present, in increasing Z, named by nameguess
  integer, parameter :: spc_given = 5 !< species table and is(:) given in the record (z(:) unused)
  integer, parameter :: lunit_bohr = 0 !< lengths (cell, Cartesian coordinates) in bohr
  integer, parameter :: lunit_ang = 1 !< lengths in angstrom
  integer, parameter :: sym_none = 0 !< no symmetry information
  integer, parameter :: sym_spg = 1 !< symmetry from the GULP space group string (spg, origin, shift)
  integer, parameter :: sym_ops = 2 !< symmetry from explicit operations (neqv, ncv, rotm, cen)

  ! bookkeeping for geometries read from GULP output files
  type :: gulpout_info
     logical :: centred = .false. !< centred conventional cell (full cell parameters given)
     real*8 :: pappl = huge(1d0) !< applied pressure (GPa)
     real*8 :: dede(3) = 0d0 !< strain derivatives dE/de1..3 (eV)
     logical :: havedede = .false.
     real*8 :: vol = 0d0 !< primitive cell volume (ang^3)
     real*8 :: stress(3) = 0d0 !< diagonal stress (GPa)
     logical :: havestress = .false.
     logical :: havefinalcell = .false.
  end type gulpout_info

  type :: rawseed
     ! atoms
     integer :: nat = 0 !< number of atoms
     real*8, allocatable :: x(:,:) !< coordinates (fractional, or Cartesian if iscart)
     integer, allocatable :: z(:) !< atomic numbers (not used with spc_given)
     character*10, allocatable :: atname(:) !< atom labels
     real*8, allocatable :: occ(:) !< site occupancies (optional)
     character*10, allocatable :: spcname(:) !< per-atom species key/name when different from the label (optional)
     logical, allocatable :: isfrac(:) !< per-atom fractional flag overriding iscart (optional)
     ! species
     integer :: spcmode = spc_label !< species policy (spc_*)
     integer :: nspc = 0 !< number of species (spc_given)
     type(species), allocatable :: spc(:) !< species table (spc_given)
     integer, allocatable :: is(:) !< species index per atom (spc_given)
     ! cell
     integer :: ndim = -1 !< 0 = molecule, 3 = periodic (GULP dimensionality; the molecule switch)
     integer :: cellmode = 0 !< 0 = none, 1 = aa/bb, 2 = rv
     real*8 :: aa(3) = 0d0 !< cell lengths
     real*8 :: bb(3) = 0d0 !< cell angles (degrees)
     real*8 :: rv(3,3) = 0d0 !< lattice vectors, rv(:,i) = vector i
     integer :: lunit = lunit_ang !< unit of lengths (lunit_*)
     logical :: iscart = .false. !< coordinates are Cartesian
     logical :: wrap = .false. !< wrap fractional coordinates into [0,1)
     ! molecular fields
     logical :: ismol = .false. !< a molecule regardless of the mol argument
     real*8 :: border = 0d0 !< molecular cell border (bohr)
     logical :: cubic = .false. !< cubic molecular cell
     logical :: havex0 = .false. !< molecular cell origin given (at zero)
     ! symmetry
     integer :: symmode = sym_none !< symmetry policy (sym_*)
     character(len=:), allocatable :: spg !< space group (HM symbol or number); empty = P1
     integer :: origin = 1 !< origin choice (1 or 2)
     logical :: haveshift = .false. !< origin shift given
     logical :: shiftfromout = .false. !< the shift was read from an output file echo
     real*8 :: shift(3) = 0d0 !< origin shift (fractional)
     integer :: neqv = 0 !< number of symmetry operations (sym_ops)
     integer :: ncv = 0 !< number of centering vectors (sym_ops)
     real*8, allocatable :: rotm(:,:,:) !< symmetry operations (sym_ops)
     real*8, allocatable :: cen(:,:) !< centering vectors (sym_ops)
     ! properties and naming
     real*8 :: energy = huge(1d0) !< energy (Hartree)
     real*8 :: pressure = huge(1d0) !< pressure (GPa)
     character(len=:), allocatable :: name !< if allocated, the seed name is file|name
     character(len=:), allocatable :: tag !< if allocated, replaces the energy in the trajectory names
     logical :: isfinal = .false. !< final (optimised) geometry of a trajectory
     type(gulpout_info) :: gulp !< GULP output bookkeeping
  end type rawseed

  integer, parameter :: mxtok = 60 !< maximum number of tokens per line

  !xx! private subroutines
  ! subroutine read_all_cif(file,mol,errmsg,nseed,mseed,seed0,dblock,ti)
  ! subroutine read_all_mol2(file,errmsg,nseed,mseed,seed0,name,ti)
  ! subroutine read_all_sdf(file,errmsg,nseed,mseed,seed0,name,ti)
  ! subroutine read_all_qeout(nseed,seed,file,mol,istruct,errmsg,ti)
  ! subroutine read_all_xyz(nseed,seed,file,mol,errmsg,seed0,ti)
  ! subroutine read_all_log(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_aimsout(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_castep_geom(nseed,seed,file,errmsg,ti)
  ! subroutine read_all_gulpin(nseed,seed,file,mol,istruct,errmsg,ti)
  ! subroutine read_all_gulpout(nseed,seed,file,mol,istruct,errmsg,ti)
  ! subroutine gulp_detect_ismol(file,isformat,ismol,ti)
  ! subroutine rawseed_to_seed(rec,seed,mol,file,isformat,errmsg)
  ! subroutine gulp_spg_ops(rec,neqv,ncv,rotm,cen,found,errmsg)
  ! subroutine rawseed_select(recs,nrec,istruct,mol,file,isformat,nseed,seed,errmsg)
  ! subroutine rawseed_add_atom(rec,z,x,atname,occ,spcname,isfrac)
  ! subroutine rawseed_grow(recs,n)
  ! subroutine read_zmat_geometry(file,nat,x,z,name,errmsg,ti)
  ! function which_out_format(file,ti)
  ! subroutine which_in_format(file,isformat,ti)
  ! function string_to_symop(str)
  ! subroutine qe_read_namelists(lu,calculation,ibrav,celldm,nat,ntyp,space_group,uniqueb,rhombohedral,origin_choice,errmsg)

contains

  !> Deallocate arrays in the seed
  module subroutine seed_end(seed)
    use param, only: isformat_r_unknown
    class(crystalseed), intent(inout) :: seed !< Crystal seed output

    seed%isused = .false.
    seed%file = ""
    seed%name = ""
    seed%libname = ""
    seed%isformat = isformat_r_unknown
    seed%nat = 0
    if (allocated(seed%x)) deallocate(seed%x)
    if (allocated(seed%is)) deallocate(seed%is)
    if (allocated(seed%atname)) deallocate(seed%atname)
    if (allocated(seed%occ)) deallocate(seed%occ)
    seed%neqlist = .false.
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
    seed%havebonds = .false.
    if (allocated(seed%nstar)) deallocate(seed%nstar)
    seed%ismolecule = .false.
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%energy = huge(1d0)
    seed%pressure = huge(1d0)

  end subroutine seed_end

  !> Check the seed fields for internal consistency. These checks guard
  !> against ill-formed seeds produced by the reader routines (programmer
  !> errors), not against bad user input. If errmsg is present, the first
  !> violation found is returned there and the routine returns; otherwise a
  !> violation aborts the program with faterr.
  module subroutine seed_check(seed,errmsg)
    use tools_io, only: ferror, faterr
    use tools_math, only: m_x2c_from_cellpar
    class(crystalseed), intent(in) :: seed !< Crystal seed
    character(len=:), allocatable, intent(out), optional :: errmsg !< error message (faterr if absent)

    character(len=:), allocatable :: msg
    integer :: i, ier
    real*8 :: xdum(3,3)

    msg = ""

    ! the seed must be in use
    if (.not.seed%isused) then
       msg = "seed is not in use"
       goto 999
    end if

    ! useabr must be in range; useabr==0 only for molecules (auto bounding box)
    if (seed%useabr < 0 .or. seed%useabr > 2) then
       msg = "unknown useabr (must be 0, 1, or 2)"
       goto 999
    end if
    if (seed%useabr == 0 .and. .not.seed%ismolecule) then
       msg = "useabr=0 (auto bounding box) is only valid for a molecule"
       goto 999
    end if

    ! The cell parameters must describe a realizable cell. This is checked
    ! here, before struct_new wipes the caller's crystal with c%init(), so a
    ! bad cell is reported without destroying the structure being edited.
    ! Uses the same math as the builder, so the two cannot disagree.
    if (seed%useabr == 1) then
       if (any(seed%aa <= 0d0)) then
          msg = "invalid cell: cell lengths must be positive"
          goto 999
       end if
       if (any(seed%bb <= 0d0) .or. any(seed%bb >= 180d0)) then
          msg = "invalid cell: cell angles must be between 0 and 180 degrees"
          goto 999
       end if
       xdum = m_x2c_from_cellpar(seed%aa,seed%bb,ier)
       if (ier /= 0) then
          msg = "invalid cell: these cell angles do not describe a realizable cell"
          goto 999
       end if
    end if

    ! atom list kind: a non-equivalent atom list requires symmetry and a crystal.
    ! A complete list (neqlist=.false.) is allowed with or without symmetry.
    if (seed%neqlist .and. (seed%havesym <= 0 .or. seed%ismolecule)) then
       msg = "non-equivalent atom list (neqlist) requires symmetry (havesym>0) and a crystal"
       goto 999
    end if

    ! molecules carry no crystallographic symmetry
    if (seed%ismolecule) then
       if (seed%havesym > 0) then
          msg = "a molecule cannot have crystallographic symmetry (havesym>0)"
          goto 999
       end if
       if (seed%checkrepeats) then
          msg = "a molecule cannot use checkrepeats"
          goto 999
       end if
    end if

    ! checkrepeats (remove symmetry-duplicate atoms before building) requires
    ! symmetry; it applies to both a complete and a non-equivalent atom list
    if (seed%checkrepeats .and. seed%havesym <= 0) then
       msg = "checkrepeats requires symmetry (havesym>0)"
       goto 999
    end if

    ! the symmetry info must be present and consistent when symmetry is claimed
    if (seed%havesym > 0) then
       if (seed%neqv < 1 .or. seed%ncv < 1) then
          msg = "havesym>0 but neqv or ncv is not positive"
          goto 999
       end if
       if (.not.allocated(seed%cen) .or. .not.allocated(seed%rotm)) then
          msg = "havesym>0 but cen or rotm is not allocated"
          goto 999
       end if
       if (size(seed%cen,2) < seed%ncv .or. size(seed%rotm,3) < seed%neqv) then
          msg = "havesym>0 but cen or rotm is too small for neqv/ncv"
          goto 999
       end if
    end if

    ! bonds are indexed over the complete cell list, so they require the
    ! complete list (neqlist=.false.); they are compatible with symmetry (the
    ! reduction in struct_new keeps the complete-list order)
    if (seed%havebonds) then
       if (seed%neqlist) then
          msg = "incompatible fields in seed: bonds and non-equivalent atom list"
          goto 999
       end if
       if (.not.allocated(seed%nstar)) then
          msg = "havebonds but nstar is not allocated"
          goto 999
       end if
       if (size(seed%nstar,1) /= seed%nat) then
          msg = "havebonds but nstar size does not match the number of atoms"
          goto 999
       end if
    end if

    ! findsym must be in range
    if (seed%findsym < -1 .or. seed%findsym > 1) then
       msg = "unknown findsym (must be -1, 0, or 1)"
       goto 999
    end if

    ! atom and species sanity
    if (seed%nat > 0) then
       if (.not.allocated(seed%x) .or. .not.allocated(seed%is) .or. .not.allocated(seed%atname)) then
          msg = "nat>0 but x, is, or atname is not allocated"
          goto 999
       end if
       if (size(seed%x,2) < seed%nat .or. size(seed%is,1) < seed%nat .or. size(seed%atname,1) < seed%nat) then
          msg = "nat>0 but x, is, or atname is too small"
          goto 999
       end if
       if (seed%nspc < 1 .or. .not.allocated(seed%spc)) then
          msg = "nat>0 but no species are defined"
          goto 999
       end if
       do i = 1, seed%nat
          if (seed%is(i) < 1 .or. seed%is(i) > seed%nspc) then
             msg = "atom with species index out of range"
             goto 999
          end if
       end do
       if (allocated(seed%occ)) then
          if (size(seed%occ,1) < seed%nat) then
             msg = "occ is allocated but too small"
             goto 999
          end if
          do i = 1, seed%nat
             if (seed%occ(i) <= 0d0 .or. seed%occ(i) > 1d0+1d-10) then
                msg = "atom with occupancy out of range (must be in (0,1])"
                goto 999
             end if
          end do
       end if
    end if

999 continue
    if (present(errmsg)) errmsg = msg
    if (len_trim(msg) > 0 .and. .not.present(errmsg)) &
       call ferror("seed_check",msg,faterr)

  end subroutine seed_check

  !> Parse a crystal environment
  module subroutine parse_crystal_env(seed,lu,oksyn,fromlib)
    use spglib, only: spg_get_hall_number_from_symbol, spg_get_symmetry_from_database
    use global, only: eval_next, dunit0, iunit, iunit_isdef, iunit_bohr
    use tools_math, only: matinv
    use tools_io, only: uin, getline, ucopy, lgetword, equal, ferror, faterr,&
       getword, lower, isinteger, isreal, string, nameguess, zatguess, equali
    use param, only: bohrtoa, isformat_r_from_input
    use types, only: realloc

    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    integer, intent(in) :: lu !< Logical unit for input
    logical, intent(out) :: oksyn !< Was there a syntax error?
    logical, intent(in), optional :: fromlib !< from library (use native units, ignore UNITS)

    character(len=:), allocatable :: word, aux, line, uword, name, errmsg
    character*255, allocatable :: sline(:)
    integer :: i, j, k, lp, nsline, luout, iat, lp2, iunit0, it
    integer :: hnum, ier
    real*8 :: rmat(3,3), scal, ascal, x(3), xdif(3), rocc
    real*8, allocatable :: occ(:)
    logical :: ok, goodspg, islib
    logical :: isset
    real*8 :: rot(3,4)

    call seed%end()
    ! library structures are stored in bohr; ignore the user's UNITS setting
    islib = .false.
    if (present(fromlib)) islib = fromlib
    if (iunit_isdef .or. islib) then
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
    allocate(seed%x(3,10),seed%is(10),seed%spc(2),seed%atname(10),occ(10))
    do while (getline(lu,line,ucopy=luout))
       lp = 1
       uword = getword(line,lp)
       word = lower(uword)
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
          ascal = 1d0/dunit0(iunit0)
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
             call ferror('parse_crystal_env','Invalid lattice vectors (matrix not invertible)',&
                faterr,line,syntax=.true.)
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
             call realloc(seed%atname,2*seed%nat)
             call realloc(occ,2*seed%nat)
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
                name = string(uword)
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
          seed%atname(seed%nat) = name

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
          end if
          seed%is(seed%nat) = it

          ! optional trailing fields: unit keyword and/or an occupancy value
          ! (a bare number, taken as the last column; default 1)
          rocc = 1d0
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
             else if (len_trim(word) > 0 .and. isreal(rocc,word)) then
                if (rocc <= 0d0 .or. rocc > 1d0) then
                   call ferror('parse_crystal_env','Occupancy must be in the range (0,1]',faterr,line,syntax=.true.)
                   return
                end if
             else if (len_trim(word) > 0) then
                call ferror('parse_crystal_env','Unknown keyword in atom input',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do
          occ(seed%nat) = rocc
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
                   call realloc(occ,2*seed%nat)
                end if
                seed%is(seed%nat) = seed%is(j)
                seed%x(:,seed%nat) = x
                occ(seed%nat) = occ(j)
             endif
          end do
       end do
       deallocate(sline)
    end if
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%spc,seed%nspc)

    ! store the occupancies (kept only if some site is partial)
    call seed_set_occ(seed,occ,seed%nat)

    ! symmetry (the crystal environment is never a molecule)
    if (goodspg) then
       seed%havesym = 1
       seed%checkrepeats = .false.
       seed%findsym = 0
       seed%neqlist = .true.
    else
       seed%havesym = 0
       seed%checkrepeats = .false.
       seed%findsym = -1
       seed%neqlist = .false.
    end if
    oksyn = .true.

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = .false.
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = ""
    seed%name = "<input>"
    seed%isformat = isformat_r_from_input

  end subroutine parse_crystal_env

  !> Parse a molecule environment
  module subroutine parse_molecule_env(seed,lu,oksyn,fromlib)
    use global, only: rborder_def, eval_next, dunit0, iunit, iunit_ang, iunit_isdef
    use tools_io, only: uin, ucopy, getline, lgetword, equal, ferror, faterr,&
       string, isinteger, isreal, nameguess, getword, zatguess, equali, lower
    use param, only: bohrtoa, isformat_r_from_input
    use types, only: realloc

    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    integer, intent(in) :: lu !< Logical unit for input
    logical, intent(out) :: oksyn !< Was there a syntax error?
    logical, intent(in), optional :: fromlib !< from library (use native units, ignore UNITS)

    character(len=:), allocatable :: uword, word, aux, line, name
    integer :: lp, lp2, luout, iat, iunit0, it, i
    real*8 :: rborder, rocc
    real*8, allocatable :: occ(:)
    logical :: ok, docube, isset, islib

    call seed%end()
    ! library molecules are stored in angstrom; ignore the user's UNITS setting
    islib = .false.
    if (present(fromlib)) islib = fromlib
    if (iunit_isdef .or. islib) then
       iunit0 = iunit_ang
    else
       iunit0 = iunit
    end if

    ok = .false.
    docube = .false.
    rborder = rborder_def
    seed%nat = 0
    seed%nspc = 0
    allocate(seed%x(3,10),seed%is(10),seed%spc(2),seed%atname(10),occ(10))
    if (lu == uin) then
       luout = ucopy
    else
       luout = -1
    endif
    do while (getline(lu,line,ucopy=luout))
       lp = 1
       uword = getword(line,lp)
       word = lower(uword)

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
             call realloc(seed%atname,2*seed%nat)
             call realloc(occ,2*seed%nat)
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
                name = string(uword)
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
          seed%atname(seed%nat) = name

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
          end if
          seed%is(seed%nat) = it

          ! optional trailing fields: unit keyword and/or an occupancy value
          ! (a bare number, taken as the last column; default 1)
          isset = .false.
          rocc = 1d0
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'ang') .or. equal(word,'angstrom')) then
                isset = .true.
                seed%x(:,seed%nat) = seed%x(:,seed%nat) / bohrtoa
             else if (equal(word,'bohr') .or. equal(word,'au')) then
                isset = .true.
             else if (len_trim(word) > 0 .and. isreal(rocc,word)) then
                if (rocc <= 0d0 .or. rocc > 1d0) then
                   call ferror('parse_molecule_input','Occupancy must be in the range (0,1]',faterr,line,syntax=.true.)
                   return
                end if
             else if (len_trim(word) > 0) then
                call ferror('parse_molecule_input','Unknown extra keyword in atomic input',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do
          occ(seed%nat) = rocc
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

    ! store the occupancies (kept only if some site is partial)
    call seed_set_occ(seed,occ,seed%nat)

    oksyn = .true.
    seed%useabr = 0

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = .true.
    seed%cubic = docube
    seed%border = rborder
    seed%havex0 = .false.
    seed%molx0 = 0d0
    seed%file = ""
    seed%name = "<input>"
    seed%isformat = isformat_r_from_input

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
    use param, only: isformat_r_from_library
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
       call seed%parse_crystal_env(lu,ok,fromlib=.true.)
    elseif (equal(word,"molecule")) then
       call seed%parse_molecule_env(lu,ok,fromlib=.true.)
    else
       call fclose(lu)
       call ferror("read_library","error parsing the crystal/molecule keyword",faterr,syntax=.true.)
       return
    endif
    seed%file = ""
    seed%libname = trim(line)
    seed%name = trim(line) // " (library)"
    seed%isformat = isformat_r_from_library

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
    use param, only: isformat_r_derived
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
    if (allocated(seed%is)) deallocate(seed%atname)
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
    do i = 1, seed%nat
       seed%x(:,i) = fr%at(iord(i))%r
       seed%is(i) = fr%at(iord(i))%is
       seed%atname(i) = fr%spc(fr%at(iord(i))%is)%name
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
    seed%isformat = isformat_r_derived

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
          seed%atname(nnat) = seed%atname(i)
          if (allocated(seed%occ)) seed%occ(nnat) = seed%occ(i)
       end if
    end do
    seed%nat = nnat
    call realloc(seed%x,3,nnat)
    call realloc(seed%is,nnat)
    call realloc(seed%atname,nnat)
    if (allocated(seed%occ)) call realloc(seed%occ,nnat)

  end subroutine strip_hydrogens

  !> Detect the file format of a file from the extension and read a
  !> crystal seed from it. If mol0 == 1, force a molecule. If mol0 ==
  !> 0, force a crystal. If mol0 == -1, let the format detection
  !> decide. If the read was successful, return an empty error
  !> message
  module subroutine read_any_file(seed,file,mol0,errmsg,ti)
    use global, only: rborder_def
    use param, only: isformat_r_cif, isformat_r_shelx, isformat_r_f21,&
       isformat_r_cube, isformat_r_bincube, isformat_r_struct, isformat_r_abinit,&
       isformat_r_elk, isformat_r_fploout,&
       isformat_r_qein, isformat_r_qeout, isformat_r_crystal, isformat_r_xyz,&
       isformat_r_wfn, isformat_r_wfx, isformat_r_fchk, isformat_r_molden,&
       isformat_r_gaussian, isformat_r_siesta, isformat_r_xsf, isformat_r_gen,&
       isformat_r_vasp, isformat_r_pwc, isformat_r_axsf, isformat_r_dat,&
       isformat_r_pgout, isformat_r_orca, isformat_r_dmain, isformat_r_aimsin,&
       isformat_r_aimsout, isformat_r_tinkerfrac, isformat_r_gjf, isformat_r_zmat,&
       isformat_r_magres, isformat_r_alamode, isformat_r_akaikkr, isformat_r_xband,&
       isformat_r_sdf, isformat_r_castepcell, isformat_r_gulpin, isformat_r_gulpout,&
       isformat_r_castepphonon, isformat_r_castepgeom, isformat_r_mol2, isformat_r_pdb
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    integer, intent(in) :: mol0
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: isformat
    logical :: mol, hastypes

    ! detect the format for this file
    errmsg = ""
    call struct_detect_read_format(file,isformat,ti=ti)

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
    if (isformat == isformat_r_cif) then
       call seed%read_cif(file," ",mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_mol2) then
       call seed%read_mol2(file,rborder_def,.false.," ",errmsg,ti=ti)

    elseif (isformat == isformat_r_sdf) then
       call seed%read_sdf(file,rborder_def,.false.,errmsg,0,ti=ti)

    elseif (isformat == isformat_r_magres) then
       call seed%read_magres(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_alamode) then
       call seed%read_alamode(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_akaikkr) then
       call seed%read_akaikkr(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_xband) then
       call seed%read_xband(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_shelx) then
       call seed%read_shelx(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_f21) then
       call seed%read_f21(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_cube) then
       call seed%read_cube(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_bincube) then
       call seed%read_bincube(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_struct) then
       call seed%read_wien(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_vasp) then
       call seed%read_vasp(file,mol,hastypes,errmsg,ti=ti)
       if (len_trim(errmsg) == 0 .and..not.hastypes) then
          errmsg = "Atom types not found in VASP file"
          return
       end if

    elseif (isformat == isformat_r_abinit) then
       call seed%read_abinit(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_elk) then
       call seed%read_elk(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_qeout) then
       call seed%read_qeout(file,mol,0,errmsg,ti=ti)

    elseif (isformat == isformat_r_gulpin) then
       call seed%read_gulpin(file,mol,0,errmsg,ti=ti)

    elseif (isformat == isformat_r_gulpout) then
       call seed%read_gulpout(file,mol,0,errmsg,ti=ti)

    elseif (isformat == isformat_r_crystal) then
       call seed%read_crystalout(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_fploout) then
       call seed%read_fploout(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_qein) then
       call seed%read_qein(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_wfn.or.&
       isformat == isformat_r_wfx.or.isformat == isformat_r_fchk.or.&
       isformat == isformat_r_molden.or.isformat == isformat_r_gaussian.or.&
       isformat == isformat_r_dat.or.isformat == isformat_r_pgout.or.&
       isformat == isformat_r_orca.or.isformat == isformat_r_gjf.or.&
       isformat == isformat_r_zmat) then
       call seed%read_mol(file,isformat,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_r_xyz) then
       call seed%read_xyz(file,mol,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_r_pdb) then
       call seed%read_pdb(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_siesta) then
       call seed%read_siesta(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_castepcell) then
       call seed%read_castep_cell(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_castepphonon) then
       call seed%read_castep_phonon(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_castepgeom) then
       call seed%read_castep_geom(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_xsf) then
       call seed%read_xsf(file,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_r_gen) then
       call seed%read_dftbp(file,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_r_pwc) then
       call seed%read_pwc(file,mol,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_r_axsf) then
       call seed%read_axsf(file,1,0d0,rborder_def,.false.,errmsg,ti=ti)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_r_dmain) then
       call seed%read_dmain(file,mol,errmsg,ti=ti)

    elseif (isformat == isformat_r_aimsin) then
       call seed%read_aimsin(file,mol,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_r_aimsout) then
       call seed%read_aimsout(file,mol,rborder_def,.false.,errmsg,ti=ti)

    elseif (isformat == isformat_r_tinkerfrac) then
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

  !> Read a structure seed from an sdf/mol file. If id is not present
  !> or zero, return the first molecule in the mol2 file, otherwise
  !> return the molecule with that number in the sequence. If mol,
  !> force a molecule/crystal system. If error, return non-empty
  !> errmsg.
  module subroutine read_sdf(seed,file,rborder,docube,errmsg,id,ti)
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character(len=:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: id
    type(thread_info), intent(in), optional :: ti

    call seed%end()
    call read_all_sdf(file,errmsg,seed0=seed,id=id,ti=ti)
    seed%border = rborder
    seed%cubic = docube

  end subroutine read_sdf

  !> Read the structure from a res or ins (shelx) file
  module subroutine read_shelx(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, isreal, isinteger,&
       lower, zatguess, fclose, getword
    use param, only: eyet, eye, bohrtoa, isformat_r_shelx, hartokjmol
    use types, only: realloc
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, ilat, ll
    logical :: ok, iscent, havecell, found
    character(len=:), allocatable :: uword, word, line, aux, word2
    real*8 :: raux, rot0(3,4), sof, pp, x0(3), xd(3)
    integer :: i, j, k, io, it, n, mm, mult
    integer :: iz
    integer :: lncv, nfvar
    real*8, allocatable :: lcen(:,:), occfac(:), fvar(:)

    real*8, parameter :: eps = 1d-5

    ! file and seed name
    call seed%end()
    seed%file = file
    seed%name = file
    seed%isformat = isformat_r_shelx
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
    if (.not.allocated(seed%atname)) allocate(seed%atname(10))
    allocate(occfac(10))

    ! free variables (fvar(1) is the overall scale); site occupation factors
    ! (sof) may reference them
    nfvar = 0
    allocate(fvar(1))
    fvar = 1d0

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
       uword = getword(line,lp)
       word = lower(uword)
       if (len_trim(word) > 4) word = word(1:4)
       if (equal(word,"titl")) then
          ! try to read the energy (output of mol-cspy, in kJ/mol per molecule)
          word2 = getword(line,lp)
          ll = min(len(word2),len(file))
          if (equal(word2(1:ll-1),file(1:ll-1))) then
             ok = isreal(raux,line,lp)
             if (ok) seed%energy = raux / hartokjmol
          end if
       elseif (equal(word,"cell")) then
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
       elseif (equal(word,"fvar")) then
          ! read the free variables (used to model sof of disordered atoms)
          do while (.true.)
             ok = isreal(raux,line,lp)
             if (.not.ok) exit
             nfvar = nfvar + 1
             if (nfvar > size(fvar,1)) call realloc(fvar,2*nfvar)
             fvar(nfvar) = raux
          end do
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
          equal(word,"free").or.equal(word,"grid").or.&
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
          equal(word,"temp").or.equal(word,"time").or.equal(word,"twin").or.&
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
             call realloc(seed%atname,2*seed%nat)
             call realloc(occfac,2*seed%nat)
          end if
          ok = isinteger(iz,line,lp)
          ok = ok .and. isreal(seed%x(1,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(2,seed%nat),line,lp)
          ok = ok .and. isreal(seed%x(3,seed%nat),line,lp)
          if (.not.ok) then
             ! not a valid atom card: could be a multi-line TITL continuation
             ! (e.g. Olex2 writes extra description lines) or an unrecognized
             ! keyword. Skip the line rather than treating it as an atom.
             seed%nat = seed%nat - 1
             cycle
          end if
          if (iz < 1 .or. iz > seed%nspc) then
             errmsg = "Atom type not found in SFAC list."
             goto 999
          end if
          seed%is(seed%nat) = iz
          seed%atname(seed%nat) = uword

          ! site occupation factor (sof): optional 4th number after the
          ! coordinates. Decode the SHELX free-variable convention m*10+p:
          ! m=0/1 fixed at p; m>1 uses p*fvar(m); m<-1 uses |p|*(1-fvar(-m)).
          ! For m<-1 both mm and pp come out negative (e.g. sof=-21 -> mm=-2,
          ! pp=-1), so the magnitude is -pp. The special-position multiplicity
          ! factor is removed below.
          occfac(seed%nat) = 1d0
          if (isreal(sof,line,lp)) then
             mm = nint(sof/10d0)
             pp = sof - 10d0*mm
             if (mm == 0 .or. mm == 1) then
                occfac(seed%nat) = pp
             elseif (mm > 1) then
                if (mm <= nfvar) then
                   occfac(seed%nat) = pp * fvar(mm)
                else
                   occfac(seed%nat) = pp
                end if
             elseif (mm < -1) then
                if (-mm <= nfvar) then
                   occfac(seed%nat) = -pp * (1d0 - fvar(-mm))
                else
                   occfac(seed%nat) = -pp
                end if
             end if
          end if
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
    call realloc(seed%atname,seed%nat)
    call realloc(occfac,seed%nat)
    call realloc(seed%cen,3,seed%ncv)

    ! use the symmetry in this file (a molecule does not apply crystallographic
    ! symmetry)
    if (.not.mol) then
       seed%havesym = 1
       seed%checkrepeats = .true.
       seed%neqlist = .true.
    else
       seed%havesym = 0
       seed%checkrepeats = .false.
       seed%neqlist = .false.
    end if
    seed%findsym = -1
    call realloc(seed%rotm,3,4,seed%neqv)
    call realloc(seed%cen,3,seed%ncv)

    ! recover the chemical occupancy from the SHELX site occupation factor. For
    ! an atom on a special position the sof also carries the ratio of the site
    ! multiplicity to the general multiplicity (neqv*ncv), which we divide out.
    do i = 1, seed%nat
       if (.not.mol) then
          ! count the distinct images (the site multiplicity)
          mult = 0
          do io = 1, seed%neqv
             do it = 1, seed%ncv
                x0 = matmul(seed%rotm(1:3,1:3,io),seed%x(:,i)) + seed%rotm(:,4,io) + seed%cen(:,it)
                found = .false.
                inner: do k = 1, io
                   do j = 1, seed%ncv
                      if (k == io .and. j >= it) exit inner
                      xd = matmul(seed%rotm(1:3,1:3,k),seed%x(:,i)) + seed%rotm(:,4,k) + seed%cen(:,j)
                      xd = x0 - xd
                      xd = xd - nint(xd)
                      if (all(abs(xd) < eps)) then
                         found = .true.
                         exit inner
                      end if
                   end do
                end do inner
                if (.not.found) mult = mult + 1
             end do
          end do
          if (mult > 0) then
             occfac(i) = occfac(i) * real(seed%neqv*seed%ncv,8) / real(mult,8)
          end if
       end if
       ! snap a nearly-full occupancy to exactly 1 (seed_set_occ does the clamp)
       if (abs(occfac(i) - 1d0) < 1d-3) occfac(i) = 1d0
    end do

    ! store the occupancies (kept only if some site is partial)
    call seed_set_occ(seed,occfac,seed%nat)

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

  !> Read the structure from a magres file, for structural and NMR info.
  module subroutine read_magres(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, lower, fclose, lgetword, equal,&
       isreal, zatguess, nameguess, getword, equalsub
    use param, only: isformat_r_magres
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok, havelat, haveatoms
    integer :: lu, lp, iz
    character(len=:), allocatable :: line, word, word2
    real*8 :: x(3)
    type(rawseed) :: rec

    ! file and seed name
    call seed%end()
    errmsg = ""

    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! locate the <atoms> tag
    ok = .false.
    do while (getline_raw(lu,line,.false.))
       line = adjustl(lower(line))
       if (equalsub(line,1,7,"<atoms>") .or. equalsub(line,1,7,"[atoms]")) then
          ok = .true.
          exit
       end if
    end do
    if (.not.ok) then
       errmsg = "<atoms> tag not found"
       goto 999
    end if

    ! initialize (angstrom, Cartesian coordinates)
    rec%spcmode = spc_z
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_ang
    rec%iscart = .true.
    havelat = .false.
    haveatoms = .false.
    do while (getline_raw(lu,line,.false.))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"units")) then
          word = lgetword(line,lp)
          word2 = lgetword(line,lp)
          ok = equal(word,"lattice") .or. equal(word,"atom")
          ok = ok .and. equal(word2,"angstrom")
          if (.not.ok) then
             errmsg = "unknown units"
             goto 999
          end if
       elseif (equal(word,"lattice")) then
          ok = isreal(rec%rv(1,1),line,lp)
          ok = ok .and. isreal(rec%rv(2,1),line,lp)
          ok = ok .and. isreal(rec%rv(3,1),line,lp)
          ok = ok .and. isreal(rec%rv(1,2),line,lp)
          ok = ok .and. isreal(rec%rv(2,2),line,lp)
          ok = ok .and. isreal(rec%rv(3,2),line,lp)
          ok = ok .and. isreal(rec%rv(1,3),line,lp)
          ok = ok .and. isreal(rec%rv(2,3),line,lp)
          ok = ok .and. isreal(rec%rv(3,3),line,lp)
          if (.not.ok) then
             errmsg = "error reading lattice vectors"
             goto 999
          end if
          havelat = .true.
       elseif (equal(word,"atom")) then
          ! get the atomic symbol
          word = getword(line,lp)
          word2 = lgetword(line,lp) ! the second symbol
          word2 = lgetword(line,lp) ! the integer ID
          iz = zatguess(word)
          if (iz <= 0) then
             errmsg = "unknown atomic symbol: " // word
             goto 999
          end if

          ! read the coordinates
          ok = isreal(x(1),line,lp)
          ok = ok .and. isreal(x(2),line,lp)
          ok = ok .and. isreal(x(3),line,lp)
          if (.not.ok) then
             errmsg = "error reading atomic coordinates"
             goto 999
          end if
          call rawseed_add_atom(rec,iz,x,word,spcname=nameguess(iz,.true.))

          ! finished
          haveatoms = .true.
       elseif (equal(word,"</atoms>") .or. equal(word,"[/atoms]")) then
          exit
       end if
    end do
    if (.not.havelat) then
       errmsg = "lattice vectors not found"
       goto 999
    end if
    if (.not.haveatoms) then
       errmsg = "atoms not found"
       goto 999
    end if

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_magres,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_magres

  !> Read the structure from an alamode input file.
  module subroutine read_alamode(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, lower, fclose, lgetword, equal,&
       isreal, zatguess, getword, getline, equalsub
    use param, only: isformat_r_alamode
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok, havelat, haveatoms
    integer :: idx, nat, nkd
    integer :: lu, lp, i
    character(len=:), allocatable :: line, word, atomstr
    real*8 :: unit
    type(rawseed) :: rec

    ! file and seed name
    call seed%end()
    errmsg = ""

    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! locate the &general tag
    ok = .false.
    do while (getline_raw(lu,line,.false.))
       line = adjustl(lower(line))
       if (equalsub(line,1,8,"&general")) then
          ok = .true.
          exit
       end if
    end do
    if (.not.ok) then
       errmsg = "&general namelist not found"
       goto 999
    end if

    ! read the &general info for nat and nkd
    errmsg = "error reading nat and nkd from &general"
    nat = 0
    nkd = 0
    do while (getline_raw(lu,line,.false.))
       line = adjustl(line)
       do while (.true.)
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,'nat')) then
             idx = index(line,"=")
             word = line(idx+1:)
             idx = index(word,';')
             if (idx > 0) word = word(:idx-1)
             read (word,*,err=999,end=999) nat
          elseif (equal(word,'nkd')) then
             idx = index(line,"=")
             word = line(idx+1:)
             idx = index(word,';')
             if (idx > 0) word = word(:idx-1)
             read (word,*,err=999,end=999) nkd
          elseif (equal(word,'kd')) then
             idx = index(line,"=")
             atomstr = line(idx+1:)
             idx = index(atomstr,';')
             if (idx > 0) atomstr = atomstr(:idx-1)
          elseif (equal(word,'/')) then
             exit
          end if

          idx = index(line,";")
          if (idx <= 0) then
             exit
          else
             line = line(idx+1:)
             if (len_trim(line) == 0) exit
          end if
       end do
    end do
    if (nat == 0) then
       errmsg = "nat not found in &general"
       goto 999
    end if
    if (nkd == 0) then
       errmsg = "nkd not found in &general"
       goto 999
    end if

    ! allocate
    rec%nat = nat
    allocate(rec%x(3,nat),rec%is(nat))

    ! rewind and re-read
    rewind(lu)
    errmsg = "error reading &cell and &position"
    havelat = .false.
    haveatoms = .false.
    do while (getline_raw(lu,line,.false.))
       line = adjustl(lower(line))
       if (equalsub(line,1,5,"&cell")) then
          ok = getline(lu,line)
          ok = ok .and. isreal(unit,line)
          if (.not.ok) goto 999

          do i = 1, 3
             lp = 1
             ok = getline(lu,line)
             ok = ok .and. isreal(rec%rv(1,i),line,lp)
             ok = ok .and. isreal(rec%rv(2,i),line,lp)
             ok = ok .and. isreal(rec%rv(3,i),line,lp)
             if (.not.ok) goto 999
          end do
          havelat = .true.

          ! fill the cell metrics (atomic units)
          rec%rv = rec%rv * unit
          rec%ndim = 3
          rec%cellmode = 2
          rec%lunit = lunit_bohr
       end if
       if (equalsub(line,1,9,"&position")) then
          do i = 1, nat
             lp = 1
             ok = getline(lu,line)
             if (.not.ok) goto 999
             read(line,*,err=999,end=999) rec%is(i), rec%x(:,i)
          end do
          haveatoms = .true.
       end if
    end do
    if (.not.havelat) then
       errmsg = "lattice vectors (&cell) not found"
       goto 999
    end if
    if (.not.haveatoms) then
       errmsg = "atomic positions (&position) not found"
       goto 999
    end if

    ! assign species
    rec%spcmode = spc_given
    rec%nspc = nkd
    allocate(rec%spc(nkd))
    lp = 1
    do i = 1, nkd
       word = getword(atomstr,lp)
       rec%spc(i)%name = word
       rec%spc(i)%z = zatguess(word)
       if (rec%spc(i)%z < 0) then
          errmsg = "unknown atom: " // word
          goto 999
       end if
    end do

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_alamode,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_alamode

  !> Read the structure from an akaikkr input file.
  module subroutine read_akaikkr(seed,file,mol,errmsg,ti)
    use tools_math, only: matinv
    use global, only: eval_next
    use tools_io, only: fopen_read, fclose, getline_raw, getword, lower,&
       isreal, isinteger
    use param, only: isformat_r_akaikkr, pi
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, i, j, k, ll, ibrav, id, ier, ncmp
    logical :: ok
    character(len=:), allocatable :: line, token
    character*3 :: brvtyp
    real*8 :: a, coa, boa, alpha, beta, gamma, cw
    real*8 :: theta, z, snrm, xc, c2x(3,3)

    real*8, parameter :: eps = 1d-6

    ! file and seed name
    call seed%end()
    seed%file = file
    seed%name = file
    seed%isformat = isformat_r_akaikkr
    errmsg = "Error reading file: " // trim(file)
    line = ""
    lp = 1

    ! open the file
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! skip go and file
    call read_next_token()
    call read_next_token()

    ! brvtyp and cell parameters
    seed%useabr = 2
    call read_next_token()
    brvtyp = lower(token)
    call read_next_token()
    ok = isreal(a,token)
    if (.not.ok) goto 999
    call read_next_token()
    ok = isreal(coa,token)
    if (.not.ok) coa = 1d0
    call read_next_token()
    ok = isreal(boa,token)
    if (.not.ok) boa = 1d0
    call read_next_token()
    ok = isreal(alpha,token)
    if (.not.ok) alpha = 90d0
    call read_next_token()
    ok = isreal(beta,token)
    if (.not.ok) beta = 90d0
    call read_next_token()
    ok = isreal(gamma,token)
    if (.not.ok) gamma = 90d0

    ! set the x2c from the Bravais type and the lattice parameters
    if (brvtyp == "fcc") then
       ibrav = 1
       boa=1d0
       coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "bcc") then
       ibrav = 2
       boa=1d0
       coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "hcp" .or. brvtyp == "hex") then
       ibrav = 3
       boa=1d0
       if(abs(coa) < eps) coa=sqrt(8d0/3d0)
       alpha=90d0
       beta=90d0
       gamma=120d0
       call calculate_x2c()
    elseif (brvtyp == "sc ") then
       ibrav = 4
       boa=1d0
       coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "bct") then
       ibrav = 5
       boa=1d0
       if(abs(coa) < eps) coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "st ") then
       ibrav = 6
       boa=1d0
       if(abs(coa) < eps) coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "fco") then
       ibrav = 7
       alpha=90d0
       beta=90d0
       gamma=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       call calculate_x2c()
    elseif (brvtyp == "bco") then
       ibrav = 8
       alpha=90d0
       beta=90d0
       gamma=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       call calculate_x2c()
    elseif (brvtyp == "bso") then
       ibrav = 9
       alpha=90d0
       beta=90d0
       gamma=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       call calculate_x2c()
    elseif (brvtyp == "so ") then
       ibrav = 10
       alpha=90d0
       beta=90d0
       gamma=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       call calculate_x2c()
    elseif (brvtyp == "bsm") then
       ibrav = 11
       alpha=90d0
       gamma=90d0
       if(abs(beta) < eps) beta=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       seed%m_x2c(1,1)=5d-1
       seed%m_x2c(2,1)=-5d-1*boa
       seed%m_x2c(3,1)=0d0
       seed%m_x2c(1,2)=5d-1
       seed%m_x2c(2,2)=5d-1*boa
       seed%m_x2c(3,2)=0d0
       seed%m_x2c(1,3)=cos(pi/180d0*beta)*coa
       seed%m_x2c(2,3)=0d0
       seed%m_x2c(3,3)=sqrt(1d0-(cos(pi/180d0*beta))**2)*coa
    elseif (brvtyp == "sm ") then
       ibrav = 12
       alpha=90d0
       gamma=90d0
       if(abs(beta) < eps) beta=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       seed%m_x2c(1,1)=1d0
       seed%m_x2c(2,1)=0d0
       seed%m_x2c(3,1)=0d0
       seed%m_x2c(1,2)=0d0
       seed%m_x2c(2,2)=boa
       seed%m_x2c(3,2)=0d0
       seed%m_x2c(1,3)=cos(pi/180d0*beta)*coa
       seed%m_x2c(2,3)=0d0
       seed%m_x2c(3,3)=sqrt(1d0-(cos(pi/180d0*beta))**2)*coa
    elseif (brvtyp == "trc") then
       ibrav = 13
       if(abs(alpha) < eps) alpha=90d0
       if(abs(beta) < eps) beta=90d0
       if(abs(gamma) < eps) gamma=90d0
       if(abs(coa) < eps) coa=1d0
       if(abs(boa) < eps) boa=1d0
       cw=1d0-(cos(pi/180d0*alpha))**2-(cos(pi/180d0*beta))**2-&
          (cos(pi/180d0*gamma))**2+2*cos(pi/180d0*alpha)*&
          cos(pi/180d0*beta)*cos(pi/180d0*gamma)
       seed%m_x2c(1,1)=1d0
       seed%m_x2c(2,1)=0d0
       seed%m_x2c(3,1)=0d0
       seed%m_x2c(1,2)=cos(pi/180d0*gamma)*boa
       seed%m_x2c(2,2)=sin(pi/180d0*gamma)*boa
       seed%m_x2c(3,2)=0d0
       seed%m_x2c(1,3)=cos(pi/180d0*beta)*coa
       seed%m_x2c(2,3)=(cos(pi/180d0*alpha)-cos(pi/180d0*beta)*&
          cos(pi/180d0*gamma))/sqrt(1d0-(cos(pi/180d0*gamma))**2)*coa
       seed%m_x2c(3,3)=sqrt(cw)/sqrt(1d0-(cos(pi/180d0*gamma))**2)*coa
    elseif (brvtyp == "rhb" .or. brvtyp == "trg") then
       ibrav = 14
       boa=1d0
       coa=1d0
       if(abs(alpha) < eps) alpha=60d0
       beta=alpha
       gamma=alpha
       theta=2d0*pi*alpha/360d0
       z=sqrt((5d-1+cos(theta))/(1d0-cos(theta)))
       snrm=sqrt(1d0+z**2)
       ! obverse setting (used in versions after 2015)
       seed%m_x2c(1,1)=(sqrt(3d0)/2d0)/snrm
       seed%m_x2c(2,1)=-(1d0/2d0)/snrm
       seed%m_x2c(3,1)=z/snrm
       seed%m_x2c(1,2)=0d0
       seed%m_x2c(2,2)=1d0/snrm
       seed%m_x2c(3,2)=z/snrm
       seed%m_x2c(1,3)=-(sqrt(3d0)/2d0)/snrm
       seed%m_x2c(2,3)=-(1d0/2d0)/snrm
       seed%m_x2c(3,3)=z/snrm
    elseif (brvtyp == "fct") then
       ibrav = 15
       boa=1d0
       if(abs(coa) < eps) coa=1d0
       alpha=90d0
       beta=90d0
       gamma=90d0
       call calculate_x2c()
    elseif (brvtyp == "aux" .or. brvtyp == "prv") then
       errmsg = "Bravais types aux and prv not supported"
       return
    end if
    seed%m_x2c = seed%m_x2c * a

    ! continue reading the input file
    ! edelt ... pmix
    do i = 1, 10
       call read_next_token()
    end do

    ! ntyp
    call read_next_token()
    ok = isinteger(seed%nspc,token)
    if (.not.ok) goto 999
    allocate(seed%spc(seed%nspc))

    ! read the atomic species
    do i = 1, seed%nspc
       ! type, number of components
       call read_next_token()
       seed%spc(i)%name = token
       call read_next_token()
       ok = isinteger(ncmp,token)
       if (.not.ok) goto 999
       do j = 1, 3 ! rmt .. mxl
          call read_next_token()
       end do

       ! use the first component as the representative
       call read_next_token() ! anclr
       ok = isinteger(seed%spc(i)%z,token)
       if (.not.ok) goto 999
       call read_next_token() ! conc

       ! skip the rest of the components
       do j = 1, ncmp-1
          call read_next_token() ! conc
          call read_next_token() ! conc
       end do

       ! wrap up
    end do

    ! natm, number of atoms, alllocate atomic info
    call read_next_token()
    ok = isinteger(seed%nat,token)
    if (.not.ok) goto 999
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
    seed%x = 0d0

    ! atmicx: atomic coordinates and types
    do i = 1, seed%nat
       do j = 1, 3
          call read_next_token()
          token = trim(adjustl(lower(token)))
          ll = len(token)
          if (token(ll:ll)=="a" .or. token(ll:ll)=="b" .or. token(ll:ll)=="c") then
             if (token(ll:ll) == "a") then
                id = 1
             elseif (token(ll:ll) == "b") then
                id = 2
             else
                id = 3
             end if
             if (ll == 1) then
                xc = 1d0
             else
                ok = eval_next(xc,token(1:ll-1))
                if (.not.ok) goto 999
             end if
             do k = 1, 3
                seed%x(k,i) = seed%x(k,i) + xc * seed%m_x2c(k,id)
             end do
          elseif (token(ll:ll)=="x" .or. token(ll:ll)=="y" .or. token(ll:ll)=="z") then
             if (ll == 1) then
                xc = 1d0
             else
                ok = eval_next(xc,token(1:ll-1))
                if (.not.ok) goto 999
             end if
             if (token(ll:ll) == "x") then
                seed%x(id,i) = seed%x(id,i) + xc * a
             elseif (token(ll:ll) == "y") then
                seed%x(id,i) = seed%x(id,i) + xc * a * boa
             else
                seed%x(id,i) = seed%x(id,i) + xc * a * coa
             end if
          else
             ok = eval_next(xc,token)
             seed%x(j,i) = seed%x(j,i) + xc * a
          end if
       end do

       ! read the atomic name and the species
       call read_next_token()
       seed%is(i) = 0
       do j = 1, seed%nspc
          if (seed%spc(j)%name == token) then
             seed%is(i) = j
             seed%atname(i) = token
             exit
          end if
       end do
       if (seed%is(i) == 0) goto 999
    end do

    ! convert to fractional coordinates
    c2x = seed%m_x2c
    call matinv(c2x,3,ier)
    if (ier /= 0) goto 999
    do i = 1, seed%nat
       seed%x(:,i) = matmul(c2x,seed%x(:,i))
    end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    ! no symmetry
    seed%havesym = 0
    seed%checkrepeats = .false.
    seed%findsym = -1

    ! rest of the seed information
    seed%isused = .true.
    seed%ismolecule = mol
    seed%cubic = .false.
    seed%border = 0d0
    seed%havex0 = .false.
    seed%molx0 = 0d0
    return

999 continue
    call fclose(lu)
    seed%isused = .false.

  contains
    subroutine read_next_token()
      logical :: ok, newl
      integer :: i, wp, len0

      token = ""
      newl = (lp > len_trim(line))
      if (.not.newl) newl = len_trim(line(lp:)) == 0

      if (newl) then
         do while(.true.)
            ! read new line
            ok = getline_raw(lu,line,.false.)
            if (.not.ok) return
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "c" .or. line(1:1) == "C") cycle
            line = adjustl(line)
            exit
         end do
         ! reset the line pointer
         lp = 1
      end if

      ! read the next token, note there are commas
      len0 = len(line)
      token = ""

      !! advance to next non-blank character; return a blank token if nothing left
      i = lp
      do while (line(i:i) == " ")
         i = i+1
         if (i > len0) then
            lp = len0
            return
         end if
      enddo

      !! if it is a comma, return a blank token
      if (line(i:i) == ",") then
         lp = min(i + 1,len0)
         return
      end if

      !! read the token; stop at blank or comma
      wp = i
      do while (i<len0 .and. line(i:i) /= " " .and. line(i:i) /= ",")
         i = i+1
      enddo
      token = line(wp:i-1)
      lp = i

      !! if stopped at blank, advance to the next non-blank; skip the comma
      do while (lp<len0 .and. line(lp:lp) == " ")
         lp = lp + 1
      end do
      if (line(lp:lp) == ",") lp = lp + 1

    end subroutine read_next_token

    subroutine calculate_x2c()

      real*8, parameter :: c0 = 0d0
      real*8, parameter :: c1 = 0.5d0
      real*8, parameter :: c2 = -0.5d0
      real*8, parameter :: c3 = 1d0
      real*8, parameter :: c5 = 8.660254037844386d-1
      real*8, parameter :: c6 = -c5
      real*8, parameter :: lats(3,3,15) = reshape((/&
         c0,c1,c1, c1,c0,c1, c1,c1,c0,  c2,c1,c1, c1,c2,c1, c1,c1,c2,&
         c1,c6,c0, c1,c5,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,&
         c2,c1,c1, c1,c2,c1, c1,c1,c2,  c3,c0,c0, c0,c3,c0, c0,c0,c3,&
         c0,c1,c1, c1,c0,c1, c1,c1,c0,  c2,c1,c1, c1,c2,c1, c1,c1,c2,&
         c1,c2,c0, c1,c1,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,&
         c1,c2,c0, c1,c1,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,&
         c3,c0,c0, c0,c3,c0, c0,c0,c3,  c0,c0,c0, c0,c0,c0, c0,c0,c0,&
         c0,c1,c1, c1,c0,c1, c1,c1,c0/),shape(lats))

      integer :: i

      do i = 1, 3
         seed%m_x2c(1,i) = lats(1,i,ibrav)
         seed%m_x2c(2,i) = lats(2,i,ibrav) * boa
         seed%m_x2c(3,i) = lats(3,i,ibrav) * coa
      end do

    end subroutine calculate_x2c

  end subroutine read_akaikkr

  !> Read the structure from an xband sys input file.
  module subroutine read_xband(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, equal, getline_raw, isinteger, isreal, getword
    use param, only: isformat_r_xband
    class(crystalseed), intent(inout)  :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lp, lu, i, iaux, id, idum
    real*8 :: a, rdum
    logical :: hada, hadx2c, hadnat, hadcoords, hadnspc, hadspc
    logical :: ok
    character(len=:), allocatable :: line, word, word2, word3, word4
    type(rawseed) :: rec

    ! file and seed name
    call seed%end()
    errmsg = "Error reading file: " // trim(file)

    ! initialize
    a = 1d0
    hada = .false.
    hadx2c = .false.
    hadnat = .false.
    hadcoords = .false.
    hadnspc = .false.
    hadspc = .false.

    ! open the file
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! read all lines in one go
    do while (getline_raw(lu,line,.false.))
       lp = 1
       word = getword(line,lp)
       word2 = getword(line,lp)
       word3 = getword(line,lp)
       word4 = getword(line,lp)

       if (equal(word,"dimension")) then
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          if (.not.equal(line,"3D")) then
             errmsg = "cannot read non-3D structures from sys files"
             goto 999
          end if
          cycle
       elseif (equal(word,"lattice").and.equal(word2,"parameter")) then
          read (lu,*,end=999,err=999) a
          hada = .true.
          cycle
       elseif (equal(word,"primitive")) then
          do i = 1, 3
             read (lu,*,end=999,err=999) rec%rv(:,i)
          end do
          hadx2c = .true.
          cycle
       elseif (equal(word,"number").and.equal(word4,"NQ")) then
          read (lu,*,end=999,err=999) rec%nat
          hadnat = .true.
          allocate(rec%x(3,rec%nat),rec%is(rec%nat))
          cycle
       elseif (equal(word,"IQ")) then
          if (.not.hadnat) goto 999
          do i = 1, rec%nat
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             read (line,*,end=999,err=999) iaux, iaux, rec%x(:,i), rdum, idum, idum, rec%is(i)
          end do
          hadcoords = .true.
          cycle
       elseif (equal(word,"number").and.equal(word3,"atom")) then
          read (lu,*,end=999,err=999) rec%nspc
          allocate(rec%spc(rec%nspc))
          hadnspc = .true.
          cycle
       elseif (equal(word,"IT")) then
          if (.not.hadnspc) goto 999
          do i = 1, rec%nspc
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999

             lp = 1
             word = getword(line,lp)
             ok = isinteger(id,word)
             word = getword(line,lp)
             ok = ok .and. isinteger(rec%spc(id)%z,word)
             rec%spc(id)%name = getword(line,lp)
             if (.not.ok) goto 999
          end do
          hadspc = .true.
          cycle
       end if
    end do
    if (.not.hada.or..not.hadx2c.or..not.hadnat.or..not.hadcoords.or.&
       .not.hadnspc.or..not.hadspc) goto 999

    ! multiply by a (Cartesian coordinates, atomic units)
    rec%rv = rec%rv * a
    rec%x = rec%x * a
    rec%spcmode = spc_given
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_bohr
    rec%iscart = .true.

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_xband,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_xband

  !> Read the structure from a fort.21 from neighcrys
  module subroutine read_f21(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, nameguess
    use param, only: bohrtoa, maxzat0, mlen, isformat_r_f21
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
    seed%isformat = isformat_r_f21

    ! read the lattice constants and the atomic coordinates
    isused = 0
    seed%useabr = 1
    seed%nat = 0
    seed%nspc = 0
    seed%neqv = 0
    seed%ncv = 0
    allocate(seed%spc(10),seed%is(10),seed%atname(10),seed%rotm(3,4,48),seed%cen(3,4))
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
             if (seed%nat > size(seed%is,1)) then
                call realloc(seed%is,2*seed%nat)
                call realloc(seed%atname,2*seed%nat)
             end if
             seed%is(seed%nat) = isused(iz)
             seed%atname(seed%nat) = str(1:10)
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
    call realloc(seed%atname,seed%nat)
    call realloc(seed%spc,seed%nspc)
    call realloc(seed%rotm,3,4,seed%neqv)
    call realloc(seed%cen,3,seed%ncv)

    ! close the file and clean up
    call fclose(lu)

    ! have symmetry, but recalculate (a molecule does not apply crystallographic
    ! symmetry)
    if (.not.mol) then
       seed%havesym = 1
       seed%neqlist = .true.
    else
       seed%havesym = 0
       seed%neqlist = .false.
    end if
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
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, nameguess, getline_raw
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: isformat_r_cube
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    integer :: i, j, nstep(3), nn, iz, it, ier
    real*8 :: x0(3), rmat(3,3), rdum, rx(3)
    logical :: ismo, ok, use0
    character(len=:), allocatable :: line
    real*8, allocatable :: rxc(:,:)

    real*8, parameter :: eps = 1d-5

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
    seed%isformat = isformat_r_cube

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

    ! Atomic positions.
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
    allocate(seed%spc(2),rxc(3,seed%nat))
    nn = seed%nat
    seed%nat = 0
    do i = 1, nn
       read (lu,*,err=999,end=999) iz, rdum, rx
       if (iz <= 0) cycle

       seed%nat = seed%nat + 1
       rxc(:,seed%nat) = rx
       seed%x(:,seed%nat) = matmul(rx - x0,rmat)
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
          seed%spc(seed%nspc)%name = nameguess(iz,.true.)
          it = seed%nspc
       end if
       seed%is(seed%nat) = it
       seed%atname(seed%nat) = seed%spc(it)%name
    end do
    if (seed%nat /= nn) then
       call realloc(seed%x,3,seed%nat)
       call realloc(seed%is,seed%nat)
       call realloc(seed%atname,seed%nat)
    end if
    call realloc(seed%spc,seed%nspc)

    errmsg = ""
999 continue
    call fclose(lu)

    ! if this is a molecule, and atoms are outside the cell, use our
    ! own cell instead of the one provided by the cube
    use0 = .false.
    if (mol) then
       do i = 1, seed%nat
          if (any(seed%x(:,i) < -eps) .or. any(seed%x(:,i) > 1d0+eps)) then
             use0 = .true.
             exit
          end if
          seed%x(:,i) = min(max(seed%x(:,i),0d0),1d0)
       end do
    end if

    if (mol) then
       ! treat this as a molecule
       if (use0) then
          seed%useabr = 0
          seed%border = rborder_def
          seed%x = rxc
       else
          seed%useabr = 2
          seed%border = 0d0
       end if
    else
       ! treat it as a crystal
       seed%useabr = 2
       seed%border = 0d0
    end if

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

  end subroutine read_cube

  !> Read the structure from a binary cube file
  module subroutine read_bincube(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, nameguess, getline_raw
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: isformat_r_bincube
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
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
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
             seed%spc(seed%nspc)%name = nameguess(iz,.true.)
             it = seed%nspc
          end if
          seed%is(seed%nat) = it
          seed%atname(seed%nat) = seed%spc(it)%name
       endif
    end do
    if (seed%nat /= nn) then
       call realloc(seed%x,3,seed%nat)
       call realloc(seed%is,seed%nat)
       call realloc(seed%atname,seed%nat)
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
    seed%isformat = isformat_r_bincube

  end subroutine read_bincube

  !> Read the crystal structure from a WIEN2k STRUCT file.
  !> Code adapted from the WIEN2k distribution.
  module subroutine read_wien(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, zatguess, fclose, equal, equali
    use types, only: realloc
    use param, only: pi, isformat_r_struct
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
    seed%isformat = isformat_r_struct

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
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
    iat = 0
    DO JATOM=1,seed%nat
       iat0 = iat
       iat = iat + 1
       if (iat > size(seed%is)) then
          call realloc(seed%x,3,2*iat)
          call realloc(seed%is,2*iat)
          call realloc(seed%atname,2*iat)
       end if
       READ(lut,1012,err=999,end=999) iatnr,seed%x(:,iat),MULTW

       istart = iat
       if (readall) then
          DO MU=1,MULTW-1
             iat = iat + 1
             if (iat > size(seed%is)) then
                call realloc(seed%x,3,2*iat)
                call realloc(seed%is,2*iat)
                call realloc(seed%atname,2*iat)
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
       seed%atname(jatom) = aname
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
    call realloc(seed%atname,iat)
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

    ! symmetry (a molecule does not apply crystallographic symmetry)
    if (seed%neqv > 0) then
       seed%findsym = 0
       if (.not.mol) then
          seed%havesym = 1
          seed%neqlist = .true.
       else
          seed%havesym = 0
          seed%neqlist = .false.
       end if
    else
       seed%findsym = -1
       seed%havesym = 0
       seed%neqlist = .false.
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
       getword, zatguess, string, isinteger, fclose
    use tools_math, only: det3sym, matinv
    use param, only: bohrtoa, isformat_r_vasp
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
    allocate(seed%x(3,seed%nat),seed%atname(seed%nat))
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
       if (hastypes) then
          seed%atname(i) = seed%spc(seed%is(i))%name
       else
          seed%atname(i) = ""
       end if
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
    seed%isformat = isformat_r_vasp

  end subroutine read_vasp

  !> Read the structure from an abinit DEN file (and similar files: ELF, LDEN, etc.)
  module subroutine read_abinit(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, nameguess, fclose
    use abinit_private, only: hdr_type, hdr_io
    use param, only: isformat_r_abinit
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, fform0
    type(hdr_type) :: hdr
    integer :: i, iz
    type(rawseed) :: rec

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

    ! cell parameters (bohr)
    rec%rv = hdr%rprimd(:,:)
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_bohr

    ! types
    rec%spcmode = spc_given
    rec%nspc = hdr%ntypat
    allocate(rec%spc(rec%nspc))
    do i = 1, rec%nspc
       iz = nint(hdr%znucltypat(i))
       rec%spc(i)%z = iz
       rec%spc(i)%name = nameguess(iz,.true.)
    end do

    ! atoms
    rec%nat = hdr%natom
    allocate(rec%x(3,rec%nat),rec%is(rec%nat))
    do i = 1, rec%nat
       rec%x(:,i) = hdr%xred(:,i)
       rec%is(i) = hdr%typat(i)
    end do

    ! abinit has symmetry in hdr%nsym/hdr%symrel, but there is no
    ! distinction between pure centering and rotation operations, and
    ! the user may not want any symmetry - let critic2 guess.
    call rawseed_to_seed(rec,seed,mol,file,isformat_r_abinit,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_abinit

  ! The following code has been adapted from the elk distribution, version 1.3.2
  ! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU General Public License.
  module subroutine read_elk(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, equal, getword,&
       zatguess, fclose, string
    use types, only: realloc
    use param, only: isformat_r_elk
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< input filename
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, atname
    integer :: lu, i, zat, j, lp, idx
    integer :: natoms
    logical :: ok
    real*8 :: x(3)
    type(rawseed) :: rec

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

    read(lu,'(3G18.10)',err=999,end=999) rec%rv(:,1)
    read(lu,'(3G18.10)',err=999,end=999) rec%rv(:,2)
    read(lu,'(3G18.10)',err=999,end=999) rec%rv(:,3)
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_bohr

    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    ok = getline_raw(lu,line,.true.)
    if (.not.ok) goto 999
    if (equal(line,'molecule')) then
       errmsg = "Isolated molecules not supported."
       goto 999
    end if

    ! species and atoms
    rec%spcmode = spc_given
    read(lu,'(I4)',err=999,end=999) rec%nspc
    allocate(rec%spc(rec%nspc),rec%is(10))
    do i = 1, rec%nspc
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
       rec%spc(i)%z = zat

       idx = index(atname,".in",.true.)
       if (idx > 1) &
          atname = trim(atname(1:idx-1))
       rec%spc(i)%name = atname

       read(lu,*,err=999,end=999) natoms
       do j = 1, natoms
          read(lu,*,err=999,end=999) x
          call rawseed_add_atom(rec,zat,x,atname)
          if (rec%nat > size(rec%is,1)) call realloc(rec%is,2*rec%nat)
          rec%is(rec%nat) = i
       end do
    end do

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_elk,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_elk

  !> Read the structure from a molecule file
  module subroutine read_mol(seed,file,fmt,rborder,docube,errmsg,alsovib,ti)
    use wfn_private, only: wfn_read_xyz_geometry, wfn_read_wfn_geometry, &
       wfn_read_wfx_geometry, wfn_read_fchk_geometry, wfn_read_molden_geometry,&
       wfn_read_log_geometry, wfn_read_dat_geometry, wfn_read_pgout_geometry,&
       wfn_read_orca_geometry, wfn_read_gjf_geometry
    use param, only: isformat_r_xyz, isformat_r_wfn, isformat_r_wfx,&
       isformat_r_fchk, isformat_r_molden, isformat_r_gaussian, isformat_r_dat,&
       isformat_r_pgout, isformat_r_orca, isformat_r_gjf, isformat_r_zmat
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: fmt !< wfn/wfx/xyz
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(out), optional :: alsovib
    type(thread_info), intent(in), optional :: ti

    integer, allocatable :: z(:)
    real*8, allocatable :: x(:,:)
    character*10, allocatable :: atname(:)
    integer :: nat
    type(rawseed) :: rec

    if (present(alsovib)) alsovib = .false.
    call seed%end()
    errmsg = ""
    nat = 0
    if (fmt == isformat_r_xyz) then
       ! xyz
       call wfn_read_xyz_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_gjf) then
       ! xyz
       call wfn_read_gjf_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_wfn) then
       ! wfn
       call wfn_read_wfn_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_wfx) then
       ! wfx
       call wfn_read_wfx_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_fchk) then
       ! fchk
       call wfn_read_fchk_geometry(file,nat,x,z,atname,errmsg,alsovib=alsovib,ti=ti)
    elseif (fmt == isformat_r_molden) then
       ! molden (psi4)
       call wfn_read_molden_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_gaussian) then
       ! Gaussian output file
       call wfn_read_log_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_dat) then
       ! psi4 output file
       call wfn_read_dat_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_pgout) then
       ! postg output file
       call wfn_read_pgout_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_orca) then
       ! orca output file
       call wfn_read_orca_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    elseif (fmt == isformat_r_zmat) then
       call read_zmat_geometry(file,nat,x,z,atname,errmsg,ti=ti)
    end if
    if (len_trim(errmsg) > 0) goto 999
    if (nat == 0) then
       errmsg = "No atomic species found."
       goto 999
    end if

    ! molecule in atomic units; species by label (case-insensitive)
    rec%spcmode = spc_label_nocase
    rec%ndim = 0
    rec%lunit = lunit_bohr
    rec%iscart = .true.
    rec%ismol = .true.
    rec%border = rborder
    rec%cubic = docube
    rec%nat = nat
    call move_alloc(x,rec%x)
    call move_alloc(z,rec%z)
    call move_alloc(atname,rec%atname)

    call rawseed_to_seed(rec,seed,.true.,file,fmt,errmsg)
999 continue

  end subroutine read_mol

  !> Read a pdb file and return a molecule or a crystal, depending
  !> on whether SCALEn is present or not.
  module subroutine read_pdb(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, zatguess, nameguess, equalsub
    use tools_math, only: matinv
    use param, only: bohrtoa, isformat_r_pdb
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: i, ier, ll, iz
    integer :: lu
    character(len=:), allocatable :: line
    logical :: readscale(3)
    real*8 :: r(3,3), x0(3), x(3)
    type(rawseed) :: rec

    ! initialize
    call seed%end()
    errmsg = ""
    readscale = .false.

    ! open the file
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    rec%spcmode = spc_z
    main: do while (getline_raw(lu,line))
       ll = len(line)
       if ((equalsub(line,1,4,"ATOM") .or. equalsub(line,1,6,"HETATM")) .and. ll >= 78) then
          iz = zatguess(line(77:78))
          if (iz <= 0 .or. iz > 118) cycle
          read (line(31:38),*,err=999,end=999) x(1)
          read (line(39:46),*,err=999,end=999) x(2)
          read (line(47:54),*,err=999,end=999) x(3)
          call rawseed_add_atom(rec,iz,x,adjustl(line(13:16)),spcname=nameguess(iz,.true.))
       elseif (line(1:5) == "SCALE") then
          if (line(6:6) == "1") then
             readscale(1) = .true.
             read(line(11:20),*) r(1,1)
             read(line(21:30),*) r(1,2)
             read(line(31:40),*) r(1,3)
             read(line(46:55),*) x0(1)
          elseif (line(6:6) == "2") then
             readscale(2) = .true.
             read(line(11:20),*) r(2,1)
             read(line(21:30),*) r(2,2)
             read(line(31:40),*) r(2,3)
             read(line(46:55),*) x0(2)
          elseif (line(6:6) == "3") then
             readscale(3) = .true.
             read(line(11:20),*) r(3,1)
             read(line(21:30),*) r(3,2)
             read(line(31:40),*) r(3,3)
             read(line(46:55),*) x0(3)
          end if
       end if
    end do main
    if (rec%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    endif

    if (.not.mol.and.all(readscale)) then
       ! read it as a crystal: the SCALE matrix gives the fractional coordinates
       rec%ndim = 3
       rec%cellmode = 2
       rec%lunit = lunit_bohr
       rec%rv = r * bohrtoa
       call matinv(rec%rv,3,ier)
       if (ier /= 0) then
          errmsg = "error inverting lattice vector matrix"
          goto 999
       end if
       do i = 1, rec%nat
          rec%x(:,i) = matmul(r,rec%x(:,i)) + x0
       end do
    else
       ! read it as a molecule (angstrom)
       rec%ndim = 0
       rec%lunit = lunit_ang
       rec%iscart = .true.
    end if

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_pdb,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_pdb

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

  !> Read a structure from a GULP input file (gin/grs). If istruct is
  !> zero, read the last configuration; otherwise, read configuration
  !> number istruct.
  module subroutine read_gulpin(seed,file,mol,istruct,errmsg,ti)
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    integer, intent(in) :: istruct
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: nseed
    type(crystalseed), allocatable :: seedaux(:)

    call seed%end()
    call read_all_gulpin(nseed,seedaux,file,mol,istruct,errmsg,ti=ti)
    if (allocated(seedaux) .and. nseed >= 1) seed = seedaux(1)

  end subroutine read_gulpin

  !> Read a structure from a GULP output file (gout/got). If istruct
  !> is zero, read the last geometry; otherwise, read geometry number
  !> istruct.
  module subroutine read_gulpout(seed,file,mol,istruct,errmsg,ti)
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    integer, intent(in) :: istruct
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: nseed
    type(crystalseed), allocatable :: seedaux(:)

    call seed%end()
    call read_all_gulpout(nseed,seedaux,file,mol,istruct,errmsg,ti=ti)
    if (allocated(seedaux) .and. nseed >= 1) seed = seedaux(1)

  end subroutine read_gulpout

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
    use param, only: bohrtoa, isformat_r_qein
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
    logical :: ok

    ! Namelist variables used below. The rest of the QE input is ignored.
    ! A real Fortran namelist read is not used because it aborts on any
    ! variable the declaration does not know about, and every QE version
    ! (and every locally patched QE) has its own extra keywords. Instead,
    ! the namelist groups are scanned by qe_read_namelists.
    character(len=80) :: calculation
    integer :: ibrav, nat, ntyp, origin_choice, space_group
    real*8 :: celldm(6)
    logical :: rhombohedral, uniqueb

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

    ! open
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    r = 0d0
    calculation = ""

    ! read the namelists (tolerant scan; see qe_read_namelists)
    call qe_read_namelists(lu,calculation,ibrav,celldm,nat,ntyp,space_group,&
       uniqueb,rhombohedral,origin_choice,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! allocate space for atoms
    seed%nat = nat
    seed%nspc = ntyp
    allocate(seed%x(3,nat),seed%is(nat),seed%spc(ntyp),seed%atname(nat))

    ! read the cards
    iunit = icrystal
    cunit = ialat
    do while (getline_raw(lu,line))
       line = lower(line)
       lp = 1
       word = getword(line,lp)
       if (equal(word,'atomic_species')) then
          errmsg = "Error reading atomic species."
          do i = 1, ntyp
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read (line,*,iostat=ios) seed%spc(i)%name
             if (ios/=0 .or. len(seed%spc(i)%name) == 0) goto 999
             seed%spc(i)%z = zatguess(seed%spc(i)%name)
             if (seed%spc(i)%z < 0) goto 999
          end do

       else if (equal(word,'atomic_positions')) then
          errmsg = "Error reading atomic positions."
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
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999

             read (line,*,iostat=ios) atm, seed%x(:,i)
             if (ios/=0 .or. len(atm) == 0) goto 999
             seed%is(i) = 0
             do j = 1, seed%nspc
                if (equal(seed%spc(j)%name,atm)) then
                   seed%is(i) = j
                   seed%atname(i) = atm
                   exit
                end if
             end do
             if (seed%is(i) == 0) goto 999
          end do
       elseif (equal(word,'cell_parameters')) then
          errmsg = "Error reading CELL_PARAMETERS."
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
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             read (line,*,iostat=ios) (r(i,j),j=1,3)
             if (ios/=0) goto 999
          end do
       endif
    end do
    errmsg = ""

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
       if (allocated(seed%atname)) deallocate(seed%atname)
       allocate(seed%x(3,nattot),seed%is(nattot),seed%atname(nattot))
       do i = 1, nattot
          seed%x(:,i) = tautot(:,i)
          seed%is(i) = ityptot(i)
          seed%atname(i) = seed%spc(ityptot(i))%name
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
       ! not wrapped into the main cell: the coordinates as given are
       ! kept until the crystal is built (struct_new records the shift)
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
    seed%isformat = isformat_r_qein

  end subroutine read_qein

  !> Read the structure from a crystal output
  module subroutine read_crystalout(seed,file,mol,errmsg,ti)
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    call seed%end()
    call read_all_crystalout(file,mol,errmsg,seed0=seed,ti=ti)

  end subroutine read_crystalout

  !> Read the structure from an FPLO output
  module subroutine read_fploout(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, isinteger, isreal,&
       zatguess, nameguess, fclose, equalsub
    use param, only: maxzat, isformat_r_fploout
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, ll, nat
    character(len=:), allocatable :: line
    integer :: iz, lp, idum1, idum2, idum3
    logical :: ok, iscell, isatoms
    character*(10) :: ats
    real*8 :: r(3,3), x(3)
    type(rawseed) :: rec

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
       elseif (equalsub(line,18,27,"Atom sites")) then
          ! number of sites
          ok = getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          ok = ok .and. isinteger(nat,line(18:))
          do i = 1, 3
             ok = ok .and. getline_raw(lu,line)
          end do
          if (.not.ok) then
             errmsg = "Error reading number of sites."
             goto 999
          end if

          ! sites (the last block wins)
          rec = rawseed()
          rec%spcmode = spc_z
          do i = 1, nat
             ok = getline_raw(lu,line)
             read (line,*) idum1, ats, idum2, idum3, x
             iz = zatguess(ats)
             if (iz < 1 .or. iz > maxzat) then
                errmsg = "Unknown atomic species: " // ats // "."
                goto 999
             end if
             call rawseed_add_atom(rec,iz,x,ats,spcname=nameguess(iz,.true.))
          end do

          ! all done
          isatoms = (rec%nat > 0)
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

    ! cell (atomic units, Cartesian coordinates wrapped into the cell)
    rec%ndim = 3
    rec%cellmode = 2
    rec%rv = transpose(r)
    rec%lunit = lunit_bohr
    rec%iscart = .true.
    rec%wrap = .true.

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_fploout,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_fploout

  !> Read the structure from a siesta STRUCT_IN/OUT input
  module subroutine read_siesta(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, nameguess, fclose
    use types, only: realloc
    use param, only: isformat_r_siesta
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, idum
    type(rawseed) :: rec

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

    ! the atoms
    rec%nspc = 0
    read (lu,*,err=999,end=999) rec%nat
    allocate(rec%x(3,rec%nat),rec%is(rec%nat),rec%spc(2))
    do i = 1, rec%nat
       read (lu,*,err=999,end=999) rec%is(i), idum, rec%x(:,i)
       if (rec%is(i) > size(rec%spc,1)) &
          call realloc(rec%spc,2*rec%is(i))
       rec%nspc = max(rec%nspc,rec%is(i))
       rec%spc(rec%is(i))%z = idum
       rec%spc(rec%is(i))%name = nameguess(idum,.true.)
    end do
    call realloc(rec%spc,rec%nspc)

    ! fill the cell metrics
    rec%spcmode = spc_given
    rec%ndim = 3
    rec%cellmode = 2
    rec%rv = transpose(r)
    rec%lunit = lunit_ang

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_siesta,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_siesta

  !> Read the structure from a CASTEP cell file
  module subroutine read_castep_cell(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, lgetword, getline_raw,&
       equal, getword, isinteger, zatguess, isreal, lower
    use param, only: bohrtoa, bohrtom, bohrtocm, bohrtonm, isformat_r_castepcell
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word
    integer :: lu, lp, idum, iz
    logical :: ok, iscart
    real*8 :: rconv, x(3)
    type(rawseed) :: rec

    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! atomic units (unit conversion applied while reading); species by label
    iscart = .false.
    rec%spcmode = spc_label_nocase
    rec%ndim = 3
    rec%lunit = lunit_bohr
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
             read(line,*,err=999,end=999) rec%aa
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) rec%bb

             ! conversion
             rec%cellmode = 1
             rec%aa = rec%aa * rconv
          elseif (equal(word,"lattice_cart")) then
             !! lattice vectors
             ! optional unit
             if (.not.get_optional_unit(rconv)) goto 999

             ! read the rest
             read(line,*,err=999,end=999) rec%rv(:,1)
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) rec%rv(:,2)
             if (.not.get_next_line()) goto 999
             read(line,*,err=999,end=999) rec%rv(:,3)
             if (.not.get_next_line()) goto 999

             ! conversion
             rec%cellmode = 2
             rec%rv = rec%rv * rconv

          elseif (equal(word,"positions_frac").or.equal(word,"positions_abs")) then
             iscart = equal(word,"positions_abs")
             if (iscart) then
                if (.not.get_optional_unit(rconv)) goto 999
             else
                if (.not.get_next_line()) goto 999
             end if

             !! positions in fractional/Cartesian coordinates (a new block replaces the atoms)
             rec%nat = 0
             if (allocated(rec%x)) deallocate(rec%x,rec%z,rec%atname)
             ok = .true.
             do while(ok)
                if (line(1:1) == "%") exit

                ! read the atomic symbol or number
                lp = 1
                word = getword(line,lp)
                if (isinteger(idum,word)) then
                   iz = idum
                else
                   iz = zatguess(word)
                   if (iz < 0) then
                      errmsg = "Unknown atomic symbol: " // word
                      goto 999
                   end if
                end if

                ! read the atomic coordinates
                ok = isreal(x(1),line,lp)
                ok = ok .and. isreal(x(2),line,lp)
                ok = ok .and. isreal(x(3),line,lp)
                if (.not.ok) goto 999
                if (iscart) x = x * rconv
                call rawseed_add_atom(rec,iz,x,word)

                ok = get_next_line()
             end do
             if (rec%nat == 0) then
                errmsg = "Error reading atoms"
                goto 999
             end if
          end if
       end if
    end do
    ! consistency checks
    if (rec%cellmode == 0) then
       errmsg = "No lattice block found"
       goto 999
    end if
    if (rec%nat == 0) then
       errmsg = "No atoms found"
       goto 999
    end if
    if (iscart .and. rec%cellmode /= 2) then
       errmsg = "Atomic Cartesian (absolute) coordinates require lattice_cart"
       goto 999
    end if
    rec%iscart = iscart

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_castepcell,errmsg)
999 continue
    call fclose(lu)

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
    use param, only: isformat_r_castepgeom
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
       if (ll > 4) then
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
       end if
    end do
    if (nat == 0) goto 999
    if (nspc == 0) goto 999

    ! second pass, actual read
    seed%useabr = 2
    seed%nat = nat
    seed%nspc = nspc
    allocate(seed%x(3,nat),seed%is(nat),seed%spc(nspc),usespc(nspc),seed%atname(nat))
    usespc = .false.
    rewind(lu)
    do i = 1, nlast-1
       read (lu,*,err=999,end=999)
    end do
    nat = 0
    do while(getline_raw(lu,line))
       line = trim(adjustl(line))
       ll = len(line)
       if (ll > 4) then
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
             seed%atname(nat) = word

             is = usen%get(lword,1)
             seed%is(nat) = is
             if (.not.usespc(is)) then
                seed%spc(is)%name = trim(word)
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
    seed%isformat = isformat_r_castepgeom

  end subroutine read_castep_geom

  !> Read the structure from a CASTEP phonon file
  module subroutine read_castep_phonon(seed,file,mol,errmsg,ti)
    use tools_io, only: nameguess, zatguess, fopen_read, fclose, getline_raw, equalsub
    use param, only: isformat_r_castepphonon
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line
    integer :: lu, i, idum, iz, nat
    character*20 :: strname
    real*8 :: x(3)
    type(rawseed) :: rec

    call seed%end()
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! read the header (angstrom, fractional coordinates)
    nat = 0
    rec%spcmode = spc_z
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_ang
    do while(getline_raw(lu,line))
       if (equalsub(line,1,15," Number of ions")) then
          read(line(16:),*,end=999,err=999) nat
       elseif (equalsub(line,1,22," Unit cell vectors (A)")) then
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) rec%rv(:,1)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) rec%rv(:,2)
          if (.not.getline_raw(lu,line)) goto 999
          read(line,*,err=999,end=999) rec%rv(:,3)
       elseif (equalsub(line,1,24," Fractional Co-ordinates")) then
          if (nat == 0) goto 999
          do i = 1, nat
             if (.not.getline_raw(lu,line)) goto 999
             read(line,*,end=999,err=999) idum, x, strname
             iz = zatguess(strname)
             if (iz <= 0) goto 999
             call rawseed_add_atom(rec,iz,x,trim(adjustl(strname)),spcname=trim(nameguess(iz,.true.)))
          end do
       end if
    end do
    if (rec%nat == 0) goto 999

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_castepphonon,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_castep_phonon

  !> Read the structure from a DMACRYS input file (dmain)
  module subroutine read_dmain(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, isreal,&
       getword, zatguess, isinteger, lower
    use tools_math, only: matinv
    use param, only: isformat_r_dmain
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, lword
    character*2 :: atsym
    logical :: havelatt, haverlscal, ok
    integer :: lu, lp, iz, lvl, nrec, ier
    real*8 :: r(3,3), rlscal, x(3)
    integer :: i
    type(rawseed) :: rec

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
    rec%spcmode = spc_z
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

             ! write down the atom (species named by the two-letter symbol)
             ok = isreal(x(1),line,lp)
             ok = ok .and. isreal(x(2),line,lp)
             ok = ok .and. isreal(x(3),line,lp)
             if (.not.ok) goto 999
             call rawseed_add_atom(rec,iz,x,word,spcname=atsym)

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
          exit
       elseif (word == "stop") then
          exit
       end if
    end do

    ! lattice vectors (angstrom, scaled) and fractional coordinates from the
    ! unscaled lattice, as DMACRYS does
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_ang
    rec%rv = transpose(r) * rlscal
    r = transpose(r)
    call matinv(r,3,ier)
    if (ier /= 0) goto 999
    do i = 1, rec%nat
       rec%x(:,i) = matmul(r,rec%x(:,i))
    end do

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_dmain,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_dmain

  !> Read the structure from a file in DFTB+ gen format.
  module subroutine read_dftbp(seed,file,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline, lower, equal, &
       getword, zatguess, nameguess, fclose
    use param, only: isformat_r_gen
    use types, only: realloc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, iz, idum, lp
    logical :: ok, molout
    character*1 :: isfrac
    character(len=:), allocatable :: line, word
    type(rawseed) :: rec

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
    read (line,*,err=999,end=999) rec%nat, isfrac
    isfrac = lower(isfrac)
    if (.not.(equal(isfrac,"f").or.equal(isfrac,"c").or.equal(isfrac,"s"))) then
       errmsg = 'Wrong coordinate selector.'
       goto 999
    end if
    allocate(rec%x(3,rec%nat),rec%is(rec%nat))

    ! atom types
    rec%nspc = 0
    allocate(rec%spc(2))
    ok = getline(lu,line)
    if (.not.ok) goto 999
    lp = 1
    word = getword(line,lp)
    iz = zatguess(word)
    do while (iz >= 0)
       rec%nspc = rec%nspc + 1
       if (rec%nspc > size(rec%spc,1)) &
          call realloc(rec%spc,2*rec%nspc)
       rec%spc(rec%nspc)%z = iz
       rec%spc(rec%nspc)%name = nameguess(iz,.true.)
       word = getword(line,lp)
       iz = zatguess(word)
    end do
    if (rec%nspc == 0) then
       errmsg = 'No atomic types found.'
       goto 999
    end if
    call realloc(rec%spc,rec%nspc)

    ! read atomic positions (Cartesian coordinates converted to bohr)
    do i = 1, rec%nat
       ok = getline(lu,line)
       if (.not.ok) goto 999
       read (line,*,err=999,end=999) idum, rec%is(i), rec%x(:,i)
    end do
    rec%spcmode = spc_given
    rec%lunit = lunit_ang
    rec%iscart = (isfrac /= "f")
    rec%border = rborder
    rec%cubic = docube

    ! read lattice vectors, if they exist
    ok = getline(lu,line)
    if (ok) then
       do i = 1, 3
          ok = getline(lu,line,.true.)
          read (line,*) r(i,:)
       end do
       if (isfrac == "c") then
          errmsg = 'Lattice plus C not supported.'
          goto 999
       end if
       rec%ndim = 3
       rec%cellmode = 2
       rec%rv = transpose(r)
       molout = .false.
    else
       ! molecule and no lattice -> set up the origin and the molecular cell
       if (isfrac == "f" .or. isfrac == "s") then
          errmsg = 'S or F coordinates but no lattice vectors.'
          goto 999
       end if
       rec%ndim = 0
       rec%ismol = .true.
       molout = .true.
    end if

    call rawseed_to_seed(rec,seed,molout,file,isformat_r_gen,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_dftbp

  !> Read the structure from an xsf file.
  module subroutine read_xsf(seed,file,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, nameguess, equal,&
       zatguess, isinteger, getword, isreal
    use param, only: isformat_r_xsf
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, name
    character*10 :: atn
    integer :: lu, lp, i, iz, nat
    real*8 :: r(3,3), x(3)
    logical :: ok, ismol, xread
    type(rawseed) :: rec

    ! open
    call seed%end()
    ismol = .false.
    errmsg = ""
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! angstrom, Cartesian coordinates
    xread = .false.
    errmsg = "Error reading file: " // trim(file)
    rec%lunit = lunit_ang
    rec%iscart = .true.
    rec%border = rborder
    rec%cubic = docube
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"primvec")) then
          do i = 1, 3
             read (lu,*,err=999,end=999) r(i,:)
          end do
          ismol = .false.
       elseif (equal(word,"primcoord").and..not.xread) then
          ! crystal: species by atomic number, named by the first label
          read (lu,*,err=999,end=999) nat
          rec%spcmode = spc_z
          do i = 1, nat
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
             ok = isreal(x(1),line,lp)
             ok = ok.and.isreal(x(2),line,lp)
             ok = ok.and.isreal(x(3),line,lp)
             if (.not.ok) then
                errmsg = 'Wrong atomic position.'
                goto 999
             end if
             call rawseed_add_atom(rec,iz,x,name)
          end do
          ismol = .false.
          xread = .true.
       elseif (equal(word,"atoms").and..not.xread) then
          ! molecule: species by label (case-insensitive)
          ismol = .true.
          rec%spcmode = spc_label_nocase
          do while (getline_raw(lu,line))
             if (len_trim(line) == 0) exit
             read (line,*,err=999,end=999) atn, x
             ok = isinteger(iz,atn)
             if (.not.ok) then
                iz = zatguess(atn)
             end if
             if (iz < 0) then
                errmsg = "Unknown atomic symbol: "//trim(atn)//"."
                goto 999
             end if
             call rawseed_add_atom(rec,iz,x,trim(atn))
          end do
          xread = .true.
       end if
    end do
    if (rec%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if

    if (.not.ismol) then
       rec%ndim = 3
       rec%cellmode = 2
       rec%rv = transpose(r)
    else
       rec%ndim = 0
       rec%ismol = .true.
    end if

    call rawseed_to_seed(rec,seed,ismol,file,isformat_r_xsf,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_xsf

  !> Read the structure from a pwc file.
  module subroutine read_pwc(seed,file,mol,errmsg,ti)
    use tools_io, only: fopen_read, fclose, zatguess
    use param, only: isformat_r_pwc
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    integer :: version, i
    character*3, allocatable :: atm(:)
    real*8 :: alat
    type(rawseed) :: rec

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

    read (lu,err=999,end=999) rec%nspc, rec%nat, alat

    ! species
    allocate(atm(rec%nspc),rec%spc(rec%nspc))
    read (lu,err=999,end=999) atm
    do i = 1, rec%nspc
       rec%spc(i)%name = trim(atm(i))
       rec%spc(i)%z = zatguess(rec%spc(i)%name)
    end do
    deallocate(atm)

    ! read the rest (atomic units, Cartesian coordinates in alat units)
    allocate(rec%x(3,rec%nat),rec%is(rec%nat))
    read (lu,err=999,end=999) rec%is
    read (lu,err=999,end=999) rec%x
    read (lu,err=999,end=999) rec%rv
    rec%rv = rec%rv * alat
    rec%x = rec%x * alat
    rec%spcmode = spc_given
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_bohr
    rec%iscart = .true.

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_pwc,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_pwc

  !> Read the structure from an axsf file (xcrysden). Read the
  !> coordinates from PRIMCOORD block and nudge them using the
  !> eigenvector on the same block scaled by the value of xnudge
  !> (bohr).
  module subroutine read_axsf(seed,file,nread0,xnudge,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isinteger, &
       string, getword, isreal, nameguess, zatguess
    use param, only: bohrtoa, isformat_r_axsf
    class(crystalseed), intent(inout) :: seed !< Crystal seed output
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: nread0
    real*8, intent(in) :: xnudge
    real*8, intent(in) :: rborder !< user-defined border in bohr
    logical, intent(in) :: docube !< if true, make the cell cubic
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, name
    integer :: lu, lp, iprim, i, iz, nat
    real*8 :: r(3,3), x(3), xd(3)
    logical :: ok, didreadr, didreadx

    type(rawseed) :: rec

    ! open
    call seed%end()
    errmsg = ""
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! atomic units, Cartesian coordinates, species by atomic number
    errmsg = "Error reading file: " // trim(file)
    rec%spcmode = spc_z
    rec%ndim = 3
    rec%cellmode = 2
    rec%lunit = lunit_bohr
    rec%iscart = .true.
    rec%border = rborder
    rec%cubic = docube
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
       elseif (equal(word,"primcoord")) then
          ok = isinteger(iprim,line,lp)
          if (iprim == nread0) then
             didreadx = .true.
             read (lu,*,err=999,end=999) nat
             do i = 1, nat
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
                ok = isreal(x(1),line,lp)
                ok = ok.and.isreal(x(2),line,lp)
                ok = ok.and.isreal(x(3),line,lp)
                if (.not.ok) then
                   errmsg = 'Wrong atomic position.'
                   goto 999
                end if
                x = x / bohrtoa

                ! read the displacement vector and apply the nudge
                ok = isreal(xd(1),line,lp)
                ok = ok.and.isreal(xd(2),line,lp)
                ok = ok.and.isreal(xd(3),line,lp)
                if (.not.ok) then
                   errmsg = 'Wrong displacement vector.'
                   goto 999
                end if
                x = x + xnudge * xd
                call rawseed_add_atom(rec,iz,x,name)
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
    if (rec%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    rec%rv = transpose(r)

    call rawseed_to_seed(rec,seed,.false.,file,isformat_r_axsf,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_axsf

  !> Read an FHIaims input file.
  module subroutine read_aimsin(seed,file,mol,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isreal, &
       getword, zatguess
    use param, only: isformat_r_aimsin
    class(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nlat, iz
    logical :: is_file_mol, ok
    character(len=:), allocatable :: line, word
    real*8 :: rlat(3,3), x(3)
    type(rawseed) :: rec

    ! open
    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! collect the information (angstrom; per-atom fractional flag)
    rec%spcmode = spc_label
    rec%lunit = lunit_ang
    rec%iscart = .true.
    rec%border = rborder
    rec%cubic = docube
    is_file_mol = .true.
    nlat = 0
    do while (getline_raw(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (len_trim(word) == 0) cycle
       if (word(1:1) == "#") cycle

       if (equal(word,'atom').or.equal(word,'atom_frac')) then
          ok = isreal(x(1),line,lp)
          ok = ok .and. isreal(x(2),line,lp)
          ok = ok .and. isreal(x(3),line,lp)
          if (.not.ok) goto 999
          ok = equal(word,'atom_frac')
          word = getword(line,lp)
          iz = zatguess(word)
          if (iz <= 0) then
             errmsg = "unknown atom type: " // word
             goto 999
          end if
          call rawseed_add_atom(rec,iz,x,word,isfrac=ok)
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

    ! handle the molecule/crystal expectation/contents of the file
    if (is_file_mol) then
       rec%ismol = .true.
       rec%ndim = 0
    else
       if (mol) then
          errmsg = "tried to load a crystal as a molecule; use the CRYSTAL keyword"
          goto 999
       end if
       rec%ndim = 3
       rec%cellmode = 2
       rec%rv = rlat
    end if

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_aimsin,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_aimsin

  !> Read an FHI aims output file.
  module subroutine read_aimsout(seed,file,mol,rborder,docube,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, equal, isreal, &
       getword, zatguess, equalsub
    use tools_math, only: matinv
    use types, only: realloc
    use hashmod, only: hash
    use param, only: bohrtoa, hartoev, eva3togpa, isformat_r_aimsout
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
    allocate(isfrac(10),seed%x(3,10),seed%is(10),seed%atname(10))

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
          call realloc(seed%atname,2*seed%nat)
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
       seed%atname(seed%nat) = word
    end do
    call realloc(isfrac,seed%nat)
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%atname,seed%nat)

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
          if (line == "  Fractional coordinates:") exit
          if (len_trim(line) == 0) cycle
          if (equalsub(line,1,4,"----")) exit
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
             seed%atname(iat) = getword(line,lp)
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
    seed%isformat = isformat_r_aimsout

  end subroutine read_aimsout

  !> Read the structure from a TINKER frac file (must have cell
  !> parameters in the second line).
  module subroutine read_tinkerfrac(seed,file,mol,errmsg,ti)
    use tools_io, only: getline_raw, fopen_read, fclose, isinteger, zatguess, nameguess
    use param, only: maxzat, isformat_r_tinkerfrac
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, i, iz, idum, nat
    logical :: ok
    character(len=:), allocatable :: line, aux
    character*30 :: atsym
    real*8 :: x(3)
    type(rawseed) :: rec

    ! open the file
    call seed%end()
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! line 1: number of atoms and name of the seed
    lp = 1
    ok = getline_raw(lu,line,.false.)
    ok = ok .and. isinteger(nat,line,lp)
    if (.not.ok) goto 999

    ! line 2: cell parameters (angstrom, fractional coordinates)
    read (lu,*,err=999,end=999) rec%aa, rec%bb
    rec%ndim = 3
    rec%cellmode = 1
    rec%lunit = lunit_ang
    rec%spcmode = spc_z
    rec%havex0 = .true.

    ! atoms
    do i = 1, nat
       read(lu,*,err=999,end=999) idum, atsym, x
       iz = zatguess(atsym)
       if (iz <= 0 .or. iz > maxzat) then
          errmsg = "Unknown atomic symbol."
          goto 999
       endif
       aux = adjustl(atsym)
       call rawseed_add_atom(rec,iz,x,aux(1:10),spcname=nameguess(iz,.true.))
    end do

    call rawseed_to_seed(rec,seed,mol,file,isformat_r_tinkerfrac,errmsg)
999 continue
    call fclose(lu)

  end subroutine read_tinkerfrac

  !> Read the structure from an xyz file. If mol0 == 1, force a
  !> molecule. If mol0 == 0, force a crystal. If mol0 == -1, let the
  !> format detection decide. If the read was successful, return an
  !> empty error message
  module subroutine read_xyz(seed,file,mol,rborder,docube,errmsg,ti)
    class(crystalseed), intent(inout) :: seed !< Output crystal seed
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: mol !< is this a molecule?
    real*8, intent(in) :: rborder
    logical, intent(in) :: docube
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed), allocatable :: seedaux(:)
    integer :: nseed, imol

    if (mol) then
       imol = 1
    else
       imol = 0
    end if
    call seed%end()
    call read_all_xyz(nseed,seedaux,file,imol,errmsg,seed0=seed,ti=ti)
    seed%border = rborder
    seed%cubic = docube

  end subroutine read_xyz

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

  !> Store the site occupancies occ(1:nat) into the seed. The values are clamped
  !> to (0,1]. seed%occ is allocated only when some site is partially occupied
  !> (occ < 1); otherwise it is left deallocated, which means all occupancies
  !> are 1 (see the seed%occ contract in crystalseedmod).
  module subroutine seed_set_occ(seed,occ,nat)
    use global, only: occ_eps
    type(crystalseed), intent(inout) :: seed
    real*8, intent(in) :: occ(:)
    integer, intent(in) :: nat

    integer :: i

    if (allocated(seed%occ)) deallocate(seed%occ)
    if (nat <= 0) return
    if (.not.any(occ(1:nat) < 1d0-occ_eps)) return
    allocate(seed%occ(nat))
    do i = 1, nat
       seed%occ(i) = max(min(occ(i),1d0),1d-10)
    end do

  end subroutine seed_set_occ

  !> Detect the format for the structure-containing file. Normally,
  !> this works by detecting the extension, but the file may be
  !> opened and searched if ambiguity is present. The format and
  !> whether the file contains a molecule or crysatl is returned.
  !> If alsofield is present, then return .true. if the file also
  !> contains a scalar field.
  module subroutine struct_detect_read_format(file,isformat,alsofield,ti)
    use param, only: isformat_r_unknown, isformat_r_cif, isformat_r_shelx,&
       isformat_r_f21, isformat_r_xyz, isformat_r_gjf,&
       isformat_r_cube, isformat_r_bincube, isformat_r_struct, isformat_r_abinit, isformat_r_elk,&
       isformat_r_wfn, isformat_r_wfx, isformat_r_fchk, isformat_r_molden,&
       isformat_r_gaussian, isformat_r_siesta, isformat_r_xsf, isformat_r_gen,&
       isformat_r_vasp, isformat_r_pwc, isformat_r_axsf, isformat_r_dat, isformat_r_pgout,&
       isformat_r_dmain, isformat_r_aimsin, isformat_r_aimsout, isformat_r_tinkerfrac,&
       isformat_r_castepcell, isformat_r_castepgeom, isformat_r_castepphonon,&
       isformat_r_qein, isformat_r_qeout, isformat_r_xband, isformat_r_gulpin, isformat_r_gulpout,&
       isformat_r_mol2, isformat_r_pdb, isformat_r_zmat, isformat_r_sdf, isformat_r_magres
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
       (index(basename,'POSCAR') > 0).or.(index(basename,'PARCHG') > 0)

    if (equal(lower(wextdot),'cif')) then
       isformat = isformat_r_cif
    elseif (equal(wextdot,'pwc')) then
       isformat = isformat_r_pwc
       alsofield_ = .true.
    elseif (equal(wextdot,'res').or.equal(wextdot,'ins').or.equal(wextdot,'16')) then
       isformat = isformat_r_shelx
    elseif (equal(wextdot,'21')) then
       isformat = isformat_r_f21
    elseif (equal(wextdot,'cube').or.equal(wextdot,'cub')) then
       isformat = isformat_r_cube
       alsofield_ = .true.
    elseif (equal(wextdot,'bincube')) then
       isformat = isformat_r_bincube
       alsofield_ = .true.
    elseif (equal(wextdot,'struct')) then
       isformat = isformat_r_struct
    elseif (equal(wextdot,'DEN').or.equal(wext_,'DEN').or.equal(wextdot,'ELF').or.equal(wext_,'ELF').or.&
       equal(wextdot,'POT').or.equal(wext_,'POT').or.equal(wextdot,'VHA').or.equal(wext_,'VHA').or.&
       equal(wextdot,'VHXC').or.equal(wext_,'VHXC').or.equal(wextdot,'VXC').or.equal(wext_,'VXC').or.&
       equal(wextdot,'GDEN1').or.equal(wext_,'GDEN1').or.equal(wextdot,'GDEN2').or.equal(wext_,'GDEN2').or.&
       equal(wextdot,'GDEN3').or.equal(wext_,'GDEN3').or.equal(wextdot,'LDEN').or.equal(wext_,'LDEN').or.&
       equal(wextdot,'KDEN').or.equal(wext_,'KDEN').or.equal(wextdot,'PAWDEN').or.equal(wext_,'PAWDEN').or.&
       equal(wextdot,'VCLMB').or.equal(wext_,'VCLMB').or.equal(wextdot,'VPSP').or.equal(wext_,'VPSP')) then
       isformat = isformat_r_abinit
       alsofield_ = .true.
    elseif (equal(wextdot,'OUT')) then
       isformat = isformat_r_elk
    elseif (equal(wextdot,'out')) then
       call which_out_format(file,isformat,ti=ti)
    elseif (equal(wextdot,'own')) then
       call which_out_format(file,isformat,ti=ti)
       if (isformat /= isformat_r_aimsout) goto 999
    elseif (equal(wextdot,'in')) then
       call which_in_format(file,isformat,ti=ti)
    elseif (equal(wextdot2,'in.next_step')) then
       call which_in_format(file,isformat,ti=ti)
       if (isformat /= isformat_r_aimsin) goto 999
    elseif (equal(wextdot,'sys')) then
       isformat = isformat_r_xband
    elseif (equal(wextdot,'pwi')) then
       isformat = isformat_r_qein
    elseif (equal(wextdot,'pwo')) then
       isformat = isformat_r_qeout
    elseif (equal(wextdot,'xyz')) then
       isformat = isformat_r_xyz
    elseif (equal(wextdot,'gjf').or.equal(wextdot,'com')) then
       isformat = isformat_r_gjf
    elseif (equal(wextdot,'zmat')) then
       isformat = isformat_r_zmat
    elseif (equal(wextdot,'sdf').or.equal(wextdot,'mol')) then
       isformat = isformat_r_sdf
    elseif (equal(wextdot,'magres')) then
       isformat = isformat_r_magres
    elseif (equal(wextdot,'pgout')) then
       isformat = isformat_r_pgout
    elseif (equal(wextdot,'wfn')) then
       isformat = isformat_r_wfn
       alsofield_ = .true.
    elseif (equal(wextdot,'wfx')) then
       isformat = isformat_r_wfx
       alsofield_ = .true.
    elseif (equal(wextdot,'log')) then
       isformat = isformat_r_gaussian
       alsofield_ = .false.
    elseif (equal(wextdot,'fchk')) then
       isformat = isformat_r_fchk
       alsofield_ = .true.
    elseif (equal(wextdot,'molden') .or. equal(wextdot2,'molden.input')) then
       isformat = isformat_r_molden
       alsofield_ = .true.
    elseif (equal(wextdot,'STRUCT_OUT').or.equal(wextdot,'STRUCT_IN')) then
       isformat = isformat_r_siesta
    elseif (equal(wextdot,'cell')) then
       isformat = isformat_r_castepcell
    elseif (equal(wextdot,'phonon')) then
       isformat = isformat_r_castepphonon
    elseif (equal(wextdot,'geom')) then
       isformat = isformat_r_castepgeom
    elseif (equal(wextdot,'dmain')) then
       isformat = isformat_r_dmain
    elseif (equal(wextdot,'xsf')) then
       isformat = isformat_r_xsf
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
       isformat = isformat_r_gen

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
       isformat = isformat_r_axsf
    elseif (equal(wextdot,'dat')) then
       isformat = isformat_r_dat
    elseif (equal(wextdot,'frac')) then
       isformat = isformat_r_tinkerfrac
    elseif (equal(lower(wextdot),'vasp')) then
       isformat = isformat_r_vasp
       alsofield_ = .false.
    elseif (isvasp) then
       isformat = isformat_r_vasp
       alsofield_ = (index(basename,'CHGCAR') > 0) .or. (index(basename,'CHG') > 0) .or. &
          (index(basename,'ELFCAR') > 0) .or. (index(basename,'AECCAR0') > 0) .or. &
          (index(basename,'AECCAR1') > 0) .or. (index(basename,'AECCAR2') > 0).or.&
          (index(basename,'PARCHG') > 0)
    elseif (equal(wextdot,'mol2')) then
       isformat = isformat_r_mol2
    elseif (equal(wextdot,'pdb')) then
       isformat = isformat_r_pdb
    elseif (equal(lower(wextdot),'gin') .or. equal(lower(wextdot),'grs')) then
       isformat = isformat_r_gulpin
    elseif (equal(lower(wextdot),'gout') .or. equal(lower(wextdot),'got')) then
       isformat = isformat_r_gulpout
    else
       goto 999
    endif
    if (present(alsofield)) alsofield = alsofield_

    if (isformat == 0) goto 999

    return
999 continue
    isformat = isformat_r_unknown

  end subroutine struct_detect_read_format

  !> Detect whether a file with format isformat contains a molecule
  !> (ismol=.true.)  or a crystal (.false.)
  module subroutine struct_detect_ismol(file,isformat,ismol,ti)
    use tools_io, only: fopen_read, fclose, getline, equali,&
       getline_raw, lgetword, equal, lower, isreal
    use param, only: isformat_r_cif, isformat_r_shelx, isformat_r_f21,&
       isformat_r_cube, isformat_r_bincube, isformat_r_struct, isformat_r_abinit, isformat_r_elk,&
       isformat_r_qein, isformat_r_qeout, isformat_r_crystal, isformat_r_fploout,&
       isformat_r_xyz, isformat_r_gjf, isformat_r_wfn,&
       isformat_r_wfx, isformat_r_fchk, isformat_r_molden, isformat_r_gaussian, isformat_r_siesta,&
       isformat_r_xsf, isformat_r_gen, isformat_r_vasp, isformat_r_pwc, isformat_r_axsf,&
       isformat_r_dat, isformat_r_pgout, isformat_r_orca, isformat_r_dmain, isformat_r_aimsin,&
       isformat_r_aimsout, isformat_r_tinkerfrac, isformat_r_castepcell, isformat_r_castepphonon,&
       isformat_r_castepgeom,&
       isformat_r_mol2, isformat_r_pdb, isformat_r_zmat, isformat_r_sdf, isformat_r_magres,&
       isformat_r_alamode, isformat_r_akaikkr, isformat_r_xband, isformat_r_gulpin,&
       isformat_r_gulpout
    character*(*), intent(in) :: file
    integer, intent(in) :: isformat
    logical, intent(out) :: ismol
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word, errmsg
    integer :: lu, ios, lp, nat, idx, i
    character*1 :: isfrac
    logical :: ok
    real*8 :: scalexx(3)

    real*8, parameter :: eps = 1d-10

    ismol = .false.
    select case (isformat)
    case (isformat_r_cif,isformat_r_pwc,isformat_r_shelx,isformat_r_f21,&
       isformat_r_struct,isformat_r_abinit,&
       isformat_r_elk,isformat_r_siesta,isformat_r_dmain,isformat_r_vasp,&
       isformat_r_axsf,isformat_r_tinkerfrac,isformat_r_qein,isformat_r_qeout,&
       isformat_r_crystal,isformat_r_fploout,isformat_r_castepcell,isformat_r_castepphonon,&
       isformat_r_castepgeom,isformat_r_magres,isformat_r_alamode,isformat_r_akaikkr,&
       isformat_r_xband)
       ! these are always crystals
       ismol = .false.

    case (isformat_r_gjf,isformat_r_pgout,isformat_r_wfn,isformat_r_wfx,&
       isformat_r_gaussian,isformat_r_fchk,isformat_r_molden,isformat_r_dat,&
       isformat_r_orca,isformat_r_mol2,isformat_r_zmat,isformat_r_sdf)
       ! these are always molecules
       ismol = .true.

    case(isformat_r_xyz)
       ! detect if this is an ASE file: check for "pbc" and "lattice *="
       ismol = .true.
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) return
       ok = getline_raw(lu,line)
       ok = ok.and.getline_raw(lu,line)
       if (.not.ok) return
       line = lower(line)

       idx = index(line,"pbc")
       if (idx > 0) then
          idx = index(line,"lattice")
          if (idx > 0) then
             ismol = .false.
             do i = idx+7,len(line)
                if (line(i:i) == " ") cycle ! skip blanks
                if (line(i:i) /= "=") then
                   ismol = .true.
                   exit
                end if
                exit
             end do
          end if
       end if
       call fclose(lu)

    case(isformat_r_pdb)
       ! In a pdb file: we have a crystal if SCALE1 and the other
       ! parameters are present. Unfortunately, these may appear in
       ! the file but be bogus, which wrecks the geometry. So we do
       ! not use them if they have unreasonable values.
       ismol = .true.
       lu = fopen_read(file,errstop=.false.,ti=ti)
       if (lu < 0) return
       do while(getline_raw(lu,line))
          lp = 1
          word = lgetword(line,lp)
          if (equal(word,"scale1")) then
             ok = isreal(scalexx(1),line,lp)
             ok = ok .and. isreal(scalexx(2),line,lp)
             ok = ok .and. isreal(scalexx(3),line,lp)
             if (ok) then
                if (abs(scalexx(1)-1d0)>eps .or. abs(scalexx(2))>eps .or. abs(scalexx(1)-1d0)>eps) then
                   ismol = .false.
                end if
             end if
             exit
          end if
       end do
       call fclose(lu)

    case(isformat_r_aimsout)
       ! Detect from the FHI-aims output file whether this is a
       ! molecule or a crystal.
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

    case (isformat_r_aimsin)
       ! FHIaims input: a crystal if we have lattice_vectors
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

    case (isformat_r_gen)
       ! A gen file: this is a crystal if we have "F" and a molecule
       ! if "C" in the first line.
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

    case (isformat_r_xsf)
       ! xsf format: a molecule if the first entry is ATOMS
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

    case (isformat_r_cube)
       ! cube format: a molecule if atoms are outside the box spanned by the grid
       ismol = ismol_cube(file,.false.,errmsg,ti=ti)

    case (isformat_r_bincube)
       ! bincube format: a molecule if atoms are outside the box spanned by the grid
       ismol = ismol_cube(file,.true.,errmsg,ti=ti)

    case (isformat_r_gulpin,isformat_r_gulpout)
       ! GULP files: a molecule if the first structure is a 0D cluster
       call gulp_detect_ismol(file,isformat,ismol,ti=ti)

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
  !> if none). If iavib is present, return the seed number for which
  !> the file can be read to obtain vibrational data. If mol0 == 1,
  !> force reading molecules, if mol0 == 0, force reading crystals, if
  !> mol0 == -1, let the routine figure it out. Returns the number of
  !> seeds read (nseed), the seeds themselves (seed) and
  !> collapse=.true. if the seeds are the steps of a geometry
  !> optimization.
  module subroutine read_seeds_from_file(file,mol0,isformat0,readlastonly,&
     nseed,seed,collapse,errmsg,iafield,iavib,ti)
    use global, only: rborder_def, doguess
    use tools_io, only: getword, equali, fopen_read, fclose
    use param, only: isformat_r_cube, isformat_r_bincube, isformat_r_xyz, isformat_r_wfn,&
       isformat_r_wfx, isformat_r_fchk, isformat_r_molden, isformat_r_gaussian, isformat_r_gjf,&
       isformat_r_zmat, isformat_r_abinit, isformat_r_cif, isformat_r_pwc, isformat_r_fploout,&
       isformat_r_crystal, isformat_r_elk, isformat_r_gen, isformat_r_qein, isformat_r_qeout,&
       isformat_r_shelx, isformat_r_siesta, isformat_r_struct, isformat_r_vasp, isformat_r_axsf,&
       isformat_r_xsf, isformat_r_castepcell, isformat_r_castepphonon, isformat_r_castepgeom,&
       isformat_r_dat, isformat_r_f21, isformat_r_unknown, isformat_r_pgout, isformat_r_orca,&
       isformat_r_dmain, isformat_r_aimsin, isformat_r_aimsout, isformat_r_tinkerfrac,&
       isformat_r_mol2, isformat_r_sdf, isformat_r_pdb, isformat_r_magres,&
       isformat_r_alamode, isformat_r_akaikkr, isformat_r_xband, isformat_r_gulpin,&
       isformat_r_gulpout
    character*(*), intent(in) :: file
    integer, intent(in) :: mol0
    integer, intent(in) :: isformat0
    logical, intent(in) :: readlastonly
    integer, intent(out) :: nseed
    type(crystalseed), allocatable, intent(inout) :: seed(:)
    logical, intent(out) :: collapse
    character(len=:), allocatable, intent(out) :: errmsg
    integer, intent(out), optional :: iafield
    integer, intent(out), optional :: iavib
    type(thread_info), intent(in), optional :: ti

    integer :: isformat, is, mol0_, i, lu
    logical :: mol, alsofield, ok, alsovib

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

    call struct_detect_read_format(file,is,alsofield,ti=ti)
    if (isformat == isformat_r_unknown) isformat = is
    if (isformat == isformat_r_unknown) then
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
    alsovib = .false.
    if (isformat == isformat_r_cif) then
       call read_all_cif(file,mol,errmsg,nseed=nseed,mseed=seed,ti=ti)
    elseif (isformat == isformat_r_mol2) then
       call read_all_mol2(file,errmsg,nseed=nseed,mseed=seed,ti=ti)
    elseif (isformat == isformat_r_sdf) then
       call read_all_sdf(file,errmsg,nseed=nseed,mseed=seed,ti=ti)
    elseif (isformat == isformat_r_magres) then
       call seed(1)%read_magres(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_alamode) then
       call seed(1)%read_alamode(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_akaikkr) then
       call seed(1)%read_akaikkr(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_xband) then
       call seed(1)%read_xband(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_pwc) then
       call seed(1)%read_pwc(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_shelx) then
       call seed(1)%read_shelx(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_f21) then
       call seed(1)%read_f21(file,mol,errmsg,ti=ti)
    else if (isformat == isformat_r_cube) then
       call seed(1)%read_cube(file,mol,errmsg,ti=ti)
    else if (isformat == isformat_r_bincube) then
       call seed(1)%read_bincube(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_struct) then
       call seed(1)%read_wien(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_vasp) then
       call read_all_vasp(nseed,seed,file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_abinit) then
       call seed(1)%read_abinit(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_elk) then
       call seed(1)%read_elk(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_qeout) then
       call read_all_qeout(nseed,seed,file,mol,-1,errmsg,ti=ti)
    elseif (isformat == isformat_r_gulpin) then
       call read_all_gulpin(nseed,seed,file,mol,-1,errmsg,ti=ti)
    elseif (isformat == isformat_r_gulpout) then
       call read_all_gulpout(nseed,seed,file,mol,-1,errmsg,ti=ti)
    elseif (isformat == isformat_r_crystal) then
       call read_all_crystalout(file,mol,errmsg,nseed=nseed,mseed=seed,alsovib=alsovib,ti=ti)
    elseif (isformat == isformat_r_fploout) then
       call seed(1)%read_fploout(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_qein) then
       call seed(1)%read_qein(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_xyz) then
       call read_all_xyz(nseed,seed,file,-1,errmsg,ti=ti)
    elseif (isformat == isformat_r_gaussian) then
       call read_all_log(nseed,seed,file,errmsg,alsovib=alsovib,ti=ti)
    elseif (isformat == isformat_r_wfn .or. isformat == isformat_r_wfx.or.&
       isformat == isformat_r_fchk.or.isformat == isformat_r_molden.or.&
       isformat == isformat_r_dat.or.isformat == isformat_r_pgout.or.&
       isformat == isformat_r_orca.or.isformat == isformat_r_gjf.or.&
       isformat == isformat_r_zmat) then
       call seed(1)%read_mol(file,isformat,rborder_def,.false.,errmsg,alsovib=alsovib,ti=ti)
    elseif (isformat == isformat_r_pdb) then
       call seed(1)%read_pdb(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_siesta) then
       call seed(1)%read_siesta(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_castepcell) then
       call seed(1)%read_castep_cell(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_castepphonon) then
       call seed(1)%read_castep_phonon(file,mol,errmsg,ti=ti)
       alsovib = .true.
    elseif (isformat == isformat_r_castepgeom) then
       call read_all_castep_geom(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_r_dmain) then
       call seed(1)%read_dmain(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_aimsin) then
       call seed(1)%read_aimsin(file,mol,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_r_aimsout) then
       call read_all_aimsout(nseed,seed,file,errmsg,ti=ti)
    elseif (isformat == isformat_r_tinkerfrac) then
       call seed(1)%read_tinkerfrac(file,mol,errmsg,ti=ti)
    elseif (isformat == isformat_r_axsf) then
       call seed(1)%read_axsf(file,1,0d0,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_r_xsf) then
       call seed(1)%read_xsf(file,rborder_def,.false.,errmsg,ti=ti)
    elseif (isformat == isformat_r_gen) then
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

    ! output iavib
    if (present(iavib)) then
       if (alsovib) then
          iavib = nseed
       else
          iavib = 0
       end if
    end if

    ! output collapse
    collapse = ((isformat == isformat_r_qeout .or. isformat == isformat_r_gaussian .or.&
       isformat == isformat_r_aimsout .or. isformat == isformat_r_castepgeom .or.&
       isformat == isformat_r_crystal .or. isformat == isformat_r_gulpout).and.nseed > 1)

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
       write (uout,'(99(A," "))') string(seed%spc(i)%z), string(seed%name)
    end do
    write (uout,*)

    write (uout,'("## atomic positions (",A,")")') string(seed%nat)
    do i = 1, seed%nat
       write (uout,'(99(A," "))') string(seed%is(i)), &
          (string(seed%x(j,i),'f',decimal=10),j=1,3)
    end do
    write (uout,*)

    write (uout,'("## cell")')
    write (uout,'("useabr = ",A)') string(seed%useabr)
    if (seed%useabr == 1) then
       write (uout,'("aa = ",3(A," "))') (string(seed%aa(i),'f',decimal=8),i=1,3)
       write (uout,'("bb = ",3(A," "))') (string(seed%bb(i),'f',decimal=5),i=1,3)
    else
       write (uout,'("x2c = ",3(A," "))') (string(seed%m_x2c(1,i),'f',decimal=8),i=1,3)
       write (uout,'("x2c = ",3(A," "))') (string(seed%m_x2c(2,i),'f',decimal=8),i=1,3)
       write (uout,'("x2c = ",3(A," "))') (string(seed%m_x2c(3,i),'f',decimal=8),i=1,3)
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
    write (uout,'("molx0 = ",3(A," "))') (string(seed%molx0(i),'f',decimal=4),i=1,3)
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
    to%libname = from%libname
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
    if (allocated(from%atname)) then
       to%atname = from%atname
    else
       if (allocated(to%atname)) deallocate(to%atname)
    end if
    if (allocated(from%occ)) then
       to%occ = from%occ
    else
       if (allocated(to%occ)) deallocate(to%occ)
    end if
    if (allocated(from%mix)) then
       to%mix = from%mix
    else
       if (allocated(to%mix)) deallocate(to%mix)
    end if
    to%neqlist = from%neqlist
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
    to%havebonds = from%havebonds
    if (allocated(from%nstar)) then
       to%nstar = from%nstar
    else
       if (allocated(to%nstar)) deallocate(to%nstar)
    end if
    to%ismolecule = from%ismolecule
    to%cubic = from%cubic
    to%border = from%border
    to%havex0 = from%havex0
    to%molx0 = from%molx0
    to%energy = from%energy
    to%pressure = from%pressure

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
    use param, only: bohrtoa, dirsep, isformat_r_vasp
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
          ofile = trim(path) // dirsep // "POTCAR"
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
       allocate(seed(nseed)%is(10),seed(nseed)%atname(10))
       do i = 1, seed(nseed)%nspc
          ok = isinteger(nn,line,lp)
          if (.not.ok) then
             errmsg = "Too many atom types"
             exit
          end if
          do j = seed(nseed)%nat+1, seed(nseed)%nat+nn
             if (j > size(seed(nseed)%is)) then
                call realloc(seed(nseed)%is,2*(seed(nseed)%nat+nn))
                call realloc(seed(nseed)%atname,2*(seed(nseed)%nat+nn))
             end if
             seed(nseed)%is(j) = i
             seed(nseed)%atname(j) = seed(nseed)%spc(i)%name
          end do
          seed(nseed)%nat = seed(nseed)%nat + nn
       end do
       allocate(seed(nseed)%x(3,seed(nseed)%nat))
       call realloc(seed(nseed)%is,seed(nseed)%nat)
       call realloc(seed(nseed)%atname,seed(nseed)%nat)

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
       seed(i)%isformat = isformat_r_vasp
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
    use param, only: tab, isformat_r_cif
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
    integer :: cols(7) ! 1=symbol, 2=x, 3=y, 4=z, 5=label, 6=symop, 7=occupancy
    integer :: nat, nop
    real*8, allocatable :: xat(:,:), occat(:)
    integer, allocatable :: zat(:)
    character*10, allocatable :: atname(:)
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
    allocate(xat(3,20),zat(20),atname(20),occat(20))
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
             ! finalize the previous seed if we have all the info
             if (indata) then
                call fill_seed(seed)
                if (present(nseed).and.present(mseed)) then
                   if (len(errmsg) == 0) then
                      nseed = nseed + 1
                      if (nseed > size(mseed,1)) call realloc_crystalseed(mseed,2*nseed)
                      mseed(nseed) = seed
                   else
                      ! keep reading
                      errmsg = ""
                      call seed%end()
                   end if
                elseif (present(seed0)) then
                   if (len(errmsg) == 0) then
                      seed0 = seed
                      return
                   else
                      ! keep reading
                      errmsg = ""
                      call seed%end()
                   end if
                else
                   errmsg = "Error reading cif file (incorrect use of the interface)"
                   return
                end if
             end if
             errmsg = ""
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
                word == "_atom_site_fract_z" .or. word == "_atom_site_occupancy") then
                ! an atom loop
                looptype = 2
                nat = 0
                if (word == "_atom_site_type_symbol") then
                   cols(1) = nhead
                elseif (word == "_atom_site_label") then
                   if (cols(1) == 0) cols(1) = nhead
                   cols(5) = nhead
                elseif (word == "_atom_site_fract_x") then
                   cols(2) = nhead
                elseif (word == "_atom_site_fract_y") then
                   cols(3) = nhead
                elseif (word == "_atom_site_fract_z") then
                   cols(4) = nhead
                elseif (word == "_atom_site_occupancy") then
                   cols(7) = nhead
                end if
             elseif (word == "_symmetry_equiv_pos_as_xyz" .or. &
                word == "_space_group_symop_operation_xyz") then
                ! a symmetry loop
                looptype = 1
                nop = 0
                if (word == "_symmetry_equiv_pos_as_xyz") then
                   cols(6) = nhead
                elseif (word == "_space_group_symop_operation_xyz") then
                   if (cols(6) == 0) cols(6) = nhead
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
                call realloc(atname,2*nat)
                call realloc(occat,2*nat)
             end if
             occat(nat) = 1d0

             ! read the fields
             lp = 1
             do i = 1, max(maxval(cols(1:5)),cols(7))
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
                   if (cols(5)==i) &
                      atname(nat) = word
                   if (cols(7)==i) then
                      ! occupancy is optional data: keep 1 if unknown (?, .) or unparseable
                      if (isreal(rdum,word)) occat(nat) = rdum
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
             do i = 1, cols(6)
                ok = isexpression_or_word(word,line,lp)
                if (ok.and.cols(6)==i) ops(nop) = word
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
       if (len(errmsg) == 0) then
          nseed = nseed + 1
          call realloc_crystalseed(mseed,nseed)
          mseed(nseed) = seed
       end if
       if (nseed == 0) then
          errmsg = "Error reading cif file: no valid data blocks were found"
       else
          errmsg = ""
       end if
    elseif (present(seed0).and.len(errmsg) == 0) then
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
      use param, only: eyet, eye
      type(crystalseed), intent(out) :: seed

      integer :: hnum
      integer :: i, j
      real*8 :: rot0(3,4)
      logical :: ok
      type(rawseed) :: rec

      ! check we have all the info
      if (.not.all(havefields)) then
         errmsg = "Error reading cif file (missing fields or data block not found)"
         return
      end if

      ! cell (angstrom), atoms (fractional), species sorted by Z, occupancies
      rec%spcmode = spc_zsorted
      rec%ndim = 3
      rec%cellmode = 1
      rec%lunit = lunit_ang
      rec%aa = aa
      rec%bb = bb
      do i = 1, nat
         call rawseed_add_atom(rec,zat(i),xat(:,i),atname(i),occ=occat(i))
      end do
      rec%name = trim(blockname)

      ! symmetry from the symops in the cif file, else from the space group label
      rec%symmode = sym_ops
      rec%neqv = 0
      rec%ncv = 1
      allocate(rec%rotm(3,4,48),rec%cen(3,4))
      rec%cen = 0d0
      if (nop > 0) then
         do j = 1, nop
            rot0 = string_to_symop(ops(j),errmsg)
            if (len_trim(errmsg) > 0) return

            if (all(abs(eyet - rot0) < 1d-12)) then
               ! the identity
               rec%neqv = rec%neqv + 1
               if (rec%neqv > size(rec%rotm,3)) &
                  call realloc(rec%rotm,3,4,2*rec%neqv)
               rec%rotm(:,:,rec%neqv) = rot0
            elseif (all(abs(eye - rot0(1:3,1:3)) < 1d-12)) then
               ! a non-zero pure translation
               ! check if I have it already
               ok = .true.
               do i = 1, rec%ncv
                  if (all(abs(rot0(:,4) - rec%cen(:,i)) < 1d-12)) then
                     ok = .false.
                     exit
                  endif
               end do
               if (ok) then
                  rec%ncv = rec%ncv + 1
                  if (rec%ncv > size(rec%cen,2)) call realloc(rec%cen,3,2*rec%ncv)
                  rec%cen(:,rec%ncv) = rot0(:,4)
               endif
            else
               ! a rotation, with some pure translation in it
               ! check if I have this rotation matrix already
               ok = .true.
               do i = 1, rec%neqv
                  if (all(abs(rec%rotm(1:3,1:3,i) - rot0(1:3,1:3)) < 1d-12)) then
                     ok = .false.
                     exit
                  endif
               end do
               if (ok) then
                  rec%neqv = rec%neqv + 1
                  rec%rotm(:,:,rec%neqv) = rot0
               endif
            endif
         end do
      else
         hnum = spg_get_hall_number_from_symbol(spg)
         if (hnum > 0) &
            call spg_get_symmetry_from_database(hnum,rec%neqv,rec%ncv,rec%rotm,rec%cen)
      end if

      call rawseed_to_seed(rec,seed,mol,file,isformat_r_cif,errmsg)

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
    use tools_io, only: fopen_read, fclose, getline, zatguess, equal, isinteger
    use param, only: maxzat, isformat_r_mol2
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
    integer :: nat, i, idum
    integer :: zat
    real*8 :: x(3)
    character*10 :: attyp, atname
    type(crystalseed) :: seed
    type(rawseed) :: rec

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

                ! fill the seed with information from this ATOM block (angstrom)
                rec = rawseed()
                rec%spcmode = spc_z
                rec%ndim = 0
                rec%lunit = lunit_ang
                rec%iscart = .true.
                rec%ismol = .true.
                rec%border = rborder_def
                if (len_trim(molname) > 0) rec%name = trim(molname)
                do i = 1, nat
                   ok = getline(lu,line)
                   read(line,*,err=999,end=999) idum, atname, x, attyp
                   zat = zatguess(attyp)
                   if (zat < 1 .or. zat > maxzat) then
                      errmsg = "error reading atom type: " // trim(attyp)
                      goto 999
                   end if
                   call rawseed_add_atom(rec,zat,x,atname,spcname=trim(attyp))
                end do
                call rawseed_to_seed(rec,seed,.true.,file,isformat_r_mol2,errmsg)
                if (len_trim(errmsg) > 0) goto 999
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

  !> Read multiple structure seeds from a sdf file. If mol, force a
  !> molecule/crystal system. If error, return non-empty errmsg.  If
  !> nseed and mseed are present, read all seeds into the mseed array
  !> and return the number of seeds in nseed. If seed0 is present,
  !> return the molecule number id. If id is not present or is empty,
  !> return the first molecule in seed0.
  subroutine read_all_sdf(file,errmsg,nseed,mseed,seed0,id,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, getline_raw, lower, isexpression_or_word,&
       isreal, isinteger, zatguess, lgetword
    use param, only: isformat_r_sdf
    integer, intent(out), optional :: nseed
    type(crystalseed), intent(inout), allocatable, optional :: mseed(:)
    type(crystalseed), intent(inout), optional :: seed0
    integer, intent(in), optional :: id
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp
    logical :: ok, havedollars
    character(len=:), allocatable :: line, word1, lword1, word2
    character*5 :: version
    integer :: i, idseed, idtarget, idum, nat
    integer :: zat
    real*8 :: x(3)
    type(crystalseed) :: seed
    type(rawseed) :: rec

    ! consistency check
    if (.not.(present(nseed).and.present(mseed)).and..not.present(seed0)) then
       errmsg = "Error reading mol/sdf file (incorrect use of the interface)"
       return
    end if

    ! initialize
    errmsg = ""
    if (present(nseed).and.present(mseed)) then
       nseed = 0
       if (allocated(mseed)) deallocate(mseed)
       allocate(mseed(10))
    else
       call seed0%end()
    end if
    idtarget = 1
    if (present(id)) then
       if (id > 0) idtarget = id
    end if

    ! main loop
    idseed = 0
    lu = fopen_read(file,ti=ti)
    main: do while (.true.)
       ! prepare the new seed
       call seed%end()
       rec = rawseed()
       rec%spcmode = spc_z
       rec%ndim = 0
       rec%lunit = lunit_ang
       rec%iscart = .true.
       rec%ismol = .true.
       rec%border = rborder_def
       havedollars = .false.
       nat = 0

       ! header
       if (.not.getline_raw(lu,line)) exit main
       if (len_trim(line) > 0) rec%name = trim(line)
       if (.not.getline_raw(lu,line)) exit main
       if (.not.getline_raw(lu,line)) exit main

       ! counts line
       ! 11 fields w 3 char, 1 field w 6 char
       if (.not.getline_raw(lu,line)) exit main
       if (len(line) >= 39) then
          version = lower(line(35:39))
       else
          version = ""
       end if
       if (version == "v2000") then
          ! read the number of atoms
          ok = isinteger(nat,line(1:3))
          if (.not.ok) goto 100

          ! read nat lines, write down coordinates
          do i = 1, nat
             if (.not.getline_raw(lu,line)) goto 100
             if (len(line) < 33) goto 100

             ! coordinates
             ok = isreal(x(1),line(1:10))
             ok = ok .and. isreal(x(2),line(11:20))
             ok = ok .and. isreal(x(3),line(21:30))
             if (.not.ok) goto 100

             ! atomic symbol and chemical species
             zat = zatguess(line(31:33))
             if (zat < 0) goto 100
             call rawseed_add_atom(rec,zat,x,line(31:33),spcname=trim(adjustl(line(31:33))))
          end do
       elseif (version == "v3000") then
          do while (.true.)
             if (.not.getline_raw(lu,line)) goto 100
             if (line(1:4) == "$$$$") then
                havedollars = .true.
                goto 100
             end if
             if (len(line) < 8) continue
             if (lower(line(1:7)) == "m  v30 ") then
                line = line(8:)
                lp = 1
                word1 = lgetword(line,lp)
                word2 = lgetword(line,lp)
                if (word1 == "begin" .and. word2 == "ctab") then
                   ! read the number of atoms
                   if (.not.getline_raw(lu,line)) goto 100
                   if (len(line) < 8) goto 100
                   line = line(8:)
                   lp = 1
                   word1 = lgetword(line,lp)

                   ! read the value
                   ok = isinteger(nat,line,lp)
                   if (.not.ok) goto 100
                elseif (word1 == "begin" .and. word2 == "atom") then
                   ! must have the number of atoms already
                   if (nat == 0) goto 100

                   ! read the atomic coordinates
                   do i = 1, nat
                      lp = 1
                      if (.not.getline_raw(lu,line)) goto 100
                      if (len(line) < 8) goto 100
                      line = line(8:)

                      ! index and atomic symbol
                      lp = 1
                      ok = isinteger(idum,line,lp)
                      ok = ok .and. isexpression_or_word(word1,line,lp)
                      if (.not.ok) goto 100

                      ! coordinates
                      ok = isreal(x(1),line,lp)
                      ok = ok .and. isreal(x(2),line,lp)
                      ok = ok .and. isreal(x(3),line,lp)
                      if (.not.ok) goto 100

                      ! atomic symbol and chemical species
                      lword1 = lower(word1)
                      if (lword1 == "a".or.lword1 == "q".or.lword1 == "*".or.lword1 == "r#") then
                         zat = 0
                      else
                         zat = zatguess(word1)
                      end if
                      if (zat < 0) goto 100
                      call rawseed_add_atom(rec,zat,x,word1,spcname=trim(adjustl(word1)))
                   end do
                   exit
                end if
             end if
          end do
       else
          errmsg = "unknown sdf/mol version"
          goto 999
       end if

       ! build the seed (mol coordinates in angstrom)
       if (rec%nat > 0) then
          call rawseed_to_seed(rec,seed,.true.,file,isformat_r_sdf,errmsg)
          if (len_trim(errmsg) > 0) goto 999
       end if

100    continue

       ! search for the next $$$$, this completes the structure
       ! chemdraw doesn't write the $$$$ at the end (?)
       if (.not.havedollars) then
          do while (getline_raw(lu,line))
             if (len(line) >= 4) then
                if (line(1:4) == "$$$$") then
                   havedollars = .true.
                   exit
                end if
             end if
          end do
       end if
       if (.not.havedollars.and.seed%nat==0 .or. seed%nspc==0) exit

       ! adopt the seed
       idseed = idseed + 1
       if (present(nseed).and.present(mseed)) then
          nseed = nseed + 1
          call realloc_crystalseed(mseed,nseed)
          mseed(nseed) = seed
       elseif (present(seed0) .and. idseed == idtarget) then
          seed0 = seed
          exit main
       end if
    end do main
    call fclose(lu)

    ! check if a seed has been read
    errmsg = "No structures found in the sdf/mol file"
    if (present(nseed).and.present(mseed)) then
       if (nseed == 0) goto 999
    elseif (present(seed0)) then
       if (.not.seed0%isused) goto 999
    end if

    ! all done
    errmsg = ""

    return
999 continue
    call fclose(lu)
    if (present(nseed)) nseed = 0
    if (present(mseed)) deallocate(mseed)

  end subroutine read_all_sdf

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
    use param, only: bohrtoa, isformat_r_qeout
    use types, only: species
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: istruct !< ID of the structure
    logical, intent(in) :: mol !< Is this a molecule?
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ideq, i, j, nee, nstr, nrec, nstruct
    character(len=:), allocatable :: line, enechar
    character*10 :: atn, sdum
    character*40 :: sene
    integer :: idum, idx
    real*8 :: alat, r(3,3), qaux, rfac, cfac, rdum
    logical :: ok, tox
    ! interim copy of seed info
    integer :: nat, nspc
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    character*10, allocatable :: atname(:)
    type(species), allocatable :: spc(:) !< Species
    real*8 :: m_x2c(3,3)
    logical :: hasx, hasis, hasspc, hasr, keep, lastkept
    type(rawseed), allocatable :: rec(:)
    type(rawseed) :: tmpl

    nseed = 0
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! first pass: count the structures and find out whether the energies
    ! are printed as "!" or "!!" lines (hybrid functionals)
    nee = 0
    nstr = 0
    do while (getline_raw(lu,line))
       if (index(line,"!") == 1) nstr = nstr + 1
       if (index(line,"!!") == 1) nee = nee + 1 ! for hybrid functionals
    end do
    if (nee > 0) then
       nstr = nee
       enechar = "!!"
    else
       enechar = "!"
    end if
    if (nstr == 0) then
       errmsg = "No valid structures found."
       goto 999
    end if
    if (allocated(seed)) deallocate(seed)

    ! records: all of them (istruct < 0) or only the requested one (slot 1)
    if (istruct < 0) then
       allocate(rec(nstr))
    else
       allocate(rec(1))
    end if
    tmpl%spcmode = spc_given
    tmpl%ndim = 3
    tmpl%cellmode = 2
    tmpl%lunit = lunit_bohr
    nrec = 0
    nstruct = 0
    lastkept = .false.
    alat = 1d0

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
          if (allocated(atname)) deallocate(atname)
          allocate(x(3,nat),is(nat),atname(nat))

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
          if (.not.allocated(is)) then
             errmsg = "Found atomic positions before the number of atoms in the QE output."
             goto 999
          end if
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
             atname(i) = atn
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

       elseif (index(line,"CELL_PARAMETERS") == 1) then
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

       elseif (index(line,"ATOMIC_POSITIONS") == 1) then
          if (.not.allocated(is)) then
             errmsg = "Found atomic positions before the number of atoms in the QE output."
             goto 999
          end if
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
             atname(i) = atn
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
       else if (index(line,enechar) == 1) then
          if (nat == 0 .or.(.not.hasx.and.nstruct == 0)) then
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

          ! a new geometry: keep it as a record if it is wanted
          if (hasx) then
             nstruct = nstruct + 1
             keep = (istruct <= 0) .or. (istruct == nstruct)
             if (keep) then
                if (istruct < 0) then
                   nrec = nrec + 1
                   call rawseed_grow(rec,nrec)
                else
                   nrec = 1
                end if
                rec(nrec) = tmpl
                rec(nrec)%nspc = nspc
                rec(nrec)%spc = spc
                rec(nrec)%is = is
                rec(nrec)%atname = atname
                rec(nrec)%nat = nat
                rec(nrec)%x = x
                rec(nrec)%rv = m_x2c
                rec(nrec)%iscart = tox
                rec(nrec)%isfinal = (nstruct == nstr)
             end if
             lastkept = keep
          end if

          ! read the energy (Ry) into the current record
          read (line,*,err=999,end=999) sdum, sdum, sdum, sdum, sene
          read (sene,*,err=999,end=999) rdum
          if (nrec > 0 .and. lastkept) rec(nrec)%energy = rdum / 2d0
          hasx = .false.
       else if (nrec > 0 .and. lastkept .and. index(line,"total   stress") > 0) then
          ! add the pressure to the current record, if available
          idx = index(line,'=')
          ok = isreal(rec(nrec)%pressure,line(idx+1:))
          if (ok) then
             rec(nrec)%pressure = rec(nrec)%pressure / 10d0 ! kbar -> GPa
          else
             rec(nrec)%pressure = huge(1d0)
          end if
       end if
    end do
    if (nrec == 0) then
       if (istruct > 0) then
          errmsg = "Structure number not found in file: " // string(istruct)
       else
          errmsg = "No valid structures found."
       end if
       goto 999
    end if

    ! build the seeds (the requested one is in slot 1 if istruct >= 0)
    call rawseed_select(rec,nrec,min(istruct,0),mol,file,isformat_r_qeout,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  end subroutine read_all_qeout

  !> Read all structures from an xyz file. Returns nseed seeds in the
  !> seed variable. mol = 1 (force molecules), 0 (force crystals), -1
  !> (detect format). If seed0 is present, read only the first
  !> structure and ignore the rest.
  subroutine read_all_xyz(nseed,seed,file,mol,errmsg,seed0,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, getline_raw, lower, zatguess,&
       isinteger, isreal, string, nameguess
    use param, only: maxzat, isformat_r_xyz
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    integer, intent(in) :: mol
    character(len=:), allocatable, intent(out) :: errmsg
    type(crystalseed), intent(inout), optional :: seed0
    type(thread_info), intent(in), optional :: ti

    integer :: lu, nat, i, iz, lp, nrec
    logical :: ok, ismol
    real*8 :: edum, r(3,3), x(3)
    character(len=:), allocatable :: line
    character*10 :: atn
    type(rawseed), allocatable :: rec(:)
    type(rawseed) :: tmpl

    nseed = 0
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! one record per frame (angstrom, Cartesian coordinates); the title
    ! line energy is stored as read (unknown units) and is not shown
    errmsg = "Error reading file: " // trim(file)
    tmpl%spcmode = spc_label_nocase
    tmpl%lunit = lunit_ang
    tmpl%iscart = .true.
    tmpl%border = rborder_def
    allocate(rec(1))
    nrec = 0
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) cycle
       read (line,*,err=999,end=999) nat
       nrec = nrec + 1
       call rawseed_grow(rec,nrec)
       rec(nrec) = tmpl

       ! interpret the title line, if possible
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       call interpret_title_line()
       rec(nrec)%energy = edum
       if (mol == 1 .or. ismol) then
          rec(nrec)%ismol = .true.
          rec(nrec)%ndim = 0
       else
          rec(nrec)%ndim = 3
          rec(nrec)%cellmode = 2
          rec(nrec)%rv = transpose(r)
       end if

       ! the atoms (angstrom); an integer symbol is the atomic number
       do i = 1, nat
          read (lu,*,err=999,end=999) atn, x
          ok = isinteger(iz,atn)
          if (ok) then
             if (iz < 0 .or. iz > maxzat) then
                errmsg = "Invalid atomic number: "//string(iz)//"."
                goto 999
             end if
             call rawseed_add_atom(rec(nrec),iz,x,atn,spcname=trim(nameguess(iz,.true.)))
          else
             iz = zatguess(atn)
             if (iz < 0) then
                errmsg = "Unknown atomic symbol: "//trim(atn)//"."
                goto 999
             end if
             call rawseed_add_atom(rec(nrec),iz,x,atn)
          end if
       end do
       if (present(seed0)) exit ! stop here if we only need 1 seed
    end do

    ! build the seeds; several seeds are named file|i
    if (nrec > 1) then
       do i = 1, nrec
          rec(i)%name = string(i)
       end do
    end if
    call rawseed_select(rec,nrec,-1,.false.,file,isformat_r_xyz,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! assign seed0
    if (present(seed0)) seed0 = seed(1)

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  contains
    ! Interpret the title line, contained in the "line" variable.
    ! Sets edum, lp, ismol, r
    subroutine interpret_title_line()
      use tools_io, only: lower
      integer :: idx
      logical :: firstpass, ok

      ! initialize
      edum = huge(1d0)
      lp = 1
      ismol = .true.
      r = 0d0

      ! check if it is a single field -> energy
      ok = isreal(edum,line,lp)
      if (.not.ok .or. len_trim(line(lp:)) > 0) then
         edum = huge(1d0)
      else
         return
      end if

      ! check if we can read the lattice parameters from the line
      line = lower(line)
      idx = index(line,"pbc")
      if (idx > 0) then
         idx = index(line,"lattice")
         if (idx > 0) then
            ismol = .false.
            firstpass = .true.
            do i = idx+7,len(line)
               if (line(i:i) == " ") cycle ! skip blanks
               if (firstpass) then
                  if (line(i:i) /= "=") then
                     ismol = .true.
                     exit
                  else
                     firstpass = .false.
                     cycle
                  end if
               end if
               if (line(i:i) == "'" .or. line(i:i) == '"') then
                  lp = i+1
                  exit
               end if
               ismol = .true.
               exit
            end do

            ! read the real numbers
            if (.not.ismol) then
               ok = isreal(r(1,1),line,lp)
               ok = ok .and. isreal(r(1,2),line,lp)
               ok = ok .and. isreal(r(1,3),line,lp)
               ok = ok .and. isreal(r(2,1),line,lp)
               ok = ok .and. isreal(r(2,2),line,lp)
               ok = ok .and. isreal(r(2,3),line,lp)
               ok = ok .and. isreal(r(3,1),line,lp)
               ok = ok .and. isreal(r(3,2),line,lp)
               ok = ok .and. isreal(r(3,3),line,lp)
               if (.not.ok) ismol = .true.
            end if
         end if
      end if

    end subroutine interpret_title_line

  end subroutine read_all_xyz

  !> Read all structures from a Gaussian output (log) file. Returns
  !> all crystal seeds.
  subroutine read_all_log(nseed,seed,file,errmsg,alsovib,ti)
    use global, only: rborder_def
    use tools_io, only: fopen_read, fclose, getline_raw, nameguess
    use param, only: isformat_r_gaussian
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file !< Input file name
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(out), optional :: alsovib
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line
    integer :: lu, idum, iz, nrec, idx
    logical :: ok, laste, lastinputor
    real*8 :: energy, x(3)
    type(rawseed), allocatable :: rec(:)
    type(rawseed) :: tmpl

    ! initialize
    nseed = 0
    errmsg = ""
    if (present(alsovib)) alsovib = .false.

    ! open file
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! one record per orientation block; the input/Z-matrix orientation
    ! wins over the standard orientation once it has been seen
    tmpl%spcmode = spc_zsorted
    tmpl%ndim = 0
    tmpl%lunit = lunit_ang
    tmpl%iscart = .true.
    tmpl%ismol = .true.
    tmpl%border = rborder_def
    allocate(rec(1))
    nrec = 0
    laste = .false.
    lastinputor = .false.
    do while (getline_raw(lu,line))
       ok = (index(line,"Input orientation:") > 0) .or. (index(line,"Z-Matrix orientation:") > 0)
       if (ok) then
          lastinputor = .true.
       elseif (.not.lastinputor) then
          ok = (index(line,"Standard orientation:") > 0)
          if (ok) lastinputor = .false.
       end if

       ! set alsovib
       if (present(alsovib)) then
          if (.not.alsovib) then
             if (index(line," Frequencies --") > 0) alsovib = .true.
          end if
       end if

       if (ok) then
          nrec = nrec + 1
          call rawseed_grow(rec,nrec)
          rec(nrec) = tmpl
          ok = getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          ok = ok .and. getline_raw(lu,line)
          if (.not.ok) goto 999
          do while (.true.)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"---------") > 0) exit
             read(line,*,err=999,end=999) idum, iz, idum, x
             if (iz < 0) iz = 0
             call rawseed_add_atom(rec(nrec),iz,x,nameguess(iz,.true.))
          end do
          laste = .false.
       elseif (nrec > 0 .and. index(line,"SCF Done") > 0) then
          idx = index(line,"=")
          if (idx > 0) then
             line = line(idx+1:)
             read (line,*) energy
             rec(nrec)%energy = energy
          end if
          laste = .true.
       end if
    end do
    if (nrec == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    if (rec(nrec)%nat == 0) then
       errmsg = "No atoms found."
       goto 999
    end if
    ! drop a trailing geometry without an energy
    if (.not.laste .and. nrec > 1) nrec = nrec - 1
    rec(nrec)%isfinal = .true.

    ! build the seeds
    call rawseed_select(rec,nrec,-1,.true.,file,isformat_r_gaussian,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

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
    use tools_io, only: fopen_read, getline_raw, fclose, lgetword, isreal, getword, zatguess
    use param, only: hartoev, eva3togpa, isformat_r_aimsout
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, nlat, i, j, idx, ll, nrec, iz
    logical :: is_file_mol, ok, isfinal
    character*1 :: cdum
    character*10 :: dum1, dum2, splbl
    character(len=:), allocatable :: line, word
    real*8 :: rlat(3,3), x(3), rdum
    type(rawseed), allocatable :: rec(:)

    nseed = 0
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

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
    ! first record: the input geometry
    allocate(rec(1))
    nrec = 1
    call init_record(rec(1))
    ok = getline_raw(lu,line)
    ok = ok .and. getline_raw(lu,line)
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) exit
       read(line,*,err=999,end=999) cdum, dum1, dum2, splbl, (x(j),j=1,3)
       word = trim(adjustl(splbl))
       iz = zatguess(word)
       if (iz <= 0) then
          errmsg = "unknown atom type: " // word
          goto 999
       end if
       call rawseed_add_atom(rec(1),iz,x,word)
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
             ok = isreal(rdum,line(idx+1:))
             if (ok) then
                rec(nrec)%energy = rdum / hartoev
             else
                rec(nrec)%energy = huge(1d0)
             end if
          elseif (index(line,'|  Pressure') > 0) then
             idx = index(line,':')
             ok = isreal(rdum,line(idx+1:))
             if (ok) then
                rec(nrec)%pressure = rdum * eva3togpa
             else
                rec(nrec)%pressure = huge(1d0)
             end if
          elseif (trim(line) == "  Updated atomic structure:") then
             isfinal = .false.
             exit
          elseif (trim(line) == "  Final atomic structure:") then
             isfinal = .true.
             exit
          end if
       end do

       ! a new record
       nrec = nrec + 1
       call rawseed_grow(rec,nrec)
       call init_record(rec(nrec))

       ! read the geometry block
       nlat = 0
       ok = getline_raw(lu,line)
       if (.not.ok) goto 999
       do while (getline_raw(lu,line))
          ll = len(line)
          if (ll >= 25) then
             if (line == "  Fractional coordinates:") exit
          end if
          if (ll >= 4) then
             if (line(1:4) == "----") exit
          end if
          if (ll == 0) cycle
          lp = 1
          word = lgetword(line,lp)
          if (word == "lattice_vector") then
             nlat = nlat + 1
             ok = isreal(rlat(1,nlat),line,lp)
             ok = ok .and. isreal(rlat(2,nlat),line,lp)
             ok = ok .and. isreal(rlat(3,nlat),line,lp)
             if (.not.ok) goto 999
          elseif (word == "atom") then
             ok = isreal(x(1),line,lp)
             ok = ok .and. isreal(x(2),line,lp)
             ok = ok .and. isreal(x(3),line,lp)
             if (.not.ok) goto 999
             word = getword(line,lp)
             iz = zatguess(word)
             if (iz <= 0) then
                errmsg = "unknown atom type: " // word
                goto 999
             end if
             call rawseed_add_atom(rec(nrec),iz,x,word)
          end if
       end do
       if (.not.is_file_mol) rec(nrec)%rv = rlat

       ! if this is the final structure and we do not have an energy, take it from last step
       if (isfinal .and. rec(nrec)%energy == huge(1d0) .and. nrec > 1)&
          rec(nrec)%energy = rec(nrec-1)%energy
    end do main
    rec(nrec)%isfinal = .true.

    ! build the seeds
    call rawseed_select(rec,nrec,-1,is_file_mol,file,isformat_r_aimsout,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! fin
    errmsg = ""
999 continue ! error condition
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  contains
    ! common settings of every record in this file
    subroutine init_record(r)
      type(rawseed), intent(inout) :: r
      r%spcmode = spc_label
      r%lunit = lunit_ang
      r%iscart = .true.
      r%ismol = is_file_mol
      r%border = rborder_def
      if (is_file_mol) then
         r%ndim = 0
      else
         r%ndim = 3
         r%cellmode = 2
         r%rv = rlat
      end if
    end subroutine init_record
  end subroutine read_all_aimsout

  !> Read all seeds from a CASTEP geom file.
  subroutine read_all_castep_geom(nseed,seed,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, getword, isinteger, isreal, zatguess
    use param, only: isformat_r_castepgeom
    integer, intent(out) :: nseed !< number of seeds
    type(crystalseed), intent(inout), allocatable :: seed(:) !< seeds on output
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word
    integer :: lu, ll, lp, idum, iz, nrec
    logical :: ok
    real*8 :: x(3)
    type(rawseed), allocatable :: rec(:)
    type(rawseed) :: tmpl

    nseed = 0
    errmsg = ""
    ! open
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! read the geometries: one record per "<-- E" block (atomic units)
    tmpl%spcmode = spc_label_nocase
    tmpl%ndim = 3
    tmpl%cellmode = 2
    tmpl%lunit = lunit_bohr
    tmpl%iscart = .true.
    allocate(rec(1))
    nrec = 0
    do while(getline_raw(lu,line))
       line = trim(adjustl(line))
       ll = len(line)
       if (ll > 4) then
          if (line(ll-4:ll) == "<-- E") then
             nrec = nrec + 1
             call rawseed_grow(rec,nrec)
             rec(nrec) = tmpl
             read(line,*,err=999,end=999) rec(nrec)%energy
          elseif (nrec == 0) then
             cycle
          elseif (line(ll-4:ll) == "<-- h") then
             read(line,*,err=999,end=999) rec(nrec)%rv(:,1)
             if (.not.getline_raw(lu,line)) goto 999
             read(line,*,err=999,end=999) rec(nrec)%rv(:,2)
             if (.not.getline_raw(lu,line)) goto 999
             read(line,*,err=999,end=999) rec(nrec)%rv(:,3)
          elseif (line(ll-4:ll) == "<-- R") then
             lp = 1
             word = getword(line,lp)
             ok = isinteger(idum,line,lp)
             ok = ok .and. isreal(x(1),line,lp)
             ok = ok .and. isreal(x(2),line,lp)
             ok = ok .and. isreal(x(3),line,lp)
             if (.not.ok) goto 999
             if (isinteger(idum,word)) then
                iz = idum
             else
                iz = zatguess(word)
                if (iz < 0) then
                   errmsg = "Unknown atomic symbol: " // word
                   goto 999
                end if
             end if
             call rawseed_add_atom(rec(nrec),iz,x,word)
          end if
       end if
    end do
    if (nrec == 0) goto 999
    rec(nrec)%isfinal = .true.

    ! build the seeds
    call rawseed_select(rec,nrec,-1,.false.,file,isformat_r_castepgeom,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  end subroutine read_all_castep_geom

  !> Read all seeds from a crystal output file.
  subroutine read_all_crystalout(file,mol,errmsg,nseed,mseed,seed0,alsovib,ti)
    use crystalseedmod, only: crystalseed
    use tools_io, only: fopen_read, getline_raw, isinteger, isreal,&
       zatguess, fclose, equali
    use tools_math, only: matinv
    use types, only: realloc
    use param, only: bohrtoa, isformat_r_crystal
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    character(len=:), allocatable, intent(out) :: errmsg
    integer, intent(out), optional :: nseed
    type(crystalseed), intent(inout), allocatable, optional :: mseed(:)
    type(crystalseed), intent(inout), optional :: seed0
    logical, intent(out), optional :: alsovib
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, ier, idx
    character(len=:), allocatable :: line
    integer :: idum, iz, lp
    real*8 :: r(3,3), x(3), energy
    logical :: ok, firstseed, readseed, inopt
    integer :: iscrystal
    character*(10) :: ats, bts
    type(crystalseed) :: seed

    ! initialize and open the file
    call seed%end()
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    if (present(nseed) .and. present(mseed)) then
       nseed = 0
       if (allocated(mseed)) deallocate(mseed)
       allocate(mseed(10))
    end if
    if (present(alsovib)) alsovib = .false.

    ! run over the file
    errmsg = "Error reading file: " // trim(file)
    r = 0d0
    iscrystal = -1
    allocate(seed%x(3,10),seed%is(10),seed%spc(2),seed%atname(10))
    seed%nat = 0
    seed%nspc = 0
    firstseed = .true.
    readseed = .false.
    inopt = .false.
    do while (getline_raw(lu,line))
       if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle

       ! set iavib
       if (present(alsovib)) then
          if (.not.alsovib) then
             if (index(line,"NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)") > 0) &
                alsovib = .true.
          end if
       end if

       ! whether this is a crystal or molecule
       if (iscrystal < 0) then
          if (index(line,"CRYSTAL CALCULATION") > 0) then
             iscrystal = 1
          elseif (index(line,"3D - CRYSTAL") > 0) then
             iscrystal = 1
          elseif (index(line,"MOLECULAR CALCULATION") > 0) then
             iscrystal = 0
          end if
       end if

       ! reading the first and last seeds, combination of the keywords:
       ! lattice = DIRECT LATTICE VECTORS CARTESIAN COMPONENTS
       ! atoms = CARTESIAN COORDINATES - PRIMITIVE CELL
       if (index(line,"DIRECT LATTICE VECTORS CARTESIAN COMPONENTS") > 0) then
          ! read the r
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
          ! read the atomic information
          do i = 1, 3
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
          end do
          line = ""
          if (firstseed) then
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
                   call realloc(seed%atname,2*seed%nat)
                end if
                read (line,*,err=999,end=999) idum, iz, ats, x
                seed%x(:,seed%nat) = x / bohrtoa
                seed%atname(seed%nat) = ats
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
             call realloc(seed%x,3,seed%nat)
             call realloc(seed%is,seed%nat)
             call realloc(seed%atname,seed%nat)
             call realloc(seed%spc,seed%nspc)

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
             seed%isformat = isformat_r_crystal

             ! done with the first seed
             firstseed = .false.
          else
             ! this is not the first seed
             do i = 1, seed%nat
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999
                if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
                read (line,*,err=999,end=999) idum, iz, ats, x
                seed%x(:,i) = x / bohrtoa
             end do
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

          ! new seed is ready
          readseed = .true.
       elseif (index(line,"OPTOPTOPTOPT") > 0) then
          inopt = .true.
       end if

       if (inopt) then
          if (index(line,"COORDINATE AND CELL OPTIMIZATION - POINT") > 0) then
             ! read cell parameters
             do i = 1, 5
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999
                if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
             end do
             read (line,*,err=999,end=999) seed%aa, seed%bb
             seed%aa = seed%aa / bohrtoa
             seed%useabr = 1

             ! read atomic positions
             do i = 1, 4
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999
                if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
             end do
             do i = 1, seed%nat
                ok = getline_raw(lu,line)
                if (.not.ok) goto 999
                if (index(line,"PROCESS") > 0 .and. index(line,"WORKING") > 0) cycle
                read (line,*,err=999,end=999) idum, bts, iz, ats, seed%x(:,i)
                seed%x(:,i) = seed%x(:,i) - floor(seed%x(:,i))
                seed%atname(i) = ats
             end do

             ! new seed is ready
             readseed = .true.
          end if
       end if

       ! read the total energy
       if (len(line) > 13) then
          if (line(1:13) == " TOTAL ENERGY") then
             do i = 1, 3
                idx = index(line,")")
                line = line(idx+1:)
             end do
             read (line,*,err=999,end=999) energy
             if (present(seed0)) seed0%energy = energy
             if (present(nseed).and.present(mseed)) mseed(nseed)%energy = energy
          end if
       end if

       ! transfer the seed
       if (readseed) then
          if (present(seed0)) seed0 = seed
          if (present(nseed) .and. present(mseed)) then
             nseed = nseed + 1
             if (nseed > size(mseed,1)) call realloc_crystalseed(mseed,2*nseed)
             mseed(nseed) = seed
          end if
          readseed = .false.
       end if
    end do

    ! rearrange seeds if optimization
    if (present(nseed).and.present(mseed)) then
       if (inopt .and. nseed > 1) then
          ! opt = last point is repeated, move to next-to-last structure but keep the energy
          energy = mseed(nseed-1)%energy
          mseed(nseed-1) = mseed(nseed)
          mseed(nseed-1)%energy = energy
          nseed = nseed - 1
       end if
       call realloc_crystalseed(mseed,nseed)
    end if

    ! consistency checks
    if (iscrystal /= 1) then
       errmsg = "Only CRYSTAL calculations supported (no MOLECULE, SLAB or POLYMER)."
       goto 999
    end if
    if (all(r == 0d0)) then
       errmsg = "Could not find lattice vectors."
       goto 999
    end if
    if (firstseed) then
       errmsg = "Could not find crystal geometries."
       goto 999
    end if

    ! done
    errmsg = ""
    call fclose(lu)
    return

999 continue
    call fclose(lu)
    if (present(nseed) .and. present(mseed)) then
       nseed = 0
       if (allocated(mseed)) deallocate(mseed)
    end if
    if (present(seed0)) call seed0%end()

  end subroutine read_all_crystalout

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
    use param, only: isformat_r_qeout, isformat_r_crystal, isformat_r_fploout,&
       isformat_r_orca, isformat_r_aimsout
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
       if (index(line,"* O   R   C   A *") > 0) then
          isformat = isformat_r_orca
          exit
       else if (index(line,"EEEEEEEEEE STARTING  DATE") > 0) then
          isformat = isformat_r_crystal
          exit
       else if (index(line,"Program PWSCF") > 0) then
          isformat = isformat_r_qeout
          exit
       elseif (index(line,"Invoking FHI-aims ...") > 0) then
          isformat = isformat_r_aimsout
          exit
       elseif (index(line,"FULL-POTENTIAL LOCAL-ORBITAL MINIMUM BASIS BANDSTRUCTURE CODE") > 0) then
          isformat = isformat_r_fploout
          exit
       end if
    end do

    ! determine whether the aims output contains a molecule or a crystal
    if (isformat == isformat_r_aimsout) then
       isformat = 0
       do while(getline_raw(lu,line))
          if (adjustl(trim(line)) == "Input geometry:") then
             isformat = isformat_r_aimsout
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
    use param, only: isformat_r_qein, isformat_r_aimsin, isformat_r_alamode,&
       isformat_r_akaikkr
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
          isformat = isformat_r_qein
          exit
       end if
       if (equal(word,"atom").or.equal(word,"atom_frac").or.equal(word,"lattice_vector")) then
          isformat = isformat_r_aimsin
          exit
       end if
       if (equal(word,"&general")) then
          isformat = isformat_r_alamode
          exit
       end if
       if (equal(word,"go")) then
          isformat = isformat_r_akaikkr
          exit
       end if
    end do

    ! If this is an FHI-aims input, determine whether it is a molecule or a crystal
    if (isformat == isformat_r_aimsin) then
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

  !> Detect whether a Gaussian cube (binary or text, according to
  !> isbin) corresponds to a molecule or a crystal. If error, return
  !> non-zero errmsg.
  function ismol_cube(file,isbin,errmsg,ti)
    use tools_io, only: fopen_read, fclose
    use tools_math, only: matinv
    character*(*), intent(in) :: file !< Input file name
    logical, intent(in) :: isbin
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    logical :: ismol_cube

    integer :: lu, nat, nstep(3)
    integer :: i, iz, ier
    real*8 :: x0(3), rmat(3,3), rdum, rx(3)

    real*8, parameter :: eps = 1d-2

    ismol_cube = .false.
    errmsg = "Error reading file: " // trim(file)
    if (isbin) then
       lu = fopen_read(file,form="unformatted",ti=ti)
    else
       lu = fopen_read(file,ti=ti)
    end if
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! ignore the title lines and read number of atoms and unit cell
    if (isbin) then
       read (lu,err=999,end=999)
       read (lu,err=999,end=999)
       read (lu,err=999,end=999) nat, x0
    else
       read (lu,*,err=999,end=999)
       read (lu,*,err=999,end=999)
       read (lu,*,err=999,end=999) nat, x0
    end if
    nat = abs(nat)

    ! grid dimensions
    if (isbin) then
       read (lu,err=999,end=999) nstep, rmat
    else
       do i = 1, 3
          read (lu,*,err=999,end=999) nstep(i), rmat(:,i)
          rmat(:,i) = rmat(:,i) * nstep(i)
       end do
    end if
    rmat = transpose(rmat)
    call matinv(rmat,3,ier)
    if (ier /= 0) then
       errmsg = "Error inverting matrix"
       goto 999
    end if

    ! atomic positions.
    do i = 1, nat
       if (isbin) then
          read (lu,err=999,end=999) iz, rdum, rx
       else
          read (lu,*,err=999,end=999) iz, rdum, rx
       end if
       if (iz <= 0) cycle
       rx = matmul(rx - x0,rmat)
       if (any(rx < -eps) .or. any(rx > 1d0+eps)) then
          ismol_cube = .true.
          errmsg = ""
          return
       end if
    end do

    errmsg = ""
999 continue
    call fclose(lu)

  end function ismol_cube

  !> Read the namelist groups of a Quantum ESPRESSO input file (unit lu,
  !> positioned at the beginning of the file) and return the handful of
  !> variables critic2 uses. A Fortran namelist read is not used here
  !> because it fails on any variable outside the declaration, and QE
  !> inputs routinely carry keywords from newer or locally patched
  !> versions of the code; in this scan, unknown variables and unknown
  !> groups are simply ignored. On output the file is left positioned on
  !> the first line after the last namelist group, ready for the card
  !> parsing. Non-zero errmsg if no namelist group was found.
  subroutine qe_read_namelists(lu,calculation,ibrav,celldm,nat,ntyp,space_group,&
     uniqueb,rhombohedral,origin_choice,errmsg)
    use tools_io, only: getline_raw, lower, isinteger, isreal
    integer, intent(in) :: lu
    character(len=80), intent(inout) :: calculation
    integer, intent(inout) :: ibrav, nat, ntyp, space_group, origin_choice
    real*8, intent(inout) :: celldm(6)
    logical, intent(inout) :: uniqueb, rhombohedral
    character(len=:), allocatable, intent(out) :: errmsg

    character(len=*), parameter :: namechar = &
       "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"

    character(len=:), allocatable :: line, grp, word
    character*1 :: ch, quote
    integer :: iline, ip, i0, len0, nline, ngrp
    logical :: ok, ingrp

    errmsg = ""
    grp = ""
    ingrp = .false.
    nline = 0
    ngrp = 0
    iline = 0
    main: do while (getline_raw(lu,line))
       iline = iline + 1

       ! remove the trailing comment, if any (only outside quotes)
       quote = " "
       len0 = len(line)
       do ip = 1, len0
          ch = line(ip:ip)
          if (quote /= " ") then
             if (ch == quote) quote = " "
          elseif (ch == "'" .or. ch == '"') then
             quote = ch
          elseif (ch == "!" .or. ch == "#") then
             line = line(1:ip-1)
             exit
          end if
       end do
       if (len_trim(line) == 0) cycle

       ! walk the line, item by item
       ip = 1
       len0 = len(line)
       do while (ip <= len0)
          ch = line(ip:ip)
          if (ch == " " .or. ch == char(9) .or. ch == ",") then
             ip = ip + 1
          elseif (ch == "&") then
             ! a group delimiter: &name opens, &end closes
             i0 = ip + 1
             ip = i0
             do while (ip <= len0)
                if (index(namechar,line(ip:ip)) == 0) exit
                ip = ip + 1
             end do
             word = lower(line(i0:ip-1))
             if (word == "end") then
                ingrp = .false.
                nline = iline
             else
                if (.not.ingrp) ngrp = ngrp + 1
                grp = word
                ingrp = .true.
             end if
          elseif (ch == "/") then
             ! end of the group
             ingrp = .false.
             nline = iline
             ip = ip + 1
          elseif (.not.ingrp) then
             ! outside a group and not a group opener: the cards start here
             exit main
          else
             ! an assignment, up to the next unquoted comma or slash
             i0 = ip
             quote = " "
             do while (ip <= len0)
                ch = line(ip:ip)
                if (quote /= " ") then
                   if (ch == quote) quote = " "
                elseif (ch == "'" .or. ch == '"') then
                   quote = ch
                elseif (ch == "," .or. ch == "/") then
                   exit
                end if
                ip = ip + 1
             end do
             call qe_nml_assign(line(i0:ip-1))
          end if
       end do
    end do main

    if (ngrp == 0) then
       errmsg = "No namelist group (&control, &system, ...) found in the file"
       return
    end if
    if (ingrp) then
       errmsg = "Unterminated namelist group &" // trim(grp) // " in the file"
       return
    end if

    ! rewind and skip the namelist section, so the caller reads the cards
    rewind(lu)
    do ip = 1, nline
       ok = getline_raw(lu,line)
       if (.not.ok) exit
    end do

  contains
    !> Interpret one "name = value" assignment belonging to namelist
    !> group grp. Names critic2 does not use are ignored.
    subroutine qe_nml_assign(str)
      character*(*), intent(in) :: str

      character(len=:), allocatable :: nam, val
      integer :: ieq, ipar, jpar, lp, idum, idx
      real*8 :: rdum
      logical :: lok

      ieq = index(str,"=")
      if (ieq <= 1) return
      nam = lower(trim(adjustl(str(1:ieq-1))))
      val = trim(adjustl(str(ieq+1:)))
      if (len_trim(nam) == 0 .or. len_trim(val) == 0) return

      ! split off the array index, if any (celldm(1), ...)
      idx = 0
      ipar = index(nam,"(")
      if (ipar > 1) then
         jpar = index(nam,")")
         if (jpar > ipar+1) then
            lp = 1
            if (.not.isinteger(idx,nam(ipar+1:jpar-1),lp)) idx = 0
         end if
         nam = nam(1:ipar-1)
      end if

      lp = 1
      if (grp == "control") then
         if (nam == "calculation") calculation = lower(qe_nml_string(val))
      elseif (grp == "system") then
         if (nam == "ibrav") then
            if (isinteger(idum,val,lp)) ibrav = idum
         elseif (nam == "nat") then
            if (isinteger(idum,val,lp)) nat = idum
         elseif (nam == "ntyp") then
            if (isinteger(idum,val,lp)) ntyp = idum
         elseif (nam == "space_group") then
            if (isinteger(idum,val,lp)) space_group = idum
         elseif (nam == "origin_choice") then
            if (isinteger(idum,val,lp)) origin_choice = idum
         elseif (nam == "celldm") then
            if (idx >= 1 .and. idx <= 6) then
               if (isreal(rdum,val,lp)) celldm(idx) = rdum
            end if
         elseif (nam == "uniqueb") then
            if (qe_nml_logical(val,lok)) uniqueb = lok
         elseif (nam == "rhombohedral") then
            if (qe_nml_logical(val,lok)) rhombohedral = lok
         end if
      end if

    end subroutine qe_nml_assign

    !> Value of a character variable, with the surrounding quotes removed.
    function qe_nml_string(str) result(res)
      character*(*), intent(in) :: str
      character(len=:), allocatable :: res

      integer :: n

      res = trim(adjustl(str))
      n = len(res)
      if (n >= 2) then
         if ((res(1:1) == "'" .and. res(n:n) == "'") .or.&
            (res(1:1) == '"' .and. res(n:n) == '"')) res = res(2:n-1)
      end if

    end function qe_nml_string

    !> Value of a logical variable. Returns false if str is not a logical.
    function qe_nml_logical(str,val) result(isok)
      character*(*), intent(in) :: str
      logical, intent(out) :: val
      logical :: isok

      character(len=:), allocatable :: s

      isok = .false.
      val = .false.
      s = lower(trim(adjustl(str)))
      if (len_trim(s) == 0) return
      if (s(1:1) == ".") s = s(2:)
      if (len_trim(s) == 0) return
      if (s(1:1) == "t") then
         val = .true.
         isok = .true.
      elseif (s(1:1) == "f") then
         isok = .true.
      end if

    end function qe_nml_logical

  end subroutine qe_read_namelists

  !> Reader for GULP input (gin, grs).  Read one or all configurations
  !> from a GULP input file. If istruct < 0, return all
  !> configurations; if istruct = 0, return the last one; otherwise,
  !> return configuration number istruct.
  subroutine read_all_gulpin(nseed,seed,file,mol,istruct,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lower, string, equal
    use param, only: bohrtoa, isformat_r_gulpin
    use global, only: rborder_def
    integer, intent(out) :: nseed
    type(crystalseed), intent(inout), allocatable :: seed(:)
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    integer, intent(in) :: istruct
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ncfg, nword, nfloat, i, j, region, icopy, nskip, z, itype
    character(len=:), allocatable :: line, pendline, keywords, w, lw, pname
    character*(mxtok*4) :: words(mxtok)
    real*8 :: floats(mxtok), unitfac, scale, x(3), vals(6), occ
    logical :: pending, lflags, nowrap, lcelllast, ok, iscart, newcfg, iscore, havekw
    type(rawseed), allocatable :: st(:)
    ! pending cell (applies to the next coordinate block); pcellmode = 0 means none
    integer :: pcellmode
    real*8 :: paa(3), pbb(3), prv(3,3)

    nseed = 0
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! first pass: keyword line plus any "keyword" option lines
    keywords = ""
    havekw = .false.
    do while (getline_raw(lu,line))
       if (len_trim(line) == 0) cycle
       if (line(1:1) == "#") cycle
       if (.not.havekw) then
          keywords = lower(trim(line))
          havekw = .true.
       else
          lw = lower(adjustl(line))
          if (.not.gulp_stem(lw,"keyw")) cycle
          call gulp_lex(line,nword,words,nfloat,floats)
          do i = 2, nword
             keywords = keywords // " " // lower(trim(words(i)))
          end do
       end if
    end do
    if (.not.havekw) then
       errmsg = "Empty GULP input file: " // trim(file)
       goto 999
    end if
    call gulp_keyword_flags(keywords,lflags,nowrap)

    ! second pass: read the structures
    rewind(lu)
    havekw = .false.
    ncfg = 0
    allocate(st(1))
    pending = .false.
    pcellmode = 0
    paa = 0d0
    pbb = 0d0
    prv = 0d0
    lcelllast = .false.
    scale = 1d0
    pname = ""
    do while (next_line(line))
       if (.not.havekw) then
          ! skip the keyword line
          if (len_trim(line) == 0) cycle
          if (line(1:1) == "#") cycle
          havekw = .true.
          cycle
       end if
       call gulp_lex(line,nword,words,nfloat,floats)
       if (nword == 0) cycle
       w = lower(trim(words(1)))
       unitfac = 1d0
       do i = 2, nword
          if (equal(lower(trim(words(i))),"au")) unitfac = bohrtoa
       end do

       if (gulp_stem(w,"titl")) then
          ! title block: n lines or until "end"
          if (nfloat >= 1) then
             nskip = nint(floats(1))
             do i = 1, nskip
                if (.not.next_line(line)) exit
             end do
          else
             do while (next_line(line))
                lw = lower(adjustl(line))
                if (gulp_stem(lw,"end")) exit
             end do
          end if
       elseif (gulp_stem(w,"keyw")) then
          cycle
       elseif (gulp_stem(w,"igno")) then
          ! ignore ... erongi
          do while (next_line(line))
             lw = lower(adjustl(line))
             if (gulp_stem(lw,"eron")) exit
          end do
       elseif (gulp_stem(w,"star") .or. equal(w,"run") .or. equal(w,"stop")) then
          exit
       elseif (gulp_stem(w,"incl")) then
          errmsg = "GULP include files are not supported: " // trim(line)
          goto 999
       elseif (gulp_stem(w,"name")) then
          if (nword >= 2) then
             pname = trim(words(2))
          else
             if (next_line(line)) then
                call gulp_lex(line,nword,words,nfloat,floats)
                if (nword >= 1) pname = trim(words(1))
             end if
          end if
       elseif (equal(w,"cell")) then
          if (nfloat >= 6) then
             vals = floats(1:6)
          else
             ok = .false.
             do while (next_line(line))
                call gulp_lex(line,nword,words,nfloat,floats)
                if (nfloat == 0) cycle
                if (nfloat < 6) then
                   errmsg = "Missing data in GULP cell specification"
                   goto 999
                end if
                vals = floats(1:6)
                ok = .true.
                exit
             end do
             if (.not.ok) goto 999
          end if
          paa = abs(vals(1:3)) * unitfac * scale
          pbb = abs(vals(4:6))
          pcellmode = 1
          lcelllast = .true.
       elseif (gulp_stem(w,"vect")) then
          do i = 1, 3
             if (.not.getline_raw(lu,line)) goto 999
             read (line,*,err=999,end=999) prv(:,i)
          end do
          prv = prv * unitfac * scale
          if (lflags) then
             if (.not.getline_raw(lu,line)) goto 999
          end if
          pcellmode = 2
          lcelllast = .true.
       elseif (gulp_stem(w,"scal")) then
          if (nfloat >= 1) then
             scale = floats(1)
          else
             if (.not.next_line(line)) goto 999
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nfloat < 1) goto 999
             scale = floats(1)
          end if
       elseif (gulp_stem(w,"frac") .or. gulp_stem(w,"cart")) then
          ! coordinate block
          iscart = gulp_stem(w,"cart")
          region = 1
          do i = 2, nword
             if (gulp_stem(lower(trim(words(i))),"regi")) then
                if (nfloat >= 1) region = nint(floats(1))
             end if
          end do
          newcfg = .true.
          if (iscart .and. region > 1 .and. .not.lcelllast .and. ncfg > 0) newcfg = .false.
          if (newcfg) then
             ncfg = ncfg + 1
             call rawseed_grow(st,ncfg)
             if (pcellmode /= 0) then
                st(ncfg)%ndim = 3
                st(ncfg)%cellmode = pcellmode
                st(ncfg)%aa = paa
                st(ncfg)%bb = pbb
                st(ncfg)%rv = prv
             else
                if (.not.iscart) then
                   errmsg = "GULP fractional coordinates given without a cell"
                   goto 999
                end if
                st(ncfg)%ndim = 0
             end if
             st(ncfg)%iscart = iscart
             st(ncfg)%tag = pname
             pname = ""
             pcellmode = 0
          else
             if (st(ncfg)%iscart .neqv. iscart) then
                errmsg = "Mixed fractional/Cartesian regions in GULP input"
                goto 999
             end if
          end if
          lcelllast = .false.

          ! read the atoms
          do while (next_line(line))
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nword + nfloat == 0) cycle
             if (nword == 0 .and. nfloat == 1) cycle
             if (nword > 0) then
                if (.not.gulp_symbol(words(1),z,itype)) then
                   ! end of block: re-dispatch this line
                   pending = .true.
                   pendline = line
                   exit
                end if
             end if
             call gulp_atom_line(nword,words,nfloat,floats,lflags,z,itype,x,occ,iscore,ok)
             if (.not.ok) then
                errmsg = "Missing data in GULP coordinate input: " // trim(line)
                goto 999
             end if
             if (.not.iscore .or. z <= 0) cycle
             if (iscart) x = x * unitfac * scale
             call gulp_add_atom(st(ncfg),z,itype,x,occ,ok)
             if (.not.ok) then
                errmsg = "Invalid site occupancy in GULP coordinate input: " // trim(line)
                goto 999
             end if
          end do
       elseif (gulp_stem(w,"spac")) then
          if (ncfg == 0) then
             errmsg = "GULP space group given before any structure"
             goto 999
          end if
          ok = .false.
          do while (next_line(line))
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nword + nfloat == 0) cycle
             if (nword > 0) then
                st(ncfg)%spg = trim(gulp_strip_comment(line))
                st(ncfg)%spg = trim(adjustl(st(ncfg)%spg(1:min(16,len(st(ncfg)%spg)))))
             elseif (nfloat == 1) then
                st(ncfg)%spg = string(nint(floats(1)))
             else
                errmsg = "Unsupported GULP space group specification: " // trim(line)
                goto 999
             end if
             ok = .true.
             exit
          end do
          if (.not.ok) goto 999
       elseif (gulp_stem(w,"orig")) then
          if (ncfg == 0) then
             errmsg = "GULP origin given before any structure"
             goto 999
          end if
          if (nfloat == 0) then
             if (.not.next_line(line)) goto 999
             call gulp_lex(line,nword,words,nfloat,floats)
          end if
          if (nfloat == 1) then
             st(ncfg)%origin = nint(floats(1))
             if (st(ncfg)%origin < 1 .or. st(ncfg)%origin > 2) then
                errmsg = "Unknown GULP origin setting: " // trim(line)
                goto 999
             end if
          elseif (nfloat == 3) then
             st(ncfg)%haveshift = .true.
             do j = 1, 3
                if (floats(j) < 1d0) then
                   st(ncfg)%shift(j) = nint(24d0 * floats(j)) / 24d0
                else
                   st(ncfg)%shift(j) = nint(floats(j)) / 24d0
                end if
             end do
          else
             errmsg = "Unknown GULP origin specification: " // trim(line)
             goto 999
          end if
       elseif (gulp_stem(w,"ditt")) then
          ! duplicate a previous configuration
          if (ncfg == 0) then
             errmsg = "GULP ditto given before any structure"
             goto 999
          end if
          icopy = ncfg
          if (nfloat >= 1) icopy = nint(floats(1))
          if (icopy < 1 .or. icopy > ncfg) then
             errmsg = "Invalid configuration in GULP ditto"
             goto 999
          end if
          ncfg = ncfg + 1
          call rawseed_grow(st,ncfg)
          st(ncfg) = st(icopy)
          if (len_trim(pname) > 0) st(ncfg)%tag = pname
          pname = ""
          pcellmode = 0
          lcelllast = .false.
       elseif (gulp_stem(w,"scel") .or. gulp_stem(w,"svec") .or. gulp_stem(w,"sfra") .or.&
          gulp_stem(w,"pcel") .or. gulp_stem(w,"pvec") .or. gulp_stem(w,"pfra")) then
          errmsg = "GULP surface/polymer (2D/1D) structures are not supported: " // trim(words(1))
          goto 999
       elseif (gulp_stem(w,"cont") .or. gulp_stem(w,"supe") .or. gulp_stem(w,"symmetry_o")) then
          errmsg = "Unsupported GULP option: " // trim(words(1))
          goto 999
       end if
    end do
    if (ncfg == 0) then
       errmsg = "No structures found in GULP input file: " // trim(file)
       goto 999
    end if

    ! conversion policies
    do i = 1, ncfg
       st(i)%symmode = sym_spg
       st(i)%wrap = .not.nowrap .and. .not.st(i)%iscart ! GULP wraps fractional input only
       if (st(i)%ndim == 0) then
          st(i)%ismol = .true.
          st(i)%border = rborder_def
       end if
    end do

    ! build the seeds
    call rawseed_select(st,ncfg,istruct,mol,file,isformat_r_gulpin,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  contains
    ! get the next logical line: honor pushback and & continuation, strip comments
    function next_line(lineo) result(oko)
      character(len=:), allocatable, intent(inout) :: lineo
      logical :: oko

      character(len=:), allocatable :: aux
      integer :: iidx

      if (pending) then
         lineo = pendline
         pending = .false.
         oko = .true.
         return
      end if
      oko = getline_raw(lu,lineo)
      if (.not.oko) return
      lineo = gulp_strip_comment(lineo)
      do while (index(lineo,"&") > 0)
         iidx = index(lineo,"&")
         lineo = lineo(1:iidx-1)
         if (.not.getline_raw(lu,aux)) exit
         lineo = lineo // " " // gulp_strip_comment(aux)
      end do

    end function next_line

  end subroutine read_all_gulpin

  !> Read one or all geometries from a GULP output file. The
  !> geometries are the input echo of every configuration followed by
  !> the final geometry of every completed optimisation. If istruct <
  !> 0, return all of them; if istruct = 0, return the last one;
  !> otherwise, return geometry number istruct.
  subroutine read_all_gulpout(nseed,seed,file,mol,istruct,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lower, isreal, isinteger
    use param, only: isformat_r_gulpout, eva3togpa, hartoev, hartokjmol, kcal2ha
    use global, only: rborder_def
    integer, intent(out) :: nseed
    type(crystalseed), intent(inout), allocatable :: seed(:)
    character*(*), intent(in) :: file
    logical, intent(in) :: mol
    integer, intent(in) :: istruct
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, nst, ncfg, icur, ifin, idx, i, j, nword, nfloat
    character(len=:), allocatable :: line, lw
    character*(mxtok*4) :: words(mxtok)
    real*8 :: floats(mxtok), rdum
    logical :: have_e0, want_prim, ok
    real*8 :: efin
    type(rawseed), allocatable :: st(:)

    nseed = 0
    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    errmsg = "Error reading file: " // trim(file)

    ! All the input echoes come before the first output block, so the
    ! echo (initial geometry) of configuration k is always st(k).
    nst = 0
    ncfg = 0
    icur = 0
    ifin = 0
    have_e0 = .true.
    want_prim = .false.
    efin = huge(1d0)
    allocate(st(1))
    do while (getline_raw(lu,line))
       if (index(line,"*  Input for Configuration =") > 0) then
          ! new configuration: initial geometry record
          ncfg = ncfg + 1
          nst = nst + 1
          call rawseed_grow(st,nst)
          st(nst)%ndim = 3
       elseif (ncfg == 0) then
          cycle
       elseif (index(line,"Dimensionality =") > 0) then
          idx = index(line,"=")
          ok = isinteger(st(ncfg)%ndim,line(idx+1:))
          if (.not.ok) goto 999
          if (st(ncfg)%ndim == 1 .or. st(ncfg)%ndim == 2) then
             errmsg = "GULP surface/polymer (2D/1D) structures are not supported"
             goto 999
          end if
       elseif (index(line,"Space group (") > 0) then
          idx = index(line,":")
          st(ncfg)%spg = trim(adjustl(line(idx+1:)))
       elseif (index(line,"Shift of the origin") > 0) then
          idx = index(line,":")
          call gulp_lex(line(idx+1:),nword,words,nfloat,floats)
          if (nfloat /= 3) goto 999
          st(ncfg)%haveshift = .true.
          st(ncfg)%shiftfromout = .true.
          st(ncfg)%shift = floats(1:3)
       elseif (index(line,"Cartesian lattice vectors (Angstroms) :") > 0 .and.&
          index(line,"Final") == 0 .and. index(line,"Reference") == 0) then
          call read_vectors(st(ncfg)%rv)
          if (st(ncfg)%cellmode == 0) st(ncfg)%cellmode = 2
       elseif (index(line,"Primitive cell parameters :") > 0 .and.&
          index(line,"Full cell parameters :") > 0) then
          st(ncfg)%gulp%centred = .true.
          do i = 1, 3
             call next_nonblank(line,ok)
             if (.not.ok) goto 999
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nfloat < 4) goto 999
             st(ncfg)%aa(i) = floats(3)
             st(ncfg)%bb(i) = floats(4)
          end do
          st(ncfg)%cellmode = 1
       elseif (index(line,"Pressure of configuration =") > 0) then
          idx = index(line,"=")
          ok = isreal(st(ncfg)%gulp%pappl,line(idx+1:))
          if (.not.ok) st(ncfg)%gulp%pappl = huge(1d0)
       elseif (index(line,"Fractional coordinates of asymmetric unit :") > 0 .or.&
          index(line,"Cartesian coordinates of cluster :") > 0) then
          st(ncfg)%iscart = (index(line,"cluster") > 0)
          call read_table(st(ncfg),.true.,ok)
          if (.not.ok) goto 999
       elseif (index(line,"Mixed fractional/Cartesian coordinates") > 0) then
          errmsg = "GULP surface/polymer (2D/1D) structures are not supported"
          goto 999
       elseif (index(line,"*  Output for configuration") > 0) then
          idx = index(line,"configuration") + len("configuration")
          ok = isinteger(icur,line(idx:))
          if (.not.ok .or. icur < 1 .or. icur > ncfg) then
             errmsg = "Invalid configuration number in GULP output"
             goto 999
          end if
          have_e0 = .false.
          want_prim = .false.
          efin = huge(1d0)
          ifin = 0
       elseif (index(line,"#  Output for point") > 0 .or. index(line,"#  Output for cell point") > 0 .or.&
          index(line,"structures optimised out of") > 0) then
          efin = huge(1d0)
          ifin = 0
       elseif (icur > 0 .and. .not.have_e0 .and. (index(line,"Total lattice energy") > 0 .or.&
          index(line,"Total lattice enthalpy") > 0 .or. index(line,"Total free energy") > 0)) then
          idx = index(line,"=")
          if (idx > 0) then
             call read_energy(line(idx+1:),rdum,ok)
             if (ok) st(icur)%energy = rdum / hartoev
             have_e0 = .true.
          else
             want_prim = .true.
          end if
       elseif (want_prim .and. index(line,"Primitive unit cell") > 0 .and. index(line,"=") > 0) then
          idx = index(line,"=")
          call read_energy(line(idx+1:),rdum,ok)
          if (ok) st(icur)%energy = rdum / hartoev
          have_e0 = .true.
          want_prim = .false.
       elseif (index(line,"Final energy =") > 0 .or. index(line,"Final enthalpy =") > 0 .or.&
          index(line,"Final free energy =") > 0) then
          idx = index(line,"=")
          call read_energy(line(idx+1:),rdum,ok)
          if (ok) then
             efin = rdum
          else
             efin = huge(1d0)
          end if
       elseif (icur > 0 .and. (index(line,"Final asymmetric unit coordinates :") > 0 .or.&
          index(line,"Final fractional coordinates of atoms :") > 0 .or.&
          index(line,"Final cartesian coordinates of atoms :") > 0)) then
          ! final geometry record
          nst = nst + 1
          call rawseed_grow(st,nst)
          st(nst) = st(icur) ! inherits the cell (fixed-cell runs), pappl, spg, occupancies
          st(nst)%isfinal = .true.
          st(nst)%energy = huge(1d0)
          if (efin /= huge(1d0)) st(nst)%energy = efin / hartoev
          st(nst)%nat = 0
          st(nst)%gulp%havedede = .false.
          st(nst)%gulp%havestress = .false.
          st(nst)%gulp%havefinalcell = .false.
          st(nst)%gulp%vol = 0d0
          call read_table(st(nst),.false.,ok)
          if (.not.ok) goto 999
          ifin = nst
          efin = huge(1d0)
       elseif (index(line,"Final fractional/Cartesian coordinates") > 0) then
          errmsg = "GULP surface/polymer (2D/1D) structures are not supported"
          goto 999
       elseif (ifin > 0 .and. index(line,"Final Cartesian lattice vectors (Angstroms) :") > 0) then
          call read_vectors(st(ifin)%rv)
          if (.not.st(ifin)%gulp%centred) then
             st(ifin)%cellmode = 2
             st(ifin)%gulp%havefinalcell = .true.
          end if
       elseif (ifin > 0 .and. index(line,"Final cell parameters") > 0) then
          ! a, b, c, alpha, beta, gamma rows, with optional strain derivatives
          call skip_rule(ok)
          if (.not.ok) goto 999
          do while (getline_raw(lu,line))
             if (gulp_is_rule(line)) exit
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nword < 1 .or. nfloat < 1) cycle
             lw = lower(trim(words(1)))
             select case (lw)
             case ("a")
                j = 1
             case ("b")
                j = 2
             case ("c")
                j = 3
             case ("alpha")
                j = 4
             case ("beta")
                j = 5
             case ("gamma")
                j = 6
             case default
                j = 0
             end select
             if (j == 0) cycle
             if (.not.st(ifin)%gulp%centred .and. .not.st(ifin)%gulp%havefinalcell) then
                if (j <= 3) then
                   st(ifin)%aa(j) = floats(1)
                else
                   st(ifin)%bb(j-3) = floats(1)
                end if
             end if
             if (j <= 3 .and. nfloat >= 2 .and. index(line,"dE/de") > 0) then
                st(ifin)%gulp%dede(j) = floats(2)
                st(ifin)%gulp%havedede = .true.
             end if
          end do
          if (.not.st(ifin)%gulp%centred .and. .not.st(ifin)%gulp%havefinalcell) then
             st(ifin)%cellmode = 1
             st(ifin)%gulp%havefinalcell = .true.
          end if
       elseif (ifin > 0 .and. index(line,"Primitive cell volume =") > 0) then
          idx = index(line,"=")
          ok = isreal(st(ifin)%gulp%vol,line(idx+1:))
          if (.not.ok) st(ifin)%gulp%vol = 0d0
       elseif (ifin > 0 .and. index(line,"Non-primitive lattice parameters :") > 0) then
          if (st(ifin)%gulp%centred) then
             do i = 1, 2
                call next_nonblank(line,ok)
                if (.not.ok) goto 999
                call gulp_lex(line,nword,words,nfloat,floats)
                if (nfloat < 3) goto 999
                if (i == 1) then
                   st(ifin)%aa = floats(1:3)
                else
                   st(ifin)%bb = floats(1:3)
                end if
             end do
             st(ifin)%cellmode = 1
             st(ifin)%gulp%havefinalcell = .true.
          end if
       elseif (ifin > 0 .and. index(line,"Final stress tensor components (GPa):") > 0) then
          call skip_rule(ok)
          if (.not.ok) goto 999
          do i = 1, 3
             call next_nonblank(line,ok)
             if (.not.ok) goto 999
             call gulp_lex(line,nword,words,nfloat,floats)
             if (nfloat < 1) goto 999
             st(ifin)%gulp%stress(i) = floats(1)
          end do
          st(ifin)%gulp%havestress = .true.
       elseif (index(line,"Job Finished at") > 0) then
          exit
       end if
    end do
    if (nst == 0) then
       errmsg = "No structures found in GULP output file: " // trim(file)
       goto 999
    end if

    ! pressures
    do i = 1, nst
       if (st(i)%ndim /= 3) then
          st(i)%pressure = huge(1d0)
       elseif (st(i)%gulp%havestress) then
          st(i)%pressure = -sum(st(i)%gulp%stress) / 3d0
       elseif (st(i)%gulp%havedede .and. st(i)%gulp%vol > 0d0 .and. st(i)%gulp%pappl /= huge(1d0)) then
          st(i)%pressure = st(i)%gulp%pappl - sum(st(i)%gulp%dede) / (3d0 * st(i)%gulp%vol) * eva3togpa
       else
          st(i)%pressure = st(i)%gulp%pappl
       end if
       st(i)%symmode = sym_spg
       if (st(i)%ndim == 0) then
          st(i)%ismol = .true.
          st(i)%border = rborder_def
       end if
    end do

    ! build the seeds
    call rawseed_select(st,nst,istruct,mol,file,isformat_r_gulpout,nseed,seed,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    errmsg = ""
999 continue
    call fclose(lu)
    if (len_trim(errmsg) > 0) then
       nseed = 0
       if (allocated(seed)) deallocate(seed)
    end if

  contains
    ! read the next non-blank line
    subroutine next_nonblank(lineo,oko)
      character(len=:), allocatable, intent(inout) :: lineo
      logical, intent(out) :: oko

      oko = .false.
      do while (getline_raw(lu,lineo))
         if (len_trim(lineo) > 0) then
            oko = .true.
            return
         end if
      end do

    end subroutine next_nonblank

    ! skip up to and including the next dashed rule
    subroutine skip_rule(oko)
      logical, intent(out) :: oko
      character(len=:), allocatable :: aux

      oko = .false.
      do while (getline_raw(lu,aux))
         if (gulp_is_rule(aux)) then
            oko = .true.
            return
         end if
      end do

    end subroutine skip_rule

    ! read three lattice vectors (rows) after a blank line
    subroutine read_vectors(rvo)
      real*8, intent(out) :: rvo(3,3)
      logical :: oko
      integer :: ii

      rvo = 0d0
      do ii = 1, 3
         call next_nonblank(line,oko)
         if (.not.oko) return
         read (line,*,err=10,end=10) rvo(:,ii)
      end do
10    continue

    end subroutine read_vectors

    ! read an energy value followed by a unit word; convert to eV
    subroutine read_energy(str,ene,oko)
      character*(*), intent(in) :: str
      real*8, intent(out) :: ene
      logical, intent(out) :: oko

      character*20 :: unit
      integer :: ios

      oko = .false.
      ene = huge(1d0)
      unit = ""
      read (str,*,iostat=ios) ene, unit
      if (ios /= 0) then
         read (str,*,iostat=ios) ene
         if (ios /= 0) return
      end if
      unit = lower(unit)
      if (unit(1:2) == "kj") then
         ene = ene * hartoev / hartokjmol
      elseif (unit(1:4) == "kcal") then
         ene = ene * kcal2ha * hartoev
      end if
      oko = .true.

    end subroutine read_energy

    ! read a coordinate table (input echo if echo=.true., final table otherwise)
    subroutine read_table(s,echo,oks)
      type(rawseed), intent(inout) :: s
      logical, intent(in) :: echo
      logical, intent(out) :: oks

      character*5 :: label
      character*2 :: ctype
      real*8 :: x(3), occ
      real*8, allocatable :: occ0(:)
      integer :: z, itype, ios, nn, nat0
      logical :: oko

      ! header: rule, two lines, rule
      oks = .false.
      call skip_rule(oko)
      if (.not.oko) return
      call skip_rule(oko)
      if (.not.oko) return

      ! keep the occupancies of the initial geometry (final tables do not list them)
      nat0 = 0
      if (.not.echo .and. allocated(s%occ)) then
         nat0 = size(s%occ,1)
         occ0 = s%occ
      end if

      s%nat = 0
      do while (getline_raw(lu,line))
         if (len_trim(line) == 0) exit
         if (gulp_is_rule(line)) cycle
         if (len(line) < 52) cycle
         read (line(1:7),*,iostat=ios) nn
         if (ios /= 0) cycle
         label = line(9:13)
         ctype = line(15:16)
         occ = 1d0
         if (echo) then
            read (line(19:27),*,iostat=ios) x(1)
            if (ios /= 0) cycle
            read (line(31:39),*,iostat=ios) x(2)
            if (ios /= 0) cycle
            read (line(43:51),*,iostat=ios) x(3)
            if (ios /= 0) cycle
            if (len(line) >= 75) then
               read (line(64:75),*,iostat=ios) occ
               if (ios /= 0) occ = 1d0
            end if
         else
            read (line(17:52),*,iostat=ios) x
            if (ios /= 0) cycle
         end if
         if (index(ctype,"s") > 0) cycle ! shells
         if (.not.gulp_symbol(label,z,itype)) then
            errmsg = "Unknown atom label in GULP output: " // trim(label)
            return
         end if
         if (z <= 0) cycle ! dummy atoms
         call gulp_add_atom(s,z,itype,x,occ,oks)
         if (.not.oks) then
            errmsg = "Invalid site occupancy in GULP output: " // trim(line)
            return
         end if
      end do
      if (.not.echo .and. nat0 == s%nat .and. s%nat > 0) s%occ(1:s%nat) = occ0(1:s%nat)
      oks = .true.

    end subroutine read_table

  end subroutine read_all_gulpout

  !> Detect whether the first structure in a GULP input or output
  !> file is a molecule (0D cluster).
  subroutine gulp_detect_ismol(file,isformat,ismol,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lower, isinteger
    use param, only: isformat_r_gulpin
    character*(*), intent(in) :: file
    integer, intent(in) :: isformat
    logical, intent(out) :: ismol
    type(thread_info), intent(in), optional :: ti

    integer :: lu, nword, nfloat, idx, ndim
    character(len=:), allocatable :: line, w
    character*(mxtok*4) :: words(mxtok)
    real*8 :: floats(mxtok)
    logical :: havecell, havekw, ok

    ismol = .false.
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) return

    if (isformat == isformat_r_gulpin) then
       havecell = .false.
       havekw = .false.
       do while (getline_raw(lu,line))
          if (len_trim(line) == 0) cycle
          if (line(1:1) == "#") cycle
          if (.not.havekw) then
             havekw = .true.
             cycle
          end if
          call gulp_lex(line,nword,words,nfloat,floats)
          if (nword == 0) cycle
          w = lower(trim(words(1)))
          if (w == "cell" .or. gulp_stem(w,"vect") .or. gulp_stem(w,"scel") .or.&
             gulp_stem(w,"svec") .or. gulp_stem(w,"pcel") .or. gulp_stem(w,"pvec")) then
             havecell = .true.
          elseif (gulp_stem(w,"frac") .or. gulp_stem(w,"sfra") .or. gulp_stem(w,"pfra")) then
             exit
          elseif (gulp_stem(w,"cart")) then
             ismol = .not.havecell
             exit
          end if
       end do
    else
       do while (getline_raw(lu,line))
          if (index(line,"Dimensionality =") > 0) then
             idx = index(line,"=")
             ok = isinteger(ndim,line(idx+1:))
             if (ok) ismol = (ndim == 0)
             exit
          end if
       end do
    end if
    call fclose(lu)

  end subroutine gulp_detect_ismol

  ! GULP reader helpers

  !> Strip a GULP comment (everything from the first #) from a line.
  function gulp_strip_comment(line) result(res)
    character*(*), intent(in) :: line
    character(len=:), allocatable :: res

    integer :: idx

    idx = index(line,"#")
    if (idx > 0) then
       res = line(1:idx-1)
    else
       res = line
    end if

  end function gulp_strip_comment

  !> True if the line is a dashed rule (starts with ----).
  function gulp_is_rule(line) result(ok)
    character*(*), intent(in) :: line
    logical :: ok

    character(len=:), allocatable :: aux

    aux = trim(adjustl(line))
    ok = .false.
    if (len(aux) < 4) return
    ok = (aux(1:4) == "----")

  end function gulp_is_rule

  !> GULP-style prefix match: word starts with stem.
  function gulp_stem(word,stem) result(ok)
    character*(*), intent(in) :: word, stem
    logical :: ok

    ok = .false.
    if (len_trim(word) < len(stem)) return
    ok = (word(1:len(stem)) == stem)

  end function gulp_stem

  !> Tokenize a GULP line into words and numbers (reproduces linepro).
  !> A token is a number if its first character (or the second, when
  !> the first is a period) is a digit or a sign and it contains fewer
  !> than two letters; a/b fractions are accepted.
  subroutine gulp_lex(line0,nword,words,nfloat,floats)
    use tools_io, only: isletter
    character*(*), intent(in) :: line0
    integer, intent(out) :: nword, nfloat
    character*(*), intent(out) :: words(:)
    real*8, intent(out) :: floats(:)

    character(len=:), allocatable :: line, tok
    integer :: i, j, n, nlet, idx, ios
    character*1 :: c
    real*8 :: rnum, rden
    logical :: isnum

    nword = 0
    nfloat = 0
    line = gulp_strip_comment(line0)
    n = len_trim(line)
    i = 1
    do while (i <= n)
       ! skip separators
       if (line(i:i) == " " .or. line(i:i) == achar(9)) then
          i = i + 1
          cycle
       end if
       j = i
       do while (j <= n)
          if (line(j:j) == " " .or. line(j:j) == achar(9)) exit
          j = j + 1
       end do
       tok = line(i:j-1)
       i = j

       ! classify
       isnum = .false.
       c = tok(1:1)
       if (c == "." .and. len(tok) > 1) c = tok(2:2)
       if (index("0123456789+-",c) > 0) then
          nlet = 0
          do j = 1, len(tok)
             if (isletter(tok(j:j))) nlet = nlet + 1
          end do
          isnum = (nlet < 2)
       end if
       if (isnum) then
          idx = index(tok,"/")
          if (idx > 0) then
             read (tok(1:idx-1),*,iostat=ios) rnum
             if (ios == 0) read (tok(idx+1:),*,iostat=ios) rden
             if (ios == 0 .and. rden /= 0d0) then
                rnum = rnum / rden
             else
                isnum = .false.
             end if
          else
             read (tok,*,iostat=ios) rnum
             if (ios /= 0) isnum = .false.
          end if
       end if
       if (isnum) then
          if (nfloat < size(floats)) then
             nfloat = nfloat + 1
             floats(nfloat) = rnum
          end if
       else
          if (nword < size(words)) then
             nword = nword + 1
             words(nword) = tok
          end if
       end if
    end do

  end subroutine gulp_lex

  !> Atomic number from a 1- or 2-character GULP element symbol
  !> (case-insensitive). D is deuterium (1); X is a dummy atom (0).
  !> Returns -1 if unknown.
  function gulp_element(sym) result(z)
    use tools_io, only: zatguess, nameguess, equali
    character*(*), intent(in) :: sym
    integer :: z

    character*2 :: s

    s = sym
    z = -1
    if (equali(trim(s),"D")) then
       z = 1
    elseif (equali(trim(s),"X")) then
       z = 0
    else
       z = zatguess(s)
       if (z > 0) then
          ! require an exact match (zatguess is lenient with the second character)
          if (.not.equali(trim(nameguess(z,.true.)),trim(s))) z = -1
       else
          z = -1
       end if
    end if

  end function gulp_element

  !> Decide whether a word is a GULP atom label and split it into
  !> element (z) and type number (itype). Reproduces worsy and ltont:
  !> the label (up to a blank or underscore) has 1 to 5 characters,
  !> its third character is not a third letter, its first n characters
  !> (n = number of letters in the label, at most 2) are an element
  !> symbol, and the digits after the symbol are the type.
  function gulp_symbol(word,z,itype) result(ok)
    use tools_io, only: isletter, isdigit
    character*(*), intent(in) :: word
    integer, intent(out) :: z, itype
    logical :: ok

    character(len=:), allocatable :: aux
    integer :: ilen, idx, nbeg, nlet, i

    ok = .false.
    z = -1
    itype = 0
    aux = trim(adjustl(word))
    idx = index(aux,"_")
    if (idx > 0) then
       ilen = idx - 1
    else
       ilen = len_trim(aux)
    end if
    if (ilen < 1 .or. ilen > 5) return
    nlet = 0
    do i = 1, ilen
       if (isletter(aux(i:i))) then
          nlet = nlet + 1
          if (i == 3 .and. nlet == 3) return
       end if
    end do
    if (nlet == 0) return
    z = gulp_element(aux(1:min(nlet,2)))
    if (z < 0) return
    nbeg = 2
    if (ilen >= 2) then
       if (isletter(aux(2:2))) nbeg = 3
    end if

    ! trailing type digits
    itype = 0
    do i = nbeg, ilen
       if (.not.isdigit(aux(i:i))) exit
       itype = itype * 10 + (ichar(aux(i:i)) - ichar("0"))
    end do
    ok = .true.

  end function gulp_symbol

  !> Interpret the numbers of a GULP coordinate line (reproduces
  !> strword): optional trailing flags (only if lflags), numeric label
  !> form (Z, or Z+100 for a shell), silent removal of three extra
  !> numbers, and x y z [q [occ [rad]]]. iscore is false for shells.
  subroutine gulp_atom_line(nword,words,nfloat0,floats,lflags,z,itype,x,occ,iscore,ok)
    use tools_io, only: lower
    integer, intent(in) :: nword
    character*(*), intent(in) :: words(:)
    integer, intent(in) :: nfloat0
    real*8, intent(in) :: floats(:)
    logical, intent(in) :: lflags
    integer, intent(inout) :: z
    integer, intent(inout) :: itype
    real*8, intent(out) :: x(3)
    real*8, intent(out) :: occ
    logical, intent(out) :: iscore
    logical, intent(out) :: ok

    integer :: nfloat, nbeg, nai
    character*4 :: w2

    ok = .false.
    iscore = .true.
    x = 0d0
    occ = 1d0
    nfloat = nfloat0
    if (lflags) nfloat = nfloat - 3
    if (nword > 0) then
       ! label given by the caller (z, itype); type word
       nbeg = 0
       if (nword >= 2) then
          w2 = lower(words(2)(1:4))
          if (w2(1:1) == "s" .or. w2(1:2) == "bs" .or. w2(1:2) == "qs" .or. w2(1:3) == "qbs") &
             iscore = .false.
       end if
    else
       if (nfloat < 1) return
       nai = nint(floats(1))
       if (nai > 100) then
          nai = nai - 100
          iscore = .false.
       end if
       z = nai
       itype = 0
       nbeg = 1
       nfloat = nfloat - 1
    end if
    if (nfloat > 6) nfloat = nfloat - 3
    if (nfloat < 3 .or. nfloat > 6) return
    x = floats(nbeg+1:nbeg+3)
    if (nfloat >= 5) occ = floats(nbeg+5)
    ok = .true.

  end subroutine gulp_atom_line

  !> Interpret the keyword line: whether optimisation flags follow
  !> the coordinates and the cell (lflags), and whether fractional
  !> coordinates are left unwrapped (nowrap).
  subroutine gulp_keyword_flags(keywords,lflags,nowrap)
    use tools_io, only: lgetword, equal
    character*(*), intent(in) :: keywords
    logical, intent(out) :: lflags, nowrap

    integer :: lp
    character(len=:), allocatable :: w
    logical :: lopt, lnoflag

    lopt = .false.
    lnoflag = .false.
    nowrap = .false.
    lp = 1
    do while (.true.)
       w = lgetword(keywords,lp)
       if (len_trim(w) == 0) exit
       if (gulp_stem(w,"opti") .or. gulp_stem(w,"grad") .or. gulp_stem(w,"harm") .or.&
          gulp_stem(w,"fit") .or. equal(w,"rfo") .or. equal(w,"mc") .or. equal(w,"md") .or.&
          gulp_stem(w,"neb") .or. gulp_stem(w,"sync")) lopt = .true.
       if (gulp_stem(w,"conp") .or. gulp_stem(w,"conv") .or. gulp_stem(w,"cell") .or.&
          gulp_stem(w,"noflag") .or. equal(w,"shell")) lnoflag = .true.
       if (gulp_stem(w,"nomod") .or. gulp_stem(w,"nowrap")) nowrap = .true.
    end do
    lflags = lopt .and. .not.lnoflag

  end subroutine gulp_keyword_flags

  !> Find the Hall numbers for a GULP space group specification (HM
  !> symbol or ITA number): hnum is the standard setting (origin
  !> choice 1; hexagonal or rhombohedral axes for R-lattice groups,
  !> decided from the cell angles bb) and hnum2 the second origin
  !> choice, or -1 if the group has only one.
  subroutine gulp_hall_number(spg,bb,hnum,hnum2,errmsg)
    use spglib, only: spg_get_hall_number_from_symbol, spg_get_hall_number_from_number,&
       spg_get_spacegroup_type, SpglibSpaceGroupType
    use tools_io, only: isinteger, equal
    character*(*), intent(in) :: spg
    real*8, intent(in) :: bb(3)
    integer, intent(out) :: hnum, hnum2
    character(len=:), allocatable, intent(out) :: errmsg

    type(SpglibSpaceGroupType) :: sa
    character(len=:), allocatable :: sym, want
    character*11 :: shsym
    integer :: lp, ispg, h
    logical :: ishex
    character*1 :: lat

    real*8, parameter :: tol = 1d-3

    errmsg = ""
    hnum = -1
    hnum2 = -1
    sym = trim(adjustl(spg))
    if (len_trim(sym) == 0) return

    lp = 1
    if (isinteger(ispg,sym,lp)) then
       if (ispg < 1 .or. ispg > 230) then
          errmsg = "Invalid GULP space group number: " // trim(sym)
          return
       end if
       hnum = spg_get_hall_number_from_number(ispg,"")
    else
       ! old-style cubic symbols without the bar (F M 3 M) are in the spglib mapping
       hnum = spg_get_hall_number_from_symbol(sym)
    end if
    if (hnum <= 0) then
       errmsg = "Unknown GULP space group: " // trim(sym)
       return
    end if
    sa = spg_get_spacegroup_type(hnum)

    ! rhombohedral groups: hexagonal or rhombohedral axes
    shsym = adjustl(sa%international_short)
    lat = shsym(1:1)
    if (lat == "R") then
       ishex = (abs(bb(1)-90d0) < tol .and. abs(bb(2)-90d0) < tol .and. abs(bb(3)-120d0) < tol)
       if (ishex) then
          want = "H"
       else
          want = "R"
       end if
       if (.not.equal(trim(sa%choice),want)) then
          h = spg_get_hall_number_from_number(sa%number,want)
          if (h <= 0) then
             errmsg = "Setting not found for GULP space group: " // trim(sym)
             return
          end if
          hnum = h
       end if
    else
       hnum2 = spg_get_hall_number_from_number(sa%number,"2")
    end if

  end subroutine gulp_hall_number

  !> True if two sets of symmetry operations (rotations + centering
  !> vectors) are the same group, up to lattice translations.
  function gulp_same_group(n1,nc1,rot1,cen1,n2,nc2,rot2,cen2) result(same)
    integer, intent(in) :: n1, nc1, n2, nc2
    real*8, intent(in) :: rot1(:,:,:), cen1(:,:), rot2(:,:,:), cen2(:,:)
    logical :: same

    integer :: i, j, k, l
    real*8 :: t(3), d(3)
    logical :: found

    real*8, parameter :: eps = 1d-5

    same = .false.
    if (n1 * nc1 /= n2 * nc2) return
    do i = 1, n1
       do j = 1, nc1
          t = rot1(:,4,i) + cen1(:,j)
          found = .false.
          do k = 1, n2
             if (any(abs(rot1(1:3,1:3,i) - rot2(1:3,1:3,k)) > eps)) cycle
             do l = 1, nc2
                d = t - rot2(:,4,k) - cen2(:,l)
                d = d - nint(d)
                if (all(abs(d) < eps)) then
                   found = .true.
                   exit
                end if
             end do
             if (found) exit
          end do
          if (.not.found) return
       end do
    end do
    same = .true.

  end function gulp_same_group

  !> Apply GULP's origin shift convention to a set of operations:
  !> t' = t + R*s - s.
  subroutine gulp_apply_shift(neqv,rotm,s)
    integer, intent(in) :: neqv
    real*8, intent(inout) :: rotm(:,:,:)
    real*8, intent(in) :: s(3)

    integer :: i

    do i = 1, neqv
       rotm(:,4,i) = rotm(:,4,i) + matmul(rotm(1:3,1:3,i),s) - s
       rotm(:,4,i) = rotm(:,4,i) - floor(rotm(:,4,i))
    end do

  end subroutine gulp_apply_shift

  !> Convert an interim structure record into a crystal seed. The
  !> record fields select the policies: species (spcmode), length
  !> units (lunit), molecule/crystal (ismol), and symmetry (symmode).
  !> See the rawseed type for the details.
  subroutine rawseed_to_seed(rec,seed,mol,file,isformat,errmsg)
    use tools_math, only: matinv, m_x2c_from_cellpar
    use tools_io, only: nameguess, string, equal, equali
    use types, only: realloc
    use param, only: bohrtoa, maxzat
    type(rawseed), intent(in) :: rec
    type(crystalseed), intent(inout) :: seed
    logical, intent(in) :: mol
    character*(*), intent(in) :: file
    integer, intent(in) :: isformat
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, j, ier, nat
    integer :: usedz(0:maxzat)
    character*10 :: key
    real*8 :: r(3,3), lfac
    logical :: found, cart, match

    errmsg = ""
    call seed%end()
    nat = rec%nat
    if (nat == 0) then
       errmsg = "No atoms found."
       return
    end if

    ! atoms
    seed%nat = nat
    allocate(seed%x(3,nat),seed%is(nat),seed%atname(nat))
    seed%x = rec%x(:,1:nat)
    if (allocated(rec%atname)) seed%atname = rec%atname(1:nat)
    if (rec%spcmode == spc_z .or. rec%spcmode == spc_zsorted) then
       do i = 1, nat
          if (rec%z(i) < 0 .or. rec%z(i) > maxzat) then
             errmsg = "Invalid atomic number: " // string(rec%z(i))
             return
          end if
       end do
    end if

    ! species
    select case (rec%spcmode)
    case (spc_given)
       seed%nspc = rec%nspc
       seed%spc = rec%spc(1:rec%nspc)
       seed%is = rec%is(1:nat)
       if (.not.allocated(rec%atname)) then
          do i = 1, nat
             seed%atname(i) = seed%spc(seed%is(i))%name
          end do
       end if
    case (spc_label,spc_label_nocase)
       ! first appearance of the species key (the label unless spcname is given)
       allocate(seed%spc(nat))
       seed%nspc = 0
       do i = 1, nat
          if (allocated(rec%spcname)) then
             key = rec%spcname(i)
          else
             key = rec%atname(i)
          end if
          seed%is(i) = 0
          do j = 1, seed%nspc
             if (rec%spcmode == spc_label) then
                match = equal(seed%spc(j)%name,key)
             else
                match = equali(seed%spc(j)%name,key)
             end if
             if (match) then
                seed%is(i) = j
                exit
             end if
          end do
          if (seed%is(i) == 0) then
             seed%nspc = seed%nspc + 1
             seed%spc(seed%nspc)%name = key
             seed%spc(seed%nspc)%z = rec%z(i)
             seed%is(i) = seed%nspc
          end if
       end do
       call realloc(seed%spc,seed%nspc)
    case (spc_z)
       ! first appearance of the atomic number, named by the species key
       allocate(seed%spc(nat))
       seed%nspc = 0
       usedz = 0
       do i = 1, nat
          if (usedz(rec%z(i)) == 0) then
             seed%nspc = seed%nspc + 1
             if (allocated(rec%spcname)) then
                seed%spc(seed%nspc)%name = rec%spcname(i)
             else
                seed%spc(seed%nspc)%name = rec%atname(i)
             end if
             seed%spc(seed%nspc)%z = rec%z(i)
             usedz(rec%z(i)) = seed%nspc
          end if
          seed%is(i) = usedz(rec%z(i))
       end do
       call realloc(seed%spc,seed%nspc)
    case (spc_zsorted)
       ! atomic numbers present, in increasing Z
       usedz = 0
       do i = 1, nat
          usedz(rec%z(i)) = 1
       end do
       seed%nspc = count(usedz > 0)
       allocate(seed%spc(seed%nspc))
       seed%nspc = 0
       do i = 0, maxzat
          if (usedz(i) == 0) cycle
          seed%nspc = seed%nspc + 1
          seed%spc(seed%nspc)%name = nameguess(i,.true.)
          seed%spc(seed%nspc)%z = i
          usedz(i) = seed%nspc
       end do
       do i = 1, nat
          seed%is(i) = usedz(rec%z(i))
       end do
    case default
       errmsg = "Unknown species policy"
       return
    end select

    ! occupancies (kept only if some site is partial)
    if (allocated(rec%occ)) call seed_set_occ(seed,rec%occ,nat)

    ! cell and coordinates
    lfac = 1d0
    if (rec%lunit == lunit_ang) lfac = bohrtoa
    if (rec%ndim == 0) then
       seed%useabr = 0
       seed%m_x2c = 0d0
       do i = 1, nat
          if (.not.atomisfrac(i)) seed%x(:,i) = seed%x(:,i) / lfac
       end do
    else
       if (rec%cellmode == 1) then
          seed%useabr = 1
          seed%aa = rec%aa / lfac
          seed%bb = rec%bb
          if (any(seed%aa <= 0d0) .or. any(seed%bb <= 0d0) .or. any(seed%bb >= 180d0)) then
             errmsg = "Invalid cell parameters"
             return
          end if
          r = m_x2c_from_cellpar(seed%aa,seed%bb,ier)
          if (ier /= 0) then
             errmsg = "Invalid cell parameters"
             return
          end if
       elseif (rec%cellmode == 2) then
          seed%useabr = 2
          seed%m_x2c = rec%rv / lfac
          r = seed%m_x2c
       else
          errmsg = "Missing cell"
          return
       end if
       cart = .false.
       do i = 1, nat
          if (.not.atomisfrac(i)) cart = .true.
       end do
       if (cart) then
          call matinv(r,3,ier)
          if (ier /= 0) then
             errmsg = "Error inverting lattice vector matrix"
             return
          end if
          do i = 1, nat
             if (.not.atomisfrac(i)) seed%x(:,i) = matmul(r,seed%x(:,i) / lfac)
          end do
       end if
       if (rec%wrap) seed%x = seed%x - floor(seed%x)
    end if

    ! molecule flags
    seed%ismolecule = mol .or. rec%ismol
    seed%cubic = rec%cubic
    seed%border = rec%border
    seed%havex0 = rec%havex0
    seed%molx0 = 0d0

    ! symmetry
    seed%havesym = 0
    seed%checkrepeats = .false.
    seed%neqlist = .false.
    seed%findsym = -1
    found = .false.
    if (rec%ndim == 3) then
       select case (rec%symmode)
       case (sym_spg)
          call gulp_spg_ops(rec,seed%neqv,seed%ncv,seed%rotm,seed%cen,found,errmsg)
          if (len_trim(errmsg) > 0) return
       case (sym_ops)
          found = (rec%neqv > 0 .and. rec%ncv > 0)
          if (found) then
             seed%neqv = rec%neqv
             seed%ncv = rec%ncv
             seed%rotm = rec%rotm(:,:,1:rec%neqv)
             seed%cen = rec%cen(:,1:rec%ncv)
          end if
       end select
    end if
    if (found) then
       ! crystals check for repeats; a molecule does not apply crystallographic symmetry
       if (.not.mol) then
          seed%havesym = 1
          seed%checkrepeats = .true.
          seed%neqlist = .true.
       end if
       seed%findsym = 0
    end if

    ! properties and flags
    if (rec%energy /= huge(1d0)) seed%energy = rec%energy
    if (rec%pressure /= huge(1d0)) seed%pressure = rec%pressure
    seed%isused = .true.
    seed%file = file
    seed%name = file
    if (allocated(rec%name)) seed%name = trim(file) // "|" // rec%name
    seed%isformat = isformat

  contains
    ! whether atom i is given in fractional coordinates
    function atomisfrac(i)
      integer, intent(in) :: i
      logical :: atomisfrac
      if (allocated(rec%isfrac)) then
         atomisfrac = rec%isfrac(i)
      else
         atomisfrac = .not.rec%iscart
      end if
    end function atomisfrac
  end subroutine rawseed_to_seed

  !> Symmetry operations from a GULP space group specification (the
  !> rec%spg string with rec%origin and rec%shift). found is true if
  !> the group is not P1; P1 is treated as no symmetry because the
  !> repeat check is O(N^2) and useless.
  subroutine gulp_spg_ops(rec,neqv,ncv,rotm,cen,found,errmsg)
    use spglib, only: spg_get_symmetry_from_database
    use tools_math, only: cellpar_from_metric
    use types, only: realloc
    type(rawseed), intent(in) :: rec
    integer, intent(out) :: neqv, ncv
    real*8, allocatable, intent(inout) :: rotm(:,:,:), cen(:,:)
    logical, intent(out) :: found
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: hnum, h1, h2, n1, nc1, n2, nc2
    real*8 :: aa(3), bb(3)
    real*8, allocatable :: rot1(:,:,:), cen1(:,:), rot2(:,:,:), cen2(:,:)
    logical :: apply

    errmsg = ""
    found = .false.
    neqv = 0
    ncv = 0
    if (.not.allocated(rec%spg)) return
    if (len_trim(rec%spg) == 0) return

    if (rec%cellmode == 1) then
       bb = rec%bb
    else
       call cellpar_from_metric(matmul(transpose(rec%rv),rec%rv),aa,bb)
    end if
    call gulp_hall_number(rec%spg,bb,h1,h2,errmsg)
    if (len_trim(errmsg) > 0) return

    ! GULP's operations for a group with two origin choices are those of
    ! origin 2 (symmet.F90): origin 1 is obtained by shifting them, and an
    ! explicit shift applies on top of the origin-2 operations. Output
    ! files print the shift both for origin 2 (where it is the shift
    ! between the two settings and is not applied) and for explicit
    ! shifts (applied); tell the two cases apart by checking whether the
    ! shifted origin-2 operations are the origin-1 group.
    hnum = h1
    apply = rec%haveshift
    if (h2 > 0) then
       if (rec%haveshift) then
          hnum = h2
          if (rec%shiftfromout) then
             call spg_get_symmetry_from_database(h1,n1,nc1,rot1,cen1)
             call spg_get_symmetry_from_database(h2,n2,nc2,rot2,cen2)
             call gulp_apply_shift(n2,rot2,rec%shift)
             apply = .not.gulp_same_group(n1,nc1,rot1,cen1,n2,nc2,rot2,cen2)
          end if
       elseif (rec%origin == 2) then
          hnum = h2
       end if
    end if
    call spg_get_symmetry_from_database(hnum,neqv,ncv,rotm,cen)
    if (neqv <= 0 .or. ncv <= 0) then
       errmsg = "Could not expand space group: " // trim(rec%spg)
       return
    end if
    if (apply) call gulp_apply_shift(neqv,rotm,rec%shift)
    call realloc(rotm,3,4,neqv)
    call realloc(cen,3,ncv)
    found = (neqv * ncv > 1)

  end subroutine gulp_spg_ops

  !> Select the structures requested by istruct from an array of
  !> interim records and convert them to seeds. istruct < 0: all; 0:
  !> last; n: number n. When several seeds are returned and the record
  !> has no name of its own, it is named file|NNNNN (E eV), or
  !> file|(fin) (E eV) for the last one if it is a final geometry; the
  !> parenthesis holds the record's tag instead when it has one, and is
  !> omitted when neither is known.
  subroutine rawseed_select(recs,nrec,istruct,mol,file,isformat,nseed,seed,errmsg)
    use tools_io, only: string
    use param, only: hartoev
    type(rawseed), intent(in) :: recs(:)
    integer, intent(in) :: nrec
    integer, intent(in) :: istruct
    logical, intent(in) :: mol
    character*(*), intent(in) :: file
    integer, intent(in) :: isformat
    integer, intent(out) :: nseed
    type(crystalseed), intent(inout), allocatable :: seed(:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, iuse, npad
    character(len=:), allocatable :: str, ename

    errmsg = ""
    if (allocated(seed)) deallocate(seed)
    if (istruct < 0) then
       nseed = nrec
       allocate(seed(nseed))
       npad = ceiling(log10(nrec-1+0.1d0))
       do i = 1, nrec
          call rawseed_to_seed(recs(i),seed(i),mol,file,isformat,errmsg)
          if (len_trim(errmsg) > 0) return
          if (nrec > 1 .and. .not.allocated(recs(i)%name)) then
             ename = ""
             if (allocated(recs(i)%tag)) then
                if (len_trim(recs(i)%tag) > 0) ename = " (" // trim(recs(i)%tag) // ")"
             elseif (recs(i)%energy /= huge(1d0)) then
                ename = " (" // trim(adjustl(string(recs(i)%energy*hartoev,'f',decimal=8))) // " eV)"
             end if
             if (i == nrec .and. recs(i)%isfinal) then
                seed(i)%name = trim(file) // "|(fin)" // ename
             else
                str = string(i,npad,pad0=.true.)
                str = string(str,length=max(5,len(str)))
                seed(i)%name = trim(file) // "|" // str // ename
             end if
          end if
       end do
    else
       if (istruct == 0) then
          iuse = nrec
       else
          iuse = istruct
       end if
       if (iuse < 1 .or. iuse > nrec) then
          errmsg = "Structure number not found in file: " // string(istruct)
          return
       end if
       nseed = 1
       allocate(seed(1))
       call rawseed_to_seed(recs(iuse),seed(1),mol,file,isformat,errmsg)
    end if

  end subroutine rawseed_select

  !> Append an atom to an interim record: atomic number z,
  !> coordinates x, label atname, and optionally the occupancy occ,
  !> a species key spcname different from the label, and the
  !> fractional flag isfrac. The optional arrays are allocated the
  !> first time the argument is given.
  subroutine rawseed_add_atom(rec,z,x,atname,occ,spcname,isfrac)
    use types, only: realloc
    type(rawseed), intent(inout) :: rec
    integer, intent(in) :: z
    real*8, intent(in) :: x(3)
    character*(*), intent(in) :: atname
    real*8, intent(in), optional :: occ
    character*(*), intent(in), optional :: spcname
    logical, intent(in), optional :: isfrac

    integer :: n

    if (.not.allocated(rec%x)) then
       allocate(rec%x(3,10),rec%z(10),rec%atname(10))
       rec%nat = 0
    end if
    rec%nat = rec%nat + 1
    n = rec%nat
    if (n > size(rec%z,1)) then
       call realloc(rec%x,3,2*n)
       call realloc(rec%z,2*n)
       call realloc(rec%atname,2*n)
       if (allocated(rec%occ)) call realloc(rec%occ,2*n)
       if (allocated(rec%spcname)) call realloc(rec%spcname,2*n)
       if (allocated(rec%isfrac)) call realloc(rec%isfrac,2*n)
    end if
    rec%x(:,n) = x
    rec%z(n) = z
    rec%atname(n) = atname
    if (present(occ)) then
       if (.not.allocated(rec%occ)) then
          allocate(rec%occ(size(rec%z,1)))
          rec%occ(1:n-1) = 1d0
       end if
       rec%occ(n) = occ
    elseif (allocated(rec%occ)) then
       rec%occ(n) = 1d0
    end if
    if (present(spcname)) then
       if (.not.allocated(rec%spcname)) then
          allocate(rec%spcname(size(rec%z,1)))
          rec%spcname(1:n-1) = rec%atname(1:n-1)
       end if
       rec%spcname(n) = spcname
    elseif (allocated(rec%spcname)) then
       rec%spcname(n) = atname
    end if
    if (present(isfrac)) then
       if (.not.allocated(rec%isfrac)) then
          allocate(rec%isfrac(size(rec%z,1)))
          rec%isfrac(1:n-1) = .not.rec%iscart
       end if
       rec%isfrac(n) = isfrac
    elseif (allocated(rec%isfrac)) then
       rec%isfrac(n) = .not.rec%iscart
    end if

  end subroutine rawseed_add_atom

  !> Append a GULP atom (element z, type number itype, occupancy occ)
  !> to a record. Returns ok = .false. if the occupancy is invalid.
  subroutine gulp_add_atom(rec,z,itype,x,occ,ok)
    use tools_io, only: nameguess, string
    type(rawseed), intent(inout) :: rec
    integer, intent(in) :: z, itype
    real*8, intent(in) :: x(3)
    real*8, intent(in) :: occ
    logical, intent(out) :: ok

    ok = (occ > 0d0 .and. occ <= 1d0 + 1d-10)
    if (.not.ok) return
    if (itype > 0) then
       call rawseed_add_atom(rec,z,x,trim(nameguess(z,.true.)) // string(itype),occ=occ)
    else
       call rawseed_add_atom(rec,z,x,nameguess(z,.true.),occ=occ)
    end if

  end subroutine gulp_add_atom


  !> Grow an array of interim records to hold at least n entries.
  subroutine rawseed_grow(recs,n)
    type(rawseed), allocatable, intent(inout) :: recs(:)
    integer, intent(in) :: n

    type(rawseed), allocatable :: aux(:)

    if (n > size(recs,1)) then
       allocate(aux(2*n))
       aux(1:size(recs,1)) = recs
       call move_alloc(aux,recs)
    end if

  end subroutine rawseed_grow


end submodule proc
