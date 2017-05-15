! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
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

! Routines for basic crystallography computations
module struct_basic
  use spglib, only: SpglibDataset
  use types, only: atom, celatom, neighstar, fragment
  implicit none

  private

  ! private to this module
  private :: lattpg
  private :: typeop
  ! private for wigner-seitz routines
  private :: equiv_tetrah
  private :: perm3
  ! other crystallography tools that are crystal-independent
  public :: search_lattice
  public :: pointgroup_info
  ! crystal seed symmetry initialization
  public :: spgs_wrap
  
  !> Minimal amount of information to generate a crystal
  type crystalseed
     ! general
     logical :: isused = .false. !< Is the seed being used?
     character(len=128) :: file = "" !< Source file, if available
     ! atoms
     integer :: nat = 0 !< Number of atoms
     integer :: usezname = 0 !< 0 = uninit; 1 = use z; 2 = use name; 3 = use both
     real*8, allocatable :: x(:,:) !< Atomic positions (crystal - fractional;molecule with useabr=0 - bohr)
     integer, allocatable :: z(:) !< Atomic numbers
     character*(10), allocatable :: name(:) !< Atomic names
     ! cell
     integer :: useabr = 0 !< 0 = uninit; 1 = use aa,bb; 2 = use crys2car
     real*8 :: aa(3) !< Cell lengths (bohr)
     real*8 :: bb(3) !< Cell angles (degrees)
     real*8 :: crys2car(3,3) !< Crystallographic to cartesian matrix
     ! symmetry
     integer :: havesym = 0 !< Have symmetry? 0 = no; 1 = yes
     integer :: findsym = -1 !< Find the symmetry? 0 = no; 1 = yes; -1 = only small
     integer :: neqv = 0 !< Number of symmetry operations
     integer :: ncv = 0 !< Number ofcentering vectors
     real*8, allocatable :: cen(:,:) !< Centering vectors
     real*8, allocatable :: rotm(:,:,:) !< Space group operations
     ! molecular fields
     logical :: ismolecule = .false. !< Is this a molecule?
     logical :: cubic = .false. !< Use a cubic cell for the molecule
     real*8 :: border = 0d0 !< border of the molecular cell (bohr)
     logical :: havex0 = .false. !< an origin of the cell for molecules has bene given
     real*8 :: molx0(3) = 0d0 !< origin of the cell for molecules
  end type crystalseed
  public :: crystalseed

  !> Crystal type
  type crystal
     ! Initialization flags
     logical :: isinit = .false. !< has the crystal structure been initialized?
     logical :: isenv = .false. !< were the atomic environments determined?
     integer :: havesym = 0 !< was the symmetry determined? (0 - nosym, 1 - full)
     logical :: isast = .false. !< have the molecular asterisms and connectivity been calculated?
     logical :: isewald = .false. !< do we have the data for ewald's sum?
     logical :: isrecip = .false. !< symmetry information about the reciprocal cell
     logical :: isnn = .false. !< information about the nearest neighbors

     ! file name for the occasional critic2 trick
     character(len=128) :: file

     !! Initialization level: isinit !!
     ! non-equivalent atoms list
     integer :: nneq = 0 !< Number of non-equivalent atoms
     type(atom), allocatable :: at(:) !< Non-equivalent atom array
     ! complete atoms list
     integer :: ncel = 0 !< Number of atoms in the main cell
     type(celatom), allocatable :: atcel(:) !< List of atoms in the main cell
     ! cell and lattice metrics
     real*8 :: aa(3) !< cell lengths (bohr)
     real*8 :: bb(3) !< cell angles (degrees)
     real*8 :: omega !< unit cell volume
     real*8 :: gtensor(3,3) !< metric tensor (3,3)
     real*8 :: ar(3) !< reciprocal cell lengths
     real*8 :: grtensor(3,3) !< reciprocal metric tensor (3,3)
     ! crystallographic/cartesian conversion matrices
     real*8 :: crys2car(3,3) !< crystallographic to cartesian matrix
     real*8 :: car2crys(3,3) !< cartesian to crystallographic matrix
     real*8 :: n2_x2c !< sqrt(3)/norm-2 of the crystallographic to cartesian matrix
     real*8 :: n2_c2x !< sqrt(3)/norm-2 of the cartesian to crystallographic matrix
     ! space-group symmetry
     type(SpglibDataset) :: spg !< spglib's symmetry dataset
     integer :: neqv !< number of symmetry operations
     integer :: neqvg !< number of symmetry operations, reciprocal space
     integer :: ncv  !< number of centering vectors
     real*8, allocatable :: cen(:,:) !< centering vectors
     real*8 :: rotm(3,4,48) !< symmetry operations
     real*8 :: rotg(3,3,48) !< symmetry operations, reciprocal space
     ! variables for molecular systems
     logical :: ismolecule = .false. !< is it a molecule?
     real*8 :: molx0(3) !< centering vector for the molecule
     real*8 :: molborder(3) !< border length (cryst coords)
     ! wigner-seitz cell 
     integer :: nws !< number of WS neighbors/faces
     integer :: ivws(3,16) !< WS neighbor lattice points
     logical :: isortho !< is the cell orthogonal?
     integer :: nvert_ws !< number of vertices of the WS cell
     integer, allocatable :: nside_ws(:) !< number of sides of WS faces
     integer, allocatable :: iside_ws(:,:) !< sides of the WS faces
     real*8, allocatable :: vws(:,:) !< vertices of the WS cell

     !! Initialization level: isenv !!
     ! atomic environment of the cell
     integer :: nenv = 0 !< Environment around the main cell
     real*8 :: dmax0_env !< Maximum environment distance
     type(celatom), allocatable :: atenv(:) !< Atoms around the main cell

     !! Initialization level: isast !!
     ! asterisms
     type(neighstar), allocatable :: nstar(:) !< Neighbor stars
     integer :: nmol = 0 !< Number of molecules in the unit cell
     type(fragment), allocatable :: mol(:) !< Molecular fragments
     logical, allocatable :: moldiscrete(:) !< Is the crystal extended or molecular?

     !! Initialization level: isewald !!
     ! ewald data
     real*8 :: rcut, hcut, eta, qsum
     integer :: lrmax(3), lhmax(3)

   contains
     procedure :: init => struct_init !< Allocate arrays and nullify variables
     procedure :: checkflags => struct_checkflags !< Check the flags for a given crystal
     procedure :: end => struct_end !< Deallocate arrays and nullify variables
     procedure :: x2c !< Convert crystallographic to cartesian
     procedure :: c2x !< Convert cartesian to crystallographic
     procedure :: distance !< Distance between points in crystallographic coordinates
     procedure :: eql_distance !< Shortest distance between lattice-translated vectors
     procedure :: shortest !< Gives the lattice-translated vector with shortest length
     procedure :: are_close !< True if a vector is at a distance less than eps of another
     procedure :: are_lclose !< True if a vector is at a distance less than eps of all latice translations of another
     procedure :: nearest_atom !< Calculate the atom nearest to a given point
     procedure :: identify_atom !< Identify an atom in the unit cell
     procedure :: identify_fragment !< Build an atomic fragment of the crystal
     procedure :: identify_fragment_from_xyz !< Build a crystal fragment from an xyz file
     procedure :: symeqv  !< Calculate the symmetry-equivalent positions of a point
     procedure :: get_mult !< Multiplicity of a point
     procedure :: get_mult_reciprocal !< Reciprocal-space multiplicity of a point
     procedure :: build_env !< Build the crystal environment (atenv)
     procedure :: find_asterisms !< Find the molecular asterisms (atomic connectivity)
     procedure :: fill_molecular_fragments !< Find the molecular fragments in the crystal
     procedure :: listatoms_cells !< List all atoms in n cells (maybe w border and molmotif)
     procedure :: listatoms_sphcub !< List all atoms in a sphere or cube
     procedure :: listmolecules !< List all molecules in the crystal
     procedure :: pointshell !< Calculate atomic shells around a point
     procedure :: sitesymm !< Determine the local-symmetry group symbol for a point
     procedure :: get_pack_ratio !< Calculate the packing ratio
     procedure :: powder !< Calculate the powder diffraction pattern
     procedure :: rdf !< Calculate the radial distribution function
     procedure :: calculate_ewald_cutoffs !< Calculate the cutoffs for Ewald's sum
     procedure :: newcell !< Change the unit cell and rebuild the crystal
     procedure :: cell_standard !< Transform the the standard cell (possibly primitive)
     procedure :: cell_niggli !< Transform to the Niggli primitive cell
     procedure :: cell_delaunay !< Transform to the Delaunay primitive cell
     procedure :: delaunay_reduction !< Perform the delaunay reduction.
     procedure :: struct_new !< Initialize the structure from minimal info (passed as arguments)
     procedure :: struct_fill !< Initialize the structure from minimal info (already in the object)
     procedure :: struct_report !< Write lots of information about the crystal structure to uout
     procedure :: struct_report_symxyz !< Write sym. ops. in crystallographic notation to uout
     procedure :: spglib_wrap !< Fill symmetry information in the crystal using spglib
     procedure :: wigner !< Calculate the WS cell and the IWS/tetrahedra
     procedure :: pmwigner !< Poor man's wigner
  end type crystal
  public :: crystal

  ! the current crystal
  type(crystal), target :: cr
  public :: cr

  ! enumeration for structure file format types
  integer, parameter, public :: isformat_unknown = 0
  integer, parameter, public :: isformat_cif = 1
  integer, parameter, public :: isformat_res = 2
  integer, parameter, public :: isformat_cube = 3
  integer, parameter, public :: isformat_struct = 4
  integer, parameter, public :: isformat_abinit = 5
  integer, parameter, public :: isformat_elk = 6
  integer, parameter, public :: isformat_qein = 7
  integer, parameter, public :: isformat_qeout = 8
  integer, parameter, public :: isformat_crystal = 9
  integer, parameter, public :: isformat_xyz = 10
  integer, parameter, public :: isformat_wfn = 11
  integer, parameter, public :: isformat_wfx = 12
  integer, parameter, public :: isformat_fchk = 13
  integer, parameter, public :: isformat_molden = 14
  integer, parameter, public :: isformat_siesta = 15
  integer, parameter, public :: isformat_xsf = 16
  integer, parameter, public :: isformat_gen = 17
  integer, parameter, public :: isformat_vasp = 18

  ! symmetry operation symbols
  integer, parameter :: ident=0 !< identifier for sym. operations
  integer, parameter :: inv=1 !< identifier for sym. operations
  integer, parameter :: c2=2 !< identifier for sym. operations
  integer, parameter :: c3=3 !< identifier for sym. operations
  integer, parameter :: c4=4 !< identifier for sym. operations
  integer, parameter :: c6=5 !< identifier for sym. operations
  integer, parameter :: s3=6 !< identifier for sym. operations
  integer, parameter :: s4=7 !< identifier for sym. operations
  integer, parameter :: s6=8 !< identifier for sym. operations
  integer, parameter :: sigma=9 !< identifier for sym. operations

  ! array initialization values
  integer, parameter :: mneq0 = 4
  integer, parameter :: mcel0 = 10
  integer, parameter :: menv0 = 100

  ! holohedry identifier
  integer, parameter, public :: holo_unk = 0 ! unknown
  integer, parameter, public :: holo_tric = 1 ! triclinic
  integer, parameter, public :: holo_mono = 2 ! monoclinic
  integer, parameter, public :: holo_ortho = 3 ! orthorhombic
  integer, parameter, public :: holo_tetra = 4 ! tetragonal
  integer, parameter, public :: holo_trig = 5 ! trigonal
  integer, parameter, public :: holo_hex = 6 ! hexagonal
  integer, parameter, public :: holo_cub = 7 ! cubic
  character(len=12), parameter, public :: holo_string(0:7) = (/ &
     "unknown     ","triclinic   ","monoclinic  ","orthorhombic",&
     "tetragonal  ","trigonal    ","hexagonal   ","cubic       "/)

  ! Laue class identifier
  integer, parameter, public :: laue_unk = 0 ! unknown
  integer, parameter, public :: laue_1 = 1 ! -1
  integer, parameter, public :: laue_2m = 2 ! 2/m
  integer, parameter, public :: laue_mmm = 3 ! mmm
  integer, parameter, public :: laue_4m = 4 ! 4/m
  integer, parameter, public :: laue_4mmm = 5 ! 4/mmm
  integer, parameter, public :: laue_3 = 6 ! -3
  integer, parameter, public :: laue_3m = 7 ! -3m
  integer, parameter, public :: laue_6m = 8 ! 6/m
  integer, parameter, public :: laue_6mmm = 9 ! 6/mmm
  integer, parameter, public :: laue_m3 = 10 ! m-3
  integer, parameter, public :: laue_m3m = 11 ! m-3m
  character(len=12), parameter, public :: laue_string(0:11) = (/ &
     "unknown","-1     ","2/m    ","mmm    ","4/m    ","4/mmm  ",&
     "-3     ","-3m    ","6/m    ","6/mmm  ","m-3    ","m-3m   "/)

contains

  !xx! crystal class methods
  !> Initialize the struct arrays
  subroutine struct_init(c)
    use param, only: eyet, eye
    class(crystal), intent(inout) :: c

    integer :: i

    ! deallocate all the arrays and reinitialize
    call c%end()

    ! allocate space for atoms
    if (.not.allocated(c%at)) allocate(c%at(mneq0))
    if (.not.allocated(c%atcel)) allocate(c%atcel(mcel0))
    if (.not.allocated(c%atenv)) allocate(c%atenv(menv0))

    ! nullify metrics
    c%aa = 0d0
    c%bb = 0d0
    c%ar = 0d0
    c%omega = 0d0
    c%gtensor = 0d0
    c%grtensor = 0d0
    c%crys2car = 0d0
    c%car2crys = 0d0
    c%n2_x2c = 0d0
    c%n2_c2x = 0d0

    ! no symmetry
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,:,1) = eyet
    c%neqvg = 1
    c%rotg = 0d0
    c%rotg(:,:,1) = eye
    c%ncv = 1
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%cen = 0d0
    c%isortho = .false.
    c%nws = 0
    c%spg%n_atoms = 0

    ! no molecule
    c%ismolecule = .false.
    c%molx0 = 0d0
    c%molborder = 0d0

    ! initialize atoms
    do i = 1, mneq0
       c%at(i)%name = ""
       c%at(i)%z = 0
       c%at(i)%zpsp = -1
       c%at(i)%qat = 0d0
       c%at(i)%rnn2 = 0d0
    end do

    ! no molecular fragments
    c%nmol = 0

    ! the crystal is not initialized until struct_fill is run
    c%isinit = .false. 
    c%isenv = .false. 
    c%havesym = 0 
    c%isast = .false. 
    c%isewald = .false. 

  end subroutine struct_init

  !> Check that the flags necessary for the operation of critic2 using
  !> the given crystal structure have been set. Several flags are
  !> available: environments available (env), wigner-seitz cell
  !> calculated (ws), asterisms determined (ast), reciprocal cell
  !> symmetry known (recip), nearest-neighbor information (nn). If
  !> crash=.true., crash with error if the requested flag is not
  !> set. If crash=.false., take the necessary steps to initialize the
  !> crystal flags that are .true.  This routine is thread-safe if
  !> crash = .true.
  subroutine struct_checkflags(c,crash,env0,ast0,recip0,nn0,ewald0)
    use tools_io, only: ferror, faterr
    class(crystal), intent(inout) :: c
    logical :: crash
    logical, intent(in), optional :: env0
    logical, intent(in), optional :: ast0
    logical, intent(in), optional :: recip0
    logical, intent(in), optional :: nn0
    logical, intent(in), optional :: ewald0

    logical :: env, ast, recip, nn, ewald
    integer :: iast
    character(len=:), allocatable :: reason
    logical :: lflag(8)

    ! initialize optional arguments
    env = .false.
    if (present(env0)) env = env0
    ast = .false.
    if (present(ast0)) ast = ast0
    recip = .false.
    if (present(recip0)) recip = recip0
    nn = .false.
    if (present(nn0)) nn = nn0
    ewald = .false.
    if (present(ewald0)) ewald = ewald0

    lflag = .false.
    if (env .and. .not.c%isenv) lflag(2) = .true.
    if (ast .and. .not.c%isast) lflag(4) = .true.
    if (recip .and. .not.c%isrecip) lflag(5) = .true.
    if (nn .and. .not.c%isnn) lflag(6) = .true.
    if (ewald .and. .not.c%isewald) lflag(7) = .true.

    if (any(lflag)) then
       if (crash) then
          if (lflag(2)) reason = "atomic environments not determined for this crystal"
          if (lflag(4)) reason = "molecular connectivity has not been calculated"
          if (lflag(5)) reason = "reciprocal cell metrics and symmetry not determined"
          if (lflag(6)) reason = "nearest-neighbor information not available"
          if (lflag(7)) reason = "ewald cutoffs not available"
          call ferror('struct_checkflags',reason,faterr)
       else
          if (lflag(4)) then
             iast = 1
          else
             iast = 0
          end if
          call c%struct_fill(lflag(2),iast,lflag(5),lflag(6),lflag(7))
       end if
    end if

  end subroutine struct_checkflags

  !> Terminate allocated arrays
  subroutine struct_end(c)
    class(crystal), intent(inout) :: c

    c%isinit = .false.
    if (allocated(c%at)) deallocate(c%at)
    if (allocated(c%atcel)) deallocate(c%atcel)
    if (allocated(c%cen)) deallocate(c%cen)
    if (allocated(c%nside_ws)) deallocate(c%nside_ws)
    if (allocated(c%iside_ws)) deallocate(c%iside_ws)
    if (allocated(c%vws)) deallocate(c%vws)
    if (allocated(c%atenv)) deallocate(c%atenv)
    if (allocated(c%nstar)) deallocate(c%nstar)
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%moldiscrete)) deallocate(c%moldiscrete)
    c%isinit = .false.
    c%isenv = .false. 
    c%havesym = 0
    c%isast = .false. 
    c%isewald = .false. 
    c%isrecip = .false. 
    c%isnn = .false. 
    c%file = ""
    c%nneq = 0
    c%ncel = 0
    c%neqv = 0
    c%neqvg = 0
    c%ncv = 0
    c%ismolecule = .false.
    c%molx0 = 0d0
    c%molborder = 0d0
    c%nws = 0
    c%nvert_ws = 0
    c%nenv = 0
    c%nmol = 0

  end subroutine struct_end

  !> Transform crystallographic to cartesian. This routine is thread-safe.
  pure function x2c(c,xx) 
    class(crystal), intent(in) :: c
    real*8, intent(in) :: xx(3) 
    real*8 :: x2c(3)

    x2c = matmul(c%crys2car,xx)

  end function x2c

  !> Transform cartesian to crystallographic. This routine is thread-safe. 
  pure function c2x(c,xx)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: c2x(3)

    c2x = matmul(c%car2crys,xx)

  end function c2x

  !> Compute the distance between points in crystallographic.  This
  !> routine is thread-safe.
  pure function distance(c,x1,x2)
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in), dimension(3) :: x1 !< First point in cryst. coordinates
    real*8, intent(in), dimension(3) :: x2 !< Second point in cryst. coordinates
    real*8 :: distance

    real*8 :: xd(3)

    xd = c%x2c(x1 - x2)
    distance = sqrt(dot_product(xd,xd))

  end function distance

  !> Compute the shortest distance between a point x1 and all
  !> lattice translations of another point x2. Input points in cryst.
  !> coordinates. This routine is thread-safe.
  pure function eql_distance(c,x1,x2)
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in), dimension(3) :: x1 !< First point in cryst. coordinates
    real*8, intent(in), dimension(3) :: x2 !< Second point in cryst. coordinates
    real*8 :: eql_distance

    real*8 :: xd(3), dist2

    xd = x1 - x2
    call c%shortest(xd,dist2)
    eql_distance = sqrt(dist2)

  end function eql_distance

  !> Given a point in crystallographic coordinates (x), find the
  !> lattice-translated copy of x with the shortest length. Returns
  !> the shortest-length vector in cartesian coordinates and 
  !> the square of the distance. This routine is thread-safe.
  pure subroutine shortest(c,x,dist2)
    class(crystal), intent(in) :: c
    real*8, intent(inout) :: x(3)
    real*8, intent(out) :: dist2

    integer :: i
    real*8 :: x0(3), rvws(3), dvws

    x = x - nint(x)
    if (c%isortho) then
       x = matmul(c%crys2car,x)
       dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    else
       x0 = x
       x = matmul(c%crys2car,x)
       dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
       do i = 1, c%nws
          rvws = c%ivws(:,i) + x0
          rvws = matmul(c%crys2car,rvws)
          dvws = rvws(1)*rvws(1)+rvws(2)*rvws(2)+rvws(3)*rvws(3)
          if (dvws < dist2) then
             x = rvws
             dist2 = dvws
          endif
       end do
    endif

  end subroutine shortest

  !> Determine if two points x0 and x1 (cryst.) are at a distance less
  !> than eps. Logical veresion of c%distance(). If d2 is present and
  !> are_close is .true., return the square of the distance in that
  !> argument.  This routine is thread-safe.
  function are_close(c,x0,x1,eps,d2)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3), x1(3)
    real*8, intent(in) :: eps
    real*8, intent(out), optional :: d2
    logical :: are_close

    real*8 :: x(3), dbound, dist2

    are_close = .false.
    x = x0 - x1
    dbound = minval(abs(x)) * c%n2_c2x
    if (dbound > eps) return
    x = matmul(c%crys2car,x)
    if (any(abs(x) > eps)) return
    dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    are_close = (dist2 < (eps*eps))
    if (present(d2) .and. are_close) d2 = dist2

  end function are_close

  !> Determine if a points x0 is at a distance less than eps from x1
  !> or any of its lattice translations. x0 and x1 are in cryst.
  !> coords. Logical version of c%ldistance(). If d2 is present and
  !> are_close is .true., return the square of the distance in that
  !> argument. This routine is thread-safe.
  function are_lclose(c,x0,x1,eps,d2)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3), x1(3)
    real*8, intent(in) :: eps
    real*8, intent(out), optional :: d2
    logical :: are_lclose

    real*8 :: x(3), xa(3), dist2, dbound
    integer :: i

    are_lclose = .false.
    x = x0 - x1
    x = x - nint(x)
    dbound = minval(abs(x)) * c%n2_c2x
    if (c%isortho .or. c%ismolecule) then
       if (dbound > eps) return
       x = matmul(c%crys2car,x)
       if (any(abs(x) > eps)) return
       dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
       are_lclose = (dist2 < (eps*eps))
       goto 999
    else
       xa = x
       if (dbound < eps) then
          x = matmul(c%crys2car,x)
          if (.not.any(abs(x) > eps)) then
             dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
             are_lclose = (dist2 < (eps*eps))
             if (are_lclose) goto 999
          end if
       end if
       do i = 1, c%nws
          x = c%ivws(:,i) + xa
          dbound = minval(abs(x)) * c%n2_c2x
          if (dbound < eps) then
             x = matmul(c%crys2car,x)
             if (.not.any(abs(x) > eps)) then
                dist2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
                are_lclose = (dist2 < (eps*eps))
                if (are_lclose) goto 999
             end if
          end if
       end do
       are_lclose = .false.
    end if
    return

999 continue
    if (present(d2) .and. are_lclose) d2 = dist2

  end function are_lclose

  !> Given the point xp in crystallographic coordinates, calculates
  !> the nearest atom. If nid /= 0, then consider only atoms of the
  !> nid type (nneq atom list). In the output, nid represents the
  !> complete list id (atcel). dist is the distance and lvec the
  !> lattice vector required to transform atcel(nid)%x to the nearest
  !> position. This routine is thread-safe.
  subroutine nearest_atom(c,xp,nid,dist,lvec)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: xp(:)
    integer, intent(inout) :: nid
    real*8, intent(out) :: dist
    integer, intent(out) :: lvec(3)

    real*8 :: temp(3), d2, d2min
    integer :: j, nin

    nin = nid
    d2min = 1d30
    do j= 1, c%ncel
       if (nin /= 0 .and. c%atcel(j)%idx /= nin) cycle
       temp = c%atcel(j)%x - xp
       call c%shortest(temp,d2)
       if (d2 < d2min) then
          nid = j
          d2min = d2
          lvec = nint(c%atcel(j)%x - xp - temp)
       end if
    end do
    dist = sqrt(d2min)

  end subroutine nearest_atom

  !> Identify an atom in the unit cell. Input: cartesian coords. Output:
  !> the non-equivalent atom index (default if lncel is false) or the
  !> complete atom index (if lncel is true). This routine is
  !> thread-safe.
  function identify_atom(c,x0,lncel0)
    use tools_io, only: ferror, faterr
    
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    logical, intent(in), optional :: lncel0
    integer :: identify_atom

    real*8 :: x(3), xd(3), dist2
    integer :: i
    logical :: lncel

    lncel = .false.
    if (present(lncel0)) lncel = lncel0

    identify_atom = 0
    x = c%c2x(x0)
    do i = 1, c%ncel
       xd = x - c%atcel(i)%x
       call c%shortest(xd,dist2)
       if (dist2 < 1d-6) then
          if (lncel) then
             identify_atom = i
          else
             identify_atom = c%atcel(i)%idx
          endif
          return
       end if
    end do
    call ferror('identify_atom','atom in fragment not in list',faterr)

  endfunction identify_atom

  !> Identify a fragment in the unit cell. Input: cartesian coords. Output:
  !> A fragment object. This routine is thread-safe.
  function identify_fragment(c,nat,x0) result(fr)
    use types, only: realloc
    class(crystal), intent(in) :: c
    integer, intent(in) :: nat
    real*8, intent(in) :: x0(3,nat)
    type(fragment) :: fr

    integer :: id, i

    fr%nat = nat
    allocate(fr%at(nat))
    do i = 1, nat
       id = identify_atom(c,x0(:,i),.true.)
       fr%at(i)%r = x0(:,i)
       fr%at(i)%x = c%c2x(x0(:,i))
       fr%at(i)%cidx = id
       fr%at(i)%idx = c%atcel(id)%idx
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(id)%x)
       fr%at(i)%z = c%at(fr%at(i)%idx)%z
    end do
    call realloc(fr%at,fr%nat)

  end function identify_fragment

  !> Identify a fragment in the unit cell from an external
  !> xyz file. An instance of a fragment object is returned.
  function identify_fragment_from_xyz(c,file) result(fr)
    use tools_io, only: fopen_read, string, ferror, faterr, fclose
    use param, only: bohrtoa
    use types, only: realloc

    class(crystal), intent(in) :: c
    character*(*) :: file
    type(fragment) :: fr

    integer :: lu, nat
    integer :: id, i
    real*8 :: x0(3)
    character(len=:), allocatable :: word

    lu = fopen_read(file)
    read(lu,*,err=999) nat
    read(lu,*,err=999) 
    
    fr%nat = nat
    allocate(fr%at(nat))
    do i = 1, nat
       word = ""
       read(lu,*,err=999) word, x0
       x0 = x0 / bohrtoa - c%molx0
       id = c%identify_atom(x0,.true.)
       fr%at(i)%r = x0 
       fr%at(i)%x = c%c2x(x0)
       fr%at(i)%cidx = id
       fr%at(i)%idx = c%atcel(id)%idx
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(id)%x)
       fr%at(i)%z = c%at(fr%at(i)%idx)%z
    end do
    call fclose(lu)
    call realloc(fr%at,fr%nat)

    return
999 continue
    call ferror('identify_fragment_from_xyz','error reading xyz file: '//string(file),faterr)

  end function identify_fragment_from_xyz

  !> Obtain symmetry equivalent positions of xp0 and write them to
  !> vec, and the multiplicity of the xp0 position in mmult. xp0 and
  !> vec are in crystallographic coordinates. irotm and icenv contain
  !> the index of the rotation matrix and centering vectors
  !> responsible for the transformation of xp into the corresponding
  !> vec. eps is the minimum distance to consider two points
  !> equivalent (in bohr). vec, irotm, icenv, and eps0 are optional. 
  subroutine symeqv(c,xp0,mmult,vec,irotm,icenv,eps0)
    use types, only: realloc
    class(crystal), intent(in) :: c !< Input crystal
    real*8, dimension(3), intent(in) :: xp0 !< input position (crystallographic)
    integer, intent(out) :: mmult !< multiplicity
    real*8, allocatable, intent(inout), optional :: vec(:,:) !< sym-eq positions (crystallographic)
    integer, allocatable, intent(inout), optional :: irotm(:) !< index of the operation
    integer, allocatable, intent(inout), optional :: icenv(:) !< index of the cent. vector
    real*8, intent(in), optional :: eps0 !< Minimum distance to consider two vectors different (bohr)

    real*8 :: avec(3,c%neqv*c%ncv)
    integer :: i, j, k, l
    integer :: mrot, mrot0
    real*8 :: tmp(3), xp(3)
    real*8 :: l2
    real*8 :: loweps, dist2, eps

    real*8, parameter :: eps_default = 1d-2

    ! the eps
    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if

    !.Translate position to main cell
    xp = xp0 - floor(xp0)

    !.Run over symmetry operations and create the list of copies
    mrot0 = 0
    do i = 1, c%neqv
       do  j = 1, c%ncv
          mrot0 = mrot0 + 1
          avec(:,mrot0) = matmul(c%rotm(1:3,1:3,i),xp) + c%rotm(:,4,i) + c%cen(:,j)
          avec(:,mrot0) = avec(:,mrot0) - floor(avec(:,mrot0))
       enddo
    enddo

    if (present(vec)) then
       if (.not.allocated(vec)) then
          allocate(vec(3,mrot0))
       elseif (size(vec,2) < mrot0) then
          call realloc(vec,3,mrot0)
       endif
    end if

    ! calculate distances and (possibly) write vec
    mrot = 0
    d: do i = 1, mrot0
       do j = 1,i-1
          if (c%are_lclose(avec(:,i),avec(:,j),eps)) cycle d
       end do
       mrot = mrot + 1
       if (present(vec)) vec(:,mrot) = avec(:,i)
    end do d

    mmult=mrot
    if(.not.present(vec)) return
    call realloc(vec,3,mmult)

    if (.not.present(irotm).or..not.present(icenv)) return
    if (.not.allocated(irotm)) then
       allocate(irotm(mmult))
    else
       call realloc(irotm,mmult)
    endif
    if (.not.allocated(icenv)) then
       allocate(icenv(mmult))
    else
       call realloc(icenv,mmult)
    endif

    ! rotation matrix identifier
    loweps = 1d-2 * eps
    l2 = loweps*loweps
    alo: do j = 1, mmult
       blo: do k = 1, c%neqv
          clo: do l = 1, c%ncv
             ! generate equivalent position in zeroth cell
             ! rotm * (xp + L) + cv + L' == vec
             ! and check against known equivalent positions
             tmp = matmul(c%rotm(:,1:3,k),xp) + c%rotm(:,4,k) + c%cen(:,l)
             tmp = tmp - vec(:,j)
             call c%shortest(tmp,dist2)
             if (dist2 < l2) then
                irotm(j) = k
                icenv(j) = l
                exit blo
             end if
          end do clo
       end do blo
    end do alo

  end subroutine symeqv

  !> Calculate the multiplicity of the point x0 (cryst. coord.)
  function get_mult(c,x0) result (mult)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer :: mult

    call c%symeqv(x0,mult)

  end function get_mult

  !> Calculate the multiplicity of the point x0 in reciprocal space
  !> (fractional coordinates). 
  function get_mult_reciprocal(c,x0) result (mult)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer :: mult

    real*8 :: xp(3), x(3), d(3)
    real*8 :: avec(3,48)
    integer :: i, j, nvec
    logical :: found

    mult = 0

    !.Translate position to main cell
    xp = x0 - floor(x0)

    !.Run over symmetry operations and create the list of copies
    nvec = 0
    do i = 1, c%neqvg
       x = matmul(c%rotg(:,:,i),xp)
       found = .false.
       do j = 1, nvec
          d = x - avec(:,j)
          d = abs(d - nint(d))
          if (all(d < 1d-5)) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          nvec = nvec + 1
          avec(:,nvec) = x
       end if
    enddo
    mult = nvec

  end function get_mult_reciprocal

  !> Build succesive shells around the target point. Each shell is
  !> formed by all the identical atoms equidistant to the target.
  !> A density cutoff of 1d-12 is used to determine atoms that contribute
  !> to the unit cell's density. Used in the structure initialization.
  !> If dmax is given, use that number as an estimate of how many cells
  !> should be included in the search for atoms. 
  subroutine build_env(c,dmax0)
    use tools_math, only: norm
    use global, only: cutrad
    use types, only: realloc
    class(crystal), intent(inout) :: c !< Input crystal
    real*8, intent(in), optional :: dmax0

    integer :: i, j, k, l(3), m
    real*8 :: xx(3), dmax, sphmax, dist
    integer :: imax, jmax, kmax

    ! allocate atenv
    if (.not.allocated(c%atenv)) allocate(c%atenv(menv0))

    ! In molecules, use only the atoms in the main cell
    if (c%ismolecule) then
       c%nenv = c%ncel
       if (c%nenv > size(c%atenv)) &
          call realloc(c%atenv,c%nenv)
       l = 0
       do m = 1, c%ncel
          xx = c%atcel(m)%x
          c%atenv(m)%x = xx
          c%atenv(m)%r = c%x2c(xx)
          c%atenv(m)%idx = c%atcel(m)%idx
          c%atenv(m)%cidx = m
          c%atenv(m)%ir = c%atcel(m)%ir
          c%atenv(m)%ic = c%atcel(m)%ic
          c%atenv(m)%lvec = c%atcel(m)%lvec + l
          c%atenv(m)%lenv = l
       end do
       return
    endif

    sphmax = norm(c%x2c((/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
    sphmax = max(sphmax,norm(c%x2c((/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm(c%x2c((/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm(c%x2c((/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))

    if (present(dmax0)) then
       dmax = dmax0
    else
       dmax = 0d0
       do i = 1, c%nneq
          if (c%at(i)%z > 0) dmax = max(dmax,cutrad(c%at(i)%z))
       end do
    end if
    c%dmax0_env = dmax
    call search_lattice(c%crys2car,dmax,imax,jmax,kmax)

    ! build environment
    c%nenv = 0
    do i = -imax, imax
       do j = -jmax, jmax
          do k = -kmax, kmax
             !.run over the ions in the (i,j,k) cell
             do m = 1, c%ncel
                l = (/i,j,k/)
                xx = c%atcel(m)%x + l
                dist = norm(c%x2c(xx - (/0.5d0,0.5d0,0.5d0/)))
                if (dist > sphmax+dmax) cycle

                c%nenv = c%nenv + 1
                if (c%nenv > size(c%atenv)) then
                   call realloc(c%atenv,2*size(c%atenv))
                endif

                ! Store the point
                c%atenv(c%nenv)%x = xx
                c%atenv(c%nenv)%r = c%x2c(xx)
                c%atenv(c%nenv)%idx = c%atcel(m)%idx
                c%atenv(c%nenv)%cidx = m
                c%atenv(c%nenv)%ir = c%atcel(m)%ir
                c%atenv(c%nenv)%ic = c%atcel(m)%ic
                c%atenv(c%nenv)%lvec = c%atcel(m)%lvec + l
                c%atenv(c%nenv)%lenv = l
             enddo  !m
          enddo  !k
       enddo  !j
    enddo  !i

    call realloc(c%atenv,c%nenv)

  end subroutine build_env

  !> Find asterisms. For every atom in the unit cell, find the atoms in the 
  !> main cell or adjacent cells that are connected to it. 
  subroutine find_asterisms(c)
    use global, only: bondfactor
    use param, only: atmcov, vsmall

    class(crystal), intent(inout) :: c

    integer :: i, j, k
    real*8 :: rvws(3), x0(3), dist, dist2
    real*8 :: d0
    integer :: lvec0(3), lvec(3)

    if (allocated(c%nstar)) deallocate(c%nstar)
    if (.not.allocated(c%nstar)) allocate(c%nstar(c%ncel))

    ! allocate the neighbor star
    do i = 1, c%ncel
       allocate(c%nstar(i)%idcon(4))
       allocate(c%nstar(i)%lcon(3,4))
    end do

    if (c%ismolecule) then
       ! run over all pairs of atoms in the molecule
       lvec = 0
       do i = 1, c%ncel
          do j = i+1, c%ncel
             d0 = bondfactor * (atmcov(c%at(c%atcel(i)%idx)%z)+atmcov(c%at(c%atcel(j)%idx)%z))
             ! use the Cartesian directly
             x0 = c%atcel(j)%r - c%atcel(i)%r
             if (any(abs(x0) > d0)) cycle
             d0 = d0 * d0
             dist2 = x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3)
             if (dist2 < d0) then
                call addpair(i,j,lvec)
                call addpair(j,i,lvec)
             end if
          end do
       end do
    else
       ! run over all pairs of atoms in the unit cell
       do i = 1, c%ncel
          do j = i, c%ncel
             x0 = c%atcel(j)%x - c%atcel(i)%x
             lvec0 = nint(x0)
             x0 = x0 - lvec0
             d0 = bondfactor * (atmcov(c%at(c%atcel(i)%idx)%z)+atmcov(c%at(c%atcel(j)%idx)%z))

             do k = 0, c%nws
                if (k == 0) then
                   rvws = x0
                   lvec = lvec0
                else
                   rvws = x0 - c%ivws(:,k)
                   lvec = lvec0 + c%ivws(:,k)
                endif
                rvws = matmul(c%crys2car,rvws)
                if (all(abs(rvws) < d0)) then
                   dist = sqrt(rvws(1)*rvws(1)+rvws(2)*rvws(2)+rvws(3)*rvws(3))
                   if (dist > vsmall .and. dist < d0) then
                      call addpair(i,j,lvec)
                      call addpair(j,i,-lvec)
                   end if
                end if
             end do
          end do
       end do
    end if

  contains
    subroutine addpair(i,j,lvec)
      use types, only: realloc
      integer :: i, j, lvec(3)

      c%nstar(i)%ncon = c%nstar(i)%ncon + 1
      if (c%nstar(i)%ncon > size(c%nstar(i)%idcon)) then
         call realloc(c%nstar(i)%idcon,2*c%nstar(i)%ncon)
         call realloc(c%nstar(i)%lcon,3,2*c%nstar(i)%ncon)
      end if
      c%nstar(i)%idcon(c%nstar(i)%ncon) = j
      c%nstar(i)%lcon(:,c%nstar(i)%ncon) = -lvec

    end subroutine addpair
  end subroutine find_asterisms

  !> List atoms in a number of cells around the main cell (nx cells),
  !> possibly with border (doborder).
  function listatoms_cells(c,nx,doborder) result(fr)
    use types, only: realloc
    class(crystal), intent(in) :: c
    integer, intent(in) :: nx(3)
    logical, intent(in) :: doborder
    type(fragment) :: fr

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr

    integer :: ix, iy, iz, i
    logical :: if1

    allocate(fr%at(1))
    fr%nat = 0

    ! All atoms in these cells
    fr%nat = 0
    do ix = 0,nx(1)-1
       do iy = 0,nx(2)-1
          do iz = 0,nx(3)-1
             do i = 1, c%ncel
                fr%nat = fr%nat + 1
                if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                fr%at(fr%nat)%cidx = i
                fr%at(fr%nat)%idx = c%atcel(i)%idx
                fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                fr%at(fr%nat)%z = c%at(c%atcel(i)%idx)%z
             end do
          end do
       end do
    end do

    ! border: pick atoms
    if (doborder) then
       do ix = -1,nx(1)
          do iy = -1,nx(2)
             do iz = -1,nx(3)
                if (ix > -1 .and. ix < nx(1) .and. iy > -1 .and. iy < nx(2) .and.&
                   iz > -1 .and. iz < nx(3)) cycle
                do i = 1, c%ncel
                   ! border
                   if1 = (ix == -1 .and. c%atcel(i)%x(1)<rthr1 .or. &
                      ix == nx(1) .and. c%atcel(i)%x(1)>rthr .or. &
                      iy == -1 .and. c%atcel(i)%x(2)<rthr1 .or. &
                      iy == nx(2) .and. c%atcel(i)%x(2)>rthr .or. &
                      iz == -1 .and. c%atcel(i)%x(3)<rthr1 .or. &
                      iz == nx(3) .and. c%atcel(i)%x(3)>rthr)
                   if (.not.if1) then
                      fr%nat = fr%nat + 1
                      if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                      fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                      fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                      fr%at(fr%nat)%cidx = i
                      fr%at(fr%nat)%idx = c%atcel(i)%idx
                      fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                      fr%at(fr%nat)%z = c%at(c%atcel(i)%idx)%z
                      cycle
                   end if
                end do
             end do
          end do
       end do
    end if
    call realloc(fr%at,fr%nat)

  end function listatoms_cells

  !> List atoms inside a sphere of radius rsph and center xsph
  !> (cryst.)  or a cube of side rcub and center xcub (cryst.). Return
  !> the list of atomic positions (Cartesian) in x, the atomic numbers
  !> in z and the number of atoms in nat. 
  function listatoms_sphcub(c,rsph,xsph,rcub,xcub) result(fr)
    use tools_io, only: ferror, faterr
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in), optional :: rsph, xsph(3)
    real*8, intent(in), optional :: rcub, xcub(3)
    type(fragment) :: fr

    integer :: ix, iy, iz, i, nn
    real*8 :: x0(3), d, rsph2
    logical :: doagain, dosph

    if (.not.(present(rsph).and.present(xsph)).and..not.(present(rcub).and.present(xcub))) &
       call ferror("listatoms_sphcub","Need sphere or cube input",faterr)
    dosph = present(rsph)

    allocate(fr%at(1))
    fr%nat = 0

    ! all atoms in a sphere
    doagain = .true.
    nn = -1
    if (dosph) rsph2 = rsph * rsph
    do while(doagain)
       doagain = .false.
       nn = nn + 1
       do ix = -nn, nn
          do iy = -nn, nn
             do iz = -nn, nn
                if (abs(ix)/=nn .and. abs(iy)/=nn .and. abs(iz)/=nn) cycle
                do i = 1, c%ncel
                   if (dosph) then
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xsph)
                      if (all(abs(x0) > rsph)) cycle
                      d = dot_product(x0,x0)
                      if (d >= rsph2) cycle
                   else
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xcub)
                      if (any(abs(x0) > rcub)) cycle
                   endif

                   ! add this atom
                   fr%nat = fr%nat + 1
                   if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                   fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                   fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                   fr%at(fr%nat)%cidx = i
                   fr%at(fr%nat)%idx = c%atcel(i)%idx
                   fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                   fr%at(fr%nat)%z = c%at(c%atcel(i)%idx)%z
                   doagain = .true.
                end do
             end do
          end do
       end do
    end do
    call realloc(fr%at,fr%nat)

  end function listatoms_sphcub

  !> Using the calculated asterisms for each atom determine the
  !> molecular in the system and whether the crystal is extended or
  !> molecular. This routine fills nmol, mol, and moldiscrete.
  subroutine fill_molecular_fragments(c)
    use fragmentmod, only: fragment_cmass
    use tools_io, only: ferror, faterr
    use types, only: realloc
    class(crystal), intent(inout) :: c

    integer :: i, j, k, l, jid, newid, newl(3)
    integer :: nat
    logical :: used(c%ncel), found, fdisc
    integer, allocatable :: id(:), lvec(:,:)
    logical, allocatable :: ldone(:)
    real*8 :: xcm(3)

    if (.not.allocated(c%nstar)) &
       call ferror('fill_molecular_fragments','no asterisms found',faterr)
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%moldiscrete)) deallocate(c%moldiscrete)

    ! initizialize 
    used = .false.
    c%nmol = 0
    allocate(c%mol(1),c%moldiscrete(1),id(10),lvec(3,10),ldone(10))
    c%moldiscrete = .true.

    ! run over atoms in the unit cell
    do i = 1, c%ncel
       if (used(i)) cycle
       
       ! increment the fragment counter
       c%nmol = c%nmol + 1
       if (c%nmol > size(c%mol)) then
          call realloc(c%mol,2*c%nmol)
          call realloc(c%moldiscrete,2*c%nmol)
          c%moldiscrete(c%nmol:2*c%nmol) = .true.
       end if

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = i
       lvec(:,1) = 0d0
       ldone(1) = .false.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)

             ! Is this atom in the fragment already? -> skip it. If it
             ! has a different lattice vector, mark the fragment as
             ! not discrete.
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                fdisc = all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) then
                if (.not.fdisc) then
                   c%moldiscrete(c%nmol) = .false.
                end if
                cycle
             end if

             ! Have we used this atom already?
             if (used(newid)) cycle

             ! Add the atom to the stack and mark it as used.
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl
             used(newid) = .true.
             ldone(nat) = .false.
          end do
       end do
       
       ! add this fragment to the list
       used(i) = .true.
       allocate(c%mol(c%nmol)%at(nat))
       c%mol(c%nmol)%nat = nat
       do j = 1, nat
          c%mol(c%nmol)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          c%mol(c%nmol)%at(j)%r = c%x2c(c%mol(c%nmol)%at(j)%x)
          c%mol(c%nmol)%at(j)%cidx = id(j)
          c%mol(c%nmol)%at(j)%idx = c%atcel(id(j))%idx
          c%mol(c%nmol)%at(j)%lvec = lvec(:,j)
          c%mol(c%nmol)%at(j)%z = c%at(c%mol(c%nmol)%at(j)%idx)%z
       end do
    end do
    call realloc(c%mol,c%nmol)
    call realloc(c%moldiscrete,c%nmol)

    ! translate all fragments to the main cell
    if (.not.c%ismolecule) then
       do i = 1, c%nmol
          xcm = fragment_cmass(c%mol(i))
          newl = floor(c%c2x(xcm))
          do j = 1, c%mol(i)%nat
             c%mol(i)%at(j)%x = c%mol(i)%at(j)%x - newl
             c%mol(i)%at(j)%r = c%x2c(c%mol(i)%at(j)%x)
             c%mol(i)%at(j)%lvec = c%mol(i)%at(j)%lvec - newl
          end do
       end do
    end if

  end subroutine fill_molecular_fragments

  !> List all molecules resulting from completing the initial fragment
  !> fri by adding adjacent atoms that are covalently bonded. Return
  !> the number of fragment (nfrag), the fragments themselves (fr),
  !> and whether the fragments are discrete (not connected to copies
  !> of themselves in a different cell). 
  subroutine listmolecules(c,fri,nfrag,fr,isdiscrete)
    use types, only: realloc
    class(crystal), intent(inout) :: c
    type(fragment), intent(in) :: fri
    integer, intent(out) :: nfrag
    type(fragment), intent(out), allocatable :: fr(:)
    logical, intent(out), allocatable :: isdiscrete(:)
    
    integer :: i, j, k, l, newid, newl(3), jid
    integer :: nat
    integer, allocatable :: id(:), lvec(:,:)
    logical, allocatable :: ldone(:)
    logical :: found, ldist
    integer :: nseed
    integer, allocatable :: idseed(:), lseed(:,:)
    logical, allocatable :: fseed(:)

    ! find the neighbor stars, if not already done
    call c%checkflags(.false.,ast0=.true.)

    ! unwrap the input fragment
    nseed = fri%nat
    allocate(idseed(nseed),lseed(3,nseed),fseed(nseed))
    do j = 1, nseed
       idseed(j) = fri%at(j)%cidx
       lseed(:,j) = fri%at(j)%lvec
    end do
    fseed = .false.

    ! allocate stuff
    nfrag = 0
    allocate(fr(1),isdiscrete(1),id(10),lvec(3,10),ldone(10))
    isdiscrete = .true.
    
    do i = 1, nseed
        if (fseed(i)) cycle

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = idseed(i)
       lvec(:,1) = lseed(:,i)
       ldone(1) = .false.
       ldist = .true.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)
 
             ! is this atom in the fragment already? -> skip it
             found = .false.
             do l = 1, nat
                found = (newid == id(l)) .and. all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) cycle

             ! is this atom in the fragment already with a different
             ! lattice vector?  -> add it to the list but not to the
             ! stack, and mark the fragment as non-discrete
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                if (found) exit
             end do
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl

             if (found) then
                ldone(nat) = .true.
                ldist = .false.
             else
                ! if it wasn't found, then add the atom to the stack
                ldone(nat) = .false.
             end if
          end do
       end do

       ! add this fragment to the list
       fseed(i) = .true.
       nfrag = nfrag + 1
       if (nfrag > size(fr)) then
          call realloc(fr,2*nfrag)
          call realloc(isdiscrete,2*nfrag)
       end if
       allocate(fr(nfrag)%at(nat))
       isdiscrete(nfrag) = ldist
       fr(nfrag)%nat = nat
       do j = 1, nat
          fr(nfrag)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          fr(nfrag)%at(j)%r = c%x2c(fr(nfrag)%at(j)%x)
          fr(nfrag)%at(j)%cidx = id(j)
          fr(nfrag)%at(j)%idx = c%atcel(id(j))%idx
          fr(nfrag)%at(j)%lvec = lvec(:,j) 
          fr(nfrag)%at(j)%z = c%at(fr(nfrag)%at(j)%idx)%z
       end do

       ! run over all atoms in the new fragment and mark those atoms in the seed
       do j = 1, nat
          do k = 1, nseed
             if (fseed(k)) cycle
             if (id(j) == idseed(k) .and. all(lvec(:,j) == lseed(:,k))) &
                fseed(k) = .true.
          end do
       end do
    end do
    call realloc(fr,nfrag)
    call realloc(isdiscrete,nfrag)

  end subroutine listmolecules

  !> Calculates the neighbor environment of a given point x0 (cryst.
  !> coords) up to shell shmax. Return the number of neighbors for
  !> each shell in nneig, the non-equivalent atom index in wat,
  !> and the distance in dist. If the argument xenv is present,
  !> return the position of a representative atom from each shell.
  subroutine pointshell(c,x0,shmax,nneig,wat,dist,xenv)
    use global, only: atomeps, atomeps2
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: shmax
    integer, intent(out) :: nneig(shmax)
    integer, intent(out) :: wat(shmax)
    real*8, intent(out) :: dist(shmax)
    real*8, intent(inout), allocatable, optional :: xenv(:,:,:)

    integer :: j, l, n
    real*8 :: aux(shmax), d2, x0c(3)
    integer :: iaux(shmax)
    real*8, allocatable :: aux2(:,:,:)

    x0c = c%x2c(x0)
    dist = 1d30
    if (present(xenv)) then
       if (allocated(xenv)) deallocate(xenv)
       allocate(xenv(3,5,shmax),aux2(3,5,shmax))
       xenv = 0d0
    endif
    nneig = 0
    wat = 0
    do j = 1, c%nenv
       d2 = dot_product(c%atenv(j)%r-x0c,c%atenv(j)%r-x0c)
       if (d2 < atomeps2) cycle
       do l = 1, shmax
          if (abs(d2 - dist(l)) < atomeps) then
             nneig(l) = nneig(l) + 1
             if (present(xenv)) then
                if (nneig(l) > size(xenv,2)) then
                   n = size(xenv)
                   call realloc(xenv,3,2*n,shmax)
                   call realloc(aux2,3,2*n,shmax)
                   xenv(:,n+1:,:) = 0d0
                endif
                xenv(:,nneig(l),l) = c%atenv(j)%x
             endif
             exit
          else if (d2 < dist(l)) then
             aux(l:shmax-1) = dist(l:shmax-1)
             dist(l+1:shmax) = aux(l:shmax-1)
             iaux(l:shmax-1) = nneig(l:shmax-1)
             nneig(l+1:shmax) = iaux(l:shmax-1)
             iaux(l:shmax-1) = wat(l:shmax-1)
             wat(l+1:shmax) = iaux(l:shmax-1)
             if (present(xenv)) then
                aux2(:,:,l:shmax-1) = xenv(:,:,l:shmax-1)
                xenv(:,:,l+1:shmax) = aux2(:,:,l:shmax-1)
             endif

             dist(l) = d2
             nneig(l) = 1
             wat(l) = c%atenv(j)%idx
             if (present(xenv)) then
                xenv(:,1,l) = c%atenv(j)%x
             endif
             exit
          end if
       end do
    end do
    dist = sqrt(dist)
    if (allocated(aux2)) deallocate(aux2)

  end subroutine pointshell

  !> Determines the site symmetry for a point x0 in cryst. coords.
  !> Two points are the same if their distance is less than eps0.
  !> Returns the site symmetry group symbol (sitesymm), the
  !> number of operations in this group (leqv) and the rotation
  !> operations (lrotm)
  function sitesymm(c,x0,eps0,leqv,lrotm)
    use tools_io, only: string
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x0(3) !< Input point in cryst. coords.
    real*8, intent(in), optional :: eps0 !< Two points are different if distance is > eps
    character*3 :: sitesymm !< point group symbol
    integer, optional :: leqv !< Number of operations in the group
    real*8, optional :: lrotm(3,3,48) !< Point group operations

    integer :: i, m
    real*8 :: dumy(3), eps2, dist2, eps, vec(3)
    integer :: type
    logical :: ok
    integer :: highest, highests
    integer :: nnsym, order
    integer :: ordersym(c%neqv), masksym(0:9)

    real*8, parameter :: eps_default = 1d-2

    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if
    eps2 = eps*eps

    ! Run over all proper symmetry elements of symmetry
    nnsym = 0
    masksym = 0
    ordersym = 0
    do i = 1, c%neqv
       ok = .false.
       do m = 1, c%ncv
          dumy = matmul(c%rotm(1:3,1:3,i),x0) 
          dumy = - dumy + x0 - c%rotm(1:3,4,i) - c%cen(:,m)
          call c%shortest(dumy,dist2)
          ok = (dist2 < eps2)
          if (ok) exit
       end do
       if (ok) then
          ! A symmetry operation at location x has been found.
          nnsym = nnsym+1
          call typeop(c%rotm(:,:,i),type,vec,order)
          ordersym(nnsym) = order
          masksym(type) = masksym(type) + 1
          if (present(leqv).and.present(lrotm)) then
             leqv = nnsym
             lrotm(:,:,leqv) = c%rotm(1:3,1:3,i)
          end if
       endif
    enddo

    ! calculate the point group
    sitesymm = ""
    if (masksym(c3) > 2) then
       ! cubic groups
       if (masksym(c4) /= 0) then
          if (masksym(inv) /= 0) then
             sitesymm='Oh'
          else
             sitesymm='O'
          endif
       elseif (masksym(s4).ne.0) then
          sitesymm= 'Td'
       elseif (masksym(inv).ne.0) then
          sitesymm = 'Th'
       else
          sitesymm = 'T'
       endif
    else
       !Compute highest order proper axis.
       highest=0
       highests=0
       if (masksym(c2) /= 0) highest=2
       if (masksym(c3) /= 0) highest=3
       if (masksym(c4) /= 0) highest=4
       if (masksym(c6) /= 0) highest=6
       if (masksym(s3) /= 0) highests=3
       if (masksym(s4) /= 0) highests=4
       if (masksym(s6) /= 0) highests=6
       if (highest == 0) then
          if (masksym(inv) /= 0) then
             sitesymm='i'
          elseif (masksym(sigma) /= 0) then
             sitesymm='Cs'
          else
             sitesymm='C1'
          endif
       elseif (masksym(c2) >= highest) then
          if (masksym(sigma) == 0) then
             sitesymm='D' // string(highest)
          elseif (masksym(inv) .eq. 1) then
             if (highest == 3) then
                sitesymm= 'D3d'
             else
                sitesymm= 'D' // string(highest) // 'h'
             endif
          else
             if (highest .eq. 3) then
                sitesymm= 'D3h'
             else
                sitesymm= 'D' // string(highest) // 'd'
             endif
          endif
       elseif (masksym(sigma) == 0) then
          if (highests /= 0) then
             sitesymm= 'S' // string(highests/2)
          else
             sitesymm= 'C' // string(highest)
          endif
       elseif (masksym(sigma) .lt. highest) then
          sitesymm= 'C' // string(highest) // 'h'
       else
          sitesymm= 'C' // string(highest) // 'v'
       endif
    endif

  end function sitesymm

  !> Calculate the packing ratio (in %) using the nearest-neighbor
  !> information. Each atom is assigned a ratio equal to half the distance
  !> to its nearest neighbor.
  function get_pack_ratio(c) result (px)
    use param, only: pi
    class(crystal), intent(inout) :: c
    real*8 :: px
    
    integer :: i

    px = 0d0
    do i = 1, c%nneq
       px = px + c%at(i)%mult * 4d0/3d0 * pi * c%at(i)%rnn2**3
    end do
    px = px / c%omega * 100d0

  end function get_pack_ratio

  !> Calculate the powder diffraction pattern. 
  !> On input, npts is the number of 2*theta points from the initial
  !> (th2ini0) to the final (th2end0) 2*theta values. Both angles are
  !> in degrees. lambda0 is the wavelength of the radiation (in
  !> angstrom). sigma is the parameter for the Gaussian broadening.
  !> fpol is the polarization correction factor (0 = unpolarized, 0.95
  !> = syncrhotron).
  !> On output, t is the 2*theta grid, ih is the intensity on the
  !> 2*theta grid, th2p is the 2*theta for the located maxima, ip is
  !> the list of maxima itensities, and hvecp is the reciprocal
  !> lattice vector corresponding to the peaks.
  subroutine powder(c,th2ini0,th2end0,npts,lambda0,fpol,&
     sigma,t,ih,th2p,ip,hvecp)
    use param, only: pi, bohrtoa, cscatt, c2scatt
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in) :: th2ini0, th2end0
    integer, intent(in) :: npts
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    real*8, intent(in) :: sigma
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)
    real*8, allocatable, intent(inout) :: th2p(:)
    real*8, allocatable, intent(inout) :: ip(:)
    integer, allocatable, intent(inout) :: hvecp(:,:)

    integer :: i, ii, np, hcell, h, k, l, iz, idx
    real*8 :: th2ini, th2end, lambda, hvec(3), kvec(3), th, sth, th2
    real*8 :: sigma2, smax, dh2, dh, dh3, sthlam, cterm, sterm
    real*8 :: ffac, as(4), bs(4), cs, c2s(4), int, mcorr, afac
    real*8 :: ipmax, ihmax
    integer :: hmax
    integer, allocatable :: multp(:)
    integer, allocatable :: io(:)
    real*8, allocatable :: th2p_(:), ip_(:)
    integer, allocatable :: hvecp_(:,:)

    integer, parameter :: mp = 20
    real*8, parameter :: ieps = 1d-5
    real*8, parameter :: theps = 1d-5

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts))
    do i = 1, npts
       t(i) = th2ini0 + real(i-1,8) / real(npts-1,8) * (th2end0-th2ini0)
    end do
    ih = 0d0
    th2ini = th2ini0 * pi / 180d0
    th2end = th2end0 * pi / 180d0

    ! allocate for peak list
    if (allocated(th2p)) deallocate(th2p)
    if (allocated(ip)) deallocate(ip)
    if (allocated(hvecp)) deallocate(hvecp)
    allocate(th2p(mp),ip(mp),multp(mp),hvecp(3,mp))

    ! cell limits, convert lambda to bohr
    lambda = lambda0 / bohrtoa
    smax = sin(th2end/2d0)
    hmax = 2*ceiling(2*smax/lambda/minval(c%ar))
    ! broadening -> gaussian
    sigma2 = sigma * sigma

    ! calculate the intensities
    np = 0
    do hcell = 1, hmax
       do h = -hcell, hcell
          do k = -hcell, hcell
             do l = -hcell, hcell
                if (abs(h)/=hcell.and.abs(k)/=hcell.and.abs(l)/=hcell) cycle
                ! reciprocal lattice vector length
                hvec = real((/h,k,l/),8)
                dh2 = dot_product(hvec,matmul(c%grtensor,hvec))
                dh = sqrt(dh2)
                dh3 = dh2 * dh

                ! the theta is not outside the spectrum range
                sth = 0.5d0 * lambda * dh
                if (abs(sth) > smax) cycle
                th = asin(sth)
                th2 = 2d0 * th
                if (th2 < th2ini .or. th2 > th2end) cycle

                ! more stuff we need
                sthlam = dh / bohrtoa / 2d0
                kvec = 2 * pi * hvec

                ! calculate the raw intensity for this (hkl)
                cterm = 0d0
                sterm = 0d0
                do i = 1, c%ncel
                   iz = c%at(c%atcel(i)%idx)%z
                   if (iz < 1 .or. iz > size(cscatt,2)) &
                      call ferror('struct_powder','invalid Z -> no atomic scattering factors',faterr)
                   as = (/cscatt(1,iz),cscatt(3,iz),cscatt(5,iz),cscatt(7,iz)/)
                   bs = (/cscatt(2,iz),cscatt(4,iz),cscatt(6,iz),cscatt(8,iz)/)
                   cs = cscatt(9,iz)
                   if (dh < 2d0) then
                      ffac = as(1)*exp(-bs(1)*dh2)+as(2)*exp(-bs(2)*dh2)+&
                         as(3)*exp(-bs(3)*dh2)+as(4)*exp(-bs(4)*dh2)+cs
                   elseif (iz == 1) then
                      ffac = 0d0
                   else
                      c2s = c2scatt(:,iz)
                      ffac = exp(c2s(1)+c2s(2)*dh+c2s(3)*dh2/10d0+c2s(4)*dh3/100d0)
                   end if
                   ffac = ffac * exp(-sthlam**2)
                   cterm = cterm + ffac * cos(dot_product(kvec,c%atcel(i)%x))
                   sterm = sterm + ffac * sin(dot_product(kvec,c%atcel(i)%x))
                end do
                int = cterm**2 + sterm**2

                ! profile correction
                ! Yinghua J. Appl. Cryst. 20 (1987) 258
                ! lpf = (1 + cos(th2)**2) / sin(th)**2
                ! int = int * lpf

                ! gdis lorentz-polarization correction; checked in diamond; 
                ! criticized by Yinghua because it is a correction for the integrated
                ! intensity
                ! lpf = 0.5*(1+cos(th2)**2) / sin(th2) / sin(th)
                ! int = int * lpf

                ! FoX-compatible
                ! lorentz correction
                mcorr = 1d0 / sin(th2)
                ! slit aperture
                mcorr = mcorr / sin(th)
                ! polarization
                afac = (1-fpol) / (1+fpol)
                mcorr = mcorr * (1+afac*(0.5d0+0.5d0*cos(2*th2))) / (1+afac)
                int = int * mcorr
                
                ! sum the peak
                if (int > ieps) then
                   ! use a gaussian profile, add the intensity
                   ih = ih + int * exp(-(t-th2*180/pi)**2 / 2d0 / sigma2)

                   ! identify the new peak
                   if (all(abs(th2p(1:np)-th2) > theps)) then
                      np = np + 1
                      if (np > size(th2p)) then
                         call realloc(th2p,2*np)
                         call realloc(ip,2*np)
                         call realloc(multp,2*np)
                         call realloc(hvecp,3,2*np)
                      end if
                      th2p(np) = th2
                      ip(np) = int
                      multp(np) = 1
                      hvecp(:,np) = (/h,k,l/)
                   else
                      do idx = 1, np
                         if (abs(th2p(idx)-th2) <= theps) exit
                      end do
                      multp(idx) = multp(idx) + 1
                      ! usually the hvec with the most positive indices is the last one
                      hvecp(:,idx) = (/h,k,l/)
                   endif
                end if
             end do
          end do
       end do
    end do
    call realloc(th2p,np)
    call realloc(ip,np)
    call realloc(multp,np)
    call realloc(hvecp,3,np)

    ! normalize the intensities to 100
    if (np == 0) &
       call ferror('struct_powder','no peaks found in the 2theta range',faterr)
    ip = ip * multp
    ipmax = maxval(ip)
    ihmax = maxval(ih)
    ih = ih / ihmax * 100
    ip = ip / ipmax * 100

    ! deallocate the multiplicities
    deallocate(multp)

    ! sort the peaks
    allocate(io(np),th2p_(np),ip_(np),hvecp_(3,np))
    do i = 1, np
       io(i) = i
    end do
    call qcksort(th2p,io,1,np)
    do ii = 1, np
       i = io(ii)
       th2p_(ii) = th2p(i)
       ip_(ii) = ip(i)
       hvecp_(:,ii) = hvecp(:,i)
    end do
    th2p = th2p_
    ip = ip_
    hvecp = hvecp_
    deallocate(th2p_,ip_,hvecp_,io)

  end subroutine powder

  !> Calculate the radial distribution function.  On input, npts is
  !> the number of bins points from the initial (0) to the final
  !> (rend) distance. On output, t is the distance grid, and ih is the
  !> value of the RDF. This routine is based on:
  !>   Willighagen et al., Acta Cryst. B 61 (2005) 29.
  !> except using the sqrt of the atomic numbers instead of the 
  !> charges.
  subroutine rdf(c,rend,npts,t,ih)
    use tools_math, only: norm
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rend
    integer, intent(in) :: npts
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)

    integer :: i, j, ibin
    real*8 :: d, hfac

    ! integer :: i, ii, np, hcell, h, k, l, iz, idx
    ! real*8 :: th2ini, th2end, lambda, hvec(3), kvec(3), th, sth, th2
    ! real*8 :: smax, dh2, dh, dh3, sthlam, cterm, sterm
    ! real*8 :: ffac, as(4), bs(4), cs, c2s(4), int, mcorr, afac
    ! real*8 :: ipmax, ihmax
    ! integer :: hmax
    ! integer, allocatable :: multp(:)
    ! integer, allocatable :: io(:)
    ! real*8, allocatable :: th2p_(:), ip_(:)
    ! integer, allocatable :: hvecp_(:,:)

    ! integer, parameter :: mp = 20
    ! real*8, parameter :: ieps = 1d-5
    ! real*8, parameter :: theps = 1d-5

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts))
    do i = 1, npts
       t(i) = real(i-1,8) / real(npts-1,8) * rend
    end do
    ih = 0d0

    ! calculate the radial distribution function for the crystal
    ! RDF(r) = sum_i=1...c%nneq sum_j=1...c%nenv sqrt(Zi*Zj) / c%nneq / rij * delta(r-rij)
    hfac = (npts-1) / rend
    do i = 1, c%nneq
       do j = 1, c%nenv
          d = norm(c%at(i)%r - c%atenv(j)%r)
          ibin = nint(d * hfac) + 1
          if (ibin <= 0 .or. ibin > npts) cycle
          ih(ibin) = ih(ibin) + sqrt(real(c%at(i)%z * c%at(c%atenv(j)%idx)%z,8))
       end do
    end do
    do i = 2, npts
       ih(i) = ih(i) / t(i)
    end do
    ih = ih / c%nneq

  end subroutine rdf

  !> Calculate real and reciprocal space sum cutoffs
  subroutine calculate_ewald_cutoffs(c)
    use tools_io, only: ferror, faterr
    use param, only: pi, rad, sqpi, tpi

    class(crystal), intent(inout) :: c

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux, q2sum
    integer :: ia, ib, ic
    real*8 :: alrmax(3)
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    ! calculate sum of charges and charges**2
    c%qsum = 0d0
    q2sum = 0d0
    do i = 1, c%nneq
       if (abs(c%at(i)%qat) < 1d-6) &
          call ferror('ewald_energy','Some of the charges are 0',faterr)
       c%qsum = c%qsum + real(c%at(i)%mult * c%at(i)%qat,8)
       q2sum = q2sum + real(c%at(i)%mult * c%at(i)%qat**2,8)
    end do

    ! determine shortest vector in real space
    aux = 0d0
    do i = 1, 3
       if (c%aa(i) > aux) then
          aux = c%aa(i)
          ia = i
       end if
    end do
    ! determine shortest vector in reciprocal space, dif. from ia
    aux = 0d0
    do i = 1, 3
       if (c%ar(i) > aux .and. i /= ia) then
          aux = c%ar(i)
          ic = i
       end if
    end do
    ! the remaining vector is ib
    ib = 1
    do i = 1, 3
       if (i /= ia .and. i /= ic) ib = i
    end do

    ! convergence parameter
    c%eta = sqrt(c%omega / pi / c%aa(ib) / sin(c%bb(ic)*rad))

    ! real space cutoff
    rcut1 = 1d0
    rcut2 = 2d0 / sgrow
    err_real = 1d30
    do while (err_real >= eeps)
       rcut2 = rcut2 * sgrow
       err_real = pi * c%ncel**2 * q2sum / c%omega * c%eta**2 * erfc(rcut2 / c%eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       c%rcut = 0.5*(rcut1+rcut2)
       err_real = pi * c%ncel**2 * q2sum / c%omega * c%eta**2 * erfc(c%rcut / c%eta)
       if (err_real > eeps) then
          rcut1 = c%rcut
       else
          rcut2 = c%rcut
       endif
    end do
    c%rcut = 0.5*(rcut1+rcut2)
    ! real space cells to explore
    alrmax = 0d0
    alrmax(1) = c%aa(2) * c%aa(3) * sin(c%bb(1)*rad)
    alrmax(2) = c%aa(1) * c%aa(3) * sin(c%bb(2)*rad)
    alrmax(3) = c%aa(1) * c%aa(2) * sin(c%bb(3)*rad)
    c%lrmax = ceiling(c%rcut * alrmax / c%omega)

    ! reciprocal space cutoff
    hcut1 = 1d0
    hcut2 = 2d0 / sgrow
    err_rec = 1d30
    do while(err_rec >= eeps)
       hcut2 = hcut2 * sgrow
       err_rec = c%ncel**2 * q2sum / sqpi / c%eta * erfc(c%eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       c%hcut = 0.5*(hcut1+hcut2)
       err_rec = c%ncel**2 * q2sum / sqpi / c%eta * erfc(c%eta * c%hcut / 2)
       if (err_rec > eeps) then
          hcut1 = c%hcut
       else
          hcut2 = c%hcut
       endif
    end do
    c%hcut = 0.5*(hcut1+hcut2)
    ! reciprocal space cells to explore
    c%lhmax = ceiling(c%aa(ia) / tpi * c%hcut)

  end subroutine calculate_ewald_cutoffs

  !> Given a crystal structure (c) and three lattice vectors in cryst.
  !> coords (x0(:,1), x0(:,2), x0(:,3)), build the same crystal
  !> structure using the unit cell given by those vectors. 
  subroutine newcell(c,x00,t0,verbose0)
    use tools_math, only: det, matinv, mnorm2
    use tools_io, only: ferror, faterr, warning, string, uout
    use param, only: pi, ctsq3
    use types, only: realloc
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x00(3,3)
    real*8, intent(in), optional :: t0(3)
    logical, intent(in), optional :: verbose0

    type(crystalseed) :: ncseed
    logical :: ok, found, verbose
    real*8 :: x0(3,3), x0inv(3,3), fvol
    real*8 :: r(3,3), g(3,3), x(3), dx(3), dd, t(3)
    integer :: i, j, k, l, m
    integer :: nr, nn
    integer :: nlat
    real*8, allocatable :: xlat(:,:)

    real*8, parameter :: eps = 1d-6

    if (c%ismolecule) &
       call ferror('newcell','NEWCELL incompatible with molecules',faterr)

    ! initialize
    x0 = x00
    dd = det(x0)

    if (abs(dd) < eps) then
       call ferror('newcell','invalid input vectors',faterr)
    elseif (dd < 0d0) then
       ! flip the cell
       x0 = -x0
       dd = det(x0)
       call ferror('newcell','det < 0; vectors flipped',warning)
    endif
    if (present(t0)) then
       t = t0
    else
       t = 0d0
    end if
    verbose = .false.
    if (present(verbose0)) verbose = verbose0

    ! check that the vectors are pure translations
    do i = 1, 3
       ok = .false.
       do j = 1, c%ncv
          ok = (c%are_lclose(x0(:,i),c%cen(:,j),1d-4))
          if (ok) exit
       end do
       if (.not.ok) &
          call ferror("newcell","Cell vector number " // string(i) // &
          " is not a pure translation",faterr)
    end do

    ! is this a smaller or a larger cell? Arrange vectors.
    if (abs(dd-1d0) < eps) then
       nr = 1
    elseif (dd > 1d0) then
       nr = nint(dd)
       if (abs(nr-dd) > eps) &
          call ferror('newcell','inconsistent determinant of lat. vectors',faterr)
    else
       nr = -nint(1d0/dd)
       if (abs(-nr-1d0/dd) > eps) &
          call ferror('newcell','inconsistent determinant of lat. vectors',faterr)
    end if

    if (verbose) then
       write (uout,'("* Transformation to a new unit cell (NEWCELL)")')
       write (uout,'("  Lattice vectors of the new cell in the old setting (cryst. coord.):")')
       write (uout,'(4X,3(A,X))') (string(x0(i,1),'f',12,7,4),i=1,3)
       write (uout,'(4X,3(A,X))') (string(x0(i,2),'f',12,7,4),i=1,3)
       write (uout,'(4X,3(A,X))') (string(x0(i,3),'f',12,7,4),i=1,3)
       write (uout,'("  Origin translation: ",3(A,X))') (string(t(i),'f',12,7,4),i=1,3)
       write (uout,*)
    end if

    ! inverse matrix
    x0inv = matinv(x0)

    ! metrics of the new cell
    r = matmul(transpose(x0),transpose(c%crys2car))
    ncseed%crys2car = transpose(r)
    ncseed%useabr = 2
    fvol = abs(det(r)) / c%omega
    if (abs(nint(fvol)-fvol) > eps .and. abs(nint(1d0/fvol)-1d0/fvol) > eps) &
       call ferror("newcell","Inconsistent newcell volume",faterr)

    ! find a star of lattice vectors and supercell centering vectors, if any
    ! first lattice vector is (0 0 0)
    allocate(xlat(3,10))
    xlat = 0d0
    nlat = 1
    do i = minval(floor(x0(1,:))),maxval(ceiling(x0(1,:)))
       do j = minval(floor(x0(2,:))),maxval(ceiling(x0(2,:)))
          do k = minval(floor(x0(3,:))),maxval(ceiling(x0(3,:)))
             x = matmul((/i, j, k/),transpose(x0inv))
             if (any(abs(x-nint(x)) > eps)) then
                ! this is a new candidate for supercell centering vector
                ! check if we have it already
                x = x - floor(x)
                found = .false.
                do l = 1, nlat
                   if (all(abs(xlat(:,l) - x) < eps)) then
                      found = .true.
                      exit
                   end if
                end do
                if (.not.found) then
                   nlat = nlat + 1
                   if (nlat > size(xlat,2)) call realloc(xlat,3,2*nlat)
                   xlat(:,nlat) = x
                end if
             endif
          end do
       end do
    end do

    ! build the new atom list
    ncseed%nat = 0
    nn = ceiling(c%ncel * abs(det(r)) / c%omega)
    allocate(ncseed%x(3,nn),ncseed%z(nn),ncseed%name(nn))
    do i = 1, nlat
       do j = 1, c%ncel
          ! candidate atom
          x = matmul(c%atcel(j)%x-t,transpose(x0inv)) + xlat(:,i)
          x = x - floor(x)

          ! check if we have it already
          ok = .true.
          do m = 1, ncseed%nat
             dx = x - ncseed%x(:,m)
             dx = abs(dx - nint(dx))
             if (all(dx < eps)) then
                ok = .false.
                exit
             end if
          end do
          if (ok) then
             ! add it to the list
             ncseed%nat = ncseed%nat + 1
             if (ncseed%nat > size(ncseed%x,2)) then
                call realloc(ncseed%x,3,2*ncseed%nat)
                call realloc(ncseed%name,2*ncseed%nat)
                call realloc(ncseed%z,2*ncseed%nat)
             end if
             ncseed%x(:,ncseed%nat) = x
             ncseed%name(ncseed%nat) = c%at(c%atcel(j)%idx)%name
             ncseed%z(ncseed%nat) = c%at(c%atcel(j)%idx)%z 
          end if
       end do
    end do
    call realloc(ncseed%x,3,ncseed%nat)
    call realloc(ncseed%name,ncseed%nat)
    call realloc(ncseed%z,ncseed%nat)
    ncseed%usezname = 3

    if (nr > 0) then
       if (ncseed%nat / c%ncel /= nr) then
          write (uout,*) "c%nneq = ", c%ncel
          write (uout,*) "ncseed%nat = ", ncseed%nat
          write (uout,*) "nr = ", nr
          call ferror('newcell','inconsistent cell # of atoms (nr > 0)',faterr)
       end if
    else
       if (c%ncel / ncseed%nat /= -nr) then
          write (uout,*) "c%nneq = ", c%ncel
          write (uout,*) "ncseed%nat = ", ncseed%nat
          write (uout,*) "nr = ", nr
          call ferror('newcell','inconsistent cell # of atoms (nr < 0)',faterr)
       end if
    endif

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = "<derived>"
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed)
    call c%struct_fill(.true.,-1,.false.,.true.,.false.)
    if (verbose) call c%struct_report()

  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive standard cell. If verbose, write
  !> information about the new crystal. If doforce = .true.,
  !> force the transformation to the primitive even if it does
  !> not lead to a smaller cell.
  subroutine cell_standard(c,toprim,doforce,verbose)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell, spg_get_dataset
    use global, only: symprec
    use tools_math, only: det, matinv
    use tools_io, only: ferror, faterr, uout
    use param, only: maxzat0, eye
    class(crystal), intent(inout) :: c
    logical, intent(in) :: toprim
    logical, intent(in) :: doforce
    logical, intent(in) :: verbose
    
    integer :: ntyp, nat
    integer :: i, iz(maxzat0), id
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types(:)
    real*8 :: rmat(3,3), t(3)

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib transformation to the standard cell
    rmat = transpose(c%crys2car)
    iz = 0
    ntyp = 0
    nat = c%ncel
    allocate(x(3,c%ncel),types(c%ncel))
    do i = 1, c%ncel
       x(:,i) = c%atcel(i)%x
       if (iz(c%at(c%atcel(i)%idx)%z) == 0) then
          ntyp = ntyp + 1
          iz(c%at(c%atcel(i)%idx)%z) = ntyp
          types(i) = ntyp
       else
          types(i) = iz(c%at(c%atcel(i)%idx)%z)
       end if
    end do

    if (toprim) then
       id = spg_standardize_cell(rmat,x,types,nat,1,1,symprec)
       if (id == 0) &
          call ferror("cell_standard","could not find primitive cell",faterr)
       rmat = transpose(rmat)
       do i = 1, 3
          rmat(:,i) = c%c2x(rmat(:,i))
       end do
    else
       id = spg_standardize_cell(rmat,x,types,nat,0,1,symprec)
       if (id == 0) &
          call ferror("cell_standard","could not find standard cell",faterr)
       rmat = transpose(rmat)
       do i = 1, 3
          rmat(:,i) = c%c2x(rmat(:,i))
       end do
    end if

    ! flip the cell?
    if (det(rmat) < 0d0) rmat = -rmat

    ! if a primitive is wanted but det is not less than 1, do not make the change
    if (all(abs(rmat - eye) < symprec)) then
       if (verbose) &
          write (uout,'("+ Cell transformation leads to the same cell: skipping."/)')
       return
    end if
    if (toprim .and. .not.(det(rmat) < 1d0-symprec) .and..not.doforce) then
       if (verbose) &
          write (uout,'("+ Cell transformation does not lead to a smaller cell: skipping."/)')
       return
    end if

    ! transform -> use the origin shift
    t = -matmul(c%spg%origin_shift,rmat)
    ! rmat = transpose(matinv(c%spg%transformation_matrix))
    call c%newcell(rmat,t,verbose)

  end subroutine cell_standard

  !> Transform to the Niggli cell. If verbose, write information
  !> about the new crystal.
  subroutine cell_niggli(c,verbose)
    use spglib, only: spg_niggli_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det
    class(crystal), intent(inout) :: c
    logical, intent(in) :: verbose
    
    real*8 :: rmat(3,3)
    integer :: id, i

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib delaunay reduction
    rmat = transpose(c%crys2car)
    id = spg_niggli_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_niggli","could not find Niggli reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det(rmat) < 0d0) rmat = -rmat

    ! transform
    call c%newcell(rmat,verbose0=verbose)

  end subroutine cell_niggli

  !> Transform to the Delaunay cell. If verbose, write information
  !> about the new crystal.
  subroutine cell_delaunay(c,verbose)
    use spglib, only: spg_delaunay_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det
    class(crystal), intent(inout) :: c
    logical, intent(in) :: verbose

    real*8 :: rmat(3,3)
    integer :: id, i

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib delaunay reduction
    rmat = transpose(c%crys2car)
    id = spg_delaunay_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_delaunay","could not find Delaunay reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det(rmat) < 0d0) rmat = -rmat

    ! transform
    call c%newcell(rmat,verbose0=verbose)

  end subroutine cell_delaunay

  !> Transforms the current basis to the Delaunay reduced basis.
  !> Return the four Delaunay vectors in crystallographic coordinates
  !> (rmat) cell, see 9.1.8 in ITC. If rmati is given, use the three
  !> vectors (cryst. coords.) as the basis for the reduction. If sco
  !> is present, use it in output for the scalar products.
  subroutine delaunay_reduction(c,rmat,rmati,sco)
    class(crystal), intent(in) :: c
    real*8, intent(out) :: rmat(3,4)
    real*8, intent(in), optional :: rmati(3,3)
    real*8, intent(out), optional :: sco(4,4)
    
    integer :: i, j, k
    real*8 :: sc(4,4)
    logical :: again

    real*8, parameter :: eps = 1d-10

    ! build the four Delaunay vectors
    if (present(rmati)) then
       do i = 1, 3
          rmat(:,i) = c%x2c(rmati(:,i))
       end do
    else
       rmat = 0d0
       do i = 1, 3
          rmat(i,i) = 1d0
          rmat(:,i) = c%x2c(rmat(:,i))
       end do
    end if
    rmat(:,4) = -(rmat(:,1)+rmat(:,2)+rmat(:,3))

    ! reduce until all the scalar products are negative or zero
    again = .true.
    sc = -1d0
    do while(again)
       do i = 1, 4
          do j = i+1, 4
             sc(i,j) = dot_product(rmat(:,i),rmat(:,j))
             sc(j,i) = sc(i,j)
          end do
       end do

       if (any(sc > eps)) then
          ai: do i = 1, 4
             aj: do j = i+1, 4
                if (sc(i,j) > eps) exit ai
             end do aj
          end do ai
          do k = 1, 4
             if (i == k .or. j == k) cycle
             rmat(:,k) = rmat(:,i) + rmat(:,k)
          end do
          rmat(:,i) = -rmat(:,i)
       else
          again = .false.
       end if
    end do
    do i = 1, 4
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    if (present(sco)) sco = sc

  end subroutine delaunay_reduction

  !> Create a new, complete crystal/molecule from a crystal seed
  subroutine struct_new(c,seed)
    use global, only: crsmall, atomeps
    use tools_math, only: crys2car_from_cellpar, car2crys_from_cellpar, matinv, &
       det, mnorm2
    use tools_io, only: ferror, faterr, zatguess, string, nameguess
    use types, only: realloc
    use param, only: pi, eyet, ctsq3, maxzat
    class(crystal), intent(inout) :: c
    type(crystalseed), intent(in) :: seed
    
    real*8 :: g(3,3), xmax(3), xmin(3), xcm(3)
    logical :: good, found, hasspg
    integer :: i, j, iat, io, it
    integer :: nnew, icpy
    real*8, allocatable :: atpos(:,:)
    integer, allocatable :: irotm(:), icenv(:)
    real*8 :: v1(3), v2(3)

    if (.not.seed%isused) &
       call ferror("struct_new","uninitialized seed",faterr)

    ! initialize the structure
    call c%init()

    ! copy the atomic information
    c%nneq = seed%nat
    if (c%nneq > 0) then
       if (allocated(c%at)) deallocate(c%at)
       allocate(c%at(c%nneq))
       if (seed%usezname == 1) then
          ! use the atomic number
          do i = 1, c%nneq
             c%at(i)%z = seed%z(i)
             if (c%at(i)%z < 0) &
                call ferror("struct_new","unknown atom with Z: " // string(c%at(i)%z),faterr)
             c%at(i)%name = nameguess(seed%z(i),.true.)
             c%at(i)%x = seed%x(:,i)
          end do
       elseif (seed%usezname == 2) then
          ! use the atomic name
          do i = 1, c%nneq
             c%at(i)%name = seed%name(i)
             c%at(i)%z = zatguess(c%at(i)%name)
             if (c%at(i)%z < 0) &
                call ferror("struct_new","unknown atom: " // string(c%at(i)%name),faterr)
             c%at(i)%x = seed%x(:,i)
          end do
       elseif (seed%usezname == 3) then
          ! use both the atomic number and the name 
          do i = 1, c%nneq
             c%at(i)%name = seed%name(i)
             c%at(i)%z = seed%z(i)
             if (c%at(i)%z < 0) &
                call ferror("struct_new","unknown atom: " // string(c%at(i)%name),faterr)
             c%at(i)%x = seed%x(:,i)
          end do
       else
          call ferror("struct_new","unknown usezname",faterr)
       end if

       ! if this is a molecule, calculate the center and encompassing cell
       if (seed%ismolecule) then
          xmax = -1d40
          xmin =  1d40
          xcm = 0d0
          do i = 1, seed%nat
             do j = 1, 3
                xmax(j) = max(seed%x(j,i)+seed%border,xmax(j))
                xmin(j) = min(seed%x(j,i)-seed%border,xmin(j))
             end do
             xcm = xcm + seed%x(:,i)
          end do
          xcm = xcm / seed%nat
          if (seed%cubic) then
             xmin = minval(xmin)
             xmax = maxval(xmax)
          end if
       end if
    else
       xmax = 1d0
       xcm = 0.5d0
       xmin = 0d0
    end if

    ! basic cell and centering information
    if (seed%useabr == 0) then
       if (.not.seed%ismolecule) &
          call ferror("struct_new","cell data unavailable",faterr)
       ! this is a molecule, for which no cell has been given
       c%aa = xmax - xmin
       c%bb = 90d0
       c%crys2car = crys2car_from_cellpar(c%aa,c%bb)
       c%car2crys = matinv(c%crys2car)
       g = matmul(transpose(c%crys2car),c%crys2car)
    elseif (seed%useabr == 1) then
       ! use aa and bb
       c%aa = seed%aa
       c%bb = seed%bb
       c%crys2car = crys2car_from_cellpar(c%aa,c%bb)
       c%car2crys = matinv(c%crys2car)
       g = matmul(transpose(c%crys2car),c%crys2car)
    elseif (seed%useabr == 2) then
       ! use crys2car
       c%crys2car = seed%crys2car
       c%car2crys = matinv(c%crys2car)
       g = matmul(transpose(c%crys2car),c%crys2car)
       do i = 1, 3
          c%aa(i) = sqrt(g(i,i))
       end do
       c%bb(1) = acos(g(2,3) / c%aa(2) / c%aa(3)) * 180d0 / pi
       c%bb(2) = acos(g(1,3) / c%aa(1) / c%aa(3)) * 180d0 / pi
       c%bb(3) = acos(g(1,2) / c%aa(1) / c%aa(2)) * 180d0 / pi
    else
       call ferror("struct_new","unknown useabr",faterr)
    end if

    ! transform the atomic coordinates in the case of a molecule, and fill 
    ! the remaining molecular fields
    if (seed%ismolecule) then
       c%ismolecule = seed%ismolecule

       if (seed%useabr == 0) then
          ! a cell has not been given

          ! center in the cell and convert to crystallographic coordinates
          do i = 1, c%nneq
             c%at(i)%x = (c%at(i)%x-xcm+0.5d0*c%aa) / c%aa
          end do

          ! Keep the (1/2,1/2,1/2) translation applied
          c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa + xcm 

          ! Set up the molecular cell. c%molborder is in fractional coordinates
          ! and gives the position of the molecular cell in each axis. By default,
          ! choose the molecular cell as the minimal encompassing cell for the molecule
          ! plus 80% of the border or 2 bohr, whichever is larger. The molecular cell 
          ! can not exceed the actual unit cell
          c%molborder = max(seed%border - max(2d0,0.8d0 * seed%border),0d0) / (xmax - xmin)
       else
          if (any(abs(c%bb - 90d0) > 1d-3)) &
             call ferror("struct_new","MOLECULE does not allow non-orthogonal cells",faterr)
          ! a cell has been given, save the origin
          if (seed%havex0) then
             c%molx0 = seed%molx0
          else
             c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa 
          endif

          ! calculate the molecular cell
          if (c%nneq > 0) then
             xmin = 1d40
             do i = 1, c%nneq
                do j = 1, 3
                   xmin(j) = min(c%at(i)%x(j),xmin(j))
                   xmin(j) = min(1d0-max(c%at(i)%x(j),1d0-xmin(j)),xmin(j))
                end do
             end do
          end if
          c%molborder = max(xmin - max(0.8d0 * xmin,2d0/c%aa),0d0)
       end if
    end if

    ! move the crystallographic coordinates to the main cell, calculate the
    ! Cartesian coordinates
    do i = 1, c%nneq
       c%at(i)%x = c%at(i)%x - floor(c%at(i)%x)
       c%at(i)%r = c%x2c(c%at(i)%x)
    end do

    ! rest of the cell metrics
    c%gtensor = g
    c%omega = sqrt(max(det(g),0d0))
    c%grtensor = matinv(g)
    do i = 1, 3
       c%ar(i) = sqrt(c%grtensor(i,i))
    end do
    c%n2_x2c = ctsq3 / mnorm2(c%crys2car)
    c%n2_c2x = ctsq3 / mnorm2(c%car2crys)

    ! calculate the wigner-seitz cell
    call c%wigner((/0d0,0d0,0d0/),nvec=c%nws,vec=c%ivws,&
       nvert_ws=c%nvert_ws,nside_ws=c%nside_ws,iside_ws=c%iside_ws,vws=c%vws)
    c%isortho = (c%nws <= 6)
    if (c%isortho) then
       do i = 1, c%nws
          c%isortho = c%isortho .and. (count(abs(c%ivws(:,i)) == 1) == 1)
       end do
    endif

    ! copy the symmetry information, if available
    if (seed%havesym > 0 .and..not.seed%ismolecule) then
       c%havesym = 1
       c%neqv = seed%neqv
       c%ncv = seed%ncv
       if (allocated(c%cen)) deallocate(c%cen)
       allocate(c%cen(3,seed%ncv))
       c%cen = seed%cen(:,1:seed%ncv)
       c%rotm = 0d0
       c%rotm(:,:,1:seed%neqv) = seed%rotm(:,:,1:seed%neqv)

       ! permute the symmetry operations to make the identity the first
       if (c%neqv > 1) then
          if (.not.all(abs(eyet - c%rotm(:,:,1)) < 1d-12)) then
             good = .false.
             do i = 1, c%neqv
                if (all(abs(eyet - c%rotm(:,:,i)) < 1d-12)) then
                   c%rotm(:,:,i) = c%rotm(:,:,1)
                   c%rotm(:,:,1) = eyet
                   good = .true.
                   exit
                end if
             end do
             if (.not.good) &
                call ferror('struct_new','identity operation not found',faterr)
          end if
       end if
    else
       c%havesym = 0
       c%neqv = 1
       c%rotm = 0d0
       c%rotm(:,:,1) = eyet
       c%ncv = 1
       if (.not.allocated(c%cen)) allocate(c%cen(3,4))
       c%cen = 0d0
    end if

    ! symmetry from spglib
    hasspg = .false.
    if (.not.seed%ismolecule .and. seed%havesym == 0 .and. &
       (seed%findsym == 1 .or. seed%findsym == -1 .and. seed%nat <= crsmall)) then
       ! symmetry was not available, and I want it
       call c%spglib_wrap(.true.,.false.)
       hasspg = .true.
    end if

    ! eliminate redundant atoms 
    nnew = 0
    do i = 1, c%nneq
       found = .false.
       if (c%at(i)%z <= maxzat) then ! skip critical points
          loio: do io = 1, c%neqv
             do it = 1, c%ncv
                v1 = matmul(c%rotm(1:3,1:3,io), c%at(i)%x) + c%rotm(:,4,io) + c%cen(:,it)
                do j = 1, nnew
                   if (c%at(j)%z > maxzat) cycle ! skip critical points
                   v2 = c%at(j)%x
                   if (c%are_lclose(v1,v2,atomeps) .and. c%at(i)%z == c%at(j)%z) then
                      found = .true.
                      icpy = j
                      exit loio
                   end if
                end do
             end do
          end do loio
       end if
       if (.not.found) then
          nnew = nnew + 1
          if (nnew > size(c%at)) then
             call realloc(c%at,2*size(c%at))
          end if
          c%at(nnew) = c%at(i)
       end if
    end do
    c%nneq = nnew
    if (c%nneq > 0) &
       call realloc(c%at,c%nneq)

    ! generate the complete atom list
    if (c%nneq > 0) then
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%nneq*c%neqv*c%ncv))
       if (c%havesym > 0) then
          allocate(atpos(3,192))
          allocate(irotm(192))
          allocate(icenv(192))
          atpos = 0
          irotm = 0
          icenv = 0
          iat = 0
          do i = 1, c%nneq
             call c%symeqv(c%at(i)%x,c%at(i)%mult,atpos,irotm,icenv,atomeps)
             do j = 1, c%at(i)%mult
                iat = iat + 1
                c%atcel(iat)%x = atpos(:,j)
                c%atcel(iat)%r = c%x2c(atpos(:,j))
                c%atcel(iat)%idx = i
                c%atcel(iat)%ir = irotm(j)
                c%atcel(iat)%ic = icenv(j)
                c%atcel(iat)%lvec = nint(atpos(:,j) - &
                   (matmul(c%rotm(1:3,1:3,irotm(j)),c%at(i)%x) + &
                   c%rotm(1:3,4,irotm(j)) + c%cen(:,icenv(j))))
             end do
          end do
          c%ncel = iat
          if (allocated(atpos)) deallocate(atpos)
          if (allocated(irotm)) deallocate(irotm)
          if (allocated(icenv)) deallocate(icenv)
       else
          c%ncel = c%nneq
          do i = 1, c%nneq
             c%at(i)%mult = 1
             c%atcel(i)%x = c%at(i)%x
             c%atcel(i)%r = c%at(i)%r
             c%atcel(i)%idx = i
             c%atcel(i)%ir = 1
             c%atcel(i)%ic = 1
             c%atcel(i)%lvec = 0
          end do
       end if
       call realloc(c%atcel,c%ncel)
    else
       c%ncel = 0
    end if

    ! symmetry is available, but I still want the space group details
    if (c%havesym > 0 .and..not.hasspg) &
       call c%spglib_wrap(.false.,.true.)

    ! the initialization is done - this crystal is ready to use
    c%file = seed%file
    c%isinit = .true.

  end subroutine struct_new

  !> This routine fills ancillary information in the crystal structure
  !> if it is not already available.  If env0, build the atomic
  !> environments. If iast0 = 1, determine the molecular asterisms;
  !> iast0 = 0, do not determine the asterisms; iast0 = -1, only for
  !> small crystals.  If recip0, determine the reciprocal cell metrics
  !> and symmetry.  If lnn0, determine the nearest-neighbor
  !> information. If ewald0, calculate the cutoffs for Ewald method.
  subroutine struct_fill(c,env0,iast0,recip0,lnn0,ewald0)
    use global, only: crsmall

    class(crystal), intent(inout) :: c
    integer :: iast0
    logical, intent(in) :: env0, recip0, lnn0, ewald0

    logical :: env, ast, recip, lnn, ewald
    integer :: i
    real*8, dimension(3) :: vec
    real*8 :: dist(1)
    integer :: nneig(1), wat(1)

    ! Handle input flag dependencies
    env = env0
    if (iast0 == 1) then
       ast = .true.
    elseif (iast0 == 0) then
       ast = .false.
    else
       ast = (c%nneq <= crsmall)
    end if
    recip = recip0
    lnn = lnn0
    ewald = ewald0

    ! nearest-neighbor shells requires environment
    if (lnn) env = .true.

    ! Build the atomic environments
    if (env) then
       call c%build_env()
       c%isenv = .true.
    end if

    ! Reciprocal cell symmetry
    if (recip) then
       ! Reciprocal space point group
       vec = 0d0
       call lattpg(c%car2crys,1,vec,c%neqvg,c%rotg)
       c%isrecip = .true.
    end if

    ! asterisms 
    if (ast) then
       call c%find_asterisms()
       c%isast = .true.
       call c%fill_molecular_fragments()
    end if

    ! nearest neighbors
    if (lnn) then
       do i = 1, c%nneq
          call c%pointshell(c%at(i)%x,1,nneig,wat,dist)
          c%at(i)%rnn2 = dist(1) / 2d0
       end do
       c%isnn = .true.
    end if

    ! preparation for ewald 
    if (ewald) then
       call c%calculate_ewald_cutoffs()
       c%isewald = .true.
    end if

  end subroutine struct_fill

  !> Write information about the crystal structure to the output.
  subroutine struct_report(c)
    use fragmentmod, only: fragment_cmass
    use global, only: iunitname0, dunit0, iunit
    use tools_math, only: gcd, norm
    use tools_io, only: uout, string, ioj_center, ioj_left, ioj_right
    use param, only: bohrtoa, maxzat
    class(crystal), intent(in) :: c

    integer, parameter :: natenvmax = 2000

    integer :: i, j, k
    integer :: nelec, holo, laue
    real*8 :: maxdv, xcm(3), x0(3)
    character(len=:), allocatable :: str1, str2
    character(len=3) :: schpg
    integer, allocatable :: nneig(:), wat(:)
    real*8, allocatable :: dist(:)

    ! Header
    if (.not.c%ismolecule) then
       write (uout,'("* Crystal structure")')
       write (uout,'("  From: ",A)') string(c%file)
       write (uout,'("  Lattice parameters (bohr): ",3(A,2X))') &
          string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
       write (uout,'("  Lattice parameters (ang): ",3(A,2X))') &
          string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
       write (uout,'("  Lattice angles (degrees): ",3(A,2X))') &
          string(c%bb(1),'f',decimal=3), string(c%bb(2),'f',decimal=3), string(c%bb(3),'f',decimal=3)
    else
       write (uout,'("* Molecular structure")')
       write (uout,'("  From: ",A)') string(c%file)
       write (uout,'("  Encompassing cell dimensions (bohr): ",3(A,2X))') &
          string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
       write (uout,'("  Encompassing cell dimensions (ang): ",3(A,2X))') &
          string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
    endif

    ! Compute unit formula, and z
    if (.not.c%ismolecule) then
       maxdv = gcd(c%at(1:c%nneq)%mult,c%nneq)
       write (uout,'("  Molecular formula: ",999(/4X,10(A,"(",A,") ")))') &
          (string(c%at(i)%name), string(nint(c%at(i)%mult/maxdv)), i=1,c%nneq)
       write (uout,'("  Number of non-equivalent atoms in the unit cell: ",A)') string(c%nneq)
       write (uout,'("  Number of atoms in the unit cell: ",A)') string(c%ncel)
    else
       write (uout,'("  Number of atoms: ",A)') string(c%ncel)
    endif
    nelec = 0
    do i = 1, c%nneq
       if (c%at(i)%z >= maxzat) cycle
       nelec = nelec + c%at(i)%z * c%at(i)%mult
    end do
    write (uout,'("  Number of electrons: ",A/)') string(nelec)
    
    ! List of atoms in crystallographic coordinates
    if (.not.c%ismolecule) then
       write (uout,'("+ List of non-equivalent atoms in the unit cell (cryst. coords.): ")')
       write (uout,'("# ",7(A,X))') string("nat",3,ioj_center), &
          string("x",14,ioj_center), string("y",14,ioj_center),&
          string("z",14,ioj_center), string("name",10,ioj_center), &
          string("mult",4,ioj_center), string("Z",4,ioj_center)
       do i=1, c%nneq
          write (uout,'(2x,7(A,X))') &
             string(i,3,ioj_center),&
             string(c%at(i)%x(1),'f',length=14,decimal=10,justify=3),&
             string(c%at(i)%x(2),'f',length=14,decimal=10,justify=3),&
             string(c%at(i)%x(3),'f',length=14,decimal=10,justify=3),& 
             string(c%at(i)%name,10,ioj_center), &
             string(c%at(i)%mult,4,ioj_center), string(c%at(i)%z,4,ioj_center)
       enddo
       write (uout,*)

       write (uout,'("+ List of atoms in the unit cell (cryst. coords.): ")')
       write (uout,'("# ",6(A,X))') string("at",3,ioj_center),&
          string("x",14,ioj_center), string("y",14,ioj_center),&
          string("z",14,ioj_center), string("name",10,ioj_center),&
          string("Z",4,ioj_center)
       do i=1,c%ncel
          write (uout,'(2x,6(A,X))') &
             string(i,3,ioj_center),&
             string(c%atcel(i)%x(1),'f',length=14,decimal=10,justify=3),&
             string(c%atcel(i)%x(2),'f',length=14,decimal=10,justify=3),&
             string(c%atcel(i)%x(3),'f',length=14,decimal=10,justify=3),& 
             string(c%at(c%atcel(i)%idx)%name,10,ioj_center),&
             string(c%at(c%atcel(i)%idx)%z,4,ioj_center)
       enddo
       write (uout,*)

       write (uout,'("+ Lattice vectors (",A,")")') iunitname0(iunit)
       do i = 1, 3
          write (uout,'(4X,A,": ",3(A,X))') string(i), (string(c%crys2car(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
       end do
       write (uout,*)
    end if

    ! List of atoms in Cartesian coordinates
    write (uout,'("+ List of atoms in Cartesian coordinates (",A,"): ")') iunitname0(iunit)
    write (uout,'("# ",6(A,X))') string("at",3,ioj_center), &
       string("x",16,ioj_center), string("y",16,ioj_center),&
       string("z",16,ioj_center), string("name",10,ioj_center),&
       string("Z",4,ioj_center)
    do i=1,c%ncel
       write (uout,'(2x,6(A,X))') &
          string(i,3,ioj_center),&
          (string((c%atcel(i)%r(j)+c%molx0(j))*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
          string(c%at(c%atcel(i)%idx)%name,10,ioj_center),&
          string(c%at(c%atcel(i)%idx)%z,4,ioj_center)
    enddo
    write (uout,*)

    ! Encompassing region for the molecule
    if (c%ismolecule) then
       write (uout,'("+ Limits of the molecular cell (in fractions of the encompassing cell).")')
       write (uout,'("  The region of the encompassing cell outside the molecular cell is")')
       write (uout,'("  assumed to represent infinity (no CPs or gradient paths in it).")')
       write (uout,'("  x-axis: ",A," -> ",A)') trim(string(c%molborder(1),'f',10,4)), trim(string(1d0-c%molborder(1),'f',10,4))
       write (uout,'("  y-axis: ",A," -> ",A)') trim(string(c%molborder(2),'f',10,4)), trim(string(1d0-c%molborder(2),'f',10,4))
       write (uout,'("  z-axis: ",A," -> ",A)') trim(string(c%molborder(3),'f',10,4)), trim(string(1d0-c%molborder(3),'f',10,4))
       write (uout,*)
    end if

    ! Cell volume
    if (.not.c%ismolecule) then
       write (uout,'("+ Cell volume (bohr^3): ",A)') string(c%omega,'f',decimal=5)
       write (uout,'("+ Cell volume (ang^3): ",A)') string(c%omega * bohrtoa**3,'f',decimal=5)
       write (uout,*)
    end if
    
    ! Write symmetry operations 
    if (.not.c%ismolecule) then
       write(uout,'("+ List of symmetry operations (",A,"):")') string(c%neqv)
       do k = 1, c%neqv
          write (uout,'(2X,"Operation ",A,":")') string(k)
          write (uout,'(2(4X,4(A,X)/),4X,4(A,X))') ((string(c%rotm(i,j,k),'f',length=9,decimal=6,justify=3), j = 1, 4), i = 1, 3)
       enddo
       write (uout,*)
    
       call c%struct_report_symxyz()

       write(uout,'("+ List of centering vectors (",A,"):")') string(c%ncv)
       do k = 1, c%ncv
          write (uout,'(2X,"Vector ",A,": ",3(A,X))') string(k), &
             (string(c%cen(i,k),'f',length=9,decimal=6), i = 1, 3)
       enddo
       write (uout,*)
    
       if (c%havesym > 0) then
          write(uout,'("+ Crystal symmetry information")')
          if (c%spg%n_atoms > 0) then
             if (len_trim(c%spg%choice) > 0) then
                write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,", setting ",A,")")') &
                   string(c%spg%international_symbol), string(c%spg%spacegroup_number),&
                   string(c%spg%choice)
             else
                write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,")")') &
                   string(c%spg%international_symbol), string(c%spg%spacegroup_number)
             end if
             write(uout,'("  Space group (Hall): ",A, " (number ",A,")")') &
                string(c%spg%hall_symbol), string(c%spg%hall_number)
             write(uout,'("  Point group (Hermann-Mauguin): ",A)') string(c%spg%pointgroup_symbol)

             call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
             write(uout,'("  Point group (Schoenflies): ",A)') string(schpg)
             write(uout,'("  Holohedry: ",A)') string(holo_string(holo))
             write(uout,'("  Laue class: ",A)') string(laue_string(laue))
          else
             write(uout,'("  Unavailable because symmetry read from external file")')
          end if
          write (uout,*)
       end if

       write (uout,'("+ Cartesian/crystallographic coordinate transformation matrices:")')
       write (uout,'("  A = car to crys (xcrys = A * xcar, ",A,"^-1)")') iunitname0(iunit)
       do i = 1, 3
          write (uout,'(4X,3(A,X))') (string(c%car2crys(i,j)/dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
       end do
       write (uout,'("  B = crys to car (xcar = B * xcrys, ",A,")")') iunitname0(iunit)
       do i = 1, 3
          write (uout,'(4X,3(A,X))') (string(c%crys2car(i,j)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
       end do
       write (uout,'("  G = metric tensor (B''*B, ",A,"^2)")') iunitname0(iunit)
       do i = 1, 3
          write (uout,'(4X,3(A,X))') (string(c%gtensor(i,j)*dunit0(iunit)**2,'f',length=16,decimal=10,justify=5),j=1,3)
       end do
       write (uout,*)
    end if

    ! Discrete molecules, if available
    if (allocated(c%nstar) .and. allocated(c%mol) .and. c%nmol > 0) then
       write (uout,'("+ List of fragments in the system")')
       write (uout,'("  Number of fragments: ",A)') string(c%nmol)
       write (uout,'("# Id  nat           Center of mass          Discrete? ")')
       do i = 1, c%nmol
          if (c%ismolecule) then
             xcm = (fragment_cmass(c%mol(i))+c%molx0) * dunit0(iunit)
          else
             xcm = c%c2x(fragment_cmass(c%mol(i)))
          end if
          write (uout,'(99(2X,A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
             (string(xcm(j),'f',10,6,3),j=1,3), string(c%moldiscrete(i))
       end do
       if (.not.c%ismolecule .and. all(c%moldiscrete(1:c%nmol))) &
          write (uout,'(/"+ This is a molecular crystal.")')
       write (uout,*)
    end if

    ! Number of atoms in the atomic environment
    if (.not.c%ismolecule) then
       write (uout,'("+ Atomic environment of the main cell")')
       write (uout,'("  Number of atoms contributing density to the main cell: ",A)') string(c%nenv)
       write (uout,*)
    end if

    ! Print out atomic environments and determine the nearest neighbor distances
    if (c%nneq <= natenvmax) then
       write (uout,'("+ Atomic environments (distances in ",A,")")') iunitname0(iunit)
       write (uout,'("# ",6(A,2X))') &
          string("id",length=4,justify=ioj_center), &
          string("atom",length=5,justify=ioj_center), &
          string("nneig",length=5,justify=ioj_center), &
          string("distance",length=11,justify=ioj_right), &
          string("nat",length=4,justify=ioj_center), &
          string("type",length=10,justify=ioj_left)
       allocate(nneig(10),wat(10),dist(10))
       do i = 1, c%nneq
          call c%pointshell(c%at(i)%x,10,nneig,wat,dist)
          do j = 1, 10
             if (j == 1) then
                str1 = string(i,length=4,justify=ioj_center)
                str2 = string(c%at(i)%name,length=5,justify=ioj_center)
             else
                str1 = string("",length=4,justify=ioj_center)
                str2 = " ... "
             end if
             if (wat(j) /= 0) then
                write (uout,'(6(2X,A))') &
                   str1, str2, &
                   string(nneig(j),length=5,justify=ioj_center), &
                   string(dist(j)*dunit0(iunit),'f',length=12,decimal=7,justify=5), &
                   string(wat(j),length=4,justify=ioj_center), &
                   string(c%at(wat(j))%name,length=10,justify=ioj_left)
             end if
          end do
       end do
       write (uout,*)
       deallocate(nneig,wat,dist)
    else
       write (uout,'("+ Atomic environments not written because of the large number ")')
       write (uout,'("  of non-equivalent atoms (",A," > ",A,"). Please, use the")') string(c%nneq), string(natenvmax)
       write (uout,'("  ENVIRON keyword to calculate the atomic environments.")')
    end if

    ! Determine nn/2 for every atom
    write (uout,'("+ List of half nearest neighbor distances (",A,")")') iunitname0(iunit)
    write (uout,'(3(2X,A))') string("id",length=4,justify=ioj_center),&
       string("atom",length=5,justify=ioj_center), &
       string("rnn/2",length=12,justify=ioj_center)
    do i = 1, c%nneq
       write (uout,'(3(2X,A))') string(i,length=4,justify=ioj_center),&
          string(c%at(i)%name,length=5,justify=ioj_center), &
          string(c%at(i)%rnn2*dunit0(iunit),'f',length=12,decimal=7,justify=4)
    end do
    write (uout,*)

    ! Determine the wigner-seitz cell
    if (.not.c%ismolecule) then
       write (uout,'("+ Vertex of the WS cell (cryst. coords.)")')
       write (uout,'(5(2X,A))') string("id",length=3,justify=ioj_right),&
          string("x",length=11,justify=ioj_center),&
          string("y",length=11,justify=ioj_center),&
          string("z",length=11,justify=ioj_center),&
          string("d ("//iunitname0(iunit)//")",length=14,justify=ioj_center)
       do i = 1, c%nvert_ws
          x0 = c%x2c(c%vws(:,i))
          write (uout,'(5(2X,A))') string(i,length=3,justify=ioj_right), &
             (string(c%vws(j,i),'f',length=11,decimal=6,justify=4),j=1,3), &
             string(norm(x0)*dunit0(iunit),'f',length=14,decimal=8,justify=4)
       enddo
       write (uout,*)

       write (uout,'("+ Faces of the WS cell")')
       write (uout,'("  Number of faces: ",A)') string(c%nws)
       do i = 1, c%nws
          write (uout,'(2X,A,": ",999(A,X))') string(i,length=2,justify=ioj_right), &
             (string(c%iside_ws(j,i),length=2),j=1,c%nside_ws(i))
       end do
       write (uout,*)

       write (uout,'("+ Lattice vectors for the Wigner-Seitz neighbors")')
       do i = 1, c%nws
          write (uout,'(2X,A,": ",99(A,X))') string(i,length=2,justify=ioj_right), &
             (string(c%ivws(j,i),length=2,justify=ioj_right),j=1,size(c%ivws,1))
       end do
       write (uout,*)
       write (uout,'("+ Is the cell orthogonal? ",L1/)') c%isortho
    end if

  end subroutine struct_report

  !> Write the list of symmetry operations to stdout, using crystallographic
  !> notation (if possible).
  subroutine struct_report_symxyz(c,strfin)
    use tools_io, only: uout, string
    use global, only: symprec
    class(crystal), intent(in) :: c
    character*255, intent(out), optional :: strfin(c%neqv)

    real*8, parameter :: rfrac(25) = (/-12d0/12d0,-11d0/12d0,-10d0/12d0,&
       -9d0/12d0,-8d0/12d0,-7d0/12d0,-6d0/12d0,-5d0/12d0,-4d0/12d0,-3d0/12d0,&
       -2d0/12d0,-1d0/12d0,0d0/12d0,1d0/12d0,2d0/12d0,3d0/12d0,4d0/12d0,&
       5d0/12d0,6d0/12d0,7d0/12d0,8d0/12d0,9d0/12d0,10d0/12d0,11d0/12d0,12d0/12d0/)
    character*6, parameter :: sfrac(25) = (/"      ","-11/12","-5/6  ",&
       "-3/4  ","-2/3  ","-7/12 ","-1/2  ","-5/12 ","-1/3  ","-1/4  ","-1/6  ",&
       "-1/12 ","      ","1/12  ","1/6   ","1/4   ","1/3   ","5/12  ","1/2   ",&
       "7/12  ","2/3   ","3/4   ","5/6   ","11/12 ","      "/)
    character*1, parameter :: xyz(3) = (/"x","y","z"/)

    logical :: ok, iszero
    integer :: i, j, k
    character*255 :: strout(c%neqv)

    do i = 1, c%neqv
       strout(i) = ""
       do j = 1, 3
          ! translation
          ok = .false.
          do k = 1, 25
             if (abs(c%rotm(j,4,i) - rfrac(k)) < symprec) then
                ok = .true.
                strout(i) = trim(strout(i)) // sfrac(k)
                iszero = (k == 13) .or. (k == 1) .or. (k == 25)
                exit
             end if
          end do
          if (.not.ok) return

          ! rotation
          do k = 1, 3
             if (abs(c%rotm(j,k,i) - 1d0) < symprec) then
                if (iszero) then
                   strout(i) = trim(strout(i)) // xyz(k)
                else
                   strout(i) = trim(strout(i)) // "+" // xyz(k)
                end if
                iszero = .false.
             elseif (abs(c%rotm(j,k,i) + 1d0) < symprec) then
                strout(i) = trim(strout(i)) // "-" // xyz(k)
                iszero = .false.
             elseif (abs(c%rotm(j,k,i)) > symprec) then
                return
             end if
          end do

          ! the comma
          if (j < 3) &
             strout(i) = trim(strout(i)) // ","
       end do
    end do

    if (present(strfin)) then
       strfin = strout
    else
       write(uout,'("+ List of symmetry operations in crystallographic notation:")')
       do k = 1, c%neqv
          write (uout,'(3X,A,": ",A)') string(k), string(strout(k))
       enddo
       write (uout,*)
    end if

  end subroutine struct_report_symxyz

  !> Use the spg library to find information about the space group.
  !> In: cell vectors (crys2car), ncel, atcel(:), at(:) Out: neqv,
  !> rotm, spg. If usenneq is .true., use nneq and at(:) instead of 
  !> ncel and atcel. If onlyspg is .true., fill only the spg field
  !> and leave the others unchanged.
  subroutine spglib_wrap(c,usenneq,onlyspg)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_dataset, spg_get_error_message
    use global, only: symprec
    use tools_io, only: string, ferror, warning, equal
    use param, only: maxzat0, eyet, eye
    use types, only: realloc
    class(crystal), intent(inout) :: c
    logical, intent(in) :: usenneq
    logical, intent(in) :: onlyspg

    real(c_double) :: lattice(3,3)
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: typ(:)
    integer :: ntyp, nat
    integer :: i, j, iz(maxzat0)
    character(len=32) :: error
    logical :: found
    real*8 :: rotm(3,3)

    ! get the dataset from spglib
    lattice = transpose(c%crys2car)
    iz = 0
    ntyp = 0
    if (usenneq) then
       nat = c%nneq
       allocate(x(3,c%nneq),typ(c%nneq))
       do i = 1, c%nneq
          x(:,i) = c%at(i)%x
          if (iz(c%at(i)%z) == 0) then
             ntyp = ntyp + 1
             iz(c%at(i)%z) = ntyp
             typ(i) = ntyp
          else
             typ(i) = iz(c%at(i)%z)
          end if
       end do
    else
       nat = c%ncel
       allocate(x(3,c%ncel),typ(c%ncel))
       do i = 1, c%ncel
          x(:,i) = c%atcel(i)%x
          if (iz(c%at(c%atcel(i)%idx)%z) == 0) then
             ntyp = ntyp + 1
             iz(c%at(c%atcel(i)%idx)%z) = ntyp
             typ(i) = ntyp
          else
             typ(i) = iz(c%at(c%atcel(i)%idx)%z)
          end if
       end do
    end if
    c%spg = spg_get_dataset(lattice,x,typ,nat,symprec)
    deallocate(x,typ)

    ! check error messages
    error = trim(spg_get_error_message(c%spg%spglib_error))
    if (.not.equal(error,"no error")) &
       call ferror("spglib_wrap","error from spglib: "//string(error),warning)

    if (onlyspg) return

    ! unpack spglib's output into pure translations and symops
    c%neqv = 1
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (.not.allocated(c%cen)) allocate(c%cen(3,1))
    c%cen = 0d0
    do i = 1, c%spg%n_operations
       rotm = transpose(c%spg%rotations(:,:,i))

       ! is this a pure translation?
       if (all(abs(rotm - eye) < symprec)) then
          found = .false.
          do j = 1, c%ncv
             if (all(abs(c%spg%translations(:,i) - c%cen(:,j)) < symprec)) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) then
             c%ncv = c%ncv + 1
             if (c%ncv > size(c%cen,2)) &
                call realloc(c%cen,3,2*c%ncv)
             c%cen(:,c%ncv) = c%spg%translations(:,i)
          end if
          cycle
       end if

       ! Do I have this rotation already?
       found = .false.
       do j = 1, c%neqv
          if (all(abs(rotm - c%rotm(1:3,1:3,j)) < symprec)) then
             found = .true.
          end if
       end do
       if (.not.found) then
          c%neqv = c%neqv + 1
          c%rotm(1:3,1:3,c%neqv) = rotm
          c%rotm(1:3,4,c%neqv) = c%spg%translations(:,i)
       end if
    end do
    call realloc(c%cen,3,c%ncv)
    c%havesym = 1

  end subroutine spglib_wrap

  !xx! Wigner-Seitz cell tools and cell partition

  !> Builds the Wigner-Seitz cell and its irreducible wedge.
  subroutine wigner(c,xorigin,nvec,vec,area0,ntetrag,tetrag,&
     nvert_ws,nside_ws,iside_ws,vws)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_int
    use global, only: fileroot
    use tools_math, only: norm, mixed, cross
    use tools_io, only: string, filepath, fopen_write, fopen_read,&
       ferror, faterr, fclose
    use param, only: dirsep
    use types, only: realloc

    interface
       subroutine doqhull(fin,fvert,fface,ithr) bind(c)
         use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_int
         character(kind=c_char) :: fin(*), fvert(*), fface(*)
         integer(kind=c_int) :: ithr
       end subroutine doqhull
    end interface

    class(crystal), intent(in) :: c
    real*8, intent(in) :: xorigin(3) !< Origin of the WS cell
    integer, intent(out), optional :: nvec !< Number of lattice point neighbors
    integer, intent(out), optional :: vec(3,16) !< Integer vectors to neighbors
    real*8, intent(out), optional :: area0(40) !< Area to neighbors
    integer, intent(out), optional :: ntetrag !< number of tetrahedra forming the irreducible WS cell
    real*8, allocatable, intent(inout), optional :: tetrag(:,:,:) !< vertices of the tetrahedra
    integer, intent(out), optional :: nvert_ws !< number of vertices of the WS cell
    integer, allocatable, intent(inout), optional :: nside_ws(:) !< number of sides of WS faces
    integer, allocatable, intent(inout), optional :: iside_ws(:,:) !< sides of the WS faces
    real*8, allocatable, intent(inout), optional :: vws(:,:) !< vertices of the WS cell

    ! three WS vertices are collinear if det < -eps
    real*8, parameter :: eps_wspesca = 1d-5 !< Criterion for tetrahedra equivalence
    real*8, parameter :: ws_eps_vol = 1d-5 !< Reject tetrahedra smaller than this.
    real*8, parameter :: eps_bary = 1d-1 !< barycenter identification
    integer, parameter :: icelmax_def = 2 !< maximum number of cells in a direction.
    integer, parameter :: icelmax_safe = 4 !< maximum number of cells in a direction (safe)
    integer, parameter :: icelmax_usafe = 10 !< maximum number of cells in a direction (ultra-safe)

    real*8 :: rnorm
    integer :: i, j, k, n, npolig, leqv, icelmax
    real*8 :: xp1(3), xp2(3), xp3(3), x0(3)
    logical, allocatable :: active(:)
    real*8, allocatable :: tvol(:), xstar(:,:)
    real*8 :: lrotm(3,3,48), sumi, xoriginc(3)
    real*8 :: area, bary(3), av(3)
    real*8 :: x2r(3,3), r2x(3,3)
    real*8, allocatable :: xws(:,:)
    integer :: nvert, iaux, lu
    integer, allocatable :: nside(:), iside(:,:)
    character(len=:), allocatable :: file1, file2, file3
    character(kind=c_char,len=1024) :: file1c, file2c, file3c
    character*3 :: pg
    real*8 :: rmat(3,4), rmati(3,3)
    ! qhull threshold for face recognition
    integer(c_int) :: ithr
    integer(c_int), parameter :: ithr_def = 5

    ! set origin
    xoriginc = c%x2c(xorigin)

    ! input for delaunay
    rmati = 0d0
    do i = 1, 3
       rmati(i,i) = 1d0
    end do

    ! cartesian/crystallographic
    x2r = c%crys2car
    r2x = c%car2crys
    rnorm = 1d0

    ! anchor for when critic2 and qhull fight each other
    ithr = ithr_def
    icelmax = icelmax_def
99  continue

    ! construct star of lattice vectors -> use Delaunay reduction
    ! see 9.1.8 in ITC.
    n = 14
    allocate(xstar(3,14))
    call c%delaunay_reduction(rmat,rmati)
    xstar(:,1)  = rmat(:,1)
    xstar(:,2)  = rmat(:,2)
    xstar(:,3)  = rmat(:,3)
    xstar(:,4)  = rmat(:,4)
    xstar(:,5)  = rmat(:,1)+rmat(:,2)
    xstar(:,6)  = rmat(:,1)+rmat(:,3)
    xstar(:,7)  = rmat(:,2)+rmat(:,3)
    xstar(:,8)  = -(rmat(:,1))
    xstar(:,9)  = -(rmat(:,2))
    xstar(:,10) = -(rmat(:,3))
    xstar(:,11) = -(rmat(:,4))
    xstar(:,12) = -(rmat(:,1)+rmat(:,2))
    xstar(:,13) = -(rmat(:,1)+rmat(:,3))
    xstar(:,14) = -(rmat(:,2)+rmat(:,3))
    do i = 1, 14
       xstar(:,i) = matmul(x2r,xstar(:,i))
    end do

    ! compute the voronoi polyhedron using libqhull
    file1 = trim(filepath) // dirsep // trim(fileroot) // "_wsstar.dat"
    file1c = trim(file1) // C_NULL_CHAR
    file2 = trim(filepath) // dirsep // trim(fileroot) // "_wsvert.dat"
    file2c = trim(file2) // C_NULL_CHAR
    file3 = trim(filepath) // dirsep // trim(fileroot) // "_wsface.dat"
    file3c = trim(file3) // C_NULL_CHAR

    ! write the file with the star vertices
    lu = fopen_write(file1,abspath0=.true.)
    write (lu,'("3"/I4)') n+1
    write (lu,'(3(F20.12,X))') 0d0,0d0,0d0
    do i = 1, n
       write (lu,'(3(F20.12,X))') xstar(:,i)
    end do
    call fclose(lu)
    deallocate(xstar)

    ! run qhull
    call doqhull(file1c,file2c,file3c,ithr)

    ! read the vertices from file number 2
    lu = fopen_read(file2,abspath0=.true.)
    read(lu,*)
    read(lu,*) nvert
    allocate(xws(3,nvert))
    do i = 1, nvert
       read(lu,*) xws(:,i)
    end do
    call fclose(lu)

    ! first pass, the number of sides of each polygon
    lu = fopen_read(file3,abspath0=.true.)
    read(lu,*)
    read(lu,*) iaux, npolig
    allocate(nside(npolig))
    do i = 1, nvert
       read(lu,*)
    end do
    do i = 1, npolig
       read(lu,*) nside(i)
    end do

    ! second pass, the side indices
    rewind(lu)
    read(lu,*)
    read(lu,*)
    do i = 1, nvert
       read(lu,*)
    end do
    allocate(iside(npolig,maxval(nside)))
    do i = 1, npolig
       read(lu,*) iaux, (iside(i,j),j=1,nside(i))
    end do
    iside = iside + 1

    call fclose(lu)

    ! save faces and vertices
    if (present(nside_ws)) then
       if (allocated(nside_ws)) deallocate(nside_ws)
       allocate(nside_ws(npolig))
       nside_ws = nside(1:npolig)
    end if
    if (present(iside_ws)) then
       if(allocated(iside_ws)) deallocate(iside_ws)
       allocate(iside_ws(maxval(nside(1:npolig)),npolig))
       iside_ws = 0
       do i = 1, npolig
          do j = 1, nside(i)
             iside_ws(j,i) = iside(i,j)
          end do
       end do
    end if
    if (present(nvert_ws)) then
       nvert_ws = nvert
    end if
    if (present(vws)) then
       if (allocated(vws)) deallocate(vws)
       allocate(vws(3,nvert))
       do i = 1, nvert
          x0 = matmul(r2x,xws(:,i))
          vws(:,i) = x0
       end do
    end if

    ! tetrahedra
    if (present(ntetrag).and.present(tetrag)) then
       ! local symmetry group
       pg = c%sitesymm(xorigin,leqv=leqv,lrotm=lrotm)

       ! build all the tetrahedra
       ntetrag = 0
       if (allocated(tetrag)) deallocate(tetrag)
       allocate(tetrag(4,3,10))
       do i = 1, npolig
          n = nside(i)
          ! calculate middle point of the face
          x0 = 0d0
          do j = 1, n
             x0 = x0 + xws(:,iside(i,j))
          end do
          x0 = x0 / n

          do j = 1, n
             xp1 = xws(:,iside(i,j))
             xp2 = xws(:,iside(i,mod(j,n)+1))
             if (ntetrag+2>size(tetrag,3)) &
                call realloc(tetrag,4,3,2*size(tetrag,3))
             ! tetrah 1
             ntetrag = ntetrag + 1
             tetrag(1,:,ntetrag) = (/ 0d0, 0d0, 0d0 /) + xoriginc
             tetrag(2,:,ntetrag) = x0 + xoriginc
             tetrag(3,:,ntetrag) = xp1 + xoriginc
             tetrag(4,:,ntetrag) = 0.5d0 * (xp1 + xp2) + xoriginc
             ! tetrah 2
             ntetrag = ntetrag + 1
             tetrag(1,:,ntetrag) = (/ 0d0, 0d0, 0d0 /) + xoriginc
             tetrag(2,:,ntetrag) = x0 + xoriginc
             tetrag(3,:,ntetrag) = xp2 + xoriginc
             tetrag(4,:,ntetrag) = 0.5d0 * (xp1 + xp2) + xoriginc
          end do
       end do

       ! output some info
       sumi = 0d0
       allocate(active(ntetrag),tvol(ntetrag))
       active = .true.
       do i = 1, ntetrag
          ! calculate volume
          xp1 = tetrag(2,:,i) - tetrag(1,:,i)
          xp2 = tetrag(3,:,i) - tetrag(1,:,i)
          xp3 = tetrag(4,:,i) - tetrag(1,:,i)
          tvol(i) = abs(mixed(xp1,xp2,xp3)) / 6d0
          if (tvol(i) < ws_eps_vol) then
             active(i) = .false.
             cycle
          end if
          sumi = sumi + tvol(i)
          do j = 1, 4
             tetrag(j,:,i) = matmul(r2x,tetrag(j,:,i))
          end do
       end do

       ! reduce tetrahedra equivalent by symmetry
       do i = 1, ntetrag
          if (.not.active(i)) cycle
          do j = i+1, ntetrag
             if (equiv_tetrah(c,xorigin,tetrag(:,:,i),tetrag(:,:,j),leqv,lrotm,eps_wspesca)) then
                active(j) = .false.
             end if
          end do
       end do

       ! rebuild the tetrahedra list
       in: do i = 1, ntetrag
          if (active(i)) cycle
          do j = i+1, ntetrag
             if (active(j)) then
                tetrag(:,:,i) = tetrag(:,:,j)
                tvol(i) = tvol(j)
                active(j) = .false.
                cycle in
             end if
          end do
          ntetrag = i - 1
          exit
       end do in
       call realloc(tetrag,4,3,ntetrag)

       ! renormalize
       tvol = tvol

       deallocate(active,tvol)
    end if

    if (present(nvec)) then
       if (.not.present(vec)) &
          call ferror('wigner','incorrect call',faterr)
       ! calculate neighbors and areas
       nvec = 0
       do i = 1, npolig
          ! lattice point
          bary = 0d0
          do j = 1, nside(i)
             bary = bary + xws(:,iside(i,j))
          end do
          bary = 2d0 * bary / nside(i)

          ! area of a convex polygon
          av = 0d0
          do j = 1, nside(i)
             k = mod(j,nside(i))+1
             av = av + cross(xws(:,iside(i,j)),xws(:,iside(i,k)))
          end do
          area = 0.5d0 * abs(dot_product(bary,av) / norm(bary)) / rnorm**2
          bary = matmul(r2x,bary)

          nvec = nvec + 1
          vec(:,nvec) = nint(bary)
          if (any(bary - nint(bary) > eps_bary)) then
             ithr = ithr - 1
             if (ithr > 1) then
                ! write (uout,'("(!!) Restarting with ithr = ",I2," icelmax = ",I2)') ithr, icelmax
                ! call ferror("wigner","wrong lattice vector",warning)
             else
                if (icelmax == icelmax_def) then
                   icelmax = icelmax_safe
                   ithr = ithr_def
                   ! write (uout,'("(!!) Restarting with ithr = ",I2," icelmax = ",I2)') ithr, icelmax
                   ! call ferror("wigner","wrong lattice vector",warning)
                elseif (icelmax == icelmax_safe) then
                   icelmax = icelmax_usafe
                   ithr = ithr_def
                   ! write (uout,'("(!!) Restarting with ithr = ",I2," icelmax = ",I2)') ithr, icelmax
                   ! call ferror("wigner","wrong lattice vector",warning)
                else
                   call ferror("wigner","wrong lattice vector",faterr)
                end if
             end if
             call unlink(file1)
             call unlink(file2)
             call unlink(file3)
             deallocate(xws,iside,nside)
             goto 99
          end if
          if (present(area0)) then
             area0(nvec) = area
          endif
       end do
    end if

    ! remove the temporary files and cleanup
    call unlink(file1)
    call unlink(file2)
    call unlink(file3)
    deallocate(xws,iside,nside)

  end subroutine wigner

  !> Partition the unit cell in tetrahedra.
  subroutine pmwigner(c,ntetrag,tetrag)
    use tools_math, only: mixed
    class(crystal), intent(in) :: c !< the crystal structure
    integer, intent(out), optional :: ntetrag !< number of tetrahedra forming the irreducible WS cell
    real*8, allocatable, intent(out), optional :: tetrag(:,:,:) !< vertices of the tetrahedra

    real*8 :: sumi, tvol(5), xp1(3), xp2(3), xp3(3)
    integer :: i

    if (allocated(tetrag)) deallocate(tetrag)
    allocate(tetrag(4,3,5))
    ntetrag = 5
    tetrag(1,:,1) = (/ 1d0, 0d0, 0d0 /)
    tetrag(2,:,1) = (/ 1d0, 0d0, 1d0 /)
    tetrag(3,:,1) = (/ 1d0, 1d0, 0d0 /)
    tetrag(4,:,1) = (/ 0d0, 0d0, 0d0 /)

    tetrag(1,:,2) = (/ 0d0, 1d0, 0d0 /)
    tetrag(2,:,2) = (/ 0d0, 1d0, 1d0 /)
    tetrag(3,:,2) = (/ 0d0, 0d0, 0d0 /)
    tetrag(4,:,2) = (/ 1d0, 1d0, 0d0 /)

    tetrag(1,:,3) = (/ 0d0, 0d0, 1d0 /)
    tetrag(2,:,3) = (/ 0d0, 0d0, 0d0 /)
    tetrag(3,:,3) = (/ 0d0, 1d0, 1d0 /)
    tetrag(4,:,3) = (/ 1d0, 0d0, 1d0 /)

    tetrag(1,:,4) = (/ 1d0, 1d0, 1d0 /)
    tetrag(2,:,4) = (/ 1d0, 1d0, 0d0 /)
    tetrag(3,:,4) = (/ 1d0, 0d0, 1d0 /)
    tetrag(4,:,4) = (/ 0d0, 1d0, 1d0 /)

    tetrag(1,:,5) = (/ 0d0, 0d0, 0d0 /)
    tetrag(2,:,5) = (/ 0d0, 1d0, 1d0 /)
    tetrag(3,:,5) = (/ 1d0, 0d0, 1d0 /)
    tetrag(4,:,5) = (/ 1d0, 1d0, 0d0 /)

    do i = 1, ntetrag
       ! calculate volume
       xp1 = c%x2c(tetrag(2,:,i) - tetrag(1,:,i))
       xp2 = c%x2c(tetrag(3,:,i) - tetrag(1,:,i))
       xp3 = c%x2c(tetrag(4,:,i) - tetrag(1,:,i))
       tvol(i) = abs(mixed(xp1,xp2,xp3)) / 6d0
       sumi = sumi + tvol(i)
    end do

  end subroutine pmwigner

  !xx! Private for wigner

  !> Private for wigner. Determines if two tetrahedra are equivalent.
  function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)

    logical :: equiv_tetrah
    type(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    real*8, dimension(0:3,3), intent(in) :: t1, t2
    integer, intent(in) :: leqv
    real*8, intent(in) :: lrotm(3,3,48), eps

    integer :: i, k, p
    real*8 :: r1(0:3,3), d1(0:3,3), dist2(0:3), eps2

    eps2 = eps*eps
    equiv_tetrah = .false.

    do i = 1, leqv
       do k = 0, 3
          r1(k,:) = matmul(lrotm(:,1:3,i),t1(k,:)-x0) + x0
       end do

       do p = 1, 6
          d1 = perm3(p,r1,t2)
          do k = 0, 3
             d1(k,:) = c%x2c(d1(k,:))
             dist2(k) = dot_product(d1(k,:),d1(k,:))
          end do
          if (all(dist2 < eps2)) then
             equiv_tetrah = .true.
             return
          end if
       end do
    end do

  end function equiv_tetrah

  !> Private for equiv_tetrah, wigner. 3! permutations.
  function perm3(p,r,t)

    real*8 :: perm3(0:3,3)
    integer, intent(in) :: p
    real*8, intent(in) :: r(0:3,3), t(0:3,3)

    select case(p)
    case(1)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(1,:)
       perm3(2,:) = r(2,:) - t(2,:)
       perm3(3,:) = r(3,:) - t(3,:)
    case(2)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(1,:)
       perm3(2,:) = r(2,:) - t(3,:)
       perm3(3,:) = r(3,:) - t(2,:)
    case(3)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(2,:)
       perm3(2,:) = r(2,:) - t(1,:)
       perm3(3,:) = r(3,:) - t(3,:)
    case(4)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(2,:)
       perm3(2,:) = r(2,:) - t(3,:)
       perm3(3,:) = r(3,:) - t(1,:)
    case(5)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(3,:)
       perm3(2,:) = r(2,:) - t(1,:)
       perm3(3,:) = r(3,:) - t(2,:)
    case(6)
       perm3(0,:) = r(0,:) - t(0,:)
       perm3(1,:) = r(1,:) - t(3,:)
       perm3(2,:) = r(2,:) - t(2,:)
       perm3(3,:) = r(3,:) - t(1,:)
    end select

  end function perm3

  !> Calculate the point group of a lattice. rmat is the matrix of
  !> lattice vectors (G = rmat * rmat'). verbose activates output to
  !> uout. ncen and xcen are the centering vectors for the lattice in
  !> crystallographic coordinates. nn and rot, if present, receive the
  !> symmetry operations for the lattice.
  subroutine lattpg(rmat,ncen,xcen,nn,rot)
    use sympg, only: nopsym, opsym, sym3d
    use tools_math, only: matinv
    use types, only: realloc
    real*8, intent(in) :: rmat(3,3)
    integer, intent(in) :: ncen
    real*8, intent(in) :: xcen(3,ncen)
    integer, intent(out), optional :: nn
    real*8, intent(out), optional :: rot(3,3,48)

    real*8 :: rmati(3,3), aal(3), gmat(3,3)
    integer :: i, na, nb, nc, npos, ia, ib, ic, it, op
    real*8  :: amax, amax2e, x(3), t(3), d2
    real*8, allocatable :: ax(:,:)
    integer, allocatable :: atZmol(:)

    ! the reciprocal-space matrix is the transpose of the inverse
    rmati = matinv(rmat)
    gmat = matmul(rmat,transpose(rmat))
    do i = 1, 3
       aal(i) = sqrt(gmat(i,i))
    end do

    ! every lattice vector must be represented in the molecule
    amax = 2d0*maxval(aal)
    amax2e = amax * amax + 1d-5
    na = int((amax/aal(1)) - 1d-7) + 2
    nb = int((amax/aal(2)) - 1d-7) + 2
    nc = int((amax/aal(3)) - 1d-7) + 2

    ! build the lattice point molecule
    npos = 0
    allocate(ax(3,10))
    do ia = -na, na
       do ib = -nb, nb
          do ic = -nc, nc
             do it = 1, ncen
                x = real((/ia,ib,ic/),8) + xcen(:,it)
                t = matmul(x,rmat)
                d2 = dot_product(t,t)
                if (d2 <= amax2e) then
                   npos = npos + 1
                   if (npos > size(ax,2)) call realloc(ax,3,2*npos)
                   ax(:,npos) = t
                end if
             end do
          end do
       end do
    end do
    call realloc(ax,3,npos)
    allocate(atZmol(npos))
    atzmol = 1

    ! Use symmetry to determine the point group. Default options
    call sym3d(rmat,ax(1,:),ax(2,:),ax(3,:),atZmol,npos,.false.)

    ! fill the reciprocal space matrices and clean up
    if (present(nn) .and. present(rot)) then
       nn = nopsym
       do op = 1, nn
          rot(:,:,op) = matmul(transpose(rmati),matmul(opsym(op,:,:),transpose(rmat)))
       end do
    end if
    deallocate(atZmol,ax)

  end subroutine lattpg

  !> Calculate the number of lattice vectors in each direction in
  !> order to be sure that the main cell is surrounded by a shell at
  !> least rmax thick. x2r is the cryst-to-car matrix for the
  !> lattice. The heuristic method has been adapted from gulp.
  subroutine search_lattice(x2r,rmax,imax,jmax,kmax)
    use param, only: pi

    real*8, intent(in) :: x2r(3,3), rmax
    integer, intent(out) :: imax, jmax, kmax

    integer :: i, nadd
    real*8 :: acell(3), bcell(3), gmat(3,3)

    gmat = matmul(transpose(x2r),x2r)
    do i = 1, 3
       acell(i) = sqrt(gmat(i,i))
    end do
    bcell(1) = acos(gmat(2,3) / acell(2) / acell(3)) * 180d0 / pi
    bcell(2) = acos(gmat(1,3) / acell(1) / acell(3)) * 180d0 / pi
    bcell(3) = acos(gmat(1,2) / acell(1) / acell(2)) * 180d0 / pi

    ! determine cell list, trick adapted from gulp
    nadd = 0
    if (bcell(1)<30 .or. bcell(2)<30 .or. bcell(3)<30 .or. bcell(1)>150 .or. bcell(2)>150 .or. bcell(3)>150) then
       nadd = 3
    else if (bcell(1)<50 .or. bcell(2)<50 .or. bcell(3)<50 .or. bcell(1)>130 .or. bcell(2)>130 .or. bcell(3)>130) then
       nadd = 2
    else if (bcell(1)<70 .or. bcell(2)<70 .or. bcell(3)<70 .or. bcell(1)>110 .or. bcell(2)>110 .or. bcell(3)>110) then
       nadd = 1
    else
       nadd = 1
    end if
    imax = ceiling(rmax / acell(1)) + nadd
    jmax = ceiling(rmax / acell(2)) + nadd
    kmax = ceiling(rmax / acell(3)) + nadd

  end subroutine search_lattice

  !> For the symmetry operation rot, compute the order and
  !> characteristic direction of the operation. The output is given in
  !> type (the type of operation: ident, c2, c3,...)  and vec (the
  !> direction of the axis or the normal to the plane, 0 for a point
  !> symmetry element). Used in the structure initialization.
  subroutine typeop(rot,type,vec,order)
    use tools_math, only: norm, eigns
    use tools_io, only: ferror, faterr
    use param, only: tpi, eye

    real*8, intent(in) :: rot(3,4) !< rotm operation
    integer, intent(out) :: type !< output type
    real*8, dimension(3), intent(out) :: vec !< output eigenvector
    integer, intent(out) :: order

    real*8, parameter :: eps = 1e-3

    integer, dimension(2,2) :: iord
    real*8 :: trace, tone
    real*8, dimension(3) :: eigen, eigeni
    real*8, dimension(3,3) :: mat, mat3
    integer :: nm, nones, nminusones
    integer :: i

    trace = rot(1,1) + rot(2,2) + rot(3,3)
    eigen = 0d0
    eigeni = 0d0
    mat = rot(1:3,1:3)
    iord = 0
    vec = 0d0

    !.Determine type of operation
    order = 1
    if (abs(trace+3) .lt. eps) then
       ! inversion
       order = 2
       type = inv
       return
    elseif (abs(trace-3) .lt.eps) then
       ! identity
       type = ident
       return
    else
       nm=3
       ! It is an axis or a plane. Diagonalize.
       call eigns(mat,eigen,eigeni)

       ! Determine actual operation.
       nones = 0
       nminusones = 0
       do i = 1, 3
          if (abs(eigen(i)-1) < eps .and. abs(eigeni(i)) < eps) then
             nones = nones + 1
             iord(1,nones) = i
          elseif (abs(eigen(i)+1) < eps .and. abs(eigeni(i)) < eps) then
             nminusones = nminusones +1
             iord(2,nminusones) = i
          endif
       enddo

       ! What is it?
       if (nones .eq. 1) then
          !.A proper axis. Order?
          tone = min(max((trace-1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(1,1))
          enddo
          if (order .eq. 2) then
             type = c2
          elseif (order .eq. 3) then
             type = c3
          elseif (order .eq. 4) then
             type = c4
          elseif (order .eq. 6) then
             type = c6
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       elseif (nones .eq. 2) then
          !A plane, Find directions.
          order = 2
          type = sigma
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
       elseif (nminusones .eq. 1) then
          !.An improper axes. Order?
          tone = min(max((trace+1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
          if (order .eq. 3) then
             type = s3
          elseif (order .eq. 4) then
             type = s4
          elseif (order .eq. 6) then
             mat3 = rot(1:3,1:3)
             mat3 = matmul(matmul(mat3,mat3),mat3)
             if (all(abs(mat3+eye) < 1d-5)) then
                type = s3
             else
                type = s6
             end if
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       else
          call ferror ('typeop', 'Sym. Element unknown', faterr)
       endif
    endif
    vec = vec / norm(vec)

  end subroutine typeop

  !> Get the holohedry and the Laue class from the Hermann-Mauguin
  !> point group label. Adapted from spglib, takes spglib HM point
  !> group labels.
  subroutine pointgroup_info(hmpg,schpg,holo,laue)
    use tools_io, only: equal
    character*(*), intent(in) :: hmpg
    character(len=3), intent(out) :: schpg
    integer, intent(out) :: holo
    integer, intent(out) :: laue

    if (equal(hmpg,"")) then
       schpg = ""
       holo = holo_unk
       laue = laue_unk
    elseif (equal(hmpg,"1")) then
       schpg = "C1"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"-1")) then
       schpg = "Ci"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"2")) then
       schpg = "C2"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"m")) then
       schpg = "Cs"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"2/m")) then
       schpg = "C2h"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"222")) then
       schpg = "D2"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mm2")) then
       schpg = "C2v"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mmm")) then
       schpg = "D2h"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"4")) then
       schpg = "C4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"-4")) then
       schpg = "S4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"4/m")) then
       schpg = "C4h"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"422")) then
       schpg = "D4"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4mm")) then
       schpg = "C4v"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"-42m")) then
       schpg = "D2d"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4/mmm")) then
       schpg = "D4h"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"3")) then
       schpg = "C3"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"-3")) then
       schpg = "C3i"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"32")) then
       schpg = "D3"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"3m")) then
       schpg = "C3v"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"-3m")) then
       schpg = "D3d"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"6")) then
       schpg = "C6"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"-6")) then
       schpg = "C3h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"6/m")) then
       schpg = "C6h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"622")) then
       schpg = "D6"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6mm")) then
       schpg = "C6v"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"-6m2")) then
       schpg = "D3h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6/mmm")) then
       schpg = "D6h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"23")) then
       schpg = "T"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"m-3")) then
       schpg = "Th"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"432")) then
       schpg = "O"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"-43m")) then
       schpg = "Td"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"m-3m")) then
       schpg = "Oh"
       holo = holo_cub
       laue = laue_m3m
    end if

  end subroutine pointgroup_info

  !> Wrapper to the spgs module. Sets the symetry in a crystal seed. 
  !> (including seed%havesym but not seed%findsym).
  subroutine spgs_wrap(seed,spg,usespgr)
    use spgs, only: spgs_ncv, spgs_cen, spgs_n, spgs_m, spgs_driver
    type(crystalseed), intent(inout) :: seed
    character*(*), intent(in) :: spg
    logical, intent(in) :: usespgr

    call spgs_driver(spg,usespgr)
    seed%ncv = spgs_ncv
    if (allocated(seed%cen)) deallocate(seed%cen)
    allocate(seed%cen(3,seed%ncv))
    seed%cen(:,1:seed%ncv) = real(spgs_cen(:,1:seed%ncv),8) / 12d0
    seed%neqv = spgs_n
    if (allocated(seed%rotm)) deallocate(seed%rotm)
    allocate(seed%rotm(3,4,spgs_n))
    seed%rotm = real(spgs_m(:,:,1:spgs_n),8)
    seed%rotm(:,4,:) = seed%rotm(:,4,:) / 12d0
    seed%havesym = 1

  end subroutine spgs_wrap

end module struct_basic
