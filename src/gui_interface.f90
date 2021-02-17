! Copyright (c) 2017 Alberto Otero de la Roza
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

!> Interface for the critic2 GUI.
module gui_interface
  use systemmod, only: system
  use crystalseedmod, only: crystalseed
  use iso_c_binding, only: c_ptr, c_null_ptr, c_float, c_char, c_int,&
     c_bool, c_null_char
  implicit none

  private

  !xx! interoperable types

  ! C-interoperable atom type
  type, bind(c) :: c_atom
     real(c_float) :: x(3) !< atom position (crystallographic)
     real(c_float) :: r(3) !< atom position (Cartesian, bohr)
     integer(c_int) :: is !< atom species
     integer(c_int) :: z !< atomic number
     character(kind=c_char,len=1) :: zsymb(3) !< atomic symbol
     character(kind=c_char,len=1) :: name(11) !< atomic name
     integer(c_int) :: idx !< index from the nneq list
     integer(c_int) :: cidx !< index from the complete list
     integer(c_int) :: flvec(3) !< lvec to the position in the fragment
     integer(c_int) :: ifrag !< which fragment this atom belongs to
     real(c_float) :: rad !< ball radius (bohr)
     real(c_float) :: rgb(4) !< color (0 to 1)
     integer(c_int) :: ncon !< number of neighbors
  end type c_atom

  ! scene type - holds all the information to render one scene.  This
  ! type is not c-interoperable. Access to individual scenes is
  ! achieved by remapping the pointers below.
  type scene
     integer :: idfile !< id of the file that generated this scene
     character(kind=c_char,len=1) :: file(512) !< name of the file
     character(kind=c_char,len=1) :: name(512) !< name of the scene
     integer :: isinit = 0 ! 0 = not init; 1 = seed; 2 = full
     logical :: readasfield = .false. ! if true, read a field from this file when initializing

     type(crystalseed) :: seed ! crystal seed for this scene
     type(system) :: sy ! system for this scene
     real*8 :: center(3) ! center of the scene (bohr)
     real(c_float) :: srad ! radius of the encompassing sphere

     logical(c_bool) :: ismolecule ! is this a molecule?

     integer(c_int) :: nf ! number of fields
     character(kind=c_char,len=1), allocatable :: fieldname(:,:) !< name of the fields

     integer(c_int) :: nat ! number of atoms
     type(c_atom), allocatable :: at(:) ! atoms

     integer(c_int), allocatable :: idcon(:,:) !< id (cidx) of the connected atom
     integer(c_int), allocatable :: lcon(:,:,:) !< lattice vector of the connected atom

     integer(c_int) :: nmol ! number of fragments
     integer(c_int), allocatable :: moldiscrete(:) ! is fragment discrete?

     real(c_float) :: avec(3,3) ! lattice vectors
     real(c_float) :: molx0(3) ! molecule centering translation
     real(c_float) :: molborder(3) ! molecular cell
  end type scene

  !xx! public interface

  ! Routines. All public routines are accessed from the GUI via critic2.h
  public :: gui_initialize
  public :: open_file
  public :: scene_initialize
  public :: set_scene_pointers
  public :: scene_set_reference_field
  public :: gui_end
  private :: realloc_scene

  ! Variables. All bind(c) variables are accessed from the GUI via critic2.h.

  ! General information for the run
  character(kind=c_char,len=1), allocatable, target :: c2home_c(:) ! Location of data files
  type(c_ptr), bind(c) :: c2home

  ! File and scene tree
  integer(c_int), bind(c) :: nfiles = 0 ! Number of files
  integer(c_int), bind(c) :: nsc = 0 ! Number of scenes
  type(scene), allocatable, target :: sc(:) ! Information about the loaded scenes
  integer :: ilastfile = 0 ! Integer pointer to the last file

  ! Error message container to pass errors to the GUI
  character(kind=c_char,len=1), target :: errmsg_c(512)
  type(c_ptr), bind(c) :: errmsg

  ! Pointers to the current scene. The correspond mostly to pointers
  ! to (or copies of) the information in an instance of the scene type
  ! above. Those variables that do not appear in the scene type have a
  ! comment after them.
  integer(c_int) :: icursc = -1 ! Id of the current scene
  logical :: scupdated = .false. ! .true. if the current scene is up to date

  integer(c_int), bind(c) :: idfile
  type(c_ptr), bind(c) :: file
  type(c_ptr), bind(c) :: name
  integer(c_int), bind(c) :: isinit

  real(c_float), bind(c) :: scenerad ! (srad, radius of the encompassing sphere)

  integer(c_int), bind(c) :: ismolecule

  integer(c_int), bind(c) :: nf
  integer(c_int), bind(c) :: iref ! Current reference field for this system
  type(c_ptr), bind(c) :: fieldname

  integer(c_int), bind(c) :: nat
  type(c_ptr), bind(c) :: at

  integer(c_int), bind(c) :: mncon ! Maximum number of atom neighbors
  type(c_ptr), bind(c) :: idcon
  type(c_ptr), bind(c) :: lcon

  integer(c_int), bind(c) :: nmol
  type(c_ptr), bind(c) :: moldiscrete

  type(c_ptr), bind(c) :: avec(3)
  type(c_ptr), bind(c) :: molx0
  type(c_ptr), bind(c) :: molborder

  !xx! private interface
  ! parameters
  real(c_float), parameter, private :: minsrad = 10.0_c_float ! minimum scene radius

contains

  !> Initialize the critic2 GUI.
  subroutine gui_initialize() bind(c)
    use c_interface_module, only: f_c_string
    use iso_fortran_env, only: input_unit, output_unit
    use iso_c_binding, only: c_loc
    use systemmod, only: systemmod_init
    use config, only: getstring, istring_datadir
    use global, only: global_init, config_write, initial_banner, critic_home
    use tools_io, only: ucopy, uout, start_clock, lualloc_init, &
       tictac, interactive, uin, filepath
    use param, only: param_init

    ! initialize parameters
    call start_clock()
    call param_init()

    ! input/output, arguments (tools_io)
    call lualloc_init()
    uin = input_unit
    uout = output_unit
    interactive = .false.
    filepath = "."

    ! set default values and initialize the rest of the modules
    call global_init("",getstring(istring_datadir))
    call systemmod_init(1)

    ! banner and compilation info; do not copy input
    call initial_banner()
    call config_write()
    call tictac('CRITIC2')
    write (uout,*)
    ucopy = -1

    ! allocate the initial scene
    if (allocated(sc)) deallocate(sc)
    allocate(sc(1))
    nsc = 0
    isinit = 0
    idfile = 0
    scenerad = minsrad
    scupdated = .true.

    ! set the global information variables
    if (allocated(c2home_c)) deallocate(c2home_c)
    allocate(c2home_c(len(critic_home)+1))
    call f_c_string(critic_home,c2home_c)
    c2home = c_loc(c2home_c)

  end subroutine gui_initialize

  !> Open one or more scenes from all files in the line. ismolecule: 0
  !> = crystal, 1 = molecule, -1 = critic2 decides.
  function open_file(file0,ismolecule) bind(c)
    use iso_c_binding, only: c_loc
    use c_interface_module, only: c_string_value, f_c_string
    use crystalseedmod, only: read_seeds_from_file, crystalseed
    type(c_ptr), intent(in) :: file0
    integer(c_int), value :: ismolecule
    integer(c_int) :: open_file

    character(len=:), allocatable :: file
    integer :: iseed, nseed, iafield
    type(crystalseed), allocatable :: seed(:)
    character(len=:), allocatable :: errmsg_

    ! transform to fortran string
    file = c_string_value(file0)

    ! read all seeds from the line
    call read_seeds_from_file(file,ismolecule,nseed,seed,errmsg_,iafield)

    if (nseed > 0) then
       nfiles = nfiles + 1
       ilastfile = ilastfile + 1

       if (nsc + nseed > size(sc,1)) &
          call realloc_scene(sc,nsc+nseed)

       if (iafield > 0) &
          sc(nsc+iafield)%readasfield  = .true.

       do iseed = 1, nseed
          ! initialize the system from the first seed
          nsc = nsc + 1
          sc(nsc)%idfile = ilastfile
          call f_c_string(trim(seed(iseed)%file),sc(nsc)%file)
          call f_c_string(trim(seed(iseed)%name),sc(nsc)%name)
          sc(nsc)%seed = seed(iseed)
          sc(nsc)%isinit = 1
          if (nsc == icursc) scupdated = .true.
       end do
       open_file = 1
    else
       open_file = 0
       call f_c_string(errmsg_,errmsg_c,512)
       errmsg = c_loc(errmsg_c)
    end if

  end function open_file

  subroutine scene_initialize(isc) bind(c)
    use tools_io, only: nameguess
    use c_interface_module, only: f_c_string
    use param, only: atmcov, jmlcol
    integer(c_int), value, intent(in) :: isc

    integer :: i, j, idx, iz, is, i1, i2, i3
    real(c_float) :: xmin(3), xmax(3)
    integer :: mncon_, id
    character(len=:), allocatable :: errmsg
    real*8 :: x(3)
    character*2 :: zsymb

    if (isc > nsc .or. isc < 1) return
    if (sc(isc)%isinit == 0 .or. sc(isc)%isinit == 2) return

    ! initialize the system
    sc(isc)%isinit = 2
    call sc(isc)%sy%new_from_seed(sc(isc)%seed)

    ! load any field
    if (sc(isc)%readasfield) then
       call sc(isc)%sy%load_field_string(sc(isc)%seed%file,id,errmsg)
       sc(isc)%sy%f(id)%file = sc(isc)%seed%file
       sc(isc)%sy%f(id)%name = sc(isc)%seed%file
    end if

    ! report
    call sc(isc)%sy%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
    sc(isc)%center = 0d0

    ! build the atom list
    mncon_ = 0
    sc(isc)%nat = sc(isc)%sy%c%ncel
    if (allocated(sc(isc)%at)) deallocate(sc(isc)%at)
    allocate(sc(isc)%at(sc(isc)%nat))
    do i = 1, sc(isc)%nat
       is = sc(isc)%sy%c%atcel(i)%is
       idx = sc(isc)%sy%c%atcel(i)%idx
       iz = sc(isc)%sy%c%spc(is)%z

       sc(isc)%at(i)%x = real(sc(isc)%sy%c%atcel(i)%x,c_float)
       sc(isc)%at(i)%r = real(sc(isc)%sy%c%atcel(i)%r,c_float)
       sc(isc)%at(i)%is = is
       sc(isc)%at(i)%idx = idx
       sc(isc)%at(i)%cidx = i
       sc(isc)%at(i)%z = iz
       zsymb = nameguess(iz,.true.)
       call f_c_string(sc(isc)%sy%c%spc(is)%name,sc(isc)%at(i)%name,11)
       call f_c_string(zsymb,sc(isc)%at(i)%zsymb,3)
       if (allocated(sc(isc)%sy%c%nstar)) then
          sc(isc)%at(i)%ncon = sc(isc)%sy%c%nstar(i)%ncon
       else
          sc(isc)%at(i)%ncon = 0
       end if
       mncon_ = max(mncon_,sc(isc)%at(i)%ncon)

       if (iz > 0) then
          sc(isc)%at(i)%rad = real(0.7d0*atmcov(iz),c_float)
          sc(isc)%at(i)%rgb(1:3) = real(jmlcol(:,iz),4) / 255.
          sc(isc)%at(i)%rgb(4) = 1.0
       else
          sc(isc)%at(i)%rad = real(0.4d0,c_float)
          sc(isc)%at(i)%rgb(1:3) = (/1d0,0.4314d0,0.7059d0/) ! hotpink1
          sc(isc)%at(i)%rgb(4) = 1.0
       end if
    end do

    ! build the fields info
    sc(isc)%nf = sc(isc)%sy%nf

    if (allocated(sc(isc)%fieldname)) deallocate(sc(isc)%fieldname)
    allocate(sc(isc)%fieldname(255,0:sc(isc)%nf))
    do i = 0, sc(isc)%nf
       call f_c_string(sc(isc)%sy%f(i)%name,sc(isc)%fieldname(:,i),255)
    end do

    ! build the fragment info
    sc(isc)%nmol = sc(isc)%sy%c%nmol
    if (allocated(sc(isc)%moldiscrete)) deallocate(sc(isc)%moldiscrete)
    allocate(sc(isc)%moldiscrete(sc(isc)%nmol))
    do i = 1, sc(isc)%sy%c%nmol
       if (sc(isc)%sy%c%mol(i)%discrete) then
          sc(isc)%moldiscrete(i) = 1
       else
          sc(isc)%moldiscrete(i) = 0
       end if
       do j = 1, sc(isc)%sy%c%mol(i)%nat
          idx = sc(isc)%sy%c%mol(i)%at(j)%cidx
          sc(isc)%at(idx)%flvec = sc(isc)%sy%c%mol(i)%at(j)%lvec
          sc(isc)%at(idx)%ifrag = i-1
       end do
    end do

    ! build the neighbor info
    if (allocated(sc(isc)%idcon)) deallocate(sc(isc)%idcon)
    if (allocated(sc(isc)%lcon)) deallocate(sc(isc)%lcon)
    allocate(sc(isc)%idcon(mncon_,sc(isc)%nat))
    allocate(sc(isc)%lcon(3,mncon_,sc(isc)%nat))
    sc(isc)%idcon = 0
    sc(isc)%lcon = 0
    do i = 1, sc(isc)%nat
       do j = 1, sc(isc)%at(i)%ncon
          sc(isc)%idcon(j,i) = sc(isc)%sy%c%nstar(i)%idcon(j)-1
          sc(isc)%lcon(:,j,i) = sc(isc)%sy%c%nstar(i)%lcon(:,j)
       end do
    end do

    ! lattice vectors
    sc(isc)%avec = real(sc(isc)%sy%c%m_x2c,c_float)
    sc(isc)%ismolecule = sc(isc)%sy%c%ismolecule
    sc(isc)%molx0 = real(sc(isc)%sy%c%molx0,c_float)
    sc(isc)%molborder = real(sc(isc)%sy%c%molborder,c_float)

    ! calculate the scene radius
    if (sc(isc)%nat == 0) then
       xmin = 0._c_float
       xmax = 0._c_float
    else
       xmin = sc(isc)%at(1)%r
       xmax = sc(isc)%at(1)%r
       do i = 2, sc(isc)%nat
          xmax = max(sc(isc)%at(i)%r,xmax)
          xmin = min(sc(isc)%at(i)%r,xmin)
       end do
       if (.not.sc(isc)%ismolecule) then
          do i1 = 0, 1
             do i2 = 0, 1
                do i3 = 0, 1
                   x = real((/i1,i2,i3/),8)
                   x = sc(isc)%sy%c%x2c(x)
                   xmax = max(real(x,c_float),xmax)
                   xmin = min(real(x,c_float),xmin)
                end do
             end do
          end do
       end if
    end if
    sc(isc)%srad = max(norm2(xmax-xmin),minsrad)

    ! this scene has been updated
    if (isc == icursc) scupdated = .true.

  end subroutine scene_initialize

  subroutine set_scene_pointers(isc) bind(c)
    use iso_c_binding, only: c_loc
    integer(c_int), value, intent(in) :: isc
    if (isc == icursc .and..not.scupdated) return

    nat = 0
    isinit = 0
    if (isc < 1 .or. isc > nsc) return

    isinit = sc(isc)%isinit
    idfile = sc(isc)%idfile
    file = c_loc(sc(isc)%file)
    name = c_loc(sc(isc)%name)

    scenerad = sc(isc)%srad

    nf = sc(isc)%nf
    iref = sc(isc)%sy%iref
    fieldname = c_loc(sc(isc)%fieldname)

    nat = sc(isc)%nat
    at = c_loc(sc(isc)%at)

    mncon = size(sc(isc)%idcon,1)
    idcon = c_loc(sc(isc)%idcon)
    lcon = c_loc(sc(isc)%lcon)

    nmol = sc(isc)%nmol
    moldiscrete = c_loc(sc(isc)%moldiscrete)

    avec(1) = c_loc(sc(isc)%avec(1,1))
    avec(2) = c_loc(sc(isc)%avec(1,2))
    avec(3) = c_loc(sc(isc)%avec(1,3))
    if (sc(isc)%ismolecule) then
       ismolecule = 1
    else
       ismolecule = 0
    end if
    molx0 = c_loc(sc(isc)%molx0)
    molborder = c_loc(sc(isc)%molborder)
    icursc = isc
    scupdated = .false.

  end subroutine set_scene_pointers

  !> Change the reference field for the system in scene isc.
  subroutine scene_set_reference_field(isc, iref) bind(c)
    integer(c_int), intent(in), value :: isc
    integer(c_int), intent(in), value :: iref

    if (isc < 1 .or. isc > nsc) return
    if (.not.sc(isc)%sy%goodfield(iref)) return
    call sc(isc)%sy%set_reference(iref,.false.)
    call sc(isc)%sy%report(.false.,.false.,.true.,.false.,.false.,.false.,.true.)
    if (isc == icursc) scupdated = .true.

  end subroutine scene_set_reference_field

  !> Clean up the GUI data structures.
  subroutine gui_end() bind(c)
    use grid1mod, only: grid1_clean_grids
    use tools_io, only: print_clock, tictac, ncomms, nwarns, uout, string

    ! deallocate scene
    if (allocated(sc)) deallocate(sc)
    nsc = 0
    nfiles = 0

    ! kill atomic grids
    call grid1_clean_grids()

    ! final message
    write (uout,'("CRITIC2 ended successfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
       string(nwarns), string(ncomms)
    call print_clock()
    call tictac('CRITIC2')

  end subroutine gui_end

  !> Adapt the size of an allocatable 1D type(scene) array
  subroutine realloc_scene(a,nnew)
    type(scene), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(scene), allocatable :: temp(:)
    integer :: l1, u1

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_scene

end module gui_interface
