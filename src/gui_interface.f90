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
  use iso_c_binding, only: c_ptr, c_null_ptr, c_float, c_char, c_int
  implicit none

  private

  !xx! interoperable types
  ! C-interoperable atom type
  type, bind(c) :: c_atom
     character(kind=c_char,len=1) :: label(11) !< Label
     integer(c_int) :: z !< atomic number
     real(c_float) :: r(3) !< center position (bohr) 
     real(c_float) :: rad !< ball radius (bohr) 
     real(c_float) :: rgb(4) !< color (0 to 1)
  end type c_atom

  !xx! private to this module
  type(c_ptr) :: window = c_null_ptr ! window pointer

  type scene
     integer :: isinit = 0 ! 0 = not init; 1 = seed; 2 = full
     type(crystalseed) :: seed ! crystal seed for this scene
     type(system) :: sy ! system for this scene
     real*8 :: center(3) ! center of the scene (bohr)
     integer(c_int) :: nat ! number of atoms
     type(c_atom), allocatable :: at(:) ! atoms
     real(c_float) :: srad ! radius of the encompassing sphere
  end type scene
  integer :: nsc = 0
  type(scene), allocatable, target :: sc(:)

  !xx! public interface
  ! routines
  public :: gui_initialize
  public :: open_file
  public :: set_scene_pointers
  public :: gui_end

  ! pointers to the current scene
  integer(c_int), bind(c) :: nat
  type(c_ptr), bind(c) :: at
  real(c_float), bind(c) :: scenerad

contains

  !> Initialize the critic2 GUI.
  subroutine gui_initialize(window_) bind(c)
    use iso_fortran_env, only: input_unit, output_unit
    use systemmod, only: systemmod_init
    use spgs, only: spgs_init
    use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
       fcflags, cc, cflags, ldflags, enable_debug, package
    use global, only: global_init, config_write, initial_banner
    use tools_io, only: ioinit, ucopy, uout, start_clock, &
       tictac, interactive, uin, filepath
    use param, only: param_init

    type(c_ptr), intent(in), value :: window_

    ! save the window pointer
    window = window_
    call assert_non_null(window,"gui_initialize","window")

    ! initialize parameters
    call start_clock()
    call param_init()

    ! input/output, arguments (tools_io)
    call ioinit()
    uin = input_unit
    uout = output_unit
    interactive = .false.
    filepath = "."

    ! set default values and initialize the rest of the modules
    call global_init("",datadir)
    call spgs_init()
    call systemmod_init(1)

    ! banner and compilation info; do not copy input
    call initial_banner()
    call config_write(package,version,atarget,adate,f77,fflags,fc,&
       fcflags,cc,cflags,ldflags,enable_debug,datadir)
    call tictac('CRITIC2')
    write (uout,*)
    ucopy = -1

    ! allocate the initial scene
    if (allocated(sc)) deallocate(sc)
    allocate(sc(1))
    nsc = 0

  end subroutine gui_initialize

  !> Open one or more scenes from all files in the line. ismolecule: 0
  !> = crystal, 1 = molecule, -1 = critic2 decides.
  subroutine open_file(line0,ismolecule) bind(c)
    use c_interface_module, only: c_string_value, f_c_string
    use iso_c_binding, only: c_int
    use crystalseedmod, only: read_seeds_from_file, crystalseed
    use param, only: pi, atmcov, jmlcol
    type(c_ptr), intent(in) :: line0
    integer(c_int), value :: ismolecule

    integer :: lp
    character(len=:), allocatable :: line
    integer :: i, idx, iz, nseed
    type(crystalseed), allocatable :: seed(:)
    real(c_float) :: xmax(3)

    ! transform to fortran string
    line = c_string_value(line0)
    
    ! read all seeds from the line
    lp = 1
    call read_seeds_from_file(line,lp,ismolecule,nseed,seed)
    
    if (nseed > 0) then
       ! initialize the system from the first seed
       nsc = 1
       sc(1)%seed = seed(1)
       sc(1)%isinit = 2
       call sc(1)%sy%new_from_seed(sc(1)%seed)
       call sc(1)%sy%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
       sc(1)%center = 0d0

       ! build the atoms
       sc(1)%nat = sc(1)%sy%c%ncel
       if (allocated(sc(1)%at)) deallocate(sc(1)%at)
       allocate(sc(1)%at(sc(1)%nat))
       do i = 1, sc(1)%nat
          idx = sc(1)%sy%c%atcel(i)%idx
          iz = sc(1)%sy%c%at(idx)%z
          sc(1)%at(i)%r = sc(1)%sy%c%atcel(i)%r
          sc(1)%at(i)%z = iz
          call f_c_string(sc(1)%sy%c%at(idx)%name,sc(1)%at(i)%label,11)
          if (atmcov(iz) > 1) then
             sc(1)%at(i)%rad = atmcov(iz)
          else
             sc(1)%at(i)%rad = 2.*atmcov(iz)
          end if
          sc(1)%at(i)%rgb(1:3) = real(jmlcol(:,iz),4) / 255.
          sc(1)%at(i)%rgb(4) = 1.0
          sc(1)%center = sc(1)%center + sc(1)%at(i)%r
       end do

       ! translate to the center of mass
       xmax = 0._c_float
       sc(1)%center = sc(1)%center / sc(1)%nat
       do i = 1, sc(1)%nat
          sc(1)%at(i)%r = sc(1)%at(i)%r - sc(1)%center
          xmax = max(abs(sc(1)%at(i)%r),xmax)
       end do
       sc(1)%srad = sqrt(dot_product(xmax,xmax))
    end if

  end subroutine open_file

  subroutine set_scene_pointers(isc) bind(c)
    use iso_c_binding, only: c_loc
    use gui_glfw
    use gui_glu
    use gui_gl
    integer(c_int), value, intent(in) :: isc

    nat = 0
    if (isc < 0 .or. isc > nsc) return

    nat = sc(isc)%nat
    at = c_loc(sc(isc)%at)
    scenerad = sc(isc)%srad

  end subroutine set_scene_pointers

  subroutine gui_end() bind(c)
    use grid1mod, only: grid1_clean_grids
    use tools_io, only: print_clock, tictac, ncomms, nwarns, uout, string

    ! deallocate scene
    if (allocated(sc)) deallocate(sc)
    nsc = 0

    ! kill atomic grids
    call grid1_clean_grids()
    
    ! final message
    write (uout,'("CRITIC2 ended succesfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
       string(nwarns), string(ncomms)
    call print_clock()
    call tictac('CRITIC2')
    
  end subroutine gui_end

  subroutine assert_non_null(ptr,routine,ptrname)
    use iso_c_binding, only: c_associated
    use tools_io, only: ferror, faterr
    character*(*), intent(in) :: routine, ptrname
    type(c_ptr), intent(in) :: ptr

    if (.not.c_associated(ptr)) &
       call ferror(trim(routine),"Unexpected NULL pointer - " // trim(ptrname),faterr)

  end subroutine assert_non_null

end module gui_interface
