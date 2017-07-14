! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Output in several common graphics formats
module graphics
  implicit none
  
  private

  type grhandle
     integer :: lu = 0 !< Logical unit for the graphics file
     integer :: lumtl = 0 !< Logical unit for the matierlas file (obj format)
     integer :: fmt = 0 !< File format
     character*(255) :: file = "" !< File name
     integer :: nball = 0 !< Number of balls
     integer :: nstick = 0 !< Number of sticks
     integer :: nsurf = 0 !< Number of surfaces
     integer :: nface = 0 !< Number of faces
     integer :: nmtl = 0 !< Number of materials
     integer :: nv = 0 !< Number of vertices
     integer :: nf = 0 !< Number of faces
     integer, allocatable :: mtlrgb(:,:) !< Material definitions
   contains
     procedure :: open => graphics_open
     procedure :: close => graphics_close
     procedure :: ball => graphics_ball
     procedure :: polygon => graphics_polygon
     procedure :: stick => graphics_stick
     procedure :: surf => graphics_surf
  end type grhandle
  public :: grhandle

  private :: graphics_init
  private :: graphics_open
  private :: graphics_close
  private :: graphics_ball
  private :: graphics_polygon
  private :: graphics_stick
  private :: graphics_surf
  private :: obj_open
  private :: obj_close
  private :: obj_ball
  private :: obj_polygon
  private :: obj_stick
  private :: obj_surf
  private :: register_texture
  private :: ply_open
  private :: ply_close
  private :: ply_ball
  private :: ply_polygon
  private :: ply_stick
  private :: ply_surf
  private :: off_open
  private :: off_close
  private :: off_ball
  private :: off_polygon
  private :: off_stick
  private :: off_surf

  ! graphics database
  logical :: isinit = .false.

  ! icosahedron models
  integer :: nvsph(0:2), nfsph(0:2)
  real*8 :: vsph(3,42,0:2)
  integer :: fsph(3,80,0:2)

  ! cylinder models
  integer :: nvcyl(0:1), nfcyl(0:1)
  real*8 :: vcyl(3,34,0:1)
  integer :: fcyl(3,64,0:1)

  ! formats
  integer, parameter :: ifmt_unk = 0
  integer, parameter :: ifmt_obj = 1
  integer, parameter :: ifmt_ply = 2
  integer, parameter :: ifmt_off = 3

contains

  ! initialize the icosahedron vertices and faces
  subroutine graphics_init()
    
    ! fill the icosahedron models
    vsph = 0d0
    fsph = 0

    ! lvl 0 -- a tetrahedron
    nvsph(0) = 4
    vsph(1:3,1:4,0) = reshape((/&
        0.57735d0,  0.57735d0,  0.57735d0,&
       -0.57735d0, -0.57735d0,  0.57735d0,&
       -0.57735d0,  0.57735d0, -0.57735d0,&
        0.57735d0, -0.57735d0, -0.57735d0/),(/3,4/))

    nfsph(0) = 4
    fsph(1:3,1:4,0) = reshape((/&
        1, 2, 3,&
        1, 2, 4,&
        1, 3, 4,&
        2, 3, 4/),(/3,4/))

    ! lvl 1 -- an icosahedron
    nvsph(1) = 12
    vsph(1:3,1:12,1) = reshape((/&
        0.85065d0,   0.00000d0,   0.52573d0,&
        0.85065d0,   0.00000d0,  -0.52573d0,&
       -0.85065d0,   0.00000d0,   0.52573d0,&
       -0.85065d0,   0.00000d0,  -0.52573d0,&
        0.00000d0,   0.52573d0,   0.85065d0,&
        0.00000d0,  -0.52573d0,   0.85065d0,&
        0.00000d0,   0.52573d0,  -0.85065d0,&
        0.00000d0,  -0.52573d0,  -0.85065d0,&
        0.52573d0,   0.85065d0,   0.00000d0,&
       -0.52573d0,   0.85065d0,   0.00000d0,&
        0.52573d0,  -0.85065d0,   0.00000d0,&
       -0.52573d0,  -0.85065d0,   0.00000d0/),(/3,12/))

    nfsph(1) = 20
    fsph(1:3,1:20,1) = reshape((/&
        1,  9,  5,&
       10,  9,  7,&
        3,  5, 10,&
        9, 10,  5,&
        3, 10,  4,&
       10,  7,  4,&
        1, 11,  2,&
        9,  1,  2,&
        7,  9,  2,&
       12,  3,  4,&
        6,  1,  5,&
       11,  1,  6,&
       12, 11,  6,&
        5,  3,  6,&
        6,  3, 12,&
        8,  7,  2,&
        8,  2, 11,&
        8, 11, 12,&
        4,  7,  8,&
       12,  4,  8/),(/3,20/))
    
    ! lvl 2 -- a subdivided icosahedron (aka Edward Carnby's head)
    nvsph(2) = 42
    vsph(1:3,1:42,2) = reshape((/&
        0.85065d0,  0.00000d0,  0.52573d0,&
        0.85065d0,  0.00000d0, -0.52573d0,&
       -0.85065d0,  0.00000d0,  0.52573d0,&
       -0.85065d0,  0.00000d0, -0.52573d0,&
        0.00000d0,  0.52573d0,  0.85065d0,&
        0.00000d0, -0.52573d0,  0.85065d0,&
        0.00000d0,  0.52573d0, -0.85065d0,&
        0.00000d0, -0.52573d0, -0.85065d0,&
        0.52573d0,  0.85065d0,  0.00000d0,&
       -0.52573d0,  0.85065d0,  0.00000d0,&
        0.52573d0, -0.85065d0,  0.00000d0,&
       -0.52573d0, -0.85065d0,  0.00000d0,&
        0.80902d0,  0.50000d0,  0.30902d0,&
        0.30902d0,  0.80902d0,  0.50000d0,&
        0.50000d0,  0.30902d0,  0.80902d0,&
        0.00000d0,  1.00000d0,  0.00000d0,&
        0.30902d0,  0.80902d0, -0.50000d0,&
       -0.30902d0,  0.80902d0, -0.50000d0,&
       -0.50000d0,  0.30902d0,  0.80902d0,&
       -0.30902d0,  0.80902d0,  0.50000d0,&
       -0.80902d0,  0.50000d0,  0.30902d0,&
       -0.80902d0,  0.50000d0, -0.30902d0,&
       -1.00000d0,  0.00000d0,  0.00000d0,&
       -0.50000d0,  0.30902d0, -0.80902d0,&
        0.80902d0, -0.50000d0,  0.30902d0,&
        0.80902d0, -0.50000d0, -0.30902d0,&
        1.00000d0,  0.00000d0,  0.00000d0,&
        0.80902d0,  0.50000d0, -0.30902d0,&
        0.50000d0,  0.30902d0, -0.80902d0,&
       -0.80902d0, -0.50000d0,  0.30902d0,&
       -0.80902d0, -0.50000d0, -0.30902d0,&
        0.50000d0, -0.30902d0,  0.80902d0,&
        0.00000d0,  0.00000d0,  1.00000d0,&
        0.30902d0, -0.80902d0,  0.50000d0,&
        0.00000d0, -1.00000d0,  0.00000d0,&
       -0.30902d0, -0.80902d0,  0.50000d0,&
       -0.50000d0, -0.30902d0,  0.80902d0,&
        0.00000d0,  0.00000d0, -1.00000d0,&
        0.50000d0, -0.30902d0, -0.80902d0,&
        0.30902d0, -0.80902d0, -0.50000d0,&
       -0.30902d0, -0.80902d0, -0.50000d0,&
       -0.50000d0, -0.30902d0, -0.80902d0/),(/3,42/))

    nfsph(2) = 80
    fsph(1:3,1:80,2) = reshape((/&
        1, 13, 15,&
       10, 16, 18,&
        3, 19, 21,&
        9, 16, 14,&
        3, 21, 23,&
       10, 18, 22,&
        1, 25, 27,&
        9, 13, 28,&
        7, 17, 29,&
       12, 30, 31,&
        6, 32, 33,&
       11, 25, 34,&
       12, 35, 36,&
        5, 19, 33,&
        6, 37, 36,&
        8, 38, 39,&
        8, 39, 40,&
        8, 40, 41,&
        4, 24, 42,&
       12, 31, 41,&
        9, 14, 13,&
        5, 15, 14,&
       13, 14, 15,&
        9, 17, 16,&
        7, 18, 17,&
       16, 17, 18,&
        5, 20, 19,&
       10, 21, 20,&
       19, 20, 21,&
       10, 20, 16,&
        5, 14, 20,&
       16, 20, 14,&
       10, 22, 21,&
        4, 23, 22,&
       21, 22, 23,&
        7, 24, 18,&
        4, 22, 24,&
       18, 24, 22,&
       11, 26, 25,&
        2, 27, 26,&
       25, 26, 27,&
        1, 27, 13,&
        2, 28, 27,&
       13, 27, 28,&
        9, 28, 17,&
        2, 29, 28,&
       17, 28, 29,&
        3, 23, 30,&
        4, 31, 23,&
       30, 23, 31,&
        1, 15, 32,&
        5, 33, 15,&
       32, 15, 33,&
        1, 32, 25,&
        6, 34, 32,&
       25, 32, 34,&
       11, 34, 35,&
        6, 36, 34,&
       35, 34, 36,&
        3, 37, 19,&
        6, 33, 37,&
       19, 37, 33,&
        3, 30, 37,&
       12, 36, 30,&
       37, 30, 36,&
        7, 29, 38,&
        2, 39, 29,&
       38, 29, 39,&
        2, 26, 39,&
       11, 40, 26,&
       39, 26, 40,&
       11, 35, 40,&
       12, 41, 35,&
       40, 35, 41,&
        7, 38, 24,&
        8, 42, 38,&
       24, 38, 42,&
        4, 42, 31,&
        8, 41, 42,&
       31, 42, 41/),(/3,80/))

    ! fill the cylinder models
    vcyl = 0d0
    fcyl = 0

    ! lvl 0 -- hexagonal
    nvcyl(0) = 14
    vcyl(1:3,1:14,0) = reshape((/&
       -0.86603d0,  0.50000d0,  0.00000d0,&
       -0.86603d0, -0.50000d0,  0.00000d0,&
       -0.00000d0, -1.00000d0,  0.00000d0,&
        0.86603d0, -0.50000d0,  0.00000d0,&
        0.86603d0,  0.50000d0,  0.00000d0,&
        0.00000d0,  1.00000d0,  0.00000d0,&
       -0.86603d0,  0.50000d0,  1.00000d0,&
       -0.86603d0, -0.50000d0,  1.00000d0,&
       -0.00000d0, -1.00000d0,  1.00000d0,&
        0.86603d0, -0.50000d0,  1.00000d0,&
        0.86603d0,  0.50000d0,  1.00000d0,&
        0.00000d0,  1.00000d0,  1.00000d0,&
        0.00000d0,  0.00000d0,  0.00000d0,&
        0.00000d0,  0.00000d0,  1.00000d0/),(/3,14/))
    
    nfcyl(0) = 24
    fcyl(1:3,1:24,0) = reshape((/&
       13,  1,  2,&
       14,  7,  8,&
        1,  2,  7,&
        8,  7,  2,&
       13,  2,  3,&
       14,  8,  9,&
        2,  3,  8,&
        9,  8,  3,&
       13,  3,  4,&
       14,  9, 10,&
        3,  4,  9,&
       10,  9,  4,&
       13,  4,  5,&
       14, 10, 11,&
        4,  5, 10,&
       11, 10,  5,&
       13,  5,  6,&
       14, 11, 12,&
        5,  6, 11,&
       12, 11,  6,&
       13,  6,  1,&
       14, 12,  7,&
        6,  1, 12,&
        7, 12,  1/),(/3,24/))
       
    ! lvl 1 -- hexadecagonal
    nvcyl(1) = 34
    vcyl(1:3,1:34,1) = reshape((/&
       -0.38268,  0.92388,  0.00000,&
       -0.70711,  0.70711,  0.00000,&
       -0.92388,  0.38268,  0.00000,&
       -1.00000, -0.00000,  0.00000,&
       -0.92388, -0.38268,  0.00000,&
       -0.70711, -0.70711,  0.00000,&
       -0.38268, -0.92388,  0.00000,&
        0.00000, -1.00000,  0.00000,&
        0.38268, -0.92388,  0.00000,&
        0.70711, -0.70711,  0.00000,&
        0.92388, -0.38268,  0.00000,&
        1.00000,  0.00000,  0.00000,&
        0.92388,  0.38268,  0.00000,&
        0.70711,  0.70711,  0.00000,&
        0.38268,  0.92388,  0.00000,&
       -0.00000,  1.00000,  0.00000,&
       -0.38268,  0.92388,  1.00000,&
       -0.70711,  0.70711,  1.00000,&
       -0.92388,  0.38268,  1.00000,&
       -1.00000, -0.00000,  1.00000,&
       -0.92388, -0.38268,  1.00000,&
       -0.70711, -0.70711,  1.00000,&
       -0.38268, -0.92388,  1.00000,&
        0.00000, -1.00000,  1.00000,&
        0.38268, -0.92388,  1.00000,&
        0.70711, -0.70711,  1.00000,&
        0.92388, -0.38268,  1.00000,&
        1.00000,  0.00000,  1.00000,&
        0.92388,  0.38268,  1.00000,&
        0.70711,  0.70711,  1.00000,&
        0.38268,  0.92388,  1.00000,&
       -0.00000,  1.00000,  1.00000,&
        0.00000,  0.00000,  0.00000,&
        0.00000,  0.00000,  1.00000/),(/3,34/))
       
    nfcyl(1) = 64
    fcyl(1:3,1:64,1) = reshape((/&
       33,  1,  2,&
       34, 17, 18,&
        1,  2, 17,&
       18, 17,  2,&
       33,  2,  3,&
       34, 18, 19,&
        2,  3, 18,&
       19, 18,  3,&
       33,  3,  4,&
       34, 19, 20,&
        3,  4, 19,&
       20, 19,  4,&
       33,  4,  5,&
       34, 20, 21,&
        4,  5, 20,&
       21, 20,  5,&
       33,  5,  6,&
       34, 21, 22,&
        5,  6, 21,&
       22, 21,  6,&
       33,  6,  7,&
       34, 22, 23,&
        6,  7, 22,&
       23, 22,  7,&
       33,  7,  8,&
       34, 23, 24,&
        7,  8, 23,&
       24, 23,  8,&
       33,  8,  9,&
       34, 24, 25,&
        8,  9, 24,&
       25, 24,  9,&
       33,  9, 10,&
       34, 25, 26,&
        9, 10, 25,&
       26, 25, 10,&
       33, 10, 11,&
       34, 26, 27,&
       10, 11, 26,&
       27, 26, 11,&
       33, 11, 12,&
       34, 27, 28,&
       11, 12, 27,&
       28, 27, 12,&
       33, 12, 13,&
       34, 28, 29,&
       12, 13, 28,&
       29, 28, 13,&
       33, 13, 14,&
       34, 29, 30,&
       13, 14, 29,&
       30, 29, 14,&
       33, 14, 15,&
       34, 30, 31,&
       14, 15, 30,&
       31, 30, 15,&
       33, 15, 16,&
       34, 31, 32,&
       15, 16, 31,&
       32, 31, 16,&
       33, 16,  1,&
       34, 32, 17,&
       16,  1, 32,&
       17, 32,  1/),(/3,64/))

    ! nullify the graphics database
    isinit = .true.

  end subroutine graphics_init

  !> Open a graphics file (file) with format fmt. Returns the logical
  !> unit and (in obj files) the mtl file.
  subroutine graphics_open(g,fmt,file)
    use tools_io, only: equal, lower
    class(grhandle), intent(inout) :: g
    character*3, intent(in) :: fmt
    character*(*), intent(in) :: file

    character*3 :: fmt0

    if (.not.isinit) then
       !$omp critical (initgraph)
       call graphics_init()
       !$omp end critical (initgraph)
    end if

    fmt0 = lower(fmt)
    if (.not.(equal(fmt0,"obj") .or. equal(fmt0,"ply") .or. equal(fmt0,"off"))) &
       return
    g%file = adjustl(file)

    if (equal(fmt0,"obj")) then
       g%fmt = ifmt_obj
       call obj_open(g)
    elseif (equal(fmt0,"ply")) then
       g%fmt = ifmt_ply
       call ply_open(g)
    elseif (equal(fmt0,"off")) then
       g%fmt = ifmt_off
       call off_open(g)
    end if

  end subroutine graphics_open

  !> Open a graphics file (lu, lumtl) with format fmt.
  subroutine graphics_close(g)
    use tools_io, only: equal
    class(grhandle), intent(inout) :: g

    if (g%fmt == ifmt_obj) then
       call obj_close(g)
    elseif (g%fmt == ifmt_ply) then 
       call ply_close(g)
    elseif (g%fmt == ifmt_off) then
       call off_close(g)
    end if
    g%lu = 0
    g%lumtl = 0
    g%fmt = ifmt_unk
    g%file = ""
    g%nball = 0
    g%nface = 0
    g%nstick = 0
    g%nsurf = 0
    g%nmtl = 0
    g%nv = 0
    g%nf = 0
    if (allocated(g%mtlrgb)) deallocate(g%mtlrgb)

  end subroutine graphics_close

  !> Write a ball to graphics file with LU lu and format fmt. The
  !> center is at x and the radius is r. rgb is the color.
  subroutine graphics_ball(g,x,rgb,r)
    use tools_io, only: equal
    class(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    if (g%fmt == ifmt_obj) then
       call obj_ball(g,x,rgb,r)
    elseif (g%fmt == ifmt_ply) then
       call ply_ball(g,x,rgb,r)
    elseif (g%fmt == ifmt_off) then
       call off_ball(g,x,rgb,r)
    end if

  end subroutine graphics_ball

  !> Write a polygon to graphics file with LU lu and format fmt.  The
  !> vertices are in x, and are assumed to be consecutive. rgb is the
  !> color.
  subroutine graphics_polygon(g,x,rgb)
    use tools_io, only: equal
    class(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    if (g%fmt == ifmt_obj) then
       call obj_polygon(g,x,rgb)
    elseif (g%fmt == ifmt_ply) then
       call ply_polygon(g,x,rgb)
    elseif (g%fmt == ifmt_off) then
       call off_polygon(g,x,rgb)
    end if

  end subroutine graphics_polygon

  !> Write a stick to graphics file with LU lu and format fmt.  The
  !> vertices are in x, and are assumed to be consecutive. rgb is the
  !> color.
  subroutine graphics_stick(g,x1,x2,rgb,r)
    use tools_io, only: equal
    class(grhandle), intent(inout) :: g
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    if (g%fmt == ifmt_obj) then
       call obj_stick(g,x1,x2,rgb,r)
    elseif (g%fmt == ifmt_ply) then
       call ply_stick(g,x1,x2,rgb,r)
    elseif (g%fmt == ifmt_off) then
       call off_stick(g,x1,x2,rgb,r)
    end if

  end subroutine graphics_stick

  !> Write a surface (srf with colors fsurf) to graphics file with LU
  !> lu and format fmt.
  subroutine graphics_surf(g,srf,fsurf)
    use tools_io, only: equal
    use surface, only: minisurf
    class(grhandle), intent(inout) :: g
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    if (g%fmt == ifmt_obj) then
       call obj_surf(g,srf,fsurf)
    elseif (g%fmt == ifmt_ply) then
       call ply_surf(g,srf,fsurf)
    elseif (g%fmt == ifmt_off) then
       call off_surf(g,srf,fsurf)
    end if

  end subroutine graphics_surf

  !> Open an obj file (and its mtl companion)
  subroutine obj_open(g)
    use tools_io, only: faterr, ferror, fopen_write
    type(grhandle), intent(inout) :: g

    integer :: idx
    character(len=:), allocatable :: filemtl, aux

    ! open the obj
    g%lu = fopen_write(g%file)

    ! name of the mtl
    filemtl = g%file
    idx = index(filemtl,'.',.true.)
    if (idx==0) call ferror("obj_open","could not parse file name",faterr)
    aux = trim(filemtl(1:idx-1)) // ".mtl"
    filemtl = aux

    ! open the mtl
    g%lumtl = fopen_write(filemtl)

    ! clear and initialize the mtl database
    g%nmtl = 0
    g%nball = 0
    g%nface = 0
    g%nstick = 0
    g%nv = 0
    if (allocated(g%mtlrgb)) deallocate(g%mtlrgb)
    allocate(g%mtlrgb(3,10))

    ! write the obj header
    write (g%lu,'("# OBJ created by critic2")')
    write (g%lu,'("mtllib ",A)') trim(filemtl)

  end subroutine obj_open

  !> Close an obj file (and its mtl companion)
  subroutine obj_close(g)
    use tools_io, only: ferror, string, fclose
    type(grhandle), intent(inout) :: g

    integer :: i

    ! write the mtl 
    do i = 1, g%nmtl
       write (g%lumtl,'("newmtl mat",A)') string(i)
       write (g%lumtl,'("Ns 96.078")')
       write (g%lumtl,'("Ka 0.0 0.0 0.0")')
       write (g%lumtl,'("Kd ",3(F12.5,X))') real(g%mtlrgb(:,i),8)/255d0
       write (g%lumtl,'("Ks 0.5 0.5 0.5")')
       write (g%lumtl,'("Ni 1.0")')
       write (g%lumtl,'("illum 2")')
       write (g%lumtl,*)
    end do

    ! clear the mtl database
    g%nmtl = 0
    g%nball = 0
    g%nface = 0
    g%nstick = 0
    g%nv = 0
    deallocate(g%mtlrgb)

    ! close both files
    call fclose(g%lu)
    call fclose(g%lumtl)

  end subroutine obj_close

  !> Write a ball to the obj file
  subroutine obj_ball(g,x,rgb,r)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: imtl, i, lvl

    ! ball detail
    if (r > 0.5) then
       lvl = 2
    elseif (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif

    ! write the ball to the obj
    g%nball = g%nball + 1
    write (g%lu,'("o ball_",A)') string(g%nball)
    write (g%lu,'("s on")')
    imtl = register_texture(g,rgb)
    write (g%lu,'("usemtl mat",A)') string(imtl)
    do i = 1, nvsph(lvl)
       write (g%lu,'("v ",3(F20.12,X))') x + r * vsph(:,i,lvl)
    end do
    do i = 1, nfsph(lvl)
       write (g%lu,'("f ",3(I10,X))') g%nv + fsph(:,i,lvl)
    end do
    g%nv = g%nv + nvsph(lvl)

  end subroutine obj_ball

  !> Write a polygon to the obj file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine obj_polygon(g,x,rgb)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: imtl, i, n
    character(len=:), allocatable :: str, aux

    ! write the face to the obj
    n = size(x,2)
    g%nface = g%nface + 1
    write (g%lu,'("o face_",A)') string(g%nface)
    write (g%lu,'("s on")')
    imtl = register_texture(g,rgb)
    write (g%lu,'("usemtl mat",A)') string(imtl)
    do i = 1, n
       write (g%lu,'("v ",3(F20.12,X))') x(:,i)
    end do
    str = ""
    do i = 1, n
       aux = str // string(g%nv+i) // " "
       str = aux
    end do
    write (g%lu,'("f ",A)') str
    g%nv = g%nv + n

  end subroutine obj_polygon
  
  !> Write a stick to the obj file
  subroutine obj_stick(g,x1,x2,rgb,r)
    use tools_io, only: ferror, string
    use tools_math, only: norm, cross
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: imtl, i, lvl
    real*8 :: dist, xd(3), v1(3), v2(3)

    ! stick detail
    if (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif

    ! kill degenerate sticks
    xd = x2 - x1
    dist = norm(xd)
    if (dist < 1d-6) return
    
    ! local coordinates of the tube
    if (abs(xd(1)) > 1d-6) then
       v1 = (/-xd(2)/xd(1), 1d0, 0d0/)
    elseif (abs(xd(2)) > 1d-6) then
       v1 = (/0d0, -xd(3)/xd(2), 1d0/)
    elseif (abs(xd(3)) > 1d-6) then
       v1 = (/1d0, 0d0, -xd(1)/xd(3)/)
    else
       return
    endif
    v2 = cross(xd,v1)
    v1 = v1 / norm(v1) * r
    v2 = v2 / norm(v2) * r

    ! write the stick to the obj
    g%nstick = g%nstick + 1
    write (g%lu,'("o stick_",A)') string(g%nstick)
    write (g%lu,'("s on")')
    imtl = register_texture(g,rgb)
    write (g%lu,'("usemtl mat",A)') string(imtl)
    do i = 1, nvcyl(lvl)
       write (g%lu,'("v ",3(F20.12,X))') x1+vcyl(1,i,lvl)*v1+&
          vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd
    end do
    do i = 1, nfcyl(lvl)
       write (g%lu,'("f ",3(I10,X))') g%nv + fcyl(:,i,lvl)
    end do
    g%nv = g%nv + nvcyl(lvl)

  end subroutine obj_stick
  
  !> Write a surface to the obj file
  subroutine obj_surf(g,srf,fsurf)
    use surface, only: minisurf
    use tools_io, only: ferror, faterr, string
    use param, only: pi
    type(grhandle), intent(inout) :: g
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, imtl
    real*8 :: x(3), maxf, minf, xrgb(3), z

    if (srf%isinit <= 1) &
       call ferror ('obj_surf','No face information in minisurf',faterr)

    ! limits of the color scale
    if (present(fsurf)) then
       if (size(fsurf) == 3) then
          xrgb = fsurf
       else
          maxf = maxval(fsurf)
          minf = minval(fsurf)
       endif
    endif

    ! write the surface to the obj
    g%nsurf = g%nsurf + 1
    write (g%lu,'("o surf_",A)') string(g%nsurf)
    write (g%lu,'("s on")')
    imtl = register_texture(g,rgb_default)
    write (g%lu,'("usemtl mat",A)') string(imtl)
    do i = 1, srf%nv
       x = srf%n + (/ srf%r(i) * sin(srf%th(i)) * cos(srf%ph(i)),&
           srf%r(i) * sin(srf%th(i)) * sin(srf%ph(i)),&
           srf%r(i) * cos(srf%th(i)) /)
       if (present(fsurf)) then
          if (size(fsurf) /= 3) then
             z = (fsurf(i) - minf) / (maxf - minf)
             xrgb(1) = sqrt(z)
             xrgb(2) = z**3
             xrgb(3) = sin(2*z*pi)
          end if
          write (g%lu,'("v ",3(E20.12,X),3(F8.5,X))') x, xrgb
       else
          write (g%lu,'("v ",3(E20.12,X))') x
       endif
    end do
    do i = 1, srf%nf
       write (g%lu,'("f ",999(I10,X))') &
          (g%nv+srf%f(i)%v(j),j=1,srf%f(i)%nv)
    end do
    g%nv = g%nv + srf%nv

  end subroutine obj_surf

  !> Register a texture using the color triplet.
  function register_texture(g,rgb) result(imtl)
    use types, only: realloc
    type(grhandle), intent(inout) :: g
    integer, intent(in) :: rgb(3)
    integer :: imtl

    integer :: i

    do i = 1, g%nmtl
       if (all(g%mtlrgb(:,i)-rgb == 0)) then
          imtl = i
          return
       end if
    end do
    g%nmtl = g%nmtl + 1
    if (g%nmtl > size(g%mtlrgb,2)) call realloc(g%mtlrgb,3,2*g%nmtl)
    g%mtlrgb(:,g%nmtl) = rgb
    imtl = g%nmtl

  end function register_texture

  !> Open a ply file
  subroutine ply_open(g)
    use tools_io, only: ferror, fopen_scratch
    type(grhandle), intent(inout) :: g

    ! open the temporary ply file
    g%lu = fopen_scratch("formatted")

    ! clear and initialize the database
    g%nv = 0
    g%nf = 0

  end subroutine ply_open

  !> Close a ply file 
  subroutine ply_close(g)
    use tools_io, only: ferror, getline_raw, fopen_write, string, fclose
    type(grhandle), intent(inout) :: g

    integer :: lu
    character(len=:), allocatable :: file, line
    logical :: ok

    ! get the temporary file name from the first line and open
    rewind(g%lu)
    lu = fopen_write(g%file)

    ! write the header
    write (lu,'("ply ")')
    write (lu,'("format ascii 1.0")')
    write (lu,'("comment PLY created by critic2")')
    write (lu,'("element vertex ",A)') string(g%nv)
    write (lu,'("property float x")')
    write (lu,'("property float y")')
    write (lu,'("property float z")')
    write (lu,'("property uchar red")')
    write (lu,'("property uchar green")')
    write (lu,'("property uchar blue")')
    write (lu,'("property uchar alpha")')
    write (lu,'("element face ",A)') string(g%nf)
    write (lu,'("property list uchar int vertex_indices")')
    write (lu,'("end_header")')

    ! transfer the vertices from the scratch file over to the ply file
    do while(getline_raw(g%lu,line))
       if (line(1:1) == "v") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! transfer the faces
    rewind(g%lu)
    ok = getline_raw(g%lu,file,.true.)
    do while(getline_raw(g%lu,line))
       if (line(1:1) == "f") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! close both
    call fclose(lu)
    call fclose(g%lu)

    ! clear the database
    g%nv = 0
    g%nf = 0
    
  end subroutine ply_close

  !> Write a ball to the ply file
  subroutine ply_ball(g,x,rgb,r)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)

    ! ball detail
    if (r > 0.5) then
       lvl = 2
    elseif (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif
    
    ! write the ball to the ply
    do i = 1, nvsph(lvl)
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') x + r * vsph(:,i,lvl), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, nfsph(lvl)
       iff = g%nv + fsph(:,i,lvl) - 1
       write (g%lu,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    g%nv = g%nv + nvsph(lvl)
    g%nf = g%nf + nfsph(lvl)

  end subroutine ply_ball
  
  !> Write a polygon to the ply file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine ply_polygon(g,x,rgb)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: i, i1, n
    real*8 :: xcm(3)

    ! center of the face
    n = size(x,2)
    xcm = 0d0
    do i = 1, n
       xcm = xcm + x(:,i)
    end do
    xcm = xcm / n

    ! write the face to the obj
    do i = 1, n
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') x(:,i), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') xcm, &
       string(rgb(1)), string(rgb(2)), string(rgb(3))

    do i = 1, n
       i1 = mod(i,n) + 1
       write (g%lu,'("f3 ",3(A,X))') string(g%nv+i-1), string(g%nv+i1-1), string(g%nv+n)
    end do
    g%nv = g%nv + n + 1
    g%nf = g%nf + n

  end subroutine ply_polygon

  !> Write a stick to the ply file
  subroutine ply_stick(g,x1,x2,rgb,r)
    use tools_io, only: ferror, string
    use tools_math, only: norm, cross
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: dist, xd(3), v1(3), v2(3)

    ! stick detail
    if (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif
    
    ! kill degenerate sticks
    xd = x2 - x1
    dist = norm(xd)
    if (dist < 1d-6) return
    
    ! local coordinates of the tube
    if (abs(xd(1)) > 1d-6) then
       v1 = (/-xd(2)/xd(1), 1d0, 0d0/)
    elseif (abs(xd(2)) > 1d-6) then
       v1 = (/0d0, -xd(3)/xd(2), 1d0/)
    elseif (abs(xd(3)) > 1d-6) then
       v1 = (/1d0, 0d0, -xd(1)/xd(3)/)
    else
       return
    endif
    v2 = cross(xd,v1)
    v1 = v1 / norm(v1) * r
    v2 = v2 / norm(v2) * r

    ! write the stick to the ply
    do i = 1, nvcyl(lvl)
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') &
          x1+vcyl(1,i,lvl)*v1+vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd,&
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, nfcyl(lvl)
       iff = g%nv + fcyl(:,i,lvl) - 1
       write (g%lu,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    g%nv = g%nv + nvcyl(lvl)
    g%nf = g%nf + nfcyl(lvl)

  end subroutine ply_stick
  
  !> Write a surface to the ply file
  subroutine ply_surf(g,srf,fsurf)
    use surface, only: minisurf
    use tools_io, only: ferror, faterr, string
    use param, only: pi
    type(grhandle), intent(inout) :: g
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, rgb(3)
    real*8 :: x(3), maxf, minf, xrgb(3), z
    
    if (srf%isinit <= 1) &
       call ferror ('ply_surf','No face information in minisurf',faterr)

    ! limits of the color scale
    rgb = rgb_default
    if (present(fsurf)) then
       if (size(fsurf) /= 3) then
          maxf = maxval(fsurf)
          minf = minval(fsurf)
       else
          rgb = nint(fsurf * 255)
       end if
    endif

    do i = 1, srf%nv
       x = srf%n + (/ srf%r(i) * sin(srf%th(i)) * cos(srf%ph(i)),&
           srf%r(i) * sin(srf%th(i)) * sin(srf%ph(i)),&
           srf%r(i) * cos(srf%th(i)) /)
       if (present(fsurf)) then
          if (size(fsurf) > 3) then
             z = (fsurf(i) - minf) / (maxf - minf)
             xrgb(1) = sqrt(z)
             xrgb(2) = z**3
             xrgb(3) = sin(2*z*pi)
             rgb = nint(xrgb * 255)
          end if
       endif
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') x,&
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, srf%nf
       write (g%lu,'("f",999(A,X))') string(srf%f(i)%nv),&
          (string(g%nv+srf%f(i)%v(j)-1),j=1,srf%f(i)%nv)
    end do
    g%nv = g%nv + srf%nv
    g%nf = g%nf + srf%nf

  end subroutine ply_surf

  !> Open an off file
  subroutine off_open(g)
    use tools_io, only: ferror, fopen_scratch
    type(grhandle), intent(inout) :: g

    ! open the temporary off file
    g%lu = fopen_scratch("formatted")

    ! clear and initialize the database
    g%nv = 0
    g%nf = 0

  end subroutine off_open

  !> Close a off file 
  subroutine off_close(g)
    use tools_io, only: ferror, getline_raw, string, fopen_write, fclose
    type(grhandle), intent(inout) :: g

    integer :: lu
    character(len=:), allocatable :: file, line
    logical :: ok

    ! get the temporary file name from the first line and open
    rewind(g%lu)
    lu = fopen_write(g%file)

    ! write the header
    write (lu,'("COFF")')
    write (lu,'(3(A,X))') string(g%nv), string(g%nf), string(g%nv+g%nf-2)

    ! transfer the vertices from the scratch file over to the off file
    do while(getline_raw(g%lu,line))
       if (line(1:1) == "v") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! transfer the faces
    rewind(g%lu)
    ok = getline_raw(g%lu,file,.true.)
    do while(getline_raw(g%lu,line))
       if (line(1:1) == "f") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! close both
    call fclose(lu)
    call fclose(g%lu)

    ! clear the database
    g%nv = 0
    g%nf = 0
    
  end subroutine off_close

  !> Write a ball to the off file
  subroutine off_ball(g,x,rgb,r)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: xrgb(3)

    ! ball detail
    if (r > 0.5) then
       lvl = 2
    elseif (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif
    
    ! write the ball to the off
    do i = 1, nvsph(lvl)
       xrgb = real(rgb,8) / 255d0
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0.0")') x + r * vsph(:,i,lvl), &
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, nfsph(lvl)
       iff = g%nv + fsph(:,i,lvl) - 1
       write (g%lu,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    g%nv = g%nv + nvsph(lvl)
    g%nf = g%nf + nfsph(lvl)

  end subroutine off_ball
  
  !> Write a polygon to the off file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine off_polygon(g,x,rgb)
    use tools_io, only: ferror, string
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: i, i1, n
    real*8 :: xcm(3)

    ! center of the face
    n = size(x,2)
    xcm = 0d0
    do i = 1, n
       xcm = xcm + x(:,i)
    end do
    xcm = xcm / n

    ! write the face to the obj
    do i = 1, n
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') x(:,i), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') xcm, &
       string(rgb(1)), string(rgb(2)), string(rgb(3))

    do i = 1, n
       i1 = mod(i,n) + 1
       write (g%lu,'("f3 ",3(A,X))') string(g%nv+i-1), string(g%nv+i1-1), string(g%nv+n)
    end do
    g%nv = g%nv + n + 1
    g%nf = g%nf + n

  end subroutine off_polygon

  !> Write a stick to the off file
  subroutine off_stick(g,x1,x2,rgb,r)
    use tools_io, only: ferror, string
    use tools_math, only: norm, cross
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: dist, xd(3), v1(3), v2(3)
    real*8 :: xrgb(3)

    ! stick detail
    if (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif
    
    ! kill degenerate sticks
    xd = x2 - x1
    dist = norm(xd)
    if (dist < 1d-6) return
    
    ! local coordinates of the tube
    if (abs(xd(1)) > 1d-6) then
       v1 = (/-xd(2)/xd(1), 1d0, 0d0/)
    elseif (abs(xd(2)) > 1d-6) then
       v1 = (/0d0, -xd(3)/xd(2), 1d0/)
    elseif (abs(xd(3)) > 1d-6) then
       v1 = (/1d0, 0d0, -xd(1)/xd(3)/)
    else
       return
    endif
    v2 = cross(xd,v1)
    v1 = v1 / norm(v1) * r
    v2 = v2 / norm(v2) * r

    ! write the stick to the off
    do i = 1, nvcyl(lvl)
       xrgb = real(rgb,8) / 255d0
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') &
          x1+vcyl(1,i,lvl)*v1+vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd,&
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, nfcyl(lvl)
       iff = g%nv + fcyl(:,i,lvl) - 1
       write (g%lu,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    g%nv = g%nv + nvcyl(lvl)
    g%nf = g%nf + nfcyl(lvl)

  end subroutine off_stick
  
  !> Write a surface to the off file
  subroutine off_surf(g,srf,fsurf)
    use surface, only: minisurf
    use tools_io, only: faterr, ferror, string
    use param, only: pi
    type(grhandle), intent(inout) :: g
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, rgb(3)
    real*8 :: x(3), z
    real*8 :: xrgb(3), maxf, minf

    if (srf%isinit <= 1) &
       call ferror ('off_surf','No face information in minisurf',faterr)

    ! limits of the color scale
    minf = 0d0
    rgb = rgb_default
    xrgb = real(rgb,8) / 255d0
    if (present(fsurf)) then
       if (size(fsurf) /= 3) then
          maxf = maxval(fsurf)
          minf = minval(fsurf)
       else
          xrgb = fsurf
       endif
    endif

    do i = 1, srf%nv
       x = srf%n + (/ srf%r(i) * sin(srf%th(i)) * cos(srf%ph(i)),&
           srf%r(i) * sin(srf%th(i)) * sin(srf%ph(i)),&
           srf%r(i) * cos(srf%th(i)) /)
       if (present(fsurf)) then
          if (size(fsurf) /= 3) then
             z = (fsurf(i) - minf) / (maxf - minf)
             xrgb(1) = sqrt(z)
             xrgb(2) = z**3
             xrgb(3) = sin(2*z*pi)
          endif
       endif
       write (g%lu,'("v",3(F20.12,X),3(A,X),"0")') x,&
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, srf%nf
       write (g%lu,'("f",999(A,X))') string(srf%f(i)%nv),&
          (string(g%nv+srf%f(i)%v(j)-1),j=1,srf%f(i)%nv)
    end do
    g%nv = g%nv + srf%nv
    g%nf = g%nf + srf%nf

  end subroutine off_surf

end module graphics
