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
  public :: graphics_init
  public :: writexyz
  public :: writecml
  public :: writegjf
  public :: graphics_open
  public :: graphics_close
  public :: graphics_ball
  public :: graphics_polygon
  public :: graphics_stick
  public :: graphics_surf
  private :: graphics_checkfmt
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
  logical :: isopen
  integer :: nball, nstick, nsurf, nface
  integer :: nmtl
  integer, allocatable :: mtlrgb(:,:)
  integer :: nv, nf

  ! icosahedron models
  integer :: nvsph(0:2), nfsph(0:2)
  real*8 :: vsph(3,42,0:2)
  integer :: fsph(3,80,0:2)

  ! cylinder models
  integer :: nvcyl(0:1), nfcyl(0:1)
  real*8 :: vcyl(3,34,0:1)
  integer :: fcyl(3,64,0:1)

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
    isopen = .false.
    nball = 0
    nface = 0
    nstick = 0
    nsurf = 0
    nmtl = 0
    nv = 0
    nf = 0

  end subroutine graphics_init

  !> write an xyz-style file from an array of atomic coordinates.
  subroutine writexyz(file,fr)
    use tools_io, only: fopen_write, nameguess, fclose
    use types, only: fragment
    use param, only: bohrtoa
    character*(*), intent(in) :: file
    type(fragment), intent(in) :: fr

    integer :: i, lu

    ! write it
    lu = fopen_write(file)
    write (lu,*) fr%nat
    write (lu,*)
    do i = 1, fr%nat
       if (fr%at(i)%z >= 0) then
          write (lu,*) nameguess(fr%at(i)%z,.true.), fr%at(i)%r * bohrtoa
       end if
    end do
    call fclose(lu)

  end subroutine writexyz

  !> Write a cml file (molecule) from an array of atomic coordinates. 
  subroutine writecml(file,fr,r,luout)
    use tools_math, only: matinv
    use tools_io, only: fopen_write, string, nameguess, fclose
    use types, only: fragment
    use param, only: pi, bohrtoa
    character*(*), intent(in) :: file
    type(fragment), intent(in) :: fr
    real*8, intent(in), optional :: r(3,3)
    integer, intent(out), optional :: luout

    integer :: i, j, lu
    real*8 :: g(3,3), aa(3), bb(3), x(3), ri(3,3)

    ! write it
    lu = fopen_write(file)
    write (lu,'("<molecule>")')

    ! crystal structure
    if (present(r)) then
       ri = matinv(r)
       g = matmul(transpose(r),r)
       do i = 1, 3
          aa(i) = sqrt(g(i,i)) 
       end do
       bb(1) = acos(g(2,3) / aa(2) / aa(3)) * 180d0 / pi
       bb(2) = acos(g(1,3) / aa(1) / aa(3)) * 180d0 / pi
       bb(3) = acos(g(1,2) / aa(1) / aa(2)) * 180d0 / pi
       aa = aa * bohrtoa
       write (lu,'(" <crystal>")')
       write (lu,'("  <scalar title=""a"" units=""units:angstrom"">",A,"</scalar>")') string(aa(1),'f',decimal=8)
       write (lu,'("  <scalar title=""b"" units=""units:angstrom"">",A,"</scalar>")') string(aa(2),'f',decimal=8)
       write (lu,'("  <scalar title=""c"" units=""units:angstrom"">",A,"</scalar>")') string(aa(3),'f',decimal=8)
       write (lu,'("  <scalar title=""alpha"" units=""units:degree"">",A,"</scalar>")') string(bb(1),'f',decimal=4)
       write (lu,'("  <scalar title=""beta"" units=""units:degree"">",A,"</scalar>")') string(bb(2),'f',decimal=4)
       write (lu,'("  <scalar title=""gamma"" units=""units:degree"">",A,"</scalar>")') string(bb(3),'f',decimal=4)
       write (lu,'("  <symmetry spaceGroup=""P 1"">")')
       write (lu,'("   <transform3>1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1</transform3>")')
       write (lu,'("  </symmetry>")')
       write (lu,'(" </crystal>")')
    end if

    write (lu,'(" <atomArray>")')
    do i = 1, fr%nat
       if (fr%at(i)%z >= 0) then
          if (present(r)) then
             x = matmul(ri,fr%at(i)%r)
             write (lu,'("<atom id=""a",A,""" elementType=""",A,""" xFract=""",A,""" yFract=""",A,""" zFract=""",A,"""/>")') &
                string(i), trim(nameguess(fr%at(i)%z,.true.)), (trim(string(x(j),'f',18,10)),j=1,3)
          else
             write (lu,'("<atom id=""a",A,""" elementType=""",A,""" x3=""",A,""" y3=""",A,""" z3=""",A,"""/>")') &
                string(i), trim(nameguess(fr%at(i)%z,.true.)), &
                (trim(string(fr%at(i)%r(j) * bohrtoa,'f',18,10)),j=1,3)
          end if
       end if
    end do
    if (present(luout)) then
       luout = lu
    else
       write (lu,'(" </atomArray>")')
       write (lu,'("</molecule>")')
       call fclose(lu)
    end if

  end subroutine writecml

  !> write an Gaussian-style input file from an array of atomic coordinates.
  subroutine writegjf(file,fr)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use types, only: fragment
    use param, only: bohrtoa

    character*(*), intent(in) :: file
    type(fragment), intent(in) :: fr

    character(len=:), allocatable :: aux
    integer :: i, lu, isum

    aux = file

    ! write it
    lu = fopen_write(aux)
    write (lu,'("#p b3lyp sto-3g")')
    write (lu,'("")')
    write (lu,'("title")')
    write (lu,'("")')

    isum = 0
    do i = 1, fr%nat
       if (fr%at(i)%z > 0) isum = isum + fr%at(i)%z
    end do
    write (lu,'("0 ",A)') string(2*modulo(isum,2)+1)
    do i = 1, fr%nat
       if (fr%at(i)%z > 0) then
          write (lu,*) nameguess(fr%at(i)%z,.true.), fr%at(i)%r * bohrtoa
       end if
    end do
    write (lu,'("")')
    call fclose(lu)

  end subroutine writegjf

  !> Open a graphics file (file) with format fmt. Returns the logical
  !> unit and (in obj files) the mtl file.
  subroutine graphics_open(fmt,file,lu,lumtl)
    use tools_io, only: faterr, ferror, equal
    character*3, intent(in) :: fmt
    character*(*), intent(in) :: file
    integer, intent(out) :: lu
    integer, intent(out), optional :: lumtl

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       if (.not.present(lumtl)) &
          call ferror("graphics_open","mtl LU argument not found",faterr)
       call obj_open(file,lu,lumtl)
    elseif (equal(fmt,"ply")) then
       call ply_open(file,lu)
    elseif (equal(fmt,"off")) then
       call off_open(file,lu)
    end if

  end subroutine graphics_open

  !> Open a graphics file (lu, lumtl) with format fmt.
  subroutine graphics_close(fmt,lu,lumtl)
    use tools_io, only: ferror, faterr, equal
    character*3, intent(in) :: fmt
    integer, intent(in) :: lu
    integer, intent(in), optional :: lumtl

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       if (.not.present(lumtl)) &
          call ferror("graphics_close","mtl LU argument not found",faterr)
       call obj_close(lu,lumtl)
    elseif (equal(fmt,"ply")) then
       call ply_close(lu)
    elseif (equal(fmt,"off")) then
       call off_close(lu)
    end if

  end subroutine graphics_close

  !> Write a ball to graphics file with LU lu and format fmt.  The
  !> center is at x and the radius is r. rgb is the color.
  subroutine graphics_ball(fmt,lu,x,rgb,r)
    use tools_io, only: equal
    character*3, intent(in) :: fmt
    integer, intent(in) :: lu
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       call obj_ball(lu,x,rgb,r)
    elseif (equal(fmt,"ply")) then
       call ply_ball(lu,x,rgb,r)
    elseif (equal(fmt,"off")) then
       call off_ball(lu,x,rgb,r)
    end if

  end subroutine graphics_ball

  !> Write a polygon to graphics file with LU lu and format fmt.  The
  !> vertices are in x, and are assumed to be consecutive. rgb is the
  !> color.
  subroutine graphics_polygon(fmt,lu,x,rgb)
    use tools_io, only: equal
    character*3, intent(in) :: fmt
    integer, intent(in) :: lu
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       call obj_polygon(lu,x,rgb)
    elseif (equal(fmt,"ply")) then
       call ply_polygon(lu,x,rgb)
    elseif (equal(fmt,"off")) then
       call off_polygon(lu,x,rgb)
    end if

  end subroutine graphics_polygon

  !> Write a stick to graphics file with LU lu and format fmt.  The
  !> vertices are in x, and are assumed to be consecutive. rgb is the
  !> color.
  subroutine graphics_stick(fmt,lu,x1,x2,rgb,r)
    use tools_io, only: equal
    character*3, intent(in) :: fmt
    integer, intent(in) :: lu
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       call obj_stick(lu,x1,x2,rgb,r)
    elseif (equal(fmt,"ply")) then
       call ply_stick(lu,x1,x2,rgb,r)
    elseif (equal(fmt,"off")) then
       call off_stick(lu,x1,x2,rgb,r)
    end if

  end subroutine graphics_stick

  !> Write a surface (srf with colors fsurf) to graphics file with LU
  !> lu and format fmt.
  subroutine graphics_surf(fmt,lu,srf,fsurf)
    use tools_io, only: equal
    use types, only: minisurf
    character*3, intent(in) :: fmt
    integer, intent(in) :: lu
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    call graphics_checkfmt(fmt)
    if (equal(fmt,"obj")) then
       call obj_surf(lu,srf,fsurf)
    elseif (equal(fmt,"ply")) then
       call ply_surf(lu,srf,fsurf)
    elseif (equal(fmt,"off")) then
       call off_surf(lu,srf,fsurf)
    end if

  end subroutine graphics_surf

  subroutine graphics_checkfmt(fmt) 
    use tools_io, only: equal, ferror, faterr
    character*3, intent(in) :: fmt
    logical :: ok

    if (.not.(equal(fmt,"obj") .or. equal(fmt,"ply") .or. equal(fmt,"off"))) &
       call ferror("graphics_checkfmt","unknown graphics format "//fmt,faterr)

  end subroutine graphics_checkfmt

  !> Open an obj file (and its mtl companion)
  subroutine obj_open(file,luobj,lumtl)
    use tools_io, only: faterr, ferror, fopen_write
    character*(*), intent(in) :: file
    integer, intent(out) :: luobj, lumtl

    integer :: idx
    character*(len(file)) :: filemtl

    if (isopen) &
       call ferror('obj_open','error opening graphics file: one is already open',faterr)
    isopen = .true.

    ! open the obj
    luobj = fopen_write(file)

    ! name of the mtl
    filemtl = file
    idx = index(filemtl,'.',.true.)
    if (idx==0) call ferror("obj_open","could not parse file name",faterr)
    filemtl = trim(filemtl(1:idx-1)) // ".mtl"

    ! open the mtl
    lumtl = fopen_write(filemtl)

    ! clear and initialize the mtl database
    nmtl = 0
    nball = 0
    nface = 0
    nstick = 0
    nv = 0
    if (allocated(mtlrgb)) deallocate(mtlrgb)
    allocate(mtlrgb(3,10))

    ! write the obj header
    write (luobj,'("# OBJ created by critic2")')
    write (luobj,'("mtllib ",A)') trim(filemtl)

  end subroutine obj_open

  !> Close an obj file (and its mtl companion)
  subroutine obj_close(luobj,lumtl)
    use tools_io, only: ferror, faterr, string, fclose
    integer, intent(in) :: luobj, lumtl

    integer :: i

    if (.not.isopen) &
       call ferror('obj_close','error: graphics file is not open',faterr)

    ! write the mtl 
    do i = 1, nmtl
       write (lumtl,'("newmtl mat",A)') string(i)
       write (lumtl,'("Ns 96.078")')
       write (lumtl,'("Ka 0.0 0.0 0.0")')
       write (lumtl,'("Kd ",3(F12.5,X))') real(mtlrgb(:,i),8)/255d0
       write (lumtl,'("Ks 0.5 0.5 0.5")')
       write (lumtl,'("Ni 1.0")')
       write (lumtl,'("illum 2")')
       write (lumtl,*)
    end do

    ! clear the mtl database
    nmtl = 0
    nball = 0
    nface = 0
    nstick = 0
    nv = 0
    deallocate(mtlrgb)

    ! close both files
    call fclose(luobj)
    call fclose(lumtl)
    isopen = .false.

  end subroutine obj_close

  !> Write a ball to the obj file
  subroutine obj_ball(luobj,x,rgb,r)
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: luobj
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: imtl, i, lvl

    if (.not.isopen) &
       call ferror('obj_ball','error: graphics file is not open',faterr)

    ! ball detail
    if (r > 0.5) then
       lvl = 2
    elseif (r > 0.2) then
       lvl = 1
    else
       lvl = 0
    endif

    ! write the ball to the obj
    nball = nball+1
    write (luobj,'("o ball_",A)') string(nball)
    write (luobj,'("s on")')
    imtl = register_texture(rgb)
    write (luobj,'("usemtl mat",A)') string(imtl)
    do i = 1, nvsph(lvl)
       write (luobj,'("v ",3(F20.12,X))') x + r * vsph(:,i,lvl)
    end do
    do i = 1, nfsph(lvl)
       write (luobj,'("f ",3(I10,X))') nv + fsph(:,i,lvl)
    end do
    nv = nv + nvsph(lvl)

  end subroutine obj_ball

  !> Write a polygon to the obj file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine obj_polygon(luobj,x,rgb)
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: luobj
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: imtl, i, n
    character(len=:), allocatable :: str, aux

    if (.not.isopen) &
       call ferror('obj_polygon','error: graphics file is not open',faterr)

    ! write the face to the obj
    n = size(x,2)
    nface = nface+1
    write (luobj,'("o face_",A)') string(nface)
    write (luobj,'("s on")')
    imtl = register_texture(rgb)
    write (luobj,'("usemtl mat",A)') string(imtl)
    do i = 1, n
       write (luobj,'("v ",3(F20.12,X))') x(:,i)
    end do
    str = ""
    do i = 1, n
       aux = str // string(nv+i) // " "
       str = aux
    end do
    write (luobj,'("f ",A)') str
    nv = nv + n

  end subroutine obj_polygon
  
  !> Write a stick to the obj file
  subroutine obj_stick(luobj,x1,x2,rgb,r)
    use tools_io, only: ferror, faterr, string
    use tools_math, only: norm, cross
    integer, intent(in) :: luobj
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: imtl, i, lvl
    real*8 :: dist, xd(3), v1(3), v2(3)

    if (.not.isopen) &
       call ferror('obj_stick','error: graphics file is not open',faterr)

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
    nstick = nstick + 1
    write (luobj,'("o stick_",A)') string(nstick)
    write (luobj,'("s on")')
    imtl = register_texture(rgb)
    write (luobj,'("usemtl mat",A)') string(imtl)
    do i = 1, nvcyl(lvl)
       write (luobj,'("v ",3(F20.12,X))') x1+vcyl(1,i,lvl)*v1+&
          vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd
    end do
    do i = 1, nfcyl(lvl)
       write (luobj,'("f ",3(I10,X))') nv + fcyl(:,i,lvl)
    end do
    nv = nv + nvcyl(lvl)

  end subroutine obj_stick
  
  !> Write a surface to the obj file
  subroutine obj_surf(luobj,srf,fsurf)
    use types, only: minisurf
    use tools_io, only: ferror, faterr, string
    use param, only: pi
    integer, intent(in) :: luobj
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, imtl
    real*8 :: x(3), maxf, minf, xrgb(3), z

    if (srf%init <= 1) &
       call ferror ('obj_surf','No face information in minisurf',faterr)
    if (.not.isopen) &
       call ferror('obj_surf','error: graphics file is not open',faterr)

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
    nsurf = nsurf + 1
    write (luobj,'("o surf_",A)') string(nsurf)
    write (luobj,'("s on")')
    imtl = register_texture(rgb_default)
    write (luobj,'("usemtl mat",A)') string(imtl)
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
          write (luobj,'("v ",3(E20.12,X),3(F8.5,X))') x, xrgb
       else
          write (luobj,'("v ",3(E20.12,X))') x
       endif
    end do
    do i = 1, srf%nf
       write (luobj,'("f ",999(I10,X))') &
          (nv+srf%f(i)%v(j),j=1,srf%f(i)%nv)
    end do
    nv = nv + srf%nv

  end subroutine obj_surf

  !> Register a texture using the color triplet.
  function register_texture(rgb) result(imtl)
    use types, only: realloc
    integer, intent(in) :: rgb(3)
    integer :: imtl

    integer :: i

    do i = 1, nmtl
       if (all(mtlrgb(:,i)-rgb == 0)) then
          imtl = i
          return
       end if
    end do
    nmtl = nmtl + 1
    if (nmtl > size(mtlrgb,2)) call realloc(mtlrgb,3,2*nmtl)
    mtlrgb(:,nmtl) = rgb
    imtl = nmtl

  end function register_texture

  !> Open a ply file
  subroutine ply_open(file,luply)
    use tools_io, only: faterr, ferror, fopen_scratch
    character*(*), intent(in) :: file
    integer, intent(out) :: luply

    if (isopen) &
       call ferror('ply_open','error opening graphics file: one is already open',faterr)
    isopen = .true.

    ! open the temporary ply file
    luply = fopen_scratch("formatted")

    ! write the target file name to the scratch file
    write (luply,'(A)') trim(file)

    ! clear and initialize the database
    nv = 0
    nf = 0

  end subroutine ply_open

  !> Close a ply file 
  subroutine ply_close(luply)
    use tools_io, only: ferror, faterr, getline_raw, fopen_write, string, fclose
    integer, intent(in) :: luply

    integer :: lu
    character(len=:), allocatable :: file, line
    logical :: ok

    if (.not.isopen) &
       call ferror('ply_close','error: graphics file is not open',faterr)

    ! get the temporary file name from the first line and open
    rewind(luply)
    ok = getline_raw(luply,file,.true.)
    lu = fopen_write(file)

    ! write the header
    write (lu,'("ply ")')
    write (lu,'("format ascii 1.0")')
    write (lu,'("comment PLY created by critic2")')
    write (lu,'("element vertex ",A)') string(nv)
    write (lu,'("property float x")')
    write (lu,'("property float y")')
    write (lu,'("property float z")')
    write (lu,'("property uchar red")')
    write (lu,'("property uchar green")')
    write (lu,'("property uchar blue")')
    write (lu,'("property uchar alpha")')
    write (lu,'("element face ",A)') string(nf)
    write (lu,'("property list uchar int vertex_indices")')
    write (lu,'("end_header")')

    ! transfer the vertices from the scratch file over to the ply file
    do while(getline_raw(luply,line))
       if (line(1:1) == "v") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! transfer the faces
    rewind(luply)
    ok = getline_raw(luply,file,.true.)
    do while(getline_raw(luply,line))
       if (line(1:1) == "f") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! close both
    call fclose(lu)
    call fclose(luply)
    isopen = .false.

    ! clear the database
    nv = 0
    nf = 0
    
  end subroutine ply_close

  !> Write a ball to the ply file
  subroutine ply_ball(luply,x,rgb,r)
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: luply
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)

    if (.not.isopen) &
       call ferror('ply_ball','error: graphics file is not open',faterr)

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
       write (luply,'("v",3(F20.12,X),3(A,X),"0")') x + r * vsph(:,i,lvl), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, nfsph(lvl)
       iff = nv + fsph(:,i,lvl) - 1
       write (luply,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    nv = nv + nvsph(lvl)
    nf = nf + nfsph(lvl)

  end subroutine ply_ball
  
  !> Write a polygon to the ply file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine ply_polygon(luply,x,rgb)
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: luply
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: i, i1, n
    character(len=:), allocatable :: str, aux
    real*8 :: xcm(3)

    if (.not.isopen) &
       call ferror('ply_polygon','error: graphics file is not open',faterr)

    ! center of the face
    n = size(x,2)
    xcm = 0d0
    do i = 1, n
       xcm = xcm + x(:,i)
    end do
    xcm = xcm / n

    ! write the face to the obj
    do i = 1, n
       write (luply,'("v",3(F20.12,X),3(A,X),"0")') x(:,i), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    write (luply,'("v",3(F20.12,X),3(A,X),"0")') xcm, &
       string(rgb(1)), string(rgb(2)), string(rgb(3))

    do i = 1, n
       i1 = mod(i,n) + 1
       write (luply,'("f3 ",3(A,X))') string(nv+i-1), string(nv+i1-1), string(nv+n)
    end do
    nv = nv + n + 1
    nf = nf + n

  end subroutine ply_polygon

  !> Write a stick to the ply file
  subroutine ply_stick(luply,x1,x2,rgb,r)
    use tools_io, only: ferror, faterr, string
    use tools_math, only: norm, cross
    integer, intent(in) :: luply
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: dist, xd(3), v1(3), v2(3)

    if (.not.isopen) &
       call ferror('ply_stick','error: graphics file is not open',faterr)

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
       write (luply,'("v",3(F20.12,X),3(A,X),"0")') &
          x1+vcyl(1,i,lvl)*v1+vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd,&
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, nfcyl(lvl)
       iff = nv + fcyl(:,i,lvl) - 1
       write (luply,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    nv = nv + nvcyl(lvl)
    nf = nf + nfcyl(lvl)

  end subroutine ply_stick
  
  !> Write a surface to the ply file
  subroutine ply_surf(luply,srf,fsurf)
    use types, only: minisurf
    use tools_io, only: ferror, faterr, string
    use param, only: pi
    integer, intent(in) :: luply
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, rgb(3)
    real*8 :: x(3), maxf, minf, xrgb(3), z
    
    if (srf%init <= 1) &
       call ferror ('ply_surf','No face information in minisurf',faterr)
    if (.not.isopen) &
       call ferror('ply_surf','error: graphics file is not open',faterr)

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
       write (luply,'("v",3(F20.12,X),3(A,X),"0")') x,&
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    do i = 1, srf%nf
       write (luply,'("f",999(A,X))') string(srf%f(i)%nv),&
          (string(nv+srf%f(i)%v(j)-1),j=1,srf%f(i)%nv)
    end do
    nv = nv + srf%nv
    nf = nf + srf%nf

  end subroutine ply_surf

  !> Open an off file
  subroutine off_open(file,luoff)
    use tools_io, only: ferror, faterr, fopen_scratch
    character*(*), intent(in) :: file
    integer, intent(out) :: luoff

    if (isopen) &
       call ferror('off_open','error opening graphics file: one is already open',faterr)
    isopen = .true.

    ! open the temporary off file
    luoff = fopen_scratch("formatted")

    ! write the target file name to the scratch file
    write (luoff,'(A)') trim(file)

    ! clear and initialize the database
    nv = 0
    nf = 0

  end subroutine off_open

  !> Close a off file 
  subroutine off_close(luoff)
    use tools_io, only: ferror, faterr, getline_raw, string, fopen_write, fclose
    integer, intent(in) :: luoff

    integer :: lu
    character(len=:), allocatable :: file, line
    logical :: ok

    if (.not.isopen) &
       call ferror('off_close','error: graphics file is not open',faterr)

    ! get the temporary file name from the first line and open
    rewind(luoff)
    ok = getline_raw(luoff,file,.true.)
    lu = fopen_write(file)

    ! write the header
    write (lu,'("COFF")')
    write (lu,'(3(A,X))') string(nv), string(nf), string(nv+nf-2)

    ! transfer the vertices from the scratch file over to the off file
    do while(getline_raw(luoff,line))
       if (line(1:1) == "v") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! transfer the faces
    rewind(luoff)
    ok = getline_raw(luoff,file,.true.)
    do while(getline_raw(luoff,line))
       if (line(1:1) == "f") then
          write (lu,'(A)') line(2:)
       endif
    end do

    ! close both
    call fclose(lu)
    call fclose(luoff)
    isopen = .false.

    ! clear the database
    nv = 0
    nf = 0
    
  end subroutine off_close

  !> Write a ball to the off file
  subroutine off_ball(luoff,x,rgb,r)
    use tools_io, only: faterr, ferror, string
    integer, intent(in) :: luoff
    real*8, intent(in) :: x(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: xrgb(3)

    if (.not.isopen) &
       call ferror('off_ball','error: graphics file is not open',faterr)

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
       write (luoff,'("v",3(F20.12,X),3(A,X),"0.0")') x + r * vsph(:,i,lvl), &
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, nfsph(lvl)
       iff = nv + fsph(:,i,lvl) - 1
       write (luoff,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    nv = nv + nvsph(lvl)
    nf = nf + nfsph(lvl)

  end subroutine off_ball
  
  !> Write a polygon to the off file. The coordinates are in array x
  !> (Cartesian), and are assumed to be in order.
  subroutine off_polygon(luoff,x,rgb)
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: luoff
    real*8, intent(in) :: x(:,:)
    integer, intent(in) :: rgb(3)

    integer :: i, i1, n
    character(len=:), allocatable :: str, aux
    real*8 :: xcm(3)

    if (.not.isopen) &
       call ferror('off_polygon','error: graphics file is not open',faterr)

    ! center of the face
    n = size(x,2)
    xcm = 0d0
    do i = 1, n
       xcm = xcm + x(:,i)
    end do
    xcm = xcm / n

    ! write the face to the obj
    do i = 1, n
       write (luoff,'("v",3(F20.12,X),3(A,X),"0")') x(:,i), &
          string(rgb(1)), string(rgb(2)), string(rgb(3))
    end do
    write (luoff,'("v",3(F20.12,X),3(A,X),"0")') xcm, &
       string(rgb(1)), string(rgb(2)), string(rgb(3))

    do i = 1, n
       i1 = mod(i,n) + 1
       write (luoff,'("f3 ",3(A,X))') string(nv+i-1), string(nv+i1-1), string(nv+n)
    end do
    nv = nv + n + 1
    nf = nf + n

  end subroutine off_polygon

  !> Write a stick to the off file
  subroutine off_stick(luoff,x1,x2,rgb,r)
    use tools_io, only: faterr, ferror, string
    use tools_math, only: norm, cross
    integer, intent(in) :: luoff
    real*8, intent(in) :: x1(3), x2(3)
    integer, intent(in) :: rgb(3)
    real*8, intent(in) :: r

    integer :: i, lvl, iff(3)
    real*8 :: dist, xd(3), v1(3), v2(3)
    real*8 :: xrgb(3)

    if (.not.isopen) &
       call ferror('off_stick','error: graphics file is not open',faterr)

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
       write (luoff,'("v",3(F20.12,X),3(A,X),"0")') &
          x1+vcyl(1,i,lvl)*v1+vcyl(2,i,lvl)*v2+vcyl(3,i,lvl)*xd,&
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, nfcyl(lvl)
       iff = nv + fcyl(:,i,lvl) - 1
       write (luoff,'("f3 ",3(A,X))') string(iff(1)), string(iff(2)), string(iff(3))
    end do
    nv = nv + nvcyl(lvl)
    nf = nf + nfcyl(lvl)

  end subroutine off_stick
  
  !> Write a surface to the off file
  subroutine off_surf(luoff,srf,fsurf)
    use types, only: minisurf
    use tools_io, only: faterr, ferror, string
    use param, only: pi
    integer, intent(in) :: luoff
    type(minisurf), intent(in) :: srf
    real*8, intent(in), optional :: fsurf(:)

    integer, parameter :: rgb_default(3) = (/128,128,128/)

    integer :: i, j, rgb(3)
    real*8 :: x(3), z
    real*8 :: xrgb(3), maxf, minf

    if (srf%init <= 1) &
       call ferror ('off_surf','No face information in minisurf',faterr)
    if (.not.isopen) &
       call ferror('off_surf','error: graphics file is not open',faterr)

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
       write (luoff,'("v",3(F20.12,X),3(A,X),"0")') x,&
          string(xrgb(1),"g"), string(xrgb(2),"g"), string(xrgb(3),"g")
    end do
    do i = 1, srf%nf
       write (luoff,'("f",999(A,X))') string(srf%f(i)%nv),&
          (string(nv+srf%f(i)%v(j)-1),j=1,srf%f(i)%nv)
    end do
    nv = nv + srf%nv
    nf = nf + srf%nf

  end subroutine off_surf

end module graphics
