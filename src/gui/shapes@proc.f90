! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! OpenGL buffers for simple shapes
submodule (shapes) proc
  use iso_c_binding
  implicit none

  ! number of vertices and faces
  integer(c_int), parameter :: sphnve(nmaxsph) = (/12/)
  integer(c_int), parameter :: sphnel(nmaxsph) = (/20/)
  integer(c_int), parameter :: sphneladd(0:nmaxsph) = (/0,20/)
  ! integer(c_int), parameter :: sphnve(nmaxsph) = (/12,42,162,642/)
  ! integer(c_int), parameter :: sphnel(nmaxsph) = (/20,80,320,1280/)
  ! integer(c_int), parameter :: sphneladd(0:nmaxsph) = (/0,20,100,420,1700/)

  ! static float cylr = cosf(30.f * PI / 180.f);
  ! static float cylnorm = sqrtf(0.5f + cylr * cylr);
  ! static float cylrn = cylr / cylnorm;
  ! static float halfn = 0.5f / cylnorm;
  ! static float halfn2 = 2 * halfn;

  ! // global variables
  ! GLuint sphVAO[nmaxsph];
  ! GLuint sphVBO;
  ! GLuint sphEBO[nmaxsph];
  ! GLuint cylVAO[nmaxcyl];
  ! GLuint cylVBO[nmaxcyl];
  ! GLuint cylEBO[nmaxcyl];
  ! const GLuint cylnve[nmaxcyl]      =    {14, 26, 38};
  ! const GLuint cylnveadd[nmaxcyl+1] = {0, 14, 40, 78};
  ! const GLuint cylnel[nmaxcyl]      =    {24, 48, 72};
  ! const GLuint cylneladd[nmaxcyl+1] = {0, 24, 72, 144};

  ! icosphere vertices (1:3) and normals (4:6)
  real(c_float), allocatable, target :: sphv(:,:)

  ! icosphere faces
  integer(c_int), allocatable, target :: sphi(:,:)

  ! // cylinder vertices (hexagonal prism)
  ! static GLfloat *cylv;
  ! static GLfloat cylv0[] = {
  !     0.0,  0.0, -0.5,     0.0,     0.0,   -1.0,
  !    cylr,  0.5, -0.5,   cylrn,   halfn, -halfn,
  !    cylr, -0.5, -0.5,   cylrn,  -halfn, -halfn,
  !     0.0, -1.0, -0.5,     0.0, -halfn2, -halfn,
  !   -cylr, -0.5, -0.5,  -cylrn,  -halfn, -halfn,
  !   -cylr,  0.5, -0.5,  -cylrn,   halfn, -halfn,
  !     0.0,  1.0, -0.5,     0.0,  halfn2, -halfn,
  !     0.0,  0.0,  0.5,     0.0,     0.0,    1.0,
  !    cylr,  0.5,  0.5,   cylrn,   halfn,  halfn,
  !    cylr, -0.5,  0.5,   cylrn,  -halfn,  halfn,
  !     0.0, -1.0,  0.5,     0.0, -halfn2,  halfn,
  !   -cylr, -0.5,  0.5,  -cylrn,  -halfn,  halfn,
  !   -cylr,  0.5,  0.5,  -cylrn,   halfn,  halfn,
  !     0.0,  1.0,  0.5,     0.0,  halfn2,  halfn,
  ! };

  ! // hexagonal prism faces
  ! static GLuint *cyli;
  ! static GLuint cyli0[] = {
  !   0,  1,   2,
  !   0,  2,   3,
  !   0,  3,   4,
  !   0,  4,   5,
  !   0,  5,   6,
  !   0,  6,   1,
  !   7,  9,   8,
  !   7, 10,   9,
  !   7, 11,  10,
  !   7, 12,  11,
  !   7, 13,  12,
  !   7,  8,  13,
  !   1,  8,   2,
  !   2,  8,   9,
  !   2,  9,   3,
  !   3,  9,  10,
  !   3, 10,   4,
  !   4, 10,  11,
  !   4, 11,   5,
  !   5, 11,  12,
  !   5, 12,   6,
  !   6, 12,  13,
  !   6, 13,   1,
  !   1, 13,   8,
  ! };

  ! test vertices
  real(c_float), allocatable, target :: testv(:,:)

contains

  !> Create and initialize the buffers for the basic shapes
  module subroutine shapes_init()
    use interfaces_opengl3
    use tools_math, only: cross

    real(c_float) :: tau, rad0, xico, zico, zero
    integer :: i
    integer(c_int) :: c_int_
    real(c_float) :: c_float_

    ! allocate icosphere vertex and face arrays
    if (allocated(sphv)) deallocate(sphv)
    if (allocated(sphi)) deallocate(sphi)
    allocate(sphv(6,sphnve(nmaxsph)))
    allocate(sphi(3,sphneladd(nmaxsph)))

    ! initialize icosahedron
    zero = 0._c_float
    tau = (1._c_float + sqrt(5._c_float))/2._c_float
    rad0 = sqrt(3._c_float - tau)
    xico = (tau - 1._c_float) / rad0
    zico = 1._c_float / rad0

    ! vertices and normals
    sphv(:,1 ) = (/-xico,  zero,  zico, -xico,  zero,  zico/)
    sphv(:,2 ) = (/ xico,  zero,  zico,  xico,  zero,  zico/)
    sphv(:,3 ) = (/-xico,  zero, -zico, -xico,  zero, -zico/)
    sphv(:,4 ) = (/ xico,  zero, -zico,  xico,  zero, -zico/)
    sphv(:,5 ) = (/ zero,  zico,  xico,  zero,  zico,  xico/)
    sphv(:,6 ) = (/ zero,  zico, -xico,  zero,  zico, -xico/)
    sphv(:,7 ) = (/ zero, -zico,  xico,  zero, -zico,  xico/)
    sphv(:,8 ) = (/ zero, -zico, -xico,  zero, -zico, -xico/)
    sphv(:,9 ) = (/ zico,  xico,  zero,  zico,  xico,  zero/)
    sphv(:,10) = (/-zico,  xico,  zero, -zico,  xico,  zero/)
    sphv(:,11) = (/ zico, -xico,  zero,  zico, -xico,  zero/)
    sphv(:,12) = (/-zico, -xico,  zero, -zico, -xico,  zero/)

    ! faces, counter-clock-wise
    sphi(:,1 ) = (/ 0_c_int,  1_c_int,  4_c_int/)
    sphi(:,2 ) = (/ 0_c_int,  4_c_int,  9_c_int/)
    sphi(:,3 ) = (/ 9_c_int,  4_c_int,  5_c_int/)
    sphi(:,4 ) = (/ 4_c_int,  8_c_int,  5_c_int/)
    sphi(:,5 ) = (/ 4_c_int,  1_c_int,  8_c_int/)
    sphi(:,6 ) = (/ 8_c_int,  1_c_int, 10_c_int/)
    sphi(:,7 ) = (/ 8_c_int, 10_c_int,  3_c_int/)
    sphi(:,8 ) = (/ 5_c_int,  8_c_int,  3_c_int/)
    sphi(:,9 ) = (/ 5_c_int,  3_c_int,  2_c_int/)
    sphi(:,10) = (/ 2_c_int,  3_c_int,  7_c_int/)
    sphi(:,11) = (/ 7_c_int,  3_c_int, 10_c_int/)
    sphi(:,12) = (/ 7_c_int, 10_c_int,  6_c_int/)
    sphi(:,13) = (/ 7_c_int,  6_c_int, 11_c_int/)
    sphi(:,14) = (/11_c_int,  6_c_int,  0_c_int/)
    sphi(:,15) = (/ 0_c_int,  6_c_int,  1_c_int/)
    sphi(:,16) = (/ 6_c_int, 10_c_int,  1_c_int/)
    sphi(:,17) = (/ 9_c_int, 11_c_int,  0_c_int/)
    sphi(:,18) = (/ 9_c_int,  2_c_int, 11_c_int/)
    sphi(:,19) = (/ 9_c_int,  5_c_int,  2_c_int/)
    sphi(:,20) = (/ 7_c_int, 11_c_int,  2_c_int/)

    ! // recursive subdivision for the other icospheres

    ! ! build the buffers for the spheres
    ! call glGenVertexArrays(nmaxsph, c_loc(sphVAO))
    ! call glGenBuffers(1, c_loc(sphVBO))
    ! call glGenBuffers(nmaxsph, c_loc(sphEBO))
    ! call glBindBuffer(GL_ARRAY_BUFFER, sphVBO)
    ! call glBufferData(GL_ARRAY_BUFFER, 6*sphnve(nmaxsph)*c_sizeof(c_float_), c_loc(sphv), GL_STATIC_DRAW)

    ! do i = 1, nmaxsph
    !    call glBindVertexArray(sphVAO(i))
    !    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphEBO(i))
    !    call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*sphnel(i)*c_sizeof(c_int_), c_loc(sphi(1,sphneladd(i-1)+1)), GL_STATIC_DRAW)

    !    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
    !       0_c_int)
    !    call glEnableVertexAttribArray(0)
    !    call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
    !       int(3 * c_sizeof(c_float_),c_int))
    !    call glEnableVertexAttribArray(1)
    ! end do

    ! call glBindBuffer(GL_ARRAY_BUFFER, 0)
    ! call glBindVertexArray(0)
    ! call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! // allocate space for the cylinders

    ! test object
    if (allocated(testv)) deallocate(testv)
    allocate(testv(3,3))
    testv(:,1) = (/-0.5_c_float, -0.5_c_float, 0.0_c_float/)
    testv(:,2) = (/ 0.5_c_float, -0.5_c_float, 0.0_c_float/)
    testv(:,3) = (/ 0.0_c_float,  0.5_c_float, 0.0_c_float/)

    call glGenVertexArrays(1, c_loc(testVAO))
    call glGenBuffers(1, c_loc(testVBO))
    call glBindVertexArray(testVAO)
    call glBindBuffer(GL_ARRAY_BUFFER, testVBO)
    call glBufferData(GL_ARRAY_BUFFER, 9*c_sizeof(c_float_), c_loc(testv), GL_STATIC_DRAW)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), &
       int(3*c_sizeof(c_float_),c_int), c_null_ptr)
    call glEnableVertexAttribArray(0)
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

  end subroutine shapes_init

  !> Terminate the shapes buffers
  module subroutine shapes_end()
    use interfaces_opengl3

    ! spheres
    call glDeleteVertexArrays(nmaxsph, c_loc(sphVAO))
    call glDeleteBuffers(1, c_loc(sphVBO))
    call glDeleteBuffers(nmaxsph, c_loc(sphEBO))
    if (allocated(sphv)) deallocate(sphv)
    if (allocated(sphi)) deallocate(sphi)

    ! glDeleteVertexArrays(nmaxcyl, cylVAO);
    ! glDeleteBuffers(nmaxcyl, cylVBO);
    ! glDeleteBuffers(nmaxcyl, cylEBO);
    ! delete [] sphi;
    ! delete [] sphv;
    ! delete [] cyli;
    ! delete [] cylv;

  end subroutine shapes_end

end submodule proc
