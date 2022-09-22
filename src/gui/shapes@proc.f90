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
    ! for (int i = 1; i < nmaxsph; i++){
    !   std::map<std::pair<int,int>,int> edgev = {};

    !   int n = sphnve[i-1] - 1;
    !   int nface = sphneladd[i]-1;
    !   for (int j = sphneladd[i-1]; j < sphneladd[i]; j++){
    !     int k1 = sphi[3*j+0];
    !     int k2 = sphi[3*j+1];
    !     int k3 = sphi[3*j+2];
    !     int nk1, nk2, nk3;

    !     // create the new vertices
    !     if (edgev[std::make_pair(k1,k2)]){
    !       nk1 = edgev[std::make_pair(k1,k2)];
    !     } else {
    !       edgev[std::make_pair(k1,k2)] = edgev[std::make_pair(k2,k1)] = nk1 = ++n;
    !       sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k1+0] + sphv[6*k2+0]);
    !       sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k1+1] + sphv[6*k2+1]);
    !       sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k1+2] + sphv[6*k2+2]);
    !     }
    !     if (edgev[std::make_pair(k1,k3)]){
    !       nk2 = edgev[std::make_pair(k1,k3)];
    !     } else {
    !       edgev[std::make_pair(k1,k3)] = edgev[std::make_pair(k3,k1)] = nk2 = ++n;
    !       sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k1+0] + sphv[6*k3+0]);
    !       sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k1+1] + sphv[6*k3+1]);
    !       sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k1+2] + sphv[6*k3+2]);
    !     }
    !     if (edgev[std::make_pair(k2,k3)]){
    !       nk3 = edgev[std::make_pair(k2,k3)];
    !     } else {
    !       edgev[std::make_pair(k2,k3)] = edgev[std::make_pair(k3,k2)] = nk3 = ++n;
    !       sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k2+0] + sphv[6*k3+0]);
    !       sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k2+1] + sphv[6*k3+1]);
    !       sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k2+2] + sphv[6*k3+2]);
    !     }

    !     // create the new faces
    !     nface++;
    !     sphi[3*nface+0] = k1;  sphi[3*nface+1] = nk1; sphi[3*nface+2] = nk2;
    !     nface++;
    !     sphi[3*nface+0] = nk1; sphi[3*nface+1] = nk3; sphi[3*nface+2] = nk2;
    !     nface++;
    !     sphi[3*nface+0] = nk1; sphi[3*nface+1] = k2;  sphi[3*nface+2] = nk3;
    !     nface++;
    !     sphi[3*nface+0] = nk2; sphi[3*nface+1] = nk3; sphi[3*nface+2] = k3;
    !   }

    !   // normalize the vertices
    !   for (int j = sphnve[i-1]; j < sphnve[i]; j++){
    !     float anorm = sqrtf(sphv[6*j+0]*sphv[6*j+0]+sphv[6*j+1]*sphv[6*j+1]+sphv[6*j+2]*sphv[6*j+2]);
    !     sphv[6*j+0] = sphv[6*j+3] = sphv[6*j+0] / anorm;
    !     sphv[6*j+1] = sphv[6*j+4] = sphv[6*j+1] / anorm;
    !     sphv[6*j+2] = sphv[6*j+5] = sphv[6*j+2] / anorm;
    !   }
    ! }

    ! build the buffers for the spheres
    call glGenVertexArrays(nmaxsph, c_loc(sphVAO))
    call glGenBuffers(1, c_loc(sphVBO))
    call glGenBuffers(nmaxsph, c_loc(sphEBO))
    call glBindBuffer(GL_ARRAY_BUFFER, sphVBO)
    call glBufferData(GL_ARRAY_BUFFER, 6*sphnve(nmaxsph)*c_sizeof(c_float_), c_loc(sphv), GL_STATIC_DRAW)

    do i = 1, nmaxsph
       call glBindVertexArray(sphVAO(i))
       call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphEBO(i))
       call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*sphnel(i)*c_sizeof(c_int_), c_loc(sphi(1,sphneladd(i-1)+1)), GL_STATIC_DRAW)

       call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
          0_c_int)
       call glEnableVertexAttribArray(0)
       call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
          int(3 * c_sizeof(c_float_),c_int))
       call glEnableVertexAttribArray(1)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! // allocate space for the cylinders
    ! cylv = new GLfloat [6*cylnveadd[nmaxcyl]];
    ! cyli = new GLuint [3*cylneladd[nmaxcyl]];
    ! memcpy(cylv,cylv0,sizeof(cylv0));
    ! memcpy(cyli,cyli0,sizeof(cyli0));

    ! // create the vertices, normals, and faces
    ! for (int ncyl = 0; ncyl < nmaxcyl; ncyl++){
    !   int n = cylnveadd[ncyl]-1;
    !   int npt = (cylnve[ncyl]-2)/2;

    !   n++;
    !   cylv[6*n+0] = cylv[6*n+1] = 0.f; cylv[6*n+2] = -0.5f;
    !   cylv[6*n+3] = cylv[6*n+4] = 0.f; cylv[6*n+5] = -1.0f;
    !   n++;
    !   cylv[6*n+0] = cylv[6*n+1] = 0.f; cylv[6*n+2] = 0.5f;
    !   cylv[6*n+3] = cylv[6*n+4] = 0.f; cylv[6*n+5] = 1.0f;


    !   float half = 0.5f;
    !   for (int j=0;j<2;j++){
    !     half = -half;
    !     for (int i=0;i<npt;i++){
    !       float angle = ((float) i) / ((float) npt) * 2.f * PI;
    !       float ca = cos(angle), sa = sin(angle);
    !       n++;
    !       cylv[6*n+0] = 0.5f * ca; cylv[6*n+1] = 0.5f * sa; cylv[6*n+2] = half;
    !       cylv[6*n+3] = ca; cylv[6*n+4] = sa; cylv[6*n+5] = 0.f;
    !     }
    !   }

    !   int nface = cylneladd[ncyl]-1;
    !   int shift1 = 2;
    !   int shift2 = 2+npt;
    !   for (int i=0; i<npt; i++){
    !     nface++;
    !     cyli[3*nface+0] = 0;
    !     cyli[3*nface+1] = (i+1)%npt+shift1;
    !     cyli[3*nface+2] = i+shift1;
    !   }
    !   for (int i=0; i<npt; i++){
    !     nface++;
    !     cyli[3*nface+0] = 1;
    !     cyli[3*nface+1] = i+shift2;
    !     cyli[3*nface+2] = (i+1)%npt+shift2;
    !   }
    !   for (int i=0; i<npt; i++){
    !     nface++;
    !     cyli[3*nface+0] = i+shift1;
    !     cyli[3*nface+1] = (i+1)%npt+shift1;
    !     cyli[3*nface+2] = (i+1)%npt+shift2;
    !   }
    !   for (int i=0; i<npt; i++){
    !     nface++;
    !     cyli[3*nface+0] = i+shift1;
    !     cyli[3*nface+1] = (i+1)%npt+shift2;
    !     cyli[3*nface+2] = i+shift2;
    !   }
    ! }


    ! // build the buffers for the cylinders
    ! glGenVertexArrays(nmaxcyl, cylVAO);
    ! glGenBuffers(nmaxcyl, cylVBO);
    ! glGenBuffers(nmaxcyl, cylEBO);

    ! for (int i=0;i<nmaxcyl;i++){
    !   glBindBuffer(GL_ARRAY_BUFFER, cylVBO[i]);
    !   glBufferData(GL_ARRAY_BUFFER, 6 * cylnve[i] * sizeof(GLfloat), &(cylv[6*cylnveadd[i]]), GL_STATIC_DRAW);

    !   glBindVertexArray(cylVAO[i]);

    !   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylEBO[i]);
    !   glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * cylnel[i] * sizeof(GLuint), &(cyli[3*cylneladd[i]]), GL_STATIC_DRAW);

    !   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
    !   glEnableVertexAttribArray(0);
    !   glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
    !   glEnableVertexAttribArray(1);
    ! }

    ! glBindBuffer(GL_ARRAY_BUFFER, 0);
    ! glBindVertexArray(0);
    ! glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  end subroutine shapes_init

  !> Terminate the shapes buffers
  module subroutine shapes_end()
    use interfaces_opengl3

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
