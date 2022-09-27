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

  ! icospheres
  real(c_float), allocatable, target :: sphv(:,:) ! vertices (1:3) and normals (4:6)
  integer(c_int), allocatable, target :: sphi(:,:) ! faces

  ! cylinders
  real(c_float), allocatable, target :: cylv(:,:) ! vertices (1:3) and normals (4:6)
  integer(c_int), allocatable, target :: cyli(:,:) ! faces

  ! test vertices
  real(c_float), allocatable, target :: testv(:,:)

contains

  !> Create and initialize the buffers for the basic shapes
  module subroutine shapes_init()
    use interfaces_opengl3
    use tools_math, only: cross
    use hashmod, only: hash
    use param, only: pi
    real(c_float) :: tau, rad0, xico, zico, zero, half, one
    real(c_float) :: cylr, pic, cylnorm, cylrn, halfn, halfn2
    integer :: i, j, k(3), nk1, nk2, nk3, n, nface
    integer :: nk1, nk2, nk3, maxsphnve
    integer(c_int) :: c_int_
    real(c_float) :: c_float_
    type(c_ptr) :: c_ptr_
    type(hash) :: ipair
    character(len=:), allocatable :: idx

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

    ! sphere faces, counter-clock-wise
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

    ! vertices and edges of the rest of the spheres
    maxsphnve = maxval(sphnve)
    do i = 1, nmaxsph-1
       n = sphnve(i)
       nface = sphneladd(i)
       do j = sphneladd(i-1)+1, sphneladd(i)
          k = sphi(:,j) + 1

          ! create the new vertices
          call ipack(k(1),k(2))
          if (ipair%iskey(idx)) then
             nk1 = ipair%get(idx,nk1)
          else
             n = n + 1
             nk1 = n
             call ipair%put(idx,n)
             sphv(1:3,n) = 0.5_c_float * (sphv(1:3,k(1)) + sphv(1:3,k(2)))
             sphv(1:3,n) = sphv(1:3,n) / norm2(sphv(1:3,n))
             sphv(4:6,n) = sphv(1:3,n)
          end if
          call ipack(k(1),k(3))
          if (ipair%iskey(idx)) then
             nk2 = ipair%get(idx,nk2)
          else
             n = n + 1
             nk2 = n
             call ipair%put(idx,n)
             sphv(1:3,n) = 0.5_c_float * (sphv(1:3,k(1)) + sphv(1:3,k(3)))
             sphv(1:3,n) = sphv(1:3,n) / norm2(sphv(1:3,n))
             sphv(4:6,n) = sphv(1:3,n)
          end if
          call ipack(k(2),k(3))
          if (ipair%iskey(idx)) then
             nk3 = ipair%get(idx,nk3)
          else
             n = n + 1
             nk3 = n
             call ipair%put(idx,n)
             sphv(1:3,n) = 0.5_c_float * (sphv(1:3,k(2)) + sphv(1:3,k(3)))
             sphv(1:3,n) = sphv(1:3,n) / norm2(sphv(1:3,n))
             sphv(4:6,n) = sphv(1:3,n)
          end if

          ! create the new faces
          nface = nface + 1
          sphi(:,nface) = (/k(1),nk1,nk2/) - 1
          nface = nface + 1
          sphi(:,nface) = (/nk1,nk3,nk2/) - 1
          nface = nface + 1
          sphi(:,nface) = (/nk1,k(2),nk3/) - 1
          nface = nface + 1
          sphi(:,nface) = (/nk2,nk3,k(3)/) - 1
       end do
    end do

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
          c_null_ptr)
       call glEnableVertexAttribArray(0)
       call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
          transfer(3_c_int * c_sizeof(c_float_),c_ptr_))
       call glEnableVertexAttribArray(1)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! allocate cylinder vertex and face arrays
    if (allocated(cylv)) deallocate(cylv)
    if (allocated(cyli)) deallocate(cyli)
    allocate(cylv(6,cylnveadd(nmaxcyl)))
    allocate(cyli(3,cylneladd(nmaxcyl)))

    ! initialize cylinders
    pic = real(pi,c_float)
    half = 0.5_c_float
    one = 1.0_c_float
    cylr = cos(30._c_float * pic / 180._c_float)
    cylnorm = sqrt(0.5_c_float + cylr * cylr)
    cylrn = cylr / cylnorm
    halfn = 0.5_c_float / cylnorm
    halfn2 = 2 * halfn

    ! vertices and normals
    cylv(:,1 ) = (/ zero, zero,-half,    zero,    zero,   -one/)
    cylv(:,2 ) = (/ cylr, half,-half,   cylrn,   halfn, -halfn/)
    cylv(:,3 ) = (/ cylr,-half,-half,   cylrn,  -halfn, -halfn/)
    cylv(:,4 ) = (/ zero, -one,-half,    zero, -halfn2, -halfn/)
    cylv(:,5 ) = (/-cylr,-half,-half,  -cylrn,  -halfn, -halfn/)
    cylv(:,6 ) = (/-cylr, half,-half,  -cylrn,   halfn, -halfn/)
    cylv(:,7 ) = (/ zero,  one,-half,    zero,  halfn2, -halfn/)
    cylv(:,8 ) = (/ zero, zero, half,    zero,    zero,    one/)
    cylv(:,9 ) = (/ cylr, half, half,   cylrn,   halfn,  halfn/)
    cylv(:,10) = (/ cylr,-half, half,   cylrn,  -halfn,  halfn/)
    cylv(:,11) = (/ zero, -one, half,    zero, -halfn2,  halfn/)
    cylv(:,12) = (/-cylr,-half, half,  -cylrn,  -halfn,  halfn/)
    cylv(:,13) = (/-cylr, half, half,  -cylrn,   halfn,  halfn/)
    cylv(:,14) = (/ zero,  one, half,    zero,  halfn2,  halfn/)

    ! cylinder faces, counter-clock-wise
    cyli(:,1 ) = (/0_c_int,  1_c_int,  2_c_int/)
    cyli(:,2 ) = (/0_c_int,  2_c_int,  3_c_int/)
    cyli(:,3 ) = (/0_c_int,  3_c_int,  4_c_int/)
    cyli(:,4 ) = (/0_c_int,  4_c_int,  5_c_int/)
    cyli(:,5 ) = (/0_c_int,  5_c_int,  6_c_int/)
    cyli(:,6 ) = (/0_c_int,  6_c_int,  1_c_int/)
    cyli(:,7 ) = (/7_c_int,  9_c_int,  8_c_int/)
    cyli(:,8 ) = (/7_c_int, 10_c_int,  9_c_int/)
    cyli(:,9 ) = (/7_c_int, 11_c_int, 10_c_int/)
    cyli(:,10) = (/7_c_int, 12_c_int, 11_c_int/)
    cyli(:,11) = (/7_c_int, 13_c_int, 12_c_int/)
    cyli(:,12) = (/7_c_int,  8_c_int, 13_c_int/)
    cyli(:,13) = (/1_c_int,  8_c_int,  2_c_int/)
    cyli(:,14) = (/2_c_int,  8_c_int,  9_c_int/)
    cyli(:,15) = (/2_c_int,  9_c_int,  3_c_int/)
    cyli(:,16) = (/3_c_int,  9_c_int, 10_c_int/)
    cyli(:,17) = (/3_c_int, 10_c_int,  4_c_int/)
    cyli(:,18) = (/4_c_int, 10_c_int, 11_c_int/)
    cyli(:,19) = (/4_c_int, 11_c_int,  5_c_int/)
    cyli(:,20) = (/5_c_int, 11_c_int, 12_c_int/)
    cyli(:,21) = (/5_c_int, 12_c_int,  6_c_int/)
    cyli(:,22) = (/6_c_int, 12_c_int, 13_c_int/)
    cyli(:,23) = (/6_c_int, 13_c_int,  1_c_int/)
    cyli(:,24) = (/1_c_int, 13_c_int,  8_c_int/)

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

    ! build the buffers for the spheres
    call glGenVertexArrays(nmaxcyl, c_loc(cylVAO))
    call glGenBuffers(nmaxcyl, c_loc(cylVBO))
    call glGenBuffers(nmaxcyl, c_loc(cylEBO))

    do i = 1, nmaxcyl
       call glBindBuffer(GL_ARRAY_BUFFER, cylVBO(i))
       call glBufferData(GL_ARRAY_BUFFER, 6*cylnve(i)*c_sizeof(c_float_), c_loc(cylv(1,cylnveadd(i-1)+1)), GL_STATIC_DRAW)
       call glBindVertexArray(cylVAO(i))
       call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylEBO(i))
       call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*cylnel(i)*c_sizeof(c_int_), c_loc(cyli(1,cylneladd(i-1)+1)), GL_STATIC_DRAW)

       call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
          c_null_ptr)
       call glEnableVertexAttribArray(0)
       call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
          transfer(3_c_int * c_sizeof(c_float_),c_ptr_))
       call glEnableVertexAttribArray(1)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! test object
    if (allocated(testv)) deallocate(testv)
    allocate(testv(6,3))
    testv(:,1) = (/ -half, -half, zero, half, half, half/)
    testv(:,2) = (/  half, -half, zero, half, half, half/)
    testv(:,3) = (/  zero,  half, zero, half, half, half/)

    call glGenVertexArrays(1, c_loc(testVAO))
    call glGenBuffers(1, c_loc(testVBO))
    call glBindVertexArray(testVAO)
    call glBindBuffer(GL_ARRAY_BUFFER, testVBO)
    call glBufferData(GL_ARRAY_BUFFER, 18*c_sizeof(c_float_), c_loc(testv), GL_STATIC_DRAW)
    ! call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int), c_null_ptr)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int), &
       transfer(0_c_int * c_sizeof(c_float_),c_ptr_))
    call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(6*c_sizeof(c_float_),c_int),&
       transfer(3_c_int * c_sizeof(c_float_),c_ptr_))

    call glEnableVertexAttribArray(0)
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

  contains
    subroutine ipack(i,j)
      use tools_io, only: string
      integer, intent(in) :: i, j

      if (i < j) then
         idx = string(i) // "_" // string(j)
      else
         idx = string(j) // "_" // string(i)
      end if

    end subroutine ipack

  end subroutine shapes_init

  !> Terminate the shapes buffers
  module subroutine shapes_end()
    use interfaces_opengl3

    ! test object
    call glDeleteVertexArrays(1, c_loc(testVAO))
    call glDeleteBuffers(1, c_loc(testVBO))
    if (allocated(testv)) deallocate(testv)

    ! spheres
    call glDeleteVertexArrays(nmaxsph, c_loc(sphVAO))
    call glDeleteBuffers(1, c_loc(sphVBO))
    call glDeleteBuffers(nmaxsph, c_loc(sphEBO))
    if (allocated(sphv)) deallocate(sphv)
    if (allocated(sphi)) deallocate(sphi)

  end subroutine shapes_end

end submodule proc
