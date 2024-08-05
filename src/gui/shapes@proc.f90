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

contains

  !> Create and initialize the buffers for the basic shapes
  module subroutine shapes_init()
    use interfaces_opengl3
    use tools_math, only: cross
    use hashmod, only: hash
    use param, only: pi
    real(c_float) :: tau, rad0, xico, zico, zero, half, one, angle
    real(c_float) :: cylr, pic, cylnorm, cylrn, halfn, halfn2, ca, sa
    integer :: i, l, j, k(3), nk1, nk2, nk3, n, nface, npt
    integer(c_int) :: shift1, shift2
    integer(c_int) :: c_int_
    real(c_float) :: c_float_
    type(c_ptr) :: c_ptr_
    type(hash) :: ipair
    character(len=:), allocatable :: idx

    ! allocate icosphere vertex and face arrays
    if (allocated(sphv)) deallocate(sphv)
    if (allocated(sphi)) deallocate(sphi)
    allocate(sphv(3,sphnve(nmaxsph)))
    allocate(sphi(3,sphneladd(nmaxsph)))

    ! initialize icosahedron
    zero = 0._c_float
    tau = (1._c_float + sqrt(5._c_float))/2._c_float
    rad0 = sqrt(3._c_float - tau)
    xico = (tau - 1._c_float) / rad0
    zico = 1._c_float / rad0

    ! vertices
    sphv(:,1 ) = (/-xico,  zero,  zico/)
    sphv(:,2 ) = (/ xico,  zero,  zico/)
    sphv(:,3 ) = (/-xico,  zero, -zico/)
    sphv(:,4 ) = (/ xico,  zero, -zico/)
    sphv(:,5 ) = (/ zero,  zico,  xico/)
    sphv(:,6 ) = (/ zero,  zico, -xico/)
    sphv(:,7 ) = (/ zero, -zico,  xico/)
    sphv(:,8 ) = (/ zero, -zico, -xico/)
    sphv(:,9 ) = (/ zico,  xico,  zero/)
    sphv(:,10) = (/-zico,  xico,  zero/)
    sphv(:,11) = (/ zico, -xico,  zero/)
    sphv(:,12) = (/-zico, -xico,  zero/)

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
    call glBufferData(GL_ARRAY_BUFFER, 3*sphnve(nmaxsph)*c_sizeof(c_float_), c_loc(sphv), GL_STATIC_DRAW)

    do i = 1, nmaxsph
       call glBindVertexArray(sphVAO(i))
       call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphEBO(i))
       call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*sphnel(i)*c_sizeof(c_int_), c_loc(sphi(1,sphneladd(i-1)+1)), GL_STATIC_DRAW)

       call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int),&
          c_null_ptr)
       call glEnableVertexAttribArray(0)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! allocate cylinder vertex and face arrays
    if (allocated(cylv)) deallocate(cylv)
    if (allocated(cyli)) deallocate(cyli)
    allocate(cylv(3,cylnveadd(nmaxcyl)))
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

    ! vertices and edges of the rest of the cylinders
    do i = 1, nmaxcyl
       n = cylnveadd(i-1)
       npt = (cylnve(i)-2)/2

       ! end caps
       n = n + 1
       cylv(1:2,n) = 0._c_float
       cylv(3,n) = -0.5_c_float
       n = n + 1
       cylv(1:2,n) = 0._c_float
       cylv(3,n) = 0.5_c_float

       half = 0.5_c_float
       do j = 1, 2
          half = -half
          do l = 1, npt
             angle = real(l-1,c_float) / real(npt,c_float) * 2._c_float * pic
             ca = cos(angle)
             sa = sin(angle)
             n = n + 1
             cylv(1,n) = 0.5_c_float * ca
             cylv(2,n) = 0.5_c_float * sa
             cylv(3,n) = half
          end do
       end do

       nface = cylneladd(i-1)
       shift1 = 2
       shift2 = 2+npt
       do j = 1, npt
          nface = nface + 1
          cyli(:,nface) = (/0, mod(j,npt)+shift1, j+shift1-1/)
       end do
       do j = 1, npt
          nface = nface + 1
          cyli(:,nface) = (/1, j+shift2-1, mod(j,npt)+shift2/)
       end do
       do j = 1, npt
          nface = nface + 1
          cyli(:,nface) = (/j+shift1-1, mod(j,npt)+shift1, mod(j,npt)+shift2/)
       end do
       do j = 1, npt
          nface = nface + 1
          cyli(:,nface) = (/j+shift1-1, mod(j,npt)+shift2, j+shift2-1/)
       end do
    end do

    ! build the buffers for the cylinder
    call glGenVertexArrays(nmaxcyl, c_loc(cylVAO))
    call glGenBuffers(nmaxcyl, c_loc(cylVBO))
    call glGenBuffers(nmaxcyl, c_loc(cylEBO))

    do i = 1, nmaxcyl
       call glBindBuffer(GL_ARRAY_BUFFER, cylVBO(i))
       call glBufferData(GL_ARRAY_BUFFER, 3*cylnve(i)*c_sizeof(c_float_), c_loc(cylv(1,cylnveadd(i-1)+1)), GL_STATIC_DRAW)
       call glBindVertexArray(cylVAO(i))
       call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylEBO(i))
       call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*cylnel(i)*c_sizeof(c_int_), c_loc(cyli(1,cylneladd(i-1)+1)), GL_STATIC_DRAW)

       call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int),&
          c_null_ptr)
       call glEnableVertexAttribArray(0)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! build the buffers for the text (direct)
    call glGenVertexArrays(1, c_loc(textVAO))
    call glGenBuffers(1, c_loc(textVBO))
    call glBindBuffer(GL_ARRAY_BUFFER, textVBO)
    call glBindVertexArray(textVAO)
    call glBufferData(GL_ARRAY_BUFFER, text_maxvert*4*c_sizeof(c_float_), c_null_ptr, GL_DYNAMIC_DRAW)
    call glEnableVertexAttribArray(0)
    call glVertexAttribPointer(0, 4, GL_FLOAT, int(GL_FALSE,c_signed_char), &
       int(4*c_sizeof(c_float_),c_int),c_null_ptr)
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

    ! build the buffers for the text (on scene)
    call glGenVertexArrays(1, c_loc(textVAOos))
    call glGenBuffers(1, c_loc(textVBOos))
    call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)
    call glBindVertexArray(textVAOos)
    call glBufferData(GL_ARRAY_BUFFER, text_maxvert*8*c_sizeof(c_float_), c_null_ptr, GL_DYNAMIC_DRAW)
    call glEnableVertexAttribArray(0)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(8*c_sizeof(c_float_),c_int),&
       c_null_ptr)
    call glEnableVertexAttribArray(1)
    call glVertexAttribPointer(1, 1, GL_FLOAT, int(GL_FALSE,c_signed_char), int(8*c_sizeof(c_float_),c_int),&
       transfer(3_c_int * c_sizeof(c_float_),c_ptr_))
    call glEnableVertexAttribArray(2)
    call glVertexAttribPointer(2, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(8*c_sizeof(c_float_),c_int),&
       transfer(4_c_int * c_sizeof(c_float_),c_ptr_))
    call glEnableVertexAttribArray(3)
    call glVertexAttribPointer(3, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(8*c_sizeof(c_float_),c_int),&
       transfer(6_c_int * c_sizeof(c_float_),c_ptr_))
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

    ! spheres
    call glDeleteVertexArrays(nmaxsph, c_loc(sphVAO))
    call glDeleteBuffers(1, c_loc(sphVBO))
    call glDeleteBuffers(nmaxsph, c_loc(sphEBO))
    if (allocated(sphv)) deallocate(sphv)
    if (allocated(sphi)) deallocate(sphi)

    ! cylinders
    call glDeleteVertexArrays(nmaxcyl, c_loc(cylVAO))
    call glDeleteBuffers(nmaxcyl, c_loc(cylVBO))
    call glDeleteBuffers(nmaxcyl, c_loc(cylEBO))
    if (allocated(cylv)) deallocate(cylv)
    if (allocated(cyli)) deallocate(cyli)

  end subroutine shapes_end

  ! !> Initialize images for the icons in the GUI.
  ! module subroutine icons_init()
  !   use interfaces_opengl3
  !   use interfaces_stb
  !   use tools_io, only: ferror, faterr
  !   use global, only: critic_home
  !   use param, only: dirsep

  !   character(kind=c_char,len=:), allocatable, target :: file
  !   type(c_ptr) :: image
  !   integer(c_int) :: idum, width, height

  !   integer(c_int) :: icon_side = 64

  !   ! load the icon
  !   file = trim(critic_home) // dirsep // "assets" // dirsep // "qe.png" // c_null_char
  !   image = stbi_load(c_loc(file), width, height, idum, 4)
  !   if (.not.c_associated(image)) &
  !      call ferror('gui_start','Could not find GUI assets: have you set CRITIC_HOME?',faterr)

  !   ! create the texture
  !   call glGenTextures(1, c_loc(treeicon1))
  !   call glBindTexture(GL_TEXTURE_2D, treeicon1)
  !   call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
  !   call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
  !   call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image)
  !   call glBindTexture(GL_TEXTURE_2D, 0)

  !   ! free
  !   call stbi_image_free(image)

  ! end subroutine icons_init

end submodule proc
