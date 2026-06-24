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

  ! cones (arrowheads)
  real(c_float), allocatable, target :: conv(:,:) ! vertices (1:3)
  integer(c_int), allocatable, target :: coni(:,:) ! faces

contains

  !> Create and initialize the buffers for the basic shapes
  module subroutine shapes_init()
    use interfaces_opengl3
    use param, only: pi
    real(c_float) :: angle, pic, ca, sa
    integer :: i, l, j, n, nface, npt
    integer(c_int) :: shift1
    integer(c_int) :: c_int_
    real(c_float) :: c_float_
    type(c_ptr) :: c_ptr_
    real(c_float), target :: quadv(3,4)
    integer(c_int), target :: quadi(3,2)
    real(c_float), target :: triv(3,3)
    integer(c_int), target :: trii(3,1)
    real(c_float), target :: impquadc(2,4)

    pic = real(pi,c_float)

    ! allocate cone vertex and face arrays
    if (allocated(conv)) deallocate(conv)
    if (allocated(coni)) deallocate(coni)
    allocate(conv(3,connveadd(nmaxcone)))
    allocate(coni(3,conneladd(nmaxcone)))

    ! vertices and faces of the cones (unit cone: base radius 0.5 at
    ! z=-0.5, apex at z=+0.5)
    do i = 1, nmaxcone
       n = connveadd(i-1)
       npt = connve(i) - 2

       ! base center and apex
       n = n + 1
       conv(1:2,n) = 0._c_float
       conv(3,n) = -0.5_c_float
       n = n + 1
       conv(1:2,n) = 0._c_float
       conv(3,n) = 0.5_c_float

       ! base ring
       do l = 1, npt
          angle = real(l-1,c_float) / real(npt,c_float) * 2._c_float * pic
          ca = cos(angle)
          sa = sin(angle)
          n = n + 1
          conv(1,n) = 0.5_c_float * ca
          conv(2,n) = 0.5_c_float * sa
          conv(3,n) = -0.5_c_float
       end do

       nface = conneladd(i-1)
       shift1 = 2
       ! base cap (center vertex 0 + base ring)
       do j = 1, npt
          nface = nface + 1
          coni(:,nface) = (/0, mod(j,npt)+shift1, j+shift1-1/)
       end do
       ! sides (apex vertex 1 + base ring)
       do j = 1, npt
          nface = nface + 1
          coni(:,nface) = (/1, j+shift1-1, mod(j,npt)+shift1/)
       end do
    end do

    ! build the buffers for the cones
    call glGenVertexArrays(nmaxcone, c_loc(coneVAO))
    call glGenBuffers(nmaxcone, c_loc(coneVBO))
    call glGenBuffers(nmaxcone, c_loc(coneEBO))

    do i = 1, nmaxcone
       call glBindBuffer(GL_ARRAY_BUFFER, coneVBO(i))
       call glBufferData(GL_ARRAY_BUFFER, 3*connve(i)*c_sizeof(c_float_), c_loc(conv(1,connveadd(i-1)+1)), GL_STATIC_DRAW)
       call glBindVertexArray(coneVAO(i))
       call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coneEBO(i))
       call glBufferData(GL_ELEMENT_ARRAY_BUFFER,3*connel(i)*c_sizeof(c_int_), c_loc(coni(1,conneladd(i-1)+1)), GL_STATIC_DRAW)

       call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int),&
          c_null_ptr)
       call glEnableVertexAttribArray(0)
    end do

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! build the buffers for the quad (unit rectangle in the z=0 plane, corners
    ! at (+/-1,+/-1,0); two triangles)
    quadv(:,1) = (/-1._c_float,-1._c_float,0._c_float/)
    quadv(:,2) = (/ 1._c_float,-1._c_float,0._c_float/)
    quadv(:,3) = (/ 1._c_float, 1._c_float,0._c_float/)
    quadv(:,4) = (/-1._c_float, 1._c_float,0._c_float/)
    quadi(:,1) = (/0_c_int,1_c_int,2_c_int/)
    quadi(:,2) = (/0_c_int,2_c_int,3_c_int/)
    call glGenVertexArrays(1, c_loc(quadVAO))
    call glGenBuffers(1, c_loc(quadVBO))
    call glGenBuffers(1, c_loc(quadEBO))
    call glBindBuffer(GL_ARRAY_BUFFER, quadVBO)
    call glBufferData(GL_ARRAY_BUFFER, 12*c_sizeof(c_float_), c_loc(quadv), GL_STATIC_DRAW)
    call glBindVertexArray(quadVAO)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, quadEBO)
    call glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6*c_sizeof(c_int_), c_loc(quadi), GL_STATIC_DRAW)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int),&
       c_null_ptr)
    call glEnableVertexAttribArray(0)
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    ! Build the buffers for the triangle (unit reference triangle in the z=0
    ! plane, corners at (0,0,0),(1,0,0),(0,1,0))
    triv(:,1) = (/0._c_float,0._c_float,0._c_float/)
    triv(:,2) = (/1._c_float,0._c_float,0._c_float/)
    triv(:,3) = (/0._c_float,1._c_float,0._c_float/)
    trii(:,1) = (/0_c_int,1_c_int,2_c_int/)
    call glGenVertexArrays(1, c_loc(triVAO))
    call glGenBuffers(1, c_loc(triVBO))
    call glGenBuffers(1, c_loc(triEBO))
    call glBindBuffer(GL_ARRAY_BUFFER, triVBO)
    call glBufferData(GL_ARRAY_BUFFER, 9*c_sizeof(c_float_), c_loc(triv), GL_STATIC_DRAW)
    call glBindVertexArray(triVAO)
    call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triEBO)
    call glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*c_sizeof(c_int_), c_loc(trii), GL_STATIC_DRAW)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int),&
       c_null_ptr)
    call glEnableVertexAttribArray(0)
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
    call glBufferData(GL_ARRAY_BUFFER, text_maxvert*10*c_sizeof(c_float_), c_null_ptr, GL_DYNAMIC_DRAW)
    call glEnableVertexAttribArray(0)
    call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(10*c_sizeof(c_float_),c_int),&
       c_null_ptr)
    call glEnableVertexAttribArray(1)
    call glVertexAttribPointer(1, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(10*c_sizeof(c_float_),c_int),&
       transfer(3_c_int * c_sizeof(c_float_),c_ptr_))
    call glEnableVertexAttribArray(2)
    call glVertexAttribPointer(2, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(10*c_sizeof(c_float_),c_int),&
       transfer(6_c_int * c_sizeof(c_float_),c_ptr_))
    call glEnableVertexAttribArray(3)
    call glVertexAttribPointer(3, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(10*c_sizeof(c_float_),c_int),&
       transfer(8_c_int * c_sizeof(c_float_),c_ptr_))
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

    ! impostor billboards: a single unit quad (triangle strip, corners in
    ! [-1,1]) is expanded per instance via instanced attributes.
    impquadc(:,1) = (/-1._c_float,-1._c_float/)
    impquadc(:,2) = (/ 1._c_float,-1._c_float/)
    impquadc(:,3) = (/-1._c_float, 1._c_float/)
    impquadc(:,4) = (/ 1._c_float, 1._c_float/)
    call glGenBuffers(1, c_loc(impQuadVBO))
    call glBindBuffer(GL_ARRAY_BUFFER, impQuadVBO)
    call glBufferData(GL_ARRAY_BUFFER, 2*4*c_sizeof(c_float_), c_loc(impquadc), GL_STATIC_DRAW)
    call glBindBuffer(GL_ARRAY_BUFFER, 0)

    ! the per-instance VAOs/VBOs (sphere/cylinder impostors and the plain meshes)
    ! are per-scene and created lazily by scene_glbuffers%init

  end subroutine shapes_init

  !> Terminate the shapes buffers
  module subroutine shapes_end()
    use interfaces_opengl3

    ! cones
    call glDeleteVertexArrays(nmaxcone, c_loc(coneVAO))
    call glDeleteBuffers(nmaxcone, c_loc(coneVBO))
    call glDeleteBuffers(nmaxcone, c_loc(coneEBO))
    if (allocated(conv)) deallocate(conv)
    if (allocated(coni)) deallocate(coni)

    ! quads
    call glDeleteVertexArrays(1, c_loc(quadVAO))
    call glDeleteBuffers(1, c_loc(quadVBO))
    call glDeleteBuffers(1, c_loc(quadEBO))

    ! triangles
    call glDeleteVertexArrays(1, c_loc(triVAO))
    call glDeleteBuffers(1, c_loc(triVBO))
    call glDeleteBuffers(1, c_loc(triEBO))

    ! shared impostor quad (the per-scene instance buffers are freed by
    ! scene_glbuffers%end)
    call glDeleteBuffers(1, c_loc(impQuadVBO))

  end subroutine shapes_end

  !> Create and configure this scene's instance VAOs/VBOs: sphere and cylinder
  !> impostors and the plain meshes (planes, polyhedra triangles, cones), each
  !> with a cached and (for the impostors and cones) a transient scratch buffer.
  !> The per-instance VAOs reference the shared template geometry (impQuadVBO and
  !> the cone/quad/tri meshes). Idempotent: a no-op if already initialized.
  module subroutine glbuffers_init(b)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout), target :: b
    real(c_float) :: c_float_
    type(c_ptr) :: c_ptr_

    if (b%isinit) return

    ! sphere impostor VAOs: cached main + transient scratch (same layout)
    call setup_sph_inst(b%sphinstVAO, b%sphinstVBO)
    call setup_sph_inst(b%sphinstVAOscr, b%sphinstVBOscr)

    ! cylinder impostor VAOs: cached main + transient scratch
    call setup_cyl_inst(b%cylinstVAO, b%cylinstVBO)
    call setup_cyl_inst(b%cylinstVAOscr, b%cylinstVBOscr)

    ! instanced plain-mesh VAOs (planes, polyhedra triangles, cones): the mesh
    ! vertex (location 0) plus a per-instance model matrix (1..4) + color (5).
    ! Cones always use the highest mesh resolution (counts are tiny); the cone
    ! scratch VAO serves the gizmo arrowheads.
    call setup_mesh_inst(b%planeinstVAO, b%planeinstVBO, quadVBO, quadEBO)
    call setup_mesh_inst(b%triinstVAO, b%triinstVBO, triVBO, triEBO)
    call setup_mesh_inst(b%coneinstVAO, b%coneinstVBO, coneVBO(nmaxcone), coneEBO(nmaxcone))
    call setup_mesh_inst(b%coneinstVAOscr, b%coneinstVBOscr, coneVBO(nmaxcone), coneEBO(nmaxcone))

    b%isinit = .true.

  contains
    !> Set up a sphere impostor VAO: shared quad corner (loc 0) + per-instance
    !> sphere attributes (loc 1..8) from the given instance buffer.
    subroutine setup_sph_inst(vao, vbo)
      integer(c_int), intent(inout), target :: vao, vbo
      integer(c_int) :: st

      call glGenVertexArrays(1, c_loc(vao))
      call glGenBuffers(1, c_loc(vbo))
      call glBindVertexArray(vao)
      call glBindBuffer(GL_ARRAY_BUFFER, impQuadVBO)
      call glEnableVertexAttribArray(0)
      call glVertexAttribPointer(0, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(2*c_sizeof(c_float_),c_int), c_null_ptr)
      call glBindBuffer(GL_ARRAY_BUFFER, vbo)
      st = int(sph_inst_nf * c_sizeof(c_float_),c_int)
      call inst_attrib(1, 3, st, 0)  ! a_center
      call inst_attrib(2, 1, st, 3)  ! a_radius
      call inst_attrib(3, 4, st, 4)  ! a_color
      call inst_attrib(4, 1, st, 8)  ! a_border
      call inst_attrib(5, 3, st, 9)  ! a_bordercolor
      call inst_attrib(6, 3, st, 12) ! a_xdelta_re
      call inst_attrib(7, 3, st, 15) ! a_xdelta_im
      call inst_attrib(8, 4, st, 18) ! a_idx
      call glBindBuffer(GL_ARRAY_BUFFER, 0)
      call glBindVertexArray(0)

    end subroutine setup_sph_inst

    !> Set up a cylinder impostor VAO: shared quad corner (loc 0) + per-instance
    !> cylinder attributes (loc 1..8) from the given instance buffer.
    subroutine setup_cyl_inst(vao, vbo)
      integer(c_int), intent(inout), target :: vao, vbo
      integer(c_int) :: st

      call glGenVertexArrays(1, c_loc(vao))
      call glGenBuffers(1, c_loc(vbo))
      call glBindVertexArray(vao)
      call glBindBuffer(GL_ARRAY_BUFFER, impQuadVBO)
      call glEnableVertexAttribArray(0)
      call glVertexAttribPointer(0, 2, GL_FLOAT, int(GL_FALSE,c_signed_char), int(2*c_sizeof(c_float_),c_int), c_null_ptr)
      call glBindBuffer(GL_ARRAY_BUFFER, vbo)
      st = int(cyl_inst_nf * c_sizeof(c_float_),c_int)
      call inst_attrib(1, 3, st, 0)  ! a_x1
      call inst_attrib(2, 3, st, 3)  ! a_x2
      call inst_attrib(3, 1, st, 6)  ! a_radius
      call inst_attrib(4, 4, st, 7)  ! a_color
      call inst_attrib(5, 1, st, 11) ! a_border
      call inst_attrib(6, 3, st, 12) ! a_bordercolor
      call inst_attrib(7, 1, st, 15) ! a_delta
      call inst_attrib(8, 3, st, 16) ! a_outward
      call glBindBuffer(GL_ARRAY_BUFFER, 0)
      call glBindVertexArray(0)

    end subroutine setup_cyl_inst
    !> Set up an instanced plain-mesh VAO: bind the mesh vertex buffer (loc 0)
    !> and element buffer, then a fresh per-instance buffer holding the model
    !> matrix columns (loc 1..4) and the color (loc 5).
    subroutine setup_mesh_inst(vao, vbo, meshvbo, meshebo)
      integer(c_int), intent(inout), target :: vao, vbo
      integer(c_int), intent(in) :: meshvbo, meshebo
      integer(c_int) :: st

      call glGenVertexArrays(1, c_loc(vao))
      call glGenBuffers(1, c_loc(vbo))
      call glBindVertexArray(vao)
      call glBindBuffer(GL_ARRAY_BUFFER, meshvbo)
      call glEnableVertexAttribArray(0)
      call glVertexAttribPointer(0, 3, GL_FLOAT, int(GL_FALSE,c_signed_char), int(3*c_sizeof(c_float_),c_int), c_null_ptr)
      call glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshebo)
      call glBindBuffer(GL_ARRAY_BUFFER, vbo)
      st = int(mesh_inst_nf * c_sizeof(c_float_),c_int)
      call inst_attrib(1, 4, st, 0)  ! model column 1
      call inst_attrib(2, 4, st, 4)  ! model column 2
      call inst_attrib(3, 4, st, 8)  ! model column 3
      call inst_attrib(4, 4, st, 12) ! model column 4
      call inst_attrib(5, 4, st, 16) ! color
      call glBindBuffer(GL_ARRAY_BUFFER, 0)
      call glBindVertexArray(0)

    end subroutine setup_mesh_inst
    !> Enable an instanced float attribute (divisor 1) at location loc with
    !> ncomp components, the given byte stride, and an offset of offf floats.
    subroutine inst_attrib(loc, ncomp, stride, offf)
      integer(c_int), intent(in) :: loc, ncomp, stride, offf

      call glEnableVertexAttribArray(loc)
      call glVertexAttribPointer(loc, ncomp, GL_FLOAT, int(GL_FALSE,c_signed_char), stride,&
         transfer(int(offf,c_intptr_t) * c_sizeof(c_float_), c_ptr_))
      call glVertexAttribDivisor(loc, 1_c_int)

    end subroutine inst_attrib

  end subroutine glbuffers_init

  !> Delete this scene's instance VAOs/VBOs. A no-op if not initialized.
  module subroutine glbuffers_end(b)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout), target :: b

    if (.not.b%isinit) return

    call glDeleteVertexArrays(1, c_loc(b%sphinstVAO))
    call glDeleteBuffers(1, c_loc(b%sphinstVBO))
    call glDeleteVertexArrays(1, c_loc(b%cylinstVAO))
    call glDeleteBuffers(1, c_loc(b%cylinstVBO))
    call glDeleteVertexArrays(1, c_loc(b%planeinstVAO))
    call glDeleteBuffers(1, c_loc(b%planeinstVBO))
    call glDeleteVertexArrays(1, c_loc(b%triinstVAO))
    call glDeleteBuffers(1, c_loc(b%triinstVBO))
    call glDeleteVertexArrays(1, c_loc(b%coneinstVAO))
    call glDeleteBuffers(1, c_loc(b%coneinstVBO))
    call glDeleteVertexArrays(1, c_loc(b%sphinstVAOscr))
    call glDeleteBuffers(1, c_loc(b%sphinstVBOscr))
    call glDeleteVertexArrays(1, c_loc(b%cylinstVAOscr))
    call glDeleteBuffers(1, c_loc(b%cylinstVBOscr))
    call glDeleteVertexArrays(1, c_loc(b%coneinstVAOscr))
    call glDeleteBuffers(1, c_loc(b%coneinstVBOscr))

    call b%detach()

  end subroutine glbuffers_end

  !> Reset the buffer handles and cached state to a fresh (uninitialized) state
  !> WITHOUT deleting any GL objects. Used after a whole-scene value copy so the
  !> copy lazily creates its own buffers instead of aliasing the source's.
  module subroutine glbuffers_detach(b)
    class(scene_glbuffers), intent(inout) :: b

    b%isinit = .false.
    b%sphinstVAO = 0
    b%sphinstVBO = 0
    b%cylinstVAO = 0
    b%cylinstVBO = 0
    b%coneinstVAO = 0
    b%coneinstVBO = 0
    b%planeinstVAO = 0
    b%planeinstVBO = 0
    b%triinstVAO = 0
    b%triinstVBO = 0
    b%sphinstVAOscr = 0
    b%sphinstVBOscr = 0
    b%cylinstVAOscr = 0
    b%cylinstVBOscr = 0
    b%coneinstVAOscr = 0
    b%coneinstVBOscr = 0
    b%inst_valid = .false.
    b%inst_last_anim = -1
    b%nsph_inst = 0
    b%ncyl_inst = 0
    b%ncone_inst = 0
    b%nplane_inst = 0
    b%ntri_inst = 0

  end subroutine glbuffers_detach

  !> Upload n packed sphere instances and draw them through the cached (or
  !> scratch) sphere impostor VAO with the currently bound shader (one instanced
  !> draw call).
  module subroutine glbuffers_draw_spheres(b,n,buf,scratch)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: n
    real(c_float), intent(in), target :: buf(sph_inst_nf,n)
    logical, intent(in) :: scratch
    integer(c_int) :: vao, vbo
    real(c_float) :: f_

    if (n <= 0) return
    if (scratch) then
       vao = b%sphinstVAOscr; vbo = b%sphinstVBOscr
    else
       vao = b%sphinstVAO; vbo = b%sphinstVBO
    end if
    call glBindVertexArray(vao)
    call glBindBuffer(GL_ARRAY_BUFFER, vbo)
    call glBufferData(GL_ARRAY_BUFFER, int(n,c_intptr_t)*sph_inst_nf*c_sizeof(f_), c_loc(buf), GL_DYNAMIC_DRAW)
    call glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0_c_int, 4_c_int, int(n,c_int))
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

  end subroutine glbuffers_draw_spheres

  !> Upload n packed cylinder instances and draw them through the cached (or
  !> scratch) cylinder impostor VAO with the currently bound shader.
  module subroutine glbuffers_draw_cylinders(b,n,buf,scratch)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: n
    real(c_float), intent(in), target :: buf(cyl_inst_nf,n)
    logical, intent(in) :: scratch
    integer(c_int) :: vao, vbo
    real(c_float) :: f_

    if (n <= 0) return
    if (scratch) then
       vao = b%cylinstVAOscr; vbo = b%cylinstVBOscr
    else
       vao = b%cylinstVAO; vbo = b%cylinstVBO
    end if
    call glBindVertexArray(vao)
    call glBindBuffer(GL_ARRAY_BUFFER, vbo)
    call glBufferData(GL_ARRAY_BUFFER, int(n,c_intptr_t)*cyl_inst_nf*c_sizeof(f_), c_loc(buf), GL_DYNAMIC_DRAW)
    call glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0_c_int, 4_c_int, int(n,c_int))
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

  end subroutine glbuffers_draw_cylinders

  !> Upload n packed mesh instances and draw them with the currently bound mesh
  !> shader, using the VAO/VBO selected by role (glb_cone/plane/tri/conescr) and
  !> the mesh's triangle count nelem.
  module subroutine glbuffers_draw_mesh(b,role,nelem,n,buf)
    use interfaces_opengl3
    use tools_io, only: ferror, faterr
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: role, nelem, n
    real(c_float), intent(in), target :: buf(mesh_inst_nf,n)
    integer(c_int) :: vao, vbo
    real(c_float) :: f_

    if (n <= 0) return
    select case (role)
    case (glb_cone)
       vao = b%coneinstVAO
       vbo = b%coneinstVBO
    case (glb_plane)
       vao = b%planeinstVAO
       vbo = b%planeinstVBO
    case (glb_tri)
       vao = b%triinstVAO
       vbo = b%triinstVBO
    case (glb_conescr)
       vao = b%coneinstVAOscr
       vbo = b%coneinstVBOscr
    case default
       call ferror('glbuffers_draw_mesh','unknown mesh role',faterr)
    end select
    call glBindVertexArray(vao)
    call glBindBuffer(GL_ARRAY_BUFFER, vbo)
    call glBufferData(GL_ARRAY_BUFFER, int(n,c_intptr_t)*mesh_inst_nf*c_sizeof(f_), c_loc(buf), GL_DYNAMIC_DRAW)
    call glDrawElementsInstanced(GL_TRIANGLES, int(3*nelem,c_int), GL_UNSIGNED_INT, c_null_ptr, int(n,c_int))
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)

  end subroutine glbuffers_draw_mesh

  !> Re-draw n already-uploaded sphere impostor instances through the cached VAO
  !> (no upload). Used when the cached buffers have not changed.
  module subroutine glbuffers_redraw_spheres(b,n)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: n

    if (n <= 0) return
    call glBindVertexArray(b%sphinstVAO)
    call glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0_c_int, 4_c_int, int(n,c_int))
    call glBindVertexArray(0)

  end subroutine glbuffers_redraw_spheres

  !> Re-draw n already-uploaded cylinder impostor instances through the cached
  !> VAO (no upload).
  module subroutine glbuffers_redraw_cylinders(b,n)
    use interfaces_opengl3
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: n

    if (n <= 0) return
    call glBindVertexArray(b%cylinstVAO)
    call glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0_c_int, 4_c_int, int(n,c_int))
    call glBindVertexArray(0)

  end subroutine glbuffers_redraw_cylinders

  !> Re-draw n already-uploaded mesh instances through the cached VAO selected by
  !> role (glb_cone/plane/tri), with nelem triangles per instance (no upload).
  module subroutine glbuffers_redraw_mesh(b,role,nelem,n)
    use interfaces_opengl3
    use tools_io, only: ferror, faterr
    class(scene_glbuffers), intent(inout) :: b
    integer, intent(in) :: role, nelem, n
    integer(c_int) :: vao

    if (n <= 0) return
    select case (role)
    case (glb_cone)
       vao = b%coneinstVAO
    case (glb_plane)
       vao = b%planeinstVAO
    case (glb_tri)
       vao = b%triinstVAO
    case default
       call ferror('glbuffers_redraw_mesh','unknown mesh role',faterr)
    end select
    call glBindVertexArray(vao)
    call glDrawElementsInstanced(GL_TRIANGLES, int(3*nelem,c_int), GL_UNSIGNED_INT, c_null_ptr, int(n,c_int))
    call glBindVertexArray(0)

  end subroutine glbuffers_redraw_mesh

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
