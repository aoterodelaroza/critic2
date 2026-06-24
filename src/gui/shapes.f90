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
module shapes
  use iso_c_binding
  implicit none

  private

  public :: shapes_init
  public :: shapes_end

  ! cone objects (arrowheads)
  integer, parameter, public :: nmaxcone = 3
  integer(c_int), target, public :: coneVAO(nmaxcone) ! cone: vertex array object
  integer(c_int), target, public :: coneVBO(nmaxcone) ! cone: vertex buffer object
  integer(c_int), target, public :: coneEBO(nmaxcone) ! cone: element buffer object

  ! cone dimensions (apex + base center + base ring of 6/12/18 vertices)
  integer(c_int), parameter, public :: connve(nmaxcone) = (/8,14,20/)
  integer(c_int), parameter, public :: connveadd(0:nmaxcone) = (/0,8,22,42/)
  integer(c_int), parameter, public :: connel(nmaxcone) = (/12,24,36/)
  integer(c_int), parameter, public :: conneladd(0:nmaxcone) = (/0,12,36,72/)

  ! quad objects (filled flat rectangles): unit quad in the z=0 plane, corners
  ! at (+/-1,+/-1,0), two triangles
  integer(c_int), target, public :: quadVAO ! quad: vertex array object
  integer(c_int), target, public :: quadVBO ! quad: vertex buffer object
  integer(c_int), target, public :: quadEBO ! quad: element buffer object
  integer(c_int), parameter, public :: quadnel = 2 ! quad: number of triangles

  ! triangle objects (filled flat triangles): unit reference triangle in the
  ! z=0 plane, corners at (0,0,0),(1,0,0),(0,1,0). A model matrix whose first
  ! two columns are (v2-v1) and (v3-v1) maps these corners onto any triangle.
  integer(c_int), target, public :: triVAO ! triangle: vertex array object
  integer(c_int), target, public :: triVBO ! triangle: vertex buffer object
  integer(c_int), target, public :: triEBO ! triangle: element buffer object
  integer(c_int), parameter, public :: trinel = 1 ! triangle: number of triangles

  ! text objects
  integer(c_int), target, public :: textVAO
  integer(c_int), target, public :: textVBO
  integer(c_int), target, public :: textVAOos
  integer(c_int), target, public :: textVBOos
  integer(c_int), parameter, public :: text_maxvert = 6 * 128

  ! impostor objects (instanced billboards). A single unit quad (corners in
  ! [-1,1], drawn as a triangle strip) is expanded per instance; the fragment
  ! shader ray-casts the actual surface. The shared unit-quad corners live here
  ! (template geometry); the per-instance VAOs/VBOs are per-scene, in the
  ! scene_glbuffers type below.
  integer(c_int), target, public :: impQuadVBO ! shared unit-quad corners
  integer(c_int), parameter, public :: sph_inst_nf = 22 ! floats per sphere instance
  integer(c_int), parameter, public :: cyl_inst_nf = 19 ! floats per cylinder instance

  ! instanced plain meshes (planes, polyhedra triangles, cones): per-instance
  ! model matrix (4 columns) + color
  integer(c_int), parameter, public :: mesh_inst_nf = 20 ! floats per mesh instance (16 model + 4 color)

  !! draw list objects
  !> spheres for the draw list
  type dl_sphere
     real(c_float) :: x(3) ! position
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     integer(c_int) :: idx(4) ! atom ID (complete atom list) + lattice vector
     complex(c_float_complex) :: xdelta(3) ! delta-vector for vibration animations
     real(c_float) :: border ! border size
     real(c_float) :: rgbborder(3) ! border color
     logical :: ghost = .false. ! invisible pick-only target (atoms hidden, bonds shown)
  end type dl_sphere
  public :: dl_sphere

  !> cylinders for the draw list
  type dl_cylinder
     real(c_float) :: x1(3) ! one end of the cylinder
     real(c_float) :: x2(3) ! other end of the cylinder
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     complex(c_float_complex) :: x1delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (end 1)
     complex(c_float_complex) :: x2delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (end 2)
     integer(c_int) :: order ! order of the bond (-1=aromatic/1.5,0=dashed,1=single,2=double,3=triple); < -1 = flat cylinder
     real(c_float) :: border ! border size
     real(c_float) :: rgbborder(3) ! border color
     real(c_float) :: arvec(3) = 0._c_float ! aromatic bonds: outward (ring-exterior) direction; 0 otherwise
     real(c_float) :: alpha = 1._c_float ! opacity (1 = opaque)
     real(c_float) :: gizwinpos(2) = 0._c_float ! window position (fractions) for window-anchored gizmo items
     logical :: gizscalewithzoom = .false. ! whether this gizmo item scales with zoom (gizmo items only)
  end type dl_cylinder
  public :: dl_cylinder

  !> strings for the draw list
  type dl_string
     real(c_float) :: x(3) ! position
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     real(c_float) :: scale ! scale (1.0 = radius)
     real(c_float) :: offset(3) ! offset of the label in pixels
     complex(c_float_complex) :: xdelta(3) ! delta-vector for vibration animations
     character(len=:), allocatable :: str ! string
     real(c_float) :: gizwinpos(2) = 0._c_float ! window position (fractions) for window-anchored gizmo items
     logical :: gizscalewithzoom = .false. ! whether this gizmo item scales with zoom (gizmo items only)
  end type dl_string
  public :: dl_string

  !> filled flat rectangles (quads) for the draw list
  type dl_plane
     real(c_float) :: x(3) ! center position
     real(c_float) :: e1(3) ! x-axis half-vector
     real(c_float) :: e2(3) ! y-axis half-vector
     real(c_float) :: rgb(3) ! color
     real(c_float) :: alpha = 1._c_float ! opacity (1 = opaque)
  end type dl_plane
  public :: dl_plane

  !> filled flat triangles for the draw list (e.g. coordination polyhedra faces)
  type dl_triangle
     real(c_float) :: x1(3) ! first vertex
     real(c_float) :: x2(3) ! second vertex
     real(c_float) :: x3(3) ! third vertex
     real(c_float) :: rgb(3) ! color
     real(c_float) :: alpha = 1._c_float ! opacity (1 = opaque)
     complex(c_float_complex) :: x1delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (vertex 1)
     complex(c_float_complex) :: x2delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (vertex 2)
     complex(c_float_complex) :: x3delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (vertex 3)
  end type dl_triangle
  public :: dl_triangle

  !> collection of objects (draw lists) belonging to a scene
  type scene_objects
     integer :: nsph ! number of spheres
     type(dl_sphere), allocatable :: sph(:) ! sphere draw list
     integer :: ncyl ! number of cylinders
     type(dl_cylinder), allocatable :: cyl(:) ! cylinder draw list
     integer :: ncylflat ! number of flat cylinders
     type(dl_cylinder), allocatable :: cylflat(:) ! flat cylinder draw list
     integer :: ncone ! number of cones (arrowheads)
     type(dl_cylinder), allocatable :: cone(:) ! cone draw list (x1=base, x2=apex)
     integer :: nplane ! number of flat rectangles
     type(dl_plane), allocatable :: plane(:) ! flat rectangle draw list
     integer :: ntriangle ! number of flat triangles
     type(dl_triangle), allocatable :: triangle(:) ! flat triangle draw list
     integer :: nstring ! number of strings
     type(dl_string), allocatable :: string(:) ! flat cylinder draw list
     ! window-anchored axes gizmo (drawn in a separate overlay pass)
     integer :: ncylgiz ! number of gizmo cylinders (axis shafts)
     type(dl_cylinder), allocatable :: cylgiz(:) ! gizmo cylinder draw list
     integer :: nconegiz ! number of gizmo cones (arrowheads)
     type(dl_cylinder), allocatable :: conegiz(:) ! gizmo cone draw list
     integer :: nstringgiz ! number of gizmo strings (axis labels)
     type(dl_string), allocatable :: stringgiz(:) ! gizmo string draw list
  end type scene_objects
  public :: scene_objects

  ! mesh-buffer role selectors for the scene_glbuffers draw_mesh/redraw_mesh methods
  integer, parameter, public :: glb_cone = 1    ! cached cone (arrowhead) mesh
  integer, parameter, public :: glb_plane = 2   ! cached plane (rectangle) mesh
  integer, parameter, public :: glb_tri = 3     ! cached triangle (polyhedron face) mesh
  integer, parameter, public :: glb_conescr = 4 ! scratch cone mesh (gizmo arrowheads)

  ! per-scene OpenGL instance buffers.
  type scene_glbuffers
     logical :: isinit = .false. ! GL buffers created?
     ! cached instance buffers (one scene's worth)
     integer(c_int) :: sphinstVAO = 0, sphinstVBO = 0    ! sphere impostors
     integer(c_int) :: cylinstVAO = 0, cylinstVBO = 0    ! cylinder impostors
     integer(c_int) :: coneinstVAO = 0, coneinstVBO = 0  ! cone mesh
     integer(c_int) :: planeinstVAO = 0, planeinstVBO = 0 ! quad mesh
     integer(c_int) :: triinstVAO = 0, triinstVBO = 0    ! triangle mesh
     ! scratch instance buffers (re-uploaded every frame: measure-selection,
     ! highlights, picking, gizmo)
     integer(c_int) :: sphinstVAOscr = 0, sphinstVBOscr = 0
     integer(c_int) :: cylinstVAOscr = 0, cylinstVBOscr = 0
     integer(c_int) :: coneinstVAOscr = 0, coneinstVBOscr = 0
     ! cached-buffer validity and instance counts
     logical :: inst_valid = .false. ! true if the cached instance buffers are current
     integer :: inst_last_anim = -1 ! animation state at the last instance-buffer build
     integer :: nsph_inst = 0   ! number of cached atom-sphere instances
     integer :: ncyl_inst = 0   ! number of cached bond/cell-cylinder instances
     integer :: ncone_inst = 0  ! number of cached cone instances
     integer :: nplane_inst = 0 ! number of cached plane instances
     integer :: ntri_inst = 0   ! number of cached triangle instances
   contains
     procedure :: init => glbuffers_init
     procedure :: end => glbuffers_end
     procedure :: detach => glbuffers_detach
     procedure :: draw_spheres => glbuffers_draw_spheres
     procedure :: draw_cylinders => glbuffers_draw_cylinders
     procedure :: draw_mesh => glbuffers_draw_mesh
     procedure :: redraw_spheres => glbuffers_redraw_spheres
     procedure :: redraw_cylinders => glbuffers_redraw_cylinders
     procedure :: redraw_mesh => glbuffers_redraw_mesh
  end type scene_glbuffers
  public :: scene_glbuffers

  ! module procedure interfaces
  interface
     module subroutine shapes_init()
     end subroutine shapes_init
     module subroutine shapes_end()
     end subroutine shapes_end
     module subroutine glbuffers_init(b)
       class(scene_glbuffers), intent(inout), target :: b
     end subroutine glbuffers_init
     module subroutine glbuffers_end(b)
       class(scene_glbuffers), intent(inout), target :: b
     end subroutine glbuffers_end
     module subroutine glbuffers_detach(b)
       class(scene_glbuffers), intent(inout) :: b
     end subroutine glbuffers_detach
     module subroutine glbuffers_draw_spheres(b,n,buf,scratch)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: n
       real(c_float), intent(in), target :: buf(sph_inst_nf,n)
       logical, intent(in) :: scratch
     end subroutine glbuffers_draw_spheres
     module subroutine glbuffers_draw_cylinders(b,n,buf,scratch)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: n
       real(c_float), intent(in), target :: buf(cyl_inst_nf,n)
       logical, intent(in) :: scratch
     end subroutine glbuffers_draw_cylinders
     module subroutine glbuffers_draw_mesh(b,role,nelem,n,buf)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: role, nelem, n
       real(c_float), intent(in), target :: buf(mesh_inst_nf,n)
     end subroutine glbuffers_draw_mesh
     module subroutine glbuffers_redraw_spheres(b,n)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: n
     end subroutine glbuffers_redraw_spheres
     module subroutine glbuffers_redraw_cylinders(b,n)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: n
     end subroutine glbuffers_redraw_cylinders
     module subroutine glbuffers_redraw_mesh(b,role,nelem,n)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: role, nelem, n
     end subroutine glbuffers_redraw_mesh
  end interface

end module shapes

