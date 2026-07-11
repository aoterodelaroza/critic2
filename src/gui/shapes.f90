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
  ! floats per on-scene-text vertex: anchor (3), eye shift (3), glyph corner
  ! (2), font-atlas UV (2), anchor vibration delta re (3) and im (3)
  integer(c_int), parameter, public :: text_vert_nf = 16

  ! impostor objects (instanced billboards). A single unit quad (corners in
  ! [-1,1], drawn as a triangle strip) is expanded per instance; the fragment
  ! shader ray-casts the actual surface. The shared unit-quad corners live here
  ! (template geometry); the per-instance VAOs/VBOs are per-scene, in the
  ! scene_glbuffers type below.
  integer(c_int), target, public :: impQuadVBO ! shared unit-quad corners
  integer, parameter, public :: maxpie = 4 ! max colored occupancy-pie sectors (rest grey)
  integer(c_int), parameter, public :: sph_inst_nf = 38 ! floats per sphere instance
  integer(c_int), parameter, public :: cyl_inst_nf = 31 ! floats per cylinder instance

  ! instanced plain meshes (planes, polyhedra triangles, cones): per-instance
  ! model matrix (4 columns) + color + per-vertex vibration deltas (real and
  ! imaginary parts for the three reference-triangle corners; zero for planes
  ! and cones, which are not animated)
  integer(c_int), parameter, public :: mesh_inst_nf = 38 ! floats per mesh instance (16 model + 4 color + 18 deltas)

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
     real(c_float) :: occ = 1._c_float ! occupancy of pie sector 1 (representative)
     real(c_float) :: occ_empty_rgb(3) = 0._c_float ! color of the empty (unoccupied) sector
     real(c_float) :: pie_cum(3) = 1._c_float ! cumulative sector boundaries t2,t3,ttot (mixed sites)
     real(c_float) :: pie_rgb(3,3) = 0._c_float ! colors of pie sectors 2,3,4 (mixed sites)
     logical :: ghost = .false. ! invisible pick-only target (atoms hidden, bonds shown)
  end type dl_sphere
  public :: dl_sphere

  !> cylinders for the draw list. Also used for cones (arrowheads), with
  !> x1 = base center and x2 = apex.
  type dl_cylinder
     real(c_float) :: x1(3) ! one end of the cylinder
     real(c_float) :: x2(3) ! other end of the cylinder
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     complex(c_float_complex) :: x1delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (end 1)
     complex(c_float_complex) :: x2delta(3) = (0._c_float,0._c_float) ! delta-vector for vibration animations (end 2)
     integer(c_int) :: order = 1 ! order of the bond (-1=aromatic/1.5,0=dashed,1=single,2=double,3=triple); < -1 = flat cylinder
     real(c_float) :: border = 0._c_float ! border size
     real(c_float) :: rgbborder(3) = 0._c_float ! border color
     real(c_float) :: arvec(3) = 0._c_float ! aromatic bonds: outward (ring-exterior) direction; 0 otherwise
     real(c_float) :: alpha = 1._c_float ! opacity (1 = opaque)
  end type dl_cylinder
  public :: dl_cylinder

  !> window-anchored overlay cylinder/cone: a dl_cylinder plus its window placement
  type, extends(dl_cylinder) :: dl_cylinder_over
     real(c_float) :: winpos(2) = 0._c_float ! window position (fractions from left and bottom)
     logical :: scalewithzoom = .false. ! whether this overlay item scales with zoom
  end type dl_cylinder_over
  public :: dl_cylinder_over

  !> strings for the draw list
  type dl_string
     real(c_float) :: x(3) ! position
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     real(c_float) :: scale ! scale (1.0 = radius)
     real(c_float) :: offset(3) ! offset of the label in pixels
     complex(c_float_complex) :: xdelta(3) ! delta-vector for vibration animations
     character(len=:), allocatable :: str ! string
  end type dl_string
  public :: dl_string

  !> window-anchored overlay string: a dl_string plus its window placement
  type, extends(dl_string) :: dl_string_over
     real(c_float) :: winpos(2) = 0._c_float ! window position (fractions from left and bottom)
     logical :: scalewithzoom = .false. ! whether this overlay item scales with zoom
  end type dl_string_over
  public :: dl_string_over

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
     type(dl_string), allocatable :: string(:) ! string draw list
     ! window-anchored overlay objects (drawn in a separate overlay pass;
     ! currently the axes gizmo, but not restricted to it)
     integer :: ncylover ! number of overlay cylinders (e.g. axes gizmo shafts)
     type(dl_cylinder_over), allocatable :: cylover(:) ! overlay cylinder draw list
     integer :: nconeover ! number of overlay cones (e.g. axes gizmo arrowheads)
     type(dl_cylinder_over), allocatable :: coneover(:) ! overlay cone draw list
     integer :: nstringover ! number of overlay strings (e.g. axes gizmo labels)
     type(dl_string_over), allocatable :: stringover(:) ! overlay string draw list
   contains
     procedure :: reset => scene_objects_reset
     procedure :: reserve => scene_objects_reserve
     procedure :: end => scene_objects_end
  end type scene_objects
  public :: scene_objects

  !> Ensure a pack scratch array is allocated with at least the requested
  !> extent (grow-only; the contents are scratch and are not preserved).
  interface ensure_pack
     module procedure ensure_pack_realc2
     module procedure ensure_pack_int1
  end interface ensure_pack
  public :: ensure_pack

  !> Append an item to a draw list, growing (or allocating) the list as needed.
  interface dl_append
     module procedure dl_append_sphere
     module procedure dl_append_cylinder
     module procedure dl_append_cylinder_over
     module procedure dl_append_string
     module procedure dl_append_string_over
     module procedure dl_append_plane
     module procedure dl_append_triangle
  end interface dl_append
  public :: dl_append

  ! mesh-buffer role selectors for the scene_glbuffers draw_mesh/redraw_mesh methods
  integer, parameter, public :: glb_cone = 1    ! cached cone (arrowhead) mesh
  integer, parameter, public :: glb_plane = 2   ! cached plane (rectangle) mesh
  integer, parameter, public :: glb_tri = 3     ! cached triangle (polyhedron face) mesh
  integer, parameter, public :: glb_conescr = 4 ! scratch cone mesh (overlay cones)

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
     ! highlights, picking, overlay)
     integer(c_int) :: sphinstVAOscr = 0, sphinstVBOscr = 0
     integer(c_int) :: cylinstVAOscr = 0, cylinstVBOscr = 0
     integer(c_int) :: coneinstVAOscr = 0, coneinstVBOscr = 0
     ! VBO storage capacities (in instances). Uploads orphan the current
     ! storage and glBufferSubData into it; the storage is re-created (2x)
     ! only when the instance count exceeds the capacity.
     integer :: sph_cap = 0, cyl_cap = 0, cone_cap = 0, plane_cap = 0, tri_cap = 0
     integer :: sphscr_cap = 0, cylscr_cap = 0, conescr_cap = 0
     ! persistent CPU pack scratch (grow-only, see ensure_pack): instance data
     ! is packed here before upload, avoiding per-frame allocations
     real(c_float), allocatable :: packsph(:,:)  ! sphere instances (also selections/highlights/pick)
     real(c_float), allocatable :: packcyl(:,:)  ! cylinder instances
     real(c_float), allocatable :: packmesh(:,:) ! mesh instances (cones, then planes, then triangles)
     ! per-scene on-scene-text buffer: the glyph vertices of the scene labels
     ! are cached here (and in the VBO) and reused while the camera and the
     ! draw lists do not change
     integer(c_int) :: textVAO = 0, textVBO = 0 ! cached label glyph buffer (text_vert_nf floats/vertex)
     integer :: text_cap = 0 ! vertex capacity of the text VBO storage
     integer, allocatable :: text_first(:) ! first vertex of each label (0-based)
     integer, allocatable :: text_count(:) ! vertex count of each label
     real(c_float), allocatable :: packtext(:,:) ! concatenated glyph vertices (text_vert_nf,:)
     logical :: text_valid = .false. ! cached text vertices are current
     real*8 :: text_build_time = -1d0 ! draw-list build time at the last text rebuild
     real(c_float) :: text_proj(4,4) = 0._c_float ! projection matrix at the last text rebuild
     real(c_float) :: text_vw3(4) = 0._c_float ! view*world third row at the last text rebuild
     real(c_float) :: text_hside = -1._c_float ! reset-zoom half-side at the last text rebuild
     real(c_float) :: text_fontsize = -1._c_float ! baked font size at the last text rebuild
     ! cached-buffer validity and instance counts
     logical :: inst_valid = .false. ! true if the cached instance buffers are current
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
     procedure :: upload_text => glbuffers_upload_text
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
     module subroutine scene_objects_reset(o)
       class(scene_objects), intent(inout) :: o
     end subroutine scene_objects_reset
     module subroutine scene_objects_reserve(o,nsph,ncyl,nstring)
       class(scene_objects), intent(inout) :: o
       integer, intent(in), optional :: nsph, ncyl, nstring
     end subroutine scene_objects_reserve
     module subroutine scene_objects_end(o)
       class(scene_objects), intent(inout) :: o
     end subroutine scene_objects_end
     module subroutine ensure_pack_realc2(arr,nf,n)
       real(c_float), allocatable, intent(inout) :: arr(:,:)
       integer, intent(in) :: nf, n
     end subroutine ensure_pack_realc2
     module subroutine ensure_pack_int1(arr,n)
       integer, allocatable, intent(inout) :: arr(:)
       integer, intent(in) :: n
     end subroutine ensure_pack_int1
     module subroutine dl_append_sphere(lst,n,it)
       type(dl_sphere), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_sphere), intent(in) :: it
     end subroutine dl_append_sphere
     module subroutine dl_append_cylinder(lst,n,it)
       type(dl_cylinder), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_cylinder), intent(in) :: it
     end subroutine dl_append_cylinder
     module subroutine dl_append_cylinder_over(lst,n,it)
       type(dl_cylinder_over), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_cylinder_over), intent(in) :: it
     end subroutine dl_append_cylinder_over
     module subroutine dl_append_string(lst,n,it)
       type(dl_string), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_string), intent(in) :: it
     end subroutine dl_append_string
     module subroutine dl_append_string_over(lst,n,it)
       type(dl_string_over), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_string_over), intent(in) :: it
     end subroutine dl_append_string_over
     module subroutine dl_append_plane(lst,n,it)
       type(dl_plane), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_plane), intent(in) :: it
     end subroutine dl_append_plane
     module subroutine dl_append_triangle(lst,n,it)
       type(dl_triangle), allocatable, intent(inout) :: lst(:)
       integer, intent(inout) :: n
       type(dl_triangle), intent(in) :: it
     end subroutine dl_append_triangle
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
     module subroutine glbuffers_upload_text(b,nvert,buf)
       class(scene_glbuffers), intent(inout) :: b
       integer, intent(in) :: nvert
       real(c_float), intent(in), target :: buf(text_vert_nf,nvert)
     end subroutine glbuffers_upload_text
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

