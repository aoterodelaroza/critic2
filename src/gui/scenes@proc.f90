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

! ---- Coordinate systems used for scene rendering ----
! 1. Object coordinates: these are the coordinates local to the object.
!
! 1 -> 2 with MODEL MATRIX: The model matrix contains the scaling,
! rotation and translation necessary to bring the object to its size,
! position, and orientation in the scene.
!
! 2. World coordinates: the coordinates corresponding to the chemical
! system. Identical to the Cartesian coordinates (%r) used by
! critic2, in bohr.
!
! 2 -> 3 with WORLD MATRIX: a translation + rotation matrix that is
! used to rotate and translate all objects at the same time.
!
! 3. Transformed world coordinates (tworld): the world coordinates
! after rotation/translation of the whole scene
!
! 3 -> 4 with VIEW MATRIX: the view matrix is calculated from the
! camera position (s%campos), the vector that indicates the direction
! of the camera (s%camfront), and the up vector (s%camup) using the
! lookat routine.
!
! 4. Eye/view coordinates: coordinate system in which the camera is at
! (0,0,0) and points in the -z direction, with the up vector being
! the +y direction.
!
! 4 -> 5 with PROJECTION MATRIX: the projection matrix is constructed
! by selecting the visible region (fustrum) and mapping that into the
! -1:1 cube.
!
! 5. Clip coordinates: homogeneous coordinates where the "visible"
! space is between -w and +w.
!
! 5 -> 6 Divide by w
!
! 6. Normalized device coordinates (NDC): the visible region is
! between -1 and +1 in all directions.
!
! 6 -> 7 Scale the -1:1 rectangle to the texture size.
!   xtex = (0.5 * xndc + 0.5) * side
!
! 7. Texture coordinates: mapping of the NDC into texture pixels,
! between 0 and the side of the texture (Tx,Ty).
!
! 7 -> 8 !! texpos to mousepos transformation !!
!   texpos is (0,ty)   (tx,ty)
!                  +---+
!                  |   |
!                  +---+
!             (0,0)    (tx,0)
!   First, transform to normalized coordinates:
!     xnorm = 2 * tx/TX - 1
!     ynorm = 1 - 2 * ty/TY
!   norm  is (-1,-1)   (1,-1)
!                  +---+
!                  |   |
!                  +---+
!             (-1,1)   (1,1)
!   Then, transform to mouse coordinaes:
!     xmouse = xmin + 1/2 * (xmax-xmin) + 1/2 * xnorm * max(xmax-xmin,ymax-ymin)
!     ymouse = ymin + 1/2 * (ymax-ymin) + 1/2 * ynorm * max(xmax-xmin,ymax-ymin)
!   where (xmin,ymin) is the top left corner of the window and
!   (xmax,ymax) is the bottom right corner of the window.
!
! 8. Mouse coordinates: coordinates of the mouse on the screen as
! determined by imgui.
!

! ---- Standard camera movements in the critic2 GUI ----
!
! 1. Zoom (default: mouse wheel). First, transform the scene center
! to tworld coordinates. In tworld coordinates, calcualte the vector
! from the scene center to the camera and reduce it or extend
! it by a factor ("ratio"). Clip the resulting vector so it is not
! shorter than min_zoom or longer than max_zoom. Move the camera to
! the new position indicated by this vector.
!
! 2. Drag (default: right button). The first time the texture is
! clicked, save the view coordinates corresponding to the mouse
! position, and the view matrix. The mouse is now moved with the right
! button down; the new texture position is also transformed to view
! coordinates.  Calculate the difference vector betwen the new and the
! old view coordinates and transform it to tworld coordinates using
! the saved view matrix. Move the camera by this vector.
!
! 3. Rotate (default: left button). The first time the texture is
! clicked, save the view and texture coordintes coordinates of the
! mouse position. When the mouse is moved, take the difference between
! the new position in view coordinates and the old position. The cross
! product of (0,0,1), which is camera front-pointing vector in view
! coordinates, with the difference vector gives the axis of rotation
! in view coordinates. For the angle, calculate the difference between
! the new and the old texture positions and divide by the side of the
! texture to normalize. The angle is this value multiplied by the
! sensitivity to rotation.
!
! 4. Rotate around perpendicular axis: the first time the scene is
! clicked, save the scene center in texture coordinates and
! the normalized difference vector between the mouse position and the
! scene center, also in texture coordinates. When the mouse is moved,
! calculate the normalized difference vector wrt the scene center
! with the new position and then the angle between this vector
! and the saved one. The camera is rotated around its front-pointing
! vector by this angle.
!

! Scene object and GL rendering utilities
submodule (scenes) proc
  implicit none

  ! some math parameters
  real(c_float), parameter :: zero = 0._c_float
  real(c_float), parameter :: one = 1._c_float
  real(c_float), parameter :: eye4(4,4) = reshape((/&
     one,zero,zero,zero,&
     zero,one,zero,zero,&
     zero,zero,one,zero,&
     zero,zero,zero,one/),shape(eye4))

  !xx! private procedures: low-level draws
  ! subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,centered)

contains

  !xx! scene

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use representations, only: reptype_atoms, reptype_unitcell, reptype_axes,&
       repflavor_atoms_ballandstick, repflavor_atoms_criticalpoints, repflavor_atoms_gradientpaths,&
       repflavor_atoms_sticks, repflavor_unitcell_basic, repflavor_axes, repflavor_NUM
    use systems, only: sys, sysc, sys_ready, ok_system
    use global, only: crsmall
    use gui_main, only: lockbehavior
    use param, only: maxzat, maxzat0
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! basic variables
    s%id = isys
    s%isinit = 1
    s%iscaminit = .false.
    s%nc = 1
    s%scenerad = 10d0
    s%scenecenter = 0d0
    s%scenexmin = 0d0
    s%scenexmax = 1d0
    s%forcebuildlists = .true.
    s%iqpt_selected = 0
    s%ifreq_selected = 0
    s%animation = 0
    s%anim_speed = anim_speed_default
    s%anim_amplitude = anim_amplitude_default

    ! resolutions
    s%atom_res = 4
    s%bond_res = 1
    s%uc_res = 1

    ! appearance default settings
    call s%set_style_defaults()

    ! measure atom sets
    s%nmsel = 0
    s%msel = 0

    ! initialize representations
    if (allocated(s%rep)) deallocate(s%rep)
    allocate(s%rep(20))
    s%nrep = 0
    if (allocated(s%icount)) deallocate(s%icount)
    allocate(s%icount(0:repflavor_NUM))
    s%icount = 0
    if (allocated(s%iord)) deallocate(s%iord)
    allocate(s%iord(20))
    s%iord = 0

    ! transient representations (none on init)
    if (allocated(s%reptrans)) deallocate(s%reptrans)
    s%nreptrans = 0
    s%reptrans_set = .false.
    s%reptrans_tag = -1

    ! atoms
    if (sys(isys)%c%ncel <= crsmall) then
       call s%add_representation(reptype_atoms,repflavor_atoms_ballandstick)
    else
       call s%add_representation(reptype_atoms,repflavor_atoms_sticks)
    end if

    ! unit cell
    if (.not.sys(isys)%c%ismolecule) &
       call s%add_representation(reptype_unitcell,repflavor_unitcell_basic)

    ! cartesian axes (shown by default for molecules)
    if (sys(isys)%c%ismolecule) &
       call s%add_representation(reptype_axes,repflavor_axes)

    ! critical points
    if (any(sys(isys)%c%spc(:)%z > maxzat)) then
       call s%add_representation(reptype_atoms,repflavor_atoms_criticalpoints)
    end if

    ! gradient paths
    if (any(sys(isys)%c%spc(:)%z == maxzat0)) then
       call s%add_representation(reptype_atoms,repflavor_atoms_gradientpaths)
    end if

    ! reset the camera later
    s%camresetdist = 1.5_c_float
    s%camratio = 1.5_c_float

    ! sort the representations next pass
    s%forcesort = .true.
    s%timelastrender = 0d0
    s%timelastbuild = 0d0
    s%timelastcamchange = 0d0

    ! locking group for the camera
    s%forceresetcam = .true.
    if (lockbehavior == 0) then ! no-lock
       s%lockedcam = 0
    elseif (lockbehavior == 2) then ! all-lock
       s%lockedcam = -1
    else ! ==1, only scf
       if (sysc(isys)%collapse < 0) then
          s%lockedcam = isys
       elseif (sysc(isys)%collapse > 0) then
          s%lockedcam = sysc(isys)%collapse
       else
          s%lockedcam = 0
       end if
    end if

  end subroutine scene_init

  !> Terminate a scene object
  module subroutine scene_end(s)
    class(scene), intent(inout), target :: s

    s%isinit = 0
    s%iscaminit = .false.
    s%id = 0
    call s%gl%end()
    call s%obj%end()
    if (allocated(s%rep)) deallocate(s%rep)
    if (allocated(s%icount)) deallocate(s%icount)
    if (allocated(s%iord)) deallocate(s%iord)
    if (allocated(s%reptrans)) deallocate(s%reptrans)
    s%nreptrans = 0
    s%reptrans_set = .false.
    s%reptrans_tag = -1
    s%nrep = 0
    s%iqpt_selected = 0
    s%ifreq_selected = 0
    s%animation = 0
    s%anim_speed = anim_speed_default
    s%anim_amplitude = anim_amplitude_default

  end subroutine scene_end

  !> Reset the camera position and direction. Sets scenerad,
  !> scenecenter, ortho_fov, persp_fov, v_center, v_up, v_pos, view,
  !> world, projection, and znear.
  module subroutine scene_reset(s)
    use param, only: pi
    use interfaces_glfw, only: glfwGetTime
    class(scene), intent(inout), target :: s

    real(c_float) :: pic
    real*8 :: hside

    ! default transformation matrices
    pic = real(pi,c_float)
    s%ortho_fov = 45._c_float
    s%persp_fov = 45._c_float

    ! world matrix
    s%world = eye4

    ! camera distance and view matrix
    hside = reset_zoom_hside(s)
    s%campos = s%scenecenter
    s%campos(3) = s%campos(3) + real(hside,c_float) / tan(0.5_c_float * s%ortho_fov * pic / 180._c_float)
    s%camfront = (/zero,zero,-one/)
    s%camup = (/zero,one,zero/)
    call s%update_view_matrix()

    ! projection matrix
    call s%update_projection_matrix()

    ! the camera has been updated
    s%timelastcamchange = glfwGetTime()
    s%forceresetcam = .false.
    s%iscaminit = .true.

  end subroutine scene_reset

  !> Reset animation parameters in the scene
  module subroutine scene_reset_animation(s)
    class(scene), intent(inout), target :: s

    s%animation = 0
    s%anim_speed = anim_speed_default
    s%anim_amplitude = anim_amplitude_default

  end subroutine scene_reset_animation

  !> Reset atom colors in the scene to the defaults
  module subroutine scene_reset_atom_colors(s)
    class(scene), intent(inout), target :: s

    integer :: irep

    if (.not.allocated(s%rep)) return

    do irep = 1, s%nrep
       call s%rep(irep)%atoms%style%reset_colors(s%rep(irep))
    end do
    s%forcebuildlists = .true.

  end subroutine scene_reset_atom_colors

  !> Build the draw lists for the current scene.
  module subroutine scene_build_lists(s)
    use representations, only: reptype_atoms, reptype_axes, reptype_symelem, axes_winfrac_def
    use interfaces_glfw, only: glfwGetTime
    use utils, only: translate
    use systems, only: sys_ready, ok_system
    use interfaces_glfw, only: glfwGetTime
    class(scene), intent(inout), target :: s

    integer :: i
    real(c_float) :: xmin(3), xmax(3), maxrad, xc(3)

    ! only build lists if system is initialized
    if (.not.ok_system(s%id,sys_ready)) return

    ! reset the draw lists: counters to zero, allocations kept and reused
    ! across rebuilds (freed only in scene_end)
    call s%obj%reset()

    ! add the items by representation. Window-anchored axes are deferred: they
    ! need the scene radius (computed below) to auto-size, and they live in the
    ! gizmo draw lists, which are excluded from the scene bounding box anyway.
    do i = 1, s%nrep
       ! update to reflect changes in the number of atoms or molecules
       call s%rep(i)%update()

       ! add draw elements (except window-anchored axes and persistent symmetry
       ! elements, which need the scene radius and are done below)
       if (s%rep(i)%type == reptype_axes .and. s%rep(i)%axes%placement == 1) cycle
       if (s%rep(i)%type == reptype_symelem) cycle
       call s%rep(i)%add_draw_elements(s%nc,s%obj,s%animation>0,s%iqpt_selected,s%ifreq_selected)
    end do

    ! reset the measure selection
    s%nmsel = 0
    s%msel = 0

    ! recalculate scene center and radius
    maxrad = 0._c_float
    do i = 1, s%nrep
       if (s%rep(i)%shown .and. s%rep(i)%type == reptype_atoms) then
          if (s%rep(i)%atoms%style%ntype > 0) then
             maxrad = max(maxrad,real(maxval(s%rep(i)%atoms%style%rad(1:s%rep(i)%atoms%style%ntype)),c_float))
          end if
       end if
    end do

    ! ghost (pick-only) spheres do not count as framing geometry: a scene whose
    ! only geometry is ghosts falls through to the default extent below
    if (count(.not.s%obj%sph(1:s%obj%nsph)%ghost) + s%obj%ncyl + s%obj%ncylflat +&
       s%obj%ncone + s%obj%nstring > 0) then
       do i = 1, 3
          xmin(i) = huge(1._c_float)
          xmax(i) = -huge(1._c_float)

          ! exclude ghost (pick-only) spheres so they do not affect framing
          xmin(i) = minval(s%obj%sph(1:s%obj%nsph)%x(i),&
             mask=.not.s%obj%sph(1:s%obj%nsph)%ghost) - maxrad
          xmin(i) = min(xmin(i),minval(s%obj%cyl(1:s%obj%ncyl)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cyl(1:s%obj%ncyl)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cylflat(1:s%obj%ncylflat)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cylflat(1:s%obj%ncylflat)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cone(1:s%obj%ncone)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cone(1:s%obj%ncone)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%obj%string(1:s%obj%nstring)%x(i)))

          xmax(i) = maxval(s%obj%sph(1:s%obj%nsph)%x(i),&
             mask=.not.s%obj%sph(1:s%obj%nsph)%ghost) + maxrad
          xmax(i) = max(xmax(i),maxval(s%obj%cyl(1:s%obj%ncyl)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cyl(1:s%obj%ncyl)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cylflat(1:s%obj%ncylflat)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cylflat(1:s%obj%ncylflat)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cone(1:s%obj%ncone)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cone(1:s%obj%ncone)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%string(1:s%obj%nstring)%x(i)))
       end do
    else
       xmin = 0._c_float
       xmax = 1._c_float
    end if

    ! new scene center and center shift in xc
    xc = s%scenecenter
    s%scenecenter = 0.5_c_float * (xmin + xmax)
    xc = s%scenecenter - xc

    ! radius and bounding box
    s%scenerad = 0.5_c_float * norm2(xmax - xmin)
    s%scenexmin = xmin
    s%scenexmax = xmax

    ! translate the scene so the center position remains unchanged
    if (s%iscaminit .and. .not.s%nextbuildlists_fixcam) &
       call translate(s%world,-xc)
    s%nextbuildlists_fixcam = .false.

    ! Now that the scene radius is known, build the deferred window-anchored
    ! axes. Auto-size the gizmo from the scene radius so it occupies a roughly
    ! constant fraction of the window regardless of the system size (and of how
    ! many cells are drawn). This is done only once, the first time the lists
    ! are built (or after the gizmo is re-anchored); the auto flag is then
    ! cleared so manual edits and later rebuilds keep the value.
    do i = 1, s%nrep
       if (s%rep(i)%type == reptype_axes .and. s%rep(i)%axes%placement == 1) then
          if (s%rep(i)%axes%scale_auto) then
             s%rep(i)%axes%scale = axes_winfrac_def * real(s%scenerad,8) / max(s%rep(i)%axes%length,1d-10)
             s%rep(i)%axes%scale_auto = .false.
          end if
          call s%rep(i)%add_draw_elements(s%nc,s%obj,s%animation>0,s%iqpt_selected,s%ifreq_selected)
       end if
    end do

    ! persistent symmetry-element sets: deferred like anchored axes because they
    ! are sized from the scene radius (and excluded from the bounding box)
    do i = 1, s%nrep
       if (s%rep(i)%type == reptype_symelem) then
          s%rep(i)%symelem%size = real(s%scenerad,8)
          s%rep(i)%symelem%cen = real(s%scenecenter,8)
          call s%rep(i)%add_draw_elements(s%nc,s%obj,s%animation>0,s%iqpt_selected,s%ifreq_selected)
       end if
    end do

    ! transient representations (e.g. hover hints): built last, after the scene
    ! bounding box is known, so they never perturb the camera or scene size
    do i = 1, s%nreptrans
       call s%reptrans(i)%add_draw_elements(s%nc,s%obj,s%animation>0,s%iqpt_selected,s%ifreq_selected)
    end do

    ! flag whether any object is anchored to the window borders (the view
    ! window uses this to re-render when the window geometry changes)
    s%hasanchoredobj = (s%obj%ncylgiz > 0 .or. s%obj%nconegiz > 0 .or. s%obj%nstringgiz > 0)

    ! rebuilding lists is done; the cached instance buffers are now stale and
    ! must be repacked/uploaded on the next render
    s%forcebuildlists = .false.
    s%gl%inst_valid = .false.
    s%isinit = 2
    s%timelastbuild = glfwGetTime()

  end subroutine scene_build_lists

  !> Draw the scene.
  module subroutine scene_render(s)
    use interfaces_glfw, only: glfwGetTime
    use interfaces_cimgui
    use interfaces_opengl3
    use shapes, only: textVAOos, textVBOos, quadnel, trinel,&
       sph_inst_nf, cyl_inst_nf, mesh_inst_nf, text_vert_nf, connel, nmaxcone,&
       glb_cone, glb_plane, glb_tri, glb_conescr, ensure_pack
    use gui_main, only: fonts, fontbakesize_large, font_large
    use systems, only: sys, sysc
    use tools_io, only: string
    use shaders, only: shader_text_onscene, shader_sphere, shader_cylinder,&
       shader_mesh, useshader, setuniform_int, setuniform_float, setuniform_vec3,&
       setuniform_mat4,&
       uniloc, u_world, u_view, u_projection, u_isortho, u_displ, u_upick,&
       u_isanchored, u_anchored_ndc, u_anchored_scale, u_textcolor
    use param, only: img, pi
    class(scene), intent(inout), target :: s

    real(c_float) :: xsel(3,4), radsel(4)
    complex(c_float_complex) :: displ
    real*8 :: deltat, fac, time
    logical :: doit, dobuild

    real(c_float), parameter :: msel_thickness = 0.1_c_float
    real(c_float), parameter :: sel_thickness = 0.2_c_float
    real(c_float), parameter :: sel_label_size = 1.2_c_float
    real*8, parameter :: freq_ref = 300d0
    real*8, parameter :: freq_min = 50d0

    ! check that the scene and system are initialized
    if (s%isinit == 0) return

    ! buffers, draw lists, camera lock, camera reset
    call scene_render_prepare(s)

    ! get the time
    time = glfwGetTime()

    ! render text with the large font
    call igPushFont(font_large)

    ! calculate the time factor. The phasor is zero unless an animation is
    ! actually running, so Animation=None shows the equilibrium geometry (the
    ! GPU applies displ to the per-instance vibration deltas).
    displ = 0._c_float
    if (s%ifreq_selected > 0 .and. s%iqpt_selected > 0 .and. s%animation > 0) then
       fac = s%anim_amplitude * sqrt(freq_ref / max(abs(sys(s%id)%c%vib%freq(s%ifreq_selected,s%iqpt_selected)),freq_min))
       if (s%animation == 1) then ! manual
          displ = cmplx(fac * exp(0.5d0 * s%anim_phase * pi * img),&
             kind=c_float_complex)
       else ! automatic
          deltat = time - s%timerefanimation
          displ = cmplx(fac * exp(deltat * s%anim_speed * img),kind=c_float_complex)
       end if
    end if

    ! decide whether this scene's cached instance buffers must be (re)built
    ! and uploaded: only when they are stale (draw lists changed). All
    ! vibration animation is applied in the vertex shaders from the
    ! per-instance deltas and the displ uniform, so the buffers (and the text
    ! cache) are animation-invariant.
    dobuild = .not.s%gl%inst_valid

    ! draw the atoms (instanced sphere impostors)
    if (s%obj%nsph > 0) then
       call setup_shader(shader_sphere)
       if (dobuild) then
          call draw_all_spheres()
       else
          call s%gl%redraw_spheres(s%gl%nsph_inst)
       end if
    end if

    ! draw the bonds and the unit-cell edges (instanced cylinder impostors).
    ! Blending is enabled because the bond alpha is plumbed through to the
    ! fragment color (translucent bonds); all cylinders are opaque (alpha=1)
    ! at present, so this has no visible effect today.
    if (s%obj%ncyl + s%obj%ncylflat > 0) then
       call setup_shader(shader_cylinder)
       call glEnable(GL_BLEND)
       call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
       if (dobuild) then
          call draw_all_cylinders()
       else
          call s%gl%redraw_cylinders(s%gl%ncyl_inst)
       end if
       call glDisable(GL_BLEND)
    end if

    ! draw the plain meshes (cones, planes, polyhedra triangles)
    if (s%obj%ncone + s%obj%nplane + s%obj%ntriangle > 0) then
       call setup_shader(shader_mesh)

       ! cones (arrowheads): opaque
       if (s%obj%ncone > 0) then
          if (dobuild) then
             call draw_all_cones()
          else
             call s%gl%redraw_mesh(glb_cone,connel(nmaxcone),s%gl%ncone_inst)
          end if
       end if

       ! flat rectangles and polyhedra triangles (translucent: keep the depth
       ! test but disable depth writes so atoms/labels behind stay visible)
       if (s%obj%nplane + s%obj%ntriangle > 0) then
          call glDisable(GL_CULL_FACE)
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          call glDepthMask(int(GL_FALSE,c_signed_char))
          if (s%obj%nplane > 0) then
             if (dobuild) then
                call draw_all_planes()
             else
                call s%gl%redraw_mesh(glb_plane,quadnel,s%gl%nplane_inst)
             end if
          end if
          if (s%obj%ntriangle > 0) then
             if (dobuild) then
                call draw_all_triangles()
             else
                call s%gl%redraw_mesh(glb_tri,trinel,s%gl%ntri_inst)
             end if
          end if
          call glDepthMask(int(GL_TRUE,c_signed_char))
          call glDisable(GL_BLEND)
          call glEnable(GL_CULL_FACE)
       end if
    end if

    ! this scene's cached instance buffers are now current
    s%gl%inst_valid = .true.

    ! draw the measure selection atoms (instanced sphere impostors)
    if (s%nmsel > 0) then
       call setup_shader(shader_sphere)
       call glEnable(GL_BLEND)
       call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
       call draw_all_mselections()
       call glDisable(GL_BLEND)
    end if

    ! render labels with on-scene text
    call setup_shader(shader_text_onscene)

    call glDisable(GL_MULTISAMPLE)
    call glEnable(GL_BLEND)
    call glBlendEquation(GL_FUNC_ADD)
    call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)

    call glActiveTexture(GL_TEXTURE0)
    call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))

    ! on-scene text quads have reversed winding (y-down glyphs flipped for the
    ! right-side-up presentation), so disable culling while drawing them
    call glDisable(GL_CULL_FACE)
    call draw_all_text()

    ! render selected atom labels with on-scene text
    if (s%nmsel > 0) &
       call draw_selection_text()
    call glEnable(GL_CULL_FACE)
    call glEnable(GL_MULTISAMPLE)
    call glDisable(GL_BLEND)

    ! highlight the highlighted/selected atoms (translucent overlay spheres;
    ! multisampling is back on so their rims are antialiased like the atoms)
    doit = .false.
    if (allocated(sysc(s%id)%highlight_rgba)) &
       doit = any(sysc(s%id)%highlight_rgba >= 0._c_float)
    if (.not.doit.and.allocated(sysc(s%id)%highlight_rgba_transient)) &
       doit = any(sysc(s%id)%highlight_rgba_transient >= 0._c_float)
    if (doit) then
       call setup_shader(shader_sphere)
       call glEnable(GL_BLEND)
       call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
       call draw_highlights()
       call glDisable(GL_BLEND)
    end if

    ! window-anchored axes gizmo, drawn on top of the scene
    if (s%obj%ncylgiz + s%obj%nconegiz + s%obj%nstringgiz > 0) &
       call render_axes_gizmo()

    ! pop the large font
    call igPopFont()

    ! clean up
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! save the rendering time
    s%timelastrender = time

  contains
    !> Bind shader ishader and set the shared scene uniforms (world, view,
    !> projection, isanchored=0, and the vibration displacement phasor; the
    !> u_displ set is a no-op for programs without the uniform). The impostor
    !> shaders also get the projection type; the sphere shader additionally
    !> gets pick mode off.
    subroutine setup_shader(ishader)
      integer, intent(in) :: ishader

      call useshader(ishader)
      call setuniform_mat4(s%world,idxi=uniloc(u_world))
      call setuniform_mat4(s%view,idxi=uniloc(u_view))
      call setuniform_mat4(s%projection,idxi=uniloc(u_projection))
      call setuniform_int(0_c_int,idxi=uniloc(u_isanchored))
      call setuniform_vec3((/real(displ,c_float),real(aimag(displ),c_float),0._c_float/),idxi=uniloc(u_displ))
      if (ishader == shader_sphere .or. ishader == shader_cylinder) &
         call setuniform_int(merge(1_c_int,0_c_int,s%isortho),idxi=uniloc(u_isortho))
      if (ishader == shader_sphere) &
         call setuniform_int(0_c_int,idxi=uniloc(u_upick))

    end subroutine setup_shader

    subroutine draw_all_spheres()
      integer :: i, n
      real(c_float), parameter :: zr(4) = 0._c_float

      call ensure_pack(s%gl%packsph,sph_inst_nf,s%obj%nsph)
      n = 0
      do i = 1, s%obj%nsph
         if (s%obj%sph(i)%ghost) cycle ! invisible pick-only target, not drawn
         n = n + 1
         call sphere_pack(s%gl%packsph(:,n),s%obj%sph(i)%x,s%obj%sph(i)%r,&
            (/s%obj%sph(i)%rgb,1._c_float/),s%obj%sph(i)%border,&
            s%obj%sph(i)%rgbborder,s%obj%sph(i)%xdelta,zr,s%obj%sph(i)%occ,&
            s%obj%sph(i)%occ_empty_rgb)
      end do
      call s%gl%draw_spheres(n,s%gl%packsph,.false.)
      s%gl%nsph_inst = n

    end subroutine draw_all_spheres

    !> Build and draw all cylinder impostors: bonds (multi-bonds
    !> expanded into laterally-offset instances; dashed and the
    !> aromatic interior expanded into short capped-cylinder dashes)
    !> and the unit-cell edges (plain, no border).
    subroutine draw_all_cylinders()
      integer :: i, n, k, j, nseg, ntot, ndash
      real(c_float) :: x1(3), x2(3), rgba(4), r, rgeo, bcol(3), bord, outw(3)
      real(c_float) :: dl(3), dir(3), blen, p, half, tc, sx1(3), sx2(3), t1, t2
      complex(c_float_complex) :: xd1(3), xd2(3), sd1(3), sd2(3)
      logical :: dsh(3)
      real(c_float), parameter :: zv(3) = 0._c_float

      ! count the instances (dashes expand into several short cylinders)
      ntot = s%obj%ncylflat
      do i = 1, s%obj%ncyl
         select case (s%obj%cyl(i)%order)
         case (-1)
            ntot = ntot + 1 + dash_count(i)
         case (0)
            ntot = ntot + dash_count(i)
         case (2)
            ntot = ntot + 2
         case (3)
            ntot = ntot + 3
         case default
            ntot = ntot + 1
         end select
      end do

      call ensure_pack(s%gl%packcyl,cyl_inst_nf,ntot)
      n = 0

      ! bonds: expand the bond order into one or more (possibly dashed) runs
      do i = 1, s%obj%ncyl
         x1 = s%obj%cyl(i)%x1
         x2 = s%obj%cyl(i)%x2
         xd1 = s%obj%cyl(i)%x1delta
         xd2 = s%obj%cyl(i)%x2delta
         r = s%obj%cyl(i)%r
         rgeo = 0.5_c_float * r
         rgba = (/s%obj%cyl(i)%rgb, s%obj%cyl(i)%alpha/)
         bord = s%obj%cyl(i)%border
         bcol = s%obj%cyl(i)%rgbborder
         outw = 0._c_float
         dsh = .false.
         select case (s%obj%cyl(i)%order)
         case (-1) ! aromatic (1.5): solid ring exterior + dashed ring interior
            nseg = 2
            dl(1) = 0.75_c_float*r
            dl(2) = -0.75_c_float*r; dsh(2) = .true.
            outw = s%obj%cyl(i)%arvec
         case (0) ! dashed
            nseg = 1
            dl(1) = 0._c_float; dsh(1) = .true.
         case (2) ! double
            nseg = 2
            dl(1) = 0.75_c_float*r; dl(2) = -0.75_c_float*r
         case (3) ! triple
            nseg = 3
            dl(1) = 1.35_c_float*r; dl(2) = 0._c_float; dl(3) = -1.35_c_float*r
         case default ! single (1) or flat (< -1)
            nseg = 1
            dl(1) = 0._c_float
         end select

         blen = norm2(x2-x1)
         dir = 0._c_float
         if (blen > 1e-7_c_float) dir = (x2-x1)/blen

         do k = 1, nseg
            if (dsh(k)) then
               ! string of short capped cylinders (dashes); each dash carries
               ! the endpoint deltas interpolated at its fractional positions,
               ! so the dashes ride the animated bond
               ndash = dash_count(i)
               p = blen / real(ndash,c_float)
               half = 0.25_c_float * p
               do j = 0, ndash-1
                  tc = (real(j,c_float) + 0.5_c_float) * p
                  sx1 = x1 + dir * (tc - half)
                  sx2 = x1 + dir * (tc + half)
                  t1 = (tc - half) / max(blen,1e-7_c_float)
                  t2 = (tc + half) / max(blen,1e-7_c_float)
                  sd1 = (1._c_float - t1) * xd1 + t1 * xd2
                  sd2 = (1._c_float - t2) * xd1 + t2 * xd2
                  n = n + 1
                  call cyl_pack(s%gl%packcyl(:,n),sx1,sx2,rgeo,rgba,bord,bcol,dl(k),outw,sd1,sd2)
               end do
            else
               n = n + 1
               call cyl_pack(s%gl%packcyl(:,n),x1,x2,rgeo,rgba,bord,bcol,dl(k),outw,xd1,xd2)
            end if
         end do
      end do

      ! unit-cell edges: plain cylinders, no border
      do i = 1, s%obj%ncylflat
         n = n + 1
         call cyl_pack(s%gl%packcyl(:,n),s%obj%cylflat(i)%x1,s%obj%cylflat(i)%x2,&
            0.5_c_float*s%obj%cylflat(i)%r,(/s%obj%cylflat(i)%rgb,1._c_float/),&
            0._c_float,zv,0._c_float,zv,s%obj%cylflat(i)%x1delta,s%obj%cylflat(i)%x2delta)
      end do

      call s%gl%draw_cylinders(n,s%gl%packcyl,.false.)
      s%gl%ncyl_inst = n

    end subroutine draw_all_cylinders

    !> Number of dashes for the dashed expansion of bond cylinder i, from the
    !> equilibrium bond length (frozen while animating).
    function dash_count(i) result(ndash)
      integer, intent(in) :: i
      integer :: ndash

      real(c_float), parameter :: dash_period = 0.4_c_float

      ndash = min(64, max(1, nint(norm2(s%obj%cyl(i)%x2 - s%obj%cyl(i)%x1)/dash_period)))

    end function dash_count

    subroutine draw_all_cones()
      integer :: i, n
      real(c_float) :: model(4,4)

      call ensure_pack(s%gl%packmesh,mesh_inst_nf,s%obj%ncone)
      n = 0
      do i = 1, s%obj%ncone
         if (norm2(s%obj%cone(i)%x2 - s%obj%cone(i)%x1) < 1e-4_c_float) cycle
         call cone_model(s%obj%cone(i)%x1,s%obj%cone(i)%x2,s%obj%cone(i)%r,model)
         n = n + 1
         call mesh_pack(s%gl%packmesh(:,n),model,(/s%obj%cone(i)%rgb,1._c_float/))
      end do
      call s%gl%draw_mesh(glb_cone,connel(nmaxcone),n,s%gl%packmesh)
      s%gl%ncone_inst = n

    end subroutine draw_all_cones

    !> Target NDC position (gizndc) and zoom-compensation factor (gizf) for a
    !> gizmo item at the window fraction winpos. The shader (isanchored) anchors the
    !> local-origin geometry at gizndc and scales it by gizf — no matrix inverse
    !> or per-item world matrix needed.
    subroutine gizmo_ndc_scale(winpos,scalewithzoom,projgiz,gizndc,gizf)
      real(c_float), intent(in) :: winpos(2)
      logical, intent(in) :: scalewithzoom
      real(c_float), intent(in) :: projgiz(4,4)
      real(c_float), intent(out) :: gizndc(3), gizf

      real(c_float) :: spanx, spany, hside

      ! winpos is given as fractions of the visible part of the render buffer
      ! (the window), from the left and bottom; map it to the NDC of the full
      ! render texture, accounting for the cropped region (viewuv0). The texture
      ! is presented right-side up, so the bottom of the window is NDC y = -1.
      spanx = 1._c_float - 2._c_float * s%viewuv0(1)
      spany = 1._c_float - 2._c_float * s%viewuv0(2)
      gizndc(1) = spanx * (2._c_float * winpos(1) - 1._c_float)
      gizndc(2) = spany * (2._c_float * winpos(2) - 1._c_float)
      gizndc(3) = 0._c_float

      ! zoom-compensation factor: when the gizmo should not scale with zoom,
      ! shrink/grow the geometry so the orthographic projection leaves its
      ! on-screen size constant. hside is the half-window size at the reset zoom,
      ! so at that zoom both modes coincide (no jump when toggling).
      if (scalewithzoom) then
         gizf = 1._c_float
      else
         hside = reset_zoom_hside(s)
         gizf = 1._c_float / (projgiz(1,1) * hside)
      end if

    end subroutine gizmo_ndc_scale

    !> Render the window-anchored axes gizmo(s). Each gizmo item carries its own
    !> window position and zoom-scale flag; the object shaders (isanchored) anchor
    !> the local-origin geometry at the requested NDC position, drawn on top with
    !> an orthographic projection. Shafts use the cylinder impostor, arrowheads
    !> the mesh shader, and labels the on-scene text shader. Consecutive items
    !> with the same window placement (in practice, the shafts/heads of one
    !> gizmo) are batched into a single instanced draw.
    subroutine render_axes_gizmo()
      real(c_float) :: projgiz(4,4), gizndc(3), gizf
      real(c_float) :: hside, siz, model(4,4)
      integer :: i, j, k, n
      integer(c_int) :: nvert
      real(c_float), allocatable, target :: vert(:,:)
      real(c_float), parameter :: zv(3) = 0._c_float
      complex(c_float_complex), parameter :: zc(3) = (0._c_float,0._c_float)

      ! the gizmo always uses an orthographic projection, with a symmetric depth
      ! range so the anchored (eye-origin-centered) geometry is never clipped
      call ortho_projection(s,projgiz,symz=.true.)

      ! clear the depth buffer so the gizmo always draws on top
      call glClear(GL_DEPTH_BUFFER_BIT)

      ! shafts: cylinder impostors, anchored, batched by window placement
      if (s%obj%ncylgiz > 0) then
         call useshader(shader_cylinder)
         call setuniform_mat4(s%world,idxi=uniloc(u_world))
         call setuniform_mat4(s%view,idxi=uniloc(u_view))
         call setuniform_mat4(projgiz,idxi=uniloc(u_projection))
         call setuniform_int(1_c_int,idxi=uniloc(u_isortho))
         call setuniform_int(1_c_int,idxi=uniloc(u_isanchored))
         call ensure_pack(s%gl%packcyl,cyl_inst_nf,s%obj%ncylgiz)
         i = 1
         do while (i <= s%obj%ncylgiz)
            j = giz_group_end(s%obj%cylgiz,s%obj%ncylgiz,i)
            call gizmo_ndc_scale(s%obj%cylgiz(i)%winpos,s%obj%cylgiz(i)%scalewithzoom,projgiz,gizndc,gizf)
            call setuniform_vec3(gizndc,idxi=uniloc(u_anchored_ndc))
            call setuniform_float(gizf,idxi=uniloc(u_anchored_scale))
            n = 0
            do k = i, j
               n = n + 1
               call cyl_pack(s%gl%packcyl(:,n),s%obj%cylgiz(k)%x1,s%obj%cylgiz(k)%x2,&
                  0.5_c_float*s%obj%cylgiz(k)%r,(/s%obj%cylgiz(k)%rgb,1._c_float/),&
                  s%obj%cylgiz(k)%border,s%obj%cylgiz(k)%rgbborder,0._c_float,zv,zc,zc)
            end do
            call s%gl%draw_cylinders(n,s%gl%packcyl,.true.)
            i = j + 1
         end do
      end if

      ! arrowheads: cone meshes, anchored, batched by window placement
      if (s%obj%nconegiz > 0) then
         call useshader(shader_mesh)
         call setuniform_mat4(s%world,idxi=uniloc(u_world))
         call setuniform_mat4(s%view,idxi=uniloc(u_view))
         call setuniform_mat4(projgiz,idxi=uniloc(u_projection))
         call setuniform_int(1_c_int,idxi=uniloc(u_isanchored))
         call ensure_pack(s%gl%packmesh,mesh_inst_nf,s%obj%nconegiz)
         i = 1
         do while (i <= s%obj%nconegiz)
            j = giz_group_end(s%obj%conegiz,s%obj%nconegiz,i)
            call gizmo_ndc_scale(s%obj%conegiz(i)%winpos,s%obj%conegiz(i)%scalewithzoom,projgiz,gizndc,gizf)
            call setuniform_vec3(gizndc,idxi=uniloc(u_anchored_ndc))
            call setuniform_float(gizf,idxi=uniloc(u_anchored_scale))
            n = 0
            do k = i, j
               if (norm2(s%obj%conegiz(k)%x2 - s%obj%conegiz(k)%x1) < 1e-4_c_float) cycle
               call cone_model(s%obj%conegiz(k)%x1,s%obj%conegiz(k)%x2,s%obj%conegiz(k)%r,model)
               n = n + 1
               call mesh_pack(s%gl%packmesh(:,n),model,(/s%obj%conegiz(k)%rgb,1._c_float/))
            end do
            call s%gl%draw_mesh(glb_conescr,connel(nmaxcone),n,s%gl%packmesh)
            i = j + 1
         end do
      end if

      ! labels: streamed per string (each carries its own color and placement)
      if (s%obj%nstringgiz > 0) then
         call useshader(shader_text_onscene)
         call setuniform_int(1_c_int,idxi=uniloc(u_isanchored))
         call setuniform_mat4(s%world,idxi=uniloc(u_world))
         call setuniform_mat4(s%view,idxi=uniloc(u_view))
         call setuniform_mat4(projgiz,idxi=uniloc(u_projection))
         call glDisable(GL_MULTISAMPLE)
         call glDisable(GL_CULL_FACE) ! y-flipped glyph quads have reversed winding
         call glEnable(GL_BLEND)
         call glBlendEquation(GL_FUNC_ADD)
         call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)
         call glActiveTexture(GL_TEXTURE0)
         call glBindVertexArray(textVAOos)
         call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))
         call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)
         do i = 1, s%obj%nstringgiz
            call gizmo_ndc_scale(s%obj%stringgiz(i)%winpos,s%obj%stringgiz(i)%scalewithzoom,projgiz,gizndc,gizf)
            call setuniform_vec3(gizndc,idxi=uniloc(u_anchored_ndc))
            call setuniform_float(gizf,idxi=uniloc(u_anchored_scale))
            call setuniform_vec3(s%obj%stringgiz(i)%rgb,idxi=uniloc(u_textcolor))
            nvert = 0
            ! the label tracks the same zoom behavior as the arrows: constant
            ! on-screen size when the gizmo does not scale with zoom, or scaling
            ! with the gizmo's orthographic projection otherwise.
            if (.not.s%obj%stringgiz(i)%scalewithzoom) then
               hside = reset_zoom_hside(s)
               siz = 2 * abs(s%obj%stringgiz(i)%scale) / fontbakesize_large / hside
            else
               siz = 2 * abs(s%obj%stringgiz(i)%scale) * projgiz(1,1) / fontbakesize_large
            end if
            call calc_text_onscene_vertices(s%obj%stringgiz(i)%str,s%obj%stringgiz(i)%x,s%obj%stringgiz(i)%r,&
               siz,nvert,vert,shift=s%obj%stringgiz(i)%offset,centered=.true.)
            call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*text_vert_nf*c_sizeof(c_float), c_loc(vert))
            call glDrawArrays(GL_TRIANGLES, 0, nvert)
         end do
         call glEnable(GL_MULTISAMPLE)
         call glEnable(GL_CULL_FACE)
         call glDisable(GL_BLEND)
         call setuniform_int(0_c_int,idxi=uniloc(u_isanchored))
      end if

    end subroutine render_axes_gizmo

    !> Last index of the run of gizmo items starting at i that share the same
    !> window placement (winpos and scalewithzoom) as item i.
    function giz_group_end(arr,ntot,i) result(j)
      type(dl_cylinder_giz), intent(in) :: arr(:)
      integer, intent(in) :: ntot, i
      integer :: j

      j = i
      do while (j < ntot)
         if (arr(j+1)%winpos(1) /= arr(i)%winpos(1) .or. arr(j+1)%winpos(2) /= arr(i)%winpos(2) .or. &
            (arr(j+1)%scalewithzoom .neqv. arr(i)%scalewithzoom)) exit
         j = j + 1
      end do

    end function giz_group_end

    ! Draw all filled rectangles. The unit quad (corners +/-1,+/-1,0)
    ! is mapped onto each rectangle by a model matrix whose first two
    ! columns are the in-plane half-edge vectors and whose translation
    ! is the rectangle center. The bound shader/uniform state must
    ! already render flat (object_type=2 in the simple shader).
    subroutine draw_all_planes()
      use tools_math, only: cross_cfloat
      integer :: i, n
      real(c_float) :: m(4,4), nrm(3)

      call ensure_pack(s%gl%packmesh,mesh_inst_nf,s%obj%nplane)
      n = 0
      do i = 1, s%obj%nplane
         nrm = cross_cfloat(s%obj%plane(i)%e1,s%obj%plane(i)%e2)
         if (norm2(nrm) > 1e-10_c_float) nrm = nrm / norm2(nrm)
         m = 0._c_float
         m(4,4) = 1._c_float
         m(1:3,1) = s%obj%plane(i)%e1
         m(1:3,2) = s%obj%plane(i)%e2
         m(1:3,3) = nrm
         m(1:3,4) = s%obj%plane(i)%x
         n = n + 1
         call mesh_pack(s%gl%packmesh(:,n),m,(/s%obj%plane(i)%rgb,s%obj%plane(i)%alpha/))
      end do
      call s%gl%draw_mesh(glb_plane,quadnel,n,s%gl%packmesh)
      s%gl%nplane_inst = n

    end subroutine draw_all_planes

    ! Draw all filled triangles. The unit reference triangle (corners
    ! (0,0,0),(1,0,0),(0,1,0)) is mapped onto each triangle by a model matrix
    ! whose first two columns are (x2-x1) and (x3-x1) and whose translation is
    ! x1. The bound shader/uniform state must already render flat (object_type=2).
    subroutine draw_all_triangles()
      use tools_math, only: cross_cfloat
      integer :: i, n
      real(c_float) :: m(4,4), nrm(3), p1(3), p2(3), p3(3)

      ! the model matrix maps the reference triangle onto the equilibrium
      ! vertices; the vibration displacement is applied in the vertex shader
      ! from the per-corner deltas
      call ensure_pack(s%gl%packmesh,mesh_inst_nf,s%obj%ntriangle)
      n = 0
      do i = 1, s%obj%ntriangle
         p1 = s%obj%triangle(i)%x1
         p2 = s%obj%triangle(i)%x2
         p3 = s%obj%triangle(i)%x3
         nrm = cross_cfloat(p2 - p1,p3 - p1)
         if (norm2(nrm) > 1e-10_c_float) nrm = nrm / norm2(nrm)
         m = 0._c_float
         m(4,4) = 1._c_float
         m(1:3,1) = p2 - p1
         m(1:3,2) = p3 - p1
         m(1:3,3) = nrm
         m(1:3,4) = p1
         n = n + 1
         call mesh_pack(s%gl%packmesh(:,n),m,(/s%obj%triangle(i)%rgb,s%obj%triangle(i)%alpha/),&
            s%obj%triangle(i)%x1delta,s%obj%triangle(i)%x2delta,s%obj%triangle(i)%x3delta)
      end do
      call s%gl%draw_mesh(glb_tri,trinel,n,s%gl%packmesh)
      s%gl%ntri_inst = n

    end subroutine draw_all_triangles

    !> Draw the measure selections
    subroutine draw_all_mselections()
      use gui_main, only: ColorMeasureSelect
      use representations, only: atomborder_def
      integer :: i, j
      real(c_float) :: x(3)
      real(c_float), parameter :: zr(4) = 0._c_float

      call ensure_pack(s%gl%packsph,sph_inst_nf,s%nmsel)
      do j = 1, s%nmsel
         i = s%msel(5,j)
         ! CPU-displaced label anchor for draw_selection_text (displ is zero
         ! when no animation is running)
         x = s%obj%sph(i)%x + real(displ * s%obj%sph(i)%xdelta,c_float)
         call sphere_pack(s%gl%packsph(:,j),s%obj%sph(i)%x,s%obj%sph(i)%r + msel_thickness,&
            ColorMeasureSelect(:,j),real(atomborder_def,c_float),(/0._c_float,0._c_float,0._c_float/),&
            s%obj%sph(i)%xdelta,zr,1._c_float,(/0._c_float,0._c_float,0._c_float/))
         radsel(j) = s%obj%sph(i)%r + msel_thickness
         xsel(:,j) = x
      end do
      call s%gl%draw_spheres(s%nmsel,s%gl%packsph,.true.)

    end subroutine draw_all_mselections

    !> Draw the highlights on the scene
    subroutine draw_highlights()
      use systems, only: sysc
      use representations, only: atomborder_def
      integer :: i, id, n
      real(c_float) :: rgba(4)
      real(c_float), parameter :: zr(4) = 0._c_float

      ! initial checks
      if (s%isinit < 2) return
      if (.not.allocated(sysc(s%id)%highlight_rgba).and.&
         .not.allocated(sysc(s%id)%highlight_rgba_transient)) return

      ! highlight the spheres
      call ensure_pack(s%gl%packsph,sph_inst_nf,s%obj%nsph)
      n = 0
      do i = 1, s%obj%nsph
         id = s%obj%sph(i)%idx(1)
         ! skip any non-atom sphere (idx < 1): it carries no cell-atom index and
         ! is not present in the per-atom highlight arrays
         if (id < 1) cycle
         rgba = -1._c_float
         if (allocated(sysc(s%id)%highlight_rgba_transient)) &
            rgba = sysc(s%id)%highlight_rgba_transient(:,id)
         if (any(rgba < 0).and.allocated(sysc(s%id)%highlight_rgba)) &
            rgba = sysc(s%id)%highlight_rgba(:,id)
         if (all(rgba >= 0)) then
            n = n + 1
            call sphere_pack(s%gl%packsph(:,n),s%obj%sph(i)%x,s%obj%sph(i)%r + sel_thickness,&
               rgba,real(atomborder_def,c_float),rgba(1:3),s%obj%sph(i)%xdelta,zr,1._c_float,&
               (/0._c_float,0._c_float,0._c_float/))
         end if
      end do
      call s%gl%draw_spheres(n,s%gl%packsph,.true.)

    end subroutine draw_highlights

    !> Draw the scene labels. The glyph vertices are cached in the per-scene
    !> text buffer and rebuilt only when their inputs change: the draw lists
    !> (timelastbuild) or any camera/projection quantity entering the
    !> on-screen glyph size (projection matrix, view*world third row,
    !> reset-zoom half-side, baked font size).
    subroutine draw_all_text()
      integer :: i, nvert, nv0
      real(c_float) :: hside, siz, x(3), vw(4,4), wclip
      logical :: rebuild

      if (s%obj%nstring <= 0) return

      ! view*world, used to get the anchor's clip-space w (= 1 in orthographic,
      ! = view depth in perspective) for projection-aware label sizing
      vw = matmul(s%view,s%world)

      ! half-window size at the reset zoom (constant on-screen size labels)
      hside = reset_zoom_hside(s)

      ! decide whether the cached glyph vertices are still valid
      rebuild = .not.s%gl%text_valid
      if (.not.rebuild) rebuild = (s%gl%text_build_time /= s%timelastbuild)
      if (.not.rebuild) rebuild = any(s%gl%text_proj /= s%projection)
      if (.not.rebuild) rebuild = any(s%gl%text_vw3 /= vw(3,:))
      if (.not.rebuild) rebuild = (s%gl%text_hside /= hside) .or.&
         (s%gl%text_fontsize /= fontbakesize_large)

      if (rebuild) then
         ! rebuild the concatenated glyph vertices for all labels and upload
         ! them once
         call ensure_pack(s%gl%text_first,s%obj%nstring)
         call ensure_pack(s%gl%text_count,s%obj%nstring)
         nvert = 0
         do i = 1, s%obj%nstring
            nv0 = nvert
            x = s%obj%string(i)%x
            if (s%obj%string(i)%scale > 0._c_float) then
               ! constant on-screen size (projection-independent)
               siz = 2 * s%obj%string(i)%scale / fontbakesize_large / hside
            else
               ! scale with zoom (projection-aware): divide by the anchor
               ! clip-space w so the label foreshortens with depth under
               ! perspective
               wclip = s%projection(4,3) * (vw(3,1)*x(1)+vw(3,2)*x(2)+vw(3,3)*x(3)+vw(3,4)) + s%projection(4,4)
               wclip = max(wclip,1e-4_c_float) ! guard anchors at/behind the camera (perspective); =1 in ortho
               siz = 2 * abs(s%obj%string(i)%scale) * s%projection(1,1) / fontbakesize_large / wclip
            end if
            call calc_text_onscene_vertices(s%obj%string(i)%str,x,s%obj%string(i)%r,&
               siz,nvert,s%gl%packtext,shift=s%obj%string(i)%offset,centered=.true.,&
               xdelta=s%obj%string(i)%xdelta)
            s%gl%text_first(i) = nv0
            s%gl%text_count(i) = nvert - nv0
         end do
         call s%gl%upload_text(nvert,s%gl%packtext)

         ! stamp the validity keys
         s%gl%text_valid = .true.
         s%gl%text_build_time = s%timelastbuild
         s%gl%text_proj = s%projection
         s%gl%text_vw3 = vw(3,:)
         s%gl%text_hside = hside
         s%gl%text_fontsize = fontbakesize_large
      end if

      ! draw each label from the cached buffer with its own color
      call glBindVertexArray(s%gl%textVAO)
      do i = 1, s%obj%nstring
         if (s%gl%text_count(i) <= 0) cycle
         call setuniform_vec3(s%obj%string(i)%rgb,idxi=uniloc(u_textcolor))
         call glDrawArrays(GL_TRIANGLES, int(s%gl%text_first(i),c_int), int(s%gl%text_count(i),c_int))
      end do
      call glBindVertexArray(0)

    end subroutine draw_all_text

    !> Draw the measure-selection numerals. At most 4 short strings, streamed
    !> through the shared on-scene text buffer every frame (not cached).
    subroutine draw_selection_text()
      integer :: j
      real(c_float) :: siz, vw(4,4), wclip
      integer(c_int) :: nvert
      real(c_float), allocatable, target :: vert(:,:)

      ! view*world, for the anchor clip-space w (projection-aware sizing; see
      ! draw_all_text). These labels always scale with zoom.
      vw = matmul(s%view,s%world)

      call glBindVertexArray(textVAOos)
      call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)
      do j = 1, s%nmsel
         call setuniform_vec3((/1._c_float,1._c_float,1._c_float/),idxi=uniloc(u_textcolor))
         wclip = s%projection(4,3) * (vw(3,1)*xsel(1,j)+vw(3,2)*xsel(2,j)+vw(3,3)*xsel(3,j)+vw(3,4)) + s%projection(4,4)
         wclip = max(wclip,1e-4_c_float) ! guard anchors at/behind the camera (perspective); =1 in ortho
         siz = sel_label_size * s%projection(1,1) / fontbakesize_large / wclip
         nvert = 0
         ! the xsel anchors are already CPU-displaced (see draw_all_mselections);
         ! xdelta is deliberately omitted here, or the numerals would displace twice
         call calc_text_onscene_vertices(string(j),xsel(:,j),radsel(j),siz,nvert,vert,centered=.true.)
         call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*text_vert_nf*c_sizeof(c_float), c_loc(vert))
         call glDrawArrays(GL_TRIANGLES, 0, nvert)
      end do
      call glBindBuffer(GL_ARRAY_BUFFER, 0)
      call glBindVertexArray(0)

    end subroutine draw_selection_text

  end subroutine scene_render

  !> Draw the scene (for object picking)
  module subroutine scene_render_pick(s)
    use interfaces_cimgui
    use interfaces_opengl3
    use systems, only: sys
    use shapes, only: sph_inst_nf, ensure_pack
    use shaders, only: shader_sphere, useshader, setuniform_int,&
       setuniform_vec3, setuniform_mat4, uniloc, u_world, u_view,&
       u_projection, u_isortho, u_upick, u_isanchored, u_displ
    use param, only: maxzat0
    class(scene), intent(inout), target :: s

    integer :: i, iz, idx, n
    real(c_float) :: ridx(4)

    ! check that the scene and system are initialized
    if (s%isinit < 2) return

    ! buffers, draw lists, camera lock, camera reset
    call scene_render_prepare(s)

    ! blending must be off: the pick colors are bit-packed indices whose alpha
    ! channel is 0.0, so any leaked blend state would zero them out
    call glDisable(GL_BLEND)

    ! set up the sphere impostor shader for picking (ray-cast so the pickable
    ! silhouette matches the visible sphere); positions are not animated
    call useshader(shader_sphere)
    call setuniform_mat4(s%world,idxi=uniloc(u_world))
    call setuniform_mat4(s%view,idxi=uniloc(u_view))
    call setuniform_mat4(s%projection,idxi=uniloc(u_projection))
    call setuniform_int(merge(1_c_int,0_c_int,s%isortho),idxi=uniloc(u_isortho))
    call setuniform_int(1_c_int,idxi=uniloc(u_upick))
    call setuniform_int(0_c_int,idxi=uniloc(u_isanchored))
    call setuniform_vec3((/0._c_float,0._c_float,0._c_float/),idxi=uniloc(u_displ))

    ! draw the atoms, each with its loop index encoded into the pick color
    if (s%obj%nsph > 0) then
       call ensure_pack(s%gl%packsph,sph_inst_nf,s%obj%nsph)
       n = 0
       do i = 1, s%obj%nsph
          ! draw the sphere, no gradient paths
          idx = s%obj%sph(i)%idx(1)
          ! skip any non-atom sphere (idx < 1): it carries no cell-atom index
          ! and must not be pickable as an atom
          if (idx < 1) cycle
          iz = sys(s%id)%c%spc(sys(s%id)%c%atcel(idx)%is)%z
          if (iz < maxzat0) then
             n = n + 1
             ridx = transfer((/i,0,0,0/),ridx)
             call sphere_pack(s%gl%packsph(:,n),s%obj%sph(i)%x,s%obj%sph(i)%r,&
                (/0._c_float,0._c_float,0._c_float,1._c_float/),0._c_float,&
                (/0._c_float,0._c_float,0._c_float/),s%obj%sph(i)%xdelta,ridx,1._c_float,&
                (/0._c_float,0._c_float,0._c_float/))
          end if
       end do
       call s%gl%draw_spheres(n,s%gl%packsph,.true.)
    end if

  end subroutine scene_render_pick

  !> Set the appearance defaults for the current scene.
  module subroutine scene_set_style_defaults(s)
    use gui_main, only: ColorSceneBg_def
    class(scene), intent(inout), target :: s

    integer :: i

    s%bgcolor = ColorSceneBg_def
    do i = 1, s%nrep
       s%rep(i)%uc%rgb = 0._c_float
    end do

  end subroutine scene_set_style_defaults

  !> Copy camera parameters from scene si to the current scene. If an
  !> integer is given instead, search the system list for the most
  !> recent render in group idx and, if found, copy the camera
  !> parameters from that scene.
  module subroutine scene_cam_copy(s,si)
    use interfaces_glfw, only: glfwGetTime
    class(scene), intent(inout), target :: s
    type(scene), intent(in), target :: si

    s%camresetdist = si%camresetdist
    s%camratio = si%camratio
    s%isortho = si%isortho
    s%ortho_fov = si%ortho_fov
    s%persp_fov = si%persp_fov
    s%campos = si%campos
    s%camfront = si%camfront
    s%camup = si%camup
    s%world = si%world
    s%view = si%view
    s%projection = si%projection
    call s%update_view_matrix()
    call s%update_projection_matrix()
    s%timelastcamchange = glfwGetTime()
    s%forceresetcam = .false.
    s%iscaminit = .true.

  end subroutine scene_cam_copy

  !> Zoom in and out towards the center of the scene by a factor equal
  !> to ratio. min_zoom and max_zoom apply.
  module subroutine scene_cam_zoom(s,ratio)
    use interfaces_glfw, only: glfwGetTime
    use utils, only: mult
    class(scene), intent(inout), target :: s
    real(c_float), intent(in) :: ratio

    real(c_float) :: xc(3), pos3(3)

    ! calculate scene center in tworld coordinates
    call mult(xc,s%world,s%scenecenter)

    ! scale the vector from camera position to scene center
    pos3 = s%campos - xc
    pos3 = pos3 - ratio * pos3
    if (norm2(pos3) < min_zoom) &
       pos3 = pos3 / norm2(pos3) * min_zoom
    if (norm2(pos3) > max_zoom * s%scenerad) &
       pos3 = pos3 / norm2(pos3) * (max_zoom * s%scenerad)

    ! move the camera
    call s%cam_move(xc + pos3)

  end subroutine scene_cam_zoom

  !> Move camera to position xc (tworld coordinates).
  module subroutine scene_cam_move(s,xc)
    use interfaces_glfw, only: glfwGetTime
    use utils, only: mult
    class(scene), intent(inout), target :: s
    real(c_float), intent(in) :: xc(3)

    s%campos = xc
    call s%update_view_matrix()
    call s%update_projection_matrix()
    s%timelastcamchange = glfwGetTime()

  end subroutine scene_cam_move

  !> Rotate the camera around the given axis (view coordinates) by angle ang
  !> (radians).
  module subroutine scene_cam_rotate(s,axis,ang)
    use interfaces_glfw, only: glfwGetTime
    use utils, only: translate, rotate, invmult
    class(scene), intent(inout), target :: s
    real(c_float), intent(in) :: axis(3)
    real(c_float), intent(in) :: ang

    real(c_float) :: lax, axis_(3)

    lax = norm2(axis)
    if (lax > 1e-10_c_float) then
       axis_ = axis / lax
       call invmult(axis_,s%world,notrans=.true.)
       call translate(s%world,s%scenecenter)
       call rotate(s%world,ang,axis_)
       call translate(s%world,-s%scenecenter)
       s%timelastcamchange = glfwGetTime()
    end if

  end subroutine scene_cam_rotate

  !> Show the representation menu (called from view). Return .true.
  !> if the scene needs to be rendered again.
  module function representation_menu(s,idparent) result(changed)
    use interfaces_cimgui
    use representations, only: reptype_atoms, reptype_unitcell, reptype_axes, reptype_symelem
    use utils, only: iw_text, iw_tooltip, iw_button, iw_checkbox, iw_menuitem, iw_inputtext
    use windows, only: stack_create_window, wintype_editrep
    use gui_main, only: ColorDangerButton, g
    use tools_io, only: string
    use tools, only: mergesort
    class(scene), intent(inout), target :: s
    integer(c_int), intent(in) :: idparent
    logical :: changed

    integer :: i, ii, id, idum
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    logical :: discol, doerase, ok
    type(ImVec2) :: szero
    integer, allocatable :: idx(:)
    integer(c_int) :: flags

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with draw_view in windows@view module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_type = 3
    integer(c_int), parameter :: ic_editbutton = 4

    ! initialization
    szero%x = 0
    szero%y = 0
    changed = .false.

    ! sort the representation array
    if (s%forcesort) then
       allocate(idx(s%nrep))
       do i = 1, s%nrep
          s%iord(i) = i
          idx(i) = s%rep(i)%iord
       end do
       call mergesort(idx,s%iord,1,s%nrep)
       deallocate(idx)
       s%forcesort = .false.
    end if

    ! representation rows
    do ii = 1, s%nrep
       i = s%iord(ii)
       if (.not.s%rep(i)%isinit) cycle

       ! close button
       doerase = .false.
       call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
       if (igTableSetColumnIndex(ic_closebutton)) then
          call igAlignTextToFramePadding()
          str1 = "##2ic_closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
          if (my_CloseButton(c_loc(str1),ColorDangerButton)) doerase = .true.
       end if

       ! view button
       if (igTableSetColumnIndex(ic_viewbutton)) then
          if (iw_checkbox("##2ic_viewbutton" // string(ic_viewbutton) // "," // string(i),s%rep(i)%shown)) &
             changed = .true.
       end if

       ! name
       discol = .not.s%rep(i)%shown
       if (igTableSetColumnIndex(ic_name)) then
          if (discol) &
             call igPushStyleColor_Vec4(ImGuiCol_Text,g%Style%Colors(ImGuiCol_TextDisabled+1))
          call iw_text(trim(s%rep(i)%name))
          if (discol) call igPopStyleColor(1)

          ! name context menu
          if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_MouseButtonRight)) then
             ! edit
             if (iw_menuitem("Edit")) &
                idum = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idparent=idparent,orraise=-1)
             call iw_tooltip("Edit this object",ttshown)

             ! duplicate
             if (iw_menuitem("Duplicate")) then
                id = s%get_new_representation_id()
                s%rep(id) = s%rep(i)
                s%rep(id)%name = trim(s%rep(i)%name) // " (copy)"
                s%icount(s%rep(i)%flavor) = s%icount(s%rep(i)%flavor) + 1
                s%icount(0) = s%icount(0) + 1
                s%rep(id)%iord = s%icount(0)
                s%forcesort = .true.
                changed = .true.
             end if
             call iw_tooltip("Make a copy of this object",ttshown)

             ! show/hide
             if (iw_menuitem("Show/Hide")) then
                s%rep(i)%shown = .not.s%rep(i)%shown
                changed = .true.
             end if
             call iw_tooltip("Toggle hide/show of this object",ttshown)

             ! rename
             str2 = "Rename" // c_null_char
             if (igBeginMenu(c_loc(str2),.true._c_bool)) then
                if (iw_inputtext("##inputrenamerep",bufsize=1023,texta=s%rep(i)%name,width=30,grabfocus=.true.,&
                   notlive=.true.,flags=ImGuiInputTextFlags_AutoSelectAll)) then
                   call igCloseCurrentPopup()
                end if
                call igEndMenu()
             end if
             call iw_tooltip("Rename this object",ttshown)

             ! delete
             if (iw_menuitem("Delete")) &
                doerase = .true.
             call iw_tooltip("Delete this object",ttshown)

             call igEndPopup()
          end if
       end if

       ! type
       if (igTableSetColumnIndex(ic_type)) then
          if (discol) &
             call igPushStyleColor_Vec4(ImGuiCol_Text,g%Style%Colors(ImGuiCol_TextDisabled+1))
          if (s%rep(i)%type == reptype_atoms) then
             str3 = "atoms" // c_null_char
          elseif (s%rep(i)%type == reptype_unitcell) then
             str3 = "cell" // c_null_char
          elseif (s%rep(i)%type == reptype_axes) then
             str3 = "axes" // c_null_char
          elseif (s%rep(i)%type == reptype_symelem) then
             str3 = "symmetry" // c_null_char
          else
             str3 = "???" // c_null_char
          end if
          flags = ImGuiSelectableFlags_SpanAllColumns
          flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
          flags = ior(flags,ImGuiSelectableFlags_DontClosePopups)
          flags = ior(flags,ImGuiSelectableFlags_AllowDoubleClick)
          ok = igSelectable_Bool(c_loc(str3),.false._c_bool,flags,szero)
          if (ok .and. igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) then
             s%rep(i)%shown = .not.s%rep(i)%shown
             changed = .true.
          end if

          if (discol) call igPopStyleColor(1)
       end if

       ! edit button
       if (igTableSetColumnIndex(ic_editbutton)) then
          if (iw_button("Edit##2ic_editbutton" // string(ic_editbutton) // "," // string(i))) then
             idum = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idparent=idparent,&
                orraise=-1)
          end if
       end if

       ! delete the representation if asked
       if (doerase) then
          call s%rep(i)%end()
          changed = .true.
       end if
    end do

  end function representation_menu

  !> Get the ID for a new representation. If necessary, reallocate the
  !> representations array.
  module function get_new_representation_id(s) result(id)
    use types, only: realloc
    use windows, only: regenerate_window_pointers
    class(scene), intent(inout), target :: s
    integer :: id

    integer :: i
    type(representation), allocatable :: aux(:)

    ! try to find an empty spot
    do i = 1, s%nrep
       if (.not.s%rep(i)%isinit) then
          id = i
          return
       end if
    end do

    ! make new representation at the end
    s%nrep = s%nrep + 1
    if (s%nrep > size(s%rep,1)) then
       allocate(aux(2*s%nrep))
       aux(1:size(s%rep,1)) = s%rep
       call move_alloc(aux,s%rep)
       call realloc(s%iord,2*s%nrep)
       call regenerate_window_pointers()
    end if
    id = s%nrep

  end function get_new_representation_id

  !> Update the projection matrix from the v_pos
  module subroutine update_projection_matrix(s)
    use interfaces_glfw, only: glfwGetTime
    use utils, only: mult, infiniteperspective
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, sc(3), znear

    pic = real(pi,c_float)

    if (s%isortho) then
       ! orthographic: the frustum half-width follows the camera distance, so
       ! moving the camera (zoom) changes the apparent size
       call ortho_projection(s,s%projection)
    else
       ! perspective: the field of view is fixed and the apparent size follows
       ! the camera distance. The near plane is placed just in front of the
       ! scene (must be strictly positive); the far plane is at infinity. The
       ! viewport is square, so the aspect ratio is 1.
       call mult(sc,s%world,s%scenecenter)
       znear = norm2(s%campos - sc) - s%scenerad
       if (znear < 0.01_c_float * s%scenerad) znear = 0.01_c_float * s%scenerad
       call infiniteperspective(s%projection,s%persp_fov * pic / 180._c_float,1._c_float,znear)
    end if

    ! the camera projection has changed: stamp the change time so consumers of
    ! the camera state (locked-camera sync, caches) see direct projection
    ! edits (e.g. the orthographic/perspective toggle), not only camera moves
    s%timelastcamchange = glfwGetTime()

  end subroutine update_projection_matrix

  !> Ratio between the window-anchored gizmo's on-screen size with "scale with
  !> zoom" off vs on.
  module function scene_gizmo_zoom_factor(s) result(f)
    use utils, only: mult
    use param, only: pi
    class(scene), intent(in) :: s
    real(c_float) :: f

    real(c_float) :: sc(3), hw2, hside

    call mult(sc,s%world,s%scenecenter)
    hw2 = tan(0.5_c_float * s%ortho_fov * real(pi,c_float) / 180._c_float) * norm2(s%campos - sc)
    hw2 = max(hw2,1e-4_c_float)
    hside = reset_zoom_hside(s)
    f = hw2 / hside

  end function scene_gizmo_zoom_factor

  !> Half of the visible window side at the reset zoom (tworld units): the
  !> quantity scene_reset uses to place the camera. Constant on-screen-size
  !> items (labels, the window-anchored gizmo) divide by it so that the two
  !> zoom-scaling modes coincide at the reset zoom.
  function reset_zoom_hside(s) result(hside)
    class(scene), intent(in) :: s
    real(c_float) :: hside

    hside = s%camresetdist * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
    hside = hside * s%camratio
    hside = max(hside,3._c_float)

  end function reset_zoom_hside

  subroutine ortho_projection(s,proj,symz)
    use utils, only: ortho, mult
    use param, only: pi
    class(scene), intent(in) :: s
    real(c_float), intent(out) :: proj(4,4)
    logical, intent(in), optional :: symz

    real(c_float) :: sc(3), hw2, zf
    logical :: sym

    sym = .false.
    if (present(symz)) sym = symz

    call mult(sc,s%world,s%scenecenter)
    hw2 = tan(0.5_c_float * s%ortho_fov * real(pi,c_float) / 180._c_float) * norm2(s%campos - sc)
    ! guard against a degenerate (zero-width or zero-depth) frustum: if the
    ! camera reaches the scene center (hw2->0) or the scene has no extent
    ! (scenerad->0), the projection matrix becomes singular and its inverse
    ! (unproject) divides by zero -> SIGFPE.
    hw2 = max(hw2,1e-4_c_float)
    zf = max((s%camresetdist * max_zoom) * s%scenerad,1e-4_c_float)
    if (sym) then
       ! symmetric depth range for the window-anchored gizmo: its geometry is
       ! centered at the eye origin (anchoring drops the eye translation), so a
       ! [0,zf] range would clip the half pointing toward the camera. Centering
       ! the range on z=0 keeps the whole widget inside the frustum.
       call ortho(proj,-hw2,hw2,-hw2,hw2,-zf,zf)
    else
       call ortho(proj,-hw2,hw2,-hw2,hw2,0._c_float,zf)
    end if

  end subroutine ortho_projection

  !> Update the view matrix from the v_pos, v_front, and v_up
  module subroutine update_view_matrix(s)
    use utils, only: lookat
    class(scene), intent(inout), target :: s

    call lookat(s%view,s%campos,s%campos+s%camfront,s%camup)

  end subroutine update_view_matrix

  !> Align the view with a given scene axis. a,b,c = 1,2,3 and x,y,z =
  !> -1,-2,-3.
  module subroutine align_view_axis(s,iaxis)
    use systems, only: sys
    use tools_math, only: cross
    use utils, only: rotate, translate
    class(scene), intent(inout), target :: s
    integer, intent(in) :: iaxis

    real*8 :: xaxis(3), oaxis(3), raxis(3)
    real(c_float) :: raxis_c(3), angle

    ! alignment axis
    if (iaxis == 1) then
       xaxis = sys(s%id)%c%m_x2c(:,1)
    elseif (iaxis == 2) then
       xaxis = sys(s%id)%c%m_x2c(:,2)
    elseif (iaxis == 3) then
       xaxis = sys(s%id)%c%m_x2c(:,3)
    elseif (iaxis == -1) then
       xaxis = (/1d0,0d0,0d0/)
    elseif (iaxis == -2) then
       xaxis = (/0d0,1d0,0d0/)
    elseif (iaxis == -3) then
       xaxis = (/0d0,0d0,1d0/)
    else
       return
    end if
    xaxis = xaxis / norm2(xaxis)

    oaxis = (/0d0,0d0,1d0/)
    raxis = cross(oaxis,xaxis)
    angle = real(asin(norm2(raxis)),c_float)

    ! reset the camera position
    call s%reset()

    ! set the world matrix
    if (angle > 1e-10_c_float) then
       raxis_c = real(raxis / norm2(raxis),c_float)
       call translate(s%world,s%scenecenter)
       call rotate(s%world,-angle,raxis_c)
       call translate(s%world,-s%scenecenter)
    end if

  end subroutine align_view_axis

  !> Add atom idx to the measure selection set. If idx(1) = 0,
  !> clear the measure selection.
  module subroutine select_atom(s,idx)
    class(scene), intent(inout), target :: s
    integer, intent(in) :: idx(5)

    integer :: i, j

    ! if no atom, clear
    if (idx(1) == 0) then
       s%nmsel = 0
       return
    end if

    ! if the atom is already selected, deselect
    do i = 1, s%nmsel
       if (idx(5) == s%msel(5,i)) then
          do j = i+1, s%nmsel
             s%msel(:,j-1) = s%msel(:,j)
          end do
          s%nmsel = s%nmsel - 1
          return
       end if
    end do

    if (s%nmsel < 4) then
       ! if the atom is not known and we have space for it, add it
       s%nmsel = s%nmsel + 1
       s%msel(:,s%nmsel) = idx
    else
       ! if we have 4 atoms selected and clicked a different atom, select only that one
       s%nmsel = 1
       s%msel(:,1) = idx
    end if

  end subroutine select_atom

  !> Add a new representation to the scene with type itype and the given flavor.
  !> If id is present, it returns the new representation id.
  module subroutine add_representation(s,itype,flavor,id)
    class(scene), intent(inout), target :: s
    integer, intent(in) :: itype
    integer, intent(in) :: flavor
    integer, intent(out), optional :: id

    integer :: id_

    id_ = s%get_new_representation_id()
    call s%rep(id_)%init(s%id,id_,itype,flavor,s%icount)
    s%forcesort = .true.
    s%forcebuildlists = .true.
    if (present(id)) id = id_

  end subroutine add_representation

  !> Add a transient representation (type itype, flavor) to the scene and
  !> return its id so the caller can configure it. Transient
  !> representations live in a separate list from the user representations:
  !> they are drawn in the scene but never appear in the object menu and
  !> are not user-controllable. They are cleared automatically each frame
  !> unless re-armed (see scene_clear_transient_representations and the
  !> per-frame logic in the main loop). A throwaway name counter is used so
  !> the scene's user-representation counters are untouched.
  module function scene_add_transient_representation(s,itype,flavor) result(id)
    use representations, only: repflavor_NUM
    use systems, only: sysc, lastchange_render
    class(scene), intent(inout), target :: s
    integer, intent(in) :: itype
    integer, intent(in) :: flavor
    integer :: id

    type(representation), allocatable :: aux(:)
    integer :: dummycount(0:repflavor_NUM)

    ! grow the transient list
    s%nreptrans = s%nreptrans + 1
    if (.not.allocated(s%reptrans)) allocate(s%reptrans(4))
    if (s%nreptrans > size(s%reptrans,1)) then
       allocate(aux(2*s%nreptrans))
       aux(1:size(s%reptrans,1)) = s%reptrans
       call move_alloc(aux,s%reptrans)
    end if
    id = s%nreptrans

    ! initialize the representation (throwaway counter, not s%icount). Transient
    ! reps are single-element and never run update(), so their symelem_style stays
    ! uninitialized and the single-element draw path is used.
    dummycount = 0
    call s%reptrans(id)%init(s%id,id,itype,flavor,dummycount)

    ! trigger a rebuild + re-render and keep the transient set alive this frame
    s%forcebuildlists = .true.
    s%reptrans_set = .true.
    if (s%id > 0) call sysc(s%id)%post_event(lastchange_render)

  end function scene_add_transient_representation

  !> Clear all transient representations from the scene. Triggers a rebuild
  !> of the draw lists (and a re-render) only when there was something to
  !> remove, so it is cheap to call every frame.
  module subroutine scene_clear_transient_representations(s)
    use systems, only: sysc, lastchange_render
    class(scene), intent(inout), target :: s

    integer :: i

    if (s%nreptrans <= 0) return
    do i = 1, s%nreptrans
       call s%reptrans(i)%end()
    end do
    s%nreptrans = 0
    s%reptrans_tag = -1
    s%forcebuildlists = .true.
    if (s%id > 0) call sysc(s%id)%post_event(lastchange_render)

  end subroutine scene_clear_transient_representations

  !> Show a transient set of standard-orientation axes at the Cartesian
  !> (bohr) point xcom, oriented by the rotation matrix rot (columns are
  !> the axis directions), with each axis of length axlen. tag identifies
  !> the content (e.g. the molecule index) so that re-hovering the same
  !> item just keeps the axes alive without rebuilding the draw lists.
  module subroutine scene_show_transient_axes(s,tag,xcom,rot,axlen)
    use representations, only: reptype_axes, repflavor_axes
    class(scene), intent(inout), target :: s
    integer, intent(in) :: tag
    real*8, intent(in) :: xcom(3)
    real*8, intent(in) :: rot(3,3)
    real*8, intent(in) :: axlen

    integer :: id

    ! already shown for this tag (same molecule): keep the axis geometry but
    ! refresh its transform, so the axes track the molecule as it is dragged
    ! (rotated/translated)
    if (s%reptrans_tag == tag .and. s%nreptrans > 0) then
       do id = 1, s%nreptrans
          if (s%reptrans(id)%type == reptype_axes) then
             s%reptrans(id)%axes%origin = xcom
             s%reptrans(id)%axes%rot = rot
             s%reptrans(id)%axes%scale = axlen / max(s%reptrans(id)%axes%length,1d-10)
          end if
       end do
       s%reptrans_set = .true.
       return
    end if

    ! (re)build the transient axes
    call s%clear_transient_representations()
    id = s%add_transient_representation(reptype_axes,repflavor_axes)
    if (id <= 0) return
    if (.not.s%reptrans(id)%isinit) return

    s%reptrans(id)%axes%placement = 0 ! drawn in world space (not a window gizmo)
    s%reptrans(id)%axes%coordtype = 2 ! origin in cartesian (bohr)
    s%reptrans(id)%axes%origin = xcom
    s%reptrans(id)%axes%kind = 0 ! cartesian base directions, reoriented by axes_rot
    s%reptrans(id)%axes%rot = rot ! orient along the molecule's principal axes
    s%reptrans(id)%axes%showlabels = .false.
    s%reptrans(id)%axes%scale = axlen / max(s%reptrans(id)%axes%length,1d-10)
    s%reptrans(id)%axes%scale_auto = .false.
    s%reptrans_tag = tag

  end subroutine scene_show_transient_axes

  !> Show a (transient) black cylinder marking the rotation axis
  !> rotdir (cartesian unit vector, bohr) through the center of mass
  !> xcom, with half-length rotlen.
  module subroutine scene_show_transient_rotaxis(s,tag,xcom,rotdir,rotlen)
    use representations, only: reptype_rotaxis, repflavor_rotaxis
    class(scene), intent(inout), target :: s
    integer, intent(in) :: tag
    real*8, intent(in) :: xcom(3)
    real*8, intent(in) :: rotdir(3)
    real*8, intent(in) :: rotlen

    integer :: id

    ! already shown for this tag: keep the geometry but refresh the transforms,
    ! so both gizmos track the molecule as it is dragged (rotated)
    if (s%reptrans_tag == tag .and. s%nreptrans > 0) then
       do id = 1, s%nreptrans
          if (s%reptrans(id)%type == reptype_rotaxis) then
             s%reptrans(id)%rotaxis%origin = xcom
             s%reptrans(id)%rotaxis%dir = rotdir
             s%reptrans(id)%rotaxis%length = rotlen
          end if
       end do
       s%reptrans_set = .true.
       return
    end if

    ! rebuild: the rotation-axis cylinder (black, through the COM)
    id = s%add_transient_representation(reptype_rotaxis,repflavor_rotaxis)
    if (id > 0) then
       if (s%reptrans(id)%isinit) then
          s%reptrans(id)%rotaxis%origin = xcom
          s%reptrans(id)%rotaxis%dir = rotdir
          s%reptrans(id)%rotaxis%length = rotlen
       end if
    end if
    s%reptrans_tag = tag

  end subroutine scene_show_transient_rotaxis

  !> Show a set of n symmetry elements as transient
  !> representations. Each element k is a plane
  !> (kind(k)=symop_kind_plane, dir = plane normal) or an axis
  !> (kind(k)=symop_kind_axis, dir = axis direction); xorig and dir
  !> are in cartesian (bohr) and the elements are sized to span the
  !> displayed system. tag identifies the requested set for dedup:
  !> when it matches the set already shown, the element geometry is
  !> refreshed in place (no list rebuild), otherwise the transient
  !> list is rebuilt from scratch. Like the other transient
  !> representations, this must be called every frame to keep the
  !> elements alive.
  module subroutine scene_show_symelems(s,tag,n,kind,xorig,dir,order)
    use representations, only: reptype_symelem, repflavor_symelem
    class(scene), intent(inout), target :: s
    integer, intent(in) :: tag
    integer, intent(in) :: n
    integer, intent(in) :: kind(n)
    real*8, intent(in) :: xorig(3,n)
    real*8, intent(in) :: dir(3,n)
    integer, intent(in) :: order(n)

    integer :: k, id

    ! nothing to show: let the main loop auto-clear the transients
    if (n <= 0) return

    ! same set already shown: refresh the element geometry, keep it alive
    if (s%reptrans_tag == tag .and. s%nreptrans == n) then
       do id = 1, n
          if (s%reptrans(id)%type /= reptype_symelem) cycle
          s%reptrans(id)%symelem%kind = kind(id)
          s%reptrans(id)%symelem%origin_transient = xorig(:,id)
          s%reptrans(id)%symelem%dir = dir(:,id)
          s%reptrans(id)%symelem%size = s%scenerad
          s%reptrans(id)%symelem%cen = real(s%scenecenter,8)
          s%reptrans(id)%symelem%order = order(id)
       end do
       s%reptrans_set = .true.
       return
    end if

    ! (re)build the transient symmetry elements
    call s%clear_transient_representations()
    do k = 1, n
       id = s%add_transient_representation(reptype_symelem,repflavor_symelem)
       if (id <= 0) cycle
       if (.not.s%reptrans(id)%isinit) cycle
       s%reptrans(id)%symelem%kind = kind(k)
       s%reptrans(id)%symelem%origin_transient = xorig(:,k)
       s%reptrans(id)%symelem%dir = dir(:,k)
       s%reptrans(id)%symelem%size = s%scenerad
       s%reptrans(id)%symelem%cen = real(s%scenecenter,8)
       s%reptrans(id)%symelem%order = order(k)
    end do
    s%reptrans_tag = tag

  end subroutine scene_show_symelems

  !xx! private procedures: low-level draws

  !> Common entry preparation for scene_render and scene_render_pick: create
  !> the GL buffers on first use, (re)build the draw lists if needed, copy the
  !> camera from the most recently moved member of the locking group, and
  !> reset the camera if required.
  subroutine scene_render_prepare(s)
    use systems, only: sysc, nsys
    class(scene), intent(inout), target :: s

    integer :: i, ifound
    real*8 :: time

    ! create this scene's GL instance buffers on first render
    if (.not.s%gl%isinit) call s%gl%init()

    ! build the draw lists if not done already or if a rebuild is forced
    if (s%isinit == 1 .or. .not.allocated(s%obj%sph) .or. s%forcebuildlists) &
       call s%build_lists()

    ! if the camera is locked, copy the camera parameters from the member
    ! of the locking group who was moved last
    if (s%lockedcam /= 0) then
       ifound = 0
       time = s%timelastcamchange
       do i = 1, nsys
          if (sysc(i)%sc%lockedcam == s%lockedcam .and. sysc(i)%sc%timelastcamchange > time) then
             ifound = i
             time = sysc(i)%sc%timelastcamchange
          end if
       end do
       if (ifound > 0) call s%cam_copy(sysc(ifound)%sc)
    end if

    ! if necessary, reset the camera
    if (s%forceresetcam) call s%reset()

  end subroutine scene_render_prepare

  !> Pack one sphere instance into column col of an instance buffer (layout
  !> must match the sphinstVAO attribute offsets set in glbuffers_init).
  subroutine sphere_pack(col,x,r,rgba,border,bcol,xdelta,ridx,occ,occempty)
    use shapes, only: sph_inst_nf
    real(c_float), intent(out) :: col(sph_inst_nf)
    real(c_float), intent(in) :: x(3), r, rgba(4), border, bcol(3)
    complex(c_float_complex), intent(in) :: xdelta(3)
    real(c_float), intent(in) :: ridx(4)
    real(c_float), intent(in) :: occ
    real(c_float), intent(in) :: occempty(3)

    col(1:3) = x
    col(4) = r
    col(5:8) = rgba
    col(9) = border
    col(10:12) = bcol
    col(13:15) = real(xdelta,c_float)
    col(16:18) = real(aimag(xdelta),c_float)
    col(19:22) = ridx
    col(23) = occ
    col(24:26) = occempty

  end subroutine sphere_pack

  !> Pack one cylinder instance into column col of an instance buffer (layout
  !> must match the cylinstVAO attribute offsets set in glbuffers_init).
  !> x1/x2 are the equilibrium endpoints; xd1/xd2 the vibration deltas.
  subroutine cyl_pack(col,x1,x2,r,rgba,border,bcol,delta,outward,xd1,xd2)
    use shapes, only: cyl_inst_nf
    real(c_float), intent(out) :: col(cyl_inst_nf)
    real(c_float), intent(in) :: x1(3), x2(3), r, rgba(4), border, bcol(3), delta, outward(3)
    complex(c_float_complex), intent(in) :: xd1(3), xd2(3)

    col(1:3) = x1
    col(4:6) = x2
    col(7) = r
    col(8:11) = rgba
    col(12) = border
    col(13:15) = bcol
    col(16) = delta
    col(17:19) = outward
    col(20:22) = real(xd1,c_float)
    col(23:25) = aimag(xd1)
    col(26:28) = real(xd2,c_float)
    col(29:31) = aimag(xd2)

  end subroutine cyl_pack

  !> Pack one mesh instance (model matrix columns + color + per-vertex
  !> vibration deltas) into column col. The deltas d1/d2/d3 correspond to the
  !> reference-triangle corners (polyhedra faces); they are zero when absent
  !> (planes, cones).
  subroutine mesh_pack(col,model,rgba,d1,d2,d3)
    use shapes, only: mesh_inst_nf
    real(c_float), intent(out) :: col(mesh_inst_nf)
    real(c_float), intent(in) :: model(4,4), rgba(4)
    complex(c_float_complex), intent(in), optional :: d1(3), d2(3), d3(3)

    col(1:16) = reshape(model,(/16/)) ! column-major: col1..col4
    col(17:20) = rgba
    col(21:38) = 0._c_float
    if (present(d1)) then
       col(21:23) = real(d1,c_float)
       col(24:26) = aimag(d1)
    end if
    if (present(d2)) then
       col(27:29) = real(d2,c_float)
       col(30:32) = aimag(d2)
    end if
    if (present(d3)) then
       col(33:35) = real(d3,c_float)
       col(36:38) = aimag(d3)
    end if

  end subroutine mesh_pack

  !> Build the model matrix that maps the unit cone mesh onto the cone from
  !> base x1 to apex x2 with base radius rad (same convention as draw_cone).
  subroutine cone_model(x1,x2,rad,model)
    use tools_math, only: cross_cfloat
    real(c_float), intent(in) :: x1(3), x2(3), rad
    real(c_float), intent(out) :: model(4,4)

    real(c_float) :: xmid(3), xdif(3), up(3), crs(3), blen, a, ca, sa, axis(3), temp(3)

    model = eye4
    xmid = 0.5_c_float * (x1 + x2)
    xdif = x2 - x1
    blen = norm2(xdif)
    if (blen < 1e-4_c_float) return
    xdif = xdif / blen
    up = (/0._c_float,0._c_float,1._c_float/)
    crs = cross_cfloat(up,xdif)
    model(1:3,4) = xmid
    if (dot_product(crs,crs) > 1e-14_c_float) then
       a = acos(dot_product(xdif,up))
       ca = cos(a)
       sa = sin(a)
       axis = crs / norm2(crs)
       temp = (1._c_float - ca) * axis
       model(1,1) = ca + temp(1) * axis(1)
       model(2,1) = temp(1) * axis(2) + sa * axis(3)
       model(3,1) = temp(1) * axis(3) - sa * axis(2)
       model(1,2) = temp(2) * axis(1) - sa * axis(3)
       model(2,2) = ca + temp(2) * axis(2)
       model(3,2) = temp(2) * axis(3) + sa * axis(1)
       model(1,3) = temp(3) * axis(1) + sa * axis(2)
       model(2,3) = temp(3) * axis(2) - sa * axis(1)
       model(3,3) = ca + temp(3) * axis(3)
    end if
    model(:,1) = model(:,1) * rad
    model(:,2) = model(:,2) * rad
    model(:,3) = model(:,3) * blen

  end subroutine cone_model

  !> Calculate the vertices for the given text and adds them to
  !> nvert/vert, on-scene version. x0 = world position of the label
  !> (equilibrium; the vibration displacement is applied in the shader from
  !> the xdelta packed with each vertex). r = radius of the associated atom.
  subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,shift,centered,xdelta)
    use interfaces_cimgui
    use shapes, only: text_vert_nf
    use gui_main, only: g, fontbakesize_large
    use types, only: realloc
    use param, only: newline, bohrtoa
    character(len=*), intent(in) :: text
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: r
    real(c_float), intent(in) :: siz
    integer(c_int), intent(inout) :: nvert
    real(c_float), allocatable, intent(inout) :: vert(:,:)
    real(c_float), intent(in), optional :: shift(3)
    logical, intent(in), optional :: centered
    complex(c_float_complex), intent(in), optional :: xdelta(3)

    integer :: i, j, nline, nvert0
    type(c_ptr) :: cptr
    type(ImFontGlyph), pointer :: glyph
    real(c_float) :: xpos, ypos, lheight, shift_(3)
    complex(c_float_complex) :: xdelta_(3)
    logical :: centered_
    real(c_float), allocatable :: xlen(:)
    integer, allocatable :: jlen(:)

    real(c_float), parameter :: rshift = 0.01_c_float

    ! initialize
    centered_ = .false.
    if (present(centered)) centered_ = centered
    shift_ = 0._c_float
    if (present(shift)) shift_ = shift / real(bohrtoa,c_float)
    xdelta_ = (0._c_float,0._c_float)
    if (present(xdelta)) xdelta_ = xdelta
    if (.not.allocated(vert)) then
       allocate(vert(text_vert_nf,100))
       nvert = 0
    end if
    nvert0 = nvert

    ! initial variables
    xpos = 0._c_float
    ypos = 0._c_float
    lheight = fontbakesize_large
    if (centered_) then
       nline = 1
       allocate(xlen(10),jlen(10))
       jlen(1) = nvert+1
       xlen(1) = 0._c_float
    end if

    ! loop over characters
    i = 0
    do while (i < len_trim(text))
       i = i + 1
       ! newline, skip line and advance one (linux)
       if (text(i:i) == newline) then
          xpos = 0._c_float
          ypos = ypos + lheight
          i = i + 1
          if (centered_) then
             nline = nline + 1
             if (nline+1 > size(xlen,1)) then
                call realloc(xlen,2*nline)
                call realloc(jlen,2*nline)
             end if
             xlen(nline) = 0._c_float
             jlen(nline) = nvert+1
          end if
          continue
       end if

       ! get the glyph
       cptr = ImFont_FindGlyph(g%Font,int(ichar(text(i:i)),c_int16_t))
       call c_f_pointer(cptr,glyph)

       ! add to the vertices
       if (nvert+6 > size(vert,2)) call realloc(vert,int(text_vert_nf),2*(nvert+6))
       do j = nvert+1, nvert+6
          vert(1:3,j) = x0
          vert(4:5,j) = shift_(1:2)
          vert(6,j) = r + rshift + shift_(3)
          vert(11:13,j) = real(xdelta_,c_float)
          vert(14:16,j) = aimag(xdelta_)
       end do

       vert(7:8,nvert+1) = (/xpos + glyph%X0, ypos + glyph%Y1/)
       vert(7:8,nvert+2) = (/xpos + glyph%X0, ypos + glyph%Y0/)
       vert(7:8,nvert+3) = (/xpos + glyph%X1, ypos + glyph%Y0/)
       vert(7:8,nvert+4) = (/xpos + glyph%X0, ypos + glyph%Y1/)
       vert(7:8,nvert+5) = (/xpos + glyph%X1, ypos + glyph%Y0/)
       vert(7:8,nvert+6) = (/xpos + glyph%X1, ypos + glyph%Y1/)

       vert(9:10,nvert+1) = (/glyph%U0, glyph%V1/)
       vert(9:10,nvert+2) = (/glyph%U0, glyph%V0/)
       vert(9:10,nvert+3) = (/glyph%U1, glyph%V0/)
       vert(9:10,nvert+4) = (/glyph%U0, glyph%V1/)
       vert(9:10,nvert+5) = (/glyph%U1, glyph%V0/)
       vert(9:10,nvert+6) = (/glyph%U1, glyph%V1/)
       nvert = nvert + 6

       ! advance xpos
       xpos = xpos + glyph%AdvanceX

       ! update
       if (centered_) then
          xlen(nline) = max(xlen(nline),xpos)
       end if
    end do
    if (centered_) then
       jlen(nline+1) = nvert+1

       do i = 1, nline
          vert(7,jlen(i):jlen(i+1)-1) = vert(7,jlen(i):jlen(i+1)-1) - 0.5_c_float * xlen(i)
       end do
       vert(8,jlen(1):nvert) = vert(8,jlen(1):nvert) - 0.5_c_float * nline * lheight
    end if
    ! scale to on-screen size. The render texture is presented right-side up
    ! (see draw_view) but ImGui glyph quads use a y-down layout, so the vertical
    ! offset is negated (folded into the scaling) to keep on-scene text upright.
    ! This reverses the quad winding, so culling must be disabled where these
    ! vertices are drawn.
    vert(7,nvert0+1:nvert) =  vert(7,nvert0+1:nvert) * siz
    vert(8,nvert0+1:nvert) = -vert(8,nvert0+1:nvert) * siz

  end subroutine calc_text_onscene_vertices

end submodule proc
