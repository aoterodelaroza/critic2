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
! the -y direction.
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

  ! space for uniforms
  integer, parameter :: iu_world = 1
  integer, parameter :: iu_view = 2
  integer, parameter :: iu_projection = 3
  integer, parameter :: iu_model = 4
  integer, parameter :: iu_object_type = 5
  integer, parameter :: iu_border = 6
  integer, parameter :: iu_bordercolor = 7
  integer, parameter :: iu_vcolor = 8
  integer, parameter :: iu_idx = 9
  integer, parameter :: iu_delta_cyl = 10
  integer, parameter :: iu_ndash_cyl = 11
  integer, parameter :: iu_NUM = 12
  integer(c_int) :: iunif(iu_NUM)

  !xx! private procedures: low-level draws
  ! subroutine draw_sphere(x0,rad,ires,rgb,index)
  ! subroutine draw_cylinder(x1,x2,rad,rgb,ires)
  ! subroutine calc_text_direct_vertices(text,x0,y0,siz,nvert,vert,centered)
  ! subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,centered)

contains

  !xx! scene

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use representations, only: reptype_atoms, reptype_unitcell,&
       repflavor_atoms_ballandstick, repflavor_atoms_sticks,&
       repflavor_unitcell_basic, repflavor_NUM
    use systems, only: sys, sysc, sys_ready, ok_system
    use global, only: crsmall
    use gui_main, only: lockbehavior
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

    ! shader default settings
    s%style = style_simple
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

    ! atoms
    if (sys(isys)%c%ncel <= crsmall) then
       call s%add_representation(reptype_atoms,repflavor_atoms_ballandstick)
    else
       call s%add_representation(reptype_atoms,repflavor_atoms_sticks)
    end if

    ! unit cell
    if (.not.sys(isys)%c%ismolecule) &
       call s%add_representation(reptype_unitcell,repflavor_unitcell_basic)

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
    if (allocated(s%rep)) deallocate(s%rep)
    if (allocated(s%icount)) deallocate(s%icount)
    if (allocated(s%iord)) deallocate(s%iord)
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
    hside = s%camresetdist * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
    hside = hside * s%camratio
    hside = max(hside,3._c_float)
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
       call s%rep(irep)%atom_style%reset_colors(s%id)
    end do
    s%forcebuildlists = .true.

  end subroutine scene_reset_atom_colors

  !> Build the draw lists for the current scene.
  module subroutine scene_build_lists(s)
    use representations, only: reptype_atoms
    use interfaces_glfw, only: glfwGetTime
    use utils, only: translate
    use systems, only: sys_ready, ok_system
    use interfaces_glfw, only: glfwGetTime
    class(scene), intent(inout), target :: s

    integer :: i
    real(c_float) :: xmin(3), xmax(3), maxrad, xc(3)

    ! only build lists if system is initialized
    if (.not.ok_system(s%id,sys_ready)) return

    ! initialize
    s%obj%nsph = 0
    if (allocated(s%obj%sph)) deallocate(s%obj%sph)
    allocate(s%obj%sph(100))
    s%obj%ncyl = 0
    if (allocated(s%obj%cyl)) deallocate(s%obj%cyl)
    allocate(s%obj%cyl(100))
    s%obj%ncylflat = 0
    if (allocated(s%obj%cylflat)) deallocate(s%obj%cylflat)
    allocate(s%obj%cylflat(10))
    s%obj%nstring = 0
    if (allocated(s%obj%string)) deallocate(s%obj%string)
    allocate(s%obj%string(10))

    ! add the items by representation
    do i = 1, s%nrep
       ! update to reflect changes in the number of atoms or molecules
       call s%rep(i)%update()

       ! add draw elements
       call s%rep(i)%add_draw_elements(s%nc,s%obj,s%animation>0,s%iqpt_selected,s%ifreq_selected)
    end do

    ! reset the measure selection
    s%nmsel = 0
    s%msel = 0

    ! recalculate scene center and radius
    maxrad = 0._c_float
    do i = 1, s%nrep
       if (s%rep(i)%shown .and. s%rep(i)%type == reptype_atoms) &
          maxrad = max(maxrad,maxval(s%rep(i)%atom_style%rad(1:s%rep(i)%atom_style%ntype)))
    end do
    if (s%obj%nsph + s%obj%ncyl + s%obj%ncylflat + s%obj%nstring > 0) then
       do i = 1, 3
          xmin(i) = huge(1._c_float)
          xmax(i) = -huge(1._c_float)

          xmin(i) = minval(s%obj%sph(1:s%obj%nsph)%x(i)) - maxrad
          xmin(i) = min(xmin(i),minval(s%obj%cyl(1:s%obj%ncyl)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cyl(1:s%obj%ncyl)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cylflat(1:s%obj%ncylflat)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%obj%cylflat(1:s%obj%ncylflat)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%obj%string(1:s%obj%nstring)%x(i)))

          xmax(i) = maxval(s%obj%sph(1:s%obj%nsph)%x(i)) + maxrad
          xmax(i) = max(xmax(i),maxval(s%obj%cyl(1:s%obj%ncyl)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cyl(1:s%obj%ncyl)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cylflat(1:s%obj%ncylflat)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%obj%cylflat(1:s%obj%ncylflat)%x2(i)))
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
    if (s%iscaminit) &
       call translate(s%world,-xc)

    ! rebuilding lists is done
    s%forcebuildlists = .false.
    s%isinit = 2
    s%timelastbuild = glfwGetTime()

  end subroutine scene_build_lists

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_glfw, only: glfwGetTime
    use interfaces_cimgui
    use interfaces_opengl3
    use shapes, only: sphVAO, cylVAO, textVAOos, textVBOos
    use gui_main, only: fonts, fontbakesize_large, font_large
    use systems, only: sys, sysc, nsys
    use utils, only: ortho, project
    use tools_math, only: eigsym, matinv_cfloat
    use tools_io, only: string
    use shaders, only: shader_phong, shader_simple, shader_text_onscene,&
       useshader, setuniform_int, setuniform_float, setuniform_vec3,&
       setuniform_vec4, setuniform_mat3, setuniform_mat4, get_uniform_location
    use param, only: img, pi
    class(scene), intent(inout), target :: s

    real(c_float) :: xsel(3,4), radsel(4)
    complex(c_float_complex) :: displ
    real*8 :: deltat, fac, time
    logical :: doit
    integer :: i, ifound

    real(c_float), parameter :: msel_thickness = 0.1_c_float
    real(c_float), parameter :: sel_thickness = 0.2_c_float
    real(c_float), parameter :: sel_label_size = 1.2_c_float
    real*8, parameter :: freq_ref = 300d0
    real*8, parameter :: freq_min = 50d0

    ! check that the scene and system are initialized
    if (s%isinit == 0) return

    ! build draw lists if not done already
    if (s%isinit == 1 .or. .not.allocated(s%obj%sph)) call s%build_lists()

    ! if necessary, rebuild draw lists
    if (s%forcebuildlists) call s%build_lists()

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

    ! get the time
    time = glfwGetTime()

    ! render text with the large font
    call igPushFont(font_large)

    ! calculate the time factor
    displ = 0._c_float
    if (s%ifreq_selected > 0 .and. s%iqpt_selected > 0) then
       fac = s%anim_amplitude * sqrt(freq_ref / max(abs(sys(s%id)%c%vib%freq(s%ifreq_selected,s%iqpt_selected)),freq_min))
       if (s%animation == 1) then ! manual
          displ = cmplx(fac * exp(0.5d0 * s%anim_phase * pi * img),&
             kind=c_float_complex)
       else ! automatic
          deltat = time - s%timerefanimation
          displ = cmplx(fac * exp(deltat * s%anim_speed * img),kind=c_float_complex)
       end if
    end if

    if (s%style == style_phong) then
       !! phong !!
       ! set up the shader and the uniforms
       call useshader(shader_phong)
       call setuniform_int(1_c_int,"uselighting")
       call setuniform_vec3(s%lightpos,"lightPos")
       call setuniform_vec3(s%lightcolor,"lightColor")
       call setuniform_float(s%ambient,"ambient")
       call setuniform_float(s%diffuse,"diffuse")
       call setuniform_float(s%specular,"specular")
       call setuniform_int(s%shininess,"shininess")
       call setuniform_mat4(s%world,"world")
       call setuniform_mat4(s%view,"view")
       call setuniform_mat4(s%projection,"projection")

       ! get all the uniforms
       iunif(iu_model) = get_uniform_location("model")
       iunif(iu_vcolor) = get_uniform_location("vColor")

       ! draw the atoms
       if (s%obj%nsph > 0) then
          call glBindVertexArray(sphVAO(s%atom_res))
          call draw_all_spheres()
       end if

       ! draw the bonds
       if (s%obj%ncyl > 0) then
          call glBindVertexArray(cylVAO(s%bond_res))
          call draw_all_cylinders()
       end if

       ! draw the flat cylinders (unit cell)
       if (s%obj%ncylflat > 0) then
          call setuniform_int(0_c_int,"uselighting")
          call glBindVertexArray(cylVAO(s%uc_res))
          call draw_all_flat_cylinders()
       end if

       ! draw the measure selection atoms
       if (s%nmsel > 0) then
          call setuniform_int(0_c_int,"uselighting")
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          call draw_all_mselections()
          call glDisable(GL_BLEND)
       end if

       ! render labels with on-scene text
       call useshader(shader_text_onscene)
       call setuniform_mat4(s%world,"world")
       call setuniform_mat4(s%view,"view")
       call setuniform_mat4(s%projection,"projection")

       call glDisable(GL_MULTISAMPLE)
       call glEnable(GL_BLEND)
       call glBlendEquation(GL_FUNC_ADD)
       call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)

       call glActiveTexture(GL_TEXTURE0)
       call glBindVertexArray(textVAOos)
       call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))
       call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)

       call draw_all_text()

       ! render selected atom labels with on-scene text
       if (s%nmsel > 0) &
          call draw_selection_text()

       ! draw the highlighted/selected atoms
       doit = .false.
       if (allocated(sysc(s%id)%highlight_rgba)) &
          doit = any(sysc(s%id)%highlight_rgba >= 0._c_float)
       if (.not.doit.and.allocated(sysc(s%id)%highlight_rgba_transient)) &
          doit = any(sysc(s%id)%highlight_rgba_transient >= 0._c_float)
       if (doit) then
          call useshader(shader_phong)
          call setuniform_int(0_c_int,"uselighting")
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          call draw_highlights()
          call glDisable(GL_BLEND)
       end if

       call glEnable(GL_MULTISAMPLE)
       call glDisable(GL_BLEND)
    else
       ! !! simple !!

       ! set up the shader
       call useshader(shader_simple)

       ! get all the uniforms
       iunif(iu_world) = get_uniform_location("world")
       iunif(iu_view) = get_uniform_location("view")
       iunif(iu_projection) = get_uniform_location("projection")
       iunif(iu_model) = get_uniform_location("model")
       iunif(iu_object_type) = get_uniform_location("object_type")
       iunif(iu_border) = get_uniform_location("rborder")
       iunif(iu_bordercolor) = get_uniform_location("bordercolor")
       iunif(iu_vcolor) = get_uniform_location("vColor")
       iunif(iu_idx) = get_uniform_location("idx")
       iunif(iu_delta_cyl) = get_uniform_location("delta_cyl")
       iunif(iu_ndash_cyl) = get_uniform_location("ndash_cyl")

       ! set the common uniforms
       call setuniform_mat4(s%world,idxi=iunif(iu_world))
       call setuniform_mat4(s%view,idxi=iunif(iu_view))
       call setuniform_mat4(s%projection,idxi=iunif(iu_projection))

       ! draw the spheres for the atoms
       call setuniform_int(0_c_int,idxi=iunif(iu_object_type))
       if (s%obj%nsph > 0) then
          call glBindVertexArray(sphVAO(s%atom_res))
          call draw_all_spheres()
       end if

       ! draw the cylinders for the bonds (inherit border from atoms)
       call setuniform_int(1_c_int,idxi=iunif(iu_object_type))
       if (s%obj%ncyl > 0) then
          call glBindVertexArray(cylVAO(s%bond_res))
          call draw_all_cylinders()
       end if
       call setuniform_int(0_c_int,idxi=iunif(iu_ndash_cyl))
       call setuniform_float(0._c_float,idxi=iunif(iu_delta_cyl))

       ! draw the flat cylinders for the unit cell
       call setuniform_float(0._c_float,idxi=iunif(iu_border))
       call setuniform_int(2_c_int,idxi=iunif(iu_object_type))
       if (s%obj%ncylflat > 0) then
          call glBindVertexArray(cylVAO(s%uc_res))
          call draw_all_flat_cylinders()
       end if
       call setuniform_int(0_c_int,idxi=iunif(iu_ndash_cyl))
       call setuniform_float(0._c_float,idxi=iunif(iu_delta_cyl))

       ! draw the measure selection atoms
       if (s%nmsel > 0) then
          call setuniform_int(0_c_int,idxi=iunif(iu_object_type))
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          call draw_all_mselections()
          call glDisable(GL_BLEND)
       end if

       ! render labels with on-scene text
       call useshader(shader_text_onscene)
       call setuniform_mat4(s%world,"world")
       call setuniform_mat4(s%view,"view")
       call setuniform_mat4(s%projection,"projection")

       call glDisable(GL_MULTISAMPLE)
       call glEnable(GL_BLEND)
       call glBlendEquation(GL_FUNC_ADD)
       call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)

       call glActiveTexture(GL_TEXTURE0)
       call glBindVertexArray(textVAOos)
       call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))
       call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)

       call draw_all_text()

       ! render selected atom labels with on-scene text
       if (s%nmsel > 0) &
          call draw_selection_text()

       ! highlight the highlighted/selected atoms
       doit = .false.
       if (allocated(sysc(s%id)%highlight_rgba)) &
          doit = any(sysc(s%id)%highlight_rgba >= 0._c_float)
       if (.not.doit.and.allocated(sysc(s%id)%highlight_rgba_transient)) &
          doit = any(sysc(s%id)%highlight_rgba_transient >= 0._c_float)
       if (doit) then
          call useshader(shader_simple)
          call setuniform_int(0_c_int,idxi=iunif(iu_object_type))
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          call draw_highlights()
          call glDisable(GL_BLEND)
       end if

       call glEnable(GL_MULTISAMPLE)
       call glDisable(GL_BLEND)
    end if

    ! pop the large font
    call igPopFont()

    ! clean up
    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! save the rendering time
    s%timelastrender = time

  contains
    subroutine draw_all_spheres()
      integer :: i
      real(c_float) :: x(3)

      do i = 1, s%obj%nsph
         x = s%obj%sph(i)%x
         if (s%animation > 0) then
            x = x + real(displ * s%obj%sph(i)%xdelta,c_float)
         end if
         call draw_sphere(x,s%obj%sph(i)%r,s%atom_res,rgb=s%obj%sph(i)%rgb,&
            border=s%obj%sph(i)%border,rgbborder=s%obj%sph(i)%rgbborder)
      end do

    end subroutine draw_all_spheres

    subroutine draw_all_cylinders()
      integer :: i
      real(c_float) :: x1(3), x2(3)

      do i = 1, s%obj%ncyl
         x1 = s%obj%cyl(i)%x1
         x2 = s%obj%cyl(i)%x2
         if (s%animation > 0) then
            x1 = x1 + real(displ * s%obj%cyl(i)%x1delta,c_float)
            x2 = x2 + real(displ * s%obj%cyl(i)%x2delta,c_float)
         end if
         call draw_cylinder(x1,x2,s%obj%cyl(i)%r,s%obj%cyl(i)%rgb,s%bond_res,&
            s%obj%cyl(i)%order,s%obj%cyl(i)%border,s%obj%cyl(i)%rgbborder)
      end do

    end subroutine draw_all_cylinders

    subroutine draw_all_flat_cylinders()
      integer :: i

      do i = 1, s%obj%ncylflat
         call draw_cylinder(s%obj%cylflat(i)%x1,s%obj%cylflat(i)%x2,&
            s%obj%cylflat(i)%r,s%obj%cylflat(i)%rgb,s%uc_res,-1,0._c_float)
      end do

    end subroutine draw_all_flat_cylinders

    !> Draw the measure selections
    subroutine draw_all_mselections()
      use gui_main, only: ColorMeasureSelect
      integer :: i, j
      real(c_float) :: x(3)

      do j = 1, s%nmsel
         i = s%msel(5,j)
         x = s%obj%sph(i)%x
         if (s%animation > 0) x = x + real(displ * s%obj%sph(i)%xdelta,c_float)
         call draw_sphere(x,s%obj%sph(i)%r + msel_thickness,s%atom_res,&
            rgba=ColorMeasureSelect(:,j))
         radsel(j) = s%obj%sph(i)%r + msel_thickness
         xsel(:,j) = x
      end do

    end subroutine draw_all_mselections

    !> Draw the highlights on the scene
    subroutine draw_highlights()
      use systems, only: sysc
      integer :: i, id
      real(c_float) :: x(3), rgba(4)

      ! initial checks
      if (s%isinit < 2) return
      if (.not.allocated(sysc(s%id)%highlight_rgba).and.&
         .not.allocated(sysc(s%id)%highlight_rgba_transient)) return

      ! highlight the spheres
      do i = 1, s%obj%nsph
         id = s%obj%sph(i)%idx(1)
         rgba = -1._c_float
         if (allocated(sysc(s%id)%highlight_rgba_transient)) &
            rgba = sysc(s%id)%highlight_rgba_transient(:,id)
         if (any(rgba < 0).and.allocated(sysc(s%id)%highlight_rgba)) &
            rgba = sysc(s%id)%highlight_rgba(:,id)
         if (all(rgba >= 0)) then
            x = s%obj%sph(i)%x
            if (s%animation > 0) x = x + real(displ * s%obj%sph(i)%xdelta,c_float)
            call draw_sphere(x,s%obj%sph(i)%r + sel_thickness,s%atom_res,rgba=rgba,rgbborder=rgba(1:3))
         end if
      end do

    end subroutine draw_highlights

    subroutine draw_all_text()
      integer :: i
      real(c_float) :: hside, siz, x(3)
      integer(c_int) :: nvert
      real(c_float), allocatable, target :: vert(:,:)

      integer :: iu

      iu = get_uniform_location("textColor")

      do i = 1, s%obj%nstring
         call setuniform_vec3(s%obj%string(i)%rgb,idxi=iu)
         nvert = 0
         if (s%obj%string(i)%scale > 0._c_float) then
            hside = s%camresetdist * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
            hside = hside * s%camratio
            hside = max(hside,3._c_float)
            siz = 2 * s%obj%string(i)%scale / fontbakesize_large / hside
         else
            siz = 2 * abs(s%obj%string(i)%scale) * s%projection(1,1) / fontbakesize_large
         end if
         x = s%obj%string(i)%x
         if (s%animation > 0) x = x + real(displ * s%obj%string(i)%xdelta,c_float)

         call calc_text_onscene_vertices(s%obj%string(i)%str,x,s%obj%string(i)%r,&
            siz,nvert,vert,shift=s%obj%string(i)%offset,centered=.true.)
         call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*10*c_sizeof(c_float), c_loc(vert))
         call glDrawArrays(GL_TRIANGLES, 0, nvert)
      end do

    end subroutine draw_all_text

    subroutine draw_selection_text()
      integer :: j
      real(c_float) :: siz
      integer(c_int) :: nvert
      real(c_float), allocatable, target :: vert(:,:)

      integer :: iu

      iu = get_uniform_location("textColor")

      do j = 1, s%nmsel
         call setuniform_vec3((/1._c_float,1._c_float,1._c_float/),idxi=iu)
         siz = sel_label_size * s%projection(1,1) / fontbakesize_large
         nvert = 0
         call calc_text_onscene_vertices(string(j),xsel(:,j),radsel(j),siz,nvert,vert,centered=.true.)
         call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*10*c_sizeof(c_float), c_loc(vert))
         call glDrawArrays(GL_TRIANGLES, 0, nvert)
      end do

    end subroutine draw_selection_text

  end subroutine scene_render

  !> Draw the scene (for object picking)
  module subroutine scene_render_pick(s)
    use interfaces_cimgui
    use interfaces_opengl3
    use systems, only: sysc, nsys
    use shapes, only: sphVAO
    use utils, only: ortho, project
    use tools_math, only: eigsym, matinv_cfloat
    use shaders, only: shader_pickindex, useshader, setuniform_int,&
       setuniform_float, setuniform_vec3, setuniform_vec4, setuniform_mat3,&
       setuniform_mat4, get_uniform_location
    class(scene), intent(inout), target :: s

    integer :: i, ifound
    real*8 :: time

    ! check that the scene and system are initialized
    if (s%isinit < 2) return

    ! build draw lists if not done already
    if (.not.allocated(s%obj%sph)) call s%build_lists()

    ! if necessary, rebuild draw lists
    if (s%forcebuildlists) call s%build_lists()

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

    ! set up the shader and the uniforms
    call useshader(shader_pickindex)
    call setuniform_mat4(s%world,"world")
    call setuniform_mat4(s%view,"view")
    call setuniform_mat4(s%projection,"projection")

    ! set the uniforms
    iunif(iu_model) = get_uniform_location("model")
    iunif(iu_idx) = get_uniform_location("idx")

    ! draw the atoms
    if (s%obj%nsph > 0) then
       call glBindVertexArray(sphVAO(s%atom_res))
       do i = 1, s%obj%nsph
          ! draw the sphere
          call draw_sphere(s%obj%sph(i)%x,s%obj%sph(i)%r,s%atom_res,idx=(/i,0,0,0/))
       end do
       call glBindVertexArray(0)
    end if

  end subroutine scene_render_pick

  !> Set the style defaults for the current scene. If style is
  !> not given, use s%style
  module subroutine scene_set_style_defaults(s,style)
    class(scene), intent(inout), target :: s
    integer(c_int), intent(in), optional :: style

    integer :: style_
    integer :: i

    style_ = s%style
    if (present(style)) style_ = style

    if (style_ == style_phong) then
       !! phong !!
       s%bgcolor = (/0.3_c_float,0.3_c_float,0.3_c_float/)
       s%lightpos = (/20._c_float,20._c_float,0._c_float/)
       s%lightcolor = (/1._c_float,1._c_float,1._c_float/)
       s%ambient = 0.2_c_float
       s%diffuse = 0.4_c_float
       s%specular = 0.6_c_float
       s%shininess = 8_c_int
       do i = 1, s%nrep
          s%rep(i)%uc_rgb = 1._c_float
       end do
    else
       !! simple !!
       s%bgcolor = (/1._c_float,1._c_float,1._c_float/)
       do i = 1, s%nrep
          s%rep(i)%uc_rgb = 0._c_float
       end do
       s%bordercolor = 0._c_float
    end if

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
    use representations, only: reptype_atoms, reptype_unitcell
    use utils, only: iw_text, iw_tooltip, iw_button, iw_checkbox, iw_menuitem
    use windows, only: stack_create_window, wintype_editrep
    use gui_main, only: ColorDangerButton, g
    use tools_io, only: string
    use tools, only: mergesort
    class(scene), intent(inout), target :: s
    integer(c_int), intent(in) :: idparent
    logical :: changed

    integer :: i, ii, id, ll, idum
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    logical :: discol, doerase, ok
    type(ImVec2) :: szero
    character(kind=c_char,len=1024), target :: txtinp
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
                str3 = "##inputrenamerep" // c_null_char
                txtinp = trim(adjustl(s%rep(i)%name)) // c_null_char
                call igSetKeyboardFocusHere(0_c_int)
                if (igInputText(c_loc(str3),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_EnterReturnsTrue,&
                   c_null_funptr,c_null_ptr)) then
                   ll = index(txtinp,c_null_char)
                   s%rep(i)%name = txtinp(1:ll-1)
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
    use utils, only: ortho, mult
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, hw2, sc(3), znear, zfar

    pic = real(pi,c_float)

    ! scene center: world to tworld
    call mult(sc,s%world,s%scenecenter)

    ! near and far planes
    znear = 0._c_float
    zfar = (s%camresetdist * max_zoom) * s%scenerad

    ! update the projection matrix
    hw2 = tan(0.5_c_float * s%ortho_fov * pic / 180._c_float) * norm2(s%campos - sc)
    call ortho(s%projection,-hw2,hw2,-hw2,hw2,znear,zfar)

  end subroutine update_projection_matrix

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
  module subroutine add_representation(s,itype,flavor)
    class(scene), intent(inout), target :: s
    integer, intent(in) :: itype
    integer, intent(in) :: flavor

    integer :: id

    id = s%get_new_representation_id()
    call s%rep(id)%init(s%id,id,itype,s%style,flavor,s%icount)
    s%forcesort = .true.
    s%forcebuildlists = .true.

  end subroutine add_representation

  !xx! private procedures: low-level draws

  !> Draw a sphere with center x0, radius rad and color rgb. Requires
  !> having the sphere VAO bound.
  subroutine draw_sphere(x0,rad,ires,rgb,rgba,idx,border,rgbborder)
    use representations, only: atomborder_def
    use interfaces_opengl3
    use shaders, only: setuniform_vec3, setuniform_vec4, setuniform_mat4, setuniform_float
    use shapes, only: sphnel
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: rad
    integer(c_int), intent(in) :: ires
    real(c_float), intent(in), optional :: rgb(3)
    real(c_float), intent(in), optional :: rgba(4)
    integer(c_int), intent(in), optional :: idx(4)
    real(c_float), intent(in), optional :: border
    real(c_float), intent(in), optional :: rgbborder(3)

    real(c_float) :: model(4,4)
    real(c_float) :: ridx(4), rgb_(4)
    real(c_float) :: border_, rgbborder_(3)

    ! the model matrix: scale and translate
    model = eye4
    model(1:3,4) = x0
    model(1,1) = rad
    model(2,2) = rad
    model(3,3) = rad

    ! border
    border_ = atomborder_def
    if (present(border)) border_ = border
    call setuniform_float(border_,idxi=iunif(iu_border))
    rgbborder_ = 0._c_float
    if (present(rgbborder)) rgbborder_ = rgbborder
    call setuniform_vec3(rgbborder_,idxi=iunif(iu_bordercolor))

    ! set the uniforms
    if (present(rgb)) then
       rgb_(1:3) = rgb
       rgb_(4) = 1._c_float
       call setuniform_vec4(rgb_,idxi=iunif(iu_vcolor))
    elseif (present(rgba)) then
       call setuniform_vec4(rgba,idxi=iunif(iu_vcolor))
    elseif (present(idx)) then
       ridx = transfer(idx,ridx)
       call setuniform_vec4(ridx,idxi=iunif(iu_idx))
    end if
    call setuniform_mat4(model,idxi=iunif(iu_model))

    ! draw
    call glDrawElements(GL_TRIANGLES, int(3*sphnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_sphere

  !> Draw a cylinder from x1 to x2 with radius rad and color rgb. ires
  !> = resolution. order = order of the bond or flat cylinder if order
  !> < 0. border = size of the border. Requires having the cylinder
  !> VAO bound.
  subroutine draw_cylinder(x1,x2,rad,rgb,ires,order,border,rgbborder)
    use interfaces_opengl3
    use tools_math, only: cross_cfloat
    use shaders, only: setuniform_vec4, setuniform_mat4, setuniform_int, setuniform_float,&
       setuniform_vec3
    use shapes, only: cylnel
    real(c_float), intent(in) :: x1(3)
    real(c_float), intent(in) :: x2(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgb(3)
    integer(c_int), intent(in) :: ires
    integer(c_int), intent(in) :: order
    real(c_float), intent(in) :: border
    real(c_float), intent(in), optional :: rgbborder(3)

    real(c_float) :: xmid(3), xdif(3), up(3), crs(3), model(4,4), blen
    real(c_float) :: a, ca, sa, axis(3), temp(3), rgb_(4), rgbborder_(3)
    integer(c_int) :: ndash

    real(c_float), parameter :: dash_length = 0.4 ! length of the dashes

    ! color of the border
    rgbborder_ = 0._c_float
    if (present(rgbborder)) rgbborder_ = rgbborder

    ! some calculations for the model matrix
    xmid = 0.5_c_float * (x1 + x2)
    xdif = x2 - x1
    blen = norm2(xdif)
    if (blen < 1e-4_c_float) return
    xdif = xdif / blen
    up = (/0._c_float,0._c_float,1._c_float/)
    crs = cross_cfloat(up,xdif)

    ! the model matrix
    model = eye4
    model(1:3,4) = xmid ! translate(m_model,xmid);
    if (dot_product(crs,crs) > 1e-14_c_float) then
       ! m_model = m_model * rotate(acos(dot(xdif,up)),crs);
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
    ! m_model = scale(m_model,glm::vec3(rad,rad,blen));
    model(:,1) = model(:,1) * rad
    model(:,2) = model(:,2) * rad
    model(:,3) = model(:,3) * blen

    ! set the uniforms
    rgb_(1:3) = rgb
    rgb_(4) = 1._c_float
    call setuniform_vec4(rgb_,idxi=iunif(iu_vcolor))
    call setuniform_mat4(model,idxi=iunif(iu_model))

    ! draw
    if (order < 0) then
       ! flat cylinder
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
    elseif (order == 0) then
       ! dashed
       ndash = nint(blen / dash_length)
       call setuniform_int(ndash,idxi=iunif(iu_ndash_cyl))
       call setuniform_float(0._c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
    elseif (order == 1) then
       ! single
       call setuniform_float(border,idxi=iunif(iu_border))
       call setuniform_vec3(rgbborder_,idxi=iunif(iu_bordercolor))
       call setuniform_float(0._c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
    elseif (order == 2) then
       ! double
       call setuniform_float(border,idxi=iunif(iu_border))
       call setuniform_vec3(rgbborder_,idxi=iunif(iu_bordercolor))
       call setuniform_float(0.75_c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
       call setuniform_float(-0.75_c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
    elseif (order == 3) then
       ! triple
       call setuniform_float(border,idxi=iunif(iu_border))
       call setuniform_vec3(rgbborder_,idxi=iunif(iu_bordercolor))
       call setuniform_float(1.35_c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
       call setuniform_float(0._c_float,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
       call setuniform_float(-1.35_c_float*rad,idxi=iunif(iu_delta_cyl))
       call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)
    end if

  end subroutine draw_cylinder

  !> Calculate the vertices for the given text and adds them to
  !> nvert/vert, direct version. (x0,y0) = position of top-left corner
  !> (pixels), siz = size in pixels. nvert/vert output vertices. If
  !> centered, center the text in x and y.
  subroutine calc_text_direct_vertices(text,x0,y0,siz,nvert,vert,centered)
    use interfaces_cimgui
    use gui_main, only: g, fontbakesize_large
    use types, only: realloc
    use param, only: newline
    character(len=*), intent(in) :: text
    real(c_float), intent(in) :: x0, y0
    real(c_float), intent(in) :: siz
    integer(c_int), intent(inout) :: nvert
    real(c_float), allocatable, intent(inout) :: vert(:,:)
    logical, intent(in), optional :: centered

    integer :: i, nline
    type(c_ptr) :: cptr
    type(ImFontGlyph), pointer :: glyph
    real(c_float) :: xpos, ypos, scale, lheight, fs
    real(c_float) :: x1, x2, y1, y2, u1, v1, u2, v2
    logical :: centered_
    real(c_float), allocatable :: xlen(:)
    integer, allocatable :: jlen(:)

    ! ibits(glyph%colored_visible_codepoint,0,1) ! colored
    ! ibits(glyph%colored_visible_codepoint,1,1) ! visible
    ! ibits(glyph%colored_visible_codepoint,2,30), ichar('R') ! codepoint
    ! glyph%AdvanceX
    ! glyph%X0, glyph%Y0, glyph%X1, glyph%Y1
    ! glyph%U0, glyph%V0, glyph%U1, glyph%V1

    ! initialize
    centered_ = .false.
    if (present(centered)) centered_ = centered
    if (.not.allocated(vert)) then
       allocate(vert(4,100))
       nvert = 0
    end if

    ! initial variables
    xpos = floor(x0)
    ypos = floor(y0)
    fs = fontbakesize_large
    scale = siz / fs
    lheight = scale * fs
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
          xpos = floor(x0)
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

       ! calculate quad and texture coordinates
       x1 = xpos + glyph%X0 * scale
       x2 = xpos + glyph%X1 * scale
       y1 = ypos + glyph%Y0 * scale
       y2 = ypos + glyph%Y1 * scale
       u1 = glyph%U0
       v1 = glyph%V1
       u2 = glyph%U1
       v2 = glyph%V0

       ! add to the vertices
       if (nvert+6 > size(vert,2)) call realloc(vert,4,2*(nvert+6))
       vert(:,nvert+1) = (/x1, y2, u1, v1/)
       vert(:,nvert+2) = (/x1, y1, u1, v2/)
       vert(:,nvert+3) = (/x2, y1, u2, v2/)
       vert(:,nvert+4) = (/x1, y2, u1, v1/)
       vert(:,nvert+5) = (/x2, y1, u2, v2/)
       vert(:,nvert+6) = (/x2, y2, u2, v1/)
       nvert = nvert + 6

       ! advance xpos
       xpos = xpos + glyph%AdvanceX * scale

       ! update
       if (centered_) then
          xlen(nline) = max(xlen(nline),xpos)
       end if
    end do
    if (centered_) then
       jlen(nline+1) = nvert+1

       do i = 1, nline
          vert(1,jlen(i):jlen(i+1)-1) = vert(1,jlen(i):jlen(i+1)-1) - 0.5_c_float * xlen(i)
       end do
       vert(2,jlen(1):nvert) = vert(2,jlen(1):nvert) - 0.5_c_float * nline * lheight
    end if

  end subroutine calc_text_direct_vertices

  !> Calculate the vertices for the given text and adds them to
  !> nvert/vert, on-scene version. x0 = world position of the label.
  !> r = radius of the associated atom.
  subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,shift,centered)
    use interfaces_cimgui
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

    integer :: i, j, nline, nvert0
    type(c_ptr) :: cptr
    type(ImFontGlyph), pointer :: glyph
    real(c_float) :: xpos, ypos, lheight, shift_(3)
    logical :: centered_
    real(c_float), allocatable :: xlen(:)
    integer, allocatable :: jlen(:)

    real(c_float), parameter :: rshift = 0.01_c_float

    ! initialize
    centered_ = .false.
    if (present(centered)) centered_ = centered
    shift_ = 0._c_float
    if (present(shift)) shift_ = shift / real(bohrtoa,c_float)
    if (.not.allocated(vert)) then
       allocate(vert(10,100))
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
       if (nvert+6 > size(vert,2)) call realloc(vert,8,2*(nvert+6))
       do j = nvert+1, nvert+6
          vert(1:3,j) = x0
          vert(4:5,j) = shift_(1:2)
          vert(6,j) = r + rshift + shift_(3)
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
    vert(7:8,nvert0+1:nvert) = vert(7:8,nvert0+1:nvert) * siz

  end subroutine calc_text_onscene_vertices

end submodule proc
