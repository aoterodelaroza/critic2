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
  ! subroutine draw_sphere(x0,rad,ires,rgb,index)
  ! subroutine draw_cylinder(x1,x2,rad,rgb,ires)
  ! subroutine calc_text_direct_vertices(text,x0,y0,siz,nvert,vert,centered)
  ! subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,centered)

contains

  !xx! scene

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use gui_main, only: nsys, sys, sysc, lockbehavior
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return

    ! basic variables
    s%id = isys
    s%isinit = 1
    s%nc = 1
    s%scenerad = 10d0
    s%scenecenter = 0d0
    s%scenexmin = 0d0
    s%scenexmax = 1d0
    s%forcebuildlists = .true.

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
    allocate(s%icount(0:reptype_NUM))
    s%icount = 0
    if (allocated(s%iord)) deallocate(s%iord)
    allocate(s%iord(20))
    s%iord = 0

    ! atoms
    s%nrep = s%nrep + 1
    call s%rep(s%nrep)%init(s,s%id,s%nrep,reptype_atoms,s%style)

    ! unit cell
    if (.not.sys(isys)%c%ismolecule) then
       s%nrep = s%nrep + 1
       call s%rep(s%nrep)%init(s,s%id,s%nrep,reptype_unitcell,s%style)
    end if

    ! reset the camera later
    s%camresetdist = 1.5_c_float
    s%camratio = 2.5_c_float

    ! sort the representations next pass
    s%forcesort = .true.
    s%timelastrender = 0d0
    s%timelastbuild = 0d0

    ! locking group for the camera
    s%forceresetcam = .true.
    if (lockbehavior == 0) then ! no-lock
       s%lockedcam = 0
    elseif (lockbehavior == 2) then ! all-lock
       s%lockedcam = 1
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
    s%id = 0
    if (allocated(s%rep)) deallocate(s%rep)
    if (allocated(s%icount)) deallocate(s%icount)
    if (allocated(s%iord)) deallocate(s%iord)
    s%nrep = 0

  end subroutine scene_end

  !> Reset the camera position and direction. Sets scenerad,
  !> scenecenter, ortho_fov, persp_fov, v_center, v_up, v_pos, view,
  !> world, projection, and znear.
  module subroutine scene_reset(s)
    use param, only: pi
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

  end subroutine scene_reset

  !> Build the draw lists for the current scene.
  module subroutine scene_build_lists(s)
    use utils, only: translate
    use gui_main, only: time
    class(scene), intent(inout), target :: s

    integer :: i
    real(c_float) :: xmin(3), xmax(3), maxrad, xc(3)

    ! initialize
    s%nsph = 0
    if (allocated(s%drawlist_sph)) deallocate(s%drawlist_sph)
    allocate(s%drawlist_sph(100))
    s%ncyl = 0
    if (allocated(s%drawlist_cyl)) deallocate(s%drawlist_cyl)
    allocate(s%drawlist_cyl(100))
    s%ncylflat = 0
    if (allocated(s%drawlist_cylflat)) deallocate(s%drawlist_cylflat)
    allocate(s%drawlist_cylflat(10))
    s%nstring = 0
    if (allocated(s%drawlist_string)) deallocate(s%drawlist_string)
    allocate(s%drawlist_string(10))

    ! add the items by representation
    do i = 1, s%nrep
       call s%rep(i)%add_draw_elements(s%nc,s%nsph,s%drawlist_sph,s%ncyl,s%drawlist_cyl,&
          s%ncylflat,s%drawlist_cylflat,s%nstring,s%drawlist_string,s%nmsel)
    end do

    ! recalculate scene center and radius
    maxrad = 0._c_float
    do i = 1, s%nrep
       if (s%rep(i)%shown .and. s%rep(i)%type == reptype_atoms) then
          maxrad = max(maxrad,maxval(s%rep(i)%atom_style(1:s%rep(i)%natom_style)%rad))
       end if
    end do
    if (s%nsph + s%ncyl + s%ncylflat + s%nstring > 0) then
       do i = 1, 3
          xmin(i) = huge(1._c_float)
          xmax(i) = -huge(1._c_float)

          xmin(i) = minval(s%drawlist_sph(1:s%nsph)%x(i)) - maxrad
          xmin(i) = min(xmin(i),minval(s%drawlist_cyl(1:s%ncyl)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cyl(1:s%ncyl)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cylflat(1:s%ncylflat)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cylflat(1:s%ncylflat)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_string(1:s%nstring)%x(i)))

          xmax(i) = maxval(s%drawlist_sph(1:s%nsph)%x(i)) + maxrad
          xmax(i) = max(xmax(i),maxval(s%drawlist_cyl(1:s%ncyl)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cyl(1:s%ncyl)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cylflat(1:s%ncylflat)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cylflat(1:s%ncylflat)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_string(1:s%nstring)%x(i)))
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

    if (s%forceresetcam) then
       ! reset the camera if requested
       call s%reset()
       s%forceresetcam = .false.
    elseif (s%lockedcam == 0) then
       ! translate the scene so the center position remains unchanged
       call translate(s%world,-xc)
    end if

    ! rebuilding lists is done
    s%forcebuildlists = .false.
    s%isinit = 2
    s%timelastbuild = time

  end subroutine scene_build_lists

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_cimgui
    use interfaces_opengl3
    use shapes, only: sphVAO, cylVAO, textVAOos, textVBOos
    use gui_main, only: fonts, fontbakesize_large, time, font_large
    use utils, only: ortho, project
    use tools_math, only: eigsym, matinv_cfloat
    use tools_io, only: string
    use shaders, only: shader_phong, shader_simple, shader_text_onscene,&
       useshader, setuniform_int, setuniform_float, setuniform_vec3,&
       setuniform_vec4, setuniform_mat3, setuniform_mat4
    class(scene), intent(inout), target :: s

    integer :: i, j
    real(c_float) :: siz, hside, xsel(3,4), radsel(4)
    integer(c_int) :: nvert
    real(c_float), allocatable, target :: vert(:,:)

    real(c_float), parameter :: rgbsel(4,4) = reshape((/&
       1._c_float,  0.4_c_float, 0.4_c_float, 0.5_c_float,&
       0.4_c_float, 1._c_float,  0.4_c_float, 0.5_c_float,&
       0.4_c_float, 0.4_c_float, 1._c_float,  0.5_c_float,&
       0.9_c_float, 0.7_c_float, 0.4_c_float, 0.5_c_float/),shape(rgbsel))
    real(c_float), parameter :: msel_thickness = 0.1_c_float
    real(c_float), parameter :: sel_label_size = 1.2_c_float

    ! check that the scene and system are initialized
    if (s%isinit == 0) return

    ! build draw lists if not done already
    if (s%isinit == 1 .or. .not.allocated(s%drawlist_sph)) call s%build_lists()

    ! render text with the large font
    call igPushFont(font_large)

    ! if necessary, rebuild draw lists
    if (s%forcebuildlists) call s%build_lists()

    if (s%style == style_phong) then
       !! phong !!
       ! set up the shader and the uniforms
       call useshader(shader_phong)
       call setuniform_int("uselighting",1_c_int)
       call setuniform_vec3("lightPos",s%lightpos)
       call setuniform_vec3("lightColor",s%lightcolor)
       call setuniform_float("ambient",s%ambient)
       call setuniform_float("diffuse",s%diffuse)
       call setuniform_float("specular",s%specular)
       call setuniform_int("shininess",s%shininess)
       call setuniform_mat4("world",s%world)
       call setuniform_mat4("view",s%view)
       call setuniform_mat4("projection",s%projection)

       ! draw the atoms
       if (s%nsph > 0) then
          call glBindVertexArray(sphVAO(s%atom_res))
          do i = 1, s%nsph
             call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r,s%atom_res,rgb=s%drawlist_sph(i)%rgb)
          end do
       end if

       ! draw the bonds
       if (s%ncyl > 0) then
          call glBindVertexArray(cylVAO(s%bond_res))
          do i = 1, s%ncyl
             call draw_cylinder(s%drawlist_cyl(i)%x1,s%drawlist_cyl(i)%x2,s%drawlist_cyl(i)%r,&
                s%drawlist_cyl(i)%rgb,s%bond_res)
          end do
       end if

       ! draw the flat cylinders (unit cell)
       if (s%ncylflat > 0) then
          call setuniform_int("uselighting",0_c_int)
          call glBindVertexArray(cylVAO(s%uc_res))
          do i = 1, s%ncylflat
             call draw_cylinder(s%drawlist_cylflat(i)%x1,s%drawlist_cylflat(i)%x2,&
                s%drawlist_cylflat(i)%r,s%drawlist_cylflat(i)%rgb,s%uc_res)
          end do
       end if

       ! draw the selected atoms
       if (s%nmsel > 0) then
          call setuniform_int("uselighting",0_c_int)
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          do i = 1, s%nsph
             do j = 1, s%nmsel
                if (all(s%drawlist_sph(i)%idx == s%msel(:,j))) then
                   call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r + msel_thickness,s%atom_res,rgba=rgbsel(:,j))
                   radsel(j) = s%drawlist_sph(i)%r + msel_thickness
                   xsel(:,j) = s%drawlist_sph(i)%x
                end if
             end do
          end do
          call glDisable(GL_BLEND)
       end if

       ! render labels with on-scene text
       call useshader(shader_text_onscene)
       call setuniform_mat4("world",s%world)
       call setuniform_mat4("view",s%view)
       call setuniform_mat4("projection",s%projection)

       call glDisable(GL_MULTISAMPLE)
       call glEnable(GL_BLEND)
       call glBlendEquation(GL_FUNC_ADD)
       call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)

       call glActiveTexture(GL_TEXTURE0)
       call glBindVertexArray(textVAOos)
       call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))
       call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)

       do i = 1, s%nstring
          call setuniform_vec3("textColor",s%drawlist_string(i)%rgb)
          nvert = 0
          if (s%drawlist_string(i)%scale > 0._c_float) then
             hside = s%camresetdist * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
             hside = hside * s%camratio
             hside = max(hside,3._c_float)
             siz = 2 * s%drawlist_string(i)%scale / fontbakesize_large / hside
          else
             siz = 2 * abs(s%drawlist_string(i)%scale) * s%projection(1,1) / fontbakesize_large
          end if
          call calc_text_onscene_vertices(s%drawlist_string(i)%str,s%drawlist_string(i)%x,s%drawlist_string(i)%r,&
             siz,nvert,vert,centered=.true.)
          call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*8*c_sizeof(c_float), c_loc(vert))
          call glDrawArrays(GL_TRIANGLES, 0, nvert)
       end do

       ! render selected atom labels with on-scene text
       if (s%nmsel > 0) then
          do j = 1, s%nmsel
             call setuniform_vec3("textColor",(/1._c_float,1._c_float,1._c_float/))
             siz = sel_label_size * s%projection(1,1) / fontbakesize_large
             nvert = 0
             call calc_text_onscene_vertices(string(j),xsel(:,j),radsel(j),siz,nvert,vert,centered=.true.)
             call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*8*c_sizeof(c_float), c_loc(vert))
             call glDrawArrays(GL_TRIANGLES, 0, nvert)
          end do
       end if
       call glEnable(GL_MULTISAMPLE)
       call glDisable(GL_BLEND)
    else
       !! simple !!
       ! set up the shader and the uniforms
       call useshader(shader_simple)
       call setuniform_mat4("world",s%world)
       call setuniform_mat4("view",s%view)
       call setuniform_mat4("projection",s%projection)

       ! draw the atoms
       call setuniform_float("rborder",s%atomborder)
       call setuniform_vec3("bordercolor",s%bordercolor)
       if (s%nsph > 0) then
          call glBindVertexArray(sphVAO(s%atom_res))
          do i = 1, s%nsph
             call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r,s%atom_res,rgb=s%drawlist_sph(i)%rgb)
          end do
       end if

       ! draw the bonds
       call setuniform_float("rborder",-1._c_float)
       if (s%ncyl > 0) then
          call glBindVertexArray(cylVAO(s%bond_res))
          do i = 1, s%ncyl
             call draw_cylinder(s%drawlist_cyl(i)%x1,s%drawlist_cyl(i)%x2,s%drawlist_cyl(i)%r,&
                s%drawlist_cyl(i)%rgb,s%bond_res)
          end do
       end if

       ! draw the flat cylinders (unit cell)
       if (s%ncylflat > 0) then
          call glBindVertexArray(cylVAO(s%uc_res))
          do i = 1, s%ncylflat
             call draw_cylinder(s%drawlist_cylflat(i)%x1,s%drawlist_cylflat(i)%x2,&
                s%drawlist_cylflat(i)%r,s%drawlist_cylflat(i)%rgb,s%uc_res)
          end do
       end if

       ! draw the selected atoms
       if (s%nmsel > 0) then
          call glBindVertexArray(sphVAO(s%atom_res))
          call glEnable(GL_BLEND)
          call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
          do i = 1, s%nsph
             do j = 1, s%nmsel
                if (all(s%drawlist_sph(i)%idx == s%msel(:,j))) then
                   call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r + msel_thickness,s%atom_res,rgba=rgbsel(:,j))
                   radsel(j) = s%drawlist_sph(i)%r + msel_thickness
                   xsel(:,j) = s%drawlist_sph(i)%x
                end if
             end do
          end do
          call glDisable(GL_BLEND)
       end if

       ! render labels with on-scene text
       call useshader(shader_text_onscene)
       call setuniform_mat4("world",s%world)
       call setuniform_mat4("view",s%view)
       call setuniform_mat4("projection",s%projection)

       call glDisable(GL_MULTISAMPLE)
       call glEnable(GL_BLEND)
       call glBlendEquation(GL_FUNC_ADD)
       call glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA)

       call glActiveTexture(GL_TEXTURE0)
       call glBindVertexArray(textVAOos)
       call glBindTexture(GL_TEXTURE_2D, transfer(fonts%TexID,1_c_int))
       call glBindBuffer(GL_ARRAY_BUFFER, textVBOos)

       do i = 1, s%nstring
          call setuniform_vec3("textColor",s%drawlist_string(i)%rgb)
          nvert = 0
          if (s%drawlist_string(i)%scale > 0._c_float) then
             hside = s%camresetdist * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
             hside = hside * s%camratio
             hside = max(hside,3._c_float)
             siz = 2 * s%drawlist_string(i)%scale / fontbakesize_large / hside
          else
             siz = 2 * abs(s%drawlist_string(i)%scale) * s%projection(1,1) / fontbakesize_large
          end if
          call calc_text_onscene_vertices(s%drawlist_string(i)%str,s%drawlist_string(i)%x,s%drawlist_string(i)%r,&
             siz,nvert,vert,centered=.true.)
          call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*8*c_sizeof(c_float), c_loc(vert))
          call glDrawArrays(GL_TRIANGLES, 0, nvert)
       end do

       ! render selected atom labels with on-scene text
       if (s%nmsel > 0) then
          do j = 1, s%nmsel
             call setuniform_vec3("textColor",(/1._c_float,1._c_float,1._c_float/))
             siz = sel_label_size * s%projection(1,1) / fontbakesize_large
             nvert = 0
             call calc_text_onscene_vertices(string(j),xsel(:,j),radsel(j),siz,nvert,vert,centered=.true.)
             call glBufferSubData(GL_ARRAY_BUFFER, 0_c_intptr_t, nvert*8*c_sizeof(c_float), c_loc(vert))
             call glDrawArrays(GL_TRIANGLES, 0, nvert)
          end do
       end if
       call glEnable(GL_MULTISAMPLE)
       call glDisable(GL_BLEND)
    end if

    ! pop the large font
    call igPopFont()

    call glBindBuffer(GL_ARRAY_BUFFER, 0)
    call glBindVertexArray(0)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! save the rendering time
    s%timelastrender = time

  end subroutine scene_render

  !> Draw the scene (for object picking
  module subroutine scene_render_pick(s)
    use interfaces_cimgui
    use interfaces_opengl3
    use shapes, only: sphVAO
    use utils, only: ortho, project
    use tools_math, only: eigsym, matinv_cfloat
    use shaders, only: shader_pickindex, useshader, setuniform_int,&
       setuniform_float, setuniform_vec3, setuniform_vec4, setuniform_mat3,&
       setuniform_mat4
    class(scene), intent(inout), target :: s

    integer :: i

    ! check that the scene and system are initialized
    if (s%isinit < 2) return

    ! build draw lists if not done already
    if (.not.allocated(s%drawlist_sph)) call s%build_lists()

    ! if necessary, rebuild draw lists
    if (s%forcebuildlists) call s%build_lists()

    ! set up the shader and the uniforms
    call useshader(shader_pickindex)
    call setuniform_mat4("world",s%world)
    call setuniform_mat4("view",s%view)
    call setuniform_mat4("projection",s%projection)

    ! draw the atoms
    if (s%nsph > 0) then
       call glBindVertexArray(sphVAO(s%atom_res))
       do i = 1, s%nsph
          call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r,s%atom_res,idx=s%drawlist_sph(i)%idx)
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
       s%bgcolor = (/0._c_float,0._c_float,0._c_float/)
       s%lightpos = (/20._c_float,20._c_float,0._c_float/)
       s%lightcolor = (/1._c_float,1._c_float,1._c_float/)
       s%ambient = 0.2_c_float
       s%diffuse = 0.4_c_float
       s%specular = 0.6_c_float
       s%shininess = 8_c_int
       do i = 1, s%nrep
          s%rep(i)%label_rgb = 1._c_float
          s%rep(i)%uc_rgb = 1._c_float
          s%rep(i)%bond_color_style=1
       end do
    else
       !! simple !!
       s%bgcolor = (/1._c_float,1._c_float,1._c_float/)
       do i = 1, s%nrep
          s%rep(i)%label_rgb = 0._c_float
          s%rep(i)%uc_rgb = 0._c_float
          s%rep(i)%bond_color_style=0
          s%rep(i)%bond_rgb = 0._c_float
       end do
       s%atomborder = 0.2_c_float
       s%bordercolor = 0._c_float
    end if

  end subroutine scene_set_style_defaults

  !> Copy camera parameters from scene si to the current scene. If an
  !> integer is given instead, search the system list for the most
  !> recent render in group idx and, if found, copy the camera
  !> parameters from that scene.
  recursive module subroutine scene_copy_cam(s,si,idx)
    use gui_main, only: nsys, sysc
    class(scene), intent(inout), target :: s
    type(scene), intent(in), target, optional :: si
    integer, intent(in), optional :: idx

    real*8 :: time
    integer :: i, ifound

    if (present(si)) then
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
    elseif (present(idx)) then
       time = 0d0
       ifound = 0
       do i = 1, nsys
          if (sysc(i)%sc%lockedcam == idx .and. sysc(i)%sc%timelastrender > time) then
             ifound = i
             time = sysc(i)%sc%timelastrender
          end if
       end do
       if (ifound > 0) call s%copy_cam(sysc(ifound)%sc)
    end if

  end subroutine scene_copy_cam

  !> Show the representation menu (called from view). Return .true.
  !> if the scene needs to be rendered again.
  module function representation_menu(s,idcaller) result(changed)
    use interfaces_cimgui
    use utils, only: iw_text, iw_tooltip, iw_button
    use windows, only: win, stack_create_window, wintype_editrep, update_window_id
    use gui_main, only: ColorDangerButton, g
    use tools_io, only: string
    use tools, only: mergesort
    class(scene), intent(inout), target :: s
    integer(c_int), intent(in) :: idcaller
    logical :: changed

    integer :: i, ii, id, ll
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    logical :: discol, doerase
    type(ImVec2) :: szero
    character(kind=c_char,len=1024), target :: txtinp
    integer, allocatable :: idx(:)

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with draw_view in windows@view module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_editbutton = 3

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

       ! update window ID
       call update_window_id(s%rep(i)%idwin)

       ! close button
       doerase = .false.
       call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
       if (igTableSetColumnIndex(ic_closebutton)) then
          str1 = "##2ic_closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
          if (my_CloseButton(c_loc(str1),ColorDangerButton)) doerase = .true.
       end if

       ! view button
       if (igTableSetColumnIndex(ic_viewbutton)) then
          str1 = "##2ic_viewbutton" // string(ic_viewbutton) // "," // string(i) // c_null_char
          if (igCheckbox(c_loc(str1), s%rep(i)%shown)) &
             changed = .true.
       end if

       ! name
       if (igTableSetColumnIndex(ic_name)) then
          discol = .not.s%rep(i)%shown
          if (discol) &
             call igPushStyleColor_Vec4(ImGuiCol_Text,g%Style%Colors(ImGuiCol_TextDisabled+1))
          call iw_text(trim(s%rep(i)%name))
          if (discol) call igPopStyleColor(1)

          ! name context menu
          if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_MouseButtonRight)) then
             ! edit
             str2 = "Edit" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                if (s%rep(i)%idwin == 0) then
                   s%rep(i)%idwin = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idcaller=idcaller)
                else
                   call igSetWindowFocus_Str(c_loc(win(s%rep(i)%idwin)%name))
                end if
             end if
             call iw_tooltip("Edit this representation",ttshown)

             ! duplicate
             str2 = "Duplicate" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = s%get_new_representation_id()
                s%rep(id) = s%rep(i)
                s%rep(id)%name = trim(s%rep(i)%name) // " (copy)"
                s%icount(s%rep(i)%type) = s%icount(s%rep(i)%type) + 1
                s%icount(0) = s%icount(0) + 1
                s%rep(id)%iord = s%icount(0)
                s%forcesort = .true.
                changed = .true.
             end if
             call iw_tooltip("Make a copy of this representation",ttshown)

             ! show/hide
             str2 = "Show/Hide" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                s%rep(i)%shown = .not.s%rep(i)%shown
                changed = .true.
             end if
             call iw_tooltip("Toggle hide/show this representation",ttshown)

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
             call iw_tooltip("Rename this representation",ttshown)

             ! delete
             str2 = "Delete" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) &
                doerase = .true.
             call iw_tooltip("Delete this representation",ttshown)

             call igEndPopup()
          end if
       end if

       ! edit button
       if (igTableSetColumnIndex(ic_editbutton)) then
          if (iw_button("Edit##2ic_editbutton" // string(ic_editbutton) // "," // string(i))) then
             if (s%rep(i)%idwin == 0) then
                s%rep(i)%idwin = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idcaller=idcaller)
             else
                call igSetWindowFocus_Str(c_loc(win(s%rep(i)%idwin)%name))
             end if
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
    use gui_main, only: sys
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
    integer, intent(in) :: idx(4)

    integer :: i, j

    ! if no atom, clear
    if (idx(1) == 0) then
       s%nmsel = 0
       return
    end if

    ! if the atom is already selected, deselect
    do i = 1, s%nmsel
       if (all(idx == s%msel(:,i))) then
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

  !xx! representation

  !> Initialize a representation. If itype is present and not _none,
  !> fill the representation with the default values for the
  !> corresponding type and set isinit = .true. and shown = .true.  sc
  !> = parent scene, isys = system ID, irep = representation ID, style
  !> = phong or simple.
  module subroutine representation_init(r,sc,isys,irep,itype,style)
    use gui_main, only: sys
    use tools_io, only: string
    class(representation), intent(inout), target :: r
    type(scene), intent(inout), target :: sc
    integer, intent(in) :: isys
    integer, intent(in) :: irep
    integer, intent(in) :: itype
    integer, intent(in) :: style

    ! common settings
    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%id = isys
    r%idrep = irep
    r%idwin = 0
    r%name = ""
    r%filter = ""
    r%errfilter = ""
    r%pertype = 1
    r%ncell = 1
    r%border = .true.
    r%onemotif = .false.
    r%atoms_display = .true.
    r%bonds_display = .true.
    r%labels_display = .false.
    r%atom_style_type = 0
    r%natom_style = 0
    r%atom_radii_reset_type = 0
    r%atom_radii_reset_scale = 0.7_c_float
    r%atom_color_reset_type = 0
    r%bond_rad = 0.2_c_float
    r%label_style = 0
    r%label_scale = 0.5_c_float
    r%label_const_size = .false._c_bool
    r%label_exclude_h = .true._c_bool
    r%uc_radius = 0.15_c_float
    r%uc_radiusinner = 0.15_c_float
    r%uc_innersteplen = 2d0
    r%uc_innerstipple = .true.
    r%uc_inner = .true.
    r%uc_coloraxes = .true.
    r%origin = 0._c_float

    ! style-dependent settings
    if (style == style_phong) then
       r%label_rgb = 1._c_float
       r%bond_color_style = 1
       r%bond_rgb = (/1._c_float,0._c_float,0._c_float/)
       r%uc_rgb = 1._c_float
    else
       r%label_rgb = 0._c_float
       r%bond_color_style = 0
       r%bond_rgb = 0._c_float
       r%uc_rgb = 0._c_float
    end if

    ! type-dependent settings
    if (itype == reptype_atoms) then
       r%isinit = .true.
       r%shown = .true.
       r%type = reptype_atoms
       r%name = "Atoms + Bonds + Labels"
       if (sys(isys)%c%ismolecule) then
          r%ncell = 0
          r%border = .false.
          r%onemotif = .false.
       else
          r%border = .true.
          r%onemotif = .true.
          r%ncell = 1
       end if
    elseif (itype == reptype_unitcell) then
       r%isinit = .true.
       r%shown = .true.
       r%type = reptype_unitcell
       r%name = "Unit cell"
    end if

    ! increment type counter and set name
    sc%icount(itype) = sc%icount(itype) + 1
    if (sc%icount(itype) > 1) then
       r%name = trim(r%name) // " " // string(sc%icount(itype))
    end if

    ! increment global counter and force sort of the parent scene
    sc%icount(0) = sc%icount(0) + 1
    r%iord = sc%icount(0)
    sc%forcesort = .true.

    ! initialize the styles
    call r%reset_atom_style()

  end subroutine representation_init

  !> Terminate a representation
  module subroutine representation_end(r)
    use windows, only: win
    class(representation), intent(inout), target :: r

    r%name = ""
    r%filter = ""
    r%errfilter = ""
    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%id = 0
    r%idrep = 0
    r%iord = 0
    if (r%idwin > 0) &
       nullify(win(r%idwin)%rep)
    r%idwin = 0

  end subroutine representation_end

  !> Reset atom styles.
  module subroutine reset_atom_style(r)
    use gui_main, only: nsys, sys
    use param, only: jmlcol, atmcov
    class(representation), intent(inout), target :: r

    integer :: i, ispc, iz

    ! set the atom style to zero
    r%natom_style = 0
    if (allocated(r%atom_style)) deallocate(r%atom_style)

    ! check the system is correct and initialized
    if (r%id < 1 .or. r%id > nsys) return

    ! fill the styles
    if (r%atom_style_type == 0) then
       ! species
       r%natom_style = sys(r%id)%c%nspc
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%nspc
          r%atom_style(i)%shown = .true.

          iz = sys(r%id)%c%spc(i)%z
          r%atom_style(i)%rgb = real(jmlcol(:,iz),c_float) / 255._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    elseif (r%atom_style_type == 1) then
       ! nneq
       r%natom_style = sys(r%id)%c%nneq
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%nneq
          r%atom_style(i)%shown = .true.

          ispc = sys(r%id)%c%at(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          r%atom_style(i)%rgb = real(jmlcol(:,iz),c_float) / 255._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    elseif (r%atom_style_type == 2) then
       ! ncel
       r%natom_style = sys(r%id)%c%ncel
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%ncel
          r%atom_style(i)%shown = .true.

          ispc = sys(r%id)%c%atcel(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          r%atom_style(i)%rgb = real(jmlcol(:,iz),c_float) / 255._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    end if

  end subroutine reset_atom_style

  !> Add the spheres, cylinder, etc. to the draw lists.
  module subroutine add_draw_elements(r,nc,nsph,drawlist_sph,ncyl,drawlist_cyl,&
     ncylflat,drawlist_cylflat,nstring,drawlist_string,nmsel)
    use gui_main, only: sys
    use tools_io, only: string, nameguess
    use hashmod, only: hash
    use param, only: bohrtoa, newline
    class(representation), intent(inout), target :: r
    integer, intent(in) :: nc(3)
    integer, intent(inout) :: nsph
    type(dl_sphere), intent(inout), allocatable :: drawlist_sph(:)
    integer, intent(inout) :: ncyl
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cyl(:)
    integer, intent(inout) :: ncylflat
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cylflat(:)
    integer, intent(inout) :: nstring
    type(dl_string), intent(inout), allocatable :: drawlist_string(:)
    integer, intent(in) :: nmsel

    type(hash) :: shown_atoms
    logical :: havefilter, step, isedge(3)
    integer :: n(3), i, j, k, imol, lvec(3), id, idaux, n0(3), n1(3), i1, i2, i3, ix(3)
    integer :: ib, ineigh, ixn(3), ix1(3), ix2(3), nstep, idx
    real(c_float) :: rgb(3), rad
    real*8 :: xx(3), x0(3), x1(3), x2(3), res, uoriginc(3)
    type(dl_sphere), allocatable :: auxsph(:)
    type(dl_cylinder), allocatable :: auxcyl(:)
    type(dl_string), allocatable :: auxstr(:)
    character(len=:), allocatable :: atcode, errmsg

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr
    integer, parameter :: uc(3,2,12) = reshape((/&
       0,0,0,  1,0,0,&
       0,0,0,  0,1,0,&
       0,0,0,  0,0,1,&
       1,1,0,  1,1,1,&
       1,0,1,  1,1,1,&
       0,1,1,  1,1,1,&
       1,0,0,  1,1,0,&
       0,1,0,  1,1,0,&
       0,1,0,  0,1,1,&
       0,0,1,  0,1,1,&
       0,0,1,  1,0,1,&
       1,0,0,  1,0,1/),shape(uc))
    integer, parameter :: ucdir(12) = (/1, 2, 3, 3, 2, 1, 2, 1, 3, 2, 1, 3/)

    ! return if not initialized
    if (.not.r%isinit) return
    if (.not.r%shown) return

    ! initialize the drawlists if not done already
    if (.not.allocated(drawlist_sph)) then
       allocate(drawlist_sph(100))
       nsph = 0
    end if
    if (.not.allocated(drawlist_cyl)) then
       allocate(drawlist_cyl(100))
       ncyl = 0
    end if
    if (.not.allocated(drawlist_cylflat)) then
       allocate(drawlist_cylflat(100))
       ncylflat = 0
    end if
    if (.not.allocated(drawlist_string)) then
       allocate(drawlist_string(100))
       nstring = 0
    end if

    if (r%type == reptype_atoms) then
       !!! atoms and bonds representation !!!

       !! first, the atoms
       ! do we have a filter?
       havefilter = (len_trim(r%filter) > 0) .and. (len_trim(r%errfilter) == 0)

       ! calculate the periodicity
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! origin shift
       if (sys(r%id)%c%ismolecule) then
          uoriginc = r%origin / bohrtoa
       else
          uoriginc = sys(r%id)%c%x2c(real(r%origin,8))
       end if

       ! run over atoms, either directly or per-molecule
       i = 0
       imol = 0
       do while(.true.)
          if (r%onemotif) then
             ! this is a new molecule if there are no molecules or this is the last atom
             ! in the previous one
             step = (imol == 0)
             if (.not.step) step = (k == sys(r%id)%c%mol(imol)%nat)
             if (step) then
                imol = imol + 1
                k = 0
             end if

             ! we are finished if we have all molecules
             if (imol > sys(r%id)%c%nmol) exit

             ! Add the new atom; only translate if the fragment is discrete
             k = k + 1
             i = sys(r%id)%c%mol(imol)%at(k)%cidx
             if (sys(r%id)%c%mol(imol)%discrete) then
                lvec = sys(r%id)%c%mol(imol)%at(k)%lvec
             else
                lvec = 0
             end if
          else
             i = i + 1
             if (i > sys(r%id)%c%ncel) exit
             lvec = 0
          end if

          ! skip hidden atoms
          if (r%atom_style_type == 0) then
             id = sys(r%id)%c%atcel(i)%is
          elseif (r%atom_style_type == 1) then
             id = sys(r%id)%c%atcel(i)%idx
          else
             id = i
          end if
          if (.not.r%atom_style(id)%shown) cycle

          ! calculate the border
          n0 = 0
          n1 = n-1
          if (r%border.and..not.r%onemotif) then
             xx = sys(r%id)%c%atcel(i)%x
             do j = 1, 3
                if (xx(j) < rthr) then
                   n1(j) = n(j)
                elseif (xx(j) > rthr1) then
                   n0(j) = -1
                end if
             end do
          end if

          ! draw the spheres and cylinders
          rgb = r%atom_style(id)%rgb
          rad = r%atom_style(id)%rad
          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   ! atoms
                   ix = (/i1,i2,i3/) + lvec
                   xx = sys(r%id)%c%atcel(i)%x + ix
                   xx = sys(r%id)%c%x2c(xx)

                   ! apply the filter
                   if (havefilter) then
                      res = sys(r%id)%eval(r%filter,errmsg,xx)
                      if (len_trim(errmsg) == 0) then
                         if (res == 0d0) cycle
                      else
                         havefilter = .false.
                         r%errfilter = errmsg
                      end if
                   end if

                   ! draw the atom, reallocate if necessary
                   if (r%atoms_display) then
                      nsph = nsph + 1
                      if (nsph > size(drawlist_sph,1)) then
                         allocate(auxsph(2*nsph))
                         auxsph(1:size(drawlist_sph,1)) = drawlist_sph
                         call move_alloc(auxsph,drawlist_sph)
                      end if

                      ! write down the sphere
                      drawlist_sph(nsph)%x = real(xx + uoriginc,c_float)
                      drawlist_sph(nsph)%r = rad
                      drawlist_sph(nsph)%rgb = rgb
                      drawlist_sph(nsph)%idx(1) = i
                      drawlist_sph(nsph)%idx(2:4) = ix
                   end if

                   ! bonds
                   if (r%bonds_display) then
                      ! add this atom to the hash
                      atcode = string(i) // "_" // string(ix(1)) // "_" // string(ix(2)) // "_" // string(ix(3))
                      call shown_atoms%put(atcode,1)

                      ! bonds
                      do ib = 1, sys(r%id)%c%nstar(i)%ncon
                         ineigh = sys(r%id)%c%nstar(i)%idcon(ib)
                         ixn = ix + sys(r%id)%c%nstar(i)%lcon(:,ib)

                         ! check if the atom has been represented already
                         atcode = string(ineigh) // "_" // string(ixn(1)) // "_" // string(ixn(2)) // "_" // string(ixn(3))
                         if (.not.shown_atoms%iskey(atcode)) cycle

                         ! draw the bond, reallocate if necessary
                         if (r%bond_color_style == 0) then
                            ncyl = ncyl + 1
                         else
                            ncyl = ncyl + 2
                         end if
                         if (ncyl > size(drawlist_cyl,1)) then
                            allocate(auxcyl(2*ncyl))
                            auxcyl(1:size(drawlist_cyl,1)) = drawlist_cyl
                            call move_alloc(auxcyl,drawlist_cyl)
                         end if

                         x1 = xx + uoriginc
                         x2 = sys(r%id)%c%atcel(ineigh)%x + ixn
                         x2 = sys(r%id)%c%x2c(x2) + uoriginc
                         if (r%bond_color_style == 0) then
                            drawlist_cyl(ncyl)%x1 = real(x1,c_float)
                            drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                            drawlist_cyl(ncyl)%r = r%bond_rad
                            drawlist_cyl(ncyl)%rgb = r%bond_rgb
                         else
                            x0 = 0.5d0 * (x1 + x2)
                            drawlist_cyl(ncyl-1)%x1 = real(x1,c_float)
                            drawlist_cyl(ncyl-1)%x2 = real(x0,c_float)
                            drawlist_cyl(ncyl-1)%r = r%bond_rad
                            drawlist_cyl(ncyl-1)%rgb = r%atom_style(id)%rgb

                            if (r%atom_style_type == 0) then
                               idaux = sys(r%id)%c%atcel(ineigh)%is
                            elseif (r%atom_style_type == 1) then
                               idaux = sys(r%id)%c%atcel(ineigh)%idx
                            else
                               idaux = ineigh
                            end if
                            drawlist_cyl(ncyl)%x1 = real(x0,c_float)
                            drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                            drawlist_cyl(ncyl)%r = r%bond_rad
                            drawlist_cyl(ncyl)%rgb = r%atom_style(idaux)%rgb
                         end if
                      end do ! ncon
                   end if

                   ! labels
                   if (r%labels_display .and. &
                      (.not.r%label_exclude_h.or.sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z/=1)) then
                      nstring = nstring + 1
                      if (nstring > size(drawlist_string,1)) then
                         allocate(auxstr(2*nstring))
                         auxstr(1:size(drawlist_string,1)) = drawlist_string
                         call move_alloc(auxstr,drawlist_string)
                      end if

                      drawlist_string(nstring)%x = real(xx + uoriginc,c_float)
                      drawlist_string(nstring)%r = rad
                      drawlist_string(nstring)%rgb = r%label_rgb
                      if (r%label_const_size) then
                         drawlist_string(nstring)%scale = r%label_scale
                      else
                         drawlist_string(nstring)%scale = -r%label_scale
                      end if
                      if (r%label_style == 0) then ! 0 = atomic symbol
                         drawlist_string(nstring)%str = trim(nameguess(sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z,.true.))
                      elseif (r%label_style == 1) then ! 1 = atom name
                         drawlist_string(nstring)%str = trim(sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%name)
                      elseif (r%label_style == 2) then ! 2 = cel-atom
                         drawlist_string(nstring)%str = string(i)
                      elseif (r%label_style == 3) then ! 3 = cel-atom + lvec
                         drawlist_string(nstring)%str = string(i) // newline // "(" // string(ix(1)) // "," //&
                            string(ix(2)) // "," // string(ix(3)) // ")"
                      elseif (r%label_style == 4) then ! 4 = neq atom
                         drawlist_string(nstring)%str = string(sys(r%id)%c%atcel(i)%idx)
                      elseif (r%label_style == 5) then ! 5 = spc
                         drawlist_string(nstring)%str = string(sys(r%id)%c%atcel(i)%is)
                      elseif (r%label_style == 6) then ! 6 = Z
                         drawlist_string(nstring)%str = string(sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z)
                      elseif (r%label_style == 7) then ! 7 = mol
                         drawlist_string(nstring)%str = string(sys(r%id)%c%idatcelmol(i))
                      elseif (r%label_style == 8) then ! 8 = wycoff
                         idx = sys(r%id)%c%atcel(i)%idx
                         drawlist_string(nstring)%str = string(sys(r%id)%c%at(idx)%mult) //&
                            string(sys(r%id)%c%at(idx)%wyc)
                      end if
                   end if
                end do ! i3
             end do ! i2
          end do ! i1
       end do ! loop over complete atom list (i)
    elseif (r%type == reptype_unitcell) then
       !!! unit cell representation !!!

       ! number of cells
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! external cell
       do i = 1, 12
          x1 = real(uc(:,1,i) * n,8) + r%origin
          x1 = sys(r%id)%c%x2c(x1)
          x2 = real(uc(:,2,i) * n,8) + r%origin
          x2 = sys(r%id)%c%x2c(x2)

          call increase_ncylflat()
          drawlist_cylflat(ncylflat)%x1 = real(x1,c_float)
          drawlist_cylflat(ncylflat)%x2 = real(x2,c_float)
          drawlist_cylflat(ncylflat)%r = r%uc_radius
          if (r%uc_coloraxes.and.i==1) then
             drawlist_cylflat(ncylflat)%rgb = (/1._c_float,0._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==2) then
             drawlist_cylflat(ncylflat)%rgb = (/0._c_float,1._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==3) then
             drawlist_cylflat(ncylflat)%rgb = (/0._c_float,0._c_float,1._c_float/)
          else
             drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
          end if
       end do

       ! draw inner cylinders
       if (r%uc_inner) then
          do i1 = 0, n(1)-1
             do i2 = 0, n(2)-1
                do i3 = 0, n(3)-1
                   do i = 1, 12
                      ix1 = (/i1,i2,i3/) + uc(:,1,i)
                      ix2 = (/i1,i2,i3/) + uc(:,2,i)

                      ! skip outer cylinders
                      isedge = (ix1 == 0 .and. ix2 == 0) .or. (ix1 == n .and. ix2 == n)
                      isedge(ucdir(i)) = .true.
                      if (all(isedge)) cycle

                      x1 = real(ix1,8) + r%origin
                      x1 = sys(r%id)%c%x2c(x1)
                      x2 = real(ix2,8) + r%origin
                      x2 = sys(r%id)%c%x2c(x2)

                      ! logical :: uc_innerstipple ! stippled lines for the inner lines
                      if (r%uc_innerstipple) then
                         nstep = ceiling(norm2(x2 - x1) / r%uc_innersteplen)
                         do j = 1, nstep
                            call increase_ncylflat()
                            drawlist_cylflat(ncylflat)%x1 = real(x1 + real(2*j-1,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%x2 = real(x1 + real(2*j,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                            drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
                         end do
                      else
                         call increase_ncylflat()
                         drawlist_cylflat(ncylflat)%x1 = real(x1 ,c_float)
                         drawlist_cylflat(ncylflat)%x2 = real(x2 ,c_float)
                         drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                         drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
                      end if
                   end do
                end do
             end do
          end do
       end if
    end if ! reptype
  contains
    subroutine increase_ncylflat()

      ncylflat = ncylflat + 1
      if (ncylflat > size(drawlist_cylflat,1)) then
         allocate(auxcyl(2*ncylflat))
         auxcyl(1:size(drawlist_cylflat,1)) = drawlist_cylflat
         call move_alloc(auxcyl,drawlist_cylflat)
      end if

    end subroutine increase_ncylflat
  end subroutine add_draw_elements

  !xx! private procedures: low-level draws

  !> Draw a sphere with center x0, radius rad and color rgb. Requires
  !> having the sphere VAO bound.
  subroutine draw_sphere(x0,rad,ires,rgb,rgba,idx)
    use interfaces_opengl3
    use shaders, only: setuniform_vec3, setuniform_vec4, setuniform_mat4
    use shapes, only: sphnel
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: rad
    integer(c_int), intent(in) :: ires
    real(c_float), intent(in), optional :: rgb(3)
    real(c_float), intent(in), optional :: rgba(4)
    integer(c_int), intent(in), optional :: idx(4)

    real(c_float) :: model(4,4)
    real(c_float) :: ridx(4), rgb_(4)

    ! the model matrix: scale and translate
    model = eye4
    model(1:3,4) = x0
    model(1,1) = rad
    model(2,2) = rad
    model(3,3) = rad

    ! draw the sphere
    if (present(rgb)) then
       rgb_(1:3) = rgb
       rgb_(4) = 1._c_float
       call setuniform_vec4("vColor",rgb_)
    elseif (present(rgba)) then
       call setuniform_vec4("vColor",rgba)
    elseif (present(idx)) then
       ridx = transfer(idx,ridx)
       call setuniform_vec4("idx",ridx)
    end if
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*sphnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_sphere

  !> Draw a cylinder from x1 to x2 with radius rad and color
  !> rgb. Requires having the cylinder VAO bound.
  subroutine draw_cylinder(x1,x2,rad,rgb,ires)
    use interfaces_opengl3
    use tools_math, only: cross_cfloat
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: cylnel
    real(c_float), intent(in) :: x1(3)
    real(c_float), intent(in) :: x2(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgb(3)
    integer(c_int), intent(in) :: ires

    real(c_float) :: xmid(3), xdif(3), up(3), crs(3), model(4,4), blen
    real(c_float) :: a, ca, sa, axis(3), temp(3), rgb_(4)

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

    ! draw the cylinder
    rgb_(1:3) = rgb
    rgb_(4) = 1._c_float
    call setuniform_vec4("vColor",rgb_)
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)

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
  subroutine calc_text_onscene_vertices(text,x0,r,siz,nvert,vert,centered)
    use interfaces_cimgui
    use gui_main, only: g, fontbakesize_large
    use types, only: realloc
    use param, only: newline
    character(len=*), intent(in) :: text
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: r
    real(c_float), intent(in) :: siz
    integer(c_int), intent(inout) :: nvert
    real(c_float), allocatable, intent(inout) :: vert(:,:)
    logical, intent(in), optional :: centered

    integer :: i, j, nline, nvert0
    type(c_ptr) :: cptr
    type(ImFontGlyph), pointer :: glyph
    real(c_float) :: xpos, ypos, lheight
    logical :: centered_
    real(c_float), allocatable :: xlen(:)
    integer, allocatable :: jlen(:)

    real(c_float), parameter :: rshift = 0.01_c_float

    ! initialize
    centered_ = .false.
    if (present(centered)) centered_ = centered
    if (.not.allocated(vert)) then
       allocate(vert(8,100))
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
          vert(4,j) = r + rshift
       end do

       vert(5:6,nvert+1) = (/xpos + glyph%X0, ypos + glyph%Y1/)
       vert(5:6,nvert+2) = (/xpos + glyph%X0, ypos + glyph%Y0/)
       vert(5:6,nvert+3) = (/xpos + glyph%X1, ypos + glyph%Y0/)
       vert(5:6,nvert+4) = (/xpos + glyph%X0, ypos + glyph%Y1/)
       vert(5:6,nvert+5) = (/xpos + glyph%X1, ypos + glyph%Y0/)
       vert(5:6,nvert+6) = (/xpos + glyph%X1, ypos + glyph%Y1/)

       vert(7:8,nvert+1) = (/glyph%U0, glyph%V1/)
       vert(7:8,nvert+2) = (/glyph%U0, glyph%V0/)
       vert(7:8,nvert+3) = (/glyph%U1, glyph%V0/)
       vert(7:8,nvert+4) = (/glyph%U0, glyph%V1/)
       vert(7:8,nvert+5) = (/glyph%U1, glyph%V0/)
       vert(7:8,nvert+6) = (/glyph%U1, glyph%V1/)
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
          vert(5,jlen(i):jlen(i+1)-1) = vert(5,jlen(i):jlen(i+1)-1) - 0.5_c_float * xlen(i)
       end do
       vert(6,jlen(1):nvert) = vert(6,jlen(1):nvert) - 0.5_c_float * nline * lheight
    end if
    vert(5:6,nvert0+1:nvert) = vert(5:6,nvert0+1:nvert) * siz

  end subroutine calc_text_onscene_vertices

end submodule proc
